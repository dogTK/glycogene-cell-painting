"""JUMP-CP: is a gene's ConA / WGA staining different from the other genes perturbed on the same plates?

Definition agreed with the user (2026-11-17): 'stronger' = the gene's wells are brighter than the wells of the other
genes on the same plates (Welch's t-test, Benjamini-Hochberg over glycogenes, FDR < 0.05); 'weaker' = darker. Comparing
with other perturbed genes cancels the shift that introducing an ORF or cutting with CRISPR causes in every well
(against negative controls, 18-25% of non-glycogenes had p < 0.05).
Wells: plate-wise robust z against negative controls; well-position correction (per modality, the median over plates
of each well position is subtracted from every feature, as in the JUMP profiling recipe, because replicate wells of a
gene sit at the same position on replicate plates); then well-level OLS adjustment for cell area, RNA, DNA and Mito
mean intensities and the cell-count proxy (as in well_level_test.py). Negative-control wells are left out of the
reference. Effect = difference of means in units of the reference SD.
Outputs per-gene tests for all genes, glycogene calls, glycogenes vs other genes (nominal p < 0.05, Fisher), and
GlycoEnzOnto over-representation of each set.
Optional: --exclude-low-infection <dir> drops the ORF wells that the JUMP-CP genetic paper excluded for low
infection efficiency or outlying values (jump-cellpainting/2024_Chandrasekaran_Morphmap
00.download-and-process-annotations: output/orf-reagents-low-infection-efficiency-and-outliers.csv.gz, matched on
plate map (input/experiment-metadata.tsv), well and broad_sample), before any normalisation (user, 2026-11-22).
Usage: python vs_other_genes_test.py <out_dir> [--exclude-low-infection <dir>]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cursor  # noqa: E402
from jump_map import load, normalise  # noqa: E402
from well_level_test import LECTIN, adjust_wells, enrich  # noqa: E402

FDR = 0.05


def tests(w):
    w = w[~w.negcon & w.gene.notna()]
    plate_idx = {p: d.index for p, d in w.groupby('plate')}
    vals = {lec: w[f'{lec}_res'] for lec in LECTIN}
    count = w.NUCLEI_NUMBER_OBJECT_NUMBER
    rows = []
    for gene, gw in w.groupby('gene'):
        ref_idx = np.concatenate([plate_idx[p] for p in gw.plate.unique()])
        ref_idx = ref_idx[~np.isin(ref_idx, gw.index)]
        row = {'gene': gene, 'n_wells': len(gw), 'n_ref': len(ref_idx)}
        for lec in LECTIN:
            x, y = vals[lec].loc[gw.index], vals[lec].loc[ref_idx]
            p = stats.ttest_ind(x, y, equal_var=False).pvalue if len(x) >= 2 else np.nan
            row[f'{lec}_diff_sd'] = (x.mean() - y.mean()) / y.std()
            row[f'{lec}_p'] = p
        row['count_diff_sd'] = (count.loc[gw.index].mean() - count.loc[ref_idx].mean()) / count.loc[ref_idx].std()
        rows.append(row)
    return pd.DataFrame(rows).set_index('gene')


def drop_low_infection(cur, wells, d):
    """Drop ORF wells listed as low infection efficiency or outliers by the JUMP-CP genetic paper."""
    cur.execute('select METADATA_JCP2022, METADATA_BROAD_SAMPLE from RAW.GLYCO_CELL_PAINTING.JUMP_ORF')
    bs = dict(cur.fetchall())
    meta = pd.read_csv(d / 'experiment-metadata.tsv', sep='\t')
    pmap = dict(zip(meta.Assay_Plate_Barcode, meta.Plate_Map_Name))
    bad = pd.read_csv(d / 'orf-reagents-low-infection-efficiency-and-outliers.csv.gz')
    bad = set(zip(bad.Metadata_plate_map_name, bad.Metadata_Well, bad.Metadata_broad_sample))
    key = list(zip(wells.plate.map(pmap), wells.well, wells.jcp.map(bs)))
    drop = np.array([k in bad for k in key])
    print(f'low infection efficiency / outlier ORF wells dropped: {drop.sum()} of {wells.plate.map(pmap).notna().sum()}'
          ' ORF wells')
    return wells[~drop].copy()


def main(out, infection_dir=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    pd.set_option('display.width', 220)
    cur = cursor()
    wells, ctrl, genes, ptype, glyco, gmt, _ = load(cur)
    if infection_dir:
        wells = drop_low_infection(cur, wells, Path(infection_dir))
    wells = normalise(wells, ctrl)
    summary, enr = [], []
    for t in ('ORF', 'CRISPR'):
        w = wells[wells.plate.map(ptype).eq(t)].copy()
        w['gene'] = w.jcp.map(genes[t])
        feats = [c for c in w.columns if c.startswith(('CELLS_', 'NUCLEI_'))]
        w[feats] = w[feats] - w.groupby('well')[feats].transform('median')  # well-position correction
        w = adjust_wells(w)
        r = tests(w)
        r['glycogene'] = r.index.isin(glyco)
        g = r[r.glycogene].copy()
        print(f'\n===== {t}: genes {len(r)}, glycogenes {len(g)}, median wells/gene {r.n_wells.median():.0f}, '
              f'median reference wells {r.n_ref.median():.0f}')
        for lec in LECTIN:
            ok = g[f'{lec}_p'].notna()
            g.loc[ok, f'{lec}_fdr'] = multipletests(g.loc[ok, f'{lec}_p'], method='fdr_bh')[1]
            sig = g[f'{lec}_fdr'] < FDR
            g[f'{lec}_call'] = np.select([sig & (g[f'{lec}_diff_sd'] > 0), sig & (g[f'{lec}_diff_sd'] < 0)],
                                         ['stronger', 'weaker'], 'no change')
            c = g[f'{lec}_call'].value_counts()
            other = r.loc[~r.glycogene, f'{lec}_p'].dropna()
            a, b = (g[f'{lec}_p'] < 0.05).sum(), (g[f'{lec}_p'] >= 0.05).sum()
            cc, dd = (other < 0.05).sum(), (other >= 0.05).sum()
            odds, p = fisher_exact([[a, b], [cc, dd]])
            summary.append((t, lec, len(g), c.get('stronger', 0), c.get('weaker', 0), c.get('no change', 0),
                            a / (a + b), cc / (cc + dd), odds, p))
            print(f'{t} {lec}: stronger {c.get("stronger", 0)}, weaker {c.get("weaker", 0)}, no change '
                  f'{c.get("no change", 0)} (of {len(g)}); nominal p<0.05 glycogenes {a / (a + b):.1%} vs other '
                  f'genes {cc / (cc + dd):.1%} (OR {odds:.2f}, p {p:.2g})')
            for call in ('stronger', 'weaker'):
                s = g[g[f'{lec}_call'] == call].sort_values(f'{lec}_diff_sd', ascending=call == 'weaker')
                print(f'   {call}: ' + ', '.join(f'{i} ({v:+.1f} SD)' for i, v in zip(s.index, s[f'{lec}_diff_sd'])))
                e = enrich(set(s.index), set(g.index), gmt).assign(perturbation=t, lectin=lec, call=call)
                enr.append(e)
                if len(s):
                    b0 = e.iloc[0]
                    print(f'   enrichment ({call}): FDR<0.05 {(e.fdr < 0.05).sum()}; best {b0.pathway} '
                          f'{b0.n_hit}/{b0.n_pathway_genes} p={b0.p:.3g} FDR={b0.fdr:.2f} [{b0.hit_genes}]')
        r.to_csv(out / f'{t}_tests_all_genes.csv')
        g.to_csv(out / f'{t}_glycogene_calls.csv')
        big = g[(g.ConA_diff_sd.abs() >= 3) | (g.WGA_diff_sd.abs() >= 3)]
        print('   glycogenes with |effect| >= 3 SD:', ', '.join(
            f'{i} (ConA {a:+.1f}, FDR {fa:.2f}; WGA {b:+.1f}, FDR {fb:.2f}; n={n})' for i, a, fa, b, fb, n in
            zip(big.index, big.ConA_diff_sd, big.ConA_fdr, big.WGA_diff_sd, big.WGA_fdr, big.n_wells)))
    pd.DataFrame(summary, columns=['perturbation', 'lectin', 'n_glycogenes', 'stronger', 'weaker', 'no_change',
                                   'nominal_glyco', 'nominal_other', 'odds_ratio', 'fisher_p']).to_csv(
        out / 'summary.csv', index=False)
    pd.concat(enr).to_csv(out / 'enrichment.csv', index=False)


if __name__ == '__main__':
    args = sys.argv[1:]
    inf = None
    if '--exclude-low-infection' in args:
        k = args.index('--exclude-low-infection')
        inf = args[k + 1]
        del args[k:k + 2]
    main(args[0], inf)
