"""JUMP-CP: is each gene's ConA / WGA staining different from negative controls? (well-level test, 2026-11-17)

Replaces the earlier 'adjusted z >= 2 across genes' call with a test against negative controls, as agreed with the
user: 'stronger' = the gene's wells are brighter than the negative-control wells (FDR < 0.05), 'weaker' = darker.
1. Each feature -> robust z per plate against the plate's negative-control wells (as in jump_map.py).
2. Well-level adjustment (per modality, over all wells including controls): OLS of the ConA (Cells ER mean
   intensity) and WGA (Cells AGP mean intensity) z on cell area, RNA, DNA and Mito mean intensities and the cell-count
   proxy (Nuclei_Number_Object_Number); the residual is the adjusted staining of the well.
3. For each gene: its wells vs the negative-control wells of the same plates (pooled over those plates; if fewer than
   8, the negative-control wells of the same batches). Welch's t-test, difference of means in units of the control
   SD, Benjamini-Hochberg over genes within modality and lectin.
4. Counts of glycogenes stronger / weaker / not different; glycogenes vs other genes (Fisher); GlycoEnzOnto
   over-representation of each set (one-sided Fisher, BH; universe = glycogenes tested in that modality).
Usage: python well_level_test.py <out_dir>
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
from jump_map import FEATS, load, normalise  # noqa: E402

COV = ['CELLS_AREASHAPE_AREA', 'CELLS_INTENSITY_MEANINTENSITY_RNA', 'CELLS_INTENSITY_MEANINTENSITY_DNA',
       'CELLS_INTENSITY_MEANINTENSITY_MITO', 'NUCLEI_NUMBER_OBJECT_NUMBER']
LECTIN = {'ConA': 'CELLS_INTENSITY_MEANINTENSITY_ER', 'WGA': 'CELLS_INTENSITY_MEANINTENSITY_AGP'}
FDR = 0.05
MIN_CTRL = 8


def adjust_wells(w):
    ok = w[COV + list(LECTIN.values())].notna().all(axis=1) & np.isfinite(w[COV + list(LECTIN.values())]).all(axis=1)
    w = w[ok].copy()
    X = np.column_stack([np.ones(len(w)), w[COV].to_numpy()])
    for lec, col in LECTIN.items():
        beta, *_ = np.linalg.lstsq(X, w[col].to_numpy(), rcond=None)
        w[f'{lec}_res'] = w[col].to_numpy() - X @ beta
    return w


def test_genes(w, batch):
    ctrl = w[w.negcon]
    by_plate = {p: d for p, d in ctrl.groupby('plate')}
    by_batch = {b: d for b, d in ctrl.groupby(ctrl.plate.map(batch))}
    rows = []
    for gene, gw in w[~w.negcon & w.gene.notna()].groupby('gene'):
        plates = gw.plate.unique()
        c = pd.concat([by_plate[p] for p in plates if p in by_plate]) if any(p in by_plate for p in plates) else ctrl.iloc[:0]
        if len(c) < MIN_CTRL:
            batches = {batch.get(p) for p in plates}
            c = pd.concat([by_batch[b] for b in batches if b in by_batch])
        row = {'gene': gene, 'n_wells': len(gw), 'n_ctrl': len(c)}
        for lec in LECTIN:
            x, y = gw[f'{lec}_res'], c[f'{lec}_res']
            t, p = stats.ttest_ind(x, y, equal_var=False) if len(x) >= 2 else (np.nan, np.nan)
            row[f'{lec}_diff_sd'] = (x.mean() - y.mean()) / y.std()
            row[f'{lec}_p'] = p
        row['count_diff_sd'] = (gw.NUCLEI_NUMBER_OBJECT_NUMBER.mean() - c.NUCLEI_NUMBER_OBJECT_NUMBER.mean()) / \
            c.NUCLEI_NUMBER_OBJECT_NUMBER.std()
        rows.append(row)
    r = pd.DataFrame(rows).set_index('gene')
    for lec in LECTIN:
        ok = r[f'{lec}_p'].notna()
        r.loc[ok, f'{lec}_fdr'] = multipletests(r.loc[ok, f'{lec}_p'], method='fdr_bh')[1]
        r[f'{lec}_call'] = np.select([(r[f'{lec}_fdr'] < FDR) & (r[f'{lec}_diff_sd'] > 0),
                                      (r[f'{lec}_fdr'] < FDR) & (r[f'{lec}_diff_sd'] < 0)], ['stronger', 'weaker'],
                                     'no change')
    return r


def enrich(hits, universe, gmt):
    rows = []
    for n, s in gmt.items():
        s = s & universe
        if len(s) < 3 or len(s) >= len(universe):
            continue
        a = len(hits & s)
        b, c = len(hits - s), len(s - hits)
        p = fisher_exact([[a, b], [c, len(universe) - a - b - c]], alternative='greater')[1] if a else 1.0
        rows.append((n, len(s), a, p, ', '.join(sorted(hits & s))))
    e = pd.DataFrame(rows, columns=['pathway', 'n_pathway_genes', 'n_hit', 'p', 'hit_genes'])
    e['fdr'] = multipletests(e.p, method='fdr_bh')[1]
    return e.sort_values('p')


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    pd.set_option('display.width', 220, 'display.max_colwidth', 90)
    cur = cursor()
    wells, ctrl, genes, ptype, glyco, gmt, lectin = load(cur)
    cur.execute('select METADATA_PLATE, METADATA_BATCH from RAW.GLYCO_CELL_PAINTING.JUMP_PLATE')
    batch = dict(cur.fetchall())
    wells = normalise(wells, ctrl)
    summary, enr = [], []
    for t in ('ORF', 'CRISPR'):
        w = wells[wells.plate.map(ptype).eq(t)].copy()
        w['gene'] = w.jcp.map(genes[t])
        nper = w[w.negcon].groupby('plate').size()
        print(f'\n===== {t}: wells {len(w)}, negcon wells {int(w.negcon.sum())}; negcon per plate median '
              f'{nper.median():.0f} (min {nper.min()}, plates without negcon '
              f'{len(set(w.plate)) - len(nper)})')
        w = adjust_wells(w)
        r = test_genes(w, batch)
        r['glycogene'] = r.index.isin(glyco)
        r.to_csv(out / f'{t}_well_level_tests.csv')
        for lec in LECTIN:
            g = r[r.glycogene]
            c = g[f'{lec}_call'].value_counts()
            ch_g = (g[f'{lec}_call'] != 'no change')
            ch_o = (r.loc[~r.glycogene, f'{lec}_call'] != 'no change')
            a, b, cc, d = ch_g.sum(), (~ch_g).sum(), ch_o.sum(), (~ch_o).sum()
            odds, p = fisher_exact([[a, b], [cc, d]])
            summary.append((t, lec, len(g), c.get('stronger', 0), c.get('weaker', 0), c.get('no change', 0),
                            a / (a + b), cc / (cc + d), odds, p))
            print(f'{t} {lec}: glycogenes {len(g)} -> stronger {c.get("stronger", 0)}, weaker {c.get("weaker", 0)}, '
                  f'no change {c.get("no change", 0)}; changed {a / (a + b):.1%} vs other genes {cc / (cc + d):.1%} '
                  f'(OR {odds:.2f}, p {p:.2g})')
            for call in ('stronger', 'weaker'):
                hits = set(g.index[g[f'{lec}_call'] == call])
                top = g[g[f'{lec}_call'] == call].sort_values(f'{lec}_diff_sd', ascending=call == 'weaker')
                print(f'   {call}: {", ".join(f"{i} ({v:+.1f})" for i, v in top[f"{lec}_diff_sd"].head(12).items())}')
                e = enrich(hits, set(g.index), gmt).assign(perturbation=t, lectin=lec, call=call)
                enr.append(e)
                sig = e[e.fdr < 0.05]
                best = e.iloc[0] if len(e) and e.iloc[0].n_hit else None
                print(f'   enrichment ({call}): FDR<0.05 {len(sig)}' + (
                    f'; best {best.pathway} {best.n_hit}/{best.n_pathway_genes} p={best.p:.3g} FDR={best.fdr:.2f}'
                    if best is not None else ''))
    pd.DataFrame(summary, columns=['perturbation', 'lectin', 'n_glycogenes', 'stronger', 'weaker', 'no_change',
                                   'changed_glyco', 'changed_other', 'odds_ratio', 'fisher_p']).to_csv(
        out / 'summary.csv', index=False)
    pd.concat(enr).to_csv(out / 'enrichment.csv', index=False)


if __name__ == '__main__':
    main(sys.argv[1])
