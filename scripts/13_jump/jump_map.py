"""JUMP-CP: which glycogene perturbations show up in ConA / WGA staining? (plan fixed in docs/progress.md, 2026-11-17)

Data: RAW.GLYCO_CELL_PAINTING.JUMP_GENETIC_WELL_FEATURES, JUMP_WELL, JUMP_ORF, JUMP_CRISPR, JUMP_PERTURBATION_CONTROL.
1. Per plate, each feature -> robust z against the plate's negative-control wells (median, 1.4826 x MAD; all wells
   if a plate has < 8 negcon wells). Per gene, median over all its wells (replicates, plates, ORF clones). ORF and
   CRISPR separately.
2. Adjusted ConA (Cells ER mean intensity) and WGA (Cells AGP mean intensity): OLS residual on cell area, RNA, DNA and
   Mito mean intensities and the cell-count proxy (Nuclei_Number_Object_Number), fitted over all genes of the
   modality, as a robust z.
3. Primary hypothesis (one-sided): ORF of TFEB, GNPTAB, GNPTG, NAGPA, M6PR, HEXB, NEU1, GLB1, GNS, FUCA2, LAMP1, LAMP2,
   CTSD raises adjusted ConA vs all other ORF genes (Welch's t-test, one-sided); difference with 95% CI.
4. Descriptive per-glycogene analysis: glycogenes with |adjusted z| >= 2 'show up'. Toxic flag: robust z of the cell-count
   proxy across all genes of the modality <= -2 (changed on 2026-11-17 from the planned z vs negative controls, which
   flagged 5,145 of 7,975 CRISPR genes because CRISPR cutting lowers counts relative to non-targeting guides; the
   adjustment, which already includes the cell count, is unchanged). GlycoEnzOnto
   over-representation of the genes that show up (per lectin, direction, modality; one-sided Fisher, BH; universe =
   glycogenes tested in that modality). Agreement with the expected direction of the lectin gene sets.
Optional: --exclude-low-infection <dir> (see vs_other_genes_test.py).
Usage: python jump_map.py <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cursor  # noqa: E402

CH = ['ER', 'AGP', 'RNA', 'DNA', 'MITO']
FEATS = [f'CELLS_INTENSITY_MEANINTENSITY_{c}' for c in CH] + ['CELLS_AREASHAPE_AREA', 'NUCLEI_NUMBER_OBJECT_NUMBER']
LYSO = ['TFEB', 'GNPTAB', 'GNPTG', 'NAGPA', 'M6PR', 'HEXB', 'NEU1', 'GLB1', 'GNS', 'FUCA2', 'LAMP1', 'LAMP2', 'CTSD']
THR = 2.0


def robust_z(x, ref):
    med = np.nanmedian(ref)
    mad = 1.4826 * np.nanmedian(np.abs(ref - med))
    return (x - med) / mad if mad > 0 else x * np.nan


def load(cur):
    cols = ', '.join(f'f.{c}' for c in FEATS)
    cur.execute(f'''select f.METADATA_PLATE, f.METADATA_WELL, w.METADATA_JCP2022, {cols}
                    from RAW.GLYCO_CELL_PAINTING.JUMP_GENETIC_WELL_FEATURES f
                    join RAW.GLYCO_CELL_PAINTING.JUMP_WELL w
                      on w.METADATA_SOURCE = f.METADATA_SOURCE and w.METADATA_PLATE = f.METADATA_PLATE
                     and w.METADATA_WELL = f.METADATA_WELL''')
    wells = pd.DataFrame(cur.fetchall(), columns=['plate', 'well', 'jcp'] + FEATS)
    wells[FEATS] = wells[FEATS].astype(float)
    cur.execute('select METADATA_JCP2022, METADATA_PERT_TYPE from RAW.GLYCO_CELL_PAINTING.JUMP_PERTURBATION_CONTROL')
    ctrl = dict(cur.fetchall())
    genes = {}
    for t in ('ORF', 'CRISPR'):
        cur.execute(f'select METADATA_JCP2022, METADATA_SYMBOL from RAW.GLYCO_CELL_PAINTING.JUMP_{t}')
        genes[t] = dict(cur.fetchall())
    cur.execute('select METADATA_PLATE, METADATA_PLATETYPE from RAW.GLYCO_CELL_PAINTING.JUMP_PLATE')
    ptype = dict(cur.fetchall())
    # Glycogenes = all GlycoEnzOnto genes (403), not the LINCS-matched 385 (user, 2026-11-21)
    cur.execute('select distinct "geneName" from RAW.GLYCODB.GLYCOENZONTO_GLYCOGENE')
    glyco = {r[0] for r in cur.fetchall()}
    cur.execute('select PATHWAY_NAME, GENE_SYMBOLS from RAW.GLYCODB.GLYCOENZONTO_PATHWAY_GMT')
    gmt = {n: set(s.split('|')) for n, s in cur.fetchall()}
    cur.execute('select LECTIN, GENE, EXPECTED_DIRECTION from RAW.GLYCO_CELL_PAINTING.LECTIN_GENE_SETS where _LOAD_ID = '
                '(select max(_LOAD_ID) from RAW.GLYCO_CELL_PAINTING.LECTIN_GENE_SETS)')
    lectin = pd.DataFrame(cur.fetchall(), columns=['lectin', 'gene', 'expected'])
    return wells, ctrl, genes, ptype, glyco, gmt, lectin


def normalise(wells, ctrl):
    wells = wells.copy()
    wells['negcon'] = wells.jcp.map(ctrl).eq('negcon')
    parts = []
    for _, p in wells.groupby('plate'):
        ref = p[p.negcon] if p.negcon.sum() >= 8 else p
        p = p.copy()
        for f in FEATS:
            p[f] = robust_z(p[f].to_numpy(), ref[f].to_numpy())
        parts.append(p)
    return pd.concat(parts)


def adjust(g):
    out = {}
    cov = ['CELLS_AREASHAPE_AREA', 'CELLS_INTENSITY_MEANINTENSITY_RNA', 'CELLS_INTENSITY_MEANINTENSITY_DNA',
           'CELLS_INTENSITY_MEANINTENSITY_MITO', 'NUCLEI_NUMBER_OBJECT_NUMBER']
    ok = g[cov + ['CELLS_INTENSITY_MEANINTENSITY_ER', 'CELLS_INTENSITY_MEANINTENSITY_AGP']].notna().all(axis=1)
    X = np.column_stack([np.ones(ok.sum()), g.loc[ok, cov].to_numpy()])
    for lectin, col in (('ConA', 'CELLS_INTENSITY_MEANINTENSITY_ER'), ('WGA', 'CELLS_INTENSITY_MEANINTENSITY_AGP')):
        y = g.loc[ok, col].to_numpy()
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = pd.Series(np.nan, index=g.index)
        r[ok] = y - X @ beta
        out[f'{lectin}_adj'] = (r - r.median()) / (1.4826 * (r - r.median()).abs().median())
        out[f'{lectin}_r2'] = 1 - np.var(y - X @ beta) / np.var(y)
    return out


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
    r = pd.DataFrame(rows, columns=['pathway', 'n_pathway_genes', 'n_hit', 'p', 'hit_genes'])
    r['fdr'] = multipletests(r.p, method='fdr_bh')[1]  # BH over every tested pathway, including those with no hit
    return r[r.n_hit > 0].sort_values('p')


def main(out, infection_dir=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    pd.set_option('display.width', 220, 'display.max_colwidth', 80)
    cur = cursor()
    wells, ctrl, genes, ptype, glyco, gmt, lectin = load(cur)
    if infection_dir:  # same ORF well exclusion as vs_other_genes_test.py (user, 2026-11-22)
        from vs_other_genes_test import drop_low_infection
        wells = drop_low_infection(cur, wells, Path(infection_dir))
    wells = normalise(wells, ctrl)
    summary = {}
    for t in ('ORF', 'CRISPR'):
        w = wells[wells.plate.map(ptype).eq(t) & ~wells.negcon].copy()
        w['gene'] = w.jcp.map(genes[t])
        w = w.dropna(subset=['gene'])
        g = w.groupby('gene')[FEATS].median()
        g['n_wells'] = w.groupby('gene').size()
        adj = adjust(g)
        g['ConA_adj'], g['WGA_adj'] = adj['ConA_adj'], adj['WGA_adj']
        g['count_z'] = g.NUCLEI_NUMBER_OBJECT_NUMBER
        c = g.count_z
        g['count_rel'] = (c - c.median()) / (1.4826 * (c - c.median()).abs().median())
        g['toxic'] = g.count_rel <= -THR
        g['glycogene'] = g.index.isin(glyco)
        g.to_csv(out / f'{t}_gene_scores.csv')
        print(f'\n===== {t}: genes {len(g)}, glycogenes {g.glycogene.sum()}, median wells/gene {g.n_wells.median():.0f}; '
              f'R2 of covariates: ConA {adj["ConA_r2"]:.2f}, WGA {adj["WGA_r2"]:.2f}; toxic genes {g.toxic.sum()}')

        if t == 'ORF':
            lyso = g.index.isin(LYSO)
            x, y = g.loc[lyso, 'ConA_adj'].dropna(), g.loc[~lyso, 'ConA_adj'].dropna()
            tt = stats.ttest_ind(x, y, equal_var=False, alternative='greater')
            se = np.sqrt(x.var() / len(x) + y.var() / len(y))
            dfree = (x.var() / len(x) + y.var() / len(y)) ** 2 / (
                (x.var() / len(x)) ** 2 / (len(x) - 1) + (y.var() / len(y)) ** 2 / (len(y) - 1))
            half = stats.t.ppf(0.975, dfree) * se
            d = x.mean() - y.mean()
            print(f'PRIMARY HYPOTHESIS (ORF lysosome/M6P genes raise adjusted ConA): n = {len(x)} vs {len(y)}; '
                  f'difference {d:.2f} (95% CI {d - half:.2f} to {d + half:.2f}); one-sided Welch p = {tt.pvalue:.3g}')
            print(g.loc[lyso, ['ConA_adj', 'WGA_adj', 'count_z', 'n_wells']].round(2).sort_values('ConA_adj',
                                                                                                ascending=False).to_string())
            summary['primary'] = dict(n_lyso=len(x), n_other=len(y), difference=d, ci_low=d - half, ci_high=d + half,
                                      p_one_sided=tt.pvalue)

        gg = g[g.glycogene]
        universe = set(gg.index)
        for lec in ('ConA', 'WGA'):
            for direction, sign in (('up', 1), ('down', -1)):
                hits = set(gg.index[(sign * gg[f'{lec}_adj'] >= THR)])
                hits_nt = set(gg.index[(sign * gg[f'{lec}_adj'] >= THR) & ~gg.toxic])
                e = enrich(hits, universe, gmt)
                e.to_csv(out / f'{t}_{lec}_{direction}_enrichment.csv', index=False)
                top = e[e.fdr < 0.05] if len(e) else e
                print(f'{t} {lec} {direction}: glycogenes {len(hits)} (non-toxic {len(hits_nt)}): '
                      f'{", ".join(sorted(hits))[:300]}')
                if len(top):
                    print('   pathways FDR < 0.05:', '; '.join(f'{n} ({a}/{m}, {f:.3f})' for n, a, m, f in
                                                             top[['pathway', 'n_hit', 'n_pathway_genes', 'fdr']].values[:5]))
            ls = lectin[lectin.lectin == lec].drop_duplicates('gene').set_index('gene')
            ls = ls[ls.index.isin(g.index)]
            sgn = np.sign(g.loc[ls.index, f'{lec}_adj'])
            # ORF raises expression: expected sign = expected_direction; CRISPR removes it: expected sign = -direction
            exp = ls.expected if t == 'ORF' else -ls.expected
            agree = (sgn == exp).sum()
            strong = ((sgn == exp) & (g.loc[ls.index, f'{lec}_adj'].abs() >= THR)).sum()
            print(f'   {lec} lectin-set genes in {t}: {len(ls)}; direction as expected {agree}/{len(ls)}; '
                  f'as expected and |z| >= 2: {strong}')
    pd.Series(summary.get('primary', {})).to_csv(out / 'primary_hypothesis.csv')


if __name__ == '__main__':
    args = sys.argv[1:]
    inf = None
    if '--exclude-low-infection' in args:
        k = args.index('--exclude-low-infection')
        inf = args[k + 1]
        del args[k:k + 2]
    main(args[0], inf)
