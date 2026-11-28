"""Per-gene Glass's delta for all five Cell Painting channels in JUMP-CP (user, 2026-11-22: are the glycogene effects on
ConA / WGA selective compared with the DNA, RNA and Mito channels?).

Same preprocessing as scripts/13_jump/vs_other_genes_test.py --exclude-low-infection:
- low infection efficiency / outlier ORF wells dropped (JUMP-CP genetic paper list),
- plate-wise robust z against negative controls (jump_map.normalise),
- well-position correction (median over plates of each well position subtracted, per perturbation),
- well-level OLS adjustment, residual = adjusted signal,
- each gene's wells vs the wells of the other genes on the same plates (negative controls excluded), Welch's t-test,
  effect = difference of means / SD of the reference wells (Glass's delta).
Adjustment models
- 'main' (as close as possible to the existing ConA / WGA model): every channel is regressed on cell area, the
  cell-count proxy and the mean intensities of the non-lectin channels other than itself:
  ConA (ER) and WGA (AGP) on area + RNA + DNA + Mito + count (identical to the existing analysis);
  DNA on area + RNA + Mito + count; RNA on area + DNA + Mito + count; Mito on area + DNA + RNA + count.
- 'symmetric' (sensitivity): every channel on area + count + the other four channels.
Outputs per perturbation: <ORF|CRISPR>_<model>_effects.csv (gene, n_wells, <channel>_delta, <channel>_p, glycogene).
Usage: python channel_effects.py <out_dir> <infection_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '13_jump'))
from common import cursor  # noqa: E402
from jump_map import load, normalise  # noqa: E402
from vs_other_genes_test import drop_low_infection  # noqa: E402

COL = {'ConA': 'CELLS_INTENSITY_MEANINTENSITY_ER', 'WGA': 'CELLS_INTENSITY_MEANINTENSITY_AGP',
       'DNA': 'CELLS_INTENSITY_MEANINTENSITY_DNA', 'RNA': 'CELLS_INTENSITY_MEANINTENSITY_RNA',
       'Mito': 'CELLS_INTENSITY_MEANINTENSITY_MITO'}
AREA, COUNT = 'CELLS_AREASHAPE_AREA', 'NUCLEI_NUMBER_OBJECT_NUMBER'
NONLECTIN = ['DNA', 'RNA', 'Mito']
MODELS = {
    'main': {ch: [AREA, COUNT] + [COL[o] for o in NONLECTIN if o != ch] for ch in COL},
    'symmetric': {ch: [AREA, COUNT] + [COL[o] for o in COL if o != ch] for ch in COL},
}
CONTROLS = {'eGFP', 'BFP', 'HcRed', 'LacZ', 'LUCIFERASE', 'Luciferase'}


def adjust(w, model):
    cols = sorted({c for cov in MODELS[model].values() for c in cov} | set(COL.values()))
    ok = w[cols].notna().all(axis=1) & np.isfinite(w[cols]).all(axis=1)
    w = w[ok].copy()
    for ch, cov in MODELS[model].items():
        X = np.column_stack([np.ones(len(w)), w[cov].to_numpy()])
        beta, *_ = np.linalg.lstsq(X, w[COL[ch]].to_numpy(), rcond=None)
        w[f'{ch}_res'] = w[COL[ch]].to_numpy() - X @ beta
    return w


def tests(w):
    w = w[~w.negcon & w.gene.notna()]
    plate_idx = {p: d.index for p, d in w.groupby('plate')}
    vals = {ch: w[f'{ch}_res'] for ch in COL}
    rows = []
    for gene, gw in w.groupby('gene'):
        ref = np.concatenate([plate_idx[p] for p in gw.plate.unique()])
        ref = ref[~np.isin(ref, gw.index)]
        row = {'gene': gene, 'n_wells': len(gw)}
        for ch in COL:
            x, y = vals[ch].loc[gw.index], vals[ch].loc[ref]
            row[f'{ch}_delta'] = (x.mean() - y.mean()) / y.std()
            row[f'{ch}_p'] = stats.ttest_ind(x, y, equal_var=False).pvalue if len(x) >= 2 else np.nan
        rows.append(row)
    return pd.DataFrame(rows).set_index('gene')


def main(out, infection_dir):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    cur = cursor()
    wells, ctrl, genes, ptype, glyco, _, _ = load(cur)
    wells = drop_low_infection(cur, wells, Path(infection_dir))
    wells = normalise(wells, ctrl)
    for t in ('ORF', 'CRISPR'):
        w = wells[wells.plate.map(ptype).eq(t)].copy()
        w['gene'] = w.jcp.map(genes[t])
        feats = [c for c in w.columns if c.startswith(('CELLS_', 'NUCLEI_'))]
        w[feats] = w[feats] - w.groupby('well')[feats].transform('median')  # well-position correction
        for model in MODELS:
            r = tests(adjust(w, model))
            r = r[~r.index.isin(CONTROLS)]
            r['glycogene'] = r.index.isin(glyco)
            r.to_csv(out / f'{t}_{model}_effects.csv')
            print(f'{t} {model}: genes {len(r)}, glycogenes {r.glycogene.sum()}')


if __name__ == '__main__':
    main(*sys.argv[1:3])
