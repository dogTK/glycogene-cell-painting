"""Well-level values for the image panels (user, 2026-11-23): Figure 3c (overexpression and knockout),
genes chosen by select_genes.py, each against negative controls of the same plate.

Reuses the main analysis (scripts/13_jump: jump_map.load / normalise, vs_other_genes_test.drop_low_infection,
well_level_test.adjust_wells) and reports, for every well of the target genes and every negative-control well on the
same plates, four stages of the ER (ConA), AGP (WGA) and DNA cell mean intensity:
  raw   : CellProfiler per-well mean of Cells_Intensity_MeanIntensity_<ch> (images as acquired, illumination-corrected
          by the JUMP pipeline)
  z     : robust z against the plate's negative-control wells (median, 1.4826 x MAD)
  z_pos : z minus the median of the same well position over all plates of the modality (well-position correction)
  adj   : residual of z_pos regressed on cell area, RNA, DNA, Mito and cell count (ConA, WGA only; the values tested)
and the Glass's delta of each target well against the other genes' wells on the same plate (the comparison used in
the manuscript), to check that the per-well numbers reproduce the published effect.
ORF wells flagged as low infection efficiency / outliers (Morphmap list) are marked and dropped from the tests, as in
the main analysis.
Output: <out>/well_values.csv (one row per target or negative-control well on the target plates),
        <out>/plate_summary.csv (per target well: value, negcon median, delta vs other genes, per stage and channel).
Usage: python well_values.py <out_dir> <panel: oe | ko>  (genes from select_genes.py)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '13_jump'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import cursor  # noqa: E402
from jump_map import load, normalise  # noqa: E402
from vs_other_genes_test import drop_low_infection  # noqa: E402
from well_level_test import adjust_wells  # noqa: E402

from panels import panel_genes  # noqa: E402
INF = Path.home() / 'research/glyco_cell_painting/data/morphmap_infection'
CH = {'ER': 'CELLS_INTENSITY_MEANINTENSITY_ER', 'AGP': 'CELLS_INTENSITY_MEANINTENSITY_AGP',
      'DNA': 'CELLS_INTENSITY_MEANINTENSITY_DNA'}


def main(out, panel):
    TARGETS = {g: m for g, m, _, _ in panel_genes(panel)}
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    cur = cursor()
    wells, ctrl, genes, ptype, glyco, gmt, _ = load(cur)
    raw = wells.copy()
    kept = drop_low_infection(cur, wells, INF)
    low_inf = set(zip(wells.plate, wells.well)) - set(zip(kept.plate, kept.well))
    z = normalise(kept, ctrl)
    rows, summ = [], []
    for t in ('ORF', 'CRISPR'):
        w = z[z.plate.map(ptype).eq(t)].copy()
        w['gene'] = w.jcp.map(genes[t])
        feats = [c for c in w.columns if c.startswith(('CELLS_', 'NUCLEI_'))]
        zz = w[feats].copy()
        w[feats] = w[feats] - w.groupby('well')[feats].transform('median')
        zpos = w[feats].copy()
        w = adjust_wells(w)
        for gene, mod in TARGETS.items():
            if mod != t:
                continue
            plates = sorted(w.loc[w.gene == gene, 'plate'].unique())
            allp = sorted(raw.loc[raw.jcp.map(genes[t]).eq(gene), 'plate'].unique())
            for p in allp:
                rp = raw[raw.plate == p]
                for _, r in rp[rp.jcp.map(genes[t]).eq(gene) | rp.jcp.map(ctrl).eq('negcon')].iterrows():
                    key = (p, r.well)
                    d = dict(target=gene, modality=t, plate=p, well=r.well, jcp=r.jcp,
                             role='target' if genes[t].get(r.jcp) == gene else 'negcon',
                             negcon_name=r.jcp if genes[t].get(r.jcp) != gene else '', low_infection=key in low_inf)
                    for ch, col in CH.items():
                        d[f'raw_{ch}'] = r[col]
                    m = (w.plate == p) & (w.well == r.well)
                    if key not in low_inf and ((zz.index.isin(w.index[m])).any() or m.any()):
                        idx = w.index[m]
                        for ch, col in CH.items():
                            d[f'z_{ch}'] = zz.loc[idx, col].iloc[0] if len(idx) else np.nan
                            d[f'zpos_{ch}'] = zpos.loc[idx, col].iloc[0] if len(idx) else np.nan
                        d['adj_ER'] = w.loc[idx, 'ConA_res'].iloc[0] if len(idx) else np.nan
                        d['adj_AGP'] = w.loc[idx, 'WGA_res'].iloc[0] if len(idx) else np.nan
                    rows.append(d)
            # Glass's delta per target well vs other genes' wells on the same plate (adjusted values)
            for p in plates:
                wp = w[w.plate == p]
                tw = wp[wp.gene == gene]
                other = wp[~wp.negcon & wp.gene.notna() & (wp.gene != gene)]
                neg = wp[wp.negcon]
                for lec, ch in (('ConA', 'ER'), ('WGA', 'AGP')):
                    col = f'{lec}_res'
                    summ.append(dict(target=gene, plate=p, channel=ch, target_adj=tw[col].mean(),
                                     negcon_adj_median=neg[col].median(), n_negcon=len(neg),
                                     other_mean=other[col].mean(), other_sd=other[col].std(),
                                     delta_vs_other=(tw[col].mean() - other[col].mean()) / other[col].std()))
    wv = pd.DataFrame(rows)
    sm = pd.DataFrame(summ)
    wv.to_csv(out / 'well_values.csv', index=False)
    sm.to_csv(out / 'plate_summary.csv', index=False)
    pd.set_option('display.width', 250)
    print(sm.round(2).to_string(index=False))
    print('\nPooled delta (mean over plates of per-plate delta):')
    print(sm.groupby(['target', 'channel']).delta_vs_other.agg(['mean', 'median', 'count']).round(2))
    print('\nNegcon wells per target plate:')
    print(wv[wv.role == 'negcon'].groupby(['target', 'plate']).size().to_string())
    print('\nLow-infection flags among target wells:', wv[(wv.role == 'target') & wv.low_infection][['target', 'plate', 'well']].to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:3])
