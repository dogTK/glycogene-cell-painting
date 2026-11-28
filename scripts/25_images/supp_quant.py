"""Supplementary figure for the image panels Figure 3c: overexpression (rows a, b) and knockout (rows c, d): quantitative check of
the images (user, 2026-11-23). Genes from select_genes.py (panels.py).

a, c: CellProfiler well-level cell mean intensity of the channel of interest (raw, as measured by the JUMP
   pipeline), divided by the median of the negative-control wells of the same plate, for the 5 replicate plates; grey =
   every eligible negative-control well (ORF: BFP, HcRed, LUCIFERASE, LacZ; CRISPR: non-targeting), coloured = the
   target well. Source: well_values.csv.
b, d: image-level cell-region mean of the channel of interest for every QC-passing site (fetch_images.py
   measurements, raw pixels), divided by the median of the control well's passing sites on the same plate; grey =
   control-well sites, coloured = target-well sites; the sites shown in Figure 3c are outlined in black.
Usage: python supp_quant.py <results/25_images dir> <out stem>  (reads <dir>/oe and <dir>/ko)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_panel import choose_site  # noqa: E402
from panels import panel_genes  # noqa: E402

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'pdf.fonttype': 42, 'axes.linewidth': 0.6})
COL = {'up': '#1f66c2', 'down': '#d9541e'}
NEG = {'ORF': {'JCP2022_915128', 'JCP2022_915129', 'JCP2022_915130', 'JCP2022_915131'}, 'KO': {'JCP2022_800002'}}
CH_LABEL = {'AGP': 'AGP (WGA + phalloidin)', 'ER': 'ER (ConA)'}


def main(r, stem):
    r = Path(r)
    fig, axs = plt.subplots(4, 3, figsize=(7.1, 7.9), gridspec_kw=dict(wspace=0.32, hspace=0.62))
    rng = np.random.default_rng(1)
    for k, panel in enumerate(('oe', 'ko')):
        draw(r / panel, panel, axs[2 * k:2 * k + 2], rng)
    for i, lab in enumerate('ABCD'):
        axs[i, 0].text(-0.42, 1.12, lab, transform=axs[i, 0].transAxes, fontsize=9, fontweight='bold')
    for ext in ('png', 'pdf'):
        fig.savefig(f'{stem}.{ext}', dpi=600 if ext == 'png' else None, bbox_inches='tight')


def draw(r, panel, axs, rng):
    wv = pd.read_csv(r / 'well_values.csv')
    sm = pd.read_csv(r / 'site_measurements.csv')
    sel = pd.read_csv(r / 'image_selection.csv')
    for j, (gene, mod, ch, direction) in enumerate(panel_genes(panel)):
        kind, col = ('ORF' if mod == 'ORF' else 'KO'), COL[direction]
        d = wv[wv.target == gene]
        plates = sorted(d.plate.unique())
        ax = axs[0, j]
        for i, p in enumerate(plates):
            n = d[(d.plate == p) & (d.role == 'negcon') & d.jcp.isin(NEG[kind]) & ~d.low_infection]
            t = d[(d.plate == p) & (d.role == 'target')].iloc[0]
            med = n[f'raw_{ch}'].median()
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(n)), n[f'raw_{ch}'] / med, s=7, color='#b5b3ac', lw=0)
            ax.scatter(i, t[f'raw_{ch}'] / med, s=22, color=col, edgecolor='black', lw=0.4, zorder=3)
        ax.axhline(1, color='#55534f', lw=0.5, ls=':')
        ax.set_xticks(range(len(plates)), [p.replace('CP-CC9-', '') for p in plates], rotation=45, ha='right',
                      fontsize=5.5)
        ax.set_ylabel(f'{CH_LABEL[ch]}\nwell mean / control median', fontsize=6)
        ax.set_title(f'{gene} {"ORF" if kind == "ORF" else "KO"}', fontsize=7, fontweight='bold', loc='left')
        ax = axs[1, j]
        sg = sm[sm.target == gene]
        for i, p in enumerate(plates):
            dp = sg[sg.plate == p]
            if dp.empty:
                continue
            q = {role: choose_site(dp[dp.role == role], ch)[1] for role in ('control', 'target')}
            med = q['control'].loc[q['control'].qc_pass, f'{ch}_cell_mean'].median()
            for role, c, off in (('control', '#b5b3ac', -0.15), ('target', col, 0.15)):
                v = q[role][q[role].qc_pass]
                ax.scatter(i + off + rng.uniform(-0.07, 0.07, len(v)), v[f'{ch}_cell_mean'] / med, s=8, color=c, lw=0)
                shown = sel[(sel.target == gene) & (sel.role == role) & (sel.plate == p)]
                for _, s in shown.iterrows():
                    vv = v[(v.well == s.well) & (v.site == s.site)]
                    ax.scatter(i + off, vv[f'{ch}_cell_mean'] / med, s=26, facecolor='none', edgecolor='black',
                               lw=0.8, zorder=4)
        ax.axhline(1, color='#55534f', lw=0.5, ls=':')
        ax.set_xticks(range(len(plates)), [p.replace('CP-CC9-', '') for p in plates], rotation=45, ha='right',
                      fontsize=5.5)
        ax.set_ylabel('site cell-region mean /\ncontrol-site median', fontsize=6)
        for a in axs[:, j]:
            for sp in ('top', 'right'):
                a.spines[sp].set_visible(False)


if __name__ == '__main__':
    main(*sys.argv[1:3])
