"""Figure 5: independent public data confirm the knockout results (user, 2026-11-23: "results に入れよ").
a Tsui et al. 2024 (A549 CRISPR knockout screen sorted by the high-mannose lectin HHL): share of Tsui hits (10% FDR,
  knockout increases HHL binding) and of the other genes whose JUMP-CP knockout increased each channel (p < 0.05,
  delta > 0); all genes and GlycoEnzOnto genes; Fisher's exact test (tsui_comparison.py).
b PERISCOPE (A549, HeLa DMEM, HeLa HPLM; WGA and phalloidin in separate channels): share of glycogenes whose |z| exceeds
  the 95th percentile of the other genes, exact 95% CI; filled = binomial p < 0.05 (periscope_channels.py).
c The genes listed in periscope_channels.py before the PERISCOPE data were analysed (all shown): JUMP-CP knockout
  (Glass's delta; AGP for WGA, ER for ConA; no separate phalloidin channel) and PERISCOPE z, WGA, phalloidin and ConA.
No title in the image.
Usage: python external_figure.py <results/24_external dir> <out stem>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from figure_style import (FS_SMALL, INK2, LW_EDGE, LW_ERR, LW_REF, LW_SEP, WIDTH, apply, despine,  # noqa: E402
                          letter, save)
R = Path.home() / 'research/glyco_cell_painting/results'
CHJ = [('ConA', 'ER\n(ConA)'), ('WGA', 'AGP\n(WGA)'), ('DNA', 'DNA'), ('RNA', 'RNA'), ('Mito', 'Mito')]
CHP = [('ConA', 'ConA'), ('WGA', 'WGA'), ('Phalloidin', 'Phalloidin'), ('DAPI', 'DAPI'), ('Mito', 'Mito')]
COND = [('A549', 'A549', '#1f66c2'), ('HeLa_DMEM', 'HeLa (DMEM)', '#d9541e'), ('HeLa_HPLM', 'HeLa (HPLM)', '#2e9e6a')]
GENES = ['MGAT1', 'MAN2A1', 'MAN1A2', 'SLC35A2', 'B4GALT1', 'B3GNT2', 'GALE', 'HSPA5', 'ST3GAL4']


def main(r, stem):
    r = Path(r)
    t = pd.read_csv(r / 'tsui/tsui_vs_jump_ko.csv')
    sel = pd.read_csv(r / 'periscope/periscope_channel_selectivity.csv')
    nm = pd.read_csv(r / 'periscope/periscope_named_genes.csv')
    jump = pd.read_csv(R / '17_channel_specificity/CRISPR_main_effects.csv').set_index('gene')
    apply()
    fig = plt.figure(figsize=(WIDTH, 4.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], width_ratios=[1, 1.25], hspace=0.42, wspace=0.34,
                          left=0.155, right=0.9, top=0.95, bottom=0.07)
    # a
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(CHJ))
    d = t[t.scope == 'all genes'].set_index('channel')
    ax.bar(x - 0.2, [d.loc[c, 'share_nonhits_up'] * 100 for c, _ in CHJ], 0.38, color='#cfcdc6', label='Other genes')
    ax.bar(x + 0.2, [d.loc[c, 'share_hits_up'] * 100 for c, _ in CHJ], 0.38, color='#1f66c2',
           label=f'Tsui hits (n = {int(d.tsui_hits.iloc[0])})')
    for i, (c, _) in enumerate(CHJ):
        p = d.loc[c, 'fisher_p']
        lab = 'p < 0.0001' if p < 1e-4 else f'p = {p:.2g}'
        ax.text(i + 0.2, d.loc[c, 'share_hits_up'] * 100 + 0.8, lab, ha='center', va='bottom', fontsize=FS_SMALL,
                rotation=90 if c != 'ConA' else 0)
    ax.set_xticks(x, [l for _, l in CHJ])
    ax.set_ylabel('Increased in JUMP-CP KO (%)')
    ax.set_ylim(0, 30)
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.24, 1.0))
    # b
    ax = fig.add_subplot(gs[0, 1])
    x = np.arange(len(CHP))
    for k, (cond, lab, col) in enumerate(COND):
        dd = sel[sel.condition == cond].set_index('channel')
        for i, (c, _) in enumerate(CHP):
            v = dd.loc[c]
            xi = i + (k - 1) * 0.24
            ax.plot([xi, xi], [v.ci_lo * 100, v.ci_hi * 100], color=col, lw=LW_ERR)
            ax.scatter(xi, v.share_above * 100, s=16, color=col if v.binom_p < 0.05 else 'white', edgecolor=col,
                       lw=LW_EDGE, zorder=3, label=lab if i == 0 else None)
    ax.axhline(5, color=INK2, lw=LW_REF, ls=':')
    ax.set_xticks(x, [l for _, l in CHP])
    ax.set_ylabel('Glycogenes above the 95th\npercentile of other genes (%)')
    ax.legend(frameon=False, loc='upper right', ncol=3, handletextpad=0.2, columnspacing=0.6)
    ax.set_ylim(0, 17)
    # c
    ax = fig.add_subplot(gs[1, :])
    cols, labels = [], []
    for ch, jch in (('WGA', 'WGA'), ('Phalloidin', None), ('ConA', 'ConA')):
        if jch:
            cols.append(('JUMP', jch)); labels.append('JUMP-CP (U2OS)')
        for cond, lab, _ in COND:
            cols.append((cond, ch)); labels.append(lab)
    m = np.full((len(GENES), len(cols)), np.nan)
    for i, g in enumerate(GENES):
        for j, (src, ch) in enumerate(cols):
            if src == 'JUMP':
                m[i, j] = jump.loc[g, f'{ch}_delta'] if g in jump.index else np.nan
            else:
                v = nm[(nm.condition == src) & (nm.gene == g)]
                m[i, j] = v[ch].iloc[0] if len(v) else np.nan
    im = ax.imshow(m.T, cmap='RdBu_r', vmin=-6, vmax=6, aspect='auto')
    for (i, j), v in np.ndenumerate(m):
        if np.isfinite(v):
            ax.text(i, j, f'{v:.1f}', ha='center', va='center', fontsize=5.5,
                    color='white' if abs(v) > 4.5 else 'black')
    ax.set_xticks(range(len(GENES)), GENES, style='italic')
    ax.set_yticks(range(len(cols)), labels, fontsize=FS_SMALL)
    for b in (3.5, 6.5):
        ax.axhline(b, color='black', lw=LW_SEP)
    for y, lab in ((1.5, 'WGA'), (5.0, 'Phalloidin'), (8.5, 'ConA')):
        ax.text(len(GENES) - 0.4, y, lab, ha='left', va='center', fontsize=FS_SMALL, fontweight='bold')
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.1)
    cb.set_label('Effect (JUMP-CP: Glass\'s Δ; PERISCOPE: z)', fontsize=FS_SMALL)
    cb.ax.tick_params(labelsize=FS_SMALL, width=0.6)
    cb.outline.set_linewidth(0.6)
    for sp in ax.spines.values():
        sp.set_linewidth(LW_SEP)
    for a in fig.axes[:2]:
        despine(a)
    letter(fig, 0.005, 0.995, 'A')
    letter(fig, 0.47, 0.995, 'B')
    letter(fig, 0.005, 0.5, 'C')
    save(fig, stem)


if __name__ == '__main__':
    main(*sys.argv[1:3])
