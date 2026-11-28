"""Figure 1 reorganised (user, 2026-11-22): "Glycogene perturbations preferentially affect lectin-containing channels".
  a  ECDF of |Glass's delta|, glycogenes vs other genes, 5 channels x ORF/CRISPR (former Fig. 4a)
  b  share of glycogenes above the other genes' 95th percentile per channel (former Fig. 4c)
  c  gene-level strip plot of the glycogenes' ConA / WGA effects (former Fig. 1b)
Supplementary Figure S3: KS D vs random gene sets (former Fig. 4b); S4: strong responders in all channels (former 4d).
Reads existing results only (scripts/17_channel_specificity outputs, scripts/13_jump vs_other_genes_infection_filtered).
Usage: python section_figures_v5.py <results dir> <out dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from matplotlib.ticker import FixedLocator, PercentFormatter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from section_figures_v4 import CONDS, NAMES, strip_panel  # noqa: E402

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5,
                     'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'pdf.fonttype': 42})
CH = ['ConA', 'WGA', 'DNA', 'RNA', 'Mito']
CH_LABEL = {'ConA': 'ER (ConA)', 'WGA': 'AGP (WGA)', 'DNA': 'DNA', 'RNA': 'RNA', 'Mito': 'Mito'}
PT = {'ORF': 'Overexpression', 'CRISPR': 'Knockout'}
GLY, OTH, INK, INK2, GRID = '#1f66c2', '#9a9890', '#111111', '#55534f', '#e6e5e1'
PT_COL = {'ORF': '#1f66c2', 'CRISPR': '#d9541e'}


def ptxt(p):
    return f'{p:.2g}' if p >= 0.001 else '{} × 10$^{{{}}}$'.format(*[int(float(v)) for v in f'{p:.0e}'.split('e')])


def ecdf_grid(fig, gs, eff, summ):
    for i, t in enumerate(('ORF', 'CRISPR')):
        r = eff[t]
        for j, ch in enumerate(CH):
            ax = fig.add_subplot(gs[i, j])
            d = r.dropna(subset=[f'{ch}_p'])
            for m, col, lw in ((~d.glycogene, OTH, 1.1), (d.glycogene, GLY, 1.3)):
                x = np.sort(d.loc[m, f'{ch}_delta'].abs())
                ax.step(np.r_[0, x], np.r_[0, np.arange(1, len(x) + 1) / len(x)], where='post', color=col, lw=lw)
            s = summ[(summ.perturbation == t) & (summ.channel == ch)].iloc[0]
            ax.set_xscale('symlog', linthresh=1, linscale=1.2)
            ax.set_xlim(0, 30)
            ax.xaxis.set_major_locator(FixedLocator([0, 0.5, 1, 3, 10, 30]))
            ax.set_xticklabels(['0', '.5', '1', '3', '10', '30'])
            ax.set_ylim(0, 1.02)
            ax.grid(color=GRID, lw=0.4)
            ax.set_axisbelow(True)
            ax.text(0.96, 0.05, f'D = {s.ks_D:.2f}\np = {ptxt(s.ks_p)}', transform=ax.transAxes, ha='right',
                    va='bottom', fontsize=6, fontweight='bold' if s.ks_p < 0.05 else 'normal')
            if i == 0:
                ax.text(0.5, 1.06, CH_LABEL[ch], transform=ax.transAxes, ha='center', va='bottom', fontsize=7)
            if j == 0:
                ax.set_ylabel(f'{PT[t]}\nCumulative fraction')
            else:
                ax.set_yticklabels([])
            if i == 1:
                ax.set_xlabel('|Glass\'s Δ| (SD)')
            else:
                ax.set_xticklabels([])
            for sp in ('top', 'right'):
                ax.spines[sp].set_visible(False)


def dot_panel(ax, summ, col, sigcol, ylab, lo=None, hi=None, ref=None, null=False):
    for dx, t in ((-0.14, 'ORF'), (0.14, 'CRISPR')):
        s = summ[summ.perturbation == t].set_index('channel').loc[CH]
        x = np.arange(len(CH)) + dx
        sig = s[sigcol] < 0.05
        if null:
            for xi, (m, q) in zip(x, zip(s.D_null_median, s.D_null_95)):
                ax.plot([xi - 0.1, xi + 0.1], [q, q], color=PT_COL[t], lw=0.6, ls=(0, (2, 1.5)))
                ax.plot([xi, xi], [m, q], color=PT_COL[t], lw=0.5, alpha=0.5)
        if lo:
            ax.errorbar(x, s[col], yerr=[s[col] - s[lo], s[hi] - s[col]], fmt='none', ecolor=PT_COL[t], lw=0.7)  # common style (figure_style.py)
        ax.scatter(x[~sig.values], s[col][~sig.values], s=18, facecolor='white', edgecolor=PT_COL[t], lw=0.7, zorder=3)
        ax.scatter(x[sig.values], s[col][sig.values], s=18, color=PT_COL[t], zorder=3)
    if ref is not None:
        ax.axhline(ref, color=INK2, lw=0.6, ls=(0, (2, 2)))
        ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.set_xticks(range(len(CH)), [CH_LABEL[c] for c in CH])
    ax.set_ylabel(ylab)
    ax.set_ylim(0, None)
    ax.grid(axis='y', color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    h = [plt.Line2D([], [], marker='o', ls='', color=PT_COL['ORF'], ms=4.5),
         plt.Line2D([], [], marker='o', ls='', color=PT_COL['CRISPR'], ms=4.5),
         plt.Line2D([], [], marker='o', ls='', mfc='white', mec=INK2, ms=4.5)]
    ax.legend(h, ['Overexpression', 'Knockout', 'p ≥ 0.05 (open)'], frameon=False, loc='upper right', ncol=3)


def heatmap(fig, st, out):
    fig2 = plt.figure(figsize=(7.1, 3.0))
    n_o, n_c = (st.perturbation == 'ORF').sum(), (st.perturbation == 'CRISPR').sum()
    gs = fig2.add_gridspec(1, 2, left=0.11, right=0.9, top=0.86, bottom=0.25, wspace=0.35, width_ratios=[n_o, max(n_c, 3)])
    norm = TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5)
    for k, t in enumerate(('ORF', 'CRISPR')):
        s = st[st.perturbation == t].sort_values(['lead', 'lead_abs'], ascending=[True, False])
        ax = fig2.add_subplot(gs[0, k])
        m = s[[f'{c}_delta' for c in CH]].to_numpy().T
        im = ax.imshow(np.clip(m, -5, 5), cmap='RdBu_r', norm=norm, aspect='auto')
        for (yi, xi), v in np.ndenumerate(m):
            if s[f'{CH[yi]}_p'].iloc[xi] < 0.05:
                ax.text(xi, yi, '*', ha='center', va='center', fontsize=8, color='white' if abs(v) > 3 else INK)
        ax.set_yticks(range(len(CH)), [CH_LABEL[c] for c in CH])
        ax.set_xticks(range(len(s)), s.index, rotation=90, fontstyle='italic', fontsize=6.3)
        ax.tick_params(length=0)
        ax.set_xlabel(PT[t], labelpad=4)
        ax.xaxis.set_label_position('top')
    cax = fig2.add_axes([0.92, 0.3, 0.012, 0.5])
    cb = fig2.colorbar(im, cax=cax)
    cb.set_label('Glass\'s Δ (SD)', fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    for ext in ('png', 'pdf'):
        fig2.savefig(out / f'FigureS4.{ext}', dpi=300)
    plt.close(fig2)


def main(r, out):
    r, out = Path(r), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    c = r / '17_channel_specificity'
    eff = {t: pd.read_csv(c / f'{t}_main_effects.csv', index_col=0) for t in ('ORF', 'CRISPR')}
    summ = pd.read_csv(c / 'channel_summary_main.csv')
    st = pd.read_csv(c / 'strong_responders_main.csv', index_col=0)
    d = r / '13_jump/vs_other_genes_infection_filtered'
    calls = {t: pd.read_csv(d / f'{t}_glycogene_calls.csv', index_col=0) for t in ('ORF', 'CRISPR')}

    fig = plt.figure(figsize=(7.1, 9.0))
    ga = fig.add_gridspec(2, 5, left=0.1, right=0.98, top=0.95, bottom=0.68, wspace=0.12, hspace=0.18)
    ecdf_grid(fig, ga, eff, summ)
    fig.legend([plt.Line2D([], [], color=GLY, lw=1.3), plt.Line2D([], [], color=OTH, lw=1.1)],
               ['Glycogenes', 'Other genes'], loc='upper right', bbox_to_anchor=(0.98, 0.995), ncol=2, frameon=False)
    gb = fig.add_gridspec(1, 1, left=0.1, right=0.98, top=0.6, bottom=0.47)
    dot_panel(fig.add_subplot(gb[0, 0]), summ, 'share_above', 'binom_p',
              'Glycogenes above the 95th\npercentile of other genes', 'share_lo', 'share_hi', ref=0.05)
    gc = fig.add_gridspec(1, 1, left=0.17, right=0.98, top=0.39, bottom=0.07)
    ax = fig.add_subplot(gc[0, 0])
    strip_panel(ax, [(f'{NAMES[t]} · {lec}', calls[t], lec) for t, lec in CONDS])
    for t in ax.texts:  # strip_panel writes its own panel letter; Figure 1 places letters below
        if t.get_text() == 'b':
            t.set_visible(False)
    for letter, (x, y) in {'a': (0.01, 0.985), 'b': (0.01, 0.625), 'c': (0.01, 0.415)}.items():
        fig.text(x, y, letter, fontsize=9, fontweight='bold', va='top')
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'Figure1.{ext}', dpi=300)
    plt.close(fig)

    fig3, ax = plt.subplots(figsize=(4.6, 2.4))
    dot_panel(ax, summ, 'ks_D', 'D_perm_p', 'KS D (glycogenes vs other genes)', null=True)
    ax.legend(ax.get_legend().legend_handles + [plt.Line2D([], [], color=INK2, lw=0.8, ls=(0, (2, 1.5)))],
              ['Overexpression', 'Knockout', 'p ≥ 0.05 (open)', '95th percentile, random sets'], frameon=False,
              loc='upper right', ncol=2)
    fig3.tight_layout()
    for ext in ('png', 'pdf'):
        fig3.savefig(out / f'FigureS3.{ext}', dpi=300)
    plt.close(fig3)
    heatmap(fig, st, out)


if __name__ == '__main__':
    main(*sys.argv[1:3])
