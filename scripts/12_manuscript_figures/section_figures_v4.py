"""Figure 1 for publication (user, 2026-11-22) and the bar chart moved to Supplementary Figure S2.

Reads existing results only; no statistic is recomputed except the threshold-free distribution comparison.
Figure 1
  a  ECDF of |Glass's delta| for GlycoEnzOnto genes vs other genes, ORF/CRISPR x ConA/WGA (all-gene tests,
     scripts/16_all_genes: <ORF|CRISPR>_all_genes_fdr.csv, control reagents already dropped). Two-sample
     Kolmogorov-Smirnov test (two-sided) per panel, written to figure1a_ks.csv.
  b  Gene-level strip plot of the glycogenes' effects (scripts/13_jump vs_other_genes_infection_filtered
     <ORF|CRISPR>_glycogene_calls.csv, FDR over glycogenes). Symmetric-log x axis (linear within +-1 SD) to spread the
     points near zero and keep extreme effects on scale. Grey small: p >= 0.05; tint: p < 0.05; full colour with black
     ring and bold name: FDR < 0.05. Named: FDR < 0.05, genes of interest, |effect| >= 3 SD.
Supplementary Figure S2: share of genes with |effect| >= 1 SD (16_all_genes/glyco_vs_other.csv), as the former Fig. 1b.
Sized for a 180 mm wide page (7.1 in); fonts 6-8 pt at print size.
Usage: python section_figures_v4.py <results dir> <out dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import FixedLocator  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from section_figures_v3 import share_panel  # noqa: E402

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5,
                     'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
GLY, OTH, INK, INK2, GRID = '#1f66c2', '#9a9890', '#111111', '#55534f', '#e6e5e1'
UP, DOWN, UP_T, DOWN_T = '#1f66c2', '#d9541e', '#9fc1ec', '#f2b394'
CONDS = [('ORF', 'ConA'), ('ORF', 'WGA'), ('CRISPR', 'ConA'), ('CRISPR', 'WGA')]
NAMES = {'ORF': 'Overexpression', 'CRISPR': 'Knockout'}
INTEREST = ['MGAT1', 'GCK', 'POMK', 'B3GNT2', 'ST3GAL4', 'FUT9', 'B4GALNT2', 'A4GNT', 'SLC35A2']
LINTHRESH = 1.0
# hand-placed names where automatic placement collides (offset points: dx, dy)
LABEL_OFFSETS = {('Overexpression · WGA', 'FUT9'): (-4, -12), ('Overexpression · WGA', 'ST6GAL1'): (-10, 9),
                 ('Overexpression · WGA', 'B4GALNT2'): (14, 9), ('Knockout · WGA', 'ST3GAL4'): (-6, -13),
                 ('Knockout · WGA', 'MGAT1'): (6, 10), ('Knockout · WGA', 'B3GNT2'): (-6, 10)}


def ecdf_panel(ax, r, lec, letter, title, show_y):
    d = r.dropna(subset=[f'{lec}_p'])
    g = d.glycogene
    x_g, x_o = np.sort(d.loc[g, f'{lec}_diff_sd'].abs()), np.sort(d.loc[~g, f'{lec}_diff_sd'].abs())
    for x, col, lw, lab in ((x_o, OTH, 1.3, f'Other genes (n = {len(x_o):,})'), (x_g, GLY, 1.5, f'Glycogenes (n = {len(x_g)})')):
        ax.step(np.r_[0, x], np.r_[0, np.arange(1, len(x) + 1) / len(x)], where='post', color=col, lw=lw, label=lab)
    ks = stats.ks_2samp(x_g, x_o)
    ax.set_xscale('symlog', linthresh=LINTHRESH, linscale=1.2)
    ax.set_xlim(0, 30)
    ax.xaxis.set_major_locator(FixedLocator([0, 0.5, 1, 3, 10, 30]))
    ax.set_xticklabels(['0', '0.5', '1', '3', '10', '30'])
    ax.set_ylim(0, 1.02)
    ax.set_xlabel('|Glass\'s Δ| (SD)')
    if show_y:
        ax.set_ylabel('Cumulative fraction of genes')
    else:
        ax.set_yticklabels([])
    ax.grid(color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    p = ks.pvalue
    ptxt = f'{p:.2g}' if p >= 0.001 else '{} × 10$^{{{}}}$'.format(*[int(float(v)) for v in f'{p:.0e}'.split('e')])
    ax.text(0.97, 0.06, f'KS D = {ks.statistic:.2f}\np = {ptxt}', transform=ax.transAxes, ha='right', va='bottom',
            fontsize=6.5, color=INK, fontweight='bold' if p < 0.05 else 'normal')
    ax.text(0.97, 0.30, title.replace(' · ', '\n'), transform=ax.transAxes, ha='right', va='bottom', fontsize=7,
            color=INK, linespacing=1.15)
    if letter:
        ax.text(-0.32 if show_y else -0.08, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight='bold',
                va='bottom')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    return [ks.statistic, p, len(x_g), len(x_o), float(np.median(x_g)), float(np.median(x_o))]


def strip_panel(ax, rows):
    rng = np.random.default_rng(3)
    for yi, (lab, g, lec) in enumerate(rows):
        x = g[f'{lec}_diff_sd'].astype(float)
        p = g[f'{lec}_p'].astype(float)
        fdr = g[f'{lec}_fdr'].astype(float)
        sig = fdr < 0.05
        nom = (p < 0.05) & ~sig
        y = yi + rng.uniform(-0.3, 0.3, len(x))
        rest = ~sig & ~nom
        ax.scatter(x[rest], y[rest.values], s=5, color='#c9c7c0', lw=0, alpha=0.8, zorder=2, rasterized=True)
        ax.scatter(x[nom], y[nom.values], s=11, c=np.where(x[nom] > 0, UP_T, DOWN_T), lw=0, zorder=3)
        ax.scatter(x[sig], y[sig.values], s=34, c=np.where(x[sig] > 0, UP, DOWN), edgecolor=INK, lw=0.7, zorder=4)
        # name a gene in a row only where it does something there
        name = g.index[sig | (g.index.isin(INTEREST) & ((p < 0.05) | (x.abs() >= 1.5))) | (x.abs() >= 3)]
        placed = []
        for gene in sorted(name, key=lambda k: -abs(x[k])):
            xi, yv = x[gene], y[g.index.get_loc(gene)]
            dx, dy = 0, 7
            for (px, py) in placed:  # simple de-overlap: alternate above/below when names are close
                if abs(np.arcsinh(px) - np.arcsinh(xi)) < 0.35 and abs(py - (yv + 0.001)) < 0.5:
                    dy = -10 if dy > 0 else 7
            dx, dy = LABEL_OFFSETS.get((lab, gene), (dx, dy))
            placed.append((xi, yv))
            ax.annotate(gene, (xi, yv), xytext=(dx, dy), textcoords='offset points', ha='center', fontsize=6.3,
                        fontstyle='italic', fontweight='bold' if sig[gene] else 'normal', color=INK, zorder=5,
                        arrowprops=dict(arrowstyle='-', lw=0.4, color=INK2, shrinkA=0, shrinkB=2))
    ax.axvline(0, color=INK2, lw=0.6, zorder=1)
    for v in (-1, 1):
        ax.axvline(v, color=GRID, lw=0.6, ls=(0, (2, 2)), zorder=1)
    ax.set_xscale('symlog', linthresh=LINTHRESH, linscale=1.2)
    ax.set_xlim(-8, 32)
    ax.xaxis.set_major_locator(FixedLocator([-5, -3, -1, 0, 1, 3, 10, 30]))
    ax.set_xticklabels(['−5', '−3', '−1', '0', '1', '3', '10', '30'])
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_xlabel('Glass\'s Δ vs wells of other genes on the same plates (SD; linear within ±1)')
    ax.grid(axis='x', color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    h = [plt.Line2D([], [], marker='o', ls='', ms=5.5, mfc=UP, mec=INK, mew=0.7),
         plt.Line2D([], [], marker='o', ls='', ms=5.5, mfc=DOWN, mec=INK, mew=0.7),
         plt.Line2D([], [], marker='o', ls='', ms=3.6, mfc=UP_T, mec='none'),
         plt.Line2D([], [], marker='o', ls='', ms=3.6, mfc=DOWN_T, mec='none'),
         plt.Line2D([], [], marker='o', ls='', ms=2.6, mfc='#c9c7c0', mec='none')]
    ax.legend(h, ['Brighter, FDR < 0.05', 'Darker, FDR < 0.05', 'Brighter, p < 0.05', 'Darker, p < 0.05', 'p ≥ 0.05'],
              frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=5, handletextpad=0.2,
              columnspacing=1.0)
    ax.text(-0.16, 1.02, 'b', transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom')


def main(r, out):
    r, out = Path(r), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    allg = {t: pd.read_csv(r / '16_all_genes' / f'{t}_all_genes_fdr.csv', index_col=0) for t in ('ORF', 'CRISPR')}
    d = r / '13_jump/vs_other_genes_infection_filtered'
    calls = {t: pd.read_csv(d / f'{t}_glycogene_calls.csv', index_col=0) for t in ('ORF', 'CRISPR')}

    fig = plt.figure(figsize=(7.1, 6.0))
    top = fig.add_gridspec(1, 4, left=0.1, right=0.98, top=0.93, bottom=0.68, wspace=0.12)
    ks_rows = []
    for j, (t, lec) in enumerate(CONDS):
        ax = fig.add_subplot(top[0, j])
        res = ecdf_panel(ax, allg[t], lec, 'a' if j == 0 else '', f'{NAMES[t]} · {lec}', j == 0)
        ks_rows.append([t, lec] + res)
        if j == 0:
            hs, ls = ax.get_legend_handles_labels()
    fig.legend([hs[1], hs[0]], ['Glycogenes', 'Other genes'], frameon=False, loc='upper right', ncol=2,
               bbox_to_anchor=(0.98, 0.995), handlelength=1.4)
    bot = fig.add_gridspec(1, 1, left=0.17, right=0.98, top=0.56, bottom=0.12)
    strip_panel(fig.add_subplot(bot[0, 0]), [(f'{NAMES[t]} · {lec}', calls[t], lec) for t, lec in CONDS])
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'Figure1.{ext}', dpi=300)
    plt.close(fig)
    pd.DataFrame(ks_rows, columns=['perturbation', 'lectin', 'ks_D', 'ks_p', 'n_glycogenes', 'n_other',
                                   'median_abs_delta_glyco', 'median_abs_delta_other']).to_csv(out / 'figure1a_ks.csv',
                                                                                               index=False)

    comp = pd.read_csv(r / '16_all_genes/glyco_vs_other.csv')
    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    share_panel(ax, comp)
    ax.set_title('', loc='left')
    for txt in ax.texts:
        txt.set_fontsize(6.5)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'FigureS2.{ext}', dpi=300)
    plt.close(fig)
    print(pd.read_csv(out / 'figure1a_ks.csv').round(4).to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:3])
