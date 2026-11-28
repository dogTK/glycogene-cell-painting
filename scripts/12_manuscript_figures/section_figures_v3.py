"""Manuscript figures after merging the glycogene-up and glycogene-down sections into one (user, 2026-11-21).

Figure 1  Glycogene overexpression (ORF) and knockout (CRISPR) in JUMP-CP: a. each glycogene's effect in the four
          combinations (colour = direction and size, marker size = -log10 FDR, ring = FDR < 0.05), b. lysosome test
Figure 2  Share of genes with p < 0.05: glycogenes vs other genes (section_figures_v2.figure_rates)
Figure 3  Drugs (LINCS) (section_figures.figure4)
'Other genes' = the wells on the same plates in which a different gene was overexpressed (ORF) or knocked out (CRISPR),
negative controls excluded (scripts/13_jump/vs_other_genes_test.py).
Usage: python section_figures_v3.py <results dir> <out dir> [<vs_other_genes dir> <ORF gene scores csv> [<all-gene comparison csv>]]
(manuscript since 2026-11-22: 13_jump/vs_other_genes_infection_filtered 13_jump/infection_filtered/ORF_gene_scores.csv)
"""
import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scipy import stats  # noqa: E402

from section_figures import DOWN, GRID, INK, INK2, LYSO, UP, figure4, lyso, save, title  # noqa: E402
from section_figures_v2 import figure_rates  # noqa: E402


LABEL_SD = 2.8
XLIM = (-6, 12)   # effects beyond the limit are drawn at the edge with their value
CMAP = LinearSegmentedColormap.from_list('div', [DOWN, '#f3c9b5', '#cfcdc6', '#bcd4f0', UP])
NORM = TwoSlopeNorm(vmin=-4, vcenter=0, vmax=4)
OFFSETS = {('Overexpression · WGA', 'POMK'): -12, ('Overexpression · WGA', 'FUT9'): -13, ('Overexpression · WGA', 'A4GNT'): 6,
           ('Overexpression · WGA', 'ST8SIA1'): 6, ('Overexpression · WGA', 'B4GALNT2'): -13,
           ('Overexpression · WGA', 'GCK'): 8, ('Overexpression · WGA', 'HS3ST5'): -15}


def strips(ax, rows):
    """One row per perturbation x lectin: each glycogene's effect, colour = direction and size of the effect,
    marker size = -log10 FDR, black ring = FDR < 0.05."""
    rng = np.random.default_rng(1)
    for yi, (lab, g, lec) in enumerate(rows):
        x = g[f'{lec}_diff_sd'].astype(float)
        fdr = g[f'{lec}_fdr'].astype(float)
        sig = g[f'{lec}_call'] != 'no change'
        xc = x.clip(*XLIM)
        y = yi + rng.uniform(-0.28, 0.28, len(x))
        size = 10 + 60 * (-np.log10(fdr.clip(lower=1e-4))).clip(upper=3)
        ax.scatter(xc[~sig], y[~sig], s=size[~sig], c=x[~sig], cmap=CMAP, norm=NORM, lw=0, alpha=0.55, zorder=2)
        ax.scatter(xc[sig], y[sig], s=size[sig] + 20, c=np.where(x[sig] > 0, UP, DOWN), edgecolor=INK, lw=1.2,
                   zorder=3)
        lab_genes = g.index[sig | (x.abs() >= LABEL_SD)]
        for gene in lab_genes:
            xi, yv = xc[gene], y[g.index.get_loc(gene)]
            txt = gene if abs(x[gene]) <= XLIM[1] else f'{gene} ({x[gene]:+.0f})'
            dy = OFFSETS.get((lab, gene), 4)
            edge = abs(x[gene]) > XLIM[1]
            ax.annotate(txt, (xi, yv), xytext=(-7 if edge else 6, dy), textcoords='offset points', fontsize=7.5,
                        ha='right' if edge else 'left', fontstyle='italic',
                        fontweight='bold' if sig[gene] else 'normal', color=INK, zorder=4)
        n_up, n_dn = (g[f'{lec}_call'] == 'stronger').sum(), (g[f'{lec}_call'] == 'weaker').sum()
        ax.text(1.06, yi, f'{n_up} stronger\n{n_dn} weaker\n(of {len(g)})', va='center', fontsize=7.5,
                color=INK2, transform=ax.get_yaxis_transform())
    ax.axvline(0, color=INK2, lw=0.8, zorder=1)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlim(XLIM[0] - 0.3, XLIM[1] + 0.3)
    ax.set_xlabel("Staining vs wells of other genes on the same plates (SD)")
    ax.grid(axis='x', color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    h = [plt.Line2D([], [], marker='o', ls='', ms=7, mfc=UP, mec=INK, mew=1.1),
         plt.Line2D([], [], marker='o', ls='', ms=7, mfc=DOWN, mec=INK, mew=1.1),
         plt.Line2D([], [], marker='o', ls='', ms=4, mfc='#cfcdc6', mec='none'),
         plt.Line2D([], [], marker='o', ls='', ms=9, mfc='#cfcdc6', mec='none')]
    ax.legend(h, ['Stronger (FDR < 0.05)', 'Weaker (FDR < 0.05)', 'FDR = 1', 'FDR ≤ 0.001'], frameon=False,
              fontsize=7.5, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=4)


def share_panel(ax, comp):
    """Share of genes with |effect| >= 1 SD, GlycoEnzOnto genes vs other genes, four conditions (Fisher p)."""
    c = comp[comp.criterion == '|effect| >= 1 SD'].set_index(['perturbation', 'lectin'])
    conds = [('ORF', 'ConA'), ('ORF', 'WGA'), ('CRISPR', 'ConA'), ('CRISPR', 'WGA')]
    names = {'ORF': 'Overexpression', 'CRISPR': 'Knockout'}
    for i, k in enumerate(conds):
        r = c.loc[k]
        for dx, v, col, lab in ((-0.19, r.glyco_share, UP, 'Glycogenes'), (0.19, r.other_share, '#b4b2ab', 'Other genes')):
            ax.bar(i + dx, v * 100, width=0.36, color=col, label=lab if i == 0 else None)
            ax.text(i + dx, v * 100 + 0.25, f'{v * 100:.1f}', ha='center', fontsize=7.5, color=INK)
        p = r.fisher_p
        ax.text(i, max(r.glyco_share, r.other_share) * 100 + 1.4, 'p = {} × 10$^{{{}}}$'.format(*[int(float(x)) for x in f'{p:.0e}'.split('e')]) if p < 0.001 else f'p = {p:.2g}',
                ha='center', fontsize=8, color=INK, fontweight='bold' if p < 0.05 else 'normal')
    ax.set_xticks(range(4), [f'{names[t]}\n{lec}'.replace('Overexpression', 'Over-\nexpression') for t, lec in conds],
                  fontsize=8)
    ax.set_ylabel('Genes with effect ≥ 1 SD in either direction (%)')
    ax.set_ylim(0, 16)
    ax.legend(frameon=False, fontsize=8, loc='upper right')
    ax.grid(axis='y', color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.set_title('b', loc='left', fontsize=11, fontweight='bold')


def figure1(d, orf_gene_level, out, with_lyso=False, comp=None):
    """with_lyso=False since 2026-11-22 (user: the lysosome test is too early for Figure 1).
    comp: all-gene comparison table (scripts/16_all_genes glyco_vs_other.csv) -> panel b (user, 2026-11-22)."""
    up = pd.read_csv(d / 'ORF_glycogene_calls.csv', index_col=0)
    down = pd.read_csv(d / 'CRISPR_glycogene_calls.csv', index_col=0)
    if comp is not None:
        fig = plt.figure(figsize=(17, 5.8))
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.36], wspace=0.32)
        ax = fig.add_subplot(gs[0])
        strips(ax, [('Overexpression · ConA', up, 'ConA'), ('Overexpression · WGA', up, 'WGA'),
                    ('Knockout · ConA', down, 'ConA'), ('Knockout · WGA', down, 'WGA')])
        ax.set_title('a', loc='left', fontsize=11, fontweight='bold')
        share_panel(fig.add_subplot(gs[1]), comp)
        save(fig, out, 'Figure1')
        return
    fig = plt.figure(figsize=(15, 5.8) if with_lyso else (12, 5.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.32], wspace=0.28) if with_lyso else fig.add_gridspec(1, 1)
    strips(fig.add_subplot(gs[0]), [('Overexpression · ConA', up, 'ConA'), ('Overexpression · WGA', up, 'WGA'),
                                    ('Knockout · ConA', down, 'ConA'), ('Knockout · WGA', down, 'WGA')])
    if not with_lyso:
        save(fig, out, 'Figure1')
        return
    ax = fig.add_subplot(gs[1])
    lyso(ax, orf_gene_level)
    lyso_genes = orf_gene_level.index.isin(LYSO)
    x = orf_gene_level.loc[lyso_genes, 'ConA_adj'].dropna()
    y = orf_gene_level.loc[~lyso_genes, 'ConA_adj'].dropna()
    tt = stats.ttest_ind(x, y, equal_var=False, alternative='greater')
    se = np.sqrt(x.var() / len(x) + y.var() / len(y))
    dfree = se ** 4 / ((x.var() / len(x)) ** 2 / (len(x) - 1) + (y.var() / len(y)) ** 2 / (len(y) - 1))
    dd, half = x.mean() - y.mean(), stats.t.ppf(0.975, dfree) * se
    ax.set_xticks([0, 1], [f'Lysosome genes\n({len(x)})', f'Other genes\n({len(y):,})'])
    ax.set_title(f'b  Pre-specified test (overexpression, ConA)\n{dd:+.2f} (95% CI {dd - half:.2f}–{dd + half:.2f}), '
                 f'one-sided p = {tt.pvalue:.3f}', loc='left', fontsize=9.5, color=INK)
    save(fig, out, 'Figure1')


def main(r, out, vs_dir='13_jump/vs_other_genes', scores='13_jump/ORF_gene_scores.csv', comp=None):
    r, out = Path(r), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    d = r / vs_dir
    figure1(d, pd.read_csv(r / scores, index_col=0), out, comp=pd.read_csv(r / comp) if comp else None)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        figure_rates(d, tmp)
        figure4(r, tmp)
        for src, dst in (('Figure3', 'Figure2'), ('Figure4', 'Figure3')):
            for ext in ('png', 'pdf'):
                shutil.copy(tmp / f'{src}.{ext}', out / f'{dst}.{ext}')
    print('saved Figure1-3 to', out)


if __name__ == '__main__':
    main(*sys.argv[1:])
