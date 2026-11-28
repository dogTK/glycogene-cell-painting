"""One figure per Results section, with 'stronger / weaker' defined by the test against other genes on the same plates
(scripts/13_jump/vs_other_genes_test.py: well-position corrected, well-level adjusted, Welch's t-test, BH over glycogenes).

Figure 1  Glycogene up (ORF): a. counts, b. ConA volcano, c. WGA volcano, d. pre-specified lysosome test
Figure 2  Glycogene down (CRISPR): a. counts, b. ConA volcano, c. WGA volcano
Figure 3  Share of genes with p < 0.05: glycogenes vs other genes
Figure 4  Drugs (LINCS): unchanged (section_figures.py figure4)
Usage: python section_figures_v2.py <results dir> <out dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from section_figures import BG, GRAY, GRID, INK, INK2, NONE, UP, DOWN, figure4, lyso, save, title  # noqa: E402

FDR = 0.05
LABEL_MIN_SD = 3.0


def counts(ax, g, t_label):
    rows = []
    for lec in ('ConA', 'WGA'):
        c = g[f'{lec}_call'].value_counts()
        rows.append((lec, c.get('stronger', 0), c.get('weaker', 0), c.get('no change', 0), len(g)))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, u, d, n, tot) in zip(y, rows):
        ax.barh(yi, u / tot * 100, color=UP, height=0.55)
        ax.barh(yi, d / tot * 100, left=u / tot * 100, color=DOWN, height=0.55)
        ax.barh(yi, n / tot * 100, left=(u + d) / tot * 100, color=NONE, height=0.55)
        ax.text(98, yi, f'{u} stronger · {d} weaker · {n} no change (of {tot})', va='center', ha='right',
                fontsize=7.5, color=INK)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_xlim(0, 100)
    ax.set_xlabel('Glycogenes (%)')
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (UP, DOWN, NONE)]
    ax.legend(handles, ['Stronger (FDR < 0.05)', 'Weaker (FDR < 0.05)', 'No change'], frameon=False, fontsize=7.5,
              loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=3)
    title(ax, t_label)


def volcano(ax, g, lec, t_label):
    x = g[f'{lec}_diff_sd']
    y = -np.log10(g[f'{lec}_fdr'].clip(lower=1e-12))
    col = np.where(g[f'{lec}_call'] == 'stronger', UP, np.where(g[f'{lec}_call'] == 'weaker', DOWN, GRAY))
    ax.scatter(x, y, s=14, c=col, lw=0, alpha=0.9)
    ax.axhline(-np.log10(FDR), color=INK2, lw=0.6, ls='--')
    ax.axvline(0, color=INK2, lw=0.6)
    lab = g[(g[f'{lec}_call'] != 'no change') | (x.abs() >= LABEL_MIN_SD)]
    for gene in lab.index:
        ax.annotate(gene, (x[gene], y[gene]), xytext=(4, 2), textcoords='offset points', fontsize=7,
                    fontstyle='italic', color=INK)
    ax.set_xlabel(f'{lec} staining vs other genes on the same plates (SD)')
    ax.set_ylabel('−log10 FDR (over glycogenes)')
    ax.grid(color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    title(ax, t_label)


def figure_up(d, orf_gene_level, out):
    g = pd.read_csv(d / 'ORF_glycogene_calls.csv', index_col=0)
    fig = plt.figure(figsize=(15, 8.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.55, 1], width_ratios=[1, 1, 0.65], hspace=0.6, wspace=0.35)
    counts(fig.add_subplot(gs[0, :2]), g, 'a  Glycogene up (ORF, U2OS): did the staining change?')
    volcano(fig.add_subplot(gs[1, 0]), g, 'ConA', 'b  ConA')
    volcano(fig.add_subplot(gs[1, 1]), g, 'WGA', 'c  WGA')
    ax = fig.add_subplot(gs[1, 2])
    lyso(ax, orf_gene_level)
    ax.set_title('d  Pre-specified test: +1.03\n(95% CI 0.26–1.79), one-sided p = 0.006', loc='left',
                 fontsize=9.5, color=INK)
    save(fig, out, 'Figure1')


def figure_down(d, out):
    g = pd.read_csv(d / 'CRISPR_glycogene_calls.csv', index_col=0)
    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.55, 1], hspace=0.6, wspace=0.3)
    counts(fig.add_subplot(gs[0, :]), g, 'a  Glycogene down (CRISPR, U2OS): did the staining change?')
    volcano(fig.add_subplot(gs[1, 0]), g, 'ConA', 'b  ConA')
    volcano(fig.add_subplot(gs[1, 1]), g, 'WGA', 'c  WGA')
    save(fig, out, 'Figure2')


def figure_rates(d, out):
    fig, ax = plt.subplots(figsize=(7, 3.6))
    groups, x = [], 0
    for t, tl in (('ORF', 'Up'), ('CRISPR', 'Down')):
        r = pd.read_csv(d / f'{t}_tests_all_genes.csv', index_col=0)
        for lec in ('ConA', 'WGA'):
            p = r[f'{lec}_p']
            ok = p.notna()
            sig = (p < 0.05) & ok
            a, b = (sig & r.glycogene).sum(), (~sig & ok & r.glycogene).sum()
            c, e = (sig & ~r.glycogene).sum(), (~sig & ok & ~r.glycogene).sum()
            pf = fisher_exact([[a, b], [c, e]])[1]
            r1, r2 = a / (a + b) * 100, c / (c + e) * 100
            ax.bar(x - 0.18, r1, width=0.34, color=UP, label='Glycogenes' if x == 0 else None)
            ax.bar(x + 0.18, r2, width=0.34, color=GRAY, label='Other genes' if x == 0 else None)
            ax.text(x, max(r1, r2) + 0.6, f'p = {pf:.2g}', ha='center', fontsize=7.5, color=INK2)
            groups.append(f'{tl}\n{lec}')
            x += 1
    ax.axhline(5, color=INK2, lw=0.6, ls=':')
    ax.text(-0.45, 5.4, 'expected by chance (5%)', ha='left', fontsize=7, color=INK2)
    ax.set_xticks(range(len(groups)), groups)
    ax.set_ylabel('Genes with p < 0.05 (%)\nvs other genes on the same plates')
    ax.set_ylim(0, 23)
    ax.legend(frameon=False, fontsize=7.5, loc='upper right')
    ax.grid(axis='y', color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    title(ax, 'Glycogenes change the staining about as often as other genes')
    save(fig, out, 'Figure3')


def main(r, out):
    r, out = Path(r), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    d = r / '13_jump/vs_other_genes'
    figure_up(d, pd.read_csv(r / '13_jump/ORF_gene_scores.csv', index_col=0), out)
    figure_down(d, out)
    figure_rates(d, out)
    figure4(r, out)
    print('saved Figure1-4 to', out)


if __name__ == '__main__':
    main(*sys.argv[1:3])
