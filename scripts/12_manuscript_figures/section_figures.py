"""One figure per Results section.

Figure 1  Glycogene up (ORF): a. stronger / weaker / no change counts for ConA and WGA, b. adjusted ConA vs WGA,
          c. pre-specified test (lysosome / M6P ORFs vs other ORFs, adjusted ConA)
Figure 2  Glycogene down (CRISPR): a. counts, b. adjusted ConA vs WGA
Figure 3  Share of genes whose staining changed (|adjusted z| >= 2): glycogenes vs other genes (Fisher's exact test)
Figure 4  Drugs (LINCS): a. share of compounds whose staining changed in the predicted direction, predicted vs other,
          before and after adjustment; b. share of toxic compounds among ConA-bright and other compounds
Usage: python section_figures.py <results dir> <out dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

UP, DOWN, NONE, GRAY = '#2a78d6', '#eb6834', '#d9d7d0', '#b4b2ab'
INK, INK2, GRID, BG = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
THR = 2.0
LYSO = ['TFEB', 'GNPTAB', 'GNPTG', 'NAGPA', 'M6PR', 'HEXB', 'NEU1', 'GLB1', 'GNS', 'FUCA2', 'LAMP1', 'LAMP2', 'CTSD']
LABELS = {'CRISPR': {'MGAT1': (4, 3), 'MGAT2': (4, -9), 'MAN2A1': (4, 3), 'OGT': (4, 3), 'B3GNT2': (4, 3)},
          'ORF': {'B3GNT2': (5, -2), 'GNPTAB': (-42, 4), 'RENBP': (5, 5), 'ST3GAL4': (5, -6), 'SLC35A1': (-48, -8),
                  'GLB1': (5, 2)}}
plt.rcParams.update({'font.family': 'Helvetica', 'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 0.8, 'figure.facecolor': BG,
                     'axes.facecolor': BG, 'savefig.facecolor': BG})


def title(ax, t):
    ax.set_title(t, loc='left', fontsize=9.5, color=INK)


def save(fig, out, name):
    fig.savefig(out / f'{name}.png', dpi=200, bbox_inches='tight')
    fig.savefig(out / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


def counts(ax, g, t_label):
    g = g[g.glycogene]
    rows = []
    for lec in ('ConA', 'WGA'):
        z = g[f'{lec}_adj'].dropna()
        rows.append((lec, (z >= THR).sum(), (z <= -THR).sum(), (z.abs() < THR).sum(), len(z)))
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
    ax.legend(handles, ['Stronger (z ≥ 2)', 'Weaker (z ≤ −2)', 'No change'], frameon=False, fontsize=7.5,
              loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=3)
    title(ax, t_label)


def scatter(ax, g, t, t_label, lim=None):
    ax.scatter(g.ConA_adj, g.WGA_adj, s=3, color=GRAY, lw=0, alpha=0.45, rasterized=True, label='All genes')
    gg = g[g.glycogene]
    ax.scatter(gg.ConA_adj, gg.WGA_adj, s=12, color=UP, lw=0, alpha=0.85, label='Glycogenes')
    for gene, off in LABELS[t].items():
        x, y = g.loc[gene, 'ConA_adj'], g.loc[gene, 'WGA_adj']
        ax.scatter([x], [y], s=26, facecolor='none', edgecolor=INK, lw=0.9)
        ax.annotate(gene, (x, y), xytext=off, textcoords='offset points', fontsize=7.5, fontstyle='italic', color=INK)
    for v in (-THR, THR):
        ax.axvline(v, color=INK2, lw=0.6, ls='--')
        ax.axhline(v, color=INK2, lw=0.6, ls='--')
    if lim:
        ax.set_xlim(*lim[0])
        ax.set_ylim(*lim[1])
    ax.set_xlabel('ConA staining, adjusted (z)')
    ax.set_ylabel('WGA staining, adjusted (z)')
    ax.grid(color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=7.5, loc='upper left', markerscale=1.5)
    title(ax, t_label)


def lyso(ax, orf):
    m = orf.index.isin(LYSO)
    x, y = orf.loc[m, 'ConA_adj'].dropna(), orf.loc[~m, 'ConA_adj'].dropna()
    rng = np.random.default_rng(0)
    ax.scatter(rng.normal(1, 0.06, len(y)), y, s=2, color=GRAY, lw=0, alpha=0.3, rasterized=True)
    ax.scatter(rng.normal(0, 0.06, len(x)), x, s=22, color=UP, lw=0)
    for xi, s in ((0, x), (1, y)):
        half = stats.t.ppf(0.975, len(s) - 1) * s.sem()
        ax.plot([xi + 0.22] * 2, [s.mean() - half, s.mean() + half], color=INK, lw=2)
        ax.plot(xi + 0.22, s.mean(), 'o', color=INK, ms=5)
    ax.set_xticks([0, 1], [f'Lysosome genes\nup ({len(x)})', f'Other genes\nup ({len(y):,})'])
    ax.set_ylim(-6, 6)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_ylabel('ConA staining, adjusted (z)')
    ax.grid(axis='y', color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    title(ax, 'c  Pre-specified test: +1.03\n(95% CI 0.26–1.79), one-sided p = 0.006')


def figure1(orf, out):
    fig = plt.figure(figsize=(14, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.6], wspace=0.45)
    counts(fig.add_subplot(gs[0]), orf, 'a  Glycogene up (ORF, U2OS)')
    scatter(fig.add_subplot(gs[1]), orf, 'ORF', 'b  Staining when each gene is up', lim=((-6, 8), (-6, 16)))
    lyso(fig.add_subplot(gs[2]), orf)
    save(fig, out, 'Figure1')


def figure2(cr, out):
    fig = plt.figure(figsize=(10, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.6)
    counts(fig.add_subplot(gs[0]), cr, 'a  Glycogene down (CRISPR, U2OS)')
    scatter(fig.add_subplot(gs[1]), cr, 'CRISPR', 'b  Staining when each gene is down')
    save(fig, out, 'Figure2')


def figure3(data, out):
    fig, ax = plt.subplots(figsize=(7, 3.6))
    groups, x = [], 0
    for t, tl in (('ORF', 'Up'), ('CRISPR', 'Down')):
        for lec in ('ConA', 'WGA'):
            g = data[t]
            z = g[f'{lec}_adj']
            ok = z.notna()
            ch = (z.abs() >= THR) & ok
            a, b = (ch & g.glycogene).sum(), (~ch & ok & g.glycogene).sum()
            c, d = (ch & ~g.glycogene).sum(), (~ch & ok & ~g.glycogene).sum()
            p = fisher_exact([[a, b], [c, d]])[1]
            r1, r2 = a / (a + b) * 100, c / (c + d) * 100
            ax.bar(x - 0.18, r1, width=0.34, color=UP, label='Glycogenes' if x == 0 else None)
            ax.bar(x + 0.18, r2, width=0.34, color=GRAY, label='Other genes' if x == 0 else None)
            ax.text(x, max(r1, r2) + 0.8, f'p = {p:.2g}', ha='center', fontsize=7.5, color=INK2)
            groups.append(f'{tl}\n{lec}')
            x += 1
    ax.set_xticks(range(len(groups)), groups)
    ax.set_ylabel('Genes whose staining changed (%)\n|adjusted z| ≥ 2')
    ax.set_ylim(0, 22)
    ax.legend(frameon=False, fontsize=7.5, loc='upper right')
    ax.grid(axis='y', color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    title(ax, 'Glycogenes change the staining about as often as other genes')
    save(fig, out, 'Figure3')


def figure4(r, out):
    f = pd.read_csv(r / '05_followup/A2_fisher_raw_vs_adjusted.csv')
    t = pd.read_csv(r / '11_toxicity/compounds_toxicity.csv')
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), gridspec_kw={'wspace': 0.45, 'width_ratios': [1.4, 0.6]})
    ax = axes[0]
    labels, x = [], 0
    for _, row in f.iterrows():
        xx = x + (0.6 if row.staining == 'adjusted' else 0)
        ax.bar(xx - 0.18, row.rate_predicted * 100, width=0.34, color=UP, label='Predicted compounds' if x == 0 else None)
        ax.bar(xx + 0.18, row.rate_other * 100, width=0.34, color=GRAY, label='Other compounds' if x == 0 else None)
        ax.text(xx, max(row.rate_predicted, row.rate_other) * 100 + 0.5, f'p = {row.fisher_p:.2g}', ha='center',
                fontsize=7, color=INK2)
        labels.append((xx, f"{row.lectin}\n{'stronger' if row.direction == 'up' else 'weaker'}"))
        x += 1
    ax.set_xticks([l[0] for l in labels], [l[1] for l in labels], fontsize=7.5)
    ax.set_ylim(0, 18)
    ax.axvline(3.8, color=INK2, lw=0.6, ls=':')
    ax.text(1.5, 16.2, 'Before adjustment', ha='center', fontsize=8, color=INK)
    ax.text(6.1, 16.2, 'After adjustment', ha='center', fontsize=8, color=INK)
    ax.set_ylabel('Compounds whose staining changed\nin the predicted direction (%)')
    ax.legend(frameon=False, fontsize=7.5, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    ax.grid(axis='y', color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    title(ax, 'a  Prediction from glycogene up/down (LINCS, A549)')
    ax = axes[1]
    rates = [t[t.conA_bright_all].toxic.mean(), t[~t.conA_bright_all].toxic.mean()]
    ns = [(int(t[t.conA_bright_all].toxic.sum()), int(t.conA_bright_all.sum())),
          (int(t[~t.conA_bright_all].toxic.sum()), int((~t.conA_bright_all).sum()))]
    ax.bar([0, 1], np.array(rates) * 100, width=0.55, color=[UP, GRAY])
    for xi, (v, (a, n)) in enumerate(zip(rates, ns)):
        ax.text(xi, v * 100 + 1.5, f'{v:.0%}\n({a}/{n})', ha='center', fontsize=7.5, color=INK2)
    ax.set_xticks([0, 1], ['ConA\nstronger', 'Other'])
    ax.set_ylim(0, 85)
    ax.set_ylabel('Toxic compounds (%)')
    ax.grid(axis='y', color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    title(ax, 'b  Cytotoxicity')
    save(fig, out, 'Figure4')


def main(r, out):
    r, out = Path(r), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    data = {t: pd.read_csv(r / f'13_jump/{t}_gene_scores.csv', index_col=0) for t in ('ORF', 'CRISPR')}
    figure1(data['ORF'], out)
    figure2(data['CRISPR'], out)
    figure3(data, out)
    figure4(r, out)
    print('saved Figure1-4 to', out)


if __name__ == '__main__':
    main(*sys.argv[1:3])
