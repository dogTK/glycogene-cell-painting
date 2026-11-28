"""Supplementary figure: compound perturbations compared with ORF and CRISPR (compounds_analysis.py output).
a share of perturbations with |z| >= 2 (robust z against same-plate negative controls, median over >= 3 wells) per
  channel, cell area and cell count; b share with cell-count z <= -2 among all perturbations and among those with
  |ER| >= 2 or |AGP| >= 2; c z values of the glycosylation inhibitors found in JUMP-CP (one row per compound, JUMP ID).
No title in the image.
Usage: python compounds_figure.py <results/26_compounds dir> <out stem>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'pdf.fonttype': 42, 'axes.linewidth': 0.6})
MOD = [('ORF', 'OE (ORF)', '#1f66c2'), ('CRISPR', 'KO (CRISPR)', '#d9541e'), ('compound', 'Compound', '#3a3a38')]
FEAT = [('ER', 'ER\n(ConA)'), ('AGP', 'AGP\n(WGA)'), ('RNA', 'RNA'), ('DNA', 'DNA'), ('Mito', 'Mito'),
        ('area', 'Cell\narea'), ('count', 'Cell\ncount')]


def main(r, stem):
    r = Path(r)
    s = pd.read_csv(r / 'modality_summary.csv').set_index('modality')
    nm = pd.read_csv(r / 'named_compounds.csv')
    fig = plt.figure(figsize=(7.1, 2.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[2.3, 1.0, 1.5], wspace=0.95, left=0.06, right=0.98, bottom=0.2, top=0.88)
    ax = fig.add_subplot(gs[0])
    x = np.arange(len(FEAT)); wdt = 0.26
    for i, (m, lab, col) in enumerate(MOD):
        ax.bar(x + (i - 1) * wdt, [s.loc[m, f'share_abs_z_ge2_{f}'] * 100 for f, _ in FEAT], wdt, color=col, label=
               f'{lab} (n = {int(s.loc[m, "n_perturbations"]):,})')
    ax.set_xticks(x, [l for _, l in FEAT], fontsize=6)
    ax.set_ylabel('Perturbations with |z| ≥ 2 (%)')
    ax.legend(frameon=False, fontsize=5.8, loc='upper left')
    ax.text(-0.12, 1.22, 'A', transform=ax.transAxes, fontsize=9, fontweight='bold')
    ax = fig.add_subplot(gs[1])
    x = np.arange(len(MOD))
    ax.bar(x - 0.2, [s.loc[m, 'share_toxic_all'] * 100 for m, _, _ in MOD], 0.38, color='#cfcdc6', label='All')
    ax.bar(x + 0.2, [s.loc[m, 'share_toxic_among_lectin_ge2'] * 100 for m, _, _ in MOD], 0.38,
           color=[c for _, _, c in MOD], label='|ER| or |AGP| ≥ 2')
    ax.set_xticks(x, ['OE', 'KO', 'Cpd'], fontsize=6)
    ax.set_ylabel('Cell count z ≤ −2 (%)')
    ax.legend(frameon=False, fontsize=5.6, loc='lower left', bbox_to_anchor=(0, 1.0), ncol=1)
    ax.text(-0.45, 1.22, 'B', transform=ax.transAxes, fontsize=9, fontweight='bold')
    ax = fig.add_subplot(gs[2])
    inh = nm[nm.group == 'glycosylation inhibitor'].drop_duplicates('Metadata_JCP2022')
    inh = inh.assign(label=inh.name.replace({'1-deoxynojirimycin': 'DNJ / DMJ'}) + ' (' + inh.n_wells.astype(int).astype(str) + ')')
    cols = ['ER', 'AGP', 'RNA', 'DNA', 'Mito', 'count']
    m = inh[cols].to_numpy()
    im = ax.imshow(m, cmap='RdBu_r', vmin=-3, vmax=3, aspect='auto')
    for (i, j), v in np.ndenumerate(m):
        ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=5.3, color='black')
    ax.set_xticks(range(len(cols)), ['ER', 'AGP', 'RNA', 'DNA', 'Mito', 'Count'], fontsize=6, rotation=45)
    ax.set_yticks(range(len(inh)), inh.label, fontsize=6)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label('z', fontsize=6); cb.ax.tick_params(labelsize=5.5)
    ax.text(-0.95, 1.22, 'C', transform=ax.transAxes, fontsize=9, fontweight='bold')
    for a in fig.axes[:2]:
        for sp in ('top', 'right'):
            a.spines[sp].set_visible(False)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{stem}.{ext}', dpi=600 if ext == 'png' else None, bbox_inches='tight')


if __name__ == '__main__':
    main(*sys.argv[1:3])
