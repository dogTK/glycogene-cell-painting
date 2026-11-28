"""Figure 4: overexpression vs knockout effects of the glycogenes tested by both (redrawn in the common style,
user 2026-11-24; numbers from all_genes_bidirectional.py: <dir>/bidirectional_<ConA|WGA>.csv,
bidirectional_correlations.csv). Same encoding and labels as all_genes_bidirectional.py: grey = p >= 0.05 in both,
light blue = p < 0.05 in one, blue = p < 0.05 in both with opposite sign, black ring = FDR < 0.05 (205 glycogenes)
in either; axes clipped at +-8 SD (value in brackets).
Usage: python orf_vs_ko_figure.py <results/16_all_genes dir> <out stem>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from figure_style import (FS_SMALL, GRID, INK, INK2, LW_EDGE, LW_GRID, LW_REF, WIDTH, apply, despine,  # noqa: E402
                          letter, save)

GLY, OTH = '#2a78d6', '#b4b2ab'
CLIP = 8
# label offsets (points) only to keep overlapping names apart; default (3, 2), right edge drawn to the left
NUDGE = {'B3GNT6': (3, 5), 'GAL3ST3': (-4, 3), 'GCK': (3, 5), 'HS3ST5': (3, -6), 'GBA': (-22, -7), 'GLB1': (3, 2),
         'UGCG': (3, -7), 'POMK': (3, 4)}


def main(r, stem):
    r = Path(r)
    apply()
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 3.6), gridspec_kw=dict(wspace=0.32, left=0.08, right=0.98,
                                                                           top=0.93, bottom=0.25))
    for ax, lec in zip(axes, ('ConA', 'WGA')):
        d = pd.read_csv(r / f'bidirectional_{lec}.csv', index_col=0)
        xc, yc = d.orf.clip(-CLIP, CLIP), d.ko.clip(-CLIP, CLIP)
        quad = np.where((d.orf_p < 0.05) & (d.ko_p < 0.05) & (np.sign(d.orf) != np.sign(d.ko)), 'both',
                        np.where((d.orf_p < 0.05) | (d.ko_p < 0.05), 'one', 'none'))
        for q, col, sz, lab in (('none', OTH, 6, 'p ≥ 0.05 in both'), ('one', '#9cc0ea', 10, 'p < 0.05 in one'),
                                ('both', GLY, 20, 'p < 0.05 in both, opposite sign')):
            m = quad == q
            ax.scatter(xc[m], yc[m], s=sz, color=col, lw=0, alpha=0.9, label=f'{lab} ({m.sum()})')
        sig = (d.orf_fdr < 0.05) | (d.ko_fdr < 0.05)
        ax.scatter(xc[sig], yc[sig], s=26, facecolor='none', edgecolor=INK, lw=LW_EDGE,
                   label='FDR < 0.05 (glycogenes) in either')
        ax.axhline(0, color=INK2, lw=LW_REF)
        ax.axvline(0, color=INK2, lw=LW_REF)
        top_opp = d.sort_values('opposite_product', ascending=False).head(6)
        lab = set(d.index[sig]) | set(top_opp.index) | set(d.index[quad == 'both'])
        lab |= set(d.index[(d.orf.abs() >= 2.5) | (d.ko.abs() >= 2)])
        for gene in lab:
            g = d.loc[gene]
            right = xc[gene] > 5
            dx, dy = NUDGE.get(gene, (-4, 3) if right else (3, 2))
            ax.annotate(gene + ('' if abs(g.orf) <= CLIP and abs(g.ko) <= CLIP else f' ({g.orf:+.0f}/{g.ko:+.0f})'),
                        (xc[gene], yc[gene]), xytext=(dx, dy), textcoords='offset points', fontsize=FS_SMALL,
                        ha='right' if dx < 0 else 'left',
                        fontstyle='italic', color=INK,
                        fontweight='bold' if quad[d.index.get_loc(gene)] == 'both' else 'normal')
        ax.set_xlabel(f'{lec}: overexpression effect (SD; clipped at ±{CLIP})')
        ax.set_ylabel(f'{lec}: knockout effect (SD)')
        ax.grid(color=GRID, lw=LW_GRID)
        ax.set_axisbelow(True)
        despine(ax)
        ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2, columnspacing=0.8, handletextpad=0.3)
    letter(fig, 0.005, 0.99, 'A')
    letter(fig, 0.5, 0.99, 'B')
    save(fig, stem)


if __name__ == '__main__':
    main(*sys.argv[1:3])
