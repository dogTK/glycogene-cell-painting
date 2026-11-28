"""Common style for the manuscript figures (user, 2026-11-24: "全体的にグラフレイアウト全てのフィグで揃ってますか？"
-> one specification for every figure; also in the parent skill manuscript-figure-style).

  size   : width 7.1 in (full page width), height <= 9 in, PNG 600 dpi (PDF vector)
  text   : Arial; 7 pt body and axis labels, 6.5 pt ticks, 6 pt small labels / legends, 9 pt bold panel letters
  lines  : axes 0.6, reference lines 0.6 (dotted), error bars and marker edges 0.7, heatmap separators 0.6, grid 0.4
  axes   : top and right spines off; light horizontal grid only where values are read off the y axis
  no title inside the image (the legend goes in the manuscript)
Usage: from figure_style import apply, despine, letter, WIDTH, ...; call apply() before making the figure.
"""
import matplotlib.pyplot as plt

WIDTH, MAX_HEIGHT, DPI = 7.1, 9.0, 600
FS, FS_TICK, FS_SMALL, FS_LETTER = 7, 6.5, 6, 9
LW_AXIS, LW_REF, LW_ERR, LW_EDGE, LW_SEP, LW_GRID = 0.6, 0.6, 0.7, 0.7, 0.6, 0.4
INK, INK2, GRID = '#111111', '#55534f', '#e6e5e1'


def apply():
    plt.rcParams.update({'font.family': 'Arial', 'font.size': FS, 'axes.labelsize': FS, 'axes.titlesize': FS,
                         'xtick.labelsize': FS_TICK, 'ytick.labelsize': FS_TICK, 'legend.fontsize': FS_SMALL,
                         'axes.linewidth': LW_AXIS, 'xtick.major.width': LW_AXIS, 'ytick.major.width': LW_AXIS,
                         'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'lines.linewidth': 1.0,
                         'patch.linewidth': LW_EDGE, 'pdf.fonttype': 42, 'savefig.dpi': DPI})


def despine(ax):
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def letter(fig, x, y, s):
    fig.text(x, y, s, fontsize=FS_LETTER, fontweight='bold', va='top', ha='left')


def save(fig, stem):
    w, h = fig.get_size_inches()
    assert w <= WIDTH + 0.01 and h <= MAX_HEIGHT + 0.01, f'figure {w:.2f} x {h:.2f} in exceeds {WIDTH} x {MAX_HEIGHT}'
    for ext in ('png', 'pdf'):
        fig.savefig(f'{stem}.{ext}', dpi=DPI if ext == 'png' else None)
