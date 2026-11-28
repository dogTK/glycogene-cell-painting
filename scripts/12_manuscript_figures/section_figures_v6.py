"""Figure 3 made shorter (user, 2026-11-23: "Fig3そもそもデカすぎる" -> option agreed: the ECDF panel goes to the
Supplementary, the image panel loses its DNA columns).
  Figure 3a  share of glycogenes above the other genes' 95th percentile per channel (former b)
  Figure 3b  gene-level strip plot of the glycogenes' ConA / WGA effects (former c, slightly shorter)
  (Figure 3c = image panel, scripts/25_images/make_panel.py; composed by scripts/25_images/compose_figure.py)
  Supplementary figure: ECDF of |Glass's delta|, glycogenes vs other genes, 5 channels x ORF / CRISPR (former a).
Same inputs and drawing functions as section_figures_v5.py; no number changes. Common style (scripts/figure_style.py,
user 2026-11-24): error bars and marker edges 0.7.
Usage: python section_figures_v6.py <results dir> <out dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from section_figures_v4 import CONDS, NAMES, strip_panel  # noqa: E402
from section_figures_v5 import GLY, OTH, dot_panel, ecdf_grid  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from figure_style import apply  # noqa: E402


def main(r, out):
    r, out = Path(r), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    c = r / '17_channel_specificity'
    eff = {t: pd.read_csv(c / f'{t}_main_effects.csv', index_col=0) for t in ('ORF', 'CRISPR')}
    summ = pd.read_csv(c / 'channel_summary_main.csv')
    d = r / '13_jump/vs_other_genes_infection_filtered'
    calls = {t: pd.read_csv(d / f'{t}_glycogene_calls.csv', index_col=0) for t in ('ORF', 'CRISPR')}

    apply()
    fig = plt.figure(figsize=(7.1, 5.3))
    gb = fig.add_gridspec(1, 1, left=0.1, right=0.98, top=0.95, bottom=0.70)
    dot_panel(fig.add_subplot(gb[0, 0]), summ, 'share_above', 'binom_p',
              'Glycogenes above the 95th\npercentile of other genes', 'share_lo', 'share_hi', ref=0.05)
    gc = fig.add_gridspec(1, 1, left=0.17, right=0.98, top=0.60, bottom=0.11)
    ax = fig.add_subplot(gc[0, 0])
    strip_panel(ax, [(f'{NAMES[t]} · {lec}', calls[t], lec) for t, lec in CONDS])
    for t in ax.texts:
        if t.get_text() in ('b', 'B'):
            t.set_visible(False)
    for letter, (x, y) in {'A': (0.01, 0.99), 'B': (0.01, 0.635)}.items():
        fig.text(x, y, letter, fontsize=9, fontweight='bold', va='top')
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'Figure3ab.{ext}', dpi=600 if ext == 'png' else None)
    plt.close(fig)

    fig = plt.figure(figsize=(7.1, 3.6))
    ga = fig.add_gridspec(2, 5, left=0.1, right=0.98, top=0.88, bottom=0.12, wspace=0.12, hspace=0.18)
    ecdf_grid(fig, ga, eff, summ)
    fig.legend([plt.Line2D([], [], color=GLY, lw=1.3), plt.Line2D([], [], color=OTH, lw=1.1)],
               ['Glycogenes', 'Other genes'], loc='upper right', bbox_to_anchor=(0.98, 0.995), ncol=2, frameon=False)
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'FigureS_ecdf.{ext}', dpi=600 if ext == 'png' else None)
    plt.close(fig)


if __name__ == '__main__':
    main(*sys.argv[1:3])
