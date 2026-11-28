"""Image panel Figure 3c (overexpression and knockout): JUMP-CP images of the genes chosen by
select_genes.py (target vs same-plate negative control) and the image-level check (user, 2026-11-23).
Run after fetch_images.py.

Site selection (rule 3 of fetch_images.py, fixed before viewing): within each well, a site passes QC if its nucleus
count and DNA focus are >= 0.5 x the well median and the saturated fraction is < 0.1% in ER, AGP and DNA; the shown
site is the passing site whose cell-region mean in the channel of interest is closest to the median of the passing
sites. Only the plate chosen by rule 1 (median replicate) is shown; all plates are summarised in the CSV.
Display: the central 540 x 540 px of the 1080 x 1080 site (crop position fixed, not chosen); for each pair
(control, target) and channel, one linear display range = 0.5th to 99.7th percentile of the pooled pixels of the two
crops; no gamma, no other processing. Merge = channel of interest in green, DNA in blue, same ranges.
Scale bar from Metadata_Pixel_Size_Microns of jump-cellpainting/datasets metadata/microscope_config.csv (commit 016e865):
source_4 Opera Phenix 0.597976 um (equal to load_data Metadata_ImageResolutionX), source_13 CV7000 0.65 um (the TIFF
resolution tags are not usable: sensor pixel for source_4, a 72-dpi placeholder for source_13).
Outputs: docs/figures/Figure3c_panel.png (600 dpi) and .pdf (composed into manuscript/figure/Figure3.png by compose_figure.py),
<out>/image_selection.csv (chosen sites with plate, well, site, source URLs, display ranges),
<out>/image_level_summary.csv (per plate: target / control cell-region mean ratio over passing sites).
Usage: python make_panel.py <results/25_images dir> <png out> <pdf out> <panel letter>  (reads <dir>/oe and <dir>/ko)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tifffile  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from panels import panel_genes  # noqa: E402

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'pdf.fonttype': 42})
CH_LABEL = {'AGP': 'AGP (WGA + phalloidin)', 'ER': 'ER (ConA)', 'DNA': 'DNA (Hoechst)'}
CROP = 540
MC = pd.read_csv(Path.home() / 'research/glyco_cell_painting/data/jump/meta/microscope_config.csv')
PIX_UM = dict(zip(MC.Metadata_Source, MC.Metadata_Pixel_Size_Microns))


def choose_site(s, ch):
    s = s.copy()
    ok = ((s.n_nuclei >= 0.5 * s.n_nuclei.median()) & (s.focus >= 0.5 * s.focus.median())
          & (s[[f'{c}_saturated' for c in ('ER', 'AGP', 'DNA')]] < 0.001).all(axis=1))
    s['qc_pass'] = ok
    p = s[ok]
    target = p[f'{ch}_cell_mean'].median()
    pick = p.iloc[(p[f'{ch}_cell_mean'] - target).abs().argsort().iloc[0]]
    return pick, s


def crop(path):
    im = tifffile.imread(path).astype(float)
    h, w = im.shape
    y0, x0 = (h - CROP) // 2, (w - CROP) // 2
    return im[y0:y0 + CROP, x0:x0 + CROP]


def scale(im, lo, hi):
    return np.clip((im - lo) / (hi - lo), 0, 1)


def summarise(r, order):
    """Image-level summary over all plates of one panel -> <r>/image_level_summary.csv."""
    sm = pd.read_csv(r / 'site_measurements.csv')
    summ = []
    for (gene, plate), d in sm.groupby(['target', 'plate']):
        ch = dict((g, c) for g, _, c in order)[gene]
        res = {}
        for role in ('target', 'control'):
            _, q = choose_site(d[d.role == role], ch)
            res[role] = q[q.qc_pass]
        row = dict(target=gene, plate=plate, channel=ch, selected_plate=bool(d.selected_plate.iloc[0]))
        for c in ('ER', 'AGP', 'DNA'):
            row[f'{c}_ratio'] = res['target'][f'{c}_cell_mean'].median() / res['control'][f'{c}_cell_mean'].median()
        row['nuclei_target'] = res['target'].n_nuclei.median()
        row['nuclei_control'] = res['control'].n_nuclei.median()
        row['sites_pass_target'] = len(res['target'])
        row['sites_pass_control'] = len(res['control'])
        summ.append(row)
    summ = pd.DataFrame(summ)
    summ.to_csv(r / 'image_level_summary.csv', index=False)
    return sm, summ


def main(root, png, pdf, letter='c'):
    root = Path(root)
    fig, axs = plt.subplots(2, 6, figsize=(7.1, 3.15), gridspec_kw=dict(wspace=0.03, hspace=0.42, left=0.06,
                                                                           right=0.995, top=0.84, bottom=0.01))
    for k, panel in enumerate(('oe', 'ko')):
        r = root / panel
        order = [(g, 'OE' if m == 'ORF' else 'KO', c) for g, m, c, _ in panel_genes(panel)]
        sm, summ = summarise(r, order)
        sel_rows = []
        for i, (gene, kind, ch) in enumerate(order):
            d = sm[(sm.target == gene) & sm.selected_plate]
            picks = {role: choose_site(d[d.role == role], ch)[0] for role in ('control', 'target')}
            pix = PIX_UM['source_4' if kind == 'OE' else 'source_13']
            ims = {role: crop(picks[role][f'local_{ch}']) for role in picks}
            pooled = np.concatenate([ims['control'].ravel(), ims['target'].ravel()])
            rng = (np.percentile(pooled, 0.5), np.percentile(pooled, 99.7))
            lab = 'ORF' if kind == 'OE' else 'KO'
            for j, role in enumerate(('control', 'target')):
                ax = axs[k, 2 * i + j]
                ax.imshow(scale(ims[role], *rng), cmap='gray', vmin=0, vmax=1, interpolation='nearest')
                ax.set_xticks([]); ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                ax.set_title('Control' if role == 'control' else lab, fontsize=6, pad=2)
                if j == 1:
                    bar = 50 / pix
                    ax.plot([CROP - 30 - bar, CROP - 30], [CROP - 30, CROP - 30], color='white', lw=1.4)
                    ax.text(CROP - 30 - bar / 2, CROP - 45, '50 µm', color='white', ha='center', va='bottom',
                            fontsize=4.8)
            axs[k, 2 * i].annotate(f'$\\it{{{gene}}}$ {lab} · {CH_LABEL[ch]}', xy=(1.015, 1.13),
                                   xycoords='axes fraction', ha='center', va='bottom', fontsize=6.3, fontweight='bold')
            for role in ('control', 'target'):
                pk = picks[role]
                sel_rows.append(dict(target=gene, role=role, plate=pk.plate, well=pk.well, site=pk.site,
                                     channel_of_interest=ch, n_nuclei=pk.n_nuclei, focus=pk.focus,
                                     cell_mean_interest=pk[f'{ch}_cell_mean'], url_interest=pk[f'url_{ch}'],
                                     display_lo_interest=rng[0], display_hi_interest=rng[1], pixel_um=pix,
                                     crop=f'central {CROP}x{CROP} px'))
        axs[k, 0].set_ylabel('Overexpression' if panel == 'oe' else 'Knockout', fontsize=6.5, fontweight='bold',
                             labelpad=3)
        pd.DataFrame(sel_rows).to_csv(r / 'image_selection.csv', index=False)
        pd.set_option('display.width', 220)
        print(summ.round(3).to_string(index=False))
        print(pd.DataFrame(sel_rows)[['target', 'role', 'plate', 'well', 'site', 'n_nuclei',
                                      'cell_mean_interest']].to_string(index=False))
    fig.text(0.005, 0.99, letter, fontsize=9, fontweight='bold', va='top')
    fig.savefig(png, dpi=600, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')


if __name__ == '__main__':
    main(*sys.argv[1:5])
