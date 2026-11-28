"""Compose a manuscript figure = existing panels a-c (unchanged) above the image panel d (make_panel.py, 600 dpi)
(user, 2026-11-23: "既存の Figure 3a–c は維持し、新しい Figure 3d として画像パネルを追加"; Figure 5d added the same way).
Figure 3: a-c = results/12_manuscript_figures/sections_v8/Figure1.png (300 dpi, section_figures_v5.py);
Figure 5: a-c = docs/figures/external_validation.png (600 dpi, scripts/24_periscope/external_figure.py).
Both inputs are resampled (Lanczos) to the width of a-c at 600 dpi (a-c dpi given as the 4th argument); no other
change. White gap of 0.05 inch between them.
Usage: python compose_figure.py <a-c png> <d png> <out png> <a-c dpi>
"""
import sys

from PIL import Image


def main(abc, d, out, dpi='300'):
    a = Image.open(abc).convert('RGB')
    b = Image.open(d).convert('RGB')
    w = round(a.width * 600 / float(dpi))
    a = a.resize((w, round(a.height * w / a.width)), Image.LANCZOS)
    b = b.resize((w, round(b.height * w / b.width)), Image.LANCZOS)
    gap = 30
    c = Image.new('RGB', (w, a.height + gap + b.height), 'white')
    c.paste(a, (0, 0))
    c.paste(b, (0, a.height + gap))
    c.save(out, dpi=(600, 600))
    print(out, c.size)


if __name__ == '__main__':
    main(*sys.argv[1:5])
