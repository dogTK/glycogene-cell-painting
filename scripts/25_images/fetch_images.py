"""Download the JUMP-CP raw images for an image panel (oe, ko; genes from select_genes.py) and measure every site (user, 2026-11-23).

Selection rules (fixed before any image was viewed):
  1. Plate: of the 5 replicate plates, the plate whose target-well adjusted value (the tested quantity) in the channel
     of interest is the median of the five.
  2. Control well, on every replicate plate: the negative-control well (ORF: BFP, HcRed, LUCIFERASE or LacZ; CRISPR:
     non-targeting only) whose raw cell mean intensity in the channel of interest is closest to the median of the
     plate's eligible negative-control wells.
  3. Sites (applied later, in make_panel.py): all 9 sites are measured; a site passes QC if its nucleus count and its
     DNA focus score (variance of the Laplacian) are each >= 0.5 x the median of the well's sites and the saturated
     pixel fraction is < 0.1% in every channel. The displayed site is the passing site whose cell-region mean in the
     channel of interest is closest to the median over the well's passing sites.
Image paths come only from the official load_data.csv of each plate (Cell Painting Gallery,
cpg0016-jump/<source>/workspace/load_data_csv/<batch>/<plate>/load_data.csv); channels by the URL_Orig<Channel> columns.
Images are stored outside the repository (~/research/glyco_cell_painting/data/jump_images/raw/), never committed.
Per-site measurements (raw pixels, no illumination correction; the same optics apply to both wells of a plate):
  nuclei: Otsu threshold of the Gaussian-smoothed (sigma 2) DNA image, connected components >= 40 px
  cell region: Otsu threshold of the smoothed ER + AGP sum
  focus: variance of the Laplacian of the DNA image
  mean intensity of ER, AGP, DNA over the cell region and over the whole image; saturated fraction (= 65535)
Output: <out>/site_measurements.csv, <out>/selection_wells.csv
Usage: python fetch_images.py <well_values.csv> <out_dir> <panel> [selected|all]  (selected: only the plate of rule 1)
"""
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

G = 's3://cellpainting-gallery/cpg0016-jump'
RAW = Path.home() / 'research/glyco_cell_painting/data/jump_images/raw'
LD = Path.home() / 'research/glyco_cell_painting/data/jump_images/load_data'
sys.path.insert(0, str(Path(__file__).resolve().parent))
from panels import SOURCE, panel_genes  # noqa: E402
NEG_ORF = {'JCP2022_915128', 'JCP2022_915129', 'JCP2022_915130', 'JCP2022_915131'}
NEG_CRISPR = {'JCP2022_800002'}
CHS = ['ER', 'AGP', 'DNA']
PLATE_BATCH = {}


def threshold_otsu(x, nbins=256):
    """Otsu's threshold (same definition as skimage.filters.threshold_otsu, histogram of 256 bins)."""
    hist, edges = np.histogram(x.ravel(), bins=nbins)
    centers = (edges[:-1] + edges[1:]) / 2
    w1 = np.cumsum(hist); w2 = np.cumsum(hist[::-1])[::-1]
    m1 = np.cumsum(hist * centers) / np.maximum(w1, 1)
    m2 = (np.cumsum((hist * centers)[::-1]) / np.maximum(w2[::-1], 1))[::-1]
    var = w1[:-1] * w2[1:] * (m1[:-1] - m2[1:]) ** 2
    return centers[:-1][np.argmax(var)]


def s3cp(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.exists():
        subprocess.run(['aws', 's3', 'cp', '--no-sign-request', '--region', 'us-east-1', '--only-show-errors',
                        src, str(dst)], check=True)
    return dst


def load_data(source, batch, plate):
    return pd.read_csv(s3cp(f'{G}/{source}/workspace/load_data_csv/{batch}/{plate}/load_data.csv',
                            LD / f'{plate}_load_data.csv'))


def measure(paths):
    im = {c: tifffile.imread(paths[c]).astype(float) for c in CHS}
    dna = ndi.gaussian_filter(im['DNA'], 2)
    nuc = dna > threshold_otsu(dna)
    lab, n = ndi.label(nuc)
    sizes = ndi.sum(nuc, lab, range(1, n + 1))
    n_nuc = int((sizes >= 40).sum())
    cs = ndi.gaussian_filter(im['ER'] + im['AGP'], 2)
    cell = cs > threshold_otsu(cs)
    focus = float(ndi.laplace(im['DNA']).var())
    r = dict(n_nuclei=n_nuc, focus=focus, cell_fraction=float(cell.mean()))
    for c in CHS:
        r[f'{c}_cell_mean'] = float(im[c][cell].mean())
        r[f'{c}_image_mean'] = float(im[c].mean())
        r[f'{c}_saturated'] = float((im[c] >= 65535).mean())
    return r


def main(wv_csv, out, panel, mode='selected'):
    TARGETS = {g: (m, SOURCE[m], c) for g, m, c, _ in panel_genes(panel)}
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    wv = pd.read_csv(wv_csv)
    plates = pd.read_csv(Path.home() / 'research/glyco_cell_painting/data/jump/meta/plate.csv.gz')
    batch = dict(zip(plates.Metadata_Plate, plates.Metadata_Batch))
    sel, jobs = [], []
    for gene, (mod, source, ch) in TARGETS.items():
        d = wv[wv.target == gene]
        tgt = d[d.role == 'target']
        med_plate = tgt.sort_values(f'adj_{ch}').iloc[len(tgt) // 2].plate
        neg_ids = NEG_ORF if mod == 'ORF' else NEG_CRISPR
        for p in sorted(tgt.plate):
            t = tgt[tgt.plate == p].iloc[0]
            n = d[(d.plate == p) & (d.role == 'negcon') & d.jcp.isin(neg_ids) & ~d.low_infection]
            c = n.iloc[(n[f'raw_{ch}'] - n[f'raw_{ch}'].median()).abs().argsort().iloc[0]]
            for role, w in (('target', t), ('control', c)):
                sel.append(dict(target=gene, modality=mod, plate=p, well=w.well, jcp=w.jcp, role=role,
                                channel_of_interest=ch, selected_plate=p == med_plate,
                                raw_value=w[f'raw_{ch}'], adj_value=w.get(f'adj_{ch}', np.nan)))
            if mode == 'selected' and p != med_plate:
                continue
            ld = load_data(source, batch[p], p)
            for role, w in (('target', t.well), ('control', c.well)):
                for _, s in ld[ld.Metadata_Well == w].iterrows():
                    urls = {cc: s[f'URL_Orig{cc}'] for cc in CHS}
                    jobs.append(dict(target=gene, plate=p, well=w, site=int(s.Metadata_Site), role=role,
                                     selected_plate=p == med_plate, source=source, urls=urls))
    todo = [(u, RAW / j['source'] / j['plate'] / Path(u).name) for j in jobs for u in j['urls'].values()]
    print('images to fetch:', len(todo), flush=True)
    with ThreadPoolExecutor(max_workers=16) as ex:
        list(ex.map(lambda a: s3cp(*a), todo))
    rows = []
    for j in jobs:
        paths = {cc: RAW / j['source'] / j['plate'] / Path(u).name for cc, u in j['urls'].items()}
        m = measure(paths)
        rows.append(dict(target=j['target'], plate=j['plate'], well=j['well'], site=j['site'], role=j['role'],
                         selected_plate=j['selected_plate'], **{f'url_{cc}': u for cc, u in j['urls'].items()},
                         **{f'local_{cc}': str(paths[cc]) for cc in CHS}, **m))
    pd.DataFrame(sel).to_csv(out / 'selection_wells.csv', index=False)
    pd.DataFrame(rows).to_csv(out / 'site_measurements.csv', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:5])
