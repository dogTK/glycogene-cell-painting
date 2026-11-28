"""Fetch the columns this study needs from the PERISCOPE gene-level profiles (Cell Painting Gallery cpg0021-periscope,
public) (user, 2026-11-23: external validation, "3 かな", "tsui もやっといて").

Source: s3://cellpainting-gallery/cpg0021-periscope/broad/workspace/profiles/<A549|HeLa>/
        *_gene_normalized_ALLBATCHES___<plate>___ALLWELLS.csv.gz (one row per gene per plate; features normalized by
        the PERISCOPE pipeline). Genome-wide CRISPR knockout, optical pooled screening; stains ConA, WGA, Phalloidin
        (separate channels), DAPI, anti-TOMM20 antibody (mitochondria; column "Mito").
Each file is streamed (aws s3 cp - | gunzip) and only these columns are kept: Metadata_Foci_Barcode_MatchedTo_GeneCode,
Cells_Intensity_MeanIntensity_{ConA, WGA, Phalloidin, DAPI_Painting, Mito}, Cells_AreaShape_Area, Nuclei_AreaShape_Area.
Output (untrusted downloaded data kept apart from scripts): ~/research/glyco_cell_painting/data/periscope/
<cell>_<plate>.parquet and manifest.csv (source key, rows, sha256 of the downloaded bytes).
Usage: python fetch_periscope.py
"""
import hashlib
import io
import subprocess
from pathlib import Path

import pandas as pd

OUT = Path.home() / 'research/glyco_cell_painting/data/periscope'
B = 's3://cellpainting-gallery/cpg0021-periscope/broad/workspace/profiles'
COLS = ['Metadata_Foci_Barcode_MatchedTo_GeneCode', 'Cells_Intensity_MeanIntensity_ConA',
        'Cells_Intensity_MeanIntensity_WGA', 'Cells_Intensity_MeanIntensity_Phalloidin',
        'Cells_Intensity_MeanIntensity_DAPI_Painting', 'Cells_Intensity_MeanIntensity_Mito', 'Cells_AreaShape_Area',
        'Nuclei_AreaShape_Area']
FILES = ([('A549', p, f'A549/20200805_A549_WG_Screen_gene_normalized_ALLBATCHES___{p}___ALLWELLS.csv.gz')
          for p in ['CP186A', 'CP186B', 'CP186C', 'CP186D', 'CP186E', 'CP186F', 'CP186G', 'CP186H', 'CP186N']]
         + [('HeLa', p, f'HeLa/20210422_6W_CP257_gene_normalized_ALLBATCHES___{p}___ALLWELLS.csv.gz')
            for p in ['CP257A', 'CP257B', 'CP257D', 'CP257F', 'CP257H', 'CP257J', 'CP257K', 'CP257L', 'CP257N']])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    man = []
    for cell, plate, key in FILES:
        out = OUT / f'{cell}_{plate}.parquet'
        if out.exists():
            print('skip', out.name, flush=True)
            continue
        raw = subprocess.run(['aws', 's3', 'cp', '--no-sign-request', '--region', 'us-east-1', f'{B}/{key}', '-'],
                             check=True, capture_output=True).stdout
        sha = hashlib.sha256(raw).hexdigest()
        df = pd.read_csv(io.BytesIO(raw), compression='gzip', usecols=COLS)
        df.insert(0, 'Metadata_Cell', cell)
        df.insert(1, 'Metadata_Plate', plate)
        df.to_parquet(out, index=False)
        man.append(dict(cell=cell, plate=plate, source=f'{B}/{key}', rows=len(df), bytes=len(raw), sha256=sha))
        print(cell, plate, len(df), 'genes', flush=True)
        del raw
    if man:
        m = pd.DataFrame(man)
        f = OUT / 'manifest.csv'
        m = pd.concat([pd.read_csv(f), m]) if f.exists() else m
        m.to_csv(f, index=False)


if __name__ == '__main__':
    main()
