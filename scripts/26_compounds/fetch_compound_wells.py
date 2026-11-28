"""Read the columns this study needs from the JUMP-CP compound well-level profiles (user, 2026-11-23: "摂動化合物 11 万
種類のやつも一通り回してサプリに突っ込んどいて").

Source: s3://cellpainting-gallery/cpg0016-jump/source_all/workspace/profiles_assembled/COMPOUND/v1.0a/profiles.parquet
(raw per-well CellProfiler aggregates of all compound plates; 804,844 wells). Parquet column projection over S3 reads
only these columns: Metadata_{Source, Plate, Well, JCP2022, Batch, PlateType, pert_type}, Cells_Intensity_MeanIntensity_
{ER, AGP, RNA, DNA, Mito}, Cells_AreaShape_Area, Nuclei_Number_Object_Number (the same columns as the ORF / CRISPR
analysis, JUMP_GENETIC_WELL_FEATURES).
Output: ~/research/glyco_cell_painting/data/jump_compound/compound_wells_raw.parquet
Usage: python fetch_compound_wells.py
"""
from pathlib import Path

import pyarrow.parquet as pq
from pyarrow import fs

SRC = 'cellpainting-gallery/cpg0016-jump/source_all/workspace/profiles_assembled/COMPOUND/v1.0a/profiles.parquet'
OUT = Path.home() / 'research/glyco_cell_painting/data/jump_compound/compound_wells_raw.parquet'
COLS = (['Metadata_Source', 'Metadata_Plate', 'Metadata_Well', 'Metadata_JCP2022', 'Metadata_Batch',
         'Metadata_PlateType', 'Metadata_pert_type']
        + [f'Cells_Intensity_MeanIntensity_{c}' for c in ('ER', 'AGP', 'RNA', 'DNA', 'Mito')]
        + ['Cells_AreaShape_Area', 'Nuclei_Number_Object_Number'])


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    s3 = fs.S3FileSystem(anonymous=True, region='us-east-1')
    t = pq.read_table(SRC, columns=COLS, filesystem=s3)
    pq.write_table(t, OUT)
    df = t.to_pandas()
    print(OUT, df.shape)
    print(df.Metadata_pert_type.value_counts(dropna=False).to_string())
    print('compounds (JCP2022):', df.Metadata_JCP2022.nunique(), 'plates:', df.Metadata_Plate.nunique(),
          'sources:', sorted(df.Metadata_Source.unique()))


if __name__ == '__main__':
    main()
