"""Fetch the JUMP-CP ORF and CRISPR well-level features this study needs and load them into Snowflake.

Sources (public, anonymous):
- metadata: jump-cellpainting/datasets commit 016e865 (plate.csv.gz, well.csv.gz, orf.csv.gz, crispr.csv.gz,
  perturbation_control.csv)
- well-level CellProfiler profiles: s3://cellpainting-gallery/cpg0016-jump/<source>/workspace/profiles/<batch>/
  <plate>/<plate>.parquet, for every plate of type ORF or CRISPR; only the columns in COLUMNS are read (Parquet
  column selection over S3). These are raw per-well aggregates (not normalised).
Writes Parquet to ~/research/glyco_cell_painting/data/jump/ and loads RAW.GLYCO_CELL_PAINTING.JUMP_GENETIC_WELL_FEATURES,
JUMP_WELL, JUMP_PLATE, JUMP_ORF, JUMP_CRISPR, JUMP_PERTURBATION_CONTROL with comments, _LOAD_ID, COLUMN_NAMES and
_LOAD_MANIFEST (same route as 01_load_to_snowflake: Parquet -> S3 -> external stage -> COPY INTO). RAW is append-only:
stops if a table exists.
Usage: python fetch_jump_genetic.py
"""
import os
import hashlib
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from pyarrow import fs

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cursor  # noqa: E402

COMMIT = '016e865fa0691244e0860943e41c7d6a88ed2580'
META_URL = f'https://raw.githubusercontent.com/jump-cellpainting/datasets/{COMMIT}/metadata'
DATA = Path.home() / 'research/glyco_cell_painting/data/jump'
SCHEMA = 'RAW.GLYCO_CELL_PAINTING'
S3 = f"s3://{os.environ.get('GLYCO_S3_BUCKET', '<your-bucket>')}/glyco_cell_painting/raw"
STAGE = os.environ.get('GLYCO_SNOWFLAKE_STAGE', '@<your_stage>') + '/glyco_cell_painting/raw'
CH = ['ER', 'AGP', 'RNA', 'DNA', 'Mito']
COLUMNS = (['Metadata_Source', 'Metadata_Plate', 'Metadata_Well', 'Cells_AreaShape_Area', 'Nuclei_AreaShape_Area',
            'Nuclei_Number_Object_Number']
           + [f'Cells_Intensity_MeanIntensity_{c}' for c in CH]
           + [f'Cytoplasm_Intensity_MeanIntensity_{c}' for c in ('ER', 'AGP')]
           + [f'Cells_Intensity_IntegratedIntensity_{c}' for c in ('ER', 'AGP')])
META_FILES = {'JUMP_PLATE': 'plate.csv.gz', 'JUMP_WELL': 'well.csv.gz', 'JUMP_ORF': 'orf.csv.gz',
              'JUMP_CRISPR': 'crispr.csv.gz', 'JUMP_PERTURBATION_CONTROL': 'perturbation_control.csv'}
COMMENTS = {
    'JUMP_GENETIC_WELL_FEATURES': 'JUMP-CP (cpg0016) ORF and CRISPR plates, U2OS: raw well-level CellProfiler '
        'aggregates for the columns this study uses (5-channel cell mean intensities, cytoplasm ER/AGP, integrated ER/AGP, '
        'cell and nucleus area, Nuclei_Number_Object_Number as a cell-count proxy: mean object number ~ (N + 1) / 2). '
        'Not normalised. One row per well',
    'JUMP_PLATE': 'JUMP-CP plate metadata (source, batch, plate, plate type), jump-cellpainting/datasets commit 016e865',
    'JUMP_WELL': 'JUMP-CP well metadata (source, plate, well -> JCP2022 perturbation id), jump-cellpainting/datasets commit 016e865',
    'JUMP_ORF': 'JUMP-CP ORF perturbations (JCP2022 id, gene symbol, vector, transcript), jump-cellpainting/datasets commit 016e865',
    'JUMP_CRISPR': 'JUMP-CP CRISPR perturbations (JCP2022 id, gene symbol), jump-cellpainting/datasets commit 016e865',
    'JUMP_PERTURBATION_CONTROL': 'JUMP-CP control perturbations (JCP2022 id, positive/negative control type), '
        'jump-cellpainting/datasets commit 016e865',
}


def meta(name):
    out = DATA / 'meta' / name
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(['curl', '-sSLf', '-o', str(out), f'{META_URL}/{name}'], check=True)
    return out


def read_plate(row, s3):
    path = (f'cellpainting-gallery/cpg0016-jump/{row.Metadata_Source}/workspace/profiles/{row.Metadata_Batch}/'
            f'{row.Metadata_Plate}/{row.Metadata_Plate}.parquet')
    try:
        return pq.read_table(path, columns=COLUMNS, filesystem=s3).to_pandas(), None
    except Exception as e:  # noqa: BLE001  (record and continue; reported at the end)
        return None, f'{path}: {e}'


def load(cur, name, df, load_id, source_url, source_file, sha, existing):
    out = DATA / 'export' / load_id / name
    out.mkdir(parents=True, exist_ok=True)
    orig = [str(c) for c in df.columns]
    df = df.copy()
    df.columns = [c.upper() for c in orig]
    for c in df.select_dtypes('object').columns:
        df[c] = df[c].astype('string')
    df.insert(0, '_LOAD_ID', load_id)
    df.to_parquet(out / f'{name}.parquet', index=False)
    subprocess.run(['aws', 's3', 'sync', str(out), f'{S3}/{load_id}/{name}', '--only-show-errors'], check=True)
    loc = f'{STAGE}/{load_id}/{name}/'
    cur.execute(f"CREATE TABLE {SCHEMA}.{name} USING TEMPLATE (SELECT ARRAY_AGG(OBJECT_CONSTRUCT(*)) "
                f"FROM TABLE(INFER_SCHEMA(LOCATION => '{loc}', FILE_FORMAT => '{SCHEMA}.PARQUET_FMT')))")
    cur.execute(f'COMMENT ON TABLE {SCHEMA}.{name} IS $${COMMENTS[name]}$$')
    cur.execute(f"COPY INTO {SCHEMA}.{name} FROM '{loc}' FILE_FORMAT = (FORMAT_NAME = '{SCHEMA}.PARQUET_FMT') "
                'MATCH_BY_COLUMN_NAME = CASE_SENSITIVE')
    cn = pd.DataFrame({'_LOAD_ID': load_id, 'TABLE_NAME': name, 'POSITION': range(len(orig)), 'ORIGINAL_NAME': orig})
    cno = DATA / 'export' / load_id / f'COLUMN_NAMES_{name}'
    cno.mkdir(parents=True, exist_ok=True)
    cn.to_parquet(cno / 'COLUMN_NAMES.parquet', index=False)
    subprocess.run(['aws', 's3', 'sync', str(cno), f'{S3}/{load_id}/COLUMN_NAMES_{name}', '--only-show-errors'],
                   check=True)
    cur.execute(f"COPY INTO {SCHEMA}.COLUMN_NAMES FROM '{STAGE}/{load_id}/COLUMN_NAMES_{name}/' "
                f"FILE_FORMAT = (FORMAT_NAME = '{SCHEMA}.PARQUET_FMT') MATCH_BY_COLUMN_NAME = CASE_SENSITIVE")
    now = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
    cur.execute(f'INSERT INTO {SCHEMA}._LOAD_MANIFEST VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)',
                (load_id, now, name, len(df), source_url, COMMIT, source_file, sha,
                 'loaded by glyco_cell_painting scripts/13_jump/fetch_jump_genetic.py'))
    cur.execute(f'SELECT COUNT(*) FROM {SCHEMA}.{name} WHERE _LOAD_ID = %s', (load_id,))
    print(f'{name}: {len(df):,} rows; in Snowflake {cur.fetchone()[0]:,}')


def main():
    cur = cursor()
    cur.execute(f'SHOW TABLES IN SCHEMA {SCHEMA}')
    existing = {r[1] for r in cur.fetchall()}
    clash = existing & (set(META_FILES) | {'JUMP_GENETIC_WELL_FEATURES'})
    if clash:
        sys.exit(f'RAW is append-only; tables already exist: {sorted(clash)}')
    load_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')

    plates = pd.read_csv(meta('plate.csv.gz'))
    gen = plates[plates.Metadata_PlateType.isin(['ORF', 'CRISPR'])]
    print('plates to read:', gen.Metadata_PlateType.value_counts().to_dict())
    s3 = fs.S3FileSystem(anonymous=True, region='us-east-1')
    with ThreadPoolExecutor(max_workers=16) as ex:
        results = list(ex.map(lambda r: read_plate(r, s3), gen.itertuples()))
    errors = [e for _, e in results if e]
    feats = pd.concat([d for d, _ in results if d is not None], ignore_index=True)
    print(f'wells {len(feats):,} from {len(results) - len(errors)} plates; failed plates {len(errors)}')
    for e in errors[:10]:
        print('  ', e[:200])
    DATA.mkdir(parents=True, exist_ok=True)
    feats.to_parquet(DATA / 'jump_genetic_well_features.parquet', index=False)

    note = ('s3://cellpainting-gallery/cpg0016-jump/<source>/workspace/profiles/<batch>/<plate>/<plate>.parquet, '
            f'{len(results) - len(errors)} ORF/CRISPR plates')
    load(cur, 'JUMP_GENETIC_WELL_FEATURES', feats, load_id, note, 'per-plate parquet (selected columns)',
         hashlib.sha256(pd.util.hash_pandas_object(feats, index=False).values.tobytes()).hexdigest(), existing)
    for name, f in META_FILES.items():
        p = meta(f)
        df = pd.read_csv(p)
        if name == 'JUMP_WELL':  # keep only wells of ORF/CRISPR plates (the full table has every compound well)
            df = df[df.Metadata_Plate.isin(set(gen.Metadata_Plate))]
        load(cur, name, df, load_id, f'{META_URL}/{f}', f, hashlib.sha256(p.read_bytes()).hexdigest(), existing)
    print('LOAD_ID', load_id)


if __name__ == '__main__':
    main()
