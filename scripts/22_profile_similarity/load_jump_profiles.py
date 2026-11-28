"""Load the JUMP-CP assembled ORF / CRISPR well-level profiles (Cell Painting Gallery, public) into Snowflake
(user, 2026-11-23: "全特徴量落としていいよ", "snowflake 入れとくのお忘れなくね").

Source (downloaded with `aws s3 cp --no-sign-request` into ~/research/glyco_cell_painting/data/jump_profiles_v1.0a/):
  s3://cellpainting-gallery/cpg0016-jump/source_all/workspace/profiles_assembled/<ORF|CRISPR>/v1.0a/
  - profiles_wellpos_cc_var_mad_outlier_featselect.parquet           (named CellProfiler features after
    well-position and cell-count correction, variance filter, MAD robust scaling, outlier removal, feature selection)
  - ORF:    ..._featselect_sphering_harmony.parquet                    (the Morphmap ORF profiles)
  - CRISPR: ..._featselect_sphering_harmony_PCA_corrected.parquet      (the Morphmap CRISPR profiles)
Tables (RAW.GLYCO_CELL_PAINTING, append-only; stops if a table exists):
  JUMP_ORF_PROFILES_FEATSELECT, JUMP_CRISPR_PROFILES_FEATSELECT, JUMP_ORF_PROFILES_HARMONY,
  JUMP_CRISPR_PROFILES_HARMONY_PCA_CORRECTED
Same route as scripts/13_jump/fetch_jump_genetic.py: Parquet -> S3 -> external stage -> CREATE TABLE USING TEMPLATE ->
COPY INTO, with a table comment, _LOAD_ID, COLUMN_NAMES (original case) and _LOAD_MANIFEST (file sha256).
Usage: python load_jump_profiles.py
"""
import os
import hashlib
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cursor  # noqa: E402

DATA = Path.home() / 'research/glyco_cell_painting/data/jump_profiles_v1.0a'
SCHEMA = 'RAW.GLYCO_CELL_PAINTING'
S3 = f"s3://{os.environ.get('GLYCO_S3_BUCKET', '<your-bucket>')}/glyco_cell_painting/raw"
STAGE = os.environ.get('GLYCO_SNOWFLAKE_STAGE', '@<your_stage>') + '/glyco_cell_painting/raw'
SRC = 's3://cellpainting-gallery/cpg0016-jump/source_all/workspace/profiles_assembled'
VERSION = 'profiles_assembled v1.0a'
FS = 'profiles_wellpos_cc_var_mad_outlier_featselect'
TABLES = {
    'JUMP_ORF_PROFILES_FEATSELECT': ('ORF', f'{FS}.parquet',
        'JUMP-CP (cpg0016) ORF plates, U2OS: well-level profiles, named CellProfiler features after well-position and '
        'cell-count correction, variance filter, MAD robust scaling, outlier removal and feature selection '
        '(Cell Painting Gallery profiles_assembled/ORF/v1.0a). One row per well'),
    'JUMP_CRISPR_PROFILES_FEATSELECT': ('CRISPR', f'{FS}.parquet',
        'JUMP-CP (cpg0016) CRISPR plates, U2OS: well-level profiles, named CellProfiler features after well-position '
        'and cell-count correction, variance filter, MAD robust scaling, outlier removal and feature selection '
        '(Cell Painting Gallery profiles_assembled/CRISPR/v1.0a). One row per well'),
    'JUMP_ORF_PROFILES_HARMONY': ('ORF', f'{FS}_sphering_harmony.parquet',
        'JUMP-CP (cpg0016) ORF plates, U2OS: well-level profiles after feature selection, sphering and Harmony batch '
        'correction (the profiles of the Morphmap paper; columns no longer map to single CellProfiler features) '
        '(Cell Painting Gallery profiles_assembled/ORF/v1.0a). One row per well'),
    'JUMP_CRISPR_PROFILES_HARMONY_PCA_CORRECTED': ('CRISPR', f'{FS}_sphering_harmony_PCA_corrected.parquet',
        'JUMP-CP (cpg0016) CRISPR plates, U2OS: well-level profiles after feature selection, sphering, Harmony, PCA '
        'and chromosome-arm correction (the profiles of the Morphmap paper; columns are components, not CellProfiler '
        'features) (Cell Painting Gallery profiles_assembled/CRISPR/v1.0a). One row per well'),
}


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load(cur, name, mod, fname, comment, load_id):
    src = DATA / mod / fname
    df = pd.read_parquet(src)
    orig = [str(c) for c in df.columns]
    df.columns = [c.upper() for c in orig]
    for c in df.select_dtypes('object').columns:
        df[c] = df[c].astype('string')
    df.insert(0, '_LOAD_ID', load_id)
    out = DATA / 'export' / load_id / name
    out.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out / f'{name}.parquet', index=False)
    subprocess.run(['aws', 's3', 'sync', str(out), f'{S3}/{load_id}/{name}', '--only-show-errors'], check=True)
    loc = f'{STAGE}/{load_id}/{name}/'
    cur.execute(f"CREATE TABLE {SCHEMA}.{name} USING TEMPLATE (SELECT ARRAY_AGG(OBJECT_CONSTRUCT(*)) "
                f"FROM TABLE(INFER_SCHEMA(LOCATION => '{loc}', FILE_FORMAT => '{SCHEMA}.PARQUET_FMT')))")
    cur.execute(f'COMMENT ON TABLE {SCHEMA}.{name} IS $${comment}$$')
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
                (load_id, now, name, len(df), f'{SRC}/{mod}/v1.0a/{fname}', VERSION, fname, sha256(src),
                 'loaded by glyco_cell_painting scripts/22_profile_similarity/load_jump_profiles.py'))
    cur.execute(f'SELECT COUNT(*) FROM {SCHEMA}.{name} WHERE _LOAD_ID = %s', (load_id,))
    print(f'{name}: {len(df):,} rows x {len(orig)} columns; in Snowflake {cur.fetchone()[0]:,}', flush=True)


def main():
    cur = cursor()
    cur.execute(f'SHOW TABLES IN SCHEMA {SCHEMA}')
    clash = {r[1] for r in cur.fetchall()} & set(TABLES)
    if clash:
        sys.exit(f'RAW is append-only; tables already exist: {sorted(clash)}')
    load_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    for name, (mod, fname, comment) in TABLES.items():
        load(cur, name, mod, fname, comment, load_id)
    print('LOAD_ID', load_id)


if __name__ == '__main__':
    main()
