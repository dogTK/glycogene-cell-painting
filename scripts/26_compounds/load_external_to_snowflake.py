"""Load the externally fetched tables of 2026-11-23 into Snowflake (RAW.GLYCO_CELL_PAINTING, append-only), by the same
route as scripts/22_profile_similarity/load_jump_profiles.py (Parquet -> S3 -> external stage -> CREATE TABLE USING
TEMPLATE -> COPY INTO; table comment, _LOAD_ID, COLUMN_NAMES, _LOAD_MANIFEST with sha256).
Tables (each loaded only if its local file exists and the table does not exist yet):
  JUMP_COMPOUND_WELL_FEATURES  data/jump_compound/compound_wells_raw.parquet (fetch_compound_wells.py)
  TSUI2024_HHL_CRISPR_SCREEN   data/tsui2024/tsui2024_genomewide_castle.csv (Supplementary Dataset 2 of Tsui et al. 2024)
  PERISCOPE_GENE_PROFILES      data/periscope/<cell>_<plate>.parquet concatenated (fetch_periscope.py), when all 18
                               plates are present
Usage: python load_external_to_snowflake.py
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

D = Path.home() / 'research/glyco_cell_painting/data'
SCHEMA = 'RAW.GLYCO_CELL_PAINTING'
S3 = f"s3://{os.environ.get('GLYCO_S3_BUCKET', '<your-bucket>')}/glyco_cell_painting/raw"
STAGE = os.environ.get('GLYCO_SNOWFLAKE_STAGE', '@<your_stage>') + '/glyco_cell_painting/raw'
EXPORT = D / 'export'


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load(cur, name, df, comment, source_url, version, source_file, sha, load_id):
    orig = [str(c) for c in df.columns]
    df = df.copy()
    df.columns = [c.upper().replace(' ', '_').replace('%', 'PCT').replace('-', '_').replace('#', 'N_') for c in orig]
    for c in df.select_dtypes(['object', 'category']).columns:
        df[c] = df[c].astype('string')
    df.insert(0, '_LOAD_ID', load_id)
    out = EXPORT / load_id / name
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
    cno = EXPORT / load_id / f'COLUMN_NAMES_{name}'
    cno.mkdir(parents=True, exist_ok=True)
    cn.to_parquet(cno / 'COLUMN_NAMES.parquet', index=False)
    subprocess.run(['aws', 's3', 'sync', str(cno), f'{S3}/{load_id}/COLUMN_NAMES_{name}', '--only-show-errors'],
                   check=True)
    cur.execute(f"COPY INTO {SCHEMA}.COLUMN_NAMES FROM '{STAGE}/{load_id}/COLUMN_NAMES_{name}/' "
                f"FILE_FORMAT = (FORMAT_NAME = '{SCHEMA}.PARQUET_FMT') MATCH_BY_COLUMN_NAME = CASE_SENSITIVE")
    now = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
    cur.execute(f'INSERT INTO {SCHEMA}._LOAD_MANIFEST VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)',
                (load_id, now, name, len(df), source_url, version, source_file, sha,
                 'loaded by glyco_cell_painting scripts/26_compounds/load_external_to_snowflake.py'))
    cur.execute(f'SELECT COUNT(*) FROM {SCHEMA}.{name} WHERE _LOAD_ID = %s', (load_id,))
    print(f'{name}: {len(df):,} rows; in Snowflake {cur.fetchone()[0]:,}', flush=True)


def main():
    cur = cursor()
    cur.execute(f'SHOW TABLES IN SCHEMA {SCHEMA}')
    existing = {r[1] for r in cur.fetchall()}
    load_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    f = D / 'jump_compound/compound_wells_raw.parquet'
    if f.exists() and 'JUMP_COMPOUND_WELL_FEATURES' not in existing:
        load(cur, 'JUMP_COMPOUND_WELL_FEATURES', pd.read_parquet(f),
             'JUMP-CP (cpg0016) compound plates (10 sources, 1,716 plates): raw well-level CellProfiler aggregates for the '
             'columns this study uses (cell mean intensity of ER, AGP, RNA, DNA, Mito; cell area; '
             'Nuclei_Number_Object_Number) with plate, well, JCP2022 id and perturbation type. Column projection of '
             'profiles_assembled/COMPOUND/v1.0a/profiles.parquet. One row per well',
             's3://cellpainting-gallery/cpg0016-jump/source_all/workspace/profiles_assembled/COMPOUND/v1.0a/profiles.parquet',
             'profiles_assembled v1.0a', 'profiles.parquet (selected columns)', sha256_file(f), load_id)
    f = D / 'tsui2024/tsui2024_genomewide_castle.csv'
    x = D / 'tsui2024/supp_data2.xlsx'
    if f.exists() and 'TSUI2024_HHL_CRISPR_SCREEN' not in existing:
        load(cur, 'TSUI2024_HHL_CRISPR_SCREEN', pd.read_csv(f),
             'Tsui et al. 2024 Nat Commun (doi:10.1038/s41467-024-53225-1) Supplementary Dataset 2: genome-wide CRISPR '
             'knockout screen in A549 sorted by binding of the high-mannose lectin HHL (top / bottom 25%), casTLE '
             'results per gene (Combo casTLE Effect > 0 = knockout increases HHL binding; 10%FDR = 1 for hits)',
             'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-53225-1/MediaObjects/41467_2024_53225_MOESM4_ESM.xlsx',
             'published 2024-11', '41467_2024_53225_MOESM4_ESM.xlsx', sha256_file(x), load_id)
    files = sorted((D / 'periscope').glob('*_CP*.parquet'))
    if len(files) == 18 and 'PERISCOPE_GENE_PROFILES' not in existing:
        df = pd.concat([pd.read_parquet(p) for p in files], ignore_index=True)
        man = pd.read_csv(D / 'periscope/manifest.csv')
        load(cur, 'PERISCOPE_GENE_PROFILES', df,
             'PERISCOPE (cpg0021, Ramezani et al. 2025): genome-wide CRISPR knockout optical pooled screens in A549 and '
             'HeLa (DMEM: CP257A, B, D, F, H; HPLM: CP257J, K, L, N); gene-level normalized profiles per plate, columns '
             'kept: cell mean intensity of ConA, WGA, Phalloidin, DAPI, Mito; cell and nucleus area. One row per gene '
             'per plate',
             's3://cellpainting-gallery/cpg0021-periscope/broad/workspace/profiles/<A549|HeLa>/*_gene_normalized_*',
             'cpg0021 profiles', '18 gene_normalized files (selected columns)',
             hashlib.sha256(''.join(man.sort_values(['cell', 'plate']).sha256).encode()).hexdigest(), load_id)
    print('LOAD_ID', load_id)


if __name__ == '__main__':
    main()
