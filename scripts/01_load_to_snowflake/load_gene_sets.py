"""Copy the hand-curated lectin gene sets (inputs/lectin_gene_sets.csv) into RAW.GLYCO_CELL_PAINTING.LECTIN_GENE_SETS.

The git CSV is the source of truth; Snowflake holds a copy per load. Each load is appended with a new _LOAD_ID and the
git commit of the CSV, so earlier versions stay readable. The CSV must be committed (clean) before loading.
Same route as load_cell_painting.py: Parquet -> S3 -> external stage -> COPY INTO.
usage: python load_gene_sets.py
"""
import os
import hashlib
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '00_compound_overlap'))
from fetch_a549_perts import connect_snowflake  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
CSV = REPO / 'inputs/lectin_gene_sets.csv'
SCHEMA = 'RAW.GLYCO_CELL_PAINTING'
TABLE = 'LECTIN_GENE_SETS'
S3 = f"s3://{os.environ.get('GLYCO_S3_BUCKET', '<your-bucket>')}/glyco_cell_painting/raw"
STAGE = os.environ.get('GLYCO_SNOWFLAKE_STAGE', '@<your_stage>') + '/glyco_cell_painting/raw'
EXPORT = Path.home() / 'research/glyco_cell_painting/data/export'
COMMENT = ('Hand-curated glycogenes whose expression is expected to change ConA or WGA staining, with the expected '
           'direction (1: gene up -> staining up, -1: gene up -> staining down). Source of truth is '
           'inputs/lectin_gene_sets.csv in the glyco_cell_painting git repo; one copy per load (_LOAD_ID)')


def git(*args):
    return subprocess.run(['git', '-C', str(REPO), *args], capture_output=True, text=True, check=True).stdout.strip()


def main():
    if git('status', '--porcelain', '--', str(CSV)):
        sys.exit('commit inputs/lectin_gene_sets.csv first')
    commit = git('log', '-1', '--format=%H', '--', str(CSV))
    load_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    d = pd.read_csv(CSV)
    d.columns = [c.upper() for c in d.columns]
    d.insert(0, '_LOAD_ID', load_id)
    out = EXPORT / load_id / TABLE
    out.mkdir(parents=True, exist_ok=True)
    d.to_parquet(out / f'{TABLE}.parquet', index=False)
    subprocess.run(['aws', 's3', 'sync', str(out), f'{S3}/{load_id}/{TABLE}', '--only-show-errors'], check=True)

    cur = connect_snowflake().cursor()
    cur.execute('use secondary roles none')
    loc = f'{STAGE}/{load_id}/{TABLE}/'
    cur.execute(f"SHOW TABLES LIKE '{TABLE}' IN SCHEMA {SCHEMA}")
    if not cur.fetchall():
        cur.execute(f"CREATE TABLE {SCHEMA}.{TABLE} USING TEMPLATE (SELECT ARRAY_AGG(OBJECT_CONSTRUCT(*)) "
                    f"FROM TABLE(INFER_SCHEMA(LOCATION => '{loc}', FILE_FORMAT => '{SCHEMA}.PARQUET_FMT')))")
        cur.execute(f'COMMENT ON TABLE {SCHEMA}.{TABLE} IS $${COMMENT}$$')
    cur.execute(f"COPY INTO {SCHEMA}.{TABLE} FROM '{loc}' FILE_FORMAT = (FORMAT_NAME = '{SCHEMA}.PARQUET_FMT') "
                'MATCH_BY_COLUMN_NAME = CASE_SENSITIVE')
    now = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
    cur.execute(f'INSERT INTO {SCHEMA}._LOAD_MANIFEST VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)',
                (load_id, now, TABLE, len(d), 'https://github.com/dogTK/glyco_cell_painting', commit,
                 'inputs/lectin_gene_sets.csv', hashlib.sha256(CSV.read_bytes()).hexdigest(),
                 'hand-curated; loaded by scripts/01_load_to_snowflake/load_gene_sets.py'))
    cur.execute(f'SELECT COUNT(*) FROM {SCHEMA}.{TABLE} WHERE _LOAD_ID = %s', (load_id,))
    print(TABLE, 'rows', cur.fetchone()[0], 'LOAD_ID', load_id, 'commit', commit[:7])


if __name__ == '__main__':
    main()
