"""List the A549 perturbations in the LINCS L1000 glycogene table (CycleGAN-predicted RNA-seq-like Level 5).

Reads RAW.LINCS.GLYCO_GENES_WIDE in Snowflake (the same table the LINCS glycogene papers use) and writes one row per
perturbation with its name, BRD id, InChIKey and number of signatures.
Usage: python fetch_a549_perts.py <out.csv>   (needs snowflake-connector-python and ~/.snowflake/config.toml)
"""
import os
import sys
import tomllib

import pandas as pd
from cryptography.hazmat.primitives import serialization
from snowflake.connector import connect


def connect_snowflake():
    c = tomllib.load(open(os.path.expanduser('~/.snowflake/config.toml'), 'rb'))['connections']['default']
    key = serialization.load_pem_private_key(open(os.path.expanduser(c['private_key_path']), 'rb').read(),
                                             password=None)
    der = key.private_bytes(serialization.Encoding.DER, serialization.PrivateFormat.PKCS8,
                            serialization.NoEncryption())
    return connect(user=c['user'], account=c['account'], private_key=der, warehouse=c.get('warehouse'),
                   role=c.get('role'), database='RAW', schema='LINCS')


def main(out):
    cur = connect_snowflake().cursor()
    cur.execute('SELECT "pertname", "pertid", "inchi_key", COUNT(*) FROM RAW.LINCS.GLYCO_GENES_WIDE '
                'WHERE "cell" = \'A549\' GROUP BY 1, 2, 3')
    df = pd.DataFrame(cur.fetchall(), columns=['pertname', 'pertid', 'inchi_key', 'n_sig'])
    df.to_csv(out, index=False)
    print('A549 perturbations', len(df), 'signatures', df.n_sig.sum())


if __name__ == '__main__':
    main(sys.argv[1])
