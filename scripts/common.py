"""Shared helpers: Snowflake connection, L1000 glycogene rank scores, compound matching, lectin gene sets."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / '00_compound_overlap'))
from fetch_a549_perts import connect_snowflake  # noqa: E402,F401

L1000_META = {'canonical_smiles', 'cell', 'cmapid', 'compound_alias', 'dose', 'inchi_key', 'pertid', 'pertname',
              'timepoint', 'sample_id', 'VALUE'}  # VALUE: a VARIANT column, not a gene


def cursor():
    cur = connect_snowflake().cursor()
    cur.execute('use secondary roles none')
    return cur


def l1000_rank_scores(cur, cell='A549', timepoint='24 h', dose='10 uM'):
    """Second-paper gene score: rank glycogenes within each signature (ties averaged; signatures with < 50 valid
    glycogenes dropped), centre by (G + 1) / 2, and take the median over each compound's signatures.
    Returns one row per pertid with ik14 and one column per glycogene (positive = up)."""
    cur.execute("select column_name from RAW.INFORMATION_SCHEMA.COLUMNS where table_schema = 'LINCS' and "
                "table_name = 'GLYCO_GENES_WIDE'")
    genes = sorted(r[0] for r in cur.fetchall() if r[0] not in L1000_META)
    cols = ', '.join(f'try_to_double("{g}") as "{g}"' for g in genes)
    cur.execute(f'select "pertid", left("inchi_key", 14), {cols} from RAW.LINCS.GLYCO_GENES_WIDE '
                'where "cell" = %s and "timepoint" = %s and "dose" = %s', (cell, timepoint, dose))
    sig = pd.DataFrame(cur.fetchall(), columns=['pertid', 'ik14'] + genes)
    vals = sig[genes].astype(float)
    ok = vals.notna().sum(axis=1) >= 50
    ranks = vals[ok].rank(axis=1)
    ranks = ranks.sub((ranks.notna().sum(axis=1) + 1) / 2, axis=0)
    ranks[['pertid', 'ik14']] = sig.loc[ok, ['pertid', 'ik14']]
    scores = ranks.groupby('pertid').agg({**{g: 'median' for g in genes}, 'ik14': 'first'}).reset_index()
    return scores, genes, len(sig), int(ok.sum())


def cp_compounds(cur, dose_recode=6):
    """Cell Painting compounds with a profile at the given dose (6 = 10 uM), with name and MOA."""
    cur.execute('''select distinct m.BROAD_ID, m.INCHIKEY14, m.PERT_INAME, m.MOA, m.TARGET
                   from RAW.GLYCO_CELL_PAINTING.CP_A549_CONSENSUS_MODZ_DMSO p
                   join RAW.GLYCO_CELL_PAINTING.CP_REPURPOSING_MOA m on p.METADATA_BROAD_SAMPLE = m.BROAD_SAMPLE
                   where p.METADATA_DOSE_RECODE = %s''', (dose_recode,))
    return pd.DataFrame(cur.fetchall(), columns=['broad_id', 'ik14', 'pert_iname', 'moa', 'target'])


def match(cp, l1000, genes):
    """Match Cell Painting compounds to L1000 by BRD id, then by InChIKey14; one row per broad_id."""
    a = cp.merge(l1000.drop(columns='ik14'), left_on='broad_id', right_on='pertid')
    b = cp.merge(l1000.drop(columns='pertid'), on='ik14')
    m = pd.concat([a, b]).drop_duplicates('broad_id')
    info = m.set_index('broad_id')[['pert_iname', 'moa', 'target']]
    return info, m.set_index('broad_id')[genes].astype(float)


def lectin_gene_sets(cur):
    """Latest load of the hand-curated lectin gene sets."""
    cur.execute('''select LECTIN, "GROUP", GENE, EXPECTED_DIRECTION from RAW.GLYCO_CELL_PAINTING.LECTIN_GENE_SETS
                   where _LOAD_ID = (select max(_LOAD_ID) from RAW.GLYCO_CELL_PAINTING.LECTIN_GENE_SETS)''')
    return pd.DataFrame(cur.fetchall(), columns=['lectin', 'group', 'gene', 'expected_direction'])
