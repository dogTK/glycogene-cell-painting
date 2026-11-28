"""Write the GlycoEnzOnto pathway gene sets as a clusterProfiler TERM2GENE table (term, gene), from
RAW.GLYCODB.GLYCOENZONTO_PATHWAY_GMT. Usage: python export_term2gene.py <out.tsv>"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cursor  # noqa: E402

with cursor() as cur:
    cur.execute('select PATHWAY_NAME, GENE_SYMBOLS from RAW.GLYCODB.GLYCOENZONTO_PATHWAY_GMT')
    rows = [(n, g) for n, s in cur.fetchall() for g in s.split('|')]
t = pd.DataFrame(rows, columns=['term', 'gene']).drop_duplicates()
Path(sys.argv[1]).parent.mkdir(parents=True, exist_ok=True)
t.to_csv(sys.argv[1], sep='\t', index=False)
print(f'{t.term.nunique()} pathways, {t.gene.nunique()} genes -> {sys.argv[1]}')
