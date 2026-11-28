"""Export the GlycoEnzOnto subcellular compartment of each glycogene as a TERM2GENE table (user, 2026-11-23:
"経路ごと以外の glycogene の軸", "やってみようか").

Source: RAW.GLYCODB.GLYCOENZONTO_GLYCOGENE, column "Compartment" (comma-separated; a gene can have several rows and
several compartments). Each listed compartment becomes a term; a gene belongs to every compartment listed for it.
Merged labels: 'Endoplasmic Reticulum Lumen' -> 'Endoplasmic Reticulum'; 'cisGolgi', 'medialGolgi', 'transGolgi',
'cisGolgi network', 'transGolgi network' -> 'Golgi'; 'Late Endosome' -> 'Endosome'.
Output: <out>/compartment_term2gene.tsv (term, gene) and compartment_counts.csv.
Usage: python export_compartment_term2gene.py <out_dir>
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cursor  # noqa: E402

MERGE = {'Endoplasmic Reticulum Lumen': 'Endoplasmic Reticulum', 'transGolgi network': 'Golgi', 'cisGolgi': 'Golgi',
         'medialGolgi': 'Golgi', 'transGolgi': 'Golgi', 'cisGolgi network': 'Golgi', 'Late Endosome': 'Endosome'}


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    cur = cursor()
    cur.execute('select distinct "geneName", "Compartment" from RAW.GLYCODB.GLYCOENZONTO_GLYCOGENE')
    rows = []
    for gene, comp in cur.fetchall():
        for c in (comp or '').split(','):
            c = MERGE.get(c.strip(), c.strip())
            if c:
                rows.append((c, gene))
    t2g = pd.DataFrame(rows, columns=['term', 'gene']).drop_duplicates().sort_values(['term', 'gene'])
    t2g.to_csv(out / 'compartment_term2gene.tsv', sep='\t', index=False)
    n = t2g.groupby('term').gene.nunique().sort_values(ascending=False)
    n.rename('n_genes').to_csv(out / 'compartment_counts.csv')
    print(n.to_string())
    print('genes with a compartment:', t2g.gene.nunique())


if __name__ == '__main__':
    main(sys.argv[1])
