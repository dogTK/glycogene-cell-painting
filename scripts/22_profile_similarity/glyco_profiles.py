"""Write the gene-level Cell Painting profiles of the GlycoEnzOnto genes for the PCA in pca_glyco.R
(user, 2026-11-23: PCA laid out as in the first paper, glycogenes only, "とりま 1 で").

Profiles: JUMP-CP profiles_assembled v1.0a ..._featselect (named CellProfiler features; also in
RAW.GLYCO_CELL_PAINTING.JUMP_<ORF|CRISPR>_PROFILES_FEATSELECT). ORF wells listed as low infection efficiency /
outliers are dropped. Gene profile = median over all wells of the gene (all reagents). Only GlycoEnzOnto genes.
Output: <out>/glyco_profiles_<ORF|CRISPR>.csv (rows genes, columns features).
Usage: python glyco_profiles.py <out_dir>
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from control_modules import D, P, R, low_infection_mask  # noqa: E402


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    glyco = set(pd.read_csv(R / '13_jump/nominal_enrichment_infection_filtered/glycoenzonto_term2gene.tsv',
                            sep='\t').gene)
    for t in ('ORF', 'CRISPR'):
        df = pd.read_parquet(P / f'{t}/profiles_wellpos_cc_var_mad_outlier_featselect.parquet')
        if t == 'ORF':
            df = df[~low_infection_mask(df)]
        meta = pd.read_csv(D / f'jump_meta/{t.lower()}.csv.gz')
        df = df.merge(meta[['Metadata_JCP2022', 'Metadata_Symbol']], on='Metadata_JCP2022')
        df = df[df.Metadata_Symbol.isin(glyco)]
        feats = [c for c in df.columns if not c.startswith('Metadata')]
        g = df.groupby('Metadata_Symbol')[feats].median().dropna(axis=1, how='any')
        g.index.name = 'gene'
        g.to_csv(out / f'glyco_profiles_{t}.csv')
        print(f'{t}: {g.shape[0]} glycogenes x {g.shape[1]} features')


if __name__ == '__main__':
    main(sys.argv[1])
