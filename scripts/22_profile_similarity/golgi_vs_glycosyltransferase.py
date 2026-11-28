"""Separate "works in the Golgi" from "is a glycosyltransferase" for the two significant compartment results
(OE PC1 (+), KO PC2 (+)) (user, 2026-11-23: "やってみようか").
Tail = top 10% of the PC scores (as in compartment_enrichment.R). Stratified 2 x 2 tables with Fisher's exact test:
within glycosyltransferases, Golgi vs not; within non-glycosyltransferases, Golgi vs not; within Golgi genes,
glycosyltransferase vs not; within non-Golgi genes, glycosyltransferase vs not. BH over the 4 tests of each
perturbation (padj_within_pert) and over all 8 (padj_all8).
Golgi: compartment_term2gene.tsv (GlycoEnzOnto Compartment); glycosyltransferase: enzyme_class_term2gene.tsv
(GO:0016757 and descendants, org.Hs.eg.db).
Usage: python golgi_vs_glycosyltransferase.py <results/22_profile_similarity dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def main(r):
    r = Path(r)
    comp = pd.read_csv(r / 'compartment/compartment_term2gene.tsv', sep='\t')
    cls = pd.read_csv(r / 'enzyme_class/enzyme_class_term2gene.tsv', sep='\t')
    golgi, gt = set(comp.gene[comp.term == 'Golgi']), set(cls.gene[cls.term == 'Glycosyltransferase'])
    rows = []
    for t, pc in [('ORF', 'PC1'), ('CRISPR', 'PC2')]:
        s = pd.read_csv(r / f'pca_glyco/pca_glyco_{t}_scores.csv')
        n10 = int(np.ceil(0.1 * len(s)))
        s['is_tail'] = s.gene.isin(set(s.sort_values(pc, ascending=False).gene[:n10]))
        s['is_golgi'], s['is_gt'] = s.gene.isin(golgi), s.gene.isin(gt)
        for stratum, sub, var in [('glycosyltransferases', s[s.is_gt], 'is_golgi'),
                                  ('non-glycosyltransferases', s[~s.is_gt], 'is_golgi'),
                                  ('Golgi genes', s[s.is_golgi], 'is_gt'), ('non-Golgi genes', s[~s.is_golgi], 'is_gt')]:
            a = (sub[var] & sub.is_tail).sum(); b = (sub[var] & ~sub.is_tail).sum()
            c = (~sub[var] & sub.is_tail).sum(); d = (~sub[var] & ~sub.is_tail).sum()
            orr, p = fisher_exact([[a, b], [c, d]])
            rows.append(dict(perturbation=t, axis=f'{pc} (+)', stratum=stratum,
                             compared='Golgi vs not' if var == 'is_golgi' else 'glycosyltransferase vs not',
                             tail_yes=f'{a}/{a + b}', tail_no=f'{c}/{c + d}', odds_ratio=orr, p=p))
    res = pd.DataFrame(rows)
    res['padj_within_pert'] = res.groupby('perturbation').p.transform(lambda x: multipletests(x, method='fdr_bh')[1])
    res['padj_all8'] = multipletests(res.p, method='fdr_bh')[1]
    res.to_csv(r / 'compartment/golgi_vs_glycosyltransferase.csv', index=False)
    print(res.round(4).to_string(index=False))


if __name__ == '__main__':
    main(sys.argv[1])
