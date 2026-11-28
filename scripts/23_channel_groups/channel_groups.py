"""Which kinds of glycogenes drive the lectin-channel selectivity of Figure 3? (user, 2026-11-23: "fig3 の 5 つのチャンネル
を経路ごとに見ると面白いけ？", "そしたら fig4 を改善できるかな", "はい！").

Input: per-gene effects in the five channels (scripts/17_channel_specificity/channel_effects.py, model 'main';
<ORF|CRISPR>_main_effects.csv: <channel>_delta (Glass's delta vs other genes' wells on the same plates), <channel>_p,
glycogene flag).

A  Pre-specified lectin gene sets (inputs/lectin_gene_sets.csv, fixed on 2026-11-15 for the LINCS study, before any
   JUMP-CP data were looked at): for each gene, the expected direction of ConA or WGA staining when the gene goes up
   (+1 / -1). For overexpression the expected sign of delta is that direction; for knockout its opposite. For each
   perturbation x set (ConA, WGA) x channel: share of the set's genes whose delta has the expected sign; two-sided
   exact binomial test against 50%; mean of delta x expected sign. BH over the 20 tests.
B  Subcellular compartment (GlycoEnzOnto, compartment_term2gene.tsv): for each perturbation x compartment (>= 10
   glycogenes tested) x channel, share of the compartment's glycogenes whose |delta| exceeds the 95th percentile of
   the non-glycogenes in that channel (as Figure 3b); 2 x 2 table against the other glycogenes, Fisher's exact test;
   BH within each perturbation.
Usage: python channel_groups.py <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

R = Path.home() / 'research/glyco_cell_painting/results'
REPO = Path(__file__).resolve().parents[2]
CH = ['ConA', 'WGA', 'DNA', 'RNA', 'Mito']


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    sets = pd.read_csv(REPO / 'inputs/lectin_gene_sets.csv')
    comp = pd.read_csv(R / '22_profile_similarity/compartment/compartment_term2gene.tsv', sep='\t')
    rows_a, rows_b = [], []
    for t in ('ORF', 'CRISPR'):
        e = pd.read_csv(R / f'17_channel_specificity/{t}_main_effects.csv').set_index('gene')
        for lec in ('ConA', 'WGA'):
            s = sets[sets.lectin == lec]
            s = s[s.gene.isin(e.index)]
            exp = s.expected_direction.to_numpy() * (1 if t == 'ORF' else -1)
            for ch in CH:
                dlt = e.loc[s.gene, f'{ch}_delta'].to_numpy()
                ok = np.isfinite(dlt)
                agree = np.sign(dlt[ok]) == exp[ok]
                k, n = int(agree.sum()), int(ok.sum())
                rows_a.append(dict(perturbation=t, gene_set=lec, channel=ch, n=n, n_expected=k, share=k / n,
                                   p=stats.binomtest(k, n, 0.5).pvalue, mean_signed_delta=float(np.mean(dlt[ok] * exp[ok])),
                                   genes_expected=', '.join(s.gene.to_numpy()[ok][agree])))
        g = e[e.glycogene]
        o = e[~e.glycogene]
        groups = {c: set(comp.gene[comp.term == c]) & set(g.index) for c in comp.term.unique()}
        groups = {c: v for c, v in groups.items() if len(v) >= 10}
        tmp = []
        for ch in CH:
            q95 = np.nanpercentile(o[f'{ch}_delta'].abs(), 95)
            above = g[f'{ch}_delta'].abs() > q95
            for c, mem in groups.items():
                inn = g.index.isin(mem)
                a, b = int((above & inn).sum()), int((~above & inn).sum())
                cc, d = int((above & ~inn).sum()), int((~above & ~inn).sum())
                orr, p = stats.fisher_exact([[a, b], [cc, d]])
                tmp.append(dict(perturbation=t, compartment=c, channel=ch, n=a + b, n_above=a, share=a / (a + b),
                                share_other_glyco=cc / (cc + d), odds_ratio=orr, p=p,
                                genes_above=', '.join(sorted(g.index[above & inn]))))
        tb = pd.DataFrame(tmp)
        tb['p_adj'] = multipletests(tb.p, method='fdr_bh')[1]
        rows_b.append(tb)
    a = pd.DataFrame(rows_a)
    a['p_adj'] = multipletests(a.p, method='fdr_bh')[1]
    b = pd.concat(rows_b)
    a.to_csv(out / 'prespecified_sets_by_channel.csv', index=False)
    b.to_csv(out / 'compartment_by_channel.csv', index=False)
    pd.set_option('display.width', 250)
    print('A. Pre-specified sets: share of genes moving in the expected direction')
    print(a.drop(columns='genes_expected').round(3).to_string(index=False))
    print('\nB. Compartment x channel: tests', len(b), '; p < 0.05:', int((b.p < 0.05).sum()), '; adjusted p < 0.05:',
          int((b.p_adj < 0.05).sum()))
    print(b[b.p < 0.05].drop(columns='genes_above').round(3).to_string(index=False))


if __name__ == '__main__':
    main(sys.argv[1])
