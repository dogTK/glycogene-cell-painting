"""Compare the JUMP-CP knockout effects with an independent lectin-binding CRISPR screen (user, 2026-11-23: external
validation, "tsui もやっといて").

External data: Tsui et al. 2024 Nat Commun (doi:10.1038/s41467-024-53225-1), Supplementary Dataset 2
(41467_2024_53225_MOESM4_ESM.xlsx from static-content.springer.com; sha256 6e58a61a...): genome-wide CRISPR knockout
in A549, cells sorted by binding of the high-mannose lectin HHL (top / bottom 25%), casTLE. Combo casTLE Effect > 0
= knockout increases HHL binding (MGAT1 +3.1, MAN1A2 +2.4); 10%FDR = 1 for 107 genes.
JUMP-CP: per-gene knockout effects in five channels (scripts/17_channel_specificity, CRISPR_main_effects.csv).
For genes present in both (all genes, and GlycoEnzOnto genes only):
  1. Spearman correlation between the casTLE effect and the JUMP delta of each channel (ConA is the matching channel:
     ConA and HHL both bind high-mannose N-glycans; the other channels are the comparison).
  2. 2 x 2 table: Tsui hit with positive effect (10% FDR, effect > 0) vs JUMP channel increased (p < 0.05, delta > 0);
     Fisher's exact test (two-sided).
Usage: python tsui_comparison.py <tsui csv (tsui2024_genomewide_castle.csv)> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = Path.home() / 'research/glyco_cell_painting/results'
CH = ['ConA', 'WGA', 'DNA', 'RNA', 'Mito']


def main(tsui_csv, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(tsui_csv).rename(columns={'Combo casTLE Effect': 'effect', '10%FDR': 'hit10'})
    t = t.dropna(subset=['Symbol']).drop_duplicates('Symbol').set_index('Symbol')
    j = pd.read_csv(R / '17_channel_specificity/CRISPR_main_effects.csv').set_index('gene')
    both = j.join(t[['effect', 'hit10', 'casTLE p-value']], how='inner')
    both['tsui_up'] = (both.hit10 == 1) & (both.effect > 0)
    rows = []
    for scope, d in [('all genes', both), ('GlycoEnzOnto genes', both[both.glycogene])]:
        for ch in CH:
            rho, p = stats.spearmanr(d.effect, d[f'{ch}_delta'], nan_policy='omit')
            up = (d[f'{ch}_p'] < 0.05) & (d[f'{ch}_delta'] > 0)
            a = int((d.tsui_up & up).sum()); b = int((d.tsui_up & ~up).sum())
            c = int((~d.tsui_up & up).sum()); e = int((~d.tsui_up & ~up).sum())
            orr, pf = stats.fisher_exact([[a, b], [c, e]])
            rows.append(dict(scope=scope, channel=ch, n=len(d), spearman_rho=rho, spearman_p=p, tsui_hits=a + b,
                             hits_with_jump_up=a, share_hits_up=a / max(a + b, 1), share_nonhits_up=c / max(c + e, 1),
                             odds_ratio=orr, fisher_p=pf,
                             genes=', '.join(sorted(d.index[d.tsui_up & up]))))
    res = pd.DataFrame(rows)
    res.to_csv(out / 'tsui_vs_jump_ko.csv', index=False)
    hits = both[both.tsui_up].sort_values('effect', ascending=False)
    hits[['effect', 'casTLE p-value', 'glycogene'] + [f'{c}_delta' for c in CH] + [f'{c}_p' for c in CH]].to_csv(
        out / 'tsui_hits_in_jump.csv')
    pd.set_option('display.width', 250)
    print(res.drop(columns='genes').round(4).to_string(index=False))
    print('\nTsui positive hits present in JUMP KO:', len(hits), '(glycogenes', int(hits.glycogene.sum()), ')')
    print(hits[['effect', 'glycogene', 'ConA_delta', 'ConA_p', 'WGA_delta', 'DNA_delta']].round(3).head(25).to_string())


if __name__ == '__main__':
    main(*sys.argv[1:3])
