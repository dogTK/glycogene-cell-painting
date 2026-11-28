"""Choose the genes shown in the image panels by fixed rules (user, 2026-11-23: "ちゃんとロジックがある遺伝子選びたい",
then "普通fig3じゃない？" -> one rule for both, both shown in Figure 3c).

Effects: Glass's delta and p of the main analysis (results/17_channel_specificity/<ORF|CRISPR>_main_effects.csv;
channel ER = ConA column, AGP = WGA column). Candidates: GlycoEnzOnto glycogenes with p < 0.05 in the channel.
Slots (as in the original request): AGP up, AGP down, ER up.
  Rule (panels 'oe' = ORF, 'ko' = CRISPR): in each slot, the candidate with the largest |delta| in that direction.
  For information only (not used for the choice): PERISCOPE WGA / ConA z (periscope_channels.py, adjusted) in the
  three conditions and whether the sign matches JUMP in all three; Tsui et al. 2024 hit (tsui_comparison.py).
Output: <out>/gene_selection.csv (every candidate per panel and slot, with the chosen one flagged).
Usage: python select_genes.py <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path.home() / 'research/glyco_cell_painting/results'
SLOTS = [('AGP', 1), ('AGP', -1), ('ER', 1)]
COL = {'AGP': 'WGA', 'ER': 'ConA'}
CONDS = ['A549', 'HeLa_DMEM', 'HeLa_HPLM']


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    gl = set(pd.read_csv(R / '13_jump/nominal_enrichment_infection_filtered/glycoenzonto_term2gene.tsv', sep='\t').gene)
    per = {c: pd.read_csv(R / f'24_external/periscope/periscope_gene_effects_{c}.csv').set_index('gene') for c in CONDS}
    tsui = set(pd.read_csv(R / '24_external/tsui/tsui_hits_in_jump.csv').iloc[:, 0])
    rows = []
    for panel, mod in (('oe', 'ORF'), ('ko', 'CRISPR')):
        e = pd.read_csv(R / f'17_channel_specificity/{mod}_main_effects.csv')
        e = e[e.gene.isin(gl)]
        for ch, sign in SLOTS:
            c = COL[ch]
            cand = e[(e[f'{c}_p'] < 0.05) & (np.sign(e[f'{c}_delta']) == sign)].copy()
            cand['panel'], cand['modality'], cand['channel'] = panel, mod, ch
            cand['direction'] = 'up' if sign > 0 else 'down'
            cand['delta'], cand['p'] = cand[f'{c}_delta'], cand[f'{c}_p']
            for k in CONDS:
                cand[f'periscope_{c}_{k}'] = cand.gene.map(per[k][c])
            same = np.ones(len(cand), bool)
            for k in CONDS:
                same &= (np.sign(cand[f'periscope_{c}_{k}']) == sign).to_numpy()
            cand['periscope_same_sign_all3'] = same
            cand['tsui_hit'] = cand.gene.isin(tsui)
            cand['eligible'] = True
            el = cand[cand.eligible]
            pick = el.loc[el.delta.abs().idxmax(), 'gene'] if len(el) else None
            cand['chosen'] = cand.gene.eq(pick)
            rows.append(cand.sort_values('delta', key=abs, ascending=False))
    keep = ['panel', 'modality', 'channel', 'direction', 'gene', 'delta', 'p'] + \
        [f'periscope_{c}_{k}' for c in ('WGA', 'ConA') for k in CONDS] + \
        ['periscope_same_sign_all3', 'tsui_hit', 'eligible', 'chosen']
    s = pd.concat(rows)
    s = s.reindex(columns=keep)
    s.to_csv(out / 'gene_selection.csv', index=False)
    pd.set_option('display.width', 250)
    print(s[['panel', 'channel', 'direction', 'gene', 'delta', 'p', 'periscope_same_sign_all3', 'tsui_hit',
             'eligible', 'chosen']].round(3).to_string(index=False))


if __name__ == '__main__':
    main(sys.argv[1])
