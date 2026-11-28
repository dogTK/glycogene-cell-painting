"""Replicate Figure 3 (lectin-channel selectivity of glycogene knockouts) in PERISCOPE (user, 2026-11-23: external
validation with PERISCOPE, "3 かな").

Input: fetch_periscope.py output (~/research/glyco_cell_painting/data/periscope/<cell>_<plate>.parquet): gene-level
profiles per plate, each feature standardised across genes by the PERISCOPE pipeline; columns
Cells_Intensity_MeanIntensity_{ConA, WGA, Phalloidin, DAPI_Painting, Mito}, Cells_AreaShape_Area, Nuclei_AreaShape_Area.
Conditions: A549 (plates CP186A-H, N), HeLa DMEM (CP257A, B, D, F, H), HeLa HPLM (CP257J, K, L, N), as in the
PERISCOPE publication files (HeLa_DMEM_*, HeLa_HPLM_*).
Per condition:
  1. Gene value = median over the plates of the condition (nontargeting row dropped).
  2. Adjustment (as the 'main' model of Figure 3, with Phalloidin now a separate channel): each channel is regressed
     on cell area, nucleus area and the other non-lectin channels (ConA, WGA: DAPI, Mito, Phalloidin; Phalloidin:
     DAPI, Mito; DAPI: Mito, Phalloidin; Mito: DAPI, Phalloidin); the residual is standardised across genes (z).
  3. Glycogenes (GlycoEnzOnto) vs other genes, per channel: two-sample KS test on |z|; share of glycogenes with |z|
     above the other genes' 95th percentile, exact binomial test against 5% and Clopper-Pearson 95% CI.
  4. Concordance with JUMP-CP knockouts (17_channel_specificity/CRISPR_main_effects.csv) for glycogenes in both:
     Spearman correlation of the same channel (ConA-ConA, WGA-WGA; JUMP DNA-DAPI, Mito-Mito).
Outputs: <out>/periscope_gene_effects_<condition>.csv, periscope_channel_selectivity.csv, periscope_vs_jump.csv,
periscope_named_genes.csv.
Usage: python periscope_channels.py <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

D = Path.home() / 'research/glyco_cell_painting/data/periscope'
R = Path.home() / 'research/glyco_cell_painting/results'
COND = {'A549': ('A549', ['CP186A', 'CP186B', 'CP186C', 'CP186D', 'CP186E', 'CP186F', 'CP186G', 'CP186H', 'CP186N']),
        'HeLa_DMEM': ('HeLa', ['CP257A', 'CP257B', 'CP257D', 'CP257F', 'CP257H']),
        'HeLa_HPLM': ('HeLa', ['CP257J', 'CP257K', 'CP257L', 'CP257N'])}
F = {'ConA': 'Cells_Intensity_MeanIntensity_ConA', 'WGA': 'Cells_Intensity_MeanIntensity_WGA',
     'Phalloidin': 'Cells_Intensity_MeanIntensity_Phalloidin', 'DAPI': 'Cells_Intensity_MeanIntensity_DAPI_Painting',
     'Mito': 'Cells_Intensity_MeanIntensity_Mito', 'area': 'Cells_AreaShape_Area', 'nuc_area': 'Nuclei_AreaShape_Area'}
COV = {'ConA': ['DAPI', 'Mito', 'Phalloidin'], 'WGA': ['DAPI', 'Mito', 'Phalloidin'], 'Phalloidin': ['DAPI', 'Mito'],
       'DAPI': ['Mito', 'Phalloidin'], 'Mito': ['DAPI', 'Phalloidin']}
CH = ['ConA', 'WGA', 'Phalloidin', 'DAPI', 'Mito']
NAMED = ['MGAT1', 'MAN2A1', 'MAN1A2', 'SLC35A2', 'B4GALT1', 'B3GNT2', 'GALE', 'HSPA5', 'ST3GAL4']
JUMP_MATCH = {'ConA': 'ConA', 'WGA': 'WGA', 'DAPI': 'DNA', 'Mito': 'Mito'}


def gene_effects(cell, plates):
    x = pd.concat([pd.read_parquet(D / f'{cell}_{p}.parquet') for p in plates])
    x = x.rename(columns={'Metadata_Foci_Barcode_MatchedTo_GeneCode': 'gene'})
    x = x[x.gene != 'nontargeting']
    g = x.groupby('gene')[list(F.values())].median().rename(columns={v: k for k, v in F.items()}).dropna()
    out = pd.DataFrame(index=g.index)
    for ch in CH:
        X = np.column_stack([np.ones(len(g)), g[['area', 'nuc_area'] + COV[ch]].to_numpy()])
        beta, *_ = np.linalg.lstsq(X, g[ch].to_numpy(), rcond=None)
        res = g[ch].to_numpy() - X @ beta
        out[ch] = (res - res.mean()) / res.std(ddof=0)
    out['n_plates'] = x.groupby('gene').size().reindex(out.index)
    return out


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    glyco = set(pd.read_csv(R / '13_jump/nominal_enrichment_infection_filtered/glycoenzonto_term2gene.tsv',
                            sep='\t').gene)
    jump = pd.read_csv(R / '17_channel_specificity/CRISPR_main_effects.csv').set_index('gene')
    sel, conc, named = [], [], []
    for cond, (cell, plates) in COND.items():
        e = gene_effects(cell, plates)
        e['glycogene'] = e.index.isin(glyco)
        e.to_csv(out / f'periscope_gene_effects_{cond}.csv')
        for ch in CH:
            xg, xo = e.loc[e.glycogene, ch].abs(), e.loc[~e.glycogene, ch].abs()
            ks = stats.ks_2samp(xg, xo)
            q95 = np.percentile(xo, 95)
            k, n = int((xg > q95).sum()), len(xg)
            bt = stats.binomtest(k, n, 0.05)
            ci = bt.proportion_ci(confidence_level=0.95, method='exact')
            sel.append(dict(condition=cond, channel=ch, n_glyco=n, n_other=len(xo), ks_D=ks.statistic, ks_p=ks.pvalue,
                            n_above=k, share_above=k / n, ci_lo=ci.low, ci_hi=ci.high, binom_p=bt.pvalue,
                            median_ratio=np.median(xg) / np.median(xo)))
        both = e[e.glycogene].join(jump, how='inner', rsuffix='_jump')
        for ch, jch in JUMP_MATCH.items():
            rho, p = stats.spearmanr(both[ch], both[f'{jch}_delta'], nan_policy='omit')
            conc.append(dict(condition=cond, periscope=ch, jump=jch, n=len(both), spearman_rho=rho, p=p))
        for gname in NAMED:
            if gname in e.index:
                named.append(dict(condition=cond, gene=gname, **{ch: e.loc[gname, ch] for ch in CH}))
    sel, conc, named = pd.DataFrame(sel), pd.DataFrame(conc), pd.DataFrame(named)
    sel.to_csv(out / 'periscope_channel_selectivity.csv', index=False)
    conc.to_csv(out / 'periscope_vs_jump.csv', index=False)
    named.to_csv(out / 'periscope_named_genes.csv', index=False)
    pd.set_option('display.width', 220)
    print(sel.round(3).to_string(index=False))
    print(conc.round(3).to_string(index=False))
    print(named.round(2).to_string(index=False))


if __name__ == '__main__':
    main(sys.argv[1])
