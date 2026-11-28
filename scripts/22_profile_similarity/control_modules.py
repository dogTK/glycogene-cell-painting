"""Do glycogenes of the same machinery have similar Cell Painting profiles? Check on modules that must cohere
(user, 2026-11-23: "やってみよう 確認から", "お願いします").

Profiles (JUMP-CP cpg0016 profiles_assembled v1.0a, also in RAW.GLYCO_CELL_PAINTING JUMP_*_PROFILES_*):
  main    : Morphmap profiles (ORF ..._sphering_harmony, CRISPR ..._sphering_harmony_PCA_corrected), all columns
  lectin  : named features (..._featselect) whose channel tokens are only ER and/or AGP (e.g. Cells_Intensity_*_ER,
            Correlation_AGP_ER), i.e. features that do not involve DNA, RNA or Mito
ORF wells listed as low infection efficiency / outliers by the JUMP genetic paper are dropped (as in the main tests).
Gene profile: median over all wells of the gene (all reagents). Similarity: cosine between gene profiles.
Modules (controls, chosen from the biology before looking at similarity):
  GPI   : GlycoEnzOnto 'glycosylphosphatidylinositol anchor biosynthetic pathway'
  LLO   : lipid-linked oligosaccharide assembly on dolichol: ALG*, DPAGT1, RFT1 (GlycoEnzOnto genes)
  OST   : oligosaccharyltransferase complex subunits STT3A, STT3B, DDOST, RPN1, RPN2, DAD1, OSTC, OST4, TMEM258,
          MAGT1, TUSC3, KRTCAP2 (half are not in GlycoEnzOnto; kept because the complex is the control)
Statistic: median pairwise cosine within the module; null = 1,000 random sets of the same size drawn from the
GlycoEnzOnto genes of the same profile set; one-sided p = (1 + #null >= observed) / 1001; also the share of null
medians below the observed (percentile).
Gene sets: all genes present (main) and only genes with a phenotype in Morphmap (any reagent active at FDR 0.1).
Usage: python control_modules.py <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

D = Path.home() / 'research/glyco_cell_painting/data'
P = D / 'jump_profiles_v1.0a'
R = Path.home() / 'research/glyco_cell_painting/results'
FILES = {('ORF', 'main'): 'ORF/profiles_wellpos_cc_var_mad_outlier_featselect_sphering_harmony.parquet',
         ('CRISPR', 'main'): 'CRISPR/profiles_wellpos_cc_var_mad_outlier_featselect_sphering_harmony_PCA_corrected.parquet',
         ('ORF', 'lectin'): 'ORF/profiles_wellpos_cc_var_mad_outlier_featselect.parquet',
         ('CRISPR', 'lectin'): 'CRISPR/profiles_wellpos_cc_var_mad_outlier_featselect.parquet'}
ACT = {'ORF': 'phenotypic-activity-wellpos_cc_var_mad_outlier_featselect_sphering_harmony.csv.gz',
       'CRISPR': 'phenotypic-activity-wellpos_cc_var_mad_outlier_featselect_sphering_harmony_PCA_corrected.csv.gz'}
CH = {'DNA', 'RNA', 'Mito', 'ER', 'AGP'}
OST = ['STT3A', 'STT3B', 'DDOST', 'RPN1', 'RPN2', 'DAD1', 'OSTC', 'OST4', 'TMEM258', 'MAGT1', 'TUSC3', 'KRTCAP2']
NPERM = 1000


def lectin_cols(cols):
    keep = []
    for c in cols:
        tok = set(c.split('_')) & CH
        if tok and tok <= {'ER', 'AGP'}:
            keep.append(c)
    return keep


def low_infection_mask(df):
    meta = pd.read_csv(D / 'morphmap_infection/experiment-metadata.tsv', sep='\t')
    pmap = dict(zip(meta.Assay_Plate_Barcode, meta.Plate_Map_Name))
    bad = pd.read_csv(D / 'morphmap_infection/orf-reagents-low-infection-efficiency-and-outliers.csv.gz')
    bad = set(zip(bad.Metadata_plate_map_name, bad.Metadata_Well, bad.Metadata_broad_sample))
    orf = pd.read_csv(D / 'jump_meta/orf.csv.gz')
    bs = dict(zip(orf.Metadata_JCP2022, orf.Metadata_broad_sample))
    key = zip(df.Metadata_Plate.map(pmap), df.Metadata_Well, df.Metadata_JCP2022.map(bs))
    return np.array([k in bad for k in key])


def gene_profiles(t, kind):
    df = pd.read_parquet(P / FILES[(t, kind)])
    feats = [c for c in df.columns if not c.startswith('Metadata')]
    if kind == 'lectin':
        feats = lectin_cols(feats)
    if t == 'ORF':
        m = low_infection_mask(df)
        print(f'  {t} {kind}: low infection / outlier wells dropped {m.sum()}', flush=True)
        df = df[~m]
    meta = pd.read_csv(D / f'jump_meta/{t.lower()}.csv.gz')
    df = df.merge(meta[['Metadata_JCP2022', 'Metadata_Symbol']], on='Metadata_JCP2022')
    df = df[df.Metadata_Symbol.notna()]
    g = df.groupby('Metadata_Symbol')[feats].median()
    g = g.dropna(axis=1, how='any')
    return g, len(feats)


def active_genes(t):
    a = pd.read_csv(D / 'morphmap_phenotypic_activity' / ACT[t])
    meta = pd.read_csv(D / f'jump_meta/{t.lower()}.csv.gz')
    a = a.merge(meta[['Metadata_JCP2022', 'Metadata_Symbol']], on='Metadata_JCP2022')
    return set(a.loc[a.below_corrected_p, 'Metadata_Symbol'])


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    t2g = pd.read_csv(R / '13_jump/nominal_enrichment_infection_filtered/glycoenzonto_term2gene.tsv', sep='\t')
    glyco = set(t2g.gene)
    gpi = set(t2g.gene[t2g.term == 'glycosylphosphatidylinositol anchor biosynthetic pathway'])
    llo = {g for g in glyco if g.startswith('ALG') or g in ('DPAGT1', 'RFT1')}
    modules = {'GPI': gpi, 'LLO': llo, 'OST': set(OST)}
    rows, pairs = [], []
    rng = np.random.default_rng(20261123)
    for t in ('ORF', 'CRISPR'):
        act = active_genes(t)
        for kind in ('main', 'lectin'):
            prof, nfeat = gene_profiles(t, kind)
            x = prof.to_numpy(float)
            x = x / np.linalg.norm(x, axis=1, keepdims=True)
            idx = {g: i for i, g in enumerate(prof.index)}
            for subset in ('all', 'active'):
                ok = (lambda g: g in idx) if subset == 'all' else (lambda g: g in idx and g in act)
                bg = [idx[g] for g in sorted(glyco) if ok(g)]
                for name, genes in modules.items():
                    mem = sorted(g for g in genes if ok(g))
                    k = len(mem)
                    if k < 3:
                        rows.append(dict(perturbation=t, profile=kind, genes=subset, module=name, n=k))
                        continue
                    mi = [idx[g] for g in mem]
                    s = x[mi] @ x[mi].T
                    iu = np.triu_indices(k, 1)
                    obs = np.median(s[iu])
                    null = np.empty(NPERM)
                    for b in range(NPERM):
                        r = rng.choice(bg, k, replace=False)
                        null[b] = np.median((x[r] @ x[r].T)[iu])
                    rows.append(dict(perturbation=t, profile=kind, genes=subset, module=name, n=k,
                                     n_features=x.shape[1], n_background=len(bg), median_cos=obs,
                                     null_median=np.median(null), null_p95=np.percentile(null, 95),
                                     percentile=(null < obs).mean() * 100, p=(1 + (null >= obs).sum()) / (NPERM + 1),
                                     members=', '.join(mem)))
                    if subset == 'all':
                        for i, j in zip(*iu):
                            pairs.append(dict(perturbation=t, profile=kind, module=name, gene1=mem[i], gene2=mem[j],
                                              cos=s[i, j]))
    res = pd.DataFrame(rows)
    res.to_csv(out / 'control_modules.csv', index=False)
    pd.DataFrame(pairs).to_csv(out / 'control_module_pairs.csv', index=False)
    pd.set_option('display.width', 250)
    print(res.drop(columns='members').round(3).to_string(index=False))


if __name__ == '__main__':
    main(sys.argv[1])
