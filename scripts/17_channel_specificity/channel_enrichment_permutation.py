"""Direct test of whether the glycogene enrichment among large responses is greater in the lectin-containing channels
(ER = ConA, AGP = WGA) than in the other channels (DNA, RNA, Mito) (user, 2026-11-24).

Input (existing outputs of channel_effects.py, not recomputed): results/17_channel_specificity/<ORF|CRISPR>_<model>_effects.csv,
gene-level Glass's delta of the 5 channels; model 'main' (the main analysis: QC, low-infection ORF wells removed,
well-position correction, ConA / WGA adjusted for cell area, RNA, DNA, Mito, cell count; DNA / RNA / Mito for cell
area, cell count and the other two non-lectin channels) and 'symmetric' (sensitivity analysis: every channel adjusted
for cell area, cell count and the other four channels). glycogene column = GlycoEnzOnto (403 genes).
Genes: those with a delta in all 5 channels (the same set for every channel and permutation).
Statistic: for each channel c, p_c = share of glycogenes whose |delta| exceeds the 95th percentile of |delta| of the
non-glycogenes (np.percentile, linear interpolation, as channel_specificity.py);
T = mean(p_ER, p_AGP) - mean(p_DNA, p_RNA, p_Mito).
Permutation: the glycogene labels are shuffled over the genes (number of glycogenes fixed), each gene keeping its five
channel values together; the non-glycogene 95th percentile is recomputed each time; 10,000 permutations, seed 20261124.
One-sided p = (1 + #(T_perm >= T_obs)) / (10,000 + 1). ORF and CRISPR separately.
Outputs (<out>/): channel_enrichment_permutation.csv (one row per perturbation x model: n genes, n glycogenes, p of each
channel, T, permutation mean / 95th percentile of T, p), permutation_T_<model>.csv.gz (all permuted T),
channel_enrichment_permutation.png / .pdf (permutation distribution and observed T; main model).
Usage: python channel_enrichment_permutation.py <results/17_channel_specificity dir> <out dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from figure_style import GRID, INK, INK2, LW_REF, WIDTH, apply, despine, letter, save  # noqa: E402

CH = ['ConA', 'WGA', 'DNA', 'RNA', 'Mito']
LECTIN, OTHER = [0, 1], [2, 3, 4]
NPERM, SEED = 10_000, 20261124


def shares(a, is_g):
    q = np.percentile(a[~is_g], 95, axis=0)
    return (a[is_g] > q).mean(axis=0)


def stat(p):
    return p[LECTIN].mean() - p[OTHER].mean()


def run(d):
    d = d.dropna(subset=[f'{c}_delta' for c in CH])
    a = d[[f'{c}_delta' for c in CH]].abs().to_numpy()
    g = d.glycogene.astype(bool).to_numpy()
    p_obs = shares(a, g)
    t_obs = stat(p_obs)
    rng = np.random.default_rng(SEED)
    t_perm = np.empty(NPERM)
    for i in range(NPERM):
        t_perm[i] = stat(shares(a, rng.permutation(g)))
    pval = (1 + (t_perm >= t_obs).sum()) / (NPERM + 1)
    return len(d), int(g.sum()), p_obs, t_obs, t_perm, pval


def main(r, out):
    r, out = Path(r), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows, perms = [], {}
    for model in ('main', 'symmetric'):
        for t in ('ORF', 'CRISPR'):
            n, ng, p_obs, t_obs, t_perm, pval = run(pd.read_csv(r / f'{t}_{model}_effects.csv'))
            rows.append(dict(model=model, perturbation=t, n_genes=n, n_glycogenes=ng,
                             **{f'p_{c}': v for c, v in zip(CH, p_obs)},
                             mean_p_lectin=p_obs[LECTIN].mean(), mean_p_other=p_obs[OTHER].mean(), T_obs=t_obs,
                             T_perm_mean=t_perm.mean(), T_perm_q95=np.percentile(t_perm, 95), n_perm=NPERM,
                             n_perm_ge_obs=int((t_perm >= t_obs).sum()), p_one_sided=pval, seed=SEED))
            perms[(model, t)] = t_perm
            print(model, t, n, ng, np.round(p_obs, 4), round(t_obs, 4), pval, flush=True)
        pd.DataFrame({t: perms[(model, t)] for t in ('ORF', 'CRISPR')}).to_csv(out / f'permutation_T_{model}.csv.gz',
                                                                                index=False)
    res = pd.DataFrame(rows)
    res.to_csv(out / 'channel_enrichment_permutation.csv', index=False)
    apply()
    fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.4), gridspec_kw=dict(wspace=0.3, left=0.08, right=0.98,
                                                                          bottom=0.2, top=0.88))
    for ax, t, lab in zip(axs, ('ORF', 'CRISPR'), ('Overexpression', 'Knockout')):
        row = res[(res.model == 'main') & (res.perturbation == t)].iloc[0]
        ax.hist(perms[('main', t)] * 100, bins=60, color='#cfcdc6', edgecolor='none')
        ax.axvline(row.T_obs * 100, color='#1f66c2', lw=1.2)
        ax.axvline(0, color=INK2, lw=LW_REF, ls=':')
        ptxt = (f'p = {row.p_one_sided:.4f}' if row.n_perm_ge_obs > 0
                else f'p ≤ {1 / (NPERM + 1):.4f} (0 of {NPERM:,})')
        ax.text(0.02, 0.95, f'{lab}\nT = {row.T_obs * 100:.1f} points\n{ptxt}', transform=ax.transAxes, ha='left',
                va='top', color=INK)
        ax.set_xlabel('T = mean(ER, AGP) − mean(DNA, RNA, Mito)\n(percentage points, permuted labels)')
        ax.set_ylabel('Permutations')
        ax.grid(axis='y', color=GRID, lw=0.4)
        ax.set_axisbelow(True)
        despine(ax)
    letter(fig, 0.005, 0.99, 'A')
    letter(fig, 0.5, 0.99, 'B')
    save(fig, str(out / 'channel_enrichment_permutation'))
    pd.set_option('display.width', 250)
    print(res.round(4).to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:3])
