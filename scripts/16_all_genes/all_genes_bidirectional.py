"""All genes and ORF/KO bidirectionality for ConA and WGA staining in JUMP-CP (user, 2026-11-22).

Reuses the per-gene tests of scripts/13_jump/vs_other_genes_test.py run with --exclude-low-infection
(<in_dir>/<ORF|CRISPR>_tests_all_genes.csv): wells position-corrected and adjusted for cell area, RNA, DNA, Mito and
the cell-count proxy; each gene's wells vs the wells of other genes on the same plates (Welch's t-test); effect =
difference of means in SD of the reference wells. Control reagents (eGFP, BFP, HcRed, LacZ, Luciferase) are dropped.
Glycogene = GlycoEnzOnto (403 genes).

1. All genes: BH-FDR over all genes per perturbation x lectin. Glycogenes vs other genes: share with p < 0.05, share
   with |effect| >= 1 SD, share with FDR < 0.05 (2 x 2 tables, Fisher's exact test). Figures: effect distributions
   (ECDF of |effect|) and volcano plots with glycogenes highlighted.
2. Bidirectionality: GlycoEnzOnto genes tested both by ORF and CRISPR (glycogenes only). Scatter ORF effect (x) vs
   CRISPR effect (y) for ConA and WGA, Spearman and Pearson correlations; FDR over these glycogenes. Labels: genes with FDR < 0.05 (all-gene
   FDR) in either perturbation with |effect| >= 1.5 SD, the largest |ORF| x |KO| products with opposite sign, and glycogenes with
   |effect| >= 2 SD in both.
Usage: python all_genes_bidirectional.py <in_dir> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

CONTROLS = {'eGFP', 'BFP', 'HcRed', 'LacZ', 'LUCIFERASE', 'Luciferase'}
LECTINS = ('ConA', 'WGA')
LABEL = {'ORF': 'Overexpression', 'CRISPR': 'Knockout'}
GLY, OTH, INK, INK2, GRID, DOWN = '#2a78d6', '#b4b2ab', '#0b0b0b', '#52514e', '#e4e3df', '#eb6834'
CLIP = 8


def load(src):
    out = {}
    for t in ('ORF', 'CRISPR'):
        r = pd.read_csv(src / f'{t}_tests_all_genes.csv', index_col=0)
        r = r[~r.index.isin(CONTROLS)]
        for lec in LECTINS:
            ok = r[f'{lec}_p'].notna()
            r.loc[ok, f'{lec}_fdr_all'] = multipletests(r.loc[ok, f'{lec}_p'], method='fdr_bh')[1]
        out[t] = r
    return out


def compare(tabs, out, say):
    rows = []
    for t, r in tabs.items():
        for lec in LECTINS:
            d = r.dropna(subset=[f'{lec}_p'])
            g = d.glycogene
            for name, hit in (('p < 0.05', d[f'{lec}_p'] < 0.05), ('|effect| >= 1 SD', d[f'{lec}_diff_sd'].abs() >= 1),
                              ('FDR < 0.05 (all genes)', d[f'{lec}_fdr_all'] < 0.05)):
                a, b, c, e = (hit & g).sum(), (~hit & g).sum(), (hit & ~g).sum(), (~hit & ~g).sum()
                orr, p = fisher_exact([[a, b], [c, e]])
                rows.append([t, lec, name, g.sum(), a, a / g.sum(), (~g).sum(), c, c / (~g).sum(), orr, p])
            up = d[(d[f'{lec}_fdr_all'] < 0.05) & (d[f'{lec}_diff_sd'] > 0)]
            dn = d[(d[f'{lec}_fdr_all'] < 0.05) & (d[f'{lec}_diff_sd'] < 0)]
            say(f'{t} {lec}: genes {len(d)}; all-gene FDR < 0.05: stronger {len(up)}, weaker {len(dn)}; glycogenes among '
                f'them: {", ".join(sorted(set(up.index[up.glycogene]) | set(dn.index[dn.glycogene]))) or "none"}')
            for lab, s in (('stronger', up.sort_values(f'{lec}_diff_sd', ascending=False)),
                           ('weaker', dn.sort_values(f'{lec}_diff_sd'))):
                say(f'   {lab} (top 12 by effect): ' + ', '.join(
                    f'{i}{"*" if gl else ""} ({v:+.1f})' for i, v, gl in
                    zip(s.index[:12], s[f'{lec}_diff_sd'][:12], s.glycogene[:12])))
    c = pd.DataFrame(rows, columns=['perturbation', 'lectin', 'criterion', 'n_glycogenes', 'glyco_hits', 'glyco_share',
                                    'n_other', 'other_hits', 'other_share', 'odds_ratio', 'fisher_p'])
    c.to_csv(out / 'glyco_vs_other.csv', index=False)
    say('\nGlycogenes vs other genes:\n' + c.assign(glyco_share=(c.glyco_share * 100).round(1),
                                                    other_share=(c.other_share * 100).round(1),
                                                    odds_ratio=c.odds_ratio.round(2),
                                                    fisher_p=c.fisher_p.map(lambda x: f'{x:.2g}')).to_string(index=False))
    return c


def figure_distributions(tabs, comp, out):
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.2))
    for j, (t, lec) in enumerate([(t, lec) for t in ('ORF', 'CRISPR') for lec in LECTINS]):
        r = tabs[t].dropna(subset=[f'{lec}_p'])
        g = r.glycogene
        ax = axes[0, j]
        for m, col, lab in ((~g, OTH, f'Other genes (n={(~g).sum():,})'), (g, GLY, f'Glycogenes (n={g.sum()})')):
            x = np.sort(r.loc[m, f'{lec}_diff_sd'].abs().values)
            ax.plot(x, 1 - np.arange(len(x)) / len(x), color=col, lw=2, label=lab)
        ax.set_xscale('symlog', linthresh=1)
        ax.set_yscale('log')
        ax.set_xlim(0, 30)
        ax.set_ylim(1e-4, 1.05)
        ax.axvline(1, color=INK2, lw=0.6, ls=':')
        cc = comp[(comp.perturbation == t) & (comp.lectin == lec)].set_index('criterion')
        txt = '\n'.join(f'{k}: {cc.loc[k, "glyco_share"] * 100:.1f}% vs {cc.loc[k, "other_share"] * 100:.1f}% '
                        f'(p = {cc.loc[k, "fisher_p"]:.2g})' for k in cc.index)
        ax.text(0.02, 0.03, txt, transform=ax.transAxes, fontsize=6.5, color=INK2, va='bottom')
        ax.set_xlabel(f'{LABEL[t]} · {lec}: |effect| (SD)')
        ax.set_ylabel('Share of genes with a larger |effect|')
        ax.legend(frameon=False, fontsize=7, loc='upper right')
        ax.set_title("ABCD"[j], loc='left', fontsize=11, fontweight='bold')
        ax.grid(color=GRID, lw=0.5)
        ax = axes[1, j]
        x = r[f'{lec}_diff_sd'].clip(-CLIP, CLIP)
        y = -np.log10(r[f'{lec}_p'].clip(lower=1e-12))
        sig = r[f'{lec}_fdr_all'] < 0.05
        ax.scatter(x[~g], y[~g], s=5, color=OTH, alpha=0.5, lw=0, rasterized=True)
        ax.scatter(x[g], y[g], s=16, color=GLY, alpha=0.9, lw=0)
        ax.scatter(x[sig], y[sig], s=22, facecolor='none', edgecolor=INK, lw=0.8)
        thr = r.loc[sig, f'{lec}_p'].max() if sig.any() else None
        if thr is not None:
            ax.axhline(-np.log10(thr), color=INK2, lw=0.6, ls='--')
        lab = r[sig].reindex(r[sig][f'{lec}_diff_sd'].abs().sort_values(ascending=False).index)[:6]
        lab = pd.concat([lab, r[g & sig]]).drop_duplicates()
        for gene, row in lab.iterrows():
            ax.annotate(gene, (np.clip(row[f'{lec}_diff_sd'], -CLIP, CLIP), -np.log10(max(row[f'{lec}_p'], 1e-12))),
                        xytext=(4, 2), textcoords='offset points', fontsize=6.5, fontstyle='italic',
                        color=GLY if row.glycogene else INK, fontweight='bold' if row.glycogene else 'normal')
        ax.set_xlabel(f'{LABEL[t]} · {lec}: effect (SD; clipped at ±{CLIP})')
        ax.set_ylabel('−log10 p')
        ax.grid(color=GRID, lw=0.5)
        ax.set_title("EFGH"[j], loc='left', fontsize=11, fontweight='bold')
        ax.text(0.02, 0.97, 'ring: FDR < 0.05 over all genes\nblue: glycogenes', transform=ax.transAxes, fontsize=6.5,
                color=INK2, va='top')
    for ax in axes.flat:
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / 'all_genes_distributions.png', dpi=200, bbox_inches='tight')
    fig.savefig(out / 'all_genes_distributions.pdf', bbox_inches='tight')
    plt.close(fig)


def bidirectional(tabs, out, say):
    """GlycoEnzOnto genes tested by both ORF and CRISPR only (user, 2026-11-22: '2 は全遺伝子いりません')."""
    o, c = tabs['ORF'], tabs['CRISPR']
    both = o.index.intersection(c.index)
    both = both[o.loc[both, 'glycogene'].values]
    say(f'\n=== ORF/KO bidirectionality (GlycoEnzOnto genes in both): {len(both)}')
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.6))
    rows = []
    for ax, lec in zip(axes, LECTINS):
        d = pd.DataFrame({'orf': o.loc[both, f'{lec}_diff_sd'], 'ko': c.loc[both, f'{lec}_diff_sd'],
                          'orf_p': o.loc[both, f'{lec}_p'], 'ko_p': c.loc[both, f'{lec}_p']}).dropna(subset=['orf', 'ko'])
        for col in ('orf', 'ko'):  # FDR over the glycogenes in this set
            d[f'{col}_fdr'] = multipletests(d[f'{col}_p'], method='fdr_bh')[1]
        d['opposite_product'] = np.where(np.sign(d.orf) != np.sign(d.ko), d.orf.abs() * d.ko.abs(), 0)
        d.to_csv(out / f'bidirectional_{lec}.csv')
        rs, ps = stats.spearmanr(d.orf, d.ko)
        rp, pp = stats.pearsonr(d.orf, d.ko)
        rows.append([lec, len(d), rs, ps, rp, pp])
        say(f'{lec}: n = {len(d)}; Spearman rho = {rs:.3f} (p = {ps:.2g}); Pearson r = {rp:.3f} (p = {pp:.2g})')
        both_sig = d[(d.orf_p < 0.05) & (d.ko_p < 0.05)]
        say('   p < 0.05 in both: ' + ', '.join(f'{i} ({r.orf:+.1f}/{r.ko:+.1f})' for i, r in both_sig.iterrows()))
        top_opp = d.sort_values('opposite_product', ascending=False).head(8)
        say('   largest opposite-sign |ORF|x|KO|: ' + ', '.join(f'{i} ({r.orf:+.1f}/{r.ko:+.1f})'
                                                             for i, r in top_opp.iterrows()))
        xc, yc = d.orf.clip(-CLIP, CLIP), d.ko.clip(-CLIP, CLIP)
        quad = np.where((d.orf_p < 0.05) & (d.ko_p < 0.05) & (np.sign(d.orf) != np.sign(d.ko)), 'both',
                        np.where((d.orf_p < 0.05) | (d.ko_p < 0.05), 'one', 'none'))
        for q, col, sz, lab in (('none', OTH, 18, 'p ≥ 0.05 in both'), ('one', '#9cc0ea', 30, 'p < 0.05 in one'),
                                ('both', GLY, 60, 'p < 0.05 in both, opposite sign')):
            m = quad == q
            ax.scatter(xc[m], yc[m], s=sz, color=col, lw=0, alpha=0.9, label=f'{lab} ({m.sum()})')
        sig = (d.orf_fdr < 0.05) | (d.ko_fdr < 0.05)
        ax.scatter(xc[sig], yc[sig], s=80, facecolor='none', edgecolor=INK, lw=1.1,
                   label='FDR < 0.05 (glycogenes) in either')
        ax.axhline(0, color=INK2, lw=0.6)
        ax.axvline(0, color=INK2, lw=0.6)
        lab = set(d.index[sig]) | set(top_opp.index[:6]) | set(d.index[quad == 'both'])
        lab |= set(d.index[(d.orf.abs() >= 2.5) | (d.ko.abs() >= 2)])
        for gene in lab:
            r = d.loc[gene]
            ax.annotate(gene + ('' if abs(r.orf) <= CLIP and abs(r.ko) <= CLIP else f' ({r.orf:+.0f}/{r.ko:+.0f})'),
                        (xc[gene], yc[gene]), xytext=(5, 3), textcoords='offset points', fontsize=7.5,
                        fontstyle='italic', color=INK, fontweight='bold' if quad[d.index.get_loc(gene)] == 'both' else 'normal')
        ax.set_xlabel(f'{lec}: overexpression effect (SD; clipped at ±{CLIP})')
        ax.set_ylabel(f'{lec}: knockout effect (SD)')
        ax.set_title("AB"[LECTINS.index(lec)], loc='left', fontsize=11, fontweight='bold')
        ax.legend(frameon=False, fontsize=7.5, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)
        ax.grid(color=GRID, lw=0.5)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    pd.DataFrame(rows, columns=['lectin', 'n', 'spearman_rho', 'spearman_p', 'pearson_r', 'pearson_p']).to_csv(
        out / 'bidirectional_correlations.csv', index=False)
    fig.tight_layout()
    fig.savefig(out / 'orf_vs_ko.png', dpi=200, bbox_inches='tight')
    fig.savefig(out / 'orf_vs_ko.pdf', bbox_inches='tight')
    plt.close(fig)


def main(src, out):
    src, out = Path(src), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    log = open(out / 'summary.txt', 'w')

    def say(*a):
        print(*a)
        print(*a, file=log)

    tabs = load(src)
    for t, r in tabs.items():
        r.to_csv(out / f'{t}_all_genes_fdr.csv')
    comp = compare(tabs, out, say)
    figure_distributions(tabs, comp, out)
    bidirectional(tabs, out, say)
    log.close()


if __name__ == '__main__':
    main(*sys.argv[1:3])
