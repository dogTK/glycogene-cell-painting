"""Are the glycogene effects selective for the ConA (ER) and WGA (AGP) channels? (user, 2026-11-22)
Input: channel_effects.py outputs (<ORF|CRISPR>_<model>_effects.csv). ORF and CRISPR are analysed separately.

1. Per channel: ECDF of |Glass's delta|, GlycoEnzOnto genes vs other genes; two-sample Kolmogorov-Smirnov test.
2. Scale-free comparison across channels (the |delta| scale differs between channels):
   - KS D (rank based, so independent of the channel's scale), against a permutation null (KS D of 1,000 random gene
     sets of the glycogene set size in the same channel; median and 95th percentile);
   - share of glycogenes whose |delta| exceeds the 95th percentile of the other genes in the same channel (5% expected
     if glycogenes behave like other genes), exact binomial test and Clopper-Pearson 95% CI;
   - ratio of median |delta| (glycogenes / other genes), bootstrap 95% CI.
3. Strong ConA / WGA responders among glycogenes (FDR < 0.05 over glycogenes in ConA or WGA, or nominal p < 0.05 with
   |delta| >= 2 SD): their delta in all five channels. Spearman correlation of the signed delta between channels across
   glycogenes.
Figure 4 (180 mm wide): a ECDFs (rows ORF / CRISPR, columns ConA, WGA, DNA, RNA, Mito); b KS D with CI;
b KS D with the permutation null; c share above the other genes' 95th percentile with exact 95% CI; d delta of the strong responders in all channels.
Usage: python channel_specificity.py <effects_dir> <out_dir> [main|symmetric]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from matplotlib.ticker import FixedLocator  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5,
                     'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'pdf.fonttype': 42})
CH = ['ConA', 'WGA', 'DNA', 'RNA', 'Mito']
CH_LABEL = {'ConA': 'ER (ConA)', 'WGA': 'AGP (WGA)', 'DNA': 'DNA', 'RNA': 'RNA', 'Mito': 'Mito'}
PT = {'ORF': 'Overexpression', 'CRISPR': 'Knockout'}
GLY, OTH, INK, INK2, GRID = '#1f66c2', '#9a9890', '#111111', '#55534f', '#e6e5e1'
PT_COL = {'ORF': '#1f66c2', 'CRISPR': '#d9541e'}
NBOOT = 2000
NPERM = 1000


def boot_ratio(x_g, x_o, rng):
    rs = [np.median(rng.choice(x_g, len(x_g))) / np.median(rng.choice(x_o, len(x_o))) for _ in range(NBOOT)]
    return np.percentile(rs, [2.5, 97.5])


def null_D(x_all, n, rng, nperm=NPERM):
    """KS D between a random set of n genes and the rest (permutation null for the glycogene set size)."""
    out = []
    for _ in range(nperm):
        idx = rng.permutation(len(x_all))
        out.append(stats.ks_2samp(x_all[idx[:n]], x_all[idx[n:]]).statistic)
    return np.array(out)


def summarise(eff, rng):
    rows = []
    for t, r in eff.items():
        for ch in CH:
            d = r.dropna(subset=[f'{ch}_p'])
            x_g = d.loc[d.glycogene, f'{ch}_delta'].abs().to_numpy()
            x_o = d.loc[~d.glycogene, f'{ch}_delta'].abs().to_numpy()
            ks = stats.ks_2samp(x_g, x_o)
            q95 = np.percentile(x_o, 95)
            k = int((x_g > q95).sum())
            bt = stats.binomtest(k, len(x_g), 0.05, alternative='greater')
            ci = stats.binomtest(k, len(x_g), 0.05).proportion_ci(confidence_level=0.95, method='exact')
            r_ci = boot_ratio(x_g, x_o, rng)
            nd = null_D(np.r_[x_g, x_o], len(x_g), rng)
            rows.append(dict(perturbation=t, channel=ch, n_glyco=len(x_g), n_other=len(x_o), ks_D=ks.statistic,
                             ks_p=ks.pvalue, D_null_median=np.median(nd), D_null_95=np.percentile(nd, 95),
                             D_perm_p=(1 + (nd >= ks.statistic).sum()) / (1 + len(nd)), median_glyco=np.median(x_g),
                             median_other=np.median(x_o), median_ratio=np.median(x_g) / np.median(x_o),
                             ratio_lo=r_ci[0], ratio_hi=r_ci[1], q95_other=q95, n_above=k, share_above=k / len(x_g),
                             share_lo=ci.low, share_hi=ci.high, binom_p=bt.pvalue))
    return pd.DataFrame(rows)


def strong(eff):
    rows = []
    for t, r in eff.items():
        g = r[r.glycogene].copy()
        for ch in CH:
            ok = g[f'{ch}_p'].notna()
            g.loc[ok, f'{ch}_fdr'] = multipletests(g.loc[ok, f'{ch}_p'], method='fdr_bh')[1]
        hit = np.zeros(len(g), bool)
        for lec in ('ConA', 'WGA'):
            hit |= (g[f'{lec}_fdr'] < 0.05) | ((g[f'{lec}_p'] < 0.05) & (g[f'{lec}_delta'].abs() >= 2))
        s = g[hit].copy()
        s['perturbation'] = t
        s['lead'] = np.where(s.ConA_delta.abs() >= s.WGA_delta.abs(), 'ConA', 'WGA')
        s['lead_abs'] = np.maximum(s.ConA_delta.abs(), s.WGA_delta.abs())
        s['max_other_abs'] = s[[f'{c}_delta' for c in ('DNA', 'RNA', 'Mito')]].abs().max(axis=1)
        s['other_p05'] = (s[[f'{c}_p' for c in ('DNA', 'RNA', 'Mito')]] < 0.05).sum(axis=1)
        rows.append(s)
    return pd.concat(rows)


def correlations(eff):
    rows = []
    for t, r in eff.items():
        g = r[r.glycogene]
        for a in CH:
            for b in CH:
                rho = stats.spearmanr(g[f'{a}_delta'], g[f'{b}_delta'], nan_policy='omit')
                rows.append([t, a, b, rho.statistic, rho.pvalue])
    return pd.DataFrame(rows, columns=['perturbation', 'channel_a', 'channel_b', 'spearman_rho', 'p'])


def ptxt(p):
    return f'{p:.2g}' if p >= 0.001 else '{} × 10$^{{{}}}$'.format(*[int(float(v)) for v in f'{p:.0e}'.split('e')])


def figure(eff, summ, st, out):
    fig = plt.figure(figsize=(7.1, 8.9))
    ga = fig.add_gridspec(2, 5, left=0.1, right=0.98, top=0.95, bottom=0.645, wspace=0.12, hspace=0.42)
    for i, t in enumerate(('ORF', 'CRISPR')):
        r = eff[t]
        for j, ch in enumerate(CH):
            ax = fig.add_subplot(ga[i, j])
            d = r.dropna(subset=[f'{ch}_delta'])
            for m, col, lw in ((~d.glycogene, OTH, 1.1), (d.glycogene, GLY, 1.3)):
                x = np.sort(d.loc[m, f'{ch}_delta'].abs())
                ax.step(np.r_[0, x], np.r_[0, np.arange(1, len(x) + 1) / len(x)], where='post', color=col, lw=lw)
            s = summ[(summ.perturbation == t) & (summ.channel == ch)].iloc[0]
            ax.set_xscale('symlog', linthresh=1, linscale=1.2)
            ax.set_xlim(0, 30)
            ax.xaxis.set_major_locator(FixedLocator([0, 0.5, 1, 3, 10, 30]))
            ax.set_xticklabels(['0', '.5', '1', '3', '10', '30'])
            ax.set_ylim(0, 1.02)
            ax.grid(color=GRID, lw=0.4)
            ax.set_axisbelow(True)
            ax.text(0.96, 0.05, f'D = {s.ks_D:.2f}\np = {ptxt(s.ks_p)}', transform=ax.transAxes, ha='right',
                    va='bottom', fontsize=6, fontweight='bold' if s.ks_p < 0.05 else 'normal')
            if i == 0:
                ax.text(0.5, 1.06, CH_LABEL[ch], transform=ax.transAxes, ha='center', va='bottom', fontsize=7)
            if j == 0:
                ax.set_ylabel(f'{PT[t]}\nCumulative fraction')
            else:
                ax.set_yticklabels([])
            if i == 1:
                ax.set_xlabel('|Glass\'s Δ| (SD)')
            for sp in ('top', 'right'):
                ax.spines[sp].set_visible(False)
    fig.legend([plt.Line2D([], [], color=GLY, lw=1.3), plt.Line2D([], [], color=OTH, lw=1.1)],
               ['Glycogenes', 'Other genes'], loc='upper right', bbox_to_anchor=(0.98, 0.995), ncol=2, frameon=False)

    gb = fig.add_gridspec(1, 2, left=0.1, right=0.98, top=0.55, bottom=0.43, wspace=0.3)
    for k, (col, lo, hi, ylab, ref, letter) in enumerate(
            (('ks_D', None, None, 'KS D (glycogenes vs other genes)', None, 'b'),
             ('share_above', 'share_lo', 'share_hi', 'Glycogenes above the 95th\npercentile of other genes',
              0.05, 'c'))):
        ax = fig.add_subplot(gb[0, k])
        for dx, t in ((-0.14, 'ORF'), (0.14, 'CRISPR')):
            s = summ[summ.perturbation == t].set_index('channel').loc[CH]
            x = np.arange(len(CH)) + dx
            if col == 'ks_D':
                sig = s.D_perm_p < 0.05
                for xi, (m, q) in zip(x, zip(s.D_null_median, s.D_null_95)):
                    ax.plot([xi - 0.1, xi + 0.1], [q, q], color=PT_COL[t], lw=0.6, ls=(0, (2, 1.5)))
                    ax.plot([xi, xi], [m, q], color=PT_COL[t], lw=0.5, alpha=0.5)
            else:
                sig = s.binom_p < 0.05
                ax.errorbar(x, s[col], yerr=[s[col] - s[lo], s[hi] - s[col]], fmt='none', ecolor=PT_COL[t], lw=0.9)
            ax.scatter(x[~sig.values], s[col][~sig.values], s=22, facecolor='white', edgecolor=PT_COL[t], lw=1.1,
                       zorder=3)
            ax.scatter(x[sig.values], s[col][sig.values], s=22, color=PT_COL[t], zorder=3, label=PT[t])
        if ref is not None:
            ax.axhline(ref, color=INK2, lw=0.6, ls=(0, (2, 2)))
            ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
        ax.set_xticks(range(len(CH)), [CH_LABEL[c] for c in CH])
        ax.set_ylabel(ylab)
        ax.set_ylim(0, None)
        ax.grid(axis='y', color=GRID, lw=0.4)
        ax.set_axisbelow(True)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        if k == 0:
            h = [plt.Line2D([], [], marker='o', ls='', color=PT_COL['ORF'], ms=4.5),
                 plt.Line2D([], [], marker='o', ls='', color=PT_COL['CRISPR'], ms=4.5),
                 plt.Line2D([], [], marker='o', ls='', mfc='white', mec=INK2, ms=4.5)]
            h.append(plt.Line2D([], [], color=INK2, lw=0.8, ls=(0, (2, 1.5))))
            fig.legend(h, ['Overexpression, p < 0.05', 'Knockout, p < 0.05', 'p ≥ 0.05 (open)',
                           'b: 95th percentile of random gene sets'], frameon=False, loc='upper center',
                       bbox_to_anchor=(0.54, 0.405), ncol=4)

    gd = fig.add_gridspec(1, 2, left=0.14, right=0.91, top=0.33, bottom=0.085, wspace=0.45,
                          width_ratios=[len(st[st.perturbation == 'ORF']), max(len(st[st.perturbation == 'CRISPR']), 3)])
    norm = TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5)
    for k, t in enumerate(('ORF', 'CRISPR')):
        s = st[st.perturbation == t].sort_values(['lead', 'lead_abs'], ascending=[True, False])
        ax = fig.add_subplot(gd[0, k])
        m = s[[f'{c}_delta' for c in CH]].to_numpy().T
        im = ax.imshow(np.clip(m, -5, 5), cmap='RdBu_r', norm=norm, aspect='auto')
        for (yi, xi), v in np.ndenumerate(m):
            p = s[f'{CH[yi]}_p'].iloc[xi]
            ax.text(xi, yi, ('*' if p < 0.05 else ''), ha='center', va='center', fontsize=8,
                    color='white' if abs(v) > 3 else INK)
        ax.set_yticks(range(len(CH)), [CH_LABEL[c] for c in CH])
        ax.set_xticks(range(len(s)), s.index, rotation=90, fontstyle='italic', fontsize=6.3)
        ax.tick_params(length=0)
        ax.set_xlabel(PT[t], labelpad=4)
        ax.xaxis.set_label_position('top')
    cax = fig.add_axes([0.93, 0.15, 0.012, 0.13])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Glass\'s Δ (SD)', fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    for letter, (x, y) in {'a': (0.01, 0.975), 'b': (0.01, 0.58), 'c': (0.49, 0.58), 'd': (0.01, 0.35)}.items():
        fig.text(x, y, letter, fontsize=9, fontweight='bold', va='top')
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'Figure4.{ext}', dpi=300)
    plt.close(fig)


def main(src, out, model='main'):
    src, out = Path(src), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    eff = {t: pd.read_csv(src / f'{t}_{model}_effects.csv', index_col=0) for t in ('ORF', 'CRISPR')}
    rng = np.random.default_rng(20261122)
    summ = summarise(eff, rng)
    summ.to_csv(out / f'channel_summary_{model}.csv', index=False)
    st = strong(eff)
    keep = ['perturbation', 'lead', 'lead_abs', 'max_other_abs', 'other_p05'] + [f'{c}_{k}' for c in CH for k in
                                                                                ('delta', 'p')]
    st[keep].to_csv(out / f'strong_responders_{model}.csv')
    cor = correlations(eff)
    cor.to_csv(out / f'channel_correlations_{model}.csv', index=False)
    pd.set_option('display.width', 250)
    print(summ.round(4).to_string(index=False))
    print(st[keep].round(2).to_string())
    print(cor[(cor.channel_a.isin(['ConA', 'WGA'])) & (cor.channel_a != cor.channel_b)].round(3).to_string(index=False))
    if model == 'main':
        figure(eff, summ, st, out)


if __name__ == '__main__':
    main(*sys.argv[1:4])
