"""How does perturbing a GlycoEnzOnto gene change cell morphology? (user, 2026-11-22: 'まず見るところから', then
'全部解析終わらしといてね'). Uses the JUMP-CP published tables only; nothing is recomputed from wells.

Inputs (downloaded, untrusted data, read only):
- Phenotypic activity: jump-cellpainting/2024_Chandrasekaran_Morphmap 03.retrieve-annotations/output
  phenotypic-activity-<ops>.csv.gz (copairs mAP of replicate retrieval vs negative controls on the same plate,
  permutation p, BH; 'below_corrected_p' = active at the paper's FDR 0.1). ORF ops = ..._sphering_harmony,
  CRISPR ops = ..._sphering_harmony_PCA_corrected. Reproduces the paper (ORF 57% of reagents, CRISPR 70%).
- Interpretable features: Zenodo 14648922 <orf|crispr>_interpretable_features.parquet (JUMP_rr). Per reagent and
  CellProfiler feature (profiles before batch correction, original feature meaning): 'Feature significance' =
  BH-corrected p of a t-test of the reagent's wells vs negative-control wells on the same plates
  (monorepo/libs/jump_rr significance.py), 'Median' = the reagent's median value. Only the top features per reagent
  and top reagents per feature are stored, so per reagent we use its 30 most significant features.
Gene level: ORF has several reagents per gene; a gene is active if any of its reagents is active.

Q1  Share of GlycoEnzOnto genes with a phenotype (and, as a reference, of the other genes; Fisher's exact test).
Q2  Per GlycoEnzOnto pathway (GMT; >= 8 genes tested): share active vs the other GlycoEnzOnto genes (Fisher, BH).
Q3  For active GlycoEnzOnto genes: which channel / feature type / compartment their top-30 features fall in,
    and which pathways go with which channel.
Usage: python glyco_morphology.py <data_dir> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
plt.rcParams.update({'font.family': 'Arial', 'pdf.fonttype': 42})
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cursor  # noqa: E402

ACT = {'ORF': 'phenotypic-activity-wellpos_cc_var_mad_outlier_featselect_sphering_harmony.csv.gz',
       'CRISPR': 'phenotypic-activity-wellpos_cc_var_mad_outlier_featselect_sphering_harmony_PCA_corrected.csv.gz'}
TOPN = 30
MIN_PATHWAY = 8
CHANNELS = ['DNA', 'ER', 'RNA', 'AGP', 'Mito']
UP, DOWN, GRAY, INK, INK2, GRID = '#2a78d6', '#eb6834', '#b4b2ab', '#0b0b0b', '#52514e', '#e4e3df'


def feature_type(f):
    for k in ('AreaShape', 'Intensity', 'Texture', 'Granularity', 'RadialDistribution', 'Correlation',
              'Neighbors', 'Location', 'Number', 'Children', 'Parent', 'ObjectSkeleton'):
        if f.startswith(k):
            return k
    return 'Other'


def main(data, out):
    data, out = Path(data), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with cursor() as cur:
        cur.execute('select distinct "geneName" from RAW.GLYCODB.GLYCOENZONTO_GLYCOGENE')
        geo = {r[0] for r in cur.fetchall()}
        cur.execute('select PATHWAY_NAME, GENE_SYMBOLS from RAW.GLYCODB.GLYCOENZONTO_PATHWAY_GMT')
        gmt = {n: set(s.split('|')) for n, s in cur.fetchall() if n != 'glycosylation-related pathway'}
        sym = {}
        for t in ('ORF', 'CRISPR'):
            cur.execute(f'select METADATA_JCP2022, METADATA_SYMBOL from RAW.GLYCO_CELL_PAINTING.JUMP_{t}')
            sym[t] = dict(cur.fetchall())
    log = open(out / 'summary.txt', 'w')

    def say(*a):
        print(*a)
        print(*a, file=log)

    say(f'GlycoEnzOnto genes: {len(geo)}; pathways used for Q2: GMT with >= {MIN_PATHWAY} tested genes')
    gene_tables, path_rows, feat_rows = {}, [], []
    for t in ('ORF', 'CRISPR'):
        a = pd.read_csv(data / 'morphmap_phenotypic_activity' / ACT[t])
        a['gene'] = a.Metadata_JCP2022.map(sym[t])
        a = a.dropna(subset=['gene'])
        g = a.groupby('gene').agg(active=('below_corrected_p', 'any'), best_map=('mean_average_precision', 'max'),
                                  best_q=('corrected_p_value', 'min'), n_reagents=('Metadata_JCP2022', 'size'))
        g['glycogene'] = g.index.isin(geo)
        gene_tables[t] = g
        # Q1
        ga, gn = g[g.glycogene].active.sum(), (~g[g.glycogene].active).sum()
        oa, on = g[~g.glycogene].active.sum(), (~g[~g.glycogene].active).sum()
        orr, p = fisher_exact([[ga, gn], [oa, on]])
        say(f'\n=== {t}: reagents {len(a)}, genes {len(g)} (active {g.active.sum()}, {g.active.mean():.1%})')
        say(f'Q1 GlycoEnzOnto genes tested {ga + gn}/{len(geo)}; active {ga} ({ga / (ga + gn):.1%}); '
            f'other genes {oa / (oa + on):.1%}; OR {orr:.2f}, Fisher p {p:.2g}')
        # Q2
        gg = g[g.glycogene]
        rows = []
        for name, genes in gmt.items():
            inp = gg.index.isin(genes)
            if inp.sum() < MIN_PATHWAY:
                continue
            a1, b1 = gg[inp].active.sum(), (~gg[inp].active).sum()
            c1, d1 = gg[~inp].active.sum(), (~gg[~inp].active).sum()
            o_, p_ = fisher_exact([[a1, b1], [c1, d1]])
            rows.append([t, name, int(inp.sum()), int(a1), a1 / inp.sum(), c1 / (c1 + d1), o_, p_,
                         ', '.join(sorted(gg[inp & gg.active.values].index))])
        r = pd.DataFrame(rows, columns=['perturbation', 'pathway', 'n_tested', 'n_active', 'share_active',
                                        'share_active_other_glycogenes', 'odds_ratio', 'p', 'active_genes'])
        r['fdr'] = multipletests(r.p, method='fdr_bh')[1]
        r = r.sort_values('p')
        path_rows.append(r)
        say(f'Q2 pathways tested {len(r)}; FDR < 0.05: {(r.fdr < 0.05).sum()}; FDR < 0.1: {(r.fdr < 0.1).sum()}')
        for _, x in r.head(8).iterrows():
            say(f'   {x.pathway}: {x.n_active}/{x.n_tested} ({x.share_active:.0%}) vs {x.share_active_other_glycogenes:.0%}'
                f' OR {x.odds_ratio:.2f} p {x.p:.2g} FDR {x.fdr:.2f}')
        # Q3
        f = pd.read_parquet(data / 'jump_rr_14648922' / f'{t.lower()}_interpretable_features.parquet',
                            columns=['JCP2022 ID', 'Compartment', 'Feature', 'Channel', 'Suffix',
                                     'Feature significance', 'Median'])
        f['gene'] = f['JCP2022 ID'].map(sym[t])
        act_reag = set(a.loc[a.below_corrected_p, 'Metadata_JCP2022'])
        f = f[f['JCP2022 ID'].isin(act_reag) & f.gene.notna()]
        f = (f.sort_values('Feature significance').groupby('JCP2022 ID', sort=False).head(TOPN))
        f['type'] = f.Feature.map(feature_type)
        f['channel'] = f.Channel.replace('', 'none').fillna('none')
        f['glycogene'] = f.gene.isin(geo)
        f['perturbation'] = t
        f['direction'] = np.where(f.Median > 0, 'up', 'down')
        feat_rows.append(f[['perturbation', 'gene', 'JCP2022 ID', 'glycogene', 'Compartment', 'Feature', 'type',
                            'channel', 'Suffix', 'Feature significance', 'Median', 'direction']])
        for col in ('channel', 'type', 'Compartment'):
            tab = pd.crosstab(f[col], f.glycogene, normalize='columns').rename(
                columns={True: 'glycogenes', False: 'other genes'})
            say(f'Q3 {t} top-{TOPN} features of active genes by {col} (share):\n'
                + (tab * 100).round(1).sort_values('glycogenes', ascending=False).to_string())
    pd.concat({k: v for k, v in gene_tables.items()}, names=['perturbation']).to_csv(out / 'gene_activity.csv')
    paths = pd.concat(path_rows)
    paths.to_csv(out / 'pathway_activity.csv', index=False)
    feats = pd.concat(feat_rows)
    feats.to_csv(out / 'top_features_active_genes.csv.gz', index=False)

    # Q3b: pathway x channel for active GlycoEnzOnto genes (share of their top features in each channel)
    pc_rows = []
    for t in ('ORF', 'CRISPR'):
        ft = feats[(feats.perturbation == t) & feats.glycogene]
        base = ft.channel.value_counts(normalize=True)
        for name, genes in gmt.items():
            sub = ft[ft.gene.isin(genes)]
            if sub.gene.nunique() < 4:
                continue
            sh = sub.channel.value_counts(normalize=True)
            for ch in CHANNELS + ['none']:
                pc_rows.append([t, name, sub.gene.nunique(), ch, sh.get(ch, 0), base.get(ch, 0)])
    pc = pd.DataFrame(pc_rows, columns=['perturbation', 'pathway', 'n_active_genes', 'channel', 'share',
                                        'share_all_active_glycogenes'])
    pc['diff'] = pc.share - pc.share_all_active_glycogenes
    pc.to_csv(out / 'pathway_channel.csv', index=False)
    say('\nQ3b pathways whose active genes lean most to one channel (share minus all active glycogenes):')
    for t in ('ORF', 'CRISPR'):
        top = pc[pc.perturbation == t].sort_values('diff', ascending=False).head(10)
        for _, x in top.iterrows():
            say(f'   {t} {x.pathway} (n={x.n_active_genes}): {x.channel} {x.share:.0%} vs {x.share_all_active_glycogenes:.0%}')

    figures(gene_tables, paths, feats, pc, out)
    log.close()


def figures(gene_tables, paths, feats, pc, out):
    # Figure A: share active, GlycoEnzOnto vs other genes
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw={'width_ratios': [0.5, 1.3, 1.1]})
    ax = axes[0]
    for i, t in enumerate(('ORF', 'CRISPR')):
        g = gene_tables[t]
        for j, (m, c, lab) in enumerate(((g.glycogene, UP, 'GlycoEnzOnto genes'), (~g.glycogene, GRAY, 'Other genes'))):
            v = g[m].active.mean() * 100
            ax.bar(i + (j - 0.5) * 0.36, v, width=0.34, color=c, label=lab if i == 0 else None)
            ax.text(i + (j - 0.5) * 0.36, v + 1.5, f'{v:.0f}%', ha='center', fontsize=8, color=INK)
            ax.text(i + (j - 0.5) * 0.36, 3, f'n={m.sum():,}', ha='center', fontsize=6.5, color='white', rotation=90,
                    va='bottom')
    ax.set_xticks([0, 1], ['Over-\nexpression', 'Knockout'])
    ax.set_ylim(0, 100)
    ax.set_ylabel('Genes with a morphological phenotype (%)')
    ax.legend(frameon=False, fontsize=7.5, loc='upper left')
    ax.set_title('A', loc='left', fontsize=11, fontweight='bold')
    # Figure B: pathways (dot: share active, size n, ring if FDR < 0.1)
    ax = axes[1]
    narrow = paths[paths.n_tested <= 60]
    sel = narrow[narrow.pathway.isin(narrow.sort_values('p').groupby('perturbation').head(10).pathway)]
    order = sel.groupby('pathway').p.min().sort_values(ascending=False).index
    y = {p: i for i, p in enumerate(order)}
    for t, c, dx in (('ORF', UP, -0.15), ('CRISPR', DOWN, 0.15)):
        s = sel[sel.perturbation == t]
        ax.scatter(s.share_active * 100, s.pathway.map(y) + dx, s=12 + s.n_tested * 3, color=c, alpha=0.85, lw=0,
                   label='Overexpression' if t == 'ORF' else 'Knockout')
    for t, c in (('ORF', UP), ('CRISPR', DOWN)):
        g = gene_tables[t]
        ax.axvline(g[g.glycogene].active.mean() * 100, color=c, lw=0.8, ls='--')
    ax.set_yticks(range(len(order)), [p[:55] for p in order], fontsize=7)
    ax.set_xlim(0, 105)
    ax.set_xlabel('Genes with a phenotype (%)   dashed: all GlycoEnzOnto genes')
    ax.legend(frameon=False, fontsize=7.5, loc='upper left', bbox_to_anchor=(0, -0.09), ncol=2, markerscale=0.6)
    ax.grid(axis='x', color=GRID, lw=0.5)
    ax.set_title('B', loc='left', fontsize=11, fontweight='bold')
    # Figure C: channel of the top features, glycogenes vs others
    ax = axes[2]
    chans = CHANNELS + ['none']
    w = 0.2
    for k, (t, gl, c, lab) in enumerate((('ORF', True, UP, 'Overexpression · GlycoEnzOnto'),
                                         ('ORF', False, '#9cc0ea', 'Overexpression · other'),
                                         ('CRISPR', True, DOWN, 'Knockout · GlycoEnzOnto'),
                                         ('CRISPR', False, '#f3b59b', 'Knockout · other'))):
        s = feats[(feats.perturbation == t) & (feats.glycogene == gl)].channel.value_counts(normalize=True)
        ax.bar(np.arange(len(chans)) + (k - 1.5) * w, [s.get(ch, 0) * 100 for ch in chans], width=w * 0.95,
               color=c, label=lab)
    ax.set_xticks(range(len(chans)), ['DNA', 'ER\n(ConA)', 'RNA', 'AGP\n(WGA+actin)', 'Mito', 'Shape /\nno channel'])
    ax.set_ylabel(f'Share of each gene\'s top {TOPN} features (%)')
    ax.legend(frameon=False, fontsize=7, loc='upper right')
    ax.grid(axis='y', color=GRID, lw=0.5)
    ax.set_title('C', loc='left', fontsize=11, fontweight='bold')
    for a in axes:
        for sp in ('top', 'right'):
            a.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / 'glyco_morphology.png', dpi=200, bbox_inches='tight')
    fig.savefig(out / 'glyco_morphology.pdf', bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main(*sys.argv[1:3])
