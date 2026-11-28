"""Compound perturbations in JUMP-CP compared with ORF and CRISPR on the same footing (user, 2026-11-23: "摂動化合物
11 万種類のやつも一通り回してサプリに突っ込んどいて", Discussion: compounds change morphology less than KO / ORF and their
effects are hard to separate from cytotoxicity).

Data: compound wells (fetch_compound_wells.py; 1,716 plates, 10 sources) and ORF / CRISPR wells
(data/jump/jump_genetic_well_features.parquet = RAW.GLYCO_CELL_PAINTING.JUMP_GENETIC_WELL_FEATURES); the same
CellProfiler columns: cell mean intensity of ER, AGP, RNA, DNA, Mito, cell area, Nuclei_Number_Object_Number (cell
count proxy). ORF wells flagged as low infection efficiency are dropped (as in the main analysis).
Identical processing for the three modalities:
  1. Per plate, robust z against the plate's negative-control wells (compounds: DMSO; ORF: BFP, HcRed, LUCIFERASE,
     LacZ; CRISPR: non-targeting, no-guide) (median, 1.4826 x MAD), as jump_map.normalise.
  2. Perturbation value = median over its wells (treatments only; controls dropped). Perturbations with >= 3 wells.
  3. Adjusted lectin channels: perturbation-level ER and AGP regressed on cell area, RNA, DNA, Mito and cell count
     within each modality (as the 'main' model); residuals scaled by 1.4826 x MAD.
Summary per modality:
  - share of perturbations with |z| >= 2 in each channel, cell area and cell count;
  - cytotoxicity: share with cell-count z <= -2 among perturbations with |ER| >= 2 or |AGP| >= 2, and Spearman rho of
    ER and AGP with cell count across perturbations;
  - R^2 of the adjustment regression (how much of the ER / AGP variation is explained by cell size, cell count and the
    other channels) and the share with |adjusted| >= 2.
Named compounds: known glycosylation inhibitors present in JUMP (PubChem InChIKey matched on the first 14 characters:
castanospermine, 1-deoxynojirimycin / 1-deoxymannojirimycin (same InChIKey skeleton; JUMP keys carry no stereo), eliglustat,
2-deoxy-D-glucose, NGI-1) and compounds whose annotated targets (RAW.GLYCO_CELL_PAINTING.CP_REPURPOSING_MOA) include a
GlycoEnzOnto gene.
Outputs (<out>/): modality_summary.csv, perturbation_values_<modality>.parquet, named_compounds.csv
Usage: python compounds_analysis.py <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '22_profile_similarity'))
from common import cursor  # noqa: E402
from control_modules import low_infection_mask  # noqa: E402

D = Path.home() / 'research/glyco_cell_painting/data'
R = Path.home() / 'research/glyco_cell_painting/results'
CH = ['ER', 'AGP', 'RNA', 'DNA', 'Mito']
COLS = {f'Cells_Intensity_MeanIntensity_{c}': c for c in CH}
COLS.update({'Cells_AreaShape_Area': 'area', 'Nuclei_Number_Object_Number': 'count'})
FEATS = CH + ['area', 'count']
NEG = {'ORF': {'JCP2022_915128', 'JCP2022_915129', 'JCP2022_915130', 'JCP2022_915131'},
       'CRISPR': {'JCP2022_800001', 'JCP2022_800002'}}


def plate_z(w, is_neg):
    out = []
    for _, p in w.groupby('plate'):
        ref = p[is_neg.loc[p.index]]
        if len(ref) < 8:
            ref = p
        p = p.copy()
        for f in FEATS:
            med = np.nanmedian(ref[f]); mad = 1.4826 * np.nanmedian(np.abs(ref[f] - med))
            p[f] = (p[f] - med) / mad if mad > 0 else np.nan
        out.append(p)
    return pd.concat(out)


def perturbation_values(w, is_neg):
    z = plate_z(w, is_neg)
    z = z[~is_neg.loc[z.index]]
    g = z.groupby('pert')[FEATS].median()
    g['n_wells'] = z.groupby('pert').size()
    g = g[g.n_wells >= 3].dropna(subset=FEATS)
    X = np.column_stack([np.ones(len(g)), g[['area', 'RNA', 'DNA', 'Mito', 'count']].to_numpy()])
    r2 = {}
    for c in ('ER', 'AGP'):
        beta, *_ = np.linalg.lstsq(X, g[c].to_numpy(), rcond=None)
        res = g[c].to_numpy() - X @ beta
        r2[c] = 1 - res.var() / g[c].var()
        g[f'{c}_adj'] = (res - np.median(res)) / (1.4826 * np.median(np.abs(res - np.median(res))))
    return g, r2


def summarise(g, r2, mod):
    row = dict(modality=mod, n_perturbations=len(g), median_wells=g.n_wells.median())
    for f in FEATS:
        row[f'share_abs_z_ge2_{f}'] = (g[f].abs() >= 2).mean()
    lec = (g.ER.abs() >= 2) | (g.AGP.abs() >= 2)
    row['n_lectin_ge2'] = int(lec.sum())
    row['share_toxic_among_lectin_ge2'] = (g.loc[lec, 'count'] <= -2).mean()
    row['share_toxic_all'] = (g['count'] <= -2).mean()
    for c in ('ER', 'AGP'):
        row[f'rho_{c}_count'] = stats.spearmanr(g[c], g['count']).statistic
        row[f'r2_{c}_adjustment'] = r2[c]
        row[f'share_abs_adj_ge2_{c}'] = (g[f'{c}_adj'].abs() >= 2).mean()
    return row


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows, vals = [], {}
    # compounds
    c = pd.read_parquet(D / 'jump_compound/compound_wells_raw.parquet').rename(columns=COLS)
    c = c.rename(columns={'Metadata_Plate': 'plate', 'Metadata_JCP2022': 'pert'})
    c['pert'] = c.pert.astype(str)
    neg = c.Metadata_pert_type.astype(str).eq('negcon')
    keep = neg | c.Metadata_pert_type.astype(str).eq('trt')
    c, neg = c[keep], neg[keep]
    g, r2 = perturbation_values(c, neg)
    vals['compound'] = g
    rows.append(summarise(g, r2, 'compound'))
    print('compound done', len(g), flush=True)
    # genetic
    wf = pd.read_parquet(D / 'jump/jump_genetic_well_features.parquet').rename(columns=COLS)
    well = pd.read_csv(D / 'jump/meta/well.csv.gz')
    plate = pd.read_csv(D / 'jump/meta/plate.csv.gz')
    wf = wf.merge(well, on=['Metadata_Source', 'Metadata_Plate', 'Metadata_Well']).merge(
        plate[['Metadata_Plate', 'Metadata_PlateType']], on='Metadata_Plate')
    ctrl = set(pd.read_csv(D / 'jump/meta/perturbation_control.csv').Metadata_JCP2022)
    for mod in ('ORF', 'CRISPR'):
        w = wf[wf.Metadata_PlateType == mod].copy()
        if mod == 'ORF':
            w = w[~low_infection_mask(w)]
        meta = pd.read_csv(D / f'jump_meta/{mod.lower()}.csv.gz')
        sym = dict(zip(meta.Metadata_JCP2022, meta.Metadata_Symbol))
        is_neg = w.Metadata_JCP2022.isin(NEG[mod])
        w = w[is_neg | (~w.Metadata_JCP2022.isin(ctrl) & w.Metadata_JCP2022.map(sym).notna())]
        is_neg = w.Metadata_JCP2022.isin(NEG[mod])
        w = w.rename(columns={'Metadata_Plate': 'plate'})
        w['pert'] = np.where(is_neg, w.Metadata_JCP2022, w.Metadata_JCP2022.map(sym))
        g, r2 = perturbation_values(w, is_neg)
        vals[mod] = g
        rows.append(summarise(g, r2, mod))
        print(mod, 'done', len(g), flush=True)
    s = pd.DataFrame(rows)
    s.to_csv(out / 'modality_summary.csv', index=False)
    for k, v in vals.items():
        v.to_parquet(out / f'perturbation_values_{k}.parquet')
    # named compounds
    cmp = pd.read_csv(D / 'jump_meta/compound.csv.gz')
    cmp['ik14'] = cmp.Metadata_InChIKey.str[:14]
    inh = pd.read_csv(D / 'jump_compound/glycosylation_inhibitors_pubchem.csv')
    inh['ik14'] = inh.inchikey.str[:14]
    glyco = set(pd.read_csv(R / '13_jump/nominal_enrichment_infection_filtered/glycoenzonto_term2gene.tsv',
                            sep='\t').gene)
    cur = cursor()
    cur.execute('select distinct INCHIKEY14, PERT_INAME, MOA, TARGET from RAW.GLYCO_CELL_PAINTING.CP_REPURPOSING_MOA')
    moa = pd.DataFrame(cur.fetchall(), columns=['ik14', 'name', 'moa', 'target'])
    moa = moa[moa.target.fillna('').str.split('|').apply(lambda t: bool(set(t) & glyco))]
    named = pd.concat([inh[['name', 'ik14']].assign(group='glycosylation inhibitor', moa='', target=''),
                       moa.assign(group='annotated glycogene target')])
    named = named.merge(cmp[['ik14', 'Metadata_JCP2022']], on='ik14').drop_duplicates(['Metadata_JCP2022', 'name'])
    gv = vals['compound']
    named = named.join(gv, on='Metadata_JCP2022', how='inner')
    named.to_csv(out / 'named_compounds.csv', index=False)
    pd.set_option('display.width', 250)
    print(s.round(3).T.to_string())
    print(named[['group', 'name', 'target', 'n_wells', 'ER', 'AGP', 'DNA', 'RNA', 'Mito', 'count', 'ER_adj',
                 'AGP_adj']].round(2).to_string(index=False))


if __name__ == '__main__':
    main(sys.argv[1])
