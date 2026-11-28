"""Build Supplementary Tables S1-S12 as one .xlsx file each (at most 2 MB per
file) from existing result files (no recomputation). Numbers are rounded to 4 significant digits.
Requires openpyxl (miniconda base python). Usage: python build_supplementary_tables.py <out dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path.home() / 'research/glyco_cell_painting/results'
REPO = Path(__file__).resolve().parents[2]
T = REPO / 'docs/tables'


def rd(df):
    out = df.copy()
    for c in out.select_dtypes('number').columns:
        out[c] = out[c].map(lambda v: float(f'{v:.4g}') if pd.notna(v) and np.isfinite(v) else v)
    return out


TABLES = {
    'S1': {'ORF': R / '13_jump/vs_other_genes_infection_filtered/ORF_glycogene_calls.csv',
           'CRISPR': R / '13_jump/vs_other_genes_infection_filtered/CRISPR_glycogene_calls.csv'},
    'S2': {'ORF_5channels': R / '17_channel_specificity/ORF_main_effects.csv',
           'CRISPR_5channels': R / '17_channel_specificity/CRISPR_main_effects.csv',
           'permutation_test': T / 'channel_enrichment_permutation.csv'},
    'S3': {'main_model': R / '17_channel_specificity/channel_summary_main.csv',
           'sensitivity_model': R / '17_channel_specificity/channel_summary_symmetric.csv'},
    'S4': {'ORF_all_genes': R / '16_all_genes/ORF_all_genes_fdr.csv',
           'CRISPR_all_genes': R / '16_all_genes/CRISPR_all_genes_fdr.csv'},
    'S5': {'gene_selection': T / 'image_panels_gene_selection.csv',
           'image_selection': T / 'image_panels_image_selection.csv',
           'selected_wells': T / 'image_panels_selection_wells.csv',
           'image_level_summary': T / 'image_panels_image_level_summary.csv',
           'plate_delta': T / 'image_panels_plate_delta.csv'},
    'S6': {'strong_responders': R / '17_channel_specificity/strong_responders_main.csv'},
    'S7': {'gene_phenotype': R / '15_morphology/gene_activity.csv',
           'pathway_phenotype': R / '15_morphology/pathway_activity.csv',
           'pathway_channel': R / '15_morphology/pathway_channel.csv'},
    'S8': {'pathway_enrichment': R / '13_jump/nominal_enrichment_infection_filtered/clusterprofiler_enrichment.csv',
           'prespecified_sets': R / '23_channel_groups/prespecified_sets_by_channel.csv'},
    'S9': {'go_bp_all': R / '18_go_all_genes/go_bp_all.csv',
           'go_bp_significant_simplified': R / '18_go_all_genes/go_bp_significant_simplified.csv'},
    'S10': {'ORF_scores': R / '22_profile_similarity/pca_glyco/pca_glyco_ORF_scores.csv',
            'ORF_loadings': R / '22_profile_similarity/pca_glyco/pca_glyco_ORF_loadings.csv',
            'CRISPR_scores': R / '22_profile_similarity/pca_glyco/pca_glyco_CRISPR_scores.csv',
            'CRISPR_loadings': R / '22_profile_similarity/pca_glyco/pca_glyco_CRISPR_loadings.csv',
            'variance': R / '22_profile_similarity/pca_glyco/pca_glyco_variance.csv',
            'localization_enrichment': R / '22_profile_similarity/pca_glyco/pca_glyco_compartment_enrichment.csv',
            'pathway_enrichment': R / '22_profile_similarity/pca_glyco/pca_glyco_enrichment.csv',
            'golgi_glycosyltransferase': R / '22_profile_similarity/compartment/golgi_vs_glycosyltransferase.csv'},
    'S11': {'tsui_vs_jump': T / 'external_tsui_vs_jump_ko.csv',
            'periscope_selectivity': T / 'external_periscope_channel_selectivity.csv',
            'periscope_nine_genes': T / 'external_periscope_named_genes.csv',
            'periscope_vs_jump': T / 'external_periscope_vs_jump.csv'},
    'S12': {'summary_by_modality': T / 'compounds_modality_summary.csv',
            'named_compounds': T / 'compounds_named.csv',
            'glycosylation_inhibitors': T / 'glycosylation_inhibitors_pubchem.csv'},
}


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for name, sheets in TABLES.items():
        f = out / f'Supplementary_Table_{name}.xlsx'
        with pd.ExcelWriter(f, engine='openpyxl') as w:
            for sh, p in sheets.items():
                df = pd.read_csv(p)
                if name == 'S9' and sh == 'go_bp_all':  # keep the file under 2 MB: terms with unadjusted P < 0.05
                    df = df[df.pvalue < 0.05]
                rd(df).to_excel(w, sheet_name=sh[:31], index=False)
        print(name, f.stat().st_size // 1024, 'KB', '(over 2 MB!)' if f.stat().st_size > 2_000_000 else '')


if __name__ == '__main__':
    main(sys.argv[1])
