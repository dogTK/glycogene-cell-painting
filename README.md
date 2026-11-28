# glycogene-cell-painting

Analysis code for the manuscript

> Koreeda T. *Large-scale Cell Painting profiling reveals enriched lectin-channel responses to glycogene overexpression.*

The study compares the effects of overexpressing (ORF) and knocking out (CRISPR) GlycoEnzOnto glycogenes with those of other genes across the five Cell Painting channels of JUMP-CP, and checks selected knockout effects against PERISCOPE and the lectin-binding screen of Tsui et al.

## Data

All input data are public and are not included in this repository. The scripts download them from the original sources:

| Data | Source |
|---|---|
| JUMP-CP ORF / CRISPR / compound profiles | Cell Painting Gallery, `cpg0016-jump` (`source_all/workspace/profiles_assembled`, v1.0a) |
| JUMP-CP metadata | `jump-cellpainting/datasets` (commit 016e865) |
| Low-infection ORF wells, phenotype calls | `jump-cellpainting/2024_Chandrasekaran_Morphmap` |
| Feature-level test summary (JUMP_rr) | Zenodo record 14648922 |
| PERISCOPE gene-level profiles | Cell Painting Gallery, `cpg0021-periscope` |
| Lectin-binding CRISPR screen | Tsui et al. 2024, *Nat Commun*, Supplementary Dataset 2 |
| Original images | Cell Painting Gallery, following each plate's `load_data` file |
| Glycogene list and annotations | GlycoEnzOnto (403 genes) |

The Cell Painting Gallery data are released under CC0 1.0.

## Layout

```
scripts/
  common.py, figure_style.py, figure_style.R    shared helpers and the common figure style
  00_compound_overlap/fetch_a549_perts.py       Snowflake connection helper used by common.py
  01_load_to_snowflake/load_gene_sets.py        load the lectin gene sets (inputs/) into Snowflake
  12_manuscript_figures/                        figure assembly (section_figures_v6.py builds on v1-v5)
  13_jump/                                      JUMP-CP well-level features, tests against other genes, enrichment
  15_morphology/                                morphological phenotypes of glycogene perturbations
  16_all_genes/                                 all genes; overexpression vs knockout
  17_channel_specificity/                       per-channel effects and the channel permutation test
  18_go_all_genes/                              GO enrichment of all genes (clusterProfiler)
  22_profile_similarity/                        gene-level profiles and principal component analysis
  23_channel_groups/                            channel selectivity by glycogene group
  24_periscope/                                 comparison with PERISCOPE and Tsui et al.
  25_images/                                    representative images and their quantification
  26_compounds/                                 compound perturbations
  27_submission/                                Supplementary Tables S1-S12
inputs/lectin_gene_sets.csv                     hand-curated lectin gene sets
environment.yml                                 conda environment
```

Figure numbers written in script comments follow the drafting history and may differ from the final numbering in the manuscript.

## Environment

Python 3.11 (`environment.yml`: pandas 2.3.3, NumPy 2.4.6, SciPy 1.17.1, statsmodels 0.15.0, PyArrow 26.0.1, tifffile 2020.6.3, Matplotlib 3.11.2, Pillow 13.0.0) and R 4.4.3 (clusterProfiler 4.14.4, enrichplot 1.26.6, GOSemSim 2.32.0, org.Hs.eg.db 3.20.0, ggplot2 3.5.1).

Intermediate tables were kept in Snowflake. The connection is read from `~/.snowflake/config.toml` (key-pair authentication); no credentials are stored in this repository. Local data and results are written under `~/research/glyco_cell_painting/`.

The scripts that load data into Snowflake through S3 read the bucket and the external stage from environment variables:

| Variable | Meaning |
|---|---|
| `GLYCO_S3_BUCKET` | S3 bucket for the staged Parquet files and archived images (bucket name only) |
| `GLYCO_SNOWFLAKE_STAGE` | Snowflake external stage on that bucket, e.g. `@DB.SCHEMA.MY_STAGE` |

## License

MIT (see `LICENSE`).
