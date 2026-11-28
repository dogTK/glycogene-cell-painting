# GO Biological Process over-representation of all genes whose ConA (ER) or WGA (AGP) staining went up or down
# in JUMP-CP (user, 2026-11-22). Reuses the existing per-gene tests (scripts/13_jump/vs_other_genes_test.py
# --exclude-low-infection: <in_dir>/<ORF|CRISPR>_tests_all_genes.csv; Glass's delta = <lec>_diff_sd, Welch p = <lec>_p).
# Sets: perturbation (ORF, CRISPR) x lectin (ConA, WGA) x direction (up: p < 0.05 & delta > 0; down: p < 0.05 & delta < 0).
# Universe: every gene with a p value in that perturbation (control reagents removed); symbols mapped by org.Hs.eg.db.
# Test: clusterProfiler compareCluster(fun = "enrichGO", ont = "BP"), hypergeometric, BH; minGSSize 10, maxGSSize 500;
# significant: BH-adjusted p < 0.05. Redundancy: clusterProfiler::simplify (Wang semantic similarity, cutoff 0.7,
# keep the term with the smallest adjusted p) on the significant terms.
# Common vs specific: per perturbation x direction, significant terms in ConA only, WGA only, both.
# Usage: Rscript go_bp_all_genes.R <in_dir> <out_dir>
suppressPackageStartupMessages({
  library(clusterProfiler); library(enrichplot); library(org.Hs.eg.db); library(ggplot2); library(dplyr); library(readr)
})
args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]; out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
controls <- c("eGFP", "BFP", "HcRed", "LacZ", "LUCIFERASE", "Luciferase")
label <- c(ORF = "Overexpression", CRISPR = "Knockout")
FDR <- 0.05

all_res <- list(); all_sets <- list(); summary_rows <- list(); all_simp <- list()
for (t in c("ORF", "CRISPR")) {
  r <- read_csv(file.path(in_dir, paste0(t, "_tests_all_genes.csv")), show_col_types = FALSE) %>%
    filter(!gene %in% controls, !is.na(ConA_p), !is.na(WGA_p))
  universe <- r$gene
  mapped <- suppressMessages(bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  cat(sprintf("%s: universe %d genes, mapped to Entrez %d\n", t, length(universe), length(unique(mapped$SYMBOL))))
  sets <- list()
  for (lec in c("ConA", "WGA")) {
    p <- r[[paste0(lec, "_p")]]; d <- r[[paste0(lec, "_diff_sd")]]
    sets[[paste(lec, label[[t]], "up", sep = "|")]] <- r$gene[p < 0.05 & d > 0]
    sets[[paste(lec, label[[t]], "down", sep = "|")]] <- r$gene[p < 0.05 & d < 0]
  }
  cc <- compareCluster(sets, fun = "enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                       universe = universe, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
  df <- as.data.frame(cc) %>% mutate(perturbation = t)
  all_res[[t]] <- df; all_sets <- c(all_sets, sets)
  cs <- cc
  cs@compareClusterResult <- cs@compareClusterResult %>% filter(p.adjust < FDR)
  if (nrow(cs@compareClusterResult) > 0) {
    all_simp[[t]] <- as.data.frame(simplify(cs, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")) %>%
      mutate(perturbation = t)
  }
  for (nm in names(sets)) {
    n_sig <- sum(df$Cluster == nm & df$p.adjust < FDR)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(set = nm, perturbation = t, n_universe = length(universe),
                                                          n_candidates = length(sets[[nm]]), n_go_sig = n_sig)
  }
}
res <- bind_rows(all_res)
write_csv(res, file.path(out_dir, "go_bp_all.csv"))
summ <- bind_rows(summary_rows)

# Redundancy reduction of the significant terms (done per perturbation inside the loop)
sig <- res %>% filter(p.adjust < FDR)
simp <- if (length(all_simp)) bind_rows(all_simp) else data.frame(Cluster = character())
write_csv(simp, file.path(out_dir, "go_bp_significant_simplified.csv"))
summ$n_go_sig_simplified <- sapply(summ$set, function(s) sum(simp$Cluster == s))
write_csv(summ, file.path(out_dir, "go_bp_summary.csv"))
print(summ, row.names = FALSE)

# Common vs specific (ConA vs WGA) per perturbation x direction, on all significant terms
cmp <- list()
for (t in label) for (dirn in c("up", "down")) {
  a <- sig$ID[sig$Cluster == paste("ConA", t, dirn, sep = "|")]
  b <- sig$ID[sig$Cluster == paste("WGA", t, dirn, sep = "|")]
  cmp[[length(cmp) + 1]] <- data.frame(perturbation = t, direction = dirn, conA_only = length(setdiff(a, b)),
                                       WGA_only = length(setdiff(b, a)), both = length(intersect(a, b)),
                                       both_terms = paste(head(sig$Description[sig$ID %in% intersect(a, b) &
                                                                                 grepl("^ConA", sig$Cluster)], 15),
                                                          collapse = "; "))
}
cmp <- bind_rows(cmp)
write_csv(cmp, file.path(out_dir, "go_bp_conA_vs_WGA.csv"))
print(cmp %>% select(-both_terms), row.names = FALSE)
cat("\nTop simplified terms per set:\n")
top <- simp %>% group_by(Cluster) %>% slice_min(p.adjust, n = 8, with_ties = FALSE) %>% ungroup() %>%
  transmute(Cluster, Description, GeneRatio, BgRatio, p.adjust = signif(p.adjust, 2))
write.table(top, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
write_csv(top, file.path(out_dir, "go_bp_top_simplified.csv"))
