# PCA of the GlycoEnzOnto genes' Cell Painting profiles, laid out as Figure 3 of the first paper (Koreeda & Honda
# 2024, Glycoconj J): a score plots, b loading plots, c heatmap of loadings, d glycopathway enrichment from PC scores
# (user, 2026-11-23: glycogenes only, no all-gene overlay; points uncoloured, Figure 1 genes named: "とりま 1 で").
# Input: glyco_profiles_<ORF|CRISPR>.csv from glyco_profiles.py (gene-level medians of the named JUMP-CP features).
# PCA: prcomp(scale. = TRUE) on the glycogenes of each perturbation (OE and KO separately), as in the first paper.
# Here the glycogenes are the observations (perturbed cells) and the CellProfiler features are the variables, so the
# loadings are features and the glycopathway enrichment uses the glycogenes' PC scores instead of loadings.
# a, b: square panels (user, 2026-11-23: "2 で"; axes are not on equal scales).
# a: PC1-PC2 and PC1-PC3 scores (outliers on the border, see clip_axis); the strong ConA/WGA responders of Figure 1 named.
# Loadings (b, c) = eigenvector x SD of the PC (rotation x sdev), i.e. the correlation of each feature with the PC, as in
# the first paper (user, 2026-11-23: "もちのろん").
# b: loadings in the PC1-PC2 and PC1-PC3 planes, the 10 features with the longest arrows in each plane.
# c: loadings of PC1-PC3 for the 50 features with the largest SD of loadings over PC1-PC3 (the first paper used 100; user, 2026-11-23: "50 にしよ"),
#    rows ordered by hierarchical clustering (Euclidean, complete); channel of each feature in the left strip.
# d: subcellular compartment (GlycoEnzOnto "Compartment", export_compartment_term2gene.py) of the glycogenes in the top
#    10% (+) and bottom 10% (-) of the scores of PC1, PC2 and PC3; compartments with >= 5 genes in the universe, every
#    tested compartment drawn as a row; size = genes hit, fill = -log10 p, black ring = BH-adjusted p < 0.05
#    (user, 2026-11-23: "d の図と差し替えるか"). Same test as compartment_enrichment.R.
# Supplementary (pca_glyco_pathway_enrichment_supp): the former d, GlycoEnzOnto pathways (3 genes or more, root
#    dropped) for the same sets, OE and KO side by side; terms with p < 0.05; black ring = adjusted p < 0.05; no gene
#    names. For all enrichments: clusterProfiler compareCluster(enricher), universe = the glycogenes in the PCA, BH.
#    The sign of a PC is arbitrary (as in prcomp).
# Usage: Rscript pca_glyco.R <dir with glyco_profiles_*.csv> <pathway term2gene.tsv> <figure out dir>
#        <compartment_term2gene.tsv>
suppressPackageStartupMessages({library(clusterProfiler); library(ggplot2); library(ggrepel); library(patchwork)
  library(dplyr); library(readr); library(tidyr)})
args <- commandArgs(trailingOnly = TRUE)
d <- args[1]; t2g_file <- args[2]; fig_dir <- args[3]
t2g_all <- read_tsv(t2g_file, show_col_types = FALSE)
comp_all <- read_tsv(args[4], show_col_types = FALSE)
CH <- c("ER", "AGP", "RNA", "DNA", "Mito")
CHCOL <- c(ER = "#1f66c2", AGP = "#d9541e", RNA = "#8e6bbf", DNA = "#3a3a38", Mito = "#2e9e6a",
           `2+ channels` = "#a9a7a0", `Shape / count` = "#c9c6bd")
NAMED <- list(ORF = c("A4GNT", "B3GNT2", "PIGC", "ST6GAL1", "FUT9", "B4GALNT2", "GCK", "B3GNT6"),
              CRISPR = c("MGAT1", "SLC35A2", "B4GALT1", "B3GNT2"))
# shorter feature names for labels: compartment abbreviated, '_256' (grey levels) dropped, '_' -> ' '
shortname <- function(f) gsub("_", " ", sub("_256$", "", sub("^Cytoplasm_", "Cyto_", sub("^Nuclei_", "Nuc_", f))))
channel <- function(f) sapply(strsplit(f, "_"), function(tok) {
  k <- intersect(CH, tok); if (length(k) == 1) k else if (length(k) > 1) "2+ channels" else "Shape / count"})
# common style of the manuscript figures (scripts/figure_style.R; user, 2026-11-24): Arial 7 pt, axis lines 0.6 pt,
# no panel box, figure 7.1 x 9 in
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "..", "figure_style.R"))
th <- theme_paper()

# Outliers (|score - median| > 10 x MAD on an axis; only HSPA5 on KO PC1) are drawn on the panel border as triangles
# so that they do not squash the other points (user, 2026-11-23: "外れ値は線上に乗せる").
clip_axis <- function(v, k = 10) {
  med <- median(v); mad <- mad(v)
  inl <- v[abs(v - med) <= k * mad]
  lo <- min(inl); hi <- max(inl); pad <- 0.06 * (hi - lo)
  lim <- c(if (any(v < lo)) lo - pad else lo - 0.05 * (hi - lo), if (any(v > hi)) hi + pad else hi + 0.05 * (hi - lo))
  list(v = pmin(pmax(v, lim[1]), lim[2]), out = v < lim[1] | v > lim[2], lim = lim,
       clipped = any(v < lo) || any(v > hi))
}

score_plot <- function(sc, var, pc, t) {
  cx <- clip_axis(sc$PC1); cy <- clip_axis(sc[[pc]])
  sc$x <- cx$v; sc$y <- cy$v; sc$out <- cx$out | cy$out
  lab <- sc %>% filter(gene %in% NAMED[[t]])
  ex <- ey <- expansion(0)
  ggplot(sc, aes(x, y)) +
    geom_hline(yintercept = 0, colour = "#8a8984", linewidth = 0.6 * PT) +
    geom_vline(xintercept = 0, colour = "#8a8984", linewidth = 0.6 * PT) +
    geom_point(data = sc %>% filter(!out), size = 0.9, colour = "#3a3a38", alpha = 0.7) +
    geom_point(data = sc %>% filter(out), size = 1.6, shape = 17, colour = "#3a3a38") +
    geom_point(data = lab, size = 1.6, shape = 21, fill = NA, colour = "#111111", stroke = 0.5) +
    geom_text_repel(data = lab, aes(label = gene), size = 2.1, fontface = "italic", min.segment.length = 0, segment.size = 0.4 * PT,
                    box.padding = 0.35, max.overlaps = Inf, seed = 1) +
    scale_x_continuous(limits = cx$lim, expand = ex) + scale_y_continuous(limits = cy$lim, expand = ey) +
    coord_cartesian(clip = "off") +
    labs(x = sprintf("PC1 (%.1f%%)", var[1] * 100), y = sprintf("%s (%.1f%%)", pc, var[as.integer(sub("PC", "", pc))] * 100)) +
    th + theme(aspect.ratio = 1)
}

loading_plot <- function(ld, var, pc) {
  ld$len <- sqrt(ld$PC1^2 + ld[[pc]]^2)
  top <- ld %>% slice_max(len, n = 10) %>% mutate(short = shortname(feature))
  rx <- range(c(0, top$PC1)); ry <- range(c(0, top[[pc]]))
  ggplot(top, aes(0, 0, xend = PC1, yend = .data[[pc]], colour = channel)) +
    geom_hline(yintercept = 0, colour = "#d8d6cf", linewidth = 0.6 * PT) +
    geom_vline(xintercept = 0, colour = "#d8d6cf", linewidth = 0.6 * PT) +
    geom_segment(arrow = arrow(length = unit(1.1, "mm")), linewidth = 0.7 * PT) +
    geom_text_repel(aes(x = PC1, y = .data[[pc]], label = short), size = 1.6, segment.size = 0.4 * PT, colour = "#3a3a38",
                    min.segment.length = 0, max.overlaps = Inf, seed = 1, show.legend = FALSE, force = 4,
                    force_pull = 0.3, box.padding = 0.3, max.time = 5, max.iter = 1e5) +
    scale_x_continuous(expand = expansion(mult = 0.3)) + scale_y_continuous(expand = expansion(mult = 0.3)) +
    scale_colour_manual(values = CHCOL, drop = TRUE, name = NULL) +
    labs(x = sprintf("PC1 loading (%.1f%%)", var[1] * 100),
         y = sprintf("%s loading (%.1f%%)", pc, var[as.integer(sub("PC", "", pc))] * 100)) +
    th + theme(legend.position = "none", aspect.ratio = 1)
}

heat <- function(ld) {
  ld$sd <- apply(ld[, c("PC1", "PC2", "PC3")], 1, sd)
  top <- ld %>% slice_max(sd, n = 50)
  ord <- hclust(dist(top[, c("PC1", "PC2", "PC3")]))$order
  top$feature <- factor(top$feature, levels = top$feature[ord])
  long <- top %>% select(feature, channel, PC1, PC2, PC3) %>% pivot_longer(PC1:PC3, names_to = "pc", values_to = "loading")
  lim <- max(abs(long$loading))
  strip <- ggplot(top, aes("", feature, fill = channel)) + geom_tile() +
    scale_fill_manual(values = CHCOL, name = NULL) + labs(x = NULL, y = NULL) + th +
    theme(axis.text.y = element_text(size = 4.4), axis.ticks = element_blank(), axis.line = element_blank(),
          legend.position = "bottom", legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 6))
  tile <- ggplot(long, aes(pc, feature, fill = loading)) + geom_tile() +
    scale_fill_gradient2(low = "#1f66c2", mid = "white", high = "#c0392b", limits = c(-lim, lim), name = "Loading") +
    labs(x = NULL, y = NULL) + th +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
          legend.position = "bottom", legend.key.height = unit(2.4, "mm"), legend.key.width = unit(6, "mm"))
  strip + tile + plot_layout(widths = c(0.12, 1), guides = "collect") & theme(legend.position = "bottom")
}

enrich <- function(sc, t) {
  universe <- sc$gene
  tg <- t2g_all %>% filter(gene %in% universe) %>% group_by(term) %>% filter(n() < length(universe)) %>% ungroup()
  n10 <- ceiling(0.1 * length(universe))
  sets <- list()
  for (pc in c("PC1", "PC2", "PC3")) {
    o <- sc[order(sc[[pc]], decreasing = TRUE), ]
    sets[[paste0(pc, " (+)")]] <- head(o$gene, n10)
    sets[[paste0(pc, " (-)")]] <- tail(o$gene, n10)
  }
  cc <- compareCluster(sets, fun = "enricher", TERM2GENE = tg, universe = universe, pvalueCutoff = 1,
                       qvalueCutoff = 1, minGSSize = 3, maxGSSize = length(universe))
  res <- as.data.frame(cc) %>% mutate(perturbation = t, n_set = n10, n_universe = length(universe))
  res
}

enrich_plot <- function(res) {
  cols <- c("PC1 (+)", "PC1 (-)", "PC2 (+)", "PC2 (-)", "PC3 (+)", "PC3 (-)")
  lev <- sub(" ", "\n", cols)
  df <- res %>% filter(pvalue < 0.05) %>%
    mutate(col = factor(sub(" ", "\n", as.character(Cluster)), levels = lev), sig = p.adjust < 0.05,
           term = gsub(" (biosynthetic |synthetic )?pathway( event)?$| process$", "", Description),
           genes = sapply(strsplit(geneID, "/"), function(g) if (length(g) <= 2) paste(g, collapse = ", ")
                          else paste0(g[1], ", ... (", length(g), ")")))
  if (!nrow(df)) return(ggplot() + annotate("text", 0, 0, label = "No pathway with p < 0.05", size = 2.4) + theme_void())
  ord <- df %>% group_by(term) %>% summarise(b = min(pvalue)) %>% arrange(desc(b))
  df$term <- factor(df$term, levels = ord$term)
  ggplot(df, aes(col, term)) +
    geom_point(aes(size = Count, fill = -log10(pvalue), colour = sig, stroke = ifelse(sig, 1.0, 0.4)), shape = 21) +
    scale_fill_gradient(low = "#dbe8f8", high = "#0d4a9c", limits = c(1.3, 3), breaks = c(1.5, 2, 2.5, 3),
                        oob = scales::squish, name = "-log10 p") +
    scale_colour_manual(values = c(`FALSE` = "#8a8984", `TRUE` = "#0b0b0b"),
                        labels = c(`FALSE` = "p < 0.05", `TRUE` = "adjusted p < 0.05"), name = NULL, drop = FALSE) +
    scale_size_continuous(range = c(1.4, 4), breaks = c(3, 10, 20), name = "Genes hit") +
    guides(fill = guide_colorbar(barwidth = unit(18, "mm"), barheight = unit(2.4, "mm"))) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = c(0.6, 0.6))) +
    labs(x = NULL, y = NULL) + th +
    theme(axis.text.x = element_text(angle = 0), legend.position = "bottom", legend.box = "vertical",
          legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 6))
}

tail_sets <- function(sc) {
  n10 <- ceiling(0.1 * nrow(sc)); sets <- list()
  for (pc in c("PC1", "PC2", "PC3")) {
    o <- sc[order(sc[[pc]], decreasing = TRUE), ]
    sets[[paste0(pc, " (+)")]] <- head(o$gene, n10)
    sets[[paste0(pc, " (-)")]] <- tail(o$gene, n10)
  }
  sets
}

enrich_comp <- function(sc, t) {
  universe <- sc$gene
  tg <- comp_all %>% filter(gene %in% universe) %>% group_by(term) %>% filter(n() >= 5, n() < length(universe)) %>% ungroup()
  cc <- compareCluster(tail_sets(sc), fun = "enricher", TERM2GENE = tg, universe = universe, pvalueCutoff = 1,
                       qvalueCutoff = 1, minGSSize = 5, maxGSSize = length(universe))
  list(res = as.data.frame(cc) %>% mutate(perturbation = t), tested = tg %>% count(term, name = "n_in_universe"))
}

comp_plot <- function(rc) {
  cols <- sub(" ", "\n", c("PC1 (+)", "PC1 (-)", "PC2 (+)", "PC2 (-)", "PC3 (+)", "PC3 (-)"))
  df <- rc$res %>% mutate(col = factor(sub(" ", "\n", as.character(Cluster)), levels = cols), sig = p.adjust < 0.05)
  df$Description <- factor(df$Description, levels = rc$tested$term[order(rc$tested$n_in_universe)])
  ggplot(df, aes(col, Description)) +
    geom_point(aes(size = Count, fill = -log10(pvalue), colour = sig, stroke = ifelse(sig, 1.0, 0.4)), shape = 21) +
    scale_fill_gradient(low = "#f1f4f9", high = "#0d4a9c", limits = c(0, 3), oob = scales::squish, name = "-log10 p") +
    scale_colour_manual(values = c(`FALSE` = "#8a8984", `TRUE` = "#0b0b0b"),
                        labels = c(`FALSE` = "adjusted p >= 0.05", `TRUE` = "adjusted p < 0.05"), name = NULL, drop = FALSE) +
    scale_size_continuous(range = c(1.2, 4.5), name = "Genes hit") +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) + th +
    guides(fill = guide_colorbar(barwidth = unit(18, "mm"), barheight = unit(2.4, "mm"))) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.key.size = unit(2.4, "mm"),
          legend.text = element_text(size = 6))
}

all_res <- list(); summ <- list(); comp_res <- list()
for (t in c("ORF", "CRISPR")) {
  x <- read_csv(file.path(d, paste0("glyco_profiles_", t, ".csv")), show_col_types = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x$gene
  m <- m[, apply(m, 2, sd) > 0]
  p <- prcomp(m, scale. = TRUE)
  var <- p$sdev^2 / sum(p$sdev^2)
  sc <- data.frame(gene = rownames(p$x), p$x[, 1:10])
  # loadings as in the first paper: eigenvector x SD of the PC (= correlation of the feature with the PC)
  ld <- data.frame(feature = rownames(p$rotation), sweep(p$rotation[, 1:10], 2, p$sdev[1:10], "*"))
  ld$channel <- channel(ld$feature)
  write_csv(sc, file.path(d, paste0("pca_glyco_", t, "_scores.csv")))
  write_csv(ld, file.path(d, paste0("pca_glyco_", t, "_loadings.csv")))
  summ[[t]] <- data.frame(perturbation = t, pc = paste0("PC", 1:10), var_explained = var[1:10])
  res <- enrich(sc, t); all_res[[t]] <- res
  rc <- enrich_comp(sc, t); comp_res[[t]] <- rc$res
  lab <- if (t == "ORF") "OE (ORF)" else "KO (CRISPR)"
  # rows are wrapped separately so the long feature names of c do not shift a and b (no shared axis alignment)
  tag <- function(p, l) p + labs(tag = l)
  top <- (tag(score_plot(sc, var, "PC2", t), "A") + score_plot(sc, var, "PC3", t)) /
    (tag(loading_plot(ld, var, "PC2"), "B") + loading_plot(ld, var, "PC3"))
  bottom <- (wrap_elements(full = heat(ld) + plot_annotation(tag_levels = list(c("C", "")))) |
    tag(comp_plot(rc), "D")) + plot_layout(widths = c(1.7, 1))
  fig <- wrap_elements(full = top) / bottom + plot_layout(heights = c(5.6, 3.4))
  ggsave(file.path(fig_dir, paste0("pca_glyco_", t, ".png")), fig, width = 7.1, height = 9, dpi = 600)
  ggsave(file.path(fig_dir, paste0("pca_glyco_", t, ".pdf")), fig & theme(text = element_text(family = "Helvetica")),
         width = 7.1, height = 9, device = "pdf")
  cat(sprintf("\n=== %s (%s): %d glycogenes, %d features; variance PC1-5: %s\n", lab, t, nrow(m), ncol(m),
              paste(sprintf("%.1f%%", var[1:5] * 100), collapse = ", ")))
  cat("Top loadings (|loading|) PC1:", paste(head(ld$feature[order(-abs(ld$PC1))], 5), collapse = "; "), "\n")
  cat("PC1 score extremes (+):", paste(head(sc$gene[order(-sc$PC1)], 8), collapse = ", "), "\n")
  cat("PC1 score extremes (-):", paste(head(sc$gene[order(sc$PC1)], 8), collapse = ", "), "\n")
  cat("Channel share of squared loadings PC1-3:\n")
  print(ld %>% pivot_longer(PC1:PC3, names_to = "pc") %>% group_by(pc, channel) %>% summarise(s = sum(value^2), .groups = "drop") %>%
          group_by(pc) %>% mutate(s = round(s / sum(s), 3)) %>% pivot_wider(names_from = channel, values_from = s) %>% as.data.frame(), row.names = FALSE)
  cat("Enrichment, p < 0.05:\n")
  print(res %>% filter(pvalue < 0.05) %>% transmute(Cluster, Description, GeneRatio, BgRatio, p = signif(pvalue, 2),
                                                     padj = signif(p.adjust, 2), geneID), row.names = FALSE)
}
write_csv(bind_rows(all_res), file.path(d, "pca_glyco_enrichment.csv"))
write_csv(bind_rows(comp_res), file.path(d, "pca_glyco_compartment_enrichment.csv"))
# supplementary: the former d (pathways), OE and KO side by side
lab2 <- c(ORF = "OE (ORF)", CRISPR = "KO (CRISPR)")
pw <- bind_rows(all_res) %>% filter(pvalue < 0.05) %>%
  mutate(col = factor(sub(" ", "\n", as.character(Cluster)),
                      levels = sub(" ", "\n", c("PC1 (+)", "PC1 (-)", "PC2 (+)", "PC2 (-)", "PC3 (+)", "PC3 (-)"))),
         panel = factor(lab2[perturbation], levels = lab2), sig = p.adjust < 0.05,
         term = gsub(" (biosynthetic |synthetic )?pathway( event)?$| process$", "", Description))
pw$term <- factor(pw$term, levels = (pw %>% group_by(term) %>% summarise(b = min(pvalue)) %>% arrange(desc(b)))$term)
ps <- ggplot(pw, aes(col, term)) +
  geom_point(aes(size = Count, fill = -log10(pvalue), colour = sig, stroke = ifelse(sig, 1.0, 0.4)), shape = 21) +
  facet_grid(. ~ panel, drop = FALSE) +
  scale_fill_gradient(low = "#dbe8f8", high = "#0d4a9c", limits = c(1.3, 3), oob = scales::squish, name = "-log10 p") +
  scale_colour_manual(values = c(`FALSE` = "#8a8984", `TRUE` = "#0b0b0b"),
                      labels = c(`FALSE` = "p < 0.05", `TRUE` = "adjusted p < 0.05"), name = NULL, drop = FALSE) +
  scale_size_continuous(range = c(1.4, 4), name = "Genes hit") +
  scale_x_discrete(drop = FALSE) + labs(x = NULL, y = NULL) + th +
  theme(strip.text = element_text(face = "bold", size = 7), strip.background = element_rect(fill = "#f0efec"),
        legend.position = "bottom", legend.key.size = unit(2.6, "mm"))
hs <- 1.6 + 0.17 * nlevels(droplevels(pw$term))
ggsave(file.path(fig_dir, "pca_glyco_pathway_enrichment_supp.png"), ps, width = 7.1, height = hs, dpi = 300)
ggsave(file.path(fig_dir, "pca_glyco_pathway_enrichment_supp.pdf"), ps + theme(text = element_text(family = "Helvetica")),
       width = 7.1, height = hs, device = "pdf")
write_csv(bind_rows(summ), file.path(d, "pca_glyco_variance.csv"))
