# GlycoEnzOnto pathway enrichment of glycogenes whose ConA / WGA staining differed from the other genes' wells at
# nominal p < 0.05, with clusterProfiler (user, 2026-11-22: "クラスタープロファイラー使ってねRの").
# Sets: perturbation (ORF, CRISPR) x lectin (ConA, WGA) x direction (brighter, darker).
# enricher() = one-sided hypergeometric test, BH; universe = GlycoEnzOnto genes tested in that perturbation;
# gene sets with 3 to (universe size - 1) tested genes (the all-gene root pathway is dropped).
# Each perturbation has its own universe, so compareCluster() is run per perturbation and the results merged.
# Usage: Rscript enrichment_clusterprofiler.R <calls_dir> <term2gene.tsv> <out_dir>
suppressPackageStartupMessages({
  library(clusterProfiler); library(enrichplot); library(ggplot2); library(dplyr); library(readr)
})
args <- commandArgs(trailingOnly = TRUE)
calls_dir <- args[1]; t2g_file <- args[2]; out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
t2g <- read_tsv(t2g_file, show_col_types = FALSE)
label <- c(ORF = "Overexpression", CRISPR = "Knockout")

res <- list(); clusters <- list()
for (t in c("ORF", "CRISPR")) {
  g <- read_csv(file.path(calls_dir, paste0(t, "_glycogene_calls.csv")), show_col_types = FALSE)
  universe <- g$gene
  t2g_t <- t2g %>% filter(gene %in% universe) %>% group_by(term) %>% filter(n() < length(universe)) %>% ungroup()
  sets <- list()
  for (lec in c("ConA", "WGA")) {
    p <- g[[paste0(lec, "_p")]]; d <- g[[paste0(lec, "_diff_sd")]]
    sets[[paste(lec, label[[t]], "brighter", sep = "\n")]] <- g$gene[!is.na(p) & p < 0.05 & d > 0]
    sets[[paste(lec, label[[t]], "darker", sep = "\n")]] <- g$gene[!is.na(p) & p < 0.05 & d < 0]
  }
  cc <- compareCluster(sets, fun = "enricher", TERM2GENE = t2g_t, universe = universe,
                       pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, maxGSSize = length(universe))
  df <- as.data.frame(cc) %>% mutate(perturbation = t, n_universe = length(universe))
  res[[t]] <- df; clusters <- c(clusters, sets)
}
all <- bind_rows(res)
write_csv(all, file.path(out_dir, "clusterprofiler_enrichment.csv"))

# Dot plot (ggplot on the compareCluster results): facets = lectin, x = perturbation x direction, y = term.
# Terms with p < 0.05 in any set; terms with the same hit genes in the same set are merged into one row.
# Size = genes hit, colour = -log10 p, black ring + bold label = BH-adjusted p < 0.1, text = hit genes.
keep <- unique(all$Description[all$pvalue < 0.05])
short <- function(x) gsub(" (biosynthetic |synthetic )?pathway( event)?$| process$", "", x)
plot_df <- all %>% filter(pvalue < 0.05) %>%
  mutate(lectin = sub("\n.*", "", Cluster),
         column = factor(sub("^[^\n]*\n", "", Cluster),
                         levels = c("Overexpression\nbrighter", "Overexpression\ndarker",
                                    "Knockout\nbrighter", "Knockout\ndarker")),
         term = short(Description), sig = p.adjust < 0.1) %>%
  group_by(lectin, column, geneID) %>%
  summarise(term = paste(unique(term), collapse = " /\n"), pvalue = min(pvalue), p.adjust = min(p.adjust),
            Count = first(Count), sig = any(sig), .groups = "drop") %>%
  arrange(lectin, column, pvalue) %>%
  mutate(term = factor(term, levels = rev(unique(term))), genes = ifelse(sig, sprintf("%s  (adjusted p = %.2f)", gsub("/", ", ", geneID), p.adjust),
                        gsub("/", ", ", geneID)))
# publication layout (7.1 in wide like Figure 1; facets named after the channels; columns OE/KO x up/down)
col_lab <- c("Overexpression\nbrighter" = "OE\nup", "Overexpression\ndarker" = "OE\ndown",
             "Knockout\nbrighter" = "KO\nup", "Knockout\ndarker" = "KO\ndown")
plot_df <- plot_df %>%
  mutate(column = factor(col_lab[as.character(column)], levels = col_lab),
         lectin = factor(ifelse(lectin == "ConA", "ER (ConA)", "AGP (WGA)"), levels = c("ER (ConA)", "AGP (WGA)")))
p <- ggplot(plot_df, aes(column, term)) +
  geom_point(aes(size = Count, fill = -log10(pvalue), colour = sig, stroke = ifelse(sig, 1.4, 0.4)), shape = 21) +
  geom_text(aes(label = genes, fontface = ifelse(sig, "bold", "plain")), hjust = 0, nudge_x = 0.13 + 0.045 * plot_df$Count + 0.08 * plot_df$sig, size = 2.0,
            colour = "#3a3a38") +
  facet_grid(. ~ lectin, drop = FALSE) +
  scale_fill_gradient(low = "#dbe8f8", high = "#0d4a9c", limits = c(1.3, 3), oob = scales::squish,
                      name = "-log10 p") +
  scale_colour_manual(values = c(`FALSE` = "#8a8984", `TRUE` = "#0b0b0b"),
                      labels = c(`FALSE` = "p < 0.05", `TRUE` = "BH-adjusted p < 0.1"), name = NULL) +
  scale_size_continuous(range = c(1.8, 5), breaks = c(1, 2, 3, 5), name = "Genes hit") +
  scale_x_discrete(drop = FALSE, expand = expansion(add = c(0.6, 3.6))) +
  labs(x = NULL, y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3.5, stroke = c(0.4, 1.4), fill = "#7fa6dc"), order = 1),
         size = guide_legend(order = 2),
         fill = guide_colorbar(barwidth = unit(25, "mm"), barheight = unit(2.6, "mm"), order = 3)) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 7.5),
        strip.background = element_rect(fill = "#f0efec"), axis.text.y = element_text(size = 6.5, colour = "#111111"),
        axis.text.x = element_text(size = 7, colour = "#111111"), legend.position = "bottom",
        legend.box = "vertical", legend.margin = margin(0, 6, 0, 6), legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 6.8))
h <- 1.9 + 0.2 * nlevels(plot_df$term)
ggsave(file.path(out_dir, "clusterprofiler_dotplot.png"), p, width = 7.1, height = h, dpi = 300)
ggsave(file.path(out_dir, "clusterprofiler_dotplot.pdf"), p + theme(text = element_text(family = "Helvetica")),
       width = 7.1, height = h, device = "pdf")  # pdf() has no Arial; Helvetica is metrically the same
lev <- names(clusters)

cat(sprintf("sets: %s\n", paste(sprintf("%s=%d", gsub("\n", " ", lev), lengths(clusters)), collapse = "; ")))
cat(sprintf("terms with p < 0.05: %d; adjusted p < 0.1: %d\n", length(keep), sum(all$p.adjust < 0.1)))
top <- all %>% filter(pvalue < 0.05) %>% arrange(pvalue) %>%
  transmute(Cluster = gsub("\n", " ", Cluster), Description, GeneRatio, BgRatio,
            pvalue = signif(pvalue, 2), p.adjust = signif(p.adjust, 2), geneID)
write.table(top, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
