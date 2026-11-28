# Dot plot of the significant GO BP terms (after clusterProfiler::simplify) from go_bp_all_genes.R.
# Facets: ConA (ER) / WGA (AGP); x: perturbation x direction (empty columns kept); y: up to TOPN simplified terms per
# set, ordered by set and adjusted p. Size: genes hit; fill: -log10 BH-adjusted p; text: first genes and the count.
# No title in the image (figure legends go to the manuscript).
# Usage: Rscript go_dotplot.R <dir with go_bp_significant_simplified.csv> [TOPN]
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
d <- args[1]; topn <- if (length(args) > 1) as.integer(args[2]) else 10
s <- read_csv(file.path(d, "go_bp_significant_simplified.csv"), show_col_types = FALSE)
cols <- c("Overexpression|up", "Overexpression|down", "Knockout|up", "Knockout|down")
nice <- c("Overexpression|up" = "OE\nup", "Overexpression|down" = "OE\ndown", "Knockout|up" = "KO\nup",
          "Knockout|down" = "KO\ndown")  # OE overexpression, KO knockout, arrows brighter / darker (legend)
df <- s %>% mutate(lectin = sub("\\|.*", "", Cluster), col = sub("^[^|]*\\|", "", Cluster)) %>%
  group_by(Cluster) %>% slice_min(p.adjust, n = topn, with_ties = FALSE) %>% ungroup() %>%
  mutate(lectin = factor(ifelse(lectin == "ConA", "ER (ConA)", "AGP (WGA)"), levels = c("ER (ConA)", "AGP (WGA)")),
         col = factor(nice[col], levels = nice[cols]),
         genes = sapply(strsplit(geneID, "/"), function(g) if (length(g) <= 3) paste(g, collapse = ", ")
                        else paste0(paste(g[1:2], collapse = ", "), ", ... (", length(g), ")"))) %>%
  arrange(lectin, col, p.adjust)
df$term <- factor(df$Description, levels = rev(unique(df$Description)))
p <- ggplot(df, aes(col, term)) +
  geom_point(aes(size = Count, fill = -log10(p.adjust)), shape = 21, colour = "#0b0b0b", stroke = 0.8) +
  geom_text(aes(label = genes), hjust = 0, nudge_x = 0.24 + 0.0075 * df$Count, size = 2.0, colour = "#3a3a38", fontface = "italic") +
  facet_grid(. ~ lectin, drop = FALSE) +
  scale_fill_gradient(low = "#bcd4f0", high = "#0d4a9c", limits = c(1.3, 3.5), oob = scales::squish,
                      name = "-log10 adjusted p") +
  scale_size_continuous(range = c(2, 6.5), name = "Genes hit") +
  scale_x_discrete(drop = FALSE, expand = expansion(add = c(0.6, 3.4))) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colorbar(barwidth = unit(28, "mm"), barheight = unit(3, "mm"))) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 8.5),
        strip.background = element_rect(fill = "#f0efec"), axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7), legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))
h <- 2.3 + 0.19 * nlevels(df$term)
ggsave(file.path(d, "go_bp_dotplot.png"), p, width = 7.1, height = h, dpi = 300)
ggsave(file.path(d, "go_bp_dotplot.pdf"), p, width = 7.1, height = h, device = "pdf")
cat("terms plotted:", nlevels(df$term), "\n")
