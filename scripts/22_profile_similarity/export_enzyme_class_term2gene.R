# Enzyme class of each GlycoEnzOnto gene from GO annotations, as a TERM2GENE table (user, 2026-11-23: after the
# compartment result, separate "where it works" from "what kind of enzyme it is": "やってみようか").
# A gene belongs to a class if it is annotated to the GO term or any descendant (org.Hs.eg.db, keytype GOALL, any
# evidence). Classes (a gene can be in several):
#   Glycosyltransferase  GO:0016757 glycosyltransferase activity
#   Glycosidase          GO:0016798 hydrolase activity, acting on glycosyl bonds
#   Sulfotransferase     GO:0008146 sulfotransferase activity
#   Sulfatase            GO:0008484 sulfuric ester hydrolase activity
#   Transporter          GO:0022857 transmembrane transporter activity
#   Nucleotide-sugar synthesis  GO:0009226 nucleotide-sugar biosynthetic process
#   Carbohydrate binding GO:0030246 carbohydrate binding
# Output: <out>/enzyme_class_term2gene.tsv (term, gene), enzyme_class_counts.csv.
# Usage: Rscript export_enzyme_class_term2gene.R <glycoenzonto term2gene.tsv> <out dir>
suppressPackageStartupMessages({library(org.Hs.eg.db); library(AnnotationDbi); library(dplyr); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
glyco <- unique(read_tsv(args[1], show_col_types = FALSE)$gene)
out <- args[2]; dir.create(out, showWarnings = FALSE, recursive = TRUE)
classes <- c(`Glycosyltransferase` = "GO:0016757", `Glycosidase` = "GO:0016798", `Sulfotransferase` = "GO:0008146",
             `Sulfatase` = "GO:0008484", `Transporter` = "GO:0022857",
             `Nucleotide-sugar synthesis` = "GO:0009226", `Carbohydrate binding` = "GO:0030246")
t2g <- bind_rows(lapply(names(classes), function(k) {
  g <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = classes[[k]], keytype = "GOALL", columns = "SYMBOL"))
  tibble(term = k, gene = intersect(unique(g$SYMBOL), glyco))
}))
write_tsv(t2g, file.path(out, "enzyme_class_term2gene.tsv"))
n <- t2g %>% count(term, name = "n_genes") %>% arrange(desc(n_genes))
write_csv(n, file.path(out, "enzyme_class_counts.csv"))
print(as.data.frame(n), row.names = FALSE)
cat("glycogenes with a class:", length(unique(t2g$gene)), "of", length(glyco), "\n")
cat("glycogenes in 2+ classes:", sum(table(t2g$gene) > 1), "\n")
