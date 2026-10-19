#!/usr/bin/env Rscript
# Export a Seurat object's expression matrix and cell metadata to plain files for summarize_h5ad.py.
#
#   Rscript export_seurat.R Daniocell2023_SeuratV4.rds <outdir>
#
# Uses raw counts when any assay stores them; otherwise the normalized values (the 'data' slot,
# log1p of counts per 10,000), which summarize_h5ad.py detects and converts back. Writes to <outdir>:
# counts_x.bin (float64), counts_i.bin (int32, 0-based row index), counts_p.bin (int32 column
# pointers) of the genes x cells sparse matrix, dims.txt (genes cells nonzeros), genes.txt,
# cells.txt, meta.tsv (cell metadata, first column 'cell') and source.txt (assay and slot used).
suppressPackageStartupMessages({ library(Matrix); library(SeuratObject) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript export_seurat.R <object.rds> <outdir>")
dir.create(args[2], showWarnings = FALSE, recursive = TRUE)
obj <- readRDS(args[1])

get_slot <- function(assay, slot) {
  m <- tryCatch({
    if (inherits(assay, "Assay5")) LayerData(assay, layer = slot) else methods::slot(assay, slot)
  }, error = function(e) NULL)
  if (is.null(m) || length(dim(m)) != 2) return(NULL)
  m
}

cat("Assays in the object (genes x cells per slot):\n")
choice <- NULL
for (slot in c("counts", "data")) {
  for (a in unique(c(DefaultAssay(obj), names(obj@assays)))) {
    m <- get_slot(obj@assays[[a]], slot)
    d <- if (is.null(m)) c(0, 0) else dim(m)
    cat(sprintf("  %-10s %-7s %7d x %d\n", a, slot, d[1], d[2]))
    if (is.null(choice) && !is.null(m) && d[1] > 0 && d[2] > 0) choice <- list(assay = a, slot = slot, m = m)
  }
}
if (is.null(choice)) stop("no assay has a non-empty counts or data matrix")
cat(sprintf("Using assay '%s', slot '%s'%s\n", choice$assay, choice$slot,
            if (choice$slot == "data") " (normalized values; converted back to counts per 10,000 downstream)" else ""))
counts <- as(choice$m, "CsparseMatrix")
cat(sprintf("%d genes x %d cells, %.0f nonzero values\n", nrow(counts), ncol(counts), length(counts@x)))

out <- function(f) file.path(args[2], f)
writeBin(as.double(counts@x), out("counts_x.bin"), size = 8)
writeBin(as.integer(counts@i), out("counts_i.bin"), size = 4)
writeBin(as.integer(counts@p), out("counts_p.bin"), size = 4)
writeLines(as.character(c(nrow(counts), ncol(counts), length(counts@x))), out("dims.txt"))
writeLines(rownames(counts), out("genes.txt"))
writeLines(colnames(counts), out("cells.txt"))
writeLines(c(choice$assay, choice$slot), out("source.txt"))
meta <- data.frame(cell = rownames(obj@meta.data), obj@meta.data, check.names = FALSE)
write.table(meta, out("meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Cell metadata columns (name: number of values, examples):\n")
for (col in colnames(obj@meta.data)) {
  u <- unique(as.character(obj@meta.data[[col]]))
  cat(sprintf("  %-28s %7d   %s\n", col, length(u), paste(head(u, 4), collapse = ", ")))
}
