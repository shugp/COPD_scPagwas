suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
# Use run.py --stage annotation --data-root PROJECT_DIRECTORY.
source(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_annotation.R"))
out <- "Myeloid/results_cleaned"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
export <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)

message("Loading original Myeloid object")
dat <- readRDS("Myeloid/Myeloid.Rds")
meta <- dat[[]]
meta$cluster <- as.character(meta[[selected]])
stopifnot(all(meta$cluster %in% names(subtypes)))
meta$subtype <- unname(subtypes[meta$cluster])
normalized <- GetAssayData(dat, assay = "RNA", layer = "data")
# Assess technical quality within library and broad lineage, preserving low-RNA granulocytes.
qc <- data.frame(cell = rownames(meta), cluster = meta$cluster, subtype = meta$subtype,
                 nCount_RNA = meta$nCount_RNA, nFeature_RNA = meta$nFeature_RNA,
                 percent.mt = meta$percent.mt, stringsAsFactors = FALSE)
qc$min_features <- 300
qc$max_mito <- 15
strata <- split(seq_len(nrow(meta)), interaction(meta$data_source, meta$orig.ident, meta$subtype, drop = TRUE))
for (idx in strata) {
  if (length(idx) < 50L) next
  lf <- log10(meta$nFeature_RNA[idx])
  qc$min_features[idx] <- max(300, 10^(median(lf, na.rm = TRUE) - 3 * mad(lf, na.rm = TRUE)))
  mt <- meta$percent.mt[idx]
  qc$max_mito[idx] <- min(15, median(mt, na.rm = TRUE) + 3 * mad(mt, na.rm = TRUE))
}
qc$low_quality <- !is.finite(qc$nCount_RNA) | !is.finite(qc$nFeature_RNA) |
  !is.finite(qc$percent.mt) | qc$nCount_RNA < 500 |
  qc$nFeature_RNA < qc$min_features | qc$nFeature_RNA >= 7000 |
  qc$percent.mt > qc$max_mito

# Require multiple non-myeloid lineage markers. Surfactant alone is widespread ambient RNA.
lineages <- list(T = c("CD3D", "CD3E", "TRAC"), B = c("MS4A1", "CD79A", "CD79B"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  Stromal = c("COL1A1", "COL1A2", "DCN"), Endothelial = c("VWF", "EMCN", "KDR"))
qc$non_myeloid <- FALSE
for (lineage in names(lineages)) {
  genes <- intersect(lineages[[lineage]], rownames(normalized))
  values <- normalized[genes, , drop = FALSE]
  detected <- Matrix::colSums(values > 0)
  score <- Matrix::colMeans(values)
  flag <- detected >= 2 & score >= 0.75
  qc[[paste0("flag_", lineage)]] <- flag
  qc$non_myeloid <- qc$non_myeloid | flag
}
# Clusters 10 and 15 have mixed lymphoid/myeloid and epithelial/myeloid profiles.
qc$ambiguous_cluster <- qc$cluster %in% c("10", "15")
qc$keep <- !qc$low_quality & !qc$non_myeloid & !qc$ambiguous_cluster
qc$exclusion_reason <- ifelse(qc$low_quality, "low_quality",
  ifelse(qc$non_myeloid, "non_myeloid_or_mixed_lineage",
    ifelse(qc$ambiguous_cluster, "ambiguous_cluster", "retained")))
export(qc, "cell_filter_audit")
export(as.data.frame(table(qc$exclusion_reason)), "filter_summary")
export(aggregate(cbind(before = rep(1, nrow(qc)), retained = as.integer(qc$keep)) ~ cluster + subtype, qc, sum),
       "cluster_filter_summary")
message("Retaining ", sum(qc$keep), " of ", nrow(qc), " cells")
dat <- subset(dat, cells = qc$cell[qc$keep])
meta <- meta[qc$keep, , drop = FALSE]
meta$subtype <- factor(meta$subtype, levels = c("Macrophage", "Monocyte", "DC", "Mast/basophil", "Neutrophil"))
counts <- GetAssayData(dat, assay = "RNA", layer = "counts")
normalized <- GetAssayData(dat, assay = "RNA", layer = "data")
stopifnot(identical(colnames(counts), rownames(meta)))
donor <- as.character(meta$orig.ident)
repeat_libraries <- meta$data_source == "GSE168191"
donor[repeat_libraries] <- sub("_rep[0-9]+\\.h5$", "", donor[repeat_libraries])
meta$donor <- paste(meta$data_source, donor, sep = ":")
donors <- unique(meta[, c("donor", "data_source", "status")])
stopifnot(!anyDuplicated(donors$donor), all(donors$status %in% c("Healthy", "COPD")))
rownames(donors) <- donors$donor
export(as.data.frame(table(meta$subtype)), "subtype_counts")
export(as.data.frame(table(meta$cluster, meta$cell_type_jin)), "cluster_jin_crosswalk")

# Preserve the original embedding for before/after comparison.
umap <- Embeddings(dat, "umap")
write.csv(umap, file.path(out, "umap_coordinates.csv"))
points <- data.frame(UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], subtype = meta$subtype)
centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ subtype, points, median)
p <- ggplot(points, aes(UMAP_1, UMAP_2, color = subtype)) +
  geom_point(size = 0.18, alpha = 0.65) +
  ggrepel::geom_text_repel(data = centers, aes(label = subtype), color = "black", seed = 42, show.legend = FALSE) +
  coord_equal() + theme_classic() + labs(color = NULL, title = "Myeloid")
ggsave(file.path(out, "myeloid_cleaned_umap.png"), p, width = 9, height = 7, dpi = 220)
ggsave(file.path(out, "myeloid_cleaned_umap.pdf"), p, width = 9, height = 7)

genes <- intersect(unique(unlist(annotation_markers)), rownames(normalized))
membership <- Matrix::sparse.model.matrix(~ 0 + subtype, meta)
denominator <- as.numeric(table(meta$subtype))
# Average expression uses the linear normalized scale, then log1p for plotting.
avg <- sweep(as.matrix(expm1(normalized[genes, , drop = FALSE]) %*% membership), 2, denominator, "/")
pct <- 100 * sweep(as.matrix((counts[genes, , drop = FALSE] > 0) %*% membership), 2, denominator, "/")
dot <- expand.grid(gene = genes, subtype = levels(meta$subtype))
dot$mean_log1p_expression <- as.vector(log1p(avg))
dot$percent_expressing <- as.vector(pct)
dot$gene <- factor(dot$gene, levels = genes)
dot$subtype <- factor(dot$subtype, levels = rev(levels(meta$subtype)))
export(dot, "annotation_marker_expression")
p <- ggplot(dot, aes(gene, subtype, size = percent_expressing, color = mean_log1p_expression)) +
  geom_point() + scale_size_area(max_size = 7, limits = c(0, 100)) +
  scale_color_gradient(low = "white", high = "#2166AC") +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid.minor = element_blank()) +
  labs(x = NULL, y = NULL, size = "Expressing (%)", color = "Mean log1p\nexpression")
ggsave(file.path(out, "annotation_markers_dotplot.pdf"), p, width = 13, height = 4.5)
ggsave(file.path(out, "annotation_markers_dotplot.png"), p, width = 13, height = 4.5, dpi = 220)

message("Testing overall Myeloid candidates")
source(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_coloc.R"))
dat$myeloid_cluster <- meta$cluster
dat$myeloid_subtype <- meta$subtype
dat$myeloid_donor <- meta$donor
Idents(dat) <- "myeloid_subtype"
dat@misc$myeloid_resolution <- selected
dat@misc$myeloid_filter <- "cell_filter_audit.csv"
dat@misc$myeloid_analysis <- "Overall Myeloid pseudobulk ~ status; five broad annotations"
message("Saving cleaned object")
temporary <- tempfile("myeloid-", tmpdir = out)
saveRDS(dat, temporary, compress = "gzip")
stopifnot(file.rename(temporary, file.path(out, "Myeloid.cleaned.rds")))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("Completed")
