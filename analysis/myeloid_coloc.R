# Overall cleaned Myeloid, one pseudobulk per biological donor, pooled disease model.
coloc <- read.csv("coloc.CSV", check.names = FALSE)
coloc <- coloc[is.finite(coloc$PP.H4.abf) & coloc$PP.H4.abf > 0.6, , drop = FALSE]
myeloid_types <- c("immune_Alveolarmacrophage", "immune_Interstitialmacrophage",
  "immune_Monocyte-derivedmacrophage", "immune_Inflammatorymonocyte",
  "immune_Monocyte", "immune_cDC2", "immune_moDC", "immune_pDC",
  "immune_Mast", "immune_Proliferating")
coloc <- coloc[coloc$cell_type %in% myeloid_types, -1L, drop = FALSE]
aliases <- c(C10orf11 = "LRMDA")
coloc$feature <- vapply(seq_len(nrow(coloc)), function(i) {
  options <- c(coloc$Gene[i], aliases[coloc$Gene[i]], coloc$gene[i])
  hit <- options[!is.na(options) & options %in% rownames(counts)]
  if (length(hit)) hit[1L] else NA_character_
}, character(1))
export(coloc, "coloc_myeloid_input")

donor_factor <- factor(meta$donor, levels = donors$donor)
membership <- Matrix::sparse.model.matrix(~ 0 + donor_factor)
pb <- as.matrix(counts %*% membership)
colnames(pb) <- donors$donor
donors$n_cells <- as.integer(table(donor_factor))
donors$library_size <- colSums(pb)
donors$included <- donors$n_cells >= 20 & donors$library_size > 0
export(donors, "donors")
pb <- pb[, donors$included, drop = FALSE]
samples <- donors[donors$included, ]
status <- factor(samples$status, levels = c("Healthy", "COPD"))
stopifnot(all(table(status) >= 3))
design <- model.matrix(~ status)
library_size <- colSums(pb)
keep_gene <- rowSums(sweep(pb, 2, library_size / 1e6, "/") >= 1) >= 3 &
  rowSums(pb) >= 15
# Cyclic loess adjusts logCPM composition using the whole expressed transcriptome.
v <- limma::voom(pb[keep_gene, , drop = FALSE], design,
                 lib.size = library_size, normalize.method = "cyclicloess")
fit <- limma::eBayes(limma::lmFit(v, design))
de <- limma::topTable(fit, coef = "statusCOPD", number = Inf,
                      sort.by = "P", confint = TRUE)
de$feature <- rownames(de)
names(de)[names(de) == "adj.P.Val"] <- "FDR_transcriptome"
export(de, "myeloid_overall_DE")
candidate <- unique(coloc[, c("gene", "Gene", "feature")])
candidate <- merge(candidate, de, by = "feature", all.x = TRUE, sort = FALSE)
candidate$FDR_coloc_candidates <- p.adjust(candidate$P.Value, "BH")
candidate$test_status <- ifelse(is.na(candidate$feature), "feature_missing",
                                ifelse(is.na(candidate$P.Value), "low_expression", "tested"))
candidate$detection_COPD <- candidate$detection_Healthy <- NA_real_
for (group in c("COPD", "Healthy")) {
  idx <- which(!is.na(candidate$feature))
  candidate[[paste0("detection_", group)]][idx] <-
    Matrix::rowMeans(counts[candidate$feature[idx], meta$status == group, drop = FALSE] > 0)
}
pcc <- read.csv("pagwas/inte_all_gene_PCC.csv")
names(pcc)[1] <- "feature"
names(pcc)[names(pcc) == "PCC"] <- "scPagwas_allcells_PCC"
candidate <- merge(candidate, pcc, by = "feature", all.x = TRUE, sort = FALSE)
candidate <- candidate[order(candidate$P.Value, na.last = TRUE), ]
export(candidate, "coloc_overall_gene_tests")
evidence <- merge(coloc, candidate, by = c("gene", "Gene", "feature"), all.x = TRUE)
evidence$PP4_gt_0.6 <- evidence$PP.H4.abf > 0.6
evidence$expression_support <- !is.na(evidence$FDR_coloc_candidates) &
  evidence$FDR_coloc_candidates < 0.05
export(evidence[order(evidence$PP.H4.abf, decreasing = TRUE), ],
       "coloc_overall_correspondence")

targets <- intersect(unique(na.omit(coloc$feature)), rownames(v$E))
expr <- as.data.frame(as.table(v$E[targets, , drop = FALSE]), stringsAsFactors = FALSE)
names(expr) <- c("feature", "donor", "log2CPM")
expr$status <- samples$status[match(expr$donor, samples$donor)]
export(expr, "coloc_overall_donor_expression")
p <- ggplot(expr, aes(status, log2CPM, color = status)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 42), size = 1.4) +
  facet_wrap(~ feature, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c(Healthy = "#377EB8", COPD = "#E41A1C")) +
  theme_bw() + theme(legend.position = "none") +
  labs(x = NULL, y = "Donor log2 CPM", title = "Overall Myeloid")
ggsave(file.path(out, "coloc_overall_expression.pdf"), p, width = 12, height = 10)
