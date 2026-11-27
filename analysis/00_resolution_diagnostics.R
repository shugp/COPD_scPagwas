# Inspect existing resolutions without changing clustering or the input object.
suppressPackageStartupMessages({library(Seurat); library(cluster)})
dat <- readRDS("Myeloid/Myeloid.Rds")
meta <- dat[[]]
out <- "Myeloid/resolution_review"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
resolutions <- grep("^RNA_snn_res\\.", colnames(meta), value = TRUE)
donor <- as.character(meta$orig.ident)
idx <- meta$data_source == "GSE168191"
donor[idx] <- sub("_rep[0-9]+\\.h5$", "", donor[idx])
meta$donor <- paste(meta$data_source, donor, sep = ":")
set.seed(42)
sampled <- sort(sample(seq_len(nrow(meta)), min(3000L, nrow(meta))))
harmony <- Embeddings(dat, "harmony")[sampled, 1:30, drop = FALSE]
distance <- dist(harmony)
rows <- list()
for (resolution in resolutions) {
  labels <- as.character(meta[[resolution]])
  silhouette <- mean(cluster::silhouette(as.integer(factor(labels[sampled])), distance)[, "sil_width"])
  counts <- as.data.frame(table(cluster = labels, cell_type_jin = meta$cell_type_jin))
  counts <- counts[counts$Freq > 0, ]
  write.csv(counts, file.path(out, paste0("jin_", resolution, ".csv")), row.names = FALSE)
  for (label in sort(unique(labels))) {
    d <- meta[labels == label, ]
    jin <- sort(table(d$cell_type_jin), decreasing = TRUE)
    rows[[paste(resolution, label)]] <- data.frame(
      resolution = resolution, cluster = label, n_cells = nrow(d),
      jin_majority = names(jin)[1], jin_purity = as.numeric(jin[1])/nrow(d),
      n_donors = length(unique(d$donor)), n_cohorts = length(unique(d$data_source)),
      max_donor_fraction = max(table(d$donor))/nrow(d),
      max_cohort_fraction = max(table(d$data_source))/nrow(d),
      harmony_silhouette = silhouette)
  }
}
write.csv(do.call(rbind, rows), file.path(out, "resolution_diagnostics.csv"), row.names = FALSE)
writeLines(c("Resolution 0.5 is the existing marker-reviewed selection, not an automatic optimum.",
  "This audit supplements cell_type_jin consistency with source/donor concentration and sampled silhouettes.",
  "Silhouettes use a seeded 3,000-cell sample; values can differ from historical diagnostics."),
  file.path(out, "interpretation.txt"))
