suppressPackageStartupMessages({
  library(Seurat)
  library(limma)
})
out <- "Myeloid/disease_pathways"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
export <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
dat <- readRDS("Myeloid/results_cleaned/Myeloid.cleaned.rds")
counts <- GetAssayData(dat, layer = "counts")
if (exists("dat")) {
  write.csv(Embeddings(dat, "umap"), "Myeloid/results_cleaned/umap_coordinates.csv")
}
meta <- read.csv("Myeloid/trajectory/cell_pseudotime_state.csv")
stopifnot(identical(meta$cell, colnames(counts)))
active <- which(meta$trajectory_scope == "monocyte_macrophage")
meta$donor <- meta$myeloid_donor
candidates <- read.csv("Myeloid/results_cleaned/coloc_overall_gene_tests.csv")
targets <- candidates$feature
axis_genes <- c("HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRB1", "HLA-DRA",
  "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DRB5", "CD74", "CIITA",
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "CTSS", "IFI30",
  "B2M", "HLA-A", "HLA-B", "HLA-C", "ERAP1", "ERAP2", "TAP1", "TAP2",
  "TAPBP", "CALR", "CANX", "PDIA3", "NLRC5", "PSMB8", "PSMB9", "PSMB10")
prepare_voom <- function(pb, design) {
  lib <- colSums(pb)
  expressed <- rowSums(sweep(pb, 2, lib / 1e6, "/") >= 1) >= 3 & rowSums(pb) >= 15
  stopifnot(qr(design)$rank == ncol(design), nrow(design) > ncol(design))
  voom(pb[expressed, , drop = FALSE], design, lib.size = lib, normalize.method = "cyclicloess")
}
group <- interaction(meta$donor[active], meta$State[active], drop = TRUE)
membership <- Matrix::sparseMatrix(i = active, j = as.integer(group), x = 1,
                                   dims = c(ncol(counts), nlevels(group)))
bulk <- as.matrix(counts %*% membership)
colnames(bulk) <- levels(group)
first <- match(levels(group), as.character(group))
samples <- meta[active[first], c("donor", "State", "status")]
samples$n_cells <- tabulate(as.integer(group), nlevels(group))
samples$sample <- colnames(bulk)
samples$included <- samples$n_cells >= 20
export(samples, "donor_state_coverage")
results <- list()
for (state in sort(unique(samples$State))) {
  keep <- samples$included & samples$State == state
  d <- droplevels(samples[keep, ])
  d$status <- factor(d$status, levels = c("Healthy", "COPD"))
  n <- table(d$status)
  if (any(n < 3)) next
  message("Within-State disease comparison: ", state)
  design <- model.matrix(~ status, d)
  v <- prepare_voom(bulk[, keep, drop = FALSE], design)
  fit <- eBayes(lmFit(v, design))
  t <- topTable(fit, coef = "statusCOPD", number = Inf, sort.by = "P", confint = TRUE)
  t$feature <- rownames(t)
  t$State <- state
  t$n_COPD <- as.integer(n["COPD"])
  t$n_Healthy <- as.integer(n["Healthy"])
  results[[state]] <- t
}
within <- do.call(rbind, results)
export(within, "within_state_allgenes")
candidate_tests <- within[within$feature %in% targets, ]
candidate_tests$FDR_candidates_all_states <- p.adjust(candidate_tests$P.Value, "BH")
export(candidate_tests, "within_state_candidates")
axis_tests <- within[within$feature %in% axis_genes, ]
axis_tests$FDR_axis_all_states <- p.adjust(axis_tests$P.Value, "BH")
export(axis_tests, "within_state_antigen_genes")

# Paired donor fixed-effect State contrasts provide signed transcriptome ranks.
d <- samples[samples$included, ]
d$donor <- factor(d$donor)
d$State <- factor(d$State)
design <- model.matrix(~ donor + State, d)
v <- prepare_voom(bulk[, samples$included, drop = FALSE], design)
fit <- eBayes(lmFit(v, design))
state_ranks <- lapply(grep("^State", colnames(design)), function(j) {
  t <- topTable(fit, coef = j, number = Inf, sort.by = "none")
  data.frame(feature = rownames(t), t = t$t, logFC = t$logFC,
             contrast = paste0(sub("^State", "", colnames(design)[j]), "_vs_State1"))
})
export(do.call(rbind, state_ranks), "state_transition_ranks")
saveRDS(list(E = v$E, samples = d, features = rownames(v$E), axis_genes = axis_genes),
        file.path(out, "donor_state_expression.rds"))

# Donor-bin repeated-measures model with a disease-by-pseudotime interaction.
breaks <- unique(quantile(meta$pseudotime[active], seq(0, 1, length.out = 9)))
bins <- cut(meta$pseudotime[active], breaks, include.lowest = TRUE)
group <- interaction(meta$donor[active], bins, drop = TRUE)
membership <- Matrix::sparseMatrix(i = active, j = as.integer(group), x = 1,
                                   dims = c(ncol(counts), nlevels(group)))
bulk <- as.matrix(counts %*% membership)
colnames(bulk) <- levels(group)
first <- match(levels(group), as.character(group))
d <- meta[active[first], c("donor", "status")]
d$bin <- as.character(bins[first])
d$n_cells <- tabulate(as.integer(group), nlevels(group))
d$pseudotime <- as.numeric(tapply(meta$pseudotime[active], group, mean))
d$sample <- colnames(bulk)
eligible <- d$n_cells >= 20
coverage <- table(d$bin[eligible], d$status[eligible])
shared <- rownames(coverage)[rowSums(coverage >= 3) == 2]
d$included <- eligible & d$bin %in% shared
export(d, "interaction_bin_coverage")
bulk <- bulk[, d$included, drop = FALSE]
d <- d[d$included, ]
d$status <- factor(d$status, levels = c("Healthy", "COPD"))
d$donor <- factor(d$donor)
design <- model.matrix(~ status * splines::ns(pseudotime, df = 3), d)
v <- prepare_voom(bulk, design)
message("Estimating repeated-donor correlation")
correlation <- duplicateCorrelation(v, design, block = d$donor)
fit <- eBayes(lmFit(v, design, block = d$donor, correlation = correlation$consensus.correlation))
saveRDS(list(fit = fit, d = d, v = v), file.path(out, "disease_interaction_fit.rds"))
interaction <- grep("^statusCOPD:", colnames(design))
t <- topTable(fit, coef = interaction, number = Inf, sort.by = "F")
t$feature <- rownames(t)
export(t, "disease_pseudotime_interaction_allgenes")
hits <- t[t$feature %in% targets, ]
hits$FDR_candidates <- p.adjust(hits$P.Value, "BH")
export(hits, "disease_pseudotime_interaction_candidates")
hits <- t[t$feature %in% axis_genes, ]
hits$FDR_axis <- p.adjust(hits$P.Value, "BH")
export(hits, "disease_pseudotime_interaction_antigen_genes")
saveRDS(list(E = v$E, samples = d, correlation = correlation$consensus.correlation),
        file.path(out, "donor_bin_expression.rds"))

observed_features <- intersect(unique(c(targets, axis_genes)), rownames(v$E))
observations <- as.data.frame(as.table(v$E[observed_features, , drop = FALSE]))
names(observations) <- c("feature", "sample", "log2CPM")
observations <- merge(observations, d, by = "sample")
export(observations, "disease_pseudotime_observations")

# Evaluate COPD-Healthy at three observed common-support positions for signed ranks.
times <- quantile(d$pseudotime, c(0.2, 0.5, 0.8))
new <- expand.grid(status = levels(d$status), pseudotime = as.numeric(times))
new$status <- factor(new$status, levels = levels(d$status))
basis <- model.frame(~ status * splines::ns(pseudotime, df = 3), d)
X <- model.matrix(~ status * splines::ns(pseudotime, df = 3),
                  data = model.frame(terms(basis), new, xlev = .getXlevels(terms(basis), basis)))
C <- t(X[c(2, 4, 6), , drop = FALSE] - X[c(1, 3, 5), , drop = FALSE])
colnames(C) <- c("early", "middle", "late")
cf <- eBayes(contrasts.fit(fit, C))
time_results <- lapply(seq_len(3), function(i) {
  x <- topTable(cf, coef = i, number = Inf, sort.by = "P", confint = TRUE)
  x$feature <- rownames(x)
  x$position <- colnames(C)[i]
  x$pseudotime <- times[i]
  x
})
export(do.call(rbind, time_results), "disease_pseudotime_position_tests")
writeLines(c(paste("consensus_donor_correlation", correlation$consensus.correlation),
             paste("common_support_bins", length(shared)),
             "State models use pooled status comparisons within each State.",
             "Interaction uses status * natural spline(pseudotime, df=3), repeated donor block.",
             "No between-cohort validation was performed."),
           file.path(out, "model_parameters.txt"))
cat("Completed disease and interaction models\n")

source(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_antigen_function.R"))
status <- system2("python", shQuote(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_gsea.py")),
                  stdout = file.path(out, "gsea.log"), stderr = file.path(out, "gsea.log"))
stopifnot(status == 0L)
status <- system2("python", shQuote(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "05_Myeloid_plots.py")),
                  stdout = file.path(out, "plots.log"), stderr = file.path(out, "plots.log"))
stopifnot(status == 0L)
source(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "06_Myeloid_result_panels.R"))
status <- system2("python", shQuote(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_result_panels.py")),
                  stdout = file.path(out, "result_panels.log"), stderr = file.path(out, "result_panels.log"))
stopifnot(status == 0L)

source(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "07_Myeloid_differential_association.R"))
source(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "08_Myeloid_pathway_pseudotime.R"))
