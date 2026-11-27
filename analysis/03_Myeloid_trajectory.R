suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
out <- "Myeloid/trajectory"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
export <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
prefix <- tempfile("myeloid-trajectory-")
  dat <- readRDS("Myeloid/results_cleaned/Myeloid.cleaned.rds")
  meta <- dat[[]]
  meta$cell <- rownames(meta)
  write.csv(meta[, c("cell", "orig.ident", "myeloid_subtype", "myeloid_cluster", "myeloid_donor", "status", "data_source")],
            paste0(prefix, "_meta.csv"), row.names = FALSE)
  write.csv(Embeddings(dat, "harmony"), paste0(prefix, "_harmony.csv"))
  write.csv(Embeddings(dat, "umap"), paste0(prefix, "_umap.csv"))
  counts <- GetAssayData(dat, layer = "counts")
  normalized <- GetAssayData(dat, layer = "data")
  status <- system2("python", shQuote(c(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_trajectory.py"), prefix, out)),
                    stdout = file.path(out, "graph.log"), stderr = file.path(out, "graph.log"))
  stopifnot(status == 0L)
  unlink(paste0(prefix, c("_meta.csv", "_harmony.csv", "_umap.csv")))
meta <- read.csv(file.path(out, "cell_pseudotime_state.csv"))
stopifnot(identical(meta$cell, colnames(counts)))
meta$donor <- meta$myeloid_donor
active <- which(meta$trajectory_scope == "monocyte_macrophage" & is.finite(meta$pseudotime))
candidates <- read.csv("Myeloid/results_cleaned/coloc_overall_gene_tests.csv")
targets <- unique(na.omit(candidates$feature))
meta$State <- factor(meta$State, levels = paste0("State", sort(as.integer(sub("State", "", unique(na.omit(meta$State)))))))

# Pseudobulk raw counts for each donor and graph State.
group <- interaction(meta$donor[active], meta$State[active], drop = TRUE)
membership <- Matrix::sparseMatrix(i = active, j = as.integer(group), x = 1,
                                   dims = c(ncol(counts), nlevels(group)))
pb <- as.matrix(counts %*% membership)
colnames(pb) <- levels(group)
first <- match(levels(group), as.character(group))
samples <- meta[active[first], c("donor", "State", "status")]
samples$n_cells <- tabulate(as.integer(group), nlevels(group))
samples$included <- samples$n_cells >= 20
coverage <- table(samples$State[samples$included])
keep_states <- names(coverage[coverage >= 3])
samples$included <- samples$included & samples$State %in% keep_states
export(samples, "state_pseudobulk_coverage")
pb <- pb[, samples$included, drop = FALSE]
samples <- droplevels(samples[samples$included, ])
samples$donor <- factor(samples$donor)
design <- model.matrix(~ donor + State, samples)
stopifnot(qr(design)$rank == ncol(design), nrow(design) > ncol(design))
lib <- colSums(pb)
expressed <- rowSums(sweep(pb, 2, lib / 1e6, "/") >= 1) >= 3 & rowSums(pb) >= 15
v <- limma::voom(pb[expressed, , drop = FALSE], design,
                 lib.size = lib, normalize.method = "cyclicloess")
fit <- limma::eBayes(limma::lmFit(v, design))
state_coef <- grep("^State", colnames(design))
stats <- limma::topTable(fit, coef = state_coef, number = Inf, sort.by = "F")
stats$feature <- rownames(stats)
names(stats)[names(stats) == "adj.P.Val"] <- "FDR_transcriptome"
export(stats, "state_allgene_tests")
hits <- merge(candidates[, c("Gene", "feature")], stats, by = "feature", all.x = TRUE)
hits$FDR_candidates <- p.adjust(hits$P.Value, "BH")
hits$test <- "donor-adjusted State omnibus"
hits <- hits[order(hits$P.Value), ]
export(hits, "candidate_state_tests")
pairwise <- do.call(rbind, lapply(state_coef, function(coef) {
  t <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none", confint = TRUE)
  t$feature <- rownames(t)
  t <- t[t$feature %in% targets, ]
  t$comparison <- paste0(sub("^State", "", colnames(design)[coef]), " vs ", levels(samples$State)[1])
  t
}))
pairwise$FDR_all_candidate_contrasts <- p.adjust(pairwise$P.Value, "BH")
export(pairwise, "candidate_state_contrasts")
heat <- pairwise
heat$State <- sub(" vs .*", "", heat$comparison)
heat$significance <- ifelse(heat$FDR_all_candidate_contrasts < 0.05, "*", "")
p <- ggplot(heat, aes(State, feature, fill = logFC)) +
  geom_tile() + geom_text(aes(label = significance), size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_minimal() + labs(x = NULL, y = NULL, fill = "Adjusted log2FC",
                        title = "State contrasts relative to State1",
                        subtitle = "Donor fixed effects; * FDR < 0.05 across candidate contrasts")
ggsave(file.path(out, "candidate_state_heatmap.png"), p, width = 8, height = 7, dpi = 220)
ggsave(file.path(out, "candidate_state_heatmap.pdf"), p, width = 8, height = 7)

state_expression <- as.data.frame(as.table(v$E[intersect(targets, rownames(v$E)), , drop = FALSE]))
names(state_expression) <- c("feature", "sample", "log2CPM")
state_expression$State <- samples$State[match(state_expression$sample, colnames(v$E))]
state_expression$donor <- as.character(samples$donor[match(state_expression$sample, colnames(v$E))])
# voom preserves matrix column names, inherited from the membership matrix.
stopifnot(!anyNA(state_expression$State))
export(state_expression, "candidate_state_donor_expression")
p <- ggplot(state_expression, aes(State, log2CPM)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 42), size = 0.7, alpha = 0.5, color = "gray40") +
  facet_wrap(~ feature, scales = "free_y", ncol = 3) + theme_classic(base_size = 10) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 50, hjust = 1)) +
  labs(x = NULL, y = "Donor-State log2 CPM")
ggsave(file.path(out, "candidate_state_expression.pdf"), p, width = 13, height = 13)

# Continuous dynamics, aggregating into eight global pseudotime quantiles within donor.
breaks <- unique(quantile(meta$pseudotime[active], seq(0, 1, length.out = 9)))
bins <- cut(meta$pseudotime[active], breaks, include.lowest = TRUE)
group <- interaction(meta$donor[active], bins, drop = TRUE)
membership <- Matrix::sparseMatrix(i = active, j = as.integer(group), x = 1,
                                   dims = c(ncol(counts), nlevels(group)))
pb <- as.matrix(counts %*% membership)
colnames(pb) <- levels(group)
first <- match(levels(group), as.character(group))
samples <- meta[active[first], c("donor", "status")]
samples$n_cells <- tabulate(as.integer(group), nlevels(group))
samples$pseudotime <- as.numeric(tapply(meta$pseudotime[active], group, mean))
samples$included <- samples$n_cells >= 20
export(samples, "pseudotime_bin_coverage")
pb <- pb[, samples$included, drop = FALSE]
samples <- samples[samples$included, ]
samples$donor <- factor(samples$donor)
design <- model.matrix(~ donor + splines::ns(pseudotime, df = 3), samples)
stopifnot(qr(design)$rank == ncol(design))
lib <- colSums(pb)
expressed <- rowSums(sweep(pb, 2, lib / 1e6, "/") >= 1) >= 3 & rowSums(pb) >= 15
v <- limma::voom(pb[expressed, , drop = FALSE], design,
                 lib.size = lib, normalize.method = "cyclicloess")
fit <- limma::eBayes(limma::lmFit(v, design))
coef <- grep("splines", colnames(design))
trend <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "F")
trend$feature <- rownames(trend)
trend <- merge(candidates[, c("Gene", "feature")], trend, by = "feature", all.x = TRUE)
trend$FDR_candidates <- p.adjust(trend$P.Value, "BH")
export(trend[order(trend$P.Value), ], "candidate_pseudotime_tests")
expr <- as.data.frame(as.table(v$E[intersect(targets, rownames(v$E)), , drop = FALSE]))
names(expr) <- c("feature", "sample", "log2CPM")
expr$pseudotime <- samples$pseudotime[match(expr$sample, colnames(v$E))]
expr$donor <- as.character(samples$donor[match(expr$sample, colnames(v$E))])
stopifnot(!anyNA(expr$pseudotime))
export(expr, "candidate_pseudotime_donor_expression")
p <- ggplot(expr, aes(pseudotime, log2CPM)) +
  geom_point(size = 0.8, alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ splines::ns(x, df = 3), color = "#2166AC") +
  facet_wrap(~ feature, scales = "free_y", ncol = 3) + theme_classic(base_size = 10) +
  labs(x = "Diffusion pseudotime", y = "Donor-bin log2 CPM")
ggsave(file.path(out, "candidate_pseudotime_curves.pdf"), p, width = 12, height = 12)

# Expression marker checks for the inferred graph States.
genes <- intersect(c("FCN1", "VCAN", "S100A8", "FCGR3A", "C1QC", "CD68", "FABP4",
                     "MARCO", "FOLR2", "TREM2", "MKI67"), rownames(normalized))
group <- droplevels(meta$State[active])
membership <- Matrix::sparseMatrix(i = active, j = as.integer(group), x = 1,
                                   dims = c(ncol(counts), nlevels(group)))
avg <- sweep(as.matrix(expm1(normalized[genes, , drop = FALSE]) %*% membership),
             2, tabulate(as.integer(group)), "/")
pct <- 100 * sweep(as.matrix((counts[genes, , drop = FALSE] > 0) %*% membership),
                   2, tabulate(as.integer(group)), "/")
dot <- expand.grid(gene = genes, State = levels(group))
dot$mean_log1p <- as.vector(log1p(avg))
dot$percent_expressing <- as.vector(pct)
export(dot, "state_marker_expression")
p <- ggplot(dot, aes(factor(gene, levels = genes), State, color = mean_log1p, size = percent_expressing)) +
  geom_point() + scale_size_area(max_size = 7) + scale_color_gradient(low = "white", high = "#2166AC") +
  theme_classic(base_size = 10) + labs(x = NULL, y = NULL, size = "Expressing (%)", color = "Mean log1p")
ggsave(file.path(out, "state_markers.pdf"), p, width = 10, height = 5)
summary <- merge(hits[, c("Gene", "feature", "P.Value", "FDR_candidates")],
                 trend[, c("feature", "P.Value", "FDR_candidates")],
                 by = "feature", suffixes = c("_State", "_pseudotime"), all = TRUE)
export(summary, "candidate_trajectory_summary")
cat("Completed donor-adjusted State and pseudotime tests\n")
