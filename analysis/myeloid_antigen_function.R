suppressPackageStartupMessages(library(limma))
out <- "Myeloid/disease_pathways"
export <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
x <- readRDS(file.path(out, "donor_state_expression.rds"))
candidates <- read.csv("Myeloid/results_cleaned/coloc_overall_gene_tests.csv")$feature
read_gmt <- function(path) {
  lines <- strsplit(readLines(path), "\t")
  setNames(lapply(lines, function(x) unique(x[-c(1, 2)])), vapply(lines, "[", "", 1))
}
sets <- read_gmt(file.path(out, "Reactome_2022.gmt"))
selected <- c("R-HSA-2132295", "R-HSA-983170", "R-HSA-1236975", "R-HSA-877300")
axis <- sets[grepl(paste(selected, collapse = "|"), names(sets))]
# Remove every coloc candidate from scores to prevent candidate self-inclusion.
score_genes <- lapply(axis, function(g) setdiff(intersect(g, rownames(x$E)), candidates))
stopifnot(all(lengths(score_genes) >= 10))
z <- t(scale(t(x$E)))
z[!is.finite(z)] <- 0
scores <- t(vapply(score_genes, function(g) colMeans(z[g, , drop = FALSE]), numeric(ncol(z))))
colnames(scores) <- colnames(x$E)
export(do.call(rbind, lapply(names(score_genes), function(n) data.frame(pathway = n, gene = score_genes[[n]]))),
       "antigen_score_genes_excluding_candidates")
score_long <- do.call(rbind, lapply(seq_len(nrow(scores)), function(i) {
  data.frame(pathway = rownames(scores)[i], score = scores[i, ],
             x$samples, row.names = NULL, check.names = FALSE)
}))
export(score_long, "antigen_donor_state_scores")
d <- x$samples
d$status <- factor(d$status, levels = c("Healthy", "COPD"))
d$donor <- factor(d$donor)
d$State <- factor(d$State)
tests <- list()
for (state in levels(d$State)) {
  idx <- which(d$State == state)
  stopifnot(all(table(d$status[idx]) >= 3))
  for (path in rownames(scores)) {
    fit <- lm(scores[path, idx] ~ status, data = d[idx, ])
    s <- summary(fit)$coefficients["statusCOPD", ]
    ci <- confint(fit)["statusCOPD", ]
    tests[[paste(state, path)]] <- data.frame(State = state, pathway = path,
      effect_COPD = s[1], ci_low = ci[1], ci_high = ci[2], p = s[4],
      n_COPD = sum(d$status[idx] == "COPD"), n_Healthy = sum(d$status[idx] == "Healthy"))
  }
}
tests <- do.call(rbind, tests)
tests$FDR <- p.adjust(tests$p, "BH")
export(tests, "antigen_score_disease_tests")
focus <- score_long[score_long$State == "State5" &
                      grepl("R-HSA-983170", score_long$pathway), ]
focus$status <- factor(focus$status, levels = c("Healthy", "COPD"))
sensitivity <- lapply(focus$donor, function(donor) {
  fit <- lm(score ~ status, data = focus[focus$donor != donor, ])
  coefficient <- summary(fit)$coefficients["statusCOPD", ]
  data.frame(excluded_donor = donor, effect_COPD = coefficient[1], p = coefficient[4])
})
export(do.call(rbind, sensitivity), "state5_mhci_donor_sensitivity")
associations <- list()
for (gene in intersect(candidates, rownames(x$E))) {
  d$gene_expression <- x$E[gene, ]
  for (path in rownames(scores)) {
    d$score <- scores[path, ]
    fit <- lm(score ~ donor + State + gene_expression, d)
    s <- summary(fit)$coefficients["gene_expression", ]
    associations[[paste(gene, path)]] <- data.frame(feature = gene, pathway = path,
      slope = s[1], p = s[4], n_donor_states = nrow(d))
  }
}
associations <- do.call(rbind, associations)
associations$FDR <- p.adjust(associations$p, "BH")
export(associations[order(associations$FDR), ], "candidate_antigen_associations")

# Residual coexpression modules, adjusting donor and State before correlation.
design <- model.matrix(~ donor + State, d)
res <- qr.resid(qr(design), t(x$E))
res <- scale(res, center = FALSE, scale = sqrt(colSums(res^2)))
res[!is.finite(res)] <- 0
modules <- list()
for (gene in intersect(candidates, colnames(res))) {
  correlations <- as.numeric(crossprod(res[, gene], res))
  names(correlations) <- colnames(res)
  correlations <- sort(correlations[names(correlations) != gene], decreasing = TRUE)
  module <- head(correlations[correlations > 0], 100)
  modules[[gene]] <- data.frame(candidate = gene, gene = names(module), partial_correlation = as.numeric(module))
}
modules <- do.call(rbind, modules)
export(modules, "candidate_coexpression_top100")
enrichment <- list()
universe <- rownames(x$E)
for (gene in unique(modules$candidate)) {
  module <- modules$gene[modules$candidate == gene]
  for (name in names(sets)) {
    overlap_set <- intersect(sets[[name]], universe)
    if (length(overlap_set) < 15 || length(overlap_set) > 500) next
    overlap <- intersect(module, overlap_set)
    enrichment[[paste(gene, name)]] <- data.frame(candidate = gene, pathway = name,
      module_size = length(module), set_size = length(overlap_set), overlap = length(overlap),
      p = phyper(length(overlap) - 1, length(overlap_set), length(universe) - length(overlap_set),
                 length(module), lower.tail = FALSE), overlap_genes = paste(overlap, collapse = ";"))
  }
}
enrichment <- do.call(rbind, enrichment)
enrichment$FDR <- p.adjust(enrichment$p, "BH")
export(enrichment[order(enrichment$FDR), ], "candidate_module_reactome_enrichment")
cat("Completed antigen scores, candidate associations and residual coexpression modules\n")
