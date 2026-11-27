suppressPackageStartupMessages({library(nlme); library(limma); library(ggplot2)})
out <- "Myeloid/disease_pathways/pathway_pseudotime"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
export <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
x <- readRDS("Myeloid/disease_pathways/donor_bin_expression.rds")
d <- x$samples
d$status <- factor(d$status, levels = c("Healthy", "COPD"))
d$source <- factor(sub(":.*$", "", d$donor))
d$donor <- factor(d$donor)
stopifnot(identical(d$sample, colnames(x$E)), !anyDuplicated(paste(d$donor, d$bin)))
candidates <- read.csv("Myeloid/results_cleaned/coloc_overall_gene_tests.csv")$feature
read_gmt <- function(file) {
  lines <- strsplit(readLines(file, warn = FALSE), "\t")
  setNames(lapply(lines, function(z) unique(z[-c(1, 2)])), vapply(lines, "[", "", 1))
}
sets <- c(read_gmt("Myeloid/disease_pathways/MSigDB_Hallmark_2020.gmt"),
          read_gmt("Myeloid/disease_pathways/Reactome_2022.gmt"))
sets <- lapply(sets, function(g) setdiff(intersect(g, rownames(x$E)), candidates))
sets <- sets[lengths(sets) >= 15 & lengths(sets) <= 500]
z <- t(scale(t(x$E)))
z[!is.finite(z)] <- 0
scores <- t(vapply(sets, function(g) colMeans(z[g, , drop = FALSE]), numeric(ncol(z))))
colnames(scores) <- colnames(x$E)
export(data.frame(pathway = names(sets), n_genes = lengths(sets)), "pathway_gene_coverage")
export(data.frame(pathway = rownames(scores), scores, check.names = FALSE), "donor_bin_pathway_scores")
knots <- as.numeric(quantile(d$pseudotime, c(1/3, 2/3)))
boundary <- range(d$pseudotime)
basis <- splines::ns(d$pseudotime, knots = knots, Boundary.knots = boundary)
d[, c("b1", "b2", "b3")] <- basis
formula <- score ~ source + status * (b1 + b2 + b3)
design <- model.matrix(formula, transform(d, score = 0))
stopifnot(qr(design)$rank == ncol(design))
export(d, "donor_bin_design")

# Exploratory screen across all measured Hallmark and Reactome pathways.
correlation <- duplicateCorrelation(scores, design, block = d$donor)
fit <- eBayes(lmFit(scores, design, block = d$donor, correlation = correlation$consensus.correlation))
interaction <- grep("^statusCOPD:", colnames(design))
disease <- c(match("statusCOPD", colnames(design)), interaction)
for (test in c("shape", "any_disease")) {
  coefs <- if (test == "shape") interaction else disease
  table <- topTable(fit, coef = coefs, number = Inf, sort.by = "F")
  table$pathway <- rownames(table)
  export(table, paste0("all_pathways_", test, "_tests"))
}

# Predefined antigen-presentation pathways, each with its own donor random intercept.
ids <- c("R-HSA-983170", "R-HSA-2132295", "R-HSA-1236975", "R-HSA-877300")
focused <- names(sets)[vapply(names(sets), function(n) any(vapply(ids, grepl, logical(1), x = n)), logical(1))]
stopifnot(length(focused) == 4)
results <- list()
curves <- differences <- observations <- list()
models <- list()
test_linear <- function(model, columns) {
  beta <- fixef(model)[columns]
  covariance <- vcov(model)[columns, columns, drop = FALSE]
  chisq <- as.numeric(crossprod(beta, solve(covariance, beta)))
  q <- length(columns)
  denominator <- unique(model$fixDF$X[columns])
  if (length(denominator) == 1L) {
    return(c(F = chisq/q, p = pf(chisq/q, q, denominator, lower.tail = FALSE),
             numDF = q, denDF = denominator))
  }
  c(F = chisq/q, p = pchisq(chisq, q, lower.tail = FALSE), numDF = q, denDF = Inf)
}

grid <- seq(max(tapply(d$pseudotime, d$status, min)),
            min(tapply(d$pseudotime, d$status, max)), length.out = 150)
grid_basis <- splines::ns(grid, knots = knots, Boundary.knots = boundary)
for (path in focused) {
  message("Focused pathway model: ", path)
  d$score <- as.numeric(scores[path, ])
  model <- lme(formula, random = ~ 1 | donor, data = d, method = "REML",
               control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
  names_coef <- names(fixef(model))
  interaction <- grep("^statusCOPD:", names_coef)
  disease <- c(match("statusCOPD", names_coef), interaction)
  shape <- test_linear(model, interaction)
  any <- test_linear(model, disease)
  results[[path]] <- data.frame(pathway = path, n_genes = length(sets[[path]]),
    n_COPD = length(unique(d$donor[d$status == "COPD"])),
    n_Healthy = length(unique(d$donor[d$status == "Healthy"])),
    shape_F = shape["F"], shape_p = shape["p"], shape_numDF = shape["numDF"], shape_denDF = shape["denDF"],
    any_disease_chisq = any["F"] * any["numDF"], any_disease_p = any["p"],
    random_intercept_sd = as.numeric(VarCorr(model)[1, "StdDev"]),
    residual_sd = model$sigma)
  models[[path]] <- model
  curve <- list()
  for (status in levels(d$status)) {
    matrices <- lapply(levels(d$source), function(source) {
      new <- data.frame(source = factor(source, levels = levels(d$source)),
                         status = factor(status, levels = levels(d$status)),
                         b1 = grid_basis[, 1], b2 = grid_basis[, 2], b3 = grid_basis[, 3])
      model.matrix(~ source + status * (b1 + b2 + b3), new)
    })
    X <- Reduce("+", matrices) / length(matrices)
    estimate <- as.numeric(X %*% fixef(model))
    se <- sqrt(rowSums((X %*% vcov(model)) * X))
    curve[[status]] <- list(X = X, estimate = estimate)
    curves[[paste(path, status)]] <- data.frame(pathway = path, status = status, pseudotime = grid,
      fitted_score = estimate, ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)
  }
  delta_X <- curve$COPD$X - curve$Healthy$X
  estimate <- as.numeric(delta_X %*% fixef(model))
  se <- sqrt(rowSums((delta_X %*% vcov(model)) * delta_X))
  differences[[path]] <- data.frame(pathway = path, pseudotime = grid,
    difference_COPD = estimate, ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)
  observations[[path]] <- data.frame(pathway = path, d[, c("sample", "donor", "source", "status", "pseudotime", "score")])
}
results <- do.call(rbind, results)
results$FDR_shape <- p.adjust(results$shape_p, "BH")
results$FDR_any_disease <- p.adjust(results$any_disease_p, "BH")
export(results, "antigen_pathway_disease_pseudotime_tests")
export(do.call(rbind, curves), "antigen_pathway_fitted_curves")
export(do.call(rbind, differences), "antigen_pathway_curve_differences")
export(do.call(rbind, observations), "antigen_pathway_donor_bin_observations")
saveRDS(models, file.path(out, "antigen_pathway_models.rds"))
saveRDS(list(scores = scores[focused, , drop = FALSE], samples = d,
             knots = knots, boundary = boundary, focused = focused), file.path(out, "focused_model_input.rds"))
writeLines(c(
  "Primary predefined tests: four antigen pathways, each with a donor random-intercept model.",
  "Fixed effects: data source + disease + natural spline(pseudotime) + disease:spline.",
  "Shape difference: joint 3-df disease:spline test.",
  "Any disease difference: joint 4-df Wald chi-square, disease main effect plus disease:spline.",
  "BH correction is applied across the four pathways for each predefined test family.",
  "All 15 coloc candidates are excluded from pathway scores.",
  "Curves average fixed effects equally across the four sources; CI is pointwise.",
  "Gene standardization is shared across disease groups. Per-group scaling is not applied.",
  paste("Pseudotime knots:", paste(knots, collapse = ", ")),
  paste("Exploratory all-pathway donor correlation:", correlation$consensus.correlation)),
  file.path(out, "methods.txt"))
print(results[, c("pathway", "shape_p", "FDR_shape", "any_disease_p", "FDR_any_disease")], row.names = FALSE)
cat("Completed pathway-level disease-pseudotime tests\n")
source(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_all_pathway_mixed_models.R"))
source(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_pathway_source_curve_sensitivity.R"))
source(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_pathway_screen_curves.R"))
status <- system2("python", shQuote(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_pathway_pseudotime.py")),
                  stdout = file.path(out, "plots.log"), stderr = file.path(out, "plots.log"))
stopifnot(status == 0L)
