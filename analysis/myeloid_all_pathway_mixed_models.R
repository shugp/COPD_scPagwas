suppressPackageStartupMessages(library(nlme))
out <- "Myeloid/disease_pathways/pathway_pseudotime"
d <- read.csv(file.path(out, "donor_bin_design.csv"))
d$status <- factor(d$status, levels = c("Healthy", "COPD"))
d$source <- factor(d$source)
d$donor <- factor(d$donor)
s <- as.matrix(read.csv(file.path(out, "donor_bin_pathway_scores.csv"), row.names = 1, check.names = FALSE))
stopifnot(identical(colnames(s), d$sample))
rows <- vector("list", nrow(s))
for (i in seq_len(nrow(s))) {
  d$score <- s[i, ]
  result <- tryCatch({
    fit <- lme(score ~ source + status*(b1+b2+b3), random = ~1|donor, data = d,
               method = "REML", control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
    shape <- grep("^statusCOPD:", names(fixef(fit)))
    beta <- fixef(fit)[shape]
    chisq <- as.numeric(crossprod(beta, solve(vcov(fit)[shape, shape], beta)))
    den <- min(fit$fixDF$X[shape])
    disease <- c(match("statusCOPD", names(fixef(fit))), shape)
    b <- fixef(fit)[disease]
    all_chisq <- as.numeric(crossprod(b, solve(vcov(fit)[disease, disease], b)))
    data.frame(pathway = rownames(s)[i], F = chisq/3,
               shape_p = pf(chisq/3, 3, den, lower.tail = FALSE),
               any_disease_p = pchisq(all_chisq, 4, lower.tail = FALSE),
               fit_status = "converged")
  }, error = function(e) data.frame(pathway = rownames(s)[i], F = NA,
                    shape_p = NA, any_disease_p = NA, fit_status = conditionMessage(e)))
  rows[[i]] <- result
  if (i %% 100 == 0) message("Completed pathway models ", i, "/", nrow(s))
}
results <- do.call(rbind, rows)
results$FDR_shape <- p.adjust(results$shape_p, "BH", n = nrow(s))
results$FDR_any_disease <- p.adjust(results$any_disease_p, "BH", n = nrow(s))
write.csv(results[order(results$shape_p), ], file.path(out, "all_pathways_individual_mixed_models.csv"), row.names = FALSE)
cat("Converged", sum(results$fit_status == "converged"), "of", nrow(results), "\n")
cat("Shape FDR discoveries", sum(results$FDR_shape < .05, na.rm = TRUE), "\n")
