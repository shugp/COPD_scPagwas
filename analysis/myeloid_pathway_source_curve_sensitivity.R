suppressPackageStartupMessages(library(nlme))
out <- "Myeloid/disease_pathways/pathway_pseudotime"
x <- readRDS(file.path(out, "focused_model_input.rds"))
d <- x$samples
all_scores <- as.matrix(read.csv(file.path(out, "donor_bin_pathway_scores.csv"), row.names = 1, check.names = FALSE))
screen <- read.csv(file.path(out, "all_pathways_individual_mixed_models.csv"))
screen_hits <- screen$pathway[screen$FDR_shape < .05]
paths <- unique(c(x$focused, screen_hits))
rows <- lapply(paths, function(path) {
  d$score <- as.numeric(all_scores[path, ])
  model <- lme(score ~ source*(b1+b2+b3) + status*(b1+b2+b3),
               random = ~1|donor, data = d, method = "REML",
               control = lmeControl(opt = "optim", msMaxIter = 200))
  names <- names(fixef(model))
  shape <- grep("^statusCOPD:|:statusCOPD$", names)
  beta <- fixef(model)[shape]
  Q <- as.numeric(crossprod(beta, solve(vcov(model)[shape, shape], beta)))
  any <- c(match("statusCOPD", names), shape)
  b <- fixef(model)[any]
  W <- as.numeric(crossprod(b, solve(vcov(model)[any, any], b)))
  data.frame(pathway = path, shape_p = pf(Q/3, 3, min(model$fixDF$X[shape]), lower.tail = FALSE),
             any_disease_p = pchisq(W, 4, lower.tail = FALSE))
})
results <- do.call(rbind, rows)
results$FDR_shape <- results$FDR_any_disease <- NA_real_
idx <- results$pathway %in% x$focused
results$FDR_shape[idx] <- p.adjust(results$shape_p[idx], "BH")
results$FDR_any_disease[idx] <- p.adjust(results$any_disease_p[idx], "BH")
results$FDR_shape_screen_sensitivity <- NA_real_
hits <- results$pathway %in% screen_hits
results$FDR_shape_screen_sensitivity[hits] <- p.adjust(results$shape_p[hits], "BH", n = nrow(screen))
write.csv(results, file.path(out, "source_curve_sensitivity.csv"), row.names = FALSE)
write.csv(results[results$pathway %in% x$focused, ],
          file.path(out, "antigen_source_curve_sensitivity.csv"), row.names = FALSE)
print(results, row.names = FALSE)
