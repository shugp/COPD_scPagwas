suppressPackageStartupMessages(library(nlme))
out <- "Myeloid/disease_pathways/pathway_pseudotime"
x <- readRDS(file.path(out, "focused_model_input.rds"))
d <- x$samples
scores <- as.matrix(read.csv(file.path(out, "donor_bin_pathway_scores.csv"), row.names = 1, check.names = FALSE))
screen <- read.csv(file.path(out, "all_pathways_individual_mixed_models.csv"))
paths <- head(screen$pathway, 4)
grid <- seq(max(tapply(d$pseudotime, d$status, min)), min(tapply(d$pseudotime, d$status, max)), length.out=150)
basis <- splines::ns(grid, knots=x$knots, Boundary.knots=x$boundary)
rows <- list()
for(path in paths) {
  d$score <- scores[path,]
  fit <- lme(score ~ source + status*(b1+b2+b3), random=~1|donor, data=d, method="REML",
             control=lmeControl(opt="optim",msMaxIter=200))
  for(status in levels(d$status)) {
    matrices <- lapply(levels(d$source),function(source) {
      new <- data.frame(source=factor(source,levels=levels(d$source)),
                        status=factor(status,levels=levels(d$status)),
                        b1=basis[,1],b2=basis[,2],b3=basis[,3])
      model.matrix(~source+status*(b1+b2+b3),new)
    })
    X <- Reduce("+",matrices)/length(matrices)
    pred <- as.numeric(X%*%fixef(fit))
    se <- sqrt(rowSums((X%*%vcov(fit))*X))
    rows[[paste(path,status)]] <- data.frame(pathway=path,status=status,pseudotime=grid,
      fitted_score=pred,ci_low=pred-1.96*se,ci_high=pred+1.96*se)
  }
}
write.csv(do.call(rbind,rows),file.path(out,"exploratory_pathway_fitted_curves.csv"),row.names=FALSE)
