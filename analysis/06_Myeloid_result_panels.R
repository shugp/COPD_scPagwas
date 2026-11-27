# Export directly comparable residual correlations for the antigen-presentation panel.
x <- readRDS("Myeloid/disease_pathways/donor_state_expression.rds")
candidates <- c("HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRB1", "B2M", "ERAP2")
partners <- c("CIITA", "CD74", "HLA-DMA", "HLA-DMB", "CTSS", "IFI30",
              "TAP1", "TAP2", "TAPBP", "HLA-A", "HLA-B", "HLA-C",
              "NLRC5", "PSMB8", "PSMB9", "STAT1", "IRF1")
genes <- intersect(c(candidates, partners), rownames(x$E))
d <- x$samples
design <- model.matrix(~ factor(donor) + factor(State), d)
residuals <- qr.resid(qr(design), t(x$E[genes, , drop = FALSE]))
correlations <- cor(residuals)
df <- nrow(residuals) - qr(design)$rank - 1
rows <- expand.grid(candidate = intersect(candidates, genes),
                    partner = intersect(partners, genes))
rows$partial_correlation <- correlations[cbind(match(rows$candidate, genes), match(rows$partner, genes))]
r <- pmin(pmax(rows$partial_correlation, -.999999), .999999)
rows$p <- 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)
rows$FDR <- p.adjust(rows$p, "BH")
write.csv(rows, "Myeloid/disease_pathways/antigen_residual_coexpression.csv", row.names = FALSE)
