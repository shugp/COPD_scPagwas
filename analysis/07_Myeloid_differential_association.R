# Export aligned independent donor-State observations for differential association.
x <- readRDS("Myeloid/disease_pathways/donor_state_expression.rds")
candidate <- read.csv("Myeloid/results_cleaned/coloc_overall_gene_tests.csv")$feature
scores <- read.csv("Myeloid/disease_pathways/antigen_donor_state_scores.csv")
samples <- x$samples
stopifnot(identical(samples$sample, colnames(x$E)),
          !anyDuplicated(paste(samples$donor, samples$State)))
samples$data_source <- sub(":.*$", "", samples$donor)
expression <- as.data.frame(t(x$E[intersect(candidate, rownames(x$E)), , drop = FALSE]))
d <- cbind(samples[, c("sample", "donor", "State", "status", "data_source")], expression)
stopifnot(!anyDuplicated(paste(scores$sample, scores$pathway)))
write.csv(d, "Myeloid/disease_pathways/differential_association_input.csv", row.names = FALSE)
status <- system2("python", shQuote(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_differential_association.py")),
                  stdout = "Myeloid/disease_pathways/differential_association.log",
                  stderr = "Myeloid/disease_pathways/differential_association.log")
stopifnot(status == 0L)

status <- system2("python", shQuote(file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "myeloid_differential_association_sensitivity.py")),
                  stdout = "Myeloid/disease_pathways/differential_association_sensitivity.log",
                  stderr = "Myeloid/disease_pathways/differential_association_sensitivity.log")
stopifnot(status == 0L)
