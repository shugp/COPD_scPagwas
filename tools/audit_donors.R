# Read the actual input object's metadata; write no additional manuscript table.
suppressPackageStartupMessages(library(SeuratObject))
dat <- readRDS("Myeloid/Myeloid.Rds")
meta <- dat[[]]
samples <- unique(meta[, c("data_source", "orig.ident", "status")])
stopifnot(!anyDuplicated(samples[, c("data_source", "orig.ident")]))
samples$donor <- as.character(samples$orig.ident)
repeated <- samples$data_source == "GSE168191"
samples$donor[repeated] <- sub("_rep[0-9]+\\.h5$", "", samples$donor[repeated])
samples$donor <- paste(samples$data_source, samples$donor, sep = ":")
donors <- unique(samples[, c("data_source", "donor", "status")])
stopifnot(!anyDuplicated(donors$donor))
cat("Input cells:", nrow(meta), "\norig.ident labels:", nrow(samples),
    "\nBiological donor units:", nrow(donors), "\n\nSample labels by status:\n")
print(table(samples$status))
cat("\nDonor units by status:\n")
print(table(donors$status))
cat("\nDonor units by source and status:\n")
print(table(donors$data_source, donors$status))
cat("\nSample-to-donor mapping (GSE168191 rep1/rep2 merged within donor):\n")
print(samples[order(samples$data_source, samples$orig.ident), ], row.names = FALSE)
