packages <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "ggrepel", "limma", "statmod", "nlme", "cluster")
if ("--scpagwas" %in% commandArgs(TRUE)) {
  packages <- c(packages, "scPagwas", "bigparallelr", "bigreadr", "bigstatsr")
}
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
versions <- vapply(packages, function(p) as.character(packageVersion(p)), character(1))
print(data.frame(package = packages, version = versions), row.names = FALSE)
