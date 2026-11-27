project_dir <- normalizePath(getwd())
setwd(project_dir)
library(Seurat)
library(scPagwas)
helper_path <- file.path(Sys.getenv("MYELOID_CODE_DIR", "analysis"), "scPagwas_fast_helpers.R")
source(helper_path)

run_scPagwas <- function() {
  input_files <- file.path(project_dir, c("inte_all.rds", "copdgwas.txt"))
  input_info <- file.info(input_files)
  if (anyNA(input_info$size)) stop("Single-cell or GWAS input file is missing.")
  cores <- suppressWarnings(as.integer(Sys.getenv("SCPAGWAS_CORES", "16")))
  if (is.na(cores) || cores < 1L) stop("SCPAGWAS_CORES must be a positive integer.")
  cores <- min(cores, bigparallelr::nb_cores(), getOption("bigstatsr.ncores.max", Inf))
  output_dir <- file.path(project_dir, "pagwas")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- file.path(project_dir, "pagwas_accelerated", "checkpoints")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  configuration <- list(
    algorithm = "factorized-regression-weighted-modules-exact-null-v3",
    implementation = unname(tools::md5sum(helper_path)),
    scPagwas = as.character(packageVersion("scPagwas")),
    inputs = input_info[, c("size", "mtime"), drop = FALSE], assay = "RNA",
    singlecell = TRUE, celltype = TRUE, iters_singlecell = 100, seed = 1L,
    pathways = Genes_by_pathway_kegg
  )
  manifest <- file.path(checkpoint_dir, "configuration.rds")
  if (file.exists(manifest)) {
    if (!identical(readRDS(manifest), configuration)) {
      stop("Checkpoint configuration differs from current inputs or parameters.")
    }
  } else {
    pagwas_save_checkpoint(configuration, manifest)
  }

  thread_vars <- c("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                   "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "R_LOCAL_CACHE")
  old_env <- Sys.getenv(thread_vars, unset = NA_character_)
  on.exit({
    Sys.unsetenv(names(old_env)[is.na(old_env)])
    if (any(!is.na(old_env))) do.call(Sys.setenv, as.list(old_env[!is.na(old_env)]))
  }, add = TRUE)
  do.call(Sys.setenv, setNames(as.list(rep("1", 4)), thread_vars[1:4]))
  old_blas <- bigparallelr::set_blas_ncores(1)
  on.exit(bigparallelr::set_blas_ncores(old_blas), add = TRUE)
  old_options <- options(default.nproc.blas = 1L, bigreadr.nThread = cores,
                         save.defaults = modifyList(getOption("save.defaults", list()),
                                                    list(compress = FALSE)))
  on.exit(options(old_options), add = TRUE)
  gwas_cache <- file.path(checkpoint_dir, "gwas_filtered.rds")
  if (!file.exists(gwas_cache)) {
    raw_gwas <- bigreadr::fread2(input_files[2])
    filtered <- getFromNamespace("GWAS_summary_input", "scPagwas")(
      Pagwas = list(), gwas_data = raw_gwas, maf_filter = 0.01)$gwas_data
    pagwas_save_checkpoint(filtered, gwas_cache)
    rm(raw_gwas, filtered)
  }
  filtered_gwas <- readRDS(gwas_cache)
  restore_functions <- install_pagwas_fast(
    checkpoint_dir,
    bootstrap_cache = file.path(checkpoint_dir, "celltype_bootstrap_results.csv"),
    prefiltered_gwas = TRUE)
  on.exit(restore_functions(), add = TRUE)

  scratch_dir <- tempfile("scPagwas-fast-", tmpdir = Sys.getenv("SCPAGWAS_SCRATCH", "/tmp"))
  if (!dir.create(scratch_dir, mode = "0700")) stop("Cannot create scratch directory.")
  old_wd <- setwd(scratch_dir)
  on.exit({
    setwd(old_wd)
    if (dir.exists(file.path(scratch_dir, "pagwas"))) {
      tables <- list.files(file.path(scratch_dir, "pagwas"),
                           pattern = "^inte_all_.*\\.(csv|txt)$", full.names = TRUE)
      if (length(tables) && !all(file.copy(tables, output_dir, overwrite = TRUE))) {
        warning("Some analysis tables remain in ", scratch_dir)
      }
    }
  }, add = TRUE)
  set.seed(1)
  message(Sys.time(), " accelerated scPagwas, cores=", cores, ", scratch=", scratch_dir)
  result <- scPagwas_main(
    Pagwas = NULL, gwas_data = filtered_gwas, Single_data = input_files[1],
    output.prefix = "inte_all", output.dirs = "pagwas",
    Pathway_list = Genes_by_pathway_kegg, n.cores = cores, assay = "RNA",
    singlecell = TRUE, iters_singlecell = 100, celltype = TRUE,
    block_annotation = block_annotation, chrom_ld = chrom_ld
  )
  # Save before plotting, using the project's durable storage.
  Pagwas <- result
  save(Pagwas, file = file.path(project_dir, "scPagwas_all.RData"),
       compress = FALSE)
  result
}

Pagwas <- run_scPagwas()
scPagwas_Visualization(
  Single_data = Pagwas, p_thre = 0.05, FigureType = "umap",
  width = 7, height = 7, lowColor = "#000957", highColor = "#EBE645",
  output.dirs = "gPAS", size = 0.5, do_plot = FALSE
)
plot_bar_positie_nagtive(
  seurat_obj = Pagwas, var_ident = "cell_type", var_group = "positiveCells",
  vec_group_colors = c("#E8D0B3", "#7EB5A6"), do_plot = FALSE,
  output.prefix = "inte_all", output.dirs = "gPAS"
)
message(Sys.time(), " completed accelerated scPagwas")
