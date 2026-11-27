# Equivalent LD entries, including last-write wins and repeated SNP names.
pagwas_sparse_ld <- function(all_snps, ld_data) {
  n <- length(all_snps)
  if (n == 1L) return(matrix(1))
  i <- match(as.character(ld_data[[1]]), all_snps)
  j <- match(as.character(ld_data[[2]]), all_snps)
  if (anyNA(i) || anyNA(j)) stop("LD contains SNPs outside the pathway.")
  values <- as.numeric(as.character(ld_data[[3]]))
  pairs <- data.table::data.table(i = pmin(i, j), j = pmax(i, j), value = values)
  pairs <- pairs[!duplicated(pairs, by = c("i", "j"), fromLast = TRUE)]
  diagonal <- rep(1, n)
  is_diagonal <- pairs$i == pairs$j
  diagonal[pairs$i[is_diagonal]] <- pairs$value[is_diagonal]
  pairs <- pairs[!is_diagonal]
  Matrix::sparseMatrix(i = c(seq_len(n), pairs$i, pairs$j),
                       j = c(seq_len(n), pairs$j, pairs$i),
                       x = c(diagonal, pairs$value, pairs$value),
                       dims = c(n, n), dimnames = list(all_snps, all_snps))
}

pagwas_factor_regression <- function(pa_block, pca_scCell_mat, data_mat,
                                     rawPathway_list, snp_gene_df,
                                     n.cores = 1, backingpath, Rns) {
  if (pa_block$n_snps <= 2L) {
    if (pa_block$n_snps == 0L) {
      pa_block$include_in_inference <- FALSE
      pa_block$x <- NULL
      return(pa_block)
    }
    return(NULL)
  }
  previous_blas <- bigparallelr::set_blas_ncores(min(n.cores, 8L))
  on.exit(bigparallelr::set_blas_ncores(previous_blas), add = TRUE)
  pathway <- unique(pa_block$block_info$pathway)
  genes <- intersect(rawPathway_list[[pathway]], rownames(data_mat))
  expression <- as.matrix(data_mat[genes, , drop = FALSE])
  if (anyNA(expression) || any(!is.finite(expression))) {
    stop("Non-finite pathway expression requires the original regression.")
  }
  if (length(genes) == 1L) {
    expression <- expression / (expression + 1e-4)
  } else {
    totals <- colSums(expression)
    nonzero <- totals != 0
    expression[, nonzero] <- sweep(expression[, nonzero, drop = FALSE],
                                   2, totals[nonzero], "/")
    expression[, !nonzero] <- 0
  }
  gene_index <- match(pa_block$snps$label, genes)
  if (anyNA(gene_index)) stop("Pathway SNP gene is absent from expression.")
  membership <- Matrix::sparseMatrix(i = seq_along(gene_index), j = gene_index,
                                     x = 1, dims = c(length(gene_index), length(genes)))
  # X = (LD^2 %*% membership) %*% expression * pathway_PCA.
  # Fit through the much smaller gene-level crossproducts.
  design <- as.matrix(pa_block$ld_matrix_squared %*% membership)
  noise <- pa_block$snps$se^2
  keep <- !is.na(pa_block$y) & !is.na(noise)
  design <- design[keep, , drop = FALSE]
  y <- pa_block$y[keep]
  covariates <- cbind(1, bigstatsr::covar_from_df(data.frame(offset(noise[keep]))))
  stopifnot(length(y) > ncol(covariates))
  decomposition <- svd(covariates, nv = 0)
  rank <- sum(decomposition$d / (sqrt(length(y)) + sqrt(ncol(covariates)) - 1) > 1e-4)
  basis <- decomposition$u[, seq_len(rank), drop = FALSE]
  residual <- design - basis %*% crossprod(basis, design)
  gram <- crossprod(residual)
  numerator <- as.numeric(crossprod(residual, y))
  score <- as.numeric(pca_scCell_mat[pathway, ])
  estimates <- numeric(ncol(expression))
  batch_size <- 4096L
  for (first in seq.int(1L, ncol(expression), by = batch_size)) {
    cells <- first:min(first + batch_size - 1L, ncol(expression))
    weights <- expression[, cells, drop = FALSE]
    denominator <- colSums(weights * (gram %*% weights))
    estimates[cells] <- as.numeric(crossprod(numerator, weights)) /
      denominator / score[cells]
  }
  estimates[score == 0 | is.na(estimates)] <- 0
  names(estimates) <- colnames(data_mat)
  estimates
}

# The installed as_matrix() truncates sparse expression to integers.
# Retain that conversion while keeping zero entries implicit.
pagwas_sparse_integer <- function(mat) {
  if (inherits(mat, "sparseMatrix")) {
    mat <- methods::as(mat, "dgCMatrix")
    mat@x <- as.numeric(as.integer(mat@x))
    return(Matrix::drop0(mat))
  }
  mat
}

pagwas_fast_pca <- function(original, original_pca) {
  pca_env <- new.env(parent = environment(original_pca))
  pca_env$colMeans <- Matrix::colMeans
  pca_env$t <- Matrix::t
  environment(original_pca) <- pca_env
  wrapper_env <- new.env(parent = environment(original))
  wrapper_env$as_matrix <- pagwas_sparse_integer
  wrapper_env$PathwayPCAtest <- original_pca
  environment(original) <- wrapper_env
  original
}

pagwas_sparse_pcc <- function(scPagwas.gPAS.score, data_mat) {
  mat <- pagwas_sparse_integer(data_mat)
  n <- ncol(mat)
  means <- Matrix::rowMeans(mat)
  y <- as.numeric(scPagwas.gPAS.score)
  covariance <- as.numeric(mat %*% y) - n * means * mean(y)
  scale_x <- sqrt(Matrix::rowSums(mat^2) - n * means^2)
  scale_y <- sqrt(sum(y^2) - n * mean(y)^2)
  result <- matrix(covariance / (scale_x * scale_y), ncol = 1,
                    dimnames = list(rownames(mat), "PCC"))
  result[is.nan(result)] <- 0
  result
}

pagwas_save_checkpoint <- function(object, path) {
  temporary <- tempfile("checkpoint-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) stop("Could not commit checkpoint ", path)
}

install_pagwas_fast <- function(checkpoint_dir) {
  ns <- asNamespace("scPagwas")
  originals <- mget(c("make_ld_matrix", "get_Pathway_sclm", "Pathway_pcascore_run",
                      "scGet_PCC", "Get_CorrectBg_p", "Pathway_block_func",
                      "scPagwas_perform_score", "link_pwpca_block"), ns)
  fast_pca <- pagwas_fast_pca(originals$Pathway_pcascore_run,
                             get("PathwayPCAtest", ns))
  cached_pca <- function(...) {
    args <- list(...)
    cache_file <- file.path(checkpoint_dir, "pathway_pca.rds")
    if (file.exists(cache_file)) {
      cache <- readRDS(cache_file)
      p <- args$Pagwas
      p[names(cache$result)] <- cache$result
      assign(".Random.seed", cache$rng, envir = .GlobalEnv)
      message("Restored pathway PCA checkpoint")
      return(p)
    }
    result <- fast_pca(...)
    pagwas_save_checkpoint(list(result = result[setdiff(names(result), "data_mat")],
                                rng = get(".Random.seed", .GlobalEnv)), cache_file)
    result
  }
  pathway_number <- 0L
  cached_regression <- function(pa_block, pca_scCell_mat, data_mat, rawPathway_list,
                                snp_gene_df, n.cores, backingpath, Rns) {
    pathway_number <<- pathway_number + 1L
    pathway <- unique(pa_block$block_info$pathway)
    cache_file <- file.path(checkpoint_dir, paste0(pathway, "_coefficients.rds"))
    started <- Sys.time()
    if (file.exists(cache_file)) {
      result <- readRDS(cache_file)
    } else {
      result <- pagwas_factor_regression(pa_block, pca_scCell_mat, data_mat,
        rawPathway_list, snp_gene_df, n.cores, backingpath, Rns)
      # Check selected real pathways against the installed estimator.
      if (pa_block$n_snps > 2 && (pathway_number <= 3 || pathway_number %% 50 == 0)) {
        cells <- unique(as.integer(seq(1, ncol(data_mat), length.out = 32)))
        reference <- originals$get_Pathway_sclm(pa_block,
          pca_scCell_mat[, cells, drop = FALSE], data_mat[, cells, drop = FALSE],
          rawPathway_list, snp_gene_df, 1L, backingpath, paste0(Rns, "_check"))
        if (!isTRUE(all.equal(reference, result[cells], tolerance = 1e-7))) {
          stop("Regression validation failed for ", pathway)
        }
        message("Validated installed regression for ", pathway)
      }
      pagwas_save_checkpoint(result, cache_file)
    }
    message(format(Sys.time(), "%F %T"), " pathway ", pathway_number,
            " ", pathway, ", SNPs=", pa_block$n_snps,
            ", regression_seconds=", round(as.numeric(difftime(Sys.time(), started, units="secs")), 2))
    result
  }
  # The package uses diag() only for the empty-LD case in this function.
  block_function <- originals$Pathway_block_func
  block_env <- new.env(parent = environment(block_function))
  block_env$diag <- function(x = 1, nrow) Matrix::Diagonal(nrow, x = x)
  environment(block_function) <- block_env

  link_function <- originals$link_pwpca_block
  link_env <- new.env(parent = environment(link_function))
  link_env$t <- Matrix::t
  environment(link_function) <- link_env
  sparse_link <- function(...) {
    result <- link_function(...)
    if (!is.null(result$x)) result$x <- as.matrix(result$x)
    result
  }

  background_function <- originals$Get_CorrectBg_p
  background_env <- new.env(parent = environment(background_function))
  background_env$apply <- function(X, MARGIN, FUN, ...) {
    if (inherits(X, "sparseMatrix") && identical(MARGIN, 1) && identical(FUN, stats::var)) {
      n <- ncol(X)
      sums <- Matrix::rowSums(X)
      return(pmax(0, (Matrix::rowSums(X^2) - sums^2 / n) / (n - 1)))
    }
    base::apply(X, MARGIN, FUN, ...)
  }
  environment(background_function) <- background_env
  replacements <- list(make_ld_matrix = pagwas_sparse_ld,
                       get_Pathway_sclm = cached_regression,
                       Pathway_pcascore_run = cached_pca,
                       scGet_PCC = pagwas_sparse_pcc,
                       scPagwas_perform_score = pagwas_fast_score,
                       Get_CorrectBg_p = background_function,
                       link_pwpca_block = sparse_link,
                       Pathway_block_func = block_function)
  for (name in names(replacements)) {
    utils::assignInNamespace(name, replacements[[name]], ns = "scPagwas")
  }
  function() {
    for (name in names(originals)) {
      utils::assignInNamespace(name, originals[[name]], ns = "scPagwas")
    }
  }
}

pagwas_fast_score <- function(Pagwas, remove_outlier = TRUE) {
  pathway_names <- colnames(Pagwas$Pathway_sclm_results)
  pathway_expression <- lapply(pathway_names, function(pathway) {
    genes <- intersect(Pagwas$Pathway_list[[pathway]], rownames(Pagwas$data_mat))
    if (!length(genes)) return(NULL)
    Matrix::colMeans(Pagwas$data_mat[genes, , drop = FALSE])
  })
  present <- !vapply(pathway_expression, is.null, logical(1))
  pathway_expression <- data.matrix(as.data.frame(pathway_expression[present]))
  pathway_names <- pathway_names[present]
  colnames(pathway_expression) <- pathway_names
  weighted_expression <- t(Pagwas$pca_scCell_mat[pathway_names, , drop = FALSE]) *
    pathway_expression
  Pagwas$Pathway_single_results <- t(Pagwas$Pathway_sclm_results[, pathway_names, drop = FALSE] *
                                      weighted_expression)
  colnames(Pagwas$Pathway_single_results) <- colnames(Pagwas$pca_scCell_mat)
  message("* Get Pathways'rankPvalue for each celltypes!")
  celltypes <- unique(Pagwas$Celltype_anno$annotation)
  rank_p <- lapply(celltypes, function(celltype) {
    print(celltype)
    cells <- Pagwas$Celltype_anno$annotation == celltype
    getFromNamespace("scGene_rankP", "scPagwas")(
      Pagwas$Pathway_single_results[, cells, drop = FALSE])$pValueHigh
  })
  Pagwas$scPathways_rankPvalue <- as.data.frame(do.call(cbind, rank_p))
  colnames(Pagwas$scPathways_rankPvalue) <- celltypes
  rownames(Pagwas$scPathways_rankPvalue) <- rownames(Pagwas$Pathway_single_results)
  message("* Get scPgwas score for each single cell")
  scores <- colSums(Pagwas$Pathway_single_results)
  gc()
  if (remove_outlier) {
    scores <- getFromNamespace("scPagwas_score_filter", "scPagwas")(scores)
  }
  Pagwas$scPagwas.gPAS.score <- setNames(scores, colnames(Pagwas$pca_scCell_mat))
  Pagwas
}
pagwas_module_context <- function(mat, seed=1, nbin=24) {
  set.seed(seed)
  avg <- Matrix::rowMeans(mat)
  avg <- avg[order(avg)]
  cut <- ggplot2::cut_number(avg + rnorm(length(avg))/1e30, n=nbin,
                             labels=FALSE, right=FALSE)
  names(cut) <- names(avg)
  list(cut=cut,pools=split(cut,cut),rng=get(".Random.seed",.GlobalEnv),
       genes=rownames(mat))
}

pagwas_module_weights <- function(context, features, ctrl=100) {
  weights <- matrix(0,length(context$genes),length(features))
  for (i in seq_along(features)) {
    genes <- intersect(features[[i]], context$genes)
    controls <- unique(unlist(lapply(genes,function(gene) {
      names(sample(context$pools[[as.character(context$cut[gene])]],
                   size=ctrl,replace=FALSE))
    }),use.names=FALSE))
    selected <- match(genes,context$genes)
    background <- match(controls,context$genes)
    weights[selected,i] <- 1/length(selected)
    weights[background,i] <- weights[background,i] - 1/length(background)
  }
  weights
}

pagwas_module_contrast <- function(mat,features,context=NULL) {
  if (is.null(context)) context <- pagwas_module_context(mat)
  assign(".Random.seed",context$rng,envir=.GlobalEnv)
  weights <- pagwas_module_weights(context,features)
  t(as.matrix(Matrix::crossprod(weights,mat)))
}

pagwas_pathway_means <- function(data_mat, pathways) {
  genes <- lapply(pathways,intersect,y=rownames(data_mat))
  sizes <- lengths(genes)
  if (any(sizes == 0L)) stop("Empty pathway in expression score.")
  membership <- Matrix::sparseMatrix(
    i=rep(seq_along(genes),sizes),
    j=match(unlist(genes,use.names=FALSE),rownames(data_mat)),
    x=rep(1/sizes,sizes),
    dims=c(length(genes),nrow(data_mat)),
    dimnames=list(names(pathways),rownames(data_mat)))
  result <- t(as.matrix(membership %*% data_mat))
  dimnames(result) <- list(colnames(data_mat),names(pathways))
  result
}

# Reuse the package's score assembly and ranking after replacing only the means.
pagwas_aggregate_score <- function(Pagwas, remove_outlier=TRUE) {
  score_function <- pagwas_fast_score
  env <- new.env(parent=environment(score_function))
  env$lapply <- function(X,FUN,...) {
    if (identical(X,colnames(Pagwas$Pathway_sclm_results))) {
      means <- pagwas_pathway_means(Pagwas$data_mat,Pagwas$Pathway_list[X])
      return(lapply(seq_along(X),function(i) means[,i]))
    }
    base::lapply(X,FUN,...)
  }
  environment(score_function) <- env
  score_function(Pagwas,remove_outlier)
}

pagwas_indexed_blocks <- function(original) {
  rewrite_queries <- function(expression) {
    if (!is.call(expression)) return(expression)
    if (identical(expression[[1]],as.name("<-")) && is.symbol(expression[[2]])) {
      target <- as.character(expression[[2]])
      if (target == "rsids") {
        expression[[3]] <- quote(indexed_snps(chrom_block$label))
        return(expression)
      }
      if (target == "rsids_gwas") {
        expression[[3]] <- quote(indexed_join(rsids,chrom))
        return(expression)
      }
    }
    for (i in seq_along(expression)[-1]) {
      expression[i] <- list(rewrite_queries(expression[[i]]))
    }
    expression
  }
  body(original) <- rewrite_queries(body(original))
  function(Pagwas=NULL,Pachrom_block_list,chrom_gwas_list,singlecell=TRUE,
           celltype=TRUE,chrom_ld,backingpath,n.cores=1) {
    gene_rows <- split(seq_len(nrow(Pagwas$snp_gene_df)),Pagwas$snp_gene_df$label)
    gwas_index <- lapply(chrom_gwas_list,function(gwas) {
      gwas <- data.table::copy(gwas)
      data.table::setkeyv(gwas,"rsid")
      gwas
    })
    env <- new.env(parent=environment(original))
    env$indexed_snps <- function(labels) {
      rows <- sort(unlist(gene_rows[labels],use.names=FALSE))
      Pagwas$snp_gene_df[rows,c("rsid","label"),drop=FALSE]
    }
    env$indexed_join <- function(rsids,chrom) {
      gwas <- gwas_index[[chrom]]
      joined <- gwas[data.table::as.data.table(rsids),on="rsid",nomatch=0L]
      columns <- c("rsid","label",setdiff(names(gwas),"rsid"))
      as.data.frame(joined[,columns,with=FALSE])
    }
    environment(original) <- env
    original(Pagwas=Pagwas,Pachrom_block_list=Pachrom_block_list,
             chrom_gwas_list=chrom_gwas_list,singlecell=singlecell,celltype=celltype,
             chrom_ld=chrom_ld,backingpath=backingpath,n.cores=n.cores)
  }
}

install_factorized_pagwas <- install_pagwas_fast
install_pagwas_fast <- function(checkpoint_dir,bootstrap_cache=NULL,
                                prefiltered_gwas=FALSE) {
  restore_factorized <- install_factorized_pagwas(checkpoint_dir)
  ns <- asNamespace("scPagwas")
  originals <- mget(c("Pathway_block_func","scPagwas_perform_score",
                      "Get_CorrectBg_p","Boot_evaluate","GWAS_summary_input",
                      "get_p_from_empi_null"),ns)
  original_module <- getFromNamespace("AddModuleScore","Seurat")
  module_context <- NULL
  module_matrix <- NULL
  module_context_for <- function(object,assay) {
    if (is.null(module_context)) {
      module_matrix <<- Seurat::GetAssayData(object,assay=assay,layer="data")
      module_context <<- pagwas_module_context(module_matrix)
    }
    module_context
  }
  optimized_module <- function(object,features,pool=NULL,nbin=24,ctrl=100,
    k=FALSE,assay=NULL,name="Cluster",seed=1,search=FALSE,slot="data",...) {
    assay <- if (is.null(assay)) Seurat::DefaultAssay(object) else assay
    if (!is.null(pool) || nbin!=24 || ctrl!=100 || k || search ||
        !identical(seed,1) || slot!="data") {
      return(original_module(object,features,pool,nbin,ctrl,k,assay,name,seed,search,slot,...))
    }
    context <- module_context_for(object,assay)
    scores <- as.data.frame(pagwas_module_contrast(module_matrix,features,context))
    colnames(scores) <- paste0(name,seq_along(features))
    rownames(scores) <- colnames(object)
    object[[colnames(scores)]] <- scores
    object
  }
  optimized_background <- function(Single_data,scPagwas.TRS.Score,iters_singlecell,
                                   n_topgenes,scPagwas_topgenes,assay="RNA") {
    context <- module_context_for(Single_data,assay)
    weights <- matrix(0,nrow(module_matrix),iters_singlecell)
    random_genes <- vector("list",iters_singlecell)
    for(i in seq_len(iters_singlecell)) {
      set.seed(i)
      random_genes[[i]] <- sample(rownames(Single_data),n_topgenes)
      assign(".Random.seed",context$rng,envir=.GlobalEnv)
      weights[,i] <- pagwas_module_weights(context,random_genes[i])
    }
    message(Sys.time()," background weights ready for ",iters_singlecell," iterations")
    controls <- t(as.matrix(Matrix::crossprod(weights,module_matrix)))
    message(Sys.time()," background module scores ready")
    genes <- intersect(rownames(Single_data),rownames(module_matrix))
    gene_matrix <- module_matrix[genes,,drop=FALSE]
    n <- ncol(gene_matrix)
    gene_variance <- pmax(0,(Matrix::rowSums(gene_matrix^2) -
                              Matrix::rowSums(gene_matrix)^2/n)/(n-1))
    names(gene_variance) <- genes
    top <- intersect(scPagwas_topgenes,genes)
    ratio <- vapply(random_genes,function(x) sum(gene_variance[x]),numeric(1)) /
      sum(gene_variance[top])
    result <- getFromNamespace("correct_background","scPagwas")(
      scPagwas.TRS.Score,controls,ratio)
    rownames(result) <- colnames(Single_data)
    result
  }
  cached_bootstrap <- function(Pagwas,bootstrap_iters=200,part=0.5) {
    if (is.null(bootstrap_cache) || !file.exists(bootstrap_cache) ||
        bootstrap_iters != 200 || part != 0.5) {
      return(originals$Boot_evaluate(Pagwas,bootstrap_iters,part))
    }
    result <- read.csv(bootstrap_cache,row.names=1,check.names=FALSE)
    # The table contains FDR-adjusted values. Restore raw values so the
    # original main function applies its usual correction exactly once.
    expected <- (result$bias_corrected_estimate + result$bootstrap_estimate)/2
    observed <- Pagwas$lm_results$parameters[result$annotation]
    if (!isTRUE(all.equal(unname(observed),expected,tolerance=1e-8))) {
      stop("Cached bootstrap model differs from current data.")
    }
    result$bp_value <- stats::pnorm(result$bt_value,lower.tail=FALSE)
    Pagwas$bootstrap_results <- result
    message("Restored completed celltype bootstrap results")
    Pagwas
  }
  passthrough_gwas <- function(Pagwas,gwas_data,...) {
    Pagwas$gwas_data <- gwas_data
    Pagwas
  }
  replacements <- list(
    Pathway_block_func=pagwas_indexed_blocks(originals$Pathway_block_func),
    scPagwas_perform_score=pagwas_aggregate_score,
    Get_CorrectBg_p=optimized_background,
    Boot_evaluate=cached_bootstrap,
    get_p_from_empi_null=pagwas_exact_null_p)
  if (prefiltered_gwas) replacements$GWAS_summary_input <- passthrough_gwas
  for(name in names(replacements)) {
    utils::assignInNamespace(name,replacements[[name]],ns="scPagwas")
  }
  utils::assignInNamespace("AddModuleScore",optimized_module,ns="Seurat")
  function() {
    utils::assignInNamespace("AddModuleScore",original_module,ns="Seurat")
    for(name in names(originals)) {
      utils::assignInNamespace(name,originals[[name]],ns="scPagwas")
    }
    restore_factorized()
  }
}

pagwas_exact_null_p <- function(v_t,v_t_null) {
  n <- length(v_t_null)
  if (!n) {
    result <- rep(1,length(v_t))
  } else if (anyNA(v_t_null)) {
    result <- rep(NA_real_,length(v_t))
  } else {
    positions <- findInterval(v_t,sort(v_t_null))
    result <- (n-positions+1)/(n+1)
  }
  names(result) <- names(v_t)
  result
}
