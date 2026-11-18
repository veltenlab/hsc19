#'Decompose RNA-seq data into 19 gene regulatory programs
#'
#'Ideally run this function on RNA-seq data of HSCs, e.g. bulk RNA-seq data of CD34+ cells (human) or LSK cells (mouse), or single cell data subset to HSCs. This function can run on bulk, pseudobulk or single cell RNA-seq.
#'
#'@details The species and gene ID format of the data (human or mouse; symbols or ENSEMBL IDs) is detected from the row names, and the row names of \code{loadings} (mouse gene symbols for \code{hsc19_loadings}) are converted accordingly with \code{\link{convertIDs}}. For each sample, the log-normalized expression of the genes shared between data and loadings is regressed on the loadings with a linear model; the coefficients are the inferred program activities. If \code{counts} is provided, a Seurat object is created and normalized internally. 
#'@param counts Count matrix, with gene IDs as row names and sample identifiers as column names. Supported gene name formats are mouse or human gene symbols (preferred), or ENSEMBL (ENSG/ENSMUSG) IDs.
#'@param s A seurat object, with \code{NormalizeData} and \code{FindVariableFeatures} run. One of \code{counts} or \code{s} needs to be provided. If \code{counts} is provided, \code{s} is ignored. The normalized data of the default assay are used.
#'@param loadings Loading matrix with genes (mouse gene symbols, as in \code{hsc19_loadings}) in rows and factors/programs in columns. Row and column names must be present and unique. Defaults to \code{hsc19_loadings}.
#'@return Returns an object of class \code{\link{decomposition}}, with the inferred activity of each program in each sample in the \code{result} slot.
#'@export
decompose <- function(counts = NULL, s = NULL, loadings = hsc19_loadings) {
  if (!is.null(counts)) {
    s <- Seurat::CreateSeuratObject(counts)
    s <- Seurat::NormalizeData(s)
    s <- Seurat::FindVariableFeatures(s)
    s <- Seurat::ScaleData(s)
  } else if (is.null(s)) {
    stop("Either a count matrix or a seurat object needs to be supplied")
  } else {
    if (is.null(VariableFeatures(s))) stop("If providing a seurat assay, data needs to be normalized and FindVariableFeatures needs to have run.")
  }

  if (!is.matrix(loadings)) {
    stop("'loadings' must be a matrix.")
  }
  if (is.null(rownames(loadings)) || is.null(colnames(loadings))) {
    stop("'loadings' must have row names (genes) and column names (factor/program names).")
  }
  if (anyDuplicated(rownames(loadings))) {
    stop("'loadings' must have unique row names.")
  }
  if (anyDuplicated(colnames(loadings))) {
    stop("'loadings' must have unique column names.")
  }

  #determine format of rownames
  if (mean(grepl("^ENSG", rownames(s))) > 0.9) {
    rid <- "hens"
  } else if (mean(grepl("^ENSMUSG", rownames(s))) > 0.9) {
    rid <- "mens"
  } else if (mean(grepl("^[^a-z]+$", rownames(s))) > 0.9) {
    rid <- "hsy"
  } else rid <- "msy"

  use_W_mat <- loadings
  #convert gene ids
  use_W_mat <- convertIDs(use_W_mat, rid)
  factor_names <- colnames(use_W_mat)

  res <- list()
  exprs <- as.matrix(Seurat::GetAssayData(s))
  usegenes = intersect(rownames(use_W_mat), rownames(exprs))
  for (i in 1:ncol(exprs)) {
    modelf <- data.frame(response = exprs[usegenes, i], use_W_mat[usegenes, , drop = FALSE], check.names = FALSE)
    mo <- lm(response ~ ., data = modelf)
    coef_i <- coef(mo)
    names(coef_i) <- c("(Intercept)", factor_names)
    res[[i]] <- coef_i
    attr(res[[i]], "R2") <- summary(mo)$r.squared
  }
  resmat <- do.call(cbind,res)
  colnames(resmat) <- colnames(exprs)
  n_factors <- ncol(use_W_mat)

  #percent variance explained
  #denominator: sum of variance of all genes
  pve_denom <- sum(apply(exprs[usegenes,],1,var))
  pve_numer <- sum(apply(use_W_mat[usegenes,], 1, function(x) var(t(resmat) %*% c(1,x))))
  pca <- prcomp(t(exprs[usegenes,]))
  #pca_hvg <- prcomp(t(exprs[VariableFeatures(s),]))
  if (ncol(s) > 50) s <- Seurat::RunPCA(s, verbose=F, npcs = min(c(50, ncol(s))))

  out <- new("decomposition")
  out@result <- resmat
  out@factor_names <- factor_names
  out@pve <- pve_numer / pve_denom
  if (ncol(s) > 50) out@pca <- s@reductions$pca
  out@s <- s
  out@pve_pca <- cumsum(pca$sdev^2 / sum(pca$sdev^2))[min(n_factors, length(pca$sdev))]
  out@usegenes <- usegenes

  cat(sprintf("Run on %d samples. On the overlapping gene set of %d genes, %d programs explain %.1f %% of the gene expression variance in your data. For comparison, %d principal components explain %.1f %%\n",
              ncol(resmat), length(out@usegenes), n_factors, 100 * out@pve, n_factors, 100 * out@pve_pca))

  return(out)
}
#' Convert gene IDs
#'
#' Convert mouse gene symbols to different formats
#'
#' @param x A matrix with (mouse) gene symbols as row names, or a character vector of (mouse) gene symbols.
#' @param rid One of hsy (human gene symbol, obtained by upper-casing), hens (human ENSEMBL IDs; symbols are upper-cased and mapped via org.Hs.eg.db) or mens (mouse ENSEMBL IDs, mapped via org.Mm.eg.db). Any other value (e.g. msy) returns \code{x} unchanged.
#'
#' @return For matrices, the matrix with converted row names; for hens and mens, genes that cannot be mapped or map ambiguously to more than one ENSEMBL ID are removed. For character vectors, the converted IDs; for hens and mens, unmapped symbols are dropped (so the result can differ in length and order from \code{x}) and symbols with several ENSEMBL IDs yield several entries.
#' @export
setGeneric("convertIDs", function(x, rid) {
  standardGeneric("convertIDs")
})

#' @describeIn convertIDs Method for matrices
#' @exportMethod convertIDs
#' @aliases convertIDs,matrix,character-method
setMethod("convertIDs", signature(x = "matrix", rid = "character"),
          function(x, rid) {
            convertIDs_matrix(x, rid)
          })

#' @describeIn convertIDs Method for character vectors
#' @exportMethod convertIDs
#' @aliases convertIDs,character,character-method
setMethod("convertIDs", signature(x = "character", rid = "character"),
          function(x, rid) {
            convertIDs_character(x, rid)
          })

convertIDs_matrix <- function(use_W_mat, rid) {

  if (rid == "hsy") {
    rownames(use_W_mat) <- toupper(rownames(use_W_mat))
  } else if (rid == "hens") {
    rownames(use_W_mat) <- toupper(rownames(use_W_mat))
    genes <-rownames(use_W_mat)
    mapped_symbols <- na.omit(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = genes,
      keytype = "SYMBOL",
      columns = c("ENSEMBL")
    ))
    use_genes <- table(mapped_symbols$SYMBOL)
    use_genes <- names(use_genes)[use_genes==1]
    use_W_mat <- use_W_mat[use_genes,]
    mapped_symbols <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = use_genes,
      keytype = "SYMBOL",
      columns = c("ENSEMBL")
    )
    rownames(use_W_mat) <- mapped_symbols$ENSEMBL
  } else if (rid == "mens") {
    genes <-rownames(use_W_mat)
    mapped_symbols <- na.omit(AnnotationDbi::select(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = genes,
      keytype = "SYMBOL",
      columns = c("ENSEMBL")
    ))
    use_genes <- table(mapped_symbols$SYMBOL)
    use_genes <- names(use_genes)[use_genes==1]
    use_W_mat <- use_W_mat[use_genes,]
    mapped_symbols <- AnnotationDbi::select(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = use_genes,
      keytype = "SYMBOL",
      columns = c("ENSEMBL")
    )
    rownames(use_W_mat) <- mapped_symbols$ENSEMBL
  }
  return(use_W_mat)

}

convertIDs_character <- function(IDs, rid) {

  if (rid == "hsy") {
    IDs <- toupper(IDs)
  } else if (rid == "hens") {
    IDs <- toupper(IDs)
    genes <-IDs
    mapped_symbols <- na.omit(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = genes,
      keytype = "SYMBOL",
      columns = c("ENSEMBL")
    ))
    IDs <- mapped_symbols$ENSEMBL
  } else if (rid == "mens") {
    genes <- IDs
    mapped_symbols <- na.omit(AnnotationDbi::select(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = genes,
      keytype = "SYMBOL",
      columns = c("ENSEMBL")
    ))
    IDs <- mapped_symbols$ENSEMBL
  }
  return(IDs)

}
