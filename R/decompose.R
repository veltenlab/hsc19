#'Decompose RNA-seq data into 18 gene regulatory programs
#'
#'Ideally run this function on RNA-seq data of HSCs, e.g. bulk RNA-seq data of CD34+ cells (human) or LSK cells (mouse), or single cell data subset to HSCs. This function can run on bulk, pseudobulk or single cell RNA-seq.
#'@param counts Count matrix, with gene IDs as row names and sample identifiers as column names. Supported gene name formats are mouse or human gene symbols (preferred), or ENSEMBL (ENSG/ENSMUSG) IDs.
#'@param s A seurat object, with \code{NormalizeData} and \code{FindVariableFeatures} run. One of \code{counts} or \code{s} needs to be provided. If \code{counts} is provided, \code{s} is ignored.
#'@return Returns an object of class \code{\link{decomposition}}
#'@export
decompose <- function(counts = NULL, s = NULL) {
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

  #determine format of rownames
  if (mean(grepl("^ENSG", rownames(s)) > 0.9)) {
    rid <- "hens"
  } else if (mean(grepl("^ENSMUSG", rownames(s))) > 0.9) {
    rid <- "mens"
  } else if (mean(grepl("^[^a-z]+$", rownames(s))) > 0.9) {
    rid <- "hsy"
  } else rid <- "msy"

  use_W_mat <- hsc19_loadings
  #convert gene ids
  use_W_mat <- convertIDs(use_W_mat, rid)

  res <- list()
  exprs <- as.matrix(Seurat::GetAssayData(s))
  usegenes = intersect(rownames(use_W_mat), rownames(exprs))
  for (i in 1:ncol(exprs)) {
    modelf <- data.frame(response = exprs[usegenes,i], use_W_mat[usegenes,])
    mo <- lm(response ~ ., data = modelf)
    res[[i]] <- coef(mo)
    attr(res[[i]], "R2") <- summary(mo)$r.squared
  }
  resmat <- do.call(cbind,res)
  colnames(resmat) <- colnames(exprs)

  #percent variance explained
  #denominator: sum of variance of all genes
  pve_denom <- sum(apply(exprs[usegenes,],1,var))
  pve_numer <- sum(apply(use_W_mat[usegenes,], 1, function(x) var(t(resmat) %*% c(1,x))))
  pca <- prcomp(t(exprs[usegenes,]))
  #pca_hvg <- prcomp(t(exprs[VariableFeatures(s),]))
  if (ncol(s) > 50) s <- Seurat::RunPCA(s, verbose=F, npcs = min(c(50, ncol(s))))

  out <- new("decomposition")
  out@result <- resmat
  out@pve <- pve_numer / pve_denom
  if (ncol(s) > 50) out@pca <- s@reductions$pca
  out@s <- s
  out@pve_pca <- cumsum(pca$sdev^2 / sum(pca$sdev^2))[19]
  out@usegenes <- usegenes

  cat(sprintf("Run on %d samples. On the overlapping gene set of %d genes, 19 GRPs explain %.1f %% of the gene expression variance in your data. For comparison, 19 principal components explain %.1f %%\n",
              ncol(resmat), length(out@usegenes), 100* out@pve, 100* out@pve_pca))

  return(out)
}
#' Convert gene IDs
#'
#' Convert mouse gene symbols to different formats
#'
#' @param x A matrix with gene symbols as row names, or a character vector gene symbols.
#' @param rid One of hsy (human gene symbol), hens or mens (ensembl IDs)
#'
#' @return The matrix with ambiguous entries removed, or the character vector
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
