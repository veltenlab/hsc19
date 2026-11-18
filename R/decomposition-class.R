

#'Decomposition class
#'
#'This class stores the results of \code{\link{decompose}}
#'
#'@slot result Matrix with the inferred activity of the GRPs per sample: one row per program, plus a row "(Intercept)"; one column per sample
#'@slot factor_names Names of the factors/programs used for the decomposition
#'@slot pca Seurat \code{DimReduc} object with the result of running PCA (\code{Seurat::RunPCA}, 50 PCs) on the variable features of the input data. Only populated if the data has more than 50 samples
#'@slot pve Fraction of the gene expression variance (summed over the genes in \code{usegenes}) explained by the programs
#'@slot s Seurat object that was used to create the decomposition
#'@slot pve_pca Fraction of the gene expression variance (in \code{usegenes}) explained by as many principal components as there are programs
#'@slot usegenes Genes (in the ID format of the input data) shared between the input data and the loadings, on which the decomposition was computed
#'@export
decomposition <- setClass(
  "decomposition",
  slots = c(
    result = "matrix",
    factor_names = "character",
    pve = "numeric",
    pve_pca = "numeric",
    pca = "DimReduc",
    s = "Seurat",
    usegenes = "character"
  )
)

#' Show a decomposition object
#'
#' Print a short summary of a \code{\link{decomposition}}.
#'
#' @param object A \code{\link{decomposition}} object.
#'
#' @return Called for its side effect of printing; returns \code{NULL} invisibly.
#' @exportMethod show
#' @aliases show,decomposition-method
setMethod(
  "show",
  "decomposition",
  function(object) {
    n_factors <- length(object@factor_names)
    cat(sprintf("Decomposition object of %d samples (%s, %s, ...), run on %d genes. %d programs explain %.1f %% of the gene expression variance in your data. For comparison, %d principal components explain %.1f %%",
        ncol(object@result), colnames(object@result)[1], colnames(object@result)[2], length(object@usegenes), n_factors, 100 * object@pve, n_factors, 100 * object@pve_pca))
  }
)

#' Get column names of a decomposition object
#'
#' @param x A \code{\link{decomposition}} object.
#' @param do.NULL,prefix Passed to \code{\link[base]{colnames}}.
#'
#' @return A character vector of column names of the \code{result} slot, i.e. the sample/cell names.
#' @exportMethod colnames
#' @aliases colnames,decomposition-method
setMethod(
  "colnames",
  signature(x = "decomposition"),
  function(x, do.NULL = TRUE, prefix = "col") {
    colnames(x@result, do.NULL = do.NULL, prefix = prefix)
  }
)

#' Set column names of a decomposition object
#'
#' @param x A \code{\link{decomposition}} object.
#' @param value A character vector of new sample/cell names. Only the \code{result} slot is renamed; the Seurat object in slot \code{s} keeps its original names.
#'
#' @return The updated \code{\link{decomposition}} object.
#' @exportMethod "colnames<-"
#' @aliases colnames<-,decomposition,character-method
setReplaceMethod(
  "colnames",
  signature(x = "decomposition", value = "character"),
  function(x, value) {
    colnames(x@result) <- value
    x
  }
)

#' Subset a decomposition object
#'
#' Keep only a subset of samples/cells in a \code{\link{decomposition}}.
#'
#' @param x A \code{\link{decomposition}} object.
#' @param include Samples/cells to keep. This is used to subset \code{x@result[, include]}
#'   and is passed to \code{cells = include} in \code{\link[Seurat:subset.Seurat]{Seurat::subset()}}.
#'
#' @return A \code{\link{decomposition}} object containing only \code{include}. The slots \code{pve}, \code{pve_pca} and \code{usegenes} are carried over unchanged and are not recomputed for the subset.
#' @seealso \code{\link[base]{subset}}
#' @exportMethod subset
#' @aliases subset,decomposition-method
setMethod(
  "subset",
  signature(x = "decomposition"),
  function(x, include) {
    x@result <- x@result[,include]
    x@s <- subset(x@s, cells = include)
    x@pca <- x@s@reductions$pca
    x
  }
)
