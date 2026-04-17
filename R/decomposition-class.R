

#'Decomposition class
#'
#'This class stores the results of \code{\link{decompose}}
#'
#'@slot result inferred activity of the 18 GRPs per sample
#'@slot factor_names Names of the factors/programs used for the decomposition
#'@slot pca result of running `prcomp` on the HVGs of the input data
#'@slot pve Percent variance explained
#'@slot s Seurat object that was used to create the decomposition
#'@slot pve_pca Percent variant explained by 18 principal components
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
  ),
  validity = function(object) {
    if (!inherits(object@pca, "prcomp")) {
      return("Slot 'pca' must be an object of class 'prcomp'")
    }
    TRUE
  }

)

#' Show a decomposition object
#'
#' Print a short summary of a \code{\link{decomposition}}.
#'
#' @param object A \code{\link{decomposition}} object.
#'
#' @return Invisibly returns \code{object}.
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
#' @param do.NULL, prefix Passed to \code{\link[base]{colnames}}.
#'
#' @return A character vector of column names.
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
#' @param value A character vector of new column names.
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
#' @return A \code{\link{decomposition}} object containing only \code{include}.
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
