

#'Decomposition class
#'
#'This class stores the results of \code{\link{decompose}}
#'
#'@slot result inferred activity of the 18 GRPs per sample
#'@slot pca result of running `prcomp` on the HVGs of the input data
#'@slot pve Percent variance explained
#'@slot s Seurat object that was used to create the decomposition
#'@slot pve_pca Percent variant explained by 18 principal components
#'@export
decomposition <- setClass(
  "decomposition",
  slots = c(
    result = "matrix",
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

setMethod(
  "show",
  "decomposition",
  function(object) {
    cat(sprintf("Decomposition object of %d samples (%s, %s, ...), run on %d genes. 19 GRPs explain %.1f %% of the gene expression variance in your data. For comparison, 19 principal components explain %.1f %%",
        ncol(object@result), colnames(object@result)[1], colnames(object@result)[2], length(object@usegenes), 100* object@pve, 100* object@pve_pca))
  }
)

setMethod(
  "colnames",
  signature(x = "decomposition"),
  function(x, do.NULL = TRUE, prefix = "col") {
    colnames(x@result, do.NULL = do.NULL, prefix = prefix)
  }
)

## Setter so `colnames(obj) <- value` works
setReplaceMethod(
  "colnames",
  signature(x = "decomposition", value = "character"),
  function(x, value) {
    colnames(x@result) <- value
    x
  }
)


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
