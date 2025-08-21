#'Show genes associated with a certain GRP
#'
#'@param factor Name of a GRP, e.g. \code{"Factor_7"}
#'@param top Whether to show genes with the highest or lowest loadings
#'@param values Whether to return numerical factor loadings
#'@return Character vector of gene symbols, sorted by absolute value of the loading, or named numeric vector
#'@export
get_genes <- function(factor, top = T, values = F) {

  a <- hsc18_loadings[,factor] * factor_annotation$multiply_with[factor_annotation$factor_id == factor]
  if (factor_annotation$multiply_with[factor_annotation$factor_id == factor] == -1) {
    cat("For purpose of display, loadings were multiplied with -1\n")
  }
  if (top) a <- a[a >0] else a <- a[a <0]
  if (values) {
    return(a[order(a, decreasing = top)])
  } else {
    return(names(a)[order(a, decreasing = top)])
  }


}

#'Show perturbations regulating a certain GRP
#'
#'@param factor Name of a GRP, e.g. \code{"Factor_7"}
#'@param top Whether to show positive or negative regulators
#'@param values Whether to return numerical factor loadings
#'@return Character vector of gene symbols, sorted by absolute value of the loading, or named numeric vector
#'@export
get_pert <- function(factor, top = T, values = F) {
a <- hsc18_perturbation_beta[,factor] * factor_annotation$multiply_with[factor_annotation$factor_id == factor]
if (factor_annotation$multiply_with[factor_annotation$factor_id == factor] == -1) {
  cat("For purpose of display, loadings were multiplied with -1\n")
}
if (top) a <- a[a >0] else a <- a[a <0]
if (values) {
  return(a[order(a, decreasing = top)])
} else {
  return(names(a)[order(a, decreasing = top)])
}

}

#'Plot perturbations regulating a certain GRP
#'
#'@param factor Name of a GRP, e.g. \code{"Factor_7"}
#'@param min.sig.gRNA Only include perturbations with this many significant gRNAs
#'@return ggplot2 object
#'@export
plot_pert <- function(factor, min.sig.gRNA = 2) {
  beta_sparse_long <- reshape2::melt(hsc18_perturbation_beta, varnames = c("guide_id", "factor_id"), value.name = "sparse")
  beta_sparse_long$guide_id <- gsub("-CRISPRi-sgRNA-", ".", beta_sparse_long$guide_id)
  beta_full_long <- reshape2::melt(beta_mat, varnames = c("guide_id", "factor_id"), value.name = "full")
  beta_full_long$guide_id <- gsub("-CRISPRi-sgRNA-", ".", beta_full_long$guide_id)
  gamma_long <- reshape2::melt(gamma_mat, varnames = c("guide_id", "factor_id"), value.name = "gamma")
  gamma_long$guide_id <- gsub("-CRISPRi-sgRNA-", ".", gamma_long$guide_id)

  factor_regulators <- merge(factor_annotation, beta_sparse_long)
  factor_regulators <- merge(factor_regulators, beta_full_long)
  factor_regulators <- merge(factor_regulators, gamma_long)
  factor_regulators$gene <- gsub("\\.\\d$", "", factor_regulators$guide_id)
  factor_genes <- plyr::ddply(factor_regulators, c("gene", "factor_id", "factor_label", "multiply_with"), plyr::summarise,
                              n.sig = sum(sparse != 0),
                              rel.sig = sum(sparse != 0) / length(sparse),
                              mean.sig = mean(sparse[sparse!=0]), mean.w = mean(full))

  ggplot2::ggplot(ggplot2::aes(y = reorder(gene, mean.w) , x = mean.w * multiply_with, color = as.factor(n.sig)), data = subset(factor_genes, n.sig >= min.sig.gRNA & factor_id == factor)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(name = "gRNA significant", values = c( "1" = "#888888","2" = "#444444", "3" = "black")) +
    ggplot2::xlab("Regulatory effect") + ggplot2::ylab("") +
    ggplot2::theme_bw( ) + ggplot2::theme(panel.grid = ggplot2::element_blank()) + ggplot2::geom_vline(xintercept = 0, linetype=3)
}
