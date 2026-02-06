#' Dataset: GRPs (Loadings)
#'
#' @format A matrix with GRPs in columns and genes in rows. Entries correspond to loadings of the genes on the GRPs.
#' @details
#' Obtained from GSFA. Corresponds to the W matrix, after setting insignificant entries (PIP<0.8) to 0
"hsc19_loadings"

#' Dataset: Effects of Perturbations on GRPs
#'
#' @format A matrix with GRPs in columns and perturbations in rows, detailing the effect of each perturbation on each GRP.
#' @details
#' Obtained from GSFA. Corresponds to the beta matrix, after setting insignificant entries (PIP<0.8) to 0
"hsc19_perturbations_beta"

#' Dataset: GSFA results, W
#'
#' @format Matrices corresponding to GSFA output (W, beta) and corresponding posterior inclusion probabilities (PIPs)
"W_mat"

#' Dataset: GSFA results, PIP
#'
#' @describeIn W_mat PIPs for W
"PIP_mat"

#' Dataset: GSFA results, Beta
#'
#' @describeIn W_mat Beta
"beta_mat"

#' Dataset: GSFA results, Gamma
#'
#' @describeIn W_mat Gamma (Posterior Inclusion Probabilities of Beta)
"gamma_mat"

#' Dataset: Gene lists from Rodriguez-Fraticelli et al., 2020
#'
#' @format List of character vectors, adapted from the supplement of Rodriguez-Fraticelli et al., Nature 2020
"alejo_lists"

#' Dataset: GRP Annotation
#'
#' @format Data frame contaning annotation of the GRPs
"factor_annotation"