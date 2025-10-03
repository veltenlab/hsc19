#in the 500 gene GRPs, c(1:16,18,20), otherwise 3:20

#'Associate GRPs with a numeric outcome
#'
#'This function takes a decomposition object and a metadata frame with a numeric outcome and, if available, other predictors. It will then a) compute statistical associations between
#'all predictors and the outcome and b) compare the performance of all predictors to that of various supervised machine learning strategies (see details)
#'
#'@param decomp Object of class \code{\link{decomposition}}, output from \code{\link{decompose}}
#'@param metadata A data frame with rownames corresponding to column names of \code{decomp} that contains the outcome of interest and any other predictors that should be tested
#'@param outcome.name Column name from \code{metadata} that contains the outcome of interest
#'@param other.predictors Names of any other columns in \code{metadata} that may serve as 0-shot predictors for the outcome
#'@param run.ML Whether to run machine-learning models.
#'@param genelists List of character vectors of gene names. If \code{run.ML=TRUE}, LASSO regression using these genes as input will be performed and compared to the other models.
#'@param replicates Number of times to repeat CV calculation with different train-test splits. Only used if \code{run.ML=TRUE}
#'@param cores Number of cores to use (if replicates > 1)
#'@param use.s Regularization to use for LASSO, one of \code{lambda.min} and \code{lambda.1se}
#'
#'@details	If \code{run.ML} is set, the function will set up a 10 fold cross validation scheme (or leave one out CV if <50 observations). On the train set, it will perform the following:
#'\itemize{
#'\item{Perform principal component regression on 18 PCs computed on highly variable genes from the underlying Seurat object}
#'\item{Perform LASSO regularized principal component regression on 50 PCs computed on highly variable genes from the underlying Seurat object}
#'\item{Perform DE testing between samples from the highest and lowest quartile of the outcome, and perform LASSO or linear regression on the top DE features}
#'\item{Perform LASSO and linear regression on the 18 GRPs from the decomposition}
#'\item{If \code{genelists} are provided, perform LASSO on the genes from the gene lists}
#'}
#'
#'@return A list with entries \code{p.values} and, if \code{run.ML} is set, \code{r.squares}, their standard deviation across replicates, and info on the different folds
#'@export
associate_numeric_outcome<- function(decomp, metadata, outcome.name, other.predictors=NULL, run.ML=T, genelists = NULL,  cores = 6, replicates = 3, use.s = "lambda.min") {

  metadata <- metadata[colnames(decomp),]
  modelf <- cbind(metadata, t(decomp@result))
  if (any(is.na(modelf[,outcome.name]))) {
    warning("Removed ", sum(is.na(modelf[,outcome.name])), "NA entries for ", outcome.name)
    modelf <- modelf[!is.na(modelf[,outcome.name]),]
  }
  p.values <- data.frame(factor_id = factor(c(sprintf("Factor_%d", 3:20),other.predictors),levels = c(sprintf("Factor_%d", 3:20),other.predictors)),
                          p_outcome = apply(modelf[,c(sprintf("Factor_%d", 3:20),other.predictors)], 2, function(x) cor.test(x, modelf[,outcome.name])$p.value),
                          sign_outcome = sign(apply(modelf[,c(sprintf("Factor_%d", 3:20),other.predictors)], 2, function(x) cor.test(x, modelf[,outcome.name])$estimate)))

  best.factor <- with(subset(p.values, factor_id %in% sprintf("Factor_%d", 3:20)), as.character(factor_id[which.min(p_outcome)]))
  if (!is.null(other.predictors)) best.other <- with(subset(p.values, factor_id %in% other.predictors), as.character(factor_id[which.min(p_outcome)]))

  if(run.ML) {
    modelf <- cbind(modelf, decomp@pca@cell.embeddings)

    cv_result <- parallel::mclapply(1:replicates, function(void) get_cv(decomp, modelf, outcome.name, other.predictors, genelists, use.s), mc.cores= cores)
    r.squares <- do.call(rbind,lapply(cv_result, "[[",1))
    best.factor <- table(unlist(lapply(cv_result, "[[",2)))
    best.other <- table(unlist(lapply(cv_result, "[[",3)))
    return(list(p.values = p.values, r.squares = colMeans(r.squares), r.squares.sd = apply(r.squares, 2,sd), best.factor = best.factor, best.other = best.other))
    # r.squares <- get_cv(decomp, modelf, outcome.name, other.predictors, genelists, use.s)
    # return(list(p.values = p.values, r.squares = r.squares))



  } else {
    return(list(p.values = p.values))
  }

}

#'Associate GRPs with a binary outcome
#'
#'This function takes a decomposition object and a metadata frame with a binary (TRUE/FALSE) outcome and, if available, other predictors. It will then a) compute statistical associations between
#'all predictors and the outcome and b) compare the performance of all predictors to that of various supervised machine learning strategies (see details)
#'
#'@param decomp Object of class \code{\link{decomposition}}, output from \code{\link{decompose}}
#'@param metadata A data frame with rownames corresponding to column names of \code{decomp} that contains the outcome of interest and any other predictors that should be tested
#'@param outcome.name Column name from \code{metadata} that contains the outcome of interest
#'@param other.predictors Names of any other columns in \code{metadata} that may serve as 0-shot predictors for the outcome
#'@param run.ML Whether to run machine-learning models.
#'@param arbitrary.models Arbitrary linear regression models, based on columns of \code{metadata}, that should also be included in the comparison
#'@param genelists List of character vectors of gene names. If \code{run.ML=TRUE}, LASSO regression using these genes as input will be performed and compared to the other models.
#'@param replicates Number of times to repeat CV calculation with different train-test splits. Only used if \code{run.ML=TRUE}
#'@param cores Number of cores to use (if replicates > 1)
#'@param use.s Regularization to use for LASSO, one of \code{lambda.min} and \code{lambda.1se}
#'@param ROCs If full ROC curves should be returned (if FALSE, only returns AUCs)
#'
#'@details	If \code{run.ML} is set, the function will set up a 10 fold cross validation scheme (or leave one out CV if <50 observations). On the train set, it will perform the following:
#'\itemize{
#'\item{Perform principal component regression on 18 PCs computed on highly variable genes from the underlying Seurat object}
#'\item{Perform LASSO regularized principal component regression on 50 PCs computed on highly variable genes from the underlying Seurat object}
#'\item{Perform DE testing between samples from the highest and lowest quartile of the outcome, and perform LASSO or linear regression on the top DE features}
#'\item{Perform LASSO and linear regression on the 18 GRPs from the decomposition}
#'\item{If \code{genelists} are provided, perform LASSO on the genes from the gene lists}
#'}
#'
#'@return A list with entries \code{p.values} and, if \code{run.ML} is set, \code{aucs}, their standard deviation across replicates, and info on the different folds
#'@export
associate_binary_outcome<- function(decomp, metadata, outcome.name, other.predictors=NULL, arbitrary.models = NULL, genelists = NULL, run.ML=T, cores = 6, replicates = 3, use.s = "lambda.min", ROCs = F) {

  metadata <- metadata[colnames(decomp),]
  modelf <- cbind(metadata, t(decomp@result))
  if (any(is.na(modelf[,outcome.name]))) {
    warning("Removed ", sum(is.na(modelf[,outcome.name])), "NA entries for ", outcome.name)
    modelf <- modelf[!is.na(modelf[,outcome.name]),]
  }
  p.values <- data.frame(factor_id = factor(c(sprintf("Factor_%d", 3:20),other.predictors),levels = c(sprintf("Factor_%d", 3:20),other.predictors)),
                         p_outcome = apply(modelf[,c(sprintf("Factor_%d", 3:20),other.predictors)], 2, function(x) wilcox.test(x[modelf[,outcome.name]],x[!modelf[,outcome.name]])$p.value),
                         sign_outcome = sign(apply(modelf[,c(sprintf("Factor_%d", 3:20),other.predictors)], 2, function(x) wilcox.test(x[modelf[,outcome.name]],x[!modelf[,outcome.name]], conf.int = TRUE)$estimate)))

  best.factor <- with(subset(p.values, factor_id %in% sprintf("Factor_%d", 3:20)), as.character(factor_id[which.min(p_outcome)]))
  if (!is.null(other.predictors)) best.other <- with(subset(p.values, factor_id %in% other.predictors), as.character(factor_id[which.min(p_outcome)]))

  if(run.ML) {
    modelf <- cbind(modelf, decomp@pca@cell.embeddings)

    cv_result <- parallel::mclapply(1:replicates, function(void) suppressWarnings(get_cv_binary(decomp, modelf, outcome.name, other.predictors, genelists, arbitrary.models, use.s, ROCs)), mc.cores= cores)
    aucs <- do.call(rbind,lapply(cv_result, "[[", 1))
    best.factor <- table(unlist(lapply(cv_result, "[[",2)))
    best.other <- table(unlist(lapply(cv_result, "[[",3)))
    if (ROCs)  {
      rocs <- lapply(cv_result, "[[", 4)
      return(list(p.values = p.values, aucs = colMeans(aucs), aucs.sd = apply(aucs, 2,sd), best.factor = best.factor, best.other = best.other, rocs = rocs))
    } else {
      return(list(p.values = p.values, aucs = colMeans(aucs), aucs.sd = apply(aucs, 2,sd), best.factor = best.factor, best.other = best.other))
    }

    # aucs <- get_cv_binary(decomp, modelf, outcome.name, other.predictors, genelists, use.s)
    # return(list(p.values = p.values, aucs = aucs))



  } else {
    return(list(p.values = p.values))
  }

}


#'Associate GRPs with a survival outcome
#'
#'This function takes a decomposition object and a metadata frame with a binary (TRUE/FALSE) outcome and, if available, other predictors. It will then a) compute statistical associations between
#'all predictors and the outcome and b) compare the performance of all predictors to that of various supervised machine learning strategies (see details)
#'
#'@param decomp Object of class \code{\link{decomposition}}, output from \code{\link{decompose}}
#'@param metadata A data frame with rownames corresponding to column names of \code{decomp} that contains the outcome of interest and any other predictors that should be tested
#'@param outcome.status.name Column name from \code{metadata} that contains the survival outcome (0 = censored, 1 = event)
#'@param outcome.time.name Column name from \code{metadata} that contains the associated follow-up time
#'@param other.predictors Names of any other columns in \code{metadata} that may serve as 0-shot predictors for the outcome
#'@param arbitrary.models Arbitrary linear regression models, based on columns of \code{metadata}, that should also be included in the comparison.
#'@param run.ML Whether to run machine-learning models.
#'@param genelists List of character vectors of gene names. If \code{run.ML=TRUE}, LASSO regression using these genes as input will be performed and compared to the other models.
#'@param replicates Number of times to repeat CV calculation with different train-test splits. Only used if \code{run.ML=TRUE}
#'@param cores Number of cores to use (if replicates > 1)
#'@param use.s Regularization to use for LASSO, one of \code{lambda.min} and \code{lambda.1se}
#'
#'@details	If \code{run.ML} is set, the function will set up a 10 fold cross validation scheme (or leave one out CV if <50 observations). On the train set, it will perform the following:
#'\itemize{
#'\item{Perform principal component regression on 18 PCs computed on highly variable genes from the underlying Seurat object}
#'\item{Perform LASSO regularized principal component regression on 50 PCs computed on highly variable genes from the underlying Seurat object}
#'\item{Perform DE testing between samples from the highest and lowest quartile of the outcome, and perform LASSO or linear regression on the top DE features}
#'\item{Perform LASSO and linear regression on the 18 GRPs from the decomposition}
#'\item{If \code{genelists} are provided, perform LASSO on the genes from the gene lists}
#'}
#'
#'@return A list with entries \code{p.values} and, if \code{run.ML} is set, \code{aucs}, their standard deviation across replicates, and info on the different folds
#'@export
associate_survival_outcome <- function(
    decomp,
    metadata,
    outcome.status.name,   # 0 = censored, 1 = event
    outcome.time.name,     # follow-up time
    other.predictors = NULL,
    arbitrary.models = NULL,
    genelists = NULL,
    run.ML = TRUE,
    cores = 6,
    replicates = 3,
    use.s = "lambda.min"
) {
  stopifnot(all(c(outcome.status.name, outcome.time.name) %in% colnames(metadata)))

  # align to decomp cells and build model frame with factors
  metadata <- metadata[colnames(decomp), , drop = FALSE]
  modelf <- cbind(metadata, t(decomp@result))

  # remove NAs in time or status
  na_mask <- is.na(modelf[, outcome.time.name]) | is.na(modelf[, outcome.status.name])
  if (any(na_mask)) {
    warning("Removed ", sum(na_mask),
            " rows with NA in ", outcome.time.name, " or ", outcome.status.name)
    modelf <- modelf[!na_mask, , drop = FALSE]
  }

  # basic check on coding
  if (!all(modelf[[outcome.status.name]] %in% c(0, 1))) {
    stop("Status variable must be coded 0/1 (0=censored, 1=event).")
  }

  # ---- Univariate Cox p-values for candidate predictors ----
  cand <- c(sprintf("Factor_%d", 3:20), other.predictors)
  cand <- cand[cand %in% colnames(modelf)]

  get_stats <- function(cl) {
    f <- as.formula(paste0("survival::Surv(", outcome.time.name, ",", outcome.status.name, ") ~ ", cl))
    fit <- try(survival::coxph(f, data = modelf, ties = "efron"), silent = TRUE)
    if (inherits(fit, "try-error")) return(c(p = NA_real_, coef = NA_real_))
    co <- try(summary(fit)$coefficients, silent = TRUE)
    if (inherits(co, "try-error") || nrow(co) < 1) return(c(p = NA_real_, coef = NA_real_))
    c(p = as.numeric(co[1, "Pr(>|z|)"]),
      coef = as.numeric(co[1, "coef"]))  # sign via coef (HR = exp(coef))
  }

  stats <- if (length(cand)) {
    vapply(cand, get_stats, FUN.VALUE = c(p = NA_real_, coef = NA_real_))
  } else {
    matrix(numeric(0), nrow = 2, dimnames = list(c("p", "coef"), NULL))
  }

  p.values <- data.frame(
    factor_id    = factor(cand, levels = c(sprintf("Factor_%d", 3:20), other.predictors)),
    p_outcome    =  as.numeric(stats["p", ]),
    sign_outcome = sign(as.numeric(stats["coef", ])),
    row.names    = cand
  )

  if (run.ML) {
    # add PCs (assumes decomp@pca@cell.embeddings columns named like PC_1, PC_2, …)
    modelf <- cbind(modelf, decomp@pca@cell.embeddings)

    # run replicated CV (returns out-of-fold C-indices from get_cv_surv)
    cv_result <- parallel::mclapply(X = 1:replicates, FUN = function(void) suppressWarnings(get_cv_surv(decomp, modelf, outcome.status.name, outcome.time.name,
                    other.predictors, genelists, arbitrary.models, use.s)),mc.cores = cores)
    cis <- do.call(rbind,lapply(cv_result, "[[",1))
    best.factor <- table(unlist(lapply(cv_result, "[[",2)))
    best.other <- table(unlist(lapply(cv_result, "[[",3)))
    # cis <- lapply(X = 1:replicates, FUN = function(void) get_cv_surv(decomp, modelf, outcome.status.name, outcome.time.name,
    #                                                                              other.predictors, genelists, arbitrary.models, use.s))
    # cis <- do.call(rbind, cis)
    out <- tryCatch( {
      list(
      p.values    = p.values,
      c.index     = colMeans(cis, na.rm = TRUE),
      c.index.sd  = apply(cis, 2, sd, na.rm = TRUE),
      best.factor = best.factor,
      best.other = best.other
      )
    }, error = function(e) {
      warning("Error in parallel execution: ", e, "\n returning parallel output...\n")
      cv_result
    }
    )
    return(out)
  } else {
    return(list(p.values = p.values))
  }
}


get_cv <- function(decomp, modelf, outcome.name, other.predictors, genelists, use.s) {
  modelf <- modelf[sample(1:nrow(modelf)),]
  #set up train-test splits
  if (nrow(modelf) > 50) {
    modelf$set <- rep(1:10, length.out = nrow(modelf))
  } else {
    modelf$set <- 1:nrow(modelf)
  }

  modelf$predicted_factors <- NA
  modelf$predicted_factors_lasso <- NA
  modelf$predicted_best_factor <- NA
  modelf$predicted_best_other <- NA
  modelf$predicted_pcr <- NA
  modelf$predicted_pcr_lasso <- NA
  modelf$predicted_de <- NA
  modelf$predicted_de_lasso <- NA
  predicted_genelists <- replicate(length(genelists), rep(NA, nrow(modelf)), simplify = F)
  names(predicted_genelists) <- names(genelists)
  best_factors <- c()
  best_others <- c()
  for (i in 1:max(modelf$set)) {
    cat(i, "\n")
    train <- subset(modelf, set != i)
    test <- subset(modelf, set == i)
    #determine the best "0shot" predictors
    p.values.internal <- data.frame(factor_id = factor(c(sprintf("Factor_%d", 3:20),other.predictors),levels = c(sprintf("Factor_%d", 3:20),other.predictors)),
                                    p_outcome = apply(train[,c(sprintf("Factor_%d", 3:20),other.predictors)], 2, function(x) cor.test(x, train[,outcome.name])$p.value))

    best.factor <- with(subset(p.values.internal, factor_id %in% sprintf("Factor_%d", 3:20)), as.character(factor_id[which.min(p_outcome)]))
    best_factors <- c(best_factors, best.factor)
    if (!is.null(other.predictors)) {
      best.other <- with(subset(p.values.internal, factor_id %in% other.predictors), as.character(factor_id[which.min(p_outcome)]))
      best_others <- c(best_others, best.other)
    }

    #0 shot predictions
    mo_best_factor <- lm(as.formula(sprintf("%s~%s", outcome.name, best.factor)), data = train)
    if (!is.null(other.predictors)) mo_best_other <- lm(as.formula(sprintf("%s~%s", outcome.name, best.other)), data = train)

    modelf$predicted_best_factor[modelf$set==i] <- predict(mo_best_factor, test)
    if (!is.null(other.predictors))  modelf$predicted_best_other[modelf$set==i] <- predict(mo_best_other, test)

    #regression with factors
    mo_factors <- lm(as.formula(sprintf("%s~%s", outcome.name, paste(sprintf("Factor_%d",3:20),collapse="+"))), data = train)
    mo_factors_lasso <- glmnet::cv.glmnet(as.matrix(train[,sprintf("Factor_%d",3:20)]), train[,outcome.name])
    modelf$predicted_factors[modelf$set==i] <- predict(mo_factors, test)
    modelf$predicted_factors_lasso[modelf$set==i] <- predict(mo_factors_lasso, as.matrix(test[,sprintf("Factor_%d",3:20)]), s = use.s)

    #PCR
    mo_pcr_lasso <- glmnet::cv.glmnet(as.matrix(train[,sprintf("PC_%d",1:30)]), train[,outcome.name])
    mo_pcr <- lm(as.formula(sprintf("%s~%s", outcome.name, paste(sprintf("PC_%d",1:18),collapse="+"))), data = train)
    modelf$predicted_pcr_lasso[modelf$set==i] <- predict(mo_pcr_lasso, as.matrix(test[,sprintf("PC_%d",1:30)]), s = use.s)
    modelf$predicted_pcr[modelf$set==i] <- predict(mo_pcr, test)

    #DE Lasso
    top <- rownames(train)[train[,outcome.name] > quantile(train[,outcome.name], 0.75)]
    bottom <- rownames(train)[train[,outcome.name] <= quantile(train[,outcome.name], 0.25)]
    forde <- subset(decomp@s, cells = c(top,bottom))
    forde$group <- colnames(forde) %in% top
    Idents(forde) <- ifelse(forde$group, "top", "bottom")
    #browser()

    de <- Seurat::FindMarkers(forde, ident.1 = "top", min.pcr = 0.2, logfc.threshold = 0.5)
    de <- de[order(de$p_val),]
    usegenes <- rownames(de)[1:50]

    ma <- t(as.matrix(GetAssayData(decomp@s)[usegenes,rownames(train)]))
    train <- cbind(train, ma)
    mo_de_lasso <- glmnet::cv.glmnet(ma, train[,outcome.name])
    mo_de <- lm(as.formula(sprintf("%s~%s", outcome.name, paste(sprintf("`%s`",usegenes[3:20]),collapse="+"))), data = train)
    ma_test <- t(as.matrix(GetAssayData(decomp@s)[usegenes,rownames(test)]))
    test <- cbind(test, ma_test)
    modelf$predicted_de_lasso[modelf$set==i] <- predict(mo_de_lasso, ma_test, s = use.s)
    modelf$predicted_de[modelf$set==i] <- predict(mo_de, test)

    #gene list lasso
    #browser()
    if (length(genelists) > 0) {

    for (gl in 1:length(genelists)) {
      usegenes <- intersect(genelists[[gl]], rownames(decomp@s))
      ma <- t(as.matrix(GetAssayData(decomp@s)[usegenes,rownames(train)]))
      mo_gllasso <- glmnet::cv.glmnet(as.matrix(ma), train[,outcome.name])
      ma_test <- t(as.matrix(GetAssayData(decomp@s)[usegenes,rownames(test)]))
      predicted_genelists[[gl]][modelf$set==i] <- predict(mo_gllasso, as.matrix(ma_test))
    }
    }





  }



  r.squares <-c(R2.best.factor = getr2(modelf[,"predicted_best_factor"], modelf[,outcome.name]),
                R2.best.other = getr2(modelf[,"predicted_best_other"], modelf[,outcome.name]),
                R2.lm.factors = getr2(modelf[,"predicted_factors"], modelf[,outcome.name]),
                R2.lasso.factors = getr2(modelf[,"predicted_factors_lasso"], modelf[,outcome.name]),
                R2.lm.de = getr2(modelf[,"predicted_de"], modelf[,outcome.name]),
                R2.lasso.de = getr2(modelf[,"predicted_de_lasso"], modelf[,outcome.name]),
                R2.lm.pcr= getr2(modelf[,"predicted_pcr"], modelf[,outcome.name]),
                R2.lasso.pcr= getr2(modelf[,"predicted_pcr_lasso"], modelf[,outcome.name]))
  if (length(genelists) > 0) {
    r.squares <-c(r.squares, sapply(predicted_genelists, function(x) getr2(x, modelf[,outcome.name])))
  }

  return(list(r.squares = r.squares, best.factor = best_factors, best.other = best_others))
}

getr2 <- function(a,b) {
  r2 <- cor(a,b)^2
  ifelse(r2>0,r2,0)
}

# helper: compute ROC AUC safely (returns NA if all-NA or single class)
get_auc <- function(probs, y) {
  if (all(is.na(probs))) return(NA_real_)
  # coerce y to factor with positive level last (pROC uses the second level as "event")
  if (is.logical(y)) y <- as.integer(y)
  if (is.numeric(y)) y <- factor(y, levels = c(0,1))
  if (!inherits(y, "factor") || length(levels(y)) != 2) stop("Outcome must be binary.")
  if (length(unique(y[!is.na(y)])) < 2) return(NA_real_)
  as.numeric(pROC::auc(y, probs, quiet = TRUE))
}

get_full_auc <- function(probs, y) {
  if (all(is.na(probs))) return(NA_real_)
  # coerce y to factor with positive level last (pROC uses the second level as "event")
  if (is.logical(y)) y <- as.integer(y)
  if (is.numeric(y)) y <- factor(y, levels = c(0,1))
  if (!inherits(y, "factor") || length(levels(y)) != 2) stop("Outcome must be binary.")
  if (length(unique(y[!is.na(y)])) < 2) return(NA_real_)
  pROC::roc(response = y, predictor = probs, quiet = TRUE)
}

get_cv_binary <- function(decomp, modelf, outcome.name, other.predictors, genelists, arbitrary.models, use.s, ROCs) {
  # ensure binary coding and no NAs in outcome
  y <- modelf[, outcome.name]
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y)) {
    if (nlevels(y) != 2) stop("Outcome factor must have exactly two levels.")
    y <- as.integer(y == levels(y)[2])
  }
  if (!all(y %in% c(0,1))) stop("Outcome must be binary (0/1), logical, or 2-level factor.")
  modelf[, outcome.name] <- y

  modelf <- modelf[sample(1:nrow(modelf)), ]

  # 10-fold-ish split (or leave-one-out if small)
  if (nrow(modelf) > 50) {
    modelf$set <- rep(1:10, length.out = nrow(modelf))
  } else {
    modelf$set <- 1:nrow(modelf)
  }

  # store probabilities (0..1)
  modelf$predicted_factors <- NA_real_
  modelf$predicted_factors_lasso <- NA_real_
  modelf$predicted_best_factor <- NA_real_
  modelf$predicted_best_other <- NA_real_
  modelf$predicted_pcr <- NA_real_
  modelf$predicted_pcr_lasso <- NA_real_
  modelf$predicted_de <- NA_real_
  modelf$predicted_de_lasso <- NA_real_
  predicted_genelists <- replicate(length(genelists), rep(NA_real_, nrow(modelf)), simplify = FALSE)
  names(predicted_genelists) <- names(genelists)
  predicted_arbitrary <- replicate(length(arbitrary.models), rep(NA_real_, nrow(modelf)), simplify = FALSE)
  names(predicted_arbitrary) <- names(arbitrary.models)

  best_factors <- c()
  best_others <- c()
  for (i in 1:max(modelf$set)) {
    cat(i, "\n")
    train <- subset(modelf, set != i)
    test  <- subset(modelf, set == i)

    # choose best single predictors via univariate logistic regression p-values
    cand.cols <- c(sprintf("Factor_%d", 3:20), other.predictors)
    cand.cols <- cand.cols[cand.cols %in% colnames(train)]
    pvals <- sapply(cand.cols, function(cl) {
      f <- as.formula(paste0(outcome.name, "~", cl))
      suppressWarnings(summary(glm(f, data = train, family = binomial()))$coefficients[2,4])
    })
    p.values.internal <- data.frame(
      factor_id = factor(cand.cols, levels = cand.cols),
      p_outcome = pvals, row.names = cand.cols
    )

    fac.only <- sprintf("Factor_%d", 3:20)
    best.factor <- with(subset(p.values.internal, factor_id %in% fac.only),
                        as.character(factor_id[which.min(p_outcome)]))
    best_factors <- c(best_factors, best.factor)


    best.other <- if (!is.null(other.predictors)) {
      with(subset(p.values.internal, factor_id %in% other.predictors),
           as.character(factor_id[which.min(p_outcome)]))
    } else NULL
    best_others <- c(best_others, best.other)
    # 0-shot predictions (single predictors) – logistic
    if (length(best.factor)) {
      mo_best_factor <- glm(as.formula(sprintf("%s~%s", outcome.name, best.factor)),
                            data = train, family = binomial())
      modelf$predicted_best_factor[modelf$set == i] <-
        predict(mo_best_factor, test, type = "response")
    }
    if (!is.null(best.other) && length(best.other)) {
      mo_best_other <- glm(as.formula(sprintf("%s~%s", outcome.name, best.other)),
                           data = train, family = binomial())
      modelf$predicted_best_other[modelf$set == i] <-
        predict(mo_best_other, test, type = "response")
    }

    # logistic with factors
    fac.cols <- sprintf("Factor_%d", 3:20)
    fac.cols <- fac.cols[fac.cols %in% colnames(train)]
    if (length(fac.cols)) {
      mo_factors <- glm(as.formula(sprintf("%s~%s", outcome.name,
                                           paste(fac.cols, collapse = "+"))),
                        data = train, family = binomial())
      modelf$predicted_factors[modelf$set == i] <- predict(mo_factors, test, type = "response")

      mo_factors_lasso <- glmnet::cv.glmnet(as.matrix(train[, fac.cols, drop = FALSE]),
                                            train[, outcome.name],
                                            family = "binomial")
      modelf$predicted_factors_lasso[modelf$set == i] <-
        as.numeric(predict(mo_factors_lasso,
                           as.matrix(test[, fac.cols, drop = FALSE]),
                           s = use.s, type = "response"))
    }

    # PCR (logistic on PCs)
    pc.cols <- sprintf("PC_%d", 1:30)
    pc.cols <- pc.cols[pc.cols %in% colnames(train)]
    if (length(pc.cols)) {
      mo_pcr_lasso <- glmnet::cv.glmnet(as.matrix(train[, pc.cols, drop = FALSE]),
                                        train[, outcome.name], family = "binomial")
      modelf$predicted_pcr_lasso[modelf$set == i] <-
        as.numeric(predict(mo_pcr_lasso,
                           as.matrix(test[, pc.cols, drop = FALSE]),
                           s = use.s, type = "response"))

      pc.lm.cols <- intersect(sprintf("PC_%d", 1:18), pc.cols)
      if (length(pc.lm.cols)) {
        mo_pcr <- glm(as.formula(sprintf("%s~%s", outcome.name,
                                         paste(pc.lm.cols, collapse = "+"))),
                      data = train, family = binomial())
        modelf$predicted_pcr[modelf$set == i] <- predict(mo_pcr, test, type = "response")
      }
    }

    # DE-based features: for binary y, 'top' = class 1, 'bottom' = class 0
    top <- rownames(train)[train[, outcome.name] == 1]
    bottom <- rownames(train)[train[, outcome.name] == 0]
    if (length(top) > 1 && length(bottom) > 1) {
      forde <- subset(decomp@s, cells = c(top, bottom))
      forde$group <- colnames(forde) %in% top
      Idents(forde) <- ifelse(forde$group, "top", "bottom")
      de <- Seurat::FindMarkers(forde, ident.1 = "top", min.pcr = 0.2, logfc.threshold = 0.5)
      de <- de[order(de$p_val), , drop = FALSE]
      usegenes <- head(rownames(de), 50)

      if (length(usegenes)) {
        ma_tr <- t(as.matrix(GetAssayData(decomp@s)[usegenes, rownames(train), drop = FALSE]))
        ma_te <- t(as.matrix(GetAssayData(decomp@s)[usegenes, rownames(test),  drop = FALSE]))

        # lasso
        mo_de_lasso <- glmnet::cv.glmnet(ma_tr, train[, outcome.name], family = "binomial")
        modelf$predicted_de_lasso[modelf$set == i] <-
          as.numeric(predict(mo_de_lasso, ma_te, s = use.s, type = "response"))

        # logistic on top DE genes (use up to 20)
        use20 <- head(usegenes, 20)
        tr_with <- cbind(train, ma_tr)
        te_with <- cbind(test,  ma_te)
        mo_de <- glm(as.formula(sprintf("%s~%s", outcome.name,
                                        paste(sprintf("`%s`", use20), collapse = "+"))),
                     data = tr_with, family = binomial())
        modelf$predicted_de[modelf$set == i] <- predict(mo_de, te_with, type = "response")
      }
    }

    # gene list lasso (binomial)
    if (length(genelists) > 0) {
      for (gl in seq_along(genelists)) {
        usegenes <- intersect(genelists[[gl]], rownames(decomp@s))
        if (length(usegenes) == 0) next
        ma_tr <- t(as.matrix(GetAssayData(decomp@s)[usegenes, rownames(train), drop = FALSE]))
        ma_te <- t(as.matrix(GetAssayData(decomp@s)[usegenes, rownames(test),  drop = FALSE]))
        mo_gllasso <- glmnet::cv.glmnet(as.matrix(ma_tr), train[, outcome.name], family = "binomial")
        predicted_genelists[[gl]][modelf$set == i] <-
          as.numeric(predict(mo_gllasso, as.matrix(ma_te), s = use.s, type = "response"))
      }
    }

    #7. arbitrary models
    if (length(arbitrary.models) > 0) {
      for (mi in seq_along(arbitrary.models)) {

        f <- as.formula(paste( outcome.name, "~",arbitrary.models[mi]))
        mo_ab <- glm(f, data = train, family = binomial())

        if (grepl("WHO_classification", arbitrary.models[mi])) {
        #handle special case: Force test factor to training levels; unseen levels become NA
        test$WHO_classification <- factor(test$WHO_classification,
                                          levels = unique(train$WHO_classification))

        # Predict only where the factor is not NA
        ok <- !is.na(test$WHO_classification)
        predicted_arbitrary[[mi]][modelf$set == i] <- rep(NA_real_, nrow(test))
        predicted_arbitrary[[mi]][modelf$set == i][ok] <- as.numeric(predict(mo_ab, newdata = test[ok, ], type = "response"))
        } else {
          predicted_arbitrary[[mi]][modelf$set == i] <- as.numeric(predict(mo_ab, newdata = test, type = "response"))
        }

      }
    }

  }

  # Compute ROC AUCs for each model's out-of-fold probabilities
  y_all <- modelf[, outcome.name]


  aucs <- c(
    AUC.best.factor   = get_auc(modelf[, "predicted_best_factor"],   y_all),
    AUC.best.other    = get_auc(modelf[, "predicted_best_other"],    y_all),
    AUC.logit.factors = get_auc(modelf[, "predicted_factors"],       y_all),
    AUC.lasso.factors = get_auc(modelf[, "predicted_factors_lasso"], y_all),
    AUC.logit.de      = get_auc(modelf[, "predicted_de"],            y_all),
    AUC.lasso.de      = get_auc(modelf[, "predicted_de_lasso"],      y_all),
    AUC.logit.pcr     = get_auc(modelf[, "predicted_pcr"],           y_all),
    AUC.lasso.pcr     = get_auc(modelf[, "predicted_pcr_lasso"],     y_all)
  )

  if (length(genelists) > 0) {
    aucs <- c(aucs, sapply(predicted_genelists, function(x) get_auc(x, y_all)))
  }
  if (length(arbitrary.models) > 0) {
    aucs <- c(aucs, sapply(predicted_arbitrary, function(x) get_auc(x, y_all)))
  }

  if (ROCs) {
    rocs <- list(
      AUC.best.factor   = get_full_auc(modelf[, "predicted_best_factor"],   y_all),
      AUC.best.other    = get_full_auc(modelf[, "predicted_best_other"],    y_all),
      AUC.logit.factors = get_full_auc(modelf[, "predicted_factors"],       y_all),
      AUC.lasso.factors = get_full_auc(modelf[, "predicted_factors_lasso"], y_all),
      AUC.logit.de      = get_full_auc(modelf[, "predicted_de"],            y_all),
      AUC.lasso.de      = get_full_auc(modelf[, "predicted_de_lasso"],      y_all),
      AUC.logit.pcr     = get_full_auc(modelf[, "predicted_pcr"],           y_all),
      AUC.lasso.pcr     = get_full_auc(modelf[, "predicted_pcr_lasso"],     y_all)
    )

    if (length(genelists) > 0) {
      rocs <- c(rocs, lapply(predicted_genelists, function(x) get_full_auc(x, y_all)))
    }
    if (length(arbitrary.models) > 0) {
      rocs <- c(rocs, lapply(predicted_arbitrary, function(x) get_full_auc(x, y_all)))
    }

    return(list(aucs = aucs, best.factor = best_factors, best.other = best_others, rocs = rocs))
  } else {
    return(list(aucs = aucs, best.factor = best_factors, best.other = best_others))
  }


}

# helper: Harrell's C-index on out-of-fold predictions
get_cindex <- function(risk, time, status) {
  ok <- !(is.na(risk) | is.na(time) | is.na(status))
  if (sum(ok) < 3) return(NA_real_)
  a <- as.numeric(survival::concordance(survival::Surv(time[ok], status[ok]) ~ risk[ok])$concordance)
  mrisk <- -1 * risk
  b <- as.numeric(survival::concordance(survival::Surv(time[ok], status[ok]) ~ mrisk[ok])$concordance)
  max(c(a,b))
}

get_cv_surv <- function(decomp, modelf, outcome.status.name, outcome.time.name,
                        other.predictors, genelists, arbitrary.models, use.s) {
  stopifnot(all(c(outcome.status.name, outcome.time.name) %in% names(modelf)))
  if (!requireNamespace("survival", quietly = TRUE)) stop("Need survival package.")
  if (!requireNamespace("glmnet", quietly = TRUE))   stop("Need glmnet package.")
  if (!requireNamespace("Seurat", quietly = TRUE))   stop("Need Seurat for DE.")

  modelf <- modelf[sample(1:nrow(modelf)), ]

  if (nrow(modelf) > 50) {
    modelf$set <- rep(1:10, length.out = nrow(modelf))
  } else {
    modelf$set <- 1:nrow(modelf)
  }

  # out-of-fold risk scores
  modelf$predicted_factors <- NA_real_
  modelf$predicted_factors_lasso <- NA_real_
  modelf$predicted_best_factor <- NA_real_
  modelf$predicted_best_other <- NA_real_
  modelf$predicted_pcr <- NA_real_
  modelf$predicted_pcr_lasso <- NA_real_
  modelf$predicted_de <- NA_real_
  modelf$predicted_de_lasso <- NA_real_
  predicted_genelists <- replicate(length(genelists), rep(NA_real_, nrow(modelf)), simplify = FALSE)
  names(predicted_genelists) <- names(genelists)
  predicted_arbitrary <- replicate(length(arbitrary.models), rep(NA_real_, nrow(modelf)), simplify = FALSE)
  names(predicted_arbitrary) <- names(arbitrary.models)

  y_time   <- modelf[[outcome.time.name]]
  y_status <- modelf[[outcome.status.name]]  # 1=event, 0=censored
  if (!all(y_status %in% c(0,1))) stop("Status must be coded 0/1 (censored/event).")

  best_factors <- c()
  best_others <- c()

  for (i in 1:max(modelf$set)) {

    cat(i, "\n")
    train <- subset(modelf, set != i)
    test  <- subset(modelf, set == i)

    SurvTrain <- survival::Surv(train[[outcome.time.name]], train[[outcome.status.name]])

    # 1) best single predictors by univariate Cox p-values
    cand.cols <- c(sprintf("Factor_%d", 3:20), other.predictors)
    cand.cols <- cand.cols[cand.cols %in% colnames(train)]
    get_p <- function(cl) {
      f <- as.formula(paste("survival::Surv(", outcome.time.name, ",", outcome.status.name, ") ~", cl))
      fit <- suppressWarnings(survival::coxph(f, data = train, ties = "efron"))
      as.numeric(summary(fit)$coefficients[1, "Pr(>|z|)"])
    }
    pvals <- if (length(cand.cols)) sapply(cand.cols, get_p) else numeric(0)
    p.df <- data.frame(factor_id = cand.cols, p_outcome = pvals)
    fac.only <- sprintf("Factor_%d", 3:20)
    best.factor <- with(subset(p.df, factor_id %in% fac.only),
                        as.character(factor_id[which.min(p_outcome)]))
    best_factors <- c(best_factors, best.factor)
    if (!is.null(other.predictors)) {
      best.other <- with(subset(p.df, factor_id %in% other.predictors),
                         as.character(factor_id[which.min(p_outcome)]))
      best_others <- c(best_others, best.other)
    }


    # 2) 0-shot Cox with best single predictors
    if (length(best.factor)) {
      f <- as.formula(paste("survival::Surv(", outcome.time.name, ",", outcome.status.name, ") ~", best.factor))
      mo_best_factor <- survival::coxph(f, data = train, ties = "efron")
      modelf$predicted_best_factor[modelf$set == i] <- as.numeric(predict(mo_best_factor, newdata = test, type = "lp"))
    }
    if (!is.null(other.predictors)) {
      f <- as.formula(paste("survival::Surv(", outcome.time.name, ",", outcome.status.name, ") ~", best.other))
      mo_best_other <- survival::coxph(f, data = train, ties = "efron")
      modelf$predicted_best_other[modelf$set == i] <- as.numeric(predict(mo_best_other, newdata = test, type = "lp"))
    }

    # 3) Cox with all factors + Lasso-Cox
    fac.cols <- sprintf("Factor_%d", 3:20)
    fac.cols <- fac.cols[fac.cols %in% colnames(train)]
    if (length(fac.cols)) {
      f <- as.formula(paste("survival::Surv(", outcome.time.name, ",", outcome.status.name, ") ~",
                            paste(fac.cols, collapse = "+")))
      mo_factors <- survival::coxph(f, data = train, ties = "efron")
      modelf$predicted_factors[modelf$set == i] <- as.numeric(predict(mo_factors, newdata = test, type = "lp"))

      x_tr <- as.matrix(train[, fac.cols, drop = FALSE])
      x_te <- as.matrix(test[,  fac.cols, drop = FALSE])
      mo_factors_lasso <- glmnet::cv.glmnet(x_tr, SurvTrain, family = "cox")
      modelf$predicted_factors_lasso[modelf$set == i] <-
        as.numeric(predict(mo_factors_lasso, x_te, s = use.s, type = "link"))
    }

    # 4) PCR (PCs in Cox + Lasso-Cox)
    pc.cols <- sprintf("PC_%d", 1:30)
    pc.cols <- pc.cols[pc.cols %in% colnames(train)]
    if (length(pc.cols)) {
      x_tr <- as.matrix(train[, pc.cols, drop = FALSE])
      x_te <- as.matrix(test[,  pc.cols, drop = FALSE])
      mo_pcr_lasso <- glmnet::cv.glmnet(x_tr, SurvTrain, family = "cox")
      modelf$predicted_pcr_lasso[modelf$set == i] <-
        as.numeric(predict(mo_pcr_lasso, x_te, s = use.s, type = "link"))

      pc.lm.cols <- intersect(sprintf("PC_%d", 1:18), pc.cols)
      if (length(pc.lm.cols)) {
        f <- as.formula(paste("survival::Surv(", outcome.time.name, ",", outcome.status.name, ") ~",
                              paste(pc.lm.cols, collapse = "+")))
        mo_pcr <- survival::coxph(f, data = train, ties = "efron")
        modelf$predicted_pcr[modelf$set == i] <- as.numeric(predict(mo_pcr, newdata = test, type = "lp"))
      }
    }

    # 5) DE-based features (UPDATED):
    # Use only events in the TRAIN fold; define groups by quartiles of time-to-death among events.
    events_train <- train[train[[outcome.status.name]] == 1, , drop = FALSE]
    if (nrow(events_train) >= 20) {  # need enough samples to form quartiles
      qt <- stats::quantile(events_train[[outcome.time.name]], probs = c(.25, .75), na.rm = TRUE)
      # shortest 25% time-to-death (early deaths)
      bottom <- rownames(events_train)[events_train[[outcome.time.name]] <= qt[1]]
      # longest 25% time-to-death (late deaths)
      top    <- rownames(events_train)[events_train[[outcome.time.name]] >= qt[2]]

      if (length(top) > 5 && length(bottom) > 5) {
        forde <- subset(decomp@s, cells = c(top, bottom))
        forde$group <- colnames(forde) %in% top
        Seurat::Idents(forde) <- ifelse(forde$group, "top", "bottom")

        de <- Seurat::FindMarkers(forde, ident.1 = "top", min.pcr = 0.2, logfc.threshold = 0.5)
        de <- de[order(de$p_val), , drop = FALSE]
        usegenes <- head(rownames(de), 50)

        if (length(usegenes)) {
          ma_tr <- t(as.matrix(GetAssayData(decomp@s)[usegenes, rownames(train), drop = FALSE]))
          ma_te <- t(as.matrix(GetAssayData(decomp@s)[usegenes, rownames(test),  drop = FALSE]))

          # Lasso-Cox on DE genes
          mo_de_lasso <- glmnet::cv.glmnet(ma_tr, SurvTrain, family = "cox")
          modelf$predicted_de_lasso[modelf$set == i] <-
            as.numeric(predict(mo_de_lasso, ma_te, s = use.s, type = "link"))

          # Cox on top 20 DE genes
          use20 <- head(usegenes, 20)
          tr_with <- cbind(train, ma_tr)
          te_with <- cbind(test,  ma_te)
          f <- as.formula(paste("survival::Surv(", outcome.time.name, ",", outcome.status.name, ") ~",
                                paste(sprintf("`%s`", use20), collapse = "+")))
          mo_de <- survival::coxph(f, data = tr_with, ties = "efron")
          modelf$predicted_de[modelf$set == i] <- as.numeric(predict(mo_de, newdata = te_with, type = "lp"))
        }
      }
    }

    # 6) Gene list Lasso-Cox
    if (length(genelists) > 0) {
      for (gl in seq_along(genelists)) {
        usegenes <- intersect(genelists[[gl]], rownames(decomp@s))
        if (length(usegenes) == 0) next
        ma_tr <- t(as.matrix(GetAssayData(decomp@s)[usegenes, rownames(train), drop = FALSE]))
        ma_te <- t(as.matrix(GetAssayData(decomp@s)[usegenes, rownames(test),  drop = FALSE]))
        mo_gllasso <- glmnet::cv.glmnet(as.matrix(ma_tr), SurvTrain, family = "cox")
        predicted_genelists[[gl]][modelf$set == i] <-
          as.numeric(predict(mo_gllasso, as.matrix(ma_te), s = use.s, type = "link"))
      }
    }

    #7. arbitrary models
    if (length(arbitrary.models) > 0) {
      for (mi in seq_along(arbitrary.models)) {

        f <- as.formula(paste("survival::Surv(", outcome.time.name, ",", outcome.status.name, ") ~",arbitrary.models[mi]))
        mo_ab <- survival::coxph(f, data = train, ties = "efron")
        predicted_arbitrary[[mi]][modelf$set == i] <- as.numeric(predict(mo_ab, newdata = test, type = "lp"))
      }
    }
  }

  # Evaluate with Harrell's C-index on out-of-fold predictions
  y_time_all   <- modelf[[outcome.time.name]]
  y_status_all <- modelf[[outcome.status.name]]
  cidx <- c(
    C.best.factor   = get_cindex(modelf[,"predicted_best_factor"],   y_time_all, y_status_all),
    C.best.other    = get_cindex(modelf[,"predicted_best_other"],    y_time_all, y_status_all),
    C.cox.factors   = get_cindex(modelf[,"predicted_factors"],       y_time_all, y_status_all),
    C.lasso.factors = get_cindex(modelf[,"predicted_factors_lasso"], y_time_all, y_status_all),
    C.cox.de        = get_cindex(modelf[,"predicted_de"],            y_time_all, y_status_all),
    C.lasso.de      = get_cindex(modelf[,"predicted_de_lasso"],      y_time_all, y_status_all),
    C.cox.pcr       = get_cindex(modelf[,"predicted_pcr"],           y_time_all, y_status_all),
    C.lasso.pcr     = get_cindex(modelf[,"predicted_pcr_lasso"],     y_time_all, y_status_all)
  )
  if (length(genelists) > 0) {
    cidx <- c(cidx, sapply(predicted_genelists, function(x) get_cindex(x, y_time_all, y_status_all)))
  }
  if (length(arbitrary.models) > 0) {
    cidx <- c(cidx, sapply(predicted_arbitrary, function(x) get_cindex(x, y_time_all, y_status_all)))
  }

  return(list(cidx = cidx, best.factor = best_factors, best.other = best_others))
}

