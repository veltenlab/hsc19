setwd("/users/lvelten/lvelten/Analysis/Aim2/hsc18/data")
alejo_lists <- as.list(read.csv("/users/lvelten/lvelten/Analysis/Aim2/CRISPRi_decomposition/alejo_genelists.csv",skip=1))
alejo_lists <- alejo_lists[1:17]
alejo_lists <- lapply(alejo_lists, function(x) x[x!=""])
factor_annotation <- read.csv("~/cluster2/lvelten/Analysis/Aim2/CRISPRi_joe_gsfa/Factor_annotation.csv")
save(alejo_lists, file="alejo_lists.rda")
save(factor_annotation, file="factor_annotation.rda")

#path to the output from GSFA
#97 genes
#fit0_guide_prtb <- readRDS("/users/lvelten/jbowness/CRISPRi500/HSC/clean_analysis/test_analysis/GSFA_inputs/fit0_guide_prtb.rds")
#500 genes
fit0_guide_prtb <- readRDS("/users/lvelten/jbowness/CRISPRi500/HSC/clean_analysis/test_analysis/GSFA_inputs/all_prtb_fit4.rds")
gibbs_PM_guide_prtb <- fit0_guide_prtb$posterior_means
lfsr_mat_guide_prtb <- fit0_guide_prtb$lfsr[, -ncol(fit0_guide_prtb$lfsr)]
perturbation_guide_names <- colnames(lfsr_mat_guide_prtb)
# load("~/fastdata/GSFA_components.rda")

PIP_mat <- gibbs_PM_guide_prtb$F_pm #these are apparently the PIPs for gene loadings
W_mat <- gibbs_PM_guide_prtb$W_pm #but the absolute value matters?
W_mat_sparse <- W_mat
W_mat_sparse[PIP_mat < 0.8] <- 0
hsc18_loadings <- W_mat_sparse[rowSums(W_mat_sparse) !=0,c(1:16,18,20)] #make sure to put this to [3:20] with the old data
save(W_mat, file="W_mat.rda")
save(PIP_mat, file="PIP_mat.rda")
save(hsc18_loadings, file="hsc18_loadings.rda")

gamma_mat <- gibbs_PM_guide_prtb$Gamma_pm #these are apparently the PIPs for gene loadings
beta_mat <- gibbs_PM_guide_prtb$beta_pm #but the absolute value matters?
beta_mat_sparse <- beta_mat
beta_mat_sparse[gamma_mat < 0.5] <- 0
hsc18_perturbation_beta <- beta_mat_sparse
save(beta_mat, file="beta_mat.rda")
save(gamma_mat, file="gamma_mat.rda")
save(hsc18_perturbation_beta, file="hsc18_perturbation_beta.rda")

