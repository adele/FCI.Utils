rm(list=ls()) # clean the environment

library(FCI.Utils)
library(pcalg)

# library(RNOmni) # useful for normal transformation

# In case you want to run in parallel
library(doFuture)
library(future.apply)
n_cores <- 12 # change to the number of cores available
#plan("multisession", workers = n_cores)
plan("multicore", workers = n_cores) # forking


##############################
# Causal Model Specification #
##############################

# Specify a DAG with latent confounders
# Note: to specify a different causal model, verify
# the models in "dag_models.R"
true.admg <- getDAG(type = "discr1_nc")
renderDAG(true.admg$amat, add_index = FALSE) # shows observed and latent variables

# Check the true PAG
true.amat.pag <- getTruePAG(true.admg$dagg, verbose = FALSE)@amat
renderAG(true.amat.pag, add_index = F)

# Generate a dataset that agrees to the causal DAG
N = 300 # sample size

# it can be "continuous", "binary", or "mixed";
type = "continuous"

# Change the seed to evaluate different datasets
# seed <- sample(1000000, 1)
# set.seed(seed)

# set.seed(861524) # seed for a faithful dataset with 300 samples
set.seed(152244) # seed for an unfaithful dataset with 300 samples

# for binary and continuous data, we can set coef_thresh, b.lower, and b.upper
# for mixed data, it is required to set f.args
dat_out <- generateDataset(true.admg$dagg, N, type,
                           coef_thresh = 0.4, # prevent coeff between -0.4, and 0.4.
                           b.lower = -0.8, b.upper = 0.8, # minimum and maximum coeffs
                           verbose=FALSE)
dat <- dat_out$dat
if (is.null(dat)) {
  stop("Error generating the dataset. Try again!")
}

head(dat)

##########################
# Assess faithfulness to #
# the true causal model  #
# (model testing)        #
##########################

alpha = 0.05
vars_names <- colnames(dat) # nodes of the graph
covs_names = c() # variables to be always conditioned on
indepTest <- mixedCITest

suffStat <- getMixedCISuffStat(dat, vars_names, covs_names)

vars_df <- dat[,vars_names, drop=FALSE]

# Computes all conditional independence tests in parallel
citestResults <- getAllCITestResults(vars_df, indepTest, suffStat)

suffStat$citestResults <- citestResults # so we don´t redo tests later

# Checking faithfulness degree
faithf_degree <- getFaithfulnessDegree(true.amat.pag, citestResults,
                                       bayesian = FALSE, alpha = alpha)
unfaith_tests <- subset(faithf_degree$f_citestResults, pf == FALSE)

if (nrow(unfaith_tests) > 0) {
  cat("Data is unfaithful to the model!\n")
  print(unfaith_tests)
} else {
  cat("Data is faithful to the model!")
}

####################
# Causal Discovery #
####################

###########################
# Running the FCI from    #
# pcalg R package         #
###########################

# Under faithfulness, FCI recovers the true PAG
labels <- colnames(suffStat$dataset)
fci_fit <- fci(suffStat, indepTest, alpha, labels)
renderAG(fci_fit@amat, add_index = FALSE)


faithf_fci <- getFaithfulnessDegree(fci_fit@amat, suffStat$citestResults,
                                    bayesian = FALSE)
faithf_fci

#######################
# Running dcFCI from  #
# dcFCI R package     #
#######################

library(dcFCI)
library(dplyr)

#library(jsonlite)
#library(rje)
#library(pcalg)

# Computes BF for all conditional independence tests
# eff_size needs to be adjusted with the min value to be considered significant
citestResults <- getAllCITestResults(vars_df, indepTest, suffStat,
                                    computeProbs = TRUE,
                                    eff_size = 0.1)

suffStat$citestResults <- citestResults # saving tests

# dcFCI is more robust under unfaithfulness
fit_dcfci <- dcFCI(suffStat, indepTest,
                   labels=labels, alpha=alpha,
                   sel_top = 1,
                   prob_sel_top = 1,
                   verbose = 1,
                   run_parallel = TRUE)


top_dcpag <- fit_dcfci$allPAGList[[1]]
renderAG(top_dcpag$amat.pag, add_index = F)

faithf_dcfci <- getFaithfulnessDegree(top_dcpag$amat.pag, citestResults)
faithf_dcfci

fit_dcfci$scoresDF
fit_dcfci$mec_score_df


##########################
# Summary of the Results #
##########################

cat("Number of unfaithful tests:", nrow(unfaith_tests), "\n")
if (nrow(unfaith_tests) > 0) {
  print(unfaith_tests)
}

# FCI's PAG:
cat("\n\n ### FCI's PAG: ### \n")
renderAG(fci_fit@amat, add_index = FALSE)
cat("SHD FCI's PAG:", shd_PAG(true.amat.pag, fci_fit@amat), "\n")
cat("Compatibility FCI's PAG:", faithf_fci$faithful_pprop, "\n")

# dcFCI's PAG:
cat("\n\n ### dcFCI's PAG: ### \n")
renderAG(top_dcpag$amat.pag, add_index = F)
cat("SHD dcFCI's PAG:", shd_PAG(true.amat.pag, top_dcpag$amat.pag), "\n")
cat("Compatibility dcFCI's PAG:", faithf_dcfci$faithful_pprop, "\n")
cat("Compatibility dcFCI's PAG (Bayesian):", faithf_dcfci$faithful_bprop, "\n")


