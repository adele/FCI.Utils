rm(list=ls()) # clean the environment

library(FCI.Utils)
library(pcalg)

# library(RNOmni) # useful for normal transformation


# If running in parallel:
run_parallel = TRUE

if (run_parallel) {
  require(doFuture)
  require(future.apply)
  n_cores <- parallel::detectCores() # or change to the number of cores required
  plan("multisession", workers = n_cores) # for running on RStudio
  # plan("multicore", workers = n_cores) # forking -- for running at the server
}


##############################
# Causal Model Specification #
##############################

# Specify a DAG with latent confounders
# Note: to specify a different causal model, verify
# the models in "dag_models.R"
true.admg <- getDAG(type = "discr1_nc") # "3anc"
renderDAG(true.admg$amat, add_index = F) # shows observed and latent variables

# Check the true PAG
true.amat.pag <- getTruePAG(true.admg$dagg, verbose = FALSE)@amat
renderAG(true.amat.pag, add_index = T)

# Check the true min sepsets
formatSepset(getPAGImpliedSepset(true.amat.pag))


###################
# Data Generation #
###################


# Generating a dataset according to the true ADMG
data_type = "continuous" #"mixed" #
N <- 1000
cur_seed <- 1544410977

done = FALSE
n_trials = 20
trial = 0
while (!done && trial <= n_trials) {
  cat("Genarating a dataset. trial:", trial, "\n")
  trial = trial + 1

  #cur_seed <- sample(1:.Machine$integer.max, 1)
  set.seed(cur_seed)

  if (data_type == "continuous") {
    dat_out <- generateDataset(true.admg$dagg, N, type=data_type,
                               coef_thresh = 0.2, # prevent coeff between -0.4, and 0.4.
                               b.lower = -0.9, b.upper = 0.9, # minimum and maximum coeffs
                               verbose=FALSE)
  } else if (data_type == "mixed") {
    all_vars = names(dagitty::coordinates(true.admg$dagg)[["x"]])
    f.args <- list()

    # 1: continuous; 2: binary; 3: multinominal
    n_vars <- ncol(true.admg$amat)
    var_levels <- sample(1:3, size=n_vars, replace = TRUE)
    # var_levels <- c(1,2,1,2,3,1,1) # setting manually for  c("A", "B", "C", "D", "E", "Ubd", "Ucd")

    for (vari in 1:length(all_vars)) {
      var_name <- all_vars[vari]
      f.args[[var_name]] <- list(levels = var_levels[vari])
    }

    dat_out <- generateDataset(adag = true.admg$dagg, N=N,
                               type=data_type, f.args=f.args,
                               coef_thresh = 0.2)
  }

  dat <- dat_out$dat
  if (!is.null(dat)) {
    done = TRUE
  }
}

if (is.null(dat)) {
  stop("Error generating the dataset. Try again!")
}

head(dat)



##########################
# Performing All CI Test #
##########################


alpha = 0.05
vars_names <- colnames(dat) # nodes of the graph
covs_names = c() # variables to be always conditioned on
indepTest <- mixedCITest

suffStat <- getMixedCISuffStat(dat, vars_names, covs_names)

vars_df <- dat[,vars_names, drop=FALSE]

# Compute all conditional independence tests in parallel
citestResults <- getAllCITestResults(
  vars_df, indepTest, suffStat,
  computeProbs = TRUE,        # Required if running dcFCI
  eff_size = 0.1,             # Minimum effect size threshold when computeProbs = TRUE
  fileid = cur_seed,          # Unique identifier for the dataset / analysis
  recover_citestResults = FALSE # Set TRUE to resume computation after a crash
)
save(citestResults,
     file=file.path(getwd(), paste0("citestResults_", cur_seed, ".RData")))

suffStat$citestResults <- citestResults


#############################################
# How faithful the data is to the true PAG? #
#############################################


true_sf_score_out <- dcFCI::getStraightforwardPAGScore(true.amat.pag,
                                                       suffStat$citestResults)
true_sf_score_out$score

cat("Dependencies that are not supported by the data -- i.e., with P(H1|D) < 0.5:\n")
print(true_sf_score_out$false_dep_tests)


cat("Dependence relation with the least evidence from the data:\n")
true_sf_score_out$min_dep_test

cat("Independencies that are not supported by the data -- i.e., with P(H0|D) < 0.5:\n")
print(true_sf_score_out$false_indep_tests)

cat("Independence relation with the least evidence from the data:\n")
true_sf_score_out$min_indep_test

mec_out <- dcFCI::getMEC(true.amat.pag, scored = TRUE,
                             citestResults = suffStat$citestResults)
mec_out$mec$mec_score


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

library(dcFCI)
fci_sf_score_out <- dcFCI::getStraightforwardPAGScore(fci_fit@amat, suffStat$citestResults)
fci_sf_score_out$score

fci_mec_out <- dcFCI::getMEC(fci_fit@amat, scored = TRUE,
                             citestResults = suffStat$citestResults)
fci_mec_out$mec$mec_score


#######################
# Running dcFCI from  #
# dcFCI R package     #
#######################

library(dcFCI)
library(dplyr)

#library(jsonlite)
#library(rje)
#library(pcalg)

# dcFCI is more robust under unfaithfulness
start_time <- Sys.time()
fit_dcfci <- dcFCI(suffStat, indepTest, labels, alpha,
                   m.max = Inf, fixedGaps = NULL, fixedEdges = NULL,
                   verbose = FALSE, # 1+: how verbose it can be
                   sel_top = 1,
                   run_parallel = TRUE,
                   allowNewTests=TRUE,
                   list.max = 1500, # may be increased for highly uncertain cases
                   pH0ThreshMin=0.3, # if pH0ThreshMin <= pH0 <= pH0ThreshMax, it will be a potential conditional independence
                   pH0ThreshMax = 1,
                   log_folder = file.path(getwd(), "tmp", "logs"))
end_time <- Sys.time()

time_taken <- end_time - start_time
cat("\n\n dcFCI runtime: ", time_taken, "\n\n")


top_dcfci_pag <- fit_dcfci$allPAGList[[1]]$amat.pag
top_dcfci_sepset <- fit_dcfci$allPAGList[[1]]$sepset

# renderAG(top_dcfci_pag, add_index = FALSE)
# formatSepset(top_dcfci_sepset)

dcfci_metrics <- getMetrics(true.amat.pag, top_dcfci_pag,
                            est.sepset=top_dcfci_sepset,
                            dat=NULL,  # specify dat if BIC should be computed
                            conservative = FALSE)
dcfci_metrics


##########################
# Summary of the Results #
##########################

cat("cur seed:", cur_seed)

# True PAG:
cat("\n\n ### True PAG: ### \n")
renderAG(fci_fit@amat, add_index = FALSE)
cat("Straightfoward PAG score:", fci_sf_score_out$score, "\n")
cat("MEC-Targeted PAG score:", fci_mec_out$mec$mec_score, "\n")

# FCI PAG:
cat("\n\n ### FCI PAG: ### \n")
renderAG(fci_fit@amat, add_index = FALSE)
cat("SHD FCI PAG:", shd_PAG(true.amat.pag, fci_fit@amat), "\n")
cat("FCI's Straightfoward PAG score:", fci_sf_score_out$score, "\n")
cat("FCI's MEC-Targeted PAG score:", fci_mec_out$mec$mec_score, "\n")


# dcFCI PAG:
cat("\n\n ### dcFCI's PAG: ### \n")
top_dcpag <- fit_dcfci$allPAGList[[top_pagid]]
renderAG(top_dcpag$amat.pag, add_index = F)
cat("SHD dcFCI's PAG:", shd_PAG(true.amat.pag, top_dcpag$amat.pag), "\n")
cat("dcFCI's Straightfoward PAG score:", dcfci_sf_score_out$score, "\n")
cat("dcFCI's MEC-Targeted PAG score:", dcfci_mec_out$mec$mec_score, "\n")


