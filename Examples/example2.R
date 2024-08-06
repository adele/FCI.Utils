rm(list=ls())
library(FCI.Utils)
library(dcFCI)
library(jsonlite)
library(lavaan)

library(doFuture)
n_cores <- 8
plan("multicore", workers = n_cores, gc=TRUE)


#g_type = "iv2"
#g_type =  "discr1_nc"

g_type =  "2discrs_nc"
g_type = "3anc"
g_type = "2ivs"
g_type =  "collfork2"

adag_out <- getDAG(g_type)
true.amat.pag <- getTruePAG(adag_out$dagg, verbose = TRUE)@amat
renderAG(true.amat.pag)

labels <- colnames(true.amat.pag)

true.sepset <- getPAGImpliedSepset(true.amat.pag)
formatSepset(true.sepset)

trueMEC <- getMEC(true.amat.pag, ag.type="pag", scored=FALSE)
trueMEC


output_folder <- paste0("../Results/", g_type, "/")
if (!file.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}
renderAG(true.amat.pag, output_folder, fileid = "truePAG", type = "png",
         add_index = TRUE)


N = 10000 # sample size
data_type = "continuous"

while(1) {
  adat_out <- generateDataset(adag = adag_out$dagg, N=N, type=data_type)
  dat <- adat_out$dat
  head(dat)
  labels <- colnames(dat)
  write.csv(dat, file=paste0("./", "dat.csv"), row.names = FALSE)

  indepTest <- mixedCITest
  alpha = 0.05
  vars_names <- labels
  covs_names = c()
  suffStat <- getMixedCISuffStat(dat = dat,
                                 vars_names = vars_names,
                                 covs_names = covs_names)

  fci_out <- runFCIHelper(indepTest, suffStat, alpha=alpha,
                          citestResults = citestResults,
                          labels=labels, fileid="fci_test",
                          output_folder="./",
                          conservative=FALSE)
  fci_posneg <- getPAGPosNegMetrics(true.amat.pag, fci_out$pag)
  fci_out$violations
  formatSepset(fci_out$sepset)

  if (fci_posneg$false_def_anc > 0) {
    dcfci_out <- runDCFCIHelper(indepTest, suffStat, alpha=alpha,
                                labels=labels, fileid="dcfci_test",
                                output_folder="./",
                                verbose=FALSE)
    dcfci_posneg <- getPAGPosNegMetrics(true.amat.pag, dcfci_out$pag)
    if (dcfci_posneg$false_discovery_rate < fci_posneg$false_discovery_rate &&
        dcfci_posneg$false_omission_rate < fci_posneg$false_omission_rate &&
        dcfci_posneg$false_omission_rate < 0.5)
      break
  }
}


renderAG(fci_out$pag)
renderAG(dcfci_out$pag)

cfci_out <- runFCIHelper(indepTest, suffStat, alpha=alpha,
                         citestResults = citestResults,
                         labels=labels, fileid="cfci_test",
                         output_folder="./",
                         conservative=TRUE)


#gets the positions of the true and FCI's PAG in the dcFCI's PAG list
true_pag_id <- getSepsetIndex(dcfci_out$fit_dcfci$best_pag_List, true.sepset)
fci_pag_id <- getSepsetIndex(dcfci_out$fit_dcfci$best_pag_List, fci_out$sepset)

formatSepset(true.sepset)
formatSepset(fci_out$sepset)
formatSepset(dcfci_out$sepset)

lapply(dcfci_out$fit_dcfci$best_pag_List, function(x) {formatSepset(x$sepset)})

dcfci_out$fit_dcfci$best_scores_df
renderAG(dcfci_out$fit_dcfci$best_pag_List[[2]]$amat.pag)
scores <- scorePAG(dcfci_out$fit_dcfci$best_pag_List[[1]]$amat.pag,
         citestResults = dcfci_out$fit_dcfci$citestResults)







