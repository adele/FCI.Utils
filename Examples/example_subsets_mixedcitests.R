rm(list=ls())

library(FCI.Utils)
library(pcalg)

library(doParallel)
registerDoParallel(cores=7)

type =  "discr2_nc"
adag_out <- getDAG(type)
truePAG <- getTruePAG(adag_out$dagg)
trueAdjM <- truePAG@amat
labels <- colnames(trueAdjM)

output_folder <- paste0("../Results/", type, "/")
renderAG(trueAdjM, output_folder, fileid = "truePAG", type = "png",
         add_index = FALSE)


N = 10000 # sample size
type = "binary"  # "continuous"

# Generating the dataset with variables as columns and observations as rows
adat_out <- generateDataset(adag = adag_out$dagg, N=N, type=type)
dat <- adat_out$dat
head(dat)

alpha = 0.05
all_vars <- colnames(dat)
covs_names = c()

citestResults <- runAllMixedCITests(dat, vars_names=all_vars,
                                    covs_names=covs_names, alpha=alpha)

indepTest <- mixedCITest
suffStat <- getMixedCISuffStat(dat = dat,
                               vars_names = all_vars,
                               covs_names = covs_names)
suffStat$citestResults <- citestResults


edgeTypesList <- NULL
metrics <- data.frame()
for (n in 3:length(all_vars)) {
  subsets <- combn(all_vars, n)
  for (j in 1:ncol(subsets)) {
    cur_var_names <- subsets[,j]
    cur_dat <- dat[,cur_var_names]

    suffStat <- getMixedCISuffStat(dat = cur_dat,
                                   vars_names = cur_var_names,
                                   covs_names = covs_names)
    suffStat$citestResults <- extractValidCITestResults(citestResults,
                                                        all_vars, cur_var_names)

    fileid <- paste0(labels, collapse="_")
    fci_out <- runFCIHelper(indepTest, suffStat, alpha=alpha,
                            labels=cur_var_names, fileid=fileid,
                            output_folder=output_folder)
    metrics <- rbind.data.frame(metrics, c(fileid, fci_out$ci_dist$dist,
                                           fci_out$violations$out))
    edgeTypesList <- getEdgeTypesList(fci_out$pag, edgeTypesList = edgeTypesList)
  }
}

colnames(metrics) <- c("fileid", "ci_dist", "violations")
metrics

edgeTypesSumm <- summarizeEdgeTypesList(edgeTypesList)
edgeTypesSumm
