rm(list=ls())

library(FCI.Utils)
library(pcalg)
library(jsonlite)

library(doFuture)
n_cores <- 8
plan("multicore", workers = n_cores, gc=TRUE)


type = "2discrs_nc" #"discr1_c"
adag_out <- getDAG(type)
truePAG <- getTruePAG(adag_out$dagg)
trueAdjM <- truePAG@amat
labels <- colnames(trueAdjM)
renderAG(trueAdjM)

output_folder <- paste0("../Results/", type, "/")
renderAG(trueAdjM, output_folder, fileid = "truePAG", type = "png",
         add_index = FALSE)


N = 50000 # sample size
type = "mixed" # "continuous" #"binary" #

# # Generating the dataset with variables as columns and observations as rows
aseed <- 59803094 # sample(1:.Machine$integer.max, 1)
set.seed(aseed)


f.args <- NULL
if (type == "mixed") {
  f.args <- list()
  if (length(labels) == 4) {
    var_levels <- c(1, 3, 2, 1)
  } else if (length(labels) == 5) {
    var_levels <- c(1, 1, 1, 3, 2)
  }

  for (vari in 1:length(labels)) {
    var_name <- labels[vari]
    f.args[[var_name]] <- list(levels = var_levels[vari])
  }
}

adat_out <- generateDataset(adag = adag_out$dagg, N=N, type=type, f.args = f.args)
dat <- adat_out$dat
head(dat)
summary(dat)

dat$A <- cut(dat$A, breaks = 2, ordered_result = TRUE)
str(dat)


alpha = 0.05
all_vars <- colnames(dat)


##############################
# Running using mixedCITests #
##############################

vars_names = all_vars
vars_df <- dat[,vars_names, drop=FALSE]
covs_names = c()

indepTest <- mixedCITest
suffStat <- getMixedCISuffStat(dat = dat,
                               vars_names = vars_names,
                               covs_names = covs_names)

###############################################
# Testing some conditional independence tests #
###############################################

# ordinal and continuous:
indepTest(1,3,NULL, suffStat)
indepTest(1,3,2, suffStat)

# ordinal and binary
indepTest(1,5,NULL, suffStat)
indepTest(1,5,2, suffStat)

# ordinal and multinominal
indepTest(1,4,NULL, suffStat)
indepTest(1,4,2, suffStat)
indepTest(1,4,c(2,3), suffStat)


################################################
# Computing all conditional independence tests #
################################################

fileid <- paste0("seed_", aseed)
citestResults <- getAllCITestResults( vars_df, indepTest, suffStat,
                                      m.max=Inf, computeProbs = TRUE,
                                      fileid=fileid,
                                      citestResults_folder=output_folder)

f_citestResults <- getFaithfulnessDegree(amat.pag = trueAdjM,
                                         citestResults = citestResults)$f_citestResults
subset(f_citestResults, bf == FALSE | pf == FALSE)

citestResults_file <- paste0(output_folder, fileid, "_citestResults.RData")
save(citestResults, file=citestResults_file)
suffStat$citestResults <- citestResults


#########################
# Running FCI Algorithm #
#########################

fileid <- paste0("mixedCI_", paste0(vars_names, collapse="_"))
fci_out <- runFCIHelper(indepTest, suffStat, alpha=alpha,
                        labels=vars_names, fileid=fileid,
                        output_folder=output_folder)

fci_out$violations$out
formatSepset(fci_out$sepset)
renderAG(fci_out$pag)


############################
# Running over all subsets #
############################

edgeTypesList <- NULL
metrics <- data.frame()
for (n in 3:length(vars_names)) {
  subsets <- as.matrix(combn(vars_names, n))
  for (j in 1:ncol(subsets)) {
    cur_var_names <- subsets[,j]
    suffStat <- getMixedCISuffStat(dat = dat,
                                   vars_names = cur_var_names,
                                   covs_names = covs_names)
    suffStat$citestResults <- extractValidCITestResults(citestResults,
                                                        vars_names, cur_var_names)

    fileid <- paste0(cur_var_names, collapse="_")
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
