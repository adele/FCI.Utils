rm(list=ls())

library(FCI.Utils)
library(pcalg)

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

# Setting up the Conditional Independence Test
if (type == "continuous") {
  # Assuming Gaussian
  indepTest <- gaussCItest
} else {
  # Assuming Binary
  indepTest <- binCItest
}

alpha = 0.05
all_vars <- colnames(dat)


edgeTypesList <- NULL
metrics <- data.frame()
for (n in 3:length(all_vars)) {
  subsets <- combn(all_vars, n)
  for (j in 1:ncol(subsets)) {
    labels <- subsets[,j]
    cur_dat <- dat[,labels]

    if (type == "continuous") {
      # Assuming Gaussian
      suffStat <- list(C = cor(cur_dat), n = nrow(cur_dat))
    } else {
      # Assuming Binary
      suffStat <- list(dm=cur_dat, adaptDF=TRUE)
    }

    fileid <- paste0(labels, collapse="_")
    fci_out <- runFCIHelper(indepTest, suffStat, alpha=alpha,
                            labels=labels, fileid=fileid,
                            output_folder=output_folder)
    metrics <- rbind.data.frame(metrics, c(fileid, fci_out$ci_dist$dist, fci_out$violations$out))
    edgeTypesList <- getEdgeTypesList(fci_out$pag, edgeTypesList = edgeTypesList)
  }
}

colnames(metrics) <- c("fileid", "ci_dist", "violations")
metrics

edgeTypesSumm <- summarizeEdgeTypesList(edgeTypesList)
edgeTypesSumm
