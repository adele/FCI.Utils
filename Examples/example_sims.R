rm(list=ls())

library(FCI.Utils)
library(pcalg)

type = "chain4"  # "discr2_nc"
adag_out <- getDAG(type=type)
true.amat.pag <- getTruePAG(adag_out$dagg)@amat
labels <- colnames(true.amat.pag)

output_folder <- paste0("../Results/", type, "/")
renderAG(true.amat.pag, output_folder, fileid = "truePAG", type = "png",
         add_index = FALSE)

true.sepset <- getPAGImpliedSepset(true.amat.pag)
formatSepset(true.sepset)


edgeTypesList <- NULL
N = 10000 # sample size
type = "binary"  # "continuous"
nsims = 10

metrics <- data.frame()

for (sim in 1:nsims) {
  # Generating the dataset with variables as columns and observations as rows
  adat_out <- generateDataset(adag = adag_out$dagg, N=N, type=type)
  dat <- adat_out$dat
  head(dat)
  labels <- colnames(dat)

  # Setting up the Conditional Independence Test
  if (type == "continuous") {
    # Assuming Gaussian
    indepTest <- gaussCItest
    suffStat <- list(C = cor(dat), n = nrow(dat))
    suffStat$test_type <- "frequentist"
  } else {
    # Assuming Binary
    indepTest <- binCItest
    suffStat <- list(dm=dat, adaptDF=TRUE)
  }


  alpha = 0.05

  fci_out <- runFCIHelper(indepTest, suffStat, alpha=alpha,
                          labels=labels, fileid=paste0("sim_", sim),
                          output_folder=output_folder)
  metrics <- rbind.data.frame(metrics, c(fci_out$ci_dist$dist, fci_out$violations$out))
  edgeTypesList <- getEdgeTypesList(fci_out$pag, edgeTypesList = edgeTypesList)
}
colnames(metrics) <- c("ci_dist", "violations")

edgeTypesSumm <- summarizeEdgeTypesList(edgeTypesList)
edgeTypesSumm
