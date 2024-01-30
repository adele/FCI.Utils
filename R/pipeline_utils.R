#' @importFrom pcalg fci
#' @export runFCIHelper
runFCIHelper <- function(indepTest, suffStat, alpha = 0.05,
                         citestResults = NULL,
                         labels=NULL, conservative=FALSE, m.max=Inf,
                         savePlots=TRUE, add_index=FALSE, saveFiles=TRUE,
                         fileid=NULL, file_type="png",
                         output_folder="./temp/") {
  if (is.null(labels)) {
    labels <- colnames(samples)
  }

  p <- length(labels)
  fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  fixedGaps = NULL
  NAdelete = FALSE
  verbose = FALSE

  # run original FCI
  fit_fci <- pcalg::fci(suffStat, indepTest = indepTest,
                     skel.method = "stable", labels = labels, m.max=m.max,
                     NAdelete = NAdelete, type = "normal", alpha = alpha,
                     verbose = verbose, conservative = conservative)

  if(savePlots) {
    renderAG(fit_fci@amat, output_folder, fileid = fileid, type = file_type,
             labels=labels, add_index = add_index)
  }

  fci_pag <- fit_fci@amat
  fci_sepset <- fixSepsetList(fit_fci@sepset)

  ci_dist <- impliedCondIndepDistance(amat.pag = fci_pag,
                                      indepTest, suffStat, alpha=alpha, verbose=TRUE)
  violations <- hasViolation(fci_pag, fci_sepset, conservative=conservative,
                             log=TRUE, verbose=TRUE)

  fci_out <- list(pag=fci_pag, sepset=fci_sepset,
              ci_dist=ci_dist, violations=violations)

  if (saveFiles) {
    save(fci_out, file=paste0(output_folder, "fci_out_", fileid, ".RData"))
  }

  return(fci_out)
}
