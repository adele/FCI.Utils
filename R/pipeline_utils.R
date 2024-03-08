#' @importFrom pcalg fci
#' @export runFCIHelper
runFCIHelper <- function(indepTest, suffStat, alpha = 0.05,
                         citestResults = NULL, labels=NULL,
                         conservative=FALSE, maj.rule=FALSE, m.max=Inf,
                         fixedEdges=NULL, fixedGaps = NULL,
                         savePlots=TRUE, add_index=FALSE, saveFiles=TRUE,
                         fileid=NULL, file_type="png",
                         output_folder="./temp/") {
  p <- length(labels)
  if (is.null(fixedEdges)) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  fixedGaps = NULL
  NAdelete = FALSE
  verbose = FALSE

  # run original FCI
  fit_fci <- pcalg::fci(suffStat, indepTest = indepTest,
                     skel.method = "stable", labels = labels, m.max=m.max,
                     NAdelete = NAdelete, type = "normal", alpha = alpha,
                     verbose = verbose, conservative = conservative,
                     maj.rule = maj.rule)

  if(savePlots) {
    renderAG(fit_fci@amat, output_folder, fileid = fileid, type = file_type,
             labels=labels, add_index = add_index)
  }

  fci_pag <- fit_fci@amat
  fci_sepset <- fixSepsetList(fit_fci@sepset)

  ci_dist <- impliedCondIndepDistance(amat.pag = fci_pag,
                                      indepTest, suffStat, alpha=alpha, verbose=TRUE)
  violations <- hasViolation(fci_pag, fci_sepset, conservative=conservative,
                             knowledge = FALSE, log=TRUE, verbose=TRUE)

  fci_out <- list(pag=fci_pag, sepset=fci_sepset,
              ci_dist=ci_dist, violations=violations)

  if (saveFiles) {
    if (conservative == F && maj.rule == F) {
      cur_fileid <- paste0("fci_out_", fileid)
      save(fci_out, file=paste0(output_folder, cur_fileid, ".RData"))
    } else if (conservative == T && maj.rule == F) {
      cur_fileid <- paste0("cfci_out_", fileid)
      cfci_out <- fci_out
      save(cfci_out, file=paste0(output_folder, cur_fileid, ".RData"))
    } else if (conservative == F && maj.rule == T) {
      cur_fileid <- paste0("mjrfci_out", fileid)
      mjrfci_out <- fci_out
      save(mjrfci_out, file=paste0(output_folder, cur_fileid, ".RData"))
    }
  }

  return(fci_out)
}
