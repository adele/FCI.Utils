# returns posterior probability of H0 and H1 that arise from use
# of the BF together with the assumption that H0 and H1 have
# equal prior probabilities of 1/2
#' @export BF2probs
BF2probs <- function(bf) { # Assumes bf01, i.e.  P(H0|D)/P(H1|D)
  pH0 = 1 - 1/(1 + bf) # equivalent to 1/(1+(1/bf))
  pH1 = 1/(1 +  bf)
  return(list(pH0=pH0, pH1=pH1))
}

getLRTStatistic <- function(pvalue) {
  qchisq(pvalue, df=1, lower.tail=FALSE)
}

#' @importFrom BFF chi2_test_BFF
#' @export pvalue2probs
pvalue2probs <- function(pvalue, n, eff_size=0.1) {
  if (pvalue < .Machine$double.eps) {
    pvalue <- .Machine$double.eps/1000
  }
  chiSqStat <- getLRTStatistic(pvalue)
  # gets BF_10 = P(D|H1)/P(D|H0)
  logBF_10_res <- BFF::chi2_test_BFF(chi2_stat = chiSqStat, df=1, n=n,
                                     maximize = FALSE, pearsons = FALSE, save = FALSE)
  logBF_10_ef <- max(logBF_10_res$log_BFF[which(logBF_10_res$effect_size > eff_size)])
  bf_10 <- exp(logBF_10_ef)
  if (bf_10 > 0) {
    probs <- BF2probs(1/bf_10)
  } else {
    probs <- list(pH0=1, pH1=0)
  }
  return(list(pH0=probs$pH0, pH1=probs$pH1))
}
