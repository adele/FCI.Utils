#####################################
# A few metrics to compare two PAGs #
#####################################

# The adjacency matrix of a PAG is encoded as following:
# [i,j] encodes the edgemark of vj in the edge vi-jj
# 0: no edge
# 1: circle
# 2: arrowhead
# 3: tail
# e.g., [i,j] = 2 and [j,i] = 1 means vi o-> vj


# returns a matrix in which [i,j]
#  = 0 implies that i is a definite non-ancestor of j
#  = 1 implies that i is a definite ancestor of j
#  = 2 implies that i is a possible ancestor of j
#' @export getAncestralMatrix
getAncestralMatrix <- function(amat.pag) {
  defTrueAncM <- getAncestralMatrixHelper(amat.pag, definite = T)
  possTrueAncM <- getAncestralMatrixHelper(amat.pag, definite = F)
  return((possTrueAncM - defTrueAncM) + possTrueAncM)
}

# returns a matrix in which [i,j] = 1 implies that
# i is a possible/definite ancestor of j
getAncestralMatrixHelper <- function(amat.pag, definite=TRUE) {
  ancM <- 0 * amat.pag
  if (definite) {
    # Removing every edge with a circle
    to_rm <- which(amat.pag == 1, arr.ind = T)
    to_rm2 <- to_rm[,c(2,1)]
    amat.pag[rbind(to_rm, to_rm2)] <- 0
  }
  labels <- colnames(ancM)
  for (vj in 1:length(labels)) {
    possAncVj <- pcalg::possAn(amat.pag, x=vj, type="pag", ds=FALSE)
    possAncVj <- setdiff(possAncVj, vj)
    if (length(possAncVj) > 0) {
      ancM[possAncVj, vj] <- 1
    }
  }
  return(ancM)
}



#' @export getPAGPosNegMetrics
getPAGPosNegMetrics <- function(amat.trueP, amat.estP) {
  trueAncM <- getAncestralMatrix(amat.trueP)
  estAncM <- getAncestralMatrix(amat.estP)
  diag(trueAncM) = NA
  diag(estAncM) = NA

  trueAdjM <- (amat.trueP != 0) * 1
  estAdjM <- (amat.estP != 0) * 1
  diag(trueAdjM) = NA
  diag(estAdjM) = NA

  # 0 mean non-edge in the true adjacency matrix
  est_non_adj <- which(estAdjM[lower.tri(estAdjM)] == 0) #
  exp_non_adj <- trueAdjM[lower.tri(trueAdjM)][est_non_adj]
  # estimated non-adjacencies that are indeed present in the true PAG
  true_non_adj <- length(which(exp_non_adj == 0))
  # estimated non-adjacencies that are not present in the true PAG
  false_non_adj <- length(which(exp_non_adj != 0))

  # 1 means definite ancestral relation
  est_def_anc <- which(estAncM == 1 & estAdjM == 1, arr.ind=T) # def. ancestral and adjacent
  exp_def_anc <- trueAncM[est_def_anc]
  # estimated def. ancestral relations that are indeed present in the true PAG
  true_def_anc <- length(which(exp_def_anc == 1))
  # estimated def. ancestral relations that are not present in the true PAG
  false_def_anc <- length(which(exp_def_anc != 1))

  # 0 means definite non-ancestral relation
  est_def_nanc <- which(estAncM == 0 & estAdjM == 1, arr.ind=T) # def. non-ancestral and adjacent
  exp_def_nanc <- trueAncM[est_def_nanc]
  # estimated def. non-ancestral relations that are indeed present in the true PAG
  true_def_nanc <- length(which(exp_def_nanc == 0))
  # estimated def. non-ancestral relations that are not present in the true PAG
  false_def_nanc <- length(which(exp_def_nanc != 0))


  # 2 means undetermined (circle)
  est_undet <- which(estAncM == 2 & estAdjM == 1, arr.ind=T)
  exp_undet <- trueAncM[est_undet]
  # estimated non-invariance (circle) relations that are indeed present in the true PAG
  true_undet <- length(which(exp_undet == 2))
  # estimated non-invariance (circle) relations that are invariances in the true PAG
  false_undet <- length(which(exp_undet != 2))

  total_pos <- false_def_anc + false_def_nanc + false_non_adj + true_def_anc + true_def_nanc + true_non_adj
  false_discovery_rate <- if (total_pos == 0) 0 else
    (false_def_anc + false_def_nanc + false_non_adj) / total_pos
  total_neg <- true_undet + false_undet
  false_omission_rate <- if (total_neg == 0) 0 else false_undet / total_neg

  ret <- list(false_discovery_rate=false_discovery_rate,
              false_omission_rate=false_omission_rate,
              true_def_anc=true_def_anc, false_def_anc=false_def_anc,
              true_def_nanc=true_def_nanc, false_def_nanc=false_def_nanc,
              true_undet=true_undet, false_undet=false_undet)
  ret
  return(ret)
}

#' @export skelDistance
skelDistance <- function(skel.trueP, skel.estP, verbose=FALSE) {
  skell.diff <- 0
  if (!is.null(skel.trueP) && !is.null(skel.estP)) {
    if (any(dim(skel.trueP) != dim(skel.estP))) {
      stop("amat.trueP and amat.estP must have same dimensions.")
    }
    skell.diff <- length(which(skel.trueP - skel.estP != 0))
  }
  return(skell.diff)
}

# Note: The SHD is not the best metric, as circles in
# definite non-colliders have a different meaning than colliders.
# For example, if A o-> B <-o C is the true PAG, we cannot say that an
# estimated PAG A o-o B o-o C has no mistakes. It would be only
# weaker if A o-o B o-o C, with A o-o C.
#' @export shd_PAG
shd_PAG <- function(amat.trueP, amat.estP, verbose=FALSE) {
  strDist <- 0
  amat.estP <- amat.estP[rownames(amat.trueP), colnames(amat.trueP)]
  if (!is.null(amat.trueP)) {
    if (any(dim(amat.estP) != dim(amat.trueP))) {
      stop("amat.trueP and amat.estP must have same dimensions.")
    }
    strDist <- length(which(amat.estP - amat.trueP != 0))
  }
  return(strDist)
}
