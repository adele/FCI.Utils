# xname = labels[x]
# yname = labels[y]
# snames <- NULL
# TODO change to impliedConditionalSepset
getMinimalSeparator <- function(pagAdjM, xname, yname, snames, definite=TRUE,
                                ignore_path_list = list(),
                                ignore_sepvar_names = list(),
                                verbose = FALSE) {
  Z <- snames
  updated = TRUE
  while (updated) {
    defconnpaths <- getMConnPaths(pagAdjM, xname, yname, Z, definite = definite)
    ndcpaths <- length(defconnpaths)
    rmids <- c()
    if (ndcpaths > 0) {
      for (i in 1:ndcpaths) {
        dcpath <- defconnpaths[[i]]
        if (sapply(ignore_path_list, function(x) {
          length(x) == length(dcpath) && all(x %in% dcpath)} )) {
          rmids <- c(rmids, i)
        }
      }
    }
    pathList <- defconnpaths
    if (!is.null(rmids)) {
      pathList <- defconnpaths[-rmids]
    }
    i = 1
    updated = FALSE
    for (curpath in pathList) {
      # there are still paths that are connecting and of definite status
      len_curpath <- length(curpath)

      if (len_curpath > 2) {
        for (k in 1:(len_curpath-2)) {
          vi_name <- curpath[k]
          vm_name <- curpath[k+1]
          vj_name <- curpath[k+2]

          triplet <- c(vi_name, vm_name, vj_name)
          vnames <- colnames(pagAdjM)
          if (isCollider(pagAdjM, triplet)) {
            vm_id <- which(vnames == triplet[2])
            devm <- vnames[pcalg::searchAM(pagAdjM, vm_id, type="de")]
            if (length(which(devm %in% Z)) == 0) {
              # no definite descendant of this collider is in Z, so the path is
              # blocked and we don't need to checked it again util Z is changed...
              #pathList <- pathList[-i]
              break
            } else {
              next
            }
          } else if (isDefiniteNonCollider(pagAdjM, triplet)) {
            if ((vm_name %in% ignore_sepvar_names)) {
              next
            } else {
              if (length(which(Z == vm_name)) == 0) {
                # we add the non-collider to Z
                Z <- c(Z, vm_name)
                updated = TRUE # since we added a new variable to Z, all paths have to be rechecked.
              }
              # since the path is definitely blocked, we go to the next path in pathList
              break
            }
          }
        }
      }
    }
  }
  return(list(sepset=Z, connpaths = pathList))
}

getIndepStats <- function(x, y, sepset, labels, citestResults, indepTest,
                          suffStat, alpha, NAdelete, verbose=FALSE) {
  xname <- labels[x]
  yname <- labels[y]

  # computing indep pvalue
  Sxy_ids_list <- list()
  cur_ps <- c()
  Sxy_list <- sepset[[xname]][[yname]]

  for (Sxy in Sxy_list) {
    if (!is.null(Sxy) && length(Sxy) == 0) {
      # this only holds when z is a collider
      Sxyids <- NULL
      Sxy_ids_list[length(Sxy_ids_list)+1] <- list(Sxyids)
    } else {
      Sxyids <- which(labels %in% Sxy)
      Sxy_ids_list[[length(Sxy_ids_list)+1]] <- Sxyids
    }

    citest_out <- getCITestResults(x, y, Sxyids, citestResults, indepTest,
                                   suffStat, NAdelete, verbose)
    citestResults <- citest_out$citestResults
    cur_ps <- c(cur_ps, citest_out$pH0)
  }
  if (verbose && length(cur_ps) > 1) {
    cat("There are more than one minimal separator...\n!")
  }

  indep_score <- mean(sapply(cur_ps, getIndepScore, alpha))

  return(list(pH0_list=cur_ps, Sxy_ids_list=Sxy_ids_list,
              indep_score=indep_score, citestResults=citestResults))
}


