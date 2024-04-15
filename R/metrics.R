# returns the number of times an implied conditional independence is not
# observed in the data.
#' @export impliedCondIndepDistance
impliedCondIndepDistance <- function(amat.pag, indepTest, suffStat, alpha,
                                     verbose=FALSE) {
  dist_indep <- 0 # number of incorrect implied independence relations
  false_indeps <- c()

  skel <- amat.pag > 0
  sepset <- getPAGImpliedSepset(amat.pag)
  if (any(is.na(sepset))) {
    return(NA)
  }

  labels <- colnames(skel)
  eids <- which(!skel & upper.tri(skel), arr.ind = TRUE) # missing edges
  eid <- 1
  while (eid <= nrow(eids)) {
    x <- eids[eid,1]
    y <- eids[eid,2]
    xname <- labels[x]
    yname <- labels[y]
    #Sxy_list <- sepset[[xname]][[yname]]
    Sxy_list <- sepset[[x]][[y]]
    if (!is.list(Sxy_list)) {
      Sxy_list <- list(Sxy_list)
    }
    for (Sxy in Sxy_list) {
      # Sxy is a minimal separator of X and Y in the MAG
      if (!is.null(Sxy) && length(Sxy) == 0) {
        indep_p <- indepTest(x, y, NULL, suffStat)
      } else {
        Sxyids <- Sxy # which(labels %in% Sxy)
        indep_p <- indepTest(x, y, Sxyids, suffStat)
      }
      if (indep_p <= alpha) { # rejects HO of independence
        Sxynames <- if (is.null(Sxy) || length(Sxy) == 0) "" else paste0(labels[Sxyids], collapse=",")
        if (verbose) {
          cat("False implied cond. indep. of ", xname, "and", yname,
              "given {", Sxynames, "} -- p-value = ", indep_p, "\n")
        }
        false_indeps <- rbind.data.frame(false_indeps, c(xname, yname, Sxynames, indep_p))
        dist_indep <- dist_indep + 1
      }
    }
    eid <- eid + 1
  }
  if (!is.null(false_indeps)) {
    colnames(false_indeps) <- c("x", "y", "S", "pvalue")
  }

  return(list(dist=dist_indep, false_indeps=false_indeps))
}


# This function checks violations in a PAG due to:
#  1) violations in the properties of MEC of ancestral graphs, which can be checked
#     by verifying if the canonical MAG, supposedly member of the PAG, is indeed
#     an ancestral graph, and if the PAG is the same as the one corresponding to
#     to the canonical MAG.
#  2) violations in the minimality of the sepsets:
#      if has be found a non-empty sepset Sxy for X and Y, then:
#        - no definite m-connecting path must exist between X and Y given Sxy; and
#        - for every S' \subset Sxy, there must exist a definite m-connecting path
#        - between X and Y given S'such that Si \in S' lies on it.
#' @param amat.pag: the adjacency matrix of a PAG.
#' @param sepset: a given list of separators, usually the one returned by the
#' FCI algorithm,  which may differ from the minimal separators implied by amat.pag.
#' @importFrom ggm makeMG isAG
#' @export hasViolation
hasViolation <- function(amat.pag, sepset, listall=TRUE, conservative=FALSE,
                         knowledge=FALSE, log=FALSE, verbose=FALSE) {

  logList <- list()
  labels <- colnames(amat.pag)
  violates <- FALSE

  ####################################################################
  # Checking Violations in the properties of MEC of ancestral graphs #
  ####################################################################

  # Here we check whether the PAG is valid by checking whether
  # the canonical MAG is ancestral.
  validPAG <- isValidPAG(amat.pag, conservative = conservative,
                         knowledge = knowledge, verbose=verbose)
  logList["validPAG"] <- validPAG

  if (!validPAG) {
    violates <- TRUE
    if (verbose) {
      cat("PAG is invalid.\n")
      cat(paste("   --> violation!\n"))
    }

    # Returns if PAG is not valid, as
    # next checks rely on  sepsetDistance and isMSeparated,
    # which are only defined for valid PAGs
    if (log) {
      return(list(out=violates, log=logList))
    } else {
      return(violates)
    }
  }


  #######################################################################
  # Checking whether the list of separators in sepset is different from #
  # the list of minimal separators implied by amat.pag                  #
  #######################################################################

  sepsetDist <- sepsetDistance(amat.pag, sepset, verbose=verbose)
  logList["sepsetDist"] <- sepsetDist

  if (sepsetDist != 0) {
    if (verbose) {
      cat("sepset does not match with list of implied minimal separators.\n")
      cat(paste("   --> violation!\n"))
    }
    violates <- TRUE
    if (!listall) {
      if (log) {
        return(list(out=violates, log=logList))
      } else {
        return(violates)
      }
    }
  }


  #####################################################################
  # Representation of the dependencies associated to separated pairs, #
  # when conditioned on proper subsets of their minimal separators.   #
  #####################################################################

  # # if no separating set has been found yet, then there is no violations
  # # regarding implied m-connecting paths have to be checked.
  # if (length(sepset) == 0) {
  #   if (log) {
  #     return(list(out=FALSE, log=logList))
  #   } else {
  #     return(FALSE)
  #   }
  # }

  checkSepsets <- c()
  for (i in 1:length(sepset)) {
    tocheck <- which(sapply(sepset[[i]], function(x) {!(is.null(x))}))
    if (length(tocheck) > 0) {
      checkSepsets <- rbind(checkSepsets, cbind(i, tocheck))
    }
  }

  # removing duplicated pairs
  checkSepsets <- as.data.frame(checkSepsets)
  i = 1
  while (i <= nrow(checkSepsets)) {
    checkSepsets[i,] <- sort(unlist(checkSepsets[i,]))
    i = i+1
  }
  checkSepsets <- checkSepsets[!duplicated(checkSepsets),]

  ssid <- 1
  logList[["m-separations"]] <- data.frame()
  logList[["def_m-connections_min"]] <- data.frame()

  # iterates oves all pairs of variables (Vi,Vj) with a sepset Sij such that |Sij| > 1
  while ((listall || !violates) && ssid <= dim(checkSepsets)[1]) {
    vi <- checkSepsets[ssid,2]
    vj <- checkSepsets[ssid,1]
    Sij_list <- sepset[[vi]][[vj]]
    if (!is.list(Sij_list)) {
      # This is a hack to support sepsets whose entries are not list
      Sij_list <- list(Sij_list)
    }
    for (Sij in Sij_list) {
      Sij <- getSepVector(Sij)

      xname <- labels[vi]
      yname <- labels[vj]
      snames <- paste0(labels[Sij], collapse=",")

      if (verbose) {
        cat(paste("Checking if {", paste0(labels[Sij], collapse={","}),
                  "} m-separates", labels[vi], "and", labels[vj],"\n"))
      }
      def_msep <- isMSeparated(amat.pag, xname, yname, labels[Sij],
                   verbose=verbose)
      logList[["m-separations"]] <- rbind.data.frame(logList[["m-separations"]],
                                          c(xname, yname, snames, def_msep))

      # The observed independence (i.e., V_i \indep V_j | Sij) has to be
      # represented by the corresponding m-separation in the PAG
      if (!def_msep) {
        violates <- TRUE
        if (verbose) {
          cat(paste("    --> violation!\n"))
        }
      } else {
        if (verbose) {
          cat(paste("    --> OK!\n"))
        }
      }

      if (length(logList[["m-separations"]]) > 0) {
        colnames(logList[["m-separations"]]) <- c("x", "y", "S", "msep")
      }

      # Since the independence given Sij is assumed to be minimal, each
      # dependence given a proper subset S' of Sij has to be
      # represented by a definite m-connecting path in the PAG containing some
      # variable in Sij \ S'
      sepmin_out <- checkSepMinimality(amat.pag, vi, vj, Sij, listall,
                                       log=log, verbose=verbose)
      if (log) {
        violates <-  violates || sepmin_out$out
        logList[["def_m-connections_min"]] <-
          rbind.data.frame(logList[["def_m-connections_min"]], sepmin_out$log)
      } else {
        violates <- violates || sepmin_out
      }
    }

    ssid = ssid + 1
  }

  if (log) {
    return(list(out=violates, log=logList))
  } else {
    return(violates)
  }
}

checkSepMinimality <- function(amat.pag, vi, vj, Sij, listall,
                               log=FALSE, verbose=FALSE) {

  logdf <- data.frame()
  labels <- colnames(amat.pag)

  violates <- FALSE
  if (length(Sij) > 0) {
    properSubsets <- getSubsets(Sij, TRUE)
    for (vs in properSubsets) {
      if (verbose) {
        cat(paste("-> Checking if there is a definite m-connecting path between",
                labels[vi], "and", labels[vj], "given {", getSepString(labels[vs]), "} \n"))
      }

      varsVminusVs <- setdiff(Sij, vs)
      labelsVminusVs <- c()
      if (length(varsVminusVs) > 0) {
        labelsVminusVs <- labels[varsVminusVs]
      }

      xname <- labels[vi]
      yname <- labels[vj]
      snames <- getSepString(labels[vs])
      curlog <- c(xname, yname, snames)

      connpaths <- getMConnPaths(amat.pag, labels[vi], labels[vj], snames,
                                 definite=TRUE, verbose=verbose)

      if (is.null(connpaths) || length(connpaths) == 0) {
        violates <- TRUE
        curlog <- c(curlog, FALSE)
        if (verbose) {
          cat(paste0("    --> violation: there is no definite m-connecting path between ",
                     labels[vi], " and ", labels[vj], " given {", snames,   "}\n"))
        }
      } else if (!any(sapply(connpaths, function(x) { any(labels[varsVminusVs] %in% x) } ))) {
          violates <- TRUE
        curlog <- c(curlog, FALSE)
        if (verbose) {
          cat(paste0("    --> violation: none of {",
                     paste0(labels[varsVminusVs], collapse = ","),
                     "} lies on a definite m-connecting path between ",
                     labels[vi], " and ", labels[vj], "\n"))
          print(connpaths)
        }
      } else {
        curlog <- c(curlog, TRUE)
        if (verbose) {
          cat(paste("    --> OK!\n"))
        }
      }
      logdf <- rbind.data.frame(logdf, curlog)
      if (violates && !listall) {
        break
      }
    }
  }

  if (length(logdf) > 0) {
    colnames(logdf) <- c("x", "y", "S", "mconn") #"v",
  }

  if (log) {
    return(list(out=violates, log=logdf))
  } else {
    return(violates)
  }
}

# Returns the number of differences between a given sepset
# and the implied minimal sepset
#' @importFrom utils combn
#' @export sepsetDistance
sepsetDistance <- function(amat.pag, sepset, verbose=FALSE) {
  distSepset <- 0
  impliedSepset <- getPAGImpliedSepset(amat.pag)
  if (any(is.na(impliedSepset))) {
    return(NA)
  }
  labels <- colnames(amat.pag)
  p <- length(labels)
  pairs <- combn(1:p,2)
  for (i in 1:ncol(pairs)) {
    v1.ind <- pairs[1, i]
    v2.ind <- pairs[2, i]
    S12_list <- impliedSepset[[v1.ind]][[v2.ind]]
    estS12 <- sepset[[v1.ind]][[v2.ind]]
    if (is.list(estS12)) {
      estS12 <- estS12[[1]]
    }
    if (is.null(S12_list) && is.null(estS12)) {
      matched = TRUE
    } else {
      matched = FALSE
      for (impliedS12 in S12_list) {
        if (length(estS12) == length(impliedS12) && all(estS12 %in% impliedS12)) {
          matched = TRUE
        }
      }
    }
    if (!matched) {
      distSepset <- distSepset + 1
      if (verbose) {
        cat("Estimated sepset for", labels[v1.ind], "and", labels[v2.ind], "is",
            paste0(labels[estS12], collapse=","),
            "but implied minimal sepsets are:",
            paste0(lapply(S12_list, function(x) {
              paste0("{", paste0(labels[x], collapse=","), "}") }), collapse = "; "), ".\n")
      }
    }
  }
  return(distSepset)
}

