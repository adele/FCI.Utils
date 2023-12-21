# returns the number of times an implied conditional independence is not
# observed in the data.
#' @export impliedCondIndepDistance
impliedCondIndepDistance <- function(amat.pag, indepTest, suffStat, alpha,
                                     verbose=FALSE) {
  dist_indep <- 0 # number of incorrect implied independence relations
  false_indeps <- c()

  skel <- amat.pag > 0
  amat.mag <- getMAG(amat.pag)$magg
  sepset <- impliedSepsets(amat.mag)

  labels <- colnames(skel)
  eids <- which(!skel & upper.tri(skel), arr.ind = TRUE) # missing edges
  eid <- 1
  while (eid <= nrow(eids)) {
    x <- eids[eid,1]
    y <- eids[eid,2]
    xname <- labels[x]
    yname <- labels[y]
    Sxy_list <- sepset[[xname]][[yname]]
    for (Sxy in Sxy_list) {
      # Sxy is a minimal separator of X and Y in the MAG
      if (!is.null(Sxy) && length(Sxy) == 0) {
        indep_p <- indepTest(x, y, NULL, suffStat)
      } else {
        Sxyids <- which(labels %in% Sxy)
        indep_p <- indepTest(x, y, Sxyids, suffStat)
      }
      if (indep_p <= alpha) { # rejects HO of independence
        if (verbose) {
          Sxynames <- if (is.null(Sxy) || length(Sxy) == 0) "" else paste0(labels[Sxyids], collapse=",")
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


# TODO rewrite this explanation
# TODO Breakdown and list different types of violations
#
# This function checks violations in a PAG due to:
#  1) violations in the properties of MEC of ancestral graphs, which can be checked
#     by verifying if the canonical MAG, supposedly member of the PAG, is indeed
#     an ancestral graph, and if the PAG is the same as the one corresponding to
#     to the canonical MAG.
#  2) preservation of m-connections observed in previous steps:
#     a) if a dependence was observed before, then it must exist a m-connecting
#        path between the pair
#     b) if has be found a non-empty sepset Sxy for X and Y, then:
#        - no definite m-connecting path must exist between X and Y given Sxy; and
#        - for every Si \in Sxy, there must exist a definite m-connecting path
#        - between X and Y given Sxy\Si such that Si in lying on it.
#' @importFrom ggm makeMG isAG
#' @export hasViolation
hasViolation <- function(pagAdjM, sepset, alpha=NULL, citestResults=NULL,
                         listall=TRUE, log=FALSE, verbose=FALSE) {

  logList <- list()
  labels <- colnames(pagAdjM)

  ##############################################################
  # 1) Violations in the properties of MEC of ancestral graphs #
  ##############################################################

  # Here we check whether the PAG is valid by checking whether
  # the canonical MAG is ancestral.
  validPAG <- isValidPAG(pagAdjM, verbose=verbose)
  logList["validPAG"] <- validPAG

  if (!validPAG) {
    if (log) {
      return(list(out=TRUE, log=logList))
    } else {
      return(TRUE)
    }
  }

  ###############################################################
  # 2) Preservation of m-connections observed in previous steps #
  ###############################################################

  obsDepend <- NULL
  if (!is.null(citestResults) && !is.null(alpha)) {
    obsDepend = citestResults[which(citestResults$pH0 < alpha),]
  }


  # a) Checking dependencies in citestResults
  #############################################

  # Verifying if the dependencies once observed
  # are still represented in pagAdjM through m-connecting paths
  violates <- FALSE
  if (!is.null(obsDepend) && nrow(obsDepend) > 0) {
    logList[["def_m-connections"]] <- data.frame()
    for (i in 1:nrow(obsDepend)) {
      x <- obsDepend[i,"X"]
      y <- obsDepend[i,"Y"]
      S <- getSepVector(obsDepend[i,"S"])

      xname <- labels[x]
      yname <- labels[y]
      snames <- paste0(labels[S], collapse={","})

      if (verbose) {
        cat(paste("Checking whether there is a definite connecting path \n",
                  "representing the observed dependence: \n",
                  xname, "and", yname,
                  " should be m-connected given S={",
                  snames,
                  "}\n"))
      }
      def_msep <- isMSeparated(pagAdjM, xname, yname, labels[S],
                            verbose=verbose)
      logList[["def_m-connections"]] <- rbind(logList[["def_m-connections"]],
                                          c(xname, yname, snames, !def_msep))

      if (def_msep) {
        violates <- TRUE
        if (verbose) {
          cat(paste("    --> violation!\n"))
        }
        if (!listall) {
          if (log) {
            return(list(out=TRUE, log=logList))
          } else {
            return(TRUE)
          }
        }
      } else {
        if (verbose) {
          cat(paste("    --> OK!\n"))
        }
      }
    }
  }

  # b) Checking dependencies implied by minimal separating sets
  ##############################################################

  # if no separating set has been found yet, then there is no violations
  # regarding implied m-connecting paths have to be checked.
  if (length(sepset) == 0) {
    if (log) {
      return(list(out=FALSE, log=logList))
    } else {
      return(FALSE)
    }
  }

  checkSepsets <- c()
  for (i in 1:length(sepset)) {
    tocheck <- sapply(sepset[[i]], function(x) {
      !(is.null(x)) # || length(x) == 0 || (length(x) == 1 && x == ""))
    })
    tocheck <- which(tocheck)
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
    Sij <- getSepVector(sepset[[vi]][[vj]])

    xname <- labels[vi]
    yname <- labels[vj]
    snames <- paste0(labels[Sij], collapse=",")

    if (verbose) {
      cat(paste("Checking if {", paste0(labels[Sij], collapse={","}),
                "} m-separates", labels[vi], "and", labels[vj],"\n"))
    }
    def_msep <- isMSeparated(pagAdjM, xname, yname, labels[Sij],
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

    # Since the independence is assumed to be minimal, the implied/observed
    # dependencies (i.e., for every S \in Sij, V_i \indep V_j | Sij \ Si) has to be
    # represented by a definite m-connecting path in the PAG containing S.
    sepmin_out <- checkSepMinimality(pagAdjM, vi, vj, Sij, listall, log=log, verbose=verbose)
    if (log) {
      violates <-  violates || sepmin_out$out
      logList[["def_m-connections_min"]] <-
        rbind.data.frame(logList[["def_m-connections_min"]], sepmin_out$log)
    } else {
      violates <- violates || sepmin_out
    }

    ssid = ssid + 1
  }

  if (log) {
    return(list(out=violates, log=logList))
  } else {
    return(violates)
  }
}

checkSepMinimality <- function(pagAdjM, vi, vj, Sij, listall,
                               log=FALSE, verbose=FALSE) {

  logdf <- data.frame()

  labels <- colnames(pagAdjM)

  violates <- FALSE
  if (length(Sij) > 0) { # (!("" %in% Sij)) {
    for (vs in Sij) {
      if (verbose)
        cat(paste("-> Checking if ", labels[vs] , "in", paste0(labels[Sij], collapse={","}),
                  "is necessary to m-separate", labels[vi], "and", labels[vj],"\n"))

      varsVminusVs <- setdiff(Sij, vs)
      labelsVminusVs <- c()
      if (length(varsVminusVs) > 0) {
        labelsVminusVs <- labels[varsVminusVs]
      }

      # Faithful separation is evaluated through the existence of a definite
      # m-connecting path between X-Y given S/{S_i} that contains Si, for every S_i in S

      xname <- labels[vi]
      yname <- labels[vj]
      snames <- paste0(labels[varsVminusVs], collapse=",")
      vname <- labels[vs] # connection should be through this path
      curlog <- c(xname, yname, snames, vname)

      connpaths <- getMConnPaths(pagAdjM, labels[vi], labels[vj], labelsVminusVs,
                                 definite=TRUE, verbose=verbose)

      if (is.null(connpaths) || length(connpaths) == 0) {
        violates <- TRUE
        curlog <- c(curlog, FALSE)
        if (verbose) {
          cat(paste0("    --> violation: there is no definite m-connecting path between ",
                     labels[vi], " and ", labels[vj], " given ", labelsVminusVs, "\n"))
        }
      } else if (!any(sapply(connpaths, function(x) { (labels[vs] %in% x) } ))) {
        violates <- TRUE
        curlog <- c(curlog, FALSE)
        if (verbose) {
          cat(paste0("    --> violation: ", labels[vs],
                     " is not in a definite m-connecting path between ",
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
    colnames(logdf) <- c("x", "y", "S", "v", "mconn")
  }

  if (log) {
    return(list(out=violates, log=logdf))
  } else {
    return(violates)
  }
}

