# Generates obs. dataset following a linear SEM, compatible with a dagitty DAG, adag
# type: defines the type of the variables, "continuous", "binary", or "mixed"
# f.args: a list indexed by the names of the variables, where each entry is another list
# with an entry named levels indicating the number of levels a discrete node takes or
# levels = 1 for continuous variable
#' @importFrom dagitty simulateLogistic simulateSEM localTests edges latents
#' @importFrom simmixdag
#' @export generateDataset
generateDataset <- function(adag, N, type="continuous", verbose=FALSE,
                            f.args = NULL, coef_thresh = 0.2,
                            b.lower = -0.8, b.upper = 0.8,
                            b.default = NULL) {
  if (!(type %in% c("continuous", "binary", "mixed")))  {
    stop("type must be continuous, binary, or mixed")
  }

  if (type == "mixed" && is.null(f.args)) {
    stop("f.args must be specified.")
  }

  lt <- NA
  done <- FALSE
  ntries <- 0
  obs.dat <- NULL
  while (!done && ntries <= 100) {
    done <- tryCatch(
      {
        if (is.null(b.default)) {
          if (b.lower < 0 && b.upper > 0) {
            b.default <- sapply(1:nrow(dagitty::edges(adag)), function(x) {
              if (runif(1) > 0.5) runif(1, b.lower, -coef_thresh)
              else runif(1, coef_thresh, b.upper)})
          } else {
            b.default <- sapply(1:nrow(dagitty::edges(adag)), function(x) {
              runif(1, b.lower, b.upper)})
          }
        }

        if(type == "binary") {
          obs.dat <- dagitty::simulateLogistic(adag, b.default = b.default, N=N, verbose=verbose)
          obs.dat <- as.data.frame(sapply(obs.dat, function(col) as.numeric(col)-1))
          lt <- dagitty::localTests(adag, obs.dat, type="cis.chisq")
          TRUE
        } else if (type == "mixed") {
          require(simMixedDAG)

          min_coef = 0
          ntries2 <- 0
          while (min_coef < coef_thresh & ntries2 <= 100) {
            param_dag_model <- parametric_dag_model(dag = adag, f.args = f.args)
            min_coef <- min(abs(unlist(lapply(param_dag_model$f.args, function(x) {x$betas} ))))
            ntries2 <- ntries2 + 1
          }

          obs.dat <- sim_mixed_dag(param_dag_model, N=N)
          lat_vars <- dagitty::latents(adag)
          lat_cols <- which(colnames(obs.dat) %in% lat_vars)
          if (length(lat_cols) > 0) {
            obs.dat <- obs.dat[, -lat_cols]
          }
          TRUE
        } else if (type == "continuous") {
          obs.dat <- dagitty::simulateSEM(adag, b.default = b.default, N=N)
          lt <- dagitty::localTests(adag, obs.dat, type="cis")
          R <- cor(obs.dat)
          valR <- matrixcalc::is.symmetric.matrix(R) &&
            matrixcalc::is.positive.definite(R, tol=1e-8)
          valR
        }
      }, error=function(cond) {
        message(paste0("ntries: ", ntries, " - ", cond))
        return(FALSE)
      })
    ntries = ntries + 1
  }
  return(list(dat=obs.dat, lt=lt))
}

#' @importFrom dagitty canonicalize
#' @export generateDatasetFromPAG
generateDatasetFromPAG <- function(apag, N, type="continuous", f.args = NULL,
                                   coef_thresh = 0.2, verbose=FALSE) {
  adag <- dagitty::canonicalize(getMAG(apag)$magg)$g
  return(generateDataset(adag, N=N, type=type, verbose=verbose, f.args = f.args,
                         coef_thresh=coef_thresh))
}


isAncestralGraph <- function(amat.mag) {
  ret <- tryCatch({
    if (!is.null(amat.mag)) {
      ug_mag <- (amat.mag == 3 & t(amat.mag == 3)) * 1
      bg_mag <- (amat.mag == 2 & t(amat.mag == 2)) * 1
      dg_mag <- (amat.mag == 2 & t(amat.mag == 3)) * 1
      mag_ggm <- ggm::makeMG(dg_mag, ug_mag, bg_mag)
      retAG <- ggm::isAG(mag_ggm)
      retAG
    } else {
      FALSE
    }
  },
  error=function(cond) {
    print(cond)
    return(FALSE)
  },
  warning=function(cond) {
    print(cond)
    return(FALSE)
  })
  return(ret)
}

#' @export getRandomMAG
getRandomMAG <- function(n_nodes, dir_edges_prob = 0.4, bidir_edges_prob = 0.2) {
  done = FALSE
  while(!done) {
    amat.mag <- matrix(0, nrow = n_nodes, ncol=n_nodes)
    colnames(amat.mag) <- rownames(amat.mag) <- LETTERS[seq( from = 1, to = n_nodes)]

    edges <- combn(1:n_nodes, 2)
    n_edges <- dim(edges)[2]
    dir_edges <- sample(1:n_edges, floor(n_edges * dir_edges_prob), replace = FALSE)
    for (i in dir_edges) {
      amat.mag[edges[1,i], edges[2,i]] <- 2
      amat.mag[edges[2,i], edges[1,i]] <- 3
    }

    bidir_edges <- sample((1:n_edges)[-dir_edges], floor(n_edges * bidir_edges_prob), replace = FALSE)
    for (i in bidir_edges) {
      amat.mag[edges[1,i], edges[2,i]] <- 2
      amat.mag[edges[2,i], edges[1,i]] <- 2
    }

    if (isAncestralGraph(amat.mag)) {
      done = TRUE
    }
  }
  return(amat.mag)
}


#' @export generateUniqueRandomPAGs
generateUniqueRandomPAGs <- function(n_graphs = 10, n_nodes = 5,
                                     dir_edges_prob = 0.2, bidir_edges_prob = 0.3,
                                     verbose=FALSE) {
  truePAGs <- list()
  stats <- c()

  while (length(truePAGs) < n_graphs) {
    amat.mag <- getRandomMAG(n_nodes, dir_edges_prob = dir_edges_prob, bidir_edges_prob = bidir_edges_prob)
    labels <- colnames(amat.mag)
    #renderAG(amat.mag)
    mec <- MAGtoMEC(amat.mag, verbose=verbose)

    if (length(which(mec$CK$ord >= 1)) > 0 || length(which(mec$NCK$ord >= 1)) > 0) {
      if (verbose) {
        cat("PAG", length(truePAGs), "with nCK1:", length(which(mec$CK$ord >= 1)),
            "and nNCK1", length(which(mec$NCK$ord >= 1)), "\n")
      }

      cur_stats <- c(nCK1 = length(which(mec$CK$ord >= 1)),
                     nNCK = length(which(mec$NCK$ord >= 1)))
      stats <- rbind(stats, cur_stats)

      amag <- pcalg::pcalg2dagitty(amat.mag, colnames(amat.mag), type="mag")
      truePAG <- getTruePAG(amag)
      amat.pag <- truePAG@amat
      #renderAG(amat.pag)
      truePAGs[[length(truePAGs) + 1]] <- amat.pag

      dupl_ids <- which(duplicated(truePAGs))
      if (length(dupl_ids) > 0) {
        truePAGs <- truePAGs[-dupl_ids]
        stats <- stats[-dupl_ids, ]
      }
    }
  }

  return(list(pags=truePAGs, stats=stats))
}


#' @export getFaithfulnessDegree
getFaithfulnessDegree <- function(amat.pag, citestResults,
                                  cutoff=0.5, alpha=0.01,
                                  bayesian=TRUE, verbose=FALSE) {
  labels <- colnames(amat.pag)
  # exp_indep <- data.frame()
  # exp_dep <- data.frame()

  f_citestResults <- c()
  for (i in 1:nrow(citestResults)) {
    cur_row <- citestResults[i, , drop=TRUE]
    snames <- labels[getSepVector(cur_row$S)]
    xname <- labels[cur_row$X]
    yname <- labels[cur_row$Y]

    def_msep <- isMSeparated(amat.pag, xname, yname, snames,
                             verbose=verbose)
    if (def_msep) {
      if (bayesian) {
        ret <- c(cur_row, type="indep",
               bf = cur_row$pH0 > cutoff,
               pf = cur_row$pvalue > alpha)
      } else {
        ret <- c(cur_row, type="indep",
                 pf = cur_row$pvalue > alpha)
      }
      f_citestResults <- rbind.data.frame(f_citestResults, ret)
    } else {
      if (bayesian) {
        ret <- c(cur_row, type="dep",
                 bf = cur_row$pH1 > cutoff,
                 pf = cur_row$pvalue <= alpha)
      } else {
        ret <- c(cur_row, type="dep",
                 pf = cur_row$pvalue <= alpha)
      }
      f_citestResults <- rbind.data.frame(f_citestResults, ret)
    }
  }

  faithful_pprop = length(which(f_citestResults$pf)) / length(f_citestResults$pf)

  if (bayesian) {
    probs = c(subset(f_citestResults, type == "indep", select = pH0, drop=TRUE),
                       subset(f_citestResults, type == "dep", select = pH1, drop=TRUE))
    faithful_bprop = length(which(f_citestResults$bf)) / length(f_citestResults$bf)
    ret <- list(probs = probs,
                f_citestResults = f_citestResults,
                faithful_bprop = faithful_bprop,
                faithful_pprop = faithful_pprop)
  } else {
    ret <- list(f_citestResults = f_citestResults,
                faithful_pprop = faithful_pprop)

  }
  return(ret)
}

#' @export getAdjFaithfulnessDegree
#' @export getAdjFaithfulnessDegree
getAdjFaithfulnessDegree <- function(apag, citestResults, test_dep_min = FALSE,
                                     ord=Inf, alpha=0.05) {
  if (!isValidPAG(apag)) {
    return(NULL)
  }

  cur_ord <- ord
  if (is.infinite(ord) || ord > ncol(apag) - 2) {
    cur_ord <- getSepsetMaxOrd(getPAGImpliedSepset(apag))
  }

  sepset <- getPAGImpliedSepset(apag)
  #print(formatSepset(sepset))

  tested_ci <- data.frame()
  for (i in 1:(ncol(apag)-1)) {
    for (j in (i+1):ncol(apag)) {
      resultsxy <- subset(citestResults, X == i & Y == j)
      resultsyx <- subset(citestResults, X == j & Y == i)
      cur_citests <- rbind(resultsxy, resultsyx)
      Sxy_list <- sepset[[i]][[j]]
      if (is.null(Sxy_list)) {
        # there is an edge
        relation="e"

        # the probability of the conjunction of all c.i. tests indicating dependency
        potsepsets_citests <- subset(citestResults, X==i & Y==j & ord <= cur_ord)
        tested_ci <- rbind(tested_ci, cbind(potsepsets_citests, type = "dep", structure="edge"))
      } else if (test_dep_min) {
        # there is one or more minimal separating sets
        for (Sxy in Sxy_list) {
          # adding the p-value associated to the minimal independence
          tested_ci <- rbind(tested_ci,
                             cbind(subset(citestResults, X==i &  Y==j & S==getSepString(Sxy)),
                                   type = "indep", structure="separ"))
          psubsets_citests <- cur_citests[ sapply(cur_citests$S, isProperSubset,
                                                  setStr=Sxy, all_subsets=TRUE), ]
          if (nrow(psubsets_citests) > 0) {
            tested_ci <- rbind(tested_ci,
                               cbind(psubsets_citests, type = "dep", structure="min"))
          }
        }
      }
    }
  }
  tested_ci <- tested_ci[order(tested_ci$ord), ]

  faithful_pprop <- NA
  if ("pvalue" %in% colnames(tested_ci)) {
    cindep_ids <- which(tested_ci$type == "indep")
    tested_ci[cindep_ids, "p_faithful"] <- tested_ci[cindep_ids, "pvalue"] > alpha
    cdep_ids <- which(tested_ci$type == "dep")
    tested_ci[cdep_ids, "p_faithful"] <- tested_ci[cdep_ids, "pvalue"] <= alpha
    faithful_pprop <- length(which(tested_ci$p_faithful))/nrow(tested_ci)
  }

  faithful_bprop <- NA
  if ("pH0" %in% colnames(tested_ci) && "pH1" %in% colnames(tested_ci)) {
    tested_ci[cindep_ids, "b_faithful"] <- tested_ci[cindep_ids, "pH0"] > 0.5

    cdep_ids <- which(tested_ci$type == "dep")
    tested_ci[cdep_ids, "b_faithful"] <- tested_ci[cdep_ids, "pH1"] > 0.5

    faithful_bprop <- length(which(tested_ci$b_faithful))/nrow(tested_ci)
  }

  ret <- list(tested_ci=tested_ci, sepset=sepset)
  if (!is.na(faithful_pprop)) {
    ret <- c(ret, faithful_pprop=faithful_pprop)
  }

  if (!is.na(faithful_bprop)) {
    ret <- c(ret, faithful_bprop=faithful_bprop)
  }

  return(ret)
}

