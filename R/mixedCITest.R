modListErrors <- function(modList) {
  for (mod in modList) {
    errMod <- tryCatch(
    {
      summary(mod)
      1
    }, error=function(cond) {
      NA
    }, warning=function(cond) {
      NA
    })
    if (is.na(errMod)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

isBinary <- function(x) {
  return(length(unique(x)) == 2 && all(names(table(x)) %in% c(0,1)))
}

#' @importFrom stats complete.cases
#' @export zicoSeqCITest
zicoSeqCITest <- function(x, y, S, suffStat) {
  labels <- colnames(suffStat$dataset)
  feature.dat <- t(suffStat$fulltaxa_df)
  meta.dat <- suffStat$dataset[, c(x, S), drop = FALSE]
  compl_ids <- complete.cases(meta.dat)
  meta.dat <- meta.dat[compl_ids, , drop=FALSE]
  feature.dat <- feature.dat[ ,compl_ids, drop=FALSE]

  xname <- labels[x]
  yname <- labels[y]
  snames <- labels[S]

  zicoseq_out <- GUniFrac::ZicoSeq(feature.dat = feature.dat,
                         meta.dat = meta.dat,
                         grp.name = xname,
                         adj.name = snames,
                         feature.dat.type = "count",
                         return.feature.dat = FALSE,
                         prev.filter = 0,
                         is.winsor = FALSE,
                         is.post.sample = TRUE,
                         outlier.pct = 0,
                         mean.abund.filter = 0,
                         max.abund.filter = 0,
                         perm.no = 1500,
                         excl.pct = 0.2,
                         stage.no = 6,
                         post.sample.no = 25,
                         ref.pct = 0.5)
  p1 <- zicoseq_out$p.raw[yname]
  return(list(p=p1, mod0=NULL, mod1=zicoseq_out))
}

#' @export lindaCITest
lindaCITest <- function(x, y, S, suffStat) {
  feature.dat <- t(suffStat$fulltaxa_df)
  meta.dat <- suffStat$dataset[, c(x, S), drop = FALSE]
  labels <- colnames(suffStat$dataset)
  xname <- labels[x]
  yname <- labels[y]

  formula <- getRHSFormulaStr(x, y, S, suffStat)

  linda_out <- MicrobiomeStat::linda(feature.dat,
                     meta.dat,
                     feature.dat.type = "count",
                     formula = formula,
                     zero.handling = "imputation",
                     alpha = 0.05,
                     prev.filter = 0,
                     adaptive = TRUE,
                     #is.winsor = FALSE,
                     outlier.pct = 0,
                     mean.abund.filter = 0)
  p1 <- linda_out$output[[xname]][yname,"pvalue"]
  return(list(p=p1, mod0=NULL, mod1=linda_out))
}


simpleZeroInflNegBinCITest <- function (x, y, S, suffStat) {
  formulae <- getFitFormulae(x, y, S, suffStat, randStr=NULL)
  if (length(suffStat$rand_varnames) == 0) {
    formulae_rand <- getFitFormulae(x, y, S, suffStat, randStr="1")
  } else {
    formulae_rand <- getFitFormulae(x, y, S, suffStat,
                                    randStr= paste0(suffStat$rand_varnames,
                                                    collapse = " + "))
  }

  p1 <- NA
  mod1 <- NULL
  if (sum(suffStat$dataset[,y]==0) >= 1) {
    tryCatch({
      mod1 <- pscl::zeroinfl(formulae_rand$formula_YXS,
                             data=suffStat$dataset, dist = "negbin")
      xids <- which(grepl(colnames(suffStat$dataset)[x],
                          rownames(summary(mod1)$coefficients$count)))
      p1 <- min(summary(mod1)$coefficients$count[xids,4])
    }, error = function(e){})
  } else{
    tryCatch({
      mod1 <- MASS::glm.nb(formulae$formula_YXS, data=suffStat$dataset)
      xids <- which(grepl(colnames(suffStat$dataset)[x],
                          rownames(summary(mod1)$coefficients)))
      p1 <- summary(mod1)$coefficients[xids,4]
    }, error = function(e){})
  }

  return(list(p=p1, mod0=NULL, mod1=mod1))
}


#' @importFrom stats glm binomial complete.cases
#' @importFrom stats anova pchisq
#' @export logisticCITest
logisticCITest <- function (x, y, S, suffStat) {
  ydat = suffStat$dataset[, y]
  xdat = suffStat$dataset[, x]

  mod0 <- NULL
  if (is.null(S) || length(S) == 0) {
    data <- data.frame(x = xdat)

    mod1 <- stats::glm(ydat ~ ., data = data, stats::binomial)
    t1 <- mod1$null.deviance - mod1$deviance
    dof1 <- length(mod1$coefficients) - 1
    p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
  } else {
    ds0 <- suffStat$dataset[, S, drop = FALSE]
    ds1 <- suffStat$dataset[, c(S, x)]
    compl_ids <- which(complete.cases(ds1))
    ydat <- ydat[compl_ids]
    ds1 <- ds1[compl_ids, , drop=FALSE]
    ds0 <- ds0[compl_ids, , drop=FALSE]

    mod1 <- stats::glm(ydat ~ ., data = ds1, stats::binomial)
    mod0 <- stats::glm(ydat ~ ., data = ds0, stats::binomial)
    a1 <- stats::anova(mod0, mod1)
    t1 <- a1[2, 4]
    dof1 <- a1[2, 3]
    p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
  }
  return(list(p=p1, mod0=mod0, mod1=mod1))
}

# returns the right hand side of a formula
getRHSFormulaStr <- function(x, y, S, suffStat) {
  x_name <- colnames(suffStat$dataset)[x]
  formula_XS_str <- paste("~ 1 + ", x_name)

  if (length(S) > 0) {
    S_names <- colnames(suffStat$dataset)[S]
    formula_XS_str <- paste(formula_XS_str, "+", paste(S_names, collapse= " + "))
  }

  return(formula_XS_str)
}

getFitFormulae <- function(x, y, S, suffStat, randStr=NULL) {
  x_name <- colnames(suffStat$dataset)[x]
  y_name <- colnames(suffStat$dataset)[y]
  formula_YXS_str <- paste(y_name, "~ 1 + ", x_name)
  formula_YS_str <- paste(y_name, "~ 1")

  if (length(S) > 0) {
    S_names <- colnames(suffStat$dataset)[S]
    formula_YXS_str <- paste(formula_YXS_str, "+", paste(S_names, collapse= " + "))
    formula_YS_str <- paste(formula_YS_str, "+", paste(S_names, collapse= " + "))
  }

  if (!is.null(randStr)) {
    formula_YXS_str <- paste0(formula_YXS_str, " | ", randStr)
    formula_YS_str <- paste0(formula_YS_str, " | ", randStr)
  }

  formula_YXS <- stats::as.formula(formula_YXS_str)
  formula_YS <- stats::as.formula(formula_YS_str)

  return(list(formula_YS=formula_YS, formula_YXS=formula_YXS))
}



simpleZeroInflNegBinCITest <- function (x, y, S, suffStat) {
  formulae <- getFitFormulae(x, y, S, suffStat, randStr="1")
  mod1 <- pscl::zeroinfl(formulae$formula_YXS,
                         data=suffStat$dataset, dist = "negbin")
  xids <- which(grepl(colnames(suffStat$dataset)[x],
                      rownames(summary(mod1)$coefficients$count)))
  p1 <- min(summary(mod1)$coefficients$count[xids,4])

  return(list(p=p1, mod0=NULL, mod1=mod1))
}



#' @importFrom MASS glm.nb
#' @importFrom pscl zeroinfl
#' @importFrom stats anova pchisq
#' @importFrom lmtest lrtest
#' @export zeroInflNegBinCITest
zeroInflNegBinCITest <- function (x, y, S, suffStat) {
  ydat = suffStat$dataset[, y]
  xdat = suffStat$dataset[, x]

  mod0 <- NULL
  if (is.null(S) || length(S) == 0) {
    if (length(which(ydat == 0)) > 0) {
      mod1 <- pscl::zeroinfl(ydat ~ xdat | xdat, dist = "negbin")
      if (modListErrors(list(mod1)) ||
          any(is.na(summary(mod1)$coefficients$count[,4]))) {
        mod1 <- pscl::zeroinfl(ydat ~ xdat | 1, dist = "negbin")
        if (modListErrors(list(mod1)) ||
            any(is.na(summary(mod1)$coefficients$count[,4]))) {
          mod1 <- MASS::glm.nb(ydat ~ xdat)
        }
      }
    } else {
      mod1 <- MASS::glm.nb(ydat ~ xdat)
    }

    if (inherits(mod1, "zeroinfl")) {
      #p1 <- log(summary(mod1)$coefficients$count[2,4])
      xids <- which(grepl(colnames(suffStat$dataset)[x],
                  rownames(summary(mod1)$coefficients$count)))
      p1 <- min(summary(mod1)$coefficients$count[xids,4])
    } else {
      t1 <- mod1$null.deviance - mod1$deviance
      dof1 <- length(mod1$coefficients) - 1
      p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
    }

  } else {
    ds0 <- suffStat$dataset[, S, drop = FALSE]
    ds1 <- suffStat$dataset[, c(S, x)]

    a1 <- tryCatch (
    {
      if (length(which(ydat == 0)) > 0) {
        mod1 <- pscl::zeroinfl(ydat ~ . | xdat, data = ds1, dist = "negbin")
        mod0 <- pscl::zeroinfl(ydat ~ . | 1, data = ds0, dist = "negbin")
        if (modListErrors(list(mod0, mod1)) ||
            any(is.na(summary(mod1)$coefficients$count[,4]))) {
          mod1 <- pscl::zeroinfl(ydat ~ . | 1, data = ds1, dist = "negbin")
          mod0 <- pscl::zeroinfl(ydat ~ . | 1, data = ds0, dist = "negbin")
        }

        if (!modListErrors(list(mod0, mod1)) &&
            !any(is.na(summary(mod1)$coefficients$count[,4]))) {
          a1 <- lmtest::lrtest(mod0, mod1)
        } else {
          mod1 <- MASS::glm.nb(ydat ~ ., data = ds1)
          mod0 <- MASS::glm.nb(ydat ~ ., data = ds0)
          if (!modListErrors(list(mod0, mod1))) {
            a1 <- stats::anova(mod0, mod1)
          } else {
            a1 <- NULL
          }
        }
      } else {
        mod1 <- MASS::glm.nb(ydat ~ ., data = ds1)
        mod0 <- MASS::glm.nb(ydat ~ ., data = ds0)
        a1 <- stats::anova(mod0, mod1)
      }
    }, error=function(cond) {
      message("Here's the original error message:")
      message(cond)
      return(NULL)
    })

    if (!is.null(a1)) {
      t1 <- a1[2, 4]
      dof1 <- a1[2, 3]
      p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
    } else {
      p1 <- NA
    }
  }

  return(list(p=p1, mod0=mod0, mod1=mod1))
}

#' @importFrom MASS glm.nb
#' @importFrom stats anova pchisq
#' @export negBinCITest
negBinCITest <- function (x, y, S, suffStat) {
  ydat = suffStat$dataset[, y]
  xdat = suffStat$dataset[, x]

  mod0 <- NULL
  if (is.null(S) || length(S) == 0) {
    mod1 <- MASS::glm.nb(ydat ~ xdat)
    t1 <- mod1$null.deviance - mod1$deviance
    dof1 <- length(mod1$coefficients) - 1
    p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
  } else {
    ds0 <- suffStat$dataset[, S, drop = FALSE]
    ds1 <- suffStat$dataset[, c(S, x)]
    mod0 <- MASS::glm.nb(ydat ~ ., data = ds0)
    mod1 <- MASS::glm.nb(ydat ~ ., data = ds1) #, start=c(coef(mod0), 0))
    a1 <- stats::anova(mod0, mod1)
    if (!is.null(a1)) {
      t1 <- a1[2, 4]
      dof1 <- a1[2, 3]
      p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
    } else {
      p1 <- NA
    }
  }
  return(list(p=p1, mod0=mod0, mod1=mod1))
}

#' @importFrom stats anova pchisq glm poisson
#' @export poissonCITest
poissonCITest <- function (x, y, S, suffStat) {
  ydat = suffStat$dataset[, y]
  xdat = suffStat$dataset[, x]

  mod0 <- NULL
  if (is.null(S) || length(S) == 0) {
    mod1 <- stats::glm(ydat ~ xdat, stats::poisson)
    t1 <- mod1$null.deviance - mod1$deviance
    dof1 <- length(mod1$coefficients) - 1
    p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
  } else {
    ds0 <- suffStat$dataset[, S, drop = FALSE]
    ds1 <- suffStat$dataset[, c(S, x)]
    mod1 <- stats::glm(ydat ~ ., data = ds1, stats::poisson)
    mod0 <- stats::glm(ydat ~ ., data = ds0, stats::poisson)
    a1 <- stats::anova(mod0, mod1)
    if (!is.null(a1)) {
      t1 <- a1[2, 4]
      dof1 <- a1[2, 3]
      p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
    } else {
      p1 <- NA
    }
  }
  return(list(p=p1, mod0=mod0, mod1=mod1))
}

#' @importFrom stats anova pchisq pf lm
#' @export gaussianCITest
gaussianCITest <- function (x, y, S, suffStat) {
  ydat = suffStat$dataset[, y]
  xdat = suffStat$dataset[, x]

  if (is.null(S) || length(S) == 0) {
    mod1 <- stats::lm(ydat ~ xdat)
    a1 <- stats::anova(mod1)
    t1 <- a1[1, 4]
    dof1 <- a1[1, 1]
    p1 <- stats::pf(t1, dof1, a1[2, 1], lower.tail = FALSE, log.p = FALSE)
  } else {
    ds0 <- suffStat$dataset[, S, drop = FALSE]
    ds1 <- suffStat$dataset[, c(S, x)]
    mod1 <- stats::lm(ydat ~ ., data = ds1)
    if (any(is.na(mod1$coefficients))) {
      p1 <- log(1)
      t1 <- 0
    } else {
      a1 <- stats::anova(mod1)
      d1 <- dim(a1)[1] - 1
      t1 <- a1[d1, 4]
      dof1 <- a1[d1, 1]
      df2 <- a1[d1 + 1, 1]
      p1 <- stats::pf(t1, dof1, df2, lower.tail = FALSE, log.p = FALSE)
    }
  }
  return(list(p=p1, mod1=mod1))
}

#' @importFrom MXM ordinal.reg
#' @importFrom stats pchisq
#' @export ordinalCITest
ordinalCITest <- function (x, y, S, suffStat) {
  ydat = suffStat$dataset[, y]
  xdat = suffStat$dataset[, x]

  if (is.null(S) || length(S) == 0) {
    mod1 <- MXM::ordinal.reg(ydat ~ xdat)
    mod0 <- MXM::ordinal.reg(ydat ~ 1)
    t1 <- mod0$devi - mod1$devi
    dof1 <- abs(length(mod1$be) - length(mod0$be))
    p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
  } else {
    ds0 <- suffStat$dataset[, S, drop = FALSE]
    ds1 <- suffStat$dataset[, c(S, x)]

    mod1 <- MXM::ordinal.reg(ydat ~ ., data = ds1)
    mod0 <- MXM::ordinal.reg(ydat ~ ., data = ds0)
    t1 <- mod0$devi - mod1$devi
    if (t1 < 0) {
      t1 <- 0
    }
    dof1 <- abs(length(mod1$be) - length(mod0$be))
    p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
  }
  return(list(p=p1, mod0=mod0, mod1=mod1))
}

#' @importFrom nnet multinom
#' @importFrom stats anova pchisq deviance coef
#' @export multinomialCITest
multinomialCITest <- function (x, y, S, suffStat) {
  ydat = suffStat$dataset[, y]
  xdat = suffStat$dataset[, x]

  if (is.null(S) || length(S) == 0) {
    mod1 <- nnet::multinom(ydat ~ xdat,  trace = FALSE)
    mod0 <- nnet::multinom(ydat ~ 1, trace = FALSE)
    a1 <- stats::anova(mod1, mod0)
    t1 <- a1[2, 6]
    dof1 <- a1[2, 5]
    p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
  } else {
    ds0 <- suffStat$dataset[, S, drop = FALSE]
    ds1 <- suffStat$dataset[, c(S, x)]

    mod1 <- nnet::multinom(ydat ~ ., data = ds1, trace = FALSE)
    mod0 <- nnet::multinom(ydat ~ ., data = ds0, trace = FALSE)
    t1 <- stats::deviance(mod0) - stats::deviance(mod1)
    dof1 <- length(stats::coef(mod1)) - length(stats::coef(mod0))
    p1 <- stats::pchisq(t1, dof1, lower.tail = FALSE, log.p = FALSE)
  }
  return(list(p=p1, mod0=mod0, mod1=mod1))
}


mixedCITestHelper <- function(x, y, S, suffStat, verbose=FALSE) {
  ydat <- suffStat$dataset[, y]
  if (isBinary(ydat)) {
    if (verbose) {
      cat("Running logistic regression for ", y, "\n")
    }
    ret <- logisticCITest(x,y,S,suffStat)
  } else if (is.numeric(ydat) && suffStat$types[y] == "mb_count") {
    if (verbose) {
      cat("Running mb_count regression for ", y)
    }

    if (is.null(suffStat$count_regr)) {
      stop(paste0("It is necessary to specify count_regr for the counting variable ", y))
    }
    if (suffStat$count_regr == "linda") {
      if (verbose) {
        cat(" using linda.\n")
      }
      ret <- lindaCITest(x,y,S,suffStat)
    } else if (suffStat$count_regr == "zicoseq") {
      if (verbose) {
        cat(" using zicoseq.\n")
      }
      ret <- zicoSeqCITest(x,y,S,suffStat)
    } else if (suffStat$count_regr == "zinb") {
      if (verbose) {
        cat("using zero-inflated nb.\n")
      }
      ret <- simpleZeroInflNegBinCITest(x, y, S, suffStat)
    }
  } else if (is.numeric(ydat) && suffStat$types[y] == "count") {
    if (verbose) {
      cat("Running count regression for ", y)
    }

    if (is.null(suffStat$count_regr)) {
      stop(paste0("It is necessary to specify count_regr for the counting variable ", y))
    }
    if (suffStat$count_regr == "simplzeroinfl") {
      if (verbose) {
        cat("using (simple) zero-inflated nb.\n")
      }
      ret <- simpleZeroInflNegBinCITest(x,y,S,suffStat)
    } else if (suffStat$count_regr == "zeroinfl") {
      if (verbose) {
        cat("using zero-inflated nb.\n")
      }
      ret <- zeroInflNegBinCITest(x,y,S,suffStat)
    } else if (suffStat$count_regr == "nb") {
      if (verbose) {
        cat("using nb.\n")
      }
      ret <- negBinCITest(x,y,S,suffStat)
    } else if (suffStat$count_regr == "poisson") {
      if (verbose) {
        cat("using poisson.\n")
      }
      ret <- poissonCITest(x,y,S,suffStat)
    } else {
      stop(paste0("No C.I. test of type ", suffStat$count_regr, " has been implemented yet!"))
    }
  } else if (is.numeric(ydat)) {
    if (verbose) {
      cat("Running Gaussian linear regression for ", y, "\n")
    }
    ret <- gaussianCITest(x,y,S,suffStat)
  } else if (is.ordered(ydat)) {
    if (verbose) {
      cat("Running ordinal regression for ", y, "\n")
    }
    ret <- ordinalCITest(x,y,S,suffStat)
  } else if ((is.factor(ydat) & !is.ordered(ydat))) {
    if (verbose) {
      cat("Running multinomial regression for ", y, "\n")
    }
    ret <- multinomialCITest(x,y,S,suffStat)
  } else {
    stop(paste0("No C.I. test for outcome ", y, " has been implemented yet!"))
  }
  return(ret)
}

#' Gets the p-value for the conditional independence of X and Y given \eqn{S \cup C},
#' where C is a fixed set of covariates.
#' @param suffStat list with the following entries:
#'     dataset: data.frame with all variables that are nodes in the graph
#'     covs: data.frame with a fixed set of covariates that will be
#'           part of the set S in all conditional independence tests
#'     citestResults: pre-computed conditional independence tests
#'                    in a dataset with columns X, Y, S, and pvalue
#'     symmetric: boolean indicating whether both I(X,Y;S) and I(Y,X;S)
#'                should be computed
#'     retall: boolean indicating whether only a p-value (retall=FALSE)
#'             or all computed statistics should be returned (retall=TRUE).
#'     comb_p_method: if "tsagris18" or NULL, then pmm = min(2* min(p1, p2), max(p1, p2)).
#'                    If "min", then pmm = min(p1, p2)
#'     verbose
#' @export mixedCITest
mixedCITest <- function(x, y, S, suffStat) {

  if (!is.null(suffStat$citestResults)) {
    SxyStr <- getSepString(sort(S))
    resultsxys <- c()
    if (!is.null(suffStat$citestResults)) {
      resultsxys <- subset(suffStat$citestResults, X == x & Y == y & S == SxyStr)
      resultsyxs <- subset(suffStat$citestResults, X == y & Y == x & S == SxyStr)
      resultsxys <- rbind(resultsxys, resultsyxs)
    }

    if (!is.null(resultsxys) && nrow(resultsxys) > 0) {
      if (suffStat$verbose) {
        cat("Returning pre-computed p-value for  X=", x, "; Y=", y,
            "given S={", paste0(S, collapse = ","), "}\n")
      }
      # the test should be the symmetric for X,Y|S and Y,X|S
      return(resultsxys[1, "pvalue"])
    }
  }

  ##########################
  # Adding covariates to S #
  ##########################
  if (suffStat$verbose) {
    cat("Testing if X=", x, "is indep of Y=", y, "given S={", paste0(S, collapse = ","), "}\n")
  }
  dat <- suffStat$dataset
  if (!is.null(suffStat$covs) && ncol(suffStat$covs) > 0) {
    dat <- cbind(suffStat$dataset, suffStat$covs)
    cov_ids <- (ncol(suffStat$dataset)+1):(ncol(suffStat$dataset)+ncol(suffStat$covs))
    if (is.null(S)) {
      S <- cov_ids
    } else {
      S <- c(S, cov_ids)
    }
  }

  suffStat2 <- suffStat
  suffStat2$dataset <- dat

  #######################################################
  # Applying tests for X,Y given S and then Y,X given S #
  #######################################################

  ret1 <- mixedCITestHelper(x,y,S,suffStat2, suffStat$verbose)
  p1 <- ret1$p
  p <- p1
  ret <- ret1

  if (!is.null(suffStat$symmetric) &&  suffStat$symmetric == TRUE) {
    ret2 <- mixedCITestHelper(y,x,S,suffStat2, suffStat$verbose)
    p2 <- ret2$p
    cat("p1:", p1, "and p2:", p2, "\n")
    if (!is.null(suffStat$comb_p_method) &&
        suffStat$comb_p_method == "min") {
      p <- min(p1, p2, na.rm = T)
    } else { # NULL or tsagris18
      minp <- min(p1, p2, na.rm = T)
      maxp <- max(p1, p2, na.rm = T)
      p <- min(2* minp, maxp, na.rm=T)
    }
    if (!is.na(p) && !is.na(p2) && p == p2) {
      ret = ret2
    }
  }

  if (!is.null(suffStat$retall) && suffStat$retall == TRUE) {
    if (suffStat$symmetric == TRUE) {
      return(list(ret=ret, ret1=ret1, ret2=ret2))
    } else {
      return(list(ret=ret))
    }
  } else {
    return(p)
  }
}

# Initializes a citestResults data frame for all pairs of nodes and possible
# conditioning sets with length up to m.max
# If csv_citestResults_file exists, then citestResults includes any
# precomputed results recorded in such a file.
#' @importFrom doFuture `%dofuture%`
#' @export initializeCITestResults
initializeCITestResults <- function(p, m.max=Inf,
                                    csv_citestResults_file=NULL,
                                    #file_type = ".csv", # TODO: accept RData
                                    computeProbs = FALSE) {

  col_names <- c("ord", "X", "Y", "S", "pvalue")
  if (computeProbs) {
    col_names <- c(col_names, "pH0", "pH1")
  }

  citestResults_prev <- NULL
  if (!is.null(csv_citestResults_file) && file.exists(csv_citestResults_file)) {
    citestResults_prev <- readCITestResultsCSVFile(csv_citestResults_file)
    if (!computeProbs && ncol(citestResults_prev) > 5) {
        citestResults_prev <- citestResults_prev[,1:5]
    } else if (computeProbs && length(colnames(citestResults_prev)) == 5) {
      citestResults_prev <- cbind(citestResults_prev, pH0=NA, pH1=NA)
    }
  }

  if (!is.null(citestResults_prev) &&
      (length(colnames(citestResults_prev)) != length(col_names) ||
      !all(colnames(citestResults_prev) == col_names))) {
    warning("Pre-computed citestResults are imcompatible.\nStarting with an empty citestResults.")
    citestResults_prev <- NULL
  }


  if (is.infinite(m.max) || m.max > p-2) {
    m.max <- p-2
  }

  pairs <- mycombn(1:p, 2)
  citestResults <-
    foreach (pair_i = 1:ncol(pairs), .combine=rbind.data.frame) %:%
      foreach (csetsize = 0:m.max, .combine=rbind.data.frame) %:%
        foreach (S_i = 1:ncol(mycombn(setdiff(1:p, pairs[,pair_i]), csetsize)),
                 .combine=rbind.data.frame) %dofuture% {
          pair <- pairs[,pair_i]
          Svars <- mycombn(setdiff(1:p, pairs[,pair_i]), csetsize)
          S <- Svars[,S_i]
          ord <- length(S)
          x = pair[1]
          y = pair[2]
          if (computeProbs) {
            data.frame(ord=ord, X=x, Y=y, S=getSepString(S),
                     pvalue=NA, pH0=NA, pH1=NA)
          } else {
            data.frame(ord=ord, X=x, Y=y, S=getSepString(S),
                       pvalue=NA)
          }
        }

  citestResults <- rbind(citestResults_prev, citestResults)
  citestResults <- citestResults[which(!duplicated(citestResults[,1:4])),]
  citestResults <- citestResults[order(citestResults$ord),]

  return(citestResults)
}

mycombn <- function(x, m) {
  if (length(x) == 1) {
    return(combn(list(x),m))
  } else {
    return(combn(x,m))
  }
}

#' @export readCITestResultsCSVFile
readCITestResultsCSVFile <- function(csvfile) {
  #citestResults <- read.csv(csvfile, header=T)
  citestResults <- read.csv(csvfile, header = T, colClasses=c("S"="character"))
  #citestResults <- citestResults[order(citestResults$ord),]
  return(citestResults)
}


#' @importFrom doFuture `%dofuture%`
#' @export getAllCITestResults


getAllCITestResults <- function(dat, indepTest, suffStat, m.max=Inf,
                                computeProbs = FALSE,
                                saveFiles = FALSE,
                                fileid = NULL,
                                citestResults_folder="./tmp/") {
  p <- ncol(dat)
  n <- nrow(dat)

  csv_citestResults_file <- NULL
  if (saveFiles) {
    if (is.null(fileid)) {
      fileid <- paste0("citestResults_", format(Sys.time(), '%Y%m%d_%H%M%S%OS.3'))
    }
    csv_citestResults_file <- paste0(citestResults_folder, "citestResults_", fileid, ".csv")
    if (!file.exists(citestResults_folder)) {
      dir.create(citestResults_folder, recursive = TRUE)
    }
  }

  citestResults <- initializeCITestResults(p, m.max, csv_citestResults_file,
                                           computeProbs=computeProbs)
  #table(citestResults$S)

  if (is.infinite(m.max) || m.max > p-2) {
    m.max <- p-2
  }

  # write.csv(citestResults[complete.cases(citestResults), ],
  #           file=csv_citestResults_file, row.names = F)

  if (saveFiles & length(which(complete.cases(citestResults))) == 0) {
    write.csv(citestResults[NULL,],
              file=csv_citestResults_file, row.names = F)
  }

  todo_citestResults <- citestResults[!complete.cases(citestResults), ]
  if (nrow(todo_citestResults) > 0) {
    new_citestResults <- foreach (i = 1:nrow(todo_citestResults),
                                  .combine=rbind.data.frame) %dofuture% {
                                    ord <- todo_citestResults[i, c("ord")]
                                    x <- todo_citestResults[i, c("X")]
                                    y <- todo_citestResults[i, c("Y")]
                                    S <- getSepVector(todo_citestResults[i, "S"])
                                    SxyStr <- getSepString(S)
                                    pvalue <- indepTest(x, y, S, suffStat = suffStat)
                                    if (computeProbs) {
                                      probs <- pvalue2probs(pvalue, n=n)
                                      pH0 <- probs$pH0
                                      pH1 <- probs$pH1
                                      ret <- data.frame(ord=ord, X=x, Y=y, S=SxyStr,
                                                        pvalue = pvalue, pH0=pH0, pH1=pH1)
                                    } else {
                                      ret <- data.frame(ord=ord, X=x, Y=y, S=SxyStr,
                                                        pvalue = pvalue)
                                    }
                                    if (saveFiles) {
                                      write.table(ret, file=csv_citestResults_file, sep=",", row.names = FALSE,
                                                  col.names = FALSE, append = TRUE)
                                    }
                                    ret
                                  }
    citestResults <- rbind(citestResults[complete.cases(citestResults), ],
                           new_citestResults)
  }

  citestResults <- citestResults[order(citestResults$ord),]
  rownames(citestResults) <- NULL

  return(citestResults)
}


getAllCITestResultsOLD <- function(dat, indepTest, suffStat, m.max=Inf,
                                computeProbs = FALSE,
                                saveFiles = FALSE,
                                fileid = NULL,
                                citestResults_folder="./tmp/") {
  p <- ncol(dat)
  n <- nrow(dat)

  csv_citestResults_file <- NULL
  if (saveFiles) {
    tmp_fileid <- NULL
    if (is.null(fileid)) {
      fileid <- format(Sys.time(), '%Y%m%d_%H%M%S%OS.3')
      tmp_fileid <- fileid
    } else {
      tmp_fileid <- paste0(fileid, "_", format(Sys.time(), '%Y%m%d_%H%M%S%OS.3'))
    }

    tmp_folder = paste0(citestResults_folder, "tmp/")
    if (!file.exists(tmp_folder)) {
      dir.create(tmp_folder, recursive = TRUE)
    }
    partial_csv_citestResults_file <- paste0(tmp_folder, "partial_citestResults_",
                                             tmp_fileid, ".csv")

    csv_citestResults_file <- paste0(citestResults_folder, "citestResults_", fileid, ".csv")
  }

  citestResults <- initializeCITestResults(p, m.max, csv_citestResults_file,
                                           computeProbs=computeProbs)
  #table(citestResults$S)

  if (is.infinite(m.max) || m.max > p-2) {
    m.max <- p-2
  }

  if (saveFiles) {
    write.csv(citestResults[NULL,], file=partial_csv_citestResults_file, row.names = F)
  }

  pairs <- mycombn(1:p, 2)
  new_citestResults <-
    foreach (pair_i = 1:ncol(pairs), .combine=rbind.data.frame) %:%
    foreach (csetsize = 0:m.max, .combine=rbind.data.frame) %:%
    foreach (S_i = 1:ncol(mycombn(setdiff(1:p, pairs[,pair_i]), csetsize)),
             .combine=rbind.data.frame) %dofuture% {
      pair <- pairs[,pair_i]
      Svars <- mycombn(setdiff(1:p, pair), csetsize)
      S <- as.numeric(Svars[,S_i, drop=FALSE])
      ord <- length(S)
      x = pair[1]
      y = pair[2]
      SxyStr <- getSepString(S)

      curid <- which(citestResults$X == x & citestResults$Y == y &
                       citestResults$S == SxyStr)

      if (length(curid) > 0) {
        if (is.na(citestResults[curid, c("pvalue")])) {
          pvalue <- indepTest(x, y, S, suffStat = suffStat)
        } else {
          pvalue <- citestResults[curid, "pvalue"]
        }

        if (computeProbs) {
          if (any(is.na(citestResults[curid, c("pH0", "pH1")]))) {
            probs <- pvalue2probs(pvalue, n=n)
            pH0 <- probs$pH0
            pH1 <- probs$pH1
          } else {
            pH0 <- citestResults[curid, "pH0"]
            pH1 <- citestResults[curid, "pH1"]
          }
          ret <- data.frame(ord=ord, X=x, Y=y, S=SxyStr,
                            pvalue = pvalue, pH0=pH0, pH1=pH1)
        } else {
          ret <- data.frame(ord=ord, X=x, Y=y, S=SxyStr,
                            pvalue = pvalue)
        }

        if (saveFiles) {
          write.table(ret, file=partial_csv_citestResults_file, row.names = FALSE,
                    col.names = FALSE, append = TRUE)
        }
        ret
      } else {
        warning("Problem with citestResults for X = ", x, " Y = ", y, " and S = ", SxyStr)
      }
   }

  citestResults <- new_citestResults[order(new_citestResults$ord),]
  rownames(citestResults) <- NULL

  if (saveFiles) {
    write.csv(citestResults, file=csv_citestResults_file, row.names = FALSE)
  }

  return(citestResults)
}


# dat contains only variables that are represented as nodes in the graph
#' @export runAllCITests
runAllCITests <- function(dat, indepTest, suffStat,
                          m.max=Inf, alpha = 0.05,
                          partial_results_file=NULL) {
  tested_independencies <- test_all_cindeps(indepTest, samples=dat,
                                            alpha=alpha, max_csetsize = m.max,
                                            suffStat=suffStat,
                                            partial_results_file=partial_results_file)
  citestResults <- convertToCITestResults(tested_independencies)
  citestResults[, c(1,2,3,5)] <- lapply(citestResults[, c(1,2,3,5)], as.numeric)

  return(citestResults)
}

#' @export runAllMixedCITests
runAllMixedCITests <- function(dat, vars_names, covs_names=c(),
                               m.max=Inf, alpha = 0.05,
                               partial_results_file=NULL) {
  indepTest <- mixedCITest
  suffStat <- getMixedCISuffStat(dat, vars_names, covs_names)
  vars_df <- dat[,vars_names, drop=FALSE]
  citestResults <- runAllCITests(vars_df, indepTest, suffStat, m.max, alpha,
                                 partial_results_file=partial_results_file)

  return(citestResults)
}


test_indeps_helper <- function(test_function, test_data, n, i, j, csetsize,
                               cur_tested_independences = NULL) {
  #start with empty set
  csetvec <- rep(0, n)

  csetvec[index(1,csetsize)] <- 1
  tested_independences <- list()

  while ( !any(is.na(csetvec) ) ) {
    runTest <- csetvec[i]==0 && csetvec[j] == 0
    cond_vars <- which(csetvec==1)
    cset<-bin.to.dec(rev(csetvec))

    if (runTest) { #only if neither x and y are cond.
      cat(i, " ", j , "|", cond_vars, "\n")
      entries <- NULL
      if (!is.null(cur_tested_independences)) {
        entries <- which(sapply(cur_tested_independences,
                              function(x) {
                                all(x$vars == sort(c(i,j))) &&
                                  length(cond_vars) == length(x$C) &&
                                  all(cond_vars %in% x$C)
                              }))
      }

      if (!is.null(entries) && length(entries) == 1) {
        test_result <- cur_tested_independences[[entries]]
      } else {
        test_result <- list()
        test_result$vars<-sort(c(i,j))
        test_result$C <- cond_vars

        #calling the test function
        test_result$p <- test_function(i, j, cond_vars, test_data$suffStat)

        #if it is bigger than the result we have independence
        test_result$independent <- ( test_result$p > test_data$p_threshold )

        #weight is always 1
        test_result$w <- 1 # TODO change for some of my scores

        #put some parameters right
        test_result$J<-test_data$J
        test_result$jset<-test_data$jset
        test_result$cset<-cset
        test_result$M<-setdiff((1:n),c(test_result$vars,test_result$C))
        test_result$mset <- getm( test_result$vars, test_result$C, n=n)
      }
      #cat(paste(test_result$M,collapse=','),'=',test_result$mset,'\n')

      #adding the test result also to tested_independences vector
      tested_independences[[length(tested_independences) + 1]] <- test_result
    } #if x and y are not in the conditioning set

    #consider next csetvec given by the following function
    csetvec<-next_colex_comb(csetvec)
  } #while csetvec != NA
  tested_independences
}

#maxcset = schedule
getTestData <- function(samples, alpha, suffStat) {
  D <- list()
  D[[1]] <- list()
  D[[1]]$data <- samples
  D[[1]]$e <- rep(0, ncol(samples)) # 0: observational; 1: interventional
  D[[1]]$N <- nrow(samples)

  test_data <- list()
  data <- D[[1]]

  # Preparing for writing indep constraints.
  jindex <- 0 # only one dataset
  test_data$jset <- bin.to.dec(rev(1*(data$e==1)))
  test_data$J <- which(data$e==1)
  test_data$names <- colnames(data$data)
  test_data$indPath <- NULL

  test_data$p_threshold<- alpha

  # setting up suffStat with the conditional indep test parameters
  test_data$suffStat <- suffStat
  return(test_data)
}


convertToCITestResults <- function(tested_independences) {
  citestResults <- data.frame()

  if (!is.null(tested_independences)) {
    for (l in tested_independences){
      vars <- sort(l$vars)
      x <- vars[1]
      y <- vars[2]
      Sxy <- getSepString(l$C)
      ord <- length(l$C)
      pvalue <- l$p
      cur_row <- c(ord=ord, X=x, Y=y, S=Sxy, pvalue=pvalue)
      citestResults <- rbind(citestResults, cur_row)
    }
  }
  colnames(citestResults) <- c("ord", "X", "Y", "S", "pvalue")
  return(citestResults)
}


#' @export getMixedCISuffStat
getMixedCISuffStat <- function(dat, vars_names, covs_names=c(), verbose=TRUE) {
  vars_df <- dat[,vars_names, drop=FALSE]
  covs_df <- dat[,covs_names, drop=FALSE]

  types <- sapply(vars_df, class)

  suffStat <- list(dataset=vars_df,
                   covs=covs_df,
                   rand_varnames = c(), # for simpl
                   n = dim(vars_df)[1],
                   retall = FALSE,
                   symmetric = TRUE,
                   comb_p_method = "tsagris18",
                   packages_list = c(), # change to packages list required for parallelization
                   types=types,
                   verbose=verbose)

  return(suffStat)
}


# n: number of observed variables
# samples=vars_df
# test_function = indepTest
#' @importFrom doFuture `%dofuture%`
#' @export test_all_cindeps
test_all_cindeps <- function(test_function, samples, alpha, suffStat,
                             max_csetsize=Inf, n=NULL,
                             partial_results_file=NULL) {
  if (!is.null(partial_results_file) && file.exists(partial_results_file)) {
    load(partial_results_file) # loading tested_independencies
    cur_tested_independences <- tested_independences
    save(cur_tested_independences, file=paste0(partial_results_file, ".tmp"))
  } else {
    cur_tested_independences <- NULL
  }

  tested_independences=list()

  if (is.null(n)) {
    n <- ncol(samples)
  }

  if (is.infinite(max_csetsize)) {
    max_csetsize <- n-2
  }

  test_data <- getTestData(samples, alpha, suffStat)

  # Function for conducting all independence tests for one data set.
  for (csetsize in index(0, max_csetsize)) { #go from simpler to more complicated tests
    for (i in 1:(n-1)) {
      tested_independences_j <-
        foreach (j = (i+1):n, .combine=rbind.data.frame, .verbose = TRUE) %dofuture% {
          curtest <- test_indeps_helper(test_function, test_data, n,
                                        i, j, csetsize,
                                        cur_tested_independences)
          return(curtest)
        }


      # #pb <- txtProgressBar(min = 1, max = length((i+1):n)+1, style = 3)
      # #progress <- function(n) setTxtProgressBar(pb, n)
      # #opts <- list(progress = progress)
      # tested_independences_j <-
      #   foreach (j = (i+1):n, #, #.options.snow = opts,
      #            .verbose = TRUE, .export = ls(globalenv()),
      #            .packages=suffStat$packages_list
      #             ) %dopar% {
      #              curtest <- test_indeps_helper(test_function, test_data, n,
      #                                            i, j, csetsize,
      #                                            cur_tested_independences)
      #              #setTxtProgressBar(pb, i)
      #              return(curtest)
      # }

      #for j
      for (test_result_list in tested_independences_j) {
        for (test_result in test_result_list) {
          tested_independences[[length(tested_independences) + 1]] <- test_result
        }
      }
      if (!is.null(partial_results_file)) {
        save(tested_independences, file=partial_results_file)
      }
    } # for i
  } # for csetsize
  tested_independences
}

fakeCITest <- function (x, y, S, suffStat) {
  # - suffStat$tested_independences: a list of independence test results.
  # Given a list of independence test results, it retrieves the correct one.
  # This function is useful in order to pass some test results for tests that were already performed to pcalg.

  if (!is.null(suffStat$tested_independences)) {
    list_indeps <- suffStat$tested_independences
    vars <- sort(c(x,y))

    for (l in list_indeps){
      if (any( vars != l$vars)){
        next
      }
      if (length(S)== 0 && length(l$C) == 0) return (l$p)
      if (length(S)== 0) next
      if (length(l$C) == 0) next
      if (length(l$C) != length(S)) next
      if (any(S != l$C)) next
      return (l$p)
    }
  }

  stop("Asking for a test that was not performed: ", x, ",", y, ", {", paste(S, collapse=","), "}")
}

bin.to.dec <- function( Base.b ) {
  #Changes a binary number into a decimal integer.
  #REMEMBER TO +1 IF USED FOR INDEXING A BINARY CPT!

  ndigits = length(Base.b)
  sum(Base.b*2^((ndigits-1):0))
}

index <- function(from, to) {
  # Indexing help function for R.
  # Similar to matlabs ':', where 3:1 = c() instead of c(3,2,1) of R.
  if ( from > to ) {
    R <- c()
  } else {
    R <- from:to
  }
  R
}

getm <- function(vars, C, n) {
  # Gives a integer representation of the M set given
  # the number of variables and conditioning set C.
  vec <- rep(1, n)
  #vec[x]<-vec[y]<-0
  vec[vars]<-0
  vec[C]<-0
  #bin.to.dec(vec)
  bin.to.dec(rev(vec))
}

next_colex_comb <- function(x) {
  #For a 0-1 vectors gives a next
  #Can be used to quickly iterate over all subset of a given size.
  #Just start with (1,1,1,….,1,0,…,0). In the end returns 0.
  #The underlying mechanism to determine the successor is to determine the lowest block of ones
  #and move its highest bit one position up.
  #The rest of the block is then moved to the low end of the word.
  j=1;
  for ( i in index(1,(length(x)-1)) ) {
    if ( x[i] == 1 ) {
      if ( x[i + 1] == 0 ) {
        x[i]<-0;x[i+1]=1;
        return(x); #switch bit to left
      } else {
        x[i]<-0;x[j]=1;j<-j+1;
      }
    }
  }
  return(NA)
}

#' TODO: make citestResults and object with both the labels and the data.frame
#' @export extractValidCITestResults
extractValidCITestResults <- function(citestResults, cur_varnames, new_varnames) {
  new_citestResults <- data.frame()
  for (i in 1:nrow(citestResults)) {
    cur_row <- citestResults[i, , drop=FALSE]

    cur_xname <- cur_varnames[cur_row$X]
    X <- which(new_varnames == cur_xname)
    if (length(X) == 1) {
      cur_yname <- cur_varnames[cur_row$Y]
      Y <- which(new_varnames == cur_yname)
      if (length(Y) == 1) {
        cur_snames = c()
        S <- c()
        if (length(getSepVector(cur_row$S)) > 0) {
          cur_snames <- cur_varnames[getSepVector(cur_row$S)]
          S <- which(new_varnames %in% cur_snames)
          S <- sort(S)
        }

        if (length(S) == length(cur_snames)) {
          sortedXY <- sort(c(X, Y))
          X <- sortedXY[1]
          Y <- sortedXY[2]
          ord <- cur_row$ord
          stats <- cur_row[,5:length(cur_row), drop=FALSE]
          #pvalue <- cur_row$pvalue

          new_citestResults <- rbind.data.frame(new_citestResults,
                                                c("ord"=ord, "X"=X, "Y"=Y,
                                                     "S"=getSepString(S), stats))
        }
      }
    }
  }
  new_citestResults[,-4] <- lapply(new_citestResults[,-4], as.numeric)

  return(new_citestResults)
}
