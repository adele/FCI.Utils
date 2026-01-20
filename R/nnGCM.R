matrixSqrt <- function(m, inverse=FALSE, tol = sqrt(.Machine$double.eps)) {
  eigDecomp <- eigen(m, symmetric=TRUE)
  Positive <- eigDecomp$values > max(tol * eigDecomp$values[1], 0)
  D = NULL
  if (inverse) {
    D = diag(1/sqrt(eigDecomp$values[Positive]))
  } else {
    D = diag(sqrt(eigDecomp$values[Positive]))
  }
  sqrtM <- eigDecomp$vectors[, Positive, drop = FALSE] %*% D %*% t(eigDecomp$vectors[, Positive, drop = FALSE])
  return(sqrtM)
}

# GCM test function for conditional independence using neural networks
# suffStat must have "dataset", and, optionally, "nruns", and "compute_MC_pvalue"
# this is a symmetric test, so there is no need to run it switching the roles of X and Y
#' @importFrom stats pchisq
#' @importFrom nnet nnet
#' @export nnGCMTest
nnGCMTest <- function(x, y, S, suffStat) {

  require(nnet)

  # Extracting variables
  ydat = suffStat$dataset[, y]
  xdat = suffStat$dataset[, x]
  n <- nrow(suffStat$dataset)
  hidden_layer_units <- suffStat$hidden_layer_units
  if (is.null(hidden_layer_units)) {
    hidden_layer_units = 10
  }

  if (is.null(S) || length(S) == 0) {
    sdat <- matrix(1, nrow = n, ncol = 1)  # intercept only
  } else {
    sdat <- as.matrix(suffStat$dataset[, S, drop = FALSE])
  }

  # Fitting neural network for X ~ Z
  if (isBinary(xdat)) {
    if (is.factor(xdat)) {
      xdat <- as.numeric(xdat) - 1
    }
    nn_X <- nnet(x = sdat, y = xdat, size = hidden_layer_units,
                 linout = FALSE, softmax = FALSE,
                 trace = FALSE, maxit = 1000)
    pred_X <- predict(nn_X, sdat)
    resid_X <- xdat - pred_X
  } else if (isMultinomial(xdat)) {
    nn_X <- nnet(x = sdat, y = class.ind(xdat), size = hidden_layer_units,
                 linout = FALSE, softmax = TRUE,
                 trace = FALSE, maxit = 1000)
    pred_X <- predict(nn_X, sdat)
    resid_X <- class.ind(xdat) - pred_X
  } else if (is.numeric(xdat)) { # x is "continuous"
    nn_X <- nnet(x = sdat, y = xdat, size = hidden_layer_units,
                 linout = TRUE, trace = FALSE, maxit = 1000)
    pred_X <- predict(nn_X, sdat)
    resid_X <- xdat - pred_X
  } else {
    stop("Unknown type for X.")
  }

  # Fitting neural network for Y ~ Z
  if (isBinary(ydat)) {
    nn_Y <- nnet(x = sdat, y = ydat, size = hidden_layer_units,
                 linout = FALSE, softmax = FALSE,
                 trace = FALSE, maxit = 1000)
    pred_Y <- predict(nn_Y, sdat)
    resid_Y <- ydat - pred_Y
  } else if (isMultinomial(ydat)) {
    nn_Y <- nnet(x = sdat, y = class.ind(ydat), size = hidden_layer_units,
                 linout = FALSE, softmax = TRUE,
                 trace = FALSE, maxit = 1000)
    pred_Y <- predict(nn_Y, sdat)
    resid_Y <- class.ind(ydat) - pred_Y
  } else if (is.numeric(ydat)) { # y is "continuous"
    nn_Y <- nnet(x = sdat, y = ydat, size = hidden_layer_units,
                 linout = TRUE, trace = FALSE, maxit = 1000)
    pred_Y <- predict(nn_Y, sdat)
    resid_Y <- ydat - pred_Y
  } else {
    stop("Unknown type for Y.")
  }

  # Computing GCM statistic
  # If X and Y continuous/binary scalar, resid_X and resid_Y are n x 1
  # If multinomial, resid_X or resid_Y are n x k matrices
  if (is.vector(resid_X)) resid_X <- as.matrix(resid_X)
  if (is.vector(resid_Y)) resid_Y <- as.matrix(resid_Y)

  n <- NROW(resid_X)
  dx <- NCOL(resid_X)
  dy <- NCOL(resid_Y)
  p <- dx * dy

  # This follows a Gaussian distribution with some covariance matrix
  R <- c()
  for (i in 1:n) {
    R <- rbind(R, as.vector(resid_X[i, ] %*% t(resid_Y[i, ])))
  }
  Sn <- colSums(R) / sqrt(n)

  # Estimating the covariance matrix of R
  CovR <- cov(R)

  # Obtaining the standardized GCM statistic, which
  # follows a standard Gaussian distribution
  if (dim(CovR)[1] == 1) {
    sigma <- sqrt(as.numeric(CovR))
    Tn <- 1/sigma * Sn
    # TODO: Talk to Aaron
    # Tn2 <- (sqrt(n) / sigma) *  colSums(R) # this is not the same
    # Tn3 <- 1 / (sigma * sqrt(n)) *  colSums(R) #  this is the same
  } else {
    CovRinv_sqrt <- matrixSqrt(CovR, inverse = TRUE)
    Tn <- as.vector(CovRinv_sqrt %*% Sn)
  }

  # Quadratic-form test
  Qn <- sum(Tn^2)
  p_value <- 1 - pchisq(Qn, df = p)

  ret1 <- list(chi2stat = Qn, df=p, pvalue = p_value)

  ret2 <- NULL
  if (!is.null(suffStat$compute_MC_pvalue) && suffStat$compute_MC_pvalue) {
    T.stat <- max(abs(Tn))

    nsim <- 1000

    # Bootstrap S_n
    Sn_sim <- t(R) %*% matrix(rnorm(n * nsim), n, nsim) / sqrt(n)   # p x nsim

    # Standardize
    if (p > 1) {
      Tn_sim <- CovRinv_sqrt %*% Sn_sim
    } else {
      Tn_sim <- 1/sigma * Sn_sim
    }

    # Max-type statistics
    Tsim <- apply(abs(Tn_sim), 2, max)

    # Monte Carlo p-value
    p_value_MC <- (sum(Tsim >= T.stat) + 1) / (nsim + 1)

    ret2 <- list(Tstat = T.stat, pvalue.MC = p_value_MC)
  }

  return(c(ret1, ret2))
}

# Runs the test suffStat$nruns times and gets the average of the MC results
nnGCMTest2 <- function(x, y, S, suffStat) {
  suffStat$compute_MC_pvalue <- TRUE

  if (is.null(suffStat$nruns)) {
    suffStat$nruns <- 10
  }

  results <- c()
  for (i in 1:suffStat$nruns) {
    nn.test <- nnGCMTest(x=x, y=y, S=S, suffStat = suffStat)
    results <- rbind(results, unlist(nn.test))
  }
  results <- as.data.frame(results)

  results <- results[order(results[,"Tstat"]), , drop=FALSE]
  pvalueMC_ave <- as.numeric(
    results[which(results[,"Tstat"] >= median(results[,"Tstat"]))[1],"pvalue.MC"])

  ret <- results[which.min(abs(results$pvalue - pvalueMC_ave)), , drop=TRUE]
  return(ret)
}
