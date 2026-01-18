rm(list=ls())

isMultinomial <-  function(x) {
  return(
    !is.ordered(x) &&
      ((is.factor(x) && length(levels(x)) > 2) ||
         !is.numeric(x) && length(unique(x)) < 10))
}

isBinary <- function(x) {
  return(
    !is.ordered(x) &&
      ((is.factor(x) && length(levels(x)) == 2) || length(unique(x)) == 2))
}

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
# suffStat must have: dataset
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

  return(list(chisq_test=ret1, MC_test=ret2))
}


nnGCMTest2 <- function(x, y, S, suffStat) {
  results <- c()
  for (i in 1:10) {
    nn.test <- nnGCMTest(x=x, y=y, S=S, suffStat = suffStat)
    results <- rbind(results,
                         c(unlist(nn.test$chisq_test), unlist(nn.test$MC_test)))
  }
  results <- results[order(results[,"Tstat"]), ]
  results <- as.data.frame(results)
  pvalueMC_ave <- as.numeric(
    results[which(results[,"Tstat"] >= median(results[,"Tstat"]))[1],"pvalue.MC"])
  results <- results[which.min(abs(results$pvalue - pvalueMC_ave)), ]
  return(results)
}

# Load package
library(nnet)

# Generate training data

results <- c()

for (i in 1:100) {
  n <- 1000

  # Simulate Z1, Z2, Z3
  data <- data.frame(
    z1 = rnorm(n, 0, 1),
    z2 = rnorm(n, 0, 1),
    z3 = rnorm(n, 0, 1)
  )

  # Non-linear, binary X
  nlin_pred_x <- sin(data$z1) + data$z2^2 - 0.5 * data$z3
  prob_x <- 1/(1+exp(- nlin_pred_x))
  data$x <- rbinom(n, size = 1, prob = prob_x)
  data$x <- as.factor(data$x)

  # # Non-linear, multinomial X with three levels
  # nlin_pred_x1 <- data$z1 + data$z2^3 + 0.7 * data$z3
  # nlin_pred_x2 <- data$z1 - data$z2 - 0.3 * data$z3^2
  # nlin_pred_x3 <- data$z1 - 0.2 * data$z2 - data$z3^2
  # nlin_pred_x4 <- 0 # reference category
  #
  # nlin_pred_x <- cbind(nlin_pred_x1, nlin_pred_x2, nlin_pred_x3, nlin_pred_x4)
  # prob_x <- t(apply(nlin_pred_x, 1, function(row) exp(row)/sum(exp(row))))
  # # Sample category for each observation
  # data$x <- apply(prob_x, 1, function(p) sample(1:4, 1, prob = p))
  # data$x <- factor(data$x, levels = 1:4)


  # Non-linear, multinomial Y with three levels,
  # fisrt level is a function of X and Z
  nlin_pred_y1 <- c()
  for (i in 1:n) {
    nlin_pred_y1 <- c(nlin_pred_y1,
                      #if (data$x[i] == 1) {
                        data$z1[i] + data$z2[i]^2 - 0.5 * data$z3[i]
                      #}  else {
                      #  data$z1[i] - 0.3 * data$z2[i]^2 + 0.5 * data$z3[i]
                      #}
    )
  }

  nlin_pred_y2 <- data$z1 - data$z2 + 0.3 * data$z3^2
  nlin_pred_y3 <- 0 # reference category

  nlin_pred_y <- cbind(nlin_pred_y1, nlin_pred_y2, nlin_pred_y3)
  prob_y <- t(apply(nlin_pred_y, 1, function(row) exp(row)/sum(exp(row))))
  # Sample category for each observation
  data$y <- apply(prob_y, 1, function(p) sample(1:3, 1, prob = p))
  data$y <- factor(data$y, levels = 1:3)

  data <- data[, c("x", "y", "z1", "z2", "z3")]
  str(data)

  suffStat <- list(dataset=data, compute_MC_pvalue=TRUE)

  # Testing X indep Y given {Z1, Z2, Z3}
  x = 1
  y = 2
  S = 3:5

  nn.test <- nnGCMTest2(x=x, y=y, S=S, suffStat = suffStat)

  results <- rbind(results, nn.test)
}

results <- as.data.frame(results)

print(length(which(results$pvalue < 0.05)))
