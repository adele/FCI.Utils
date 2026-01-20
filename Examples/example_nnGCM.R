rm(list=ls()) # clean the environment

library(FCI.Utils)
library(pcalg)

# library(RNOmni) # useful for normal transformation


# If running in parallel:
run_parallel = TRUE

if (run_parallel) {
  require(doFuture)
  require(future.apply)
  n_cores <- parallel::detectCores() # or change to the number of cores required
  plan("multisession", workers = n_cores) # for running on RStudio
  # plan("multicore", workers = n_cores) # forking -- for running at the server
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
                      if (data$x[i] == 1) {
                        data$z1[i] + data$z2[i]^2 - 0.5 * data$z3[i]
                      }  else {
                        data$z1[i] - 0.3 * data$z2[i]^2 + 0.5 * data$z3[i]
                      }
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


  vars_names <- colnames(data) # nodes of the graph
  covs_names = c() # variables to be always conditioned on
  indepTest <- mixedCITest
  suffStat <- getMixedCISuffStat(dat = data, vars_names = vars_names,
                                 covs_names = covs_names,
                                 method = "nnGCM", verbose = TRUE)
  suffStat$retall = TRUE  # Required if computing probs for dcFCI
  suffStat$nruns = 5
  suffStat$compute_MC_pvalue = FALSE

  # Testing X indep Y given {Z1, Z2, Z3}
  ret <- indepTest(1,2,3:5,suffStat)

  # vars_df <- data[,vars_names, drop=FALSE]
  #
  # # Compute all conditional independence tests in parallel
  # citestResults <- getAllCITestResults(
  #   vars_df, indepTest, suffStat, m.max = 2,
  #   computeProbs = TRUE,        # Required if running dcFCI
  #   eff_size = 0.05,            # Minimum effect size threshold when computeProbs = TRUE
  #   fileid = "",          # Unique identifier for the dataset / analysis
  #   recover_citestResults = FALSE # Set TRUE to resume computation after a crash
  # )

  results <- rbind(results, ret$ret)
}

results <- as.data.frame(results)

print(length(which(results$pvalue < 0.05)))
