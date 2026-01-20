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



# Initializes a citestResults data frame for all pairs of nodes and possible
# conditioning sets with length up to m.max
# If csv_citestResults_file exists, then citestResults includes any
# precomputed results recorded in such a file.
# @importFrom doFuture `%dofuture%` `%:%`
#' @import doFuture
#' @import foreach
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
    futile.logger::flog.info(paste("Reading results from: ",
                                   csv_citestResults_file), name = "citests_log")

    citestResults_prev <- readCITestResultsCSVFile(csv_citestResults_file)
    if (!computeProbs && ncol(citestResults_prev) > 5) {
      citestResults_prev <- citestResults_prev[,1:5]
    } else if (computeProbs && length(colnames(citestResults_prev)) == 5) {
      citestResults_prev <- cbind(citestResults_prev, pH0=NA, pH1=NA)
    }
    futile.logger::flog.info(paste("Loaded citestResults with ",
                                   nrow(citestResults_prev), "rows."), name = "citests_log")
  }

  if (!is.null(citestResults_prev) &&
      (length(colnames(citestResults_prev)) != length(col_names) ||
       !all(colnames(citestResults_prev) == col_names))) {
    futile.logger::flog.info(paste("Pre-computed citestResults are imcompatible."), name = "citests_log")
    citestResults_prev <- NULL
  }

  if (is.null(citestResults_prev)) {
    futile.logger::flog.info(paste("Starting with an empty citestResults."), name = "citests_log")
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

  futile.logger::flog.info(paste("Total number of citests: ", nrow(citestResults)), name = "citests_log")

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


get_last_modified_file <- function(folder_path, pattern) {
  # List only files with pattern
  files <- list.files(folder_path, pattern=pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)  # no files in folder

  # Get file info and find the most recently modified file
  file_info <- file.info(files)
  last_file <- rownames(file_info)[which.max(file_info$mtime)]
  return(last_file)
}


#' @importFrom futile.logger flog.appender appender.file flog.info
#' @importFrom doFuture `%dofuture%`
#' @export getAllCITestResults
getAllCITestResults <- function(dat, indepTest, suffStat, m.max=Inf,
                                computeProbs = FALSE,
                                eff_size=0.05,
                                fileid = "",
                                save_files = TRUE,
                                recover_citestResults = FALSE, # restart from the most recently generated files
                                results_folder= file.path(getwd(), "tmp","citestsResults"), # folder where results will be recovered from and/or saved to
                                log_folder = file.path(getwd(), "tmp", "logs")) {

  if  (suffStat$method == "nnGCM" && computeProbs && !suffStat$retall) {
    stop("Set suffStat$retall = TRUE when using nnGCM and BFF.")
  }


  # Set up logging to file
  if (!file.exists(log_folder))
    dir.create(log_folder, recursive = TRUE)

  futile.logger::flog.appender(futile.logger::appender.file(
    paste0(file.path(log_folder, "citests_log_"), format(Sys.time(), '%Y%m%d_%H%M%S%OS.3'), ".txt")),
    name = "citests_log")


  p <- ncol(dat)
  n <- nrow(dat)

  if (save_files && !file.exists(results_folder)) {
    dir.create(results_folder, recursive = TRUE)
  }

  csv_citestResults_file <- NULL
  if (recover_citestResults) {
    # if there is a citestResults file, then loads the results and starts from there.
    file_pattern <- paste0("^citestResults_", fileid, ".*\\.csv$")
    csv_citestResults_file <- get_last_modified_file(results_folder, pattern = file_pattern)
    futile.logger::flog.info(paste("Recovering results from: ", csv_citestResults_file), name = "citests_log")
  }

  citestResults <- initializeCITestResults(p, m.max, csv_citestResults_file,
                                           computeProbs=computeProbs)
  #table(citestResults$S)

  if (is.infinite(m.max) || m.max > p-2) {
    m.max <- p-2
  }

  if (save_files) {
    # A new file is created with the current time
    csv_citestResults_file <- file.path(results_folder,
                                        paste0("citestResults_", fileid, "_",
                                               format(Sys.time(), '%Y%m%d_%H%M%S%OS.3'), ".csv"))
    write.csv(citestResults[complete.cases(citestResults), ],
              file=csv_citestResults_file, row.names = F)
    futile.logger::flog.info(paste("A new csv_citestResults_file has been created: ",
                                   csv_citestResults_file), name = "citests_log")
  }


  if (computeProbs) {
    # Note: in the case where probs are computed, we mark tests that failed by
    # setting pH0 = pH1 = 0.5, so we can select only those with pH0 = pH1 = NA
    done_citestResults <- citestResults[!is.na(citestResults$pH0), ]
    todo_citestResults <- citestResults[is.na(citestResults$pH0), ]
  } else {
    # Note: Previously computed tests that returned an NA pvalue will always be recomputed.
    done_citestResults <- citestResults[!is.na(citestResults$pvalue), ]
    todo_citestResults <- citestResults[is.na(citestResults$pvalue), ]
  }

  futile.logger::flog.info(paste("Number of already computed citests: ", nrow(done_citestResults)),
                           name = "citests_log")
  futile.logger::flog.info(paste("Number of citests remaining to compute: ", nrow(todo_citestResults)),
                           name = "citests_log")


  if (nrow(todo_citestResults) > 0) {
    new_citestResults <- foreach (i = 1:nrow(todo_citestResults),
                                  .combine=rbind.data.frame) %dofuture% {
                                    ord <- todo_citestResults[i, c("ord")]
                                    x <- todo_citestResults[i, c("X")]
                                    y <- todo_citestResults[i, c("Y")]
                                    S <- getSepVector(todo_citestResults[i, "S"])
                                    SxyStr <- getSepString(S)
                                    test_out <- indepTest(x, y, S, suffStat = suffStat)
                                    chiSqStat = df = NULL # default will be used for LR GLM-based CI Test
                                    if (suffStat$retall) {
                                      pvalue <- test_out$p
                                      if (suffStat$method == "nnGCM" && computeProbs) {
                                        pvalue <- test_out$ret$pvalue
                                        chi2stat <- test_out$ret$chi2stat
                                        df <- test_out$ret$df
                                      }
                                    } else {
                                      pvalue <- test_out
                                    }
                                    if (computeProbs) {
                                      if (is.na(pvalue)) {
                                        ret <- data.frame(ord=ord, X=x, Y=y, S=SxyStr,
                                                          pvalue = NA, pH0=0.5, pH1=0.5)
                                      } else {
                                        probs <- pvalue2probs(pvalue, n=n, eff_size=eff_size,
                                                              chiSqStat = chi2stat, df=df)
                                        pH0 <- probs$pH0
                                        pH1 <- probs$pH1
                                        ret <- data.frame(ord=ord, X=x, Y=y, S=SxyStr,
                                                          pvalue = pvalue, pH0=pH0, pH1=pH1)
                                      }
                                    } else {
                                      ret <- data.frame(ord=ord, X=x, Y=y, S=SxyStr,
                                                        pvalue = pvalue)
                                    }
                                    if (save_files) {
                                      write.table(ret, file=csv_citestResults_file, sep=",",
                                                  row.names = FALSE, col.names = FALSE, append = TRUE)
                                    }
                                    futile.logger::flog.info(
                                      paste0("Finished citest: ", paste0(ret, collapse = "; ")),
                                      name = "citests_log")
                                    ret
                                  }
    citestResults <- rbind(done_citestResults, new_citestResults)
  }

  citestResults <- citestResults[order(citestResults$ord),]
  rownames(citestResults) <- NULL

  return(citestResults)
}

#' TODO: make citestResults an object with both the labels and the data.frame
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


#' @export getNumberCITests
getNumberCITests <- function(p, m.max=Inf, verbose=FALSE) {
  if (is.infinite(m.max)) {
    m.max = p-2
  }

  ntests = 0
  n_pairs <- ncol(mycombn(1:p, 2))
  n_tests_pair <- 0
  for (csetsize in 0:m.max) {
    ncs <- ncol(mycombn(1:(p-2), csetsize))
    n_tests_pair <- n_tests_pair + ncs
    if (verbose) {
      print(paste0("ord: ", csetsize, " -- ", ncs, " conditioning set(s)."))
    }
  }
  ntests <- n_tests_pair * n_pairs
  return(ntests)
}

# for (p in c(4,5,10,20)) {
#   print(getNumberCITests(p, verbose = FALSE))
# }


