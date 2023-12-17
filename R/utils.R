sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}


# triplet is an array of three nodes
#' @export isCollider
isCollider <- function(amat, triplet) {
  if (length(triplet) != 3) {
    print("A triplet must have exactly three nodes")
    return(NULL)
  }
  vi <- triplet[1]
  vm <- triplet[2]
  vj <- triplet[3]
  return(amat[vi, vm] == 2 && amat[vj, vm] == 2)
}

# triplet is an array of three nodes
#' @export isDefiniteNonCollider
isDefiniteNonCollider <- function(amat, triplet) {
  if (length(triplet) != 3) {
    print("A triplet must have exactly three nodes")
    return(NULL)
  }
  vi <- triplet[1]
  vm <- triplet[2]
  vj <- triplet[3]
  return ( (amat[vi, vm] == 3 || amat[vj, vm] == 3) ||  # vm has a tail
             (amat[vi, vm] == 1 && amat[vj, vm] == 1 && amat[vi, vj] == 0) ) # vi -o vm o- vj
}

#' if definite=TRUE, it only returns FALSE if the triplet is of definite status
#' and it is either a collider such that no descendants of it in in Z, or a
#' definite non-collider.
#' @importFrom pcalg searchAM
#' @export isMConTriplet
isMConTriplet <- function(amat, vi_name, vm_name, vj_name, Z,
                          definite=TRUE, verbose=FALSE) {
  vnames <- colnames(amat)
  connecting <- TRUE
  triplet <- c(vi_name, vm_name, vj_name)
  if (isCollider(amat, triplet)) {
    vm_id <- which(vnames == triplet[2])
    devm <- vnames[pcalg::searchAM(amat, vm_id, type="de")]
    if (length(which(devm %in% Z)) == 0) {
      # no descendant of this collider is in Z
      connecting <- FALSE
    }
  } else if (isDefiniteNonCollider(amat, triplet)) {
    if (length(which(Z == vm_name)) > 0) {
      # this non-collider is in Z
      connecting <- FALSE
    }
  } else {
    # triplet has not a definite status, so it is not definitely connecting if definite=TRUE
    # otherwise, can be either connecting or not connecting (i.e., output is TRUE)
    connecting <- !definite
  }
  return(connecting) # it only returns FALSE
}

# x is a node
#' @export getAdjNodes
getAdjNodes <- function(amat, x) {
  adjNodes <- names(which(amat[,x] != 0))
  return(adjNodes)
}

# TODO it was getMinimalSeparators
#' @importFrom dagitty impliedConditionalIndependencies
#' @export impliedSepsets
impliedSepsets <- function(magg) {
  sepset <- list()
  impliedCI <- dagitty::impliedConditionalIndependencies(magg)
  if (!is.null(impliedCI)) {
    sepset <- list()
    ci_ind = 1
    while (ci_ind <= length(impliedCI)) {
      x <- impliedCI[[ci_ind]]$X
      y <- impliedCI[[ci_ind]]$Y
      Sxy <- impliedCI[[ci_ind]]$Z

      if (is.null(sepset[[x]][[y]]))  {
        sepset[[x]][[y]] <- sepset[[y]][[x]] <- list()
        sepset[[x]][[y]][[1]] <- Sxy
        sepset[[y]][[x]][[1]] <- Sxy
      } else {
        nid <- length(sepset[[x]][[y]]) + 1
        sepset[[x]][[y]][[nid]] <- Sxy
        sepset[[y]][[x]][[nid]] <- Sxy
      }
      ci_ind <- ci_ind + 1
    }
  }
  return(sepset)
}

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


# type options are {"all", "all_shortest", "one_shortest"}
#' @export getMConnPaths
getMConnPaths <- function(amat, xnames, ynames, snames, definite=FALSE,
                          type="all", maxLength = Inf, verbose=FALSE) {
  connpaths <- list()
  for (x in xnames) {
    for (y in ynames) {
      pathList <- list(x)
      i = 1
      while (length(pathList) > 0) {
        #cat("New iteration...\n")
        #print(pathList)

        # there are still paths starting from X that are connecting and of definite status
        curpath <- pathList[[i]]
        len_curpath <- length(curpath)
        lastnode <- curpath[len_curpath]
        if (lastnode == y && len_curpath == 2 && type == "one_shortest") {
          return(curpath) # it is an edge X - Y
        }

        if (length(curpath) > 2) {
          vi_name <- curpath[length(curpath)-2]
          vm_name <- curpath[length(curpath)-1]
          vj_name <- lastnode
          if (!isMConTriplet(amat, vi_name, vm_name, vj_name, snames,
                             definite=definite, verbose=verbose)) {
            # removes the path from the list, as it will never be
            # a definite connecting path between X and Y given S
            pathList <- pathList[-i]
            next
          }
        }

        if (length(curpath) > maxLength) {
          break  # this goes to the next pair of x and y
        }

        # curpath starts with X and is definite connecting...
        if (lastnode == y) {
          if (verbose) {
            cat(paste0("connecting path between {",
                       paste0(xnames, collapse=","), "} and {",
                       paste0(ynames, collapse=","), "} given {",
                       paste0(snames, collapse=","), "} : ",
                       paste0(curpath, collapse=","), "\n"))
          }

          if (type == "one_shortest") {
            return(curpath)
          } else if (type == "all_shortest") {
            if (length(connpaths) == 0 ||
                length(curpath) == length(connpaths[[length(connpaths)]])) {
              # curpath is the first or of the same size than the previous paths
              connpaths[[length(connpaths)+1]] <- curpath
              pathList <- pathList[-i]
              next
            } else {
              break # this goes to the next pair of x and y
            }
          } else {
            connpaths[[length(connpaths)+1]] <- curpath
            pathList <- pathList[-i]
            next # this goes to the next path between x and y
          }
        }

        # not sure if this is a path to Y, so we continue traversing the path
        lastnode_i <- which(colnames(amat) == lastnode)
        adjnodes <- setdiff(getAdjNodes(amat, lastnode), curpath)
        n_adjn <- length(adjnodes)
        if (n_adjn == 0) {
          pathList <- pathList[-i] # path does not reach y
        } else {
          # checking if Y is one of the adj nodes
          if (y %in% adjnodes) {
            # replaces the current path with the one that has additionally the node Y
            pathList[[i]] <- c(curpath, y)
            adjnodes <- setdiff(adjnodes, y)
            n_adjn <- n_adjn - 1
          } else {
            pathList <- pathList[-i]
          }
          if (n_adjn > 0) {
            temp <- matrix(rep(curpath, n_adjn), ncol=n_adjn)
            temp <- rbind(temp, adjnodes)
            pathList <- append(pathList, split(temp, rep(1:ncol(temp), each = nrow(temp))))
          }
        }
      }
    }
  }
  # All paths between xnames and ynames are m-separated by snames
  return(connpaths)
}

# if TRUE, then X and Y are definitely m-separated by S
#' @export isMSeparated
isMSeparated <- function(amat, xnames, ynames, snames, verbose=FALSE) {
  # mhat separated (from old zhang's calculus) would have definite=FALSE
  connpaths <- getMConnPaths(amat, xnames, ynames, snames=snames,
                             type = "one_shortest",
                             definite=TRUE, verbose=verbose)
  return(is.null(connpaths) || length(connpaths) == 0)
}

isPossMConnected <- function(amat, xnames, ynames, snames,  verbose=FALSE) {
  connpaths <- getMConnPaths(amat, xnames, ynames, snames=snames,
                             definite=FALSE, verbose=verbose)
  return(!is.null(connpaths) && length(connpaths) > 0)
}

# type can be "pag" or "dag"
#' @importFrom pcalg pcalg2dagitty pag2magAM
#' @importFrom dagitty toMAG
#' @export getMAG
getMAG <- function(amat, type="pag") {
  if (type == "pag") {
    amat.mag <- pcalg::pag2magAM(amat, 1)
    #plotAG(amat.mag)
    magg <- pcalg::pcalg2dagitty(amat.mag, colnames(amat.mag), type="mag")
    #plot(magg)
  } else {
    adagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
    magg <- dagitty::toMAG(adagg)
    amat.mag <- dagitty2amat(magg, type="mag")
  }
  return(list(amat.mag = amat.mag, magg=magg))
}

#' @importFrom dagitty dseparated
#' @export dagittyCIOracle
dagittyCIOracle <- function(x, y, S, suffStat) {
  g <- suffStat$g
  labels <- names(g)
  if (dagitty::dseparated(g, labels[x], labels[y], labels[S])) {
    return(1)
  } else {
    return(0)
  }
}

# receives a dagitty g of type "dag" or "mag" and returns
# the true PAG as an pcalg fci object
#' @importFrom dagitty toMAG graphType
#' @importFrom pcalg fci
#' @export getTruePAG
getTruePAG <- function(g, verbose = FALSE) {
  indepTest <- dagittyCIOracle
  if (dagitty::graphType(g) == "dag") {
    g <- dagitty::toMAG(g)
  }
  suffStat <- list(g=g)
  truePag <- pcalg::fci(suffStat,
                        indepTest = indepTest,
                        labels= names(suffStat$g), alpha = 0.9999)
  return(truePag)
}

# returns an pcalg amat (adjacency matrix) of type amat.pag (same type for MAGs),
# where:
# 0: No edge
# 1: Circle
# 2: Arrowhead
# 3: Tail
#' @importFrom dagitty edges
#' @export dagitty2amat
dagitty2amat <- function(adagg, type="mag") {
  e <- NULL
  edg <- dagitty::edges(adagg)
  node_names <- names(adagg)
  ans_mat <- matrix(
    data = 0, nrow = length(node_names),
    ncol = length(node_names),
    dimnames = list(node_names, node_names)
  )

  diredg <- subset(edg, e == "->")

  ans_mat[as.matrix(diredg[c("w", "v")])] <- 3
  ans_mat[as.matrix(diredg[c("v", "w")])] <- 2

  bidiredg <-  subset(edg, e == "<->")
  ans_mat[as.matrix(bidiredg[c("w", "v")])] <- 2
  ans_mat[as.matrix(bidiredg[c("v", "w")])] <- 2

  return(ans_mat)
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


#' @importFrom igraph graph_from_adjacency_matrix
#' @export getIgraphMAG
getIgraphMAG <- function(amat.mag) {
  adjM <- amat.mag * 0
  for (i in 1:nrow(amat.mag)) {
    for (j in 1:ncol(amat.mag)) {
      if (amat.mag[i,j] == 2 && amat.mag[j,i] == 3) {
        # i -> j
        adjM[i,j] = 1
      } else if (amat.mag[i,j] == 3 && amat.mag[j,i] == 2) {
        # j -> i
        adjM[j,i] = 1
      } else if (amat.mag[i,j] == 2 && amat.mag[j,i] == 2) {
        # j -> i
        adjM[i,j] = adjM[j,i] = 1
      }
    }
  }
  return(igraph::graph_from_adjacency_matrix(adjM))
}


# types can be: png, pdf, or svg
#' @importFrom rsvg rsvg_png
#' @importFrom DOT dot
#' @export renderAG
renderAG <- function(amat, output_folder=NULL, fileid=NULL, type="png",
                     width=NULL, height=NULL, labels=NULL, add_index=TRUE) {
  if (is.null(labels)) {
    labels <- colnames(amat)
  }

  if (is.null(output_folder)) {
    output_folder = "./tmp/"
  }

  if (!file.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  if (is.null(fileid)) {
    fileid <- "ag"
  }

  if (type == "png") {
    png_filename <- paste0(output_folder, fileid, ".png")
    if (is.null(width)) {
      width = 1024
    }
    if (is.null(height)) {
      height = 1024
    }
  } else if (type == "pdf") {
    pdf_filename <- paste0(output_folder, fileid, ".pdf")
    if (is.null(width)) {
      width = 5
    }
    if (is.null(height)) {
      height = 5
    }
  }

  dot_filename <- paste0(output_folder, fileid, ".dot")
  svg_filename <- paste0(output_folder, fileid, ".svg")

  graphFile <- file(dot_filename, "w")
  cat('digraph graphname {', file=graphFile)
  cat('node [shape = oval];\n', file=graphFile)

  # For each row:
  for (i in 1:(nrow(amat)-1)) {
    # For each column:
    for (j in (i+1):ncol(amat)) {
      label_i <- labels[i]
      label_j <- labels[j]

      if (add_index) {
        label_i <- paste(label_i, "_", i, sep="")
        label_j <- paste(label_j, "_", j, sep="")
      }

      if (amat[i,j] > 0) {
        cat(label_i, "->", label_j, "[color=black, dir=both,",  file=graphFile)
        if (amat[i,j] == 1) {
          cat("arrowhead=odot, ", file=graphFile)
        } else if (amat[i,j] == 2) {
          cat("arrowhead=normal, ", file=graphFile)
        } else if (amat[i,j] == 3) {
          cat("arrowhead=none, ", file=graphFile)
        }

        if (amat[j,i] == 1) {
          cat("arrowtail=odot", file=graphFile)
        } else if (amat[j,i] == 2) {
          cat("arrowtail=normal", file=graphFile)
        } else if (amat[j,i] == 3) {
          cat("arrowtail=none", file=graphFile)
        }
        cat("];\n", file=graphFile)
      }
    }
  }
  cat("}\n", file=graphFile)
  close(graphFile)

  DOT::dot(paste(readLines(dot_filename), collapse=" "), file=svg_filename)

  if (type == "png") {
    rsvg::rsvg_png(svg_filename, png_filename, width = 1024, height = 1024)
  } else if (type == "pdf") {
    rsvg::rsvg_pdf(svg_filename, pdf_filename, width = 5, height = 5)
  }
}

getSepString <- function(S) {
  sepStr <- ""
  if (!is.null(S) && length(S) > 0) {
    sepStr <- paste0(S, collapse=",")
  }
  return(sepStr)
}


getSepVector <- function(sepStr) {
  if (is.null(sepStr)) {
    return(c())
  } else if (is.numeric(sepStr)) {
    return(sepStr)
  }
  return(as.numeric(unlist(strsplit(sepStr, ","))))
}


isValidPAG <- function(pagAdjM, verbose=FALSE) {
  amag <- getMAG(pagAdjM)
  ug_mag <- (amag$amat.mag == 3 & t(amag$amat.mag == 3)) * 1
  bg_mag <- (amag$amat.mag == 2 & t(amag$amat.mag == 2)) * 1
  dg_mag <- (amag$amat.mag == 2 & t(amag$amat.mag == 3)) * 1
  mag_ggm <- ggm::makeMG(dg_mag, ug_mag, bg_mag)
  isAGret <- tryCatch({
    ggm::isAG(mag_ggm)
  },
  error=function(cond) {
    print(cond)
    return(FALSE)
  },
  warning=function(cond) {
    print(cond)
    return(FALSE)
  })

  if (!isAGret) {
    if (verbose) {
      cat(paste("PAG is invalid -- canonical MAG is not ancestral.\n"))
    }
    return(FALSE)
  } else {
    # Here we check whether the PAG is valid by checking whether
    # we can perfectly recovery of the original PAG when
    # the canonical MAG is used as a C.I. oracle
    recPAG <- getTruePAG(amag$magg, verbose = TRUE)
    # plotAG(amag$amat.mag)
    # plot(recPAG)
    # plotAG(pagAdjM)
    if (any(recPAG@amat[colnames(pagAdjM), colnames(pagAdjM)] - pagAdjM != 0)) {
      if (verbose) {
        cat(paste("PAG is invalid! -- it is not the same as the MEC of its canonical MAG.\n"))
      }
      return(FALSE)
    }
    return(TRUE)
  }
}

#' @export getEdgeTypesList
getEdgeTypesList <- function(amat, edgeTypesList=NULL, labels=NULL) {
  if (is.null(edgeTypesList) || length(edgeTypesList) == 0) {
    edgeTypesList <- list()
  }
  if (is.null(labels)) {
    labels <- colnames(amat)
  }
  for (i in 1:(nrow(amat)-1)) {
    for (j in (i+1):ncol(amat)) {
      cur_entry <- paste0(labels[c(i,j)], collapse=",")
      if (is.null(edgeTypesList[[cur_entry]])) {
        edgeTypesList[[cur_entry]] <- list("0"=c(0), "edgeTypes"=matrix(0, 3, 3))
      }
      if (amat[i,j] == 0) {
        edgeTypesList[[cur_entry]][["0"]] <-
          edgeTypesList[[cur_entry]][["0"]] + 1
      } else {
        edgeTypesList[[cur_entry]]$edgeTypes[amat[j,i],amat[i,j]] <-
          edgeTypesList[[cur_entry]]$edgeTypes[amat[j,i],amat[i,j]] + 1
      }
    }
  }
  return(edgeTypesList)
}

#' @export summarizeEdgeTypesList
summarizeEdgeTypesList <- function(edgeTypesList) {
  relations <- names(edgeTypesList)
  dat_rels <- c()
  for (rel in relations) {
    cur_item <- edgeTypesList[[rel]]
    zero <-  as.numeric(cur_item$`0`)
    values <- as.numeric(cur_item$edgeTypes)
    cur_sum <- sum(as.numeric(cur_item$edgeTypes)) + zero
    dat_rels <- rbind(dat_rels, c(rel, zero, values, cur_sum))
  }

  dat_rels <- cbind.data.frame(dat_rels[,1], apply(dat_rels[,-1], 2, as.numeric))
  colnames(dat_rels) <- c("relation", "0", "o-o", "<-o", "-o",
                          "o->", "<->", "->",
                          "o-", "<-", "-", "sum")
  return(dat_rels)
}

#' @export getPAGImpliedSepset
getPAGImpliedSepset <- function(amat, amat.type="pag", ret_df=FALSE) {
  mag_out <- getMAG(amat)
  magg <- mag_out$magg
  return(getImpliedSepset(magg, labels, ret_df = ret_df))
}

# This returns the sepset given a dagitty MAG
getImpliedSepset <- function(magg, labels, ret_df=FALSE) {
  p <- length(labels)
  seq_p <- seq_len(p) # creates a sequence of integers from 1 to p
  sepset <- lapply(seq_p, function(.) vector("list",p)) # a list of lists [p x p]

  impliedCI <- dagitty::impliedConditionalIndependencies(magg)
  if (!is.null(impliedCI)) {
    ci_ind = 1
    while (ci_ind <= length(impliedCI)) {
      x <- which(labels %in% impliedCI[[ci_ind]]$X)
      y <- which(labels %in% impliedCI[[ci_ind]]$Y)
      Sxy <- which(labels %in% impliedCI[[ci_ind]]$Z)

      sepset[[x]][[y]] <- sepset[[y]][[x]] <- Sxy

      ci_ind <- ci_ind + 1
    }
  }

  if (!ret_df) {
    return(sepset)
  } else {
    sepset_df <- c()
    for (i in 1:length(sepset)) {
      for (j in i:length(sepset)) {
        if (!is.null(sepset[[i]][[j]])) {
          cur_ord <- length(sepset[[i]][[j]])
          sepset_df <- rbind.data.frame(sepset_df,
                                        c(ord=cur_ord, X=i, Y=j,
                                          S=getSepString(sepset[[i]][[j]])))
        }
      }
    }

    if (!is.null(sepset_df)) {
      colnames(sepset_df) <- c("ord", "X", "Y", "S")
    }

    return(sepset_df)
  }
}
