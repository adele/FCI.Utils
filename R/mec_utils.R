getZeroOrderTriplets <- function(amat, dupl=FALSE, verbose=FALSE) {
  X = Y = Z = NULL

  ind <- which(amat > 0, arr.ind = TRUE)
  labels <- colnames(amat)
  c0_trplts <- c()
  nc0_trplts <- c()
  for (i in seq_len(nrow(ind))) {
    x <- ind[i, 1]
    z <- ind[i, 2]
    allY <- setdiff(which(amat[z, ] != 0), x) # neighbors of Z, except X
    for (y in allY) {
      if (!dupl) {
        if (!is.null(c0_trplts) &&
            (nrow(subset(c0_trplts, X == x & Y == y & Z == z)) > 0 ||
             nrow(subset(c0_trplts, X == y & Y == x & Z == z)) > 0))
          next
        if (!is.null(nc0_trplts) &&
            (nrow(subset(nc0_trplts, X == x & Y == y & Z == z)) > 0 ||
             nrow(subset(nc0_trplts, X == y & Y == x & Z == z)) > 0))
          next
      }
      if (amat[x, y] == 0) { #X and Y are unshielded
        if (amat[x,z] == 2 && amat[y,z] == 2) {
          if (verbose) cat("\n Adding c0 triplet ",labels[x],"->", labels[z], "<-", labels[y], "\n")
          c0_trplts <- rbind(c0_trplts, data.frame("X"=x, "Z"=z, "Y"=y))
        } else {
          if (verbose) cat("\n Adding nc0 triplet ",labels[x],"-", labels[z], "-", labels[y], "\n")
          nc0_trplts <- rbind(nc0_trplts, data.frame("X"=x, "Z"=z, "Y"=y))
        }
      }
    }
  }
  return(list(C0=c0_trplts, NC0=nc0_trplts))
}

# returns the skeleton, and lists of colliders and non-colliders with order
#' @importFrom FCI.Utils getSepString
#' @importFrom jsonlite toJSON fromJSON
#' @export MAGtoMEC
MAGtoMEC <- function(amat.mag, verbose=FALSE) {
  X = Z = NULL
  skel <- 1 * (amat.mag > 0)
  cnc0 <- getZeroOrderTriplets(amat.mag, dupl = T, verbose)

  processList <- list()
  ncK_trplts <- c()
  if (!is.null(cnc0$NC0)) {
    path <- apply(cnc0$NC0, 1, getSepString)
    ncK_trplts <- cbind(cnc0$NC0, "X0"=cnc0$NC0$X, "Y0"=cnc0$NC0$Y,
                        ord=0, "SepOrd"=getSepString(c()), "path"=path)
    nc0.i <- 1
    while (nc0.i <= nrow(cnc0$NC0)) {
      x <- cnc0$NC0[nc0.i,"X"]
      z <- cnc0$NC0[nc0.i,"Z"]
      y <- cnc0$NC0[nc0.i,"Y"]
      x0 <- x
      y0 <- y
      curord <- 0

      if (!is.null(cnc0$C0)) {
        i = 1
        potQ <- subset(cnc0$C0, X == x & Z == z)
        if (nrow(potQ) > 0) {
          q.i <- 1
          while (q.i <= nrow(potQ)) {
            q <- potQ[q.i,"Y"]
            if (amat.mag[y, q] > 0) {
              processList[[length(processList)+1]] <- list("X"=z, "Z"=q, "Y"=y,
                                                           "X0"=x0, "Y0"=y0,
                                                           "ord"=curord+1,
                                                           "SepOrd"=c(z),
                                                           path=c(x,z,q,y))
              i = i + 1
            }
            q.i <- q.i + 1
          }
        }
      }
      nc0.i <- nc0.i + 1
    }
  }

  cK_trplts <- c()
  if (!is.null(cnc0$C0)) {
    path <- apply(cnc0$C0, 1, getSepString)
    cK_trplts <- cbind(cnc0$C0, "X0"=cnc0$C0$X, "Y0"=cnc0$C0$Y,
                       ord=0, "SepOrd"=getSepString(c()), "path"=path)
  }

  while (length(processList) > 0) {
    if (verbose) {
      cat("length of processList: ", length(processList), "\n")
    }
    trplt <- processList[[1]]
    x <- trplt$X
    z <- trplt$Z
    y <- trplt$Y
    x0 <- trplt$X0
    y0 <- trplt$Y0
    curpath <- trplt$path
    curord <- trplt$ord
    sepOrd <- trplt$SepOrd # a vector

    processList <- processList[-1]
    if (amat.mag[x,z] == 2 && amat.mag[y,z] == 2) {
      if (verbose) {
        cat("collider triplet", toJSON(trplt), "\n")
      }

      # <X,Z,Y> is a collider along the discriminanting path between X0 and Y (=Y0)
      cK_trplts <- rbind(cK_trplts, data.frame("X"=x, "Z"=z, "Y"=y, "X0"=x0, "Y0"=y0,
                                               "ord"=curord,
                                               "SepOrd"=getSepString(sepOrd),
                                               "path"=getSepString(curpath)))
      potQ <- rbind(subset(ncK_trplts, X == x & Z == z),
                    subset(ncK_trplts, X == y & Z == z))  # x = y and y = x
      if (nrow(potQ) > 0) {
        q.i <- 1
        curX0 = x0
        curY0 = y0
        cur_curord = curord
        cur_sepOrd = sepOrd
        cur_curpath = curpath
        while (q.i <= nrow(potQ)) {
          x0 = curX0
          y0 = curY0
          curord = cur_curord
          sepOrd = cur_sepOrd
          curpath = cur_curpath
          if (x != potQ[q.i,"X"]) {
            ycopy = y
            y = x
            x = ycopy
          }
          q <- potQ[q.i,"Y"] # <X, Z, Q> is a non-collider
          if (amat.mag[y, q] > 0 && !(q %in% curpath)) {
            y0 <- q
            potX0_ids <- curpath # [1:which(curpath==z)]
            sep_X0s <- which(amat.mag[potX0_ids, y0] == 0)
            if (length(sep_X0s) > 0) {
              potX0_ids <- potX0_ids[sep_X0s]
              if (length(potX0_ids) > 1) {
                potX0_ids = potX0_ids[length(potX0_ids)]
              }
              x0 <- potX0_ids
              curpath <- c(curpath[which(curpath == x0):which(curpath==z)], y, q)
              curord <- length(curpath) - 3
              sepOrd <- curpath[2:((length(curpath)-2))]

              # if (amat.mag[x, q] == 0) {
              #   y0 <- q
              #   x0 <- x #  <X, Z, Q> is the first non-collider out of the discriminating path <X, Z, Y, Q>
              #   curord = 0
              #   sepOrd = c()
              #   curpath <- curpath[which(curpath == x0):length(curpath)]
              # }
              #processList[[length(processList)+1]] <- c(z, y, q, "X0"=x0, "Y0"=y0)

              processList[[length(processList)+1]] <- #list("X"=z, "Z"=q, "Y"=y,
                list("X"=z, "Z"=y, "Y"=q,
                     "X0"=x0, "Y0"=y0,
                     "ord"= curord,# curord+1,
                     "SepOrd"= sepOrd, #c(sepOrd, z),
                     "path"= curpath) #c(curpath[1:which(curpath==z)], y, q))
              #"path"= c(curpath[1:which(curpath==z)], q, y))
            }
          }
          q.i <- q.i + 1
        }
      }
    } else {
      if (verbose) {
        cat("non-collider triplet", toJSON(trplt), "\n")
      }
      ncK_trplts <- rbind(ncK_trplts,
                          data.frame("X"=x, "Z"=z, "Y"=y, "X0"=x0, "Y0"=y0,
                                     "ord"=curord, "SepOrd"=getSepString(sepOrd),
                                     "path"=getSepString(curpath)))
      potQ <- subset(cK_trplts, X == x & Z == z)
      if (nrow(potQ) > 0) {
        q.i <- 1
        while (q.i <= nrow(potQ)) {
          q <- potQ[q.i,"Y"]
          if (amat.mag[y, q] > 0 && !(q %in% curpath)) {
            processList[[length(processList)+1]] <- list("X"=z, "Z"=q, "Y"=y,
                                                         "X0"=x0, "Y0"=y0,
                                                         "ord"=curord+1,
                                                         "SepOrd"=c(sepOrd, z),
                                                         "path"= c(curpath[1:which(curpath==z)], q, y))
          }
          q.i <- q.i + 1
        }
      }
    }
  }

  # remove duplicated paths
  # NOTE: each line in this matrix correspond to a path between X0 and Y0, in which the variables
  # in SepOrd are colliders in the beginning of the path such that, when condition on, makes the
  # variable in z being either a collider (if in the ck_trplts) or a non-collider (if in the nck_trplts)
  if (!is.null(cK_trplts)) {
    ck_trplts_str <- apply(cK_trplts[,c("X", "Z", "Y", "X0", "Y0")], 1, function(x) {
      paste0(c(min(x[1],x[3]), x[2], max(x[1],x[3]),
               min(c(x[4], x[5])), max(c(x[4], x[5]))), collapse=",")})
    cK_trplts <- cK_trplts[!duplicated(ck_trplts_str),]
    cK_trplts <- cK_trplts[order(cK_trplts$X0, cK_trplts$Y0), ]
  }

  if (!is.null(ncK_trplts)) {
    ncK_trplts_str <- apply(ncK_trplts[,c("X", "Z", "Y", "X0", "Y0")], 1, function(x) {
      paste0(c(min(x[1],x[3]), x[2], max(x[1],x[3]),
               min(c(x[4], x[5])), max(c(x[4], x[5]))), collapse=",")})
    ncK_trplts <- ncK_trplts[!duplicated(ncK_trplts_str),]
    ncK_trplts <- ncK_trplts[order(ncK_trplts$X0, ncK_trplts$Y0), ]
  }
  return(list(skel=skel, CK=cK_trplts, NCK=ncK_trplts))
}


#' @importFrom jsonlite toJSON
#' @export updateMECTripletIsolators
updateMECTripletIsolators <- function(mec, amat.ag, verbose=FALSE) {
  ncK_trplts <- mec$NCK
  cK_trplts <- mec$CK
  skel <- mec$skel

  labels <- colnames(amat.ag)

  if (!is.null(cK_trplts)) {
    for (i in 1:nrow(cK_trplts)) {
      x = cK_trplts[i, "X0"]
      y = cK_trplts[i, "Y0"]
      z = cK_trplts[i, "Z"]
      pathvars = getSepVector(cK_trplts[i, "path"])
      if (length(pathvars) > 0) {
        pathvars = labels[pathvars]
      }
      ignore_path_list=list(pathvars)
      sepset_out = getImpliedConditionalSepset(amat.ag, labels[x], labels[y], labels[z],
                                               definite=TRUE, ignore_path_list,
                                               ignore_sepvar_names = NULL,
                                               verbose=verbose)
      done = FALSE
      if (!any(is.na(sepset_out$sepset))) {
        connpaths <- getMConnPaths(amat.ag, labels[x], labels[y], sepset_out$sepset, definite=TRUE)
        nconnpaths <- length(connpaths)
        if (nconnpaths > 0) {
          if (any(sapply(connpaths, function(x) {
            length(x) == length(pathvars) && all(pathvars %in% x) }))) {
            cur_sepset <- sepset_out$sepset
            if (length(cur_sepset) > 0) {
              cur_sepset <- which(labels %in% sepset_out$sepset)
            }
            cK_trplts[i, "dep_sepset"] <-  getSepString(cur_sepset)
            cK_trplts[i, "nconnpaths"] <- nconnpaths
            cK_trplts[i, "connpaths"] <- toJSON(connpaths)
            done = TRUE
          }
        }
      }
      if (!done) {
        stop("There is a problem with MEC separators!")
      }
    }
  }

  if (!is.null(ncK_trplts)) {
    for (i in 1:nrow(ncK_trplts)) {
      x = ncK_trplts[i, "X0"]
      y = ncK_trplts[i, "Y0"]
      z = ncK_trplts[i, "Z"]
      pathvars = getSepVector(ncK_trplts[i, "path"])
      if (length(pathvars) > 0) {
        pathvars = labels[pathvars]
      }
      ignore_path_list=list(pathvars)
      ignore_sepvar_names <- labels[z]
      sepset_out = getImpliedConditionalSepset(amat.ag, labels[x], labels[y], NULL,
                                               definite=TRUE, ignore_path_list, ignore_sepvar_names,
                                               verbose=verbose) #TODO
      done = FALSE
      if (!any(is.na(sepset_out$sepset))) {
        connpaths <- getMConnPaths(amat.ag, labels[x], labels[y], sepset_out$sepset, definite=TRUE)
        nconnpaths <- length(connpaths)
        if (nconnpaths > 0) {
          if (any(sapply(connpaths, function(x) {
            length(x) == length(pathvars) && all(pathvars %in% x) }))) {
            cur_sepset <- sepset_out$sepset
            if (length(cur_sepset) > 0) {
              cur_sepset <- which(labels %in% sepset_out$sepset)
            }
            ncK_trplts[i, "dep_sepset"] <-  getSepString(cur_sepset)
            ncK_trplts[i, "connpaths"] <- toJSON(connpaths)
            ncK_trplts[i, "nconnpaths"] <- nconnpaths
            done = TRUE
          }
        }
      }
      if (!done) {
        stop("There is a problem with MEC separators for NC!")
      }
    }
  }

  return(list(skel=skel, CK=cK_trplts, NCK=ncK_trplts))
}
