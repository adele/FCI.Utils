# Generates obs. dataset following a linear SEM, compatible with a dagitty DAG, adag
# type: defines the type of the variables, "continuous", "binary", or "mixed"
# f.args: a list indexed by the names of the variables, where each entry is another list
# with an entry names levels indicating the number of levels a discrete node takes or
# levels = 1 for continuous variable
#' @importFrom dagitty simulateLogistic simulateSEM localTests
#' @export generateDataset
generateDataset <- function(adag, N, type="continuous", verbose=FALSE,
                            f.args = NULL, coef_thresh = 0.2) {
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
        if(type == "binary") {
          obs.dat <- dagitty::simulateLogistic(adag, N=N, verbose=verbose)
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
          obs.dat <- dagitty::simulateSEM(adag, N=N)
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

# a path with 1 bidirected edge
# A -> B <-> C <- D
getDAG1BE <- function() {
  allvars <- c("A", "B", "C", "D", "Ubc")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["A","B"] <- 0; amat["B","A"] <- 1; # A -> B
  amat["D","C"] <- 0; amat["C","D"] <- 1; # D -> C
  amat["Ubc","B"] <- 0; amat["B","Ubc"] <- 1; # Ubc -> B
  amat["Ubc","C"] <- 0; amat["C","Ubc"] <- 1; # Ubc -> C

  lat <- c("Ubc")
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}


# a path with 4 bidirected edged
# A <-> B <-> C <-> D; A <-> D
getDAG4BE <- function() {
  allvars <- c("A", "B", "C", "D", "Uab", "Ubc", "Ucd", "Uad")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["Uab","A"] <- 0; amat["A","Uab"] <- 1; # Uab -> A
  amat["Uab","B"] <- 0; amat["B","Uab"] <- 1; # Uab -> B
  amat["Ubc","B"] <- 0; amat["B","Ubc"] <- 1; # Ubc -> B
  amat["Ubc","C"] <- 0; amat["C","Ubc"] <- 1; # Ubc -> C
  amat["Ucd","C"] <- 0; amat["C","Ucd"] <- 1; # Ucd -> C
  amat["Ucd","D"] <- 0; amat["D","Ucd"] <- 1; # Ucd -> D
  amat["Uad","A"] <- 0; amat["A","Uad"] <- 1; # Uad -> A
  amat["Uad","D"] <- 0; amat["D","Uad"] <- 1; # Uad -> D

  lat <- c("Uab", "Ubc", "Ucd", "Uad")
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}


# An edge Vi -> Vj is represented by a 0 in (Vi,Vj) and a 1 in (Vj, Vi)
getDAGPdSep <- function() {
  allvars <- c("X", "Z", "W", "Y", "Uzw")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["X","Z"] <- 0; amat["Z","X"] <- 1; # x -> z
  amat["Z","Y"] <- 0; amat["Y","Z"] <- 1; # z -> y
  amat["W","Y"] <- 0; amat["Y","W"] <- 1; # w -> y
  amat["Uzw","Z"] <- 0; amat["Z","Uzw"] <- 1; # uzw -> z
  amat["Uzw","W"] <- 0; amat["W","Uzw"] <- 1; # uzw -> w

  lat <- c("Uzw")
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}


# A -> B -> C; B <- D -> C
getDAGIV <- function() {
  allvars <- c("A", "B", "C", "D")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["A","B"] <- 0; amat["B","A"] <- 1; # a -> b
  amat["B","C"] <- 0; amat["C","B"] <- 1; # b -> c
  amat["D","B"] <- 0; amat["B","D"] <- 1; # d -> b
  amat["D","C"] <- 0; amat["C","D"] <- 1; # d -> c

  lat <- c()
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}


getDAG3Anc <- function() {
  allvars <- c("A", "B", "C", "D", "E", "Ubd", "Ucd")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["E","A"] <- 1; amat["A","E"] <- 0; # A -> E
  amat["B","A"] <- 1; amat["A","B"] <- 0; # A -> B
  amat["E","B"] <- 1; amat["B","E"] <- 0; # B -> E
  amat["C","B"] <- 1; amat["B","C"] <- 0; # B -> C
  amat["E","D"] <- 1; amat["D","E"] <- 0; # D -> E

  amat["D","Ubd"] <- 1; amat["Ubd","D"] <- 1; # Ubd -> D
  amat["B","Ubd"] <- 1; amat["Ubd","B"] <- 1; # Ubd -> B
  amat["D","Ucd"] <- 1; amat["Ucd","D"] <- 1; # Ucd -> D
  amat["C","Ucd"] <- 1; amat["Ucd","C"] <- 1; # Ucd -> C

  lat <- c("Ubd", "Ucd")
  dag <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dag) <- lat

  return(list(amat=amat, lat=lat, dagg=dag))
}



# D -> C <- B; A <- C
getDAG3Colliders <- function() {
  allvars <- c("A", "B", "C", "D")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["A","C"] <- 0; amat["C","A"] <- 1; # a -> c
  amat["B","C"] <- 0; amat["C","B"] <- 1; # b -> c
  amat["D","C"] <- 0; amat["C","D"] <- 1; # d -> c

  lat <- c()
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}


# D -> C <- B; A <- C -> B
getDAGDescCollider <- function() {
  allvars <- c("A", "B", "C", "D")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["A","C"] <- 0; amat["C","A"] <- 1; # a -> c
  amat["B","C"] <- 0; amat["C","B"] <- 1; # b -> c
  amat["C","D"] <- 0; amat["D","C"] <- 1; # c -> d

  lat <- c()
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}

# D -> C <- E; A <- C -> B
getDAGTwoDescCollider <- function() {
  allvars <- c("A", "B", "C", "D", "E")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["C","A"] <- 0; amat["A","C"] <- 1; # c -> a
  amat["C","B"] <- 0; amat["B","C"] <- 1; # c -> b
  amat["D","C"] <- 0; amat["C","D"] <- 1; # d -> c
  amat["E","C"] <- 0; amat["C","E"] <- 1; # e -> c

  lat <- c()
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}

# X -> Z <- W -> Y; X <- UXW -> W; Z <- UZY -> Y
getDAGCollFork <- function() {
  allvars <- c("W", "X", "Y", "Z", "Uxw", "Uzy")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["X","Z"] <- 0; amat["Z","X"] <- 1; # x -> z
  amat["W","Z"] <- 0; amat["Z","W"] <- 1; # w -> z
  amat["W","Y"] <- 0; amat["Y","W"] <- 1; # w -> y
  amat["Uzy","Z"] <- 0; amat["Z","Uzy"] <- 1; # uzy -> z
  amat["Uzy","Y"] <- 0; amat["Y","Uzy"] <- 1; # uzy -> y
  amat["Uxw","X"] <- 0; amat["X","Uxw"] <- 1; # uxw -> x
  amat["Uxw","W"] <- 0; amat["W","Uxw"] <- 1; # uxw -> w

  lat <- c("Uxw", "Uzy")
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}




# A -> B -> C -> D
getDAGChain <- function() {
  allvars <- c("A", "B", "C", "D")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["A","B"] <- 0; amat["B","A"] <- 1; # a -> b
  amat["B","C"] <- 0; amat["C","B"] <- 1; # b -> c
  amat["C","D"] <- 0; amat["D","C"] <- 1; # c -> d

  lat <- c()
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}



getDAGTriplet <- function(collider=TRUE) {
  allvars <- c("X", "Y", "Z")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["X","Z"] <- 0; amat["Z","X"] <- 1; # x -> z

  if (collider) {
    amat["Y","Z"] <- 0; amat["Z","Y"] <- 1; # y -> z
  } else {
    amat["Z","Y"] <- 0; amat["Y","Z"] <- 1; # z -> y
  }

  lat <- c()
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}



getDAGDiscrPaths2 <- function() {
  allvars <- c("X", "Y", "A", "B","C", "D","E", "Z",
               "Uab", "Ubc", "Uad", "Ude", "Uey")

  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars


  amat["X","Z"] <- 0; amat["Z","X"] <- 1; # x -> z
  amat["Z","Y"] <- 0; amat["Y","Z"] <- 1; # z -> y

  amat["X","A"] <- 0; amat["A","X"] <- 1; # x -> a
  amat["A","Y"] <- 0; amat["Y","A"] <- 1; # a -> y
  amat["Uab","A"] <- 0; amat["A","Uab"] <- 1; # uab -> a
  amat["Uab","B"] <- 0; amat["B","Uab"] <- 1; # uab -> b
  amat["Ubc","B"] <- 0; amat["B","Ubc"] <- 1; # ubc -> b
  amat["Ubc","C"] <- 0; amat["C","Ubc"] <- 1; # ubc -> c
  amat["B","Y"] <- 0; amat["Y","B"] <- 1; # b -> y
  amat["C","Y"] <- 0; amat["Y","C"] <- 1; # c -> y

  amat["Uad","A"] <- 0; amat["A","Uad"] <- 1; # uad -> a
  amat["Uad","D"] <- 0; amat["D","Uad"] <- 1; # uad -> d
  amat["Ude","D"] <- 0; amat["D","Ude"] <- 1; # ude -> d
  amat["Ude","E"] <- 0; amat["E","Ude"] <- 1; # ude -> e
  amat["D","Y"] <- 0; amat["Y","D"] <- 1; # d -> y
  amat["Uey","E"] <- 0; amat["E","Uey"] <- 1; # uey -> e
  amat["Uey","Y"] <- 0; amat["Y","Uey"] <- 1; # uey -> y


  lat <- c("Uab", "Ubc", "Uad", "Ude", "Uey")

  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}

# A - > X; B -> X -> Y; X <- C <-> Y;
getDAG2IVs <- function() {
  allvars <- c("X", "Y", "A", "B", "C", "Ucy")

  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars

  amat["A","X"] <- 0; amat["X","A"] <- 1; # A -> X
  amat["B","X"] <- 0; amat["X","B"] <- 1; # B -> X
  amat["C","X"] <- 0; amat["X","C"] <- 1; # C -> X
  amat["X","Y"] <- 0; amat["Y","X"] <- 1; # b -> y

  amat["Ucy","C"] <- 0; amat["C","Ucy"] <- 1; # Ucy -> C
  amat["Ucy","Y"] <- 0; amat["Y","Ucy"] <- 1; # Ucy -> Y

  lat <- c("Ucy")

  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}


# W -> B <-> A <- X; B -> Y; A -> Y
getDAGIV2 <- function() {
  allvars <- c("W", "X", "Y", "A", "B", "Uab")

  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars

  amat["W","B"] <- 0; amat["B","W"] <- 1; # w -> b
  amat["X","A"] <- 0; amat["A","X"] <- 1; # x -> a
  amat["A","Y"] <- 0; amat["Y","A"] <- 1; # a -> y
  amat["Uab","A"] <- 0; amat["A","Uab"] <- 1; # uab -> a
  amat["Uab","B"] <- 0; amat["B","Uab"] <- 1; # uab -> b
  amat["B","Y"] <- 0; amat["Y","B"] <- 1; # b -> y

  lat <- c("Uab")

  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}

getDAGCollFork2 <- function() {
  allvars <- c("A", "B", "C", "D", "E", "Uab", "Uad")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars

  # A <-> B
  amat["Uab","A"] <- 0; amat["A","Uab"] <- 1; # Uab -> A
  amat["Uab","B"] <- 0; amat["B","Uab"] <- 1; # Uab -> B

  amat["C","B"] <- 0; amat["B","C"] <- 1; # C -> B
  amat["C","E"] <- 0; amat["E","C"] <- 1; # C -> E
  amat["B","D"] <- 0; amat["D","B"] <- 1; # B -> D
  amat["B","E"] <- 0; amat["E","B"] <- 1; # B -> E

  # A <-> D
  amat["Uad","A"] <- 0; amat["A","Uad"] <- 1; # Uad -> A
  amat["Uad","D"] <- 0; amat["D","Uad"] <- 1; # Uad -> D

  lat <- c("Uab", "Uad")
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
}


getDAG2DiscrPaths <- function() {
  allvars <- c("A", "B", "C", "D", "E", "Ubc", "Uce")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars

  amat["A","B"] <- 0; amat["B","A"] <- 1; # A -> B
  amat["B","D"] <- 0; amat["D","B"] <- 1; # B -> D
  amat["B","E"] <- 0; amat["E","B"] <- 1; # B -> E
  amat["C","D"] <- 0; amat["D","C"] <- 1; # C -> D

  # B <-> C
  amat["Ubc","B"] <- 0; amat["B","Ubc"] <- 1; # Ubc -> B
  amat["Ubc","C"] <- 0; amat["C","Ubc"] <- 1; # Ubc -> C

  # C <-> E
  amat["Uce","C"] <- 0; amat["C","Uce"] <- 1; # Uce -> C
  amat["Uce","E"] <- 0; amat["E","Uce"] <- 1; # Uce -> E

  lat <- c("Ubc", "Uce")
  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat
  dagitty::coordinates(dagg) <-
    list( x=c(A=1, B=1, C=1, D=0, E= 2, Ubc=1, Uce=1.5),
          y=c(A=0, B=1, C=2, D=2.5, E=2.5, Ubc=1.5, Uce=2.25) )

  return(list(amat=amat, lat=lat, dagg=dagg))
}



# A discriminating path between X and Y
# if the path is discriminating for "B",
#    if collider == TRUE, then
#       X -> Y; A -> Y; X -> A <-> B <-> Y
#    else,
#       X -> Y; A -> Y; X -> A <-> B -> Y
# if the path is discriminating for "C",
#    if collider == TRUE, then
#       X -> Y; A -> Y; B -> Y, X -> A <-> B <-> C <-> Y
#    else,
#       X -> Y; A -> Y; B -> Y, X -> A <-> B <-> C -> Y
getDAGDiscrPath <- function(collider = TRUE, discr_var="B") {
  allvars <- c("X", "Y", "A", "B", "Uab")
  if (discr_var == "B" && collider) {
    allvars <- c(allvars, "Uby")
  }
  if (discr_var == "C") {
    allvars <- c(allvars, c("C", "Ubc"))
    if (collider) {
      allvars <- c(allvars, "Ucy")
    }
  }

  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars

  amat["X","A"] <- 0; amat["A","X"] <- 1; # x -> a
  amat["A","Y"] <- 0; amat["Y","A"] <- 1; # a -> y
  amat["Uab","A"] <- 0; amat["A","Uab"] <- 1; # uab -> a
  amat["Uab","B"] <- 0; amat["B","Uab"] <- 1; # uab -> b

  if (discr_var == "C") {
    amat["Ubc","B"] <- 0; amat["B","Ubc"] <- 1; # ubc -> b
    amat["Ubc","C"] <- 0; amat["C","Ubc"] <- 1; # ubc -> c
  }

  if (collider && discr_var == "B") {
    amat["Uby","B"] <- 0; amat["B","Uby"] <- 1; # uby -> b
    amat["Uby","Y"] <- 0; amat["Y","Uby"] <- 1; # uby -> y
  }

  if ((!collider && discr_var == "B") || (discr_var == "C")) {
    amat["B","Y"] <- 0; amat["Y","B"] <- 1; # b -> y
  }

  if (collider && discr_var == "C") {
    amat["Ucy","C"] <- 0; amat["C","Ucy"] <- 1; # ucy -> c
    amat["Ucy","Y"] <- 0; amat["Y","Ucy"] <- 1; # ucy -> y
  }

  if (!collider && discr_var == "C") {
    amat["C","Y"] <- 0; amat["Y","C"] <- 1; # c -> y
  }

  lat <- c("Uab")
  if (collider && discr_var == "B") {
    lat <- c(lat, "Uby")
  }

  if (discr_var == "C") {
    lat <- c(lat, "Ubc")
    if (collider) {
      lat <- c(lat, "Ucy")
    }
  }

  dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(dagg) <- lat

  return(list(amat=amat, lat=lat, dagg=dagg))
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


# dag types:
#    fork: X <- Z -> Y
#    collider: X -> Z <- Y
#    chain4: X -> Z -> Y -> W
#    iv: A -> B -> C; B <- D -> C
#    collfork: X -> Z <- W -> Y; X <-> W; Z <-> Y
#    iv2: X -> Y; A -> Y; X -> A <-> B <-> Y
#    discr1_c: X -> Y; A -> Y; X -> A <-> B <-> Y
#    discr1_nc: X -> Y; A -> Y; X -> A <-> B -> Y
#    discr2_c: X -> Y; A -> Y; B->Y; X -> A <-> B <-> C <-> Y
#    discr2_nc: X -> Y; A -> Y; B->Y; X -> A <-> B <-> C -> Y
#.   2descColl: D -> C <- E; A <- C -> B
#    descColl: A -> C <- B; D <- C
#' @export getDAG
getDAG <- function(type="fork") {
  if (type == "fork") {
    return(getDAGTriplet(collider = FALSE))
  } else if (type == "collider") {
    return(getDAGTriplet(collider = TRUE))
  } else if (type == "chain4") {
    return(getDAGChain())
  } else if (type == "iv") {
    return(getDAGIV())
  } else if (type == "iv2") {
      return(getDAGIV2())
  } else if (type == "2ivs") {
      return(getDAG2IVs())
  } else if (type == "collfork") {
    return(getDAGCollFork())
  } else if (type == "collfork2") {
    return(getDAGCollFork2())
  } else if (type == "3Colls") {
    return(getDAG3Colliders())
  } else if (type == "2descColl") {
    return(getDAGTwoDescCollider())
  } else if (type == "descColl") {
    return(getDAGDescCollider())
  } else if (type == "discr1_c") {
    return(getDAGDiscrPath(collider = TRUE, discr_var = "B"))
  } else if (type == "discr1_nc") {
    return(getDAGDiscrPath(collider = FALSE, discr_var = "B"))
  } else if (type == "discr2_c") {
    return(getDAGDiscrPath(collider = TRUE, discr_var = "C"))
  } else if (type == "discr2_nc") {
    return(getDAGDiscrPath(collider = FALSE, discr_var = "C"))
  } else if (type == "2discrs_nc") {
    return(getDAG2DiscrPaths())
  } else if (type == "3anc") {
    return(getDAG3Anc())
  } else if (type == "pdsep_g") { # same as iv
    return(getDAGPdSep())
  } else if (type == "1be") {
    return(getDAG1BE())
  } else if (type == "4be") {
    return(getDAG4BE())
  }
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
    min_bscore = min(c(subset(f_citestResults, type == "indep", select = pH0, drop=TRUE),
                       subset(f_citestResults, type == "dep", select = pH1, drop=TRUE)))
    faithful_bprop = length(which(f_citestResults$bf)) / length(f_citestResults$bf)
    ret <- list(min_bscore = min_bscore,
                f_citestResults = f_citestResults,
                faithful_bprop = faithful_bprop,
                faithful_pprop = faithful_pprop)
  } else {
    ret <- list(f_citestResults = f_citestResults,
                faithful_pprop = faithful_pprop)

  }
  return(ret)
}
