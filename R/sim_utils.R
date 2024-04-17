# Generates obs. dataset following a linear SEM, compatible with a dagitty DAG, adag
# Type defines the type of the variables, which can be either "continuous" or "binary"
#' @importFrom dagitty simulateLogistic simulateSEM localTests
#' @export generateDataset
generateDataset <- function(adag, N, type="continuous", verbose=FALSE) {
  if (!(type %in% c("continuous", "binary")))  {
    stop("type must be either continuous or binary")
  }

  done <- FALSE
  while (!done) {
    done <- tryCatch(
      {
        if(type == "binary") {
          obs.dat <- dagitty::simulateLogistic(adag, N=N, verbose=FALSE)
          obs.dat <- as.data.frame(sapply(obs.dat, function(col) as.numeric(col)-1))
          lt <- dagitty::localTests(adag, obs.dat, type="cis.chisq")
        } else if (type == "continuous") {
          obs.dat <- dagitty::simulateSEM(adag, N=N)
          lt <- dagitty::localTests(adag, obs.dat, type="cis")
        }
        TRUE
      }, error=function(cond) {
        message(cond)
        return(FALSE)
      })
  }
  return(list(dat=obs.dat, lt=lt))
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


getDAG2DiscrPaths <- function() {
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
  } else if (type == "collfork") {
    return(getDAGCollFork())
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
  } else if (type == "pdsep_g") {
    return(getDAGPdSep())
  } else if (type == "1be") {
    return(getDAG1BE())
  } else if (type == "4be") {
    return(getDAG4BE())
  }
}


