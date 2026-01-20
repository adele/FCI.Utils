rm(list=ls())

size = "x_large"
seed = 2
type = "dense"

sizes = c("small", "med", "large") #, "x_large")
seeds = c(1,2,3)


data_folder <- "../../SCARY_CLeaR2023/SCARY-main/CSV_unfaithful/sparse/"
if (type == "dense") {
  data_folder <- "../../SCARY_CLeaR2023/SCARY-main/CSV_unfaithful/dense/"
}

output_folder <- "../../SCARY_results/sparse/"
if (type == "dense") {
  output_folder <- "../../SCARY_results/dense/"
}
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = T)
}


for (size in sizes) {
  for (seed in seeds) {
    cat("Processing size:", size, ", seed:", seed, "\n")

    data_id <- paste0(size, "_mixed_unfaithful_", seed)
    if (type == "dense") {
      data_id <- paste0("dense_", size, "_mixed_unfaithful_", seed)
    }

    cur_output_folder <- paste0(output_folder, data_id, "/")
    if (!dir.exists(cur_output_folder)) {
      dir.create(cur_output_folder, recursive = T)
    }

    cur_lat_dirs <- list.dirs(cur_output_folder, recursive = F)
    done = (length(cur_lat_dirs) > 0)

    while(!done) {
      data <- read.csv(paste0(data_folder, data_id, "_data.csv"), header = F)
      #str(data)
      full_data_file <- paste0(cur_output_folder, "full_data.csv")
      if (!file.exists(full_data_file)) {
        write.csv(data, file=full_data_file)
      }

      amat <- read.csv(paste0(data_folder,  data_id, "_graph.csv"), header = F)
      row.names(amat) <- colnames(amat)
      amat <- t(amat)
      renderDAG(amat, add_index = T) # shows all variables

      p <- ncol(amat)
      nlat <- 0.2*p
      # lat_cols <- c(33,36,37,41,43,7,20,42,46,47) # sparse, x_large, 1

      # lat_cols <- c(3, 9) # dense, small, 1
      # lat_cols <- c(7, 9) # dense, small, 2
      # lat_cols <- c(4, 5) # dense, small, 3

      # lat_cols <- c(2,5,13) # dense, med, 3

      # lat_cols <- c(11, 4, 19, 22, 3) # dense, large, 1
      # lat_cols <- c(2, 6, 18, 9, 16) # dense, large, 2
      # lat_cols <- c(1, 5, 17,22, 8) # dense, large, 3

      # lat_cols <- sample(size = nlat, x = 1:p, replace = FALSE)

      lat <- c(colnames(amat)[sort(lat_cols)])

      dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
      dagitty::latents(dagg) <- lat
      dagitty::impliedConditionalIndependencies(dagg)
      # Check the true PAG
      true.amat.pag <- getTruePAG(dagg, verbose = TRUE, m.max = 2)@amat
      renderAG(true.amat.pag)

      ntails <- length(which(true.amat.pag == 3))
      cat("Number of tails:", ntails, "\n")
      nbidir <- length(intersect(which(t(true.amat.pag) == 2),
                                 which(true.amat.pag == 2)))
      cat("Number of bidirected edges:", nbidir, "\n")

      if (ntails  > 1) {
        # Check the true min sepsets
        true.sepset <- getPAGImpliedSepset(true.amat.pag)
        formatSepset(true.sepset)
        max.ord <- max(unlist(lapply(true.sepset, function(x) {
          lapply(x, function(y) { ( lapply(y, function(w) {
            length(w)} )   ) } )  })))
        cat("Max.ord: ", max.ord, "\n")
        if (max.ord %in% c(1,2)) {
          cur_lat_id <- paste0("lat_", paste0(lat, collapse = "."))
          cur_lat_folder <- paste0(cur_output_folder, cur_lat_id, "/")
          if (!dir.exists(cur_lat_folder)) {
            dir.create(cur_lat_folder, recursive = T)
          }

          true_objs <- list()
          true_objs$amat <- amat
          true_objs$lat <- lat
          true_objs$dagg <- dagg
          true_objs$true.amat.pag <- true.amat.pag
          true_objs$true.sepset <- true.sepset

          obs_data <- data[,-which(colnames(data) %in% lat)]
          write.csv(obs_data, file=paste0(cur_lat_folder, "obs_data.csv"), row.names = F)
          save(true_objs, file=paste0(cur_lat_folder, "true_objs.RData"))

          renderAG(true.amat.pag, add_index = T, output_folder = cur_lat_folder,
                   fileid = paste0(data_id, "_", cur_lat_id))
          done = TRUE
        }
      }
    }
  }
}

