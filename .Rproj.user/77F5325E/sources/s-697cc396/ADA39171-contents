Ctilde <- function(coords_ptt,
                   z,
                   params = list(sig_sq = 2,
                                 a = 5,
                                 c = 0.5,
                                 kappa = 0.9)) {
  ##################
  # to get started #
  ##################
  # data size
  n <- nrow(coords_ptt)
  H = R <- matrix(0, n, n)

  # select directions
  directions <- unique(as.vector(z))
  K <- length(directions)

  # index set
  nd <- nchar(strsplit(coords_ptt$partition[1], ",")[[1]][1])
  ptts <- sort(unique(coords_ptt$partition))
  idx <- list()
  for (m in ptts) {
    idx[[m]] <- which(coords_ptt$partition == m)
  }

  # assign rownames
  if (is.null(rownames(z))) {
    rownames(z) <- ptts
  }

  # parent partition and idx
  pptts_wind = pptts_list <- list()
  for (m in ptts) {
    # level 1: partition
    tmp <- data.frame(x = m) %>%
      tidyr::separate(x, c("row", "col", "time"),
                      sep = ",", convert = TRUE)

    pptts_wind[1:K] <- list(NULL)
    names(pptts_wind) <- directions

    for (h in directions) {
      # level 2: wind direction
      ppttstmp <- c(parentS(tmp, coords_ptt, h, nd),
                    parentT(tmp, coords_ptt, nd))
      # $ppartition and $idx exist only when there is at least one parents
      if (sum(is.na(ppttstmp)) < 2) {
        pptts_wind_inf <- list()
        pptts_wind_inf[['ppartition']] <- ppttstmp
        pptts_wind_inf[['pidx']] <- c(idx[[ppttstmp[1]]],
                                      idx[[ppttstmp[2]]])                       # change: important to keep ppartition order
        pptts_wind[[h]] <- pptts_wind_inf
      }
    }
    pptts_list[[m]] <- pptts_wind
  }

  a <- params$a
  c <- params$c
  kappa <- params$kappa
  sig_sq <- params$sig_sq

  ##########
  # Ctilde #
  ##########
  for (l in ptts) {
    idx_l <- idx[[l]]
    n_l <- length(idx_l)
    h_l <- z[l,]
    if (length(h_l) > 1) {
      idx_pl <- NULL
      for (s in 1:length(h_l)) {
        idx_pl <- c(idx_pl, pptts_list[[l]][[h_l[s]]]$pidx)
      }
    } else {
      idx_pl <- pptts_list[[l]][[h_l]]$pidx
    }
    coords_tmp <- rbind(coords_ptt[idx_l, c("easting", "northing", "time")],
                        coords_ptt[idx_pl, c("easting", "northing", "time")])
    spatdist <- dist(coords_tmp[,c("easting", "northing")], method = "euclidean")
    timedist <- dist(coords_tmp[,c("time")], method = "manhattan")
    invaup1 <- 1/(a*as.matrix(timedist)+1)
    Cor <- invaup1*exp(-c*as.matrix(spatdist)*(invaup1^(kappa/2)))
    if (is.null(idx_pl)) {
      R[idx_l, idx_l] <- sig_sq*Cor[1:n_l, 1:n_l]
    } else {
      CinvC <- Cor[1:n_l, -c(1:n_l)] %*%
        solve(Cor[-c(1:n_l), -c(1:n_l)] + diag(0.0005, length(idx_pl)))
      R[idx_l, idx_l] <- sig_sq*(Cor[1:n_l, 1:n_l] - CinvC%*%Cor[-c(1:n_l), 1:n_l])
      H[idx_l, idx_pl] <- CinvC
    }
  }

  invImH <- solve(diag(1, n) - H)
  return(invImH %*% R %*% t(invImH))
}
