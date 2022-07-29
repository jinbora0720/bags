rbag <- function(coords,
                 n_partition,
                 breaks_partition = list(breaks_easting = NULL,
                                         breaks_northing = NULL,
                                         breaks_time = NULL),
                 z = NULL,
                 params = list(tau_sq = 0.01,
                               sig_sq = 2,
                               a = 5,
                               c = 0.5,
                               kappa = 0.9,
                               beta = NULL),
                 seed = 123) {

  ##################
  # to get started #
  ##################
  # set seed
  set.seed(seed)

  # dimension of data
  n <- nrow(coords)
  d <- ncol(coords)
  if (d == 2) {
    coords <- cbind(coords, 0)
  }
  coords <- data.frame(coords)
  colnames(coords) <- c("easting", "northing", "time")
  rownames(coords) <- NULL
  if (is.null(z)) {
    directions <- "W"
  } else {
    directions <- unique(z)
  }
  K <- length(directions)

  # number of partitions
  if (is.null(breaks_partition$breaks_easting)) {
    n_easting <- n_partition[1]
  } else {
    breaks_easting <- breaks_partition$breaks_easting
    n_easting <- length(breaks_easting) - 1
  }
  if (is.null(breaks_partition$breaks_northing)) {
    n_northing <- n_partition[2]
  } else {
    breaks_northing <- breaks_partition$breaks_northing
    n_northing <- length(breaks_northing) - 1
  }
  if (d == 2) {
    n_time <- 1
  } else {
    if (is.null(breaks_partition$breaks_time)) {
      n_time <- n_partition[3]
    } else {
      breaks_time <- breaks_partition$breaks_time
      n_time <- length(breaks_time) - 1
    }
  }

  a <- params$a
  c <- params$c
  kappa <- params$kappa
  if (!is.null(params$beta)) {
    beta <- params$beta
    p <- length(beta)

    # covariates
    X <- matrix(rnorm(n*p, 0, 0.1), nrow = n, ncol = p)
  }
  sig_sq <- params$sig_sq
  tau_sq <- params$tau_sq

  # assign partitions
  nd <- floor(log10(max(n_northing, n_easting, n_time))) + 1
  coords_ptt <- coords
  if (is.null(breaks_partition$breaks_northing)) {                              # row first col later
    coords_ptt$row <- (n_northing+1) -
      as.numeric(cut_interval(coords_ptt$northing,
                              n = n_northing,
                              labels = 1:n_northing))                           # 1 from the top

  } else {
    coords_ptt$row <- (n_northing+1) -
      as.numeric(cut(coords_ptt$northing,
                     breaks = breaks_northing,
                     labels = 1:n_northing,
                     include.lowest = T))                                       # 1 from the top
  }
  if (is.null(breaks_partition$breaks_easting)) {
    coords_ptt$col <- as.numeric(cut_interval(coords_ptt$easting,
                                              n = n_easting,
                                              labels = 1:n_easting))
  } else {
    coords_ptt$col <- as.numeric(cut(coords_ptt$easting,
                                     breaks = breaks_easting,
                                     labels = 1:n_easting,
                                     include.lowest = T))
  }
  if (is.null(breaks_partition$breaks_time)) {
    coords_ptt$time_d <-
      as.numeric(as.character(coords_ptt$time*(n_time-1) + 1))
  } else {
    coords_ptt$time_d <- as.numeric(cut(coords_ptt$time,
                                        breaks = breaks_time,
                                        labels = 1:n_time,
                                        include.lowest = T))
  }

  coords_ptt$partition <- sprintf(ptt_format(nd),
                                  coords_ptt$row,
                                  coords_ptt$col,
                                  coords_ptt$time_d)

  # index set
  ptts <- sort(unique(coords_ptt$partition))
  idx <- list()
  for (m in ptts) {
    idx[[m]] <- which(coords_ptt$partition == m)
  }

  # parent partition and idx
  pptts_wind = pptts_list <- list()
  for (m in ptts) {
    # level 1: partition
    tmp <- data.frame(x = m) %>%
      tidyr::separate(x,
                      c("row", "col", "time"),
                      sep = ",", convert = TRUE)

    pptts_wind[1:K] <- list(NULL)
    names(pptts_wind) <- directions

    for (h in directions) {
      # level 2: wind directions
      ppttstmp <- c(parentS(tmp, coords_ptt, h, nd),
                    parentT(tmp, coords_ptt, nd))
      # $ppartition and $idx exist only when there is at least one parents
      if (sum(is.na(ppttstmp)) < 2) {
        pptts_wind_inf <- list()
        pptts_wind_inf[['ppartition']] <- ppttstmp
        pptts_wind_inf[['pidx']] <- c(idx[[ppttstmp[1]]],
                                      idx[[ppttstmp[2]]])
        pptts_wind[[h]] <- pptts_wind_inf
      }
    }
    pptts_list[[m]] <- pptts_wind
  }

  # wind directions
  if (is.null(z)) z <- rep(directions[1], length(ptts))
  names(z) <- ptts

  #################
  # generate data #
  #################
  # create w
  w <- rep(0, n)
  for (m in ptts) {
    idx_m <- idx[[m]]
    n_m <- length(idx_m)
    h_m <- z[m]
    idx_p <- pptts_list[[m]][[h_m]]$pidx

    coords_mp <- rbind(coords_ptt[idx_m, c("easting", "northing", "time")],
                       coords_ptt[idx_p, c("easting", "northing", "time")])
    rownames(coords_mp) <- NULL
    spatdist <- dist(coords_mp[,c("easting", "northing")], method = "euclidean")
    timedist <- dist(coords_mp[,c("time")], method = "manhattan")
    invaup1 <- 1/(a*as.matrix(timedist)+1)
    Cor <- invaup1*exp(-c*as.matrix(spatdist)*(invaup1^(kappa/2)))
    if (is.null(idx_p)) {
      LR_m <- t(chol(Cor))
      w[idx_m] <- sqrt(sig_sq)*LR_m %*% rnorm(n_m)
    } else {
      CinvC <- Cor[1:n_m, -c(1:n_m)] %*%
        solve(Cor[-c(1:n_m), -c(1:n_m)] + diag(0.0005, nrow(Cor) - n_m))
      LR_m <- t(chol(Cor[1:n_m, 1:n_m] - CinvC %*% Cor[-c(1:n_m), 1:n_m]))
      H_m <- CinvC
      w[idx_m] <- H_m %*% w[idx_p] + sqrt(sig_sq)*LR_m %*% rnorm(n_m)
    }
  }

  # observations y
  y <- w + sqrt(tau_sq)*rnorm(n)
  if (!is.null(params$beta)) {
    y <- y + X%*%beta
  }

  # save data
  if (is.null(params$beta)) {
    data <- data.frame(coords_ptt, y = y, w = w)
  } else {
    data <- data.frame(coords_ptt, y = y, w = w, X = X)
  }

  return(data)
}
