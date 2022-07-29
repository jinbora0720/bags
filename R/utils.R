# To mute NOTE by package check "no visible binding for global variable 'x'"
utils::globalVariables(c("x"))
#-------------------------------------------------------------------------------

ptt_format <- function(nd) {
  paste0("%0", nd, ".0f,","%0", nd, ".0f,", "%0", nd, ".0f")
}
#-------------------------------------------------------------------------------

# find parent partition in space
# NA if no parent partition in space
parentS <- function(coords, coords_ptt,
                    wind_dir = c("W", "NW", "S", "SW", "N", "NE", "E", "SE"),
                    nd) {
  # coords should have one row, col, and time
  # coords_ptt should have row, col, and time

  if (wind_dir == "W") {
    prow <- coords$row
    candd <- coords_ptt[coords_ptt$row == prow, c("row","col")]
    if (sum(candd$col < coords$col) == 0) {
      pcol <- NA
    } else {
      pcol <- max(candd$col[candd$col < coords$col])
    }
  }

  if (wind_dir == "E") {
    prow <- coords$row
    candd <- coords_ptt[coords_ptt$row == prow, c("row","col")]
    if (sum(candd$col > coords$col) == 0) {
      pcol <- NA
    } else {
      pcol <- min(candd$col[candd$col > coords$col])
    }
  }

  if (wind_dir == "NW") {
    diff <- coords$row - coords$col
    candd <- coords_ptt[as.numeric(coords_ptt$row) -
                          as.numeric(coords_ptt$col) == diff, c("row", "col")]
    if (sum(candd$col < coords$col)*sum(candd$row < coords$row) == 0) {
      pcol <- NA; prow <- NA
    } else {
      prow <- max(candd$row[candd$row < coords$row])
      pcol <- max(candd$col[candd$col < coords$col])
    }
  }

  if (wind_dir == "SE") {
    diff <- coords$row - coords$col
    candd <- coords_ptt[as.numeric(coords_ptt$row) -
                          as.numeric(coords_ptt$col) == diff, c("row", "col")]
    if (sum(candd$col > coords$col)*sum(candd$row > coords$row) == 0) {
      pcol <- NA; prow <- NA
    } else {
      prow <- min(candd$row[candd$row > coords$row])
      pcol <- min(candd$col[candd$col > coords$col])
    }
  }

  if (wind_dir == "S") {
    pcol <- coords$col
    candd <- coords_ptt[coords_ptt$col == pcol, c("row","col")]
    if (sum(candd$row > coords$row) == 0) {
      prow <- NA
    } else {
      prow <- min(candd$row[candd$row > coords$row])
    }
  }

  if (wind_dir == "N") {
    pcol <- coords$col
    candd <- coords_ptt[coords_ptt$col == pcol, c("row","col")]
    if (sum(candd$row < coords$row) == 0) {
      prow <- NA
    } else {
      prow <- max(candd$row[candd$row < coords$row])
    }
  }

  if (wind_dir == "SW") {
    total <- coords$row + coords$col
    candd <- coords_ptt[as.numeric(coords_ptt$row) +
                          as.numeric(coords_ptt$col) == total, c("row", "col")]
    if (sum(candd$col < coords$col)*sum(candd$row > coords$row) == 0) {
      pcol <- NA; prow <- NA
    } else {
      prow <- min(candd$row[candd$row > coords$row])
      pcol <- max(candd$col[candd$col < coords$col])
    }
  }

  if (wind_dir == "NE") {
    total <- coords$row + coords$col
    candd <- coords_ptt[as.numeric(coords_ptt$row) +
                          as.numeric(coords_ptt$col) == total, c("row", "col")]
    if (sum(candd$col > coords$col)*sum(candd$row < coords$row) == 0) {
      pcol <- NA; prow <- NA
    } else {
      prow <- max(candd$row[candd$row < coords$row])
      pcol <- min(candd$col[candd$col > coords$col])
    }
  }

  return(ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                sprintf(ptt_format(nd), prow, pcol, coords$time)))
}
#-------------------------------------------------------------------------------

# find parent partition in time
# NA if no parent partition in time
parentT <- function(coords, coords_ptt, nd) {
  # coords should have one row, col, and time (always in day)
  # coords_ptt should have row, col, time, and possibly time_d (in day)

  prow <- coords$row
  pcol <- coords$col
  if (is.null(coords_ptt$time_d)) {
    candd <- coords_ptt[which(coords_ptt$row == prow &
                                coords_ptt$col == pcol), "time"]
  } else {
    candd <- coords_ptt[which(coords_ptt$row == prow &
                                coords_ptt$col == pcol), "time_d"]
  }
  if (sum(candd < coords$time) == 0) {
    ptime <- NA
  } else{
    ptime <- max(candd[candd < coords$time])
  }

  return(ifelse(is.na(ptime), NA, sprintf(ptt_format(nd), prow, pcol, ptime)))
}
#-------------------------------------------------------------------------------

# find nearest parent partitions along axes and time
parent4NN <- function(coords, coords_ptt, nd) {
  # coords should have one row, col, and time
  # coords_ptt should have row, col, and time_d

  # spatial
  # left and right
  prow <- coords$row
  candd <- coords_ptt[(coords_ptt$row == prow & coords_ptt$time_d == coords$time),
                      c("row","col")]
  if (sum(candd$col < coords$col) == 0) {
    pcol <- NA
  } else {
    pcol <- max(candd$col[candd$col < coords$col])
  }
  parentl <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                    sprintf(ptt_format(nd), prow, pcol, coords$time))

  if (sum(candd$col > coords$col) == 0) {
    pcol <- NA
  } else {
    pcol <- min(candd$col[candd$col > coords$col])
  }
  parentr <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                    sprintf(ptt_format(nd), prow, pcol, coords$time))

  # up and down
  pcol <- coords$col
  candd <- coords_ptt[(coords_ptt$col == pcol & coords_ptt$time_d == coords$time),
                      c("row","col")]
  if (sum(candd$row > coords$row) == 0) {
    prow <- NA
  } else {
    prow <- min(candd$row[candd$row > coords$row])
  }
  parentd <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                    sprintf(ptt_format(nd), prow, pcol, coords$time))

  if (sum(candd$row < coords$row) == 0) {
    prow <- NA
  } else {
    prow <- max(candd$row[candd$row < coords$row])
  }
  parentu <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                    sprintf(ptt_format(nd), prow, pcol, coords$time))

  # temporal future and past
  candd <- coords_ptt[(coords_ptt$col == coords$col & coords_ptt$row == coords$row),
                      "time_d"]
  if (sum(candd < coords$time) == 0) {
    ptime <- NA
  } else {
    ptime <- max(candd[candd < coords$time])
  }
  parentp <- ifelse(is.na(ptime), NA,
                    sprintf(ptt_format(nd), coords$row, coords$col, ptime))

  if (sum(candd > coords$time) == 0) {
    ptime <- NA
  } else {
    ptime <- min(candd[candd > coords$time])
  }
  parentf <- ifelse(is.na(ptime), NA,
                    sprintf(ptt_format(nd), coords$row, coords$col, ptime))

  parents <- c(parentl, parentr, parentd, parentu, parentp, parentf)
  return(parents[!is.na(parents)])
}
#-------------------------------------------------------------------------------

# find nearest parent partitions along axes and time
parent8NN <- function(coords, coords_ptt, nd) {
  # coords should have one row, col, and time
  # coords_ptt should have row, col, and time_d

  # spatial
  # left and right
  prow <- coords$row
  candd <- coords_ptt[(coords_ptt$row == prow & coords_ptt$time_d == coords$time),
                      c("row","col")]
  if (sum(candd$col < coords$col) == 0) {
    pcol <- NA
  } else {
    pcol <- max(candd$col[candd$col < coords$col])
  }
  parentl <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                    sprintf(ptt_format(nd), prow, pcol, coords$time))

  if (sum(candd$col > coords$col) == 0) {
    pcol <- NA
  } else {
    pcol <- min(candd$col[candd$col > coords$col])
  }
  parentr <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                    sprintf(ptt_format(nd), prow, pcol, coords$time))

  # up and down
  pcol <- coords$col
  candd <- coords_ptt[(coords_ptt$col == pcol & coords_ptt$time_d == coords$time),
                      c("row","col")]
  if (sum(candd$row > coords$row) == 0) {
    prow <- NA
  } else {
    prow <- min(candd$row[candd$row > coords$row])
  }
  parentd <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                    sprintf(ptt_format(nd), prow, pcol, coords$time))

  if (sum(candd$row < coords$row) == 0) {
    prow <- NA
  } else {
    prow <- max(candd$row[candd$row < coords$row])
  }
  parentu <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                    sprintf(ptt_format(nd), prow, pcol, coords$time))

  # diagonal
  diff <- coords$row - coords$col
  candd <- coords_ptt[(as.numeric(coords_ptt$row) - as.numeric(coords_ptt$col)) == diff &
                        coords_ptt$time_d == coords$time, c("row", "col")]
  if (sum(candd$col < coords$col)*sum(candd$row < coords$row) == 0) {
    pcol <- NA; prow <- NA
  } else {
    prow <- max(candd$row[candd$row < coords$row])
    pcol <- max(candd$col[candd$col < coords$col])
  }
  parentnw <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                     sprintf(ptt_format(nd), prow, pcol, coords$time))

  if (sum(candd$col > coords$col)*sum(candd$row > coords$row) == 0) {
    pcol <- NA; prow <- NA
  } else {
    prow <- min(candd$row[candd$row > coords$row])
    pcol <- min(candd$col[candd$col > coords$col])
  }
  parentse <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                     sprintf(ptt_format(nd), prow, pcol, coords$time))

  total <- coords$row + coords$col
  candd <- coords_ptt[(as.numeric(coords_ptt$row) + as.numeric(coords_ptt$col)) == total &
                        coords_ptt$time_d == coords$time, c("row", "col")]
  if (sum(candd$col < coords$col)*sum(candd$row > coords$row) == 0) {
    pcol <- NA; prow <- NA
  } else {
    prow <- min(candd$row[candd$row > coords$row])
    pcol <- max(candd$col[candd$col < coords$col])
  }
  parentsw <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                     sprintf(ptt_format(nd), prow, pcol, coords$time))

  if (sum(candd$col > coords$col)*sum(candd$row < coords$row) == 0) {
    pcol <- NA; prow <- NA
  } else {
    prow <- max(candd$row[candd$row < coords$row])
    pcol <- min(candd$col[candd$col > coords$col])
  }
  parentne <- ifelse(is.na(prow) + is.na(pcol) > 0, NA,
                     sprintf(ptt_format(nd), prow, pcol, coords$time))

  # temporal future and past
  candd <- coords_ptt[(coords_ptt$col == coords$col & coords_ptt$row == coords$row),
                      "time_d"]
  if (sum(candd < coords$time) == 0) {
    ptime <- NA
  } else {
    ptime <- max(candd[candd < coords$time])
  }
  parentp <- ifelse(is.na(ptime), NA,
                    sprintf(ptt_format(nd), coords$row, coords$col, ptime))

  if (sum(candd > coords$time) == 0) {
    ptime <- NA
  } else {
    ptime <- min(candd[candd > coords$time])
  }
  parentf <- ifelse(is.na(ptime), NA,
                    sprintf(ptt_format(nd), coords$row, coords$col, ptime))

  parents <- c(parentl, parentr, parentd, parentu,
               parentnw, parentse, parentsw, parentne,
               parentp, parentf)
  return(parents[!is.na(parents)])
}
#-------------------------------------------------------------------------------

# transform psi to theta
logloglogit <- function(psi, unifprior = FALSE,
                        la = NULL, ua = NULL, lc = NULL, uc = NULL) {
  if (unifprior) {
    return(c(log((psi[1]-la)/(ua-psi[1])),
             log((psi[2]-lc)/(uc-psi[2])),
             log(psi[3]/(1-psi[3]))))
  } else {
    return(c(log(psi[1]), log(psi[2]), log(psi[3]/(1-psi[3]))))
  }
}
#-------------------------------------------------------------------------------

# beta update
beta_update <- function(y, X, w, tau_sq, mu0, invV0) {
  V <- solve(invV0 + crossprod(X)/tau_sq)
  LV <- t(chol(V))
  return(V%*%(invV0%*%mu0 + crossprod(X, y-w)/tau_sq) + LV%*%rnorm(length(mu0)))
}
#-------------------------------------------------------------------------------

# tau_sq update
tau_sq_update <- function(res, w, at, bt) {
  n <- length(w)
  eps <- res - w
  1/rgamma(1, at + n/2, rate = bt + sum(eps^2)/2)
}

