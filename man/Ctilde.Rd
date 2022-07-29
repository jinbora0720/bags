\name{Ctilde}
\alias{Ctilde}
\title{The nonstationary covariance matrix induced by G-BAGs}
\usage{
Ctilde(coords_ptt, z, params = list(sig_sq = 2, a = 5, c = 0.5, kappa = 0.9))
}
\arguments{
\item{coords_ptt}{a data frame of \code{n} rows whose columns contain \code{easting},
\code{northing}, \code{time}, \code{row}, \code{col}, \code{time_d}, and \code{partition}.}
\item{z}{character matrix of \code{M} (number of partitions) rows. Each row is the
select direction(s) at each partition.}
\item{params}{list of select parameter values.
\code{sig_sq} is a partial sill.
\code{a} is a temporal decay, \code{c} is a spatial decay, and \code{kappa} is space-time
interaction in a stationary nonseparable spatiotemporal covariance function
\deqn{cor(w(l), w(l + (h,u))) = exp(-c||h||/(a|u|+1)^{\kappa/2})/(a|u|+1)} in
\href{https://doi.org/10.1198/016214502760047113}{Gneiting (2002)}.}
}
\value{
\code{n x n} covariance matrix
}
\description{
Compute a covariance matrix among reference locations
conditional on a select DAG configuration \code{z}
}
\examples{
\dontrun{
  library(tidyverse)
  theme_set(theme_bw())

  # set seed
  seed <- 20220720

  # create coordinates on a grid
  ngrid_e <- 30 # number of grid on the easting axis
  ngrid_n <- 30 # number of grid on the northing axis
  ngrid_t <- 2  # number of grid on the time axis

  grid_e <- seq(0, 1, length = ngrid_e)
  grid_n <- seq(0, 1, length = ngrid_n)
  grid_t <- seq(0, 1, length = ngrid_t)
  coords <- expand.grid(easting = grid_e, northing = grid_n, time = grid_t) \%>\%
    arrange(time, easting, northing)
  n <- nrow(coords)

  # set the number of partitions in each coordinate axis
  n_easting <- 20
  n_northing <- 1
  n_time <- 2

  # set the true direction
  z_true <- rep("W", times = n_easting*n_northing*n_time)

  # set parameters
  a <- 1
  c <- 7
  kappa <- 0
  sig_sq <- 1
  tau_sq <- 0.1
  beta <- -1

  # generate data
  data <- rbag(coords = coords,
               n_partition = c(n_easting, n_northing, n_time),
               breaks_partition = list(breaks_easting = NULL,
                                       breaks_northing = NULL,
                                       breaks_time = NULL),
               z = z_true,
               params = list(tau_sq = tau_sq,
                             sig_sq = sig_sq,
                             a = a, c = c, kappa = kappa,
                             beta = beta),
               seed = seed)

  # extract coordinate-relevant columns
  coords_ptt <- data \%>\%
    select(-y, -w, -starts_with("X"))

  # compute G-BAG derived nonstationary covariance matrix
  nonstCov <- Ctilde(coords_ptt = coords_ptt,
                     z = as.matrix(z_true),
                     params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))

  # plot the covariance with a randomly selected point
  i_ref <- 200
  coords_ptt \%>\%
    mutate(cov = nonstCov[i_ref,]) \%>\%
    ggplot() +
    geom_contour_filled(aes(easting, northing, z = cov)) +
    scale_fill_viridis_d() +
    geom_point(data = coords_ptt[i_ref,],
               aes(easting, northing), size = 2, col = "red") +
    facet_grid( ~ factor(time, label = paste("Time =", sort(unique(time))))) +
    labs(x = "Easting", y = "Northing", fill = "Covariance",
         title = "G-BAG induced covariance: nonstationarity")
}
}
\references{
Jin, B., Peruzzi, M., & Dunson, D. B. (2021).
Bag of DAGs: Flexible Nonstationary Modeling of Spatiotemporal Dependence.
\href{https://doi.org/10.48550/arXiv.2112.11870}{arXiv.2112.11870}
}
\author{
Bora Jin \href{mailto:bora.jin@duke.edu}{bora.jin@duke.edu}
}
