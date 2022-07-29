\name{bags-package}
\alias{bags-package}
\alias{bags}
\docType{package}
\title{
\packageTitle{bags}
}
\description{
Fits Bayesian spatial or spatiotemporal univariate regression models with latent
Gaussian Bag of Directed Acyclic Graphs Processes (G-BAGs)
}
\details{
The function \code{\link{bag}} conducts posterior sampling using
Markov chain Monte Carlo for univariate Bayesian regression models with
Gaussian errors whose spatial or spatiotemporal random effects are modeled with
G-BAGs.

The function \code{\link{rbag}} generates random data from prior sampling of
univariate Bayesian regression models based on G-BAGs with Gaussian errors.
}
\references{
Jin, B., Peruzzi, M., & Dunson, D. B. (2021).
Bag of DAGs: Flexible Nonstationary Modeling of Spatiotemporal Dependence.
\href{https://doi.org/10.48550/arXiv.2112.11870}{arXiv.2112.11870}
}
\author{
Bora Jin \href{mailto:bora.jin@duke.edu}{bora.jin@duke.edu}
}
\keyword{ package }
\seealso{
\code{\link{bag}}, \code{\link{rbag}}
}
