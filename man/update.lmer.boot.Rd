\name{update.lmer.boot}
\alias{update.lmer.boot}
\title{Update a bootstrapped mixed-effect model with additional replications}
\usage{
  \method{update}{lmer.boot} (object, R, parallel = NULL,
    ...)
}
\arguments{
  \item{object}{An object of class \code{lmer.boot}}

  \item{R}{The number of additional bootstrap replications}

  \item{parallel}{Controls parallel execution, as in
  \code{\link{lmer.boot}}}

  \item{...}{Additional arguments (unused)}
}
\value{
  An object of class \code{lmer.boot}, with both the old
  and new samples
}
\description{
  \code{update} performs additional replications, using the
  same sampling scheme as the original object.
}

