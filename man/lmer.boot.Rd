\name{lmer.boot}
\alias{lmer.boot}
\title{Bootstrap resampling of a mixed-effect model}
\usage{
  lmer.boot(m, R, type = c("parametric", "residuals"), ...)
}
\arguments{
  \item{m}{A fitted model of type
  \code{\link[lme4]{mer-class}}}

  \item{R}{The number of bootstrap replicates}

  \item{type}{The resampling scheme to use.  Defaults to
  parametric.}

  \item{...}{Additional arguments, passed to
  \code{\link[boot]{boot}}.}
}
\value{
  An object of class \code{lmer.boot}, which is an object
  of class \code{\link[boot]{boot}} with some additional
  information appended.
}
\description{
  Generates \code{R} bootstrap replicates of the parameter
  estimates of a mixed-effect regression model.  Two types
  of resampling are supported: a fully parametric
  bootstrap, and a semiparametric bootstrap that resamples
  residuals and random effects.
}
\references{
  Carpenter, J. R., Goldstein, H., and Rasbash, J.  (2003).
  A novel bootstrap procedure for assessing the
  relationship between class size and achievement. Journal
  of the Royal Statistical Society: Series C (Applied
  Statistics), 52(4):431-443.
}

