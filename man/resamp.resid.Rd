\name{resamp.resid}
\alias{resamp.resid}
\title{Resample residuals of a mixed-effect model}
\usage{
  resamp.resid(m)
}
\arguments{
  \item{m}{A fitted model of type
  \code{\link[lme4]{mer-class}}}
}
\value{
  A numeric vector of the same length as the original
  response
}
\description{
  Generates a new response vector for a fitted mixed-effect
  model by resampling residuals at all levels (including
  all random effects). All residuals are adjusted to match
  their (co)variance to the variance and covariance
  estimates of the model.
}
\seealso{
  \code{\link{ranef.reflate}}
}

