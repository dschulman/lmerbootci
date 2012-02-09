\name{ranef.reflate}
\alias{ranef.reflate}
\title{Shrinkage-corrected estimates of random effects}
\usage{
  ranef.reflate(m)
}
\arguments{
  \item{m}{A fitted model of type
  \code{\link[lme4]{mer-class}}}
}
\value{
  A list of data frames, one for each grouping factor for
  the random effects.
}
\description{
  Give estimates of the random effects from a fitted model
  produced by \code{\link[lme4]{lmer}} or related
  functions.
}
\details{
  Relative to the shrinkage estimates given by
  \code{\link[lme4]{ranef}}, these have been transformed so
  that the (co)variance of the estimates matches the
  estimated (co)variance of the random-effect terms. They
  should be more suitable for use in bootstrapping or other
  resampling.
}
\references{
  Carpenter, J. R., Goldstein, H., and Rasbash, J.  (2003).
  A novel bootstrap procedure for assessing the
  relationship between class size and achievement. Journal
  of the Royal Statistical Society: Series C (Applied
  Statistics), 52(4):431-443.
}

