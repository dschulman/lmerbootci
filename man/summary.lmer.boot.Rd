\name{summary.lmer.boot}
\alias{summary.lmer.boot}
\title{Summarize a bootstraped mixed-effect model}
\usage{
  \method{summary}{lmer.boot} (object, ...)
}
\arguments{
  \item{object}{An object of class \code{lmer.boot}}

  \item{...}{Additional arguments}
}
\value{
  An object of class \code{summary.lmer.boot}
}
\description{
  The \code{summary} method for the results of
  \code{\link{lmer.boot}}, which produces formatted
  summaries of bootstrapped standard errors, bias, and
  percentile-based p values for fixed and random effects.
}
\seealso{
  \code{\link{lmer.boot}}
}

