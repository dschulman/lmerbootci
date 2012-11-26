\name{lmer.bootlr}
\alias{lmer.bootlr}
\title{Perform a bootstrapped likelihood-ratio test}
\usage{
  lmer.bootlr(m0, m1, R,
    type = c("parametric", "residuals"), parallel = NULL)
}
\arguments{
  \item{m0}{The inner model (type
  \code{\link[lme4]{mer-class}})}

  \item{m1}{The outer model (type
  \code{\link[lme4]{mer-class}})}

  \item{R}{The number of bootstrap replicates}

  \item{type}{The resampling scheme to use.  Defaults to
  parametric.}

  \item{parallel}{Specifies whether or how to use parallel
  execution.  See below for details.  Defaults to no
  parallelization.}
}
\value{
  an object of class \code{lmer.bootlr}, whose \code{print}
  and \code{summary} methods can be used to show results.
}
\description{
  Given a set of nested models, generated \code{R}
  bootstrap replicates from the inner (smaller, null
  hypothesis) model, and performs a likelihood-ratio test,
  using the likelihood ratios of the replicates as a
  reference distribution.  Both fully parametric and
  residual resampling are supported.
}
\details{
  Note that this function currently does not check if the
  models are properly nested, and may not give sensible
  results in other cases.

  For models fit with restricted maximum likelihood (REML),
  this still uses the likelihood (rather than the REML
  criterion), following the \code{lme4} package. Models are
  refit to the bootstrap replicates using whichever method,
  REML or ML, was used originally.

  The \code{parallel} argument can be used two ways: first,
  it can accept an integer >=2, in which case it will start
  up a cluster and run it for the duration of the call.
  Alternatively, you can give it a \code{cluster} object,
  as returned by \code{\link[parallel]{makeCluster}} or
  related functions.

  Using a pre-existing cluster is faster if you are doing
  multiple runs, since it avoids overhead of repeated
  startups.  However, you must make sure this package is
  loaded on the cluster: something like
  \code{clusterEvalQ(cl, library(lmerbootci))} should do
  it.
}

