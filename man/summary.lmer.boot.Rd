\name{summary.lmer.boot}
\alias{summary.lmer.boot}
\title{Summarize a bootstraped mixed-effect model}
\usage{
  \method{summary}{lmer.boot} (object,
    ci.type = c("perc", "norm", "basic", "stud", "pca"),
    ...)
}
\arguments{
  \item{object}{An object of class \code{lmer.boot}}

  \item{ci.type}{The type of confidence interval, as in
  \code{\link[boot]{boot.ci}}}

  \item{...}{Additional arguments}
}
\value{
  An object of class \code{summary.lmer.boot}
}
\description{
  The \code{summary} method for the results of
  \code{\link{lmer.boot}}, which produces formatted
  summaries of bootstrapped bias, confidence intervals and
  p values for fixed effects and random effects
  covariances.
}
\details{
  Studentized and bias-corrected accelerated (BCA)
  confidence intervals are not supported.  I think that
  studentized will be in the future, but there may be more
  fundamental problems with BCA. P values correspond to the
  largest confidence interval that excludes zero, and
  therefore will change with the selected type of
  confidence interval.
}
\seealso{
  \code{\link{lmer.boot}}
}

