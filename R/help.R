#' Bootstrapped confidence intervals for linear mixed-effect regression
#'
#' lmerbootci provides easy methods of performing bootstrap
#' resampling on linear mixed-effect regression models fit with the lme4
#' package.  It implements two types of resampling: a fully parametric
#' bootstrap, and a semiparametric bootstrap which resamples residuals.
#'
#' The core method is \code{\link{lmer.boot}}: given a model object
#' returned by \code{\link[lme4]{lmer}} or related methods, it
#' produces a bootstrapped sample suitable for use with 
#' \code{\link[boot]{boot.ci}} and other methods of the boot package.
#'
#' @docType package
#' @name lmerbootci
#' @aliases lmerbootci lmerbootci-package
#' @importMethodsFrom Matrix t
#' @import parallel
NULL
