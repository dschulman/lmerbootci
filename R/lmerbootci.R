#' Shrinkage-corrected estimates of random effects
#'
#' Give estimates of the random effects from a fitted
#' model produced by \code{\link[lme4]{lmer}} or related 
#' functions.
#'
#' Relative to the shrinkage estimates given by 
#' \code{\link[lme4]{ranef}}, these have been transformed 
#' so that the (co)variance of the estimates matches 
#' the estimated (co)variance of the random-effect terms.
#' They should be more suitable for use in bootstrapping or
#' other resampling.
#' 
#' @param m A fitted model of type \code{\link[lme4]{mer-class}}
#' @return A list of data frames, one for each grouping factor
#'   for the random effects.
#' @export
#' @references Carpenter, J. R., Goldstein, H., and Rasbash, J. 
#'   (2003). A novel bootstrap procedure for assessing the 
#'   relationship between class size and achievement. Journal of the
#'   Royal Statistical Society: Series C (Applied Statistics), 
#'   52(4):431-443.
ranef.reflate <- function(m) {
  mapply(function(vc.est, resid) {
    u <- as.matrix(resid)
    vc.emp <- cov(u)
    a <- solve(chol(vc.emp)) %*% chol(vc.est)
    as.data.frame(u %*% a)
  }, vc.est=VarCorr(m), resid=ranef(m), SIMPLIFY=F)
}

sample.rows <- function(df, size=nrow(df), replace=F, prob=NULL) {
  df[sample.int(nrow(df), size, replace, prob),]
}

#' Resample residuals of a mixed-effect model
#'
#' Generates a new response vector for a fitted mixed-effect model
#' by resampling residuals at all levels (including all random effects).
#' All residuals are adjusted to match their (co)variance to the
#' variance and covariance estimates of the model.
#'
#' @param m A fitted model of type \code{\link[lme4]{mer-class}}
#' @return A numeric vector of the same length as the original response
#' @seealso \code{\link{ranef.reflate}}
#' @export
resamp.resid <- function(m) {
  resid.scaled <- attr(VarCorr(m), 'sc') * resid(m) / sd(resid(m))
  resid.boot <- sample(resid.scaled, replace=T)
  ranef.boot <- lapply(ranef.reflate(m), sample.rows, replace=T)
  ranef.boot <- unlist(ranef.boot, use.name=F)
  eta.fix <- as.vector(model.matrix(m) %*% fixef(m))
  eta.ranef <- t(m@Zt) %*% ranef.boot
  as.vector(eta.fix + eta.ranef + resid.boot)
}

as.named.vector <- function(x, sep=':') {
  structure(
    as.vector(x),
    names=apply(expand.grid(dimnames(x)), 1, paste, collapse=sep))
}

extract.estimates <- function(m) {
  c(fixef(m),
    unlist(lapply(VarCorr(m), function (group) {
      sds <- attr(group, 'stddev')
      cors <- attr(group, 'correlation')
      c(sds, as.named.vector(cors)[lower.tri(cors)])
    })),
    Residual=attr(VarCorr(m), 'sc'))
}

#' Bootstrap resampling of a mixed-effect model
#'
#' Generates \code{R} bootstrap replicates of the parameter estimates
#' of a mixed-effect regression model.  Two types of resampling are
#' supported: a fully parametric bootstrap, and a semiparametric
#' bootstrap that resamples residuals and random effects.
#'
#' @param m A fitted model of type \code{\link[lme4]{mer-class}}
#' @param R The number of bootstrap replicates
#' @param type The resampling scheme to use.  Defaults to parametric. 
#' @param ... Additional arguments, passed to \code{\link[boot]{boot}}.
#' @return An object of class \code{\link[boot]{boot}}
#' @references Carpenter, J. R., Goldstein, H., and Rasbash, J. 
#'   (2003). A novel bootstrap procedure for assessing the 
#'   relationship between class size and achievement. Journal of the
#'   Royal Statistical Society: Series C (Applied Statistics), 
#'   52(4):431-443.
#' @export
lmer.boot <- function(m, R, type=c('parametric','residuals'), ...) {
  type <- match.arg(type)
  boot(
    data=model.response(model.frame(m)),
    statistic=function(data) extract.estimates(refit(m, data)),
    sim='parametric',
    ran.gen=switch(type,
      parametric=function(data, mle) simulate(mle),
      residuals=function(data, mle) resamp.resid(mle)),
    mle=m,
    R=R,
    ...)
}
