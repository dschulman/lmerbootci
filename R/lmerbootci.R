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

family.mer <- function(object, ...) {
  family <- object@call$family
  if (is.null(family))
    family <- gaussian()
  if (is.call(family))
    family <- eval(family, envir=parent.frame(2))
  if (is.symbol(family)) 
    family <- as.character(family)
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame(2))
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) 
    stop("'family' not recognized")
  family
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
  f <- family(m)
  ranef.boot <- lapply(ranef.reflate(m), sample.rows, replace=T)
  ranef.boot <- unlist(ranef.boot, use.name=F)
  eta.fix <- as.vector(model.matrix(m) %*% fixef(m))
  eta.ranef <- t(m@Zt) %*% ranef.boot
  if (f$family=='gaussian' && f$link=='identity') {
    resid.scaled <- attr(VarCorr(m), 'sc') * resid(m) / sd(resid(m))
    resid.boot <- sample(resid.scaled, replace=T)
    as.vector(eta.fix + eta.ranef + resid.boot)
  } else if (f$family=='binomial' || f$family=='poisson') {
    mu <- f$linkinv(as.numeric(eta.fix + eta.ranef))
    resid.adjust <- scale(residuals(m), center=T, scale=F)
    resid.boot <- sample(residuals(m), replace=T)
    resid.depearson <- resid.boot*sqrt(f$var(mu)/m@pWt)
    y.boot <- as.vector(mu + resid.depearson)
    if (f$family=='poisson') {
      pmax(0, round(y.boot))
    } else {
      resp <- model.response(m@frame)
      if (! is.matrix(resp)) {
        pmax(0, pmin(1, round(y.boot)))
      } else {
        trials <- rowSums(resp)
        success <- pmax(0, pmin(trials, round(y.boot*trials)))
        cbind(success, trials - success)
      }
    }
  } else
    stop("residual resampling not implemented for this glmm family")
}

setup.cluster <- function(parallel) {
  cl <- NULL
  if (inherits(parallel, 'cluster'))
    cl <- parallel
  else if (!is.null(parallel)) {
    stopifnot(is.numeric(parallel) && parallel>1)
    cl <- makeCluster(parallel, 'PSOCK')
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, library(lmerbootci))
  }
  cl
}

#' Perform a bootstrapped likelihood-ratio test
#' 
#' Given a set of nested models, generated \code{R} bootstrap 
#' replicates from the inner (smaller, null hypothesis) model, and 
#' performs a likelihood-ratio test, using the likelihood ratios of the
#' replicates as a reference distribution.  Both fully parametric
#' and residual resampling are supported.
#' 
#' Note that this function currently does not check if the models
#' are properly nested, and may not give sensible results in other
#' cases.
#' 
#' For models fit with restricted maximum likelihood (REML), this still
#' uses the likelihood (rather than the REML criterion), following the
#' \code{lme4} package. Models are refit to the bootstrap replicates using 
#' whichever method, REML or ML, was used originally.
#' 
#' The \code{parallel} argument can be used two ways: first, it can
#' accept an integer >=2, in which case it will start up a cluster
#' and run it for the duration of the call.  Alternatively, you can
#' give it a \code{cluster} object, as returned by 
#' \code{\link[parallel]{makeCluster}} or related functions.
#'
#' Using a pre-existing cluster is faster if you are doing multiple
#' runs, since it avoids overhead of repeated startups.  However, you
#' must make sure this package is loaded on the cluster: something like
#' \code{clusterEvalQ(cl, library(lmerbootci))} should do it.
#' 
#' @param m0 The inner model (type \code{\link[lme4]{mer-class}})
#' @param m1 The outer model (type \code{\link[lme4]{mer-class}})
#' @param R The number of bootstrap replicates
#' @param type The resampling scheme to use.  Defaults to parametric.
#' @param parallel Specifies whether or how to use parallel execution.
#'    See below for details.  Defaults to no parallelization.
#' @return an object of class \code{lmer.bootlr}, whose \code{print} and
#'    \code{summary} methods can be used to show results.
#' @export
#' @S3method summary lmer.bootlr
#' @S3method print summary.lmer.bootlr
lmer.bootlr <- function(m0, m1, R, type=c('parametric','residuals'), 
                        parallel=NULL) {
  type <- match.arg(type)
  cl <- setup.cluster(parallel)
  starttime <- proc.time()
  lr <- 2 * (logLik(m1, REML=F) - logLik(m0, REML=F))
  newresp <- switch(type, 
                    parametric=simulate, 
                    residuals=resamp.resid)
  f <- function(x) {
    y <- newresp(m0)
    2 * (logLik(refit(m1, y), REML=F) - logLik(refit(m0, y), REML=F))
  }
  if (is.null(cl))
    lrref <- sapply(integer(R), f)
  else
    lrref <- parSapply(cl, integer(R), f)
  elapsed <- proc.time() - starttime
  structure(
    list(lr=lr, lrref=lrref, R=R, type=type, time=elapsed,
         m0=m0, m1=m1),
    class='lmer.bootlr')
}

summary.lmer.bootlr <- function(object, ...) {
  b <- object
  df <- data.frame(
    AIC=c(AIC(b$m0), AIC(b$m1)),
    BIC=c(BIC(b$m0), BIC(b$m1)),
    logLik=c(logLik(b$m0, REML=F), logLik(b$m1, REML=F)),
    chiSq=c(NA, b$lr),
    p=c(NA, mean(b$lrref > b$lr)),
    row.names=c('m0','m1'))
  class(df) <- c('summary.lmer.bootlr', class(df))
  attr(df, 'R') <- b$R
  attr(df, 'type') <- b$type
  attr(df, 'time') <- b$time
  df
}

print.summary.lmer.bootlr <- function(x, digits=max(3, getOption('digits')-3), ...) {
  cat(sprintf('%s Bootstrap likelihood-ratio test (%d replicates)\n', 
              tocaps(attr(x, 'type')),
              attr(x, 'R')))
  printCoefmat(x, digits, has.Pvalue=T, na.print='', ...)
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
    Residual=na.omit(attr(VarCorr(m), 'sc')))
}

#' Bootstrap resampling of a mixed-effect model
#'
#' Generates \code{R} bootstrap replicates of the parameter estimates
#' of a mixed-effect regression model.  Two types of resampling are
#' supported: a fully parametric bootstrap, and a semiparametric
#' bootstrap that resamples residuals and random effects.
#'
#' The \code{parallel} argument can be used two ways: first, it can
#' accept an integer >=2, in which case it will start up a cluster
#' and run it for the duration of the call.  Alternatively, you can
#' give it a \code{cluster} object, as returned by 
#' \code{\link[parallel]{makeCluster}} or related functions.
#'
#' Using a pre-existing cluster is faster if you are doing multiple
#' runs, since it avoids overhead of repeated startups.  However, you
#' must make sure this package is loaded on the cluster: something like
#' \code{clusterEvalQ(cl, library(lmerbootci))} should do it.
#'
#' @param m A fitted model of type \code{\link[lme4]{mer-class}}
#' @param R The number of bootstrap replicates
#' @param type The resampling scheme to use.  Defaults to parametric. 
#' @param parallel Specifies whether or how to use parallel execution.
#'    See below for details.  Defaults to no parallelization.
#' @param ... Additional arguments, passed to \code{\link[boot]{boot}}.
#' @return An object of class \code{lmer.boot}, which is an object of class
#'   \code{\link[boot]{boot}} with some additional information appended.
#' @references Carpenter, J. R., Goldstein, H., and Rasbash, J. 
#'   (2003). A novel bootstrap procedure for assessing the 
#'   relationship between class size and achievement. Journal of the
#'   Royal Statistical Society: Series C (Applied Statistics), 
#'   52(4):431-443.
#' @export
lmer.boot <- function(m, R, type=c('parametric','residuals'), 
                      parallel=NULL, ...) {
  type <- match.arg(type)
  cl <- setup.cluster(parallel)
  elapsed <- system.time(
    b <- boot(
      data=model.response(model.frame(m)),
      statistic=function(data) extract.estimates(refit(m, data)),
      sim='parametric',
      ran.gen=switch(type,
                     parametric=function(data, mle) simulate(mle),
                     residuals=function(data, mle) resamp.resid(mle)),
      mle=m,
      R=R,
      # next two arguments are kind of a hack
      # relying on boot just wanting to know that they are there,
      # then ignoring them and using the value of "cl" instead
      parallel=if (is.null(cl)) 'no' else 'snow',
      ncpus=if (is.null(cl)) 1 else 2,
      cl=cl,
      ...))
  b$type <- type
  b$time <- elapsed
  class(b) <- c('lmer.boot', class(b))
  b
}

boot.2tailp <- function(t, t0, R, type=c('perc','norm','basic')) {
  index.to.p <- function(x) pmax(1, 2*pmin(x, R-x))/R
  switch(
    match.arg(type),
    perc=index.to.p(colSums(t>0)),
    norm=pnorm(0, mean=abs(2*t0-colMeans(t)), sd=apply(t, 2, sd))*2,
    basic=index.to.p(rowSums(apply(t, 1, '<', 2*t0))))
}

coefmat <- function(b, index, ci.type='perc', pvals=T) {
  if (length(pvals)==1)
    pvals <- rep(pvals, length(index))
  t0 <- b$t0[index]
  t <- b$t[,index,drop=F]
  R <- b$R
  cis <- sapply(index, function(i) {
    bci <- boot.ci(b, type=ci.type, index=i)
    switch(ci.type,
           norm=bci$normal[1,2:3],
           perc=bci$percent[1,4:5],
           basic=bci$basic[1,4:5],
           stud=bci$student[1,4:5],
           bca=bci$bca[1,4:5])
  })
  data.frame(
    Estimate=t0,
    Bias=colMeans(t) - t0,
    CI.low=cis[1,], CI.high=cis[2,],
    p.boot=ifelse(pvals, boot.2tailp(t, t0, R, ci.type), NA),
    row.names=names(t0))
}

#' Summarize a bootstraped mixed-effect model
#' 
#' The \code{summary} method for the results of \code{\link{lmer.boot}}, which
#' produces formatted summaries of bootstrapped bias, confidence intervals and
#' p values for fixed effects and random effects covariances.
#'
#' Studentized and bias-corrected accelerated (BCA) confidence 
#' intervals are not supported.  I think that studentized will be in the
#' future, but there may be more fundamental problems with BCA.
#' P values correspond to the largest confidence interval that excludes zero,
#' and therefore will change with the selected type of confidence interval.
#'
#' @param object An object of class \code{lmer.boot}
#' @param ci.type The type of confidence interval, as in 
#'   \code{\link[boot]{boot.ci}}
#' @param ... Additional arguments
#' @return An object of class \code{summary.lmer.boot}
#' @seealso \code{\link{lmer.boot}}
#' @method summary lmer.boot
#' @S3method summary lmer.boot
#' @S3method print summary.lmer.boot
summary.lmer.boot <- function(object, ci.type=c('perc','norm','basic','stud','pca'), ...) {
  ci.type <- match.arg(ci.type)
  if (ci.type=='stud')
    stop('Studentized confidence intervals are not yet supported')
  if (ci.type=='bca')
    stop('Bias-corrected confidence intervals are not supported')
  b <- object
  m <- b$mle
  ind.fixef <- 1:length(fixef(m))
  ind.ranef <- (length(fixef(m))+1):length(b$t0)
  pvals.ranef <- unlist(lapply(VarCorr(m), function(re) {
    n <- nrow(re)
    c(rep(F, n), rep(T, choose(n, 2)))
  }))
  extra <- length(ind.ranef) - length(pvals.ranef)
  pvals.ranef <- c(pvals.ranef, rep(F, extra))
  structure(
    list(R=b$R, type=b$type, ci.type=ci.type, time=b$time,
         fixef=coefmat(b, ind.fixef, ci.type),
         ranef=coefmat(b, ind.ranef, ci.type, pvals.ranef)),
    class='summary.lmer.boot')
}

tocaps <- function(x) {
  paste(toupper(substring(x, 1, 1)), substring(x, 2), sep='')
}

print.summary.lmer.boot <- function(x, digits=max(3, getOption('digits')-3), ...) {
  cat(sprintf('%s Bootstrap (%d replicates)\n', tocaps(x$type), x$R))
  cat(sprintf('%s confidence intervals\n', switch(x$ci.type,
    perc='Percentile',
    basic='Basic',
    norm='Normal')))
  cat('\nFixed Effects:\n')
  printCoefmat(x$fixef, digits, has.Pvalue=T, ...)
  cat('\nRandom Effects and Residuals:\n')
  printCoefmat(x$ranef, digits, has.Pvalue=T, ...)
}

#' Update a bootstrapped mixed-effect model with additional replications
#'
#' \code{update} performs additional replications, using the same sampling
#' scheme as the original object.
#'
#' @param object An object of class \code{lmer.boot}
#' @param R The number of additional bootstrap replications
#' @param parallel Controls parallel execution, as in \code{\link{lmer.boot}}
#' @param ... Additional arguments (unused)
#' @return An object of class \code{lmer.boot}, with both the old and new samples
#' @method update lmer.boot
#' @S3method update lmer.boot
update.lmer.boot <- function(object, R, parallel=NULL, ...) {
  obj2 <- lmer.boot(m=object$mle, R=R, type=object$type, parallel=parallel)
  result <- c(object, obj2)
  result$type <- object$type
  result$time <- object$time + obj2$time
  class(result) <- class(object)
  result
}