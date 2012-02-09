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

resamp.resid <- function(m) {
  resid.scaled <- attr(VarCorr(m), 'sc') * resid(m) / sd(resid(m))
  resid.boot <- sample(resid.scaled, replace=T)
  ranef.boot <- lapply(ranef.reflate(m), sample.rows, replace=T)
  ranef.boot <- unlist(ranef.boot, use.name=F)
  eta.fix <- as.vector(model.matrix(m) %*% fixef(m))
  eta.ranef <- t(m@Zt) %*% ranef.boot
  as.vector(eta.fix + eta.ranef + resid.boot)
}

lmer.boot <- function(m, R, type=c('parametric','residuals'), ...) {
  type <- match.arg(type)
  boot(
    data=model.response(model.frame(m)),
    statistic=function(data) { 
      m2 <- refit(m, data)
      c(fixef(m2), sapply(VarCorr(m2), diag))
    },
    sim='parametric',
    ran.gen=switch(type,
      parametric=function(data, mle) simulate(mle),
      residuals=function(data, mle) resamp.resid(mle)),
    mle=m,
    R=R,
    ...)
}