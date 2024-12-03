#'
#' Fit AFT models with non-linear cumulative survival function
#'
#' @param formula,family,data same arguments as in [mgcv::gam].
#' @param ... further arguments to be passed to [mgcv::gam].
#' @name aft_nl
#' @rdname aft_nl
#' @export
#'
aft_nl <- function(formula, data = list(),val = list(),outer.args = list(),
                   censored=NULL, min.lambda = 0.005, trace=FALSE,ps = TRUE, ...){

  dat <- .prep.df(formula, data, val,censored, m=outer.args$m, k=outer.args$k, pspline = ps)

  form_comp <- .compile_formula(dat$formula)

  sm <- aftNested:::smooth.construct.misi.smooth.spec(object = form_comp,
                                                data = dat$X, knots = list(),...)
  eff <- eff_misi(sm$xt$si$X, basis = sm$xt$basis)

  lambda <- initial.lambda(sm)
  init.par <- initialize.coef(sm, lambda[1])
  out_env <- new.env()
  lopt <- lambda.optim(eff, sm, delta = data[[censored]], vpen = 1e4, mpen = 1e2,
                       out_env,val = dat$val, delta.val = dat$val[[censored]],
                       init.par = init.par,trace=trace)
  out_l <- optim(par = lambda, fn = lopt, lower = min.lambda, method = "L-BFGS-B")

  out <- out_env$out
  out$lambda <- out_l$par
  out$gcv <- out_l$value
  out$convergence <- out_l$convergence
  out$iterations <- out_l$counts

  class(out) <- c("aftnl", class(out))

  return( out )

}


