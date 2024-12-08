#' compute single index
#'
#' @param object an aft.nl object
#' @param X a model matrix
#'
#' @export
#'
get.si <- function(sm, coefs , X=NULL, deriv = 0){

  si <- sm$xt$si
  na <- length(si$alpha)
  alpha <- coefs[1:na]
  beta <- coefs[-(1:na)]

  if(is.null(X)){X <- si$X}

  eff <- eff_misi(X, basis = sm$xt$basis)
  eff <- eff$eval(param = c(alpha,beta))

  # unparametrize the single index
  B <- si$B
  a0 <- si$a0
  if( is.null(a0) ){
    a0 <- alpha * 0
  }
  alpha <- drop(B %*% (alpha + a0))
  alpha <- c(exp(alpha[1]), -alpha[-1])/exp(alpha[1])

  out <- list(xba = drop(X %*% alpha) + drop((si$xm) %*% alpha),
              f = eff$f,
              s = log.Stx(eff, deriv = deriv),
              h = log.htx(eff, deriv = deriv),
              alpha = alpha)
}
