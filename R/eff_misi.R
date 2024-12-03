#'
#' Build single index effect
#'
#' @param Xi matrix to be projected via single index vector \code{alpha}.
#' @param basis function which takes \code{si = Xi\%*\%alpha} as input and returns model
#'                  matrix and its derivatives w.r.t. \code{si}.
#' @name eff_misi
#' @rdname eff_misi
#' @export eff_misi
#'
eff_misi <- function(Xi, basis, a0 = NULL){

  force(Xi); force(basis); force(a0)
  unif <- basis$unif

  eval <- function(param, deriv = 0, kde=NULL, kernel="gaussian"){
    na <- ncol( Xi )
    nb <- length(param) - na

    # Single index and spline coefficients
    alpha <- param[1:na]
    alpha[1] <- exp(alpha[1])
    alpha[-1] <- - alpha[ -1 ]
    beta <- param[ -(1:na) ]
    if(basis$drop.coef) {betat <- exp(beta)
    } else {betat <- c(beta[1], exp(beta[-1]))}

    if( is.null(a0) ){
      a0 <- alpha * 0
    }

    # Project covariates on single index vector
    ax1 <- drop( Xi %*% (alpha + a0) )
    if(unif & is.null(kde)){
      # kde <- .ckde(ax1, h = (4/(3*nrow(Xi)))^(1/5), Xi_fix = Xi, kernel = "gauss")
      if(kernel=="gaussian"){
        ss <- seq(1, nrow(Xi), length.out=min(nrow(Xi), 100))
        kde <- .ckde(ax1[ss], h = (4/(3*nrow(Xi[ss,])))^(1/5), Xi_fix = Xi[ss,], kernel = "gauss")
      }
      if(kernel=="triweight"){
        kde <- .ckde(ax1, h = (4/(3*nrow(Xi)))^(1/5), Xi_fix = Xi, kernel = "tri")
      }
      ax <- kde(ax1,Xi = Xi, deriv = deriv)
    } else if (unif & !is.null(kde)){
      ax <- kde(ax1,Xi = Xi, deriv = deriv)
    } else {ax <- list(d0 = ax1,
                       d1 = Xi)}

    # Build P-spline basis and its derivatives
    store <- basis$evalX(x = ax$d0, deriv = deriv)
    store$Xi <- Xi
    store$g1 <- c(alpha[1], rep(-1, na-1))
    if( deriv >= 1 ){

      if(basis$drop.coef) {btb <- betat
      } else {btb <- c(1,betat[-1])}

      store$f1 <- drop( store$X1 %*% betat )
      store$u1 <- ax$d1
      store$fb <- sweep(store$X0, 2, btb, FUN="*")
      store$f1b <- sweep(store$X1, 2, btb, FUN="*")

      if( deriv >= 2 ){
        store$f2 <- drop( store$X2 %*% betat )
        store$fbb <- store$fb[,-1]
        store$u2 <- ax$d2
        # store$g2 <- diag(c(alpha[1], rep(0,na)))
        if( deriv >= 3 ){
          store$f3 <- drop( store$X3 %*% betat )
          store$g3 <- store$g2
        }# end deriv 3
      }# end deriv 2
    }# end deriv 1

    o <- eff_misi(Xi = Xi, basis = basis)
    o$f <- drop( store$X0 %*% betat )
    o$ax <- ax$d0
    o$param <- param
    o$a0 <- a0
    o$na <- na
    o$store <- store
    o$deriv <- deriv
    o$kde <- kde
    o$drop.coef <- basis$drop.coef

    return( o )

  }

  out <- structure(list("eval" = eval,
                        unif = unif), class = c("misi", "nested"))

  return( out )

}









