#'
#' Log-survival function
#'
#' @rdname log.Stx
#' @export log.Stx
#'
log.Stx <- function(o, deriv = 0, param = NULL){

  if( is.null(param)){
    param <- o$param
    if( is.null(param) ){ stop("param vector not provided!") }
  } else if (!is.null(param)) {
    o <- o$eval(param = param, deriv = deriv)
  }
  if (o$deriv<=deriv) {
    o <- o$eval(param = param, deriv = deriv)
  }

  d0 <- o$f; d1 <- o$store$f1; d2 <- o$store$f2
  Xi <- o$store$Xi; X <- o$store$X0; X1 <- o$store$X1
  u1 <- o$store$u1; u2 <- o$store$u2;
  g1 <- o$store$g1

  alpha <- param[1:o$na]
  alpha <- c(exp(alpha[1]), -alpha[-1]) # reparametrise alpha
  beta <- param[-(1:o$na)]
  if(o$drop.coef) {betat <- exp(beta)
  } else {betat <- c(beta[1], exp(beta[-1]))}

  n <- nrow(Xi)

  out <- list()
  out$s <- -exp(d0)

  if(deriv >=1){
    if(o$drop.coef) {btb <- betat
    } else {btb <- c(1,betat[-1])}

    aQ <- sweep(u1, 2, g1, "*")
    bM <- sweep(X, 2,btb , "*")
    bM1 <- sweep(X1, 2,btb , "*")

    out$dhda <- crossprod(out$s*d1,aQ)
    out$dhdb <- crossprod(out$s, bM)

    if(deriv >=2){

      if(o$unif){
        daa <- diag(rep(0, o$na))
        daa[lower.tri(daa, diag = TRUE)] <- crossprod(out$s*d1, u2)
        daa[upper.tri(daa, diag=FALSE)] <- t(daa)[upper.tri(daa, diag = FALSE)]
        daa <- daa*g1[1,]%o%g1[1,]
      } else {daa <- 0}

      h_aa <- crossprod((out$s * (d1^2 +d2)) *aQ, aQ) - daa
      h_aa[1,1] <- h_aa[1,1] + out$dhda[1]

      h_bb <- crossprod(out$s * bM,bM)
      if(o$drop.coef) {
        diag(h_bb) <- diag(h_bb) + (out$dhdb)
      } else {diag(h_bb)[-1] <- diag(h_bb)[-1] + (out$dhdb)[-1]}

      h_ba <- crossprod(out$s * (d1* bM + bM1), aQ)

      out$dhab <- rbind(cbind(h_aa, t(h_ba)),
                        cbind(h_ba, h_bb))
    }
  }
  return( out )
}
