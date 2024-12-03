#'
#' Log-Hazard function
#'
#' @rdname log.htx
#' @export log.htx
#'
DHessDrho.htx <- function(o, delta = NULL, llk, DbDr){

  param <- o$param

  if(o$deriv < 4){o <- o$eval(param = param, deriv = 4)}

  if(is.null(delta)){
    delta <- rep(TRUE, length(o$f))
  } else{delta <- as.logical(delta)}

  # 1-4th order derivatives of the smooth
  d1 <- o$store$f1[delta]; d2 <- o$store$f2[delta]
  d3 <- o$store$f3[delta]; d4 <- o$store$f4[delta]

  # model matrices
  Xi <- o$store$Xi[delta,]; X <- o$store$X0[delta,]
  X1 <- o$store$X1[delta,]; X2 <- o$store$X2[delta,]

  n <- nrow(Xi)
  na <- o$na
  alpha <- param[1:na]; beta <- param[-(1:na)]
  betat <- c(beta[1], exp(beta[-1]))
  btb <- c(1,betat[-1]) #d betatilde / d beta


  # Derivatives of alpha and beta wrt log-smoothing parameters
  dArho = DbDr[1:na, , drop = FALSE]
  dBrho = DbDr[-(1:na), , drop = FALSE]

  f3 <- d3 + d4/d1 - d2/d1^2*(3*d3 - 2*d2^2/d1)
  g3 <- 2*n/alpha[1]^3
  df2dbt <- X1%*%(-d3/d1^2 + 2*(d2/d1)^2) + X2%*%(1-2*d2/d1^2) + X3/d1
  df1dbt <- X1 + (X2*d1 - X1*d2)/(d1^2)
  df1dbtbt <- -2/d1^3 * (t(X2) %*% X1) - t(X1/d1^3) %*% (X2*d1 - 2*d1*d2*X1)

  der3 <- list()

  coun <- 1
  for(jj in 1:na){      # AXX
    XJ <- Xi[ , jj]
    for(kk in jj:na){   # AAX
      XJK <- XJ * Xi[ , kk]
      for(ll in kk:na){ # AAA
        der3[[coun]] <-  sum(f3 * XJK * Xi[ , ll])
        if(all.equal(rep(1,3), c(jj,kk,ll))){der3[[coun]] <- der3[[coun]] + g3}
        coun <- coun + 1
      }
      for(ll in 1:nb){  # AAB
        der3[[coun]] <- sum(df2dbt[,ll]*btb[ll] * XJK)
        coun <- coun + 1
      }
    }
    for(kk in 1:nb){    # ABB
      XJK <- XJ * X[ , kk]
      for(ll in kk:nb){
        df1dbtbt <- -2/d1^3 * (X2[,ll] * X1[,kk]) - t(X1[,ll]/d1^3) * (X2[,kk]*d1 - 2*d1*d2*X1[,kk])
        der3[[coun]] <- sum(df1dbtbt*btb[kk]*btb[ll]*XJ)
        if(kk==ll){
          der3[[coun]] <- der3[[coun]] + sum(df1dbt[,kk] * XJ)
        }
        coun <- coun + 1
      }
    }
  }

  for(jj in 1:nb){      # BBB
    XJ <- leee * X[ , jj]
    for(kk in jj:nb){
      XJK <- XJ * X[ , kk]
      for(ll in kk:nb){
        der3[[coun]] <- sum( XJK * X[ , ll] ) # Same as t(leee) %*% (X[ , jj] * X[ , kk] * X[ , ll])
        coun <- coun + 1
      }
    }
  }
  der3 <- do.call("c", der3)

  return( out )
}
