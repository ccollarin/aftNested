######
# Penalty on variance of a single index vector
#
pen_var_si <- function(o, v, deriv = 0){

  na <- o$na
  a <- o$param[1:na]
  a[1] <- exp(a[1])
  a[-1] <- - a[ -1 ]
  a <- a + o$a0

  x <- o$store$Xi

  p <- length(a)
  ax <- x%*%a

  # Empirical variance
  vhat <- sum(ax^2)/length(ax) - mean(ax)^2

  # Loss
  l0 <- (vhat - v)^2

  l1 <- l2 <- l3 <- NULL
  if( deriv ){
    n <- nrow(x)

    S <- cov(x) * (n-1) / n
    Sa <- (S %*% a)

    l1 <- 4 * (vhat - v) * Sa * o$store$g1

    if( deriv > 1 ){
      l2 <- 8 * tcrossprod(Sa, Sa) + 4 * (vhat - v) * S
      l2 <- l2*tcrossprod(o$store$g1, o$store$g1)
      l2[1,1] <- l2[1,1] + l1[1]
    }
  }

  return( list("d0" = l0, "d1" = l1, "d2" = l2, "d3" = l3) )

}

