initialize.coef <- function(sm, lambda){
  #initialized alpha
  alpha <- sm$xt$si$alpha
  alpha[1] <- exp(alpha[1]); alpha[-1] <- -alpha[-1]
  na <- length(alpha)

  # compute the initial beta
  ss <- function(bt, y, X, S){
    if(sm$drop.coef){
      bet <- exp(bt)
    } else{ bet <- c(bt[1], exp(bt[-1]))}
    res <- y-X%*%bet
    pen <- t(bt)%*%S%*%bt
    sum(res^2)/length(y) + pen
  }
  yy <- log(-log(1-rank(sm$xt$si$X%*%alpha)/(nrow(sm$X)+1)))
  kk <- ncol(sm$X)-na
  S <- sm$S[[1]][-(1:na), -(1:na)]
  b1 <- optim(fn = ss, par = c(1, rnorm(kk-1)), y = yy,
              X = sm$X[ , -(1:na)], S = lambda*S,
              method = "L-BFGS-B", lower = -50, upper = 100)
  # c(sm$xt$si$alpha, runif(ncol(sm$X)-na))
  c(sm$xt$si$alpha, b1$par)
}

