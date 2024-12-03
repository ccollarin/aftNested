# gcv.aft_nl <- function(o,X){
gcv.aft_nl <- function(o, llk, kde){
  # y <- X[,1]
  # xa <- X[,-1]%*%o$coefficients[2:ncol(X)]/exp(o$coefficients[1])
  #
  # pp <- rank(y-xa)/length(y)
  # ff <- get.si(o, X = X)
  # qq <- 1-exp(ff$s)
  # sum(pp*log(pp/qq))

  llk(o$coefficients, deriv = 0, kde = kde)
}
