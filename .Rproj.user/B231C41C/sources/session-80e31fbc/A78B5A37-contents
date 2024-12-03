initial.lambda <- function(sm){
  Xtmp <- sm$X
  Xtmp[,1:length(sm$xt$si$alpha)] <- sm$xt$si$X

  lambda.out <- initial.sp(sm$X[, sm$first.para:sm$last.para],
                           list(sm$S[[1]][sm$first.para:sm$last.para,sm$first.para:sm$last.para]),
                           1)
  lambdas <- c(outer = lambda.out)
  if(length(sm$S)>1){
    si <- sm$xt$si
    lambda.int <- sapply(1:length(si$S), function(ii){
      fplp <- si$fplp[[ii]]
      initial.sp(si$X[, fplp],
                 list(si$S[[ii]][fplp, fplp]),
                 1)
    })
    lambdas["inner"] <- lambda.int
  }
  lambdas
}
