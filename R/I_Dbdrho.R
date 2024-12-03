.Dbdrho <- function(o, llk, lambda){

  sm <- o$smooth
  H1 <- solve(attr(llk, "hessian"))
  Spen <- Reduce("+", lapply(1:length(sm$S), function(ii) lambda[ii]*sm$S[[ii]])) # Outer and inner ridge penalties
  theta <- o$coefficients

  lambda * (H1 %*% (Spen %*% theta))
}
