df.aft_nl <- function(o,llk){

  H <- attr(llk(o$coefficients),"hessian")
  Spen <- Reduce("+", lapply(1:length(o$smooth$S), function(ii) o$lambda[ii]*o$smooth$S[[ii]]))

  pdf <- diag(solve(-H)%*%(-H+2*Spen))
  df <- sum(pdf)
  df
}
