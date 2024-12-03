
.wrap_nested_basis <- function(b, P,Xth,xm, add_slope, unif, drop.coef){

  force(b); force(P);force(Xth); force(unif); force(drop.coef);force(xm)

  evalX <- function(x, deriv){
    if(unif & (min(x)<0|max(x)>1)) stop("x is not unif(0,1)")
    withCallingHandlers({

      o <- b$evalX(x = x, deriv = deriv) # Get raw basis & derivatives
      if(add_slope){ # Add a slope in the last column: 1st derivative is 1, rest 0.
        o$X0 <- cbind(o$X0, x)
        if(deriv){
          o$X1 <- cbind(o$X1, 1)
          if(deriv >= 2){
            for(ii in 2:deriv){
              o[[paste0("X", ii)]] <- cbind(o[[paste0("X", ii)]], 0)
            }
          }
        }
      }
      o <- lapply(o, function(X) X %*% P)       # Reparametrise
      if(drop.coef) {
        o <- lapply(o, function(X) X[,-1])       # Reparametrise
        o$X0 <- sweep(o$X0, MARGIN = 2, xm , "-")
        Xbo <- sweep((Xth$X0%*%P)[,-1],  MARGIN = 2, xm , "-")
        Xbo1 <- (Xth$X1%*%P)[,-1]
      }else{
        Xbo <- (Xth$X0%*%P)
        Xbo1 <- (Xth$X1%*%P)
        }
      if(!unif){

        o <- linextr(x = x, b = o, th = b$krange, # Linearly extrapolate: NOTE method = "smooth" won't work if add_slope == TRUE
                     Xbo = Xbo, Xbo1 = Xbo1, method = "simple")
      }
      o
    }, warning = function(w) {
      if (length(grep("there is \\*no\\* information about some basis coefficients", conditionMessage(w)))){
        invokeRestart("muffleWarning")
      }
    })
  }

  out <- list("evalX" = evalX,
              drop.coef = drop.coef)

  return(out)
}
