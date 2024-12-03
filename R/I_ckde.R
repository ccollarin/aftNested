.ckde <- function(x_fix,Xi_fix = NULL, h=NULL, kernel="triweigth"){
  force(x_fix);force(h);force(Xi_fix); force(kernel)
  n <- length(x_fix)

  if(is.null(h)){
    h <- (4/(3*n))^(1/5)*sd(x_fix)
  }

  eval <- function(x = NULL, Xi = NULL, deriv = 0){
      if(is.null(Xi)) Xi <- Xi_fix
      if(is.null(x)) x <- x_fix

      if(substr(kernel, 1,1) == "g"){
        return(.gauss_kce_cpp(x = x, Xi = Xi,deriv = deriv,
                x_fix = x_fix, Xi_fix = Xi_fix, h=h))
      }
      if(substr(kernel, 1,1) == "t"){
        return(.triweight_kce_cpp(x = x, Xi = Xi,deriv = deriv,
                x_fix = x_fix, Xi_fix = Xi_fix, h=h))

      }
  }

  return(eval)
}

  # eval <- function(x = NULL, Xi = NULL, deriv = 0){
  #   o <- list()
  #   if(is.null(Xi)) Xi <- Xi_fix
  #   if(is.null(x)) x <- x_fix
  #   o$d0 <- sapply(x, function(a) sum(pnorm((a-x_fix)/h))/n)
  #
  #   if(deriv>0){
  #     # allocate storage
  #     o$d1 <- matrix(nrow=length(x), ncol=p)
  #
  #     if(deriv>1){
  #       o$d2 <- matrix(nrow=length(x), ncol=p*(p+1)/2)
  #     }
  #
  #     for(ir in 1:length(x)){
  #       zz <- 1
  #       dist <- x[ir]-x_fix
  #       u1 <- dnorm(dist/h)
  #
  #       for (ic in 1:p) {
  #         dist_Xic <- Xi[ir,ic]-Xi_fix[,ic]
  #         o$d1[ir,ic] <- crossprod(u1, dist_Xic)/(n*h)
  #
  #         if(deriv>1){
  #           for (jc in ic:p) {
  #             dist_Xjc <- Xi[ir,jc]-Xi_fix[,jc]
  #             o$d2[ir,zz] <- crossprod(u1 * dist, dist_Xic * dist_Xjc)/(n*h^3)
  #             zz <- zz+1
  #           } # end for jc
  #         }#end if deriv>1
  #       }# end loog ic
  #     }# end loop ir
  #   }# end if deriv>0
  #
  #   return(o)
  # }
      # if(deriv>1){
      #   o$d2 <- do.call("rbind",
      #                   lapply(1:length(x),
      #                          function(ii){
      #                            dist <- x[ii]-x_fix
      #                            dist_Xi <-do.call("cbind",
      #                                              lapply(1:p, function(jj) Xi[ii,jj]-Xi_fix[,jj]))
      #                            colSums(drop(dnorm(dist/h))*dist_Xi)/(n*h)
      #                            # sum(drop(dnorm(dist/h)))*Xi[ii,]/(n*h)
      #                          }))
      # }
