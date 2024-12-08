lambda.optim <- function(eff, sm,delta, vpen,mpen, env,val,delta.val,
                         trace=FALSE,init.par, ...){

  force(eff); force(sm); force(vpen);force(mpen)
  force(delta);force(env);force(val); force(delta.val)
  env$coef <- init.par
  env$out <- NULL

  eff.val <- eff_misi(.predict.matrix.si(sm, val, get.xa = TRUE)$Xi,
                      basis = sm$xt$basis)

  .optim_lambda <- function(rho){
    lambda <- exp(rho)
    llk <- llk_afc(eff, sm = sm, delta = delta, vpen = vpen, mpen = mpen,
                   lambda = lambda, ...)
    llk.val <- llk_afc(eff.val, sm = sm, delta = delta.val, mpen=mpen,
                       vpen = vpen, lambda = lambda, ...)

    out <- nlm(f = llk, p = env$coef, iterlim = 10000)

    out$lambda <- lambda
    sm$coefficients <- out$coefficients <- out$estimate; out$estimate <- NULL
    sm$xt$si$alpha <- out$coefficients[1:length(sm$xt$si$alpha)]
    out$smooth <- sm
    out$formula <- formula
    env$coef <- out$coefficients
    env$out <- out

    # loss <- gcv.aft_nl(out, X = .predict.matrix.si(sm, val, get.xa = TRUE)$Xi[delta.val,])
    loss <- gcv.aft_nl(out, llk.val, kde = eff$eval(sm$coefficients)$kde)
    attr(loss, "gradient") <- attr(loss, "hessian") <-NULL

    # Spen <- Reduce("+", lapply(1:length(sm$S), function(ii) lambda[ii]*sm$S[[ii]]))
    # evS <- eigen(Spen, only.values = TRUE)$values
    # H <- attr(llk(param = sm$coefficients), "hessian")
    # evH <- eigen(H, only.values = TRUE)$values
    #
    # laml <- c(out$minimum + 0.5 *( sum(log(evS[evS>1e-16])) -
    #   0.5 * sum(log(evH[evS>1e-16])) + sm$null.space.dim * log(2*pi) ))
    # print(laml)
    if(trace){
      cat(paste0("iter = ", out$iterations,
                 "\nGCV: ", loss, "  ||   llk = ", env$out$minimum,
                 "\nlambda = ", paste(lambda, collapse = "\t "), "\n\n"))
    }

    return(loss)
  }
  return(.optim_lambda)
}
