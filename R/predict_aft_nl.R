#' Prediction from fitted AFT.NL model
#'
#' @param object a fitted aft.nl object
#' @param type When this has value ``response'' the predicted time values are returned.
#'
#' @export
#'
predict.aft_nl <- function(object, type = "response"){
  if(type == "response"){
    eff <- eff_misi(object$smooth$xt$si$X, basis = object$smooth$xt$basis)
    eff <- eff$eval(object$coefficients, deriv = 2)
    X <- object$smooth$xt$si$X

    si <- get.si(object)

    dats <- data.frame(shat = exp(log.Stx(eff, deriv=0)$s),
                       x = si$xba)
    fit_x <- gam(x~s(shat), data = dats)

    exp(scale((fit_x$fitted.values - X[,-1]%*%si$alpha[-1]), scale=FALSE))
  }
}
