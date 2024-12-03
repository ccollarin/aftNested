#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(order_cpp)]]
List order_(NumericVector x) {
  if (is_true(any(duplicated(x)))) {
    Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(x).sort();

  IntegerVector ord = match(sorted, x)-1;
  IntegerVector rank = match(x, sorted)-1;

  return List::create(Named("order", ord),
                      Named("rank", rank));
}
