#include <RcppArmadillo.h>
#include "order_cpp.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.gauss_kce_cpp)]]
List gauss_kce(NumericVector x, arma::mat Xi,int deriv,
          NumericVector x_fix, arma::mat Xi_fix, double h ) {

  int m = Xi_fix.n_rows;
  int n = Xi.n_rows;
  int p = Xi.n_cols;
  float ccut = R::qnorm(1, 0,h,true, false );
  // float ccut = R::qnorm((1-1e-14), 0,h,true, false );
  NumericVector dist_i(m);
  mat dist_Xic(m,p);
  colvec u1(m);
  colvec dist(m);

  colvec d0(n);
  mat d1, d2;

  if (deriv) {
    d1 = mat(n, p);
    if (deriv > 1) {
      d2 = mat(n, p * (p + 1) / 2);
    }
  }

  List ind_fix = order_(x_fix);
  IntegerVector ord_fix = ind_fix[0];
  IntegerVector rank_fix = ind_fix[1];

  List ind = order_(x);
  IntegerVector ord = ind["order"];
  IntegerVector rank = ind["rank"];

  int start_i = 0;
  for (int ii = 0; ii < n; ++ii) {

    bool upstart = FALSE;
    int zz = 0;
    double xii = x[ord[ii]];

    for(int ix = start_i; ix < m; ++ix){

      double xfii = x_fix[ord_fix[ix]];
      double dd = (xii-xfii)/h;

      if((dd<ccut) & (dd>-ccut)){
        dist_i[zz] = dd;
        if(!upstart){
          start_i = ix;
          upstart = TRUE;
        }
        if(deriv>0){
          for (int ic =0; ic < p; ++ic) {
            dist_Xic(zz, ic) = (Xi(ord[ii],ic)-Xi_fix(ord_fix[ix],ic))/(m*h);
          }
        }
        zz += 1;
      } else if(upstart){
        break;
      }
    }
    if(zz==0){ // compute all distances (to get correct derivatives)
      cout<<"in smaller\n";
      start_i = 0;
      for(int ix = start_i; ix < m; ++ix){

        double xfii = x_fix[ord_fix[ix]];
        double dd = (xii-xfii)/h;

        dist_i[zz] = dd;
        if(deriv>0){
          for (int ic =0; ic < p; ++ic) {
            dist_Xic(zz, ic) = (Xi(ord[ii],ic)-Xi_fix(ord_fix[ix],ic))/(m*h);
          }
        }
        zz += 1;
      }
    }//end if zz=0

    d0(ord[ii]) = sum(Rcpp::pnorm(dist_i[Rcpp::seq(0, zz-1)]))/m + (double)start_i/(double)m;

    if(deriv>0){

      int kk = 0;

      dist = dist_i[Rcpp::seq(0, zz-1)];
      u1 = Rcpp::dnorm(dist_i[Rcpp::seq(0, zz-1)]);

      for (int ic =0; ic < p; ++ic) {
        d1(ord[ii], ic) = dot(dist_Xic.rows(0,zz-1).col(ic),  u1);
        // d1(ii, ic) = dot(dist_Xic.col(ic),  u1);
        if(deriv>1){
          for (int jc=ic; jc<p; ++jc) {
            d2(ord[ii],kk) = dot(u1 % dist,
               dist_Xic.rows(0,zz-1).col(ic) % dist_Xic.rows(0,zz-1).col(jc))*m;
              // d2(ii,kk) = dot(u1 % dist, dist_Xic.col(ic) % dist_Xic.col(jc))*m;
            kk += 1;
          } // end for jc
        }//end if deriv>1
      }// end loop ic
    } // end deriv>0
  }// end loop ii

  return List::create(Named("d0", d0), Named("d1", d1), Named("d2", d2));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.triweight_kce_cpp)]]
List triweight_kce(NumericVector x, arma::mat Xi,int deriv,
                   NumericVector x_fix, arma::mat Xi_fix, double h) {
  int n = Xi.n_rows;
  int m_fix = Xi_fix.n_rows;
  int p = Xi.n_cols;
  colvec u1(m_fix);
  colvec dist(m_fix);

  colvec d0(n);
  mat d1, d2;

  if (deriv) {
    d1 = mat(n, p);
    if (deriv > 1) {
      d2 = mat(n, p * (p + 1) / 2);
    }
  }

  List ind_fix = order_(x_fix);
  IntegerVector ord_fix = ind_fix[0];
  IntegerVector rank_fix = ind_fix[1];

  List ind = order_(x);
  IntegerVector ord = ind["order"];
  IntegerVector rank = ind["rank"];

  NumericVector dist_i(m_fix);
  mat dist_Xic(m_fix,p);

  double start_i = 0;
  int stop_i;
  // int m;
  // bool upstop = FALSE;

  for (int ii = 0; ii < n; ++ii) {
    bool upstart = FALSE;
    double xii = x[ord[ii]];

    int zz = 0;
    for(int ix = start_i; ix < m_fix; ++ix){
      double xfii = x_fix[ord_fix[ix]];
      double dd = (xii-xfii)/h;

      if((dd<1) & (dd>-1)){
        dist_i[zz] = dd;
        if(!upstart){
          start_i = ix;
          upstart = TRUE;
        }
        if(deriv>0){
          for (int ic =0; ic < p; ++ic) {
            dist_Xic(zz, ic) = (Xi(ord[ii],ic)-Xi_fix(ord_fix[ix],ic))/(m_fix*h);
          }
        }
        zz += 1;
      } else if(upstart){
        stop_i = ix+1;
        break;
      }
    }
    NumericVector dq;
    if(zz==0){
      cout<<"z è 0";
      d0(ord[ii]) = (double)stop_i/(double)m_fix; // might be changed with sum(x_fix<xii)/m_fix
    }else{
      dq = dist_i[Rcpp::seq(0, zz-1)]*dist_i[Rcpp::seq(0, zz-1)];
      d0(ord[ii]) = sum(0.5 + 35/32*(dist_i[Rcpp::seq(0, zz-1)]*
        (1+dq*(-1+dq*(0.6-dq/7)))))/m_fix + (double)start_i/(double)m_fix;
    }

    if(deriv>0){
      int kk = 0;

      u1 = 35/32*(1-dq)*(1-dq)*(1-dq);
      // cout<<u1<<"\n";

      for (int ic = 0; ic < p; ++ic) {
        if(zz==0){
          d1(ord[ii], ic) = 0;
        } else{
          d1(ord[ii], ic) = dot(dist_Xic.rows(0,zz-1).col(ic),  u1);
        }

        if(deriv>1){
          for (int jc=ic; jc<p; ++jc) {
            if(zz==0){
              d2(ord[ii],kk) = 0;
            } else{
              dist = dist_i[Rcpp::seq(0, zz-1)];
              d2(ord[ii],kk) = dot(u1 % dist,
                 dist_Xic.rows(0,zz-1).col(ic) % dist_Xic.rows(0,zz-1).col(jc))*m_fix;
            }
            kk += 1;
          } // end for jc
        }//end if deriv>1
      }// end loop ic
    }
  }

  return List::create(Named("d0", d0), Named("d1", d1), Named("d2", d2));
}

// else if(zz==1){
//   double t2 = dist_Xic(0,ic) * dist_Xic(0,jc);
//   d2(ord[ii],kk) = u1[0] * dist_i[0]*t2*m_fix;
// }
