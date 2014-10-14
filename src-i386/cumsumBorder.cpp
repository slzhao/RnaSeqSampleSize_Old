#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]                                                             
NumericVector cumsumBorder(NumericVector x, double border) {
  NumericVector res(1);
  double acc = 0;
  for (int i=0; i < x.size(); ++i) {
    res=i-1;
    acc += x[i];
    if (acc > border)  return res;
  }
  return 0;
}