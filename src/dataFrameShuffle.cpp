#include <Rcpp.h>

//' @useDynLib ciber

// [[Rcpp::export]]
Rcpp::DataFrame dataFrameShuffle(Rcpp::DataFrame df)
{
  	std::shuffle(df.begin(),df.end());
    return df;
}
