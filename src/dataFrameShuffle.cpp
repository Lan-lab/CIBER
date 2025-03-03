#include <Rcpp.h>
#include <random>

//' @useDynLib ciber

// [[Rcpp::export]]
Rcpp::DataFrame dataFrameShuffle(Rcpp::DataFrame df)
{
    std::mt19937 g(rd());
    std::shuffle(df.begin(),df.end(),g);
    return df;
}
