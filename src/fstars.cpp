#include <Rcpp.h>
#include <iostream>
#include <vector>
using namespace Rcpp;



// [[Rcpp::export]]
void get_factors_stars(List dimensions, bool curvilinear = false) {
  std::cout << Rcpp::as<int>(Rcpp::as<List>(dimensions[0])["from"]) << std::endl;
  std::cout << Rcpp::as<int>(Rcpp::as<List>(dimensions[0])["to"]) << std::endl;
  std::cout << Rcpp::as<int>(Rcpp::as<List>(dimensions[0])["offset"]) << std::endl;
  std::cout << Rcpp::as<int>(Rcpp::as<List>(dimensions[0])["delta"]) << std::endl;
}
