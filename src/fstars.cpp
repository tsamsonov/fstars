#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <proj.h>

// [[Rcpp::export]]
void get_factors_stars(Rcpp::List dimensions, bool curvilinear = false) {
  Rcpp::List dim = dimensions[0];
  int from = dim["from"];
  int to = dim["to"];
  double offset = dim["offset"];
  double delta = dim["delta"];
  std::string wkt = Rcpp::as<Rcpp::List>(dim["refsys"])["wkt"];

  std::cout << from << std::endl;
  std::cout << to << std::endl;
  std::cout << offset << std::endl;
  std::cout << delta << std::endl;
}

// [[Rcpp::export]]
int test_proj() {
  PJ_CONTEXT *C;
  PJ *P;
  PJ* P_for_GIS;
  PJ_COORD a, b;

  /* or you may set C=PJ_DEFAULT_CTX if you are sure you will     */
  /* use PJ objects from only one thread                          */
  // C = proj_context_create();

  C = PJ_DEFAULT_CTX;

  proj_context_use_proj4_init_rules(PJ_DEFAULT_CTX, 1);

  P = proj_create_crs_to_crs (C,
                              "EPSG:4326",
                              "+proj=eck4", /* or EPSG:32632 */
  NULL);

  // if (0==P) {
  //   fprintf(stderr, "Oops\n");
  //   return 1;
  // }

  /* a coordinate union representing Copenhagen: 55d N, 12d E    */
  /* Given that we have used proj_normalize_for_visualization(), the order of
   /* coordinates is longitude, latitude, and values are expressed in degrees. */
   a = proj_coord (12, 55, 0, 0);

   /* transform to UTM zone 32, then back to geographical */
   b = proj_trans (P, PJ_FWD, a);
   printf ("easting: %.3f, northing: %.3f\n", b.enu.e, b.enu.n);
   b = proj_trans (P, PJ_INV, b);
   printf ("longitude: %g, latitude: %g\n", b.lp.lam, b.lp.phi);

   /* Clean up */
   proj_destroy (P);
   proj_context_destroy (C); /* may be omitted in the single threaded case */
   return 0;
}
