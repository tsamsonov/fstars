#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <proj.h>

using namespace std;

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
   PJ_CONTEXT *C = PJ_DEFAULT_CTX;
   PJ_COORD a, b, c, d;

   /* or you may set C=PJ_DEFAULT_CTX if you are sure you will     */
   /* use PJ objects from only one thread                          */
   // C = proj_context_create();

   // proj_context_use_proj4_init_rules(PJ_DEFAULT_CTX, 1);

   const char* prj = "+proj=moll";
   PJ* crs = proj_create(C, (string(prj) + " +type=crs").c_str());
   PJ* P = proj_create(C, prj);

   auto ell = proj_get_ellipsoid(C, crs);

   double out_semi_major_metre;
   double out_semi_minor_metre;
   int out_is_semi_minor_computed;
   double out_inv_flattening;

   proj_ellipsoid_get_parameters(C, ell,
                                 &out_semi_major_metre,
                                 &out_semi_minor_metre,
                                 &out_is_semi_minor_computed,
                                 &out_inv_flattening);
                                 cout << "ELLIPSOID:" << endl;
                                 cout << "out_semi_major_metre: "  << out_semi_major_metre << endl;
                                 cout << "out_semi_minor_metre: " << out_semi_minor_metre << endl;
                                 cout << "out_is_semi_minor_computed: " << out_is_semi_minor_computed << endl;
                                 cout << "out_inv_flattening: " << out_inv_flattening << endl << endl;

   cout << "COORDINATES:" << endl;
   double lat = 60, lon = 100;

   /* a coordinate union representing Copenhagen: 55d N, 12d E    */
   /* Given that we have used proj_normalize_for_visualization(), the order of
   /* coordinates is longitude, latitude, and values are expressed in degrees. */

   a.lp.lam = proj_torad(lon);
   a.lp.phi = proj_torad(lat);

   /* transform to UTM zone 32, then back to geographical */
   b = proj_trans(P, PJ_FWD, a);
   cout << "easting: " << b.enu.e << ", northing: " << b.enu.n << endl;

   c = proj_trans(P, PJ_INV, b);
   cout << "longitude: " << c.lp.lam << ", latitude: " << c.lp.phi << endl << endl;

   PJ_FACTORS pf = proj_factors(P, a);

   cout << "PROJECTION FACTORS:" << endl;
   cout << "meridional_scale: "  << pf.meridional_scale << endl;
   cout << "parallel_scale: " << pf.parallel_scale << endl;
   cout << "areal_scale: " << pf.areal_scale << endl;

   cout << "angular_distortion: " << pf.angular_distortion << endl;
   cout << "meridian_parallel_angle: " << pf.meridian_parallel_angle << endl;
   cout << "meridian_convergence: " << pf.meridian_convergence << endl;

   cout << "tissot_semimajor: " << pf.tissot_semimajor << endl;
   cout << "tissot_semiminor: " << pf.tissot_semiminor << endl;

   cout << "dx_dlam: " << pf.dx_dlam << endl;
   cout << "dx_dphi: " << pf.dx_dphi << endl;
   cout << "dy_dlam: " << pf.dy_dlam << endl;
   cout << "dy_dphi: " << pf.dy_dphi << endl << endl;

   /* Clean up */
   proj_destroy (P);
   proj_context_destroy (C); /* may be omitted in the single threaded case */
   return 0;
}
