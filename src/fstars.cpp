#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <proj.h>

using namespace std;

union Values {
   Rcpp::NumericVector rect; // rectilinear raster
   Rcpp::NumericMatrix curv; // curvilinear raster
};

struct Dimension {
   int from;
   int to;
   double offset;
   double delta;
   const char* refsys;
   bool point;
   // Values values;
   Dimension(Rcpp::List dim) {
      from = dim["from"];
      to = dim["to"];
      offset = dim["offset"];
      delta = dim["delta"];
      refsys = dim["refsys"];
      point = dim["point"];
      // values.curv = Rcpp::as<Rcpp::NumericMatrix>(dim["values"]);
   }
};

struct Bilinear {
   Rcpp::NumericMatrix a00, a01, a10, a11;
   Bilinear(int ni, int nj) {
      a00 = Rcpp::NumericMatrix(ni, nj);
      a01 = Rcpp::NumericMatrix(ni, nj);
      a10 = Rcpp::NumericMatrix(ni, nj);
      a11 = Rcpp::NumericMatrix(ni, nj);
   }
};

Bilinear bilinear_coef(Rcpp::NumericMatrix  matrix) {
   int ni = matrix.nrow() - 1;
   int nj = matrix.ncol() - 1;
   Bilinear coef(ni, nj);

   double z00, z10, z01, z11;

   for (auto i = 0; i < ni; ++i) {
      for (auto j = 0; j < nj; ++j) {
         z00 = matrix(i, j);
         z10 = matrix(i+1, j);
         z01 = matrix(i, j+1);
         z11 = matrix(i+1, j+1);
         coef.a00(i, j) = z00;
         coef.a01(i, j) = z01 - z00;
         coef.a10(i, j) = z10 - z00;
         coef.a11(i, j) = z00 - z10 - z01 + z11;
      }
   }

   return coef;
}

double interpolate_ij(const Bilinear& coef, const double& di, const double& dj) {
   int i = floor(di);
   int j = floor(dj);
   double y = di - i;
   double x = dj - j;
   return coef.a00(i, j) + coef.a01(i, j) * x + coef.a10(i, j) * y + coef.a11(i, j) * x * y;
}

// [[Rcpp::export]]
double cpp_interpolate_xy(Rcpp::NumericMatrix matrix, Rcpp::List dimensions,
                          const double& x, const double& y) {
   vector<Dimension> dims;
   for (Rcpp::List dim: dimensions) {
      dims.emplace_back(dim);
   }

   auto coef = bilinear_coef(matrix);
   double dj = (y - dims[1].offset) / dims[1].delta;
   double di = (x - dims[0].offset) / dims[0].delta;

   if (di < 0 or dj < 0 or di >= coef.a00.nrow() or dj >= coef.a00.ncol()) {
      cout << di << ' ' << dj << endl;
      return NA_REAL;
   }
   else
      return interpolate_ij(coef, di, dj);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix filter_matrix(Rcpp::NumericMatrix  matrix,
                                  Rcpp::NumericMatrix  kernel) {
   int nk = kernel.nrow();
   int nl = kernel.ncol();

   int ni = matrix.nrow() - nk + 1;
   int nj = matrix.ncol() - nl + 1;

   Rcpp::NumericMatrix ishift(nk, nl);
   Rcpp::NumericMatrix jshift(nk, nl);

   int istart = nk / 2;
   int jstart = nl / 2;

   int di = -istart;
   int dj;

   double ksum = 0;

   for (auto k = 0; k < nk; ++k) {
      dj = -jstart;
      for (auto l = 0; l < nl; ++l) {
         ishift(k, l) = di;
         jshift(k, l) = dj;
         ksum += kernel(k, l);
         dj++;
      }
      di++;
   }

   Rcpp::NumericMatrix nodata(matrix.nrow(), matrix.ncol());

   for (auto i = 0; i < matrix.nrow(); ++i)
      for (auto j = 0; j < matrix.ncol(); ++j)
         nodata(i, j) = to_string(matrix(i, j)) == "nan";

   Rcpp::NumericMatrix res(ni, nj);
   double val, penalty;
   int idx, jdx, ikdx, jldx;

   for (auto i = 0; i < ni; ++i) {
      for (auto j = 0; j < nj; ++j) {
         idx = i + istart;
         jdx = j + jstart;
         if (nodata(idx, jdx) == 1) {
            res(i, j) = NA_REAL;
         } else {
            res(i, j) = 0;
            penalty = 0;
            for (auto k = 0; k < nk; ++k) {
               for (auto l = 0; l < nl; ++l) {
                  ikdx = idx + ishift(k, l);
                  jldx = jdx + jshift(k, l);
                  if (nodata(ikdx, jldx) == 1) {
                     penalty += kernel(k, l);
                  } else {
                     res(i, j) += matrix(ikdx, jldx) * kernel(k, l);
                  }
               }
            }
            res(i, j) = res(i, j) / (ksum - penalty);
         }
      }
   }

   return res;
}

// [[Rcpp::export]]
Rcpp::List get_factors_stars(Rcpp::List dimensions, std::string CRS, bool curvilinear = false) {

   Rcpp::List idim = dimensions[1];
   int ifrom = idim["from"];
   int ito = idim["to"];
   double ioffset = idim["offset"];
   double idelta = idim["delta"];
   int ni = ito - ifrom + 1;

   Rcpp::List jdim = dimensions[0];
   int jfrom = jdim["from"];
   int jto = jdim["to"];
   double joffset = jdim["offset"];
   double jdelta = jdim["delta"];
   int nj = jto - jfrom + 1;

   int n = ni * nj;

   const char *prj = CRS.c_str();

   PJ_CONTEXT *C = PJ_DEFAULT_CTX;
   PJ* P = proj_create(C, prj);

   PJ_FACTORS pf;
   PJ_COORD pt, gpt;

   vector<double> meridional_scale(n, 0);
   vector<double> parallel_scale(n, 0);
   vector<double> areal_scale(n, 0);
   vector<double> angular_distortion(n, 0);
   vector<double> meridian_parallel_angle(n, 0);
   vector<double> meridian_convergence(n, 0);
   vector<double> parallel_convergence(n, 0);
   vector<double> tissot_semimajor(n, 0);
   vector<double> tissot_semiminor(n, 0);
   vector<double> tissot_orientation(n, 0);
   vector<double> dx_dlam(n, 0);
   vector<double> dx_dphi(n, 0);
   vector<double> dy_dlam(n, 0);
   vector<double> dy_dphi(n, 0);

   int idx = 0;

   double sintheta, tsin, tcos, dir, a2, b2, k2, h2, dx, dy, theta;

   for (auto i = 0; i < ni; ++i) {
      for (auto j = 0; j < nj; ++j) {
         idx = nj * i + j;
         pt.enu.n = ioffset + i * idelta;
         pt.enu.e = joffset + j * jdelta;
         gpt = proj_trans(P, PJ_INV, pt);
         pf = proj_factors(P, gpt);

         meridional_scale[idx] = pf.meridional_scale;
         parallel_scale[idx] = pf.parallel_scale;
         areal_scale[idx] = pf.areal_scale;
         angular_distortion[idx] = pf.angular_distortion;
         meridian_convergence[idx] = pf.meridian_convergence;
         meridian_parallel_angle[idx] = pf.meridian_parallel_angle;
         tissot_semimajor[idx] = pf.tissot_semimajor;
         tissot_semiminor[idx] = pf.tissot_semiminor;
         dx_dlam[idx] = pf.dx_dlam;
         dx_dphi[idx] = pf.dx_dphi;
         dy_dlam[idx] = pf.dy_dlam;
         dy_dphi[idx] = pf.dy_dphi;

         parallel_convergence[idx] = atan2(pf.dy_dlam, pf.dx_dlam);

         // Tissot ellipse orientation

         // theta = (pf.meridian_convergence > 0) ? M_PI_2 + pf.meridian_convergence : pf.meridian_parallel_angle;
         tsin = sin(2.0 * pf.meridian_parallel_angle);
         tcos = cos(2.0 * pf.meridian_parallel_angle);
         k2 =  pf.parallel_scale * pf.parallel_scale;
         h2 =  pf.meridional_scale * pf.meridional_scale;
         // a2 = pf.tissot_semimajor * pf.tissot_semimajor;
         // b2 = pf.tissot_semiminor * pf.tissot_semiminor;
         dx = h2 * tsin;
         dy = k2 + h2 * tcos;

         dir =  0.5 * atan2(dx, dy);
         tissot_orientation[idx] = (pf.meridian_convergence > 0) ? 0.5 * atan2(dx, dy) : -0.5 * atan2(dx, dy);

         // tissot_orientation[idx] = sqrt((a2 - h2) / (h2 - b2)) * pf.tissot_semiminor / pf.tissot_semimajor;


      }
   }

   return Rcpp::List::create(Rcpp::Named("meridional_scale") = meridional_scale,
                             Rcpp::Named("parallel_scale") = parallel_scale,
                             Rcpp::Named("areal_scale") = areal_scale,
                             Rcpp::Named("angular_distortion") = angular_distortion,
                             Rcpp::Named("meridian_parallel_angle") = meridian_parallel_angle,
                             Rcpp::Named("meridian_convergence") = meridian_convergence,
                             Rcpp::Named("parallel_convergence") = parallel_convergence,
                             Rcpp::Named("tissot_semimajor") = tissot_semimajor,
                             Rcpp::Named("tissot_semiminor") = tissot_semiminor,
                             Rcpp::Named("tissot_orientation") = tissot_orientation,
                             Rcpp::Named("dx_dlam") = dx_dlam,
                             Rcpp::Named("dx_dphi") = dx_dphi,
                             Rcpp::Named("dy_dlam") = dy_dlam,
                             Rcpp::Named("dy_dphi") = dy_dphi);

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
