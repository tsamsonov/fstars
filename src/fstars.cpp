#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <map>
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include <proj.h>

using namespace std;

enum Method {
   DIRECT,
   AFFINE,
   TRAPEZOIDAL
};

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
   double x = di - i;
   double y = dj - j;
   return coef.a00(i, j) + coef.a01(i, j) * x + coef.a10(i, j) * y + coef.a11(i, j) * x * y;
}

// [[Rcpp::export]]
double rcpp_interpolate_xy(Rcpp::NumericMatrix matrix, Rcpp::List dimensions,
                          const double& x, const double& y) {
   vector<Dimension> dims;
   for (Rcpp::List dim: dimensions) {
      dims.emplace_back(dim);
   }

   auto coef = bilinear_coef(matrix);
   double dj = (y - dims[1].offset) / dims[1].delta;
   double di = (x - dims[0].offset) / dims[0].delta;

   if (di < 0 or dj < 0 or di >= coef.a00.nrow() or dj >= coef.a00.ncol()) {
      // cout << di << ' ' << dj << endl;
      return NA_REAL;
   }
   else
      return interpolate_ij(coef, di, dj);
}

vector<PJ_FACTORS> get_factors(const vector<Dimension>& dims,
                               const std::string& CRS,
                               const bool& curvilinear = false) {

   // cout << "TRYING FACTORS" << endl;

   int ifrom = dims[0].from;
   int ito = dims[0].to;
   double ioffset = dims[0].offset;
   double idelta = dims[0].delta;
   int ni = ito - ifrom + 1;

   int jfrom = dims[1].from;
   int jto = dims[1].to;
   double joffset = dims[1].offset;
   double jdelta = dims[1].delta;
   int nj = jto - jfrom + 1;

   int n = ni * nj;

   const char *prj = CRS.c_str();

   PJ_CONTEXT *C = PJ_DEFAULT_CTX;
   PJ* P = proj_create(C, prj);

   PJ_FACTORS pf;
   PJ_COORD pt, gpt;

   int idx = 0;

   vector<PJ_FACTORS> factors(n);
   PJ_FACTORS f;

   // cout << "READY TO CYCLE" << endl;

   for (auto i = 0; i < ni; ++i) {
      for (auto j = 0; j < nj; ++j) {
         idx = ni * j + i;

         pt.enu.e = ioffset + i * idelta;
         pt.enu.n = joffset + j * jdelta;

         gpt = proj_trans(P, PJ_INV, pt);

         // cout << idx << endl;
         // cout << i << endl;
         // cout << j << endl;
         // cout << pt.enu.e << endl;
         // cout << pt.enu.n << endl;
         // cout << gpt.lp.lam << endl;
         // cout << gpt.lp.phi << endl << endl;

         // if (idx == 10277) {
         //    cout << gpt.lp.lam << endl;
         //    cout << gpt.lp.phi << endl;
         // }

         factors[idx] = proj_factors(P, gpt);
      }
   }

   return factors;
}

// [[Rcpp::export]]
Rcpp::List rcpp_get_factors(Rcpp::List dimensions, const std::string& CRS,
                            const bool& curvilinear = false) {
   vector<Dimension> dims;
   for (Rcpp::List dim: dimensions) {
      dims.emplace_back(dim);
   }

   auto factors = get_factors(dims, CRS, curvilinear);

   int ni = dims[0].to - dims[0].from + 1;
   int nj = dims[1].to - dims[1].from + 1;
   int n = ni * nj;

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

   double tsin, tcos, k2, h2, dx, dy;

   int idx;

   for (auto i = 0; i < ni; ++i) {
      for (auto j = 0; j < nj; ++j) {
         idx = nj * i + j;

         meridional_scale[idx] = factors[idx].meridional_scale;
         parallel_scale[idx] = factors[idx].parallel_scale;
         areal_scale[idx] = factors[idx].areal_scale;
         angular_distortion[idx] = factors[idx].angular_distortion;
         meridian_parallel_angle[idx] = factors[idx].meridian_parallel_angle;
         meridian_convergence[idx] = factors[idx].meridian_convergence;
         tissot_semimajor[idx] = factors[idx].tissot_semimajor;
         tissot_semiminor[idx] = factors[idx].tissot_semiminor;
         dx_dlam[idx] = factors[idx].dx_dlam;
         dx_dphi[idx] = factors[idx].dx_dphi;
         dy_dlam[idx] = factors[idx].dy_dlam;
         dy_dphi[idx] = factors[idx].dy_dphi;

         // ADDITIONAL

         // Parallel convergence
         parallel_convergence[idx] = atan2(factors[idx].dy_dlam, factors[idx].dx_dlam);

         // Tissot ellipse orientation
         tsin = sin(2.0 * meridian_parallel_angle[idx]);
         tcos = cos(2.0 * meridian_parallel_angle[idx]);
         k2 =  parallel_scale[idx] * parallel_scale[idx];
         h2 =  meridional_scale[idx] * meridional_scale[idx];
         dx = h2 * tsin;
         dy = k2 + h2 * tcos;

         tissot_orientation[idx] = 0.5 * atan2(dx, dy);
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

double east_to_geo(const double& az) {
   if (az >= 0 && az <= M_PI_2) {
      return M_PI_2 - az;
   } else if (az > M_PI_2){
      return M_PI + M_PI + M_PI_2 - az;
   } else {
      return M_PI_2 - az;
   }
}

double north_to_geo(const double& az) {
   if (az < 0) {
      return M_PI + M_PI + az;
   } else {
      return az;
   }
}

double to_closest(const double& az, const double& az_ref) {
   auto dif = az - az_ref;
   if (abs(dif) < M_PI_2) {
      return az;
   } else if (dif < 0) {
      return az + M_PI;
   } else {
      return az - M_PI;
   }
}

double to_near(const double& az, const double& ref) {

}

std::pair<Rcpp::NumericMatrix, Rcpp::NumericMatrix> get_xy_kernel(const int& i,
                                                                  const int& j,
                                                                  const Rcpp::NumericMatrix& ishift,
                                                                  const Rcpp::NumericMatrix& jshift,
                                                                  const vector<Dimension>& dims,
                                                                  const vector<PJ_FACTORS>& pf,
                                                                  const double& dfactor = 1.0,
                                                                  const Method& = DIRECT) {
   auto nrow = ishift.nrow();
   auto ncol = ishift.ncol();

   Rcpp::NumericMatrix x(nrow, ncol);
   Rcpp::NumericMatrix y(nrow, ncol);

   auto idx = j * (dims[0].to - dims[0].from + 1) + i;

   auto lambdaScale = pf[idx].meridional_scale * pf[idx].tissot_semimajor / dfactor;
   auto phiScale = pf[idx].parallel_scale * pf[idx].tissot_semimajor / dfactor;
   auto parallel_convergence = atan2(pf[idx].dy_dlam, pf[idx].dx_dlam);

   int di, dj;
   double D, A, a, mu;

   for (auto k = 0; k < nrow; ++k) {
      for (auto l = 0; l < ncol; ++l) {
         di = ishift(k, l);
         dj = jshift(k, l);
         if (di == 0 and dj == 0) {
            x(k, l) = 0;
            y(k, l) = 0;
         } else if (dj == 0) {
            x(k, l) = di * dims[0].delta * phiScale * cos(parallel_convergence);
            y(k, l) = dj * dims[1].delta * phiScale * sin(parallel_convergence);
         } else if (di == 0) {
            x(k, l) = di * dims[0].delta * lambdaScale * sin(pf[idx].meridian_convergence);
            y(k, l) = dj * dims[1].delta * lambdaScale * cos(pf[idx].meridian_convergence);
         } else {
            D = sqrt(pow(di * dims[0].delta, 2) + pow(dj * dims[1].delta, 2));
            A = east_to_geo(atan2(dj * dims[1].delta, di * dims[0].delta));
            a = to_closest(north_to_geo(atan(phiScale * sin(pf[idx].meridian_parallel_angle) * tan(A) /
                           (lambdaScale + phiScale * cos(pf[idx].meridian_parallel_angle) * tan(A)))), A);
            mu = sqrt(pow(lambdaScale, 2) * pow(cos(A), 2) +
                      lambdaScale * phiScale * cos(pf[idx].meridian_parallel_angle) * sin(2 * A) +
                      pow(phiScale, 2) * pow(sin(A), 2));

            // cout << D << ' ' << A << ' ' << mu * D << ' ' << a << endl;

            x(k, l) = D * mu * sin(a - pf[idx].meridian_convergence);
            y(k, l) = D * mu * cos(a - pf[idx].meridian_convergence);
         }

         x(k, l) += i * dims[0].delta + dims[0].offset;
         y(k, l) += j * dims[1].delta + dims[1].offset;
      }
   }
   // cout << endl;

   return std::pair(x, y);

}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_filter_matrix(const Rcpp::NumericMatrix&  matrix,
                                       const Rcpp::NumericMatrix&  kernel,
                                       const Rcpp::List& dimensions,
                                       const std::string& CRS,
                                       const bool& curvilinear = false,
                                       const bool& adaptive = false) {

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
   double penalty;
   int idx, jdx, ikdx, jldx;

   if (adaptive) {
      // cout << "ADAPTIVE" << endl;
      vector<Dimension> dims;
      for (Rcpp::List dim: dimensions) {
         dims.emplace_back(dim);
      }
      auto factors = get_factors(dims, CRS, curvilinear);
      double max_semimajor = 0.0;
      for (auto f : factors) {
         if (f.tissot_semimajor > max_semimajor){
            max_semimajor = f.tissot_semimajor;
         }
      }

      // cout << "FACTORS" << endl;
      auto coef = bilinear_coef(matrix);
      // cout << "BILINEARS" << endl;
      double di, dj, value;

      // cout << "START FILTERING" << endl;

      for (auto i = 0; i < ni; ++i) {
         for (auto j = 0; j < nj; ++j) {
            idx = i + istart;
            jdx = j + jstart;
            if (nodata(idx, jdx) == 1) {
               res(i, j) = NA_REAL;
            } else {
               res(i, j) = 0;
               penalty = 0;
               auto [x, y] = get_xy_kernel(idx, jdx, ishift, jshift, dims, factors, max_semimajor);

               for (auto k = 0; k < nk; ++k) {
                  for (auto l = 0; l < nl; ++l) {
                     di = (x(k, l) - dims[0].offset) / dims[0].delta;
                     dj = (y(k, l) - dims[1].offset) / dims[1].delta;

                     // cout << di << ' ' << dj << endl;

                     int idi = floor(di);
                     int idj = floor(dj);

                     if (nodata(idi, idj) + nodata(idi, idj + 1) + nodata(idi + 1, idj) + nodata(idi + 1, idj + 1) > 0) {
                        penalty += kernel(k, l);
                     } else {
                        value = interpolate_ij(coef, di, dj);
                        res(i, j) += value * kernel(k, l);
                     }

                  }
               }

               // cout << endl;


               res(i, j) = res(i, j) / (ksum - penalty);
            }
         }
      }
   } else {
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

   }

   return res;
}

// [[Rcpp::export]]
int test_proj() {
   PJ_CONTEXT *C = PJ_DEFAULT_CTX;
   PJ_COORD a, b, c;

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
