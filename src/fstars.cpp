#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <map>
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include <chrono>
#include <thread>
#include <proj.h>

using namespace std;

const double M_7_PI_4 = 7 * M_PI_4;

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

   Dimension(const Rcpp::List& dim) {
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
   Bilinear(const int& ni, const int& nj) {
      a00 = Rcpp::NumericMatrix(ni, nj);
      a01 = Rcpp::NumericMatrix(ni, nj);
      a10 = Rcpp::NumericMatrix(ni, nj);
      a11 = Rcpp::NumericMatrix(ni, nj);
   }
};


Bilinear bilinear_coef(const Rcpp::NumericMatrix&  matrix) {
   int ni = matrix.nrow();
   int nj = matrix.ncol();
   Bilinear coef(ni, nj);
   double z00, z10, z01, z11;

   for (auto i = 0; i < ni-1; ++i) {
      for (auto j = 0; j < nj-1; ++j) {
         z00 = matrix(i,   j);
         z10 = matrix(i+1, j);
         z01 = matrix(i,   j+1);
         z11 = matrix(i+1, j+1);
         coef.a00(i, j) = z00;
         coef.a01(i, j) = z01 - z00;
         coef.a10(i, j) = z10 - z00;
         coef.a11(i, j) = z00 - z10 - z01 + z11;
      }
   }

   for (auto i = 0; i < ni-1; ++i) {
      z00 = matrix(i,   nj-1);
      z10 = matrix(i+1, nj-1);
      z01 = matrix(i,   nj-1);
      z11 = matrix(i+1, nj-1);
      coef.a00(i, nj-1) = z00;
      coef.a01(i, nj-1) = z01 - z00;
      coef.a10(i, nj-1) = z10 - z00;
      coef.a11(i, nj-1) = z00 - z10 - z01 + z11;
   }

   for (auto j = 0; j < nj-1; ++j) {
      z00 = matrix(ni-1, j);
      z10 = matrix(ni-1, j);
      z01 = matrix(ni-1, j+1);
      z11 = matrix(ni-1, j+1);
      coef.a00(ni-1, j) = z00;
      coef.a01(ni-1, j) = z01 - z00;
      coef.a10(ni-1, j) = z10 - z00;
      coef.a11(ni-1, j) = z00 - z10 - z01 + z11;
   }

   coef.a00(ni-1, nj-1) = matrix(ni-1, nj-1);
   coef.a01(ni-1, nj-1) = 0;
   coef.a10(ni-1, nj-1) = 0;
   coef.a11(ni-1, nj-1) = 0;

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
double rcpp_interpolate_xy(const Rcpp::NumericMatrix& matrix, const Rcpp::List& dimensions,
                          const double& x, const double& y) {
   vector<Dimension> dims;
   for (Rcpp::List dim: dimensions) {
      dims.emplace_back(dim);
   }

   auto coef = bilinear_coef(matrix);
   double dj = (y - dims[1].offset) / dims[1].delta;
   double di = (x - dims[0].offset) / dims[0].delta;

   if (di < 0 or dj < 0 or di >= coef.a00.nrow() or dj >= coef.a00.ncol()) {
      return NA_REAL;
   }
   else
      return interpolate_ij(coef, di, dj);
}

double rcpp_interpolate_xy(const Rcpp::NumericMatrix& matrix,
                           const vector<Dimension>& dims,
                           const Bilinear& coef,
                           const double& x, const double& y) {
   double dj = (y - dims[1].offset) / dims[1].delta;
   double di = (x - dims[0].offset) / dims[0].delta;

   if (di < 0 or dj < 0 or di >= coef.a00.nrow() or dj >= coef.a00.ncol()) {
      return NA_REAL;
   }
   else
      return interpolate_ij(coef, di, dj);
}

class ZevenbergenSurface {
   double A, B, C, D, E, F, G, H, I;
   double z1, z2, z3, z4, z5, z6, z7, z8, z9;
   double res;

   public:
      ZevenbergenSurface(const Rcpp::NumericMatrix& matrix,
                         const vector<Dimension>& dims, const Bilinear& coef,
                         const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& y,
                         const double& scale) {
         z1 = rcpp_interpolate_xy(matrix, dims, coef, x(0, 0), y(0, 0));
         z2 = rcpp_interpolate_xy(matrix, dims, coef, x(1, 0), y(1, 0));
         z3 = rcpp_interpolate_xy(matrix, dims, coef, x(2, 0), y(2, 0));
         z4 = rcpp_interpolate_xy(matrix, dims, coef, x(0, 1), y(0, 1));
         z5 = rcpp_interpolate_xy(matrix, dims, coef, x(1, 1), y(1, 1));
         z6 = rcpp_interpolate_xy(matrix, dims, coef, x(2, 1), y(2, 1));
         z7 = rcpp_interpolate_xy(matrix, dims, coef, x(0, 2), y(0, 2));
         z8 = rcpp_interpolate_xy(matrix, dims, coef, x(1, 2), y(1, 2));
         z9 = rcpp_interpolate_xy(matrix, dims, coef, x(2, 2), y(2, 2));

         res = scale * abs(dims[0].delta);

         cout << res<< endl;

         A = ((z1 + z3 + z7 + z9)/4.f - (z2 + z4 + z6 + z8)/2.f + z5) / pow(res, 4);
         B = ((z1 + z3 - z7 - z9)/4.f - (z2 - z8)/2.f) / pow(res, 3);
         C = ((-z1 + z3 - z7 + z9)/4.f + (z2 - z6)/2.f) / pow(res, 3);
         D = ((z4 + z6)/2.f - z5) / pow(res, 2);
         E = ((z2 + z8)/2.f - z5) / pow(res, 2);
         F = (-z1 + z3 + z7 - z9) / (4.f * pow(res, 2));
         G = (-z4 + z6) / (2.f * res);
         H = (z2 - z8) / (2.f * res);
         I = z5;
      }

      double slope() {
         return atan(sqrt(G*G + H*H));
      }

      double hillshade(double azimuth = M_7_PI_4, double height = M_PI_4, double zFactor = 1.0) {
         double Lx = cos(height) * cos(azimuth);
         double Ly = cos(height) * sin(azimuth);
         double Lz = sin(height);

         double Nx = G * zFactor;
         double Ny = H * zFactor;
         double Nz = 1;

         double L = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
         double N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);

         return (Lx*Nx + Ly*Ny + Lz*Nz) / (L*N);
      }
};

vector<PJ_FACTORS> get_factors(const vector<Dimension>& dims,
                               const std::string& CRS,
                               const bool& curvilinear = false) {
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

   PJ_COORD pt, gpt;

   int idx = 0;

   vector<PJ_FACTORS> factors(n);

   for (auto i = 0; i < ni; ++i) {
      for (auto j = 0; j < nj; ++j) {
         idx = ni * j + i;
         pt.enu.e = ioffset + i * idelta;
         pt.enu.n = joffset + j * jdelta;
         gpt = proj_trans(P, PJ_INV, pt);
         factors[idx] = proj_factors(P, gpt);
      }
   }

   return factors;
}

// [[Rcpp::export]]
Rcpp::List rcpp_get_factors(const Rcpp::List& dimensions, const std::string& CRS,
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

std::pair<Rcpp::NumericMatrix, Rcpp::NumericMatrix> get_shifts(const int&  ksize) {

   Rcpp::NumericMatrix ishift(ksize, ksize);
   Rcpp::NumericMatrix jshift(ksize, ksize);

   int istart = ksize / 2;
   int jstart = ksize / 2;

   int di = -istart;
   int dj;

   for (auto k = 0; k < ksize; ++k) {
      dj = -jstart;
      for (auto l = 0; l < ksize; ++l) {
         ishift(k, l) = di;
         jshift(k, l) = dj;
         dj++;
      }
      di++;
   }

   return std::pair(ishift, jshift);
}

std::tuple<Rcpp::NumericMatrix, Rcpp::NumericMatrix, double> get_xy_kernel(const int& i,
                                                                           const int& j,
                                                                           const int& ksize,
                                                                           const vector<Dimension>& dims,
                                                                           const vector<PJ_FACTORS>& pf,
                                                                           const bool& fixed = false,
                                                                           const double& dfactor = 1.0,
                                                                           const Method& = DIRECT) {

   auto idx = j * (dims[0].to - dims[0].from + 1) + i;

   int nk = fixed ? ksize : 2 * ceil(ksize * pf[idx].tissot_semimajor / 2.0) + 1;
   int nl = fixed ? ksize : 2 * ceil(ksize * pf[idx].tissot_semimajor / 2.0) + 1;

   auto [ishift, jshift] = get_shifts(nk);

   double lambdaScale = ksize * pf[idx].meridional_scale / nk;
   double phiScale = ksize * pf[idx].parallel_scale / nl;
   double parallel_convergence = atan2(pf[idx].dy_dlam, pf[idx].dx_dlam);
   auto nrow = ishift.nrow();
   auto ncol = ishift.ncol();

   Rcpp::NumericMatrix x(nrow, ncol);
   Rcpp::NumericMatrix y(nrow, ncol);

   double D, A, a, mu, di, dj;

   for (auto k = 0; k < nrow; ++k) {
      for (auto l = 0; l < ncol; ++l) {
         di = ishift(k, l);
         dj = jshift(k, l);
         if (di == 0 and dj == 0) {
            x(k, l) = 0;
            y(k, l) = 0;
         } else if (dj == 0) {
            x(k, l) = di * dims[0].delta * phiScale * cos(parallel_convergence);
            y(k, l) = di * dims[1].delta * phiScale * sin(parallel_convergence);
         } else if (di == 0) {
            x(k, l) = dj * dims[0].delta * lambdaScale * sin(pf[idx].meridian_convergence);
            y(k, l) = dj * dims[1].delta * lambdaScale * cos(pf[idx].meridian_convergence);
         } else {
            D = sqrt(pow(di * dims[0].delta, 2) + pow(dj * dims[1].delta, 2));
            A = east_to_geo(atan2(dj * dims[1].delta, di * dims[0].delta));
            a = to_closest(north_to_geo(atan(phiScale * sin(pf[idx].meridian_parallel_angle) * tan(A) /
                           (lambdaScale + phiScale * cos(pf[idx].meridian_parallel_angle) * tan(A)))), A);
            mu = sqrt(pow(lambdaScale, 2) * pow(cos(A), 2) +
                      lambdaScale * phiScale * cos(pf[idx].meridian_parallel_angle) * sin(2 * A) +
                      pow(phiScale, 2) * pow(sin(A), 2));

            x(k, l) = D * mu * sin(a - pf[idx].meridian_convergence);
            y(k, l) = D * mu * cos(a - pf[idx].meridian_convergence);
         }
      }
   }

   double scale = 1.0;

   if (fixed) {
      auto abs_compare = [](double a, double b){ return std::abs(a) < std::abs(b); };
      double xmax = std::abs(*std::max_element(x.begin(), x.end(), abs_compare));
      double ymax = std::abs(*std::max_element(y.begin(), y.end(), abs_compare));
      int k2 = ksize / 2;
      scale = k2 * dims[0].delta / std::max(xmax, ymax);
      cout << k2 << ' ' << dims[0].delta << ' ' << std::max(xmax, ymax) << ' ' << scale << endl;
   }

   x = x * scale + i * dims[0].delta + dims[0].offset;
   y = y * scale + j * dims[1].delta + dims[1].offset;

   return std::tuple(x, y, scale);

}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_filter_matrix(const Rcpp::NumericMatrix&  matrix,
                                       const Rcpp::List& dimensions,
                                       const std::string& CRS,
                                       const int& ksize,
                                       const std::vector<std::string>& stats,
                                       const bool& curvilinear = false,
                                       const bool& adaptive = false,
                                       const bool& fixed = false) {
   int ni = matrix.nrow();
   int nj = matrix.ncol();

   int nk = ksize;
   int nl = ksize;

   auto [ishift, jshift] = get_shifts(ksize);

   Rcpp::NumericMatrix nodata(ni, nj);
   for (auto i = 0; i < ni; ++i)
      for (auto j = 0; j < nj; ++j)
         nodata(i, j) = to_string(matrix(i, j)) == "nan";

   Rcpp::NumericMatrix is_bilinear(ni, nj);
   for (auto i = 0; i < ni-1; ++i)
      for (auto j = 0; j < nj-1; ++j)
         is_bilinear(i, j) = nodata(i, j) + nodata(i, j + 1) + nodata(i + 1, j) + nodata(i + 1, j + 1) < 1;
   for (auto i = 0; i < ni-1; ++i)
         is_bilinear(i, nj-1) = nodata(i, nj-1) + nodata(i + 1, nj-1) < 1;
   for (auto j = 0; j < nj-1; ++j)
      is_bilinear(ni-1, j) = nodata(ni-1, j) + nodata(ni-1, j+1) < 1;
   is_bilinear(ni-1, nj-1) = nodata(ni-1, nj-1) < 1;

   Rcpp::NumericMatrix res(ni, nj);
   int ikl, jkl, idi, idj;

   if (adaptive) {
      vector<Dimension> dims;
      for (Rcpp::List dim: dimensions) {
         dims.emplace_back(dim);
      }

      auto factors = get_factors(dims, CRS, curvilinear);

      auto coef = bilinear_coef(matrix);
      double di, dj;

      std::vector<double> values;

      for (auto i = 0; i < ni; ++i) {
         for (auto j = 0; j < nj; ++j) {
            if (nodata(i, j) == 1) {
               res(i, j) = NA_REAL;
            } else {
               auto [x, y, scale] = get_xy_kernel(i, j, ksize, dims, factors, fixed);

               nk = x.nrow();
               nl = x.ncol();

               for (auto k = 0; k < nk; ++k) {
                  for (auto l = 0; l < nl; ++l) {

                     di = (x(k, l) - dims[0].offset) / dims[0].delta;
                     dj = (y(k, l) - dims[1].offset) / dims[1].delta;

                     cout << di << ' ' << dj << endl;

                     if (di >= 0 and dj >= 0 and di < ni and dj < nj) {
                        idi = floor(di);
                        idj = floor(dj);

                        if (is_bilinear(idi, idj) == 1) {
                           values.push_back(interpolate_ij(coef, di, dj));
                        }
                     }
                  }
               }

               cout << endl;

               if (double n = values.size(); n > 0) {
                  for (auto stat : stats) {
                     if (stat == "mean") {
                        res(i, j) = std::accumulate(values.begin(), values.end(), 0) / n;
                     } else if (stat == "slope") {
                        auto surf = ZevenbergenSurface(matrix, dims, coef, x, y, scale);
                        res(i, j) = surf.slope();
                     } else if (stat == "hill") {
                        auto surf = ZevenbergenSurface(matrix, dims, coef, x, y, scale);
                        res(i, j) = surf.hillshade();
                     } else {
                        res(i, j) = std::accumulate(values.begin(), values.end(), 0) / n; // use mean by default
                     }
                  }

               } else {
                  res(i, j) = NA_REAL;
               }

               values.clear();

            }
         }
      }
   } else {
      std::vector<double> values;
      values.reserve(nk * nl);
      for (auto i = 0; i < ni; ++i) {
         for (auto j = 0; j < nj; ++j) {
            if (nodata(i, j) == 1) {
               res(i, j) = NA_REAL;
            } else {
               for (auto k = 0; k < nk; ++k) {
                  for (auto l = 0; l < nl; ++l) {
                     ikl = i + ishift(k, l);
                     jkl = j + jshift(k, l);
                     if (ikl >= 0 and jkl >= 0 and ikl < ni and jkl < nj) {
                        if (nodata(ikl, jkl) == 0) {
                           values.push_back(matrix(ikl, jkl));
                        }
                     }
                  }
               }
               res(i, j) = std::accumulate(values.begin(), values.end(), 0) / (double) (values.size());

               values.clear();
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
