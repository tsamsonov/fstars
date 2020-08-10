#include <vector>
#include <cmath>
#include <iostream>
#include <proj.h>

const double M_7_PI_4 = 7 * M_PI_4;

using namespace std;

enum SurfaceType {
  EVANS,
  ZEVENBERGEN,
  FLORINSKY,
  UNKNOWN
};

static SurfaceType to_surface_type(const char* type) {
  if (strncmp(type, "ZEVENBERGEN", 11) == 0)
    return ZEVENBERGEN;
  else if (strncmp(type, "EVANS", 5) == 0)
    return EVANS;
  else if (strncmp(type, "FLORINSKY", 9) == 0)
    return FLORINSKY;
  else return UNKNOWN;
}

class Surface {

  double A, B, C, D, E, F, G, H, I;
  SurfaceType type;

public:
  Surface(const vector<double>& z, const double& res, SurfaceType stype = ZEVENBERGEN,
          const double& phi = 0, const double& lam = 0, PJ *P = nullptr, bool warn = false) {
    type = stype;
    switch(type) {
    case EVANS:
      D = (z[0] + z[2] + z[3] + z[5] + z[6] + z[8] - 2 * (z[1] + z[5] + z[7])) / (3 * pow(res, 2));
      E = (z[0] + z[1] + z[2] + z[6] + z[7] + z[8] - 2 * (z[3] + z[4] + z[5])) / (3 * pow(res, 2));
      F = (z[2] + z[6] - z[0] - z[8]) / (4 * pow(res, 2));
      G = (z[2] + z[5] + z[8] - z[0] - z[3] - z[6]) / (6 * res);
      H = (z[0] + z[1] + z[2] - z[6] - z[7] - z[8]) / (6 * res);
      I = z[4];
      break;
    case UNKNOWN:
    case ZEVENBERGEN:
      A = ((z[0] + z[2] + z[6] + z[8])/4.f - (z[1] + z[3] + z[5] + z[7])/2.f + z[4]) / pow(res, 4);
      B = ((z[0] + z[2] - z[6] - z[8])/4.f - (z[1] - z[7])/2.0) / pow(res, 3);
      C = ((-z[0] + z[2] - z[6] + z[8])/4.f + (z[1] - z[5])/2.0) / pow(res, 3);
      D = ((z[3] + z[5])/2.f - z[4]) / pow(res, 2);
      E = ((z[1] + z[7])/2.f - z[4]) / pow(res, 2);
      F = (-z[0] + z[2] + z[6] - z[8]) / (4.0 * pow(res, 2));
      G = (-z[3] + z[5]) / (2.0 * res);
      H = (z[1] - z[7]) / (2.0 * res);
      I = z[4];
      break;
    case FLORINSKY:
      PJ_COORD p1, p2;

      p1.lpz.lam = lam - res;
      p1.lpz.phi = phi - res;
      p2.lpz.lam = lam;
      p2.lpz.phi = phi - res;
      double a = proj_lp_dist(P, p1, p2);
      if (warn)
        cout << p1.lpz.lam << ' ' << p1.lpz.phi << ' ' << p2.lpz.lam << ' ' << p2.lpz.phi << ' -> ' << a << endl;

      p1.lpz.lam = lam - res;
      p1.lpz.phi = phi;
      p2.lpz.lam = lam;
      p2.lpz.phi = phi;
      double b = proj_lp_dist(P, p1, p2);
      if (warn)
        cout << p1.lpz.lam << ' ' << p1.lpz.phi << ' ' << p2.lpz.lam << ' ' << p2.lpz.phi << ' -> ' << b << endl;

      p1.lpz.lam = lam - res;
      p1.lpz.phi = phi + res;
      p2.lpz.lam = lam;
      p2.lpz.phi = phi + res;
      double c = proj_lp_dist(P, p1, p2);
      if (warn)
        cout << p1.lpz.lam << ' ' << p1.lpz.phi << ' ' << p2.lpz.lam << ' ' << p2.lpz.phi << ' -> ' << c << endl;

      p1.lpz.lam = lam - res;
      p1.lpz.phi = phi - res;
      p2.lpz.lam = lam - res;
      p2.lpz.phi = phi;
      double d = proj_lp_dist(P, p1, p2);
      if (warn)
        cout << p1.lpz.lam << ' ' << p1.lpz.phi << ' ' << p2.lpz.lam << ' ' << p2.lpz.phi << ' -> ' << d << endl;

      p1.lpz.lam = lam - res;
      p1.lpz.phi = phi;
      p2.lpz.lam = lam - res;
      p2.lpz.phi = phi + res;
      double e = proj_lp_dist(P, p1, p2);
      if (warn)
        cout << p1.lpz.lam << ' ' << p1.lpz.phi << ' ' << p2.lpz.lam << ' ' << p2.lpz.phi << ' -> ' << e << endl;

      // if (warn)
      //   cout << a << ' ' << b << ' ' << c << ' ' << d << ' ' << e << endl;

      // double a, b, c, d, e;
      // a = b = c = d = e = 1;

      double a2 = a * a;
      double b2 = b * b;
      double c2 = c * c;
      double d2 = d * d;
      double e2 = e * e;
      double a4 = a2 * a2;
      double b4 = b2 * b2;
      double c4 = c2 * c2;

      double AA = (z[2] - z[0]) * a2 * c * d * (d + e) +
        (z[5] - z[3]) * b * (a2 * d2 + c2 * e2) +
        (z[8] - z[6]) * a * c2 * e * (d + e);
      double BB = (a2 * c2 * (d + e) * (d + e) +
                  b2 * (a2 * d2 + c2 * e2)) * 2;

      G = AA / BB;

      double CC =  (z[0] + z[2]) * (d2 * (a4 + b4 + b2 * c2) + c2 * e2 * (a2 - b2));
      double DD =  (z[3] + z[5]) * (d2 * (a4 + c4 + b2 * c2) - e2 * (a4 + c4 + a2 * b2));
      double EE =  (z[6] + z[8]) * (e2 * (b4 + c4 + a2 * b2) - a2 * d2 * (b2 - c2));
      double FF = ((z[1] - z[4]) * (a4 - 2 * b2 * c2) + (3 * z[1] - z[4]) * c4 + (z[1] - 3 * z[4]) * b4) * d2;
      double GG = ((z[4] - z[7]) * (c4 - 2 * a2 * b2) + (3 * z[4] - z[7]) * b4 + (z[4] - 3 * z[7]) * a4) * e2;
      double HH = (z[1] * c2 * e2 * (a2 - b2) + z[7] * a2 * d2* (b2 - c2)) * 2;
      double II = d * e * (d + e) * (a4 + b4 + c4) * 3;

      H = (CC - DD - EE + FF + GG - HH) / II;

      break;
    }
  }

  double slope() {
    return 180 * atan(sqrt(G*G + H*H)) / M_PI;
  }

  double aspect() {
    double a = atan2(G, H);
    if (a > 0){
      if(a <= 0.5 * M_PI){
        a = (0.5 * M_PI - a);
      } else {
        a = (2.5 * M_PI - a);
      }
    } else {
      a = (0.5 * M_PI - a);
    }
    return a;
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

  double profc() {
    return -2 * (D*G*G + E*H*H + F*G*H) / (G*G + H*H);
  }

  double planc() {
    return 2 * (D*H*H + E*G*G - F*G*H) / (G*G + H*H);
  }
};
