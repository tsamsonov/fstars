#include <vector>
#include <cmath>

const double M_7_PI_4 = 7 * M_PI_4;

using namespace std;

class ZevenbergenSurface {
  double A, B, C, D, E, F, G, H, I;

public:
  ZevenbergenSurface(const vector<double>& z, const double& res) {
    A = ((z[0] + z[2] + z[6] + z[8])/4.f - (z[1] + z[3] + z[5] + z[7])/2.f + z[4]) / pow(res, 4);
    B = ((z[0] + z[2] - z[6] - z[8])/4.f - (z[1] - z[7])/2.0) / pow(res, 3);
    C = ((-z[0] + z[2] - z[6] + z[8])/4.f + (z[1] - z[5])/2.0) / pow(res, 3);
    D = ((z[3] + z[5])/2.f - z[4]) / pow(res, 2);
    E = ((z[1] + z[7])/2.f - z[4]) / pow(res, 2);
    F = (-z[0] + z[2] + z[6] - z[8]) / (4.0 * pow(res, 2));
    G = (-z[3] + z[5]) / (2.0 * res);
    H = (z[1] - z[7]) / (2.0 * res);
    I = z[4];
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
