#include "zbessel.h"
#include "zbessel.hh"

// Define visibility macro for different compilers
#ifdef _MSC_VER
  #define EXPORT_SYMBOL  // Empty for MSVC, use header declarations
#elif defined(__GNUC__) || defined(__clang__)
  #define EXPORT_SYMBOL [[gnu::visibility("default")]]
#else
  #define EXPORT_SYMBOL
#endif

extern "C" {

EXPORT_SYMBOL
int zbesh(double zr, double zi, double fnu, int kode, int m,
          int n, double *cyr, double *cyi, int *nz) {
  return zbessel::zbesh(zr, zi, fnu, kode, m, n, cyr, cyi, nz);
}

EXPORT_SYMBOL
int zbesi(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz) {
  return zbessel::zbesi(zr, zi, fnu, kode, n, cyr, cyi, nz);
}

EXPORT_SYMBOL
int zbesj(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz) {
  return zbessel::zbesj(zr, zi, fnu, kode, n, cyr, cyi, nz);
}

EXPORT_SYMBOL
int zbesk(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz) {
  return zbessel::zbesk(zr, zi, fnu, kode, n, cyr, cyi, nz);
}

EXPORT_SYMBOL
int zbesy(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz, double *cwrkr, double *cwrki) {
  return zbessel::zbesy(zr, zi, fnu, kode, n, cyr, cyi, nz, cwrkr, cwrki);
}

EXPORT_SYMBOL
int zairy(double zr, double zi, int id, int kode, double *air, double *aii,
          int *nz) {
  return zbessel::zairy(zr, zi, id, kode, air, aii, nz);
}

EXPORT_SYMBOL
int zbiry(double zr, double zi, int id, int kode, double *bir, double *bii) {
  return zbessel::zbiry(zr, zi, id, kode, bir, bii);
}

}  // extern "C"
