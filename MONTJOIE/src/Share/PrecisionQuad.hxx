#ifndef MONTJOIE_FILE_PRECISION_QUAD_HXX

// including quadmath library
#ifdef MONTJOIE_USE_FLOAT128
#include <quadmath.h>
#endif

namespace Montjoie
{
  double toDouble(const __float128& a);
  int toInteger(const __float128& a);  

}

namespace std
{
  ostream& operator<<(ostream& out, const __float128& x);
  istream& operator>>(istream& in, __float128& x);
  istream& operator>>(istream& in, complex<__float128>& z);
  
  //__float128 abs(const __float128& x);
  __float128 floor(const __float128& x);
  __float128 ceil(const __float128& x);
  __float128 round(const __float128& x);
  __float128 isnan(const __float128& x);
  __float128 isinf(const __float128& x);
  __float128 sqrt(const __float128& x);
  __float128 exp(const __float128& x);
  __float128 log(const __float128& x);
  __float128 log10(const __float128& x);
  __float128 sin(const __float128& x);
  __float128 cos(const __float128& x);
  __float128 tan(const __float128& x);
  __float128 asin(const __float128& x);
  __float128 acos(const __float128& x);
  __float128 atan(const __float128& x);
  __float128 sinh(const __float128& x);
  __float128 cosh(const __float128& x);
  __float128 tanh(const __float128& x);
  __float128 fmod(const __float128& x, const __float128& y);
  __float128 pow(const __float128& x, const __float128& y);
  __float128 atan2(const __float128& x, const __float128& y);
  __float128 ldexp(const __float128& x, int n);
  __float128 frexp(const __float128& x, int* n);
  __float128 trunc(const __float128& x);
  __float128 tgamma(const __float128& x);
  

  bool isnan(const complex<__float128>& x);
  __float128 sign(const __float128& a);

  __float128 abs(const complex<__float128>& x);
  __float128 ComplexAbs(const complex<__float128>& x);

  complex<__float128> sqrt(const complex<__float128>& x);
  __float128 Sqrt(const __float128& x);
  complex<__float128> Sqrt(const complex<__float128>& x);

  __float128 pow(const __float128& x, int n);
  complex<__float128> pow(const complex<__float128>& z, int n);
  complex<__float128> pow(const complex<__float128>& z, const complex<double>& x);
  complex<__float128> pow(const complex<__float128>& z, double x);
  complex<__float128> pow(const complex<__float128>& z, const __float128& x);
  complex<__float128> pow(const complex<__float128>& z, const complex<__float128>& x);

  __float128 arg(const complex<__float128>& z);
  
  complex<__float128> cos(const complex<__float128>& z);
  complex<__float128> sin(const complex<__float128>& z);
  complex<__float128> tan(const complex<__float128>& z);
  complex<__float128> acos(const complex<__float128>& z);
  complex<__float128> asin(const complex<__float128>& z);
  complex<__float128> atan(const complex<__float128>& z);
  
  complex<__float128> cosh(const complex<__float128>& z);
  complex<__float128> sinh(const complex<__float128>& z);
  complex<__float128> tanh(const complex<__float128>& z);
  complex<__float128> acosh(const complex<__float128>& z);
  complex<__float128> asinh(const complex<__float128>& z);
  complex<__float128> atanh(const complex<__float128>& z);

  complex<__float128> exp(const complex<__float128>& x);
  complex<__float128> expi(const __float128& x);

  complex<__float128> log(const complex<__float128>& x);
  complex<__float128> log10(const complex<__float128>& x);

} // end namespace

#define MONTJOIE_FILE_PRECISION_QUAD_HXX
#endif


