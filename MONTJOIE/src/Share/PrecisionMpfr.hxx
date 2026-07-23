#ifndef MONTJOIE_FILE_PRECISION_MPFR_HXX

#include <gmpxx.h>
#include "Share/mpreal.h"

namespace Montjoie
{
  // Montjoie functions
  using namespace mpfr;

  double toDouble(const mpreal& a);
  int toInteger(const mpreal& a);  

}

namespace std
{
  // standard functions
  istream& operator>>(istream& in, complex<mpfr::mpreal>& z);
  
  template<class T>
  T jn(int n, const T& x);

  template<class T>
  T yn(int n, const T& x);
  
  template<class T>
  T tgamma(const T& x);

  mpfr::mpreal sign(const mpfr::mpreal& a);
  
  mpfr::mpreal real(const mpfr::mpreal& f);
  
  bool isnan(const complex<mpfr::mpreal>& x);
  
  mpfr::mpreal abs(const complex<mpfr::mpreal>& x);
  
  mpfr::mpreal ComplexAbs(const complex<mpfr::mpreal>& x);

  mpfr::mpreal Sqrt(const mpfr::mpreal& x);
  
  complex<mpfr::mpreal> Sqrt(const complex<mpfr::mpreal>& x);
}

#define MONTJOIE_FILE_PRECISION_MPFR_HXX
#endif


