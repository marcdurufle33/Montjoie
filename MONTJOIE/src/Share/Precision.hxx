#ifndef MONTJOIE_FILE_PRECISION_HXX

#include <complex>


#ifdef MONTJOIE_WITH_MPFR
#include "PrecisionMpfr.hxx"
#endif

#ifdef MONTJOIE_WITH_FLOAT128
#define MONTJOIE_USE_FLOAT128
#endif

#ifdef MONTJOIE_USE_FLOAT128
#include "PrecisionQuad.hxx"
#endif

namespace Montjoie
{
  using namespace std;
  
  // definition of real and complex numbers 
  // (simple, double precision or multiple precision)

#ifdef MONTJOIE_WITH_MPFR
  // MPFR -> multiple precision floating library
  // it can be used to obtain very accurate results on small cases
  // be careful ...
  typedef mpreal Real_wp; //!< float
#else

#ifdef MONTJOIE_WITH_FLOAT128
  // using library quadmath provided by gcc
  typedef __float128 Real_wp;
#else

#ifdef MONTJOIE_WITH_FLOAT80
  // using library quadmath provided by gcc
  typedef long double Real_wp;
#else
  //! real number (wp = working precision)
  typedef double Real_wp;
#endif

#endif

#endif
  
  //! complex number at working precision
  typedef complex<Real_wp> Complex_wp;
  
  // Real_wp = reel,Complex_wp = complexe 

  //! pi variable at working precision
  extern Real_wp pi_wp;

  //! complex number I
  extern complex<Real_wp> Iwp;

  //! machine precision (1e-7 in simple, 2e-16 in double)
  extern Real_wp epsilon_machine;
  
  // some conversion functions
  int toInteger(const double& a);
  int toInteger(const float& a);
  double toDouble(const double& a);
  int toInteger(const long double& a);
  double toDouble(const long double& a);
  
  template<class real>
  real GetPrecisionMachine(const real& x);
  
  template<class T0, class T1>
  void to_complex(const T0& x, T1& y);
  
  template<class T0, class T1>
  void to_complex(const complex<T0>& x, T1& y);

  template<class T0, class T1>
  void to_complex(const T0& x, complex<T1>& y);

  template<class T0, class T1>
  void to_complex(const complex<T0>& x, complex<T1>& y);

  void SetImaginaryI(Real_wp&);
  void SetImaginaryI(Complex_wp&);
  
  void InitPrecisionConstants();
}

namespace Seldon
{
  typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;  
}

namespace std
{
  template<class T>
  T imagpart(const T& x);

  template<class T>
  T imagpart(const complex<T>& x);  

  int sign(const int& a);
  float sign(const float& a);
  double sign(const double& a);
  long double sign(const long double& a);

  double Sqrt(double x);
  float Sqrt(float x);
  complex<double> Sqrt(const complex<double>& x);
  complex<float> Sqrt(const complex<float>& x);

  template<class T>
  T abs_real(const T& x);

  template<class T>
  T abs_real(const complex<T>& x);  

  template<class T>
  T abs_cplx(const T& x);

  template<class T>
  T abs_cplx(const complex<T>& x);  

  bool isnan(const complex<double>& x);

  template<class T>
  bool operator <(const complex<T>& x, const complex<T>& y);
  
  template<class T>
  bool operator >(const complex<T>& x, const complex<T>& y);
  
} // end namespace

#define MONTJOIE_FILE_PRECISION_HXX
#endif


