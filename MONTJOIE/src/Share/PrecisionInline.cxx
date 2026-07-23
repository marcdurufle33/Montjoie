#ifndef MONTJOIE_FILE_PRECISION_INLINE_CXX

#ifdef MONTJOIE_WITH_MPFR
#include "PrecisionMpfrInline.cxx"
#endif

#ifdef MONTJOIE_USE_FLOAT128
#include "PrecisionQuadInline.cxx"
#endif


namespace std
{
  inline bool isnan(const complex<double>& x)
  {
    if (isnan(real(x)))
      return true;
    
    if (isnan(imag(x)))
      return true;
    
    return false;
  }
    
  inline bool operator==(const complex<double>& x, int n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return true;
    
    return false;
  }

  inline bool operator!=(const complex<double>& x, int n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return false;
    
    return true;
  }

  // returns x + y
  inline complex<double> operator+(int n, const complex<double>& y)
  {
    return complex<double>(n+real(y), imag(y));
  }


  // returns x + y
  inline complex<double> operator+(const complex<double>& y, int n)
  {
    return complex<double>(n+real(y), imag(y));
  }

  // returns x - y
  inline complex<double> operator-(int x, const complex<double>& y)
  {
    return complex<double>(x-real(y), -imag(y));
  }

  // returns y - x
  inline complex<double> operator-(const complex<double>& y, int x)
  {
    return complex<double>(real(y)-x, imag(y));
  }

  // returns x*y
  inline complex<double> operator*(int x, const complex<double>& y)
  {
    return complex<double>(x*real(y), x*imag(y));
  }

  
  // returns x*y
  inline complex<double> operator*(const complex<double>& y, int x)
  {
    return complex<double>(x*real(y), x*imag(y));
  }

  // returns x/y
  inline complex<double> operator/(int x, const complex<double>& y)
  {
    double x_div_module = double(x) / (square(real(y)) + square(imag(y)));
    return complex<double>(real(y)*x_div_module, -imag(y)*x_div_module);
  }

  
  // returns y/x
  inline complex<double> operator/(const complex<double>& y, int x)
  {
    double coef = double(1)/double(x);
    return complex<double>(real(y)*coef, imag(y)*coef);
  }

    /*****************
   * Complex stuff *
   *****************/

  
  //////////////////////////
  // Comparison operators //

  inline bool operator==(const complex<long double>& x, int n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return true;
    
    return false;
  }


  inline bool operator==(const complex<long double>& x, double n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return true;
    
    return false;
  }


  inline bool operator==(const complex<long double>& x, const complex<double>& y)
  {
    if ((real(x) == real(y)) && (imag(x) == imag(y)))
      return true;
    
    return false;
  }


  inline bool operator!=(const complex<long double>& x, int n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return false;
    
    return true;
  }


  inline bool operator!=(const complex<long double>& x, double n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return false;
    
    return true;
  }


  inline bool operator!=(const complex<long double>& x, const complex<double>& y)
  {
    if ((real(x) == real(y)) && (imag(x) == imag(y)))
      return false;
    
    return true;
  }


  // Comparison operators //
  //////////////////////////


  ////////////////
  // Operator + //


  // returns x + y
  inline complex<long double> operator+(int n, const complex<long double>& y)
  {
    return complex<long double>(n+real(y), imag(y));
  }


  // returns x + y
  inline complex<long double> operator+(const complex<long double>& y, int n)
  {
    return complex<long double>(n+real(y), imag(y));
  }


  // returns x + y
  inline complex<long double> operator+(double x, const complex<long double>& y)
  {
    return complex<long double>(x+real(y), imag(y));
  }


  // returns x + y
  inline complex<long double> operator+(const complex<long double>& y, double x)
  {
    return complex<long double>(x+real(y), imag(y));
  }


  // returns x + y
  inline complex<long double> operator+(const complex<double>& x, const complex<long double>& y)
  {
    return complex<long double>(real(x)+real(y), imag(x)+imag(y));
  }


  // returns x + y
  inline complex<long double> operator+(const complex<long double>& y, const complex<double>& x)
  {
    return complex<long double>(real(x)+real(y), imag(x)+imag(y));
  }


  // Operator + //
  ////////////////
  
  
  ////////////////
  // Operator - //

  
  // returns x - y
  inline complex<long double> operator-(int x, const complex<long double>& y)
  {
    return complex<long double>(x-real(y), -imag(y));
  }


  // returns x - y
  inline complex<long double> operator-(double x, const complex<long double>& y)
  {
    return complex<long double>(x-real(y), -imag(y));
  }
  

  // returns x - y
  inline complex<long double> operator-(const complex<double>& x, const complex<long double>& y)
  {
    return complex<long double>(real(x)-real(y), imag(x) - imag(y));
  }


  // returns y - x
  inline complex<long double> operator-(const complex<long double>& y, int x)
  {
    return complex<long double>(real(y)-x, imag(y));
  }


  // returns y - x
  inline complex<long double> operator-(const complex<long double>& y, double x)
  {
    return complex<long double>(real(y)-x, imag(y));
  }


  // returns y - x
  inline complex<long double> operator-(const complex<long double>& y, const complex<double>& x)
  {
    return complex<long double>(real(y)-real(x), imag(y) - imag(x));
  }


  // Operator - //
  /////////////////


  /////////////////
  // Operator * //
  
  
  // returns x*y
  inline complex<long double> operator*(double x, const complex<long double>& y)
  {
    return complex<long double>(x*real(y), x*imag(y));
  }
  
  
  // returns x*y
  inline complex<long double> operator*(const complex<long double>& y, double x)
  {
    return complex<long double>(x*real(y), x*imag(y));
  }

  
  // returns x*y
  inline complex<long double> operator*(int x, const complex<long double>& y)
  {
    return complex<long double>(x*real(y), x*imag(y));
  }

  
  // returns x*y
  inline complex<long double> operator*(const complex<long double>& y, int x)
  {
    return complex<long double>(x*real(y), x*imag(y));
  }


  // returns x*y
  inline complex<long double> operator*(const complex<double>& x, const complex<long double>& y)
  {
    return complex<long double>(real(x)*real(y)-imag(x)*imag(y),
                                 imag(x)*real(y) + real(x)*imag(y));
  }


  // returns x*y
  inline complex<long double> operator*(const complex<long double>& y, const complex<double>& x)
  {
    return complex<long double>(real(x)*real(y)-imag(x)*imag(y),
                                 imag(x)*real(y) + real(x)*imag(y));
  }


  // Operator * //
  /////////////////
  
  
  /////////////////
  // Operator /  //
  

  // returns x/y
  inline complex<long double> operator/(double x, const complex<long double>& y)
  {
    typedef long double long_d;
    long double x_div_module = long_d(x) / (square(real(y)) + square(imag(y)));
    return complex<long double>(real(y)*x_div_module, -imag(y)*x_div_module);
  }
  
  
  // returns y/x
  inline complex<long double> operator/(const complex<long double>& y, double x)
  {
    typedef long double long_d;
    long double coef = long_d(1)/long_d(x);
    return complex<long double>(real(y)*coef, imag(y)*coef);
  }

  
  // returns x/y
  inline complex<long double> operator/(int x, const complex<long double>& y)
  {
    typedef long double long_d;
    long double x_div_module = long_d(x) / (square(real(y)) + square(imag(y)));
    return complex<long double>(real(y)*x_div_module, -imag(y)*x_div_module);
  }

  
  // returns y/x
  inline complex<long double> operator/(const complex<long double>& y, int x)
  {
    return complex<long double>(real(y)/x, imag(y)/x);
  }


  // returns x/y
  inline complex<long double> operator/(const complex<double>& x, const complex<long double>& y)
  {
    return complex<long double>(real(x), imag(x)) / y;
  }


  // returns y/x
  inline complex<long double> operator/(const complex<long double>& y, const complex<double>& x)
  {
    return y / complex<long double>(real(x), imag(x));
  }

  
  // Operator /  //
  /////////////////


  /*********************************
   * Functions for complex numbers *
   *********************************/

  
  //! returns true if x is not a number
  inline bool isnan(const complex<long double>& x)
  {
    if (isnan(real(x)))
      return true;
    
    if (isnan(imag(x)))
      return true;
    
    return false;
  }
  

  //! returns modulus of x  
  inline long double ComplexAbs(const complex<long double>& x)
  {
    return abs(x);
  }
  

  inline long double Sqrt(const long double& x)
  {
    return sqrt(x);
  }

  inline complex<long double> Sqrt(const complex<long double>& x)
  {
    return sqrt(x);
  }

  
  //! returns power of z
  inline complex<long double> pow(const complex<long double>& z, int n)
  {
    if (n == 0)
      return complex<long double>(1, 0);
    else if (n < 0)
      return complex<long double>(1, 0)/pow(z, -n);
    else
      {
        // using powers of 2
        complex<long double> pow_z(z), res(1, 0);
        while (n > 1)
          {
            if (n%2 == 1)
              res *= pow_z;
            
            n = n/2;
            pow_z = pow_z*pow_z;
          }
        
        res *= pow_z;
        return res;
      }
  }

  
  //! returns power of z
  inline complex<long double> pow(const complex<long double>& z, const complex<double>& x)
  {
    return pow(z, complex<long double>(x));
  }


  //! returns power of z
  inline complex<long double> pow(const complex<long double>& z, double x)
  {
    typedef long double long_d;
    return pow(z, long_d(x));
  }  
  
}

namespace Montjoie
{

  //! conversion to integer
  inline int toInteger(const double& a)
  {
    return int(a);
  }
  
  
  //! conversion to integer
  inline int toInteger(const float& a)
  {
    return int(a);
  }
  
  
  //! conversion to double
  inline double toDouble(const double& a)
  {
    return a;
  }


  inline double toDouble(const long double& a)
  {
    return double(a);
  }


  inline int toInteger(const long double& a)
  {
    return int(a);
  }


  //! returns machine epsilon (1e-7 for single precision, 1e-15 for double precision)
  template<class real>
  inline real GetPrecisionMachine(const real& x)
  {
    
    //     estimate unit roundoff in quantities of size x.
    
    real  a,b,c,d,eps;
    
    //     this program should function properly on all systems
    //     satisfying the following two assumptions,
    //        1.  the base used in representing floating point
    //            numbers is not a power of three.
    //        2.  the quantity  a  in statement 10 is represented to 
    //            the accuracy used in floating point variables
    //            that are stored in memory.
    //     the statement number 10 and the go to 10 are intended to
    //     force optimizing compilers to generate code satisfying 
    //     assumption 2.
    //     under these assumptions, it should be true that,
    //            a  is not exactly equal to four-thirds,
    //            b  has a zero for its last bit or digit,
    //            c  is not exactly equal to one,
    //            eps  measures the separation of 1.0 from
    //                 the next larger floating point number.
    //     the developers of eispack would appreciate being informed
    //     about any systems where these assumptions do not hold.
    
    //     this version dated 4/6/83.
    
    a= 4e0; d=3e0;
    a = a/d;
    eps = real(0);
    while (eps==real(0))
      {
	b = a - 1e0;
	c = b + b + b;
	eps = abs(c-1e0);
      }
    
    real dzeps = eps*abs(x);//  DISP(dzeps);
    return dzeps;
    
  }
 
  template<class T0, class T1>
  inline void to_complex(const T0& x, T1& y) { y = T1(x); }
  
  template<class T0, class T1>
  inline void to_complex(const complex<T0>& x, T1& y) { y = realpart(x); }

  template<class T0, class T1>
  inline void to_complex(const T0& x, complex<T1>& y) { y = complex<T1>(T1(x), T1(0)); }

  template<class T0, class T1>
  inline void to_complex(const complex<T0>& x, complex<T1>& y) { y = x; }

  inline void SetImaginaryI(Real_wp& x)
  {
    x = 1.0;
  }
  
  inline void SetImaginaryI(Complex_wp& x)
  {
    x = Iwp;
  }


  //! Initialization of constants (pi, I, and epsilon)
  inline void InitPrecisionConstants()
  {
    pi_wp = acos(Real_wp(-1));
    Iwp = complex<Real_wp>(0, Real_wp(1));
    epsilon_machine = GetPrecisionMachine<Real_wp>(Real_wp(1));
  }  

} // namespace Montjoie


namespace std
{
  template<class T>
  inline T imagpart(const T& x)
  {
    return T(0);
  }

  template<class T>
  inline T imagpart(const complex<T>& x)
  {
    return imag(x);
  }

  inline int sign(const int& a)
  {
    if (a < 0)
      return -1;
    else if (a > 0)
      return 1;

    return 0;
  }

  inline float sign(const float& a)
  {
    if (a < 0)
      return -1.0f;
    else
      return 1.0f;
  }

  inline double sign(const double& a)
  {
    if (a < 0)
      return -1.0;
    else
      return 1.0;
  }

  inline long double sign(const long double& a)
  {
    if (a < 0)
      return -1.0;
    else
      return 1.0;
  }

  inline double Sqrt(double x)
  {
    return sqrt(x);
  }

  inline float Sqrt(float x)
  {
    return sqrt(x);
  }

  inline complex<double> Sqrt(const complex<double>& x)
  {
    return sqrt(x);
  }

  inline complex<float> Sqrt(const complex<float>& x)
  {
    return sqrt(x);
  }

  template<class T>
  inline T abs_real(const T& x)
  {
    return abs(x);
  }

  template<class T>
  inline T abs_real(const complex<T>& x)
  {
    if (real(x) < 0)
      return -real(x);
    
    return real(x);
  }

  template<class T>
  inline T abs_cplx(const T& x)
  {
    return abs(x);
  }


  template<class T>
  inline T abs_cplx(const complex<T>& x)
  {
    return abs(x);
  }


  template<class T>
  inline bool operator <(const complex<T>& x, const complex<T>& y)
  {
    if (real(x) < real(y))
      return true;
    
    if (real(x) > real(y))
      return false;
    
    if (imag(x) < imag(y))
      return true;

    return false;
  }


  template<class T>
  inline bool operator >(const complex<T>& x, const complex<T>& y)
  {
    if (real(x) > real(y))
      return true;

    if (real(x) < real(y))
      return false;
    
    if (imag(x) > imag(y))
      return true;

    return false;
  }
  
}
 

#define MONTJOIE_FILE_PRECISION_INLINE_CXX
#endif
