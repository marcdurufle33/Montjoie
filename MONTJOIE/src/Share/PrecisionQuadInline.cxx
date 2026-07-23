#ifndef MONTJOIE_FILE_PRECISION_QUAD_INLINE_CXX

namespace Montjoie
{
  inline double toDouble(const __float128& a)
  {
    return double(a);
  }


  inline int toInteger(const __float128& a)
  {
    return int(a);
  }

}

namespace Seldon
{

#ifdef SELDON_WITH_MPI
  inline MPI_Datatype GetMpiDataType(const Vector<__float128>&)
  {
    return MPI_REAL16;
  }
  
  inline MPI_Datatype GetMpiDataType(const Vector<complex<__float128> >&)
  {
    return MPI_REAL16;
  }

  inline MPI_Datatype GetMpiDataType(const __float128&)
  {
    return MPI_REAL16;
  }
  
  inline MPI_Datatype GetMpiDataType(const complex<__float128>&)
  {
    return MPI_REAL16;
  }
#endif

}

namespace std
{

  /**************
   * Real stuff *
   **************/


  /*inline __float128 abs(const __float128& x)
  {
    return fabsq(x);
    }*/


  inline __float128 floor(const __float128& x)
  {
    return floorq(x);
  }


  inline __float128 ceil(const __float128& x)
  {
    return ceilq(x);
  }


  inline __float128 round(const __float128& x)
  {
    return roundq(x);
  }


  inline __float128 isnan(const __float128& x)
  {
    return x != x;
  }


  inline __float128 isinf(const __float128& x)
  {
    return fabsq(x) > 1.18973149535723176508575932662800702e4932Q;
  }


  inline __float128 sqrt(const __float128& x)
  {
    return sqrtq(x);
  }


  inline __float128 exp(const __float128& x)
  {
    return expq(x);
  }


  inline __float128 log(const __float128& x)
  {
    return logq(x);
  }


  inline __float128 log10(const __float128& x)
  {
    return log10q(x);
  }


  inline __float128 sin(const __float128& x)
  {
    return sinq(x);
  }


  inline __float128 cos(const __float128& x)
  {
    return cosq(x);
  }


  inline __float128 tan(const __float128& x)
  {
    return tanq(x);
  }


  inline __float128 asin(const __float128& x)
  {
    return asinq(x);
  }


  inline __float128 acos(const __float128& x)
  {
    return acosq(x);
  }


  inline __float128 atan(const __float128& x)
  {
    return atanq(x);
  }


  inline __float128 sinh(const __float128& x)
  {
    return sinq(x);
  }


  inline __float128 cosh(const __float128& x)
  {
    return cosq(x);
  }


  inline __float128 tanh(const __float128& x)
  {
    return tanq(x);
  }


  inline __float128 fmod(const __float128& x, const __float128& y)
  {
    return fmodq(x, y);
  }


  inline __float128 pow(const __float128& x, const __float128& y)
  {
    return powq(x, y);
  }


  inline __float128 atan2(const __float128& x, const __float128& y)
  {
    return atan2q(x, y);
  }


  inline __float128 ldexp(const __float128& x, int n)
  {
    return ldexpq(x, n);
  }

  
  inline __float128 frexp(const __float128& x, int* n)
  {
    return frexpq(x, n);
  }


  inline __float128 trunc(const __float128& x)
  {
    return truncq(x);
  }


  inline __float128 tgamma(const __float128& x)
  {
    return tgammaq(x);
  }


  inline __float128 sign(const __float128& a)
  {
    if (a < 0)
      return -__float128(1);
    else
      return __float128(1);
  }

  
  /*****************
   * Complex stuff *
   *****************/

  
  //////////////////////////
  // Comparison operators //

  inline bool operator==(const complex<__float128>& x, int n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return true;
    
    return false;
  }


  inline bool operator==(const complex<__float128>& x, double n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return true;
    
    return false;
  }


  inline bool operator==(const complex<__float128>& x, const __float128& n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return true;
    
    return false;
  }


  inline bool operator==(const complex<__float128>& x, const complex<double>& y)
  {
    if ((real(x) == real(y)) && (imag(x) == imag(y)))
      return true;
    
    return false;
  }


  inline bool operator!=(const complex<__float128>& x, int n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return false;
    
    return true;
  }


  inline bool operator!=(const complex<__float128>& x, double n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return false;
    
    return true;
  }


  inline bool operator!=(const complex<__float128>& x, const __float128& n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return false;
    
    return true;
  }


  inline bool operator!=(const complex<__float128>& x, const complex<double>& y)
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
  inline complex<__float128> operator+(int n, const complex<__float128>& y)
  {
    return complex<__float128>(n+real(y), imag(y));
  }


  // returns x + y
  inline complex<__float128> operator+(const complex<__float128>& y, int n)
  {
    return complex<__float128>(n+real(y), imag(y));
  }


  // returns x + y
  inline complex<__float128> operator+(double x, const complex<__float128>& y)
  {
    return complex<__float128>(x+real(y), imag(y));
  }


  // returns x + y
  inline complex<__float128> operator+(const complex<__float128>& y, double x)
  {
    return complex<__float128>(x+real(y), imag(y));
  }


  // returns x + y
  inline complex<__float128> operator+(const complex<double>& x, const complex<__float128>& y)
  {
    return complex<__float128>(real(x)+real(y), imag(x)+imag(y));
  }


  // returns x + y
  inline complex<__float128> operator+(const complex<__float128>& y, const complex<double>& x)
  {
    return complex<__float128>(real(x)+real(y), imag(x)+imag(y));
  }


  // Operator + //
  ////////////////
  
  
  ////////////////
  // Operator - //

  
  // returns x - y
  inline complex<__float128> operator-(int x, const complex<__float128>& y)
  {
    return complex<__float128>(x-real(y), -imag(y));
  }


  // returns x - y
  inline complex<__float128> operator-(double x, const complex<__float128>& y)
  {
    return complex<__float128>(x-real(y), -imag(y));
  }
  

  // returns x - y
  inline complex<__float128> operator-(const complex<double>& x, const complex<__float128>& y)
  {
    return complex<__float128>(real(x)-real(y), imag(x) - imag(y));
  }


  // returns y - x
  inline complex<__float128> operator-(const complex<__float128>& y, int x)
  {
    return complex<__float128>(real(y)-x, imag(y));
  }


  // returns y - x
  inline complex<__float128> operator-(const complex<__float128>& y, double x)
  {
    return complex<__float128>(real(y)-x, imag(y));
  }


  // returns y - x
  inline complex<__float128> operator-(const complex<__float128>& y, const complex<double>& x)
  {
    return complex<__float128>(real(y)-real(x), imag(y) - imag(x));
  }


  // Operator - //
  /////////////////


  /////////////////
  // Operator * //
  
  
  // returns x*y
  inline complex<__float128> operator*(double x, const complex<__float128>& y)
  {
    return complex<__float128>(x*real(y), x*imag(y));
  }
  
  
  // returns x*y
  inline complex<__float128> operator*(const complex<__float128>& y, double x)
  {
    return complex<__float128>(x*real(y), x*imag(y));
  }

  
  // returns x*y
  inline complex<__float128> operator*(int x, const complex<__float128>& y)
  {
    return complex<__float128>(x*real(y), x*imag(y));
  }

  
  // returns x*y
  inline complex<__float128> operator*(const complex<__float128>& y, int x)
  {
    return complex<__float128>(x*real(y), x*imag(y));
  }


  // returns x*y
  inline complex<__float128> operator*(const complex<double>& x, const complex<__float128>& y)
  {
    return complex<__float128>(real(x)*real(y)-imag(x)*imag(y),
                                 imag(x)*real(y) + real(x)*imag(y));
  }


  // returns x*y
  inline complex<__float128> operator*(const complex<__float128>& y, const complex<double>& x)
  {
    return complex<__float128>(real(x)*real(y)-imag(x)*imag(y),
                                 imag(x)*real(y) + real(x)*imag(y));
  }


  // Operator * //
  /////////////////
  
  
  /////////////////
  // Operator /  //
  

  // returns x/y
  inline complex<__float128> operator/(double x, const complex<__float128>& y)
  {
    __float128 x_div_module = __float128(x) / (square(real(y)) + square(imag(y)));
    return complex<__float128>(real(y)*x_div_module, -imag(y)*x_div_module);
  }
  
  
  // returns y/x
  inline complex<__float128> operator/(const complex<__float128>& y, double x)
  {
    __float128 coef = __float128(1)/__float128(x);
    return complex<__float128>(real(y)*coef, imag(y)*coef);
  }

  
  // returns x/y
  inline complex<__float128> operator/(int x, const complex<__float128>& y)
  {
    __float128 x_div_module = __float128(x) / (square(real(y)) + square(imag(y)));
    return complex<__float128>(real(y)*x_div_module, -imag(y)*x_div_module);
  }

  
  // returns y/x
  inline complex<__float128> operator/(const complex<__float128>& y, int x)
  {
    __float128 coef = __float128(1)/__float128(x);
    return complex<__float128>(real(y)*coef, imag(y)*coef);
  }


  // returns x/y
  inline complex<__float128> operator/(const complex<double>& x, const complex<__float128>& y)
  {
    return complex<__float128>(real(x), imag(x)) / y;
  }


  // returns y/x
  inline complex<__float128> operator/(const complex<__float128>& y, const complex<double>& x)
  {
    return y / complex<__float128>(real(x), imag(x));
  }


  // returns z / x
/*inline complex<__float128> operator/(const complex<__float128>& z, const __float128& x)
  {
    return complex<__float128>(real(z)/x, imag(z)/x);
  }
  
  
  // returns x / z
  inline complex<__float128> operator/(const __float128& x, const complex<__float128>& z)
  {
    __float128 x_div_module = x / (square(real(z)) + square(imag(z)));
    return complex<__float128>(real(z)*x_div_module, -imag(z)*x_div_module);
  }
*/
  
  // Operator /  //
  /////////////////


  /*********************************
   * Functions for complex numbers *
   *********************************/

  
  //! returns true if x is not a number
  inline bool isnan(const complex<__float128>& x)
  {
    if (isnan(real(x)))
      return true;
    
    if (isnan(imag(x)))
      return true;
    
    return false;
  }
  

  //! returns modulus of x
  inline __float128 abs(const complex<__float128>& x)
  {
    return cabsq(reinterpret_cast<const __complex128&>(x));
    //return sqrt(square(real(x)) + square(imag(x)));
  }
  
  
  //! returns modulus of x  
  inline __float128 ComplexAbs(const complex<__float128>& x)
  {
    return cabsq(reinterpret_cast<const __complex128&>(x));
    //return sqrt(square(real(x)) + square(imag(x)));
  }
  

  //! returns square root of x
  inline complex<__float128> sqrt(const complex<__float128>& x)
  {
    /*__float128 a(real(x)), b(imag(x));
    __float128 mod_x = sqrt(a*a + b*b);
    if (b < 0)
      return complex<__float128>(sqrt((a+mod_x)/2), -sqrt((mod_x-a)/2));
    else
    return complex<__float128>(sqrt((a+mod_x)/2), sqrt((mod_x-a)/2));*/
    
    __complex128 z = csqrtq(reinterpret_cast<const __complex128&>(x));
    return complex<__float128>(crealq(z), cimagq(z));
  }

  inline __float128 Sqrt(const __float128& x)
  {
    return sqrt(x);
  }

  inline complex<__float128> Sqrt(const complex<__float128>& x)
  {
    return sqrt(x);
  }


  //! returns power of z
  inline __float128 pow(const __float128& z, int n)
  {
    return pow(z, __float128(n));
  }
  
  
  //! returns power of z
  inline complex<__float128> pow(const complex<__float128>& z, int n)
  {
    if (n == 0)
      return complex<__float128>(1, 0);
    else if (n < 0)
      return complex<__float128>(1, 0)/pow(z, -n);
    else
      {
        // using powers of 2
        complex<__float128> pow_z(z), res(1, 0);
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
  inline complex<__float128> pow(const complex<__float128>& z, const complex<double>& x)
  {
    complex<__float128> xc(real(x), imag(x));
    __complex128 y = cpowq(reinterpret_cast<const __complex128&>(z),
                           reinterpret_cast<const __complex128&>(xc));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  //! returns power of z
  inline complex<__float128> pow(const complex<__float128>& z, double x)
  {
    __complex128 y = cpowq(reinterpret_cast<const __complex128&>(z), __float128(x));
    return complex<__float128>(crealq(y), cimagq(y));
  }  


  //! returns power of z
  inline complex<__float128> pow(const complex<__float128>& z, const __float128& x)
  {
    __complex128 y = cpowq(reinterpret_cast<const __complex128&>(z), x);
    return complex<__float128>(crealq(y), cimagq(y));
  }


  //! returns power of z
  inline complex<__float128> pow(const complex<__float128>& z, const complex<__float128>& x)
  {
    __complex128 y = cpowq(reinterpret_cast<const __complex128&>(z),
                           reinterpret_cast<const __complex128&>(x));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline __float128 arg(const complex<__float128>& z)
  {
    return cargq(reinterpret_cast<const __complex128&>(z));
  }
  

  inline complex<__float128> cos(const complex<__float128>& z)
  {
    __complex128 y = ccosq(reinterpret_cast<const __complex128&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> sin(const complex<__float128>& z)
  {
    __complex128 y = csinq(reinterpret_cast<const __complex128&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> tan(const complex<__float128>& z)
  {
    __complex128 y = ctanq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> acos(const complex<__float128>& z)
  {
    __complex128 y = cacosq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> asin(const complex<__float128>& z)
  {
    __complex128 y = casinq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> atan(const complex<__float128>& z)
  {
    __complex128 y = catanq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }
  

  inline complex<__float128> cosh(const complex<__float128>& z)
  {
    __complex128 y = ccoshq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> sinh(const complex<__float128>& z)
  {
    __complex128 y = csinhq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> tanh(const complex<__float128>& z)
  {
    __complex128 y = ctanhq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> acosh(const complex<__float128>& z)
  {
    __complex128 y = cacoshq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> asinh(const complex<__float128>& z)
  {
    __complex128 y = casinhq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> atanh(const complex<__float128>& z)
  {
    __complex128 y = catanhq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> exp(const complex<__float128>& z)
  {
    __complex128 y = cexpq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> expi(const __float128& x)
  {
    __complex128 y = cexpiq(x);
    return complex<__float128>(crealq(y), cimagq(y));
  }
  

  inline complex<__float128> log(const complex<__float128>& z)
  {
    __complex128 y = clogq(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }


  inline complex<__float128> log10(const complex<__float128>& z)
  {
    __complex128 y = clog10q(reinterpret_cast<__complex128 const&>(z));
    return complex<__float128>(crealq(y), cimagq(y));
  }
  
}
 

#define MONTJOIE_FILE_PRECISION_QUAD_INLINE_CXX
#endif
