#ifndef MONTJOIE_FILE_PRECISION_MPFR_INLINE_CXX

namespace std
{
  
  //! returns real part of x
  inline mpfr::mpreal real(const mpfr::mpreal& x)
  {
    return x;
  }


  /*********************************
   * Operators for complex numbers *
   *********************************/
  
  
  ////////////////
  // Operator = //
  
  
  //! z = double
  inline void SetComplexReal(double x, complex<mpfr::mpreal>& z)
  {
    z = complex<mpfr::mpreal>(x, 0);
  }
  
  //! affectation x = n
  //inline complex<mpfr::mpreal>& operator=(const complex<mpfr::mpreal>& x, int n)
  // {
  //x = complex<mpfr::mpreal>(n, 0);
  //return x;
  //}

  
  // Operator = 
  /////////////////

  
  //////////////////////////
  // Comparison operators //


  inline bool operator==(const complex<mpfr::mpreal>& x, int n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return true;
    
    return false;
  }


  inline bool operator==(const complex<mpfr::mpreal>& x, double n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return true;
    
    return false;
  }


  inline bool operator==(const complex<mpfr::mpreal>& x, const mpfr::mpreal& n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return true;
    
    return false;
  }


  inline bool operator==(const complex<mpfr::mpreal>& x, const complex<double>& y)
  {
    if ((real(x) == real(y)) && (imag(x) == imag(y)))
      return true;
    
    return false;
  }


  inline bool operator!=(const complex<mpfr::mpreal>& x, int n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return false;
    
    return true;
  }


  inline bool operator!=(const complex<mpfr::mpreal>& x, double n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return false;
    
    return true;
  }


  inline bool operator!=(const complex<mpfr::mpreal>& x, const mpfr::mpreal& n)
  {
    if ((real(x) == n) && (imag(x) == 0))
      return false;
    
    return true;
  }


  inline bool operator!=(const complex<mpfr::mpreal>& x, const complex<double>& y)
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
  inline complex<mpfr::mpreal> operator+(int n, const complex<mpfr::mpreal>& y)
  {
    return complex<mpfr::mpreal>(n+real(y), imag(y));
  }


  // returns x + y
  inline complex<mpfr::mpreal> operator+(const complex<mpfr::mpreal>& y, int n)
  {
    return complex<mpfr::mpreal>(n+real(y), imag(y));
  }


  // returns x + y
  inline complex<mpfr::mpreal> operator+(double x, const complex<mpfr::mpreal>& y)
  {
    return complex<mpfr::mpreal>(x+real(y), imag(y));
  }


  // returns x + y
  inline complex<mpfr::mpreal> operator+(const complex<mpfr::mpreal>& y, double x)
  {
    return complex<mpfr::mpreal>(x+real(y), imag(y));
  }


  // returns x + y
  inline complex<mpfr::mpreal> operator+(const complex<double>& x, const complex<mpfr::mpreal>& y)
  {
    return complex<mpfr::mpreal>(real(x)+real(y), imag(x)+imag(y));
  }


  // returns x + y
  inline complex<mpfr::mpreal> operator+(const complex<mpfr::mpreal>& y, const complex<double>& x)
  {
    return complex<mpfr::mpreal>(real(x)+real(y), imag(x)+imag(y));
  }


  // Operator + //
  ////////////////
  
  
  ////////////////
  // Operator - //

  
  // returns x - y
  inline complex<mpfr::mpreal> operator-(int x, const complex<mpfr::mpreal>& y)
  {
    return complex<mpfr::mpreal>(x-real(y), -imag(y));
  }


  // returns x - y
  inline complex<mpfr::mpreal> operator-(double x, const complex<mpfr::mpreal>& y)
  {
    return complex<mpfr::mpreal>(x-real(y), -imag(y));
  }
  

  // returns x - y
  inline complex<mpfr::mpreal> operator-(const complex<double>& x, const complex<mpfr::mpreal>& y)
  {
    return complex<mpfr::mpreal>(real(x)-real(y), imag(x) - imag(y));
  }


  // returns y - x
  inline complex<mpfr::mpreal> operator-(const complex<mpfr::mpreal>& y, int x)
  {
    return complex<mpfr::mpreal>(real(y)-x, imag(y));
  }


  // returns y - x
  inline complex<mpfr::mpreal> operator-(const complex<mpfr::mpreal>& y, double x)
  {
    return complex<mpfr::mpreal>(real(y)-x, imag(y));
  }


  // returns y - x
  inline complex<mpfr::mpreal> operator-(const complex<mpfr::mpreal>& y, const complex<double>& x)
  {
    return complex<mpfr::mpreal>(real(y)-real(x), imag(y) - imag(x));
  }


  // Operator - //
  /////////////////


  /////////////////
  // Operator * //
  
  
  // returns x*y
  inline complex<mpfr::mpreal> operator*(double x, const complex<mpfr::mpreal>& y)
  {
    return complex<mpfr::mpreal>(x*real(y), x*imag(y));
  }
  
  
  // returns x*y
  inline complex<mpfr::mpreal> operator*(const complex<mpfr::mpreal>& y, double x)
  {
    return complex<mpfr::mpreal>(x*real(y), x*imag(y));
  }

  
  // returns x*y
  inline complex<mpfr::mpreal> operator*(int x, const complex<mpfr::mpreal>& y)
  {
    return complex<mpfr::mpreal>(x*real(y), x*imag(y));
  }

  
  // returns x*y
  inline complex<mpfr::mpreal> operator*(const complex<mpfr::mpreal>& y, int x)
  {
    return complex<mpfr::mpreal>(x*real(y), x*imag(y));
  }


  // returns x*y
  inline complex<mpfr::mpreal> operator*(const complex<double>& x, const complex<mpfr::mpreal>& y)
  {
    return complex<mpfr::mpreal>(real(x)*real(y)-imag(x)*imag(y),
                                 imag(x)*real(y) + real(x)*imag(y));
  }


  // returns x*y
  inline complex<mpfr::mpreal> operator*(const complex<mpfr::mpreal>& y, const complex<double>& x)
  {
    return complex<mpfr::mpreal>(real(x)*real(y)-imag(x)*imag(y),
                                 imag(x)*real(y) + real(x)*imag(y));
  }


  // Operator * //
  /////////////////
  
  
  /////////////////
  // Operator /  //
  

  // returns x/y
  inline complex<mpfr::mpreal> operator/(double x, const complex<mpfr::mpreal>& y)
  {
    mpfr::mpreal x_div_module = mpfr::mpreal(x) / (square(real(y)) + square(imag(y)));
    return complex<mpfr::mpreal>(real(y)*x_div_module, -imag(y)*x_div_module);
  }
  
  
  // returns y/x
  inline complex<mpfr::mpreal> operator/(const complex<mpfr::mpreal>& y, double x)
  {
    mpfr::mpreal coef = mpfr::mpreal(1)/mpfr::mpreal(x);
    return complex<mpfr::mpreal>(real(y)*coef, imag(y)*coef);
  }

  
  // returns x/y
  inline complex<mpfr::mpreal> operator/(int x, const complex<mpfr::mpreal>& y)
  {
    mpfr::mpreal x_div_module = mpfr::mpreal(x) / (square(real(y)) + square(imag(y)));
    return complex<mpfr::mpreal>(real(y)*x_div_module, -imag(y)*x_div_module);
  }

  
  // returns y/x
  inline complex<mpfr::mpreal> operator/(const complex<mpfr::mpreal>& y, int x)
  {
    mpfr::mpreal coef = mpfr::mpreal(1)/mpfr::mpreal(x);
    return complex<mpfr::mpreal>(real(y)*coef, imag(y)*coef);
  }


  // returns x/y
  inline complex<mpfr::mpreal> operator/(const complex<double>& x, const complex<mpfr::mpreal>& y)
  {
    return complex<mpfr::mpreal>(real(x), imag(x)) / y;
  }


  // returns y/x
  inline complex<mpfr::mpreal> operator/(const complex<mpfr::mpreal>& y, const complex<double>& x)
  {
    return y / complex<mpfr::mpreal>(real(x), imag(x));
  }


  // returns z / x
  inline complex<mpfr::mpreal> operator/(const complex<mpfr::mpreal>& z, const mpfr::mpreal& x)
  {
    return complex<mpfr::mpreal>(real(z)/x, imag(z)/x);
  }
  
  
  // returns x / z
  inline complex<mpfr::mpreal> operator/(const mpfr::mpreal& x, const complex<mpfr::mpreal>& z)
  {
    mpfr::mpreal x_div_module = x / (square(real(z)) + square(imag(z)));
    return complex<mpfr::mpreal>(real(z)*x_div_module, -imag(z)*x_div_module);
  }

  
  // Operator /  //
  /////////////////


  /*********************************
   * Functions for complex numbers *
   *********************************/

  
  //! returns true if x is not a number
  inline bool isnan(const complex<mpfr::mpreal>& x)
  {
    if (mpfr::isnan(real(x)))
      return true;
    
    if (mpfr::isnan(imag(x)))
      return true;
    
    return false;
  }
  

  //! returns modulus of x
  inline mpfr::mpreal abs(const complex<mpfr::mpreal>& x)
  {
    return sqrt(square(real(x)) + square(imag(x)));
  }
  
  
  //! returns modulus of x  
  inline mpfr::mpreal ComplexAbs(const complex<mpfr::mpreal>& x)
  {
    return sqrt(square(real(x)) + square(imag(x)));
  }
  

  //! returns square root of x
  inline complex<mpfr::mpreal> sqrt(const complex<mpfr::mpreal>& x)
  {
    mpfr::mpreal a(real(x)), b(imag(x));
    mpfr::mpreal mod_x = sqrt(a*a + b*b);
    if (b < 0)
      return complex<mpfr::mpreal>(sqrt((a+mod_x)/2), -sqrt((mod_x-a)/2));
    else
      return complex<mpfr::mpreal>(sqrt((a+mod_x)/2), sqrt((mod_x-a)/2));
  }

  inline mpfr::mpreal Sqrt(const mpfr::mpreal& x)
  {
    return sqrt(x);
  }

  inline complex<mpfr::mpreal> Sqrt(const complex<mpfr::mpreal>& x)
  {
    return sqrt(x);
  }

  
  //! returns power of z
  inline complex<mpfr::mpreal> pow(const complex<mpfr::mpreal>& z, int n)
  {
    if (n == 0)
      return complex<mpfr::mpreal>(1, 0);
    else if (n < 0)
      return complex<mpfr::mpreal>(1, 0)/pow(z, -n);
    else
      {
        // using powers of 2
        complex<mpfr::mpreal> pow_z(z), res(1, 0);
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
  inline complex<mpfr::mpreal> pow(const complex<mpfr::mpreal>& z, double x)
  {
    return pow(z, mpfr::mpreal(x));
  }


  //! returns power of z
  inline complex<mpfr::mpreal> pow(const complex<mpfr::mpreal>& z, complex<double> x)
  {
    return pow(z, complex<mpfr::mpreal>(real(x), imag(x)));
  }

}

namespace Montjoie
{

  /******************************
   * Functions for real numbers *
   ******************************/


  //! conversion to double
  inline double toDouble(const mpreal& a)
  {
    return double(a);
  }
  
  //! conversion to integer
  inline int toInteger(const mpreal& a)
  {
    return int(long(a));
  }

}


namespace std
{
  
  template<class T>
  inline T jn(int n, const T& x)
  {
    return besseljn(n, x);
  }

  template<class T>
  inline T yn(int n, const T& x)
  {
    return besselyn(n, x);
  }
  
  template<class T>
  inline T tgamma(const T& x)
  {
    return gamma(x);
  }


  inline mpfr::mpreal sign(const mpfr::mpreal& a)
  {
    if (a < 0)
      return -mpfr::mpreal(1);
    else
      return mpfr::mpreal(1);
  }
  
}
 

#define MONTJOIE_FILE_PRECISION_MPFR_INLINE_CXX
#endif
