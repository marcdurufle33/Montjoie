#ifndef MONTJOIE_FILE_BESSEL_FUNCTIONS_INTERFACE_CXX


#ifndef MONTJOIE_WITH_MULTIPLE
extern "C"
{
  //! fortran routine for bessel functions
  void cbesj_(const void*, const double*, int*, int*, void*, int*, int*);
  void cbesy_(const void*, const double*, int*, int*, void*, int*, void*, int*);
  void cbesk_(const void*, const double*, int*, int*, void*, int*, int*);
  
  //! fortran routine for Hankel functions
  void cbesh_(void*, double*, int*, int*, int*, void*,int*, int*); 
  void lambert_(const double*, double*);
}
#endif

#ifdef MONTJOIE_WITH_MATHLIB
extern "C"
{
  void wwhitm_sub_(void*, void*, void*);
  void wwhitk_sub_(void*, void*, void*);
}
#endif

#ifdef MONTJOIE_WITH_GSL
#include <gsl/gsl_sf_hyperg.h>
#endif
 
namespace Montjoie
{
  
#ifdef MONTJOIE_WITH_MULTIPLE
  Real_wp GetKn(int n, const Real_wp& x)
  {
    cout << "Bessel functions not available in multiple precision " << endl;
    cout << "Recompile Montjoie in double precision " << endl;
    abort();
    return Real_wp(0);
  }


  //! returns bessel function J_n(x) (first kind)  
  Real_wp GetJn(int n, const Real_wp& x)
  {
#ifdef MONTJOIE_WITH_MPFR
    Real_wp res;
    mpfr_jn(res.mpfr_ptr(), n, x.mpfr_srcptr(), mpreal::get_default_rnd());
    return res;
#else
    cout << "Not implemented" << endl;
    abort();
    return Real_wp(0);
#endif
  }
  

  //! returns bessel function Y_n(x) (second kind)
  Real_wp GetYn(int n, const Real_wp& x)
  {
#ifdef MONTJOIE_WITH_MPFR
    Real_wp res;
    mpfr_yn(res.mpfr_ptr(), n, x.mpfr_srcptr(), mpreal::get_default_rnd());
    return res;
#else
    cout << "Not implemented" << endl;
    abort();
    return Real_wp(0);
#endif

  }
  

  //! computes bessel function J_n(x) and its derivative
  void GetDeriveJn(int n, const Real_wp& x, Real_wp& jn_, Real_wp& djn_)
  {
    Real_wp jn_prev = GetJn(n-1, x);
    jn_ = GetJn(n, x);
    if (abs(x) < epsilon_machine)
      {
        Real_wp jn_next = GetJn(n+1, x);
        djn_ = 0.5*(jn_prev - jn_next);
      }
    else
      djn_ = jn_prev - n/x*jn_;
  }


  //! computes bessel functions J_n(x) for integers and complex argument
  void ComputeBesselFunctions(int begin_order, int n,
                              const Complex_wp& z, VectComplex_wp& Jn)
  {
    if (n <= 0)
      {
        Jn.Clear();
        return;
      }
    
    if (imag(z) == Real_wp(0))
      {
        Real_wp x = real(z), j0_, j1_, jn_next;
        Jn.Reallocate(n);
        if (x == Real_wp(0))
          {
            Jn.Fill(0);
            if ((-begin_order >= 0) && (-begin_order < n))
              SetComplexOne(Jn(-begin_order));
            
            return;
          }
        
        j0_ = GetJn(begin_order, x);
        Jn(0) = Complex_wp(j0_, 0);
        if (n > 1)
          {
            j1_ = GetJn(begin_order+1, x);
            Jn(1) = Complex_wp(j1_, 0);
          }
                
        for (int i = 2; i < n; i++)
          {
            jn_next = (2*(begin_order+i)-2)*j1_/x - j0_;
            Jn(i) = Complex_wp(jn_next, 0);
            j0_ = j1_;
            j1_ = jn_next;
          }
        
        return;
      }
    
    cout << "Bessel functions not available in multiple precision and complex numbers " << endl;
    cout << "Recompile Montjoie in double precision " << endl;
    abort();
  }

  void ComputeBesselFunctions(int begin_order, int n,
                              const Complex_wp& z, VectComplex_wp& Jn, VectComplex_wp& Yn)
  {
    if (n <= 0)
      {
        Jn.Clear();
        Yn.Clear();
        return;
      }
    
    if (imag(z) == Real_wp(0))
      {
        Real_wp x = real(z), j0_, j1_, y0_, y1_, jn_next, yn_next;
        Jn.Reallocate(n); Yn.Reallocate(n);
        if (x == Real_wp(0))
          {
            Jn.Fill(0); Yn.Fill(0);
            if ((-begin_order >= 0) && (-begin_order < n))
              SetComplexOne(Jn(-begin_order));
            
            return;
          }

        j0_ = GetJn(begin_order, x);
        Jn(0) = Complex_wp(j0_, 0); 
        y0_ = GetYn(begin_order, x);
        Yn(0) = Complex_wp(y0_, 0); 
        if (n > 1)
          {
            j1_ = GetJn(begin_order+1, x);
            Jn(1) = Complex_wp(j1_, 0);
            y1_ = GetYn(begin_order+1, x);
            Yn(1) = Complex_wp(y1_, 0);
          }
                
        for (int i = 2; i < n; i++)
          {
            jn_next = (2*(begin_order+i)-2)*j1_/x - j0_;
            yn_next = (2*(begin_order+i)-2)*y1_/x - y0_;
            Jn(i) = Complex_wp(jn_next, 0);
            Yn(i) = Complex_wp(yn_next, 0);
            j0_ = j1_; y0_ = y1_;
            j1_ = jn_next; y1_ = yn_next;
          }
        
        return;
      }
    
    cout << "Bessel functions not available in multiple precision and complex numbers " << endl;
    cout << "Recompile Montjoie in double precision " << endl;
    abort();
  }
  
  void ComputeBesselFunctions(int begin_order, int n, const Real_wp& x, VectReal_wp& Jn)
  {
    if (n <= 0)
      {
        Jn.Clear();
        return;
      }
    
    Real_wp j0_, j1_, jn_next, zero(0);
    Jn.Reallocate(n);
    if (x == zero)
      {
        Jn.Fill(zero);
        if ((-begin_order >= 0) && (-begin_order < n))
          SetComplexOne(Jn(-begin_order));
        
        return;
      }
    
    j0_ = GetJn(begin_order, x);
    Jn(0) = j0_;
    if (n > 1)
      {
        j1_ = GetJn(begin_order+1, x);
        Jn(1) = j1_;
      }
    
    for (int i = 2; i < n; i++)
      {
        jn_next = (2*(begin_order+i)-2)*j1_/x - j0_;
        Jn(i) = jn_next;
        j0_ = j1_;
        j1_ = jn_next;
      }    
  }

  void ComputeDeriveBesselFunctions(int begin_order, int n, const Real_wp& x,
                                    VectReal_wp& Jn, VectReal_wp& dJn)
  {
    if (n <= 1)
      {
        Jn.Clear();
        dJn.Clear();
        return;
      }
    
    ComputeBesselFunctions(begin_order, n, x, Jn);
    
    dJn.Reallocate(n);
    if (x == Real_wp(0))
      {
        Real_wp one(1);
        dJn.Fill(0);
        if ((-begin_order+1 >= 0) && (-begin_order+1 < n))
          dJn(-begin_order+1) = one/2;

        if ((-begin_order-1 >= 0) && (-begin_order-1 < n))
          dJn(-begin_order-1) = -one/2;
        
        return;
      }
    
    dJn(0) = -Jn(1) + begin_order*Jn(0)/x;
    for (int i = 1; i < n-1; i++)
      dJn(i) = 0.5*(Jn(i-1) - Jn(i+1));
    
    dJn(n-1) = Jn(n-2) - (begin_order+n-1)*Jn(n-1)/x;
  }
  
  void ComputeBesselAndHankelFunctions(const Real_wp& begin_order, int n, const Complex_wp& z,
                                       VectComplex_wp& Jn, VectComplex_wp& Hn)
  {
    if (n <= 0)
      {
        Jn.Clear();
        Hn.Clear();
        return;
      }
    
    int b = toInteger(begin_order);
    bool int_order = true;
    if (begin_order != Real_wp(b))
      int_order = false;
    
    if ((imag(z) == Real_wp(0)) && (int_order))
      {
        Real_wp x = real(z), j0_, j1_, y0_, y1_, jn_next, yn_next;
        Jn.Reallocate(n); Hn.Reallocate(n);
        if (x == Real_wp(0))
          {
            Jn.Fill(0); Hn.Fill(0);
            if ((-b >= 0) && (-b < n))
              {
                SetComplexOne(Jn(-b));
                SetComplexOne(Hn(-b));
              }
            
            return;
          }

        j0_ = GetJn(b, x);
        Jn(0) = Complex_wp(j0_, 0); 
        y0_ = GetYn(b, x);
        Hn(0) = Complex_wp(j0_, y0_); 
        if (n > 1)
          {
            j1_ = GetJn(b+1, x);
            Jn(1) = Complex_wp(j1_, 0);
            y1_ = GetYn(b+1, x);
            Hn(1) = Complex_wp(j1_, y1_);
          }
                
        for (int i = 2; i < n; i++)
          {
            jn_next = (2*(b+i)-2)*j1_/x - j0_;
            yn_next = (2*(b+i)-2)*y1_/x - y0_;
            Jn(i) = Complex_wp(jn_next, 0);
            Hn(i) = Complex_wp(jn_next, yn_next);
            j0_ = j1_; y0_ = y1_;
            j1_ = jn_next; y1_ = yn_next;
          }
        
        return;
      }

    cout << "Bessel functions not available in multiple precision " << endl;
    cout << "Recompile Montjoie in double precision " << endl;
    abort();
  }

  void ComputeDeriveBesselAndHankel(const Real_wp& begin_order, int n, const Complex_wp& z,
                                    VectComplex_wp& Jn, VectComplex_wp& Hn,
                                    VectComplex_wp& dJn, VectComplex_wp& dHn)
  {
    if (n <= 1)
      {
        Jn.Clear(); Hn.Clear();
        dJn.Clear(); dHn.Clear();
        return;
      }

    ComputeBesselAndHankelFunctions(begin_order, n, z, Jn, Hn);
    
    dJn.Reallocate(n); dHn.Reallocate(n);
    
    // for the first order, we use the following formula
    // dvp{Z_n}{x} = -Z_{n+1}(x) + \frac{n}{x} Z_n(x)  with Z = J or H
    dJn(0) = -Jn(1) + Jn(0) * Real_wp(begin_order)/z ;
    dHn(0) = -Hn(1) + Hn(0) * Real_wp(begin_order)/z ;
    
    // for the middle orders, we use
    // dvp{Z_n}{x} = 1/2 ( Z_{n-1} - Z_{n+1}) with Z = J or H
    for (int i = 1; i < (n-1); i++)
      {
	dJn(i) = (Jn(i-1) - Jn(i+1))*0.5;
	dHn(i) = (Hn(i-1) - Hn(i+1))*0.5;
      }
    
    // for the last order, we use
    // dvp{Z_n}{x} = Z_{n-1}(x) - \frac{n}{x} Z_n(x)  with Z = J or H
    Real_wp last_order = begin_order + n-1;
    dJn(n-1) = Jn(n-2) - Jn(n-1) * last_order/z;
    dHn(n-1) = Hn(n-2) - Hn(n-1) * last_order/z;    
  }
  
  void ComputeSphericalBessel(int n, const Real_wp& x, VectReal_wp& Jn)
  {
    if (n <= 0)
      {
        Jn.Clear();
        return;
      }

    if (x == 0)
      {
        Jn.Reallocate(n);
        Jn.Fill(0);
        Jn(0) = 1.0;
        
        return;
      }
    
    Jn.Reallocate(n);
    Jn(0) = sin(x)/x;
    
    // forward recurrence, backward recurrence would be more accurate...
    if (n >= 1)
      Jn(1) = (Jn(0) - cos(x))/x;
    
    for (int i = 2; i < n; i++)
      Jn(i) = (2*i-1)*Jn(i-1)/x - Jn(i-2);
  }
    
  void ComputeDeriveSphericalBessel(int n, const Real_wp& x, VectReal_wp& Jn, VectReal_wp& dJn)
  {
    if (n <= 1)
      {
        Jn.Clear();
        dJn.Clear();
        return;
      }

    if (x == 0)
      {
        Jn.Reallocate(n);
        dJn.Reallocate(n);
        Jn.Fill(0);
        dJn.Fill(0);
        Jn(0) = 1.0;
        if (n > 1)
          dJn(1) = 1.0/3;
        
        return;
      }
    
    ComputeSphericalBessel(n, x, Jn);
    
    dJn.Reallocate(n);
    dJn(0) = -Jn(1);
    
    for (int i = 1; i < n-1; i++)
      dJn(i) = (i*Jn(i-1) - (i+1)*Jn(i+1))/(2*i+1);
    
    dJn(n-1) = Jn(n-2) - n*Jn(n-1)/x;
  }

  void ComputeSphericalHankel(int n, const Real_wp& x, VectComplex_wp& Hn)
  {
    if (n <= 0)
      {
        Hn.Clear();
        return;
      }

    if (x == 0)
      {
        // cancelling y^n because of its singularity
        Hn.Reallocate(n);
        Hn.Fill(0);
        Hn(0) = Complex_wp(1, 0);
        
        return;
      }
    
    Hn.Reallocate(n);
    Real_wp jn0, yn0, jn1(0), yn1(0), jn_next, yn_next;
    jn0 = sin(x)/x; yn0 = -cos(x)/x;
    Hn(0) = Complex_wp(jn0, yn0);
    
    // forward recurrence, backward recurrence would be more accurate for Jn
    if (n >= 1)
      {
        jn1 = jn0/x + yn0;
        yn1 = yn0/x - jn0;
        Hn(1) = Complex_wp(jn1, yn1);
      }
    
    for (int i = 2; i < n; i++)
      {
        jn_next = (2*i-1)*jn1/x - jn0;
        yn_next = (2*i-1)*yn1/x - yn0;
        jn0 = jn1; jn1 = jn_next;
        yn0 = yn1; yn1 = yn_next;
        Hn(i) = Complex_wp(jn1, yn1);
      }    
  }
  
  void ComputeDeriveSphericalHankel(int n, const Real_wp& x, VectComplex_wp& Hn,
                                    VectComplex_wp& dHn)
  {
    if (n <= 1)
      {
        Hn.Clear();
        dHn.Clear();
        return;
      }

    if (x == 0)
      {
        Hn.Reallocate(n);
        dHn.Reallocate(n);
        Hn.Fill(0);
        dHn.Fill(0);
        Hn(0) = Complex_wp(1.0, 0.0);
        if (n > 1)
          dHn(1) = Complex_wp(1,0)/3;
        
        return;
      }
    
    ComputeSphericalHankel(n, x, Hn);
    
    dHn.Reallocate(n);
    dHn(0) = -Hn(1);
    
    for (int i = 1; i < n-1; i++)
      dHn(i) = (i*Hn(i-1) - (i+1)*Hn(i+1))/(2*i+1);
    
    dHn(n-1) = Hn(n-2) - n*Hn(n-1)/x;
  }
  
  void ComputeDeriveSphericalBesselHankel(const Real_wp& begin_order, int n, const Complex_wp& z,
                                          VectComplex_wp& Jn, VectComplex_wp& Hn,
                                          VectComplex_wp& dJn, VectComplex_wp& dHn)
  {
    if (n <= 1)
      {
        Jn.Clear(); Hn.Clear();
        dJn.Clear(); Hn.Clear();
      }
    
    if (begin_order == Real_wp(0.5))
      {
        Complex_wp zero, one;
        SetComplexZero(zero); SetComplexOne(one);

        Hn.Reallocate(n); Jn.Reallocate(n);
        dHn.Reallocate(n); dJn.Reallocate(n);
        
        if (z == zero)
          {
            Hn.Fill(zero); Jn.Fill(zero);
            dHn.Fill(zero); dJn.Fill(zero);
            Hn(0) = one; Jn(0) = one;
            if (n > 1)
              {
                dJn(1) = one/Real_wp(3);
                dHn(1) = dJn(1);
              }
            
            return;
          }
        
        Complex_wp jn0, yn0, jn1(zero), yn1(zero), jn_next, yn_next;
        jn0 = sin(z)/z; yn0 = -cos(z)/z;
        Jn(0) = jn0; Hn(0) = jn0 + Iwp*yn0;
        
        // forward recurrence, backward recurrence would be more accurate for Jn
        if (n >= 1)
          {
            jn1 = jn0/z + yn0;
            yn1 = yn0/z - jn0;
            Jn(1) = jn1;
            Hn(1) = jn1 + Iwp*yn1;
          }
        
        for (int i = 2; i < n; i++)
          {
            jn_next = (2*i-1)*jn1/z - jn0;
            yn_next = (2*i-1)*yn1/z - yn0;
            jn0 = jn1; jn1 = jn_next;
            yn0 = yn1; yn1 = yn_next;
            Jn(i) = jn1;
            Hn(i) = jn1 + Iwp*yn1;
          }    
        
        dJn(0) = -Jn(1);
        dHn(0) = -Hn(1);
        
        for (int i = 1; i < n-1; i++)
          {
            dJn(i) = (i*Jn(i-1) - (i+1)*Jn(i+1))/(2*i+1);
            dHn(i) = (i*Hn(i-1) - (i+1)*Hn(i+1))/(2*i+1);
          }
    
        dJn(n-1) = Jn(n-2) - n*Jn(n-1)/z;
        dHn(n-1) = Hn(n-2) - n*Hn(n-1)/z;
        
        return;
      }
    
    cout << "Bessel functions not available in multiple precision " << endl;
    cout << "Recompile Montjoie in double precision " << endl;
    abort();
  }
  
  void ComputeDeriveRiccatiBessel(const Real_wp& begin_order, int n, const Complex_wp& z,
                                  VectComplex_wp& Jn, VectComplex_wp& Hn,
                                  VectComplex_wp& Jn_prime, VectComplex_wp& Hn_prime,
                                  VectComplex_wp& Psi_n, VectComplex_wp& Xi_n,
                                  VectComplex_wp& Psi_prime, VectComplex_wp& Xi_prime)
  {
    ComputeDeriveSphericalBesselHankel(begin_order, n, z,
                                       Jn, Hn, Jn_prime, Hn_prime);
    
    Psi_n.Reallocate(n); Xi_n.Reallocate(n);
    Psi_prime.Reallocate(n); Xi_prime.Reallocate(n);
    for (int i = 0; i < n; i++)
      {
        Psi_n(i) = z*Jn(i); Xi_n(i) = z*Hn(i);
        Psi_prime(i) = z*Jn_prime(i) + Jn(i);
        Xi_prime(i) = z*Hn_prime(i) + Hn(i);
      }
  }
#else

  //! returns modified bessel function K_n(x)
  double GetKn(int n, const double& x)
  {
    complex<double> z = x;
    int kode = 1; int nz = 0, ierr = 0, m=1;
    double begin_order = n;
    if (n < 0)
      begin_order = -begin_order;
    
    complex<double> kn;
    cbesk_(reinterpret_cast<void*>(&z), &begin_order, &kode, &m,
	   reinterpret_cast<void*>(&kn), &nz, &ierr);
    
    return real(kn);
  }
  
  
  //! returns bessel function J_n(x) (first kind)
  double GetJn(int n, const double& x)
  {
    return jn(n, x);
  }


  //! returns bessel function Y_n(x) (second kind)
  double GetYn(int n, const double& x)
  {
    return yn(n, x);
  }


  //! computes bessel function J_n(x) and its derivative
  void GetDeriveJn(int n, const double& x, double& jn_, double& djn_)
  {
    double jn_prev = jn(n-1, x);
    jn_ = jn(n, x);
    if (abs(x) < epsilon_machine)
      {
        double jn_next = jn(n+1, x);
        djn_ = 0.5*(jn_prev - jn_next);
      }
    else
      djn_ = jn_prev - n/x*jn_;
  }
  
  
  //! Computation of \f$ J_p(z) \f$ where p goes from begin_order until begin_order + n
  /*!
    \param[in] begin_order first order to compute
    \param[in] n number of bessel functions to compute
    \param[in] z complex argument of bessel functions
    \param[out] Jn bessel functions evaluated to the asked points
   */
  void ComputeBesselFunctions(int begin_order, int n,
                              const complex<double>& z, ZVect& Jn)
  {
    if (n <= 0)
      {
        Jn.Clear();
        return;
      }
    
    // if not allocated
    Jn.Reallocate(n);
    
    int kode = 1; int nz = 0,ierr = 0; 
    int n0 = begin_order; int n1 = begin_order+n-1;
    double n0_real(n0);
    if ((n0 < 0)&&(n1 <= 0))
      {
	ZVect Jn_(n);
	int j = n0;
	n0 = -n1;
	n1 = -j;
	n0_real = double(n0);
	// we call fortran routines to compute bessel and hankel functions
	cbesj_(reinterpret_cast<const void*>(&z), &n0_real, &kode, &n,
	       reinterpret_cast<void*>(Jn_.GetData()), &nz, &ierr);
	
	for (int i = n0; i <= n1; i++)
	  {
	    if ((n1+n0-i)%2 == 0)
	      Jn(i-n0) = Jn_(n0+n-1-i);
	    else
	      Jn(i-n0) = -Jn_(n0+n-1-i);
	  }
      }
    else if (n0 < 0)
      {
	ZVect Jn_(n);
	n0 = -n0+1; double zero = 0.0;
	// we call fortran routines to compute bessel and hankel functions
	cbesj_(reinterpret_cast<const void*>(&z), &zero, &kode, &n0,
	       reinterpret_cast<void*>(Jn_.GetData()), &nz, &ierr);
	
	for (int i = 0; i < n0; i++)
	  {
	    if ((n0-1-i)%2 == 0)
	      Jn(i) = Jn_(n0-1-i);
	    else
	      Jn(i) = -Jn_(n0-1-i);
	  }
	zero = 1.0;
	// we call fortran routines to compute bessel and hankel functions
	cbesj_(reinterpret_cast<const void*>(&z), &zero, &kode, &n1,
	       reinterpret_cast<void*>(Jn_.GetData()), &nz, &ierr);
	
	for (int i = 0; i < n1; i++)
	  Jn(n0+i) = Jn_(i);
	
      }
    else
      {
	// we call fortran routines to compute bessel and hankel functions
	cbesj_(reinterpret_cast<const void*>(&z), &n0_real, &kode, &n,
	       reinterpret_cast<void*>(Jn.GetData()), &nz, &ierr);
      }
  }
    

  //! Computation of \f$ J_p(z) \f$ where p goes from begin_order until begin_order + n
  /*!
    \param[in] begin_order first order to compute
    \param[in] n number of bessel functions to compute
    \param[in] z complex argument of bessel functions
    \param[out] Jn bessel functions evaluated to the asked points
   */
  void ComputeBesselFunctions(int begin_order, int n,
                              const complex<double>& z, ZVect& Jn, ZVect& Yn)
  {
    if (n <= 0)
      {
        Jn.Clear();
        Yn.Clear();
        return;
      }

    // if not allocated
    Jn.Reallocate(n);
    Yn.Reallocate(n);
    int kode = 1; int nz = 0,ierr = 0; 
    int n0 = begin_order; int n1 = begin_order+n-1;
    double n0_real(n0); Vector<Complex_wp> cwrk(n);
    if ((n0 < 0)&&(n1 <= 0))
      {
	ZVect Jn_(n), Yn_(n);
	int j = n0;
	n0 = -n1;
	n1 = -j;
	n0_real = double(n0);
	// we call fortran routines to compute bessel functions
	cbesj_(reinterpret_cast<const void*>(&z), &n0_real, &kode, &n,
	       reinterpret_cast<void*>(Jn_.GetData()), &nz, &ierr);

	cbesy_(reinterpret_cast<const void*>(&z), &n0_real, &kode, &n,
	       reinterpret_cast<void*>(Yn_.GetData()), &nz,
               reinterpret_cast<void*>(cwrk.GetData()), &ierr);
	
	for (int i = n0; i <= n1; i++)
	  {
	    if ((n1+n0-i)%2 == 0)
	      {
                Jn(i-n0) = Jn_(n0+n-1-i);
                Yn(i-n0) = Yn_(n0+n-1-i);
              }
	    else
	      {
                Jn(i-n0) = -Jn_(n0+n-1-i);
                Yn(i-n0) = -Yn_(n0+n-1-i);
              }
	  }
      }
    else if (n0 < 0)
      {
	ZVect Jn_(n), Yn_(n);
	n0 = -n0+1; double zero = 0.0;
	// we call fortran routines to compute bessel and hankel functions
	cbesj_(reinterpret_cast<const void*>(&z), &zero, &kode, &n0,
	       reinterpret_cast<void*>(Jn_.GetData()), &nz, &ierr);

	cbesy_(reinterpret_cast<const void*>(&z), &zero, &kode, &n0,
	       reinterpret_cast<void*>(Yn_.GetData()), &nz,
               reinterpret_cast<void*>(cwrk.GetData()), &ierr);
	
	for (int i = 0; i < n0; i++)
	  {
	    if ((n0-1-i)%2 == 0)
	      {
                Jn(i) = Jn_(n0-1-i);
                Yn(i) = Yn_(n0-1-i);
              }
	    else
	      {
                Jn(i) = -Jn_(n0-1-i);
                Yn(i) = -Yn_(n0-1-i);
              }
	  }
        
	zero = 1.0;
	// we call fortran routines to compute bessel and hankel functions
	cbesj_(reinterpret_cast<const void*>(&z), &zero, &kode, &n1,
	       reinterpret_cast<void*>(Jn_.GetData()), &nz, &ierr);

	cbesy_(reinterpret_cast<const void*>(&z), &zero, &kode, &n1,
	       reinterpret_cast<void*>(Yn_.GetData()), &nz,
               reinterpret_cast<void*>(cwrk.GetData()), &ierr);
	
	for (int i = 0; i < n1; i++)
	  {
            Jn(n0+i) = Jn_(i);
            Yn(n0+i) = Yn_(i);
          }	
      }
    else
      {
	// we call fortran routines to compute bessel and hankel functions
	cbesj_(reinterpret_cast<const void*>(&z), &n0_real, &kode, &n,
	       reinterpret_cast<void*>(Jn.GetData()), &nz, &ierr);

	cbesy_(reinterpret_cast<const void*>(&z), &n0_real, &kode, &n,
	       reinterpret_cast<void*>(Yn.GetData()), &nz,
               reinterpret_cast<void*>(cwrk.GetData()), &ierr);
      }
  }

  
  //! Computation of \f$ J_p(z) \f$ where p goes from begin_order until begin_order + n
  //! z is real
  /*!
    \param[in] begin_order first order to compute
    \param[in] n number of bessel functions to evaluate
    \param[in] x real argument of bessel functions
    \param[out] Jn bessel functions evaluated to the asked points
   */
  void ComputeBesselFunctions(int begin_order, int n, const double& x, DVect& Jn)
  {
    complex<double> z = x; ZVect Jn_cplx;
    ComputeBesselFunctions(begin_order, n, z, Jn_cplx);
    // if not allocated
    Jn.Reallocate(n);
    for (int i = 0; i < n; i++)
      Jn(i) = real(Jn_cplx(i));
    
    return;
  }
  

  //! Computation of \f$ J_p(z) \f$ where p goes from begin_order until begin_order + n
  //! z is real
  /*!
    \param[in] begin_order first order to compute
    \param[in] n number of bessel functions to evaluate
    \param[in] x real argument of bessel functions
    \param[out] Jn bessel functions evaluated to the asked points
   */
  void ComputeDeriveBesselFunctions(int begin_order, int n, const double& x, DVect& Jn, DVect& dJn)
  {
    if (n <= 1)
      {
        Jn.Clear();
        dJn.Clear();
        return;
      }
    
    if (x == 0)
      {
        Jn.Reallocate(n);
        dJn.Reallocate(n);
        Jn.Fill(0);
        dJn.Fill(0);
        if ((begin_order <= 0) && (begin_order >= -n+1))
          Jn(-begin_order) = 1.0;
        
        if ((begin_order <= 1) && (begin_order >= -n+2))
          dJn(-begin_order+1) = 0.5;
        
        return;
      }
    
    complex<double> z = x; ZVect Jn_cplx;
    ComputeBesselFunctions(begin_order, n, z, Jn_cplx);
    // if not allocated
    Jn.Reallocate(n);
    dJn.Reallocate(n);
    for (int i = 0; i < n; i++)
      Jn(i) = real(Jn_cplx(i));
    
    dJn(0) = -Jn(1) + begin_order*Jn(0)/x;
    for (int i = 1; i < n-1; i++)
      dJn(i) = 0.5*(Jn(i-1) - Jn(i+1));
    
    if (n > 1)
      dJn(n-1) = Jn(n-2) - (begin_order+n-1)*Jn(n-1)/x;
      
    return;
  }

  
  //! Computation of \f$ J_p(z) \mbox{ and } H_p^{(1)}(z) \f$
  //! where p goes from begin_order until begin_order + n
  /*!
    \param[in] begin_order first order to compute
    \param[in] n number of bessel functions to evaluate
    \param[in] z complex argument of bessel functions
    \param[out] Jn Bessel functions evaluated to the asked points
    \param[out] Hn Hankel functions evaluated to the asked points
   */
  void ComputeBesselAndHankelFunctions(const double& begin_order, int n, const complex<double>& z,
                                       ZVect& Jn, ZVect& Hn)
  {
    if (n <= 0)
      {
        Jn.Clear();
        Hn.Clear();
        return;
      }

    // if not allocated
    Jn.Reallocate(n); Hn.Reallocate(n);
    
    int kode = 1; int nz = 0, ierr = 0;
    // we call fortran routines to compute bessel and hankel functions
    cbesj_(reinterpret_cast<const void*>(&z), &begin_order, &kode, &n,
	   reinterpret_cast<void*>(Jn.GetData()), &nz, &ierr);
    
    ZVect Yn_(n);
    cbesy_(reinterpret_cast<const void*>(&z), &begin_order, &kode, &n,
	   reinterpret_cast<void*>(Yn_.GetData()), &nz,
	   reinterpret_cast<void*>(Hn.GetData()), &ierr);
    
    for (int i = 0; i < n; i++)
      {
	if (imag(z) == 0)
	  {
	    Jn(i) = real(Jn(i));
	    Yn_(i) = real(Yn_(i));
	  }
	
	Hn(i) = Jn(i) + Iwp*Yn_(i);
      }
  }
  
  
  //! Computation of J_p(z), H_p(z), J'_p(z) and H_p^{'(1)}(z)
  void ComputeDeriveBesselAndHankel(const double& begin_order, int n, const complex<double>& z,
                                    ZVect& Jn, ZVect& Hn,
                                    ZVect& Jn_prime, ZVect& Hn_prime)
  {
    if (n <= 1)
      {
        Jn.Clear(); Hn.Clear();
        Jn_prime.Clear(); Hn_prime.Clear();
        return;
      }
    
    ComputeBesselAndHankelFunctions(begin_order, n, z, Jn, Hn);
    
    // if not allocated
    Jn_prime.Reallocate(n); Hn_prime.Reallocate(n);
    
    // for the first order, we use the following formula
    // dvp{Z_n}{x} = -Z_{n+1}(x) + \frac{n}{x} Z_n(x)  with Z = J or H
    Jn_prime(0) = -Jn(1) + Jn(0) * double(begin_order)/z ;
    Hn_prime(0) = -Hn(1) + Hn(0) * double(begin_order)/z ;
    
    // for the middle orders, we use
    // dvp{Z_n}{x} = 1/2 ( Z_{n-1} - Z_{n+1}) with Z = J or H
    for (int i = 1; i < (n-1); i++)
      {
	Jn_prime(i) = (Jn(i-1) - Jn(i+1))*0.5;
	Hn_prime(i) = (Hn(i-1) - Hn(i+1))*0.5;
      }
    
    // for the last order, we use
    // dvp{Z_n}{x} = Z_{n-1}(x) - \frac{n}{x} Z_n(x)  with Z = J or H
    double last_order = begin_order + n-1;
    Jn_prime(n-1) = Jn(n-2) - Jn(n-1) * last_order/z;
    Hn_prime(n-1) = Hn(n-2) - Hn(n-1) * last_order/z;    
  }
  
  
  //! Compute spherical bessel function j_n(x)
  void ComputeSphericalBessel(int n, const double& x, Vector<double>& Jn)
  {
    if (n <= 0)
      {
        Jn.Clear();
        return;
      }

    if (x == 0)
      {
        Jn.Reallocate(n);
        Jn.Fill(0);
        Jn(0) = 1.0;
        
        return;
      }
    
    complex<double> z(x, 0);
    
    Vector<complex<double> > Jn_(n);
    int kode = 1; int nz = 0, ierr = 0;
    double begin_order = 0.5;
    cbesj_(reinterpret_cast<void*>(&z), &begin_order, &kode, &n,
	   reinterpret_cast<void*>(Jn_.GetData()), &nz, &ierr);
    
    double coef = sqrt(0.5*M_PI/x);
    Jn.Reallocate(n);
    for (int i = 0; i < n; i++)
      Jn(i) = coef * real(Jn_(i));    
  }
  

  //! Compute spherical bessel function j_n(x) and derivatives j_n'(x)
  void ComputeDeriveSphericalBessel(int n, const double& x,
                                    Vector<double>& Jn, Vector<double>& dJn)
  {
    if (n <= 1)
      {
        Jn.Clear();
        dJn.Clear();
        return;
      }

    if (x == 0)
      {
        Jn.Reallocate(n);
        dJn.Reallocate(n);
        Jn.Fill(0);
        dJn.Fill(0);
        Jn(0) = 1.0;
        if (n > 1)
          dJn(1) = 1.0/3;
        
        return;
      }

    complex<double> z(x, 0);
    
    Vector<complex<double> > Jn_(n);
    int kode = 1; int nz = 0, ierr = 0;
    double begin_order = 0.5;
    cbesj_(reinterpret_cast<void*>(&z), &begin_order, &kode, &n,
	   reinterpret_cast<void*>(Jn_.GetData()), &nz, &ierr);
    
    double half_divX = 0.5/x;
    double coef = sqrt(M_PI*half_divX);
    double last_order_divX = half_divX*(2*n-1);
    double jn_prime;
    Jn.Reallocate(n);
    dJn.Reallocate(n);
    for (int i = 0; i < n; i++)
      {
        if (i == 0)
          jn_prime = -real(Jn_(1)) + real(Jn_(0))*half_divX;
        else if (i == n-1)
          jn_prime = real(Jn_(n-2)) - real(Jn_(n-1))*last_order_divX;
        else
          jn_prime = 0.5*(real(Jn_(i-1)) - real(Jn_(i+1)));
        
        Jn(i) = coef * real(Jn_(i));    
        dJn(i) = -Jn(i)*half_divX + coef*jn_prime;
      }
  }
  

  //! Compute spherical Hankel function h_n(x)
  void ComputeSphericalHankel(int n, const double& x, Vector<complex<double> >& Hn)
  {
    if (n <= 0)
      {
        Hn.Clear();
        return;
      }

    if (x == 0)
      {
        // cancelling y^n because of its singularity
        Hn.Reallocate(n);
        Hn.Fill(0);
        Hn(0) = 1.0;
        
        return;
      }
    
    complex<double> z(x, 0);
    
    Vector<complex<double> > Jn_(n), Yn_(n), cwrk(n);
    int kode = 1; int nz = 0, ierr = 0;
    double begin_order = 0.5;
    cbesj_(reinterpret_cast<void*>(&z), &begin_order, &kode, &n,
	   reinterpret_cast<void*>(Jn_.GetData()), &nz, &ierr);
    
    cbesy_(reinterpret_cast<void*>(&z), &begin_order, &kode, &n,
	   reinterpret_cast<void*>(Yn_.GetData()), &nz,
           reinterpret_cast<void*>(cwrk.GetData()), &ierr);
    
    double coef = sqrt(0.5*M_PI/x);
    Hn.Reallocate(n);
    for (int i = 0; i < n; i++)
      Hn(i) = coef * complex<double>(real(Jn_(i)), real(Yn_(i)));    
  }
  

  //! Compute spherical bessel function h_n(x) and derivatives h_n'(x)
  void ComputeDeriveSphericalHankel(int n, const double& x, Vector<complex<double> >& Hn,
                                    Vector<complex<double> >& dHn)
  {
    if (n <= 1)
      {
        Hn.Clear();
        dHn.Clear();
        return;
      }

    if (x == 0)
      {
        Hn.Reallocate(n);
        dHn.Reallocate(n);
        Hn.Fill(0);
        dHn.Fill(0);
        Hn(0) = 1.0;
        if (n > 1)
          dHn(1) = 1.0/3;
        
        return;
      }

    complex<double> z(x, 0);
    
    Vector<complex<double> > Jn_(n), Yn_(n), cwrk(n);
    int kode = 1; int nz = 0, ierr = 0;
    double begin_order = 0.5;
    cbesj_(reinterpret_cast<void*>(&z), &begin_order, &kode, &n,
	   reinterpret_cast<void*>(Jn_.GetData()), &nz, &ierr);

    cbesy_(reinterpret_cast<void*>(&z), &begin_order, &kode, &n,
	   reinterpret_cast<void*>(Yn_.GetData()), &nz,
           reinterpret_cast<void*>(cwrk.GetData()), &ierr);
    
    double half_divX = 0.5/x;
    double coef = sqrt(M_PI*half_divX);
    double last_order_divX = half_divX*(2*n-1);
    double jn_prime, yn_prime;
    Hn.Reallocate(n);
    dHn.Reallocate(n);
    for (int i = 0; i < n; i++)
      {
        if (i == 0)
          {
            jn_prime = -real(Jn_(1)) + real(Jn_(0))*half_divX;
            yn_prime = -real(Yn_(1)) + real(Yn_(0))*half_divX;
          }
        else if (i == n-1)
          {
            jn_prime = real(Jn_(n-2)) - real(Jn_(n-1))*last_order_divX;
            yn_prime = real(Yn_(n-2)) - real(Yn_(n-1))*last_order_divX;
          }
        else
          {
            jn_prime = 0.5*(real(Jn_(i-1)) - real(Jn_(i+1)));
            yn_prime = 0.5*(real(Yn_(i-1)) - real(Yn_(i+1)));
          }
        
        Hn(i) = coef * complex<double>(real(Jn_(i)), real(Yn_(i)));    
        dHn(i) = -Hn(i)*half_divX + coef*complex<double>(jn_prime, yn_prime);
      }
  }

  
  //! Compute spherical bessel function j_n(z) and h_n^{(1)}(z) and derivatives
  void ComputeDeriveSphericalBesselHankel(const double& begin_order, int n,
                                          const complex<double>& z,
                                          ZVect& Jn,
                                          ZVect& Hn,
                                          ZVect& Jn_prime,
                                          ZVect& Hn_prime)
  {
    ComputeDeriveBesselAndHankel(begin_order, n, z, Jn, Hn, Jn_prime, Hn_prime);
    
    // we replace hankel function of the first kind by
    // spherical hankel function of the first kind
    // h_n^(1) (z) = sqrt(pi/(2z)) H_{n+1/2}^{(1)} (z) 
    // and the same thing for bessel function
    // j_n (z) = sqrt(pi/(2z)) J_{n+1/2} (z)
    
    // and we replace H_{n+1/2}^{(1)'}(z), J_{n+1/2}' (z)
    // by derivatives of the spherical hankel and bessel functions
    for (int i = 0; i < n; i++)
      {
	Jn(i) = sqrt(0.5*M_PI/z) * Jn(i);
	Hn(i) = sqrt(0.5*M_PI/z) * Hn(i);
	Jn_prime(i) = -Jn(i)*0.5/z + sqrt(0.5*M_PI/z) * Jn_prime(i);
	Hn_prime(i) = -Hn(i)*0.5/z + sqrt(0.5*M_PI/z) * Hn_prime(i);
      }
  }
  
  
  //! Compute spherical bessel function j_n(z) and h_n^{(1)}(z) and derivatives, and associated Riccati functions
  void ComputeDeriveRiccatiBessel(const double& begin_order, int n,
                                  const complex<double>& z,
                                  ZVect& Jn, ZVect& Hn,
                                  ZVect& Jn_prime, ZVect& Hn_prime,
                                  ZVect& Psi_n, ZVect& Xi_n,
                                  ZVect& Psi_prime,ZVect& Xi_prime)
  {
    ComputeDeriveBesselAndHankel(begin_order, n, z, Jn, Hn, Jn_prime, Hn_prime);
    
    // if not allocated
    Xi_n.Reallocate(n); Psi_n.Reallocate(n);
    Xi_prime.Reallocate(n); Psi_prime.Reallocate(n);
    for (int i = 0;  i < n; i++)
      {
	Jn(i) = sqrt(0.5*M_PI/z) * Jn(i);
	Psi_n(i) = z * Jn(i);
	Hn(i) = sqrt(0.5*M_PI/z) * Hn(i);
	Xi_n(i) = z * Hn(i);
	Psi_prime(i) = Jn(i)*0.5 + sqrt(M_PI*0.5*z) * Jn_prime(i);
	Xi_prime(i) = Hn(i)*0.5 + sqrt(M_PI*0.5*z) * Hn_prime(i);
	Jn_prime(i) = -Jn(i)*0.5/z + sqrt(0.5*M_PI/z) * Jn_prime(i);
	Hn_prime(i) = -Hn(i)*0.5/z + sqrt(0.5*M_PI/z) * Hn_prime(i);
      }
  }
#endif
  
  
  //! evaluation of the maximum mode to evaluate in Bessel expansion
  /*!
    \param[in] ka wave number multiplied by the radius of the circle/sphere
    \param[in] epsilon_machine precision wished
    \returns Nmax of Bessel expansion
  */
  template<class T>
  int ComputeOrder(const T& ka, const T& epsilon_machine)
  {
    T epsilon = epsilon_machine*pow(abs(ka), T(1.0/0.3));
    T v = abs(ka);
    T x_maxi = 1.0, x_mini = 0.0, x = 1.0, b, f, a, g, crit;
    while (abs(x_mini-x_maxi) >= epsilon_machine )
      {
	x = (x_mini + x_maxi)/2;
	b = sqrt(1-x*x);
	f = log(1.0+b) - log(x) -b;
	a = 0.5*sqrt(x)/sqrt(b*b*b);
	g = a*exp(-v*f/x);
	
	crit = g-epsilon;
	if (crit<=0)
	  x_mini = x;
	else
	  x_maxi = x;
      }
    return toInteger(floor(v/x-0.5));
  }
  
  
  //! computation of \f$ I^n \f$
  /*!
    \param[in] n power of I
   */
  Complex_wp ComputePowerI(int n)
  {
    if (n >= 0)
      {
	if (n%4 == 0)
	  return Complex_wp(1.0 , 0.0);
	else if (n%4 == 1)
	  return Complex_wp(0.0 , 1.0);
	else if (n%4 == 2)
	  return Complex_wp(-1.0 , 0.0);
	else
	  return Complex_wp(0.0 , -1.0);
      }
    else
      {
	if ( (-n)%4 == 0)
	  return Complex_wp(1.0 , 0.0);
	else if ( (-n)%4 == 1)
	  return Complex_wp(0.0 , -1.0);
	else if ( (-n)%4 == 2)
	  return Complex_wp(-1.0 , 0.0);
	else
	  return Complex_wp(0.0 , 1.0);
      }
    return Complex_wp(0);
  }

#ifdef MONTJOIE_WITH_GSL
  Real_wp GetLambertW0(const Real_wp& x)
  {
    return Real_wp(gsl_sf_lambert_W0(toDouble(x)));
  }
#else
  Real_wp GetLambertW0(const Real_wp& x)
  {
#ifndef MONTJOIE_WITH_MULTIPLE

    // using Fortran function in Bessel.f
    Real_wp res;
    lambert_(&x, &res);
    return res;
#else
    cout << "Lambert function not available in multiple precision" << endl;
    abort();
    return Real_wp(0);
#endif
  }
#endif

  Real_wp GetWhittakerM(const Real_wp& Kappa, const Real_wp& mu, const Real_wp& x) 
  {
#ifdef MONTJOIE_WITH_GSL
    // equivalent of whitm in mpmath
    Real_wp res = exp(-x/2)*pow(x, mu+0.5)*gsl_sf_hyperg_1F1(mu-Kappa+0.5, 1.0+2*mu, x);
    return res;
#else
    cout << "Whittaker function not available without Gsl" << endl;
    abort();
    return Real_wp(0);
#endif
  }
  
  Real_wp GetWhittakerW(const Real_wp& Kappa, const Real_wp& mu, const Real_wp& x)
  {
#ifdef MONTJOIE_WITH_GSL
    // equivalent of whitw in mpmath
    Real_wp res = exp(-x/2)*pow(x, mu+0.5)*gsl_sf_hyperg_U(mu-Kappa+0.5, 1.0+2*mu, x);
    return res;
#else
    cout << "Whittaker function not available without Gsl" << endl;
    abort();
    return Real_wp(0);
#endif
  }

  /*  Complex_wp GetWhittakerM(const Complex_wp& Kappa, const Complex_wp& mu, const Complex_wp& z)
  {
#ifdef MONTJOIE_WITH_MATHLIB
    Complex_wp res
      = wwhitm_sub_(reinterpret_cast<const void*>(&z), reinterpret_cast<const void*>(&Kappa),
		    reinterpret_cast<const void*>(&mu));

    return res;
#else
    cout << "Whittaker function not available without Mathlib" << endl;
    abort();
    return Complex_wp(0);
#endif
  }
  
  Complex_wp GetWhittakerW(const Complex_wp& Kappa, const Complex_wp& mu, const Complex_wp& x)
  {
#ifdef MONTJOIE_WITH_MATHLIB
#else
    cout << "Whittaker function not available without Mathlib" << endl;
    abort();
    return Complex_wp(0);
#endif

}*/
  
} // end namespace

#define MONTJOIE_FILE_BESSEL_FUNCTIONS_INTERFACE_CXX
#endif
