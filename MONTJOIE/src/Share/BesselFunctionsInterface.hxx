#ifndef MONTJOIE_FILE_BESSEL_FUNCTIONS_INTERFACE_HXX

namespace Montjoie
{
  Real_wp GetKn(int n, const Real_wp& x);  
  Real_wp GetJn(int n, const Real_wp& x);
  Real_wp GetYn(int n, const Real_wp& x);
  
  void GetDeriveJn(int n, const Real_wp& x, Real_wp&, Real_wp&);

  void ComputeBesselFunctions(int begin_order, int n,
                              const Complex_wp& z, VectComplex_wp& Jn);

  void ComputeBesselFunctions(int begin_order, int n,
                              const Complex_wp& z, VectComplex_wp& Jn, VectComplex_wp& Yn);
  
  void ComputeBesselFunctions(int begin_order, int n, const Real_wp& x, VectReal_wp& Jn);

  void ComputeDeriveBesselFunctions(int begin_order, int n, const Real_wp& x,
                                    VectReal_wp& Jn, VectReal_wp& dJn);
  
  void ComputeBesselAndHankelFunctions(const Real_wp& begin_order, int n, const Complex_wp& z,
                                       VectComplex_wp& Jn, VectComplex_wp& Hn);

  void ComputeDeriveBesselAndHankel(const Real_wp& begin_order, int n, const Complex_wp& z,
                                    VectComplex_wp& Jn, VectComplex_wp& Hn,
                                    VectComplex_wp& dJn, VectComplex_wp& dHn);
  
  void ComputeSphericalBessel(int n, const Real_wp& x, VectReal_wp& Jn);
  void ComputeDeriveSphericalBessel(int n, const Real_wp& x, VectReal_wp& Jn, VectReal_wp& dJn);

  void ComputeSphericalHankel(int n, const Real_wp& x, VectComplex_wp& Hn);
  void ComputeDeriveSphericalHankel(int n, const Real_wp& x, VectComplex_wp& Hn,
                                    VectComplex_wp& dHn);  
  
  void ComputeDeriveSphericalBesselHankel(const Real_wp& begin_order, int n, const Complex_wp& z,
                                          VectComplex_wp& Jn, VectComplex_wp& Hn,
                                          VectComplex_wp& Jn_prime, VectComplex_wp& Hn_prime);
  
  void ComputeDeriveRiccatiBessel(const Real_wp& begin_order, int n, const Complex_wp& z,
                                  VectComplex_wp& Jn, VectComplex_wp& Hn,
                                  VectComplex_wp& Jn_prime, VectComplex_wp& Hn_prime,
                                  VectComplex_wp& Psi_n, VectComplex_wp& Xi_n,
                                  VectComplex_wp& Psi_prime, VectComplex_wp& Xi_prime);

  template<class T>
  int ComputeOrder(const T& ka, const T& epsilon_machine);
  
  Complex_wp ComputePowerI(int n);

  Real_wp GetLambertW0(const Real_wp& x);

  Real_wp GetWhittakerM(const Real_wp& Kappa, const Real_wp& mu, const Real_wp& x);
  Real_wp GetWhittakerW(const Real_wp& Kappa, const Real_wp& mu, const Real_wp& x);

  //Complex_wp GetWhittakerM(const Complex_wp& Kappa, const Complex_wp& mu, const Complex_wp& x);
  //Complex_wp GetWhittakerW(const Complex_wp& Kappa, const Complex_wp& mu, const Complex_wp& x);

} // namespace Montjoie

#define MONTJOIE_FILE_BESSEL_FUNCTIONS_INTERFACE_HXX
#endif

