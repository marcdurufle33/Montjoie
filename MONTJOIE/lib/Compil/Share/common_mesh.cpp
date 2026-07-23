#include "MontjoieFlag.hxx"

// pour un InitMontjoie qui marche bien (avec initialisation de Mesh)
#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Share/MontjoieTypes.cxx"
#include "Share/CommonMontjoie.cxx"
#include "Share/Precision.cxx"

#ifdef MONTJOIE_WITH_MPFR
#include "Share/PrecisionMpfr.cxx"
#endif

#ifdef MONTJOIE_USE_FLOAT128
#include "Share/PrecisionQuad.cxx"
#endif

#endif

namespace Seldon
{
  SELDON_EXTERN template void Copy(list<TinyVector<Real_wp, 3> >&, Vector<TinyVector<Real_wp, 3> >&);
  SELDON_EXTERN template void Copy(list<int>&, Vector<int>&);
  SELDON_EXTERN template void Copy(list<Real_wp>&, Vector<Real_wp>&);  
  SELDON_EXTERN template void Copy(const Vector<Real_wp>&, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void Copy(const Vector<Complex_wp >&, Vector<Complex_wp, VectSparse>&);
  
  SELDON_EXTERN template void CopyVector(const Vector<Real_wp>&, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void Linspace(const Montjoie::Real_wp&, const Montjoie::Real_wp&, int, Vector<Montjoie::Real_wp>&);

}

namespace Montjoie
{  
  SELDON_EXTERN template void ComputePadeCoefficientsSqrt(const Real_wp&, int, complex<Real_wp>&, Vector<complex<Real_wp> >&, Vector<complex<Real_wp> >&);
  
}
