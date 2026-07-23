#include "MontjoieFlag.hxx"

#include "Share/MontjoieCommonHeader.hxx"
#include "Share/MontjoieCommonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Share/FFT.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class VirtualFftObject<Complex_wp>;
  SELDON_EXTERN template class ManualFftObject<Complex_wp>;
  SELDON_EXTERN template class FftInterface<Complex_wp>;

  SELDON_EXTERN template class VirtualFftRealObject<Real_wp>;
  
#ifdef MONTJOIE_WITH_GSL
  SELDON_EXTERN template class GslFftObject<Complex_wp>;
  SELDON_EXTERN template class GslFftRealObject<Real_wp>;
#endif

#ifdef MONTJOIE_WITH_FFTW
  SELDON_EXTERN template class FftwFftObject<Complex_wp>;
  SELDON_EXTERN template class FftwFftRealObject<Real_wp>;
#endif

#ifdef SELDON_WITH_MKL
  SELDON_EXTERN template class MklFftObject<Complex_wp>;
  SELDON_EXTERN template class MklFftRealObject<Real_wp>;
#endif

}
