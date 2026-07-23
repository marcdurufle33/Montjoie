#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Source/SourceSpatiale.cxx"
#include "vector/Vector.cxx"
#include "vector/Functions_Arrays.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void Vector<Montjoie::ParamOutputClass>::Resize(size_t);
 
  SELDON_EXTERN template void Vector<Montjoie::MeshInterpolation<Montjoie::Dimension2> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::MeshInterpolation<Montjoie::Dimension3> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::GridInterpolationFull<Montjoie::Dimension1> >::Resize(size_t);  
  SELDON_EXTERN template void Vector<Montjoie::GridInterpolationFull<Montjoie::Dimension2> >::Resize(size_t);  
  SELDON_EXTERN template void Vector<Montjoie::GridInterpolationFull<Montjoie::Dimension3> >::Resize(size_t);  

  SELDON_EXTERN template void Vector<Montjoie::PhysicalVaryingMedia<Montjoie::Dimension2, Real_wp>* >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::PhysicalVaryingMedia<Montjoie::Dimension2, Complex_wp>* >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::PhysicalVaryingMedia<Montjoie::Dimension3, Real_wp>* >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::PhysicalVaryingMedia<Montjoie::Dimension3, Complex_wp>* >::Resize(size_t);

  SELDON_EXTERN template void RemoveDuplicate(Vector<Montjoie::VirtualSourceField<Complex_wp, Montjoie::Dimension2>* >&);
  SELDON_EXTERN template void RemoveDuplicate(Vector<Montjoie::VirtualSourceField<Real_wp, Montjoie::Dimension2>* >&);
  SELDON_EXTERN template void RemoveDuplicate(Vector<Montjoie::VirtualSourceField<Complex_wp, Montjoie::Dimension3>* >&);
  SELDON_EXTERN template void RemoveDuplicate(Vector<Montjoie::VirtualSourceField<Real_wp, Montjoie::Dimension3>* >&);

}

namespace Montjoie
{
  SELDON_EXTERN template class VirtualSourceField<Real_wp, Dimension2>;
  SELDON_EXTERN template class VirtualSourceField<Real_wp, Dimension3>;
  SELDON_EXTERN template class VirtualSourceField<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VirtualSourceField<Complex_wp, Dimension3>;

  SELDON_EXTERN template class GaussianSource<Dimension2>;
  SELDON_EXTERN template class GaussianSource<Dimension3>;

  SELDON_EXTERN template class GaussianSourceField<Real_wp, Dimension2>;
  SELDON_EXTERN template class GaussianSourceField<Real_wp, Dimension3>;
  SELDON_EXTERN template class GaussianSourceField<Complex_wp, Dimension2>;
  SELDON_EXTERN template class GaussianSourceField<Complex_wp, Dimension3>;
  
  SELDON_EXTERN template class UniformSourceField<Real_wp, Dimension2>;
  SELDON_EXTERN template class UniformSourceField<Real_wp, Dimension3>;
  SELDON_EXTERN template class UniformSourceField<Complex_wp, Dimension2>;
  SELDON_EXTERN template class UniformSourceField<Complex_wp, Dimension3>;

  SELDON_EXTERN template class IncidentWaveField<Real_wp, Dimension2>;
  SELDON_EXTERN template class IncidentWaveField<Real_wp, Dimension3>;
  SELDON_EXTERN template class IncidentWaveField<Complex_wp, Dimension2>;
  SELDON_EXTERN template class IncidentWaveField<Complex_wp, Dimension3>;

  SELDON_EXTERN template class PlaneWaveIncidentField<Real_wp, Dimension2>;
  SELDON_EXTERN template class PlaneWaveIncidentField<Real_wp, Dimension3>;
  SELDON_EXTERN template class PlaneWaveIncidentField<Complex_wp, Dimension2>;
  SELDON_EXTERN template class PlaneWaveIncidentField<Complex_wp, Dimension3>;

  SELDON_EXTERN template class PlaneWaveIncidentFieldComplex<Real_wp, Dimension2>;
  SELDON_EXTERN template class PlaneWaveIncidentFieldComplex<Real_wp, Dimension3>;
  SELDON_EXTERN template class PlaneWaveIncidentFieldComplex<Complex_wp, Dimension2>;
  SELDON_EXTERN template class PlaneWaveIncidentFieldComplex<Complex_wp, Dimension3>;

  SELDON_EXTERN template class LayeredPlaneWaveIncidentField<Real_wp, Dimension2>;
  SELDON_EXTERN template class LayeredPlaneWaveIncidentField<Real_wp, Dimension3>;
  SELDON_EXTERN template class LayeredPlaneWaveIncidentField<Complex_wp, Dimension2>;
  SELDON_EXTERN template class LayeredPlaneWaveIncidentField<Complex_wp, Dimension3>;

  SELDON_EXTERN template class HankelIncidentField<Real_wp, Dimension2>;
  SELDON_EXTERN template class HankelIncidentField<Real_wp, Dimension3>;
  SELDON_EXTERN template class HankelIncidentField<Complex_wp, Dimension2>;
  SELDON_EXTERN template class HankelIncidentField<Complex_wp, Dimension3>;

  SELDON_EXTERN template class GaussianBeamIncidentField<Real_wp, Dimension2>;
  SELDON_EXTERN template class GaussianBeamIncidentField<Real_wp, Dimension3>;
  SELDON_EXTERN template class GaussianBeamIncidentField<Complex_wp, Dimension2>;
  SELDON_EXTERN template class GaussianBeamIncidentField<Complex_wp, Dimension3>;

  SELDON_EXTERN template class IncidentWaveSourceField<Real_wp, Dimension2>;
  SELDON_EXTERN template class IncidentWaveSourceField<Real_wp, Dimension3>;
  SELDON_EXTERN template class IncidentWaveSourceField<Complex_wp, Dimension2>;
  SELDON_EXTERN template class IncidentWaveSourceField<Complex_wp, Dimension3>;

  SELDON_EXTERN template class VirtualProjectorFEM_Base<Real_wp>;
  SELDON_EXTERN template class VirtualProjectorFEM_Base<Complex_wp>;

  SELDON_EXTERN template class VirtualProjectorFEM<Real_wp, Dimension2>;
  SELDON_EXTERN template class VirtualProjectorFEM<Real_wp, Dimension3>;
  SELDON_EXTERN template class VirtualProjectorFEM<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VirtualProjectorFEM<Complex_wp, Dimension3>;

  SELDON_EXTERN template class IncidentWaveProjector<Real_wp, Dimension2>;
  SELDON_EXTERN template class IncidentWaveProjector<Real_wp, Dimension3>;
  SELDON_EXTERN template class IncidentWaveProjector<Complex_wp, Dimension2>;
  SELDON_EXTERN template class IncidentWaveProjector<Complex_wp, Dimension3>;

  SELDON_EXTERN template class VirtualSourceFEM<Real_wp, Dimension2>;
  SELDON_EXTERN template class VirtualSourceFEM<Real_wp, Dimension3>;
  SELDON_EXTERN template class VirtualSourceFEM<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VirtualSourceFEM<Complex_wp, Dimension3>;

  SELDON_EXTERN template class VolumetricSource_Base<Real_wp, Dimension2>;
  SELDON_EXTERN template class VolumetricSource_Base<Real_wp, Dimension3>;
  SELDON_EXTERN template class VolumetricSource_Base<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VolumetricSource_Base<Complex_wp, Dimension3>;

}
