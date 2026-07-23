#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Output/MeshInterpolationFEM.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class MeshInterpolationFEM<Dimension2>;
  SELDON_EXTERN template class MeshInterpolationFEM<Dimension3>;


  SELDON_EXTERN template void MeshInterpolationFEM<Dimension2>::ComputeEnHnOnBoundary(const VarProblem<Dimension2>&, const Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&, bool, bool) const;

  SELDON_EXTERN template void MeshInterpolationFEM<Dimension2>::ComputeEnHnOnBoundary(const VarProblem<Dimension2>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, Vector<Complex_wp>&, bool, bool) const;


  SELDON_EXTERN template void MeshInterpolationFEM<Dimension3>::ComputeEnHnOnBoundary(const VarProblem<Dimension3>&, const Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&, bool, bool) const;

  SELDON_EXTERN template void MeshInterpolationFEM<Dimension3>::ComputeEnHnOnBoundary(const VarProblem<Dimension3>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, Vector<Complex_wp>&, bool, bool) const;

  
  SELDON_EXTERN template void MeshInterpolationFEM<Dimension2>::AddSourceBoundary(const VarProblem<Dimension2>&, const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void MeshInterpolationFEM<Dimension2>::AddSourceBoundary(const VarProblem<Dimension2>&, const VectComplex_wp&, VectComplex_wp&) const;


  SELDON_EXTERN template void MeshInterpolationFEM<Dimension3>::AddSourceBoundary(const VarProblem<Dimension3>&, const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void MeshInterpolationFEM<Dimension3>::AddSourceBoundary(const VarProblem<Dimension3>&, const VectComplex_wp&, VectComplex_wp&) const;
  
}
