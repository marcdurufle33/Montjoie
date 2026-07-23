#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Source/DefineSourceElliptic.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void VarSourceProblem_Base::InitGaussianParameter(GaussianSource<Dimension2>&, const VectString&, int&) const;
  SELDON_EXTERN template void VarSourceProblem_Base::InitGaussianParameter(GaussianSource<Dimension3>&, const VectString&, int&) const;

  SELDON_EXTERN template void VarSourceProblem_Base::InitRandomGaussianParameter(GaussianSource<Dimension2>&, const VectString&, int) const;
  SELDON_EXTERN template void VarSourceProblem_Base::InitRandomGaussianParameter(GaussianSource<Dimension3>&, const VectString&, int) const;

  SELDON_EXTERN template class VarSourceProblem_Dim<Dimension2>;
  SELDON_EXTERN template class VarSourceProblem_Dim<Dimension3>;

  SELDON_EXTERN template class VarSourceProblem_Cplx<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarSourceProblem_Cplx<Real_wp, Dimension3>;
  SELDON_EXTERN template class VarSourceProblem_Cplx<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VarSourceProblem_Cplx<Complex_wp, Dimension3>;

  SELDON_EXTERN template class VarSourceProblem_Fem<Dimension2>;
  SELDON_EXTERN template class VarSourceProblem_Fem<Dimension3>;

  SELDON_EXTERN template class VarSourceProblem<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarSourceProblem<Real_wp, Dimension3>;

  SELDON_EXTERN template class VarSourceProblem<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VarSourceProblem<Complex_wp, Dimension3>;

  SELDON_EXTERN template void VarSourceProblem_Fem<Dimension2>::
  AddDiracSourceGen(const Real_wp&, const VectR2&, const Vector<VectReal_wp>&, const IVect&, Vector<VectReal_wp>&, Vector<VirtualSourceFEM<Real_wp, Dimension2>* >&) const;

  SELDON_EXTERN template void VarSourceProblem_Fem<Dimension3>::
  AddDiracSourceGen(const Real_wp&, const VectR3&, const Vector<VectReal_wp>&, const IVect&, Vector<VectReal_wp>&, Vector<VirtualSourceFEM<Real_wp, Dimension3>* >&) const;


  SELDON_EXTERN template void VarSourceProblem_Fem<Dimension2>::
  AddDiracSourceGen(const Complex_wp&, const VectR2&, const Vector<VectComplex_wp>&, const IVect&, Vector<VectComplex_wp>&, Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >&) const;

  SELDON_EXTERN template void VarSourceProblem_Fem<Dimension3>::
  AddDiracSourceGen(const Complex_wp&, const VectR3&, const Vector<VectComplex_wp>&, const IVect&, Vector<VectComplex_wp>&, Vector<VirtualSourceFEM<Complex_wp, Dimension3>* >&) const;

  SELDON_EXTERN template void VarSourceProblem_Fem<Dimension2>::
  AddDiracSourceGen(const Real_wp&, Vector<VectReal_wp>&, Vector<VirtualSourceFEM<Real_wp, Dimension2>* >&) const;

  SELDON_EXTERN template void VarSourceProblem_Fem<Dimension2>::
  AddDiracSourceGen(const Complex_wp&, Vector<VectComplex_wp>&, Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >&) const;


  SELDON_EXTERN template void VarSourceProblem_Fem<Dimension3>::
  AddDiracSourceGen(const Real_wp&, Vector<VectReal_wp>&, Vector<VirtualSourceFEM<Real_wp, Dimension3>* >&) const;

  SELDON_EXTERN template void VarSourceProblem_Fem<Dimension3>::
  AddDiracSourceGen(const Complex_wp&, Vector<VectComplex_wp>&, Vector<VirtualSourceFEM<Complex_wp, Dimension3>* >&) const;

  SELDON_EXTERN template void VarSourceProblem_Fem<Dimension2>::
  AddDiracSourceAxisym(const Real_wp&, const Vector<bool>&, int, Vector<VectReal_wp>&, const VectR3&, Vector<VirtualSourceFEM<Real_wp, Dimension2>* >&) const;

  SELDON_EXTERN template void VarSourceProblem_Fem<Dimension2>::
  AddDiracSourceAxisym(const Complex_wp&, const Vector<bool>&, int, Vector<VectComplex_wp>&, const VectR3&, Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >&) const;
  
}
