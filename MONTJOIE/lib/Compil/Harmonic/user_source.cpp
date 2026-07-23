#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Source/UserSource.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class UserDefinedSource<Real_wp, Dimension2>;
  SELDON_EXTERN template class UserDefinedSource<Complex_wp, Dimension2>;
  SELDON_EXTERN template class UserDefinedSource<Real_wp, Dimension3>;
  SELDON_EXTERN template class UserDefinedSource<Complex_wp, Dimension3>;

  SELDON_EXTERN template class InitialUserFunction<Real_wp, Dimension2>;
  SELDON_EXTERN template class InitialUserFunction<Real_wp, Dimension3>;

  SELDON_EXTERN template class VirtualTimeSource<Real_wp>;

  SELDON_EXTERN template void ComputeVariableUserIndex(const VarGeometryProblem<Dimension2>&, const IVect&, int, int, int, const Vector<VectR2>&, bool, bool,
						       Vector<Vector<Real_wp> >&, Vector<Vector<TinyVector<Real_wp, 2> > >& grad_rho,
						       Vector<Vector<TinyMatrix<Real_wp, Symmetric, 2, 2> > >&, const Real_wp&, const Real_wp&);

  SELDON_EXTERN template void ComputeVariableUserIndex(const VarGeometryProblem<Dimension3>&, const IVect&, int, int, int, const Vector<VectR3>&, bool, bool,
						       Vector<Vector<Real_wp> >&, Vector<Vector<TinyVector<Real_wp, 3> > >& grad_rho,
						       Vector<Vector<TinyMatrix<Real_wp, Symmetric, 3, 3> > >&, const Real_wp&, const Real_wp&);

  SELDON_EXTERN template void ComputeVariableUserIndex(const VarGeometryProblem<Dimension2>&, const IVect&, int, int, int, const Vector<VectR2>&, bool, bool,
						       Vector<Vector<Complex_wp> >&, Vector<Vector<TinyVector<Complex_wp, 2> > >& grad_rho,
						       Vector<Vector<TinyMatrix<Complex_wp, Symmetric, 2, 2> > >&, const Complex_wp&, const Complex_wp&);

  SELDON_EXTERN template void ComputeVariableUserIndex(const VarGeometryProblem<Dimension3>&, const IVect&, int, int, int, const Vector<VectR3>&, bool, bool,
						       Vector<Vector<Complex_wp> >&, Vector<Vector<TinyVector<Complex_wp, 3> > >& grad_rho,
						       Vector<Vector<TinyMatrix<Complex_wp, Symmetric, 3, 3> > >&, const Complex_wp&, const Complex_wp&);
  
}
