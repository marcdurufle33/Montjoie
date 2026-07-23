#include "MontjoieFlag.hxx"

#include "Elliptic/Helmholtz/MontjoieHelmholtzHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtzInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/VarLaplace.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"
#include "Computation/ProdMatVectScalarH1.cxx"
#include "Computation/ProdMatVectH1.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class LaplaceEquation<Dimension2>;
  SELDON_EXTERN template class LaplaceEquation<Dimension3>;
  SELDON_EXTERN template class LaplaceEquationDG<Dimension2>;
  SELDON_EXTERN template class LaplaceEquationDG<Dimension3>;

  SELDON_EXTERN template void LaplaceEquationDG<Dimension2>::GetNabc(TinyMatrix<Real_wp, General, 3, 3>&, TinyVector<Real_wp, 2>&, int, int, int, GlobalGenericMatrix<Real_wp> const&, int, EllipticProblem<LaplaceEquationDG<Dimension2> > const&, ElementReference<Dimension2, 1> const&);
  
  SELDON_EXTERN template void LaplaceEquationDG<Dimension3>::GetNabc(TinyMatrix<Real_wp, General, 4, 4>&, TinyVector<Real_wp, 3>&, int, int, int, GlobalGenericMatrix<Real_wp> const&, int, EllipticProblem<LaplaceEquationDG<Dimension3> > const&, ElementReference<Dimension3, 1> const&);

  SELDON_EXTERN template void LaplaceEquationDG<Dimension2>::GetPenalDG(TinyMatrix<Real_wp, General, 3, 3>&, TinyVector<Real_wp, 2>&, int, int, int, GlobalGenericMatrix<Real_wp> const&, int, int, EllipticProblem<LaplaceEquationDG<Dimension2> > const&, ElementReference<Dimension2, 1> const&);

  SELDON_EXTERN template void LaplaceEquationDG<Dimension3>::GetPenalDG(TinyMatrix<Real_wp, General, 4, 4>&, TinyVector<Real_wp, 3>&, int, int, int, GlobalGenericMatrix<Real_wp> const&, int, int, EllipticProblem<LaplaceEquationDG<Dimension3> > const&, ElementReference<Dimension3, 1> const&);

  SELDON_EXTERN template void LaplaceEquationDG<Dimension2>::GetNabc(TinyMatrix<Complex_wp, General, 3, 3>&, TinyVector<Real_wp, 2>&, int, int, int, GlobalGenericMatrix<Complex_wp> const&, int, EllipticProblem<LaplaceEquationDG<Dimension2> > const&, ElementReference<Dimension2, 1> const&);
  
  SELDON_EXTERN template void LaplaceEquationDG<Dimension3>::GetNabc(TinyMatrix<Complex_wp, General, 4, 4>&, TinyVector<Real_wp, 3>&, int, int, int, GlobalGenericMatrix<Complex_wp> const&, int, EllipticProblem<LaplaceEquationDG<Dimension3> > const&, ElementReference<Dimension3, 1> const&);

  SELDON_EXTERN template void LaplaceEquationDG<Dimension2>::GetPenalDG(TinyMatrix<Complex_wp, General, 3, 3>&, TinyVector<Real_wp, 2>&, int, int, int, GlobalGenericMatrix<Complex_wp> const&, int, int, EllipticProblem<LaplaceEquationDG<Dimension2> > const&, ElementReference<Dimension2, 1> const&);

  SELDON_EXTERN template void LaplaceEquationDG<Dimension3>::GetPenalDG(TinyMatrix<Complex_wp, General, 4, 4>&, TinyVector<Real_wp, 3>&, int, int, int, GlobalGenericMatrix<Complex_wp> const&, int, int, EllipticProblem<LaplaceEquationDG<Dimension3> > const&, ElementReference<Dimension3, 1> const&);

  SELDON_EXTERN template void LaplaceEquationDG<Dimension2>::GetTensorMass(EllipticProblem<LaplaceEquationDG<Dimension2> > const&, int, int, GlobalGenericMatrix<Real_wp> const&, int, TinyMatrix<Real_wp, General, 3, 3>&);
  SELDON_EXTERN template void LaplaceEquationDG<Dimension3>::GetTensorMass(EllipticProblem<LaplaceEquationDG<Dimension3> > const&, int, int, GlobalGenericMatrix<Real_wp> const&, int, TinyMatrix<Real_wp, General, 4, 4>&);

  SELDON_EXTERN template class VarHarmonic<LaplaceEquation<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<LaplaceEquation<Dimension3> >;
  SELDON_EXTERN template class VarHarmonic<LaplaceEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<LaplaceEquationDG<Dimension3> >;

  SELDON_EXTERN template class VarLaplace_Base<LaplaceEquation<Dimension2> >;
  SELDON_EXTERN template class VarLaplace_Base<LaplaceEquation<Dimension3> >;
  SELDON_EXTERN template class VarLaplace_Base<LaplaceEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarLaplace_Base<LaplaceEquationDG<Dimension3> >;

  SELDON_EXTERN template class VarComputationRCS_Laplace<Dimension2>;
  SELDON_EXTERN template class VarComputationRCS_Laplace<Dimension3>;

  SELDON_EXTERN template class IncidentWaveProjector_Laplace<Dimension2>;
  SELDON_EXTERN template class IncidentWaveProjector_Laplace<Dimension3>;

  SELDON_EXTERN template class DiffractedWaveSource_Laplace<Dimension2>;
  SELDON_EXTERN template class DiffractedWaveSource_Laplace<Dimension3>;

  SELDON_EXTERN template class TotalWaveSource_Laplace<Dimension2>;
  SELDON_EXTERN template class TotalWaveSource_Laplace<Dimension3>;

  SELDON_EXTERN template class IncidentWaveProjector_LaplaceDG<Dimension2>;
  SELDON_EXTERN template class IncidentWaveProjector_LaplaceDG<Dimension3>;

  SELDON_EXTERN template class DiffractedWaveSource_LaplaceDG<Dimension2>;
  SELDON_EXTERN template class DiffractedWaveSource_LaplaceDG<Dimension3>;

  SELDON_EXTERN template class TotalWaveSource_LaplaceDG<Dimension2>;
  SELDON_EXTERN template class TotalWaveSource_LaplaceDG<Dimension3>;

  SELDON_EXTERN template class EllipticProblem<LaplaceEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<LaplaceEquation<Dimension3> >;
  SELDON_EXTERN template class EllipticProblem<LaplaceEquationDG<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<LaplaceEquationDG<Dimension3> >;

  SELDON_EXTERN template void GetExtrapolationU(EllipticProblem<LaplaceEquation<Dimension2> > const&, int, TinyVector<Vector<Real_wp>, 2> const&, TinyVector<Vector<Real_wp>, 2>&, TinyVector<Vector<Real_wp>, 2>&, ElementReference<Dimension2, 1> const&, TinyVector<Vector<Real_wp>, 2>&, bool);
  SELDON_EXTERN template void GetExtrapolationU(EllipticProblem<LaplaceEquation<Dimension3> > const&, int, TinyVector<Vector<Real_wp>, 3> const&, TinyVector<Vector<Real_wp>, 3>&, TinyVector<Vector<Real_wp>, 3>&, ElementReference<Dimension3, 1> const&, TinyVector<Vector<Real_wp>, 3>&, bool);

  SELDON_EXTERN template void GetExtrapolationU(EllipticProblem<LaplaceEquationDG<Dimension2> > const&, int, TinyVector<Vector<Real_wp>, 2> const&, TinyVector<Vector<Real_wp>, 2>&, TinyVector<Vector<Real_wp>, 2>&, ElementReference<Dimension2, 1> const&, TinyVector<Vector<Real_wp>, 2>&, bool);
  SELDON_EXTERN template void GetExtrapolationU(EllipticProblem<LaplaceEquationDG<Dimension3> > const&, int, TinyVector<Vector<Real_wp>, 3> const&, TinyVector<Vector<Real_wp>, 3>&, TinyVector<Vector<Real_wp>, 3>&, ElementReference<Dimension3, 1> const&, TinyVector<Vector<Real_wp>, 3>&, bool);

  SELDON_EXTERN template void MltAdd_SquareHex_VectorialH1(Real_wp const&, int, FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension2> > const&, Vector<Real_wp> const&, Real_wp const&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd_SquareHex_VectorialH1(Real_wp const&, int, FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension3> > const&, Vector<Real_wp> const&, Real_wp const&, Vector<Real_wp>&);

  SELDON_EXTERN template void MltAdd_SquareHex_ScalarH1(Real_wp const&, int, FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension2> > const&, Vector<Real_wp> const&, Real_wp const&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd_SquareHex_ScalarH1(Real_wp const&, int, FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension3> > const&, Vector<Real_wp> const&, Real_wp const&, Vector<Real_wp>&);
  
  SELDON_EXTERN template void MltAddScalarH1(const Real_wp&, const GlobalGenericMatrix<Real_wp>&, const SeldonTranspose&, int, const FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension2> >&, const VectReal_wp&, const Real_wp&, VectReal_wp&, bool);
  SELDON_EXTERN template void MltAddScalarH1(const Real_wp&, const GlobalGenericMatrix<Real_wp>&, const SeldonTranspose&, int, const FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension3> >&, const VectReal_wp&, const Real_wp&, VectReal_wp&, bool);
  
  SELDON_EXTERN template void AddElementaryFluxesDG(VirtualMatrix<Real_wp>&, GlobalGenericMatrix<Real_wp> const&, EllipticProblem<LaplaceEquation<Dimension2> > const&, int, int);
  SELDON_EXTERN template void AddElementaryFluxesDG(VirtualMatrix<Real_wp>&, GlobalGenericMatrix<Real_wp> const&, EllipticProblem<LaplaceEquation<Dimension3> > const&, int, int);

}
