#include "MontjoieFlag.hxx"

#include "Elliptic/Helmholtz/MontjoieHelmholtzHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtzInline.hxx"

#include "Elliptic/Helmholtz/CellulairePeriodique.hxx"

#include "Hyperbolic/Acoustic/MontjoieAcousticHeader.hxx"
#include "Hyperbolic/Acoustic/MontjoieAcousticInline.hxx"

#include "Hyperbolic/Acoustic/VarAdvection.hxx"
#include "Hyperbolic/Acoustic/VarAdvectionInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/VarHelmholtz.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class VarHelmholtzIndex_Base<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarHelmholtzIndex_Base<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VarHelmholtzIndex_Base<Real_wp, Dimension3>;
  SELDON_EXTERN template class VarHelmholtzIndex_Base<Complex_wp, Dimension3>;
  
  SELDON_EXTERN template class HelmholtzEquation<Dimension2>;
  SELDON_EXTERN template class HelmholtzEquation<Dimension3>;

  SELDON_EXTERN template class VarHarmonic<HelmholtzEquation<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<HelmholtzEquation<Dimension3> >;
  SELDON_EXTERN template class VarHarmonic<HelmholtzEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<HelmholtzEquationDG<Dimension3> >;

  SELDON_EXTERN template class VarHelmholtz_Eq<LaplaceEquation<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtz_Eq<LaplaceEquation<Dimension3> >;
  SELDON_EXTERN template class VarHelmholtz_Eq<LaplaceEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtz_Eq<LaplaceEquationDG<Dimension3> >;

  SELDON_EXTERN template class VarHelmholtz_Eq<HarmonicAdvectionEquation<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtz_Eq<HarmonicAdvectionEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtz_Eq<AdvectionEquationDG<Dimension2> >;
  //SELDON_EXTERN template class VarHelmholtz_Eq<HarmonicDoubleAdvectionEquation<Dimension2> >;
  //SELDON_EXTERN template class VarHelmholtz_Eq<HarmonicDoubleAdvectionEquationDG<Dimension2> >;

  SELDON_EXTERN template class VarHarmonic<HelmholtzEquation_CellulairePeriodique<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<HelmholtzEquation_CellulairePeriodique<Dimension3> >;

  SELDON_EXTERN template class VarHelmholtz_Eq<HelmholtzEquation_CellulairePeriodique<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtz_Eq<HelmholtzEquation_CellulairePeriodique<Dimension3> >;
  
  SELDON_EXTERN template class VarHelmholtz_Base<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarHelmholtz_Base<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VarHelmholtz_Base<Real_wp, Dimension3>;
  SELDON_EXTERN template class VarHelmholtz_Base<Complex_wp, Dimension3>;

  SELDON_EXTERN template void VarHelmholtz_Base<Complex_wp, Dimension2>::GetPenalizationTauIPP(const VarProblem<Dimension2>&, int, Complex_wp&) const;
  SELDON_EXTERN template void VarHelmholtz_Base<Complex_wp, Dimension3>::GetPenalizationTauIPP(const VarProblem<Dimension3>&, int, Complex_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension2>::GetPenalizationTauIPP(const VarProblem<Dimension2>&, int, Complex_wp&) const;
  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension3>::GetPenalizationTauIPP(const VarProblem<Dimension3>&, int, Complex_wp&) const;
  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension2>::GetPenalizationTauIPP(const VarProblem<Dimension2>&, int, Real_wp&) const;
  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension3>::GetPenalizationTauIPP(const VarProblem<Dimension3>&, int, Real_wp&) const;
  
  SELDON_EXTERN template void VarHelmholtz_Base<Complex_wp, Dimension2>::GetCoefficientIPP(const VarProblem<Dimension2>&, Complex_wp&, Real_wp&,
                                                                                           Complex_wp&, Complex_wp&, Complex_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Complex_wp, Dimension3>::GetCoefficientIPP(const VarProblem<Dimension3>&, Complex_wp&, Real_wp&,
                                                                                           Complex_wp&, Complex_wp&, Complex_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension2>::GetCoefficientIPP(const VarProblem<Dimension2>&, Real_wp&, Real_wp&,
                                                                                        Real_wp&, Real_wp&, Real_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension3>::GetCoefficientIPP(const VarProblem<Dimension3>&, Real_wp&, Real_wp&,
                                                                                        Real_wp&, Real_wp&, Real_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension2>::GetCoefficientIPP(const VarProblem<Dimension2>&, Real_wp&, Real_wp&,
                                                                                        Complex_wp&, Complex_wp&, Complex_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension3>::GetCoefficientIPP(const VarProblem<Dimension3>&, Real_wp&, Real_wp&,
                                                                                        Complex_wp&, Complex_wp&, Complex_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Complex_wp, Dimension2>::ModifyCoefficientBC_HDG(const VarProblem<Dimension2>&, const VarBoundaryCondition_Dim<Complex_wp, Dimension2>&, const Mesh<Dimension2>&, int, int, const Complex_wp&, const Real_wp&, const Complex_wp&, int, int, int, Complex_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Complex_wp, Dimension3>::ModifyCoefficientBC_HDG(const VarProblem<Dimension3>&, const VarBoundaryCondition_Dim<Complex_wp, Dimension3>&, const Mesh<Dimension3>&, int, int, const Complex_wp&, const Real_wp&, const Complex_wp&, int, int, int, Complex_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension2>::ModifyCoefficientBC_HDG(const VarProblem<Dimension2>&, const VarBoundaryCondition_Dim<Real_wp, Dimension2>&, const Mesh<Dimension2>&, int, int, const Complex_wp&, const Real_wp&, const Complex_wp&, int, int, int, Complex_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension3>::ModifyCoefficientBC_HDG(const VarProblem<Dimension3>&, const VarBoundaryCondition_Dim<Real_wp, Dimension3>&, const Mesh<Dimension3>&, int, int, const Complex_wp&, const Real_wp&, const Complex_wp&, int, int, int, Complex_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension2>::ModifyCoefficientBC_HDG(const VarProblem<Dimension2>&, const VarBoundaryCondition_Dim<Real_wp, Dimension2>&, const Mesh<Dimension2>&, int, int, const Real_wp&, const Real_wp&, const Real_wp&, int, int, int, Real_wp&) const;
  
  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension3>::ModifyCoefficientBC_HDG(const VarProblem<Dimension3>&, const VarBoundaryCondition_Dim<Real_wp, Dimension3>&, const Mesh<Dimension3>&, int, int, const Real_wp&, const Real_wp&, const Real_wp&, int, int, int, Real_wp&) const;

  SELDON_EXTERN template void VarHelmholtz_Base<Complex_wp, Dimension2>::ComputeElementaryMatrixHelmDG(int iquad, IVect&, VirtualMatrix<Complex_wp>& mat_interac,
                                                                                                       CondensationBlockSolver_Base<Complex_wp>&,
                                                                                                       const GlobalGenericMatrix<Complex_wp>& nat_mat,
                                                                                                       const VarProblem<Dimension2>&,
                                                                                                       const VarBoundaryCondition_Dim<Complex_wp, Dimension2>&,
                                                                                                       const ElementReference<Dimension2, 1>& Fb);

    SELDON_EXTERN template void VarHelmholtz_Base<Complex_wp, Dimension3>::ComputeElementaryMatrixHelmDG(int iquad, IVect&, VirtualMatrix<Complex_wp>& mat_interac,
                                                                                                       CondensationBlockSolver_Base<Complex_wp>&,
                                                                                                       const GlobalGenericMatrix<Complex_wp>& nat_mat,
                                                                                                       const VarProblem<Dimension3>&,
                                                                                                       const VarBoundaryCondition_Dim<Complex_wp, Dimension3>&,
                                                                                                       const ElementReference<Dimension3, 1>& Fb);
  
  SELDON_EXTERN template class VarHelmholtz_Eq<HelmholtzEquation<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtz_Eq<HelmholtzEquation<Dimension3> >;
  SELDON_EXTERN template class VarHelmholtz_Eq<HelmholtzEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtz_Eq<HelmholtzEquationDG<Dimension3> >;

  SELDON_EXTERN template class VarHelmholtz_Cplx<HelmholtzEquation<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtz_Cplx<HelmholtzEquation<Dimension3> >;
  SELDON_EXTERN template class VarHelmholtz_Cplx<HelmholtzEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtz_Cplx<HelmholtzEquationDG<Dimension3> >;

  SELDON_EXTERN template class CondensationBlockSolver_Helm<Real_wp, HelmholtzEquation<Dimension2> >;
  SELDON_EXTERN template class CondensationBlockSolver_Helm<Real_wp, HelmholtzEquation<Dimension3> >;
  SELDON_EXTERN template class CondensationBlockSolver_Helm<Complex_wp, HelmholtzEquation<Dimension2> >;
  SELDON_EXTERN template class CondensationBlockSolver_Helm<Complex_wp, HelmholtzEquation<Dimension3> >;

  SELDON_EXTERN template class CondensationBlockSolver_Helm<Real_wp, LaplaceEquation<Dimension2> >;
  SELDON_EXTERN template class CondensationBlockSolver_Helm<Real_wp, LaplaceEquation<Dimension3> >;
  SELDON_EXTERN template class CondensationBlockSolver_Helm<Complex_wp, LaplaceEquation<Dimension2> >;
  SELDON_EXTERN template class CondensationBlockSolver_Helm<Complex_wp, LaplaceEquation<Dimension3> >;

  SELDON_EXTERN template class CondensationBlockSolver_HelmDG<Real_wp, Complex_wp, Dimension2>;
  SELDON_EXTERN template class CondensationBlockSolver_HelmDG<Real_wp, Complex_wp, Dimension3>;
  SELDON_EXTERN template class CondensationBlockSolver_HelmDG<Complex_wp, Complex_wp, Dimension2>;
  SELDON_EXTERN template class CondensationBlockSolver_HelmDG<Complex_wp, Complex_wp, Dimension3>;

  SELDON_EXTERN template class CondensationBlockSolver_HelmDG<Real_wp, Real_wp, Dimension2>;
  SELDON_EXTERN template class CondensationBlockSolver_HelmDG<Real_wp, Real_wp, Dimension3>;
  SELDON_EXTERN template class CondensationBlockSolver_HelmDG<Complex_wp, Real_wp, Dimension2>;
  SELDON_EXTERN template class CondensationBlockSolver_HelmDG<Complex_wp, Real_wp, Dimension3>;

  SELDON_EXTERN template class HelmholtzEquationDG<Dimension2>;
  SELDON_EXTERN template class HelmholtzEquationDG<Dimension3>;

  SELDON_EXTERN template class VarComputationRCS_Helm<Dimension2>;
  SELDON_EXTERN template class VarComputationRCS_Helm<Dimension3>;

  SELDON_EXTERN template class TransparencySolver_Helm<Dimension2>;
  SELDON_EXTERN template class TransparencySolver_Helm<Dimension3>;

  SELDON_EXTERN template class EllipticProblem<HelmholtzEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HelmholtzEquation<Dimension3> >;
  SELDON_EXTERN template class EllipticProblem<HelmholtzEquationDG<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HelmholtzEquationDG<Dimension3> >;

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension2>
  ::AddElementaryFluxesHelmDG(VirtualMatrix<Real_wp>&, const GlobalGenericMatrix<Real_wp>&,
			      const EllipticProblem<LaplaceEquationDG<Dimension2> >&, int, int);

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension2>
  ::AddElementaryFluxesHelmDG(VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&,
			      const EllipticProblem<LaplaceEquationDG<Dimension2> >&, int, int);

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension3>
  ::AddElementaryFluxesHelmDG(VirtualMatrix<Real_wp>&, const GlobalGenericMatrix<Real_wp>&,
			      const EllipticProblem<LaplaceEquationDG<Dimension3> >&, int, int);

  SELDON_EXTERN template void VarHelmholtz_Base<Real_wp, Dimension3>
  ::AddElementaryFluxesHelmDG(VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&,
			      const EllipticProblem<LaplaceEquationDG<Dimension3> >&, int, int);

  SELDON_EXTERN template void VarHelmholtz_Base<Complex_wp, Dimension2>
  ::AddElementaryFluxesHelmDG(VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&,
			      const EllipticProblem<HelmholtzEquationDG<Dimension2> >&, int, int);

  SELDON_EXTERN template void VarHelmholtz_Base<Complex_wp, Dimension3>
  ::AddElementaryFluxesHelmDG(VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&,
			      const EllipticProblem<HelmholtzEquationDG<Dimension3> >&, int, int);

  SELDON_EXTERN template void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>&, GlobalGenericMatrix<Complex_wp> const&, EllipticProblem<HelmholtzEquation<Dimension2> > const&, int, int);
  SELDON_EXTERN template void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>&, GlobalGenericMatrix<Complex_wp> const&, EllipticProblem<HelmholtzEquation<Dimension3> > const&, int, int);
  
}
