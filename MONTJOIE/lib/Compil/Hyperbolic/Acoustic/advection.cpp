#include "MontjoieFlag.hxx"

#include "Algebra/MontjoieAlgebraHeader.hxx"
#include "Algebra/MontjoieAlgebraInline.hxx"
#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Computation/ElementaryMatrixH1.hxx"

#include "Hyperbolic/Acoustic/MontjoieAcousticHeader.hxx"
#include "Hyperbolic/Acoustic/MontjoieAcousticInline.hxx"

#include "Hyperbolic/Acoustic/VarAdvection.hxx"
#include "Hyperbolic/Acoustic/VarAdvectionInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Acoustic/VarAdvection.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class HarmonicAdvectionEquation<Dimension2>;
  SELDON_EXTERN template class HarmonicAdvectionEquationDG<Dimension2>;
  SELDON_EXTERN template class AdvectionEquationDG<Dimension2>;
  //SELDON_EXTERN template class HarmonicDoubleAdvectionEquation<Dimension2>;
  //SELDON_EXTERN template class HarmonicDoubleAdvectionEquationDG<Dimension2>;

  SELDON_EXTERN template class VarAdvection_Base<HarmonicAdvectionEquation<Dimension2> >;
  SELDON_EXTERN template class VarAdvection_Base<HarmonicAdvectionEquationDG<Dimension2> >;
  //SELDON_EXTERN template class VarAdvection_Base<HarmonicDoubleAdvectionEquation<Dimension2> >;
  //SELDON_EXTERN template class VarAdvection_Base<HarmonicDoubleAdvectionEquationDG<Dimension2> >;

  SELDON_EXTERN template class EllipticProblem<HarmonicAdvectionEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HarmonicAdvectionEquationDG<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<AdvectionEquationDG<Dimension2> >;
  //SELDON_EXTERN template class EllipticProblem<HarmonicDoubleAdvectionEquation<Dimension2> >;
  //SELDON_EXTERN template class EllipticProblem<HarmonicDoubleAdvectionEquationDG<Dimension2> >;

  SELDON_EXTERN template class VarHarmonic<HarmonicAdvectionEquation<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<HarmonicAdvectionEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<AdvectionEquationDG<Dimension2> >;
  //SELDON_EXTERN template class VarHarmonic<HarmonicDoubleAdvectionEquation<Dimension2> >;
  //SELDON_EXTERN template class VarHarmonic<HarmonicDoubleAdvectionEquationDG<Dimension2> >;

  SELDON_EXTERN template class AdvectionEquation_Base<Complex_wp, Dimension2>;

  SELDON_EXTERN template void AdvectionEquation_Base<Complex_wp, Dimension2>::GetGradPhiTensor(const EllipticProblem<HarmonicAdvectionEquationDG<Dimension2> >&, int, int, const GlobalGenericMatrix<Complex_wp>&, int, TinyMatrix<TinyVector<Complex_wp, 2>, General, 1, 1>&, TinyMatrix<TinyVector<Complex_wp, 2>, General, 1, 1>&); 
  
  SELDON_EXTERN template void AdvectionEquation_Base<Complex_wp, Dimension2>::GetTensorMass(const EllipticProblem<HarmonicAdvectionEquationDG<Dimension2> >&, int, int, const GlobalGenericMatrix<Complex_wp>&, int, TinyMatrix<Complex_wp, General, 1, 1>&);

  SELDON_EXTERN template void AdvectionEquation_Base<Complex_wp, Dimension2>::GetGradPhiTensor(const EllipticProblem<HarmonicAdvectionEquation<Dimension2> >&, int, int, const GlobalGenericMatrix<Complex_wp>&, int, TinyMatrix<TinyVector<Complex_wp, 2>, General, 1, 1>&, TinyMatrix<TinyVector<Complex_wp, 2>, General, 1, 1>&); 
  
  SELDON_EXTERN template void AdvectionEquation_Base<Complex_wp, Dimension2>::GetTensorMass(const EllipticProblem<HarmonicAdvectionEquation<Dimension2> >&, int, int, const GlobalGenericMatrix<Complex_wp>&, int, TinyMatrix<Complex_wp, General, 1, 1>&);

  //SELDON_EXTERN template void DoubleAdvectionEquation_Base<Complex_wp, Dimension2>::GetGradPhiTensor(const EllipticProblem<HarmonicDoubleAdvectionEquationDG<Dimension2> >&, int, int, const GlobalGenericMatrix<Complex_wp>&, int, TinyMatrix<TinyVector<Complex_wp, 2>, General, 1, 1>&, TinyMatrix<TinyVector<Complex_wp, 2>, General, 1, 1>&); 
  
  //SELDON_EXTERN template void DoubleAdvectionEquation_Base<Complex_wp, Dimension2>::GetTensorMass(const EllipticProblem<HarmonicDoubleAdvectionEquationDG<Dimension2> >&, int, int, const GlobalGenericMatrix<Complex_wp>&, int, TinyMatrix<Complex_wp, General, 1, 1>&);

  //SELDON_EXTERN template void DoubleAdvectionEquation_Base<Complex_wp, Dimension2>::GetGradPhiTensor(const EllipticProblem<HarmonicDoubleAdvectionEquation<Dimension2> >&, int, int, const GlobalGenericMatrix<Complex_wp>&, int, TinyMatrix<TinyVector<Complex_wp, 2>, General, 1, 1>&, TinyMatrix<TinyVector<Complex_wp, 2>, General, 1, 1>&); 
  
  //SELDON_EXTERN template void DoubleAdvectionEquation_Base<Complex_wp, Dimension2>::GetTensorMass(const EllipticProblem<HarmonicDoubleAdvectionEquation<Dimension2> >&, int, int, const GlobalGenericMatrix<Complex_wp>&, int, TinyMatrix<Complex_wp, General, 1, 1>&);

}
