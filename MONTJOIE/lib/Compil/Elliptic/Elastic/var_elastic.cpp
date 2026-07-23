#include "MontjoieFlag.hxx"

#include "Elliptic/Elastic/MontjoieElasticHeader.hxx"
#include "Elliptic/Elastic/MontjoieElasticInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Elastic/VarElastic.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"
#include "Computation/ProdMatVectH1.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarElastic_Dim<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarElastic_Dim<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VarElastic_Dim<Real_wp, Dimension3>;
  SELDON_EXTERN template class VarElastic_Dim<Complex_wp, Dimension3>;

  SELDON_EXTERN template class VarElastic_Eq<ElasticEquation<Dimension2> >;
  SELDON_EXTERN template class VarElastic_Eq<ElasticEquation<Dimension3> >;
  SELDON_EXTERN template class VarElastic_Eq<HarmonicElasticEquation<Dimension2> >;
  SELDON_EXTERN template class VarElastic_Eq<HarmonicElasticEquation<Dimension3> >;

  SELDON_EXTERN template class VarElastic_Eq<ElasticEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarElastic_Eq<ElasticEquationDG<Dimension3> >;
  SELDON_EXTERN template class VarElastic_Eq<HarmonicElasticEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarElastic_Eq<HarmonicElasticEquationDG<Dimension3> >;

  SELDON_EXTERN template class VarElastic_Cplx<HarmonicElasticEquation<Dimension2> >;
  SELDON_EXTERN template class VarElastic_Cplx<HarmonicElasticEquation<Dimension3> >;
  SELDON_EXTERN template class VarElastic_Cplx<HarmonicElasticEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarElastic_Cplx<HarmonicElasticEquationDG<Dimension3> >;

  SELDON_EXTERN template class VarElastic_Real<ElasticEquation<Dimension2> >;
  SELDON_EXTERN template class VarElastic_Real<ElasticEquation<Dimension3> >;
  SELDON_EXTERN template class VarElastic_Real<ElasticEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarElastic_Real<ElasticEquationDG<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<ElasticEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<ElasticEquation<Dimension3> >;
  SELDON_EXTERN template class EllipticProblem<HarmonicElasticEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HarmonicElasticEquation<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<ElasticEquationDG<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<ElasticEquationDG<Dimension3> >;
  SELDON_EXTERN template class EllipticProblem<HarmonicElasticEquationDG<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HarmonicElasticEquationDG<Dimension3> >;

  SELDON_EXTERN template class ImpedanceABC<Complex_wp, HarmonicElasticEquation<Dimension2> >;
  SELDON_EXTERN template class ImpedanceABC<Real_wp, ElasticEquation<Dimension2> >;
  SELDON_EXTERN template class ImpedanceABC<Complex_wp, ElasticEquation<Dimension2> >;

  SELDON_EXTERN template class ImpedanceABC<Complex_wp, HarmonicElasticEquation<Dimension3> >;
  SELDON_EXTERN template class ImpedanceABC<Real_wp, ElasticEquation<Dimension3> >;
  SELDON_EXTERN template class ImpedanceABC<Complex_wp, ElasticEquation<Dimension3> >;

}
