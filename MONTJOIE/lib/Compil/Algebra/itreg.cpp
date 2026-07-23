#include "MontjoieFlag.hxx"

#include "Regularization/ItregHeader.hxx"
#include "Regularization/ItregInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "Regularization/ForwardOperator.cxx"
#include "Regularization/RegularizationMethod.cxx"
#include "Regularization/Landweber.cxx"
#include "Regularization/NewtonInner.cxx"
#include "Regularization/NewtonTikh.cxx"
#include "Regularization/IrgnmCg.cxx"
#include "Regularization/MinimizationMethod.cxx"

#endif

namespace itreg
{
  SELDON_EXTERN template class ForwardOperator<Real_wp>;
  
  SELDON_EXTERN template class NewtonCG_with_G_P<Real_wp>;

  SELDON_EXTERN template class IRGNM_CG<Real_wp>;
  SELDON_EXTERN template class IRGNM_CG_with_G<Real_wp>;
  SELDON_EXTERN template class IRGNM_CG_Reortho<Real_wp>;
  SELDON_EXTERN template class NewtonCG_with_G<Real_wp>;
  SELDON_EXTERN template class IRGNM_CG_Precond<Real_wp>;
  SELDON_EXTERN template class IRGNM_CG_Precond_with_G<Real_wp>;
  SELDON_EXTERN template class IRGNM_CG_Stefan1<Real_wp>;
  SELDON_EXTERN template class IRGNM_CG_Stefan2<Real_wp>;

  SELDON_EXTERN template class LandweberP<Real_wp>;
  SELDON_EXTERN template class ModLandweberP<Real_wp>;

  SELDON_EXTERN template class ScaledLandweber<Real_wp>;
  SELDON_EXTERN template class Landweber<Real_wp>;
  SELDON_EXTERN template class ModLandweber<Real_wp>;

  SELDON_EXTERN template class NewtonInnerP<Real_wp>;
  SELDON_EXTERN template class NewtonLwP<Real_wp>;
  SELDON_EXTERN template class NewtonNuP<Real_wp>;

  SELDON_EXTERN template class NewtonInner<Real_wp>;
  SELDON_EXTERN template class NewtonLandweber<Real_wp>;
  SELDON_EXTERN template class NewtonNuMethod<Real_wp>;
  SELDON_EXTERN template class NewtonCG<Real_wp>;
  SELDON_EXTERN template class NewtonCGReOrtho<Real_wp>;
  SELDON_EXTERN template class NewtonCGPrecond<Real_wp>;
  SELDON_EXTERN template class SvdOrtho<Real_wp>;

  SELDON_EXTERN template class IRGNMP<Real_wp>;
  SELDON_EXTERN template class IRGNM_with_G_P<Real_wp>;
  SELDON_EXTERN template class LM2P<Real_wp>;
  SELDON_EXTERN template class IterTikhP<Real_wp>;
  SELDON_EXTERN template class NLTikhP<Real_wp>;

  SELDON_EXTERN template class NewtonTikh<Real_wp>;
  SELDON_EXTERN template class IRGNM<Real_wp>;
  SELDON_EXTERN template class IRGNM_with_G<Real_wp>;
  SELDON_EXTERN template class LevMarquardt<Real_wp>;
  SELDON_EXTERN template class LM2<Real_wp>;
  SELDON_EXTERN template class IterTikh<Real_wp>;
  SELDON_EXTERN template class NLTikh<Real_wp>;
  SELDON_EXTERN template class IRGNMFrozen<Real_wp>;
  SELDON_EXTERN template class LevMarquardtFrozen<Real_wp>;
  SELDON_EXTERN template class LevMarquardtGenuine<Real_wp>;

  SELDON_EXTERN template class RegMethodP<Real_wp>;
  SELDON_EXTERN template class RegularizationMethod<Real_wp>;
  
  SELDON_EXTERN template class GslRegMethodP<Real_wp>;
  SELDON_EXTERN template class MinimizationGslReg<Real_wp>;
  
}
