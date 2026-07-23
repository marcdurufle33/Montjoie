#ifndef ITREG_FILE_ITREG_HEADER_HXX

#include "Share/MontjoieCommonHeader.hxx"

namespace itreg
{
  using namespace Seldon;
}

#include "Solver/NonLinearLeastSquares.hxx"
#include "Solver/NonLinearEquations.hxx"

#include "ForwardOperator.hxx"
#include "RegularizationMethod.hxx"
#include "Landweber.hxx"
#include "NewtonInner.hxx"
#include "NewtonTikh.hxx"
#include "IrgnmCg.hxx"
#include "MinimizationMethod.hxx"

#define ITREG_FILE_ITREG_HEADER_HXX
#endif
