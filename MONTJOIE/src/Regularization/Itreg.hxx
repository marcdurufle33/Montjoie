#ifndef ITREG_FILE_ITREG_HXX

#include "ItregHeader.hxx"
#include "ItregInline.hxx"

#include "Share/MontjoieCommon.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "ForwardOperator.cxx"
#include "RegularizationMethod.cxx"
#include "Landweber.cxx"
#include "NewtonInner.cxx"
#include "NewtonTikh.cxx"
#include "IrgnmCg.cxx"

#include "Solver/NonLinearLeastSquares.cxx"
#include "Solver/NonLinearEquations.cxx"

#include "MinimizationMethod.cxx"

#endif


#define ITREG_FILE_ITREG_HXX
#endif
