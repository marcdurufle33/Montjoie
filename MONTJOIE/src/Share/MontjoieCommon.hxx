#ifndef MONTJOIE_FILE_MONTJOIE_COMMON_HXX

#include "MontjoieCommonHeader.hxx"
#include "MontjoieCommonInline.hxx"

// including Seldon and other Seldon-like methods
#include "Algebra/MontjoieAlgebra.hxx"

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "MontjoieCommonCompil.hxx"
#endif

// Files in directory Share : 

#ifndef SELDON_WITH_COMPILED_LIBRARY

// FFT stuff
#include "Share/FFT.cxx"

// other useful functions
#include "Share/CommonMontjoie.cxx"
#include "Share/Timer.cxx"
#include "Share/RandomGenerator.cxx"

// basic polynoms classes (univariate -> one variable x
//  multivariate -> any numbers of variables (for example x,y,z))
#include "Share/UnivariatePolynomial.cxx"
#include "Share/MultivariatePolynomial.cxx"

// interface for Bessel functions
#include "Share/BesselFunctionsInterface.cxx"
#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_COMMON_HXX
#endif

