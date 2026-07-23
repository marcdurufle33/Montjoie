#ifndef MONTJOIE_FILE_MONTJOIE_HELMHOLTZ_HXX

#include "MontjoieHelmholtzHeader.hxx"
#include "MontjoieHelmholtzInline.hxx"

#include "Harmonic/MontjoieHarmonic.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY

// sources
#include "Elliptic/Helmholtz/DefineSourceHelmholtz.cxx"

// general functions to solve Laplace equations
#include "Elliptic/Helmholtz/VarLaplace.cxx"

#include "Elliptic/Helmholtz/LaplacianModalSource.cxx"

#ifdef MONTJOIE_WITH_THIN_SLOT_MODEL
// thin slots models for 2-D computations
#include "Elliptic/Helmholtz/ThinSlotHelmholtzModel.cxx"
#endif

#ifdef MONTJOIE_WITH_TRANSMISSION
#include "Elliptic/Helmholtz/TransmissionModelHelmholtz.cxx"
#endif

#ifdef MONTJOIE_WITH_INTEGRAL_EQUATION
// integral equations in 3-D (Galerkin, no multipole, very basic)
#include "Elliptic/Helmholtz/IntegralEquationHelmholtz.cxx"
#endif

// general functions to solve Helmholtz equations
#include "Elliptic/Helmholtz/VarHelmholtz.cxx"
#include "Elliptic/Helmholtz/ImpedanceHelmholtz.cxx"


#ifdef MONTJOIE_WITH_ONE_DIM
// 1-D pathologic case
#include "Elliptic/Helmholtz/Helmholtz1D.cxx"
#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_HELMHOLTZ_HXX
#endif
