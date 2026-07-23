#ifndef MONTJOIE_FILE_MONTJOIE_HELMHOLTZ_HEADER_HXX

#include "Harmonic/MontjoieHarmonicHeader.hxx"

#ifdef MONTJOIE_WITH_ONE_DIM
// 1-D pathologic case
#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
#endif

#ifdef MONTJOIE_WITH_THIN_SLOT_MODEL
// thin slots models for 2-D computations
#include "Elliptic/Helmholtz/ThinSlotHelmholtzModel.hxx"
#endif

#ifdef MONTJOIE_WITH_INTEGRAL_EQUATION
// integral equations in 3-D (Galerkin, no multipole, very basic)
#include "Elliptic/Helmholtz/IntegralEquationHelmholtz.hxx"
#endif

// general functions to solve Helmholtz equations
#include "Elliptic/Helmholtz/VarHelmholtz.hxx"
#include "Elliptic/Helmholtz/ProdMatVectHelmholtz.hxx"
#include "Elliptic/Helmholtz/ImpedanceHelmholtz.hxx"

#ifdef MONTJOIE_WITH_TRANSMISSION
#include "Elliptic/Helmholtz/TransmissionModelHelmholtz.hxx"
#endif

// modal source
#include "Elliptic/Helmholtz/LaplacianModalSource.hxx"

// sources
#include "Elliptic/Helmholtz/DefineSourceHelmholtz.hxx"

// declaration of classes used to solve Laplace equation
#include "Elliptic/Helmholtz/VarLaplace.hxx"

#define MONTJOIE_FILE_MONTJOIE_HELMHOLTZ_HEADER_HXX
#endif
