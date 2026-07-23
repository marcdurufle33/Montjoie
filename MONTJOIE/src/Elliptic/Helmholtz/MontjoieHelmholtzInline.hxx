#ifndef MONTJOIE_FILE_MONTJOIE_HELMHOLTZ_INLINE_HXX

#include "Harmonic/MontjoieHarmonicInline.hxx"

// sources
#include "Elliptic/Helmholtz/DefineSourceHelmholtzInline.cxx"
#include "Elliptic/Helmholtz/LaplacianModalSourceInline.cxx"

#include "Elliptic/Helmholtz/VarLaplaceInline.cxx"

#ifdef MONTJOIE_WITH_TRANSMISSION
#include "Elliptic/Helmholtz/TransmissionModelHelmholtzInline.cxx"
#endif

#include "Elliptic/Helmholtz/VarHelmholtzInline.cxx"
#include "Elliptic/Helmholtz/ProdMatVectHelmholtzInline.cxx"
#include "Elliptic/Helmholtz/ImpedanceHelmholtzInline.cxx"


#ifdef MONTJOIE_WITH_ONE_DIM
#include "Elliptic/Helmholtz/Helmholtz1D_Inline.cxx"
#endif

#define MONTJOIE_FILE_MONTJOIE_HELMHOLTZ_INLINE_HXX
#endif
