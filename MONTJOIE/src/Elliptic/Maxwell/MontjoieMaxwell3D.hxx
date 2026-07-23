#ifndef MONTJOIE_FILE_MONTJOIE_MAXWELL_3D_HXX

#include "MontjoieMaxwell3D_Header.hxx"
#include "MontjoieMaxwell3D_Inline.hxx"

#include "Harmonic/MontjoieHarmonic.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "Computation/ProdMatVectHcurl3D.cxx"

#include "Elliptic/Maxwell/3D/HarmonicMaxwell3D.cxx"
#include "Elliptic/Maxwell/3D/ComputationModeMaxwell.cxx"
#include "Computation/ProdMatVectHcurl2D.cxx"
#include "Computation/ElementaryMatrixHcurl2D.cxx"
#include "Elliptic/Maxwell/3D/DefineSourceHarmonicMaxwell3D.cxx"
#include "Elliptic/Maxwell/3D/StaticMaxwell3D.cxx"


#ifdef MONTJOIE_WITH_TRANSMISSION
#include "Elliptic/Maxwell/3D/TransmissionModelMaxwell3D.cxx"
#endif

#ifdef SELDON_WITH_PRECONDITIONING  
#include "Elliptic/Maxwell/3D/PrecondMaxwell.cxx"
#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_MAXWELL_3D_HXX
#endif
