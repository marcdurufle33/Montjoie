#ifndef MONTJOIE_FILE_MONTJOIE_MAXWELL_3D_INLINE_HXX

#include "Harmonic/MontjoieHarmonicInline.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtzInline.hxx"

#ifdef MONTJOIE_WITH_TRANSMISSION
#include "Elliptic/Maxwell/3D/TransmissionModelMaxwell3D_Inline.cxx"
#endif

#include "Elliptic/Maxwell/3D/DefineSourceHarmonicMaxwell3D_Inline.cxx"
#include "Elliptic/Maxwell/3D/ComputationModeMaxwellInline.cxx"

#include "Elliptic/Maxwell/3D/HarmonicMaxwell3D_Inline.cxx"
#include "Elliptic/Maxwell/3D/DefineSourceHarmonicMaxwell3D_Inline.cxx"
#include "Elliptic/Maxwell/3D/StaticMaxwell3D_Inline.cxx"

#ifdef SELDON_WITH_PRECONDITIONING  
#include "Elliptic/Maxwell/3D/PrecondMaxwellInline.cxx"
#endif

#define MONTJOIE_FILE_MONTJOIE_MAXWELL_3D_INLINE_HXX
#endif
