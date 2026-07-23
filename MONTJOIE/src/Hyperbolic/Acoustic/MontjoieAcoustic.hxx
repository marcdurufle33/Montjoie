#ifndef MONTJOIE_FILE_MONTJOIE_ACOUSTIC_HXX

#include "MontjoieAcousticHeader.hxx"
#include "MontjoieAcousticInline.hxx"

#include "Instationary/MontjoieUnsteady.hxx"

// including Laplace equation
#include "Elliptic/Helmholtz/MontjoieLaplace.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Acoustic/TimeAcoustic.cxx"
#include "Hyperbolic/Acoustic/ProdMatVectAcoustic.cxx"

#ifdef MONTJOIE_WITH_TIME_REVERSAL
#include "Hyperbolic/Acoustic/TimeReversal.cxx"
#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_ACOUSTIC_HXX
#endif

