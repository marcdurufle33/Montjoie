#ifndef MONTJOIE_FILE_MONTJOIE_UNSTEADY_HXX

#include "MontjoieUnsteadyHeader.hxx"
#include "MontjoieUnsteadyInline.hxx"

#include "Harmonic/MontjoieHarmonic.hxx"
#include "MontjoieTime.hxx"

/***********************************************
 * Generic classes for time-dependent problems *
 ***********************************************/

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "Instationary/MontjoieUnsteadyCompil.hxx"
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "Source/TimeSource.cxx"

#include "Instationary/ContinuousMassMatrix.cxx"
#include "Instationary/DiscontinuousMassMatrix.cxx"
#include "Instationary/ContinuousUnsteadyMassMatrix.cxx"
#include "Instationary/DiscontinuousUnsteadyMassMatrix.cxx"

#include "Instationary/EvaluationCflWave.cxx"
#include "Instationary/VarInstationary.cxx"
#include "Instationary/LocalTimeSchemes.cxx"

#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_UNSTEADY_HXX
#endif
