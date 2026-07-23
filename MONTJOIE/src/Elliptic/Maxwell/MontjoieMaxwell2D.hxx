#ifndef MONTJOIE_FILE_MONTJOIE_MAXWELL_2D_HXX

#include "MontjoieMaxwell2D_Header.hxx"
#include "MontjoieMaxwell2D_Inline.hxx"

#include "Harmonic/MontjoieHarmonic.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/2D/HarmonicMaxwell2D.cxx"
#include "Elliptic/Maxwell/2D/StaticMaxwell2D.cxx"
#include "Elliptic/Maxwell/2D/DefineSourceHarmonicMaxwell2D.cxx"

#include "Computation/ElementaryMatrixHcurl2D.cxx"
#include "Computation/ProdMatVectHcurl2D.cxx"
#endif

#define MONTJOIE_FILE_MONTJOIE_MAXWELL_2D_HXX
#endif
