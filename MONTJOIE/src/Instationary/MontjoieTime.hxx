#ifndef MONTJOIE_FILE_MONTJOIE_TIME_HXX

#include "Quadrature/MontjoieQuadrature.hxx"
#include "Solver/MontjoieSolver.hxx"

#include "MontjoieTimeHeader.hxx"
#include "MontjoieTimeInline.hxx"

#include "Mesh/MontjoieMesh.hxx"

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "Instationary/MontjoieTimeCompil.hxx"
#endif

#include "Output/CommonInputOutput.cxx"
#include "Solver/NewtonSolver.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
// time-schemes
#include "Instationary/TimeSchemes.cxx"
#include "Instationary/ExplicitTimeSchemes.cxx"
//#include "Instationary/LocalTimeSchemes.cxx"
#include "Instationary/ImplicitTimeSchemes.cxx"
#include "Instationary/ModifiedEquation.cxx"

#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_TIME_HXX
#endif
