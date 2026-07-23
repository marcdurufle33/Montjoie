#ifndef MONTJOIE_FILE_MONTJOIE_SOLVER_HXX

#include "Share/MontjoieCommon.hxx"

#include "Solver/MontjoieSolverHeader.hxx"
#include "Solver/MontjoieSolverInline.hxx"

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "Solver/MontjoieSolverCompil.hxx"
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Solver/SolveSystem.cxx"
#include "Solver/NonLinearEquations.cxx"
#include "Solver/NonLinearLeastSquares.cxx"
#include "Solver/NewtonSolver.cxx"
#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_SOLVER_HXX
#endif
