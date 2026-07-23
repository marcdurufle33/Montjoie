#ifndef MONTJOIE_FILE_MONTJOIE_LAPLACE_HXX

#include "Elliptic/Helmholtz/MontjoieLaplaceHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieLaplaceInline.hxx"

// including classes for Helmholtz equation
#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
// general functions to solve Laplace equations
#include "Elliptic/Helmholtz/VarLaplace.cxx"
#endif

#define MONTJOIE_FILE_MONTJOIE_LAPLACE_HXX
#endif
