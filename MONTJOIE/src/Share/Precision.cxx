#ifndef MONTJOIE_FILE_PRECISION_CXX

#ifdef MONTJOIE_WITH_MPFR
#include "PrecisionMpfr.cxx"
#endif

#ifdef MONTJOIE_WITH_FLOAT128
#include "PrecisionQuad.cxx"
#endif

#define MONTJOIE_FILE_PRECISION_CXX
#endif
