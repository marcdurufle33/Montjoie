#ifndef MONTJOIE_FILE_MONTJOIE_QUADRATURE_HXX

#include "Quadrature/MontjoieQuadratureHeader.hxx"
#include "Quadrature/MontjoieQuadratureInline.hxx"

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#include "Share/MontjoieCommon.hxx"

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "Quadrature/MontjoieQuadratureCompil.hxx"
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "Quadrature/GaussLobattoPoints.cxx"
#include "Quadrature/GaussLobatto.cxx"
#include "Quadrature/GaussJacobi.cxx"

#include "Quadrature/TriangleQuadrature.cxx"
#include "Quadrature/QuadrangleQuadrature.cxx"

#ifdef MONTJOIE_WITH_THREE_DIM
#include "Quadrature/TetrahedronQuadrature.cxx"
#include "Quadrature/PyramidQuadrature.cxx"
#include "Quadrature/HexahedronQuadrature.cxx"
#endif

#endif
#endif

#define MONTJOIE_FILE_MONTJOIE_QUADRATURE_HXX
#endif
