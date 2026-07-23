#ifndef MONTJOIE_FILE_MONTJOIE_QUADRATURE_HEADER_HXX

#include "Share/MontjoieCommonHeader.hxx"

// classes for computing Gauss-Lobatto points
#include "Quadrature/GaussLobattoPoints.hxx"
#include "Quadrature/GaussLobatto.hxx"
// classes for computing Gauss-Jacobi points
#include "Quadrature/GaussJacobi.hxx"

#include "Quadrature/TriangleQuadrature.hxx"
#include "Quadrature/QuadrangleQuadrature.hxx"

#ifdef MONTJOIE_WITH_THREE_DIM
#include "Quadrature/TetrahedronQuadrature.hxx"
#include "Quadrature/PyramidQuadrature.hxx"
#include "Quadrature/HexahedronQuadrature.hxx"
#endif

#define MONTJOIE_FILE_MONTJOIE_QUADRATURE_HEADER_HXX
#endif
