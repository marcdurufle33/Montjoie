#ifndef MONTJOIE_FILE_MONTJOIE_FINITE_ELEMENT_HEADER_HXX

#include "Mesh/MontjoieMeshHeader.hxx"

// definition of effect of transformation Fi for H1 elements
#include "FiniteElement/FiniteElementH1.hxx"
// base class for all finite elements
#include "FiniteElement/ElementReference.hxx"

#ifdef MONTJOIE_WITH_NODAL_HCURL
#include "FiniteElement/FiniteElementHcurl2D.hxx"
#ifdef MONTJOIE_WITH_THREE_DIM
#include "FiniteElement/FiniteElementHcurl3D.hxx"
#endif
#endif

#ifdef MONTJOIE_WITH_NODAL_HDIV
#include "FiniteElement/FiniteElementHdiv.hxx"
#endif

// 2-D nodal finite elements
#include "FiniteElement/FaceReference.hxx"
#include "FiniteElement/Quadrangle/QuadrangleReference.hxx"
#include "FiniteElement/Triangle/TriangleReference.hxx"

#ifdef MONTJOIE_WITH_THREE_DIM
// base class for 3-D finite element
#include "FiniteElement/VolumeReference.hxx"
// base class for hexahedron finite elements
#include "FiniteElement/Hexahedron/HexahedronReference.hxx"
// base class for pyramidal finite elements
#include "FiniteElement/Pyramid/PyramidReference.hxx"
 
// alternative class for pyramids using Knabbner transformation (pyramid split into two tetrahedra)
// #include "FiniteElement/Pyramid/PyramidReferenceSplit.hxx"

// base class for wedge finite elements
#include "FiniteElement/Wedge/WedgeReference.hxx"

// base class for tetrahedronl finite elements
#include "FiniteElement/Tetrahedron/TetrahedronReference.hxx"
#endif

/***********************
 * 2-D Finite Elements *
 ***********************/

#ifdef MONTJOIE_WITH_NODAL_H1
// Quadrilateral finite elements with Gauss-Legendre points for quadrature (Gauss-Lobatto for dofs)
#include "FiniteElement/Quadrangle/QuadrangleGauss.hxx"
// Quadrilateral Gauss-Lobatto nodal finite element 
#include "FiniteElement/Quadrangle/QuadrangleLobatto.hxx"
#include "FiniteElement/Quadrangle/QuadrangleRadau.hxx"
// classical Pk finite element (can be either regular points or hesthaven's like points)
#include "FiniteElement/Triangle/TriangleClassical.hxx"
// quasi-lumped triangles (sparse factorization of mass matrix)
#include "FiniteElement/Triangle/TriangleQuasiLumped.hxx"
#endif

#ifdef MONTJOIE_WITH_HP_H1
// hp finite elements for triangles and quadrangles
#include "FiniteElement/Quadrangle/QuadrangleHierarchic.hxx"
#include "FiniteElement/Triangle/TriangleHierarchic.hxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_DG
// Quadrilateral finite elements with Gauss-Legendre points for dofs and quadrature
#include "FiniteElement/Quadrangle/QuadrangleDgGauss.hxx"
// mass-lumped triangle with lobatto points on edges
#include "FiniteElement/Triangle/TriangleDgMassLumped.hxx"
#endif

#ifdef MONTJOIE_WITH_ORTHO_DG
// orthogonal functions on triangle 
#include "FiniteElement/Triangle/TriangleDgOrtho.hxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HCURL
// edge finite element Rk for first family of Nedelec on triangles (see Graglia)
#include "FiniteElement/Triangle/TriangleHcurlFirstFamily.hxx"
// edge finite element Pk for second family, with nodal points
#include "FiniteElement/Triangle/TriangleHcurlSecondFamily.hxx"
// first family of nedelec on quadrilaterals
#include "FiniteElement/Quadrangle/QuadrangleHcurlFirstFamily.hxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlGaussFirstFamily.hxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlOptimalFirstFamily.hxx"
// second family of nedelec on quadrilaterals
#include "FiniteElement/Quadrangle/QuadrangleHcurlLobatto.hxx"
#endif

#ifdef MONTJOIE_WITH_HP_HCURL
// optimal hp edge finite element Rk for first family of Nedelec on quadrangles
#include "FiniteElement/Quadrangle/QuadrangleHcurlOptimalHpFirstFamily.hxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlHpFirstFamily.hxx"
// edge finite element Pk for second family of Nedelec on triangles (see Solin)
#include "FiniteElement/Triangle/TriangleHcurlHierarchic.hxx"
// optimal hp edge finite element Rk for first family of Nedelec on triangles
#include "FiniteElement/Triangle/TriangleHcurlOptimalHpFirstFamily.hxx"
#endif


#ifdef MONTJOIE_WITH_NODAL_HDIV
// first family of nedelec on quadrilaterals
#include "FiniteElement/Quadrangle/QuadrangleHdivFirstFamily.hxx"
#include "FiniteElement/Quadrangle/QuadrangleHdivOptimalFirstFamily.hxx"
// first family of nedelec on triangles
#include "FiniteElement/Triangle/TriangleHdivFirstFamily.hxx"
#endif

#ifdef MONTJOIE_WITH_HP_HDIV
// optimal hp edge finite element Rk for first family of Nedelec on quadrangles
#include "FiniteElement/Quadrangle/QuadrangleHdivOptimalHpFirstFamily.hxx"
#include "FiniteElement/Quadrangle/QuadrangleHdivHpFirstFamily.hxx"
// optimal hp edge finite element Rk for first family of Nedelec on triangles
#include "FiniteElement/Triangle/TriangleHdivOptimalHpFirstFamily.hxx"
#endif

/***********************
 * 3-D Finite Elements *
 ***********************/


#ifdef MONTJOIE_WITH_THREE_DIM

#ifdef MONTJOIE_WITH_NODAL_H1
// nodal finite elements in 3-D for hexahedra, pyramids, prisms and tetrahedra
#include "FiniteElement/Hexahedron/HexahedronGauss.hxx"
#include "FiniteElement/Hexahedron/HexahedronLobatto.hxx"
#include "FiniteElement/Pyramid/PyramidClassical.hxx"
// #include "FiniteElement/Pyramid/PyramidClassicalSplit.hxx"
#include "FiniteElement/Wedge/WedgeClassical.hxx"
#include "FiniteElement/Tetrahedron/TetrahedronClassical.hxx"
// #include "FiniteElement/Tetrahedron/TetrahedronMassLumped.hxx"
// #include "FiniteElement/Tetrahedron/TetrahedronQuasiLumped.hxx"
#endif

#ifdef MONTJOIE_WITH_HP_H1
// hp finite elements in 3-D by using tensorized points (doesn't work on general meshes)
#include "FiniteElement/Hexahedron/HexahedronHierarchic.hxx"
#include "FiniteElement/Pyramid/PyramidHierarchic.hxx"
#include "FiniteElement/Wedge/WedgeHierarchic.hxx"
#include "FiniteElement/Tetrahedron/TetrahedronHierarchic.hxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HCURL
// first family of nedelec on hexahedron
#include "FiniteElement/Hexahedron/HexahedronHcurlFirstFamily.hxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlOptimalFirstFamily.hxx"
// second family of nedelec on hexahedron
#include "FiniteElement/Hexahedron/HexahedronHcurlLobatto.hxx"
// first family of nedelec on wedges
#include "FiniteElement/Wedge/WedgeHcurlFirstFamily.hxx"
#include "FiniteElement/Wedge/WedgeHcurlOptimalFirstFamily.hxx"
// first family of nedelec on pyramids
#include "FiniteElement/Pyramid/PyramidHcurlFirstFamily.hxx"
#include "FiniteElement/Pyramid/PyramidHcurlOptimalFirstFamily.hxx"
// first family of nedelec on tets
#include "FiniteElement/Tetrahedron/TetrahedronHcurlFirstFamily.hxx"
// second family of nedelec on tets
//#include "FiniteElement/Tetrahedron/TetrahedronHcurlSecondFamily.hxx"
#include "FiniteElement/Tetrahedron/TetrahedronHcurlLobatto.hxx"
#endif

#ifdef MONTJOIE_WITH_HP_HCURL
// optimal first family of nedelec on hexahedra
#include "FiniteElement/Hexahedron/HexahedronHcurlOptimalHpFirstFamily.hxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlHpFirstFamily.hxx"
// first family of nedelec on pyramids
#include "FiniteElement/Pyramid/PyramidHcurlOptimalHpFirstFamily.hxx"
#include "FiniteElement/Pyramid/PyramidHcurlHpFirstFamily.hxx"
// first family of nedelec on wedges
#include "FiniteElement/Wedge/WedgeHcurlOptimalHpFirstFamily.hxx"
#include "FiniteElement/Wedge/WedgeHcurlHpFirstFamily.hxx"
// optimal hp first family of nedelec on tets
#include "FiniteElement/Tetrahedron/TetrahedronHcurlOptimalHpFirstFamily.hxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_DG
// nodal finite elements for DG
#include "FiniteElement/Hexahedron/HexahedronDgGauss.hxx"
#include "FiniteElement/Wedge/WedgeDgClassical.hxx"
// #include "FiniteElement/Hexahedron/HexahedronDgJacobi.hxx"
#endif

#ifdef MONTJOIE_WITH_ORTHO_DG
// orthogonal functions for DG
#include "FiniteElement/Hexahedron/HexahedronDgGauss.hxx"
#include "FiniteElement/Pyramid/PyramidDgOrtho.hxx"
#include "FiniteElement/Wedge/WedgeDgOrtho.hxx"
#include "FiniteElement/Tetrahedron/TetrahedronDgOrtho.hxx"
#endif

#ifdef MONTJOIE_WITH_LEGENDRE_DG
// use of P_r for each element (even hexahedron, pyramid and wedge) with orthogonal functions
#include "FiniteElement/Hexahedron/HexahedronDgLegendre.hxx"
#include "FiniteElement/Pyramid/PyramidDgLegendre.hxx"
#include "FiniteElement/Wedge/WedgeDgLegendre.hxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HDIV
// first family of nedelec on hexahedron
#include "FiniteElement/Hexahedron/HexahedronHdivFirstFamily.hxx"
#include "FiniteElement/Hexahedron/HexahedronHdivOptimalFirstFamily.hxx"
// first family of nedelec on wedges
#include "FiniteElement/Wedge/WedgeHdivFirstFamily.hxx"
#include "FiniteElement/Wedge/WedgeHdivOptimalFirstFamily.hxx"
// first family of nedelec on pyramids
#include "FiniteElement/Pyramid/PyramidHdivFirstFamily.hxx"
#include "FiniteElement/Pyramid/PyramidHdivOptimalFirstFamily.hxx"
// first family of nedelec on tets
#include "FiniteElement/Tetrahedron/TetrahedronHdivFirstFamily.hxx"
#endif

#ifdef MONTJOIE_WITH_HP_HDIV
#include "FiniteElement/Triangle/TriangleDgOrtho.hxx"
#include "FiniteElement/Quadrangle/QuadrangleDgOrtho.hxx"
// optimal first family of nedelec on hexahedra
#include "FiniteElement/Hexahedron/HexahedronHdivOptimalHpFirstFamily.hxx"
#include "FiniteElement/Hexahedron/HexahedronHdivHpFirstFamily.hxx"
// first family of nedelec on pyramids
#include "FiniteElement/Pyramid/PyramidHdivOptimalHpFirstFamily.hxx"
#include "FiniteElement/Pyramid/PyramidHdivHpFirstFamily.hxx"
// first family of nedelec on wedges
#include "FiniteElement/Wedge/WedgeHdivOptimalHpFirstFamily.hxx"
#include "FiniteElement/Wedge/WedgeHdivHpFirstFamily.hxx"
// optimal hp first family of nedelec on tets
#include "FiniteElement/Tetrahedron/TetrahedronHdivOptimalHpFirstFamily.hxx"
#endif
 
#endif

// class storing various finite elements
#include "Harmonic/VarFiniteElement.hxx"


#define MONTJOIE_FILE_MONTJOIE_FINITE_ELEMENT_HEADER_HXX
#endif

