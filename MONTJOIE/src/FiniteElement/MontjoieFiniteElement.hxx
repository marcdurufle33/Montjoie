#ifndef MONTJOIE_FILE_MONTJOIE_FINITE_ELEMENT_HXX

#include "MontjoieFiniteElementHeader.hxx"
#include "MontjoieFiniteElementInline.hxx"

#include "Mesh/MontjoieMesh.hxx"

/***********************
 * 2-D Finite Elements *
 ***********************/

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/MontjoieFiniteElementCompil.hxx"
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "FiniteElement/FiniteElementH1.cxx"
#include "FiniteElement/ElementReference.cxx"


// 2-D nodal finite elements
#include "FiniteElement/FaceReference.cxx"
#include "FiniteElement/Triangle/TriangleReference.cxx"
#include "FiniteElement/Quadrangle/QuadrangleReference.cxx"

#ifdef MONTJOIE_WITH_NODAL_H1
// Lobatto finite element
#include "FiniteElement/Quadrangle/QuadrangleGauss.cxx"
#include "FiniteElement/Quadrangle/QuadrangleLobatto.cxx"
#include "FiniteElement/Quadrangle/QuadrangleRadau.cxx"
// classical Pk finite element (can be either regular points or hesthaven's like points)
#include "FiniteElement/Triangle/TriangleClassical.cxx"
#include "FiniteElement/Triangle/TriangleQuasiLumped.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_H1
#include "FiniteElement/Quadrangle/QuadrangleHierarchic.cxx"
#include "FiniteElement/Triangle/TriangleHierarchic.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_DG
// mass-lumped triangle with lobatto points on edges
#include "FiniteElement/Triangle/TriangleDgMassLumped.cxx"
// gauss case
#include "FiniteElement/Quadrangle/QuadrangleDgGauss.cxx"
#endif

#ifdef MONTJOIE_WITH_ORTHO_DG
#include "FiniteElement/Triangle/TriangleDgOrtho.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HCURL
#include "FiniteElement/FiniteElementHcurl2D.cxx"
// edge finite element Rk for first family of Nedelec on triangles (see Graglia)
#include "FiniteElement/Triangle/TriangleHcurlFirstFamily.cxx"
// edge finite element Pk for second family, with nodal points
#include "FiniteElement/Triangle/TriangleHcurlSecondFamily.cxx"
// first family of nedelec on quadrilaterals
#include "FiniteElement/Quadrangle/QuadrangleHcurlFirstFamily.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlGaussFirstFamily.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlOptimalFirstFamily.cxx"
// second family of nedelec on quadrilaterals
#include "FiniteElement/Quadrangle/QuadrangleHcurlLobatto.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_HCURL
// optimal hp edge finite element Rk for first family of Nedelec on quadrangles
#include "FiniteElement/Quadrangle/QuadrangleHcurlOptimalHpFirstFamily.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlHpFirstFamily.cxx"
// edge finite element Pk for second family of Nedelec on triangles (see Solin)
#include "FiniteElement/Triangle/TriangleHcurlHierarchic.cxx"
// optimal hp edge finite element Rk for first family of Nedelec on triangles
#include "FiniteElement/Triangle/TriangleHcurlOptimalHpFirstFamily.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HDIV
#include "FiniteElement/FiniteElementHdiv.cxx"
// first family of nedelec on quadrilaterals
#include "FiniteElement/Quadrangle/QuadrangleHdivFirstFamily.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHdivOptimalFirstFamily.cxx"
// first family of nedelec on triangles
#include "FiniteElement/Triangle/TriangleHdivFirstFamily.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_HDIV
// optimal hp edge finite element Rk for first family of Nedelec on quadrangles
#include "FiniteElement/Quadrangle/QuadrangleHdivOptimalHpFirstFamily.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHdivHpFirstFamily.cxx"
// optimal hp edge finite element Rk for first family of Nedelec on triangles
#include "FiniteElement/Triangle/TriangleHdivOptimalHpFirstFamily.cxx"
#endif


/***********************
 * 3-D Finite Elements *
 ***********************/


#ifdef MONTJOIE_WITH_THREE_DIM

// base class for 3-D finite element
#include "FiniteElement/VolumeReference.cxx"
// base class for tetrahedronl finite elements
#include "FiniteElement/Tetrahedron/TetrahedronReference.cxx"
// base class for hexahedronl finite elements
#include "FiniteElement/Hexahedron/HexahedronReference.cxx"
// base class for pyramidal finite elements
#include "FiniteElement/Pyramid/PyramidReference.cxx"
//#include "FiniteElement/Pyramid/PyramidReferenceSplit.cxx"
// base class for wedge finite elements
#include "FiniteElement/Wedge/WedgeReference.cxx"


#ifdef MONTJOIE_WITH_NODAL_H1
// nodal Lobatto hexahedronl elements
#include "FiniteElement/Hexahedron/HexahedronGauss.cxx"
#include "FiniteElement/Hexahedron/HexahedronLobatto.cxx"
// classical Pk finite element on tets
#include "FiniteElement/Tetrahedron/TetrahedronClassical.cxx"
#include "FiniteElement/Pyramid/PyramidClassical.cxx"
//#include "FiniteElement/Pyramid/PyramidClassicalSplit.cxx"
#include "FiniteElement/Wedge/WedgeClassical.cxx"
//#include "FiniteElement/Tetrahedron/TetrahedronMassLumped.cxx"
//#include "FiniteElement/Tetrahedron/TetrahedronQuasiLumped.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_H1
#include "FiniteElement/Hexahedron/HexahedronHierarchic.cxx"
#include "FiniteElement/Tetrahedron/TetrahedronHierarchic.cxx"
#include "FiniteElement/Pyramid/PyramidHierarchic.cxx"
#include "FiniteElement/Wedge/WedgeHierarchic.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HCURL
#include "FiniteElement/FiniteElementHcurl3D.cxx"
// first family of nedelec on hexahedra
#include "FiniteElement/Hexahedron/HexahedronHcurlFirstFamily.cxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlOptimalFirstFamily.cxx"
// first family of nedelec on hexahedra
#include "FiniteElement/Hexahedron/HexahedronHcurlLobatto.cxx"
// first family of nedelec on wedges
#include "FiniteElement/Wedge/WedgeHcurlFirstFamily.cxx"
#include "FiniteElement/Wedge/WedgeHcurlOptimalFirstFamily.cxx"
// first family of nedelec on pyramids
#include "FiniteElement/Pyramid/PyramidHcurlFirstFamily.cxx"
#include "FiniteElement/Pyramid/PyramidHcurlOptimalFirstFamily.cxx"
// first family of nedelec on tets
#include "FiniteElement/Tetrahedron/TetrahedronHcurlFirstFamily.cxx"
// second family of nedelec on tets
//#include "FiniteElement/Tetrahedron/TetrahedronHcurlSecondFamily.cxx"
#include "FiniteElement/Tetrahedron/TetrahedronHcurlLobatto.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_HCURL
// optimal first family of nedelec on hexahedronls
#include "FiniteElement/Hexahedron/HexahedronHcurlOptimalHpFirstFamily.cxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlHpFirstFamily.cxx"
// first family of nedelec on pyramids
#include "FiniteElement/Pyramid/PyramidHcurlOptimalHpFirstFamily.cxx"
#include "FiniteElement/Pyramid/PyramidHcurlHpFirstFamily.cxx"
// first family of nedelec on wedges
#include "FiniteElement/Wedge/WedgeHcurlOptimalHpFirstFamily.cxx"
#include "FiniteElement/Wedge/WedgeHcurlHpFirstFamily.cxx"
// optimal hp first family of nedelec on tets
#include "FiniteElement/Tetrahedron/TetrahedronHcurlOptimalHpFirstFamily.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_DG
// Gauss-Lobatto points
#include "FiniteElement/Hexahedron/HexahedronDgGauss.cxx"
// Jacobi points
// #include "FiniteElement/Hexahedron/HexahedronDgJacobi.cxx"
// DG classical Pk finite element on tets
#include "FiniteElement/Wedge/WedgeDgClassical.cxx"
#endif

#ifdef MONTJOIE_WITH_ORTHO_DG
#include "FiniteElement/Hexahedron/HexahedronDgGauss.cxx"
#include "FiniteElement/Tetrahedron/TetrahedronDgOrtho.cxx"
#include "FiniteElement/Pyramid/PyramidDgOrtho.cxx"
#include "FiniteElement/Wedge/WedgeDgOrtho.cxx"
#endif

#ifdef MONTJOIE_WITH_LEGENDRE_DG
#include "FiniteElement/Hexahedron/HexahedronDgLegendre.cxx"
#include "FiniteElement/Pyramid/PyramidDgLegendre.cxx"
#include "FiniteElement/Wedge/WedgeDgLegendre.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HDIV
// first family of nedelec on hexahedron
#include "FiniteElement/Hexahedron/HexahedronHdivFirstFamily.cxx"
#include "FiniteElement/Hexahedron/HexahedronHdivOptimalFirstFamily.cxx"
// first family of nedelec on wedges
#include "FiniteElement/Wedge/WedgeHdivFirstFamily.cxx"
#include "FiniteElement/Wedge/WedgeHdivOptimalFirstFamily.cxx"
// first family of nedelec on pyramids
#include "FiniteElement/Pyramid/PyramidHdivFirstFamily.cxx"
#include "FiniteElement/Pyramid/PyramidHdivOptimalFirstFamily.cxx"
// first family of nedelec on tets
#include "FiniteElement/Tetrahedron/TetrahedronHdivFirstFamily.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_HDIV
#include "FiniteElement/Triangle/TriangleDgOrtho.cxx"
#include "FiniteElement/Quadrangle/QuadrangleDgOrtho.cxx"
// optimal first family of nedelec on hexahedra
#include "FiniteElement/Hexahedron/HexahedronHdivOptimalHpFirstFamily.cxx"
#include "FiniteElement/Hexahedron/HexahedronHdivHpFirstFamily.cxx"
// first family of nedelec on pyramids
#include "FiniteElement/Pyramid/PyramidHdivOptimalHpFirstFamily.cxx"
#include "FiniteElement/Pyramid/PyramidHdivHpFirstFamily.cxx"
// first family of nedelec on wedges
#include "FiniteElement/Wedge/WedgeHdivOptimalHpFirstFamily.cxx"
#include "FiniteElement/Wedge/WedgeHdivHpFirstFamily.cxx"
// optimal hp first family of nedelec on tets
#include "FiniteElement/Tetrahedron/TetrahedronHdivOptimalHpFirstFamily.cxx"
#endif

#endif

// class storing various finite elements
#include "Harmonic/VarFiniteElement.cxx"

#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_FINITE_ELEMENT_HXX
#endif
