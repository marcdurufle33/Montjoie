#ifndef MONTJOIE_FILE_MONTJOIE_FINITE_ELEMENT_INLINE_HXX

#include "Mesh/MontjoieMeshInline.hxx"

#include "FiniteElement/ElementReferenceInline.cxx"
#include "FiniteElement/FaceReferenceInline.cxx"

#ifdef MONTJOIE_WITH_NODAL_H1
#include "FiniteElement/Quadrangle/QuadrangleGaussInline.cxx"
#include "FiniteElement/Quadrangle/QuadrangleLobattoInline.cxx"
#include "FiniteElement/Quadrangle/QuadrangleRadauInline.cxx"
#include "FiniteElement/Triangle/TriangleClassicalInline.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_H1
#include "FiniteElement/Quadrangle/QuadrangleHierarchicInline.cxx"
#include "FiniteElement/Triangle/TriangleHierarchicInline.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_DG
#include "FiniteElement/Quadrangle/QuadrangleDgGaussInline.cxx"
#endif

#ifdef MONTJOIE_WITH_ORTHO_DG
#include "FiniteElement/Triangle/TriangleDgOrthoInline.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HCURL
#include "FiniteElement/Triangle/TriangleHcurlFirstFamilyInline.cxx"
#include "FiniteElement/Triangle/TriangleHcurlSecondFamilyInline.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlFirstFamilyInline.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlGaussFirstFamilyInline.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlOptimalFirstFamilyInline.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlLobattoInline.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_HCURL
#include "FiniteElement/Quadrangle/QuadrangleHcurlOptimalHpFirstFamilyInline.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHcurlHpFirstFamilyInline.cxx"
#include "FiniteElement/Triangle/TriangleHcurlOptimalHpFirstFamilyInline.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HDIV
#include "FiniteElement/Quadrangle/QuadrangleHdivFirstFamilyInline.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHdivOptimalFirstFamilyInline.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_HDIV
#include "FiniteElement/Quadrangle/QuadrangleHdivOptimalHpFirstFamilyInline.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHdivHpFirstFamilyInline.cxx"
#include "FiniteElement/Triangle/TriangleHdivOptimalHpFirstFamilyInline.cxx"
#endif

// 3-D Finite elements

#ifdef MONTJOIE_WITH_THREE_DIM

#include "FiniteElement/VolumeReferenceInline.cxx"

#include "FiniteElement/Tetrahedron/TetrahedronReferenceInline.cxx"
#include "FiniteElement/Pyramid/PyramidReferenceInline.cxx"
#include "FiniteElement/Wedge/WedgeReferenceInline.cxx"
#include "FiniteElement/Hexahedron/HexahedronReferenceInline.cxx"

#ifdef MONTJOIE_WITH_NODAL_H1
#include "FiniteElement/Hexahedron/HexahedronGaussInline.cxx"
#include "FiniteElement/Hexahedron/HexahedronLobattoInline.cxx"
#include "FiniteElement/Tetrahedron/TetrahedronClassicalInline.cxx"
#include "FiniteElement/Pyramid/PyramidClassicalInline.cxx"
#include "FiniteElement/Wedge/WedgeClassicalInline.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_H1
#include "FiniteElement/Hexahedron/HexahedronHierarchicInline.cxx"
#include "FiniteElement/Tetrahedron/TetrahedronHierarchicInline.cxx"
#include "FiniteElement/Pyramid/PyramidHierarchicInline.cxx"
#include "FiniteElement/Wedge/WedgeHierarchicInline.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HCURL
#include "FiniteElement/Hexahedron/HexahedronHcurlFirstFamilyInline.cxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlOptimalFirstFamilyInline.cxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlLobattoInline.cxx"
#include "FiniteElement/Wedge/WedgeHcurlFirstFamilyInline.cxx"
#include "FiniteElement/Wedge/WedgeHcurlOptimalFirstFamilyInline.cxx"
#include "FiniteElement/Pyramid/PyramidHcurlFirstFamilyInline.cxx"
#include "FiniteElement/Pyramid/PyramidHcurlOptimalFirstFamilyInline.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_HCURL
#include "FiniteElement/Hexahedron/HexahedronHcurlOptimalHpFirstFamilyInline.cxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlHpFirstFamilyInline.cxx"
#include "FiniteElement/Pyramid/PyramidHcurlOptimalHpFirstFamilyInline.cxx"
#include "FiniteElement/Pyramid/PyramidHcurlHpFirstFamilyInline.cxx"
#include "FiniteElement/Wedge/WedgeHcurlOptimalHpFirstFamilyInline.cxx"
#include "FiniteElement/Wedge/WedgeHcurlHpFirstFamilyInline.cxx"
#include "FiniteElement/Tetrahedron/TetrahedronHcurlOptimalHpFirstFamilyInline.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_DG
#include "FiniteElement/Hexahedron/HexahedronDgGaussInline.cxx"
// #include "FiniteElement/Hexahedron/HexahedronDgJacobiInline.cxx"
#include "FiniteElement/Wedge/WedgeDgClassicalInline.cxx"
#endif

#ifdef MONTJOIE_WITH_ORTHO_DG
#include "FiniteElement/Hexahedron/HexahedronDgGaussInline.cxx"
#include "FiniteElement/Tetrahedron/TetrahedronDgOrthoInline.cxx"
#include "FiniteElement/Pyramid/PyramidDgOrthoInline.cxx"
#include "FiniteElement/Wedge/WedgeDgOrthoInline.cxx"
#endif

#ifdef MONTJOIE_WITH_LEGENDRE_DG
#include "FiniteElement/Hexahedron/HexahedronDgLegendreInline.cxx"
#include "FiniteElement/Pyramid/PyramidDgLegendreInline.cxx"
#include "FiniteElement/Wedge/WedgeDgLegendreInline.cxx"
#endif

#ifdef MONTJOIE_WITH_NODAL_HDIV
#include "FiniteElement/Hexahedron/HexahedronHdivFirstFamilyInline.cxx"
#include "FiniteElement/Hexahedron/HexahedronHdivOptimalFirstFamilyInline.cxx"
#include "FiniteElement/Wedge/WedgeHdivFirstFamilyInline.cxx"
#include "FiniteElement/Wedge/WedgeHdivOptimalFirstFamilyInline.cxx"
#include "FiniteElement/Pyramid/PyramidHdivFirstFamilyInline.cxx"
#include "FiniteElement/Pyramid/PyramidHdivOptimalFirstFamilyInline.cxx"
#endif

#ifdef MONTJOIE_WITH_HP_HDIV
#include "FiniteElement/Triangle/TriangleDgOrthoInline.cxx"
#include "FiniteElement/Pyramid/PyramidHdivOptimalHpFirstFamilyInline.cxx"
#include "FiniteElement/Pyramid/PyramidHdivHpFirstFamilyInline.cxx"
#include "FiniteElement/Wedge/WedgeHdivOptimalHpFirstFamilyInline.cxx"
#include "FiniteElement/Wedge/WedgeHdivHpFirstFamilyInline.cxx"
#endif

#endif

// class storing various finite elements
#include "Harmonic/VarFiniteElementInline.cxx"

#define MONTJOIE_FILE_MONTJOIE_FINITE_ELEMENT_INLINE_HXX
#endif
