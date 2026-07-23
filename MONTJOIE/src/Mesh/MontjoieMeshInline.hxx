#ifndef MONTJOIE_FILE_MONTJOIE_MESH_INLINE_HXX

#include "Quadrature/MontjoieQuadratureInline.hxx"
#include "Solver/MontjoieSolverInline.hxx"

#include "Mesh/MeshElementInline.cxx"
#include "Mesh/MeshBaseInline.cxx"
#ifdef MONTJOIE_WITH_ONE_DIM
#include "Mesh/Mesh1DInline.cxx"
#endif
#include "Mesh/NumberMapInline.cxx"
#include "Mesh/NumberMeshInline.cxx"
#include "Mesh/MeshBoundariesInline.cxx"
#include "Mesh/Mesh2DBoundariesInline.cxx"
#include "Mesh/Mesh2DInline.cxx"

#ifdef MONTJOIE_WITH_THREE_DIM
#include "Mesh/Mesh3DInline.cxx"
#include "Mesh/Mesh3DBoundariesInline.cxx"
#endif

#include "Output/CommonInputOutputInline.cxx"

#include "FiniteElement/PointsReferenceInline.cxx"
#include "FiniteElement/ProjectionOperatorInline.cxx"

#include "FiniteElement/ElementReferenceBaseInline.cxx"
#include "FiniteElement/ElementReferenceInline.cxx"
#include "FiniteElement/Edge/EdgeReferenceInline.cxx"
#include "FiniteElement/ElementGeomReferenceInline.cxx"
#include "FiniteElement/FaceGeomReferenceInline.cxx"

#include "FiniteElement/Triangle/TriangleGeomReferenceInline.cxx"
#include "FiniteElement/Quadrangle/QuadrangleGeomReferenceInline.cxx"


#ifdef MONTJOIE_WITH_THREE_DIM
#include "FiniteElement/VolumeGeomReferenceInline.cxx"
#include "FiniteElement/Tetrahedron/TetrahedronGeomReferenceInline.cxx"
#include "FiniteElement/Pyramid/PyramidGeomReferenceInline.cxx"
#include "FiniteElement/Wedge/WedgeGeomReferenceInline.cxx"
#include "FiniteElement/Hexahedron/HexahedronGeomReferenceInline.cxx"
#endif

#define MONTJOIE_FILE_MONTJOIE_MESH_INLINE_HXX
#endif
