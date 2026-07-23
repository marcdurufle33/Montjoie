#ifndef MONTJOIE_FILE_MONTJOIE_MESH_HXX

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#include "Quadrature/MontjoieQuadrature.hxx"
#include "Solver/MontjoieSolver.hxx"

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/MontjoieMeshCompil.hxx"
#endif

#include "Output/CommonInputOutput.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "Solver/NonLinearEquations.cxx"
#include "Solver/NonLinearLeastSquares.cxx"
#include "FiniteElement/ProjectionOperator.cxx"

#include "FiniteElement/ElementReferenceBase.cxx"
#include "FiniteElement/ElementReference.cxx"

#include "FiniteElement/Edge/EdgeReference.cxx"
#include "FiniteElement/PointsReference.cxx"

#include "FiniteElement/ElementGeomReference.cxx"

#include "FiniteElement/FaceGeomReference.cxx"
#include "FiniteElement/Triangle/TriangleGeomReference.cxx"
#include "FiniteElement/Quadrangle/QuadrangleGeomReference.cxx"

#ifdef MONTJOIE_WITH_THREE_DIM

// base class for 3-D finite element
#include "FiniteElement/VolumeGeomReference.cxx"
#include "FiniteElement/Tetrahedron/TetrahedronGeomReference.cxx"
#include "FiniteElement/Hexahedron/HexahedronGeomReference.cxx"
#include "FiniteElement/Pyramid/PyramidGeomReference.cxx"
#include "FiniteElement/Wedge/WedgeGeomReference.cxx"

#endif

#include "Mesh/MeshElement.cxx"
#include "Mesh/MeshBase.cxx"
#include "Mesh/NumberMap.cxx"
#include "Mesh/NumberMesh.cxx"
#include "Mesh/MeshBoundaries.cxx"

#ifdef MONTJOIE_WITH_ONE_DIM
#include "Mesh/Mesh1D.cxx"
#endif

#include "Mesh/Mesh2DBoundaries.cxx"
#include "Mesh/Mesh2D.cxx"

#ifdef MONTJOIE_WITH_THREE_DIM
#include "Mesh/Mesh3DBoundaries.cxx"
#include "Mesh/Mesh3D.cxx"
#endif

#include "Mesh/FunctionsMesh.cxx"

#ifdef SELDON_WITH_MPI
#include "Mesh/ParallelMeshFunctions.cxx"
#endif

#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_MESH_HXX
#endif
