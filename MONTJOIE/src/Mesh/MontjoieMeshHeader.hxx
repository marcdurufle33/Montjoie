#ifndef MONTJOIE_FILE_MONTJOIE_MESH_HEADER_HXX

/****************
 * Common Files *
 ****************/

#include "Quadrature/MontjoieQuadratureHeader.hxx"

// output functions
#include "Output/CommonInputOutput.hxx"

// linear and non-linear solvers
#include "Solver/MontjoieSolverHeader.hxx"


/********
 * Mesh *
 ********/

// numbering stuff for meshes
#include "Mesh/NumberMesh.hxx"

// declaration of points and matrices on the reference element (unit cube, unit tet ...)
#include "FiniteElement/PointsReference.hxx"
// interpolation from nodal points to another sets of points
#include "FiniteElement/ProjectionOperator.hxx"

#include "FiniteElement/ElementReferenceBase.hxx"
#include "FiniteElement/ElementReference.hxx"
#include "FiniteElement/Edge/EdgeReference.hxx"

// geometric elements
#include "FiniteElement/ElementGeomReference.hxx"
#include "FiniteElement/FaceGeomReference.hxx"
#include "FiniteElement/Quadrangle/QuadrangleGeomReference.hxx"
#include "FiniteElement/Triangle/TriangleGeomReference.hxx"

#ifdef MONTJOIE_WITH_THREE_DIM
#include "FiniteElement/VolumeGeomReference.hxx"
#include "FiniteElement/Hexahedron/HexahedronGeomReference.hxx"
#include "FiniteElement/Pyramid/PyramidGeomReference.hxx"
#include "FiniteElement/Wedge/WedgeGeomReference.hxx"
#include "FiniteElement/Tetrahedron/TetrahedronGeomReference.hxx"
#endif

// declarations of edges, faces, hex, tets ...
#include "Mesh/MeshElement.hxx"

// base class for 2-D and 3-D meshes
#include "Mesh/MeshBase.hxx"
// basic methods for handling boundaries and curved boundaries of a mesh
#include "Mesh/MeshBoundaries.hxx"

#ifdef MONTJOIE_WITH_ONE_DIM
// main class for 1-D mesh
#include "Mesh/Mesh1D.hxx"
#endif

// methods for handling curves for 2-D meshes
#include "Mesh/Mesh2DBoundaries.hxx"
// main class for 2-D mesh
#include "Mesh/Mesh2D.hxx"

#ifdef MONTJOIE_WITH_THREE_DIM
// methods for handling curves for 3-D meshes
#include "Mesh/Mesh3DBoundaries.hxx"
// main class for 3-D mesh
#include "Mesh/Mesh3D.hxx"
#endif

// functions for meshes
#include "Mesh/FunctionsMesh.hxx"

#ifdef SELDON_WITH_MPI
#include "Mesh/ParallelMeshFunctions.hxx"
#endif

#define MONTJOIE_FILE_MONTJOIE_MESH_HEADER_HXX
#endif
