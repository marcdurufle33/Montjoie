#ifndef MONTJOIE_FILE_MONTJOIE_HARMONIC_HEADER_HXX


/****************
 * Common Files *
 ****************/

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "Output/MontjoieOutputHeader.hxx"

// functions to assemble finite element matrix
#include "Computation/AssembleMatrix.hxx"

// some functions that interface linear solvers with montjoie
#include "Solver/SolveHarmonic.hxx"


/***********************************
 * Classes for stationary problems *
 ***********************************/


// generic definition of systems of equations solved by Montjoie
#include "Harmonic/GenericEquation.hxx"

// declaration of matrix-free matrices
// with these classes, matrix-vector product is done on the fly (without storing the matrix)
#include "Computation/FemMatrixFreeClass.hxx"

// declaration of preconditioners
#include "Solver/Preconditioner.hxx"

// eigenvalues computation
#include "Solver/EigenvaluesHarmonic.hxx"

// class for sources
#include "Source/SourceSpatiale.hxx"
#include "Source/ModalSource.hxx"
#include "Source/DefineSourceElliptic.hxx"
#include "Source/UserSource.hxx"

// class for outputs
#include "Output/MeshInterpolationFEM.hxx"
#include "Output/OutputHarmonic.hxx"

// boundary conditions
#include "Harmonic/BoundaryConditionHarmonic.hxx"

// other models
#ifdef MONTJOIE_WITH_TRANSMISSION
#include "Harmonic/TransmissionModel.hxx"
#endif

#include "Harmonic/GeneralizedImpedanceModel.hxx"

#ifdef MONTJOIE_WITH_WIRES
#include "Harmonic/ModelWires.hxx"
#endif

// transparent condition
#include "Harmonic/TransparencyCondition.hxx"

// treatment of physical indices
#include "Elliptic/Maxwell/PhysicalConstant.hxx"
#include "Elliptic/PhysicalProperty.hxx"

// base classes of EllipticProblem
#include "Harmonic/VarProblemBase.hxx"
#include "Harmonic/VarAxisymProblem.hxx"
#include "Harmonic/VarGeometryProblem.hxx"
#include "Harmonic/DistributedProblem.hxx"
#include "Harmonic/VarProblem.hxx"
#include "Harmonic/VarHarmonicBase.hxx"
#include "Harmonic/VarHarmonic.hxx"

// multi-frequency problems
#include "MultiFrequencyProblem.hxx"

#ifdef MONTJOIE_WITH_NODAL_H1
// elementary matrices and matrix vector product
#include "Computation/ElementaryMatrixH1.hxx"
#include "Computation/ProdMatVectH1.hxx"
#endif

#ifdef MONTJOIE_WITH_ONE_DIM
#include "Elliptic/VarProblem1D.hxx"
#endif

#define MONTJOIE_FILE_MONTJOIE_HARMONIC_HEADER_HXX
#endif
