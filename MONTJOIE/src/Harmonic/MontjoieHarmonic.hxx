#ifndef MONTJOIE_FILE_MONTJOIE_HARMONIC_HXX

#include "MontjoieHarmonicHeader.hxx"
#include "MontjoieHarmonicInline.hxx"

#include "FiniteElement/MontjoieFiniteElement.hxx"
#include "Output/MontjoieOutput.hxx"

/****************
 * Common Files *
 ****************/

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/MontjoieHarmonicCompil.hxx"
#endif


// template files that may need additional instantiations
#include "Elliptic/PhysicalProperty.cxx"

#include "Harmonic/GenericEquation.cxx"

#ifdef MONTJOIE_WITH_NODAL_H1
// elementary matrices
#include "Computation/ElementaryMatrixH1.cxx"
#include "Computation/ProdMatVectH1.cxx"
#endif

#include "Harmonic/VarHarmonic.cxx"

// then templates files that are completely instantiated in lib/Compil/Harmonic
#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "Computation/AssembleMatrix.cxx"

#include "Solver/SolveHarmonic.cxx"

// declaration of matrix-free matrices
// with these classes, matrix-vector product is done on the fly
#include "Computation/FemMatrixFreeClass.cxx"

/***********************************
 * Classes for stationary problems *
 ***********************************/

#include "Elliptic/Maxwell/PhysicalConstant.cxx"

// declaration of preconditioners
#include "Solver/Preconditioner.cxx"

// class for sources
#include "Source/SourceSpatiale.cxx"
#include "Source/DefineSourceElliptic.cxx"
#include "Source/ModalSource.cxx"
#include "Source/UserSource.cxx"

// class for outputs
#include "Output/MeshInterpolationFEM.cxx"
#include "Output/OutputHarmonic.cxx"

// class for boundary conditions
#include "Harmonic/BoundaryConditionHarmonic.cxx"

#ifdef MONTJOIE_WITH_TRANSMISSION
#include "Harmonic/TransmissionModel.cxx"
#endif

#include "Harmonic/GeneralizedImpedanceModel.cxx"

#ifdef MONTJOIE_WITH_WIRES
// wires models 
#include "Harmonic/ModelWires.cxx"
#endif

// transparent condition
#include "Harmonic/TransparencyCondition.cxx"

// base classes for elliptic problems
#include "Harmonic/VarProblemBase.cxx"
#include "Harmonic/VarAxisymProblem.cxx"
#include "Harmonic/VarGeometryProblem.cxx"
#include "Harmonic/DistributedProblem.cxx"
#include "Harmonic/VarProblem.cxx"
#include "Harmonic/VarHarmonicBase.cxx"

// multi-frequency problems
#include "MultiFrequencyProblem.cxx"

// eigenvalues computation
#include "Solver/EigenvaluesHarmonic.cxx"

#ifdef MONTJOIE_WITH_ONE_DIM
#include "Elliptic/VarProblem1D.cxx"
#endif

#endif

#endif

#endif
