#ifndef MONTJOIE_FILE_MONTJOIE_HARMONIC_INLINE_HXX

// Common files
#include "FiniteElement/MontjoieFiniteElementInline.hxx"
#include "Output/MontjoieOutputInline.hxx"

#include "Solver/PreconditionerInline.cxx"
#include "Computation/AssembleMatrixInline.cxx"
#include "Solver/SolveHarmonicInline.cxx"

// generic definition of systems of equations solved by Montjoie
#include "Harmonic/GenericEquationInline.cxx"

// declaration of matrix-free matrices
// with these classes, matrix-vector product is done on the fly
#include "Computation/FemMatrixFreeClassInline.cxx"

// source
#include "Source/DefineSourceEllipticInline.cxx"
#include "Source/SourceSpatialeInline.cxx"
#include "Source/ModalSourceInline.cxx"
#include "Source/UserSourceInline.cxx"

// outputs
#include "Output/OutputHarmonicInline.cxx"

// boundary conditions
#include "Harmonic/BoundaryConditionHarmonicInline.cxx"

#ifdef MONTJOIE_WITH_TRANSMISSION
#include "Harmonic/TransmissionModelInline.cxx"
#endif

#include "Harmonic/GeneralizedImpedanceModelInline.cxx"

#include "Harmonic/TransparencyConditionInline.cxx"

#include "Elliptic/PhysicalPropertyInline.cxx"

#include "Harmonic/VarProblemBaseInline.cxx"
#include "Harmonic/VarAxisymProblemInline.cxx"
#include "Harmonic/VarGeometryProblemInline.cxx"
#include "Harmonic/DistributedProblemInline.cxx"
#include "Harmonic/VarProblemInline.cxx"
#include "Harmonic/VarHarmonicBaseInline.cxx"
#include "Harmonic/VarHarmonicInline.cxx"

#include "Solver/EigenvaluesHarmonicInline.cxx"

#ifdef MONTJOIE_WITH_ONE_DIM
#include "Elliptic/VarProblem1D_Inline.cxx"
#endif

#define MONTJOIE_FILE_MONTJOIE_HARMONIC_INLINE_HXX
#endif
