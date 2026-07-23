#ifndef MONTJOIE_FILE_MONTJOIE_ALGEBRA_HXX

#include "MontjoieAlgebraHeader.hxx"
#include "MontjoieAlgebraInline.hxx"

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "MontjoieAlgebraCompil.hxx"
#endif

// Inclusion of Seldon and files related to algebra stuff

// Seldon library is included
#include "Seldon.hxx"
// and complex matrices
#include "SeldonComplexMatrix.hxx"
//  and solvers
#include "SeldonSolver.hxx"

#ifndef SELDON_WITH_MPI
#include "computation/solver/DistributedSolver.cxx"
#endif

#include "Algebra/TinyArray3D.cxx"
#include "Algebra/DiagonalMatrix.hxx"
#include "Algebra/MatrixBlockDiagonal.cxx"
//#include "Algebra/OptMatrixBlockDiagonal.cxx"

#include "matrix_sparse/BandMatrix.cxx"
#include "vector/Vector2.cxx"
#include "vector/Vector3.cxx"

#include "Algebra/TinyBandMatrix.cxx"
#include "Algebra/TinySymmetricTensor.cxx"
#include "Algebra/TinyBlockSolver1D.cxx"

// Seldon-like classes
#include "Algebra/SymmetricSkyLineMatrix.cxx"
#include "Algebra/GeneralSkyLineMatrix.cxx"
#include "Algebra/FactorisationLU.cxx"
#include "Algebra/Eigenvalue.cxx"

#ifdef SELDON_WITH_MPI
// distributed vectors and matrices
#include "Algebra/DistributedBlockDiagonalMatrix.cxx"
#include "Algebra/DistributedDenseMatrix.cxx"
#endif

#include "Algebra/FunctionMatrixExtraction.cxx"

#include "Share/MontjoieTypes.cxx"

// including Scalapack interface
#ifndef SELDON_WITH_COMPILED_LIBRARY
#ifdef SELDON_WITH_SCALAPACK
#include "Algebra/Scalapack.cxx"
#endif
#endif

#include "Share/CommonMontjoie.cxx"
#include "Share/Precision.cxx"

#endif


#define MONTJOIE_FILE_MONTJOIE_ALGEBRA_HXX
#endif
