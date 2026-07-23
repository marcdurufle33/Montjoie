#ifndef MONTJOIE_FILE_MONTJOIE_ALGEBRA_HEADER_HXX

#ifdef SELDON_WITH_COMPILED_LIBRARY
// for fast compilation, all the flags are defined
#include "MontjoieFlag.hxx"
#endif

// NewAlloc as default allocator in order to avoid problems
// with vectors of complex types
#define SELDON_DEFAULT_ALLOCATOR NewAlloc
// seldon will call abort() when encountering an exception
#define SELDON_WITH_ABORT
// no call of srand by Seldon
#define SELDON_WITHOUT_REINIT_RANDOM
// use of virtual methods in Seldon for a fast compilation
#define SELDON_WITH_VIRTUAL
// use of preconditioning
#define SELDON_WITH_PRECONDITIONING
// to avoid too many instantiations in Seldon
#define SELDON_WITH_REDUCED_TEMPLATE

// including packages of STL not included in Seldon
#include <list>
#include <iomanip>
#include <unistd.h>
#include <stdint.h>

// To display a variable (with its name and rank processor)
#ifndef DISP
// #define DISP(x) if (rank_global_processor ==0) { cout << #x ": " << (x) << endl; } else { cerr << #x ": " << (x) << endl; }
#ifdef SELDON_WITH_MPI
#define DISP(x) std::cout << rank_global_processor << " " << #x ": " << (x) << std::endl
#endif
#endif

// default precision for mpfr reals
#define MONTJOIE_DEFAULT_PRECISION 500

// size for mpfr reals (size of array _mpfr_d)
#define MONTJOIE_SIZE_MPFR (MONTJOIE_DEFAULT_PRECISION-1)/64+1

//! Seed used for random generation
extern unsigned int Seed;

//! Define if we have ALREADY put on the seed
extern bool is_seed_defined;

//! global processor rank
extern int rank_global_processor;

//! global variables for parallel computation
extern double size_all_communications;

#ifdef _OPENMP
#include "omp.h"
#endif

#ifdef SELDON_WITH_MPI
#include "mpi.h"
#endif

#if defined(MONTJOIE_WITH_MPFR) || defined(MONTJOIE_WITH_FLOAT128) || defined(MONTJOIE_WITH_FLOAT80)
// double is not used
#define MONTJOIE_WITH_MULTIPLE
#endif

// including definition of Real_wp and precision-machine constants
#include "Share/Precision.hxx"
  
// Inclusion of Seldon and files related to algebra stuff

// Seldon library is included
#include "SeldonHeader.hxx"
#include "SeldonSolverHeader.hxx"
#include "SeldonComplexMatrixHeader.hxx"

#include "matrix_sparse/BandMatrix.hxx"
#include "vector/Vector2.hxx"
#include "vector/Vector3.hxx"

namespace Montjoie
{
  // seldon library is used directly without writing Seldon::
  using namespace Seldon;
}

#include "Algebra/SeldonRedirection.hxx"

#include "Algebra/TinyArray3D.hxx"

#ifndef SELDON_WITH_MPI
#include "computation/solver/DistributedSolver.hxx"
#endif

#include "Algebra/DiagonalMatrix.hxx"
#include "Algebra/MatrixBlockDiagonal.hxx"
//#include "Algebra/OptMatrixBlockDiagonal.hxx"

#include "Algebra/SymmetricSkyLineMatrix.hxx"
#include "Algebra/GeneralSkyLineMatrix.hxx"
#include "Algebra/FactorisationLU.hxx"
#include "Algebra/Eigenvalue.hxx"
#include "Algebra/FunctionMatrixExtraction.hxx"
#include "Algebra/TinyBandMatrix.hxx"
#include "Algebra/TinySymmetricTensor.hxx"
#include "Algebra/TinyBlockSolver1D.hxx"

#ifdef SELDON_WITH_MPI
#include "Algebra/DistributedBlockDiagonalMatrix.hxx"
#include "Algebra/DistributedDenseMatrix.hxx"
#endif

// including Scalapack interface
#ifdef SELDON_WITH_SCALAPACK
#include "Algebra/Scalapack.hxx"
#endif

// list of typedef (Real_wp, Complex_wp, vectors, matrices, etc)
#include "Share/MontjoieTypes.hxx"

#include "Share/CommonMontjoie.hxx"

#define MONTJOIE_FILE_MONTJOIE_ALGEBRA_HEADER_HXX
#endif
