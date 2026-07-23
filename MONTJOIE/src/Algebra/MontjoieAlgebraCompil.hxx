#ifndef MONTJOIE_FILE_MONTJOIE_ALGEBRA_COMPIL_HXX

// files in Seldon directory
#include "Compil/Seldon/Array3D.cpp"
#include "Compil/Seldon/BandMatrix.cpp"
#include "Compil/Seldon/Blas.cpp"
#include "Compil/Seldon/CholeskySolver.cpp"
#include "Compil/Seldon/Cholmod.cpp"
#include "Compil/Seldon/Common.cpp"
#include "Compil/Seldon/EigenvalueSolver.cpp"
#include "Compil/Seldon/FunctionsMatrixArray.cpp"
#include "Compil/Seldon/FunctionsMatrixComplex.cpp"
#include "Compil/Seldon/FunctionsMatrix.cpp"
#include "Compil/Seldon/FunctionsMatrixDense.cpp"
#include "Compil/Seldon/FunctionsMatrixVector.cpp"
#include "Compil/Seldon/FunctionsMatVectComplex.cpp"

#ifdef SELDON_WITH_PRECONDITIONING
#include "Compil/Seldon/IlutPreconditioning.cpp"
#endif

#include "Compil/Seldon/IOMatrixMarket.cpp"
#include "Compil/Seldon/IterativeSolver.cpp"
#include "Compil/Seldon/Lapack.cpp"
#include "Compil/Seldon/MatrixComplexConversions.cpp"
#include "Compil/Seldon/MatrixComplexSparse.cpp"
#include "Compil/Seldon/MatrixConversion.cpp"
#include "Compil/Seldon/MatrixHermSymTriang.cpp"
#include "Compil/Seldon/MatrixPacked.cpp"
#include "Compil/Seldon/MatrixPointers.cpp"
#include "Compil/Seldon/MatrixSparse.cpp"
#include "Compil/Seldon/PermutationScalingMatrix.cpp"
#include "Compil/Seldon/RelaxationMatrixVector.cpp"
#include "Compil/Seldon/SparseDirectSolver.cpp"
#include "Compil/Seldon/SparseSeldonSolver.cpp"
#include "Compil/Seldon/Vector.cpp"

#ifdef SELDON_WITH_MUMPS
#include "Compil/Seldon/Mumps.cpp"
#endif

#ifdef SELDON_WITH_MPI
#include "Compil/Seldon/DistributedMatrix.cpp"
#include "Compil/Seldon/DistributedVector.cpp"
#endif

// files in Algebra directory
#include "Compil/Algebra/additional_function.cpp"
#include "Compil/Algebra/matrix_blockdiagonal.cpp"
#include "Compil/Algebra/skyline_matrix.cpp"
#include "Compil/Algebra/tiny_band_matrix.cpp"
#include "Compil/Algebra/tiny_vector.cpp"

#define MONTJOIE_FILE_MONTJOIE_ALGEBRA_COMPIL_HXX
#endif

