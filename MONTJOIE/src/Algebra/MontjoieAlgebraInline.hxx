#ifndef MONTJOIE_FILE_MONTJOIE_ALGEBRA_INLINE_HXX

// including all files containing inlined fonctions
#include "Algebra/TinyArray3D_Inline.cxx"
#include "Algebra/SeldonRedirectionInline.cxx"

#include "SeldonInline.hxx"
#include "SeldonSolverInline.hxx"
#include "SeldonComplexMatrixInline.hxx"

#include "matrix_sparse/BandMatrixInline.cxx"

#include "Algebra/TinyArray3D_Inline.cxx"
#include "Algebra/DiagonalMatrixInline.cxx"
#include "Algebra/MatrixBlockDiagonalInline.cxx"
#include "Algebra/GeneralSkyLineMatrixInline.cxx"
#include "Algebra/SymmetricSkyLineMatrixInline.cxx"
#include "Algebra/FunctionMatrixExtractionInline.cxx"
#include "Algebra/TinyBandMatrixInline.cxx"
#include "Algebra/TinySymmetricTensorInline.cxx"

#ifdef SELDON_WITH_MPI
#include "Algebra/DistributedBlockDiagonalMatrixInline.cxx"
#endif

// Common files
#include "Share/PrecisionInline.cxx"
#include "Share/MontjoieTypesInline.cxx"
#include "Share/CommonMontjoieInline.cxx"

#define MONTJOIE_FILE_MONTJOIE_ALGEBRA_INLINE_HXX
#endif
