#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/VarProblemBase.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template void VarProblem_Base::GetVelocityOnElements(VectReal_wp&, const Mesh<Dimension2>&);
  SELDON_EXTERN template void VarProblem_Base::GetVelocityOnElements(VectReal_wp&, const Mesh<Dimension3>&);

  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(FemMatrixFreeClass_Base<Real_wp>&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const GlobalGenericMatrix<Real_wp>&, int , int, CondensationBlockSolver_Fem<Real_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(FemMatrixFreeClass_Base<Real_wp>&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const GlobalGenericMatrix<Real_wp>&, int , int, CondensationBlockSolver_Fem<Real_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(FemMatrixFreeClass_Base<Complex_wp>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(FemMatrixFreeClass_Base<Complex_wp>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(FemMatrixFreeClass_Base<Complex_wp>&, DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(FemMatrixFreeClass_Base<Complex_wp>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const GlobalGenericMatrix<Real_wp>&, int , int, CondensationBlockSolver_Fem<Real_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const GlobalGenericMatrix<Real_wp>&, int , int, CondensationBlockSolver_Fem<Real_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(DistributedMatrix<Real_wp, General, BlockDiagRow>&, const GlobalGenericMatrix<Real_wp>&, int , int, CondensationBlockSolver_Fem<Real_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&, const GlobalGenericMatrix<Real_wp>&, int , int, CondensationBlockSolver_Fem<Real_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(DistributedMatrix<Complex_wp, General, BlockDiagRow>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(Matrix<Complex_wp, Symmetric, DiagonalRow>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(Matrix<Real_wp, Symmetric, DiagonalRow>&, const GlobalGenericMatrix<Real_wp>&, int , int, CondensationBlockSolver_Fem<Real_wp>*, bool);

  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(FemMatrixFreeClass_Base<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(FemMatrixFreeClass_Base<Real_wp>&, const GlobalGenericMatrix<Real_wp>&, int , int, CondensationBlockSolver_Fem<Real_wp>*, bool);

  SELDON_EXTERN template void VarComputationProblem::AddMatrixFEM(FemMatrixFreeClass_Base<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&, int , int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  SELDON_EXTERN template void VarComputationProblem::AddMatrixFEM(FemMatrixFreeClass_Base<Real_wp>&, const GlobalGenericMatrix<Real_wp>&, int , int, CondensationBlockSolver_Fem<Real_wp>*, bool);

  SELDON_EXTERN template void VarComputationProblem::AddMatrixFEM(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const GlobalGenericMatrix<Real_wp>&,
                                                                  int, int, CondensationBlockSolver_Fem<Real_wp>*, bool);
  
  SELDON_EXTERN template void VarComputationProblem::AddMatrixFEM(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const GlobalGenericMatrix<Complex_wp>&,
                                                                  int, int, CondensationBlockSolver_Fem<Complex_wp>*, bool);
  
  SELDON_EXTERN template void VarComputationProblem::ComputeDiagonalMatrix(Vector<Real_wp>&, const GlobalGenericMatrix<Real_wp>&, bool);
  SELDON_EXTERN template void VarComputationProblem::ComputeDiagonalMatrix(Vector<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&, bool);

  SELDON_EXTERN template void VarComputationProblem::ComputeDiagonalMatrix(Vector<Real_wp>&, const FemMatrixFreeClass_Base<Real_wp>&, const GlobalGenericMatrix<Real_wp>&, bool);
  SELDON_EXTERN template void VarComputationProblem::ComputeDiagonalMatrix(Vector<Complex_wp>&, const FemMatrixFreeClass_Base<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&, bool);

  SELDON_EXTERN template void VarComputationProblem::ComputeDiagonalMatrix(Vector<Complex_wp>&, const Matrix<Complex_wp, General, ArrayRowSparse>&, const GlobalGenericMatrix<Complex_wp>&);
  SELDON_EXTERN template void VarComputationProblem::ComputeDiagonalMatrix(Vector<Complex_wp>&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const GlobalGenericMatrix<Complex_wp>&);
  SELDON_EXTERN template void VarComputationProblem::ComputeDiagonalMatrix(Vector<Complex_wp>&, const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const GlobalGenericMatrix<Complex_wp>&);
  SELDON_EXTERN template void VarComputationProblem::ComputeDiagonalMatrix(Vector<Complex_wp>&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const GlobalGenericMatrix<Complex_wp>&);

  SELDON_EXTERN template void VarComputationProblem::ComputeDiagonalMatrix(Vector<Real_wp>&, const Matrix<Real_wp, General, ArrayRowSparse>&, const GlobalGenericMatrix<Real_wp>&);
  SELDON_EXTERN template void VarComputationProblem::ComputeDiagonalMatrix(Vector<Real_wp>&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const GlobalGenericMatrix<Real_wp>&);

  SELDON_EXTERN template void VarComputationProblem::ComputeBlockDiagonalMatrix(DistributedMatrix<Real_wp, General, BlockDiagRow>&, int, const FemMatrixFreeClass_Base<Real_wp>&, const GlobalGenericMatrix<Real_wp>&);
  SELDON_EXTERN template void VarComputationProblem::ComputeBlockDiagonalMatrix(DistributedMatrix<Complex_wp, General, BlockDiagRow>&, int, const FemMatrixFreeClass_Base<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&);
  SELDON_EXTERN template void VarComputationProblem::ComputeBlockDiagonalMatrix(DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&, int, const FemMatrixFreeClass_Base<Real_wp>&, const GlobalGenericMatrix<Real_wp>&);
  SELDON_EXTERN template void VarComputationProblem::ComputeBlockDiagonalMatrix(DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&, int, const FemMatrixFreeClass_Base<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&);
  
}
