#ifndef MONTJOIE_FILE_CONTINUOUS_MASS_MATRIX_INLINE_CXX

namespace Montjoie
{
  
  inline int ContinuousMassMatrix_Base::GetMatrixType() const
  {
    return type_matrix;
  }
  
  
  inline void ContinuousMassMatrix_Base::SetDiagonal(const VectReal_wp& D)
  {
    type_matrix = FemMassMatrix::DIAGONAL;
    diagonal = D; 
  }
  
  
  inline void ContinuousMassMatrix_Base::SetInverseDiagonal(const VectReal_wp& invD)
  {
    type_matrix = FemMassMatrix::DIAGONAL;
    invDiagonal = invD;
  }
  
  
  inline VectReal_wp& ContinuousMassMatrix_Base::GetDiagonal()
  {
    return diagonal; 
  }
  
  
  inline const VectReal_wp& ContinuousMassMatrix_Base::GetDiagonal() const
  {
    return diagonal; 
  }
  
  
  inline VectReal_wp& ContinuousMassMatrix_Base::GetInverseDiagonal()
  {
    return invDiagonal;
  }
  
  
  inline const VectReal_wp& ContinuousMassMatrix_Base::GetInverseDiagonal() const
  {
    return invDiagonal;
  }

  
  inline DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& 
  ContinuousMassMatrix_Base::GetBlockDiagonal()
  {
    return block_diagonal; 
  }
  
  
  inline DistributedMatrix<Real_wp, General, BlockDiagRow>& 
  ContinuousMassMatrix_Base::GetUnsymmetricBlockDiagonal()
  {
    return block_diagonal_unsym; 
  }


  inline DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& 
  ContinuousMassMatrix_Base::GetInverseBlockDiagonal()
  {
    return invBlock_diagonal; 
  }
  
  
  inline DistributedMatrix<Real_wp, General, BlockDiagRow>& 
  ContinuousMassMatrix_Base::GetInverseUnsymmetricBlockDiagonal()
  {
    return invBlock_diagonal_unsym; 
  }

  
  inline DistributedMatrix<Real_wp, General, ArrayRowSparse>& 
  ContinuousMassMatrix_Base::GetSparseMatrix()
  {
    return mat_sparse_unsym;
  }


  inline DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& 
  ContinuousMassMatrix_Base::GetSymSparseMatrix()
  {
    return mat_sparse;
  }

}

#define MONTJOIE_FILE_CONTINUOUS_MASS_MATRIX_INLINE_CXX
#endif
