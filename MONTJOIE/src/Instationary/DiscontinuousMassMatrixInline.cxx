#ifndef MONTJOIE_FILE_DISCONTINUOUS_MASS_MATRIX_INLINE_CXX

namespace Montjoie
{

  //! sets the threshold used to convert a dense matrix to a sparse matrix
  inline void DiscontinuousMassMatrix_Base::SetThreshold(const Real_wp& eps)
  {
    threshold = eps;
  }
  
  
  //! returns the threshold used to convert a dense matrix to a sparse matrix
  inline Real_wp DiscontinuousMassMatrix_Base::GetThreshold() const
  {
    return threshold;
  }
  

  //! sets the stopping criterion used by conjugate gradient
  inline void DiscontinuousMassMatrix_Base::SetStoppingCriterion(const Real_wp& eps)
  {
    stopping_criterion = eps;
  }
  
  
  //! returns the stopping criterion used by conjugate gradient
  inline Real_wp DiscontinuousMassMatrix_Base::GetStoppingCriterion() const
  {
    return stopping_criterion;
  }


  //! returns the number of elements
  inline int DiscontinuousMassMatrix_Base::GetNbElt() const
  {
    return var_problem.GetNbElt();
  }
  
  
  //! returns the type of matrix 
  inline int DiscontinuousMassMatrix_Base::GetMatrixType() const
  {
    return type_matrix; 
  }
  

  //! sets the unknown number
  inline void DiscontinuousMassMatrix_Base::SetUnknownNumber(int m)
  {
    num_unknown = m; 
  }
  

  //! sets mass and damping coefficient
  inline void DiscontinuousMassMatrix_Base::SetCoefficient(const Real_wp& a, const Real_wp& b)
  {
    coef_mass_ = a;
    coef_damp_ = b;
  }
  

  //! returns the diagonal (if the mass matrix is diagonal)
  inline VectReal_wp& DiscontinuousMassMatrix_Base::GetDiagonal()
  {
    return diagonal;
  }
  

  //! returns the diagonal (if the mass matrix is diagonal)
  inline const VectReal_wp& DiscontinuousMassMatrix_Base::GetDiagonal() const
  {
    return diagonal; 
  }
  

  //! returns the inverse of diagonal (if the mass matrix is diagonal)
  inline VectReal_wp& DiscontinuousMassMatrix_Base::GetInverseDiagonal()
  {
    return invDiagonal; 
  }
  

  //! returns the inverse diagonal (if the mass matrix is diagonal)
  inline const VectReal_wp& DiscontinuousMassMatrix_Base::GetInverseDiagonal() const
  {
    return invDiagonal;
  }


  //! return the offset for each block
  inline const IVect& DiscontinuousMassMatrix_Base::GetOffsetSizeBlocks() const
  {
    return block_size;
  }
  
  
  //! returns the distribution of level times
  inline MatrixVectorProductLevel& DiscontinuousMassMatrix_Base::GetLevelTime()
  {
    return level_time;
  }
  
}

#define MONTJOIE_FILE_DISCONTINUOUS_MASS_MATRIX_INLINE_CXX
#endif
