#ifndef MONTJOIE_FILE_CONTINUOUS_UNSTEADY_MASS_MATRIX_INLINE_CXX

namespace Montjoie
{

  template<class TypeEquation>
  inline VirtualMassMatrix::VirtualMassMatrix(HyperbolicProblem<TypeEquation>& var)
    : var_time(var)
  {
    InitDefaultValues();
  }


  //! default constructor
  template<class TypeEquation>
  inline ContinuousUnsteadyMassMatrix_Base::
  ContinuousUnsteadyMassMatrix_Base(HyperbolicProblem<TypeEquation>& vars)
    : VirtualMassMatrix(vars), var_computation(vars.var_harmonic), var_problem(vars.var_harmonic),
      var_comm(vars.var_harmonic), var_boundary(vars.var_harmonic), var_time(vars)
  {
    InitDefaultValues();
  }

  
  inline const VectReal_wp& ContinuousUnsteadyMassMatrix_Base::GetDiagonalDh() const
  {
    return Dh->GetDiagonal(); 
  }


  inline const VectReal_wp& ContinuousUnsteadyMassMatrix_Base::GetInverseDiagonalDh() const
  {
    return Dh->GetInverseDiagonal(); 
  }


  inline const VectReal_wp& ContinuousUnsteadyMassMatrix_Base::GetDiagonalSh() const
  {
    return Sh->GetDiagonal(); 
  }
  
  
  inline const VectReal_wp& ContinuousUnsteadyMassMatrix_Base::GetDiagonalDhMinusdtSh() const
  {
    return DhMinusdtSh->GetDiagonal();
  }
  
  
  inline const VectReal_wp& ContinuousUnsteadyMassMatrix_Base::GetInverseDiagonalDhPlusdtSh() const
  {
    return DhPlusdtSh->GetInverseDiagonal();
  }
  

  inline const VectReal_wp& ContinuousUnsteadyMassMatrix_Base::GetDiagonalBh() const
  {
    return Bh->GetDiagonal(); 
  }

  
  inline const VectReal_wp& ContinuousUnsteadyMassMatrix_Base::GetDiagonalBhMinusdtSh() const
  {
    return BhMinusdtSh->GetDiagonal(); 
  }
  
  
  inline const VectReal_wp& ContinuousUnsteadyMassMatrix_Base::GetInverseDiagonalBhPlusdtSh() const
  {
    return BhPlusdtSh->GetInverseDiagonal();
  }
  
  
  inline void ContinuousUnsteadyMassMatrix_Base::SetDiagonalDh(const VectReal_wp& D)
  {
    Dh->SetDiagonal(D); 
  }
  
  
  inline void ContinuousUnsteadyMassMatrix_Base::SetDiagonalSh(const VectReal_wp& D)
  {
    Sh->SetDiagonal(D); 
  }
  
  
  inline void ContinuousUnsteadyMassMatrix_Base::SetInverseDiagonalDhPlusdtSh(const VectReal_wp& D)
  {
    DhPlusdtSh->SetInverseDiagonal(D);
  }


  template<class Dimension> template<class TypeEquation>
  inline ContinuousUnsteadyMassMatrix_Dim<Dimension>
  ::ContinuousUnsteadyMassMatrix_Dim(HyperbolicProblem<TypeEquation>& vars)
    : ContinuousUnsteadyMassMatrix_Base(vars)
  {
  }
    
  
  template<class Dimension> template<class TypeEquation>
  inline ContinuousUnsteadyMassMatrix<Dimension>
  ::ContinuousUnsteadyMassMatrix(HyperbolicProblem<TypeEquation>& vars)
    : ContinuousUnsteadyMassMatrix_Dim<Dimension>(vars), var_problem(vars.var_harmonic)
  {
  }
  
}

#define MONTJOIE_FILE_CONTINUOUS_UNSTEADY_MASS_MATRIX_INLINE_CXX
#endif
