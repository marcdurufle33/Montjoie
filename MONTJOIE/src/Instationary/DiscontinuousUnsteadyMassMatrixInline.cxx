#ifndef MONTJOIE_FILE_DISCONTINUOUS_UNSTEADY_MASS_MATRIX_INLINE_CXX

namespace Montjoie
{

  //! default constructor
  template<class TypeEquation>
  inline DiscontinuousUnsteadyMassMatrix_Base::
  DiscontinuousUnsteadyMassMatrix_Base(HyperbolicProblem<TypeEquation>& vars)
    : VirtualMassMatrix(vars), var_computation(vars.var_harmonic), var_problem(vars.var_harmonic),
      var_boundary(vars.var_harmonic), var_time(vars), var_comm(vars.var_harmonic)
  {
    InitDefaultValues();
  }

  
  inline const VectReal_wp& DiscontinuousUnsteadyMassMatrix_Base::GetDiagonalDh() const
  {
    return Dh->GetDiagonal();
  }


  inline const VectReal_wp& DiscontinuousUnsteadyMassMatrix_Base::GetInverseDiagonalDh() const
  {
    return Dh->GetInverseDiagonal();
  }


  inline const VectReal_wp& DiscontinuousUnsteadyMassMatrix_Base::GetDiagonalSh() const
  {
    return Sh->GetDiagonal();
  }
  
  
  inline const VectReal_wp& DiscontinuousUnsteadyMassMatrix_Base::GetDiagonalDhMinusdtSh() const
  {
    return DhMinusdtSh->GetDiagonal(); 
  }
  
  
  inline const VectReal_wp& DiscontinuousUnsteadyMassMatrix_Base
  ::GetInverseDiagonalDhPlusdtSh() const
  {
    return DhPlusdtSh->GetInverseDiagonal();
  }
  

  inline const VectReal_wp& DiscontinuousUnsteadyMassMatrix_Base::GetDiagonalBh() const
  {
    return Bh->GetDiagonal();
  }

  
  inline const VectReal_wp& DiscontinuousUnsteadyMassMatrix_Base
  ::GetDiagonalBhMinusdtSh() const
  {
    return BhMinusdtSh->GetDiagonal();
  }
  
  
  inline const VectReal_wp& DiscontinuousUnsteadyMassMatrix_Base
  ::GetInverseDiagonalBhPlusdtSh() const
  {
    return BhPlusdtSh->GetInverseDiagonal();
  }
  

  inline DiscontinuousMassMatrix_Base& DiscontinuousUnsteadyMassMatrix_Base::GetOperatorDh()
  {
    return *Dh;
  }


  inline DiscontinuousMassMatrix_Base& DiscontinuousUnsteadyMassMatrix_Base::GetOperatorBh()
  {
    return *Bh;
  }


  inline DiscontinuousMassMatrix_Base& DiscontinuousUnsteadyMassMatrix_Base::GetOperatorSh()
  {
    return *Sh;
  }


  inline DiscontinuousMassMatrix_Base& DiscontinuousUnsteadyMassMatrix_Base::GetOperatorShVec()
  {
    return *ShVec;
  }


  //! returns operator C_h (for HDG formulation only)
  inline DiscontinuousMassMatrix_Base& DiscontinuousUnsteadyMassMatrix_Base::GetOperatorCh()
  {
    return *ChSurf;
  }
  
  
  //! constructor
  template<class Dimension> template<class TypeEquation>
  inline DiscontinuousUnsteadyMassMatrix<Dimension>
  ::DiscontinuousUnsteadyMassMatrix(HyperbolicProblem<TypeEquation>& var)
    : DiscontinuousUnsteadyMassMatrix_Base(var), var_problem(var.var_harmonic), var_time(var)
  {
  }



}

#define MONTJOIE_FILE_DISCONTINUOUS_UNSTEADY_MASS_MATRIX_INLINE_CXX
#endif

