#ifndef MONTJOIE_FILE_GENERALIZED_IMPEDANCE_MODEL_INLINE_CXX

namespace Montjoie
{

  //! constructor with a given problem
  template<class TypeEquation>
  inline VarGeneralizedImpedance_Base
  ::VarGeneralizedImpedance_Base(EllipticProblem<TypeEquation>& var)
    : var_problem(var), var_boundary(var)
  {
    InitDefaultValues();
  }

  inline void VarGeneralizedImpedance_Base::InitGIBC()
  {
  }
  
  
  inline void VarGeneralizedImpedance_Base::InitReferences(int N)
  {
    ref_cond.Reallocate(N+1); ref_cond.Fill(0);
  }
  
  
  inline int VarGeneralizedImpedance_Base::GetNbVolumeDofOnSurface() const
  {
    return nb_ddl_volume_on_surface;
  }
  
  
  inline bool VarGeneralizedImpedance_Base::UsePadeApproximationSquareRoot() const
  {
    return pade_approx_sqrt; 
  }
  
  
  inline bool VarGeneralizedImpedance_Base::UseUnsymmetricImplementation() const
  {
    return use_unsymmetric_algo;
  }
  
  
  inline const Vector<int>& VarGeneralizedImpedance_Base::GetLocalVolumeDofNumber() const
  {
    return DdlVol;
  }

  
  inline void VarGeneralizedImpedance_Base
  ::AddGibcTerms(const Real_wp& alpha, const GlobalGenericMatrix<Real_wp>& nat_mat,
		 VirtualMatrix<Real_wp>& mat_sp, int offset_row, int offset_col)
  {
  } 

  
  inline void VarGeneralizedImpedance_Base
  ::AddGibcTerms(const Complex_wp& alpha, const GlobalGenericMatrix<Complex_wp>& nat_mat,
		 VirtualMatrix<Complex_wp>& mat_sp, int offset_row, int offset_col)
  {
  } 


  template<class TypeEquation>
  inline VarGeneralizedImpedance<TypeEquation>
  ::VarGeneralizedImpedance(EllipticProblem<TypeEquation>& var)
    : VarGeneralizedImpedance_Base(var), var_problem(var)
  {
  }
  
  
  template<class TypeEquation>
  inline void VarGeneralizedImpedance<TypeEquation>::TreatGibc(const IVect& Epart)
  {
    FindDofsGibc(var_problem, Epart);
  }
  
}

#define MONTJOIE_FILE_GENERALIZED_IMPEDANCE_MODEL_INLINE_CXX
#endif
