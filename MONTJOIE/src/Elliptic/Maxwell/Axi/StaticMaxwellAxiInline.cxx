#ifndef MONTJOIE_FILE_STATIC_MAXWELL_AXI_INLINE_CXX

namespace Montjoie
{

  template<class TypeEquation>
  inline VolumetricSource<StaticMaxwellEquation_HcurlAxi>
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_MaxwellAxi<Real_wp>(var, param)
  {
  }
  

  inline bool StaticMaxwellEquation_HcurlAxi::SymmetricGlobalMatrix()
  {
    return false;
  }
  
  
  inline bool StaticMaxwellEquation_HcurlAxi::SymmetricElementaryMatrix()
  {
    return false; 
  }
  
  
  inline Real_wp StaticMaxwellEquation_HcurlAxi::GetCoefficientMassMatrix(const VarProblem_Base& vars)
  {
    return -1.0;
  }
  
  
  inline Real_wp StaticMaxwellEquation_HcurlAxi::GetCoefficientStiffnessMatrix(const VarProblem_Base& vars)
  {
    return 1.0;
  }
  
  
  inline Real_wp StaticMaxwellEquation_HcurlAxi::GetCoefMode(const VarProblem_Base& vars)
  {
    return 1.0;
  }
  

  inline bool StaticMaxwellEquation_Axi_DG::SymmetricGlobalMatrix()
  {
    return true; 
  }
  
  
  inline bool StaticMaxwellEquation_Axi_DG::SymmetricElementaryMatrix()
  {
    return false; 
  }
  
}

#define MONTJOIE_FILE_STATIC_MAXWELL_AXI_INLINE_CXX
#endif

