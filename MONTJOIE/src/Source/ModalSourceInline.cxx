#ifndef MONTJOIE_FILE_MODAL_SOURCE_INLINE_CXX

namespace Montjoie
{

  //! constructor with a given problem
  template<class TypeEquation>
  inline ModalSourceBoundary_Base
  ::ModalSourceBoundary_Base(const EllipticProblem<TypeEquation>& var)
  {
    InitDefaultValues(var.mesh.GetNbReferences());
  }


  inline int ModalSourceBoundary_Base::GetNbModesToCompute() const
  {
    return nb_modes_to_compute;
  }
  
  
  inline const Real_wp& ModalSourceBoundary_Base::GetEigenvalueShift() const
  {
    return shift_eigenval_mode; 
  }
  
  
  inline const VectReal_wp& ModalSourceBoundary_Base::GetCoefficientPlane() const
  {
    return coef_plane; 
  }
  
  
  inline const IVect& ModalSourceBoundary_Base::GetModeNumberToCombine() const
  {
    return number_mode_combined;
  }
  
  
  inline int ModalSourceBoundary_Base::GetLateralBoundaryCondition() const
  {
    return lateral_condition; 
  }
  
  
  inline int ModalSourceBoundary_Base::GetVolumeReference() const
  {
    return reference_volume; 
  }
  
  
  //! constructor with a given problem
  template<class T, class Dimension> template<class TypeEquation>
  inline ModalSourceBoundary_Dim<T, Dimension>
  ::ModalSourceBoundary_Dim(const EllipticProblem<TypeEquation>& var)
    : ModalSourceBoundary_Base(var), var_problem(var)
  {
  }

  template<class T, class Dimension>
  inline ModalSourceBoundary_Dim<T, Dimension>::~ModalSourceBoundary_Dim()
  {
  }
  
}

#define MONTJOIE_FILE_MODAL_SOURCE_INLINE_CXX
#endif

