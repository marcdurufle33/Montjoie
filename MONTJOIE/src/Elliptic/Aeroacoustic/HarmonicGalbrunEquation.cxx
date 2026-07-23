#ifndef MONTJOIE_FILE_HARMONIC_GALBRUN_EQUATION_CXX

namespace Montjoie
{

  /*****************
   * VarGalbrun_Eq *
   *****************/


  template<class TypeEquation>
  VarGalbrun_Eq<TypeEquation>::VarGalbrun_Eq()
    : VarHarmonic<TypeEquation>(),
      VarGalbrun_Dim<typename TypeEquation::Complexe,
		     typename TypeEquation::Dimension>(this->GetLeafClass())
  {
  }
  

  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {    
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
    VarGalbrun_Dim<Complexe, Dimension>::SetInputData(description_field, parameters);
  }

    
  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>::InitIndices(int n)
  {
    VarGalbrun_Dim<Complexe, Dimension>::InitIndices(n);
  }


  template<class TypeEquation>
  int VarGalbrun_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return VarGalbrun_Dim<Complexe, Dimension>::GetNbPhysicalIndices();
  }


  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    VarGalbrun_Dim<Complexe, Dimension>::SetIndices(i, parameters);
  }


  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    VarGalbrun_Dim<Complexe, Dimension>::SetPhysicalIndex(name_media, i, parameters);
  }


  template<class TypeEquation>
  string VarGalbrun_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return VarGalbrun_Dim<Complexe, Dimension>::GetPhysicalIndexName(m);
  }

  
  template<class TypeEquation>
  bool VarGalbrun_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return VarGalbrun_Dim<Complexe, Dimension>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  bool VarGalbrun_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return VarGalbrun_Dim<Complexe, Dimension>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  Real_wp VarGalbrun_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return VarGalbrun_Dim<Complexe, Dimension>::GetVelocityOfMedia(ref);
  }


  template<class TypeEquation>
  Real_wp VarGalbrun_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return VarGalbrun_Dim<Complexe, Dimension>::GetVelocityOfInfinity();
  }

  
  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
		      Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
		      IVect& num_ref, IVect& num_index, IVect& num_component,
		      Vector<bool>& compute_grad, Vector<bool>& compute_hess)
  {
    return VarGalbrun_Dim<Complexe, Dimension>
      ::GetVaryingIndices(rho_real, num_ref, num_index, num_component,
			  compute_grad, compute_hess);
  }

  
  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>
  ::ComputeMeshAndFiniteElement(const string & name_elt, bool split_mesh)
  {
    VarGalbrun_Dim<Complexe, Dimension>::CheckBoundaryCondition(this->mesh.GetBoundaryCondition());
    VarHarmonic<TypeEquation>::ComputeMeshAndFiniteElement(name_elt, split_mesh);
  }
  
  
  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>::ComputePhysicalCoefficients()
  {
    TypeEquation::SetIndexToCompute(*this);
    
    VarHarmonic<TypeEquation>::ComputePhysicalCoefficients();      
    VarGalbrun_Dim<Complexe, Dimension>::ComputePhysicalCoefficients();
  }


  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>::PerformOtherInitializations()
  {
    SetComplexOne(this->coefficient_volumic_source);
  }
  

  template<class TypeEquation> template<class T, class Dim>
  void VarGalbrun_Eq<TypeEquation>
  ::ModifyVolumetricSourceGen(int i, int j, const typename Dim::R_N& x,
                              const VirtualSourceField<T, Dim>& fsrc,
                              Vector<T>& f) const
  {
    abort();
  }

  
  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>
  ::ModifyVolumetricSourceGen(int i, int j, const typename Dimension::R_N& x,
                              const VirtualSourceField<Complexe, Dimension>& fsrc,
                              Vector<Complexe>& f) const
  {
    VarGalbrun_Dim<Complexe, Dimension>::ModifyVolumetricSource(i, j, x, fsrc, f);
  }


  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>
  ::ModifyVolumetricSource(int i, int j, const R2& pt,
                           const VirtualSourceField<Real_wp, Dimension2>& src,
                           Vector<Real_wp>& f) const
  {
    ModifyVolumetricSourceGen(i, j, pt, src, f);
  }
  

  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>
  ::ModifyVolumetricSource(int i, int j, const R2& pt,
                           const VirtualSourceField<Complex_wp, Dimension2>& src,
                           Vector<Complex_wp>& f) const
  {
    ModifyVolumetricSourceGen(i, j, pt, src, f);
  }


  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>
  ::ModifyVolumetricSource(int i, int j, const R3& pt,
                           const VirtualSourceField<Real_wp, Dimension3>& src,
                           Vector<Real_wp>& f) const
  {
    ModifyVolumetricSourceGen(i, j, pt, src, f);
  }
  
    
  template<class TypeEquation>
  void VarGalbrun_Eq<TypeEquation>
  ::ModifyVolumetricSource(int i, int j, const R3& pt,
                           const VirtualSourceField<Complex_wp, Dimension3>& src,
                           Vector<Complex_wp>& f) const
  {
    ModifyVolumetricSourceGen(i, j, pt, src, f);
  }
  
}

#define MONTJOIE_FILE_HARMONIC_GALBRUN_EQUATION_CXX
#endif
