#ifndef MONTJOIE_FILE_MAXWELL_AXISYM_HARMONIC_INLINE_CXX

namespace Montjoie
{
 
  /***********************
   * HarmonicMaxwell_Axi *
   ***********************/
  

  //! default constructor
  template<class Complexe> template<class TypeEquation>
  inline HarmonicMaxwell_Axi<Complexe>::HarmonicMaxwell_Axi(EllipticProblem<TypeEquation>& var)
    : VarAxisymProblem(var), var_problem(var), var_boundary(var)
  {
    // for a vectorial equation -1 and +1 are preferred modes
    name_file_outputJ = "HyCircle.dat";
    nb_subdiv_outputJ = 0; output_Hy_teta = true;
    
    epsilon0 = 1.0; mu0 = 1.0; invMu0 = 1.0;
  }
  
   
  //! returns the number of physical media
  template<class Complexe>
  inline int HarmonicMaxwell_Axi<Complexe>::GetNbPhysicalIndices() const
  {
    return ref_epsilon.GetM();
  }
  
  
  //! returns true if physical properties of media are varying
  template<class Complexe>
  inline bool HarmonicMaxwell_Axi<Complexe>::IsVaryingMedia(int i) const
  {
    if (ref_epsilon(i).IsVarying()||ref_mu(i).IsVarying()||ref_sigma(i).IsVarying())
      return true;
    
    return false;
  }
  
  
  //! returns the velocity of waves in physical media whose reference is ref
  template<class Complexe>
  inline Real_wp HarmonicMaxwell_Axi<Complexe>::GetVelocityOfMedia(int ref) const
  {
    return 1.0/sqrt(abs(this->ref_epsilon(ref)(0,0) * this->ref_mu(ref)(0, 0)));
  }

  
  //! returns the velocity of waves at infinity
  template<class Complexe>
  inline Real_wp HarmonicMaxwell_Axi<Complexe>::GetVelocityOfInfinity() const
  {
    return 1.0/sqrt(abs(this->epsilon0 * this->mu0));
  }


  template<>
  inline void HarmonicMaxwell_Axi<Real_wp>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
                      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
                      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                      Vector<bool>& compute_hess)
  {
    GetVaryingIndex(rho_real, num_ref, num_index, num_component, compute_grad, compute_hess);
  }


  template<>
  inline void HarmonicMaxwell_Axi<Complex_wp>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
                      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
                      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                      Vector<bool>& compute_hess)
  {
    GetVaryingIndex(rho_complex, num_ref, num_index, num_component, compute_grad, compute_hess);
  }
  

  /*************************
   * HarmonicMaxwellAxi_Eq *
   *************************/


  template<class TypeEquation>
  inline HarmonicMaxwellAxi_Eq<TypeEquation>::HarmonicMaxwellAxi_Eq()
    : HarmonicMaxwell_Axi<typename TypeEquation::Complexe>(this->GetLeafClass())
  {
    this->list_number_mode.Reallocate(2);
    this->list_number_mode(0) = -1; this->list_number_mode(1) = 1;
    
    this->threshold_rhs = epsilon_machine;
    this->do_not_store_modes = true;
  }
  
  
  template<class TypeEquation>
  inline void HarmonicMaxwellAxi_Eq<TypeEquation>::SetWaveVector(const R3& k)
  {
    this->kwave3D = k; 
  }
  

  template<class TypeEquation>
  inline const R3& HarmonicMaxwellAxi_Eq<TypeEquation>::GetWaveVector() const
  {
    return this->kwave3D; 
  }
  

  template<class TypeEquation>  
  inline const R3& HarmonicMaxwellAxi_Eq<TypeEquation>::GetPhaseOrigin() const
  {
    return HarmonicMaxwell_Axi<Complexe>::GetPhaseOrigin();
  }

  
  template<class TypeEquation>
  inline int HarmonicMaxwellAxi_Eq<TypeEquation>::GetNbModesSource() const
  {
    return 1;
  }  
  
  
  template<class TypeEquation>
  inline bool HarmonicMaxwellAxi_Eq<TypeEquation>::UseNumericalIntegration(int iquad) const
  {
    return true;
  }

  
  //! returns true if the problem is symmetric
  template<class TypeEquation>
  inline bool HarmonicMaxwellAxi_Eq<TypeEquation>::IsSymmetricProblem(bool eigen) const
  {
    if (this->FirstOrderFormulation())
      return false;
    
    return VarHarmonic<TypeEquation>::IsSymmetricProblem(eigen);
  }
  
  
  template<class TypeEquation>
  inline int HarmonicMaxwellAxi_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return HarmonicMaxwell_Axi<Complexe>::GetNbPhysicalIndices();
  }
  
  
  template<class TypeEquation>
  inline bool HarmonicMaxwellAxi_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return HarmonicMaxwell_Axi<Complexe>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  inline bool HarmonicMaxwellAxi_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return HarmonicMaxwell_Axi<Complexe>::IsVaryingMedia(i);
  }

  
  template<class TypeEquation>
  inline Real_wp HarmonicMaxwellAxi_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return HarmonicMaxwell_Axi<Complexe>::GetVelocityOfMedia(ref);
  }
  
  
  template<class TypeEquation>
  inline Real_wp HarmonicMaxwellAxi_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return HarmonicMaxwell_Axi<Complexe>::GetVelocityOfInfinity();
  }
        
  
  template<class TypeEquation>
  inline void HarmonicMaxwellAxi_Eq<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
    HarmonicMaxwell_Axi<Complexe>::SetInputData(description_field, parameters);
  }
    

  template<class TypeEquation>
  inline void HarmonicMaxwellAxi_Eq<TypeEquation>::CheckInputMesh()
  {
    HarmonicMaxwell_Axi<Complexe>::CheckSectionMeshAxi();
  }
  
    
  template<class TypeEquation>
  inline void HarmonicMaxwellAxi_Eq<TypeEquation>::InitIndices(int n)
  {
    HarmonicMaxwell_Axi<Complexe>::InitIndices(n);
  }
  
    
  template<class TypeEquation>
  inline void HarmonicMaxwellAxi_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    HarmonicMaxwell_Axi<Complexe>::SetIndices(i, parameters);
  }
  
  
  template<class TypeEquation>
  inline void HarmonicMaxwellAxi_Eq<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    HarmonicMaxwell_Axi<Complexe>::SetPhysicalIndex(name_media, i, parameters);
  }

  
  template<class TypeEquation>
  inline string HarmonicMaxwellAxi_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return HarmonicMaxwell_Axi<Complexe>::GetPhysicalIndexName(m);
  }

  
  template<class TypeEquation>
  inline void HarmonicMaxwellAxi_Eq<TypeEquation>::FinalizeComputationVaryingIndices()
  {
    HarmonicMaxwell_Axi<Complexe>::FinalizeComputationVaryingIndices();
  }
  
  
  template<class TypeEquation>
  inline void HarmonicMaxwellAxi_Eq<TypeEquation>::AllocateMassMatrices()
  {
    HarmonicMaxwell_Axi<Complexe>::AllocateMassMatrices();
  }
  
  
  template<class TypeEquation>
  inline void HarmonicMaxwellAxi_Eq<TypeEquation>
  ::ComputeLocalMassMatrix(int i, int N, bool linear_sparse,
                           SetPoints<Dimension2>& PointsElem,
                           SetMatrices<Dimension2>& MatricesElem,
                           IVect& OrderFace, const ElementGeomReference<Dimension2>& Fb)
  {
    VarHarmonic<TypeEquation>::ComputeLocalMassMatrix(i, N, linear_sparse, PointsElem,
                                                      MatricesElem, OrderFace, Fb);
    
    HarmonicMaxwell_Axi<Complexe>::ComputeLocalMassMatrix(i);
  }
  

  template<class TypeEquation>
  inline void HarmonicMaxwellAxi_Eq<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
                      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
                      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                      Vector<bool>& compute_hess)
  {
    HarmonicMaxwell_Axi<Complexe>::GetVaryingIndices(rho_complex, rho_real, num_ref, num_index,
						     num_component, compute_grad, compute_hess);
  }
  

  /************************************
   * HarmonicMaxwellEquation_HcurlAxi *
   ************************************/


  inline bool HarmonicMaxwellEquation_HcurlAxi::SymmetricGlobalMatrix()
  {
    return true; 
  }
  
  
  inline bool HarmonicMaxwellEquation_HcurlAxi::SymmetricElementaryMatrix()
  {
    return false;
  }
  
  
  template<class TypeEquation>
  inline void HarmonicMaxwellEquation_HcurlAxi
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
		      int i, const ElementReference_Dim<Dimension2>& Fb)
  {
    var.ComputeHcurlMassMatrix(i, Fb);
  }
  
  
  //! time-harmonic -> -omega
  inline Real_wp HarmonicMaxwellEquation_HcurlAxi
  ::GetCoefficientMassMatrix(const VarProblem_Base& var)
  {
    return -var.GetOmega();
  }
  
  
  //! time-harmonic -> 1/omega
  inline Real_wp HarmonicMaxwellEquation_HcurlAxi
  ::GetCoefficientStiffnessMatrix(const VarProblem_Base& vars)
  {
    return 1.0/vars.GetOmega();
  }
  
  
  //! not used
  inline Real_wp HarmonicMaxwellEquation_HcurlAxi::GetCoefMode(const VarProblem_Base& vars)
  {
    return 1.0;
  }


  template<class TypeEquation>
  inline ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complex_wp, Dimension2>(var),
      var_maxwell(var)
  {
  }

  inline bool ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>::PresenceGradient() const
  {    
    if ((modified_formulation) && (number_mode != 0))
      return true;

    return false;
  }
  
  
  template<class TypeEquation>
  inline ImpedanceGeneric<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>
  ::ImpedanceGeneric(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complex_wp, Dimension2>(var)
  {
  }
  
  
  template<class TypeEquation>
  inline VarComputationRCS_MaxwellAxi::VarComputationRCS_MaxwellAxi(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_Axi(var), var_maxwell(var)
  {
  }
  

  template<class TypeEquation>
  inline VarComputationRCS<HarmonicMaxwellEquation_HcurlAxi>::VarComputationRCS(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_MaxwellAxi(var)
  {
  }

}

#define MONTJOIE_FILE_MAXWELL_AXISYM_HARMONIC_INLINE_CXX
#endif
