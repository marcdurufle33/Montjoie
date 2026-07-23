#ifndef MONTJOIE_FILE_VAR_HARMONIC_BASE_INLINE_CXX

namespace Montjoie
{
  //! constructor of VarHarmonicClass, set default values
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline VarHarmonic_Base<Complexe, Dimension>::VarHarmonic_Base(EllipticProblem<TypeEquation>& var)
    : VarProblem<Dimension>(var),
      VarBoundaryCondition<Complexe, Dimension>(var),
      VarOutputProblem<Dimension>(var),
      VarSourceProblem<Complexe, Dimension>(var)
    , var_gibc_base(var.var_gibc)

#ifdef MONTJOIE_WITH_TRANSMISSION
    , var_transmission_base(var.var_transmission)
#endif

    , output_rcs_base(var.output_rcs_param)
  {
  }


  //! true if the problem is solved in complex numbers
  template<class Complexe, class Dimension>
  inline bool VarHarmonic_Base<Complexe, Dimension>::IsComplexProblem() const
  {
    Complexe z(0);
    return IsComplexNumber(z);
  }


  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>
  ::ComputeMassMatrix(bool compute_rho, bool delete_points)
  {
    VarProblem<Dimension>::ComputeMassMatrix(compute_rho, delete_points);
  }
  

  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>
  ::ComputeStoreCoefficientsPML(int i1, int num, const VectR_N& AllPoints)
  {
    VarBoundaryCondition_Dim<Complexe, Dimension>::ComputeStoreCoefficientsPML(i1, num, AllPoints);
  }
  

  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>
  ::SetComputationFarPoints(VectR_N& Points, Real_wp dt)
  {
    this->output_rcs_base.SetOutsidePoints(Points);
    this->output_rcs_base.SetTimeStep(dt);
  }

  
  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>::InitVarGrid()
  {
    VarOutputProblem<Dimension>::InitVarGrid();
    
    if (output_rcs_base.GetNbPointsOutside() > 0)
      output_rcs_base.InitComputationRCS();
  }
  
  
  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>::FindElementsInsidePML()
  {
    VarBoundaryCondition<Complexe, Dimension>::FindElementsInsidePML();
  }
  

  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>::TreatGibc(const IVect& Epart)
  {
    var_gibc_base.TreatGibc(Epart); 
  }
  
  
  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>::InitGibcReferences(int N)
  {
    var_gibc_base.InitReferences(N); 
  }


#ifdef MONTJOIE_WITH_TRANSMISSION
  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>::PartMeshTransmission()
  {
    var_transmission_base.PartMeshTransmission();
  }
  
  
  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>::TreatTransmission(const IVect& epart)
  {
    var_transmission_base.TreatTransmission(epart);
  }
  
  
  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>
  ::SendTransmissionDofs(const IVect& num, int& nb0, int& nb1, IVect& info0, IVect& info1) 
  {
    var_transmission_base.SendTransmissionDofs(num, nb0, nb1, info0, info1);
  }

  
  template<class Complexe, class Dimension>
  inline void VarHarmonic_Base<Complexe, Dimension>
  ::DistributeTransmissionDofs(int nb0, int nb1, IVect& info0, IVect& info1) 
  { 
    var_transmission_base.DistributeTransmissionDofs(nb0, nb1, info0, info1);
  }
#endif

  
}

#define MONTJOIE_FILE_VAR_HARMONIC_BASE_INLINE_CXX
#endif

