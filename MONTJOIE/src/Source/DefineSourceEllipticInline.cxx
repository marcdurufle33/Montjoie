#ifndef MONTJOIE_FILE_DEFINE_SOURCE_ELLIPTIC_INLINE_CXX

namespace Montjoie
{
  
  /*************************
   * VarSourceProblem_Base *
   *************************/
  
  
  //! default constructor
  template<class TypeEquation>
  inline VarSourceProblem_Base::VarSourceProblem_Base(EllipticProblem<TypeEquation>& var)
    : var_computation(var), var_boundary(var), var_problem(var)
  {
    InitDefaultValues();
  }
  
  
  //! returns true if only one right hand side needs to be computed
  inline bool VarSourceProblem_Base::OnlyOneSource() const
  {
    return (GetNbRhs() == 1);
  }

  
  //! returns the number of right hand sides
  inline int VarSourceProblem_Base::GetNbRhs() const
  {
    return type_source.GetM();
  }


  //! returns the type of the source n
  inline int VarSourceProblem_Base::GetSourceType(int n) const
  {
    return type_source(n);
  }


  //! returns the threshold used to neglect some modes
  inline Real_wp VarSourceProblem_Base::GetThresholdSource() const
  {
    return threshold_rhs;
  }


  //! returns the type of incident field
  inline int VarSourceProblem_Base::GetIncidentFieldType(int n) const
  {
    return incident_source(n);
  }
  
  
  //! returns the number of lines for source n
  inline int VarSourceProblem_Base::GetNbParameterSource(int n) const
  {
    return source_space_param(n).GetM();
  }
  

  //! returns the parameters for the line num of source n
  inline const VectString& VarSourceProblem_Base::GetParameterSource(int n, int num) const
  {
    return source_space_param(n)(num);
  }
  
  
  inline void VarSourceProblem_Base::GetMemoryUsed(map<string, size_t>& var) const
  {
  }

  
  inline void VarSourceProblem_Base
  ::ModifyVolumetricSource(int i, int j, const R2&,
			   const VirtualSourceField<Real_wp, Dimension2>&,
			   Vector<Real_wp>&) const 
  {
  }
  
  inline void VarSourceProblem_Base
  ::ModifyVolumetricSource(int i, int j, const R2&,
			   const VirtualSourceField<Complex_wp, Dimension2>&,
			   Vector<Complex_wp>&) const 
  {
  }
  
  inline void VarSourceProblem_Base
  ::ModifyVolumetricSource(int i, int j, const R3&,
			   const VirtualSourceField<Real_wp, Dimension3>&,
			   Vector<Real_wp>&) const 
  {
  }
  
  inline void VarSourceProblem_Base
  ::ModifyVolumetricSource(int i, int j, const R3&,
			   const VirtualSourceField<Complex_wp, Dimension3>&,
			   Vector<Complex_wp>&) const 
  {
  }
  
  
  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline VarSourceProblem_Dim<Dimension>::VarSourceProblem_Dim(EllipticProblem<TypeEquation>& var)
    : var_problem(var), var_source(var)
  {
  }


  /*************************
   * VarSourceProblem_Cplx *
   *************************/


  //! default constructor
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline VarSourceProblem_Cplx<Complexe, Dimension>
  ::VarSourceProblem_Cplx(EllipticProblem<TypeEquation>& var)
    : VarSourceProblem_Dim<Dimension>(var), var_problem(var), var_source(var)
  {
    SetComplexOne(coefficient_volumic_source);
  }

  
  //! returns coefficient for volume source
  template<class Complexe, class Dimension>
  inline Complexe VarSourceProblem_Cplx<Complexe, Dimension>
  ::GetCoefficientVolumeSource() const
  {
    return coefficient_volumic_source;
  }


  template<class Complexe, class Dimension> template<class T>
  inline IncidentWaveField<T, Dimension>*
  VarSourceProblem_Cplx<Complexe, Dimension>::GetIncidentFieldGen(int n, const T& val) const
  {
    cout << "Not implemented for any number" << endl;
    abort();
    return NULL;
  }
    
  
  template<class Complexe, class Dimension> template<class T>
  inline IncidentWaveProjector<T, Dimension>* 
  VarSourceProblem_Cplx<Complexe, Dimension>::GetIncidentWaveProjectorGen(int n, const T&) const
  {
    cout << "Not implemented for any number" << endl;
    abort();
    return NULL;
  }

  
  template<class Complexe, class Dimension>
  inline IncidentWaveField<Real_wp, Dimension>* VarSourceProblem_Cplx<Complexe, Dimension>
  ::GetIncidentField(int n, const Real_wp& val) const
  {
    return GetIncidentFieldGen(n, val);
  }
  
  
  template<class Complexe, class Dimension>
  inline IncidentWaveProjector<Real_wp, Dimension>* VarSourceProblem_Cplx<Complexe, Dimension>
  ::GetIncidentWaveProjector(int n, const Real_wp& val) const
  {
    return GetIncidentWaveProjectorGen(n, val);
  }

  
  template<class Complexe, class Dimension>
  inline IncidentWaveField<Complex_wp, Dimension>* VarSourceProblem_Cplx<Complexe, Dimension>
  ::GetIncidentField(int n, const Complex_wp& val) const
  {
    return GetIncidentFieldGen(n, val);
  }
  
  
  template<class Complexe, class Dimension>
  inline IncidentWaveProjector<Complex_wp, Dimension>* VarSourceProblem_Cplx<Complexe, Dimension>
  ::GetIncidentWaveProjector(int n, const Complex_wp& val) const
  {
    return GetIncidentWaveProjectorGen(n, val);
  }


  /************************
   * VarSourceProblem_Fem *
   ************************/

  
  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline VarSourceProblem_Fem<Dimension>::VarSourceProblem_Fem(EllipticProblem<TypeEquation>& var)
    : VarSourceProblem_Base(var), var_problem(var), var_boundary(var), var_source(var)
  {
  }
  

  template<class Dimension>
  inline void VarSourceProblem_Fem<Dimension>
  ::AddDiracSource(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
		   Vector<VirtualSourceFEM<Real_wp, Dimension>* >& f) const
  {
    AddDiracSourceGen(alpha, b_src, f);
  }
    

  template<class Dimension>  
  inline void VarSourceProblem_Fem<Dimension>
  ::AddDiracSource(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
		   Vector<VirtualSourceFEM<Complex_wp, Dimension>* >& f) const
  {
    AddDiracSourceGen(alpha, b_src, f);
  }
  

  /********************
   * VarSourceProblem *
   ********************/

  
  //! default constructor
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline VarSourceProblem<Complexe, Dimension>
  ::VarSourceProblem(EllipticProblem<TypeEquation>& var) : VarSourceProblem_Fem<Dimension>(var),
							   VarSourceProblem_Cplx<Complexe, Dimension>(var),
							   var_problem(var)
  {
  }


  template<class Complexe, class Dimension> template<class T>
  inline void VarSourceProblem<Complexe, Dimension>
  ::AddIncidentWaveGen(const T&, Vector<Vector<T> >&)
  {
    cout << "Not implemented for any type of number" << endl;
    abort();
  }

  
  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::AddIncidentWave(const Real_wp& alpha, Vector<VectReal_wp>& U)
  {
    AddIncidentWaveGen(alpha, U);
  }
  

  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::AddIncidentWave(const Complex_wp& alpha, Vector<VectComplex_wp>& U)
  {
    AddIncidentWaveGen(alpha, U);
  }


  template<class Complexe, class Dimension> template<class T>
  inline void VarSourceProblem<Complexe, Dimension>
  ::ComputeRightHandSideGen(Vector<T>&, bool assemble)
  {
    cout << "Not implemented for any type of number" << endl;
    abort();
  }

  
  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::ComputeRightHandSide(VectReal_wp& b_source, bool assemble)
  {
    ComputeRightHandSideGen(b_source, assemble);
  }   
  

  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::ComputeRightHandSide(VectComplex_wp& b_source, bool assemble)  
  {
    ComputeRightHandSideGen(b_source, assemble);
  }
  

  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::AddVolumeProjection(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
			Vector<VirtualProjectorFEM<Real_wp, Dimension>* >& f) const
  {
    this->AddVolumetricProjection(alpha, b_src, f);
  }


  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::AddVolumeProjection(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
			Vector<VirtualProjectorFEM<Complex_wp, Dimension>* >& f) const
  {
    this->AddVolumetricProjection(alpha, b_src, f);
  }


  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::AddSurfaceSource(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
		     Vector<VirtualSourceFEM<Real_wp, Dimension>* >& f) const
  {
    this->AddSurfacicSource(alpha, b_src, f);
  }


  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::AddSurfaceSource(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
		     Vector<VirtualSourceFEM<Complex_wp, Dimension>* >& f) const
  {
    this->AddSurfacicSource(alpha, b_src, f);
  }


  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::SetSurfaceProjection(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
			 Vector<VirtualSourceFEM<Real_wp, Dimension>* >& f) const
  {
    this->SetSurfacicProjection(alpha, b_src, f);
  }


  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::SetSurfaceProjection(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
			 Vector<VirtualSourceFEM<Complex_wp, Dimension>* >& f) const
  {
    this->SetSurfacicProjection(alpha, b_src, f);
  }


  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::AddVolumeSource(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
		    Vector<VirtualSourceFEM<Real_wp, Dimension>* >& f) const
  {
    this->AddVolumetricSource(alpha, b_src, f);
  }


  template<class Complexe, class Dimension>
  inline void VarSourceProblem<Complexe, Dimension>
  ::AddVolumeSource(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
		    Vector<VirtualSourceFEM<Complex_wp, Dimension>* >& f) const
  {
    this->AddVolumetricSource(alpha, b_src, f);
  }
    
} // end namespace

#define MONTJOIE_FILE_DEFINE_SOURCE_ELLIPTIC_INLINE_CXX
#endif  
