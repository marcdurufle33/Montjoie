#ifndef MONTJOIE_FILE_VAR_LAPLACE_INLINE_CXX

namespace Montjoie
{

  /*******************
   * VarLaplace_Base *
   *******************/
  

  template<class TypeEquation>
  inline void VarLaplace_Base<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
                      Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real, IVect& num_ref,
                      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                      Vector<bool>& compute_hess)
  {
    VarHelmholtz_Base<Complexe, Dimension>::GetVaryingIndices(rho_real, num_ref, num_index,
                                                              num_component, compute_grad, compute_hess);
  }
  
  
  /*********************
   * VarComputationRCS *
   *********************/


  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline VarComputationRCS_Laplace<Dimension>
  ::VarComputationRCS_Laplace(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_Base<Dimension>(var), var_problem(var), var_output(var)
  {
    InitDefaultValues();
  }

  
  template<class Dimension>
  inline void VarComputationRCS_Laplace<Dimension>
  ::ComputeRCS(const VectReal_wp& )
  {
  }
  
  
  template<class Dimension>
  inline void VarComputationRCS_Laplace<Dimension>
  ::ComputeRCS(const VectComplex_wp& ) 
  {
  }


  template<class Dimension>
  inline const string& VarComputationRCS_Laplace<Dimension>
  ::GetFileNameNormL2() const
  {
    return file_name_normL2; 
  }
  
  
  template<class Dimension>
  inline void VarComputationRCS_Laplace<Dimension>
  ::SetFileNameNormL2(const string& nom)
  {
    file_name_normL2 = nom; 
  }
  
  
  template<class Dimension>
  inline void VarComputationRCS_Laplace<Dimension>
  ::SetL2NormOutput(bool compute_L2_norm_, const Real_wp& t_begin_calculL2_,
		    const string& file_name_normL2_)
  {
    compute_L2_norm = compute_L2_norm_;
    t_begin_calculL2 = t_begin_calculL2_;
    file_name_normL2 = file_name_normL2_;
  }
  
  
  template<class Dimension> template<class TypeEquation>
  inline VarComputationRCS<LaplaceEquation<Dimension> >
  ::VarComputationRCS(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_Laplace<Dimension>(var)
  {
  }
  

#ifdef MONTJOIE_WITH_TRANSMISSION
  template<class TypeEquation>
  inline VarTransmission<LaplaceEquation<Dimension2> > 
  ::VarTransmission(EllipticProblem<TypeEquation>& var)
    : VarTransmission_Helm<Real_wp, Dimension2>(var) 
  {
  }
#endif
  

  template<class T, class Dimension>
  inline CondensationBlockSolver<T, LaplaceEquation<Dimension> >
  ::CondensationBlockSolver(EllipticProblem<LaplaceEquation<Dimension> >& var)
    : CondensationBlockSolver_Helm<T, LaplaceEquation<Dimension> >(var)
  {
  }


  template<class T, class Dimension>
  inline CondensationBlockSolver<T, LaplaceEquationDG<Dimension> >
  ::CondensationBlockSolver(EllipticProblem<LaplaceEquationDG<Dimension> >& var)
    : CondensationBlockSolver_HelmDG<T, Real_wp, Dimension>(var)
  {
  }

  /********************
   * Plane wave stuff *
   ********************/


  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline IncidentWaveProjector_Laplace<Dimension>
  ::IncidentWaveProjector_Laplace(const EllipticProblem<TypeEquation>& var,
				  IncidentWaveField<Real_wp, Dimension>& u_inc)
    : IncidentWaveProjector<Real_wp, Dimension>(var, u_inc), var_source(var),
      var_boundary(var), var_problem(var)     
  {
  }    
  

  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline DiffractedWaveSource_Laplace<Dimension>
  ::DiffractedWaveSource_Laplace(const EllipticProblem<TypeEquation>& var,
				 IncidentWaveField<Real_wp, Dimension>& u_inc)
    : VirtualSourceFEM<Real_wp, Dimension>(var), incident_wave(u_inc), var_source(var),
      var_boundary(var), var_problem(var), var_laplace(var)
  {
    scalar_eq = true;
    dirichlet_cond = true;
  }    
  

  //! constructor with given problem
  template<class Dimension> template<class TypeEquation>
  inline DiffractedWaveSource<LaplaceEquation<Dimension> >
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Real_wp, Dimension>& fsrc_)
    : DiffractedWaveSource_Laplace<Dimension>(vars_helm, fsrc_) 
  {
  }
  

  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline TotalWaveSource_Laplace<Dimension>
  ::TotalWaveSource_Laplace(const EllipticProblem<TypeEquation>& var,
			    IncidentWaveField<Real_wp, Dimension>& u_inc)
    : VirtualSourceFEM<Real_wp, Dimension>(var), incident_wave(u_inc), var_source(var),
      var_boundary(var), var_problem(var), var_laplace(var)
  {
    scalar_eq = true;
  }    
  

  //! constructor with given problem
  template<class Dimension> template<class TypeEquation>
  inline TotalWaveSource<LaplaceEquation<Dimension> >
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Real_wp, Dimension>& fsrc)
    : TotalWaveSource_Laplace<Dimension>(vars_helm, fsrc) 
  {
  }


  //! constructor with given problem
  template<class Dimension>
  VolumetricSource<LaplaceEquationDG<Dimension> >
  ::VolumetricSource(const EllipticProblem<LaplaceEquationDG<Dimension> >& var,
		     const Vector<VectString>& param)
    : VolumetricSource_HelmDG<Real_wp, Dimension>(var, param)
  {    
  }

  
  template<class Dimension> 
  inline IncidentWaveProjector<Real_wp, Dimension>*
  EllipticProblem<LaplaceEquation<Dimension> >
  ::GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Real_wp, Dimension>& u_inc) const
  {
    IncidentWaveProjector_Laplace<Dimension>* proj;
    proj = new IncidentWaveProjector_Laplace<Dimension>(*this, u_inc);
    proj->Init(0.0);
    return proj;
  }

  template<class Dimension>
  inline void EllipticProblem<LaplaceEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelm(i, num_dof, mat_elem, nat_mat, *this, *this,
					this->GetReferenceElementH1(i)); 
  }


  template<class Dimension>
  inline void EllipticProblem<LaplaceEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelm(i, num_dof, mat_elem, nat_mat, *this, *this,
					this->GetReferenceElementH1(i)); 
  }  
  
  template<class Dimension>
  inline void EllipticProblem<LaplaceEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);    
  }


  template<class Dimension>
  inline void EllipticProblem<LaplaceEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }

  template<class Dimension>
  inline FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<LaplaceEquation<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Real_wp, LaplaceEquation<Dimension> >(var_)
  {
  };


  template<class Dimension>
  inline FemMatrixFreeClass<Complex_wp, LaplaceEquation<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<LaplaceEquation<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Complex_wp, LaplaceEquation<Dimension> >(var_)
  {
  };

   
  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline IncidentWaveProjector_LaplaceDG<Dimension>
  ::IncidentWaveProjector_LaplaceDG(const EllipticProblem<TypeEquation>& var,
				    IncidentWaveField<Real_wp, Dimension>& u_inc)
    : IncidentWaveProjector_Laplace<Dimension>(var, u_inc)
  {
  }    
  

  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline DiffractedWaveSource_LaplaceDG<Dimension>
  ::DiffractedWaveSource_LaplaceDG(const EllipticProblem<TypeEquation>& var,
				   IncidentWaveField<Real_wp, Dimension>& u_inc)
    : VirtualSourceFEM<Real_wp, Dimension>(var), incident_wave(u_inc), var_source(var),
      var_boundary(var), var_problem(var), var_laplace(var)
  {
    scalar_eq = true;
  }    
  
  
  //! constructor with given problem
  template<class Dimension> template<class TypeEquation>
  inline DiffractedWaveSource<LaplaceEquationDG<Dimension> >
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Real_wp, Dimension>& fsrc_)
    : DiffractedWaveSource_LaplaceDG<Dimension>(vars_helm, fsrc_) 
  {
  }


  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline TotalWaveSource_LaplaceDG<Dimension>
  ::TotalWaveSource_LaplaceDG(const EllipticProblem<TypeEquation>& var,
			      IncidentWaveField<Real_wp, Dimension>& u_inc)
    : VirtualSourceFEM<Real_wp, Dimension>(var), incident_wave(u_inc), var_source(var),
      var_boundary(var), var_problem(var), var_laplace(var)
  {
    scalar_eq = true;
  }    
  

  //! constructor with given problem
  template<class Dimension> template<class TypeEquation>
  inline TotalWaveSource<LaplaceEquationDG<Dimension> >
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Real_wp, Dimension>& fsrc)
    : TotalWaveSource_LaplaceDG<Dimension>(vars_helm, fsrc) 
  {
  }


  template<class Dimension>
  inline IncidentWaveProjector<Real_wp, Dimension>*
  EllipticProblem<LaplaceEquationDG<Dimension> >
  ::GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Real_wp, Dimension>& u_inc) const
  {
    IncidentWaveProjector_LaplaceDG<Dimension>* proj;
    proj = new IncidentWaveProjector_LaplaceDG<Dimension>(*this, u_inc);
    proj->Init(0.0);
    return proj;
  }


  template<class Dimension>
  inline void EllipticProblem<LaplaceEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>& solver,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelmDG(i, num_dof, mat_elem, solver, nat_mat, *this, *this,
                                        this->GetReferenceElementH1(i)); 
  }


  template<class Dimension>
  inline void EllipticProblem<LaplaceEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>& solver,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelmDG(i, num_dof, mat_elem, solver, nat_mat, *this, *this,
					this->GetReferenceElementH1(i)); 
  }  

  template<class Dimension>
  inline void EllipticProblem<LaplaceEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    this->AddElementaryFluxesHelmDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }


  template<class Dimension>
  inline void EllipticProblem<LaplaceEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    this->AddElementaryFluxesHelmDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }
  

  template<class Dimension>
  inline FemMatrixFreeClass<Real_wp, LaplaceEquationDG<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<LaplaceEquationDG<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Real_wp, LaplaceEquationDG<Dimension> >(var_)
  {
  };
    
  
  template<class Dimension>
  inline FemMatrixFreeClass<Complex_wp, LaplaceEquationDG<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<LaplaceEquationDG<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Complex_wp, LaplaceEquationDG<Dimension> >(var_)
  {
  };
  
   
  /*************
   * Impedance *
   *************/


  template<class Complexe, class Dimension> template<class TypeEquation>
  inline ImpedanceABC_Laplace<Complexe, Dimension>
  ::ImpedanceABC_Laplace(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complexe, Dimension>(var),  var_helm(var), var_problem(var)
  {
  }
  
  
  template<class Complexe, class Dimension>
  inline ImpedanceABC<Complexe, LaplaceEquation<Dimension> >
  ::ImpedanceABC(const EllipticProblem<LaplaceEquation<Dimension> >& var)
    : ImpedanceABC_Laplace<Complexe, Dimension>(var)
  {
  }


  template<class Complexe, class Dimension>
  inline ImpedanceABC<Complexe, LaplaceEquationDG<Dimension> >
  ::ImpedanceABC(const EllipticProblem<LaplaceEquationDG<Dimension> >& var)
    : ImpedanceABC_Laplace<Complexe, Dimension>(var)
  {
  }
  
}

#define MONTJOIE_FILE_VAR_LAPLACE_INLINE_CXX
#endif
