#ifndef MONTJOIE_FILE_VAR_INSTATIONARY_INLINE_CXX

namespace Montjoie
{
  //! default constructor
  template<class TypeEquation>
  inline VarInstationary_Base::VarInstationary_Base(HyperbolicProblem<TypeEquation>& var)
    : dim_N(TypeEquation::TypeEquationStationary::Dimension::dim_N),
      var_problem(var.var_harmonic), var_computation(var.var_harmonic),
      var_comm(var.var_harmonic), var_output(var.var_harmonic), var_boundary(var.var_harmonic),
      var_source(var.var_harmonic)
  {
    InitDefaultValues();
  }


  inline Real_wp VarInstationary_Base::GetFinalTimeSource() const
  {
    return tlimit_source;
  }
  
  
  inline void VarInstationary_Base::SetFinalTimeSource(const Real_wp& t)
  {
    tlimit_source = t;
  }
  
  
  inline Vector<Real_wp, VectSparse>& VarInstationary_Base::GetSpaceSource()
  {
    return sparse_vector_source;
  }
  
  
  //! returns true if the system can be split with two unknowns (like E and H for Maxwell)
  inline bool VarInstationary_Base::SplitSystem() const
  {
    return false; 
  }
  

  inline Real_wp VarInstationary_Base::GetInitialTime() const
  {
    return initial_time; 
  }
  
  
  inline Real_wp VarInstationary_Base::GetFinalTime() const
  {
    return final_time;
  }
  
  
  //! is the time scheme implicit ?
  inline bool VarInstationary_Base::IsImplicitScheme() const
  {
    return var_time_scheme.IsImplicitScheme();
  }
  
  
  //! returns true if the time scheme is a scheme adapted for first-order formulation
  inline bool VarInstationary_Base::FirstOrderScheme() const
  {
    return var_time_scheme.FirstOrderScheme();
  }
  
  
  //! returns the time scheme used in the simulation
  inline int VarInstationary_Base::GetTimeSchemeType() const
  {
    return var_time_scheme.GetTimeSchemeType();
  }
  
  
  //! returns the order of used time scheme
  inline int VarInstationary_Base::GetTimeSchemeOrder() const
  {
    return var_time_scheme.GetOrder();
  }
  
  
  //! returns the time step
  inline Real_wp VarInstationary_Base::GetTimeStep() const
  {
    return deltat; 
  }


  //! sets the time step
  inline void VarInstationary_Base::SetTimeStep(const Real_wp& dt)
  {
    deltat = dt;
  }
  
  
  //! returns the object containing different levels for local time stepping
  inline MatrixVectorProductLevel& VarInstationary_Base::GetTimeLevelDistribution()
  {
    return level_time_scheme;
  }


  //! returns the object containing different levels for local time stepping
  inline const MatrixVectorProductLevel& VarInstationary_Base::GetTimeLevelDistribution() const
  {
    return level_time_scheme;
  }
  

  inline int VarInstationary_Base::GetOperatorLevel() const
  {
    return global_level_operator; 
  }
  

  inline void VarInstationary_Base::ForceComputationTimeStep(bool to_be_computed)
  {
    time_step_to_be_computed = to_be_computed;
  }


  inline int VarInstationary_Base::GetCFL_ComputationalMethod() const
  {
    return method_computation_cfl;
  }
  
  
  inline const string& VarInstationary_Base::GetFileStoringLocalDt() const
  {
    return file_with_local_dt;
  }
  

  inline void VarInstationary_Base::SetLevel(int level)
  {
    global_level_operator = level; 
  }
  
  
  inline void VarInstationary_Base::ComputeExtrapolationVectorial(const VectReal_wp&, const VectReal_wp&)
  {
  }
  
  
  inline void VarInstationary_Base::ComputeExtrapolationScalar(const VectReal_wp&, const VectReal_wp&)
  {
  }    
  
  
  //! returns the number of unknowns of the global system
  inline int VarInstationary_Base::GetNumberOfUnknowns() const
  {
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      return var_problem.GetNbDof() - var_problem.nb_unknowns_hdg*var_problem.GetNbMeshDof();
    
    return var_problem.GetNbDof();
  }
  
  
  //! returns the size of the scalar solution
  inline int VarInstationary_Base::GetNbScalarUnknowns() const
  {
    return var_comm.GetNbMainUnknownDof()*var_problem.nb_unknowns_scal;
  }
  
  
  //! applies the stiffness matrix 
  /*!
    This method computes 
    Yh = beta Yh + alpha K Uh
    where K is the stiffness matrix involved in second-order formulation :
    M d^2 / dt^2 (U) + S dU/dt + K U = 0
  */
  inline void VarInstationary_Base
  ::ApplyOperatorKh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
		    const Real_wp& beta, VectReal_wp& Yh)
  {
    Glob_mat_Kh->MltAddVector(-alpha, Uh, beta, Yh, false);
  }
  
  
  //! applies the mass matrix Dh, Prod_Uh = beta Prod_Uh + alpha Dh Uh
  inline void VarInstationary_Base
  ::ApplyOperatorDh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
		    const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    Glob_mat_Dh->ApplyOperatorDh(alpha, tn, Uh, beta, Prod_Uh);
  }
  
  
  //! applies the damping matrix Sh, Prod_Uh = beta Prod_Uh + alpha Sh Uh
  inline void VarInstationary_Base
  ::ApplyOperatorSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
		    const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    Glob_mat_Dh->ApplyOperatorSh(alpha, tn, Uh, beta, Prod_Uh);
  }
  
  
  //! applies the damping matrix ShV, Prod_Uh = beta Prod_Uh + alpha ShV Uh
  inline void VarInstationary_Base::
  ApplyOperatorShVectorial(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			   const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    Glob_mat_Dh->ApplyOperatorShVectorial(alpha, tn, Uh, beta, Prod_Uh);
  }
  
  
  //! applies the matrix Dh - dt/2 Sh, Prod_Uh = beta Prod_Uh + alpha (Dh - dt/2 Sh) Uh
  inline void VarInstationary_Base::
  ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			   const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    Glob_mat_Dh->ApplyOperatorDhMinusdtSh(alpha, tn, Uh, beta, Prod_Uh);
  }
  
  
  //! solves the matrix Dh + dt/2 Sh, Y = (Dh + dt/2 Sh)^{-1} Y
  inline void VarInstationary_Base::
  SolveOperatorDhPlusdtSh(VectReal_wp& Y)
  {
    this->Assemble(Y);
    Glob_mat_Dh->SolveOperatorDhPlusdtSh(Y);
  }
  
  
  //! solves the matrix Dh + dt/2 Sh, Y = Y + alpha (Dh + dt/2 Sh)^{-1} X
  inline void VarInstationary_Base::
  SolveOperatorDhPlusdtSh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y)
  {
    VectReal_wp X2(X);
    this->Assemble(X2);
    Glob_mat_Dh->SolveOperatorDhPlusdtSh(alpha, X2, Y);
  }
  
  
  //! solves the matrix Dh, Y = Dh^{-1} Y
  inline void VarInstationary_Base::SolveOperatorDh(VectReal_wp& Y)
  {
    this->Assemble(Y);
    Glob_mat_Dh->SolveOperatorDh(Y);
  }
  
  //! solves the matrix Dh, Y = Dh^{-1} Y
  inline void VarInstationary_Base::SolveOperatorDh(VectComplex_wp& Y)
  {
    cout << "not implemented" << endl;
    abort();
    // this->Assemble(Y);
    //Glob_mat_Dh->SolveOperatorDh(Y);
  }
  
  //! solves the matrix L where Dh = L L^t, Y = L^{-1} Y or Y = L^{-T} Y
  inline void VarInstationary_Base::SolveCholeskyDh(const SeldonTranspose& transA, VectReal_wp& Y)
  {
    Glob_mat_Dh->SolveCholeskyDh(transA, Y);
  }
  
  
  //! solves the matrix Dh, Y = Y + alpha Dh^{-1} X
  inline void VarInstationary_Base::
  SolveOperatorDh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y)
  {
    VectReal_wp X2(X);
    this->Assemble(X2);
    Glob_mat_Dh->SolveOperatorDh(alpha, X2, Y);
  }
  
  
  //! applies the scalar stiffness matrix RhS
  /*!
    The method performs the operation :
    Yh = beta Yh + alpha RhS Uh
    where the evolution system can be written as :
    Dh dU/dt + Sh U + RhS V = 0
    Bh dV/dt + ShV V + RhV U = 0
  */
  inline void VarInstationary_Base
  ::ApplyOperatorRhScalar(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
			  const Real_wp& beta, VectReal_wp& Yh, bool extrapol)
  {
    if (Glob_mat_RhS.GetM() > 0)
      MltAdd(alpha, Glob_mat_RhS, Uh, beta, Yh);
    else
      this->MltAddStiffnessScalar(alpha, global_level_operator, Uh, beta, Yh);
  }
  
  
  //! applies the vectorial stiffness matrix RhV
  /*!
    The method performs the operation :
    Yh = beta Yh + alpha RhV Uh
    where the evolution system can be written as :
    Dh dU/dt + Sh U + RhS V = 0
    Bh dV/dt + ShV V + RhV U = 0
  */  
  inline void VarInstationary_Base
  ::ApplyOperatorRhVectorial(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
			     const Real_wp& beta, VectReal_wp& Yh, bool extrapol)
  {
    if (Glob_mat_RhV.GetM() > 0)
      MltAdd(alpha, Glob_mat_RhV, Uh, beta, Yh);
    else
      this->MltAddStiffnessVectorial(alpha, global_level_operator, Uh, beta, Yh);
  }
  
  
  //! applies the mass matrix Bh, Prod_Uh = beta Prod_Uh + alpha Bh Uh
  inline void VarInstationary_Base::
  ApplyOperatorBh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
                  const Real_wp& beta, VectReal_wp& Yh)
  {
    Glob_mat_Dh->ApplyOperatorBh(alpha, t, Uh, beta, Yh);
  }
  
  
  //! applies the matrix Bh - dt/2 ShV, Prod_Uh = beta Prod_Uh + alpha (Bh - dt/2 ShV) Uh
  inline void VarInstationary_Base::
  ApplyOperatorBhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
			   const Real_wp& beta, VectReal_wp& Yh)
  {
    Glob_mat_Dh->ApplyOperatorBhMinusdtSh(alpha, t, Uh, beta, Yh);
  }
  
  
  //! solves the matrix Bh + dt/2 ShV, Y = (Bh + dt/2 ShV)^{-1} Y
  inline void VarInstationary_Base::
  SolveOperatorBhPlusdtSh(VectReal_wp& Y)
  {
    Glob_mat_Dh->SolveOperatorBhPlusdtSh(Y);
  }
  
  
  //! solves the matrix Bh + dt/2 ShV, Y = Y + alpha (Bh + dt/2 ShV)^{-1} X
  inline void VarInstationary_Base::
  SolveOperatorBhPlusdtSh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y)
  {
    this->Glob_mat_Dh->SolveOperatorBhPlusdtSh(alpha, X, Y);
  }
  
  
  //! solves the matrix Bh, Y = Bh^{-1} Y
  inline void VarInstationary_Base::SolveOperatorBh(VectReal_wp& Y)
  {
    this->Glob_mat_Dh->SolveOperatorBh(Y);
  }
  
  
  //! solves the matrix Dh, Y = Y + alpha Bh^{-1} X
  inline void VarInstationary_Base::
  SolveOperatorBh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y)
  {
    this->Glob_mat_Dh->SolveOperatorBh(alpha, X, Y);
  }
  
  
  //! assembling the vector Vh between processors (values on shared dofs are summed)
  inline void VarInstationary_Base::Assemble(VectReal_wp& Vh) const
  {
    if (var_problem.FormulationDG() == ElementReference_Base::CONTINUOUS)
      this->var_comm.AddDomains(Vh);
  }
  
  
  //! evaluation of Y = G(tn,X)
  inline void VarInstationary_Base
  ::EvaluateFunction(const Real_wp& tn, const VectReal_wp& X, VectReal_wp& Y,
		     bool invert_mass, bool source)
  {
    this->EvaluateDerivativeFunction(tn, 0, X, Y, invert_mass, source);
  }
  
  
  //! evaluation of Y = G(tn,X)
  inline void VarInstationary_Base::
  EvaluateFunction(const Real_wp& tn, const Real_wp& alpha, 
                   int level, const VectReal_wp& X, VectReal_wp& Y, bool invert_mass, bool source)
  {
    this->EvaluateDerivativeFunction(tn, 0, alpha, level, X, Y, invert_mass, source);
  }
  
  
  //! method called during time iterations for outputs 
  /*!
    \param[in] nb_iter iterate number
    \param[in] tn time
    \param[in] Y solution vector
  */
  inline void VarInstationary_Base
  ::GiveIterate(int nb_iter, const Real_wp& tn, Vector<Real_wp>& Y)
  {
    this->WriteSnapshot(nb_iter, tn, Y, display_computational_time);
  }
  
  
  //! for a vectorial iterate provided by time scheme
  inline void VarInstationary_Base
  ::GiveVectorialIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y)
  {
  }
  
  
  //! last iterate provided by time scheme
  inline void VarInstationary_Base
  ::GiveFinalIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Un)
  {
    GiveIterate(nb_iter, tn, Un);
  }
  
  
  //! number of iterations to reach end of the experience, and time step dt
  inline void VarInstationary_Base
  ::GiveNumberIterations(const Real_wp& dt,int nb_max_iter)
  {
  }  

  
  template<class Dimension> template<class TypeEquation>
  inline VarInstationary_Dim<Dimension>::VarInstationary_Dim(HyperbolicProblem<TypeEquation>& var)
    : VarInstationary_Base(var), var_problem(var.var_harmonic), var_output(var.var_harmonic),
      var_source(var.var_harmonic), var_boundary(var.var_harmonic)
  {
  }
  
  
  template<class Dimension> template<class TypeEquation>
  inline VarInstationary_Fem<Dimension>::VarInstationary_Fem(HyperbolicProblem<TypeEquation>& var)
    : VarInstationary_Dim<Dimension>(var), var_problem(var.var_harmonic)
  {
  }
  
  
  template<class TypeEquation>
  inline VarInstationary<TypeEquation>::VarInstationary()
    : VarInstationary_Fem<Dimension>(GetLeafClass())
  {
  }
  
  
  template<class TypeEquation>
  inline HyperbolicProblem<TypeEquation>& VarInstationary<TypeEquation>::GetLeafClass()
  {
    return static_cast<HyperbolicProblem<TypeEquation>& >(*this);
  }

  template<class TypeEquation>
  inline VirtualMassMatrix* VarInstationary<TypeEquation>::GetNewMassMatrix()
  {
    if (var_harmonic.FormulationDG() != ElementReference_Base::CONTINUOUS)
      return new DiscontinuousUnsteadyMassMatrix<Dimension>(GetLeafClass());
    
    return new ContinuousUnsteadyMassMatrix<Dimension>(GetLeafClass());
  }
  

  template<class TypeEquation>
  inline VirtualProjectorFEM<Real_wp, typename TypeEquation::TypeEquationStationary::Dimension>*
  VarInstationary<TypeEquation>::GetNewInitialCondition()
  {
    return new InitialUserFunction<Real_wp, Dimension>(var_harmonic, this->param_initial_condition);
  }


  template<class TypeEquation>
  inline void VarInstationary<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    var_harmonic.SetInputData(description_field, parameters);
    VarInstationary_Fem<Dimension>::SetInputData(description_field, parameters);
  }

  template<class TypeEquation>
  inline Real_wp VarInstationary<TypeEquation>::EvaluateCFL(bool init_computation)
  {
    return Montjoie::EvaluateCFL(this->GetLeafClass(),
                                 this->var_time_scheme.GetStabilityFunction(), init_computation);
  }
  

  template<class TypeEquation>
  inline VarInstationary_Dim<typename TypeEquation::TypeEquationStationary::Dimension>*
  VarInstationary<TypeEquation>::GetNewHyperbolicProblem() const
  {
    return new HyperbolicProblem<TypeEquation>();
  }
  
}

#define MONTJOIE_FILE_VAR_INSTATIONARY_INLINE_CXX
#endif
