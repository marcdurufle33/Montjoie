#ifndef MONTJOIE_FILE_CAUCHY_POISSON_CXX

namespace Montjoie
{
  //! definition of Cauchy-Poisson equation
  template<class Dimension>
  class CauchyPoissonEquation : public AcousticEquation<Dimension>
  {
  public :
    typedef LaplaceEquation<Dimension> TypeEquation_Harmonic; //!< related stationary problem
    
  };
  
  
  //! class to solve Cauchy-Poisson equation
  template<class Dimension>
  class HyperbolicProblem<CauchyPoissonEquation<Dimension> >
    : public TimeAcoustic_Eq<CauchyPoissonEquation<Dimension> >
  {
  public :
    Real_wp teta; //!< teta parameter for teta scheme (implicit scheme)
    IVect TypeConditionImpedance; //!< absorbing or free-surface condition or other ?
    IVect ddl_Sh; //!< dofs with free-surface condition
    IVect ddl_Sigh; //!< dofs with absorbing condition
    VectReal_wp value_Sh; 
    //!< diagonal matrix \f$ \int \varphi_i \varphi_j \f$ for free-surface condition
    VectReal_wp value_Sigh; 
    //!< diagonal matrix \f$ \int \varphi_i \varphi_j \f$ for absorbing condition
    //! surface condition
    static const int CONDITION_WAVE = 1, CONDITION_ABSORBING = 2;
    
    //! default constructor
    HyperbolicProblem()
      : TimeAcoustic_Eq<CauchyPoissonEquation<Dimension> >()
    { 
      teta = 0.25;
      TypeConditionImpedance.Reallocate(50);
      TypeConditionImpedance.Fill(0);
    }
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    bool IsImplicitScheme() const;
    
    void TreatBoundaryConditions();
    
    // pour lancer les iterations
    void RunTimeIterations();
    
  };
  

  //! modification of parameters of class CauchyPoisson according data file
  template<class Dimension>
  void HyperbolicProblem<CauchyPoissonEquation<Dimension> >
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    TimeAcoustic_Eq<CauchyPoissonEquation<Dimension> >
      ::SetInputData(description_field, parameters);
    
    if (!description_field.compare("SurfaceFluid"))
      {
	// free-surface boundary
	for (int i = 0; i < parameters.GetM(); i++)
	  TypeConditionImpedance(to_num<int>(parameters(i))) = CONDITION_WAVE;
      }
    else if (!description_field.compare("SurfaceAbsorbing"))
      {
	// absorbing boundary
	for (int i = 0; i < parameters.GetM(); i++)
	  TypeConditionImpedance(to_num<int>(parameters(i))) = CONDITION_ABSORBING;
      }
  }
  

  //! always implicit
  template<class Dimension>
  bool HyperbolicProblem<CauchyPoissonEquation<Dimension> >::IsImplicitScheme() const
  {
    return true;
  }
  

  //! boundary conditions, computation of matrix 
  //! coming from absorbing condition and free-surface condition
  template<class Dimension>
  void HyperbolicProblem<CauchyPoissonEquation<Dimension> >::TreatBoundaryConditions()
  {
    // this->var_harmonic.mat_glob_iterative.nature_matrix.SetCoefStiffness(teta);
    for (int i = 0; i < this->var_harmonic.ref_mu.GetM(); i++)
      {
	this->var_harmonic.ref_mu(i).Mlt(teta);
	this->var_harmonic.ref_invMu(i).Mlt(1.0/teta);
      }
    
    ImpedanceFunction_Base<Real_wp, Dimension, 1> fimped(this->var_harmonic);
    
    int N = this->var_harmonic.mesh_num.GetNbDof();
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> mat_abc(N, N);
    Real_wp invDt2 = 1.0/(this->deltat*this->deltat);
    for (int ref = 0; ref < this->var_harmonic.Param_condition.GetM(); ref++)
      {

	if (TypeConditionImpedance(ref) == CONDITION_WAVE)
	  Mlt(invDt2, this->var_harmonic.Param_condition(ref));
	else if (TypeConditionImpedance(ref) == CONDITION_ABSORBING)
	  {
	    Mlt(0.5/this->deltat, this->var_harmonic.Param_condition(ref));
	    this->var_harmonic.mesh.SetBoundaryCondition(ref, BoundaryConditionEnum::LINE_INSIDE);
	  }
      }
    
    GlobalGenericMatrix<Real_wp> nat_mat;
    this->var_harmonic.AddMatrixImpedanceBoundary(1.0, this->var_harmonic.mesh.GetBoundaryCondition(),
						  BoundaryConditionEnum::LINE_IMPEDANCE,
						  nat_mat, mat_abc, 0, 0, fimped);
    
    ddl_Sh.Clear(); value_Sh.Clear();
    for (int i = 0; i < N; i++)
      {
	Real_wp val = mat_abc(i,i);
	if (abs(val) > this->var_harmonic.epsilon_machine)
	  {
	    ddl_Sh.PushBack(i);
	    value_Sh.PushBack(val);
	  }
      }
    
    for (int ref = 0; ref < this->var_harmonic.Param_condition.GetM(); ref++)
      {
	if (TypeConditionImpedance(ref) == CONDITION_WAVE)
	  this->var_harmonic.mesh.Cond_curve(ref) = this->var_harmonic.LINE_INSIDE;
	else if (TypeConditionImpedance(ref) == CONDITION_ABSORBING)
	  this->var_harmonic.mesh.Cond_curve(ref) = this->var_harmonic.LINE_IMPEDANCE;
      }
    
    mat_abc.Clear(); mat_abc.Reallocate(N,N);
    this->var_harmonic.AddMatrixImpedanceBoundary(1.0, this->var_harmonic.mesh.GetBoundaryCondition(), 
						  BoundaryConditionEnum::LINE_IMPEDANCE,
						  nat_mat, mat_abc, 0, 0, fimped);
    
    ddl_Sigh.Clear(); value_Sigh.Clear();
    for (int i = 0; i < N; i++)
      {
	Real_wp val = mat_abc(i,i);
	if (abs(val) > epsilon_machine)
	  {
	    ddl_Sigh.PushBack(i);
	    value_Sigh.PushBack(val);
	  }
      }
    
    for (int ref = 0; ref < this->var_harmonic.Param_condition.GetM(); ref++)
      {
	if (TypeConditionImpedance(ref) == CONDITION_WAVE)
	  this->var_harmonic.mesh.Cond_curve(ref) = BoundaryConditionEnum::LINE_IMPEDANCE;
	else if (TypeConditionImpedance(ref) == CONDITION_ABSORBING)
	  this->var_harmonic.mesh.Cond_curve(ref) = BoundaryConditionEnum::LINE_IMPEDANCE;
      }
    // this->var_harmonic.ComputeMassMatrix();
  }
  
  //! time iterations
  template<class Dimension>
  void HyperbolicProblem<CauchyPoissonEquation<Dimension> >::RunTimeIterations()
  {
    int N = this->var_harmonic.mesh_num.GetNbDof(); Real_wp zero(0);
    // allocation des vecteurs d'iterations, et mise a zero
    VectReal_wp Un(N), Unm1(N), Xn(N), Ah_Un(N);
    Un.Fill(zero); Unm1.Fill(zero); Xn.Fill(zero); Ah_Un.Fill(zero);
    
    // initialisation des iterations en temps
    Real_wp t = this->initial_time;
    Real_wp dt = this->deltat; // invDt2 = 1.0/(dt*dt);
    int nb_iter = 0;
    
    // on lance la boucle en temps
    // (Sigma_h/ (2 dt) + Sh / dt^2 + teta K_h) U^{n+1} + ( (1-2 teta) K_h - 2 Sh / dt^2) U^n
    // + (teta K_h + S_h/ dt^2 - Sigma_h / (2 dt) ) U^{n-1} = F_h^n
    // on note la matrice A_h = (Sigma_h / (2 dt) + Sh / dt^2 + teta K_h) 
    // (c'est mat_glob_iterative)
    // on a alors
    // A_h U^{n+1} + [(1-2teta)/teta (A_h - Sigma_h / (2 dt) ) + 1/ teta Sh / dt^2 ] U^n 
    // + (A_h - Sigma_h / dt) U^{n-1}
    DISP(t); DISP(this->final_time);    
    while (t < this->final_time)
      {
	// on stocke les instantanes si demande
	this->WriteSnapshot(nb_iter, t, Un);
	
	// calcul de la source
	Ah_Un.Zero();
	if (t < this->tlimit_source)
	  this->AddPrimitiveSourceAtTime(1.0, t, 1, Ah_Un);
	
	// on rajoute (2 teta-1)/teta A_h U^n - A_h U^{n-1}
	Copy(Unm1, Xn); Add( (1-2.0*teta)/teta, Un, Xn);
	this->Glob_mat_Kh->MltAddVector(-1.0, Xn, 1.0, Ah_Un);
	
	// on rajoute 1/ teta S_h / dt^2
	Real_wp alpha = 1.0/teta;
	for (int i = 0; i < ddl_Sh.GetM(); i++)
	  Ah_Un(ddl_Sh(i)) += alpha*value_Sh(i)*Un(ddl_Sh(i));
	
	// on rajoute -(2 teta - 1)/ teta Sigma_h / (2 dt) U^n  + Sigma_h / dt U^{n-1}
	alpha = (1-2.0*teta)/teta;
	for (int i = 0; i < ddl_Sigh.GetM(); i++)
	  {
	    int num_dof = ddl_Sigh(i);
	    Ah_Un(num_dof) += value_Sigh(i)*(alpha*Un(num_dof) + 2.0*Unm1(num_dof));
	  }
	
	// DISP(this->var_harmonic.print_level);
	// on resout A_h X = B
	Xn.Zero(); 
	//ComputeSolution(Ah_Un, Xn);
	abort();
	
	// on met a jour les vecteurs
	Copy(Un, Unm1); Copy(Xn, Un);
        
        nb_iter++;
	t = nb_iter*dt + this->initial_time;
      }
  }

} // end namespace

#define MONTJOIE_FILE_CAUCHY_POISSON_CXX
#endif
