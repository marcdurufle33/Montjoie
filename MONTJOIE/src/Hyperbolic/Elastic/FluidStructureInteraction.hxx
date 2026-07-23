#ifndef MONTJOIE_FILE_FLUID_STRUCTURE_INTERACTION_HXX

namespace Montjoie
{
  
  //! base class for coupling two equations in adjacent domains
  template<class Dimension>
  class FluidStructureInteraction_Base : public InputDataProblem_Base
  {
  private:
    VarProblem<Dimension>& var_fluid;
    VarProblem<Dimension>& var_solid;
    
  protected:
    //! references for the interface
    IVect ref_interface;
    //! references for the fluid
    IVect ref_fluid;
    //! 3-D faces for fluid and solid
    IVect NumElementFluid, NumElementSolid, NumLocFluid, NumLocSolid;
    //! for parallel computations
    IVect ProcElementFluid, ProcElementSolid;
    Vector<IVect> ConnecSolid, ConnecFluid;
    //! number of degrees of freedom
    int nodl;

  public:
    FluidStructureInteraction_Base(VarProblem<Dimension>& var_f,
				   VarProblem<Dimension>& var_s);

    int GetNbDof() const;
    
    void ComputeMeshAndFiniteElement(const string& name_element, bool first = true);
    
    const ElementReference<Dimension, 1>&
    GetFluidData(int i, int& proc_fluid,
		 int& rot2, IVect& NodleAcous, int& num_loc_acous) const;
    
    int GetFluidRotation(int i, int& proc_fluid) const;

    void InitIndices(int nb);
    void ConstructAll(const string& input_file, const string& name_element, bool first = true);
    
  };
  

  //! Base class for fluid-structure coupling
  template<class Complexe, class Dimension>
  class FluidStructureInteraction_Elas : public FluidStructureInteraction_Base<Dimension>
  {
    typedef Vector<Complexe> VectComplexe;
    typedef typename Dimension::R_N R_N;
    
  private:
    VarHarmonic_Base<Complexe, Dimension>& var_fluid;
    VarHarmonic_Base<Complexe, Dimension>& var_solid;

  protected:
    IVect NbPointsQuadFluid, NbPointsQuadSolid;

  public:
    FluidStructureInteraction_Elas(VarHarmonic_Base<Complexe, Dimension>& var_f,
				   VarHarmonic_Base<Complexe, Dimension>& var_s);
    
    void ConstructAll(const string& input_file, const string& name_element, bool first = true);
    
    void SendUnFluid(Vector<VectComplexe>& UnQuad);
    void SendUnSolid(Vector<VectComplexe>& UnQuad);
    
    void ApplyCouplingOperatorB(const Complexe& alpha, const Vector<Complexe>& X,
                                const Complexe& beta, Vector<Complexe>& Y);
    
    void ApplyCouplingTransposeB(const Complexe& alpha, const Vector<Complexe>& X,
                                 const Complexe& beta, Vector<Complexe>& Y);
    
    void AddOperatorB(const Complexe& alpha, VirtualMatrix<Complexe>& A);
    
  };
  
  
  //! object for solving elasto-acoustic system in time-harmonic domain
  template<class TypeEqFluid, class TypeEqSolid>
  class FluidStructureInteractionHarmonic :
    public FluidStructureInteraction_Elas<typename TypeEqFluid::Complexe,
					 typename TypeEqFluid::Dimension>
  {
  public :
    typedef typename TypeEqFluid::Complexe Complexe;
    typedef typename TypeEqFluid::Dimension Dimension;
    typedef typename Dimension::R_N R_N;
    
    //! object for the solution of acoustic equation in the fluid
    EllipticProblem<TypeEqFluid> var_fluid;
    
    //! object for the solution of elastic equation in the solid
    EllipticProblem<TypeEqSolid> var_solid;

  public:
    FluidStructureInteractionHarmonic();
    
    void SetInputData(const string& keyword, const Vector<string>& param);
    void RunAll(const string& input_file, const string& name_element);
    
  };
  

  //! object for the solution of fluid-structure interaction in time-domain
  /*!
    we are considering coupling between acoustic equation and elastodynamic 
    rho du/dt - div (sigma) = f_u
    dsigma/dt - C eps(u) = 0
    1/(c^2 rho_f) dP/dt + div v = f_p
    rho_f dv/dt + grad P = 0
    where rho : density of the solid, C : elasticity tensor of the solid
    rho_f : density of the fluid, c : speed of sound in the fluid
    Unknowns are the pressure P, displacement of the fluid v, displacement 
    of the solid u, stress tensor sigma
    with the following transmission conditions
    u.n = v.n (normal component of displacement is continuous)
    sigma n . n = -P 
    
    The obtained variational formulation is therefore :
    M_u dU/dt + R_u Sigma + B* P = f_u
    M_Sigma dSigma/dt - R_sigma U = 0
    M_p dP/dt - R_p V - B U = f_p
    M_v dV/dt + R_v P = 0
    with M_u, M_Sigma, M_p, M_v standard mass matrices (tensor C is contained in M_Sigma)
    stifness matrices are equal to
    (R_u) Sigma = \int_Omega  Sigma : grad u^T
    (R_Sigma) U = \int_Omega  grad u : Sigma^T
    (R_p) V = \int_Omega  v \cdot grad P^T
    (R_v) P = \int_Omega  grad P \cdot v^T
    coupling matrix B is equal to :
    B U = \int_\Gamma U.n P^T
    where Gamma is the interface between the fluid and solid
    the normal n is outward from the solid to the fluid
    
    cf. Diaz's thesis for details about the equations
   */
  template<class TypeEqFluid, class TypeEqSolid>
  class FluidStructureInteraction :
    public FluidStructureInteraction_Elas<typename TypeEqFluid::TypeEquationStationary::Complexe,
					 typename TypeEqFluid::TypeEquationStationary::Dimension>,
    public VirtualOdeSystem<Real_wp>
  {
  public :
    typedef typename TypeEqFluid::TypeEquationStationary::Dimension Dimension;
    typedef typename Dimension::R_N R_N;
    
    //! object for the resolution of acoustic equation in the fluid
    HyperbolicProblem<TypeEqFluid> var_fluid;
    
    //! object for the resolution of elastic equation in the solid
    HyperbolicProblem<TypeEqSolid> var_solid;
    
  public:
    inline FluidStructureInteraction()
      : FluidStructureInteraction_Elas<Real_wp, Dimension>(var_fluid.var_harmonic, var_solid.var_harmonic) {}
    
    void SetInputData(const string& keyword, const Vector<string>& param);    
    
    Real_wp GetFrequency() const;
    Real_wp GetTimeStep() const;
    Real_wp GetInitialTime() const;
    Real_wp GetFinalTime() const;
    
    void GiveIterate(int n, const Real_wp& t, VectReal_wp& Y);
    
    void ConstructAll(const string& input_file, const string& name_element);
    
    void SolveMassPressure(VectReal_wp& X);
    void SolveMassDisplacementFluid(VectReal_wp& X);
    void SolveMassDisplacementSolid(VectReal_wp& X);
    void SolveMassStress(VectReal_wp& X);
    
    void ApplyStiffnessPressure(const Real_wp& alpha, const VectReal_wp& X,
                                const Real_wp& beta, VectReal_wp& Y);

    void ApplyStiffnessDisplacementFluid(const Real_wp& alpha, const VectReal_wp& X,
                                         const Real_wp& beta, VectReal_wp& Y);

    void ApplyStiffnessDisplacementSolid(const Real_wp& alpha, const VectReal_wp& X,
                                         const Real_wp& beta, VectReal_wp& Y);

    void ApplyStiffnessStress(const Real_wp& alpha, const VectReal_wp& X,
                              const Real_wp& beta, VectReal_wp& Y);
    
    inline void SetDirichletCondition(const Real_wp& tn, int nb_deriv,
				      VectReal_wp& X, Real_wp alpha = 0) {}
    
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
				    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert = true, bool source = true);

    void EvaluateFunction(const Real_wp& tn, const VectReal_wp& X, VectReal_wp& Y,
                          bool invert_mass = true, bool source = true);
    
  };
  
  
  //! second-order leap frog scheme for vibro-acoustic equation
  class FluidStructureLeapFrogScheme
  {
  public :
    // time step
    Real_wp dt;
    
    // vectors
    VectReal_wp Pn, Sn, Un_half, Vn_half;
    VectReal_wp ProdPn, ProdSn, ProdUn, ProdVn;
    
    // initialization of all the variables
    template<class GenericPb>
    void SetInitialCondition(const Real_wp& t0, const Real_wp& dt_, GenericPb& var);
    
    // scheme is advanced of a time step
    template<class GenericPb>
    void Advance(const Real_wp& t, int n, GenericPb& var);
      
  };


  template<class Dimension>
  class AcousticAcousticInteraction :
    public FluidStructureInteraction_Base<Dimension>
  {
  public:
    Vector<int> DdlLambda;
    
    //! object for the resolution of acoustic equation in the fluid
    HyperbolicProblem<AcousticEquation<Dimension> > var_fluid;
    
    //! object for the resolution of elastic equation in the solid
    HyperbolicProblem<AcousticEquation<Dimension> > var_solid;

  public:
    inline AcousticAcousticInteraction()
      : FluidStructureInteraction_Base<Dimension>(var_fluid.var_harmonic, var_solid.var_harmonic) {}
    
    void SetInputData(const string& keyword, const Vector<string>& param);
    void ConstructAll(const string& input_file, const string& name_element);

    void AddOperatorB(const Real_wp& alpha, VirtualMatrix<Real_wp>& mat_sp);
    
    void ComputeMatrices(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& M,
			 DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& K);
    
  };
  
}

#define MONTJOIE_FILE_FLUID_STRUCTURE_INTERACTION_HXX
#endif
