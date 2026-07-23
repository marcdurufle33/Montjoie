#ifndef MONTJOIE_FILE_VIBRO_ACOUSTIC_SESSION_HXX

namespace Montjoie
{
  
  //! object for the resolution of vibro-acoustic equations
  /*!
    we are considering coupling between acoustic equation and Reissner-Mindlin equations :
    1/(rho_f c^2) d^2 P/dt^2 - div( 1/rho_f grad P) = f
    -1/rho_f dP/dn = du/dt on the plate
    rho_f : density of the fluid
    c^2 : sound speed in the fluid
    u : displacement in z-coordinate)
    P : primitive of the pressure
    FV =>  d^2/dt^2 \int_\Omega  1/(rho_f c^2) P phi  + \int_\Omega 1/rho_f grad P \cdot grad phi
                             + d/dt \int_\Gamma u [ phi] = \int_\Omega f phi
    Gamma : surface where the plate is located (surface z = cte)
    [phi] = phi^-  -  phi^+
    Evolution equation => Mh d^2 P /dt^2  + Kh P + Ch^* d Lambda/dt = F
    where Lambda = [teta_x, teta_y, u]
    
    Lambda is satisfying Reissner-Mindlin equations :
    rho delta^3/12 d^2 theta/dt^2 - delta^3/12 div(C eps(theta) ) + delta mu (grad(u) + theta) = f
    rho delta d^2 u/dt^2 - delta div(mu (grad(u) + theta) ) = f + [ dP/dt]
    
    Evolution equation : M^L d^2 Lambda/dt^2 + K^L Lambda - Ch dP/dt = 0        
    cf. Grob's thesis for details about the equations
   */
  template<class TypeEqPlate, class TypeEqVol>
  class VibroAcousticSession : public VirtualOdeSystem<Real_wp>, public InputDataProblem_Base,
			       public VirtualMatrix<Real_wp>
  {
  public :    
    //! celerity of acoustic waves in the fluid
    Real_wp rho, c2;    
    //! references for the plate
    IVect ref_plate;
    //! 3-D faces corresponding to 2-D mesh elements
    /*! 
      NumFaceMinus contains the number of 3-D faces  which are on the side z < 0
      (these faces are on the plate)
      while NumFacePlus is related to the side z > 0
      Those arrays are not used anymore
    */
    IVect NumFaceMinus, NumFacePlus;
    //! 3-D dofs corresponding to 2-D dofs
    /*!     
      ddl_vol_minus contains the number of volume dofs on the plate, but
      on the side z < 0, while ddl_vol_plus contains volumes dofs for z > 0
    */
    IVect ddl_vol_minus, ddl_vol_plus;
    //! dof numbers on the plate for acoustic and mecanic
    /*!
      surface dof numbers when using the order of approximation for acoustic equation
      or for Reissner-Mindlin (i.e. mecanic) equation
     */
    Vector<IVect> NodleAcous, NodleMeca;
    //! omega_k J_i(\xi_k) on the plate (with acoustic dofs)
    Vector<VectReal_wp> JacobianAcous;
    //! dofs on the plate (acoustic) on Dirichlet boundaries
    IVect DirichletAcous;
    //! projection from mecanic unknowns to acoustic unknowns
    FiniteElementInterpolator proj_meca;
    
    //! object for the resolution of acoustic equation in the volume
    HyperbolicProblem<TypeEqVol> var_volume;
    //! object for the resolution of Reissner-Mindlin equation on the plate
    HyperbolicProblem<TypeEqPlate> var_plate;
    
#ifdef SELDON_WITH_MPI
    IVect NumberProc_DofVol;
    Vector<IVect> DofVolMinusOnProc, DofVolPlusOnProc;
    IVect NumberProc_DofSurf;
    Vector<IVect> DofSurfMinusOnProc, DofSurfPlusOnProc;
#endif
    
    //! mass matrix for mecanic unknowns and volume unknowns (but restricted to the surface)
    VectReal_wp MassLambda;
    int nodl_acous;
    IVect NumberProcMatchingAcous;
    Vector<IVect> MatchingDofAcous;
    VectReal_wp invDiagonalSchur;
    IVect OverlapDofsUz;
    
    void SetInputData(const string& keyword, const Vector<string>& param);
    
    void ComputeMeshAndFiniteElement();
    void ConstructAll(const string&);
    void DistributeTransmissionDofs();
    
    void ExchangeSurfaceToVolume(VectReal_wp& U, VectReal_wp& Uminus, VectReal_wp& Uplus);
    void ExchangeVolumeToSurface(VectReal_wp& Uminus, VectReal_wp& Uplus,
                                 VectReal_wp& U, VectReal_wp&);
    void AssembleAcoustic(VectReal_wp& X) const;

    void ComputeDiagonalSchur(VectReal_wp& invMh, VectReal_wp&, VectReal_wp&,
                              bool assemble = true);
    
    template<class MatrixSparse>
    void ComputeCbarre(MatrixSparse& Ch);
    
    Real_wp GetCelerity() const;
    Real_wp GetFrequency() const;
    int GetOrder() const;
    Real_wp GetMeshSize() const;
    Real_wp GetTimeStep() const;
    Real_wp GetInitialTime() const;
    Real_wp GetFinalTime() const;
    
    int GetNbMecanicDofSurface() const;
    int GetNbAcousticDofSurface() const;
    int GetM() const;
    int GetNbDofVolume() const;
    
    VectReal_wp& GetAcousticMass();
    VectReal_wp& GetAcousticDamping();
    const VectReal_wp& GetMecanicMass() const;    
    
    void InitTimeIterations();

    int GetNbDof() const;
    
    void GiveNumberIterations(const Real_wp&, int);
    void GiveIterate(int n, const Real_wp& t, VectReal_wp& Y);
    void GiveFinalIterate(int n, const Real_wp& t, VectReal_wp& Y);

    inline void SetDirichletCondition(const Real_wp& tn, int nb_deriv,
				      VectReal_wp& X, Real_wp alpha = 0) {}
    
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
				    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert = true, bool source = true);

    void EvaluateFunction(const Real_wp& tn, const VectReal_wp& X, VectReal_wp& Y,
                          bool invert_mass = true, bool source = true);
    
    void MltAddCh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y) const;
    void MltAddChTranspose(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y) const;
    
    void SolveSchur(VectReal_wp & X, VectReal_wp & B);

    void MltAddVector(const Real_wp& alpha,
		      const VectReal_wp& X, const Real_wp& beta, VectReal_wp& Yvec) const;

    void MltAddVector(const Complex_wp& alpha,
		      const VectComplex_wp& X, const Complex_wp& beta, VectComplex_wp& Yvec) const;

    void MltAddVector(const Real_wp& alpha, const SeldonTranspose&,
		      const VectReal_wp& X, const Real_wp& beta, VectReal_wp& Yvec) const;

    void MltAddVector(const Complex_wp& alpha, const SeldonTranspose&,
		      const VectComplex_wp& X, const Complex_wp& beta, VectComplex_wp& Yvec) const;
    
    void MltVector(const VectReal_wp& X, VectReal_wp& Yvec) const;
    void MltVector(const VectComplex_wp& X, VectComplex_wp& Yvec) const;

    void MltVector(const SeldonTranspose&, const VectReal_wp& X, VectReal_wp& Yvec) const;
    void MltVector(const SeldonTranspose&, const VectComplex_wp& X, VectComplex_wp& Yvec) const;
    
  };
  
  
  //! second-order leap frog scheme for vibro-acoustic equation
  /*!
    This scheme is written as :
    Mh (P^(n+1) - 2 P^n + P^(n-1)) / dt^2  + Sh (P^(n+1) - P^(n-1))/(2 dt)
    + Kh P^n + Ch (Lambda^(n+1) - Lambda^(n-1)) / (2 dt) = F^n
    
    M^L (Lambda^(n+1) - 2 Lambda^n + Lambda^(n-1)) / dt^2 + K^L Lambda^n
    - Ch (P^(n+1) - P^(n-1)) / (2 dt) = 0
    This semi-explicit scheme is solved with a Schur-complement method
    to reduce the linear system to a system with only Lambda^n+1
   */
  class VibroAcousticLeapFrogScheme
  {
  public :
    // time step
    Real_wp dt;
    
    // vectors
    VectReal_wp Ph_n, Ph_nm1, KhPn;
    VectReal_wp Lambda_np1, Lambda_n, Lambda_nm1;
    VectReal_wp Gvol_minus, Gvol_plus, Gsurf_minus, Gsurf_plus;
    VectReal_wp KlambdaLn, Fuz, Fvol, U_np1;
    VectReal_wp ChLambdaMinus_n, ChLambdaPlus_n, ChLambdaMinus_nm1, ChLambdaPlus_nm1;
    VectReal_wp Vh_n, Vh_np1;
    bool first_order;
    
    VibroAcousticLeapFrogScheme();
    
    // initialization of all the variables
    template<class GenericPb>
    void SetInitialCondition(const Real_wp& t0, const Real_wp& dt_, GenericPb& var);
    
    // scheme is advanced of a time step
    template<class GenericPb>
    void Advance(const Real_wp& t, int n, GenericPb& var);
      
  };
  
}

#define MONTJOIE_FILE_VIBRO_ACOUSTIC_SESSION_HXX
#endif
