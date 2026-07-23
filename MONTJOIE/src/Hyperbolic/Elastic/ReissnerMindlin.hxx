#ifndef MONTJOIE_FILE_REISSNER_MINDLIN_HXX

namespace Montjoie
{
  
  //! class to specify resolution of Reissner-Mindlin equation for H1 finite element
  /*!
    Reissner-Mindlin equation has the following form :
    rho delta^2/12 d^2 theta/dt^2 - delta^2/12 div(C' eps(theta)) + mu (grad(u) + theta) = f
    rho d^2 u/dt^2 - div( mu (grad(u) + theta)) = f
    
    cf. Grob's thesis for details about the equations
   */
  template<class T>
  class ReissnerMindlinEquation_Base : public GenericEquation<T>
  {
  public :
    typedef Dimension2 Dimension;
    
    enum {nb_unknowns = 3, nb_unknowns_hdg=0,
	  nb_components_en = 1, nb_components_hn = 2, nb_unknowns_scal = 3, nb_unknowns_vec = 6};
   
    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return false; }
    
    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var,
                                     int i, int j,
				     const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                     Vector1& dU, Vector2& dV);

    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector2& V);

    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
                                     int i, int j,
				     const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                     Vector1& dU, Vector2& V);

    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
                                     int i, int j,
				     const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                     Vector1& U, Vector2& dV);

    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                  int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                  int ref, MatStiff& Cgrad_grad);

    template<class TypeEquation, class T0, class Matrix1>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, Matrix1& mass);

    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int i, int j,
                                 const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
  };
  
  
  //! class to specify resolution of Reissner-Mindlin equation for H1 finite element
  class ReissnerMindlinEquation : public ReissnerMindlinEquation_Base<Real_wp>
  {
  public :

    template<class TypeEquation, class T0, class Matrix1>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, Matrix1& mass);
    
    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector2& V);
    
  };


  //! class to specify resolution of Reissner-Mindlin equation for H1 finite element
  class HarmonicReissnerMindlinEquation
    : public ReissnerMindlinEquation_Base<Complex_wp>
  {
  public :
    template<class TypeEquation, class T0, class Matrix1>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, Matrix1& mass);
    
    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector2& V);
    
  };
    

  //! base class to solve stationary Reissner-Mindlin equations
  class ReissnerMindlin_Base
  {
  public :
    //! indices rho
    Vector<ScalarPhysicalIndice<Dimension2, Real_wp> > ref_delta;    
    //! indices rho
    Vector<ScalarPhysicalIndice<Dimension2, Real_wp> > ref_rho;
    //! angles d'orthotropie
    Vector<ScalarPhysicalIndice<Dimension2, Real_wp> > ref_theta_orthotrope;
    //! Young coefficient
    Vector<VectorPhysicalIndice<Dimension2, 2, Real_wp> > ref_Eyoung;
    //! Poisson coefficient
    Vector<VectorPhysicalIndice<Dimension2, 2, Real_wp> > ref_nu_poisson;
    //! shear modulus
    Vector<VectorPhysicalIndice<Dimension2, 3, Real_wp> > ref_shear_modulus;
    //! Timoshenko coefficient
    Vector<VectorPhysicalIndice<Dimension2, 2, Real_wp> > ref_timoshenko;
    
    Vector<TinySymmetricTensor<Real_wp, 2> > ref_tensorC;
    Vector<TinyMatrix<Real_wp, Symmetric, 2, 2> > ref_Grot;
    
    Real_wp omega_bar, u_bar;

  private:
    DistributedProblem<Dimension2>& var_problem;
    
  public:
    template<class TypeEquation>
    ReissnerMindlin_Base(EllipticProblem<TypeEquation>&);
    
    void InitIndices(int n);
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    int GetNbPhysicalIndices() const;
    bool IsVaryingMedia(int i) const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
			   IVect& num_ref, IVect& num_index, IVect& num_component,
			   Vector<bool>& compute_grad, Vector<bool>& compute_hess);
    
    template<class Complexe>
    void UpdateShiftAdimensionalization(Complexe&, Complexe&);
    
    template<class Complexe>
    void UpdateEigenvaluesAdimensionalization(Vector<Complexe>& lambda, Vector<Complexe>& lambda_imag,
					      Matrix<Complexe, General, ColMajor>& eigen_vec);
    
    void ComputePhysicalCoefficients();
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
  };

  
  template<class TypeEquation>
  class ReissnerMindlin_Eq : public ReissnerMindlin_Base, public VarHarmonic<TypeEquation>
  {
    typedef typename TypeEquation::Complexe Complexe;
    
  public:
    ReissnerMindlin_Eq();
    
    void InitIndices(int n);
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    int GetNbPhysicalIndices() const;
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_cplx,
			   Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
			   IVect& num_ref, IVect& num_index, IVect& num_component,
			   Vector<bool>& compute_grad, Vector<bool>& compute_hess);
    
    void UpdateShiftAdimensionalization(Real_wp&, Real_wp&);
    
    void UpdateEigenvaluesAdimensionalization(Vector<Real_wp>& lambda, Vector<Real_wp>& lambda_imag,
					      Matrix<Real_wp, General, ColMajor>& eigen_vec);

    void UpdateShiftAdimensionalization(Complex_wp&, Complex_wp&);
    
    void UpdateEigenvaluesAdimensionalization(Vector<Complex_wp>& lambda, Vector<Complex_wp>& lambda_imag,
					      Matrix<Complex_wp, General, ColMajor>& eigen_vec);
    
    void ComputePhysicalCoefficients();
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
  };


  //! class used to solve real Reissner-Mindlin equation with H1 elements
  template<>
  class EllipticProblem<ReissnerMindlinEquation>
    : public ReissnerMindlin_Eq<ReissnerMindlinEquation>
  {    
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };
  
  
  template<>
  class FemMatrixFreeClass<Real_wp, ReissnerMindlinEquation>
    : public FemMatrixFreeClass_Eq<Real_wp, ReissnerMindlinEquation>
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<ReissnerMindlinEquation>& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };

  
  //! class used to solve harmonic Reissner-Mindlin equation with H1 elements
  template<>
  class EllipticProblem<HarmonicReissnerMindlinEquation>
    : public ReissnerMindlin_Eq<HarmonicReissnerMindlinEquation>
  {    
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };
  
  
  //! base class for Reissner-Mindlin equation in time-domain
  class TimeReissnerMindlinEquation
  {
  public :
    typedef ReissnerMindlinEquation TypeEquationStationary; //!< stationary equation
    
  };
  
  
  //! base class to solve Reissner-Mindlin equation in time-domain
  template<class TypeEquation> 
  class TimeReissnerMindlin_Eq : public VarInstationary<TypeEquation>
  {
  };
  
  
  //! class used to solve real Reissner-Mindlin equation with H1 elements
  template<>
  class HyperbolicProblem<TimeReissnerMindlinEquation>
    : public TimeReissnerMindlin_Eq<TimeReissnerMindlinEquation>
  {    
  };
  
}

#define MONTJOIE_FILE_REISSNER_MINDLIN_HXX
#endif
