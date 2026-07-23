#ifndef MONTJOIE_FILE_STATIC_MAXWELL_3D_HXX

namespace Montjoie
{
  
  //! class to specify the resolution of stationary Maxwell equations
  /*!
    The Maxwell equation solved is
    \f$ \varepsilon \, E + \mbox{rot}( \frac{1}{\mu} \mbox{rot} E) \, = \, f \f$
   */    
  class StaticMaxwellEquation_3D : public MaxwellEquation3D_Base<Real_wp>
  {
  };


  //! class to specify the resolution of stationary Maxwell equations with HDG formulation
  class StaticMaxwellEquationHdg_3D : public MaxwellEquationHdg_3D<Real_wp>
  {
  };


  //! class to specify the resolution of static maxwell equation with dg method
  /*!
    The maxwell equation solved is
    \f$ \varepsilon \, E + \mbox{rot}( H) \, = \, f \f$
    \f$ \mu H - \mbox{rot} E \, = \, g \f$
   */
  class StaticMaxwellEquation_3D_DG
    : public GenericEquation<Real_wp>
  {
  public :
    typedef Dimension3 Dimension;
    
    static const bool FirstOrderFormulation = true;
    static bool store_dfjm1;
    
    enum { nb_unknowns = 6, nb_unknowns_hdg=3,
	   nb_components_en = 3, nb_components_hn = 3, nb_unknowns_vec = 3, nb_unknowns_scal = 3};
    
    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();

    template<class TypeEquation>    
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
				  int i, const ElementReference_Dim<Dimension3>&);
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& var,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			      int ref, MatMass& mass);

    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref,
				Vector1& U, Vector1& V);
    
    template<class TypeEquation, class T0, class Vector1>
    static void GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
				    const GlobalGenericMatrix<T0>& nat_mat,
				    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
				 int num_elem, int jloc,
				 const GlobalGenericMatrix<T0>& nat_mat, int ref,
				 MatStiff& Ephi_grad, MatStiff& Dgrad_phi);
    
    template<class TypeEquation,
	     class T0, class Vector1, class Vector2>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
				     int ref, Vector1& Vn, Vector2& Un);
    
    template<class TypeEquation,
	     class T0, class Vector1, class Vector2>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
				     int ref, Vector1& Un, Vector2& Vn);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, R3& normale, int ref, int iquad, int k,
			const GlobalGenericMatrix<T0>& nat_mat, int ref2, const GenericPb& vars,
			const ElementReference<Dimension3, 1>& Fb);

    template<class Matrix1, class GenericPb, class T0>
    static void GetPenalDG(Matrix1& Nabc, R3& normale, int iquad, int k, int nf,
			   const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
			   const GenericPb& vars, const ElementReference<Dimension3, 1>& Fb);
      
    template<class Vector1, class Vector2, class GenericPb, class T0>
    static void MltPenalDG(const R3& normale, const Vector1& Vn, Vector2& Un,
			   int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
			   int ref, int ref2, const GenericPb& vars,
			   const ElementReference<Dimension3, 1>& Fb);
    
    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(R3& normale, int ref, const Vector1& Vn, Vector1& Un,
			int num_elem1, int k, const GlobalGenericMatrix<T0>& nat_mat,
			int ref2, const EllipticProblem<TypeEquation>& vars,
			const ElementReference<Dimension3, 1>& Fb);
    
  };

  
  //! base class to solve static Maxwell equations
  template<class TypeEquation>
  class StaticMaxwell_3D : public HarmonicMaxwell_3D_Eq<TypeEquation>
  {
  public :
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension3, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension3, Real_wp>* >& rho_real, IVect& num_ref,
                           IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
    
  };


  //! base class to compute incident field in time-domain and Maxwell's equations
  class IncidentWaveProjector_StaticMaxwell3D : public IncidentWaveProjector<Real_wp, Dimension3>
  {
  protected :
    R3 k_wave; //!< wave vector
    R3 polar, polarH; //!< polarization vector
    Real_wp c0; //!< wave velocity
    Real_wp offset_arg; //!< f will be evaluated for k.x - c t - offset_arg
    
    VirtualTimeSource<Real_wp>* wave_pulse;
    
  private:
    const VarSourceProblem_Base& var_source;
    const VarBoundaryCondition_Base& var_boundary;
    const VarProblem<Dimension3>& var_problem;

  public :
    template<class TypeEquation>
    IncidentWaveProjector_StaticMaxwell3D(const EllipticProblem<TypeEquation>& vars_helm,
                                          IncidentWaveField<Real_wp, Dimension3>&);
    
    void Init(const Real_wp& t0);
    
    void EvaluateFunction(int i, int j, const R3& x, VectReal_wp& f);    
    void EvaluateFunctionGradient(int i, int j, const R3& x, VectReal_wp& f, VectReal_wp& df);
    
  };
  
  
  //! source for the scattering of an incident wave
  /*!
    the diffracted field is computed
   */
  class DiffractedWaveSource_StaticMaxwell3D : public VirtualSourceFEM<Real_wp, Dimension3>
  {
  protected :
    R3 k_wave; //!< wave vector
    R3 polar, polarH; //!< polarization vector
    Real_wp c0; //!< wave velocity
    bool element_insidePML; //!< element of a PML layer ?
    Real_wp t; //!< current time
    Real_wp offset_arg; //!< f will be evaluated for k.x - c t - offset_arg
    int num_deriv;
    bool scalar_eq;
    bool dirichlet_cond;
    
  private:
    IncidentWaveField<Real_wp, Dimension3>& incident_wave;
    VirtualTimeSource<Real_wp>* wave_pulse;
    
    const VarSourceProblem_Base& var_source;
    const VarBoundaryCondition_Base& var_boundary;
    const VarProblem<Dimension3>& var_problem;
    const HarmonicMaxwell_3D<Real_wp>& var_maxwell;

  public :
    template<class TypeEquation>
    DiffractedWaveSource_StaticMaxwell3D(const EllipticProblem<TypeEquation>& vars_helm,
                                         IncidentWaveField<Real_wp, Dimension3>&);

    Real_wp Init(const Real_wp&, const Real_wp&, int p, int n, bool);

    bool PresenceDirichlet() const;
    void EvaluateFunction(int i, int j, const R3& x, VectReal_wp& f);    
    
    void InitElement(int num_elem, const VectR3& s);
    bool IsNonNullVolumetricSource(const VectR3& s);
    
    void EvaluateVolumetricSource(int i, int j, const R3& x, VectReal_wp& f);
    
    bool IsNonNullGradientSource(const VectR3& s);
    
    void EvaluateGradientSource(int i, int j, const R3& x, VectReal_wp& f);    
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
				 const SetMatrices<Dimension3>& MatricesElem, VectReal_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  
  
  //! source for the scattering of an incident wave
  template<>
  class DiffractedWaveSource<StaticMaxwellEquation_3D>
    : public DiffractedWaveSource_StaticMaxwell3D
  {
  public :
    //! constructor with given problem
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Real_wp, Dimension3>& fsrc_);

  };

  
  //! source for the scattering of an incident wave
  /*!
    the total field is computed
   */
  class TotalWaveSource_StaticMaxwell3D : public VirtualSourceFEM<Real_wp, Dimension3>
  {
  protected:
    R3 k_wave; //!< wave vector
    R3 polar, polarH; //!< polarization vector
    Real_wp c0; //!< wave velocity
    bool element_insidePML; //!< element of a PML layer ?
    Real_wp t; //!< current time
    Real_wp offset_arg; //!< f will be evaluated for k.x - c t - offset_arg
    int num_deriv;
    bool scalar_eq;
    
  private:
    IncidentWaveField<Real_wp, Dimension3>& incident_wave;
    VirtualTimeSource<Real_wp>* wave_pulse;
    
    const VarSourceProblem_Base& var_source;
    const VarBoundaryCondition_Base& var_boundary;
    const VarProblem<Dimension3>& var_problem;
    const HarmonicMaxwell_3D<Real_wp>& var_maxwell;

  public :
    template<class TypeEquation>
    TotalWaveSource_StaticMaxwell3D(const EllipticProblem<TypeEquation>& vars_helm,
                                    IncidentWaveField<Real_wp, Dimension3>&);
    
    Real_wp Init(const Real_wp&, const Real_wp&, int p, int n, bool);
    
    void EvaluateFunction(int i, int j, const R3& x, VectReal_wp& f);
    
    void InitElement(int num_elem, const VectR3& s);
    bool IsNonNullVolumetricSource(const VectR3& s);
    
    void EvaluateVolumetricSource(int i, int j, const R3& x, VectReal_wp& f);
    
    bool IsNonNullGradientSource(const VectR3& s);
    
    void EvaluateGradientSource(int i, int j, const R3& x, VectReal_wp& f);
    
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
				 const SetMatrices<Dimension3>& MatricesElem, VectReal_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  
  
  //! source for the scattering of an incident wave
  template<>
  class TotalWaveSource<StaticMaxwellEquation_3D>
    : public TotalWaveSource_StaticMaxwell3D
  {
  public :
    template<class TypeEquation>
    TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Real_wp, Dimension3>& fsrc);
  };


  //! class for first-order absorbing boundary condition
  template<class Complexe>
  class ImpedanceABC<Complexe, StaticMaxwellEquation_3D>
    : public ImpedanceABC_Maxwell3D<Complexe, Real_wp>
  {
  public:
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);

  };
  
  
  //! modal source for 3-D static maxwell equation
  /* template<class TypeElement, class TypeSource>
  class ModalSourceBoundary<TypeElement,TypeSource,StaticMaxwellEquation_3D>
    : public ModalSourceBoundary_Maxwell<TypeElement,TypeSource,StaticMaxwellEquation_3D>
  {
  public :
    
    //! default constructor
    ModalSourceBoundary(const EllipticProblem<TypeElement, StaticMaxwellEquation_3D>& var, TypeSource& fsrc_)
      : ModalSourceBoundary_Maxwell<TypeElement,TypeSource,StaticMaxwellEquation_3D>(var, fsrc_) { }
    
      };*/

  template<class T>
  class CondensationBlockSolver<T, StaticMaxwellEquation_3D>
    : public CondensationBlockSolver_Maxwell3D<T, StaticMaxwellEquation_3D>
  {
  public:
    CondensationBlockSolver(EllipticProblem<StaticMaxwellEquation_3D>& var);
    
  };

  
  //! class used to solve static Maxwell equation with edge finite element
  template <>
  class EllipticProblem<StaticMaxwellEquation_3D>
    : public StaticMaxwell_3D<StaticMaxwellEquation_3D>
  {

  public:
    IncidentWaveProjector<Real_wp, Dimension3>*
    GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Real_wp, Dimension3>& u_inc) const;
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
    void ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
                               Vector<VectReal_wp>& curl_quadrature,
                               int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
                               bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const;

    void ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
                               Vector<VectComplex_wp>& curl_quadrature,
                               int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
                               bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const;
    
    void ComputeEnHnNodal(Vector<VectReal_wp >& u_nodal, Vector<VectReal_wp >& grad_nodal,
                          int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
                          Vector<VectReal_wp >& En_nodal, Vector<VectReal_wp >& Hn_nodal) const;

    void ComputeEnHnNodal(Vector<VectComplex_wp >& u_nodal, Vector<VectComplex_wp >& grad_nodal,
                          int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
                          Vector<VectComplex_wp >& En_nodal, Vector<VectComplex_wp >& Hn_nodal) const;
    
  };


  template<>
  class FemMatrixFreeClass<Real_wp, StaticMaxwellEquation_3D>
    : public FemMatrixFreeClass_Eq<Real_wp, StaticMaxwellEquation_3D>
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<StaticMaxwellEquation_3D>& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };

  
  template<>
  class FemMatrixFreeClass<Complex_wp, StaticMaxwellEquation_3D>
    : public FemMatrixFreeClass_Eq<Complex_wp, StaticMaxwellEquation_3D>
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<StaticMaxwellEquation_3D>& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };

  
  //! class used to solve static Maxwell equation with DG method
  template <>
  class EllipticProblem<StaticMaxwellEquation_3D_DG>
    : public StaticMaxwell_3D<StaticMaxwellEquation_3D_DG>
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
  };  


  template<>
  class VolumetricSource<StaticMaxwellEquationHdg_3D>
    : public VolumetricSource_MaxwellHdg3D<Real_wp>
  {
  public :    
    VolumetricSource(const EllipticProblem<StaticMaxwellEquationHdg_3D>&,
		     const Vector<VectString>& param);
    
  };

  template<class T>
  class CondensationBlockSolver<T, StaticMaxwellEquationHdg_3D>
    : public CondensationBlockSolver_MaxwellHdg3D<T, Real_wp>
  {
  public:
    CondensationBlockSolver(EllipticProblem<StaticMaxwellEquationHdg_3D>& var);
    
  };

  
  template<>
  class FemMatrixFreeClass<Real_wp, StaticMaxwellEquationHdg_3D>
    : public FemMatrixFreeClass_Eq<Real_wp, StaticMaxwellEquationHdg_3D>
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<StaticMaxwellEquationHdg_3D>& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };

  
  template<>
  class FemMatrixFreeClass<Complex_wp, StaticMaxwellEquationHdg_3D>
    : public FemMatrixFreeClass_Eq<Complex_wp, StaticMaxwellEquationHdg_3D>
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<StaticMaxwellEquationHdg_3D>& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };

  
  //! class used to solve static Maxwell equation with HDG method
  template <>
  class EllipticProblem<StaticMaxwellEquationHdg_3D>
    : public StaticMaxwell_3D<StaticMaxwellEquationHdg_3D>
  {
  public:
    EllipticProblem();
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };  

}

#define MONTJOIE_FILE_STATIC_MAXWELL_3D_HXX
#endif

