#ifndef MONTJOIE_FILE_VAR_ADVECTION_HXX

namespace Montjoie
{

  template<class TypeEquation>
  class VarAdvection_Base : public VarHelmholtz_Cplx<TypeEquation>
  {
  public :
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeEquation::Dimension Dimension;

    void ModifyOutputUnknown(Vector<Real_wp>& val_u,
			     Vector<Real_wp>& grad_u, int i,
			     const GridInterpolation<Dimension>& var_interp,
			     int iquad, bool compute_grad) const {}

    void ModifyOutputUnknown(Vector<Complex_wp>& val_u,
			     Vector<Complex_wp>& grad_u, int i,
			     const GridInterpolation<Dimension>& var_interp,
			     int iquad, bool compute_grad) const {}
    
    void ModifyOutputUnknown(Vector<VectReal_wp>& Unodal, Vector<VectReal_wp>& GradNodal,
			     int i, bool u_component, bool grad_component) const {}
    
    void ModifyOutputUnknown(Vector<VectComplex_wp>& Unodal, Vector<VectComplex_wp>& GradNodal,
			     int i, bool u_component, bool grad_component) const {}
    
  };
  
  
  template<class T, class Dim>
  class AdvectionEquation_Base :  public GenericEquation<T>
  {
  public :
    typedef Dim Dimension;
    
    enum {nb_unknowns = 1, nb_unknowns_scal = 1,
	  nb_components_en = 1, nb_components_hn = 1, nb_unknowns_vec = 0, nb_unknowns_hdg=1};
    
    static const bool FirstOrderFormulation = true;

    // Inline functions
    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();
    static bool ComputeDFjm1();

    static const Vector<TinyMatrix<T, Symmetric,
				   Dimension::dim_N, Dimension::dim_N> >&
    GetStiffCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int num_elem, int ref);
    
    static const Vector<T>&
    GetMassCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int num_elem, int ref);
    
    // Other functions
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                              int ref, MatMass& mass);

  };

  template<class Dimension>
  class AdvectionEquation : public AdvectionEquation_Base<Real_wp, Dimension>
  {
  public :
  };


  template<class Dimension>
  class HarmonicAdvectionEquation : public AdvectionEquation_Base<Complex_wp, Dimension>
  {
  public :    
  };
  
  //! class for absorbing boundary condition and advection equation
  template<class Complexe, class Dimension>
  class ImpedanceABC<Complexe, AdvectionEquation<Dimension> > 
    : public ImpedanceFunction_Base<Complexe, Dimension>
  {

  private:
    const VarHelmholtz_Base<Complexe, Dimension>& var_helm;
    const VarProblem<Dimension>& var_problem;

  public :
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);

    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<Complexe>& nat_mat, int ref,
                              const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat);
    
  };


  //! class for absorbing boundary condition and advection equation
  template<class Dimension>
  class ImpedanceABC<Complex_wp, HarmonicAdvectionEquation<Dimension> > 
    : public ImpedanceFunction_Base<Complex_wp, Dimension>
  {
  private:
    const VarHelmholtz_Base<Complex_wp, Dimension>& var_helm;
    const VarProblem<Dimension>& var_problem;

  public :
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);

    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
                              const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat);
    
  };
  

  //! class to solve advection equation with H1 elements
  template<class Dimension>
  class EllipticProblem<AdvectionEquation<Dimension> >
    : public VarLaplace_Base<AdvectionEquation<Dimension> >
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


  //! class to solve advection equation with H1 elements in time-harmonic domain
  template<class Dimension>
  class EllipticProblem<HarmonicAdvectionEquation<Dimension> >
    : public VarAdvection_Base<HarmonicAdvectionEquation<Dimension> >
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


  template<class T, class Dim>
  class AdvectionEquationDG_Base : public AdvectionEquation_Base<T, Dim>
  {
  public :
    
    template<class Matrix1, class T0, class GenericPb>
    static void GetNabc(Matrix1& Nabc, typename Dim::R_N& normale,
			int ref, int iquad, int k,const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
			const GenericPb& vars, const ElementReference<Dim, 1>& Fb);
    
    template<class Matrix1, class T0, class GenericPb>
    static void GetPenalDG(Matrix1& Nabc, typename Dim::R_N& normale,
			   int iquad, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
			   const GenericPb& vars, const ElementReference<Dim, 1>& Fb);
    
  };

  template<class Dimension>
  class AdvectionEquationDG : public AdvectionEquationDG_Base<Real_wp, Dimension>
  {
  public :
  };


  template<class Dimension>
  class HarmonicAdvectionEquationDG : public AdvectionEquationDG_Base<Complex_wp, Dimension>
  {
  public :
    
  };


  //! class to solve advection equation with DG elements
  template<class Dimension>
  class EllipticProblem<AdvectionEquationDG<Dimension> >
    : public VarLaplace_Base<AdvectionEquationDG<Dimension> >
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


  //! class to solve advection equation with DG elements
  template<class Dimension>
  class EllipticProblem<HarmonicAdvectionEquationDG<Dimension> >
    : public VarAdvection_Base<HarmonicAdvectionEquationDG<Dimension> >
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


  template<class Dimension>
  class TimeAdvectionEquation
  {
  public :
    typedef AdvectionEquation<Dimension> TypeEquationStationary; //!< stationary equation
    
  };

  template<class Dimension>
  class TimeAdvectionEquationDG : public TimeAdvectionEquation<Dimension>
  {
  public :
    typedef AdvectionEquationDG<Dimension> TypeEquationStationary; //!< stationary equation
    
  };
  
    
  //! class to solve advection equation with H1 elements
  template<class Dimension>
  class HyperbolicProblem<TimeAdvectionEquation<Dimension> >
    : public VarInstationary<TimeAdvectionEquation<Dimension> >
  {
  };

  
  //! class to solve advection equation with discontinuous Galerkin
  template<class Dimension>
  class HyperbolicProblem<TimeAdvectionEquationDG<Dimension> >
    : public VarInstationary<TimeAdvectionEquationDG<Dimension> >
  {
  };  



  template<class T, class Dim>
  class DoubleAdvectionEquation_Base :  public GenericEquation<T>
  {
  public :
    typedef Dim Dimension;
    
    enum {nb_unknowns = 2, nb_unknowns_scal = 2,
	  nb_components_en = 2, nb_components_hn = 2, nb_unknowns_vec = 0, nb_unknowns_hdg = 1};

    static const bool FirstOrderFormulation = true;
    
    // Inline functions
    static bool ComputeDFjm1();

    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();

    static const Vector<TinyMatrix<T, Symmetric,
				   Dimension::dim_N, Dimension::dim_N> >&
    GetStiffCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int num_elem, int ref);
    
    static const Vector<T>&
    GetMassCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int num_elem, int ref);

    // other functions
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                              int ref, MatMass& mass);

    template<class Matrix1, class T0, class GenericPb>
    static void GetNabc(Matrix1& Nabc, typename Dim::R_N& normale,
			int ref, int iquad, int k,const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
			const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);

  };


  template<class Dimension>
  class DoubleAdvectionEquation : public DoubleAdvectionEquation_Base<Real_wp, Dimension>
  {
  public :
  };


  template<class Dimension>
  class HarmonicDoubleAdvectionEquation : public DoubleAdvectionEquation_Base<Complex_wp, Dimension>
  {
  public :
    
  };
  
  //! class for absorbing boundary condition and advection equation
  template<class Complexe, class Dimension>
  class ImpedanceABC_DoubleAdvec
    : public ImpedanceFunction_Base<Complexe, Dimension>
  {
  private:
    const VarHelmholtz_Base<Complexe, Dimension>& var_helm;
    const VarProblem<Dimension>& var_problem;

  public :
    template<class TypeEquation>
    ImpedanceABC_DoubleAdvec(const EllipticProblem<TypeEquation>& var);

    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<Complexe>& nat_mat, int ref,
                              const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat);
    
  };


  //! class for absorbing boundary condition and advection equation
  template<class Complexe, class Dimension>
  class ImpedanceABC<Complexe, DoubleAdvectionEquation<Dimension> >
    : public ImpedanceABC_DoubleAdvec<Complexe, Dimension>
  { 
  public:
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);
  };


  //! class for absorbing boundary condition and advection equation
  template<class Dimension>
  class ImpedanceABC<Complex_wp, HarmonicDoubleAdvectionEquation<Dimension> > 
    : public ImpedanceABC_DoubleAdvec<Complex_wp, Dimension>
  {
  public:
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);
    
  };

  
  //! class to solve advection equation with H1 elements
  template<class Dimension>
  class EllipticProblem<DoubleAdvectionEquation<Dimension> >
    : public VarLaplace_Base<DoubleAdvectionEquation<Dimension> >
  {
    public :
    EllipticProblem();
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };  

  //! class to solve advection equation with H1 elements in time-harmonic domain
  template<class Dimension>
  class EllipticProblem<HarmonicDoubleAdvectionEquation<Dimension> >
    : public VarAdvection_Base<HarmonicDoubleAdvectionEquation<Dimension> >
  {
  public :
    EllipticProblem();
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };  

 
  template<class T, class Dim>
  class DoubleAdvectionEquationDG_Base : public DoubleAdvectionEquation_Base<T, Dim>
  {
  public :
    
    template<class Matrix1, class T0, class GenericPb>
    static void GetPenalDG(Matrix1& Nabc, typename Dim::R_N& normale,
			   int iquad, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
			   const GenericPb& vars, const ElementReference<Dim, 1>& Fb);
    
  };

  template<class Dimension>
  class DoubleAdvectionEquationDG : public DoubleAdvectionEquationDG_Base<Real_wp, Dimension>
  {
  public :
  };


  template<class Dimension>
  class HarmonicDoubleAdvectionEquationDG : public DoubleAdvectionEquationDG_Base<Complex_wp, Dimension>
  {
  public :
    
  };

  
  //! class to solve advection equation with DG elements
  template<class Dimension>
  class EllipticProblem<DoubleAdvectionEquationDG<Dimension> >
    : public VarLaplace_Base<DoubleAdvectionEquationDG<Dimension> >
  {
    public :
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


  //! class to solve advection equation with DG elements
  template<class Dimension>
  class EllipticProblem<HarmonicDoubleAdvectionEquationDG<Dimension> >
    : public VarAdvection_Base<HarmonicDoubleAdvectionEquationDG<Dimension> >
  {
    public :
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


  template<class Dimension>
  class TimeDoubleAdvectionEquation
  {
  public :
    typedef DoubleAdvectionEquation<Dimension> TypeEquationStationary; //!< stationary equation

  };

  template<class Dimension>
  class TimeDoubleAdvectionEquationDG : public TimeDoubleAdvectionEquation<Dimension>
  {
  public :
    typedef DoubleAdvectionEquationDG<Dimension> TypeEquationStationary; //!< stationary equation
    
  };
  
    
  //! class to solve advection equation with H1 elements
  template<class Dimension>
  class HyperbolicProblem<TimeDoubleAdvectionEquation<Dimension> >
    : public VarInstationary<TimeDoubleAdvectionEquation<Dimension> >
  {
    public :
    
  };

  
  //! class to solve advection equation with discontinuous Galerkin
  template<class Dimension>
  class HyperbolicProblem<TimeDoubleAdvectionEquationDG<Dimension> >
    : public VarInstationary<TimeDoubleAdvectionEquationDG<Dimension> >
  {
    public :
  };  

}

#define MONTJOIE_FILE_VAR_ADVECTION_HXX
#endif
