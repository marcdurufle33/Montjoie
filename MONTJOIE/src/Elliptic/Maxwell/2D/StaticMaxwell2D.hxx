#ifndef MONTJOIE_FILE_STATIC_MAXWELL_2D_HXX

namespace Montjoie
{

  //! class to describe static maxwell equation in 2-D
  class StaticMaxwellEquation_2D : public GenericEquation<Real_wp>
  {
  public :
    typedef Dimension2 Dimension;
    
    enum {nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_hdg=1,
	  nb_components_en = 1, nb_components_hn = 1, nb_unknowns_vec = 1,
    type_element = 2};
    
    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return false; }
    
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
				  int num_elem, const ElementReference_Dim<Dimension2>& Fb);
    
    template<class TypeEquation>
    static void ComputeMassMatrixQuad(EllipticProblem<TypeEquation>& var,
				      int num_elem, const QuadrangleHcurlFirstFamily& Fb);
    
    template<class TypeEquation, class T0>
    static void GetPenalizationCoef(const EllipticProblem<TypeEquation>& var,
				    int i, const GlobalGenericMatrix<T0>& nat_mat, int ref, T0& coef);
    
  };
  
  
  //! class to describe static maxwell equation in 2-D, with DG method
  class StaticMaxwellEquation_2D_DG : public GenericEquation<Real_wp>
  {
  public :
    typedef Dimension2 Dimension;
    
    static const bool FirstOrderFormulation = true;
    
    enum {nb_unknowns = 3, nb_unknowns_hdg=1,
	  nb_components_en = 1, nb_components_hn = 1, nb_unknowns_scal = 2, nb_unknowns_vec = 1};
    
    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return false; }

    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int i, const ElementReference_Dim<Dimension>&);
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& var,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V);
    
    template<class TypeEquation, class T0, class Vector1>
    static void GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
                                    const GlobalGenericMatrix<T0>& nat_mat,
                                    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
				 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
				 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& Vn, Vector2& Un);

    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& Un, Vector2& Vn);
    
    template<class Matrix1, class R_N, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, R_N& normale, int ref, int iquad, int k, 
                        const GlobalGenericMatrix<T0>& nat_mat, int ref2,
                        const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb);

    template<class Matrix1, class R_N, class GenericPb, class T0>
    static void GetPenalDG(Matrix1& Nabc, R_N& normale, int iquad, int k, int nf,
                           const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb);

    template<class Vector1, class Vector2, class GenericPb, class T0>
    static void MltPenalDG(const R2& normale, const Vector1& Vn, Vector2& Un, int i, int k,
			   int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb);

    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(const R2& normale, int ref, const Vector1& Vn, Vector1& Un,
                        int num_elem1, int npoint, const GlobalGenericMatrix<T0>& nat_mat, int ref2,
                        const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension2, 1>& Fb);
    
  };
  
  //! base class for 2-D static Maxwell equation
  template<class TypeEquation>
  class StaticMaxwell_2D : public HarmonicMaxwell_2D_Eq<TypeEquation>
  {
    typedef typename TypeEquation::Complexe Complexe;
    
  public:
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
			   IVect& num_ref, IVect& num_index, IVect& num_component,
			   Vector<bool>& compute_grad, Vector<bool>& compute_hess);
    
  };


  //! class for first-order absorbing boundary condition
  template<class T>
  class ImpedanceABC<T, StaticMaxwellEquation_2D>
    : public ImpedanceFunction_Base<T, Dimension2>
  {
  private:
    const VarProblem<Dimension2>& var_problem;
    const HarmonicMaxwell_2D<Real_wp>& var_maxwell;
    
  public :
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);
    
    //! impedance
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<T>& nat_mat, int ref,
                              const SetPoints<Dimension2>& Pts, const SetMatrices<Dimension2>& Mat);  
    
  };
  
  
  //! class to solve 2-D static Maxwell equations with edge finite element
  template<>
  class EllipticProblem<StaticMaxwellEquation_2D> : public StaticMaxwell_2D<StaticMaxwellEquation_2D>
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);

  };

  //! class to solve 2-D static Maxwell equations with DG method
  template<>
  class EllipticProblem<StaticMaxwellEquation_2D_DG>
    : public StaticMaxwell_2D<StaticMaxwellEquation_2D_DG>
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
  
}

#define MONTJOIE_FILE_STATIC_MAXWELL_2D_HXX
#endif

