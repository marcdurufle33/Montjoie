#ifndef MONTJOIE_FILE_COMPUTATION_MODE_MAXWELL_HXX

namespace Montjoie
{
  
  //! equation used to compute modes of 3-D, computation of 2-D eigenvalues
  template<class T>
  class ModeMaxwellEquation : public GenericEquation<T>
  {
  public :
    typedef Dimension2 Dimension;
    
    enum {nb_unknowns = 2, nb_unknowns_vec = 2, nb_unknowns_scal = 2,
          nb_components_en = 1, nb_components_hn = 1, nb_unknowns_hdg=0,
          type_element = 2};
    
    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return false; }

    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int num_elem, const ElementReference_Dim<Dimension2>& Fb);

    template<class TypeEquation>
    static void ComputeMassMatrixQuad(EllipticProblem<TypeEquation>& var,
				      int num_elem, const QuadrangleHcurlFirstFamily& Fb);
    
  };
  
  
  //! base class to compute modes in 3-D guides
  template<class TypeEquation>
  class VarModeMaxwell : public VarHarmonic<TypeEquation>
  {
  public :
    typedef typename TypeEquation::Complexe Complexe;
    
    //! array used to compute elementary matrices
    Vector<Vector<Real_wp> > Glob_matMass_Dh, Glob_matMassPoids;
    //! array used to compute elementary matrices
    Vector<Vector<TinyMatrix<Real_wp, Symmetric, 2, 2> > > Glob_matMass_Bh;
    //! physical indices epsilon, mu, sigma
    Vector<ScalarPhysicalIndice<Dimension2, Complexe> > ref_epsilon, ref_mu, ref_sigma;

    template<class T0, class Prop>
    void GetMassMatrixCoef(int num_elem, int j, const GlobalGenericMatrix<T0>& nat_mat,
			   int m, int n, TinyMatrix<T0, Prop, 2, 2>& mass) const;
    
    template<class T0>
    void GetStiffMatrixCoef(int num_elem, int j,
			    const GlobalGenericMatrix<T0>& nat_mat, int m, int n, T0& stiff) const;
    
    void InitIndices(int n);
    int GetNbPhysicalIndices() const;
    
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);    
    string GetPhysicalIndexName(int m) const;
    
    void UpdateShiftAdimensionalization(Real_wp& shift, Real_wp& shift_imag);
    void UpdateShiftAdimensionalization(Complex_wp& shift, Complex_wp& shift_imag);
        
    void UpdateEigenvaluesAdimensionalization(Vector<Real_wp>& eigen_values,
                                              Vector<Real_wp>& lambda_imag,
					      Matrix<Real_wp, General, ColMajor>& eigen_vectors);

    void UpdateEigenvaluesAdimensionalization(Vector<Complex_wp>& eigen_values,
                                              Vector<Complex_wp>& lambda_imag,
					      Matrix<Complex_wp, General, ColMajor>& eigen_vectors);

    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
			   IVect& num_ref, IVect& num_index, IVect& num_component,
			   Vector<bool>& compute_grad, Vector<bool>& compute_hess);
    
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;
    
    void AllocateMassMatrices();
    
  };
  

  //! class to compute Maxwell modes in 3-D guides
  template<>
  class EllipticProblem<ModeMaxwellEquation<Real_wp> >
    : public VarModeMaxwell<ModeMaxwellEquation<Real_wp> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };


  //! class to compute Maxwell modes in 3-D guides
  template<>
  class EllipticProblem<ModeMaxwellEquation<Complex_wp> >
    : public VarModeMaxwell<ModeMaxwellEquation<Complex_wp> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };


#ifdef MONTJOIE_WITH_THREE_DIM
  //! class to compute modes in a section for 3-D Maxwell's equations  
  class ModalSourceMaxwell3D_Base
  {
  protected:
    //! reference for each boundary in order to solve Poisson problem
    IVect ref_mode_poisson; 
    //! values of constant potentials for Poisson problem
    VectReal_wp coef_mode_poisson;
    
    int mode_to_compute;
    enum{MODE_POISSON, MODE_MAXWELL};

  private:
    const DistributedProblem<Dimension3>& var_problem;
    const ModalSourceBoundary_Base& modal_source;
    
  public:
    template<class TypeEquation>
    ModalSourceMaxwell3D_Base(const EllipticProblem<TypeEquation>& var,
			      const ModalSourceBoundary_Base& mode);

    void GetCoefficientBoundaryCondition(const VectReal_wp& eigenval, VectReal_wp& coef_eigenval,
					 const Real_wp& imped, int boundary_condition);
    
    void GetCoefficientBoundaryCondition(const VectReal_wp& eigenval, VectComplex_wp& coef_eigenval,
					 const Complex_wp& imped, int boundary_condition);
    
    void EvaluateModePoisson(const Mesh<Dimension3>& result_mesh, int proc_mode,
			     const Vector<int>& IndexBoundary, const Vector<int>& IndexElement,
			     const IVect& RotElement, const Vector<VectR2>& LocalQuadPoints,
			     Vector<VectR3>& eval_mode);
    
    void EvaluateModeMaxwell(const Mesh<Dimension3>& result_mesh, int proc_mode,
			     const Vector<int>& IndexBoundary, const Vector<int>& IndexElement,
			     const IVect& RotElement, const Vector<VectR2>& LocalQuadPoints,
			     Vector<Vector<VectR3> >& eval_mode, VectReal_wp& eigenval);
      
    virtual void CopyIndices(VarHelmholtz_Base<Real_wp, Dimension2>& var_eig) = 0;    
    virtual void CopyIndices(VarModeMaxwell<ModeMaxwellEquation<Real_wp> >& var_eig) = 0;
    
  };

  
  //! class to compute modes in a section for 3-D Maxwell's equations  
  template<class Complexe>
  class ModalSourceMaxwell3D : public ModalSourceBoundary_Dim<Complexe, Dimension3>,
			       public ModalSourceMaxwell3D_Base
  {
  protected:
    Complexe voltage;
    
  private:
    const DistributedProblem<Dimension3>& var_problem;
    const HarmonicMaxwell_3D<Complexe>& var_maxwell;
    
  public:
    template<class TypeEquation>
    ModalSourceMaxwell3D(const EllipticProblem<TypeEquation>& var);
    
    void ReadSourceParameters(const IVect& boundary_condition, const IVect& ref,
			      const VectString& source_param);

    void EvaluateMode(const Mesh<Dimension3>& mesh, int proc_mode, const Vector<int>& IndexBoundary, 
		      const Vector<int>&, const Vector<int>& RotElement,
		      const Vector<VectR2>& LocalQuadPoints,
		      Vector<Vector<Vector<Complexe> > >& eval);
    
    void CopyIndices(VarHelmholtz_Base<Real_wp, Dimension2>& var_eig);
    void CopyIndices(VarModeMaxwell<ModeMaxwellEquation<Real_wp> >& var_eig);
    
  };
#endif

}

#define MONTJOIE_FILE_COMPUTATION_MODE_MAXWELL_HXX
#endif

