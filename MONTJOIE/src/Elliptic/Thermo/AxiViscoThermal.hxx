#ifndef MONTJOIE_FILE_AXI_VISCO_THERMAL_HXX

namespace Montjoie
{
  // visco-thermal equation in axisymmetric domain
  template<class T>
  class ViscoThermalEquationAxi : public GenericEquation<T>
  {
  public:
    typedef Dimension2 Dimension;

    enum{nb_unknowns = 4, nb_components_en = 1, nb_components_hn = 1,
         nb_unknowns_scal = 4, nb_unknowns_vec = 0 , nb_unknowns_hdg = 0};

    static const bool FormulationDG = false;
    static const bool FirstOrderFormulation = false;
    static const bool DiscontinuousDiMatrix = false;
    static const bool TensorStiffnessSymmetric = false;

    static inline IVect GetOtherElementType() { IVect num(1); num(0) = 1; return num; }

    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return true; }

    // providing C in Cgrad_grad
    template<class TypeEquation, class NatureMat, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                  int num_elem, int jloc, const NatureMat& nat_mat,
                                  int ref, MatStiff& Cgrad_grad);

    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc,
                                 const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);

    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                              int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);

  };

  
  //! class for axisymmetric visco-thermal equations
  template<class TypeEquation>
  class VarViscoThermal_Axi : public VarAxisymProblem, public VarHarmonic<TypeEquation>
  {
    const VarProblem<Dimension2>& var_problem;
    
  public:    
    //! radius on quadrature points
    Vector<Vector<Real_wp> > Glob_radius;

    Real_wp c_p, T0, P0, rho0, mu, lambda, zeta;
    bool alternative_ipp;
    
  public:
    VarViscoThermal_Axi(EllipticProblem<TypeEquation>& var);

    inline VarViscoThermal_Axi() {c_p = 1.0; T0 = 1.0; P0 = 1.0; rho0 = 1.0; mu = 1.0; lambda = 1.0; zeta = 1.0;};

    void AllocateMassMatrices();
    
    void ComputeLocalMassMatrix(int i, int N, bool linear_sparse,
				SetPoints<Dimension2>& PointsElem,
                                SetMatrices<Dimension2>& MatricesElem,
                                IVect& OrderFace, const ElementGeomReference<Dimension2>& Fb);
    
    inline bool UseNumericalIntegration(int i) const { return true; }
    void UpdateDirichlet(int n);    
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;

    inline void InitIndices(int n) {}
    inline int GetNbPhysicalIndices() const{ return 100;}
    inline void SetIndices(int i, const VectString& parameters) {}
    inline bool IsVaryingMedia(int i) const { return false; }
    inline bool IsVaryingMedia(int m, int i) const { return false; }
    inline Real_wp GetVelocityOfMedia(int ref) const{ return 1.0;}
    inline Real_wp GetVelocityOfInfinity() const{return 1.0;}

    inline void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
                                  IVect& num_ref, IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                                  Vector<bool>& compute_hess) {}

    inline int GetNbModesSource() const { return 1; }
    inline const R3& GetPhaseOrigin() const { return R3(0, 0, 0); }
    
  };


  //! volumetric source
  template<>
  class VolumetricSource<ViscoThermalEquationAxi<Complex_wp> > : public VirtualSourceFEM<Complex_wp, Dimension2>
  {
  private:
    const EllipticProblem<ViscoThermalEquationAxi<Complex_wp> >& var_problem;
    
  protected:
    VirtualSourceField<Complex_wp, Dimension3>* fsrc;
    Vector<bool> unif_source_surf, gaussian_source_surf;
    Vector<Vector<Complex_wp> > coef_source_surf;
    
  public :
    template<class TypeEquation>
    VolumetricSource(const EllipticProblem<TypeEquation>& var,
                     const Vector<VectString>& param);
    
    ~VolumetricSource();

    void SetSurfaceSource(int ref, VirtualSourceField<Complex_wp, Dimension3>* f);
    
    bool IsNonNullSurfacicSource(int ref);
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
				const SetMatrices<Dimension2>& MatricesElem, Vector<Complex_wp>& f);

    void EvaluateFunction(int i, int j, const R2& x, Vector<Complex_wp>& f);    
    void EvaluateMode(const R2& x, Vector<Complex_wp>& f);
    
  };
  
  
  //---------------------------
  // EllipticProblem
  //---------------------------
  template<>
  class EllipticProblem<ViscoThermalEquationAxi<Complex_wp> >
    : public VarViscoThermal_Axi<ViscoThermalEquationAxi<Complex_wp> >
  {
  public:
    inline EllipticProblem<ViscoThermalEquationAxi<Complex_wp> >()
      : VarViscoThermal_Axi<ViscoThermalEquationAxi<Complex_wp> > (*this)
    {
    }       
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                                 CondensationBlockSolver_Base<Real_wp>&,
                                 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                                 CondensationBlockSolver_Base<Complex_wp>&,
                                 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
  };

}

#define MONTJOIE_FILE_AXI_VISCO_THERMAL_HXX
#endif
