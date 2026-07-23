#ifndef MONTJOIE_FILE_VAR_PROBLEM_1D_HXX

namespace Montjoie
{

  //! base class for PML damping function
  class PmlDampingFunction
  {
  protected:
    //! damping factor in PML
    Real_wp vsigma, exponent_pml;
    //! thickness of PML
    Real_wp thickness;
    //! PML is outside the interval [xmin, xmax]
    Real_wp xmin, xmax;
    //! type of damping in PML
    int type_function;
    //! maximum velocity in PML layers
    Real_wp max_velocity;
    //! available damping functions
    enum {PML_PARABOLE, PML_SHIFTED_PARABOLE, PML_CONSTANT, PML_EXPONENT};
    //! offset for damping in PML
    Real_wp offset_damping;

    void EvaluateFunctionTau(const Real_wp& dx, const Real_wp& dsig, const Real_wp& a,
			     Real_wp& zeta, Real_wp& zeta_prime) const;
    
  public:
    PmlDampingFunction();
    
    void SetInputData(const Vector<string>& parameters);

    void SetMaximumVelocity(const Real_wp& v);
    void SetPmlArea(const Real_wp& x0, const Real_wp& xN, const Real_wp& d);
    
    void GetDampingTau(const Real_wp& point, Real_wp& zeta, Real_wp& zeta_prime) const;
    
  };
  

  //! class for 1-D linear problems
  class VarProblem_1D : public VarFiniteElement<Dimension1>, public InputDataProblem_Base,
			public VarComputationProblem_Base, public VarPhysicalProblem
  {
  protected:
    int type_coordinate;
    enum { CARTESIAN, POLAR, SPHERICAL};
    
    int nb_unknowns;
    int condition_left_side; //!< left boundary condition
    int condition_right_side; //!< rigt boundary condition

    GridInterpolation<Dimension1> var_section; //!< interpolation grid
    ParamOutputClass output_section_param; //!< parameters for output in interpolation grid
    
    IVect OffsetDofV; //!< offsets for dofs for the unknown v (first-order formulation)
    int nodl; //!< total number of degrees of freedom
    
    Real_wp omega, omega2; //!< pulsation
    
    Vector<VectString> mesh_data; //!< parameters for the 1-D mesh
    
    IVect Dirichlet_dof; //!< dofs with Dirichlet condition
    
    int print_level; //!< verbose level    
    Real_wp threshold_matrix; //!< values below threshold are skipped
    bool mixed_formulation; //!< if true, the first-order formulation is used
    //! envelope to compute ?
    bool calcul_enveloppe;

    //! extremities of the physical domain
    Real_wp xmin_, xmax_;
    VectReal_wp CoorDofs;
    
    //! global number of dofs
    int nodl_mesh;
    
    //! for parallel stuff
    IVect OverlapDofNumber_Subdomain;

  public :
    //! PML damping function
    PmlDampingFunction pml_damping;
    
    Mesh<Dimension1> mesh; //!< 1-D mesh

#ifdef SELDON_WITH_MPI
    // processors sharing the same mesh
    MPI_Comm comm_group_mode;
#endif

    
  public :
    VarProblem_1D();
    
    // Inline methods
    void GetMiomega(Real_wp& m_iomega);
    void GetMiomega(Complex_wp& m_iomega);

    void SetOmega(const Real_wp& om);
    
    bool FirstOrderFormulation() const;
    void SetFirstOrderFormulation(bool first = true);
    Real_wp GetXmin() const;
    Real_wp GetXmax() const;
    
    int GetNbDofBoundaries(int i);
    int GetNbDof() const;
    int GetNbRows() const;
    int GetNbElt() const;
    int GetPrintLevel() const;
    
    VectString& GetMeshData(int n);    
    void SetMeshData(const Vector<VectString>& param);
    
    GridInterpolation<Dimension1>& GetSectionGrid();
    const GridInterpolation<Dimension1>& GetSectionGrid() const;
    ParamOutputClass& GetOutputParameters();
    const ParamOutputClass& GetOutputParameters() const;

    const VectReal_wp& GetCoordinateDof() const;
    const Real_wp& GetCoordinateDof(int i) const;
    const IVect& GetOverlappedDofNumber() const;

    template<class MatrixSparse>
    void SetPeriodicCondition(MatrixSparse& A);
    
    template<class MatrixSparse>
    void SetDirichletCondition(MatrixSparse& A);

    int GetLeftBoundaryCondition() const;
    int GetRightBoundaryCondition() const;

    bool EnvelopeToCompute() const;
    
    void SetBoundaryCondition(int left_cond, int right_cond);

    void SetPulsation(const Real_wp& om);

    void GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
				 Vector<int>& intern_node) const;
    
    // Other methods
    void SetInputData(const string& description_field, const VectString& parameters);
    virtual int ReadBoundaryCondition(int& nb, const VectString& parameters, int side);
    
    void ComputeQuadraturePoints(VectReal_wp& pts);
    void GetIntegralMass(Vector<Real_wp>& M);
    
    template<class MatrixB>
    void GetStiffnessMatrix(MatrixB& mat_stiff);

    template<class MatrixB>
    void GetGradientMatrix(MatrixB& mat_stiff);

    template<class VectorB>
    void GetLocalMassMatrix(VectorB& mat_mass);
    
    void EvaluateFunctionTauPML(const Real_wp& dx, const Real_wp& dsig, const Real_wp& a,
				Real_wp& zeta, Real_wp& zeta_prime) const;
    
    void GetDampingTauPML(const Real_wp& point, Real_wp& zeta, Real_wp& zeta_prime) const;
    
    void ComputeMeshAndFiniteElement(const string& name_elt, bool split_mesh = true);

    template<class T>
    void AddDomains(Vector<T>& X, int nb_u = -1) const;
    
    template<class T>
    void AddVolumetricSource(Vector<T> & b_source, VirtualSourceField<T, Dimension1>& f);
    
    /* template<class T>
       void AddSurfacicSource(Vector<T>& b_source, VirtualSourceFEM<T, Dimension1>& f);

       template<class T>
       void AddVolumetricProjection(const T& alpha, Vector<T>&, VirtualProjectorFEM<T, Dimension1>& proj); */
    
    template<class T>
    void ComputeRightHandSide(Vector<T>& bsrc);
    
    void TreatDirichletCondition();
    
    template<class T>
    void AddMatrixFEM(VirtualMatrix<T>& mat_direct, GlobalGenericMatrix<T>& nat_mat);
    
    template<class T>
    void WriteDatas(Vector<T> & U0);
    
    template<class T>
    void ComputeInterpolationU(const Vector<T> & U0,
                               const GridInterpolation<Dimension1> & var_interp,
                               Vector<T>& trace_u, bool discont = false, bool grad = false) const;
    
    template<class T> 
    T GetInterpolate(const Vector<T>& U0, int iquad, const Real_wp& point_loc,
		     VectReal_wp& val_phi, bool discont = false, bool grad = false) const;

    // Virtual methods
    virtual void ReadImpedance(int&, const VectString& param, int side) = 0;
    virtual Real_wp GetVelocityOfMedia(int ref) const = 0;
    
  };

}

#define MONTJOIE_FILE_VAR_PROBLEM_1D_HXX
#endif
