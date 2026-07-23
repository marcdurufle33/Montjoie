#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_HXX

namespace Montjoie
{
  //! base class for transmission conditions in Montjoie
  /*!
    For a surface Gamma, one may impose transmission conditions 
    (e.g. [u] = f and [du/dn] = g where [] is the jump across Gamma)    
   */
  template<class Dimension>
  class VarTransmission_Base
  {
  protected :
    typedef Vector<typename Dimension::R_N> VectR_N;
    typedef typename Dimension::MatrixN_N MatrixN_N;

    //! list of references for the initial surface Gamma
    IVect transmission_references;
    //! order for transmission condition
    int OrdreTransmission ;
  
    //! radius of the circle/cylinder where the condition is set
    //! (if Gamma is a circle/cylinder)
    Real_wp R0 ;
    //! number of sections (exact problem is a cyclic problem)
    //! delta = 2 pi/nb_sections
    int nb_sections;
    //! reference attributed to Gamma-
    int RefGammaMoins;
    //! reference attributed to Gamma+
    int RefGammaPlus;
    //! reference for Omega- (so that we can orient normales)
    int RefOmegaMoins;
    //! reference for Omega+ (so that we can orient normales)
    int RefOmegaPlus;
    
    //! reference attributed to the additional layer 
    /*!
      if exact solution is required, one may ask to add layers in a small layer
      placed between the two domains Omega- and Omega+. ref_omega_layer
      is the reference of the elements contained in this layer.
    */
    int ref_omega_layer;
    //! if true, a layer of small elements will be inserted between Omega- and Omega+
    bool exact_solution;
    //! Gamma- is shifted by -Decalage*delta, and Gamma+ by +Decalage*delta
    /*!
      If Decalage is set to 0, the two surfaces are the same (the vertex number and edge numbers
      are different so that the solution is discontinuous)
      If Decalage is set to 0.5, the surfaces are distant with a distance equal to delta
     */
    Real_wp Decalage;
    //! thickness of the layer
    Real_wp delta ;
    //! mesh of the initial surface Gamma
    Mesh<Dimension> mesh_surf;
    //! 2-D mesh of the surface (if it is a plate)
    Mesh<Dimension2> mesh_plate;
    MeshNumbering<Dimension2> mesh_plate_num;
    //! if true, Gamma is a plate
    bool gamma_is_plate;
    //! if true, Gamma is a sphere
    bool gamma_is_sphere;
    //! face numbers for Gamma- and Gamma+ (exact matching between faces)
    IVect NumFaceMinus, NumFacePlus;
    //! dof numbers associated with Gamma- and Gamma+ (exact matching)
    IVect DdlVolMinus, DdlVolPlus, DdlSurfLambda;
    //! global number of dofs on the surface
    int nb_ddl_volume_on_surface;
    // local to global numbering for surface
    Vector<IVect> NodleSurf, NodleAcous, NodleMeca, NodleVolMinus, NodleVolPlus;
    //! for each dof on gamma, the processor associated with the associated volume dof 
    //! (the processor can be different for Omega- and Omega+)
    IVect ProcVolMinus, ProcVolPlus;
    //! processors that communicate with the current processor
    //! dofs for lambda are shared with other processors
    IVect MatchingNumber_Subdomain;
    //! for each processor p, the dof numbers that are shared with processor p
    Vector<IVect> MatchingDofOrig_Subdomain;
    //! weighted jacobian on the surface
    Vector<VectReal_wp> JacobianAcous;
    //! dofs on the boundary
    IVect DirichletAcous;
    IVect OldNumFaceSurf, OldNumVertexSurf, IndexVert; // temporary variable;
    //! if true extremities of Gamma are duplicated (useful for periodic conditions)
    bool duplicate_extremity;
    
    //! geometrical quantitites for Gamma
    Vector<VectR_N> NormaleNodes, PosNodes;
    Vector<VectReal_wp> DsNodes;
    Vector<Vector<Vector<MatrixN_N> > > MatDfjNodes;
    
    //! variables used to retrieve jacobian matrices
    int nb_elt_local, nb_nodes_local;
    IVect nb_pts_quad_local;
    VectReal_wp coef_dfj_quad_local;

  private:
    void InitDefaultValues();
    
    DistributedProblem<Dimension>& var_problem;
    VarBoundaryCondition_Base& var_boundary;

  public:
    
    template<class TypeEquation>
    VarTransmission_Base(EllipticProblem<TypeEquation>&);
    
    // Inline methods
    virtual void InitTransmission();
    virtual void UpdateNumberOfDofs(int& nodl, int& nodl_scalar);
    virtual void UpdateGlobalDofs();
    virtual bool IsSymmetricCondition() const;
    
    const MatrixN_N& GetDfjMatrixQuadraturePoint(int i, int n, int k) const;

    const IVect& GetDofNumberOmegaMinus() const;
    const IVect& GetDofNumberOmegaPlus() const;
    
    int GetNbVolumeDofOnSurface() const;
    int GetNbInterfaces() const;
    
    void SetReferences(const IVect& ref);
    const Mesh<Dimension2>& GetMeshPlate() const;
    void SetMeshPlate(const Mesh<Dimension2>& mesh_p);
    void SetMeshPlateNumbering(const MeshNumbering<Dimension2>& mesh_n);
    void DuplicateExtremities(bool ext = true);

    Vector<IVect>& GetNodleAcous();
    IVect& GetDirichletAcous();
    IVect& GetProcVolMinus();
    IVect& GetProcVolPlus();
    Vector<VectReal_wp>& GetJacobianAcous();
    Vector<IVect>& GetNodleVolMinus();
    Vector<IVect>& GetNodleVolPlus();
    
    // other methods
    virtual size_t GetMemorySize() const;
    virtual void SetInputData(const string&, const VectString& param);
    
    void SetModifiedColNumbers(int offset);
    void SetModifiedRowNumbers(int offset);
    
    // methods to handle dofs on the two sides of the transmission layer
    void PartMeshTransmission();        
    virtual void TreatTransmission(const IVect& Epart);

    virtual void SendTransmissionDofs(const IVect&, int& nb_surfaces_involved, int& nb_loc_ddl_on_surface,
				      IVect& InfoSurfInvolved, IVect& InfoSurfMatching);

    virtual void DistributeTransmissionDofs(int nb_surfaces_involved, int nb_loc_ddl_on_surface,
					    IVect& InfoSurfInvolved, IVect& InfoSurfMatching);
    
    void Clear();
    
    // overloaded method
    virtual void AddTransmissionTerms(const Real_wp& alpha,
				      const GlobalGenericMatrix<Real_wp>& nat_mat,
				      VirtualMatrix<Real_wp>& mat_sp, int, int);

    virtual void AddTransmissionTerms(const Complex_wp& alpha,
				      const GlobalGenericMatrix<Complex_wp>& nat_mat,
				      VirtualMatrix<Complex_wp>& mat_sp, int, int);
    
    void CopyArray(IVect& ddl_vol_minus, IVect& ddl_vol_plus, 
		   Vector<IVect>& NodleAcous_, Vector<IVect>& NodleMeca_,
		   Vector<VectReal_wp>& JacobianAcous_, IVect& DirichletAcous_);
    
  };

  
  template<class TypeEquation>
  class VarTransmission : public VarTransmission_Base<typename TypeEquation::Dimension>
  {
  public:
    VarTransmission(EllipticProblem<TypeEquation>& var);
    
  };
  
} // end namespace

#define MONTJOIE_FILE_TRANSMISSION_MODEL_HXX
#endif

  
