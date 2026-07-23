#ifndef MONTJOIE_FILE_GENERALIZED_IMPEDANCE_MODEL_HXX

namespace Montjoie
{
  //! base class for GIBC (Generalized Impedance Boundary Condition) in Montjoie
  /*!
    For a surface Gamma, one may impose impedance boundary condition
    (e.g.  du/dn + T(u) = 0)
   */
  class VarGeneralizedImpedance_Base
  {
  protected :
    //! list of references for the initial surface Gamma
    IVect gibc_references;
    //! if ref_cond(ref) = 1, reference ref is in Gamma
    IVect ref_cond;
    
    //! order for GIBC
    int order_gibc;
    
    int nodl_mesh; 
    //! dof numbers associated with Gamma
    IVect DdlVol; // DdlSurfLambda;
    Vector<IVect> NodleVol;
    
    //! global number of dofs on the surface
    int nb_ddl_volume_on_surface;
    
    //! processors that communicate with the current processor
    //! dofs for lambda are shared with other processors
    IVect MatchingNumber_Subdomain;
    
    //! for each processor p, the dof numbers that are shared with processor p
    Vector<IVect> MatchingDofOrig_Subdomain;

    //! using unsymmetric implementation (matrix will be unsymmetric) ?
    bool use_unsymmetric_algo;

    bool pade_approx_sqrt; Real_wp epsilon_pade;
    VectComplex_wp CoefAl, CoefBl; Complex_wp CoefC0;

  private:
    VarProblem_Base& var_problem;
    VarBoundaryCondition_Base& var_boundary;

    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    VarGeneralizedImpedance_Base(EllipticProblem<TypeEquation>& var);
    
    // Inline methods
    virtual void InitGIBC();

    void InitReferences(int N);
    int GetNbVolumeDofOnSurface() const;
    bool UsePadeApproximationSquareRoot() const;
    bool UseUnsymmetricImplementation() const;

    const Vector<int>& GetLocalVolumeDofNumber() const;
    
    // Other methods
    virtual size_t GetMemorySize() const;
    void SetInputData(const string&, const VectString& param);
    void Clear();
    
    void SetModifiedColNumbers(int offset);
    void SetModifiedRowNumbers(int offset);

    // methods to retrieve dofs on the boundary
    template<class Dimension>
    void FindLocalDofs(int i, int ne, int num_loc, VarProblem<Dimension>& var,
		       IVect& IndexDof, const ElementReference_Dim<Dimension>& Fb);

    template<class Dimension>
    void FindDofsGibc(VarProblem<Dimension>& var, const IVect& Epart);

    // overloaded virtual methods
    virtual void TreatGibc(const IVect& Epart) = 0;
    
    virtual void AddGibcTerms(const Real_wp& alpha, const GlobalGenericMatrix<Real_wp>& nat_mat,
			      VirtualMatrix<Real_wp>& mat_sp, int offset_row = 0, int offset_col = 0);

    virtual void AddGibcTerms(const Complex_wp& alpha, const GlobalGenericMatrix<Complex_wp>& nat_mat,
			      VirtualMatrix<Complex_wp>& mat_sp, int offset_row = 0, int offset_col = 0);
    
  };

  
  template<class TypeEquation>
  class VarGeneralizedImpedance : public VarGeneralizedImpedance_Base
  {
  private:
    EllipticProblem<TypeEquation>& var_problem;
    
  public:
    VarGeneralizedImpedance(EllipticProblem<TypeEquation>& var);
    
    void TreatGibc(const IVect& Epart);

  };
  
} // end namespace

#define MONTJOIE_FILE_GENERALIZED_IMPEDANCE_MODEL_HXX
#endif

  
