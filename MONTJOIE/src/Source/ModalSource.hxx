#ifndef MONTJOIE_FILE_MODAL_SOURCE_HXX

namespace Montjoie
{
  //! Base class for modal source 
  class ModalSourceBoundary_Base
  {
  protected :
    //! references of the surface of the mode
    IVect is_ref_on_mode;
    //! number of modes to compute (computation of eigenvalues)
    int nb_modes_to_compute;
    //! shift used to compute eigenvalues
    Real_wp shift_eigenval_mode;
    //! Boundary condition associated to lateral surfaces
    int lateral_condition;
    //! Boundary condition associated with the mode
    int boundary_condition_mode;
    //! coefficients (a,b,c,d) for definition of a plane a x + b y + c z + d = 0
    /*!
      Used when the mode is supposed to be on a plane
     */
    VectReal_wp coef_plane;
    
    //! the number of modes to combine to obtain the final source
    int nb_modes_to_combine;
    //! list of mode numbers to be combined
    IVect number_mode_combined;
    //! a reference of the volume in contact with the section
    int reference_volume;
    //! keyword associated with the mode
    string keyword_mode;
    
  private :
    void InitDefaultValues(int nb_ref_mesh);
    
  public:
    template<class TypeEquation>
    ModalSourceBoundary_Base(const EllipticProblem<TypeEquation>&);

    int GetNbModesToCompute() const;
    const Real_wp& GetEigenvalueShift() const;
    const VectReal_wp& GetCoefficientPlane() const;
    const IVect& GetModeNumberToCombine() const;
    int GetLateralBoundaryCondition() const;
    int GetVolumeReference() const;
    
#ifdef MONTJOIE_WITH_TWO_DIM  
    void GetMeshAndQuadraturePoints(const DistributedProblem<Dimension2>& var,
				    Mesh<Dimension2>& result_mesh, int proc_mode,
				    Vector<int>& IndexBoundary, Vector<int>&,
				    Vector<int>&, Vector<VectReal_wp>& LocalQuadPoints);
#endif
    
#ifdef MONTJOIE_WITH_THREE_DIM  
    void GetMeshAndQuadraturePoints(const DistributedProblem<Dimension3>& var,
				    Mesh<Dimension3>& result_mesh, int proc_mode,
				    Vector<int>& IndexBoundary, Vector<int>&,
				    Vector<int>&, Vector<VectR2>& LocalQuadPoints);

    void ConstructProblem(const Mesh<Dimension3>& result_mesh,
			  const DistributedProblem<Dimension3>& var_problem,
			  DistributedProblem<Dimension2>& var_eig, IVect& RotationFaceElement,
			  const string& name_element, R3& vec_u, R3& vec_v) const;

    void ComputeEigenvalues(EigenProblemMontjoie<Real_wp>& eigen_solver,
			    DistributedProblem<Dimension2>& var_eig, VectReal_wp& eigen_values,
			    Matrix<Real_wp, General, ColMajor>& eigen_vectors,
			    int mode_eig = -1, Real_wp coef_shift = 1.0, bool init = true) const;

    void ComputeProjector(const DistributedProblem<Dimension3>& var_problem,
			  const DistributedProblem<Dimension2>& var_eig,
			  const IVect& IndexBoundary, const IVect& IndexElement,
			  const Vector<VectR2>& LocalQuadPoints,
			  Vector<FiniteElementInterpolator>& proj,
			  TinyVector<Vector<Matrix<int> >, 2>& RotationPoints) const;
    
#endif
    
  };
    
  
  //! Class for modal source with dimension-dependent functions
  /*!
    on a boundary, source associated to a Dirichlet condition u = f
    or a Neumann condition du/dn = g, or Robin du/dn + lambda u = f
    f is a mode (sinus or cosinus in 2-D)
   */
  template<class T, class Dimension>
  class ModalSourceBoundary_Dim : public ModalSourceBoundary_Base
  {
  protected :
    typedef typename Dimension::DimensionBoundary DimensionB;
    
    //! coefficients for each mode
    Vector<T> coef_mode_combined;
    
  private :
    const DistributedProblem<Dimension>& var_problem;
    
  public :    
    template<class TypeEquation>
    ModalSourceBoundary_Dim(const EllipticProblem<TypeEquation>& var);

    virtual ~ModalSourceBoundary_Dim();
    
    virtual void ReadSourceParameters(const IVect& boundary_condition, const IVect& ref,
				      const VectString& source_param);
    
    virtual void ComputeMode(int num_mode, const IVect& ref,
			     const VectString& param, Vector<Vector<Vector<T> > >& eval);
    
    virtual void EvaluateMode(const Mesh<Dimension>& mesh, int proc_mode, const Vector<int>& IndexBoundary, 
			      const Vector<int>&, const Vector<int>&,
			      const Vector<typename DimensionB::VectR_N>& LocalQuadPoints,
			      Vector<Vector<Vector<T> > >& eval);
    
  };

} // end namespace

#define MONTJOIE_FILE_MODAL_SOURCE_HXX
#endif
