#ifndef MONTJOIE_FILE_OUTPUT_HARMONIC_HXX

namespace Montjoie
{
  //! base class for outputs
  class VarOutputProblem_Base : public InputDataProblem_Base
  {
  public:
    Vector<ParamOutputClass> output_grid_param; //!< output parameters
    Vector<ParamOutputClass> output_mesh_param; //!< output parameters
    ParamOutputClass output_points_outside_param; //! parameters for points outside the domain    
    bool grid_to_be_computed;
    
  protected:
    bool display_solution; //!< solution ready to be displayed ?
    bool keep_grid_localization_array;
    
    //! component number to be displayed on grids (-1 -> all components)
    int choice_field_output;

    //! format of the output files (binary, ascii)
    int format_output_file;
    //! extension to use for results obtained on predefined grids (.dat or .vtk)
    string output_file_extension;
    //! precision for output
    int precision_output_file;
    
    //! subdivided mesh relies on Gauss-Lobatto points ?
    bool output_lobatto_points; 
    
    //! parameters for writing directly the output (without interpolation)
    ParamOutputClass output_solution_param;
    //! if true, the solution is directly written in .elb (without interpolation)
    bool output_solution_file;
    bool write_solution_on_quad_points;
    string file_name_solution_quad, file_name_weights_quad;
    
    //! if true the points are moved (for an output on a surface mesh)
    //! such that we can see the deformation of a plate for example
    bool move_points_on_output_mesh;
    //! which component of the solution is used to move the points ?
    int component_used_to_translate;
    //! coefficient used to move the points
    Real_wp coefficient_used_to_translate;
    
  public :
    //! directory where output files are written
    string DOSSIER_output;
    //! if true, the output files must be written, otherwise they can be bufferized
    bool last_output_solution;
    
  private:
    void InitDefaultValues();

    VarProblem_Base& var_problem;
    DistributedProblem_Base& var_comm;
    VarBoundaryCondition_Base& var_boundary;
    VarSourceProblem_Base& var_source;

  public :

    template<class TypeEquation>
    VarOutputProblem_Base(EllipticProblem<TypeEquation>&);
    
    // Inline methods
    bool OutputWrittenInDoublePrecision();
    void InitOutput(const Real_wp& t0);
    void SetOutputComponent(int t);
    int GetOutputComponent() const;
    void KeepGridLocalizationArrays();

  protected:
    template<class T>
    static void AddContribInterpol(const T& coef, const TinyVector<Real_wp, 1>& phi,
                                   int p, Vector<T>& val);
    
    template<class T>
    static void AddContribInterpol(const T& coef, const R2& phi, int p, Vector<T>& val);

    template<class T>
    static void AddContribInterpol(const T& coef, const R3& phi, int p, Vector<T>& val);

  public:
    // other methods
    int GetIndexOutputFiles(int type);

    template<class T>
    void WriteDatas(Vector<T> & U0);
    
    virtual void WriteDatas(Vector<VectReal_wp> & U0);
    virtual void WriteDatas(Vector<VectComplex_wp> & U0);

    template<class T>
    void GetModalOutput(int n, int nm, FftInterface<Complex_wp>& fft_interface,
			Vector<T>& val_v, Vector<T>& val_grad_v, 
			int off_u, Vector<T>& val_u, int off_v, Vector<T>& val_grad_u,
                        bool compute_grad) const;    
    
    void SetInputData(const string& description_field, const VectString& parameters);

    template<class T>
    void WriteOutputFileReshaped(const Vector<Vector<T> >& U0, const string& name_file) const;

    template<class T>
    void ReshapeVector(const Vector<T>& U, Vector<Vector<T> >& Uvec, int nb_u = -1) const;
    
    template<class T>
    void WriteOutputFile(const Vector<T>& U0, int num_file, int nb_u = -1) const;
    
    template<class T>
    void WriteOutputFile(const Vector<T>& U0, int num_file, int nb_u, const Vector<int>&) const;

    template<class T>
    void WriteOutputFile(const Vector<T>& U0, const string& name_file, int nb_u, const Vector<int>&) const;
    
    template<class T>
    void WriteOutputFile(const Vector<T>& U0, const string& name_file, int nb_u = -1) const;

    void ChangeTimeSnapshot(int, const Real_wp&);
    
    template<class T>
    void WriteOutputFile(const Vector<Vector<T> > & U0, int num_file, int nb_uloc = -1) const;

    void CopyInputData(const VarOutputProblem_Base& var);    

    inline void ApplyInversePoint(FftInterface<Complex_wp>& fft, int ix, const VectReal_wp& u, Real_wp& v) const {}

    inline void ApplyInversePoint(FftInterface<Complex_wp>& fft, int ix, int iy, const VectReal_wp& u, Real_wp& v) const {}
    
    inline void ApplyInversePoint(FftInterface<Complex_wp>& fft, int ix, int iy, int iz, const VectReal_wp& u, Real_wp& v) const {}

    inline void ApplyInversePoint(FftInterface<Complex_wp>& fft, int ix, const VectComplex_wp& u, Complex_wp& v) const { fft.ApplyInversePoint(ix, u, v); }

    inline void ApplyInversePoint(FftInterface<Complex_wp>& fft, int ix, int iy, const VectComplex_wp& u, Complex_wp& v) const { fft.ApplyInversePoint(ix, iy, u, v); }
    
    inline void ApplyInversePoint(FftInterface<Complex_wp>& fft, int ix, int iy, int iz, const VectComplex_wp& u, Complex_wp& v) const { fft.ApplyInversePoint(ix, iy, iz, u, v); }

    // Vtk outputs
    virtual void GetVtkVolumeOutput(int nb_u, int nb_du, int dim_N, bool cplx,
				    IVect&, Vector<IVect>&, Vector<string>&,
				    IVect&, Vector<IVect>&, Vector<string>&) const;
    
    virtual void GetVtkSurfaceOutput(int nb_u, int nb_du, int dim_N, bool cplx,
				     IVect&, Vector<IVect>&, Vector<string>&,
				     IVect&, Vector<IVect>&, Vector<string>&) const;

    virtual void WriteOutputFileReshaped(const Vector<VectReal_wp>& U0, int num_file) const = 0;
    virtual void WriteOutputFileReshaped(const Vector<VectComplex_wp>& U0, int num_file) const = 0;    

    virtual void WriteSnapshot(int& nb_iter, const Real_wp& t, VectReal_wp& Uh, int num_file) const = 0;
    virtual void WriteSnapshot(int& nb_iter, const Real_wp& t, VectComplex_wp& Uh, int num_file) const = 0;

    virtual void WriteUquadrature(const Vector<VectReal_wp>& U0,
                                  const string& file_name_u, const string& file_name_poids) const = 0;

    virtual void WriteUquadrature(const Vector<VectComplex_wp>& U0,
                                  const string& file_name_u, const string& file_name_poids) const = 0;
    
    virtual void InitVarGrid() = 0;
    virtual void ComputeVarGrid() = 0;
    
  };
  
  
  //! base class for outputs of finite element solutions
  template<class Dimension>
  class VarOutputProblem_Dim : public VarOutputProblem_Base
  {
  public:
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    typedef typename Dimension::MatrixN_N MatrixN_N;

    //! localization of display points in the mesh
    GridInterpolation<Dimension> all_points_display;
    //! predefined display grids
    Vector<GridInterpolationFull<Dimension> > var_grid;
    //! display can be done on meshes
    Vector<MeshInterpolation<Dimension> > var_mesh_interp;
    
  private :
    VarBoundaryCondition_Base& var_boundary;
    VarSourceProblem_Dim<Dimension>& var_source;
    VarSourceProblem_Base& var_source_base;
    VarProblem<Dimension>& var_problem;    
    VarComputationRCS_Base<Dimension>& output_rcs_param;
    
  public :

    template<class TypeEquation>
    VarOutputProblem_Dim(EllipticProblem<TypeEquation>&);

  protected:
    // internal methods
    void SetNbCyclicModes();
    
  public :
    // public methods
    void GetMemoryUsed(map<string, size_t>& var) const;    
    void SetInputData(const string& description_field, const VectString& parameters);
    void InitVarGrid();
    void ComputeVarGrid();

    void InitVarMeshVolumetric(MeshInterpolation<Dimension>& var_interp);
    void InitVarMeshSurfacic(MeshInterpolation<Dimension>& var_interp);
    
    template<class T>
    void ComputeInterpolationU_GradU(const Vector<Vector<T> > & U0,
				     const GridInterpolation<Dimension> & var_interp,
				     const GridInterpolationFull<Dimension>& var_gr,
				     Vector<Vector<T> >& trace_u, Vector<Vector<T> >& trace_grad_u,
				     int add_total_field, const Real_wp& t,
				     Vector<Vector<T> >& trace_u_diff, Vector<Vector<T> >& trace_grad_diff,
				     bool compute_grad) const;

    template<class T>
    void ComputeInterpolationU_GradU(const Vector<T> & U0, int num_grid, bool compute_grad,
				     Vector<Vector<T> >& trace_u, Vector<Vector<T> >& trace_grad_u,
                                     int nb_u = 1) const;
    
    template<class T>
    void ComputeInterpolationNodalU(const Vector<Vector<Vector<T> > > & Unodal,
				    const Vector<Vector<Vector<T> > >& GradNodal,
				    const GridInterpolationFull<Dimension>& var_gr,
				    int nb_comp_all, Vector<Vector<T> >& trace_u) const;
    
    template<class T>
    void ComputeInterpolationNodalUloc(const Vector<Vector<Vector<T> > > & Unodal,
				       const Vector<Vector<Vector<T> > >& GradNodal,
				       const GridInterpolation<Dimension> & var_interp,
				       Vector<T>& trace_vec, int nb_comp_all,
				       const IVect& list_points, int nnz) const;
    
    void CopyInputData(const VarOutputProblem_Base& var);

    template<class T>
    void WriteOutputFileGen(const Vector<Vector<T> >& U0,
			    const Vector<MeshInterpolation<Dimension> >& var_interp,
			    const Vector<ParamOutputClass>& grid, int num_file, int nb_u = -1) const;
    
    template<class T, class TypeInterpolation>
    void WriteOutputFileGen(const Vector<Vector<T> > & U0, const TypeInterpolation & var_interp,
			    int num_file, const string & file_name, 
			    int add_total_field, const string& file_name_diff,
			    const ParamOutputClass& param, int append) const;
    
    VarComputationRCS_Base<Dimension>& GetParameterOutputRCS();

    template<class T>
    void WriteOutputFileGen(const Vector<Vector<T> >& U0,
			    const VarComputationRCS_Base<Dimension>& var_rcs,
			    const ParamOutputClass& param, int num_file) const;
    
    virtual void SetComputationFarPoints(VectR_N& Points, Real_wp dt) = 0;

    template<class T, class TypeInterpolation>
    void WriteOutputFileNodalGen(const Vector<Vector<T> > & U0, const TypeInterpolation & var_interp,
				 const Vector<ParamOutputClass>& grid, int num_file) const;
    
    template<class T, class TypeInterpolation>
    void WriteOutputFileGen(const Vector<Vector<T> > & U0, const TypeInterpolation & var_interp,
			    const Vector<ParamOutputClass>& grid, int num_file) const;

    void WriteOutputFileReshaped(const Vector<VectReal_wp>& U0, int num_file) const;
    void WriteOutputFileReshaped(const Vector<VectComplex_wp>& U0, int num_file) const;
    
    template<class T, class VarInterpolation>
    void WriteOutputVector(const VarInterpolation & var_interp, Vector<T> & output_vector,
			   Vector<Vector<T> >& output_u, Vector<Vector<T> >& output_du,
			   int num_file, const string & file_name,
			   int type, int m, const ParamOutputClass& param, int append) const;
    
    template<class T>
    void WriteSnapshot(int nb_iter, const Real_wp& t, Vector<T>& Uh,
		       const Vector<MeshInterpolation<Dimension> >& var_interp,
		       bool increment) const;
    
    template<class T>
    void WriteSnapshot(int nb_iter, const Real_wp& t, Vector<T>& Uh,
		       const ParamOutputClass& output_param,
		       const Vector<GridInterpolationFull<Dimension> >& var_interp,
		       bool increment) const;
    
    virtual void ComputeInterpolationUloc(const Vector<Vector<Real_wp> > & U0,
					  const GridInterpolation<Dimension> & var_interp,
					  Vector<Real_wp>& trace_vec, Vector<Real_wp>& trace_grad_vec,
					  const IVect& list_points, int nnz, bool compute_grad) const = 0;
    
    virtual void ComputeInterpolationUloc(const Vector<Vector<Complex_wp> > & U0,
					  const GridInterpolation<Dimension> & var_interp,
					  Vector<Complex_wp>& trace_vec, Vector<Complex_wp>& trace_grad_vec,
					  const IVect& list_points, int nnz, bool compute_grad) const = 0;
    
  public:
    void WriteSnapshot(int& nb_iter, const Real_wp& t, VectReal_wp& Uh, int num_file) const;
    void WriteSnapshot(int& nb_iter, const Real_wp& t, VectComplex_wp& Uh, int num_file) const;
    
    virtual void ComputeNodalUgradU(const Vector<VectReal_wp>& U0,
				    Vector<Vector<VectReal_wp> >& Unodal,
				    Vector<Vector<VectReal_wp> >& GradNodal,
				    bool u_component, bool grad_component,
				    bool hess_component = false,
				    bool true_unknown = false, int nb_u = -1) const = 0;
    
    virtual void ComputeNodalUgradU(const Vector<VectComplex_wp>& U0,
				    Vector<Vector<VectComplex_wp> >& Unodal,
				    Vector<Vector<VectComplex_wp> >& GradNodal,
				    bool u_component, bool grad_component,
				    bool hess_component = false,
				    bool true_unknown = false, int nb_u = -1) const = 0;

    virtual void ComputeQuadratureUgradU(const Vector<VectReal_wp>& U0,
					 Vector<Vector<VectReal_wp > >& Uquad,
					 Vector<Vector<VectReal_wp> >& GradQuad,
					 bool u_component, bool grad_component,
					 bool true_unknown, int nb_u = -1) const = 0;
    
    virtual void ComputeQuadratureUgradU(const Vector<VectComplex_wp>& U0,
					 Vector<Vector<VectComplex_wp> >& Uquad,
					 Vector<Vector<VectComplex_wp> >& GradQuad,
					 bool u_component, bool grad_component,
					 bool true_unknown, int nb_u = -1) const = 0;

    template<class T>
    void WriteUquadratureGen(const Vector<Vector<T> >& U0,
                             const string& file_name_u, const string& file_name_poids) const;
    template<class T>
    void ComputeValueNodalBoundary(const Vector<T>& Uloc, Vector<T>& Uboundary, int num_elem, int num_loc) const;
    
    virtual void WriteUquadrature(const Vector<VectReal_wp>& U0,
                                  const string& file_name_u, const string& file_name_poids) const;

    virtual void WriteUquadrature(const Vector<VectComplex_wp>& U0,
                                  const string& file_name_u, const string& file_name_poids) const;
    
  };


  //! base class implementing some methods for the outputs
  template<class Dimension>
  class VarOutputProblem : public VarOutputProblem_Dim<Dimension>
  {

  private:
    VarProblem<Dimension>& var_problem;
    VarSourceProblem_Fem<Dimension>& var_source;
    VarBoundaryCondition_Base& var_boundary;
    
  protected :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    
    template<class T>
    void ComputeInterpolationUlocGen(const Vector<Vector<T> > & U0,
				     const GridInterpolation<Dimension> & var_interp,
				     Vector<T>& trace_vec, Vector<T>& trace_grad_vec,
				     const IVect& list_points, int nnz, bool compute_grad) const;

    // Inline methods
    virtual void ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u, int i,
				     const GridInterpolation<Dimension>& var_interp,
				     int iquad, bool compute_grad) const;
    
    virtual void ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u, int i,
				     const GridInterpolation<Dimension>& var_interp,
				     int iquad, bool compute_grad) const;

    virtual void ModifyOutputUnknown(Vector<VectReal_wp>&, Vector<VectReal_wp>&,
				     int, bool, bool) const;
    
    virtual void ModifyOutputUnknown(Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
				     int, bool, bool) const;

    virtual void ComputeInterpolationUloc(const Vector<Vector<Real_wp> > & U0,
					  const GridInterpolation<Dimension> & var_interp,
					  Vector<Real_wp>& trace_vec, Vector<Real_wp>& trace_grad_vec,
					  const IVect& list_points, int nnz, bool compute_grad) const;
    
    virtual void ComputeInterpolationUloc(const Vector<Vector<Complex_wp> > & U0,
					  const GridInterpolation<Dimension> & var_interp,
					  Vector<Complex_wp>& trace_vec, Vector<Complex_wp>& trace_grad_vec,
					  const IVect& list_points, int nnz, bool compute_grad) const;
  public :
    template<class TypeEquation>
    VarOutputProblem(EllipticProblem<TypeEquation>&);
    
    void ComputeNodalUgradU(const Vector<VectReal_wp>& U0,
			    Vector<Vector<VectReal_wp> >& Unodal,
			    Vector<Vector<VectReal_wp> >& GradNodal,
			    bool u_component, bool grad_component,
			    bool hess_component = false,
			    bool true_unknown = false, int nb_u = -1) const;

    void ComputeNodalUgradU(const Vector<VectComplex_wp>& U0,
			    Vector<Vector<VectComplex_wp> >& Unodal,
			    Vector<Vector<VectComplex_wp> >& GradNodal,
			    bool u_component, bool grad_component,
			    bool hess_component = false,
			    bool true_unknown = false, int nb_u = -1) const;

    void ComputeQuadratureUgradU(const Vector<VectReal_wp>& U0,
				 Vector<Vector<VectReal_wp > >& Uquad,
				 Vector<Vector<VectReal_wp> >& GradQuad,
				 bool u_component, bool grad_component,
				 bool true_unknown, int nb_u = -1) const;

    void ComputeQuadratureUgradU(const Vector<VectComplex_wp>& U0,
				 Vector<Vector<VectComplex_wp> >& Uquad,
				 Vector<Vector<VectComplex_wp> >& GradQuad,
				 bool u_component, bool grad_component,
				 bool true_unknown, int nb_u = -1) const;
    
    // other methods
  protected:
    template<class T>
    void ComputeNodalUgradU_Gen(const Vector<Vector<T> >& U0,
				Vector<Vector<Vector<T> > >& Unodal,
				Vector<Vector<Vector<T> > >& GradNodal,
				bool u_component, bool grad_component,
				bool hess_component = false,
				bool true_unknown = false, int nb_u = -1) const;
    
    template<class T>
    void ComputeQuadratureUgradU_Gen(const Vector<Vector<T> >& U0,
				     Vector<Vector<Vector<T> > >& Uquad,
				     Vector<Vector<Vector<T> > >& GradQuad,
				     bool u_component, bool grad_component,
				     bool true_unknown, int nb_u = -1) const;
    
  };
    
} // namespace Montjoie

#define MONTJOIE_FILE_OUTPUT_HARMONIC_HXX
#endif
