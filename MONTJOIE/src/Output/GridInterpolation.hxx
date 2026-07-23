#ifndef MONTJOIE_FILE_GRID_INTERPOLATION_HXX

namespace Montjoie
{
  
  //! localization of all points of display grids on the mesh
  /*
    this class regroups all the points you need to localize
    on the mesh. In order to add points, use the method append.
  */
  template<class Dimension>
  class GridInterpolation_Base
  {    
  protected :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_N MatrixN_N;
    
    //! global coordinates
    VectR_N GlobalCoord;
    bool dfjm1_store; //!< DF_i^{-1} has to be computed ?
    VectReal_wp TetaInterp; //!< theta variable (for axisymmetric domain)
    
    //! local coordinates (inside the reference element)
    VectR_N CoorInterp;
    //! values of DF_i^{-1}
    Vector<MatrixN_N> dfjm1;
    
    //! ElementInterp(i) : element number where the point i is
    //! -1 if no element contains the point
    IVect ElementInterp;
    //! distance of each point to the boundary of the element
    //! if negative, the point is slightly outside of the element
    VectReal_wp Epsilon_InsideElement;
    
    // parameters in the case of cyclic domain
    IVect CyclicSection; //!< section number for each point (cyclic domain)
    int nb_sections_cyclic ; //!< number of sections (cyclic domain)
    int nb_sections_x, nb_sections_y, nb_sections_z; //!< for periodic domains
    
    //! number of interpolation points
    int nbPoints_grid;
    Real_wp xmin; //!< minimal x-coordinate
    Real_wp xmax; //!< maximal x-coordinate
    Real_wp ymin; //!< minimal y-coordinate
    Real_wp ymax; //!< maximal y-coordinate
    Real_wp zmin; //!< minimal z-coordinate
    Real_wp zmax; //!< maximal z-coordinate
    
    Vector<bool> point_treated; //!< to avoid relocalize points already localized
    int nb_subdiv_grid_x; //!< number of grid points along x
    int nb_subdiv_grid_y; //!< number of grid points along y
    int nb_subdiv_grid_z; //!< number of grid points along z
    Real_wp step_subdiv_grid_x; //!< space step along x
    Real_wp step_subdiv_grid_y; //!< space step along y
    Real_wp step_subdiv_grid_z; //!< space step along z
    //! ListeBoxGrid_Element(i) is the list of cells intersecting element i
    Vector<IVect> ListeBoxGrid_Element;

#ifdef SELDON_WITH_MPI
    MPI_Comm comm_;
#endif
    
  public :
    
    GridInterpolation_Base();

    size_t GetMemorySize() const;

    int GetNbElt() const;
    int GetNbGlobalPoints() const;
    int GetElementNumber(int i) const;
    void SetElementNumber(int i, int n);
    
    const R_N& GetLocalCoordinate(int i) const;
    const R_N& GetGlobalCoordinate(int i) const;
    const VectR_N& GetGlobalCoordinate() const;
    VectR_N& GetGlobalCoordinate();
    void SetLocalCoordinate(int i, const R_N&);
    void SetGlobalCoordinate(int i, const R_N&);
    
    int GetNbTheta() const;
    const Real_wp& GetTheta(int i) const;
    const VectReal_wp& GetTheta() const;
    VectReal_wp& GetTheta();
    void SetTheta(int i, const Real_wp&);

    bool StoreDFjm1() const;
    const MatrixN_N& GetDFjm1(int i) const;
    void SetDFjm1(int, const MatrixN_N&);
    
    int GetNbSectionNumber() const;
    int GetSectionNumber(int i) const;
    void SetSectionNumber(int, int);
    
    Real_wp GetEpsilonInsideElement(int i) const;
    void SetEpsilonInsideElement(int i, const Real_wp&);
    
    const Real_wp& GetXmin() const;
    const Real_wp& GetXmax() const;
    const Real_wp& GetYmin() const;
    const Real_wp& GetYmax() const;
    const Real_wp& GetZmin() const;
    const Real_wp& GetZmax() const;

    void SetXmin(const Real_wp& x);
    void SetXmax(const Real_wp& x);
    void SetYmin(const Real_wp& y);
    void SetYmax(const Real_wp& y);
    void SetZmin(const Real_wp& z);
    void SetZmax(const Real_wp& z);
    
    void ResizeElements(int n);

    void SetNbSubdivisions(int nx, int ny);
    void SetSubdivisionStep(const Real_wp& dx, const Real_wp& dy);

    void SetNbSubdivisions(int nx, int ny, int nz);
    void SetSubdivisionStep(const Real_wp& dx, const Real_wp& dy, const Real_wp& dz);
    
    void SetNbCyclicSections(int);
    void SetNbCartesianSections(int, int, int);
    
    int GetNbCyclicSections() const;
    void GetNbCartesianSections(int&, int&, int&) const;
    
#ifdef SELDON_WITH_MPI
    void SetCommunicator(const MPI_Comm& comm);
    const MPI_Comm& GetCommunicator() const;
#endif
    
    // appends Points_ND to the list of points
    void Append(const VectR_N& Points_ND, const VectReal_wp&, const IVect&, IVect& number);

    void CompressGrid(Vector<GridInterpolationFull<Dimension> >& grid);
    
    void ClearPrelocalizationArrays();
    
  };
  
  
  //! displaying details of class GridInterpolation_Base
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolation_Base<Dimension>& g);
  
  
  //! class for basic predefined grids
  template<class Dimension>
  class GridInterpolationFull_Base
  {
    template<class Dim>
    friend class GridInterpolation_Base;
    
  protected :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3

    //! pointer to the localization of the points
    GridInterpolation<Dimension>* var_interp;
    
    //! number of the points in object var_interp
    IVect list_points;
    
    //! numbers of global points stored in this processor
    IVect IndexPoints;

    //! total number of grid points
    int nb_global_grid_points;
    
    //! for the master processor, we store point numbers for each processor
    Vector<IVect> list_points_proc;
    
    int nbPoints_x; //!< number of points along x
    int nbPoints_y; //!< number of points along y
    int nbPoints_z; //!< number of points along z
    
    //! boundaries are automatically computed to fit the mesh ?
    bool grid_interval_to_be_computed;
    Real_wp xmin, xmin0; //!< minimal x-coordinate
    Real_wp xmax, xmax0; //!< maximal x-coordinate
    Real_wp ymin, ymin0; //!< minimal y-coordinate
    Real_wp ymax, ymax0; //!< maximal y-coordinate
    Real_wp zmin, zmin0; //!< minimal z-coordinate
    Real_wp zmax, zmax0; //!< maximal z-coordinate
    Real_wp z0_adim; //!< caracteristic length
    
    //! type of grid
    int type_output_file;
    R3 center; //!< center
    //! file where the points may be read
    string PointsInputFile;

  public :
    VectR3 GlobalPoints3D;
    VectR2 GlobalPoints2D;

    GridInterpolationFull_Base();

    size_t GetMemorySize() const;
    
    bool GridIntervalToBeComputed() const;
    
    int GetOutputType() const;
    void SetOutputType(int i);
        
    int GetNbPoints() const;
    int GetNbGlobalPoints() const;
    int GetPointNumber(int) const;
    const IVect& GetPointNumber() const;
    const Vector<IVect>& GetPointNumberAllProc() const;
    
    int GetNbPointsX() const;
    int GetNbPointsY() const;
    int GetNbPointsZ() const;
    void SetNbPointsX(int i);
    void SetNbPointsY(int i);
    
    const R3& GetCenter() const;
    const Real_wp& GetXmin() const;
    const Real_wp& GetXmax() const;
    const Real_wp& GetYmin() const;
    const Real_wp& GetYmax() const;
    const Real_wp& GetZmin() const;
    const Real_wp& GetZmax() const;
    
    void SetXmin(const Real_wp& x);
    void SetXmax(const Real_wp& x);
    void SetYmin(const Real_wp& x);
    void SetYmax(const Real_wp& x);
    void SetZmin(const Real_wp& x);
    void SetZmax(const Real_wp& x);
    
    void SetXmin0(const Real_wp& x);
    void SetXmax0(const Real_wp& x);
    void SetYmin0(const Real_wp& x);
    void SetYmax0(const Real_wp& x);
    void SetZmin0(const Real_wp& x);
    void SetZmax0(const Real_wp& x);

    void AdimSpace(const Real_wp& z);
    
  };
  
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolationFull_Base<Dimension>& g);
  
  
#ifdef MONTJOIE_WITH_ONE_DIM
  //! localization of points for 1-D meshes
  template<>
  class GridInterpolation<Dimension1>
  {
  protected :
    //! global coordinate for each point of the display grid
    VectReal_wp Points;
    //! element number for each point of the display grid
    IVect ElemInterp;
    //! local coordinate for each point of the display grid
    VectReal_wp CoorInterp;
    int nb_global_points;
    
  public :
    VectReal_wp TetaInterp, PhiInterp;
    Vector<IVect> ProcNumberPoints; //!< global point numbers for each proc
    
    GridInterpolation();
    
    size_t GetMemorySize() const;
    
    const Real_wp& GetXmin();
    const Real_wp& GetXmax();
    
    void Init(const Real_wp& xmin, const Real_wp& xmax, int N);
    void Init(const VectReal_wp& coord);
    void AddPoint(const Real_wp& x);
    
    int GetNbPointsGrid() const;
    int GetElementNumber(int i) const;
    Real_wp GetLocalCoordinate(int i) const;
    Real_wp GetGlobalCoordinate(int i) const;
    const VectReal_wp& GetGlobalCoordinate() const;
    
    void ReallocatePoints(int);
    void SetElementNumber(int i, int);
    void SetLocalCoordinate(int i, const Real_wp&);
    void SetGlobalCoordinate(int i, const Real_wp&);

    void Append(const VectReal_wp& Points_ND, const VectReal_wp& teta,
                const VectReal_wp& phi, IVect& number);    
    
    void LocalizePoints(Mesh<Dimension1>& mesh);    
    
    int GetGlobalNumberPoints() const;
    void SetGlobalNumberPoints(int);
    
  };
  
  
  //! basic grids in 1-D
  /*
    it deals with outputs on 1-D regular grids, lines or points
  */
  template<>
  class GridInterpolationFull<Dimension1>
    : public GridInterpolationFull_Base<Dimension1>
  {
  public :
    enum {LINE, POINT, PLANE_AXI, LINE_AXI, POINT_AXI, CIRCLE_AXI,
          THREE_PLANES_AXI, VOLUME_AXI, POINTS_FILE_AXI, PLANE};
    
    GridInterpolationFull();
    
    void SetInputData(const string& description_field,
		      const VectString& parameters);
    
    static int GetType(const string& description);
    static int GetDimension(int type);

    void Init(int type, const Real_wp& x0, const Real_wp& xN, int N);
    
    void SetPlaneAxiOutput(const R3&, const R3&, const R3&, int, int);
    void SetLineAxiOutput(const R3& ptA, const R3& ptB, int nbx);
    void SetPointAxiOutput(const R3& ptA);
    
    void SetCircleAxiOutput(const R3& ptA, const R3& normale, 
			    const Real_wp&, const Real_wp&, int nbx);

    void SetPointsFileAxiOutput(const string&);
    
    void SetThreePlanesAxiOutput(const R3&, const R3&, const R3&,
				 int nbx_, int nby_, int nbz_);

    void SetVolumeAxiOutput(const R3&, const R3&,
			    int nbx_, int nby_, int nbz_);

    void GenerateGridPoints(VectReal_wp& Points1D, VectReal_wp& TetaPoints);    
    void GenerateGridPoints(VectReal_wp& Points1D, VectReal_wp& TetaPoints, VectReal_wp& PhiPoints);    

    void InitGrid(GridInterpolation<Dimension1>& grid, int dim = 3);

    void ReadText(istream& FileStream);
    
    void Read(istream& FileStream, bool double_prec = false);
    
    void ReadVtk(istream& FileStream, bool&, int&);    
    
    void Write(ostream& FileStream, bool double_prec = false) const;
    
    template<class real>
    void WriteBinary(ostream& FileStream, real& a) const;
    
    void WriteText(ostream& FileStream, int double_prec = 0) const;

    void WriteVtk(ostream& FileStream, bool double_prec = false, bool ascii = false) const;

    void WriteGrid(const string& file_name);
    
  };
  
  ostream& operator<<(ostream& out, const GridInterpolationFull<Dimension1>& g);

#endif
  
  
#ifdef MONTJOIE_WITH_TWO_DIM
  //! localization of points for 2-D meshes
  template<>
  class GridInterpolation<Dimension2> : public GridInterpolation_Base<Dimension2>
  {
  public :
    GridInterpolation();
    
    void InitInterpolationGrid(const Mesh<Dimension2>& mesh);
    
    void LocalizePoints(const Mesh<Dimension2>& mesh);
    void PreLocalizePoints(IVect& NumBoxGrid_Point, Vector<IVect>& ListPoints_Grid,
			   IVect& NbPoints_Grid, bool compute_list);
    
    void LocalizePoints(const Mesh<Dimension2>& mesh, const VectR2& Points2D,
			IVect& ElementInt, VectR2& CoorInt, VectReal_wp& EpsilonInt,
			Vector<bool>& point_treat, Vector<Matrix2_2>& mat_dfjm1,
			bool store_df = false);

  };
  
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolation<Dimension2>& g);
  
  
  //! class for predefined grids in 2-D
  /*
    it deals with outputs on 2-D regular grids, lines, points or circles
   */
  template<>
  class GridInterpolationFull<Dimension2>
    : public GridInterpolationFull_Base<Dimension2>
  {
  public :
    
    enum {PLANE, PLANE_AXI, LINE, LINE_AXI, POINT, POINT_AXI, CIRCLE,
	  CIRCLE_AXI, THREE_PLANES_AXI, VOLUME_AXI, POINTS_FILE, POINTS_FILE_AXI}; 

  protected :
    template<class real>
    void WriteBinary(ostream& FileStream, real& a) const;

    template<class real>
    void ReadBinary(istream& FileStream, real& a);
    
  public :
    GridInterpolationFull();
    
    void SetInputData(const string& description_field,
		      const VectString& parameters);
    
    static int GetType(const string& description);
    static int GetDimension(int type);
    
    void SetPlaneOutput(const Real_wp& xmin_, const Real_wp& xmax_,
			const Real_wp& ymin_, const Real_wp& ymax_,
			int nbx_, int nby_);

    void SetPlaneAxiOutput(const R3&, const R3&, const R3&, int, int);

    void SetLineOutput(const Real_wp& xmin_, const Real_wp& xmax_,
		       const Real_wp& ymin_, const Real_wp& ymax_, int nbx_);

    void SetLineAxiOutput(const R3& ptA, const R3& ptB, int nbx);

    void SetPointOutput(const R2& ptA);
    void SetPointsOutput(const VectR2& ptA);
    
    void SetPointAxiOutput(const R3& ptA);
    
    void SetCircleOutput(const R2& center, const Real_wp& radius, int nbx);
    void SetCircleAxiOutput(const R3& ptA, const R3& normale, 
			    const Real_wp&, const Real_wp&, int nbx);
    
    void SetPointsFileOutput(const string&);
    void SetPointsFileAxiOutput(const string&);
    
    void SetThreePlanesAxiOutput(const R3&, const R3&, const R3&,
				 int nbx_, int nby_, int nbz_);

    void SetVolumeAxiOutput(const R3&, const R3&,
			    int nbx_, int nby_, int nbz_);
    
    void GenerateGridPoints(const GridInterpolation<Dimension2>& grid, IVect& Index,
			    VectR2& Points2D, VectReal_wp& TetaPoints, IVect& NumPartie);

    void GenerateGridPoints(VectR2& Points2D, VectReal_wp& TetaPoints);
    
    bool TranslatePoint(const R3& point, bool axi, const GridInterpolation<Dimension2>& grid,
			R2& pt2D, Real_wp& theta, int& n);
        
    void InitGrid(GridInterpolation<Dimension2>& grid);
    
    // writing header in binary
    void Write(ostream& FileStream, bool double_prec = false) const;
    
    // writing header in ascii
    void WriteText(ostream& FileStream, int double_prec = 0) const;
    
    void WriteVtk(ostream& FileStream, bool double_prec = false,
		  bool ascii = false) const;

    void Read(istream& FileStream, bool double_prec = false);
    
    void ReadText(istream& FileStream);

    void ReadVtk(istream& FileStream, bool&, int&);
    
    void WriteGrid(const string& file_name);
  };
  
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolationFull<Dimension2>& g);
  
#endif
#if defined(MONTJOIE_WITH_THREE_DIM)
  //! localization of points for 3-D meshes
  template <>
  class GridInterpolation<Dimension3> : public GridInterpolation_Base<Dimension3>
  {
  public :

    GridInterpolation();
      
    void InitInterpolationGrid(const Mesh<Dimension3>& mesh);
    void PreLocalizePoints(IVect& NumBoxGrid_Point, Vector<IVect>& ListPoints_Grid,
			   IVect& NbPoints_Grid, bool compute_list);

    void LocalizePoints(const Mesh<Dimension3>& mesh);
    
    void LocalizePoints(const Mesh<Dimension3>& mesh, const VectR3& Points2D,
			IVect& ElementInt, VectR3& CoorInt, VectReal_wp& EpsilonInt,
			Vector<bool>& point_treat, Vector<Matrix3_3>& mat_dfjm1,
			bool store_df = false);
    
  };
  
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolation<Dimension3>& g);
  
  
  //! basic for predefined grids in 3-D
  /*
    it deals with outputs on 3-D regular grids, planes, lines, or points
   */
  template<>
  class GridInterpolationFull<Dimension3>
    : public GridInterpolationFull_Base<Dimension3>
  {
  public :
    enum {VOLUME, THREE_PLANES, PLANE, LINE,
	  POINT, POINTS_FILE, CIRCLE};

  protected :
    template<class real>
    void WriteBinary(ostream& FileStream, real& a) const;

    template<class real>
    void ReadBinary(istream& FileStream, real& a);
    
  public :    
    GridInterpolationFull();
    
    void SetInputData(const string& description_field,
		      const VectString& parameters);
    
    static int GetType(const string& description);
    static int GetDimension(int type);
    
    void SetPlaneOutput(const R3& ptO, const R3& ptA, const R3& ptB,
			int nbx, int nby);
    
    void SetLineOutput(const R3& ptA, const R3& ptB, int nbx);
    void SetPointOutput(const R3& ptA);
    
    void SetCircleOutput(const R3& ptA, const R3& normale,
			 const Real_wp& rx, const Real_wp& ry, int nbx);
    
    void SetPointsFileOutput(const string& name);
    void SetPointsOutput(const VectR3& ptA);
    
    void SetThreePlanesOutput(const R3& ptO, const R3& ptMin, const R3& ptMax,
			      int nbx_, int nby_, int nbz_);

    void SetVolumeOutput(const R3& ptMin, const R3& ptMax,
			 int nbx_, int nby_, int nbz_);
  
    void GenerateGridPoints(const GridInterpolation<Dimension3>&, IVect&,
			    VectR3& Points3D, IVect&);

    void GenerateGridPoints(VectR3& Points3D);
    
    void InitGrid(GridInterpolation<Dimension3>& grid);

    bool TranslatePoint(const R3& point, const GridInterpolation<Dimension3>& grid,
			R3& pt3D, int& num_partie);
    
    // writing header in binary
    void Write(ostream& FileStream, bool double_prec = false) const;
    
    // writing header in ascii
    void WriteText(ostream& FileStream, int double_prec = 0) const;
    
    void WriteVtk(ostream& FileStream, bool double_prec = false,
		  bool ascii = false) const;
    
    void Read(istream& FileStream, bool double_prec = false);
    
    void ReadText(istream& FileStream);

    void ReadVtk(istream& FileStream, bool&, int&);
    
    void WriteGrid(const string& file_name);
  };
  
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolationFull<Dimension3>& g);

#endif

  template<class T, class Dimension>
  void WriteVtk(const Vector<T>& val, const GridInterpolationFull<Dimension>& grid,
		const string& file_name, bool double_prec = false, bool ascii = false);

  template<class T, class Dimension>
  void ReadVtk(Vector<T>& val, GridInterpolationFull<Dimension>& grid,
	       const string& file_name, bool double_prec, bool ascii);

  void WriteBinaryTypeData_DoubleOrFloat(const VectReal_wp& output_vector,
					 ostream& file_out, bool double_prec);
  
  void WriteBinaryTypeData_DoubleOrFloat(const VectComplex_wp& output_vector,
					 ostream& file_out, bool double_prec);

  void WriteTypeData_DoubleOrFloat(const VectReal_wp& output_vector,
				   ostream& file_out, int type_prec);

  void WriteTypeData_DoubleOrFloat(const VectComplex_wp& output_vector,
				   ostream& file_out, int type_prec);

  template<class T, class Dimension>
  void WriteMatlab(const Vector<T>& val,
		   GridInterpolationFull<Dimension>& grid,
		   const string& file_name,
		   int double_prec = OutputTypeEnum::SINGLE_PRECISION, bool ascii = false);

  template<class T, class Dimension>
  void WriteMatlab(const Vector<T>& val, int nb_true,
		   const Vector<GridInterpolationFull<Dimension> >& grid, int type,
		   const string& file_name,
		   int double_prec = OutputTypeEnum::SINGLE_PRECISION, bool ascii = false);

  template<class T, class Dimension>
  void ReadMatlab(Vector<T>& val,
		  GridInterpolationFull<Dimension>& grid,
		  const string& file_name, bool ascii);
  
  template<class T, class Dimension>
  void ReadMatlab(Vector<T>& output_vector,
		  Vector<GridInterpolationFull<Dimension> >& grid,
		  const string& file_name, bool ascii);
  
} // namespace Montjoie

#define MONTJOIE_FILE_GRID_INTERPOLATION_HXX
#endif

