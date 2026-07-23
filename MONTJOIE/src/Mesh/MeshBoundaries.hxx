#ifndef MONTJOIE_FILE_MESH_BOUNDARIES_HXX

namespace Montjoie
{
  
  template<class Dimension, class T = Real_wp>
  class PmlRegionParameter
  {
  };


  //! parameters specifying parameters of the PML region to add
  template<class Dimension, class T = Real_wp>
  class PmlRegionParameter_Base
  {
  protected:
    // PML layers to add ?
    int pml_add_axis_x; //!< pml layers to add in x-direction ?
    int pml_add_axis_y; //!< pml layers to add in y-direction ?
    int pml_add_axis_z; //!< pml layers to add in z-direction ?  
    T thicknessPML; //!< PML thickness
    T radiusPML; //!< radius of boundary before adding PML layers
    TinyVector<T, Dimension::dim_N> originPML;
    int nb_layers_pml; //! number of layers in PML domain
    Vector<int> ref_pml;
    
  public:
    // PML layers
    enum {PML_NO, PML_BOTH_SIDES, PML_POSITIVE_SIDE, PML_NEGATIVE_SIDE, PML_RADIAL};

    PmlRegionParameter_Base();
    
    // inline methods
    const T& GetThicknessPML() const;
    const TinyVector<T, Dimension::dim_N>& GetOriginRadialPML() const;
    const T& GetRadiusPML() const;
    void SetThicknessPML(const T&);
    int GetTypeAdditionPML(int) const;    
    int GetNbLayersPML() const;
    void SetRadiusPML(const T&);
    void SetOriginRadialPML(const TinyVector<T, Dimension::dim_N>&);
    
    PmlRegionParameter<Dimension, T>& GetLeafClass();
    
    // other methods
    void SetAdditionPML(int, int, int, int);
    void SetParameters(const Vector<string>& parameters);
    
    void AddPML(int, Mesh<Dimension>&);
  };


  //! base class for curved boundaries
  template<class Dimension, class T = Real_wp>
  class MeshBoundaries_Base
  {
    template<class DimB, class T0>
    friend class Mesh_Base;
    
  protected :
    typedef TinyVector<T, Dimension::dim_N> R_N;
    typedef Vector<R_N> VectR_N;
    
    
    //! list of PML areas to add
    Vector<PmlRegionParameter<Dimension, T> > pml_areas;
    
    /**********
     * Curves *
     **********/
    
    //! type of curve for each referenced boundary (circle, sphere, etc)
    IVect Type_curve;
    //! boundary condition for each referenced boundary (dirichlet, neumann, etc)
    IVect Cond_curve;
    //! owning body for each referenced boundary (used for transparent condition)
    IVect Body_curve;
    VectBool Param_known; //!< true if parameters of curve are known    
    
    // Param_curve(i) : parameters for curve i
    // for example center and radius in case of a circle
    // Param_condition(i) : parameters for curve i
    // (parameters related to a boundary conditions)
    //! parameters of the curve (radius and center of circle by example)
    Vector<Vector<T> > Param_curve;
  
    //! interpolation on curved edges
    Globatto<Real_wp> lob_curve;

    
    /***********************
     * Boundary conditions *
     ***********************/
    
    
    //! list of references with a periodic condition
    Vector<TinyVector<int, 2> > periodicity_references;
    
    //! type of periodicity
    IVect type_coord_periodicity;
    R_N periodic_center_polar; Real_wp periodic_alpha_polar;
    
    Mesh<Dimension, T>& GetLeafClass();
    const Mesh<Dimension, T>& GetLeafClass() const;
    
  public :
    
    MeshBoundaries_Base();
    
    /* Basic methods */

    size_t GetMemorySize() const;    
    int GetNbReferences() const;
    int GetCurveType(int) const;
    const IVect& GetCurveType() const;
    int GetBoundaryCondition(int) const; 
    const IVect& GetBoundaryCondition() const;
    int GetBodyNumber(int) const;
    const IVect& GetBodyNumber() const;
    void GetCurveParameter(int, Vector<T>&) const;
    const Vector<T>& GetCurveParameter(int) const;
    const Vector<Vector<T> >& GetCurveParameter() const;
    const Globatto<Real_wp>& GetCurveFunctions1D() const;
    int GetNbPmlAreas() const;
    void ReallocatePmlAreas(int n);
    PmlRegionParameter<Dimension, T>& GetPmlArea(int i);
    const PmlRegionParameter<Dimension, T>& GetPmlArea(int i) const;
    
    template<class Dim>
    void SetBoundaryCondition(const Mesh<Dim,T>&); 
    void SetBoundaryCondition(int, int); 
    void SetCurveType(int, int); 
    void SetBodyNumber(int, int);
    void SetCurveParameter(int, const Vector<T>&);
    void CopyCurves(const Mesh<Dimension>&);
    void CopyInputData(const Mesh<Dimension>&);
    
    /* Convenient methods */
    
    // return an unused reference number
    int GetNewReference();
    void GetNewReference(int&, int&);
    // resizing arrays in order to deal with a curve referenced ref_max
    void ResizeNbReferences(int ref_max);
    // changing type of curve (straight, circle, spline...)
    // and associated boundary condition of a reference
    void SetNewReference(int ref, int type, int cond);
    
    /* Periodicity */
    
    int GetNbPeriodicReferences() const;
    int GetPeriodicReference(int, int) const;
    void AddPeriodicCondition(const TinyVector<int, 2>&, int);
    
    int GetPeriodicityTypeReference(int i) const;
    
    Real_wp GetPeriodicAlpha() const;
    void SetPeriodicAlpha(const Real_wp& alpha);
    
    void ClearPeriodicCondition();
    
    /* Methods for curved boundaries */
    
    // returns safety coefficient for a curved element, returns 0 for a straight element
    T GetSafetyCoef(int i) const;
    void SetInputData(const string& description_field, const VectString& parameters);
    
    bool AreCurveEqual(int ref1, int ref2) const;
    int GetGeometryOrder() const;
    Real_wp GetInterpolate1D(int i, const Real_wp& x);
    const Real_wp& GetNodalPointInsideEdge(int i) const; 

    int FindEdgeWithExtremities(int n1_coarse, int n2_coarse) const;
    int FindEdgeRefWithExtremities(int n1_coarse, int n2_coarse) const;
    
    /* Methods to clear arrays */
    
    void ClearConnectivity();
    void ClearCurves();
    void Clear();
    
  };
    
} // namespace Montjoie

#define MONTJOIE_FILE_MESH_BOUNDARIES_HXX
#endif
