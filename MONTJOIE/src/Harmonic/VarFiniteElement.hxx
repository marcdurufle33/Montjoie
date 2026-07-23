#ifndef MONTJOIE_FILE_VAR_FINITE_ELEMENT_HXX

namespace Montjoie
{
  
  template<class Dimension, int type>
  class VarFiniteElementEnum_Base
  {
  protected :
    static Vector<map<string, int> > list_element;
    
  public :
    static int GetIdentityNumber(const string& name_elt, int t);
    
  };

  template<class Dimension, int type>
  Vector<map<string, int> > VarFiniteElementEnum_Base<Dimension, type>::list_element;

  template<class Dimension, int type>
  class VarFiniteElementEnum : public VarFiniteElementEnum_Base<Dimension, type>
  {
  public :
    static void InitStaticData();
    static ElementReference<Dimension, type>* GetNewReferenceElement(int id);
    
  };

#ifdef MONTJOIE_WITH_NODAL_H1

#ifdef MONTJOIE_WITH_ONE_DIM
  template<>
  class VarFiniteElementEnum<Dimension1, 1> : public VarFiniteElementEnum_Base<Dimension1, 1>
  {
  public :
    enum {_None, _EdgeLobatto, _EdgeRadau, _EdgeGauss, _EdgeHierarchic};

    static void InitStaticData();
    static ElementReference<Dimension1, 1>* GetNewReferenceElement(int id);
    
  };
#endif

#ifdef MONTJOIE_WITH_TWO_DIM
  template<>
  class VarFiniteElementEnum<Dimension2, 1> : public VarFiniteElementEnum_Base<Dimension2, 1>
  {
  public :
    enum {_None, _TriangleClassical, _QuadrangleGauss, _TriangleLobatto, _QuadrangleLobatto,
	  _TriangleHierarchic, _QuadrangleHierarchic, _QuadrangleDgGauss,
	  _TriangleDgOrtho, _QuadrangleRadau};

    static void InitStaticData();
    static ElementReference<Dimension2, 1>* GetNewReferenceElement(int id);
    static ElementReference<Dimension1, 1>* GetNewSurfaceElement(int r, int rg, int rq, int id);
    
  };
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
  template<>
  class VarFiniteElementEnum<Dimension3, 1> : public VarFiniteElementEnum_Base<Dimension3, 1>
  {
  public :
    enum {_None, _TetrahedronClassical, _PyramidClassical, _WedgeClassical, _HexahedronGauss,
	  _PyramidLobatto, _WedgeLobatto, _HexahedronLobatto,
	  _TetrahedronHierarchic, _PyramidHierarchic, _WedgeHierarchic, _HexahedronHierarchic,	  
	  _WedgeDgClassical, _HexahedronDgGauss, _TetrahedronDgOrtho, _PyramidDgOrtho, _WedgeDgOrtho,
	  _PyramidDgLegendre, _WedgeDgLegendre, _HexahedronDgLegendre};
    
    static void InitStaticData();
    static ElementReference<Dimension3, 1>* GetNewReferenceElement(int id);
    static ElementReference<Dimension2, 1>* GetNewSurfaceElement(int r, int rg, int rq, int id);    
    
  };
#endif

#endif


#ifdef MONTJOIE_WITH_NODAL_HCURL

#ifdef MONTJOIE_WITH_TWO_DIM
  template<>
  class VarFiniteElementEnum<Dimension2, 2> : public VarFiniteElementEnum_Base<Dimension2, 2>
  {
  public :
    enum {_None, _TriangleHcurlFirstFamily, _QuadrangleHcurlFirstFamily, _QuadrangleHcurlGaussFirstFamily,
	  _TriangleHcurlOptimalFirstFamily, _QuadrangleHcurlOptimalFirstFamily,
	  _TriangleHcurlSecondFamily, _QuadrangleHcurlLobatto, _TriangleHcurlHierarchic,
	  _QuadrangleHcurlHpFirstFamily, _TriangleHcurlOptimalHpFirstFamily, 
	  _QuadrangleHcurlOptimalHpFirstFamily};
    
    static void InitStaticData();
    static ElementReference<Dimension2, 2>* GetNewReferenceElement(int id);
    static ElementReference<Dimension1, 1>* GetNewSurfaceElement(int r, int rg, int rq, int id);
    
  };
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
  template<>
  class VarFiniteElementEnum<Dimension3, 2> : public VarFiniteElementEnum_Base<Dimension3, 2>
  {
  public :
    enum {_None, _TetrahedronHcurlFirstFamily, _PyramidHcurlFirstFamily,
          _WedgeHcurlFirstFamily, _HexahedronHcurlFirstFamily,
	  _TetrahedronHcurlOptimalFirstFamily, _PyramidHcurlOptimalFirstFamily,
          _WedgeHcurlOptimalFirstFamily, _HexahedronHcurlOptimalFirstFamily,
	  _TetrahedronHcurlLobatto, _PyramidHcurlLobatto, _WedgeHcurlLobatto, _HexahedronHcurlLobatto,
	  _TetrahedronHcurlOptimalHpFirstFamily, _PyramidHcurlOptimalHpFirstFamily,
	  _WedgeHcurlOptimalHpFirstFamily, _HexahedronHcurlOptimalHpFirstFamily, 
	  _PyramidHcurlHpFirstFamily, _WedgeHcurlHpFirstFamily, _HexahedronHcurlHpFirstFamily};
    
    static void InitStaticData();
    static ElementReference<Dimension3, 2>* GetNewReferenceElement(int id);    
    static ElementReference<Dimension2, 2>* GetNewSurfaceElement(int r, int rg, int rq, int id);
    
  };
#endif

#endif


#ifdef MONTJOIE_WITH_NODAL_HDIV

#ifdef MONTJOIE_WITH_TWO_DIM
  template<>
  class VarFiniteElementEnum<Dimension2, 3> : public VarFiniteElementEnum_Base<Dimension2, 3>
  {
  public :
    enum {_None, _TriangleHdivFirstFamily, _QuadrangleHdivFirstFamily,
	  _TriangleHdivOptimalFirstFamily, _QuadrangleHdivOptimalFirstFamily,
	  _QuadrangleHdivHpFirstFamily,
	  _TriangleHdivOptimalHpFirstFamily, _QuadrangleHdivOptimalHpFirstFamily};
    
    static void InitStaticData();
    static ElementReference<Dimension2, 3>* GetNewReferenceElement(int id);
    static ElementReference<Dimension1, 1>* GetNewSurfaceElement(int r, int rg, int rq, int id);
    
  };
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
  template<>
  class VarFiniteElementEnum<Dimension3, 3> : public VarFiniteElementEnum_Base<Dimension3, 3>
  {
  public :
    enum {_None, _TetrahedronHdivFirstFamily, _PyramidHdivFirstFamily,
          _WedgeHdivFirstFamily, _HexahedronHdivFirstFamily,
	  _TetrahedronHdivOptimalFirstFamily, _PyramidHdivOptimalFirstFamily,
          _WedgeHdivOptimalFirstFamily, _HexahedronHdivOptimalFirstFamily,
	  _TetrahedronHdivOptimalHpFirstFamily, _PyramidHdivOptimalHpFirstFamily,
          _WedgeHdivOptimalHpFirstFamily, _HexahedronHdivOptimalHpFirstFamily,
	  _PyramidHdivHpFirstFamily, _WedgeHdivHpFirstFamily, _HexahedronHdivHpFirstFamily};
    
    static void InitStaticData();
    static ElementReference<Dimension3, 3>* GetNewReferenceElement(int id);    
    static ElementReference<Dimension2, 1>* GetNewSurfaceElement(int r, int rg, int rq, int id);
    
  };
#endif

#endif


  template<class Dimension>
  class VarFiniteElement;
  
  class VarFiniteElement_Base
  {
  protected:
    int nb_elt_, nb_surf_, nb_surf_neigh_;
    
     bool use_piola_transform;
    
  public :
    VarFiniteElement_Base();
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
#ifdef MONTJOIE_WITH_ONE_DIM
    void AddFiniteElement1D(const string& name_elt, int order, const Vector<bool>&,
			    Mesh<Dimension1>& mesh);
#endif
    
#ifdef MONTJOIE_WITH_THREE_DIM
    void GetTypeIntegrationBoundary(const Mesh<Dimension3>& mesh, const MeshNumbering<Dimension3>& mesh_num,
                                    const Vector<MeshNumbering<Dimension3>* >& other_mesh_num,
				    int& type_integration_edge,
				    int& type_integration_tri, int& type_integration_quad);
#endif
    
#ifdef MONTJOIE_WITH_TWO_DIM
    void GetTypeIntegrationBoundary(const Mesh<Dimension2>& mesh, const MeshNumbering<Dimension2>& mesh_num,
                                    const Vector<MeshNumbering<Dimension2>* >& other_mesh_num,
				    int& type_integration_edge,
				    int& type_integration_tri, int& type_integration_quad);
#endif
    
  };

  //! class handling finite element
  template<class Dimension>
  class VarFiniteElement : public VarFiniteElement_Base
  {
    friend class VarFiniteElement_Base;

  private :
    // volume finite elements
    Vector<int> id_stored_reference_element, type_stored_reference_element;
    Vector<ElementReference_Dim<Dimension>* > stored_reference_element;
    Vector<ElementReference_Dim<Dimension>* > reference_element; //!< array of finite elements
    
    // surface finite elements
    typedef ElementReference_Dim<typename Dimension::DimensionBoundary> ElementReferenceTrace;
    Vector<int> id_stored_surface_element, type_stored_surface_element;
    Vector<ElementReferenceTrace*> stored_surface_element;
    Vector<ElementReferenceTrace*> surface_element; //!< array of finite elements    
    
    //! finite elements associated with other processors
    Vector<ElementReference_Dim<Dimension>* > reference_neighbor_element; 

  protected:
    static int GetIdentityNumber(const string& name, int t, int type);

    ElementReferenceTrace* GetSurfaceFiniteElementPtr(int r, int id_elt, int type_elt) const;
    ElementReference_Dim<Dimension>* GetReferenceFiniteElementPtr(int r, int id_elt, int type_elt) const;

    ElementReferenceTrace* GetNewSurfaceElement(int r, int rgeom, int rquad, int id_elt, int type_elt);

    ElementReference_Dim<Dimension>* GetNewReferenceElement(int id_elt, int type_elt);

    void ConstructFiniteElement(int r, int rgeom, int rquad, int type_elt,
                                ElementReference_Dim<Dimension>* elt);
    
  public :
    ~VarFiniteElement();
    
    void ClearFiniteElement();
    
    void GetMemoryUsed(map<string, size_t>& var) const;

    void AddSurfaceFiniteElement(const string& name_elt, TinyVector<IVect, 4>& order,
				 Vector<bool>& change_elt, Mesh<Dimension>& mesh,
				 MeshNumbering<Dimension>& mesh_num, int r_over_quad,
				 int discontinuous_formulation, int type = 1);

    void AddSurfaceFiniteElement(const string& name_elt, TinyVector<IVect, 4>& order, int type,
                                 const Vector<string>& other_name,
                                 Vector<TinyVector<IVect, 4> >& other_order, Vector<int>& other_type,
                                 Vector<bool>& change_elt, Mesh<Dimension>& mesh,
                                 MeshNumbering<Dimension>& mesh_num,
                                 Vector<MeshNumbering<Dimension>* >& other_mesh_num,
                                 Vector<int>& r_over_quad, int discontinuous_formulation);
    
    void AddFiniteElement(const string& name_elt, TinyVector<IVect, 4>& order,
			  Vector<bool>& change_elt, Mesh<Dimension>& mesh,
			  MeshNumbering<Dimension>& mesh_num, int r_over_quad,
			  int discontinuous_formulation, int type = 1);

    void AddFiniteElement(const string& name_elt, TinyVector<IVect, 4>& order, int type,
                          const Vector<string>& other_name,
                          Vector<TinyVector<IVect, 4> >& other_order,
                          const Vector<int>& other_type,
			  Vector<bool>& change_elt, Mesh<Dimension>& mesh,
			  MeshNumbering<Dimension>& mesh_num,
                          Vector<MeshNumbering<Dimension>* >& other_mesh_num,
                          Vector<int>& r_over_quad,
			  int discontinuous_formulation,
                          IVect& other_discont_form);

    void UpdateInterpolationElement(Mesh<Dimension>& mesh, MeshNumbering<Dimension>& mesh_num);
    
    void CopyFiniteElement(const VarFiniteElement<Dimension>&,
                           const Mesh<Dimension>& mesh,
                           const MeshNumbering<Dimension>& mesh_num,
                           const Vector<MeshNumbering<Dimension>* >& other_mesh_num,
                           int type, const Vector<int>& other_type);
    
    void SetFiniteElement(const string&, const Vector<bool>&, int type, int n = 0);
    
    void GetReferenceElement(Vector<const ElementReference_Dim<Dimension>* >& elt) const;
    void GetSurfaceFiniteElement(Vector<const ElementReferenceTrace*>& elt) const;
    
    void SetAxisymGeometry(bool axisym);
    
    // Inline methods
    const ElementReference_Dim<Dimension>& GetReferenceElement(int i, int n = 0) const;
    ElementReference_Dim<Dimension>& GetReferenceElement(int i, int n = 0);

    const ElementReference<Dimension, 1>& GetReferenceElementH1(int i, int n = 0) const;
    ElementReference<Dimension, 1>& GetReferenceElementH1(int i, int n = 0);

    const ElementReference<Dimension, 2>& GetReferenceElementHcurl(int i, int n = 0) const;
    ElementReference<Dimension, 2>& GetReferenceElementHcurl(int i, int n = 0);

    const ElementReference<Dimension, 3>& GetReferenceElementHdiv(int i, int n = 0) const;
    ElementReference<Dimension, 3>& GetReferenceElementHdiv(int i, int n = 0);

    const ElementReferenceTrace& GetSurfaceFiniteElement(int i, int n = 0) const;
    ElementReferenceTrace& GetSurfaceFiniteElement(int i, int n = 0);

    const ElementReference<typename Dimension::DimensionBoundary, 1>& GetSurfaceFiniteElementH1(int i, int n = 0) const;
    ElementReference<typename Dimension::DimensionBoundary, 1>& GetSurfaceFiniteElementH1(int i, int n = 0);

    const ElementReference<typename Dimension::DimensionBoundary, 2>& GetSurfaceFiniteElementHcurl(int i, int n = 0) const;
    ElementReference<typename Dimension::DimensionBoundary, 2>& GetSurfaceFiniteElementHcurl(int i, int n = 0);

    const ElementReference_Dim<Dimension>& GetNeighborReferenceElement(int i, int n = 0) const;
    ElementReference_Dim<Dimension>& GetNeighborReferenceElement(int i, int n = 0);

  };

  
}

#define MONTJOIE_FILE_VAR_FINITE_ELEMENT_HXX
#endif
