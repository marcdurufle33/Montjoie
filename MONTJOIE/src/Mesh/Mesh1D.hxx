#ifndef MONTJOIE_FILE_MESH1D_HXX

namespace Montjoie
{
  
  //! Class implementing 1-D mesh 
  /*!
    It contains both the geometrical mesh and the dof numbering
   */
  template<class T>
  class Mesh<Dimension1, T>
  {
  protected :
    //! Nodle(i,j) returns the degree of freedom number j of the element i
    Matrix<int> Nodle;
    //! x_i of each point of the 1D mesh
    Vector<T> Vertices;
    //! elements of the mesh (intervals [x_i, x_{i+1}]
    Vector<Edge<Dimension1> > elements;
    //! order of discretization
    int order;
    //! number of degrees of freedom
    int nodl;
    //! if true, the mesh is periodic
    bool periodic_ext;
    //! type of addition
    int pml_add_type;
    //! number of layers in PML domain
    int nb_layers_pml, ref_pml;
    //! PML thickness
    T thicknessPML;
    
    //! level of refinement for each required vertex
    IVect level_refinement_vertex;
    //! position of vertex around which refinement is required
    Vector<T> position_refinement_vertex;
    //! ratio of refinement for each required vertex
    Vector<T> ratio_refinement_vertex;

    //! for each layer, number of cells to be added and associated reference
    IVect nb_layers_to_add, reference_layers_to_add;
    //! for each layer, final position (xmax)
    Vector<T> position_layers_to_add;

  public :
    //! pml can be added on the right, left or both sides
    enum {PML_NO, PML_NEGATIVE_SIDE, PML_POSITIVE_SIDE, PML_BOTH_SIDES};
    //! global dof numbers
    IVect GlobDofNumber_Subdomain;
    //! verbose level
    int print_level;
    
    Mesh();
    
    const Vector<T>& Vertex() const;
    const T& Vertex(int i) const;
    T& Vertex(int i);
    
    const Edge<Dimension1>& Element(int i) const;
    Edge<Dimension1>& Element(int i);
    
    template<class Vector1>
    void GetVerticesElement(int num_face, Vector1& points) const;
    
    const T& GetXmin() const;
    const T& GetXmax() const;
    
    int GetNumberDof(int num_elt, int num_dof) const;
    
    int GetOrder() const;
    int GetNbDof() const;
    int GetNbElt() const;
    int GetNbVertices() const;
    void SetOrder(int);
    void SetPeriodicExtremity(bool per = true);
    int GetOrderElement(int) const;

    T GetMeshSize() const;
    
    // Read the mesh
    void Read(const string & file_mesh);
    void Write(const string & file_name);
    
    size_t GetMemorySize() const;

    void ConstructMesh(const Vector<string>& parameters, Real_wp&, Real_wp&);
    
    // interval subdivised regularly
    void CreateRegularMesh(const T& xmin, const T& xmax, int nbPoints_x, int ref_domain);
    void CreateNonRegularMesh(const Vector<T>& Points, const IVect& ref_domain);
    void CreateLayeredMesh(const Vector<T>& pos, const IVect& n, const IVect& ref);
    
    void AddLayersInput(const T& pos, int N, int ref);
    void ClearLayersInput();
    void AddRefinementVertex(const T& pos, int lvl, const T& ratio);
    void ClearRefinementVertex();
    
    // Number the mesh, fill the array Nodle !
    void NumberMesh(bool dg_form = false);
    
    void SubdivideMesh(Mesh<Dimension1, T>& fine_mesh,
		       const IVect& factor_refinement) const;

    void CreateSubmesh(Mesh<Dimension1, T>& sub_mesh, int n0, int n1) const;
    
    void Clear();
    
    /* Methods for PML */
    
    const T& GetThicknessPML() const;
    void SetThicknessPML(const T& delta);
    int GetTypeAdditionPML(int) const;
    int GetNbLayersPML() const;
    void SetAdditionPML(int type_add, int nb_layers, int ref = 0);

    void AddPMLElements(int nb_div);    
    
  };
  
  template<class T>
  ostream& operator <<(ostream& out, const Mesh<Dimension1, T>& mesh);
  
} // namespace Montjoie  

#define MONTJOIE_FILE_MESH1D_HXX
#endif
