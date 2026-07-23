#ifndef MONTJOIE_FILE_NUMBER_MAP_HXX

namespace Montjoie
{
  
  //! this class is useful to number the mesh 
  class NumberMap
  {
  protected :
    //! number of dofs on each vertex
    int nb_dof_vertex;
    //! number of dofs on each edge
    IVect nb_dof_edge;
    //! number of dofs on each quad
    IVect nb_dof_quad; 
    //! number of dofs on each triangle
    IVect nb_dof_tri;
    //! number of dofs on each hex
    IVect nb_dof_hexa;
    //! number of dofs on each tet
    IVect nb_dof_tetra;
    //! number of dofs on each pyramid
    IVect nb_dof_pyramid;
    //! number of dofs on each pyramid
    IVect nb_dof_wedge;
    //! rotation of surfacic dofs for quadrangular faces
    Vector<Matrix<int> > FacesDof_Rotation_Quad;
    //! rotation of surfacic dofs for triangular faces
    Vector<Matrix<int> > FacesDof_Rotation_Tri;
    //! signs of dofs after rotation
    Vector<Matrix<bool> > SignDof_Rotation_Quad, SignDof_Rotation_Tri;
    // What becomes the dof numbers of an edge, if the edge is flipped
    // 1,2..k -> k ... 1 ?
    //! symmetry of dofs on edges
    Vector<IVect> EdgesDof_Symmetry;
    Vector<Vector<bool> > EdgesDof_SkewSymmetry;
    bool only_hexahedral; //!< true if only hexahedral elements
    int dg_element; //!< type of dg formulation used
    
    //! quadrature points on each edge
    Vector<VectReal_wp> points_quadrature_edge;
    //! quadrature points on each triangle
    Vector<VectR2> points_quadrature_tri;
    //! quadrature points on each quadrangle
    Vector<VectR2> points_quadrature_quad;
    //! 1/2 omega_i
    Vector<VectReal_wp> half_weight_edge;
    //! 1/2 omega_i
    Vector<VectReal_wp> half_weight_tri;
    //! 1/2 omega_i
    Vector<VectReal_wp> half_weight_quad;

    //! quadrature points of an edge after rotation
    Vector<Matrix<int> > EdgesQuad_Rotation;
    //! quadrature points of quadrilaterals after rotation
    Vector<Matrix<int> > FacesQuad_Rotation_Quad;
    //! quadrature points of triangles after rotation
    Vector<Matrix<int> > FacesQuad_Rotation_Tri;
    
    //! true if the dofs on the face are not invariant by rotation
    /*!
      In that case, you have to provide a dense matrix for projection 
      on local degrees of freedom when the global face is oriented differently from local face
    */
    Vector<bool> general_combination_tri, general_combination_quad;
    //! dense operator involved for triangular faces
    //Vector< Vector<Matrix<Real_wp> > > coef_combination_tri, invCoef_combination_tri;
    //! dense operator involved for quadrilateral faces
    //Vector< Vector<Matrix<Real_wp> > > coef_combination_quad, invCoef_combination_quad;

    //! sparse operator involved for triangular faces
    Vector< Vector<Matrix<Real_wp, General, RowSparse> > >
    coef_combination_tri, invCoef_combination_tri;
    
    //! sparse operator involved for quadrilateral faces
    Vector< Vector<Matrix<Real_wp, General, RowSparse> > >
    coef_combination_quad, invCoef_combination_quad;
    
  public:
    static IVect NbNodesGmshEntity;
    static IVect NbVerticesGmshEntity;
    static IVect OrderGmshEntity;
    static IVect DimensionGmshEntity;
    
    NumberMap();

    size_t GetMemorySize() const;
    
    /*************************
     * Local to Global stuff *
     *************************/
    
    template<class Dimension>
    IVect GetDofNumberOnElement(const MeshNumbering<Dimension>& mesh, int) const;
    
    template<class Dimension, class T, int nb_unknowns>
    void GetLocalUnknownVector(const MeshNumbering<Dimension>& mesh,
			       const Vector<T>& B, int i,
			       TinyVector<Vector<T>, nb_unknowns>& Uloc) const;

    template<class Dimension, class T>
    void GetLocalUnknownVector(const MeshNumbering<Dimension>& mesh,
			       const Vector<T>& B, int i,
			       Vector<Vector<T> >& Uloc) const;

    template<class Dimension, class T, int nb_unknowns>
    void AddLocalUnknownVector(const MeshNumbering<Dimension>& mesh,
			       const T& alpha, TinyVector<Vector<T>, nb_unknowns>& Uloc,
			       int i, Vector<T>& C) const;
        
    template<class Dimension, class T>
    void ModifyLocalUnknownVector(const MeshNumbering<Dimension>&,
				  Vector<Vector<T> >& U, int iquad) const;

    template<class Dimension, class T>
    void ModifyLocalUnknownVector(const MeshNumbering<Dimension>&,
				  Vector<T>& U, int iquad) const;

    template<class Dimension, class T, int nb_unknowns>
    void ModifyLocalUnknownVector(const MeshNumbering<Dimension>&,
				  Vector<TinyVector<T, nb_unknowns> >& U, int iquad) const;

    template<class Dimension, class T>
    void ModifyLocalComponentVector(const MeshNumbering<Dimension>&,
				    Vector<Vector<T> >& Eloc, int iquad) const;

    template<class Dimension, class T>
    void ModifyLocalComponentVector(const MeshNumbering<Dimension>&,
				    Vector<T>& Eloc, int iquad) const;
    
    template<class Dimension, class T>
    void GetGlobalUnknownVector(const MeshNumbering<Dimension>&,
				Vector<Vector<T> >& U, int iquad) const;
    
    template<class Dimension, class T>
    void ModifyLocalRowMatrix(const MeshNumbering<Dimension>& mesh, 
                              VirtualMatrix<T>& mat, int iquad,
                              int nb_u = 0, int off_row = 0) const;
    
    template<class Dimension, class T>
    void ModifyLocalColumnMatrix(const MeshNumbering<Dimension>& mesh, 
                                 VirtualMatrix<T>& mat, int iquad,
                                 int nb_u = 0, int off_col = 0) const;
    
    
    /*********************************************
     * Treatment of dofs inside triangular faces *
     *********************************************/
    
  protected :
    template<class Dimension, class T, class Prop, class Storage>
    void ModifyLocalDenseRowMatrix(const MeshNumbering<Dimension>&,
                                   Matrix<T, Prop, Storage>& mat, int iquad,
                                   int nb_u = 0, int off_row = 0) const;
    
    template<class Dimension, class T, class Prop, class Storage>
    void ModifyLocalDenseColumnMatrix(const MeshNumbering<Dimension>& mesh, 
                                      Matrix<T, Prop, Storage>& mat, int iquad,
                                      int nb_u = 0, int off_col = 0) const;

    template<class Dimension, class T, class Prop, class Storage>
    void ModifyLocalSparseRowMatrix(const MeshNumbering<Dimension>&,
                                    Matrix<T, Prop, Storage>& mat, int iquad,
                                    int nb_u = 0, int off_row = 0) const;
    
    template<class Dimension, class T, class Prop, class Storage>
    void ModifyLocalSparseColumnMatrix(const MeshNumbering<Dimension>& mesh, 
                                       Matrix<T, Prop, Storage>& mat, int iquad,
                                       int nb_u = 0, int off_col = 0) const;
    
    template<class Vector1>
    inline void ProjectToGlobalDofs(Vector1& U, int rot, int r, int offset,
				    const Edge<Dimension2>& f) const {}
    
    template<class Vector1>
    inline void TransposeProjectToGlobalDofs(Vector1& U, int rot, int r, int offset,
					     const Edge<Dimension2>& f) const {}
    
    template<class Vector1>
    inline void ProjectToLocalDofs(Vector1& U, int rot, int r, int offset,
				   const Edge<Dimension2>& f) const {}

    template<class Matrix1>
    inline void ModifyRowToGlobalDofs(Matrix1& A, int rot, int r, int offset,
				      const Edge<Dimension2>& f) const {}
    
    template<class Matrix1>
    inline void ModifyColumnToGlobalDofs(Matrix1& A, int rot, int r, int offset,
					 const Edge<Dimension2>& f) const {}
    
    template<class Vector1>
    void ProjectToGlobalDofs(Vector1& U, int rot, int r, int offset,
			     const Face<Dimension3>& f) const;

    template<class Vector1>
    void TransposeProjectToGlobalDofs(Vector1& U, int rot, int r, int offset,
				      const Face<Dimension3>& f) const;

    template<class Vector1>
    void ProjectToLocalDofs(Vector1& U, int rot, int r, int offset,
			    const Face<Dimension3>& f) const;

    template<class Matrix1>
    void ModifyRowToGlobalDofs(Matrix1& A, int rot, int r, int offset,
			       const Face<Dimension3>& f) const;

    template<class Matrix1>
    void ModifyColumnToGlobalDofs(Matrix1& A, int rot, int r, int offset,
				  const Face<Dimension3>& f) const;
    
    template<class Vector1>
    void ApplyOperator(const Matrix<Real_wp, General, RowSparse>& A,
                       Vector1& U, int offset, int nb_dof) const;

    template<class Vector1>
    void ApplyTransposeOperator(const Matrix<Real_wp, General, RowSparse>& A,
				Vector1& U, int offset, int nb_dof) const;

    template<class Matrix1>
    void ApplyOperatorRow(const Matrix<Real_wp, General, RowSparse>& A,
                          Matrix1& U, int offset, int nb_dof) const;

    template<class Matrix1>
    void ApplyOperatorColumn(const Matrix<Real_wp, General, RowSparse>& A,
                             Matrix1& U, int offset, int nb_dof) const;

  public:

    template<class T>
    void ApplyOperatorFace(const SeldonTranspose& trans,
			   int r, int rot, const Edge<Dimension2>& e,
			   const Vector<T>& U, Vector<T>& Un) const;
    
    template<class T>
    void ApplyOperatorFace(const SeldonTranspose& trans,
			   int r, int rot, const Face<Dimension3>& f,
			   const Vector<T>& U, Vector<T>& V) const;
    

    /*******************************
     * Counting degrees of freedom *
     *******************************/    
    
    int GetNbDofVertex(int order) const;
    int GetNbDofEdge(int order) const;
    int GetNbDofTriangle(int order) const;
    int GetNbDofQuadrangle(int order) const;
    int GetNbDofTetrahedron(int order) const;
    int GetNbDofPyramid(int order) const;
    int GetNbDofWedge(int order) const;
    int GetNbDofHexahedron(int order) const;
    
    void SetNbDofVertex(int order, int);
    void SetNbDofEdge(int order, int);
    void SetNbDofTriangle(int order, int);
    void SetNbDofQuadrangle(int order, int);
    void SetNbDofTetrahedron(int order, int);
    void SetNbDofPyramid(int order, int);
    void SetNbDofWedge(int order, int);
    void SetNbDofHexahedron(int order, int);

    template<class Dimension>
    int GetNbDofElement(int order, const Edge<Dimension>&) const;

    template<class Dimension>
    int GetNbDofElement(int order, const Face<Dimension>&) const;
    
    int GetNbDofElement(int order, const Volume&) const;

    int GetNbDofBoundary(int order, const Edge<Dimension2>&) const;
    int GetNbDofBoundary(int order, const Face<Dimension3>&) const;

    void SetOrder(const Mesh<Dimension2>& mesh, int order);
    
#ifdef MONTJOIE_WITH_THREE_DIM
    void SetOrder(const Mesh<Dimension3>& mesh, int order);
#endif

    /*******************************************
     * Symmetry of degrees of freedom on edges *
     *******************************************/
    
    int GetSymmetricEdgeDof(int order, int k) const;
    bool IsSkewSymmetricEdgeDof(int order, int k) const;

    void SetOppositeEdgesDofSymmetry(int r, int n);
    void SetEqualEdgesDofSymmetry(int r, int n);
    void SetEdgesDofSymmetry(int r, int i, int j, bool skew_sym);
    void SetSkewSymmetryEdgesDof(int r, const Vector<bool>&);
    void SetAllEdgesDofToSkewSymmetric(int r);
    void SetOddEdgesDofToSkewSymmetric(int r);
    void SetEvenEdgesDofToSkewSymmetric(int r);
    
    /**********************************************
     * Degrees of freedom on faces after rotation *
     **********************************************/
        
    bool DofInvariantByRotation(int order, const Face<Dimension3>&) const;    
    
    void SetFacesDofRotationTri(int, const Matrix<int>& );
    void SetFacesDofRotationQuad(int, const Matrix<int>& );
    void SetFacesDofRotationTri(int, const Vector<Matrix<Real_wp> >& );
    void SetFacesDofRotationQuad(int, const Vector<Matrix<Real_wp> >& );
    const Vector<Matrix<Real_wp, General, RowSparse> >& GetFacesDofRotationTri(int) const;
    
    void SetSignDofRotationTri(int, const Matrix<bool>& );
    void SetSignDofRotationQuad(int, const Matrix<bool>& );
    
    bool IsNegativeDof(int rot, int r, int k, const Face<Dimension3>&) const;
    int GetRotationFaceDof(int rot, int r, int k, const Edge<Dimension2>&) const;
    int GetRotationFaceDof(int rot, int r, int k, const Face<Dimension3>&) const;
        
    /*********************************
     * Methods for quadrature points *
     *********************************/
    
    int GetNbPointsQuadBoundary(int order, const Edge<Dimension2>&) const;
    int GetNbPointsQuadBoundary(int order, const Face<Dimension3>&) const;
    
    Real_wp GetQuadraturePoint(int rf, int k, const Edge<Dimension2>&) const;

#ifdef MONTJOIE_WITH_THREE_DIM
    const R2& GetQuadraturePoint(int rf, int k, const Face<Dimension3>& face) const;
#endif

    void ConstructQuadrature2D(const TinyVector<IVect, 4>& order, int type_quad = 0);
    void ConstructQuadrature3D(const TinyVector<IVect, 4>& order, int type_tri = 0,
                               int type_quad = 0);
    
    const VectReal_wp& GetEdgeQuadrature(int order) const;
    const VectR2& GetTriangleQuadrature(int order) const;
    const VectR2& GetQuadrangleQuadrature(int order) const;
    
    void SetPointsQuadratureTriangle(int order, const VectR2&);
    void SetPointsQuadratureQuadrangle(int order, const VectR2&);
    void SetPointsQuadratureEdge(int order, const VectReal_wp&);
    
    void SetFluxWeightEdge(int order, const VectReal_wp& );
    const VectReal_wp& GetFluxWeight(int order, const Edge<Dimension2>& elt) const;

    void SetFluxWeightTriangle(int order, const VectReal_wp& );
    void SetFluxWeightQuadrangle(int order, const VectReal_wp& );
    const VectReal_wp& GetFluxWeight(int order, const Face<Dimension3>& elt) const;
    
    const Matrix<int>& GetRotationQuadraturePoints(int r, const Edge<Dimension2>&) const;
    const Matrix<int>& GetRotationQuadraturePoints(int r, const Face<Dimension3>&) const;
    const Matrix<int>& GetRotationQuadraturePoints(int r, int nb_vertices) const;
    

    /***************
     * Gmsh format *
     ***************/
    
    static void GetGmshTriangularNumbering(int r, IVect& num);
    static void GetGmshQuadrilateralNumbering(int r, IVect& num);

#ifdef MONTJOIE_WITH_THREE_DIM
    static void GetGmshTetrahedralNumbering(int r, IVect& num);
    static void GetGmshPyramidalNumbering(int r, IVect& num);
    static void GetGmshPrismaticNumbering(int r, IVect& num);
    static void GetGmshHexahedralNumbering(int r, IVect& num);
#endif
    
    template<class Dimension>
    static int GetGmshEntityNumber(const Edge<Dimension>& e, int order);

    template<class Dimension>
    static int GetGmshEntityNumber(const Face<Dimension>& e, int order);

    static int GetGmshEntityNumber(const Volume& e, int order);
    
    
    /*****************
     * Other methods *
     *****************/
    
    int FormulationDG() const;
    void SetFormulationDG(int);
    
    void Clear();
        
    friend ostream& operator <<(ostream& out, const NumberMap& nmap);
  };
  
}

#define MONTJOIE_FILE_NUMBER_MAP_HXX
#endif
