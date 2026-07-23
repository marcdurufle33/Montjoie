#ifndef MONTJOIE_FILE_NUMBER_MAP_CXX

namespace Montjoie
{  

#ifndef SELDON_WITH_COMPILED_LIBRARY
  IVect NumberMap::NbNodesGmshEntity;
  IVect NumberMap::NbVerticesGmshEntity;
  IVect NumberMap::OrderGmshEntity;
  IVect NumberMap::DimensionGmshEntity;
#endif


  //! returns the number of dofs strictly inside the element
  template<class Dimension>
  int NumberMap::GetNbDofElement(int order, const Face<Dimension>& elt) const
  {
    int nv = elt.GetNbVertices();
    if (nv == 3)
      {
	if (nb_dof_tri.GetM() > order)
	  return nb_dof_tri(order);
	
	return 0;
      }
    
    if (nb_dof_quad.GetM() > order)
      return nb_dof_quad(order);
    
    return 0;
  }
  
  
  //! returns the number of dofs strictly inside the element
  int NumberMap::GetNbDofElement(int order, const Volume& elt) const
  {
    int nv = elt.GetNbVertices();
    if (nv == 4)
      {
	if (nb_dof_tetra.GetM() > order)
	  return nb_dof_tetra(order);
	
	return 0;
      }
    else if (nv == 5)
      {
	if (nb_dof_pyramid.GetM() > order)
	  return nb_dof_pyramid(order);
	
	return 0;
      }
    if (nv == 6)
      {
	if (nb_dof_wedge.GetM() > order)
	  return nb_dof_wedge(order);
	
	return 0;
      }
    
    if (nb_dof_hexa.GetM() > order)
      return nb_dof_hexa(order);
    
    return 0;
  }

  
  //! sets the number of dofs per edge for a given order
  void NumberMap::SetNbDofEdge(int order, int nb_dof)
  {
    if (nb_dof_edge.GetM() <= order)
      {
	int m = nb_dof_edge.GetM();
	nb_dof_edge.Resize(order+1);
	for (int i = m; i <= order; i++)
	  nb_dof_edge(i) = 0;
      }
    
    nb_dof_edge(order) = nb_dof;
  }
  

  //! sets the number of dofs per triangle for a given order  
  void NumberMap::SetNbDofTriangle(int order, int nb_dof)
  {
    if (nb_dof_tri.GetM() <= order)
      {
	int m = nb_dof_tri.GetM();
	nb_dof_tri.Resize(order+1);
	for (int i = m; i <= order; i++)
	  nb_dof_tri(i) = 0;
      }
    
    nb_dof_tri(order) = nb_dof;
  }
  
  
  //! sets the number of dofs per quadrangle for a given order
  void NumberMap::SetNbDofQuadrangle(int order, int nb_dof)
  {
    if (nb_dof_quad.GetM() <= order)
      {
	int m = nb_dof_quad.GetM();
	nb_dof_quad.Resize(order+1);
	for (int i = m; i <= order; i++)
	  nb_dof_quad(i) = 0;
      }
    
    nb_dof_quad(order) = nb_dof;
  }
  
  
  //! sets the number of dofs per tetrahedron for a given order
  void NumberMap::SetNbDofTetrahedron(int order, int nb_dof)
  {
    if (nb_dof_tetra.GetM() <= order)
      {
	int m = nb_dof_tetra.GetM();
	nb_dof_tetra.Resize(order+1);
	for (int i = m; i <= order; i++)
	  nb_dof_tetra(i) = 0;
      }
    
    nb_dof_tetra(order) = nb_dof;
  }
  
  
  //! sets the number of dofs per pyramid for a given order
  void NumberMap::SetNbDofPyramid(int order, int nb_dof)
  {
    if (nb_dof_pyramid.GetM() <= order)
      {
	int m = nb_dof_pyramid.GetM();
	nb_dof_pyramid.Resize(order+1);
	for (int i = m; i <= order; i++)
	  nb_dof_pyramid(i) = 0;
      }
    
    nb_dof_pyramid(order) = nb_dof;
  }
  
  
  //! sets the number of dofs per prism for a given order
  void NumberMap::SetNbDofWedge(int order, int nb_dof)
  {
    if (nb_dof_wedge.GetM() <= order)
      {
	int m = nb_dof_wedge.GetM();
	nb_dof_wedge.Resize(order+1);
	for (int i = m; i <= order; i++)
	  nb_dof_wedge(i) = 0;
      }
    
    nb_dof_wedge(order) = nb_dof;
  }
  
  
  //! sets the number of dofs per hexahedron for a given order
  void NumberMap::SetNbDofHexahedron(int order, int nb_dof)
  {
    if (nb_dof_hexa.GetM() <= order)
      {
	int m = nb_dof_hexa.GetM();
	nb_dof_hexa.Resize(order+1);
	for (int i = m; i <= order; i++)
	  nb_dof_hexa(i) = 0;
      }
    
    nb_dof_hexa(order) = nb_dof;
  }
  
  
  //! sets the quadrature points to use for a triangular face
  void NumberMap::SetPointsQuadratureTriangle(int order, const VectR2& xi)
  {
    if (points_quadrature_tri.GetM() <= order)
      {
        points_quadrature_tri.Resize(order+1);
        FacesQuad_Rotation_Tri.Resize(order+1);
      }
    
    points_quadrature_tri(order) = xi;
#ifdef MONTJOIE_WITH_THREE_DIM
    MeshNumbering<Dimension3>::
      GetRotationTriangularFace(xi, FacesQuad_Rotation_Tri(order));
#endif
  }
  
  
  //! sets the quadrature points to use for a quadrangular face
  void NumberMap::SetPointsQuadratureQuadrangle(int order, const VectR2& xi)
  {
    if (points_quadrature_quad.GetM() <= order)
      {
        points_quadrature_quad.Resize(order+1);
        FacesQuad_Rotation_Quad.Resize(order+1);
      }
        
    points_quadrature_quad(order) = xi;
#ifdef MONTJOIE_WITH_THREE_DIM
    MeshNumbering<Dimension3>::
      GetRotationQuadrilateralFace(xi, FacesQuad_Rotation_Quad(order));
#endif
  }

  
  //! sets the half of quadrature weights to use for an edge
  void NumberMap::SetFluxWeightEdge(int order, const VectReal_wp& PoidsFlux)
  {
    if (half_weight_edge.GetM() <= order)
      half_weight_edge.Resize(order+1);
    
    half_weight_edge(order) = PoidsFlux;
  }
  
  
  //! sets the half of quadrature weights to use for a triangle
  void NumberMap::SetFluxWeightTriangle(int order, const VectReal_wp& PoidsFlux)
  {
    if (half_weight_tri.GetM() <= order)
      half_weight_tri.Resize(order+1);
    
    half_weight_tri(order) = PoidsFlux;
  }
  
  
  //! sets the half of quadrature weights to use for a quadrangle
  void NumberMap::SetFluxWeightQuadrangle(int order, const VectReal_wp& PoidsFlux)
  {
    if (half_weight_quad.GetM() <= order)
      half_weight_quad.Resize(order+1);
    
    half_weight_quad(order) = PoidsFlux;
  }
  
  
  size_t NumberMap::GetMemorySize() const
  {
    size_t taille = 4*sizeof(int);
    taille += nb_dof_edge.GetMemorySize() + nb_dof_quad.GetMemorySize() + nb_dof_tri.GetMemorySize();
    taille += nb_dof_hexa.GetMemorySize() + nb_dof_tetra.GetMemorySize()
      + nb_dof_pyramid.GetMemorySize() + nb_dof_wedge.GetMemorySize();
    taille += Seldon::GetMemorySize(FacesDof_Rotation_Quad) + Seldon::GetMemorySize(FacesDof_Rotation_Tri);
    taille += Seldon::GetMemorySize(SignDof_Rotation_Quad) + Seldon::GetMemorySize(SignDof_Rotation_Tri);
    taille += Seldon::GetMemorySize(EdgesDof_Symmetry) + Seldon::GetMemorySize(EdgesDof_SkewSymmetry);

    taille += Seldon::GetMemorySize(points_quadrature_edge) + Seldon::GetMemorySize(points_quadrature_tri);
    taille += Seldon::GetMemorySize(points_quadrature_quad) + Seldon::GetMemorySize(half_weight_edge);
    taille += Seldon::GetMemorySize(half_weight_tri) + Seldon::GetMemorySize(half_weight_quad);

    taille += Seldon::GetMemorySize(EdgesQuad_Rotation) + Seldon::GetMemorySize(FacesQuad_Rotation_Quad);
    taille += Seldon::GetMemorySize(FacesQuad_Rotation_Tri);
    
    taille += general_combination_tri.GetMemorySize() + general_combination_quad.GetMemorySize();
    taille += Seldon::GetMemorySize(coef_combination_tri) + Seldon::GetMemorySize(invCoef_combination_tri);
    taille += Seldon::GetMemorySize(coef_combination_quad) + Seldon::GetMemorySize(invCoef_combination_quad);
    
    return taille;
  }

  
  // filling dof scheme for a classical nodal finite element of order r
  void NumberMap::SetOrder(const Mesh<Dimension2>& mesh, int r)
  {
    SetNbDofVertex(r, 1);
    SetNbDofEdge(r, r-1);
    SetNbDofTriangle(r, (r-1)*(r-2)/2);
    SetOppositeEdgesDofSymmetry(r, r-1);
    SetNbDofQuadrangle(r, (r-1)*(r-1));
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(r);
    this->SetPointsQuadratureEdge(r, gauss.Points());
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  // filling dof scheme for a classical nodal finite element of order r
  void NumberMap::SetOrder(const Mesh<Dimension3>& mesh, int r)
  {
    SetNbDofVertex(r, 1);
    SetNbDofEdge(r, r-1);
    SetNbDofTriangle(r, (r-1)*(r-2)/2);
    SetOppositeEdgesDofSymmetry(r, r-1);
    SetNbDofQuadrangle(r, (r-1)*(r-1));
    SetNbDofTetrahedron(r, (r-1)*(r-2)*(r-3)/6);
    SetNbDofPyramid(r, (r-1)*(r-2)*(2*r-3)/6);
    SetNbDofWedge(r, (r-1)*(r-1)*(r-2)/2);
    SetNbDofHexahedron(r, (r-1)*(r-1)*(r-1));
    
    if (r > 2)
      {
	Vector<R2> points_nodal2d_inside_tri((r-1)*(r-2)/2);
	for (int j = 3*r; j < (r+1)*(r+2)/2; j++)
	  points_nodal2d_inside_tri(j-3*r) = mesh.triangle_reference.PointsNodalND(j);
	
	Matrix<int> FacesDof_Rotation_Tri;
	MeshNumbering<Dimension3>::
	  GetRotationTriangularFace(points_nodal2d_inside_tri, FacesDof_Rotation_Tri);
	
	this->SetFacesDofRotationTri(r, FacesDof_Rotation_Tri);
      }
    
    if (r > 1)
      {
	Vector<R2> points_nodal2d_inside_quad((r-1)*(r-1));
	for (int j = 4*r; j < (r+1)*(r+1); j++)
	  points_nodal2d_inside_quad(j-4*r) = mesh.quadrangle_reference.PointsNodalND(j);
	
	Matrix<int> FacesDof_Rotation_Quad;
	MeshNumbering<Dimension3>::
	  GetRotationQuadrilateralFace(points_nodal2d_inside_quad, FacesDof_Rotation_Quad);
	
	this->SetFacesDofRotationQuad(r, FacesDof_Rotation_Quad);
      }
  }
#endif
  
  
  //! 2-D quadrature rules are set as standard Gauss rules
  void NumberMap::ConstructQuadrature2D(const TinyVector<IVect, 4>& order, int type_quad)
  {
    Globatto<Real_wp> gauss;
    VectReal_wp half_weights;
    for (int i = 0; i < order(0).GetM(); i++)
      {
        int r = order(0)(i);
        
        gauss.ConstructQuadrature(r, type_quad);
        half_weights = gauss.Weights();
        Mlt(Real_wp(1)/2, half_weights);
        
        this->SetPointsQuadratureEdge(r, gauss.Points());
        this->SetFluxWeightEdge(r, half_weights);
      }
  }
  

  //! 3-D quadrature rules are set as standard Gauss rules
  void NumberMap::ConstructQuadrature3D(const TinyVector<IVect, 4>& order,
					int type_tri, int type_quad)
  {
    VectR2 points; VectReal_wp half_weights;
    Matrix<int> FacesPoints_Rotation;
    for (int i = 0; i < order(0).GetM(); i++)
      {
        int r = order(0)(i);
        
        TriangleQuadrature::ConstructQuadrature(2*r, points, half_weights, type_tri);
        Mlt(Real_wp(0.5), half_weights);
            
        this->SetPointsQuadratureTriangle(r, points);
        this->SetFluxWeightTriangle(r, half_weights);
      }
    
    for (int i = 0; i < order(1).GetM(); i++)
      {
        int r = order(1)(i);            
        QuadrangleQuadrature::ConstructQuadrature(2*r, points, half_weights, type_quad);
        Mlt(Real_wp(0.5), half_weights);
        
        this->SetPointsQuadratureQuadrangle(r, points);
        this->SetFluxWeightQuadrangle(r, half_weights);            
      }
    
  }


  //! sets the quadrature points to use for an edge
  void NumberMap::SetPointsQuadratureEdge(int order, const VectReal_wp& xi)
  {
    if (points_quadrature_edge.GetM() <= order)
      {
        points_quadrature_edge.Resize(order+1);
        EdgesQuad_Rotation.Resize(order+1);
      }
    
    points_quadrature_edge(order) = xi;
    EdgesQuad_Rotation(order).Reallocate(2, xi.GetM());
    for (int i = 0; i < xi.GetM(); i++)
      {
        EdgesQuad_Rotation(order)(0, i) = i;
        EdgesQuad_Rotation(order)(1, i) = xi.GetM() - 1 - i;
      }
  }
  
  
  //! After symmetry, the dof numbers are reversed without sign
  void NumberMap::SetOppositeEdgesDofSymmetry(int i, int m)
  {
    EdgesDof_SkewSymmetry.Resize(nb_dof_edge.GetM());
    EdgesDof_Symmetry.Resize(nb_dof_edge.GetM());
    int r = m - 1;
    if (r >= 0)
      {
	EdgesDof_Symmetry(i).Reallocate(r+1);
	EdgesDof_SkewSymmetry(i).Reallocate(r+1);
	for (int j = 0; j <= r; j++)
	  {
	    EdgesDof_Symmetry(i)(j) = r-j;
	    EdgesDof_SkewSymmetry(i)(j) = false;
          }
      }    
  }
  
  
  //! After symmetry, the dof numbers are the same, with opposite sign for even dofs
  void NumberMap::SetEqualEdgesDofSymmetry(int i, int m)
  {
    EdgesDof_SkewSymmetry.Resize(nb_dof_edge.GetM());
    EdgesDof_Symmetry.Resize(nb_dof_edge.GetM());
    int r = m - 1;
    if (r >= 0)
      {
	EdgesDof_Symmetry(i).Reallocate(r+1);
	EdgesDof_SkewSymmetry(i).Reallocate(r+1);
        for (int j = 0; j <= r; j++)
	  {
	    EdgesDof_Symmetry(i)(j) = j;
	    EdgesDof_SkewSymmetry(i)(j) = (j%2 == 0);
          }
      }
  }

  
  //! provides dof numbers after rotation of a quadrilateral face
  void NumberMap::SetFacesDofRotationQuad(int order, const Matrix<int>& num)
  {
    if (FacesDof_Rotation_Quad.GetM() <= order)
      {
	FacesDof_Rotation_Quad.Resize(order+1);
	SignDof_Rotation_Quad.Resize(order+1);
       }

    if (general_combination_quad.GetM() <= order)
      general_combination_quad.Resize(order+1);
    
    FacesDof_Rotation_Quad(order) = num;
    SignDof_Rotation_Quad(order).Reallocate(num.GetM(), num.GetN());
    SignDof_Rotation_Quad(order).Fill(false);
    general_combination_quad(order) = false;
  }

  
  //! provides dof numbers after rotation of a triangular face
  void NumberMap::SetFacesDofRotationTri(int order, const Matrix<int>& num)
  {
    if (FacesDof_Rotation_Tri.GetM() <= order)
      {
	FacesDof_Rotation_Tri.Resize(order+1);
	SignDof_Rotation_Tri.Resize(order+1);        
      }

    if (general_combination_tri.GetM() <= order)
      general_combination_tri.Resize(order+1);
    
    FacesDof_Rotation_Tri(order) = num;
    SignDof_Rotation_Tri(order).Reallocate(num.GetM(), num.GetN());
    SignDof_Rotation_Tri(order).Fill(false);
    general_combination_tri(order) = false;
  }
  
  
  //! provides linear operator to apply to dofs of a quadrilateral face after rotation
  void NumberMap::SetFacesDofRotationQuad(int order, const Vector<Matrix<Real_wp> >& coef)
  {
    if (coef_combination_quad.GetM() <= order)
      {
	coef_combination_quad.Resize(order+1);
	invCoef_combination_quad.Resize(order+1);
        general_combination_quad.Resize(order+1);
      }
    
    Matrix<Real_wp> invCoef;
    Real_wp threshold = 1e-10;
    coef_combination_quad(order).Reallocate(8);
    invCoef_combination_quad(order).Reallocate(8);
    for (int i = 0; i < 8; i++)
      {
        invCoef = coef(i);        
        GetInverse(invCoef);
        ConvertToSparse(coef(i), coef_combination_quad(order)(i), threshold);
        ConvertToSparse(invCoef, invCoef_combination_quad(order)(i), threshold);
      }

    general_combination_quad(order) = true;
  }


  //! provides linear operator to apply to dofs of a triangular face after rotation
  void NumberMap::SetFacesDofRotationTri(int order, const Vector<Matrix<Real_wp> >& coef)
  {
    if (coef_combination_tri.GetM() <= order)
      {
	coef_combination_tri.Resize(order+1);
	invCoef_combination_tri.Resize(order+1);
        general_combination_tri.Resize(order+1);
      }

    Matrix<Real_wp> invCoef;
    Real_wp threshold = 1e-10;
    coef_combination_tri(order).Reallocate(6);
    invCoef_combination_tri(order).Reallocate(6);
    for (int i = 0; i < 6; i++)
      {
        invCoef = coef(i);        
        GetInverse(invCoef);
        ConvertToSparse(coef(i), coef_combination_tri(order)(i), threshold);
        ConvertToSparse(invCoef, invCoef_combination_tri(order)(i), threshold);
      }
    
    general_combination_tri(order) = true;
  }
   
  
  //! Applies the linear operator of the face to a vector U
  template<class Vector1>
  void NumberMap::ApplyOperator(const Matrix<Real_wp, General, RowSparse>& A,
                                Vector1& U, int offset, int nb_dof) const
  {
    Vector1 V(nb_dof), Av(nb_dof);
    for (int i = 0; i < nb_dof; i++)
      V(i) = U(offset+i);
    
    Mlt(SeldonTrans, A, V, Av);
    for (int i = 0; i < nb_dof; i++)
      U(offset+i) = Av(i);
  }

  
  //! Applies the transpose of linear operator of the face to a vector U
  template<class Vector1>
  void NumberMap::ApplyTransposeOperator(const Matrix<Real_wp, General, RowSparse>& A,
                                         Vector1& U, int offset, int nb_dof) const
  {
    Vector1 V(nb_dof), Av(nb_dof);
    for (int i = 0; i < nb_dof; i++)
      V(i) = U(offset+i);
    
    Mlt(A, V, Av);
    
    for (int i = 0; i < nb_dof; i++)
      U(offset+i) = Av(i);
  }


  //! applies the linear operator for a single face
  template<class T>
  void NumberMap::ApplyOperatorFace(const SeldonTranspose& trans,
				    int r, int rot, const Edge<Dimension2>& e,
				    const Vector<T>& U, Vector<T>& V) const
  {
    if (rot)
      {
	V = U;
	return;
      }
    
    for (int k = 0; k < U.GetM(); k++)
      {
	if (EdgesDof_SkewSymmetry(r)(k))
	  V(k) = -U(k);
	else
	  V(k) = U(k);
      }    
  }


  //! applies the linear operator for a single face
  template<class T>
  void NumberMap::ApplyOperatorFace(const SeldonTranspose& trans,
				    int r, int rot, const Face<Dimension3>& f,
				    const Vector<T>& U, Vector<T>& V) const
  {
    if (f.GetNbVertices() == 3)
      {
	if (general_combination_tri(r))
	  Mlt(trans, coef_combination_tri(r)(rot), U, V);
	else
	  for (int k = 0; k < U.GetM(); k++)
	    {
	      if (SignDof_Rotation_Tri(r)(rot, k))
		V(k) = -U(k);
	      else
		V(k) = U(k);
	    }
      }
    else
      {
	if (general_combination_quad(r))
	  Mlt(trans, coef_combination_quad(r)(rot), U, V);
	else
	  for (int k = 0; k < U.GetM(); k++)
	    {
	      if (SignDof_Rotation_Quad(r)(rot, k))
		V(k) = -U(k);
	      else
		V(k) = U(k);
	    }
      }
  }

  
  //! Applies the linear operator of the face to rows of the matrix A
  template<class Matrix1>
  void NumberMap::ApplyOperatorRow(const Matrix<Real_wp, General, RowSparse>& A,
                                   Matrix1& U, int offset, int nb_dof) const
  {
    if (IsSymmetricMatrix(U))
      {
        abort();
      }
    
    int n = U.GetN();
    Matrix1 V(nb_dof, n);
    Vector<typename Matrix1::value_type> Vcol(nb_dof), Av(nb_dof);
    for (int i = 0; i < nb_dof; i++)
      for (int j = 0; j < n; j++)
        V(i, j) = U(offset+i, j);
    
    for (int j = 0; j < n; j++)
      {
        for (int k = 0; k < nb_dof; k++)
          Vcol(k) = V(k, j);
        
        Mlt(SeldonTrans, A, Vcol, Av);
        for (int i = 0; i < nb_dof; i++)
          U(offset+i, j) = Av(i);
      }
  }
  
  
  //! Applies the linear operator of the face to columns of the matrix A
  template<class Matrix1>
  void NumberMap::ApplyOperatorColumn(const Matrix<Real_wp, General, RowSparse>& A,
                                      Matrix1& U, int offset, int nb_dof) const
  {
    if (IsSymmetricMatrix(U))
      {
        abort();
      }
    
    int n = U.GetM();
    Matrix1 V(n, nb_dof);
    Vector<typename Matrix1::value_type> Vrow(nb_dof), Av(nb_dof);
    for (int i = 0; i < nb_dof; i++)
      for (int j = 0; j < n; j++)
        V(j, i) = U(j, offset+i);
    
    for (int j = 0; j < n; j++)
      {
        for (int k = 0; k < nb_dof; k++)
          Vrow(k) = V(j, k);
        
        Mlt(SeldonTrans, A, Vrow, Av);
        for (int i = 0; i < nb_dof; i++)
          U(j, offset+i) = Av(i);
      }
  }  
  
  
  //! Retrieves local components of the solution on the given element
  template<class Dimension, class T, int nb_unknowns>
  void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension>& mesh_num,
			const Vector<T>& B, int i,
			TinyVector<Vector<T>, nb_unknowns>& Eloc) const
  {
    int N = mesh_num.GetNbLocalDof(i);
    const IVect& Nodle = mesh_num.Element(i).GetNodle();
    for (int m = 0; m < nb_unknowns; m++)
      {
	Eloc(m).Reallocate(N);
	Eloc(m).Zero();
      }
    
    int Nvol = mesh_num.GetNbDof();
    int nb_u = min(int(B.GetM()/Nvol), nb_unknowns);
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      {
        int offset = Nodle(0);
        for (int j = 0; j < N; j++)
          {
            for (int m = 0; m < nb_u; m++)
              Eloc(m)(j) = B(offset + j + m*Nvol);
          }
      }
    else
      {        
        for (int j = 0; j < N; j++)
          {
            int num_dof = Nodle(j);
	    if (num_dof >= 0)
	      for (int m = 0; m < nb_u; m++)
		Eloc(m)(j) = B(num_dof+m*Nvol);
          }
	
        for (int m = 0; m < nb_u; m++)
          ModifyLocalComponentVector(mesh_num, Eloc(m), i);        
      }
  }


  //! Retrieves local components of the solution on the given element
  template<class Dimension, class T>
  void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension>& mesh_num,
			const Vector<T>& B, int i,
			Vector<Vector<T> >& Eloc) const
  {
    int N = mesh_num.GetNbLocalDof(i);
    const IVect& Nodle = mesh_num.Element(i).GetNodle();
    int nb_unknowns = Eloc.GetM();
    for (int m = 0; m < nb_unknowns; m++)
      {
	Eloc(m).Reallocate(N);
	Eloc(m).Zero();
      }
    
    int Nvol = mesh_num.GetNbDof();
    int nb_u = min(int(B.GetM()/Nvol), nb_unknowns);
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      {
        int offset = Nodle(0);
        for (int j = 0; j < N; j++)
          {
            for (int m = 0; m < nb_u; m++)
              Eloc(m)(j) = B(offset + j + m*Nvol);
          }
      }
    else
      {        
        for (int j = 0; j < N; j++)
          {
            int num_dof = Nodle(j);
            if (num_dof >= 0)
	      for (int m = 0; m < nb_u; m++)
		Eloc(m)(j) = B(num_dof+m*Nvol);
          }
        
	ModifyLocalComponentVector(mesh_num, Eloc, i);        
      }
  }
  
  
  //! Adds a local vector to the global solution
  template<class Dimension, class T, int nb_unknowns>
  void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension>& mesh_num,
			const T& alpha, TinyVector<Vector<T>, nb_unknowns>& Eloc,
			int i, Vector<T>& C) const
  {
    int N = mesh_num.GetNbLocalDof(i);
    const IVect& Nodle = mesh_num.Element(i).GetNodle();
    int Nvol = mesh_num.GetNbDof();
    int nb_u = min(int(C.GetM()/Nvol), nb_unknowns);
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      {
        int offset = Nodle(0);
        for (int j = 0; j < N; j++)
          {
            for (int m = 0; m < nb_u; m++)
              C(offset + j + m*Nvol) += alpha*Eloc(m)(j);
          }
      }
    else
      {
        for (int m = 0; m < nb_u; m++)
          ModifyLocalUnknownVector(mesh_num, Eloc(m), i);
        
        for (int j = 0; j < N; j++)
          {
            int num_dof = Nodle(j);
            if (num_dof >= 0)
	      for (int m = 0; m < nb_u; m++)
		C(num_dof+m*Nvol) += alpha*Eloc(m)(j);
          }
      }
  }
  
  
  //! Modification of a local vector so that it can be directly added to the global one
  template<class Dimension, class T>
  void NumberMap::
  ModifyLocalUnknownVector(const MeshNumbering<Dimension>& mesh_num,
			   Vector<Vector<T> >& U, int iquad) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;

    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    const ElementNumbering& elt = mesh_num.Element(iquad);
    const typename Dimension::Element_Rn& elt_ = mesh.Element(iquad);
    const IVect& NegativeNumDof = elt.GetNegativeDofNumber();
    
    Real_wp one(1);
    for (int m = 0; m < U.GetM(); m++)
      for (int k = 0; k < NegativeNumDof.GetM(); k++)
	U(m)(NegativeNumDof(k)) *= -one;
   
    for (int num_loc = 0; num_loc < elt_.GetNbBoundary(); num_loc++)
      {
        int r = elt.GetOrderFace(num_loc); 
                
        if (r > 0)
          {
            int rot = elt_.GetOrientationFace(num_loc);
            int offset = elt.GetOffsetFace(num_loc);
            int nf = elt_.numBoundary(num_loc);
	    for (int m = 0; m < U.GetM(); m++)
	      this->ProjectToGlobalDofs(U(m), rot, r, offset, mesh.Boundary(nf));
          }
      }
  }


  //! Modification of a local vector so that it can be directly added to the global one
  template<class Dimension, class T>
  void NumberMap::
  ModifyLocalUnknownVector(const MeshNumbering<Dimension>& mesh_num,
			   Vector<T>& U, int iquad) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;

    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    const ElementNumbering& elt = mesh_num.Element(iquad);
    const typename Dimension::Element_Rn& elt_ = mesh.Element(iquad);
    const IVect& NegativeNumDof = elt.GetNegativeDofNumber();
    
    Real_wp one(1);
    for (int k = 0; k < NegativeNumDof.GetM(); k++)
      U(NegativeNumDof(k)) *= -one;
      
    for (int num_loc = 0; num_loc < elt_.GetNbBoundary(); num_loc++)
      {
        int r = elt.GetOrderFace(num_loc); 
                
        if (r > 0)
          {
            int rot = elt_.GetOrientationFace(num_loc);
            int offset = elt.GetOffsetFace(num_loc);
            int nf = elt_.numBoundary(num_loc);
	    this->ProjectToGlobalDofs(U, rot, r, offset, mesh.Boundary(nf));
          }
      }
  }


  //! Modification of a local vector so that it can be directly added to the global one
  template<class Dimension, class T, int nb_unknowns>
  void NumberMap::
  ModifyLocalUnknownVector(const MeshNumbering<Dimension>& mesh_num,
			   Vector<TinyVector<T, nb_unknowns> >& U, int iquad) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;

    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    const ElementNumbering& elt = mesh_num.Element(iquad);
    const typename Dimension::Element_Rn& elt_ = mesh.Element(iquad);
    const IVect& NegativeNumDof = elt.GetNegativeDofNumber();
    
    Real_wp one(1);
    for (int k = 0; k < NegativeNumDof.GetM(); k++)
      U(NegativeNumDof(k)) *= -one;
   
    for (int num_loc = 0; num_loc < elt_.GetNbBoundary(); num_loc++)
      {
        int r = elt.GetOrderFace(num_loc); 
                
        if (r > 0)
          {
            int rot = elt_.GetOrientationFace(num_loc);
            int offset = elt.GetOffsetFace(num_loc);
            int nf = elt_.numBoundary(num_loc);
	    this->ProjectToGlobalDofs(U, rot, r, offset, mesh.Boundary(nf));
          }
      }
  }
  
  
  //! Modification of a local vector so that it can be directly added to the global one
  template<class Dimension, class T>
  void NumberMap::
  ModifyLocalComponentVector(const MeshNumbering<Dimension>& mesh_num,
			     Vector<Vector<T> >& U, int iquad) const
  {
    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    const ElementNumbering& elt = mesh_num.Element(iquad);
    const typename Dimension::Element_Rn& elt_ = mesh.Element(iquad);

    const IVect& NegativeNumDof = elt.GetNegativeDofNumber();

    for (int m = 0; m < U.GetM(); m++)    
      for (int k = 0; k < NegativeNumDof.GetM(); k++)
	U(m)(NegativeNumDof(k)) *= -1.0;

    for (int num_loc = 0; num_loc < elt_.GetNbBoundary(); num_loc++)
      {
        int r = elt.GetOrderFace(num_loc); 
        
        if (r > 0)
          {
            int rot = elt_.GetOrientationFace(num_loc);
            int offset = elt.GetOffsetFace(num_loc);
            int nf = elt_.numBoundary(num_loc);
	    for (int m = 0; m < U.GetM(); m++)
	      this->TransposeProjectToGlobalDofs(U(m), rot, r, offset, mesh.Boundary(nf));
          }
      }
  }


  //! Modification of a local vector so that it can be directly added to the global one
  template<class Dimension, class T>
  void NumberMap::
  ModifyLocalComponentVector(const MeshNumbering<Dimension>& mesh_num,
			     Vector<T>& U, int iquad) const
  {
    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    const ElementNumbering& elt = mesh_num.Element(iquad);
    const typename Dimension::Element_Rn& elt_ = mesh.Element(iquad);

    const IVect& NegativeNumDof = elt.GetNegativeDofNumber();

    for (int k = 0; k < NegativeNumDof.GetM(); k++)
      U(NegativeNumDof(k)) *= -1.0;
    
    for (int num_loc = 0; num_loc < elt_.GetNbBoundary(); num_loc++)
      {
        int r = elt.GetOrderFace(num_loc); 
        
        if (r > 0)
          {
            int rot = elt_.GetOrientationFace(num_loc);
            int offset = elt.GetOffsetFace(num_loc);
            int nf = elt_.numBoundary(num_loc);
	    this->TransposeProjectToGlobalDofs(U, rot, r, offset, mesh.Boundary(nf));
          }
      }
  }
  
  
  //! projection from a global vector to a local vector
  template<class Dimension, class T>
  void NumberMap::GetGlobalUnknownVector(const MeshNumbering<Dimension>& mesh_num,
					 Vector<Vector<T> >& U, int iquad) const
  {
    if (dg_element)
      return;

    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    const ElementNumbering& elt = mesh_num.Element(iquad);
    const typename Dimension::Element_Rn& elt_ = mesh.Element(iquad);
    const IVect& NegativeNumDof = elt.GetNegativeDofNumber();
    
    for (int m = 0; m < U.GetM(); m++)    
      for (int k = 0; k < NegativeNumDof.GetM(); k++)
	U(m)(NegativeNumDof(k)) *= -1.0;

    for (int num_loc = 0; num_loc < elt_.GetNbBoundary(); num_loc++)
      {
        int r = elt.GetOrderFace(num_loc); 
        
        if (r > 0)
          {
            int rot = elt_.GetOrientationFace(num_loc);
            int offset = elt.GetOffsetFace(num_loc);
            int nf = elt_.numBoundary(num_loc);
            for (int m = 0; m < U.GetM(); m++)    
	      this->ProjectToLocalDofs(U(m), rot, r, offset, mesh.Boundary(nf));
          }
      }
  }
  

  //! Modifies rows of the local elementary matrix
  //! so that it can be directly added to the global matrix
  template<class Dimension, class T, class Prop, class Storage>
  void NumberMap::ModifyLocalDenseRowMatrix(const MeshNumbering<Dimension>& mesh_num,
                                            Matrix<T, Prop, Storage>& mat,
                                            int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    bool sym = IsSymmetricMatrix(mat);
    int n = mat.GetN();
    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    const ElementNumbering& elt = mesh_num.Element(iquad);
    const typename Dimension::Element_Rn& elt_ = mesh.Element(iquad);
    const IVect& NegativeNumDof = elt.GetNegativeDofNumber();
    int nb_dof = elt.GetNbDof();
    if (nb_u == 0)
      nb_u = mat.GetM()/nb_dof;

    if (sym)
      for (int k = 0; k < NegativeNumDof.GetM(); k++)
        for (int m = 0; m < nb_u; m++)
          {
            int i = NegativeNumDof(k) + m*nb_dof;

            for (int j = off_row + i; j < n; j++)
              mat(off_row + i, j) *= -1.0;
            
          }
    else
      for (int k = 0; k < NegativeNumDof.GetM(); k++)
        for (int m = 0; m < nb_u; m++)
          {
            int i = NegativeNumDof(k) + m*nb_dof;
            for (int j = 0; j < n; j++)
              mat(off_row + i, j) *= -1.0;
          }

    for (int num_loc = 0; num_loc < elt_.GetNbBoundary(); num_loc++)
      {
        int r = elt.GetOrderFace(num_loc); 
        
        if (r > 0)
          {
            int rot = elt_.GetOrientationFace(num_loc);
            int offset = elt.GetOffsetFace(num_loc);
            int nf = elt_.numBoundary(num_loc);
            this->ModifyRowToGlobalDofs(mat, rot, r, off_row + offset, mesh.Boundary(nf));
          }
      }
  }
  

  //! Modifies columns of the local elementary matrix
  //! so that it can be directly added to the global matrix
  template<class Dimension, class T, class Prop, class Storage>
  void NumberMap::ModifyLocalDenseColumnMatrix(const MeshNumbering<Dimension>& mesh_num,
                                               Matrix<T, Prop, Storage>& mat, int iquad, int nb_u, int off_col) const
  {
    bool sym = IsSymmetricMatrix(mat);
    int n = mat.GetM();
    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    const ElementNumbering& elt = mesh_num.Element(iquad);
    const typename Dimension::Element_Rn& elt_ = mesh.Element(iquad);
    const IVect& NegativeNumDof = elt.GetNegativeDofNumber();
    int nb_dof = elt.GetNbDof();
    if (nb_u == 0)
      nb_u = mat.GetN()/nb_dof;
    
    if (sym)
      for (int k = 0; k < NegativeNumDof.GetM(); k++)
	for (int m = 0; m < nb_u; m++)
          {
            int j = NegativeNumDof(k) + m*nb_dof;
            for (int i = 0; i <= off_col + j; i++)
              mat(i, off_col + j) *= -1.0;
          }
    else
      for (int k = 0; k < NegativeNumDof.GetM(); k++)
        for (int m = 0; m < nb_u; m++)
          {
            int j = NegativeNumDof(k) + m*nb_dof;
            for (int i = 0; i < n; i++)
              mat(i, off_col + j) *= -1.0;
          }

    for (int num_loc = 0; num_loc < elt_.GetNbBoundary(); num_loc++)
      {
        int r = elt.GetOrderFace(num_loc); 
        
        if (r > 0)
          {
            int rot = elt_.GetOrientationFace(num_loc);
            int offset = elt.GetOffsetFace(num_loc);
            int nf = elt_.numBoundary(num_loc);
            this->ModifyColumnToGlobalDofs(mat, rot, r, off_col + offset, mesh.Boundary(nf));
          }
      }
  }


  //! Modifies rows of the local elementary matrix
  //! so that it can be directly added to the global matrix
  template<class Dimension, class T, class Prop, class Storage>
  void NumberMap::ModifyLocalSparseRowMatrix(const MeshNumbering<Dimension>& mesh_num,
                                             Matrix<T, Prop, Storage>& mat,
                                             int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    //int n = mat.GetN();
    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    const ElementNumbering& elt = mesh_num.Element(iquad);
    const typename Dimension::Element_Rn& elt_ = mesh.Element(iquad);
    const IVect& NegativeNumDof = elt.GetNegativeDofNumber();
    int nb_dof = elt.GetNbDof();
    if (nb_u == 0)
      nb_u = mat.GetM()/nb_dof;

    for (int k = 0; k < NegativeNumDof.GetM(); k++)
      for (int m = 0; m < nb_u; m++)
        {
          int i = NegativeNumDof(k) + m*nb_dof;
          int size_row = mat.GetRowSize(off_row+i);
          for (int j = 0; j < size_row; j++)
            mat.Value(off_row + i, j) *= -1.0;          
        }
    
    for (int num_loc = 0; num_loc < elt_.GetNbBoundary(); num_loc++)
      {
        int r = elt.GetOrderFace(num_loc); 
        
        if (r > 0)
          {
            cout << "Not implemented" << endl;
            abort();
          }
      }
  }
  

  //! Modifies columns of the local elementary matrix
  //! so that it can be directly added to the global matrix
  template<class Dimension, class T, class Prop, class Storage>
  void NumberMap::ModifyLocalSparseColumnMatrix(const MeshNumbering<Dimension>& mesh_num,
                                                Matrix<T, Prop, Storage>& mat, int iquad, int nb_u, int off_col) const
  {
    int n = mat.GetM();
    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    const ElementNumbering& elt = mesh_num.Element(iquad);
    const typename Dimension::Element_Rn& elt_ = mesh.Element(iquad);
    const IVect& NegativeNumDof = elt.GetNegativeDofNumber();
    int nb_dof = elt.GetNbDof();
    if (nb_u == 0)
      nb_u = mat.GetN()/nb_dof;

    if (NegativeNumDof.GetM() > 0)
      {
        Vector<bool> is_col_neg(mat.GetN());
        is_col_neg.Fill(false);    
        for (int k = 0; k < NegativeNumDof.GetM(); k++)
          for (int m = 0; m < nb_u; m++)
            {
              int j = NegativeNumDof(k) + m*nb_dof;
              is_col_neg(off_col+j) = true;
            }

        for (int i = 0; i < n; i++)
          for (int j0 = 0; j0 < mat.GetRowSize(i); j0++)
            if (is_col_neg(mat.Index(i, j0)))
              mat.Value(i, j0) *= -1.0;
      }        

    for (int num_loc = 0; num_loc < elt_.GetNbBoundary(); num_loc++)
      {
        int r = elt.GetOrderFace(num_loc); 
        
        if (r > 0)
          {
            cout << "Not implemented" << endl;
            abort();
          }
      }
  }


  //! Modifies rows of the local elementary matrix
  //! so that it can be directly added to the global matrix
  template<class Dimension, class T>
  void NumberMap::ModifyLocalRowMatrix(const MeshNumbering<Dimension>& mesh_num,
                                       VirtualMatrix<T>& mat,
                                       int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;

    // we search the true type of the matrix
    if (mat.IsSymmetric())
      {
        try
          {
            Matrix<T, Symmetric, RowSymPacked>& A = dynamic_cast<Matrix<T, Symmetric, RowSymPacked>& >(mat);
            ModifyLocalDenseRowMatrix(mesh_num, A, iquad, nb_u, off_row);
          }
        catch (const std::bad_cast&)
          {
            try
              {
                Matrix<T, Symmetric, ArrayRowSymSparse>& A = dynamic_cast<Matrix<T, Symmetric, ArrayRowSymSparse>& >(mat);
                ModifyLocalSparseRowMatrix(mesh_num, A, iquad, nb_u, off_row);
              }
            catch (const std::bad_cast&)
              {
                try
                  {
                    //Matrix<T, Symmetric, DiagonalRow>& A = dynamic_cast<Matrix<T, Symmetric, DiagonalRow>& >(mat);
                    dynamic_cast<Matrix<T, Symmetric, DiagonalRow>& >(mat);
                    // nothing to do
                  }
                catch (const std::bad_cast&)
                  {
                    cout << "Not implemented " << endl;
                    abort();
                  }
              }
          }
      }
    else
      {
        try
          {
            Matrix<T>& A = dynamic_cast<Matrix<T>& >(mat);
            ModifyLocalDenseRowMatrix(mesh_num, A, iquad, nb_u, off_row);
          }
        catch (const std::bad_cast&)
          {
            try
              {
                Matrix<T, General, ArrayRowSparse>& A = dynamic_cast<Matrix<T, General, ArrayRowSparse>& >(mat);
                ModifyLocalSparseRowMatrix(mesh_num, A, iquad, nb_u, off_row);
              }
            catch (const std::bad_cast&)
              {
                cout << "Not implemented " << endl;
                abort();
              }
          }
      }
  }


  //! Modifies columns of the local elementary matrix
  //! so that it can be directly added to the global matrix
  template<class Dimension, class T>
  void NumberMap::ModifyLocalColumnMatrix(const MeshNumbering<Dimension>& mesh_num,
                                          VirtualMatrix<T>& mat,
                                          int iquad, int nb_u, int off_col) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;

    // we search the true type of the matrix
    if (mat.IsSymmetric())
      {
        try
          {
            Matrix<T, Symmetric, RowSymPacked>& A = dynamic_cast<Matrix<T, Symmetric, RowSymPacked>& >(mat);
            ModifyLocalDenseColumnMatrix(mesh_num, A, iquad, nb_u, off_col);
          }
        catch (const std::bad_cast&)
          {
            try
              {
                Matrix<T, Symmetric, ArrayRowSymSparse>& A = dynamic_cast<Matrix<T, Symmetric, ArrayRowSymSparse>& >(mat);
                ModifyLocalSparseColumnMatrix(mesh_num, A, iquad, nb_u, off_col);
              }
            catch (const std::bad_cast&)
              {
                try
                  {
                    //Matrix<T, Symmetric, DiagonalRow>& A = dynamic_cast<Matrix<T, Symmetric, DiagonalRow>& >(mat);
                    dynamic_cast<Matrix<T, Symmetric, DiagonalRow>& >(mat);
                    // nothing to do
                  }
                catch (const std::bad_cast&)
                  {
                    cout << "Not implemented " << endl;
                    abort();
                  }
              }
          }
      }
    else
      {
        try
          {
            Matrix<T>& A = dynamic_cast<Matrix<T>& >(mat);
            ModifyLocalDenseColumnMatrix(mesh_num, A, iquad, nb_u, off_col);
          }
        catch (const std::bad_cast&)
          {
            try
              {
                Matrix<T, General, ArrayRowSparse>& A = dynamic_cast<Matrix<T, General, ArrayRowSparse>& >(mat);
                ModifyLocalSparseColumnMatrix(mesh_num, A, iquad, nb_u, off_col);
              }
            catch (const std::bad_cast&)
              {
                cout << "Not implemented " << endl;
                abort();
              }
          }
      }
  }

  
  template<>
  void NumberMap::ModifyLocalRowMatrix(const MeshNumbering<Dimension2>& mesh_num,
                                       VirtualMatrix<TinyVector<Real_wp, 1> >& mat,
                                       int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    try
      {
        Matrix<TinyVector<Real_wp, 1> >& A = dynamic_cast<Matrix<TinyVector<Real_wp, 1> >& >(mat);
        ModifyLocalDenseRowMatrix(mesh_num, A, iquad, nb_u, off_row);
      }
    catch (const std::bad_cast&)
      {
        cout << "Not implemented " << endl;
        abort();
      }    
  }


  template<>
  void NumberMap::ModifyLocalRowMatrix(const MeshNumbering<Dimension2>& mesh_num,
                                       VirtualMatrix<TinyVector<Real_wp, 2> >& mat,
                                       int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    try
      {
        Matrix<TinyVector<Real_wp, 2> >& A = dynamic_cast<Matrix<TinyVector<Real_wp, 2> >& >(mat);
        ModifyLocalDenseRowMatrix(mesh_num, A, iquad, nb_u, off_row);
      }
    catch (const std::bad_cast&)
      {
        cout << "Not implemented " << endl;
        abort();
      }    
  }

  
  template<>
  void NumberMap::ModifyLocalRowMatrix(const MeshNumbering<Dimension3>& mesh_num,
                                       VirtualMatrix<TinyVector<Real_wp, 1> >& mat,
                                       int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    try
      {
        Matrix<TinyVector<Real_wp, 1> >& A = dynamic_cast<Matrix<TinyVector<Real_wp, 1> >& >(mat);
        ModifyLocalDenseRowMatrix(mesh_num, A, iquad, nb_u, off_row);
      }
    catch (const std::bad_cast&)
      {
        cout << "Not implemented " << endl;
        abort();
      }    
  }


  template<>
  void NumberMap::ModifyLocalRowMatrix(const MeshNumbering<Dimension3>& mesh_num,
                                       VirtualMatrix<TinyVector<Real_wp, 2> >& mat,
                                       int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    try
      {
        Matrix<TinyVector<Real_wp, 2> >& A = dynamic_cast<Matrix<TinyVector<Real_wp, 2> >& >(mat);
        ModifyLocalDenseRowMatrix(mesh_num, A, iquad, nb_u, off_row);
      }
    catch (const std::bad_cast&)
      {
        cout << "Not implemented " << endl;
        abort();
      }    
  }  

  
  template<>
  void NumberMap::ModifyLocalRowMatrix(const MeshNumbering<Dimension3>& mesh_num,
                                       VirtualMatrix<TinyVector<Real_wp, 3> >& mat,
                                       int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    try
      {
        Matrix<TinyVector<Real_wp, 3> >& A = dynamic_cast<Matrix<TinyVector<Real_wp, 3> >& >(mat);
        ModifyLocalDenseRowMatrix(mesh_num, A, iquad, nb_u, off_row);
      }
    catch (const std::bad_cast&)
      {
        cout << "Not implemented " << endl;
        abort();
      }    
  }


  template<>
  void NumberMap::ModifyLocalColumnMatrix(const MeshNumbering<Dimension2>& mesh_num,
                                          VirtualMatrix<TinyVector<Real_wp, 1> >& mat,
                                          int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    try
      {
        Matrix<TinyVector<Real_wp, 1> >& A = dynamic_cast<Matrix<TinyVector<Real_wp, 1> >& >(mat);
        ModifyLocalDenseColumnMatrix(mesh_num, A, iquad, nb_u, off_row);
      }
    catch (const std::bad_cast&)
      {
        cout << "Not implemented " << endl;
        abort();
      }    
  }


  template<>
  void NumberMap::ModifyLocalColumnMatrix(const MeshNumbering<Dimension2>& mesh_num,
                                          VirtualMatrix<TinyVector<Real_wp, 2> >& mat,
                                          int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    try
      {
        Matrix<TinyVector<Real_wp, 2> >& A = dynamic_cast<Matrix<TinyVector<Real_wp, 2> >& >(mat);
        ModifyLocalDenseColumnMatrix(mesh_num, A, iquad, nb_u, off_row);
      }
    catch (const std::bad_cast&)
      {
        cout << "Not implemented " << endl;
        abort();
      }    
  }
  
  
  template<>
  void NumberMap::ModifyLocalColumnMatrix(const MeshNumbering<Dimension3>& mesh_num,
                                       VirtualMatrix<TinyVector<Real_wp, 1> >& mat,
                                       int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    try
      {
        Matrix<TinyVector<Real_wp, 1> >& A = dynamic_cast<Matrix<TinyVector<Real_wp, 1> >& >(mat);
        ModifyLocalDenseColumnMatrix(mesh_num, A, iquad, nb_u, off_row);
      }
    catch (const std::bad_cast&)
      {
        cout << "Not implemented " << endl;
        abort();
      }    
  }


  template<>
  void NumberMap::ModifyLocalColumnMatrix(const MeshNumbering<Dimension3>& mesh_num,
                                          VirtualMatrix<TinyVector<Real_wp, 2> >& mat,
                                          int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    try
      {
        Matrix<TinyVector<Real_wp, 2> >& A = dynamic_cast<Matrix<TinyVector<Real_wp, 2> >& >(mat);
        ModifyLocalDenseColumnMatrix(mesh_num, A, iquad, nb_u, off_row);
      }
    catch (const std::bad_cast&)
      {
        cout << "Not implemented " << endl;
        abort();
      }    
  }
  

  template<>
  void NumberMap::ModifyLocalColumnMatrix(const MeshNumbering<Dimension3>& mesh_num,
                                          VirtualMatrix<TinyVector<Real_wp, 3> >& mat,
                                          int iquad, int nb_u, int off_row) const
  {
    if (dg_element == ElementReference_Base::DISCONTINUOUS)
      return;
    
    try
      {
        Matrix<TinyVector<Real_wp, 3> >& A = dynamic_cast<Matrix<TinyVector<Real_wp, 3> >& >(mat);
        ModifyLocalDenseColumnMatrix(mesh_num, A, iquad, nb_u, off_row);
      }
    catch (const std::bad_cast&)
      {
        cout << "Not implemented " << endl;
        abort();
      }    
  }

  
  //! clear the details about how to number a mesh
  void NumberMap::Clear()
  {
    nb_dof_edge.Clear();
    nb_dof_quad.Clear(); nb_dof_tri.Clear(); nb_dof_pyramid.Clear();
    nb_dof_hexa.Clear(); nb_dof_tetra.Clear(); nb_dof_wedge.Clear();
    FacesDof_Rotation_Quad.Clear(); FacesDof_Rotation_Tri.Clear();
    FacesQuad_Rotation_Quad.Clear(); FacesQuad_Rotation_Tri.Clear();
    EdgesDof_Symmetry.Clear(); EdgesQuad_Rotation.Clear();
    only_hexahedral = false; dg_element = false;
  }
  
  
  //! Node numbers for triangles in Gmsh
  void ConstructGmshTriangularNumbering(int r, Matrix<int>& NumNodes)
  {
    int r1 = 0, r2 = r;
    NumNodes.Reallocate(r+1, r+1);
    NumNodes.Fill(-1);
    int node = 0;
    // recursive numbering
    while (r1 <= r2)
      {
	// three vertices
	NumNodes(r1, r1) = node++;
	
	if (r1 < r2)
	  {
	    NumNodes(r2, r1) = node++;
	    NumNodes(r1, r2) = node++;
	  }
	
	// three edges
	for (int j = r1+1; j < r2; j++)	  
	  NumNodes(j, r1) = node++;
	
	for (int j = r1+1; j < r2; j++)
	  NumNodes(r1+r2-j, j) = node++;
	
	for (int j = r1+1; j < r2; j++)
	  NumNodes(r1, r1+r2-j) = node++;
	
	// next sub-triangle
	r1++; r2 -= 2;
      }    
  }  
  
  
  //! Equivalence between Gmsh numbers and Montjoie numbers
  void NumberMap::GetGmshTriangularNumbering(int r, IVect& num)
  {
    Matrix<int> coor, NumNodes2D, NumGmsh;
    MeshNumbering<Dimension2>::ConstructTriangularNumbering(r, NumNodes2D, coor);
    ConstructGmshTriangularNumbering(r, NumGmsh);
    num.Reallocate((r+1)*(r+2)/2);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	num(NumGmsh(i, j)) = NumNodes2D(i, j);
    
  }
  
  
  //! Node numbers for quadrilaterals in Gmsh
  void ConstructGmshQuadrilateralNumbering(int r, Matrix<int>& NumNodes)
  {
    int r1 = 0, r2 = r;
    NumNodes.Reallocate(r+1, r+1);
    NumNodes.Fill(-1);
    int node = 0;
    // recursive numbering
    while (r1 <= r2)
      {
	// four vertices
	NumNodes(r1, r1) = node++;
	
	if (r1 < r2)
	  {
	    NumNodes(r2, r1) = node++;
	    NumNodes(r2, r2) = node++;
	    NumNodes(r1, r2) = node++;
	  }
	
	// four edges
	for (int j = r1+1; j < r2; j++)	  
	  NumNodes(j, r1) = node++;

	for (int j = r1+1; j < r2; j++)	  
	  NumNodes(r2, j) = node++;
	
	for (int j = r1+1; j < r2; j++)
	  NumNodes(r1+r2-j, r2) = node++;
	
	for (int j = r1+1; j < r2; j++)
	  NumNodes(r1, r1+r2-j) = node++;
	
	// next sub-quadrilateral
	r1++; r2--;
      }
  }  

  
  //! Equivalence between Gmsh numbers and Montjoie numbers
  void NumberMap::GetGmshQuadrilateralNumbering(int r, IVect& num)
  {
    Matrix<int> coor, NumNodes2D, NumGmsh;
    MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(r, NumNodes2D, coor);
    ConstructGmshQuadrilateralNumbering(r, NumGmsh);
    num.Reallocate((r+1)*(r+1));
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	num(NumGmsh(i, j)) = NumNodes2D(i, j);
  }
  

  //! returns the global numbers of dofs of the element
  template<class Dimension>
  IVect NumberMap
  ::GetDofNumberOnElement(const MeshNumbering<Dimension>& mesh, int i) const
  {
    if (dg_element)
      {
	IVect Nodle(mesh.GetNbLocalDof(i));
	Nodle.Fill();
        int offset = mesh.Element(i).GetNumberDof(0);
	for (int i = 0; i < Nodle.GetM(); i++)
	  Nodle(i) += offset;
	
	return Nodle;
      }
    
    return mesh.Element(i).GetNodle();
  }
  
  
  //! returns the Gmsh entity number for an edge
  template<class Dimension>
  int NumberMap::GetGmshEntityNumber(const Edge<Dimension>& e, int order)
  {
    switch (order)
      {
      case 1 :
	return 1;
      case 2 :
	return 8;
      case 3 :
	return 26;
      case 4 :
	return 27;
      case 5 :
	return 28;
      case 6 :
	return 62;
      case 7 :
	return 63;
      case 8 :
	return 64;
      case 9 :
	return 65;
      case 10 :
	return 66;
      default :
	abort();
      }
    
    return -1;
  }

  
  //! returns the Gmsh entity number for a face (triangle or quadrangle)
  template<class Dimension>
  int NumberMap::GetGmshEntityNumber(const Face<Dimension>& e, int order)
  {
    if (e.GetNbVertices() == 3)
      {
	switch (order)
	  {
	  case 1 :
	    return 2;
	  case 2 :
	    return 9;
	  case 3 :
	    return 21;
	  case 4 :
	    return 23;
	  case 5 :
	    return 25;
	  case 6 :
	    return 42;
	  case 7 :
	    return 43;
	  case 8 :
	    return 44;
	  case 9 :
	    return 45;
	  case 10 :
	    return 46;
	  default :
	    abort();
	  }
      }
    else if (e.GetNbVertices() == 4)
      {
	switch (order)
	  {
	  case 1 :
	    return 3;
	  case 2 :
	    return 10;
	  case 3 :
            return 36;
          case 4 :
            return 37;
          case 5 :
            return 38;
	  case 6 :
	    return 47;
	  case 7 :
	    return 48;
	  case 8 :
	    return 49;
	  case 9 :
	    return 50;
	  case 10 :
	    return 51;
          default :
	    abort();
	  }
      }
    
    return -1;
  }

  
  //! returns the Gmsh entity number for a volume (tetra, pyramid, prism or hexa)
  int NumberMap::GetGmshEntityNumber(const Volume& e, int order)
  {
    if (e.GetNbVertices() == 4)
      {
	switch (order)
	  {
	  case 1 :
	    return 4;
	  case 2 :
	    return 11;
	  case 3 :
	    return 29;
	  case 4 :
	    return 30;
	  case 5 :
	    return 31;
	  case 6 :
	    return 71;
	  case 7 :
	    return 72;
	  case 8 :
	    return 73;
	  case 9 :
	    return 74;
	  case 10 :
	    return 75;
	  default :
	    abort();
	  }
      }
    else if (e.GetNbVertices() == 5)
      {
	switch (order)
	  {
	  case 1 :
	    return 7;
	  case 2 :
	    return 14;
	  case 3 :
	    return 109;
	  case 4 :
	    return 110;
	  case 5 :
	    return 111;
	  default :
	    abort();
	  }
      }
    else if (e.GetNbVertices() == 6)
      {
	switch (order)
	  {
	  case 1 :
	    return 6;
	  case 2 :
	    return 13;
	  case 3 :
	    return 90;
	  case 4 :
	    return 91;
	  case 5 :
	    return 112;
	  default :
	    abort();
	  }
      }
    else if (e.GetNbVertices() == 8)
      {
	switch (order)
	  {
	  case 1 :
	    return 5;
	  case 2 :
	    return 12;
	  case 3 :
	    return 92;
	  case 4 :
	    return 93;
	  case 5 :
	    return 94;
	  case 6 :
	    return 95;
	  case 7 :
	    return 96;
	  case 8 :
	    return 97;
	  case 9 :
	    return 98;
	  default :
	    abort();
	  }
      }
    
    return -1;
  }
  
  
  //! displays details about class NumberMap
  ostream& operator <<(ostream& out, const NumberMap& nmap)
  {
    out<<nmap.nb_dof_vertex<<" dof by vertex"<<endl;
    out<<nmap.nb_dof_edge<<" dof by edge"<<endl;
    out<<nmap.nb_dof_quad<<" dof by quadrilateral"<<endl;
    out<<nmap.nb_dof_tri<<" dof by triangle"<<endl;
    out<<nmap.nb_dof_hexa<<" dof by hexahedron"<<endl;
    out<<nmap.nb_dof_tetra<<" dof by tetrahedron"<<endl;
    return out;
  }

  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! Node numbers for tetrahedra in Gmsh
  void ConstructGmshTetrahedralNumbering(int r, Array3D<int>& NumNodes)
  {
    Matrix<int> NumTri;
    NumNodes.Reallocate(r+1, r+1, r+1);
    NumNodes.Fill(-1);

    int r1 = 0, r2 = r;    
    int node = 0;
    while ( r1 <= r2)
      {
	// four vertices
	NumNodes(r1, r1, r1) = node++;
	
	if (r1 < r2)
	  {
	    NumNodes(r2, r1, r1) = node++;
	    NumNodes(r1, r2, r1) = node++;
	    NumNodes(r1, r1, r2) = node++;
	  }
	
	// nodes associated with edges
	// first edge 0 - 1 
	for (int i = r1+1; i < r2; i++)
	  NumNodes(i, r1, r1) = node++;
	
	// second edge 1 - 2
	for (int i = r1+1; i < r2; i++)
	  NumNodes(r1+r2-i, i, r1) = node++;
	
	// third edge 2 - 0
	for (int i = r1+1; i < r2; i++)
	  NumNodes(r1, r1+r2-i, r1) = node++;
	
	// fourth edge 3 - 0
	for (int i = r1+1; i < r2; i++)
	  NumNodes(r1, r1, r1+r2-i) = node++;
	
	// fifth edge 3 - 2
	for (int i = r1+1; i < r2; i++)
	  NumNodes(r1, i, r1+r2-i) = node++;
	
	// sixth edge = 3 - 1
	for (int i = r1+1; i < r2; i++)
	  NumNodes(i, r1, r1+r2-i) = node++;
	
	// nodes associated with faces
	if (r2 > r1+2)
	  {
	    ConstructGmshTriangularNumbering(r2-r1-3, NumTri);
	    	    
	    // first face = 0-1-2
	    for (int i = r1+1; i < r2; i++)
	      for (int j = r1+1; j < r2+r1-i; j++)
		NumNodes(j, i, r1) = node + NumTri(i-r1-1, j-r1-1);
			    
	    node += (r2-r1-2)*(r2-r1-1)/2;

	    // second face = 0-1-3
	    for (int i = r1+1; i < r2; i++)
	      for (int j = r1+1; j < r2+r1-i; j++)
		NumNodes(i, r1, j) = node + NumTri(i-r1-1, j-r1-1);
	    
	    node += (r2-r1-2)*(r2-r1-1)/2;	    

	    // third face = 0-2-3
	    for (int i = r1+1; i < r2; i++)
	      for (int j = r1+1; j < r2+r1-i; j++)
		NumNodes(r1, j, i) = node + NumTri(i-r1-1, j-r1-1);
	    
	    node += (r2-r1-2)*(r2-r1-1)/2;	    
	    
	    // fourth face = 1-2-3
	    for (int i = r1+1; i < r2; i++)
	      for (int j = r1+1; j < r2+r1-i; j++)
		NumNodes(i, j, 2*r1+r2-i-j) = node + NumTri(i-r1-1, j-r1-1);
	    
	    node += (r2-r1-2)*(r2-r1-1)/2;
	  }
	
	// next sub-tetrahedral
	r1++; r2 -= 3;
      }
    
  }
  

  //! Equivalence between Gmsh numbers and Montjoie numbers
  void NumberMap::GetGmshTetrahedralNumbering(int r, IVect& num)
  {
    Matrix<int> coor; Array3D<int> NumNodes3D, NumGmsh;
    MeshNumbering<Dimension3>::ConstructTetrahedralNumbering(r, NumNodes3D, coor);
    ConstructGmshTetrahedralNumbering(r, NumGmsh);
    num.Reallocate((r+1)*(r+2)*(r+3)/6);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	for (int k = 0; k <= r-i-j; k++)
	  num(NumGmsh(i, j, k)) = NumNodes3D(i, j, k);    
  }
  
  
  //! Equivalence between Gmsh numbers and Montjoie numbers
  void NumberMap::GetGmshPyramidalNumbering(int r, IVect& num)
  {
    Matrix<int> coor; Array3D<int> NumNodes3D;
    MeshNumbering<Dimension3>::ConstructPyramidalNumbering(r, NumNodes3D, coor);
    num.Reallocate((r+1)*(r+2)*(2*r+3)/6); num.Fill();
    // nodes associated with edges
    // first edge 0 - 1 
    int node = 5;
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(i, 0, 0);    
    
    // second edge 0 - 3
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(0, i, 0);    

    // third edge 0 - 4
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(0, 0, i);
    
    // fourth edge 1 - 2
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(r, i, 0);
    
    // fifth edge 1 - 4
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(r-i, 0, i);

    // sixth edge 2 - 3
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(r-i, r, 0);
    
    // seventh edge 2 - 4
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(r-i, r-i, i);
    
    // eigth edge 3 - 4
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(0, r-i, i);
    
    // quadrilateral base
    for (int j = 1; j < r; j++)
      for (int i = 1; i < r; i++)
	num(node++) = NumNodes3D(i, j, 0);
    
    // triangular faces
    for (int j = 1; j < r; j++)
      for (int i = 1; i < r-j; i++)
	num(node++) = NumNodes3D(i, 0, j);

    for (int j = 1; j < r; j++)
      for (int i = 1; i < r-j; i++)
	num(node++) = NumNodes3D(r-j, i, j);

    for (int j = 1; j < r; j++)
      for (int i = 1; i < r-j; i++)
	num(node++) = NumNodes3D(i, r-j, j);

    for (int j = 1; j < r; j++)
      for (int i = 1; i < r-j; i++)
	num(node++) = NumNodes3D(0, i, j);
    
    // interior
    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)
        for (int k = 1; k < r-max(i,j); k++)
          num(node++) = NumNodes3D(i, j, k);
    
  }
  

  //! Equivalence between Gmsh numbers and Montjoie numbers  
  void NumberMap::GetGmshPrismaticNumbering(int r, IVect& num)
  {
    Matrix<int> coor; Array3D<int> NumNodes3D;
    MeshNumbering<Dimension3>::ConstructPrismaticNumbering(r, NumNodes3D, coor);
    num.Reallocate((r+1)*(r+2)*(r+1)/2); num.Fill();
    // nodes associated with edges
    // first edge 0 - 1 
    int node = 6;
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(i, 0, 0);    
    
    // second edge 0 2
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(0, i, 0);    

    // third edge 0 3
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(0, 0, i);

    // fourth edge 1 2
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(r-i, i, 0);    
    
    // fifth edge 1 4
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(r, 0, i);    
    
    // sixth edge 2 5
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(0, r, i);
    
    // seventh edge 3 4
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(i, 0, r);        
    
    // eigth edge 3 5
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(0, i, r);    
    
    // ninth edge 4 5
    for (int i = 1; i < r; i++)
      num(node++) = NumNodes3D(r-i, i, r);    
    
    // quadrilateral faces
    for (int j = 1; j < r; j++)
      for (int i = 1; i < r; i++)
	num(node++) = NumNodes3D(i, 0, j);

    for (int j = 1; j < r; j++)
      for (int i = 1; i < r; i++)
	num(node++) = NumNodes3D(0, i, j);

    for (int j = 1; j < r; j++)
      for (int i = 1; i < r; i++)
	num(node++) = NumNodes3D(r-i, i, j);
    
    // triangular faces
    for (int j = 1; j < r; j++)
      for (int i = 1; i < r-j; i++)
	num(node++) = NumNodes3D(i, j, 0);

    for (int j = 1; j < r; j++)
      for (int i = 1; i < r-j; i++)
	num(node++) = NumNodes3D(i, j, r);

    // interior
    for (int i = 1; i < r; i++)
      for (int j = 1; j < r-i; j++)
        for (int k = 1; k < r; k++)
          num(node++) = NumNodes3D(i, j, k);
  }
  
  
  //! Gmsh node numbers for hexahedra
  void ConstructGmshHexahedralNumbering(int r, Array3D<int>& NumNodes)
  {
    int r1 = 0, r2 = r;
    NumNodes.Reallocate(r+1, r+1, r+1);
    NumNodes.Fill(-1);
    int node = 0;
    Matrix<int> NumQuad;
    // recursive numbering
    while (r1 <= r2)
      {
	// nodes associated with the eight vertices
	NumNodes(r1, r1, r1) = node++;
	
	if (r1 < r2)
	  {
	    NumNodes(r2, r1, r1) = node++;
	    NumNodes(r2, r2, r1) = node++;
	    NumNodes(r1, r2, r1) = node++;
	    NumNodes(r1, r1, r2) = node++;
	    NumNodes(r2, r1, r2) = node++;
	    NumNodes(r2, r2, r2) = node++;
	    NumNodes(r1, r2, r2) = node++;

	    // nodes associated with edges
	    // first edge 0 - 1 
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(i, r1, r1) = node++;    
	    
	    // second edge 0 - 3
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(r1, i, r1) = node++;    
	    
	    // third edge 0 - 4
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(r1, r1, i) = node++;
	    
	    // fourth edge 1 2
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(r2, i, r1) = node++;
	    
	    // fifth edge 1 5
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(r2, r1, i) = node++;    
	    
	    // sixth edge 3 2
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(r1+r2-i, r2, r1) = node++;    
	    
	    // seventh edge 2 6
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(r2, r2, i) = node++;    
	    
	    // eigth edge 3 7
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(r1, r2, i) = node++;
	    
	    // ninth edge 4 5
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(i, r1, r2) = node++;        
	    
	    // tenth edge 4 7
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(r1, i, r2) = node++;
	    
	    // eleventh edge 5 6
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(r2, i, r2) = node++;        
	    
	    // twelfth edge 7 6
	    for (int i = r1+1; i < r2; i++)
	      NumNodes(r1+r2-i, r2, r2) = node++;
	  }
	
	if (r2 > r1+1)
	  {
	    ConstructGmshQuadrilateralNumbering(r2-r1-2, NumQuad);
	    
	    // quadrilateral faces
	    for (int j = r1+1; j < r2; j++)
	      for (int i = r1+1; i < r2; i++)
		NumNodes(j, i, r1) = node + NumQuad(i-r1-1, j-r1-1);
	    
	    node += (r2-r1-1)*(r2-r1-1);
	    
	    for (int j = r1+1; j < r2; j++)
	      for (int i = r1+1; i < r2; i++)
		NumNodes(i, r1, j) = node + NumQuad(i-r1-1, j-r1-1);
	    
	    node += (r2-r1-1)*(r2-r1-1);
	    
	    for (int j = r1+1; j < r2; j++)
	      for (int i = r1+1; i < r2; i++)
		NumNodes(r1, j, i) = node + NumQuad(i-r1-1, j-r1-1);
	    
	    node += (r2-r1-1)*(r2-r1-1);
	    
	    for (int j = r1+1; j < r2; j++)
	      for (int i = r1+1; i < r2; i++)
		NumNodes(r2, i, j) = node + NumQuad(i-r1-1, j-r1-1);
	    
	    node += (r2-r1-1)*(r2-r1-1);
	    
	    for (int j = r1+1; j < r2; j++)
	      for (int i = r1+1; i < r2; i++)
		NumNodes(r1+r2-i, r2, j) = node + NumQuad(i-r1-1, j-r1-1);
	    
	    node += (r2-r1-1)*(r2-r1-1);
	    
	    for (int j = r1+1; j < r2; j++)
	      for (int i = r1+1; i < r2; i++)
		NumNodes(i, j, r2) = node + NumQuad(i-r1-1, j-r1-1);
	    
	    node += (r2-r1-1)*(r2-r1-1);
	  }
	
	// next sub-hexahedral
	r1++; r2--;
      }    
  }
  

  //! Equivalence between Gmsh numbers and Montjoie numbers  
  void NumberMap::GetGmshHexahedralNumbering(int r, IVect& num)
  {
    Matrix<int> coor; Array3D<int> NumNodes3D, NumGmsh;
    MeshNumbering<Dimension3>::ConstructHexahedralNumbering(r, NumNodes3D, coor);
    ConstructGmshHexahedralNumbering(r, NumGmsh);
    num.Reallocate((r+1)*(r+1)*(r+1)); num.Fill();
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	for (int k = 0; k <= r; k++)
	  num(NumGmsh(i, j, k)) = NumNodes3D(i, j, k);    
  }
#endif
  
} // namespace Montjoie

#define MONTJOIE_FILE_NUMBER_MAP_CXX
#endif
