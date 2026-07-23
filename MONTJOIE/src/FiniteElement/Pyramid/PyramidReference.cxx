#ifndef MONTJOIE_FILE_PYRAMID_REFERENCE_CXX

namespace Montjoie
{
  //! returns the size of memory used by the object
  template<int type>
  size_t PyramidReference<type>::GetMemorySize() const
  {
    size_t taille = ElementReference<Dimension3, type>::GetMemorySize();
    taille += sizeof(Real_wp)*(points1d_z.GetM()+weights1d_z.GetM());
    taille += Fb_geom.GetMemorySize();
    return taille;
  }

  
  //! how to number mesh
  template<int type> void PyramidReference<type>::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    nmap.SetFormulationDG(dg);
    if (dg == ElementReference_Base::DISCONTINUOUS)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofQuadrangle(this->order, 0);
	nmap.SetNbDofTriangle(this->order, 0);
	nmap.SetNbDofPyramid(this->order, this->nb_dof_loc);        
	return;
      }
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, (this->order+2)*(this->order+1)/2);
	nmap.SetNbDofQuadrangle(this->order, (this->order+1)*(this->order+1));
	nmap.SetNbDofPyramid(this->order, 0);

	// rotation of dofs on faces
        Matrix<int> FacesDof_Rotation_Tri;
        MeshNumbering<Dimension3>::
          GetRotationTriangularFace(Fb_geom.PointsDof2D_tri(), FacesDof_Rotation_Tri);
        
        nmap.SetFacesDofRotationTri(this->order, FacesDof_Rotation_Tri);
        
        Matrix<int>& NumQuad2D = this->elt_geom.GetNumQuad2D();
        Matrix<int> FacesDof_Rotation_Quad;
        MeshNumbering<Dimension3>::
          GetRotationQuadrilateralFace(NumQuad2D, FacesDof_Rotation_Quad);
    
        nmap.SetFacesDofRotationQuad(this->order, FacesDof_Rotation_Quad);
	
	return;
      }
    
    nmap.SetNbDofVertex(this->order, 1);
    nmap.SetNbDofEdge(this->order, this->order-1);
    nmap.SetNbDofTriangle(this->order, (this->order-1)*(this->order-2)/2);
    nmap.SetNbDofQuadrangle(this->order, (this->order-1)*(this->order-1));
    nmap.SetNbDofPyramid(this->order, this->nb_dof_loc - 3*this->order*this->order - 2);
    
    nmap.SetOppositeEdgesDofSymmetry(this->order, this->order-1);
    
    // rotation of dofs on faces
    VectR2 points_inside(this->nb_dof_tri - 3*this->order);
    for (int i = 0; i < points_inside.GetM(); i++)
      points_inside(i) = Fb_geom.PointsDof2D_tri()(i+3*this->order);
    
    Matrix<int> FacesDof_Rotation_Tri;
    MeshNumbering<Dimension3>::
      GetRotationTriangularFace(points_inside, FacesDof_Rotation_Tri);
    
    nmap.SetFacesDofRotationTri(this->order, FacesDof_Rotation_Tri);

    Matrix<int> NumNodes2D_InsideQuad(this->order-1, this->order-1);
    for (int i = 1;  i < this->order; i++)
      for (int j = 1;  j < this->order; j++)
	NumNodes2D_InsideQuad(i-1, j-1) = (i-1)*(this->order-1) + j-1;
    
    Matrix<int> FacesDof_Rotation_Quad;
    MeshNumbering<Dimension3>::
      GetRotationQuadrilateralFace(NumNodes2D_InsideQuad, FacesDof_Rotation_Quad);
    
    nmap.SetFacesDofRotationQuad(this->order, FacesDof_Rotation_Quad);
  }
  
  
  //! constructing finite element
  template<int type>
  void PyramidReference<type>::
  ConstructFiniteElement(int r, int rgeom, int rquad,
                         int type_quad, int rsurf_tri, int rsurf_quad,
                         int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (rgeom == 0)
      rgeom = r;
    
    if (rquad == 0)
      rquad = r;
    
    if (type_quad == -1)
      {
        if (ElementReference_Base::use_warburton_trick)
          type_quad = PyramidQuadrature::QUADRATURE_JACOBI1;
        else
          type_quad = PyramidQuadrature::QUADRATURE_JACOBI2;
        
        //type_quad = PyramidQuadrature::QUADRATURE_NON_PRODUCT;
      }
    
    this->order = r;
    this->SetGeometryOrder(rgeom);
    this->order_quad = rquad;
    
    // quadrature rules
    ConstructQuadrature(this->order_quad, type_quad,
                        rsurf_tri, rsurf_quad, type_surf_tri, type_surf_quad);
    
    // 2-D basis functions
    TriangleGeomReference& function_tri = Fb_geom.GetTriangularSurfaceFiniteElement();
    if (type_interpolation == PyramidGeomReference::REGULAR_BASIS)
      function_tri.type_interpolation_nodal = TriangleGeomReference::REGULAR_BASIS;
    else
      function_tri.type_interpolation_nodal = TriangleGeomReference::LOBATTO_BASIS;
    
    function_tri.ConstructFiniteElement(rgeom);
    
    // orthogonal functions
    Fb_geom.ComputeOrthogonalFunctions(max(r, rgeom));
    
    // shape functions (used for Fi)
    Fb_geom.ConstructNodalShapeFunctions(rgeom);
    this->FacesNode = Fb_geom.GetNodalNumber();
    
    this->nb_dof_loc = this->GetNbPointsNodalElt();
    
    Fb_geom.SetPointsDof2D_tri(Fb_geom.PointsNodal2D_tri());
    Fb_geom.SetPointsDof2D_quad(Fb_geom.PointsNodal2D_quad());    
  }
  
  
  //! construction of quadrature points and nodal points
  template<int type> 
  void PyramidReference<type>::
  ConstructQuadrature(int r, int type_quad, int rtri, int rquad,
                      int type_surf_tri, int type_surf_quad)
  {
    if (rtri == 0)
      rtri = r;
    
    if (rquad == 0)
      rquad = r;
    
    if (type_surf_tri == -1)
      type_surf_tri = TriangleQuadrature::QUADRATURE_GAUSS;
    
    if (type_surf_quad == -1)
      type_surf_quad = Globatto<Real_wp>::QUADRATURE_GAUSS;
        
    // 1-D quadrature points
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(rquad, type_surf_quad);
    
    VectReal_wp points1d, weights1d;
    points1d = gauss.Points();
    weights1d = gauss.Weights();
    
    VectReal_wp weights3d; VectR3 points3d;
    PyramidQuadrature::ConstructQuadrature(r, points3d, weights3d,
                                           this->points1d_z, this->weights1d_z, type_quad);
    
    Matrix<int> coor;
    Matrix<int>& NumQuad2D = this->elt_geom.GetNumQuad2D();
    MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(rquad, NumQuad2D, coor);
    
    // quadrature points on the triangle
    VectR2 points2d_tri; VectReal_wp weights2d_tri;
    TriangleQuadrature::
      ConstructQuadrature(2*rtri, points2d_tri, weights2d_tri, type_surf_tri);
        
    // quadrature points on the quadrilateral base
    VectR2 points2d_quad; VectReal_wp weights2d_quad;
    points2d_quad.Reallocate((rquad+1)*(rquad+1));
    weights2d_quad.Reallocate((rquad+1)*(rquad+1));
    
    for (int i = 0; i <= rquad; i++)
      for (int j = 0; j <= rquad; j++)
	{
	  points2d_quad(NumQuad2D(i, j)).Init(points1d(i), points1d(j));
	  weights2d_quad(NumQuad2D(i,j)) = weights1d(i)*weights1d(j);
	}
    
    // quadrature points on all the faces
    this->nb_points_quadrature_boundaries
      = points2d_quad.GetM() + 4*points2d_tri.GetM();
    
    this->nb_points_quadrature_inside = points3d.GetM();
    
    // we add quadrature points of the face at the end of points3d
    // Points of integration for the first face z = 0 (square)
    points3d
      .Resize(this->nb_points_quadrature_inside + this->nb_points_quadrature_boundaries);
    
    Real_wp two(2), one(1); int ind = this->nb_points_quadrature_inside;
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d(ind++) = R3(two*points2d_quad(i)(0)-one,
			   two*points2d_quad(i)(1)-one, Real_wp(0));
    
    // Points of integration for the second face z-y = 1 (triangle)
    for (int i = 0; i < points2d_tri.GetM(); i++)
      points3d(ind++) = R3(two*points2d_tri(i)(0) + points2d_tri(i)(1) - one,
			   points2d_tri(i)(1) - one, points2d_tri(i)(1));
    
    // Points of integration for the third face z+x = 1 (triangle)
    for (int i = 0; i < points2d_tri.GetM(); i++)
      points3d(ind++) = R3(one - points2d_tri(i)(1),
			   two*points2d_tri(i)(0) + points2d_tri(i)(1) - one,
			   points2d_tri(i)(1));
    
    // Points of integration for the fourth face z+y = 1 (triangle)
    for (int i = 0; i < points2d_tri.GetM(); i++)
      points3d(ind++) = R3(two*points2d_tri(i)(0) + points2d_tri(i)(1) - one,
			   one - points2d_tri(i)(1), points2d_tri(i)(1));
    
    // Points of integration for the fifth face z-x = 1 (triangle)
    for (int i = 0; i < points2d_tri.GetM(); i++)
      points3d(ind++) = R3(points2d_tri(i)(1) - Real_wp(1),
			   two*points2d_tri(i)(0) + points2d_tri(i)(1) - one,
			   points2d_tri(i)(1));
    
    // constructing array this->num_quad_points_surf
    // to get numbers of quadrature points of the faces
    ind = this->nb_points_quadrature_inside;
    this->num_quad_points_surf.Reallocate(5);
    this->num_quad_points_surf(0).Reallocate(points2d_quad.GetM());
    for (int n = 1; n < 5; n++)
      this->num_quad_points_surf(n).Reallocate(points2d_tri.GetM());

    // first face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      this->num_quad_points_surf(0)(i) = ind++;
    
    // second face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      this->num_quad_points_surf(1)(i) = ind++;
    
    // third face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      this->num_quad_points_surf(2)(i) = ind++;
    
    // fourth face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      this->num_quad_points_surf(3)(i) = ind++;
      
    // fifth face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      this->num_quad_points_surf(4)(i) = ind++;
    
    this->offset_faceSh.Reallocate(6);
    this->offset_faceSh(0) = 0;
    for (int i = 1; i < 5; i++)
      this->offset_faceSh(i) = points2d_quad.GetM() + (i-1)*points2d_tri.GetM();
    
    this->offset_faceSh(5) = this->nb_points_quadrature_boundaries;
    
    this->PoidsFlux(0) = weights2d_quad;
    Mlt(Real_wp(0.5), this->PoidsFlux(0));
    for (int k = 1; k < 5; k++)
      {
        this->PoidsFlux(k) = weights2d_tri;
        Mlt(Real_wp(0.5), this->PoidsFlux(k));
      }

    this->SetPoints1D(points1d); this->SetWeights1D(weights1d);

    Fb_geom.SetPoints2D_tri(points2d_tri);
    Fb_geom.SetWeights2D_tri(weights2d_tri);
    Fb_geom.SetPoints2D_quad(points2d_quad);
    Fb_geom.SetWeights2D_quad(weights2d_quad);
    
    this->SetPointsND(points3d);
    this->SetWeightsND(weights3d);

  }
  
  //! displays informations class PyramidReference<type>
  template<int type>
  ostream& operator <<(ostream& out,const PyramidReference<type>& e)
  {
    out << static_cast<const ElementReference<Dimension3, type>& >(e);
    out<<"Type interpolation "<<e.type_interpolation<<endl;
    return out;
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_REFERENCE_CXX
#endif
