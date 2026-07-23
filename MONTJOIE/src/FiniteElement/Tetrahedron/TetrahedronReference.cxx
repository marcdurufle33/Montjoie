#ifndef MONTJOIE_FILE_TETRAHEDRON_REFERENCE_CXX

namespace Montjoie
{
  //! returns the size of memory used by the object
  template<int type>
  size_t TetrahedronReference<type>::GetMemorySize() const
  {
    size_t taille = ElementReference<Dimension3, type>::GetMemorySize();
    taille += Fb_geom.GetMemorySize();
    return taille;
  }
  
  
  /****************************
   * Initialization functions *
   ****************************/
  
  
  //! how to number mesh
  template<int type>
  void TetrahedronReference<type>::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    nmap.SetFormulationDG(dg);
    if (dg == ElementReference_Base::DISCONTINUOUS)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofQuadrangle(this->order, 0);
	nmap.SetNbDofTriangle(this->order, 0);
	nmap.SetNbDofTetrahedron(this->order, this->nb_dof_loc);        
	return;
      }
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, (this->order+2)*(this->order+1)/2);
	nmap.SetNbDofQuadrangle(this->order, (this->order+1)*(this->order+1));
	nmap.SetNbDofTetrahedron(this->order, 0);

	// rotation of dofs on faces
        Matrix<int> FacesDof_Rotation_Tri;
        MeshNumbering<Dimension3>::
          GetRotationTriangularFace(Fb_geom.PointsDof2D_tri(), FacesDof_Rotation_Tri);
        
        //DISP(FacesDof_Rotation_Tri);
        nmap.SetFacesDofRotationTri(this->order, FacesDof_Rotation_Tri);
        
	return;
      }
    
    nmap.SetNbDofVertex(this->order, 1);
    nmap.SetNbDofEdge(this->order, this->order-1);
    nmap.SetNbDofQuadrangle(this->order, (this->order-1)*(this->order-1));
    nmap.SetNbDofTriangle(this->order, (this->order-2)*(this->order-1)/2);
    nmap.SetNbDofTetrahedron(this->order, (this->order-3)*(this->order-2)*(this->order-1)/6);
    
    nmap.SetOppositeEdgesDofSymmetry(this->order, this->order-1);
    
    // rotation of dofs on faces
    VectR2 points_inside(this->nb_dof_tri - 3*this->order);
    for (int i = 0; i < points_inside.GetM(); i++)
      points_inside(i) = Fb_geom.PointsDof2D_tri()(i+3*this->order);
    
    Matrix<int> FacesDof_Rotation_Tri;
    MeshNumbering<Dimension3>::
      GetRotationTriangularFace(points_inside, FacesDof_Rotation_Tri);
    
    nmap.SetFacesDofRotationTri(this->order, FacesDof_Rotation_Tri);
  }
  
  
  //! constructing finite element
  template<int type>
  void TetrahedronReference<type>
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad, int rsurf_tri,
			   int rsurf_quad, int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (rgeom == 0)
      rgeom = r;
    
    if (rquad == 0)
      rquad = r;
    
    if (type_quad == -1)
      type_quad = TetrahedronQuadrature::QUADRATURE_GAUSS;
    
    this->order = r;
    this->SetGeometryOrder(rgeom);
    this->order_quad = rquad;
    
    // quadrature rules
    ConstructQuadrature(rquad, type_quad);
    
    // orthogonal functions
    Fb_geom.ComputeOrthogonalFunctions(max(r, rgeom));
    
    // shape functions (used for Fi)
    Fb_geom.ConstructNodalShapeFunctions(rgeom);
    this->FacesNode = Fb_geom.GetNodalNumber();
    
    this->nb_dof_tri = Fb_geom.PointsNodal2D_tri().GetM();
    this->nb_dof_quad = 0;
    this->nb_dof_loc = this->GetNbPointsNodalElt();
    
    Fb_geom.SetPointsDof2D_tri(Fb_geom.PointsNodal2D_tri());
    Fb_geom.SetPointsDof2D_quad(Fb_geom.PointsNodal2D_quad());    
  }
  
  
  //! construction of quadrature points and nodal points
  template<int type>
  void TetrahedronReference<type>::ConstructQuadrature(int r, int type_quad)
  {
    // 1-D points
    VectReal_wp points1d, weights1d;
    ComputeGaussLegendre(points1d, weights1d, r);
    
    // 2-D points
    VectR2 points2d_tri; VectReal_wp weights2d_tri;
    TriangleQuadrature::ConstructQuadrature(2*r, points2d_tri, weights2d_tri);
    
    int p = 2*r;
    if (r >= 5)
      p++;
    
    // quadrature points on tetrahedron
    VectR3 points3d; VectReal_wp weights3d;
    TetrahedronQuadrature::ConstructQuadrature(p, points3d, weights3d, type_quad);
    
    this->nb_points_quadrature_boundaries = 4*points2d_tri.GetM();
    this->nb_points_quadrature_inside = points3d.GetM();
      
    // we add quadrature points on the faces
    points3d
      .Resize(this->nb_points_quadrature_inside + this->nb_points_quadrature_boundaries);
    
    int nb = this->nb_points_quadrature_inside;
    this->num_quad_points_surf.Reallocate(4);
    for (int n = 0; n < 4; n++)
      this->num_quad_points_surf(n).Reallocate(points2d_tri.GetM());

    // Points of integration for the first face z = 0
    for (int i = 0; i < points2d_tri.GetM(); i++)
      {
	this->num_quad_points_surf(0)(i) = nb;
	points3d(nb++).Init(points2d_tri(i)(0), points2d_tri(i)(1), Real_wp(0));
      }
    
    // Points of integration for the second face y = 0
    for (int i = 0; i < points2d_tri.GetM(); i++)
      {
	this->num_quad_points_surf(1)(i) = nb;
	points3d(nb++).Init(points2d_tri(i)(0), Real_wp(0), points2d_tri(i)(1));
      }
    
    // Points of integration for the third face x = 0
    for (int i = 0; i < points2d_tri.GetM(); i++)
      {
	this->num_quad_points_surf(2)(i) = nb;
	points3d(nb++).Init(Real_wp(0), points2d_tri(i)(0), points2d_tri(i)(1));
      }
    
    // Points of integration for the fourth face x+y+z = 1
    for (int i = 0; i < points2d_tri.GetM(); i++)
      {
	this->num_quad_points_surf(3)(i) = nb;
	points3d(nb++).Init(Real_wp(1) - points2d_tri(i)(0) - points2d_tri(i)(1),
			    points2d_tri(i)(0), points2d_tri(i)(1));
      }
    
    this->offset_faceSh.Reallocate(5); this->offset_faceSh.Fill();
    Mlt(int(points2d_tri.GetM()), this->offset_faceSh);
    
    for (int k = 0; k < 4; k++)
      {
        this->PoidsFlux(k) = weights2d_tri;
        Mlt(Real_wp(0.5), this->PoidsFlux(k));
      }
    
    this->SetPoints1D(points1d); this->SetWeights1D(weights1d);
    Fb_geom.SetPoints2D_tri(points2d_tri);
    Fb_geom.SetWeights2D_tri(weights2d_tri);
    
    this->SetPointsND(points3d);
    this->SetWeightsND(weights3d);

  }
  
  
  //! displays informations class TetrahedronReference<type>
  template<int type> ostream& operator <<(ostream& out, const TetrahedronReference<type>& e)
  {
    out<<static_cast<const VolumeReference<type>&>(e);
    return out;
  }
  
} // end namespace

#define MONTJOIE_FILE_TETRAHEDRON_REFERENCE_CXX
#endif
