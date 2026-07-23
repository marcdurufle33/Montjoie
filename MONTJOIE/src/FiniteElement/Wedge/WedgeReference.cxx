#ifndef MONTJOIE_FILE_WEDGE_REFERENCE_CXX

namespace Montjoie
{
  //! returns the size of memory used by the object
  template<int type>
  size_t WedgeReference<type>::GetMemorySize() const
  {
    size_t taille = ElementReference<Dimension3, type>::GetMemorySize();
    taille += sizeof(Real_wp)*(points1d_z.GetM()+weights1d_z.GetM());
    taille += sizeof(int)*(NumQuadTri.GetDataSize()
                           +CoordinateQuad.GetDataSize());
    taille += Fb_geom.GetMemorySize();
    return taille;
  }

  
  //! how to number mesh
  template<int type>
  void WedgeReference<type>::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    nmap.SetFormulationDG(dg);
    if (dg == ElementReference_Base::DISCONTINUOUS)
      {
	nmap.SetNbDofWedge(this->order, this->nb_dof_loc);
	return;
      }
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, (this->order+2)*(this->order+1)/2);
	nmap.SetNbDofQuadrangle(this->order, (this->order+1)*(this->order+1));
	nmap.SetNbDofWedge(this->order, 0);

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
    nmap.SetNbDofWedge(this->order, (this->order-1)*(this->order-1)*(this->order-2)/2);
    
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
  void WedgeReference<type>::
  ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                         int rsurf_tri, int rsurf_quad,
                         int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (rgeom == 0)
      rgeom = r;
    
    if (rquad == 0)
      rquad = r;

    if (type_quad == -1)
      type_quad = TriangleQuadrature::QUADRATURE_GAUSS;
    
    if (gauss_z == -1)
      gauss_z = Globatto<Real_wp>::QUADRATURE_GAUSS;
    
    this->order = r;
    this->SetGeometryOrder(rgeom);
    this->order_quad = rquad;
    
    // quadrature
    ConstructQuadrature(rquad, type_quad, gauss_z, rsurf_tri,
                        rsurf_quad, type_surf_tri, type_surf_quad);
    
    // nodal functions
    Fb_geom.ConstructNodalShapeFunctions(rgeom);
    this->FacesNode = Fb_geom.GetNodalNumber();
    
    Fb_geom.ConstructOrthogonalFunctions(max(r, rgeom));
    this->nb_dof_loc = this->GetNbPointsNodalElt();

    this->SetPointsDofND(this->PointsNodalND());    
    Fb_geom.SetPointsDof2D_tri(Fb_geom.PointsNodal2D_tri());
    Fb_geom.SetPointsDof2D_quad(Fb_geom.PointsNodal2D_quad());
  }
  
  
  //! construction of quadrature points and nodal points
  template<int type>
  void WedgeReference<type>::
  ConstructQuadrature(int r, int type_quad, int gauss_z, int rtri, int rquad,
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

    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(rquad, type_surf_quad);

    VectReal_wp points1d, weights1d;
    points1d = gauss.Points();
    weights1d = gauss.Weights();
    
    if (gauss_z == Globatto<Real_wp>::QUADRATURE_GAUSS)
      ComputeGaussLegendre(this->points1d_z, this->weights1d_z, r);
    else
      ComputeGaussLobatto(this->points1d_z, this->weights1d_z, r);
    
    // quadrature points on triangle (for 3-D internal points)
    VectR2 points2d; VectReal_wp weights2d;
    TriangleQuadrature::ConstructQuadrature(2*r, points2d, weights2d, type_quad);
    
    // for triangular faces
    VectR2 points2d_tri; VectReal_wp weights2d_tri;
    TriangleQuadrature::ConstructQuadrature(2*rtri, points2d_tri,
                                            weights2d_tri, type_surf_tri);

    // quadrature points on quadrangle
    Matrix<int> coor; VectR2 points2d_quad; VectReal_wp weights2d_quad;
    Matrix<int>& NumQuad2D = this->elt_geom.GetNumQuad2D();        
    MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(rquad, NumQuad2D, coor);
    points2d_quad.Reallocate((rquad+1)*(rquad+1));
    weights2d_quad.Reallocate((rquad+1)*(rquad+1));    
    for (int i = 0; i <= rquad; i++)
      for (int j = 0; j <= rquad; j++)
	{
	  points2d_quad(NumQuad2D(i, j)).Init(points1d(i), points1d(j));
	  weights2d_quad(NumQuad2D(i,j)) = weights1d(i)*weights1d(j);
	}
        
    int Nq = points2d.GetM()*(r+1);
    this->nb_points_quadrature_inside = Nq;
    VectR3 points3d; VectReal_wp weights3d;
    points3d.Reallocate(Nq + 2*points2d_tri.GetM() + 3*points2d_quad.GetM());
    weights3d.Reallocate(Nq);
    CoordinateQuad.Reallocate(points3d.GetM(), 3);
    NumQuadTri.Reallocate(points2d.GetM(), r+1);
    int k = 0;
    for (int i = 0; i < points2d.GetM(); i++)    
      for (int j = 0; j < r+1; j++)  
	{
	  points3d(k)(0) = points2d(i)(0);
	  points3d(k)(1) = points2d(i)(1);
	  points3d(k)(2) = this->points1d_z(j);
	  weights3d(k) = weights2d(i)*this->weights1d_z(j);
	  NumQuadTri(i, j) = k;
	  CoordinateQuad(k, 0) = i;
	  CoordinateQuad(k, 1) = j;
          CoordinateQuad(k, 2) = -1;
	  k++;
	}

    // quadrature points on all the faces
    this->nb_points_quadrature_boundaries
      = 3*points2d_quad.GetM() + 2*points2d_tri.GetM();
   
    // Points of integration for the first face z = 0 (triangle)
    for (int i = 0; i < points2d_tri.GetM(); i++)
      points3d(Nq+i) = R3(points2d_tri(i)(0), points2d_tri(i)(1), Real_wp(0));
    
    Nq += points2d_tri.GetM();
    // Points of integration for the second face y = 0 (quadrangle)
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d(Nq+i) = R3(points2d_quad(i)(0), 0, points2d_quad(i)(1));
    
    Nq += points2d_quad.GetM();
    // Points of integration for the second face x+y = 1 (quadrangle)
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d(Nq + i) = R3(Real_wp(1)-points2d_quad(i)(0),
			    points2d_quad(i)(0), points2d_quad(i)(1));
    
    Nq += points2d_quad.GetM();
    // Points of integration for the fourth face x = 0 (quadrangle)
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d(Nq + i) = R3(0, points2d_quad(i)(0), points2d_quad(i)(1));
    
    Nq += points2d_quad.GetM();
    // Points of integration for the fifth face z = 1 (triangle)
    for (int i = 0; i < points2d_tri.GetM(); i++)
      points3d(Nq + i) = R3(points2d_tri(i)(0), points2d_tri(i)(1), Real_wp(1));
    
    
    int ind = this->nb_points_quadrature_inside;
    this->num_quad_points_surf.Reallocate(5);
    this->num_quad_points_surf(0).Reallocate(points2d_tri.GetM());
    this->num_quad_points_surf(4).Reallocate(points2d_tri.GetM());
    for (int n = 1; n < 4; n++)
      this->num_quad_points_surf(n).Reallocate(points2d_quad.GetM());

    // first face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      {
        CoordinateQuad(ind, 0) = i;
        CoordinateQuad(ind, 1) = 0;
        CoordinateQuad(ind, 2) = 0;
        this->num_quad_points_surf(0)(i) = ind++;
      }

    
    // second face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      {
        CoordinateQuad(ind, 0) = coor(i, 0);
        CoordinateQuad(ind, 1) = coor(i, 1);
        CoordinateQuad(ind, 2) = 1;
        this->num_quad_points_surf(1)(i) = ind++;
      }
    
    // third face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      {
        CoordinateQuad(ind, 0) = coor(i, 0);
        CoordinateQuad(ind, 1) = coor(i, 1);
        CoordinateQuad(ind, 2) = 2;
        this->num_quad_points_surf(2)(i) = ind++;
      }
    
    // fourth face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      {
        CoordinateQuad(ind, 0) = coor(i, 0);
        CoordinateQuad(ind, 1) = coor(i, 1);
        CoordinateQuad(ind, 2) = 3;
	this->num_quad_points_surf(3)(i) = ind++;
      }
      
    // fifth face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      {
        CoordinateQuad(ind, 0) = i;
        CoordinateQuad(ind, 1) = r;
        CoordinateQuad(ind, 2) = 4;
        this->num_quad_points_surf(4)(i) = ind++;
      }
    
    this->offset_faceSh.Reallocate(6);
    this->offset_faceSh(0) = 0;
    for (int i = 1; i < 4; i++)
      this->offset_faceSh(i) = points2d_tri.GetM() + (i-1)*points2d_quad.GetM();
    
    this->offset_faceSh(4) = points2d_tri.GetM() + 3*points2d_quad.GetM();
    
    this->nb_points_quadrature_boundaries
      = 3*points2d_quad.GetM() + 2*points2d_tri.GetM();
    this->offset_faceSh(5) = this->nb_points_quadrature_boundaries;
    
    this->PoidsFlux(0) = weights2d_tri;
    Mlt(Real_wp(0.5), this->PoidsFlux(0));
    this->PoidsFlux(4) = weights2d_tri;
    Mlt(Real_wp(0.5), this->PoidsFlux(4));
    for (int k = 1; k < 4; k++)
      {
        this->PoidsFlux(k) = weights2d_quad;
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

  
  //! displays informations class WedgeReference<type>
  template<int type> ostream& operator <<(ostream& out, const WedgeReference<type>& e)
  {
    out<<static_cast<const VolumeReference<type>&>(e);
    return out;
  }
  
} // end namespace

#define MONTJOIE_FILE_WEDGE_REFERENCE_CXX
#endif
