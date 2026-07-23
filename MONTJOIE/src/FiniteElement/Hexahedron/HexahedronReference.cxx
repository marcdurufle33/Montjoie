#ifndef MONTJOIE_FILE_HEXAHEDRON_REFERENCE_CXX

namespace Montjoie
{

  //! returns size of memory used by the object
  template<int type>
  size_t HexahedronReference<type>::GetMemorySize() const
  {
    size_t taille = ElementReference<Dimension3, type>::GetMemorySize();
    taille += lob_basis.GetMemorySize()+lob_quad.GetMemorySize();
    taille += NumQuad3D.GetMemorySize();
    taille += CoordinateQuad3D.GetMemorySize();
    taille += CoordinateQuad2D.GetMemorySize();
    taille += Fb_geom.GetMemorySize();
    return taille;
  }
  
  
  //! how to number mesh
  template<int type> void HexahedronReference<type>::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    // dg formulation
    nmap.SetFormulationDG(dg);

    if (dg == ElementReference_Base::DISCONTINUOUS)
      {
	// dofs are only inside hexahedra
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, 0);
	nmap.SetNbDofQuadrangle(this->order, 0);
	nmap.SetNbDofHexahedron(this->order, this->nb_dof_loc);
	return;
      }
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofQuadrangle(this->order, (this->order+1)*(this->order+1));
	nmap.SetNbDofHexahedron(this->order, 0);

	// rotation of dofs on faces
        Matrix<int>& NumQuad2D = this->elt_geom.GetNumQuad2D();
        Matrix<int> FacesDof_Rotation_Quad;
        MeshNumbering<Dimension3>::
          GetRotationQuadrilateralFace(NumQuad2D, FacesDof_Rotation_Quad);

        nmap.SetFacesDofRotationQuad(this->order, FacesDof_Rotation_Quad);
        
	return;
      }
    
    nmap.SetNbDofVertex(this->order, 1);
    nmap.SetNbDofEdge(this->order, this->order-1);
    nmap.SetNbDofQuadrangle(this->order, (this->order-1)*(this->order-1));
    if (serendip_element)
      nmap.SetNbDofQuadrangle(this->order, 0);
        
    nmap.SetNbDofHexahedron(this->order, (this->order-1)*(this->order-1)*(this->order-1));
    	
    nmap.SetOppositeEdgesDofSymmetry(this->order, this->order-1);
    
    // rotation of dofs on faces
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
  void HexahedronReference<type>
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                           int rsurf_tri, int rsurf, int type_surf_tri,
                           int type_quad_surf, int gauss_z)
  {
    if (rgeom == 0)
      rgeom = r;
    
    if (rquad == 0)
      rquad = r;
    
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
        
    this->order = r;
    this->SetGeometryOrder(rgeom);
    this->order_quad = rquad;
    
    ConstructQuadrature(rquad, type_quad, rsurf, type_quad_surf);
    
    Fb_geom.ConstructNodalShapeFunctions(rgeom);
    
    this->FacesNode = Fb_geom.GetNodalNumber();
    this->nb_dof_tri = 0;
    this->nb_dof_quad = Fb_geom.PointsNodal2D_quad().GetM();
    
    this->SetPointsDofND(this->PointsNodalND());
    
    this->nb_dof_loc = this->GetNbPointsNodalElt();
    
    Fb_geom.SetPointsDof2D_tri(Fb_geom.PointsNodal2D_tri());
    Fb_geom.SetPointsDof2D_quad(Fb_geom.PointsNodal2D_quad());
    
    Fb_geom.ComputeCurvedTransformation();
    Fb_geom.ComputeCoefficientTransformation();
  }
  
  
  //! construction of quadrature points and nodal points
  template<int type>
  void HexahedronReference<type>::ConstructQuadrature(int r, int type_quad,
						      int rsurf, int type_surf)
  {
    if (rsurf == 0)
      rsurf = r;
    
    if (type_surf == -1)
      type_surf = type_quad;
    
    lob_quad.ConstructQuadrature(r, type_quad);
    this->SetPoints1D(lob_quad.Points());
    this->SetWeights1D(lob_quad.Weights());
    
    Matrix<int>& NumQuad2D = this->elt_geom.GetNumQuad2D();
    MeshNumbering<Dimension3>::ConstructHexahedralNumbering(r, NumQuad3D, CoordinateQuad3D);    
    MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(r, NumQuad2D, CoordinateQuad2D);

    const VectReal_wp& weights1d = this->Weights1D();
    const VectReal_wp& points1d = this->Points1D();
    VectR2 points2d_quad;
    VectReal_wp weights2d_quad;
    weights2d_quad.Reallocate((r+1)*(r+1));
    points2d_quad.Reallocate((r+1)*(r+1));
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  weights2d_quad(NumQuad2D(i, j)) = weights1d(i)*weights1d(j);
	  points2d_quad(NumQuad2D(i, j)).Init(points1d(i), points1d(j));
	}
    
    // quadrature points
    VectR3 points3d; VectReal_wp weights3d;
    points3d.Reallocate((r+1)*(r+1)*(r+1));
    this->nb_points_quadrature_inside = points3d.GetM();
    this->nb_points_quadrature_boundaries = 6*points2d_quad.GetM();
    
    // points and weights of integration
    weights3d.Reallocate(points3d.GetM());
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	for (int k = 0; k <= r; k++)
	  {
	    points3d(NumQuad3D(i, j, k))
              .Init(points1d(i), points1d(j), points1d(k));
	    
	    weights3d(NumQuad3D(i, j, k))
              = weights1d(i)*weights1d(j)*weights1d(k);
	  }

    this->num_quad_points_surf.Reallocate(6);
    for (int n = 0; n < 6; n++)
      this->num_quad_points_surf(n).Reallocate(points2d_quad.GetM());

    if (type_quad == Globatto<Real_wp>::QUADRATURE_LOBATTO)
      {
        // quadrature points on first face
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            this->num_quad_points_surf(0)(NumQuad2D(i, j)) = NumQuad3D(0, i, j);
        
        // quadrature points on second face
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            this->num_quad_points_surf(1)(NumQuad2D(i, j)) = NumQuad3D(i, 0, j);
        
        // quadrature points on third face
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            this->num_quad_points_surf(2)(NumQuad2D(i, j)) = NumQuad3D(i, j, 0);
        
        // quadrature points on fourth face
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            this->num_quad_points_surf(3)(NumQuad2D(i, j)) = NumQuad3D(i, j, r);
        
        // quadrature points on fifth face
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            this->num_quad_points_surf(4)(NumQuad2D(i, j)) = NumQuad3D(i, r, j);
        
        // quadrature points on sixth face
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            this->num_quad_points_surf(5)(NumQuad2D(i, j)) = NumQuad3D(r, i, j);
      }
    else
      {
        points3d.Resize(points3d.GetM() + 6*points2d_quad.GetM());
        int nb = this->nb_points_quadrature_inside;
        // Points of integration for the first face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(0.0, points2d_quad(i)(0), points2d_quad(i)(1));
        
        // Points of integration for the second face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(points2d_quad(i)(0), 0.0, points2d_quad(i)(1));
        
        // Points of integration for the third face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(points2d_quad(i)(0), points2d_quad(i)(1), 0.0);
        
        // Points of integration for the fourth face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(points2d_quad(i)(0), points2d_quad(i)(1), 1.0);
        
        // Points of integration for the fifth face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(points2d_quad(i)(0), 1.0, points2d_quad(i)(1));
        
        // Points of integration for the sixth face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(1.0, points2d_quad(i)(0), points2d_quad(i)(1));
        
        nb = this->nb_points_quadrature_inside;
        for (int num_loc = 0; num_loc < 6; num_loc++)
          for (int i = 0; i < points2d_quad.GetM(); i++)
            this->num_quad_points_surf(num_loc)(i) = nb++;
        
      }
    
    this->offset_faceSh.Reallocate(7); this->offset_faceSh.Fill();
    Mlt(int(points2d_quad.GetM()), this->offset_faceSh);
    
    for (int k = 0; k < 6; k++)
      {
        this->PoidsFlux(k) = weights2d_quad;
        Mlt(Real_wp(0.5), this->PoidsFlux(k));
      }
    
    if ( (rsurf != r) || (type_surf != type_quad))
      {
	Globatto<Real_wp> lob_surf;
	lob_surf.ConstructQuadrature(rsurf, type_surf);
	
	MeshNumbering<Dimension2>::
          ConstructQuadrilateralNumbering(rsurf, NumQuad2D, CoordinateQuad2D);

	weights2d_quad.Reallocate((rsurf+1)*(rsurf+1));
	points2d_quad.Reallocate((rsurf+1)*(rsurf+1));
	for (int i = 0; i <= rsurf; i++)
	  for (int j = 0; j <= rsurf; j++)
	    {
	      weights2d_quad(NumQuad2D(i, j)) = lob_surf.Weights(i)*lob_surf.Weights(j);
	      points2d_quad(NumQuad2D(i, j)).Init(lob_surf.Points(i), lob_surf.Points(j));
	    }
	
	this->SetPoints1D(lob_surf.Points());
	this->SetWeights1D(lob_surf.Weights());

	this->nb_points_quadrature_boundaries = 6*points2d_quad.GetM();
        points3d.Resize(this->nb_points_quadrature_inside + 6*points2d_quad.GetM());
        int nb = this->nb_points_quadrature_inside;
        // Points of integration for the first face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(0.0, points2d_quad(i)(0), points2d_quad(i)(1));
        
        // Points of integration for the second face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(points2d_quad(i)(0), 0.0, points2d_quad(i)(1));
        
        // Points of integration for the third face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(points2d_quad(i)(0), points2d_quad(i)(1), 0.0);
        
        // Points of integration for the fourth face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(points2d_quad(i)(0), points2d_quad(i)(1), 1.0);
        
        // Points of integration for the fifth face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(points2d_quad(i)(0), 1.0, points2d_quad(i)(1));
        
        // Points of integration for the sixth face
        for (int i = 0; i < points2d_quad.GetM(); i++)
          points3d(nb++).Init(1.0, points2d_quad(i)(0), points2d_quad(i)(1));

	this->num_quad_points_surf.Reallocate(6);
	for (int n = 0; n < 6; n++)
	  this->num_quad_points_surf(n).Reallocate(points2d_quad.GetM());
        
        nb = this->nb_points_quadrature_inside;
        for (int num_loc = 0; num_loc < 6; num_loc++)
          for (int i = 0; i < points2d_quad.GetM(); i++)
            this->num_quad_points_surf(num_loc)(i) = nb++;
	
	this->offset_faceSh.Reallocate(7); this->offset_faceSh.Fill();
	Mlt(int(points2d_quad.GetM()), this->offset_faceSh);
	
	for (int k = 0; k < 6; k++)
	  {
	    this->PoidsFlux(k) = weights2d_quad;
	    Mlt(Real_wp(0.5), this->PoidsFlux(k));
	  }
      }    
    
    Fb_geom.SetPoints2D_quad(points2d_quad);
    Fb_geom.SetWeights2D_quad(weights2d_quad);
    
    this->SetPointsND(points3d);
    this->SetWeightsND(weights3d);
  }
  
  
  //! displays informations class HexahedronReference<type>
  template<int type>
  ostream& operator <<(ostream& out, const HexahedronReference<type>& e)
  {
    out << static_cast<const VolumeReference<type>&>(e);
    return out;
  }
  
} // end namespace

#define MONTJOIE_FILE_HEXAHEDRON_REFERENCE_CXX
#endif
