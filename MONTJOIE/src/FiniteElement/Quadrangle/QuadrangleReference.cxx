#ifndef MONTJOIE_FILE_QUADRANGLE_REFERENCE_CXX

namespace Montjoie
{
  //! default constructor
  template<int type>
  QuadrangleReference<type>::QuadrangleReference() : ElementReference<Dimension2, type>(Fb_geom)
  {
  }
  
  
  //! returns the size of memory used by the object
  template<int type>
  size_t QuadrangleReference<type>::GetMemorySize() const
  {
    size_t taille = ElementReference<Dimension2, type>::GetMemorySize();
    taille += lob_quad.GetMemorySize()+lob_basis.GetMemorySize();
    taille += sizeof(int)*(NumQuad2D.GetDataSize()+CoordinateQuad.GetDataSize());
    taille += Fb_geom.GetMemorySize();
    return taille;
  }


  //! returns a reference to the object storing shape functions
  template<int type>
  QuadrangleGeomReference& QuadrangleReference<type>::GetGeometricElement()
  {
    return Fb_geom;
  }
  
  
  //! returns a reference to the object storing shape functions
  template<int type>
  const QuadrangleGeomReference& QuadrangleReference<type>::GetGeometricElement() const
  {
    return Fb_geom;
  }
  
    
  /****************************
   * Initialization functions *
   ****************************/
  
  
  //! constructing finite element
  template<int type>
  void QuadrangleReference<type>::
  ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad, int rsurf, int type_surf)
  {
    if (rgeom == 0)
      rgeom = r;
    
    if (rquad == 0)
      rquad = r;
    
    if (rsurf == 0)
      rsurf = rquad;
    
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
    
    this->order = r;
    this->SetGeometryOrder(rgeom);
    this->order_quad = rquad;

    // 1-D finite element 
    ConstructFiniteElement1D(r, r, rquad, type_surf);
        
    Fb_geom.ConstructNodalShapeFunctions(rgeom);
    this->EdgesNodal = this->GetNodalNumber();
    ConstructQuadrature(rquad, type_quad, rsurf);
    
    Fb_geom.ComputeCurvedTransformation();
    Fb_geom.ComputeCoefficientTransformation();
  }


  //! constructs 1-D finite element (restriction of the 2-D finite element to an edge)
  template<int type>
  void QuadrangleReference<type>
  ::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {
    if (this->element_surface != NULL)
      delete this->element_surface;

    EdgeGaussReference* edge = new EdgeGaussReference();
    edge->ConstructFiniteElement(r, rgeom, rquad, type_quad);    
    this->element_surface = edge;    
  }
  
  
  //! construction of quadrature points in 2-D
  /*!
    \param[in] order_ order of approximation
  */
  template<int type>
  void QuadrangleReference<type>::ConstructQuadrature(int r, int type_quad, int rsurf)
  {    
    if (rsurf == 0)
      rsurf = r;
    
    // quadrature formulas on edges
    lob_quad.ConstructQuadrature(rsurf, type_quad);
    
    this->SetWeights1D(lob_quad.Weights());
    this->SetPoints1D(lob_quad.Points());
    
    // quadrature formulas inside
    MeshNumbering<Dimension2>::
      ConstructQuadrilateralNumbering(r, NumQuad2D, CoordinateQuad);
    
    lob_quad.ConstructQuadrature(r, type_quad);
    this->nb_points_quadrature_inside = (r+1)*(r+1);
    VectR2 points2d; VectReal_wp weights2d;
    weights2d.Reallocate(this->nb_points_quadrature_inside);
    points2d.Reallocate(this->nb_points_quadrature_inside);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  weights2d(NumQuad2D(i,j)) = lob_quad.Weights(i)*lob_quad.Weights(j);
	  points2d(NumQuad2D(i,j)).Init(lob_quad.Points(i), lob_quad.Points(j));
	}

    this->num_quad_points_surf.Reallocate(4);
    for (int n = 0; n < 4; n++)
      this->num_quad_points_surf(n).Reallocate(rsurf+1);

    if (type_quad == Globatto<Real_wp>::QUADRATURE_LOBATTO)
      {
        if (r != rsurf)
          {
            abort();
          }
        
        // the quadrature points  on edges are Gauss-Lobatto points
        // First Edge
        for (int i = 0; i <= r; i++)
          this->num_quad_points_surf(0)(i) = NumQuad2D(i,0);
        
        // Second Edge
        for (int i = 0; i <= r; i++)
          this->num_quad_points_surf(1)(i) = NumQuad2D(r,i);
        
        // Third Edge
        for (int i = 0; i <= r; i++)
          this->num_quad_points_surf(2)(i) = NumQuad2D(r-i,r);
        
        // Fourth Edge
        for (int i = 0; i <= r; i++)
          this->num_quad_points_surf(3)(i) = NumQuad2D(0,r-i);
        
      }
    else
      {
        // Gauss points on edges
        points2d.Resize(points2d.GetM() + 4*(rsurf+1));
        
        int nb = (r+1)*(r+1);
        // Points of integration for the first edge
        for (int i = 0; i <= rsurf; i++)
          {
            this->num_quad_points_surf(0)(i) = nb;
            points2d(nb++).Init(this->Points1D(i), 0.0);
          }
        
        // Points of integration for the second edge
        for (int i = 0; i <= rsurf; i++)
          {
            this->num_quad_points_surf(1)(i) = nb;
            points2d(nb++).Init(1.0, this->Points1D(i));
          }
        
        // Points of integration for the third edge
        for (int i = 0; i <= rsurf; i++)
          {
            this->num_quad_points_surf(2)(i) = nb;
            points2d(nb++).Init(1.0-this->Points1D(i), 1.0);
          }
        
        // Points of integration for the fourth edge
        for (int i = 0; i <= rsurf; i++)
          {
            this->num_quad_points_surf(3)(i) = nb;
            points2d(nb++).Init(0.0, 1.0-this->Points1D(i));
          }
        
      }
    
    this->nb_points_quadrature_boundaries = 4*(rsurf+1);
    this->offset_faceSh.Reallocate(5); this->offset_faceSh.Fill();
    Mlt(rsurf+1, this->offset_faceSh);
    
    for (int k = 0; k < 4; k++)
      {
        this->PoidsFlux(k) = this->Weights1D();
        Mlt(Real_wp(0.5), this->PoidsFlux(k));
      }

    this->SetWeightsND(weights2d);
    this->SetPointsND(points2d);
  }
  
  
  //! method useful to number mesh
  template<int type>
  void QuadrangleReference<type>::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    int r = this->order;
    nmap.SetFormulationDG(dg);
    if (dg == ElementReference_Base::DISCONTINUOUS)
      {
	nmap.SetNbDofVertex(r, 0);
	nmap.SetNbDofEdge(r, 0);
	nmap.SetNbDofQuadrangle(r, this->nb_dof_loc);

	return;
      }
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, this->order+1);
	nmap.SetNbDofQuadrangle(this->order, 0);

	nmap.SetOppositeEdgesDofSymmetry(r, r+1);
	
	return;
      }
    
    nmap.SetNbDofVertex(r, 1);
    nmap.SetNbDofEdge(r, r-1);
    nmap.SetNbDofQuadrangle(r, (r-1)*(r-1));
    
    nmap.SetOppositeEdgesDofSymmetry(r, r-1);
  }

  

  //! displays informations class QuadrangleReference<type>
  template<int type> ostream& operator <<(ostream& out, const QuadrangleReference<type>& e)
  {
    out<<static_cast<const ElementReference<Dimension2, type>&>(e);
    return out;
  }
}
  
#define MONTJOIE_FILE_QUADRANGLE_REFERENCE_CXX
#endif
