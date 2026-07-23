#ifndef MONTJOIE_FILE_TRIANGLE_REFERENCE_CXX

namespace Montjoie
{
  
  //! default constructor
  template<int type>
  TriangleReference<type>::TriangleReference() : ElementReference<Dimension2, type>(Fb_geom)
  {
    type_interpolation = TriangleGeomReference::LOBATTO_BASIS;
  }
  

  //! returns the size of memory used by the object
  template<int type>
  size_t TriangleReference<type>::GetMemorySize() const
  {
    size_t taille = ElementReference<Dimension2, type>::GetMemorySize();
    taille += Fb_geom.GetMemorySize();
    return taille;
  }


  //! returns a reference to the object storing shape functions
  template<int type>
  TriangleGeomReference& TriangleReference<type>::GetGeometricElement()
  {
    return Fb_geom;
  }
  
  
  //! returns a reference to the object storing shape functions
  template<int type>
  const TriangleGeomReference& TriangleReference<type>::GetGeometricElement() const
  {
    return Fb_geom;
  }


  //! method useful to number mesh
  template<int type>
  void TriangleReference<type>::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    nmap.SetFormulationDG(dg);
    if (dg == ElementReference_Base::DISCONTINUOUS)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, this->nb_dof_loc);
	
	return;
      }
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, this->order+1);
	nmap.SetNbDofTriangle(this->order, 0);
	
	nmap.SetOppositeEdgesDofSymmetry(this->order, this->order+1);
	    
	return;
      }
    
    nmap.SetNbDofVertex(this->order, 1);
    nmap.SetNbDofEdge(this->order, this->order-1);
    nmap.SetNbDofQuadrangle(this->order, (this->order-1)*(this->order-1));
    nmap.SetNbDofTriangle(this->order, (this->order-1)*(this->order-2)/2);
    
    nmap.SetOppositeEdgesDofSymmetry(this->order, this->order-1);
  }
  
  
  //! construction of finite element
  template<int type>
  void TriangleReference<type>::ConstructFiniteElement(int r, int rgeom, int rquad,
                                                       int type_quad, int rsurf, int type_surf)
  {
    if (rgeom == 0)
      rgeom = r;
    
    if (rquad == 0)
      rquad = r;
    
    if (type_quad == -1)
      type_quad = TriangleQuadrature::QUADRATURE_GAUSS;
    
    if (rsurf == 0)
      rsurf = rquad;
    
    if (type_surf == -1)
      type_surf = Globatto<Real_wp>::QUADRATURE_GAUSS;
    
    this->order = r;
    this->SetGeometryOrder(rgeom);
    this->order_quad = rquad;

    // 1-D finite element 
    ConstructFiniteElement1D(r, r, rquad, type_surf);
    
    // quadrature rules
    ConstructQuadrature(rquad, type_quad, rsurf, type_surf);
    
    // orthogonal functions
    Fb_geom.ComputeOrthogonalFunctions(max(r, rgeom));

    // nodal shape functions
    Fb_geom.ConstructNodalShapeFunctions(rgeom);
    this->EdgesNodal = this->GetNodalNumber();
    
  } 


  //! constructs 1-D finite element (restriction of the 2-D finite element to an edge)
  template<int type>
  void TriangleReference<type>
  ::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {
    if (this->element_surface != NULL)
      delete this->element_surface;

    EdgeGaussReference* edge = new EdgeGaussReference();
    edge->ConstructFiniteElement(r, rgeom, rquad, type_quad);    
    this->element_surface = edge;    
  }

  
  //! construction of the quadrature formulas over a triangle
  /*!
    \param[in] r order of quadrature
  */
  template<int type>
  void TriangleReference<type>::ConstructQuadrature(int r, int type_quad, int rsurf, int type_surf)
  {
    // quadrature points on triangles
    VectR2 points2d; VectReal_wp weights2d;
    TriangleQuadrature::ConstructQuadrature(2*r, points2d, weights2d, type_quad); 
    int N = points2d.GetM();
    
    // quadrature points for boundary integrals
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(rsurf, type_surf);
    this->SetPoints1D(gauss.Points());
    this->SetWeights1D(gauss.Weights());
        
    points2d.Resize(N + 3*(rsurf+1));
    this->num_quad_points_surf.Reallocate(3);
    int nb = N;
    this->num_quad_points_surf(0).Reallocate(rsurf+1);
    // Points of integration for the first edge
    for (int i = 0; i <= rsurf; i++)
      {
	points2d(nb).Init(this->Points1D(i), 0);
	this->num_quad_points_surf(0)(i) = nb++;
      }
    
    // Points of integration for the second edge
    this->num_quad_points_surf(1).Reallocate(rsurf+1);
    for (int i = 0; i <= rsurf; i++)
      {
	points2d(nb).Init(Real_wp(1)-this->Points1D(i), this->Points1D(i));
	this->num_quad_points_surf(1)(i) = nb++;
      }
    
    // Points of integration for the third edge
    this->num_quad_points_surf(2).Reallocate(rsurf+1);
    for (int i = 0; i <= rsurf; i++)
      {
	points2d(nb).Init(0, Real_wp(1)-this->Points1D(i));
	this->num_quad_points_surf(2)(i) = nb++;    
      }
    
    this->nb_points_quadrature_inside = N;    
    this->nb_points_quadrature_boundaries = 3*(rsurf+1);
    
    this->offset_faceSh.Reallocate(4); this->offset_faceSh.Fill();
    Mlt(rsurf+1, this->offset_faceSh);
    
    for (int k = 0; k < 3; k++)
      {
        this->PoidsFlux(k) = this->Weights1D();
        Mlt(Real_wp(0.5), this->PoidsFlux(k));
      }
    
    this->SetWeightsND(weights2d);
    this->SetPointsND(points2d);
  }


  //! displays informations about class TriangleReference<type>
  template<int type>
  ostream& operator <<(ostream& out, const TriangleReference<type>& e)
  {
    out << static_cast<const ElementReference<Dimension2, type>& >(e);
    return out;
  }
  
}


#define MONTJOIE_FILE_TRIANGLE_REFERENCE_CXX
#endif
