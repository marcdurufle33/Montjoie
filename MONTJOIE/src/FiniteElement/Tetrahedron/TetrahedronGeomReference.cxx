#ifndef MONTJOIE_FILE_TETRAHEDRON_GEOM_REFERENCE_CXX

namespace Montjoie
{
  
  /****************************
   * Initialization functions *
   ****************************/
  

  //! default constructor
  TetrahedronGeomReference::TetrahedronGeomReference() : ElementGeomReference<Dimension3>()
  {
    type_interpolation_nodal = LOBATTO_BASIS;
    
    // normales of the four faces
    this->normale.Reallocate(4);
    this->normale(0).Init(0,0,-1);
    this->normale(1).Init(0,-1,0);
    this->normale(2).Init(-1,0,0);
    this->normale(3).Init(1,1,1);    

    this->tangente_loc_x.Reallocate(4);
    this->tangente_loc_y.Reallocate(4);
    this->tangente_loc_x(0).Init(1,0,0);
    this->tangente_loc_y(0).Init(0,1,0);
    this->tangente_loc_x(1).Init(1,0,0);
    this->tangente_loc_y(1).Init(0,0,1);
    this->tangente_loc_x(2).Init(0,1,0);
    this->tangente_loc_y(2).Init(0,0,1);
    this->tangente_loc_x(3).Init(-1,1,0);
    this->tangente_loc_y(3).Init(-1,0,1);
    
    this->hybrid_type_elt = 0;
    this->nb_vertices_elt = 4;
    this->nb_edges_elt = 6;
    this->nb_boundaries_elt = 4;
    
    this->is_local_face_quad.Reallocate(4);
    this->is_local_face_quad.Fill(false); 
  }
  

  size_t TetrahedronGeomReference::GetMemorySize() const
  {
    size_t taille = ElementGeomReference<Dimension3>::GetMemorySize();
    taille += function_tri.GetMemorySize();
    taille += coefFi_curve.GetMemorySize();
    taille += EdgesNodal.GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += Seldon::GetMemorySize(OddJacobiPolynom);
    taille += Seldon::GetMemorySize(EvenJacobiPolynom);
    taille += NumOrtho3D.GetMemorySize();
    taille += InvWeightPolynomial.GetMemorySize();
    taille += CoefLegendre.GetMemorySize();
    taille += CoefOddJacobi.GetMemorySize();
    taille += CoefEvenJacobi.GetMemorySize();
    taille += InverseVDM.GetMemorySize();
    return taille;
  }


  //! constructing finite element
  void TetrahedronGeomReference
  ::ConstructFiniteElement(int rgeom)
  {
    this->order_geom = rgeom;
    
    // orthogonal functions
    ComputeOrthogonalFunctions(rgeom);
    
    // shape functions (used for Fi)
    ConstructNodalShapeFunctions(rgeom);
  }
  
    
  //! construction of nodal points, which are used for interpolation  
  void TetrahedronGeomReference::ConstructNodalShapeFunctions(int r)
  {
    // clearing previous points
    this->points_nodal1d.Clear(); this->points_nodal2d_quad.Clear();
    this->points_nodal2d_tri.Clear();
    
    MeshNumbering<Dimension2>::
      ConstructTriangularNumbering(r, this->NumNodes2D_tri, this->CoordinateNodes2D_tri);
    
    MeshNumbering<Dimension2>::
      ConstructQuadrilateralNumbering(r, this->NumNodes2D_quad, this->CoordinateNodes2D_quad);

    // numbering on 3-D tetrahedron
    MeshNumbering<Dimension3>::
      ConstructTetrahedralNumbering(r, this->NumNodes3D, this->CoordinateNodes);
    
    // nodal shape functions on triangle
    function_tri.ConstructFiniteElement(r);
    
    ConstructLobattoPoints(r, this->points_nodal1d, this->points_nodal2d_tri,
                           this->points_nodal_nd);    
    
    InvWeightPolynomial.Reallocate(this->points_nodal_nd.GetM());
    NumOrtho3D.Reallocate(r+1, r+1, r+1); NumOrtho3D.Fill(-1);
    int nb = 0;
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	for (int k = 0; k <= r-i-j; k++)
	  {
	    NumOrtho3D(i, j, k) = nb;
	    InvWeightPolynomial(nb) = CoefLegendre(i)*CoefOddJacobi(i, j)*CoefEvenJacobi(i+j, k);
	    nb++;
	  }
    
    ComputeLagrangianFunctions(r);
    
    // nodes on edges
    EdgesNodal.Reallocate(r+1, 6);
    for (int i = 0; i < 6; i++)
      {
	EdgesNodal(0, i) = MeshNumbering<Dimension3>::FirstExtremityEdge(0, i);
	EdgesNodal(r, i) = MeshNumbering<Dimension3>::SecondExtremityEdge(0, i);
	for (int k = 1; k < r; k++)
	  EdgesNodal(k, i) = i*(r-1) + 3 + k;
      }
    
    this->FacesNodal.Reallocate(this->points_nodal2d_tri.GetM(), 4);
    // nodes on faces
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  int node = this->NumNodes2D_tri(i,j);
	  this->FacesNodal(node,0) = this->NumNodes3D(i, j, 0);
	  this->FacesNodal(node,1) = this->NumNodes3D(i, 0, j);
	  this->FacesNodal(node,2) = this->NumNodes3D(0, i, j);
	  this->FacesNodal(node,3) = this->NumNodes3D(r-i-j, i, j);
	}
    
    if (r > 2)
      {
        VectR2 points_nodal2d_inside_tri((r-1)*(r-2)/2);
        for (int i = 0; i < points_nodal2d_inside_tri.GetM(); i++)
          points_nodal2d_inside_tri(i) = this->points_nodal2d_tri(i+3*r);
        
        MeshNumbering<Dimension3>::
          GetRotationTriangularFace(points_nodal2d_inside_tri, this->FacesNodal_Rotation_Tri);
      }
    
    VectR3 grad_phi;
    this->GradientPhi_Nodal.Reallocate(this->points_nodal_nd.GetM(), this->points_nodal_nd.GetM());
    for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
      {
	ComputeGradientPhiNodalRef(this->points_nodal_nd(j), grad_phi);
	
        for (int k = 0; k < this->points_nodal_nd.GetM(); k++)
	  this->GradientPhi_Nodal(k, j) = grad_phi(k);
      }
    
    ComputeCurvedTransformation();
    ComputeCoefficientTransformation();
  }
  
  
  //! nodal points are set on Pascal's tetrahedron, regularly spaced points  
  void TetrahedronGeomReference::
  ConstructRegularPoints(int r, VectReal_wp& points_nodal1d_, VectR2& points_nodal2d_tri_,
                         const Matrix<int>& NumNodes_tri,
                         VectR3& points_nodal3d_, const Array3D<int>& NumNodes)   
  {
    // number of nodal points on a triangle : (r+1)*(r+2)/2
    int nb_points_nodal_tri = (r+2)*(r+1)/2;
    
    // number of nodal points on a tetrahedron : (r+1)(r+2)(r+3)/6
    int nb_points_nodal_elt = (r+1)*(r+2)*(r+3)/6;
    
    // we construct equally spaced points on the triangle
    points_nodal1d_.Reallocate(r+1);
    for (int i = 0; i <= r; i++)
      points_nodal1d_(i) = Real_wp(i)/r;
        
    // 2-D nodal points
    points_nodal2d_tri_.Reallocate(nb_points_nodal_tri);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	points_nodal2d_tri_(NumNodes_tri(i,j)).Init(points_nodal1d_(i), points_nodal1d_(j));
        
    points_nodal3d_.Reallocate(nb_points_nodal_elt);        
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	for (int k = 0; k <= r-i-j; k++)
	  points_nodal3d_(NumNodes(i,j,k)).Init(points_nodal1d_(i),
                                                points_nodal1d_(j), points_nodal1d_(k));
    
  }
  
  
  //! Construction of Hesthaven's points  
  void TetrahedronGeomReference::
  ConstructLobattoPoints(int r, VectReal_wp& points_nodal1d_, VectR2& points_nodal2d_tri_,
                         VectR3& points_nodal3d_)
  {
    points_nodal3d_.Reallocate((r+1)*(r+2)*(r+3)/6);
    
    // points associated with vertices
    points_nodal3d_(0) = R3(Real_wp(0), Real_wp(0), Real_wp(0));
    points_nodal3d_(1) = R3(Real_wp(1), Real_wp(0), Real_wp(0));
    points_nodal3d_(2) = R3(Real_wp(0), Real_wp(1), Real_wp(0));
    points_nodal3d_(3) = R3(Real_wp(0), Real_wp(0), Real_wp(1));
    
    TriangleGeomReference::ConstructLobattoPoints(r, 0, points_nodal1d_,
						  points_nodal2d_tri_);
    
    int ind = 4; R3 point;
    for (int ne = 0; ne < 6; ne++)
      {
	int n1 = MeshNumbering<Dimension3>::FirstExtremityEdge(0, ne);
	int n2 = MeshNumbering<Dimension3>::SecondExtremityEdge(0, ne);
	
	for (int k = 1; k < r; k++)
	  {
	    point.Fill(0);
	    Add(Real_wp(1)-points_nodal1d_(k), points_nodal3d_(n1), point);
	    Add(points_nodal1d_(k), points_nodal3d_(n2), point);
	    points_nodal3d_(ind) = point; ind++;
	  }
      }
    
    // points associated with triangular faces
    int nb_points_nodal_inside_tri = 0;
    VectR2 points_nodal2d_inside_tri;
    if (r > 2)
      {
	points_nodal2d_inside_tri.Reallocate((r-1)*(r-2)/2);
        int node = 3*r;
	for (int j = 1; j < r; j++)
	  for (int k = 1; k < r-j; k++)
	    {
	      points_nodal2d_inside_tri(nb_points_nodal_inside_tri) = points_nodal2d_tri_(node);
	      nb_points_nodal_inside_tri++; node++;
	    }
      }
    
    // first face z = 0
    for (int i = 0;  i < nb_points_nodal_inside_tri; i++)
      points_nodal3d_(ind++).Init(points_nodal2d_inside_tri(i)(0),
                                  points_nodal2d_inside_tri(i)(1), Real_wp(0));
    
    // second face y = 0
    for (int i = 0;  i < nb_points_nodal_inside_tri; i++)
      points_nodal3d_(ind++).Init(points_nodal2d_inside_tri(i)(0),
                                  Real_wp(0), points_nodal2d_inside_tri(i)(1));
    
    // third face z = 0
    for (int i = 0;  i < nb_points_nodal_inside_tri; i++)
      points_nodal3d_(ind++).Init(Real_wp(0), points_nodal2d_inside_tri(i)(0),
                                  points_nodal2d_inside_tri(i)(1));
    
    // fourth face z = 0
    for (int i = 0;  i < nb_points_nodal_inside_tri; i++)
      points_nodal3d_(ind++)
        .Init(Real_wp(1)-points_nodal2d_inside_tri(i)(0) - points_nodal2d_inside_tri(i)(1),
              points_nodal2d_inside_tri(i)(0), points_nodal2d_inside_tri(i)(1));
    
    R3 center(0.25, 0.25, 0.25);
    int n1 = 0, n4 = 0, n6 = 0, n12 = 0, n24 = 0;
    ind = 2*r*r+2; 
    // inside points
    switch (r)
      {
      case 1 :
      case 2 :
      case 3 :
	{
	  // no inside node
	}
	break;
      case 4 :
	{
	  n1 = 1;
	  points_nodal3d_(ind) = center;
	}
	break;
      case 5 :
	{
	  n4 = 1;
	  points_nodal3d_(ind).Init(0.1834903473, 0.1834903473, 0.1834903473);
	}
	break;
      case 6 :
	{
	  n4 = 1; n6 = 1;
	  points_nodal3d_(ind).Init(0.1402705081, 0.1402705081, 0.1402705081);
	  points_nodal3d_(ind+4).Init(0.1457947417, 0.1457947417, 0.3542052583);
	}
	break;
      case 7 :
	{
	  n4 = 2; n12 = 1;
	  points_nodal3d_(ind).Init(0.1144606542, 0.1144606542, 0.1144606542);
	  points_nodal3d_(ind+4).Init(0.2917002822, 0.2917002822, 0.2917002822);
	  points_nodal3d_(ind+8).Init(0.1208429970, 0.1208429970, 0.2812936703);
	}
	break;
      case 8 :
	{
	  n1 = 1; n4 = 1; n6 = 1; n12 = 2;
	  points_nodal3d_(ind) = center;
	  points_nodal3d_(ind+1).Init(0.0991203900, 0.0991203900, 0.0991203900);
	  points_nodal3d_(ind+5).Init(0.1079468963, 0.1079468963, 0.3920531037);
	  points_nodal3d_(ind+11).Init(0.1047451941, 0.1047451941, 0.2323149656);
	  points_nodal3d_(ind+23).Init(0.2419418605, 0.2419418605, 0.1099065340);
	}
	break;
      case 9 :
	{
	  n4 = 2; n12 = 4;
	  points_nodal3d_(ind).Init(0.0823287303, 0.0823287303, 0.0823287303);
	  points_nodal3d_(ind+4).Init(0.2123055477, 0.2123055477, 0.2123055477);
	  points_nodal3d_(ind+8).Init(0.0873980781, 0.0873980781, 0.1954980564);
	  points_nodal3d_(ind+20).Init(0.0916714679, 0.0916714679, 0.3347047619);
	  points_nodal3d_(ind+32).Init(0.2040338880, 0.2040338880, 0.0923029247);
	  points_nodal3d_(ind+44).Init(0.3483881173, 0.3483881173, 0.0956734931);
	}
	break;
      case 10 :
	{
	  n4 = 3; n6 = 2; n12 = 3; n24 = 1;
	  points_nodal3d_(ind).Init(0.0678316144, 0.0678316144, 0.0678316144);
	  points_nodal3d_(ind+4).Init(0.1805746957, 0.1805746957, 0.1805746957);
	  points_nodal3d_(ind+8).Init(0.3051527124, 0.3051527124, 0.3051527124);
	  points_nodal3d_(ind+12).Init(0.1835663764, 0.1835663764, 0.3164336236);
	  points_nodal3d_(ind+18).Init(0.0780456199, 0.0780456199, 0.4219543801);
	  points_nodal3d_(ind+24).Init(0.1734555313, 0.1734555313, 0.0768158196);
	  points_nodal3d_(ind+36).Init(0.0724033935, 0.0724033935, 0.1658367169);
	  points_nodal3d_(ind+48).Init(0.0768451848, 0.0768451848, 0.2889363346);
	  points_nodal3d_(ind+60).Init(0.0809350066, 0.1776946096, 0.3019928872);
	}
	break;
      default :
	{
	  cout << "order not implemented " << endl;
	  abort();
	}
      }
    
    int nb_couples = n1 + n4 + n6 + n12 + n24;
    IVect NbPoints_couple(nb_couples);
    for (int i = 0; i < n1; i++)
      NbPoints_couple(i) = 1;
    
    for (int i = 0; i < n4; i++)
      NbPoints_couple(n1+i) = 4;

    for (int i = 0; i < n6; i++)
      NbPoints_couple(n1+n4+i) = 6;
    
    for (int i = 0; i < n12; i++)
      NbPoints_couple(n1+n4+n6+i) = 12;
    
    for (int i = 0; i < n24; i++)
      NbPoints_couple(n1+n4+n6+n12+i) = 24;
    
    Real_wp x, y, z;
    // now we generate the other points by symmetry
    for (int i = 0; i < nb_couples; i++)
      {
	int np = NbPoints_couple(i);
	
	if (np==4)
	  {
	    x = points_nodal3d_(ind)(0);
	    y = Real_wp(1) - Real_wp(3)*x;
	    points_nodal3d_(ind).Init(x, x, y);
	    points_nodal3d_(ind+1).Init(x, x, x);
	    points_nodal3d_(ind+2).Init(x, y, x);
	    points_nodal3d_(ind+3).Init(y, x, x);
	  }
	else if (np==6)
	  {
	    x = points_nodal3d_(ind)(0);
	    y = 0.5*(Real_wp(1) - Real_wp(2)*x);
	    points_nodal3d_(ind).Init(x, x, y);
	    points_nodal3d_(ind+1).Init(x, y, x);
	    points_nodal3d_(ind+2).Init(y, x, x);
	    points_nodal3d_(ind+3).Init(y, y, x);
	    points_nodal3d_(ind+4).Init(y, x, y);
	    points_nodal3d_(ind+5).Init(x, y, y);
	  }
	else if (np==12)
	  {
	    x = points_nodal3d_(ind)(0);
	    y = points_nodal3d_(ind)(2);
	    z = Real_wp(1) - Real_wp(2)*x -y;
	    points_nodal3d_(ind).Init(x, x, y);
	    points_nodal3d_(ind+1).Init(x, y, x);
	    points_nodal3d_(ind+2).Init(y, x, x);
	    points_nodal3d_(ind+3).Init(x, x, z);
	    points_nodal3d_(ind+4).Init(x, z, x);
	    points_nodal3d_(ind+5).Init(z, x, x);
	    points_nodal3d_(ind+6).Init(x, y, z);
	    points_nodal3d_(ind+7).Init(x, z, y);
	    points_nodal3d_(ind+8).Init(y, x, z);
	    points_nodal3d_(ind+9).Init(z, x, y);
	    points_nodal3d_(ind+10).Init(y, z, x);
	    points_nodal3d_(ind+11).Init(z, y, x);
	  }
	else if (np==24)
	  {
	    Real_wp l0 = points_nodal3d_(ind)(0);
	    Real_wp l1 = points_nodal3d_(ind)(1);
	    Real_wp l2 = points_nodal3d_(ind)(2);
	    Real_wp l3 = Real_wp(1) - l1 - l2 - l0;
	    points_nodal3d_(ind+1).Init(l1, l2, l3);
	    points_nodal3d_(ind+2).Init(l1, l3, l2);
	    points_nodal3d_(ind+3).Init(l2, l1, l3);
	    points_nodal3d_(ind+4).Init(l2, l3, l1);
	    points_nodal3d_(ind+5).Init(l3, l1, l2);
	    points_nodal3d_(ind+6).Init(l3, l2, l1);
	    points_nodal3d_(ind+7).Init(l1, l0, l3);
	    points_nodal3d_(ind+8).Init(l1, l3, l0);
	    points_nodal3d_(ind+9).Init(l0, l3, l1);
	    points_nodal3d_(ind+10).Init(l0, l1, l3);
	    points_nodal3d_(ind+11).Init(l3, l1, l0);
	    points_nodal3d_(ind+12).Init(l3, l0, l1);
	    points_nodal3d_(ind+13).Init(l0, l2, l3);
	    points_nodal3d_(ind+14).Init(l0, l3, l2);
	    points_nodal3d_(ind+15).Init(l2, l0, l3);
	    points_nodal3d_(ind+16).Init(l2, l3, l0);
	    points_nodal3d_(ind+17).Init(l3, l0, l2);
	    points_nodal3d_(ind+18).Init(l3, l2, l0);
	    points_nodal3d_(ind+19).Init(l1, l0, l2);
	    points_nodal3d_(ind+20).Init(l1, l2, l0);
	    points_nodal3d_(ind+21).Init(l0, l2, l1);
	    points_nodal3d_(ind+22).Init(l2, l1, l0);
	    points_nodal3d_(ind+23).Init(l2, l0, l1);
	  }
	
	ind += np;
      }    
  }
  

  //! constructing orthogonal functions    
  void TetrahedronGeomReference::ComputeOrthogonalFunctions(int r)
  {
    // generating Jacobi polynomials P_m^{2i,0}
    GetJacobiPolynomial(LegendrePolynom, r, Real_wp(0), Real_wp(0));
    OddJacobiPolynom.Reallocate(r+1);
    EvenJacobiPolynom.Reallocate(r+1);
    for (int i = 0; i <= r; i++)
      {
	GetJacobiPolynomial(OddJacobiPolynom(i), r+1-i, Real_wp(2*i+1), Real_wp(0));
	GetJacobiPolynomial(EvenJacobiPolynom(i), r+1-i, Real_wp(2*i+2), Real_wp(0));
      }
    
    // we compute weight of polynomials in order to obtain orthonormal expansion
    CoefLegendre.Reallocate(r+1); CoefLegendre.Fill(0);
    CoefOddJacobi.Reallocate(r+1, r+1); CoefOddJacobi.Fill(0);
    CoefEvenJacobi.Reallocate(r+1, r+1); CoefEvenJacobi.Fill(0);
    VectReal_wp xi, omega;
    ComputeGaussLegendre(xi, omega, r+1);
    for (int q = 0; q <= r+1; q++)
      {
	VectReal_wp Pn;
	EvaluateJacobiPolynomial(LegendrePolynom, r, 2.0*xi(q) - 1.0, Pn);
	for (int i = 0; i <= r; i++)
	  CoefLegendre(i) += omega(q)*Pn(i)*Pn(i);
	
	for (int i = 0; i <= r; i++)
	  {
	    EvaluateJacobiPolynomial(OddJacobiPolynom(i), r-i, 2.0*xi(q) - 1.0, Pn);
	    for (int j = 0; j <= r-i; j++)
	      CoefOddJacobi(i, j) += omega(q)*Pn(j)*Pn(j)*pow(1.0-xi(q), 2*i+1); 
	  }

	for (int i = 0; i <= r; i++)
	  {
	    EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, 2.0*xi(q) - 1.0, Pn);
	    for (int j = 0; j <= r-i; j++)
	      CoefEvenJacobi(i, j) += omega(q)*Pn(j)*Pn(j)*pow(1.0-xi(q), 2*i+2); 
	  }
      }
    
    // we keep 1/sqrt(omega)
    for (int i = 0; i <= r; i++)
      {
	CoefLegendre(i) = 1.0/sqrt(CoefLegendre(i));
	for (int j = 0; j <= r-i; j++)
	  {
	    CoefOddJacobi(i, j) = 1.0/sqrt(CoefOddJacobi(i, j));
	    CoefEvenJacobi(i, j) = 1.0/sqrt(CoefEvenJacobi(i, j));
	  }
      }
  }
  
  
  //! computation of shape functions  
  void TetrahedronGeomReference::ComputeLagrangianFunctions(int r)
  {
    int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    InverseVDM.Reallocate(nb_points_nodal_elt, nb_points_nodal_elt); 
    InverseVDM.Fill(Real_wp(0));
            
    VectReal_wp phi;
    for (int k = 0; k < nb_points_nodal_elt; k++)
      {
	ComputeValuesPhiOrthoRef(r, NumOrtho3D, InvWeightPolynomial, this->points_nodal_nd(k), phi);
	for (int j = 0; j < nb_points_nodal_elt; j++)
	  InverseVDM(j, k) = phi(j);
      }
    
    GetInverse(InverseVDM);
  }

  
  //! computation of coefficients to compute quickly Fi and DFi  
  void TetrahedronGeomReference::ComputeCurvedTransformation()
  {
    int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    
    coefFi_curve.Reallocate(nb_points_nodal_elt, nb_points_nodal_elt);
    coefFi_curve.Fill(0); VectReal_wp val_phi;
    // use of transformation Fi defined in Solin book for inside points
    TinyVector<Real_wp, 4> lambda; Real_wp zeta, coef, zetas, coefs, l1, l2, l3;
    for (int i = 2*this->order_geom*this->order_geom+2; i < nb_points_nodal_elt; i++)
      {
        lambda(1) = this->points_nodal_nd(i)(0);
	lambda(2) = this->points_nodal_nd(i)(1);
	lambda(3) = this->points_nodal_nd(i)(2);
	lambda(0) = 1.0 - lambda(1) - lambda(2) -lambda(3);
	
	// part due to vertices
	for (int k = 0; k < 4; k++)
	  coefFi_curve(i, k) = lambda(k);
	
	// part due to edges
	for (int j = 0; j < 6; j++)
	  {
	    int n1 = MeshNumbering<Dimension3>::FirstExtremityEdge(0,j);
	    int n2 = MeshNumbering<Dimension3>::SecondExtremityEdge(0,j);
	    zeta = 0.5*(lambda(n2)-lambda(n1)) + 0.5; 
	    coef = lambda(n1)*lambda(n2)/(zeta*(1.0-zeta));
	    // contribution of points of the edge
	    function_tri.ComputeValuesNodalPhi1D(zeta, val_phi);
	    for (int k = 0; k <= this->order_geom; k++)
	      coefFi_curve(i, EdgesNodal(k, j) ) += val_phi(k)*coef;
	
	    // contribution of the extremities of the edge
	    coefFi_curve(i, EdgesNodal(0,j) ) -= (1.0-zeta)*coef;
	    coefFi_curve(i, EdgesNodal(this->order_geom,j) ) -= zeta*coef;
	  }
	
	// part due to faces
	for (int j = 0; j < 4; j++)
	  {
	    int n1 = MeshNumbering<Dimension3>::FirstExtremityFace(0,j);
	    int n2 = MeshNumbering<Dimension3>::SecondExtremityFace(0,j);
	    int n3 = MeshNumbering<Dimension3>::ThirdExtremityFace(0,j);
	    
	    l2 = (2.0*lambda(n2)-lambda(n1)-lambda(n3) + 1.0)/3;
	    l3 = (2.0*lambda(n3)-lambda(n1)-lambda(n2) + 1.0)/3;
	    l1 = 1.0 - l2 - l3;
	    coef = lambda(n1)*lambda(n2)*lambda(n3)/(l1*l2*l3);
	    
	    // part due to first edge of the face
	    int ne = MeshNumbering<Dimension3>::FirstEdgeFace(0, j);
	    zetas = 0.5*(l2-l1) + 0.5; 
	    coefs = l1*l2/(zetas*(1.0-zetas));
	    function_tri.ComputeValuesNodalPhi1D(zetas, val_phi);
	    if (MeshNumbering<Dimension3>::FirstExtremityEdge(0,ne) == n1)
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, EdgesNodal(k, ne)) -= coef*coefs*val_phi(k);
	    else
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, EdgesNodal(this->order_geom-k, ne)) -= coef*coefs*val_phi(k);
	    
	    // contribution of the extremities
	    coefFi_curve(i, n1) += (1.0-zetas)*coef*coefs;
	    coefFi_curve(i, n2) += zetas*coef*coefs;

	    // part due to second edge of the face
	    ne = MeshNumbering<Dimension3>::SecondEdgeFace(0, j);
	    zetas = 0.5*(l3-l2) + 0.5;
	    coefs = l2*l3/(zetas*(1.0-zetas));
	    function_tri.ComputeValuesNodalPhi1D(zetas, val_phi);
	    if (MeshNumbering<Dimension3>::FirstExtremityEdge(0,ne) == n2)
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, EdgesNodal(k, ne)) -= coef*coefs*val_phi(k);
	    else
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, EdgesNodal(this->order_geom-k, ne)) -= coef*coefs*val_phi(k);
	    
	    // contribution of the extremities
	    coefFi_curve(i, n2) += (1.0-zetas)*coef*coefs;
	    coefFi_curve(i, n3) += zetas*coef*coefs;
	    
	    // part due to third edge of the face
	    ne = MeshNumbering<Dimension3>::ThirdEdgeFace(0, j);
	    zetas = 0.5*(l1-l3) + 0.5;
	    coefs = l3*l1/(zetas*(1.0-zetas));
	    function_tri.ComputeValuesNodalPhi1D(zetas, val_phi);
	    if (MeshNumbering<Dimension3>::FirstExtremityEdge(0,ne) == n3)
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, EdgesNodal(k, ne)) -= coef*coefs*val_phi(k);
	    else
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, EdgesNodal(this->order_geom-k, ne)) -= coef*coefs*val_phi(k);
	    
	    // contribution of the extremities
	    coefFi_curve(i, n3) += (1.0-zetas)*coef*coefs;
	    coefFi_curve(i, n1) += zetas*coef*coefs;
	    
	    // contribution of the vertices of the face
	    coefFi_curve(i, n1) -= l1*coef;
	    coefFi_curve(i, n2) -= l2*coef;
	    coefFi_curve(i, n3) -= l3*coef;
	    
	    // contribution of the parametrization of the face
	    function_tri.ComputeValuesPhiNodalRef(R2(l2, l3), val_phi);
	    for (int k = 0; k < this->points_nodal2d_tri.GetM(); k++)
	      coefFi_curve(i, this->FacesNodal(k, j)) += coef*val_phi(k);
	    
	  }
      }

    const TriangleGeomReference& tri = this->GetTriangularSurfaceFiniteElement();
    this->InitTriangularCurvature(tri);
  }
  

  //! computation of coefficients to compute quickly Fi and DFi    
  void TetrahedronGeomReference::ComputeCoefficientTransformation()
  {
    ComputeCurvedTransformation();
    VectR3 all_points = this->points_nodal_nd;
    int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    all_points.PushBack(this->points3d); all_points.PushBack(this->points_dof3d);
    
    int nb_points_all = all_points.GetM();
    coefFi_curve.Resize(nb_points_all, nb_points_nodal_elt);
    
    // for other points, we use lagrangian interpolation
    VectReal_wp val_phi;
    for (int i = nb_points_nodal_elt; i < nb_points_all; i++)
      {
	ComputeValuesPhiNodalRef(all_points(i), val_phi);
	for (int j = 0; j < nb_points_nodal_elt; j++)
	  coefFi_curve(i, j) = val_phi(j);
      }
    
  }


  void TetrahedronGeomReference
  ::GetGradient3D_FromGradient2D(int num_loc, int i, const VectReal_wp& nabla_nx,
				 const VectReal_wp& nabla_ny, const VectReal_wp& nabla_nz,
				 R3& d_nx, R3& d_ny, R3& d_nz) const
  {
    switch (num_loc)
      {
      case 0:
	d_nx(0) = nabla_nx(2*i); d_nx(1) = nabla_nx(2*i+1); d_nx(2) = 0;
	d_ny(0) = nabla_ny(2*i); d_ny(1) = nabla_ny(2*i+1); d_ny(2) = 0;
	d_nz(0) = nabla_nz(2*i); d_nz(1) = nabla_nz(2*i+1); d_nz(2) = 0;
	break;
      case 1:
	d_nx(0) = nabla_nx(2*i); d_nx(1) = 0; d_nx(2) = nabla_nx(2*i+1);
	d_ny(0) = nabla_ny(2*i); d_ny(1) = 0; d_ny(2) = nabla_ny(2*i+1);
	d_nz(0) = nabla_nz(2*i); d_nz(1) = 0; d_nz(2) = nabla_nz(2*i+1);
	break;
      case 2:
	d_nx(0) = 0; d_nx(1) = nabla_nx(2*i); d_nx(2) = nabla_nx(2*i+1);
	d_ny(0) = 0; d_ny(1) = nabla_ny(2*i); d_ny(2) = nabla_ny(2*i+1);
	d_nz(0) = 0; d_nz(1) = nabla_nz(2*i); d_nz(2) = nabla_nz(2*i+1);
	break;
      case 3:
	// for the last face,  the transformation from unit triangle
	// to the face reads :
	// x = (1- xt - yt)
	// y = xt
	// z = yt
	// the jacobian matrix is equal to [-1 -1; 1 0; 0 1]
	// DF^*{-1} is equal to 1/3 [-1 -1; 2 -1; -1 2]
	// we use the rule \nabla n_i = DF^*{-1} \hat{\nabla} n_i
	d_nx(0) = (-nabla_nx(2*i)-nabla_nx(2*i+1))/3;
	d_ny(0) = (-nabla_ny(2*i)-nabla_ny(2*i+1))/3;
	d_nz(0) = (-nabla_nz(2*i)-nabla_nz(2*i+1))/3;
        
	d_nx(1) = (2.0*nabla_nx(2*i)-nabla_nx(2*i+1))/3;
	d_ny(1) = (2.0*nabla_ny(2*i)-nabla_ny(2*i+1))/3;
	d_nz(1) = (2.0*nabla_nz(2*i)-nabla_nz(2*i+1))/3;
        
	d_nx(2) = (-nabla_nx(2*i)+2.0*nabla_nx(2*i+1))/3;
	d_ny(2) = (-nabla_ny(2*i)+2.0*nabla_ny(2*i+1))/3;
	d_nz(2) = (-nabla_nz(2*i)+2.0*nabla_nz(2*i+1))/3;
      }
  }
  
  
  //! returns a tangential vector on a face of the tetrahedron from a 2-D vector on the square/triangle
  R3 TetrahedronGeomReference::GetTangentialVector(int num_loc, const R2& vec_u) const
  {
    switch (num_loc)
      {
      case 0:
	return R3(vec_u(0), vec_u(1), Real_wp(0));
      case 1:
	return R3(vec_u(0), Real_wp(0), vec_u(1));
      case 2:
	return R3(Real_wp(0), vec_u(0), vec_u(1));
      case 3:
	// for the last face,  the transformation from unit triangle
	// to the face reads :
	// x = (1- xt - yt)
	// y = xt
	// z = yt
	// the jacobian matrix is equal to [-1 -1; 1 0; 0 1]
	// DF^*{-1} is equal to 1/3 [-1 -1; 2 -1; -1 2]
	return R3((-vec_u(0)-vec_u(1))/3,
		  (2.0*vec_u(0)-vec_u(1))/3,
		  (-vec_u(0)+2.0*vec_u(1))/3);
      }
    
    return R3();
  }
  
  
  //! transpose operation of GetTangentialVector
  R2 TetrahedronGeomReference::TransposeTangentialVector(int num_loc, const R3& vec_u) const
  {
    switch (num_loc)
      {
      case 0:
	return R2(vec_u(0), vec_u(1));
      case 1:
	return R2(vec_u(0), vec_u(2));
      case 2:
	return R2(vec_u(1), vec_u(2));
      case 3:
        return R2((-vec_u(0)+2.0*vec_u(1)-vec_u(2))/3,
                  (-vec_u(0)-vec_u(1)+2.0*vec_u(2))/3);
        
      }

    return R2();
  }


  //! transpose operation of GetTangentialVector
  R2_Complex_wp TetrahedronGeomReference::TransposeTangentialVector(int num_loc, const R3_Complex_wp& vec_u) const
  {
    switch (num_loc)
      {
      case 0:
	return R2_Complex_wp(vec_u(0), vec_u(1));
      case 1:
	return R2_Complex_wp(vec_u(0), vec_u(2));
      case 2:
	return R2_Complex_wp(vec_u(1), vec_u(2));
      case 3:
        return R2_Complex_wp((-vec_u(0)+2.0*vec_u(1)-vec_u(2))/3,
                  (-vec_u(0)-vec_u(1)+2.0*vec_u(2))/3);
        
      }

    return R2_Complex_wp();
  }


  //! returns the minimal length of the element
  Real_wp TetrahedronGeomReference
  ::GetMinimalSize(const VectR3& s) const
  {
    Real_wp h = Distance(s(0), s(1));
    h = min(h, Distance(s(0), s(2)));
    h = min(h, Distance(s(0), s(3)));
    h = min(h, Distance(s(1), s(2)));
    h = min(h, Distance(s(1), s(3)));
    h = min(h, Distance(s(2), s(3)));
    return h;
  }
  

  //! returns true if the point is outside the unit tetrahedron
  /*!
    \param[in] Xn local coordinates of the point
    \param[in] epsilon threshold 
    \return true if the point is outside
   */
  bool TetrahedronGeomReference::
  OutsideReferenceElement(const VectR3& s, const R3& Xn, const Real_wp& epsilon) const
  {
    Real_wp Rmax = max(abs(s(0)(0)), abs(s(0)(1))); Rmax = max(Rmax, abs(s(0)(2)));
    Rmax = max(Rmax, abs(s(1)(0))); Rmax = max(Rmax, abs(s(1)(1))); Rmax = max(Rmax, abs(s(1)(2)));
    Rmax = max(Rmax, abs(s(2)(0))); Rmax = max(Rmax, abs(s(2)(1))); Rmax = max(Rmax, abs(s(2)(2)));
    Rmax = max(Rmax, abs(s(3)(0))); Rmax = max(Rmax, abs(s(3)(1))); Rmax = max(Rmax, abs(s(3)(2)));
    
    Real_wp dx = s(0).Distance(s(1)), dy = s(0).Distance(s(2)), dz = s(0).Distance(s(3));
    if ((Rmax == Real_wp(0)) || (dx == Real_wp(0)) || (dy == Real_wp(0)) || (dz == Real_wp(0)))
      return true;
    
    if (dx*Xn(0) < -epsilon*Rmax)
      return true;

    if (dy*Xn(1) < -epsilon*Rmax)
      return true;

    if (dz*Xn(2) < -epsilon*Rmax)
      return true;

    if (min(dz, min(dx, dy))*(Xn(0) + Xn(1) + Xn(2) - Real_wp(1)) > epsilon*Rmax)
      return true;
    
    return false;
  }
  
  
  //! returns distance of the local point to the boundary of the element
  /*!
    \param[in] pointloc local coordinates of the point
    \return distance to the boundary (negative if the point is outside)
   */
  Real_wp TetrahedronGeomReference::
  GetDistanceToBoundary(const R3& pointloc) const
  {
    // distance = min(\lambda_1,\lambda_2,\lambda_3,\lambda_4)
    Real_wp distance = pointloc(0); distance = min(distance, pointloc(1));
    distance = min(distance, pointloc(2));
    distance = min(distance, Real_wp(Real_wp(1)-pointloc(0)));
    distance = min(distance, Real_wp(Real_wp(1)-pointloc(1)));
    distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)));    
    distance = min(distance, Real_wp(Real_wp(1)-pointloc(0)-pointloc(1)-pointloc(2)));
    return distance;
  }
  
  
  Real_wp TetrahedronGeomReference::GetDistanceToBoundary(const R3& pointloc, int n) const
  {
    cout << "Not implemented" << endl;
    abort();
    return Real_wp(0);
  }


  //! project the local point on the boundary, if outside the unit tetrahedron  
  int TetrahedronGeomReference::ProjectPointOnBoundary(R3& pointloc) const
  {
    pointloc(0) = max( Real_wp(0), pointloc(0) );
    pointloc(1) = max( Real_wp(0), pointloc(1) );
    pointloc(2) = max( Real_wp(0), pointloc(2) );
    Real_wp one(1);
    if ((pointloc(0)<=pointloc(1))&&(pointloc(0)<=pointloc(2)))
      pointloc(0) = min(pointloc(0), Real_wp(one-pointloc(1)-pointloc(2)));
    else if ((pointloc(1)<=pointloc(0))&&(pointloc(1)<=pointloc(2)))
      pointloc(1) = min(pointloc(1), Real_wp(one-pointloc(0)-pointloc(2)));
    else
      pointloc(2) = min(pointloc(2), Real_wp(one-pointloc(0)-pointloc(1)));
    
    return -1;
  }
  
  
  //! returns 3-D coordinates of a point on a face
  /*!
    \param[in] num_loc local face number in the quad
    \param[in] point_loc coordinate on the face
    \param[out] res local coordinates of the point in the unit square
  */
  void TetrahedronGeomReference::
  GetLocalCoordOnBoundary(int num_loc, const R2& point_loc, R3& res) const
  {
    switch(num_loc)
      {
      case 0:
	res.Init(point_loc(0),point_loc(1),0);break;
      case 1:
	res.Init(point_loc(0),0,point_loc(1));break;
      case 2:
	res.Init(0,point_loc(0),point_loc(1));break;
      case 3:
	res.Init(Real_wp(1)-point_loc(0)-point_loc(1),
		 point_loc(0), point_loc(1));break;
      }
  }


  //! computes 2-D coordinates of a point on an element
  /*!
    \param[in] num_loc local face number in the quad
    \param[out] point_loc coordinate on the face
    \param[int] pt local coordinates of the point in the unit tetrahedron
  */
  void TetrahedronGeomReference::
  GetLocalCoordOnBoundary(int num_loc, const R3& pt, R2& point_loc) const
  {
    switch(num_loc)
      {
      case 0:
	point_loc.Init(pt(0), pt(1)); break;
      case 1:
	point_loc.Init(pt(0), pt(2)); break;
      case 2:
	point_loc.Init(pt(1), pt(2)); break;
      case 3:
	point_loc.Init(pt(1), pt(2)); break;
      }
  }
    
  
  /****************
   * Fj transform *
   ****************/
  
  
  //! transformation Fi in the case of curved tetrahedron
  void TetrahedronGeomReference::
  FjCurve(const SetPoints<Dimension3>& PTReel, const R3& pointloc, R3& res) const
  {
    res.Fill(0); VectReal_wp phi;
    ComputeValuesPhiNodalRef(pointloc, phi);
    for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
      Add(phi(node), PTReel.GetPointNodal(node), res);
  }
  
  
  //! transformation DFi in the case of curved tetrahedron
  void TetrahedronGeomReference::
  DFjCurve(const SetPoints<Dimension3>& PTReel,
	   const R3& pointloc, Matrix3_3& res) const
  {
    R3 grad; VectR3 grad_phi;
    res.Fill(0); ComputeGradientPhiNodalRef(pointloc, grad_phi);
    for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
      {
	grad = grad_phi(node);
	res(0,0) += PTReel.GetPointNodal(node)(0)*grad(0);
	res(1,0) += PTReel.GetPointNodal(node)(1)*grad(0);
	res(2,0) += PTReel.GetPointNodal(node)(2)*grad(0);
	res(0,1) += PTReel.GetPointNodal(node)(0)*grad(1);
	res(1,1) += PTReel.GetPointNodal(node)(1)*grad(1);
	res(2,1) += PTReel.GetPointNodal(node)(2)*grad(1);
	res(0,2) += PTReel.GetPointNodal(node)(0)*grad(2);
	res(1,2) += PTReel.GetPointNodal(node)(1)*grad(2);
	res(2,2) += PTReel.GetPointNodal(node)(2)*grad(2);
      }
  }
  
  
  //! computation of the polynomial decomposition of the jacobian |DF|
  void TetrahedronGeomReference
  ::ComputeCoefJacobian(const VectR3& s, VectReal_wp& CoefJacobian) const
  {
    CoefJacobian.Reallocate(1);
    Matrix3_3 DF; R3 pt; DFjLinear(s, pt, DF);
    CoefJacobian(0) = Det(DF);
    /*
    CoefComatrix.Reallocate(1);
    GetInverse(DF); Transpose(DF);    
    CoefComatrix(0) = DF;
    Mlt(CoefJacobian(0), CoefComatrix(0)); */
  }
  
  
  /**********************
   * FjElem and DFjElem *
   **********************/
  
  
  //! computes res = Fi(point) for all reference points
  /*!
    \param[in] s list of vertices of the element
    \param[out] res "reference points" after transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  void TetrahedronGeomReference::
  FjElem(const VectR3& s, SetPoints<Dimension3>& res,
	 const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
	if (!this->quadrature_equal_nodal)
	  FjElemQuadratureCurve(s, res, mesh, nquad);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  FjElemDofCurve(s, res, mesh, nquad);
      }
    else
      {
	FjElemNodalLinear(s, res);
	if (!this->quadrature_equal_nodal)
	  FjElemQuadratureLinear(s, res);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  FjElemDofLinear(s, res);
      }
  }
  
  
  //! computes res = Fi(point) for all nodal points
  void TetrahedronGeomReference::
  FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
	      const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
    else
      FjElemNodalLinear(s, res);
  }
  
  
  //! computes res = Fi(point) for all quadrature points
  void TetrahedronGeomReference::
  FjElemQuadrature(const VectR3& s, SetPoints<Dimension3>& res,
		   const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  FjElemQuadratureCurve(s, res, mesh, nquad);
      }
    else
      FjElemQuadratureLinear(s, res);
  }
  
  
  //! computes res = Fi(point) for all dof points
  void TetrahedronGeomReference::
  FjElemDof(const VectR3& s, SetPoints<Dimension3>& res,
	    const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else
	  FjElemDofCurve(s, res, mesh, nquad);
      }
    else
      FjElemDofLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all reference points
  /*!
    \param[in] s list of vertices of the element
    \param[in] PTReel "reference points" after transformation Fi
    \param[out] res jacobian matrices DFi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  void TetrahedronGeomReference::
  DFjElem(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	  SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh,
	  int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (!this->quadrature_equal_nodal)
	  DFjElemQuadratureCurve(s, PTReel, res, mesh, nquad);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  DFjElemDofCurve(s, PTReel, res, mesh, nquad);
      }
    else
      {
	DFjElemNodalLinear(s, res);
	if (!this->quadrature_equal_nodal)
	  DFjElemQuadratureLinear(s, res);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  DFjElemDofLinear(s, res);
      }
  }
  
  
  //! computes res = DFi(point) for all nodal points
  void TetrahedronGeomReference::
  DFjElemNodal(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	       SetMatrices<Dimension3>& res,
	       const Mesh<Dimension3>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      DFjElemNodalCurve(s,PTReel,res,mesh,nquad);
    else
      DFjElemNodalLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all quadrature points
  void TetrahedronGeomReference::
  DFjElemQuadrature(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		    SetMatrices<Dimension3>& res,
		    const Mesh<Dimension3>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  DFjElemQuadratureCurve(s, PTReel, res, mesh, nquad);
      }
    else
      DFjElemQuadratureLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all dof points
  void TetrahedronGeomReference::
  DFjElemDof(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	     SetMatrices<Dimension3>& res,
	     const Mesh<Dimension3>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else
	  DFjElemDofCurve(s, PTReel, res, mesh, nquad);
      }
    else
      DFjElemDofLinear(s, res);
  }
  
  
  //! nodal points in the case of straight tetrahedron
  void TetrahedronGeomReference::
  FjElemNodalLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
      {
	res_n.Fill(0);
	Real_wp l1 = Real_wp(1) - this->points_nodal_nd(i)(0) - this->points_nodal_nd(i)(1) 
	  - this->points_nodal_nd(i)(2);
	
	Add(l1, s(0), res_n);
	Add(this->points_nodal_nd(i)(0), s(1), res_n);
	Add(this->points_nodal_nd(i)(1), s(2), res_n);
	Add(this->points_nodal_nd(i)(2), s(3), res_n);
	
	res.SetPointNodal(i, res_n);
      }
  }
  
  
  //! quadrature points in the case of straight tetrahedron
  void TetrahedronGeomReference::
  FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int i = 0; i < this->points3d.GetM(); i++)
      {
	res_n.Fill(0);
	Real_wp l1
          = Real_wp(1) - this->points3d(i)(0) - this->points3d(i)(1) - this->points3d(i)(2);
	Add(l1, s(0), res_n);
	Add(this->points3d(i)(0), s(1), res_n);
	Add(this->points3d(i)(1), s(2), res_n);
	Add(this->points3d(i)(2), s(3), res_n);
	    
	res.SetPointQuadrature(i, res_n);
      }
    
  }
  
  
  //! dof points in the case of straight tetrahedron
  void TetrahedronGeomReference::
  FjElemDofLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int i = 0; i < this->points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	Real_wp l1 = Real_wp(1) - this->points_dof3d(i)(0) - this->points_dof3d(i)(1)
	  - this->points_dof3d(i)(2);
	Add(l1, s(0), res_n); Add(this->points_dof3d(i)(0), s(1), res_n);
	Add(this->points_dof3d(i)(1),s(2),res_n); Add(this->points_dof3d(i)(2),s(3),res_n);
	    
	res.SetPointDof(i, res_n);
      }
  }
  
  
  //! transformation DFi for nodal points in the case of straight tetrahedron
  void TetrahedronGeomReference::
  DFjElemNodalLinear(const VectR3& s, SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    DFjLinear(s, s(0), mat);
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    res.FillNodal(mat);
  }
  
  
  //! transformation DFi for quadrature points in the case of straight tetrahedron
  void TetrahedronGeomReference::
  DFjElemQuadratureLinear(const VectR3& s, SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    DFjLinear(s, s(0), mat);
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    res.FillQuadrature(mat);
  }
  
  
  //! transformation DFi for dof points in the case of straight tetrahedron
  void TetrahedronGeomReference::
  DFjElemDofLinear(const VectR3& s, SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    DFjLinear(s, s(0), mat);
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    res.FillDof(mat);
  }
  
  
  //! transformation Fi for nodal points in the case of curved tetrahedron
  /*!
    \param[in] s four vertices of the element
    \param[out] res references points after the transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number in the mesh
   */
  void TetrahedronGeomReference::
  FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
		   const Mesh<Dimension3>& mesh, int nquad, const Volume& elt) const
  {
    ElementGeomReference<Dimension3>::FjElemNodalCurve(s, res, mesh, nquad);
    
    R3 res_n;
    // for other nodal points, we use coefFi_curve
    for (int node = 2*this->order_geom*this->order_geom+2;
         node < this->points_nodal_nd.GetM(); node++)
      {
	res_n.Zero();
	for (int k = 0; k < 2*this->order_geom*this->order_geom+2; k++)
	  Add(coefFi_curve(node, k), res.GetPointNodal(k), res_n);
	
	res.SetPointNodal(node, res_n);
      }
  }
  
  
  //! transformation Fi for quadrature points in the case of curved tetrahedron
  void TetrahedronGeomReference::
  FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
			const Mesh<Dimension3>& mesh, int nquad) const
  {
    R3 res_n; int offset = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int i = 0; i < this->points3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointQuadrature(i,res_n);
      }
  }
  
  
  //! transformation Fi for dof points in the case of curved tetrahedron
  void TetrahedronGeomReference::
  FjElemDofCurve(const VectR3& s, SetPoints<Dimension3>& res,
		 const Mesh<Dimension3>& mesh, int nquad) const
  {
    R3 res_n; int offset = this->points_nodal_nd.GetM() + this->points3d.GetM();
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int i = 0; i < this->points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  Add(coefFi_curve(offset+i, node) ,res.GetPointNodal(node), res_n);
	
	res.SetPointDof(i,res_n);
      }
  }
  
  
  //! transformation DFi for nodal points in the case of curved tetrahedron
  void TetrahedronGeomReference::
  DFjElemNodalCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		    SetMatrices<Dimension3>& res,
		    const Mesh<Dimension3>& mesh, int nquad) const
  {
    Matrix3_3 tmp;
    R3 grad; grad.Fill(0);
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
      {
	tmp.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  {
	    grad = this->GradientPhi_Nodal(node, i);
	    tmp(0,0) += PTReel.GetPointNodal(node)(0)*grad(0);
	    tmp(1,0) += PTReel.GetPointNodal(node)(1)*grad(0);
            tmp(2,0) += PTReel.GetPointNodal(node)(2)*grad(0);
	    
            tmp(0,1) += PTReel.GetPointNodal(node)(0)*grad(1);
	    tmp(1,1) += PTReel.GetPointNodal(node)(1)*grad(1);
            tmp(2,1) += PTReel.GetPointNodal(node)(2)*grad(1);
            
            tmp(0,2) += PTReel.GetPointNodal(node)(0)*grad(2);
	    tmp(1,2) += PTReel.GetPointNodal(node)(1)*grad(2);
            tmp(2,2) += PTReel.GetPointNodal(node)(2)*grad(2);
	  }
        
	res.SetPointNodal(i,tmp);
      } 
  }
  
  
  //! transformation DFi for quadrature points in the case of curved tetrahedron
  void TetrahedronGeomReference::
  DFjElemQuadratureCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			 SetMatrices<Dimension3>& res,
			 const Mesh<Dimension3>& mesh, int nquad) const
  {
    Matrix3_3 res_n; int offset = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int i = 0; i < this->points3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointQuadrature(i, res_n);
      }
  }
  
  
  //! transformation DFi for dof points in the case of curved tetrahedron
  void TetrahedronGeomReference::
  DFjElemDofCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		  SetMatrices<Dimension3>& res,
		  const Mesh<Dimension3>& mesh, int nquad) const
  {
    Matrix3_3 res_n; int offset = this->points_nodal_nd.GetM() + this->points3d.GetM();
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int i = 0; i < this->points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointDof(i,res_n);
      }
  }
  
  
  /*******************
   * Other functions *
   *******************/
  
  
  void TetrahedronGeomReference::ComputeValuesNodalPhi1D(const Real_wp&, VectReal_wp&) const
  {
    cout << "not used " << endl;
    abort();
  }

  
  //! returns \f$ \varphi_{node}(pointloc) \f$
  /*!
    where \f$ \varphi_{node} \f$ is the basis function
    associated with the nodal point \f$ \xi_{node} \f$
  */  
  void TetrahedronGeomReference
  ::ComputeValuesPhiNodalRef(const R3 &pointloc, VectReal_wp& phi) const
  {
    VectReal_wp psi; phi.Reallocate(this->points_nodal_nd.GetM()); phi.Fill(0);
    ComputeValuesPhiOrthoRef(this->order_geom, NumOrtho3D, InvWeightPolynomial, pointloc, psi);
    Mlt(InverseVDM, psi, phi);
  }
  
  
  //! Evaluating gradient of shape functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of shape functions on pointloc
  */  
  void TetrahedronGeomReference::
  ComputeGradientPhiNodalRef(const R3& pointloc, VectR3& grad_phi) const
  {
    VectReal_wp psi, phi; VectR3 grad_psi; 
    psi.Reallocate(this->points_nodal_nd.GetM()); psi.Fill(0);
    phi.Reallocate(this->points_nodal_nd.GetM()); phi.Fill(0);
    grad_phi.Reallocate(this->points_nodal_nd.GetM());
    ComputeGradientPhiOrthoRef(this->order_geom, NumOrtho3D,
                               InvWeightPolynomial, pointloc, grad_psi);
    for (int i = 0; i < psi.GetM(); i++)
      psi(i) = grad_psi(i)(0);
	
    Mlt(InverseVDM, psi, phi);
    
    for (int i = 0; i < psi.GetM(); i++)
      {
	grad_phi(i)(0) = phi(i);
	psi(i) = grad_psi(i)(1);
      }
    
    Mlt(InverseVDM, psi, phi);
    for (int i = 0; i < psi.GetM(); i++)
      {
	grad_phi(i)(1) = phi(i);
	psi(i) = grad_psi(i)(2);
      }
    
    Mlt(InverseVDM, psi, phi);
    for (int i = 0; i < psi.GetM(); i++)
      grad_phi(i)(2) = phi(i);
  }
  

  //! Evaluating orthogonal functions at pointloc  
  void TetrahedronGeomReference::
  ComputeValuesPhiOrthoRef(int r, const Array3D<int>& NumOrtho, const VectReal_wp& InvWeight,
                           const R3& ptloc, VectReal_wp& phi) const
  {
    // orthogonal polynomials :
    // L_i(2x/(1-y-z) - 1) P_j^{2i+1}(2y/(1-z)-1) (1-y-z)^i P_k^{2i+2j+2}(2z-1) (1-z)^j
    // = L_i(a) P_j^(2i+1)(b) ((1-b)/2)^i P_k^(2i+2j+2)(c) ((1-c)/2)^(i+j)
    // where (a, b, c) is on the cube [-1,1]^3 and (x,y,z) on the unit tetrahedron
    
    phi.Reallocate(InvWeight.GetM());
    Real_wp x = ptloc(0), y = ptloc(1), z = ptloc(2);
    /* for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	for (int k = 0; k <= r-i-j; k++)
	  phi(NumOrtho(i, j, k)) = pow(x, i)*pow(y, j)*pow(z, k);
	  
	  return;*/
    
    // coordinates on the cube [-1, 1]
    Real_wp a(-1), b(-1), c(-1);
    if (abs(1.0-y-z) > epsilon_machine)
      a = 2.0*x/(1.0-y-z) - 1.0;
    
    if (abs(1.0-z) > epsilon_machine)
      b = 2.0*y/(1.0-z) - 1.0;
    
    c = 2.0*z - 1.0;
    
    VectReal_wp Px;
    EvaluateJacobiPolynomial(LegendrePolynom, r, a, Px);
    Vector<VectReal_wp> Py(r+1);
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(OddJacobiPolynom(i), r-i, b, Py(i));
    
    Vector<VectReal_wp> Pz(r+1);
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, c, Pz(i));
    
    Real_wp valx, valy;
    Real_wp pow_oneMinusYZ = 1.0, pow_oneMinusZ = 1.0;
    for (int i = 0; i <= r; i++)
      {
        // computation of valx = L_i(a) (1-y-z)^i
        if (i == 0)
          valx = 1.0;
        else if (i == 1)
          valx = 2.0*x + y + z - 1.0;
        else
          valx = pow_oneMinusYZ*Px(i);
        
        pow_oneMinusZ = 1.0;
        for (int j = 0; j <= r-i; j++)
          {
            // computation of valy = P_j^(2i+1)(a) (1-z)^j
            if (j == 0)
              valy = 1.0;
            else if (j == 1)
              valy = 2.0*y + z - 1.0 - OddJacobiPolynom(i)(0,0)*(1.0-z);
            else
              valy = pow_oneMinusZ*Py(i)(j);
            
            for (int k = 0; k <= r-i-j; k++)
              {
                int node = NumOrtho(i, j, k);
                phi(node) = valx*valy*Pz(i+j)(k);
              }

            pow_oneMinusZ *= 1.0-z;
          }
        
        pow_oneMinusYZ *= 1.0-y-z;
      }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < phi.GetM(); i++)
      phi(i) *= InvWeight(i);
  }
  

  //! evaluating gradient of orthogonal functions at pointloc  
  void TetrahedronGeomReference::
  ComputeGradientPhiOrthoRef(int r, const Array3D<int>& NumOrtho, const VectReal_wp& InvWeight,
                             const R3& ptloc, VectR3& grad_phi) const
  {
    // orthogonal polynomials :
    // L_i(2x/(1-y-z) - 1) P_j^{2i+1}(2y/(1-z)-1) (1-y)^i P_k^{2i+2j+2}(2z-1) (1-z)^(i+j)
    // = L_i(a) P_j^(2i+1)(b) ((1-b)/2)^i P_k^(2i+2j+2)(c) ((1-c)/2)^(i+j)
    // where (a, b, c) is on the cube [-1,1]^3 and (x,y,z) on the unit tetrahedron
    
    // derivative with respect to a, b, c is equal to :
    // L'_i(a) P_j^(2i+1)(b) ((1-b)/2)^i P_k^(2i+2j+2)(c) ((1-c)/2)^(i+j)
    // L_i(a) ((1-b)/2)^(i-1) ( (1-b)/2 P'_j^(2i+1)(b) - i/2 P_j^(2i+1)(b) ) 
    //         P_k^(2i+2j+2)(c) ((1-c)/2)^(i+j)
    // L_i(a) P_j^(2i+1)(b) ((1-b)/2)^i ((1-c)/2)^(i+j-1)
    //        ( (1-c)/2 P'_k^(2i+2j+2)(c) - (i+j)/2 P_k^(2i+2j+2)(c)) 
    grad_phi.Reallocate(InvWeight.GetM());
    Real_wp x = ptloc(0), y = ptloc(1), z = ptloc(2);
    // coordinates on the cube [-1, 1]
    Real_wp a(-1), b(-1), c(-1);
    Real_wp da_dx(0), da_dy(0), db_dy(0), db_dz(0);
    if (abs(1.0-y-z) > epsilon_machine)
      {
	a = 2.0*x/(1.0-y-z) - 1.0;
	da_dx = 2.0/(1.0-y-z);
	da_dy = 2.0*x/square(1.0-y-z);
      }
    
    if (abs(1.0-z) > epsilon_machine)
      {
	b = 2.0*y/(1.0-z) - 1.0;
	db_dy = 2.0/(1.0-z);
	db_dz = 2.0*y/square(1.0-z);
      }
    
    c = 2.0*z - 1.0;
    
    VectReal_wp Px, dPx;
    EvaluateJacobiPolynomial(LegendrePolynom, r, a, Px, dPx);
    Vector<VectReal_wp> Py(r+1), dPy(r+1);
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(OddJacobiPolynom(i), r-i, b, Py(i), dPy(i));
    
    Vector<VectReal_wp> Pz(r+1), dPz(r+1);
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, c, Pz(i), dPz(i));
    
    Real_wp valx, valy, valz;
    Real_wp dvalx_dx, dvalx_dy, dvaly_dy, dvaly_dz, dvalz_dz;
    Real_wp pow_oneMinusYZ = 1.0, pow_oneMinusZ = 1.0;    
    Real_wp powm1_oneMinusYZ = 0.0, powm1_oneMinusZ = 0.0;    
    for (int i = 0; i <= r; i++)
      {
        if (i == 0)
          {
            valx = 1.0;
            dvalx_dx = 0;
            dvalx_dy = 0;
          }
        else if (i == 1)
          {
            valx = 2.0*x + y + z - 1.0;
            dvalx_dx = 2.0;
            dvalx_dy = 1.0;
          }
        else
          {
            valx = pow_oneMinusYZ*Px(i);
            dvalx_dx = pow_oneMinusYZ*dPx(i)*da_dx;
            dvalx_dy = pow_oneMinusYZ*dPx(i)*da_dy - i*powm1_oneMinusYZ*Px(i);
          }
        
        powm1_oneMinusZ = 0; pow_oneMinusZ = 1.0;
        for (int j = 0; j <= r-i; j++)
          {
            if (j == 0)
              {
                valy = 1.0;
                dvaly_dy = 0;
                dvaly_dz = 0;
              }
            else if (j == 1)
              {
                valy = 2.0*y + z - 1.0 - OddJacobiPolynom(i)(0,0)*(1.0-z);
                dvaly_dy = 2.0;
                dvaly_dz = 1.0 + OddJacobiPolynom(i)(0,0);
              }
            else
              {
                valy = pow_oneMinusZ*Py(i)(j);
                dvaly_dy = pow_oneMinusZ*dPy(i)(j)*db_dy;
                dvaly_dz = pow_oneMinusZ*dPy(i)(j)*db_dz - j*powm1_oneMinusZ*Py(i)(j);
              }
                   
            for (int k = 0; k <= r-i-j; k++)
              {
                int node = NumOrtho(i, j, k);
                valz = Pz(i+j)(k);
                dvalz_dz = 2.0*dPz(i+j)(k);
                grad_phi(node).Init(dvalx_dx*valy*valz, (dvalx_dy*valy + valx*dvaly_dy)*valz, 
                                    dvalx_dy*valy*valz + valx*dvaly_dz*valz + valx*valy*dvalz_dz);
              }
            
            powm1_oneMinusZ = pow_oneMinusZ;
            pow_oneMinusZ *= 1.0-z;
          }
        
        powm1_oneMinusYZ = pow_oneMinusYZ;
	pow_oneMinusYZ *= 1.0-y-z;
      }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < grad_phi.GetM(); i++)
      grad_phi(i) *= InvWeight(i);
    
  }
    
    
  //! displays informations class TetrahedronGeomReference
  ostream& operator <<(ostream& out, const TetrahedronGeomReference& e)
  {
    out << static_cast<const ElementGeomReference<Dimension3>& >(e);
    return out;
  }
  
} // end namespace

#define MONTJOIE_FILE_TETRAHEDRON_GEOM_REFERENCE_CXX
#endif
