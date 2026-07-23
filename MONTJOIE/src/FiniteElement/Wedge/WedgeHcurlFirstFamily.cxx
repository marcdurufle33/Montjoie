#ifndef MONTJOIE_FILE_WEDGE_HCURL_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  WedgeHcurlFirstFamily::WedgeHcurlFirstFamily() : WedgeReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;

    type_integration_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
  }
  
    
  //! how to number mesh
  void WedgeHcurlFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return WedgeReference<2>::ConstructNumberMap(nmap, dg);
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, this->order*(this->order+2));
	nmap.SetNbDofQuadrangle(this->order, 2*this->order*(this->order+1));
	nmap.SetNbDofWedge(this->order, 0);

	// rotation of dofs on faces
	element_quad_surf->FindHcurlRotationQuad(nmap, 0);
	element_tri_surf->FindHcurlRotationTri(nmap, 0);
	
	return;
      }
    
    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofEdge(order, order);
    // nothing on quadrangular face !
    nmap.SetNbDofQuadrangle(order, 2*order*(order-1));
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order-1));
    
    // dofs inside the tetrahedron
    nmap.SetNbDofWedge(order, nb_dof_loc-nb_dof_boundaries);

    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
    
    // FindHcurlSignEdge(nmap);
    // rotation of dofs on faces
    FindHcurlLinearCombinationRotation(nmap, *element_tri_surf, *element_quad_surf);
  }


  size_t WedgeHcurlFirstFamily::GetMemorySize() const
  {
    size_t taille = WedgeReference<2>::GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    taille += lobz_hcurl.GetMemorySize() + lobz_h1.GetMemorySize();
    taille += GL_GX.GetMemorySize();
    taille += element_tri_h1.GetMemorySize();
    return taille;
  }

  
  //! construction of finite element
  void WedgeHcurlFirstFamily::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						     int rsurf_tri, int rsurf_quad,
						     int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    // type_quad = TriangleQuadrature::QUADRATURE_TENSOR;
    WedgeReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad,
                                              r, r, TriangleQuadrature::QUADRATURE_GAUSS,
                                              Globatto<Real_wp>::QUADRATURE_LOBATTO,
					      Globatto<Real_wp>::QUADRATURE_LOBATTO);
    
    // computation of basis functions
    ConstructFunctions();
    this->nb_points_dof_inside = this->nb_dof_loc;
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    // construction of Value_PhiVec, Curl_Phi and elementary matrices
    ConstructHcurlElementaryMatrix();
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace(true);
    
    TriangleHcurlFirstFamily* Fb_tri = new TriangleHcurlFirstFamily();
    QuadrangleHcurlFirstFamily* Fb_quad = new QuadrangleHcurlFirstFamily();
    Fb_tri->ConstructFiniteElement(order);
    Fb_tri->type_nodal_basis = Fb_tri->NODAL_GAUSS;        
    Fb_quad->ConstructFiniteElement(order);
    element_tri_surf = Fb_tri;
    element_quad_surf = Fb_quad;

    VectR2 points_dof2d_tri = Fb_tri->PointsDofND();
    this->SetPointsDof2D_tri(points_dof2d_tri);
    
    VectR2 points_dof2d_quad = Fb_quad->PointsDofND();
    this->SetPointsDof2D_quad(points_dof2d_quad);
    
    this->num_dof_points_surf.Reallocate(5);
    for (int n = 0; n < 5; n++)
      {
	this->num_dof_points_surf(n).Reallocate(this->GetNbDofBoundary(n));
	for (int i = 0; i < this->GetNbDofBoundary(n); i++)
	  this->num_dof_points_surf(n)(i) = FacesDof(i, n);
      }    

  }
  
  
  //! construction of basis functions
  void WedgeHcurlFirstFamily::ConstructFunctions()
  {
    int r = order;
    nb_dof_loc = 3*order*(order+2)*(order+1)/2;
    
    nb_dof_boundaries = 9*order + 2*order*(order-1) + 6*order*(order-1);
    nb_dof_quad = 4*order + 2*order*(order-1);
    nb_dof_tri = 3*order + order*(order-1);
    
    // triangular finite elements
    TriangleHcurlFirstFamily* element_tri_hcurl = new TriangleHcurlFirstFamily();
    element_tri_surf = element_tri_hcurl;

    element_tri_h1.ConstructFiniteElement(order);
    element_tri_hcurl->ConstructFiniteElement(order);
    
    // Gauss and Gauss-Lobatto points
    lobz_h1.ConstructQuadrature(order-1, Globatto<Real_wp>::QUADRATURE_GAUSS);
    lobz_hcurl.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_LOBATTO);

    VectReal_wp points_gauss, poids_gauss;
    ComputeGaussLegendre(points_gauss, poids_gauss, order);
    
    GL_GX.Reallocate(r+1, r+1);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	GL_GX(i, j) = lobz_hcurl.EvaluatePhi(i, points_gauss(j));
    
    VectReal_wp points_dof1d; VectR3 points_dof3d;
    points_dof1d = lobz_h1.Points();
    
    points_dof3d.Reallocate(nb_dof_loc);
    tangente_dof.Reallocate(nb_dof_loc);
    
    CoordinateDofs.Reallocate(nb_dof_loc, 3);
    CoordinateDofs.Fill(-2);
    // dofs on edges
    for (int i = 0; i < order; i++)
      {
        CoordinateDofs(i, 0) = i;
        CoordinateDofs(i, 1) = 0;
        CoordinateDofs(i, 2) = 1;
	tangente_dof(i).Init(1, 0, 0);
	points_dof3d(i).Init(points_dof1d(i), 0, 0);
        
        CoordinateDofs(r+i, 0) = r+i;
        CoordinateDofs(r+i, 1) = 0;
        CoordinateDofs(r+i, 2) = 1;
	tangente_dof(r+i).Init(-1, 1, 0);
	points_dof3d(r+i).Init(1.0-points_dof1d(i), points_dof1d(i), 0);

        CoordinateDofs(2*r+i, 0) = 2*r + r - 1 - i;
        CoordinateDofs(2*r+i, 1) = 0;
        CoordinateDofs(2*r+i, 2) = -1;
	tangente_dof(2*r+i).Init(0, 1, 0);
	points_dof3d(2*r+i).Init(0, points_dof1d(i), 0);

        CoordinateDofs(3*r+i, 0) = 0;
        CoordinateDofs(3*r+i, 1) = i;
        CoordinateDofs(3*r+i, 2) = 2;
	tangente_dof(3*r+i).Init(0, 0, 1);
	points_dof3d(3*r+i).Init(0, 0, points_dof1d(i));

        CoordinateDofs(4*r+i, 0) = 1;
        CoordinateDofs(4*r+i, 1) = i;
        CoordinateDofs(4*r+i, 2) = 2;
	tangente_dof(4*r+i).Init(0, 0, 1);
	points_dof3d(4*r+i).Init(1, 0, points_dof1d(i));

        CoordinateDofs(5*r+i, 0) = 2;
        CoordinateDofs(5*r+i, 1) = i;
        CoordinateDofs(5*r+i, 2) = 2;
	tangente_dof(5*r+i).Init(0, 0, 1);
	points_dof3d(5*r+i).Init(0, 1, points_dof1d(i));

        CoordinateDofs(6*r+i, 0) = i;
        CoordinateDofs(6*r+i, 1) = r;
        CoordinateDofs(6*r+i, 2) = 1;
	tangente_dof(6*r+i).Init(1, 0, 0);
	points_dof3d(6*r+i).Init(points_dof1d(i), 0, 1);
        
        CoordinateDofs(7*r+i, 0) = r+i;
        CoordinateDofs(7*r+i, 1) = r;
        CoordinateDofs(7*r+i, 2) = 1;
	tangente_dof(7*r+i).Init(-1, 1, 0);
	points_dof3d(7*r+i).Init(1.0-points_dof1d(i), points_dof1d(i), 1);

        CoordinateDofs(8*r+i, 0) = 2*r + r - 1 - i;
        CoordinateDofs(8*r+i, 1) = r;
        CoordinateDofs(8*r+i, 2) = -1;
	tangente_dof(8*r+i).Init(0, 1, 0);
	points_dof3d(8*r+i).Init(0, points_dof1d(i), 1);
      }
    
    // dofs on faces
    int offset = 9*r;

    const VectR2& tangente_tri = element_tri_surf->GetOrientationDofs();    
    // face z = 0
    for (int i = 3*r; i < element_tri_surf->GetNbDof(); i++)
      {
        CoordinateDofs(offset, 0) = i;
        CoordinateDofs(offset, 1) = 0;
        CoordinateDofs(offset, 2) = 1;
	Real_wp tx = tangente_tri(i)(0);
	Real_wp ty = tangente_tri(i)(1);
	tangente_dof(offset).Init(tx, ty, 0);
	Real_wp x = element_tri_surf->PointsDofND(i)(0);
	Real_wp y = element_tri_surf->PointsDofND(i)(1);
	points_dof3d(offset).Init(x, y, 0);
        offset++;
      }
    
    // face y = 0
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r-1; j++)
        {
          CoordinateDofs(offset, 0) = i;
          CoordinateDofs(offset, 1) = j+1;
          CoordinateDofs(offset, 2) = 1;
	  tangente_dof(offset).Init(1, 0, 0);
	  points_dof3d(offset).Init(points_dof1d(i), 0, lobz_hcurl.Points(j+1));
          offset++;
        }
    
    for (int i = 0; i < r-1; i++)
      for (int j = 0; j < r; j++)
        {
          CoordinateDofs(offset, 0) = 3 + i;
          CoordinateDofs(offset, 1) = j;
          CoordinateDofs(offset, 2) = 2;
	  tangente_dof(offset).Init(0, 0, 1);
	  points_dof3d(offset).Init(lobz_hcurl.Points(i+1), 0, points_dof1d(j));
          offset++;
        }

    // face x+y = 1
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r-1; j++)
        {
          CoordinateDofs(offset, 0) = r+i;
          CoordinateDofs(offset, 1) = j+1;
          CoordinateDofs(offset, 2) = 1;
	  tangente_dof(offset).Init(-1, 1, 0);
	  points_dof3d(offset).Init(1.0-points_dof1d(i), points_dof1d(i), lobz_hcurl.Points(j+1));
          offset++;
        }
    
    for (int i = 0; i < r-1; i++)
      for (int j = 0; j < r; j++)
        {
          CoordinateDofs(offset, 0) = 3 + (r-1) + i;
          CoordinateDofs(offset, 1) = j;
          CoordinateDofs(offset, 2) = 2;
	  tangente_dof(offset).Init(0, 0, 1);
	  points_dof3d(offset)
            .Init(1.0-lobz_hcurl.Points(i+1), lobz_hcurl.Points(i+1), points_dof1d(j));
          offset++;
        }

    // face x = 0
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r-1; j++)
        {
          CoordinateDofs(offset, 0) = 2*r + r-1-i;
          CoordinateDofs(offset, 1) = j+1;
          CoordinateDofs(offset, 2) = -1;
	  tangente_dof(offset).Init(0, 1, 0);
	  points_dof3d(offset).Init(0, points_dof1d(i), lobz_hcurl.Points(j+1));
          offset++;
        }
    
    for (int i = 0; i < r-1; i++)
      for (int j = 0; j < r; j++)
        {
          CoordinateDofs(offset, 0) = 3 + 2*(r-1) + r-2-i;
          CoordinateDofs(offset, 1) = j;
          CoordinateDofs(offset, 2) = 2;
	  tangente_dof(offset).Init(0, 0, 1);
	  points_dof3d(offset).Init(0, lobz_hcurl.Points(i+1), points_dof1d(j));
          offset++;
        }
    
    // face z = 1
    for (int i = 3*r; i < element_tri_surf->GetNbDof(); i++)
      {
        CoordinateDofs(offset, 0) = i;
        CoordinateDofs(offset, 1) = order;
        CoordinateDofs(offset, 2) = 1;
	Real_wp tx = tangente_tri(i)(0);
	Real_wp ty = tangente_tri(i)(1);
	tangente_dof(offset).Init(tx, ty, 0);
	Real_wp x = element_tri_surf->PointsDofND(i)(0);
	Real_wp y = element_tri_surf->PointsDofND(i)(1);
	points_dof3d(offset).Init(x, y, 1.0);
        offset++;
      }
    
    // interior
    for (int i = 3*r; i < element_tri_surf->GetNbDof(); i++) 
      for (int j = 0; j <= r-2; j++)
        {
          CoordinateDofs(offset, 0) = i;
          CoordinateDofs(offset, 1) = j+1;
          CoordinateDofs(offset, 2) = 1;
	  Real_wp tx = tangente_tri(i)(0);
	  Real_wp ty = tangente_tri(i)(1);
	  tangente_dof(offset).Init(tx, ty, 0);
	  Real_wp x = element_tri_surf->PointsDofND(i)(0);
	  Real_wp y = element_tri_surf->PointsDofND(i)(1);
	  points_dof3d(offset).Init(x, y, lobz_hcurl.Points(j+1));
          offset++;
        }
    
    for (int i = 3*r; i < element_tri_h1.GetNbPointsNodalElt(); i++) 
      for (int j = 0; j <= r-1; j++)
        {
          CoordinateDofs(offset, 0) = i;
          CoordinateDofs(offset, 1) = j;
          CoordinateDofs(offset, 2) = 2;
	  tangente_dof(offset).Init(0, 0, 1);
	  Real_wp x = element_tri_h1.PointsNodalND(i)(0);
	  Real_wp y = element_tri_h1.PointsNodalND(i)(1);
	  points_dof3d(offset).Init(x, y, points_dof1d(j));
          offset++;
        }
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof3d);

    /*
    for (int i = 0; i < points_dof3d.GetM(); i++)
      {
	VectR3 phi;
	ComputeValuesPhiRef(points_dof3d(i), phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  {
	    Real_wp phin = DotProd(phi(j), tangente_dof(i));
	    Real_wp val_ex = 0;
	    if (i == j)
	      val_ex = 1.0;
	    
	    if (abs(phin - val_ex) > 1e5*epsilon_machine)
	      {
		DISP(i); DISP(j);
		DISP(phin); DISP(val_ex);
		abort();
	      }
	  }
	  } */
  }
  

  //! Integration against basis functions on a face
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc face number
    (\omega, \xi_k) are here Gauss points and not Gauss-Lobatto points
  */  
  template<class Vector1, class Vector2>
  void WedgeHcurlFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    if (num_loc%4 != 0)
      {
	const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
	int Nquad = (order+1)*(order+1);
	Vector1 fx(3*Nquad), fy(3*Nquad);
	FillZero(fx); FillZero(fy);
	// we express feval on Gauss-Lobatto points
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int node = 3*NumQuad2D(i, j);
	      for (int m = 0; m <= order; m++)
		{
		  int p = 3*NumQuad2D(m, j);
		  fx(node) += GL_GX(i, m)*feval(p);
		  fx(node+1) += GL_GX(i, m)*feval(p+1);
		  fx(node+2) += GL_GX(i, m)*feval(p+2);
		}
	    }
	
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int node = 3*NumQuad2D(i, j);
	      for (int m = 0; m <= order; m++)
		{
		  int p = 3*NumQuad2D(i, m);
		  fy(node) += GL_GX(j, m)*fx(p);
		  fy(node+1) += GL_GX(j, m)*fx(p+1);
		  fy(node+2) += GL_GX(j, m)*fx(p+2);
		}
	    }
	
	// then we use integration with Gauss-Lobatto points
	ComputeIntegralSurfaceRef(fy, res, num_loc);
      }
    else
      ComputeIntegralSurfaceRef(feval, res, num_loc);
  }

    
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void WedgeHcurlFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    R2 pt2D(point_loc(0), point_loc(1));
    Real_wp z = point_loc(2);
    VectReal_wp phiG(order), phiGL(order+1);
    for (int i = 0; i < order; i++)
      phiG(i) = lobz_h1.EvaluatePhi(i, z);
    
    for (int i = 0; i <= order; i++)
      phiGL(i) = lobz_hcurl.EvaluatePhi(i, z);
    
    VectReal_wp phi_scal; VectR2 phi_vec;
    element_tri_h1.ComputeValuesPhiNodalRef(pt2D, phi_scal);
    element_tri_surf->ComputeValuesPhiRef(pt2D, phi_vec);
    
    for (int i = 0; i < nb_dof_loc; i++)
      {
        int node_tri = CoordinateDofs(i, 0);
        int k = CoordinateDofs(i, 1);
        int direction = CoordinateDofs(i, 2);
        if (direction == 2)
          {
            phi(i).Init(0, 0, phi_scal(node_tri)*phiG(k));
          }
        else
          {
            if (direction == 1)
              phi(i).Init(phi_vec(node_tri)(0)*phiGL(k), phi_vec(node_tri)(1)*phiGL(k), 0);
            else
              phi(i).Init(-phi_vec(node_tri)(0)*phiGL(k), -phi_vec(node_tri)(1)*phiGL(k), 0);
          }
      }
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] curl curl of basis functions on point_loc
  */
  void WedgeHcurlFirstFamily::ComputeCurlPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    R2 pt2D(point_loc(0), point_loc(1));
    Real_wp z = point_loc(2);
    VectReal_wp phiG(order), phiGL(order+1), dphiGL(order+1);
    for (int i = 0; i < order; i++)
      phiG(i) = lobz_h1.EvaluatePhi(i, z);
    
    for (int i = 0; i <= order; i++)
      {
        phiGL(i) = lobz_hcurl.EvaluatePhi(i, z);
        dphiGL(i) = lobz_hcurl.EvaluatePhiGrad(i, z);
      }
    
    VectReal_wp phi_scal; VectR2 phi_vec;
    VectR2 grad_scal; VectReal_wp curl_vec;
    element_tri_h1.ComputeValuesPhiNodalRef(pt2D, phi_scal);
    element_tri_surf->ComputeValuesPhiRef(pt2D, phi_vec);

    element_tri_h1.ComputeGradientPhiNodalRef(pt2D, grad_scal);
    element_tri_surf->ComputeCurlPhiRef(pt2D, curl_vec);
    
    for (int i = 0; i < nb_dof_loc; i++)
      {
        int node_tri = CoordinateDofs(i, 0);
        int k = CoordinateDofs(i, 1);
        int direction = CoordinateDofs(i, 2);
        if (direction == 2)
          {
            res(i).Init(phiG(k)*grad_scal(node_tri)(1), -phiG(k)*grad_scal(node_tri)(0), 0);
          }
        else
          {
            res(i).Init(-dphiGL(k)*phi_vec(node_tri)(1), dphiGL(k)*phi_vec(node_tri)(0), 
                        phiGL(k)*curl_vec(node_tri));
            
            if (direction != 1)
              Mlt(-1.0, res(i));
          }
      }
  }
  

  //! displays details of class WedgeHcurlHpFirstFamily
  ostream& operator <<(ostream& out, const WedgeHcurlFirstFamily& e)
  {
    out<<static_cast<const WedgeReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_WEDGE_HCURL_FIRST_FAMILY_CXX
#endif
