#ifndef MONTJOIE_FILE_WEDGE_HCURL_OPTIMAL_HP_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  WedgeHcurlOptimalHpFirstFamily::WedgeHcurlOptimalHpFirstFamily() : WedgeReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  
  
  //! how to number mesh
  void WedgeHcurlOptimalHpFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return WedgeReference<2>::ConstructNumberMap(nmap, dg);

    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofEdge(order, order);
    // nothing on quadrangular face !
    nmap.SetNbDofQuadrangle(order, 2*order*order);
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order-1));
    
    // dofs inside the tetrahedron
    nmap.SetNbDofWedge(order, nb_dof_loc-nb_dof_boundaries);

    nmap.SetEqualEdgesDofSymmetry(order, order+1);
   
    //FindHcurlSignEdge(*this, nmap);
    // rotation of dofs on faces
    FindHcurlLinearCombinationRotation(nmap, *element_tri_surf, *element_quad_surf);
    // FindHcurlLinearCombinationRotation(*this, nmap);
    
  }


  size_t WedgeHcurlOptimalHpFirstFamily::GetMemorySize() const
  {
    size_t taille = WedgeReference<2>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize() + JacobiPolynom.GetMemorySize();    
    taille += NumDofsX.GetMemorySize() + NumDofsY.GetMemorySize() + NumDofsZ.GetMemorySize();
    taille += NumQuad2D_gauss.GetMemorySize() + G_GX.GetMemorySize();
    return taille;
  }

  
  //! construction of finite element
  void WedgeHcurlOptimalHpFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf_tri, int rsurf_quad,
			   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (rquad <= r)
      rquad = r+1;
    
    // type_quad = TriangleQuadrature::QUADRATURE_TENSOR;
    WedgeReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad,
					      r, r+1, -1, -1, Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    // computation of basis functions
    ConstructFunctions();
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    // construction of Value_PhiVec, Curl_Phi and elementary matrices
    ConstructHcurlElementaryMatrix();
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace(false);

    //Matrix<Real_wp, Symmetric, RowSymPacked>
    // Mh(nb_dof_loc, nb_dof_loc), Kh(nb_dof_loc, nb_dof_loc);
    //Mh.Fill(0); Kh.Fill(0); Matrix3_3 coef; coef.SetIdentity();
    //AddConstantMassMatrix(0, 0, coef, Mh);
    //AddConstantStiffnessMatrix(0, 0, coef, Kh);
    //Mh.Write("MhWedCurl"+to_str(this->order)+".dat");
    //Kh.Write("KhWedCurl"+to_str(this->order)+".dat");

  }
  
  
  //! construction of basis functions
  void WedgeHcurlOptimalHpFirstFamily::ConstructFunctions()
  {
    nb_dof_loc = order*(order+2)*(3*order+7)/2;
    
    nb_dof_boundaries = 9*order + 2*order*(order-1) + 6*order*order;
    nb_dof_quad = 4*order + 2*order*order;
    nb_dof_tri = 3*order + order*(order-1);
    
    //points_dof3d.Reallocate(nb_points_quadrature_inside);
    //for (int i = 0; i < nb_points_quadrature_inside; i++)
    //points_dof3d(i) = points3d(i);

    TriangleGeomReference element_tri_h1;
    TriangleHcurlOptimalFirstFamily element_tri_hcurl;
    
    element_tri_h1.ConstructFiniteElement(order);
    element_tri_hcurl.ConstructFiniteElement(order);
    
    VectReal_wp points_lob, poids_lob;
    ComputeGaussLobatto(points_lob, poids_lob, order+1);
    
    VectReal_wp points_dof1d; VectR3 points_dof3d;
    
    points_dof1d.Reallocate(order);
    for (int i = 0; i < order; i++)
      points_dof1d(i) = points_lob(i+1);
    
    /* points_dof3d.Reallocate(nb_dof_loc);
    tangente_dof.Reallocate(nb_dof_loc);
    
    // dofs on edges
    int r = order;
    for (int i = 0; i < order; i++)
      {
	tangente_dof(i).Init(1, 0, 0);
	points_dof3d(i).Init(points_dof1d(i), 0, 0);
        
	tangente_dof(r+i).Init(-1, 1, 0);
	points_dof3d(r+i).Init(1.0-points_dof1d(i), points_dof1d(i), 0);

	tangente_dof(2*r+i).Init(0, 1, 0);
	points_dof3d(2*r+i).Init(0, points_dof1d(i), 0);

	tangente_dof(3*r+i).Init(0, 0, 1);
	points_dof3d(3*r+i).Init(0, 0, points_dof1d(i));

	tangente_dof(4*r+i).Init(0, 0, 1);
	points_dof3d(4*r+i).Init(1, 0, points_dof1d(i));

	tangente_dof(5*r+i).Init(0, 0, 1);
	points_dof3d(5*r+i).Init(0, 1, points_dof1d(i));

	tangente_dof(6*r+i).Init(1, 0, 0);
	points_dof3d(6*r+i).Init(points_dof1d(i), 0, 1);
        
	tangente_dof(7*r+i).Init(-1, 1, 0);
	points_dof3d(7*r+i).Init(1.0-points_dof1d(i), points_dof1d(i), 1);

	tangente_dof(8*r+i).Init(0, 1, 0);
	points_dof3d(8*r+i).Init(0, points_dof1d(i), 1);
      }
    
    // dofs on faces
    int offset = 9*r;

    const Vector<R2>& tangente_tri = element_tri_hcurl.GetOrientationDofs();
    // face z = 0
    for (int i = 3*r; i < element_tri_hcurl.GetNbDof(); i++)
      {
	Real_wp tx = tangente_tri(i)(0);
	Real_wp ty = tangente_tri(i)(1);
	tangente_dof(offset).Init(tx, ty, 0);
	Real_wp x = element_tri_hcurl.PointsDofND(i)(0);
	Real_wp y = element_tri_hcurl.PointsDofND(i)(1);
	points_dof3d(offset).Init(x, y, 0);
        offset++;
      }
    
    // face y = 0
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
        {
	  tangente_dof(offset).Init(1, 0, 0);
	  points_dof3d(offset).Init(points_dof1d(i), 0, points_lob(j+1));          
          offset++;
        }
    
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
        {
	  tangente_dof(offset).Init(0, 0, 1);
	  points_dof3d(offset).Init(points_lob(i+1), 0, points_dof1d(j));          
          offset++;
        }

    // face x+y = 1
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
        {
	  tangente_dof(offset).Init(-1, 1, 0);
	  points_dof3d(offset).Init(1.0-points_dof1d(i), points_dof1d(i), points_lob(j+1));    
          offset++;
        }
    
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
        {
	  tangente_dof(offset).Init(0, 0, 1);
	  points_dof3d(offset).Init(1.0-points_lob(i+1), points_lob(i+1), points_dof1d(j));
          offset++;
        }

    // face x = 0
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
        {
	  tangente_dof(offset).Init(0, 1, 0);
	  points_dof3d(offset).Init(0, points_dof1d(i), points_lob(j+1));
          offset++;
        }
    
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
        {
	  tangente_dof(offset).Init(0, 0, 1);
	  points_dof3d(offset).Init(0, points_lob(i+1), points_dof1d(j));
          offset++;
        }
    
    // face z = 1
    for (int i = 3*r; i < element_tri_hcurl.GetNbDof(); i++)
      {
	Real_wp tx = tangente_tri(i)(0);
	Real_wp ty = tangente_tri(i)(1);
	tangente_dof(offset).Init(tx, ty, 0);
	Real_wp x = element_tri_hcurl.PointsDofND(i)(0);
	Real_wp y = element_tri_hcurl.PointsDofND(i)(1);
	points_dof3d(offset).Init(x, y, 1.0);
        offset++;
      }
    
    // interior
    for (int i = 3*r; i < element_tri_hcurl.GetNbDof(); i++) 
      for (int j = 0; j <= r-1; j++)
        {
	  Real_wp tx = tangente_tri(i)(0);
	  Real_wp ty = tangente_tri(i)(1);
	  tangente_dof(offset).Init(tx, ty, 0);
	  Real_wp x = element_tri_hcurl.PointsDofND(i)(0);
	  Real_wp y = element_tri_hcurl.PointsDofND(i)(1);
	  points_dof3d(offset).Init(x, y, points_lob(j+1));
          offset++;
        }
    
    for (int i = 3*(r+1); i < element_tri_h1.GetNbPointsNodalElt(); i++) 
      for (int j = 0; j <= r-1; j++)
        {
	  tangente_dof(offset).Init(0, 0, 1);
	  Real_wp x = element_tri_h1.PointsNodalND(i)(0);
	  Real_wp y = element_tri_h1.PointsNodalND(i)(1);
	  points_dof3d(offset).Init(x, y, points_dof1d(j));
          offset++;
        }
    */
    
    points_dof3d = this->PointsND();
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDof2D_tri(this->Points2D_tri());
    this->SetPointsDof2D_quad(this->Points2D_quad());
    this->SetPointsDofND(points_dof3d);

    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;

    GetJacobiPolynomial(LegendrePolynom, order, Real_wp(0), Real_wp(0));
    GetJacobiPolynomial(JacobiPolynom, order, Real_wp(1), Real_wp(1));

    TriangleHcurlOptimalHpFirstFamily* Fb_tri = new TriangleHcurlOptimalHpFirstFamily();
    QuadrangleHcurlOptimalHpFirstFamily* Fb_quad = new QuadrangleHcurlOptimalHpFirstFamily();
    Fb_tri->ConstructFiniteElement(order);
    Fb_quad->ConstructFiniteElement(order);
    
    element_tri_surf = Fb_tri;
    element_quad_surf = Fb_quad;

    NumDofsX.Reallocate(order, order, order);
    NumDofsY.Reallocate(order, order, order);
    NumDofsZ.Reallocate(order, order, order);
    NumDofsX.Fill(-1); NumDofsY.Fill(-1); NumDofsZ.Fill(-1);
    int offset = nb_dof_boundaries;
    for (int diag = 0; diag <= order-1; diag++)
      for (int i = 0; i < diag; i++)
	for (int j = 0; j < diag-i; j++)
	  for (int k = 0; k <= diag; k++)
	    if ((i+j) == diag-1 || k == diag)
	      {
		NumDofsX(i, j, k) = offset++;
		NumDofsY(i, j, k) = offset++;
		NumDofsZ(i, j, k) = offset++;
	      }
    
    Matrix<int> coor;
    MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(order, NumQuad2D_gauss, coor);
    Globatto<Real_wp> lob_xy;
    lob_xy.ConstructQuadrature(order+1, lob_xy.QUADRATURE_GAUSS);
    VectReal_wp points_gauss, weights_gauss;
    ComputeGaussLegendre(points_gauss, weights_gauss, order);
    G_GX.Reallocate(order+2, order+1);
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order; j++)
	G_GX(i, j) = lob_xy.EvaluatePhi(i, points_gauss(j));
    
  }
  
  
  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are located
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void WedgeHcurlOptimalHpFirstFamily
  ::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    //WedgeReference<2>::ComputeProjectionDofRef(feval, contrib);
    
    const VectReal_wp& weights3d = this->WeightsND();    
    Vector1 feval_weight = feval;
    for (int i = 0; i < feval_weight.GetM(); i+=3)
      {
       int j = i/3;
       feval_weight(i) *= weights3d(j);
       feval_weight(i+1) *= weights3d(j);
       feval_weight(i+2) *= weights3d(j);
     }
    
    contrib.Reallocate(nb_dof_loc);
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }
  

  void WedgeHcurlOptimalHpFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {
    bool tri = (num_loc%4 == 0);
    for (int i = 0; i < this->order; i++)
      if (i%2 == 0)
	{
	  contrib(2*this->order+i) = -contrib(2*this->order+i);
 	  if (!tri)
	    contrib(3*this->order+i) = -contrib(3*this->order+i);
	}
  }


  void WedgeHcurlOptimalHpFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    bool tri = (num_loc%4 == 0);
    for (int i = 0; i < this->order; i++)
      if (i%2 == 0)
	{
	  contrib(2*this->order+i) = -contrib(2*this->order+i);
	  if (!tri)
	    contrib(3*this->order+i) = -contrib(3*this->order+i);
	}
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
  void WedgeHcurlOptimalHpFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    if (num_loc%4 != 0)
      {
	const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
	Vector1 fx(3*(order+1)*(order+2)), fy(3*(order+2)*(order+2));
	FillZero(fx); FillZero(fy);
	// we express feval on Gauss-Lobatto points
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int node = 3*(i*(order+1) + j);
	      for (int m = 0; m <= order; m++)
		{
		  int p = 3*NumQuad2D_gauss(m, j);
		  fx(node) += G_GX(i, m)*feval(p);
		  fx(node+1) += G_GX(i, m)*feval(p+1);
		  fx(node+2) += G_GX(i, m)*feval(p+2);
		}
	    }
	
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      int node = 3*NumQuad2D(i, j);
	      for (int m = 0; m <= order; m++)
		{
		  int p = 3*(i*(order+1) + m);
		  fy(node) += G_GX(j, m)*fx(p);
		  fy(node+1) += G_GX(j, m)*fx(p+1);
		  fy(node+2) += G_GX(j, m)*fx(p+2);
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
    \param[out] phi values of basis functions on point_loc
  */
  void WedgeHcurlOptimalHpFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    int node = 0;
        
    VectReal_wp P1, P2, P3, Px, Py, Pz, Jx, Jy, Jz;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*x-1, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*y-1, Py);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*z-1, Pz);

    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*x-1, Jx);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*y-1, Jy);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*z-1, Jz);

    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*x+y-1, P1);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, y-x, P2);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*y+x-1, P3);
    
    // edges
    // triangle 1 - edge 1
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = (1-z)*(1-y)*P1(i);
        phi(node)(1) = (1-z)*x*P1(i);
        phi(node)(2) = 0;
        node++;
      }
    
    // triangle 1 - edge 2
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = -(1-z)*y*P2(i);
        phi(node)(1) = (1-z)*x*P2(i);
        phi(node)(2) = 0;
        node++;
      }
	  
    // triangle 1 - edge 3
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = (1-z)*y*P3(i);
        phi(node)(1) = (1-z)*(1-x)*P3(i);
        phi(node)(2) = 0;
        node++;
      }
    
    // vertical edge 1
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = 0;
        phi(node)(1) = 0;
        phi(node)(2) = (1-x-y)*Pz(i);
        node++;
      }
    
    // vertical edge 2
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = 0;
        phi(node)(1) = 0;
        phi(node)(2) = x*Pz(i);
        node++;
      }
    
    // vertical edge 3
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = 0;
        phi(node)(1) = 0;
        phi(node)(2) = y*Pz(i);
        node++;
      }
	  
    // triangle 2 - edge 1
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = z*(1-y)*P1(i);
        phi(node)(1) = z*x*P1(i);
        phi(node)(2) = 0;
        node++;
      }
    
    // triangle 2 - edge 2
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = -z*y*P2(i);
        phi(node)(1) = z*x*P2(i);
        phi(node)(2) = 0;
        node++;
      }
    
    // triangle 2 - edge 3
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = z*y*P3(i);
        phi(node)(1) = z*(1-x)*P3(i);
        phi(node)(2) = 0;
        node++;
      }
    
    const TriangleHcurlOptimalHpFirstFamily& Fb_tri
      = static_cast<const TriangleHcurlOptimalHpFirstFamily& >(*element_tri_surf);
    
    const Matrix<int>& NumDofsX_tri = Fb_tri.GetNumDofsX();
    const Matrix<int>& NumDofsY_tri = Fb_tri.GetNumDofsY();
    
    const QuadrangleHcurlOptimalHpFirstFamily& Fb_quad
      = static_cast<const QuadrangleHcurlOptimalHpFirstFamily& >(*element_quad_surf);
    
    const Matrix<int>& NumDofsX_quad = Fb_quad.GetNumDofsX();
    const Matrix<int>& NumDofsY_quad = Fb_quad.GetNumDofsY();

    Real_wp Sxiyj;
    // face triangulaire 1
    int offset = 6*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {		
          Sxiyj = P1(i)*P3(j);
	  
	  node = offset + NumDofsX_tri(i, j);
          phi(node)(0) = (1-z)*x*y*Sxiyj;
          phi(node)(1) = (1-z)*x*(1.0-x)*Sxiyj;
          phi(node)(2) = 0;
	  
	  node = offset + NumDofsY_tri(i, j);
          phi(node)(0) = (1-z)*y*(1.0-y)*Sxiyj;
          phi(node)(1) = (1-z)*y*x*Sxiyj;
          phi(node)(2) = 0;
	}	
    
    
    Real_wp Sxizj, Syizj;
    // face quadrangulaire 1
    offset += Fb_tri.GetNbDof() - 4*order;
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1; j++)
        {            
          Sxizj = Px(i)*Jz(j);
          
	  node = offset + NumDofsX_quad(i, j);
          phi(node)(0) = z*(1-z)*(1-y)*Sxizj;
          phi(node)(1) = z*(1-z)*x*Sxizj;
          phi(node)(2) = 0;
          
          Sxizj = Jx(i)*Pz(j);
          
	  node = offset + NumDofsY_quad(i, j);
          phi(node)(0) = 0;
          phi(node)(1) = 0;
          phi(node)(2) = (1-x-y)*x*Sxizj;
        }
    
    // face quadrangulaire 2
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1; j++)
        {  
          Syizj = Py(i)*Jz(j);
          
	  node = offset + NumDofsX_quad(i, j);
          phi(node)(0) = -z*(1-z)*y*Syizj;
          phi(node)(1) = z*(1-z)*x*Syizj;
          phi(node)(2) = 0;
			  
          Syizj = Jy(i)*Pz(j);
          
	  node = offset + NumDofsY_quad(i, j);
          phi(node)(0) = 0;
          phi(node)(1) = 0;
          phi(node)(2) = x*y*Syizj;
        }
    
    // face quadrangulaire 3
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1; j++)
        {  
          Syizj = Py(i)*Jz(j);
          
	  node = offset + NumDofsX_quad(i, j);
          phi(node)(0) = z*(1-z)*y*Syizj;
          phi(node)(1) = z*(1-z)*(1-x)*Syizj;
          phi(node)(2) = 0;
          
          Syizj = Jy(i)*Pz(j);
          
	  node = offset + NumDofsY_quad(i, j);
          phi(node)(0) = 0;
          phi(node)(1) = 0;
          phi(node)(2) = (1-x-y)*y*Syizj;
        }
    
    // face triangulaire 2
    offset += Fb_quad.GetNbDof() - 3*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {		
          Sxiyj = P1(i)*P3(j);
          
	  node = offset + NumDofsX_tri(i, j);
          phi(node)(0) = z*x*y*Sxiyj;
          phi(node)(1) = z*x*(1.0-x)*Sxiyj;
          phi(node)(2) = 0;
          
	  node = offset + NumDofsY_tri(i, j);
          phi(node)(0) = z*y*(1.0-y)*Sxiyj;
          phi(node)(1) = z*y*x*Sxiyj;
          phi(node)(2) = 0;
	}
    
    Real_wp Sijk;
    // interior
    for (int k = 0; k <= order-1; k++)
      for (int i = 0; i <= order-2; i++)
        for (int j = 0; j <= order-2; j++)
          if (i+j <= order-2)
            {
              Sijk = Px(i)*Py(j)*Pz(k);
	      
	      node = NumDofsX(i, j, k);
              phi(node)(0) = -z*(1-z)*x*y*Sijk;
              phi(node)(1) = -z*(1-z)*x*(1-x)*Sijk;
              phi(node)(2) = 0;
              
	      node = NumDofsY(i, j, k);
              phi(node)(0) = -(1-x-y)*z*(1-z)*y*Sijk;
              phi(node)(1) = (1-x-y)*z*(1-z)*x*Sijk;
              phi(node)(2) = 0;
              
	      node = NumDofsZ(i, j, k);
              phi(node)(0) = 0;
              phi(node)(1) = 0;
              phi(node)(2) = x*y*(1-x-y)*Sijk;
	    }	
    
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void WedgeHcurlOptimalHpFirstFamily::ComputeCurlPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    int node = 0;
    
    VectReal_wp Px, Py, Pz, dPx, dPy, dPz, P1, P2, P3, dP1, dP2, dP3;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*x-1, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*y-1, Py, dPy);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*z-1, Pz, dPz);
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*x+y-1, P1, dP1);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, y-x, P2, dP2);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*y+x-1, P3, dP3);
    
    VectReal_wp Jx, Jy, Jz, dJx, dJy, dJz;
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*x-1, Jx, dJx);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*y-1, Jy, dJy);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*z-1, Jz, dJz);
    
    // edges
    // triangle 1 - edge 1
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = x*P1(i);
        res(node)(1) = -(1-y)*P1(i);
        res(node)(2) = -(1-z)*(dP1(i)-2*x*dP1(i)-2*P1(i)-y*dP1(i));
        node++;
      }
    
    // triangle 1 - edge 2
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = x*P2(i);
        res(node)(1) = y*P2(i);
        res(node)(2) = -(1-z)*(x*dP2(i)-2*P2(i)-y*dP2(i));
        node++;
      }
    
    // triangle 1 - edge 3
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = (1-x)*P3(i);
        res(node)(1) = -y*P3(i);
        res(node)(2) = -(1-z)*(-dP3(i)+2*y*dP3(i)+x*dP3(i)+2*P3(i));
        node++;
      }
    
    // vertical edge 1
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = -Pz(i);
        res(node)(1) = Pz(i);
        res(node)(2) = 0;
        node++;
      }
    
    // vertical edge 2
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = 0;
        res(node)(1) = -Pz(i);
        res(node)(2) = 0;
        node++;
      }
    
    // vertical edge 3
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = Pz(i);
        res(node)(1) = 0;
        res(node)(2) = 0;
        node++;
      }
    
    // triangle 2 - edge 1
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = -x*P1(i);
        res(node)(1) = (1-y)*P1(i);
        res(node)(2) = -z*(dP1(i)-2*x*dP1(i)-2*P1(i)-y*dP1(i));
        node++;
      }
    
    // triangle 2 - edge 2
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = -x*P2(i);
        res(node)(1) = -y*P2(i);
        res(node)(2) = -z*(x*dP2(i)-2*P2(i)-y*dP2(i));
        node++;
      }
    
    // triangle 2 - edge 3
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = -(1-x)*P3(i);
        res(node)(1) = y*P3(i);
        res(node)(2) = -z*(-dP3(i)+2*y*dP3(i)+x*dP3(i)+2*P3(i));
        node++;
      }

    const TriangleHcurlOptimalHpFirstFamily& Fb_tri
      = static_cast<const TriangleHcurlOptimalHpFirstFamily& >(*element_tri_surf);
    
    const Matrix<int>& NumDofsX_tri = Fb_tri.GetNumDofsX();
    const Matrix<int>& NumDofsY_tri = Fb_tri.GetNumDofsY();
    
    const QuadrangleHcurlOptimalHpFirstFamily& Fb_quad
      = static_cast<const QuadrangleHcurlOptimalHpFirstFamily& >(*element_quad_surf);
    
    const Matrix<int>& NumDofsX_quad = Fb_quad.GetNumDofsX();
    const Matrix<int>& NumDofsY_quad = Fb_quad.GetNumDofsY();
    
    Real_wp Sxiyj; R3 grad_zPiPj, vec_u, grad_lambda, curl_u;
    Real_wp lambda;
    // face triangulaire 1
    int offset = 6*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        { 
          Sxiyj = P1(i)*P3(j);
          grad_zPiPj(0) = (1.0-z)*(2.0*dP1(i)*P3(j) + P1(i)*dP3(j));
          grad_zPiPj(1) = (1.0-z)*(dP1(i)*P3(j) + 2.0*P1(i)*dP3(j));
          grad_zPiPj(2) = -Sxiyj;
          
          grad_lambda = x*grad_zPiPj;
          grad_lambda(0) += (1.0-z)*Sxiyj;
          
          lambda = x*(1.0-z)*Sxiyj;
          vec_u.Init(y, 1.0-x, 0);
          
          TimesProd(grad_lambda, vec_u, curl_u);
          
	  node = offset + NumDofsX_tri(i, j);
          res(node)(0) = curl_u(0);
          res(node)(1) = curl_u(1);
          res(node)(2) = -2.0*lambda + curl_u(2);
          
          
          grad_lambda = y*grad_zPiPj;
          grad_lambda(1) += (1.0-z)*Sxiyj;
          
          lambda = y*(1.0-z)*Sxiyj;
          vec_u.Init(1.0-y, x, 0);
          
          TimesProd(grad_lambda, vec_u, curl_u);
          
          node = offset + NumDofsY_tri(i, j);
          res(node)(0) = curl_u(0);
          res(node)(1) = curl_u(1);
          res(node)(2) = 2.0*lambda + curl_u(2);
          
          node++;
        }		 
    
    Real_wp  Sxizj, Syizj;
    // face quadrangulaire 1
    offset += Fb_tri.GetNbDof() - 4*order;
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1; j++)
        {  
          Sxizj = Px(i)*Jz(j);
          
          lambda = z*(1.0-z)*Sxizj;
          grad_lambda.Init(2.0*z*(1.0-z)*dPx(i)*Jz(j),
                           0, (1.0-2.0*z)*Sxizj + 2.0*z*(1.0-z)*Px(i)*dJz(j));
          
          vec_u.Init((1.0-y), x, 0);
          TimesProd(grad_lambda, vec_u, curl_u);
          
	  node = offset + NumDofsX_quad(i, j);
          res(node)(0) = curl_u(0);
          res(node)(1) = curl_u(1);
          res(node)(2) = 2.0*lambda + curl_u(2);
          
          Sxizj = Jx(i)*Pz(j);
          
          grad_lambda.Init((1.0-2.0*x-y)*Sxizj + 2.0*(1.0-x-y)*x*dJx(i)*Pz(j),
                           -x*Sxizj, 0);
          
	  node = offset + NumDofsY_quad(i, j);
          res(node)(0) = grad_lambda(1);
          res(node)(1) = -grad_lambda(0);
          res(node)(2) = 0;
        }
    
    // face quadrangulaire 2
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1; j++)
        {  
          Syizj = Py(i)*Jz(j);
          
          lambda = z*(1.0-z)*Syizj;
          grad_lambda.Init(0, 2.0*z*(1.0-z)*dPy(i)*Jz(j),
                           (1.0-2.0*z)*Syizj + 2.0*z*(1.0-z)*Py(i)*dJz(j));
          
          vec_u.Init(-y, x, 0);
          TimesProd(grad_lambda, vec_u, curl_u);
          
	  node = offset + NumDofsX_quad(i, j);
          res(node)(0) = curl_u(0);
          res(node)(1) = curl_u(1);
          res(node)(2) = 2.0*lambda + curl_u(2);          
          
          Syizj = Jy(i)*Pz(j);
          grad_lambda.Init(y*Syizj, x*Syizj + 2.0*x*y*dJy(i)*Pz(j), 0);
          
	  node = offset + NumDofsY_quad(i, j);
          res(node)(0) = grad_lambda(1);
          res(node)(1) = -grad_lambda(0);
          res(node)(2) = 0;
        }
    
    // face quadrangulaire 3
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1; j++)
        {  
          Syizj = Py(i)*Jz(j);
          
          lambda = z*(1.0-z)*Syizj;
          grad_lambda.Init(0, 2.0*z*(1.0-z)*dPy(i)*Jz(j),
                           (1.0-2.0*z)*Syizj + 2.0*z*(1.0-z)*Py(i)*dJz(j));
          
          vec_u.Init(y, 1.0-x, 0);
          TimesProd(grad_lambda, vec_u, curl_u);
          
	  node = offset + NumDofsX_quad(i, j);
          res(node)(0) = curl_u(0);
          res(node)(1) = curl_u(1);
          res(node)(2) = -2.0*lambda + curl_u(2);
          
          Syizj = Jy(i)*Pz(j);
          grad_lambda.Init(-y*Syizj, (1.0-x-2.0*y)*Syizj + 2.0*(1.0-x-y)*y*dJy(i)*Pz(j), 0);
          
	  node = offset + NumDofsY_quad(i, j);
          res(node)(0) = grad_lambda(1);
          res(node)(1) = -grad_lambda(0);
          res(node)(2) = 0;
        }
    
    // face triangulaire 2
    offset += Fb_quad.GetNbDof() - 3*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {		
          Sxiyj = P1(i)*P3(j);
          grad_zPiPj(0) = z*(2.0*dP1(i)*P3(j) + P1(i)*dP3(j));
          grad_zPiPj(1) = z*(dP1(i)*P3(j) + 2.0*P1(i)*dP3(j));
          grad_zPiPj(2) = Sxiyj;
          
          grad_lambda = x*grad_zPiPj;
          grad_lambda(0) += z*Sxiyj;
          
          lambda = x*z*Sxiyj;
          vec_u.Init(y, 1.0-x, 0);
          
          TimesProd(grad_lambda, vec_u, curl_u);
          
	  node = offset + NumDofsX_tri(i, j);
          res(node)(0) = curl_u(0);
          res(node)(1) = curl_u(1);
          res(node)(2) = -2.0*lambda + curl_u(2);
          
          grad_lambda = y*grad_zPiPj;
          grad_lambda(1) += z*Sxiyj;
          
          lambda = y*z*Sxiyj;
          vec_u.Init(1.0-y, x, 0);
          
          TimesProd(grad_lambda, vec_u, curl_u);
	  
	  node = offset + NumDofsY_tri(i, j);
          res(node)(0) = curl_u(0);
          res(node)(1) = curl_u(1);
          res(node)(2) = 2.0*lambda + curl_u(2);
        }
    
    // interior
    for (int k = 0; k <= order-1; k++)
      for (int i = 0; i <= order-2; i++)
        for (int j = 0; j <= order-2; j++)
          if (i+j <= order-2)
            {
              Sxiyj = Px(i)*Py(j);
	      
	      node = NumDofsX(i, j, k);
              res(node)(0) = (-1+x)*x*Px(i)*Py(j)*(-Pz(k)+2*z*Pz(k)+2*pow(z,2)*dPz(k)-2*z*dPz(k));
              res(node)(1) = y*x*Px(i)*Py(j)*(-Pz(k)+2*z*Pz(k)+2*pow(z,2)*dPz(k)-2*z*dPz(k));
              res(node)(2) = -(-1+z)*z*Pz(k)*(-Px(i)*Py(j)-2*x*dPx(i)*Py(j)
                                              +2*pow(x,2)*dPx(i)*Py(j)+2*y*x*dPy(j)*Px(i)
                                              +3*x*Px(i)*Py(j));
              
              node = NumDofsY(i, j, k);
              res(node)(0) = (1-x-y)*x*Px(i)*Py(j)*(-Pz(k)+2*z*Pz(k)+2*pow(z,2)*dPz(k)-2*z*dPz(k));
              res(node)(1) = (1-x-y)*y*Px(i)*Py(j)*(-Pz(k)+2*z*Pz(k)+2*pow(z,2)*dPz(k)-2*z*dPz(k));
              res(node)(2) = -(1-z)*z*Pz(k)*(2*pow(y,2)*dPy(j)*Px(i)+3*x*Px(i)*Py(j)
                                             +3*y*Px(i)*Py(j)-2*Px(i)*Py(j)-2*x*dPx(i)*Py(j)
                                             +2*y*x*dPy(j)*Px(i)+2*pow(x,2)*dPx(i)*Py(j)
                                             +2*x*y*dPx(i)*Py(j)-2*y*dPy(j)*Px(i));
	      
              node = NumDofsZ(i, j, k);
              res(node)(0) = -x*Px(i)*Pz(k)*(2*x*y*dPy(j)+x*Py(j)-2*y*dPy(j)
                                             +2*pow(y,2)*dPy(j)-Py(j)+2*y*Py(j));
              res(node)(1) = y*Py(j)*Pz(k)*(2*x*y*dPx(i)+2*pow(x,2)*dPx(i)+2*x*Px(i)-Px(i)
                                            -2*x*dPx(i)+y*Px(i));
              res(node)(2) = 0;
            }	
    
  }
  

  //! displays details of class WedgeHcurlOptimalHpFirstFamily
  ostream& operator <<(ostream& out, const WedgeHcurlOptimalHpFirstFamily& e)
  {
    out<<static_cast<const WedgeReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_WEDGE_HCURL_OPTIMAL_HP_FIRST_FAMILY_CXX
#endif
