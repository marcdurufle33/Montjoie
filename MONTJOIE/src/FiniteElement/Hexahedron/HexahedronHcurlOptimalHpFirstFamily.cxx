#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_OPTIMAL_HP_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronHcurlOptimalHpFirstFamily::HexahedronHcurlOptimalHpFirstFamily()
    : HexahedronReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  

  //! returns size of memory used by the object
  size_t HexahedronHcurlOptimalHpFirstFamily::GetMemorySize() const
  {
    size_t taille = HexahedronReference<2>::GetMemorySize();
    taille += NumDofsX.GetMemorySize() + NumDofsY.GetMemorySize() + NumDofsZ.GetMemorySize();
    taille += NumQuad2D_gauss.GetMemorySize() + G_GX.GetMemorySize();
    taille += LegendrePolynom.GetMemorySize() + JacobiPolynom.GetMemorySize();    
    taille += ValLeg.GetMemorySize() + ValBubbleJac.GetMemorySize() 
      + DerBubbleJac.GetMemorySize();

    taille += ch1_node.GetMemorySize() + ch2_node.GetMemorySize() + ch3_node.GetMemorySize()
      + ch1_loc.GetMemorySize() + ch2_loc.GetMemorySize() + ch3_loc.GetMemorySize()
      + rh_loc.GetMemorySize();
    
    for (int k = 0; k < 6; k++)
      {
        taille += sh_loc(k).GetMemorySize() + sh_locX(k).GetMemorySize() + sh_locY(k).GetMemorySize()
          + sh_curl_loc(k).GetMemorySize() + sh_curl_locX(k).GetMemorySize() + sh_curl_locY(k).GetMemorySize()
          + sh_nodeX(k).GetMemorySize() + sh_nodeY(k).GetMemorySize() +
          sh_curl_nodeX(k).GetMemorySize() + sh_curl_nodeY(k).GetMemorySize();
      }
    
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        {
          taille += sparse_mass_matrix(i, j).GetMemorySize() +
            sparse_stiff_matrix(i, j).GetMemorySize();
        }
    
    taille += sparse_matrix_chol.GetMemorySize();
    
    return taille;
  }

  
  //! how to number mesh
  void HexahedronHcurlOptimalHpFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return HexahedronReference<2>::ConstructNumberMap(nmap, dg);

    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofEdge(order, order);
    // nothing on quadrangular face !
    nmap.SetNbDofQuadrangle(order, 2*order*order);
    
    // dofs inside the tetrahedron
    nmap.SetNbDofHexahedron(order, nb_dof_loc-nb_dof_boundaries);

    nmap.SetEqualEdgesDofSymmetry(order, order+1);
   
    element_quad_surf->FindHcurlRotationQuad(nmap);
  }

  
  //! construction of finite element
  void HexahedronHcurlOptimalHpFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf_tri, int rsurf_quad,
			   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (rquad <= r)
      rquad = r+1;
    
    // type_quad = HexahedronQuadrature::QUADRATURE_JACOBI1;    
    HexahedronReference<2>::
      ConstructFiniteElement(r, rgeom, rquad, Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    // computation of basis functions
    ConstructFunctions();
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructMassMatrix();
    ConstructStiffnessMatrix();

    /* VectR3 phi;
    this->MatProjectionDof.Reallocate(this->GetNbPointsDof(), this->nb_dof_loc);
    for (int i = 0; i < this->GetNbPointsDof(); i++)
      {
	this->ComputeValuesPhiRef(this->PointsDofND(i), phi);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  this->MatProjectionDof(i, j) = DotProd(phi(j), this->tangente_dof(i));
      }
    
      GetQR(this->MatProjectionDof, this->tauProjectionDof);*/
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace(false);

  }
  
  
  //! construction of basis functions
  void HexahedronHcurlOptimalHpFirstFamily::ConstructFunctions()
  {
    nb_dof_loc = 3*order*(order+2)*(order+2);
    nb_dof_boundaries = 12*order + 12*order*order;
    nb_dof_quad = 2*order*(order+2);
    nb_dof_tri = 0;
    
    VectReal_wp points_gauss, poids_gauss, points_lob, poids_lob;
    ComputeGaussLegendre(points_gauss, poids_gauss, order-1);
    ComputeGaussLobatto(points_lob, poids_lob, order+1);
    
    VectR3 points_dof3d;
    points_dof3d = this->PointsND();
    
    /* points_dof3d.Reallocate(nb_dof_loc);
    tangente_dof.Reallocate(nb_dof_loc);
    int nb = 0;
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
	for (int k = 0; k <= order+1; k++)
	  {
	    points_dof3d(nb).Init(points_gauss(i), points_lob(j), points_lob(k));
	    tangente_dof(nb).Init(1.0, 0, 0);
	    points_dof3d(nb+1).Init(points_lob(j), points_gauss(i), points_lob(k));
	    tangente_dof(nb+1).Init(0, 1.0, 0);
	    points_dof3d(nb+2).Init(points_lob(j), points_lob(k), points_gauss(i));
	    tangente_dof(nb+2).Init(0.0, 0, 1.0);
	    nb += 3;
	  }
    */
    
    this->SetPointsDof2D_quad(this->Points2D_quad());
    this->SetPointsDofND(points_dof3d);

    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;

    GetJacobiPolynomial(LegendrePolynom, order, Real_wp(0), Real_wp(0));
    GetJacobiPolynomial(JacobiPolynom, order, Real_wp(1), Real_wp(1));

    QuadrangleHcurlOptimalHpFirstFamily* Fb_quad = new QuadrangleHcurlOptimalHpFirstFamily();
    Fb_quad->ConstructFiniteElement(order);

    element_quad_surf = Fb_quad;
    
    NumDofsX.Reallocate(order, order+2, order+2);
    NumDofsY.Reallocate(order+2, order, order+2);
    NumDofsZ.Reallocate(order+2, order+2, order);
    
    // dofs on edges
    for (int i = 0; i < order; i++)
      {
        NumDofsX(i, 0, 0) = i;
        NumDofsY(order+1, i, 0) = order + i;
        NumDofsX(i, order+1, 0) = 2*order + i;
        NumDofsY(0, i, 0) = 3*order + i;
        NumDofsZ(0, 0, i) = 4*order + i;
        NumDofsZ(order+1, 0, i) = 5*order + i;
        NumDofsZ(order+1, order+1, i) = 6*order + i;
        NumDofsZ(0, order+1, i) = 7*order + i;
        NumDofsX(i, 0, order+1) = 8*order + i;
        NumDofsY(order+1, i, order+1) = 9*order + i;
        NumDofsX(i, order+1, order+1) = 10*order + i;
        NumDofsY(0, i, order+1) = 11*order + i;
      }
    
    const Matrix<int>& NumDofsX_quad = Fb_quad->GetNumDofsX();
    const Matrix<int>& NumDofsY_quad = Fb_quad->GetNumDofsY();

    // dofs on faces
    int N = Fb_quad->GetNbDof();
    for (int i = 0; i < order; i++)
      for (int j = 0; j < order; j++)
        {
          int node = 8*order + NumDofsX_quad(i, j);
          NumDofsY(0, i, j+1) = node;
          node = 8*order + NumDofsY_quad(j, i);
          NumDofsZ(0, j+1, i) = node;
          
          node = 4*order + N + NumDofsX_quad(i, j);
          NumDofsX(i, 0, j+1) = node;
          node = 4*order + N + NumDofsY_quad(j, i);
          NumDofsZ(j+1, 0, i) = node;
          
          node = 2*N + NumDofsX_quad(i, j);
          NumDofsX(i, j+1, 0) = node;
          node = 2*N + NumDofsY_quad(j, i);
          NumDofsY(j+1, i, 0) = node;

          node = -4*order + 3*N + NumDofsX_quad(i, j);
          NumDofsX(i, j+1, order+1) = node;
          node = -4*order + 3*N + NumDofsY_quad(j, i);
          NumDofsY(j+1, i, order+1) = node;

          node = -8*order + 4*N + NumDofsX_quad(i, j);
          NumDofsX(i, order+1, j+1) = node;
          node = -8*order + 4*N + NumDofsY_quad(j, i);
          NumDofsZ(j+1, order+1, i) = node;

          node = -12*order + 5*N + NumDofsX_quad(i, j);
          NumDofsY(order+1, i, j+1) = node;
          node = -12*order + 5*N + NumDofsY_quad(j, i);
          NumDofsZ(order+1, j+1, i) = node;
        }
    
    // dofs inside the hexahedron
    int offset = nb_dof_boundaries;
    for (int diag = 0; diag <= order-1; diag++)
      for (int i = 0; i <= diag; i++)
	for (int j = 0; j <= diag; j++)
	  for (int k = 0; k <= diag; k++)
	    if ((i == diag) || (j == diag) || (k == diag))
	      {
		NumDofsX(i, j+1, k+1) = offset++;
		NumDofsY(i+1, j, k+1) = offset++;
		NumDofsZ(i+1, j+1, k) = offset++;
	      }

    Matrix<int> coor;
    MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(order, NumQuad2D_gauss, coor);
    Globatto<Real_wp> lob_xy;
    lob_xy.ConstructQuadrature(order+1, lob_xy.QUADRATURE_GAUSS);
    
    VectReal_wp weights_gauss;
    ComputeGaussLegendre(points_gauss, weights_gauss, order);
    G_GX.Reallocate(order+2, order+1);
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order; j++)
	G_GX(i, j) = lob_xy.EvaluatePhi(i, points_gauss(j));
    
    // storing P_j^(0,0) ( \xi_i)  and (1-x)*x*P_j^{1,1}(\xi_i)
    int r = order_quad;
    ValLeg.Reallocate(order, r+3);
    ValBubbleJac.Reallocate(order+2, r+3);
    DerBubbleJac.Reallocate(order+2, r+3);
    ValBubbleJac.Fill(0);
    DerBubbleJac.Fill(0);
    VectReal_wp Px, Jx, dJx;
    const VectReal_wp& points1d = this->Points1D();
    for (int i = 0; i <= r; i++)
      {
        EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*points1d(i)-1.0, Px);
        EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2.0*points1d(i)-1.0, Jx, dJx);
        for (int j = 0; j < order; j++)
          ValLeg(j, i) = Px(j);
        
        ValBubbleJac(0, i) = 1.0-points1d(i); DerBubbleJac(0, i) = -1.0;
        ValBubbleJac(order+1, i) = points1d(i); DerBubbleJac(order+1, i) = 1.0;
        for (int j = 0; j < order; j++)
          {
            ValBubbleJac(j+1, i) = Jx(j)*points1d(i)*(1.0-points1d(i));
            DerBubbleJac(j+1, i)
              = (1.0 - 2.0*points1d(i))*Jx(j) + 2.0*dJx(j)*points1d(i)*(1.0-points1d(i));
          }        
      }
    
    // values on extremities are stored at the end
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, Real_wp(-1), Px);
    for (int j = 0; j < order; j++)
      ValLeg(j, r+1) = Px(j);
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, Real_wp(1), Px);
    for (int j = 0; j < order; j++)
      ValLeg(j, r+2) = Px(j);
    
    ValBubbleJac(0, r+1) = 1.0; ValBubbleJac(order+1, r+1) = 0.0;
    DerBubbleJac(0, r+1) = -1.0; DerBubbleJac(order+1, r+1) = 1.0;
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, Real_wp(-1), Jx);
    for (int j = 0; j < order; j++)
      DerBubbleJac(j+1, r+1) = Jx(j);
    
    ValBubbleJac(0, r+2) = 0.0; ValBubbleJac(order+1, r+2) = 1.0;
    DerBubbleJac(0, r+2) = -1.0; DerBubbleJac(order+1, r+2) = 1.0;
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, Real_wp(1), Jx);
    for (int j = 0; j < order; j++)
      DerBubbleJac(j+1, r+2) = -Jx(j);
    
    // adding quadrature points of the face in CoordinateQuad3D
    CoordinateQuad3D.Resize(this->GetNbPointsQuadrature(), 3);
    N = nb_points_quadrature_inside;
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        {
          int p = N + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = r+1;
          CoordinateQuad3D(p, 1) = i;
          CoordinateQuad3D(p, 2) = j;

          p = N + (r+1)*(r+1) + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = i;
          CoordinateQuad3D(p, 1) = r+1;
          CoordinateQuad3D(p, 2) = j;

          p = N + 2*(r+1)*(r+1) + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = i;
          CoordinateQuad3D(p, 1) = j;
          CoordinateQuad3D(p, 2) = r+1;


          p = N + 3*(r+1)*(r+1) + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = i;
          CoordinateQuad3D(p, 1) = j;
          CoordinateQuad3D(p, 2) = r+2;


          p = N + 4*(r+1)*(r+1) + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = i;
          CoordinateQuad3D(p, 1) = r+2;
          CoordinateQuad3D(p, 2) = j;

          p = N + 5*(r+1)*(r+1) + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = r+2;
          CoordinateQuad3D(p, 1) = i;
          CoordinateQuad3D(p, 2) = j;
        }

  }
  
  
  //! constructing mass matrix
  void HexahedronHcurlOptimalHpFirstFamily::ConstructMassMatrix()
  {
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2, Ch3;
    ComputeProjectionMatrix(Ch1, Ch2, Ch3, NumQuad3D, this->Points1D(),
                            ValLeg, ValBubbleJac);

    // conversion to csr format for more efficiency
    Copy(Ch1, ch1_loc);
    Copy(Ch2, ch2_loc);
    Copy(Ch3, ch3_loc);
    
    int order_geom = this->GetGeometryOrder();
    const VectReal_wp& points_nodal1d = this->PointsNodal1D();
    
    // now projection on nodes
    Matrix<Real_wp> ValLegNode(order, order_geom+1);
    Matrix<Real_wp> ValBubbleJacNode(order+2, order_geom+1);
    Matrix<Real_wp> DerBubbleJacNode(order+2, order_geom+1);
    VectReal_wp Px, Jx, dJx;
    for (int i = 0; i <= order_geom; i++)
      {
        EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*points_nodal1d(i)-1.0, Px);
        EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2.0*points_nodal1d(i)-1.0, Jx, dJx);
        for (int j = 0; j < order; j++)
          ValLegNode(j, i) = Px(j);
        
        ValBubbleJacNode(0, i) = 1.0-points_nodal1d(i);
        ValBubbleJacNode(order+1, i) = points_nodal1d(i);
        DerBubbleJacNode(0, i) = -1.0;
        DerBubbleJacNode(order+1, i) = 1.0;
        for (int j = 0; j < order; j++)
          {
            ValBubbleJacNode(j+1, i) = Jx(j)*points_nodal1d(i)*(1.0-points_nodal1d(i));
            DerBubbleJacNode(j+1, i) = (1.0-2.0*points_nodal1d(i))*Jx(j)
              + 2.0*dJx(j)*points_nodal1d(i)*(1.0-points_nodal1d(i));
          }
      }
    
    ComputeProjectionMatrix(Ch1, Ch2, Ch3, Fb_geom.GetNumNodes3D(), points_nodal1d,
                            ValLegNode, ValBubbleJacNode);
    
    Copy(Ch1, ch1_node);
    Copy(Ch2, ch2_node);
    Copy(Ch3, ch3_node);
    
    // projection on surface nodal points
    TinyVector<Matrix<Real_wp, General, ArrayRowSparse>, 6> Sh;
    int r = order;
    int rq = order_quad;
    
    // projection on surface of Legendre polynomials
    for (int n = 0; n < 6; n++)
      Sh(n).Reallocate(nb_dof_loc, 3*(r+2)*(r+2));
    
    Real_wp threshold = 1e3*epsilon_machine;
    for (int k = 0; k <= order-1; k++)
      for (int i = 0; i <= order+1; i++)
        for (int j = 0; j <= order+1; j++)
          {
            int dof = NumDofsX(k, i, j);
            int node = i*(r+2) + j;
            Sh(0).AddInteraction(dof, 3*node, ValLeg(k, rq+1));
            Sh(5).AddInteraction(dof, 3*node, ValLeg(k, rq+2));
            if (abs(ValBubbleJac(i, rq+1)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(1).AddInteraction(dof, 3*node, ValBubbleJac(i, rq+1));
              }

            if (abs(ValBubbleJac(i, rq+2)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(4).AddInteraction(dof, 3*node, ValBubbleJac(i, rq+2));
              }

            if (abs(ValBubbleJac(j, rq+1)) > threshold)
              {
                node = k*(r+2) + i;
                Sh(2).AddInteraction(dof, 3*node, ValBubbleJac(j, rq+1));
              }

            if (abs(ValBubbleJac(j, rq+2)) > threshold)
              {
                node = k*(r+2) + i;
                Sh(3).AddInteraction(dof, 3*node, ValBubbleJac(j, rq+2));
              }
            
            dof = NumDofsY(i, k, j);
            node = i*(r+2) + j;
            Sh(1).AddInteraction(dof, 3*node+1, ValLeg(k, rq+1));
            Sh(4).AddInteraction(dof, 3*node+1, ValLeg(k, rq+2));
            if (abs(ValBubbleJac(i, rq+1)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(0).AddInteraction(dof, 3*node+1, ValBubbleJac(i, rq+1));
              }

            if (abs(ValBubbleJac(i, rq+2)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(5).AddInteraction(dof, 3*node+1, ValBubbleJac(i, rq+2));
              }

            if (abs(ValBubbleJac(j, rq+1)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(2).AddInteraction(dof, 3*node+1, ValBubbleJac(j, rq+1));
              }

            if (abs(ValBubbleJac(j, rq+2)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(3).AddInteraction(dof, 3*node+1, ValBubbleJac(j, rq+2));
              }

            dof = NumDofsZ(i, j, k);
            node = i*(r+2) + j;
            Sh(2).AddInteraction(dof, 3*node+2, ValLeg(k, rq+1));
            Sh(3).AddInteraction(dof, 3*node+2, ValLeg(k, rq+2));
            if (abs(ValBubbleJac(i, rq+1)) > threshold)
              {
                node = j*(r+2) + k;
                Sh(0).AddInteraction(dof, 3*node+2, ValBubbleJac(i, rq+1));
              }

            if (abs(ValBubbleJac(i, rq+2)) > threshold)
              {
                node = j*(r+2) + k;
                Sh(5).AddInteraction(dof, 3*node+2, ValBubbleJac(i, rq+2));
              }

            if (abs(ValBubbleJac(j, rq+1)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(1).AddInteraction(dof, 3*node+2, ValBubbleJac(j, rq+1));
              }

            if (abs(ValBubbleJac(j, rq+2)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(4).AddInteraction(dof, 3*node+2, ValBubbleJac(j, rq+2));
              }            
          }
    
    for (int n = 0; n < 6; n++)
      Copy(Sh(n), sh_loc(n));
    
    const Matrix<int>& NumNodes2D = Fb_geom.GetQuadrangularSurfaceFiniteElement().GetNumNodes2D();
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    ComputeSurfaceProjection(sh_locX, sh_locY, this->Points1D(), ValLeg, ValBubbleJac, NumQuad2D);
    ComputeSurfaceProjection(sh_nodeX, sh_nodeY, points_nodal1d,
                             ValLegNode, ValBubbleJacNode, NumNodes2D);
    
    // projection on surface of Legendre polynomials
    for (int n = 0; n < 6; n++)
      {
        Sh(n).Clear();
        Sh(n).Reallocate(nb_dof_loc, 5*(r+2)*(r+2));
      }
    
    for (int k = 0; k <= order-1; k++)
      for (int i = 0; i <= order+1; i++)
        for (int j = 0; j <= order+1; j++)
          {
            int dof = NumDofsX(k, i, j);
            int node = i*(r+2) + j;
            Sh(0).AddInteraction(dof, 5*node, ValLeg(k, rq+1));
            Sh(5).AddInteraction(dof, 5*node, ValLeg(k, rq+2));
            if (abs(ValBubbleJac(i, rq+1)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(1).AddInteraction(dof, 5*node, ValBubbleJac(i, rq+1));
              }

            if (abs(ValBubbleJac(i, rq+2)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(4).AddInteraction(dof, 5*node, ValBubbleJac(i, rq+2));
              }

            if (abs(ValBubbleJac(j, rq+1)) > threshold)
              {
                node = k*(r+2) + i;
                Sh(2).AddInteraction(dof, 5*node, ValBubbleJac(j, rq+1));
              }

            if (abs(ValBubbleJac(j, rq+2)) > threshold)
              {
                node = k*(r+2) + i;
                Sh(3).AddInteraction(dof, 5*node, ValBubbleJac(j, rq+2));
              }

            if (abs(DerBubbleJac(i, rq+1)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(1).AddInteraction(dof, 5*node+3, DerBubbleJac(i, rq+1));
              }

            if (abs(DerBubbleJac(i, rq+2)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(4).AddInteraction(dof, 5*node+3, DerBubbleJac(i, rq+2));
              }

            if (abs(DerBubbleJac(j, rq+1)) > threshold)
              {
                node = k*(r+2) + i;
                Sh(2).AddInteraction(dof, 5*node+3, DerBubbleJac(j, rq+1));
              }

            if (abs(DerBubbleJac(j, rq+2)) > threshold)
              {
                node = k*(r+2) + i;
                Sh(3).AddInteraction(dof, 5*node+3, DerBubbleJac(j, rq+2));
              }
            
            dof = NumDofsY(i, k, j);
            node = i*(r+2) + j;
            Sh(1).AddInteraction(dof, 5*node+1, ValLeg(k, rq+1));
            Sh(4).AddInteraction(dof, 5*node+1, ValLeg(k, rq+2));
            if (abs(ValBubbleJac(i, rq+1)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(0).AddInteraction(dof, 5*node+1, ValBubbleJac(i, rq+1));
              }

            if (abs(ValBubbleJac(i, rq+2)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(5).AddInteraction(dof, 5*node+1, ValBubbleJac(i, rq+2));
              }

            if (abs(ValBubbleJac(j, rq+1)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(2).AddInteraction(dof, 5*node+1, ValBubbleJac(j, rq+1));
              }

            if (abs(ValBubbleJac(j, rq+2)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(3).AddInteraction(dof, 5*node+1, ValBubbleJac(j, rq+2));
              }
            
            if (abs(DerBubbleJac(i, rq+1)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(0).AddInteraction(dof, 5*node+3, DerBubbleJac(i, rq+1));
              }

            if (abs(DerBubbleJac(i, rq+2)) > threshold)
              {
                node = k*(r+2) + j;
                Sh(5).AddInteraction(dof, 5*node+3, DerBubbleJac(i, rq+2));
              }

            if (abs(DerBubbleJac(j, rq+1)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(2).AddInteraction(dof, 5*node+4, DerBubbleJac(j, rq+1));
              }

            if (abs(DerBubbleJac(j, rq+2)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(3).AddInteraction(dof, 5*node+4, DerBubbleJac(j, rq+2));
              }
            
            dof = NumDofsZ(i, j, k);
            node = i*(r+2) + j;
            Sh(2).AddInteraction(dof, 5*node+2, ValLeg(k, rq+1));
            Sh(3).AddInteraction(dof, 5*node+2, ValLeg(k, rq+2));
            if (abs(ValBubbleJac(i, rq+1)) > threshold)
              {
                node = j*(r+2) + k;
                Sh(0).AddInteraction(dof, 5*node+2, ValBubbleJac(i, rq+1));
              }

            if (abs(ValBubbleJac(i, rq+2)) > threshold)
              {
                node = j*(r+2) + k;
                Sh(5).AddInteraction(dof, 5*node+2, ValBubbleJac(i, rq+2));
              }

            if (abs(ValBubbleJac(j, rq+1)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(1).AddInteraction(dof, 5*node+2, ValBubbleJac(j, rq+1));
              }

            if (abs(ValBubbleJac(j, rq+2)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(4).AddInteraction(dof, 5*node+2, ValBubbleJac(j, rq+2));
              }            

            if (abs(DerBubbleJac(i, rq+1)) > threshold)
              {
                node = j*(r+2) + k;
                Sh(0).AddInteraction(dof, 5*node+4, DerBubbleJac(i, rq+1));
              }

            if (abs(DerBubbleJac(i, rq+2)) > threshold)
              {
                node = j*(r+2) + k;
                Sh(5).AddInteraction(dof, 5*node+4, DerBubbleJac(i, rq+2));
              }

            if (abs(DerBubbleJac(j, rq+1)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(1).AddInteraction(dof, 5*node+4, DerBubbleJac(j, rq+1));
              }

            if (abs(DerBubbleJac(j, rq+2)) > threshold)
              {
                node = i*(r+2) + k;
                Sh(4).AddInteraction(dof, 5*node+4, DerBubbleJac(j, rq+2));
              }            
          }
    
    for (int n = 0; n < 6; n++)
      Copy(Sh(n), sh_curl_loc(n));
    
    ComputeSurfaceCurlProjection(sh_curl_locX, sh_curl_locY, this->Points1D(),
                                 DerBubbleJac, ValLeg, ValBubbleJac, NumQuad2D);

    ComputeSurfaceCurlProjection(sh_curl_nodeX, sh_curl_nodeY, points_nodal1d,
                                 DerBubbleJacNode, ValLegNode, ValBubbleJacNode, NumNodes2D);
    
  }
  

  //! computing the projection from degrees of freedom to a set of points
  /*!
    \param[out] Ch1 part of projection operator for interpolation along x
    \param[out] Ch2 part of projection operator for interpolation along y
    \param[out] Ch3 part of projection operator for interpolation along z
    \param[in] Num numbering of final points
    \param[in] points 1-D points, we compute projection on (points_i, points_j, points_k)
    \param[in] ValLegendre ValLegendre(i, j) = P_i^{0,0}(2 points_j - 1)
    \param[in] ValBubble ValBubble(i, j) = (1-x) if i = 0, x if i = 1,
                                    x (1 - x) P_i^{1,1}(2 x - 1) otherwise with x = points_j
    The operator Ch1^* Ch2^* Ch3^* projects degrees freedom into values on specified points    
   */
  void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeProjectionMatrix(Matrix<Real_wp, General, ArrayRowSparse>& Ch1,
                            Matrix<Real_wp, General, ArrayRowSparse>& Ch2,
                            Matrix<Real_wp, General, ArrayRowSparse>& Ch3,
                            const Array3D<int>& Num, const VectReal_wp& points,
                            Matrix<Real_wp>& ValLegendre, Matrix<Real_wp>& ValBubble)
  {
    int r = order;
    int rq = points.GetM() - 1;
    int N = r*(r+2);
    
    Ch3.Reallocate(nb_dof_loc, 3*r*(r+2)*(rq+1));
    // first step, from dofs, we compute components on the intermediary basis
    // y P_i(x), (1-y) P_i(x),  y (1-y) P_i(x)  (for x-dof)
    // x P_i(y), (1-x) P_i(y),  x (1-x) P_i(y)  (for y-dof)
    // x P_i(z), (1-x) P_i(z),  x (1-x) P_i(z)  (for z-dof)
    
    Real_wp vloc, threshold = 1e3*epsilon_machine;
    for (int k = 0; k <= rq; k++)
      for (int i = 0; i < r; i++)
        for (int j = 0; j <= r+1; j++)
          for (int l = 0; l <= r+1; l++)
            {
              int dof = NumDofsX(i, j, l);
              int node = 3*(N*k + j*r + i);
              vloc = ValBubble(l, k);
              if (abs(vloc) > threshold)
                Ch3.AddInteraction(dof, node, vloc);
              
              dof = NumDofsY(j, i, l);
              node = 3*(N*k + j*r + i) + 1;
              vloc = ValBubble(l, k);
              if (abs(vloc) > threshold)
                Ch3.AddInteraction(dof, node, vloc);
              
              dof = NumDofsZ(j, l, i);
              node = 3*(N*k + j*r + i) + 2;
              vloc = ValBubble(l, k);
              if (abs(vloc) > threshold)
                Ch3.AddInteraction(dof, node, vloc);              
            }
    
    // now we compute components on P_i(x), P_i(y) and P_i(z)
    Ch2.Reallocate(3*r*(r+2)*(rq+1), 3*r*(rq+1)*(rq+1));
    for (int k = 0; k <= rq; k++)
      for (int m = 0; m <= rq; m++)
        for (int i = 0; i < r; i++)
          {
            // y(1-y) P_j^{1,1}(y) P_i(x)
            for (int j = 0; j <= r+1; j++)
              {
                int node = 3*(N*k + j*r + i);
                int ny = 3*(r*(rq+1)*k + m*r + i);
                vloc = ValBubble(j, m);
                if (abs(vloc) > threshold)
                  Ch2.AddInteraction(node, ny, vloc);
              }
            
            // x(1-x) P_j^{1,1}(x) P_i(y)
            for (int j = 0; j <= r+1; j++)
              {
                int node = 3*(N*k + j*r + i) + 1;
                int ny = 3*(r*(rq+1)*k + m*r + i) + 1;
                vloc = ValBubble(j, m);
                if (abs(vloc) > threshold)
                  Ch2.AddInteraction(node, ny, vloc);
              }

            // x(1-x) P_j^{1,1}(x) P_i(z)
            for (int j = 0; j <= r+1; j++)
              {
                int node = 3*(N*k + j*r + i) + 2;
                int ny = 3*(r*(rq+1)*k + m*r + i) + 2;
                vloc = ValBubble(j, m);
                if (abs(vloc) > threshold)
                  Ch2.AddInteraction(node, ny, vloc);
              }           
          }
    
    // then value on quadrature points
    Ch1.Reallocate(3*r*(rq+1)*(rq+1), 3*(rq+1)*(rq+1)*(rq+1));
    for (int k = 0; k <= rq; k++)
      for (int m = 0; m <= rq; m++)
        for (int n = 0; n <= rq; n++)
          for (int i = 0; i < r; i++)
            {
              int ny = 3*(r*(rq+1)*k + m*r + i);
              int nx = 3*Num(n, m, k);
              vloc = ValLegendre(i, n);
              if (abs(vloc) > threshold)
                Ch1.AddInteraction(ny, nx, vloc);

              ny = 3*(r*(rq+1)*k + m*r + i) + 1;
              nx = 3*Num(m, n, k) + 1;
              vloc = ValLegendre(i, n);
              if (abs(vloc) > threshold)
                Ch1.AddInteraction(ny, nx, vloc);

              ny = 3*(r*(rq+1)*k + m*r + i) + 2;
              nx = 3*Num(m, k, n) + 2;
              vloc = ValLegendre(i, n);
              if (abs(vloc) > threshold)
                Ch1.AddInteraction(ny, nx, vloc);
            }
  }

  
  //! computing the projection from dofs on the face to other points of the face
  /*!
    \param[out] shX part of projection operator for interpolation along x
    \param[out] shY part of projection operator for interpolation along y
    \param[in] Num numbering of final points
    \param[in] points 1-D points, we compute projection on (points_i, points_j)
    \param[in] ValLegendre ?
    \param[in] ValBubble ?
    The operator shY^* shX^* does the projection
   */
  void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeSurfaceProjection(TinyVector<Matrix<Real_wp, General, RowSparse>, 6>& shX,
                             TinyVector<Matrix<Real_wp, General, RowSparse>, 6>& shY,
                             const VectReal_wp& points, const Matrix<Real_wp>& ValLegendre,
                             const Matrix<Real_wp>& ValBubble, const Matrix<int>& Num)
  {
    Matrix<Real_wp, General, ArrayRowSparse> Sh1, Sh2;
    
    // face 0-5 
    int r = order;
    int rq = points.GetM() - 1;
    Sh1.Reallocate(3*(r+2)*(r+2), 3*(rq+1)*(r+2));
    Sh2.Reallocate(3*(rq+1)*(r+2), 3*(rq+1)*(rq+1));
    Real_wp vloc, threshold = 1e3*epsilon_machine;
    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= r+1; j++)
            {
              int node = i*(r+2) + j;
              int n1 = i*(rq+1) + k;
              vloc = ValBubble(j, k);
              if (abs(vloc) > threshold)
                Sh1.AddInteraction(3*node, 3*n1, vloc);

              if (i < r)
                {
                  vloc = ValBubble(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(3*node+1, 3*n1+1, vloc);
                  
                }
              
              if (j < r)
                {
                  vloc = ValLegendre(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(3*node+2, 3*n1+2, vloc);
                }
            }
      }

    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= rq; j++)
            {
              int node = i*(rq+1) + j;
              int n2 = Num(k, j);
              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(3*node, 3*n2, vloc);
              
              if (i < r)
                {
                  vloc = ValLegendre(i, k);
                  if (abs(vloc) > threshold)
                    Sh2.AddInteraction(3*node+1, 3*n2+1, vloc);
                }
              
              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(3*node+2, 3*n2+2, vloc);
            }
      }
    
    Copy(Sh1, shX(0));
    Copy(Sh2, shY(0));
    Copy(Sh1, shX(5));
    Copy(Sh2, shY(5));
    
    Sh1.Clear(); Sh2.Clear();
    Sh1.Reallocate(3*(r+2)*(r+2), 3*(rq+1)*(r+2));
    Sh2.Reallocate(3*(rq+1)*(r+2), 3*(rq+1)*(rq+1));
    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= r+1; j++)
            {
              int node = i*(r+2) + j;
              int n1 = i*(rq+1) + k;
              vloc = ValBubble(j, k);
              if (abs(vloc) > threshold)
                Sh1.AddInteraction(3*node+1, 3*n1+1, vloc);

              if (i < r)
                {
                  vloc = ValBubble(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(3*node, 3*n1, vloc);
                }
              
              if (j < r)
                {
                  vloc = ValLegendre(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(3*node+2, 3*n1+2, vloc);
                }
            }
      }

    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= rq; j++)
            {
              int node = i*(rq+1) + j;
              int n2 = Num(k, j);
              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(3*node+1, 3*n2+1, vloc);
              
              if (i < r)
                {
                  vloc = ValLegendre(i, k);
                  if (abs(vloc) > threshold)
                    Sh2.AddInteraction(3*node, 3*n2, vloc);
                }
              
              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(3*node+2, 3*n2+2, vloc);
            }
      }
    
    Copy(Sh1, shX(1));
    Copy(Sh2, shY(1));
    Copy(Sh1, shX(4));
    Copy(Sh2, shY(4));
    
    Sh1.Clear(); Sh2.Clear();
    Sh1.Reallocate(3*(r+2)*(r+2), 3*(rq+1)*(r+2));
    Sh2.Reallocate(3*(rq+1)*(r+2), 3*(rq+1)*(rq+1));
    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= r+1; j++)
            {
              int node = i*(r+2) + j;
              int n1 = i*(rq+1) + k;
              vloc = ValBubble(j, k);
              if (abs(vloc) > threshold)
                Sh1.AddInteraction(3*node+2, 3*n1+2, vloc);

              if (i < r)
                {
                  vloc = ValBubble(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(3*node, 3*n1, vloc);
                }
              
              if (j < r)
                {
                  vloc = ValLegendre(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(3*node+1, 3*n1+1, vloc);
                }              
            }
      }

    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= rq; j++)
            {
              int node = i*(rq+1) + j;
              int n2 = Num(k, j);
              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(3*node+2, 3*n2+2, vloc);
              
              if (i < r)
                {
                  vloc = ValLegendre(i, k);
                  if (abs(vloc) > threshold)
                    Sh2.AddInteraction(3*node, 3*n2, vloc);
                }
              
              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(3*node+1, 3*n2+1, vloc);
            }
      }
    
    Copy(Sh1, shX(2));
    Copy(Sh2, shY(2));
    Copy(Sh1, shX(3));
    Copy(Sh2, shY(3));
  }
  
  
  //! computing the projection from nodal points on the face to other points of the face
  //! (curl is computed)
  /*!
    \param[out] shX part of projection operator for interpolation along x
    \param[out] shY part of projection operator for interpolation along y
    \param[in] Num numbering of final points
    \param[in] points 1-D points, we compute projection on (points_i, points_j)
    \param[in] ValLegendre ?
    \param[in] ValBubble ?
    \param[in] DerBubble ?
    The operator shY^* shX^* does the projection
   */
  void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeSurfaceCurlProjection(TinyVector<Matrix<Real_wp, General, RowSparse>, 6>& shX,
                                 TinyVector<Matrix<Real_wp, General, RowSparse>, 6>& shY,
                                 const VectReal_wp& points,
                                 const Matrix<Real_wp>& DerBubble,
                                 const Matrix<Real_wp>& ValLegendre,
                                 const Matrix<Real_wp>& ValBubble, const Matrix<int>& Num)
  {
    Matrix<Real_wp, General, ArrayRowSparse> Sh1, Sh2;
    
    // face 0-5 
    int r = order;
    int rq = points.GetM() - 1;
    Sh1.Reallocate(5*(r+2)*(r+2), 6*(rq+1)*(r+2));
    Sh2.Reallocate(6*(rq+1)*(r+2), 3*(rq+1)*(rq+1));
    Real_wp vloc, threshold = 1e3*epsilon_machine;
    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= r+1; j++)
            {
              int node = i*(r+2) + j;
              int n1 = i*(rq+1) + k;
              vloc = DerBubble(j, k);
              if (abs(vloc) > threshold)
                Sh1.AddInteraction(5*node, 6*n1+1, vloc);

              vloc = ValBubble(j, k);
              if (abs(vloc) > threshold)
                Sh1.AddInteraction(5*node, 6*n1, -vloc);
              
              if (i < r)
                {
                  vloc = DerBubble(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node+1, 6*n1+3, -vloc);
                  
                  vloc = ValBubble(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node+3, 6*n1+2, vloc);
                }
              
              if (j < r)
                {
                  vloc = ValLegendre(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node+2, 6*n1+5, vloc);

                  vloc = ValLegendre(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node+4, 6*n1+4, -vloc);
                }
            }
      }
    
    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= rq; j++)
            {
              int node = i*(rq+1) + j;
              int n2 = Num(k, j);
              
              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+1, 3*n2+1, vloc);

              vloc = DerBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node, 3*n2+2, vloc);
              
              if (i < r)
                {
                  vloc = ValLegendre(i, k);
                  if (abs(vloc) > threshold)
                    Sh2.AddInteraction(6*node+3, 3*n2, vloc);
                  
                  vloc = ValLegendre(i, k);
                  if (abs(vloc) > threshold)
                    Sh2.AddInteraction(6*node+2, 3*n2+2, vloc);
                }
              
              vloc = DerBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+5, 3*n2, vloc);

              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+4, 3*n2+1, vloc);
            }
      }

    Copy(Sh1, shX(0));
    Copy(Sh2, shY(0));
    Copy(Sh1, shX(5));
    Copy(Sh2, shY(5));
    
    Sh1.Clear(); Sh2.Clear();
    Sh1.Reallocate(5*(r+2)*(r+2), 6*(rq+1)*(r+2));
    Sh2.Reallocate(6*(rq+1)*(r+2), 3*(rq+1)*(rq+1));
    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= r+1; j++)
            {
              int node = i*(r+2) + j;
              int n1 = i*(rq+1) + k;
              vloc = DerBubble(j, k);
              if (abs(vloc) > threshold)
                Sh1.AddInteraction(5*node+1, 6*n1+3, -vloc);

              vloc = ValBubble(j, k);
              if (abs(vloc) > threshold)
                Sh1.AddInteraction(5*node+1, 6*n1+2, vloc);
              
              if (i < r)
                {
                  vloc = ValBubble(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node+3, 6*n1, -vloc);
                  
                  vloc = DerBubble(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node, 6*n1+1, vloc);
                }
              
              if (j < r)
                {
                  vloc = ValLegendre(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node+2, 6*n1+4, -vloc);

                  vloc = ValLegendre(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node+4, 6*n1+5, vloc);
                }
            }
      }
    
    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= rq; j++)
            {
              int node = i*(rq+1) + j;
              int n2 = Num(k, j);
              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+3, 3*n2, vloc);
              
              vloc = DerBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+2, 3*n2+2, vloc);
              
              if (i < r)
                {
                  vloc = ValLegendre(i, k);
                  if (abs(vloc) > threshold)
                    Sh2.AddInteraction(6*node, 3*n2+2, vloc);
                  
                  vloc = ValLegendre(i, k);
                  if (abs(vloc) > threshold)
                    Sh2.AddInteraction(6*node+1, 3*n2+1, vloc);
                }
              
              vloc = DerBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+4, 3*n2+1, vloc);

              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+5, 3*n2, vloc);
            }
      }
        
    Copy(Sh1, shX(1));
    Copy(Sh2, shY(1));
    Copy(Sh1, shX(4));
    Copy(Sh2, shY(4));
    
    Sh1.Clear(); Sh2.Clear();
    Sh1.Reallocate(5*(r+2)*(r+2), 6*(rq+1)*(r+2));
    Sh2.Reallocate(6*(rq+1)*(r+2), 3*(rq+1)*(rq+1));
    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= r+1; j++)
            {
              int node = i*(r+2) + j;
              int n1 = i*(rq+1) + k;
              
              vloc = DerBubble(j, k);
              if (abs(vloc) > threshold)
                Sh1.AddInteraction(5*node+2, 6*n1+5, vloc);
              
              vloc = ValBubble(j, k);
              if (abs(vloc) > threshold)
                Sh1.AddInteraction(5*node+2, 6*n1+4, -vloc);
              
              if (i < r)
                {
                  vloc = ValBubble(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node+3, 6*n1+1, vloc);
                  
                  vloc = DerBubble(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node, 6*n1, -vloc);
                }
              
              if (j < r)
                {
                  vloc = ValLegendre(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node+1, 6*n1+2, vloc);

                  vloc = ValLegendre(j, k);
                  if (abs(vloc) > threshold)
                    Sh1.AddInteraction(5*node+4, 6*n1+3, -vloc);
                }
            }
      }

    for (int k = 0; k <= rq; k++)
      {
        for (int i = 0; i <= r+1; i++)
          for (int j = 0; j <= rq; j++)
            {
              int node = i*(rq+1) + j;
              int n2 = Num(k, j);
              
              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+5, 3*n2, vloc);
              
              vloc = DerBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+4, 3*n2+1, vloc);
              
              if (i < r)
                {
                  vloc = ValLegendre(i, k);
                  if (abs(vloc) > threshold)
                    Sh2.AddInteraction(6*node+1, 3*n2+1, vloc);
                  
                  vloc = ValLegendre(i, k);
                  if (abs(vloc) > threshold)
                    Sh2.AddInteraction(6*node, 3*n2+2, vloc);
                }
              
              vloc = DerBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+2, 3*n2+2, vloc);

              vloc = ValBubble(i, k);
              if (abs(vloc) > threshold)
                Sh2.AddInteraction(6*node+3, 3*n2, vloc);
            }
      }
    
    Copy(Sh1, shX(2));
    Copy(Sh2, shY(2));
    Copy(Sh1, shX(3));
    Copy(Sh2, shY(3));
  }

  
  //! constructing stiffness matrix
  void HexahedronHcurlOptimalHpFirstFamily::ConstructStiffnessMatrix()
  {
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    const Matrix<Real_wp>& dG = lob_quad.GradPhi();
    
    Matrix<Real_wp, General, ArrayRowSparse>
      Rh(3*nb_points_quadrature_inside, 3*nb_points_quadrature_inside);
    
    for (int i = 0; i <= order_quad; i++)
      for (int j = 0; j <= order_quad; j++)
	for (int k = 0; k <= order_quad; k++)
	  for (int m = 0; m <= order_quad; m++)
	    {
	      int col = 3*NumQuad3D(i,j,k);
	      int row = 3*NumQuad3D(i,j,m)+1;
	      Rh.AddInteraction(row, col, -dG(m,k));
              
	      row = 3*NumQuad3D(i,m,k)+2;
	      Rh.AddInteraction(row, col, dG(m,j));
	      	      
	      col = 3*NumQuad3D(i,j,k)+1;
	      row = 3*NumQuad3D(i,j,m);
	      Rh.AddInteraction(row, col, dG(m,k));
	      	      
	      row = 3*NumQuad3D(m,j,k)+2;
	      Rh.AddInteraction(row, col, -dG(m,i));
	      	      
	      col = 3*NumQuad3D(i,j,k)+2;
	      row = 3*NumQuad3D(i,m,k);
	      Rh.AddInteraction(row, col, -dG(m,j));
	      	      
	      row = 3*NumQuad3D(m,j,k)+1;
	      Rh.AddInteraction(row, col, dG(m,i));
            }
    
    Seldon::Copy(Rh, rh_loc);
    
    // computation of constant mass and stiffness matrix
    Matrix<Real_wp> IntPiPj(order, order), IntPiLj(order, order+2), IntLiLj(order+2, order+2);
    Matrix<Real_wp> IntPiDLj(order, order+2), IntLiDLj(order+2, order+2),
      IntDLiDLj(order+2, order+2);
    
    IntPiPj.Fill(0); IntPiLj.Fill(0); IntLiLj.Fill(0);
    IntPiDLj.Fill(0); IntLiDLj.Fill(0); IntDLiDLj.Fill(0);
    const VectReal_wp& weights1d = this->Weights1D();
    for (int k = 0; k < weights1d.GetM(); k++)
      {
        for (int i = 0; i < order; i++)
          {
            for (int j = 0; j < order; j++)
              IntPiPj(i, j) += weights1d(k)*ValLeg(i, k)*ValLeg(j, k);
            
            for (int j = 0; j <= order+1; j++)
              {
                IntPiLj(i, j) += weights1d(k)*ValLeg(i, k)*ValBubbleJac(j, k);
                IntPiDLj(i, j) += weights1d(k)*ValLeg(i, k)*DerBubbleJac(j, k);
              }
          }
        
        for (int i = 0; i <= order+1; i++)
          for (int j = 0; j <= order+1; j++)
            {
              IntLiLj(i, j) += weights1d(k)*ValBubbleJac(i, k)*ValBubbleJac(j, k);
              IntLiDLj(i, j) += weights1d(k)*ValBubbleJac(i, k)*DerBubbleJac(j, k);
              IntDLiDLj(i, j) += weights1d(k)*DerBubbleJac(i, k)*DerBubbleJac(j, k);
            }
      }
    
    // sparse mass matrix
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++)
        sparse_mass_matrix(i, j).Reallocate(nb_dof_loc, nb_dof_loc);

    mass_matrix_chol.Reallocate(nb_dof_loc, nb_dof_loc);
    mass_matrix_chol.Fill(0);
        
    Real_wp vloc;
    VectReal_wp value(nb_dof_loc); IVect col(nb_dof_loc);
    Real_wp threshold = 1e3*epsilon_machine;
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
        for (int k = 0; k <= order+1; k++)
          {
            int irow = NumDofsX(i, j, k);
            int nb = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsX(i2, j2, k2);
                    vloc = IntPiPj(i, i2)*IntLiLj(j, j2)*IntLiLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        if (irow <= jcol)
                          mass_matrix_chol(irow, jcol) = vloc;
                        
                        nb++;
                      }
                  }
            
            sparse_mass_matrix(0, 0).AddInteractionRow(irow, nb, col, value);
            
            nb = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsY(j2, i2, k2);
                    vloc = IntPiLj(i, j2)*IntPiLj(i2, j)*IntLiLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        nb++;
                      }
                  }
            
            sparse_mass_matrix(0, 1).AddInteractionRow(irow, nb, col, value);

            nb = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsZ(j2, k2, i2);
                    vloc = IntPiLj(i, j2)*IntPiLj(i2, k)*IntLiLj(j, k2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        nb++;
                      }
                  }
            
            sparse_mass_matrix(0, 2).AddInteractionRow(irow, nb, col, value);
            
            irow = NumDofsY(j, i, k);
            nb = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsY(j2, i2, k2);
                    vloc = IntPiPj(i, i2)*IntLiLj(j, j2)*IntLiLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        if (irow <= jcol)
                          mass_matrix_chol(irow, jcol) = vloc;
                        
                        nb++;
                      }
                  }
            
            sparse_mass_matrix(1, 1).AddInteractionRow(irow, nb, col, value);

            nb = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsZ(j2, k2, i2);
                    vloc = IntPiLj(i, k2)*IntPiLj(i2, k)*IntLiLj(j, j2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;                        
                        nb++;
                      }
                  }
            
            sparse_mass_matrix(1, 2).AddInteractionRow(irow, nb, col, value);
            
            irow = NumDofsZ(j, k, i);
            nb = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsZ(j2, k2, i2);
                    vloc = IntPiPj(i, i2)*IntLiLj(j, j2)*IntLiLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        if (irow <= jcol)
                          mass_matrix_chol(irow, jcol) = vloc;
                        
                        nb++;
                      }
                  }
            
            sparse_mass_matrix(2, 2).AddInteractionRow(irow, nb, col, value);
          }
    
    Transpose(sparse_mass_matrix(0, 1), sparse_mass_matrix(1, 0));
    Transpose(sparse_mass_matrix(0, 2), sparse_mass_matrix(2, 0));
    Transpose(sparse_mass_matrix(1, 2), sparse_mass_matrix(2, 1));
    
    GetCholesky(mass_matrix_chol);
    ConvertToSparse(mass_matrix_chol, sparse_matrix_chol, threshold);
    mass_matrix_chol.Clear();

    // sparse stiffness matrix
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++)
        sparse_stiff_matrix(i, j).Reallocate(nb_dof_loc, nb_dof_loc);
    
    IVect col2(nb_dof_loc);
    VectReal_wp value2(nb_dof_loc);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
        for (int k = 0; k <= order+1; k++)
          {
            int irow = NumDofsY(j, i, k);
            int irow2 = NumDofsZ(j, k, i);
            int nb = 0, nb2 = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsY(j2, i2, k2);
                    int jcol2 = NumDofsZ(j2, k2, i2);
                    vloc = IntPiPj(i, i2)*IntLiLj(j, j2)*IntDLiDLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = -IntPiDLj(i, k2)*IntLiLj(j, j2)*IntPiDLj(i2, k);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol2;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = -IntLiLj(j, j2)*IntPiDLj(i2, k)*IntPiDLj(i, k2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                    
                    vloc = IntLiLj(j, j2)*IntDLiDLj(k, k2)*IntPiPj(i, i2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol2;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                  }
            
            sparse_stiff_matrix(0, 0).AddInteractionRow(irow, nb, col, value);
            sparse_stiff_matrix(0, 0).AddInteractionRow(irow2, nb2, col2, value2);
            
            nb = 0; nb2 = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsX(i2, j2, k2);
                    int jcol2 = NumDofsZ(j2, k2, i2);
                    vloc = -IntPiLj(i2, j)*IntPiLj(i, j2)*IntDLiDLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = IntLiDLj(j, j2)*IntPiLj(i, k2)*IntPiDLj(i2, k);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol2;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = IntPiLj(i2, j)*IntLiDLj(j2, k)*IntPiDLj(i, k2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                    
                    vloc = -IntLiDLj(j, j2)*IntLiDLj(k2, k)*IntPiPj(i, i2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol2;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                  }
            
            sparse_stiff_matrix(0, 1).AddInteractionRow(irow, nb, col, value);
            sparse_stiff_matrix(0, 1).AddInteractionRow(irow2, nb2, col2, value2);
            
            nb = 0; nb2 = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsX(i2, j2, k2);
                    int jcol2 = NumDofsY(j2, i2, k2);
                    vloc = IntPiLj(i2, j)*IntPiDLj(i, j2)*IntLiDLj(k2, k);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = -IntLiDLj(j, j2)*IntPiPj(i, i2)*IntLiDLj(k2, k);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol2;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = -IntPiLj(i2, j)*IntDLiDLj(k, j2)*IntPiLj(i, k2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                    
                    vloc = IntLiDLj(j, j2)*IntPiDLj(i2, k)*IntPiLj(i, k2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol2;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                  }
            
            sparse_stiff_matrix(0, 2).AddInteractionRow(irow, nb, col, value);
            sparse_stiff_matrix(0, 2).AddInteractionRow(irow2, nb2, col2, value2);
            
            irow = NumDofsX(i, j, k);
            irow2 = NumDofsZ(j, k, i);
            nb = 0; nb2 = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsX(i2, j2, k2);
                    int jcol2 = NumDofsZ(j2, k2, i2);
                    vloc = IntPiPj(i, i2)*IntLiLj(j, j2)*IntDLiDLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = -IntPiDLj(i, j2)*IntLiLj(j, k2)*IntPiDLj(i2, k);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol2;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = -IntPiDLj(i2, j)*IntLiLj(j2, k)*IntPiDLj(i, k2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                    
                    vloc = IntDLiDLj(j, j2)*IntLiLj(k, k2)*IntPiPj(i, i2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol2;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                  }
            
            sparse_stiff_matrix(1, 1).AddInteractionRow(irow, nb, col, value);
            sparse_stiff_matrix(1, 1).AddInteractionRow(irow2, nb2, col2, value2);
            
            nb = 0; nb2 = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsX(i2, j2, k2);
                    int jcol2 = NumDofsY(j2, i2, k2);
                    vloc = -IntPiPj(i, i2)*IntLiDLj(j, j2)*IntLiDLj(k2, k);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = IntPiDLj(i, j2)*IntPiLj(i2, j)*IntLiDLj(k2, k);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol2;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = IntPiDLj(i2, j)*IntLiDLj(k, j2)*IntPiLj(i, k2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                    
                    vloc = -IntDLiDLj(j, j2)*IntPiLj(i2, k)*IntPiLj(i, k2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol2;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                  }
            
            sparse_stiff_matrix(1, 2).AddInteractionRow(irow, nb, col, value);
            sparse_stiff_matrix(1, 2).AddInteractionRow(irow2, nb2, col2, value2);
            
            irow = NumDofsX(i, j, k);
            irow2 = NumDofsY(j, i, k);
            nb = 0; nb2 = 0;
            for (int i2 = 0; i2 < order; i2++)
              for (int j2 = 0; j2 <= order+1; j2++)
                for (int k2 = 0; k2 <= order+1; k2++)
                  {
                    int jcol = NumDofsX(i2, j2, k2);
                    int jcol2 = NumDofsY(j2, i2, k2);
                    vloc = IntPiPj(i, i2)*IntDLiDLj(j, j2)*IntLiLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = -IntPiDLj(i, j2)*IntPiDLj(i2, j)*IntLiLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col(nb) = jcol2;
                        value(nb) = vloc;
                        nb++;
                      }
                    
                    vloc = -IntPiDLj(i2, j)*IntPiDLj(i, j2)*IntLiLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                    
                    vloc = IntDLiDLj(j, j2)*IntPiPj(i, i2)*IntLiLj(k, k2);
                    if (abs(vloc) > threshold)
                      {
                        col2(nb2) = jcol2;
                        value2(nb2) = vloc;
                        nb2++;
                      }
                  }
            
            sparse_stiff_matrix(2, 2).AddInteractionRow(irow, nb, col, value);
            sparse_stiff_matrix(2, 2).AddInteractionRow(irow2, nb2, col2, value2);
          }    
    
    Transpose(sparse_stiff_matrix(0, 1), sparse_stiff_matrix(1, 0));
    Transpose(sparse_stiff_matrix(0, 2), sparse_stiff_matrix(2, 0));
    Transpose(sparse_stiff_matrix(1, 2), sparse_stiff_matrix(2, 1));
  }

  
  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are located
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    //HexahedronReference<2>::ComputeProjectionDofRef(feval, contrib);
    
    Vector1 feval_weight = feval;
    const VectReal_wp& weights3d = this->WeightsND();
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

  
  void HexahedronHcurlOptimalHpFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {
    for (int i = 0; i < this->order; i++)
      if (i%2 == 0)
	{
	  contrib(2*this->order+i) = -contrib(2*this->order+i);
 	  contrib(3*this->order+i) = -contrib(3*this->order+i);
	}
  }


  void HexahedronHcurlOptimalHpFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    for (int i = 0; i < this->order; i++)
      if (i%2 == 0)
	{
	  contrib(2*this->order+i) = -contrib(2*this->order+i);
 	  contrib(3*this->order+i) = -contrib(3*this->order+i);
	}
  }

  
  //! projection from values on quadrature points to dof components
  /*!
    \param[in] Equad values on quadrature points
    \param[out] Edof dof components
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily::
  ProjectQuadratureToDofGen(const Vector1& Equad, Vector2& Edof) const
  {
    Vector1 feval_weight = Equad;
    const VectReal_wp& weights3d = this->WeightsND();
    for (int i = 0; i < feval_weight.GetM(); i+=3)
      {
       int j = i/3;
       feval_weight(i) *= weights3d(j);
       feval_weight(i+1) *= weights3d(j);
       feval_weight(i+2) *= weights3d(j);
     }
    
    Edof.Reallocate(nb_dof_loc);
    ApplyCh(feval_weight, Edof);
    SolveMassMatrix(Edof);
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
  void HexahedronHcurlOptimalHpFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    Vector1 fx(3*(order+1)*(order+2)), fy(3*(order+2)*(order+2));
    FillZero(fx); FillZero(fy);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
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
  

  //! computation of U on nodal points (on reference element)
  /*!
    \param[in] u_loc components of U on dofs
    \param[out] Uloc_node values of U on nodal points
  */
  template<class Vector1,class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily::
  ComputeNodalValuesGen(const Vector1& u_loc, Vector2& u_node) const
  {
    Vector2 Uy(ch3_node.GetN()), Ux(ch2_node.GetN());
    Mlt(SeldonTrans, ch3_node, u_loc, Uy);
    Mlt(SeldonTrans, ch2_node, Uy, Ux);
    Mlt(SeldonTrans, ch1_node, Ux, u_node);
  }
  
  
  //! computation of u on nodal points of a face
  /*!
    \param[in] Uh dof components of u
    \param[out] Vh values of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1,class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily::ComputeValueBoundaryGen(const Vector1& Uh,
								    Vector2& Vh, int num_loc) const
  {
    Vector2 Ux(sh_loc(0).GetN()), Uy(sh_nodeX(0).GetN());
    
    Seldon::Mlt(SeldonTrans, sh_loc(num_loc), Uh, Ux);
    Seldon::Mlt(SeldonTrans, sh_nodeX(num_loc), Ux, Uy);
    Seldon::Mlt(SeldonTrans, sh_nodeY(num_loc), Uy, Vh);

  }
  
  
  //! computation of curl(u) on nodal points of a face
  /*!
    \param[in] Uh dof components of u
    \param[out] Vh curl of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1,class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily::ComputeCurlBoundaryGen(const Vector1& Uh,
								   Vector2& Vh, int num_loc) const
  {
    Vector2 Ux(sh_curl_loc(0).GetN()), Uy(sh_curl_nodeX(0).GetN());
    
    Seldon::Mlt(SeldonTrans, sh_curl_loc(num_loc), Uh, Ux);
    Seldon::Mlt(SeldonTrans, sh_curl_nodeX(num_loc), Ux, Uy);
    Seldon::Mlt(SeldonTrans, sh_curl_nodeY(num_loc), Uy, Vh);
  }
  
    
  //! x is overwritten by M^-1 x where M is the mass matrix    
  template<class Vector1>
  void HexahedronHcurlOptimalHpFirstFamily::SolveMassMatrixGen(Vector1& x) const
  {
    Seldon::SolveCholesky(SeldonNoTrans, sparse_matrix_chol, x);
    Seldon::SolveCholesky(SeldonTrans, sparse_matrix_chol, x);
  }
  

  //! x is overwritten by L^-1 x or L^-T where M = L L^t is the mass matrix    
  template<class Vector1>
  void HexahedronHcurlOptimalHpFirstFamily::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& x) const
  {
    Seldon::SolveCholesky(TransA, sparse_matrix_chol, x);
  }
    
  
  //! x is overwritten by M x where M is the mass matrix  
  template<class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily::MltMassMatrixGen(Vector2& x) const
  {
    Seldon::MltCholesky(SeldonTrans, sparse_matrix_chol, x);
    Seldon::MltCholesky(SeldonNoTrans, sparse_matrix_chol, x);
  }
    
  
  //! Integration against curl of basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \nabla \times \varphi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Rh Uh
    where (Rh)_{i,j} = \nabla \times \varphi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Uquad(3*nb_points_quadrature_inside);
    ApplyRhQuadrature(Uh, Uquad);
    ApplyCh(Uquad, Vh);
  }
    
  
  //! computation of curl of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh curl of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j curl phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = curl phi_i(xi_j)
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Uquad(3*nb_points_quadrature_inside);
    ApplyChTranspose(Uh, Uquad);
    ApplyRhQuadratureTranspose(Uquad, Vh);
  }
  

  //! Integration against curl of basis functions associated with quadrature points
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vh Uh_i = \int_K f curl(\psi_i) dx
    This operation is equivalent to a matrix vector product
    Vh = Rh Uh
    where (Rh)_{i,j} = \nabla \times \psi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily::ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(rh_loc, Uh, Vh);
  }
  
    
  //! Computation of curl on quadrature points from values on quadrature points
  /*!
    \param[in] Uh values of u on quadrature points
    \param[out] Vh curl of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Rh* Uh
    where (Rh)_{i,j} = \nabla \times \psi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  }
  
    
  //! Integration against basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Ch Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Ux(ch1_loc.GetM()), Uy(ch2_loc.GetM());
    Seldon::Mlt(ch1_loc, Uh, Ux);
    Seldon::Mlt(ch2_loc, Ux, Uy);
    Seldon::Mlt(ch3_loc, Uy, Vh);
  }
  
    
  //! computation of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Ch* Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Uy(ch3_loc.GetN()), Ux(ch2_loc.GetN());
    Seldon::Mlt(SeldonTrans, ch3_loc, Uh, Uy);
    Seldon::Mlt(SeldonTrans, ch2_loc, Uy, Ux);
    Seldon::Mlt(SeldonTrans, ch1_loc, Ux, Vh);
  }
  

  //! computation of u on quadrature points of a boundary
  /*!
    \param[in] num_loc boundary number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vector2 Ux(sh_loc(0).GetN()), Uy(sh_locX(0).GetN());
    
    Seldon::Mlt(SeldonTrans, sh_loc(num_loc), Uh, Ux);
    Seldon::Mlt(SeldonTrans, sh_locX(num_loc), Ux, Uy);
    Seldon::Mlt(SeldonTrans, sh_locY(num_loc), Uy, Vh);
  }
  
  
  //! integration against basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    T0 one; SetComplexOne(one);
    Vector2 Ux(sh_locY(0).GetM()), Uy(sh_locX(0).GetM());
    
    Seldon::Mlt(sh_locY(num_loc), Uh, Ux);
    Seldon::Mlt(sh_locX(num_loc), Ux, Uy);
    Seldon::MltAdd(alpha, sh_loc(num_loc), Uy, one, Vh);
  }
  
  
  //! computation of curl of u on quadrature points of a boundary
  /*!
    \param[in] num_loc boundary number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh curl of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = curl \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vector2 Ux(sh_curl_loc(0).GetN()), Uy(sh_curl_locX(0).GetN());
    
    Seldon::Mlt(SeldonTrans, sh_curl_loc(num_loc), Uh, Ux);
    Seldon::Mlt(SeldonTrans, sh_curl_locX(num_loc), Ux, Uy);
    Seldon::Mlt(SeldonTrans, sh_curl_locY(num_loc), Uy, Vh);
  }
  
  
  //! integration against curl of basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f curl(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
  */
  template<class T0, class Vector1, class Vector2>
  void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    T0 one; SetComplexOne(one);
    Vector2 Ux(sh_curl_locY(0).GetM()), Uy(sh_curl_locX(0).GetM());
    
    Seldon::Mlt(sh_curl_locY(num_loc), Uh, Ux);
    Seldon::Mlt(sh_curl_locX(num_loc), Ux, Uy);
    Seldon::MltAdd(alpha, sh_curl_loc(num_loc), Uy, one, Vh);
  }
  
  
  //! we add constant mass matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass tensor
    \param[inout] A matrix modified
    A(m:, n:) = A(m:, n:) + M
    where M is the mass matrix, M_ij = \int mass \varphi_j \varphi_i dx
   */
  template<class T, class Prop>
  void HexahedronHcurlOptimalHpFirstFamily
  ::AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& mass, VirtualMatrix<T>& A) const
  {
    for (int i = 0; i < nb_dof_loc; i++)
      {
        for (int p = 0; p < 3; p++)
          for (int q = 0; q < 3; q++)
            for (int jloc = 0; jloc < sparse_mass_matrix(p, q).GetRowSize(i); jloc++)
              {
                int j = sparse_mass_matrix(p, q).Index(i, jloc);
		A.AddInteraction(m+i, n+j, sparse_mass_matrix(p, q).Value(i, jloc)*mass(p, q));
              }
      }
  }
  
    
  //! we add constant stiffness matrix \int_K  C curl(phi_j) curl(phi_i)
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] C tensor
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + \sum_p  C(p, q) S^{p,q}
    where S^{p,q} is a stiffness matrix equal to :
    (S^{p,q })_{i, j} = \int_K curl(\phi_j)_q curl(\phi_i)_p  dx
   */  
  template<class T, class Prop>
  void HexahedronHcurlOptimalHpFirstFamily
  ::AddConstantStiffnessMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const
  {
    for (int i = 0; i < nb_dof_loc; i++)
      {
        for (int p = 0; p < 3; p++)
          for (int q = 0; q < 3; q++)
            for (int jloc = 0; jloc < sparse_stiff_matrix(p, q).GetRowSize(i); jloc++)
              {
                int j = sparse_stiff_matrix(p, q).Index(i, jloc);
		A.AddInteraction(m+i, n+j, sparse_stiff_matrix(p, q).Value(i, jloc)*C(p, q));
              }
      }
  }
  

  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] phi values of basis functions on point
  */
  void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeValuesPhiRef(const R3& point_loc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    int node = 0;

    VectReal_wp Px, Py, Pz, Jx, Jy, Jz, Jx2, Jy2, Jz2;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*x-1, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*y-1, Py);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*z-1, Pz);

    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*x-1, Jx);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*y-1, Jy);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*z-1, Jz);
    
    Jx2.Reallocate(order+2);
    Jy2.Reallocate(order+2);
    Jz2.Reallocate(order+2);
    
    Jx2(0) = 1.0-x; Jx2(order+1) = x;
    Jy2(0) = 1.0-y; Jy2(order+1) = y;
    Jz2(0) = 1.0-z; Jz2(order+1) = z;
    for (int i = 0; i < order; i++)
      {
        Jx2(i+1) = x*(1.0-x)*Jx(i);
        Jy2(i+1) = y*(1.0-y)*Jy(i);
        Jz2(i+1) = z*(1.0-z)*Jz(i);
      }
    
    for (int k = 0; k <= order-1; k++)
      for (int i = 0; i <= order+1; i++)
        for (int j = 0; j <= order+1; j++)
          {  
            node = NumDofsX(k, i, j);
	    phi(node)(0) = Px(k)*Jy2(i)*Jz2(j);
            phi(node)(1) = 0;
            phi(node)(2) = 0;
            
	    node = NumDofsY(i, k, j);
            phi(node)(0) = 0;
            phi(node)(1) = Py(k)*Jx2(i)*Jz2(j);
            phi(node)(2) = 0;
	    
            node = NumDofsZ(i, j, k);
            phi(node)(0) = 0;
            phi(node)(1) = 0;
            phi(node)(2) = Jx2(i)*Jy2(j)*Pz(k);
	  }	
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point
  */
  void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeCurlPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    int node = 0;
    
    VectReal_wp Px, Py, Pz, Jx, Jy, Jz, dJx, dJy, dJz;
    VectReal_wp Jx2, Jy2, Jz2, dJx2, dJy2, dJz2;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*x-1, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*y-1, Py);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*z-1, Pz);
    
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*x-1, Jx, dJx);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*y-1, Jy, dJy);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*z-1, Jz, dJz);
    
    Jx2.Reallocate(order+2); dJx2.Reallocate(order+2);
    Jy2.Reallocate(order+2); dJy2.Reallocate(order+2);
    Jz2.Reallocate(order+2); dJz2.Reallocate(order+2);
    
    Jx2(0) = 1.0-x; Jx2(order+1) = x; dJx2(0) = -1.0; dJx2(order+1) = 1.0;
    Jy2(0) = 1.0-y; Jy2(order+1) = y; dJy2(0) = -1.0; dJy2(order+1) = 1.0;
    Jz2(0) = 1.0-z; Jz2(order+1) = z; dJz2(0) = -1.0; dJz2(order+1) = 1.0;
    for (int i = 0; i < order; i++)
      {
        Jx2(i+1) = x*(1.0-x)*Jx(i);
        Jy2(i+1) = y*(1.0-y)*Jy(i);
        Jz2(i+1) = z*(1.0-z)*Jz(i);

        dJx2(i+1) = (1.0 - 2.0*x)*Jx(i) + 2.0*x*(1.0-x)*dJx(i);
        dJy2(i+1) = (1.0 - 2.0*y)*Jy(i) + 2.0*y*(1.0-y)*dJy(i);
        dJz2(i+1) = (1.0 - 2.0*z)*Jz(i) + 2.0*z*(1.0-z)*dJz(i);
      }
    
    for (int k = 0; k <= order-1; k++)
      for (int i = 0; i <= order+1; i++)
        for (int j = 0; j <= order+1; j++)
          { 
	    node = NumDofsX(k, i, j); 
            res(node)(0) = 0;
            res(node)(1) = Px(k)*dJz2(j)*Jy2(i);
            res(node)(2) = -Px(k)*dJy2(i)*Jz2(j);
	    
	    node = NumDofsY(i, k, j);
            res(node)(0) = -Py(k)*Jx2(i)*dJz2(j);
            res(node)(1) = 0;
            res(node)(2) = Py(k)*dJx2(i)*Jz2(j);
	    
	    node = NumDofsZ(i, j, k);
            res(node)(0) = Pz(k)*Jx2(i)*dJy2(j);
            res(node)(1) = -Pz(k)*dJx2(i)*Jy2(j);
            res(node)(2) = 0;
	  }      
  }
  
       
  //! computation of val_phi \f$ = \varphi_{node} (\xi_{kg})  \f$
  /*!
    \param[in] kg quadrature point number
    \param[out] phi values of basis functions
    (xi_i) are Gauss-Legendre quadrature points
  */  
  void HexahedronHcurlOptimalHpFirstFamily::GetValuePhiOnQuadraturePoint(int kg, VectR3& phi) const
  {
    int kx = CoordinateQuad3D(kg, 0);
    int ky = CoordinateQuad3D(kg, 1);
    int kz = CoordinateQuad3D(kg, 2);
    
    Real_wp zero(0);
    phi.Reallocate(nb_dof_loc);
    for (int k = 0; k <= order-1; k++)
      for (int i = 0; i <= order+1; i++)
        for (int j = 0; j <= order+1; j++)
          {  
            int node = NumDofsX(k, i, j);
            phi(node).Init(ValLeg(k, kx)*ValBubbleJac(i, ky)*ValBubbleJac(j, kz), zero, zero);
            node = NumDofsY(i, k, j);
            phi(node).Init(zero, ValLeg(k, ky)*ValBubbleJac(i, kx)*ValBubbleJac(j, kz), zero);
            node = NumDofsZ(i, j, k);
            phi(node).Init(zero, zero, ValLeg(k, kz)*ValBubbleJac(i, kx)*ValBubbleJac(j, ky));
          }
    
  }
  
  
  //! computation of curl_phi = \f$ \nabla \times \phi_{node} (\xi_{kg}) \f$ 
  /*!
    \param[in] kg quadrature point number
    \param[out] curl_phi curl of basis functions 
  */ 
  void HexahedronHcurlOptimalHpFirstFamily
  ::GetCurlPhiOnQuadraturePoint(int kg, VectR3& curl_phi) const
  {
    int kx = CoordinateQuad3D(kg, 0);
    int ky = CoordinateQuad3D(kg, 1);
    int kz = CoordinateQuad3D(kg, 2);
    
    curl_phi.Reallocate(nb_dof_loc);
    FillZero(curl_phi);
    
    // dofs inside
    for (int k = 0; k <= order-1; k++)
      for (int i = 0; i <= order+1; i++)
        for (int j = 0; j <= order+1; j++)
          {  
            int node = NumDofsX(k, i, j);
            curl_phi(node)(1) = ValLeg(k, kx)*ValBubbleJac(i, ky)*DerBubbleJac(j, kz);
            curl_phi(node)(2) = -ValLeg(k, kx)*DerBubbleJac(i, ky)*ValBubbleJac(j, kz);

            node = NumDofsY(i, k, j);
            curl_phi(node)(0) = -ValLeg(k, ky)*ValBubbleJac(i, kx)*DerBubbleJac(j, kz);
            curl_phi(node)(2) = ValLeg(k, ky)*DerBubbleJac(i, kx)*ValBubbleJac(j, kz);
            
            node = NumDofsZ(i, j, k);
            curl_phi(node)(0) = ValLeg(k, kz)*ValBubbleJac(i, kx)*DerBubbleJac(j, ky);
            curl_phi(node)(1) = -ValLeg(k, kz)*DerBubbleJac(i, kx)*ValBubbleJac(j, ky);
          }
  }
  

  //! displays details of class HexahedronHcurlOptimalHpFirstFamily
  ostream& operator <<(ostream& out, const HexahedronHcurlOptimalHpFirstFamily& e)
  {
    out<<static_cast<const HexahedronReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_OPTIMAL_HP_FIRST_FAMILY_CXX
#endif
