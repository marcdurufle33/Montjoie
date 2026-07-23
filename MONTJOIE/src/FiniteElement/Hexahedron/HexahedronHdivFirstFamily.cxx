#ifndef MONTJOIE_FILE_HEXAHEDRON_HDIV_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronHdivFirstFamily::HexahedronHdivFirstFamily() : HexahedronReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  

  //! returns size of memory used by the object
  size_t HexahedronHdivFirstFamily::GetMemorySize() const
  {
    size_t taille = HexahedronReference<3>::GetMemorySize();
    taille += gauss.GetMemorySize();
    taille += NumDofs_X.GetMemorySize() + NumDofs_Y.GetMemorySize()
                           +NumDofs_Z.GetMemorySize()+CoordinateDofs.GetMemorySize();
    
    taille += G_GL.GetMemorySize()+GL_G.GetMemorySize()+dGL_G.GetMemorySize()
      +dGL_GL.GetMemorySize()+stiff1d.GetMemorySize();
    
    taille += rh_loc.GetMemorySize()+ch_loc.GetMemorySize()+ch1_loc.GetMemorySize();
    taille += ch2_loc.GetMemorySize()+ch3_loc.GetMemorySize();
    taille += rh_exact.GetMemorySize()+ch1_node.GetMemorySize()+ch2_node.GetMemorySize();
    taille += ch3_node.GetMemorySize() + gauss_exact.GetMemorySize();
    taille += G_GX.GetMemorySize()+GL_GX.GetMemorySize()+dGX_GX.GetMemorySize()
      +weights3d_Exact.GetMemorySize()+WeightsMassG.GetMemorySize()+invWeightsMassG.GetMemorySize()
      +invSqrtWeightsMassG.GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    return taille;
  }

  
  //! in order to number a mesh
  void HexahedronHdivFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return HexahedronReference<3>::ConstructNumberMap(nmap, dg);
    
    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, 0);
    nmap.SetNbDofQuadrangle(order, order*order);
    
    nmap.SetNbDofHexahedron(order, nb_dof_loc - 6*order*order);

    Matrix<int> FacesDof_Rotation_Quad;
    MeshNumbering<Dimension3>::
      GetRotationQuadrilateralFace(this->PointsDof2D_quad(), FacesDof_Rotation_Quad);

    nmap.SetFacesDofRotationQuad(this->order, FacesDof_Rotation_Quad);
    
    Matrix<bool> OppositeSigneDof(8, order*order);
    OppositeSigneDof.Fill(false);
    for (int i = 4; i < 8; i++)
      for (int j = 0; j < order*order; j++)
	OppositeSigneDof(i, j) = true;
    
    nmap.SetSignDofRotationQuad(order, OppositeSigneDof);
  }

  
  //! constructing finite element
  void HexahedronHdivFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf_tri, int rsurf_quad,
			   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    HexahedronReference<3>::ConstructFiniteElement(r, rgeom, rquad, type_quad); 
    
    lob_basis.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    gauss.ConstructQuadrature(order-1, Globatto<Real_wp>::QUADRATURE_GAUSS);
    gauss_exact.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    ConstructFunctions();
    this->nb_points_dof_inside = this->nb_dof_loc;
    
    ConstructStiffnessMatrix();
    ConstructMassMatrix();

    ConstructHdivElementaryMatrix();
    
    FindDofsOnFace();
    
    this->Fb_geom.ComputeCoefficientTransformation();

    QuadrangleGauss* Fb_quad = new QuadrangleGauss();
    Fb_quad->ConstructFiniteElement(order);
    element_quad_surf = Fb_quad;
    
    this->num_dof_points_surf.Reallocate(6);
    for (int n = 0; n < 6; n++)
      {
	this->num_dof_points_surf(n).Reallocate(this->GetNbDofBoundary(n));
	for (int i = 0; i < this->GetNbDofBoundary(n); i++)
	  this->num_dof_points_surf(n)(i) = FacesDof(i, n);
      }    
  }

  
  //! construction of orthogonal functions generating the finite element space
  void HexahedronHdivFirstFamily::ConstructOrthogonalBasis(int r)
  {
    order = r;
    nb_dof_loc = 3*(order+1)*order*order;
    GetJacobiPolynomial(LegendrePolynom, r+1, Real_wp(0), Real_wp(0));
  }
  
  
  //! evaluating orthogonal functions generating the finite element space
  void HexahedronHdivFirstFamily::ComputeValuesPhiOrthoRef(const R3& x, VectR3& phi)
  {
    VectReal_wp Px, Py, Pz;
    int r = order;
    EvaluateJacobiPolynomial(LegendrePolynom, r, 2.0*x(0) - 1.0, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, r, 2.0*x(1) - 1.0, Py);
    EvaluateJacobiPolynomial(LegendrePolynom, r, 2.0*x(2) - 1.0, Pz);
    
    phi.Reallocate(3*(order+1)*order*order);
    int num = 0;
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-1; j++)
        for (int k = 0; k <= r-1; k++)
          {
            phi(num++).Init(Px(i)*Py(j)*Pz(k), 0.0, 0.0);
            phi(num++).Init(0.0, Py(i)*Px(j)*Pz(k), 0.0);
            phi(num++).Init(0.0, 0.0, Pz(i)*Px(j)*Py(k));
          }
          
  }
  

  //! construction of basis functions
  void HexahedronHdivFirstFamily::ConstructFunctions()
  {
    nb_dof_quad = order*order;
    nb_dof_loc = 3*order*order*(order+1);
    nb_dof_tri = 0;
    
    int nb_dof_inside_quad = nb_dof_quad;
    nb_dof_boundaries = 6*nb_dof_inside_quad;

    FacesDof.Reallocate(nb_dof_quad, 6);
    NumDofs_X.Reallocate(order+1, order, order);
    NumDofs_Y.Reallocate(order, order+1, order);
    NumDofs_Z.Reallocate(order, order, order+1);
    
    VectR2 points_dof2d_quad; VectR3 points_dof3d;
    points_dof2d_quad.Reallocate(order*order);
    for (int i = 0; i < order; i++)
      for (int j = 0; j < order; j++)
	points_dof2d_quad(i*order + j).Init(gauss.Points(i), gauss.Points(j));

    // dofs on first face
    int offset = 0;
    for (int i = 0; i < order ; i++)
      for (int j = 0 ; j < order ; j++)
	NumDofs_X(0, i, j) = offset + i*order + j;
    
    // dofs on second face
    offset += nb_dof_quad;
    for (int i = 0; i < order ; i++)
      for (int j = 0 ; j < order ; j++)
	NumDofs_Y(i, 0, j) = offset + i*order + j;

    // dofs on third face
    offset += nb_dof_quad;
    for (int i = 0; i < order ; i++)
      for (int j = 0 ; j < order ; j++)
	NumDofs_Z(i, j, 0) = offset + i*order + j;
    
    // dofs on fourth face
    offset += nb_dof_quad;
    for (int i = 0; i < order ; i++)
      for (int j = 0 ; j < order ; j++)
	NumDofs_Z(i, j, order) = offset + i*order + j;
    
    // dofs on fifth face
    offset += nb_dof_quad;
    for (int i = 0; i < order ; i++)
      for (int j = 0 ; j < order ; j++)
	NumDofs_Y(i, order, j) = offset + i*order + j;
    
    // dofs on sixth face
    offset += nb_dof_quad;
    for (int i = 0; i < order ; i++)
      for (int j = 0 ; j < order ; j++)
	NumDofs_X(order, i, j) = offset + i*order + j;
    
    // dofs inside the hexahedron
    offset += nb_dof_quad;
    for (int i = 1; i < order; i++)
      for (int j = 0; j < order; j++)
	for (int k = 0; k < order; k++)
	  {
	    NumDofs_X(i, j, k) = offset++;
	    NumDofs_Y(j, i, k) = offset++;
	    NumDofs_Z(j, k, i) = offset++;
	  }
    
    points_dof3d.Reallocate(nb_dof_loc);
    CoordinateDofs.Reallocate(nb_dof_loc, 4);
    
    // DISP(NumDofs_X); DISP(NumDofs_Y); DISP(NumDofs_Z);
    for (int i = 0 ; i <= order ; i++)
      for (int j = 0 ; j < order ; j++)
	for (int k = 0 ; k < order ; k++)
	  {
	    points_dof3d(NumDofs_X(i, j, k))
              .Init(lob_basis.Points(i), gauss.Points(j), gauss.Points(k));
	    points_dof3d(NumDofs_Y(j, i, k))
              .Init(gauss.Points(j), lob_basis.Points(i), gauss.Points(k));
	    points_dof3d(NumDofs_Z(k, j, i))
              .Init(gauss.Points(k), gauss.Points(j), lob_basis.Points(i));
	    
	    CoordinateDofs(NumDofs_X(i, j, k), 0) = i;
	    CoordinateDofs(NumDofs_X(i, j, k), 1) = j;
	    CoordinateDofs(NumDofs_X(i, j, k), 2) = k;
	    CoordinateDofs(NumDofs_X(i, j, k), 3) = 0;
	    CoordinateDofs(NumDofs_Y(j, i, k), 0) = j;
	    CoordinateDofs(NumDofs_Y(j, i, k), 1) = i;
	    CoordinateDofs(NumDofs_Y(j, i, k), 2) = k;
	    CoordinateDofs(NumDofs_Y(j, i, k), 3) = 1;
	    CoordinateDofs(NumDofs_Z(k, j, i), 0) = k;
	    CoordinateDofs(NumDofs_Z(k, j, i), 1) = j;
	    CoordinateDofs(NumDofs_Z(k, j, i), 2) = i;
	    CoordinateDofs(NumDofs_Z(k, j, i), 3) = 2;
	  }
    
    this->SetPointsDof2D_quad(points_dof2d_quad);
    this->SetPointsDofND(points_dof3d);
  }
  
  
  //! computation of some arrays, G_GL, GL_G and dGL_G
  //! needed for the computation of stiffness matrix
  void HexahedronHdivFirstFamily::ConstructStiffnessMatrix()
  {
    G_GL.Reallocate(order, order+1);
    GL_G.Reallocate(order+1, order);
    dGL_G.Reallocate(order+1, order);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	{
	  G_GL(i, j) = gauss.EvaluatePhi(i, lob_basis.Points(j));
	  GL_G(j, i) = lob_basis.EvaluatePhi(j, gauss.Points(i));
	  dGL_G(j, i) = lob_basis.EvaluatePhiGrad(j, gauss.Points(i));
	}
    
    lob_basis.ComputeGradPhi(1e3*epsilon_machine);
    dGL_GL = lob_basis.GradPhi();
    
    stiff1d.Reallocate(order+1, order+1);
    stiff1d.Fill(0);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          stiff1d(i, j) += this->Weights1D(k)*dGL_GL(i, k)*dGL_GL(j, k);

    G_GX.Reallocate(order, order+1); GL_GX.Reallocate(order+1, order+1);
    dGX_GX.Reallocate(order+1, order+1);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	G_GX(i,j) = gauss.EvaluatePhi(i, gauss_exact.Points(j));
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  GL_GX(i,j) = lob_basis.EvaluatePhi(i, gauss_exact.Points(j));
	  dGX_GX(i,j) = gauss_exact.EvaluatePhiGrad(i, gauss_exact.Points(j));
	}
    
    int N = nb_points_quadrature_inside;
    weights3d_Exact.Reallocate(N);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  weights3d_Exact(NumQuad3D(i,j,k)) = gauss_exact.Weights(i)
	    *gauss_exact.Weights(j)*gauss_exact.Weights(k);
    
  }
  
  
  //! mass matrix
  void HexahedronHdivFirstFamily::ConstructMassMatrix()
  {    
    Matrix<Real_wp> G_Geom(order, order+1), GL_Geom(order+1, order+1);
    int order_geom = this->GetGeometryOrder();
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
        GL_Geom(i, j) = lob_basis.EvaluatePhi(i, this->PointsNodal1D(j));
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order_geom; j++)
        G_Geom(i, j) = gauss.EvaluatePhi(i, this->PointsNodal1D(j));    
  }

  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point
  */
  void HexahedronHdivFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc); res.Fill(R3(0,0,0));
    VectReal_wp phiGx(order), phiGy(order), phiGz(order),
      phiGLx(order+1), phiGLy(order+1), phiGLz(order+1);
    for (int i = 0; i < order; i++)
      {
	phiGx(i) = gauss.EvaluatePhi(i, point_loc(0));
	phiGy(i) = gauss.EvaluatePhi(i, point_loc(1));
	phiGz(i) = gauss.EvaluatePhi(i, point_loc(2));
      }

    for (int i = 0; i <= order; i++)
      {
	phiGLx(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiGLy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
	phiGLz(i) = lob_basis.EvaluatePhi(i, point_loc(2));
      }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j < order; j++)
	for (int k = 0; k < order; k++)
	  {
	    int nx = NumDofs_X(i, j, k);
	    res(nx)(0) = phiGLx(i)*phiGy(j)*phiGz(k);
	    
	    int ny = NumDofs_Y(j, i, k);
	    res(ny)(1) = -phiGx(j)*phiGLy(i)*phiGz(k);
	    
	    int nz = NumDofs_Z(k, j, i);
	    res(nz)(2) = phiGx(k)*phiGy(j)*phiGLz(i);
	  }
  }

  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] div divergence of basis functions on point
  */
  void HexahedronHdivFirstFamily::ComputeDivPhiRef(const R3& point_loc, VectReal_wp& div) const
  {
    div.Reallocate(nb_dof_loc); div.Fill(0);
    VectReal_wp phiGx(order), phiGy(order), phiGz(order);
    VectReal_wp dphiGLx(order+1), dphiGLy(order+1), dphiGLz(order+1);
    for (int i = 0; i < order; i++)
      {
	phiGx(i) = gauss.EvaluatePhi(i, point_loc(0));
	phiGy(i) = gauss.EvaluatePhi(i, point_loc(1));
	phiGz(i) = gauss.EvaluatePhi(i, point_loc(2));
      }

    for (int i = 0; i <= order; i++)
      {
	dphiGLx(i) = lob_basis.EvaluatePhiGrad(i, point_loc(0));
	dphiGLy(i) = lob_basis.EvaluatePhiGrad(i, point_loc(1));
	dphiGLz(i) = lob_basis.EvaluatePhiGrad(i, point_loc(2));
      }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j < order; j++)
	for (int k = 0; k < order; k++)
	  {
	    int nx = NumDofs_X(i, j, k);
	    div(nx) = dphiGLx(i)*phiGy(j)*phiGz(k);
	    
	    int ny = NumDofs_Y(j, i, k);
	    div(ny) = -dphiGLy(i)*phiGx(j)*phiGz(k);

	    int nz = NumDofs_Z(k, j, i);
	    div(nz) = dphiGLz(i)*phiGx(k)*phiGy(j);
	  }
    
  }

  
  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are placed
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void HexahedronHdivFirstFamily::
  ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    contrib.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      {
	if (CoordinateDofs(i, 3) == 1)
	  contrib(i) = -feval(3*i + CoordinateDofs(i,3));
	else
	  contrib(i) = feval(3*i + CoordinateDofs(i,3));
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
  void HexahedronHdivFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    int Nquad = (order+1)*(order+1);
    Vector1 fx(3*Nquad), fy(3*Nquad);
    FillZero(fx); FillZero(fy);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
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


  void HexahedronHdivFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
    if (num_loc%2 == 0)
      for (int i = 0; i < contrib.GetM(); i++)
	contrib(i) = -contrib(i);
  }
  

  void HexahedronHdivFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    if (num_loc %= 2)
      for (int i = 0; i < contrib.GetM(); i++)
	contrib(i) = -contrib(i);
  }
  

  //! displays details of class HexahedronHdivFirstFamily
  ostream& operator <<(ostream& out, const HexahedronHdivFirstFamily& e)
  {
    out << static_cast<const HexahedronReference<3>& >(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_HEXAHEDRON_HDIV_FIRST_FAMILY_CXX
#endif
