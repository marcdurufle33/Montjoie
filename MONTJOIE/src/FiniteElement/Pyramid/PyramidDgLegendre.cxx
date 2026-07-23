#ifndef MONTJOIE_FILE_PYRAMID_DG_LEGENDRE_CXX

namespace Montjoie
{
  
  //! default constructor
  PyramidDgLegendre::PyramidDgLegendre() : PyramidReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;

    this->discontinuous_element = true;
    use_quadrature_for_sh = true;
    use_quadrature_for_rh = true;
  }

  
  size_t PyramidDgLegendre::GetMemorySize() const
  {
    size_t taille = PyramidReference<1>::GetMemorySize();
    for (int k = 0; k < 5; k++)
      {
        taille += ShLoc(k).GetMemorySize() + DerivShLoc(k).GetMemorySize();
        taille += DerivShSurf(k).GetMemorySize();
      }

    for (int k = 0; k < 4; k++)
      taille += ShSurf(k).GetMemorySize();
    
    taille += ChX.GetMemorySize() + ChY.GetMemorySize() + ChZ.GetMemorySize();
    taille += RhLoc.GetMemorySize();
    taille += DerivDxtildeDx.GetMemorySize() + DerivDxtildeDz.GetMemorySize()
      + DerivDytildeDy.GetMemorySize() + DerivDytildeDz.GetMemorySize();
    
    taille += NumFct3D.GetMemorySize() + InvWeightFunction.GetMemorySize();    
    taille += gauss_x.GetMemorySize() + gauss_z.GetMemorySize();
    taille += Seldon::GetMemorySize(ProjOperatorDxtildeQuad) + Seldon::GetMemorySize(ProjOperatorDerivShQuad);
    for (int k = 0; k < ProjOperatorDerivShTri.GetM(); k++)
      for (int p = 0; p < 4; p++)
        taille += ProjOperatorDerivShTri(k)(p).GetMemorySize();
    
    for (int k = 0; k < ProjOperatorDxtildeTri.GetM(); k++)
      for (int p = 0; p < 4; p++)
        taille += ProjOperatorDxtildeTri(k)(p).GetMemorySize();
    
    return taille;
  }

    
  //! constructing finite element
  void PyramidDgLegendre::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						 int rsurf_tri, int rsurf_quad,
						 int type_surf_tri, int type_surf_quad, int gauss_z)
  {    
    PyramidReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    ConstructFunctions();
    ConstructMassMatrix();
    ConstructStiffnessMatrix();
    
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructElementaryMatrix(*this);
  }
  
  
  //! constructing basis functions
  void PyramidDgLegendre::ConstructFunctions()
  {
    nb_dof_loc = (order+1)*(order+2)*(order+3)/6;
    InvWeightFunction.Reallocate(nb_dof_loc);
    nb_dof_loc = 0;
    NumFct3D.Reallocate(order+1, order+1, order+1);
    NumFct3D.Fill(-1);
    const Matrix<Real_wp>& CoefEvenJacobi = Fb_geom.GetCoefEvenJacobi();
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefLegendre();
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        for (int k = 0; k <= order-i-j; k++)
          {
	    NumFct3D(i, j, k) = nb_dof_loc;
	    InvWeightFunction(nb_dof_loc) = CoefLegendre(i)*CoefLegendre(j)*CoefEvenJacobi(i+j, k);
            nb_dof_loc++;
	  }
    
    nb_dof_boundaries = 0;
    
    VectR3 points_dof3d;
    points_dof3d.Reallocate(nb_points_quadrature_inside);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      points_dof3d(i) = this->PointsND(i);
    
    this->SetPointsDofND(points_dof3d);

    this->nb_points_dof_inside = nb_points_quadrature_inside;
    this->num_dof_points_surf.Reallocate(5);

    TriangleDgOrtho* tri = new TriangleDgOrtho();
    tri->ConstructFiniteElement(order);
    element_tri_surf = tri;

    QuadrangleDgGauss* quad = new QuadrangleDgGauss();
    quad->ConstructFiniteElement(order);
    element_quad_surf = quad;
  }

  
  //! constructing mass matrix
  void PyramidDgLegendre::ConstructMassMatrix()
  {
    Array3D<int> NumPrism3D(order+1, order+1, order+1);
    NumPrism3D.Fill(-1); int node = 0;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	for (int k = 0; k <= order; k++)
	  NumPrism3D(i, j, k) = node++;
    
    // first interpolation along z :
    // v^z_i,j,k = \sum_m P_m^(2max(i,j)+2)(xi_k) ((1-xi_k)/2)^(max(i,j)) u_{i,j,m}
    Array3D<Real_wp> Pv(order+1, order+1, order+1);
    Pv.Fill(0); VectReal_wp Pm;
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    const Matrix<Real_wp>& CoefEvenJacobi = Fb_geom.GetCoefEvenJacobi();
    for (int k = 0; k <= order; k++)
      {
	Real_wp xi_k = 2.0*points1d_z(k) - 1.0;
	for (int i = 0; i <= order; i++)
	  {
	    EvaluateJacobiPolynomial(EvenJacobiPolynom(i), order-i, xi_k, Pm);
	    for (int j = 0; j < Pm.GetM(); j++)
	      Pv(i, j, k) = Pm(j)*CoefEvenJacobi(i, j)*pow(0.5*(1.0-xi_k), i);
	  }
      }
    
    // CopyPvZ(Pv);
    int nb_points_z = (order+1)*(order+2)*(order+1)/2;
    Matrix<Real_wp, General, ArrayRowSparse> ch(nb_dof_loc, nb_points_z);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = NumPrism3D(i, j, k);
	    for (int m = 0; m <= order-i-j; m++)
              ch.AddInteraction(NumFct3D(i, j, m), node, Pv(i+j, m, k));
          }
    
    // conversion to csr
    Copy(ch, ChZ);
    
    // second interpolation along y :
    // v^y_i,j,k = \sum_m L_m(xi_j) v^z_{i,m,k}
    Pv.Fill(0);
    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefLegendre();
    const VectReal_wp& points1d = this->Points1D();
    for (int k = 0; k <= order; k++)
      {
	Real_wp xi_k = 2.0*points1d(k) - 1.0;
	EvaluateJacobiPolynomial(LegendrePolynom, order, xi_k, Pm);
	for (int j = 0; j < Pm.GetM(); j++)
	  Pv(0, j, k) = Pm(j)*CoefLegendre(j);
      }
    
    ch.Clear();
    ch.Reallocate(nb_points_z, nb_points_quadrature_inside);
    int rp1 = order+1;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = rp1*(rp1*i + j) + k;
	    for (int m = 0; m <= order-i; m++)
              ch.AddInteraction(NumPrism3D(i, m, k), node, Pv(0, m, j));
          }
    
    // conversion to csr
    Copy(ch, ChY);

    // third interpolation along x :
    // v^y_i,j,k = \sum_m L_m(xi_i) v^y_{m,j,k}    
    ch.Clear();
    ch.Reallocate(nb_points_quadrature_inside, nb_points_quadrature_inside);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    for (int m = 0; m <= order; m++)
              ch.AddInteraction(rp1*(m*rp1 + j) + k, rp1*(i*rp1 + j) + k, Pv(0, m, i));
          }
    
    // conversion to csr
    Copy(ch, ChX);

  }

  
  //! construction of stiffness matrix
  void PyramidDgLegendre::ConstructStiffnessMatrix()
  {
    const VectReal_wp& points1d = this->Points1D();
    const VectR3& points3d = this->PointsND();
    // Computation of gradient from values on quadrature points
    gauss_x.AffectPoints(points1d);
    gauss_z.AffectPoints(points1d_z);
    gauss_x.ComputeGradPhi(1e3*epsilon_machine);
    gauss_z.ComputeGradPhi(1e3*epsilon_machine);
    
    int N = (order+1)*(order+1)*(order+1);
    Matrix<Real_wp, General, ArrayRowSparse> rh;
    rh.Reallocate(N, 3*N);
    int rp1 = order+1;
    DerivDxtildeDx.Reallocate(points3d.GetM());
    DerivDxtildeDz.Reallocate(points3d.GetM());
    DerivDytildeDy.Reallocate(points3d.GetM());
    DerivDytildeDz.Reallocate(points3d.GetM());
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = rp1*(i*rp1 + j) + k;
	    for (int m = 0; m <= order; m++)
	      {
		int node2 = 3*(rp1*(m*rp1 + j) + k);
                rh.AddInteraction(node, node2, gauss_x.GradPhi(i, m));

		node2 = 3*(rp1*(i*rp1 + m) + k) + 1;
		rh.AddInteraction(node, node2, gauss_x.GradPhi(j, m));

		node2 = 3*(rp1*(i*rp1 + j) + m) + 2;
                rh.AddInteraction(node, node2, gauss_z.GradPhi(k, m));
	      }
	    
	    Real_wp x = points3d(node)(0), y = points3d(node)(1), z = points3d(node)(2);
	    DerivDxtildeDx(node) = 0.5/(1.0-z);
	    DerivDxtildeDz(node) = 0.5*x/square(1.0-z);
	    DerivDytildeDy(node) = 0.5/(1.0-z);
	    DerivDytildeDz(node) = 0.5*y/square(1.0-z);
	  }
    
    Copy(rh, RhLoc);
    
    // construction of ShLoc and ShLocTranspose
    int Nquad = (order+1)*(order+1);
    rh.Clear();
    rh.Reallocate(N, 5*Nquad);
    VectReal_wp ValPhiX0(order+1), ValPhiX1(order+1), ValPhiZ0(order+1);
    VectReal_wp DerivPhiX0(order+1), DerivPhiX1(order+1), DerivPhiZ0(order+1);
    
    for (int i = 0; i <= order; i++)
      {
	ValPhiX0(i) = gauss_x.EvaluatePhi(i, Real_wp(0));
	ValPhiX1(i) = gauss_x.EvaluatePhi(i, Real_wp(1));
	ValPhiZ0(i) = gauss_z.EvaluatePhi(i, Real_wp(0));

	DerivPhiX0(i) = gauss_x.EvaluatePhiGrad(i, Real_wp(0));
	DerivPhiX1(i) = gauss_x.EvaluatePhiGrad(i, Real_wp(1));
	DerivPhiZ0(i) = gauss_z.EvaluatePhiGrad(i, Real_wp(0));
      }
    
    // CopyShLoc(ValPhiX0, ValPhiX1, ValPhiZ0);
    
    Matrix<Real_wp, General, ArrayRowSparse> sh0, sh1, sh2, sh3, sh4;
    sh0.Reallocate(nb_points_quadrature_inside, Nquad);
    sh1.Reallocate(nb_points_quadrature_inside, Nquad);
    sh2.Reallocate(nb_points_quadrature_inside, Nquad);
    sh3.Reallocate(nb_points_quadrature_inside, Nquad);
    sh4.Reallocate(nb_points_quadrature_inside, Nquad);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
 	  {
	    int nv = rp1*(i*rp1 + j) + k;
	    int ns = NumQuad2D(i, j);
            sh0.AddInteraction(nv, ns, ValPhiZ0(k));
	    
	    ns = NumQuad2D(i, k);
            sh1.AddInteraction(nv, ns, ValPhiX0(j));
	    
	    ns = NumQuad2D(j, k);
            sh2.AddInteraction(nv, ns, ValPhiX1(i));
	    
	    ns = NumQuad2D(i, k);
            sh3.AddInteraction(nv, ns, ValPhiX1(j));
	    
	    ns = NumQuad2D(j, k);
            sh4.AddInteraction(nv, ns, ValPhiX0(i));
	  }
    
    Copy(sh0, ShLoc(0)); Copy(sh1, ShLoc(1)); Copy(sh2, ShLoc(2));
    Copy(sh3, ShLoc(3)); Copy(sh4, ShLoc(4));
    
    sh0.Clear(); sh1.Clear(); sh2.Clear(); sh3.Clear(); sh4.Clear();
    sh0.Reallocate(N, 2*Nquad); sh1.Reallocate(N, 2*Nquad);
    sh2.Reallocate(N, 2*Nquad); sh3.Reallocate(N, 2*Nquad); sh4.Reallocate(N, 2*Nquad);
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int nv = rp1*(i*rp1+j)+k;
	    int ns = NumQuad2D(i, j);
            sh0.AddInteraction(nv, 2*ns, ValPhiZ0(k));
            sh0.AddInteraction(nv, 2*ns+1, DerivPhiZ0(k));
            
	    ns = NumQuad2D(i, k);
            sh1.AddInteraction(nv, 2*ns, ValPhiX0(j));
            sh1.AddInteraction(nv, 2*ns+1, DerivPhiX0(j));
	    
	    ns = NumQuad2D(j, k);
            sh2.AddInteraction(nv, 2*ns, ValPhiX1(i));
            sh2.AddInteraction(nv, 2*ns+1, DerivPhiX1(i));
	    
	    ns = NumQuad2D(i, k);
            sh3.AddInteraction(nv, 2*ns, ValPhiX1(j));
            sh3.AddInteraction(nv, 2*ns+1, DerivPhiX1(j));

	    ns = NumQuad2D(j, k);
            sh4.AddInteraction(nv, 2*ns, ValPhiX0(i));
            sh4.AddInteraction(nv, 2*ns+1, DerivPhiX0(i));
	  }

    Copy(sh0, DerivShLoc(0)); Copy(sh1, DerivShLoc(1));
    Copy(sh2, DerivShLoc(2)); Copy(sh3, DerivShLoc(3)); Copy(sh4, DerivShLoc(4));
    
    // now computing ShSurf and ShSurfTranspose
    int Ntri = this->Points2D_tri().GetM();
    for (int k = 0; k < 4; k++)
      ShSurf(k).Reallocate(Nquad, Ntri); 

    sh0.Clear(); sh1.Clear(); sh2.Clear(); sh3.Clear(); sh4.Clear();
    sh0.Reallocate(2*Nquad, 3*Nquad);
    sh1.Reallocate(2*Nquad, 3*Ntri); sh3.Reallocate(2*Nquad, 3*Ntri); 
    sh2.Reallocate(2*Nquad, 3*Ntri); sh4.Reallocate(2*Nquad, 3*Ntri);

    for (int i = 0; i < Nquad; i++)
      {
        int m = nb_points_quadrature_inside + i;
        Real_wp x = points3d(m)(0), y = points3d(m)(1), z = points3d(m)(2);
        Real_wp xt = 0.5*x + 0.5, yt = 0.5*y + 0.5;
        
        // face z = 0
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_x.EvaluatePhi(j, xt)*gauss_x.EvaluatePhi(k, yt);
	      int col = NumQuad2D(j, k);
              sh0.AddInteraction(2*col, 3*i, gauss_x.EvaluatePhiGrad(j, xt)*gauss_x.EvaluatePhi(k, yt));
              sh0.AddInteraction(2*col, 3*i+1, gauss_x.EvaluatePhi(j, xt)*gauss_x.EvaluatePhiGrad(k, yt));
              sh0.AddInteraction(2*col+1, 3*i+2, vloc);
            }

        DerivDxtildeDx(m) = 0.5/(1.0-z);
        DerivDxtildeDz(m) = 0.5*x/square(1.0-z);
        DerivDytildeDy(m) = 0.5/(1.0-z);
        DerivDytildeDz(m) = 0.5*y/square(1.0-z);
      }
    
    for (int i = 0; i < Ntri; i++)
      {
	int m = nb_points_quadrature_inside + Nquad + i;
	Real_wp x = points3d(m)(0), y = points3d(m)(1), z = points3d(m)(2);
	Real_wp xt = 0.5*x/(1.0-z) + 0.5, yt = 0.5*y/(1.0-z) + 0.5, zt = z;
	
	// face y = -(1-z)
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_x.EvaluatePhi(j, xt)*gauss_z.EvaluatePhi(k, zt);
	      int col = NumQuad2D(j, k);
              ShSurf(0)(col, i) = vloc;
              sh1.AddInteraction(2*col, 3*i, gauss_x.EvaluatePhiGrad(j, xt)*gauss_z.EvaluatePhi(k, zt));
              sh1.AddInteraction(2*col, 3*i+2, gauss_x.EvaluatePhi(j, xt)*gauss_z.EvaluatePhiGrad(k, zt));
              sh1.AddInteraction(2*col+1, 3*i+1, vloc);
	    }

        DerivDxtildeDx(m) = 0.5/(1.0-z);
        DerivDxtildeDz(m) = 0.5*x/square(1.0-z);
        DerivDytildeDy(m) = 0.5/(1.0-z);
        DerivDytildeDz(m) = 0.5*y/square(1.0-z);

	// face x = (1-z)
	m = nb_points_quadrature_inside + Ntri + Nquad + i;
	x = points3d(m)(0); y = points3d(m)(1); z = points3d(m)(2);
	xt = 0.5*x/(1.0-z) + 0.5; yt = 0.5*y/(1.0-z) + 0.5; zt = z;
	
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_x.EvaluatePhi(j, yt)*gauss_z.EvaluatePhi(k, zt);
	      int col = NumQuad2D(j, k);
              ShSurf(1)(col, i) = vloc;
              sh2.AddInteraction(2*col, 3*i+1, gauss_x.EvaluatePhiGrad(j, yt)*gauss_z.EvaluatePhi(k, zt));
              sh2.AddInteraction(2*col, 3*i+2, gauss_x.EvaluatePhi(j, yt)*gauss_z.EvaluatePhiGrad(k, zt));
              sh2.AddInteraction(2*col+1, 3*i, vloc);
	    }

        DerivDxtildeDx(m) = 0.5/(1.0-z);
        DerivDxtildeDz(m) = 0.5*x/square(1.0-z);
        DerivDytildeDy(m) = 0.5/(1.0-z);
        DerivDytildeDz(m) = 0.5*y/square(1.0-z);

        // face y = (1-z)
	m = nb_points_quadrature_inside + 2*Ntri + Nquad + i;
	x = points3d(m)(0); y = points3d(m)(1); z = points3d(m)(2);
	xt = 0.5*x/(1.0-z) + 0.5; yt = 0.5*y/(1.0-z) + 0.5; zt = z;
	
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_x.EvaluatePhi(j, xt)*gauss_z.EvaluatePhi(k, zt);
	      int col = NumQuad2D(j, k);
              ShSurf(2)(col, i) = vloc;
              sh3.AddInteraction(2*col, 3*i, gauss_x.EvaluatePhiGrad(j, xt)*gauss_z.EvaluatePhi(k, zt));
              sh3.AddInteraction(2*col, 3*i+2, gauss_x.EvaluatePhi(j, xt)*gauss_z.EvaluatePhiGrad(k, zt));
              sh3.AddInteraction(2*col+1, 3*i+1, vloc);
	    }

        DerivDxtildeDx(m) = 0.5/(1.0-z);
        DerivDxtildeDz(m) = 0.5*x/square(1.0-z);
        DerivDytildeDy(m) = 0.5/(1.0-z);
        DerivDytildeDz(m) = 0.5*y/square(1.0-z);
        
        // face x = -(1-z)
	m = nb_points_quadrature_inside + 3*Ntri + Nquad + i;
	x = points3d(m)(0); y = points3d(m)(1); z = points3d(m)(2);
	xt = 0.5*x/(1.0-z) + 0.5; yt = 0.5*y/(1.0-z) + 0.5; zt = z;
		
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_x.EvaluatePhi(j, yt)*gauss_z.EvaluatePhi(k, zt);
	      int col = NumQuad2D(j, k);
              ShSurf(3)(col, i) = vloc;
              sh4.AddInteraction(2*col, 3*i+1, gauss_x.EvaluatePhiGrad(j, yt)*gauss_z.EvaluatePhi(k, zt));
              sh4.AddInteraction(2*col, 3*i+2, gauss_x.EvaluatePhi(j, yt)*gauss_z.EvaluatePhiGrad(k, zt));
              sh4.AddInteraction(2*col+1, 3*i, vloc);
	    }
        
        DerivDxtildeDx(m) = 0.5/(1.0-z);
        DerivDxtildeDz(m) = 0.5*x/square(1.0-z);
        DerivDytildeDy(m) = 0.5/(1.0-z);
        DerivDytildeDz(m) = 0.5*y/square(1.0-z);
      }
    
    Copy(sh0, DerivShSurf(0)); Copy(sh1, DerivShSurf(1));
    Copy(sh2, DerivShSurf(2)); Copy(sh3, DerivShSurf(3)); Copy(sh4, DerivShSurf(4));
  }


  //! Integration against basis functions
  /*!
    \param[in] U vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] V V_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    V = Ch U
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void PyramidDgLegendre::ApplyChGen(const Vector1& U, Vector2& V) const
  {
    //if (!UseQuadraturePointsForSh())
    //return PyramidReference::ApplyCh(U, V);
    
    Vector1 Ux(ChX.GetM()), Uy(ChY.GetM());
    Ux.Fill(0); Uy.Fill(0);
    Mlt(ChX, U, Ux);
    Mlt(ChY, Ux, Uy);
    Mlt(ChZ, Uy, V);
  }
   
  
  //! computation of u on quadrature points
  /*!
    \param[in] U components of u on degrees of freedom
    \param[out] V values of u on quadrature points
    This operation is equivalent to a matrix vector product
    V = Ch* U
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void PyramidDgLegendre::ApplyChTransposeGen(const Vector1& U, Vector2& V) const
  {
    //if (!UseQuadraturePointsForSh())
    //return PyramidReference::ApplyChTranspose(U, V);
    
    Vector1 Uy(ChY.GetN()), Uz(ChZ.GetN());
    Uy.Fill(0); Uz.Fill(0);
    Mlt(SeldonTrans, ChZ, U, Uz);
    Mlt(SeldonTrans, ChY, Uz, Uy);
    Mlt(SeldonTrans, ChX, Uy, V);    
  }
  
  
  //! Integration against gradient of basis functions associated with quadrature points
  /*!
    \param[in] U vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] V V_i = \int_K f \nabla \psi_i dx
    This operation is equivalent to a matrix vector product
    V = Rh U
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void PyramidDgLegendre::ApplyRhQuadratureGen(const Vector1& U, Vector2& V) const
  {
    // applying chaine rule to get derivatives on unit tetrahedron
    Vector1 dU(U.GetM());
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	dU(3*i) = DerivDxtildeDx(i)*U(3*i) + DerivDxtildeDz(i)*U(3*i+2);
	dU(3*i+1) = U(3*i+1)*DerivDytildeDy(i) + U(3*i+2)*DerivDytildeDz(i);
	dU(3*i+2) = U(3*i+2);
      }
    
    Mlt(RhLoc, dU, V);
  }
  
  
  //! Computation of gradient on quadrature points from values on quadrature points
  /*!
    \param[in] U values of u on quadrature points
    \param[out] V gradient of u on quadrature points
    This operation is equivalent to a matrix vector product
    V = Rh* U
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void PyramidDgLegendre::ApplyRhQuadratureTransposeGen(const Vector1& U, Vector2& V) const
  {    
    Mlt(SeldonTrans, RhLoc, U, V);
    // applying chaine rule to get derivatives on unit tetrahedron
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	V(3*i+2) += V(3*i+1)*DerivDytildeDz(i) + V(3*i)*DerivDxtildeDz(i);
	V(3*i+1) *= DerivDytildeDy(i);
	V(3*i) *= DerivDxtildeDx(i);
      }
  }
  
    
  //! integration against quadrature functions on a face
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc face number
    \param[in] U vector containing \omega_k f(\xi_k) 
    \param[out] V V_i = V_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    \psi_i are functions associated with quadrature points
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void PyramidDgLegendre
  ::ApplyShQuadratureGen(const T0& alpha, int num_loc, const Vector1& U, Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    Vector1 Uquad(Nquad);
    if ((r == 0) || (r == order))
      {
	if (num_loc == 0)
	  Copy(U, Uquad);
	else
	  Mlt(ShSurf(num_loc-1), U, Uquad);
      }
    else
      {
	if (num_loc == 0)
	  Mlt(ProjOperatorQuadOrder(r), U, Uquad);
	else
	  Mlt(ProjOperatorTriOrder(r), U, Uquad);
      }

    typename Vector2::value_type one(1);
    MltAdd(alpha, ShLoc(num_loc), Uquad, one, V); 
  }

    
  //! computation of gradient of u on quadrature points of a boundary
  /*!
    \param[in] num_loc boundary number
    \param[in] U values of u on quadrature points
    \param[out] V values of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    V = Sh* U
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
    and \varphi_i are functions associated with quadrature points
  */
  template<class Vector1, class Vector2>
  void PyramidDgLegendre
  ::ApplyShQuadratureTransposeGen(int num_loc, const Vector1& U, Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    Vector1 Uquad(Nquad);

    Mlt(SeldonTrans, ShLoc(num_loc), U, Uquad);
    
    if ((r == 0) || (r == order))
      {
	if (num_loc == 0)
	  Copy(Uquad, V);
	else
	  Mlt(SeldonTrans, ShSurf(num_loc-1), Uquad, V);
      }
    else
      {
	if (num_loc == 0)
	  Mlt(SeldonTrans, ProjOperatorQuadOrder(r), Uquad, V);
	else
	  Mlt(SeldonTrans, ProjOperatorTriOrder(r), Uquad, V);
      }
  }
  
  
  //! integration against gradient of quadrature functions on a face
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc face number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] V V_i = V_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    \psi_i are functions associated with quadrature points
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void PyramidDgLegendre
  ::ApplyNablaShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh,
			      Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    //Vector1 Uquad(Nquad);
    Vector1 Uquad, U = Uh;
    Uquad.Reallocate(2*Nquad);

    if ((r == 0) || (r == order))
      {
        int offset = nb_points_quadrature_inside + this->offset_faceSh(num_loc);
        for (int i = 0; i < GetNbQuadBoundary(num_loc); i++)
          {
            U(3*i) = DerivDxtildeDx(offset+i)*U(3*i) + DerivDxtildeDz(offset+i)*U(3*i+2);
            U(3*i+1) = U(3*i+1)*DerivDytildeDy(offset+i) + U(3*i+2)*DerivDytildeDz(offset+i);
          }
        
        Mlt(DerivShSurf(num_loc), U, Uquad);
      }
    else
      {
        if (num_loc == 0)
          {
            const Vector<TinyVector<Real_wp, 4> >& Dx = ProjOperatorDxtildeQuad(r);
            for (int i = 0; i < Dx.GetM(); i++)
              {
                U(3*i) = Dx(i)(0)*U(3*i) + Dx(i)(1)*U(3*i+2);
                U(3*i+1) = Dx(i)(2)*U(3*i+1) + Dx(i)(3)*U(3*i+2);
              }
            
            Mlt(ProjOperatorDerivShQuad(r), U, Uquad);
          }
        else
          {
            const Vector<TinyVector<Real_wp, 4> >& Dx = ProjOperatorDxtildeTri(r)(num_loc-1);
            for (int i = 0; i < Dx.GetM(); i++)
              {
                U(3*i) = Dx(i)(0)*U(3*i) + Dx(i)(1)*U(3*i+2);
                U(3*i+1) = Dx(i)(2)*U(3*i+1) + Dx(i)(3)*U(3*i+2);
              }
            
            Mlt(ProjOperatorDerivShTri(r)(num_loc-1), U, Uquad);
          }
      }

    typename Vector2::value_type one(1);
    MltAdd(alpha, DerivShLoc(num_loc), Uquad, one, V); 
  }
  
  
  //! computation of gradient of u on quadrature points of a boundary
  /*!
    \param[in] num_loc boundary number
    \param[in] U values of u on quadrature points
    \param[out] V gradient of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    V = Sh* U
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
    and \varphi_i functions associated with quadrature points
  */
  template<class Vector1, class Vector2>
  void PyramidDgLegendre
  ::ApplyNablaShQuadratureTransposeGen(int num_loc, const Vector1& U, Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    // Vector1 Uquad(Nquad);
    Vector1 Uquad;
    Uquad.Reallocate(2*Nquad);
    
    Mlt(SeldonTrans, DerivShLoc(num_loc), U, Uquad);
    
    if ((r == 0) || (r == order))
      {
        Mlt(SeldonTrans, DerivShSurf(num_loc), Uquad, V);
        
        int offset = nb_points_quadrature_inside + this->offset_faceSh(num_loc);
        for (int i = 0; i < GetNbQuadBoundary(num_loc); i++)
          {
            V(3*i+2) += V(3*i+1)*DerivDytildeDz(offset+i) + V(3*i)*DerivDxtildeDz(offset+i);
            V(3*i+1) *= DerivDytildeDy(offset+i);
            V(3*i) *= DerivDxtildeDx(offset+i);
          }
      }
    else
      {
	if (num_loc == 0)
          {
            Mlt(SeldonTrans, ProjOperatorDerivShQuad(r), Uquad, V);
            
            const Vector<TinyVector<Real_wp, 4> >& Dx = ProjOperatorDxtildeQuad(r);
            for (int i = 0; i < Dx.GetM(); i++)
              {
                V(3*i+2) += Dx(i)(3)*V(3*i+1) + Dx(i)(1)*V(3*i);
                V(3*i+1) *= Dx(i)(2);
                V(3*i) *= Dx(i)(0);
              }
          }
        else
          {
            Mlt(SeldonTrans, ProjOperatorDerivShTri(r)(num_loc-1), Uquad, V);
            
            const Vector<TinyVector<Real_wp, 4> >& Dx = ProjOperatorDxtildeTri(r)(num_loc-1);
            for (int i = 0; i < Dx.GetM(); i++)
              {
                V(3*i+2) += Dx(i)(3)*V(3*i+1) + Dx(i)(1)*V(3*i);
                V(3*i+1) *= Dx(i)(2);
                V(3*i) *= Dx(i)(0);
              }
          }
      }
  }

  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  template<class Vector1, class Vector2>
  void PyramidDgLegendre::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    Vector1 feval_weight = feval;
    const VectReal_wp& weights3d = this->WeightsND();
    for (int i = 0; i < feval_weight.GetM(); i++)
      feval_weight(i) *= weights3d(i);
    
    ApplyCh(feval_weight, contrib);
  }


  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of basis functions on point
    */  
  void PyramidDgLegendre::ComputeValuesPhiRef(const R3& pointloc, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    Real_wp a(0), b(0);
    if (abs(1.0 - z) > epsilon_machine)
      {
	a = x/(1.0-z);
	b = y/(1.0-z);
      }
    
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    
    VectReal_wp Px, Py, pow_oneMinusZ(order+1);
    EvaluateJacobiPolynomial(LegendrePolynom, order, a, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order, b, Py);
    Vector<VectReal_wp> Pz(order+1);
    pow_oneMinusZ(0) = 1.0;
    for (int i = 0; i <= order; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), order-i, 2.0*z-1, Pz(i));
    
    for (int i = 0; i < order; i++)
      pow_oneMinusZ(i+1) = pow_oneMinusZ(i)*(1.0-z);
    
    Real_wp val;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        {
          int m = i + j;
          if (m == 0)
            val = 1.0;
          else if (m == 1)
            {
              if (i == 0)
                val = y;
              else
                val = x;
            }
          else
            val = Px(i)*Py(j)*pow_oneMinusZ(m);
          
          for (int k = 0; k <= order-m; k++)
            {
	      int node = NumFct3D(i, j, k);
              phi(node) = val*Pz(m)(k);
	    }
        }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < phi.GetM(); i++)
      phi(i) *= InvWeightFunction(i);    
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on point
  */
  void PyramidDgLegendre::ComputeGradientPhiRef(const R3& pointloc, VectR3& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc);
    
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    Real_wp a(0), b(0), da_dx(0), da_dz(0), db_dy(0), db_dz(0), one(1);
    if (abs(1-z) >  1e3*epsilon_machine)
      {
	a = x/(one-z);
	b = y/(one-z);
	da_dx = one/(one-z); da_dz = x/square(one-z);
	db_dy = one/(one-z); db_dz = y/square(one-z);
      }

    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    
    VectReal_wp Px, Py, dPx, dPy, pow_oneMinusZ(order+1);
    EvaluateJacobiPolynomial(LegendrePolynom, order, a, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, order, b, Py, dPy);
    Vector<VectReal_wp> Pz(order+1), dPz(order+1);
    pow_oneMinusZ(0) = 1.0;
    for (int i = 0; i <= order; i++)      
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), order-i, 2.0*z-1, Pz(i), dPz(i));
    
    for (int i = 0; i < order; i++)
      pow_oneMinusZ(i+1) = pow_oneMinusZ(i)*(1.0-z);
    
    Real_wp val, dval_dx, dval_dz, dval_dy, valz;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        {
          int m = i + j;
          if (m == 0)
            {
              val = 1.0;
              dval_dx = 0.0; dval_dy = 0; dval_dz = 0;
            }
          else if (m == 1)
            {
              if (i == 0)
                {
                  val = y;
                  dval_dx = 0.0; dval_dy = 1.0; dval_dz = 0;
                }
              else
                {
                  val = x;
                  dval_dx = 1.0; dval_dy = 0; dval_dz = 0;
                }
            }
          else
            {
              val = Px(i)*Py(j)*pow_oneMinusZ(m);
              dval_dx = dPx(i)*Py(j)*da_dx*pow_oneMinusZ(m);
              dval_dy = dPy(j)*Px(i)*db_dy*pow_oneMinusZ(m);
              dval_dz = dPx(i)*Py(j)*da_dz*pow_oneMinusZ(m)
                + Px(i)*dPy(j)*db_dz*pow_oneMinusZ(m) - m*Px(i)*Py(j)*pow_oneMinusZ(m-1); 
            }
          
          for (int k = 0; k <= order-m; k++)
            {
	      int node = NumFct3D(i, j, k);
              valz = Pz(m)(k);
              grad_phi(node).Init(dval_dx*valz, dval_dy*valz, dval_dz*valz + 2.0*val*dPz(m)(k));
            }
        }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < grad_phi.GetM(); i++)
      grad_phi(i) *= InvWeightFunction(i);
    
  }
  
  
  //! computation of projection between finite element of different orders  
  void PyramidDgLegendre::
  ComputeTriangularInterpolationProjector(const IVect& order_elt,
                                          const ElementReference<Dimension2, 1>& Fb,
                                          const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();    
    VectReal_wp phi;
    ProjOperatorTriOrder.Reallocate(rmax+1);
    ProjOperatorDxtildeTri.Reallocate(rmax+1);
    ProjOperatorDerivShTri.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
        int Nquad = (order+1)*(order+1);
        int Ntri = Pts(i).GetM();
	ProjOperatorTriOrder(r).Reallocate(Nquad, Ntri);
        for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Real_wp yt = Pts(i)(j)(1), xt = Pts(i)(j)(0)/(1.0-yt);
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = NumQuad2D(i1, i2);
		  ProjOperatorTriOrder(r)(k, j) = gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhi(i2, yt);
		}
	  }
        
        // face y = -1 + z
        ProjOperatorDxtildeTri(r)(0).Reallocate(Ntri);
        Matrix<Real_wp, General, ArrayRowSparse> sh(2*Nquad, 3*Ntri);		
        for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Real_wp zt = Pts(i)(j)(1), xt = Pts(i)(j)(0)/(1.0-zt);
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
                {
                  int k = NumQuad2D(i1, i2);
                  sh.AddInteraction(2*k, 3*j, gauss_x.EvaluatePhiGrad(i1, xt)*gauss_z.EvaluatePhi(i2, zt));
                  sh.AddInteraction(2*k, 3*j+2, gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhiGrad(i2, zt));
                  sh.AddInteraction(2*k+1, 3*j+1, gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhi(i2, zt));
                }
            
            Real_wp x = 2.0*Pts(i)(j)(0) + Pts(i)(j)(1) - 1.0;
            Real_wp y = Pts(i)(j)(1) - 1.0;
            Real_wp z = Pts(i)(j)(1);
            ProjOperatorDxtildeTri(r)(0)(j)(0) = 0.5/(1.0-z);
            ProjOperatorDxtildeTri(r)(0)(j)(1) = 0.5*x/square(1.0-z);
            ProjOperatorDxtildeTri(r)(0)(j)(2) = 0.5/(1.0-z);
            ProjOperatorDxtildeTri(r)(0)(j)(3) = 0.5*y/square(1.0-z);
          }
        
        Copy(sh, ProjOperatorDerivShTri(r)(0));
        
        // face x = 1 - z
        ProjOperatorDxtildeTri(r)(1).Reallocate(Ntri);
        sh.Clear(); sh.Reallocate(2*Nquad, 3*Ntri);		
        for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Real_wp zt = Pts(i)(j)(1), xt = Pts(i)(j)(0)/(1.0-zt);
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
                {
                  int k = NumQuad2D(i1, i2);
                  sh.AddInteraction(2*k, 3*j+1, gauss_x.EvaluatePhiGrad(i1, xt)*gauss_z.EvaluatePhi(i2, zt));
                  sh.AddInteraction(2*k, 3*j+2, gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhiGrad(i2, zt));
                  sh.AddInteraction(2*k+1, 3*j, gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhi(i2, zt));
                }
            
            Real_wp x = 1.0 - Pts(i)(j)(1);
            Real_wp y = 2.0*Pts(i)(j)(0) + Pts(i)(j)(1) - 1.0;
            Real_wp z = Pts(i)(j)(1);
            ProjOperatorDxtildeTri(r)(1)(j)(0) = 0.5/(1.0-z);
            ProjOperatorDxtildeTri(r)(1)(j)(1) = 0.5*x/square(1.0-z);
            ProjOperatorDxtildeTri(r)(1)(j)(2) = 0.5/(1.0-z);
            ProjOperatorDxtildeTri(r)(1)(j)(3) = 0.5*y/square(1.0-z);
          }
        
        Copy(sh, ProjOperatorDerivShTri(r)(1));

        // face y = 1 - z
        ProjOperatorDxtildeTri(r)(2).Reallocate(Ntri);
        sh.Clear(); sh.Reallocate(2*Nquad, 3*Ntri);		
        for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Real_wp zt = Pts(i)(j)(1), xt = Pts(i)(j)(0)/(1.0-zt);
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
                {
                  int k = NumQuad2D(i1, i2);
                  sh.AddInteraction(2*k, 3*j, gauss_x.EvaluatePhiGrad(i1, xt)*gauss_z.EvaluatePhi(i2, zt));
                  sh.AddInteraction(2*k, 3*j+2, gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhiGrad(i2, zt));
                  sh.AddInteraction(2*k+1, 3*j+1, gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhi(i2, zt));
                }
            
            Real_wp x = 2.0*Pts(i)(j)(0) + Pts(i)(j)(1) - 1.0;
            Real_wp y = 1.0 - Pts(i)(j)(1);
            Real_wp z = Pts(i)(j)(1);
            ProjOperatorDxtildeTri(r)(2)(j)(0) = 0.5/(1.0-z);
            ProjOperatorDxtildeTri(r)(2)(j)(1) = 0.5*x/square(1.0-z);
            ProjOperatorDxtildeTri(r)(2)(j)(2) = 0.5/(1.0-z);
            ProjOperatorDxtildeTri(r)(2)(j)(3) = 0.5*y/square(1.0-z);
          }
        
        Copy(sh, ProjOperatorDerivShTri(r)(2));

        // face x = -1 + z
        ProjOperatorDxtildeTri(r)(3).Reallocate(Ntri);
        sh.Clear(); sh.Reallocate(2*Nquad, 3*Ntri);		
        for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Real_wp zt = Pts(i)(j)(1), xt = Pts(i)(j)(0)/(1.0-zt);
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
                {
                  int k = NumQuad2D(i1, i2);
                  sh.AddInteraction(2*k, 3*j+1, gauss_x.EvaluatePhiGrad(i1, xt)*gauss_z.EvaluatePhi(i2, zt));
                  sh.AddInteraction(2*k, 3*j+2, gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhiGrad(i2, zt));
                  sh.AddInteraction(2*k+1, 3*j, gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhi(i2, zt));
                }
            
            Real_wp x = -1.0 + Pts(i)(j)(1);
            Real_wp y = 2.0*Pts(i)(j)(0) + Pts(i)(j)(1) - 1.0;
            Real_wp z = Pts(i)(j)(1);
            ProjOperatorDxtildeTri(r)(3)(j)(0) = 0.5/(1.0-z);
            ProjOperatorDxtildeTri(r)(3)(j)(1) = 0.5*x/square(1.0-z);
            ProjOperatorDxtildeTri(r)(3)(j)(2) = 0.5/(1.0-z);
            ProjOperatorDxtildeTri(r)(3)(j)(3) = 0.5*y/square(1.0-z);
          }
        
        Copy(sh, ProjOperatorDerivShTri(r)(3));
        
      }
  }
  

  //! computation of projection between finite element of different orders    
  void PyramidDgLegendre::
  ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                            const ElementReference<Dimension2, 1>& Fb,
                                            const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    ProjOperatorQuadOrder.Reallocate(rmax+1);
    ProjOperatorDxtildeQuad.Reallocate(rmax+1);
    ProjOperatorDerivShQuad.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorQuadOrder(r).Reallocate((order+1)*(order+1), Pts(i).GetM());
        ProjOperatorDxtildeQuad(r).Reallocate(Pts(i).GetM());
        Matrix<Real_wp, General, ArrayRowSparse> sh(2*(order+1)*(order+1), 3*Pts(i).GetM());
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = NumQuad2D(i1, i2);
		  Real_wp yt = Pts(i)(j)(1), xt = Pts(i)(j)(0);
		  ProjOperatorQuadOrder(r)(k, j) = gauss_x.EvaluatePhi(i1, xt)*gauss_x.EvaluatePhi(i2, yt);
                  sh.AddInteraction(2*k, 3*j, gauss_x.EvaluatePhiGrad(i1, xt)*gauss_x.EvaluatePhi(i2, yt));
                  sh.AddInteraction(2*k, 3*j+1, gauss_x.EvaluatePhi(i1, xt)*gauss_x.EvaluatePhiGrad(i2, yt));
                  sh.AddInteraction(2*k+1, 3*j+2, gauss_x.EvaluatePhi(i1, xt)*gauss_x.EvaluatePhi(i2, yt));
		}	    

            Real_wp x = 2.0*Pts(i)(j)(0) - 1.0;
            Real_wp y = 2.0*Pts(i)(j)(1) - 1.0;
            Real_wp z = 0;
            ProjOperatorDxtildeQuad(r)(j)(0) = 0.5/(1.0-z);
            ProjOperatorDxtildeQuad(r)(j)(1) = 0.5*x/square(1.0-z);
            ProjOperatorDxtildeQuad(r)(j)(2) = 0.5/(1.0-z);
            ProjOperatorDxtildeQuad(r)(j)(3) = 0.5*y/square(1.0-z);
	  }
        
        Copy(sh, ProjOperatorDerivShQuad(r));
      }
  }
  
}

#define MONTJOIE_FILE_PYRAMID_DG_LEGENDRE_CXX
#endif
