#ifndef MONTJOIE_FILE_TETRAHEDRON_DG_ORTHO_CXX

namespace Montjoie
{
  //! default constructor
  TetrahedronDgOrtho::TetrahedronDgOrtho() : TetrahedronReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;

    this->discontinuous_element = true;
    use_quadrature_for_sh = true;
    use_quadrature_for_rh = true;
  }


  //! how to number mesh
  void TetrahedronDgOrtho::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg != ElementReference_Base::HDG)
      return TetrahedronReference<1>::ConstructNumberMap(nmap, dg);
    
    nmap.SetFormulationDG(dg);
    
    nmap.SetNbDofVertex(this->order, 0);
    nmap.SetNbDofEdge(this->order, 0);
    nmap.SetNbDofTriangle(this->order, (this->order+2)*(this->order+1)/2);
    nmap.SetNbDofQuadrangle(this->order, (this->order+1)*(this->order+1));
    nmap.SetNbDofTetrahedron(this->order, 0);
    
    // rotation on triangular faces
    ElementReference<Dimension2, 1>::
      FindH1RotationTri(this->order, this->Points2D_tri(), this->Weights2D_tri(),
                        element_tri_surf->GetValuePhi(), nmap);
    
  }

  
  size_t TetrahedronDgOrtho::GetMemorySize() const
  {
    size_t taille = TetrahedronReference<1>::GetMemorySize();
    for (int k = 0; k < 4; k++)
      {
        taille += ShLoc(k).GetMemorySize() + DerivShLoc(k).GetMemorySize();
        taille += ShSurf(k).GetMemorySize();
        taille += DerivShSurf(k).GetMemorySize();
      }
    
    taille += ChX.GetMemorySize() + ChY.GetMemorySize() + ChZ.GetMemorySize();
    taille += RhLoc.GetMemorySize() + InvWeightFct.GetMemorySize(); 
    taille += points1d_y.GetMemorySize() + points1d_z.GetMemorySize();
    taille += DerivDxtildeDx.GetMemorySize() + DerivDxtildeDy.GetMemorySize()
      + DerivDytildeDy.GetMemorySize() + DerivDytildeDz.GetMemorySize();
    
    taille += gauss_x.GetMemorySize() + gauss_y.GetMemorySize() + gauss_z.GetMemorySize();
    for (int k = 0; k < ProjOperatorSh.GetM(); k++)
      for (int p = 0; p < 4; p++)
        taille += ProjOperatorSh(k)(p).GetMemorySize();

    for (int k = 0; k < ProjOperatorDerivSh.GetM(); k++)
      for (int p = 0; p < 4; p++)
        taille += ProjOperatorDerivSh(k)(p).GetMemorySize();
    
    for (int k = 0; k < ProjOperatorDxtilde.GetM(); k++)
      for (int p = 0; p < 4; p++)
        taille += ProjOperatorDxtilde(k)(p).GetMemorySize();
    
    return taille;
  }
  
  
  //! constructing finite element
  void TetrahedronDgOrtho::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						  int rsurf_tri, int rsurf_quad,
						  int type_surf_tri, int type_surf_quad, int gauss_z)
  {    
    TetrahedronReference<1>::
      ConstructFiniteElement(r, rgeom, rquad, TetrahedronQuadrature::QUADRATURE_TENSOR);
    
    ConstructFunctions();
    ConstructMassMatrix();
    ConstructStiffnessMatrix();
    
    this->Fb_geom.ComputeCoefficientTransformation();
    
    this->ConstructElementaryMatrix(*this);
  }
  
  
  //! construction of mass matrix
  void TetrahedronDgOrtho::ConstructMassMatrix()
  {
    // computation of ChX, ChY, ChZ
    
    // first interpolation along z :
    // v^z_i,j,k = \sum_m P_m^(2i+2j+2)(xi_k) ((1-xi_k)/2)^(i+j) u_{i,j,m}
    Array3D<Real_wp> Pv(order+1, order+1, order+1);
    Pv.Fill(0); VectReal_wp Pn;
    const Vector<Matrix<Real_wp> >& OddJacobiPolynom = Fb_geom.GetOddJacobiPolynomial();
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefOddJacobi();
    const Matrix<Real_wp>& CoefEvenJacobi = Fb_geom.GetCoefEvenJacobi();
    for (int k = 0; k <= order; k++)
      {
	Real_wp xi_k = 2.0*points1d_z(k) - 1.0;
	for (int i = 0; i <= order; i++)
	  {
	    EvaluateJacobiPolynomial(EvenJacobiPolynom(i), order-i, xi_k, Pn);
	    for (int j = 0; j < Pn.GetM(); j++)
	      Pv(i, j, k) = Pn(j)*CoefEvenJacobi(i, j)*pow(0.5*(1.0-xi_k), i);
	  }
      }
    
    int nb_points_z = (order+1)*(order+1)*(order+2)/2;
    Matrix<Real_wp, General, ArrayRowSparse> ch(nb_dof_loc, nb_points_z);
    Array3D<int> NumNodesZ(order+1, order+1, order+1); int node = 0;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	for (int k = 0; k <= order; k++)
	  {
	    for (int m = 0; m <= order-i-j; m++)
              ch.AddInteraction(NumFct3D(i, j, m), node, Pv(i+j, m, k));
            
	    NumNodesZ(i, j, k) = node++;
	  }
    
    // conversion to csr
    Copy(ch, ChZ);
    
    // second interpolation along y :
    // v^y_i,j,k = \sum_m P_m^(2i+1)(xi_j) ((1-xi_j)/2)^i v^z_{i,j,m}
    Pv.Fill(0);
    for (int k = 0; k <= order; k++)
      {
	Real_wp xi_k = 2.0*points1d_y(k) - 1.0;
	for (int i = 0; i <= order; i++)
	  {
	    EvaluateJacobiPolynomial(OddJacobiPolynom(i), order-i, xi_k, Pn);
	    for (int j = 0; j < Pn.GetM(); j++)
	      Pv(i, j, k) = Pn(j)*CoefOddJacobi(i, j)*pow(0.5*(1.0-xi_k), i);
	  }
      }
    
    int nb_points_y = (order+1)*(order+1)*(order+1);
    ch.Clear();
    ch.Reallocate(nb_points_z, nb_points_y);
    Array3D<int> NumNodesY(order+1, order+1, order+1); node = 0;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    for (int m = 0; m <= order-i; m++)
              ch.AddInteraction(NumNodesZ(i, m, k), node, Pv(i, m, j));
	
	    NumNodesY(i, j, k) = node++;
	  }
    
    // conversion to csr
    Copy(ch, ChY);

    // third interpolation along x :
    // v^y_i,j,k = \sum_m L_m(xi_i) v^y_{m,j,k}
    Pv.Fill(0);
    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefLegendre();
    const VectReal_wp& points1d = this->Points1D();
    for (int k = 0; k <= order; k++)
      {
	Real_wp xi_k = 2.0*points1d(k) - 1.0;
	EvaluateJacobiPolynomial(LegendrePolynom, order, xi_k, Pn);
	for (int j = 0; j < Pn.GetM(); j++)
	  Pv(0, j, k) = Pn(j)*CoefLegendre(j);
      }
    
    ch.Clear();
    ch.Reallocate(nb_points_y, nb_points_y);
    int rp1 = order+1;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    for (int m = 0; m <= order; m++)
              ch.AddInteraction(NumNodesY(m, j, k), rp1*(i*rp1+j)+k, Pv(0, m, i));
	  }
    
    // conversion to csr
    Copy(ch, ChX);
    
  }
  
  
  //! construction of stiffness matrix
  void TetrahedronDgOrtho::ConstructStiffnessMatrix()
  {
    const VectReal_wp& points1d = this->Points1D();
    const VectR3& points3d = this->PointsND();
    
    // Computation of gradient from values on quadrature points
    gauss_x.AffectPoints(points1d);
    gauss_y.AffectPoints(points1d_y);
    gauss_z.AffectPoints(points1d_z);
    gauss_x.ComputeGradPhi(1e3*epsilon_machine);
    gauss_y.ComputeGradPhi(1e3*epsilon_machine);
    gauss_z.ComputeGradPhi(1e3*epsilon_machine);
    
    int N = (order+1)*(order+1)*(order+1);
    Matrix<Real_wp, General, ArrayRowSparse> rh;
    rh.Reallocate(N, 3*N);
    DerivDxtildeDx.Reallocate(points3d.GetM());
    DerivDxtildeDy.Reallocate(points3d.GetM());
    DerivDytildeDy.Reallocate(points3d.GetM());
    DerivDytildeDz.Reallocate(points3d.GetM());
    int rp1 = order+1;
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
                rh.AddInteraction(node, node2, gauss_y.GradPhi(j, m));

		node2 = 3*(rp1*(i*rp1 + j) + m) + 2;
                rh.AddInteraction(node, node2, gauss_z.GradPhi(k, m));
	      }
	    
	    Real_wp x = points3d(node)(0), y = points3d(node)(1), z = points3d(node)(2);
	    DerivDxtildeDx(node) = 1.0/(1.0-y-z);
	    DerivDxtildeDy(node) = x/square(1.0-y-z);
	    DerivDytildeDy(node) = 1.0/(1.0-z);
	    DerivDytildeDz(node) = y/square(1.0-z);
	  }
    
    Copy(rh, RhLoc);
    
    // construction of ShLoc
    int Nquad = (order+1)*(order+1);
    rh.Clear();
    
    Matrix<Real_wp, General, ArrayRowSparse> sh0, sh1, sh2, sh3;
    sh0.Reallocate(N, Nquad); sh1.Reallocate(N, Nquad);
    sh2.Reallocate(N, Nquad); sh3.Reallocate(N, Nquad);

    VectReal_wp ValPhiX0, ValPhiX1, ValPhiY0, ValPhiZ0;
    VectReal_wp DerivPhiX0, DerivPhiX1, DerivPhiY0, DerivPhiZ0;
    
    ValPhiX0.Reallocate(order+1); ValPhiX1.Reallocate(order+1);
    ValPhiY0.Reallocate(order+1); ValPhiZ0.Reallocate(order+1);
    DerivPhiX0.Reallocate(order+1); DerivPhiX1.Reallocate(order+1);
    DerivPhiY0.Reallocate(order+1); DerivPhiZ0.Reallocate(order+1);
    
    for (int i = 0; i <= order; i++)
      {
	ValPhiX0(i) = gauss_x.EvaluatePhi(i, Real_wp(0));
	ValPhiX1(i) = gauss_x.EvaluatePhi(i, Real_wp(1));
	ValPhiY0(i) = gauss_y.EvaluatePhi(i, Real_wp(0));
	ValPhiZ0(i) = gauss_z.EvaluatePhi(i, Real_wp(0));

	DerivPhiX0(i) = gauss_x.EvaluatePhiGrad(i, Real_wp(0));
	DerivPhiX1(i) = gauss_x.EvaluatePhiGrad(i, Real_wp(1));
	DerivPhiY0(i) = gauss_y.EvaluatePhiGrad(i, Real_wp(0));
	DerivPhiZ0(i) = gauss_z.EvaluatePhiGrad(i, Real_wp(0));
      }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int nv = rp1*(i*rp1+j)+k;
	    int ns = rp1*i + j;
            sh0.AddInteraction(nv, ns, ValPhiZ0(k));
	    
	    ns = rp1*i + k;
            sh1.AddInteraction(nv, ns, ValPhiY0(j));
	    
	    ns = rp1*j + k;
            sh2.AddInteraction(nv, ns, ValPhiX0(i));
	    
	    ns = rp1*j + k;
            sh3.AddInteraction(nv, ns, ValPhiX1(i));
	  }
    
    Copy(sh0, ShLoc(0)); Copy(sh1, ShLoc(1));
    Copy(sh2, ShLoc(2)); Copy(sh3, ShLoc(3));

    sh0.Clear(); sh1.Clear(); sh2.Clear(); sh3.Clear();
    sh0.Reallocate(N, 2*Nquad); sh1.Reallocate(N, 2*Nquad);
    sh2.Reallocate(N, 2*Nquad); sh3.Reallocate(N, 2*Nquad);

    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int nv = rp1*(i*rp1+j)+k;
	    int ns = rp1*i + j;
            sh0.AddInteraction(nv, 2*ns, ValPhiZ0(k));
            sh0.AddInteraction(nv, 2*ns+1, DerivPhiZ0(k));
            
	    ns = rp1*i + k;
            sh1.AddInteraction(nv, 2*ns, ValPhiY0(j));
            sh1.AddInteraction(nv, 2*ns+1, DerivPhiY0(j));
	    
	    ns = rp1*j + k;
            sh2.AddInteraction(nv, 2*ns, ValPhiX0(i));
            sh2.AddInteraction(nv, 2*ns+1, DerivPhiX0(i));
	    
	    ns = rp1*j + k;
            sh3.AddInteraction(nv, 2*ns, ValPhiX1(i));
            sh3.AddInteraction(nv, 2*ns+1, DerivPhiX1(i));
	  }

    Copy(sh0, DerivShLoc(0)); Copy(sh1, DerivShLoc(1));
    Copy(sh2, DerivShLoc(2)); Copy(sh3, DerivShLoc(3));
        
    // now computing ShSurf
    int Ntri = this->Points2D_tri().GetM();
    for (int k = 0; k < 4; k++)
      ShSurf(k).Reallocate(Nquad, Ntri); 
    
    sh0.Clear(); sh1.Clear(); sh2.Clear(); sh3.Clear();
    sh0.Reallocate(2*Nquad, 3*Ntri); sh1.Reallocate(2*Nquad, 3*Ntri); 
    sh2.Reallocate(2*Nquad, 3*Ntri); sh3.Reallocate(2*Nquad, 3*Ntri); 
    
    for (int i = 0; i < Ntri; i++)
      {	
	// face z = 0
        int m = nb_points_quadrature_inside + i;
	Real_wp x = points3d(m)(0), y = points3d(m)(1), z = points3d(m)(2);
	Real_wp xt = x/(1.0-y-z), yt = y/(1.0-z), zt = z;

        for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_x.EvaluatePhi(j, xt)*gauss_y.EvaluatePhi(k, yt);
	      int col = rp1*j + k;
              ShSurf(0)(col, i) = vloc;
              sh0.AddInteraction(2*col, 3*i, gauss_x.EvaluatePhiGrad(j, xt)*gauss_y.EvaluatePhi(k, yt));
              sh0.AddInteraction(2*col, 3*i+1, gauss_x.EvaluatePhi(j, xt)*gauss_y.EvaluatePhiGrad(k, yt));
              sh0.AddInteraction(2*col+1, 3*i+2, vloc);
	    }

        DerivDxtildeDx(m) = 1.0/(1.0-y-z);
        DerivDxtildeDy(m) = x/square(1.0-y-z);
        DerivDytildeDy(m) = 1.0/(1.0-z);
        DerivDytildeDz(m) = y/square(1.0-z);

	// face y = 0
	m = nb_points_quadrature_inside + Ntri + i;
	x = points3d(m)(0); y = points3d(m)(1); z = points3d(m)(2);
	xt = x/(1.0-y-z); yt = y/(1.0-z); zt = z;
	
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_x.EvaluatePhi(j, xt)*gauss_z.EvaluatePhi(k, zt);
	      int col = rp1*j + k;
              ShSurf(1)(col, i) = vloc;
              sh1.AddInteraction(2*col, 3*i, gauss_x.EvaluatePhiGrad(j, xt)*gauss_z.EvaluatePhi(k, zt));
              sh1.AddInteraction(2*col, 3*i+2, gauss_x.EvaluatePhi(j, xt)*gauss_z.EvaluatePhiGrad(k, zt));
              sh1.AddInteraction(2*col+1, 3*i+1, vloc);
	    }

        DerivDxtildeDx(m) = 1.0/(1.0-y-z);
        DerivDxtildeDy(m) = x/square(1.0-y-z);
        DerivDytildeDy(m) = 1.0/(1.0-z);
        DerivDytildeDz(m) = y/square(1.0-z);

	// face x = 0
	m = nb_points_quadrature_inside + 2*Ntri + i;
	x = points3d(m)(0); y = points3d(m)(1); z = points3d(m)(2);
	xt = x/(1.0-y-z); yt = y/(1.0-z); zt = z;
	
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_y.EvaluatePhi(j, yt)*gauss_z.EvaluatePhi(k, zt);
	      int col = rp1*j + k;
              ShSurf(2)(col, i) = vloc;
              sh2.AddInteraction(2*col, 3*i+1, gauss_y.EvaluatePhiGrad(j, yt)*gauss_z.EvaluatePhi(k, zt));
              sh2.AddInteraction(2*col, 3*i+2, gauss_y.EvaluatePhi(j, yt)*gauss_z.EvaluatePhiGrad(k, zt));
              sh2.AddInteraction(2*col+1, 3*i, vloc);
	    }
        
        DerivDxtildeDx(m) = 1.0/(1.0-y-z);
        DerivDxtildeDy(m) = x/square(1.0-y-z);
        DerivDytildeDy(m) = 1.0/(1.0-z);
        DerivDytildeDz(m) = y/square(1.0-z);
        
	// face x+y+z = 1
	m = nb_points_quadrature_inside + 3*Ntri + i;
	x = points3d(m)(0); y = points3d(m)(1); z = points3d(m)(2);
	xt = x/(1.0-y-z); yt = y/(1.0-z); zt = z;
	
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_y.EvaluatePhi(j, yt)*gauss_z.EvaluatePhi(k, zt);
	      int col = rp1*j + k;
              ShSurf(3)(col, i) = vloc;
              sh3.AddInteraction(2*col, 3*i+1, gauss_y.EvaluatePhiGrad(j, yt)*gauss_z.EvaluatePhi(k, zt));
              sh3.AddInteraction(2*col, 3*i+2, gauss_y.EvaluatePhi(j, yt)*gauss_z.EvaluatePhiGrad(k, zt));
              sh3.AddInteraction(2*col+1, 3*i, vloc);
	    }
        
        DerivDxtildeDx(m) = 1.0/(1.0-y-z);
        DerivDxtildeDy(m) = x/square(1.0-y-z);
        DerivDytildeDy(m) = 1.0/(1.0-z);
        DerivDytildeDz(m) = y/square(1.0-z);
      }

    Copy(sh0, DerivShSurf(0)); Copy(sh1, DerivShSurf(1));
    Copy(sh2, DerivShSurf(2)); Copy(sh3, DerivShSurf(3));
    
  }


  //! constructing basis functions
  void TetrahedronDgOrtho::ConstructFunctions()
  {
    int r = order;
    InvWeightFct.Reallocate((r+1)*(r+2)*(r+3)/6);
    NumFct3D.Reallocate(r+1, r+1, r+1); NumFct3D.Fill(-1);
    int nb = 0;
    const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefOddJacobi();
    const Matrix<Real_wp>& CoefEvenJacobi = Fb_geom.GetCoefEvenJacobi();
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefLegendre();
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        for (int k = 0; k <= r-i-j; k++)
          {
            NumFct3D(i, j, k) = nb;
            InvWeightFct(nb) = CoefLegendre(i)*CoefOddJacobi(i, j)*CoefEvenJacobi(i+j, k);
            nb++;
          }            
    
    nb_dof_loc = InvWeightFct.GetM();
    VectR3 points_dof3d;
    points_dof3d.Reallocate(nb_points_quadrature_inside);
    for (int i = 0; i < points_dof3d.GetM(); i++)
      points_dof3d(i) = this->PointsND(i);
    
    this->SetPointsDofND(points_dof3d);
    this->nb_points_dof_inside = points_dof3d.GetM();
    
    // retrieving points1d_z and points1d_y
    points1d_z.Reallocate(order+1);
    points1d_y.Reallocate(order+1);
    for (int i = 0; i <= order; i++)
      {
        points1d_z(i) = this->PointsND(i)(2);
        points1d_y(i) = this->PointsND(i*(order+1))(1)/(1.0-points1d_z(0));
      }
    
    TriangleDgOrtho* tri = new TriangleDgOrtho();
    tri->ConstructFiniteElement(order);
    element_tri_surf = tri;
  }
  

  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  template<class Vector1, class Vector2>
  void TetrahedronDgOrtho::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    Vector1 feval_weight = feval;
    const VectReal_wp& weights3d = this->WeightsND();
    for (int i = 0; i < feval_weight.GetM(); i++)
      feval_weight(i) *= weights3d(i);
    
    ApplyCh(feval_weight, contrib);
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
  void TetrahedronDgOrtho::ApplyChGen(const Vector1& U, Vector2& V) const
  {
    //Vector1 Ux(ChX.GetM()), Uy(ChY.GetM());
    typedef typename Vector1::value_type Complexe;
    Vector1& Ux = TetrahedronDgOrthoVariables<Complexe>::Ux(order);
    Vector1& Uy = TetrahedronDgOrthoVariables<Complexe>::Uy(order);
    Ux.Reallocate(ChX.GetM());
    Uy.Reallocate(ChY.GetM());

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
  void TetrahedronDgOrtho::ApplyChTransposeGen(const Vector1& U, Vector2& V) const
  {
    //Vector1 Uy(ChY.GetN()), Uz(ChZ.GetN());
    typedef typename Vector1::value_type Complexe;
    Vector1& Uy2 = TetrahedronDgOrthoVariables<Complexe>::Uy2(order);
    Vector1& Uz = TetrahedronDgOrthoVariables<Complexe>::Uz(order);
    
    Uy2.Reallocate(ChY.GetN());
    Uz.Reallocate(ChZ.GetN());
    
    Mlt(SeldonTrans, ChZ, U, Uz);
    Mlt(SeldonTrans, ChY, Uz, Uy2);
    Mlt(SeldonTrans, ChX, Uy2, V);    
  }
  
  
  //! Integration against gradient of basis functions associated with quadrature points
  /*!
    \param[in] U vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] V Vh_i = \int_K f \nabla \psi_i dx
    This operation is equivalent to a matrix vector product
    V = Rh U
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void TetrahedronDgOrtho::ApplyRhQuadratureGen(const Vector1& U, Vector2& V) const
  {
    Vector1 dU(U.GetM());
    // applying chaine rule to get derivatives on unit tetrahedron
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	dU(3*i) = DerivDxtildeDx(i)*U(3*i)
          + DerivDxtildeDy(i)*U(3*i+1) + DerivDxtildeDy(i)*U(3*i+2);
	
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
  void TetrahedronDgOrtho::ApplyRhQuadratureTransposeGen(const Vector1& U, Vector2& V) const
  {    
    Mlt(SeldonTrans, RhLoc, U, V);
    
    // applying chaine rule to get derivatives on unit tetrahedron
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	V(3*i+2) += V(3*i+1)*DerivDytildeDz(i) + V(3*i)*DerivDxtildeDy(i);
	V(3*i+1) = V(3*i+1)*DerivDytildeDy(i) + V(3*i)*DerivDxtildeDy(i);
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
  void TetrahedronDgOrtho::ApplyShQuadratureGen(const T0& alpha, int num_loc,
						const Vector1& U, Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    //Vector1 Uquad(Nquad);
    typedef typename Vector1::value_type Complexe;
    Vector1& Uquad = TetrahedronDgOrthoVariables<Complexe>::Uquad(order);
    Uquad.Reallocate(Nquad);
    
    if ( (r == 0) || (r == order))
      {
	Mlt(ShSurf(num_loc), U, Uquad);	
      }
    else
      {
	Mlt(ProjOperatorSh(r)(num_loc), U, Uquad);
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
  void TetrahedronDgOrtho::ApplyShQuadratureTransposeGen(int num_loc, const Vector1& U,
							 Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    //Vector1 Uquad(Nquad);
    typedef typename Vector1::value_type Complexe;
    Vector1& Uquad = TetrahedronDgOrthoVariables<Complexe>::Uquad(order);
    Uquad.Reallocate(Nquad);
    
    Mlt(SeldonTrans, ShLoc(num_loc), U, Uquad);
    
    if ( (r == 0) || (r == order))
      {	
	Mlt(SeldonTrans, ShSurf(num_loc), Uquad, V);
      }
    else
      {
	Mlt(SeldonTrans, ProjOperatorSh(r)(num_loc), Uquad, V);
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
  void TetrahedronDgOrtho::ApplyNablaShQuadratureGen(const T0& alpha, int num_loc,
						     const Vector1& Uh, Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    //Vector1 Uquad(Nquad);
    //typedef typename Vector1::value_type Complexe;
    Vector1 Uquad, U = Uh;
    Uquad.Reallocate(2*Nquad);
    
    
    if ( (r == 0) || (r == order))
      {
	int Ntri = this->Points2D_tri().GetM();
        int offset = nb_points_quadrature_inside + num_loc*Ntri;
        for (int i = 0; i < Ntri; i++)
          {
            U(3*i) = DerivDxtildeDx(offset+i)*U(3*i) + DerivDxtildeDy(offset+i)*U(3*i+1)
              + DerivDxtildeDy(offset+i)*U(3*i+2);
            
            U(3*i+1) = U(3*i+1)*DerivDytildeDy(offset+i) + U(3*i+2)*DerivDytildeDz(offset+i);
          }
        
	Mlt(DerivShSurf(num_loc), U, Uquad);	
      }
    else
      {
        const Vector<TinyVector<Real_wp, 4> >& Dx = ProjOperatorDxtilde(r)(num_loc);
        for (int i = 0; i < Dx.GetM(); i++)
          {
            U(3*i) = Dx(i)(0)*U(3*i) + Dx(i)(1)*U(3*i+1) + Dx(i)(1)*U(3*i+2);
            U(3*i+1) = Dx(i)(2)*U(3*i+1) + Dx(i)(3)*U(3*i+2);
          }
        
	Mlt(ProjOperatorDerivSh(r)(num_loc), U, Uquad);
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
  void TetrahedronDgOrtho
  ::ApplyNablaShQuadratureTransposeGen(int num_loc, const Vector1& U, Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    //Vector1 Uquad(Nquad);
    //typedef typename Vector1::value_type Complexe;
    Vector1 Uquad;
    Uquad.Reallocate(2*Nquad);
    
    Mlt(SeldonTrans, DerivShLoc(num_loc), U, Uquad);
    
    if ( (r == 0) || (r == order))
      {	
	Mlt(SeldonTrans, DerivShSurf(num_loc), Uquad, V);
        
	int Ntri = this->Points2D_tri().GetM();
        int offset = nb_points_quadrature_inside + num_loc*Ntri;
        for (int i = 0; i < Ntri; i++)
          {
            V(3*i+2) += V(3*i+1)*DerivDytildeDz(offset+i) + V(3*i)*DerivDxtildeDy(offset+i);
            V(3*i+1) = V(3*i+1)*DerivDytildeDy(offset+i) + V(3*i)*DerivDxtildeDy(offset+i);
            V(3*i) *= DerivDxtildeDx(offset+i);
          }
        
      }
    else
      {
	Mlt(SeldonTrans, ProjOperatorDerivSh(r)(num_loc), Uquad, V);
        
        const Vector<TinyVector<Real_wp, 4> >& Dx = ProjOperatorDxtilde(r)(num_loc);
        for (int i = 0; i < Dx.GetM(); i++)
          {
            V(3*i+2) += Dx(i)(3)*V(3*i+1) + Dx(i)(1)*V(3*i);
            V(3*i+1) = Dx(i)(2)*V(3*i+1) + Dx(i)(1)*V(3*i);
            V(3*i) *= Dx(i)(0);
          }
        
      }
  }
  
  
  //! computation of projection between finite element of different orders    
  void TetrahedronDgOrtho::
  ComputeTriangularInterpolationProjector(const IVect& order_elt,
                                          const ElementReference<Dimension2, 1>& Fb,
                                          const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi;
    ProjOperatorSh.Reallocate(rmax+1);
    ProjOperatorDxtilde.Reallocate(rmax+1);
    ProjOperatorDerivSh.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
        int Nquad = (order+1)*(order+1);
        int Ntri = Pts(i).GetM();
        // face z = 0
	ProjOperatorSh(r)(0).Reallocate(Nquad, Ntri);
        ProjOperatorDxtilde(r)(0).Reallocate(Ntri);
        Matrix<Real_wp, General, ArrayRowSparse> sh(2*Nquad, 3*Ntri);
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Real_wp yt = Pts(i)(j)(1), xt = Pts(i)(j)(0)/(1.0-yt);
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = i1*(order+1) + i2;		  
		  ProjOperatorSh(r)(0)(k, j) = gauss_x.EvaluatePhi(i1, xt)*gauss_y.EvaluatePhi(i2, yt);
                  sh.AddInteraction(2*k, 3*j, gauss_x.EvaluatePhiGrad(i1, xt)*gauss_y.EvaluatePhi(i2, yt));
                  sh.AddInteraction(2*k, 3*j+1, gauss_x.EvaluatePhi(i1, xt)*gauss_y.EvaluatePhiGrad(i2, yt));
                  sh.AddInteraction(2*k+1, 3*j+2, gauss_x.EvaluatePhi(i1, xt)*gauss_y.EvaluatePhi(i2, yt));
		}
            
            Real_wp x = Pts(i)(j)(0), y = Pts(i)(j)(1), z = 0;
            ProjOperatorDxtilde(r)(0)(j)(0) = 1.0/(1.0-y-z);
            ProjOperatorDxtilde(r)(0)(j)(1) = x/square(1.0-y-z);
            ProjOperatorDxtilde(r)(0)(j)(2) = 1.0/(1.0-z);
            ProjOperatorDxtilde(r)(0)(j)(3) = y/square(1.0-z);
	  }
        
        Copy(sh, ProjOperatorDerivSh(r)(0));
        
        // face y = 0
	ProjOperatorSh(r)(1).Reallocate((order+1)*(order+1), Pts(i).GetM());
        ProjOperatorDxtilde(r)(1).Reallocate(Ntri);
        sh.Clear(); sh.Reallocate(2*Nquad, 3*Ntri);
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Real_wp zt = Pts(i)(j)(1), xt = Pts(i)(j)(0)/(1.0-zt);
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = i1*(order+1) + i2;
		  ProjOperatorSh(r)(1)(k, j) = gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhi(i2, zt);
                  sh.AddInteraction(2*k, 3*j, gauss_x.EvaluatePhiGrad(i1, xt)*gauss_z.EvaluatePhi(i2, zt));
                  sh.AddInteraction(2*k, 3*j+2, gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhiGrad(i2, zt));
                  sh.AddInteraction(2*k+1, 3*j+1, gauss_x.EvaluatePhi(i1, xt)*gauss_z.EvaluatePhi(i2, zt));
		}
            
            Real_wp x = Pts(i)(j)(0), y = 0, z = Pts(i)(j)(1);
            ProjOperatorDxtilde(r)(1)(j)(0) = 1.0/(1.0-y-z);
            ProjOperatorDxtilde(r)(1)(j)(1) = x/square(1.0-y-z);
            ProjOperatorDxtilde(r)(1)(j)(2) = 1.0/(1.0-z);
            ProjOperatorDxtilde(r)(1)(j)(3) = y/square(1.0-z);
	  }
        
        Copy(sh, ProjOperatorDerivSh(r)(1));

	ProjOperatorSh(r)(2).Reallocate((order+1)*(order+1), Pts(i).GetM());
        ProjOperatorDxtilde(r)(2).Reallocate(Ntri);
        sh.Clear(); sh.Reallocate(2*Nquad, 3*Ntri);
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Real_wp zt = Pts(i)(j)(1), yt = Pts(i)(j)(0)/(1.0-zt);
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = i1*(order+1) + i2;
		  ProjOperatorSh(r)(2)(k, j) = gauss_y.EvaluatePhi(i1, yt)*gauss_z.EvaluatePhi(i2, zt);
                  sh.AddInteraction(2*k, 3*j+1, gauss_y.EvaluatePhiGrad(i1, yt)*gauss_z.EvaluatePhi(i2, zt));
                  sh.AddInteraction(2*k, 3*j+2, gauss_y.EvaluatePhi(i1, yt)*gauss_z.EvaluatePhiGrad(i2, zt));
                  sh.AddInteraction(2*k+1, 3*j, gauss_y.EvaluatePhi(i1, yt)*gauss_z.EvaluatePhi(i2, zt));
		}
            
            Real_wp x = 0, y = Pts(i)(j)(0), z = Pts(i)(j)(1);
            ProjOperatorDxtilde(r)(2)(j)(0) = 1.0/(1.0-y-z);
            ProjOperatorDxtilde(r)(2)(j)(1) = x/square(1.0-y-z);
            ProjOperatorDxtilde(r)(2)(j)(2) = 1.0/(1.0-z);
            ProjOperatorDxtilde(r)(2)(j)(3) = y/square(1.0-z);
	  }
        
        Copy(sh, ProjOperatorDerivSh(r)(2));
        
	ProjOperatorSh(r)(3).Reallocate((order+1)*(order+1), Pts(i).GetM());
        ProjOperatorDxtilde(r)(3).Reallocate(Ntri);
        sh.Clear(); sh.Reallocate(2*Nquad, 3*Ntri);
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Real_wp zt = Pts(i)(j)(1), yt = Pts(i)(j)(0)/(1.0-zt);
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = i1*(order+1) + i2;		  
		  ProjOperatorSh(r)(3)(k, j) = gauss_y.EvaluatePhi(i1, yt)*gauss_z.EvaluatePhi(i2, zt);
                  sh.AddInteraction(2*k, 3*j+1, gauss_y.EvaluatePhiGrad(i1, yt)*gauss_z.EvaluatePhi(i2, zt));
                  sh.AddInteraction(2*k, 3*j+2, gauss_y.EvaluatePhi(i1, yt)*gauss_z.EvaluatePhiGrad(i2, zt));
                  sh.AddInteraction(2*k+1, 3*j, gauss_y.EvaluatePhi(i1, yt)*gauss_z.EvaluatePhi(i2, zt));
		}
            
            Real_wp y = Pts(i)(j)(0), z = Pts(i)(j)(1), x = 1.0-y-z;
            ProjOperatorDxtilde(r)(3)(j)(0) = 1.0/(1.0-y-z);
            ProjOperatorDxtilde(r)(3)(j)(1) = x/square(1.0-y-z);
            ProjOperatorDxtilde(r)(3)(j)(2) = 1.0/(1.0-z);
            ProjOperatorDxtilde(r)(3)(j)(3) = y/square(1.0-z);
	  }
        
        Copy(sh, ProjOperatorDerivSh(r)(3));
      }
  }
      
}

#define MONTJOIE_FILE_TETRAHEDRON_DG_ORTHO_CXX
#endif
