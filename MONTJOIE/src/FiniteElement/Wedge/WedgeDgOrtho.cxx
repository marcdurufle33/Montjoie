#ifndef MONTJOIE_FILE_WEDGE_DG_ORTHO_CXX

namespace Montjoie
{

  //! default constructor
  WedgeDgOrtho::WedgeDgOrtho() : WedgeReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;

    this->discontinuous_element = true;
    this->optimized_mass_matrix = true;
    this->linear_sparse_mass_matrix = true;
    this->sparse_mass_matrix = true;
    
    use_quadrature_for_sh = true;
    use_quadrature_for_rh = true;
  }


  //! how to number mesh
  void WedgeDgOrtho::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg != ElementReference_Base::HDG)
      return WedgeReference<1>::ConstructNumberMap(nmap, dg);
    
    nmap.SetFormulationDG(dg);
    
    nmap.SetNbDofVertex(this->order, 0);
    nmap.SetNbDofEdge(this->order, 0);
    nmap.SetNbDofTriangle(this->order, (this->order+2)*(this->order+1)/2);
    nmap.SetNbDofQuadrangle(this->order, (this->order+1)*(this->order+1));
    nmap.SetNbDofPyramid(this->order, 0);
    
    // rotation on triangular faces
    ElementReference<Dimension2, 1>::
      FindH1RotationTri(this->order, this->Points2D_tri(), this->Weights2D_tri(),
                        element_tri_surf->GetValuePhi(), nmap);
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    Matrix<int> FacesDof_Rotation_Quad;
    MeshNumbering<Dimension3>::
      GetRotationQuadrilateralFace(NumQuad2D, FacesDof_Rotation_Quad);
    
    nmap.SetFacesDofRotationQuad(this->order, FacesDof_Rotation_Quad);
  }


  size_t WedgeDgOrtho::GetMemorySize() const
  {
    size_t taille = WedgeReference<1>::GetMemorySize();
    
    for (int i = 0; i < ShLoc.GetM(); i++)
      {
        taille += ShLoc(i).GetMemorySize();
        taille += DerivShLoc(i).GetMemorySize();
      }

    taille += ChX.GetMemorySize() + ChY.GetMemorySize();
    taille += RhLoc.GetMemorySize() + ShQuad.GetMemorySize();
    taille += NumQuad3D.GetMemorySize() + NumOrtho3D.GetMemorySize(); 
    taille += points1d_y.GetMemorySize() + weights1d_y.GetMemorySize();
    taille += DerivDxtildeDx.GetMemorySize() + DerivDxtildeDy.GetMemorySize();
    taille += InvWeight2D.GetMemorySize() + ChTri.GetMemorySize();
    taille += NumOrtho2D.GetMemorySize() + NumDofTri.GetMemorySize() 
      + CoordinateDofTri.GetMemorySize() + gauss.GetMemorySize();;
    taille += function_tri.GetMemorySize();
    taille += gauss_x.GetMemorySize() + gauss_y.GetMemorySize();
    
    taille += ShTri.GetMemorySize();
    taille += DerivShTri.GetMemorySize() + DerivShQuad.GetMemorySize()
      + DerivShGauss.GetMemorySize();
    
    taille += Seldon::GetMemorySize(ProjOperatorQuadOrderY);
    
    taille += Seldon::GetMemorySize(ProjOperatorDxtildeTri);
    taille += Seldon::GetMemorySize(ProjOperatorDxtildeQuad);
    
    taille += Seldon::GetMemorySize(ProjOperatorDerivQuadOrder);
    taille += Seldon::GetMemorySize(ProjOperatorDerivTriOrder);
    taille += Seldon::GetMemorySize(ProjOperatorDerivQuadOrderY);

    taille += alpha.GetMemorySize() + beta.GetMemorySize() + gamma.GetMemorySize();
    taille += C.GetMemorySize() + CmO.GetMemorySize() + CpO.GetMemorySize();
    return taille;
  }

    
  //! constructing finite element
  void WedgeDgOrtho::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					    int rsurf_tri, int rsurf_quad,
					    int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    WedgeReference<1>::
      ConstructFiniteElement(r, rgeom, rquad, TriangleQuadrature::QUADRATURE_TENSOR);
        
    ConstructFunctions();    
    
    ConstructMassMatrix();
    ConstructStiffnessMatrix();
    
    this->Fb_geom.ComputeCoefficientTransformation();
    
    this->ConstructElementaryMatrix(*this);
  }
  
        
  //! construction of mass matrix
  void WedgeDgOrtho::ConstructMassMatrix()
  {
    // computation of ChX, ChY
        
    // first interpolation along y :
    // v^y_i,j,k = \sum_m P_m^(2i+1)(xi_j) ((1-xi_j)/2)^i v^z_{i,m,k}
    Array3D<Real_wp> Pv(order+1, order+1, order+1); Pv.Fill(0);
    VectReal_wp Pm;
    const Vector<Matrix<Real_wp> >& tri_OddJacobiPolynom = function_tri.GetOddJacobiPolynomial();
    const Matrix<Real_wp>& tri_CoefOddJacobi = function_tri.GetCoefficientOddJacobi();
    for (int k = 0; k <= order; k++)
      {
	Real_wp xi_k = 2.0*points1d_y(k) - 1.0;
	for (int i = 0; i <= order; i++)
	  {
	    EvaluateJacobiPolynomial(tri_OddJacobiPolynom(i), order-i, xi_k, Pm);
	    for (int j = 0; j < Pm.GetM(); j++)
	      Pv(i, j, k) = Pm(j)*tri_CoefOddJacobi(i, j)*pow(0.5*(1.0-xi_k), i);
	  }
      }
    
    int nb_points_y = (order+1)*(order+1)*(order+1);
    Matrix<Real_wp, General, ArrayRowSparse> ch;
    ch.Clear();
    ch.Reallocate(nb_dof_loc, nb_points_y);
    int rp1 = order+1;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = rp1*(i*rp1 + j) + k;
	    for (int m = 0; m <= order-i; m++)
              ch.AddInteraction(NumOrtho3D(i, m, k), node, Pv(i, m, j));
            
	  }
    
    // conversion to csr
    Copy(ch, ChY);

    // third interpolation along x :
    // v_i,j,k = \sum_m L_m(xi_i) v^y_{m,j,k}
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
    ch.Reallocate(nb_points_y, nb_points_y);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    for (int m = 0; m <= order; m++)
              ch.AddInteraction(rp1*(rp1*m + j) + k, rp1*(rp1*i + j) + k, Pv(0, m, i));            
	  }
    
    // conversion to csr
    Copy(ch, ChX);
    
    // useful coefficients
    alpha.Reallocate(order+1, 2*order+2); alpha.Fill(0);
    beta.Reallocate(order+1, 2*order+2); beta.Fill(0);
    gamma.Reallocate(order+1, 2*order+2); gamma.Fill(0);

    const Matrix<Real_wp>& tri_LegendrePolynom = function_tri.GetLegendrePolynomial();
    const VectReal_wp& tri_CoefLegendre = function_tri.GetCoefficientLegendre();
    for (int i = 0; i <= order; i++)
      {
	alpha(0, 0) = tri_LegendrePolynom(0,0);
	alpha(0, 2*i+1) = tri_OddJacobiPolynom(i)(0,0);
	
	gamma(0, 0) = tri_CoefLegendre(0)/tri_CoefLegendre(1);
	gamma(0, 2*i+1) = tri_CoefOddJacobi(i,0)/tri_CoefOddJacobi(i,1);
	
	for (int n = 1; n <= order-i; n++)
	  {     
	    alpha(n, 0) = tri_LegendrePolynom(n,0);
	    alpha(n, 2*i+1) = tri_OddJacobiPolynom(i)(n,0);
	    
	    beta(n, 0) = tri_LegendrePolynom(n,1)
              *tri_CoefLegendre(n)/tri_CoefLegendre(n-1);
	    beta(n, 2*i+1) = tri_OddJacobiPolynom(i)(n,1)
              *tri_CoefOddJacobi(i,n)/tri_CoefOddJacobi(i,n-1);
	    
	    gamma(n, 0) = tri_CoefLegendre(n)/tri_CoefLegendre(n+1);
	    gamma(n, 2*i+1) = tri_CoefOddJacobi(i,n)/tri_CoefOddJacobi(i,n+1);
	  }
      }
    
    VectReal_wp Pn(order+2),PnmO(order+2),PnpO(order+2);
    C.Reallocate(order+1, order+1, order+1); C.Fill(0);
    CpO.Reallocate(order+1, order+1, order+1); CpO.Fill(0);
    CmO.Reallocate(order+1, order+1, order+1); CmO.Fill(0);

    VectReal_wp xi, omega;
    ComputeGaussLegendre(xi, omega, order+1);
    for (int k = 0; k <= order+1; k++)
      for (int i = 0; i <= order; i++)
	{
	  EvaluateJacobiPolynomial(tri_OddJacobiPolynom(i), order-i, 2*xi(k)-1, Pn);
	  EvaluateJacobiPolynomial(tri_OddJacobiPolynom(i+1), order-1-i, 2*xi(k)-1, PnpO);
	  if (i==0)
	    PnmO = 0;
	  else 
	    {
	      EvaluateJacobiPolynomial(tri_OddJacobiPolynom(i-1),
                                       order+1-i, 2*xi(k)-1, PnmO);
	      for (int j = 0; j < PnmO.GetM(); j++)
		PnmO(j) *= tri_CoefOddJacobi(i-1,j);
	    }
	  
	  for (int j = 0; j < Pn.GetM(); j++)
	    Pn(j) *= tri_CoefOddJacobi(i,j);
	  
	  for (int j = 0; j < PnpO.GetM(); j++)
	    PnpO(j) *= tri_CoefOddJacobi(i+1,j);
	  
	  for (int j1 = 0; j1 <= order-i; j1++)
	    {
	      for (int j2 = 0; j2 <= order-i; j2++)
		C(i,j1,j2) += omega(k)*pow(1-xi(k),2*i+2)*Pn(j1)*Pn(j2);
	      
	      for (int j2 = 0; j2 <= order-1-i; j2++)
		CpO(i,j1,j2) += omega(k)*pow(1-xi(k),2*i+3)*Pn(j1)*PnpO(j2);
	      
	      if (i > 0)
		for (int j2 = 0; j2 <= order+1-i; j2++)
		  CmO(i,j1,j2) += omega(k)*pow(1-xi(k),2*i+1)*Pn(j1)*PnmO(j2);
	    }
	}
  }
 

  //! construction of stiffness matrix
  void WedgeDgOrtho::ConstructStiffnessMatrix()
  {
    const VectReal_wp& points1d = this->Points1D();
    const VectR3& points3d = this->PointsND();

    // Computation of gradient from values on quadrature points
    gauss_x.AffectPoints(points1d);
    gauss_y.AffectPoints(points1d_y);
    gauss_x.ComputeGradPhi(1e3*epsilon_machine);
    gauss_y.ComputeGradPhi(1e3*epsilon_machine);
    
    int N = (order+1)*(order+1)*(order+1);
    Matrix<Real_wp, General, ArrayRowSparse> rh;
    rh.Reallocate(N, 3*N);
    DerivDxtildeDx.Reallocate(points3d.GetM());
    DerivDxtildeDy.Reallocate(points3d.GetM());
    int rp1 = order+1;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = rp1*(rp1*i + j) + k;
	    for (int m = 0; m <= order; m++)
	      {
		int node2 = 3*(rp1*(rp1*m + j) + k);
                rh.AddInteraction(node, node2, gauss_x.GradPhi(i, m));

		node2 = 3*(rp1*(rp1*i + m) + k) + 1;
                rh.AddInteraction(node, node2, gauss_y.GradPhi(j, m));

		node2 = 3*(rp1*(rp1*i + j) + m) + 2;
                rh.AddInteraction(node, node2, gauss_x.GradPhi(k, m));
	      }
	    
	    Real_wp x = points3d(node)(0), y = points3d(node)(1);
	    DerivDxtildeDx(node) = 1.0/(1.0-y);
	    DerivDxtildeDy(node) = x/square(1.0-y);
	  }
    
    Copy(rh, RhLoc);
    
    // construction of ShLoc
    int Nquad = (order+1)*(order+1);
    rh.Clear();
    Matrix<Real_wp, General, ArrayRowSparse> sh0, sh1, sh2, sh3, sh4;
    sh0.Reallocate(N, Nquad); sh1.Reallocate(N, Nquad); sh2.Reallocate(N, Nquad);
    sh3.Reallocate(N, Nquad); sh4.Reallocate(N, Nquad);
    VectReal_wp ValPhiX0(order+1), ValPhiX1(order+1), ValPhiY0(order+1);
    VectReal_wp DerivPhiX0(order+1), DerivPhiX1(order+1), DerivPhiY0(order+1);
    
    for (int i = 0; i <= order; i++)
      {
	ValPhiX0(i) = gauss_x.EvaluatePhi(i, Real_wp(0));
	ValPhiX1(i) = gauss_x.EvaluatePhi(i, Real_wp(1));
	ValPhiY0(i) = gauss_y.EvaluatePhi(i, Real_wp(0));

	DerivPhiX0(i) = gauss_x.EvaluatePhiGrad(i, Real_wp(0));
	DerivPhiX1(i) = gauss_x.EvaluatePhiGrad(i, Real_wp(1));
	DerivPhiY0(i) = gauss_y.EvaluatePhiGrad(i, Real_wp(0));
      }
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int nv = rp1*(rp1*i +j) + k;
	    int ns = NumQuad2D(i, j);
	    sh0.AddInteraction(nv, ns, ValPhiX0(k));
            
	    ns = NumQuad2D(i, k);
            sh1.AddInteraction(nv, ns, ValPhiY0(j));
	    
	    ns = NumQuad2D(j, k);
            sh2.AddInteraction(nv, ns, ValPhiX1(i));
	    
	    ns = NumQuad2D(j, k);
            sh3.AddInteraction(nv, ns, ValPhiX0(i));

	    ns = NumQuad2D(i, j);
            sh4.AddInteraction(nv, ns, ValPhiX1(k));
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
            sh0.AddInteraction(nv, 2*ns, ValPhiX0(k));
            sh0.AddInteraction(nv, 2*ns+1, DerivPhiX0(k));
            
	    ns = NumQuad2D(i, k);
            sh1.AddInteraction(nv, 2*ns, ValPhiY0(j));
            sh1.AddInteraction(nv, 2*ns+1, DerivPhiY0(j));
	    
	    ns = NumQuad2D(j, k);
            sh2.AddInteraction(nv, 2*ns, ValPhiX1(i));
            sh2.AddInteraction(nv, 2*ns+1, DerivPhiX1(i));
	    
	    ns = NumQuad2D(j, k);
            sh3.AddInteraction(nv, 2*ns, ValPhiX0(i));
            sh3.AddInteraction(nv, 2*ns+1, DerivPhiX0(i));
            
            ns = NumQuad2D(i, j);
            sh4.AddInteraction(nv, 2*ns, ValPhiX1(k));
            sh4.AddInteraction(nv, 2*ns+1, DerivPhiX1(k));            
	  }

    Copy(sh0, DerivShLoc(0)); Copy(sh1, DerivShLoc(1));
    Copy(sh2, DerivShLoc(2)); Copy(sh3, DerivShLoc(3)); Copy(sh4, DerivShLoc(4));
    
    // computing ShQuad, ShQuadTranspose
    rh.Reallocate(Nquad, Nquad);
    Matrix<Real_wp> ValPhiY(order+1, order+1);
    Matrix<Real_wp> DerivPhiY(order+1, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
          ValPhiY(i, j) = gauss_y.EvaluatePhi(i, gauss_x.Points(j));
          DerivPhiY(i, j) = gauss_y.EvaluatePhiGrad(i, gauss_x.Points(j));
        }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  int ns = NumQuad2D(i, j);
	  for (int m = 0; m <= order; m++)
            rh.AddInteraction(NumQuad2D(m, j), ns, ValPhiY(m, i));
	}
    
    Copy(rh, ShQuad);
    rh.Clear();
    
    // now computing ShSurf and ShSurfTranspose
    int Ntri = this->Points2D_tri().GetM();
    ShTri.Reallocate(Nquad, Ntri); 
    sh0.Clear(); sh1.Clear(); sh2.Clear();
    sh0.Reallocate(2*Nquad, 3*Ntri); sh1.Reallocate(2*Nquad, 3*Nquad); 
    sh2.Reallocate(2*Nquad, 3*Nquad);
         
    for (int i = 0; i < Ntri; i++)
      {
	int m = nb_points_quadrature_inside + i;
	Real_wp x = points3d(m)(0), y = points3d(m)(1);
	Real_wp xt = x/(1.0-y), yt = y;
	
	// face z = 0
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_x.EvaluatePhi(j, xt)*gauss_y.EvaluatePhi(k, yt);
	      int col = NumQuad2D(j, k);
	      ShTri(col, i) = vloc;
              sh0.AddInteraction(2*col, 3*i, gauss_x.EvaluatePhiGrad(j, xt)*gauss_y.EvaluatePhi(k, yt));
              sh0.AddInteraction(2*col, 3*i+1, gauss_x.EvaluatePhi(j, xt)*gauss_y.EvaluatePhiGrad(k, yt));
              sh0.AddInteraction(2*col+1, 3*i+2, vloc);
	    }	
        
        DerivDxtildeDx(m) = 1.0/(1.0-y);
        DerivDxtildeDy(m) = x/square(1.0-y);

        m += Ntri + 3*this->Points2D_quad().GetM();
        DerivDxtildeDx(m) = 1.0/(1.0-y);
        DerivDxtildeDy(m) = x/square(1.0-y);
        
      }    

    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
        {
          int i = NumQuad2D(i1, i2);
          for (int j1 = 0; j1 <= order; j1++)
            for (int j2 = 0; j2 <= order; j2++)
              {
                Real_wp delta = 0;
                if (i2 == j2)
                  delta = 1.0;
                
                int col = NumQuad2D(j1, j2);
                sh1.AddInteraction(2*col+1, 3*i, ValPhiY(j1, i1)*delta);
                sh1.AddInteraction(2*col, 3*i+1, DerivPhiY(j1, i1)*delta);
                sh1.AddInteraction(2*col, 3*i+2, ValPhiY(j1, i1)*gauss_x.GradPhi(j2, i2));
              }
          
          int col = NumQuad2D(i1, i2);
          sh2.AddInteraction(2*col+1, 3*i+1, 1.0);
          for (int j1 = 0; j1 <= order; j1++)
            {
              col = NumQuad2D(j1, i2);
              sh2.AddInteraction(2*col, 3*i, gauss_x.GradPhi(j1, i1));

              col = NumQuad2D(i1, j1);
              sh2.AddInteraction(2*col, 3*i+2, gauss_x.GradPhi(j1, i2));
            }

          int m = nb_points_quadrature_inside + Ntri + i;
          Real_wp x = gauss_x.Points(i1), y = 0.0;
          DerivDxtildeDx(m) = 1.0/(1.0-y);
          DerivDxtildeDy(m) = x/square(1.0-y);
          
          m = nb_points_quadrature_inside + Ntri + Nquad + i;
          y = gauss_x.Points(i1); x = 1.0 - y;
          DerivDxtildeDx(m) = 1.0/(1.0-y);
          DerivDxtildeDy(m) = x/square(1.0-y);

          m = nb_points_quadrature_inside + Ntri + 2*Nquad + i;
          y = gauss_x.Points(i1); x = 0;
          DerivDxtildeDx(m) = 1.0/(1.0-y);
          DerivDxtildeDy(m) = x/square(1.0-y);
          
        }
    
    Copy(sh0, DerivShTri); Copy(sh1, DerivShQuad); Copy(sh2, DerivShGauss);    
  }


  //! construction of basis functions
  void WedgeDgOrtho::ConstructFunctions()
  {    
    function_tri.ConstructFiniteElement(order);
    
    VectR3 points_dof3d;
    points_dof3d.Reallocate(nb_points_quadrature_inside);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      points_dof3d(i) = this->PointsND(i);
    
    this->SetPointsDofND(points_dof3d);
    
    this->nb_points_dof_inside = nb_points_quadrature_inside;
    nb_dof_loc = (order+1)*(order+2)*(order+1)/2;
    nb_dof_boundaries = 0;
    
    gauss.ConstructQuadrature(order, gauss.QUADRATURE_GAUSS);
    NumOrtho2D = function_tri.GetNumOrtho2D();
    InvWeight2D = function_tri.GetInverseWeightPolynomial();
    
    NumDofTri.Reallocate((order+1)*(order+2)/2, order+1);
    NumOrtho3D.Reallocate(order+1, order+1, order+1);
    
    int node = 0;
    CoordinateDofTri.Reallocate(nb_dof_loc, 2);
    for (int k = 0; k <= order; k++)
      for (int i = 0; i <= order; i++)
        for (int j = 0; j <= order-i; j++)
          {
            NumDofTri(NumOrtho2D(i, j), k) = node;
            CoordinateDofTri(node, 0) = NumOrtho2D(i, j);
            CoordinateDofTri(node, 1) = k;
            NumOrtho3D(i, j, k) = node++;
          }

    VectReal_wp phi; node = 0; R2 pointA; R3 pt3D;
    ChTri.Reallocate((order+1)*(order+2)/2, (order+1)*(order+1));
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        {
          pt3D = this->PointsND((order+1)*(i*(order+1) + j));
          pointA.Init(pt3D(0), pt3D(1));
          function_tri.ComputeValuesPhiOrthoRef(order, NumOrtho2D, InvWeight2D, pointA, phi);
          for (int k = 0; k < phi.GetM(); k++)
	    ChTri(k, node) = phi(k);
          
          node++;
        }        
    
    ComputeGaussJacobi(points1d_y, weights1d_y, order, Real_wp(1), Real_wp(0));
    
    TriangleDgOrtho* tri = new TriangleDgOrtho();
    tri->ConstructFiniteElement(order);
    element_tri_surf = tri;

    QuadrangleDgGauss* quad = new QuadrangleDgGauss();
    quad->ConstructFiniteElement(order);
    element_quad_surf = quad;
  }
  
  
  //! computation of mass matrix knowing decomposition of Jacobian for non-affine elements
  /*!
    \param[out] Mm computed mass matrix
    \param[in] CoefJacobian decomposition of jacobian |DF|
   */
  template<class Matrix1>
  void WedgeDgOrtho
  ::PerformComputationMassMatrix(Matrix1 & Mm, const VectReal_wp & CoefJacobian) const
  {
    Mm.Reallocate(nb_dof_loc,nb_dof_loc); Mm.Fill(0);
    Real_wp A = CoefJacobian(0);
    Real_wp B1 = CoefJacobian(1), B2 = CoefJacobian(2), B3 = CoefJacobian(3);
    Real_wp C1 = CoefJacobian(4), C2 = CoefJacobian(5);
    Real_wp D = CoefJacobian(6);
    
    for (int i3 = 0; i3 <= order; i3++)
      {
        Real_wp poids = gauss.Weights(i3);
        Real_wp z = gauss.Points(i3);
        for (int i1 = 0; i1 <= order; i1++)
          for (int i2 = 0; i2 <= order-i1; i2++)
            {
              int node = NumOrtho2D(i1, i2);
              int i = NumDofTri(node, i3);
              // diagonal coefficient
              Mm.Get(i, i) += (A + 0.5*B2*(alpha(i2, 2*i1+1)+1.0)
                               + z*(B3 + 0.5*C2*(alpha(i2, 2*i1+1)+1.0) + z*D))*poids;
              if (i2 < order-i1)
                {
                  int j = NumDofTri(NumOrtho2D(i1, i2+1), i3);
                  Mm.Get(i, j) += (0.5*B2 + 0.5*z*C2)*gamma(i2, 2*i1+1)*poids;
                }

              if (i2 > 0)
                {
                  int j = NumDofTri(NumOrtho2D(i1, i2-1), i3);
                  Mm.Get(i, j) += (0.5*B2 + 0.5*z*C2)*beta(i2, 2*i1+1)*poids;
                }

              for (int j2 = 0; j2 <= order+1-i1; j2++)
                {
                  if (i1 <= order-j2)
                    {
                      int j = NumDofTri(NumOrtho2D(i1, j2), i3);
                      Mm.Get(i, j) += (0.5*B1 + 0.5*z*C1)*C(i1, i2, j2)*poids;
                    }
                  
                  if (i1 < order-j2)
                    {
                      int j = NumDofTri(NumOrtho2D(i1+1, j2), i3);
                      Mm.Get(i, j) += (0.5*B1 + 0.5*z*C1)*gamma(i1, 0)*CpO(i1, i2, j2)*poids;
                    }

                  if (i1 > 0)
                    {
                      int j = NumDofTri(NumOrtho2D(i1-1, j2), i3);
                      Mm.Get(i, j) += (0.5*B1 + 0.5*z*C1)*beta(i1, 0)*CmO(i1, i2, j2)*poids;
                    }
                }
            }
      }
  }
  
  
  //! x is overwritten by M^-1 x where M is the mass matrix  
  template<class Vector1>
  void WedgeDgOrtho::SolveMassMatrixGen(Vector1& x) const
  {
    // mass matrix with gauss weights
    const VectReal_wp& weights1d = this->Weights1D();
    for (int k = 0; k <= order; k++)
      {
        Real_wp coef = 1.0/weights1d(k);
        for (int i = 0; i < NumDofTri.GetM(); i++)
          x(NumDofTri(i, k)) *= coef;
      }
  }
  

  //! y is overwritten by M y where M is the mass matrix    
  template<class Vector1>
  void WedgeDgOrtho::MltMassMatrixGen(Vector1& y) const
  {
    // mass matrix with gauss weights
    const VectReal_wp& weights1d = this->Weights1D();
    for (int k = 0; k <= order; k++)
      {
        Real_wp coef = weights1d(k);
        for (int i = 0; i < NumDofTri.GetM(); i++)
          y(NumDofTri(i, k)) *= coef;
      }
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
  void WedgeDgOrtho::ApplyChGen(const Vector1& U, Vector2& V) const
  {
    Vector1 Ux(ChX.GetM());
    Ux.Fill(0);
    Mlt(ChX, U, Ux);
    Mlt(ChY, Ux, V);
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
  void WedgeDgOrtho::ApplyChTransposeGen(const Vector1& U, Vector2& V) const
  {
    Vector1 Uy(ChY.GetN());
    Uy.Fill(0);
    Mlt(SeldonTrans, ChY, U, Uy);
    Mlt(SeldonTrans, ChX, Uy, V);    
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
  void WedgeDgOrtho::ApplyRhQuadratureGen(const Vector1& U, Vector2& V) const
  {
    // applying chaine rule to get derivatives on unit prism
    Vector1 dU(U.GetM());
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	dU(3*i) = DerivDxtildeDx(i)*U(3*i) + DerivDxtildeDy(i)*U(3*i+1);
	dU(3*i+1) = U(3*i+1);
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
  void WedgeDgOrtho::ApplyRhQuadratureTransposeGen(const Vector1& U, Vector2& V) const
  {    
    Mlt(SeldonTrans, RhLoc, U, V);
    
    // applying chaine rule to get derivatives on unit prism
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	V(3*i+1) += V(3*i)*DerivDxtildeDy(i);
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
  void WedgeDgOrtho::ApplyShQuadratureGen(const T0& alpha, int num_loc, const Vector1& U,
					  Vector2& V, int r) const
  {
    int Nquad = ShTri.GetM();
    Vector1 Uquad(Nquad);
    switch (num_loc)
      {
      case 0 :
      case 4 :
        if ((r == 0) || (r == order))
	  Mlt(ShTri, U, Uquad);
	else
	  Mlt(ProjOperatorTriOrder(r), U, Uquad);
	
        break;
      case 1 :
	if ((r == 0) || (r == order))
	  Copy(U, Uquad);
	else
	  Mlt(ProjOperatorQuadOrder(r), U, Uquad);
        
	break;
      case 2 :
      case 3 :
        if ((r == 0) || (r == order))
	  Mlt(ShQuad, U, Uquad);
	else
	  Mlt(ProjOperatorQuadOrderY(r), U, Uquad);
        
	break;
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
  void WedgeDgOrtho::ApplyShQuadratureTransposeGen(int num_loc, const Vector1& U,
						   Vector2& V, int r) const
  {
    int Nquad = ShTri.GetM();
    Vector1 Uquad(Nquad);
    
    Mlt(SeldonTrans, ShLoc(num_loc), U, Uquad);

    switch (num_loc)
      {
      case 0 :
      case 4 :
        if ((r == 0) || (r == order))
	  Mlt(SeldonTrans, ShTri, Uquad, V);
	else
	  Mlt(SeldonTrans, ProjOperatorTriOrder(r), Uquad, V);
        break;
      case 1 :
        if ((r == 0) || (r == order))
	  Copy(Uquad, V);
	else
	  Mlt(SeldonTrans, ProjOperatorQuadOrder(r), Uquad, V);

        break;
      case 2 :
      case 3 :
        if ((r == 0) || (r == order))
	  Mlt(SeldonTrans, ShQuad, Uquad, V);
	else
	  Mlt(SeldonTrans, ProjOperatorQuadOrderY(r), Uquad, V);
	
        break;
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
  void WedgeDgOrtho::ApplyNablaShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh,
					       Vector2& V, int r) const
  {
    int Nquad = ShTri.GetM();
    Vector1 Uquad(2*Nquad), U = Uh;
    bool variable_order = true;
    if ((r == 0) || (r == order))
      variable_order = false;
    
    if (variable_order)
      {
        if (num_loc%4 == 0)
          {
            const VectR2& Dx = ProjOperatorDxtildeTri(r);
            for (int i = 0; i < Dx.GetM(); i++)
              U(3*i) = Dx(i)(0)*U(3*i) + Dx(i)(1)*U(3*i+1);        
          }
        else
          {
            const VectR2& Dx = ProjOperatorDxtildeQuad(r)(num_loc-1);
            for (int i = 0; i < Dx.GetM(); i++)
              U(3*i) = Dx(i)(0)*U(3*i) + Dx(i)(1)*U(3*i+1);
          }
      }
    else
      {
        int offset = nb_points_quadrature_inside + this->offset_faceSh(num_loc);
        for (int i = 0; i < this->GetNbQuadBoundary(num_loc); i++)
          U(3*i) = DerivDxtildeDx(offset+i)*U(3*i) + DerivDxtildeDy(offset+i)*U(3*i+1);        
      }
    
    switch (num_loc)
      {
      case 0 :
      case 4 :
        if (!variable_order)
	  Mlt(DerivShTri, U, Uquad);
	else
	  Mlt(ProjOperatorDerivTriOrder(r), U, Uquad);
	
        break;
      case 1 :
	if (!variable_order)
	  Mlt(DerivShGauss, U, Uquad);
	else
	  Mlt(ProjOperatorDerivQuadOrder(r), U, Uquad);
        
	break;
      case 2 :
      case 3 :
        if (!variable_order)
	  Mlt(DerivShQuad, U, Uquad);
	else
	  Mlt(ProjOperatorDerivQuadOrderY(r), U, Uquad);
        
	break;
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
  void WedgeDgOrtho
  ::ApplyNablaShQuadratureTransposeGen(int num_loc, const Vector1& U, Vector2& V, int r) const
  {
    int Nquad = ShTri.GetM();
    Vector1 Uquad(2*Nquad);
    bool variable_order = true;
    if ((r == 0) || (r == order))
      variable_order = false;
        
    Mlt(SeldonTrans, DerivShLoc(num_loc), U, Uquad);
    
    switch (num_loc)
      {
      case 0 :
      case 4 :
        if (!variable_order)
	  Mlt(SeldonTrans, DerivShTri, Uquad, V);
	else
	  Mlt(SeldonTrans, ProjOperatorDerivTriOrder(r), Uquad, V);
        break;
      case 1 :
        if (!variable_order)
	  Mlt(SeldonTrans, DerivShGauss, Uquad, V);
	else
	  Mlt(SeldonTrans, ProjOperatorDerivQuadOrder(r), Uquad, V);

        break;
      case 2 :
      case 3 :
        if (!variable_order)
	  Mlt(SeldonTrans, DerivShQuad, Uquad, V);
	else
	  Mlt(SeldonTrans, ProjOperatorDerivQuadOrderY(r), Uquad, V);
	
        break;
      }
    
    if (variable_order)
      {
        if (num_loc%4 == 0)
          {
            const VectR2& Dx = ProjOperatorDxtildeTri(r);
            for (int i = 0; i < Dx.GetM(); i++)
              {
                V(3*i+1) += V(3*i)*Dx(i)(1);
                V(3*i) *= Dx(i)(0);
              } 
          }
        else
          {
            const VectR2& Dx = ProjOperatorDxtildeQuad(r)(num_loc-1);
            for (int i = 0; i < Dx.GetM(); i++)
              {
                V(3*i+1) += V(3*i)*Dx(i)(1);
                V(3*i) *= Dx(i)(0);
              } 
          }
      }
    else
      {
        int offset = nb_points_quadrature_inside + this->offset_faceSh(num_loc);
        for (int i = 0; i < this->GetNbQuadBoundary(num_loc); i++)
          {
            V(3*i+1) += V(3*i)*DerivDxtildeDy(offset+i);
            V(3*i) *= DerivDxtildeDx(offset+i);
          }
      }
    
  }

  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  template<class Vector1, class Vector2>
  void WedgeDgOrtho::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    Vector1 feval_weight = feval;
    const VectReal_wp& weights3d = this->WeightsND();
    for (int i = 0; i < feval_weight.GetM(); i++)
      feval_weight(i) *= weights3d(i);
    
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }


  //! computation of projection between finite element of different orders    
  void WedgeDgOrtho
  ::ComputeTriangularInterpolationProjector(const IVect& order_elt,
                                            const ElementReference<Dimension2, 1>& Fb,
                                            const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi;
    ProjOperatorTriOrder.Reallocate(rmax+1);
    ProjOperatorDxtildeTri.Reallocate(rmax+1);
    ProjOperatorDerivTriOrder.Reallocate(rmax+1);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
        int Nquad = (order+1)*(order+1);
        int Ntri = Pts(i).GetM();
	ProjOperatorTriOrder(r).Reallocate(Nquad, Ntri);
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = NumQuad2D(i1, i2);
		  Real_wp yt = Pts(i)(j)(1), xt = Pts(i)(j)(0)/(1.0-yt);
		  ProjOperatorTriOrder(r)(k, j) = gauss_x.EvaluatePhi(i1, xt)*gauss_y.EvaluatePhi(i2, yt);
		}
	  }
        
        ProjOperatorDxtildeTri(r).Reallocate(Ntri);
        Matrix<Real_wp, General, ArrayRowSparse> sh(2*Nquad, 3*Ntri);		
        for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Real_wp yt = Pts(i)(j)(1), xt = Pts(i)(j)(0)/(1.0-yt);
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
                {
                  int k = NumQuad2D(i1, i2);
                  sh.AddInteraction(2*k, 3*j, gauss_x.EvaluatePhiGrad(i1, xt)*gauss_y.EvaluatePhi(i2, yt));
                  sh.AddInteraction(2*k, 3*j+1, gauss_x.EvaluatePhi(i1, xt)*gauss_y.EvaluatePhiGrad(i2, yt));
                  sh.AddInteraction(2*k+1, 3*j+2, gauss_x.EvaluatePhi(i1, xt)*gauss_y.EvaluatePhi(i2, yt));
                }
            
            Real_wp x = Pts(i)(j)(0);
            Real_wp y = Pts(i)(j)(1);
            ProjOperatorDxtildeTri(r)(j)(0) = 1.0/(1.0-y);
            ProjOperatorDxtildeTri(r)(j)(1) = x/square(1.0-y);
          }
        
        Copy(sh, ProjOperatorDerivTriOrder(r));

      }
  }
  

  //! computation of projection between finite element of different orders    
  void WedgeDgOrtho
  ::ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                              const ElementReference<Dimension2, 1>& Fb,
                                              const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi;
    ProjOperatorQuadOrder.Reallocate(rmax+1);
    ProjOperatorQuadOrderY.Reallocate(rmax+1);
    ProjOperatorDxtildeQuad.Reallocate(rmax+1);
    ProjOperatorDerivQuadOrder.Reallocate(rmax+1);
    ProjOperatorDerivQuadOrderY.Reallocate(rmax+1);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
        int Nquad = (order+1)*(order+1);
        int Ntri = Pts(i).GetM();
	ProjOperatorQuadOrder(r).Reallocate(Nquad, Ntri);
	ProjOperatorQuadOrderY(r).Reallocate(Nquad, Ntri);
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = NumQuad2D(i1, i2);
		  Real_wp yt = Pts(i)(j)(1), xt = Pts(i)(j)(0);
		  ProjOperatorQuadOrder(r)(k, j) = gauss_x.EvaluatePhi(i1, xt)*gauss_x.EvaluatePhi(i2, yt);
		  ProjOperatorQuadOrderY(r)(k, j) = gauss_y.EvaluatePhi(i1, xt)*gauss_x.EvaluatePhi(i2, yt);
		}	    
	  }
        
        ProjOperatorDxtildeQuad(r)(0).Reallocate(Ntri);
        ProjOperatorDxtildeQuad(r)(1).Reallocate(Ntri);
        ProjOperatorDxtildeQuad(r)(2).Reallocate(Ntri);
        Matrix<Real_wp, General, ArrayRowSparse>
          sh_quad(2*Nquad, 3*Ntri), sh_gauss(2*Nquad, 3*Ntri);
        for (int j = 0; j < Pts(i).GetM(); j++)
	  {
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = NumQuad2D(i1, i2);
		  Real_wp zt = Pts(i)(j)(1), xt = Pts(i)(j)(0);
                  sh_gauss.AddInteraction(2*k, 3*j, gauss_x.EvaluatePhiGrad(i1, xt)*gauss_x.EvaluatePhi(i2, zt));
                  sh_gauss.AddInteraction(2*k+1, 3*j+1,
                                          gauss_x.EvaluatePhi(i1, xt)*gauss_x.EvaluatePhi(i2, zt));
                  sh_gauss.AddInteraction(2*k, 3*j+2,
                                          gauss_x.EvaluatePhi(i1, xt)*gauss_x.EvaluatePhiGrad(i2, zt));

                  sh_quad.AddInteraction(2*k+1, 3*j, gauss_y.EvaluatePhi(i1, xt)*gauss_x.EvaluatePhi(i2, zt));
                  sh_quad.AddInteraction(2*k, 3*j+1, gauss_y.EvaluatePhiGrad(i1, xt)*gauss_x.EvaluatePhi(i2, zt));
                  sh_quad.AddInteraction(2*k, 3*j+2, gauss_y.EvaluatePhi(i1, xt)*gauss_x.EvaluatePhiGrad(i2, zt));
                }
            
            Real_wp x = Pts(i)(j)(0), y = 0;
            ProjOperatorDxtildeQuad(r)(0)(j)(0) = 1.0/(1.0-y);
            ProjOperatorDxtildeQuad(r)(0)(j)(1) = x/square(1.0-y);

            x = 1.0-Pts(i)(j)(0); y = Pts(i)(j)(0);
            ProjOperatorDxtildeQuad(r)(1)(j)(0) = 1.0/(1.0-y);
            ProjOperatorDxtildeQuad(r)(1)(j)(1) = x/square(1.0-y);

            x = 0; y = Pts(i)(j)(0);
            ProjOperatorDxtildeQuad(r)(2)(j)(0) = 1.0/(1.0-y);
            ProjOperatorDxtildeQuad(r)(2)(j)(1) = x/square(1.0-y);
          }
        
        Copy(sh_gauss, ProjOperatorDerivQuadOrder(r));
        Copy(sh_quad, ProjOperatorDerivQuadOrderY(r));
      }
  }


  //! computation of mass matrix knowing decomposition of Jacobian for non-affine elements
  /*!
    \param[out] A computed mass matrix
    \param[in] CoefJacobian decomposition of jacobian |DF|
   */
  void WedgeDgOrtho::ComputeMassMatrix(Matrix<Real_wp, Symmetric, SymColSkyLine> & A,
                                       const VectReal_wp & CoefJacobian) const
  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> Mm;
    ComputeMassMatrix(Mm, CoefJacobian);    
    
    Copy(Mm, A);
  }
  
  
  //! computation of mass matrix knowing decomposition of Jacobian for non-affine elements
  /*!
    \param[out] Mm computed mass matrix
    \param[in] CoefJacobian decomposition of jacobian |DF|
   */
  void WedgeDgOrtho::ComputeMassMatrix(Matrix<Real_wp, Symmetric, ArrayRowSymSparse> & Mm,
                                       const VectReal_wp & CoefJacobian) const
  {
    PerformComputationMassMatrix(Mm, CoefJacobian);
    
    for (int i = 0; i < Mm.GetM(); i++)
      for (int j = 0; j < Mm.GetRowSize(i); j++)
        if (Mm.Index(i, j) > i)
          Mm.Value(i, j) *= 0.5;
    
    Real_wp M_max = Norm1(Mm);
    Mm.RemoveSmallEntry(10*epsilon_machine*M_max);
  }
  
  
  //! computation of mass matrix knowing decomposition of Jacobian for non-affine elements
  /*!
    \param[out] Mm computed mass matrix
    \param[in] CoefJacobian decomposition of jacobian |DF|
   */
  void WedgeDgOrtho::ComputeMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & Mm,
                                       const VectReal_wp & CoefJacobian) const
  {
    PerformComputationMassMatrix(Mm, CoefJacobian);
    for (int i = 0; i < Mm.GetM(); i++)
      for (int j = i+1; j < Mm.GetM(); j++)
        Mm(i, j) *= 0.5;        
  }


  //! computation of mass matrix knowing weighted jacobian on each quadrature point
  /*!
    \param[out] A computed mass matrix
    \param[in] coef values of jacobian |DF| on quadrature points
   */
  void WedgeDgOrtho::IntegrateMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
                                         const VectReal_wp & coef) const
  {
    A.Reallocate(nb_dof_loc, nb_dof_loc);
    A.Fill(0);
    int Ntri = InvWeight2D.GetM();
    Real_wp vloc;
    for (int i = 0; i < Ntri; i++)
      for (int j = 0; j < Ntri; j++)
        for (int k = 0; k <= order; k++)
          {
            int irow = NumDofTri(i, k);
            int icol = NumDofTri(j, k);
            vloc = 0;
            for (int m = 0; m < NumQuadTri.GetM(); m++)
              vloc += ChTri(i, m)*ChTri(j, m)*coef(NumQuadTri(m, k));
            
            A(irow, icol) = vloc;
          }
  }


  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of basis functions on pointloc
  */  
  void WedgeDgOrtho::ComputeValuesPhiRef(const R3& pointloc, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    R2 pt2D(pointloc(0), pointloc(1));
    VectReal_wp phi_tri;
    function_tri.ComputeValuesPhiOrthoRef(order, NumOrtho2D, InvWeight2D, pt2D, phi_tri);
    for (int k = 0; k <= order; k++)
      {
        Real_wp val = gauss.EvaluatePhi(k, pointloc(2));
        for (int i = 0; i < NumDofTri.GetM(); i++)
          phi(NumDofTri(i, k)) = val*phi_tri(i);
      }        
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on pointloc
  */
  void WedgeDgOrtho::ComputeGradientPhiRef(const R3& pointloc, VectR3& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc);
    R2 pt2D(pointloc(0), pointloc(1));
    VectReal_wp phi_tri; VectR2 grad_phi_tri;
    function_tri.ComputeValuesPhiOrthoRef(order, NumOrtho2D, InvWeight2D, pt2D, phi_tri);
    function_tri.ComputeGradientPhiOrthoRef(order, NumOrtho2D, InvWeight2D, pt2D, grad_phi_tri);
    for (int k = 0; k <= order; k++)
      {
        Real_wp val = gauss.EvaluatePhi(k, pointloc(2));
        Real_wp dval = gauss.EvaluatePhiGrad(k, pointloc(2));
        for (int i = 0; i < NumDofTri.GetM(); i++)
          grad_phi(NumDofTri(i, k))
            .Init(val*grad_phi_tri(i)(0), val*grad_phi_tri(i)(1), dval*phi_tri(i));
      }
  }
          
}

#define MONTJOIE_FILE_WEDGE_DG_ORTHO_CXX
#endif
