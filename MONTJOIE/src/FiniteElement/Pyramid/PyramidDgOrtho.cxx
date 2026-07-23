#ifndef MONTJOIE_FILE_PYRAMID_DG_ORTHO_CXX

// #include "OptPyramidDgOrtho.cxx"

namespace Montjoie
{
  
  //! default constructor
  PyramidDgOrtho::PyramidDgOrtho() : PyramidReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;

    this->discontinuous_element = true;
    
    use_quadrature_for_sh = true;
    use_quadrature_for_rh = true;
    linear_sparse_mass_matrix = true;
  }
  

  //! how to number mesh
  void PyramidDgOrtho::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg != ElementReference_Base::HDG)
      return PyramidReference<1>::ConstructNumberMap(nmap, dg);
    
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

  
  size_t PyramidDgOrtho::GetMemorySize() const
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
    taille += RhLoc.GetMemorySize() + InvWeightFct.GetMemorySize(); 
    taille += DerivDxtildeDx.GetMemorySize() + DerivDxtildeDz.GetMemorySize()
      + DerivDytildeDy.GetMemorySize() + DerivDytildeDz.GetMemorySize();
    
    taille += beta.GetMemorySize() + gamma.GetMemorySize();
    taille += Seldon::GetMemorySize(CmO1) +  Seldon::GetMemorySize(CpO1)
      + Seldon::GetMemorySize(CmO2) + Seldon::GetMemorySize(CpO2) + Seldon::GetMemorySize(CmO1mO2)
      + Seldon::GetMemorySize(CpO1mO2) + Seldon::GetMemorySize(CmO1pO2)
      + Seldon::GetMemorySize(CpO1pO2);
    
    taille += NumFct3D.GetMemorySize() + InvWeightFct.GetMemorySize();    
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
  void PyramidDgOrtho::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					      int rsurf_tri, int rsurf_quad,
					      int type_surf_tri, int type_surf_quad, int gauss_z)
  {    
    PyramidReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    ConstructFunctions();
    ConstructMassMatrix();
    ConstructStiffnessMatrix();
    
    this->Fb_geom.ComputeCoefficientTransformation();
    
    this->ConstructElementaryMatrix(*this);      
  }
  

  //! construction of mass matrix
  void PyramidDgOrtho::ConstructMassMatrix()
  {
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
    int rp1 = order+1;
    int nb_points_z = (order+1)*(order+1)*(order+1);
    Matrix<Real_wp, General, ArrayRowSparse> ch(nb_dof_loc, nb_points_z);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = rp1*(i*rp1 + j) + k;
	    for (int m = 0; m <= order-max(i,j); m++)
              ch.AddInteraction(NumFct3D(i, j, m), node, Pv(max(i,j), m, k));

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
    
    // CopyPvX(Pv);
    ch.Clear();
    ch.Reallocate(nb_points_z, nb_points_z);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = rp1*(i*rp1 + j) + k;
	    for (int m = 0; m <= order; m++)
              ch.AddInteraction(rp1*(i*rp1 + m) + k, node, Pv(0, m, j));
          }
    
    // conversion to csr
    Copy(ch, ChY);

    // third interpolation along x :
    // v^y_i,j,k = \sum_m L_m(xi_i) v^y_{m,j,k}    
    ch.Clear();
    ch.Reallocate(nb_points_z, nb_points_z);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    for (int m = 0; m <= order; m++)
              ch.AddInteraction(rp1*(m*rp1 + j) + k, rp1*(i*rp1 + j) + k, Pv(0, m, i));
          }
    
    // conversion to csr
    Copy(ch, ChX);
    
    // useful coefficients
    beta.Reallocate(order+1); beta.Fill(0);
    gamma.Reallocate(order+1); gamma.Fill(0);
    
    for (int i = 0; i <= order; i++)
      for (int n = 0; n <= order; n++)
	{     
	  if (n > 0) 
	    beta(n) = LegendrePolynom(n,1)*CoefLegendre(n)/CoefLegendre(n-1);
	  
          if (n < order)
	    gamma(n) = CoefLegendre(n)/CoefLegendre(n+1);
	}
    
    VectReal_wp Pn(order+2), PnmO1(order+2), PnpO1(order+2), PnmO2(order+2), PnpO2(order+2);
    VectReal_wp PnmO1mO2(order+2), PnmO1pO2(order+2), PnpO1mO2(order+2), PnpO1pO2(order+2);
    
    CpO1.Reallocate(order+1, order+1);
    CpO2.Reallocate(order+1, order+1);
    CmO1.Reallocate(order+1, order+1);
    CmO2.Reallocate(order+1, order+1);
    
    CpO1pO2.Reallocate(order+1, order+1);
    CpO1mO2.Reallocate(order+1, order+1);
    CmO1pO2.Reallocate(order+1, order+1);
    CmO1mO2.Reallocate(order+1, order+1);

    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  CpO1(i,j).Reallocate(order+1, order+1); CpO1(i,j).Fill(0);
	  CpO2(i,j).Reallocate(order+1, order+1); CpO2(i,j).Fill(0);
	  CmO1(i,j).Reallocate(order+1, order+1); CmO1(i,j).Fill(0);
	  CmO2(i,j).Reallocate(order+1, order+1); CmO2(i,j).Fill(0);
	  
	  CpO1pO2(i,j).Reallocate(order+1, order+1); CpO1pO2(i,j).Fill(0);
	  CpO1mO2(i,j).Reallocate(order+1, order+1); CpO1mO2(i,j).Fill(0);
	  CmO1pO2(i,j).Reallocate(order+1, order+1); CmO1pO2(i,j).Fill(0);
	  CmO1mO2(i,j).Reallocate(order+1, order+1); CmO1mO2(i,j).Fill(0);
	}
    
    VectReal_wp xi, omega;
    ComputeGaussLegendre(xi, omega, order+1);
    for (int k = 0; k <= order+1; k++)
      for (int i = 0; i <= order; i++)
        for (int j = 0; j <= order; j++)
	  {
	    EvaluateJacobiPolynomial(EvenJacobiPolynom(max(i,j)), order, 2*xi(k)-1, Pn);
	    
	    if (i<order) 
	      EvaluateJacobiPolynomial(EvenJacobiPolynom(max(i+1,j)), order, 2*xi(k)-1, PnpO1);
	    if (j<order) 
	      EvaluateJacobiPolynomial(EvenJacobiPolynom(max(i,j+1)), order, 2*xi(k)-1, PnpO2);
	    if (i>0) 
	      EvaluateJacobiPolynomial(EvenJacobiPolynom(max(i-1,j)), order, 2*xi(k)-1, PnmO1);
	    if (j>0) 
	      EvaluateJacobiPolynomial(EvenJacobiPolynom(max(i,j-1)), order, 2*xi(k)-1, PnmO2);

	    if ((i<order)&&(j<order)) 
	      EvaluateJacobiPolynomial(EvenJacobiPolynom(max(i+1,j+1)),
                                       order, 2*xi(k)-1, PnpO1pO2);
	    if ((i<order)&&(j>0)) 
	      EvaluateJacobiPolynomial(EvenJacobiPolynom(max(i+1,j-1)),
                                       order, 2*xi(k)-1, PnpO1mO2);
	    if ((i>0)&&(j<order)) 
	      EvaluateJacobiPolynomial(EvenJacobiPolynom(max(i-1,j+1)),
                                       order, 2*xi(k)-1, PnmO1pO2);
	    if ((i>0)&&(j>0)) 
	      EvaluateJacobiPolynomial(EvenJacobiPolynom(max(i-1,j-1)),
                                       order, 2*xi(k)-1, PnmO1mO2);
	    
	    for (int m = 0; m < Pn.GetM(); m++)
	      {
		Pn(m) *= CoefEvenJacobi(max(i,j),m);
		
		if (max(i+1,j) <= order) 
		  PnpO1(m) *= CoefEvenJacobi(max(i+1,j),m);
		if (max(i,j+1) <= order)
		  PnpO2(m) *= CoefEvenJacobi(max(i,j+1),m);
		if (max(i-1,j) >= 0) 
		  PnmO1(m) *= CoefEvenJacobi(max(i-1,j),m);
		if (max(i,j-1) >= 0)
		  PnmO2(m) *= CoefEvenJacobi(max(i,j-1),m);
		
		if ((i<order)&&(j<order)) 
		  PnpO1pO2(m) *= CoefEvenJacobi(max(i+1,j+1),m);
		if ((i<order)&&(j>0)) 
		  PnpO1mO2(m) *= CoefEvenJacobi(max(i+1,j-1),m);
		if ((i>0)&&(j<order)) 
		  PnmO1pO2(m) *= CoefEvenJacobi(max(i-1,j+1),m);
		if ((i>0)&&(j>0)) 
		  PnmO1mO2(m) *= CoefEvenJacobi(max(i-1,j-1),m);
	      }

	    for (int i3 = 0; i3 <= order; i3++)
	      for (int j3 = 0; j3 <= order; j3++)
		{
		  if (i<order) 
		    CpO1(i,j)(i3,j3) += omega(k)*pow(1-xi(k),max(i,j)
                                                     +max(i+1,j)+2)*Pn(i3)*PnpO1(j3);
		  if (i>0) 
		    CmO1(i,j)(i3,j3) += omega(k)*pow(1-xi(k),max(i,j)
                                                     +max(i-1,j)+2)*Pn(i3)*PnmO1(j3);
		  if (j<order) 
		    CpO2(i,j)(i3,j3) += omega(k)*pow(1-xi(k),max(i,j)
                                                     +max(i,j+1)+2)*Pn(i3)*PnpO2(j3);
		  if (j>0) 
		    CmO2(i,j)(i3,j3) += omega(k)*pow(1-xi(k),max(i,j)
                                                     +max(i,j-1)+2)*Pn(i3)*PnmO2(j3);
		  
		  if ((i<order)&&(j<order)) 
		    CpO1pO2(i,j)(i3,j3) += omega(k)*pow(1-xi(k),max(i,j)
                                                        +max(i+1,j+1)+2)*Pn(i3)*PnpO1pO2(j3);
		  if ((i<order)&&(j>0)) 
		    CpO1mO2(i,j)(i3,j3) += omega(k)*pow(1-xi(k),max(i,j)
                                                        +max(i+1,j-1)+2)*Pn(i3)*PnpO1mO2(j3);
		  if ((i>0)&&(j<order)) 
		    CmO1pO2(i,j)(i3,j3) += omega(k)*pow(1-xi(k),max(i,j)
                                                        +max(i-1,j+1)+2)*Pn(i3)*PnmO1pO2(j3);
		  if ((i>0)&&(j>0)) 
		    CmO1mO2(i,j)(i3,j3) += omega(k)*pow(1-xi(k),max(i,j)
                                                        +max(i-1,j-1)+2)*Pn(i3)*PnmO1mO2(j3);
		}
	  }
  }

  
  //! construction of stiffness matrix
  void PyramidDgOrtho::ConstructStiffnessMatrix()
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
	    int node = rp1*(i*rp1+j)+k;
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
    
    // construction of ShLoc
    int Nquad = (order+1)*(order+1);
    rh.Clear();
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
	// face y = -(1-z)
	int m = nb_points_quadrature_inside + Nquad + i;
	Real_wp x = points3d(m)(0), y = points3d(m)(1), z = points3d(m)(2);
	Real_wp xt = 0.5*x/(1.0-z) + 0.5, yt = 0.5*y/(1.0-z) + 0.5, zt = z;
	
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


  //! constructing basis functions
  void PyramidDgOrtho::ConstructFunctions()
  {
    NumFct3D.Reallocate(order+1, order+1, order+1);
    bool renum_ortho = false;
    if (order <= 10)
      {
        ifstream file_in("src/FiniteElement/Pyramid/ReorderingPyramidOrtho.txt");
        string ligne, line_search;
        line_search = string("Reordering  of  the  dof  ")
          + "for  the  reference  pyramid  of  order  p=" + Seldon::to_str(order);
        
        if (!file_in.is_open())
          {
	    // cout<<" File containing reordering not found"<<endl;
	  }
        else
          {
            renum_ortho = true;
            
            while (!file_in.eof())
              {
                getline(file_in, ligne);
                DeleteSpaceAtExtremityOfString(ligne);
                
                if (!ligne.compare(line_search))
                  {
                    int nb;
                    for (int k = 0; k <= order; k++)
                      for (int i = 0; i <= order-k; i++)
                        for (int j = 0; j <= order-k; j++)
                          {
                            file_in >> nb;
                            NumFct3D(i, j, k) = nb-1;                        
                          }
                  }
              }
          }
      }
    
    if (!renum_ortho)
      {
	int nb = 0;
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    for (int k = 0; k <= order-max(i,j); k++)
	      {
		NumFct3D(i, j, k) = nb;
		nb++;
	      }
      }
    
    InvWeightFct.Reallocate((order+1)*(order+2)*(2*order+3)/6);
    const Matrix<Real_wp>& CoefEvenJacobi = Fb_geom.GetCoefEvenJacobi();
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefLegendre();
    for (int k = 0; k <= order; k++)
      for (int i = 0; i <= order-k; i++)
        for (int j = 0; j <= order-k; j++)
          InvWeightFct(NumFct3D(i, j, k)) = CoefLegendre(i)*CoefLegendre(j)
            *CoefEvenJacobi(max(i,j), k);
        
    nb_dof_loc = InvWeightFct.GetM();

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
  
  
  //! computation of mass matrix knowing decomposition of Jacobian for non-affine elements
  /*!
    \param[out] Mm computed mass matrix
    \param[in] CoefJacobian decomposition of jacobian |DF|
   */
  template<class Matrix1>
  void PyramidDgOrtho
  ::PerformComputationMassMatrix(Matrix1 & Mm, const VectReal_wp & CoefJacobian) const
  {
    Mm.Reallocate(nb_dof_loc,nb_dof_loc); Mm.Fill(0);
    
    Real_wp A1 = CoefJacobian(0);
    Real_wp B1 = CoefJacobian(1), B2 = CoefJacobian(2);
    Real_wp C1 = CoefJacobian(3);
    
    int i;
    
    for (int i3 = 0; i3 <= order; i3++)
      for (int i2 = 0; i2 <= order-i3; i2++)
        for (int i1 = 0; i1 <= order-i3; i1++)
          { 
            i = NumFct3D(i1,i2,i3);
	    
            Mm.Get(i,i) += A1;
	    
            for (int j3 = 0; j3 <= order-max(i1-1,i2); j3++)
              if (i1>0)
                Mm.Get(i,NumFct3D(i1-1,i2,j3)) += B1*beta(i1)*CmO1(i1,i2)(i3,j3);
            
            for (int j3 = 0; j3 <= order-max(i1+1,i2); j3++)
              if (i1<order-j3)
                Mm.Get(i,NumFct3D(i1+1,i2,j3)) += B1*gamma(i1)*CpO1(i1,i2)(i3,j3);
	    
            for (int j3 = 0; j3 <= order-max(i1,i2-1); j3++)
              if (i2>0)
                Mm.Get(i,NumFct3D(i1,i2-1,j3)) += B2*beta(i2)*CmO2(i1,i2)(i3,j3);
            
            for (int j3 = 0; j3 <= order-max(i1,i2+1); j3++)
              if (i2<order-j3)
                Mm.Get(i,NumFct3D(i1,i2+1,j3)) += B2*gamma(i2)*CpO2(i1,i2)(i3,j3);	    
	    
            for (int j3 = 0; j3 <= order-max(i1-1,i2-1); j3++)
              if (i1>0)
                if (i2>0)
                  Mm.Get(i,NumFct3D(i1-1,i2-1,j3)) += C1*beta(i1)*beta(i2)*CmO1mO2(i1,i2)(i3,j3);
	    
            for (int j3 = 0; j3 <= order-max(i1-1,i2+1); j3++)
              if (i1>0)
                if (i2<order-j3)
                  Mm.Get(i,NumFct3D(i1-1,i2+1,j3)) += C1*beta(i1)*gamma(i2)*CmO1pO2(i1,i2)(i3,j3);
	    
            for (int j3 = 0; j3 <= order-max(i1+1,i2-1); j3++)
              if (i1<order-j3)
                if (i2>0)
                  Mm.Get(i,NumFct3D(i1+1,i2-1,j3)) += C1*gamma(i1)*beta(i2)*CpO1mO2(i1,i2)(i3,j3);
	    
            for (int j3 = 0; j3 <= order-max(i1+1,i2+1); j3++)
              if (i1<order-j3)
                if (i2<order-j3)
                  Mm.Get(i,NumFct3D(i1+1,i2+1,j3)) += C1*gamma(i1)*gamma(i2)*CpO1pO2(i1,i2)(i3,j3);
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
  void PyramidDgOrtho::ApplyChGen(const Vector1& U, Vector2& V) const
  {
    //Vector1 Ux(ChX.GetM()), Uy(ChY.GetM());
    typedef typename Vector1::value_type Complexe;
    Vector1& Ux = PyramidDgOrthoVariables<Complexe>::Ux(order);
    Vector1& Uy = PyramidDgOrthoVariables<Complexe>::Uy(order);
    
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
  void PyramidDgOrtho::ApplyChTransposeGen(const Vector1& U, Vector2& V) const
  {
    //Vector1 Uy(ChY.GetN()), Uz(ChZ.GetN());
    typedef typename Vector1::value_type Complexe;
    Vector1& Uy = PyramidDgOrthoVariables<Complexe>::Uy2(order);
    Vector1& Uz = PyramidDgOrthoVariables<Complexe>::Uz(order);
    
    Uy.Reallocate(ChY.GetN());
    Uz.Reallocate(ChZ.GetN());
    
    Mlt(SeldonTrans, ChZ, U, Uz);
    Mlt(SeldonTrans, ChY, Uz, Uy);
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
  void PyramidDgOrtho::ApplyRhQuadratureGen(const Vector1& U, Vector2& V) const
  {
    // applying chaine rule to get derivatives on symmetric pyramid
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
  void PyramidDgOrtho::ApplyRhQuadratureTransposeGen(const Vector1& U, Vector2& V) const
  {    
    Mlt(SeldonTrans, RhLoc, U, V);
    
    // applying chaine rule to get derivatives on symmetric pyramid
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
  void PyramidDgOrtho::ApplyShQuadratureGen(const T0& alpha, int num_loc, const Vector1& U,
                                         Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    //Vector1 Uquad(Nquad);
    typedef typename Vector1::value_type Complexe;
    Vector1& Uquad = PyramidDgOrthoVariables<Complexe>::Uquad(order);
    Uquad.Reallocate(Nquad);

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
  void PyramidDgOrtho::ApplyShQuadratureTransposeGen(int num_loc, const Vector1& U,
						     Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    // Vector1 Uquad(Nquad);
    typedef typename Vector1::value_type Complexe;
    Vector1& Uquad = PyramidDgOrthoVariables<Complexe>::Uquad(order);
    Uquad.Reallocate(Nquad);
    
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
  void PyramidDgOrtho::ApplyNablaShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh,
						 Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    //Vector1 Uquad(Nquad);
    //typedef typename Vector1::value_type Complexe;
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
  void PyramidDgOrtho::ApplyNablaShQuadratureTransposeGen(int num_loc, const Vector1& U,
							  Vector2& V, int r) const
  {
    int Nquad = ShSurf(0).GetM();
    // Vector1 Uquad(Nquad);
    //typedef typename Vector1::value_type Complexe;
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
  void PyramidDgOrtho::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    Vector1 feval_weight = feval;
    const VectReal_wp& weights3d = this->WeightsND();
    for (int i = 0; i < feval_weight.GetM(); i++)
      feval_weight(i) *= weights3d(i);
    
    ApplyCh(feval_weight, contrib);
  }

  
  //! computation of projection between finite element of different orders    
  void PyramidDgOrtho::
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
    ProjOperatorTriOrder.Reallocate(rmax+1);
    ProjOperatorDxtildeTri.Reallocate(rmax+1);
    ProjOperatorDerivShTri.Reallocate(rmax+1);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();    
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
  void PyramidDgOrtho::
  ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                            const ElementReference<Dimension2, 1>& Fb,
                                            const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    ProjOperatorQuadOrder.Reallocate(rmax+1);
    ProjOperatorDxtildeQuad.Reallocate(rmax+1);
    ProjOperatorDerivShQuad.Reallocate(rmax+1);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();    
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


  //! computation of mass matrix knowing decomposition of Jacobian for non-affine elements
  /*!
    \param[out] A computed mass matrix
    \param[in] CoefJacobian decomposition of jacobian |DF|
   */
  void PyramidDgOrtho::ComputeMassMatrix(Matrix<Real_wp, Symmetric, SymColSkyLine> & A,
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
  void PyramidDgOrtho::ComputeMassMatrix(Matrix<Real_wp, Symmetric, ArrayRowSymSparse> & Mm,
                                         const VectReal_wp & CoefJacobian) const
  {
    PerformComputationMassMatrix(Mm, CoefJacobian);
    
    for (int i = 0; i < Mm.GetM(); i++)
      for (int j = 0; j < Mm.GetRowSize(i); j++)
        if (Mm.Index(i, j) > i)
          Mm.Value(i, j) *= 0.5;
    
    Real_wp max_M = Norm1(Mm);
    Mm.RemoveSmallEntry(10.0*epsilon_machine*max_M);
  }
  
  
  //! computation of mass matrix knowing decomposition of Jacobian for non-affine elements
  /*!
    \param[out] Mm computed mass matrix
    \param[in] CoefJacobian decomposition of jacobian |DF|
  */
  void PyramidDgOrtho::ComputeMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & Mm,
                                         const VectReal_wp & CoefJacobian) const
  {
    PerformComputationMassMatrix(Mm, CoefJacobian);
    for (int i = 0; i < Mm.GetM(); i++)
      for (int j = i+1; j < Mm.GetM(); j++)
        Mm(i, j) *= 0.5;        
  }
  
}

#define MONTJOIE_FILE_PYRAMID_DG_ORTHO_CXX
#endif
