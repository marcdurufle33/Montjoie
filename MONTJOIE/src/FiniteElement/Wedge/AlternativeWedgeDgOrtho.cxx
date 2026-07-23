#ifndef MONTJOIE_FILE_WEDGE_DG_ORTHO_CXX

namespace Montjoie
{
  //! default constructor
  WedgeDgOrtho::WedgeDgOrtho() : WedgeReference()
  {
    quadrature_equal_nodal = false;
    dof_equal_nodal = true;
    dof_equal_quadrature = false;

    this->linear_sparse_mass_matrix = true;
    this->use_quadrature_for_sh = true;
  }

  
  //! constructions of basis functions
  void WedgeDgOrtho::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					    int rsurf_tri, int rsurf_quad,
					    int type_surf_tri, int type_surf_quad, int gauss_z);
  {
    order = r;    
        
    function_tri.type_interpolation = function_tri.LOBATTO_ELEMENT;
    function_tri.ConstructFiniteElement(r);
    ConstructNodalPoints();
    ConstructQuadrature(r);
    ComputeLagrangianFunctions();

    points_dof3d = points_nodal3d;
    ConstructFunctions();    
    
    VectR3 grad_phi;
    GradientPhi_Nodal.Reallocate(points_nodal3d.GetM(), points_nodal3d.GetM());
    for (int j = 0; j < points_nodal3d.GetM(); j++)
      {
        ComputeGradientPhiNodalRef(points_nodal3d(j), grad_phi);
        
        for (int k = 0; k < points_nodal3d.GetM(); k++)
          GradientPhi_Nodal(k, j) = grad_phi(k);
      }
    
    ConstructMassMatrix();
    ConstructStiffnessMatrix();
    ComputeCoefficientTransformation();
    
    this->ConstructStiffnessMatrixDG(*this);
    
  }
  
  
  //! construction of quadrature points and nodal points
  void WedgeDgOrtho::ConstructQuadrature(int order_)
  {
    // quadrature points on the faces
    Globatto<Dimension2> face_gauss_tri;
    face_gauss_tri.ConstructQuadrature(order, face_gauss_tri.QUADRATURE_TRIANGLE_GAUSS);
    points2d_tri = face_gauss_tri.Points2D();
    weights2d_tri = face_gauss_tri.Weights2D();

    ComputeGaussLegendre(points1d, weights1d, order);
    ComputeGaussJacobi(points1d_y, weights1d_y, order, 1.0, 0.0);
    points2d_quad.Reallocate((order+1)*(order+1));
    weights2d_quad.Reallocate((order+1)*(order+1));
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  points2d_quad(NumNodes2D_quad(i, j)).Init(points1d(i), points1d(j));
	  weights2d_quad(NumNodes2D_quad(i,j)) = weights1d(i)*weights1d(j);
	}
    
    // tensorized quadrature points on the volume
    // quadrature points based on Gauss-Jacobi rules (with tensorization)
    Matrix<int> coor;
    Mesh<Dimension3>::ConstructHexahedralNumbering(order, NumQuad3D, coor);
    points3d.Reallocate((order+1)*(order+1)*(order+1));
    weights3d.Reallocate(points3d.GetM());
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)      
	  {
	    int node = NumQuad3D(i, j, k);
	    points3d(node)(0) = points1d(i)*(1.0-points1d_y(j));
	    points3d(node)(1) = points1d_y(j);
	    points3d(node)(2) = points1d(k);
	    weights3d(node) = weights1d(i)*weights1d_y(j)*weights1d(k);
	  }
    
    nb_points_quadrature_inside = points3d.GetM();    
    // quadrature points on all the faces
    nb_points_quadrature_boundaries = 3*points2d_quad.GetM() + 2*points2d_tri.GetM();
   
    // Points of integration for the first face z = 0 (triangle)
    for (int i = 0; i < points2d_tri.GetM(); i++)
      points3d.PushBack(R3(points2d_tri(i)(0), points2d_tri(i)(1), Real_wp(0)));
    
    // Points of integration for the second face y = 0 (quadrangle)
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d.PushBack(R3(points2d_quad(i)(0), 0, points2d_quad(i)(1)));

    // Points of integration for the second face x+y = 1 (quadrangle)
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d.PushBack(R3(Real_wp(1)-points2d_quad(i)(0),
                           points2d_quad(i)(0), points2d_quad(i)(1)));
    
    // Points of integration for the fourth face x = 0 (quadrangle)
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d.PushBack(R3(0, points2d_quad(i)(0), points2d_quad(i)(1)));
    
    // Points of integration for the fifth face z = 1 (triangle)
    for (int i = 0; i < points2d_tri.GetM(); i++)
      points3d.PushBack(R3(points2d_tri(i)(0), points2d_tri(i)(1), Real_wp(1)));
    
    
    int ind = nb_points_quadrature_inside;
    this->num_quad_points_surf.Reallocate(5);
    this->num_quad_points_surf(0).Reallocate(points2d_tri.GetM());
    this->num_quad_points_surf(4).Reallocate(points2d_tri.GetM());
    for (int n = 1; n < 4; n++)
      this->num_quad_points_surf(n).Reallocate(points2d_quad.GetM());

    // first face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      this->num_quad_points_surf(0)(i) = ind++;
    
    // second face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      this->num_quad_points_surf(1)(i) = ind++;
    
    // third face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      this->num_quad_points_surf(2)(i) = ind++;
    
    // fourth face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      this->num_quad_points_surf(3)(i) = ind++;
      
    // fifth face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      this->num_quad_points_surf(4)(i) = ind++;
    
    offset_faceSh.Reallocate(6);
    offset_faceSh(0) = 0;
    for (int i = 1; i < 4; i++)
      offset_faceSh(i) = points2d_tri.GetM() + (i-1)*points2d_quad.GetM();
    
    offset_faceSh(4) = points2d_tri.GetM() + 3*points2d_quad.GetM();
    
    nb_points_quadrature_boundaries = 3*points2d_quad.GetM() + 2*points2d_tri.GetM();
    offset_faceSh(5) = nb_points_quadrature_boundaries;

  }
  
      
  //! construction of mass matrix
  void WedgeDgOrtho::ConstructMassMatrix()
  {
    Value_Phi.Reallocate(nb_dof_loc, points3d.GetM());
    VectReal_wp phi;
    for (int i = 0; i < points3d.GetM(); i++)
      {
	ComputeValuesPhiRef(points3d(i), phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  Value_Phi(j,i) = phi(j);
      }
    
    // computation of ChX, ChY, ChZ
    
    // first interpolation along z :
    // v^z_i,j,k = \sum_m L_m(xi_k) u_{i,j,m}
    Array3D<Real_wp> Pv(order+1, order+1, order+1);
    Pv.Fill(0); VectReal_wp Pm;
    for (int k = 0; k <= order; k++)
      {
	Real_wp xi_k = 2.0*points1d(k) - 1.0;
	EvaluateJacobiPolynomial(LegendrePolynom, order, xi_k, Pm);
	for (int j = 0; j < Pm.GetM(); j++)
	  Pv(0, j, k) = Pm(j)*CoefLegendre(j);
      }
    
    Matrix<Real_wp, General, ArrayRowSparse>
      ch(nb_dof_loc, nb_dof_loc), ch_trans(nb_dof_loc, nb_dof_loc);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = NumOrtho3D(i, j, k);
	    for (int m = 0; m <= order; m++)
	      {
		ch.AddInteraction(NumOrtho3D(i, j, m), node, Pv(0, m, k));
		ch_trans.AddInteraction(node, NumOrtho3D(i, j, m), Pv(0, m, k));
	      }
	  }
    
    // conversion to csr
    Copy(ch, ChZ); Copy(ch_trans, ChTransposeZ);
    
    // second interpolation along y :
    // v^y_i,j,k = \sum_m P_m^(2i+1)(xi_j) ((1-xi_j)/2)^i v^z_{i,m,k}
    Pv.Fill(0);
    for (int k = 0; k <= order; k++)
      {
	Real_wp xi_k = 2.0*points1d_y(k) - 1.0;
	for (int i = 0; i <= order; i++)
	  {
	    EvaluateJacobiPolynomial(function_tri.OddJacobiPolynom(i), order-i, xi_k, Pm);
	    for (int j = 0; j < Pm.GetM(); j++)
	      Pv(i, j, k) = Pm(j)*function_tri.CoefOddJacobi(i, j)*pow(0.5*(1.0-xi_k), i);
	  }
      }
    
    int nb_points_y = (order+1)*(order+1)*(order+1);
    ch.Clear(); ch_trans.Clear();
    ch.Reallocate(nb_dof_loc, nb_points_y);
    ch_trans.Reallocate(nb_points_y, nb_dof_loc);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = NumQuad3D(i, j, k);
	    for (int m = 0; m <= order-i; m++)
	      {
		ch.AddInteraction(NumOrtho3D(i, m, k), node, Pv(i, m, j));
		ch_trans.AddInteraction(node, NumOrtho3D(i, m, k), Pv(i, m, j));
	      }
	  }
    
    // conversion to csr
    Copy(ch, ChY); Copy(ch_trans, ChTransposeY);

    // third interpolation along x :
    // v_i,j,k = \sum_m L_m(xi_i) v^y_{m,j,k}
    Pv.Fill(0);
    for (int k = 0; k <= order; k++)
      {
	Real_wp xi_k = 2.0*points1d(k) - 1.0;
	EvaluateJacobiPolynomial(LegendrePolynom, order, xi_k, Pm);
	for (int j = 0; j < Pm.GetM(); j++)
	  Pv(0, j, k) = Pm(j)*CoefLegendre(j);
      }
    
    ch.Clear(); ch_trans.Clear();
    ch.Reallocate(nb_points_y, nb_points_y);
    ch_trans.Reallocate(nb_points_y, nb_points_y);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    for (int m = 0; m <= order; m++)
	      {
		ch.AddInteraction(NumQuad3D(m, j, k), NumQuad3D(i,j,k), Pv(0, m, i));
		ch_trans.AddInteraction(NumQuad3D(i,j,k), NumQuad3D(m, j, k), Pv(0, m, i));
	      }
	  }
    
    // conversion to csr
    Copy(ch, ChX); Copy(ch_trans, ChTransposeX);

    // usefull coefficients
    alpha.Reallocate(order+1,2*order+2); alpha.Fill(0);
    beta.Reallocate(order+1,2*order+2); beta.Fill(0);
    gamma.Reallocate(order+1,2*order+2); gamma.Fill(0);
    
    for (int i = 0; i <= order; i++)
      {
	alpha(0, 0) = function_tri.LegendrePolynom(0,0);
	alpha(0,2*i+1) = function_tri.OddJacobiPolynom(i)(0,0);
	
	gamma(0,0) = function_tri.CoefLegendre(0)/function_tri.CoefLegendre(1);
	gamma(0,2*i+1) = function_tri.CoefOddJacobi(i,0)/function_tri.CoefOddJacobi(i,1);
	
	for (int n = 1; n <= order-i; n++)
	  {     
	    alpha(n,0) = function_tri.LegendrePolynom(n,0);
	    alpha(n,2*i+1) = function_tri.OddJacobiPolynom(i)(n,0);
	    
	    beta(n,0) = function_tri.LegendrePolynom(n,1)
              *function_tri.CoefLegendre(n)/function_tri.CoefLegendre(n-1);
	    beta(n,2*i+1) = function_tri.OddJacobiPolynom(i)(n,1)
              *function_tri.CoefOddJacobi(i,n)/function_tri.CoefOddJacobi(i,n-1);
	    
	    gamma(n,0) = function_tri.CoefLegendre(n)/function_tri.CoefLegendre(n+1);
	    gamma(n,2*i+1) = function_tri.CoefOddJacobi(i,n)/function_tri.CoefOddJacobi(i,n+1);
	  }
      }

    VectReal_wp Pn(order+2),PnmO(order+2),PnpO(order+2);
    C.Reallocate(order+1,order+1,order+1); C.Fill(0);
    CpO.Reallocate(order+1,order+1,order+1); CpO.Fill(0);
    CmO.Reallocate(order+1,order+1,order+1); CmO.Fill(0);

    VectReal_wp xi, omega;
    ComputeGaussLegendre(xi, omega, order+1);
    for (int k = 0; k <= order+1; k++)
      for (int i = 0; i <= order; i++)
	{
	  EvaluateJacobiPolynomial(function_tri.OddJacobiPolynom(i), order-i, 2*xi(k)-1, Pn);
	  EvaluateJacobiPolynomial(function_tri.OddJacobiPolynom(i+1), order-1-i, 2*xi(k)-1, PnpO);
	  if (i==0)
	    PnmO = 0;
	  else 
	    {
	      EvaluateJacobiPolynomial(function_tri.OddJacobiPolynom(i-1),
                                       order+1-i, 2*xi(k)-1, PnmO);
	      for (int j = 0; j < PnmO.GetM(); j++)
		PnmO(j) *= function_tri.CoefOddJacobi(i-1,j);
	    }
	  
	  for (int j = 0; j < Pn.GetM(); j++)
	    Pn(j) *= function_tri.CoefOddJacobi(i,j);
	  
	  for (int j = 0; j < PnpO.GetM(); j++)
	    PnpO(j) *= function_tri.CoefOddJacobi(i+1,j);
	  
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
    Gradient_Phi.Reallocate(nb_dof_loc, points3d.GetM());
    VectR3 grad_phi;
    for (int i = 0; i < points3d.GetM(); i++)
      {
	ComputeGradientPhiRef(points3d(i), grad_phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  Gradient_Phi(j, i) = grad_phi(j);
      }        
    
    /*
    VectReal_wp phi_xm, phi_xp, phi_ym, phi_yp, phi_zm, phi_zp; R3 pt;
    for (int j = 0; j < nb_points_quadrature_inside; j++)
      {
	Real_wp h = pow(epsilon_machine, 1.0/3.0);
	ComputeGradientPhiRef(points3d(j), grad_phi);
	
	pt = points3d(j); pt(0) -= h;
	ComputeValuesPhiRef(pt, phi_xm);

	pt = points3d(j); pt(0) += h;
	ComputeValuesPhiRef(pt, phi_xp);

	pt = points3d(j); pt(1) -= h;
	ComputeValuesPhiRef(pt, phi_ym);

	pt = points3d(j); pt(1) += h;
	ComputeValuesPhiRef(pt, phi_yp);

	pt = points3d(j); pt(2) -= h;
	ComputeValuesPhiRef(pt, phi_zm);

	pt = points3d(j); pt(2) += h;
	ComputeValuesPhiRef(pt, phi_zp);
	
	for (int i = 0; i < nb_dof_loc; i++)
	  {
	    R3 grad_num;
	    grad_num(0) = (phi_xp(i) - phi_xm(i))/(2.0*h);
	    grad_num(1) = (phi_yp(i) - phi_ym(i))/(2.0*h);
	    grad_num(2) = (phi_zp(i) - phi_zm(i))/(2.0*h);
	    
	    if (grad_num.Distance(grad_phi(i)) > sqrt(epsilon_machine)*10.0)
	      abort();
	  }	
      }
    */
    
    // Computation of gradient from values on quadrature points
    Globatto<Dimension1> gauss_x, gauss_y;
    gauss_x.AffectPoints(points1d);
    gauss_y.AffectPoints(points1d_y);
    gauss_x.ComputeGradPhi(1e3*epsilon_machine);
    gauss_y.ComputeGradPhi(1e3*epsilon_machine);
    
    int N = (order+1)*(order+1)*(order+1);
    Matrix<Real_wp, General, ArrayRowSparse> rh, rh_trans;
    rh.Reallocate(N, 3*N);
    rh_trans.Reallocate(3*N, N);
    DerivDxtildeDx.Reallocate(N); DerivDxtildeDy.Reallocate(N);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = NumQuad3D(i, j, k);
	    for (int m = 0; m <= order; m++)
	      {
		int node2 = 3*NumQuad3D(m, j, k);
		rh_trans.AddInteraction(node2, node, gauss_x.GradPhi(i, m));
		rh.AddInteraction(node, node2, gauss_x.GradPhi(i, m));

		node2 = 3*NumQuad3D(i, m, k)+1;
		rh_trans.AddInteraction(node2, node, gauss_y.GradPhi(j, m));
		rh.AddInteraction(node, node2, gauss_y.GradPhi(j, m));

		node2 = 3*NumQuad3D(i, j, m)+2;
		rh_trans.AddInteraction(node2, node, gauss_x.GradPhi(k, m));
		rh.AddInteraction(node, node2, gauss_x.GradPhi(k, m));
	      }
	    
	    Real_wp x = points3d(node)(0), y = points3d(node)(1);
	    DerivDxtildeDx(node) = 1.0/(1.0-y);
	    DerivDxtildeDy(node) = x/square(1.0-y);
	  }
    
    Copy(rh, RhLoc); Copy(rh_trans, RhLocTranspose);
    
    // construction of ShLoc and ShLocTranspose
    int Nquad = (order+1)*(order+1);
    rh.Clear(); rh_trans.Clear();
    rh.Reallocate(N, 5*Nquad);
    rh_trans.Reallocate(5*Nquad, N);
    VectReal_wp ValPhiX0(order+1), ValPhiX1(order+1), ValPhiY0(order+1);
    
    for (int i = 0; i <= order; i++)
      {
	ValPhiX0(i) = gauss_x.EvaluatePhi(i, 0.0);
	ValPhiX1(i) = gauss_x.EvaluatePhi(i, 1.0);
	ValPhiY0(i) = gauss_y.EvaluatePhi(i, 0.0);
      }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int nv = NumQuad3D(i, j, k);
	    int ns = NumNodes2D_quad(i, j);
	    rh_trans.AddInteraction(ns, nv, ValPhiX0(k));
	    rh.AddInteraction(nv, ns, ValPhiX0(k));
	    
	    ns = NumNodes2D_quad(i, k) + Nquad;
	    rh_trans.AddInteraction(ns, nv, ValPhiY0(j));
	    rh.AddInteraction(nv, ns, ValPhiY0(j));
	    
	    ns = NumNodes2D_quad(j, k) + 2*Nquad;
	    rh_trans.AddInteraction(ns, nv, ValPhiX1(i));
	    rh.AddInteraction(nv, ns, ValPhiX1(i));
	    
	    ns = NumNodes2D_quad(j, k) + 3*Nquad;
	    rh_trans.AddInteraction(ns, nv, ValPhiX0(i));
	    rh.AddInteraction(nv, ns, ValPhiX0(i));

	    ns = NumNodes2D_quad(i, j) + 4*Nquad;
	    rh_trans.AddInteraction(ns, nv, ValPhiX1(k));
	    rh.AddInteraction(nv, ns, ValPhiX1(k));
	  }
    
    Copy(rh, ShLoc); Copy(rh_trans, ShLocTranspose);
    rh.Clear(); rh_trans.Clear();
    
    // computing ShQuad, ShQuadTranspose
    rh.Reallocate(Nquad, Nquad); rh_trans.Reallocate(Nquad, Nquad);
    Matrix<Real_wp> ValPhiY(order+1, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	ValPhiY(i, j) = gauss_y.EvaluatePhi(i, gauss_x.Points(j));
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  int ns = NumNodes2D_quad(i, j);
	  for (int m = 0; m <= order; m++)
	    {
	      rh_trans.AddInteraction(ns, NumNodes2D_quad(m, j), ValPhiY(m, i));
	      rh.AddInteraction(NumNodes2D_quad(m, j), ns, ValPhiY(m, i));
	    }
	}
    
    Copy(rh, ShQuad); Copy(rh_trans, ShQuadTranspose);
    rh.Clear(); rh_trans.Clear();
    
    // now computing ShSurf and ShSurfTranspose
    int Ntri = points2d_tri.GetM();
    for (int k = 0; k < 2; k++)
      {
	ShSurfTranspose(k).Reallocate(Ntri, Nquad); 
	ShSurf(k).Reallocate(Nquad, Ntri); 
      }
    
    for (int i = 0; i < points2d_tri.GetM(); i++)
      {
	int m = nb_points_quadrature_inside + i;
	Real_wp x = points3d(m)(0), y = points3d(m)(1), z = points3d(m)(2);
	Real_wp xt = x/(1.0-y), yt = y, zt = z;
	
	// face z = 0
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_x.EvaluatePhi(j, xt)*gauss_y.EvaluatePhi(k, yt);
	      int col = NumNodes2D_quad(j, k);
	      ShSurfTranspose(0)(i, col) = vloc;
	      ShSurf(0)(col, i) = vloc;
	    }
	
	m = nb_points_quadrature_inside + Ntri + 3*Nquad + i;
	x = points3d(m)(0); y = points3d(m)(1); z = points3d(m)(2);
	xt = x/(1.0-z); yt = y; zt = z;
	
	// face z = 1
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      Real_wp vloc = gauss_x.EvaluatePhi(j, xt)*gauss_y.EvaluatePhi(k, yt);
	      int col = NumNodes2D_quad(j, k);
	      ShSurfTranspose(1)(i, col) = vloc;
	      ShSurf(1)(col, i) = vloc;
	    }
      }    
  }

  
  void WedgeDgOrtho::ConstructFunctions()
  {
        
    NumOrtho3D.Reallocate(order+1, order+1, order+1); NumOrtho3D.Fill(-1);
    
    bool renum_ortho = false;
    if (order <= 10)
      {
	ifstream file_in("src/FiniteElement/Wedge/ReorderingWedgeOrtho.txt");
	string ligne, line_search;
	line_search = string("Reordering  of  the  dof  ")
	  + "for  the  reference  wedge  of  order  p=" + Seldon::to_str(order);
	
	if (!file_in.is_open())
	  {
	    // cout<<" File containing reordering not found"<<endl;
	  }
	else
	  {
	    while (!file_in.eof())
	      {
		getline(file_in, ligne);
		DeleteSpaceAtExtremityOfString(ligne);
		
		if (!ligne.compare(line_search))
		  {
		    renum_ortho = true;
		    int nb;
		    for (int i = 0; i <= order; i++)
		      for (int j = 0; j <= order-i; j++)
			for (int k = 0; k <= order; k++)
			  {
			    file_in >> nb;
			    NumOrtho3D(i, j, k) = nb-1;
			  }
		  }
	      }
	  }
      }
    
    if (!renum_ortho)
      {
	int nb = 0;
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order-i; j++)
	    for (int k = 0; k <= order; k++)
	      {
		NumOrtho3D(i, j, k) = nb++;
	      }
      }

    nb_dof_loc = points_nodal3d.GetM();
    nb_dof_boundaries = 0;
    points_dof3d = points_nodal3d;
  }
  
  
  //! construction of mass matrix
  void WedgeDgOrtho::ComputeMassMatrix(Matrix<Real_wp, Symmetric, SymColSkyLine> & A,
                                       const VectReal_wp & CoefJacobian) const
  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> Mm;
    ComputeMassMatrix(Mm, CoefJacobian);    
    
    Copy(Mm, A);
  }
  
  
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
  
  
  void WedgeDgOrtho::ComputeMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & Mm,
                                       const VectReal_wp & CoefJacobian) const
  {
    PerformComputationMassMatrix(Mm, CoefJacobian);
    for (int i = 0; i < Mm.GetM(); i++)
      for (int j = i+1; j < Mm.GetM(); j++)
        Mm(i, j) *= 0.5;        
  }

  
  //! construction of mass matrix	
  template<class Matrix1>
  void WedgeDgOrtho
  ::PerformComputationMassMatrix(Matrix1 & Mm, const VectReal_wp & CoefJacobian) const
  {
    Mm.Reallocate(nb_dof_loc,nb_dof_loc); Mm.Fill(0);
    
    Real_wp A1 = CoefJacobian(0);
    Real_wp B1 = CoefJacobian(1), B2 = CoefJacobian(2), B3 = CoefJacobian(3);
    Real_wp C1 = CoefJacobian(4), C2 = CoefJacobian(5);
    Real_wp D1 = CoefJacobian(6);
    
    int i;
    
    for	(int i3 = 0; i3 <= order; i3++)
      for (int i2 = 0; i2 <= order; i2++)
        for (int i1 = 0; i1 <= order-i2; i1++)
          { 
            i = NumOrtho3D(i1,i2,i3);
	    
            Mm(i,i) += A1 + B2/2*(1+alpha(i2,2*i1+1)) + B3/2 + C2/4*(1+alpha(i2,2*i1+1)) + D1/4;
            if (i3>0)
              Mm(i,i) += D1/4*gamma(i3-1,0)*beta(i3,0);
            if (i3<order)
              Mm(i,i ) += D1/4*gamma(i3,0)*beta(i3+1,0);
	    
            if(i2>0)
              Mm(i,NumOrtho3D(i1,i2-1,i3)) += B2/2*beta(i2,2*i1+1) + C2/4*beta(i2,2*i1+1);  
            if (i2<order-i1)
              Mm(i,NumOrtho3D(i1,i2+1,i3)) += B2/2*gamma(i2,2*i1+1) + C2/4*gamma(i2,2*i1+1);
	    
            if (i3<order)
              {
                Mm(i,NumOrtho3D(i1,i2,i3+1))
                  += B3/2*gamma(i3,0) + C2/4*(1+alpha(i2,2*i1+1))*gamma(i3,0) + D1/2*gamma(i3,0);
                if (i2>0)
                  Mm(i,NumOrtho3D(i1,i2-1,i3+1)) += C2/4*beta(i2,2*i1+1)*gamma(i3,0);  
                if (i2<order-i1)
                  Mm(i,NumOrtho3D(i1,i2+1,i3+1)) += C2/4*gamma(i2,2*i1+1)*gamma(i3,0);
              }
            if (i3>0)
              {
                Mm(i,NumOrtho3D(i1,i2,i3-1))
                  += B3/2*beta(i3,0) + C2/4*(1+alpha(i2,2*i1+1))*beta(i3,0) + D1/2*beta(i3,0);
                if (i2>0)
                  Mm(i,NumOrtho3D(i1,i2-1,i3-1)) += C2/4*beta(i2,2*i1+1)*beta(i3,0);
                if (i2<order-i1)
                  Mm(i,NumOrtho3D(i1,i2+1,i3-1)) += C2/4*gamma(i2,2*i1+1)*beta(i3,0);
              }
	    
            if (i3<order-1)
              Mm(i,NumOrtho3D(i1,i2,i3+2)) += D1/4*gamma(i3,0)*gamma(i3+1,0);
            if (i3>1)
              Mm(i,NumOrtho3D(i1,i2,i3-2)) += D1/4*beta(i3,0)*beta(i3-1,0);
            
            for (int j2 = 0; j2 <= order-i1; j2++)
              {	
                Mm(i,NumOrtho3D(i1,j2,i3)) += B1/2*C(i1,i2,j2) + C1/4*C(i1,i2,j2);
		
                if (i3<order)
                  Mm(i,NumOrtho3D(i1,j2,i3+1)) += C1/4*gamma(i3,0)*C(i1,i2,j2);
                if (i3>0)
                  Mm(i,NumOrtho3D(i1,j2,i3-1)) += C1/4*beta(i3,0)*C(i1,i2,j2);
		
              }
	    
            for (int j2 = 0; j2 <= order+1-i1; j2++)
              if (i1>0)
                {
                  Mm(i,NumOrtho3D(i1-1,j2,i3))
                    += B1/2*beta(i1,0)*CmO(i1,i2,j2) + C1/4*beta(i1,0)*CmO(i1,i2,j2);
                  if (i3<order)
                    Mm(i,NumOrtho3D(i1-1,j2,i3+1)) += C1/4*beta(i1,0)*gamma(i3,0)*CmO(i1,i2,j2);
                  if (i3>0)
                    Mm(i,NumOrtho3D(i1-1,j2,i3-1)) += C1/4*beta(i1,0)*beta(i3,0)*CmO(i1,i2,j2);
                }
            
            for (int j2 = 0; j2 <= order-1-i1; j2++)
              {
                Mm(i,NumOrtho3D(i1+1,j2,i3))
                  += B1/2*gamma(i1,0)*CpO(i1,i2,j2) + C1/4*gamma(i1,0)*CpO(i1,i2,j2);
		
                if (i3<order)
                  Mm(i,NumOrtho3D(i1+1,j2,i3+1)) += C1/4*gamma(i1,0)*gamma(i3,0)*CpO(i1,i2,j2);
                if (i3>0)
                  Mm(i,NumOrtho3D(i1+1,j2,i3-1)) += C1/4*gamma(i1,0)*beta(i3,0)*CpO(i1,i2,j2);
              }
          }    
  }
    
  
  template<class Vector1>
  inline void WedgeDgOrtho::SolveMassMatrix(Vector1& x) const
  {
    // mass matrix equal to identity, we don't change x
  }
  
  
  template<class Vector1>
  inline void WedgeDgOrtho::MltMassMatrix(Vector1& y) const
  {
    // mass matrix equal to identity
  }
  
  
  inline Real_wp WedgeDgOrtho::GetMassMatrix(int i, int j) const
  {
    if (i == j)
      return Real_wp(1);
    
    return Real_wp(0);
  }
  
  
  template<class Vector1, class Vector2>
  void WedgeDgOrtho::ApplyCh(const Vector1& U, Vector2& V) const
  {
    //if (!UseQuadraturePointsForSh())
    //return WedgeReference::ApplyCh(U, V);
    
    Vector1 Ux(ChX.GetM()), Uy(ChY.GetM());
    Ux.Fill(0); Uy.Fill(0);
    Mlt(ChX, U, Ux);
    Mlt(ChY, Ux, Uy);
    Mlt(ChZ, Uy, V);
  }
   
  
  template<class Vector1, class Vector2>
  void WedgeDgOrtho::ApplyChTranspose(const Vector1& U, Vector2& V) const
  {
    //if (!UseQuadraturePointsForSh())
    //return WedgeReference::ApplyChTranspose(U, V);
    
    Vector1 Uy(ChTransposeY.GetM()), Uz(ChTransposeZ.GetM());
    Uy.Fill(0); Uz.Fill(0);
    Mlt(ChTransposeZ, U, Uz);
    Mlt(ChTransposeY, Uz, Uy);
    Mlt(ChTransposeX, Uy, V);    
    //DISP(V);
  }
  
  
  template<class Vector1, class Vector2>
  void WedgeDgOrtho::ApplyRhQuadrature(Vector1& U, Vector2& V) const
  {
    // applying chaine rule to get derivatives on unit tetrahedron
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	U(3*i) = DerivDxtildeDx(i)*U(3*i) + DerivDxtildeDy(i)*U(3*i+1);
      }
    
    Mlt(RhLoc, U, V);
  }
  
  
  template<class Vector1, class Vector2>
  void WedgeDgOrtho::ApplyRhQuadratureTranspose(const Vector1& U, Vector2& V) const
  {    
    Mlt(RhLocTranspose, U, V);
    // applying chaine rule to get derivatives on unit tetrahedron
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	V(3*i+1) += V(3*i)*DerivDxtildeDy(i);
	V(3*i) *= DerivDxtildeDx(i);
      }
  }
  
  
  template<class T0, class Vector1, class Vector2>
  void WedgeDgOrtho::ApplyShQuadrature(const T0& alpha, const Vector1& U, Vector2& V) const
  {
    int Nquad = ShSurf(0).GetM();
    int Ntri = ShSurf(0).GetN();
    Vector1 Utri(Ntri), Uquad(Nquad), Uquad_bis(Nquad);
    Uquad.Fill(0); Uquad_bis.Fill(0);
    Vector1 Usurf(5*Nquad);
    for (int i = 0; i < Nquad; i++)
      Usurf(Nquad+i) = U(Ntri + i);
    
    for (int i = 0; i < Nquad; i++)
      Uquad(i) = U(Ntri+Nquad+i);
    
    Mlt(ShQuad, Uquad, Uquad_bis);
    for (int i = 0; i < Nquad; i++)
      Usurf(2*Nquad+i) = Uquad_bis(i);
    
    for (int i = 0; i < Nquad; i++)
      Uquad(i) = U(Ntri+2*Nquad+i);
    
    Mlt(ShQuad, Uquad, Uquad_bis);
    for (int i = 0; i < Nquad; i++)
      Usurf(3*Nquad+i) = Uquad_bis(i);
    
    for (int num_loc = 0; num_loc <= 4; num_loc+=4)
      {
	int offset = 0;
	if (num_loc == 4)
	  offset = Ntri + 3*Nquad;
	
	for (int i = 0; i < Ntri; i++)
	  Utri(i) = U(offset+i);
	
	Mlt(ShSurf(0), Utri, Uquad);
	offset = Nquad*num_loc;
	for (int i = 0; i < Nquad; i++)
	  Usurf(offset+i) = Uquad(i);
      }
    
    Seldon::MltAdd(alpha, ShLoc, Usurf, T0(1), V);
  }
  
  
  template<class Vector1, class Vector2>
  void WedgeDgOrtho::ApplyShQuadratureTranspose(const Vector1& U, Vector2& V) const
  {
    Vector1 Usurf(ShLocTranspose.GetM()); Usurf.Fill(0);
    Mlt(ShLocTranspose, U, Usurf);
    
    int Nquad = ShSurf(0).GetM();
    int Ntri = ShSurf(0).GetN();
    Vector1 Utri(Ntri), Uquad(Nquad), Uquad_bis(Nquad);
    Utri.Fill(0); Uquad_bis.Fill(0);
    for (int i = 0; i < Nquad; i++)
      V(Ntri + i) = Usurf(Nquad+i);
    
    for (int i = 0; i < Nquad; i++)
      Uquad(i) = Usurf(2*Nquad+i);
    
    Mlt(ShQuadTranspose, Uquad, Uquad_bis);
    for (int i = 0; i < Nquad; i++)
      V(Ntri+Nquad+i) = Uquad_bis(i);
    
    for (int i = 0; i < Nquad; i++)
      Uquad(i) = Usurf(3*Nquad+i);
    
    Mlt(ShQuadTranspose, Uquad, Uquad_bis);
    for (int i = 0; i < Nquad; i++)
      V(Ntri+2*Nquad+i) = Uquad_bis(i);
    
    for (int num_loc = 0; num_loc <= 4; num_loc+=4)
      {
	int offset = Nquad*num_loc;	
	for (int i = 0; i < Nquad; i++)
	  Uquad(i) = Usurf(offset+i);
	
	Mlt(ShSurfTranspose(0), Uquad, Utri);
	offset = 0;
	if (num_loc == 4)
	  offset = Ntri + 3*Nquad;
        
	for (int i = 0; i < Ntri; i++)
	  V(offset+i) = Utri(i);
      }
    // DISP(V);
  }
  
  
  void WedgeDgOrtho::ComputeValuesPhiRef(const R3& pointloc, VectReal_wp& phi) const
  {
    WedgeReference::ComputeValuesPhiOrthoRef(pointloc, phi);
  }
  
  
  void WedgeDgOrtho::ComputeGradientPhiRef(const R3& pointloc, VectR3& grad_phi) const
  {
    WedgeReference::ComputeGradientPhiOrthoRef(pointloc, grad_phi);
  }
  
  
  //! returns \f$ \varphi_i(\xi_j) \f$
  /*! 
    num_dof = i : degree of freedom number
    num_points = j : quadrature point number
  */
  inline Real_wp WedgeDgOrtho::GetValuePhi(int num_dof, int num_point) const
  {
    return Value_Phi(num_dof, num_point);
  }
  
  
  //! returns \f$ \varphi^{2D}_k(\xi_j) \f$
  inline Real_wp WedgeDgOrtho::GetValuePhiOnBoundary(int node, int num_point, int num_loc) const
  {
    return Value_Phi(node, this->num_quad_points_surf(num_loc)(num_point));
  }

  
  //! computes \f$ \nabla \varphi_{numdof} (\xi_{numpoint}) \f$
  inline const R3& WedgeDgOrtho::
  GetGradientPhi(int num_dof, int num_point) const
  {
    return Gradient_Phi(num_dof, num_point);
  }
    
  
  //! returns \f$ \varphi_k(\xi_j) \f$
  /*!
    where \f$ \varphi \mbox{ and } \xi \f$
    are restriction of functions on boundary and quadrature point on the boundary
  */
  inline const R3& WedgeDgOrtho::
  GetGradientPhiOnBoundary(int node, int num_point, int num_loc) const
  {
    return GetGradientPhi(node, this->num_quad_points_surf(num_loc)(num_point));
  }
  
  
  template<class Vector1, class Vector2>
  void WedgeDgOrtho::ComputeProjectionDofRef(const Vector1& feval, Vector2& contrib) const
  {
    contrib.Reallocate(nb_dof_loc);
    contrib.Zero();
    //MltAdd(1.0, SeldonTrans, InverseVDM, feval, 0.0, contrib); 
    abort();
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \hat{\varphi}_i \f$
  /*!
    \param[in] feval values of the function f on quadrature points
    \param[out] res integral agains each basis function
    \warning feval(i) is assumed to be already multiplied by the integration weight
   */
  template<class Vector1,class Vector2>
  void WedgeDgOrtho::ComputeIntegralRef(const Vector1 & feval, Vector2& res) const
  {
    res.Reallocate(nb_dof_loc);
    res.Fill(0);
    for (int j = 0; j < nb_dof_loc; j++)
      for (int i = 0; i < nb_points_quadrature_inside; i++)
	res(j) += feval(i)*GetValuePhi(j,i);
    
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \hat{\nabla} \hat{\varphi} \f$
  /*!
    \param[in] feval values of the function f on quadrature points
    \param[out] res integral against gradient of each basis function
    \warning feval(i) is assumed to be already multiplied by the integration weight
   */
  template<class Vector1,class Vector2>
  void WedgeDgOrtho::
  ComputeIntegralGradientRef(const Vector1 & feval, Vector2& res) const
  {
    res.Reallocate(nb_dof_loc);
    res.Fill(0);
    for (int j = 0; j < nb_dof_loc; j++)
      for (int i = 0; i < nb_points_quadrature_inside; i++)
	res(j) += DotProd(Gradient_Phi(j,i), feval(i));
  }
  
  
  //! computation of \f$ \int_\Sigma f \phi_i ds \f$
  /*!
    \param[in] feval values of the function f on quadrature points of the face
    \param[out] res integral against each basis function
    \param[in] num_loc local face number
    \warning feval(i) is assumed to be already multiplied by the integration weight
   */
  template<class Vector1,class Vector2>
  void WedgeDgOrtho::
  ComputeIntegralSurfaceRef(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    res.Reallocate(nb_dof_loc);
    res.Fill(0);
    for (int j = 0; j < points_nodal2d_tri.GetM(); j++)
      {
	int num_dof = FacesDof(j, num_loc);
	for (int i = 0; i < points2d_tri.GetM(); i++)
	  res(num_dof) += feval(i)*GetValuePhiOnBoundary(j, i, num_loc);
      }
  }
  
  
  template<class Vector1,class Vector2>
  void WedgeDgOrtho::
  ComputeGaussIntegralSurfaceRef(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    res.Reallocate(nb_dof_loc);
    res.Fill(0);
    for (int j = 0; j < points_nodal2d_tri.GetM(); j++)
      {
	int num_dof = FacesDof(j, num_loc);
	for (int i = 0; i < points2d_tri.GetM(); i++)
	  res(num_dof) += feval(i)*GetValuePhiOnBoundary(j, i, num_loc);
      }
  }
  
  
  //! not impemented
  template<class Vector1,class Vector2>
  void WedgeDgOrtho::
  ComputeIntegralSurfaceGradientRef(const Vector1 & feval, Vector2& res,
				    int num_loc) const
  {
    res.Fill(0);
  }

}

#define MONTJOIE_FILE_WEDGE_DG_ORTHO_CXX
#endif
