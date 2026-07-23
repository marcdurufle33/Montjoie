#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_ELASTIC_CXX

namespace Montjoie
{

  /*****************************
   * VarTransmission_Elas_Base *
   *****************************/
  
  //! initialization before adding transmission terms
  template<class Dimension>
  void VarTransmission_Elas_Base<Dimension>::InitTransmission()
  {
    if (this->transmission_references.GetM() <= 0)
      return;
    
    var_boundary.ProcRowNumbers_Impedance.Reallocate(var_problem.GetNbDof());
    var_boundary.ProcColumnNumbers_Impedance.Reallocate(var_problem.GetNbDof());
    var_boundary.ProcRowNumbers_Impedance.Fill(0);
    var_boundary.ProcColumnNumbers_Impedance.Fill(0);
    
    var_boundary.NewRowNumbers_Impedance.Reallocate(var_problem.GetNbDof());
    var_boundary.NewColumnNumbers_Impedance.Reallocate(var_problem.GetNbDof());
    var_boundary.NewRowNumbers_Impedance.Fill(-1);
    var_boundary.NewColumnNumbers_Impedance.Fill(-1);

    // interpolation of DF_i from nodal points to quadrature points
    int type_quad = var_problem.GetReferenceElement(0).GetTypeIntegrationEdge();
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    int r = mesh_num.GetOrder();
    MeshInterpolation<Dimension> mesh_interp;
    mesh_interp.SetGaussQuadrature(r, type_quad);
    
    typedef typename Dimension::DimensionBoundary DimensionB;
    int nb_mat = 0;
    for (int i = 0; i < this->MatDfjNodes.GetM(); i++)
      for (int n = 0; n < this->MatDfjNodes(i).GetM(); n++)
        if (this->MatDfjNodes(i)(n).GetM() > 0)
	  nb_mat++;
    
    Vector<const ElementGeomReference<DimensionB>* > elt_geom_surf(nb_mat);
    nb_mat = 0;
    for (int i = 0; i < this->MatDfjNodes.GetM(); i++)
      for (int n = 0; n < this->MatDfjNodes(i).GetM(); n++)
        if (this->MatDfjNodes(i)(n).GetM() > 0)
          {
	    const ElementGeomReference<Dimension>& Fb_geom = var_problem.GetShapeElement(i);
	    elt_geom_surf(nb_mat) = &Fb_geom.GetSurfaceFiniteElement(n);
	    nb_mat++;
	  }
    
    FiniteElementInterpolator proj_surf;
    proj_surf.InitProjection(elt_geom_surf, mesh_interp.GetSubdivisionStep(),
                             mesh_interp.PointsReferenceSurface());
    
    for (int i = 0; i < this->MatDfjNodes.GetM(); i++)
      for (int n = 0; n < this->MatDfjNodes(i).GetM(); n++)
        if (this->MatDfjNodes(i)(n).GetM() > 0)
          {
            Vector<typename Dimension::MatrixN_N> mat_quad;
	    int nf = var_problem.mesh.Element(i).numBoundary(n);
            int type_elt = var_problem.mesh.Boundary(nf).GetHybridType();
	    proj_surf.Project(this->MatDfjNodes(i)(n), mat_quad, type_elt);
            this->MatDfjNodes(i)(n) = mat_quad;
	  } 
  }


  //! adding terms due to transmission conditions
  /*!
    Transmission conditions (developments are performed in the thesis of Alienor Burel).
    The following transmission conditions are implemented
    delta < t(u) > = A [u] + delta B J \partial_1 <u>
    [t(u)] = delta ( rho_i \partial_t^2 - J A J \partial_1^2) <u> - J B \partial_1 [u]
    where
    t(u) = C grad(u) n 
    A = (mu_i 0; 0 lambda_i + 2 mu_i), B = (mu_i 0; 0 lambda_i), J = (0 1; 1 0)
    [u] = u+ - u-
    <u> = 1/2 (u+  + u-)
    
    The variational formulation associated with these transmission conditions is equal to :
    - \int_\Gamma t(u)- v- dx + \int_\Gamma t(u)+  v+ dx
    = \int_\Gamma <t(u)> [v] + [t(u)] <v> dx
    = \delta rho_i \partial_t^2 \int_\Gamma <u> <v> dx
    \delta \int_\Gamma J A J <\partial_1 u> <\partial_1 v> dx
    - \int_\Gamma BJ <u> \partial_1 [v] dx 
    - \int_\Gamma BJ <v> \partial_1 [u] dx 
    + 1/\delta \int_\Gamma A [u] [v] dx
   */
  template<> template<class Complexe>
  void VarTransmission_Elas_Base<Dimension2>
  ::AddTransmissionTermGen(const Complexe& rho_layer, const Complexe& lambda_layer,
			   const Complexe& mu_layer, const Complexe& sigma_layer,
			   const Complexe& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			   VirtualMatrix<Complexe>& mat_sp, int offset_row, int offset_col)
  {
    if (this->transmission_references.GetM() <= 0)
      return;

    const MeshNumbering<Dimension2>& mesh_num = var_problem.GetMeshNumbering(0);
    Complexe a0, b0, a1, b1, m_omega2, m_iomega;

    // diagonal of A and B
    a0 = mu_layer; a1 = lambda_layer+2.0*mu_layer;
    b0 = mu_layer; b1 = lambda_layer;
    var_problem.GetMiomega(m_iomega);
    var_problem.GetMomega2(m_omega2); // 1 for real numbers, -omega^2 for complex numbers
    
    TinyMatrix<Complexe, General, 2, 2> coef_u_phi, coef_jump;
    
    // formulation without Lagrange multiplier
    this->SetModifiedColNumbers(0);
    this->SetModifiedRowNumbers(0);

    bool ref_case = false;
    
    // DISP(ref_case);
    if (ref_case)
      {
	// on recupere les ddls de Gamma
	int r = mesh_num.GetOrder();
	
	int nb_edges = this->NumFaceMinus.GetM();
	Vector<IVect> NodleMinus(nb_edges), NodlePlus(nb_edges);
	int Nvol = mesh_num.GetNbDof();
	const ElementReference<Dimension2, 1>& Fb = var_problem.GetReferenceElementH1(0);
	for (int i = 0; i < this->NumFaceMinus.GetM(); i++)
	  {
	    int nf_m = this->NumFaceMinus(i);
	    int ne_m = var_problem.mesh.BoundaryRef(nf_m).numElement(0);
	    int num_loc_m = var_problem.mesh.Element(ne_m).GetPositionBoundary(nf_m);
	    
	    bool same_rot = var_problem.mesh.Element(ne_m).GetOrientationEdge(num_loc_m);
	    NodleMinus(i).Reallocate(2*(r+1));
	    for (int j = 0; j <= r; j++)
	      {
		int num_dof_loc = Fb.GetLocalNumber(num_loc_m, j);
		if (!same_rot)
		  num_dof_loc = Fb.GetLocalNumber(num_loc_m, r-j);
		
		int num_dof = mesh_num.Element(ne_m).GetNumberDof(num_dof_loc);
		NodleMinus(i)(2*j) = num_dof;
		NodleMinus(i)(2*j+1) = num_dof + Nvol;
	      }
	    
	    int nf_p = this->NumFacePlus(i);
	    int ne_p = var_problem.mesh.BoundaryRef(nf_p).numElement(0);
	    int num_loc_p = var_problem.mesh.Element(ne_p).GetPositionBoundary(nf_p);
	    
	    same_rot = var_problem.mesh.Element(ne_p).GetOrientationEdge(num_loc_p);
	    NodlePlus(i).Reallocate(2*(r+1));
	    for (int j = 0; j <= r; j++)
	      {
		int num_dof_loc = Fb.GetLocalNumber(num_loc_p, j);
		if (!same_rot)
		  num_dof_loc = Fb.GetLocalNumber(num_loc_p, r-j);
		
		int num_dof = mesh_num.Element(ne_p).GetNumberDof(num_dof_loc);
		NodlePlus(i)(2*j) = num_dof;
		NodlePlus(i)(2*j+1) = num_dof + Nvol;
	      }
	    
	    //DISP(i);
	    //DISP(var.mesh.BoundaryRef(nf_m)); DISP(var.mesh.BoundaryRef(nf_p));
	    //DISP(NodleMinus(i)); DISP(NodlePlus(i));
	  }	
	
	// calcul des fonctions de base
	const ElementReference<Dimension1, 1>& edge = Fb.GetSurfaceFiniteElement(0);
	//DISP(edge.Points);
	Matrix<Real_wp> ValPhi(r+1, r+1), GradPhi(r+1, r+1);
	for (int k = 0; k <= r; k++)
	  {
	    for (int j = 0; j <= r; j++)
	      {
		ValPhi(j, k) = edge.GetValuePhi1D(j, edge.Points(k));
		GradPhi(j, k) = edge.GetGradientPhi1D(j, edge.Points(k));
	      }
	  }
	
	//DISP(ValPhi); DISP(GradPhi);
	Vector<Matrix<R2> > ValPhiGlob(nb_edges), DerPhiGlob(nb_edges);
	for (int i = 0; i < nb_edges; i++)
	  {
	    int nf_m = this->NumFaceMinus(i);
	    int ne_m = var_problem.mesh.BoundaryRef(nf_m).numElement(0);
	    int num_loc_m = var_problem.mesh.Element(ne_m).GetPositionBoundary(nf_m);
	    
	    bool same_rot = var_problem.mesh.Element(ne_m).GetOrientationEdge(num_loc_m);    

	    ValPhiGlob(i).Reallocate(2*(r+1), r+1);
	    DerPhiGlob(i).Reallocate(2*(r+1), r+1);
	    R2 normale, tau; Real_wp ds;
	    for (int k = 0; k <= r; k++)
	      {
		// on recupere tau et n
		ds = this->DsNodes(i)(k);
		normale = this->NormaleNodes(i)(k);
		tau.Init(normale(1), -normale(0));
		
		for (int j = 0; j <= r; j++)
		  {
		    ValPhiGlob(i)(2*j, k)(0) = ValPhi(j, k)*tau(0);
		    ValPhiGlob(i)(2*j, k)(1) = ValPhi(j, k)*normale(0);
		    
		    ValPhiGlob(i)(2*j+1, k)(0) = ValPhi(j, k)*tau(1);
		    ValPhiGlob(i)(2*j+1, k)(1) = ValPhi(j, k)*normale(1);
		    
		    DerPhiGlob(i)(2*j, k)(0) = GradPhi(j, k)/ds*tau(0);
		    DerPhiGlob(i)(2*j, k)(1) = GradPhi(j, k)/ds*normale(0);
		    
		    DerPhiGlob(i)(2*j+1, k)(0) = GradPhi(j, k)/ds*tau(1);
		    DerPhiGlob(i)(2*j+1, k)(1) = GradPhi(j, k)/ds*normale(1);
		    
		    if (!same_rot)
		      {
			Mlt(-1.0, DerPhiGlob(i)(2*j, k));
			Mlt(-1.0, DerPhiGlob(i)(2*j+1, k));
		      }		    
		    
		    //DISP(i); DISP(j); DISP(k); DISP(DerPhiGlob(i)(2*j, k));
		    //DISP(DerPhiGlob(i)(2*j+1, k));
		  }
	      }
	  }
	
	// adding terms into the matrix
	// loop on edges of Gamma+
	Matrix<Complexe> mat_PP(2*(r+1), 2*(r+1));
	Matrix<Complexe> mat_MM(2*(r+1), 2*(r+1));
	Matrix<Complexe> mat_PM(2*(r+1), 2*(r+1));
	Matrix<Complexe> mat_MP(2*(r+1), 2*(r+1));
	for (int n = 0; n < nb_edges; n++)
	  {
	    mat_PP.Fill(0); mat_MM.Fill(0);
	    mat_PM.Fill(0); mat_MP.Fill(0);
	    Real_wp ds; Complexe vloc;
	    TinyMatrix<Complexe, General, 2, 2> A;
	    TinyVector<Complexe, 2> vec_u;
	    for (int k = 0; k <= r; k++)
	      {
		ds = this->DsNodes(n)(k)*edge.Weights(k);
		
		// term - \eta \rho \omega^2 <u> <v>
		Complexe coef = (rho_layer*m_omega2*nat_mat.GetCoefMass()
				 + sigma_layer*m_iomega*nat_mat.GetCoefDamping());
		coef *= 0.25*this->delta*ds;
		
		for (int i = 0; i <= 2*r+1; i++)
		  for (int j = 0; j <= 2*r+1; j++)
		    {
		      vloc = DotProd(ValPhiGlob(n)(i, k), ValPhiGlob(n)(j, k))*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) += vloc;
		      mat_MP(i, j) += vloc;
		      mat_PM(i, j) += vloc;
		    }
		
		// term 1/eta A [u] [v]
		A.Fill(0); A(0, 0) = a0; A(1, 1) = a1;
		coef = ds*coef_ch/this->delta*nat_mat.GetCoefStiffness();
		
		for (int i = 0; i <= 2*r+1; i++)
		  for (int j = 0; j <= 2*r+1; j++)
		    {
		      Mlt(A, ValPhiGlob(n)(j, k), vec_u);
		      vloc = DotProd(ValPhiGlob(n)(i, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) += vloc;
		      mat_MP(i, j) -= vloc;
		      mat_PM(i, j) -= vloc;
		    }
		
		// term (0 mu, lambda 0) <du/ds> [v]
		A.Fill(0); A(0, 1) = mu_layer; A(1, 0) = lambda_layer;
		coef = -0.5*ds*coef_sh*nat_mat.GetCoefStiffness();
		for (int i = 0; i <= 2*r+1; i++)
		  for (int j = 0; j <= 2*r+1; j++)
		    {
		      Mlt(A, DerPhiGlob(n)(j, k), vec_u);
		      vloc = DotProd(ValPhiGlob(n)(i, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) -= vloc;
		      mat_MP(i, j) += vloc;
		      mat_PM(i, j) -= vloc;
		    }
		
		// term (0 mu, lambda 0) <dv/ds> [u]
		for (int i = 0; i <= 2*r+1; i++)
		  for (int j = 0; j <= 2*r+1; j++)
		    {
		      Mlt(A, DerPhiGlob(n)(i, k), vec_u);
		      vloc = DotProd(ValPhiGlob(n)(j, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) -= vloc;
		      mat_MP(i, j) -= vloc;
		      mat_PM(i, j) += vloc;
		    }
		
		// term (lambda+2 mu, 0; 0, mu) <du/ds> <dv/ds>
		A.Fill(0); A(0, 0) = a1; A(1, 1) = a0;
		coef = 0.25*ds*coef_sh*this->delta*nat_mat.GetCoefStiffness();
		
		for (int i = 0; i <= 2*r+1; i++)
		  for (int j = 0; j <= 2*r+1; j++)
		    {
		      Mlt(A, DerPhiGlob(n)(j, k), vec_u);
		      vloc = DotProd(DerPhiGlob(n)(i, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) += vloc;
		      mat_MP(i, j) += vloc;
		      mat_PM(i, j) += vloc;
		    }
	      }
	    
	    // on rajoute les matrices elementaires
	    Complexe val;
	    for (int i = 0; i <= 2*r+1; i++)
	      for (int j = 0; j <= 2*r+1; j++)
		{
		  to_complex(mat_PP(i, j)*alpha, val);
		  mat_sp.AddInteraction(NodlePlus(n)(i), NodlePlus(n)(j), val);
		  
		  to_complex(mat_MM(i, j)*alpha, val);
		  mat_sp.AddInteraction(NodleMinus(n)(i), NodleMinus(n)(j), val);
		  
		  to_complex(mat_PM(i, j)*alpha, val);
		  mat_sp.AddInteraction(NodleMinus(n)(i), NodlePlus(n)(j), val);

		  to_complex(mat_MP(i, j)*alpha, val);
		  mat_sp.AddInteraction(NodlePlus(n)(i), NodleMinus(n)(j), val);
		}
	  }
	  
      }
    
    if (ref_case)
      return;
    
    VarBoundaryCondition<Complexe, Dimension2>& var_b 
      = static_cast<VarBoundaryCondition<Complexe, Dimension2>& >(var_boundary);
    
    // coefficients specific to the 2-D case
    TinyMatrix<Complexe, General, 2, 2>
      coef_du_dphi, coef_BJ;
    
    // coef_u_phi = -delta rho omega^2
    coef_u_phi.SetIdentity();
    coef_u_phi *= 0.25*this->delta*(rho_layer*m_omega2*nat_mat.GetCoefMass()
                                    + sigma_layer*m_iomega*nat_mat.GetCoefDamping());
    
    // coef_jump = 1/delta
    coef_jump(0, 0) = a0;
    coef_jump(1, 1) = a1;
    coef_jump *= coef_ch/this->delta*nat_mat.GetCoefStiffness();
    
    // coef_du_dphi = delta JAJ
    coef_du_dphi(0, 0) = 0.25*this->delta*a1;
    coef_du_dphi(1, 1) = 0.25*this->delta*a0;
    
    // coef_u_dphi = B J
    coef_BJ(0, 1) = 0.5*b1;
    coef_BJ(1, 0) = 0.5*b0;
    
    coef_du_dphi *= coef_sh*nat_mat.GetCoefStiffness();
    coef_BJ *= coef_sh*nat_mat.GetCoefStiffness();
    
    IVect ref_plus(var_problem.mesh.GetNbReferences()+1);
    IVect ref_minus(var_problem.mesh.GetNbReferences()+1);
    ref_plus.Fill(0); ref_minus.Fill(0);
    ref_plus(this->RefGammaPlus) = 1; ref_minus(this->RefGammaMoins) = 1;

    // term  - 1/4 delta rho omega^2 u- phi- +  A/delta u- phi-  
    //       + 1/4 delta JAJ d1(u-) d1(phi-) + 1/2 BJ u- d_1 phi- + 1/2 BJ phi- d_1 u-
    ImpedanceTransmission_Elas<Complexe, Dimension2> fimpedance(var_problem, *this);
    fimpedance.curvilinear = curvilinear_case;
    fimpedance.coef_u = coef_u_phi + coef_jump;
    fimpedance.coef_du_dphi = coef_du_dphi;
    fimpedance.coef_u_dphi = -coef_BJ;
    fimpedance.coef_phi_du = -coef_BJ;
    
    var_b.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
				     offset_row, offset_col, fimpedance, false, false);
    
    // term  - 1/4 delta rho omega^2 u+ phi- -  A/delta u+ phi-  
    // term + 1/4 delta JAJ d1(u+) d1(phi-) + 1/2 BJ u+ d_1 phi-  - 1/2 BJ phi- d_1 u+
    fimpedance.coef_u = coef_u_phi - coef_jump;
    fimpedance.coef_u_dphi = coef_BJ;
    fimpedance.coef_phi_du = -coef_BJ;
    
    var_b.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
				     offset_row, offset_col, fimpedance, true, false);
    
    // term  - 1/4 delta rho omega^2 u- phi+ -  A/delta u- phi+
    // term + 1/4 delta JAJ d1(u-) d1(phi+) - 1/2 BJ u- d_1 phi+  +  1/2 BJ phi+ d_1 u-
    fimpedance.coef_u_dphi = -coef_BJ;
    fimpedance.coef_phi_du = coef_BJ;
    
    var_b.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
				     offset_row, offset_col, fimpedance, false, true);
    
    // term  - 1/4 delta rho omega^2 u+ phi+ +  A/delta u+ phi+
    // term + 1/4 delta JAJ d1(u+) d1(phi+) - 1/2 BJ u+ d_1 phi+  -  1/2 BJ phi+ d_1 u+
    ImpedanceTransmission_Elas<Complexe, Dimension2> fimpedance_p(var_problem, *this);
    fimpedance_p.curvilinear = curvilinear_case;
    fimpedance_p.coef_du_dphi = coef_du_dphi;
    fimpedance_p.coef_u_dphi = coef_BJ;
    fimpedance_p.coef_phi_du = coef_BJ;
    fimpedance_p.coef_u = coef_u_phi + coef_jump;
    
    // changement de signe a cause de s qui change de signe de gamma- vers gamma+
    if (curvilinear_case)
      {
	Mlt(-1.0, fimpedance_p.coef_u_dphi);
	Mlt(-1.0, fimpedance_p.coef_phi_du);
      }
    
    var_b.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
				     offset_row, offset_col, fimpedance_p, false, false);
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! transmission terms in 3-D
  template<> template<class Complexe>
  void VarTransmission_Elas_Base<Dimension3>
  ::AddTransmissionTermGen(const Complexe& rho_layer, const Complexe& lambda_layer,
			   const Complexe& mu_layer, const Complexe& sigma_layer,
			   const Complexe& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			   VirtualMatrix<Complexe>& mat_sp, int offset_row, int offset_col)
  {
    if (this->transmission_references.GetM() <= 0)
      return;

    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    Complexe a0, b0, a1, b1, m_omega2, m_iomega;
    
    // diagonal of A and B
    a0 = mu_layer; a1 = lambda_layer+2.0*mu_layer;
    b0 = mu_layer; b1 = lambda_layer;
    var_problem.GetMiomega(m_iomega);
    var_problem.GetMomega2(m_omega2); // 1 for real numbers, -omega^2 for complex numbers
    
    // formulation without Lagrange multiplier
    this->SetModifiedColNumbers(0);
    this->SetModifiedRowNumbers(0);

    bool ref_case = curvilinear_case;
    //bool ref_case = false;
    if (ref_case)
      {
      	// on recupere les ddls de Gamma
	const ElementReference<Dimension3, 1>& Fb = var_problem.GetReferenceElementH1(0);
	
        Vector<R3> PointsDof(mesh_num.GetNbDof());
        Vector<R3> points;
        for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
          {
            var_problem.mesh.GetDofsElement(i, points, Fb.GetGeometricElement());
            for (int j = 0; j < mesh_num.Element(i).GetNbDof(); j++)
              {
                int jglob = mesh_num.Element(i).GetNumberDof(j);
                PointsDof(jglob) = points(j);
              }
          }
        
	int r = Fb.GetOrder();
        Matrix<int> FacesDof_Rotation, NumNodes2D, coor_nodes;
        MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(r, NumNodes2D, coor_nodes);
        MeshNumbering<Dimension3>::GetRotationQuadrilateralFace(NumNodes2D, FacesDof_Rotation);
        //DISP(FacesDof_Rotation);
        
	int nb_faces = this->NumFaceMinus.GetM();
	Vector<IVect> NodleMinus(nb_faces), NodlePlus(nb_faces);
	int Nvol = mesh_num.GetNbDof();
	for (int i = 0; i < this->NumFaceMinus.GetM(); i++)
	  {
	    int nf_m = this->NumFaceMinus(i);
	    int ne_m = var_problem.mesh.BoundaryRef(nf_m).numElement(0);
	    int num_loc_m = var_problem.mesh.Element(ne_m).GetPositionBoundary(nf_m);

	    int rot = var_problem.mesh.Element(ne_m).GetOrientationFace(num_loc_m);
	    NodleMinus(i).Reallocate(3*(r+1)*(r+1));
	    for (int j = 0; j < (r+1)*(r+1); j++)
	      {
                int jrot = FacesDof_Rotation(rot, j);
		int num_dof_loc = Fb.GetLocalNumber(num_loc_m, jrot); 
                
		int num_dof = mesh_num.Element(ne_m).GetNumberDof(num_dof_loc);
		NodleMinus(i)(3*j) = num_dof;
		NodleMinus(i)(3*j+1) = num_dof + Nvol;
                NodleMinus(i)(3*j+2) = num_dof + 2*Nvol;
	      }
	    
	    int nf_p = this->NumFacePlus(i);
	    int ne_p = var_problem.mesh.BoundaryRef(nf_p).numElement(0);
	    int num_loc_p = var_problem.mesh.Element(ne_p).GetPositionBoundary(nf_p);
	    
	    rot = var_problem.mesh.Element(ne_p).GetOrientationFace(num_loc_p);
	    NodlePlus(i).Reallocate(3*(r+1)*(r+1));
            for (int j = 0; j < (r+1)*(r+1); j++)
	      {
                int jrot = FacesDof_Rotation(rot, j);
		int num_dof_loc = Fb.GetLocalNumber(num_loc_p, jrot);
		
		int num_dof = mesh_num.Element(ne_p).GetNumberDof(num_dof_loc);
		NodlePlus(i)(3*j) = num_dof;
		NodlePlus(i)(3*j+1) = num_dof + Nvol;
                NodlePlus(i)(3*j+2) = num_dof + 2*Nvol;
	      }
	    
            //DISP(i); DISP(rot);
	    //DISP(var_problem.mesh.BoundaryRef(nf_m)); DISP(var_problem.mesh.BoundaryRef(nf_p));
	    //DISP(NodleMinus(i)); DISP(NodlePlus(i));

            /* for (int j = 0; j < (r+1)*(r+1); j++)
              {
                if (j < 4)
                  {
                    // we check vertices of the face
                    int nA = var_problem.mesh.BoundaryRef(nf_m).numVertex(j);
                    R3 ptA = var_problem.mesh.Vertex(nA);
                    if (PointsDof(NodleMinus(i)(3*j)) != ptA)
                      {
                        cout << "Dof incorrect" << endl;
                        abort();
                      }
                    
                    nA = var_problem.mesh.BoundaryRef(nf_p).numVertex(j);
                    ptA = var_problem.mesh.Vertex(nA);
                    if (PointsDof(NodlePlus(i)(3*j)) != ptA)
                      {
                        cout << "Dof incorrect" << endl;
                        abort();
                      }                                        
                  }
                
                if (PointsDof(NodleMinus(i)(3*j)) != PointsDof(NodlePlus(i)(3*j)))
                  {
                    cout << "Dof incorrect" << endl;
                    abort();
                  }
              }
            */
	  }	
	
	// calcul des fonctions de base
	QuadrangleLobatto edge;
        edge.ConstructFiniteElement(r);
	//DISP(edge.Points);
	Matrix<Real_wp> ValPhi((r+1)*(r+1), (r+1)*(r+1));
        Matrix<R2> GradPhi((r+1)*(r+1), (r+1)*(r+1));
        for (int k = 0; k < edge.GetNbDof(); k++)
	  {
            VectReal_wp val_phi; VectR2 grad_phi;
            edge.ComputeValuesPhiRef(edge.PointsND(k), val_phi);
            edge.ComputeGradientPhiRef(edge.PointsND(k), grad_phi);
	    for (int j = 0; j < edge.GetNbDof(); j++)
	      {
		ValPhi(j, k) = val_phi(j);
		GradPhi(j, k) = grad_phi(j);
	      }
	  }
	
	//DISP(ValPhi); DISP(GradPhi);
	Vector<Matrix<R3> > ValPhiGlob(nb_faces), DerPhiGlobX(nb_faces), DerPhiGlobY(nb_faces);
	for (int i = 0; i < nb_faces; i++)
	  {
	    //int nf_m = this->NumFaceMinus(i);
	    //int ne_m = var_problem.mesh.BoundaryRef(nf_m).numElement(0);
	    //int num_loc_m = var_problem.mesh.Element(ne_m).GetPositionBoundary(nf_m);
	    
	    ValPhiGlob(i).Reallocate(3*(r+1)*(r+1), (r+1)*(r+1));
	    DerPhiGlobX(i).Reallocate(3*(r+1)*(r+1), (r+1)*(r+1));
            DerPhiGlobY(i).Reallocate(3*(r+1)*(r+1), (r+1)*(r+1));
	    R3 normale, tauX, tauY, dF_dx, dF_dy, grad_chap, gradP;
            Real_wp dphi1, dphi2; Matrix3_3 mat_dfj, dfjm1;
	    for (int k = 0; k < edge.GetNbDof(); k++)
	      {
                dF_dx.Fill(0); dF_dy.Fill(0);
                for (int j = 0; j < edge.GetNbDof(); j++)
                  {
                    Add(GradPhi(j, k)(0), this->PosNodes(i)(j), dF_dx);
                    Add(GradPhi(j, k)(1), this->PosNodes(i)(j), dF_dy);
                  }
                
                TimesProd(dF_dx, dF_dy, normale);
                Mlt(1.0/Norm2(normale), normale);
                
                // jacobian matrix
                mat_dfj.Fill(0);
                mat_dfj(0, 0) = dF_dx(0); mat_dfj(1, 0) = dF_dx(1); mat_dfj(2, 0) = dF_dx(2);
                mat_dfj(0, 1) = dF_dy(0); mat_dfj(1, 1) = dF_dy(1); mat_dfj(2, 1) = dF_dy(2);
                mat_dfj(0, 2) = normale(0); mat_dfj(1, 2) = normale(1); mat_dfj(2, 2) = normale(2);
                GetInverse(mat_dfj, dfjm1);
                
		// on recupere tau et n
		normale = this->NormaleNodes(i)(k);
                GetVectorPlane(normale, tauX, tauY);
                tauX.Init(1, 0, 0); tauY.Init(0, 1, 0);
                
		for (int j = 0; j < edge.GetNbDof(); j++)
		  {
		    ValPhiGlob(i)(3*j, k)(0) = ValPhi(j, k)*tauX(0);
                    ValPhiGlob(i)(3*j, k)(1) = ValPhi(j, k)*tauY(0);
		    ValPhiGlob(i)(3*j, k)(2) = ValPhi(j, k)*normale(0);
		    
		    ValPhiGlob(i)(3*j+1, k)(0) = ValPhi(j, k)*tauX(1);
                    ValPhiGlob(i)(3*j+1, k)(1) = ValPhi(j, k)*tauY(1);
		    ValPhiGlob(i)(3*j+1, k)(2) = ValPhi(j, k)*normale(1);

		    ValPhiGlob(i)(3*j+2, k)(0) = ValPhi(j, k)*tauX(2);
                    ValPhiGlob(i)(3*j+2, k)(1) = ValPhi(j, k)*tauY(2);
		    ValPhiGlob(i)(3*j+2, k)(2) = ValPhi(j, k)*normale(2);
		    
                    grad_chap.Init(GradPhi(j, k)(0), GradPhi(j, k)(1), 0);
                    MltTrans(dfjm1, grad_chap, gradP);
                    dphi1 = DotProd(gradP, tauX); dphi2 = DotProd(gradP, tauY);
                    
		    DerPhiGlobX(i)(3*j, k)(0) = dphi1*tauX(0);
                    DerPhiGlobX(i)(3*j, k)(1) = dphi1*tauY(0);
		    DerPhiGlobX(i)(3*j, k)(2) = dphi1*normale(0);

		    DerPhiGlobX(i)(3*j+1, k)(0) = dphi1*tauX(1);
                    DerPhiGlobX(i)(3*j+1, k)(1) = dphi1*tauY(1);
		    DerPhiGlobX(i)(3*j+1, k)(2) = dphi1*normale(1);

		    DerPhiGlobX(i)(3*j+2, k)(0) = dphi1*tauX(2);
                    DerPhiGlobX(i)(3*j+2, k)(1) = dphi1*tauY(2);
		    DerPhiGlobX(i)(3*j+2, k)(2) = dphi1*normale(2);
		    
		    DerPhiGlobY(i)(3*j, k)(0) = dphi2*tauX(0);
                    DerPhiGlobY(i)(3*j, k)(1) = dphi2*tauY(0);
		    DerPhiGlobY(i)(3*j, k)(2) = dphi2*normale(0);

		    DerPhiGlobY(i)(3*j+1, k)(0) = dphi2*tauX(1);
                    DerPhiGlobY(i)(3*j+1, k)(1) = dphi2*tauY(1);
		    DerPhiGlobY(i)(3*j+1, k)(2) = dphi2*normale(1);

		    DerPhiGlobY(i)(3*j+2, k)(0) = dphi2*tauX(2);
                    DerPhiGlobY(i)(3*j+2, k)(1) = dphi2*tauY(2);
		    DerPhiGlobY(i)(3*j+2, k)(2) = dphi2*normale(2);

		    //DISP(i); DISP(j); DISP(k); DISP(DerPhiGlob(i)(2*j, k));
		    //DISP(DerPhiGlob(i)(2*j+1, k));
		  }
	      }
	  }
	
	// adding terms into the matrix
	// loop on edges of Gamma+
	Matrix<Complexe> mat_PP(3*(r+1)*(r+1), 3*(r+1)*(r+1));
	Matrix<Complexe> mat_MM(3*(r+1)*(r+1), 3*(r+1)*(r+1));
	Matrix<Complexe> mat_PM(3*(r+1)*(r+1), 3*(r+1)*(r+1));
	Matrix<Complexe> mat_MP(3*(r+1)*(r+1), 3*(r+1)*(r+1));
	for (int n = 0; n < nb_faces; n++)
	  {
	    mat_PP.Fill(0); mat_MM.Fill(0);
	    mat_PM.Fill(0); mat_MP.Fill(0);
	    Real_wp ds; Complexe vloc;
	    TinyMatrix<Complexe, General, 3, 3> A;
	    TinyVector<Complexe, 3> vec_u;
	    for (int k = 0; k < (r+1)*(r+1); k++)
	      {
		ds = this->DsNodes(n)(k)*edge.WeightsND(k);
		
		// term - \eta \rho \omega^2 <u> <v>
		Complexe coef = (rho_layer*m_omega2*nat_mat.GetCoefMass()
				 + sigma_layer*m_iomega*nat_mat.GetCoefDamping());
		coef *= 0.25*this->delta*ds;
                
                for (int i = 0; i < 3*edge.GetNbDof(); i++)
		  for (int j = 0; j < 3*edge.GetNbDof(); j++)
		    {
		      vloc = DotProd(ValPhiGlob(n)(i, k), ValPhiGlob(n)(j, k))*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) += vloc;
		      mat_MP(i, j) += vloc;
		      mat_PM(i, j) += vloc;
		    }
		
		// term 1/eta A [u] [v]
		A.Fill(0); A(0, 0) = a0; A(1, 1) = a0; A(2, 2) = a1;
		coef = ds*coef_ch/this->delta*nat_mat.GetCoefStiffness();
		
		for (int i = 0; i < 3*edge.GetNbDof(); i++)
		  for (int j = 0; j < 3*edge.GetNbDof(); j++)
		    {
		      Mlt(A, ValPhiGlob(n)(j, k), vec_u);
		      vloc = DotProd(ValPhiGlob(n)(i, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) += vloc;
		      mat_MP(i, j) -= vloc;
		      mat_PM(i, j) -= vloc;
                    }
                
		// term (0 0 mu, 0 0 0, lambda 0 0) <du/dx> [v]
                A.Fill(0); A(0, 2) = mu_layer; A(2, 0) = lambda_layer;
		coef = 0.5*ds*coef_sh*nat_mat.GetCoefStiffness();
		for (int i = 0; i < 3*edge.GetNbDof(); i++)
		  for (int j = 0; j < 3*edge.GetNbDof(); j++)
		    {
		      Mlt(A, DerPhiGlobX(n)(j, k), vec_u);
		      vloc = DotProd(ValPhiGlob(n)(i, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) -= vloc;
		      mat_MP(i, j) += vloc;
		      mat_PM(i, j) -= vloc;
		    }
		
		// term (0 0 mu, 0 0 0, lambda 0 0) <dv/dx> [u]
		for (int i = 0; i < 3*edge.GetNbDof(); i++)
		  for (int j = 0; j < 3*edge.GetNbDof(); j++)
		    {
		      Mlt(A, DerPhiGlobX(n)(i, k), vec_u);
		      vloc = DotProd(ValPhiGlob(n)(j, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) -= vloc;
		      mat_MP(i, j) -= vloc;
		      mat_PM(i, j) += vloc;
		    }

		// term (0 0 0, 0 0 mu, 0 lambda 0) <du/dy> [v]
		A.Fill(0); A(1, 2) = mu_layer; A(2, 1) = lambda_layer;
		coef = 0.5*ds*coef_sh*nat_mat.GetCoefStiffness();
		for (int i = 0; i < 3*edge.GetNbDof(); i++)
		  for (int j = 0; j < 3*edge.GetNbDof(); j++)
		    {
		      Mlt(A, DerPhiGlobY(n)(j, k), vec_u);
		      vloc = DotProd(ValPhiGlob(n)(i, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) -= vloc;
		      mat_MP(i, j) += vloc;
		      mat_PM(i, j) -= vloc;
                    }
		
		// term (0 0 0, 0 0 mu, 0 lambda 0) <dv/dy> [u]
		for (int i = 0; i < 3*edge.GetNbDof(); i++)
		  for (int j = 0; j < 3*edge.GetNbDof(); j++)
		    {
		      Mlt(A, DerPhiGlobY(n)(i, k), vec_u);
		      vloc = DotProd(ValPhiGlob(n)(j, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) -= vloc;
		      mat_MP(i, j) -= vloc;
		      mat_PM(i, j) += vloc;
                    }
		
		// term (0 lambda 0, mu 0 0, 0 0 0) <du/dx> <dv/dy>
                A.Fill(0); A(0, 1) = lambda_layer; A(1, 0) = mu_layer;
		coef = 0.25*ds*coef_sh*this->delta*nat_mat.GetCoefStiffness();
		for (int i = 0; i < 3*edge.GetNbDof(); i++)
		  for (int j = 0; j < 3*edge.GetNbDof(); j++)
		    {
		      Mlt(A, DerPhiGlobX(n)(j, k), vec_u);
		      vloc = DotProd(DerPhiGlobY(n)(i, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) += vloc;
		      mat_MP(i, j) += vloc;
		      mat_PM(i, j) += vloc;
		    }
		
		// term (0 lambda 0, mu 0 0, 0 0 0) <dv/dx> <du/dy>
		for (int i = 0; i < 3*edge.GetNbDof(); i++)
		  for (int j = 0; j < 3*edge.GetNbDof(); j++)
		    {
		      Mlt(A, DerPhiGlobX(n)(i, k), vec_u);
		      vloc = DotProd(DerPhiGlobY(n)(j, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) += vloc;
		      mat_MP(i, j) += vloc;
		      mat_PM(i, j) += vloc;
                    }
		
		// term (mu 0 0, 0 lambda+2 mu 0, 0 0 mu) <du/dy> <dv/dy>
		A.Fill(0); A(0, 0) = a0; A(1, 1) = a1; A(2, 2) = a0;
		coef = 0.25*ds*coef_sh*this->delta*nat_mat.GetCoefStiffness();
		
		for (int i = 0; i < 3*edge.GetNbDof(); i++)
		  for (int j = 0; j < 3*edge.GetNbDof(); j++)
		    {
		      Mlt(A, DerPhiGlobY(n)(j, k), vec_u);
		      vloc = DotProd(DerPhiGlobY(n)(i, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) += vloc;
		      mat_MP(i, j) += vloc;
		      mat_PM(i, j) += vloc;
		    }

		// term (lambda + 2 mu 0 0, 0 mu 0, 0 0 mu) <du/dx> <dv/dx>
		A.Fill(0); A(0, 0) = a1; A(1, 1) = a0; A(2, 2) = a0;
		coef = 0.25*ds*coef_sh*this->delta*nat_mat.GetCoefStiffness();
		
		for (int i = 0; i < 3*edge.GetNbDof(); i++)
		  for (int j = 0; j < 3*edge.GetNbDof(); j++)
		    {
		      Mlt(A, DerPhiGlobX(n)(j, k), vec_u);
		      vloc = DotProd(DerPhiGlobX(n)(i, k), vec_u)*coef;
		      mat_PP(i, j) += vloc;
		      mat_MM(i, j) += vloc;
		      mat_MP(i, j) += vloc;
		      mat_PM(i, j) += vloc;
		    }
	      }
	    
	    // on rajoute les matrices elementaires
	    Complexe val;
	    for (int i = 0; i < 3*edge.GetNbDof(); i++)
	      for (int j = 0; j < 3*edge.GetNbDof(); j++)
		{
		  to_complex(mat_PP(i, j)*alpha, val);
		  if (abs(val) > var_problem.GetThresholdMatrix())
                    mat_sp.AddInteraction(NodlePlus(n)(i), NodlePlus(n)(j), val);
		  
		  to_complex(mat_MM(i, j)*alpha, val);
		  if (abs(val) > var_problem.GetThresholdMatrix())
                    mat_sp.AddInteraction(NodleMinus(n)(i), NodleMinus(n)(j), val);
		  
		  to_complex(mat_PM(i, j)*alpha, val);
		  if (abs(val) > var_problem.GetThresholdMatrix())
                    mat_sp.AddInteraction(NodleMinus(n)(i), NodlePlus(n)(j), val);

		  to_complex(mat_MP(i, j)*alpha, val);
		  if (abs(val) > var_problem.GetThresholdMatrix())
                    mat_sp.AddInteraction(NodlePlus(n)(i), NodleMinus(n)(j), val);
		}
          }        
      }
    
    if (ref_case)
      return;

    VarBoundaryCondition<Complexe, Dimension3>& var_b 
      = static_cast<VarBoundaryCondition<Complexe, Dimension3>& >(var_boundary);

    TinyMatrix<Complexe, General, 3, 3> coef_u_phi, coef_jump;
    
    // coefficients specific to the 3-D case
    TinyMatrix<Complexe, General, 3, 3> coef_B2_J2, coef_B1_J1, coef_J3_B3,
      coef_J1_A_J1, coef_J2_A_J2;
    
    // coef_u_phi = -delta rho omega^2
    coef_u_phi.SetIdentity();
    coef_u_phi *= 0.25*this->delta*(rho_layer*m_omega2*nat_mat.GetCoefMass()
                                    + sigma_layer*m_iomega*nat_mat.GetCoefDamping());
    
    // coef_jump = 1/delta A
    coef_jump(0, 0) = a0;
    coef_jump(1, 1) = a0;
    coef_jump(2, 2) = a1;
    coef_jump *= coef_ch/this->delta*nat_mat.GetCoefStiffness();
    
    // other coefficients
    coef_B1_J1(1, 2) = 0.5*b0;
    coef_B1_J1(2, 1) = 0.5*b1;
    coef_B2_J2(0, 2) = 0.5*b0;
    coef_B2_J2(2, 0) = 0.5*b1;
    
    coef_J1_A_J1(0, 0) = 0.25*this->delta*a0;
    coef_J1_A_J1(1, 1) = 0.25*this->delta*a1;
    coef_J1_A_J1(2, 2) = 0.25*this->delta*a0;
    
    coef_J2_A_J2(0, 0) = 0.25*this->delta*a1;
    coef_J2_A_J2(1, 1) = 0.25*this->delta*a0;
    coef_J2_A_J2(2, 2) = 0.25*this->delta*a0;
    
    coef_J3_B3(0, 1) = 0.25*this->delta*b1;
    coef_J3_B3(1, 0) = 0.25*this->delta*b0;
    
    coef_B1_J1 *= coef_sh*nat_mat.GetCoefStiffness();
    coef_B2_J2 *= coef_sh*nat_mat.GetCoefStiffness();
    coef_J1_A_J1 *= coef_sh*nat_mat.GetCoefStiffness();
    coef_J2_A_J2 *= coef_sh*nat_mat.GetCoefStiffness();
    coef_J3_B3 *= coef_sh*nat_mat.GetCoefStiffness();

    IVect ref_plus(var_problem.mesh.GetNbReferences()+1);
    IVect ref_minus(var_problem.mesh.GetNbReferences()+1);
    ref_plus.Fill(0); ref_minus.Fill(0);
    ref_plus(this->RefGammaPlus) = 1; ref_minus(this->RefGammaMoins) = 1;
    
    // term  - 1/4 delta rho omega^2 u- phi- +  A/delta u- phi-  
    //       + 1/4 delta J1 A J1 d1(u-) d1(phi-) + 1/2 B2 J2 u- d_1 phi-
    //       + 1/2 B2 J2 phi- d_1 u-
    //       + 1/4 delta J2 A J2 d2(u-) d2(phi-) + 1/2 B1 J1 u- d_2 phi- 
    //       + 1/2 B1 J1 phi- d_2 u-
    //       + 1/4 delta J3 B3 d1(u-) d2(phi-) + 1/4 delta J3 B3 d1(phi-) d2(u-) 
    ImpedanceTransmission_Elas<Complexe, Dimension3> fimpedance(var_problem, *this);
    fimpedance.curvilinear = curvilinear_case;
    fimpedance.coef_u = coef_u_phi + coef_jump;
    fimpedance.coef_d1u_d1phi = coef_J2_A_J2;
    fimpedance.coef_d2u_d2phi = coef_J1_A_J1;
    fimpedance.coef_u_d1_phi = coef_B2_J2;
    fimpedance.coef_phi_d1_u = coef_B2_J2;
    fimpedance.coef_u_d2_phi = coef_B1_J1;
    fimpedance.coef_phi_d2_u = coef_B1_J1;
    fimpedance.coef_d1u_d2phi = coef_J3_B3;
    fimpedance.coef_d1phi_d2u = coef_J3_B3;
    var_b.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
				     offset_row, offset_col, fimpedance, false, false);
    
    // term  - 1/4 delta rho omega^2 u+ phi- -  A/delta u+ phi-  
    // term + 1/4 delta J1 A J1 d1(u+) d1(phi-) + 1/2 B2 J2 u+ d_1 phi- 
    // - 1/2 B2 J2 phi- d_1 u+
    // term + 1/4 delta J2 A J2 d2(u+) d2(phi-) + 1/2 B1 J1 u+ d_2 phi-  
    // - 1/2 B1 J1 phi- d_2 u+
    //      + 1/4 delta J3 B3 d1(u+) d2(phi-) + 1/4 delta J3 B3 d1(phi-) d2(u+) 
    fimpedance.coef_u = coef_u_phi - coef_jump;
    fimpedance.coef_u_d1_phi = coef_B2_J2;
    fimpedance.coef_phi_d1_u = -coef_B2_J2;
    fimpedance.coef_u_d2_phi = coef_B1_J1;
    fimpedance.coef_phi_d2_u = -coef_B1_J1;
    var_b.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
				     offset_row, offset_col, fimpedance, true, false);
    
    // term  - 1/4 delta rho omega^2 u- phi+ -  A/delta u- phi+
    // term + 1/4 delta J1 A J1 d1(u-) d1(phi+) - 1/2 B2 J2 u- d_1 phi+ 
    //  +  1/2 B2 J2 phi+ d_1 u-
    // term + 1/4 delta J2 A J2 d2(u-) d2(phi+) - 1/2 B1 J1 u- d_2 phi+ 
    //  +  1/2 B1 J1 phi+ d_2 u-
    //      + 1/4 delta J3 B3 d1(u-) d2(phi+) + 1/4 delta J3 B3 d1(phi+) d2(u-) 
    fimpedance.coef_u_d1_phi = -coef_B2_J2;
    fimpedance.coef_phi_d1_u = coef_B2_J2;
    fimpedance.coef_u_d2_phi = -coef_B1_J1;
    fimpedance.coef_phi_d2_u = coef_B1_J1;
    var_b.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
				     offset_row, offset_col, fimpedance, false, true);
    
    ImpedanceTransmission_Elas<Complexe, Dimension3> fimpedance_p(var_problem, *this);
    fimpedance_p.coef_d1u_d1phi = coef_J2_A_J2;
    fimpedance_p.coef_d2u_d2phi = coef_J1_A_J1;
    fimpedance_p.coef_d1u_d2phi = coef_J3_B3;
    fimpedance_p.coef_d1phi_d2u = coef_J3_B3;
        
    // term  - 1/4 delta rho omega^2 u+ phi+ +  A/delta u+ phi+
    // term + 1/4 delta J1 A J1 d1(u+) d1(phi+) - 1/2 B2 J2 u+ d_1 phi+
    // -  1/2 B2 J2 phi+ d_1 u+
    // term + 1/4 delta J2 A J2 d2(u+) d2(phi+) - 1/2 B1 J1 u+ d_2 phi+ 
    // -  1/2 B1 J1 phi+ d_2 u+
    //      + 1/4 delta J3 B3 d1(u+) d2(phi+) + 1/4 delta J3 B3 d1(phi+) d2(u+) 
    fimpedance_p.coef_u_d1_phi = -coef_B2_J2;
    fimpedance_p.coef_phi_d1_u = -coef_B2_J2;
    fimpedance_p.coef_u_d2_phi = -coef_B1_J1;
    fimpedance_p.coef_phi_d2_u = -coef_B1_J1;
    fimpedance_p.coef_u = coef_u_phi + coef_jump;
    
    if (curvilinear_case)
      {
	Mlt(-1.0, fimpedance_p.coef_u_d1_phi);
	Mlt(-1.0, fimpedance_p.coef_phi_d1_u);
	Mlt(-1.0, fimpedance_p.coef_u_d2_phi);
	Mlt(-1.0, fimpedance_p.coef_phi_d2_u);
      }
    
    var_b.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
				     offset_row, offset_col, fimpedance_p, false, false);
  }
#endif

  
  /************************
   * VarTransmission_Elas *
   ************************/


  //! we read data file
  template<class Complexe, class Dimension>
  void VarTransmission_Elas<Complexe, Dimension>::
  SetInputData(const string& description_field, const VectString& parameters)
  {
    VarTransmission_Elas_Base<Dimension>::SetInputData(description_field, parameters);
    
    if (!description_field.compare("ReferenceTransmission"))
      {
	// ReferenceTransmssion ref PARAM rho lambda mu Decalage delta
	// e.g. ReferenceTransmission = 1 PARAM 0.5 2.3 0.8 0.5 0.01
        int pos = -1;
        int ref_max = 0;
        for (int i = 0; i < parameters.GetM(); i++)
          if ( !  parameters(i).compare("PARAM")) 
            pos = i;
        
        if (pos <= 0)
          {
	    cout << "Needs references before keyword PARAM in ReferenceTransmission" << endl;
            abort();
          }
        
        IVect ref_(pos);
        for (int i = 0; i < pos; i++)
          {
            ref_(i) = to_num<int>(parameters(i));
            ref_max = max(ref_max, ref_(i));
          }
        
        if (ref_max > var_problem.mesh.GetNbReferences())
          var_problem.mesh.ResizeNbReferences(ref_max);
        
	Vector<Vector<Complexe> >& Param_condition = var_boundary.GetParamCondition();
        if (ref_max > Param_condition.GetM())
          Param_condition.Resize(ref_max+1);
        
        if (parameters.GetM() > pos)
          {
            // we store these values in parameters_condition
            Vector<Complexe> parameters_condition(parameters.GetM() - pos-1);
            
            for (int i = pos+1; i < parameters.GetM(); i++)
              parameters_condition(i - pos-1) = to_num<Complexe>(parameters(i));
	    
	    if (parameters_condition.GetM() <= 4)
	      {
		cout << "In SetInputData of VarTransmission_Elas" << endl;
		cout << "ReferenceTransmission needs more parameters, for instance :" << endl;
		cout << "ReferenceTransmission = ref PARAM rho lambda mu Decalage delta" << endl;
		cout << "Current parameters are : " << endl << parameters_condition << endl;
		abort();
	      }
	    
            // dans l'ordre, rho, lambda, mu, Decalage, delta
            rho_layer = parameters_condition(0) ;
            lambda_layer = parameters_condition(1) ;
            mu_layer = parameters_condition(2) ;
            this->Decalage = realpart(parameters_condition(3)) ;
            this->delta = realpart(parameters_condition(4));
	    if (parameters_condition.GetM() > 5)
	      sigma_layer = parameters_condition(5);
          }
        else
          {
            abort();
          }
        
        this->transmission_references = ref_;
      }
    else if (!description_field.compare("CurvilinearTransmission"))
      {
	if (parameters(0) == "YES")
	  this->curvilinear_case = true;
      }
  }
  

#ifdef MONTJOIE_WITH_TWO_DIM  
  //! adding transmission terms (for any dimension)
  template<>
  void VarTransmission_Elas<Real_wp, Dimension2>::
  AddTransmissionTerms(const Real_wp& alpha, const GlobalGenericMatrix<Real_wp>& nat_mat,
                       VirtualMatrix<Real_wp>& mat_sp, int offset_row, int offset_col)
  {
    this->AddTransmissionTermGen(rho_layer, lambda_layer, mu_layer, sigma_layer,
				 alpha, nat_mat, mat_sp, offset_row, offset_col);
  }


  //! adding transmission terms (for any dimension)
  template<>
  void VarTransmission_Elas<Real_wp, Dimension2>::
  AddTransmissionTerms(const Complex_wp& alpha, const GlobalGenericMatrix<Complex_wp>& nat_mat,
                       VirtualMatrix<Complex_wp>& mat_sp, int offset_row, int offset_col)
  {
    cout << "not implemented" << endl;
    abort();
  }


  //! adding transmission terms (for any dimension)
  template<>
  void VarTransmission_Elas<Complex_wp, Dimension2>::
  AddTransmissionTerms(const Real_wp& alpha, const GlobalGenericMatrix<Real_wp>& nat_mat,
                       VirtualMatrix<Real_wp>& mat_sp, int offset_row, int offset_col)
  {
    cout << "not possible" << endl;
    abort();
  }


  //! adding transmission terms (for any dimension)
  template<>
  void VarTransmission_Elas<Complex_wp, Dimension2>::
  AddTransmissionTerms(const Complex_wp& alpha, const GlobalGenericMatrix<Complex_wp>& nat_mat,
                       VirtualMatrix<Complex_wp>& mat_sp, int offset_row, int offset_col)
  {
    this->AddTransmissionTermGen(rho_layer, lambda_layer, mu_layer, sigma_layer,
				 alpha, nat_mat, mat_sp, offset_row, offset_col);
  }
#endif


#ifdef MONTJOIE_WITH_THREE_DIM  
  //! adding transmission terms (for any dimension)
  template<>
  void VarTransmission_Elas<Real_wp, Dimension3>::
  AddTransmissionTerms(const Real_wp& alpha, const GlobalGenericMatrix<Real_wp>& nat_mat,
                       VirtualMatrix<Real_wp>& mat_sp, int offset_row, int offset_col)
  {
    this->AddTransmissionTermGen(rho_layer, lambda_layer, mu_layer, sigma_layer,
				 alpha, nat_mat, mat_sp, offset_row, offset_col);
  }


  //! adding transmission terms (for any dimension)
  template<>
  void VarTransmission_Elas<Real_wp, Dimension3>::
  AddTransmissionTerms(const Complex_wp& alpha, const GlobalGenericMatrix<Complex_wp>& nat_mat,
                       VirtualMatrix<Complex_wp>& mat_sp, int offset_row, int offset_col)
  {
    cout << "not implemented" << endl;
    abort();
  }


  //! adding transmission terms (for any dimension)
  template<>
  void VarTransmission_Elas<Complex_wp, Dimension3>::
  AddTransmissionTerms(const Real_wp& alpha, const GlobalGenericMatrix<Real_wp>& nat_mat,
                       VirtualMatrix<Real_wp>& mat_sp, int offset_row, int offset_col)
  {
    cout << "not possible" << endl;
    abort();
  }


  //! adding transmission terms (for any dimension)
  template<>
  void VarTransmission_Elas<Complex_wp, Dimension3>::
  AddTransmissionTerms(const Complex_wp& alpha, const GlobalGenericMatrix<Complex_wp>& nat_mat,
                       VirtualMatrix<Complex_wp>& mat_sp, int offset_row, int offset_col)
  {
    this->AddTransmissionTermGen(rho_layer, lambda_layer, mu_layer, sigma_layer,
				 alpha, nat_mat, mat_sp, offset_row, offset_col);
  }
#endif


  /***********************************
   * ImpedanceTransmission_Elas_Base *
   ***********************************/

  
  template<class Dimension>
  ImpedanceTransmission_Elas_Base<Dimension>
  ::ImpedanceTransmission_Elas_Base(VarProblem<Dimension>& var,
				    VarTransmission_Elas_Base<Dimension>& var_t)
    : var_problem(var), var_transmission(var_t)
  { 
    curvilinear = false;
    cond_on_gamma_plus = false;
  }
  

  template<>
  void ImpedanceTransmission_Elas_Base<Dimension2>
  ::SetFrenetBasis(const R2& normale, int k)
  {
    VecTau(k).Init(normale(1), -normale(0));
    VecNormale(k) = normale;
  }


  template<>
  void ImpedanceTransmission_Elas_Base<Dimension3>
  ::SetFrenetBasis(const R3& normale, int k)
  {
    GetVectorPlane(normale, VecTauX(k), VecTauY(k));
    VecNormale(k) = normale;
  }
  

  //! evaluation of impedance for tangential derivatives
  template<class Dimension>
  void ImpedanceTransmission_Elas_Base<Dimension>
  ::ConstructFrenetBasis(int i, int num_elem, int num_edge, int num_loc, int k,
			 const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat_)
  {
    SetMatrices<Dimension>& Mat = const_cast<SetMatrices<Dimension>& >(Mat_);
    
    if (k == 0)
      {
	typename Dimension::MatrixN_N mat_dfj, dfjm1;
	typename Dimension::R_N normale_fj; Real_wp dsj(0);	
	for (int k2 = 0; k2 < Pts.GetNbPointsQuadratureBoundary(); k2++)
	  {
	    mat_dfj = var_transmission.GetDfjMatrixQuadraturePoint(num_elem, num_loc, k2);
	    GetInverse(mat_dfj, dfjm1);
	    var_problem.mesh.GetNormale(dfjm1, normale_fj, dsj, num_elem, num_loc);
	    Mat.SetPointQuadratureBoundary(k2, mat_dfj);
	    Mat.SetNormaleQuadratureBoundary(k2, normale_fj);
	    Mat.SetDsQuadratureBoundary(k2, dsj);
	  }
	
	VecNormale.Reallocate(Pts.GetNbPointsQuadratureBoundary());
        if (Dimension::dim_N == 3)
          {
            VecTauX.Reallocate(Pts.GetNbPointsQuadratureBoundary());
            VecTauY.Reallocate(Pts.GetNbPointsQuadratureBoundary());
          }
        else
          VecTau.Reallocate(Pts.GetNbPointsQuadratureBoundary());
      }
   
    SetFrenetBasis(Mat.GetNormaleQuadratureBoundary(k), k);
  }
  

  //! evaluation of f(phi, \nabla phi) in 2-D
  /*!
    This term is involved in the integral \int_\Gamma f(\varphi_i, \nabla \varphi_i) \varphi_j ds
    computed in function AddMatrixImpedanceBoundary of class VarBoundaryCondition
    n is the component number of varphi_i
   */
  template<> template<class T>
  void ImpedanceTransmission_Elas_Base<Dimension2>
  ::ApplyImpedancePhiGen(int n, int j, const TinyVector<Real_wp, 1>& phi_boundary,
			 const R2& dphi, Vector<T>& feval_phi)
  {
    ImpedanceTransmission_Elas<T, Dimension2>& imped
      = static_cast<ImpedanceTransmission_Elas<T, Dimension2>& >(*this);
    
    if (curvilinear)
      {
        TinyVector<T, 2> vec_u, vec_v, vec_w;
	R2 phi_nt, dphi_tn_ds;
	
        // on calcule phi.n et phi.t
        phi_nt(0) = phi_boundary(0)*VecTau(j)(n);
        phi_nt(1) = phi_boundary(0)*VecNormale(j)(n);
        
        // on multiplie par coef_u
        MltTrans(imped.coef_u, phi_nt, vec_u);
	
        // on calcule dphi_i/ds
        Real_wp dphi_ds = DotProd(dphi, VecTau(j));
	
        // il faut en deduire dphi_n/ds et dphi_tau/ds
        dphi_tn_ds(0) = dphi_ds*VecTau(j)(n);
        dphi_tn_ds(1) = dphi_ds*VecNormale(j)(n);
	
        // on multiplie par coef_u_dphi
        MltTrans(imped.coef_u_dphi, dphi_tn_ds, vec_v);
	
        vec_w = vec_u + vec_v;
		
        // maintenant on repasse w en (x, y) (produit avec transpose(P))
        feval_phi(0) = VecTau(j)(0)*vec_w(0) + VecNormale(j)(0)*vec_w(1);
        feval_phi(1) = VecTau(j)(1)*vec_w(0) + VecNormale(j)(1)*vec_w(1);	    
		
        /*TinyMatrix<Complexe, General, Dimension::dim_N, Dimension::dim_N> P, Atmp, Ptrans, coef_ub, coef_u_dphi_t;
          P(0, 0) = VecTau(j)(0); P(0, 1) = VecTau(j)(1);
          P(1, 0) = VecNormale(j)(0); P(1, 1) = VecNormale(j)(1);
          Transpose(P, Ptrans);
          Transpose(coef_u_dphi, coef_u_dphi_t);
	  
          Mlt(coef_u_dphi_t, P, Atmp);
          Mlt(Ptrans, Atmp, coef_ub);
          DISP(P); DISP(coef_u_dphi); DISP(coef_ub);
          
          Complexe feval_x = coef_ub(0, n)*dphi_ds;
          Complexe feval_y = coef_ub(1, n)*dphi_ds;
          DISP(feval_phi); DISP(feval_x); DISP(feval_y);
          if (abs(feval_phi(0) - feval_x) > 1e-6)
          {
          cout << "Value different" << endl;
          abort();
          }
	  
          if (abs(feval_phi(1) - feval_y) > 1e-6)
          {
          cout << "Value different" << endl;
          abort();
          }*/
        
        return;
      }
    
    feval_phi(0) = imped.coef_u(n, 0)*phi_boundary(0) + imped.coef_u_dphi(n, 0)*dphi(0);
    feval_phi(1) = imped.coef_u(n, 1)*phi_boundary(0) + imped.coef_u_dphi(n, 1)*dphi(0);
  }
  

#ifdef MONTJOIE_WITH_THREE_DIM
  //! evaluation of f(phi, \nabla phi) in 3-D
  /*!
    This term is involved in the integral \int_\Gamma f(\varphi_i, \nabla \varphi_i) \varphi_j ds
    computed in function AddMatrixImpedanceBoundary of class VarBoundaryCondition
    n is the component number of varphi_i
   */
  template<> template<class T>
  void ImpedanceTransmission_Elas_Base<Dimension3>
  ::ApplyImpedancePhiGen(int n, int j, const TinyVector<Real_wp, 1>& phi_boundary,
			 const R3& dphi, Vector<T>& feval_phi)
  {
    ImpedanceTransmission_Elas<T, Dimension3>& imped
      = static_cast<ImpedanceTransmission_Elas<T, Dimension3>& >(*this);

    if (curvilinear)
      {
        TinyVector<T, 3> vec_u, vec_vx, vec_vy, vec_w;
	R3 phi_nt, dphi_tn_dtauX, dphi_tn_dtauY;
        
        // on calcule phi dans la base (tau_x, tau_y, n)
        phi_nt(0) = phi_boundary(0)*VecTauX(j)(n);
        phi_nt(1) = phi_boundary(0)*VecTauY(j)(n);
        phi_nt(2) = phi_boundary(0)*VecNormale(j)(n);
	
        // on multiplie par coef_u
        MltTrans(imped.coef_u, phi_nt, vec_u);
        
        // on calcule dphi_dtauX, dphi_dtauY
        Real_wp dphi_dtauX = DotProd(dphi, VecTauX(j));
        Real_wp dphi_dtauY = DotProd(dphi, VecTauY(j));
        	
        // on en deduit dphi_tauX/dtauX, etc
        dphi_tn_dtauX(0) = dphi_dtauX*VecTauX(j)(n);
        dphi_tn_dtauX(1) = dphi_dtauX*VecTauY(j)(n);
        dphi_tn_dtauX(2) = dphi_dtauX*VecNormale(j)(n);

        dphi_tn_dtauY(0) = dphi_dtauY*VecTauX(j)(n);
        dphi_tn_dtauY(1) = dphi_dtauY*VecTauY(j)(n);
        dphi_tn_dtauY(2) = dphi_dtauY*VecNormale(j)(n);
	
        // on multiplie par coef_u_dphi
        MltTrans(imped.coef_u_d1_phi, dphi_tn_dtauX, vec_vx);
        MltTrans(imped.coef_u_d2_phi, dphi_tn_dtauY, vec_vy);
        
        vec_w = vec_u + vec_vx + vec_vy;
        
        // maintenant on repasse w en (x, y, z) (produit avec transpose(P))
        feval_phi(0) = VecTauX(j)(0)*vec_w(0) + VecTauY(j)(0)*vec_w(1) + VecNormale(j)(0)*vec_w(2);
        feval_phi(1) = VecTauX(j)(1)*vec_w(0) + VecTauY(j)(1)*vec_w(1) + VecNormale(j)(1)*vec_w(2);	    
        feval_phi(2) = VecTauX(j)(2)*vec_w(0) + VecTauY(j)(2)*vec_w(1) + VecNormale(j)(2)*vec_w(2);
        
        return;
      }
    
    feval_phi(0) = imped.coef_u(n, 0)*phi_boundary(0) 
      + imped.coef_u_d1_phi(n, 0)*dphi(0) + imped.coef_u_d2_phi(n, 0)*dphi(1);
    
    feval_phi(1) = imped.coef_u(n, 1)*phi_boundary(0) 
      + imped.coef_u_d1_phi(n, 1)*dphi(0) + imped.coef_u_d2_phi(n, 1)*dphi(1);
    
    feval_phi(2) = imped.coef_u(n, 2)*phi_boundary(0) 
      + imped.coef_u_d1_phi(n, 2)*dphi(0) + imped.coef_u_d2_phi(n, 2)*dphi(1);
  }
#endif
  

  //! evaluation of f(phi, \nabla phi) in 2-D
  /*!
    This term is involved in the integral 
    \int_\Gamma f(\varphi_i, \nabla \varphi_i) \nabla \varphi_j ds
    computed in function AddMatrixImpedanceBoundary of class VarBoundaryCondition
   */  
  template<> template<class T>
  void ImpedanceTransmission_Elas_Base<Dimension2>
  ::ApplyImpedanceGradGen(int n, int j, const TinyVector<Real_wp, 1>& phi,
			  const R2& dphi, Vector<T>& geval_phi)
  {    
    ImpedanceTransmission_Elas<T, Dimension2>& imped
      = static_cast<ImpedanceTransmission_Elas<T, Dimension2>& >(*this);

    if (curvilinear)
      {
        TinyVector<T, 2> vec_u, vec_v, vec_w;
        TinyVector<Real_wp, 2> phi_nt, dphi_tn_ds; 
	
        // on calcule phi.n et phi.t
        phi_nt(0) = phi(0)*VecTau(j)(n);
        phi_nt(1) = phi(0)*VecNormale(j)(n);
	
        // on multiplie par coef_phi_du
        Mlt(imped.coef_phi_du, phi_nt, vec_u);
	
        // on calcule dphi_i/ds, l'autre vaut 0
        Real_wp dphi_ds = DotProd(dphi, VecTau(j));
	
        // il faut en deduire dphi_n/ds et dphi_tau/ds
        dphi_tn_ds(0) = dphi_ds*VecTau(j)(n);
        dphi_tn_ds(1) = dphi_ds*VecNormale(j)(n);
	
        // on multiplie par coef_du_dphi
        MltTrans(imped.coef_du_dphi, dphi_tn_ds, vec_v);
	
        vec_w = vec_u + vec_v;
	
        // maintenant on fait w \cdot du/ds
        // mais qu'il faut remettre dans la base cartesienne
        // w_tau d(u \cdot tau)/ds + w_n d(u \cdot n)/ds
        // soit w_tau ( du/ds \cdot tau + n/R \cdot u)  et w_n ( du/ds \cdot n - tau/R \cdot u)
        // ou du/ds = grad u \cdot tau
	
        // d'abord le gradient de u_x
        geval_phi(0) = vec_w(0)*VecTau(j)(0)*VecTau(j)(0) + vec_w(1)*VecNormale(j)(0)*VecTau(j)(0);
        geval_phi(1) = vec_w(0)*VecTau(j)(0)*VecTau(j)(1) + vec_w(1)*VecNormale(j)(0)*VecTau(j)(1);
	
        // ensuite le gradient de u_y
        geval_phi(2) = vec_w(0)*VecTau(j)(1)*VecTau(j)(0) + vec_w(1)*VecNormale(j)(1)*VecTau(j)(0);
        geval_phi(3) = vec_w(0)*VecTau(j)(1)*VecTau(j)(1) + vec_w(1)*VecNormale(j)(1)*VecTau(j)(1);
	
        /* TinyMatrix<Complexe, General, 2, 2> P, Atmp, Ptrans, coef_ub;
           P(0, 0) = VecTau(j)(0); P(0, 1) = VecTau(j)(1);
           P(1, 0) = VecNormale(j)(0); P(1, 1) = VecNormale(j)(1);
           Transpose(P, Ptrans);
	   
           Mlt(coef_phi_du, P, Atmp);
           Mlt(Ptrans, Atmp, coef_ub);
           DISP(P); DISP(coef_phi_du); DISP(coef_ub);
	   
           Complexe feval_xx = coef_ub(0, n)*phi(0)*VecTau(j)(0);
           Complexe feval_xy = coef_ub(0, n)*phi(0)*VecTau(j)(1);
           Complexe feval_yx = coef_ub(1, n)*phi(0)*VecTau(j)(0);
           Complexe feval_yy = coef_ub(1, n)*phi(0)*VecTau(j)(1);
	   
           DISP(j); DISP(VecTau(j)); DISP(VecNormale(j));
           DISP(feval_phi);
           
           if (abs(feval_phi(0) - feval_xx) > 1e-6)
           {
           cout << "Value grad different" << endl;
           abort();
           }
           
           if (abs(feval_phi(1) - feval_xy) > 1e-6)
           {
           cout << "Value grad different" << endl;
           abort();
           }
           
           if (abs(feval_phi(2) - feval_yx) > 1e-6)
           {
           cout << "Value grad different" << endl;
           abort();
           }
           
           if (abs(feval_phi(3) - feval_yy) > 1e-6)
           {
           cout << "Value grad different" << endl;
           abort();
           }
        */
	
        return;
      }
    
    geval_phi(1) = 0.0; geval_phi(3) = 0.0;
    geval_phi(0) = imped.coef_phi_du(0, n)*phi(0) + imped.coef_du_dphi(n, 0)*dphi(0);
    geval_phi(2) = imped.coef_phi_du(1, n)*phi(0) + imped.coef_du_dphi(n, 1)*dphi(0);
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! evaluation of f(phi, \nabla phi) in 3-D
  /*!
    This term is involved in the integral 
    \int_\Gamma f(\varphi_i, \nabla \varphi_i) \nabla \varphi_j ds
    computed in function AddMatrixImpedanceBoundary of class VarBoundaryCondition
   */  
  template<> template<class T>
  void ImpedanceTransmission_Elas_Base<Dimension3>
  ::ApplyImpedanceGradGen(int n, int j, const TinyVector<Real_wp, 1>& phi,
			  const R3& dphi, Vector<T>& feval_phi)
  {    
    ImpedanceTransmission_Elas<T, Dimension3>& imped
      = static_cast<ImpedanceTransmission_Elas<T, Dimension3>& >(*this);

    if (curvilinear)
      {
        TinyVector<T, 3> vec_ux, vec_uy, vec_vx, vec_vy, vec_v12, vec_v21, vec_wx, vec_wy;
        TinyVector<Real_wp, 3> phi_nt, dphi_tn_dtauX, dphi_tn_dtauY; 
        
        // on calcule phi.n et phi.t
        phi_nt(0) = phi(0)*VecTauX(j)(n);
        phi_nt(1) = phi(0)*VecTauY(j)(n);
        phi_nt(2) = phi(0)*VecNormale(j)(n);
        
        Mlt(imped.coef_phi_d1_u, phi_nt, vec_ux);
        Mlt(imped.coef_phi_d2_u, phi_nt, vec_uy);
        
        // on calcule dphi_dtauX, dphi_dtauY
        Real_wp dphi_dtauX = DotProd(dphi, VecTauX(j));
        Real_wp dphi_dtauY = DotProd(dphi, VecTauY(j));
        
        // on en deduit dphi_tauX/dtauX, etc
        dphi_tn_dtauX(0) = dphi_dtauX*VecTauX(j)(n);
        dphi_tn_dtauX(1) = dphi_dtauX*VecTauY(j)(n);
        dphi_tn_dtauX(2) = dphi_dtauX*VecNormale(j)(n);

        dphi_tn_dtauY(0) = dphi_dtauY*VecTauX(j)(n);
        dphi_tn_dtauY(1) = dphi_dtauY*VecTauY(j)(n);
        dphi_tn_dtauY(2) = dphi_dtauY*VecNormale(j)(n);
        
        // on multiplie par coef_du_dphi
        MltTrans(imped.coef_d1u_d1phi, dphi_tn_dtauX, vec_vx);
        MltTrans(imped.coef_d2u_d2phi, dphi_tn_dtauY, vec_vy);
        
        MltTrans(imped.coef_d1u_d2phi, dphi_tn_dtauY, vec_v12);
        Mlt(imped.coef_d1phi_d2u, dphi_tn_dtauX, vec_v21);
        
        vec_wx = vec_ux + vec_vx + vec_v12;
        vec_wy = vec_uy + vec_vy + vec_v21;
        
        // gradient de u_x
        vec_ux(0) = vec_wx(0)*VecTauX(j)(0) + vec_wx(1)*VecTauY(j)(0) + vec_wx(2)*VecNormale(j)(0);
        vec_uy(0) = vec_wy(0)*VecTauX(j)(0) + vec_wy(1)*VecTauY(j)(0) + vec_wy(2)*VecNormale(j)(0);
        feval_phi(0) = VecTauX(j)(0)*vec_ux(0) + VecTauY(j)(0)*vec_uy(0);
        feval_phi(1) = VecTauX(j)(1)*vec_ux(0) + VecTauY(j)(1)*vec_uy(0);
        feval_phi(2) = VecTauX(j)(2)*vec_ux(0) + VecTauY(j)(2)*vec_uy(0);
        
        // gradient de u_y
        vec_ux(1) = vec_wx(0)*VecTauX(j)(1) + vec_wx(1)*VecTauY(j)(1) + vec_wx(2)*VecNormale(j)(1);
        vec_uy(1) = vec_wy(0)*VecTauX(j)(1) + vec_wy(1)*VecTauY(j)(1) + vec_wy(2)*VecNormale(j)(1);
        feval_phi(3) = VecTauX(j)(0)*vec_ux(1) + VecTauY(j)(0)*vec_uy(1);
        feval_phi(4) = VecTauX(j)(1)*vec_ux(1) + VecTauY(j)(1)*vec_uy(1);
        feval_phi(5) = VecTauX(j)(2)*vec_ux(1) + VecTauY(j)(2)*vec_uy(1);

        // gradient de u_z
        vec_ux(2) = vec_wx(0)*VecTauX(j)(2) + vec_wx(1)*VecTauY(j)(2) + vec_wx(2)*VecNormale(j)(2);
        vec_uy(2) = vec_wy(0)*VecTauX(j)(2) + vec_wy(1)*VecTauY(j)(2) + vec_wy(2)*VecNormale(j)(2);
        feval_phi(6) = VecTauX(j)(0)*vec_ux(2) + VecTauY(j)(0)*vec_uy(2);
        feval_phi(7) = VecTauX(j)(1)*vec_ux(2) + VecTauY(j)(1)*vec_uy(2);
        feval_phi(8) = VecTauX(j)(2)*vec_ux(2) + VecTauY(j)(2)*vec_uy(2);
        
        return;
      }
    
    feval_phi(2) = 0.0; feval_phi(5) = 0.0; feval_phi(8) = 0.0;
    feval_phi(0) = imped.coef_phi_d1_u(0, n)*phi(0)
      + imped.coef_d1u_d1phi(n, 0)*dphi(0) + imped.coef_d1u_d2phi(n, 0)*dphi(1);
    
    feval_phi(1) = imped.coef_phi_d2_u(0, n)*phi(0)
      + imped.coef_d2u_d2phi(n, 0)*dphi(1) + imped.coef_d1phi_d2u(0, n)*dphi(0);
    
    feval_phi(3) = imped.coef_phi_d1_u(1, n)*phi(0)
      + imped.coef_d1u_d1phi(n, 1)*dphi(0) + imped.coef_d1u_d2phi(n, 1)*dphi(1);
    
    feval_phi(4) = imped.coef_phi_d2_u(1, n)*phi(0)
      + imped.coef_d2u_d2phi(n, 1)*dphi(1) + imped.coef_d1phi_d2u(1, n)*dphi(0);
    
    feval_phi(6) = imped.coef_phi_d1_u(2, n)*phi(0)
      + imped.coef_d1u_d1phi(n, 2)*dphi(0) + imped.coef_d1u_d2phi(n, 2)*dphi(1);
    
    feval_phi(7) = imped.coef_phi_d2_u(2, n)*phi(0)
      + imped.coef_d2u_d2phi(n, 2)*dphi(1) + imped.coef_d1phi_d2u(2, n)*dphi(0);
    
  }
#endif
  
    
  /******************************
   * ImpedanceTransmission_Elas *
   ******************************/

  
  template<class T, class Dimension>
  ImpedanceTransmission_Elas<T, Dimension>
  ::ImpedanceTransmission_Elas(VarProblem<Dimension>& var,
			       VarTransmission_Elas_Base<Dimension>& var_t)
    : ImpedanceTransmission_Elas_Base<Dimension>(var, var_t),
      ImpedanceFunction_Base<T, Dimension>(var)
  {
  }

  
  //! evaluation of impedance for tangential derivatives
  template<class T, class Dimension>
  void ImpedanceTransmission_Elas<T, Dimension>
  ::EvaluateImpedanceGrad(int i, int num_elem, int num_edge, int num_loc, int k,
			  const GlobalGenericMatrix<T>& nat_mat, int ref,
			  const SetPoints<Dimension>& Pts,
			  const SetMatrices<Dimension>& Mat)
  {
    this->ConstructFrenetBasis(i, num_elem, num_edge, num_loc, k, Pts, Mat);
  }
  
  
  //! evaluation of f(phi, \nabla phi)
  /*!
    This term is involved in the integral \int_\Gamma f(\varphi_i, \nabla \varphi_i) \varphi_j ds
    computed in function AddMatrixImpedanceBoundary of class VarBoundaryCondition
    n is the component number of varphi_i
   */
  template<class T, class Dimension>
  void ImpedanceTransmission_Elas<T, Dimension>
  ::ApplyImpedancePhi_H1(int n, int j, int, const TinyVector<Real_wp, 1>& phi_boundary,
		      const typename Dimension::R_N& dphi, Vector<T>& feval_phi)
  {
    this->ApplyImpedancePhiGen(n, j, phi_boundary, dphi, feval_phi);
  }
  
  
  //! evaluation of f(phi, \nabla phi)
  /*!
    This term is involved in the integral 
    \int_\Gamma f(\varphi_i, \nabla \varphi_i) \nabla \varphi_j ds
    computed in function AddMatrixImpedanceBoundary of class VarBoundaryCondition
   */  
  template<class T, class Dimension>
  void ImpedanceTransmission_Elas<T, Dimension>
  ::ApplyImpedanceGrad(int n, int j, int, const TinyVector<Real_wp, 1>& phi,
		       const typename Dimension::R_N& dphi, Vector<T>& feval_phi)
  {
    this->ApplyImpedanceGradGen(n, j, phi, dphi, feval_phi);
  }
  
}// end namespace

#define MONTJOIE_FILE_TRANSMISSION_MODEL_ELASTIC_CXX
#endif

  
