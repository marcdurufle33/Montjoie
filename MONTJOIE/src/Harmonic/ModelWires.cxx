#ifndef MONTJOIE_FILE_MODEL_WIRES_CXX
 
namespace Montjoie
{
#ifdef MONTJOIE_WITH_TWO_DIM 
  //! adding dofs for wires
  template<class TypeElement, class TypeEquation>
  void VarHarmonic<TypeElement,TypeEquation>::IncrementNumberDofs_Wires(Dimension2 & dim)
  {
    //On ajoute des degres de liberte 
    // cout<<"ATTENTION AU nb_fils"<<endl;    
    int nb_fils = Tab_Fil.GetM();
    if (nb_fils <= 0)
      return;
        
    num_ddl_fil.Reallocate(nb_fils);
    for (int i = 0; i < nb_fils; i++)
      {
	num_ddl_fil(i) = this->nodl++;
	// cout<<"Numero des ddl ajoutes"<<endl;
	// DISP(num_ddl_fil(i));
      }
    
    if (Tab_Fil(0).Fct_Add != GlobalVariables::FCT_WIRE_AUCUNE)
      {
	num_ddl_fctfil.Reallocate(nb_fils);
	
	for (int i = 0; i < nb_fils; i++)
	  num_ddl_fctfil(i) = this->nodl++;
      }
  }
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
  //! incrementation of dofs in 3-D case
  template<class TypeElement, class TypeEquation>
  void VarHarmonic<TypeElement,TypeEquation>::IncrementNumberDofs_Wires(Dimension3 & dim)
  {
    int nb_fils = Tab_Fil.GetM();
    if (this->lagrange_wire)
      for (int i = 0; i < nb_fils; i++)
	Tab_Fil(i).IncrementDofs(this->nodl);
    
    
    for (int i = 0; i < nb_fils; i++)
      if (Tab_Fil(i).Fct_Add != GlobalVariables::FCT_WIRE_AUCUNE)
	Tab_Fil(i).Increment_AdditionalFunctions(this->nodl);
    
    // DISP(this->nodl); exit(0);
  }
  
  //! modification of matrix to consider wires
  template<class TypeElement, class TypeEquation> template<class MatrixSparse>
  void VarHarmonic<TypeElement,TypeEquation>::TreatWires(MatrixSparse& mat_sp, Dimension3& dim)
  {
    int nb_fils = this->Tab_Fil.GetM();
    if (nb_fils <= 0)
      return;
    
    if (holland_wire)
      return TreatWire_Holland(mat_sp);
    
    // computation of different operators 
    Vector<Matrix<Real_wp, General, ArrayRowSparse>>
      Cfict, Badd_stiff, Badd_mass, beta_stiff, beta_mass;
    this->var_deriv->ComputeOperator_FictitiousDomain(Cfict, Badd_stiff, Badd_mass,
                                                      beta_stiff, beta_mass);
    
    Real_wp value;
    for (int num_fil = 0; num_fil < nb_fils; num_fil++)
      {
	Wire<Dimension3>& wire = this->Tab_Fil(num_fil);
	// loop on each degree of freedom of the wire
	for (int i1 = 0; i1 < wire.GetNbDof(); i1++)
	  {
	    // adding C and C^t
	    int i = this->Tab_Fil(num_fil).GetNumberDof(i1);
	    if (this->lagrange_wire)
	      for (int k = 0; k < Cfict(num_fil).GetRowSize(i1); k++)
		{
		  int j = Cfict(num_fil).Index(i1,k);
		  value = Cfict(num_fil).Value(i1,k);
		  mat_sp.AddInteraction(i, j, value);
		  mat_sp.AddInteraction(j, i, value);
		}
	    // WriteMatrix_Matlab("Cfict.dat", Cfict(num_fil));
	    
	    if (Tab_Fil(num_fil).Fct_Add != GlobalVariables::FCT_WIRE_AUCUNE)
	      {
		// adding B and B^t
		i = this->Tab_Fil(num_fil).GetNumberDof_FctAdd(i1);
		for (int k = 0; k < Badd_mass(num_fil).GetRowSize(i1); k++)
		  {
		    int j = Badd_mass(num_fil).Index(i1,k);
		    value = -this->omega2*Badd_mass(num_fil).Value(i1,k);
		    mat_sp.AddInteraction(i, j, value);
		    mat_sp.AddInteraction(j, i, value);
		  }
		
		for (int k = 0; k < Badd_stiff(num_fil).GetRowSize(i1); k++)
		  {
		    int j = Badd_stiff(num_fil).Index(i1,k);
		    value = Badd_stiff(num_fil).Value(i1,k);
		    mat_sp.AddInteraction(i, j, value);
		    mat_sp.AddInteraction(j, i, value);
		  }
		// WriteMatrix_Matlab("Bmass.dat",Badd_mass(num_fil));
		// WriteMatrix_Matlab("Bstiff.dat",Badd_stiff(num_fil));
		
		// adding beta
		for (int k = 0; k < beta_mass(num_fil).GetRowSize(i1); k++)
		  {
		    int j = this->Tab_Fil(num_fil)
                      .GetNumberDof_FctAdd(beta_mass(num_fil).Index(i1,k));
		    value = -this->omega2*beta_mass(num_fil).Value(i1,k);
		    mat_sp.AddInteraction(i, j, value);
		    mat_sp.AddInteraction(j, i, value);
		  }
		
		for (int k = 0; k < beta_stiff(num_fil).GetRowSize(i1); k++)
		  {
		    int j = this->Tab_Fil(num_fil)
                      .GetNumberDof_FctAdd(beta_stiff(num_fil).Index(i1,k));
		    value = beta_stiff(num_fil).Value(i1,k);
		    mat_sp.AddInteraction(i, j, value);
		    mat_sp.AddInteraction(j, i, value);
		  }
		// WriteMatrix_Matlab("beta_mass.dat",beta_mass(num_fil));
		// WriteMatrix_Matlab("beta_stiff.dat",beta_stiff(num_fil));
		
	      }
	  }
      }
  }
  
  //! computation of matrices coming from wire model
  template<class TypeElement, class TypeEquation> template<class VectMatrix>
  void VarHarmonic<TypeElement,TypeEquation>
  ::ComputeOperator_FictitiousDomain(VectMatrix& Cfict, VectMatrix& Badd_stiff, 
                                     VectMatrix& Badd_mass, VectMatrix& beta_stiff,
                                     VectMatrix& beta_mass)
  {
    int nb_fils = this->Tab_Fil.GetM();
    // allocation of output matrices
    Cfict.Reallocate(nb_fils); Badd_stiff.Reallocate(nb_fils); Badd_mass.Reallocate(nb_fils);
    beta_mass.Reallocate(nb_fils); beta_stiff.Reallocate(nb_fils);
    for (int i = 0; i < nb_fils; i++)
      {
	int nb_dof = this->Tab_Fil(i).GetNbDof();
	Cfict(i).Reallocate(nb_dof, this->mesh.nodl);
	Badd_stiff(i).Reallocate(nb_dof, this->mesh.nodl);
	Badd_mass(i).Reallocate(nb_dof, this->mesh.nodl);
	beta_stiff(i).Reallocate(nb_dof, nb_dof);
	beta_mass(i).Reallocate(nb_dof, nb_dof);
      }
    
    
    ////////////////////////
    // INTEGRALS ON GAMMA //
    
    // For each wire, we generate quadrature points
    VectR3 Mesh_Glob_Wire, Direction_Glob_Wire; VectReal_wp Weights_Glob_Wire;
    Vector<Matrix<Real_wp, General, ArrayRowSparse>> Psi_Glob_Wire(nb_fils);
    VectR3 PointsQuad, DirectionQuad; VectReal_wp WeightsQuad;
    
    for(int i = 0; i < nb_fils; i++)
      {
	this->Tab_Fil(i).Gen_Mesh_Wire(PointsQuad, WeightsQuad, DirectionQuad, Psi_Glob_Wire(i));
	Append(Mesh_Glob_Wire, PointsQuad); Append(Weights_Glob_Wire, WeightsQuad);
	Append(Direction_Glob_Wire, DirectionQuad);
      }
    
    // DISP(Mesh_Glob_Wire); DISP(Direction_Glob_Wire); DISP(Weights_Glob_Wire);
    // DISP(Psi_Glob_Wire);
    // checking that the surface of the wire is well computed
    Real_wp poids(0);
    for (int i = 0; i < Weights_Glob_Wire.GetM(); i++)
      poids += Weights_Glob_Wire(i);
    
    if ((rank_proc == 0)||(this->print_level >= 10))
      if (this->print_level >= 2)
	cout<<"Le fil a une surface de "<<poids<<endl;
    
    int taille = Mesh_Glob_Wire.GetM(); 
    IVect num(taille);

    // and now we localize the quadrature points on the mesh
    GridInterpolation<Dimension>& loc = this->all_points_display;
    loc.Append(Mesh_Glob_Wire, num);
    loc.LocalizePoints(this->mesh);
    
    ArrayVectReal_wp val; Real_wp value, phi_scal_nu, troncat;
    // loop on wires
    int old_taille = 0; taille = 0;
    for (int i = 0; i < nb_fils; i++)
      {
	old_taille = taille;
	Wire<Dimension3>& wire = this->Tab_Fil(i);
	taille += wire.GetNbPointsQuadrature();
	// DISP(old_taille); DISP(taille);
	
	// boucle sur les points d'integration concernes
	for (int j = old_taille; j < taille; j++)
	  { 
	    // calcul de val: phi_i(x_j) pour tout i, et x_j = loc.CoorInterp(num(j))
	    int num_elem = loc.ElementInterp(num(j));
	    this->ComputeValuesPhi(loc.CoorInterp(num(j)),
                                   val, loc.dfjm1(num(j)), num_elem);
	    
	    troncat = wire.GetFunction_Troncature1D(wire.GetEps(), wire.Fct_Add);
	    troncat /= wire.GetEps();
	    if (num_elem < 0)
	      {
		cout<<"le fil n'est pas dans le maillage, changez la position du fil"<<endl;
		exit(0);
	      }
	    
	    // loop on basis functions of the 3-D mesh
	    for (int k = 0; k < val.GetM(); k++)
	      {
		// global number of the 3-D basis function
		int num_glob1 = this->mesh.GetNumberDof(num_elem, k);
		phi_scal_nu = this->var_deriv->
                  DotProd_WireDirection(val(k), Direction_Glob_Wire(j));
		
		// DISP(k); DISP(phi_scal_nu); DISP(Direction_Glob_Wire(j));
		// loop on 1-D basis functions of the wire
		for (int m = 0; m < Psi_Glob_Wire(i).GetRowSize(j-old_taille); m++)
		  {
		    int num_glob_fil = Psi_Glob_Wire(i).Index(j-old_taille,m);
		    
		    value = Psi_Glob_Wire(i).Value(j-old_taille,m)
                      *phi_scal_nu*Weights_Glob_Wire(j);
		    
                    // DISP(num_glob1); DISP(num_glob_fil); DISP(value);
		    
		    
		    if (num_glob_fil >= 0)
		      Cfict(i).AddInteraction(num_glob_fil, num_glob1, value);
		    
		    if (wire.Fct_Add != GlobalVariables::FCT_WIRE_AUCUNE)
		      {
			// part -\int_\Gamma 1 / epsilon  f \Xi \psi 
			value *= -troncat;
			if (num_glob_fil >= 0)
			  Badd_stiff(i).AddInteraction(num_glob_fil, num_glob1, value);
		      }
		  }
	      }
	  }
      }
    
    // INTEGRALS ON GAMMA //
    ////////////////////////
    
    ///////////////////////
    // VOLUMIC INTEGRALS //
    
    // in this case, the wire is always straight
    Mesh<Dimension3> & mesh = this->mesh;
    TypeElement* Fb = this->FaceBasis;
    
    // same integration for all the elements
    Globatto<Dimension3> gauss;
    gauss.ConstructQuadrature(this->Tab_Fil(0).Order_Quad_Wire);
    
    // basis functions on reference points
    MatrixFullR3 Value_Phi(Fb->nb_dof_loc, gauss.nb_points_quadrature),
      Gradient_Phi(Fb->nb_dof_loc, gauss.nb_points_quadrature);
    
    for (int i = 0; i < Fb->nb_dof_loc; i++)
      for (int j = 0; j < gauss.nb_points_quadrature; j++)
	{
	  Fb->GetValue_Phi(i, gauss.Points3D(j), Value_Phi(i,j));
	  Fb->GetGradient_Phi(i, gauss.Points3D(j), Gradient_Phi(i,j));
	}
    
    SetPoints<Dimension3> PointsElem(*this->var_deriv);
    SetMatrices<Dimension3> MatricesElem(*this->var_deriv);
    VectR3 s(8); R3 pt_glob, axis, proj_on_axis; Matrix3_3 mat_dfj, dfjm1;
    int num_seg; Real_wp coor_loc_seg, radius; 
    
    VectReal_wp val_phi(3); R3 vec_er, vec_ez, grad_phi;
    Real_wp grad_troncat1D, fct_log;
    R3 grad_fct_log, grad_troncat; Real_wp vol(0);
    int old_percent = 0, new_percent = 0;
    for (int num_fil = 0; num_fil < nb_fils; num_fil++)
      {
	Wire<Dimension3>& wire = this->Tab_Fil(num_fil);
	VectReal_wp val_phi_wire(wire.GetNbDof_Element(0)),
          grad_phi_wire(wire.GetNbDof_Element(0));
	
	if (Tab_Fil(num_fil).Fct_Add != GlobalVariables::FCT_WIRE_AUCUNE)
	  for (int i = 0; i < mesh.nb_elt; i++)
	    {
	      // barre de progression
	      if (this->print_level >= 1)
		{
		  new_percent = toInteger(round(Real_wp(i)/(mesh.nb_elt-1)*80));
		  for (int percent = old_percent; percent < new_percent; percent++)
		    { cout<<"#"; cout.flush(); }
		  
		  old_percent = new_percent;
		}
	      
	      bool non_null = false;
	      // on effectue une localisation, pour savoir qui interagit
	      mesh.GetVerticesElement(i,s);
	      Fb->FjElemNodal(s, PointsElem, mesh, i);
	      
	      Real_wp epsilon = wire.GetEps();
	      
	      // distance of each point of the axis
	      for (int j = 0; j < Fb->nb_points_nodal_elt; j++)
		{
		  // localisation of nodal point on the wire
		  // radius is the distance from the wire
		  // coor_loc_seg and num_seg are local coordinates and edge number of the wire
		  // axis is the tangent vector to the wire at this point
		  // proj_on_axis is the projection of the nodal point on the wire
		  bool test_inside = wire.FjInv(PointsElem.GetPointNodal(j), num_seg, coor_loc_seg,
						radius, axis, proj_on_axis);
		  
		  // DISP(j);DISP(radius); DISP(PointsElem.GetPointNodal(j));
		  if ((radius < wire.Rayon_Fct_Wire)&&(test_inside))
		    non_null = true;
		}
	      
	      if (non_null)
		{
		  for (int j = 0; j < gauss.nb_points_quadrature; j++)
		    {
		      Fb->Fj(s, PointsElem, gauss.Points3D(j), pt_glob, mesh, i);
		      bool test_inside = wire.FjInv(pt_glob, num_seg, coor_loc_seg,
						    radius, axis, proj_on_axis);
		      
		      vec_ez = axis;

		      if ((radius < wire.Rayon_Fct_Wire)&&(test_inside)&&(radius > epsilon))
			{
			  Fb->DFj(s, PointsElem, gauss.Points3D(j), mat_dfj, mesh, i);
			  Inverse(mat_dfj, dfjm1);
			  Real_wp poids = Det(mat_dfj)*gauss.Weights3D(j);
			  
			  // troncature function and gradient (troncat and grad_troncat)
			  troncat = wire.GetFunction_Troncature1D(radius, wire.Fct_Add);
			  grad_troncat1D = wire.GetGradient_Troncature1D(radius);
			  Subtract(pt_glob, proj_on_axis, vec_er);
			  Mlt(1.0/Norm2(vec_er), vec_er);
			  grad_troncat = vec_er;
                          // DISP(troncat); DISP(grad_troncat); DISP(pt_glob);
			  
                          Mlt(grad_troncat1D, grad_troncat);
			  
			  // logarithmic function
			  fct_log = log(radius/epsilon);
			  grad_fct_log = vec_er; Mlt(1.0/radius, grad_fct_log);
			  
			  if (num_seg >= 0)
			    {
			      for (int m = 0; m < wire.GetNbDof_Element(num_seg); m++)
				{
				  val_phi_wire(m) = wire.GetValue_Phi(m, coor_loc_seg);
				  // gradient is oriented by axis
				  grad_phi_wire(m)
                                    = wire.GetGradient_Phi(num_seg, m, coor_loc_seg);
				}
			      // DISP(poids); DISP(val_phi_wire); DISP(vec_ez);
			      vol += poids*square(val_phi_wire(0)*troncat);
			      
			      for (int k = 0; k < Fb->nb_dof_loc; k++)
				{
				  int num_glob_phi = mesh.GetNumberDof(i,k);
				  // scalar product axis with basis functions
				  Fb->ApplyDF_ValuePhi(Value_Phi(k,j), mat_dfj, dfjm1, val_phi,
						       this->var_deriv, this->FaceBasis,
                                                       mesh, i, k);
				  
                                  Fb->ApplyDF_GradientPhi(Gradient_Phi(k,j), mat_dfj,
                                                          dfjm1, grad_phi,
							  this->var_deriv, this->FaceBasis,
                                                          mesh, i, k);
				  
				  for (int m = 0; m < wire.GetNbDof_Element(num_seg); m++)
				    {
				      int num_glob_fil = wire.GetNumberDof(num_seg, m);
				      // part \int grad(\Xi \psi) \cdot grad \varphi log(r/eps)  
				      value = DotProd(grad_troncat, grad_phi)
                                        *val_phi_wire(m) + DotProd(vec_ez, grad_phi)
                                        *grad_phi_wire(m)*troncat;
                                      
                                      // value = DotProd(grad_troncat, grad_phi)*val_phi_wire(m);
				      value *= fct_log; // value = 0.0;
				      // part - \int phi grad( ln(r/eps) ) \cdot \grad( \Xi) \psi
				      value -=  DotProd(grad_fct_log, grad_troncat)
                                        *val_phi(0)*val_phi_wire(m);
				      
                                      value *= poids;
				      if (num_glob_fil >= 0)
					Badd_stiff(num_fil).
                                          AddInteraction(num_glob_fil, num_glob_phi, value);
				      
				      // mass part
				      value = troncat*val_phi(0)*fct_log*val_phi_wire(m)*poids;
				      if (num_glob_fil >= 0)
					Badd_mass(num_fil).
                                          AddInteraction(num_glob_fil, num_glob_phi, value);
				    }
				}
			      
			      // computation of beta
			      // DISP(pt_glob); DISP(radius); DISP(fct_log); DISP(troncat);
			      // DISP(coor_loc_seg); DISP(val_phi_wire);
			      for (int m = 0; m < wire.GetNbDof_Element(num_seg); m++)
				{
				  int num_glob_fil = wire.GetNumberDof(num_seg, m);
				  for (int n = 0; n < wire.GetNbDof_Element(num_seg); n++)
				    {
				      int num_glob2_fil = wire.GetNumberDof(num_seg, n);
				      // stiffness part
				      value = DotProd(grad_troncat, grad_troncat)*val_phi_wire(m)
                                        *val_phi_wire(n) + troncat*troncat*grad_phi_wire(m)
                                        *grad_phi_wire(n);
				      
                                      value *= poids*fct_log*fct_log;
				      if ((num_glob_fil >= 0)&&(num_glob2_fil >= 0))
					beta_stiff(num_fil).
                                          AddInteraction(num_glob_fil, num_glob2_fil, value);
				      
				      // mass part
				      value = val_phi_wire(m)*val_phi_wire(n)
                                        *square(troncat*fct_log)*poids;
				      
                                      if ((num_glob_fil >= 0)&&(num_glob2_fil >= 0))
					beta_mass(num_fil).
                                          AddInteraction(num_glob_fil, num_glob2_fil, value);
				    }
				}
			    }
			}
		    }
		}
	    }
      }
    // DISP(vol); DISP(beta_mass(0)(0,0));
    
    // VOLUMIC INTEGRALS //
    ///////////////////////
  }
  
  //! modification of u to take into account additional functions
  template<class TypeElement, class TypeEquation>
  void VarHarmonic<TypeElement,TypeEquation>
  ::AddComponent_AdditionalFunction(const VectComplexe& U0, int num_fil,
                                    const R3& pt_glob, VectComplexe& val_u) const
  {
    const Wire<Dimension3>& wire = this->Tab_Fil(num_fil);
    int num_seg; Real_wp coor_loc_seg(0), radius; R3 axis, proj_on_axis;
    bool test_inside = wire.FjInv(pt_glob, num_seg, coor_loc_seg,
				  radius, axis, proj_on_axis);
    Real_wp troncat(0);
    if ((radius > wire.GetEps())&&(radius < wire.Rayon_Fct_Wire)&&(test_inside))
      {
	troncat = wire.GetFunction_Troncature1D(radius, wire.Fct_Add);
	troncat *= log(radius/wire.GetEps());
      }
    else
      return;
    
    if (num_seg >= 0)
      for (int i = 0; i < wire.GetNbDof_Element(num_seg); i++)
	{
	  int dof = wire.GetNumberDof(num_seg, i);
	  if (dof >= 0)
	    {
	      int num_glob = wire.GetNumberDof_FctAdd(dof);
	      Real_wp phi = wire.GetValue_Phi(i, coor_loc_seg);
	      val_u(0) += phi*troncat*U0(num_glob);
	    }
	}
  }
  
  //! modification of matrix for Holland model
  template<class TypeElement, class TypeEquation> template<class MatrixSparse>
  void VarHarmonic<TypeElement,TypeEquation>::TreatWire_Holland(MatrixSparse& mat_sp)
  {
    Vector<Matrix<Real_wp, General, ArrayRowSparse>> Glob_mat_Bholland;
    // computing all the operators for Holland model
    this->var_deriv->ComputeHollandOperator(Glob_mat_Bholland);
    
    // now we fill the matrix
    int nb_fils = this->Tab_Fil.GetM();
    for (int num_fil = 0; num_fil < nb_fils; num_fil++)
      {
	Wire<Dimension3>& wire = this->Tab_Fil(num_fil);
	Real_wp Lapp = wire.GetLapp(), value(0);
	for (int i1 = 0; i1 < wire.GetNbDof(); i1++)
	  {
	    int i = wire.GetNumberDof(i1);
	    for (int j1 = 0; j1 < Glob_mat_Bholland(num_fil).GetRowSize(i1); j1++)
	      {
		int j = Glob_mat_Bholland(num_fil).Index(i1,j1);
		value = -Glob_mat_Bholland(num_fil).Value(i1,j1);
		// DISP(i); DISP(j); DISP(value);
		mat_sp.AddInteraction(i, j, value);
		mat_sp.AddInteraction(j, i, value);
	      }
	    
	    // mass matrix of wire
	    value = wire.GetMass(i1)*Lapp;
	    mat_sp.AddInteraction(i, i, value);
	    // DISP(i); DISP(value);
	    
	    // stiffness matrix of wire
	    for (int j1 = 0; j1 < wire.mat_stiffness.GetRowSize(i1); j1++)
	      {
		int j = wire.GetNumberDof(wire.mat_stiffness.Index(i1,j1));
		value = wire.mat_stiffness.Value(i1,j1)*Lapp/this->omega2;
		mat_sp.AddInteraction(i, j, value);
		// DISP(i); DISP(j); DISP(value);
	      }
	  }
      }
  }
  
  //! computation of matrix involved in Holland model
  template<class TypeElement, class TypeEquation> template<class VectMatrix>
  void VarHarmonic<TypeElement,TypeEquation>::ComputeHollandOperator(VectMatrix& Glob_mat_Bholland)
  {
    int nb_fils = this->Tab_Fil.GetM();
    if ((this->holland_wire)&&(nb_fils > 0))
      {
	// in this case, the wire is always straight
	Mesh<Dimension3> & mesh = this->mesh;
	TypeElement* Fb = this->FaceBasis;
	
	// same integration for all the elements
	Globatto<Dimension3> gauss;
	gauss.ConstructQuadrature(this->Tab_Fil(0).Order_Quad_Wire);
	// basis functions on reference points
	MatrixFullR3 Value_Phi(Fb->nb_dof_loc, gauss.nb_points_quadrature);
	for (int i = 0; i < Fb->nb_dof_loc; i++)
	  for (int j = 0; j < gauss.nb_points_quadrature; j++)
	    Fb->GetValue_Phi(i, gauss.Points3D(j), Value_Phi(i,j));
	
	SetPoints<Dimension3> PointsElem(*this->var_deriv);
	SetMatrices<Dimension3> MatricesElem(*this->var_deriv);
	VectR3 s(8); R3 pt_glob, pt_loc, axis; Matrix3_3 mat_dfj, dfjm1; Real_wp value;
	int num_seg; Real_wp coor_loc_seg; R3 proj_on_axis;
	
	int Nvol = mesh.nodl;
	Glob_mat_Bholland.Reallocate(nb_fils); VectReal_wp val_phi(3);
	for (int num_fil = 0; num_fil < nb_fils; num_fil++)
	  {
	    Wire<Dimension3>& wire = this->Tab_Fil(num_fil);
	    VectReal_wp val_phi_wire(wire.GetNbDof_Element(0));
	    // construction of mass and stiffness matrices of the wire
	    wire.ConstructStiffnessMatrix();
	    
	    // computation of \int \delta_app 
	    // in order to normalize (ie, we want to have \int \delta_app = 1
	    Real_wp norme_Lapp(0), eval_Lapp(0), radius;
	    for (int i = 0; i < gauss.Points.GetM(); i++)
	      {
		radius = wire.Rayon_Fct_Wire*gauss.Points(i);
		Real_wp delta_app = wire.GetDeltaApp_Holland(radius);
		Real_wp poids = gauss.Weights(i)*wire.Rayon_Fct_Wire;
		norme_Lapp += delta_app*radius*poids;
		
		radius = wire.GetEps()*(1.0-gauss.Points(i)) + wire.Rayon_Fct_Wire*gauss.Points(i);
		delta_app = wire.GetDeltaApp_Holland(radius);
		poids = gauss.Weights(i)*(wire.Rayon_Fct_Wire-wire.GetEps());
		eval_Lapp += delta_app*radius*log(radius/wire.GetEps())*poids;
	      }
	    eval_Lapp /= norme_Lapp; norme_Lapp *= 2.0*pi_wp;
	    
	    // matrice B stockee en creux, de maniere bourine
	    Glob_mat_Bholland(num_fil).Reallocate(wire.GetNbDof(), Nvol);
	    Real_wp Lepsilon = 0.0, poids_Lapp = 0.0;
	    for (int i = 0; i < mesh.nb_elt; i++)
	      {
		bool non_null = false;
		// on effectue une localisation, pour savoir qui interagit
		mesh.GetVerticesElement(i,s);
		Fb->FjElemNodal(s, PointsElem, mesh, i);
		
		// distance of each point of the axis
		for (int j = 0; j < Fb->nb_points_nodal_elt; j++)
		  {
		    // localisation of nodal point on the wire
		    // radius is the distance from the wire
		    // coor_loc_seg and num_seg are local coordinates and edge number of the wire
		    // axis is the tangent vector to the wire at this point
		    // proj_on_axis is the projection of the nodal point on the wire
		    bool test_inside = wire.FjInv(PointsElem.GetPointNodal(j), num_seg,
                                                  coor_loc_seg,
						  radius, axis, proj_on_axis);
		    
		    // DISP(i); DISP(j); DISP(PointsElem.GetPointNodal(j));
		    // DISP(radius); DISP(wire.Rayon_Fct_Wire); DISP(axis); DISP(proj_on_axis);
		    
		    if ((radius < wire.Rayon_Fct_Wire)&&(test_inside))
		      non_null = true;
		  }		
		
		if (non_null)
		  {
		    for (int j = 0; j < gauss.nb_points_quadrature; j++)
		      {
			Fb->Fj(s, PointsElem, gauss.Points3D(j), pt_glob, mesh, i);
			// localisation of the quadrature point on the wire
			bool test_inside = wire.FjInv(pt_glob, num_seg, coor_loc_seg,
						      radius, axis, proj_on_axis);
			
			// DISP(i); DISP(j); DISP(pt_glob); DISP(num_seg); DISP(coor_loc_seg);
			// DISP(radius); DISP(wire.Rayon_Fct_Wire); DISP(axis); DISP(proj_on_axis);
			// if ((radius < wire.Rayon_Fct_Wire)
                        // &&(Labs > -wire.Rayon_Fct_Wire)&&(Labs < (1+wire.Rayon_Fct_Wire)))
			if ((radius < wire.Rayon_Fct_Wire)&&(test_inside))
			  {
			    Fb->DFj(s, PointsElem, gauss.Points3D(j), mat_dfj, mesh, i);
			    Inverse(mat_dfj, dfjm1);
			    // poids_Lapp = \int_\Omega \delta_app 
			    // Lepsilon = 1/(2 pi) \int_\Omega log(r/a) \delta_app 
			    Real_wp delta_app = wire.GetDeltaApp_Holland(radius);
			    Real_wp poids = Det(mat_dfj)*gauss.Weights3D(j)*delta_app;
			    poids_Lapp += poids;
			    if (radius > wire.GetEps())
			      Lepsilon += poids*log(radius/wire.GetEps());
			    
			    // DISP(num_seg); DISP(coor_loc_seg); // DISP(radius);
			    for (int m = 0; m < wire.GetNbDof_Element(num_seg); m++)
			      val_phi_wire(m) = wire.GetValue_Phi(m, coor_loc_seg);
			    
			    // DISP(val_phi_wire); // DISP(axis);
			    for (int k = 0; k < Fb->nb_dof_loc; k++)
			      {
				int num_glob_phi = mesh.GetNumberDof(i,k);
				// scalar product axis with basis functions
				Fb->ApplyDF_ValuePhi(Value_Phi(k,j), mat_dfj, dfjm1, val_phi,
						     this->var_deriv, this->FaceBasis, mesh, i, k);
				
				Real_wp phi_scal_nu
                                  = this->var_deriv->DotProd_WireDirection(val_phi, axis);
				
				// DISP(k); DISP(phi_scal_nu);
				// loop on basis functions of the wire
				for (int m = 0; m < wire.GetNbDof_Element(num_seg); m++)
				  {
				    int num_glob_fil = wire.GetNumberDof(num_seg, m);
				    value = poids*phi_scal_nu*val_phi_wire(m);
				    if (num_glob_fil >= 0)
				      Glob_mat_Bholland(num_fil).
                                        AddInteraction(num_glob_fil, num_glob_phi, value);
				  }
				
			      } // end if quadrature point near the wire
			    
			  }
			
		      } // end loop on quadrature points of the element
		  
		  }  // end if non null integrals
		
	      } // end loop on elements
	    Real_wp length_wire = wire.GetLength();
	    poids_Lapp /= length_wire;
	    Lepsilon /= length_wire*poids_Lapp;
	    
	    // DISP(Lepsilon); DISP(poids_Lapp); DISP(length_wire);
	    // DISP(norme_Lapp); DISP(eval_Lapp);
	    //Mlt(1.0/poids_Lapp, Glob_mat_Bholland(num_fil));
	    Mlt(1.0/norme_Lapp, Glob_mat_Bholland(num_fil));
	    // WriteMatrix_Matlab("mat_holl.dat",Glob_mat_Bholland(num_fil));
	  } // end loop on wires
      }
  }
  
#endif
    
#ifdef MONTJOIE_WITH_TWO_DIM 
  //! modification of matrix to take into account wires
  template<class TypeElement, class TypeEquation> template<class MatrixSparse>
  void VarHarmonic<TypeElement,TypeEquation>::TreatWires(MatrixSparse& mat_sp, Dimension2& dim)
  {
    
    VectR_N s;
    SetPoints<Dimension2> PointsElem(*this->var_deriv);
    SetPoints<Dimension2> PointsElem2(*this->var_deriv);
    Matrix2_2 dfj,inv_dfj;
    R2 pt_glob;
    Complex_wp Int_e = 0;

    ArrayVectReal_wp phi(this->GetNbMaxDof()),gradphi(this->GetNbMaxDof());
    
    
    // Je suppose que je me suis donne un tableau de fils Wire Tab_Fil[nb_fils];
    int taille = 0; int nb_fils = Tab_Fil.GetM();
    if (nb_fils <= 0)
      return;
    
    for(int i = 0; i < nb_fils; i++)
      {
	taille += Tab_Fil(i).GetPrec_Mesh_Wire();
      }
    
    IVect num(taille);
    
    
    // Pour chaque fil je genere un ensemble de points de quadrature
    VectR_N Mesh_Glob_Wire;
    
    for(int i = 0; i < nb_fils; i++)
      {
	Tab_Fil(i).Gen_Mesh_Wire();
	// DISP(i); DISP(Tab_Fil(i).Mesh_Wire);
      }
    
    // Maintenant je concatene les maillages de fil obtenus
    
    for(int i=0; i<nb_fils; i++)
      {
	// DISP(Tab_Fil(i).Mesh_Wire);
	Append(Mesh_Glob_Wire, Tab_Fil(i).Mesh_Wire);
      }
    
    
    
    GridInterpolation<Dimension2>& loc = this->all_points_display;
    loc.Append(Mesh_Glob_Wire, num);
    loc.LocalizePoints(this->mesh);
    
    // DISP(Mesh_Glob_Wire); DISP(loc.CoorInterp); DISP(loc.ElementInterp);
    Vector<Vector<double> > val;
    Matrix2_2 mat;
     
    // boucle sur tous les fils
    taille = 0; int old_taille = taille;
    for (int i = 0; i < nb_fils; i++)
      {
	
	Real_wp epsilon = Tab_Fil(i).GetEps();
	Real_wp val_log(0);
	


	VectReal_wp val_logphi(this->nodl); val_logphi.Zero();
	//Real_wp rayon_troncat = Tab_Fil(i).Rayon_Fct_Wire;
	
	
	old_taille = taille;
	taille += Tab_Fil(i).GetPrec_Mesh_Wire();
	
	// boucle sur les points d'integration concernes
	for(int j = old_taille; j < taille; j++)
	  { 
	    // calcul de val: phi_i(x_j) pour tout i, et x_j = loc.CoorInterp(num(j))
	    int num_elem = loc.ElementInterp(num(j));
	    this->ComputeValuesPhi(loc.CoorInterp(num(j)), pt_glob, val, mat, num_elem);

	    if (num_elem < 0)
	      {
		cout<<"le fil n'est pas dans le maillage, changez la position du fil"<<endl;
		exit(0);
	      }
	    
	    VectReal_wp val_interac(nb_fils);
	    // boucle sur toutes les fonctions de base
	    for (int k = 0; k < val.GetM(); k++)
	      {
		// fonctions de base 2-D
		int num_glob1 = this->mesh.GetNumberDof(num_elem, k);

		// fonction de base du fil : num_ddl_fil(i)
		// interaction ?

		Real_wp poids = 2.0*pi_wp/(Tab_Fil(i).GetPrec_Mesh_Wire());
		val_interac(i) = val(k)(0)*poids;
		
		
		
		//val_logphi(num_glob1) 
                // += -val_interac(i)*exp(-10*pow(epsilon/rayon_troncat,2));		
		val_logphi(num_glob1) += -val_interac(i);

		
		// DISP(num_glob1); DISP(num_ddl_fil(i)); DISP(val_interac(i));
		mat_sp.AddInteraction(num_glob1, num_ddl_fil(i), val_interac(i));
		mat_sp.AddInteraction(num_ddl_fil(i), num_glob1, val_interac(i));
	      }
	 
	    
	    
	  }
	
	if (Tab_Fil(0).Fct_Add == GlobalVariables::FCT_WIRE_AUCUNE)
	  return;
	


	////////////////////////////////////////////////////////////////////
	////    On prend en compte la fonction de base additionnelle   /////
	////////////////////////////////////////////////////////////////////


	Globatto<Dimension2> gauss;

	// points de Gauss tensorises pour les quads
	int type_quadrature = gauss.QUADRATURE_GAUSS;

	// points de quadrature pour les triangles
	if (this->mesh.nb_triangles > 0)
	  type_quadrature = gauss.QUADRATURE_TRIANGLE_GAUSS;
	
	// on construit les points
	gauss.ConstructQuadrature(Tab_Fil(i).Order_Quad_Wire, type_quadrature);


	// Boucle sur tous les elements du maillage
	Vector<bool> Element_Intersec(this->mesh.nb_elt); Element_Intersec.Fill(false);

	Vector<bool> Dof_Intersec(this->nodl); Dof_Intersec.Fill(false);
	
	
	for (int l = 0; l < this->mesh.nb_elt; l++)
	  {
	    int nb_dof_loc = this->GetNbLocalDof(l);
	    VectR_N points;
	    this->mesh.GetDofsElement(l, points, *(this->FaceBasis) );
	    for (int m1 = 0; m1 < nb_dof_loc; m1++)
	      {
		if (points(m1).Distance(Tab_Fil(i).GetCtr()) < Tab_Fil(i).Rayon_Fct_Wire)
		  {

		    Dof_Intersec(this->mesh.GetNumberDof(l,m1)) = true;
		    Element_Intersec(l) = true;		  
		  
		  }
		
	      }	    	    
	    			
	    
	  }
	
	  ///////////////////////////////////////////////////////////////
	  ///// Quadrature sans traitement special de la singularite
	  //////////////////////////////////////////////////////////////

       


	for (int l = 0; l < this->mesh.nb_elt; l++)
	    if (Element_Intersec(l))
	    {
	      int nb_dof_loc = this->GetNbLocalDof(l);
	      this->mesh.GetVerticesElement(l,s);
	      this->FjElem(s,PointsElem,this->mesh,l);
	
	      
	      
	      // boucle sur les points de quadrature
	      for(int n1 = 0; n1< (gauss.Points2D).GetLength(); n1++)
		{
		  this->Fj(s, PointsElem, gauss.Points2D(n1),pt_glob,this->mesh,l);
		  R_N diff;
		  
		  
		  
		  Subtract(pt_glob, Tab_Fil(i).GetCtr(), diff);
		  Real_wp radius = Norm2(diff);
		  
		  
		 
		  if (radius > epsilon)
		    {
		      
		      this->DFj(s,PointsElem,gauss.Points2D(n1),dfj,this->mesh,l);
		      Inverse(dfj, inv_dfj);
		      Real_wp jacob = Det(dfj);

		      Real_wp poids = gauss.Weights2D(n1)*jacob;
		      
		      this->ComputeValuesPhi(gauss.Points2D(n1), pt_glob, phi, inv_dfj, l);
		      this->ComputeValuesGradientPhi(gauss.Points2D(n1), gradphi, inv_dfj, l);
	
		      VectReal_wp ValuePhi(nb_dof_loc);
		      VectR_N GradientPhi(nb_dof_loc);

		      for (int m1 = 0; m1 < nb_dof_loc; m1++)
			{
			  ValuePhi(m1) = phi(m1)(0);
			  for (int p = 0; p < this->dim_N; p++)
			    GradientPhi(m1)(p) = gradphi(m1)(p);
			}
		      
		      Real_wp fct_log = log(radius/epsilon);		      		
		      		      
		      R_N grad_fct_log = diff;
		      Mlt(1.0/(radius*radius), grad_fct_log); 
		      		      
		      Real_wp troncat; R_N troncat_grad;		      		      
		      
		            		      
		      //////////////////////////////////////////////////////////////
		      //// 1ere  troncature en plateau
		      /////////////////////////////////////////////////////////////
		      /*
		      
		      if(radius>2*rayon_troncat/3)
			{
			  troncat = 0;
			  troncat_grad = 0;
			  
			}
		      else if (radius<rayon_troncat/3)
			{
			  troncat = 1;
			  troncat_grad = 0;
			  
			}
		      else
			{
			  
			  Real_wp val_int1 = (3*radius-rayon_troncat)/rayon_troncat;
			  Real_wp val_int2 = 18/rayon_troncat*(pow(val_int1,2)-val_int1);
			  troncat = 2*pow(val_int1,3)-3*pow(val_int1,2)+1;
			  troncat_grad = diff;
			  Mlt(1.0/radius, troncat_grad);
			  Mlt(val_int2, troncat_grad);
		
			  }
		      
		      */      		      
		      //////////////////////////////////////////////////////////////
		      //// 2eme troncature en plateau
		      /////////////////////////////////////////////////////////////
		      /*
		      
		      if(radius>2*rayon_troncat/3)
			{
			  troncat = 0;
			  troncat_grad = 0;
			  
			}
		      else if (radius<rayon_troncat/3)
			{
			  troncat = 1;
			  troncat_grad = 0;
			  
			}
		      else
			{
			  
			  Real_wp val_int1 = (3*radius-rayon_troncat)/rayon_troncat;
			  Real_wp val_int2 = 3/rayon_troncat*
                          (-5*6*pow(val_int1,4)+15*4*pow(val_int1,3)-10*3*pow(val_int1,2));
			  troncat = 1-6*pow(val_int1,5)+15*pow(val_int1,4)-10*pow(val_int1,3);
			  troncat_grad = diff;
			  Mlt(1.0/radius, troncat_grad);
			  Mlt(val_int2, troncat_grad);
		
			  }
		      */
		      //////////////////////////////////////////////////////////////
		      ////  troncature par une gaussienne
		      /////////////////////////////////////////////////////////////

		      /*
		      troncat = exp(-10*pow(radius/rayon_troncat,2));
		      troncat_grad = diff;
		      Mlt(-20*troncat/pow(rayon_troncat,2),troncat_grad);
		      */
		      
		      this->Tab_Fil(i).Compute_Troncat( radius, troncat);
		      this->Tab_Fil(i).Compute_Troncat_Grad( diff,troncat_grad);
		      
		      //interaction du log avec lui meme
		      val_log += -this->omega2*poids*troncat*troncat*fct_log*fct_log
		      	+ DotProd(troncat_grad, troncat_grad)*fct_log*fct_log*poids;
		      
		      
		      
		      // boucle sur les dofs de chaque element
		      // Interaction du log avec les autres fct de forme
		      for (int m1 = 0; m1 < nb_dof_loc; m1++)
			{
			  int num_dof = this->mesh.GetNumberDof(l,m1);
			  
			  
			  
			   
			  val_logphi(num_dof) += poids
                            *(-this->omega2*ValuePhi(m1)*troncat*fct_log +
                              DotProd(GradientPhi(m1),troncat_grad)*fct_log
                              -ValuePhi(m1)*DotProd(grad_fct_log,troncat_grad));
			  
			  
			  
			  
			  
			  
  
			}  	
		      
		    }
		}
	      
	      
	      
	    } // end loop on elements
	
	
	// ajout de l'element sur la diagonale
	
	mat_sp.AddInteraction(num_ddl_fctfil(i), num_ddl_fctfil(i), val_log);
	
	// et les autres
	for (int m = 0; m < this->nodl; m++)
	  if (val_logphi(m) != 0)
	    {
	      
	      mat_sp.AddInteraction(m, num_ddl_fctfil(i), val_logphi(m));
	      mat_sp.AddInteraction(num_ddl_fctfil(i), m, val_logphi(m));
	    }
	
	
	////////////////////////////////////////////////////////////
	///    Fin de la partie fonction de base additionnelle   ///
	////////////////////////////////////////////////////////////
	
	
      }
    
  }


  //! cut-off function
  void Wire<Dimension2>::Compute_Troncat(const Real_wp& radius, Real_wp& troncat)
  {
  
    Real_wp rayon_troncat = this->Rayon_Fct_Wire;

    if(this->Fct_Add==GlobalVariables::FCT_WIRE_GAUSSIENNE)
	troncat = exp(-10*pow(radius/rayon_troncat,2));	      
     
    else if(this->Fct_Add==GlobalVariables::FCT_WIRE_PLATEAU_H2)
      {
	if(radius>2*rayon_troncat/3)
	  {
	    troncat = 0;
	  }
	else if (radius<rayon_troncat/3)
	  {
	    troncat = 1; 
	  }
	else
	  {
	    Real_wp val_int1 = (3*radius-rayon_troncat)/rayon_troncat;
	    troncat = 2*pow(val_int1,3)-3*pow(val_int1,2)+1;
	  }	      
	
      }
    else if(this->Fct_Add==GlobalVariables::FCT_WIRE_PLATEAU_H3)
      {
	
	if(radius>2*rayon_troncat/3)
	  {
	    troncat = 0;					  
	  }
	else if (radius<rayon_troncat/3)
	  {
	    troncat = 1;					  
	  }
	else
	  {			  
	    Real_wp val_int1 = (3*radius-rayon_troncat)/rayon_troncat;
	    troncat = 1-6*pow(val_int1,5)+15*pow(val_int1,4)-10*pow(val_int1,3);
	  }			
      }
    else
      {
	// A COMPLETER EN ENVOYANT UNE ERREUR
      }
  }

  
  //! gradient of cut-off function
  void Wire<Dimension2>::Compute_Troncat_Grad(const R2& diff,R2& troncat_grad)
  {
    Real_wp rayon_troncat = this->Rayon_Fct_Wire;
    Real_wp radius = Norm2(diff);
    if(this->Fct_Add==GlobalVariables::FCT_WIRE_GAUSSIENNE)
      {
	Real_wp troncat = exp(-10*pow(radius/rayon_troncat,2));
	troncat_grad = diff;
	Mlt(-20*troncat/pow(rayon_troncat,2),troncat_grad);
      }
    else if(this->Fct_Add==GlobalVariables::FCT_WIRE_PLATEAU_H2)
      {	      
	
	if(radius>2*rayon_troncat/3)
	  {			 
	    troncat_grad = 0;			  
	  }
	else if (radius<rayon_troncat/3)
	  {
	    troncat_grad = 0;			  
	  }
	else
	  {			  
	    Real_wp val_int1 = (3*radius-rayon_troncat)/rayon_troncat;
	    Real_wp val_int2 = 18/rayon_troncat*(pow(val_int1,2)-val_int1);
	    troncat_grad = diff;
	    Mlt(1.0/radius, troncat_grad);
	    Mlt(val_int2, troncat_grad);		
	  }
      }
    else if(this->Fct_Add==GlobalVariables::FCT_WIRE_PLATEAU_H3)
      {
	if(radius>2*rayon_troncat/3)
	  {		
	    troncat_grad = 0;			  
	  }
	else if (radius<rayon_troncat/3)
	  {
	    troncat_grad = 0;			  
	  }
	else
	  {
	    
	    Real_wp val_int1 = (3*radius-rayon_troncat)/rayon_troncat;
	    Real_wp val_int2 = 3/rayon_troncat
              *(-5*6*pow(val_int1,4)+15*4*pow(val_int1,3)-10*3*pow(val_int1,2));
            
	    troncat_grad = diff;
	    Mlt(1.0/radius, troncat_grad);
	    Mlt(val_int2, troncat_grad);
	    
	  }

      }
    else
      {
	// A COMPLETER EN ENVOYANT UNE ERREUR
      }
  }
  
  
  

  
  
#endif
  
}

#define MONTJOIE_FILE_MODEL_WIRES_CXX
#endif
