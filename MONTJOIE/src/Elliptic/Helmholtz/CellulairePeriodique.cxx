#ifndef MONTJOIE_FILE_CELLULAIRE_PERIODIQUE_CXX

namespace Montjoie
{
  
  template<class TypeEquation>
  VarHelmholtz_CellulairePeriodique<TypeEquation>::
  VarHelmholtz_CellulairePeriodique() : VarHelmholtz_Cplx<TypeEquation>()
  {
    troncature = 10;
    R0 = 1;
    PeriodeY = 1 ;
    PeriodeX = 1.0;
    SymmetricCase = false; 
    source_valeur_absolue = false ;
    ParametrePenalisation = -1e-8;
    CasDirichlet = false ;
    eps_max_infty = 1.0;
    mu_max_infty = 1.0;
    sigma_max_infty = 0.0;
    mu_infty = 1.0;
    file_name_constant = "constantes.txt";
  }
  
  
  // sucharge de SetInputData pour prendre en compte les nouveaux champs R0 et 
  // troncature (troncature de l operateur Dirichlet To Neumann) 
  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique <TypeEquation>
  ::SetInputData(const string & description_field,  const VectString& parameters)
  {
    VarHelmholtz_Cplx<TypeEquation>::SetInputData(description_field, parameters);
    
    if (!description_field.compare("RayonMoyen"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHelmholtz_CellulairePeriodique" << endl;
	    cout << "RayonMoyen needs more parameters, for instance :" << endl;
	    cout << "RayonMoyen = R0" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	// cout<<parameters(0)<<endl ;
	R0 = to_num<Real_wp>(parameters(0));
	DISP(R0) ;
      }
    else if (!description_field.compare("Troncature"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHelmholtz_CellulairePeriodique" << endl;
	    cout << "Troncature needs more parameters, for instance :" << endl;
	    cout << "Troncature = Nmax" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        //cout<<parameters(0)<<endl ;
        troncature = to_num<int>(parameters(0));
        DISP(troncature) ;
      }
    else if (!description_field.compare("PenalisationLaplacien"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHelmholtz_CellulairePeriodique" << endl;
	    cout << "PenalisationLaplacien needs more parameters, for instance :" << endl;
	    cout << "PenalisationLaplacien = epsilon" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        ParametrePenalisation = to_num<Real_wp>(parameters(0));
      }
    else if (!description_field.compare("CelluleSymetrique"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHelmholtz_CellulairePeriodique" << endl;
	    cout << "CelluleSymetrique needs more parameters, for instance :" << endl;
	    cout << "CelluleSymetrique = FALSE" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        //cout<<parameters(0)<<endl ;
        if (!parameters(0).compare("FALSE"))
          SymmetricCase = false;
        else if (!parameters(0).compare("false"))
          SymmetricCase = false;
        else if (!parameters(0).compare("TRUE"))
          SymmetricCase = true;
        else if (!parameters(0).compare("true"))
          SymmetricCase = true;
        else
          {
            cout << " probleme : CelluleSymetrique = true ou  CelluleSymetrique =false "  << endl ;
            abort(); 
          }
        
        DISP(SymmetricCase) ;
        
      }
    else if (!description_field.compare("CasDirichlet"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHelmholtz_CellulairePeriodique" << endl;
	    cout << "CasDirichlet needs more parameters, for instance :" << endl;
	    cout << "CasDirichlet = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (!parameters(0).compare("YES"))
          CasDirichlet = true ;
        else
          CasDirichlet = false ;
      }
    
  }


  // surcharge de InitIndices
  // prise en compte d un nouveau champ ref_rho2 
  // ref_rho est nul (artifice pour resoudre div(mu grad(U)) =f
  // ref_rho2 contient les valeurs de rho 
  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique<TypeEquation>::InitIndices(int n)
  {
    if (n <= 1)
      {
	cout<<" ATTENTION VOUS N'AVEZ MIS QU'UN SEUL MILIEU 0!!!!!!!!!"<<endl;
	// throw WrongIndices("Number of domains should be greater than 2");
      }
    
    VarHelmholtz_Cplx<TypeEquation>::InitIndices(n);
    
    this->ref_epsilon_max.Reallocate(n+1);
    this->ref_mu_max.Reallocate(n+1);
    this->ref_sigma_max.Reallocate(n+1);
    
    Complex_wp one(1, 0);
    for (int i = 0; i <= n; i++)
      {
	this->ref_epsilon_max(i).SetConstant(one);
	this->ref_mu_max(i).SetIdentity();
	this->ref_sigma_max(i).Zero();
      }
  }
  
  
  // surcharge de setIndices
  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique<TypeEquation>
  ::SetIndices(int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_epsilon_max.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_epsilon_max.GetM() << endl;
        abort();
      }
    
    // l'utilisateur donne epsilon, mu et sigma (pour Maxwell)
    this->ref_epsilon_max(i).SetInputData(nb, parameters, parameters(0));
    this->ref_mu_max(i).SetInputData(nb, parameters, parameters(0));
    this->ref_sigma_max(i).SetInputData(nb, parameters, parameters(0));
    
    // on met rho = penal ou 0 (avec penal une petite valeur, par exemple 1e-6)
    if ( SymmetricCase == false)
      {
        this->ref_rho(i).SetIdentity();
        if ( this->CasDirichlet) 
          this->ParametrePenalisation = 0;
        
        this->ref_rho(i).Mlt(ParametrePenalisation/this->GetSquareOmega()) ;
        this->ref_invRho(i) = this->ref_rho(i);
        this->ref_invRho(i).GetInverse();
      }
    else
      {
	Complex_wp zero(0, 0);
        this->ref_rho(i).SetConstant(zero);
        this->ref_invRho(i).SetConstant(zero);
      }
    
  }
  

  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
   */
  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_epsilon_max.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_epsilon_max.GetM() << endl;
        abort();
      }

    if (name_media == "epsilon")
      this->ref_epsilon_max(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "mu")
      this->ref_mu_max(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "sigma")
      this->ref_sigma_max(i).SetInputData(nb, parameters, parameters(0));
    else
      {
        cout << "Unknown media : " << name_media << endl;
        abort();
      }
  }
  
    
  //! returns the name associated with the physical index num
  template<class TypeEquation>
  string VarHelmholtz_CellulairePeriodique<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("epsilon");
      case 1: return string("mu");
      case 2: return string("sigma");
      }

    return string();
  }  
  
  
  //calcul des constantes B1 B2
  // B1 = (omega^2  R0 / \mu_inf) (int_{omega 1/2}  rho_inf -rho)
  // B2 = (1/mu_inf) (int_{omega 1/2} (mu_inf/R0) -mu(dU/dZ +1/R0))
  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique<TypeEquation>
  ::CalculB1B2(const VectComplexe & U0, Complexe& B1, Complexe& B2)
  {
    Complexe sum, sum2, sum3;
    sum = 0 ;
    sum2 = 0 ;
    sum3 = 0 ;
    VectComplexe Uquad, grad_Uquad;
    TinyVector<VectComplexe, 1> Uloc;
    typename Dimension::MatrixN_N dfjm1;
    TinyVector<Complexe, Dimension::dim_N> hat_grad_u, grad_u;
    Real_wp jacob(0), poids(0);
	
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	// on repere sur quel materiaux est situé l'élément
	// l integral n est faite que sur les materiaux non references par 1 (omega 1/2)
	int ref = this->mesh.Element(i).GetReference();
	if (ref <= 1)
	  continue;
	
	const ElementReference<Dimension, 1>& Fb = this->GetReferenceElementH1(i);
	int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
	TinyMatrix<Complexe, Symmetric, Dimension::dim_N, Dimension::dim_N> mu;
	Complexe rho2 ;
	
	// on recupere U0 localement
	this->GetLocalUnknownVector(U0, i, Uloc);
	
	// calcul de grad u
	grad_Uquad.Reallocate(Dimension::dim_N*nb_points_quadrature);
	if (Fb.UseQuadraturePointsForRh())
	  {
	    Uquad.Reallocate(nb_points_quadrature);
	    Fb.ApplyChTranspose(Uloc(0), Uquad);
	    Fb.ApplyRhQuadratureTranspose(Uquad, grad_Uquad);
	  }
	else
	  Fb.ApplyRhTranspose(Uloc(0), grad_Uquad);
	
	// omega 1/2 + =  omega 1/2 inter R>0 (on a utilise les symetries du probleme)
	// terme a calculer : (1/mu_inf) (int_{omega 1/2 +} (mu_inf/R0) -mu(dPhi_j/dZ +1/R0))
	// int_{T} dPhi_j/dZ = int_{T^} grad{Phi_j}. [1 0] 
	//                   = int_{Tchapeau}  det(BT) gard(Phi_j_chapeau) . B^-1 [1 0]
	// [1 0 ] =temp1
	// B^-1 = matriceInverse
	// boucle sur les points de quadrature
	for (int k = 0; k < nb_points_quadrature; k++)
	  {
	    // coordonnnes du point k
	    //on recupere les valeurs de mu et de rho2 ....(oui puisque rho =0)
	    mu = this->ref_mu(ref).GetCoefficient(*this, i, k);
	    Complexe epsilon = this->ref_epsilon_max(ref).GetCoefficient(*this, i, k);
	    Complexe sigma = this->ref_sigma_max(ref).GetCoefficient(*this, i, k);
	    rho2 = epsilon + Iwp*sigma/this->GetOmega();
	    poids = Fb.WeightsND(k);
	    
	    if (this->mesh.IsElementAffine(i))
	      {
		jacob = this->Glob_jacobian(i)(0);
		dfjm1 = this->Glob_DFjm1(i)(0);
	      }
	    else
	      {
		dfjm1 = this->Glob_DFjm1(i)(k);
		jacob = this->Glob_jacobian(i)(k)/poids;
	      }
	    
	    CopyVector(grad_Uquad, k, hat_grad_u);
	    MltTrans(dfjm1, hat_grad_u, grad_u);
	    
	    //DISP(grad_u);
	    // terme contribuant int_{T}mu * dPhi_j/dZ
	    sum += grad_u(0)* mu(0,0) * poids; 
	    // terme contribuant a B1 mais ne dependant pas de U (int_T mu)
	    sum2 += mu(0,0)* poids/(this->R0)*jacob;
	    // terme conttribuant a B1
	    sum3 +=  rho2 * poids*jacob;
	    
	  } //end loop points de quadrature
	
      } //end loop elements du maillage
    
    // on multiplie par 2 parce que l on a fait les calcul sur omega + 1/2 et non pas omega 1/2 
    // par parite de la solution il suffit de multiplier par 2 
    DISP(sum) ;
    DISP(sum2) ;
    DISP(sum3); 
    
    sum = (-2.0 * (sum +sum2) + this->mu0/this->R0) / this->mu0 ;
    
    B2 = sum ;
    // meme explication pour la multiplication par 2 
    B1  = ( this->eps_max_infty  - 2. * sum3 )* this->omega2 * this->R0/ this->mu0  ;
    
  }


  /////////////////////////////////////////////////////////////////////
  // Pour le 3D
  /////////////////////////////////////////////////////////////////////
  


  // D11 = (1/mu_inf T ) (int_{omega 1/2 +} (mu_inf) (dPhi_j/dX +1)) -1
  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique<TypeEquation>::
  CalculD11D22A3D(const VectComplexe & U0, Complexe& D11, Complexe& D22, Complexe& A3D)
  {
    cout << "calcul de D11" << endl ;
    Complex_wp sumD11,  sumD22, sum2D11D22, sumA3D;
    sumD11 = 0 ;
    sumD22 = 0 ;
    sum2D11D22 = 0 ;
    sumA3D = 0;
    TinyVector<Vector<Complexe>, 1> Uloc;
    VectComplexe GradUloc, Uquad;
    typename Dimension::MatrixN_Nsym_Complex_wp mu;
    Complex_wp rho2 ;
    TinyVector<Complexe, Dimension::dim_N> GradChap, GradU;
    MatrixN_N DFjm1 ;  
    Real_wp jacob_weighted;
    
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	int ref = this->mesh.Element(i).GetReference();
        const ElementReference<Dimension, 1>& Fb = this->GetReferenceElementH1(i);
	
	// on repere sur quel materiaux est situé l'élément
	// l integral n est faite que sur les materiaux non references par 1 (omega 1/2)
	if (ref <= 1)
	  continue;
	
	// omega 1/2 + =  omega 1/2 inter R>0 (on a utilise les symetries du probleme)
	// terme a calculer : (1/mu_inf T ) (int_{omega 1/2 +} (mu_inf) (dPhi_j/dX +1)) -1
	int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
	GradUloc.Resize(Dimension::dim_N * nb_points_quadrature);
	this->GetLocalUnknownVector(U0, i, Uloc);
	
	if (Fb.UseQuadraturePointsForRh())
	  {
	    Uquad.Reallocate(nb_points_quadrature);
	    Fb.ApplyChTranspose(Uloc(0), Uquad);
	    Fb.ApplyRhQuadratureTranspose(Uquad, GradUloc);
	  }
	else
	  Fb.ApplyRhTranspose(Uloc(0), GradUloc);
	
	for (int k = 0; k < nb_points_quadrature; k++)
	  {
	    Real_wp poids = Fb.WeightsND(k);
	    
	    if (this->mesh.IsElementAffine(i))
	      {
		DFjm1 = this->Glob_DFjm1(i)(0);
		jacob_weighted = this->Glob_jacobian(i)(0)*poids;
	      }
	    else
	      {
		DFjm1 = this->Glob_DFjm1(i)(k);
		jacob_weighted = this->Glob_jacobian(i)(k);
	      }
	    
	    CopyVector(GradUloc, k, GradChap);
	    
	    MltTrans(DFjm1, GradChap, GradU) ;
	    // coordonnnes du point k
	    //on recupere les valeurs de mu et de rho2 ....(oui puisque rho =0)
	    mu = this->ref_mu(ref).GetCoefficient(*this, i, k);
	    
	    //cout << "apres le mu" << endl ;
	    // terme contribuant int_{T}mu * dPhi_j/dX
	    sumD11 += GradU(0) * mu(0,0) * poids;
	    
	    // terme contribuant a D11  mais ne dependant pas de U (int_T mu)
	    sum2D11D22 += mu(0,0)* jacob_weighted;
	    
	    // terme contribuant int_{T}mu * dPhi_j/dY
	    sumD22 += GradU(1)* mu(0,0) * poids; 
	    //  terme contribuant int_{T} dPhi_j/dZ
	    sumA3D += GradU(2) * poids; 
	  }   
	
      } // end loop elements du maillage
    
    // dans le cas symetrique on multiplie par 2 parce que l on a fait les calcul
    // sur omega + 1/2 et non pas omega 1/2 
    // par parite de la solution il suffit de multiplier par 2 
    DISP( sum2D11D22) ;
    DISP(sumD11) ;
    DISP(sumD22) ;
    DISP(sumA3D) ;
    if ( SymmetricCase )
      {
        sumD11 =   2.0 * (sumD11 +sum2D11D22) / (this->mu0 * this->PeriodeY * this->PeriodeX) ;
        sumD22 =   2.0 * (sumD22 +sum2D11D22) / (this->mu0 * this->PeriodeY * this->PeriodeX) ;
        sumA3D =   2.0 * (sumA3D ) / (this->mu0 * this->PeriodeY * this->PeriodeX) ;
      }
    else 
      {
        sumD11 =    (sumD11 +sum2D11D22) / (this->mu0 * this->PeriodeY * this->PeriodeX) ;
        sumD22 =    (sumD22 +sum2D11D22) / (this->mu0 * this->PeriodeY * this->PeriodeX) ;
        sumA3D =    (sumA3D ) / (this->mu0 * this->PeriodeY * this->PeriodeX) ;
      }
    
    sumD11 += -1.0 ;
    sumD22 += -1.0 ;
    sumA3D += -1.0 ;
    
    D11 = sumD11  ;
    D22 = sumD22  ;
    A3D = sumA3D  ;
  }
  


  //calcul de la constante A
  // A = 2/R0 * ( int_{R =Re} U(theta,Re)  dtheta - Re)  
  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique<TypeEquation>
  ::CalculA(const VectComplexe & U0, Complexe& A, bool saut)
  {
    Real_wp Rexterieur ;
    Complexe sum; SetComplexZero(sum);
    
    SetPoints<Dimension> PointsElem1;
    SetMatrices<Dimension> MatricesElem1;
    TinyVector<VectComplexe, 1> Uloc;
    
    Rexterieur = 0  ;
    // premiere boucle pour trouver Re
    for (int e = 0; e < this->mesh.GetNbBoundaryRef(); e++)
      {
	int ref_line_e =  this->mesh.BoundaryRef(e).GetReference();
	if (this->mesh.GetBoundaryCondition(ref_line_e)== BoundaryConditionEnum::LINE_ABSORBING )
	  {
            int iquad_e = this->mesh.BoundaryRef(e).numElement(0);
            int num_loc_e = this->mesh.Element(iquad_e).GetPositionBoundary(e);
                        
	    int n1 = this->mesh.Element(iquad_e).numVertex(num_loc_e);
	    Rexterieur = abs(this->mesh.Vertex(n1)(1));
            
	    break ;
	  }
      }
    
    // on remplit petit a petit  int_{R =Re} U(theta,Re)
    // il faut donc calculer int_{arrete de Gamma+} phi_j(theta,Re) * U_j
    // boucle sur les aretes referencees
    bool soustraction ;
    for (int e = 0; e < this->mesh.GetNbBoundaryRef(); e++)
      {
	// on recupere la reference de l arete
        int ref_line_e =this->mesh.BoundaryRef(e).GetReference();
	
        // si e fait partie des aretes de gamma +
	if (this->mesh.GetBoundaryCondition(ref_line_e) == BoundaryConditionEnum::LINE_ABSORBING )
	  {
            int iquad_e = this->mesh.BoundaryRef(e).numElement(0);
            int num_loc_e = this->mesh.Element(iquad_e).GetPositionBoundary(e);
            int num_sommet = this->mesh.Element(iquad_e).numVertex(num_loc_e);
	    Real_wp coordonnees = (this->mesh.Vertex(num_sommet)(1));
            if ((coordonnees < 0) && saut )
              soustraction = true  ;
            else
              soustraction = false ;
	    
	    const ElementReference<Dimension, 1>& Fb = this->GetReferenceElementH1(iquad_e);
	    // poids et elements de surfaces
	    Real_wp  poidsm,  dsm;
	    // coordonnees des sommets 
	    VectR_N s1;
	    R_N sommet ;
	    Complexe sum2 ;
	    sum2 = 0 ;
	    Complex_wp  phi_im;
	    // sommet de l element iquad_e
	    this->mesh.GetVerticesElement(iquad_e, s1);
	    
	    // on recupere les transformations elementaires
	    Fb.FjElem(s1, PointsElem1, this->mesh, iquad_e);
	    Fb.DFjElem(s1, PointsElem1, MatricesElem1, this->mesh, iquad_e);
	    
	    Fb.FjSurfaceElem(s1, PointsElem1, this->mesh, iquad_e, num_loc_e);
	    Fb.DFjSurfaceElem(s1, PointsElem1, MatricesElem1, this->mesh, iquad_e, num_loc_e);
	    
	    // on recupere U0 localement
	    this->GetLocalUnknownVector(U0, iquad_e, Uloc);
	    
	    // evaluation de u sur les points de quadrature
	    int nb_quad = Fb.GetNbQuadBoundary(num_loc_e);
	    VectComplexe Uquad(nb_quad);
	    Fb.ApplyShTranspose(num_loc_e, Uloc(0), Uquad);
	    
	    // boucle sur les poinrs de quadrature    
	    for( int m = 0; m < nb_quad; m++)
	      {
		poidsm = Fb.WeightsQuadratureBoundary(m, num_loc_e) ;
		dsm =  MatricesElem1.GetDsQuadratureBoundary(m);
		
		sum2 +=  poidsm * dsm * Uquad(m) ;
	      }
	    if (soustraction)
	      sum -= sum2 ;
	    else
	      sum +=sum2 ;
	  } //fin if line absorbing
        
      } //fin boucle sur les aretes
    
    if (this-> CasDirichlet)
      {
        if ( this->source_valeur_absolue) 
          {
            if (saut)
              A = (sum)/this->R0 ;
           
            else
              A = 0.5 * (sum - 2 * Rexterieur)/this->R0; 
          }
        else
          {
            if (saut)
              A = (sum - 2 * Rexterieur)/this->R0 ;
            else
              A = 0.5 * (sum )/this->R0; 
          }
      }
    else
      {
        A = (2.0 * (sum - Rexterieur))/this->R0 ;
      }
    
    DISP( this->source_valeur_absolue) ;
    DISP(saut) ;
    DISP(A) ;
 
    
  }
  
  
  //   // D22 = (1/mu_inf T ) (int_{omega 1/2 +} (mu_inf) (dPhi_j/dY +1)) -1
  //   template<class TypeEquation>
  //   void VarHelmholtz_CellulairePeriodique<TypeEquation>
  // ::CalculD22(const VectComplexe & U0)
  //   {
  
  //     Vector<Vector<TinyVector<Complexe,1> > > Unodal;
  //     Vector<Vector<TinyVector<Complexe, 3> > > GradNodal;
  //     //var_leaf.ComputeU_GradU_Nodes( U0, Unodal, GradNodal, false, true);
  //     var_leaf.ComputeNodalUgradU( U0, Unodal, GradNodal, false, true);
  //     Complex_wp sum, sum2;
  //     int ref ;
  //     sum = 0 ;
  //     sum2 = 0 ;
  //     // omega 1/2 + =  omega 1/2 inter R>0 (on a utilise les symetries du probleme)
  //     // terme a calculer : (1/mu_inf T ) (int_{omega 1/2 +} (mu_inf) (dPhi_j/dX +1)) -1
  //    cout<<"Calcul de D22 "<<endl;
  //     int nb_points_quadrature;
  //     typename Dimension::MatrixN_Nsym_Complex_wp mu;
  //     Complex_wp rho2 ;
  
  //     for (int i = 0; i < this->mesh.GetNbElt(); i++)
  //       {
    //     const ElementReference<Dimension, 1>& Fb = this->GetReferenceElement(i);
  // 	ref = this->mesh.Element(i).GetReference();
  
  // 	// on repere sur quel materiaux est situé l'élément
  // 	// l integral n est faite que sur les materiaux non references par 1 (omega 1/2)
  // 	if (ref > 1)
  // 	  {
  
  
  // 	    nb_points_quadrature = this->GetNbPointsQuadratureInside(i);
  
  //             //boucle sur les points de quadrature
  //             // on separe le cas affine de l autre cas .... GRRR 
  //             if  (this->mesh.IsElementAffine(i)==true)
  //               {
  
  //                 for (int k = 0; k < nb_points_quadrature; k++)
  //                   {
  //                     // coordonnnes du point k
  //                     //on recupere les valeurs de mu et de rho2 ....(oui puisque rho =0)
  //                     mu = this->ref_mu(ref).GetCoefficient(var_leaf, i, k);
  
  
  //                     // terme contribuant int_{T}mu * dPhi_j/dY
  //                    /sum += GradNodal(i)(k)(1)* mu(0,0)
  //         * this->Glob_jacobian(i)(0)*Fb.WeightsND(k); 
  
  //                     // terme contribuant a D11  mais ne dependant pas de U (int_T mu)
  //                     sum2 +=mu(0,0)* this->Glob_jacobian(i)(0) *Fb.WeightsND(k) ;
  
  //                   }
  //               }
  //             else
  //               {
  //                 for (int k = 0; k < nb_points_quadrature; k++)
  //                   {
  //                     //cout << " non afine" << endl ;
  //                     // coordonnnes du point k
  //                     //on recupere les valeurs de mu et de rho2 ....(oui puisque rho =0)
  //                     mu = this->ref_mu(ref).GetCoefficient(var_leaf, i, k);
  
  
  //                     // terme contribuant int_{T}mu * dPhi_j/dY
  //                      sum += GradNodal(i)(k)(1)* mu(0,0) * this->Glob_jacobian(i)(k); 
  
  //                     // terme contribuant a D11  mais ne dependant pas de U (int_T mu)
  //                     sum2 +=mu(0,0)* this->Glob_jacobian(i)(k);
  
  //                   }   
  //               }
  
  //             //end loop points de quadrature
  //         // }//end loop degres de liberte
  //           }//end if ref >1
  //       }//end loop elements du maillage
  //     // om multiplie par 2 parce que l on a fait les calcul sur
  //      omega + 1/2 et non pas omega 1/2 
  //     // par parite de la solution il suffit de multiplier par 2 
  //      cout << "sum doit etre negative " << endl ;
  //     DISP(sum) ;
  //     DISP(sum2) ;
  
  //     if (SymmetricCase == true)
  //       sum =   2.0 * (sum +sum2) / (this->mu0 * this->PeriodeY * this->PeriodeX) ;
  //     else
  //       sum =  (sum +sum2) / (this->mu0 * this->PeriodeY * this->PeriodeX) ;
  
  //     sum += -1.0 ;
  
  
  
  //   }
  
  
  
  
  
  //   // A3D = (1/mu_inf T ) (int_{omega 1/2 +} (dPhi_j/dZ) -1
  //   template<class TypeEquation>
  //   void VarHelmholtz_CellulairePeriodique<TypeEquation>
  // ::CalculA3D(const VectComplexe & U0)
  //   {
  //     // TinyVector<Matrix<Complexe>, 1> Unodal;
  //     //TinyVector<Matrix<Complexe>, 3> GradNodal;
  //     Vector<Vector<TinyVector<Complexe,1> > > Unodal;
  //     Vector<Vector<TinyVector<Complexe, 3> > > GradNodal;
  //     //var_leaf.ComputeU_GradU_Nodes( U0, Unodal, GradNodal, false, true);
  //     var_leaf.ComputeNodalUgradU( U0, Unodal, GradNodal, false, true);
  //     Complex_wp sum, sum2;
  //     int ref ;
  //     sum = 0 ;
  //     sum2 = 0 ;
  //     // omega 1/2 + =  omega 1/2 inter R>0 (on a utilise les symetries du probleme)
  //     // terme a calculer : (1/mu_inf T ) (int_{omega 1/2 +} (mu_inf) (dPhi_j/dX +1)) -1
  //    cout<<"Calcul de A3D "<<endl;
  //     int nb_points_quadrature;
  //     typename Dimension::MatrixN_Nsym_Complex_wp mu;
  //     Complex_wp rho2 ;
  
  //     for (int i = 0; i < this->mesh.GetNbElt(); i++)
  //       {
    //     const ElementReference<Dimension, 1>& Fb = this->GetReferenceElement(i);
  // 	ref = this->mesh.Element(i).GetReference();
  
  // 	// on repere sur quel materiaux est situé l'élément
  // 	// l integral n est faite que sur les materiaux non references par 1 (omega 1/2)
  // 	if (ref > 1)
  // 	  {
  
  
  // 	    nb_points_quadrature = this->GetNbPointsQuadratureInside(i);
  
  //             //boucle sur les points de quadrature
  //             // on separe le cas affine de l autre cas .... GRRR 
  //             if  (this->mesh.IsElementAffine(i)==true)
  //               {
  
  //                 for (int k = 0; k < nb_points_quadrature; k++)
  //                   {
  //                     // coordonnnes du point k
  //                     // terme contribuant int_{T} dPhi_j/dZ
  //                     sum += GradNodal(i)(k)(2) * this->Glob_jacobian(i)(0)*Fb.WeightsND(k); 
  
  
  
  //                   }
  //               }
  //             else
  //               {
  //                 for (int k = 0; k < nb_points_quadrature; k++)
  //                   {
  
  //                     // terme contribuant int_{T} dPhi_j/dZ
  //                     sum += GradNodal(i)(k)(2) * this->Glob_jacobian(i)(k); 
  
  //                   }   
  //               }
  
  //             //end loop points de quadrature
  //         // }//end loop degres de liberte
  //           }//end if ref >1
  //       }//end loop elements du maillage
  //     // om multiplie par 2 parce que l on a fait les calcul sur
  //    omega + 1/2 et non pas omega 1/2 
  //     // par parite de la solution il suffit de multiplier par 2 
  //     DISP(sum) ;
  //     if (SymmetricCase ==true )
  
  //       sum =   2.0 * (sum ) / (this->mu0 * this->PeriodeY * this->PeriodeX) ;
  //     else
  //       sum = (sum ) / (this->mu0 * this->PeriodeY * this->PeriodeX) ;
  
  //     sum += -1.0 ;
  //     this->A3D = sum  ;
  //     DISP(this->A3D) ;
  
  
  //   }
  
  

  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique <TypeEquation>
  ::AddBoundaryConditionTerms(VirtualMatrix<Real_wp>& mat_sp,
			      const GlobalGenericMatrix<Real_wp>& nat_mat,
			      int offset_row, int offset_col)
  {
  }


  // surcharge de l operateur DtN
  // double boucle sur les aretes .. un peu lent
  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique <TypeEquation>
  ::AddBoundaryConditionTerms(VirtualMatrix<Complex_wp>& mat_sp,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat,
			      int offset_row, int offset_col)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(this->comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    int order_integration = this->mesh_num.GetOrder();
    
    for (int ref = 1; ref <= this->mesh.GetNbReferences(); ref++)
      {
	if (this->mesh.GetBoundaryCondition(ref) != BoundaryConditionEnum::LINE_ABSORBING)
	  continue;
	
	// on recupere le maillage surfacique de reference ref_cond
	// et on calcule les poids d'integration
	MeshInterpolationFEM<Dimension> var_surf;
	var_surf.SetGaussQuadrature(order_integration);
	IVect ref_cond_(this->mesh.GetNbReferences()+1); ref_cond_.Fill(0);
	ref_cond_(ref) = 1;
	
	Mesh<Dimension> mesh_subdiv;
	var_surf.InitProjectionSurface(this->mesh);
	var_surf.ComputeSurfaceMesh(ref_cond_, this->mesh, mesh_subdiv, *this);
	
	// on compte le nombre de faces concernees
	int nb_faces_concerned = var_surf.GetNbBoundary();
	
	int nb_total_faces_concerned = nb_faces_concerned;
	int old_percent = 0, new_percent = 0;
	nb_faces_concerned = 0;
	
	Vector<Vector<Complex_wp> > feval_phi;
	VectComplexe feval, contrib_prime, contrib;
	// boucle principale sur les faces
	for (int e = 0; e < var_surf.GetNbBoundary(); e++)
	  {
	    int iquad_e  = var_surf.GetElementNumberOfSurface(e);
	    int num_loc_e = var_surf.GetLocalPositionOfSurface(e);
	    
	    // displays a progress bar
	    new_percent = toInteger(round(Real_wp(nb_faces_concerned)
					  /(nb_total_faces_concerned-1)*80));
	    if (rank_proc == 0)
	      if (this->print_level >= 2)
		for (int percent = old_percent; percent < new_percent; percent++)
		  { cout<<"#"; cout.flush(); }
	    
	    nb_faces_concerned++;
	    old_percent = new_percent;
	    
	    // seconde boucle sur les faces
	    for (int eprime = 0; eprime < var_surf.GetNbBoundary(); eprime++)
	      {
		int iquad_eprime  = var_surf.GetElementNumberOfSurface(eprime);
		int num_loc_eprime = var_surf.GetLocalPositionOfSurface(eprime);
		
		const ElementReference<Dimension, 1>& Fb = this->GetReferenceElementH1(iquad_e);
		const ElementReference<Dimension, 1>& Fb_prime = this->GetReferenceElementH1(iquad_eprime);
		int nb_quad_e = var_surf.GetNbPointsQuadrature(e);
		int nb_quad_eprime = var_surf.GetNbPointsQuadrature(eprime);
		feval.Reallocate(nb_quad_e);
		feval_phi.Reallocate(Fb.GetNbDof());
		contrib_prime.Reallocate(Fb_prime.GetNbDof());
		contrib.Reallocate(Fb.GetNbDof());
		for (int i = 0; i < Fb.GetNbDof(); i++)
		  feval_phi(i).Reallocate(nb_quad_eprime);
		
		Real_wp  poidsm, poidsn;
		R_N point_n, point_m ;
		// int Troncature = 10 ;
		Complex_wp sum;
		Real_wp coef_sum = 2 * pi_wp;
		Real_wp invPerY = 1.0/PeriodeY;
		Real_wp invPerX = 1.0/PeriodeX;
		
		for( int n = 0; n < nb_quad_eprime; n++)
		  {
		    point_n = var_surf.PointsQuadrature(eprime, n);
		    // poidsn contient le ds
		    poidsn = var_surf.WeightsQuadrature(eprime, n);
		    
		    for( int m = 0; m < nb_quad_e; m++)
		      {
			point_m = var_surf.PointsQuadrature(e, m);
			// poidsm contient le ds
			poidsm = var_surf.WeightsQuadrature(e, m);
			
			sum = 0;
			if (point_n.GetM() == 2)
			  {
			    for(int l = - this->troncature ; l < this->troncature+1  ; l++)
			      sum +=  2 * pi_wp * abs(l) * exp(Iwp * (point_m(0)-point_n(0))
							       * 2.0 * pi_wp *Real_wp(l))    ;
			    
			    sum= sum * this->mu0 * poidsn * poidsm;
			  }
			else
			  {
			    for (int p = - this->troncature ; p < this->troncature+1  ; p++)
			      for ( int q = - this->troncature ; q < this->troncature+1  ; q++)
				{
				  
				  sum +=  coef_sum * sqrt(square(p*invPerX) + square(q*invPerY))
				    * exp(Iwp *  2.0 * pi_wp * ((point_m(0)-point_n(0))  *Real_wp(p)*invPerX
								+(point_m(1)-point_n(1)) *Real_wp(q)*invPerY));
				  
				}
			    
			    
			    
			    sum = sum * this->mu0 * poidsn * poidsm;
			    sum = sum/(PeriodeY*PeriodeX) ;
			  }
			
			feval(m) = sum;
		      }
		    
		    Fb.ComputeGaussIntegralSurfaceRef(feval, contrib, num_loc_e);
		    for (int i = 0; i < Fb.GetNbDof(); i++)
		      feval_phi(i)(n) = contrib(i);
		  }
		
		Complex_wp coef = nat_mat.GetCoefStiffness();
		Vector<Complex_wp> val_interac(Fb_prime.GetNbDof());
		IVect col_interac(Fb_prime.GetNbDof());
		for (int i = 0; i < Fb.GetNbDof(); i++)
		  if (Norm2(feval_phi(i)) > epsilon_machine)
		    {
		      Fb_prime.ComputeGaussIntegralSurfaceRef(feval_phi(i), contrib_prime, num_loc_eprime);
		      int c_interac = 0;
		      int irow = this->mesh_num.Element(iquad_e).GetNumberDof(i);
		      
		      for (int j = 0; j < Fb_prime.GetNbDof(); j++)
			if (abs(contrib_prime(j)) > epsilon_machine)
			  {
			    to_complex(coef*contrib_prime(j), val_interac(c_interac));
			    col_interac(c_interac) = offset_col + this->mesh_num
			      .Element(iquad_eprime).GetNumberDof(j);
			    c_interac++;
			  }
		      
		      if (c_interac > 0)
			mat_sp.AddInteractionRow(offset_row + irow, c_interac, col_interac, val_interac);
		    }
		
	      }
	  }
	cout << ' ' << endl ;
      }
    
  } // end computeImpedanceBoundary
  

  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique<TypeEquation>::
  ComputeRightHandSide(VectReal_wp & b_rhs, bool assemble)
  {
  }

  
  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique<TypeEquation>::
  ComputeRightHandSide(VectComplex_wp & b_rhs, bool assemble)
  {
    Vector<VectComplexe> vec_rhs;
    vec_rhs.SetData(1, &b_rhs);
    
    ComputeRightHandSide(vec_rhs, assemble);
    
    vec_rhs.Nullify();
  }
  
  
  template<class TypeEquation>
  void VarHelmholtz_CellulairePeriodique<TypeEquation>::
  ComputeRightHandSide(Vector<VectComplexe> & b_rhs, bool assemble)
  {
    const EllipticProblem<TypeEquation>& var_leaf = this->GetLeafClass();

    Vector<VirtualSourceFEM<Complexe, Dimension>* > f_vec(1);
    if (this->type_source_cellule == SRC_MU)
      {
        MuSource<Dimension> mu_source(var_leaf);
	f_vec(0) = &mu_source;
        
        cout << "ouhouh je suis dans le ComputeRightHandSide "<< endl ;
        var_leaf.ComputeGenericSource(b_rhs, f_vec, assemble);
      }
    else if (this->type_source_cellule == SRC_MU_X)
      {
	MuSource<Dimension> mu_source(var_leaf);
	f_vec(0) = &mu_source;
	
        cout << "ouhouh je suis dans le ComputeRightHandSide SRC_MU_X"<< endl ;
        mu_source.Source_x =true ;
        DISP(mu_source.Source_x) ;
        var_leaf.ComputeGenericSource(b_rhs, f_vec, assemble);
        
      }
    else if (this->type_source_cellule == SRC_MU_Y)
      {
        MuSource<Dimension> mu_source(var_leaf);
	f_vec(0) = &mu_source;
        
        cout << "ouhouh je suis dans le ComputeRightHandSide "<< endl ;
        mu_source.Source_x =false ;
        var_leaf.ComputeGenericSource(b_rhs, f_vec, assemble);
        
      }
    else if (this->type_source_cellule == SRC_ABSORBING)
      {
        SourceAbsorbante<Dimension> sourceAbsorbante(var_leaf);
        f_vec(0) = &sourceAbsorbante;
	
        cout << "ouhouh je suis dans le ComputeRightHandSide "<< endl ;
        if ( this->source_valeur_absolue)
            sourceAbsorbante.valeurAbsolue = true ;
        else
           sourceAbsorbante.valeurAbsolue = false;
            
        var_leaf.ComputeGenericSource(b_rhs, f_vec, assemble);
      }
  }  


  template<class Dimension>
  void EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelm(i, num_dof, mat_elem, nat_mat, *this, *this,
					this->GetReferenceElementH1(i)); 
  }

  
  /************
   * MuSource *
   ************/
  
  
  template<class Dimension> template<class TypeEquation>
  MuSource<Dimension>::MuSource(const EllipticProblem<TypeEquation>& var)
    : VirtualSourceFEM<Complex_wp, Dimension>(var), vars(var)
  {
  }
  
  
  template<class Dimension>
  bool MuSource<Dimension>::IsNonNullGradientSource(const VectR_N& s) { return true; }
  
  
  template<class Dimension>
  void MuSource<Dimension>::EvaluateGradientSource(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    //   cout << "ouhouh je suis dans EvaluateSourceGradient "<< endl ;
    f.Fill(0);
    // on recupere la matrice mu
    // int ref = this->vars.mesh.elements(i).GetReference();
    int ref = this->vars.mesh.Element(i).GetReference();
    //DISP(this->vars->R0) ;
    
    typename Dimension::MatrixN_Nsym_Complex_wp mu;
    mu = this->vars.ref_mu(ref).GetCoefficient(this->vars, i, j);
    
    // int  num_dof = this->vars.mesh.GetNumberDof(i,j) ;
    // cas 2d
    if (f.GetM() == 2)
      {
	f(0) = -mu(0,0)/ this->vars.R0 ; 
	f(1) = 0;
      }
    // Cas 3D cartesien
    else
      {
	if(Source_x == true)
	  {
	    f(0) = -mu(0,0) ;
	  }
	if(Source_x == false)
	  {
	    f(1) = -mu(0,0) ;                
	  }
      }
    
  }


  /********************
   * SourceAbsorbante *
   ********************/
  
  
  template<class Dimension> template<class TypeEquation>
  SourceAbsorbante<Dimension>
  ::SourceAbsorbante(const EllipticProblem<TypeEquation>& var)
    : VirtualSourceFEM<Complex_wp, Dimension>(var), vars(var)
  {
    ReferenceLineAbsorbing = BoundaryConditionEnum::LINE_ABSORBING; 
    valeurAbsolue =false ; 
  }
  
  
  template<class Dimension>
  bool SourceAbsorbante<Dimension>::IsNonNullSurfacicSource ( int ref )
  { 
    if (this->vars.mesh.GetBoundaryCondition(ref) == ReferenceLineAbsorbing)
      {
	return true;
      }
    else
      {
	return false;
      }
  }
  
  
  template<class Dimension>
  void SourceAbsorbante<Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
			   const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f)
  {
    // f.Reallocate(1);
    R_N X = PointsElem.GetPointQuadratureBoundary(k) ;
    
    if (X(Dimension::dim_N -1) > 0)
      f(0) =   1 ;
    if (X(Dimension::dim_N-1) < 0)
      {
	if (valeurAbsolue) 
	  f(0) = 1 ;
          else
            f(0) = -1 ;
      }    
  }
  
}//end montjoie

#define MONTJOIE_FILE_CELLULAIRE_PERIODIQUE_CXX
#endif
