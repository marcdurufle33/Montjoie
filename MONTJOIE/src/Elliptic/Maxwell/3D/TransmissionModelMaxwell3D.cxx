#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_MAXWELL3D_CXX

namespace Montjoie
{
  
  /*****************************
   * VarTransmission_Maxwell3D *
   *****************************/


  //! default constructor
  void VarTransmission_Maxwell3D::InitDefaultValues()
  {
    // parametres pour un couche homogene
    epsilon_layer = 1.0;
    sigma_layer = 0.0;
    mu_layer = 1.0;
    penalization_layer = 0.0;
    type_model_thin_layer = MODEL_PARAMETRIC;
    
    homogeneous_thin_layer = false;
    resistive_model = false;
    
    // parametres pour une couche periodique
    A = 0;
    B = 0;
    D1.Zero();
    D2.Zero();
    alpha = 0.5;
  }


  void VarTransmission_Maxwell3D::UpdateNumberOfDofs(int& nodl, int& nodl_scalar)
  {
    if (!homogeneous_thin_layer)
      {
        nodl += this->DdlVolMinus.GetM();
        nodl_scalar += this->DdlVolMinus.GetM();
        if (type_model_thin_layer == MODEL_CONDUCTIVE)
          {
            nodl += this->DdlVolMinus.GetM();
            nodl_scalar += this->DdlVolMinus.GetM();
          }
      }
  }


  void VarTransmission_Maxwell3D::UpdateGlobalDofs()
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int nb_proc = 1;
#endif

    MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    // updating MatchingDofOrig_Subdomain in order to add dofs for lambda
    if ((this->transmission_references.GetM() > 0) 
        && (!this->homogeneous_thin_layer) && (nb_proc > 1))
      {
        int nb_lambda = 1;
        if (type_model_thin_layer == MODEL_CONDUCTIVE)
          nb_lambda = 2;
        
        // we add global dof numbers of lambda on the array GlobDofNumber_Subdomain
        int Nsurf = this->DdlSurfLambda.GetM();
        int Nsurf_glob = this->GetNbVolumeDofOnSurface();
        mesh_num.GlobDofNumber_Subdomain.Resize(var_problem.GetNbDof());
        for (int i = 0; i < this->DdlSurfLambda.GetM(); i++)
          mesh_num.GlobDofNumber_Subdomain(mesh_num.GetNbDof() + i)
            = var_problem.GetNbGlobalMeshDof() + this->DdlSurfLambda(i);

        if (nb_lambda == 2)
          for (int i = 0; i < this->DdlSurfLambda.GetM(); i++)
            mesh_num.GlobDofNumber_Subdomain(mesh_num.GetNbDof() + Nsurf + i)
            = var_problem.GetNbGlobalMeshDof() + Nsurf_glob + this->DdlSurfLambda(i);            
        
        // we are checking if new processors are involved
        IVect proc_used(nb_proc);
        proc_used.Fill(-1);

        IVect& MatchingNumber_Pb = var_problem.GetProcMatchingNeighbor();
        Vector<IVect>& MatchingDofOrig_Pb = var_problem.GetOriginalMatchingDofNeighbor();
        for (int i = 0; i < MatchingNumber_Pb.GetM(); i++)
          proc_used(MatchingNumber_Pb(i)) = i;
        
        int nb_new = 0;
        for (int i = 0; i < this->MatchingNumber_Subdomain.GetM(); i++)
          if (proc_used(this->MatchingNumber_Subdomain(i)) == -1)
            {
              proc_used(this->MatchingNumber_Subdomain(i)) = nb_proc + i;
              nb_new++;
            }
	
        // resizing arrays if needed
        if (nb_new > 0)
          {
            nb_new += MatchingNumber_Pb.GetM();
            IVect OldMatchingNumber(MatchingNumber_Pb);
            Vector<IVect> OldMatchingOrig(MatchingDofOrig_Pb);
            
            MatchingNumber_Pb.Reallocate(nb_new);
            MatchingDofOrig_Pb.Reallocate(nb_new);
            nb_new = 0;
            for (int i = 0; i < proc_used.GetM(); i++)
              if (proc_used(i) != -1)
                {
                  if (proc_used(i) < nb_proc)
                    {
                      MatchingNumber_Pb(nb_new) = OldMatchingNumber(proc_used(i));
                      MatchingDofOrig_Pb(nb_new) = OldMatchingOrig(proc_used(i));
                    }
                  else
                    {
                      MatchingNumber_Pb(nb_new)
                        = this->MatchingNumber_Subdomain(proc_used(i)-nb_proc);
                      
                      MatchingDofOrig_Pb(nb_new).Clear();
                    }
                  
                  nb_new++;
                }
          }
        
        // Filling MatchingDofOrig with lambda dofs
        proc_used.Fill(-1);
        for (int i = 0; i < MatchingNumber_Pb.GetM(); i++)
          proc_used(MatchingNumber_Pb(i)) = i;
        
        int offset_glob = mesh_num.GetNbDof();
        for (int i = 0; i < this->MatchingNumber_Subdomain.GetM(); i++)
          {
            int k = proc_used(this->MatchingNumber_Subdomain(i));
            int offset = MatchingDofOrig_Pb(k).GetM();
            int nb_dof = offset + nb_lambda*this->MatchingDofOrig_Subdomain(i).GetM();
            MatchingDofOrig_Pb(k).Resize(nb_dof);
            for (int j = 0; j < this->MatchingDofOrig_Subdomain(i).GetM(); j++)
              MatchingDofOrig_Pb(k)(offset + j)
                = offset_glob + this->MatchingDofOrig_Subdomain(i)(j);
            
            int Ndof = this->MatchingDofOrig_Subdomain(i).GetM();
            if (nb_lambda == 2)
              for (int j = 0; j < this->MatchingDofOrig_Subdomain(i).GetM(); j++)
		MatchingDofOrig_Pb(k)(offset + Ndof + j)
		  = offset_glob + Nsurf + this->MatchingDofOrig_Subdomain(i)(j);
          }
      }
  }
  
  
  bool VarTransmission_Maxwell3D::IsSymmetricCondition() const
  {
    if (this->transmission_references.GetM() <= 0)
      return true;
    
    if (this->exact_solution)
      return true;
    
    if (this->Decalage*this->delta != 0)
      return false;
    
    return true;
  }
  
  
  //! reading a line of the data file
  void VarTransmission_Maxwell3D::
  SetInputData(const string& description_field, const VectString& parameters)
  {
    VarTransmission_Base<Dimension3>::SetInputData(description_field, parameters);
    
    if (!description_field.compare("HomogeneousThinLayer") )
      {
        // HomogeneousThinLayer corresponds to a test case with transmission conditions
        // to solve the continuous problem [E x n ] = 0, [H x n] = 0
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarTransmission_Maxwell" << endl;
	    cout << "HomogeneousThinLayer needs at least one parameter, for instance :" << endl;
	    cout << "HomogeneousThinLayer = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        homogeneous_thin_layer = false;
        if (!parameters(0).compare("YES"))
          homogeneous_thin_layer = true;
      }
    else if (!description_field.compare("ResistiveModel") )
      {
        // no longer used (AddThinLayerGitc)
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarTransmission_Maxwell" << endl;
	    cout << "ResistiveModel needs at least one parameter, for instance :" << endl;
	    cout << "ResistiveModel = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        resistive_model = false;
        if (!parameters(0).compare("YES"))
          resistive_model = true;
      }     
    else if (!description_field.compare("ReferenceTransmission"))
      {
        int pos = -1;
        int ref_max = 0;
        string keyword;
        for (int i = 0; i < parameters.GetM(); i++)
          if ( ( parameters(i) == "PARAM") || ( parameters(i) == "UNIFORM_SYM")
               || (parameters(i) == "UNIFORM_UNSYM") || (parameters(i) == "CHUN")
               || (parameters(i) == "RESISTIVE") || (parameters(i) == "RESISTIVE_SIGMA") 
               || (parameters(i) == "CONDUCTIVE") )
            {
              keyword = parameters(i);
              pos = i;
            }
         
        if (pos <= 0)
          {
            cout << "Expecting keyword PARAM, UNIFORM_SYM, UNIFORM_UNSYM,"
                 <<" RESISTIVE, RESISTIVE_SIGMA or CHUN" << endl;
            abort();
          }
         
        IVect ref_(pos);
        for (int i = 0; i < pos; i++)
          {
            ref_(i) = to_num<int>(parameters(i));
            ref_max = max (ref_max,ref_(i));
          }
         
        if (ref_max > var_problem.mesh.GetNbReferences())
          var_problem.mesh.ResizeNbReferences(ref_max);
         
	Vector<Vector<Complex_wp> >& Param_condition = var_boundary.GetParamCondition();
        if (ref_max > Param_condition.GetM())
          Param_condition.Resize(ref_max+1);
         
        if (parameters.GetM() > pos)
          {
            // we store these values in parameters_condition
            Vector<Complex_wp> parameters_condition(parameters.GetM() - pos-1);
            for (int i = pos+1; i < parameters.GetM(); i++)
              parameters_condition(i - pos-1) = to_num<Complex_wp>(parameters(i));
             
            if (homogeneous_thin_layer)
              {
                if (parameters_condition.GetM() <= 4)
                  {
                    cout << "In SetInputData of VarTransmission_Maxwell" << endl;
                    cout
                      << "ReferenceTransmission needs more parameters, for instance :" << endl;
                    cout
                      << "ReferenceTransmission = ref PARAM delta Decalage epsilon mu sigma"
                      << endl;
		     
                    cout << "Current parameters are : " << endl << parameters_condition << endl;
                    abort();
                  }

                // dans l'ordre delta, Decalage, epsilon, mu, sigma
                this->delta = real(parameters_condition(0));
                this->Decalage = real(parameters_condition(1)) ;
                epsilon_layer = parameters_condition(2);
                mu_layer = parameters_condition(3);
                sigma_layer = parameters_condition(4);                
                 
                if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
                  sigma_layer *= PhysicalConstant::impedance0;                 
              }
            else
              {
                if ( (keyword == "UNIFORM_SYM") || (keyword == "UNIFORM_UNSYM")
                     || (keyword == "RESISTIVE") || (keyword == "RESISTIVE_SIGMA")
                     || (keyword == "CHUN") || (keyword == "CONDUCTIVE") )
                  {
                    if (parameters_condition.GetM() <= 5)
                      {
                        cout << "In SetInputData of VarTransmission_Maxwell" << endl;
                        cout << "ReferenceTransmission needs more parameters, for instance :" 
                             << endl;
			 
                        cout << "ReferenceTransmission = ref UNIFORM_SYM "
                             << "delta Decalage epsilon mu sigma coef_penal" << endl;
                        cout << "Current parameters are : " << endl
                             << parameters_condition << endl;
			 
                        abort();
                      }

                    // modele de Chun : prendre Decalage = 0.5
                    // couche mince uniforme, on donne epsilon et mu de la couche
                    this->delta = real(parameters_condition(0)) ;
                    this->Decalage = real(parameters_condition(1)) ;
                    epsilon_layer = parameters_condition(2);
                    mu_layer = parameters_condition(3);
                    sigma_layer = parameters_condition(4);
                    penalization_layer = parameters_condition(5);
                     
                    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
                      sigma_layer *= PhysicalConstant::impedance0;
                     
                    if (keyword == "UNIFORM_SYM")
                      type_model_thin_layer = MODEL_DELOURME;
                    else if (keyword == "UNIFORM_UNSYM")
                      type_model_thin_layer = MODEL_PERON;
                    else if (keyword == "RESISTIVE")
                      type_model_thin_layer = MODEL_PERON_RESISTIVE;
                    else if (keyword == "RESISTIVE_SIGMA")
                      type_model_thin_layer = MODEL_PERON_RESISTIVE_SIGMA;
                    else if (keyword == "CHUN")
                      type_model_thin_layer = MODEL_CHUN;
                    else if (keyword == "CONDUCTIVE")
                      type_model_thin_layer = MODEL_CONDUCTIVE;
                  }
                else
                  {
                    if (parameters_condition.GetM() <= 7)
                      {
                        cout << "In SetInputData of VarTransmission_MaxwellAxi" << endl;
                        cout << "ReferenceTransmission needs more parameters, for instance :"
                             << endl;
			 
                        cout << "ReferenceTransmission = ref PARAM A B D1_11 D1_22 "
                             << "D2_11 D2_22 Decalage delta" << endl;
                        cout << "Current parameters are : " << endl
                             << parameters_condition << endl;
                        abort();
                      }

                    this->type_model_thin_layer = MODEL_PARAMETRIC;
                    this->Decalage = real(parameters_condition(6)) ;
                    this->delta = real(parameters_condition(7)) ;
                     
                    // les parametres sont dans l ordre suivant : 
                    // A,B,D1(0,0),D1(1,1),D2(0,0),D2(1,1),Decalage,delta
                    A =  parameters_condition(0) + 2 * this->Decalage ;
                    B =  parameters_condition(1) + 2 * this->Decalage ;
                     
                    D1(1,1) =  parameters_condition(2) + 2 * this->Decalage ;
                    D1(2,2) =  parameters_condition(3) + 2 * this->Decalage ;
                    D2(1,1) =  parameters_condition(4) + 2 * this->Decalage ;
                    D2(2,2) =  parameters_condition(5) + 2 * this->Decalage ;
                     
                    if ((abs(D1(1, 1)-D1(2, 2)) > 1e-10) || (abs(D2(1, 1)-D2(2, 2)) > 1e-10))
                      {
                        // cas orthotrope non gere
                        cout << "Anisotropic case not handled in transmission " << endl;
                        abort();
                      }
                    else
                      {
                        // cas isotrope
                        D1(0,0) = D1(1,1);  D2(0,0) = D2(1,1); 
                      }

                    DISP(this->Decalage); DISP(this->delta); DISP(parameters_condition);
                    DISP(D1) ; DISP(D2) ;DISP(A) ;DISP(B) ;
                  }                 
              }
          }
        else
          {
            abort();
          }
         
        for (int i = 0; i < pos; i++)
          var_problem.mesh.SetBoundaryCondition(ref_(i), BoundaryConditionEnum::LINE_TRANSMISSION);
         
        this->transmission_references = ref_;
      }
    else if (!description_field.compare("AlphaTransmission"))
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of VarTransmission_Maxwell" << endl;
            cout << "AlphaTransmission needs at least one parameter, for instance :" << endl;
            cout << "AlphaTransmission = 0.0" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
         
        alpha = to_num<Real_wp>(parameters(0));
      }
  }
  
  
  //! initialization of computations of transmission terms
  void VarTransmission_Maxwell3D::InitTransmission()
  {
    if (this->transmission_references.GetM() <= 0)
      return;
    
    if (this->exact_solution)
      return;

#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int nb_proc = 1, rank_proc(0);
#endif
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    
    if (!homogeneous_thin_layer)
      {
        Complex_wp eps0 = epsilon_layer;
        Complex_wp eps_minus = var_maxwell.ref_epsilon(this->RefOmegaMoins).GetConstant()(0, 0);
        Complex_wp mu_minus = var_maxwell.ref_mu(this->RefOmegaMoins).GetConstant()(0, 0);
        Complex_wp sigma_minus = var_maxwell.ref_sigma(this->RefOmegaMoins).GetConstant()(0, 0);

        Complex_wp mu0 = mu_layer;
        Complex_wp eps_plus = var_maxwell.ref_epsilon(this->RefOmegaPlus).GetConstant()(0, 0);
        Complex_wp mu_plus = var_maxwell.ref_mu(this->RefOmegaPlus).GetConstant()(0, 0);
        Complex_wp sigma_plus = var_maxwell.ref_sigma(this->RefOmegaPlus).GetConstant()(0, 0);
        Complex_wp alpha = this->Decalage;
        
        switch (type_model_thin_layer)
          {
          case MODEL_DELOURME:
            eps0 += Iwp*sigma_layer/var_problem.GetOmega();
            eps_minus += Iwp*sigma_minus/var_problem.GetOmega();
            eps_plus += Iwp*sigma_plus/var_problem.GetOmega();
            
            A = 1.0/eps0 - 0.5/eps_minus - 0.5/eps_plus + 2.0*alpha;
            B = 1.0/mu0 - 0.5/mu_minus - 0.5/mu_plus + 2.0*alpha;
            D1.SetIdentity();
            D2.SetIdentity();
            D1 *= mu0 - 0.5*(mu_minus + mu_plus) + 2.0*alpha;
            D2 *= eps0 - 0.5*(eps_minus + eps_plus) + 2.0*alpha;        
            DISP(A); DISP(B); DISP(D1); DISP(D2);
            DISP(eps0); DISP(eps_minus); DISP(eps_plus);
            DISP(mu0); DISP(mu_minus); DISP(mu_plus);
            break;
          case MODEL_PERON:
            eps0 += Iwp*sigma_layer/var_problem.GetOmega();
            eps_minus += Iwp*sigma_minus/var_problem.GetOmega();
            eps_plus += Iwp*sigma_plus/var_problem.GetOmega();
            
            A = 1.0/eps0 - 1.0/eps_plus;
            B = 1.0/mu0 - 1.0/mu_plus;
            D1.SetIdentity();
            D2.SetIdentity();
            D1 *= mu0 - mu_plus;
            D2 *= eps0 - eps_plus;
            break;
          case MODEL_PERON_RESISTIVE :
            eps0 += Iwp*sigma_layer/var_problem.GetOmega();
            eps_minus += Iwp*sigma_minus/var_problem.GetOmega();
            eps_plus += Iwp*sigma_plus/var_problem.GetOmega();
            
            A = 1.0/eps0;
            B = 0.0;
            D1.Zero();
            D2.Zero();
            break;
          case MODEL_PERON_RESISTIVE_SIGMA :
            eps_minus += Iwp*sigma_minus/var_problem.GetOmega();
            eps_plus += Iwp*sigma_plus/var_problem.GetOmega();
            
            A = 1.0/eps0 - 1.0/eps_plus;
            B = 1.0/mu0 - 1.0/mu_plus;
            D1.SetIdentity();
            D2.SetIdentity();
            D1 *= mu0 - mu_plus;
            D2 *= eps0 - eps_plus;
            break;
          case MODEL_CHUN :
            eps0 += Iwp*sigma_layer/var_problem.GetOmega();
            
            A = 1.0/eps0;
            B = 1.0/mu0;
            D1.SetIdentity();
            D2.SetIdentity();
            D1 *= mu0;
            D2 *= eps0;
            break;
          case MODEL_CONDUCTIVE :
            {
              eps_minus += Iwp*sigma_minus/var_problem.GetOmega();
              eps_plus += Iwp*sigma_plus/var_problem.GetOmega();
              Complex_wp sigma_tilde = sigma_layer * this->delta*this->delta;
              Complex_wp gamma = exp(3.0*Iwp*pi_wp/4)*sqrt(var_problem.GetOmega()*mu0*sigma_tilde);
              A1 = -0.5/var_problem.GetSquareOmega()*(1.0/eps_minus + 1.0/eps_plus);
              B1 = 2.0*mu0/gamma*tanh(gamma/2) - 0.5*(mu_plus+mu_minus);
              A2 = 0.25*A1; 
              B2 = mu0/(2.0*gamma*tanh(gamma/2)) - 0.125*(mu_plus+mu_minus);
              A3 = -0.25/var_problem.GetSquareOmega()*(1.0/eps_plus - 1.0/eps_minus);
              B3 = -0.25*(mu_plus - mu_minus);
              //DISP(A1); DISP(B1); DISP(A2); DISP(B2); DISP(A3); DISP(B3);
            }
            break;
          }
        
        // rajout de la penalisation
        for (int k = 0; k < 3; k++)
          D1(k, k) += this->delta*penalization_layer;
        
        // since lambda is the trace of H on the boundary
        // Neumann condition on E becomes Dirichlet condition on H
        // we retrieve these Dirichlet dofs
        // DofOnNeumann(i) = true if the dof is located on a Neumann boundary
        Vector<int> DofOnNeumann(var_problem.GetNbDof()); DofOnNeumann.Fill(0);
        for (int i = 0; i < var_problem.mesh.GetNbBoundaryRef(); i++)
          {
            int ref = var_problem.mesh.BoundaryRef(i).GetReference();
            if (var_problem.mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEUMANN)
              {
                int num_face = i;
                int ne = var_problem.mesh.Boundary(num_face).numElement(0);
                if (ne >= 0)
                  {
                    int num_loc = var_problem.mesh.Element(ne).GetPositionBoundary(num_face);
		    const ElementReference<Dimension3, 2>& Fb = var_problem.GetReferenceElementHcurl(ne);

		    for (int j = 0; j < Fb.GetNbDofBoundary(num_loc); j++)
		      DofOnNeumann(mesh_num.Element(ne).
				   GetNumberDof(Fb.GetLocalNumber(num_loc, j))) = 1;
		  }
              }
          }
	
#ifdef SELDON_WITH_MPI
        var_problem.ReduceDistributedVector(DofOnNeumann, MPI_MAX, 1);
#endif

        // counting the number of Dirichlet dofs for lambda
        int nb_dir = var_boundary.GetNbDirichletDof();
        Vector<bool> is_lambda_dir(this->DdlVolMinus.GetM());
        is_lambda_dir.Fill(false);
        if (nb_proc > 1)
          for (int i = 0; i < this->DdlVolMinus.GetM(); i++)
            {
              if (this->ProcVolMinus(i) == rank_proc)
                if (DofOnNeumann(this->DdlVolMinus(i)) == 1)
                  is_lambda_dir(i) = true;
              
              if (this->ProcVolPlus(i) == rank_proc)
                if (DofOnNeumann(this->DdlVolPlus(i)) == 1)
                  is_lambda_dir(i) = true;
              
              if (is_lambda_dir(i))
                nb_dir++;
            }
        else
          for (int i = 0; i < this->DdlVolMinus.GetM(); i++)
            if ((DofOnNeumann(this->DdlVolMinus(i)) == 1)
                || (DofOnNeumann(this->DdlVolPlus(i)) == 1))
              {
                is_lambda_dir(i) = true;
                nb_dir++;
                if (type_model_thin_layer == MODEL_CONDUCTIVE)
                  nb_dir++;
              }
        
        int offset = mesh_num.GetNbDof();
        if (nb_dir > var_boundary.GetNbDirichletDof())
          {
            int offset2 = offset + this->DdlVolPlus.GetM();
            // updating Dirichlet dofs
            IVect Dir_dof = var_boundary.GetDirichletDofNumber();
            Dir_dof.Resize(nb_dir);
            nb_dir = var_boundary.GetNbDirichletDof();
            for (int i = 0; i < this->DdlVolMinus.GetM(); i++)
              if (is_lambda_dir(i))
                {
                  Dir_dof(nb_dir) = offset + i;
                  nb_dir++;
                  if (type_model_thin_layer == MODEL_CONDUCTIVE)
                    {
                      Dir_dof(nb_dir) = offset2 + i;
                      nb_dir++;
                    }
                }
            
            var_boundary.SetDirichletDofs(nb_dir, Dir_dof);
          }
      }
    
    var_boundary.NewColumnNumbers_Impedance.Reallocate(var_problem.GetNbDof()) ;
    var_boundary.NewRowNumbers_Impedance.Reallocate(var_problem.GetNbDof());
    var_boundary.NewRowNumbers_Impedance.Fill(-1) ;
    var_boundary.NewColumnNumbers_Impedance.Fill(-1) ;

    var_boundary.ProcColumnNumbers_Impedance.Reallocate(var_problem.GetNbDof()) ;
    var_boundary.ProcRowNumbers_Impedance.Reallocate(var_problem.GetNbDof());
    var_boundary.ProcRowNumbers_Impedance.Fill(rank_proc) ;
    var_boundary.ProcColumnNumbers_Impedance.Fill(rank_proc) ;
    
    //DISP(A3); DISP(B3); DISP(A2); DISP(B2); DISP(A1); DISP(B1);
    
    {
      //A1 = 0.0; B1 = 1.0; A2 = 0; B2 = 0; A3 = 0; B3 = 0;
      /*
	this->DdlVolMinus.WriteText("num_minus.dat");
	this->DdlVolPlus.WriteText("num_plus.dat");
	IVect InverseDof(var.GetNbDof());
	int Nsurf = this->DdlVolMinus.GetM();
	InverseDof.Fill(-1);
	for (int i = 0; i < this->DdlVolMinus.GetM(); i++)
        InverseDof(this->DdlVolMinus(i)) = i;
      
	SetPoints<Dimension> PointsElem;
	SetMatrices<Dimension> MatricesElem;
	VectR3 s; R3 pt, pt_loc, normale; Complexe val;
	Matrix3_3 dfj, dfjm1; Real_wp poids; TinyVector<Complexe, 3> vec_u, vec_v;
	Matrix<Complexe, General, ArrayRowSparse> MatNTimesPhi(Nsurf, Nsurf), MatMassPhi(Nsurf, Nsurf), MatCurlPhi(Nsurf, Nsurf);
      
	for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
        {
	int ref = var.mesh.BoundaryRef(i).GetReference();
	if (ref == this->RefGammaMoins)
	{
	int iquad = var.mesh.BoundaryRef(i).numElement(0);
	int r = var.mesh_num.GetOrderElement(iquad);
	const ElementReference<Dimension, 2>& Fb = var.GetReferenceElement(iquad);
              
	var.mesh.GetVerticesElement(iquad, s);
	Fb.FjElem(s, PointsElem, var.mesh, iquad);
	Fb.DFjElem(s, PointsElem, MatricesElem, var.mesh, iquad);
              
	int num_face = i;
	int num_loc = var.mesh.Element(iquad).GetPositionBoundary(num_face);
	int nb_points_quadrature = Fb.GetNbQuadBoundary(num_loc);
              
	int nb_dof_face = Fb.GetNbDofBoundary(num_loc);
              
	Fb.FjSurfaceElem(s, PointsElem, var.mesh, iquad, num_loc);
              
	Fb.DFjSurfaceElem(s, PointsElem, MatricesElem,
	var.mesh, iquad, num_loc);
              
	VectR3 phi(Fb.GetNbDof()), curl_phi(Fb.GetNbDof());
	for (int k = 0; k < nb_points_quadrature; k++)
	{
	pt_loc = Fb.PointsND(Fb.GetQuadNumber(num_loc, k));
	pt = PointsElem.GetPointQuadratureBoundary(k);
	normale = MatricesElem.GetNormaleQuadratureBoundary(k);
	dfj = MatricesElem.GetPointQuadratureBoundary(k);
	poids = MatricesElem.GetDsQuadratureBoundary(k)*Fb.WeightsQuadratureBoundary(k, num_loc);
	GetInverse(dfj, dfjm1);
                  
	Fb.ComputeValuesPhi(pt_loc, phi, dfjm1, var.mesh_num, iquad);

	Fb.ComputeValuesGradientPhi(pt_loc, curl_phi, dfjm1, var.mesh_num, iquad);
                  
	//DISP(k); DISP(pt_loc); DISP(phi); DISP(curl_phi); DISP(poids); DISP(normale);
	for (int iloc = 0; iloc < nb_dof_face; iloc++)
	for (int jloc = 0; jloc < nb_dof_face; jloc++)
	{
	int i = Fb.GetLocalNumber(num_loc, iloc);
	int dofi = var.mesh_num.Element(iquad).GetNumberDof(i);
	dofi = InverseDof(dofi);

	int j = Fb.GetLocalNumber(num_loc, jloc);
	int dofj = var.mesh_num.Element(iquad).GetNumberDof(j);
	dofj = InverseDof(dofj);
                        
	TimesProd(normale, phi(j), vec_u);
	val = poids*DotProd(vec_u, phi(i));
	MatNTimesPhi.AddInteraction(dofi, dofj, val);
                        
	TimesProd(normale, phi(i), vec_v);
	val = poids*DotProd(vec_u, vec_v);
	MatMassPhi.AddInteraction(dofi, dofj, val);
                        
	val = poids*DotProd(normale, curl_phi(i))*DotProd(normale, curl_phi(j));
	MatCurlPhi.AddInteraction(dofi, dofj, val);
	}
	}
	}
          
	MatNTimesPhi.WriteText("MatNTimesPhi.dat");
	MatMassPhi.WriteText("MatMass.dat");
	MatCurlPhi.WriteText("MatStiff.dat");
        }
      */
    }
  }
  

  void VarTransmission_Maxwell3D
  ::AddTransmissionTerms(const Real_wp& alpha,
			 const GlobalGenericMatrix<Real_wp>& nat_mat,
			 VirtualMatrix<Real_wp>& mat_sp, int, int)
  {
    cout << "Not implemented" << endl;
    abort();
  }

  
  //! adding terms of variational formulation due to transmission conditions
  void VarTransmission_Maxwell3D
  ::AddTransmissionTerms(const Complex_wp& alpha,
			 const GlobalGenericMatrix<Complex_wp>& nat_mat,
			 VirtualMatrix<Complex_wp>& mat_sp, int offset_row, int offset_col)
  {
    if (transmission_references.GetM() <= 0)
      return;
    
    IVect ref_plus(var_problem.mesh.GetNbReferences()+1);
    IVect ref_minus(var_problem.mesh.GetNbReferences()+1);
    ref_plus.Fill(0); ref_minus.Fill(0);
    ref_plus(this->RefGammaPlus) = 1; ref_minus(this->RefGammaMoins) = 1;

    if (homogeneous_thin_layer)
      {
        this->SetModifiedRowNumbers(0);
        this->SetModifiedColNumbers(0);
        
	// here the integrals are computed over Gamma+,
	// therefore the normale is opposite to normale n
        
	// we add term alpha \int n/mu \times rot E^e \cdot U^e
	// + alpha \int  n/mu \times rot U^e \cdot E^e
	fimpedance.coef_val = -this->alpha; fimpedance.coef_grad = -this->alpha;
	//fimpedance.coef_u = var.alpha_penalization;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, false, var_problem);
	
	// we add term -alpha \int n/mu \times rot E^e \cdot U^c
	fimpedance.coef_val = 0.0; fimpedance.coef_grad = this->alpha;
	//fimpedance.coef_u = -var.alpha_penalization;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, true, var_problem);
	
	// we add term -alpha \int n/mu \times rot U^e \cdot E^c
	fimpedance.coef_val = this->alpha; fimpedance.coef_grad = 0.0;
	//fimpedance.coef_u = 0.0;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	
	
	// integrals over Gamma-
	
	// we add term (1-alpha) \int n/mu \times rot E^c \cdot U^c
	// and (1-alpha) \int  n/mu \times rot U^c \cdot E^c
	fimpedance.coef_val = -(1.0-this->alpha); fimpedance.coef_grad = -(1.0-this->alpha);
	//fimpedance.coef_u = var.alpha_penalization;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, false, var_problem);
	
	// we add term -(1-alpha) \int n/mu \times rot E^c \cdot U^e
	fimpedance.coef_val = 0.0; fimpedance.coef_grad = (1.0-this->alpha);
	//fimpedance.coef_u = -var.alpha_penalization;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, true, var_problem);
	
	// we add term  -(1-alpha) \int n/mu \times rot U^c \cdot E^e
	fimpedance.coef_val = (1.0-this->alpha); fimpedance.coef_grad = 0.0;
	//fimpedance.coef_u = 0.0;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);                        
	
        return;        
      }

    if (type_model_thin_layer == MODEL_CONDUCTIVE)
      {
        this->SetModifiedRowNumbers(0);
        
	// integrals over Gamma+
	
	// terme : - \int_Gamma+ (n \times \lambda_1) phi^+ 
	fimpedance.C1 = -1.0;
	fimpedance.C2.Zero();
	fimpedance.C3 = 0 ;
        
	this->SetModifiedColNumbers(1);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	
	// terme : - 1/2 \int_Gamma+ (n \times \lambda_2) phi^+
	fimpedance.C1 = -0.5;
	fimpedance.C2.Zero();
	fimpedance.C3 = 0 ;
        
	this->SetModifiedColNumbers(2);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	
	// terme : -\int_\Gamma+ E^+ \times n \cdot \mu_1
	fimpedance.C1 = 1.0;
	fimpedance.C2.Zero();
	fimpedance.C3 = 0 ;
        
	this->SetModifiedRowNumbers(1);
	this->SetModifiedColNumbers(0);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, true, var_problem);
	
	// terme : +\int_\Gamma+ E^- \times n \cdot \mu_1
	fimpedance.C1 = -1.0;
	fimpedance.C2.Zero();
	fimpedance.C3 = 0 ;
        
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, true, var_problem);
	
	// terme : -1/2 \int_\Gamma+ E^+ \times n \cdot \mu_2
	fimpedance.C1 = 0.5;
	fimpedance.C2.Zero() ;
	fimpedance.C3 = 0 ;
        
	this->SetModifiedRowNumbers(2);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, true, var_problem);
	
	// terme : -1/2 \int_\Gamma+ E^- \times n \cdot \mu_2
	fimpedance.C1 = 0.5;
	fimpedance.C2.Zero() ;
	fimpedance.C3 = 0 ;
        
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
					       offset_row, offset_col, fimpedance, true, true, var_problem);
	
	// terme : delta \int_\Gamma+  A1 curl_\Gamma \lambda_1 curl_\Gamma \mu_1
	fimpedance.C1 = 0;
	fimpedance.C2.Zero();
	fimpedance.C3 = this->delta*A1;
        
	this->SetModifiedRowNumbers(1);
	this->SetModifiedColNumbers(1);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, true, var_problem);
	
	// terme : -delta \int_Gamma+ B1 \lambda_1 \cdot \mu_1
	fimpedance.C1 = 0;
	fimpedance.C2.SetDiagonal(-this->delta*B1);
	fimpedance.C3 = 0;
        
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, true, var_problem);
	
	// terme : delta \int_\Gamma+ A2 curl_\Gamma \lambda_2 curl_\Gamma mu_2
	fimpedance.C1 = 0;
	fimpedance.C2.Zero();
	fimpedance.C3 = this->delta*A2;
        
	this->SetModifiedRowNumbers(2);
	this->SetModifiedColNumbers(2);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, true, var_problem);
            
	// terme : -delta \int_\Gamma+ B2 \lambda_2 \cdot \mu_2
	fimpedance.C1 = 0;
	fimpedance.C2.SetDiagonal(-this->delta*B2);
	fimpedance.C3 = 0;
        
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, true, var_problem);
	
	// terme : delta \int_\Gamma+ A3 curl_\Gamma \lambda_1 curl_\Gamma mu_2
	fimpedance.C1 = 0;
	fimpedance.C2.Zero();
	fimpedance.C3 = this->delta*A3;
        
	this->SetModifiedRowNumbers(2);
	this->SetModifiedColNumbers(1);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, true, var_problem);
	
	// terme : -delta \int_\Gamma+ B3 \lambda_1 \cdot \mu_2
	fimpedance.C1 = 0;
	fimpedance.C2.SetDiagonal(-this->delta*B3);
	fimpedance.C3 = 0;
        
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, true, var_problem);

	// terme : delta \int_\Gamma+ A3 curl_\Gamma \lambda_2 curl_\Gamma mu_1
	fimpedance.C1 = 0;
	fimpedance.C2.Zero();
	fimpedance.C3 = this->delta*A3;
	
	this->SetModifiedRowNumbers(1);
	this->SetModifiedColNumbers(2);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, true, var_problem);
	
	// terme : -delta \int_\Gamma+ B3 \lambda_2 \cdot \mu_1
	fimpedance.C1 = 0;
	fimpedance.C2.SetDiagonal(-this->delta*B3);
	fimpedance.C3 = 0;
        
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, true, var_problem);
	
	// integrals over Gamma-
	
	// terme : \int_Gamma+ (n times lambda_1) phi^-
	fimpedance.C1 = -1.0;
	fimpedance.C2.Zero() ;
	fimpedance.C3 = 0 ;
        
	this->SetModifiedColNumbers(1);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	
	// terme : - 1/2 \int_Gamma+ (n \times \lambda_2) phi^-
	fimpedance.C1 = 0.5;
	fimpedance.C2.Zero() ;
	fimpedance.C3 = 0 ;
        
	this->SetModifiedColNumbers(2);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	
        return;
      }
    
    // coefficient 1/2 C delta
    Real_wp coef_curve = 0.5*this->delta;
    coef_curve = 0;
    
    // Integral over Gamma+
    
    // terme :  - \int_gamma+ (n times E^+) mu
    fimpedance.C1 = -1.0+coef_curve;
    fimpedance.C2.Zero() ;
    fimpedance.C3 = 0 ;
    
    this->SetModifiedRowNumbers(1);           
    var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
					    offset_row, offset_col, fimpedance, false, true, var_problem);
    //exit(0);

    // terme : \int_{gamma+} (n times E^-) mu
    fimpedance.C1 = 1.0+coef_curve;
    fimpedance.C2.Zero() ;
    fimpedance.C3 = 0 ;
    
    this->SetModifiedRowNumbers(1);
    this->SetModifiedColNumbers(0);
    var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
					    offset_row, offset_col, fimpedance, true, true, var_problem);
    
    // terme : \int_Gamma+ (n times lambda) phi^+ 
    fimpedance.C1 = 1.0-coef_curve;
    fimpedance.C2.Zero() ;
    fimpedance.C3 = 0 ;
    
    this->SetModifiedColNumbers(1);
    var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
					    offset_row, offset_col, fimpedance, true, false, var_problem);
    
    // termes : this->delta A /omega^2  int_Gamma+ curl lambda curl mu
    // - \delta int D1 lambda mu
    fimpedance.C1 = 0.0; //DISP(D1); DISP(this->delta);
    fimpedance.C2 = -this->delta * D1 ;     
    fimpedance.C3 = this->delta * A /var_problem.GetSquareOmega();
    //fimpedance.C3 = 0.0;
    
    this->SetModifiedRowNumbers(1);
    this->SetModifiedColNumbers(1);        
    var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
					    offset_row, offset_col, fimpedance, true, true, var_problem);
    
    // Termes : delta B int_Gamma+ curl E^+ curl phi^+
    // - int_Gamma+ omega^2 D2 delta E^+ phi^+
    fimpedance.C1 = 0 ;
    fimpedance.C2 = -(square(1.0-this->alpha) * this->delta *var_problem.GetSquareOmega()) *   D2  ;
    fimpedance.C3 = this->delta *square(1.0- this->alpha) *  B ;
    
    var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
					    offset_row, offset_col, fimpedance, false, false, var_problem);
    
    // termes : -delta B int_Gamma- curl E^- curl phi+  
    // + int_Gamma- omega^2 D2 delta E^- phi+ (changement de lignes avec -)
    fimpedance.C1 = 0 ;
    fimpedance.C2 = -(this->delta  * this->alpha*(1.0-this->alpha)
		      * var_problem.GetSquareOmega()) *   D2  ;
    fimpedance.C3 = this->alpha*(1.0-this->alpha)* this->delta * B ;
    
    this->SetModifiedColNumbers(0);
    var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
					    offset_row, offset_col, fimpedance, true, false, var_problem);
        
    // Integrals over Gamma-
    
    // terme 6 : int_{Gamma-} (n \times lambda} phi^- 
    fimpedance.C1 = 1.0+coef_curve;
    fimpedance.C2.Zero() ;
    fimpedance.C3 = 0 ;
    
    this->SetModifiedColNumbers(1);
    var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
					    offset_row, offset_col, fimpedance, true, false, var_problem);
    
    // termes 10 et 14  : delta B int_Gamma- curl E^- curl phi^-
    // - int_Gamma- omega^2 D2 delta E^- phi^-
    fimpedance.C1 = 0 ;
    fimpedance.C2 = (-this->alpha*this->alpha *this->delta *var_problem.GetSquareOmega()) *   D2  ;
    fimpedance.C3 = this->alpha*this->alpha * this->delta * B ;
        
    var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
					    offset_row, offset_col, fimpedance, false, false, var_problem);
      
    // Termes : delta B int_Gamma+ curl E^+ curl phi- 
    // - int_Gamma+ omega^2 D2 delta E^+ phi-
    fimpedance.C1 = 0 ;
    fimpedance.C2 = -(this->delta *var_problem.GetSquareOmega()* this->alpha*(1.0-this->alpha)) *D2;
    fimpedance.C3 = this->delta* this->alpha*(1.0-this->alpha)  * B ;
    
    this->SetModifiedColNumbers(0);
    var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
					    offset_row, offset_col, fimpedance, true, false, var_problem);
  }


  /**********************************
   * ImpedanceTransmissionMaxwell3D *
   **********************************/
  

  //! computation of impedance
  void ImpedanceTransmissionMaxwell3D::
  EvaluateImpedancePhi(int i, int iquad, int num_edge, int num_loc, int k,
		       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
		       const SetPoints<Dimension3>& Pts, const SetMatrices<Dimension3>& Mat)
  {
    if (k == 0)
      feval.Reallocate(Mat.GetNbPointsQuadratureBoundary());
    
    if (var_transmission.homogeneous_thin_layer)
      {
        feval(k).Zero();

        R3 normale = Mat.GetNormaleQuadratureBoundary(k);
        TinyMatrix<Complex_wp, Symmetric, 3, 3> ProjTangent;
        GetTangentialProjector(normale, ProjTangent);
        feval(k) = Real_wp(-0.5)*ProjTangent;
      }
    else
      {
        R3 normale = Mat.GetNormaleQuadratureBoundary(k);
        
        TinyMatrix<Complex_wp, General, 3, 3> ProjTimesN;
        ProjTimesN(1,0) = normale(2) ; 
        ProjTimesN(2,0) = -normale(1) ;
        ProjTimesN(0,1) = -normale(2) ; 
        ProjTimesN(2,1) = normale(0) ;
        ProjTimesN(0,2) = normale(1) ;
        ProjTimesN(1,2) = -normale(0) ;
        
        feval(k) = -dot(ProjTimesN, dot(C2, ProjTimesN)) + C1 * ProjTimesN;
      }    
  }
  
  
  //! computation of impedance for curl
  void ImpedanceTransmissionMaxwell3D::
  EvaluateImpedanceGrad(int i, int iquad, int num_edge, int num_loc, int k,
			const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_domain,
			const SetPoints<Dimension3>& Pts, const SetMatrices<Dimension3>& Mat)
  {
    if (k == 0)
      feval_grad.Reallocate(Pts.GetNbPointsQuadratureBoundary());

    R3 normale = Mat.GetNormaleQuadratureBoundary(k);
    if (var_transmission.homogeneous_thin_layer)
      {
        if (k == 0)
          invMu_normale.Reallocate(Pts.GetNbPointsQuadratureBoundary());
        
        int refD = var_transmission.RefOmegaPlus;
        Complex_wp mu = var_maxwell.ref_mu(refD)(0, 0); mu = 1.0;
        invMu_normale(k) = normale; Mlt(1.0/mu, invMu_normale(k));
        feval_grad(k).Zero();
      }
    else
      {
        TinyMatrix<Complex_wp, Symmetric, 3, 3> ProjNormale;
        GetNormalProjector(normale, ProjNormale);
        
        feval_grad(k) = C3 * ProjNormale;    
      }
  }
  

  //! applying impedance to phi
  void ImpedanceTransmissionMaxwell3D
  ::ApplyImpedancePhi_Hcurl(int n, int j, int offset, const TinyVector<Real_wp, 3>& phi_boundary,
                            const TinyVector<Real_wp, 3>& dphi, VectComplex_wp& feval_phi)
  {
    R3_Complex_wp fvec, fphi;
    if (var_transmission.homogeneous_thin_layer)
      {
        TimesProd(dphi, invMu_normale(j), fvec);
        Mlt(coef_val, fvec);
        
        Mlt(feval(j), phi_boundary, fphi);
        Add(coef_u, fphi, fvec);
      }
    else
      Mlt(feval(j), phi_boundary, fvec) ;    
    
    CopyVector(fvec, offset, feval_phi);
  }
  

  //! applying impedance to phi  
  void ImpedanceTransmissionMaxwell3D
  ::ApplyImpedanceCurl(int n, int j, int offset, const TinyVector<Real_wp, 3>& phi,
		       const TinyVector<Real_wp, 3>& curl_phi, VectComplex_wp& g_phi)
  {
    R3_Complex_wp fvec;
    if (var_transmission.homogeneous_thin_layer)
      {
        TimesProd(invMu_normale(j), phi, fvec);
        Mlt(coef_grad, fvec);
       }
    else
      Mlt(feval_grad(j), curl_phi, fvec) ;
    
    CopyVector(fvec, offset, g_phi);
  }      
  
} // end namespace

#define MONTJOIE_FILE_TRANSMISSION_MODEL_MAXWELL3D_CXX
#endif

  
