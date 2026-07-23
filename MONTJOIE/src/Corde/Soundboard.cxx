#ifndef MONTJOIE_FILE_SOUNDBOARD_CXX

#include "Oscillator.cxx"

namespace Montjoie
{
  
  //! Default constructor
  SourceChi::SourceChi()
  {      
    x0.Init(0.3,0.3);
    r0 = 0.1;
    alpha =-log(1e-6);
    beta = alpha / (pi_wp*r0*r0);
    ampl = 1.0;
  }
  
  
  //! computation of the source 
  void SourceChi::EvaluateFunction(const R2& x,  Vector<Real_wp>& f) const
  {
    // the source is oriented along u_z and of amplitude
    // ampl beta exp( - alpha (r/r0)^2)
    // beta is here chosen such that the integral of f is equal to 1 (if ampl = 1)
    Real_wp r = x.Distance(x0);
    f(0) = 0;
    f(1) = 0;
    f(2) = ampl*beta*exp(-alpha*square(r/r0) );
    //DISP(x0); DISP(r); DISP(beta); DISP(ampl); DISP(alpha); DISP(r0); DISP(f);
  }
  
  
  //! modification of (x0, y0), r0 and ampl with the data file
  void SourceChi::SetInputData(const string& description_field, const VectString& parameters,
                               ParameterMultistring& param_corde)
  {
    if (!description_field.compare("MecanicInitialCondition"))
      {
        // if the first parameter is a note (C2, D3, etc)
        // we use the datas contained in param_corde
        int nb = 3;
        if (parameters.GetM() < 3)
          {
  	    cout << "In SetInputData of SourceChi" << endl;
  	    cout << "Not enough parameters, possible formats are :" << endl;
  	    cout << "MecanicInitialCondition = Note r0" << endl;
  	    cout << "MecanicInitialCondition = Note r0 amplitude" << endl;
  	    cout << "MecanicInitialCondition = x0 y0 r0" << endl;
  	    cout << "MecanicInitialCondition = x0 y0 r0 amplitude" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
          }
	
        if (isalpha(parameters(1)[0]))
          {
            param_corde.GetParamAttache(parameters(1), x0(0), x0(1));
            nb = 2;
          }
        else
          x0.Init(to_num<Real_wp>(parameters(1)), to_num<Real_wp>(parameters(2)));
	
        if (parameters.GetM() <= nb)
  	  {
  	    cout << "In SetInputData of SourceChi" << endl;
  	    cout << "Radius r0 is missing :" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
  	  }
	
        r0 = to_num<Real_wp>(parameters(nb++));
        beta = alpha / (pi_wp*r0*r0);
    	if(parameters.GetM() > nb)
    	  ampl = to_num<Real_wp>(parameters(nb++));
    	else
    	  ampl = 1.0;
	
    	//DISP(x0);
    	//DISP(r0);
        //DISP(alpha);
    	//DISP(ampl);
      }
  }
  
  
  /**************
   * Soundboard *
   **************/
 

  //! constructeur par defaut
  Soundboard::Soundboard()
  {    
    nb_eigenvalues_per_run = 0;
    size_buffer_mecanic_sismo = 100;
    save_eigenvectors_asked = false;
    load_eigenvectors_asked = false;
    path_eigenvectors = "./";
    
#ifdef SELDON_WITH_SCALAPACK
    size_block_schur = 64;
#endif
  }
 
  
  //! set parameters read in the data file
  /*!
    \param[in] description_field keyword of the line of the data file
    \param[in] parameters parameters associated with the keyword
   */
  void Soundboard::SetInputData(const string & description_field, const VectString & parameters)
  {
    Oscillator_Base::SetInputData(description_field, parameters);
    if(!description_field.compare("ModalDamping"))
      {
        if (parameters.GetM() <= 2)
  	  {
  	    cout << "In SetInputData of Soundboard" << endl;
  	    cout << "ModalDamping needs three parameters, for instance :" << endl;
  	    cout << "ModalDamping = a b c" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
  	  }
        
        for(int i = 0; i < 3; i++)
          VectParamAmo(i) = to_num<Real_wp>(parameters(i));
      }
    else if (description_field == "NbEigenvaluesToComputePerRun")
      {
        if (parameters.GetM() <= 0)
  	  {
  	    cout << "In SetInputData of Soundboard" << endl;
  	    cout << "NbEigenvaluesToComputePerRun needs one parameter, for instance :" << endl;
  	    cout << "NbEigenvaluesToComputePerRun = N" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
  	  }

      nb_eigenvalues_per_run = to_num<int>(parameters(0));
      }
    else if (description_field == "SaveEquivalentImpedance")
      {
        if (parameters.GetM() <= 0)
  	  {
  	    cout << "In SetInputData of Soundboard" << endl;
  	    cout << "SaveEquivalentImpedance needs four parameters, for instance :" << endl;
  	    cout << "SaveEquivalentImpedance = file1 file2 file3 file4" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
  	  }
	
      write_oscillator_param_bool = true;
      mass_file_name = parameters(0);
      lambda_file_name = parameters(1);
      amo_file_name = parameters(2);
      vect_second_membre_file_name = parameters(3);
      }
    else if (description_field == "SaveEigenvectors")
      {
        if (parameters.GetM() <= 0)
  	  {
  	    cout << "In SetInputData of Soundboard" << endl;
  	    cout << "SaveEigenvectors needs one parameter, for instance :" << endl;
  	    cout << "SaveEigenvectors = YES" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
  	  }
	
        if (parameters(0) == "YES")
          save_eigenvectors_asked = true;
        else
          save_eigenvectors_asked = false;
      }
    else if (description_field == "LoadEigenvectors")
      {
        if (parameters.GetM() <= 0)
  	  {
  	    cout << "In SetInputData of Soundboard" << endl;
  	    cout << "LoadEigenvectors needs one parameter, for instance :" << endl;
  	    cout << "LoadEigenvectors = YES" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
  	  }
        
        if (parameters(0) == "YES")
          {
            load_eigenvectors_asked = true;
          }
        else
          load_eigenvectors_asked = false;
      }
    else if (description_field == "PathEigenvectors")
      {
        if (parameters.GetM() <= 0)
  	  {
  	    cout << "In SetInputData of Soundboard" << endl;
  	    cout << "PathEigenvectors needs one parameter, for instance :" << endl;
  	    cout << "PathEigenvectors = path" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
  	  }
        
        path_eigenvectors = parameters(0);
      }
#ifdef SELDON_WITH_SCALAPACK
    else if (description_field == "SizeSchurBlock")
      {
        size_block_schur = to_num<int>(parameters(0));
      }
#endif
  }
  
  
  //! Avance d'un pas de temps
  void Soundboard::ActualiseLambda(Real_wp coef, Real_wp pulse)
  {
    chrono.Start(VirtualTimer::COMM);
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = VarSoundboard.comm_group_mode;
    
    // broadcasting coef to all the processors
    Vector<int64_t> xtmp;
    MpiBcast(comm, &coef, xtmp, 1, 0);
#endif
    
    chrono.Stop(VirtualTimer::COMM);
    

    // avance d'un pas
    VectReal_wp Lambdaold(Lambda);
    VectReal_wp LambdaPointold(LambdaPoint);
    Lambda.Fill(Real_wp(0));
    LambdaPoint.Fill(Real_wp(0));
	
    for (int g = 0; g < nev; g++)
      {
        Lambda(g) = MatriceResolExacte0(g) * Lambdaold(g)
          + MatriceResolExacte1(g) * LambdaPointold(g) + coef*Lambda0(g);
        
        LambdaPoint(g) = MatriceResolExacte0Point(g) * Lambdaold(g)
          + MatriceResolExacte1Point(g) * LambdaPointold(g) + coef*Lambda0Point(g);
      }
    
    if (RhsMecanicModal.GetM() == nev)
      {
        for (int g = 0; g < nev; g++)
          {
            Lambda(g) += RhsMecanicModal(g)*pulse*MatriceResolExacteSecondMembre(g);
            LambdaPoint(g) += RhsMecanicModal(g)*pulse*MatriceResolExacteSecondMembrePoint(g);
          }
      }
  }
  
  
  //! Construction of the 2-D mesh and computation of eigenmodes
  /*!
    \param[in] DOSSIER directory where outputs are written
    \param[in] source_coupling_string distribution function of the bridge on the soundboard
    \param[in] lines_input_file data file stored in a vector of strings
    \param[in] remove_file if true, the output files are removed before new values are appended
   */
  void Soundboard::RunSoundboardDiagonalization(string DOSSIER, SourceChi& source_coupling_string,
                                                Vector<string>& lines_input_file, bool remove_file)
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = VarSoundboard.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif

    VarSoundboard.InitIndices(PhysicalConstant::nb_max_indices);
    VarSoundboard.SetTypeEquation("REISSNER_MINDLIN");
    
    // lecture fichier de donnees
    ReadInputFile(lines_input_file, VarSoundboard);
    ReadInputFile(lines_input_file, *this);
    
    if (rank_proc == 0)
      {
        output_lambda.  Init( DOSSIER + "SismoLambda.bin", size_buffer_lambda, remove_file);
        output_lambda.SetBinary(true);
        output_lambda.SetDoublePrecision(false);
      }
    
    // calcul maillage et element fini
    VarSoundboard.ComputeMeshAndFiniteElement(string("TRIANGLE_LOBATTO"), true);
    if(nb_proc == 1)
      VarSoundboard.mesh.Write("Soundboard.mesh");
    
    // autres initialisations
    VarSoundboard.PerformOtherInitializations();
    VarSoundboard.ComputeMassMatrix();
    PerformDiagonalization(VarSoundboard, path_eigenvectors, path_eigenvectors,
                           lines_input_file, save_eigenvectors_asked, load_eigenvectors_asked,
                           source_coupling_string);
    
  }
  
  
  //! computation of eigenmodes
  /*!
    \param[in] var_soundboard object describing the soundboard (mesh + materials)
    \param[in] prefix prefixe utilise pour les reprises (chargement)
    \param[in] prefix_save prefixe utilise pour les reprises (sauvegarde)
    \param[in] lines_data_file data file stored in a vector of strings
    \param[in] save_reprise si vrai, les donnees sont sauvegardees
    pour assurer une reprise ulterieure
    \param[in] load_reprise si vrai, les donnees sont chargees d'une reprise precedente
    \param[in] source_coupling_string fonction de distribution du chevalet sur la table d'harmonie
    \param[in] display_eigen si vrai, les vecteurs propres sont interpoles et ecrits sur le disque
  */
  void Soundboard::PerformDiagonalization(EllipticProblem<ReissnerMindlinEquation> & var_soundboard,
                                          const string& prefix, const string& prefix_save,
                                          const Vector<string>& lines_data_file,
                                          bool save_reprise, bool load_reprise,
                                          SourceChi& source_coupling_string, bool display_eigen)
  {
#ifdef SELDON_WITH_MPI
    chrono.SetCommunicator(var_soundboard.comm_group_mode);
    int nb_proc; MPI_Comm_size(var_soundboard.comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(var_soundboard.comm_group_mode, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif
    
    // le solver 
    Dimension2 dim; Symmetric property;
    All_LinearSolver glob_solver(var_soundboard);
    ReadInputFile(lines_data_file,  glob_solver);
    
    //var_soundboard.print_level = 3;
    //glob_solver.SetPrintLevel(3);
    const MeshNumbering<Dimension2>& mesh_num = var_soundboard.GetMeshNumbering(0);
    glob_solver.SelectOptimalLinearSolver(mesh_num.GetOrder(),
					  mesh_num.GetNbDof(), dim, property);
    
    // objet declare le pb aux valeurs propres
    EigenProblemMontjoie<Real_wp> eigen_solver(var_soundboard, glob_solver);
    
    nFEM = mesh_num.GetNbDof();
    
    ReadInputFile(lines_data_file, eigen_solver);
    nev = eigen_solver.GetNbAskedEigenvalues();
    if (nb_eigenvalues_per_run > 0)
      eigen_solver.SetNbAskedEigenvalues(nb_eigenvalues_per_run);
    
    var_soundboard.print_level = 3;
    eigen_solver.SetPrintLevel(2);
    
    string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
    if (load_reprise)
      {
        eigen_values.Read(prefix+"EigenvalPianoSave"+suffix);
        eigen_vectors.Read(prefix+"EigenvecPianoSave"+suffix);
        if (nev < eigen_values.GetM())
          {
            eigen_values.Resize(nev);
            eigen_vectors.Resize(eigen_vectors.GetM(), nev);
          }
        else
          {
            if (nev > eigen_values.GetM())
              {
                if (rank_proc == 0)
                  cout << "Continuing computations with " << eigen_values.GetM()
                       << " eigenvalues" << endl;
              }            
            nev = eigen_values.GetM();
          }
      }
    else
      {
        if (nev > 0)
          {
            var_soundboard.SetStorageFiniteElementMatrix(var_soundboard.MATRIX_STORED);
            
            FemMatrixFreeClass_Base<Real_wp>* Kh;
            FemMatrixFreeClass_Base<Real_wp>* Mh;
            
            Mh = var_soundboard.GetNewIterativeMatrix(Real_wp(0));
            Kh = var_soundboard.GetNewIterativeMatrix(Real_wp(0));
            
            eigen_solver.InitMatrix(*Kh, *Mh);
            
            if (eigen_solver.GetNbAskedEigenvalues() != nev)
              {
                // on alloue les tableaux resultats
                eigen_values.Reallocate(nev);
                eigen_vectors.Reallocate(eigen_solver.GetM(), nev);
                
                // variables temporaires
                glob_chrono.Reset(VirtualTimer::PML);
                glob_chrono.Start(VirtualTimer::PML);
                VectReal_wp lambda_r, lambda_i;
                Matrix<Real_wp, General, ColMajor> eigen_vec;
                Real_wp omega_cur = 0;
                Real_wp shift = eigen_solver.GetShiftValue(), shift_imag = 0;
                int n = 0; Real_wp coef = 1.0+1e-8;
                while (n < nev)
                  {
                    // on calcule les valeurs propres pour le shift choisi
                    var_soundboard.UpdateShiftAdimensionalization(shift, shift_imag);
                    eigen_solver.SetShiftValue(shift); eigen_solver.SetImagShiftValue(shift_imag);
                    
                    GetEigenvaluesEigenvectors(eigen_solver, lambda_r, lambda_i, eigen_vec);
                    
                    var_soundboard.UpdateEigenvaluesAdimensionalization(lambda_r,
                                                                        lambda_i, eigen_vec);
                    
                    if (rank_proc == 0)
                      cout << "Temporary eigenvalues :" << endl << lambda_r << endl << endl;
                    
                    // on stocke toutes les nouvelles valeurs propres
                    for (int i = 0; i < lambda_r.GetM()-4; i++)
                      {
                        if ((lambda_r(i) > coef*omega_cur) && (n < nev))
                          {
                            eigen_values(n) = lambda_r(i);
                            for (int j = 0; j < eigen_vec.GetM(); j++)
                              eigen_vectors(j, n) = eigen_vec(j, i);
                            
                            omega_cur = lambda_r(i);
                            n++;
                          }
                      }
                    
                    // on passe aux valeurs propres suivantes
                    shift = omega_cur + 0.35*abs(omega_cur - lambda_r(0));
                  }
                
                glob_chrono.Stop(VirtualTimer::PML);
                if (rank_proc == 0)
                  cout << "Time to compute eigenvalues and eigenvectors = " 
                       << glob_chrono.GetSeconds(VirtualTimer::PML) << endl;
              }
            else
              {
                Real_wp shift = eigen_solver.GetShiftValue(),
                  shift_imag = eigen_solver.GetImagShiftValue();
                
                var_soundboard.UpdateShiftAdimensionalization(shift, shift_imag);
                eigen_solver.SetShiftValue(shift); eigen_solver.SetImagShiftValue(shift_imag);
                
                // on remplit eigen_values, lambda_imag, eigen_vectors;
                // renvoie lambda et u tels que K u = lambda M u
                GetEigenvaluesEigenvectors(eigen_solver, eigen_values, lambda_imag, eigen_vectors);
                
                var_soundboard.UpdateEigenvaluesAdimensionalization(eigen_values, lambda_imag,
                                                                    eigen_vectors);
              }
            
            delete Kh; delete Mh;
            
            if (rank_proc == 0)
              {
                DISP(eigen_solver.GetStoppingCriterion());
                DISP(eigen_values);
                // writing eigenvalue on the disk
                eigen_values.WriteText(eigen_solver.GetEigenvalFileName());
              }
            
            //if (nb_proc == 2)
            //for (int i = 0; i < eigen_vectors.GetM(); i++)
            //  eigen_vectors(i, 1) *= -1.0;
            
            Real_wp Lmax = eigen_values(eigen_values.GetM() - 1);
            
            bool check_ortho = false;
            if (check_ortho)
              {        
                // checking eigenvalues and eigenvectors
                int N = var_soundboard.GetNbDof();
                MatrixDirect Mh(N, N), Kh(N, N);
                VectReal_wp Dh(N);
                FemMatrixFreeClass<Real_wp, ReissnerMindlinEquation> free_mat(var_soundboard);
                GlobalGenericMatrix<Real_wp> nat_mat;
                nat_mat.SetCoefMass(1);
                nat_mat.SetCoefStiffness(0);
                var_soundboard.AddMatrixWithBC(free_mat, Mh, nat_mat);
                
                nat_mat.SetCoefMass(0);
                nat_mat.SetCoefStiffness(1);
                var_soundboard.AddMatrixWithBC(free_mat, Kh, nat_mat);
                //Mh.WriteText("Mh.dat");
                //Kh.WriteText("Kh.dat");
                
                for (int i = 0; i < N; i++)
                  Dh(i) = Mh(i, i);
                
                VectReal_wp x(N), y(N), Mx(N);
                Real_wp threshold = 1e-10;
                for (int i = 0; i < eigen_values.GetM(); i++)
                  {
                    for (int j = 0; j < N; j++)
                      x(j) = eigen_vectors(j, i);
                    
                    Mlt(Kh, x, y);
                    for (int j = 0; j < N; j++)
                      y(j) /= Dh(j);
                    
                    Real_wp vmax = 0;
                    Real_wp Li = 0;
                    for (int j = 0; j < N; j++)
                      {
                        if (abs(x(j)) > vmax)
                          {
                            vmax = abs(x(j));
                            Li = y(j)/x(j);
                          }
                      }
                    
                    if (abs(Li - eigen_values(i))/Lmax > 10000*threshold)
                      {
                        cout << "not a good eigenvalue, error = "
                             << abs(Li - eigen_values(i))/Lmax << endl;
                        //abort();
                      }
                    
                    Mlt(1.0/Li, y);
                    if (Li > threshold*Lmax)
                      {
                        bool eigen_vec_ok = true;
                        Real_wp err_max = 0;
                        for (int j = 0; j < N; j++)
                          {
                            if (abs(x(j) - y(j)) > 10000*threshold)
                              {
                                eigen_vec_ok = false;
                                err_max = max(err_max, abs(x(j) - y(j)));
                              }
                            
                            if (!eigen_vec_ok)
                              {
                                cout << "not a good eigenvector, error = " << err_max << endl;
                                //abort();                    
                              }
                          }
                      }
                  }
                
                // checking orthogonality of eigenvectors
                for (int i = 0; i < eigen_values.GetM(); i++)
                  {
                    for (int i2 = 0; i2 < eigen_values.GetM(); i2++)
                      {
                        for (int j = 0; j < N; j++)
                          {
                            x(j) = eigen_vectors(j, i);
                            y(j) = eigen_vectors(j, i2);
                          }
                        
                        Mlt(Mh, x, Mx);
                        Real_wp scal_prod = DotProd(Mx, y);
                        if (i == i2)
                          {
                            cout << "Norme au carre de x = " << scal_prod << endl;
                          }
                        else
                          {
                            if (abs(scal_prod) > threshold)
                              {
                                cout << "eigenvectors not orthogonal" << endl;
                                abort();
                              }
                          }
                      }
                  }
              }
            
            VectReal_wp eigen_mode(3*nFEM);
            Matrix<Real_wp, General, ColMajor> eigen_vec2(eigen_vectors);
            eigen_vectors.Reallocate(3*nFEM, nev);
            Vector<Vector<Real_wp> > eigen_mode_v;
            eigen_mode_v.SetData(1, &eigen_mode);
            const IVect& local_col_numbers = eigen_solver.GetLocalColumnNumbers();
            for (int j = 0; j < nev; j++)
              {
                eigen_mode.Fill(0);
                // original eigenvectors are computed
                for (int i = 0; i < eigen_vec2.GetM(); i++)
                  eigen_mode(local_col_numbers(i)) = eigen_vec2(i, j);
                
                if (nb_proc > 1)
                  var_soundboard.AddDomains(eigen_mode);
                
                //display_eigen = true;
                if (display_eigen)
                  {
                    if ((j == 0) && (nb_proc == 1))
                      eigen_mode.Write("eigenmode.dat");
                    
                    string numero = NumberToString(j);
                    for (int k = 0; k < var_soundboard.output_mesh_param.GetM(); k++)
                      var_soundboard.output_mesh_param(k).
                        SetFileName(2, GetBaseString(var_soundboard.output_mesh_param(k)
                                                     .GetTotalFieldFile())
                                    + "." + numero.substr(numero.size()-3,3) + string(".vtk"));
                    
                    for (int k = 0; k < var_soundboard.output_grid_param.GetM(); k++)
                      var_soundboard.output_grid_param(k).
                        SetFileName(2, GetBaseString(var_soundboard.output_grid_param(k)
                                                     .GetTotalFieldFile())
                                    + numero + string(".dat"));
                    
                    var_soundboard.WriteOutputFile(eigen_mode_v, 2);
                  }
                
                for (int i = 0; i < 3*nFEM; i++)
                  eigen_vectors(i, j) = eigen_mode(i);                
              }
            
            eigen_mode_v.Nullify();
            
          }
        if (save_reprise)
          {
            eigen_values.Write(prefix_save+"EigenvalPianoSave"+suffix);
            eigen_vectors.Write(prefix_save+"EigenvecPianoSave"+suffix);
          }
      }
    
    // on remplit les vecteurs contenant les coefs pour l'EDO
    vect_lambda.Reallocate(nev);
    vect_masse.Reallocate(nev);
    vect_amortissement.Reallocate(nev);
    for (int g = 0; g < nev; g++)
      {
        vect_lambda(g) = eigen_values(g);
        vect_masse(g) = 1.0;
        Real_wp freqtemp = sqrt(vect_lambda(g))/(2.0*pi_wp);
        vect_amortissement(g) = VectParamAmo(0)*freqtemp*freqtemp
          + VectParamAmo(1)*freqtemp+VectParamAmo(2);
      }


    
    // pour faire une boucle sur les elts
    int Ntot = mesh_num.GetNbDof(); // nb de ddl du maillage;
    
    // matrice de masse sans materiaux (sans rho delta^3/12 et rho delta)
    VectReal_wp DiagMasse(3*mesh_num.GetNbDof());
    DiagMasse.Fill(0);
    for (int i = 0; i < var_soundboard.mesh.GetNbElt(); i++)
      {
        int nb_dof = mesh_num.Element(i).GetNbDof();
        if (var_soundboard.Glob_jacobian(i).GetM() > 1)
          {
            for (int j = 0; j < nb_dof; j++)
              DiagMasse(mesh_num.Element(i).GetNumberDof(j))
                += var_soundboard.Glob_jacobian(i)(j);
          }
        else
          {
            const VectReal_wp& weights = var_soundboard.WeightsND(i);
            Real_wp jacob = var_soundboard.Glob_jacobian(i)(0);
            for (int j = 0; j < nb_dof; j++)
              DiagMasse(mesh_num.Element(i).GetNumberDof(j)) += jacob*weights(j);
          }
      }
    
    var_soundboard.AddDomains(DiagMasse);

    Vector<VectReal_wp> chi(1); //un vecteur de vecteurs avec un seul vecteur (bon...)    
    chi(0).Reallocate(Ntot*ReissnerMindlinEquation::nb_unknowns);
    
    Vector<VectString> param_s;
    VolumetricSource<ReissnerMindlinEquation> f(var_soundboard, param_s);
    f.SetVolumeSourceFunction(source_coupling_string);
    
    //DISP(nFEM); DiagMasse.WriteText("Dh.dat");
    /* pour visualiser les points d'attache
      ParameterMultistring param_corde;
    Vector<string> param(4);
    for (int num = 0; num < 88; num++)
    //for (int num = 0; num < 1; num++)
      {        
        param(0) = "USER";
        param(1) = param_corde.GetNoteString(num);
        param(2) = "0.1";
        param(3) = "1.0";
        //DISP(param(1));
        
        //f.fsrc = source_coupling_string;
        f.fsrc.SetInputData(string("MecanicInitialCondition"), param,
                            param_corde);
        
        chi(0).Fill(0);
        var_soundboard.AddVolumetricSource(chi, f);
        
        //chi(0).WriteText("MhSol.dat");
        for (int i = 0; i < nFEM; i++)
          chi(0)(2*nFEM+i) /= DiagMasse(i);
        
        string numero = to_str(num);
        for (int k = 0; k < var_soundboard.output_mesh_param.GetM(); k++)
          var_soundboard.output_mesh_param(k).
            SetFileName(2, GetBaseString(var_soundboard.output_mesh_param(k).GetTotalFieldFile())
                        + "." + numero + string(".bb"));
        
        for (int k = 0; k < var_soundboard.output_grid_param.GetM(); k++)
          var_soundboard.output_grid_param(k).
            SetFileName(2, GetBaseString(var_soundboard.output_grid_param(k).GetTotalFieldFile())
                        + numero + string(".dat"));
        
        var_soundboard.WriteOutputFile(chi, 2);
      }
    */

    Vector<VirtualSourceFEM<Real_wp, Dimension2>* > f_vec(1);
    f_vec(0) = &f;
    chi(0).Fill(0);
    var_soundboard.AddVolumetricSource(Real_wp(1), chi, f_vec);
    
    // VectSecondMembre = J_h
    VectSecondMembre.Reallocate(Ntot*ReissnerMindlinEquation::nb_unknowns);
    for (int g = 0; g < Ntot*ReissnerMindlinEquation::nb_unknowns; g++) 
      VectSecondMembre(g) = chi(0)(g);

    f.NullifyVolumeSourceFunction();
    
    // calcul du second membre mecanique
    var_soundboard.ComputeRightHandSide(chi, false);
    RhsMecanic = chi(0);
    //chi(0).Write("rhs.dat");
    //DISP(Norm2(RhsMecanic));
    
    // matrice de masse avec materiaux
    VectReal_wp DiagMasseMateriau(3*mesh_num.GetNbDof());
    DiagMasseMateriau.Fill(0);
    var_soundboard.GetMassMatrix(DiagMasseMateriau);
    
    // on veut obtenir M f avec M la matrice de masse avec materiaux
    for (int i = 0; i < nFEM; i++)
      {
        Real_wp invMass = 1.0/DiagMasse(i);
        RhsMecanic(i) *= invMass*DiagMasseMateriau(i);
        RhsMecanic(nFEM + i) *= invMass*DiagMasseMateriau(i);
        RhsMecanic(2*nFEM + i) *= invMass*DiagMasseMateriau(i);
      }
    
    // actu de Pu
    Pu.Reallocate(nFEM, nev);
    for (int j = 0; j < nev; j++)
      {
        for (int i = 0; i < nFEM; i++)
          Pu(i, j) = eigen_vectors(2*nFEM + i, j);
      }
    
    // calcul de psi_j(x_0) pour chaque point d'un MecanicSismoPoint
    int nb_points_sismo = 0;    
    VectR2 Points;
    for (int i = 0; i < var_soundboard.var_grid.GetM(); i++)
      if (var_soundboard.var_grid(i).GetOutputType() == var_soundboard.var_grid(i).POINT)
        Points.PushBack(R2(var_soundboard.var_grid(i).GetXmin(),
                           var_soundboard.var_grid(i).GetYmin())); 
    
    if (Points.GetM() > 0)
      {
        nb_points_sismo = Points.GetM();
        ValLambdaSismoPoint.Reallocate(nb_points_sismo, nev);
        ValLambdaSismoPoint.Fill(0);
        
        int nb = 0; VectReal_wp val_psi(nev);
        Vector<TinyVector<Real_wp, 1> > phi;
        for (int i = 0; i < var_soundboard.var_grid.GetM(); i++)
          if (var_soundboard.var_grid(i).GetOutputType() == var_soundboard.var_grid(i).POINT)
            {
              int npoint = -1;
              if (var_soundboard.var_grid(i).GetNbPoints() > 0)
                npoint = var_soundboard.var_grid(i).GetPointNumber(0);
	      
              val_psi.Fill(0);
              if (npoint >= 0)
                {
                  int num_elem = var_soundboard.all_points_display.GetElementNumber(npoint);
                  R2 pt_loc = var_soundboard.all_points_display.GetLocalCoordinate(npoint);
                  if (num_elem >= 0)
                    {
                      var_soundboard.GetReferenceElementH1(num_elem).
                        ComputeValuesPhi(pt_loc, phi, var_soundboard.all_points_display.GetDFjm1(npoint),
                                         mesh_num, num_elem);
                      
                      for (int j = 0; j < phi.GetM(); j++)
                        {
                          int num_dof = mesh_num.Element(num_elem).GetNumberDof(j);
                          for (int k = 0; k < nev; k++)
                            val_psi(k) += phi(j)(0)*eigen_vectors(2*nFEM+num_dof, k);
                        }
                    }
                }
              
              AssembleModalVector(val_psi);
              for (int k = 0; k < nev; k++)
                ValLambdaSismoPoint(nb, k) = val_psi(k);
              
              nb++;
            }
        
      }
    
    if (rank_proc == 0)
      cout << "sismo OK" << endl;
    
    //exit(0);
    //abort();
  }
  
  
  //! sums the modal vector among processors involved in the computation
  void Soundboard::AssembleModalVector(VectReal_wp& x)
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = VarSoundboard.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc > 1)
      {
        VectReal_wp rhs(x); Vector<int64_t> xtmp;
        MpiAllreduce(comm, rhs, xtmp, x, rhs.GetM(), MPI_SUM);	
      }
#endif
  }
  
  
  //! Construction of modal right hand sides
  void Soundboard::ConstruitVectSecondMembreModal()
  {
    // avoir déjà les vecteurs propres
    // VectSecondMembreModal = 	tr(P_h^U) J_h
    VectSecondMembreModal.Reallocate(nev); VectSecondMembreModal.Fill(Real_wp(0));
    Mlt(SeldonTrans, eigen_vectors, VectSecondMembre, VectSecondMembreModal);
    
    this->AssembleModalVector(VectSecondMembreModal);

    if (write_oscillator_param_bool)
      {
        //        vect_lambda.WriteText("RaideursDiago.txt");
        //        vect_masse.WriteText("MassesDiago.txt");
        //        vect_amortissement.WriteText("DissipDiago.txt");
        //        VectSecondMembreModal.WriteText("VectSecondMembreModalDiago.txt");
        
        vect_lambda.WriteText(lambda_file_name);
        vect_masse.WriteText(mass_file_name);
        vect_amortissement.WriteText(amo_file_name);
        VectSecondMembreModal.WriteText(vect_second_membre_file_name);
      }

    RhsMecanicModal.Reallocate(nev); RhsMecanicModal.Fill(Real_wp(0));
    Mlt(SeldonTrans, eigen_vectors, RhsMecanic, RhsMecanicModal);
    
    this->AssembleModalVector(RhsMecanicModal);
    //DISP(RhsMecanicModal);
    // VectSecondMembreModal.Write("vectPar"
    // +to_str(rank_proc)+"N"+to_str(nb_proc)+".dat");	
  }
  
  
  //! Resolution du complement de Schur et calcul de scal_b
  Real_wp Soundboard::ConstruitLambdaSchur()
  {
    LambdaVolLambda.Reallocate(nev);  LambdaVolLambda.Fill(0);
    LambdaJ.Reallocate(nev);
    Copy(Lambda0, LambdaJ);
    SolveSchur(LambdaJ);
    // LambdaJ= S_p^{-1} FJ = S_p^{-1} Lambda0
    Real_wp scal_b = 0.0;
    for (int g = 0; g < nev; g++)
      scal_b += LambdaJ(g)*VectSecondMembreModal(g);     
    return(scal_b);
  }

  
  //! construction of needed arrays for time iterations
  Real_wp Soundboard::ConstructAll()
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = VarSoundboard.comm_group_mode;
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int rank_proc(0);
#endif

    ConstruitVectSecondMembreModal();
    Real_wp scal_b = ConstruitMatricesResolExacte();
    ConstruitLambda0etPoint();
    if(deltat_output_lambda!= 0){
      
      if(rank_proc == 0)
      {
        VectReal_wp freq_propres(nev), alpha(nev);
        for(int g = 0; g < nev; g++)
        {
          Real_wp invmasse =1.0/vect_masse(g);
          freq_propres(g) = sqrt(vect_lambda(g)*invmasse)/(2.0*pi_wp);
          alpha(g) = vect_amortissement(g)*invmasse;
        }
        output_lambda.AddVect(freq_propres);
        output_lambda.AddVect(alpha);
      }
    }
    return(scal_b);    
  }

  //! initial condition
  void Soundboard::Initialise(Real_wp rhodelta)
  {
    // construction de la donnee ini sur les FEM
    // pas besoin de multiplier par M, on a F tel que Mtilde^-1 F est le truc sur les EF
    // il suffit donc d'ajuster le rho*delta delta parce que la projection sur U^0 
    // utilise la mass matrice FEM et non pas physique
		
    // on utilise VectSecondMembre, on le multiplie par rho*delta tr(P)
    // passage en base modale
    //		Real_wp rho = VarSoundboard.ref_rho(1);
    Mlt(SeldonTrans, eigen_vectors, VectSecondMembre, Lambda);
    Mlt(rhodelta, Lambda);
    
    this->AssembleModalVector(Lambda);
  }
  
  
  //! avance d'un pas (calcul de Lambda^{n+1/2})
  void Soundboard::ActualiseLambda(Real_wp coef, VectReal_wp & SecondMembreP)
  {
    chrono.Start(VirtualTimer::COMM);
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = VarSoundboard.comm_group_mode;
    
    // broadcasting coef to all the processors
    Vector<int64_t> xtmp;
    MpiBcast(comm, &coef, xtmp, 1, 0);
#endif

    chrono.Stop(VirtualTimer::COMM);


    VectReal_wp Lambdaold(Lambda);
    VectReal_wp LambdaPointold(LambdaPoint);
    Lambda.Fill(Real_wp(0));
    LambdaPoint.Fill(Real_wp(0));
    for (int g = 0; g < nev; g++)
    {
      Lambda(g) = MatriceResolExacte0(g) * Lambdaold(g)
          + MatriceResolExacte1(g) * LambdaPointold(g)
          + MatriceResolExacteSecondMembre(g)*SecondMembreP(g)+ coef*Lambda0(g);
        
      LambdaPoint(g) = MatriceResolExacte0Point(g) * Lambdaold(g)
          + MatriceResolExacte1Point(g) * LambdaPointold(g)
	        + MatriceResolExacteSecondMembrePoint(g)*SecondMembreP(g)+ coef*Lambda0Point(g);        
    }
  }
  
  
  //! On met a jour Lambda^{n+1/2} apres calcul de P^{n+1}
  void Soundboard::ActualiseLambdaSchur(Real_wp coef, VectReal_wp& rhs, VectReal_wp& Fvol)
  {

    chrono.Start(VirtualTimer::COMM);
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = VarSoundboard.comm_group_mode;

    // broadcasting coef to all the processors
    Vector<int64_t> xtmp;
    MpiBcast(comm, &coef, xtmp, 1, 0);
#endif

    chrono.Stop(VirtualTimer::COMM);
    chrono.Start(VirtualTimer::SCHEME);
    
    VectReal_wp Lambdaold(Lambda);
    VectReal_wp LambdaPointold(LambdaPoint);
    VectReal_wp LambdaTmp(nev);
    
    // Lambda^n+1/2 = LambdaTilde^n+1/2 + F Lambda^J
    for (int g = 0; g < nev; g++)
      Lambda(g) = LambdaVolLambda(g)+ coef*LambdaJ(g);
    
    //on calcule 1/2 P^T Ch Mh^-1 Ch^T P Lambda^n+1/2 = (rhs + F Lambda0 - Lambda^{n+1/2}) / R_{dt}
    for (int g = 0; g < nev; g++)
      LambdaTmp(g) = (rhs(g) + coef*Lambda0(g) - Lambda(g)) / MatriceResolExacteSecondMembre(g);
    
    // on met a jour Fvol pour contenir P^T Ch (P^{n+1/2} + P^{n-1/2})/2
    for (int g = 0; g < nev; g++)
      Fvol(g) -= LambdaTmp(g);
    
    // on calcule LambdaPoint^n+1/2
    for (int g = 0; g < nev; g++)
    {
      LambdaPoint(g) = MatriceResolExacte0Point(g) * Lambdaold(g)
          + MatriceResolExacte1Point(g) * LambdaPointold(g)
	        + MatriceResolExacteSecondMembrePoint(g)*Fvol(g) + coef*Lambda0Point(g);        
    }           

    chrono.Stop(VirtualTimer::SCHEME);
  }
  
  
  //! Calcul de Lambda^{n+1/2}
  template<class VibroPb>
  void Soundboard::CalculeLambda(const Real_wp& coef, const Real_wp& pulse,
                                 const VectReal_wp& Gvol_minus,
				 const VectReal_wp& Gvol_plus, const VectReal_wp& Mh,
                                 VibroPb& var_vibro, VectReal_wp &LambdaNew)
  {		    
    VectReal_wp Lambdaold(Lambda);
    VectReal_wp LambdaPointold(LambdaPoint);

    // calcul de Fvol = R_{∆t} tr(Ph) Cbarre Gvol 
    // où Gvol = (P^n - P^{n-1})/∆t -∆t/2 Mh^-1 ( K_h P^n - Ch^T Pu Lambda^{n-1/2} (ordre 2)
    int Nacous = var_vibro.nodl_acous;
    VectReal_wp Gvol(Nacous);
    for (int i = 0; i < Nacous; i++)
      Gvol(i) = Gvol_minus(i) - Gvol_plus(i);

    chrono.Start(VirtualTimer::FLUX);    
    VectReal_wp tmp; tmp.Reallocate(nFEM); tmp.Fill(0.0);
    VectReal_wp  Fvol; Fvol.Reallocate(nev); Fvol.Fill(0.0);
    var_vibro.MltAddCh(1.0, Gvol, tmp);
    
    chrono.Stop(VirtualTimer::FLUX); chrono.Start(VirtualTimer::PROD);
    
    Mlt(SeldonTrans, Pu, tmp, Fvol);

    chrono.Stop(VirtualTimer::PROD); chrono.Start(VirtualTimer::COMM);
    // on assemble Fvol sur tous les processeurs
    this->AssembleModalVector(Fvol);
    
    chrono.Stop(VirtualTimer::COMM); chrono.Start(VirtualTimer::SCHEME);
    VectReal_wp LambdaTmp(nev);
    LambdaTmp.Fill(0.0);
    for (int g = 0; g < nev; g++)
      LambdaTmp(g) = MatriceResolExacte0(g) * Lambda(g) + MatriceResolExacte1(g) * LambdaPoint(g);
    
    // computation of Lambda^n+1/2
    LambdaNew.Reallocate(nev);
    for (int g = 0; g < nev; g++)
      LambdaNew(g) = (Fvol(g) + RhsMecanicModal(g)*pulse)*MatriceResolExacteSecondMembre(g)
        + LambdaTmp(g) + coef * Lambda0(g);
		
    // SolveSchur
    // on inverse par I + dt/2 R_{dt} Pu^T Cbarre Mh^-1 Cbarre^T
    Copy(LambdaNew, Lambda);
    SolveSchur(Lambda);
    
    // on calcule 1/2 P^T Ch Mh^-1 Ch^T P Lambda^n+1/2 = (LambdaNew - Lambda^{n+1/2}) / R_{dt}
    for (int g = 0; g < nev; g++)
      LambdaTmp(g) = (LambdaNew(g) - Lambda(g)) / MatriceResolExacteSecondMembre(g);
    
    // on met a jour Fvol pour contenir P^T Ch (P^{n+1/2} + P^{n-1/2})/2
    for (int g = 0; g < nev; g++)
      Fvol(g) -= LambdaTmp(g);
    
    // on calcule LambdaPoint^n+1/2
    for (int g = 0; g < nev; g++)
    {
      LambdaPoint(g) = MatriceResolExacte0Point(g) * Lambdaold(g)
          + MatriceResolExacte1Point(g) * LambdaPointold(g)
          + MatriceResolExacteSecondMembrePoint(g)*(Fvol(g) + pulse*RhsMecanicModal(g))
          + coef*Lambda0Point(g);        
    }       

    chrono.Stop(VirtualTimer::SCHEME);
  }
  
  
  //! Calcul de scal_a
  template<class VibroPb>
  Real_wp Soundboard::ComputeScal_a(const VectReal_wp& Gvol_minus, 
                                    const VectReal_wp& Gvol_plus, 
                                    VectReal_wp& rhs, VectReal_wp& Fvol, VibroPb& var_vibro)
  {
    // calcul de Fvol = R_{∆t} tr(Ph) Cbarre Gvol 
    // où Gvol = (P^n - P^{n-1})/∆t -∆t/2 Mh^-1 ( K_h P^n - Ch^T Pu Lambda^{n-1/2} (ordre 2)
    int Nacous = var_vibro.nodl_acous;
    VectReal_wp Gvol(Nacous);
    for (int i = 0; i < Nacous; i++)
      Gvol(i) = Gvol_minus(i) - Gvol_plus(i);
    
    chrono.Start(VirtualTimer::FLUX);
    VectReal_wp tmp; tmp.Reallocate(nFEM); tmp.Fill(0.0);
    Fvol.Reallocate(nev); Fvol.Fill(0.0);
    var_vibro.MltAddCh(1.0, Gvol, tmp);
    
    chrono.Stop(VirtualTimer::FLUX); chrono.Start(VirtualTimer::PROD);
    Mlt(SeldonTrans, Pu, tmp, Fvol);
    
    chrono.Stop(VirtualTimer::PROD); chrono.Start(VirtualTimer::COMM);
    // on assemble Fvol sur tous les processeurs
    this->AssembleModalVector(Fvol);
    
    chrono.Stop(VirtualTimer::COMM); chrono.Start(VirtualTimer::SCHEME);
    VectReal_wp LambdaTmp(nev);
    LambdaTmp.Fill(0.0);
    for (int g = 0; g < nev; g++)
      LambdaTmp(g) = MatriceResolExacte0(g) * Lambda(g) + MatriceResolExacte1(g) * LambdaPoint(g);

    // computation of Lambda^n+1/2
    LambdaVolLambda.Reallocate(nev);
    for (int g = 0; g < nev; g++)
      LambdaVolLambda(g) = Fvol(g)*MatriceResolExacteSecondMembre(g) + LambdaTmp(g);
		
    // SolveSchur
    // on inverse par I + dt/2 R_{dt} Pu^T Cbarre Mh^-1 Cbarre^T
    // pour obtenir LambdaTilde^n+1/2
    rhs = LambdaVolLambda;
    SolveSchur(LambdaVolLambda);
    
    // on calcule scal_a
    Real_wp    scal_a = 0.0;    
    for (int g = 0; g < nev; g++)
    {
      Real_wp Force = (LambdaVolLambda(g) - Lambda(g))*VectSecondMembreModal(g);
      scal_a += Force;
    }
    
    chrono.Stop(VirtualTimer::SCHEME);
    return scal_a;
  }

  
  //! Calcul de la factorisation du complement de Schur
  template<class VibroVar>
  void Soundboard::CalculeLU(VectReal_wp& Mh, const VibroVar& var)
  { 
    // square root of R_{\Delta t}
    RDeltaDemi.Reallocate(nev);
    for (int g = 0; g < nev; g++)
    {
      Real_wp tmp = MatriceResolExacteSecondMembre(g);
      if (tmp < -epsilon_machine)
      {
        cout << "Matrice Second Membre negative pour g = "<< g<<endl;
        abort();
      }
      else if (tmp < epsilon_machine)
      {
        RDeltaDemi(g) = 0;
      }
      else
      {
        RDeltaDemi(g) = sqrt(abs(tmp));
      }
    }

    // computation of LUSchurP column per column
    // LuSchurP = Id + 1/2 R_{∆t}^{1/2} tr(P_h) Ch Mh^-1 tr(Cbarre) P_h R_{∆t}^{1/2}	
    int Nacous = var.nodl_acous;
    VectReal_wp Ones(nev), LuOnes(nev);
    Ones.Fill(0); LuOnes.Fill(0);
    VectReal_wp Pu_Ones(nFEM), Y(nFEM); Pu_Ones.Fill(0);
    VectReal_wp X(Nacous); X.Fill(0);

#ifdef SELDON_WITH_SCALAPACK
    /* LUSchurP.Init(global_blacs_handle, nev, nev,
                  size_block_schur, size_block_schur);
    Real_wp val;
    LUSchurP.Fill(0);
    
    Xloc_schur.Init(global_blacs_handle, nev, 1,
                    size_block_schur, size_block_schur);
    
                    Xloc_schur.Fill(0); */
    LUSchurP.Reallocate(nev, nev); LUSchurP.SetIdentity();
#else
    LUSchurP.Reallocate(nev, nev); LUSchurP.SetIdentity();
#endif
    
    for (int g = 0; g < nev; g++)
    {
	    Ones.Fill(0);
	    Ones(g) = RDeltaDemi(g);
	
    	// multiplication by Pu
    	Mlt(Pu, Ones, Pu_Ones);
    	
    	// multiplication by Ch^*
    	X.Fill(0);
    	var.MltAddChTranspose(1.0, Pu_Ones, X);
    	
    	var.AssembleAcoustic(X);
    	
    	// inversion by Mh
    	for (int i = 0; i < X.GetM(); i++)
    	  X(i) *= Mh(i);
    	
    	// multiplication by Ch
    	Y.Fill(0);
    	var.MltAddCh(1.0, X, Y);
    	
    	// and by Pu*
    	Mlt(SeldonTrans, Pu, Y, LuOnes);
    	
    	// and RDeltaDemi
    	for (int i = 0; i < nev; i++)
    	  LuOnes(i) *= RDeltaDemi(i);
    	
    	// assembling
    	this->AssembleModalVector(LuOnes);

#ifdef SELDON_WITH_SCALAPACK	
        /*        for (int i = 0; i < nev; i++)
          {
            val =0.0;
            if (i == g)
              val = 1.0;
            
            val += 0.5*LuOnes(i);
            LUSchurP.SetGlobal(i, g, val);
            } */

    	for (int i = g; i < nev; i++)
    	  LUSchurP(i, g) += 0.5*LuOnes(i);
#else
    	for (int i = g; i < nev; i++)
    	  LUSchurP(i, g) += 0.5*LuOnes(i);
#endif
    }
    
    //LUSchurP.WriteText("LUSchur"+to_str(rank_proc)+"N"+to_str(nb_proc)+".dat");
    
    // calcul de sa facto de Cholesky et stockage dans LUSchurP
    GetCholesky(LUSchurP);		
  }
	
  
  //! on ecrase X par Y = (I + dt/2 R_{dt} Pu^T Cbarre Mh^{-1} Cbarre^T Pu)^-1 X
  void Soundboard::SolveSchur(VectReal_wp & X)
  {
	
    // la facto de Cholesky etant base sur la matrice symetrique :
    // I + dt/2 sqrt(R_dt) Pu^T Cbarre Mh^{-1} Cbarre^T Pu sqrt(R_dt)
    // on notera que la solution Y peut etre obtenue comme
    // Y = sqrt(R_{dt}) [ I + dt/2 sqrt(R_dt) Pu^T Cbarre Mh^{-1} Cbarre^T Pu sqrt(R_dt) ]
    //     sqrt(R_{dt})^-1 X
    // on divise donc par sqrt(R_{dt})
    for (int g = 0; g < nev; g++)
      X(g) /= RDeltaDemi(g);
    
    // puis on resout par [ I + dt/2 sqrt(R_dt) Pu^T Cbarre Mh^{-1} Cbarre^T Pu sqrt(R_dt) ]
#ifdef SELDON_WITH_SCALAPACK
    /* DistributeVector(X, Xloc_schur);
    
    SolveCholesky(LUSchurP, Xloc_schur);
    
    GatherVector(Xloc_schur, X); */

    SolveCholesky(SeldonNoTrans, LUSchurP, X); // descente
    SolveCholesky(SeldonTrans, LUSchurP, X); // remontee
#else
    SolveCholesky(SeldonNoTrans, LUSchurP, X); // descente
    SolveCholesky(SeldonTrans, LUSchurP, X); // remontee
#endif

    // et on multiplie par sqrt(R_{dt})
    for (int g = 0; g < nev; g++)
      X(g)*=RDeltaDemi(g);

    //paf
  }
  
  
  //! sortie des instantanes (deplacement u et rotations theta)
  void Soundboard::WriteOutputFEM(string toto, EllipticProblem<ReissnerMindlinEquation> & var_soundboard)
  {
    abort();
    // rajouter temps dans les arguments pour que ca marche, merci
    /*bool snapshot = false;
    for (int i = 0; i < var_soundboard.output_grid_param.GetM(); i++)
      if (var_soundboard.output_grid_param(i).SnapshotToStore(temps))
        snapshot = true;
    
    if (snapshot)
      {
        VectReal_wp LambdaFEM(3*nFEM);
        LambdaFEM.Fill(Real_wp(0));
        Mlt(eigen_vectors,Lambda,LambdaFEM);
        var_soundboard.WriteOutputFile(LambdaFEM, toto);
        }*/
  }
  
  
  //! initialisation des sorties
  void Soundboard::InitOutput(EllipticProblem<ReissnerMindlinEquation> & var_soundboard)
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = VarSoundboard.comm_group_mode;
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int rank_proc(0);
#endif

    string file_name; 
    for (int i = 0; i < var_soundboard.output_grid_param.GetM(); i++)
      if (var_soundboard.output_grid_param(i).GetNature()
          == GridInterpolationFull<Dimension2>::POINT)
        file_name = var_soundboard.output_grid_param(i).GetDiffractedFieldFile();
        
    if ( file_name.size() > 0)
      if (rank_proc == 0)
        output_mecanic_sismo.Init(file_name, size_buffer_mecanic_sismo, false);
  }
  

  //! sortie des instantanes (deplacement u et rotations theta)
  void Soundboard::WriteOutputFEM(int nt, Real_wp temps,
				  EllipticProblem<ReissnerMindlinEquation> & var_soundboard)
  {    
    
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = VarSoundboard.comm_group_mode;
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int rank_proc(0);
#endif

    bool snapshot = false; string file_name; 
    for (int i = 0; i < var_soundboard.output_grid_param.GetM(); i++)
      if (var_soundboard.output_grid_param(i).SnapshotToStore(temps))
      {
        if (var_soundboard.output_grid_param(i).GetNature()
              == GridInterpolationFull<Dimension2>::POINT)
        {
          file_name = var_soundboard.output_grid_param(i).GetDiffractedFieldFile();
          var_soundboard.output_grid_param(i).IncrementSnapshot();
        }
        else
          snapshot = true;
      }

    for (int i = 0; i < var_soundboard.output_mesh_param.GetM(); i++)
      if (var_soundboard.output_mesh_param(i).SnapshotToStore(temps))
        snapshot = true;
    
    if ( file_name.size() > 0)
      if (rank_proc == 0)
      {
        VectReal_wp val(1+ValLambdaSismoPoint.GetM());
        val.Fill(0);
        val(0) = temps;
        for (int nb = 0; nb < val.GetM()-1; nb++)
          for (int k = 0; k < nev; k++)
            val(nb+1) += ValLambdaSismoPoint(nb, k)*Lambda(k);
          
        output_mecanic_sismo.AddVect(val);
      }
    
    if (snapshot)
    {
      VectReal_wp LambdaFEM(3*nFEM);
      LambdaFEM.Fill(Real_wp(0));
      Mlt(eigen_vectors, Lambda, LambdaFEM);
      var_soundboard.WriteSnapshot(nt, temps, LambdaFEM, 1);
    }
  }


  //! sortie des instantanes (deplacement u et rotations theta)	
  void Soundboard::WriteOutputVect(VectReal_wp& Vect, string toto,
                                   EllipticProblem<ReissnerMindlinEquation> & var_soundboard)
  {    
    VectReal_wp VectFEM(3*nFEM);
    VectFEM.Fill(Real_wp(0));
    Mlt(eigen_vectors, Vect, VectFEM);
    var_soundboard.WriteOutputFile(VectFEM, toto);
  }


  //#endif
}

#define MONTJOIE_FILE_SOUNDBOARD_CXX
#endif 
