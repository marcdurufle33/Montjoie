#ifndef MONTJOIE_FILE_MEGA_NEWTON_SOLVER_CXX

namespace Montjoie
{

  void LinearSolverNewtonSparse::Factorize(Matrix<Real_wp, General, ArrayRowSparse>& A)
  {
    mat_sparse_lu.Factorize(A, true, true);
  }

  void LinearSolverNewtonSparse::Factorize(Matrix<Real_wp, General, ArrowCol>& A)
  {
    cout << "Incompatible solver" << endl;
    abort();
  }

  void LinearSolverNewtonSparse::Solve(VectReal_wp& b)
  {
    mat_sparse_lu.Solve(b);
  }

  void LinearSolverNewtonArrow::Factorize(Matrix<Real_wp, General, ArrayRowSparse>& A)
  {
    cout << "Not implemented" << endl;
    abort();

    //Copy(A, mat_arrow_lu); A.Clear();
    //mat_arrow_lu.Factorize();
  }

  void LinearSolverNewtonArrow::Factorize(Matrix<Real_wp, General, ArrowCol>& A)
  {
    mat_arrow_lu = A; A.Clear();
    mat_arrow_lu.Factorize();
  }

  void LinearSolverNewtonArrow::Solve(VectReal_wp& b)
  {
    mat_arrow_lu.Solve(b);
  }


  template<int d>
  void LinearSolverNewtonTinyBand<d>::Factorize(Matrix<Real_wp, General, ArrayRowSparse>& A)
  {
    Copy(A, mat_band_lu); A.Clear();
    mat_band_lu.Factorize();
  }

  template<int d>
  void LinearSolverNewtonTinyBand<d>::Factorize(Matrix<Real_wp, General, ArrowCol>& A)
  {
    cout << "Not implemented" << endl;
    abort();

    //Copy(A, mat_band_lu); A.Clear();
    //mat_band_lu.Factorize();
  }

  template<int d>
  void LinearSolverNewtonTinyBand<d>::Solve(VectReal_wp& b)
  {
    mat_band_lu.Solve(b);
  }

  template<int d, int m>
  void LinearSolverNewtonTinyArrow<d, m>::Factorize(Matrix<Real_wp, General, ArrayRowSparse>& A)
  {
    Copy(A, mat_arrow_lu); A.Clear();
    mat_arrow_lu.Factorize();
  }

  template<int d, int m>
  void LinearSolverNewtonTinyArrow<d, m>::Factorize(Matrix<Real_wp, General, ArrowCol>& A)
  {
    cout << "Not implemented" << endl;
    abort();

    //Copy(A, mat_arrow_lu); A.Clear();
    //mat_arrow_lu.Factorize();
  }

  template<int d, int m>
  void LinearSolverNewtonTinyArrow<d, m>::Solve(VectReal_wp& b)
  {
    mat_arrow_lu.Solve(b);
  }

  //! constructeur par defaut
  MegaNewtonSolver::MegaNewtonSolver()
  {
    linear_solver = NULL;
    tol = 1e-6;
    nbitermax = 5;
    erreur = 0.0;
    scale_newton_test = true;
    force_reevaluation_jacobian = false;
    threshold_norme_scheme = 1e-5;
    t_diff = -1;
    scheme_is_linear = false;

#ifdef SELDON_WITH_PASTIX
    // si on utilise Pastix, on raffine la solution
    mat_sparse_lu.RefineSolution();
#endif
    
    size_band_lu = 0;
    size_last_row_lu = 0;
    type_matrix_lu = LU_SPARSE;
  }
  

  MegaNewtonSolver::~MegaNewtonSolver()
  {
    if (linear_solver != NULL)
      delete linear_solver;
  }

  
  //! returns a reference to the vector containing F(X)
  VectReal_wp& MegaNewtonSolver::GetScheme()
  {
    return scheme;
  }
  
  
  // A supprimer apres debug
  /*
  template<class GenericPb>
  void MegaNewtonSolver::RemplitTest(GenericPb & var)
  {
    int size = var.size;
    DiffMatrix.Reallocate(size,size);
    for (int i = 0; i < var.nb_corde; i++)
      {
  for (int ne = 0; ne < var.piano_strings(i).mesh.GetNbElt(); ne++)
    for (int pp = 0; pp < var.piano_strings(i).mesh.GetOrder()+1; pp++)
      {
        int p = var.piano_strings(i).mesh.GetNumberDof(ne,pp);
        for (int k = 0; k < var.dimension; k++)
    {
      var.SetDiffStringString(i, k, p, p, i, k, p, DiffMatrix);
      var.SetDiffLMString(i, 10, i, k, p, DiffMatrix);
      var.SetDiffStringLM(i, k, p, -10, i, DiffMatrix);
      var.SetDiffHammerString(20, i, k, p, DiffMatrix);
      var.SetDiffStringHammer(i, k, p, -20, DiffMatrix);
    }
      }
  
  var.SetDiffLMLM(i, 5, var.nb_corde-i-1, DiffMatrix);
      }
    
    var.SetDiffHammerHammer(12, DiffMatrix);
  }
  */
  
  
  //! on change les parametres de Newton avec une ligne du fichier de donnees 
  /*!
    \param[in] description_field mot-cle du fichier de donnees
    \param[in] parameters liste des parametres associes au mot-cle
  */
  void MegaNewtonSolver
  ::SetInputData(const string & description_field, const Vector<string> &parameters)
  {
    if (!description_field.compare("NewtonParam"))
      {
  if (parameters.GetM() <= 1)
    {
      cout << "In SetInputData of MegaNewtonSolver" << endl;
      cout << "NewtonParam needs at least two parameters, for instance :" << endl;
      cout << "NewtonParam = tol nb_iter_max" << endl;
      cout << "Current parameters are : " << endl << parameters << endl;
      abort();
    }
  
  tol = to_num<Real_wp>(parameters(0)); // critere d'arret de Newton
  nbitermax = to_num<int>(parameters(1)); // nombre max d'iterations
  
  // critere pour determiner si Newton a converge ou pas
  if (parameters.GetM() < 3)
    tol_convergence = tol;
  else
            tol_convergence = to_num<Real_wp>(parameters(2));
        
  // si on met YES_SCALE => scaling des equations
  // si on met NO_SCALE => les equations sont laissees telles quelles
  if (parameters.GetM() < 4)
    scale_newton_test = true;
  else
    {
      if (parameters(3) == "NO_SCALE")
        scale_newton_test = false;
    }
  
  // le residu vaut |F(x)| / norme_ref
  // avec norme_ref = max(|F(x_0)|, threshold_norme_scheme)
  // ca permet de ne pas s'arreter au debut lorsque la solution est toute petite
  if (parameters.GetM() >= 5)
    threshold_norme_scheme = to_num<Real_wp>(parameters(4));
      }
    else if (!description_field.compare("OutputDiffMatrix"))
      {

  if (parameters.GetM() <= 0)
    {
      cout << "In SetInputData of MegaNewtonSolver" << endl;
      cout << "OUtputDiffMatrix needs one parameter, for instance :" << endl;
      cout << "OutputDiffMatrix = 1e-5" << endl;
      cout << "Current parameters are : " << endl << parameters << endl;
      abort();
    }
        t_diff = to_num<Real_wp>(parameters(0));
      }
    else if (!description_field.compare("ForceReevaluationJacobian"))
      {
  if (parameters.GetM() <= 0)
    {
      cout << "In SetInputData of MegaNewtonSolver" << endl;
      cout << "ForceReevaluationJacobian needs one parameter, for instance :" << endl;
      cout << "ForceReevaluationJacobian = YES" << endl;
      cout << "Current parameters are : " << endl << parameters << endl;
      abort();
    }
  
  if (parameters(0) == "YES")
    force_reevaluation_jacobian = true;
  else
    force_reevaluation_jacobian = false;
      }
  }
  
  
  //! computes and factorizes the differential matrix
  template<class GenericPb>
  void MegaNewtonSolver
  ::ComputeAndFactorizeDiff(GenericPb& var, Matrix<Real_wp>& x, bool compute_df)
  {
    if (compute_df)
    {
      //if (type_matrix_lu == LU_SPARSE)
      if (type_matrix_lu != LU_ARROW)
      {
        var.ComputeDiff(DiffMatrix, x);
        ScaleLeftMatrix(DiffMatrix, scale_newton);
      }
      else
      {
        var.ComputeDiff(mat_arrow, x);
        ScaleLeftMatrix(mat_arrow, scale_newton);
      }
    }
    
    if (var.print_level >= 6)
      chrono_string.Start(VirtualTimer::FACTO);
    
    //if (type_matrix_lu == LU_SPARSE)
    if (type_matrix_lu != LU_ARROW)
    {
      linear_solver->Factorize(DiffMatrix);
    }
    else
      linear_solver->Factorize(mat_arrow);
    
    if (var.print_level >= 6)
      chrono_string.Stop(VirtualTimer::FACTO);
  }


  //! computes and factorizes the differential matrix
  template<class GenericPb>
  void MegaNewtonSolver
  ::UpdateAndFactorizeDiff(GenericPb& var, Matrix<Real_wp>& x)
  {
    //if (type_matrix_lu == LU_SPARSE)
    if (type_matrix_lu != LU_ARROW)
    {
      DiffMatrix = RowDiffMatrix;
      var.UpdateDiff(DiffMatrix, x);
      ScaleLeftMatrix(DiffMatrix, scale_newton);
      linear_solver->Factorize(DiffMatrix);
    }
    else
    {
      var.UpdateDiff(mat_arrow, x);
      ScaleLeftMatrix(mat_arrow, scale_newton);
      linear_solver->Factorize(mat_arrow);
    }
  }

  
  //! calcul de DF^-1 b, b est remplace par la solution
  void MegaNewtonSolver::SolveDifferential(VectReal_wp& b, int print_level)
  {
    if (print_level >= 6)
      chrono_string.Start(VirtualTimer::SOLVE);
    
    linear_solver->Solve(b);
    
    if (print_level >= 6)
      chrono_string.Stop(VirtualTimer::SOLVE);
  }   
  
  
  //! multiplication by scaling
  void MegaNewtonSolver::ApplyScaling(VectReal_wp& b)
  {
    for (int i = 0; i < scale_newton.GetM(); i++)
      b(i) *= scale_newton(i);
  }
  
  
  //! Initialisation des iterations de Newton
  /*!
    \param[in] var definition du systeme non-lineaire qu'on veut resoudre
    \param[in] x matrice a trois colonnes x(i, 0) = x^{n-1}
    x(i, 1) = x^n et x(i, 2) = x^{n+1}
    
    On alloue ici les tableaux, la jacobienne est calculee et factorisee
  */
  template<class GenericPb>
  void MegaNewtonSolver::Init(GenericPb & var, Matrix<Real_wp> & x)
  {
    int size = x.GetM();
    
    // on desalloue toutes les differentielles eventuellement utilisees anterieurement
    DiffMatrix.Clear();
    mat_arrow.Clear();
    
    // allocation des tableaux
    scheme.Reallocate(size);
    b.Reallocate(size);
    save_scheme.Reallocate(size);
    erreur_relative.Reallocate(size);
    x_memory.Reallocate(size);
    x_prec.Reallocate(size);
    rhs.Reallocate(size);
    scheme.Fill(0);
    b.Fill(0);
    save_scheme.Fill(0);
    erreur_relative.Fill(0);
    x_memory.Fill(0);
    x_prec.Fill(0);
    rhs.Fill(0);

    // the linear solver is constructed
    if (linear_solver != NULL)
      delete linear_solver;
    
    linear_solver = NULL;
    if (type_matrix_lu == LU_SPARSE)
    {
      linear_solver = new LinearSolverNewtonSparse();
      cout << "A direct solver has been selected" << endl;
    }
    else if (type_matrix_lu == LU_ARROW)
    {
      linear_solver = new LinearSolverNewtonArrow();
      cout << "An arrow solver has been selected" << endl;
    }
    else
    {
      // tiny arrow matrices
      if (size_last_row_lu == 0)
        switch (size_band_lu)
        {
          case 1: linear_solver = new LinearSolverNewtonTinyBand<1>(); break;
          case 2: linear_solver = new LinearSolverNewtonTinyBand<2>(); break;
          case 3: linear_solver = new LinearSolverNewtonTinyBand<3>(); break;
          case 4: linear_solver = new LinearSolverNewtonTinyBand<4>(); break;
          case 6: linear_solver = new LinearSolverNewtonTinyBand<6>(); break;
        }
      else if (size_last_row_lu == 1)
        switch (size_band_lu)
        {
          case 1: linear_solver = new LinearSolverNewtonTinyArrow<1, 1>(); break;
          case 2: linear_solver = new LinearSolverNewtonTinyArrow<2, 1>(); break;
          case 3: linear_solver = new LinearSolverNewtonTinyArrow<3, 1>(); break;
          case 4: linear_solver = new LinearSolverNewtonTinyArrow<4, 1>(); break;
          case 6: linear_solver = new LinearSolverNewtonTinyArrow<6, 1>(); break;
        }
      else if (size_last_row_lu == 2)
        switch (size_band_lu)
        {
          case 1: linear_solver = new LinearSolverNewtonTinyArrow<1, 2>(); break;
          case 2: linear_solver = new LinearSolverNewtonTinyArrow<2, 2>(); break;
          case 3: linear_solver = new LinearSolverNewtonTinyArrow<3, 2>(); break;
          case 4: linear_solver = new LinearSolverNewtonTinyArrow<4, 2>(); break;
          case 6: linear_solver = new LinearSolverNewtonTinyArrow<6, 2>(); break;
        }
      else if (size_last_row_lu == 3)
        switch (size_band_lu)
        {
          case 1: linear_solver = new LinearSolverNewtonTinyArrow<1, 3>(); break;
          case 2: linear_solver = new LinearSolverNewtonTinyArrow<2, 3>(); break;
          case 3: linear_solver = new LinearSolverNewtonTinyArrow<3, 3>(); break;
          case 4: linear_solver = new LinearSolverNewtonTinyArrow<4, 3>(); break;
          case 6: linear_solver = new LinearSolverNewtonTinyArrow<6, 3>(); break;
        }

      if (linear_solver == NULL)
      {
        cout << "An arrow solver has been selected" << endl;
        type_matrix_lu = LU_ARROW;
        linear_solver = new LinearSolverNewtonArrow();
      }
      else
        cout << "A tiny arrow solver has been selected" << endl;
    }
        
    // calcul de la jacobienne DF(x)
    //if (type_matrix_lu == LU_SPARSE)
    if (type_matrix_lu != LU_ARROW)
    {
      DiffMatrix.Reallocate(size, size);
      RowDiffMatrix.Reallocate(size, size);
      var.ComputeDiff(DiffMatrix, x);
      RowDiffMatrix = DiffMatrix; // on sauvegarde la diff matrix non scalée et non factorisée
    }
    else
    {
      mat_arrow.Reallocate(size, size, size_band_lu, size_band_lu,
                             size_last_row_lu, size_last_row_lu);
      mat_arrow.Fill(Real_wp(0));
      var.ComputeDiff(mat_arrow, x); 
    }

    // on calcule le scaling si necessaire
    if (scale_newton_test)
    {
      if (type_matrix_lu == LU_SPARSE)
        GetRowSum(scale_newton, DiffMatrix);
      else
        GetRowSum(scale_newton, mat_arrow);
        
      //DISP(size_last_row_lu);
      //mat_arrow2_lu.WriteText("arrow2.dat");
      //scale_newton.WriteText("scale.dat");
      for (int i = 0; i < scale_newton.GetM(); i++)
      {
        if (scale_newton(i) == 0.0)
        {
          cout << "Differential matrix not invertible" << endl;
          abort();
        }
        else
          scale_newton(i) = 1.0/scale_newton(i);
      }
  
      //if (type_matrix_lu == LU_SPARSE)
      if (type_matrix_lu != LU_ARROW)
        ScaleLeftMatrix(DiffMatrix, scale_newton);
      else
        ScaleLeftMatrix(mat_arrow, scale_newton);
    }
    else
    {
      scale_newton.Reallocate(size);
      scale_newton.Fill(1.0);
    }
    
    //cout << "On compare les matrices, pour le moment non" << endl;
    //this->CompareDiffMatrix(var, x);
    
    // factorisation de la jacobienne
    ComputeAndFactorizeDiff(var, x, false);

    cout << " INIT OK" << endl;
    // fichier stockant les residus
    if (var.print_level >= 10)
      file_out_residual_newton.open("residu_newton.txt");
  }
  
  
  //! resolution du system F(x) = 0
  /*!
    \param[in] var definition du systeme non-lineaire qu'on veut resoudre
    \param[in] x matrice a trois colonnes x(i, 0) = x^{n-1}
    x(i, 1) = x^n et x(i, 2) = x^{n+1}
    en pratique seul x^{n+1} est l'inconnue recherchee,
    les vecteurs x^{n-1} et x^n sont juste des parametres de F
    
  */
  template<class GenericPb>
  void MegaNewtonSolver::Solve(GenericPb & var, Matrix<Real_wp> & x)
  { 
    // var contient deux membres : CalculeSchemeLM (F) et CalculeDiffLM (dF)
    // Solve remplace x par la sol de F(x) = 0
    erreur = 1.0;
    int nb = 0;

    if (nbitermax == 1)
    {
      scheme_is_linear = true;
    }
    
    // on evalue scheme = F(x)
    var.ComputeRightHandSide(rhs, x, true);
    var.ComputeScheme(scheme, rhs, x);
    
    for (int i = 0; i < scheme.GetM(); i++)
      scheme(i) *= scale_newton(i);
    
    // on va comparer |F(x)| avec max(|F(x_0)|, threshold_norme_scheme)
    // le threshold sert a eviter les problemes lorsque la solution x est
    // toute petite, ce qui arrive au debut des iterations en temps
    Real_wp norme_ref = max(threshold_norme_scheme, Norm2(scheme));
    
    // on memorise F(x_0) et x_0
    for (int i = 0; i < x.GetM(); i++)
    {
      save_scheme(i) = scheme(i);
      x_memory(i) = x(i, 2);
    }
    
    // on effectue les iterations de Newton
    Newton(var, x, nb, norme_ref);

    if(abs(var.t_courant - t_diff) < 1e-10)
    {    
      CompareDiffMatrix(var, x);
      abort();
    }
    if ( (!scheme_is_linear) && (nb == nbitermax) )
    {
      // DISP("On recommence");
      // l'algo de Newton a atteint le nombre maximal d'iterations
      // on considere alors qu'il faut remettre a jour la jacobienne
      // on recalcule la jacobienne et on recommence du debut
      for (int i = 0; i < x.GetM(); i++)
      {
        x(i, 2) = x_memory(i);
        scheme(i) = save_scheme(i);
      }
  
      ComputeAndFactorizeDiff(var, x);
        
      nb = 0;
      Newton(var, x, nb, norme_ref);
    }
    
    if (var.print_level >= 10)
    {
      file_out_residual_newton << nb << " " << erreur << " " << norme_ref << endl;
    }
    //  scheme.WriteText(var.DOSSIER + "scheme" + to_str(var.nb_iter_courant) + ".txt");
    // Si on n'a pas converge, on arrete le calcul
    if ( (!scheme_is_linear) && (erreur > tol_convergence*norme_ref) )
    {
      cout << "help !! pb de convergence "<<endl;
      DISP(erreur/norme_ref);
      DISP(tol_convergence);
      DISP(nbitermax);
      scheme.Write("scheme.dat");
      // on ecrit la matrice jacobienne exacte et numerique
      CompareDiffMatrix(var, x);
        
      // les eventuels buffers de sortie sont vides
      var.FlushBuffer();
      abort();
    }
    
    //int test_input; cout << "waiting" << endl; cin >> test_input;
    
    // on affiche le nombre d'iterations effectuees
    
    int n = toInteger(round(var.t_courant/(var.t_display)));
    if ((abs(var.t_courant/(var.t_display) - n) < var.Deltat/2))
    {
      if(nb==0)
        cout << "Newton a fait " << nbitermax << " iterations" <<endl;
      else
      {
        if(nb==1)
          cout << "Newton a fait " << nb << " iteration" <<endl;
        else 
          cout << "Newton a fait " << nb << " iterations" <<endl;
      }
    }
    
    // Debug
    //	  int nodl = var.piano_strings(0).mesh.GetNbDof()-1;
    //	  Real_wp perte = var.GetLM(0,2,x)*(var.GetValue(0,0,nodl,2,x)
    //                     -var.GetValue(0,0,nodl,0,x))*0.5*var.invDeltat;
    //      WriteOnTheGo(var.DOSSIER  + "pertes.txt", perte);
    //
    
    // on passe au temps suivant x^{n-1} = x^n et x^n = x^{n+1}
    // Mise a jour des valeurs de x
    for(int n = 0; n < x.GetM(); n++)
    {
      x(n, 0) = x(n, 1);
      x(n, 1) = x(n, 2);
    }

    if ((var.hammer.with_auxiliary_variable) && (!force_reevaluation_jacobian))
    {
      //Pour mettre à jour uniquement les termes de la Jacobienne qui ont changé
      //depuis l'itération précédante
      UpdateAndFactorizeDiff(var, x);
    }


  }
  
  
  //! Iterations de Newton pour resoudre F(x) = 0
  /*!
    \param[in] var definition du systeme non-lineaire qu'on veut resoudre
    \param[in] x matrice a trois colonnes x(i, 0) = x^{n-1}
    x(i, 1) = x^n et x(i, 2) = x^{n+1}
    en pratique seul x^{n+1} est l'inconnue recherchee,
    les vecteurs x^{n-1} et x^n sont juste des parametres de F
    
    \param[out] nb nombre d'iterations necessaires pour obtenir la solution
    \param[in] norme_ref on s'arrete lorsque |F(x)| / norme_ref <= tol
    \param[in] reevaluate_jacobian si vrai, la jacobienne est reevaluee a chaque iteration
    si non, on garde la meme jacobienne jusqu'a la fin
    x contient le guess initial en entree, la solution de F(x) = 0 en sortie
  */
  template<class GenericPb>
  void MegaNewtonSolver::Newton(GenericPb & var, Matrix<Real_wp> & x, int &nb,
     const Real_wp& norme_ref, bool reevaluate_jacobian)
  {

    Real_wp err_prec = 1e20;
    erreur = err_prec*0.99;
    
    // on verifie que le newton ne stagne pas, ou que le residu est assez petit
    //DISP(tol*norme_ref); DISP(norme_ref);
    while ((erreur < err_prec) && (erreur > tol*norme_ref) && (nb < nbitermax))
    {          
      nb++;
      err_prec = erreur;          
      for (int r = 0; r < x.GetM(); r++)
      {
        x_prec(r) = x(r, 2);
        b(r) = -scheme(r);
      }
  
      if (force_reevaluation_jacobian)
        reevaluate_jacobian = true;
  
      if (reevaluate_jacobian)
      {
        // calcul et factorisation de la jacobienne
        ComputeAndFactorizeDiff(var, x);
      }  
  
      // on calcule b = - DF^{-1} F(x)
      SolveDifferential(b, var.print_level);
        
      // nouvel itere x - DF^{-1} F(x)
      for(int n = 0; n < x.GetM(); n++)
        x(n, 2) += b(n);
  
      // on force les contraintes s'il y en a
      var.ForceConstraint(x);
        
	// on calcul F(x) avec le nouvel itere
	var.ComputeScheme(scheme, rhs, x);
	
	for (int i = 0; i < scheme.GetM(); i++)
	  scheme(i) *= scale_newton(i);
	
	// 	 scheme.WriteText(var.DOSSIER + "Newscheme" + to_str(nb) + ".txt");
	//for(int n=0;n<x.GetM();n++){
	//		  if(abs(scheme(n))<=1e-14){
	//			// magnifique wunderbar wonderful
	//			erreur_relative(n)=old_scheme(n)-scheme(n);
	//		  }
	//		  else 
	//		  {
	//			erreur_relative(n)=(old_scheme(n)-scheme(n))/scheme(n);
	//		  }
	//		}
	//		Y = Y + alpha * X
	//		Add(alpha, X, Y);
	//Add(-1.0,scheme,old_scheme);
	//		DISP(Norm2(scheme));
	//		DISP(Norm2(old_scheme));
	//		erreur=Norm2(old_scheme) / Norm2(scheme);
	// min(Norm2(scheme),Norm2(old_scheme-scheme));
	//		DISP(Norm2(erreur_relative));
	// DISP(Norm2(scheme));
	
	// erreur = | F(x) |
	erreur = Norm2(scheme);
	
	if (((reevaluate_jacobian) && (!force_reevaluation_jacobian)) || (var.print_level >= 10))
	  cout << "Erreur dans le Newton a l'iteration " << nb << " = " << erreur << endl;
	
	/*if ((nb >= 2) && (erreur > 1e-8))
	  {
	  x.WriteText("x_example.dat");
	  CompareDiffMatrix(var, x);
	  abort();
	  }
	*/
	//DISP(erreur);
	//	DISP(err_prec);
	// Norm2(erreur_relative);
	
      }
    
    // si on a |F(x_{k-1})| < |F(x_k)| on retourne x_{k-1}
    if (erreur > err_prec)
    {
      for (int i = 0; i < x.GetM(); i++)
        x(i, 2) = x_prec(i);
  
      if ( (!scheme_is_linear) && (erreur > tol_convergence*norme_ref) )
      {
        nb = nbitermax;
      }
    }
  }
  
  
  //! on ecrit la jacobienne sur le disque
  void MegaNewtonSolver::EcritDiffMatrix(Real_wp t, string DOSSIER)
  {
    // Ecrit la matrice dans un fichier, au format de matlab
    string name_file = DOSSIER + "Diff"+to_str(t)+".txt";
    remove(name_file.data());
    //if (type_matrix_lu == LU_SPARSE)
    if (type_matrix_lu != LU_ARROW)
      DiffMatrix.WriteText(name_file);
    else
      mat_arrow.WriteText(name_file);
  }
  
  
  //! on ecrit F(x) sur le disque
  void MegaNewtonSolver::EcritScheme(Real_wp t, string DOSSIER)
  {
    string name_file = DOSSIER+"/Scheme"+to_str(t)+".txt";
    remove(name_file.data());
    ofstream objetfichier;
    objetfichier.open(name_file.data(), ios::out);//on ouvre le fichier en ecriture
    objetfichier <<scheme << endl;//*
    objetfichier.close();
  }
  
  
  //! comparaison de la jacobienne exacte et de la jacobienne numerique
  //! Les deux matrices sont ecrites sur le disque
  template<class GenericPb>
  void MegaNewtonSolver::CompareDiffMatrix(GenericPb & var,Matrix<Real_wp> & val)
  {
    // fill with random values
    // Matrix<Real_wp> valpouet;
    //    valpouet.Reallocate((val.GetM())/2*3,3);
    //    valpouet.FillRand();
    //    int p=0;
    //    int q=0;
    //    for (int i=0; i<(valpouet.GetM()-3)/3; i++) {
    //    val(p,0)=valpouet(q,0);val(p,1)=valpouet(q,1);val(p,2)=valpouet(q,2);
    //    p++;
    //    q++;
    //    val(p,0)=valpouet(q,0);val(p,1)=valpouet(q,1);val(p,2)=valpouet(q,2);
    //    p++;
    //    q++;
    //    q++;
    //    }
    //    val(p,0)=valpouet(q,0);val(p,1)=valpouet(q,1);val(p,2)=valpouet(q,2);
    //    val(p+1,0)=valpouet(q+1,0);val(p+1,1)=valpouet(q+1,1);val(p+1,2)=valpouet(q+1,2);
    //    DISP(val);
    //    val*=1e-12;
    //    val.WriteText(var.DOSSIER + "/VarRand.txt");
    // end random
    int size = val.GetM();
    //if (type_matrix_lu == LU_SPARSE)
    if (type_matrix_lu != LU_ARROW)
      {
  DiffMatrix.Reallocate(size, size);
  var.ComputeDiff(DiffMatrix, val);
  ScaleLeftMatrix(DiffMatrix, scale_newton);
      }
    else
      {
  mat_arrow.Reallocate(size, size, size_band_lu, size_band_lu,
                             size_last_row_lu, size_last_row_lu);
  
  var.ComputeDiff(mat_arrow, val);
  ScaleLeftMatrix(mat_arrow, scale_newton);
      }
    
    EcritDiffMatrix(var.t_courant, var.DOSSIER);
    
    Matrix<Real_wp, General, ArrayRowSparse> FalseDiffMatrix(val.GetM(), val.GetM());
    
    Real_wp h = pow(epsilon_machine, Real_wp(1.0)/Real_wp(3.0)); DISP(h);
    VectReal_wp scheme_minus(val.GetM()), scheme_plus(val.GetM());
    //      Real_wp a(0), b(0);
    for(int p = 0; p < val.GetM(); p++)
      {
  //boucle sur les variables
  //variations et calcul du scheme
  Real_wp souvenir = val(p, 2);
  val(p, 2) += h;
  var.ComputeScheme(scheme_plus, rhs, val);
        
  val(p, 2) -= 2.0*h;
  var.ComputeScheme(scheme_minus, rhs, val);
        
  for(int q = 0; q < val.GetM(); q++)
    {
      Real_wp DF = (scale_newton(q)*scheme_plus(q)
        - scale_newton(q)*scheme_minus(q)) / Real_wp(2.0*h);
      
      if (abs(DF) > epsilon_machine)
        FalseDiffMatrix.AddInteraction(q, p, DF);
      
      if (abs(DF - DiffMatrix(q, p)) > 1e-5)
        {
    DISP(p); DISP(q);
    cout << "Value of the numerical jacobian : " << DF << endl;
    cout << "Value computed with ComputeDiff/AddDiff : " << DiffMatrix(q, p) << endl;
    DISP(DF-DiffMatrix(q, p));
        }
      
    }
  
  val(p, 2) = souvenir;
      } //for p
    
    string name_file = var.DOSSIER + "FalseDiff"+to_str(var.t_courant)+".txt";
    remove(name_file.data());
    FalseDiffMatrix.WriteText(name_file);
  }
  
}

#define MONTJOIE_FILE_MEGA_NEWTON_SOLVER_CXX
#endif
