void ModeEs_Solver::ResolutionNonLinear_DD_Dirichlet(Vector<Complex_wp>& betas,
                                  Vector<VectComplex_wp>& E3, Vector<VectComplex_wp>& H3,
                                  Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                                  Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Kr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Cr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rir2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sir2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Dr2h,
                                  VectR2& vec_Pos)
{
  int N = this->mesh_num.GetNbDof();

  ModeEsNonLinear_DD_DirichletEigenProblem var_eig(N, *this, vec_Mh, vec_Mr2h, vec_Kh, vec_Kr2h,
                                      vec_Ch, vec_Cr2h, vec_Rjh, vec_Rjr2h, vec_Rih, vec_Rir2h,
                                      vec_Sjh, vec_Sjr2h, vec_Sih, vec_Sir2h, vec_Dh, vec_Dr2h, vec_Pos,
                                      this->DofKeptDir, this->IndexDirichlet,
                                      this->use_split_formulation);

  if (this->use_explicit_matrix)
    var_eig.SetExplicitMatrix(); // en calculant les matrices

  var_eig.SetExactPreconditioning();
  var_eig.SetPrintLevel(4);

  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(nb_asked_eigenval);

  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, shift);
  SlepcParamNep& param = var_eig.GetSlepcParameters();
  if (type_eigensolver == NL_RII)
    param.SetEigensolverType(param.RII);
  else if (type_eigensolver == NL_SLP)
    param.SetEigensolverType(param.SLP);
  else
    param.SetEigensolverType(param.NLEIGS);

  Real_wp reShift = realpart(this->shift);
  Real_wp imShift = imagpart(this->shift);
  Real_wp rayonIntervalle = this->radRegion;

  param.SetIntervalRegion(reShift-rayonIntervalle, reShift+rayonIntervalle, imShift-rayonIntervalle, imShift+rayonIntervalle);

  // param.SetIntervalRegion(0.1, 10.0, 0.001, 10.0);
  param.EnableCommandLineOptions();
  if (this->use_default_petsc_solver)
    param.SetDefaultPetscSolver();

  Vector<Complex_wp> lambda, lambda_imag;
  Matrix<Complex_wp, General, ColMajor> eigen_vec;

  DISP(var_eig.GetM());
  FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);

  DISP(lambda); DISP(lambda_imag);

  int nbFoundModes = eigen_vec.GetN();
  betas.Reallocate(nbFoundModes);
  E3.Reallocate(nbFoundModes);
  H3.Reallocate(nbFoundModes);

  for (int k = 0; k < nbFoundModes; k++)
    {
      int nb_dof_E3 = DofKeptDir.GetM();
      int nb_dof_H3 = N;
      VectComplex_wp x(2*nb_dof_H3), X;
      VectComplex_wp y;
      GetCol(eigen_vec, k, X);

      VectComplex_wp Es, Hs, E1, E2, H1, H2;
      var_eig.ComputeEsHs(lambda(k), X, Es, Hs, E1, E2, H1, H2);

      X.Resize(nb_dof_E3 + nb_dof_H3);
      VectComplex_wp Y(nb_dof_E3 + nb_dof_H3); Y.Zero();
      ModeEsNonLinear2NDirichletEigenProblem var_eig_model(N, *this,
                                          this->DofKeptDir, this->IndexDirichlet,
                                          false);

      var_eig_model.MltOperator(lambda(k), SeldonNoTrans, X, Y);
      cout << "VP " << k << ": " << lambda(k) << " ; ||T(VP) X|| : " << Norm2(Y) << " ; erreur relative : " << Norm2(Y) / Norm2(X) << endl;

      var_eig.ExtractDirichlet(X, x);
      E3(k).Reallocate(nb_dof_H3); E3(k).Zero();
      H3(k).Reallocate(nb_dof_H3); H3(k).Zero();
      for (int i = 0; i < nb_dof_H3; i++)
        {
          E3(k).Get(i) = x(i);
          H3(k).Get(i) = x(nb_dof_H3 + i);
        }
      betas.Get(k) = lambda(k);
      WriteOutputFileScalar(E3(k), "ModeE3_"+ to_str(k));
      WriteOutputFileScalar(H3(k), "ModeH3_"+ to_str(k));
      WriteOutputFileScalar(Es, "ModeEs_"+ to_str(k));
    }
}

void ModeEs_Solver::ResolutionNonLinear_Dk_Dirichlet(Vector<Complex_wp>& betas,
                                  Vector<VectComplex_wp>& E3, Vector<VectComplex_wp>& H3,
                                  Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                                  Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Kr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Cr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rir2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sir2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Dr2h,
                                  VectR2& vec_Pos)
{
  int N = this->mesh_num.GetNbDof();

  ModeEsNonLinear_Dk_DirichletEigenProblem var_eig(N, *this, vec_Mh, vec_Mr2h, vec_Kh, vec_Kr2h,
                                      vec_Ch, vec_Cr2h, vec_Rjh, vec_Rjr2h, vec_Rih, vec_Rir2h,
                                      vec_Sjh, vec_Sjr2h, vec_Sih, vec_Sir2h, vec_Dh, vec_Dr2h, vec_Pos,
                                      this->DofKeptDir, this->IndexDirichlet,
                                      this->use_split_formulation);

  if (this->use_explicit_matrix)
    var_eig.SetExplicitMatrix(); // en calculant les matrices

  var_eig.SetExactPreconditioning();
  var_eig.SetPrintLevel(4);

  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(nb_asked_eigenval);

  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, shift);
  SlepcParamNep& param = var_eig.GetSlepcParameters();
  if (type_eigensolver == NL_RII)
    param.SetEigensolverType(param.RII);
  else if (type_eigensolver == NL_SLP)
    param.SetEigensolverType(param.SLP);
  else
    param.SetEigensolverType(param.NLEIGS);

  Real_wp reShift = realpart(this->shift);
  Real_wp imShift = imagpart(this->shift);
  Real_wp rayonIntervalle = this->radRegion;

  param.SetIntervalRegion(reShift-rayonIntervalle, reShift+rayonIntervalle, imShift-rayonIntervalle, imShift+rayonIntervalle);

  // param.SetIntervalRegion(0.1, 10.0, 0.001, 10.0);
  param.EnableCommandLineOptions();
  if (this->use_default_petsc_solver)
    param.SetDefaultPetscSolver();

  Vector<Complex_wp> lambda, lambda_imag;
  Matrix<Complex_wp, General, ColMajor> eigen_vec;

  DISP(var_eig.GetM());
  FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);

  DISP(lambda); DISP(lambda_imag);

  int nbFoundModes = eigen_vec.GetN();
  betas.Reallocate(nbFoundModes);
  E3.Reallocate(nbFoundModes);
  H3.Reallocate(nbFoundModes);

  for (int k = 0; k < nbFoundModes; k++)
    {
      int nb_dof_E3 = DofKeptDir.GetM();
      int nb_dof_H3 = N;
      VectComplex_wp x(2*nb_dof_H3), X;
      VectComplex_wp y;
      GetCol(eigen_vec, k, X);

      VectComplex_wp Es, Hs, E1, E2, H1, H2;
      var_eig.ComputeEsHs(lambda(k), X, Es, Hs, E1, E2, H1, H2);

      X.Resize(nb_dof_E3 + nb_dof_H3);
      VectComplex_wp Y(nb_dof_E3 + nb_dof_H3); Y.Zero();
      ModeEsNonLinear2NDirichletEigenProblem var_eig_model(N, *this,
                                          this->DofKeptDir, this->IndexDirichlet,
                                          false);

      var_eig_model.MltOperator(lambda(k), SeldonNoTrans, X, Y);
      cout << "VP " << k << ": " << lambda(k) << " ; ||T(VP) X|| : " << Norm2(Y) << " ; erreur relative : " << Norm2(Y) / Norm2(X) << endl;

      var_eig.ExtractDirichlet(X, x);
      E3(k).Reallocate(nb_dof_H3); E3(k).Zero();
      H3(k).Reallocate(nb_dof_H3); H3(k).Zero();
      for (int i = 0; i < nb_dof_H3; i++)
        {
          E3(k).Get(i) = x(i);
          H3(k).Get(i) = x(nb_dof_H3 + i);
        }
      betas.Get(k) = lambda(k);
      WriteOutputFileScalar(E3(k), "ModeE3_"+ to_str(k));
      WriteOutputFileScalar(H3(k), "ModeH3_"+ to_str(k));
      WriteOutputFileScalar(Es, "ModeEs_"+ to_str(k));
      WriteOutputFileScalar(Hs, "ModeHs_"+ to_str(k));
    }
}

void ModeEs_Solver::ResolutionNonLinear2NDirichlet(Vector<Complex_wp>& betas,
                                  Vector<VectComplex_wp>& E3, Vector<VectComplex_wp>& H3)
{
  int N = this->mesh_num.GetNbDof();

  ModeEsNonLinear2NDirichletEigenProblem var_eig(N, *this,
                                      this->DofKeptDir, this->IndexDirichlet,
                                      this->use_split_formulation);

  if (this->use_explicit_matrix)
    var_eig.SetExplicitMatrix(); // en calculant les matrices

  var_eig.SetExactPreconditioning();
  var_eig.SetPrintLevel(4);

  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(nb_asked_eigenval);

  SlepcParamNep& param = var_eig.GetSlepcParameters();
  if (type_eigensolver == NL_RII)
    param.SetEigensolverType(param.RII);
  else if (type_eigensolver == NL_SLP)
    param.SetEigensolverType(param.SLP);
  else
    param.SetEigensolverType(param.NLEIGS);

  param.EnableCommandLineOptions();
  if (this->use_default_petsc_solver)
    param.SetDefaultPetscSolver();

  // cout << "Penser à bien cibler  (shift + région de recherche +-0.1) !!!" << endl;

  /*
    On centre la recherche sur un shift (il devrait donc être différent entre chaque subdivision d'intervalle)
  */
  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, shift);
  /*
    Définition de l'intervalle complexe
    Pas obligé de le faire avec un rayon, là c'est ça qui va permettre de
      subdiviser l'ensemble de recherche de la valeur propre.
  */
  Real_wp reShift = realpart(this->shift);
  Real_wp imShift = imagpart(this->shift);
  Real_wp rayonIntervalle = this->radRegion;
  param.SetIntervalRegion(reShift-rayonIntervalle, reShift+rayonIntervalle, imShift-rayonIntervalle, imShift+rayonIntervalle);
  /*
    C'est là qu'on recherche la valeur et le vecteur propre.
  */
  Vector<Complex_wp> lambda, lambda_imag;
  Matrix<Complex_wp, General, ColMajor> eigen_vec;
  DISP(var_eig.GetM());
  FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);

  /*
    Du coup, ce que je voudrais, c'est faire une recherche parallèle
    avec subdivision du rectangle de recherche en 2*2 ou 3*3 (ou ce qui te semble le plus pratique / efficace).
    Dès qu'un programme s'arrête dans un coeur, a priori ça veut dire que la VP est dans la subdivision pour lequel il s'est arrêté
      (mais ce serait prudent de remplir un fichier "historique de recherche", au cas où il se montrait capricieux).
    Réitérer jusqu'à ce que l'erreur soit assez petite.

    // Pour le calcul d'erreur, à la fin de chaque itération, on peut utiliser ces lignes :
    int nb_dof_E3 = DofKeptDir.GetM();
    int nb_dof_H3 = N;
    VectComplex_wp X; GetCol(eigen_vec, 0, X);
    X.Resize(nb_dof_E3 + nb_dof_H3);
    VectComplex_wp Y(nb_dof_E3 + nb_dof_H3); Y.Zero();

    var_eig.MltOperator(lambda(0), SeldonNoTrans, X, Y);
    // L'erreur à tester : Norm2(Y) / Norm2(X), selon un threshold de l'ordre de 10^(-12) (à voir)
  */




  DISP(lambda); DISP(lambda_imag);

  int nbFoundModes = eigen_vec.GetN();
  betas.Reallocate(nbFoundModes);
  E3.Reallocate(nbFoundModes);
  H3.Reallocate(nbFoundModes);

  for (int k = 0; k < nbFoundModes; k++)
    {
      int nb_dof_H3 = N;
      VectComplex_wp x(2*nb_dof_H3), X;
      VectComplex_wp y;
      GetCol(eigen_vec, k, X);

      VectComplex_wp Es, Hs, E1, E2, H1, H2;
      var_eig.ComputeEsHs(lambda(k), X, Es, Hs, E1, E2, H1, H2);

      var_eig.ExtractDirichlet(X, x);
      E3(k).Reallocate(nb_dof_H3); E3(k).Zero();
      H3(k).Reallocate(nb_dof_H3); H3(k).Zero();
      for (int i = 0; i < nb_dof_H3; i++)
        {
          E3(k).Get(i) = x(i);
          H3(k).Get(i) = x(nb_dof_H3 + i);
        }
      betas.Get(k) = lambda(k);

      cout << "VP " << k << ": " << lambda(k) <<  endl;
      // WriteOutputFileScalar(Es, "ModeEs_"+ to_str(k));
      // WriteOutputFileScalar(Hs, "ModeHs_"+ to_str(k));
      WriteOutputFileScalar(E3(k), "ModeE3_"+ to_str(k));
      WriteOutputFileScalar(H3(k), "ModeH3_"+ to_str(k));
      WriteOutputFileScalar(Es, "ModeEs_"+ to_str(k));
      WriteOutputFileScalar(Hs, "ModeHs_"+ to_str(k));
    }
}

void ModeEs_Solver::ResolutionNonLinear_DD_CLA(Vector<Complex_wp>& betas,
                                  Vector<VectComplex_wp>& E3, Vector<VectComplex_wp>& H3,
                                  Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                                  Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Kr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Cr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Rir2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjr2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Sir2h,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_Dr2h,
                                  VectR2& vec_Pos,
                                  VectReal_wp& Mh_surf,
                                  Matrix<Complexe, Prop, Storage>& Kh_surf,
                                  Matrix<Complexe, Prop, Storage>& Rh_surf)
{
  cout << "À réimplémenter" << endl;
  int N = this->mesh_num.GetNbDof();

  ModeEsNonLinear_DD_CLA_EigenProblem var_eig(N, *this, vec_Mh, vec_Mr2h, vec_Kh, vec_Kr2h,
                                      vec_Ch, vec_Cr2h, vec_Rjh, vec_Rjr2h, vec_Rih, vec_Rir2h,
                                      vec_Sjh, vec_Sjr2h, vec_Sih, vec_Sir2h, vec_Dh, vec_Dr2h, vec_Pos,
                                      this->use_split_formulation);

  if (this->use_explicit_matrix)
    var_eig.SetExplicitMatrix(); // en calculant les matrices

  var_eig.SetExactPreconditioning();
  var_eig.SetPrintLevel(4);

  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(nb_asked_eigenval);

  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, shift);
  SlepcParamNep& param = var_eig.GetSlepcParameters();
  if (type_eigensolver == NL_RII)
    param.SetEigensolverType(param.RII);
  else if (type_eigensolver == NL_SLP)
    param.SetEigensolverType(param.SLP);
  else
    param.SetEigensolverType(param.NLEIGS);

  Real_wp reShift = realpart(this->shift);
  Real_wp imShift = imagpart(this->shift);
  Real_wp rayonIntervalle = this->radRegion;

  param.SetIntervalRegion(reShift-rayonIntervalle, reShift+rayonIntervalle, imShift-rayonIntervalle, imShift+rayonIntervalle);

  // param.SetIntervalRegion(0.1, 10.0, 0.001, 10.0);
  param.EnableCommandLineOptions();
  if (this->use_default_petsc_solver)
    param.SetDefaultPetscSolver();

  Vector<Complex_wp> lambda, lambda_imag;
  Matrix<Complex_wp, General, ColMajor> eigen_vec;

  DISP(var_eig.GetM());
  FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);

  DISP(lambda); DISP(lambda_imag);

  int nbFoundModes = eigen_vec.GetN();
  betas.Reallocate(nbFoundModes);
  E3.Reallocate(nbFoundModes);
  H3.Reallocate(nbFoundModes);

  for (int k = 0; k < nbFoundModes; k++)
    {
      int nb_dof = N;
      VectComplex_wp x(2*nb_dof), X;
      VectComplex_wp y;
      GetCol(eigen_vec, k, X);

      VectComplex_wp Es, Hs, E1, E2, H1, H2;
      var_eig.ComputeEsHs(lambda(k), X, Es, Hs, E1, E2, H1, H2);

      // X.Resize(nb_dof_E3 + nb_dof_H3);
      // VectComplex_wp Y(nb_dof_E3 + nb_dof_H3); Y.Zero();
      // ModeEsNonLinear2NDirichletEigenProblem var_eig_model(N, *this,
      //                                     this->DofKeptDir, this->IndexDirichlet,
      //                                     false);
      //
      // var_eig_model.MltOperator(lambda(k), SeldonNoTrans, X, Y);
      cout << "VP " << k << ": " << lambda(k) << endl;
      // cout << "VP " << k << ": " << lambda(k) << " ; ||T(VP) X|| : " << Norm2(Y) << " ; erreur relative : " << Norm2(Y) / Norm2(X) << endl;

      E3(k).Reallocate(nb_dof); E3(k).Zero();
      H3(k).Reallocate(nb_dof); H3(k).Zero();
      for (int i = 0; i < nb_dof; i++)
        {
          E3(k).Get(i) = X(i);
          H3(k).Get(i) = X(nb_dof + i);
        }
      betas.Get(k) = lambda(k);
      WriteOutputFileScalar(E3(k), "ModeE3_"+ to_str(k));
      WriteOutputFileScalar(H3(k), "ModeH3_"+ to_str(k));
      WriteOutputFileScalar(Es, "ModeEs_"+ to_str(k));
    }
}
