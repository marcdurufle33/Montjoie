void ModeEz_Solver::ResolutionNonLinearFullDirichlet(Vector<Complex_wp>& betas,
                                  Vector<VectComplex_wp>& Ez, Vector<VectComplex_wp>& Hz,
                                  Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                                  Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                                  Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch)
{
  int N = this->mesh_num.GetNbDof();

  ModeEzNonLinearFullDirichletEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch,
                                        this->DofKeptNonDir, this->IndexNonDirichlet,
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

  param.SetIntervalRegion(0.1, 10.0, 0.001, 10.0);
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
  Ez.Reallocate(nbFoundModes);
  Hz.Reallocate(nbFoundModes);

  for (int k = 0; k < nbFoundModes; k++)
    {
      int nb_dof_mesh = N;
      VectComplex_wp x(2*nb_dof_mesh), X;
      VectComplex_wp y;
      GetCol(eigen_vec, k, X);
      var_eig.ExtractDirichlet(X, x);
      Ez(k).Reallocate(nb_dof_mesh); Ez(k).Zero();
      Hz(k).Reallocate(nb_dof_mesh); Hz(k).Zero();
      for (int i = 0; i < nb_dof_mesh; i++)
        {
          Ez(k).Get(i) = x(i);
          if (!only_Ez)
            {
              Hz(k).Get(i) = x(nb_dof_mesh + i);
            }
        }
      betas.Get(k) = lambda(k);
    }
}

void ModeEz_Solver::ResolutionNonLinearDirichlet(Vector<Complex_wp>& betas,
                                  Vector<VectComplex_wp>& Ez, Vector<VectComplex_wp>& Hz,
                                  Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                                  Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                                  Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch)
{
  int N = this->mesh_num.GetNbDof();

  ModeEzNonLinearDirichletEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch,
                                        this->DofKeptNonDir, this->IndexNonDirichlet,
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

  param.SetIntervalRegion(0.1, 10.0, 0.001, 10.0);
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
  Ez.Reallocate(nbFoundModes);
  Hz.Reallocate(nbFoundModes);

  for (int k = 0; k < nbFoundModes; k++)
    {
      int nb_dof_Hz = N;
      VectComplex_wp x(2*nb_dof_Hz), X;
      VectComplex_wp y;
      GetCol(eigen_vec, k, X);
      var_eig.ExtractDirichlet(X, x);
      Ez(k).Reallocate(nb_dof_Hz); Ez(k).Zero();
      Hz(k).Reallocate(nb_dof_Hz); Hz(k).Zero();
      for (int i = 0; i < nb_dof_Hz; i++)
        {
          Ez(k).Get(i) = x(i);
          if (!only_Ez)
            {
              Hz(k).Get(i) = x(nb_dof_Hz + i);
            }
        }
      betas.Get(k) = lambda(k);
    }
}

void ModeEz_Solver::ResolutionNonLinearBGT(Vector<Complex_wp>& betas,
                                  Vector<VectComplex_wp>& Ez, Vector<VectComplex_wp>& Hz,
                                  Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                                  Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                                  Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch,
                                  VectReal_wp& Mh_surf,
                                  Matrix<Complexe, Prop, Storage>& Kh_surf,
                                  Matrix<Complexe, Prop, Storage>& Rh_surf)
{
  int N = this->mesh_num.GetNbDof();

  ModeEzNonLinearBGTEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch,
                                          Kh_surf, Mh_surf, Rh_surf, this->use_split_formulation);

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

  param.SetIntervalRegion(0.1, 10.0, 0.001, 10.0);
  param.EnableCommandLineOptions();
  if (this->use_default_petsc_solver)
    param.SetDefaultPetscSolver();

  cout << "On va commencer !!!" << endl;

  Vector<Complex_wp> lambda, lambda_imag;
  Matrix<Complex_wp, General, ColMajor> eigen_vec;

  DISP(var_eig.GetM());
  FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);

  DISP(lambda); DISP(lambda_imag);

  int nbFoundModes = eigen_vec.GetN();
  betas.Reallocate(nbFoundModes);
  Ez.Reallocate(nbFoundModes);
  Hz.Reallocate(nbFoundModes);

  for (int k = 0; k < nbFoundModes; k++)
    {
      int nb_dof = N;
      VectComplex_wp X;
      VectComplex_wp Y(2*nb_dof); Y.Zero();
      GetCol(eigen_vec, k, X);

      var_eig.MltOperator(lambda(k), SeldonNoTrans, X, Y);
      cout << "VP " << k << ": " << lambda(k) << " ; norme : " << Norm2(Y) / Norm2(X) << endl;

      Ez(k).Reallocate(nb_dof); Ez(k).Zero();
      Hz(k).Reallocate(nb_dof); Hz(k).Zero();
      for (int i = 0; i < nb_dof; i++)
        {
          Ez(k).Get(i) = X(i);
          if (!only_Ez)
            {
              Hz(k).Get(i) = X(nb_dof + i);
            }
        }
      betas.Get(k) = lambda(k);
    }
}
