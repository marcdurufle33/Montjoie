// initialisation du probleme
ModeEsPolynomialDirichletEigenProblem
::ModeEsPolynomialDirichletEigenProblem(int N, ModeEs_Solver& var0,
                              Vector<Matrix<Complexe, Prop, Storage> >& Mh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ch,
                              Vector<Matrix<Complexe, Prop, Storage> >& Dh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Th,
                              Vector<Matrix<Complexe, Prop, Storage> >& Uh,
                              Vector<Matrix<Complexe, Symmetric, DiagonalRow> >& DTau0,
                              Vector<Matrix<Complexe, Symmetric, DiagonalRow> >& ktr2,
                              Vector<int>& DofDir, Vector<int>& IndexDir)
  :
  var(var0), vec_Mh(Mh), vec_Kh(Kh), vec_Ch(Ch), vec_Dh(Dh), vec_Th(Th), vec_Uh(Uh), vec_D_tau0(DTau0), vec_ktrSquared(ktr2), DofKeptDir(DofDir), IndexDirichlet(IndexDir)
{
    // pour Es, condition de Dirichlet
    nb_dof_Es = DofKeptDir.GetM();

    // pour Hs, condition de Neumann
    nb_dof_Hs = N;

    // taille du pb aux valeurs propres
    nb_rows = nb_dof_Es;
    nb_dof_all = N;
    nb_rows += nb_dof_Hs;
    nb_dof_all += N;

    vec_twisted.Reallocate(var.ref_tau.GetM());
    vec_twisted.Zero();
    for (int ref = 0 ; ref < var.ref_tau.GetM() ; ref++)
      {
        vec_twisted(ref) = (var.ref_tau(ref) != 0.0);
      }
    DISP(vec_twisted);

    adim_coef = 10000.0;
    DISP(adim_coef);

    this->Init(nb_rows);
    DISP(nb_rows); DISP(nb_dof_Es); DISP(nb_dof_Hs);

    ComputeMatrix();
}

void ModeEsPolynomialDirichletEigenProblem::RemoveDirichletRowCol(const Matrix<Complexe, Prop, Storage>& A,
                                                        Matrix<Complexe, Prop, Storage>& B)
{
  B.Clear();
  B.Reallocate(nb_dof_Es, nb_dof_Es);
  for (int i = 0; i < nb_dof_Es; i++)
    {
      int i0 = DofKeptDir(i);
      int nb_val = 0;
      for (int j = 0; j < A.GetRowSize(i0); j++)
        if (IndexDirichlet(A.Index(i0, j)) >= 0)
          nb_val++;

      B.ReallocateRow(i, nb_val);
      nb_val = 0;
      for (int j = 0; j < A.GetRowSize(i0); j++)
        if (IndexDirichlet(A.Index(i0, j)) >= 0)
          {
            B.Index(i, nb_val) = IndexDirichlet(A.Index(i0, j));
            B.Value(i, nb_val) = A.Value(i0, j);
            nb_val++;
          }
    }
}

void ModeEsPolynomialDirichletEigenProblem::ExtractDirichlet(const VectComplex_wp& Xdir, VectComplex_wp& X)
{
  X.Zero();
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    X(DofKeptDir(i)) = Xdir(i);

  for (int i = 0; i < nb_dof_Hs; i++)
    X(nb_dof_Hs + i) = Xdir(nb_dof_Es + i);
}

void ModeEsPolynomialDirichletEigenProblem::ExpandDirichlet(const VectComplex_wp& X, VectComplex_wp& Xdir)
{
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    Xdir(i) = X(DofKeptDir(i));

  for (int i = 0; i < nb_dof_Hs; i++)
    Xdir(nb_dof_Es + i) = X(nb_dof_Hs + i);
}


void ModeEsPolynomialDirichletEigenProblem::ComputeInvDBetaMatrix(const Complex_wp& L)
{
  int N = nb_dof_Hs;
  Complex_wp beta = L*var.omega;
  for (int ref = 0; ref < vec_D_tau0.GetM(); ref++)
    if (vec_D_tau0(ref).GetM() > 0)
      {
        for (int i = 0 ; i < N ; i++)
        {
          vec_Inv_D_tauBeta(ref).Get(i,i) = 1.0/(vec_D_tau0(ref)(i,i) + beta * beta);
        }
      }
}

void ModeEsPolynomialDirichletEigenProblem::ComputeA0()
{
  A0.Clear();
  A0.Reallocate(nb_rows, nb_rows);

  // Real_wp contrib_Eigenvalue = 1.0;
  Complex_wp contrib_Eigenvalue = adim_coef * 1.0 / (var.omega * var.omega);
  // Real_wp contrib_Eigenvalue = 1.0 / (var.omega * var.omega * var.omega * var.omega);

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    {
      if (vec_Kh(ref).GetM() > 0)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = - eps_sigma * mu_tilde;
          DISP(ref);
          DISP(eps_sigma);
          DISP(mu_tilde);
          DISP(k2);
          Real_wp tau = var.ref_tau(ref);

          // on traite les lignes avec V2
          for (int iDir = 0; iDir < nb_dof_Es; iDir++)
            {
              int i = DofKeptDir(iDir); // i = ligne globale
              int nM_Hs = vec_Mh(ref).GetRowSize(i);
              int nK_Hs = vec_Kh(ref).GetRowSize(i);
              int nD_Hs = vec_Dh(ref).GetRowSize(i);
              int nT_Hs = vec_Th(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;
              int num = 0;
              col.Reallocate(4*nM_Hs); val.Reallocate(4*nM_Hs);
              // parties avec Uh
              for (int j = 0; j < nM_Hs; j++)
                {
                  int j_Glob = vec_Mh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                  if (jDir >= 0)
                    {
                      col(num) = jDir;
                      // partie - eps_sigma Mh V2t
                      val(num) = eps_sigma * k2 * k2 * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                      num++;

                      if (vec_twisted(ref))
                        {
                          col(num) = jDir;
                          // partie - eps_sigma Mh V2t
                          val(num) = eps_sigma * k2 * contrib_Eigenvalue * vec_Mh(ref).Value(i, j) * ktr2;
                          num++;

                          col(num) = jDir;
                          // partie - eps_sigma Mh V2t
                          val(num) = 4.0 * eps_sigma * k2 * tau * tau * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                          num++;
                        }
                    }

                  if (vec_twisted(ref))
                    {
                      col(num) = nb_dof_Es + j_Glob;
                      // partie 2 k2 tau Mh W2
                      val(num) = - 2.0 * k2 * k2 * tau * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                      num++;
                    }
                }

              A0.AddInteractionRow(iDir, num, col, val);

              num = 0;
              col.Reallocate(2 * nK_Hs); val.Reallocate(2 * nK_Hs);
              // partie avec Kh
              for (int j = 0; j < nK_Hs; j++)
                {
                  int j_Glob = vec_Kh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                  if (jDir >= 0)
                    {
                      col(num) = jDir;
                      // partie eps_sigma Kh V2
                      val(num) = - eps_sigma * k2 * contrib_Eigenvalue * vec_Kh(ref).Value(i, j);
                      num++;

                      if (vec_twisted(ref))
                        {
                          col(num) = jDir;
                          // partie eps_sigma Kh V2
                          val(num) = - eps_sigma * contrib_Eigenvalue * vec_Kh(ref).Value(i, j) * ktr2;
                          num++;
                        }
                    }
                }

              A0.AddInteractionRow(iDir, num, col, val);


              if (vec_twisted(ref))
                {
                  num = 0;
                  col.Reallocate(2 * nD_Hs); val.Reallocate(2 * nD_Hs);
                  // partie avec Dh
                  for (int j = 0; j < nD_Hs; j++)
                    {
                      int j_Glob = vec_Dh(ref).Index(i, j);
                      int jDir = IndexDirichlet(j_Glob);
                      Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                      if (jDir >= 0)
                        {
                          col(num) = jDir;
                          // partie eps_sigma tau^2 Dh V2
                          val(num) = - eps_sigma * k2 * tau * tau * contrib_Eigenvalue * vec_Dh(ref).Value(i, j);
                          num++;

                          col(num) = jDir;
                          // partie eps_sigma tau^2 Dh V2
                          val(num) = - eps_sigma * tau * tau * contrib_Eigenvalue * vec_Dh(ref).Value(i, j) * ktr2;
                          num++;
                        }
                    }

                  A0.AddInteractionRow(iDir, num, col, val);

                  num = 0;
                  col.Reallocate(nT_Hs); val.Reallocate(nT_Hs);
                  // partie avec Dh
                  for (int j = 0; j < nT_Hs; j++)
                    {
                      int j_Glob = vec_Th(ref).Index(i, j);
                      int jDir = IndexDirichlet(j_Glob);
                      if (jDir >= 0)
                        {
                          col(num) = jDir;
                          // partie eps_sigma tau^2 Dh V2
                          val(num) = 2.0 * eps_sigma * tau * tau * k2 * contrib_Eigenvalue * vec_Th(ref).Value(i, j);
                          num++;
                        }
                    }

                  A0.AddInteractionRow(iDir, num, col, val);
                }
            }

          // on traite les lignes avec W2
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              int nM_Hs = vec_Mh(ref).GetRowSize(i);
              int nK_Hs = vec_Kh(ref).GetRowSize(i);
              int nD_Hs = vec_Dh(ref).GetRowSize(i);
              int nT_Hs = vec_Th(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;

              int num = 0;
              col.Reallocate(4*nM_Hs); val.Reallocate(4*nM_Hs);
              // parties avec Uh
              for (int j = 0; j < nM_Hs; j++)
                {
                  int j_Glob = vec_Mh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);

                  col(num) = nb_dof_Es + j_Glob;
                  // partie - mu_tilde Mh W2t
                  val(num) = mu_tilde * k2 * k2 * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                  num++;

                  if (vec_twisted(ref))
                    {
                      col(num) = nb_dof_Es + j_Glob;
                      // partie - mu_tilde Mh W2t
                      val(num) = mu_tilde * k2 * contrib_Eigenvalue * vec_Mh(ref).Value(i, j) * ktr2;
                      num++;

                      col(num) = nb_dof_Es + j_Glob;
                      // partie - mu_tilde Mh W2t
                      val(num) = 4.0 * k2 * tau * tau * mu_tilde * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                      num++;
                    }

                  if ((jDir >= 0) && vec_twisted(ref))
                    {
                      col(num) = jDir;
                      // partie - 2 k2 tau Mh, V2
                      val(num) = 2.0 * k2 * k2 * tau * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                      num++;
                    }
                }

              A0.AddInteractionRow(nb_dof_Es + i, num, col, val);

              num = 0;
              col.Reallocate(2*nK_Hs); val.Reallocate(2*nK_Hs);
              // partie avec Kh
              for (int j = 0; j < nK_Hs; j++)
                {
                  int j_Glob = vec_Kh(ref).Index(i, j);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                  // int jDir = IndexDirichlet(j_Glob);
                  // Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob, j_Glob);
                  col(num) = nb_dof_Es + j_Glob;
                  // partie mu_tilde Kh W2
                  val(num) = - mu_tilde * k2 * contrib_Eigenvalue * vec_Kh(ref).Value(i, j);
                  num++;

                  if (vec_twisted(ref))
                    {
                      col(num) = nb_dof_Es + j_Glob;
                      // partie mu_tilde Kh W2
                      val(num) = - mu_tilde * contrib_Eigenvalue * vec_Kh(ref).Value(i, j) * ktr2;
                      num++;
                    }
                }

              A0.AddInteractionRow(nb_dof_Es + i, num, col, val);

              if (vec_twisted(ref))
                {
                  num = 0;
                  col.Reallocate(2 * nD_Hs); val.Reallocate(2 * nD_Hs);
                  // partie avec Dh
                  for (int j = 0; j < nD_Hs; j++)
                    {
                      int j_Glob = vec_Dh(ref).Index(i, j);
                      Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                      // int jDir = IndexDirichlet(j_Glob);
                      // Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob, j_Glob);
                      col(num) = nb_dof_Es + j_Glob;
                      // partie mu_tilde tau^2 Dh W2
                      val(num) = - mu_tilde * k2 * tau * tau * contrib_Eigenvalue * vec_Dh(ref).Value(i, j);
                      num++;

                      col(num) = nb_dof_Es + j_Glob;
                      // partie mu_tilde tau^2 Dh W2
                      val(num) = - mu_tilde * tau * tau * contrib_Eigenvalue * vec_Dh(ref).Value(i, j) * ktr2;
                      num++;
                    }

                  A0.AddInteractionRow(nb_dof_Es + i, num, col, val);


                  num = 0;
                  col.Reallocate(nT_Hs); val.Reallocate(nT_Hs);
                  // partie avec Dh
                  for (int j = 0; j < nT_Hs; j++)
                    {
                      int j_Glob = vec_Th(ref).Index(i, j);
                      // int jDir = IndexDirichlet(j_Glob);
                      // Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob, j_Glob);
                      col(num) = nb_dof_Es + j_Glob;
                      // partie mu_tilde tau^2 Dh W2
                      val(num) = 2.0 * k2 * mu_tilde * tau * tau * contrib_Eigenvalue * vec_Th(ref).Value(i, j);
                      num++;
                    }

                  A0.AddInteractionRow(nb_dof_Es + i, num, col, val);
                }
            }
        }
    }
}

void ModeEsPolynomialDirichletEigenProblem::ComputeA1()
{
  A1.Clear();
  A1.Reallocate(nb_rows, nb_rows);

  // Real_wp contrib_Eigenvalue = var.omega;
  Complex_wp contrib_Eigenvalue = adim_coef * 1.0 / var.omega;
  // Real_wp contrib_Eigenvalue = 1.0 / (var.omega * var.omega * var.omega);

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    {
      if (vec_Kh(ref).GetM() > 0)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = - eps_sigma * mu_tilde;
          Real_wp tau = var.ref_tau(ref);

          // on traite les lignes avec V2
          for (int iDir = 0; iDir < nb_dof_Es; iDir++)
            {
              int i = DofKeptDir(iDir); // i = ligne globale
              int nC_Hs = vec_Ch(ref).GetRowSize(i);
              int nU_Hs = vec_Uh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;
              int num = 0;
              col.Reallocate(2*nC_Hs); val.Reallocate(2*nC_Hs);
              // partie avec Ch
              for (int j = 0; j < nC_Hs; j++)
                {
                  int j_Glob = vec_Ch(ref).Index(i, j);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);

                  col(num) = nb_dof_Es + j_Glob;
                  // partie - Iwp Ch W2t
                  val(num) = - Iwp * k2 * contrib_Eigenvalue * vec_Ch(ref).Value(i, j);
                  num++;

                  if (vec_twisted(ref))
                    {
                      col(num) = nb_dof_Es + j_Glob;
                      // partie - Iwp Ch W2t
                      val(num) = - Iwp * contrib_Eigenvalue * vec_Ch(ref).Value(i, j) * ktr2;
                      num++;
                    }
                }

              A1.AddInteractionRow(iDir, num, col, val);

              if (vec_twisted(ref))
                {
                  num = 0;
                  col.Reallocate(3*nU_Hs); val.Reallocate(3*nU_Hs);
                  // parties avec Uh
                  for (int j = 0; j < nU_Hs; j++)
                    {
                      int j_Glob = vec_Uh(ref).Index(i, j);
                      int jDir = IndexDirichlet(j_Glob);
                      Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);

                      if (jDir >= 0)
                        {
                          col(num) = jDir;
                          // partie 2 i tau eps_sigma Uh V2t
                          val(num) = 2.0 * Iwp * tau * eps_sigma * contrib_Eigenvalue * vec_Uh(ref).Value(i, j) * ktr2;
                          num++;

                          col(num) = jDir;
                          // partie 2 i tau eps_sigma Uh V2t
                          val(num) = 2.0 * Iwp * tau * eps_sigma * k2 * contrib_Eigenvalue * vec_Uh(ref).Value(i, j);
                          num++;
                        }

                      col(num) = nb_dof_Es + j_Glob;
                      // partie -2 Iwp k2 tau^2 Uh W2
                      val(num) = - 2.0 * Iwp * k2 * tau * tau * contrib_Eigenvalue * vec_Uh(ref).Value(i, j);
                      num++;
                    }
                  A1.AddInteractionRow(iDir, num, col, val);
                }
            }

          // on traite les lignes avec W2
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              int nC_Hs = vec_Ch(ref).GetRowSize(i);
              int nU_Hs = vec_Uh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;

              int num = 0;
              col.Reallocate(2*nC_Hs); val.Reallocate(2*nC_Hs);
              // partie avec Ch
              for (int j = 0; j < nC_Hs; j++)
                {
                  int j_Glob = vec_Ch(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                  if (j >= 0)
                    {
                      if (vec_twisted(ref))
                        {
                          col(num) = jDir;
                          // partie Iwp Ch V2t
                          val(num) = Iwp * contrib_Eigenvalue * vec_Ch(ref).Value(i, j) * ktr2;
                          num++;
                        }
                      col(num) = jDir;
                      // partie Iwp Ch V2t
                      val(num) = Iwp * k2 * contrib_Eigenvalue * vec_Ch(ref).Value(i, j);
                      num++;
                    }
                }

              A1.AddInteractionRow(nb_dof_Es + i, num, col, val);

              if (vec_twisted(ref))
                {
                  num = 0;
                  col.Reallocate(3*nU_Hs); val.Reallocate(3*nU_Hs);
                  // parties avec Uh
                  for (int j = 0; j < nU_Hs; j++)
                    {
                      int j_Glob = vec_Uh(ref).Index(i, j);
                      int jDir = IndexDirichlet(j_Glob);

                      Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);

                      col(num) = nb_dof_Es + j_Glob;
                      // partie 2 i mu_tilde tau Uh W2t
                      val(num) = 2.0 * Iwp * tau * mu_tilde * contrib_Eigenvalue * vec_Uh(ref).Value(i, j) * ktr2;
                      num++;

                      col(num) = nb_dof_Es + j_Glob;
                      // partie 2 i mu_tilde tau Uh W2t
                      val(num) = 2.0 * Iwp * tau * mu_tilde * k2 * contrib_Eigenvalue * vec_Uh(ref).Value(i, j);
                      num++;

                      if (jDir >= 0)
                        {
                          col(num) = jDir;
                          // partie 2 i k2 tau^2 Uh V2
                          val(num) = 2.0 * Iwp * k2 * tau * tau * contrib_Eigenvalue * vec_Uh(ref).Value(i, j);
                          num++;
                        }
                    }

                  A1.AddInteractionRow(nb_dof_Es + i, num, col, val);
              }
            }
        }
    }
}

void ModeEsPolynomialDirichletEigenProblem::ComputeA2()
{
  A2.Clear();
  A2.Reallocate(nb_rows, nb_rows);

  // Real_wp contrib_Eigenvalue = var.omega * var.omega;
  Complex_wp contrib_Eigenvalue = adim_coef * 1.0;
  // Real_wp contrib_Eigenvalue = 1.0 / (var.omega * var.omega);

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    {
      if (vec_Kh(ref).GetM() > 0)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = - eps_sigma * mu_tilde;
          Real_wp tau = var.ref_tau(ref);


          // on traite les lignes avec V2
          for (int iDir = 0; iDir < nb_dof_Es; iDir++)
            {
              int i = DofKeptDir(iDir); // i = ligne globale
              int nM_Hs = vec_Mh(ref).GetRowSize(i);
              int nK_Hs = vec_Kh(ref).GetRowSize(i);
              int nD_Hs = vec_Dh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;
              int num = 0;
              col.Reallocate(3*nM_Hs); val.Reallocate(3*nM_Hs);
              // parties avec Uh
              for (int j = 0; j < nM_Hs; j++)
                {
                  int j_Glob = vec_Mh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                  if (jDir >= 0)
                    {
                      col(num) = jDir;
                      // partie - 2 eps_sigma k^2 Mh V2
                      val(num) = - 2.0 * eps_sigma * k2 * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                      num++;

                      if (vec_twisted(ref))
                        {
                          col(num) = jDir;
                          // partie - eps_sigma Mh V2t
                          val(num) = - eps_sigma * contrib_Eigenvalue * vec_Mh(ref).Value(i, j) * ktr2;
                          num++;
                        }
                    }

                  if (vec_twisted(ref))
                    {
                      col(num) = nb_dof_Es + j_Glob;
                      // partie 2 k2 tau Mh W2
                      val(num) = 2.0 * k2 * tau * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                      num++;
                    }
                }

              A2.AddInteractionRow(iDir, num, col, val);

              num = 0;
              col.Reallocate(nK_Hs); val.Reallocate(nK_Hs);
              // partie avec Kh
              for (int j = 0; j < nK_Hs; j++)
                {
                  int j_Glob = vec_Kh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  if (jDir >= 0)
                    {
                      col(num) = jDir;
                      // partie eps_sigma Kh V2
                      val(num) = eps_sigma * contrib_Eigenvalue * vec_Kh(ref).Value(i, j);
                      num++;
                    }
                }

              A2.AddInteractionRow(iDir, num, col, val);
              if (vec_twisted(ref))
                {
                  num = 0;
                  col.Reallocate(nD_Hs); val.Reallocate(nD_Hs);
                  // partie avec Dh
                  for (int j = 0; j < nD_Hs; j++)
                    {
                      int j_Glob = vec_Dh(ref).Index(i, j);
                      int jDir = IndexDirichlet(j_Glob);
                      if (jDir >= 0)
                        {
                          col(num) = jDir;
                          // partie eps_sigma tau^2 Dh V2
                          val(num) = eps_sigma * tau * tau * contrib_Eigenvalue * vec_Dh(ref).Value(i, j);
                          num++;
                        }
                    }

                  A2.AddInteractionRow(iDir, num, col, val);
                }
            }

          // on traite les lignes avec W2
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              int nM_Hs = vec_Mh(ref).GetRowSize(i);
              int nK_Hs = vec_Kh(ref).GetRowSize(i);
              int nD_Hs = vec_Dh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;

              int num = 0;
              col.Reallocate(3*nM_Hs); val.Reallocate(3*nM_Hs);
              // parties avec Uh
              for (int j = 0; j < nM_Hs; j++)
                {
                  int j_Glob = vec_Mh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);

                  col(num) = nb_dof_Es + j_Glob;
                  // partie - mu_tilde Mh W2t
                  val(num) = - 2.0 * mu_tilde * k2 * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                  num++;

                  if (vec_twisted(ref))
                    {
                      col(num) = nb_dof_Es + j_Glob;
                      // partie - mu_tilde Mh W2t
                      val(num) = - mu_tilde * contrib_Eigenvalue * vec_Mh(ref).Value(i, j) * ktr2;
                      num++;
                    }

                  if ((jDir >= 0) && vec_twisted(ref))
                    {
                      col(num) = jDir;
                      // partie - 2 k2 tau Mh, V2
                      val(num) = - 2.0 * k2 * tau * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                      num++;
                    }
                }

              A2.AddInteractionRow(nb_dof_Es + i, num, col, val);

              num = 0;
              col.Reallocate(nK_Hs); val.Reallocate(nK_Hs);
              // partie avec Kh
              for (int j = 0; j < nK_Hs; j++)
                {
                  int j_Glob = vec_Kh(ref).Index(i, j);
                  // int jDir = IndexDirichlet(j_Glob);
                  // Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob, j_Glob);
                  col(num) = nb_dof_Es + j_Glob;
                  // partie mu_tilde Kh W2
                  val(num) = mu_tilde * contrib_Eigenvalue * vec_Kh(ref).Value(i, j);
                  num++;
                }

              A2.AddInteractionRow(nb_dof_Es + i, num, col, val);

              if (vec_twisted(ref))
                {
                  num = 0;
                  col.Reallocate(nD_Hs); val.Reallocate(nD_Hs);
                  // partie avec Dh
                  for (int j = 0; j < nD_Hs; j++)
                    {
                      int j_Glob = vec_Dh(ref).Index(i, j);
                      // int jDir = IndexDirichlet(j_Glob);
                      // Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob, j_Glob);
                      col(num) = nb_dof_Es + j_Glob;
                      // partie mu_tilde tau^2 Dh W2
                      val(num) = mu_tilde * tau * tau * contrib_Eigenvalue * vec_Dh(ref).Value(i, j);
                      num++;
                    }

                  A2.AddInteractionRow(nb_dof_Es + i, num, col, val);
                }
            }
        }
    }
}

void ModeEsPolynomialDirichletEigenProblem::ComputeA3()
{
  A3.Clear();
  A3.Reallocate(nb_rows, nb_rows);

  // Real_wp contrib_Eigenvalue = var.omega * var.omega * var.omega;
  Complex_wp contrib_Eigenvalue = adim_coef * var.omega;
  // Real_wp contrib_Eigenvalue = 1.0 / var.omega;

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    {
      if (vec_Kh(ref).GetM() > 0)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = - eps_sigma * mu_tilde;
          Real_wp tau = var.ref_tau(ref);

          // on traite les lignes avec V2
          for (int iDir = 0; iDir < nb_dof_Es; iDir++)
            {
              int i = DofKeptDir(iDir); // i = ligne globale
              int nC_Hs = vec_Ch(ref).GetRowSize(i);
              int nU_Hs = vec_Uh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;
              int num = 0;
              col.Reallocate(nC_Hs); val.Reallocate(nC_Hs);
              // partie avec Ch
              for (int j = 0; j < nC_Hs; j++)
                {
                  int j_Glob = vec_Ch(ref).Index(i, j);
                  col(num) = nb_dof_Es + j_Glob;
                  // partie Iwp Ch W2
                  val(num) = Iwp * contrib_Eigenvalue * vec_Ch(ref).Value(i, j);
                  num++;
                }

              A3.AddInteractionRow(iDir, num, col, val);

              if (vec_twisted(ref))
                {
                  num = 0;
                  col.Reallocate(nU_Hs); val.Reallocate(nU_Hs);
                  // parties avec Uh
                  for (int j = 0; j < nU_Hs; j++)
                    {
                      int j_Glob = vec_Uh(ref).Index(i, j);
                      int jDir = IndexDirichlet(j_Glob);
                      if (jDir >= 0)
                        {
                          col(num) = jDir;
                          // partie -2 i tau eps_sigma Uh V2
                          val(num) = -2.0 * Iwp * tau * eps_sigma * contrib_Eigenvalue * vec_Uh(ref).Value(i, j);
                          num++;
                        }
                    }

                  A3.AddInteractionRow(iDir, num, col, val);
                }
            }

          // on traite les lignes avec W2
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              int nC_Hs = vec_Ch(ref).GetRowSize(i);
              int nU_Hs = vec_Uh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;

              int num = 0;
              col.Reallocate(nC_Hs); val.Reallocate(nC_Hs);
              // partie avec Ch
              for (int j = 0; j < nC_Hs; j++)
                {
                  int j_Glob = vec_Ch(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  if (j >= 0)
                    {
                      col(num) = jDir;
                      // partie - i Ch V2
                      val(num) = - Iwp * contrib_Eigenvalue * vec_Ch(ref).Value(i, j);
                      num++;
                    }
                }

              A3.AddInteractionRow(nb_dof_Es + i, num, col, val);

              if (vec_twisted(ref))
                {
                  num = 0;
                  col.Reallocate(nU_Hs); val.Reallocate(nU_Hs);
                  // parties avec Uh
                  for (int j = 0; j < nU_Hs; j++)
                    {
                      int j_Glob = vec_Uh(ref).Index(i, j);
                      int jDir = IndexDirichlet(j_Glob);

                      col(num) = nb_dof_Es + j_Glob;
                      // partie - 2 i mu_tilde tau Uh W2
                      val(num) = - 2.0 * Iwp * tau * mu_tilde * contrib_Eigenvalue * vec_Uh(ref).Value(i, j);
                      num++;
                    }

                  A3.AddInteractionRow(nb_dof_Es + i, num, col, val);
                }
            }

        }
    }
}

void ModeEsPolynomialDirichletEigenProblem::ComputeA4()
{
  A4.Clear();
  A4.Reallocate(nb_rows, nb_rows);

  // Real_wp contrib_Eigenvalue = var.omega * var.omega * var.omega * var.omega;
  Complex_wp contrib_Eigenvalue = adim_coef * var.omega * var.omega;
  // Real_wp contrib_Eigenvalue = 1.0;

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    {
      if (vec_Kh(ref).GetM() > 0)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);

          // on traite les lignes avec V2
          for (int iDir = 0; iDir < nb_dof_Es; iDir++)
            {
              int i = DofKeptDir(iDir); // i = ligne globale
              int nM_Hs = vec_Mh(ref).GetRowSize(i);
              Vector<int> col; Vector<Complex_wp> val;
              int num = 0;
              col.Reallocate(nM_Hs); val.Reallocate(nM_Hs);
              // parties avec Mh
              for (int j = 0; j < nM_Hs; j++)
                {
                  int j_Glob = vec_Mh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  if (jDir >= 0)
                    {
                      col(num) = jDir;
                      // partie eps_sigma M_h V2
                      val(num) = eps_sigma * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                      num++;
                    }
                }
              A4.AddInteractionRow(iDir, num, col, val);
            }

          // on traite les lignes avec W2
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              int nM_Hs = vec_Mh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;

              int num = 0;
              col.Reallocate(nM_Hs); val.Reallocate(nM_Hs);
              // parties avec Uh
              for (int j = 0; j < nM_Hs; j++)
                {
                  int j_Glob = vec_Mh(ref).Index(i, j);

                  col(num) = nb_dof_Es + j_Glob;
                  // partie mu_tilde Mh W2
                  val(num) = mu_tilde * contrib_Eigenvalue * vec_Mh(ref).Value(i, j);
                  num++;
                }

              A4.AddInteractionRow(nb_dof_Es + i, num, col, val);
            }

        }
    }
}

void ModeEsPolynomialDirichletEigenProblem::ComputeMatrix()
{
  ComputeA0();
  ComputeA1();
  ComputeA2();
  ComputeA3();
  ComputeA4();
}
