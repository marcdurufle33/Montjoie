// initialisation du probleme
ModeEsNonLinear_DD_DirichletEigenProblem
::ModeEsNonLinear_DD_DirichletEigenProblem(int N, ModeEs_Solver& var0,
                              Vector<Matrix<Complex_wp, Prop, Storage> >& Mh,
                              Vector<Matrix<Complex_wp, Prop, Storage> >& Mr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ch,
                              Vector<Matrix<Complexe, Prop, Storage> >& Cr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Rjh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Rjr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Rih,
                              Vector<Matrix<Complexe, Prop, Storage> >& Rir2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Sjh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Sjr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Sih,
                              Vector<Matrix<Complexe, Prop, Storage> >& Sir2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Dh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Dr2h,
                              VectR2& Pos,
                              Vector<int>& DofDir, Vector<int>& IndexDir, bool use_split)
  :
  var(var0), vec_Mh(Mh), vec_Mr2h(Mr2h), vec_Kh(Kh), vec_Kr2h(Kr2h), vec_Ch(Ch), vec_Cr2h(Cr2h),
  vec_Rjh(Rjh), vec_Rjr2h(Rjr2h), vec_Rih(Rih), vec_Rir2h(Rir2h), vec_Sjh(Sjh), vec_Sjr2h(Sjr2h),
  vec_Sih(Sih), vec_Sir2h(Sir2h), vec_Dh(Dh), vec_Dr2h(Dr2h), vec_Pos(Pos), DofKeptDir(DofDir), IndexNonDirichlet(IndexDir)
{
  // pour Es, condition de Dirichlet
  nb_dof_E3 = DofKeptDir.GetM();

  // pour Hs, condition de Neumann
  nb_dof_H3 = N;

  countIteration = 0;

  // taille du pb aux valeurs propres
  nb_rows = nb_dof_E3;
  nb_dof_all = N;
  nb_rows += nb_dof_H3;
  nb_dof_all += N;

  // construction de la matrice R2
  R2.Reallocate(N, N); R2.Zero();
  for (int i = 0 ; i < N ; i++)
    {
      Vector<int> col(1); Vector<Complex_wp> val(1);
      col(0) = i;
      val(0) = vec_Pos(i)(0)*vec_Pos(i)(0) + vec_Pos(i)(1)*vec_Pos(i)(1);
      R2.AddInteractionRow(i, 1, col, val);
    }

  // construction des numerotations
  cout << "N = " << N << endl;
  nb_couches = vec_Mh.GetM();
  nb_dof_couche.Reallocate(nb_couches);
  nb_dof_coucheE3.Reallocate(nb_couches);
  nb_dof_couche.Zero(); nb_dof_coucheE3.Zero();
  IndexDofE3.Reallocate(nb_couches);
  NumGlobE3.Reallocate(nb_couches);
  IndexDof.Reallocate(nb_couches);
  NumGlob.Reallocate(nb_couches);

  for (int ref = 0; ref < nb_couches; ref++)
    if (vec_Mh(ref).GetM() == N)
      {
        IndexDof(ref).Reallocate(N);
        IndexDof(ref).Fill(-1);

        IndexDofE3(ref).Reallocate(N);
        IndexDofE3(ref).Fill(-1);
      }

  // IndexDof(ref)(j) renvoie le numero local du ddl global j dans la couche ref
  // NumGlob(ref)(i) renvoie le numero global j du ddl local i dans la couche ref
  // IndexDofE3, NumGlobE3 pour l'inconnue E3 (ddl Dirichlets enleves)
  for (int i = 0; i < var.mesh.GetNbElt(); i++)
    {
      int ref = var.mesh.Element(i).GetReference();
      int nb_dof_loc = var.mesh_num.Element(i).GetNbDof();
      for (int j = 0; j < nb_dof_loc; j++)
        {
          int num_dof = var.mesh_num.Element(i).GetNumberDof(j);
          if (IndexDof(ref)(num_dof) == -1) // si c'est un nouveau
            {
              IndexDof(ref)(num_dof) = nb_dof_couche(ref);
              NumGlob(ref).PushBack(num_dof);
              nb_dof_couche(ref)++;
            }

          if ((IndexNonDirichlet(num_dof) >= 0) && (IndexDofE3(ref)(num_dof) == -1)) // si c'est un nouveau
            {
              IndexDofE3(ref)(num_dof) = nb_dof_coucheE3(ref);
              NumGlobE3(ref).PushBack(num_dof);
              nb_dof_coucheE3(ref)++;
            }
        }
    }
  DISP(nb_dof_E3); DISP(nb_dof_H3);

  // twisted = (var.torsion != 0.0);
  // DISP(twisted);
  index_used.Reallocate(var.ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < var.mesh.GetNbElt(); i++)
    index_used(var.mesh.Element(i).GetReference()) = true;

  for (int ref = 0 ; ref < nb_couches ; ref++)
    if (index_used(ref))
      {
        DISP(ref);
        nb_rows += nb_dof_coucheE3(ref) + nb_dof_couche(ref);
      }

  this->Init(nb_rows);
  DISP(nb_rows);

  // on calcule les poles du pb aux valeurs propres
  VectComplex_wp s;
  split_reference.Clear();
  this->numer_pol_split.Clear();
  this->denom_pol_split.Clear();
  if (use_split)
    {
      Vector<Complexe> num, denom;
      // A1_sum.Reallocate(nb_dof_H3, nb_dof_H3);
      // A1_sumH.Reallocate(nb_dof_H3, nb_dof_H3);

      this->numer_pol_split.PushBack(num);
      this->denom_pol_split.PushBack(denom);
    }

  ofstream file_out;
  if (use_split)
    {
      file_out.open("coef.dat");
      file_out.precision(15);
    }

  DistributedMatrix<Complexe, Prop, Storage> A;
  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
  {
    if (vec_Kh(ref).GetM() > 0)
      {
        Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
        Complex_wp mu_tilde = - Iwp * var.omega * var.ref_mu(ref);
        Complex_wp k2 = - eps_sigma*mu_tilde;
        Complex_wp k2_c = (var.ref_epsilon(ref) + Iwp * var.ref_sigma(ref) / var.omega) * var.ref_mu(ref);
        Complex_wp z = sqrt(-k2) / var.omega;
        // Real_wp tau = var.ref_tau(ref);
        Matrix<Complexe, Prop, Storage> temp;
        s.PushBack(z); s.PushBack(-z);

        if (use_split)
          {

            cout << "Remplissage de f_i" << endl;


            Vector<Complex_wp> num(1), denom(3);
            // f_2(beta) = 1/D(beta)
            num(0) = 1.0;
            denom(0) = var.omega*var.omega; denom(1) = 0.0; denom(2) = -k2;
            // f_2(beta) = k2/D(beta)
            // num(0) = 1.0;
            // denom(0) =  1.0 / k2_c; denom(1) = 0.0; denom(2) = -1.0;
            split_reference.PushBack(ref);
            this->numer_pol_split.PushBack(num);
            this->denom_pol_split.PushBack(denom);

            // f_3(beta) = i beta/D(beta)
            num.Reallocate(2);
            num(1) = 0.0; num(0) = Iwp * var.omega;
            this->numer_pol_split.PushBack(num);
            this->denom_pol_split.PushBack(denom);

            // f_4(beta) = i beta k2L/D(beta)^2
            denom.Reallocate(5);
            // denom(4) = k2*k2; denom(3) = 0.0; denom(2) = - 2.0*k2*var.omega*var.omega; denom(1) = 0.0; denom(0) = var.omega*var.omega*var.omega*var.omega;
            // num(0) = Iwp * var.omega * k2;
            // denom(4) = k2; denom(3) = 0.0; denom(2) = - 2.0*var.omega*var.omega; denom(1) = 0.0; denom(0) = ((var.omega*var.omega)/k2)*var.omega*var.omega;
            denom(4) = k2; denom(3) = 0.0; denom(2) = - 2.0*var.omega*var.omega; denom(1) = 0.0; denom(0) = var.omega*var.omega / k2_c;
            num(0) = Iwp * var.omega; num(1) = 0.0;
            this->numer_pol_split.PushBack(num);
            this->denom_pol_split.PushBack(denom);

            // f_5(beta) = k2L/D(beta)^2
            num.Reallocate(1);
            // num(0) = k2;
            num(0) = 1.0;
            // f_5(beta) = k4L/D(beta)^2
            // denom(4) = 1.0; denom(3) = 0.0; denom(2) = - 2.0/k2_c; denom(1) = 0.0; denom(0) = 1.0 / (k2_c * k2_c);
            this->numer_pol_split.PushBack(num);
            this->denom_pol_split.PushBack(denom);
          }
      }
  }

  DISP(split_reference.GetM());
  this->nb_split_matrix = this->numer_pol_split.GetM();
  DISP(this->nb_split_matrix);

  this->SetSingularities(s); //DISP(s);
}

void ModeEsNonLinear_DD_DirichletEigenProblem::RemoveDirichletRowCol(const Matrix<Complexe, Prop, Storage>& A,
                                                        Matrix<Complexe, Prop, Storage>& B)
{
  cout << "RemoveDirichletRowCol" << endl;
  abort();
}

void ModeEsNonLinear_DD_DirichletEigenProblem::ExtractDirichlet(const VectComplex_wp& Xdir, VectComplex_wp& X)
{
  X.Zero();
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    X(DofKeptDir(i)) = Xdir(i);

  for (int i = 0; i < nb_dof_H3; i++)
    X(nb_dof_H3 + i) = Xdir(nb_dof_E3 + i);
}

void ModeEsNonLinear_DD_DirichletEigenProblem::ExpandDirichlet(const VectComplex_wp& X, VectComplex_wp& Xdir)
{
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    Xdir(i) = X(DofKeptDir(i));

  for (int i = 0; i < nb_dof_H3; i++)
    Xdir(nb_dof_E3 + i) = X(nb_dof_H3 + i);
}

void ModeEsNonLinear_DD_DirichletEigenProblem::ExtractDirichletAuxVariable(const VectComplex_wp& EHdir, const int& ref, VectComplex_wp& E3, VectComplex_wp& H3, VectComplex_wp& Fref, VectComplex_wp& Gref)
{
  int N = this->nb_dof_H3;
  Vector<Complex_wp> EH(nb_dof_all);
  ExtractDirichlet(EHdir, EH);
  for (int i = 0 ; i < N ; i++)
    {
      E3(i) = EH(i);
      H3(i) = EH(N + i);
    }

  Fref.Zero();
  Gref.Zero();

  int inc = nb_dof_E3 + nb_dof_H3;
  if (!index_used(ref))
    {
      cout << "La référence : " << ref << " n'existe pas !" << endl;
      abort();
    }
  else
    {
      for (int i = 0 ; i < ref ; i++)
        {
          if (index_used(i))
            inc += nb_dof_coucheE3(i) + nb_dof_couche(i);
        }
    }

  for (int i = 0; i < NumGlobE3(ref).GetM(); i++)
    {
      Fref(NumGlobE3(ref)(i)) = EHdir(inc + i);
    }

  for (int i = 0; i < NumGlob(ref).GetM(); i++)
    {
      Gref(NumGlob(ref)(i)) = EHdir(inc + nb_dof_coucheE3(ref) + i);
    }
}

void ModeEsNonLinear_DD_DirichletEigenProblem::ExpandDirichletAuxVariable(const int& ref, const VectComplex_wp& E3, const VectComplex_wp& H3, const VectComplex_wp& Fref, const VectComplex_wp& Gref, VectComplex_wp& EHdir)
{
  EHdir.Reallocate(nb_rows); EHdir.Zero();
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    EHdir(i) = E3(DofKeptDir(i));
  for (int i = 0; i < nb_dof_H3; i++)
    EHdir(nb_dof_E3 + i) = H3(i);

  int inc = nb_dof_E3 + nb_dof_H3;
  if (!index_used(ref))
    {
      cout << "La référence : " << ref << " n'existe pas !" << endl;
      abort();
    }
  else
    {
      for (int i = 0 ; i < ref ; i++)
        {
          if (index_used(i))
            inc += nb_dof_coucheE3(i) + nb_dof_couche(i);
        }
    }

  for (int i = 0; i < NumGlobE3(ref).GetM(); i++)
    {
      EHdir(inc + i) = Fref(NumGlobE3(ref)(i));
    }

  for (int i = 0; i < NumGlob(ref).GetM(); i++)
    {
      EHdir(inc + nb_dof_coucheE3(ref) + i) = Gref(NumGlob(ref)(i));
    }
}


void ModeEsNonLinear_DD_DirichletEigenProblem::ComputeOperator(const Complex_wp& L)
{

}

void ModeEsNonLinear_DD_DirichletEigenProblem::Mult(const Matrix<Complexe, Prop, Storage>& A, const  Matrix<Complexe, Symmetric, DiagonalRow>& D, Matrix<Complexe, Prop, Storage>& M)
{
  int N = nb_dof_H3;
  M.Reallocate(N,N);
  for (int i = 0 ; i < N ; i++)
  {
    int n = A.GetRowSize(i);
    Vector<int> col(n); Vector<Complex_wp> val(n);
    int num = 0;
    for (int j0 = 0 ; j0 < n ; j0++)
    {
      int j = A.Index(i, j0);
      col(num) = j;
      val(num) = A.Value(i,j0)*D(j,j);
      num++;
    }
    M.AddInteractionRow(i, num, col, val);
  }
}

void ModeEsNonLinear_DD_DirichletEigenProblem::ComputeInvDBetaMatrix(const Complex_wp& L)
{
}

void ModeEsNonLinear_DD_DirichletEigenProblem::ComputeChgmtVariable(const Complex_wp& L)
{
}


void ModeEsNonLinear_DD_DirichletEigenProblem::ComputeEsHs(const Complex_wp& L,
                                                  const Vector<Complex_wp>& eigenvector,
                                                  Vector<Complex_wp>& Es,
                                                  Vector<Complex_wp>& Hs,
                                                  Vector<Complex_wp>& E1,
                                                  Vector<Complex_wp>& E2,
                                                  Vector<Complex_wp>& H1,
                                                  Vector<Complex_wp>& H2)
{
  int N = var.mesh_num.GetNbDof();
  Matrix<Complexe, General, ArrayRowSparse> MM, MMe;
  Complex_wp beta = L * var.omega;
  Complex_wp one(1, 0);

  Vector<Complex_wp> grE3vXi;
  Vector<Complex_wp> grH3sXi;

  Es.Reallocate(N); Es.Zero();
  Hs.Reallocate(N); Hs.Zero();
  // grE3vXi.Reallocate(N); grE3vXi.Zero();
  // grH3sXi.Reallocate(N); grH3sXi.Zero();
  E1.Reallocate(N); E1.Zero();
  E2.Reallocate(N); E2.Zero();
  H1.Reallocate(N); H1.Zero();
  H2.Reallocate(N); H2.Zero();
  MM.Reallocate(N,N); MM.Zero();
  MMe.Reallocate(N,N); MMe.Zero();

  // int inc = nb_dof_E3 + nb_dof_H3;
  for (int ref = 0 ; ref < nb_couches ; ref++)
    if (index_used(ref))
      {
        Vector<Complex_wp> E3(N), H3(N), F(N), G(N);
        ExtractDirichletAuxVariable(eigenvector, ref, E3, H3, F, G);

        Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
        Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
        Complex_wp k2 = - eps_sigma*mu_tilde;
        Complex_wp DBeta = - k2 + beta * beta;
        Complex_wp tau = var.torsion;


        MltAdd(one, vec_Mh(ref), F, one, Es);
        MltAdd(- k2 * tau * tau / DBeta, vec_Mr2h(ref), F, one, Es);
        MltAdd(mu_tilde * tau / DBeta, vec_Rjh(ref), G, one, Es);
        MltAdd(- (mu_tilde / DBeta) * tau * tau * tau * (k2 / DBeta), vec_Rjr2h(ref), G, one, Es);
        MltAdd(- 4.0 * (mu_tilde / DBeta) * tau * tau * tau * (k2 / DBeta), vec_Mr2h(ref), G, one, Es);
        MltAdd(Iwp * beta * tau / DBeta, vec_Sjh(ref), F, one, Es);
        MltAdd(-(Iwp * beta / DBeta) * tau * tau * tau * (k2 / DBeta), vec_Sjr2h(ref), F, one, Es);

        Add(1.0, vec_Mh(ref), MM);
        Add(1.0, vec_Mh(ref), MMe);
      }

    // IL FAUT REALLOUER TOUTE LA COLONNE !!!!
    for (int i = 0; i < N; i++)
      if (IndexNonDirichlet(i) == -1)
        {
          MMe.ReallocateRow(i, 1); // cette ligne ne contient qu'une entrée non nulle
          MMe.Index(i, 0) = i; // la première entrée non nulle est à la colonne (réelle) i (sur la diagonale donc)
          MMe.Value(i, 0) = 1.0; // et elle vaut 1

          Es(i) = 0.0;
        }

    SparseDistributedSolver<Complex_wp> mat_lu;
    mat_lu.Factorize(MMe);
    mat_lu.Solve(Es);

    // mat_lu.Factorize(MM);
    // mat_lu.Solve(Hs);

    // mat_lu.Solve(grE3vXi);
    // mat_lu.Solve(grH3sXi);
    //
    // var.WriteOutputFileScalar(grE3vXi, "GrE3xi");
    // var.WriteOutputFileScalar(grH3sXi, "GrH3xi");
}


// calcul de Y = T(L) X
void ModeEsNonLinear_DD_DirichletEigenProblem::MltOperator(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  cout << "MltOperator : " << ++countIteration << endl;
  // abort();
}

// calcul explicite de T(L) (utilise si on calcule explicitement les matrices)
void ModeEsNonLinear_DD_DirichletEigenProblem
::ComputeOperatorExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "ComputeOperatorExplicit" << endl;
  abort();
}


void ModeEsNonLinear_DD_DirichletEigenProblem::ComputeJacobian(const Complex_wp& L)
{}

// calcul de Y = T'(L) X
void ModeEsNonLinear_DD_DirichletEigenProblem::MltJacobian(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  cout << "MltJacobian" << endl;
  abort();
}

// calcul explicite de T'(L)
void ModeEsNonLinear_DD_DirichletEigenProblem
::ComputeJacobianExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "Pas écrit : ComputeJacobianExplicit" << endl;
  abort();
}

void ModeEsNonLinear_DD_DirichletEigenProblem::ComputeExplicitPreconditioning(DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  mat_lu.Factorize(A);
}

// factorisation de T(L)
void ModeEsNonLinear_DD_DirichletEigenProblem::ComputePreconditioning(const Complex_wp& L)
{
  DistributedMatrix<Complex_wp, Prop, Storage> A;
  ComputeOperatorExplicit(L, A);

  mat_lu.Factorize(A);
}


void ModeEsNonLinear_DD_DirichletEigenProblem::ComputePreconditioning(const VectComplex_wp& L, const VectComplex_wp& coef)
{
  cout << "Pas encore bien écrit : ComputePreconditioning" << endl;
  abort();
}

void ModeEsNonLinear_DD_DirichletEigenProblem::ComputeSplitPreconditioning(const Vector<int>& numL, const VectComplex_wp& coef)
{
  DistributedMatrix<Complex_wp, Prop, Storage> A, B;

  A.Reallocate(nb_rows, nb_rows);
  for (int k = 0; k < numL.GetM(); k++)
    {
      ComputeOperatorSplitExplicit(numL(k), B);
      Add(coef(k), B, A);
    }

  mat_lu.Factorize(A);
}

// application du preconditionneur Y = T(L)^{-1} X
void ModeEsNonLinear_DD_DirichletEigenProblem::ApplyPreconditioning(const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
{
  // cout << "Pas encore bien écrit : ApplyPreconditioning" << endl;
  // abort();
  //cout << "appel apply" << Norm2(X) << endl;
  mat_lu.Solve(Y, X);
}


void ModeEsNonLinear_DD_DirichletEigenProblem
::ComputeOperatorSplitExplicit(int num, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  if (num == 0)
    {
      A.Clear();
      A.Reallocate(nb_rows, nb_rows);

      int inc = nb_dof_E3 + nb_dof_H3;
      for (int ref = 0 ; ref < nb_couches ; ref++)
        if (index_used(ref))
          {
            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);

            for (int iDir = 0; iDir < nb_dof_E3; iDir++)
              {
                int i = DofKeptDir(iDir); // i = ligne globale
                int nM = vec_Mh(ref).GetRowSize(i);

                Vector<int> col; Vector<Complex_wp> val;
                int num = 0;
                col.Reallocate(nM); val.Reallocate(nM);
                // parties avec Uh
                for (int j = 0; j < nM; j++)
                  {
                    int j_Glob = vec_Mh(ref).Index(i, j);
                    int jDir = IndexNonDirichlet(j_Glob);
                    if (jDir >= 0)
                      {
                        int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                        col(num) = inc + j_loc_E3;
                        // partie - eps_sigma Mh V2t
                        val(num) = eps_sigma * vec_Mh(ref).Value(i, j);
                        num++;
                      }
                  }
                A.AddInteractionRow(iDir, num, col, val);
              }

            for (int i = 0; i < nb_dof_H3; i++)
              {
                int nM = vec_Mh(ref).GetRowSize(i);

                Vector<int> col; Vector<Complex_wp> val;
                int num = 0;
                col.Reallocate(nM); val.Reallocate(nM);
                // partie avec Kh
                for (int j = 0; j < nM; j++)
                  {
                    int j_Glob = vec_Mh(ref).Index(i, j);
                    int j_loc = IndexDof(ref)(j_Glob);
                    col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
                    // partie eps_sigma Kh V2
                    val(num) = mu_tilde * vec_Mh(ref).Value(i, j);
                    num++;
                  }
                A.AddInteractionRow(nb_dof_E3 + i, num, col, val);
              }

            // Variables auxiliaires
            for (int i = 0 ; i < NumGlobE3(ref).GetM() ; i++)
              {
                int i_Glob = NumGlobE3(ref)(i);
                int iDir = IndexNonDirichlet(i_Glob);
                Vector<int> col(2); Vector<Complex_wp> val(2);
                col(0) = inc + i;
                val(0) = 1.0;
                col(1) = iDir;
                val(1) = -1.0;
                A.AddInteractionRow(inc + i, 2, col, val);
              }
            for (int i = 0 ; i < NumGlob(ref).GetM() ; i++)
              {
                int i_Glob = NumGlob(ref)(i);
                Vector<int> col(2); Vector<Complex_wp> val(2);
                col(0) = inc + nb_dof_coucheE3(ref) + i;
                val(0) = 1.0;
                col(1) = nb_dof_E3 + i_Glob;
                val(1) = -1.0;
                A.AddInteractionRow(inc + nb_dof_coucheE3(ref) + i, 2, col, val);
              }

            inc += nb_dof_coucheE3(ref) + nb_dof_couche(ref);
          }
    }
  else if (num%4 == 1)
    {
      int ref = this->split_reference(num/4);
      Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
      Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
      Complex_wp k2 = -eps_sigma*mu_tilde;
      Complex_wp tau = var.torsion;

      int inc = nb_dof_E3 + nb_dof_H3;
      for (int i = 0 ; i < ref ; i++)
        {
          if (index_used(i))
            inc += nb_dof_coucheE3(i) + nb_dof_couche(i);
        }

      A.Clear();
      A.Reallocate(nb_rows, nb_rows);

      // on traite les lignes avec V2
      for (int iDir = 0; iDir < nb_dof_E3; iDir++)
        {
          int i = DofKeptDir(iDir); // i = ligne globale
          int nMr2 = vec_Mr2h(ref).GetRowSize(i);
          int nRj = vec_Rjh(ref).GetRowSize(i);
          int nRi = vec_Rih(ref).GetRowSize(i);
          int nK = vec_Kh(ref).GetRowSize(i);
          int nD = vec_Dh(ref).GetRowSize(i);

          Vector<int> col; Vector<Complex_wp> val;
          int num = 0;
          col.Reallocate(nMr2); val.Reallocate(nMr2);
          // parties avec Uh
          for (int j = 0; j < nMr2; j++)
            {
              int j_Glob = vec_Mr2h(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie - eps_sigma Mh V2t
                  val(num) = - eps_sigma * k2 * tau * tau * vec_Mr2h(ref).Value(i, j);
                  // val(num) = - eps_sigma * tau * tau * vec_Mr2h(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nRj); val.Reallocate(nRj);
          // partie avec Kh
          for (int j = 0; j < nRj; j++)
            {
              int j_Glob = vec_Rjh(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie eps_sigma Kh V2
              val(num) = - k2 * tau * vec_Rjh(ref).Value(i, j);
              // val(num) = - tau * vec_Rjh(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nRi); val.Reallocate(nRi);
          // partie avec Kh
          for (int j = 0; j < nRi; j++)
            {
              int j_Glob = vec_Rih(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie eps_sigma Kh V2
              val(num) = - k2 * tau * vec_Rih(ref).Value(i, j);
              // val(num) = - tau * vec_Rih(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nK); val.Reallocate(nK);
          // partie avec Dh
          for (int j = 0; j < nK; j++)
            {
              int j_Glob = vec_Kh(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie - eps_sigma Mh V2t
                  val(num) = eps_sigma * vec_Kh(ref).Value(i, j);
                  // val(num) = -1.0/mu_tilde * vec_Kh(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nD); val.Reallocate(nD);
          // partie avec Dh
          for (int j = 0; j < nD; j++)
            {
              int j_Glob = vec_Dh(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma tau^2 Dh V2
                  val(num) = - eps_sigma * tau * tau * vec_Dh(ref).Value(i, j);
                  // val(num) = 1.0 / mu_tilde * tau * tau * vec_Dh(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(iDir, num, col, val);
        }

      // on traite les lignes avec W2
      for (int i = 0; i < nb_dof_H3; i++)
        {
          int nMr2 = vec_Mr2h(ref).GetRowSize(i);
          int nRj = vec_Rjh(ref).GetRowSize(i);
          int nRi = vec_Rih(ref).GetRowSize(i);
          int nK = vec_Kh(ref).GetRowSize(i);
          int nD = vec_Dh(ref).GetRowSize(i);

          Vector<int> col; Vector<Complex_wp> val;

          int num = 0;
          col.Reallocate(nMr2); val.Reallocate(nMr2);
          // parties avec Uh
          for (int j = 0; j < nMr2; j++)
            {
              int j_Glob = vec_Mr2h(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie - mu_tilde Mh W2t
              val(num) = - mu_tilde * k2 * tau * tau * vec_Mr2h(ref).Value(i, j);
              // val(num) = - mu_tilde * tau * tau * vec_Mr2h(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nRj); val.Reallocate(nRj);
          // partie avec Dh
          for (int j = 0; j < nRj; j++)
            {
              int j_Glob = vec_Rjh(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma tau^2 Dh V2
                  val(num) = k2 * tau * vec_Rjh(ref).Value(i, j);
                  // val(num) = tau * vec_Rjh(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nRi); val.Reallocate(nRi);
          // partie avec Dh
          for (int j = 0; j < nRi; j++)
            {
              int j_Glob = vec_Rih(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma tau^2 Dh V2
                  val(num) = k2 * tau * vec_Rih(ref).Value(i, j);
                  // val(num) = tau * vec_Rih(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nK); val.Reallocate(nK);
          // partie avec Kh
          for (int j = 0; j < nK; j++)
            {
              int j_Glob = vec_Kh(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              // int jDir = IndexNonDirichlet(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie mu_tilde Kh W2
              val(num) = mu_tilde * vec_Kh(ref).Value(i, j);
              // val(num) = -1.0/eps_sigma * vec_Kh(ref).Value(i, j);
              num++;
            }

          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nD); val.Reallocate(nD);
          // partie avec Dh
          for (int j = 0; j < nD; j++)
            {
              int j_Glob = vec_Dh(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie mu_tilde tau^2 Dh W2
              val(num) = - mu_tilde * tau * tau * vec_Dh(ref).Value(i, j);
              // val(num) = 1.0/eps_sigma * tau * tau * vec_Dh(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);
        }

      // Variables auxiliaires
      for (int i = 0 ; i < NumGlobE3(ref).GetM() ; i++)
        {
          int i_Glob = NumGlobE3(ref)(i);
          Vector<int> col(1); Vector<Complex_wp> val(1);
          col(0) = inc + i;
          val(0) = - 2.0*k2*tau*tau * (vec_Pos(i_Glob)(0)*vec_Pos(i_Glob)(0) + vec_Pos(i_Glob)(1)*vec_Pos(i_Glob)(1));
          A.AddInteractionRow(inc + i, 1, col, val);
        }
      for (int i = 0 ; i < NumGlob(ref).GetM() ; i++)
        {
          int i_Glob = NumGlob(ref)(i);
          Vector<int> col(1); Vector<Complex_wp> val(1);
          col(0) = inc + nb_dof_coucheE3(ref) + i;
          val(0) = - 2.0*k2*tau*tau * (vec_Pos(i_Glob)(0)*vec_Pos(i_Glob)(0) + vec_Pos(i_Glob)(1)*vec_Pos(i_Glob)(1));
          A.AddInteractionRow(inc + nb_dof_coucheE3(ref) + i, 1, col, val);
        }
    }
  else if (num%4 == 2) // terme après f_3 = i beta /D_L(beta)
    {
      int ref = this->split_reference(num/4);
      Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
      Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
      Complex_wp k2 = -eps_sigma*mu_tilde;
      Complex_wp tau = var.torsion;

      int inc = nb_dof_E3 + nb_dof_H3;
      for (int i = 0 ; i < ref ; i++)
        {
          if (index_used(i))
            inc += nb_dof_coucheE3(i) + nb_dof_couche(i);
        }

      A.Clear();
      A.Reallocate(nb_rows, nb_rows);

      // on traite les lignes avec V2
      for (int iDir = 0; iDir < nb_dof_E3; iDir++)
        {
          int i = DofKeptDir(iDir); // i = ligne globale
          int nSj = vec_Sjh(ref).GetRowSize(i);
          int nSi = vec_Sih(ref).GetRowSize(i);
          int nC = vec_Ch(ref).GetRowSize(i);

          Vector<int> col; Vector<Complex_wp> val;
          int num = 0;
          col.Reallocate(nSj); val.Reallocate(nSj);
          // partie avec Kh
          for (int j = 0; j < nSj; j++)
            {
              int j_Glob = vec_Sjh(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma Kh V2
                  val(num) = eps_sigma * tau * vec_Sjh(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nSi); val.Reallocate(nSi);
          // partie avec Kh
          for (int j = 0; j < nSi; j++)
            {
              int j_Glob = vec_Sih(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma Kh V2
                  val(num) = - eps_sigma * tau * vec_Sih(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nC); val.Reallocate(nC);
          // partie avec Dh
          for (int j = 0; j < nC; j++)
            {
              int j_Glob = vec_Ch(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie - eps_sigma Mh V2t
              val(num) = vec_Ch(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(iDir, num, col, val);
        }

      // on traite les lignes avec W2
      for (int i = 0; i < nb_dof_H3; i++)
        {
          int nSj = vec_Sjh(ref).GetRowSize(i);
          int nSi = vec_Sih(ref).GetRowSize(i);
          int nC = vec_Ch(ref).GetRowSize(i);

          Vector<int> col; Vector<Complex_wp> val;

          int num = 0;

          num = 0;
          col.Reallocate(nSj); val.Reallocate(nSj);
          // partie avec Dh
          for (int j = 0; j < nSj; j++)
            {
              int j_Glob = vec_Sjh(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie eps_sigma tau^2 Dh V2
              val(num) = mu_tilde * tau * vec_Sjh(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nSi); val.Reallocate(nSi);
          // partie avec Dh
          for (int j = 0; j < nSi; j++)
            {
              int j_Glob = vec_Sih(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie eps_sigma tau^2 Dh V2
              val(num) = - mu_tilde * tau * vec_Sih(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nC); val.Reallocate(nC);
          // partie avec Kh
          for (int j = 0; j < nC; j++)
            {
              int j_Glob = vec_Ch(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma Kh V2
                  val(num) = - vec_Ch(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);
        }
    }
  else if (num%4 == 3) // terme après f_3 = i beta /D_L(beta)
    {
      int ref = this->split_reference(num/4);
      Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
      Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
      Complex_wp k2 = -eps_sigma*mu_tilde;
      Complex_wp tau = var.torsion;

      int inc = nb_dof_E3 + nb_dof_H3;
      for (int i = 0 ; i < ref ; i++)
        {
          if (index_used(i))
            inc += nb_dof_coucheE3(i) + nb_dof_couche(i);
        }

      A.Clear();
      A.Reallocate(nb_rows, nb_rows);

      // on traite les lignes avec V2
      for (int iDir = 0; iDir < nb_dof_E3; iDir++)
        {
          int i = DofKeptDir(iDir); // i = ligne globale
          int nSjr2 = vec_Sjr2h(ref).GetRowSize(i);
          int nSir2 = vec_Sir2h(ref).GetRowSize(i);
          int nCr2 = vec_Cr2h(ref).GetRowSize(i);
          int nSi = vec_Sih(ref).GetRowSize(i);

          Vector<int> col; Vector<Complex_wp> val;
          int num = 0;
          col.Reallocate(nSjr2); val.Reallocate(nSjr2);
          // partie avec Kh
          for (int j = 0; j < nSjr2; j++)
            {
              int j_Glob = vec_Sjr2h(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma Kh V2
                  val(num) = - eps_sigma * tau * tau * tau * vec_Sjr2h(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nSir2); val.Reallocate(nSir2);
          // partie avec Kh
          for (int j = 0; j < nSir2; j++)
            {
              int j_Glob = vec_Sir2h(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma Kh V2
                  val(num) = eps_sigma * tau * tau * tau * vec_Sir2h(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nCr2); val.Reallocate(nCr2);
          // partie avec Dh
          for (int j = 0; j < nCr2; j++)
            {
              int j_Glob = vec_Cr2h(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie - eps_sigma Mh V2t
              val(num) = - tau * tau * vec_Cr2h(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nSi); val.Reallocate(nSi);
          // partie avec Dh
          for (int j = 0; j < nSi; j++)
            {
              int j_Glob = vec_Sih(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie - eps_sigma Mh V2t
              val(num) = - 4.0 * tau * tau * vec_Sih(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(iDir, num, col, val);
        }

      // on traite les lignes avec W2
      for (int i = 0; i < nb_dof_H3; i++)
        {
          int nSjr2 = vec_Sjr2h(ref).GetRowSize(i);
          int nSir2 = vec_Sir2h(ref).GetRowSize(i);
          int nCr2 = vec_Cr2h(ref).GetRowSize(i);
          int nSi = vec_Sih(ref).GetRowSize(i);

          Vector<int> col; Vector<Complex_wp> val;

          int num = 0;

          num = 0;
          col.Reallocate(nSjr2); val.Reallocate(nSjr2);
          // partie avec Dh
          for (int j = 0; j < nSjr2; j++)
            {
              int j_Glob = vec_Sjr2h(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie eps_sigma tau^2 Dh V2
              val(num) = - mu_tilde * tau * tau * tau * vec_Sjr2h(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nSir2); val.Reallocate(nSir2);
          // partie avec Dh
          for (int j = 0; j < nSir2; j++)
            {
              int j_Glob = vec_Sir2h(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie eps_sigma tau^2 Dh V2
              val(num) = mu_tilde * tau * tau * tau * vec_Sir2h(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nCr2); val.Reallocate(nCr2);
          // partie avec Kh
          for (int j = 0; j < nCr2; j++)
            {
              int j_Glob = vec_Cr2h(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma Kh V2
                  val(num) = tau * tau * vec_Cr2h(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nSi); val.Reallocate(nSi);
          // partie avec Kh
          for (int j = 0; j < nSi; j++)
            {
              int j_Glob = vec_Sih(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma Kh V2
                  val(num) = 4.0 * tau * tau * vec_Sih(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);
        }
    }
  else
    {
      int ref = this->split_reference(num/4-1);
      Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
      Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
      Complex_wp k2 = -eps_sigma*mu_tilde;
      Complex_wp tau = var.torsion;

      int inc = nb_dof_E3 + nb_dof_H3;
      for (int i = 0 ; i < ref ; i++)
        {
          if (index_used(i))
            inc += nb_dof_coucheE3(i) + nb_dof_couche(i);
        }

      A.Clear();
      A.Reallocate(nb_rows, nb_rows);

      // on traite les lignes avec V2
      for (int iDir = 0; iDir < nb_dof_E3; iDir++)
        {
          int i = DofKeptDir(iDir); // i = ligne globale
          int nRjr2 = vec_Rjr2h(ref).GetRowSize(i);
          int nMr2 = vec_Mr2h(ref).GetRowSize(i);
          int nRir2 = vec_Rir2h(ref).GetRowSize(i);
          int nKr2 = vec_Kr2h(ref).GetRowSize(i);
          int nRi = vec_Rih(ref).GetRowSize(i);
          int nDr2 = vec_Dr2h(ref).GetRowSize(i);

          Vector<int> col; Vector<Complex_wp> val;
          int num = 0;
          col.Reallocate(nMr2); val.Reallocate(nMr2);


          num = 0;
          col.Reallocate(nRjr2); val.Reallocate(nRjr2);
          // partie avec Kh
          for (int j = 0; j < nRjr2; j++)
            {
              int j_Glob = vec_Rjr2h(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie eps_sigma Kh V2
              val(num) = k2 * tau * tau * tau * vec_Rjr2h(ref).Value(i, j);
              // val(num) = tau * tau * tau * vec_Rjr2h(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nMr2); val.Reallocate(nMr2);
          // partie avec Kh
          for (int j = 0; j < nMr2; j++)
            {
              int j_Glob = vec_Mr2h(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie eps_sigma Kh V2
              val(num) = 4.0 * k2 * tau * tau * tau * vec_Mr2h(ref).Value(i, j);
              // val(num) = 4.0 * tau * tau * tau * vec_Mr2h(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nRir2); val.Reallocate(nRir2);
          // partie avec Kh
          for (int j = 0; j < nRir2; j++)
            {
              int j_Glob = vec_Rir2h(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie eps_sigma Kh V2
              val(num) = k2 * tau * tau * tau * vec_Rir2h(ref).Value(i, j);
              // val(num) = tau * tau * tau * vec_Rir2h(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nKr2); val.Reallocate(nKr2);
          // parties avec Uh
          for (int j = 0; j < nKr2; j++)
            {
              int j_Glob = vec_Kr2h(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie - eps_sigma Mh V2t
                  val(num) = - eps_sigma * tau * tau * vec_Kr2h(ref).Value(i, j);
                  // val(num) = tau * tau / mu_tilde * vec_Kr2h(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nRi); val.Reallocate(nRi);
          // parties avec Uh
          for (int j = 0; j < nRi; j++)
            {
              int j_Glob = vec_Rih(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie - eps_sigma Mh V2t
                  val(num) = - 4.0 * eps_sigma * tau * tau * vec_Rih(ref).Value(i, j);
                  // val(num) = 4.0 * tau * tau / mu_tilde * vec_Rih(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(iDir, num, col, val);

          num = 0;
          col.Reallocate(nDr2); val.Reallocate(nDr2);
          // parties avec Uh
          for (int j = 0; j < nDr2; j++)
            {
              int j_Glob = vec_Dr2h(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie - eps_sigma Mh V2t
                  val(num) = eps_sigma * tau * tau * tau * tau * vec_Dr2h(ref).Value(i, j);
                  // val(num) = - tau * tau * tau * tau / mu_tilde * vec_Dr2h(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(iDir, num, col, val);
        }

      // on traite les lignes avec W2
      for (int i = 0; i < nb_dof_H3; i++)
        {
          int nRjr2 = vec_Rjr2h(ref).GetRowSize(i);
          int nMr2 = vec_Mr2h(ref).GetRowSize(i);
          int nRir2 = vec_Rir2h(ref).GetRowSize(i);
          int nKr2 = vec_Kr2h(ref).GetRowSize(i);
          int nRi = vec_Rih(ref).GetRowSize(i);
          int nDr2 = vec_Dr2h(ref).GetRowSize(i);

          Vector<int> col; Vector<Complex_wp> val;

          int num = 0;

          num = 0;
          col.Reallocate(nRjr2); val.Reallocate(nRjr2);
          // partie avec Dh
          for (int j = 0; j < nRjr2; j++)
            {
              int j_Glob = vec_Rjr2h(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma tau^2 Dh V2
                  val(num) = - k2 * tau * tau * tau * vec_Rjr2h(ref).Value(i, j);
                  // val(num) = - tau * tau * tau * vec_Rjr2h(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nMr2); val.Reallocate(nMr2);
          // partie avec Dh
          for (int j = 0; j < nMr2; j++)
            {
              int j_Glob = vec_Mr2h(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma tau^2 Dh V2
                  val(num) = - 4.0 * k2 * tau * tau * tau * vec_Mr2h(ref).Value(i, j);
                  // val(num) = - 4.0 * tau * tau * tau * vec_Mr2h(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nRir2); val.Reallocate(nRir2);
          // partie avec Dh
          for (int j = 0; j < nRir2; j++)
            {
              int j_Glob = vec_Rir2h(ref).Index(i, j);
              int jDir = IndexNonDirichlet(j_Glob);
              if (jDir >= 0)
                {
                  int j_loc_E3 = IndexDofE3(ref)(j_Glob);
                  col(num) = inc + j_loc_E3;
                  // partie eps_sigma tau^2 Dh V2
                  val(num) = - k2 * tau * tau * tau * vec_Rir2h(ref).Value(i, j);
                  // val(num) = - tau * tau * tau * vec_Rir2h(ref).Value(i, j);
                  num++;
                }
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nKr2); val.Reallocate(nKr2);
          // parties avec Uh
          for (int j = 0; j < nKr2; j++)
            {
              int j_Glob = vec_Kr2h(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie - mu_tilde Mh W2t
              val(num) = - mu_tilde * tau * tau * vec_Kr2h(ref).Value(i, j);
              // val(num) = tau * tau / eps_sigma * vec_Kr2h(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nRi); val.Reallocate(nRi);
          // parties avec Uh
          for (int j = 0; j < nRi; j++)
            {
              int j_Glob = vec_Rih(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie - mu_tilde Mh W2t
              val(num) = - 4.0 * mu_tilde * tau * tau * vec_Rih(ref).Value(i, j);
              // val(num) = 4.0 * tau * tau / eps_sigma * vec_Rih(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);

          num = 0;
          col.Reallocate(nDr2); val.Reallocate(nDr2);
          // parties avec Uh
          for (int j = 0; j < nDr2; j++)
            {
              int j_Glob = vec_Dr2h(ref).Index(i, j);
              int j_loc = IndexDof(ref)(j_Glob);
              col(num) = inc + nb_dof_coucheE3(ref) + j_loc;
              // partie - mu_tilde Mh W2t
              val(num) = mu_tilde * tau * tau * tau * tau * vec_Dr2h(ref).Value(i, j);
              // val(num) = - tau * tau * tau * tau / eps_sigma * vec_Dr2h(ref).Value(i, j);
              num++;
            }
          A.AddInteractionRow(nb_dof_E3 + i, num, col, val);
        }

      // Variables auxiliaires
      for (int i = 0 ; i < NumGlobE3(ref).GetM() ; i++)
        {
          int i_Glob = NumGlobE3(ref)(i);
          Vector<int> col(1); Vector<Complex_wp> val(1);
          Complex_wp r2 = vec_Pos(i_Glob)(0)*vec_Pos(i_Glob)(0) + vec_Pos(i_Glob)(1)*vec_Pos(i_Glob)(1);
          col(0) = inc + i;
          val(0) = k2*tau*tau*tau*tau * r2 * r2;
          // val(0) = tau*tau*tau*tau * r2 * r2;
          A.AddInteractionRow(inc + i, 1, col, val);
        }
      for (int i = 0 ; i < NumGlob(ref).GetM() ; i++)
        {
          int i_Glob = NumGlob(ref)(i);
          Vector<int> col(1); Vector<Complex_wp> val(1);
          Complex_wp r2 = vec_Pos(i_Glob)(0)*vec_Pos(i_Glob)(0) + vec_Pos(i_Glob)(1)*vec_Pos(i_Glob)(1);
          col(0) = inc + nb_dof_coucheE3(ref) + i;
          val(0) = k2*tau*tau*tau*tau * r2 * r2;
          // val(0) = tau*tau*tau*tau * r2 * r2;
          A.AddInteractionRow(inc + nb_dof_coucheE3(ref) + i, 1, col, val);
        }
    }
}

void ModeEsNonLinear_DD_DirichletEigenProblem
::MltOperatorSplit(int num, const SeldonTranspose& trans, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{

  Vector<Complex_wp> EH(nb_dof_all), FG(nb_dof_all), Y(nb_dof_all);
  // DISP(nb_dof_all);
  Complex_wp one(1, 0);

  // Y.Zero();
  //
  int N = this->nb_dof_H3;
  Vector<Complex_wp> E3(N), H3(N), F(N), G(N), Y1(N), Y2(N), Var1(N), Var2(N);
  // ExtractDirichletAuxVariable(Xdir, 2, E3, H3, F, G);
  Y1.Zero(); Y2.Zero(); Var1.Zero(), Var2.Zero();

  if (num == 0)
    {
      Ydir.Reallocate(nb_rows); Ydir.Zero();

      int inc = nb_dof_E3 + nb_dof_H3;
      for (int ref = 0 ; ref < nb_couches ; ref++)
        if (index_used(ref))
          {
            ExtractDirichletAuxVariable(Xdir, ref, E3, H3, F, G);

            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);

            // Équation principale
            MltAdd(eps_sigma, vec_Mh(ref), F, one, Y1);
            MltAdd(mu_tilde, vec_Mh(ref), G, one, Y2);

            // Variables auxiliaires
            for (int i = 0; i < NumGlobE3(ref).GetM(); i++)
              {
                Ydir(inc + i) = Xdir(inc + i) - Xdir(NumGlobE3(ref)(i)); // <- F - E3
              }
            for (int i = 0; i < NumGlob(ref).GetM(); i++)
              {
                Ydir(inc + nb_dof_coucheE3(ref) + i) = Xdir(inc + nb_dof_coucheE3(ref) + i) - Xdir(nb_dof_E3 + NumGlob(ref)(i)); // <- G - H3
              }
            inc += nb_dof_coucheE3(ref) + nb_dof_couche(ref);
          }

      Y.Zero();
      for (int i = 0 ; i < N ; i++)
        {
          Y(i) = Y1(i);
          Y(N + i) = Y2(i);
        }
      ExpandDirichlet(Y, Ydir); // cette méthode ne réinitialise pas le vecteur, elle n'écrase donc pas le variables auxiliaires
    }
  else if (num%4 == 1) // terme après f_2 = 1/D_L(beta)
    {
      int ref = this->split_reference(num/4);
      ExtractDirichletAuxVariable(Xdir, ref, E3, H3, F, G);
      Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
      Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
      Complex_wp k2 = -eps_sigma*mu_tilde;
      Complex_wp tau = var.torsion;

      // split : 1/D(Beta)

      // Équations principales
      MltAdd(- eps_sigma * k2 * tau * tau, vec_Mr2h(ref), F, one, Y1);
      MltAdd(- k2 * tau, vec_Rjh(ref), G, one, Y1);
      MltAdd(- k2 * tau, vec_Rih(ref), G, one, Y1);
      MltAdd(eps_sigma, vec_Kh(ref), F, one, Y1);
      MltAdd(-eps_sigma*tau*tau, vec_Dh(ref), F, one, Y1);

      MltAdd(- mu_tilde * k2 * tau * tau, vec_Mr2h(ref), G, one, Y2);
      MltAdd(k2 * tau, vec_Rjh(ref), F, one, Y2);
      MltAdd(k2 * tau, vec_Rih(ref), F, one, Y2);
      MltAdd(mu_tilde, vec_Kh(ref), G, one, Y2);
      MltAdd(-mu_tilde*tau*tau, vec_Dh(ref), G, one, Y2);

      // Variables auxiliaires
      Mlt(R2, F, Var1); Mlt(- 2.0*k2*tau*tau, Var1);
      Mlt(R2, G, Var2); Mlt(- 2.0*k2*tau*tau, Var2);

      // split : k2/D(Beta)

      // Équations principales
      // MltAdd(- eps_sigma * tau * tau, vec_Mr2h(ref), F, one, Y1);
      // MltAdd(- tau, vec_Rjh(ref), G, one, Y1);
      // MltAdd(- tau, vec_Rih(ref), G, one, Y1);
      // MltAdd(-1.0/mu_tilde, vec_Kh(ref), F, one, Y1);
      // MltAdd(tau*tau / mu_tilde, vec_Dh(ref), F, one, Y1);
      //
      // MltAdd(- mu_tilde * tau * tau, vec_Mr2h(ref), G, one, Y2);
      // MltAdd(tau, vec_Rjh(ref), F, one, Y2);
      // MltAdd(tau, vec_Rih(ref), F, one, Y2);
      // MltAdd(- 1.0/eps_sigma, vec_Kh(ref), G, one, Y2);
      // MltAdd(tau*tau/eps_sigma, vec_Dh(ref), G, one, Y2);
      //
      // // Variables auxiliaires
      // Mlt(R2, F, Var1); Mlt(- 2.0*tau*tau, Var1);
      // Mlt(R2, G, Var2); Mlt(- 2.0*tau*tau, Var2);

      ExpandDirichletAuxVariable(ref, Y1, Y2, Var1, Var2, Ydir);
    }
  else if (num%4 == 2) // terme après f_3 = i beta /D_L(beta)
    {
      int ref = this->split_reference(num/4);
      ExtractDirichletAuxVariable(Xdir, ref, E3, H3, F, G);
      Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
      Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
      Complex_wp k2 = -eps_sigma*mu_tilde;
      Complex_wp tau = var.torsion;

      MltAdd(eps_sigma * tau, vec_Sjh(ref), F, one, Y1);
      MltAdd(-eps_sigma * tau, vec_Sih(ref), F, one, Y1);
      MltAdd(one, vec_Ch(ref), G, one, Y1);

      MltAdd(mu_tilde * tau, vec_Sjh(ref), G, one, Y2);
      MltAdd(-mu_tilde * tau, vec_Sih(ref), G, one, Y2);
      MltAdd(-one, vec_Ch(ref), F, one, Y2);

      Vector<Complex_wp> ver0(N); ver0.Zero();
      ExpandDirichletAuxVariable(ref, Y1, Y2, ver0, ver0, Ydir);
    }
  else if (num%4 == 3) // terme après f_4 = i beta k2L /D_L(beta^2
    {
      int ref = this->split_reference(num/4);
      ExtractDirichletAuxVariable(Xdir, ref, E3, H3, F, G);
      Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
      Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
      Complex_wp k2 = -eps_sigma*mu_tilde;
      Complex_wp tau = var.torsion;

      MltAdd(- eps_sigma * tau * tau * tau, vec_Sjr2h(ref), F, one, Y1);
      MltAdd(eps_sigma * tau * tau * tau, vec_Sir2h(ref), F, one, Y1);
      MltAdd(- tau * tau, vec_Cr2h(ref), G, one, Y1);
      MltAdd(- 4.0 * tau * tau, vec_Sih(ref), G, one, Y1);

      MltAdd(- mu_tilde * tau * tau * tau, vec_Sjr2h(ref), G, one, Y2);
      MltAdd(mu_tilde * tau * tau * tau, vec_Sir2h(ref), G, one, Y2);
      MltAdd(tau * tau, vec_Cr2h(ref), F, one, Y2);
      MltAdd(4.0 * tau * tau, vec_Sih(ref), F, one, Y2);

      Vector<Complex_wp> ver0(N); ver0.Zero();
      ExpandDirichletAuxVariable(ref, Y1, Y2, ver0, ver0, Ydir);
    }
  else // terme après f_5 = k2L /D_L(beta)^2
    {
      int ref = this->split_reference(num/4-1);
      ExtractDirichletAuxVariable(Xdir, ref, E3, H3, F, G);
      Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
      Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
      Complex_wp k2 = -eps_sigma*mu_tilde;
      Complex_wp tau = var.torsion;

      // split : k2/D(Beta)^2

      // Équation principale
      MltAdd(k2 * tau * tau * tau, vec_Rjr2h(ref), G, one, Y1);
      MltAdd(4.0 * k2 * tau * tau * tau, vec_Mr2h(ref), G, one, Y1);
      MltAdd(k2 * tau * tau * tau, vec_Rir2h(ref), G, one, Y1);
      MltAdd(- eps_sigma * tau * tau, vec_Kr2h(ref), F, one, Y1);
      MltAdd(- 4.0 * eps_sigma * tau * tau, vec_Rih(ref), F, one, Y1);
      MltAdd(eps_sigma * tau * tau * tau * tau, vec_Dr2h(ref), F, one, Y1);

      MltAdd(- k2 * tau * tau * tau, vec_Rjr2h(ref), F, one, Y2);
      MltAdd(- 4.0 * k2 * tau * tau * tau, vec_Mr2h(ref), F, one, Y2);
      MltAdd(- k2 * tau * tau * tau, vec_Rir2h(ref), F, one, Y2);
      MltAdd(- mu_tilde * tau * tau, vec_Kr2h(ref), G, one, Y2);
      MltAdd(- 4.0 * mu_tilde * tau * tau, vec_Rih(ref), G, one, Y2);
      MltAdd(mu_tilde * tau * tau * tau * tau, vec_Dr2h(ref), G, one, Y2);

      // Variables auxiliaires
      Mlt(R2, F, Var1); Mlt(R2, Var1, Var1); Mlt(k2*tau*tau*tau*tau, Var1);
      Mlt(R2, G, Var2); Mlt(R2, Var2, Var2); Mlt(k2*tau*tau*tau*tau, Var2);

      // split : k4/D(Beta)^2

      // // Équation principale
      // MltAdd(tau * tau * tau, vec_Rjr2h(ref), G, one, Y1);
      // MltAdd(4.0 * tau * tau * tau, vec_Mr2h(ref), G, one, Y1);
      // MltAdd(tau * tau * tau, vec_Rir2h(ref), G, one, Y1);
      // MltAdd(tau * tau/mu_tilde, vec_Kr2h(ref), F, one, Y1);
      // MltAdd(4.0 * tau * tau/mu_tilde, vec_Rih(ref), F, one, Y1);
      // MltAdd(-tau * tau * tau * tau/mu_tilde, vec_Dr2h(ref), F, one, Y1);
      //
      // MltAdd(- tau * tau * tau, vec_Rjr2h(ref), F, one, Y2);
      // MltAdd(- 4.0 * tau * tau * tau, vec_Mr2h(ref), F, one, Y2);
      // MltAdd(- tau * tau * tau, vec_Rir2h(ref), F, one, Y2);
      // MltAdd(tau * tau / eps_sigma, vec_Kr2h(ref), G, one, Y2);
      // MltAdd(4.0 * tau * tau / eps_sigma, vec_Rih(ref), G, one, Y2);
      // MltAdd(- tau * tau * tau * tau/eps_sigma, vec_Dr2h(ref), G, one, Y2);
      //
      // // Variables auxiliaires
      // Mlt(R2, F, Var1); Mlt(R2, Var1, Var1); Mlt(tau*tau*tau*tau, Var1);
      // Mlt(R2, G, Var2); Mlt(R2, Var2, Var2); Mlt(tau*tau*tau*tau, Var2);

      ExpandDirichletAuxVariable(ref, Y1, Y2, Var1, Var2, Ydir);
    }
}
