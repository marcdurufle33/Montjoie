// initialisation du probleme
ModeEsNonLinear2NDirichletEigenProblem
::ModeEsNonLinear2NDirichletEigenProblem(int N, ModeEs_Solver& var0,
                            Vector<int>& DofDir, Vector<int>& IndexDir,
                            bool use_split)
  :
  var(var0), DofKeptDir(DofDir), IndexNonDirichlet(IndexDir)
{
  if (use_split)
    {
      cout << "Formulation split inconnue pour le modèle 2N." << endl;
      abort();
    }
  // pour Es, condition de Dirichlet
  nb_dof_E3 = DofKeptDir.GetM();

  // pour Hs, condition de Neumann
  nb_dof_H3 = N;

  countIteration = 0;

  cout << "Initiation des matrices" << endl;
  var.ComputeFemMatricesNonLin2N(var.shift, vec_MhB, vec_KhB, vec_ChB, vec_RjhB, vec_RihB, vec_SjhB, vec_SihB, vec_DhB);

  // taille du pb aux valeurs propres
  nb_rows = nb_dof_E3;
  nb_dof_all = N;
  nb_rows += nb_dof_H3;
  nb_dof_all += N;

  // construction des numerotations
  cout << "N = " << N << endl;

  index_used.Reallocate(var.ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < var.mesh.GetNbElt(); i++)
    index_used(var.mesh.Element(i).GetReference()) = true;

  nb_couches = index_used.GetM();
  nb_dof_couche.Reallocate(nb_couches);
  nb_dof_coucheE3.Reallocate(nb_couches);
  nb_dof_couche.Zero(); nb_dof_coucheE3.Zero();
  IndexDofE3.Reallocate(nb_couches);
  NumGlobE3.Reallocate(nb_couches);
  IndexDof.Reallocate(nb_couches);
  NumGlob.Reallocate(nb_couches);

  for (int ref = 0; ref < nb_couches; ref++)
    if (index_used(ref))
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

  this->Init(nb_rows);
  DISP(nb_rows);

  // on calcule les poles du pb aux valeurs propres
  VectComplex_wp s;

  DistributedMatrix<Complexe, Prop, Storage> A;
  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
  {
    if (index_used(ref))
      {
        Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
        Complex_wp mu_tilde = - Iwp * var.omega * var.ref_mu(ref);
        Complex_wp k2 = - eps_sigma*mu_tilde;
        // Complex_wp k2_c = (var.ref_epsilon(ref) + Iwp * var.ref_sigma(ref) / var.omega) * var.ref_mu(ref);
        Complex_wp z = sqrt(-k2) / var.omega;
        // Real_wp tau = var.ref_tau(ref);
        // Matrix<Complexe, Prop, Storage> temp;
        s.PushBack(z); s.PushBack(-z);
      }
  }

  // this->SetSingularities(s); //DISP(s);
}

void ModeEsNonLinear2NDirichletEigenProblem::RemoveDirichletRowCol(const Matrix<Complexe, Prop, Storage>& A,
                                                        Matrix<Complexe, Prop, Storage>& B)
{
  cout << "RemoveDirichletRowCol" << endl;
  abort();
}

void ModeEsNonLinear2NDirichletEigenProblem::ExtractDirichlet(const VectComplex_wp& Xdir, VectComplex_wp& X)
{
  X.Zero();
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    X(DofKeptDir(i)) = Xdir(i);

  for (int i = 0; i < nb_dof_H3; i++)
    X(nb_dof_H3 + i) = Xdir(nb_dof_E3 + i);
}

void ModeEsNonLinear2NDirichletEigenProblem::ExpandDirichlet(const VectComplex_wp& X, VectComplex_wp& Xdir)
{
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    Xdir(i) = X(DofKeptDir(i));

  for (int i = 0; i < nb_dof_H3; i++)
    Xdir(nb_dof_E3 + i) = X(nb_dof_H3 + i);
}


void ModeEsNonLinear2NDirichletEigenProblem::ComputeOperator(const Complex_wp& L)
{

}

// void ModeEsNonLinear2NDirichletEigenProblem::Mult(const Matrix<Complexe, Prop, Storage>& A, const  Matrix<Complexe, Symmetric, DiagonalRow>& D, Matrix<Complexe, Prop, Storage>& M)
// {
//   int N = nb_dof_H3;
//   M.Reallocate(N,N);
//   for (int i = 0 ; i < N ; i++)
//   {
//     int n = A.GetRowSize(i);
//     Vector<int> col(n); Vector<Complex_wp> val(n);
//     int num = 0;
//     for (int j0 = 0 ; j0 < n ; j0++)
//     {
//       int j = A.Index(i, j0);
//       col(num) = j;
//       val(num) = A.Value(i,j0)*D(j,j);
//       num++;
//     }
//     M.AddInteractionRow(i, num, col, val);
//   }
// }


void ModeEsNonLinear2NDirichletEigenProblem::ComputeEsHs(const Complex_wp& L,
                                                  const Vector<Complex_wp>& eigenvector,
                                                  Vector<Complex_wp>& Es,
                                                  Vector<Complex_wp>& Hs,
                                                  Vector<Complex_wp>& E1,
                                                  Vector<Complex_wp>& E2,
                                                  Vector<Complex_wp>& H1,
                                                  Vector<Complex_wp>& H2)
{
    var.ComputeFemMatricesNonLin2N(L, vec_MhB, vec_KhB, vec_ChB, vec_RjhB, vec_RihB, vec_SjhB, vec_SihB, vec_DhB);
    Vector<Matrix<Complexe, Prop, Storage>> vec_Mh;
    var.ComputeMassMatrix(vec_Mh);
    int N = var.mesh_num.GetNbDof();
    Matrix<Complexe, General, ArrayRowSparse> MM, MMe;
    Complex_wp beta = L * var.omega;
    Complex_wp one(1, 0);

    // On récupère E3 et H3
    VectComplex_wp x(2*N), E3(N), H3(N);
    E3.Zero(); H3.Zero();
    ExtractDirichlet(eigenvector, x);
    for (int i = 0; i < nb_dof_H3; i++)
      {
        E3.Get(i) = x(i);
        H3.Get(i) = x(nb_dof_H3 + i);
      }

    // Vector<Complex_wp> grE3vXi;
    // Vector<Complex_wp> grH3sXi;

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
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = - eps_sigma*mu_tilde;
          Complex_wp DBeta = - k2 + beta * beta;
          Complex_wp tau = var.torsion;

          MltAdd(one, vec_MhB(ref), E3, one, Es);
          MltAdd(mu_tilde * tau/DBeta, vec_RjhB(ref), H3, one, Es);
          MltAdd(Iwp * beta * tau/DBeta, vec_SjhB(ref), E3, one, Es);

          MltAdd(one, vec_MhB(ref), H3, one, Hs);
          MltAdd(-eps_sigma * tau/DBeta, vec_RjhB(ref), E3, one, Hs);
          MltAdd(Iwp * beta * tau/DBeta, vec_SjhB(ref), H3, one, Hs);

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

    mat_lu.Factorize(MM);
    mat_lu.Solve(Hs);
}


// calcul de Y = T(L) X
void ModeEsNonLinear2NDirichletEigenProblem::MltOperator(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  DISP(L);
  this->currentL = L;
  countIteration++;
  if (countIteration%10 ==0)
    cout << "MltOperator " << countIteration << endl;
  var.ComputeFemMatricesNonLin2N(L, vec_MhB, vec_KhB, vec_ChB, vec_RjhB, vec_RihB, vec_SjhB, vec_SihB, vec_DhB);

  Complex_wp one(1, 0);
  Complex_wp beta = L*var.omega;

  Vector<Complex_wp> X(nb_dof_all), Y(nb_dof_all);
  ExtractDirichlet(Xdir, X);
  Y.Zero();

  // le vecteur calculé est décomposé en deux : Ez, puis Hz
  int N = nb_dof_H3;
  Vector<Complex_wp> E3(N), H3(N), Y1(N), Y2(N);
  Y1.Zero(); Y2.Zero();
  //cout << "Je passe à : MltOperator" << endl;
  for (int i = 0 ; i < N; i++)
    {
      E3(i) = X(i);
      H3(i) = X(N + i);
    }

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    if (index_used(ref))
      {
        Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
        Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
        Complex_wp k2 = -eps_sigma*mu_tilde;
        Complex_wp Delta = -k2 + beta*beta;
        Complex_wp tau = var.torsion;

        MltAdd(eps_sigma, vec_MhB(ref), E3, one, Y1);
        MltAdd(mu_tilde, vec_MhB(ref), H3, one, Y2);
        MltAdd(Iwp * beta * eps_sigma * tau / Delta, vec_SjhB(ref), E3, one, Y1);
        MltAdd(Iwp * beta * mu_tilde * tau / Delta, vec_SjhB(ref), H3, one, Y2);
        MltAdd(-k2 * tau / Delta, vec_RjhB(ref), H3, one, Y1);
        MltAdd(k2 * tau / Delta, vec_RjhB(ref), E3, one, Y2);
        MltAdd(-k2 * tau / Delta, vec_RihB(ref), H3, one, Y1);
        MltAdd(k2 * tau / Delta, vec_RihB(ref), E3, one, Y2);
        MltAdd(-Iwp * beta * eps_sigma * tau / Delta, vec_SihB(ref), E3, one, Y1);
        MltAdd(-Iwp * beta * mu_tilde * tau / Delta, vec_SihB(ref), H3, one, Y2);
        MltAdd(Iwp * beta / Delta, vec_ChB(ref), H3, one, Y1);
        MltAdd(-Iwp * beta / Delta, vec_ChB(ref), E3, one, Y2);
        MltAdd(eps_sigma / Delta, vec_KhB(ref), E3, one, Y1);
        MltAdd(mu_tilde / Delta, vec_KhB(ref), H3, one, Y2);
        MltAdd(-eps_sigma*tau*tau / Delta, vec_DhB(ref), E3, one, Y1);
        MltAdd(-mu_tilde*tau*tau / Delta, vec_DhB(ref), H3, one, Y2);
      }

  for (int i = 0 ; i < N ; i++)
    {
      Y(i) = Y1(i);
      Y(N + i) = Y2(i);
    }

  ExpandDirichlet(Y, Ydir);
}

// calcul explicite de T(L) (utilise si on calcule explicitement les matrices)
void ModeEsNonLinear2NDirichletEigenProblem
::ComputeOperatorExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "ComputeOperatorExplicit" << endl;
  abort();
}


void ModeEsNonLinear2NDirichletEigenProblem::ComputeJacobian(const Complex_wp& L)
{}

// calcul de Y = T'(L) X
void ModeEsNonLinear2NDirichletEigenProblem::MltJacobian(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  cout << "MltJacobian" << endl;
  abort();
}

// calcul explicite de T'(L)
void ModeEsNonLinear2NDirichletEigenProblem
::ComputeJacobianExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "Pas écrit : ComputeJacobianExplicit" << endl;
  abort();
}

void ModeEsNonLinear2NDirichletEigenProblem::ComputeExplicitPreconditioning(DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "ComputeExplicitPreconditioning not implemented." << endl;
  abort();
  mat_lu.Factorize(A);
}

// factorisation de T(L)
void ModeEsNonLinear2NDirichletEigenProblem::ComputePreconditioning(const Complex_wp& L)
{
  cout << "ComputePreconditioning not implemented." << endl;
  abort();
  DistributedMatrix<Complex_wp, Prop, Storage> A;
  ComputeOperatorExplicit(L, A);

  mat_lu.Factorize(A);
}


void ModeEsNonLinear2NDirichletEigenProblem::ComputePreconditioning(const VectComplex_wp& L, const VectComplex_wp& coef)
{
  cout << "ComputePreconditioning " << endl;

  Matrix<Complex_wp, Prop, Storage> A;

  Complex_wp tau = var.torsion;
  VectComplex_wp beta(L.GetM());

  A.Reallocate(nb_rows, nb_rows);

  DISP(coef.GetM());
  for (int k = 0; k < coef.GetM(); k++)
    {
      beta(k) = var.omega*L(k);
      var.ComputeFemMatricesNonLin2N(L(k), vec_MhB, vec_KhB, vec_ChB, vec_RjhB, vec_RihB, vec_SjhB, vec_SihB, vec_DhB);
      for (int ref = 0; ref < index_used.GetM(); ref++)
        if (index_used(ref))
          {
            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
            Complex_wp k2 = - eps_sigma * mu_tilde;
            Complex_wp Delta = - k2 + beta(k)*beta(k);

            // on traite les lignes avec E3
            for (int i0 = 0; i0 < nb_dof_E3; i0++)
              {
                int i = DofKeptDir(i0);
                int n1 = vec_MhB(ref).GetRowSize(i);
                int n2 = vec_SjhB(ref).GetRowSize(i);
                int n3 = vec_RjhB(ref).GetRowSize(i);
                int n4 = vec_RihB(ref).GetRowSize(i);
                int n5 = vec_SihB(ref).GetRowSize(i);
                int n6 = vec_ChB(ref).GetRowSize(i);
                int n7 = vec_KhB(ref).GetRowSize(i);
                int n8 = vec_DhB(ref).GetRowSize(i);

                Vector<int> col(n1+n2+n3+n4+n5+n6+n7+n8); Vector<Complex_wp> val(n1+n2+n3+n4+n5+n6+n7+n8);
                int num = 0;

                // partie eps_sigma DBeta MhB
                for (int j0 = 0; j0 < n1; j0++)
                  {
                    int j1 = vec_MhB(ref).Index(i, j0);
                    int j = IndexNonDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = vec_MhB(ref).Value(i, j0)*eps_sigma;
                        num++;
                      }
                  }

                // partie i beta eps_sigma tau Sjh
                for (int j0 = 0; j0 < n2; j0++)
                  {
                    int j1 = vec_SjhB(ref).Index(i, j0);
                    int j = IndexNonDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = Iwp * beta(k) * eps_sigma * tau / Delta * vec_SjhB(ref).Value(i, j0);
                        num++;
                      }
                  }

                // partie - k2 tau Rjh
                for (int j0 = 0; j0 < n3; j0++)
                  {
                    int j = vec_RjhB(ref).Index(i, j0);
                    col(num) = nb_dof_E3 + j;
                    val(num) = -k2 * tau / Delta * vec_RjhB(ref).Value(i, j0);
                    num++;
                  }

                // partie - k2 tau Rih
                for (int j0 = 0; j0 < n4; j0++)
                  {
                    int j = vec_RihB(ref).Index(i, j0);
                    col(num) = nb_dof_E3 + j;
                    val(num) = -k2 * tau / Delta * vec_RihB(ref).Value(i, j0);
                    num++;
                  }

                // partie -i beta eps_sigma tau Sjh
                for (int j0 = 0; j0 < n5; j0++)
                  {
                    int j1 = vec_SihB(ref).Index(i, j0);
                    int j = IndexNonDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = -Iwp * beta(k) * eps_sigma * tau / Delta * vec_SihB(ref).Value(i, j0);
                        num++;
                      }
                  }

                // partie i beta Ch
                for (int j0 = 0; j0 < n6; j0++)
                  {
                    int j = vec_ChB(ref).Index(i, j0);
                    col(num) = nb_dof_E3 + j;
                    val(num) = Iwp * beta(k) / Delta * vec_ChB(ref).Value(i, j0);
                    num++;
                  }

                // partie eps_sigma Kh
                for (int j0 = 0; j0 < n7; j0++)
                  {
                    int j1 = vec_KhB(ref).Index(i, j0);
                    int j = IndexNonDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = eps_sigma / Delta * vec_KhB(ref).Value(i, j0);
                        num++;
                      }
                  }

                // partie - eps_sigma tau^2 Dh
                for (int j0 = 0; j0 < n8; j0++)
                  {
                    int j1 = vec_DhB(ref).Index(i, j0);
                    int j = IndexNonDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = -eps_sigma * tau * tau / Delta * vec_DhB(ref).Value(i, j0);
                        num++;
                      }
                  }

                A.AddInteractionRow(i0, num, col, val);
              }

            // on traite les lignes avec H3
            for (int i = 0; i < nb_dof_H3; i++)
              {
                int n1 = vec_MhB(ref).GetRowSize(i);
                int n2 = vec_SjhB(ref).GetRowSize(i);
                int n3 = vec_RjhB(ref).GetRowSize(i);
                int n4 = vec_RihB(ref).GetRowSize(i);
                int n5 = vec_SihB(ref).GetRowSize(i);
                int n6 = vec_ChB(ref).GetRowSize(i);
                int n7 = vec_KhB(ref).GetRowSize(i);
                int n8 = vec_DhB(ref).GetRowSize(i);

                Vector<int> col(n1+n2+n3+n4+n5+n6+n7+n8); Vector<Complex_wp> val(n1+n2+n3+n4+n5+n6+n7+n8);
                int num = 0;

                // partie eps_sigma DBeta MhB
                for (int j0 = 0; j0 < n1; j0++)
                  {
                    int j = vec_MhB(ref).Index(i, j0);
                    col(num) = nb_dof_E3 + j;
                    val(num) = vec_MhB(ref).Value(i, j0)*mu_tilde;
                    num++;
                  }

                // partie i beta mu_tilde tau Sjh
                for (int j0 = 0; j0 < n2; j0++)
                  {
                    int j = vec_SjhB(ref).Index(i, j0);
                    col(num) = nb_dof_E3 + j;
                    val(num) = Iwp * beta(k) * mu_tilde * tau / Delta * vec_SjhB(ref).Value(i, j0);
                    num++;
                  }

                // partie k2 tau Rjh
                for (int j0 = 0; j0 < n3; j0++)
                  {
                    int j1 = vec_RjhB(ref).Index(i, j0);
                    int j = IndexNonDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = k2 * tau / Delta * vec_RjhB(ref).Value(i, j0);
                        num++;
                      }
                  }

                // partie k2 tau Rih
                for (int j0 = 0; j0 < n4; j0++)
                  {
                    int j1 = vec_RihB(ref).Index(i, j0);
                    int j = IndexNonDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = k2 * tau / Delta * vec_RihB(ref).Value(i, j0);
                        num++;
                      }
                  }

                // partie -i beta mu_tilde tau Sih
                for (int j0 = 0; j0 < n5; j0++)
                  {
                    int j = vec_SihB(ref).Index(i, j0);
                    col(num) = nb_dof_E3 + j;
                    val(num) = - Iwp * beta(k) * mu_tilde * tau / Delta * vec_SihB(ref).Value(i, j0);
                    num++;
                  }

                // partie -i beta Ch
                for (int j0 = 0; j0 < n6; j0++)
                  {
                    int j1 = vec_ChB(ref).Index(i, j0);
                    int j = IndexNonDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = -Iwp * beta(k) / Delta * vec_ChB(ref).Value(i, j0);
                        num++;
                      }
                  }

                // partie mu_tilde Kh
                for (int j0 = 0; j0 < n7; j0++)
                  {
                    int j = vec_KhB(ref).Index(i, j0);
                    col(num) = nb_dof_E3 + j;
                    val(num) = mu_tilde / Delta * vec_KhB(ref).Value(i, j0);
                    num++;
                  }

                // partie - mu_tilde tau^2 Dh
                for (int j0 = 0; j0 < n8; j0++)
                  {
                    int j = vec_DhB(ref).Index(i, j0);
                    col(num) = nb_dof_E3 + j;
                    val(num) = - mu_tilde * tau * tau / Delta * vec_DhB(ref).Value(i, j0);
                    num++;
                  }

                A.AddInteractionRow(nb_dof_E3 + i, num, col, val);
              }
          }
    }

  mat_lu.Factorize(A);
}

void ModeEsNonLinear2NDirichletEigenProblem::ComputeSplitPreconditioning(const Vector<int>& numL, const VectComplex_wp& coef)
{
  cout << "ComputeSplitPreconditioning not implemented." << endl;
  abort();
}

// application du preconditionneur Y = T(L)^{-1} X
void ModeEsNonLinear2NDirichletEigenProblem::ApplyPreconditioning(const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
{
  cout << "ApplyPreconditioning." << endl;
  // cout << "Pas encore bien écrit : ApplyPreconditioning" << endl;
  // abort();
  //cout << "appel apply" << Norm2(X) << endl;
  // if (precond)
  //   {
  //     DISP(currentL);
  //     VectComplex_wp L(1); L(0) = currentL;
  //     VectComplex_wp coef(1); coef(0) = 1.0;
  //     ComputePreconditioning(L, coef);
  //   }
  // else
  //   {
  //     precond = true;
  //   }
  mat_lu.Solve(Y, X);
}


void ModeEsNonLinear2NDirichletEigenProblem
::ComputeOperatorSplitExplicit(int num, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "ComputeOperatorSplitExplicit not implemented." << endl;
  abort();
}

void ModeEsNonLinear2NDirichletEigenProblem
::MltOperatorSplit(int num, const SeldonTranspose& trans, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  cout << "MltOperatorSplit not implemented." << endl;
  abort();
}
