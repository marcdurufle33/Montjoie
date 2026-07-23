// initialisation du probleme
ModeEsNonLinearEigenProblem
::ModeEsNonLinearEigenProblem(int N, ModeEs_Solver& var0,
                              Vector<Matrix<Complexe, Prop, Storage> >& Mh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ch,
                              Vector<Matrix<Complexe, Prop, Storage> >& Dh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Th,
                              Vector<Matrix<Complexe, Prop, Storage> >& Uh,
                              Vector<Matrix<Complexe, Symmetric, DiagonalRow> >& DTau0,
                              Vector<Matrix<Complexe, Symmetric, DiagonalRow> >& ktr2,
                              Vector<int>& DofDir, Vector<int>& IndexDir, bool use_split)
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

  // twisted = (var.torsion != 0.0);
  // DISP(twisted);
  index_used.Reallocate(var.ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < var.mesh.GetNbElt(); i++)
    index_used(var.mesh.Element(i).GetReference()) = true;

  vec_twisted.Reallocate(var.ref_tau.GetM());
  vec_twisted.Zero();
  for (int ref = 0 ; ref < var.ref_tau.GetM() ; ref++)
    {
      vec_twisted(ref) = (var.ref_tau(ref) != 0.0);
      // vec_twisted(ref) = true;
    }


  this->Init(nb_rows);
  DISP(nb_rows); DISP(nb_dof_Es); DISP(nb_dof_Hs);

  // on calcule les poles du pb aux valeurs propres
  VectComplex_wp s;
  split_reference.Clear();
  this->numer_pol_split.Clear();
  this->denom_pol_split.Clear();
  if (use_split)
    {
      Vector<Complexe> num, denom;
      A1_sum.Reallocate(nb_dof_Hs, nb_dof_Hs);
      A1_sumH.Reallocate(nb_dof_Hs, nb_dof_Hs);

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
        Complex_wp z = sqrt(-k2) / var.omega;
        Real_wp tau = var.ref_tau(ref);
        Matrix<Complexe, Prop, Storage> temp;
        s.PushBack(z); s.PushBack(-z);

        if (use_split)
          {

            cout << "Remplissage de f_i" << endl;
            // remplissage de la matrice dont la fonction f_i est 1
            Add(eps_sigma, vec_Mh(ref), A1_sum);
            Add(mu_tilde, vec_Mh(ref), A1_sumH);
            // if (twisted)
            //   {
            //     Add(4 * k2 * tau * tau * eps_sigma, vec_Mh(ref), A1_sum);
            //     Add(4 * k2 * tau * tau * mu_tilde, vec_Mh(ref), A1_sumH);
            //     Add(2 * k2 * tau * tau * eps_sigma, vec_Th(ref), A1_sum);
            //     Add(2 * k2 * tau * tau * mu_tilde, vec_Th(ref), A1_sumH);
            //     Mult(vec_Kh(ref), vec_ktrSquared(ref), temp);
            //     Add(-eps_sigma, temp, A1_sum);
            //     Add(-mu_tilde, temp, A1_sumH);
            //     Mult(vec_Dh(ref), vec_ktrSquared(ref), temp);
            //     Add(-eps_sigma * tau * tau, temp, A1_sum);
            //     Add(-mu_tilde * tau * tau, temp, A1_sumH);
            //
            //   }


            Vector<Complex_wp> num(1), denom(3);
            denom(0) = var.omega*var.omega;
            denom(1) = 0.0;
            denom(2) = -k2;

            string name;

            // f_2(beta) = 1/D(beta)
            num(0) = 1.0;
            split_reference.PushBack(ref);
            this->numer_pol_split.PushBack(num);
            this->denom_pol_split.PushBack(denom);
            name = "A2_" + to_str(ref) + ".mtx";
            DISP(this->numer_pol_split.GetM());
            file_out << " 2 " << real(num(0)) << " " << imag(num(0)) << " " << var.omega << " "
                     << real(denom(0)) << " " << imag(denom(0)) << " " << name << endl;
            ComputeOperatorSplitExplicit(this->numer_pol_split.GetM()-1, A);
            WriteMatrixMarket(A, name);

            // f_3(beta) = beta/D(beta)
            num.Reallocate(2);
            num(1) = 0.0; num(0) = var.omega;
            split_reference.PushBack(ref);
            this->numer_pol_split.PushBack(num);
            this->denom_pol_split.PushBack(denom);
            name = "A3_" + to_str(ref) + ".mtx";
            file_out << " 3 " << real(num(0)) << " " << imag(num(0)) << " " << var.omega << " "
                     << real(denom(0)) << " " << imag(denom(0)) << " " << name << endl;
            ComputeOperatorSplitExplicit(this->numer_pol_split.GetM()-1, A);
            WriteMatrixMarket(A, name);

            if (vec_twisted(ref))
              {
                // f_4(beta) = beta/D(beta)^2
                denom.Reallocate(5);
                denom(4) = k2*k2; denom(3) = 0.0; denom(2) = - 2*k2*var.omega*var.omega; denom(1) = 0.0; denom(0) = var.omega*var.omega*var.omega*var.omega;
                split_reference.PushBack(ref);
                this->numer_pol_split.PushBack(num);
                this->denom_pol_split.PushBack(denom);
                name = "A4_" + to_str(ref) + ".mtx";
                file_out << " 4 " << real(num(0)) << " " << imag(num(0)) << " " << var.omega << " "
                         << real(denom(0)) << " " << imag(denom(0)) << " " << name << endl;
                ComputeOperatorSplitExplicit(this->numer_pol_split.GetM()-1, A);
                WriteMatrixMarket(A, name);

                // f_5(beta) = 1/D(beta)^2
                num.Reallocate(1);
                num(0) = 1.0;
                split_reference.PushBack(ref);
                this->numer_pol_split.PushBack(num);
                this->denom_pol_split.PushBack(denom);
                name = "A5_" + to_str(ref) + ".mtx";
                file_out << " 5 " << real(num(0)) << " " << imag(num(0)) << " " << var.omega << " "
                         << real(denom(0)) << " " << imag(denom(0)) << " " << name << endl;
                ComputeOperatorSplitExplicit(this->numer_pol_split.GetM()-1, A);
                WriteMatrixMarket(A, name);
              }
          }
      }
    DISP(split_reference.GetM());
    cout << "Matériaux torsadés" << endl;
    vec_Inv_D_tauBeta.Reallocate(vec_D_tau0.GetM());
    vec_ChgmtVariables.Reallocate(vec_D_tau0.GetM());
    if (vec_D_tau0(ref).GetM() > 0)
      {
        vec_Inv_D_tauBeta(ref).Reallocate(N, N);
        vec_ChgmtVariables(ref).Reallocate(N, N);
      }
  }

  if (use_split)
    {
      file_out << "1 1.0 0.0 " << var.omega << " 1.0 0.0 A1.mtx" << endl;
      ComputeOperatorSplitExplicit(0, A);
      WriteMatrixMarket(A, "A1.mtx");
    }
  file_out.close();

  this->nb_split_matrix = this->numer_pol_split.GetM();
  DISP(this->nb_split_matrix);
  //DISP(this->nb_split_matrix);

  // test for the split formulation
  /* Complex_wp L = Complex_wp(1.6, 1.2);

  DistributedMatrix<Complexe, Prop, Storage> B;
  B.Reallocate(nb_rows, nb_rows);
  for (int i = 0; i < this->nb_split_matrix; i++)
    {
      Complex_wp numer(1, 0), denom(1, 0);

      if (this->numer_pol_split(i).GetM() > 0)
        {
          int n = 0;
          numer = this->numer_pol_split(i)(n);
          while (n < this->numer_pol_split(i).GetM()-1)
            {
              n++;
              numer = L*numer + this->numer_pol_split(i)(n);
            }
        }

      if (this->denom_pol_split(i).GetM() > 0)
        {
          int n = 0;
          denom = this->denom_pol_split(i)(n);
          while (n < this->denom_pol_split(i).GetM()-1)
            {
              n++;
              denom = L*denom + this->denom_pol_split(i)(n);
            }
        }

      Complex_wp coef = numer / denom;
      ComputeOperatorSplitExplicit(i, A);
      Add(coef, A, B);
    }

  B.WriteText("As.dat");
  ComputeOperatorExplicit(L, A);
  A.WriteText("A.dat");
  exit(0);*/

  this->SetSingularities(s); //DISP(s);
}

void ModeEsNonLinearEigenProblem::RemoveDirichletRowCol(const Matrix<Complexe, Prop, Storage>& A,
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

void ModeEsNonLinearEigenProblem::ExtractDirichlet(const VectComplex_wp& Xdir, VectComplex_wp& X)
{
  X.Zero();
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    X(DofKeptDir(i)) = Xdir(i);

  for (int i = 0; i < nb_dof_Hs; i++)
    X(nb_dof_Hs + i) = Xdir(nb_dof_Es + i);
}

void ModeEsNonLinearEigenProblem::ExpandDirichlet(const VectComplex_wp& X, VectComplex_wp& Xdir)
{
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    Xdir(i) = X(DofKeptDir(i));

  for (int i = 0; i < nb_dof_Hs; i++)
    Xdir(nb_dof_Es + i) = X(nb_dof_Hs + i);
}

int ModeEsNonLinearEigenProblem::GetMatrixIndex(int num)
{
  bool stop = false;
  for (int ref = 0 ; ref < index_used.GetM() ; ref++)
    {
      if (index_used(ref) && !stop)
        {
          if (vec_twisted(ref))
            {
              if (num >= 4)
                {
                  num -= 4;
                }
              else
                {
                  stop = true;
                }
            }
          else if (!vec_twisted(ref))
            {
              if (num >= 2)
                {
                  num -= 2;
                }
              else
                {
                  stop = true;
                }
            }
        }
    }
    return num;
}

void ModeEsNonLinearEigenProblem::ComputeOperator(const Complex_wp& L)
{
}

void ModeEsNonLinearEigenProblem::Mult(const Matrix<Complexe, Prop, Storage>& A, const  Matrix<Complexe, Symmetric, DiagonalRow>& D, Matrix<Complexe, Prop, Storage>& M)
{
  int N = nb_dof_Hs;
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

void ModeEsNonLinearEigenProblem::ComputeInvDBetaMatrix(const Complex_wp& L)
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

void ModeEsNonLinearEigenProblem::ComputeChgmtVariable(const Complex_wp& L)
{
  int N = nb_dof_Hs;
  Complex_wp beta = L*var.omega;
  for (int ref = 0; ref < vec_D_tau0.GetM(); ref++)
    if (vec_D_tau0(ref).GetM() > 0)
      {
        Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
        Complex_wp mu_tilde = - Iwp * var.omega * var.ref_mu(ref);
        Complex_wp k2 = Iwp*var.omega*eps_sigma*var.ref_mu(ref);
        for (int i = 0 ; i < N ; i++)
        {
          vec_ChgmtVariables(ref).Get(i,i) = (-k2 + beta * beta)/(vec_D_tau0(ref)(i,i) + beta * beta);
        }
      }
}

// calcul de Y = T(L) X
void ModeEsNonLinearEigenProblem::MltOperator(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  if (isinf(abs(L)) || isnan(abs(L)))
    {
      DISP(L);
      abort();
    }

  Complex_wp one(1, 0);
  Complex_wp beta = L*var.omega;

  Vector<Complex_wp> X(nb_dof_all), Y(nb_dof_all);
  ExtractDirichlet(Xdir, X);
  Y.Zero();

  // calcul des matrices diagonales pour les "variables intermédiaires"
  ComputeInvDBetaMatrix(L);

  // le vecteur calculé est décomposé en deux : Es, puis Hs
  int N = nb_dof_Hs;
  Vector<Complex_wp> Es(N), Hs(N), Y1(N), Y2(N);
  Y1.Zero(); Y2.Zero();
  //cout << "Je passe à : MltOperator" << endl;
  for (int i = 0 ; i < N; i++)
    {
      Es(i) = X(i);
      // V1(i) = Es(i) / (vec_D_tau0(ref).Get(i,i) + beta*beta);
      // V2(i) = V1(i) / (vec_D_tau0(ref).Get(i,i) + beta*beta);
      Hs(i) = X(N + i);
      // W1(i) = Hs(i) / (vec_D_tau0(ref).Get(i,i) + beta*beta);
      // W2(i) = W1(i) / (vec_D_tau0(ref).Get(i,i) + beta*beta);
    }

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    if (vec_Kh(ref).GetM() > 0)
      {
        Vector<Complex_wp> F1(N), F2(N), G1(N), G2(N);
        Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
        Complex_wp mu_tilde = - Iwp * var.omega * var.ref_mu(ref);
        Complex_wp k2 = Iwp*var.omega*eps_sigma*var.ref_mu(ref);
        Complex_wp D0 = -k2 + beta*beta;

        // Mlt(vec_Inv_D_tauBeta(ref), Es, V1);
        // Mlt(vec_Inv_D_tauBeta(ref), V1, V2);
        // Mlt(vec_Inv_D_tauBeta(ref), Hs, W1);
        // Mlt(vec_Inv_D_tauBeta(ref), W1, W2);

        Mlt(vec_ChgmtVariables(ref), Es, F1);
        Mlt(vec_ChgmtVariables(ref), F1, F2);
        Mlt(vec_ChgmtVariables(ref), Hs, G1);
        Mlt(vec_ChgmtVariables(ref), G1, G2);

        // pb sans torsion
        // MltAdd(eps_sigma, vec_Mh(ref), Es, one, Y1);
        // MltAdd(- Iwp * var.omega * var.ref_mu(ref), vec_Mh(ref), Hs, one, Y2);
        // MltAdd(eps_sigma/D0, vec_Kh(ref), Es, one, Y1);
        // MltAdd((- Iwp * var.omega * var.ref_mu(ref))/D0, vec_Kh(ref), Hs, one, Y2);
        // MltAdd(Iwp * beta / D0, vec_Ch(ref), Hs, one, Y1);
        // MltAdd(-Iwp * beta / D0, vec_Ch(ref), Es, one, Y2);

        // pb avec torsion
        // MltAdd(eps_sigma * D0, vec_Mh(ref), V1, one, Y1);
        // MltAdd(- 2.0 * Iwp * beta * var.torsion * eps_sigma , vec_Uh(ref), V1, one, Y1);
        // MltAdd(2.0 * k2 * var.torsion * D0, vec_Mh(ref), W2, one, Y1);
        // MltAdd(-2.0 * Iwp * beta * k2 * var.torsion * var.torsion, vec_Uh(ref), W2, one, Y1);
        // MltAdd(4.0 * k2 * var.torsion * var.torsion * eps_sigma, vec_Mh(ref), V2, one, Y1);
        // MltAdd(2.0 * k2 * var.torsion * var.torsion * eps_sigma, vec_Th(ref), V2, one, Y1);
        // MltAdd(eps_sigma, vec_Kh(ref), V1, one, Y1);
        // MltAdd(eps_sigma * var.torsion * var.torsion, vec_Dh(ref), V1, one, Y1);
        // MltAdd(Iwp * beta, vec_Ch(ref), W1, one, Y1);
        //
        // MltAdd(mu_tilde * D0, vec_Mh(ref), W1, one, Y2);
        // MltAdd(- 2.0 * Iwp * beta * var.torsion * mu_tilde , vec_Uh(ref), W1, one, Y2);
        // MltAdd(- 2.0 * k2 * var.torsion * D0, vec_Mh(ref), V2, one, Y2);
        // MltAdd(2.0 * Iwp * beta * k2 * var.torsion * var.torsion, vec_Uh(ref), V2, one, Y2);
        // MltAdd(4.0 * k2 * var.torsion * var.torsion * mu_tilde, vec_Mh(ref), W2, one, Y2);
        // MltAdd(2.0 * k2 * var.torsion * var.torsion * mu_tilde, vec_Th(ref), W2, one, Y2);
        // MltAdd(mu_tilde, vec_Kh(ref), W1, one, Y2);
        // MltAdd(mu_tilde * var.torsion * var.torsion, vec_Dh(ref), W1, one, Y2);
        // MltAdd(-Iwp * beta, vec_Ch(ref), V1, one, Y2);

        MltAdd(eps_sigma, vec_Mh(ref), F1, one, Y1);
        MltAdd(- 2.0 * Iwp * beta * var.ref_tau(ref) * eps_sigma/D0 , vec_Uh(ref), F1, one, Y1);
        MltAdd(2.0 * k2 * var.ref_tau(ref) / D0, vec_Mh(ref), G2, one, Y1);
        MltAdd(-2.0 * Iwp * beta * k2 * var.ref_tau(ref) * var.ref_tau(ref)/ (D0*D0), vec_Uh(ref), G2, one, Y1);
        MltAdd(4.0 * k2 * var.ref_tau(ref) * var.ref_tau(ref) * eps_sigma/(D0*D0), vec_Mh(ref), F2, one, Y1);
        MltAdd(2.0 * k2 * var.ref_tau(ref) * var.ref_tau(ref) * eps_sigma/(D0*D0), vec_Th(ref), F2, one, Y1);
        MltAdd(eps_sigma/D0, vec_Kh(ref), F1, one, Y1);
        MltAdd(eps_sigma * var.ref_tau(ref) * var.ref_tau(ref)/D0, vec_Dh(ref), F1, one, Y1);
        MltAdd(Iwp * beta/D0, vec_Ch(ref), G1, one, Y1);

        MltAdd(mu_tilde, vec_Mh(ref), G1, one, Y2);
        MltAdd(- 2.0 * Iwp * beta * var.ref_tau(ref) * mu_tilde/D0 , vec_Uh(ref), G1, one, Y2);
        MltAdd(- 2.0 * k2 * var.ref_tau(ref) / D0, vec_Mh(ref), F2, one, Y2);
        MltAdd(2.0 * Iwp * beta * k2 * var.ref_tau(ref) * var.ref_tau(ref)/(D0*D0), vec_Uh(ref), F2, one, Y2);
        MltAdd(4.0 * k2 * var.ref_tau(ref) * var.ref_tau(ref) * mu_tilde/(D0*D0), vec_Mh(ref), G2, one, Y2);
        MltAdd(2.0 * k2 * var.ref_tau(ref) * var.ref_tau(ref) * mu_tilde/(D0*D0), vec_Th(ref), G2, one, Y2);
        MltAdd(mu_tilde/D0, vec_Kh(ref), G1, one, Y2);
        MltAdd(mu_tilde * var.ref_tau(ref) * var.ref_tau(ref)/D0, vec_Dh(ref), G1, one, Y2);
        MltAdd(-Iwp * beta/D0, vec_Ch(ref), F1, one, Y2);
      }

  for (int i = 0 ; i < N ; i++)
    {
      Y(i) = Y1(i);
      Y(N + i) = Y2(i);
    }

  ExpandDirichlet(Y, Ydir);
}

// calcul explicite de T(L) (utilise si on calcule explicitement les matrices)
void ModeEsNonLinearEigenProblem
::ComputeOperatorExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  A.Reallocate(nb_rows, nb_rows);

  /*DISP(L);
  VectComplex_wp Ones(nb_rows), Aones(nb_rows);
  ofstream file_out("B.dat"); file_out.precision(15);
  for (int i = 0; i < nb_rows; i++)
    {
      Ones.Zero(); Ones(i) = 1.0;
      MltOperator(L, SeldonNoTrans, Ones, Aones);
      for (int j = 0; j < nb_rows; j++)
        if (abs(Aones(j)) > 1e-12)
          file_out << j+1 << " " << i+1 << " " << Aones(j) << "\n";
    }

    file_out.close(); */

  Complex_wp beta = L*var.omega;

  // calcul des matrices diagonales pour les "variables intermédiaires"
  ComputeInvDBetaMatrix(L);

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    if (vec_Mh(ref).GetM() > 0)
      {
        Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
        Complex_wp mu_tilde = - Iwp * var.omega * var.ref_mu(ref);
        Complex_wp k2 = - eps_sigma * mu_tilde;
        Complex_wp D0 = - k2 + beta*beta;
        Real_wp tau = var.ref_tau(ref);

        // on traite les lignes avec Es
        for (int i0 = 0; i0 < nb_dof_Es; i0++)
          {
            int i = DofKeptDir(i0);
            int n1 = vec_Mh(ref).GetRowSize(i);
            int n2 = vec_Kh(ref).GetRowSize(i);
            int n3 = vec_Ch(ref).GetRowSize(i);
            int n4 = vec_Uh(ref).GetRowSize(i);
            int n5 = vec_Th(ref).GetRowSize(i);
            int n6 = vec_Dh(ref).GetRowSize(i);

            Vector<int> col(n1+n2+n3+n4+n5+n6); Vector<Complex_wp> val(n1+n2+n3+n4+n5+n6);
            int num = 0;

            // parties avec Mh
            for (int j0 = 0; j0 < n1; j0++)
              {
                int j1 = vec_Mh(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j1,j1);
                if (j >= 0)
                  {
                    col(num) = j;
                    val(num) = 0;
                    // partie eps_sigma * D0 Mh V1
                    val(num) += eps_sigma * D0 * vec_Mh(ref).Value(i, j0) * InvDtauBeta;
                    // partie 4 k2 tau^2 eps_sigma Mh V2
                    val(num) += 4.0 * k2 * tau * tau * eps_sigma * vec_Mh(ref).Value(i, j0) * InvDtauBeta * InvDtauBeta;
                    num++;
                  }

                col(num) = nb_dof_Es + j1;
                // partie 2 k2 torsion D0 Mh W2
                val(num) = 2.0 * k2 * tau * D0 * vec_Mh(ref).Value(i, j0) * InvDtauBeta * InvDtauBeta;
                num++;
              }

            A.AddInteractionRow(i0, num, col, val);

            num = 0;
            // parties avec Kh
            for (int j0 = 0; j0 < n2; j0++)
              {
                int j1 = vec_Kh(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                if (j >= 0)
                  {
                    Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j1,j1);
                    col(num) = j;
                    // eps_sigma Kh V1
                    val(num) = eps_sigma * vec_Kh(ref).Value(i, j0) * InvDtauBeta;
                    num++;
                  }
              }

            // partie avec Ch
            for (int j0 = 0; j0 < n3; j0++)
              {
                int j = vec_Ch(ref).Index(i, j0);
                Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j,j);
                col(num) = nb_dof_Es + j;
                // partie Iwp * beta Ch W1
                val(num) = Iwp * beta * vec_Ch(ref).Value(i, j0) * InvDtauBeta;
                num++;
              }

            A.AddInteractionRow(i0, num, col, val);

            num = 0;
            // parties avec Uh
            for (int j0 = 0; j0 < n4; j0++)
              {
                int j1 = vec_Uh(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j1,j1);
                if (j >= 0)
                  {
                    col(num) = j;
                    // partie - 2 Iwp beta tau eps_sigma Uh V1
                    val(num) = -2.0 * Iwp * beta * tau * eps_sigma * vec_Uh(ref).Value(i, j0) * InvDtauBeta;
                    num++;
                  }

                col(num) = nb_dof_Es + j1;
                // partie -2 Iwp beta k2 tau^2 Uh W2
                val(num) = - 2.0 * Iwp * beta * k2 * tau * tau * vec_Uh(ref).Value(i, j0) * InvDtauBeta * InvDtauBeta;
                num++;
              }

            A.AddInteractionRow(i0, num, col, val);

            num = 0;
            // parties avec Th
            for (int j0 = 0; j0 < n5; j0++)
              {
                int j1 = vec_Th(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                if (j >= 0)
                  {
                    Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j1,j1);
                    col(num) = j;
                    // partie 2 k2 tau^2 eps_sigma Th V2
                    val(num) = 2.0 * k2 * tau * tau * eps_sigma * vec_Th(ref).Value(i, j0) * InvDtauBeta * InvDtauBeta;
                    num++;
                  }
              }

            A.AddInteractionRow(i0, num, col, val);

            num = 0;
            // parties avec Dh
            for (int j0 = 0; j0 < n6; j0++)
              {
                int j1 = vec_Dh(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                if (j >= 0)
                  {
                    Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j1,j1);
                    col(num) = j;
                    // partie eps_sigma tau^2 Dh V1
                    val(num) = tau * tau * eps_sigma * vec_Dh(ref).Value(i, j0) * InvDtauBeta;
                    num++;
                  }
              }

            A.AddInteractionRow(i0, num, col, val);
          }

        // on traite les lignes avec Hs
        for (int i = 0; i < nb_dof_Hs; i++)
          {
            int n1 = vec_Mh(ref).GetRowSize(i);
            int n2 = vec_Kh(ref).GetRowSize(i);
            int n3 = vec_Ch(ref).GetRowSize(i);
            int n4 = vec_Uh(ref).GetRowSize(i);
            int n5 = vec_Th(ref).GetRowSize(i);
            int n6 = vec_Dh(ref).GetRowSize(i);

            Vector<int> col(n1+n2+n3+n4+n5+n6); Vector<Complex_wp> val(n1+n2+n3+n4+n5+n6);

            int num = 0;

            // parties avec Mh
            for (int j0 = 0; j0 < n1; j0++)
              {
                int j1 = vec_Mh(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j1,j1);
                col(num) = nb_dof_Es + j1;
                val(num) = 0;
                // partie mu_tilde * D0 Mh W1
                val(num) += mu_tilde * D0 * vec_Mh(ref).Value(i, j0) * InvDtauBeta;
                // partie 4 k2 tau^2 mu_tilde Mh W2
                val(num) += 4.0 * k2 * tau * tau * mu_tilde * vec_Mh(ref).Value(i, j0) * InvDtauBeta * InvDtauBeta;
                num++;

                if (j >= 0)
                  {
                    col(num) = j;
                    // partie -2 k2 torsion D0 Mh V2
                    val(num) = -2.0 * k2 * tau * D0 * vec_Mh(ref).Value(i, j0) * InvDtauBeta * InvDtauBeta;
                    num++;
                  }
              }

            A.AddInteractionRow(nb_dof_Es + i, num, col, val);

            num = 0;

            // parties avec Kh
            for (int j0 = 0; j0 < n2; j0++)
              {
                int j = vec_Kh(ref).Index(i, j0);
                Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j,j);
                col(num) = nb_dof_Es + j;
                // mu_tilde Kh W1
                val(num) = mu_tilde * vec_Kh(ref).Value(i, j0) * InvDtauBeta;
                num++;
              }

            // partie avec Ch
            for (int j0 = 0; j0 < n3; j0++)
              {
                int j1 = vec_Ch(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j1,j1);
                if (j >= 0)
                  {
                    col(num) = j;
                    // partie -Iwp * beta Ch V1
                    val(num) = -Iwp * beta * vec_Ch(ref).Value(i, j0) * InvDtauBeta;
                    num++;
                  }
              }

            A.AddInteractionRow(nb_dof_Es + i, num, col, val);

            num = 0;
            // parties avec Uh
            for (int j0 = 0; j0 < n4; j0++)
              {
                int j1 = vec_Uh(ref).Index(i, j0);
                int j = IndexDirichlet(j1);

                Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j1,j1);
                col(num) = nb_dof_Es + j1;
                // partie - 2 Iwp beta tau mu_tilde Uh W1
                val(num) = -2.0 * Iwp * beta * tau * mu_tilde * vec_Uh(ref).Value(i, j0) * InvDtauBeta;
                num++;

                if (j >= 0)
                  {
                    col(num) = j;
                    // partie 2 Iwp beta k2 tau^2 Uh V2
                    val(num) = 2.0 * Iwp * beta * k2 * tau * tau * vec_Uh(ref).Value(i, j0) * InvDtauBeta * InvDtauBeta;
                    num++;
                  }
              }

            A.AddInteractionRow(nb_dof_Es + i, num, col, val);

            num = 0;
            // parties avec Th
            for (int j0 = 0; j0 < n5; j0++)
              {
                int j = vec_Th(ref).Index(i, j0);
                Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j,j);
                col(num) = nb_dof_Es + j;
                // partie 2 k2 tau^2 mu_tilde Th W2
                val(num) = 2.0 * k2 * tau * tau * mu_tilde * vec_Th(ref).Value(i, j0) * InvDtauBeta * InvDtauBeta;
                num++;
              }

            A.AddInteractionRow(nb_dof_Es + i, num, col, val);

            num = 0;
            // parties avec Dh
            for (int j0 = 0; j0 < n6; j0++)
              {
                int j = vec_Dh(ref).Index(i, j0);
                Complex_wp InvDtauBeta = vec_Inv_D_tauBeta(ref).Get(j,j);
                col(num) = nb_dof_Es + j;
                // partie mu_tilde tau^2 Dh W1
                val(num) = tau * tau * mu_tilde * vec_Dh(ref).Value(i, j0) * InvDtauBeta;
                num++;
              }

            A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          }


      }
  //A.WriteText("B.dat");
  //exit(0);
}


void ModeEsNonLinearEigenProblem::ComputeJacobian(const Complex_wp& L)
{}

// calcul de Y = T'(L) X
void ModeEsNonLinearEigenProblem::MltJacobian(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  Vector<Complex_wp> X(nb_dof_all), Y(nb_dof_all);
  ExtractDirichlet(Xdir, X);
  Y.Zero();

  Complex_wp one(1, 0);
  Complex_wp beta = L*var.omega;

  ComputeInvDBetaMatrix(L);

  // le vecteur calculé est décomposé en deux : Es, puis Hs
  int N = this->nb_dof_Hs;
  Vector<Complex_wp> Es(N), Hs(N), Y1(N), Y2(N);
  Y1.Zero(); Y2.Zero();
  //cout << "Je passe à : MltJacobian" << endl;
  for (int i = 0 ; i < N ; i++)
    {
      Es(i) = X(i);
      Hs(i) = X(N + i);
    }

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    if (vec_Mh(ref).GetM() > 0)
      {

        Vector<Complex_wp> V1(N), V2(N), V2Star(N), V3(N), W1(N), W2(N), W2Star(N), W3(N);
        Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
        Complex_wp mu_tilde = - Iwp * var.omega * var.ref_mu(ref);
        Complex_wp k2 = Iwp*var.omega*eps_sigma*var.ref_mu(ref);
        Complex_wp D0 = -k2 + beta*beta;
        Complex_wp dD = 2.0*var.omega*beta;

        Mlt(vec_Inv_D_tauBeta(ref), Es, V1);
        Mlt(vec_Inv_D_tauBeta(ref), V1, V2);
        Mlt(vec_ktrSquared(ref), V2, V2Star);
        Mlt(vec_Inv_D_tauBeta(ref), V2, V3);
        Mlt(vec_Inv_D_tauBeta(ref), Hs, W1);
        Mlt(vec_Inv_D_tauBeta(ref), W1, W2);
        Mlt(vec_ktrSquared(ref), W2, W2Star);
        Mlt(vec_Inv_D_tauBeta(ref), W2, W3);

        // // pb sans torsion
        // MltAdd(-eps_sigma/(Delta*Delta)*dDelta, vec_Kh(ref), Es, one, Y1);
        // MltAdd(Iwp * var.omega * var.ref_mu(ref)/(Delta * Delta) * dDelta, vec_Kh(ref), Hs, one, Y2);
        // MltAdd((Iwp * var.omega * Delta - Iwp * beta * dDelta) / (Delta * Delta), vec_Ch(ref), Hs, one, Y1);
        // MltAdd((-Iwp * var.omega * Delta + Iwp * beta * dDelta) / (Delta * Delta), vec_Ch(ref), Es, one, Y2);

        // pb avec torsion
        MltAdd(-2.0*beta*var.omega*eps_sigma, vec_Mh(ref), V2Star, one, Y1);
        MltAdd(-2.0 * Iwp * eps_sigma * var.omega * var.ref_tau(ref), vec_Uh(ref), V1, one, Y1);
        MltAdd(4.0 * beta * beta * Iwp * eps_sigma * var.omega * var.ref_tau(ref), vec_Uh(ref), V2, one, Y1);
        MltAdd(4.0 * k2 * var.ref_tau(ref) * beta * var.omega, vec_Mh(ref), W2, one, Y1);
        MltAdd(-8.0 * k2 * var.ref_tau(ref) * beta * var.omega * D0, vec_Mh(ref), W3, one, Y1);
        MltAdd(Iwp * var.omega, vec_Ch(ref), W1, one, Y1);
        MltAdd(-2.0 * Iwp * var.omega * k2 * var.ref_tau(ref) * var.ref_tau(ref), vec_Uh(ref), W2, one, Y1);
        MltAdd(-2.0 * Iwp * beta * beta * var.omega, vec_Ch(ref), W2, one, Y1);
        MltAdd(8.0 * Iwp * beta * beta * var.omega * k2 * var.ref_tau(ref) * var.ref_tau(ref), vec_Uh(ref), W3, one, Y1);
        MltAdd(-2.0 * beta * var.omega * eps_sigma, vec_Kh(ref), V2, one, Y1);
        MltAdd(-8.0 * k2 * var.ref_tau(ref) * var.ref_tau(ref) * beta * var.omega * eps_sigma, vec_Th(ref), V3, one, Y1);
        MltAdd(-16.0 * k2 * var.ref_tau(ref) * var.ref_tau(ref) * beta * var.omega * eps_sigma, vec_Mh(ref), V3, one, Y1);
        MltAdd(-2.0 * eps_sigma * var.ref_tau(ref) * var.ref_tau(ref) * beta * var.omega, vec_Dh(ref), V2, one, Y1);


        MltAdd(-2.0*beta*var.omega*mu_tilde, vec_Mh(ref), W2Star, one, Y2);
        MltAdd(-2.0 * Iwp * mu_tilde * var.omega * var.ref_tau(ref), vec_Uh(ref), W1, one, Y2);
        MltAdd(4.0 * beta * beta * Iwp * mu_tilde * var.omega * var.ref_tau(ref), vec_Uh(ref), W2, one, Y2);
        MltAdd(-4.0 * k2 * var.ref_tau(ref) * beta * var.omega, vec_Mh(ref), V2, one, Y2);
        MltAdd(8.0 * k2 * var.ref_tau(ref) * beta * var.omega * D0, vec_Mh(ref), V3, one, Y2);
        MltAdd(-Iwp * var.omega, vec_Ch(ref), V1, one, Y2);
        MltAdd(2.0 * Iwp * var.omega * k2 * var.ref_tau(ref) * var.ref_tau(ref), vec_Uh(ref), V2, one, Y2);
        MltAdd(2.0 * Iwp * beta * beta * var.omega, vec_Ch(ref), V2, one, Y2);
        MltAdd(-8.0 * Iwp * beta * beta * var.omega * k2 * var.ref_tau(ref) * var.ref_tau(ref), vec_Uh(ref), V3, one, Y2);
        MltAdd(-2.0 * beta * var.omega * mu_tilde, vec_Kh(ref), W2, one, Y2);
        MltAdd(-8.0 * k2 * var.ref_tau(ref) * var.ref_tau(ref) * beta * var.omega * mu_tilde, vec_Th(ref), W3, one, Y2);
        MltAdd(-16.0 * k2 * var.ref_tau(ref) * var.ref_tau(ref) * beta * var.omega * mu_tilde, vec_Mh(ref), W3, one, Y2);
        MltAdd(-2.0 * mu_tilde * var.ref_tau(ref) * var.ref_tau(ref) * beta * var.omega, vec_Dh(ref), W2, one, Y2);
      }

  for (int i = 0 ; i < N ; i++)
    {
      Y(i) = Y1(i);
      Y(N + i) = Y2(i);
    }

  ExpandDirichlet(Y, Ydir);
}

// calcul explicite de T'(L)
void ModeEsNonLinearEigenProblem
::ComputeJacobianExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "Pas écrit : ComputeJacobianExplicit" << endl;
  abort();
  A.Reallocate(nb_rows, nb_rows);

  Complex_wp beta = L*var.omega;

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    if (vec_Mh(ref).GetM() > 0)
      {
        Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
        Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
        Complex_wp dDelta = 2.0*var.omega*beta;
        Complex_wp coef_11 = -eps_sigma/(Delta*Delta)*dDelta;
        Complex_wp coef_12 = (Iwp * var.omega * Delta - Iwp * beta * dDelta) / (Delta * Delta);

        // on traite les lignes avec Es
        for (int i0 = 0; i0 < nb_dof_Es; i0++)
          {
            int i = DofKeptDir(i0);
            int n2 = vec_Kh(ref).GetRowSize(i);
            int n3 = vec_Ch(ref).GetRowSize(i);

            Vector<int> col(n2+n3); Vector<Complex_wp> val(n2+n3);
            int num = 0;
            for (int j0 = 0; j0 < n2; j0++)
              {
                int j1 = vec_Kh(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                if (j >= 0)
                  {
                    col(num) = j;
                    val(num) = vec_Kh(ref).Value(i, j0)*coef_11;
                    num++;
                  }
              }

            for (int j0 = 0; j0 < n3; j0++)
              {
                int j = vec_Ch(ref).Index(i, j0);
                col(num) = nb_dof_Es + j;
                val(num) = vec_Ch(ref).Value(i, j0)*coef_12;
                num++;
              }

            A.AddInteractionRow(i0, num, col, val);
          }

        Complex_wp coef_22 = Iwp * var.omega * var.ref_mu(ref)/(Delta * Delta) * dDelta;
        Complex_wp coef_21 = (-Iwp * var.omega * Delta + Iwp * beta * dDelta) / (Delta * Delta);

        // on traite les lignes avec Hs
        for (int i = 0; i < nb_dof_Hs; i++)
          {
            int n2 = vec_Kh(ref).GetRowSize(i);
            int n3 = vec_Ch(ref).GetRowSize(i);

            Vector<int> col(n2+n3); Vector<Complex_wp> val(n2+n3);
            int num = 0;

            for (int j0 = 0; j0 < n2; j0++)
              {
                int j = vec_Kh(ref).Index(i, j0);
                col(num) = nb_dof_Es + j;
                val(num) = vec_Kh(ref).Value(i, j0)*coef_22;
                num++;
              }

            for (int j0 = 0; j0 < n3; j0++)
              {
                int j1 = vec_Ch(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                if (j >= 0)
                  {
                    col(num) = j;
                    val(num) = vec_Ch(ref).Value(i, j0)*coef_21;
                    num++;
                  }
              }

            A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          }
      }
}

void ModeEsNonLinearEigenProblem::ComputeExplicitPreconditioning(DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  mat_lu.Factorize(A);
}

// factorisation de T(L)
void ModeEsNonLinearEigenProblem::ComputePreconditioning(const Complex_wp& L)
{
  DistributedMatrix<Complex_wp, Prop, Storage> A;
  ComputeOperatorExplicit(L, A);

  mat_lu.Factorize(A);
}


void ModeEsNonLinearEigenProblem::ComputePreconditioning(const VectComplex_wp& L, const VectComplex_wp& coef)
{
  cout << "Pas encore bien écrit : Preconditionning" << endl;
  abort();
  Matrix<Complex_wp, Prop, Storage> A;

  VectComplex_wp beta(L.GetM());
  for (int k = 0; k < L.GetM(); k++)
    beta(k) = var.omega*L(k);

  A.Reallocate(nb_rows, nb_rows);

  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    if (vec_Mh(ref).GetM() > 0)
      {
        Complex_wp coef_11, coef_12, coef_22, coefKh_22, coef_21, coefM_11;
        SetComplexZero(coef_11); SetComplexZero(coef_12); SetComplexZero(coef_22);
        SetComplexZero(coefKh_22); SetComplexZero(coef_21); SetComplexZero(coefM_11);
        for (int k = 0; k < coef.GetM(); k++)
          {
            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta(k)*beta(k);
            coefM_11 += coef(k)*eps_sigma;
            coef_11 += coef(k)*eps_sigma/Delta;
            coef_12 += coef(k)*Iwp * beta(k) / Delta;
            coef_22 += - coef(k)*Iwp * var.omega * var.ref_mu(ref);
            coefKh_22 += - Iwp * coef(k)*var.omega * var.ref_mu(ref)/Delta;
            coef_21 -= coef(k)*Iwp * beta(k) / Delta;
          }

        // on traite les lignes avec Es
        for (int i0 = 0; i0 < nb_dof_Es; i0++)
          {
            int i = DofKeptDir(i0);
            int n1 = vec_Mh(ref).GetRowSize(i);
            int n2 = vec_Kh(ref).GetRowSize(i);
            int n3 = vec_Ch(ref).GetRowSize(i);

            Vector<int> col(n1+n2+n3); Vector<Complex_wp> val(n1+n2+n3);
            int num = 0;

            // partie eps_sigma Mh
            for (int j0 = 0; j0 < n1; j0++)
              {
                int j1 = vec_Mh(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                if (j >= 0)
                  {
                    col(num) = j;
                    val(num) = vec_Mh(ref).Value(i, j0)*coefM_11;
                    num++;
                  }
              }

            // partie eps_sigma/Delta Kh
            for (int j0 = 0; j0 < n2; j0++)
              {
                int j1 = vec_Kh(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                if (j >= 0)
                  {
                    col(num) = j;
                    val(num) = vec_Kh(ref).Value(i, j0)*coef_11;
                    num++;
                  }
              }

            // partie i beta/Delta Ch
            for (int j0 = 0; j0 < n3; j0++)
              {
                int j = vec_Ch(ref).Index(i, j0);
                col(num) = nb_dof_Es + j;
                val(num) = vec_Ch(ref).Value(i, j0)*coef_12;
                num++;
              }

            A.AddInteractionRow(i0, num, col, val);
          }

        // on traite les lignes avec Hs
        for (int i = 0; i < nb_dof_Hs; i++)
          {
            int n1 = vec_Mh(ref).GetRowSize(i);
            int n2 = vec_Kh(ref).GetRowSize(i);
            int n3 = vec_Ch(ref).GetRowSize(i);

            Vector<int> col(n1+n2+n3); Vector<Complex_wp> val(n1+n2+n3);
            int num = 0;

            // partie -i omega mu Mh
            for (int j0 = 0; j0 < n1; j0++)
              {
                int j = vec_Mh(ref).Index(i, j0);
                col(num) = nb_dof_Es + j;
                val(num) = vec_Mh(ref).Value(i, j0)*coef_22;
                num++;
              }

            // partie -i omega mu/Delta Kh
            for (int j0 = 0; j0 < n2; j0++)
              {
                int j = vec_Kh(ref).Index(i, j0);
                col(num) = nb_dof_Es + j;
                val(num) = vec_Kh(ref).Value(i, j0)*coefKh_22;
                num++;
              }

            // partie -i beta/Delta Ch
            for (int j0 = 0; j0 < n3; j0++)
              {
                int j1 = vec_Ch(ref).Index(i, j0);
                int j = IndexDirichlet(j1);
                if (j >= 0)
                  {
                    col(num) = j;
                    val(num) = vec_Ch(ref).Value(i, j0)*coef_21;
                    num++;
                  }
              }

            A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          }
      }

  mat_lu.Factorize(A);
}

void ModeEsNonLinearEigenProblem::ComputeSplitPreconditioning(const Vector<int>& numL, const VectComplex_wp& coef)
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
void ModeEsNonLinearEigenProblem::ApplyPreconditioning(const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
{
  //cout << "appel apply" << Norm2(X) << endl;
  mat_lu.Solve(Y, X);
}


void ModeEsNonLinearEigenProblem
::ComputeOperatorSplitExplicit(int num, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  int ref = this->split_reference(num-1);
  if (vec_twisted(ref))
    {
      if (num == 0)
        {
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              int i0 = IndexDirichlet(i);
              int nb_val = 0;
              for (int j = 0; j < A1_sum.GetRowSize(i); j++)
                if (IndexDirichlet(A1_sum.Index(i, j)) >= 0)
                  nb_val++;

              if (i0 >= 0)
                {
                  A.ReallocateRow(i0, nb_val);
                  nb_val = 0;
                  for (int j = 0; j < A1_sum.GetRowSize(i); j++)
                    if (IndexDirichlet(A1_sum.Index(i, j)) >= 0)
                      {
                        A.Index(i0, nb_val) = IndexDirichlet(A1_sum.Index(i, j));
                        A.Value(i0, nb_val) = A1_sum.Value(i, j);
                        nb_val++;
                      }
                }

              nb_val = A1_sumH.GetRowSize(i);
              A.ReallocateRow(nb_dof_Es + i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(nb_dof_Es + i, j) = nb_dof_Es + A1_sumH.Index(i, j);
                  A.Value(nb_dof_Es + i, j) = A1_sumH.Value(i, j);
                }
            }
        }
      else if (GetMatrixIndex(num) == 1)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = - eps_sigma * mu_tilde;
          Real_wp tau = var.ref_tau(ref);
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);


          // on traite les lignes avec V2
          for (int iDir = 0; iDir < nb_dof_Es; iDir++)
            {
              int i = DofKeptDir(iDir); // i = ligne globale
              int nM_Hs = vec_Mh(ref).GetRowSize(i);
              int nK_Hs = vec_Kh(ref).GetRowSize(i);
              int nD_Hs = vec_Dh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;
              int num = 0;
              col.Reallocate(2*nM_Hs); val.Reallocate(2*nM_Hs);
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
                      val(num) = - eps_sigma * vec_Mh(ref).Value(i, j) * ktr2;
                      num++;
                    }

                  col(num) = nb_dof_Es + j_Glob;
                  // partie 2 k2 tau Mh W2
                  val(num) = 2.0 * k2 * tau * vec_Mh(ref).Value(i, j);
                  num++;
                }

              A.AddInteractionRow(iDir, num, col, val);

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
                      val(num) = eps_sigma * vec_Kh(ref).Value(i, j);
                      num++;
                    }
                }

              A.AddInteractionRow(iDir, num, col, val);

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
                      val(num) = eps_sigma * tau * tau * vec_Dh(ref).Value(i, j);
                      num++;
                    }
                }

              A.AddInteractionRow(iDir, num, col, val);

            }

          // on traite les lignes avec W2
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              int nM_Hs = vec_Mh(ref).GetRowSize(i);
              int nK_Hs = vec_Kh(ref).GetRowSize(i);
              int nD_Hs = vec_Dh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;

              int num = 0;
              col.Reallocate(2*nM_Hs); val.Reallocate(2*nM_Hs);
              // parties avec Uh
              for (int j = 0; j < nM_Hs; j++)
                {
                  int j_Glob = vec_Mh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);

                  col(num) = nb_dof_Es + j_Glob;
                  // partie - mu_tilde Mh W2t
                  val(num) = - mu_tilde * vec_Mh(ref).Value(i, j) * ktr2;
                  num++;

                  if (jDir >= 0)
                    {
                      col(num) = jDir;
                      // partie - 2 k2 tau Mh, V2
                      val(num) = - 2.0 * k2 * tau * vec_Mh(ref).Value(i, j);
                      num++;
                    }
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);

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
                  val(num) = mu_tilde * vec_Kh(ref).Value(i, j);
                  num++;
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);

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
                  val(num) = mu_tilde * tau * tau * vec_Dh(ref).Value(i, j);
                  num++;
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);
            }






          // int ref = this->split_reference(num/4);
          // Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          // Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          // Complex_wp k2 = - eps_sigma * mu_tilde;
          // Real_wp tau = var.torsion;
          // A.Clear();
          // A.Reallocate(nb_rows, nb_rows);
          //
          // // on traite les lignes avec V2
          // for (int iDir = 0; iDir < nb_dof_Es; iDir++)
          //   {
          //     int i = DofKeptDir(iDir); // i = ligne globale
          //     int nC_Hs = vec_Ch(ref).GetRowSize(i);
          //     int nU_Hs = vec_Uh(ref).GetRowSize(i);
          //
          //     Vector<int> col; Vector<Complex_wp> val;
          //     int num = 0;
          //     col.Reallocate(nC_Hs); val.Reallocate(nC_Hs);
          //     // partie avec Ch
          //     for (int j = 0; j < nC_Hs; j++)
          //       {
          //         int j_Glob = vec_Ch(ref).Index(i, j);
          //         Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
          //         col(num) = nb_dof_Es + j_Glob;
          //         // partie - Iwp Ch W2t
          //         val(num) = - Iwp * vec_Ch(ref).Value(i, j) * ktr2;
          //         num++;
          //       }
          //
          //     A.AddInteractionRow(iDir, num, col, val);
          //
          //     num = 0;
          //     col.Reallocate(2*nU_Hs); val.Reallocate(2*nU_Hs);
          //     // parties avec Uh
          //     for (int j = 0; j < nU_Hs; j++)
          //       {
          //         int j_Glob = vec_Uh(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //         Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
          //         if (jDir >= 0)
          //           {
          //             col(num) = jDir;
          //             // partie 2 i tau eps_sigma Uh V2t
          //             val(num) = 2.0 * Iwp * tau * eps_sigma * vec_Uh(ref).Value(i, j) * ktr2;
          //             num++;
          //           }
          //
          //         col(num) = nb_dof_Es + j_Glob;
          //         // partie -2 Iwp k2 tau^2 Uh W2
          //         val(num) = - 2.0 * Iwp * k2 * tau * tau * vec_Uh(ref).Value(i, j);
          //         num++;
          //       }
          //
          //     A.AddInteractionRow(iDir, num, col, val);
          //   }
          //
          // // on traite les lignes avec W2
          // for (int i = 0; i < nb_dof_Hs; i++)
          //   {
          //     int nC_Hs = vec_Ch(ref).GetRowSize(i);
          //     int nU_Hs = vec_Uh(ref).GetRowSize(i);
          //
          //     Vector<int> col; Vector<Complex_wp> val;
          //
          //     int num = 0;
          //     col.Reallocate(nC_Hs); val.Reallocate(nC_Hs);
          //     // partie avec Ch
          //     for (int j = 0; j < nC_Hs; j++)
          //       {
          //         int j_Glob = vec_Ch(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //         Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
          //         if (j >= 0)
          //           {
          //             col(num) = jDir;
          //             // partie Iwp Ch V2t
          //             val(num) = Iwp * vec_Ch(ref).Value(i, j) * ktr2;
          //             num++;
          //           }
          //       }
          //
          //     A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          //
          //     num = 0;
          //     col.Reallocate(2*nU_Hs); val.Reallocate(2*nU_Hs);
          //     // parties avec Uh
          //     for (int j = 0; j < nU_Hs; j++)
          //       {
          //         int j_Glob = vec_Uh(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //
          //         Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
          //         col(num) = nb_dof_Es + j_Glob;
          //         // partie 2 i mu_tilde tau Uh W2t
          //         val(num) = 2.0 * Iwp * tau * mu_tilde * vec_Uh(ref).Value(i, j) * ktr2;
          //         num++;
          //
          //         if (jDir >= 0)
          //           {
          //             col(num) = jDir;
          //             // partie 2 i k2 tau^2 Uh V2
          //             val(num) = 2.0 * Iwp * k2 * tau * tau * vec_Uh(ref).Value(i, j);
          //             num++;
          //           }
          //       }
          //
          //     A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          //   }
        }
      else if (GetMatrixIndex(num) == 2)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = - eps_sigma * mu_tilde;
          Real_wp tau = var.ref_tau(ref);
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);

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
                  val(num) = Iwp * vec_Ch(ref).Value(i, j);
                  num++;
                }

              A.AddInteractionRow(iDir, num, col, val);

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
                      val(num) = -2.0 * Iwp * tau * eps_sigma * vec_Uh(ref).Value(i, j);
                      num++;
                    }
                }

              A.AddInteractionRow(iDir, num, col, val);
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
                      val(num) = - Iwp * vec_Ch(ref).Value(i, j);
                      num++;
                    }
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);

              num = 0;
              col.Reallocate(nU_Hs); val.Reallocate(nU_Hs);
              // parties avec Uh
              for (int j = 0; j < nU_Hs; j++)
                {
                  int j_Glob = vec_Uh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);

                  col(num) = nb_dof_Es + j_Glob;
                  // partie - 2 i mu_tilde tau Uh W2
                  val(num) = - 2.0 * Iwp * tau * mu_tilde * vec_Uh(ref).Value(i, j);
                  num++;
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);
            }

          // int ref = this->split_reference(num/4);
          // Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          // Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          // Complex_wp k2 = - eps_sigma * mu_tilde;
          // Real_wp tau = var.torsion;
          // A.Clear();
          // A.Reallocate(nb_rows, nb_rows);
          //
          //
          // // on traite les lignes avec V2
          // for (int iDir = 0; iDir < nb_dof_Es; iDir++)
          //   {
          //     int i = DofKeptDir(iDir); // i = ligne globale
          //     int nM_Hs = vec_Mh(ref).GetRowSize(i);
          //     int nK_Hs = vec_Kh(ref).GetRowSize(i);
          //     int nD_Hs = vec_Dh(ref).GetRowSize(i);
          //
          //     Vector<int> col; Vector<Complex_wp> val;
          //     int num = 0;
          //     col.Reallocate(2*nM_Hs); val.Reallocate(2*nM_Hs);
          //     // parties avec Uh
          //     for (int j = 0; j < nM_Hs; j++)
          //       {
          //         int j_Glob = vec_Mh(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //         Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
          //         if (jDir >= 0)
          //           {
          //             col(num) = jDir;
          //             // partie - eps_sigma Mh V2t
          //             val(num) = - eps_sigma * vec_Mh(ref).Value(i, j) * ktr2;
          //             num++;
          //           }
          //
          //         col(num) = nb_dof_Es + j_Glob;
          //         // partie 2 k2 tau Mh W2
          //         val(num) = 2.0 * k2 * tau * vec_Mh(ref).Value(i, j);
          //         num++;
          //       }
          //
          //     A.AddInteractionRow(iDir, num, col, val);
          //
          //     num = 0;
          //     col.Reallocate(nK_Hs); val.Reallocate(nK_Hs);
          //     // partie avec Kh
          //     for (int j = 0; j < nK_Hs; j++)
          //       {
          //         int j_Glob = vec_Kh(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //         if (jDir >= 0)
          //           {
          //             col(num) = jDir;
          //             // partie eps_sigma Kh V2
          //             val(num) = eps_sigma * vec_Kh(ref).Value(i, j);
          //             num++;
          //           }
          //       }
          //
          //     A.AddInteractionRow(iDir, num, col, val);
          //
          //     num = 0;
          //     col.Reallocate(nD_Hs); val.Reallocate(nD_Hs);
          //     // partie avec Dh
          //     for (int j = 0; j < nD_Hs; j++)
          //       {
          //         int j_Glob = vec_Dh(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //         if (jDir >= 0)
          //           {
          //             col(num) = jDir;
          //             // partie eps_sigma tau^2 Dh V2
          //             val(num) = eps_sigma * tau * tau * vec_Dh(ref).Value(i, j);
          //             num++;
          //           }
          //       }
          //
          //     A.AddInteractionRow(iDir, num, col, val);
          //
          //   }
          //
          // // on traite les lignes avec W2
          // for (int i = 0; i < nb_dof_Hs; i++)
          //   {
          //     int nM_Hs = vec_Mh(ref).GetRowSize(i);
          //     int nK_Hs = vec_Kh(ref).GetRowSize(i);
          //     int nD_Hs = vec_Dh(ref).GetRowSize(i);
          //
          //     Vector<int> col; Vector<Complex_wp> val;
          //
          //     int num = 0;
          //     col.Reallocate(2*nM_Hs); val.Reallocate(2*nM_Hs);
          //     // parties avec Uh
          //     for (int j = 0; j < nM_Hs; j++)
          //       {
          //         int j_Glob = vec_Mh(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //         Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
          //
          //         col(num) = nb_dof_Es + j_Glob;
          //         // partie - mu_tilde Mh W2t
          //         val(num) = - mu_tilde * vec_Mh(ref).Value(i, j) * ktr2;
          //         num++;
          //
          //         if (jDir >= 0)
          //           {
          //             col(num) = jDir;
          //             // partie - 2 k2 tau Mh, V2
          //             val(num) = - 2.0 * k2 * tau * vec_Mh(ref).Value(i, j);
          //             num++;
          //           }
          //       }
          //
          //     A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          //
          //     num = 0;
          //     col.Reallocate(nK_Hs); val.Reallocate(nK_Hs);
          //     // partie avec Kh
          //     for (int j = 0; j < nK_Hs; j++)
          //       {
          //         int j_Glob = vec_Kh(ref).Index(i, j);
          //         // int jDir = IndexDirichlet(j_Glob);
          //         // Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob, j_Glob);
          //         col(num) = nb_dof_Es + j_Glob;
          //         // partie mu_tilde Kh W2
          //         val(num) = mu_tilde * vec_Kh(ref).Value(i, j);
          //         num++;
          //       }
          //
          //     A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          //
          //     num = 0;
          //     col.Reallocate(nD_Hs); val.Reallocate(nD_Hs);
          //     // partie avec Dh
          //     for (int j = 0; j < nD_Hs; j++)
          //       {
          //         int j_Glob = vec_Dh(ref).Index(i, j);
          //         // int jDir = IndexDirichlet(j_Glob);
          //         // Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob, j_Glob);
          //         col(num) = nb_dof_Es + j_Glob;
          //         // partie mu_tilde tau^2 Dh W2
          //         val(num) = mu_tilde * tau * tau * vec_Dh(ref).Value(i, j);
          //         num++;
          //       }
          //
          //     A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          //   }
        }
      else if (GetMatrixIndex(num) == 3)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = - eps_sigma * mu_tilde;
          Real_wp tau = var.ref_tau(ref);
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);

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
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                  col(num) = nb_dof_Es + j_Glob;
                  // partie - Iwp Ch W2t
                  val(num) = - Iwp * vec_Ch(ref).Value(i, j) * ktr2;
                  num++;
                }

              A.AddInteractionRow(iDir, num, col, val);

              num = 0;
              col.Reallocate(2*nU_Hs); val.Reallocate(2*nU_Hs);
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
                      val(num) = 2.0 * Iwp * tau * eps_sigma * vec_Uh(ref).Value(i, j) * ktr2;
                      num++;
                    }

                  col(num) = nb_dof_Es + j_Glob;
                  // partie -2 Iwp k2 tau^2 Uh W2
                  val(num) = - 2.0 * Iwp * k2 * tau * tau * vec_Uh(ref).Value(i, j);
                  num++;
                }

              A.AddInteractionRow(iDir, num, col, val);
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
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                  if (j >= 0)
                    {
                      col(num) = jDir;
                      // partie Iwp Ch V2t
                      val(num) = Iwp * vec_Ch(ref).Value(i, j) * ktr2;
                      num++;
                    }
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);

              num = 0;
              col.Reallocate(2*nU_Hs); val.Reallocate(2*nU_Hs);
              // parties avec Uh
              for (int j = 0; j < nU_Hs; j++)
                {
                  int j_Glob = vec_Uh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);

                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                  col(num) = nb_dof_Es + j_Glob;
                  // partie 2 i mu_tilde tau Uh W2t
                  val(num) = 2.0 * Iwp * tau * mu_tilde * vec_Uh(ref).Value(i, j) * ktr2;
                  num++;

                  if (jDir >= 0)
                    {
                      col(num) = jDir;
                      // partie 2 i k2 tau^2 Uh V2
                      val(num) = 2.0 * Iwp * k2 * tau * tau * vec_Uh(ref).Value(i, j);
                      num++;
                    }
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);
            }






          // int ref = this->split_reference(num/4);
          // Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          // Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          // Complex_wp k2 = - eps_sigma * mu_tilde;
          // Real_wp tau = var.torsion;
          // A.Clear();
          // A.Reallocate(nb_rows, nb_rows);
          // cout << "Allocation" << endl;
          //
          // // on traite les lignes avec V2
          // for (int iDir = 0; iDir < nb_dof_Es; iDir++)
          //   {
          //     int i = DofKeptDir(iDir); // i = ligne globale
          //     int nM_Hs = vec_Mh(ref).GetRowSize(i);
          //     Vector<int> col; Vector<Complex_wp> val;
          //     int num = 0;
          //     col.Reallocate(nM_Hs); val.Reallocate(nM_Hs);
          //     // parties avec Uh
          //     for (int j = 0; j < nM_Hs; j++)
          //       {
          //         int j_Glob = vec_Mh(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //         if (jDir >= 0)
          //           {
          //             col(num) = jDir;
          //             // partie eps_sigma M_h V2
          //             val(num) = eps_sigma * vec_Mh(ref).Value(i, j);
          //             num++;
          //           }
          //       }
          //     A.AddInteractionRow(iDir, num, col, val);
          //   }
          //
          // // on traite les lignes avec W2
          // for (int i = 0; i < nb_dof_Hs; i++)
          //   {
          //     int nM_Hs = vec_Mh(ref).GetRowSize(i);
          //
          //     Vector<int> col; Vector<Complex_wp> val;
          //
          //     int num = 0;
          //     col.Reallocate(nM_Hs); val.Reallocate(nM_Hs);
          //     // parties avec Uh
          //     for (int j = 0; j < nM_Hs; j++)
          //       {
          //         int j_Glob = vec_Mh(ref).Index(i, j);
          //         Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
          //
          //         col(num) = nb_dof_Es + j_Glob;
          //         // partie mu_tilde Mh W2
          //         val(num) = mu_tilde * vec_Mh(ref).Value(i, j);
          //         num++;
          //       }
          //
          //     A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          //   }
        }
      else
        {
          DISP(num);
          DISP(ref);
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = - eps_sigma * mu_tilde;
          Real_wp tau = var.ref_tau(ref);
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);

          // on traite les lignes avec V2
          for (int iDir = 0; iDir < nb_dof_Es; iDir++)
            {
              int i = DofKeptDir(iDir); // i = ligne globale
              int nM_Hs = vec_Mh(ref).GetRowSize(i);
              int nT_Hs = vec_Th(ref).GetRowSize(i);
              int nK_Hs = vec_Kh(ref).GetRowSize(i);
              int nD_Hs = vec_Dh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;
              int num = 0;
              col.Reallocate(nM_Hs); val.Reallocate(nM_Hs);
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
                      val(num) = 4 * k2 * tau * tau * eps_sigma * vec_Mh(ref).Value(i, j);
                      num++;
                    }
                }

              A.AddInteractionRow(iDir, num, col, val);

              num = 0;
              col.Reallocate(nT_Hs); val.Reallocate(nT_Hs);
              // parties avec Uh
              for (int j = 0; j < nT_Hs; j++)
                {
                  int j_Glob = vec_Th(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);
                  if (jDir >= 0)
                    {
                      col(num) = jDir;
                      // partie - eps_sigma Mh V2t
                      val(num) = 2 * k2 * tau * tau * eps_sigma * vec_Th(ref).Value(i, j);
                      num++;
                    }
                }

              A.AddInteractionRow(iDir, num, col, val);

              num = 0;
              col.Reallocate(nK_Hs); val.Reallocate(nK_Hs);
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
                      val(num) = - eps_sigma * vec_Kh(ref).Value(i, j) * ktr2;
                      num++;
                    }
                }

              A.AddInteractionRow(iDir, num, col, val);

              num = 0;
              col.Reallocate(nD_Hs); val.Reallocate(nD_Hs);
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
                      val(num) = eps_sigma * tau * tau * vec_Dh(ref).Value(i, j) * ktr2;
                      num++;
                    }
                }

              A.AddInteractionRow(iDir, num, col, val);

            }

          // on traite les lignes avec W2
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              int nM_Hs = vec_Mh(ref).GetRowSize(i);
              int nT_Hs = vec_Th(ref).GetRowSize(i);
              int nK_Hs = vec_Kh(ref).GetRowSize(i);
              int nD_Hs = vec_Dh(ref).GetRowSize(i);

              Vector<int> col; Vector<Complex_wp> val;

              int num = 0;
              col.Reallocate(nM_Hs); val.Reallocate(nM_Hs);
              // parties avec Uh
              for (int j = 0; j < nM_Hs; j++)
                {
                  int j_Glob = vec_Mh(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);

                  col(num) = nb_dof_Es + j_Glob;
                  // partie - mu_tilde Mh W2t
                  val(num) = 4 * k2 * tau * tau * mu_tilde * vec_Mh(ref).Value(i, j);
                  num++;
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);

              num = 0;
              col.Reallocate(nT_Hs); val.Reallocate(nT_Hs);
              // parties avec Uh
              for (int j = 0; j < nT_Hs; j++)
                {
                  int j_Glob = vec_Th(ref).Index(i, j);
                  int jDir = IndexDirichlet(j_Glob);

                  col(num) = nb_dof_Es + j_Glob;
                  // partie - mu_tilde Mh W2t
                  val(num) = 2 * k2 * tau * tau * mu_tilde * vec_Th(ref).Value(i, j);
                  num++;
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);

              num = 0;
              col.Reallocate(nK_Hs); val.Reallocate(nK_Hs);
              // partie avec Kh
              for (int j = 0; j < nK_Hs; j++)
                {
                  int j_Glob = vec_Kh(ref).Index(i, j);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob,j_Glob);
                  // int jDir = IndexDirichlet(j_Glob);
                  // Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob, j_Glob);
                  col(num) = nb_dof_Es + j_Glob;
                  // partie mu_tilde Kh W2
                  val(num) = - mu_tilde * vec_Kh(ref).Value(i, j) * ktr2;
                  num++;
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);

              num = 0;
              col.Reallocate(nD_Hs); val.Reallocate(nD_Hs);
              // partie avec Dh
              for (int j = 0; j < nD_Hs; j++)
                {
                  int j_Glob = vec_Dh(ref).Index(i, j);
                  // int jDir = IndexDirichlet(j_Glob);
                  Complex_wp ktr2 = vec_ktrSquared(ref).Get(j_Glob, j_Glob);
                  col(num) = nb_dof_Es + j_Glob;
                  // partie mu_tilde tau^2 Dh W2
                  val(num) = - mu_tilde * tau * tau * vec_Dh(ref).Value(i, j) * ktr2;
                  num++;
                }

              A.AddInteractionRow(nb_dof_Es + i, num, col, val);
            }





          // int ref = this->split_reference(num/4-1);
          // Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          // Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          // Complex_wp k2 = - eps_sigma * mu_tilde;
          // Real_wp tau = var.torsion;
          // A.Clear();
          // A.Reallocate(nb_rows, nb_rows);
          //
          // // on traite les lignes avec V2
          // for (int iDir = 0; iDir < nb_dof_Es; iDir++)
          //   {
          //     int i = DofKeptDir(iDir); // i = ligne globale
          //     int nC_Hs = vec_Ch(ref).GetRowSize(i);
          //     int nU_Hs = vec_Uh(ref).GetRowSize(i);
          //
          //     Vector<int> col; Vector<Complex_wp> val;
          //     int num = 0;
          //     col.Reallocate(nC_Hs); val.Reallocate(nC_Hs);
          //     // partie avec Ch
          //     for (int j = 0; j < nC_Hs; j++)
          //       {
          //         int j_Glob = vec_Ch(ref).Index(i, j);
          //         col(num) = nb_dof_Es + j_Glob;
          //         // partie Iwp Ch W2
          //         val(num) = Iwp * vec_Ch(ref).Value(i, j);
          //         num++;
          //       }
          //
          //     A.AddInteractionRow(iDir, num, col, val);
          //
          //     num = 0;
          //     col.Reallocate(nU_Hs); val.Reallocate(nU_Hs);
          //     // parties avec Uh
          //     for (int j = 0; j < nU_Hs; j++)
          //       {
          //         int j_Glob = vec_Uh(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //         if (jDir >= 0)
          //           {
          //             col(num) = jDir;
          //             // partie -2 i tau eps_sigma Uh V2
          //             val(num) = -2.0 * Iwp * tau * eps_sigma * vec_Uh(ref).Value(i, j);
          //             num++;
          //           }
          //       }
          //
          //     A.AddInteractionRow(iDir, num, col, val);
          //   }
          //
          // // on traite les lignes avec W2
          // for (int i = 0; i < nb_dof_Hs; i++)
          //   {
          //     int nC_Hs = vec_Ch(ref).GetRowSize(i);
          //     int nU_Hs = vec_Uh(ref).GetRowSize(i);
          //
          //     Vector<int> col; Vector<Complex_wp> val;
          //
          //     int num = 0;
          //     col.Reallocate(nC_Hs); val.Reallocate(nC_Hs);
          //     // partie avec Ch
          //     for (int j = 0; j < nC_Hs; j++)
          //       {
          //         int j_Glob = vec_Ch(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //         if (j >= 0)
          //           {
          //             col(num) = jDir;
          //             // partie - i Ch V2
          //             val(num) = - Iwp * vec_Ch(ref).Value(i, j);
          //             num++;
          //           }
          //       }
          //
          //     A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          //
          //     num = 0;
          //     col.Reallocate(nU_Hs); val.Reallocate(nU_Hs);
          //     // parties avec Uh
          //     for (int j = 0; j < nU_Hs; j++)
          //       {
          //         int j_Glob = vec_Uh(ref).Index(i, j);
          //         int jDir = IndexDirichlet(j_Glob);
          //
          //         col(num) = nb_dof_Es + j_Glob;
          //         // partie - 2 i mu_tilde tau Uh W2
          //         val(num) = - 2.0 * Iwp * tau * mu_tilde * vec_Uh(ref).Value(i, j);
          //         num++;
          //       }
          //
          //     A.AddInteractionRow(nb_dof_Es + i, num, col, val);
          //   }
        }
    }
  else
    {
      DISP(num);
      if (num == 0)
        {
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              // partie (-i omega epsilon + sigma) M_h Ez
              int i0 = IndexDirichlet(i);
              int nb_val = 0;
              for (int j = 0; j < A1_sum.GetRowSize(i); j++)
                if (IndexDirichlet(A1_sum.Index(i, j)) >= 0)
                  nb_val++;

              if (i0 >= 0)
                {
                  A.ReallocateRow(i0, nb_val);
                  nb_val = 0;
                  for (int j = 0; j < A1_sum.GetRowSize(i); j++)
                    if (IndexDirichlet(A1_sum.Index(i, j)) >= 0)
                      {
                        A.Index(i0, nb_val) = IndexDirichlet(A1_sum.Index(i, j));
                        A.Value(i0, nb_val) = A1_sum.Value(i, j);
                        nb_val++;
                      }
                }

              // partie -i omega mu M_h H_z
              nb_val = A1_sumH.GetRowSize(i);
              A.ReallocateRow(nb_dof_Es + i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(nb_dof_Es + i, j) = nb_dof_Es + A1_sumH.Index(i, j);
                  A.Value(nb_dof_Es + i, j) = A1_sumH.Value(i, j);
                }
            }
        }
      else if (GetMatrixIndex(num) == 1)
        {
          DISP(num);
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp coef_mu = -Iwp*var.omega*var.ref_mu(ref);
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              // partie (-i omega epsilon + sigma)/Delta K_h Ez
              int i0 = IndexDirichlet(i);
              int nb_val = 0;
              for (int j = 0; j < vec_Kh(ref).GetRowSize(i); j++)
                if (IndexDirichlet(vec_Kh(ref).Index(i, j)) >= 0)
                  nb_val++;

              if (i0 >= 0)
                {
                  A.ReallocateRow(i0, nb_val);
                  nb_val = 0;
                  for (int j = 0; j < vec_Kh(ref).GetRowSize(i); j++)
                    if (IndexDirichlet(vec_Kh(ref).Index(i, j)) >= 0)
                      {
                        A.Index(i0, nb_val) = IndexDirichlet(vec_Kh(ref).Index(i, j));
                        A.Value(i0, nb_val) = eps_sigma*vec_Kh(ref).Value(i, j);
                        nb_val++;
                      }
                }

              // partie -i omega mu K_h H_z
              nb_val = vec_Kh(ref).GetRowSize(i);
              A.ReallocateRow(nb_dof_Es + i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(nb_dof_Es + i, j) = nb_dof_Es + vec_Kh(ref).Index(i, j);
                  A.Value(nb_dof_Es + i, j) = coef_mu*vec_Kh(ref).Value(i, j);
                }
            }
        }
      else
        {
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          for (int i = 0; i < nb_dof_Hs; i++)
            {
              // partie i beta/Delta C_h Hz
              int i0 = IndexDirichlet(i);
              if (i0 >= 0)
                {
                  int nb_val = vec_Ch(ref).GetRowSize(i);
                  A.ReallocateRow(i0, nb_val);
                  for (int j = 0; j < nb_val; j++)
                    {
                      A.Index(i0, j) = nb_dof_Es + vec_Ch(ref).Index(i, j);
                      A.Value(i0, j) = Iwp * vec_Ch(ref).Value(i, j);
                    }
                }

              // partie -i beta/Delta C_h E_z
              int nb_val = 0;
              for (int j = 0; j < vec_Ch(ref).GetRowSize(i); j++)
                if (IndexDirichlet(vec_Ch(ref).Index(i, j)) >= 0)
                  nb_val++;

              A.ReallocateRow(nb_dof_Es + i, nb_val);
              nb_val = 0;
              for (int j = 0; j < vec_Ch(ref).GetRowSize(i); j++)
                if (IndexDirichlet(vec_Ch(ref).Index(i, j)) >= 0)
                  {
                    A.Index(nb_dof_Es + i, nb_val) = IndexDirichlet(vec_Ch(ref).Index(i, j));
                    A.Value(nb_dof_Es + i, nb_val) = - Iwp * vec_Ch(ref).Value(i, j);
                    nb_val++;
                  }
            }
        }
    }
}

void ModeEsNonLinearEigenProblem
::MltOperatorSplit(int num, const SeldonTranspose& trans, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  // DISP(num);
  // DISP(GetMatrixIndex(num));
  // DISP(num%4);
  // DISP(num%4==GetMatrixIndex(num));
  // if(num%4!=GetMatrixIndex(num))
  //   abort();
  int ref = this->split_reference(num-1);
  Vector<Complex_wp> X(nb_dof_all), Y(nb_dof_all);
  Complex_wp one(1, 0);
  ExtractDirichlet(Xdir, X);
  Y.Zero();

  int N = this->nb_dof_Hs;
  Vector<Complex_wp> F2(N), G2(N), Y1(N), Y2(N), F2t(N), G2t(N);
  Y1.Zero(); Y2.Zero();

  for (int i = 0 ; i < N ; i++)
    {
      F2(i) = X(i);
      G2(i) = X(N + i);
    }

  if (vec_twisted(ref))
    {
      // cout << "--- Torsadé ---" << endl;
      // DISP(num);
      // DISP(ref);
      // DISP(GetMatrixIndex(num));
      // DISP(var.ref_tau(ref));
      if (num == 0)
        {
          Mlt(A1_sum, F2, Y1);
          Mlt(A1_sumH, G2, Y2);
        }
      else if (GetMatrixIndex(num) == 1) // terme après f_2
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = -eps_sigma*mu_tilde;
          Matrix<Complexe, Prop, Storage> temp;

          Mlt(vec_ktrSquared(ref), F2, F2t);
          Mlt(vec_ktrSquared(ref), G2, G2t);

          MltAdd(- eps_sigma, vec_Mh(ref), F2t, one, Y1);
          MltAdd(- mu_tilde, vec_Mh(ref), G2t, one, Y2);
          MltAdd(2*k2*var.ref_tau(ref), vec_Mh(ref), G2, one, Y1);
          MltAdd(- 2*k2*var.ref_tau(ref), vec_Mh(ref), F2, one, Y2);
          MltAdd(eps_sigma, vec_Kh(ref), F2, one, Y1);
          MltAdd(mu_tilde, vec_Kh(ref), G2, one, Y2);
          MltAdd(eps_sigma*var.ref_tau(ref)*var.ref_tau(ref), vec_Dh(ref), F2, one, Y1);
          MltAdd(mu_tilde*var.ref_tau(ref)*var.ref_tau(ref), vec_Dh(ref), G2, one, Y2);
        }
      else if (GetMatrixIndex(num) == 2) // terme après f_3
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          MltAdd(- 2 * Iwp * eps_sigma * var.ref_tau(ref), vec_Uh(ref), F2, one, Y1);
          MltAdd(- 2 * Iwp * mu_tilde * var.ref_tau(ref), vec_Uh(ref), G2, one, Y2);
          MltAdd(Iwp, vec_Ch(ref), G2, one, Y1);
          MltAdd(- Iwp, vec_Ch(ref), F2, one, Y2);
        }
      else if (GetMatrixIndex(num) == 3) // terme après f_4
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = -eps_sigma*mu_tilde;
          Matrix<Complexe, Prop, Storage> temp;

          Mlt(vec_ktrSquared(ref), F2, F2t);
          Mlt(vec_ktrSquared(ref), G2, G2t);

          MltAdd(2 * Iwp * eps_sigma * var.ref_tau(ref), vec_Uh(ref), F2t, one, Y1);
          MltAdd(2 * Iwp * mu_tilde * var.ref_tau(ref), vec_Uh(ref), G2t, one, Y2);
          MltAdd(-2 * Iwp * k2 * var.ref_tau(ref) * var.ref_tau(ref), vec_Uh(ref), G2, one, Y1);
          MltAdd(2 * Iwp * k2 * var.ref_tau(ref) * var.ref_tau(ref), vec_Uh(ref), F2, one, Y2);
          MltAdd(- Iwp, vec_Ch(ref), G2t, one, Y1);
          MltAdd(Iwp, vec_Ch(ref), F2t, one, Y2);
        }
      else if (GetMatrixIndex(num) == 0) // terme après f_5
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp k2 = -eps_sigma*mu_tilde;
          Real_wp tau = var.ref_tau(ref);


          Mlt(vec_ktrSquared(ref), F2, F2t);
          Mlt(vec_ktrSquared(ref), G2, G2t);

          MltAdd(4 * k2 * tau * tau * eps_sigma, vec_Mh(ref), F2, one, Y1);
          MltAdd(4 * k2 * tau * tau * mu_tilde, vec_Mh(ref), G2, one, Y2);
          MltAdd(2 * k2 * tau * tau * eps_sigma, vec_Th(ref), F2, one, Y1);
          MltAdd(2 * k2 * tau * tau * mu_tilde, vec_Th(ref), G2, one, Y2);
          MltAdd(-eps_sigma, vec_Kh(ref), F2t, one, Y1);
          MltAdd(-mu_tilde, vec_Kh(ref), G2t, one, Y2);
          MltAdd(-eps_sigma * tau * tau, vec_Dh(ref), F2t, one, Y1);
          MltAdd(-mu_tilde * tau * tau, vec_Dh(ref), G2t, one, Y2);
        }
    }
  else
    {
      // cout << "--- Non torsadé ---" << endl;
      // DISP(num);
      // DISP(ref);
      // DISP(GetMatrixIndex(num));
      if (num == 0)
        {
          Mlt(A1_sum, F2, Y1);
          Mlt(A1_sumH, G2, Y2);
        }
      else if (GetMatrixIndex(num) == 1)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);          Complex_wp coef_mu = -Iwp*var.omega*var.ref_mu(ref);

          Mlt(vec_Kh(ref), F2, Y1); Mlt(eps_sigma, Y1);
          Mlt(vec_Kh(ref), G2, Y2); Mlt(coef_mu, Y2);
        }
      else
        {
          Mlt(vec_Ch(ref), G2, Y1); Mlt(Iwp, Y1);
          Mlt(vec_Ch(ref), F2, Y2); Mlt(-Iwp, Y2);
        }
    }


  for (int i = 0 ; i < N ; i++)
    {
      Y(i) = Y1(i);
      Y(N + i) = Y2(i);
    }

  ExpandDirichlet(Y, Ydir);
}
