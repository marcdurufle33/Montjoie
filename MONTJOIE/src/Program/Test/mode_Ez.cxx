// initialisation du probleme
ModeEzNonLinearEigenProblem
::ModeEzNonLinearEigenProblem(int N, ModeEz_Solver& var0,
                              Vector<Matrix<Complexe, Prop, Storage> >& Mh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ch,
                              Vector<int>& DofDir, Vector<int>& IndexDir, bool use_split)
  :
  var(var0), vec_Mh(Mh), vec_Kh(Kh), vec_Ch(Ch), DofKeptDir(DofDir), IndexDirichlet(IndexDir)
{
  // pour Ez, condition de Dirichlet
  nb_dof_Ez = DofKeptDir.GetM();

  // pour Hz, condition de Neumann
  nb_dof_Hz = N;

  // taille du pb aux valeurs propres
  nb_rows = nb_dof_Ez;
  nb_dof_all = N;
  if (!only_Ez)
    {
      nb_rows += nb_dof_Hz;
      nb_dof_all += N;
    }

  this->Init(nb_rows);
  DISP(nb_rows); DISP(nb_dof_Ez); DISP(nb_dof_Hz);

  // on calcule les poles du pb aux valeurs propres
  VectComplex_wp s; split_reference.Clear();
  this->numer_pol_split.Clear();
  this->denom_pol_split.Clear();
  if (use_split)
    {
      Vector<Complexe> num, denom;
      Mh_sum.Reallocate(nb_dof_Hz, nb_dof_Hz);
      if (!only_Ez)
        Mh_sumH.Reallocate(nb_dof_Hz, nb_dof_Hz);

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
    if (vec_Kh(ref).GetM() > 0)
      {
        Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
        Complex_wp Delta = Iwp*var.omega*eps_sigma*var.ref_mu(ref);
        Complex_wp z = sqrt(Delta) / var.omega;
        Complex_wp coef_mu = -Iwp*var.omega*var.ref_mu(ref);
        s.PushBack(z); s.PushBack(-z);

        if (use_split)
          {
            Add(eps_sigma, vec_Mh(ref), Mh_sum);
            if (!only_Ez)
              Add(coef_mu, vec_Mh(ref), Mh_sumH);

            split_reference.PushBack(ref);

            Vector<Complex_wp> num(1), denom(3);
            denom(0) = var.omega*var.omega;
            denom(1) = 0.0;
            denom(2) = -Delta;

            if (only_Ez)
              {
                num(0) = eps_sigma;

                this->numer_pol_split.PushBack(num);
                this->denom_pol_split.PushBack(denom);

                string name = "Kh" + to_str(ref) + ".mtx";
                file_out << " 0 " << real(num(0)) << " " << imag(num(0)) << " " << var.omega << " "
                         << real(denom(2)) << " " << imag(denom(2)) << " " << name << endl;

                RemoveDirichletRowCol(vec_Kh(ref), A);
                WriteMatrixMarket(A, name);
              }
            else
              {
                num(0) = 1.0;

                this->numer_pol_split.PushBack(num);
                this->denom_pol_split.PushBack(denom);

                string name = "Kh" + to_str(ref) + ".mtx";
                file_out << " 0 " << real(num(0)) << " " << imag(num(0)) << " " << var.omega << " "
                         << real(denom(2)) << " " << imag(denom(2)) << " " << name << endl;
                ComputeOperatorSplitExplicit(this->numer_pol_split.GetM()-1, A);
                WriteMatrixMarket(A, name);

                num.Reallocate(2);
                num(1) = 0.0; num(0) = Iwp*var.omega;
                this->numer_pol_split.PushBack(num);
                this->denom_pol_split.PushBack(denom);

                name = "Ch" + to_str(ref) + ".mtx";
                file_out << " 1 " << real(num(0)) << " " << imag(num(0)) << " " << var.omega << " "
                         << real(denom(2)) << " " << imag(denom(2)) << " " << name << endl;
                ComputeOperatorSplitExplicit(this->numer_pol_split.GetM()-1, A);
                WriteMatrixMarket(A, name);
              }
          }
      }

  if (use_split)
    {
      file_out << "0 1.0 0.0 0.0 1.0 0.0 Mh0.mtx" << endl;
      ComputeOperatorSplitExplicit(0, A);
      WriteMatrixMarket(A, "Mh0.mtx");
      file_out.close();
    }

  this->nb_split_matrix = this->numer_pol_split.GetM();
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

void ModeEzNonLinearEigenProblem::RemoveDirichletRowCol(const Matrix<Complexe, Prop, Storage>& A,
                                                        Matrix<Complexe, Prop, Storage>& B)
{
  B.Clear();
  B.Reallocate(nb_dof_Ez, nb_dof_Ez);
  for (int i = 0; i < nb_dof_Ez; i++)
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

void ModeEzNonLinearEigenProblem::ExtractDirichlet(const VectComplex_wp& Xdir, VectComplex_wp& X)
{
  X.Zero();
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    X(DofKeptDir(i)) = Xdir(i);

  if (!only_Ez)
    for (int i = 0; i < nb_dof_Hz; i++)
      X(nb_dof_Hz + i) = Xdir(nb_dof_Ez + i);
}

void ModeEzNonLinearEigenProblem::ExpandDirichlet(const VectComplex_wp& X, VectComplex_wp& Xdir)
{
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    Xdir(i) = X(DofKeptDir(i));

  if (!only_Ez)
    for (int i = 0; i < nb_dof_Hz; i++)
      Xdir(nb_dof_Ez + i) = X(nb_dof_Hz + i);
}

void ModeEzNonLinearEigenProblem::ComputeOperator(const Complex_wp& L)
{
}


// calcul de Y = T(L) X
void ModeEzNonLinearEigenProblem::MltOperator(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
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

  if (only_Ez)
    {
      for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
        if (vec_Kh(ref).GetM() > 0)
          {
            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;

            MltAdd(eps_sigma, vec_Mh(ref), X, one, Y);
            MltAdd(eps_sigma/Delta, vec_Kh(ref), X, one, Y);
          }
    }
  else
    {
      // le vecteur calculé est décomposé en deux : Ez, puis Hz
      int N = nb_dof_Hz;
      Vector<Complex_wp> Ez(N), Hz(N), Y1(N), Y2(N);
      Y1.Zero(); Y2.Zero();
      //cout << "Je passe à : MltOperator" << endl;
      for (int i = 0 ; i < N; i++)
        {
          Ez(i) = X(i);
          Hz(i) = X(N + i);
        }

      for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
        if (vec_Kh(ref).GetM() > 0)
          {
            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;

            MltAdd(eps_sigma, vec_Mh(ref), Ez, one, Y1);
            MltAdd(- Iwp * var.omega * var.ref_mu(ref), vec_Mh(ref), Hz, one, Y2);
            MltAdd(eps_sigma/Delta, vec_Kh(ref), Ez, one, Y1);
            MltAdd((- Iwp * var.omega * var.ref_mu(ref))/Delta, vec_Kh(ref), Hz, one, Y2);
            MltAdd(Iwp * beta / Delta, vec_Ch(ref), Hz, one, Y1);
            MltAdd(-Iwp * beta / Delta, vec_Ch(ref), Ez, one, Y2);
          }

      for (int i = 0 ; i < N ; i++)
        {
          Y(i) = Y1(i);
          Y(N + i) = Y2(i);
        }
    }

  ExpandDirichlet(Y, Ydir);
}

// calcul explicite de T(L) (utilise si on calcule explicitement les matrices)
void ModeEzNonLinearEigenProblem
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
  if (only_Ez)
    {

      for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
        if (vec_Mh(ref).GetM() > 0)
          {
            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
            Complex_wp coef = eps_sigma / Delta;

            // on effectue eps_sigma Mh + coef Kh en eliminant les noeuds Dirichlet
            for (int i0 = 0; i0 < nb_dof_Ez; i0++)
              {
                int i = DofKeptDir(i0);
                int n1 = vec_Mh(ref).GetRowSize(i);
                int n2 = vec_Kh(ref).GetRowSize(i);
                Vector<int> col(n1+n2); Vector<Complex_wp> val(n1+n2); int num = 0;
                for (int j0 = 0; j0 < n1; j0++)
                  {
                    int j1 = vec_Mh(ref).Index(i, j0);
                    int j = IndexDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = vec_Mh(ref).Value(i, j0)*eps_sigma;
                        num++;
                      }
                  }

                for (int j0 = 0; j0 < n2; j0++)
                  {
                    int j1 = vec_Kh(ref).Index(i, j0);
                    int j = IndexDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = vec_Kh(ref).Value(i, j0)*coef;
                        num++;
                      }
                  }

                A.AddInteractionRow(i0, num, col, val);
              }
          }
    }
  else
    {
      for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
        if (vec_Mh(ref).GetM() > 0)
          {
            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
            Complex_wp coef_11 = eps_sigma/Delta;
            Complex_wp coef_12 = Iwp * beta / Delta;

            // on traite les lignes avec Ez
            for (int i0 = 0; i0 < nb_dof_Ez; i0++)
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
                        val(num) = vec_Mh(ref).Value(i, j0)*eps_sigma;
                        num++;
                      }
                  }

                A.AddInteractionRow(i0, num, col, val);
                num = 0;

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
                    col(num) = nb_dof_Ez + j;
                    val(num) = vec_Ch(ref).Value(i, j0)*coef_12;
                    num++;
                  }

                A.AddInteractionRow(i0, num, col, val);
              }

            Complex_wp coef_22 = - Iwp * var.omega * var.ref_mu(ref);
            Complex_wp coefKh_22 = - Iwp * var.omega * var.ref_mu(ref)/Delta;
            Complex_wp coef_21 = -Iwp * beta / Delta;

            // on traite les lignes avec Hz
            for (int i = 0; i < nb_dof_Hz; i++)
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
                    col(num) = nb_dof_Ez + j;
                    val(num) = vec_Mh(ref).Value(i, j0)*coef_22;
                    num++;
                  }

                A.AddInteractionRow(nb_dof_Ez+i, num, col, val);
                num = 0;

                // partie -i omega mu/Delta Kh
                for (int j0 = 0; j0 < n2; j0++)
                  {
                    int j = vec_Kh(ref).Index(i, j0);
                    col(num) = nb_dof_Ez + j;
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

                A.AddInteractionRow(nb_dof_Ez + i, num, col, val);
              }
          }
    }
  //A.WriteText("B.dat");
  //exit(0);
}


void ModeEzNonLinearEigenProblem::ComputeJacobian(const Complex_wp& L)
{}

// calcul de Y = T'(L) X
void ModeEzNonLinearEigenProblem::MltJacobian(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  Vector<Complex_wp> X(nb_dof_all), Y(nb_dof_all);
  ExtractDirichlet(Xdir, X);
  Y.Zero();

  Complex_wp one(1, 0);
  Complex_wp beta = L*var.omega;
  if (only_Ez)
    {
      for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
        if (vec_Mh(ref).GetM() > 0)
          {
            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
            Complex_wp dDelta = 2.0*var.omega*beta;

            MltAdd(-eps_sigma/(Delta*Delta)*dDelta, vec_Kh(ref), X, one, Y);
          }
    }
  else
    {
      // le vecteur calculé est décomposé en deux : Ez, puis Hz
      int N = this->nb_dof_Hz;
      Vector<Complex_wp> Ez(N), Hz(N), Y1(N), Y2(N);
      Y1.Zero(); Y2.Zero();
      //cout << "Je passe à : MltJacobian" << endl;
      for (int i = 0 ; i < N ; i++)
        {
          Ez(i) = X(i);
          Hz(i) = X(N + i);
        }

      for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
        if (vec_Mh(ref).GetM() > 0)
          {

            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
            Complex_wp dDelta = 2.0*var.omega*beta;

            MltAdd(-eps_sigma/(Delta*Delta)*dDelta, vec_Kh(ref), Ez, one, Y1);
            MltAdd(Iwp * var.omega * var.ref_mu(ref)/(Delta * Delta) * dDelta, vec_Kh(ref), Hz, one, Y2);
            MltAdd((Iwp * var.omega * Delta - Iwp * beta * dDelta) / (Delta * Delta), vec_Ch(ref), Hz, one, Y1);
            MltAdd((-Iwp * var.omega * Delta + Iwp * beta * dDelta) / (Delta * Delta), vec_Ch(ref), Ez, one, Y2);
          }

      for (int i = 0 ; i < N ; i++)
        {
          Y(i) = Y1(i);
          Y(N + i) = Y2(i);
        }
    }

  ExpandDirichlet(Y, Ydir);
}

// calcul explicite de T'(L)
void ModeEzNonLinearEigenProblem
::ComputeJacobianExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  A.Reallocate(nb_rows, nb_rows);

  Complex_wp beta = L*var.omega;
  if (only_Ez)
    {

      for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
        if (vec_Mh(ref).GetM() > 0)
          {
            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
            Complex_wp dDelta = 2.0*var.omega*beta;
            Complex_wp coef = -eps_sigma / (Delta*Delta)*dDelta;

            for (int i0 = 0; i0 < nb_dof_Ez; i0++)
              {
                int i = DofKeptDir(i0);
                int n2 = vec_Kh(ref).GetRowSize(i);
                Vector<int> col(n2); Vector<Complex_wp> val(n2); int num = 0;
                for (int j0 = 0; j0 < n2; j0++)
                  {
                    int j1 = vec_Kh(ref).Index(i, j0);
                    int j = IndexDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = vec_Kh(ref).Value(i, j0)*coef;
                        num++;
                      }
                  }

                A.AddInteractionRow(i0, num, col, val);
              }
          }
    }
  else
    {
      for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
        if (vec_Mh(ref).GetM() > 0)
          {
            Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
            Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
            Complex_wp dDelta = 2.0*var.omega*beta;
            Complex_wp coef_11 = -eps_sigma/(Delta*Delta)*dDelta;
            Complex_wp coef_12 = (Iwp * var.omega * Delta - Iwp * beta * dDelta) / (Delta * Delta);

            // on traite les lignes avec Ez
            for (int i0 = 0; i0 < nb_dof_Ez; i0++)
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
                    col(num) = nb_dof_Ez + j;
                    val(num) = vec_Ch(ref).Value(i, j0)*coef_12;
                    num++;
                  }

                A.AddInteractionRow(i0, num, col, val);
              }

            Complex_wp coef_22 = Iwp * var.omega * var.ref_mu(ref)/(Delta * Delta) * dDelta;
            Complex_wp coef_21 = (-Iwp * var.omega * Delta + Iwp * beta * dDelta) / (Delta * Delta);

            // on traite les lignes avec Hz
            for (int i = 0; i < nb_dof_Hz; i++)
              {
                int n2 = vec_Kh(ref).GetRowSize(i);
                int n3 = vec_Ch(ref).GetRowSize(i);

                Vector<int> col(n2+n3); Vector<Complex_wp> val(n2+n3);
                int num = 0;

                for (int j0 = 0; j0 < n2; j0++)
                  {
                    int j = vec_Kh(ref).Index(i, j0);
                    col(num) = nb_dof_Ez + j;
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

                A.AddInteractionRow(nb_dof_Ez + i, num, col, val);
              }
          }
    }
}

void ModeEzNonLinearEigenProblem::ComputeExplicitPreconditioning(DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  mat_lu.Factorize(A);
}

// factorisation de T(L)
void ModeEzNonLinearEigenProblem::ComputePreconditioning(const Complex_wp& L)
{
  DistributedMatrix<Complex_wp, Prop, Storage> A;
  ComputeOperatorExplicit(L, A);

  mat_lu.Factorize(A);
}


void ModeEzNonLinearEigenProblem::ComputePreconditioning(const VectComplex_wp& L, const VectComplex_wp& coef)
{
  Matrix<Complex_wp, Prop, Storage> A;

  VectComplex_wp beta(L.GetM());
  for (int k = 0; k < L.GetM(); k++)
    beta(k) = var.omega*L(k);

  A.Reallocate(nb_rows, nb_rows);

  if (only_Ez)
    {

      for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
        if (vec_Mh(ref).GetM() > 0)
          {
            Complex_wp coef1, coef2;
            SetComplexZero(coef1); SetComplexZero(coef2);
            for (int k = 0; k < coef.GetM(); k++)
              {
                Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
                Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta(k)*beta(k);
                coef1 += eps_sigma*coef(k); coef2 += coef(k) * eps_sigma / Delta;
              }

            // on effectue eps_sigma Mh + coef Kh en eliminant les noeuds Dirichlet
            for (int i0 = 0; i0 < nb_dof_Ez; i0++)
              {
                int i = DofKeptDir(i0);
                int n1 = vec_Mh(ref).GetRowSize(i);
                int n2 = vec_Kh(ref).GetRowSize(i);
                Vector<int> col(n1+n2); Vector<Complex_wp> val(n1+n2); int num = 0;
                for (int j0 = 0; j0 < n1; j0++)
                  {
                    int j1 = vec_Mh(ref).Index(i, j0);
                    int j = IndexDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = vec_Mh(ref).Value(i, j0)*coef1;
                        num++;
                      }
                  }

                for (int j0 = 0; j0 < n2; j0++)
                  {
                    int j1 = vec_Kh(ref).Index(i, j0);
                    int j = IndexDirichlet(j1);
                    if (j >= 0)
                      {
                        col(num) = j;
                        val(num) = vec_Kh(ref).Value(i, j0)*coef2;
                        num++;
                      }
                  }

                A.AddInteractionRow(i0, num, col, val);
              }
          }
    }
  else
    {
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

            // on traite les lignes avec Ez
            for (int i0 = 0; i0 < nb_dof_Ez; i0++)
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
                    col(num) = nb_dof_Ez + j;
                    val(num) = vec_Ch(ref).Value(i, j0)*coef_12;
                    num++;
                  }

                A.AddInteractionRow(i0, num, col, val);
              }

            // on traite les lignes avec Hz
            for (int i = 0; i < nb_dof_Hz; i++)
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
                    col(num) = nb_dof_Ez + j;
                    val(num) = vec_Mh(ref).Value(i, j0)*coef_22;
                    num++;
                  }

                // partie -i omega mu/Delta Kh
                for (int j0 = 0; j0 < n2; j0++)
                  {
                    int j = vec_Kh(ref).Index(i, j0);
                    col(num) = nb_dof_Ez + j;
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

                A.AddInteractionRow(nb_dof_Ez + i, num, col, val);
              }
          }
    }

  mat_lu.Factorize(A);
}

void ModeEzNonLinearEigenProblem::ComputeSplitPreconditioning(const Vector<int>& numL, const VectComplex_wp& coef)
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
void ModeEzNonLinearEigenProblem::ApplyPreconditioning(const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
{
  //cout << "appel apply" << Norm2(X) << endl;
  mat_lu.Solve(Y, X);
}


void ModeEzNonLinearEigenProblem
::ComputeOperatorSplitExplicit(int num, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  if (only_Ez)
    {
      if (num == 0)
        RemoveDirichletRowCol(Mh_sum, A);
      else
        {
          int ref = this->split_reference(num-1);
          RemoveDirichletRowCol(vec_Kh(ref), A);
        }
    }
  else
    {
      if (num == 0)
        {
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          for (int i = 0; i < nb_dof_Hz; i++)
            {
              // partie (-i omega epsilon + sigma) M_h Ez
              int i0 = IndexDirichlet(i);
              int nb_val = 0;
              for (int j = 0; j < Mh_sum.GetRowSize(i); j++)
                if (IndexDirichlet(Mh_sum.Index(i, j)) >= 0)
                  nb_val++;

              if (i0 >= 0)
                {
                  A.ReallocateRow(i0, nb_val);
                  nb_val = 0;
                  for (int j = 0; j < Mh_sum.GetRowSize(i); j++)
                    if (IndexDirichlet(Mh_sum.Index(i, j)) >= 0)
                      {
                        A.Index(i0, nb_val) = IndexDirichlet(Mh_sum.Index(i, j));
                        A.Value(i0, nb_val) = Mh_sum.Value(i, j);
                        nb_val++;
                      }
                }

              // partie -i omega mu M_h H_z
              nb_val = Mh_sumH.GetRowSize(i);
              A.ReallocateRow(nb_dof_Ez + i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(nb_dof_Ez + i, j) = nb_dof_Ez + Mh_sumH.Index(i, j);
                  A.Value(nb_dof_Ez + i, j) = Mh_sumH.Value(i, j);
                }
            }
        }
      else if (num%2 == 1)
        {
          int ref = this->split_reference(num/2);
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp coef_mu = -Iwp*var.omega*var.ref_mu(ref);
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          for (int i = 0; i < nb_dof_Hz; i++)
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
              A.ReallocateRow(nb_dof_Ez + i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(nb_dof_Ez + i, j) = nb_dof_Ez + vec_Kh(ref).Index(i, j);
                  A.Value(nb_dof_Ez + i, j) = coef_mu*vec_Kh(ref).Value(i, j);
                }
            }
        }
      else
        {
          int ref = this->split_reference(num/2-1);
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          for (int i = 0; i < nb_dof_Hz; i++)
            {
              // partie i beta/Delta C_h Hz
              int i0 = IndexDirichlet(i);
              if (i0 >= 0)
                {
                  int nb_val = vec_Ch(ref).GetRowSize(i);
                  A.ReallocateRow(i0, nb_val);
                  for (int j = 0; j < nb_val; j++)
                    {
                      A.Index(i0, j) = nb_dof_Ez + vec_Ch(ref).Index(i, j);
                      A.Value(i0, j) = vec_Ch(ref).Value(i, j);
                    }
                }

              // partie -i beta/Delta C_h E_z
              int nb_val = 0;
              for (int j = 0; j < vec_Ch(ref).GetRowSize(i); j++)
                if (IndexDirichlet(vec_Ch(ref).Index(i, j)) >= 0)
                  nb_val++;

              A.ReallocateRow(nb_dof_Ez + i, nb_val);
              nb_val = 0;
              for (int j = 0; j < vec_Ch(ref).GetRowSize(i); j++)
                if (IndexDirichlet(vec_Ch(ref).Index(i, j)) >= 0)
                  {
                    A.Index(nb_dof_Ez + i, nb_val) = IndexDirichlet(vec_Ch(ref).Index(i, j));
                    A.Value(nb_dof_Ez + i, nb_val) = -vec_Ch(ref).Value(i, j);
                    nb_val++;
                  }
            }
        }
    }
}

void ModeEzNonLinearEigenProblem
::MltOperatorSplit(int num, const SeldonTranspose& trans, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  Vector<Complex_wp> X(nb_dof_all), Y(nb_dof_all);
  ExtractDirichlet(Xdir, X);
  Y.Zero();

  if (only_Ez)
    {
      if (num == 0)
        Mlt(trans, Mh_sum, X, Y);
      else
        {
          int ref = this->split_reference(num-1);
          Mlt(trans, vec_Kh(ref), X, Y);
        }
    }
  else
    {
      int N = this->nb_dof_Hz;
      Vector<Complex_wp> Ez(N), Hz(N), Y1(N), Y2(N);
      Y1.Zero(); Y2.Zero();

      for (int i = 0 ; i < N ; i++)
        {
          Ez(i) = X(i);
          Hz(i) = X(N + i);
        }

      if (num == 0)
        {
          Mlt(Mh_sum, Ez, Y1);
          Mlt(Mh_sumH, Hz, Y2);
        }
      else if (num%2 == 1)
        {
          int ref = this->split_reference(num/2);
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);          Complex_wp coef_mu = -Iwp*var.omega*var.ref_mu(ref);

          Mlt(vec_Kh(ref), Ez, Y1); Mlt(eps_sigma, Y1);
          Mlt(vec_Kh(ref), Hz, Y2); Mlt(coef_mu, Y2);
        }
      else
        {
          int ref = this->split_reference(num/2-1);
          Mlt(vec_Ch(ref), Hz, Y1); 
          Mlt(vec_Ch(ref), Ez, Y2); Y2 = -Y2;
        }

      for (int i = 0 ; i < N ; i++)
        {
          Y(i) = Y1(i);
          Y(N + i) = Y2(i);
        }
    }

  ExpandDirichlet(Y, Ydir);
}
