void ModeEs_Solver::ComputeTwistedMaterialsDiagonalMatrix_Bis(const Vector<Matrix<Real_wp, Symmetric, DiagonalRow>>& vec_diag_tau0_Bis,
                                            Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                                            Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_ktrSquared,
                                            Vector<int>& nbRefNoeud)
{
  cout << "Construction des matrices de matériaux torsadés" << endl;
  // à refaire intégralement lors de la prise en compte des torsades multiples
  vec_D_tau0_Bis.Reallocate(this->ref_epsilon.GetM());
  vec_ktrSquared.Reallocate(this->ref_epsilon.GetM());
  int N = this->mesh_num.GetNbDof();
  nbRefNoeud.Reallocate(N);
  nbRefNoeud.Zero();

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  for (int ref = 0; ref < index_used.GetM(); ref++)
    if (index_used(ref))
      {
        DISP(ref);
        vec_D_tau0_Bis(ref).Reallocate(N,N);
        vec_ktrSquared(ref).Reallocate(N,N);
        Complex_wp k2 = Iwp * this->omega * this->ref_mu(ref) * (-Iwp * this->omega * this->ref_epsilon(ref) + this->ref_sigma(ref));
        for (int i = 0 ; i < N ; i++)
          {
            Real_wp val = vec_diag_tau0_Bis(ref)(i,i);
            if (val > 0.9) // tous les cas intéressants sont >= 1, on met 0.9 pour ne pas conter les éventuels 0 non comptabilisés numériquement
              {
                vec_D_tau0_Bis(ref).Get(i,i) = -k2 * val;
                vec_ktrSquared(ref).Get(i,i) = k2 * (val - 1.0);
                nbRefNoeud.Get(i) += 1;
              }
            else
              {
                vec_D_tau0_Bis(ref).Get(i,i) = 0.0;
                vec_ktrSquared(ref).Get(i,i) = 0.0;
              }
          }
      }
}


void ModeEs_Solver::ComputeMatrixV2ToV1(const Complex_wp& L,
                          const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                          const Vector<int>& nbRefNoeud,
                          Matrix<Complexe, Symmetric, DiagonalRow>& Dtau)
{
  int N = this->mesh_num.GetNbDof();
  Complex_wp beta = L * this->omega;
  Dtau.Reallocate(N,N);

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  for (int i = 0 ; i < N ; i++)
    {
      Dtau.Get(i,i) = 0.0;
      for (int ref = 0; ref < index_used.GetM(); ref++)
        if (index_used(ref))
          {
            Matrix<Complexe, Symmetric, DiagonalRow> D = vec_D_tau0_Bis(ref);
            if (D(i,i) != 0.0)
              Dtau.Get(i,i) += (D.Get(i,i) + beta*beta) / double(nbRefNoeud.Get(i));
              // Dtau.Get(i,i) += (D.Get(i,i) + beta*beta);
            // DISP(i+1);
            // cout << D(i,i) << ", beta : " << beta << endl;
            // DISP(Dtau2.Get(i,i));
            // DISP(nbRefNoeud.Get(i));
          }
    }
}

void ModeEs_Solver::ComputeMatrixV2ToEs(const Complex_wp& L,
                          const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                          const Vector<int>& nbRefNoeud,
                          Matrix<Complexe, Symmetric, DiagonalRow>& Dtau2)
{
  int N = this->mesh_num.GetNbDof();
  Complex_wp beta = L * this->omega;
  Dtau2.Reallocate(N,N);

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  for (int i = 0 ; i < N ; i++)
    {
      Dtau2.Get(i,i) = 0.0;
      for (int ref = 0; ref < index_used.GetM(); ref++)
        if (index_used(ref))
          {
            Matrix<Complexe, Symmetric, DiagonalRow> D = vec_D_tau0_Bis(ref);
            if (D(i,i) != 0.0)
              Dtau2.Get(i,i) += (D.Get(i,i) + beta*beta)*(D.Get(i,i) + beta*beta);
              // Dtau2.Get(i,i) += (D.Get(i,i) + beta*beta)*(D.Get(i,i) + beta*beta) / nbRefNoeud.Get(i);
            // DISP(i+1);
            // cout << D(i,i) << ", beta : " << beta << endl;
            // DISP(Dtau2.Get(i,i));
            // DISP(nbRefNoeud.Get(i));
          }
    }
}


void ModeEs_Solver::ComputeEsHs(const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                const Vector<int>& nbRefNoeud,
                const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                const Vector<Matrix<Complexe, Prop, Storage>>& vec_Uh,
                const Vector<Matrix<Complexe, Prop, Storage>> &vec_Th,
                const Vector<Matrix<Complexe, Symmetric, DiagonalRow>> &vec_ktrSquared,
                const Complex_wp& L,
                const Vector<Complex_wp>& V2, const Vector<Complex_wp>& W2,
                Vector<Complex_wp>& Es, Vector<Complex_wp>& Hs)
{
  int N = this->mesh_num.GetNbDof();
  Vector<Complex_wp> V1(N);
  Vector<Complex_wp> W1(N);
  Matrix<Complexe, Symmetric, DiagonalRow> DTau;
  Matrix<Complexe, General, ArrayRowSparse> MM;
  Complex_wp beta = L * this->omega;

  ComputeMatrixV2ToV1(L, vec_D_tau0_Bis, nbRefNoeud, DTau);

  Es.Reallocate(N); Es.Zero();
  Hs.Reallocate(N); Hs.Zero();
  MM.Reallocate(N,N);

  Mlt(DTau, V2, V1);
  Mlt(DTau, W2, W1);

  for (int ref = 0; ref < this->ref_epsilon.GetM(); ref++)
    {
      if (vec_Mh(ref).GetM() > 0)
        {
          Complex_wp eps_sigma = -Iwp*this->omega*this->ref_epsilon(ref) + this->ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*this->omega*this->ref_mu(ref);
          Complex_wp k2 = -eps_sigma*mu_tilde;
          Complex_wp one(1,0);
          Vector<Complex_wp> V2t(N);
          Vector<Complex_wp> W2t(N);
          Mlt(vec_ktrSquared(ref), V2, V2t);
          Mlt(vec_ktrSquared(ref), W2, W2t);

          MltAdd(- k2 + beta * beta, vec_Mh(ref), V1, one, Es);
          MltAdd(- k2 + beta * beta, vec_Mh(ref), W1, one, Hs);
          MltAdd(- Iwp * beta * this->ref_tau(ref), vec_Uh(ref), V1, one, Es);
          MltAdd(- Iwp * beta * this->ref_tau(ref), vec_Uh(ref), W1, one, Hs);
          MltAdd(mu_tilde * this->ref_tau(ref), vec_Th(ref), W1, one, Es);
          MltAdd(- eps_sigma * this->ref_tau(ref), vec_Th(ref), V1, one, Hs);
          MltAdd(- 2.0 * mu_tilde * this->ref_tau(ref), vec_Mh(ref), W2t, one, Es);
          MltAdd(2.0 * eps_sigma * this->ref_tau(ref), vec_Mh(ref), V2t, one, Hs);

          Add(1.0, vec_Mh(ref), MM);
        }
    }
    SparseDistributedSolver<Complex_wp> mat_lu;
    mat_lu.Factorize(MM);
    mat_lu.Solve(Es);
    mat_lu.Solve(Hs);
}

void ModeEs_Solver::ComputeMatrixF2ToF1(const Complex_wp& L,
                          const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                          const Vector<int>& nbRefNoeud,
                          Matrix<Complexe, Symmetric, DiagonalRow>& Dtau)
{
  int N = this->mesh_num.GetNbDof();
  Complex_wp beta = L * this->omega;
  Dtau.Reallocate(N,N);

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  for (int i = 0 ; i < N ; i++)
    {
      Dtau.Get(i,i) = 0.0;
      for (int ref = 0; ref < index_used.GetM(); ref++)
        if (index_used(ref))
          {
            Complex_wp eps_sigma = -Iwp*this->omega*this->ref_epsilon(ref) + this->ref_sigma(ref);
            Complex_wp mu_tilde = -Iwp*this->omega*this->ref_mu(ref);
            Complex_wp k2 = -eps_sigma*mu_tilde;
            Matrix<Complexe, Symmetric, DiagonalRow> D = vec_D_tau0_Bis(ref);
            if (D(i,i) != 0.0)
              Dtau.Get(i,i) += ((D.Get(i,i) + beta*beta) / (-k2 + beta*beta)) / double(nbRefNoeud.Get(i));
              // Dtau.Get(i,i) += (D.Get(i,i) + beta*beta);
            // DISP(i+1);
            // cout << D(i,i) << ", beta : " << beta << endl;
            // DISP(Dtau2.Get(i,i));
            // DISP(nbRefNoeud.Get(i));
          }
    }
}

void ModeEs_Solver::FromF2G2ToEsHs(const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                const Vector<int>& nbRefNoeud,
                const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                const Vector<Matrix<Complexe, Prop, Storage>>& vec_Uh,
                const Vector<Matrix<Complexe, Prop, Storage>> &vec_Th,
                const Vector<Matrix<Complexe, Symmetric, DiagonalRow>> &vec_ktrSquared,
                const Complex_wp& L,
                const Vector<Complex_wp>& F2, const Vector<Complex_wp>& G2,
                Vector<Complex_wp>& Es, Vector<Complex_wp>& Hs)
{
  int N = this->mesh_num.GetNbDof();
  Matrix<Complexe, General, ArrayRowSparse> MM, MMe;
  Complex_wp beta = L * this->omega;

  // ComputeMatrixF2ToF1(L, vec_D_tau0_Bis, nbRefNoeud, DTau);

  Es.Reallocate(N); Es.Zero();
  Hs.Reallocate(N); Hs.Zero();
  MM.Reallocate(N,N);
  MMe.Reallocate(N,N);


  for (int ref = 0; ref < this->ref_epsilon.GetM(); ref++)
    {
      if (vec_Mh(ref).GetM() > 0)
        {
          Complex_wp eps_sigma = -Iwp*this->omega*this->ref_epsilon(ref) + this->ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*this->omega*this->ref_mu(ref);
          Complex_wp k2 = -eps_sigma*mu_tilde;
          Complex_wp DBeta = - k2 + beta * beta;
          Complex_wp one(1,0);
          Vector<Complex_wp> F1(F2);
          Vector<Complex_wp> G1(G2);
          Vector<Complex_wp> temp(N); temp.Zero();
          Vector<Complex_wp> F2t(N);
          Vector<Complex_wp> G2t(N);

          Mlt(vec_ktrSquared(ref), F2, F2t);
          Mlt(vec_ktrSquared(ref), G2, G2t);
          Add(-1.0 / DBeta, F2t, F1);
          Add(-1.0 / DBeta, G2t, G1);
          //
          // for (int i = 0 ; i < N ; i++)
          //   {
          //     if (vec_ktrSquared(ref)(i,i) == 0)
          //      {
          //        F1.Get(i) = 0.0;
          //      }
          //   }

          WriteOutputFileScalar(F1, "F1_"+ to_str(ref));
          WriteOutputFileScalar(G1, "G1_"+ to_str(ref));


          MltAdd(one, vec_Mh(ref), F1, one, Es);
          MltAdd(one, vec_Mh(ref), G1, one, Hs);
          MltAdd(- Iwp * beta * this->ref_tau(ref) / DBeta, vec_Uh(ref), F1, one, Es);
          MltAdd(- Iwp * beta * this->ref_tau(ref) / DBeta, vec_Uh(ref), G1, one, Hs);
          MltAdd(mu_tilde * this->ref_tau(ref) / DBeta, vec_Th(ref), G1, one, Es);
          MltAdd(- eps_sigma * this->ref_tau(ref) / DBeta, vec_Th(ref), F1, one, Hs);
          MltAdd(- 2.0 * mu_tilde * this->ref_tau(ref) / (DBeta * DBeta), vec_Mh(ref), G2t, one, Es);
          MltAdd(2.0 * eps_sigma * this->ref_tau(ref) / (DBeta * DBeta), vec_Mh(ref), F2t, one, Hs);

          Add(1.0, vec_Mh(ref), MM);
          Add(1.0, vec_Mh(ref), MMe);
        }
    }

    // IL FAUT REALLOUER TOUTE LA COLONNE !!!!
    for (int i = 0; i < N; i++)
      if (IndexDirichlet(i) == -1)
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
