// initialisation du probleme
ModeEsPoly1DirichletEigenProblem
::ModeEsPoly1DirichletEigenProblem(int N, ModeEs_Solver& var0,
                              Vector<Matrix<Complex_wp, Prop, Storage> >& Mh,
                              Vector<Matrix<Complex_wp, Prop, Storage> >& M1h,
                              Vector<Matrix<Complex_wp, Prop, Storage> >& M2h,
                              Vector<Matrix<Complex_wp, Prop, Storage> >& Mr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Rih,
                              Vector<Matrix<Complexe, Prop, Storage> >& Sih,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ti1h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ti2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Tj1h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Tj2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Di1h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Di2h,
                              Vector<int>& DofDir, Vector<int>& IndexDir)
  :
  var(var0), vec_Mh(Mh), vec_M1h(M1h), vec_M2h(M2h), vec_Mr2h(Mr2h), vec_Rih(Rih), vec_Sih(Sih), vec_Ti1h(Ti1h), vec_Ti2h(Ti2h), vec_Tj1h(Tj1h), vec_Tj2h(Tj2h), vec_Di1h(Di1h), vec_Di2h(Di2h), DofKeptDir(DofDir), IndexNonDirichlet(IndexDir)
{
    // pour Es, condition de Dirichlet
    nb_dof_E3 = DofKeptDir.GetM();

    // pour Hs, condition de Neumann
    nb_dof_H3 = N;

    // taille du pb aux valeurs propres
    nb_rows = nb_dof_E3;
    nb_dof_all = N;
    nb_rows += nb_dof_H3;
    nb_dof_all += N;

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

    DISP(nb_couches);

    for (int ref = 0 ; ref < nb_couches ; ref++)
      if (vec_Mh(ref).GetM() == N)
        {
          nb_rows += 5 * nb_dof_couche(ref) + nb_dof_coucheE3(ref);
        }

    DISP(nb_rows);
    this->Init(nb_rows);

    ComputeMatrix();
}

void ModeEsPoly1DirichletEigenProblem::RemoveDirichletRowCol(const Matrix<Complexe, Prop, Storage>& A,
                                                        Matrix<Complexe, Prop, Storage>& B)
{
  B.Clear();
  B.Reallocate(nb_dof_E3, nb_dof_E3);
  for (int i = 0; i < nb_dof_E3; i++)
    {
      int i0 = DofKeptDir(i);
      int nb_val = 0;
      for (int j = 0; j < A.GetRowSize(i0); j++)
        if (IndexNonDirichlet(A.Index(i0, j)) >= 0)
          nb_val++;

      B.ReallocateRow(i, nb_val);
      nb_val = 0;
      for (int j = 0; j < A.GetRowSize(i0); j++)
        if (IndexNonDirichlet(A.Index(i0, j)) >= 0)
          {
            B.Index(i, nb_val) = IndexNonDirichlet(A.Index(i0, j));
            B.Value(i, nb_val) = A.Value(i0, j);
            nb_val++;
          }
    }
}

void ModeEsPoly1DirichletEigenProblem::ExtractDirichlet(const VectComplex_wp& Xdir, VectComplex_wp& X)
{
  X.Zero();
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    X(DofKeptDir(i)) = Xdir(i);

  for (int i = 0; i < nb_dof_H3; i++)
    X(nb_dof_H3 + i) = Xdir(nb_dof_E3 + i);
}

void ModeEsPoly1DirichletEigenProblem::ExpandDirichlet(const VectComplex_wp& X, VectComplex_wp& Xdir)
{
  for (int i = 0; i < DofKeptDir.GetM(); i++)
    Xdir(i) = X(DofKeptDir(i));

  for (int i = 0; i < nb_dof_H3; i++)
    Xdir(nb_dof_E3 + i) = X(nb_dof_H3 + i);
}



void ModeEsPoly1DirichletEigenProblem::ComputeA0()
{
  cout << "Création de la matrice A0" << endl;
  A0.Reallocate(nb_rows, nb_rows);
  A0.Zero();

  int inc = nb_dof_E3 + nb_dof_H3;
  for (int ref = 0; ref < nb_couches; ref++)
    {
      if (nb_dof_couche(ref) > 0)
        {
          Complex_wp epsilon = var.ref_epsilon(ref);
          Complex_wp mu = var.ref_mu(ref);
          Complex_wp sigma = var.ref_sigma(ref);
          Real_wp tau = var.torsion;

          Complex_wp eps_tilde = - Iwp * epsilon * var.omega + sigma;
          Complex_wp mu_tilde = - Iwp * mu * var.omega;
          Complex_wp k2L = - eps_tilde * mu_tilde;

          int inc_F0 = 0;
          int inc_F1 = nb_dof_coucheE3(ref);
          int inc_F2 = nb_dof_coucheE3(ref) + nb_dof_couche(ref); // les bords de F1 doivent être pris en compte ! E3 = 0 sur le bord =/=> grad E3 = 0
          int inc_G0 = nb_dof_coucheE3(ref) + 2*nb_dof_couche(ref);
          int inc_G1 = nb_dof_coucheE3(ref) + 3*nb_dof_couche(ref);
          int inc_G2 = nb_dof_coucheE3(ref) + 4*nb_dof_couche(ref);

          Matrix<Complexe, Prop, Storage> Mh = vec_Mh(ref);
          Matrix<Complexe, Prop, Storage> M1h = vec_M1h(ref);
          Matrix<Complexe, Prop, Storage> M2h = vec_M2h(ref);
          Matrix<Complexe, Prop, Storage> Mr2h = vec_Mr2h(ref);
          Matrix<Complexe, Prop, Storage> Rih = vec_Rih(ref);
          Matrix<Complexe, Prop, Storage> Ti1h = vec_Ti1h(ref);
          Matrix<Complexe, Prop, Storage> Ti2h = vec_Ti2h(ref);
          Matrix<Complexe, Prop, Storage> Tj1h = vec_Tj1h(ref);
          Matrix<Complexe, Prop, Storage> Tj2h = vec_Tj2h(ref);
          Matrix<Complexe, Prop, Storage> Di1h = vec_Di1h(ref);
          Matrix<Complexe, Prop, Storage> Di2h = vec_Di2h(ref);

          for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
            {
              int m = NumGlob(ref)(ddl_loc); // numéro global avec les noeuds de Dirichlet
              int m_loc = IndexDof(ref)(m);
              int m_loc_E3 = IndexDofE3(ref)(m);

              int size_row = Mh.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Mh.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_E3 = IndexDofE3(ref)(n);

                  // variables auxiliaires
                  if ((IndexNonDirichlet(m) != -1) && (IndexNonDirichlet(n) != -1))
                    {
                      A0.Get(inc + inc_F0 + m_loc_E3, IndexNonDirichlet(n)) += Mh.Value(m, j);
                      A0.Get(inc + inc_F0 + m_loc_E3, inc + inc_F0 + n_loc_E3) += k2L * Mh.Value(m, j);
                    }
                  A0.Get(inc + inc_G0 + m_loc, nb_dof_E3 + n) += Mh.Value(m, j);
                  A0.Get(inc + inc_G0 + m_loc, inc + inc_G0 + n_loc) += k2L * Mh.Value(m, j);

                  A0.Get(inc + inc_F1 + m_loc, inc + inc_F1 + n_loc) += k2L * Mh.Value(m, j);
                  A0.Get(inc + inc_F2 + m_loc, inc + inc_F2 + n_loc) += k2L * Mh.Value(m, j);

                  A0.Get(inc + inc_G1 + m_loc, inc + inc_G1 + n_loc) += k2L * Mh.Value(m, j);
                  A0.Get(inc + inc_G2 + m_loc, inc + inc_G2 + n_loc) += k2L * Mh.Value(m, j);

                  // équation principale
                  if ((IndexNonDirichlet(m) != -1) && (IndexNonDirichlet(n) != -1))
                    A0.Get(IndexNonDirichlet(m), inc + inc_F0 + n_loc_E3) += - eps_tilde * k2L * Mh.Value(m, j);
                  A0.Get(nb_dof_E3 + m, inc + inc_G0 + n_loc) += - mu_tilde * k2L * Mh.Value(m, j);
                }

              size_row = Mr2h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Mr2h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_E3 = IndexDofE3(ref)(n);

                  if ((IndexNonDirichlet(m) != -1) && (IndexNonDirichlet(n) != -1))
                    A0.Get(inc + inc_F0 + m_loc_E3, inc + inc_F0 + n_loc_E3) += tau * tau * k2L * Mr2h.Value(m, j);
                  A0.Get(inc + inc_G0 + m_loc, inc + inc_G0 + n_loc) += tau * tau * k2L * Mr2h.Value(m, j);

                  A0.Get(inc + inc_F1 + m_loc, inc + inc_F1 + n_loc) += tau * tau * k2L * Mr2h.Value(m, j);
                  A0.Get(inc + inc_F2 + m_loc, inc + inc_F2 + n_loc) += tau * tau * k2L * Mr2h.Value(m, j);

                  A0.Get(inc + inc_G1 + m_loc, inc + inc_G1 + n_loc) += tau * tau * k2L * Mr2h.Value(m, j);
                  A0.Get(inc + inc_G2 + m_loc, inc + inc_G2 + n_loc) += tau * tau * k2L * Mr2h.Value(m, j);
                }

              size_row = Tj1h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Tj1h.Index(m, j);
                  if (IndexNonDirichlet(n) != -1)
                    A0.Get(inc + inc_F1 + m_loc, IndexNonDirichlet(n)) += Tj1h.Value(m, j);
                  A0.Get(inc + inc_G1 + m_loc, nb_dof_E3 + n) += Tj1h.Value(m, j);
                }

              size_row = Tj2h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Tj2h.Index(m, j);
                  if (IndexNonDirichlet(n) != -1)
                    A0.Get(inc + inc_F2 + m_loc, IndexNonDirichlet(n)) += Tj2h.Value(m, j);
                  A0.Get(inc + inc_G2 + m_loc, nb_dof_E3 + n) += Tj2h.Value(m, j);
                }

              size_row = M1h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = M1h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A0.Get(IndexNonDirichlet(m), inc + inc_G1 + n_loc) += - k2L * tau * M1h.Value(m, j);
                  A0.Get(nb_dof_E3 + m, inc + inc_F1 + n_loc) += k2L * tau * M1h.Value(m, j);
                }

              size_row = M2h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = M2h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A0.Get(IndexNonDirichlet(m), inc + inc_G2 + n_loc) += - k2L * tau * M2h.Value(m, j);
                  A0.Get(nb_dof_E3 + m, inc + inc_F2 + n_loc) += k2L * tau * M2h.Value(m, j);
                }

              size_row = Rih.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Rih.Index(m, j);
                  int n_loc = IndexDof(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A0.Get(IndexNonDirichlet(m), inc + inc_G0 + n_loc) += - k2L * tau * Rih.Value(m, j);
                  A0.Get(nb_dof_E3 + m, inc + inc_F0 + n_loc) += k2L * tau * Rih.Value(m, j);
                }

              size_row = Ti1h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Ti1h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A0.Get(IndexNonDirichlet(m), inc + inc_F1 + n_loc) += eps_tilde * Ti1h.Value(m, j);
                  A0.Get(nb_dof_E3 + m, inc + inc_G1 + n_loc) += mu_tilde * Ti1h.Value(m, j);
                }

              size_row = Ti2h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Ti2h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A0.Get(IndexNonDirichlet(m), inc + inc_F2 + n_loc) += eps_tilde * Ti2h.Value(m, j);
                  A0.Get(nb_dof_E3 + m, inc + inc_G2 + n_loc) += mu_tilde * Ti2h.Value(m, j);
                }

              size_row = Di1h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Di1h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A0.Get(IndexNonDirichlet(m), inc + inc_F2 + n_loc) += eps_tilde * tau * tau * Di1h.Value(m, j);
                  A0.Get(nb_dof_E3 + m, inc + inc_G2 + n_loc) += mu_tilde * tau * tau * Di1h.Value(m, j);
                }

              size_row = Di2h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Di2h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A0.Get(IndexNonDirichlet(m), inc + inc_F1 + n_loc) += - eps_tilde * tau * tau * Di2h.Value(m, j);
                  A0.Get(nb_dof_E3 + m, inc + inc_G1 + n_loc) += - mu_tilde * tau * tau * Di2h.Value(m, j);
                }
            }
          inc += nb_dof_coucheE3(ref) + 5 * nb_dof_couche(ref); //on arrive au premier indice de la couche suivante dans la matrice
        }
    }
}

void ModeEsPoly1DirichletEigenProblem::ComputeA1()
{
  cout << "Création de la matrice A1" << endl;
  A1.Reallocate(nb_rows, nb_rows);
  A1.Zero();

  int inc = nb_dof_E3 + nb_dof_H3;
  for (int ref = 0; ref < nb_couches; ref++)
    {
      if (nb_dof_couche(ref) > 0)
        {
          Complex_wp epsilon = var.ref_epsilon(ref);
          Complex_wp mu = var.ref_mu(ref);
          Complex_wp sigma = var.ref_sigma(ref);
          Real_wp tau = var.torsion;

          Complex_wp eps_tilde = - Iwp * epsilon * var.omega + sigma;
          Complex_wp mu_tilde = - Iwp * mu * var.omega;
          Complex_wp k2L = - eps_tilde * mu_tilde;

          int inc_F0 = 0;
          int inc_F1 = nb_dof_coucheE3(ref);
          int inc_F2 = nb_dof_coucheE3(ref) + nb_dof_couche(ref); // les bords de F1 doivent être pris en compte ! E3 = 0 sur le bord =/=> grad E3 = 0
          int inc_G0 = nb_dof_coucheE3(ref) + 2*nb_dof_couche(ref);
          int inc_G1 = nb_dof_coucheE3(ref) + 3*nb_dof_couche(ref);
          int inc_G2 = nb_dof_coucheE3(ref) + 4*nb_dof_couche(ref);

          Matrix<Complexe, Prop, Storage> M1h = vec_M1h(ref);
          Matrix<Complexe, Prop, Storage> M2h = vec_M2h(ref);
          Matrix<Complexe, Prop, Storage> Sih = vec_Sih(ref);
          Matrix<Complexe, Prop, Storage> Ti1h = vec_Ti1h(ref);
          Matrix<Complexe, Prop, Storage> Ti2h = vec_Ti2h(ref);

          for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
            {
              int m = NumGlob(ref)(ddl_loc); // numéro global avec les noeuds de Dirichlet

              int size_row = M1h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = M1h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A1.Get(IndexNonDirichlet(m), inc + inc_F2 + n_loc) += Iwp * var.omega * eps_tilde * tau * M1h.Value(m, j);
                  A1.Get(nb_dof_E3 + m, inc + inc_G2 + n_loc) += Iwp * var.omega * mu_tilde * tau * M1h.Value(m, j);
                }

              size_row = M2h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = M2h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A1.Get(IndexNonDirichlet(m), inc + inc_F1 + n_loc) += - Iwp * var.omega * eps_tilde * tau * M2h.Value(m, j);
                  A1.Get(nb_dof_E3 + m, inc + inc_G1 + n_loc) += - Iwp * var.omega * mu_tilde * tau * M2h.Value(m, j);
                }

              size_row = Sih.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Sih.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_E3 = IndexDofE3(ref)(n);

                  // équation principale
                  if ((IndexNonDirichlet(m) != -1) && (IndexNonDirichlet(n) != -1))
                    A1.Get(IndexNonDirichlet(m), inc + inc_F0 + n_loc_E3) += - Iwp * var.omega * eps_tilde * tau * Sih.Value(m, j);
                  A1.Get(nb_dof_E3 + m, inc + inc_G0 + n_loc) += - Iwp * var.omega * mu_tilde * tau * Sih.Value(m, j);
                }

              size_row = Ti1h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Ti1h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_E3 = IndexDofE3(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A1.Get(IndexNonDirichlet(m), inc + inc_G2 + n_loc) += - Iwp * var.omega * Ti1h.Value(m, j);
                  A1.Get(nb_dof_E3 + m, inc + inc_F2 + n_loc) += Iwp * var.omega * Ti1h.Value(m, j);
                }

              size_row = Ti2h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Ti2h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_E3 = IndexDofE3(ref)(n);

                  // équation principale
                  if (IndexNonDirichlet(m) != -1)
                    A1.Get(IndexNonDirichlet(m), inc + inc_G1 + n_loc) += Iwp * var.omega * Ti2h.Value(m, j);
                  A1.Get(nb_dof_E3 + m, inc + inc_F1 + n_loc) += - Iwp * var.omega * Ti2h.Value(m, j);
                }
            }
          inc += nb_dof_coucheE3(ref) + 5 * nb_dof_couche(ref); //on arrive au premier indice de la couche suivante dans la matrice
        }
    }
}

void ModeEsPoly1DirichletEigenProblem::ComputeA2()
{
  cout << "Création de la matrice A2" << endl;
  A2.Reallocate(nb_rows, nb_rows);
  A2.Zero();

  int inc = nb_dof_E3 + nb_dof_H3;
  for (int ref = 0; ref < nb_couches; ref++)
    {
      if (nb_dof_couche(ref) > 0)
        {
          Complex_wp epsilon = var.ref_epsilon(ref);
          Complex_wp mu = var.ref_mu(ref);
          Complex_wp sigma = var.ref_sigma(ref);
          Real_wp tau = var.torsion;

          Complex_wp eps_tilde = - Iwp * epsilon * var.omega + sigma;
          Complex_wp mu_tilde = - Iwp * mu * var.omega;
          Complex_wp k2L = - eps_tilde * mu_tilde;

          int inc_F0 = 0;
          int inc_F1 = nb_dof_coucheE3(ref);
          int inc_F2 = nb_dof_coucheE3(ref) + nb_dof_couche(ref); // les bords de F1 doivent être pris en compte ! E3 = 0 sur le bord =/=> grad E3 = 0
          int inc_G0 = nb_dof_coucheE3(ref) + 2*nb_dof_couche(ref);
          int inc_G1 = nb_dof_coucheE3(ref) + 3*nb_dof_couche(ref);
          int inc_G2 = nb_dof_coucheE3(ref) + 4*nb_dof_couche(ref);

          Matrix<Complexe, Prop, Storage> Mh = vec_Mh(ref);

          for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
            {
              int m = NumGlob(ref)(ddl_loc); // numéro global avec les noeuds de Dirichlet
              int m_loc = IndexDof(ref)(m);
              int m_loc_E3 = IndexDofE3(ref)(m);

              int size_row = Mh.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Mh.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_E3 = IndexDofE3(ref)(n);

                  // variables auxiliaires
                  if ((IndexNonDirichlet(m) != -1) && (IndexNonDirichlet(n) != -1))
                    A2.Get(inc + inc_F0 + m_loc_E3, inc + inc_F0 + n_loc_E3) += - var.omega * var.omega * Mh.Value(m, j);
                  A2.Get(inc + inc_G0 + m_loc, inc + inc_G0 + n_loc) += - var.omega * var.omega * Mh.Value(m, j);

                  A2.Get(inc + inc_F1 + m_loc, inc + inc_F1 + n_loc) += - var.omega * var.omega * Mh.Value(m, j);
                  A2.Get(inc + inc_F2 + m_loc, inc + inc_F2 + n_loc) += - var.omega * var.omega * Mh.Value(m, j);

                  A2.Get(inc + inc_G1 + m_loc, inc + inc_G1 + n_loc) += - var.omega * var.omega * Mh.Value(m, j);
                  A2.Get(inc + inc_G2 + m_loc, inc + inc_G2 + n_loc) += - var.omega * var.omega * Mh.Value(m, j);

                  // équation principale
                  if ((IndexNonDirichlet(m) != -1) && (IndexNonDirichlet(n) != -1))
                    A2.Get(IndexNonDirichlet(m), inc + inc_F0 + n_loc_E3) += var.omega * var.omega * eps_tilde * Mh.Value(m, j);
                  A2.Get(nb_dof_E3 + m, inc + inc_G0 + n_loc) += var.omega * var.omega * mu_tilde * Mh.Value(m, j);
                }
            }
          inc += nb_dof_coucheE3(ref) + 5 * nb_dof_couche(ref); //on arrive au premier indice de la couche suivante dans la matrice
        }
    }
}


void ModeEsPoly1DirichletEigenProblem::ComputeMatrix()
{
  ComputeA0();
  ComputeA1();
  ComputeA2();
}

void ModeEsPoly1DirichletEigenProblem::ComputeEsHs(const Complex_wp& L,
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

  Vector<Complex_wp> grE3vXi;
  Vector<Complex_wp> grH3sXi;

  Es.Reallocate(N); Es.Zero();
  Hs.Reallocate(N); Hs.Zero();
  grE3vXi.Reallocate(N); grE3vXi.Zero();
  grH3sXi.Reallocate(N); grH3sXi.Zero();
  E1.Reallocate(N); Es.Zero();
  E2.Reallocate(N); Hs.Zero();
  H1.Reallocate(N); Es.Zero();
  H2.Reallocate(N); Hs.Zero();
  MM.Reallocate(N,N); MM.Zero();
  MMe.Reallocate(N,N); MMe.Zero();

  // On va résoudre MMe Es = EE (potentiellement en rajoutant les zéros), et MM Hs = HH
  //  La boucle construit la matrice MM, et les vecteurs EE et HH
  int inc = nb_dof_E3 + nb_dof_H3;
  for (int ref = 0; ref < nb_couches; ref++)
    {
      if (nb_dof_couche(ref) > 0)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp mu_tilde = -Iwp*var.omega*var.ref_mu(ref);
          Complex_wp DBeta = eps_sigma*mu_tilde + beta * beta;

          int inc_F0 = 0;
          int inc_F1 = nb_dof_coucheE3(ref);
          int inc_F2 = nb_dof_coucheE3(ref) + nb_dof_couche(ref); // les bords de F1 doivent être pris en compte ! E3 = 0 sur le bord =/=> grad E3 = 0
          int inc_G0 = nb_dof_coucheE3(ref) + 2*nb_dof_couche(ref);
          int inc_G1 = nb_dof_coucheE3(ref) + 3*nb_dof_couche(ref);
          int inc_G2 = nb_dof_coucheE3(ref) + 4*nb_dof_couche(ref);

          Matrix<Complexe, Prop, Storage> Mh = vec_Mh(ref);
          Matrix<Complexe, Prop, Storage> M1h = vec_M1h(ref);
          Matrix<Complexe, Prop, Storage> M2h = vec_M2h(ref);


          for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
            {
              int m = NumGlob(ref)(ddl_loc); // numéro global avec les noeuds de Dirichlet

              int size_row = Mh.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Mh.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_E3 = IndexDofE3(ref)(n);

                  if (IndexNonDirichlet(n) != -1)
                    Es(m) += DBeta * Mh.Value(m, j) * eigenvector(inc + inc_F0 + n_loc_E3);
                  Hs(m) += DBeta * Mh.Value(m, j) * eigenvector(inc + inc_G0 + n_loc);
                }

              size_row = M1h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = M1h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_E3 = IndexDofE3(ref)(n);

                  Es(m) += mu_tilde * var.torsion * M1h.Value(m, j) * eigenvector(inc + inc_G1 + n_loc);
                  Es(m) += Iwp * beta * var.torsion * M1h.Value(m, j) * eigenvector(inc + inc_F2 + n_loc);
                  grH3sXi(m) += mu_tilde * var.torsion * M1h.Value(m, j) *  eigenvector(inc + inc_G1 + n_loc);
                  grE3vXi(m) += Iwp * beta * var.torsion * M1h.Value(m, j) * eigenvector(inc + inc_F2 + n_loc);
                  Hs(m) += - eps_sigma * var.torsion * M1h.Value(m, j) * eigenvector(inc + inc_F1 + n_loc);
                  Hs(m) += Iwp * beta * var.torsion * M1h.Value(m, j) * eigenvector(inc + inc_G2 + n_loc);
                }

              size_row = M2h.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = M2h.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_E3 = IndexDofE3(ref)(n);

                  Es(m) += mu_tilde * var.torsion * M2h.Value(m, j) * eigenvector(inc + inc_G2 + n_loc);
                  Es(m) += - Iwp * beta * var.torsion * M2h.Value(m, j) * eigenvector(inc + inc_F1 + n_loc);
                  grH3sXi(m) += mu_tilde * var.torsion * M2h.Value(m, j) *  eigenvector(inc + inc_G2 + n_loc);
                  grE3vXi(m) += - Iwp * beta * var.torsion * M2h.Value(m, j) * eigenvector(inc + inc_F1 + n_loc);
                  Hs(m) += - eps_sigma * var.torsion * M2h.Value(m, j) * eigenvector(inc + inc_F2 + n_loc);
                  Hs(m) += - Iwp * beta * var.torsion * M2h.Value(m, j) * eigenvector(inc + inc_G1 + n_loc);
                }
            }

          Add(1.0, vec_Mh(ref), MM);
          Add(1.0, vec_Mh(ref), MMe);

          inc += nb_dof_coucheE3(ref) + 5 * nb_dof_couche(ref); //on arrive au premier indice de la couche suivante dans la matrice
        }
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

    mat_lu.Solve(grE3vXi);
    mat_lu.Solve(grH3sXi);

    var.WriteOutputFileScalar(grE3vXi, "GrE3xi");
    var.WriteOutputFileScalar(grH3sXi, "GrH3xi");
}
