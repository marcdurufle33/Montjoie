void ModeEs_Solver::AssemblageDirichletMatricesVerifDroit(Matrix<Complex_wp, General, ArrayRowSparse>& K_tilde,
                        Matrix<Complex_wp, General, ArrayRowSparse>& M_tilde,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch,
                        const Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                        const Vector<Matrix<Complexe, Prop, Storage> >& vec_Th,
                        const Vector<Matrix<Complexe, Prop, Storage> >& vec_Uh,
                        const VectR2 vec_Pos)
{
  // construction des numerotations
  int N = this->mesh_num.GetNbDof();
  int nb_dof_Es = DofKeptDir.GetM();
  int nb_dof_Hs = N;
  cout << "N = " << N << endl;
  DISP(vec_Pos.GetM());
  int nb_couches = vec_Mh.GetM();
  Vector<int> nb_dof_couche(nb_couches), nb_dof_coucheEs(nb_couches);
  nb_dof_couche.Zero(); nb_dof_coucheEs.Zero();
  Vector<Vector<int> > IndexDofEs(nb_couches), NumGlobEs(nb_couches);
  Vector<Vector<int> > IndexDof(nb_couches), NumGlob(nb_couches);
  for (int ref = 0; ref < nb_couches; ref++)
    if (vec_Mh(ref).GetM() == N)
      {
        IndexDof(ref).Reallocate(N);
        IndexDof(ref).Fill(-1);

        IndexDofEs(ref).Reallocate(N);
        IndexDofEs(ref).Fill(-1);
      }

  // IndexDof(ref)(j) renvoie le numero local du ddl global j dans la couche ref
  // NumGlob(ref)(i) renvoie le numero global j du ddl local i dans la couche ref
  // IndexDofEs, NumGlobEz pour l'inconnue Es (ddl Dirichlets enleves)
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    {
      int ref = this->mesh.Element(i).GetReference();
      int nb_dof_loc = this->mesh_num.Element(i).GetNbDof();
      for (int j = 0; j < nb_dof_loc; j++)
        {
          int num_dof = this->mesh_num.Element(i).GetNumberDof(j);
          if (IndexDof(ref)(num_dof) == -1) // si c'est un nouveau
            {
              IndexDof(ref)(num_dof) = nb_dof_couche(ref);
              NumGlob(ref).PushBack(num_dof);
              nb_dof_couche(ref)++;
            }

          if ((IndexDirichlet(num_dof) >= 0) && (IndexDofEs(ref)(num_dof) == -1)) // si c'est un nouveau
            {
              IndexDofEs(ref)(num_dof) = nb_dof_coucheEs(ref);
              NumGlobEs(ref).PushBack(num_dof);
              nb_dof_coucheEs(ref)++;
            }
        }
    }


  // Initialisation
  int ordreMatrices = nb_dof_Es + nb_dof_Hs; // lignes pour Es et Hs
  /* lignes pour V_ref et W_ref et leurs variables temporaires
  (pour outrepasser le caractère quadratique de l'équation) */
  for (int ref = 0 ; ref < nb_couches ; ref++)
    ordreMatrices += 4 * nb_dof_couche(ref) + 4 * nb_dof_coucheEs(ref);


  K_tilde.Reallocate(ordreMatrices, ordreMatrices);
  M_tilde.Reallocate(ordreMatrices, ordreMatrices);
  DISP(ordreMatrices);
  K_tilde.Zero(); M_tilde.Zero();

  DISP(this->ref_tau);
  cout << "Création des matrices" << endl;
  int inc = nb_dof_Es + nb_dof_Hs; // premier indice de la première couche dans la matrice
  for (int ref = 0 ; ref < nb_couches ; ref++)
    if (nb_dof_couche(ref) > 0)
      {
        Complex_wp epsilon = this->ref_epsilon(ref);
        Complex_wp mu = this->ref_mu(ref);
        Complex_wp sigma = this->ref_sigma(ref);
        Complex_wp eps_tilde = - Iwp * this->omega * epsilon + sigma;
        Complex_wp mu_tilde = - Iwp * this->omega * mu;
        Complex_wp k2 = this->omega * this->omega * epsilon * mu + Iwp * this->omega * sigma * mu;
        Real_wp tau = this->ref_tau(ref);
        // conversion en non-symetrique
        Matrix<Complexe, Prop, Storage> Kh, Mh, Ch, Dh, Th, Uh;
        // Matrix<Real_wp, General, ArrayRowSparse> Kh, Mh;
        Copy(vec_Kh(ref), Kh);
        Copy(vec_Mh(ref), Mh);
        Copy(vec_Ch(ref), Ch);
        Copy(vec_Dh(ref), Dh);
        Copy(vec_Th(ref), Th);
        Copy(vec_Uh(ref), Uh);


        int inc_tempV1 = 0;
        int inc_V1 = nb_dof_coucheEs(ref);
        int inc_tempV2 = 2*nb_dof_coucheEs(ref);
        int inc_V2 = 3*nb_dof_coucheEs(ref);
        int inc_tempW1 = 4 * nb_dof_coucheEs(ref);
        int inc_W1 = 4 * nb_dof_coucheEs(ref) + nb_dof_couche(ref);
        int inc_tempW2 = 4 * nb_dof_coucheEs(ref) + 2*nb_dof_couche(ref);
        int inc_W2 = 4 * nb_dof_coucheEs(ref) + 3*nb_dof_couche(ref);

        // int inc_tempV1 = nb_dof_coucheEs(ref) + nb_dof_couche(ref);
        // int inc_V1 = 0;
        // int inc_W1 = nb_dof_coucheEs(ref);
        // int inc_tempW1 = 2 * nb_dof_coucheEs(ref) + nb_dof_couche(ref);

        for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
          {
            int m = NumGlob(ref)(ddl_loc); // numéro global avec les noeuds de Dirichlet
            int ddl_loc_Es = IndexDofEs(ref)(m);
            R2 pos = vec_Pos(m);
            Real_wp r2 = pos(0)*pos(0) + pos(1)*pos(1);// à compléter
            if (IndexDirichlet(m) != -1) // si on n'est pas sur un noeud de Dirichlet
              {
                // équation avec la variable temporaire : V_11 = beta omega V_1
                K_tilde.Get(inc + inc_tempV1 + ddl_loc_Es, inc + inc_tempV1 + ddl_loc_Es) += 1.0;
                M_tilde.Get(inc + inc_tempV1 + ddl_loc_Es, inc + inc_V1 + ddl_loc_Es) += this->omega;
                // équation avec la variable temporaire : V_21 = beta omega V_2
                K_tilde.Get(inc + inc_tempV2 + ddl_loc_Es, inc + inc_tempV2 + ddl_loc_Es) += 1.0;
                M_tilde.Get(inc + inc_tempV2 + ddl_loc_Es, inc + inc_V2 + ddl_loc_Es) += this->omega;


                // equation en Es : - Delta_tau,k V_1 + Es = 0
                //                   k2(1+tau^2 r^2) V_1 + Es = beta^2 omega^2 V1 = beta omega V_11
                // partie Ez dans l'equation en V_1
                K_tilde.Get(inc + inc_V1 + ddl_loc_Es, IndexDirichlet(m)) += 1.0;
                // partie k2(1+tau^2 r^2) V_1 dans l'equation en V_1
                K_tilde.Get(inc + inc_V1 + ddl_loc_Es, inc + inc_V1 + ddl_loc_Es) += k2 * (1 + tau * tau * r2);
                // partie beta omega V_11 dans l'equation en V1
                M_tilde.Get(inc + inc_V1 + ddl_loc_Es, inc + inc_tempV1 + ddl_loc_Es) += this->omega;

                // equation en V2 : - Delta_tau,k V_2 + V1 = 0
                //                   k2(1+tau^2 r^2) V_2 + V1 = beta^2 omega^2 V2 = beta omega V_21
                // partie Es dans l'equation en V_2
                K_tilde.Get(inc + inc_V2 + ddl_loc_Es, inc + inc_V1 + ddl_loc_Es) += 1.0;
                // partie k2(1+tau^2 r^2) V_2 dans l'equation en V_2
                K_tilde.Get(inc + inc_V2 + ddl_loc_Es, inc + inc_V2 + ddl_loc_Es) += k2 * (1 + tau * tau * r2);
                // partie beta omega V_21 dans l'equation en V2
                M_tilde.Get(inc + inc_V2 + ddl_loc_Es, inc + inc_tempV2 + ddl_loc_Es) += this->omega;
              }

            // équation avec la variable temporaire : W_11 = beta omega W_1
            K_tilde.Get(inc + inc_tempW1 + ddl_loc, inc + inc_tempW1 + ddl_loc) += 1.0;
            M_tilde.Get(inc + inc_tempW1 + ddl_loc, inc + inc_W1 + ddl_loc) += this->omega;
            // équation avec la variable temporaire : W_21 = beta omega W_2
            K_tilde.Get(inc + inc_tempW2 + ddl_loc, inc + inc_tempW2 + ddl_loc) += 1.0;
            M_tilde.Get(inc + inc_tempW2 + ddl_loc, inc + inc_W2 + ddl_loc) += this->omega;


            // equation en Hs : - Delta_tau,k W_1 + Hs = 0
            //                   k2(1+tau^2 r^2) W_1 + Hs = beta^2 omega^2 W1 = beta omega W_11
            // partie Hs dans l'equation en W1
            K_tilde.Get(inc + inc_W1 + ddl_loc, nb_dof_Es + m) += 1.0;
            // partie k2(1+tau^2 r^2) W_1 dans l'equation en W1
            K_tilde.Get(inc + inc_W1 + ddl_loc, inc + inc_W1 + ddl_loc) += k2 * (1 + tau * tau * r2);
            // partie beta omega W_11 dans l'equation en W1
            M_tilde.Get(inc + inc_W1 + ddl_loc, inc + inc_tempW1 + ddl_loc) += this->omega;



            // equation en W2 : - Delta_tau,k W_2 + W1 = 0
            //                   k2(1+tau^2 r^2) W_2 + W1 = beta^2 omega^2 W2 = beta omega W_21
            // partie W1 dans l'equation en W2
            K_tilde.Get(inc + inc_W2 + ddl_loc, inc + inc_W1 + ddl_loc) += 1.0;
            // partie k2(1+tau^2 r^2) W_2 dans l'equation en W2
            K_tilde.Get(inc + inc_W2 + ddl_loc, inc + inc_W2 + ddl_loc) += k2 * (1.0 + tau * tau * r2);
            // partie beta omega W22 dans l'equation en W2
            M_tilde.Get(inc + inc_W2 + ddl_loc, inc + inc_tempW2 + ddl_loc) += this->omega;



            // équations principales (en V1 et V2)
            // (- i omega epsilon + sigma) M_h E_z + (-i omega epsilon + sigma) K_h V_k = i beta omega C_h W_k
            int size_row = Kh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Kh.Index(m, j);
                int n_loc = IndexDof(ref)(n);
                int n_loc_Es = IndexDofEs(ref)(n);

                if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                  K_tilde.Get(IndexDirichlet(m), inc + inc_V1 + n_loc_Es) += eps_tilde * Kh.Value(m, j);
                K_tilde.Get(nb_dof_Es + m, inc + inc_W1 + n_loc) += mu_tilde * Kh.Value(m, j);
                if (m == 7)
                { DISP(ref); DISP(vec_Kh(ref).Index(m, j)); DISP(n_loc+inc); DISP( K_tilde.Get(m, n_loc + inc)); DISP(Kh.Value(m, j)); }
              }

            size_row = Mh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Mh.Index(m, j);
                int n_loc = IndexDof(ref)(n);
                int n_loc_Es = IndexDofEs(ref)(n);
                if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                  {
                    // - eps_tilde k^2 M_h V1 = - eps_tilde beta^2 omega^2 M_h V1 = - eps_tilde beta omega M_h V11
                    K_tilde.Get(IndexDirichlet(m), inc + inc_V1 + n_loc_Es) += - eps_tilde * k2 * Mh.Value(m, j);
                    M_tilde.Get(IndexDirichlet(m), inc + inc_tempV1 + n_loc_Es) += - eps_tilde * this->omega * Mh.Value(m, j);
                    // - 2 k^4 tau M_h W2 = - 2 k^2 beta^2 omega^2 tau M_h W2 =  2 k^2 beta omega tau M_h W21
                    K_tilde.Get(IndexDirichlet(m), inc + inc_W2 + n_loc) += - 2.0 * k2 * k2 * tau * Mh.Value(m, j);
                    M_tilde.Get(IndexDirichlet(m), inc + inc_tempW2 + n_loc) += - 2.0 * k2 * tau * this->omega * Mh.Value(m, j);
                    // 4 k^4 tau^2 eps_tilde Mh V2 = 0
                    K_tilde.Get(IndexDirichlet(m), inc + inc_V2 + n_loc_Es) += 4.0 * k2 * tau * tau * eps_tilde * Mh.Value(m, j);

                    K_tilde.Get(nb_dof_Es + m, inc + inc_V2 + n_loc_Es) += 2.0 * k2 * k2 * tau * Mh.Value(m, j);
                    M_tilde.Get(nb_dof_Es + m, inc + inc_tempV2 + n_loc_Es) += 2.0 * k2 * tau * this->omega * Mh.Value(m, j);
                  }
                K_tilde.Get(nb_dof_Es + m, inc + inc_W1 + n_loc) += - mu_tilde * k2 * Mh.Value(m, j);
                M_tilde.Get(nb_dof_Es + m, inc + inc_tempW1 + n_loc) += - mu_tilde * this->omega * Mh.Value(m, j);
                K_tilde.Get(nb_dof_Es + m, inc + inc_W2 + n_loc) += 4.0 * k2 * tau * tau * mu_tilde * Mh.Value(m, j);
              }

            size_row = Ch.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Ch.Index(m, j);
                int n_loc = IndexDof(ref)(n);
                int n_loc_Es = IndexDofEs(ref)(n);

                if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                  {
                    M_tilde.Get(IndexDirichlet(m), n_loc + inc_W1 + inc) += - Iwp * this->omega * Ch.Value(m, j);
                    M_tilde.Get(nb_dof_Es + m, n_loc_Es + inc_V1 + inc) += Iwp * this->omega * Ch.Value(m, j);
                  }
              }

            size_row = Uh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Uh.Index(m, j);
                int n_loc = IndexDof(ref)(n);
                int n_loc_Es = IndexDofEs(ref)(n);

                if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                  {
                    // 0 = 2 i beta omega tau eps_tilde U_h V1
                    M_tilde.Get(IndexDirichlet(m), inc + inc_V1 + n_loc_Es) += 2.0 * this->omega * Iwp * tau * eps_tilde * Uh.Value(m,j);
                    // 0 = 2 i beta omega k^2 tau^2 U_h W2
                    M_tilde.Get(IndexDirichlet(m), inc + inc_W2 + n_loc) += 2.0 * this->omega * Iwp * k2 * tau * tau * Uh.Value(m,j);

                    M_tilde.Get(nb_dof_Es + m, inc + inc_V2 + n_loc_Es) += -2.0 * this->omega * Iwp * k2 * tau * tau * Uh.Value(m,j);
                  }
                M_tilde.Get(nb_dof_Es + m, inc + inc_W1 + n_loc) += 2.0 * this->omega * Iwp * tau * mu_tilde * Uh.Value(m,j);
              }

            size_row = Th.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Th.Index(m, j);
                int n_loc = IndexDof(ref)(n);
                int n_loc_Es = IndexDofEs(ref)(n);

                if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                  {
                    // 2 k^2 tau^2 eps_tilde T_h V2 = 0
                    K_tilde.Get(IndexDirichlet(m), inc + inc_V2 + n_loc_Es) += 2.0 * k2 * tau * tau * eps_tilde * Th.Value(m,j);
                  }
                K_tilde.Get(nb_dof_Es + m, inc + inc_W2 + n_loc) += 2.0 * k2 * tau * tau * mu_tilde * Th.Value(m,j);
              }

            size_row = Dh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Dh.Index(m, j);
                int n_loc = IndexDof(ref)(n);
                int n_loc_Es = IndexDofEs(ref)(n);

                if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                  {
                    // eps_tilde tau^2 D_h V1 = 0
                    K_tilde.Get(IndexDirichlet(m), inc + inc_V1 + n_loc_Es) += tau * tau * eps_tilde * Dh.Value(m,j);
                  }
                K_tilde.Get(nb_dof_Es + m, inc + inc_W1 + n_loc) += tau * tau * mu_tilde * Dh.Value(m,j);
              }
          }


        inc += 4 * nb_dof_coucheEs(ref) + 4 * nb_dof_couche(ref); //on arrive au premier indice de la couche suivante dans la matrice
        // inc += 2 * nb_dof_coucheEs(ref) + 2 * nb_dof_couche(ref);
        //DISP(inc);
      }
}
