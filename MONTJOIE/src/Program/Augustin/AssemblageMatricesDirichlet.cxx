
void ModeEz_Solver::AssemblageDirichletMatricesOnlyEz(Matrix<Complex_wp, General, ArrayRowSparse>& K_tilde,
                        Matrix<Complex_wp, General, ArrayRowSparse>& M_tilde,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh)
{
  int N = this->mesh_num.GetNbDof();
  int nb_dof_Ez = DofKeptNonDir.GetM();
  // construction des numerotations
  cout << "N = " << N << endl;
  int nb_couches = vec_Mh.GetM();
  Vector<int> nb_dof_coucheEz(nb_couches);
  nb_dof_coucheEz.Zero();
  Vector<Vector<int> > IndexDofEz(nb_couches), NumGlobEz(nb_couches);
  for (int ref = 0; ref < nb_couches; ref++)
    if (vec_Mh(ref).GetM() == N)
      {
        IndexDofEz(ref).Reallocate(N);
        IndexDofEz(ref).Fill(-1);
      }

  // IndexDof(ref)(j) renvoie le numero local du ddl global j dans la couche ref
  // NumGlob(ref)(i) renvoie le numero global j du ddl local i dans la couche ref
  // IndexDofEz, NumGlobEz pour l'inconnue Ez (ddl Dirichlets enleves)
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    {
      int ref = this->mesh.Element(i).GetReference();
      int nb_dof_loc = this->mesh_num.Element(i).GetNbDof();
      for (int j = 0; j < nb_dof_loc; j++)
        {
          int num_dof = this->mesh_num.Element(i).GetNumberDof(j);
          if ((IndexNonDirichlet(num_dof) >= 0) && (IndexDofEz(ref)(num_dof) == -1)) // si c'est un nouveau
            {
              IndexDofEz(ref)(num_dof) = nb_dof_coucheEz(ref);
              NumGlobEz(ref).PushBack(num_dof);
              nb_dof_coucheEz(ref)++;
            }
        }
    }

  int ordreMatrices = nb_dof_Ez; // lignes pour Ez
  for (int ref = 0 ; ref < nb_couches ; ref++)
    ordreMatrices += nb_dof_coucheEz(ref); // lignes pour V_ref

  K_tilde.Reallocate(ordreMatrices, ordreMatrices);
  M_tilde.Reallocate(ordreMatrices, ordreMatrices);
  DISP(ordreMatrices);
  K_tilde.Zero(); M_tilde.Zero();

  cout << "Création des matrices" << endl;
  int inc = nb_dof_Ez; // premier indice de la première couche dans la matrice
  for (int ref = 0 ; ref < nb_couches ; ref++)
    if (nb_dof_coucheEz(ref) > 0)
    {
      Complex_wp epsilon = this->ref_epsilon(ref);
      Complex_wp mu = this->ref_mu(ref);
      Complex_wp sigma = this->ref_sigma(ref);
      // conversion en non-symetrique
      Matrix<Complexe, Prop, Storage> Kh, Mh;
      // Matrix<Real_wp, General, ArrayRowSparse> Kh, Mh;
      Copy(vec_Kh(ref), Kh);
      Copy(vec_Mh(ref), Mh);
      for (int ddl_loc = 0 ; ddl_loc < nb_dof_coucheEz(ref) ; ddl_loc++)
        {
          int m = NumGlobEz(ref)(ddl_loc);
          // partie -Ez dans l'equation en V_k
          K_tilde.Get(inc + ddl_loc, IndexNonDirichlet(m)) -= 1.0;

          // partie (-omega^2 epsilon - i omega sigma) mu V_k dans l'equation en V_k
          K_tilde.Get(inc + ddl_loc, inc + ddl_loc) -= this->omega * this->omega * epsilon * mu + Iwp * this->omega * sigma * mu;

          // partie -beta^2 omega^2 V_k dans l'equation en Vk
          M_tilde.Get(inc + ddl_loc, inc + ddl_loc) -= this->omega * this->omega;
          // cout << K_tilde.Get(inc + ddl_loc, m) << " ," << K_tilde.Get(inc + ddl_loc + 1, m) << endl;

          // équation en Ez
          int size_row = Kh.GetRowSize(m);
          for (int j = 0 ; j < size_row ; j++)
            {
              int n = Kh.Index(m, j);
              int n_loc = IndexDofEz(ref)(n);
              if (IndexNonDirichlet(n) != -1)
                K_tilde.Get(IndexNonDirichlet(m), n_loc + inc) += (- Iwp * this->omega * epsilon + sigma) * Kh.Value(m, j);
              //if (m == 7)
              // { DISP(ref); DISP(vec_Kh(ref).Index(m, j)); DISP(n_loc+inc); DISP( K_tilde.Get(m, n_loc + inc)); }
            }

          size_row = Mh.GetRowSize(m);
          for (int j = 0 ; j < size_row ; j++)
            {
              int n = Mh.Index(m, j);
              if (IndexNonDirichlet(n) != -1)
                K_tilde.Get(IndexNonDirichlet(m), IndexNonDirichlet(n)) += (- Iwp * this->omega * epsilon + sigma) * Mh.Value(m, j);
            }
        }

      inc += nb_dof_coucheEz(ref); //on arrive au premier indice de la couche suivante dans la matrice
      //DISP(inc);
    }
}

void ModeEz_Solver::AssemblageDirichletMatricesEzHz(Matrix<Complex_wp, General, ArrayRowSparse>& K_tilde,
                        Matrix<Complex_wp, General, ArrayRowSparse>& M_tilde,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch)
{
  int N = this->mesh_num.GetNbDof();
  int nb_dof_Ez = DofKeptNonDir.GetM();
  int nb_dof_Hz = N;
  // construction des numerotations
  cout << "N = " << N << endl;
  int nb_couches = vec_Mh.GetM();
  Vector<int> nb_dof_couche(nb_couches), nb_dof_coucheEz(nb_couches);
  nb_dof_couche.Zero(); nb_dof_coucheEz.Zero();
  Vector<Vector<int> > IndexDofEz(nb_couches), NumGlobEz(nb_couches);
  Vector<Vector<int> > IndexDof(nb_couches), NumGlob(nb_couches);
  for (int ref = 0; ref < nb_couches; ref++)
    if (vec_Mh(ref).GetM() == N)
      {
        IndexDof(ref).Reallocate(N);
        IndexDof(ref).Fill(-1);

        IndexDofEz(ref).Reallocate(N);
        IndexDofEz(ref).Fill(-1);
      }

  // IndexDof(ref)(j) renvoie le numero local du ddl global j dans la couche ref
  // NumGlob(ref)(i) renvoie le numero global j du ddl local i dans la couche ref
  // IndexDofEz, NumGlobEz pour l'inconnue Ez (ddl Dirichlets enleves)
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

          if ((IndexNonDirichlet(num_dof) >= 0) && (IndexDofEz(ref)(num_dof) == -1)) // si c'est un nouveau
            {
              IndexDofEz(ref)(num_dof) = nb_dof_coucheEz(ref);
              NumGlobEz(ref).PushBack(num_dof);
              nb_dof_coucheEz(ref)++;
            }
        }
    }

  // Initialisation
  int ordreMatrices = nb_dof_Ez + nb_dof_Hz; // lignes pour Ez et Hz
  /* lignes pour V_ref et W_ref et leurs variables temporaires
  (pour outrepasser le caractère quadratique de l'équation) */
  for (int ref = 0 ; ref < nb_couches ; ref++)
    ordreMatrices += 2 * nb_dof_couche(ref) + 2 * nb_dof_coucheEz(ref);


  K_tilde.Reallocate(ordreMatrices, ordreMatrices);
  M_tilde.Reallocate(ordreMatrices, ordreMatrices);
  DISP(ordreMatrices);
  K_tilde.Zero(); M_tilde.Zero();

  cout << "Création des matrices" << endl;
  int inc = nb_dof_Ez + nb_dof_Hz; // premier indice de la première couche dans la matrice
  for (int ref = 0 ; ref < nb_couches ; ref++)
    if (nb_dof_couche(ref) > 0)
      {
        Complex_wp epsilon = this->ref_epsilon(ref);
        Complex_wp mu = this->ref_mu(ref);
        Complex_wp sigma = this->ref_sigma(ref);
        // conversion en non-symetrique
        Matrix<Complexe, Prop, Storage> Kh, Mh, Ch;
        // Matrix<Real_wp, General, ArrayRowSparse> Kh, Mh;
        Copy(vec_Kh(ref), Kh);
        Copy(vec_Mh(ref), Mh);
        Copy(vec_Ch(ref), Ch);

        int inc_tempVk = 0;
        int inc_Vk = nb_dof_coucheEz(ref);
        int inc_tempWk = 2 * nb_dof_coucheEz(ref);
        int inc_Wk = 2 * nb_dof_coucheEz(ref) + nb_dof_couche(ref);

        for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
          {
            int m = NumGlob(ref)(ddl_loc); // numéro global avec les noeuds de Dirichlet
            int ddl_loc_Ez = IndexDofEz(ref)(m);
            if (IndexNonDirichlet(m) != -1) // si on n'est pas sur un noeud de Dirichlet
              {
                // équation avec la variable temporaire : TV_k = beta omega V_k
                K_tilde.Get(inc + inc_tempVk + ddl_loc_Ez, inc + inc_tempVk + ddl_loc_Ez) = 1.0;
                M_tilde.Get(inc + inc_tempVk + ddl_loc_Ez, inc + inc_Vk + ddl_loc_Ez) = this->omega;

                // equation en Vk : - Delta_k V_k + E_z = 0
                //                   (omega^2 epsilon mu + i omega sigma mu) V_k + Ez = beta^2 omega^2 Vk = beta omega TV_k
                // partie Ez dans l'equation en V_k
                K_tilde.Get(inc + inc_Vk + ddl_loc_Ez, IndexNonDirichlet(m)) = 1.0;
                // partie (omega^2 epsilon mu + i omega sigma) mu V_k dans l'equation en V_k
                K_tilde.Get(inc + inc_Vk + ddl_loc_Ez, inc + inc_Vk + ddl_loc_Ez) = this->omega * this->omega * epsilon * mu + Iwp * this->omega * sigma * mu;
                // partie beta omega tempV_k dans l'equation en Vk
                M_tilde.Get(inc + inc_Vk + ddl_loc_Ez, inc + inc_tempVk + ddl_loc_Ez) = this->omega;
              }

            // equation en TWk : TW_k = beta omega W_k
            // équation avec la variable temporaire
            K_tilde.Get(inc + inc_tempWk + ddl_loc, inc + inc_tempWk + ddl_loc) = 1.0;
            M_tilde.Get(inc + inc_tempWk + ddl_loc, inc + inc_Wk + ddl_loc) = this->omega;

            // equation en Wk : - Delta_k W_k + H_z = 0
            //                   (omega^2 epsilon mu + i omega sigma mu) W_k + Ez = beta^2 omega^2 Wk = beta omega TW_k
            // partie Hz dans l'equation en W_k
            K_tilde.Get(inc + inc_Wk + ddl_loc, nb_dof_Ez + m) = 1.0;
            // partie (omega^2 epsilon mu + i omega sigma) mu W_k dans l'equation en W_k
            K_tilde.Get(inc + inc_Wk + ddl_loc, inc + inc_Wk + ddl_loc) = this->omega * this->omega * epsilon * mu + Iwp * this->omega * sigma * mu;
            // partie beta omega tempW_k dans l'equation en Wk
            M_tilde.Get(inc + inc_Wk + ddl_loc, inc + inc_tempWk + ddl_loc) = this->omega;


            // équations en Ez et Hz
            // (- i omega epsilon + sigma) M_h E_z + (-i omega epsilon + sigma) K_h V_k = i beta omega C_h W_k
            int size_row = Kh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Kh.Index(m, j);
                int n_loc = IndexDof(ref)(n);
                int n_loc_Ez = IndexDofEz(ref)(n);

                if ((IndexNonDirichlet(m) != -1) && (IndexNonDirichlet(n) != -1))
                  K_tilde.Get(IndexNonDirichlet(m), inc + inc_Vk + n_loc_Ez) = (- Iwp * this->omega * epsilon + sigma) * Kh.Value(m, j);
                K_tilde.Get(nb_dof_Ez + m, inc + inc_Wk + n_loc) = (- Iwp * this->omega * mu) * Kh.Value(m, j);
                //if (m == 7)
                // { DISP(ref); DISP(vec_Kh(ref).Index(m, j)); DISP(n_loc+inc); DISP( K_tilde.Get(m, n_loc + inc)); }
              }

            size_row = Mh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Mh.Index(m, j);
                if ((IndexNonDirichlet(m) != -1) && (IndexNonDirichlet(n) != -1))
                  K_tilde.Get(IndexNonDirichlet(m), IndexNonDirichlet(n)) += (- Iwp * this->omega * epsilon + sigma) * Mh.Value(m, j);
                K_tilde.Get(nb_dof_Ez + m, nb_dof_Ez + n) += (- Iwp * this->omega * mu) * Mh.Value(m, j);
              }

            size_row = Ch.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Ch.Index(m, j);
                int n_loc = IndexDof(ref)(n);
                int n_loc_Ez = IndexDofEz(ref)(n);

                if ((IndexNonDirichlet(m) != -1) && (IndexNonDirichlet(n) != -1))
                  {
                    M_tilde.Get(IndexNonDirichlet(m), n_loc + inc_Wk + inc) = - Iwp * this->omega * Ch.Value(m, j);
                    M_tilde.Get(nb_dof_Ez + m, n_loc_Ez + inc_Vk + inc) = Iwp * this->omega * Ch.Value(m, j);
                  }
              }
          }


        inc += 2 * nb_dof_coucheEz(ref) + 2 * nb_dof_couche(ref); //on arrive au premier indice de la couche suivante dans la matrice
        //DISP(inc);
      }

}
