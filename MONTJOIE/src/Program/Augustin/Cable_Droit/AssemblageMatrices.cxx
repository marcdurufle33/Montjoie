void ModeEz_Solver::AssemblageMatricesOnlyEz(Matrix<Complex_wp, General, ArrayRowSparse>& K_tilde,
                        Matrix<Complex_wp, General, ArrayRowSparse>& M_tilde,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                        const VectReal_wp& Mh_surf,
                        const Matrix<Complexe, Prop, Storage>& Kh_surf)
{
  int N = this->mesh_num.GetNbDof();
  int nbVarAuxiliairesBords = 0;
  int nb_dof_border = NumGlobBorder.GetM();

  // On cherche le nb de variables auxiliaires
  if ((type_CLA == BGT1 && linearization_step == 2) || (type_CLA == BGT2 && linearization_step == 1))
    nbVarAuxiliairesBords = 1;
  else if (type_CLA == BGT2 && linearization_step == 2)
    nbVarAuxiliairesBords = 2;
  else if (type_CLA == BGT2 && linearization_step == 3)
    nbVarAuxiliairesBords = 4;

  // construction des numerotations
  cout << "N = " << N << endl;
  int nb_couches = vec_Mh.GetM();
  Vector<int> nb_dof_couche(nb_couches);
  nb_dof_couche.Zero();
  Vector<Vector<int> > IndexDof(nb_couches), NumGlob(nb_couches);
  Complex_wp rhoPlus = 0;
  for (int ref = 0; ref < nb_couches; ref++)
    if (vec_Mh(ref).GetM() == N)
      {
        IndexDof(ref).Reallocate(N);
        IndexDof(ref).Fill(-1);

        Complex_wp rho = this->ref_epsilon(ref)*this->ref_mu(ref);
        if (realpart(rho) > realpart(rhoPlus))
            rhoPlus = rho;
      }
  DISP(rhoPlus);
  // rhoPlus = 2.0;

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
        }
    }

  cout << "ordreMatrices = " << N << " + ";
  int ordreMatrices = N; // lignes pour Ez
  for (int ref = 0 ; ref < nb_couches ; ref++)
    ordreMatrices += nb_dof_couche(ref); // lignes pour V_ref
  cout << ordreMatrices - N  << " + " << nbVarAuxiliairesBords * nb_dof_border << " = ";


  ordreMatrices += nbVarAuxiliairesBords * nb_dof_border;
  cout << ordreMatrices << endl;



  K_tilde.Reallocate(ordreMatrices, ordreMatrices);
  M_tilde.Reallocate(ordreMatrices, ordreMatrices);
  K_tilde.Zero(); M_tilde.Zero();

  cout << "Création des matrices" << endl;
  cout << "Partie Omega" << endl;
  int inc = N; // premier indice de la première couche dans la matrice
  int inc_water = 0;
  for (int ref = 0 ; ref < nb_couches ; ref++)
    if (nb_dof_couche(ref) > 0)
      {
        Complex_wp epsilon = this->ref_epsilon(ref);
        Complex_wp mu = this->ref_mu(ref);
        Complex_wp sigma = this->ref_sigma(ref);
        // conversion en non-symetrique
        Matrix<Complexe, Prop, Storage> Kh, Mh;
        // Matrix<Real_wp, General, ArrayRowSparse> Kh, Mh;
        Copy(vec_Kh(ref), Kh);
        Copy(vec_Mh(ref), Mh);
        for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
          {
            int m = NumGlob(ref)(ddl_loc);

            // équation en Ez
            int size_row = Kh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Kh.Index(m, j);
                int n_loc = IndexDof(ref)(n);
                K_tilde.Get(m, inc + n_loc) += (- Iwp * this->omega * epsilon + sigma) * Kh.Value(m, j);
              }

            size_row = Mh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Mh.Index(m, j);
                K_tilde.Get(m, n) += (- Iwp * this->omega * epsilon + sigma) * Mh.Value(m, j);
              }

            // partie Ez dans l'equation en V_k
            K_tilde.Get(inc + ddl_loc, m) = 1.0;
            // partie (-omega^2 epsilon - i omega sigma) mu V_k dans l'equation en V_k
            K_tilde.Get(inc + ddl_loc, inc + ddl_loc) = this->omega * this->omega * epsilon * mu
                                                        + Iwp * this->omega * sigma * mu;
            // partie -beta^2 omega^2 V_k dans l'equation en Vk
            M_tilde.Get(inc + ddl_loc, inc + ddl_loc) = this->omega * this->omega;
          }
        if (ref == this->ref_water)
          inc_water = inc;

        inc += nb_dof_couche(ref); //on arrive au premier indice de la couche suivante dans la matrice
      }

  cout << "Partie D Omega" << endl;
  // partie sur le bord
  if (this->type_CLA != NEUMANN)
    {
      DISP(rhoInfini);
      Complex_wp alpha = this->omega * Sqrt(rhoPlus - rhoInfini);
      DISP(alpha);
      Real_wp R = this->ray;
      Complex_wp eps_tilde = - Iwp * this->omega * this->ref_epsilon(this->ref_water) + this->ref_sigma(this->ref_water);
      for (int ddl_loc_border = 0 ; ddl_loc_border < nb_dof_border ; ddl_loc_border++)
        {
          int m = NumGlobBorder(ddl_loc_border); // m est le numéro global
          int loc_Vl = IndexDof(this->ref_water)(m);

          if (this->type_CLA == BGT1)
            {
              if (this->linearization_step == 0)
                {
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * (alpha + 1/(2*R)) * Mh_surf(m);
                }
              else if (this->linearization_step == 1)
                {
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * (alpha/2.0 - rhoInfini * this->omega * this->omega/(2.0*alpha) + 1/(2*R)) * Mh_surf(m);
                  M_tilde.Get(m, inc_water + loc_Vl) -= eps_tilde * Mh_surf(m) / (2.0 * alpha) * this->omega * this->omega;
                }
              else if (this->linearization_step == 2)
                {
                  int mV = inc + ddl_loc_border;
                  // partie heut gauche (eq de u en fonction de u)
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * (alpha/4.0 - rhoInfini * this->omega * this->omega/(4.0*alpha) + 1/(2*R)) * Mh_surf(m);
                  M_tilde.Get(m, inc_water + loc_Vl) -= eps_tilde * Mh_surf(m) / (4.0 * alpha) * this->omega * this->omega;
                  // partie haut droite (eq de u en fonction de v)
                  K_tilde.Get(m, mV) -= eps_tilde * rhoInfini * this->omega * this->omega * Mh_surf(m);
                  M_tilde.Get(m, mV) -= eps_tilde * Mh_surf(m) * this->omega * this->omega;
                  // partie bas gauche (eq de v en fonction de u)
                  K_tilde.Get(mV, inc_water + loc_Vl) += 1.0;
                  // partie bas droite (eq de v en fonction de v)
                  K_tilde.Get(mV, mV) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  M_tilde.Get(mV, mV) += 1.0 / alpha * this->omega * this->omega;
                }
            }
          else if (this->type_CLA == BGT2)
            {
              if (this->linearization_step == 0)
                {
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * ((3/(2*R) - (9/(8*R) + R * rhoInfini * this->omega * this->omega)/(1.0 + R * alpha)) * Mh_surf(m));
                  M_tilde.Get(m, inc_water + loc_Vl) -= eps_tilde * R * Mh_surf(m) / (1.0 + R * alpha) * this->omega * this->omega;

                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * R * Kh_surf.Value(m, j) / (2.0 * (1.0 + R * alpha));
                      //if (m == 7)
                      // { DISP(ref); DISP(vec_Kh(ref).Index(m, j)); DISP(n_loc+inc); DISP( K_tilde.Get(m, n_loc + inc)); }
                    }
                }
              else if (this->linearization_step == 1)
                {
                  int mV = inc + ddl_loc_border;
                  // partie heut gauche (eq de u en fonction de u)
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * 3.0/(2*R) * Mh_surf(m);
                  // partie haut droite (eq de u en fonction de v)
                  K_tilde.Get(m, mV) -= eps_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nV = inc + IndexLocalBorder(n);
                      K_tilde.Get(m, nV) += eps_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(m, mV) -= eps_tilde * R * Mh_surf(m) * this->omega * this->omega;
                  // partie bas gauche (eq de v en fonction de u)
                  K_tilde.Get(mV, inc_water + loc_Vl) += 1.0;
                  // partie bas droite (eq de v en fonction de v)
                  K_tilde.Get(mV, mV) -= 1.0 + alpha * R / 2.0 - R * rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV, mV) += R / (2.0*alpha) * this->omega * this->omega;
                }
              else if (this->linearization_step == 2)
                {
                  int mV = inc + ddl_loc_border;
                  int mW = inc + nb_dof_border + ddl_loc_border;
                  // partie heut gauche (eq de u en fonction de u)
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * 3.0/(2*R) * Mh_surf(m);
                  // partie haut milieu (eq de u en fonction de v)
                  K_tilde.Get(m, mV) -= eps_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nV = inc + IndexLocalBorder(n);
                      K_tilde.Get(m, nV) += eps_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(m, mV) -= eps_tilde * R * Mh_surf(m) * this->omega * this->omega;
                  // partie milieu gauche (eq de v en fonction de u)
                  K_tilde.Get(mV, inc_water + loc_Vl) += 1.0;
                  // partie milieu milieu (eq de v en fonction de v)
                  K_tilde.Get(mV, mV) -= 1.0 + alpha * R / 4.0 - R * rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV, mV) += R / (4.0*alpha) * this->omega * this->omega;
                  // partie milieu droite (eq de v en fonction de w)
                  K_tilde.Get(mV, mW) += R * rhoInfini * this->omega * this->omega;
                  M_tilde.Get(mV, mW) += R * this->omega * this->omega;
                  // partie bas milieu (eq de w en fonction de v)
                  K_tilde.Get(mW, mV) += 1.0;
                  // partie bas droite (eq de v en fonction de w)
                  K_tilde.Get(mW, mW) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  M_tilde.Get(mW, mW) += 1.0 / alpha * this->omega * this->omega;
                }
              else if (this->linearization_step == 3)
                {
                  int mU1 = inc + ddl_loc_border;
                  int mV = inc + nb_dof_border + ddl_loc_border;
                  int mW = inc + 2 * nb_dof_border + ddl_loc_border;
                  int mX = inc + 3 * nb_dof_border + ddl_loc_border;
                  // partie 1 1 (eq de u en fonction de u)
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * 3.0/(2*R) * Mh_surf(m);
                  // partie 1 2 (eq de u en fonction de u')
                  K_tilde.Get(m, mU1) -= eps_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nU1 = inc + IndexLocalBorder(n);
                      K_tilde.Get(m, nU1) += eps_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(m, mU1) -= eps_tilde * R * Mh_surf(m) * this->omega * this->omega;

                  // partie 2 1 (eq de u' en fonction de u)
                  K_tilde.Get(mU1, inc_water + loc_Vl) += 1.0;
                  // partie 2 2 (eq de u' en fonction de u')
                  K_tilde.Get(mU1, mU1) -= 1.0 + alpha * R / 8.0 - R * rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mU1, mU1) += R / (8.0*alpha) * this->omega * this->omega;
                  // partie 2 3 (eq de u' en fonction de v)
                  K_tilde.Get(mU1, mV) += R * rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mU1, mV) += R / 2.0 * this->omega * this->omega;
                  // partie 2 4 (eq de u' en fonction de v)
                  K_tilde.Get(mU1, mW) += R * rhoInfini * this->omega * this->omega;
                  M_tilde.Get(mU1, mW) += R * this->omega * this->omega;

                  // partie 3 2 (eq de v en fonction de u')
                  K_tilde.Get(mV, mU1) += 1.0;
                  // partie 3 3 (eq de v en fonction de v)
                  K_tilde.Get(mV, mV) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  // M_tilde.Get(mV, mV) += R / alpha * this->omega * this->omega;
                  M_tilde.Get(mV, mV) += 1.0 / alpha * this->omega * this->omega;

                  // partie 4 2 (eq de w en fonction de u')
                  // K_tilde.Get(mW, m) += 1.0;
                  K_tilde.Get(mW, mU1) += 1.0;
                  // partie 4 4 (eq de w en fonction de w)
                  K_tilde.Get(mW, mW) -= alpha / 2.0 - rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW, mW) += 1.0 / (2.0*alpha) * this->omega * this->omega;
                  // partie 4 5 (eq de w en fonction de x)
                  K_tilde.Get(mW, mX) += 2.0 * rhoInfini * this->omega * this->omega;
                  M_tilde.Get(mW, mX) += 2.0 * this->omega * this->omega;

                  // partie 5 4 (eq de x en fonction de w)
                  K_tilde.Get(mX, mW) += 1.0;
                  // partie 5 5 (eq de x en fonction de x)
                  K_tilde.Get(mX, mX) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  M_tilde.Get(mX, mX) += 1.0 / alpha * this->omega * this->omega;
                }
            }
        }
    }
}

void ModeEz_Solver::AssemblageMatricesEzHz(Matrix<Complex_wp, General, ArrayRowSparse>& K_tilde,
                        Matrix<Complex_wp, General, ArrayRowSparse>& M_tilde,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                        const Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch,
                        const Matrix<Complexe, Prop, Storage>& Rh_surf,
                        const VectReal_wp& Mh_surf,
                        const Matrix<Complexe, Prop, Storage>& Kh_surf)
{
  int N = this->mesh_num.GetNbDof();
  int nb_dof = N;
  int nbVarAuxiliairesBords = 0;
  int nb_dof_border = NumGlobBorder.GetM();

  // On cherche le nb de variables auxiliaires
  if ((type_CLA == BGT1 && linearization_step == 2) || (type_CLA == BGT2 && linearization_step == 1))
    nbVarAuxiliairesBords = 1;
  else if (type_CLA == BGT2 && linearization_step == 2)
    nbVarAuxiliairesBords = 2;
  else if (type_CLA == BGT2 && linearization_step == 3)
    nbVarAuxiliairesBords = 4;
  else if (type_CLA == BGT2 && linearization_step == 4)
    nbVarAuxiliairesBords = 8;
  else if (type_CLA == BGT2 && linearization_step == 5)
    nbVarAuxiliairesBords = 16;
  else if (type_CLA == BGT2 && linearization_step == 6)
    nbVarAuxiliairesBords = 32;

  // construction des numerotations
  cout << "N = " << N << endl;
  int nb_couches = vec_Mh.GetM();
  Vector<int> nb_dof_couche(nb_couches);
  nb_dof_couche.Zero();
  Vector<Vector<int> > IndexDof(nb_couches), NumGlob(nb_couches);
  Complex_wp rhoPlus = 0;
  for (int ref = 0; ref < nb_couches; ref++)
    if (vec_Mh(ref).GetM() == N)
      {
        IndexDof(ref).Reallocate(N);
        IndexDof(ref).Fill(-1);

        Complex_wp rho = this->ref_epsilon(ref)*this->ref_mu(ref);
        if (realpart(rho) > realpart(rhoPlus))
            rhoPlus = this->ref_epsilon(ref)*this->ref_mu(ref) + Iwp * this->ref_sigma(ref) * this->ref_mu(ref) / this->omega;;
      }
  DISP(rhoPlus);

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
        }
    }


  cout << "ordreMatrices = " << 2*N << " + ";
  int ordreMatrices = 2*N; // lignes pour Ez et Hz
  /* lignes pour V_ref et W_ref et leurs variables temporaires
  (pour outrepasser le caractère quadratique de l'équation) */
  for (int ref = 0 ; ref < nb_couches ; ref++)
    ordreMatrices += 4 * nb_dof_couche(ref);
  cout << ordreMatrices - N  << " + " << 2 * nbVarAuxiliairesBords * nb_dof_border << " = ";

  ordreMatrices += 2 * nbVarAuxiliairesBords * nb_dof_border;
  cout << ordreMatrices << endl;


  K_tilde.Reallocate(ordreMatrices, ordreMatrices);
  M_tilde.Reallocate(ordreMatrices, ordreMatrices);
  K_tilde.Zero(); M_tilde.Zero();

  cout << "Création des matrices" << endl;
  cout << "Partie Omega" << endl;
  int inc = 2*N; // premier indice de la première couche dans la matrice
  int inc_water = 0;
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
        int inc_Vk = nb_dof_couche(ref);
        int inc_tempWk = 2 * nb_dof_couche(ref);
        int inc_Wk = 3 * nb_dof_couche(ref);

        for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
          {
            int m = NumGlob(ref)(ddl_loc);

            // équation avec la variable temporaire : TV_k = beta omega V_k
            K_tilde.Get(inc + inc_tempVk + ddl_loc, inc + inc_tempVk + ddl_loc) = 1.0;
            M_tilde.Get(inc + inc_tempVk + ddl_loc, inc + inc_Vk + ddl_loc) = this->omega;

            // equation en Vk : - Delta_k V_k + E_z = 0
            //                   (omega^2 epsilon mu + i omega sigma mu) V_k + Ez = beta^2 omega^2 Vk = beta omega TV_k
            // partie Ez dans l'equation en V_k
            K_tilde.Get(inc + inc_Vk + ddl_loc, m) = 1.0;
            // partie (omega^2 epsilon mu + i omega sigma) mu V_k dans l'equation en V_k
            K_tilde.Get(inc + inc_Vk + ddl_loc, inc + inc_Vk + ddl_loc) = this->omega * this->omega * epsilon * mu + Iwp * this->omega * sigma * mu;
            // partie beta omega tempV_k dans l'equation en Vk
            M_tilde.Get(inc + inc_Vk + ddl_loc, inc + inc_tempVk + ddl_loc) = this->omega;

            // equation en TWk : TW_k = beta omega W_k
            // équation avec la variable temporaire
            K_tilde.Get(inc + inc_tempWk + ddl_loc, inc + inc_tempWk + ddl_loc) = 1.0;
            M_tilde.Get(inc + inc_tempWk + ddl_loc, inc + inc_Wk + ddl_loc) = this->omega;

            // equation en Wk : - Delta_k W_k + H_z = 0
            //                   (omega^2 epsilon mu + i omega sigma mu) W_k + Ez = beta^2 omega^2 Wk = beta omega TW_k
            // partie Hz dans l'equation en W_k
            K_tilde.Get(inc + inc_Wk + ddl_loc, nb_dof + m) = 1.0;
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

                K_tilde.Get(m, inc + inc_Vk + n_loc) = (- Iwp * this->omega * epsilon + sigma) * Kh.Value(m, j);
                K_tilde.Get(nb_dof + m, inc + inc_Wk + n_loc) = (- Iwp * this->omega * mu) * Kh.Value(m, j);
              }

            size_row = Mh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Mh.Index(m, j);
                K_tilde.Get(m, n) += (- Iwp * this->omega * epsilon + sigma) * Mh.Value(m, j);
                K_tilde.Get(nb_dof + m, nb_dof + n) += (- Iwp * this->omega * mu) * Mh.Value(m, j);
              }

            size_row = Ch.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Ch.Index(m, j);
                int n_loc = IndexDof(ref)(n);

                M_tilde.Get(m, n_loc + inc_Wk + inc) = - Iwp * this->omega * Ch.Value(m, j);
                M_tilde.Get(nb_dof + m, n_loc + inc_Vk + inc) = Iwp * this->omega * Ch.Value(m, j);
              }
          }
        if (ref == this->ref_water)
          inc_water = inc;

        inc += 4 * nb_dof_couche(ref); //on arrive au premier indice de la couche suivante dans la matrice
        //DISP(inc);
      }

  cout << "Partie D Omega" << endl;
  Real_wp R = this->ray;
  for (int ddl_loc_border = 0 ; ddl_loc_border < nb_dof_border ; ddl_loc_border++)
    {
      int m = NumGlobBorder(ddl_loc_border); // m est le numéro global
      int size_row = Rh_surf.GetRowSize(m);
      for (int j = 0 ; j < size_row ; j++)
        {
          int n = Rh_surf.Index(m, j);
          int nW = inc_water + 3 * nb_dof_couche(this->ref_water) + IndexLocalBorder(n);
          int nV = inc_water + nb_dof_couche(this->ref_water) + IndexLocalBorder(n);
          M_tilde.Get(m, nW) += - Iwp * this->omega * Rh_surf.Value(m, j);
          M_tilde.Get(nb_dof + m, nV) += Iwp * this->omega * Rh_surf.Value(m, j);
        }
    }

  if (this->type_CLA != NEUMANN)
    {
      DISP(rhoInfini);
      rhoPlus = 2.1198386293541 + 0.212979996978107 * Iwp;
      Complex_wp alpha = this->omega * Sqrt(rhoPlus - rhoInfini);
      // Complex_wp alpha = this->omega * Sqrt(rhoPlus - real(rhoInfini)) + Iwp * 0.1;
      DISP(alpha);
      Real_wp R = this->ray;
      Complex_wp eps_tilde = - Iwp * this->omega * this->ref_epsilon(this->ref_water) + this->ref_sigma(this->ref_water);
      Complex_wp mu_tilde = - Iwp * this->omega * this->ref_mu(this->ref_water);
      for (int ddl_loc_border = 0 ; ddl_loc_border < nb_dof_border ; ddl_loc_border++)
        {
          int m = NumGlobBorder(ddl_loc_border); // m est le numéro global
          int loc_Vl = nb_dof_couche(this->ref_water) + IndexDof(this->ref_water)(m);
          int loc_Wl = 3 * nb_dof_couche(this->ref_water) + IndexDof(this->ref_water)(m);

          if (this->type_CLA == BGT1)
            {
              if (this->linearization_step == 0)
                {
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * (alpha + 1/(2*R)) * Mh_surf(m);

                  K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * (alpha + 1/(2*R)) * Mh_surf(m);
                }
              else if (this->linearization_step == 1)
                {
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * (alpha/2.0 - rhoInfini * this->omega * this->omega/(2.0*alpha) + 1/(2*R)) * Mh_surf(m);
                  M_tilde.Get(m, inc_water + loc_Vl) -= eps_tilde * Mh_surf(m) / (2.0 * alpha) * this->omega * this->omega;

                  K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * (alpha/2.0 - rhoInfini * this->omega * this->omega/(2.0*alpha) + 1/(2*R)) * Mh_surf(m);
                  M_tilde.Get(nb_dof + m, inc_water + loc_Wl) -= mu_tilde * Mh_surf(m) / (2.0 * alpha) * this->omega * this->omega;
                }
              else if (this->linearization_step == 2)
                {
                  int mV1 = inc + ddl_loc_border;
                  int mW1 = inc + nb_dof_border + ddl_loc_border;
                  // partie heut gauche (eq de u en fonction de u)
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * (alpha/4.0 - rhoInfini * this->omega * this->omega/(4.0*alpha) + 1/(2*R)) * Mh_surf(m);
                  M_tilde.Get(m, inc_water + loc_Vl) -= eps_tilde * Mh_surf(m) / (4.0 * alpha) * this->omega * this->omega;
                  // partie haut droite (eq de u en fonction de v)
                  K_tilde.Get(m, mV1) -= eps_tilde * rhoInfini * this->omega * this->omega * Mh_surf(m);
                  M_tilde.Get(m, mV1) -= eps_tilde * Mh_surf(m) * this->omega * this->omega;
                  // partie bas gauche (eq de v en fonction de u)
                  K_tilde.Get(mV1, inc_water + loc_Vl) += 1.0;
                  // partie bas droite (eq de v en fonction de v)
                  K_tilde.Get(mV1, mV1) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  M_tilde.Get(mV1, mV1) += 1.0 / alpha * this->omega * this->omega;

                  // partie heut gauche (eq de u en fonction de u)
                  K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * (alpha/4.0 - rhoInfini * this->omega * this->omega/(4.0*alpha) + 1/(2*R)) * Mh_surf(m);
                  M_tilde.Get(nb_dof + m, inc_water + loc_Wl) -= mu_tilde * Mh_surf(m) / (4.0 * alpha) * this->omega * this->omega;
                  // partie haut droite (eq de u en fonction de v)
                  K_tilde.Get(nb_dof + m, mW1) -= mu_tilde * rhoInfini * this->omega * this->omega * Mh_surf(m);
                  M_tilde.Get(nb_dof + m, mW1) -= mu_tilde * Mh_surf(m) * this->omega * this->omega;
                  // partie bas gauche (eq de v en fonction de u)
                  K_tilde.Get(mW1, inc_water + loc_Wl) += 1.0;
                  // partie bas droite (eq de v en fonction de v)
                  K_tilde.Get(mW1, mW1) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  M_tilde.Get(mW1, mW1) += 1.0 / alpha * this->omega * this->omega;
                }
            }
          else if (this->type_CLA == BGT2)
            {
              if (this->linearization_step == 0)
                {
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * ((3/(2*R) - (9/(8*R) + R * rhoInfini * this->omega * this->omega)/(1.0 + R * alpha)) * Mh_surf(m));
                  M_tilde.Get(m, inc_water + loc_Vl) -= eps_tilde * R * Mh_surf(m) / (1.0 + R * alpha) * this->omega * this->omega;

                  K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * ((3/(2*R) - (9/(8*R) + R * rhoInfini * this->omega * this->omega)/(1.0 + R * alpha)) * Mh_surf(m));
                  M_tilde.Get(nb_dof + m, inc_water + loc_Wl) -= mu_tilde * R * Mh_surf(m) / (1.0 + R * alpha) * this->omega * this->omega;

                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      // int n = Kh_surf.Index(m, j);
                      K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * R * Kh_surf.Value(m, j) / (2.0 * (1.0 + R * alpha));

                      K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * R * Kh_surf.Value(m, j) / (2.0 * (1.0 + R * alpha));
                      //if (m == 7)
                      // { DISP(ref); DISP(vec_Kh(ref).Index(m, j)); DISP(n_loc+inc); DISP( K_tilde.Get(m, n_loc + inc)); }
                    }
                }
              else if (this->linearization_step == 1)
                {
                  int mV1 = inc + ddl_loc_border;
                  int mW1 = inc + nb_dof_border + ddl_loc_border;
                  // partie heut gauche (eq de u en fonction de u)
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * 3.0/(2*R) * Mh_surf(m);
                  // partie haut droite (eq de u en fonction de v)
                  K_tilde.Get(m, mV1) -= eps_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nV = inc + IndexLocalBorder(n);
                      K_tilde.Get(m, nV) += eps_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(m, mV1) -= eps_tilde * R * Mh_surf(m) * this->omega * this->omega;
                  // partie bas gauche (eq de v en fonction de u)
                  K_tilde.Get(mV1, inc_water + loc_Vl) += 1.0;
                  // partie bas droite (eq de v en fonction de v)
                  K_tilde.Get(mV1, mV1) -= 1.0 + alpha * R / 2.0 - R * rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV1, mV1) += R / (2.0*alpha) * this->omega * this->omega;

                  // partie heut gauche (eq de u en fonction de u)
                  K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * 3.0/(2*R) * Mh_surf(m);
                  // partie haut droite (eq de u en fonction de v)
                  K_tilde.Get(nb_dof + m, mW1) -= mu_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nW = inc + nb_dof_border + IndexLocalBorder(n);
                      K_tilde.Get(nb_dof + m, nW) += mu_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(nb_dof + m, mW1) -= mu_tilde * R * Mh_surf(m) * this->omega * this->omega;
                  // partie bas gauche (eq de v en fonction de u)
                  K_tilde.Get(mW1, inc_water + loc_Wl) += 1.0;
                  // partie bas droite (eq de v en fonction de v)
                  K_tilde.Get(mW1, mW1) -= 1.0 + alpha * R / 2.0 - R * rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW1, mW1) += R / (2.0*alpha) * this->omega * this->omega;
                }
              else if (this->linearization_step == 2)
                {

                  int mV1 = inc + ddl_loc_border;
                  int mV2 = inc + nb_dof_border + ddl_loc_border;
                  int mW1 = inc + 2*nb_dof_border + ddl_loc_border;
                  int mW2 = inc + 3*nb_dof_border + ddl_loc_border;

                  // partie heut gauche (eq de u en fonction de u)
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * 3.0/(2*R) * Mh_surf(m);
                  // partie haut milieu (eq de u en fonction de v)
                  K_tilde.Get(m, mV1) -= eps_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nV = inc + IndexLocalBorder(n);
                      K_tilde.Get(m, nV) += eps_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(m, mV1) -= eps_tilde * R * Mh_surf(m) * this->omega * this->omega;

                  // partie milieu gauche (eq de v en fonction de u)
                  K_tilde.Get(mV1, inc_water + loc_Vl) += 1.0;
                  // partie milieu milieu (eq de v en fonction de v)
                  K_tilde.Get(mV1, mV1) -= 1.0 + alpha * R / 4.0 - R * rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV1, mV1) += R / (4.0*alpha) * this->omega * this->omega;
                  // partie milieu droite (eq de v en fonction de w)
                  K_tilde.Get(mV1, mV2) += R * rhoInfini * this->omega * this->omega;
                  M_tilde.Get(mV1, mV2) += R * this->omega * this->omega;
                  // partie bas milieu (eq de w en fonction de v)
                  K_tilde.Get(mV2, mV1) += 1.0;
                  // partie bas droite (eq de v en fonction de w)
                  K_tilde.Get(mV2, mV2) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  M_tilde.Get(mV2, mV2) += 1.0 / alpha * this->omega * this->omega;

                  // partie heut gauche (eq de u en fonction de u)
                  K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * 3.0/(2*R) * Mh_surf(m);
                  // partie haut milieu (eq de u en fonction de v)
                  K_tilde.Get(nb_dof + m, mW1) -= mu_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nW = inc + nb_dof_border + IndexLocalBorder(n);
                      K_tilde.Get(nb_dof + m, nW) += mu_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(nb_dof + m, mW1) -= mu_tilde * R * Mh_surf(m) * this->omega * this->omega;
                  // partie milieu gauche (eq de v en fonction de u)
                  K_tilde.Get(mW1, inc_water + loc_Wl) += 1.0;
                  // partie milieu milieu (eq de v en fonction de v)
                  K_tilde.Get(mW1, mW1) -= 1.0 + alpha * R / 4.0 - R * rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW1, mW1) += R / (4.0*alpha) * this->omega * this->omega;
                  // partie milieu droite (eq de v en fonction de w)
                  K_tilde.Get(mW1, mW2) += R * rhoInfini * this->omega * this->omega;
                  M_tilde.Get(mW1, mW2) += R * this->omega * this->omega;
                  // partie bas milieu (eq de w en fonction de v)
                  K_tilde.Get(mW2, mW1) += 1.0;
                  // partie bas droite (eq de v en fonction de w)
                  K_tilde.Get(mW2, mW2) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  M_tilde.Get(mW2, mW2) += 1.0 / alpha * this->omega * this->omega;
                }
              else if (this->linearization_step == 3)
                {
                  int mV1 = inc + ddl_loc_border;
                  int mV2 = inc + nb_dof_border + ddl_loc_border;
                  int mV3 = inc + 2 * nb_dof_border + ddl_loc_border;
                  int mV4 = inc + 3 * nb_dof_border + ddl_loc_border;
                  int mW1 = inc + 4 * nb_dof_border + ddl_loc_border;
                  int mW2 = inc + 5 * nb_dof_border + ddl_loc_border;
                  int mW3 = inc + 6 * nb_dof_border + ddl_loc_border;
                  int mW4 = inc + 7 * nb_dof_border + ddl_loc_border;
                  // partie 1 1 (eq de u en fonction de u)
                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * 3.0/(2*R) * Mh_surf(m);
                  // partie 1 2 (eq de u en fonction de u')
                  K_tilde.Get(m, mV1) -= eps_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nV = inc + IndexLocalBorder(n);
                      K_tilde.Get(m, nV) += eps_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(m, mV1) -= eps_tilde * R * Mh_surf(m) * this->omega * this->omega;

                  // partie 2 1 (eq de u' en fonction de u)
                  K_tilde.Get(mV1, inc_water + loc_Vl) += 1.0;
                  // partie 2 2 (eq de u' en fonction de u')
                  K_tilde.Get(mV1, mV1) -= 1.0 + alpha * R / 8.0 - R * rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mV1, mV1) += R / (8.0*alpha) * this->omega * this->omega;
                  // partie 2 3 (eq de u' en fonction de v)
                  K_tilde.Get(mV1, mV2) += R * rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV1, mV2) += R / 2.0 * this->omega * this->omega;
                  // partie 2 4 (eq de u' en fonction de v)
                  K_tilde.Get(mV1, mV3) += R * rhoInfini * this->omega * this->omega;
                  M_tilde.Get(mV1, mV3) += R * this->omega * this->omega;

                  // partie 3 2 (eq de v en fonction de u')
                  K_tilde.Get(mV2, mV1) += 1.0;
                  // partie 3 3 (eq de v en fonction de v)
                  K_tilde.Get(mV2, mV2) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  // M_tilde.Get(mV2, mV2) += R / alpha * this->omega * this->omega;
                  M_tilde.Get(mV2, mV2) += 1.0 / alpha * this->omega * this->omega;

                  // partie 4 2 (eq de w en fonction de u')
                  // K_tilde.Get(mV3, m) += 1.0;
                  K_tilde.Get(mV3, mV1) += 1.0;
                  // partie 4 4 (eq de w en fonction de w)
                  K_tilde.Get(mV3, mV3) -= alpha / 2.0 - rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV3, mV3) += 1.0 / (2.0*alpha) * this->omega * this->omega;
                  // partie 4 5 (eq de w en fonction de x)
                  K_tilde.Get(mV3, mV4) += 2.0 * rhoInfini * this->omega * this->omega;
                  M_tilde.Get(mV3, mV4) += 2.0 * this->omega * this->omega;

                  // partie 5 4 (eq de x en fonction de w)
                  K_tilde.Get(mV4, mV3) += 1.0;
                  // partie 5 5 (eq de x en fonction de x)
                  K_tilde.Get(mV4, mV4) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  M_tilde.Get(mV4, mV4) += 1.0 / alpha * this->omega * this->omega;


                  // partie 1 1 (eq de u en fonction de u)
                  K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * 3.0/(2*R) * Mh_surf(m);
                  // partie 1 2 (eq de u en fonction de u')
                  K_tilde.Get(nb_dof + m, mW1) -= mu_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nW = inc + nb_dof_border + IndexLocalBorder(n);
                      K_tilde.Get(nb_dof + m, nW) += mu_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(nb_dof + m, mW1) -= mu_tilde * R * Mh_surf(m) * this->omega * this->omega;

                  // partie 2 1 (eq de u' en fonction de u)
                  K_tilde.Get(mW1, inc_water + loc_Wl) += 1.0;
                  // partie 2 2 (eq de u' en fonction de u')
                  K_tilde.Get(mW1, mW1) -= 1.0 + alpha * R / 8.0 - R * rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mW1, mW1) += R / (8.0*alpha) * this->omega * this->omega;
                  // partie 2 3 (eq de u' en fonction de v)
                  K_tilde.Get(mW1, mW2) += R * rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW1, mW2) += R / 2.0 * this->omega * this->omega;
                  // partie 2 4 (eq de u' en fonction de v)
                  K_tilde.Get(mW1, mW3) += R * rhoInfini * this->omega * this->omega;
                  M_tilde.Get(mW1, mW3) += R * this->omega * this->omega;

                  // partie 3 2 (eq de v en fonction de u')
                  K_tilde.Get(mW2, mW1) += 1.0;
                  // partie 3 3 (eq de v en fonction de v)
                  K_tilde.Get(mW2, mW2) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  // M_tilde.Get(mW2, mW2) += R / alpha * this->omega * this->omega;
                  M_tilde.Get(mW2, mW2) += 1.0 / alpha * this->omega * this->omega;

                  // partie 4 2 (eq de w en fonction de u')
                  // K_tilde.Get(mW3, m) += 1.0;
                  K_tilde.Get(mW3, mW1) += 1.0;
                  // partie 4 4 (eq de w en fonction de w)
                  K_tilde.Get(mW3, mW3) -= alpha / 2.0 - rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW3, mW3) += 1.0 / (2.0*alpha) * this->omega * this->omega;
                  // partie 4 5 (eq de w en fonction de x)
                  K_tilde.Get(mW3, mW4) += 2.0 * rhoInfini * this->omega * this->omega;
                  M_tilde.Get(mW3, mW4) += 2.0 * this->omega * this->omega;

                  // partie 5 4 (eq de x en fonction de w)
                  K_tilde.Get(mW4, mW3) += 1.0;
                  // partie 5 5 (eq de x en fonction de x)
                  K_tilde.Get(mW4, mW4) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                  M_tilde.Get(mW4, mW4) += 1.0 / alpha * this->omega * this->omega;
                }
              else if (this->linearization_step == 4)
                {
                  int mV1 = inc + 0 * nb_dof_border + ddl_loc_border;
                  int mV2 = inc + 1 * nb_dof_border + ddl_loc_border;
                  int mV3 = inc + 2 * nb_dof_border + ddl_loc_border;
                  int mV4 = inc + 3 * nb_dof_border + ddl_loc_border;
                  int mV5 = inc + 4 * nb_dof_border + ddl_loc_border;
                  int mV6 = inc + 5 * nb_dof_border + ddl_loc_border;
                  int mV7 = inc + 6 * nb_dof_border + ddl_loc_border;
                  int mV8 = inc + 7 * nb_dof_border + ddl_loc_border;
                  int mW1 = inc + 8 * nb_dof_border + ddl_loc_border;
                  int mW2 = inc + 9 * nb_dof_border + ddl_loc_border;
                  int mW3 = inc + 10 * nb_dof_border + ddl_loc_border;
                  int mW4 = inc + 11 * nb_dof_border + ddl_loc_border;
                  int mW5 = inc + 12 * nb_dof_border + ddl_loc_border;
                  int mW6 = inc + 13 * nb_dof_border + ddl_loc_border;
                  int mW7 = inc + 14 * nb_dof_border + ddl_loc_border;
                  int mW8 = inc + 15 * nb_dof_border + ddl_loc_border;

                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * 3.0/(2*R) * Mh_surf(m);
                  K_tilde.Get(m, mV1) -= eps_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nV = inc + IndexLocalBorder(n);
                      K_tilde.Get(m, nV) += eps_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(m, mV1) -= eps_tilde * R * Mh_surf(m) * this->omega * this->omega;

                  K_tilde.Get(mV1, inc_water + loc_Vl) += 1.0;
                  K_tilde.Get(mV1, mV1) += -1.0 - alpha * R / 16.0 + R * rhoInfini * this->omega * this->omega / (16.0*alpha);
                  M_tilde.Get(mV1, mV1) += R / (16.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV1, mV2) += R * rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV1, mV2) += R / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV1, mV4) += R * rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV1, mV4) += R / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mV1, mV8) += R * rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mV1, mV8) += R / 8.0 * this->omega * this->omega;

                  K_tilde.Get(mV4, mV1) += 1.0;
                  K_tilde.Get(mV4, mV4) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV4, mV4) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV4, mV5) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV4, mV5) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV2, mV1) += 1.0;
                  K_tilde.Get(mV2, mV2) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mV2, mV2) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV2, mV3) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV2, mV3) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV2, mV7) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV2, mV7) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mV3, mV2) += 1.0;
                  K_tilde.Get(mV3, mV3) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV3, mV3) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV3, mV6) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV3, mV6) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV8, mV1) += 1.0;
                  K_tilde.Get(mV8, mV8) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV8, mV8) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV7, mV2) += 1.0;
                  K_tilde.Get(mV7, mV7) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV7, mV7) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV6, mV3) += 1.0;
                  K_tilde.Get(mV6, mV6) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV6, mV6) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV5, mV4) += 1.0;
                  K_tilde.Get(mV5, mV5) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV5, mV5) += 1.0 / (2.0*alpha) * this->omega * this->omega;


                  K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * 3.0/(2*R) * Mh_surf(m);
                  K_tilde.Get(nb_dof + m, mW1) -= mu_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nW = inc + nb_dof_border + IndexLocalBorder(n);
                      K_tilde.Get(nb_dof + m, nW) += mu_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(nb_dof + m, mW1) -= mu_tilde * R * Mh_surf(m) * this->omega * this->omega;

                  K_tilde.Get(mW1, inc_water + loc_Wl) += 1.0;
                  K_tilde.Get(mW1, mW1) += -1.0 - alpha * R / 16.0 + R * rhoInfini * this->omega * this->omega / (16.0*alpha);
                  M_tilde.Get(mW1, mW1) += R / (16.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW1, mW2) += R * rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW1, mW2) += R / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW1, mW4) += R * rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW1, mW4) += R / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mW1, mW8) += R * rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mW1, mW8) += R / 8.0 * this->omega * this->omega;

                  K_tilde.Get(mW4, mW1) += 1.0;
                  K_tilde.Get(mW4, mW4) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW4, mW4) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW4, mW5) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW4, mW5) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW2, mW1) += 1.0;
                  K_tilde.Get(mW2, mW2) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mW2, mW2) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW2, mW3) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW2, mW3) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW2, mW7) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW2, mW7) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mW3, mW2) += 1.0;
                  K_tilde.Get(mW3, mW3) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW3, mW3) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW3, mW6) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW3, mW6) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW8, mW1) += 1.0;
                  K_tilde.Get(mW8, mW8) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW8, mW8) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW7, mW2) += 1.0;
                  K_tilde.Get(mW7, mW7) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW7, mW7) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW6, mW3) += 1.0;
                  K_tilde.Get(mW6, mW6) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW6, mW6) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW5, mW4) += 1.0;
                  K_tilde.Get(mW5, mW5) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW5, mW5) += 1.0 / (2.0*alpha) * this->omega * this->omega;
                }
              else if (this->linearization_step == 5)
                {
                  int mV1 = inc + 0 * nb_dof_border + ddl_loc_border;
                  int mV2 = inc + 1 * nb_dof_border + ddl_loc_border;
                  int mV3 = inc + 2 * nb_dof_border + ddl_loc_border;
                  int mV4 = inc + 3 * nb_dof_border + ddl_loc_border;
                  int mV5 = inc + 4 * nb_dof_border + ddl_loc_border;
                  int mV6 = inc + 5 * nb_dof_border + ddl_loc_border;
                  int mV7 = inc + 6 * nb_dof_border + ddl_loc_border;
                  int mV8 = inc + 7 * nb_dof_border + ddl_loc_border;
                  int mV9 = inc + 8 * nb_dof_border + ddl_loc_border;
                  int mV10 = inc + 9 * nb_dof_border + ddl_loc_border;
                  int mV11 = inc + 10 * nb_dof_border + ddl_loc_border;
                  int mV12 = inc + 11 * nb_dof_border + ddl_loc_border;
                  int mV13 = inc + 12 * nb_dof_border + ddl_loc_border;
                  int mV14 = inc + 13 * nb_dof_border + ddl_loc_border;
                  int mV15 = inc + 14 * nb_dof_border + ddl_loc_border;
                  int mV16 = inc + 15 * nb_dof_border + ddl_loc_border;
                  int mW1 = inc + 16 * nb_dof_border + ddl_loc_border;
                  int mW2 = inc + 17 * nb_dof_border + ddl_loc_border;
                  int mW3 = inc + 18 * nb_dof_border + ddl_loc_border;
                  int mW4 = inc + 19 * nb_dof_border + ddl_loc_border;
                  int mW5 = inc + 20 * nb_dof_border + ddl_loc_border;
                  int mW6 = inc + 21 * nb_dof_border + ddl_loc_border;
                  int mW7 = inc + 22 * nb_dof_border + ddl_loc_border;
                  int mW8 = inc + 23 * nb_dof_border + ddl_loc_border;
                  int mW9 = inc + 24 * nb_dof_border + ddl_loc_border;
                  int mW10 = inc + 25 * nb_dof_border + ddl_loc_border;
                  int mW11 = inc + 26 * nb_dof_border + ddl_loc_border;
                  int mW12 = inc + 27 * nb_dof_border + ddl_loc_border;
                  int mW13 = inc + 28 * nb_dof_border + ddl_loc_border;
                  int mW14 = inc + 29 * nb_dof_border + ddl_loc_border;
                  int mW15 = inc + 30 * nb_dof_border + ddl_loc_border;
                  int mW16 = inc + 31 * nb_dof_border + ddl_loc_border;

                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * 3.0/(2*R) * Mh_surf(m);
                  K_tilde.Get(m, mV1) -= eps_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nV = inc + IndexLocalBorder(n);
                      K_tilde.Get(m, nV) += eps_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(m, mV1) -= eps_tilde * R * Mh_surf(m) * this->omega * this->omega;

                  K_tilde.Get(mV1, inc_water + loc_Vl) += 1.0;
                  K_tilde.Get(mV1, mV1) += -1.0 - alpha * R / 32.0 + R * rhoInfini * this->omega * this->omega / (32.0*alpha);
                  M_tilde.Get(mV1, mV1) += R / (32.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV1, mV2) += R * rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV1, mV2) += R / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV1, mV4) += R * rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV1, mV4) += R / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mV1, mV8) += R * rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mV1, mV8) += R / 8.0 * this->omega * this->omega;
                  K_tilde.Get(mV1, mV16) += R * rhoInfini * this->omega * this->omega / 16.0;
                  M_tilde.Get(mV1, mV16) += R / 16.0 * this->omega * this->omega;

                  K_tilde.Get(mV8, mV1) += 1.0;
                  K_tilde.Get(mV8, mV8) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV8, mV8) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV8, mV9) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV8, mV9) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV4, mV1) += 1.0;
                  K_tilde.Get(mV4, mV4) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mV4, mV4) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV4, mV5) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV4, mV5) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV4, mV13) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV4, mV13) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mV2, mV1) += 1.0;
                  K_tilde.Get(mV2, mV2) += - alpha  / 16.0 + rhoInfini * this->omega * this->omega / (16.0*alpha);
                  M_tilde.Get(mV2, mV2) += 1.0 / (16.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV2, mV3) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV2, mV3) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV2, mV7) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV2, mV7) += 1.0 / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mV2, mV15) += rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mV2, mV15) += 1.0 / 8.0 * this->omega * this->omega;

                  K_tilde.Get(mV7, mV2) += 1.0;
                  K_tilde.Get(mV7, mV7) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV7, mV7) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV7, mV10) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV7, mV10) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV3, mV2) += 1.0;
                  K_tilde.Get(mV3, mV3) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mV3, mV3) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV3, mV6) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV3, mV6) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV3, mV14) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV3, mV14) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mV6, mV3) += 1.0;
                  K_tilde.Get(mV6, mV6) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV6, mV6) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV6, mV11) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV6, mV11) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV5, mV4) += 1.0;
                  K_tilde.Get(mV5, mV5) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV5, mV5) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV5, mV12) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV5, mV12) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV16, mV1) += 1.0;
                  K_tilde.Get(mV16, mV16) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV16, mV16) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV15, mV2) += 1.0;
                  K_tilde.Get(mV15, mV15) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV15, mV15) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV14, mV3) += 1.0;
                  K_tilde.Get(mV14, mV14) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV14, mV14) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV13, mV4) += 1.0;
                  K_tilde.Get(mV13, mV13) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV13, mV13) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV12, mV5) += 1.0;
                  K_tilde.Get(mV12, mV12) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV12, mV12) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV11, mV6) += 1.0;
                  K_tilde.Get(mV11, mV11) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV11, mV11) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV10, mV7) += 1.0;
                  K_tilde.Get(mV10, mV10) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV10, mV10) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV9, mV8) += 1.0;
                  K_tilde.Get(mV9, mV9) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV9, mV9) += 1.0 / (2.0*alpha) * this->omega * this->omega;


                  K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * 3.0/(2*R) * Mh_surf(m);
                  K_tilde.Get(nb_dof + m, mW1) -= mu_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nW = inc + nb_dof_border + IndexLocalBorder(n);
                      K_tilde.Get(nb_dof + m, nW) += mu_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(nb_dof + m, mW1) -= mu_tilde * R * Mh_surf(m) * this->omega * this->omega;

                  K_tilde.Get(mW1, inc_water + loc_Wl) += 1.0;
                  K_tilde.Get(mW1, mW1) += -1.0 - alpha * R / 32.0 + R * rhoInfini * this->omega * this->omega / (32.0*alpha);
                  M_tilde.Get(mW1, mW1) += R / (32.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW1, mW2) += R * rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW1, mW2) += R / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW1, mW4) += R * rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW1, mW4) += R / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mW1, mW8) += R * rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mW1, mW8) += R / 8.0 * this->omega * this->omega;
                  K_tilde.Get(mW1, mW16) += R * rhoInfini * this->omega * this->omega / 16.0;
                  M_tilde.Get(mW1, mW16) += R / 16.0 * this->omega * this->omega;

                  K_tilde.Get(mW8, mW1) += 1.0;
                  K_tilde.Get(mW8, mW8) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW8, mW8) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW8, mW9) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW8, mW9) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW4, mW1) += 1.0;
                  K_tilde.Get(mW4, mW4) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mW4, mW4) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW4, mW5) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW4, mW5) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW4, mW13) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW4, mW13) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mW2, mW1) += 1.0;
                  K_tilde.Get(mW2, mW2) += - alpha  / 16.0 + rhoInfini * this->omega * this->omega / (16.0*alpha);
                  M_tilde.Get(mW2, mW2) += 1.0 / (16.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW2, mW3) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW2, mW3) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW2, mW7) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW2, mW7) += 1.0 / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mW2, mW15) += rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mW2, mW15) += 1.0 / 8.0 * this->omega * this->omega;

                  K_tilde.Get(mW7, mW2) += 1.0;
                  K_tilde.Get(mW7, mW7) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW7, mW7) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW7, mW10) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW7, mW10) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW3, mW2) += 1.0;
                  K_tilde.Get(mW3, mW3) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mW3, mW3) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW3, mW6) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW3, mW6) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW3, mW14) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW3, mW14) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mW6, mW3) += 1.0;
                  K_tilde.Get(mW6, mW6) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW6, mW6) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW6, mW11) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW6, mW11) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW5, mW4) += 1.0;
                  K_tilde.Get(mW5, mW5) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW5, mW5) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW5, mW12) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW5, mW12) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW16, mW1) += 1.0;
                  K_tilde.Get(mW16, mW16) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW16, mW16) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW15, mW2) += 1.0;
                  K_tilde.Get(mW15, mW15) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW15, mW15) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW14, mW3) += 1.0;
                  K_tilde.Get(mW14, mW14) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW14, mW14) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW13, mW4) += 1.0;
                  K_tilde.Get(mW13, mW13) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW13, mW13) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW12, mW5) += 1.0;
                  K_tilde.Get(mW12, mW12) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW12, mW12) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW11, mW6) += 1.0;
                  K_tilde.Get(mW11, mW11) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW11, mW11) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW10, mW7) += 1.0;
                  K_tilde.Get(mW10, mW10) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW10, mW10) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW9, mW8) += 1.0;
                  K_tilde.Get(mW9, mW9) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW9, mW9) += 1.0 / (2.0*alpha) * this->omega * this->omega;
                }
              else if (this->linearization_step == 6)
                {
                  int mV1 = inc + 0 * nb_dof_border + ddl_loc_border;
                  int mV2 = inc + 1 * nb_dof_border + ddl_loc_border;
                  int mV3 = inc + 2 * nb_dof_border + ddl_loc_border;
                  int mV4 = inc + 3 * nb_dof_border + ddl_loc_border;
                  int mV5 = inc + 4 * nb_dof_border + ddl_loc_border;
                  int mV6 = inc + 5 * nb_dof_border + ddl_loc_border;
                  int mV7 = inc + 6 * nb_dof_border + ddl_loc_border;
                  int mV8 = inc + 7 * nb_dof_border + ddl_loc_border;
                  int mV9 = inc + 8 * nb_dof_border + ddl_loc_border;
                  int mV10 = inc + 9 * nb_dof_border + ddl_loc_border;
                  int mV11 = inc + 10 * nb_dof_border + ddl_loc_border;
                  int mV12 = inc + 11 * nb_dof_border + ddl_loc_border;
                  int mV13 = inc + 12 * nb_dof_border + ddl_loc_border;
                  int mV14 = inc + 13 * nb_dof_border + ddl_loc_border;
                  int mV15 = inc + 14 * nb_dof_border + ddl_loc_border;
                  int mV16 = inc + 15 * nb_dof_border + ddl_loc_border;
                  int mV17 = inc + 16 * nb_dof_border + ddl_loc_border;
                  int mV18 = inc + 17 * nb_dof_border + ddl_loc_border;
                  int mV19 = inc + 18 * nb_dof_border + ddl_loc_border;
                  int mV20 = inc + 19 * nb_dof_border + ddl_loc_border;
                  int mV21 = inc + 20 * nb_dof_border + ddl_loc_border;
                  int mV22 = inc + 21 * nb_dof_border + ddl_loc_border;
                  int mV23 = inc + 22 * nb_dof_border + ddl_loc_border;
                  int mV24 = inc + 23 * nb_dof_border + ddl_loc_border;
                  int mV25 = inc + 24 * nb_dof_border + ddl_loc_border;
                  int mV26 = inc + 25 * nb_dof_border + ddl_loc_border;
                  int mV27 = inc + 26 * nb_dof_border + ddl_loc_border;
                  int mV28 = inc + 27 * nb_dof_border + ddl_loc_border;
                  int mV29 = inc + 28 * nb_dof_border + ddl_loc_border;
                  int mV30 = inc + 29 * nb_dof_border + ddl_loc_border;
                  int mV31 = inc + 30 * nb_dof_border + ddl_loc_border;
                  int mV32 = inc + 31 * nb_dof_border + ddl_loc_border;
                  int mW1 = inc + 32 * nb_dof_border + ddl_loc_border;
                  int mW2 = inc + 33 * nb_dof_border + ddl_loc_border;
                  int mW3 = inc + 34 * nb_dof_border + ddl_loc_border;
                  int mW4 = inc + 35 * nb_dof_border + ddl_loc_border;
                  int mW5 = inc + 36 * nb_dof_border + ddl_loc_border;
                  int mW6 = inc + 37 * nb_dof_border + ddl_loc_border;
                  int mW7 = inc + 38 * nb_dof_border + ddl_loc_border;
                  int mW8 = inc + 39 * nb_dof_border + ddl_loc_border;
                  int mW9 = inc + 40 * nb_dof_border + ddl_loc_border;
                  int mW10 = inc + 41 * nb_dof_border + ddl_loc_border;
                  int mW11 = inc + 42 * nb_dof_border + ddl_loc_border;
                  int mW12 = inc + 43 * nb_dof_border + ddl_loc_border;
                  int mW13 = inc + 44 * nb_dof_border + ddl_loc_border;
                  int mW14 = inc + 45 * nb_dof_border + ddl_loc_border;
                  int mW15 = inc + 46 * nb_dof_border + ddl_loc_border;
                  int mW16 = inc + 47 * nb_dof_border + ddl_loc_border;
                  int mW17 = inc + 48 * nb_dof_border + ddl_loc_border;
                  int mW18 = inc + 49 * nb_dof_border + ddl_loc_border;
                  int mW19 = inc + 50 * nb_dof_border + ddl_loc_border;
                  int mW20 = inc + 51 * nb_dof_border + ddl_loc_border;
                  int mW21 = inc + 52 * nb_dof_border + ddl_loc_border;
                  int mW22 = inc + 53 * nb_dof_border + ddl_loc_border;
                  int mW23 = inc + 54 * nb_dof_border + ddl_loc_border;
                  int mW24 = inc + 55 * nb_dof_border + ddl_loc_border;
                  int mW25 = inc + 56 * nb_dof_border + ddl_loc_border;
                  int mW26 = inc + 57 * nb_dof_border + ddl_loc_border;
                  int mW27 = inc + 58 * nb_dof_border + ddl_loc_border;
                  int mW28 = inc + 59 * nb_dof_border + ddl_loc_border;
                  int mW29 = inc + 60 * nb_dof_border + ddl_loc_border;
                  int mW30 = inc + 61 * nb_dof_border + ddl_loc_border;
                  int mW31 = inc + 62 * nb_dof_border + ddl_loc_border;
                  int mW32 = inc + 63 * nb_dof_border + ddl_loc_border;

                  K_tilde.Get(m, inc_water + loc_Vl) += eps_tilde * 3.0/(2*R) * Mh_surf(m);
                  K_tilde.Get(m, mV1) -= eps_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  int size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nV = inc + IndexLocalBorder(n);
                      K_tilde.Get(m, nV) += eps_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(m, mV1) -= eps_tilde * R * Mh_surf(m) * this->omega * this->omega;

                  K_tilde.Get(mV1, inc_water + loc_Vl) += 1.0;
                  K_tilde.Get(mV1, mV1) += -1.0 - alpha * R / 64.0 + R * rhoInfini * this->omega * this->omega / (64.0*alpha);
                  M_tilde.Get(mV1, mV1) += R / (64.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV1, mV2) += R * rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV1, mV2) += R / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV1, mV4) += R * rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV1, mV4) += R / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mV1, mV8) += R * rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mV1, mV8) += R / 8.0 * this->omega * this->omega;
                  K_tilde.Get(mV1, mV16) += R * rhoInfini * this->omega * this->omega / 16.0;
                  M_tilde.Get(mV1, mV16) += R / 16.0 * this->omega * this->omega;
                  K_tilde.Get(mV1, mV32) += R * rhoInfini * this->omega * this->omega / 32.0;
                  M_tilde.Get(mV1, mV32) += R / 32.0 * this->omega * this->omega;

                  K_tilde.Get(mV16, mV1) += 1.0;
                  K_tilde.Get(mV16, mV16) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV16, mV16) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV16, mV17) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV16, mV17) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV8, mV1) += 1.0;
                  K_tilde.Get(mV8, mV8) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mV8, mV8) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV8, mV9) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV8, mV9) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV8, mV25) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV8, mV25) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mV4, mV1) += 1.0;
                  K_tilde.Get(mV4, mV4) += - alpha  / 16.0 + rhoInfini * this->omega * this->omega / (16.0*alpha);
                  M_tilde.Get(mV4, mV4) += 1.0 / (16.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV4, mV5) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV4, mV5) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV4, mV13) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV4, mV13) += 1.0 / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mV4, mV29) += rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mV4, mV29) += 1.0 / 8.0 * this->omega * this->omega;

                  K_tilde.Get(mV2, mV1) += 1.0;
                  K_tilde.Get(mV2, mV2) += - alpha  / 32.0 + rhoInfini * this->omega * this->omega / (32.0*alpha);
                  M_tilde.Get(mV2, mV2) += 1.0 / (32.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV2, mV3) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV2, mV3) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV2, mV7) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV2, mV7) += 1.0 / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mV2, mV15) += rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mV2, mV15) += 1.0 / 8.0 * this->omega * this->omega;
                  K_tilde.Get(mV2, mV31) += rhoInfini * this->omega * this->omega / 16.0;
                  M_tilde.Get(mV2, mV31) += 1.0 / 16.0 * this->omega * this->omega;

                  K_tilde.Get(mV15, mV2) += 1.0;
                  K_tilde.Get(mV15, mV15) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV15, mV15) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV15, mV18) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV15, mV18) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV7, mV2) += 1.0;
                  K_tilde.Get(mV7, mV7) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mV7, mV7) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV7, mV10) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV7, mV10) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV7, mV26) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV7, mV26) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mV3, mV2) += 1.0;
                  K_tilde.Get(mV3, mV3) += - alpha  / 16.0 + rhoInfini * this->omega * this->omega / (16.0*alpha);
                  M_tilde.Get(mV3, mV3) += 1.0 / (16.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV3, mV6) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV3, mV6) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV3, mV14) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV3, mV14) += 1.0 / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mV3, mV30) += rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mV3, mV30) += 1.0 / 8.0 * this->omega * this->omega;

                  K_tilde.Get(mV14, mV3) += 1.0;
                  K_tilde.Get(mV14, mV14) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV14, mV14) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV14, mV19) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV14, mV19) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV6, mV3) += 1.0;
                  K_tilde.Get(mV6, mV6) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mV6, mV6) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV6, mV11) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV6, mV11) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV6, mV27) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV6, mV27) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mV13, mV4) += 1.0;
                  K_tilde.Get(mV13, mV13) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV13, mV13) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV13, mV20) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV13, mV20) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV5, mV4) += 1.0;
                  K_tilde.Get(mV5, mV5) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mV5, mV5) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV5, mV12) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV5, mV12) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mV5, mV28) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mV5, mV28) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mV12, mV5) += 1.0;
                  K_tilde.Get(mV12, mV12) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV12, mV12) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV12, mV21) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV12, mV21) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV11, mV6) += 1.0;
                  K_tilde.Get(mV11, mV11) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV11, mV11) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV11, mV22) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV11, mV22) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV10, mV7) += 1.0;
                  K_tilde.Get(mV10, mV10) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV10, mV10) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV10, mV23) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV10, mV23) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV9, mV8) += 1.0;
                  K_tilde.Get(mV9, mV9) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mV9, mV9) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mV9, mV24) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mV9, mV24) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mV32, mV1) += 1.0;
                  K_tilde.Get(mV32, mV32) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV32, mV32) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV31, mV2) += 1.0;
                  K_tilde.Get(mV31, mV31) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV31, mV31) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV30, mV3) += 1.0;
                  K_tilde.Get(mV30, mV30) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV30, mV30) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV29, mV4) += 1.0;
                  K_tilde.Get(mV29, mV29) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV29, mV29) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV28, mV5) += 1.0;
                  K_tilde.Get(mV28, mV28) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV28, mV28) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV27, mV6) += 1.0;
                  K_tilde.Get(mV27, mV27) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV27, mV27) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV26, mV7) += 1.0;
                  K_tilde.Get(mV26, mV26) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV26, mV26) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV25, mV8) += 1.0;
                  K_tilde.Get(mV25, mV25) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV25, mV25) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV24, mV9) += 1.0;
                  K_tilde.Get(mV24, mV24) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV24, mV24) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV23, mV10) += 1.0;
                  K_tilde.Get(mV23, mV23) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV23, mV23) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV22, mV11) += 1.0;
                  K_tilde.Get(mV22, mV22) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV22, mV22) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV21, mV12) += 1.0;
                  K_tilde.Get(mV21, mV21) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV21, mV21) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV20, mV13) += 1.0;
                  K_tilde.Get(mV20, mV20) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV20, mV20) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV19, mV14) += 1.0;
                  K_tilde.Get(mV19, mV19) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV19, mV19) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV18, mV15) += 1.0;
                  K_tilde.Get(mV18, mV18) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV18, mV18) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mV17, mV16) += 1.0;
                  K_tilde.Get(mV17, mV17) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mV17, mV17) += 1.0 / (2.0*alpha) * this->omega * this->omega;


                  K_tilde.Get(nb_dof + m, inc_water + loc_Wl) += mu_tilde * 3.0/(2*R) * Mh_surf(m);
                  K_tilde.Get(nb_dof + m, mW1) -= mu_tilde * (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                  size_row = Kh_surf.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh_surf.Index(m, j);
                      int nW = inc + nb_dof_border + IndexLocalBorder(n);
                      K_tilde.Get(nb_dof + m, nW) += mu_tilde * R * Kh_surf.Value(m, j) / 2.0;
                    }
                  M_tilde.Get(nb_dof + m, mW1) -= mu_tilde * R * Mh_surf(m) * this->omega * this->omega;

                  K_tilde.Get(mW1, inc_water + loc_Wl) += 1.0;
                  K_tilde.Get(mW1, mW1) += -1.0 - alpha * R / 64.0 + R * rhoInfini * this->omega * this->omega / (64.0*alpha);
                  M_tilde.Get(mW1, mW1) += R / (64.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW1, mW2) += R * rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW1, mW2) += R / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW1, mW4) += R * rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW1, mW4) += R / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mW1, mW8) += R * rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mW1, mW8) += R / 8.0 * this->omega * this->omega;
                  K_tilde.Get(mW1, mW16) += R * rhoInfini * this->omega * this->omega / 16.0;
                  M_tilde.Get(mW1, mW16) += R / 16.0 * this->omega * this->omega;
                  K_tilde.Get(mW1, mW32) += R * rhoInfini * this->omega * this->omega / 32.0;
                  M_tilde.Get(mW1, mW32) += R / 32.0 * this->omega * this->omega;

                  K_tilde.Get(mW16, mW1) += 1.0;
                  K_tilde.Get(mW16, mW16) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW16, mW16) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW16, mW17) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW16, mW17) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW8, mW1) += 1.0;
                  K_tilde.Get(mW8, mW8) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mW8, mW8) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW8, mW9) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW8, mW9) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW8, mW25) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW8, mW25) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mW4, mW1) += 1.0;
                  K_tilde.Get(mW4, mW4) += - alpha  / 16.0 + rhoInfini * this->omega * this->omega / (16.0*alpha);
                  M_tilde.Get(mW4, mW4) += 1.0 / (16.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW4, mW5) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW4, mW5) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW4, mW13) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW4, mW13) += 1.0 / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mW4, mW29) += rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mW4, mW29) += 1.0 / 8.0 * this->omega * this->omega;

                  K_tilde.Get(mW2, mW1) += 1.0;
                  K_tilde.Get(mW2, mW2) += - alpha  / 32.0 + rhoInfini * this->omega * this->omega / (32.0*alpha);
                  M_tilde.Get(mW2, mW2) += 1.0 / (32.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW2, mW3) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW2, mW3) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW2, mW7) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW2, mW7) += 1.0 / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mW2, mW15) += rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mW2, mW15) += 1.0 / 8.0 * this->omega * this->omega;
                  K_tilde.Get(mW2, mW31) += rhoInfini * this->omega * this->omega / 16.0;
                  M_tilde.Get(mW2, mW31) += 1.0 / 16.0 * this->omega * this->omega;

                  K_tilde.Get(mW15, mW2) += 1.0;
                  K_tilde.Get(mW15, mW15) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW15, mW15) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW15, mW18) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW15, mW18) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW7, mW2) += 1.0;
                  K_tilde.Get(mW7, mW7) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mW7, mW7) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW7, mW10) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW7, mW10) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW7, mW26) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW7, mW26) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mW3, mW2) += 1.0;
                  K_tilde.Get(mW3, mW3) += - alpha  / 16.0 + rhoInfini * this->omega * this->omega / (16.0*alpha);
                  M_tilde.Get(mW3, mW3) += 1.0 / (16.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW3, mW6) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW3, mW6) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW3, mW14) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW3, mW14) += 1.0 / 4.0 * this->omega * this->omega;
                  K_tilde.Get(mW3, mW30) += rhoInfini * this->omega * this->omega / 8.0;
                  M_tilde.Get(mW3, mW30) += 1.0 / 8.0 * this->omega * this->omega;

                  K_tilde.Get(mW14, mW3) += 1.0;
                  K_tilde.Get(mW14, mW14) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW14, mW14) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW14, mW19) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW14, mW19) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW6, mW3) += 1.0;
                  K_tilde.Get(mW6, mW6) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mW6, mW6) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW6, mW11) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW6, mW11) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW6, mW27) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW6, mW27) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mW13, mW4) += 1.0;
                  K_tilde.Get(mW13, mW13) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW13, mW13) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW13, mW20) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW13, mW20) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW5, mW4) += 1.0;
                  K_tilde.Get(mW5, mW5) += - alpha  / 8.0 + rhoInfini * this->omega * this->omega / (8.0*alpha);
                  M_tilde.Get(mW5, mW5) += 1.0 / (8.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW5, mW12) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW5, mW12) += 1.0 / 2.0 * this->omega * this->omega;
                  K_tilde.Get(mW5, mW28) += rhoInfini * this->omega * this->omega / 4.0;
                  M_tilde.Get(mW5, mW28) += 1.0 / 4.0 * this->omega * this->omega;

                  K_tilde.Get(mW12, mW5) += 1.0;
                  K_tilde.Get(mW12, mW12) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW12, mW12) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW12, mW21) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW12, mW21) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW11, mW6) += 1.0;
                  K_tilde.Get(mW11, mW11) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW11, mW11) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW11, mW22) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW11, mW22) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW10, mW7) += 1.0;
                  K_tilde.Get(mW10, mW10) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW10, mW10) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW10, mW23) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW10, mW23) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW9, mW8) += 1.0;
                  K_tilde.Get(mW9, mW9) += - alpha  / 4.0 + rhoInfini * this->omega * this->omega / (4.0*alpha);
                  M_tilde.Get(mW9, mW9) += 1.0 / (4.0*alpha) * this->omega * this->omega;
                  K_tilde.Get(mW9, mW24) += rhoInfini * this->omega * this->omega / 2.0;
                  M_tilde.Get(mW9, mW24) += 1.0 / 2.0 * this->omega * this->omega;

                  K_tilde.Get(mW32, mW1) += 1.0;
                  K_tilde.Get(mW32, mW32) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW32, mW32) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW31, mW2) += 1.0;
                  K_tilde.Get(mW31, mW31) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW31, mW31) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW30, mW3) += 1.0;
                  K_tilde.Get(mW30, mW30) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW30, mW30) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW29, mW4) += 1.0;
                  K_tilde.Get(mW29, mW29) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW29, mW29) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW28, mW5) += 1.0;
                  K_tilde.Get(mW28, mW28) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW28, mW28) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW27, mW6) += 1.0;
                  K_tilde.Get(mW27, mW27) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW27, mW27) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW26, mW7) += 1.0;
                  K_tilde.Get(mW26, mW26) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW26, mW26) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW25, mW8) += 1.0;
                  K_tilde.Get(mW25, mW25) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW25, mW25) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW24, mW9) += 1.0;
                  K_tilde.Get(mW24, mW24) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW24, mW24) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW23, mW10) += 1.0;
                  K_tilde.Get(mW23, mW23) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW23, mW23) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW22, mW11) += 1.0;
                  K_tilde.Get(mW22, mW22) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW22, mW22) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW21, mW12) += 1.0;
                  K_tilde.Get(mW21, mW21) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW21, mW21) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW20, mW13) += 1.0;
                  K_tilde.Get(mW20, mW20) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW20, mW20) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW19, mW14) += 1.0;
                  K_tilde.Get(mW19, mW19) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW19, mW19) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW18, mW15) += 1.0;
                  K_tilde.Get(mW18, mW18) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW18, mW18) += 1.0 / (2.0*alpha) * this->omega * this->omega;

                  K_tilde.Get(mW17, mW16) += 1.0;
                  K_tilde.Get(mW17, mW17) += - alpha  / 2.0 + rhoInfini * this->omega * this->omega / (2.0*alpha);
                  M_tilde.Get(mW17, mW17) += 1.0 / (2.0*alpha) * this->omega * this->omega;
                }
            }
        }
    }
}
