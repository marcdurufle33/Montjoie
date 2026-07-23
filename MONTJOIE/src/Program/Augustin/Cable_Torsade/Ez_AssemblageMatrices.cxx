
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
          M_tilde.Get(m, nW) += - Iwp * this->omega * R * Rh_surf.Value(m, j);
          M_tilde.Get(nb_dof + m, nV) += Iwp * this->omega * R * Rh_surf.Value(m, j);
        }
    }

  if (this->type_CLA != NEUMANN)
    {
      DISP(rhoInfini);
      Complex_wp alpha = this->omega * Sqrt(rhoPlus - rhoInfini);
      DISP(alpha);
      Real_wp R = this->ray;
      Complex_wp eps_tilde = - Iwp * this->omega * this->ref_epsilon(this->ref_water) + this->ref_sigma(this->ref_water);
      Complex_wp mu_tilde = - Iwp * this->omega * this->ref_mu(this->ref_water);
      for (int ddl_loc_border = 0 ; ddl_loc_border < nb_dof_border ; ddl_loc_border++)
        {
          int m = NumGlobBorder(ddl_loc_border); // m est le numéro global
          int loc_Vl = nb_dof_couche(this->ref_water) + IndexDof(this->ref_water)(m);
          int loc_Wl = 3 * nb_dof_couche(this->ref_water) +IndexDof(this->ref_water)(m);

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
            }
        }
    }
}
