// initialisation du probleme
ModeEzNonLinearBGTEigenProblem
::ModeEzNonLinearBGTEigenProblem(int N, ModeEz_Solver& var0,
                              Vector<Matrix<Complexe, Prop, Storage> >& Mh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ch,
                              Matrix<Complexe, Prop, Storage>& Kh_surf,
                              VectReal_wp& Mh_surf,
                              Matrix<Complexe, Prop, Storage>& Rh_surf,
                              bool use_split)
  :
  var(var0), vec_Mh(Mh), vec_Kh(Kh), vec_Ch(Ch), Kh_bord(Kh_surf), Mh_bord(Mh_surf), Rh_bord(Rh_surf)
{

  // pour Ez et Hz, condition de BGT1
  nb_dof = N;

  // taille du pb aux valeurs propres
  nb_rows = N;
  int multi = 1;
  if (!var.only_Ez)
    {
      nb_rows += N;
      multi = 2;
    }

  this->Init(nb_rows);
  DISP(nb_rows); DISP(nb_dof);

  int Nref = 0;
  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    if (vec_Kh(ref).GetM() > 0)
      {
        Nref++;
      }

  Vector<FN> fct; PetscScalar coef0[1];
  if (var.type_CLA == var.BGT1)
    {
      this->SetSplitMatrices(multi * Nref + 2);
      fct.Reallocate(multi * Nref + 2);
      cout << "BGT1 : " << multi * Nref + 2 << " fonctions split." << endl;
    }
  else if (var.type_CLA == var.BGT2)
    {
      this->SetSplitMatrices(multi * Nref + 3);
      fct.Reallocate(multi * Nref + 3);
      cout << "BGT2 : " << multi * Nref + 3 << " fonctions split." << endl;
    }

  matMh_bord.Reallocate(nb_dof, nb_dof);
  for (int i = 0 ; i < nb_dof ; i++)
    {
      matMh_bord.Get(i,i) = Mh_surf(i);
    }

  // on calcule les poles du pb aux valeurs propres
  // VectComplex_wp s; split_reference.Clear();Condition = BGT1 Linear 2

  Mh_sum.Reallocate(nb_dof, nb_dof);
  if (!var.only_Ez)
    Mh_sumH.Reallocate(nb_dof, nb_dof);

  // f_0 = 1
  FNCreate(MPI_COMM_WORLD, &fct(0));
  FNSetType(fct(0), FNRATIONAL); coef0[0] = 1.0;
  FNRationalSetNumerator(fct(0), 1, coef0);
  FNRationalSetDenominator(fct(0), 0, NULL);

  DistributedMatrix<Complexe, Prop, Storage> A;
  int count = 1;
  for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
    {
      if (vec_Kh(ref).GetM() > 0)
        {
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp Delta = Iwp*var.omega*eps_sigma*var.ref_mu(ref);
          Complex_wp coef_mu = -Iwp*var.omega*var.ref_mu(ref);
          DISP(ref);


          // Première matrice A0
          Add(eps_sigma, vec_Mh(ref), Mh_sum);
          if (!var.only_Ez)
            Add(coef_mu, vec_Mh(ref), Mh_sumH);

          split_reference.PushBack(ref);

          PetscScalar num[1], DBeta[3];
          num[0] = 1.0;
          DBeta[0] = var.omega*var.omega;
          DBeta[1] = 0.0;
          DBeta[2] = -Delta;
          FNCreate(MPI_COMM_WORLD, &fct(count));
          FNSetType(fct(count), FNRATIONAL);
          FNRationalSetNumerator(fct(count), 1, num);
          FNRationalSetDenominator(fct(count), 3, DBeta);

          count++;

          if (!var.only_Ez)
            {
              PetscScalar num2[2];
              num2[1] = 0.0; num2[0] = Iwp*var.omega;

              FNCreate(MPI_COMM_WORLD, &fct(count));
              FNSetType(fct(count), FNRATIONAL);
              FNRationalSetNumerator(fct(count), 2, num2);
              FNRationalSetDenominator(fct(count), 3, DBeta);

              count++;
            }
        }
      }

    if (var.type_CLA == var.BGT1)
      {
        FN fDBeta, fSqrt, fSBeta, fUn, f1erTerme, f2eTerme;
        PetscScalar DBetaWater[3], un[1];
        // polynôme
        FNCreate(MPI_COMM_WORLD, &fDBeta);
        FNSetType(fDBeta, FNRATIONAL);
        DBetaWater[0] = var.omega * var.omega;
        DBetaWater[1] = 0.0;
        DBetaWater[2] = - var.rhoInfini * var.omega * var.omega;
        FNRationalSetNumerator(fDBeta, 3, DBetaWater);
        FNRationalSetDenominator(fDBeta, 0, NULL);
        // racine
        FNCreate(MPI_COMM_WORLD, &fSqrt);
        FNSetType(fSqrt, FNSQRT);
        // S(Beta)
        FNCreate(MPI_COMM_WORLD, &fSBeta);
        FNSetType(fSBeta, FNCOMBINE);
        FNCombineSetChildren(fSBeta, FN_COMBINE_COMPOSE, fDBeta, fSqrt);
        // eps_tilde
        FNCreate(MPI_COMM_WORLD, &fUn);
        FNSetType(fUn, FNRATIONAL);
        un[0] = 1.0;
        FNRationalSetNumerator(fUn, 1, un);
        FNRationalSetDenominator(fUn, 0, NULL);
        // premier terme
        FNCreate(MPI_COMM_WORLD, &f1erTerme);
        FNSetType(f1erTerme, FNCOMBINE);
        FNCombineSetChildren(f1erTerme, FN_COMBINE_DIVIDE, fUn, fSBeta);
        // deuxième terme
        FNCreate(MPI_COMM_WORLD, &f2eTerme);
        FNSetType(f2eTerme, FNCOMBINE);
        FNCombineSetChildren(f2eTerme, FN_COMBINE_DIVIDE, fUn, fDBeta);
        FNSetScale(f2eTerme, 1.0, 1.0 / (2.0 * var.ray));
        // ajout de la constante
        FNCreate(MPI_COMM_WORLD, &fct(multi * Nref + 1));
        FNSetType(fct(multi * Nref + 1), FNCOMBINE);
        FNCombineSetChildren(fct(multi * Nref + 1), FN_COMBINE_ADD, f1erTerme, f2eTerme);
        cout << "Fonctions split" << endl;
      }
    else if (var.type_CLA == var.BGT2)
      {
        FN fDBeta, fPoly2, fSqrt, fSBeta, fUn, fProd, fDenom1, fDenom2, fDenom3, f1erTerme, f2eTerme, f3eTerme, fTempD, fTempS, fSommeTemp;
        PetscScalar coef2[3], coef3[3], coef4[1], coef5[1], coef6[1];
        PetscScalar DBetaWater[3], un[1];
        // D(Beta)
        FNCreate(MPI_COMM_WORLD, &fDBeta);
        FNSetType(fDBeta, FNRATIONAL);
        DBetaWater[0] = var.omega * var.omega;
        DBetaWater[1] = 0.0;
        DBetaWater[2] = - var.rhoInfini * var.omega * var.omega;
        FNRationalSetNumerator(fDBeta, 3, DBetaWater);
        FNRationalSetDenominator(fDBeta, 0, NULL);
        // racine
        FNCreate(MPI_COMM_WORLD, &fSqrt);
        FNSetType(fSqrt, FNSQRT);
        // 1
        FNCreate(MPI_COMM_WORLD, &fUn);
        FNSetType(fUn, FNRATIONAL);
        un[0] = 1.0;
        FNRationalSetNumerator(fUn, 1, un);
        FNRationalSetDenominator(fUn, 0, NULL);
        // S(Beta)
        FNCreate(MPI_COMM_WORLD, &fSBeta);
        FNSetType(fSBeta, FNCOMBINE);
        FNCombineSetChildren(fSBeta, FN_COMBINE_COMPOSE, fDBeta, fSqrt);
        // S(Beta) * D(Beta)
        FNCreate(MPI_COMM_WORLD, &fProd);
        FNSetType(fProd, FNCOMBINE);
        FNCombineSetChildren(fProd, FN_COMBINE_MULTIPLY, fDBeta, fSBeta);
        FNSetScale(fProd, 1.0, var.ray);
        // D(Beta) + R S(Beta) * D(Beta)
        FNCreate(MPI_COMM_WORLD, &fDenom1);
        FNSetType(fDenom1, FNCOMBINE);
        FNCombineSetChildren(fDenom1, FN_COMBINE_ADD, fDBeta, fProd);
        // fct split devant K
        FNCreate(MPI_COMM_WORLD, &fct(multi * Nref + 2));
        FNSetType(fct(multi * Nref + 2), FNCOMBINE);
        FNCombineSetChildren(fct(multi * Nref + 2), FN_COMBINE_DIVIDE, fUn, fDenom1);
        FNSetScale(fct(multi * Nref + 2), 1.0, var.ray/2.0);

        // premier terme devant M
        FNCreate(MPI_COMM_WORLD, &f1erTerme);
        FNSetType(f1erTerme, FNCOMBINE);
        FNCombineSetChildren(f1erTerme, FN_COMBINE_DIVIDE, fUn, fDenom1);
        FNSetScale(f1erTerme, 1.0, 3.0 / (8.0*var.ray));
        //R D(Beta)
        FNDuplicate(fDBeta, MPI_COMM_WORLD, &fTempD);
        FNSetScale(fTempD, 1.0, var.ray);
        // S(Beta) + R * D(Beta)
        FNCreate(MPI_COMM_WORLD, &fDenom2);
        FNSetType(fDenom2, FNCOMBINE);
        FNCombineSetChildren(fDenom2, FN_COMBINE_ADD, fTempD, fSBeta);
        // deuxieme terme devant M
        FNCreate(MPI_COMM_WORLD, &f2eTerme);
        FNSetType(f2eTerme, FNCOMBINE);
        FNCombineSetChildren(f2eTerme, FN_COMBINE_DIVIDE, fUn, fDenom2);
        FNSetScale(f2eTerme, 1.0, 3.0 / 2.0);
        //R S(Beta)
        FNDuplicate(fSBeta, MPI_COMM_WORLD, &fTempS);
        FNSetScale(fTempS, 1.0, var.ray);
        // 1 + R * S(Beta)
        FNCreate(MPI_COMM_WORLD, &fDenom3);
        FNSetType(fDenom3, FNCOMBINE);
        FNCombineSetChildren(fDenom3, FN_COMBINE_ADD, fTempS, fUn);
        // troisième terme devant M
        FNCreate(MPI_COMM_WORLD, &f3eTerme);
        FNSetType(f3eTerme, FNCOMBINE);
        FNCombineSetChildren(f3eTerme, FN_COMBINE_DIVIDE, fUn, fDenom3);
        FNSetScale(f3eTerme, 1.0, var.ray);
        // somme des deux premiers termes
        FNCreate(MPI_COMM_WORLD, &fSommeTemp);
        FNSetType(fSommeTemp, FNCOMBINE);
        FNCombineSetChildren(fSommeTemp, FN_COMBINE_ADD, f1erTerme, f2eTerme);
        // fct split devant M
        FNCreate(MPI_COMM_WORLD, &fct(multi * Nref + 1));
        FNSetType(fct(multi * Nref + 1), FNCOMBINE);
        FNCombineSetChildren(fct(multi * Nref + 1), FN_COMBINE_ADD, fSommeTemp, f3eTerme);
      }

  this->SetSlepcFunction_NonRational(fct);
}

void ModeEzNonLinearBGTEigenProblem::ComputeOperator(const Complex_wp& L)
{
  cout << "Not implemented : ComputeOperator" << endl;
  abort();
}


// calcul de Y = T(L) X
void ModeEzNonLinearBGTEigenProblem::MltOperator(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  cout << "Not implemented : MltOperator" << endl;
  abort();
}

// calcul explicite de T(L) (utilise si on calcule explicitement les matrices)
void ModeEzNonLinearBGTEigenProblem
::ComputeOperatorExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "Not implemented : ComputeOperatorExplicit" << endl;
  abort();
}


void ModeEzNonLinearBGTEigenProblem::ComputeJacobian(const Complex_wp& L)
{
  cout << "Not implemented : ComputeJacobian" << endl;
  abort();
}

// calcul de Y = T'(L) X
void ModeEzNonLinearBGTEigenProblem::MltJacobian(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  cout << "Not implemented : MltJacobian" << endl;
  abort();
}

// calcul explicite de T'(L)
void ModeEzNonLinearBGTEigenProblem
::ComputeJacobianExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "Not implemented : ComputeJacobianExplicit" << endl;
  abort();
}

void ModeEzNonLinearBGTEigenProblem::ComputeExplicitPreconditioning(DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "ComputeExplicitPreconditioning" << endl;
  mat_lu.Factorize(A);
}

// factorisation de T(L)
void ModeEzNonLinearBGTEigenProblem::ComputePreconditioning(const Complex_wp& L)
{
  cout << "ComputePreconditioning" << endl;
  DistributedMatrix<Complex_wp, Prop, Storage> A;
  ComputeOperatorExplicit(L, A);

  mat_lu.Factorize(A);
}


void ModeEzNonLinearBGTEigenProblem::ComputePreconditioning(const VectComplex_wp& L, const VectComplex_wp& coef)
{
  cout << "Not implemented : ComputePreconditioning" << endl;
  abort();
}

void ModeEzNonLinearBGTEigenProblem::ComputeSplitPreconditioning(const Vector<int>& numL, const VectComplex_wp& coef)
{
  cout << "ComputeSplitPreconditioning" << endl;
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
void ModeEzNonLinearBGTEigenProblem::ApplyPreconditioning(const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
{
  cout << "ApplyPreconditioning" << endl;
  //cout << "appel apply" << Norm2(X) << endl;
  mat_lu.Solve(Y, X);
}


void ModeEzNonLinearBGTEigenProblem
::ComputeOperatorSplitExplicit(int num, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "ComputeOperatorSplitExplicit" << endl;
  int nbRef = split_reference.GetM();
  if (var.only_Ez)
    {
      if (num == 0)
        {
          A.Clear();
          A.Reallocate(nb_dof, nb_dof);
          for (int i = 0; i < nb_dof; i++)
            {
              int nb_val = Mh_sum.GetRowSize(i);
              A.ReallocateRow(i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(i, j) = Mh_sum.Index(i, j);
                  A.Value(i, j) = Mh_sum.Value(i, j);
                }
            }
        }
      else if (num == nbRef + 1)
        {
          Complex_wp eps_tilde = - Iwp * var.omega * var.ref_epsilon(var.ref_water) + var.ref_sigma(var.ref_water);
          A.Clear();
          A.Reallocate(nb_dof, nb_dof);
          for (int i = 0; i < nb_dof; i++)
            {
              A.ReallocateRow(i, 1);
              A.Index(i, 0) = i;
              A.Value(i, 0) = Mh_bord(i) * eps_tilde;
            }
        }
      else if (num == nbRef + 2)
        {
          Complex_wp eps_tilde = - Iwp * var.omega * var.ref_epsilon(var.ref_water) + var.ref_sigma(var.ref_water);
          A.Clear();
          A.Reallocate(nb_dof, nb_dof);
          for (int i = 0; i < nb_dof; i++)
            {
              int nb_val = Kh_bord.GetRowSize(i);
              A.ReallocateRow(i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(i, j) = Kh_bord.Index(i, j);
                  A.Value(i, j) = Kh_bord.Value(i, j) * eps_tilde;
                }
            }
        }
      else
        {
          int ref = this->split_reference(num-1);
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          A.Clear();
          A.Reallocate(nb_dof, nb_dof);
          for (int i = 0; i < nb_dof; i++)
            {
              int nb_val = vec_Kh(ref).GetRowSize(i);
              A.ReallocateRow(i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(i, j) = vec_Kh(ref).Index(i, j);
                  A.Value(i, j) = vec_Kh(ref).Value(i, j) * eps_sigma;
                }
            }
        }
    }
  else
    {
      if (num == 0)
        {
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          for (int i = 0; i < nb_dof; i++)
            {
              // partie (-i omega epsilon + sigma) M_h Ez
              int nb_val = Mh_sum.GetRowSize(i);
              A.ReallocateRow(i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(i, j) = Mh_sum.Index(i, j);
                  A.Value(i, j) = Mh_sum.Value(i, j);
                }

              // partie -i omega mu M_h H_z
              nb_val = Mh_sumH.GetRowSize(i);
              A.ReallocateRow(nb_dof + i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(nb_dof + i, j) = nb_dof + Mh_sumH.Index(i, j);
                  A.Value(nb_dof + i, j) = Mh_sumH.Value(i, j);
                }
            }
        }
      else if (num == 2 * nbRef + 1)
        {
          Complex_wp eps_tilde = - Iwp * var.omega * var.ref_epsilon(var.ref_water) + var.ref_sigma(var.ref_water);
          Complex_wp mu_tilde = - Iwp * var.omega * var.ref_mu(var.ref_water);
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          for (int i = 0; i < nb_dof; i++)
            {
              A.ReallocateRow(i, 1);
              A.Index(i, 0) = i;
              A.Value(i, 0) = Mh_bord(i) * eps_tilde;

              A.ReallocateRow(nb_dof + i, 1);
              A.Index(nb_dof + i, 0) = nb_dof + i;
              A.Value(nb_dof + i, 0) = Mh_bord(i) * mu_tilde;
            }
        }
      else if (num == 2 * nbRef + 2)
        {
          Complex_wp eps_tilde = - Iwp * var.omega * var.ref_epsilon(var.ref_water) + var.ref_sigma(var.ref_water);
          Complex_wp mu_tilde = - Iwp * var.omega * var.ref_mu(var.ref_water);
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          for (int i = 0; i < nb_dof; i++)
            {
              int nb_val = Kh_bord.GetRowSize(i);
              A.ReallocateRow(i, nb_val);
              A.ReallocateRow(nb_dof + i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(i, j) = Kh_bord.Index(i, j);
                  A.Value(i, j) = Kh_bord.Value(i, j) * eps_tilde;

                  A.Index(nb_dof + i, j) = nb_dof + Kh_bord.Index(i, j);
                  A.Value(nb_dof + i, j) = Kh_bord.Value(i, j) * mu_tilde;
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
          for (int i = 0; i < nb_dof; i++)
            {
              int nb_val = vec_Kh(ref).GetRowSize(i);
              A.ReallocateRow(i, nb_val);
              A.ReallocateRow(nb_dof + i, nb_val);
              for (int j = 0; j < nb_val; j++)
                {
                  A.Index(i, j) = vec_Kh(ref).Index(i, j);
                  A.Value(i, j) = eps_sigma*vec_Kh(ref).Value(i, j);

                  A.Index(nb_dof + i, j) = nb_dof + vec_Kh(ref).Index(i, j);
                  A.Value(nb_dof + i, j) = coef_mu*vec_Kh(ref).Value(i, j);
                }
            }
        }
      else
        {
          int ref = this->split_reference(num/2-1);
          A.Clear();
          A.Reallocate(nb_rows, nb_rows);
          if (ref != var.ref_water)
            {
              for (int i = 0; i < nb_dof; i++)
                {
                  // partie i beta/Delta C_h Hz
                  int nb_val = vec_Ch(ref).GetRowSize(i);
                  A.ReallocateRow(i, nb_val);
                  A.ReallocateRow(nb_dof + i, nb_val);
                  for (int j = 0; j < nb_val; j++)
                    {
                      A.Index(i, j) = nb_dof + vec_Ch(ref).Index(i, j);
                      A.Value(i, j) = vec_Ch(ref).Value(i, j);

                      A.Index(nb_dof + i, j) = vec_Ch(ref).Index(i, j);
                      A.Value(nb_dof + i, j) = -vec_Ch(ref).Value(i, j);
                    }
                }
            }
          else
            {
              for (int i = 0; i < nb_dof; i++)
                {
                  int n1 = vec_Ch(ref).GetRowSize(i);
                  int n2 = Rh_bord.GetRowSize(i);
                  Vector<int> col(n1+n2);
                  Vector<Complex_wp> val(n1+n2);
                  int num = 0;
                  for (int j0 = 0; j0 < n1; j0++)
                    {
                      int j = vec_Ch(ref).Index(i, j0);
                      col(num) = nb_dof + j;
                      val(num) = vec_Ch(ref).Value(i, j0);
                      num++;
                    }
                  A.AddInteractionRow(i, num, col, val);
                  num = 0;

                  for (int j0 = 0; j0 < n2; j0++)
                    {
                      int j = Rh_bord.Index(i, j0);
                      col(num) = nb_dof + j;
                      val(num) = Rh_bord.Value(i, j0);
                      num++;
                    }
                  A.AddInteractionRow(i, num, col, val);
                  num = 0;

                  for (int j0 = 0; j0 < n1; j0++)
                    {
                      int j = vec_Ch(ref).Index(i, j0);
                      col(num) = j;
                      val(num) = -vec_Ch(ref).Value(i, j0);
                      num++;
                    }
                  A.AddInteractionRow(nb_dof + i, num, col, val);
                  num = 0;

                  for (int j0 = 0; j0 < n2; j0++)
                    {
                      int j = Rh_bord.Index(i, j0);
                      col(num) = j;
                      val(num) = -Rh_bord.Value(i, j0);
                      num++;
                    }
                  A.AddInteractionRow(nb_dof + i, num, col, val);
                }
            }
        }
    }
}

void ModeEzNonLinearBGTEigenProblem
::MltOperatorSplit(int num, const SeldonTranspose& trans, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
{
  cout << "MltOperatorSplit" << endl;
  Y.Reallocate(nb_rows);
  Y.Zero();

  int nbRef = split_reference.GetM();

  if (var.only_Ez)
    {
      if (num == 0)
        Mlt(trans, Mh_sum, X, Y);
      else if (num == nbRef + 1)
        {
          Complex_wp eps_tilde = - Iwp * var.omega * var.ref_epsilon(var.ref_water) + var.ref_sigma(var.ref_water);
          // int ref = this->split_reference(num-1);
          // Mlt(trans, vec_Kh(ref), X, Y);
          Mlt(matMh_bord, X, Y);
          Mlt(eps_tilde, Y);
        }
      else if (num == nbRef + 2)
        {
          Complex_wp eps_tilde = - Iwp * var.omega * var.ref_epsilon(var.ref_water) + var.ref_sigma(var.ref_water);
          // int ref = this->split_reference(num-1);
          // Mlt(trans, vec_Kh(ref), X, Y);
          Mlt(Kh_bord, X, Y);
          Mlt(eps_tilde, Y);
        }
      else
        {
          int ref = this->split_reference(num-1);
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Mlt(trans, vec_Kh(ref), X, Y);
          Mlt(eps_sigma, Y);
        }
    }
  else
    {
      int N = this->nb_dof;
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
      else if (num == 2 * nbRef + 1)
        {
          Complex_wp eps_tilde = - Iwp * var.omega * var.ref_epsilon(var.ref_water) + var.ref_sigma(var.ref_water);
          Complex_wp mu_tilde = - Iwp * var.omega * var.ref_mu(var.ref_water);
          Mlt(matMh_bord, Ez, Y1);
          Mlt(eps_tilde, Y1);
          Mlt(matMh_bord, Hz, Y2);
          Mlt(mu_tilde, Y2);
        }
      else if (num == 2 * nbRef + 2)
        {
          Complex_wp eps_tilde = - Iwp * var.omega * var.ref_epsilon(var.ref_water) + var.ref_sigma(var.ref_water);
          Complex_wp mu_tilde = - Iwp * var.omega * var.ref_mu(var.ref_water);
          Mlt(Kh_bord, Ez, Y1);
          Mlt(eps_tilde, Y1);
          Mlt(Kh_bord, Hz, Y2);
          Mlt(mu_tilde, Y2);
        }
      else if (num%2 == 1)
        {
          int ref = this->split_reference(num/2);
          Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
          Complex_wp coef_mu = -Iwp*var.omega*var.ref_mu(ref);

          Mlt(vec_Kh(ref), Ez, Y1); Mlt(eps_sigma, Y1);
          Mlt(vec_Kh(ref), Hz, Y2); Mlt(coef_mu, Y2);
        }
      else
        {
          int ref = this->split_reference(num/2-1);
          Mlt(vec_Ch(ref), Hz, Y1);
          Mlt(vec_Ch(ref), Ez, Y2);
          if (ref == var.ref_water)
            {
              Vector<Complex_wp> tempE(N), tempH(N);
              tempE.Zero();
              tempH.Zero();
              Mlt(Rh_bord, Ez, tempE);
              Mlt(Rh_bord, Hz, tempH);
              Y1 += tempH;
              Y2 += tempE;
            }
          Y2 = -Y2;
        }

      for (int i = 0 ; i < N ; i++)
        {
          Y(i) = Y1(i);
          Y(N + i) = Y2(i);
        }
    }
}
