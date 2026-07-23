// initialisation du probleme
ModeUNonLinearEigenProblem
::ModeUNonLinearEigenProblem(int N, ModeU_Solver& var0,
                            Vector<Matrix<Complexe, Prop, Storage> >& Mh,
                            Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                            Matrix<Complexe, Prop, Storage>& Kh_surf,
                            VectReal_wp& Mh_surf,
                            bool use_split)
  :
  var(var0), vec_Mh(Mh), vec_Kh(Kh), Kh_bord(Kh_surf), Mh_bord(Mh_surf)
{
  // cout << "Not implemented : Constructor" << endl;
  // abort();
  nb_dof = N;
  type_CLA = var.type_CLA;
  this->Init(nb_dof);
  DISP(nb_dof);

  Vector<FN> fct; PetscScalar coef0[1], coef1[3];
  if (type_CLA == BGT1)
    {
      this->SetSplitMatrices(3);
      fct.Reallocate(3);
    }
  else if ((type_CLA == BGT2) || (type_CLA == NEUMANN))
    {
      this->SetSplitMatrices(4);
      fct.Reallocate(4);
    }

  A0.Reallocate(nb_dof, nb_dof);
  for (int ref = 0; ref < var.ref_rho.GetM(); ref++)
    if (vec_Kh(ref).GetM() > 0)
      {
        Add(1.0, vec_Kh(ref), A0);
        Add(- var.omega * var.omega * var.ref_rho(ref), vec_Mh(ref), A0);
      }
  FNCreate(MPI_COMM_WORLD, &fct(0));
  FNSetType(fct(0), FNRATIONAL); coef0[0] = 1.0;
  FNRationalSetNumerator(fct(0), 1, coef0);
  FNRationalSetDenominator(fct(0), 0, NULL);

  A1.Reallocate(nb_dof, nb_dof);
  for (int ref = 0; ref < var.ref_rho.GetM(); ref++)
    if (vec_Kh(ref).GetM() > 0)
      {
        Add(1.0, vec_Mh(ref), A1);
      }
  FNCreate(MPI_COMM_WORLD, &fct(1));
  FNSetType(fct(1), FNRATIONAL); coef1[0] = var.omega * var.omega; coef1[1] = 0.0; coef1[2] = 0.0;
  FNRationalSetNumerator(fct(1), 3, coef1);
  FNRationalSetDenominator(fct(1), 0, NULL);

  A2.Reallocate(nb_dof, nb_dof);
  for (int i = 0 ; i < nb_dof ; i++)
    {
      A2.Get(i,i) = Mh_surf(i);
    }
  if (type_CLA == BGT1)
    {
      FN fPoly, fSqrt, fTemp, fCst;
      PetscScalar coef2[3], coef4[1];
      // polynôme
      FNCreate(MPI_COMM_WORLD, &fPoly);
      FNSetType(fPoly, FNRATIONAL); coef2[0] = 1; coef2[1] = 0.0; coef2[2] = - var.rhoInfini;
      FNRationalSetNumerator(fPoly, 3, coef2);
      FNRationalSetDenominator(fPoly, 0, NULL);
      // racine
      FNCreate(MPI_COMM_WORLD, &fSqrt);
      FNSetType(fSqrt, FNSQRT);
      // 1/2R
      FNCreate(MPI_COMM_WORLD, &fCst);
      FNSetType(fCst, FNRATIONAL); coef4[0] = 1/(2.0 * var.ray);
      FNRationalSetNumerator(fCst, 1, coef4);
      FNRationalSetDenominator(fCst, 0, NULL);
      // première composition
      FNCreate(MPI_COMM_WORLD, &fTemp);
      FNSetType(fTemp, FNCOMBINE);
      FNCombineSetChildren(fTemp, FN_COMBINE_COMPOSE, fPoly, fSqrt);
      FNSetScale(fTemp, 1.0, var.omega);
      // ajout de la constante
      FNCreate(MPI_COMM_WORLD, &fct(2));
      FNSetType(fct(2), FNCOMBINE);
      FNCombineSetChildren(fct(2), FN_COMBINE_ADD, fTemp, fCst);
    }
  else if (type_CLA == BGT2)
    {
      // fct.Reallocate(4);

      FN fPoly, fPoly2, fSqrt, fTemp, fUn, fDenom, fCst1, fCst2, fNum, fQuotient;
      PetscScalar coef2[3], coef3[3], coef4[1], coef5[1], coef6[1];
      // polynôme
      FNCreate(MPI_COMM_WORLD, &fPoly);
      FNSetType(fPoly, FNRATIONAL); coef2[0] = var.omega * var.omega; coef2[1] = 0.0; coef2[2] = - var.rhoInfini * var.omega * var.omega;
      FNRationalSetNumerator(fPoly, 3, coef2);
      FNRationalSetDenominator(fPoly, 0, NULL);
      // racine
      FNCreate(MPI_COMM_WORLD, &fSqrt);
      FNSetType(fSqrt, FNSQRT);
      FNSetScale(fSqrt, 1.0, var.ray);
      // 1
      FNCreate(MPI_COMM_WORLD, &fUn);
      FNSetType(fUn, FNRATIONAL); coef4[0] = 1.0;
      FNRationalSetNumerator(fUn, 1, coef4);
      FNRationalSetDenominator(fUn, 0, NULL);
      // première composition (= R*S(\beta))
      FNCreate(MPI_COMM_WORLD, &fTemp);
      FNSetType(fTemp, FNCOMBINE);
      FNCombineSetChildren(fTemp, FN_COMBINE_COMPOSE, fPoly, fSqrt);
      // ajout de la constante
      FNCreate(MPI_COMM_WORLD, &fDenom);
      FNSetType(fDenom, FNCOMBINE);
      FNCombineSetChildren(fDenom, FN_COMBINE_ADD, fTemp, fUn);
      // inverse
      FNCreate(MPI_COMM_WORLD, &fct(3));
      FNSetType(fct(3), FNCOMBINE);
      FNCombineSetChildren(fct(3), FN_COMBINE_DIVIDE, fUn, fDenom);
      FNSetScale(fct(3), 1.0, var.ray/2.0);
      // FNSetScale(fct(3), 1.0, 1/(2.0*var.ray));

      FNCreate(MPI_COMM_WORLD, &fPoly2);
      // FNSetType(fPoly2, FNRATIONAL); coef3[0] = -var.ray * var.omega * var.omega; coef2[1] = 0.0; coef2[2] = var.ray * rhoInfini * var.omega * var.omega;
      FNSetType(fPoly2, FNRATIONAL); coef3[0] = -var.ray * var.omega * var.omega; coef3[1] = 0.0; coef3[2] = var.ray * var.rhoInfini * var.omega * var.omega;
      FNRationalSetNumerator(fPoly2, 3, coef3);
      FNRationalSetDenominator(fPoly2, 0, NULL);
      // 9/(8 R)
      FNCreate(MPI_COMM_WORLD, &fCst1);
      FNSetType(fCst1, FNRATIONAL); coef5[0] = 9.0 / (8.0 * var.ray);
      FNRationalSetNumerator(fCst1, 1, coef5);
      FNRationalSetDenominator(fCst1, 0, NULL);
      // 3 / (2 R)
      FNCreate(MPI_COMM_WORLD, &fCst2);
      FNSetType(fCst2, FNRATIONAL); coef6[0] = 3.0 / (2.0 * var.ray);
      FNRationalSetNumerator(fCst2, 1, coef6);
      FNRationalSetDenominator(fCst2, 0, NULL);
      // numérateur
      FNCreate(MPI_COMM_WORLD, &fNum);
      FNSetType(fNum, FNCOMBINE);
      FNCombineSetChildren(fNum, FN_COMBINE_ADD, fCst1, fPoly2);
      // quotient
      FNCreate(MPI_COMM_WORLD, &fQuotient);
      FNSetType(fQuotient, FNCOMBINE);
      FNCombineSetChildren(fQuotient, FN_COMBINE_DIVIDE, fNum, fDenom);
      FNSetScale(fQuotient, 1.0, -1.0);
      // f2
      FNCreate(MPI_COMM_WORLD, &fct(2));
      FNSetType(fct(2), FNCOMBINE);
      FNCombineSetChildren(fct(2), FN_COMBINE_ADD, fCst2, fQuotient);

      A3.Reallocate(nb_dof, nb_dof);
      Add(1.0, Kh_bord, A3);
    }

  this->SetSlepcFunction_NonRational(fct);
}

void ModeUNonLinearEigenProblem::ComputeOperator(const Complex_wp& L)
{
  cout << "Not implemented : ComputeOperator" << endl;
  abort();
}


// calcul de Y = T(L) X
void ModeUNonLinearEigenProblem::MltOperator(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
{

  if (isinf(abs(L)) || isnan(abs(L)))
    {
      DISP(L);
      abort();
    }

  if (type_CLA == BGT1)
    {
      Complex_wp one(1, 0);
      Complex_wp beta = L*var.omega;
      Real_wp R = var.ray;

      Y.Zero();

      for (int ref = 0; ref < var.ref_rho.GetM(); ref++)
        if (vec_Kh(ref).GetM() > 0)
          {
            MltAdd(one, vec_Kh(ref), X, one, Y);
            MltAdd(beta*beta - var.omega * var.omega * var.ref_rho(ref), vec_Mh(ref), X, one, Y);
          }
      for (int i = 0 ; i < nb_dof ; i++)
        {
          Y.Get(i) += (Sqrt(beta*beta - var.rhoInfini * var.omega * var.omega) + 1/(2*R)) * Mh_bord.Get(i) * X.Get(i);
        }
    }
  else if (type_CLA == BGT2)
    {
      Complex_wp one(1, 0);
      Complex_wp beta = L*var.omega;
      Real_wp R = var.ray;

      Y.Zero();

      for (int ref = 0; ref < var.ref_rho.GetM(); ref++)
        if (vec_Kh(ref).GetM() > 0)
          {
            MltAdd(one, vec_Kh(ref), X, one, Y);
            MltAdd(beta*beta - var.omega * var.omega * var.ref_rho(ref), vec_Mh(ref), X, one, Y);
          }
      MltAdd(R / (2.0 * (1.0 + R * Sqrt(beta*beta - var.rhoInfini * var.omega * var.omega))), Kh_bord, X, one, Y);
      for (int i = 0 ; i < nb_dof ; i++)
        {
          Y.Get(i) += (3.0/(2.0*R) - (9.0/(8.0*R) - R*(beta*beta - var.rhoInfini * var.omega * var.omega))/(1.0 + R * Sqrt(beta*beta - var.rhoInfini * var.omega * var.omega)) ) * Mh_bord.Get(i) * X.Get(i);
        }
    }
}

// calcul explicite de T(L) (utilise si on calcule explicitement les matrices)
void ModeUNonLinearEigenProblem
::ComputeOperatorExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "Not implemented : ComputeOperatorExplicit" << endl;
  abort();
  // A.Reallocate(nb_rows, nb_rows);
  //
  // /*DISP(L);
  // VectComplex_wp Ones(nb_rows), Aones(nb_rows);
  // ofstream file_out("B.dat"); file_out.precision(15);
  // for (int i = 0; i < nb_rows; i++)
  //   {
  //     Ones.Zero(); Ones(i) = 1.0;
  //     MltOperator(L, SeldonNoTrans, Ones, Aones);
  //     for (int j = 0; j < nb_rows; j++)
  //       if (abs(Aones(j)) > 1e-12)
  //         file_out << j+1 << " " << i+1 << " " << Aones(j) << "\n";
  //   }
  //
  //   file_out.close(); */
  //
  // Complex_wp beta = L*var.omega;
  // if (only_Ez)
  //   {
  //
  //     for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
  //       if (vec_Mh(ref).GetM() > 0)
  //         {
  //           Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
  //           Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
  //           Complex_wp coef = eps_sigma / Delta;
  //
  //           // on effectue eps_sigma Mh + coef Kh en eliminant les noeuds Dirichlet
  //           for (int i0 = 0; i0 < nb_dof_Ez; i0++)
  //             {
  //               int i = DofKeptDir(i0);
  //               int n1 = vec_Mh(ref).GetRowSize(i);
  //               int n2 = vec_Kh(ref).GetRowSize(i);
  //               Vector<int> col(n1+n2); Vector<Complex_wp> val(n1+n2); int num = 0;
  //               for (int j0 = 0; j0 < n1; j0++)
  //                 {
  //                   int j1 = vec_Mh(ref).Index(i, j0);
  //                   int j = IndexDirichlet(j1);
  //                   if (j >= 0)
  //                     {
  //                       col(num) = j;
  //                       val(num) = vec_Mh(ref).Value(i, j0)*eps_sigma;
  //                       num++;
  //                     }
  //                 }
  //
  //               for (int j0 = 0; j0 < n2; j0++)
  //                 {
  //                   int j1 = vec_Kh(ref).Index(i, j0);
  //                   int j = IndexDirichlet(j1);
  //                   if (j >= 0)
  //                     {
  //                       col(num) = j;
  //                       val(num) = vec_Kh(ref).Value(i, j0)*coef;
  //                       num++;
  //                     }
  //                 }
  //
  //               A.AddInteractionRow(i0, num, col, val);
  //             }
  //         }
  //   }
  // else
  //   {
  //     for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
  //       if (vec_Mh(ref).GetM() > 0)
  //         {
  //           Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
  //           Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
  //           Complex_wp coef_11 = eps_sigma/Delta;
  //           Complex_wp coef_12 = Iwp * beta / Delta;
  //
  //           // on traite les lignes avec Ez
  //           for (int i0 = 0; i0 < nb_dof_Ez; i0++)
  //             {
  //               int i = DofKeptDir(i0);
  //               int n1 = vec_Mh(ref).GetRowSize(i);
  //               int n2 = vec_Kh(ref).GetRowSize(i);
  //               int n3 = vec_Ch(ref).GetRowSize(i);
  //
  //               Vector<int> col(n1+n2+n3); Vector<Complex_wp> val(n1+n2+n3);
  //               int num = 0;
  //
  //               // partie eps_sigma Mh
  //               for (int j0 = 0; j0 < n1; j0++)
  //                 {
  //                   int j1 = vec_Mh(ref).Index(i, j0);
  //                   int j = IndexDirichlet(j1);
  //                   if (j >= 0)
  //                     {
  //                       col(num) = j;
  //                       val(num) = vec_Mh(ref).Value(i, j0)*eps_sigma;
  //                       num++;
  //                     }
  //                 }
  //
  //               A.AddInteractionRow(i0, num, col, val);
  //               num = 0;
  //
  //               // partie eps_sigma/Delta Kh
  //               for (int j0 = 0; j0 < n2; j0++)
  //                 {
  //                   int j1 = vec_Kh(ref).Index(i, j0);
  //                   int j = IndexDirichlet(j1);
  //                   if (j >= 0)
  //                     {
  //                       col(num) = j;
  //                       val(num) = vec_Kh(ref).Value(i, j0)*coef_11;
  //                       num++;
  //                     }
  //                 }
  //
  //               // partie i beta/Delta Ch
  //               for (int j0 = 0; j0 < n3; j0++)
  //                 {
  //                   int j = vec_Ch(ref).Index(i, j0);
  //                   col(num) = nb_dof_Ez + j;
  //                   val(num) = vec_Ch(ref).Value(i, j0)*coef_12;
  //                   num++;
  //                 }
  //
  //               A.AddInteractionRow(i0, num, col, val);
  //             }
  //
  //           Complex_wp coef_22 = - Iwp * var.omega * var.ref_mu(ref);
  //           Complex_wp coefKh_22 = - Iwp * var.omega * var.ref_mu(ref)/Delta;
  //           Complex_wp coef_21 = -Iwp * beta / Delta;
  //
  //           // on traite les lignes avec Hz
  //           for (int i = 0; i < nb_dof_Hz; i++)
  //             {
  //               int n1 = vec_Mh(ref).GetRowSize(i);
  //               int n2 = vec_Kh(ref).GetRowSize(i);
  //               int n3 = vec_Ch(ref).GetRowSize(i);
  //
  //               Vector<int> col(n1+n2+n3); Vector<Complex_wp> val(n1+n2+n3);
  //               int num = 0;
  //
  //               // partie -i omega mu Mh
  //               for (int j0 = 0; j0 < n1; j0++)
  //                 {
  //                   int j = vec_Mh(ref).Index(i, j0);
  //                   col(num) = nb_dof_Ez + j;
  //                   val(num) = vec_Mh(ref).Value(i, j0)*coef_22;
  //                   num++;
  //                 }
  //
  //               A.AddInteractionRow(nb_dof_Ez+i, num, col, val);
  //               num = 0;
  //
  //               // partie -i omega mu/Delta Kh
  //               for (int j0 = 0; j0 < n2; j0++)
  //                 {
  //                   int j = vec_Kh(ref).Index(i, j0);
  //                   col(num) = nb_dof_Ez + j;
  //                   val(num) = vec_Kh(ref).Value(i, j0)*coefKh_22;
  //                   num++;
  //                 }
  //
  //               // partie -i beta/Delta Ch
  //               for (int j0 = 0; j0 < n3; j0++)
  //                 {
  //                   int j1 = vec_Ch(ref).Index(i, j0);
  //                   int j = IndexDirichlet(j1);
  //                   if (j >= 0)
  //                     {
  //                       col(num) = j;
  //                       val(num) = vec_Ch(ref).Value(i, j0)*coef_21;
  //                       num++;
  //                     }
  //                 }
  //
  //               A.AddInteractionRow(nb_dof_Ez + i, num, col, val);
  //             }
  //         }
  //   }
  //A.WriteText("B.dat");
  //exit(0);
}


void ModeUNonLinearEigenProblem::ComputeJacobian(const Complex_wp& L)
{
  cout << "Not implemented : ComputeJacobian" << endl;
  abort();
}

// calcul de Y = T'(L) X
void ModeUNonLinearEigenProblem::MltJacobian(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& Xdir, Vector<Complex_wp>& Ydir)
{
  cout << "Not implemented : MltJacobian" << endl;
  abort();
  // Vector<Complex_wp> X(nb_dof_all), Y(nb_dof_all);
  // ExtractDirichlet(Xdir, X);
  // Y.Zero();
  //
  // Complex_wp one(1, 0);
  // Complex_wp beta = L*var.omega;
  // if (only_Ez)
  //   {
  //     for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
  //       if (vec_Mh(ref).GetM() > 0)
  //         {
  //           Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
  //           Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
  //           Complex_wp dDelta = 2.0*var.omega*beta;
  //
  //           MltAdd(-eps_sigma/(Delta*Delta)*dDelta, vec_Kh(ref), X, one, Y);
  //         }
  //   }
  // else
  //   {
  //     // le vecteur calculé est décomposé en deux : Ez, puis Hz
  //     int N = this->nb_dof_Hz;
  //     Vector<Complex_wp> Ez(N), Hz(N), Y1(N), Y2(N);
  //     Y1.Zero(); Y2.Zero();
  //     //cout << "Je passe à : MltJacobian" << endl;
  //     for (int i = 0 ; i < N ; i++)
  //       {
  //         Ez(i) = X(i);
  //         Hz(i) = X(N + i);
  //       }
  //
  //     for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
  //       if (vec_Mh(ref).GetM() > 0)
  //         {
  //
  //           Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
  //           Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
  //           Complex_wp dDelta = 2.0*var.omega*beta;
  //
  //           MltAdd(-eps_sigma/(Delta*Delta)*dDelta, vec_Kh(ref), Ez, one, Y1);
  //           MltAdd(Iwp * var.omega * var.ref_mu(ref)/(Delta * Delta) * dDelta, vec_Kh(ref), Hz, one, Y2);
  //           MltAdd((Iwp * var.omega * Delta - Iwp * beta * dDelta) / (Delta * Delta), vec_Ch(ref), Hz, one, Y1);
  //           MltAdd((-Iwp * var.omega * Delta + Iwp * beta * dDelta) / (Delta * Delta), vec_Ch(ref), Ez, one, Y2);
  //         }
  //
  //     for (int i = 0 ; i < N ; i++)
  //       {
  //         Y(i) = Y1(i);
  //         Y(N + i) = Y2(i);
  //       }
  //   }
  //
  // ExpandDirichlet(Y, Ydir);
}

// calcul explicite de T'(L)
void ModeUNonLinearEigenProblem
::ComputeJacobianExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  cout << "Not implemented : ComputeJacobianExplicit" << endl;
  abort();
  // A.Reallocate(nb_rows, nb_rows);
  //
  // Complex_wp beta = L*var.omega;
  // if (only_Ez)
  //   {
  //
  //     for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
  //       if (vec_Mh(ref).GetM() > 0)
  //         {
  //           Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
  //           Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
  //           Complex_wp dDelta = 2.0*var.omega*beta;
  //           Complex_wp coef = -eps_sigma / (Delta*Delta)*dDelta;
  //
  //           for (int i0 = 0; i0 < nb_dof_Ez; i0++)
  //             {
  //               int i = DofKeptDir(i0);
  //               int n2 = vec_Kh(ref).GetRowSize(i);
  //               Vector<int> col(n2); Vector<Complex_wp> val(n2); int num = 0;
  //               for (int j0 = 0; j0 < n2; j0++)
  //                 {
  //                   int j1 = vec_Kh(ref).Index(i, j0);
  //                   int j = IndexDirichlet(j1);
  //                   if (j >= 0)
  //                     {
  //                       col(num) = j;
  //                       val(num) = vec_Kh(ref).Value(i, j0)*coef;
  //                       num++;
  //                     }
  //                 }
  //
  //               A.AddInteractionRow(i0, num, col, val);
  //             }
  //         }
  //   }
  // else
  //   {
  //     for (int ref = 0; ref < var.ref_epsilon.GetM(); ref++)
  //       if (vec_Mh(ref).GetM() > 0)
  //         {
  //           Complex_wp eps_sigma = -Iwp*var.omega*var.ref_epsilon(ref) + var.ref_sigma(ref);
  //           Complex_wp Delta = -Iwp*var.omega*eps_sigma*var.ref_mu(ref) + beta*beta;
  //           Complex_wp dDelta = 2.0*var.omega*beta;
  //           Complex_wp coef_11 = -eps_sigma/(Delta*Delta)*dDelta;
  //           Complex_wp coef_12 = (Iwp * var.omega * Delta - Iwp * beta * dDelta) / (Delta * Delta);
  //
  //           // on traite les lignes avec Ez
  //           for (int i0 = 0; i0 < nb_dof_Ez; i0++)
  //             {
  //               int i = DofKeptDir(i0);
  //               int n2 = vec_Kh(ref).GetRowSize(i);
  //               int n3 = vec_Ch(ref).GetRowSize(i);
  //
  //               Vector<int> col(n2+n3); Vector<Complex_wp> val(n2+n3);
  //               int num = 0;
  //               for (int j0 = 0; j0 < n2; j0++)
  //                 {
  //                   int j1 = vec_Kh(ref).Index(i, j0);
  //                   int j = IndexDirichlet(j1);
  //                   if (j >= 0)
  //                     {
  //                       col(num) = j;
  //                       val(num) = vec_Kh(ref).Value(i, j0)*coef_11;
  //                       num++;
  //                     }
  //                 }
  //
  //               for (int j0 = 0; j0 < n3; j0++)
  //                 {
  //                   int j = vec_Ch(ref).Index(i, j0);
  //                   col(num) = nb_dof_Ez + j;
  //                   val(num) = vec_Ch(ref).Value(i, j0)*coef_12;
  //                   num++;
  //                 }
  //
  //               A.AddInteractionRow(i0, num, col, val);
  //             }
  //
  //           Complex_wp coef_22 = Iwp * var.omega * var.ref_mu(ref)/(Delta * Delta) * dDelta;
  //           Complex_wp coef_21 = (-Iwp * var.omega * Delta + Iwp * beta * dDelta) / (Delta * Delta);
  //
  //           // on traite les lignes avec Hz
  //           for (int i = 0; i < nb_dof_Hz; i++)
  //             {
  //               int n2 = vec_Kh(ref).GetRowSize(i);
  //               int n3 = vec_Ch(ref).GetRowSize(i);
  //
  //               Vector<int> col(n2+n3); Vector<Complex_wp> val(n2+n3);
  //               int num = 0;
  //
  //               for (int j0 = 0; j0 < n2; j0++)
  //                 {
  //                   int j = vec_Kh(ref).Index(i, j0);
  //                   col(num) = nb_dof_Ez + j;
  //                   val(num) = vec_Kh(ref).Value(i, j0)*coef_22;
  //                   num++;
  //                 }
  //
  //               for (int j0 = 0; j0 < n3; j0++)
  //                 {
  //                   int j1 = vec_Ch(ref).Index(i, j0);
  //                   int j = IndexDirichlet(j1);
  //                   if (j >= 0)
  //                     {
  //                       col(num) = j;
  //                       val(num) = vec_Ch(ref).Value(i, j0)*coef_21;
  //                       num++;
  //                     }
  //                 }
  //
  //               A.AddInteractionRow(nb_dof_Ez + i, num, col, val);
  //             }
  //         }
  //   }
}

void ModeUNonLinearEigenProblem::ComputeExplicitPreconditioning(DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  mat_lu.Factorize(A);
}

// factorisation de T(L)
void ModeUNonLinearEigenProblem::ComputePreconditioning(const Complex_wp& L)
{
  DistributedMatrix<Complex_wp, Prop, Storage> A;
  ComputeOperatorExplicit(L, A);

  mat_lu.Factorize(A);
}


void ModeUNonLinearEigenProblem::ComputePreconditioning(const VectComplex_wp& L, const VectComplex_wp& coef)
{
  cout << "Not implemented : ComputePreconditioning" << endl;
  abort();
  // Matrix<Complex_wp, Prop, Storage> A;
  //
  // VectComplex_wp beta(L.GetM());
  // for (int k = 0; k < L.GetM(); k++)
  //   beta(k) = var.omega*L(k);
  //
  // A.Reallocate(nb_dof, nb_dof);
  //
  // for (int ref = 0; ref < var.ref_rho.GetM(); ref++)
  //   if (vec_Mh(ref).GetM() > 0)
  //     {
  //       Complex_wp coefK, coefM;
  //       SetComplexZero(coefK); SetComplexZero(coefM);
  //       for (int k = 0; k < coef.GetM(); k++)
  //         {
  //           coefK += coef(k);
  //           coefM += coef(k) * (beta(k)*beta(k) - var.omega * var.omega * var.ref_rho(ref));
  //         }
  //
  //       // on effectue Kh + (beta^2 - w^2 rho) Mh
  //       for (int i = 0; i < nb_dof; i++)
  //         {
  //           int n1 = vec_Kh(ref).GetRowSize(i);
  //           int n2 = vec_Mh(ref).GetRowSize(i);
  //           Vector<int> col(n1+n2); Vector<Complex_wp> val(n1+n2); int num = 0;
  //           for (int j = 0; j < n1; j++)
  //             {
  //               int j1 = vec_Kh(ref).Index(i, j);
  //               if (j1 >= 0)
  //                 {
  //                   col(num) = j1;
  //                   val(num) = vec_Kh(ref).Value(i, j)*coefK;
  //                   num++;
  //                 }
  //             }
  //
  //           for (int j = 0; j < n2; j++)
  //             {
  //               int j1 = vec_Kh(ref).Index(i, j);
  //               if (j1 >= 0)
  //                 {
  //                   col(num) = j1;
  //                   val(num) = vec_Mh(ref).Value(i, j)*coefM;
  //                   num++;
  //                 }
  //             }
  //
  //           A.AddInteractionRow(i, num, col, val);
  //         }
  //     }
  //
  // cout << "REDONNER LE VRAI RHOINFINI !!!" << endl;
  // Complex_wp rhoInfini = 1.0;
  // Complex_wp coefMs; SetComplexZero(coefMs);
  // for (int k = 0; k < coef.GetM(); k++)
  //   {
  //     coefMs += coef(k) * (Sqrt(beta(k)*beta(k) - var.omega * var.omega * rhoInfini) + 1 / (2 * var.ray));
  //   }
  //
  // for (int i = 0; i < nb_dof; i++)
  //   {
  //     A.Get(i, i) += Mh_bord(i) * coefMs;
  //   }
  //
  // mat_lu.Factorize(A);
}

void ModeUNonLinearEigenProblem::ComputeSplitPreconditioning(const Vector<int>& numL, const VectComplex_wp& coef)
{
  // cout << "Not implemented : ComputeSplitPreconditioning" << endl;
  // abort();
  DistributedMatrix<Complex_wp, Prop, Storage> A, B;

  A.Reallocate(nb_dof, nb_dof);
  for (int k = 0; k < numL.GetM(); k++)
    {
      ComputeOperatorSplitExplicit(numL(k), B);
      Add(coef(k), B, A);
    }

  mat_lu.Factorize(A);
}

// application du preconditionneur Y = T(L)^{-1} X
void ModeUNonLinearEigenProblem::ApplyPreconditioning(const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
{
  //cout << "appel apply" << Norm2(X) << endl;
  mat_lu.Solve(Y, X);
}


void ModeUNonLinearEigenProblem
::ComputeOperatorSplitExplicit(int num, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
{
  // cout << "Not implemented : ComputeOperatorSplitExplicit" << endl;
  // abort();

  if (num == 0)
    {
      A.Clear();
      A.Reallocate(nb_dof, nb_dof);
      for (int i = 0; i < nb_dof; i++)
        {
          int nb_val = A0.GetRowSize(i);
          A.ReallocateRow(i, nb_val);
          for (int j = 0; j < nb_val; j++)
            {
              A.Index(i, j) = A0.Index(i, j);
              A.Value(i, j) = A0.Value(i, j);
            }
        }
    }
  else if (num == 1)
    {
      A.Clear();
      A.Reallocate(nb_dof, nb_dof);
      for (int i = 0; i < nb_dof; i++)
        {
          int nb_val = A1.GetRowSize(i);
          A.ReallocateRow(i, nb_val);
          for (int j = 0; j < nb_val; j++)
            {
              A.Index(i, j) = A1.Index(i, j);
              A.Value(i, j) = A1.Value(i, j);
            }
        }
    }
  else if (num == 2)
    {
      A.Clear();
      A.Reallocate(nb_dof, nb_dof);
      for (int i = 0; i < nb_dof; i++)
        {
          int nb_val = A2.GetRowSize(i);
          A.ReallocateRow(i, nb_val);
          for (int j = 0; j < nb_val; j++)
            {
              A.Index(i, j) = A2.Index(i, j);
              A.Value(i, j) = A2.Value(i, j);
            }
        }
    }
  else
    {
      A.Clear();
      A.Reallocate(nb_dof, nb_dof);
      for (int i = 0; i < nb_dof; i++)
        {
          int nb_val = A3.GetRowSize(i);
          A.ReallocateRow(i, nb_val);
          for (int j = 0; j < nb_val; j++)
            {
              A.Index(i, j) = A3.Index(i, j);
              A.Value(i, j) = A3.Value(i, j);
            }
        }
    }
}

void ModeUNonLinearEigenProblem
::MltOperatorSplit(int num, const SeldonTranspose& trans, const Vector<Complex_wp>& U, Vector<Complex_wp>& Y)
{
  // cout << "Not implemented : MltOperatorSplit" << endl;
  // abort();
  //
  Y.Zero();

  if (num == 0)
    {
      Mlt(A0, U, Y);
    }
  else if (num == 1)
    {
      Mlt(A1, U, Y);
    }
  else if (num == 2)
    {
      Mlt(A2, U, Y);
    }
  else
    {
      Mlt(A3, U, Y);
    }

}
