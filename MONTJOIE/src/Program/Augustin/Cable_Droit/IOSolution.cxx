
void ModeEz_Solver::SetInputData(const string& keyword, const VectString& parameters)
{
  // on appelle SetInputData pour mesh et mesh_num
  mesh.SetInputData(keyword, parameters);
  mesh_num.SetInputData(keyword, parameters);

  // autres mots-cles
  if (keyword == "FileMesh")
    name_mesh = PATH_AUGUSTIN + to_str("/Datas/") + parameters(0);
  else if (keyword == "Frequency")
    this->omega = 2.0*pi_wp*to_num<Real_wp>(parameters(0)) + to_num<Real_wp>(parameters(1));
  else if (keyword == "MateriauDielec")
    {
      // on garde le meme format que pour mode_maxwell (parameters(1) pas utilise)
      int ref = to_num<int>(parameters(0));
      ref_epsilon(ref) = to_num<Complex_wp>(parameters(2));
      ref_mu(ref) = to_num<Complex_wp>(parameters(3));
      ref_sigma(ref) = to_num<Complex_wp>(parameters(4));
    }
  else if (keyword == "Ray")
    {
      ray = to_num<Real_wp>(parameters(0));
    }
  else if ((keyword == "ReferenceDirichlet") || (keyword == "ReferenceBord"))
    {
      // liste des bords de type Dirichlet
      ref_bord.Reallocate(parameters.GetM());
      for (int i = 0; i < parameters.GetM(); i++)
        ref_bord(i) = to_num<int>(parameters(i));
    }
  else if (keyword == "ReferencePeriodic")
    {
      // liste des bords de type periodic
      ref_periodic(0) = to_num<int>(parameters(0));
      ref_periodic(1) = to_num<int>(parameters(1));

      // on rajoute la condition de maillage periodique
      mesh.AddPeriodicCondition(ref_periodic, BoundaryConditionEnum::PERIODIC_THETA);
      mesh_num.SetFormulationForPeriodicCondition(mesh_num.STRONG_PERIODIC);
    }
  else if (keyword == "NumberModePeriodic")
    {
      num_mode_periodic = to_num<int>(parameters(0));
    }
  else if (keyword == "LaplaceEigenvalue")
    {
      if (parameters(0) == "YES")
        compute_laplace_eigen = true;
      else
        compute_laplace_eigen = false;
    }
  else if (keyword == "SismoPlane")
    {
      // output on a plane
      int i = var_grid.GetM();
      var_grid.Resize(var_grid.GetM()+1);
      var_grid(i).SetInputData(keyword, parameters);
    }
  else if ( keyword == "SismoPointsFile")
    {
      // output on points given in a file
      int i = var_grid.GetM();
      var_grid.Resize(var_grid.GetM()+1);
      var_grid(i).SetInputData(keyword, parameters);
    }
  else if (keyword == "SismoLine")
    {
      // output on a plane
      int i = var_grid.GetM();
      var_grid.Resize(var_grid.GetM()+1);
      var_grid(i).SetInputData(keyword, parameters);
    }
  else if ( keyword == "SismoCircle")
    {
      // output on a circle
      int i = var_grid.GetM();
      var_grid.Resize(var_grid.GetM()+1);
      var_grid(i).SetInputData(keyword, parameters);
    }
  else if (keyword == "Condition")
    {
      cout << "Condition : " << parameters(0);

      if ((parameters(0) == "Neumann") || (parameters(0) == "Dirichlet") || (parameters(0) == "Full_Dirichlet"))
        {
          cout << " " << parameters(1) << endl;
        }
      else
        {
          cout << " " << parameters(1);
          if (parameters(1) == "Linear")
            cout << " " << parameters(2) << endl;
          else
            cout << endl;
        }

      if (parameters(0) == "Dirichlet")
        {
          type_CLA = DIRICHLET;
        }
      else if (parameters(0) == "Full_Dirichlet")
        {
          type_CLA = FULL_DIRICHLET;
        }
      else if (parameters(0) == "Neumann")
        {
          type_CLA = NEUMANN;
        }
      else if (parameters(0) == "BGT1")
        type_CLA = BGT1;
      else if (parameters(0) == "BGT2")
        type_CLA = BGT2;
      else
        {
          cout << "Unknown condition " << parameters(0) << endl;
          abort();
        }

      if (parameters(1) == "Linear")
        {
          linear = true;
          if ((type_CLA == BGT1) || ((type_CLA == BGT2)))
            linearization_step = to_num<int>(parameters(2));
          type_eigensolver = LINEAR_ARPACK;
        }
      else if (parameters(1) == "Non_Linear")
        {
          linear = false;
          type_eigensolver = NL_EIGS;
          use_split_formulation = true;
        }
      else
        {
          cout << "Unknown condition resolution " << parameters(1) << endl;
          abort();
        }
    }
  else if (keyword == "Eigensolver")
    {
      cout << "Choix du solveur : " << parameters(0) << endl;
      if (parameters(0) == "Linear" || parameters(0) == "Linear_Arpack")
        type_eigensolver = LINEAR_ARPACK;
      else if (parameters(0) == "Linear_Slepc")
        type_eigensolver = LINEAR_SLEPC;
      else if (parameters(0) == "Rii")
        type_eigensolver = NL_RII;
      else if (parameters(0) == "Slp")
        type_eigensolver = NL_SLP;
      else if (parameters(0) == "Nleigs")
        {
          type_eigensolver = NL_EIGS;
          linear = false;
        }
      else if (parameters(0) == "Direct")
        type_eigensolver = DIRECT;
      else
        {
          cout << "Unknown eigensolver " << parameters(0) << endl;
          abort();
        }

      for (int k = 1; k < parameters.GetM(); k++)
        {
          if (parameters(k) == "Split")
            use_split_formulation = true;
          if (parameters(k) == "Explicit")
            use_explicit_matrix = true;
          if (parameters(k) == "Default")
            use_default_petsc_solver = true;
        }
    }
  else if (keyword == "ShiftEigenvalue")
    {
      shift = to_num<Complex_wp>(parameters(0));
    }
  else if (keyword == "Intervalle")
    {
      minRe = to_num<Real_wp>(parameters(0));
      maxRe = to_num<Real_wp>(parameters(1));
      minIm = to_num<Real_wp>(parameters(2));
      maxIm = to_num<Real_wp>(parameters(3));
      cout << "Intervalle de recherche : [" << minRe << ";" << maxRe << "] + i [" << minIm << ";" << maxIm << "]" << endl;
    }
  else if (keyword == "NumberEigenvalue")
    {
      nb_asked_eigenval = to_num<int>(parameters(0));
    }
  else if (keyword == "CoordOtherFields")
    {
      if (parameters(0) == "Cartesian")
        type_coord_other_fields = CARTESIAN;
      else if (parameters(0) == "Polar")
        type_coord_other_fields = POLAR;
      else if (parameters(0) == "None")
        type_coord_other_fields = NONE;
      else
        {
          cout << "Unknown type of cooordonates " << parameters(0) << ". The other fields won't be computed." << endl;
        }
    }
  else if (keyword == "OnlyEz")
    {
      if (parameters(0) == "YES")
        only_Ez = true;
    }
  else if (keyword == "ReferenceConductor")
    {
      int N = parameters.GetM()/2;
      ref_conductor.Reallocate(N);
      ref_vol_conductor.Reallocate(N);
      phase_conductor.Reallocate(N);
      Complex_wp alpha = exp(Iwp*2.0*pi_wp/double(N));
      Complex_wp phase(1, 0);
      for (int i = 0; i < N; i++)
        {
          ref_conductor(i) = to_num<int>(parameters(2*i));
          ref_vol_conductor(i) = to_num<int>(parameters(2*i+1));
          phase_conductor(i) = phase;
          phase *= alpha;
        }
    }
  else if (keyword == "ReferenceIsolant")
    {
      int N = parameters.GetM();
      ref_vol_isolant.Reallocate(N);
      for (int i = 0; i < N; i++)
          ref_vol_isolant(i) = to_num<int>(parameters(i));
    }
}


// version simplifiee de InitGrid pour localiser les points du maillage sur les sorties
void ModeEz_Solver::InitGrid()
{
  // boundaries of the grid are set to boundaries of the domain if required
  for (int i = 0; i < var_grid.GetM(); i++)
    {
      var_grid(i).SetXmin0(this->mesh.GetXmin());
      var_grid(i).SetXmax0(this->mesh.GetXmax());
      var_grid(i).SetYmin0(this->mesh.GetYmin());
      var_grid(i).SetYmax0(this->mesh.GetYmax());
      var_grid(i).SetZmin0(this->mesh.GetZmin());
      var_grid(i).SetZmax0(this->mesh.GetZmax());
    }

  all_points_display.SetXmin(this->mesh.GetXmin());
  all_points_display.SetXmax(this->mesh.GetXmax());
  all_points_display.SetYmin(this->mesh.GetYmin());
  all_points_display.SetYmax(this->mesh.GetYmax());
  all_points_display.SetZmin(this->mesh.GetZmin());
  all_points_display.SetZmax(this->mesh.GetZmax());

  // predefined grids are generated
  // all_points_display is an union of all the predefined grids
  // so that the localization step is done once for the grid all_points_display
  for (int i = 0; i < var_grid.GetM(); i++)
    var_grid(i).InitGrid(all_points_display);

  // pre-localization of the elements of the mesh so that
  // each new point to search can be fastly found
  // implementation of this method is in GridInterpolation.cxx
  all_points_display.InitInterpolationGrid(this->mesh);

  // implementation of this method in file GridInterpolation.cxx
  all_points_display.LocalizePoints(this->mesh);

  // releasing memory used to prelocalize points
  all_points_display.ClearPrelocalizationArrays();

  // the grids are compressed
  // all_points_display.CompressGrid(var_grid);
}

template<class T>
void ModeEz_Solver::WriteOutputFile(const Vector<T>& x, const Vector<T>& y, int num_mode, Complex_wp beta)
{
  GridInterpolation<Dimension2>& var_interp = all_points_display;

  // calcul de la phase
  Complexe phase;
  SetComplexOne(phase);
  if (this->mesh_num.GetNbPeriodicDof() > 0)
    {
      Real_wp teta = this->mesh.GetPeriodicAlpha();
      to_complex(exp(Iwp*teta*double(this->num_mode_periodic)), phase);
      //DISP(teta); DISP(phase);
    }

  bool compute_grad = false;
  if (type_coord_other_fields != NONE)
    compute_grad = true;

  cout << "Write file" << endl;

  Vector<TinyVector<Real_wp, 1> > phi; VectR2 grad_phi, s;
  // boucle sur les grilles
  for (int n = 0; n < this->var_grid.GetM(); n++)
    {
      Vector<T> trace_Ez;
      Vector<T> trace_Hz;
      // si on prend les coordonnées polaires
      Vector<Complex_wp> trace_Hteta;
      Vector<Complex_wp> trace_Er;
      // si on prend les coordonnées cartésiennes
      Vector<Complex_wp> trace_Ex;
      Vector<Complex_wp> trace_Ey;
      Vector<Complex_wp> trace_Hx;
      Vector<Complex_wp> trace_Hy;
      GridInterpolationFull<Dimension2>& var_gr = this->var_grid(n); // grille de rendu final (pas maillage)
      const IVect& list_points = var_gr.GetPointNumber(); // points de la grille
      int nnz = list_points.GetM(); // nombre de points
      trace_Ez.Reallocate(nnz); trace_Ez.Zero();
      trace_Hz.Reallocate(nnz); trace_Hz.Zero();
      trace_Hteta.Reallocate(nnz); trace_Hteta.Zero();
      trace_Er.Reallocate(nnz); trace_Er.Zero();
      trace_Ex.Reallocate(nnz); trace_Ex.Zero();
      trace_Ey.Reallocate(nnz); trace_Ey.Zero();
      trace_Hx.Reallocate(nnz); trace_Hx.Zero();
      trace_Hy.Reallocate(nnz); trace_Hy.Zero();
      // boucle sur les points de la grille
      for (int i1 = 0; i1 < nnz; i1++)
        {
          int i = list_points(i1); // point i
          int iquad = var_interp.GetElementNumber(i); // polygone correspondant dans le maillage
          R2 point_loc = var_interp.GetLocalCoordinate(i); // coordonnées locales
          if ((iquad >= 0) && (iquad < this->mesh.GetNbElt()))
            {
              const IVect& num_ddl = this->mesh_num.Element(iquad).GetNodle();
              this->mesh.GetVerticesElement(iquad, s);
              const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(iquad);

              // we compute the values of basis functions at point_loc
              Fb.ComputeValuesPhi(point_loc, phi, var_interp.GetDFjm1(i),
                                  this->mesh_num, iquad);

              if (compute_grad)
                Fb.ComputeValuesGradientPhi(point_loc, grad_phi, var_interp.GetDFjm1(i),
                                            this->mesh_num, iquad);

              int ref = this->mesh.Element(iquad).GetReference();
              Complex_wp epsilon = this->ref_epsilon(ref);
              Complex_wp mu = this->ref_mu(ref);
              Complex_wp sigma = this->ref_sigma(ref);
              Complex_wp epsilon_tilde = epsilon + Iwp * sigma / this->omega;
              Complex_wp Delta = this->omega*this->omega*(- epsilon_tilde*mu + beta*beta);
              Complex_wp tempA = - Iwp * this->omega * epsilon_tilde / Delta;
              Complex_wp tempB = - Iwp * beta * this->omega / Delta;
              Complex_wp tempC = - Iwp * this->omega * mu / Delta;
              // cout << "Eps : " << epsilon << endl;
              // cout << "Mu : " << mu << endl;
              // cout << "Sigma : " << sigma << endl;
              // cout << "Eps_tild : " << epsilon_tilde << endl;
              // cout << "Delta : " << Delta << endl;
              // cout << "TempHt : " << tempHt << endl;
              // cout << "TempEr : " << tempEr << endl;

              // attention : pour un cas general il faut calculer Ex et Ey
              // et Hx, Hy
              // On trouve le vect r unitaire
              // DISP(var_interp.GetGlobalCoordinate(i));
              double x1 = var_interp.GetGlobalCoordinate(i)(0);
              double x2 = var_interp.GetGlobalCoordinate(i)(1);
              // DISP(point_loc); DISP(x1); DISP(x2);
              double norm = sqrt(x1*x1 + x2*x2);
              x1 /= norm;
              x2 /= norm;
              // DISP(norm); DISP(x1); DISP(x2);

              T valEz; SetComplexZero(valEz);
              T valHz; SetComplexZero(valHz);
              Complex_wp valEr; SetComplexZero(valEr);
              Complex_wp valHt; SetComplexZero(valHt);
              Complex_wp valEx; SetComplexZero(valEx);
              Complex_wp valEy; SetComplexZero(valEy);
              Complex_wp valHx; SetComplexZero(valHx);
              Complex_wp valHy; SetComplexZero(valHy);
              for (int i = 0; i < num_ddl.GetM(); i++)
                if (num_ddl(i) >= 0)
                  {
                    Complexe coef; SetComplexOne(coef);
                    if (this->is_dof_with_phase(iquad)(i))
                      coef = phase;

                    // si on se trouve dans l'élément
                    valEz += phi(i)(0)*coef*x(num_ddl(i));
                    if (!only_Ez)
                      valHz += phi(i)(0)*coef*y(num_ddl(i));

                    if (type_coord_other_fields == CARTESIAN)
                      {
                        if (only_Ez)
                          {
                            valEx += grad_phi(i)(0) * coef * x(num_ddl(i)) * tempB;
                            valEy += grad_phi(i)(1) * coef * x(num_ddl(i)) * tempB;
                            valHx -= grad_phi(i)(1) * coef * x(num_ddl(i)) * tempA;
                            valHy += grad_phi(i)(0) * coef * x(num_ddl(i)) * tempA;
                          }
                        else
                          {
                            valEx += grad_phi(i)(0) * coef * x(num_ddl(i)) * tempB;
                            valEy += grad_phi(i)(1) * coef * x(num_ddl(i)) * tempB;
                            valHx -= grad_phi(i)(1) * coef * x(num_ddl(i)) * tempA;
                            valHy += grad_phi(i)(0) * coef * x(num_ddl(i)) * tempA;

                            valEx += grad_phi(i)(1) * coef * y(num_ddl(i)) * tempC;
                            valEy -= grad_phi(i)(0) * coef * y(num_ddl(i)) * tempC;
                            valHx += grad_phi(i)(0) * coef * y(num_ddl(i)) * tempB;
                            valHy += grad_phi(i)(1) * coef * y(num_ddl(i)) * tempB;
                          }
                      }
                    else if (type_coord_other_fields == POLAR)
                      {
                        double dPhi_dr = x1 * grad_phi(i)(0) + x2 * grad_phi(i)(1);
                        valHt += dPhi_dr * coef * x(num_ddl(i)) * tempA;
                        valEr += dPhi_dr * coef * x(num_ddl(i)) * tempB;
                      }
                  }

              trace_Ez(i1) = valEz;
              trace_Hz(i1) = valHz;
              if (type_coord_other_fields == CARTESIAN)
                {
                  trace_Ex(i1) = valEx;
                  trace_Ey(i1) = valEy;
                  trace_Hx(i1) = valHx;
                  trace_Hy(i1) = valHy;
                }
              else if (type_coord_other_fields == POLAR)
                {
                  trace_Hteta(i1) = valHt;
                  trace_Er(i1) = valEr;
                }
            }
        }

      // on ecrit au format lisible par Python (loadND)
      string file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeEz" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
      cout << file_name << endl;
      WriteMatlab(trace_Ez, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
      if (!only_Ez)
        {
          file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeHz" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
          WriteMatlab(trace_Hz, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
        }

      if (type_coord_other_fields == CARTESIAN)
        {
          file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeEx" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
          WriteMatlab(trace_Ex, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
          file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeEy" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
          WriteMatlab(trace_Ey, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
          file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeHx" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
          WriteMatlab(trace_Hx, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
          file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeHy" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
          WriteMatlab(trace_Hy, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
        }
      else if (type_coord_other_fields == POLAR)
        {
          file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeEr" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
          WriteMatlab(trace_Er, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
          file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeHt" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
          WriteMatlab(trace_Hteta, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
        }
    }
}


template<class T>
void ModeEz_Solver::WriteOutputFileScalar(const Vector<T>& u, const string& name)
{
  GridInterpolation<Dimension2>& var_interp = all_points_display;

  // calcul de la phase
  Complexe phase;
  SetComplexOne(phase);
  if (this->mesh_num.GetNbPeriodicDof() > 0)
    {
      Real_wp teta = this->mesh.GetPeriodicAlpha();
      to_complex(exp(Iwp*teta*double(this->num_mode_periodic)), phase);
      //DISP(teta); DISP(phase);
    }

  Vector<TinyVector<Real_wp, 1> > phi; VectR2 grad_phi, s;
  // boucle sur les grilles
  for (int n = 0; n < this->var_grid.GetM(); n++)
    {
      Vector<T> trace_u, trace_du_dx, trace_du_dy;
      GridInterpolationFull<Dimension2>& var_gr = this->var_grid(n); // grille de rendu final (pas maillage)
      const IVect& list_points = var_gr.GetPointNumber(); // points de la grille
      int nnz = list_points.GetM(); // nombre de points
      trace_u.Reallocate(nnz); trace_u.Zero();
      trace_du_dx.Reallocate(nnz); trace_du_dx.Zero();
      trace_du_dy.Reallocate(nnz); trace_du_dy.Zero();
      // boucle sur les points de la grille
      for (int i1 = 0; i1 < nnz; i1++)
        {
          int i = list_points(i1); // point i
          int iquad = var_interp.GetElementNumber(i); // polygone correspondant dans le maillage
          R2 point_loc = var_interp.GetLocalCoordinate(i); // coordonnées locales
          if ((iquad >= 0) && (iquad < this->mesh.GetNbElt()))
            {
              const IVect& num_ddl = this->mesh_num.Element(iquad).GetNodle();
              this->mesh.GetVerticesElement(iquad, s);
              const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(iquad);

              // we compute the values of basis functions at point_loc
              Fb.ComputeValuesPhi(point_loc, phi, var_interp.GetDFjm1(i),
                                  this->mesh_num, iquad);

              Fb.ComputeValuesGradientPhi(point_loc, grad_phi, var_interp.GetDFjm1(i),
                                          this->mesh_num, iquad);

              Complex_wp val_u(0, 0);
              R2_Complex_wp grad_u;
              for (int i = 0; i < num_ddl.GetM(); i++)
                if (num_ddl(i) >= 0)
                  {
                    Complexe coef; SetComplexOne(coef);
                    if (this->is_dof_with_phase(iquad)(i))
                      coef = phase;

                    // si on se trouve dans l'élément
                    val_u += phi(i)(0)*coef*u(num_ddl(i));
                    Add(coef*u(num_ddl(i)), grad_phi(i), grad_u);
                  }

              trace_u(i1) = val_u;
              trace_du_dx(i1) = grad_u(0);
              trace_du_dy(i1) = grad_u(1);
            }
        }

      // on ecrit au format lisible par Python (loadND)
      string root = GetBaseString(name);
      string file_name = root + "_G" + to_str(n) + ".dat";
      string file_name_dx = root + "_G" + to_str(n) + "_DX.dat";
      string file_name_dy = root + "_G" + to_str(n) + "_DY.dat";
      WriteMatlab(trace_u, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
      WriteMatlab(trace_du_dx, var_gr, file_name_dx, OutputTypeEnum::DOUBLE_PRECISION, false);
      WriteMatlab(trace_du_dy, var_gr, file_name_dy, OutputTypeEnum::DOUBLE_PRECISION, false);
    }
}
