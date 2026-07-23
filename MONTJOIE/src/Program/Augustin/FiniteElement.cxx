

// calcul des matrices elements finis
// À COMMENTER POUR BIEN COMPRENDRE
void ModeEz_Solver::ComputeFemMatrices(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                        VectComplex_wp& vec_Vk,
                        Vector<VectComplex_wp >& vec_Vk2,
                        int ref_surf,
                        VectReal_wp& Mh_surf,
                        Matrix<Complexe, Prop, Storage>& Kh_surf,
                        Matrix<Complexe, Prop, Storage>& Rh_surf)
{
  // on alloue les matrices
  int N = this->mesh_num.GetNbDof();
  vec_Mh.Reallocate(this->ref_epsilon.GetM());
  vec_Kh.Reallocate(this->ref_epsilon.GetM());
  vec_Ch.Reallocate(this->ref_epsilon.GetM());
  Lk_conductor.Reallocate(this->ref_epsilon.GetM());
  Lk_conductor.Zero();
  vec_Lk.Reallocate(this->ref_epsilon.GetM());
  vec_Lk.Zero();
  vec_Vk.Reallocate(N);
  vec_Vk.Zero();
  vec_Vk2.Reallocate(this->ref_epsilon.GetM());

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  for (int ref = 0; ref < index_used.GetM(); ref++)
    if (index_used(ref))
      {
        vec_Mh(ref).Reallocate(N, N);
        vec_Kh(ref).Reallocate(N, N);
        if (!only_Ez)
          vec_Ch(ref).Reallocate(N, N);
        vec_Vk2(ref).Reallocate(N);
      }

  VectR2 s;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;

  // calcul de la phase
  Complexe phase;
  SetComplexOne(phase);
  if (this->mesh_num.GetNbPeriodicDof() > 0)
    {
      Real_wp teta = this->mesh.GetPeriodicAlpha();
      to_complex(exp(Iwp*teta*double(this->num_mode_periodic)), phase);
    }

  Real_wp threshold = 1e-15;

  // boucle sur les elements
  VectReal_wp Bzero, Bmass; Matrix2_2 dfjm1, Ctmp, Cref;
  VectR2 Dzero, Ezero; Vector<Matrix2_2> Astiff, Cstiff;
  TinyVector<bool, 4> null_term(true, false, true, true);
  Matrix<Real_wp> mat_elem_Mh, mat_elem_Kh, mat_elem_Ch;
  Cref(0, 1) = -1.0; Cref(1, 0) = 1.0;
  DISP(this->mesh.GetNbElt());
  num_conduc_ddl.Reallocate(N); num_conduc_ddl.Fill(-1);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(i);
      int ref = this->mesh.Element(i).GetReference();
      // on recupere le numero du conducteur si present
      int num_conductor = -1;
      for (int k = 0; k < ref_vol_conductor.GetM(); k++)
        if (ref_vol_conductor(k) == ref)
          num_conductor = k;

      // calcul des matrices jacobiennes sur les points de quadrature
      this->mesh.GetVerticesElement(i, s);
      Fb.FjElemQuadrature(s, PointsElem, this->mesh, i);
      Fb.DFjElemQuadrature(s, PointsElem, MatricesElem, this->mesh, i);

      int Nquad = Fb.GetNbPointsQuadratureInside();
      Bmass.Reallocate(Nquad); Astiff.Reallocate(Nquad);
      Cstiff.Reallocate(Nquad);

      // on recupere le numero de conducteur pour chaque ddl
      VectReal_wp phi;
      int nb_dof = Fb.GetNbDof();
      IVect num_ddl = this->mesh_num.Element(i).GetNodle();
      if (num_conductor >= 0)
        for (int k = 0; k < nb_dof; k++)
          if (num_ddl(k) >= 0)
            num_conduc_ddl(num_ddl(k)) = num_conductor;

      // boucle sur les points de quadrature
      for (int j = 0; j < Nquad; j++)
        {
          // pour le calcul de Mh (matrice de masse)
          Real_wp weight = Fb.WeightsND(j);
          Real_wp jacob = Det(MatricesElem.GetPointQuadrature(j));
          Bmass(j) = jacob*weight;

          // calcul de Vk
          // pour avoir phi_i(xi_j) appeler Fb.GetValuePhiOnQuadraturePoint(j, phi);

          vec_Lk(ref) += jacob*weight;
          Fb.GetValuePhiOnQuadraturePoint(j, phi);
          for (int k = 0; k < nb_dof; k++)
            if (num_ddl(k) >= 0)
              vec_Vk2(ref)(num_ddl(k)) += jacob*weight*phi(k);

          if (num_conductor >= 0)
            {
              Lk_conductor(num_conductor) += jacob*weight;
              Fb.GetValuePhiOnQuadraturePoint(j, phi);
              for (int k = 0; k < nb_dof; k++)
                if (num_ddl(k) >= 0)
                  vec_Vk(num_ddl(k)) += jacob*weight*phi(k);
            }

          // pour le calcul de Kh (matrice de rigidite)
          GetInverse(MatricesElem.GetPointQuadrature(j), dfjm1);
          MltTrans(dfjm1, dfjm1, Astiff(j));
          Mlt(jacob*weight, Astiff(j));

          // pour le calcul de Ch (matrice de couplage)
          if (!only_Ez)
            {
              MltTrans(Cref, dfjm1, Ctmp);
              Mlt(dfjm1, Ctmp, Cstiff(j));
              Mlt(jacob*weight, Cstiff(j));
            }
        }

      // calcul des matrices de masse et rigidite elementaires
      mat_elem_Mh.Reallocate(nb_dof, nb_dof); mat_elem_Mh.Zero();
      mat_elem_Kh.Reallocate(nb_dof, nb_dof); mat_elem_Kh.Zero();
      if (!only_Ez)
        {
          mat_elem_Ch.Reallocate(nb_dof, nb_dof);
          mat_elem_Ch.Zero();
        }

      Fb.AddVariableMassMatrix(0, 0, Bmass, mat_elem_Mh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Astiff, Dzero, Ezero, null_term, mat_elem_Kh);
      if (!only_Ez)
        Fb.AddVariableElemMatrix(0, 0, Bzero, Cstiff, Dzero, Ezero, null_term, mat_elem_Ch);

      // assemblage des matrices
      IVect permut(nb_dof); permut.Fill();
      Sort(num_ddl, permut);

      Vector<int> col(nb_dof); Vector<Complexe> val(nb_dof); int nb_val = 0;
      col.Fill(-1); val.Zero();
      for (int j = 0; j < nb_dof; j++)
        if (num_ddl(j) >= 0)
          {
            nb_val = 0; Complexe coef1; SetComplexOne(coef1);
            if (this->is_dof_with_phase(i)(permut(j)))
              coef1 *= conjugate(phase);

            //DISP(i); DISP(j); DISP(coef1); DISP(num_ddl(j));
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Mh(permut(j), permut(k))) >= threshold))
                {
                  Complex_wp coef = coef1;
                  if (this->is_dof_with_phase(i)(permut(k)))
                    coef *= phase;

                  col(nb_val) = num_ddl(k);
                  val(nb_val) = coef*mat_elem_Mh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Mh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Kh(permut(j), permut(k))) >= threshold))
                {
                  Complexe coef = coef1;
                  if (this->is_dof_with_phase(i)(permut(k)))
                    coef *= phase;

                  col(nb_val) = num_ddl(k);
                  val(nb_val) = coef*mat_elem_Kh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Kh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            if (!only_Ez)
              {
                nb_val = 0;
                for (int k = 0; k < nb_dof; k++)
                  if ((num_ddl(k) >= 0) && (abs(mat_elem_Ch(permut(j), permut(k))) >= threshold))
                    {
                      Complexe coef = coef1;
                      if (this->is_dof_with_phase(i)(permut(k)))
                        coef *= phase;

                      col(nb_val) = num_ddl(k);
                      val(nb_val) = coef*mat_elem_Ch(permut(j), permut(k));
                      nb_val++;
                    }

                vec_Ch(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);
              }
          }
    }

  for (int ref = 0; ref < index_used.GetM(); ref++)
    if (index_used(ref))
      {
        if (!only_Ez)
          vec_Ch(ref).RemoveSmallEntry(threshold);
        for (int i = 0; i < N; i++)
          {
            Real_wp Lk = vec_Lk(ref);
            vec_Vk2(ref)(i) /= Lk;
          }
      }

  DISP(Lk_conductor);
  DISP(vec_Lk);
  // on divise les Vk par l'aire pour avoir V_i = 1/|L_k| \int_{L_k} \varphi_i dx
  for (int i = 0; i < N; i++)
    if (num_conduc_ddl(i) >= 0)
      {
        Real_wp Lk = Lk_conductor(num_conduc_ddl(i));
        vec_Vk(i) /= Lk;
      }

  // boucle sur les aretes (matrices surfaciques)
  Mh_surf.Reallocate(N);
  Mh_surf.Zero();
  Kh_surf.Clear(); Kh_surf.Reallocate(N, N);
  Kh_surf.Zero();
  Rh_surf.Clear(); Rh_surf.Reallocate(N, N);
  Rh_surf.Zero();
  // vec_theta.Reallocate(N);
  // vec_theta.Zero();
  cout << "Calcul matrices surfaciques" << endl;
  for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
    {
      int ref = this->mesh.BoundaryRef(i).GetReference();
      if (ref != ref_surf)
        continue;

      int num_elem = this->mesh.BoundaryRef(i).numElement(0);
      IVect num_ddl = this->mesh_num.Element(num_elem).GetNodle();
      int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
      const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(num_elem);

      this->mesh.GetVerticesElement(num_elem, s);
      Fb.FjElemQuadrature(s, PointsElem, this->mesh, num_elem);
      Fb.DFjElemQuadrature(s, PointsElem, MatricesElem, this->mesh, num_elem);

      // restriction of DF_i on surface
      Fb.FjSurfaceElem(s, PointsElem, mesh, num_elem, num_loc);
      Fb.DFjSurfaceElem(s, PointsElem, MatricesElem, mesh, num_elem, num_loc);

      // on calcule les gradients surfaciques des fonctions de base sur les points de quadrature
      // et la matrice de masse surfacique
      int nb_points_quadrature_edge = Fb.GetNbQuadBoundary(num_loc);
      int Ns = Fb.GetNbDofBoundary(num_loc);
      Matrix<R2> grad_phi_surf(Ns, nb_points_quadrature_edge); // gradient fcts de base sur la surface
      Matrix<Real_wp> dPhi_ds(Ns, nb_points_quadrature_edge); // gradient fcts de base sur la surface
      VectR2 grad_phi; R2 grad;
      for (int k = 0; k < nb_points_quadrature_edge; k++)
        {
          // jacobien surfacique
          Real_wp ds = MatricesElem.GetDsQuadratureBoundary(k);
          // poids d'integration sur l'intervalle unite
          Real_wp poids = Fb.WeightsQuadratureBoundary(k, num_loc);
          int num_dof_loc = Fb.GetLocalNumber(num_loc, k);
          int num_dof = num_ddl(num_dof_loc);
          Mh_surf(num_dof) += poids*ds;

          // // point de quadrature
          // // DISP(PointsElem.GetPointQuadratureBoundary(k));
          // R2 pointQuadGlob = PointsElem.GetPointQuadratureBoundary(k);
          // vec_theta(num_dof) = atan2(pointQuadGlob);

          // normale sur le point
          R2 normale = MatricesElem.GetNormaleQuadratureBoundary(k);

          // tangente sur le point
          R2 tangente(-normale(1), normale(0));

          // on calcule grad(phi) sur le point de quadrature
          GetInverse(MatricesElem.GetPointQuadratureBoundary(k), dfjm1);
          int num_point = Fb.GetQuadNumber(num_loc, k);
          Fb.GetGradientPhiOnQuadraturePoint(num_point, grad_phi);
          for (int i0 = 0; i0 < Ns; i0++)
            {
              num_dof_loc = Fb.GetLocalNumber(num_loc, i0);
              MltTrans(dfjm1, grad_phi(num_dof_loc), grad);
              Real_wp grad_dot_n = DotProd(grad, normale);
              grad -= grad_dot_n*normale;
              grad_phi_surf(i0, k) = grad;
              dPhi_ds(i0, k) = DotProd(grad, tangente);
            }
        }

      // on calcule la matrice de rigidite surfacique
      for (int k = 0; k < nb_points_quadrature_edge; k++)
        {
          // jacobien surfacique
          Real_wp ds = MatricesElem.GetDsQuadratureBoundary(k);
          // poids d'integration sur l'intervalle unite
          Real_wp poids = Fb.WeightsQuadratureBoundary(k, num_loc);

          for (int i0 = 0; i0 < Ns; i0++)
            for (int j = 0; j < Ns; j++)
              {
                int num_dof_loc_i0 = Fb.GetLocalNumber(num_loc, i0);
                int num_dof_i0 = num_ddl(num_dof_loc_i0);
                int num_dof_loc_j = Fb.GetLocalNumber(num_loc, j);
                int num_dof_j = num_ddl(num_dof_loc_j);
                Real_wp vloc = DotProd(grad_phi_surf(i0, k), grad_phi_surf(j, k));
                Kh_surf.AddInteraction(num_dof_i0, num_dof_j, vloc*poids*ds);
              }

          for (int j = 0; j < Ns; j++)
            {
              int num_dof_loc_k = Fb.GetLocalNumber(num_loc, k);
              int num_dof_k = num_ddl(num_dof_loc_k);
              int num_dof_loc_j = Fb.GetLocalNumber(num_loc, j);
              int num_dof_j = num_ddl(num_dof_loc_j);


              Real_wp grad_vloc = dPhi_ds(j, k);
              Rh_surf.AddInteraction(num_dof_k, num_dof_j, grad_vloc*poids*ds);
            }
        }
    }
}
