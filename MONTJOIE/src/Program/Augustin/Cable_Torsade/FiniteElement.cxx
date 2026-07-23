// calcul des matrices elements finis
// À COMMENTER POUR BIEN COMPRENDRE
void ModeEs_Solver::ComputeFemMatricesNonLin(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                        Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Kr2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Cr2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjh,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjr2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Rir2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjh,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjr2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Sir2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Dr2h,
                        VectR2& vec_Pos
                      )
{
  // on alloue les matrices
  int N = this->mesh_num.GetNbDof();
  vec_Mh.Reallocate(this->ref_epsilon.GetM());
  vec_Mr2h.Reallocate(this->ref_epsilon.GetM());
  vec_Kh.Reallocate(this->ref_epsilon.GetM());
  vec_Kr2h.Reallocate(this->ref_epsilon.GetM());
  vec_Ch.Reallocate(this->ref_epsilon.GetM());
  vec_Cr2h.Reallocate(this->ref_epsilon.GetM());
  vec_Rjh.Reallocate(this->ref_epsilon.GetM());
  vec_Rjr2h.Reallocate(this->ref_epsilon.GetM());
  vec_Rih.Reallocate(this->ref_epsilon.GetM());
  vec_Rir2h.Reallocate(this->ref_epsilon.GetM());
  vec_Sjh.Reallocate(this->ref_epsilon.GetM());
  vec_Sjr2h.Reallocate(this->ref_epsilon.GetM());
  vec_Sih.Reallocate(this->ref_epsilon.GetM());
  vec_Sir2h.Reallocate(this->ref_epsilon.GetM());
  vec_Dh.Reallocate(this->ref_epsilon.GetM());
  vec_Dr2h.Reallocate(this->ref_epsilon.GetM());
  vec_Pos.Reallocate(N); vec_Pos.Fill(-100.0);

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  for (int ref = 0; ref < index_used.GetM(); ref++)
    if (index_used(ref))
      {
        vec_Mh(ref).Reallocate(N, N); vec_Mh(ref).Zero();
        vec_Mr2h(ref).Reallocate(N, N); vec_Mr2h(ref).Zero();
        vec_Kh(ref).Reallocate(N, N); vec_Kh(ref).Zero();
        vec_Kr2h(ref).Reallocate(N, N); vec_Kr2h(ref).Zero();
        vec_Ch(ref).Reallocate(N, N); vec_Ch(ref).Zero();
        vec_Cr2h(ref).Reallocate(N, N); vec_Cr2h(ref).Zero();
        vec_Rjh(ref).Reallocate(N, N); vec_Rjh(ref).Zero();
        vec_Rjr2h(ref).Reallocate(N, N); vec_Rjr2h(ref).Zero();
        vec_Rih(ref).Reallocate(N, N); vec_Rih(ref).Zero();
        vec_Rir2h(ref).Reallocate(N, N); vec_Rir2h(ref).Zero();
        vec_Sjh(ref).Reallocate(N, N); vec_Sjh(ref).Zero();
        vec_Sjr2h(ref).Reallocate(N, N); vec_Sjr2h(ref).Zero();
        vec_Sih(ref).Reallocate(N, N); vec_Sih(ref).Zero();
        vec_Sir2h(ref).Reallocate(N, N); vec_Sir2h(ref).Zero();
        vec_Dh(ref).Reallocate(N, N); vec_Dh(ref).Zero();
        vec_Dr2h(ref).Reallocate(N, N); vec_Dr2h(ref).Zero();
      }

  VectR2 s;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;


  Real_wp threshold = 1e-15;

  // boucle sur les elements
  VectReal_wp Bzero, Bmass, BmassR2; Matrix2_2 dfjm1, Ctmp, Cref, Dtmp, Dref;
  VectR2 Radv, Sadv, RRadv, SRadv, Ezero; Vector<Matrix2_2> Kstiff, KRstiff, Cstiff, CRstiff, Dstiff, DRstiff, Czero;
  Cref(0, 1) = -1.0; Cref(1, 0) = 1.0; // pour "retourner" les dérivées
  TinyVector<bool, 4> null_stiff(true, false, true, true);
  TinyVector<bool, 4> null_advi(true, true, false, true);
  TinyVector<bool, 4> null_advj(true, true, true, false);
  Matrix<Real_wp> mat_elem_Mh, mat_elem_Mr2h, mat_elem_Kh, mat_elem_Kr2h, mat_elem_Ch, mat_elem_Cr2h, mat_elem_Rjh,
                  mat_elem_Rjr2h, mat_elem_Rih, mat_elem_Rir2h, mat_elem_Sjh, mat_elem_Sjr2h, mat_elem_Sih, mat_elem_Sir2h, mat_elem_Dh, mat_elem_Dr2h;
  DISP(this->mesh.GetNbElt());
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(i);

      int ref = this->mesh.Element(i).GetReference();

      // calcul des matrices jacobiennes sur les points de quadrature
      this->mesh.GetVerticesElement(i, s);
      Fb.FjElem(s, PointsElem, this->mesh, i);
      Fb.DFjElem(s, PointsElem, MatricesElem, this->mesh, i);

      int Nquad = Fb.GetNbPointsQuadratureInside();
      Bmass.Reallocate(Nquad);
      BmassR2.Reallocate(Nquad);
      Radv.Reallocate(Nquad);
      Sadv.Reallocate(Nquad);
      Kstiff.Reallocate(Nquad);
      Cstiff.Reallocate(Nquad);
      Dstiff.Reallocate(Nquad);
      RRadv.Reallocate(Nquad);
      SRadv.Reallocate(Nquad);
      KRstiff.Reallocate(Nquad);
      CRstiff.Reallocate(Nquad);
      DRstiff.Reallocate(Nquad);

      // on recupere le numero de conducteur pour chaque ddl
      VectReal_wp phi;
      int nb_dof = Fb.GetNbDof();
      IVect num_ddl = this->mesh_num.Element(i).GetNodle();

      // boucle sur les points de quadrature
      for (int j = 0; j < Nquad; j++)
        {
          // pour le calcul de Mh (matrice de masse)
          Real_wp weight = Fb.WeightsND(j);
          Real_wp jacob = Det(MatricesElem.GetPointQuadrature(j));
          Real_wp x_quadGlob = PointsElem.GetPointQuadrature(j)(0);
          Real_wp y_quadGlob = PointsElem.GetPointQuadrature(j)(1);
          R2 PRef(x_quadGlob, y_quadGlob);
          R2 VRef(-y_quadGlob, x_quadGlob);
          Real_wp r2 = x_quadGlob*x_quadGlob + y_quadGlob*y_quadGlob;
          // R2 TRef(0.0, 0.0);
          // R2 URef(0.0, 0.0);
          Bmass(j) = jacob*weight;
          BmassR2(j) = jacob*weight*r2;

          // pour le calcul de Kh (matrice de rigidite)
          GetInverse(MatricesElem.GetPointQuadrature(j), dfjm1);

          MltTrans(dfjm1, dfjm1, Kstiff(j));
          Mlt(jacob*weight, Kstiff(j));

          MltTrans(dfjm1, dfjm1, KRstiff(j));
          Mlt(jacob*weight*r2, KRstiff(j));

          // pour le calcul de Ch (matrice de couplage)
          MltTrans(Cref, dfjm1, Ctmp);
          Mlt(dfjm1, Ctmp, Cstiff(j));
          Mlt(jacob*weight, Cstiff(j));

          Mlt(dfjm1, Ctmp, CRstiff(j));
          Mlt(jacob*weight*r2, CRstiff(j));

          // pour le calcul de Dh
          // calcul de l'application sur les gradients
          Dref(0, 0) = -y_quadGlob * y_quadGlob;
          Dref(1, 0) = x_quadGlob * y_quadGlob;
          Dref(0, 1) = x_quadGlob * y_quadGlob;
          Dref(1, 1) = -x_quadGlob * x_quadGlob;

          MltTrans(Dref, dfjm1, Dtmp);
          Mlt(dfjm1, Dtmp, Dstiff(j));
          Mlt(jacob*weight, Dstiff(j));

          Mlt(dfjm1, Dtmp, DRstiff(j));
          Mlt(jacob*weight*r2, DRstiff(j));

          // (xi scalaire grad phi_i) phi_j
          Mlt(dfjm1, PRef, Radv(j));
          Mlt(jacob*weight, Radv(j));

          Mlt(dfjm1, PRef, RRadv(j));
          Mlt(jacob*weight*r2, RRadv(j));

          // (xi vectoriel grad phi_i) phi_j
          Mlt(dfjm1, VRef, Sadv(j));
          Mlt(jacob*weight, Sadv(j));

          Mlt(dfjm1, VRef, SRadv(j));
          Mlt(jacob*weight*r2, SRadv(j));
        }

      if (Fb.GetNbPointsDof() != Fb.GetNbDof())
        {
          cout << "Element fini pas nodal " << endl;
          abort();
        }

      // boucle sur les ddls (supposes nodaux)
      for (int k = 0; k < Fb.GetNbPointsDof(); k++)
        {
          if (num_ddl(k) >= 0)
            {
              vec_Pos(num_ddl(k))(0) = PointsElem.GetPointDof(k)(0);
              vec_Pos(num_ddl(k))(1) = PointsElem.GetPointDof(k)(1);
            }
        }

      // calcul des matrices de masse et rigidite elementaires
      mat_elem_Mh.Reallocate(nb_dof, nb_dof); mat_elem_Mh.Zero();
      mat_elem_Mr2h.Reallocate(nb_dof, nb_dof); mat_elem_Mr2h.Zero();
      mat_elem_Kh.Reallocate(nb_dof, nb_dof); mat_elem_Kh.Zero();
      mat_elem_Kr2h.Reallocate(nb_dof, nb_dof); mat_elem_Kr2h.Zero();
      mat_elem_Ch.Reallocate(nb_dof, nb_dof); mat_elem_Ch.Zero();
      mat_elem_Cr2h.Reallocate(nb_dof, nb_dof); mat_elem_Cr2h.Zero();
      mat_elem_Rjh.Reallocate(nb_dof, nb_dof); mat_elem_Rjh.Zero();
      mat_elem_Rjr2h.Reallocate(nb_dof, nb_dof); mat_elem_Rjr2h.Zero();
      mat_elem_Rih.Reallocate(nb_dof, nb_dof); mat_elem_Rih.Zero();
      mat_elem_Rir2h.Reallocate(nb_dof, nb_dof); mat_elem_Rir2h.Zero();
      mat_elem_Sjh.Reallocate(nb_dof, nb_dof); mat_elem_Sjh.Zero();
      mat_elem_Sjr2h.Reallocate(nb_dof, nb_dof); mat_elem_Sjr2h.Zero();
      mat_elem_Sih.Reallocate(nb_dof, nb_dof); mat_elem_Sih.Zero();
      mat_elem_Sir2h.Reallocate(nb_dof, nb_dof); mat_elem_Sir2h.Zero();
      mat_elem_Dh.Reallocate(nb_dof, nb_dof); mat_elem_Dh.Zero();
      mat_elem_Dr2h.Reallocate(nb_dof, nb_dof); mat_elem_Dr2h.Zero();


      Fb.AddVariableMassMatrix(0, 0, Bmass, mat_elem_Mh);
      Fb.AddVariableMassMatrix(0, 0, BmassR2, mat_elem_Mr2h);

      Fb.AddVariableElemMatrix(0, 0, Bzero, Kstiff, Ezero, Ezero, null_stiff, mat_elem_Kh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, KRstiff, Ezero, Ezero, null_stiff, mat_elem_Kr2h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Cstiff, Ezero, Ezero, null_stiff, mat_elem_Ch);
      Fb.AddVariableElemMatrix(0, 0, Bzero, CRstiff, Ezero, Ezero, null_stiff, mat_elem_Cr2h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Dstiff, Ezero, Ezero, null_stiff, mat_elem_Dh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, DRstiff, Ezero, Ezero, null_stiff, mat_elem_Dr2h);

      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Radv, Ezero, null_advi, mat_elem_Rih);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, RRadv, Ezero, null_advi, mat_elem_Rir2h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Sadv, Ezero, null_advi, mat_elem_Sih);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, SRadv, Ezero, null_advi, mat_elem_Sir2h);

      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Ezero, Radv, null_advj, mat_elem_Rjh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Ezero, RRadv, null_advj, mat_elem_Rjr2h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Ezero, Sadv, null_advj, mat_elem_Sjh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Ezero, SRadv, null_advj, mat_elem_Sjr2h);


      // assemblage des matrices
      IVect permut(nb_dof); permut.Fill();
      Sort(num_ddl, permut);

      Vector<int> col(nb_dof); Vector<Complexe> val(nb_dof); int nb_val = 0;
      col.Fill(-1); val.Zero();
      for (int j = 0; j < nb_dof; j++)
        if (num_ddl(j) >= 0)
          {
            nb_val = 0;

            //DISP(i); DISP(j); DISP(coef1); DISP(num_ddl(j));
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Mh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Mh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Mh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Mr2h(permut(j), permut(k))) >= threshold))
                {

                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Mr2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Mr2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Rjh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Rjh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Rjh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Rjr2h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Rjr2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Rjr2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Rih(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Rih(permut(j), permut(k));
                  nb_val++;
                }

            vec_Rih(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Rir2h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Rir2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Rir2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Sih(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Sih(permut(j), permut(k));
                  nb_val++;
                }

            vec_Sih(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Sir2h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Sir2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Sir2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Sjh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Sjh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Sjh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Sjr2h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Sjr2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Sjr2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Kh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Kh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Kh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Kr2h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Kr2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Kr2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Ch(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Ch(permut(j), permut(k));
                  nb_val++;
                }

            vec_Ch(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Cr2h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Cr2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Cr2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Dh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Dh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Dh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Dr2h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Dr2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Dr2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);
          }
    }
}


void ModeEs_Solver::ComputeFemMatricesSurf(int ref_surf,
                        VectReal_wp& Mh_surf,
                        Matrix<Complexe, Prop, Storage>& Kh_surf,
                        Matrix<Complexe, Prop, Storage>& Rh_surf)
{
  // on alloue les matrices
  int N = this->mesh_num.GetNbDof();

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  VectR2 s;
  Matrix2_2 dfjm1;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;



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


void ModeEs_Solver::ComputeFemMatricesLin(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                        Vector<Matrix<Complex_wp, Prop, Storage> >& vec_M1h,
                        Vector<Matrix<Complex_wp, Prop, Storage> >& vec_M2h,
                        Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Ti1h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Ti2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Tj1h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Tj2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Di1h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Di2h,
                        VectR2& vec_Pos
                      )
{
  // on alloue les matrices
  int N = this->mesh_num.GetNbDof();
  vec_Mh.Reallocate(this->ref_epsilon.GetM());
  vec_M1h.Reallocate(this->ref_epsilon.GetM());
  vec_M2h.Reallocate(this->ref_epsilon.GetM());
  vec_Mr2h.Reallocate(this->ref_epsilon.GetM());
  vec_Rih.Reallocate(this->ref_epsilon.GetM());
  vec_Sih.Reallocate(this->ref_epsilon.GetM());
  vec_Ti1h.Reallocate(this->ref_epsilon.GetM());
  vec_Ti2h.Reallocate(this->ref_epsilon.GetM());
  vec_Tj1h.Reallocate(this->ref_epsilon.GetM());
  vec_Tj2h.Reallocate(this->ref_epsilon.GetM());
  vec_Di1h.Reallocate(this->ref_epsilon.GetM());
  vec_Di2h.Reallocate(this->ref_epsilon.GetM());
  vec_Pos.Reallocate(N); vec_Pos.Fill(-100.0);

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  for (int ref = 0; ref < index_used.GetM(); ref++)
    if (index_used(ref))
      {
        vec_Mh(ref).Reallocate(N, N);
        vec_M1h(ref).Reallocate(N, N);
        vec_M2h(ref).Reallocate(N, N);
        vec_Mr2h(ref).Reallocate(N, N);
        vec_Rih(ref).Reallocate(N, N);
        vec_Sih(ref).Reallocate(N, N);
        vec_Ti1h(ref).Reallocate(N, N);
        vec_Ti2h(ref).Reallocate(N, N);
        vec_Tj1h(ref).Reallocate(N, N);
        vec_Tj2h(ref).Reallocate(N, N);
        vec_Di1h(ref).Reallocate(N, N);
        vec_Di2h(ref).Reallocate(N, N);
      }

  VectR2 s;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;


  Real_wp threshold = 1e-15;

  // boucle sur les elements
  VectReal_wp Bzero, Bmass, Bmass1, Bmass2, BmassR2; Matrix2_2 dfjm1;
  VectR2 Dzero, Radv, Sadv, T1adv, T2adv, D1adv, D2adv, Ezero; Vector<Matrix2_2> Czero;
  // TinyVector<bool, 4> null_stiff(true, false, true, true);
  TinyVector<bool, 4> null_advi(true, true, false, true);
  TinyVector<bool, 4> null_advj(true, true, true, false);
  Matrix<Real_wp> mat_elem_Mh, mat_elem_M1h, mat_elem_M2h, mat_elem_Mr2h, mat_elem_Rh, mat_elem_Sh, mat_elem_Ti1h, mat_elem_Ti2h, mat_elem_Tj1h, mat_elem_Tj2h, mat_elem_D1h, mat_elem_D2h;
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(i);

      int ref = this->mesh.Element(i).GetReference();

      // calcul des matrices jacobiennes sur les points de quadrature
      this->mesh.GetVerticesElement(i, s);
      Fb.FjElem(s, PointsElem, this->mesh, i);
      Fb.DFjElem(s, PointsElem, MatricesElem, this->mesh, i);

      int Nquad = Fb.GetNbPointsQuadratureInside();
      Bmass.Reallocate(Nquad);
      Bmass1.Reallocate(Nquad);
      Bmass2.Reallocate(Nquad);
      BmassR2.Reallocate(Nquad);
      Radv.Reallocate(Nquad);
      Sadv.Reallocate(Nquad);
      T1adv.Reallocate(Nquad);
      T2adv.Reallocate(Nquad);
      D1adv.Reallocate(Nquad);
      D2adv.Reallocate(Nquad);

      // on recupere le numero de conducteur pour chaque ddl
      VectReal_wp phi;
      int nb_dof = Fb.GetNbDof();
      IVect num_ddl = this->mesh_num.Element(i).GetNodle();

      // boucle sur les points de quadrature
      for (int j = 0; j < Nquad; j++)
        {
          // pour le calcul de Mh (matrice de masse)
          Real_wp weight = Fb.WeightsND(j);
          Real_wp jacob = Det(MatricesElem.GetPointQuadrature(j));
          Real_wp x_quadGlob = PointsElem.GetPointQuadrature(j)(0);
          Real_wp y_quadGlob = PointsElem.GetPointQuadrature(j)(1);
          R2 PRef(x_quadGlob, y_quadGlob);
          R2 VRef(-y_quadGlob, x_quadGlob);
          R2 T1(1.0, 0.0);
          R2 T2(0.0, 1.0);
          // R2 TRef(0.0, 0.0);
          // R2 URef(0.0, 0.0);
          Bmass(j) = jacob*weight;
          Bmass1(j) = jacob*weight*x_quadGlob;
          Bmass2(j) = jacob*weight*y_quadGlob;
          BmassR2(j) = jacob*weight*(x_quadGlob*x_quadGlob + y_quadGlob*y_quadGlob);

          // pour le calcul de Kh (matrice de rigidite)
          GetInverse(MatricesElem.GetPointQuadrature(j), dfjm1);

          // (xi scalaire grad phi_i) phi_j
          Mlt(dfjm1, PRef, Radv(j));
          Mlt(jacob*weight, Radv(j));

          // (xi vectoriel grad phi_i) phi_j
          Mlt(dfjm1, VRef, Sadv(j));
          Mlt(jacob*weight, Sadv(j));

          // (dx1 phi_i) phi_j
          Mlt(dfjm1, T1, T1adv(j));
          Mlt(jacob*weight, T1adv(j));

          // (dx2 phi_i) phi_j
          Mlt(dfjm1, T2, T2adv(j));
          Mlt(jacob*weight, T2adv(j));

          // x1 (xi vectoriel grad phi_i) phi_j
          Mlt(dfjm1, VRef, D1adv(j));
          Mlt(jacob*weight*x_quadGlob, D1adv(j));

          // x2 (xi vectoriel grad phi_i) phi_j
          Mlt(dfjm1, VRef, D2adv(j));
          Mlt(jacob*weight*y_quadGlob, D2adv(j));

        }

      if (Fb.GetNbPointsDof() != Fb.GetNbDof())
        {
          cout << "Element fini pas nodal " << endl;
          abort();
        }

      // boucle sur les ddls (supposes nodaux)
      for (int k = 0; k < Fb.GetNbPointsDof(); k++)
        {
          if (num_ddl(k) >= 0)
            {
              vec_Pos(num_ddl(k))(0) = PointsElem.GetPointDof(k)(0);
              vec_Pos(num_ddl(k))(1) = PointsElem.GetPointDof(k)(1);
            }
        }

      // calcul des matrices de masse et rigidite elementaires
      mat_elem_Mh.Reallocate(nb_dof, nb_dof); mat_elem_Mh.Zero();
      mat_elem_M1h.Reallocate(nb_dof, nb_dof); mat_elem_M1h.Zero();
      mat_elem_M2h.Reallocate(nb_dof, nb_dof); mat_elem_M2h.Zero();
      mat_elem_Mr2h.Reallocate(nb_dof, nb_dof); mat_elem_Mr2h.Zero();
      mat_elem_Rh.Reallocate(nb_dof, nb_dof); mat_elem_Rh.Zero();
      mat_elem_Sh.Reallocate(nb_dof, nb_dof); mat_elem_Sh.Zero();
      mat_elem_Ti1h.Reallocate(nb_dof, nb_dof); mat_elem_Ti1h.Zero();
      mat_elem_Ti2h.Reallocate(nb_dof, nb_dof); mat_elem_Ti2h.Zero();
      mat_elem_Tj1h.Reallocate(nb_dof, nb_dof); mat_elem_Tj1h.Zero();
      mat_elem_Tj2h.Reallocate(nb_dof, nb_dof); mat_elem_Tj2h.Zero();
      mat_elem_D1h.Reallocate(nb_dof, nb_dof); mat_elem_D1h.Zero();
      mat_elem_D2h.Reallocate(nb_dof, nb_dof); mat_elem_D2h.Zero();

      Fb.AddVariableMassMatrix(0, 0, Bmass, mat_elem_Mh);
      Fb.AddVariableMassMatrix(0, 0, Bmass1, mat_elem_M1h);
      Fb.AddVariableMassMatrix(0, 0, Bmass2, mat_elem_M2h);
      Fb.AddVariableMassMatrix(0, 0, BmassR2, mat_elem_Mr2h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Radv, Ezero, null_advi, mat_elem_Rh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Sadv, Ezero, null_advi, mat_elem_Sh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, T1adv, Ezero, null_advi, mat_elem_Ti1h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, T2adv, Ezero, null_advi, mat_elem_Ti2h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Ezero, T1adv, null_advj, mat_elem_Tj1h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Ezero, T2adv, null_advj, mat_elem_Tj2h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, D1adv, Ezero, null_advi, mat_elem_D1h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, D2adv, Ezero, null_advi, mat_elem_D2h);

      // assemblage des matrices
      IVect permut(nb_dof); permut.Fill();
      Sort(num_ddl, permut);

      Vector<int> col(nb_dof); Vector<Complexe> val(nb_dof); int nb_val = 0;
      col.Fill(-1); val.Zero();
      for (int j = 0; j < nb_dof; j++)
        if (num_ddl(j) >= 0)
          {
            nb_val = 0;

            //DISP(i); DISP(j); DISP(coef1); DISP(num_ddl(j));
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Mh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Mh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Mh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_M1h(permut(j), permut(k))) >= threshold))
                {

                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_M1h(permut(j), permut(k));
                  nb_val++;
                }

            vec_M1h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_M2h(permut(j), permut(k))) >= threshold))
                {

                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_M2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_M2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Mr2h(permut(j), permut(k))) >= threshold))
                {

                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Mr2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Mr2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Rh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Rh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Rih(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Sh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Sh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Sih(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Ti1h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Ti1h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Ti1h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Ti2h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Ti2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Ti2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Tj1h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Tj1h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Tj1h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Tj2h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Tj2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Tj2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_D1h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_D1h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Di1h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_D2h(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_D2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Di2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);
          }
    }
}


void ModeEs_Solver::ComputeFemMatricesPoly2(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                        Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjh,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjh,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                        Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                        VectR2& vec_Pos
                      )
{
  // on alloue les matrices
  int N = this->mesh_num.GetNbDof();
  vec_Mh.Reallocate(this->ref_epsilon.GetM());
  vec_Mr2h.Reallocate(this->ref_epsilon.GetM());
  vec_Kh.Reallocate(this->ref_epsilon.GetM());
  vec_Ch.Reallocate(this->ref_epsilon.GetM());
  vec_Rjh.Reallocate(this->ref_epsilon.GetM());
  vec_Rih.Reallocate(this->ref_epsilon.GetM());
  vec_Sjh.Reallocate(this->ref_epsilon.GetM());
  vec_Sih.Reallocate(this->ref_epsilon.GetM());
  vec_Dh.Reallocate(this->ref_epsilon.GetM());
  vec_Pos.Reallocate(N); vec_Pos.Fill(-100.0);

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  for (int ref = 0; ref < index_used.GetM(); ref++)
    if (index_used(ref))
      {
        vec_Mh(ref).Reallocate(N, N); vec_Mh(ref).Zero();
        vec_Mr2h(ref).Reallocate(N, N); vec_Mr2h(ref).Zero();
        vec_Kh(ref).Reallocate(N, N); vec_Kh(ref).Zero();
        vec_Ch(ref).Reallocate(N, N); vec_Ch(ref).Zero();
        vec_Rjh(ref).Reallocate(N, N); vec_Rjh(ref).Zero();
        vec_Rih(ref).Reallocate(N, N); vec_Rih(ref).Zero();
        vec_Sjh(ref).Reallocate(N, N); vec_Sjh(ref).Zero();
        vec_Sih(ref).Reallocate(N, N); vec_Sih(ref).Zero();
        vec_Dh(ref).Reallocate(N, N); vec_Dh(ref).Zero();
      }

  VectR2 s;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;


  Real_wp threshold = 1e-15;

  // boucle sur les elements
  VectReal_wp Bzero, Bmass, BmassR2; Matrix2_2 dfjm1, Ctmp, Cref, Dtmp, Dref;
  VectR2 Radv, Sadv, Ezero; Vector<Matrix2_2> Kstiff, Cstiff, Dstiff, Czero;
  Cref(0, 1) = -1.0; Cref(1, 0) = 1.0; // pour "retourner" les dérivées
  TinyVector<bool, 4> null_stiff(true, false, true, true);
  TinyVector<bool, 4> null_advi(true, true, false, true);
  TinyVector<bool, 4> null_advj(true, true, true, false);
  Matrix<Real_wp> mat_elem_Mh, mat_elem_Mr2h, mat_elem_Kh, mat_elem_Ch, mat_elem_Rjh, mat_elem_Rih, mat_elem_Sjh, mat_elem_Sih, mat_elem_Dh;
  DISP(this->mesh.GetNbElt());
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(i);

      int ref = this->mesh.Element(i).GetReference();

      // calcul des matrices jacobiennes sur les points de quadrature
      this->mesh.GetVerticesElement(i, s);
      Fb.FjElem(s, PointsElem, this->mesh, i);
      Fb.DFjElem(s, PointsElem, MatricesElem, this->mesh, i);

      int Nquad = Fb.GetNbPointsQuadratureInside();
      Bmass.Reallocate(Nquad);
      BmassR2.Reallocate(Nquad);
      Radv.Reallocate(Nquad);
      Sadv.Reallocate(Nquad);
      Kstiff.Reallocate(Nquad);
      Cstiff.Reallocate(Nquad);
      Dstiff.Reallocate(Nquad);

      // on recupere le numero de conducteur pour chaque ddl
      VectReal_wp phi;
      int nb_dof = Fb.GetNbDof();
      IVect num_ddl = this->mesh_num.Element(i).GetNodle();

      // boucle sur les points de quadrature
      for (int j = 0; j < Nquad; j++)
        {
          // pour le calcul de Mh (matrice de masse)
          Real_wp weight = Fb.WeightsND(j);
          Real_wp jacob = Det(MatricesElem.GetPointQuadrature(j));
          Real_wp x_quadGlob = PointsElem.GetPointQuadrature(j)(0);
          Real_wp y_quadGlob = PointsElem.GetPointQuadrature(j)(1);
          R2 PRef(x_quadGlob, y_quadGlob);
          R2 VRef(-y_quadGlob, x_quadGlob);
          // R2 TRef(0.0, 0.0);
          // R2 URef(0.0, 0.0);
          Bmass(j) = jacob*weight;
          BmassR2(j) = jacob*weight*(x_quadGlob*x_quadGlob + y_quadGlob*y_quadGlob);

          // pour le calcul de Kh (matrice de rigidite)
          GetInverse(MatricesElem.GetPointQuadrature(j), dfjm1);

          MltTrans(dfjm1, dfjm1, Kstiff(j));
          Mlt(jacob*weight, Kstiff(j));

          // pour le calcul de Ch (matrice de couplage)
          MltTrans(Cref, dfjm1, Ctmp);
          Mlt(dfjm1, Ctmp, Cstiff(j));
          Mlt(jacob*weight, Cstiff(j));


          // pour le calcul de Dh
          // calcul de l'application sur les gradients
          Dref(0, 0) = -y_quadGlob * y_quadGlob;
          Dref(1, 0) = x_quadGlob * y_quadGlob;
          Dref(0, 1) = x_quadGlob * y_quadGlob;
          Dref(1, 1) = -x_quadGlob * x_quadGlob;

          MltTrans(Dref, dfjm1, Dtmp);
          Mlt(dfjm1, Dtmp, Dstiff(j));
          Mlt(jacob*weight, Dstiff(j));

          // (xi scalaire grad phi_i) phi_j
          Mlt(dfjm1, PRef, Radv(j));
          Mlt(jacob*weight, Radv(j));

          // (xi vectoriel grad phi_i) phi_j
          Mlt(dfjm1, VRef, Sadv(j));
          Mlt(jacob*weight, Sadv(j));

        }

      if (Fb.GetNbPointsDof() != Fb.GetNbDof())
        {
          cout << "Element fini pas nodal " << endl;
          abort();
        }

      // boucle sur les ddls (supposes nodaux)
      for (int k = 0; k < Fb.GetNbPointsDof(); k++)
        {
          if (num_ddl(k) >= 0)
            {
              vec_Pos(num_ddl(k))(0) = PointsElem.GetPointDof(k)(0);
              vec_Pos(num_ddl(k))(1) = PointsElem.GetPointDof(k)(1);
            }
        }

      // calcul des matrices de masse et rigidite elementaires
      mat_elem_Mh.Reallocate(nb_dof, nb_dof); mat_elem_Mh.Zero();
      mat_elem_Mr2h.Reallocate(nb_dof, nb_dof); mat_elem_Mr2h.Zero();
      mat_elem_Kh.Reallocate(nb_dof, nb_dof); mat_elem_Kh.Zero();
      mat_elem_Ch.Reallocate(nb_dof, nb_dof); mat_elem_Ch.Zero();
      mat_elem_Rjh.Reallocate(nb_dof, nb_dof); mat_elem_Rjh.Zero();
      mat_elem_Rih.Reallocate(nb_dof, nb_dof); mat_elem_Rih.Zero();
      mat_elem_Sjh.Reallocate(nb_dof, nb_dof); mat_elem_Sjh.Zero();
      mat_elem_Sih.Reallocate(nb_dof, nb_dof); mat_elem_Sih.Zero();
      mat_elem_Dh.Reallocate(nb_dof, nb_dof); mat_elem_Dh.Zero();

      Fb.AddVariableMassMatrix(0, 0, Bmass, mat_elem_Mh);
      Fb.AddVariableMassMatrix(0, 0, BmassR2, mat_elem_Mr2h);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Kstiff, Ezero, Ezero, null_stiff, mat_elem_Kh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Cstiff, Ezero, Ezero, null_stiff, mat_elem_Ch);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Dstiff, Ezero, Ezero, null_stiff, mat_elem_Dh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Radv, Ezero, null_advi, mat_elem_Rih);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Sadv, Ezero, null_advi, mat_elem_Sih);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Ezero, Radv, null_advj, mat_elem_Rjh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Ezero, Sadv, null_advj, mat_elem_Sjh);

      // assemblage des matrices
      IVect permut(nb_dof); permut.Fill();
      Sort(num_ddl, permut);

      Vector<int> col(nb_dof); Vector<Complexe> val(nb_dof); int nb_val = 0;
      col.Fill(-1); val.Zero();
      for (int j = 0; j < nb_dof; j++)
        if (num_ddl(j) >= 0)
          {
            nb_val = 0;

            //DISP(i); DISP(j); DISP(coef1); DISP(num_ddl(j));
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Mh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Mh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Mh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Mr2h(permut(j), permut(k))) >= threshold))
                {

                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Mr2h(permut(j), permut(k));
                  nb_val++;
                }

            vec_Mr2h(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Rjh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Rjh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Rjh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Rih(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Rih(permut(j), permut(k));
                  nb_val++;
                }

            vec_Rih(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Sih(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Sih(permut(j), permut(k));
                  nb_val++;
                }

            vec_Sih(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Sjh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Sjh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Sjh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Kh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Kh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Kh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Ch(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Ch(permut(j), permut(k));
                  nb_val++;
                }

            vec_Ch(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Dh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Dh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Dh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);
          }
    }
}


void ModeEs_Solver::ComputeFemMatricesNonLin2N(const Complex_wp& L,
                                              Vector<Matrix<Complex_wp, Prop, Storage> >& vec_MhB,
                                              Vector<Matrix<Complexe, Prop, Storage> >& vec_KhB,
                                              Vector<Matrix<Complexe, Prop, Storage> >& vec_ChB,
                                              Vector<Matrix<Complexe, Prop, Storage> >& vec_RjhB,
                                              Vector<Matrix<Complexe, Prop, Storage> >& vec_RihB,
                                              Vector<Matrix<Complexe, Prop, Storage> >& vec_SjhB,
                                              Vector<Matrix<Complexe, Prop, Storage> >& vec_SihB,
                                              Vector<Matrix<Complexe, Prop, Storage> >& vec_DhB
                                            )
{
  if (isinf(abs(L)) || isnan(abs(L)))
    {
      DISP(L);
      abort();
    }
  Complex_wp beta = L*this->omega;

  // on alloue les matrices
  int N = this->mesh_num.GetNbDof();
  vec_MhB.Reallocate(this->ref_epsilon.GetM());
  vec_KhB.Reallocate(this->ref_epsilon.GetM());
  vec_ChB.Reallocate(this->ref_epsilon.GetM());
  vec_RjhB.Reallocate(this->ref_epsilon.GetM());
  vec_RihB.Reallocate(this->ref_epsilon.GetM());
  vec_SjhB.Reallocate(this->ref_epsilon.GetM());
  vec_SihB.Reallocate(this->ref_epsilon.GetM());
  vec_DhB.Reallocate(this->ref_epsilon.GetM());

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  for (int ref = 0; ref < index_used.GetM(); ref++)
    if (index_used(ref))
      {
        vec_MhB(ref).Reallocate(N, N); vec_MhB(ref).Zero();
        vec_KhB(ref).Reallocate(N, N); vec_KhB(ref).Zero();
        vec_ChB(ref).Reallocate(N, N); vec_ChB(ref).Zero();
        vec_RjhB(ref).Reallocate(N, N); vec_RjhB(ref).Zero();
        vec_RihB(ref).Reallocate(N, N); vec_RihB(ref).Zero();
        vec_SjhB(ref).Reallocate(N, N); vec_SjhB(ref).Zero();
        vec_SihB(ref).Reallocate(N, N); vec_SihB(ref).Zero();
        vec_DhB(ref).Reallocate(N, N); vec_DhB(ref).Zero();
      }

  VectR2 s;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;


  Real_wp threshold = 1e-15;

  // boucle sur les elements
  VectComplex_wp Bzero, Bmass; Matrix2_2 dfjm1, Ctmp, Cref, Dtmp, Dref;
  VectR2_Complex_wp Radv, Sadv, Ezero; Vector<Matrix2_2_Complex_wp> Kstiff, Cstiff, Dstiff, Czero;
  Cref(0, 1) = -1.0; Cref(1, 0) = 1.0; // pour "retourner" les dérivées
  TinyVector<bool, 4> null_stiff(true, false, true, true);
  TinyVector<bool, 4> null_advi(true, true, false, true);
  TinyVector<bool, 4> null_advj(true, true, true, false);
  Matrix<Complex_wp> mat_elem_Mh, mat_elem_Kh, mat_elem_Ch, mat_elem_Rjh, mat_elem_Rih, mat_elem_Sjh, mat_elem_Sih, mat_elem_Dh;
  // DISP(this->mesh.GetNbElt());
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(i);

      int ref = this->mesh.Element(i).GetReference();

      Complex_wp eps_sigma = -Iwp*this->omega*this->ref_epsilon(ref) + this->ref_sigma(ref);
      Complex_wp mu_tilde = -Iwp*this->omega*this->ref_mu(ref);
      Complex_wp k2 = -eps_sigma*mu_tilde;

      // calcul des matrices jacobiennes sur les points de quadrature
      this->mesh.GetVerticesElement(i, s);
      Fb.FjElem(s, PointsElem, this->mesh, i);
      Fb.DFjElem(s, PointsElem, MatricesElem, this->mesh, i);

      int Nquad = Fb.GetNbPointsQuadratureInside();
      Bmass.Reallocate(Nquad);
      Radv.Reallocate(Nquad);
      Sadv.Reallocate(Nquad);
      Kstiff.Reallocate(Nquad);
      Cstiff.Reallocate(Nquad);
      Dstiff.Reallocate(Nquad);

      // on recupere le numero de conducteur pour chaque ddl
      VectReal_wp phi;
      int nb_dof = Fb.GetNbDof();
      IVect num_ddl = this->mesh_num.Element(i).GetNodle();

      // boucle sur les points de quadrature
      for (int j = 0; j < Nquad; j++)
        {
          // pour le calcul de Mh (matrice de masse)
          Real_wp weight = Fb.WeightsND(j);
          Real_wp jacob = Det(MatricesElem.GetPointQuadrature(j));
          Real_wp x_quadGlob = PointsElem.GetPointQuadrature(j)(0);
          Real_wp y_quadGlob = PointsElem.GetPointQuadrature(j)(1);
          R2 PRef(x_quadGlob, y_quadGlob);
          R2 VRef(-y_quadGlob, x_quadGlob);

          Real_wp r2 = x_quadGlob * x_quadGlob + y_quadGlob * y_quadGlob;

          Complex_wp D_Beta = beta*beta - k2;
          Complex_wp D_tauBeta = beta*beta - k2 * (1 + this->torsion * this->torsion * r2);

          Complex_wp coefB = D_Beta / D_tauBeta;

          // R2 TRef(0.0, 0.0);
          // R2 URef(0.0, 0.0);
          Bmass(j) = jacob*weight * coefB;

          // pour le calcul de Kh (matrice de rigidite)
          GetInverse(MatricesElem.GetPointQuadrature(j), dfjm1);

          MltTrans(dfjm1, dfjm1, Kstiff(j));
          Mlt(jacob*weight * coefB, Kstiff(j));

          // pour le calcul de Ch (matrice de couplage)
          MltTrans(Cref, dfjm1, Ctmp);
          Mlt(dfjm1, Ctmp, Cstiff(j));
          Mlt(jacob*weight * coefB, Cstiff(j));


          // pour le calcul de Dh
          // calcul de l'application sur les gradients
          Dref(0, 0) = -y_quadGlob * y_quadGlob;
          Dref(1, 0) = x_quadGlob * y_quadGlob;
          Dref(0, 1) = x_quadGlob * y_quadGlob;
          Dref(1, 1) = -x_quadGlob * x_quadGlob;

          MltTrans(Dref, dfjm1, Dtmp);
          Mlt(dfjm1, Dtmp, Dstiff(j));
          Mlt(jacob*weight * coefB, Dstiff(j));

          // (xi scalaire grad phi_i) phi_j
          Mlt(dfjm1, PRef, Radv(j));
          Mlt(jacob*weight * coefB, Radv(j));

          // (xi vectoriel grad phi_i) phi_j
          Mlt(dfjm1, VRef, Sadv(j));
          Mlt(jacob*weight * coefB, Sadv(j));

        }

      if (Fb.GetNbPointsDof() != Fb.GetNbDof())
        {
          cout << "Element fini pas nodal " << endl;
          abort();
        }

      // calcul des matrices de masse et rigidite elementaires
      mat_elem_Mh.Reallocate(nb_dof, nb_dof); mat_elem_Mh.Zero();
      mat_elem_Kh.Reallocate(nb_dof, nb_dof); mat_elem_Kh.Zero();
      mat_elem_Ch.Reallocate(nb_dof, nb_dof); mat_elem_Ch.Zero();
      mat_elem_Rjh.Reallocate(nb_dof, nb_dof); mat_elem_Rjh.Zero();
      mat_elem_Rih.Reallocate(nb_dof, nb_dof); mat_elem_Rih.Zero();
      mat_elem_Sjh.Reallocate(nb_dof, nb_dof); mat_elem_Sjh.Zero();
      mat_elem_Sih.Reallocate(nb_dof, nb_dof); mat_elem_Sih.Zero();
      mat_elem_Dh.Reallocate(nb_dof, nb_dof); mat_elem_Dh.Zero();

      Fb.AddVariableMassMatrix(0, 0, Bmass, mat_elem_Mh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Kstiff, Ezero, Ezero, null_stiff, mat_elem_Kh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Cstiff, Ezero, Ezero, null_stiff, mat_elem_Ch);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Dstiff, Ezero, Ezero, null_stiff, mat_elem_Dh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Radv, Ezero, null_advi, mat_elem_Rih);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Sadv, Ezero, null_advi, mat_elem_Sih);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Ezero, Radv, null_advj, mat_elem_Rjh);
      Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Ezero, Sadv, null_advj, mat_elem_Sjh);

      // assemblage des matrices
      IVect permut(nb_dof); permut.Fill();
      Sort(num_ddl, permut);

      Vector<int> col(nb_dof); Vector<Complexe> val(nb_dof); int nb_val = 0;
      col.Fill(-1); val.Zero();
      for (int j = 0; j < nb_dof; j++)
        if (num_ddl(j) >= 0)
          {
            nb_val = 0;

            //DISP(i); DISP(j); DISP(coef1); DISP(num_ddl(j));
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Mh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Mh(permut(j), permut(k));
                  nb_val++;
                }

            vec_MhB(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Rjh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Rjh(permut(j), permut(k));
                  nb_val++;
                }

            vec_RjhB(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Rih(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Rih(permut(j), permut(k));
                  nb_val++;
                }

            vec_RihB(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Sih(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Sih(permut(j), permut(k));
                  nb_val++;
                }

            vec_SihB(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Sjh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Sjh(permut(j), permut(k));
                  nb_val++;
                }

            vec_SjhB(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Kh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Kh(permut(j), permut(k));
                  nb_val++;
                }

            vec_KhB(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Ch(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Ch(permut(j), permut(k));
                  nb_val++;
                }

            vec_ChB(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            nb_val = 0;
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Dh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Dh(permut(j), permut(k));
                  nb_val++;
                }

            vec_DhB(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);
          }
    }
}

void ModeEs_Solver::ComputeMassMatrix(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh)
{
  // on alloue les matrices
  int N = this->mesh_num.GetNbDof();
  vec_Mh.Reallocate(this->ref_epsilon.GetM());

  Vector<bool> index_used(this->ref_epsilon.GetM());
  index_used.Fill(false);
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    index_used(this->mesh.Element(i).GetReference()) = true;

  for (int ref = 0; ref < index_used.GetM(); ref++)
    if (index_used(ref))
      {
        vec_Mh(ref).Reallocate(N, N); vec_Mh(ref).Zero();
      }

  VectR2 s;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;


  Real_wp threshold = 1e-15;

  // boucle sur les elements
  VectComplex_wp Bzero, Bmass;
  Matrix<Complex_wp> mat_elem_Mh;
  // DISP(this->mesh.GetNbElt());
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(i);

      int ref = this->mesh.Element(i).GetReference();

      // calcul des matrices jacobiennes sur les points de quadrature
      this->mesh.GetVerticesElement(i, s);
      Fb.FjElem(s, PointsElem, this->mesh, i);
      Fb.DFjElem(s, PointsElem, MatricesElem, this->mesh, i);

      int Nquad = Fb.GetNbPointsQuadratureInside();
      Bmass.Reallocate(Nquad);

      // on recupere le numero de conducteur pour chaque ddl
      VectReal_wp phi;
      int nb_dof = Fb.GetNbDof();
      IVect num_ddl = this->mesh_num.Element(i).GetNodle();

      // boucle sur les points de quadrature
      for (int j = 0; j < Nquad; j++)
        {
          // pour le calcul de Mh (matrice de masse)
          Real_wp weight = Fb.WeightsND(j);
          Real_wp jacob = Det(MatricesElem.GetPointQuadrature(j));

          // R2 TRef(0.0, 0.0);
          // R2 URef(0.0, 0.0);
          Bmass(j) = jacob*weight;
        }

      if (Fb.GetNbPointsDof() != Fb.GetNbDof())
        {
          cout << "Element fini pas nodal " << endl;
          abort();
        }

      // calcul des matrices de masse et rigidite elementaires
      mat_elem_Mh.Reallocate(nb_dof, nb_dof); mat_elem_Mh.Zero();

      Fb.AddVariableMassMatrix(0, 0, Bmass, mat_elem_Mh);

      // assemblage des matrices
      IVect permut(nb_dof); permut.Fill();
      Sort(num_ddl, permut);

      Vector<int> col(nb_dof); Vector<Complexe> val(nb_dof); int nb_val = 0;
      col.Fill(-1); val.Zero();
      for (int j = 0; j < nb_dof; j++)
        if (num_ddl(j) >= 0)
          {
            nb_val = 0;

            //DISP(i); DISP(j); DISP(coef1); DISP(num_ddl(j));
            for (int k = 0; k < nb_dof; k++)
              if ((num_ddl(k) >= 0) && (abs(mat_elem_Mh(permut(j), permut(k))) >= threshold))
                {
                  col(nb_val) = num_ddl(k);
                  val(nb_val) = mat_elem_Mh(permut(j), permut(k));
                  nb_val++;
                }

            vec_Mh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);
          }
    }
}
