#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_HP_HCURL

#include "Hyperbolic/Aeroacoustic/MontjoieTimeAeroacoustic.hxx"
#include "Hyperbolic/Maxwell/MontjoieTimeMaxwell3D.hxx"

using namespace Montjoie;

template<class TypeEquation, class Vector1>
void ModifyVelocity(EllipticProblem<TypeEquation>& var, const Vector1& PointsQuadrature)
{
  typedef typename TypeEquation::Dimension Dimension;
  typedef typename Dimension::R_N R_N;
  
  int nb_elt = PointsQuadrature.GetM();
  var.eval_flow.Reallocate(nb_elt);  
  var.grad_flow.Reallocate(nb_elt);
  var.div_flow.Reallocate(nb_elt);
  var.eval_c0.Reallocate(nb_elt);
  var.grad_c0.Reallocate(nb_elt);
  R_N v0;
  // loop on quadrature points
  for (int i = 0; i < nb_elt; i++)
    {
      int Nquad = PointsQuadrature(i).GetM();
      var.eval_flow(i).Reallocate(Nquad);  
      var.grad_flow(i).Reallocate(Nquad);
      var.div_flow(i).Reallocate(Nquad);
      var.eval_c0(i).Reallocate(Nquad);
      var.grad_c0(i).Reallocate(Nquad);
      for (int j = 0; j < Nquad; j++)
	{
          Real_wp x = PointsQuadrature(i)(j)(0);
	  Real_wp y = PointsQuadrature(i)(j)(1);
          Real_wp vx = 0.5*cos(pi_wp*x/8.0);
          Real_wp vy = 0.5*cos(pi_wp*y/8.0);
          Real_wp dvx_dx = -0.5*pi_wp/8.0*sin(pi_wp*x/8.0);
          Real_wp dvy_dy = -0.5*pi_wp/8.0*sin(pi_wp*y/8.0);
          Real_wp dvx_dy = 0, dvy_dx = 0;
          //vx = 0; vy = 0; dvx_dx = 0; dvy_dy = 0;
          // flow
	  v0.Init(vx, vy);
	  var.eval_flow(i)(j) = v0;;
	  
          // and gradient of flow
          var.grad_flow(i)(j).Zero();
          var.grad_flow(i)(j)(0, 0) = dvx_dx;
          var.grad_flow(i)(j)(0, 1) = dvx_dy;
          var.grad_flow(i)(j)(1, 0) = dvy_dx;
          var.grad_flow(i)(j)(1, 1) = dvy_dy;
	  
	  var.eval_c0(i)(j) = 1.0;
	  var.grad_c0(i)(j).Zero();
	  var.div_flow(i)(j) = dvx_dx + dvy_dy;
	}
    }
  
}

template<class TypeEquation>
void RunAll(HyperbolicProblem<TypeEquation>& var, const string& input_file, const string& name_element)
{
  //const Mesh<Dimension>& mesh = var.var_harmonic.mesh;
  
  var.var_harmonic.InitIndices(100);
  
  ReadInputFile(input_file, var);
  
  var.var_harmonic.SetFirstOrderFormulation(true);
  var.var_harmonic.ComputeMeshAndFiniteElement(name_element);
  
  var.var_harmonic.PerformOtherInitializations();
  
  var.ComputeRightHandSide();
  
  //var.var_harmonic.ComputeMassMatrix(false, false);
  var.var_harmonic.ComputeMassMatrix();
  var.var_harmonic.ComputeQuasiPeriodicPhase();
  
  // on modifie l'indice du domaine 1
  //ModifyVelocity(var.var_harmonic, var.var_harmonic.Glob_PointsQuadrature);
  
  // on efface le tableau PointsQuadrature
  var.var_harmonic.Glob_PointsQuadrature.Clear();
  
  var.var_harmonic.mesh.Write("test.mesh");
  
  DistributedMatrix<Real_wp, General, ArrayRowSparse> Sh, Kh;
  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> Mh;
  VectReal_wp Dh;
  
  GlobalGenericMatrix<Real_wp> nat_mat;
  nat_mat.SetCoefMass(0);
  nat_mat.SetCoefStiffness(1.0);
  nat_mat.SetCoefDamping(0.0);
  
  var.var_harmonic.AddMatrixWithBC(Kh, nat_mat);
  
  nat_mat.SetCoefMass(1.0);
  nat_mat.SetCoefStiffness(0.0);
  nat_mat.SetCoefDamping(0.0);
  
  var.var_harmonic.AddMatrixWithBC(Mh, nat_mat);

  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefStiffness(0.0);
  nat_mat.SetCoefDamping(1.0);
  
  var.var_harmonic.AddMatrixWithBC(Sh, nat_mat);
  Mh.WriteText("Mh.dat");

  Kh.WriteText("Kh.dat");
  Sh.WriteText("Sh.dat");
  
  // pas de temps locaux
  /*VectReal_wp local_dt;
  ComputeLocalTimeStep(var, local_dt);
  
  IVect level(var.var_harmonic.mesh.GetNbElt());
  level.Fill(0);
  Real_wp dt_min = local_dt(0);
  for (int i = 0; i < level.GetM(); i++)
    {
      dt_min = min(dt_min, local_dt(i));
      if (local_dt(i) < var.deltat)
	level(i) = 1;
    }
  
  // ordre ?
  bool test_loop = true;
  int r = 0;
  while (test_loop)
    {
      if (var.deltat/(2.0*r+1.0) < dt_min)
	test_loop = false;
      else
	r++;
    }
  DISP(r);
  */
  
  var.InitTimeIterations();
  int N = var.var_harmonic.GetNbDof();
  Dh.Reallocate(N);
  for (int i = 0; i < N; i++)
    Dh(i) = Mh(i, i);
    
  // calcul de P
  /* typedef typename TypeEquation::TypeEquationStationary TypeEquationStationary;
  DistributedMatrix<Real_wp, General, ArrayRowSparse> P(N, N);
  int Nvol = mesh.GetNbDof();
  for (int i = 0; i < level.GetM(); i++)
    if (level(i) > 0)
      {
	int offset1 = mesh.Element(i).GetNumberDof(0);
	int offset2;
        if (i == mesh.GetNbElt()-1)
          offset2 = Nvol;
        else
          offset2 = mesh.Element(i+1).GetNumberDof(0);
	
	for (int j = offset1; j < offset2; j++)
	  for (int m = 0; m < TypeEquationStationary::nb_unknowns; m++)
	    P.AddInteraction(j + m*Nvol, j+m*Nvol, 1.0);
      }
  
  // et de KhP
  DistributedMatrix<Real_wp, General, ArrayRowSparse> KhP;
  Mlt(Kh, P, KhP);      
  */
  
  // calcul de vecteurs propres
  DistributedMatrix<Real_wp, General, ArrayRowSparse> Ah(N, N);
  Copy(Kh, Ah);
  Add(1.0, Sh, Ah);

  cout << "Tapez 1 pour utiliser la methode alternative de penalisation" << endl;
  int type_penal;
  cin >> type_penal;
  
  int nb_dof_u = N;
  
  if (type_penal == 1)
    {
      Ah.Reallocate(2*N, 2*N);
      Mh.Resize(2*N, 2*N);
      Dh.Resize(2*N);
      nb_dof_u = 2*N;
      for (int i = 0; i < N; i++)
        {
          int size_row_Kh = Kh.GetRowSize(i);
          int size_row_Sh = Sh.GetRowSize(i);
          
          Ah.ReallocateRow(i, size_row_Kh+size_row_Sh);
          for (int j = 0; j < size_row_Kh; j++)
            {
              Ah.Index(i, j) = Kh.Index(i, j);
              Ah.Value(i, j) = Kh.Value(i, j);
            }
          
          for (int j = 0; j < size_row_Sh; j++)
            {
              Ah.Index(i, size_row_Kh+j) = N + Sh.Index(i, j);
              Ah.Value(i, size_row_Kh+j) = Sh.Value(i, j);
            }
          
          Mh.ReallocateRow(i+N, 1);
          Ah.ReallocateRow(i+N, 1);
          Mh.Index(i+N, 0) = i+N;
          Mh.Value(i+N, 0) = 1.0;
          Ah.Index(i+N, 0) = i;
          Ah.Value(i+N, 0) =  -1.0;
          Dh(i+N) = 1.0;
        }
      
      //Mh.WriteText("Mh.dat");
      //Kh.WriteText("Kh.dat");
      //Sh.WriteText("Sh.dat");
      //Ah.WriteText("Ah.dat");
    }  
  
  cout << "Entrez le nombre de valeurs propres a calculer" <<endl;
  int nb_eigenval;
  cin >> nb_eigenval;
  
  if (nb_eigenval > 0)
    {
      SparseEigenProblem<double, DistributedMatrix<double, General, ArrayRowSparse> > var_eig;
      var_eig.SetStoppingCriterion(1e-12);
      var_eig.SetNbAskedEigenvalues(nb_eigenval);
      var_eig.SetPrintLevel(1);
      
      var_eig.InitMatrix(Ah, Mh);
      var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
      var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES,
                              complex<double>(0, var.var_harmonic.GetOmega()), var_eig.SORTED_MODULUS);
      
      VectReal_wp lambda, lambda_imag;
      Matrix<Real_wp> eigen_vec;
      GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
      
      DISP(lambda); DISP(lambda_imag);
      lambda.WriteText("eigenval_real.dat");
      lambda_imag.WriteText("eigenval_imag.dat");
      
      // the modes are written on files
      Vector<Real_wp> eigen_mode(eigen_vec.GetM());
      Vector<Vector<Real_wp> > eigen_mode_v;
      eigen_mode_v.SetData(1, &eigen_mode);
      for (int i = 0; i < nb_eigenval; i++)
        {
          for (int j = 0; j < eigen_vec.GetM(); j++)
            eigen_mode(j) = eigen_vec(j, i);
          
          string numero = NumberToString(i);
          for (int j = 0; j < var.var_harmonic.output_mesh_param.GetM(); j++)
            var.var_harmonic.output_mesh_param(j).
              SetFileName(2, GetBaseString(var.var_harmonic.output_mesh_param(j).GetTotalFieldFile())
                          + "." + numero.substr(numero.size()-3,3) + string(".bb"));
          
          for (int j = 0; j < var.var_harmonic.output_grid_param.GetM(); j++)
            var.var_harmonic.output_grid_param(j).
              SetFileName(2, GetBaseString(var.var_harmonic.output_grid_param(j).GetTotalFieldFile())
                          + numero + string(".dat"));
          
          var.var_harmonic.WriteOutputFile(eigen_mode_v, 2);
        }
      
      eigen_mode_v.Nullify();
      return;
    }
  
  // schema en temps
  int nb_iter = (var.GetFinalTime() - var.GetInitialTime())/var.GetTimeStep();  
  VectReal_wp Un(nb_dof_u), Unm1(nb_dof_u), Unp1(nb_dof_u), KhUn(nb_dof_u);
  VectReal_wp w(nb_dof_u), Q0(nb_dof_u), Q1(nb_dof_u), Q2(nb_dof_u), PUn(nb_dof_u), ImPUn(nb_dof_u);
  Un.Fill(0); Unm1.Fill(0); Unp1.Fill(0); KhUn.Fill(0);
  Q0.Fill(0); Q1.Fill(0); Q2.Fill(0); PUn.Fill(0); ImPUn.Fill(0);
  //ShUnm1.Fill(0);
  //Real_wp dt_loc = var.GetTimeStep()/(2.0*r+1);
  //DISP(dt_loc); Dh.Write("Dh.dat");
  for (int nt = 0; nt <= nb_iter+1; nt++)
    {
      Real_wp t = var.GetInitialTime() + nt*var.GetTimeStep();
      var.WriteSnapshot(nt, t, Un);
      
      // schema saute-mouton
      Mlt(-1.0, Kh, Un, KhUn);
      if (t < var.GetFinalTimeSource())
	var.AddPrimitiveTimeSource(1.0, t, 0, KhUn);
      
      MltAdd(-1.0, Sh, Unm1, 1.0, KhUn);
      for (int i = 0; i < nb_dof_u; i++)
	KhUn(i) /= Dh(i);
      
      for (int i = 0; i < nb_dof_u; i++)
	Unp1(i) = Unm1(i) + 2.0*var.GetTimeStep()*KhUn(i);
      
      Copy(Un, Unm1);
      Copy(Unp1, Un);
      
      // schema Runge-Kutta
      /*
      Mlt(-1.0, Ah, Un, KhUn);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(1.0, t, 0, KhUn);
      
      for (int i = 0; i < nb_dof_u; i++)
	KhUn(i) /= Dh(i);
      
      Copy(Un, Unp1);
      Add(var.GetTimeStep()/6.0, KhUn, Unp1);
      
      for (int i = 0; i < nb_dof_u; i++)
	Q0(i) = Un(i) + 0.5*var.GetTimeStep()*KhUn(i);
      
      Mlt(-1.0, Ah, Q0, KhUn);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(1.0, t+0.5*var.GetTimeStep(), 0, KhUn);
      
      for (int i = 0; i < nb_dof_u; i++)
	KhUn(i) /= Dh(i);

      Add(var.GetTimeStep()/3.0, KhUn, Unp1);

      for (int i = 0; i < nb_dof_u; i++)
	Q0(i) = Un(i) + 0.5*var.GetTimeStep()*KhUn(i);
      
      Mlt(-1.0, Ah, Q0, KhUn);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(1.0, t+0.5*var.GetTimeStep(), 0, KhUn);
      
      for (int i = 0; i < nb_dof_u; i++)
	KhUn(i) /= Dh(i);

      Add(var.GetTimeStep()/3.0, KhUn, Unp1);
      
      for (int i = 0; i < nb_dof_u; i++)
	Q0(i) = Un(i) + var.GetTimeStep()*KhUn(i);
      
      Mlt(-1.0, Ah, Q0, KhUn);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(1.0, t+var.GetTimeStep(), 0, KhUn);
      
      for (int i = 0; i < nb_dof_u; i++)
	KhUn(i) /= Dh(i);

      Add(var.GetTimeStep()/6.0, KhUn, Unp1);
      Copy(Unp1, Un);
      */
      // schema saute-mouton avec pas de temps local      
      //Un.Write("Un.dat");
      /*Mlt(P, Un, PUn);
      Copy(Un, ImPUn); Add(-1.0, PUn, ImPUn); //ImPUn.Write("ImPUn.dat");
      
      Mlt(-1.0, Kh, ImPUn, w);
      MltAdd(-1.0, Sh, Unm1, 1.0, w);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(1.0, t, 0, w);
      
      //w.Write("w.dat");
      //DISP(t); DISP(var.tlimit_source);
      
      Copy(Un, Q0);
      Mlt(-1.0, KhP, Q0, Q1); //Q1.Write("KhP_Q0.dat");
      for (int i = 0; i < Q1.GetM(); i++)
	Q1(i) = dt_loc*(w(i) + Q1(i))/Dh(i);
      
      //Q1.Write("Q1.dat");
      
      for (int k = 1; k <= 2*r; k++)
	{
	  Mlt(-1.0, KhP, Q1, KhUn);
	  if (k%2 == 0)
	    for (int i = 0; i < Q1.GetM(); i++)
	      Q2(i) = Q0(i) + 2.0*dt_loc*(KhUn(i) + w(i))/Dh(i);
	  else
	    for (int i = 0; i < Q1.GetM(); i++)
	      Q2(i) = Q0(i) + 2.0*dt_loc*KhUn(i)/Dh(i);
	  
	  //Q2.Write("Qn"+to_str(k)+".dat");
	  Copy(Q1, Q0);
	  Copy(Q2, Q1);
	}
      
      for (int i = 0; i < Un.GetM(); i++)
	Unp1(i) = Unm1(i) + 2.0*Q1(i);
      */
      //Unp1.Write("Unp1.dat");
      //int test_input; cout << "we wait " << endl; cin >> test_input;
      // schema saute-mouton standard
      
      /*Mlt(Kh, Un, KhUn);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(-1.0, t, 0, KhUn);
      
      Mlt(Sh, Unm1, Q0);
      for (int i = 0; i < N; i++)
	Unp1(i) = Unm1(i) + 2.0*var.GetTimeStep()*(-Q0(i) - KhUn(i))/Dh(i);
      */
      
      
      // schema splitte Runge-Kutta
      /*
      Mlt(Sh, Un, Q0);
      for (int i = 0; i < N; i++)
	Q1(i) = Un(i) - 0.5*var.GetTimeStep()*Q0(i)/Dh(i);
  
      Mlt(-1.0, Kh, Q1, KhUn);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(1.0, t, 0, KhUn);
      
      for (int i = 0; i < N; i++)
	KhUn(i) /= Dh(i);
      
      Copy(Q1, Unp1);
      Add(var.GetTimeStep()/6.0, KhUn, Unp1);
      
      for (int i = 0; i < N; i++)
	Q0(i) = Q1(i) + 0.5*var.GetTimeStep()*KhUn(i);
      
      Mlt(-1.0, Kh, Q0, KhUn);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(1.0, t+0.5*var.GetTimeStep(), 0, KhUn);
      
      for (int i = 0; i < N; i++)
	KhUn(i) /= Dh(i);

      Add(var.GetTimeStep()/3.0, KhUn, Unp1);

      for (int i = 0; i < N; i++)
	Q0(i) = Q1(i) + 0.5*var.GetTimeStep()*KhUn(i);
      
      Mlt(-1.0, Kh, Q0, KhUn);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(1.0, t+0.5*var.GetTimeStep(), 0, KhUn);
      
      for (int i = 0; i < N; i++)
	KhUn(i) /= Dh(i);

      Add(var.GetTimeStep()/3.0, KhUn, Unp1);
      
      for (int i = 0; i < N; i++)
	Q0(i) = Q1(i) + var.GetTimeStep()*KhUn(i);
      
      Mlt(-1.0, Kh, Q0, KhUn);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(1.0, t+var.GetTimeStep(), 0, KhUn);
      
      for (int i = 0; i < N; i++)
	KhUn(i) /= Dh(i);

      Add(var.GetTimeStep()/6.0, KhUn, Unp1);
      
      Mlt(Sh, Unp1, Q0);
      for (int i = 0; i < N; i++)
	Unp1(i) = Unp1(i) - 0.5*var.GetTimeStep()*Q0(i)/Dh(i);
      */       
      // schema splitte saute-mouton
      /*
      Mlt(Sh, Unm1, Q0);
      for (int i = 0; i < N; i++)
	Q1(i) = Unm1(i) - 2.0*var.GetTimeStep()*Q0(i)/Dh(i);
      
      //Copy(Un, Q1);
      
      Mlt(Kh, Un, KhUn);
      if (t < var.tlimit_source)
	var.AddPrimitiveTimeSource(-1.0, t, 0, KhUn);
      
      for (int i = 0; i < N; i++)
	Q2(i) = Q1(i) - 2.0*var.GetTimeStep()*KhUn(i)/Dh(i);
      
      //Copy(Q2, Unp1);
      Mlt(Sh, Q2, Q0);
      for (int i = 0; i < N; i++)
	Unp1(i) = Q2(i) - var.GetTimeStep()*Q0(i)/Dh(i);
      
      Copy(Un, Unm1);
      Copy(Unp1, Un);
      */           
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 2)
    {
      cout << "Entrez le fichier de donnees" << endl;
      abort();
    }

  string input_file(argv[1]);

  string type_element, type_equation;
  getElement_Equation(input_file, type_element, type_equation);

  HyperbolicProblem<TimeAeroAcousticEquation<Dimension2> > var;
  RunAll(var, input_file, type_element);

  return FinalizeMontjoie();
}
