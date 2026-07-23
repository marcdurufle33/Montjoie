#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_HCURL_AXI
#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

//#define MONTJOIE_WITH_INTEGRAL_EQUATION
#define MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
#define MONTJOIE_WITH_TRANSMISSION

#include "Elliptic/Maxwell/MontjoieMaxwellAxi.hxx"

using namespace Montjoie;

template<class TypeEquation>
void AddSurfaceOutput(const MeshInterpolation<Dimension2>& var_interp,
                      const VectComplex_wp& x_sol, 
                      const EllipticProblem<TypeEquation>& var,
                      VectR2& PosPoints, VectComplex_wp& ValHy)
{
  typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;

  Vector<VectComplex_wp> U0_vec(3);
  int N1 = var.offset_dof_unknown(1);
  int N2 = var.offset_dof_unknown(2) - N1;
  U0_vec(0).SetData(N1, const_cast<Complex_wp*>(&x_sol(0)));
  U0_vec(1).SetData(N2, const_cast<Complex_wp*>(&x_sol(N1)));
  if (var.FirstOrderFormulation())
    U0_vec(2).SetData(x_sol.GetM()-(N1+N2), const_cast<Complex_wp*>(&x_sol(N1+N2)));
  
  // on calcule H
  Vector<VectComplex_wp> EvalH_Nodal;
  var.EvaluateH_MixedFormulation(U0_vec, EvalH_Nodal);

  U0_vec(0).Nullify();
  U0_vec(1).Nullify();
  if (var.FirstOrderFormulation())
    U0_vec(2).Nullify();

  VectR2 s; R2 pt_loc;
  Matrix2_2 dfj, dfjm1;
  SetPoints<Dimension2> PointsElem;
  
  // boucle sur les aretes
  int r = var_interp.GetNbSubdivisions();
  const VectReal_wp& step = var_interp.GetSubdivisionStep();
  Vector<TinyVector<Complex_wp, 3> > valH;
  VectR2 pt_glob;
  for (int i = 0; i < var_interp.GetNbBoundary(); i++)
    {
      int num_elem = var_interp.GetElementNumberOfSurface(i);
      int num_loc = var_interp.GetLocalPositionOfSurface(i);
      
      const ElementReference<Dimension2, 2>& Fb_hcurl = var.GetReferenceElementHcurl(num_elem);
      const ElementReference<Dimension2, 1>& Fb_h1 = var.GetReferenceElementH1(num_elem, 1);
      var.mesh.GetVerticesElement(num_elem, s);
      Fb_h1.FjElemNodal(s, PointsElem, var.mesh, num_elem);

      valH.Reallocate(r+1);
      pt_glob.Reallocate(r+1);
      for (int j = 0; j <= r; j++)
	{
	  Real_wp t_loc = step(j);
	  Fb_h1.GetLocalCoordOnBoundary(num_loc, t_loc, pt_loc);
	  
	  Fb_h1.Fj(s, PointsElem, pt_loc, pt_glob(j), var.mesh, num_elem);
	  Fb_h1.DFj(s, PointsElem, pt_loc, dfj, var.mesh, num_elem);
	  GetInverse(dfj, dfjm1);
	  var.ComputeHpolar(EvalH_Nodal, num_elem, pt_loc,
			    pt_glob(j), dfjm1, valH(j), Fb_hcurl, Fb_h1);
	}
      
      for (int j = 0; j <= r; j++)
        {
          ValHy(i*(r+1) + j) += valH(j)(1);
          PosPoints(i*(r+1) + j) = pt_glob(j);
        }
    }
}


#ifdef SELDON_WITH_MPI
template<class T>
void GatherVector(Vector<T>& x, int root, MPI_Comm& comm)
{
  Vector<int64_t> xtmp;
  int rank_proc; MPI_Comm_rank(comm, &rank_proc);
  int nb_proc; MPI_Comm_size(comm, &nb_proc);
  if (rank_proc != root)
    {
      int n = x.GetM();
      MPI_Ssend(&n, 1, MPI_INTEGER, root, 8, comm);
      MpiSsend(comm, x, xtmp, n, root, 9);
    }
  else
    {
      IVect num(nb_proc);
      num.Zero();
      MPI_Status status;
      int nb = x.GetM();
      for (int p = 0; p < nb_proc; p++)
        if (p != root)
          {
            int n;
            MPI_Recv(&n, 1, MPI_INTEGER, p, 8, comm, &status);
            num(p) = n;
            nb += n;
          }
      
      int nb_old = x.GetM();
      x.Resize(nb);
      nb = nb_old;
      for (int p = 0; p < nb_proc; p++)
        if (p != root)
          {            
            MpiRecv(comm, &x(nb), xtmp, num(p), p, 9, status);
            nb += num(p);
          }            
    }
}
#endif


template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& name_file,
	    string& name_element, const string& name_equation)
{
  typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;

  MontjoieTimer var_chrono;
  var_chrono.SetMessage("GlobalComputation", "perform the overall simulation");
  var_chrono.Start("GlobalComputation");
  
  var.InitIndices(PhysicalConstant::nb_max_indices);
  var.SetTypeEquation("none");
  var.SetTypeElement(name_element);
  DISP(name_element);
  
  Vector<string> lines_data_file;
  
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(var.comm_group_mode, &rank_proc);
  int nb_proc; MPI_Comm_size(var.comm_group_mode, &nb_proc);
  ReadLinesFile(name_file, lines_data_file, var.comm_group_mode);
#else
  int rank_proc(0), nb_proc(1);
  ReadLinesFile(name_file, lines_data_file);
#endif

  ReadInputFile(lines_data_file, var);
  
#ifdef MONTJOIE_WITH_INTEGRAL_EQUATION
  VarIntegralEquationAxi<TypeEquation> var_integral(var);
  ReadInputFile(lines_data_file, var_integral);
#endif

#ifdef SELDON_WITH_MPI
  var_chrono.SetCommunicator(var.comm_group_mode);
#endif

  var_chrono.SetMessage("MeshGeneration", "construct and number the mesh");
  var_chrono.Start("MeshGeneration");
  
  // we read mesh and construct reference element (for finite element method)
  var.ComputeMeshAndFiniteElement(name_element);

  if (nb_proc == 1)
    var.mesh.Write("test.mesh");
  
  var_chrono.Stop("MeshGeneration");
 
  All_LinearSolver glob_solver(var);
  VarComputationRCS<TypeEquation> rcs_param(var);
  ReadInputFile(lines_data_file, rcs_param);

  Dimension2 dim; Symmetric property;
  MeshNumbering<Dimension2>& mesh_num = var.GetMeshNumbering(0);
  glob_solver.SelectOptimalLinearSolver(mesh_num.GetOrder(), var.GetNbDof(), dim, property);
  
  ReadInputFile(lines_data_file, glob_solver);

  if (rank_proc == 0)
    cout << "Frequency = " << var.GetFrequency() << endl;

  var_chrono.SetMessage("Initialisation", "construct geometric quantities and other arrays");
  var_chrono.Start("Initialisation");
  
  // right hand side, solution and vector used to display solution
  VectComplex_wp source_rhs, output_vector;
  VectR3_Complex_wp RCS_Vector;
  var.InitRcs(rcs_param);
  if (var.NumberOfModesToBeComputed())
    var.ComputeListMode(rcs_param);
  
  var.InitBesselArray(rcs_param);
  var.ComputeMassMatrix();
  
  var.PerformOtherInitializations();

  GlobalGenericMatrix<Complex_wp> nat_mat;
  
#ifdef MONTJOIE_WITH_INTEGRAL_EQUATION
  // initialization for integral equation
  if (var_integral.coupling_integral_equation != var_integral.FEM_ONLY)
    {
      var_integral.GetDofsBoundary(var.LINE_ABSORBING, var_integral.mesh_ie);
      var_integral.InitComputation_IntegralEquation();
      
      if (var_integral.type_output_file_ie/2 == 0)
	var_integral.ComputeSectionGrid_IE(var_integral.mesh_ie);
      
      if (var_integral.type_output_file_ie/2 == 1)
	{
	  var_integral.ConstructSurfacicMesh_IE(var.LINE_ABSORBING, var_integral.mesh_ie, var_integral.mesh_refined_ie);
	  var_integral.mesh_boundary_inside.Write(var_integral.name_file_mesh3D_ie);
	}
    }
#endif
  
  // output on surfaces
  MeshInterpolation<Dimension2> var_interp;
  Vector<Vector<Complex_wp> > ValHy;
  VectR2 PosPoints;
  if (var.ref_outputJ.GetM() > 0)
    {
      IVect ref_cond(var.mesh.GetNbReferences()+1);  
      ref_cond.Fill(0);  
      for (int i = 0; i < var.ref_outputJ.GetM(); i++)
        ref_cond(var.ref_outputJ(i)) = 1;
      
      Mesh<Dimension2> mesh_subdiv;
      int r = var.nb_subdiv_outputJ;
      
      var_interp.SetRegularSubdivisions(r);
      var_interp.ComputeSurfaceMesh(ref_cond, var.mesh, mesh_subdiv);
      //mesh_subdiv.Write("mesh_subdiv.mesh");
      
      var_interp.InitProjectionSurface(var.mesh);
      //DISP(var_interp.PointsReferenceSurface());
      ValHy.Reallocate(var.GetNbModes());
      for (int n = 0; n < var.GetNbModes(); n++)
        {
          ValHy(n).Reallocate(var_interp.GetNbBoundary()*(r+1));
          ValHy(n).Zero();
        }
      
      PosPoints.Reallocate(var_interp.GetNbBoundary()*(r+1));
    }
  
  int nb_source = var.GetNbRightHandSide(rcs_param);
  Vector<Complex_wp> x_sol(var.GetNbDof());
  Matrix<Complex_wp, General, ColMajor> M_sol(var.GetNbDof(), nb_source);
  M_sol.Fill(0); x_sol.Fill(0);
  
  var_chrono.Stop("Initialisation");
  var_chrono.Stop("GlobalComputation");

  var_chrono.SetMessage("Factorisation", "factorize finite element matrices");
  var_chrono.SetMessage("RightHandSide", "compute the right hand sides");
  var_chrono.SetMessage("Resolution", "solve the linear systems");
  var_chrono.SetMessage("RCS", "compute the radar cross section");
  var_chrono.SetMessage("Output", "post-process the solution");
  var_chrono.SetMessage("ComputationNorm", "compute the one-norm of matrices");
  
  // loop on the modes
  var.last_output_solution = false;
  VectReal_wp normE, normH; int prev_mode = 1000000; Real_wp norme_max(0);
  Vector<int> list_mode;
  for (int n = 0; n < var.GetNbModes(); n++)
    {
      var_chrono.Start("GlobalComputation");
      
      // mode number, that we consider
      list_mode.PushBack(var.GetModeNumber(n));
      var.SetCurrentModeNumber(var.GetModeNumber(n));
      var.UpdateDirichlet(n);
      if (rank_proc == 0)
        cout << "We solve mode number " << var.GetCurrentModeNumber() << endl;
      
      if (n == var.GetNbModes()-1)
	var.last_output_solution = true;

#ifdef MONTJOIE_WITH_INTEGRAL_EQUATION
      if (var_integral.coupling_integral_equation == var_integral.FEM_ONLY)
	{
#endif
	  var_chrono.Start("Factorisation");
	  
	  if (abs(prev_mode) != abs(var.GetCurrentModeNumber()))
	    glob_solver.PerformFactorizationStep(nat_mat);
	  
	  var_chrono.Stop("Factorisation");
	  
	  var_chrono.Start("RightHandSide");
	  
	  for (int k = 0; k < nb_source; k++)
	    {
	      if (rcs_param.GetRcsType() == rcs_param.MONOSTATIC_RCS)
		{
		  // modification of wave vector and polarization
		  R3 polar3D, kwave3D;
		  rcs_param.GetWaveVectorMonostatic(k, var.GetOmega(), kwave3D, polar3D);
		  var.SetPolarization(polar3D); var.SetWaveVector(kwave3D);
		}
	      
	      // computation of right hand side finite element
	      var.ComputeRightHandSide(x_sol);
	      if ((var.print_level >= 6) && (k == 0))
		{
                  if (nb_proc == 1)
                    x_sol.Write("rhs" + to_str(n) + ".dat");
                  else
                    x_sol.Write("rhs" + to_str(n) + "_P" + to_str(rank_proc) + ".dat");
                }
	      
	      SetCol(x_sol, k, M_sol);
	    }
	  
	  var_chrono.Stop("RightHandSide");
	  var_chrono.Start("Resolution");
	  
	  // we solve the linear system LU x = b
	  glob_solver.ComputeSolution(M_sol, nat_mat);
	  	  
	  if (rank_proc == 0)
	    cout<<"Solution computed successfully"<<endl;

	  if (var.print_level >= 6)
	    {
	      GetCol(M_sol, 0, x_sol);
              if (nb_proc == 1)
                x_sol.Write("sol" + to_str(n) + ".dat");
              else
                x_sol.Write("sol" + to_str(n) + "_P" + to_str(rank_proc) + ".dat");
	    }
	  
	  var_chrono.Stop("Resolution");
	  var_chrono.Start("RCS");
	  
          rcs_param.ComputeRCS(M_sol, RCS_Vector);
	  
	  var_chrono.Stop("RCS");

	  if (nb_source == 1)
	    {
	      var_chrono.Start("Output");
	      
	      GetCol(M_sol, 0, x_sol);
	      var.WriteDatas(x_sol);
	      
	      if (var.ref_outputJ.GetM() > 0)
		AddSurfaceOutput(var_interp, x_sol, var, PosPoints, ValHy(n));
	      
	      var_chrono.Stop("Output");
	    }
	  
	  var_chrono.Start("ComputationNorm");
	  
	  if (abs(prev_mode) != abs(var.GetCurrentModeNumber()))
	    var_chrono.DisplayAll();
	  
	  if (var.last_output_solution)
	    {
	      // exiting the loop over modes
	      break;
	    }
	  
	  // we retrieve the infinite norm of the solution
	  Real_wp x_inf = Norm1(M_sol);
#ifdef SELDON_WITH_MPI
	  Real_wp max_x = x_inf; Vector<int64_t> xtmp;
	  MpiAllreduce(var.comm_group_mode, &max_x, xtmp, &x_inf, 1, MPI_MAX);    
#endif
	  
	  if (rank_proc == 0)
	    cout << "1-Norm of x / norme_max = " << x_inf/norme_max << endl;
	  
	  norme_max = max(norme_max, x_inf);
	  
	  if ((var.NumberOfModesToBeComputed()) && (n >= 2))
	    {
	      if ((norme_max > 0) && (x_inf/norme_max < var.GetModeThreshold()))
		var.last_output_solution = true;
	    }
	  
	  var_chrono.Stop("ComputationNorm");
	  
	  
          //glob_solver.ClearFactorization();
#ifdef MONTJOIE_WITH_INTEGRAL_EQUATION
	}
      else if (var_integral.coupling_integral_equation == var_integral.BIE_ONLY)
	{
	  Matrix<Complex_wp> mat_integral; VectComplex_wp rhs_integral; IVect pivot;
	  var_integral.ComputeMatrixIntegralEquation(mat_integral);
	  //mat_integral.Write(string("Mh_m")+to_str(var.GetCurrentModeNumber())+".dat");
	  GetLU(mat_integral, pivot);
	  
	  var_integral.ComputeRightHandSideIntegralEquation(rhs_integral, 0);
	  
	  SolveLU(mat_integral, pivot, rhs_integral);
	  
	  // and get currents J
	  var_integral.GetComponentsSolution_IntegralEquation(rhs_integral, x_sol);
		
	  // write results on files
	  var_integral.WriteOutputFile_CurrentsJ(x_sol, var_integral.name_file_Jt,
						 var_integral.name_file_Jb, var_integral.name_file_Jnorme, true);
	}
#endif

      var_chrono.Stop("GlobalComputation");	
      prev_mode = var.GetCurrentModeNumber();
    } // end loop over modes
  
  if (var.ref_outputJ.GetM() > 0)
    {
      if (rank_proc == 0)
        list_mode.WriteText("mode_number.dat");
      
      if (var.output_Hy_teta)
        {
          // sorting points with theta parameter
          VectReal_wp ValTheta(PosPoints.GetM());
          Real_wp radius, theta;
          for (int i = 0; i < PosPoints.GetM(); i++)
            {
              CartesianToPolar(PosPoints(i)(0), PosPoints(i)(1), radius, theta);
              ValTheta(i) = theta;
            }
          
          /* RemoveDuplicate(ValTheta, ValHy);
          
#ifdef SELDON_WITH_MPI
          if (nb_proc > 1)
            {
              GatherVector(ValTheta, 0, var.comm_group_mode);
              GatherVector(ValHy, 0, var.comm_group_mode);
              if (rank_proc == 0)
                RemoveDuplicate(ValTheta, ValHy);
            }
#endif
          
          if (rank_proc == 0)
            {
              ofstream file_out(var.name_file_outputJ.data());
              file_out.precision(14);
              for (int i = 0; i < ValHy.GetM(); i++)
                file_out << ValTheta(i) << " " << real(ValHy(i)) << " " << imag(ValHy(i)) << endl;
              
              file_out.close();
            }
          */
        }
      else
        {
          // Gathering points and values of Hy on the nodal points of all the edges
          VectReal_wp PosX(PosPoints.GetM()), PosY(PosPoints.GetM());
          for (int i = 0; i < PosPoints.GetM(); i++)
            {
              PosX(i) = PosPoints(i)(0);
              PosY(i) = PosPoints(i)(1);
            }
          
#ifdef SELDON_WITH_MPI
          if (nb_proc > 1)
            {
              GatherVector(PosX, 0, var.comm_group_mode);
              GatherVector(PosY, 0, var.comm_group_mode);
              for (int n = 0; n < var.GetNbModes(); n++)
                GatherVector(ValHy(n), 0, var.comm_group_mode);
            }
#endif          
          
          if (rank_proc == 0)
            {
              // for root processor, PosPoints is recomposed
              PosPoints.Reallocate(PosX.GetM());
              for (int i = 0; i < PosPoints.GetM(); i++)
                PosPoints(i).Init(PosX(i), PosY(i));
              
              // PosNodes contains only the two vertices of each edge
              int r = var.nb_subdiv_outputJ;
              PosX.Clear(); PosY.Clear();
              VectR2 PosNodes(PosPoints.GetM()/(r+1)*2);
              for (int i = 0; i < PosPoints.GetM()/(r+1); i++)
                {
                  PosNodes(2*i).Init(PosPoints(i*(r+1))(1), PosPoints(i*(r+1))(0));
                  PosNodes(2*i+1).Init(PosPoints(i*(r+1)+r)(1), PosPoints(i*(r+1)+r)(0));
                }
              
              // these vertices are sorted with respect to y then x
              // permutation array and reciprocal array are stored
              IVect permut(PosNodes.GetM()); permut.Fill();
              Sort(PosNodes, permut);
              
              IVect inv_permut(permut.GetM());
              for (int i = 0; i < permut.GetM(); i++)
                inv_permut(permut(i)) = i;
              
              // constructing ListPoints that will contain an ordered set of points along the boundary
              Vector<bool> PtUsed(PosNodes.GetM()); PtUsed.Fill(false);
              int num = 0; bool test_loop = true;
              IVect ListPoints(PosPoints.GetM()); int nb = 0;
              PtUsed(num) = true;
              while (test_loop)
                {
                  int p = permut(num); // position in the original array
                  // new_p is the other extremity of the edge
                  int new_p = -1;
                  if (p%2 == 0)
                    {
                      new_p = p+1;
                      for (int i = 0; i <= r; i++)
                        ListPoints(nb++) = (p/2)*(r+1)+i;
                    }
                  else
                    {
                      new_p = p-1;
                      for (int i = 0; i <= r; i++)
                        ListPoints(nb++) = (p/2)*(r+1)+r-i;
                    }
                  
                  // then we try to find a point of another edge that is equal
                  // to the selected extremity of the edge
                  num = inv_permut(new_p);
                  if ((num > 0) && (PosNodes(num-1) == PosNodes(num)) && (!PtUsed(num-1)))
                    {
                      num--;
                      PtUsed(num) = true;
                    }
                  else if ((num < PosNodes.GetM()-1)
                           && (PosNodes(num+1) == PosNodes(num)) && (!PtUsed(num+1)))
                    {
                      num++;
                      PtUsed(num) = true;
                    }
                  else
                    test_loop = false;                  
                }
              
              // now writing all the points without duplicate
              ofstream file_out(var.name_file_outputJ.data());
              file_out.precision(14);

              int nb_points = nb; nb = 0; int jmin=0;
              for (int i = 0; i < nb_points/(r+1); i++)
                {
                  if (i == 0)
                    jmin = 0;
                  else
                    {
                      jmin = 1;
                      nb++;
                    }

                  for (int j = jmin; j <= r; j++)
                    {
                      int p = ListPoints(nb); nb++;
                      file_out << PosPoints(p)(0) << " " << PosPoints(p)(1) << " ";
                      for (int n = 0; n < var.GetNbModes(); n++)
                        file_out << real(ValHy(n)(p)) << " " << imag(ValHy(n)(p)) << ' ';
                      
                      file_out << '\n';
                    }
                }
              
              file_out.close();
            }
        }
    }
  
  cout << " we destroy the variables " << endl;
  
} // end method RunAll


int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      
      EllipticProblem<HarmonicMaxwellEquation_HcurlAxi> Vars;      
      if (Vars.GetRankProcMode() == 0)
        cout << "Maxwell Solver with Axi-symmetric Hcurl quadrangles" << endl; 
      
      RunAll(Vars, file_name_data, type_element, type_equation);
      if (Vars.GetRankProcMode() == 0)
        cout << " we destroy the variables " << endl;
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"maxwell_axi.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
