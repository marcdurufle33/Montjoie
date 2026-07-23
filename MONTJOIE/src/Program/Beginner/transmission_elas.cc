#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TRANSMISSION

#define MONTJOIE_WITH_NODAL_H1

#include "Hyperbolic/Elastic/MontjoieTimeElastic.hxx"

using namespace Montjoie;

#ifdef MONTJOIE_WITH_TRANSMISSION

template<class TypeEquation>
void ComputeTransmissionMatrix(VirtualMatrix<Real_wp>& mat_sp,
                               EllipticProblem<TypeEquation>& var,
                               GlobalGenericMatrix<Real_wp>& nat_mat)
{
  var.var_transmission.AddTransmissionTerms(1.0, nat_mat, mat_sp, 0, 0);
}

template<class TypeEquation>
void RunAll(HyperbolicProblem<TypeEquation>& var, const string& input_file,
	    const string& name_element, const string& name_equation)
{
  typedef typename TypeEquation::TypeEquationStationary::Dimension Dimension;
  int dim_N = Dimension::dim_N;
  var.var_harmonic.InitIndices(PhysicalConstant::nb_max_indices);
  
  // The input file is read, see file Reading_InputFiles.cxx
  Vector<string> lines_data_file;
  
#ifdef SELDON_WITH_MPI
  int nb_proc; MPI_Comm_size(var.var_harmonic.comm_group_mode, &nb_proc);
  ReadLinesFile(input_file, lines_data_file, var.var_harmonic.comm_group_mode);
#else
  int nb_proc(1);
  ReadLinesFile(input_file, lines_data_file);
#endif

  ReadInputFile(lines_data_file, var);
  ReadInputFile(lines_data_file, *var.glob_solver);
  
  // mesh is read, finite element constructed
  var.var_harmonic.mesh_num.SetSameNumberPeriodicDofs();
  var.var_harmonic.ComputeMeshAndFiniteElement(name_element);
  
  // additional dofs
  var.var_harmonic.PerformOtherInitializations();
  
  // right hand side
  var.ComputeRightHandSide();
  
  // mass matrix for instationary case
  var.var_harmonic.ComputeMassMatrix();
  var.var_harmonic.ComputeQuasiPeriodicPhase();
  
  if (nb_proc == 1)
    var.var_harmonic.mesh.Write("test.mesh");
  
  // stiffness matrix
  DistributedMatrix<Real_wp, General, ArrayRowSparse> Ah, Mh, Kh, Sh, Ch;
  DistributedMatrix<Real_wp, General, ArrayRowSparse> MhTrue, MhAdd;
  int N = var.var_harmonic.mesh_num.GetNbDof();
  Real_wp zero(0), one(1);
  DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym> MhBlock;
  int type_schema;
  All_MatrixLU<Real_wp> mat_lu;
  IVect DdlVol, IndexDdl;
  VectReal_wp invDh;
  
  if (var.var_time_scheme.GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_ORDER2)
    {
      cout << "Choisissez le schema a utiliser " << endl;
      cout << "1- Schema totalement explicite " << endl;
      cout << "2- Schema quasi explicite (terme 1/eta Sh implicite)" << endl;
      cout << "3- Schema semi-explicite (terme 1/eta Sh + Ch + part en eta de Kh implicite)" << endl;
      cout << "4- Schema semi-implicite (tous les termes de surface)" << endl;
      cin >> type_schema;
      //type_schema = 1;
      cout << "Schema choisi :  " << type_schema << endl;
      
      Mh.Reallocate(dim_N*N, dim_N*N);
      Kh.Reallocate(dim_N*N, dim_N*N);
      Sh.Reallocate(dim_N*N, dim_N*N);
      Ch.Reallocate(dim_N*N, dim_N*N);
      
      cout <<" Calcul de Mh " << endl;
      
      GlobalGenericMatrix<Real_wp> nat_mat;
      nat_mat.SetCoefMass(1.0);
      nat_mat.SetCoefStiffness(0.0);
      nat_mat.SetCoefDamping(0.0);
      var.var_harmonic.AddMatrixWithBC(Mh, nat_mat);
      
      cout <<" Calcul de Kh " << endl;
      
      nat_mat.SetCoefMass(0.0);
      nat_mat.SetCoefStiffness(1.0);
      var.var_harmonic.var_transmission.coef_ch = 0.0;
      var.var_harmonic.var_transmission.coef_sh = 0.0;      
      var.var_harmonic.AddMatrixWithBC(Kh, nat_mat);

      cout <<" Calcul de Ch " << endl;
      
      var.var_harmonic.var_transmission.coef_ch = 1.0;
      var.var_harmonic.var_transmission.coef_sh = 0.0;
      ComputeTransmissionMatrix(Ch, var.var_harmonic, nat_mat);

      cout <<" Calcul de Sh " << endl;
 
      var.var_harmonic.var_transmission.coef_ch = 0.0;
      var.var_harmonic.var_transmission.coef_sh = 1.0;
      ComputeTransmissionMatrix(Sh, var.var_harmonic, nat_mat);
      
      if (false)
	{
	  DistributedMatrix<Real_wp, General, ArrayRowSparse> KhTest(dim_N*N, dim_N*N);
	  var.var_harmonic.var_transmission.coef_ch = 1.0;
	  var.var_harmonic.var_transmission.coef_sh = 1.0;      
	  var.var_harmonic.AddMatrixWithBC(KhTest, nat_mat);
	  KhTest.WriteText("KhAll.dat");
	  Kh.WriteText("Kh.dat");
          Mh.WriteText("Mh.dat");
	  Ch.WriteText("Ch.dat");
	  Sh.WriteText("Sh.dat");
	}
      
      cout <<" Calcul de MhAdd " << endl;
      
      // modification of mass matrix Mh for semi-implicit schemes
      MhTrue = Mh;
      MhAdd.Reallocate(dim_N*N, dim_N*N);
      if (type_schema == 2)
	{
	  // forming Mh + dt^2/4 Ch
	  MhAdd = Ch;
	}
      else if ((type_schema == 3) || (type_schema == 4))
	{	  	  
	  // forming Mh + dt^2/4 (Ch+Sh)
	  MhAdd = Ch; Add(one, Sh, MhAdd);
	  
	  // stiffness part on the surface
	  const IVect& DdlVolMinus = var.var_harmonic.var_transmission.GetDofNumberOmegaMinus();
	  const IVect& DdlVolPlus = var.var_harmonic.var_transmission.GetDofNumberOmegaPlus();
	  DdlVol.Reallocate(DdlVolMinus.GetM()*2);
	  int nb = 0;
	  for (int i = 0; i < DdlVolMinus.GetM(); i++)
	    DdlVol(nb++) = DdlVolMinus(i);

	  for (int i = 0; i < DdlVolPlus.GetM(); i++)
	    DdlVol(nb++) = DdlVolPlus(i);
	  
	  IndexDdl.Reallocate(N);
	  IndexDdl.Fill(-1);
	  for (int i = 0; i < DdlVol.GetM(); i++)
	    IndexDdl(DdlVol(i)) = i;
	  
	  if (type_schema == 4)
	    {
	      for (int row_p = 0; row_p < dim_N; row_p++)
		for (int i = 0; i < DdlVol.GetM(); i++)
		  {
		    int row = DdlVol(i);
		    int iglob = row_p*N + row;
		    for (int j = 0; j < Kh.GetRowSize(iglob); j++)
		      {
			int col = Kh.Index(iglob, j);
			Real_wp val = Kh.Value(iglob, j);
			if (val != zero)
			  {
			    if (IndexDdl(col%N) >= 0)
			      MhAdd.AddInteraction(row + row_p*N, col, val);
			  }
		      }
		  }	  
	    }
	}

      cout <<" Calcul de Kh (somme) " << endl;

      Add(one, Ch, Kh); Add(one, Sh, Kh);

      if (var.GetTimeStep() > 1.0)
	{
	  // on evalue la CFL du schema
	  SparseEigenProblem<double, Matrix<double, Symmetric, ArrayRowSymSparse> > pb_eig;
	  
	  Matrix<double, Symmetric, ArrayRowSymSparse> K, M;
	  Matrix<double, General, ArrayRowSparse> KhTrue;
	  
	  KhTrue = Kh;
	  Add(-one, MhAdd, KhTrue);
	  
	  Copy(MhTrue, M);
	  Copy(KhTrue, K);
	  
	  pb_eig.SetStoppingCriterion(1e-4);
	  pb_eig.SetNbAskedEigenvalues(4);
	  pb_eig.SetNbArnoldiVectors(20);
	  pb_eig.SetComputationalMode(pb_eig.REGULAR_MODE);
	  pb_eig.SetTypeSpectrum(pb_eig.LARGE_EIGENVALUES, 0.0);
	  
	  pb_eig.SetPrintLevel(6);
	  pb_eig.SetCholeskyFactoForMass(true);
	  pb_eig.InitMatrix(K, M);
	  
	  Vector<Real_wp> eigen_values, eigen_imag;
	  Matrix<Real_wp, General, ColMajor> eigen_vec;
	  GetEigenvaluesEigenvectors(pb_eig, eigen_values, eigen_imag, eigen_vec);
    
	  DISP(eigen_values); DISP(eigen_imag);
	  Real_wp cfl = 0;    
	  for (int i = 0; i < pb_eig.GetNbAskedEigenvalues(); i++)
	    cfl = max(cfl, sqrt(eigen_values(i)*eigen_values(i) + eigen_imag(i)*eigen_imag(i)));
	  
	  cfl = 2.0/sqrt(cfl);
	  cout << "CFL of the scheme = " << cfl << endl;
	  exit(0);
	}
      
      Add(0.25*square(var.GetTimeStep()), MhAdd, Mh);
      
      Ch.Clear(); Sh.Clear();
      MhTrue.Clear(); MhAdd.Clear(); 
      
      if (type_schema >= 3)
	{
	  // non optimal version : we invert the global matrix Mh
	  // mat_lu.Factorize(Mh);
	  
	  // optimal version : inversion of Mh restricted to the surface
	  // extracting part on the surface
	  DistributedMatrix<Real_wp, General, ArrayRowSparse> MhSurf;
	  MhSurf.Reallocate(dim_N*DdlVol.GetM(), dim_N*DdlVol.GetM());
	  int Ns = DdlVol.GetM();
	  for (int row_p = 0; row_p < dim_N; row_p++)
	    for (int i = 0; i < DdlVol.GetM(); i++)
	      {
		int row = row_p*N + DdlVol(i);
		for (int j = 0; j < Mh.GetRowSize(row); j++)
		  {
		    int col = Mh.Index(row, j)%N;
		    int col_p = Mh.Index(row, j)/N;
		    Real_wp val = Mh.Value(row, j);
		    if ((IndexDdl(col) >= 0) && (val != zero))
		      MhSurf.AddInteraction(row_p*Ns + i, col_p*Ns + IndexDdl(col), val);
		  }
	      }
	  
	  // factorising this part
	  mat_lu.Factorize(MhSurf);
	  
	  // inverting diagonal part of Mh
	  invDh.Reallocate(dim_N*N);
	  invDh.Fill(one);
	  for (int i = 0; i < N; i++)
	    if (IndexDdl(i) == -1)
	      {
		invDh(i) = 1.0/Mh(i, i);
		invDh(i+N) = 1.0/Mh(i+N, i+N);
                if (dim_N == 3)
                  invDh(i+2*N) = 1.0/Mh(i+2*N, i+2*N);
	      }
	}
      else
	{
	  ConvertToBlockDiagonal(Mh, MhBlock);
	  GetInverse(MhBlock);
	}
    }
  else
    {
      var.ComputeStiffnessMatrix();
      var.ComputeMassMatrix();
    }

  cout <<" Initialisation des iterations en temps " << endl;
 
  var.InitTimeIterations();

  if (var.var_time_scheme.GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_ORDER2)
    {       
      // compressing stiffness matrix
      cout <<" Compression de Kh " << endl;
      DistributedMatrix<Real_wp, General, RowSparse> KhSparse;
      Copy(Kh, KhSparse);
      Kh.Clear();

      cout <<" On commence les iterations en temps " << endl;
      
      VectReal_wp Un(dim_N*N), Un_prev(dim_N*N), Un_next(dim_N*N), ProdUn(dim_N*N);
      Un.Fill(zero); Un_prev.Fill(zero); Un_next.Fill(zero); ProdUn.Fill(zero);

      int Ns = DdlVol.GetM();
      VectReal_wp ProdSurf(dim_N*Ns);
      
      int nb_max_iter = toInteger(round((var.GetFinalTime()-var.GetInitialTime())/var.GetTimeStep()));
      for (int n_time = 0; n_time < nb_max_iter; n_time++)
	{
	  Real_wp t = var.GetInitialTime() + n_time*var.GetTimeStep();
	  var.GiveIterate(n_time, t, Un);
	  
	  // advancing time scheme
	  // adding contribution of Kh
	  Mlt(KhSparse, Un, ProdUn);
	  
	  // source term
	  var.AddScalarSourceAtTime(-one, t, 1, ProdUn);
          
          #pragma omp parallel for
	  for (int i = 0; i < Un.GetM(); i++)
	    Un_next(i) = 2.0*Un(i) - Un_prev(i);
	  
	  if (type_schema >= 3)
	    {
	      // non-optimal version : resolution with the global matrix Mh
	      // mat_lu.Solve(ProdUn);
	      
	      // optimal version, extracting the components on the surface
	      for (int i = 0; i < DdlVol.GetM(); i++)
		{
		  ProdSurf(i) = ProdUn(DdlVol(i));
		  ProdSurf(i+Ns) = ProdUn(DdlVol(i) + N);
                  if (dim_N == 3)
                    ProdSurf(i+2*Ns) = ProdUn(DdlVol(i) + 2*N);
		}
	      
	      // inverting them
	      mat_lu.Solve(ProdSurf);
	      
	      for (int i = 0; i < DdlVol.GetM(); i++)
		{
		  ProdUn(DdlVol(i)) = ProdSurf(i);
		  ProdUn(DdlVol(i) + N) = ProdSurf(i+Ns);
                  if (dim_N == 3)
                    ProdUn(DdlVol(i) + 2*N) = ProdSurf(i+2*Ns);
		}
	      
	      // and treating diagonal components (outside the surface)
	      for (int i = 0; i < ProdUn.GetM(); i++)
		ProdUn(i) *= invDh(i);

	      Add(-var.GetTimeStep()*var.GetTimeStep(), ProdUn, Un_next);
	    }
	  else
	    MltAdd(-var.GetTimeStep()*var.GetTimeStep(), MhBlock, ProdUn, one, Un_next);
	  
	  // updating iterates
	  Copy(Un, Un_prev); Copy(Un_next, Un);
	}
    }
  else
    var.RunTimeIterations();
}

#else

template<class TypeEquation>
void RunAll(HyperbolicProblem<TypeEquation>& var, const string& input_file, const string& name_element,
	    const string& name_equation)
{
  var.RunAll(input_file, name_element, name_equation);
}

#endif

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

      //HyperbolicProblem<TimeElasticEquation<Dimension2> > Vars;
      HyperbolicProblem<TimeElasticEquation<Dimension3> > Vars;

      if (Vars.var_harmonic.GetRankProcMode() == 0)
	cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;
      
      if (Vars.var_harmonic.GetRankProcMode() == 0)
	cout<<"Elastic Solver with DG Lobatto quadrilaterals"<<endl; 
      
      RunAll(Vars, file_name_data, type_element, type_equation);
      if (Vars.var_harmonic.GetRankProcMode() == 0)
	cout<<" we destroy the variables "<<endl;
      
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"main.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
}
