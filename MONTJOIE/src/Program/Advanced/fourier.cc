#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_LAPLACE
#define MONTJOIE_WITH_HELMHOLTZ

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"

using namespace Montjoie;


template<class TypeEquation>
void InitializeComputation(HyperbolicProblem<TypeEquation>& var_time, const string& input_file,
			   const string& name_element)
{
  // on recupere l'objet calculant le laplacien
  typedef typename TypeEquation::TypeEquationStationary TypeEquationStationary;
  EllipticProblem<TypeEquationStationary>& var = var_time.var_harmonic;
  
  // initialisation tableaux pour les indices physiques (rho, mu, sigma)
  var.InitIndices(PhysicalConstant::nb_max_indices);
  
  Vector<string> lines_data_file;
  
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  ReadLinesFile(input_file, lines_data_file, MPI_COMM_WORLD);
#else
  int rank_proc(0);
  ReadLinesFile(input_file, lines_data_file);
#endif

  // on lit le fichier de donnees
  ReadInputFile(lines_data_file, var_time);

  var_time.glob_solver = var.GetNewLinearSolver();
  ReadInputFile(lines_data_file, *var_time.glob_solver);
  
  // on construit le maillage et element fini
  var.ComputeMeshAndFiniteElement(name_element);
  
  if ((rank_proc == 0)||(var.print_level >= 10))
    if (var.print_level >= 2)
      cout << rank_proc << "Mesh and finite element constructed "<<endl;
  
  if (var_time.FirstOrderScheme())
    var.SetFirstOrderFormulation(true);
  
  // on trouve ddls avec Dirichlet
  var.PerformOtherInitializations();
  var.TreatDirichletCondition();    
  
  // on calcule la source spatiale
  var_time.ComputeRightHandSide();
  
  if ((rank_proc == 0)||(var.print_level >= 10))
    if (var.print_level >= 2)
      cout << rank_proc << " right hand side computed "<<endl;
    
  // calcul des donnes geometriques (jacobien Ji, et matrice jacobienne DFi)
  var.ComputeMassMatrix();
  var.ComputeQuasiPeriodicPhase();
  if ((rank_proc == 0)||(var.print_level >= 10))
    if (var.print_level >= 2)
      cout<<rank_proc<<" mass matrix OK "<<endl;
  
  // calcul de la matrice de rigidite
  var_time.ComputeStiffnessMatrix();
  
  if ((rank_proc == 0)||(var.print_level >= 10))
    if (var.print_level >= 3)
      cout<<rank_proc<<" Mass matrix for unsteady problems "<<endl;
  
  var_time.ComputeMassMatrix();    

  if ((rank_proc == 0)||(var.print_level >= 10))
    if (var.print_level >= 3)
      cout<<rank_proc<<" Mass matrix computed "<<endl;
    
#ifdef SELDON_WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
    
  if ((rank_proc == 0)||(var.print_level >= 10))
    if (var.print_level >= 5)
      cout<<rank_proc<<" Initialization iterations ... "<<endl;
  
  // last initialization before time iterations
  var_time.InitTimeIterations();
  
  // mesh is written on the disk
  if (var.GetNbProcPerMode() == 1)
    var.mesh.Write(string("test.mesh"));
  
#ifdef SELDON_WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  glob_chrono.Start(VirtualTimer::ALL);
  
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc == 2)
    {
      string file_name(argv[1]);
      
      string type_element, type_equation;
      getElement_Equation(file_name, type_element, type_equation);
      
      // on resout l'acoustique en 2-D
      HyperbolicProblem<AcousticEquation<Dimension2> > var;
      
      //VarComputationRCS<LaplaceEquation<Dimension2> > rcs_param(var.var_harmonic);
      //ReadInputFile(file_name, rcs_param);
      
      // on initialise les calculs 
      InitializeComputation(var, file_name, type_element);
      
      // on utilise un schema d'ordre 4
      ModifiedEquationSystemIterator<Real_wp> scheme;
      scheme.SetOrder(4);
      
      // on compte le nombre d'iterations
      int nb_max_iter = toInteger(ceil(abs(var.GetFinalTime() - var.GetInitialTime())/var.GetTimeStep()) );
      
      // conditions initiales nulles
      //int N = var.var_harmonic.GetNbDof();
      int Nscal = var.GetNbScalarUnknowns();
      int Nvec = var.GetNbVectorialUnknowns();

      VectReal_wp Uh(Nscal), Vh(Nvec);
      Real_wp t = var.GetInitialTime(), dt = var.GetTimeStep();
      Uh.Fill(0); Vh.Fill(0);
      var.GiveNumberIterations(dt, nb_max_iter);
      scheme.SetInitialConditionS(t, dt, Uh, Vh, var);
      
      // on calcule la transfo de fourier et on la stocke dans Uchap_r, Uchap_i
      VectReal_wp Uchap_r(Nscal), Uchap_i(Nscal); Uchap_r.Fill(0); Uchap_i.Fill(0);
      
      // boucle principale en temps
      for (int n_time = 0; n_time < nb_max_iter; n_time++)
	{
	  t = var.GetInitialTime() + dt*(n_time+1);
	  // on met a jour la transfo de Fourier
	  Complex_wp phase = exp(-Iwp*var.var_harmonic.GetOmega()*t );
	  Add(real(phase)*dt, scheme.GetIterate(), Uchap_r);
	  Add(imag(phase)*dt, scheme.GetIterate(), Uchap_i);
	  
	  // on avance le schema
	  scheme.Advance(t, n_time, var);
	}
      
      // et on imprime le resultat
      var.var_harmonic.WriteOutputFile(Uchap_r, "Uchap_real.dat");
      var.var_harmonic.WriteOutputFile(Uchap_i, "Uchap_imag.dat");
      
      VectComplex_wp Uchap(Nscal);
      for (int i = 0; i < Nscal; i++)
        Uchap(i) = Complex_wp(Uchap_r(i), Uchap_i(i));
      
      //rcs_param.ComputeRCS(var.var_harmonic, Uchap);
      
    }
  else
    {
      cout << "Entrez un fichier de donnees " << endl;
      abort();
    }
  
  FinalizeMontjoie();
}
