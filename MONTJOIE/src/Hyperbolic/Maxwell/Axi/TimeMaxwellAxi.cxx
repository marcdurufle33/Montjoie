#ifndef MONTJOIE_FILE_TIME_MAXWELL_AXI_CXX

namespace Montjoie
{
  
  template<class TypeEquation> template<class Matrix1>
  void TimeMaxwell_Axi<TypeEquation>
  ::SetMinus_DofTeta(const VectBool& IsDofTeta, Matrix1& A) const
  {
    for (int j1 = 0; j1 < IsDofTeta.GetM(); j1++)
      if (IsDofTeta(j1))
	{
	  for (int k1 = 0; k1 < A.GetRowSize(j1); k1++)
	    {
	      int k = A.Index(j1,k1);
	      if (IsDofTeta(k))
		A.value(j1,k1) *= -1.0;
	    }
	}
  }
  
  template<class TypeEquation> 
  void TimeMaxwell_Axi<TypeEquation>::RunAll(const string& input_file, const string& name_element,
					     const string& name_equation)
  {
    this->var_harmonic.InitIndices(10);

    // The input file is read, see file Reading_InputFiles.cxx
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(this->var_harmonic.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(this->var_harmonic.comm_group_mode, &nb_proc);
    ReadLinesFile(input_file, lines_data_file, this->var_harmonic.comm_group_mode);
#else
    int rank_proc(0), nb_proc(1);
    ReadLinesFile(input_file, lines_data_file);
#endif

    ReadInputFile(lines_data_file, *this);

    this->glob_solver = this->var_harmonic.GetNewLinearSolver();
    ReadInputFile(lines_data_file, *this->glob_solver);
 
    int print_level = this->var_harmonic.print_level;   
    if ((rank_proc == 0)||(print_level >= 10))
      if (print_level >= 0)
	{
	  cout << "Proc " << rank_proc;
	  cout << " out of " << nb_proc << " input file read" << endl;
	}
    
    // Mesh is read, finite element constructed
    this->var_harmonic.SetSameNumberPeriodicDofs();
    // implementation of this method is done in VarProblem.cxx
    this->var_harmonic.ComputeMeshAndFiniteElement(name_element);
    
    if (print_level >= 6)
      cout << rank_proc << "Mesh and finite element constructed "<<endl;
    
    if (this->FirstOrderScheme())
      this->var_harmonic.SetFirstOrderFormulation(true);
    else
      {
	if (abs(this->var_harmonic.delta_penalization) != Real_wp(0))
	  {
	    cout<<"For second order scheme, you can't have penalization terms "<<endl;
	    abort();
	  }
      }

    // Dirichlet condition (finding dofs associated to Dirichlet condition)
    // method implemented in TreatBoundaryConditions.cxx
    if (print_level >= 5)
      cout << rank_proc << " Dirichlet condition ... " << endl;
    
    this->var_harmonic.PerformOtherInitializations();
    this->var_harmonic.TreatDirichletCondition();    

    // computation of the right-hand side (source terms)
    // this method is implemented in this file VarInstationary.cxx
    this->ComputeRightHandSide();
    
    if ((rank_proc == 0)||(print_level >= 10))
      if (print_level >= 5)
	cout << rank_proc << " right hand side computed " << endl;
    
    Real_wp rhs_max = Norm2(this->sparse_vector_source);
    this->var_harmonic.SetCurrentModeNumber(0);
    bool test_calcul_mode = true;
    while (test_calcul_mode)
      {
	cout<<"Resolution of mode "<<this->var_harmonic.GetCurrentModeNumber()<<endl;
	
	// computation of DF_i^{-1} and other variables (J_i, normales)
	// method is implemented in Harmonic_MassMatrix.cxx
	// this->var_harmonic.ompute = true;
	this->var_harmonic.ComputeMassMatrix();
	if ((rank_proc == 0)||(print_level >= 10))
	  if (print_level >= 5)
	    cout << rank_proc << " mass matrix OK " << endl;
	
	// computation of stiffness matrix
	GlobalGenericMatrix<Real_wp> nat_mat;
	nat_mat.SetCoefMass(0.0); nat_mat.SetCoefStiffness(1.0);
	this->Glob_mat_Kh->Clear();
	this->var_harmonic.AddMatrixFEM(*this->Glob_mat_Kh, nat_mat);
	
	if (print_level >= 5)
	  cout << rank_proc << " Mass matrix for unsteady problems " << endl;
	
	this->ComputeMassMatrix();    
	
	if (print_level >= 5)
	  cout << rank_proc << " Mass matrix computed " << endl;
	
#ifdef SELDON_WITH_MPI
	MPI_Barrier(this->var_harmonic.comm_group_mode);
#endif
	
	if (print_level >= 5)
	  cout << rank_proc << " Initialization iterations ... " << endl;
	
	// last initialization before time iterations
	this->InitTimeIterations();
	
	// mesh is written on the disk
	if (nb_proc == 1)
	  this->var_harmonic.mesh.Write(string("test.mesh"));
	
#ifdef SELDON_WITH_MPI
        MPI_Barrier(this->var_harmonic.comm_group_mode);
#endif
	
	if (print_level >= 5)
	  cout << rank_proc << " Launching iterations ... " << endl;
	
	glob_chrono.Start(VirtualTimer::ALL);
	
	// Running time iteraions
	this->RunTimeIterations();
	
	// on passe au mode suivant 
	int number_mode = this->var_harmonic.GetCurrentModeNumber();
	this->var_harmonic.SetCurrentModeNumber(number_mode+1);
	this->ComputeRightHandSide();
	Real_wp rhs_current = Norm2(this->sparse_vector_source);
	DISP(rhs_current/rhs_max); DISP(rhs_current);
	if (this->var_harmonic.NumberOfModesToBeComputed())
	  {
	    rhs_max = max(rhs_max, rhs_current);
	    if (rhs_current/rhs_max <= this->var_harmonic.GetThresholdSource())
	      test_calcul_mode = false;
	  }
	else
	  {
	    //if (TypeEquationStationary::number_mode > this->var_harmonic.number_maximum_mode)
	    //test_calcul_mode = false;
	    abort();
	  }
	
      }
  }
}

#define MONTJOIE_FILE_TIME_MAXWELL_AXI_CXX
#endif
