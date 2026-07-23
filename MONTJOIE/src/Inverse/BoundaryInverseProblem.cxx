#ifndef MONTJOIE_FILE_BOUNDARY_INVERSE_PROBLEM_CXX

#include "BoundaryInverseProblem.hxx"

namespace Montjoie
{
  
  //! list of exact functions that must be recovered
  template<class Dimension>
  typename BoundaryInverseProblem_Base<Dimension>::type_ptr_fctC
  BoundaryInverseProblem_Base<Dimension>::exact_fct_to_find[42];
  
  
  /*******************************
   * BoundaryInverseProblem_Base *
   *******************************/
  
  
  //! Default constructor
  template<class Dimension>
  BoundaryInverseProblem_Base<Dimension>
  ::BoundaryInverseProblem_Base(VarHarmonic_Base<Complex_wp, Dimension>& var,
				VarMigration_Base<Dimension>& var_mig)
    : var_helm(var), var_migration(var_mig)
  {
    measured_part = ALL_PART;
    nb_observables = 0;
    NdisplayC = 100;
    
    nb_comp_source_invert = 0;
    nb_angles_measure = 0;
    nb_angles_to_invert = 0;
    
    number_facto_lu = 0;
    number_solve = 0;
    write_experimental_data = false;

    file_name_experimental_data = "y_ref";    
  }
  
  
  //! Destructor
  template<class Dimension>
  BoundaryInverseProblem_Base<Dimension>::~BoundaryInverseProblem_Base()
  {
    for (int k = 0; k < linear_solver.GetM(); k++)
      delete linear_solver(k);
    
    for (int k = 0; k < c_discrete.GetM(); k++)
      if (c_discrete(k) != NULL)
	delete c_discrete(k);
  }
  
  
  //! changes the number of parameters to find
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::ResizeNbParameters(int n)
  {
    c_discrete.Reallocate(n);
    param_c_name.Reallocate(n);
    param_c_type.Reallocate(n);
    c_exact.Reallocate(n);
    for (int k = 0; k < n; k++)
      {
	c_discrete(k) = NULL;
	c_exact(k) = NULL;
      }
  }


  //! sets the exact parameter k
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::SetExactParameter(int k, const string& name, type_ptr_fctC fct)
  {
    c_exact(k) = fct;
    param_c_name(k) = name;
    param_c_type(k) = 0;
  }
  
  
  //! returns parameters associated with the iterative regularization method
  template<class Dimension>
  const Vector<string>& BoundaryInverseProblem_Base<Dimension>::GetRegularizationParameter() const
  {
    return regularization_parameter;
  }
  
  
  //! reads a line of the Montjoie data file
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::SetInputData(const string& keyword, const Vector<string>& param)
  {
    MultiFrequencyProblem::SetInputData(keyword, param);
    var_migration.SetInputData(keyword, param);
    
    if (keyword == "ApproximationC")
      approx_c_parameter = param;
    else if (keyword == "RegularizationMethod")
      regularization_parameter = param;
    else if (keyword == "DirichletC")
      dirichlet_c.PushBack(param);
    else if (keyword == "ReferenceC")
      {
	ref_volume.Reallocate(param.GetM());
	for (int k = 0; k < ref_volume.GetM(); k++)
	  ref_volume(k) = to_num<int>(param(k));
      }
    else if (keyword == "DisplayC")
      NdisplayC = to_num<int>(param(0));
    else if (keyword == "WriteExperimentalData")
      {
	if (param(0) == "YES")
	  write_experimental_data = true;
	else
	  write_experimental_data = false;

        if (param.GetM() > 1)
          file_name_experimental_data = param(1);
      }
    else if (keyword == "InputFrequencySource")
      {
        coef_input_frequency.ReadText(param(0));
        if (coef_input_frequency.GetM() != this->omega.GetM())
          {
            cout << "Number of frequencies = " << this->omega.GetM() << endl;
            cout << "Number of coefficients in file " << param(0) << " = " << coef_input_frequency.GetM() << endl;
            cout << "Size differ => correct the problem" << endl;
            abort();
          }
      }
    else if (keyword == "InitialGuess")
      {
	if (isalpha(param(0)[0]))
          {
            C_init_guess.Clear(); VectReal_wp Cval;
            for (int k = 0; k < param.GetM(); k++)
              {
                Cval.Read(param(k));
                C_init_guess.PushBack(Cval);
              }
          }
        else
          {
            C_init_guess.Reallocate(param.GetM());
            for (int k = 0; k < param.GetM(); k++)
              C_init_guess(k) = to_num<Real_wp>(param(k));
          }
      }
    else if (keyword == "ParamScaling")
      {
	param_scaling.Reallocate(param.GetM());
	for (int k = 0; k < param.GetM(); k++)
	  param_scaling(k) = to_num<Real_wp>(param(k));	
      }
    else if (keyword == "TypeInversion")
      this->SetInversionParameters(param);
    else if (keyword == "TypeMeasurement")
      {
        if (param(0) == "Real")
          measured_part = REAL_PART;
        else if (param(0) == "Imag")
          measured_part = IMAG_PART;
        else if (param(0) == "Modulus")
          measured_part = MODULUS_PART;
        else if (param(0) == "Complex")
          measured_part = ALL_PART;        
      }
    else if (keyword == "MeasureAngle")
      {
        this->nb_angles_measure = 1;
        coef_combination_measure.Reallocate(2);
        Real_wp teta = to_num<Real_wp>(param(0))*pi_wp/180;
        coef_combination_measure(0) = cos(teta);
        coef_combination_measure(1) = sin(teta);
      }
  }
  
  
  //! returns true if the exact solution has been stored
  template<class Dimension>
  bool BoundaryInverseProblem_Base<Dimension>::ExactSolutionKnown() const
  {
    return true;
  }
  
  
  //! returns the distance to the exact solution
  template<class Dimension>
  Real_wp BoundaryInverseProblem_Base<Dimension>::DistanceToExactSolution(const VectReal_wp& x)
  {
    return Norm2(x - C_store_exact) / Norm2(x);
  }
  

  //! From the vector cn0, we extract the different parameters c
  //! A shallow copy is made (you need to nullify cn after)
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::ExtractParametersC(const VectReal_wp& cn0, Vector<VectReal_wp>& cn)
  {
    cn.Reallocate(c_discrete.GetM());
    int offset = 0;
    for (int k = 0; k < c_discrete.GetM(); k++)
      {
	cn(k).SetData(c_discrete(k)->GetNbDof(),
		      const_cast<Real_wp*>(&cn0(offset)));
	
	offset += c_discrete(k)->GetNbDof();
      }
  }


  //! From the vector cn0, we extract the different parameters c
  //! the values are duplicated
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::CopyParametersC(const VectReal_wp& cn0, Vector<VectReal_wp>& cn)
  {
    cn.Reallocate(c_discrete.GetM());
    int offset = 0;
    for (int k = 0; k < c_discrete.GetM(); k++)
      {
	cn(k).Reallocate(c_discrete(k)->GetNbDof());
	for (int i = 0; i < cn(k).GetM(); i++)
	  cn(k)(i) = cn0(offset + i);
	
	offset += c_discrete(k)->GetNbDof();
      }
  }
  
  
  //! Computes the interpolation error between the exact c and interpolated c
  template<class Dimension>
  Real_wp BoundaryInverseProblem_Base<Dimension>::GetInterpolationError()
  {
    Vector<VectReal_wp> cn;
    CopyParametersC(C_store_exact, cn);
    
    Real_wp err(0), sum(0);
    for (int p = 0; p < cn.GetM(); p++)
      {
	for (int i = 0; i < weights_basis(p).GetM(); i++)
	  cn(p)(i) *= weights_basis(p)(i);
	
	c_discrete(p)->Init(cn(p));
	
	// c is written on the output grids defined in var_laplace
	for (int n = 0; n < var_helm.var_grid.GetM(); n++)
	  {
	    GridInterpolationFull<Dimension>& var_gr = var_helm.var_grid(n);
	    GridInterpolation<Dimension>& var_interp = var_helm.all_points_display;
	    const IVect& list_points = var_gr.GetPointNumber();
	    //const Vector<IVect>& list_points_proc = var_gr.GetPointNumberAllProc();
	    
	    int nnz = list_points.GetM();
	    int Nglob = var_gr.GetNbGlobalPoints();
	    VectReal_wp interp_index(Nglob), exact_index(Nglob);
	    interp_index.Zero(); exact_index.Zero();
	    for (int i1 = 0; i1 < nnz; i1++)
	      {
		int i = list_points(i1);
		if (i < 0)
		  continue;
		
		R_N pt_glob = var_interp.GetGlobalCoordinate(i);
		int num_elem = var_interp.GetElementNumber(i);
		R_N pt_loc = var_interp.GetLocalCoordinate(i);
		Real_wp cn_ex = (*c_exact(p))(pt_glob);
		Real_wp cn_app = c_discrete(p)->Evaluate(pt_glob, num_elem, pt_loc);
		
		err += square(cn_app - cn_ex);
		sum += square(cn_ex);
	      }
	  }
      }

    return sqrt(err/sum);
  }

  
  //! writes the parameters c in output files
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::WriteDatas(const VectReal_wp& cn0)
  {
    Vector<VectReal_wp> cn;
    CopyParametersC(cn0, cn);
    
    for (int p = 0; p < cn.GetM(); p++)
      {
	for (int i = 0; i < weights_basis(p).GetM(); i++)
	  cn(p)(i) *= weights_basis(p)(i);
	
	c_discrete(p)->Init(cn(p));
	
	// c is written on the output grids defined in var_laplace
	for (int n = 0; n < var_helm.var_grid.GetM(); n++)
	  {
	    GridInterpolationFull<Dimension>& var_gr = var_helm.var_grid(n);
	    GridInterpolation<Dimension>& var_interp = var_helm.all_points_display;
	    const IVect& list_points = var_gr.GetPointNumber();
	    const Vector<IVect>& list_points_proc = var_gr.GetPointNumberAllProc();
	    
	    int nnz = list_points.GetM();
	    int Nglob = var_gr.GetNbGlobalPoints();
	    VectReal_wp interp_index(Nglob), exact_index(Nglob);
	    interp_index.Zero(); exact_index.Zero();
	    for (int i1 = 0; i1 < nnz; i1++)
	      {
		int i = list_points(i1);
		if (i < 0)
		  continue;
		
		R_N pt_glob = var_interp.GetGlobalCoordinate(i);
		int num_elem = var_interp.GetElementNumber(i);
		R_N pt_loc = var_interp.GetLocalCoordinate(i);
		Real_wp cn_ex = (*c_exact(p))(pt_glob);
		Real_wp cn_app = c_discrete(p)->Evaluate(pt_glob, num_elem, pt_loc);
		
		int iglob = list_points_proc(0)(i1);
		interp_index(iglob) = cn_app;
		exact_index(iglob) = cn_ex;
	      }
	    
	    cout << " We write the data " << endl;
	    
	    WriteMatlab(exact_index, var_gr, string("Cexact_P") + to_str(p) + "_G" + to_str(n) + ".dat", 
			var_helm.OutputWrittenInDoublePrecision());
	    
	    WriteMatlab(interp_index, var_gr, string("Capp_P") + to_str(p) + "_G" + to_str(n) + ".dat", 
			var_helm.OutputWrittenInDoublePrecision());
	  }
      }

    int rank_proc = 0;
#ifdef SELDON_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#endif
    
    if (rank_proc == 0)
      {
        if (cn0.GetM() <= 50)
          cout << "Solution obtained : " << cn0 << endl;
        
        
        glob_chrono.Stop(VirtualTimer::ALL);
        cout << "Temps de resolution = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << "s" << endl;
        cout << "Number of factorizations = " << number_facto_lu << " \t Number of solves with single rhs = " << number_solve << endl;
        
        glob_chrono.Start(VirtualTimer::ALL);
      }
  }
  
    
  //! initialization of montjoie object var_helm
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::Construct(const Vector<string>& lines_data_file, const string& name_elt,
	      const string name_eq)
  {
    name_element = name_elt; name_equation = name_eq;
    
    // the input file is read
    ReadInputFile(lines_data_file, *this);

    if (this->omega.GetM() <= 0)
      {
	cout << "RangeFrequency absent" << endl;
	abort();
      }
    
    if (c_discrete.GetM() <= 0)
      {
        cout << "TypeInversion absent" << endl;
        abort();
      }
        
    this->DistributeFrequency(coef_input_frequency);

#ifdef SELDON_WITH_MPI
    var_helm.comm_group_mode = this->comm_group_mode;
#endif
    
    linear_solver.Reallocate(this->omega.GetM());
    
    var_helm.KeepGridLocalizationArrays();
    var_helm.ConstructAll(var_migration.input_file_simu, name_element, name_equation, linear_solver(0), true, false);
    
    linear_solver(0)->SetPrintLevel(-1);
    for (int k = 1; k < linear_solver.GetM(); k++)
      {
	linear_solver(k) = var_helm.GetNewLinearSolver();
	ReadInputFile(var_migration.input_file_simu, *linear_solver(k));
	
	linear_solver(k)->SetPrintLevel(-1);
      }
    
    this->InitIndices();
    
    int nb_proc(1);
#ifdef SELDON_WITH_MPI
    MPI_Comm_size(this->comm_freq, &nb_proc);
#endif
    
    if (nb_proc == 1)
      var_helm.mesh.Write("test.mesh");
    
    // computation of Ji, DFi on the mesh
    points_quad = var_helm.Glob_PointsQuadrature;
    var_helm.Glob_PointsQuadrature.Clear();
    
    // the source is computed
    var_migration.ComputeRightHandSide(var_helm, source_rhs);
    
    // the discretization for c is constructed
    for (int k = 0; k < c_discrete.GetM(); k++)
      c_discrete(k) = VariableParameterND_Base<Real_wp, Dimension>::
	GetNewApproximation(var_helm, approx_c_parameter);
    
    // we keep only elements with the good reference
    int ref_max = ref_volume.GetNormInf();
    for (int i = 0; i < var_helm.mesh.GetNbElt(); i++)
      ref_max = max(ref_max, var_helm.mesh.Element(i).GetReference());
    
    Vector<bool> is_ref_kept(ref_max+1);
    is_ref_kept.Fill(false);
    for (int i = 0; i < ref_volume.GetM(); i++)
      is_ref_kept(ref_volume(i)) = true;
    
    elt_to_keep.Reallocate(var_helm.mesh.GetNbElt());
    elt_to_keep.Fill(false);
    for (int i = 0; i < var_helm.mesh.GetNbElt(); i++)
      if (is_ref_kept(var_helm.mesh.Element(i).GetReference()))
	elt_to_keep(i) = true;
    
    //DISP(elt_to_keep);
    for (int k = 0; k < c_discrete.GetM(); k++)
      c_discrete(k)->KeepElements(elt_to_keep);
    
    // we find dofs of c for which values of c are enfored
    ConstructDirichlet();
    //DISP(DofDirichletC); DISP(ValueDirichletC);
    
    // normalization for basis functions of c
    weights_basis.Reallocate(c_discrete.GetM());
    for (int k = 0; k < c_discrete.GetM(); k++)
      c_discrete(k)->GetWeightFunction(weights_basis(k));
    
    //weights_basis.Reallocate(c_discrete->GetNbDof());
    //weights_basis.Fill(Real_wp(1));

    if (param_scaling.GetM() <= 0)
      {
	param_scaling.Reallocate(c_discrete.GetM());
	param_scaling.Fill(1);
      }
    
    for (int k = 0; k < c_discrete.GetM(); k++)
      for (int i = 0; i < weights_basis(k).GetM(); i++)
	weights_basis(k)(i) = param_scaling(k)/(sqrt(weights_basis(k)(i)));
    
    // the exact value of c (projection in the discrete space) is computed
    VectReal_wp evalC_exact;
    C_store_exact.Clear();
    for (int k = 0; k < c_discrete.GetM(); k++)
      {
	c_discrete(k)->ProjectFunction(*c_exact(k), evalC_exact);      
	C_store_exact.PushBack(evalC_exact);
      }
    
    int offset = 0;
    for (int k = 0; k < c_discrete.GetM(); k++)
      {
	for (int i = 0; i < weights_basis(k).GetM(); i++)
	  C_store_exact(offset + i) /= weights_basis(k)(i);
	
	offset += weights_basis(k).GetM();
      }
    
    if (nb_proc == 1)
      {
        ofstream file_out("PointsInterp.dat");
        file_out.precision(15);
        const typename Dimension::VectR_N& PointsInterp = c_discrete(0)->GetInterpolationPoints();
        for (int i = 0; i < PointsInterp.GetM(); i++)
          {
            PrintNoBrackets(file_out, PointsInterp(i));
            file_out << '\n';
          }
        
        file_out.close();
      }
    
    int nb_pts_quad = 0;
    for (int i = 0; i < points_quad.GetM(); i++)
      if (elt_to_keep(i))
	nb_pts_quad += points_quad(i).GetM();
    
    // computation of CollMatrix that is needed
    // to compute the interpolation of c on the quadrature points
    int nb_basis_fct = c_discrete(0)->GetNbDof();
    CollMatrix.Reallocate(nb_pts_quad, nb_basis_fct);    
    Vector<Real_wp, VectSparse> vec_phi;
    // here CollMatrix(i, j) = phi_j(points_quad(i))
    int irow = 0;
    for (int i = 0; i < points_quad.GetM(); i++)
      if (elt_to_keep(i))
	for (int j = 0; j < points_quad(i).GetM(); j++)
	  {
	    const ElementReference<Dimension, 1>& Fb = var_helm.GetReferenceElementH1(i);
	    const typename Dimension::VectR_N& points_loc = Fb.PointsQuadND();
	    c_discrete(0)->ComputeValuesPhiRef(i, j, points_quad(i)(j),
					       points_loc(j), vec_phi);
	    
	    for (int j = 0; j < vec_phi.GetM(); j++)
	      CollMatrix.Get(irow, vec_phi.Index(j)) = vec_phi.Value(j);
	    
	    irow++;
	  }
    
    var_helm.SetPrintLevel(-1);
  }

  
  //! provides the additional exact values 
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::ProvideExactValues(const VectReal_wp& exact_cte)
  {
    int offset = c_discrete.GetM();
    VectReal_wp all_val;
    for (int k = 0; k < nb_comp_source_invert+nb_angles_to_invert; k++)
      all_val.PushBack(exact_cte(offset++));
    
    C_store_exact.PushBack(all_val);
    
    // C_store_exact.WriteText("Cexact_app.dat");
  }
  
    
  //! finds dofs with Dirichlet condition on c
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::ConstructDirichlet()
  {
    // arrays are allocated
    DofDirichletC.Reallocate(c_discrete.GetM());
    ValueDirichletC.Reallocate(c_discrete.GetM());
    for (int k = 0; k < DofDirichletC.GetM(); k++)
      {
	DofDirichletC(k).Clear();
	ValueDirichletC(k).Clear();
      }
    
    if (dirichlet_c.GetM() <= 0)
      return;
    
    // the input file is read
    for (int k = 0; k < dirichlet_c.GetM(); k++)
      {
	const Vector<string>& param = dirichlet_c(k);
	
	int position_word = -1;
	// we look for the first parameter, with begins by a letter
	// This parameter is the description of the index
	for (int i = 0; i < param.GetM(); i++)
	  if (isalpha(param(i)[0]))
	    {
	      position_word = i;
	      break;
	    }
	
	if (position_word == -1)
	  {
	    cout << "You need to specify the index" << endl;
	    DISP(param);
	    abort();
	  }

	// list of references for this Dirichlet condition
	// is_ref_dirichlet(k) will be true if reference k is a Dirichlet condition
	Vector<bool> is_ref_dirichlet(var_helm.mesh.GetNbReferences()+1);
	is_ref_dirichlet.Fill(false);
	Vector<int> ref(position_word);
	for (int i = 0; i < ref.GetM(); i++)
	  {
	    ref(i) = to_num<int>(param(i));
	    is_ref_dirichlet(ref(i)) = true;
	  }
	
	// we find element numbers and local boundary position
	int num_index = -1;
	for (int p = 0; p < param_c_name.GetM(); p++)
	  if (param_c_name(p) == param(position_word))
	    num_index = p;
	
	if (num_index == -1)
	  {
	    cout << "Unknown parameter " << param(position_word) << endl;
	    cout << "Available parameters : " << param_c_name(k) << endl;
	    abort();
	  }

	SameMeshIndex<Real_wp, Dimension>& c
	  = dynamic_cast<SameMeshIndex<Real_wp, Dimension>& >(*c_discrete(num_index));
	
	Vector<bool> is_dof_dirichlet(c.GetNbDof());
	is_dof_dirichlet.Fill(false);
	int nb_dir = 0;
	for (int i = 0; i < var_helm.mesh.GetNbBoundaryRef(); i++)
	  {
	    int ref = var_helm.mesh.BoundaryRef(i).GetReference();
	    if (is_ref_dirichlet(ref))
	      {
		int num_elem = var_helm.mesh.BoundaryRef(i).numElement(0);
		int num_face = i;
		int num_loc = var_helm.mesh.Element(num_elem).GetPositionBoundary(num_face);
		Vector<int> Nodle = c.GetDofsOnBoundary(num_elem, num_loc);
		for (int j = 0; j < Nodle.GetM(); j++)
		  {
		    if (!is_dof_dirichlet(Nodle(j)))
		      {
			is_dof_dirichlet(Nodle(j)) = true;
			nb_dir++;
		      }
		  }
	      }
	  }
	
	// list of dof numbers
	IVect ListeDof(nb_dir);
	nb_dir = 0;
	for (int i = 0; i < is_dof_dirichlet.GetM(); i++)
	  if (is_dof_dirichlet(i))
	    ListeDof(nb_dir++) = i;
	
	// dof numbers are written
	//ListeDof.WriteText("DofDirichlet" + to_str(k) + ".dat");
	
	Vector<Real_wp> ListeValue(nb_dir);
	if (param(position_word+1) == "UNIFORM")
	  {
	    Real_wp coef = to_num<Real_wp>(param(position_word+2));
	    ListeValue.Fill(coef);
	  }
	else
	  ListeValue.Read(param(position_word+1));
	
	DofDirichletC(num_index).PushBack(ListeDof);
	ValueDirichletC(num_index).PushBack(ListeValue);
      }
  }
  
  
  //! impose fixed values (Dirichlet) for parameters (usually values on the boundary)
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::ImposeDirichletCondition(Vector<VectReal_wp>& cn)
  {
    for (int k = 0; k < cn.GetM(); k++)
      for (int i = 0; i < DofDirichletC(k).GetM(); i++)
	cn(k)(DofDirichletC(k)(i)) = ValueDirichletC(k)(i);        
  }
  
  
  //! initial guess
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::FindInitGuess(VectReal_wp& cn_init)
  {
    // initial guess, cn = 1
    cn_init.Reallocate(C_store_exact.GetM());
    cn_init.Fill(Real_wp(1));

    for (int k = 0; k < min(cn_init.GetM(), C_init_guess.GetM()); k++)
      cn_init(k) = C_init_guess(k);
    
    Vector<VectReal_wp> cn;
    ExtractParametersC(cn_init, cn);
    ImposeDirichletCondition(cn);
    
    for (int k = 0; k < cn.GetM(); k++)
      for (int i = 0; i < weights_basis(k).GetM(); i++)
	cn(k)(i) /= weights_basis(k)(i);

    cn.Nullify();
    
    //cn_init = C_store_exact;
    
    glob_chrono.Start(VirtualTimer::ALL);
  }
  
  
  //! c is given to object var
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::SetVariableIndex(const Vector<VectR_N>& pts_quad, const Vector<VectReal_wp>& eval_c,
		     VarHarmonic_Base<Complex_wp, Dimension>& var)
  {
    int ref_max = var_helm.GetNbPhysicalIndices();
    IVect NbElement_ByRef(ref_max+1);
    NbElement_ByRef.Zero();
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      {
	int ref = var.mesh.Element(i).GetReference();
	NbElement_ByRef(ref)++;
      }
    
    Vector<bool> is_ref_vol(ref_max+1);
    is_ref_vol.Fill(false);
    Vector<PhysicalVaryingMedia<Dimension, Real_wp>* > rho_real;
    Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* > rho_cplx;
    for (int n = 0; n < ref_volume.GetM(); n++)
      {
	int ref = ref_volume(n);
	is_ref_vol(ref) = true;

	rho_real.Clear(); rho_cplx.Clear();
	this->GetVariableIndex(ref, var, rho_real, rho_cplx);
	
	for (int i = 0; i < rho_real.GetM(); i++)
	  rho_real(i)->Reallocate(NbElement_ByRef(ref), false);

	for (int i = 0; i < rho_cplx.GetM(); i++)
	  rho_cplx(i)->Reallocate(NbElement_ByRef(ref), false);
      }
    
    int num = 0;
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      {
	int ref = var.mesh.Element(i).GetReference();
	if (ref <= ref_max)
	  if (is_ref_vol(ref))
	    {
	      rho_real.Clear(); rho_cplx.Clear();
	      this->GetVariableIndex(ref, var, rho_real, rho_cplx);
	      
	      for (int n = 0; n < rho_real.GetM(); n++)
		rho_real(n)->ReallocateOnElement(var.ElementRho(i), pts_quad(i).GetM(), false);

	      for (int n = 0; n < rho_cplx.GetM(); n++)
		rho_cplx(n)->ReallocateOnElement(var.ElementRho(i), pts_quad(i).GetM(), false);
	      
	      for (int j = 0; j < pts_quad(i).GetM(); j++)
		this->SetVariableIndexPoint(ref, i, var.ElementRho(i), j, num, eval_c,
					    rho_real, rho_cplx);
	      
	      num += pts_quad(i).GetM();
	    }	
      }    
  }

  
  //! modifies physical indexes with values contained in cn0
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::FillVariableC(const VectReal_wp& cn0)
  {
    Vector<VectReal_wp> cn;
    CopyParametersC(cn0, cn);
    for (int p = 0; p < cn.GetM(); p++)
      for (int i = 0; i < cn(p).GetM(); i++)
	cn(p)(i) *= weights_basis(p)(i);    
    
    ImposeDirichletCondition(cn);

    // c is evaluated on quadrature points
    eval_c_quad.Reallocate(c_discrete.GetM());
    for (int p = 0; p < c_discrete.GetM(); p++)
      {
	c_discrete(p)->Init(cn(p));
	
	eval_c_quad(p).Reallocate(CollMatrix.GetM());
	Mlt(CollMatrix, cn(p), eval_c_quad(p));
      }
    
    SetVariableIndex(this->points_quad, eval_c_quad, var_helm);
  }
  
  
  //! modifies the source with parameters contained in cn0
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::ModifySource(const VarHarmonic_Base<Complex_wp, Dimension>& var,
                 const VectReal_wp& cn0,
                 Matrix<Complex_wp, General, ColMajor>& source)
  {
    int offset = 0;
    for (int k = 0; k < c_discrete.GetM(); k++)
      offset += c_discrete(k)->GetNbDof();
    
    for (int k = nb_comp_source_invert-1; k >= 0; k--)
      {
        Real_wp Fk = cn0(offset + k);
        int Nvol = var.offset_dof_unknown(1);
        for (int i = 0; i < Nvol; i++)
          for (int j = 0; j < source.GetN(); j++)
            source(i + Nvol*k, j) = source(i, j)*Fk; 
      }
  }

  
  //! computes coefficient used to combine measurements
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::InitCombinationMeasure(const VectReal_wp& cn0)
  {
    int offset = 0;
    for (int k = 0; k < c_discrete.GetM(); k++)
      offset += c_discrete(k)->GetNbDof();
    
    offset += nb_comp_source_invert;
    if (nb_angles_to_invert == 1)
      {
        coef_combination_measure.Reallocate(2);
        Real_wp teta = cn0(offset)*pi_wp/180;
        coef_combination_measure(0) = cos(teta);
        coef_combination_measure(1) = sin(teta);
      }
  }
  

  //! computes the solution for a given set of parameters cn
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::Init(const VectReal_wp& cn0)
  {
    FillVariableC(cn0);
    InitCombinationMeasure(cn0);
    
    // Finite element matrix is computed and factorized
    GlobalGenericMatrix<Complex_wp> nat_mat;
    sol_u.Reallocate(this->omega.GetM());
    u_observable.Reallocate(this->omega.GetM());
    for (int k = 0; k < this->omega.GetM(); k++)
      {
	var_helm.SetOmega(this->omega(k));
	
	// Computing geometry quantities
	var_helm.ComputeMassMatrix(false);

	// LU factorization is updated
	if (this->do_not_store_facto)
	  {
	    if ((k == 1) && (this->same_pattern_LU))
	      linear_solver(0)->UseOldPattern(true);
	    
	    linear_solver(0)->PerformFactorizationStep(nat_mat);
	  }
	else
	  linear_solver(k)->PerformFactorizationStep(nat_mat);
        
	// the solution is computed
	if (this->do_not_store_facto)
	  {
	    sol_u(0) = source_rhs;
            if (nb_comp_source_invert > 0)
              this->ModifySource(var_helm, cn0, sol_u(0));

	    linear_solver(0)->ComputeSolution(sol_u(0), nat_mat);
	  }
	else
	  {
	    sol_u(k) = source_rhs;
            if (nb_comp_source_invert > 0)
              this->ModifySource(var_helm, cn0, sol_u(k));
            
	    linear_solver(k)->ComputeSolution(sol_u(k), nat_mat);
	  }
        
	// interpolation on measure points
	if (this->do_not_store_facto)
	  var_migration.ProjectSimulationData(var_helm, sol_u(0), u_observable(k));
	else
	  var_migration.ProjectSimulationData(var_helm, sol_u(k), u_observable(k));
	
	// multiplication by weights
        if (coef_input_frequency.GetM() > 0)
          u_observable(k) *= coef_input_frequency(k);
	
        ScaleLeftMatrix(u_observable(k), weights_measure);

	number_facto_lu++; number_solve += sol_u(0).GetN();
      }
  }
  

  //! computes the observable vector
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::Evaluate(VectReal_wp& y)
  {
    y.Reallocate(nb_observables);
    int num = 0;
    for (int k = 0; k < this->omega.GetM(); k++)
      AppendData(u_observable(k), num, y);
  }

  
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::MltParamDeriv(const Real_wp& alpha, const Vector<VectReal_wp>& c_quad,
		  Vector<Matrix<Complex_wp, General, ColMajor> >& y)
  {
    int nb_rhs = sol_u(0).GetN();
    y.Reallocate(this->omega.GetM());
    VectComplex_wp u_glob(var_helm.GetNbDof()), y_col(var_helm.GetNbDof());
    for (int k = 0; k < this->omega.GetM(); k++)
      {
	y(k).Reallocate(var_helm.GetNbDof(), nb_rhs);
	for (int n = 0; n < nb_rhs; n++)
	  {
	    y_col.Zero();
	    GetCol(sol_u(k), n, u_glob);
	    int offset = 0;
	    Real_wp omega2 = this->omega(k)*this->omega(k);
	    for (int i = 0; i < var_helm.mesh.GetNbElt(); i++)
	      if (elt_to_keep(i))
		{
		  this->MltParamDerivLocal(k, omega2, u_glob, i, 
					   alpha, c_quad, offset, y_col);
		  
		  offset += points_quad(i).GetM();
		}
	    
	    var_helm.ImposeNullDirichletCondition(y_col);
	    SetCol(y_col, n, y(k));
	  }
      }
  }


  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::MltParamDerivOmega(const Real_wp& alpha, int k, const Vector<VectReal_wp>& c_quad,
		       const Matrix<Complex_wp, General, ColMajor>& x,
		       Matrix<Complex_wp, General, ColMajor>& y)
  {
    int nb_rhs = source_rhs.GetN();
    VectComplex_wp u_glob(var_helm.GetNbDof()), y_col(var_helm.GetNbDof());
    y.Reallocate(var_helm.GetNbDof(), nb_rhs);
    for (int n = 0; n < nb_rhs; n++)
      {
	y_col.Zero();
	GetCol(x, n, u_glob);
	int offset = 0;
	Real_wp omega2 = this->omega(k)*this->omega(k);
	for (int i = 0; i < var_helm.mesh.GetNbElt(); i++)
	  if (elt_to_keep(i))
	    {
	      this->MltParamDerivLocal(k, omega2, u_glob, i, 
				       alpha, c_quad, offset, y_col);
		  
	      offset += points_quad(i).GetM();
	    }
	
	var_helm.ImposeNullDirichletCondition(y_col);
	SetCol(y_col, n, y);
      }
  }
  
  
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::MltParamDerivTrans(const Real_wp& alpha, const Vector<Matrix<Complex_wp, General, ColMajor> >& y,
		       Vector<VectReal_wp>& ProdC)
  {
    int nb_pts_quad = 0;
    for (int i = 0; i < points_quad.GetM(); i++)
      if (elt_to_keep(i))
	nb_pts_quad += points_quad(i).GetM();
    
    ProdC.Reallocate(c_discrete.GetM());
    for (int p = 0; p < c_discrete.GetM(); p++)
      {
	ProdC(p).Reallocate(nb_pts_quad);
	ProdC(p).Zero();
      }
    
    int nb_rhs = sol_u(0).GetN();
    VectComplex_wp u_glob(var_helm.GetNbDof()), y_col(var_helm.GetNbDof());
    for (int k = 0; k < this->omega.GetM(); k++)
      {
	for (int n = 0; n < nb_rhs; n++)
	  {
	    GetCol(y(k), n, y_col);
	    var_helm.ImposeNullDirichletCondition(y_col);
	    
	    GetCol(sol_u(k), n, u_glob);
	    int offset = 0;
	    Real_wp omega2 = this->omega(k)*this->omega(k);
	    for (int i = 0; i < var_helm.mesh.GetNbElt(); i++)
	      if (elt_to_keep(i))
		{
		  this->MltParamDerivTransLoc(k, omega2, u_glob, i,
					      alpha, offset, y_col, ProdC);
		  
		  offset += points_quad(i).GetM();
		}
	  }
      }
  }
  
  
  //! Applies jacobian matrix of F to vector h, res = DF(cn) h
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::Derivative(const VectReal_wp& h0, VectReal_wp& res)
  {
    if (this->do_not_store_facto)
      {
	cout << "Not implemented " << endl;
	abort();
      }
    
    Vector<VectReal_wp> h;
    CopyParametersC(h0, h);

    Vector<VectReal_wp> rhs(h.GetM());
    Vector<Matrix<Complex_wp, General, ColMajor> > rhs_cplx;
    for (int num_index = 0; num_index < h.GetM(); num_index++)
      {
	for (int i = 0; i < h(num_index).GetM(); i++)
	  h(num_index)(i) *= weights_basis(num_index)(i);
	
	for (int i = 0; i < DofDirichletC(num_index).GetM(); i++)
	  SetComplexZero(h(num_index)(DofDirichletC(num_index)(i)));
	
	// evaluation of h on quadrature points
	rhs(num_index).Reallocate(CollMatrix.GetM());
	Mlt(CollMatrix, h(num_index), rhs(num_index));
      }
    
    // multiplication by derivative of parameters
    MltParamDeriv(Real_wp(1), rhs, rhs_cplx);
    
    // we add terms due to source parameters
    for (int k = 0; k < this->omega.GetM(); k++)
      {
        int offset = 0;
        for (int p = 0; p < c_discrete.GetM(); p++)
          offset += c_discrete(p)->GetNbDof();
        
        for (int m = nb_comp_source_invert-1; m >= 0; m--)
          {
            Real_wp hm = h0(offset + m);
            int Nvol = var_helm.offset_dof_unknown(1);
            for (int i = 0; i < Nvol; i++)
              for (int j = 0; j < rhs_cplx(k).GetN(); j++)
                rhs_cplx(k)(i + Nvol*m, j) += source_rhs(i, j)*hm;
          }
      }

    res.Reallocate(nb_observables);
    res.Zero();

    // solving linear systems
    GlobalGenericMatrix<Complex_wp> nat_mat;
    Matrix<Complex_wp, General, ColMajor> data_simu, data_comb;
    int num = 0;
    for (int k = 0; k < this->omega.GetM(); k++)
      {
	linear_solver(k)->ComputeSolution(rhs_cplx(k), nat_mat);
	
	// value on quadrature points
	var_migration.ProjectSimulationData(var_helm, rhs_cplx(k), data_simu);

	// multiplication by weights
        if (coef_input_frequency.GetM() > 0)
          data_simu *= coef_input_frequency(k);
	
        ScaleLeftMatrix(data_simu, weights_measure);
	
        // result is appended in res
        CombineDeriveData(h0, k, data_simu, data_comb);
        ScatterData(k, data_comb, num, res);
        
	number_solve += rhs_cplx(k).GetN();
      }
  }

  
  //! Applies adjoint of jacobian matrix of F to vector G, res = DF*(cn) h
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::Adjoint(const VectReal_wp& G, VectReal_wp& res)
  {
    if (this->do_not_store_facto)
      {
	cout << "Not implemented " << endl;
	abort();
      }
    
    // solving linear systems
    Vector<Matrix<Complex_wp, General, ColMajor> > rhs_cplx(this->omega.GetM());
    GlobalGenericMatrix<Complex_wp> nat_mat;
    Matrix<Complex_wp, General, ColMajor> data_simu, data_comb;
    int nb_pts = var_migration.GetNbPointsQuadrature();
    int nb_rhs = sol_u(0).GetN();
    int nb_unknowns_measured = var_helm.nb_unknowns_scal - nb_angles_measure;
    data_comb.Reallocate(nb_pts, nb_rhs*nb_unknowns_measured);
    data_simu.Reallocate(nb_pts, nb_rhs*var_helm.nb_unknowns_scal);
    int num = 0;
    res.Reallocate(C_store_exact.GetM());
    for (int k = 0; k < this->omega.GetM(); k++)
      {
	ExtractData(k, G, num, data_comb);
        CombineAdjointData(k, data_comb, data_simu, res);
	Conjugate(data_simu);
	
        if (coef_input_frequency.GetM() > 0)
          data_simu *= coef_input_frequency(k);
        
	for (int i = 0; i < data_simu.GetM(); i++)
	  for (int j = 0; j < data_simu.GetN(); j++)
	    data_simu(i, j) /= weights_measure(i);
	
	var_migration.ComputeIntegralResidu(var_helm, data_simu, rhs_cplx(k));
	linear_solver(k)->ComputeSolution(SeldonTrans, rhs_cplx(k), nat_mat);
	Conjugate(rhs_cplx(k));
	
	number_solve += rhs_cplx(k).GetN();
      }
    
    // multiplication by mass matrix
    Vector<VectReal_wp> aux;
    Vector<VectReal_wp> res_n;
    ExtractParametersC(res, res_n);
    
    // we add terms due to source parameters
    for (int k = 0; k < this->omega.GetM(); k++)
      {
        int offset = 0;
        for (int p = 0; p < c_discrete.GetM(); p++)
          offset += c_discrete(p)->GetNbDof();
        
        for (int m = nb_comp_source_invert-1; m >= 0; m--)
          {
            int Nvol = var_helm.offset_dof_unknown(1);
            for (int i = 0; i < Nvol; i++)
              for (int j = 0; j < rhs_cplx(k).GetN(); j++)
                res(offset+m) += realpart(rhs_cplx(k)(i + Nvol*m, j))*realpart(source_rhs(i, j))
                  + imagpart(rhs_cplx(k)(i + Nvol*m, j))*imagpart(source_rhs(i, j));
          }
      }
    
    // transpose of MltParamDeriv operation
    MltParamDerivTrans(Real_wp(1), rhs_cplx, aux);

    for (int num_index = 0; num_index < aux.GetM(); num_index++)
      {
	// transpose of CollMatrix
	Mlt(SeldonTrans, CollMatrix, aux(num_index), res_n(num_index));
	
	for (int i = 0; i < res_n(num_index).GetM(); i++)
	  res_n(num_index)(i) *= weights_basis(num_index)(i);

	for (int i = 0; i < DofDirichletC(num_index).GetM(); i++)
	  SetComplexZero(res_n(num_index)(DofDirichletC(num_index)(i)));    
      }
    
    res_n.Nullify();
  }


  //! Computes both f(x) and df(x)
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::EvaluateJacobian(const VectReal_wp& cn0, VectReal_wp& f, Matrix<Real_wp>& df)
  {
    int nb_ddl = cn0.GetM();
    FillVariableC(cn0);
    InitCombinationMeasure(cn0);
    
    VectReal_wp rhs, zero_vec; Vector<VectReal_wp> rhs_vec(c_discrete.GetM());
    VectComplex_wp df_col; Vector<Real_wp, VectSparse> rhs_sparse;
    Matrix<Complex_wp, General, ColMajor> rhs_cplx, rhs_store;
        
    GlobalGenericMatrix<Complex_wp> nat_mat;
    u_observable.Reallocate(this->omega.GetM());

    df.Reallocate(nb_observables, nb_ddl);
    df.Zero();
    
    // main loop on frequencies
    Matrix<Complex_wp, General, ColMajor> data_simu, x_sol, data_comb;
    Vector<int> num(nb_ddl); num.Zero();
    for (int k = 0; k < this->omega.GetM(); k++)
      {
	var_helm.SetOmega(this->omega(k));
	
	// Computing geometry quantities
	var_helm.ComputeMassMatrix(false);

	// LU factorization is updated
	if ((k == 1) && (this->same_pattern_LU))
	  linear_solver(0)->UseOldPattern(true);
	
	linear_solver(0)->PerformFactorizationStep(nat_mat);

	// the solution is computed
	x_sol = source_rhs;
        if (nb_comp_source_invert > 0)
          this->ModifySource(var_helm, cn0, x_sol);
        
	linear_solver(0)->ComputeSolution(x_sol, nat_mat);
	number_facto_lu++; number_solve += x_sol.GetN();

	// observation is computed
	var_migration.ProjectSimulationData(var_helm, x_sol, u_observable(k));
	
	// multiplication by weights
        if (coef_input_frequency.GetM() > 0)
          u_observable(k) *= coef_input_frequency(k);

	ScaleLeftMatrix(u_observable(k), weights_measure);

	// loop over columns of df
	int nb_rhs = 1, num_ddl = 0;
	for (int num_index = 0; num_index < weights_basis.GetM(); num_index++)
	  for (int k2 = 0; k2 < weights_basis(num_index).GetM(); k2++)
	    {
	      GetCol(CollMatrix, k2, rhs_sparse);
	      rhs.Reallocate(CollMatrix.GetM());
	      rhs.Zero();
	      for (int j = 0; j < rhs_sparse.GetM(); j++)
		rhs(rhs_sparse.Index(j)) = rhs_sparse.Value(j);
	      
	      rhs *= weights_basis(num_index)(k2);

	      zero_vec.Reallocate(rhs.GetM()); zero_vec.Zero();
	      
	      rhs_vec.Reallocate(weights_basis.GetM());
	      rhs_vec.Fill(zero_vec);
	      
	      rhs_vec(num_index) = rhs;
	      
	      MltParamDerivOmega(Real_wp(1), k, rhs_vec, x_sol, rhs_cplx);
	      nb_rhs = rhs_cplx.GetN();

	      if (rhs_store.GetM() == 0)
		rhs_store.Reallocate(rhs_cplx.GetM(), nb_rhs*nb_ddl);

	      for (int j = 0; j < nb_rhs; j++)
		{
		  GetCol(rhs_cplx, j, df_col);
		  SetCol(df_col, nb_rhs*num_ddl+j, rhs_store); 
		}

	      num_ddl++;
	    }

        // we add terms due to source parameters
        int offset = 0;
        for (int p = 0; p < c_discrete.GetM(); p++)
          offset += c_discrete(p)->GetNbDof();
        
        for (int m = nb_comp_source_invert-1; m >= 0; m--)
          {
            int Nvol = var_helm.offset_dof_unknown(1);
            for (int i = 0; i < Nvol; i++)
              for (int j = 0; j < nb_rhs; j++)
                rhs_store(i + Nvol*m, nb_rhs*(offset+m) + j) = source_rhs(i, j);
          }
        
	// we solve all needed right hand sides
	linear_solver(0)->ComputeSolution(rhs_store, nat_mat);
	number_solve += rhs_store.GetN();

	var_migration.ProjectSimulationData(var_helm, rhs_store, data_simu);

        if (coef_input_frequency.GetM() > 0)
          data_simu *= coef_input_frequency(k);
        
	ScaleLeftMatrix(data_simu, weights_measure);

	// df is updated
        if (this->nb_angles_measure == 1)
          {
            offset += this->nb_comp_source_invert;
            Real_wp coef_teta = pi_wp/180;
            Real_wp c1 = this->coef_combination_measure(0);
            Real_wp c2 = this->coef_combination_measure(1);
            data_comb.Reallocate(data_simu.GetM(), data_simu.GetN()/2);
            for (int i = 0; i < data_comb.GetM(); i++)
              for (int j = 0; j < data_comb.GetN(); j++)
                data_comb(i, j) = data_simu(i, 2*j)*c1 + data_simu(i, 2*j+1)*c2;

            if (this->nb_angles_to_invert == 1)
              for (int i = 0; i < data_comb.GetM(); i++)
                for (int j = 0; j < nb_rhs; j++)
                  {
                    Complex_wp ux = u_observable(k)(i, 2*j), uy = u_observable(k)(i, 2*j+1);
                    data_comb(i, nb_rhs*offset + j) = (-c2*ux + c1*uy)*coef_teta;
                  }
          }
        else
          data_comb = data_simu;
        
	int nb_u = data_comb.GetN() / rhs_store.GetN();
	for (int i = 0; i < data_comb.GetN(); i++)
	  for (int j = 0; j < data_comb.GetM(); j++)
	    {
	      int nb_rhs_d = nb_rhs*nb_u;
	      int col = i/nb_rhs_d;
	      switch (measured_part)
		{
		case ALL_PART:
		  df(num(col), col) = realpart(data_comb(j, i));
		  df(num(col)+1, col) = imagpart(data_comb(j, i));
		  num(col) += 2;
		  break;
		case REAL_PART:
		  df(num(col), col) = realpart(data_comb(j, i));
		  num(col)++;
		  break;
		case IMAG_PART:
		  df(num(col), col) = imagpart(data_comb(j, i));
		  num(col)++;
		  break;
		case MODULUS_PART:
		  {
		    Real_wp ur = realpart(u_observable(k)(j, i%nb_rhs_d)),
		      ui = imagpart(u_observable(k)(j, i%nb_rhs_d));
		    
		    Real_wp dur = realpart(data_comb(j, i)), dui = imagpart(data_comb(j, i));
		    df(num(col), col) = (dur*ur + dui*ui) / abs(u_observable(k)(j, i%nb_rhs_d));
		    num(col)++;
		  }
		}

	      col++;
	    }
      }
    
    Evaluate(f);
  }
  
  
  //! Computes the exact measures y, parameter c will
  //! solve the non-linear least-squares problem F(c) = y
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::EvaluateDifferently(VectReal_wp& y)
  {
    // object related to experimental data is constructed
    VarHarmonic_Base<Complex_wp, Dimension>* var_exp_ptr;
    var_exp_ptr = dynamic_cast<VarHarmonic_Base<Complex_wp, Dimension>* >(var_helm.GetNewEllipticProblem());
    VarHarmonic_Base<Complex_wp, Dimension>& var_exp = *var_exp_ptr;    
    All_LinearSolver* solver_exp;

    int rank_proc = 0;
    
#ifdef SELDON_WITH_MPI
    int nb_proc = 1;
    var_exp.comm_group_mode = var_helm.comm_group_mode;
    MPI_Comm_rank(this->comm_freq, &rank_proc);
    MPI_Comm_size(this->comm_freq, &nb_proc);
#endif

    // y will store the solution evaluated with exact c
    bool experimental_data = false;
    if (GetExtension(var_migration.input_file_exp) != "ini")
      {
        experimental_data = true;
        y.Read(var_migration.input_file_exp);
      }
    else
      {
        var_exp.KeepGridLocalizationArrays();
        var_exp.ConstructAll(var_migration.input_file_exp, name_element, name_equation,
                             solver_exp, true, false);

        var_exp.SetPrintLevel(-1);
        solver_exp->SetPrintLevel(-1);    
      }

    // quadrature points on the surface where du/dn is measured are constructed
    var_migration.ComputeQuadrature(var_helm, var_exp);
    
    var_migration.GetWeightsQuadrature(weights_measure);
    
    for (int i = 0; i < weights_measure.GetM(); i++)
      weights_measure(i) = sqrt(weights_measure(i));
    
    // number of observables
    int nb_rhs = source_rhs.GetN(), nb_points = weights_measure.GetM();
    int nb_unknowns_measured = var_helm.nb_unknowns_scal - nb_angles_measure;
    nb_observables = this->omega.GetM()*nb_rhs*nb_points*nb_unknowns_measured;
    if (measured_part == ALL_PART)
      nb_observables *= 2;
        	
    if (experimental_data)
      {
#ifdef SELDON_WITH_MPI
        if (nb_proc > 1)
          {
            VectReal_wp y0(y);
            y.Resize(nb_observables);
            int offset = this->offset_omega*nb_rhs*nb_points*nb_unknowns_measured;
            if (measured_part == ALL_PART)
              offset *= 2;

            for (int i = 0; i < nb_observables; i++)
              y(i) = y0(offset+i);
          }
#endif

        if (y.GetM() != nb_observables)
          {
            cout << "The experiment data file " << var_migration.input_file_exp
                 << " does contain " << y.GetM() << " values " << endl;
            cout << "But " << nb_observables << " are expected " << endl;
            abort();
          }
        
        delete var_exp_ptr;
        return;
      }

    int ref_max = this->ref_volume.GetNormInf();
    for (int i = 0; i < var_exp.mesh.GetNbElt(); i++)
      ref_max = max(ref_max, var_exp.mesh.Element(i).GetReference());
    
    Vector<bool> is_ref_kept(ref_max+1);
    is_ref_kept.Fill(false);
    for (int i = 0; i < ref_volume.GetM(); i++)
      is_ref_kept(ref_volume(i)) = true;

    Vector<bool> elt_to_keep_ref(var_exp.mesh.GetNbElt());
    elt_to_keep_ref.Fill(false);
    for (int i = 0; i < var_exp.mesh.GetNbElt(); i++)
      if (is_ref_kept(var_exp.mesh.Element(i).GetReference()))
	elt_to_keep_ref(i) = true;
    
    // evaluating the exact parameters c on quadrature points
    int nb_pts = 0;
    Vector<typename Dimension::VectR_N>& pts_quad = var_exp.Glob_PointsQuadrature;
    for (int i = 0; i < pts_quad.GetM(); i++)
      if (elt_to_keep_ref(i))
	nb_pts += pts_quad(i).GetM();

    int num = 0;    
    Vector<VectReal_wp> eval_c(c_exact.GetM());
    for (int p = 0; p < c_exact.GetM(); p++)
      {
	eval_c(p).Reallocate(nb_pts);
	num = 0;
	for (int i = 0; i < pts_quad.GetM(); i++)
	  if (elt_to_keep_ref(i))
	    for (int j = 0; j < pts_quad(i).GetM(); j++)
	      eval_c(p)(num++) = (*c_exact(p))(pts_quad(i)(j));
      }
    
    SetVariableIndex(pts_quad, eval_c, var_exp);    
    
    // vector containing measures is allocated
    y.Reallocate(nb_observables);
    y.Zero();

    // the source is computed
    BoundaryInverseProblemMultiFreqParam<Dimension> param_freq(*this, y);
    var_migration.ComputeRightHandSide(var_exp, param_freq.source_rhs_ref);
    
    // modification of the source with exact parameters (if required)
    if (nb_comp_source_invert > 0)
      this->ModifySource(var_exp, C_store_exact, param_freq.source_rhs_ref);

    DISP(MaxAbs(param_freq.source_rhs_ref));
    InitCombinationMeasure(C_store_exact);
    
    // all the frequencies are solved
    this->SolveAllFreq(var_exp, solver_exp, param_freq);
    
    if (write_experimental_data)
      WriteMeasuredY(y, file_name_experimental_data);
    
    delete solver_exp; delete var_exp_ptr;
    
    if (rank_proc == 0)
      cout << "Interpolation error when projecting c_exact on finite element mesh = "
	   << this->GetInterpolationError() << endl;
    
    if (NdisplayC == -1)
      {
        if (rank_proc == 0)
	  cout << "Computation ended successfully" << endl;
	
#ifdef SELDON_WITH_MPI
        MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
        
        exit(0);
      }
  }

  
  //! Adds observables to y 
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::AppendData(const Matrix<Complex_wp, General, ColMajor>& data0,
	       int& num, VectReal_wp& y)
  {
    // linear combination if required
    Matrix<Complex_wp, General, ColMajor> data;
    CombineData(data0, data);
    
    // then we extract the measured part
    for (int i = 0; i < data.GetN(); i++)
      for (int j = 0; j < data.GetM(); j++)
	{
	  switch (measured_part)
            {
            case ALL_PART:
              y(num) = realpart(data(j, i));
              y(num+1) = imagpart(data(j, i));
              num += 2;
              break;
            case REAL_PART:
              y(num) = realpart(data(j, i)); num++;
              break;
            case IMAG_PART:
              y(num) = imagpart(data(j, i)); num++;
              break;
            case MODULUS_PART:
              y(num) = abs(data(j, i)); num++;
              break;
            }
	}
  }

  
  //! linear combination of data
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::CombineData(const Matrix<Complex_wp, General, ColMajor>& data0,
                Matrix<Complex_wp, General, ColMajor>& data)
  {
    // default choice : no combination
    data = data0;
  }
  

  //! linear combination of data (and derivative)
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::CombineDeriveData(const VectReal_wp& h0, int k,
                      const Matrix<Complex_wp, General, ColMajor>& data_simu,
                      Matrix<Complex_wp, General, ColMajor>& data_comb)
  {
    // default choice : no combination
    data_comb = data_simu;
  }


  //! linear combination of data (and adjoint of derivative)
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::CombineAdjointData(int k, const Matrix<Complex_wp, General, ColMajor>& data_simu,
                       Matrix<Complex_wp, General, ColMajor>& data_comb, VectReal_wp& h0)
  {
    // default choice : no combination
    data_comb = data_simu;
  }

  
  //! Adds derivatives of observables to y
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::ScatterData(int k, const Matrix<Complex_wp, General, ColMajor>& data_simu,
                int& num, VectReal_wp& res)
  {
    for (int i = 0; i < data_simu.GetN(); i++)
      for (int j = 0; j < data_simu.GetM(); j++)
        {
          switch (measured_part)
            {
            case ALL_PART:
              res(num) = realpart(data_simu(j, i));
              res(num+1) = imagpart(data_simu(j, i)); num += 2;
              break;
            case REAL_PART:
              res(num) = realpart(data_simu(j, i)); num++;
              break;
            case IMAG_PART:
              res(num) = imagpart(data_simu(j, i)); num++;
              break;
            case MODULUS_PART:
              Real_wp ur = realpart(u_observable(k)(j, i)),
                ui = imagpart(u_observable(k)(j, i));
              
              Real_wp dur = realpart(data_simu(j, i)), dui = imagpart(data_simu(j, i));
              res(num) = (dur*ur + dui*ui) / abs(u_observable(k)(j, i));
              num++;
            }
        }
  }
  

  //! adjoint of ScatterData
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>
  ::ExtractData(int k, const VectReal_wp& y,
		int& num, Matrix<Complex_wp, General, ColMajor>& data)
  {
    for (int i = 0; i < data.GetN(); i++)
      for (int j = 0; j < data.GetM(); j++)
	{
	  switch (measured_part)
            {
            case ALL_PART:
              data(j, i) = Complex_wp(y(num), y(num+1));
              num += 2;
              break;
            case REAL_PART:
              data(j, i) = y(num); num++;
              break;
            case IMAG_PART:
              data(j, i) = Complex_wp(0, y(num)); num++;
              break;
            case MODULUS_PART:
              Real_wp ur = realpart(u_observable(k)(j, i)),
                ui = imagpart(u_observable(k)(j, i));
              
              Real_wp moduleU = abs(u_observable(k)(j, i));
              data(j, i) = Complex_wp(ur, ui)/moduleU*y(num);
              num++;
              break;
            }
	}
  }


  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::WriteMeasuredY(const VectReal_wp& y, const string& output_name)
  {
    int nb_proc = 1;
    
#ifdef SELDON_WITH_MPI
    int rank_proc = 0;
    MPI_Comm_rank(comm_freq, &rank_proc);
    MPI_Comm_size(comm_freq, &nb_proc);
#endif

    if (nb_proc == 1)
      {
        y.Write(output_name);
        return;
      }
    
#ifdef SELDON_WITH_MPI
    IVect nb_observ_per_proc;
    if (rank_proc == 0)
      nb_observ_per_proc.Reallocate(nb_proc);
    
    MPI_Gather(&nb_observables, 1, MPI_INTEGER, nb_observ_per_proc.GetData(), 1, MPI_INTEGER, 0, comm_freq);
    
    Vector<int64_t> xtmp;
    MPI_Status status;
    if (rank_proc == 0)
      {
        VectReal_wp y_all, y_proc;
        int nb_all = 0;
        for (int i = 0; i < nb_proc; i++)
          nb_all += nb_observ_per_proc(i);
        
        y_all.Reallocate(nb_all);
        int offset = 0;
        for (int i = 0; i < nb_proc; i++)
          if (nb_observ_per_proc(i) > 0)
            {
              y_proc.Reallocate(nb_observ_per_proc(i));
              if (i == 0)
                y_proc = y;
              else
                MpiRecv(comm_freq, y_proc, xtmp, y_proc.GetM(), i, 17, status);
              
              for (int k = 0; k < nb_observ_per_proc(i); k++)
                y_all(offset + k) = y_proc(k);
              
              offset += nb_observ_per_proc(i);
            }
        
        y_all.Write(output_name);
      }
    else
      {
        if (nb_observables > 0)
          MpiSsend(comm_freq, const_cast<VectReal_wp&>(y), xtmp, y.GetM(), 0, 17);
      }
#endif
  }
  
  
  template<class Dimension>
  void BoundaryInverseProblem_Base<Dimension>::GiveIterate(int n, const VectReal_wp& x)
  {
    if (n%NdisplayC == 0)
      {
	WriteDatas(x);
	if (write_experimental_data)
	  {
	    VectReal_wp y;
	    this->Evaluate(y);
            this->WriteMeasuredY(y, "y_n"+to_str(n)+".dat");
          }
      }
  }


  /****************************************
   * BoundaryInverseProblemMultiFreqParam *
   ****************************************/
  
  
  //! constructor with a given problem
  template<class Dimension>
  BoundaryInverseProblemMultiFreqParam<Dimension>
  ::BoundaryInverseProblemMultiFreqParam(BoundaryInverseProblem_Base<Dimension>& var, VectReal_wp& y0)
    : var_migration(var.GetMigrationProblem()), var_inverse(var), y(y0)
  {
    num = 0;
  }
  
  
  //! fills the source for the local frequency k
  template<class Dimension>
  void BoundaryInverseProblemMultiFreqParam<Dimension>
  ::SetSourceFrequency(int k, int global_k, int m, int mg, const VarComputationProblem& var_exp,
                       Matrix<Complex_wp, General, ColMajor>& sol_exp)
  {
    sol_exp = source_rhs_ref;
  }

  
  //! treats the solution for the local frequency k
  template<class Dimension>
  void BoundaryInverseProblemMultiFreqParam<Dimension>
  ::GiveSolutionFrequency(int k, int global_k, int m, int mg, 
                          const VarComputationProblem& var_exp,
                          Matrix<Complex_wp, General, ColMajor>& sol_exp)
  {
    VarHarmonic_Base<Complex_wp, Dimension>& var_harmonic
      = const_cast<VarHarmonic_Base<Complex_wp, Dimension>& >(dynamic_cast<const VarHarmonic_Base<Complex_wp, Dimension>& >(var_exp));
    
    // values on quadrature points are retrieved
    var_migration.ProjectExperimentData(var_harmonic, sol_exp, data_exp);
    
    // multiplication by weights
    if (var_inverse.coef_input_frequency.GetM() > 0)
      data_exp *= var_inverse.coef_input_frequency(k);
    
    ScaleLeftMatrix(data_exp, var_inverse.weights_measure);
    
    // measured data is appended to vector y
    var_inverse.AppendData(data_exp, num, y);
    
    const VarProblem_Base& var_problem = var_exp.GetVarProblemBase();
    const VarOutputProblem_Base& var_output = dynamic_cast<const VarOutputProblem_Base&>(var_exp);
    if (var_inverse.write_experimental_data)
      {	    
        VectComplex_wp x_sol(var_exp.GetNbRows());
        for (int n = 0; n < sol_exp.GetN(); n++)
          {
            GetCol(sol_exp, n, x_sol);
            var_output.WriteOutputFile(x_sol, "SolOm" + to_str(global_k) + "P" + to_str(n) + "_G", var_problem.nb_unknowns_scal);
          }
      }  
  }
  
  
  /*********************************************************
   * BoundaryInverseProblem<HelmholtzEquation<Dimension> > *
   *********************************************************/

  
  //! Default constructor
  template<class Dimension>
  BoundaryInverseProblem<HelmholtzEquation<Dimension> >::BoundaryInverseProblem()
    : BoundaryInverseProblem_Base<Dimension>(var_helm, var_migration)
  {
    type_convention = HELIO;
  }
  
  
  //! sets the parameters to be inverted
  template<class Dimension>
  void BoundaryInverseProblem<HelmholtzEquation<Dimension> >::
  SetInversionParameters(const VectString& param)
  {
    if (param.GetM() <= 0)
      {
        cout << "Please provide at least a parameter to invert" << endl;
        abort();
      }

    type_convention = NONE;
    this->ResizeNbParameters(param.GetM());
    for (int k = 0; k < param.GetM(); k++)
      this->SetExactParameter(k, param(k), this->exact_fct_to_find[k]);
  }


  //! computation of density and speed of sound prior to inversion
  template<class Dimension>
  void BoundaryInverseProblem<HelmholtzEquation<Dimension> >::InitIndices()
  {
    ref_density.Reallocate(var_helm.ref_rho.GetM());
    ref_c.Reallocate(var_helm.ref_rho.GetM());
    ref_gamma.Reallocate(var_helm.ref_rho.GetM());
    
    ref_rho.Reallocate(var_helm.ref_rho.GetM());
    ref_mu.Reallocate(var_helm.ref_rho.GetM());
    ref_sigma.Reallocate(var_helm.ref_rho.GetM());
    ref_M.Reallocate(var_helm.ref_rho.GetM());
    ref_beta.Reallocate(var_helm.ref_rho.GetM());
    ref_beta.Zero();
    for (int ref = 0; ref < var_helm.ref_rho.GetM(); ref++)
      {
	Real_wp rho = realpart(var_helm.ref_rho(ref).GetConstant());
	Real_wp mu = realpart(var_helm.ref_mu(ref).GetConstant()(0, 0));
	ref_density(ref) = 1.0/mu;
	ref_c(ref) = sqrt(mu/rho);
        
	if (type_convention == NATURAL)
          {
            ref_rho(ref) = rho;
            ref_mu(ref) = mu;
            ref_sigma(ref) = realpart(var_helm.ref_sigma(ref).GetConstant());
            if (var_helm.include_flow_term)
              for (int m = 0; m < Dimension::dim_N; m++)
                ref_M(ref)(m) = realpart(var_helm.ref_v0(ref).GetConstant()(m));
            
            if (var_helm.include_double_gradM_flow)
              ref_beta(ref) = realpart(var_helm.ref_beta(ref).GetConstant());
          }
        else
          {
            ref_gamma(ref) = realpart(var_helm.ref_sigma(ref).GetConstant()) / (Real_wp(2)*rho);
            if (var_helm.include_flow_term)
              for (int m = 0; m < Dimension::dim_N; m++)
                ref_M(ref)(m) = realpart(var_helm.ref_v0(ref).GetConstant()(m)) / rho;
          }
      }
    
    int ref_max = this->ref_volume.GetNormInf();
    bool varying_rho = false, varying_mu = false, varying_sigma = false;
    bool varying_flow = false, varying_beta = false;
    Vector<bool> is_elt_ref(ref_max+1);
    is_elt_ref.Fill(false);
    for (int kref = 0; kref < this->ref_volume.GetM(); kref++)      
      {
	int ref = this->ref_volume(kref);
	is_elt_ref(ref) = true;
	if (var_helm.ref_rho(ref).IsVarying())
	  varying_rho = true;

	if (var_helm.ref_mu(ref).IsVarying())
	  varying_mu = true;

	if (var_helm.ref_sigma(ref).IsVarying())
	  varying_sigma = true;

	if (var_helm.include_flow_term)
	  if (var_helm.ref_v0(ref).IsVarying())
	    varying_flow = true;
	
	if (var_helm.include_double_gradM_flow)
	  if (var_helm.ref_beta(ref).IsVarying())
	    varying_beta = true;
      }
    
    if (type_convention == NATURAL)
      {
	if (varying_rho)
	  variable_rho.Reallocate(var_helm.mesh.GetNbElt());

	if (varying_mu)
	  variable_mu.Reallocate(var_helm.mesh.GetNbElt());

	if (varying_sigma)
	  variable_sigma.Reallocate(var_helm.mesh.GetNbElt());

	if (varying_flow)
	  variable_M.Reallocate(var_helm.mesh.GetNbElt());

	if (varying_beta)
	  variable_beta.Reallocate(var_helm.mesh.GetNbElt());

	for (int i = 0; i < var_helm.mesh.GetNbElt(); i++)
	  {
	    int ref = var_helm.mesh.Element(i).GetReference();
	    int i1 = var_helm.ElementRho(i);
	    if ((ref < is_elt_ref.GetM()) && (is_elt_ref(ref)))
	      {
		if (varying_rho)
		  {
		    variable_rho(i).Reallocate(var_helm.ref_rho(ref).GetVaryingIndex().GetValue()(i1).GetM());
		    for (int j = 0; j < variable_rho(i).GetM(); j++)
		      variable_rho(i)(j) = realpart(var_helm.ref_rho(ref).GetVaryingIndex().GetValue()(i1)(j));
		  }
		
		if (varying_mu)
		  {
		    variable_mu(i).Reallocate(var_helm.ref_mu(ref).GetVaryingIndex(0, 0).GetValue()(i1).GetM());
		    for (int j = 0; j < variable_mu(i).GetM(); j++)
		      variable_mu(i)(j) = realpart(var_helm.ref_mu(ref).GetVaryingIndex(0, 0).GetValue()(i1)(j));
		  }
		
		if (varying_sigma)
		  {
		    variable_sigma(i).Reallocate(var_helm.ref_sigma(ref).GetVaryingIndex().GetValue()(i1).GetM());
		    for (int j = 0; j < variable_sigma(i).GetM(); j++)
		      variable_sigma(i)(j) = realpart(var_helm.ref_sigma(ref).GetVaryingIndex().GetValue()(i1)(j));
		  }

		if (varying_beta)
		  {
		    variable_beta(i).Reallocate(var_helm.ref_beta(ref).GetVaryingIndex().GetValue()(i1).GetM());
		    for (int j = 0; j < variable_beta(i).GetM(); j++)
		      variable_beta(i)(j) = realpart(var_helm.ref_beta(ref).GetVaryingIndex().GetValue()(i1)(j));
		  }
		if (varying_flow)
		  {
		    int nb_pts = var_helm.ref_v0(ref).GetVaryingIndex(0).GetValue()(i1).GetM();
		    variable_M(i).Reallocate(nb_pts);
		    for (int j = 0; j < nb_pts; j++)
		      for (int d = 0; d < Dimension::dim_N; d++)
			variable_M(i)(j)(d) = realpart(var_helm.ref_v0(ref)
						       .GetVaryingIndex(d).GetValue()(i1)(j));
		  }
	      }	    
	  }
      }
    else
      {
	if ((varying_rho) || (varying_mu))
	  variable_c.Reallocate(var_helm.mesh.GetNbElt());

	if (varying_mu)
	  variable_density.Reallocate(var_helm.mesh.GetNbElt());

	if ((varying_sigma) || (varying_rho))
	  variable_gamma.Reallocate(var_helm.mesh.GetNbElt());

	if (var_helm.include_flow_term)
	  if ((varying_flow) || (varying_rho))
	    variable_M.Reallocate(var_helm.mesh.GetNbElt());
	
	for (int i = 0; i < var_helm.mesh.GetNbElt(); i++)
	  {
	    int ref = var_helm.mesh.Element(i).GetReference();
	    int i1 = var_helm.ElementRho(i);
	    if ((ref < is_elt_ref.GetM()) && (is_elt_ref(ref)))
	      {
		int nb_pts = 0;
		if (varying_rho)
		  nb_pts = var_helm.ref_rho(ref).GetVaryingIndex().GetValue()(i1).GetM();
		else if (varying_mu)
		  nb_pts = var_helm.ref_mu(ref).GetVaryingIndex(0, 0).GetValue()(i1).GetM();
		else if (varying_sigma)
		  nb_pts = var_helm.ref_sigma(ref).GetVaryingIndex().GetValue()(i1).GetM();
		else if (varying_flow)
		  nb_pts = var_helm.ref_v0(ref).GetVaryingIndex(0).GetValue()(i1).GetM();
		
		if ((varying_rho) || (varying_mu))
		  variable_c(i).Reallocate(nb_pts);

		if (varying_mu)
		  variable_density(i).Reallocate(nb_pts);
		
		if ((varying_sigma) || (varying_rho))
		  variable_gamma(i).Reallocate(nb_pts);

		if (var_helm.include_flow_term)
		  if ((varying_flow) || (varying_rho))
		    variable_M(i).Reallocate(nb_pts);
		
		for (int j = 0; j < nb_pts; j++)
		  {
		    Real_wp rho = realpart(var_helm.ref_rho(ref).GetCoefficient(var_helm, i, j));
		    Real_wp mu = realpart(var_helm.ref_mu(ref).GetCoefficient(var_helm, i, j)(0, 0));
		    Real_wp sigma = realpart(var_helm.ref_sigma(ref).GetCoefficient(var_helm, i, j));
		    if ((varying_rho) || (varying_mu))
		      variable_c(i)(j) = sqrt(mu/rho);
		    
		    if (varying_mu)
		      variable_density(i)(j) = Real_wp(1) / mu;
		    
		    if ((varying_sigma) || (varying_rho))
		      variable_gamma(i) = sigma / (Real_wp(2)*rho);

		    if (var_helm.include_flow_term)
		      if ((varying_flow) || (varying_rho))
			{
			  TinyVector<Complex_wp, Dimension::dim_N> flow
			    = var_helm.ref_v0(ref).GetCoefficient(var_helm, i, j);
			  
			  for (int d = 0; d < Dimension::dim_N; d++)
			    variable_M(i)(j)(d) = realpart(flow(d)) / rho;
			}		    
		  }
	      }
	  }
      }
    
    // mass and/or stiffness terms ?
    presence_mass = false;
    presence_stiff = false;
    presence_flow = false;
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO : case PARAM_SIGMA: presence_mass = true; break;
	  case PARAM_MU : presence_stiff = true; break;
          case PARAM_MX : case PARAM_MY : case PARAM_MZ :
            presence_mass = true; presence_stiff = true; presence_flow = true; break;
            
	  case PARAM_C : case PARAM_GAMMA: presence_mass = true; break;
	  case PARAM_DENSITY : presence_mass = true; presence_stiff = true; break;
          case PARAM_UX : case PARAM_UY : case PARAM_UZ :
            presence_mass = true; presence_stiff = true; presence_flow = true; break;
	  }
      }    
    
    if ((var_helm.include_flow_term) && (type_convention == HELIO))
      {
	presence_flow = true;      
	presence_stiff = true;
      }
  }


  //! returns the integer associated with a parameter name
  template<class Dimension>
  int BoundaryInverseProblem<HelmholtzEquation<Dimension> >
  ::GetIdNumberParameter(const string& name) const
  {
    if (name == "rho")
      return PARAM_RHO;
    else if (name == "mu")
      return PARAM_MU;
    else if (name == "mx")
      return PARAM_MX;
    else if (name == "my")
      return PARAM_MY;
    else if (name == "mz")
      return PARAM_MZ;
    else if (name == "sigma")
      return PARAM_SIGMA;
    else if (name == "density")
      return PARAM_DENSITY;
    else if (name == "c")
      return PARAM_C;
    else if (name == "ux")
      return PARAM_UX;
    else if (name == "uy")
      return PARAM_UY;
    else if (name == "uz")
      return PARAM_UZ;
    else if (name == "gamma")
      return PARAM_GAMMA;
    else
      {
	cout << "Unknown parameter" << endl;
	abort();
      }

    return -1;
  }

  
  //! sets the exact parameter c
  template<class Dimension>
  void BoundaryInverseProblem<HelmholtzEquation<Dimension> >::
  SetExactParameter(int k, const string& name, type_ptr_fctC fct)
  {
    BoundaryInverseProblem_Base<Dimension>::SetExactParameter(k, name, fct);
    this->param_c_type(k) = GetIdNumberParameter(name);
    
    if (this->param_c_type(k) >= PARAM_DENSITY)
      {
        if (type_convention == NATURAL)
          {
            cout << "Use natural or helio convention, not both" << endl;
            abort();
          }
        
        type_convention = HELIO;
      }
    else
      {
        if (type_convention == HELIO)
          {
            cout << "Use natural or helio convention, not both" << endl;
            abort();
          }
        
        type_convention = NATURAL;
      }
  }
  

  //! returns the variable indexes in arrays rho_real and rho_cplx
  template<class Dimension>
  void BoundaryInverseProblem<HelmholtzEquation<Dimension> >
  ::GetVariableIndex(int ref,
		     VarHarmonic_Base<Complex_wp, Dimension>& var,
		     Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
		     Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_cplx)
  {
    EllipticProblem<HelmholtzEquation<Dimension> >& var_
      = dynamic_cast<EllipticProblem<HelmholtzEquation<Dimension> >& >(var);
    
    bool rho_added = false, mu_added = false, sigma_added = false;
    bool flow_added = false, beta_added = false;
    for (int num_index = 0; num_index < this->param_c_type.GetM(); num_index++)
      {
	if (this->param_c_type(num_index) == PARAM_RHO)
	  {
	    if (!rho_added)
	      {
		PhysicalVaryingMedia<Dimension, Complex_wp>& rho 
		  = var_.ref_rho(ref).GetVaryingIndex();
		
		var_.ref_rho(ref).SetMediaType(rho.USER);
		rho_cplx.PushBack(&rho);
		rho_added = true;
	      }
	  }
	else if (this->param_c_type(num_index) == PARAM_MU)
	  {
	    if (!mu_added)
	      {
		PhysicalVaryingMedia<Dimension, Complex_wp>& mu 
		  = var_.ref_mu(ref).GetVaryingIndex(0, 0);
		
		var_.ref_mu(ref).SetMediaType(mu.USER);
		rho_cplx.PushBack(&mu);
		mu_added = true;
	      }
	  }
	else if ((this->param_c_type(num_index) == PARAM_MX)
                 || (this->param_c_type(num_index) == PARAM_MY)
                 || (this->param_c_type(num_index) == PARAM_MZ) )
	  {
	    if (!flow_added)
	      for (int k = 0; k < Dimension::dim_N; k++)
		{
		  PhysicalVaryingMedia<Dimension, Complex_wp>& mx
		    = var_.ref_v0(ref).GetVaryingIndex(k);
		  
		  var_.ref_v0(ref).SetMediaType(mx.USER);
		  rho_cplx.PushBack(&mx);
		}
	    
	    flow_added = true;
	  }
	else if (this->param_c_type(num_index) == PARAM_SIGMA)
	  {
	    if (!sigma_added)
	      {
		PhysicalVaryingMedia<Dimension, Complex_wp>& sigma
		  = var_.ref_sigma(ref).GetVaryingIndex();
		
		var_.ref_sigma(ref).SetMediaType(sigma.USER);
		rho_cplx.PushBack(&sigma);
		sigma_added = true;
	      }
	  }
	else if (this->param_c_type(num_index) == PARAM_DENSITY)
	  {
	    if (!rho_added)
	      {
		PhysicalVaryingMedia<Dimension, Complex_wp>& rho 
		  = var_.ref_rho(ref).GetVaryingIndex();
		
		var_.ref_rho(ref).SetMediaType(rho.USER);
		rho_cplx.PushBack(&rho);
		rho_added = true;
	      }
	    
	    if (!mu_added)
	      {
		PhysicalVaryingMedia<Dimension, Complex_wp>& mu 
		  = var_.ref_mu(ref).GetVaryingIndex(0, 0);
		
		var_.ref_mu(ref).SetMediaType(mu.USER);
		rho_cplx.PushBack(&mu);
		mu_added = true;
	      }
	    
	    if (!sigma_added)
	      {
		PhysicalVaryingMedia<Dimension, Complex_wp>& sigma 
		  = var_.ref_sigma(ref).GetVaryingIndex();
		
		var_.ref_sigma(ref).SetMediaType(sigma.USER);
		rho_cplx.PushBack(&sigma);
		sigma_added = true;
	      }

	    if (var_helm.include_flow_term)
              {
                if ((var_helm.include_double_gradM_flow) && (!beta_added))
                  {
                    PhysicalVaryingMedia<Dimension, Complex_wp>& beta 
                      = var_.ref_beta(ref).GetVaryingIndex();
                    
                    var_.ref_beta(ref).SetMediaType(beta.USER);
                    rho_cplx.PushBack(&beta);
		    beta_added = true;
                  }
                
		if (!flow_added)
		  for (int k = 0; k < Dimension::dim_N; k++)
		    {
		      PhysicalVaryingMedia<Dimension, Complex_wp>& mx
			= var_.ref_v0(ref).GetVaryingIndex(k);
		      
		      var_.ref_v0(ref).SetMediaType(mx.USER);
		      rho_cplx.PushBack(&mx);
		    }
		
		flow_added = true;
              }
	  }
	else if (this->param_c_type(num_index) == PARAM_C)
	  {
	    if (!rho_added)
	      {
		PhysicalVaryingMedia<Dimension, Complex_wp>& rho 
		  = var_.ref_rho(ref).GetVaryingIndex();
		
		var_.ref_rho(ref).SetMediaType(rho.USER);
		rho_cplx.PushBack(&rho);
		rho_added = true;
	      }
	    
	    if (!sigma_added)
	      {
		PhysicalVaryingMedia<Dimension, Complex_wp>& sigma 
		  = var_.ref_sigma(ref).GetVaryingIndex();
		
		var_.ref_sigma(ref).SetMediaType(sigma.USER);
		rho_cplx.PushBack(&sigma);
		sigma_added = true;
	      }

	    if (var_helm.include_flow_term)
              {
		if ((var_helm.include_double_gradM_flow) && (!beta_added))
                  {
                    PhysicalVaryingMedia<Dimension, Complex_wp>& beta 
                      = var_.ref_beta(ref).GetVaryingIndex();
                    
                    var_.ref_beta(ref).SetMediaType(beta.USER);
                    rho_cplx.PushBack(&beta);
		    beta_added = true;
                  }
		
                if (!flow_added)
		  for (int k = 0; k < Dimension::dim_N; k++)
		    {
		      PhysicalVaryingMedia<Dimension, Complex_wp>& mx
			= var_.ref_v0(ref).GetVaryingIndex(k);
		      
		      var_.ref_v0(ref).SetMediaType(mx.USER);
		      rho_cplx.PushBack(&mx);
		    }
		
		flow_added = true;
              }
	  }
	else if (this->param_c_type(num_index) == PARAM_GAMMA)
	  {
	    if (!sigma_added)
	      {
		PhysicalVaryingMedia<Dimension, Complex_wp>& sigma
		  = var_.ref_sigma(ref).GetVaryingIndex();
		
		var_.ref_sigma(ref).SetMediaType(sigma.USER);
		rho_cplx.PushBack(&sigma);
		sigma_added = true;
	      }
	  }
	else if ((this->param_c_type(num_index) == PARAM_UX)
                 || (this->param_c_type(num_index) == PARAM_UY)
                 || (this->param_c_type(num_index) == PARAM_UZ) )
	  {
	    if (!flow_added)
	      for (int k = 0; k < Dimension::dim_N; k++)
		{
		  PhysicalVaryingMedia<Dimension, Complex_wp>& mx
		    = var_.ref_v0(ref).GetVaryingIndex(k);
		  
		  var_.ref_v0(ref).SetMediaType(mx.USER);
		  rho_cplx.PushBack(&mx);
		}
	    
	    flow_added = true;
	    if ((var_helm.include_double_gradM_flow) && (!beta_added))
              {
                PhysicalVaryingMedia<Dimension, Complex_wp>& beta
                  = var_.ref_beta(ref).GetVaryingIndex();
	    
                var_.ref_beta(ref).SetMediaType(beta.USER);
                rho_cplx.PushBack(&beta);
		beta_added = true;
              }
	  }
	else
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
      }
  }
  
  
  //! Fills variables indexes on element i1, point j from values contained in eval_c
  template<class Dimension>
  void BoundaryInverseProblem<HelmholtzEquation<Dimension> >
  ::SetVariableIndexPoint(int ref, int i, int i1, int j, int offset, const Vector<VectReal_wp>& eval_c,
			  Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
			  Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_cplx)
  {
    R_N flow = ref_M(ref); Real_wp rho = ref_density(ref), c = ref_c(ref);
    Real_wp gamma = ref_gamma(ref);
    if ((variable_M.GetM() > 0) && (variable_M(i).GetM() > 0))
      flow = variable_M(i)(j);

    if ((variable_density.GetM() > 0) && (variable_density(i).GetM() > 0))
      rho = variable_density(i)(j);

    if ((variable_c.GetM() > 0) && (variable_c(i).GetM() > 0))
      c = variable_c(i)(j);
    
    if ((variable_gamma.GetM() > 0) && (variable_gamma(i).GetM() > 0))
      gamma = variable_gamma(i)(j);
    
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	if (this->param_c_type(p) == PARAM_DENSITY)
	  rho = eval_c(p)(offset+j);
	else if (this->param_c_type(p) == PARAM_C)
	  c = eval_c(p)(offset+j);
	else if (this->param_c_type(p) == PARAM_GAMMA)
	  gamma = eval_c(p)(offset+j);
        else if ((this->param_c_type(p) == PARAM_MX)
                 || (this->param_c_type(p) == PARAM_MY)
                 || (this->param_c_type(p) == PARAM_MZ))
          {
            int comp = this->param_c_type(p) - PARAM_MX;
            flow(comp) = eval_c(p)(offset+j);
          }
        else if ((this->param_c_type(p) == PARAM_UX)
                 || (this->param_c_type(p) == PARAM_UY)
                 || (this->param_c_type(p) == PARAM_UZ))
          {
            int comp = this->param_c_type(p) - PARAM_UX;
            flow(comp) = eval_c(p)(offset+j);
          }
      }
    
    int num = 0;
    int p_flow = -1, p_beta = -1;
    bool rho_added = false, mu_added = false, sigma_added = false;
    bool flow_added = false, beta_added = false;
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	if (this->param_c_type(p) == PARAM_RHO)
	  {
	    if (!rho_added)
	      {
		rho_cplx(num)->SetCoefficient(i1, j, eval_c(p)(offset+j));
		num++; rho_added = true;
	      }
	  }
	else if (this->param_c_type(p) == PARAM_MU)
	  {
	    if (!mu_added)
	      {
		rho_cplx(num)->SetCoefficient(i1, j, eval_c(p)(offset+j));
		num++; mu_added = true;
	      }
	  }
        else if ((this->param_c_type(p) == PARAM_MX) 
                 || (this->param_c_type(p) == PARAM_MY)
                 || (this->param_c_type(p) == PARAM_MZ))
          {
            if (!flow_added)
	      {
		p_flow = num;
		num += Dimension::dim_N;
		flow_added = true;
	      }
          }
	else if (this->param_c_type(p) == PARAM_SIGMA)
	  {
	    if (!sigma_added)
	      {
		rho_cplx(num)->SetCoefficient(i1, j, eval_c(p)(offset+j));
		num++; sigma_added = true;
	      }
	  }
        else if (this->param_c_type(p) == PARAM_DENSITY)
	  {
	    if (!rho_added)
	      {
		rho_cplx(num)->SetCoefficient(i1, j, Complex_wp(1, 0)/(rho*c*c));
		num++; rho_added = true;
	      }
	    
	    if (!mu_added)
	      {
		rho_cplx(num)->SetCoefficient(i1, j, Complex_wp(1, 0)/rho);
		num++; mu_added = true;
	      }

	    if (!sigma_added)
	      {
		rho_cplx(num)->SetCoefficient(i1, j, Complex_wp(2, 0)*gamma/(rho*c*c));
		num++; sigma_added = true;
	      }
            
            if (var_helm.include_flow_term)
              {
                if ((var_helm.include_double_gradM_flow) && (!beta_added))
                  p_beta = num++;
                
                if (!flow_added)
		  {
		    p_flow = num;
		    num += Dimension::dim_N;
		  }
		
		beta_added = true; flow_added = true;
              }
	  }
	else if (this->param_c_type(p) == PARAM_C)
	  {
	    if (!rho_added)
	      {
		rho_cplx(num)->SetCoefficient(i1, j, Complex_wp(1, 0)/(rho*c*c));
		num++; rho_added = true;
	      }
	    
	    if (!sigma_added)
	      {
		rho_cplx(num)->SetCoefficient(i1, j, Complex_wp(2, 0)*gamma/(rho*c*c));
		num++; sigma_added = true;
	      }
	    
            if (var_helm.include_flow_term)
              {
                if ((var_helm.include_double_gradM_flow) && (!beta_added))
                  p_beta = num++;
                
                if (!flow_added)
		  {
		    p_flow = num;
		    num += Dimension::dim_N;
		  }
		
		beta_added = true; flow_added = true;
              }
	  }
	else if (this->param_c_type(p) == PARAM_GAMMA)
	  {
	    if (!sigma_added)
	      {
		rho_cplx(num)->SetCoefficient(i1, j, Complex_wp(2, 0)*gamma/(rho*c*c));
		num++; sigma_added = true;
	      }
	  }
        else if ((this->param_c_type(p) == PARAM_UX) 
                 || (this->param_c_type(p) == PARAM_UY)
                 || (this->param_c_type(p) == PARAM_UZ))
          {
            if (!flow_added)
	      {
		p_flow = num;
		num += Dimension::dim_N;
		flow_added = true;
	      }

            if ((var_helm.include_double_gradM_flow) && (!beta_added))
	      {
		p_beta = num++;
		beta_added = true;
	      }
          }
	else
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
      }	

    if (type_convention == NATURAL)
      {
        if (p_flow >= 0)
          for (int m = 0; m < Dimension::dim_N; m++)
            rho_cplx(p_flow+m)->SetCoefficient(i1, j, flow(m));
      }
    else
      {
        if (p_flow >= 0)
          for (int m = 0; m < Dimension::dim_N; m++)
            rho_cplx(p_flow+m)->SetCoefficient(i1, j, flow(m)/(rho*c*c));
        
        if (p_beta >= 0)
          rho_cplx(p_beta)->SetCoefficient(i1, j, rho*c*c);
      }
  }
  
  
  //! Retrieves values of M, beta, rho and c of element i
  template<class Dimension>
  void BoundaryInverseProblem<HelmholtzEquation<Dimension> >
  ::GetPhysicalIndexesElem(int i, int offset, int nb_pts, TinyVector<VectReal_wp, Dimension::dim_N>& flow,
                           VectReal_wp& beta, VectReal_wp& rho, VectReal_wp& c, VectReal_wp& gamma)
  {
    int ref = var_helm.mesh.Element(i).GetReference();
    if (var_helm.include_flow_term)
      {
        for (int m = 0; m < Dimension::dim_N; m++)
          {
            flow(m).Reallocate(nb_pts);
	    if ((variable_M.GetM() > 0) && (variable_M(i).GetM() > 0))
	      {
		for (int j = 0; j < nb_pts; j++)
		  flow(m)(j) = variable_M(i)(j)(m);
	      }
	    else
	      flow(m).Fill(ref_M(ref)(m));
	  }
	
	if ((variable_beta.GetM() > 0) && (variable_beta(i).GetM() > 0))
	  beta = variable_beta(i);
	else
	  {
	    beta.Reallocate(nb_pts);
	    beta.Fill(ref_beta(ref));
	  }
      }
    
    if (type_convention == HELIO)
      {
	if ((variable_gamma.GetM() > 0) && (variable_gamma(i).GetM() > 0))
	  gamma = variable_gamma(i);
	else
	  {
	    gamma.Reallocate(nb_pts);
	    gamma.Fill(ref_gamma(ref));
	  }
	
	if ((variable_density.GetM() > 0) && (variable_density(i).GetM() > 0))
	  rho = variable_density(i);
	else
	  {
	    rho.Reallocate(nb_pts);
	    rho.Fill(ref_density(ref));
	  }
	
	if ((variable_c.GetM() > 0) && (variable_c(i).GetM() > 0))
	  c = variable_c(i);
	else
	  {
	    c.Reallocate(nb_pts);
	    c.Fill(ref_c(ref));
	  }	    
      }    
    else
      {
	if ((variable_rho.GetM() > 0) && (variable_rho(i).GetM() > 0))
	  rho = variable_rho(i);
	else
	  {
	    rho.Reallocate(nb_pts);
	    rho.Fill(ref_rho(ref));
	  }

	if ((variable_mu.GetM() > 0) && (variable_mu(i).GetM() > 0))
	  c = variable_mu(i);
	else
	  {
	    c.Reallocate(nb_pts);
	    c.Fill(ref_mu(ref));
	  }

	if ((variable_sigma.GetM() > 0) && (variable_sigma(i).GetM() > 0))
	  gamma = variable_sigma(i);
	else
	  {
	    gamma.Reallocate(nb_pts);
	    gamma.Fill(ref_sigma(ref));
	  }
      }    

    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
          {
          case PARAM_MX: 
          case PARAM_MY:
          case PARAM_MZ:
            {
              int comp = this->param_c_type(p) - PARAM_MX;
              for (int k = 0; k < nb_pts; k++)
                flow(comp)(k) = this->eval_c_quad(p)(offset+k);
            }
            break;
          case PARAM_RHO:
          case PARAM_DENSITY:
            {
              for (int k = 0; k < nb_pts; k++)
                rho(k) = this->eval_c_quad(p)(offset+k);
            }
            break;
          case PARAM_C:
          case PARAM_MU:
            {
              for (int k = 0; k < nb_pts; k++)
                c(k) = this->eval_c_quad(p)(offset+k);
            }
            break;
          case PARAM_UX: 
          case PARAM_UY:
          case PARAM_UZ:
            {
              int comp = this->param_c_type(p) - PARAM_UX;
              for (int k = 0; k < nb_pts; k++)
                flow(comp)(k) = this->eval_c_quad(p)(offset+k);
            }
            break;
          case PARAM_GAMMA:
          case PARAM_SIGMA:
            {
              for (int k = 0; k < nb_pts; k++)
                gamma(k) = this->eval_c_quad(p)(offset+k);
            }
            break;
          }
      }
  }

  
  //! local product by derivative with respect to parameters 
  template<class Dimension>
  void BoundaryInverseProblem<HelmholtzEquation<Dimension> >
  ::MltParamDerivLocal(int nk, const Real_wp& omega2, VectComplex_wp& u_glob,
		       int i, const Real_wp& alpha, const Vector<VectReal_wp>& c_quad,
		       int offset, VectComplex_wp& y_col)
  {
    // we count the number of quadrature points
    const ElementReference<Dimension, 1>& Fb = var_helm.GetReferenceElementH1(i);
    int nb_pts = Fb.GetNbPointsQuadratureInside();
    for (int num_loc = 0; num_loc < var_helm.mesh.Element(i).GetNbBoundary(); num_loc++)
      nb_pts += Fb.GetNbQuadBoundary(num_loc);

    // values of physical indexes are filled
    VectReal_wp beta, rho, c, gamma;
    TinyVector<VectReal_wp, Dimension::dim_N> flow;

    GetPhysicalIndexesElem(i, offset, nb_pts, flow, beta, rho, c, gamma);

    // local values of u are retrieved
    TinyVector<Complex_wp, Dimension::dim_N> vec_u, vec_v;
    VectComplex_wp gradUh_quad, Uh_quad, feval, grad_feval, Uh_val;
    TinyVector<VectComplex_wp, 1> Uh_dof;
    var_helm.GetLocalUnknownVector(u_glob, i, Uh_dof);
    
    Uh_val = Uh_dof(0);
    // values of u (and/or grad u) on quadrature points are computed
    // Uh_quad contains u(xi_k) and gradUh_quad contains \nabla u(\xi_k)
    int Nquad = Fb.GetNbPointsQuadratureInside();
    Uh_quad.Reallocate(Nquad);
    if (presence_stiff)
      gradUh_quad.Reallocate(Nquad*Dimension::dim_N);
    
    if (Fb.UseQuadraturePointsForRh())
      {
	Fb.ApplyChTranspose(Uh_dof(0), Uh_quad);
	if (presence_stiff)
	  Fb.ApplyRhQuadratureTranspose(Uh_quad, gradUh_quad);
      }
    else
      {
	if (presence_mass)
	  Fb.ApplyChTranspose(Uh_dof(0), Uh_quad);
	
	if (presence_stiff)
	  Fb.ApplyRhTranspose(Uh_dof(0), gradUh_quad);
      }
    
    // we retrieve weight of integration and matrices DF_i^{-1}
    VectReal_wp poids(Nquad);
    Vector<typename Dimension::MatrixN_N> dfjm1;
    if (presence_stiff)
      dfjm1.Reallocate(Nquad);
    
    bool affine = var_helm.mesh.IsElementAffine(i);
    if (affine)
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(0)*Fb.WeightsND(k);
	  if (presence_stiff)
	    {
	      dfjm1(k) = var_helm.Glob_DFjm1(i)(0);
	      Mlt(1.0/var_helm.Glob_jacobian(i)(0), dfjm1(k));
	    }
	} 
    else
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(k);
	  if (presence_stiff)
	    {
	      dfjm1(k) = var_helm.Glob_DFjm1(i)(k);
	      Mlt(Fb.WeightsND(k)/var_helm.Glob_jacobian(i)(k), dfjm1(k));
	    }
	}
    
    feval.Reallocate(Nquad);
    feval.Zero();
    
    if (presence_stiff)
      {
	grad_feval.Reallocate(Nquad*Dimension::dim_N);
	grad_feval.Zero();

        // application of transformation DF_i^{-T}
        for (int k = 0; k < Nquad; k++)
          {
            CopyVector(gradUh_quad, k, vec_u);
            MltTrans(dfjm1(k), vec_u, vec_v);
            CopyVector(vec_v, k, gradUh_quad);
          }
      }

    VectComplex_wp MdotGradU;
    if (presence_flow)
      {
        MdotGradU.Reallocate(Nquad);
        MdotGradU.Zero();
        for (int k = 0; k < Nquad; k++)
          {
            CopyVector(gradUh_quad, k, vec_u);
            for (int m = 0; m < Dimension::dim_N; m++)
              MdotGradU(k) += flow(m)(k)*vec_u(m);
          }
      }
    
    // we update feval and grad_feval by using the parameters
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO:
	    {
	      Real_wp coef = -omega2;
	      for (int k = 0; k < Nquad; k++)
		feval(k) += coef*c_quad(p)(offset+k)*Uh_quad(k);
	    }
	    break;
	  case PARAM_MU:
	    {
	      for (int k = 0; k < Nquad; k++)
		{
		  CopyVector(gradUh_quad, k, vec_u);
		  vec_u *= c_quad(p)(offset+k);
                  AddVector(vec_u, Dimension::dim_N*k, grad_feval);
		}
	    }
	    break;
          case PARAM_SIGMA:
            {
	      Complex_wp coef = -Iwp*this->omega(nk);
	      for (int k = 0; k < Nquad; k++)
		feval(k) += coef*c_quad(p)(offset+k)*Uh_quad(k);              
            }
            break;
          case PARAM_MX:
          case PARAM_MY:
          case PARAM_MZ:
            {
              int comp = this->param_c_type(p) - PARAM_MX;
              Complex_wp coef = -Iwp*this->omega(nk);
              if (var_helm.include_flow_term)
                for (int k = 0; k < Nquad; k++)
                  {
                    int off_k = Dimension::dim_N*k;
                    feval(k) += coef*c_quad(p)(offset+k)*gradUh_quad(off_k+comp);
                    grad_feval(off_k+comp) -= coef*c_quad(p)(offset+k)*Uh_quad(k);
                    
                    if (var_helm.include_double_gradM_flow)
                      {
                        for (int m = 0; m < Dimension::dim_N; m++)
                          {
                            grad_feval(off_k+comp) -= beta(k)*flow(m)(k)*gradUh_quad(off_k+m)*c_quad(p)(offset+k);
                            grad_feval(off_k+m) -= beta(k)*flow(m)(k)*gradUh_quad(off_k+comp)*c_quad(p)(offset+k);
                          } 
                      }
                  }
            }
            break;            
	  case PARAM_DENSITY:
	    {
	      Real_wp coef = omega2;
	      for (int k = 0; k < Nquad; k++)
		{
                  Real_wp coefd = c_quad(p)(offset+k) / square(rho(k)*c(k));
		  feval(k) += coef*coefd*Uh_quad(k);
                  
		  CopyVector(gradUh_quad, k, vec_u);
		  vec_u *= -c_quad(p)(offset+k)/(rho(k)*rho(k));
		  AddVector(vec_u, Dimension::dim_N*k, grad_feval);
                  
                  feval(k) += Real_wp(2)*Iwp*this->omega(nk)*gamma(k)*coefd*Uh_quad(k);
                  
                  if (var_helm.include_flow_term)
                    {
                      int off_k = Dimension::dim_N*k;
                      Complex_wp coefb = Iwp*this->omega(nk)*coefd;
                      for (int m = 0; m < Dimension::dim_N; m++)
                        {
                          feval(k) += coefb*gradUh_quad(off_k+m)*flow(m)(k);
                          grad_feval(off_k+m) -= coefb*Uh_quad(k)*flow(m)(k);
                          
                          if (var_helm.include_double_gradM_flow)
                            grad_feval(off_k+m) += coefd*MdotGradU(k)*flow(m)(k);
                        } 
                    }
		}
	    }
	    break;
	  case PARAM_C:
	    {
	      Real_wp coef = omega2;
	      for (int k = 0; k < Nquad; k++)
                {
                  Complex_wp coefc = 2*c_quad(p)(offset+k) / (rho(k)*c(k)*c(k)*c(k));
                  feval(k) += coef*coefc*Uh_quad(k);

                  feval(k) += Real_wp(2)*Iwp*this->omega(nk)*gamma(k)*coefc*Uh_quad(k);
                  if (var_helm.include_flow_term)
                    {
                      int off_k = Dimension::dim_N*k;
                      Complex_wp coefb = Iwp*this->omega(nk)*coefc;
                      for (int m = 0; m < Dimension::dim_N; m++)
                        {
                          feval(k) += coefb*gradUh_quad(off_k+m)*flow(m)(k);
                          grad_feval(off_k+m) -= coefb*Uh_quad(k)*flow(m)(k);
                          
                          if (var_helm.include_double_gradM_flow)
                            grad_feval(off_k+m) += coefc*MdotGradU(k)*flow(m)(k);
                        } 
                    }
                }
	    }
	    break;
          case PARAM_GAMMA:
            {
	      Complex_wp coef = -Real_wp(2)*Iwp*this->omega(nk);
	      for (int k = 0; k < Nquad; k++)
		feval(k) += coef*c_quad(p)(offset+k)*Uh_quad(k)/(rho(k)*c(k)*c(k));     
            }
            break;
          case PARAM_UX:
          case PARAM_UY:
          case PARAM_UZ:
            {
              int comp = this->param_c_type(p) - PARAM_UX;
              Complex_wp coef = -Iwp*this->omega(nk);
              if (var_helm.include_flow_term)
                for (int k = 0; k < Nquad; k++)
                  {
                    int off_k = Dimension::dim_N*k;
                    Complex_wp coefb = coef / (rho(k)*c(k)*c(k));
                    feval(k) += coefb*c_quad(p)(offset+k)*gradUh_quad(off_k+comp);
                    grad_feval(off_k+comp) -= coefb*c_quad(p)(offset+k)*Uh_quad(k);
                    
                    if (var_helm.include_double_gradM_flow)
                      {
                        coefb = c_quad(p)(offset+k) / (rho(k)*c(k)*c(k));
                        for (int m = 0; m < Dimension::dim_N; m++)
                          {
                            grad_feval(off_k+comp) -= flow(m)(k)*gradUh_quad(off_k+m)*coefb;
                            grad_feval(off_k+m) -= flow(m)(k)*gradUh_quad(off_k+comp)*coefb;
                          } 
                      }
                  }
            }
            break;
	  }
      }
    
    // applying transformation DF_i and coefficient -alpha
    for (int k = 0; k < Nquad; k++)
      {
        feval(k) *= -alpha*poids(k);
        if (presence_stiff)
          {
            CopyVector(grad_feval, k, vec_u);
            Mlt(dfjm1(k), vec_u, vec_v);
            vec_v *= -alpha*poids(k);
            CopyVector(vec_v, k, grad_feval);
          }
      }
	
    // y_col is updated with feval and grad_feval
    if (Fb.UseQuadraturePointsForRh())
      {
	if (presence_stiff)
	  {
	    Fb.ApplyRhQuadrature(grad_feval, Uh_dof(0));
	    feval = feval + Uh_dof(0);
	  }
	
	Fb.ApplyCh(feval, Uh_dof(0));
      }
    else
      {	
	if (presence_mass)
	  Fb.ApplyCh(feval, Uh_dof(0));
	else
	  Uh_dof(0).Zero();
	
	if (presence_stiff)
	  {
	    Uh_quad.Reallocate(Uh_dof(0).GetM());
	    Fb.ApplyRh(grad_feval, Uh_quad);
	    Uh_dof(0) = Uh_dof(0) + Uh_quad;
	  }	
      }
    
    // contribution of absorbing boundary condition
    const MeshNumbering<Dimension>& mesh_num_helm = var_helm.GetMeshNumbering(0);
    Complex_wp coef = Real_wp(0.5)*alpha*Iwp*sqrt(omega2);
    for (int num_loc = 0; num_loc < var_helm.mesh.Element(i).GetNbBoundary(); num_loc++)
      {
	int nf = var_helm.mesh.Element(i).numBoundary(num_loc);
	int ref_line = var_helm.mesh.Boundary(nf).GetReference();
	if (var_helm.mesh.GetBoundaryCondition(ref_line) == BoundaryConditionEnum::LINE_ABSORBING)
	  {
	    nb_pts = mesh_num_helm.GetNbPointsQuadratureBoundary(nf);
	    poids.Reallocate(nb_pts);
	    for (int k = 0; k < nb_pts; k++)
	      poids(k) = var_helm.Glob_dsj(nf)(k)*Fb.WeightsQuadratureBoundary(k, num_loc);
	    
	    Uh_quad.Reallocate(nb_pts);
	    feval.Reallocate(nb_pts); feval.Zero();
	    Fb.ApplyShTranspose(num_loc, Uh_val, Uh_quad);
	    int of = Fb.GetNbPointsQuadratureInside();
	    for (int k = 0; k < num_loc; k++)
	      of += mesh_num_helm.GetNbPointsQuadratureBoundary(var_helm.mesh.Element(i).numBoundary(k));
	    
	    for (int p = 0; p < this->param_c_type.GetM(); p++)
	      {
		switch (this->param_c_type(p))
		  {
		  case PARAM_RHO:
		    {
		      for (int k = 0; k < nb_pts; k++)
                        {
                          Complex_wp zeta = Iwp*gamma(of+k)/this->omega(nk);
                          feval(k) += coef*Uh_quad(k)*poids(k)*
                            sqrt(c(of+k)/(rho(of+k) + zeta))*c_quad(p)(offset+of+k);
			}
		    }
		    break;
		  case PARAM_MU:
		    {
		      for (int k = 0; k < nb_pts; k++)
			{
			  Complex_wp zeta = Iwp*gamma(of+k)/this->omega(nk);
			  feval(k) += coef*Uh_quad(k)*poids(k)*
                            sqrt((rho(of+k)+zeta)/c(of+k))*c_quad(p)(offset+of+k);
			}
		    }
		    break;
		  case PARAM_SIGMA:
		    {
		      for (int k = 0; k < nb_pts; k++)
                        {
                          Complex_wp zeta = Iwp*gamma(of+k)/this->omega(nk);
                          feval(k) -= alpha/2*Uh_quad(k)*poids(k)*
                            sqrt(c(of+k)/(rho(of+k) + zeta))*c_quad(p)(offset+of+k);
			}
		    }
		    break;
		  case PARAM_DENSITY:
		    {
		      for (int k = 0; k < nb_pts; k++)
			{
                          Complex_wp vloc = Real_wp(2)*coef/(square(rho(of+k))*c(of+k))
                            * sqrt(Real_wp(1) + Real_wp(2)*Iwp*gamma(of+k)/this->omega(nk));
			  feval(k) -= vloc*Uh_quad(k)*poids(k)*c_quad(p)(offset+of+k);
			}
		    }
		    break;
		  case PARAM_C:
		    {
		      for (int k = 0; k < nb_pts; k++)
			{
                          Complex_wp vloc = Real_wp(2)*coef/(square(c(of+k))*rho(of+k))
                            * sqrt(Real_wp(1) + Real_wp(2)*Iwp*gamma(of+k)/this->omega(nk));
			  feval(k) -= vloc*Uh_quad(k)*poids(k)*c_quad(p)(offset+of+k);
			}
		    }
		    break;
		  case PARAM_GAMMA:
		    {
		      for (int k = 0; k < nb_pts; k++)
			{
                          Complex_wp vloc = Real_wp(1)/(rho(of+k)*c(of+k))
                            / sqrt(Real_wp(1) + Real_wp(2)*Iwp*gamma(of+k)/this->omega(nk));
			  feval(k) -= vloc*Uh_quad(k)*poids(k)*c_quad(p)(offset+of+k);
			}
		    }
		    break;
                    
		  }
	      }
	    
	    Fb.ApplySh(Complex_wp(1, 0), num_loc, feval, Uh_dof(0), Fb.GetOrder());
	  }
      }
    
    var_helm.AddLocalUnknownVector(Complex_wp(1, 0), Uh_dof, i, y_col);   
  }


  //! local product by transpose of derivative with respect to parameters   
  template<class Dimension>
  void BoundaryInverseProblem<HelmholtzEquation<Dimension> >
  ::MltParamDerivTransLoc(int nk, const Real_wp& omega2, VectComplex_wp& u_glob,
			  int i, const Real_wp& alpha, int offset, 
			  const VectComplex_wp& y_col, Vector<VectReal_wp>& ProdC)
  {
    // we count the number of quadrature points
    const ElementReference<Dimension, 1>& Fb = var_helm.GetReferenceElementH1(i);
    int nb_pts = Fb.GetNbPointsQuadratureInside();
    for (int num_loc = 0; num_loc < var_helm.mesh.Element(i).GetNbBoundary(); num_loc++)
      nb_pts += Fb.GetNbQuadBoundary(num_loc);
    
    // we retrieve values of physical indexes
    VectReal_wp beta, rho, c, gamma;
    TinyVector<VectReal_wp, Dimension::dim_N> flow;
    GetPhysicalIndexesElem(i, offset, nb_pts, flow, beta, rho, c, gamma);

    // values of Uh and Vh are retrieved on element i
    VectComplex_wp Uh_quad, Vh_quad, gradUh_quad, gradVh_quad;
    TinyVector<VectComplex_wp, 1> Uh_dof, Vh_dof;
    var_helm.GetLocalUnknownVector(u_glob, i, Uh_dof);
    
    var_helm.GetLocalUnknownVector(y_col, i, Vh_dof);
    
    // values of u and v are computed on quadrature points (and gradient if needed)
    int Nquad = Fb.GetNbPointsQuadratureInside();
    Uh_quad.Reallocate(Nquad);
    Vh_quad.Reallocate(Nquad);
    
    if (presence_stiff)
      {
	gradUh_quad.Reallocate(Dimension::dim_N*Nquad);
	gradVh_quad.Reallocate(Dimension::dim_N*Nquad);
      }
    
    if (Fb.UseQuadraturePointsForRh())
      {
	Fb.ApplyChTranspose(Uh_dof(0), Uh_quad);
	Fb.ApplyChTranspose(Vh_dof(0), Vh_quad);	
	if (presence_stiff)
	  {
	    Fb.ApplyRhQuadratureTranspose(Uh_quad, gradUh_quad);
	    Fb.ApplyRhQuadratureTranspose(Vh_quad, gradVh_quad);
	  }
      }
    else
      {
	if (presence_mass)
	  {
	    Fb.ApplyChTranspose(Uh_dof(0), Uh_quad);
	    Fb.ApplyChTranspose(Vh_dof(0), Vh_quad);
	  }
	
	if (presence_stiff)
	  {
	    Fb.ApplyRhTranspose(Uh_dof(0), gradUh_quad);
	    Fb.ApplyRhTranspose(Vh_dof(0), gradVh_quad);
	  }
      }
    
    // we retrieve weight of integration
    VectReal_wp poids(Nquad);
    Vector<typename Dimension::MatrixN_N> dfjm1;
    if (presence_stiff)
      dfjm1.Reallocate(Nquad);

    bool affine = var_helm.mesh.IsElementAffine(i);
    if (affine)
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(0)*Fb.WeightsND(k);
	  if (presence_stiff)
	    {
	      dfjm1(k) = var_helm.Glob_DFjm1(i)(0);
	      Mlt(1.0/var_helm.Glob_jacobian(i)(0), dfjm1(k));
	    }
	} 
    else
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(k);
	  if (presence_stiff)
	    {
	      dfjm1(k) = var_helm.Glob_DFjm1(i)(k);
	      Mlt(Fb.WeightsND(k)/var_helm.Glob_jacobian(i)(k), dfjm1(k));
	    }
	}

    TinyVector<Complex_wp, Dimension::dim_N> vec_u, vec_v;
    if (presence_stiff)
      {
        // application of transformation DF_i^{-T}
        for (int k = 0; k < Nquad; k++)
          {
            CopyVector(gradUh_quad, k, vec_u);
            MltTrans(dfjm1(k), vec_u, vec_v);
            CopyVector(vec_v, k, gradUh_quad);

            CopyVector(gradVh_quad, k, vec_u);
            MltTrans(dfjm1(k), vec_u, vec_v);
            CopyVector(vec_v, k, gradVh_quad);
          }
      }

    VectComplex_wp MdotGradU, MdotGradV;
    if (presence_flow)
      {
        MdotGradU.Reallocate(Nquad); MdotGradV.Reallocate(Nquad);
        MdotGradU.Zero(); MdotGradV.Zero();
        for (int k = 0; k < Nquad; k++)
          {
            CopyVector(gradUh_quad, k, vec_u);
            CopyVector(gradVh_quad, k, vec_v);

            for (int m = 0; m < Dimension::dim_N; m++)
              {
                MdotGradU(k) += flow(m)(k)*vec_u(m);
                MdotGradV(k) += flow(m)(k)*vec_v(m);
              }
          }
      }
    
    // we update ProdC
    TinyVector<Real_wp, Dimension::dim_N> vec_ur, vec_ui, vec_vr, vec_vi;
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO:
	    {
	      Real_wp coef = alpha*omega2;
	      for (int k = 0; k < Nquad; k++)
		ProdC(p)(offset+k) += coef*poids(k)*(real(Uh_quad(k))*real(Vh_quad(k))
						     + imag(Uh_quad(k))*imag(Vh_quad(k)));
	    }
	    break;
	  case PARAM_MU:
	    {
	      Real_wp coef = -alpha;
	      for (int k = 0; k < Nquad; k++)
		{
		  CopyVector(gradUh_quad, k, vec_ur, vec_ui);
		  CopyVector(gradVh_quad, k, vec_vr, vec_vi);
		  ProdC(p)(offset+k) += coef*poids(k)*(DotProd(vec_ur, vec_vr) + DotProd(vec_ui, vec_vi));
		}
	    }
	    break;
          case PARAM_SIGMA:
            {
              Real_wp coef = alpha*this->omega(nk);
              for (int k = 0; k < Nquad; k++)
                ProdC(p)(offset+k) += coef*poids(k)*(real(Uh_quad(k))*imag(Vh_quad(k))
                                                     - imag(Uh_quad(k))*real(Vh_quad(k)));
            }
            break;
          case PARAM_MX:
          case PARAM_MY:
          case PARAM_MZ:
            {
              int comp = this->param_c_type(p) - PARAM_MX;
              Real_wp coef = alpha*this->omega(nk);
              for (int k = 0; k < Nquad; k++)
                {
                  int off_k = Dimension::dim_N*k;
                  ProdC(p)(offset+k) += coef*poids(k)*
                    (real(gradUh_quad(off_k+comp))*imag(Vh_quad(k))
                     - imag(gradUh_quad(off_k+comp))*real(Vh_quad(k))
                     + real(gradVh_quad(off_k+comp))*imag(Uh_quad(k))
                     - imag(gradVh_quad(off_k+comp))*real(Uh_quad(k)));

                  if (var_helm.include_double_gradM_flow)
                    ProdC(p)(offset+k) += alpha*poids(k)*beta(k)*
                      (real(MdotGradU(k))*real(gradVh_quad(off_k+comp))
                       + imag(MdotGradU(k))*imag(gradVh_quad(off_k+comp))
                       + real(MdotGradV(k))*real(gradUh_quad(off_k+comp))
                       + imag(MdotGradV(k))*imag(gradUh_quad(off_k+comp)));
                }
            }
            break;
	  case PARAM_DENSITY:
	    {
	      Real_wp coef = -alpha*omega2;
	      for (int k = 0; k < Nquad; k++)
		{
                  Real_wp coefd = poids(k)/square(rho(k)*c(k));
		  ProdC(p)(offset+k) += coef*coefd*
		    (real(Uh_quad(k))*real(Vh_quad(k))
		     + imag(Uh_quad(k))*imag(Vh_quad(k)));
		  
		  CopyVector(gradUh_quad, k, vec_ur, vec_ui);
		  CopyVector(gradVh_quad, k, vec_vr, vec_vi);
		  ProdC(p)(offset+k) += alpha*poids(k)/(rho(k)*rho(k))*(DotProd(vec_ur, vec_vr) + DotProd(vec_ui, vec_vi));
                  
                  ProdC(p)(offset+k) -= 2*alpha*this->omega(nk)*coefd*gamma(k)
                    *(real(Uh_quad(k))*imag(Vh_quad(k))
                      - imag(Uh_quad(k))*real(Vh_quad(k)));
                  
                  if (var_helm.include_flow_term)
                    ProdC(p)(offset+k) -= alpha*this->omega(nk)*poids(k) / square(rho(k)*c(k))
                      *(real(MdotGradU(k))*imag(Vh_quad(k))
                        - imag(MdotGradU(k))*real(Vh_quad(k))
                        + real(MdotGradV(k))*imag(Uh_quad(k))
                        - imag(MdotGradV(k))*real(Uh_quad(k)));

                  if (var_helm.include_double_gradM_flow)
                    ProdC(p)(offset+k) -= alpha*poids(k) / square(rho(k)*c(k))
                      *(real(MdotGradU(k))*real(MdotGradV(k))
                        + imag(MdotGradU(k))*imag(MdotGradV(k)));
                }
	    }
	    break;
	  case PARAM_C:
	    {
	      Real_wp coef = -alpha*omega2;
	      for (int k = 0; k < Nquad; k++)
		{
                  Real_wp coefc = 2*poids(k)/(rho(k)*c(k)*c(k)*c(k));
		  ProdC(p)(offset+k) += coef*coefc*
		    (real(Uh_quad(k))*real(Vh_quad(k))
		     + imag(Uh_quad(k))*imag(Vh_quad(k)));

                  ProdC(p)(offset+k) -= 2*alpha*this->omega(nk)*coefc*gamma(k)
                    *(real(Uh_quad(k))*imag(Vh_quad(k))
                      - imag(Uh_quad(k))*real(Vh_quad(k)));

                  if (var_helm.include_flow_term)
                    ProdC(p)(offset+k) -= alpha*this->omega(nk)*coefc
                      *(real(MdotGradU(k))*imag(Vh_quad(k))
                        - imag(MdotGradU(k))*real(Vh_quad(k))
                        + real(MdotGradV(k))*imag(Uh_quad(k))
                        - imag(MdotGradV(k))*real(Uh_quad(k)));
                  
                  if (var_helm.include_double_gradM_flow)
                    ProdC(p)(offset+k) -= alpha*coefc
                      *(real(MdotGradU(k))*real(MdotGradV(k))
                        + imag(MdotGradU(k))*imag(MdotGradV(k)));
		}
	    }
	    break;
          case PARAM_GAMMA:
            {
              Real_wp coef = 2*alpha*this->omega(nk);
              for (int k = 0; k < Nquad; k++)
                ProdC(p)(offset+k) += coef*poids(k)/(rho(k)*c(k)*c(k))
                  *(real(Uh_quad(k))*imag(Vh_quad(k))
                    - imag(Uh_quad(k))*real(Vh_quad(k)));
            }
            break;
          case PARAM_UX:
          case PARAM_UY:
          case PARAM_UZ:
            {
              int comp = this->param_c_type(p) - PARAM_UX;
              Real_wp coef = alpha*this->omega(nk);
              for (int k = 0; k < Nquad; k++)
                {
                  int off_k = Dimension::dim_N*k;
                  ProdC(p)(offset+k) += coef*poids(k)/(rho(k)*c(k)*c(k))*
                    (real(gradUh_quad(off_k+comp))*imag(Vh_quad(k))
                     - imag(gradUh_quad(off_k+comp))*real(Vh_quad(k))
                     + real(gradVh_quad(off_k+comp))*imag(Uh_quad(k))
                     - imag(gradVh_quad(off_k+comp))*real(Uh_quad(k)));

                  if (var_helm.include_double_gradM_flow)
                    ProdC(p)(offset+k) += alpha*poids(k)/(rho(k)*c(k)*c(k))*
                      (real(MdotGradU(k))*real(gradVh_quad(off_k+comp))
                       + imag(MdotGradU(k))*imag(gradVh_quad(off_k+comp))
                       + real(MdotGradV(k))*real(gradUh_quad(off_k+comp))
                       + imag(MdotGradV(k))*imag(gradUh_quad(off_k+comp)));
                }
            }
            break;
	  }
      }

    // contribution of absorbing boundary condition
    const MeshNumbering<Dimension>& mesh_num_helm = var_helm.GetMeshNumbering(0);
    Complex_wp coef = Real_wp(0.5)*alpha*Iwp*sqrt(omega2);
    for (int num_loc = 0; num_loc < var_helm.mesh.Element(i).GetNbBoundary(); num_loc++)
      {
	int nf = var_helm.mesh.Element(i).numBoundary(num_loc);
	int ref_line = var_helm.mesh.Boundary(nf).GetReference();
	if (var_helm.mesh.GetBoundaryCondition(ref_line) == BoundaryConditionEnum::LINE_ABSORBING)
	  {
	    nb_pts = mesh_num_helm.GetNbPointsQuadratureBoundary(nf);
	    poids.Reallocate(nb_pts);
	    for (int k = 0; k < nb_pts; k++)
	      poids(k) = var_helm.Glob_dsj(nf)(k)*Fb.WeightsQuadratureBoundary(k, num_loc);
	    
	    Uh_quad.Reallocate(nb_pts);
	    Vh_quad.Reallocate(nb_pts);
	    Fb.ApplyShTranspose(num_loc, Uh_dof(0), Uh_quad);
	    Fb.ApplyShTranspose(num_loc, Vh_dof(0), Vh_quad);
	    int of = Fb.GetNbPointsQuadratureInside();
	    for (int k = 0; k < num_loc; k++)
	      of += mesh_num_helm.GetNbPointsQuadratureBoundary(var_helm.mesh.Element(i).numBoundary(k));
	    
	    for (int p = 0; p < this->param_c_type.GetM(); p++)
	      {
		switch (this->param_c_type(p))
		  {
		  case PARAM_RHO:
		    {
		      for (int k = 0; k < nb_pts; k++)
			{
                          Complex_wp zeta = Iwp*gamma(of+k)/this->omega(nk);
			  Complex_wp vloc = coef*Uh_quad(k)*poids(k)*sqrt(c(of+k)/(rho(of+k)+zeta));
			  ProdC(p)(offset+of+k) += real(Vh_quad(k))*real(vloc) + imag(Vh_quad(k))*imag(vloc);
			}
		    }
		    break;
		  case PARAM_MU:
		    {
		      for (int k = 0; k < nb_pts; k++)
			{
			  Complex_wp zeta = Iwp*gamma(of+k)/this->omega(nk);
			  Complex_wp vloc = coef*Uh_quad(k)*poids(k)*sqrt((rho(of+k)+zeta)/c(of+k));
			  ProdC(p)(offset+of+k) += real(Vh_quad(k))*real(vloc) + imag(Vh_quad(k))*imag(vloc);
			}
		    }
		    break;
		  case PARAM_SIGMA:
		    {
		      for (int k = 0; k < nb_pts; k++)
			{
                          Complex_wp zeta = Iwp*gamma(of+k)/this->omega(nk);
			  Complex_wp vloc = -alpha/2*Uh_quad(k)*poids(k)*sqrt(c(of+k)/(rho(of+k)+zeta));
			  ProdC(p)(offset+of+k) += real(Vh_quad(k))*real(vloc) + imag(Vh_quad(k))*imag(vloc);
			}
		    }
		    break;
		  case PARAM_DENSITY:
		    {
		      for (int k = 0; k < nb_pts; k++)
			{
                          Complex_wp vloc = Real_wp(2)*poids(k)*Uh_quad(k)*coef/(square(rho(of+k))*c(of+k))
                            * sqrt(Real_wp(1) + Real_wp(2)*Iwp*gamma(of+k)/this->omega(nk));
			  ProdC(p)(offset+of+k) -= real(Vh_quad(k))*real(vloc) + imag(Vh_quad(k))*imag(vloc);
			}
		    }
		    break;
		  case PARAM_C:
		    {
		      for (int k = 0; k < nb_pts; k++)
			{
                          Complex_wp vloc = Real_wp(2)*poids(k)*Uh_quad(k)*coef/(square(c(of+k))*rho(of+k))
                            * sqrt(Real_wp(1) + Real_wp(2)*Iwp*gamma(of+k)/this->omega(nk));
			  ProdC(p)(offset+of+k) -= real(Vh_quad(k))*real(vloc) + imag(Vh_quad(k))*imag(vloc);
			}
		    }
		    break;
		  case PARAM_GAMMA:
		    {
		      for (int k = 0; k < nb_pts; k++)
			{
                          Complex_wp vloc = alpha*poids(k)*Uh_quad(k)/(c(of+k)*rho(of+k))
                            / sqrt(Real_wp(1) + Real_wp(2)*Iwp*gamma(of+k)/this->omega(nk));
			  ProdC(p)(offset+of+k) -= real(Vh_quad(k))*real(vloc) + imag(Vh_quad(k))*imag(vloc);
			}
		    }
		    break;
		  }
	      }
	  }
      }
  }



  /************************************************
   * BoundaryInverseProblem<HelmholtzEquationAxi> *
   ************************************************/

  
  //! Default constructor
  BoundaryInverseProblem<HelmholtzEquationAxi>::BoundaryInverseProblem()
    : BoundaryInverseProblem_Base<Dimension2>(var_helm, var_migration)
  {
    var_migration.axisymmetric_problem = true;
  }
  
  
  //! computation of density and speed of sound prior to inversion
  void BoundaryInverseProblem<HelmholtzEquationAxi>::InitIndices()
  {
    ref_density.Reallocate(var_helm.ref_rho.GetM());
    ref_c.Reallocate(var_helm.ref_rho.GetM());
    ref_rho.Reallocate(var_helm.ref_rho.GetM());
    ref_mu.Reallocate(var_helm.ref_rho.GetM());
    ref_sigma.Reallocate(var_helm.ref_rho.GetM());
    for (int ref = 0; ref < var_helm.ref_rho.GetM(); ref++)
      {
	Complex_wp rho = var_helm.ref_rho(ref).GetConstant();
	Complex_wp mu = var_helm.ref_mu(ref).GetConstant()(0, 0);
	ref_density(ref) = 1.0/mu;
	ref_c(ref) = sqrt(mu/rho);
	ref_rho(ref) = rho;
	ref_mu(ref) = mu;
	ref_sigma(ref) = var_helm.ref_sigma(ref).GetConstant();
      }
  }
  

  void BoundaryInverseProblem<HelmholtzEquationAxi>::
  SetInversionParameters(const VectString& param)
  {
    if (param.GetM() <= 0)
      {
        cout << "Please provide at least a parameter to invert" << endl;
        abort();
      }

    this->ResizeNbParameters(param.GetM());
    for (int k = 0; k < param.GetM(); k++)
      this->SetExactParameter(k, param(k), this->exact_fct_to_find[k]);
  }

  
  //! sets the exact parameter c
  void BoundaryInverseProblem<HelmholtzEquationAxi>::
  SetExactParameter(int k, const string& name, type_ptr_fctC fct)
  {
    BoundaryInverseProblem_Base<Dimension2>::SetExactParameter(k, name, fct);
    if (name == "rho")
      this->param_c_type(k) = PARAM_RHO;
    else if (name == "mu")
      this->param_c_type(k) = PARAM_MU;
    else if (name == "density")
      this->param_c_type(k) = PARAM_DENSITY;
    else if (name == "c")
      this->param_c_type(k) = PARAM_C;
    else
      {
	cout << "Unknown parameter" << endl;
	abort();
      }
  }
  

  //! returns the variable indexes in arrays rho_real and rho_cplx
  void BoundaryInverseProblem<HelmholtzEquationAxi>
  ::GetVariableIndex(int ref,
		     VarHarmonic_Base<Complex_wp, Dimension2>& var,
		     Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
		     Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_cplx)
  {
    EllipticProblem<HelmholtzEquationAxi>& var_
      = dynamic_cast<EllipticProblem<HelmholtzEquationAxi>& >(var);
    
    for (int num_index = 0; num_index < this->param_c_type.GetM(); num_index++)
      {
	if (this->param_c_type(num_index) == PARAM_RHO)
	  {
	    PhysicalVaryingMedia<Dimension2, Complex_wp>& rho 
	      = var_.ref_rho(ref).GetVaryingIndex();
	
	    var_.ref_rho(ref).SetMediaType(rho.USER);
	    rho_cplx.PushBack(&rho);
	  }
	else if (this->param_c_type(num_index) == PARAM_MU)
	  {
	    PhysicalVaryingMedia<Dimension2, Complex_wp>& mu 
	      = var_.ref_mu(ref).GetVaryingIndex(0, 0);
	    
	    var_.ref_mu(ref).SetMediaType(mu.USER);
	    rho_cplx.PushBack(&mu);
	  }
	else if (this->param_c_type(num_index) == PARAM_DENSITY)
	  {
	    PhysicalVaryingMedia<Dimension2, Complex_wp>& rho 
	      = var_.ref_rho(ref).GetVaryingIndex();
	    
	    PhysicalVaryingMedia<Dimension2, Complex_wp>& mu 
	      = var_.ref_mu(ref).GetVaryingIndex(0, 0);
	    
	    var_.ref_rho(ref).SetMediaType(rho.USER);
	    var_.ref_mu(ref).SetMediaType(mu.USER);
	    rho_cplx.PushBack(&rho);
	    rho_cplx.PushBack(&mu);
	  }
	else if (this->param_c_type(num_index) == PARAM_C)
	  {
	    PhysicalVaryingMedia<Dimension2, Complex_wp>& rho 
	      = var_.ref_rho(ref).GetVaryingIndex();
	    
	    var_.ref_rho(ref).SetMediaType(rho.USER);
	    rho_cplx.PushBack(&rho);
	  }
	else
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
      }
  }
  
  
  //! Fills variables indexes on element i1, point j from values contained in eval_c
  void BoundaryInverseProblem<HelmholtzEquationAxi>
	::SetVariableIndexPoint(int ref, int i, int i1, int j, int offset, const Vector<VectReal_wp>& eval_c,
			  Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
			  Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_cplx)
  {
    int pd = -1, pc = -1;
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	if (this->param_c_type(p) == PARAM_DENSITY)
	  pd = p;
	else if (this->param_c_type(p) == PARAM_C)
	  pc = p;
      }
    
    int num = 0;
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	if (this->param_c_type(p) == PARAM_RHO)
	  {
	    rho_cplx(num)->SetCoefficient(i1, j, eval_c(p)(offset+j));
	    num++;
	  }
	else if (this->param_c_type(p) == PARAM_MU)
	  {
	    rho_cplx(num)->SetCoefficient(i1, j, eval_c(p)(offset+j));
	    num++;
	  }
	else if (this->param_c_type(p) == PARAM_DENSITY)
	  {
	    Complex_wp c = ref_c(ref), rho = eval_c(p)(offset+j);
	    if (pc >= 0)
	      c = eval_c(pc)(offset+j);
	    
	    rho_cplx(num)->SetCoefficient(i1, j, Complex_wp(1, 0)/(rho*c*c));
	    rho_cplx(num+1)->SetCoefficient(i1, j, Complex_wp(1, 0)/rho);
	    num += 2;
	  }
	else if (this->param_c_type(p) == PARAM_C)
	  {
	    Complex_wp c = eval_c(p)(offset+j), rho = ref_density(ref);
	    if (pd >= 0)
	      rho = eval_c(pd)(offset+j);
	    
	    rho_cplx(num)->SetCoefficient(i1, j, Complex_wp(1, 0)/(rho*c*c));
	    num++;
	  }
      }	
  }
  
  
  //! local product by derivative with respect to parameters 
  void BoundaryInverseProblem<HelmholtzEquationAxi>
  ::MltParamDerivLocal(int nk, const Real_wp& omega2, VectComplex_wp& u_glob,
		       int i, const Real_wp& alpha, const Vector<VectReal_wp>& c_quad,
		       int offset, VectComplex_wp& y_col)
  {
    bool mass = false, stiff = false;
    int pd = -1, pc = -1, prho = -1, pmu = -1;
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO : mass = true; prho = p; break;
	  case PARAM_MU : stiff = true; pmu = p; break;
	  case PARAM_DENSITY : mass = true; stiff = true; pd = p; break;
	  case PARAM_C : mass = true; pc = p; break;
	  }
      }
    
    TinyVector<Complex_wp, 2> vec_u, vec_v;
    VectComplex_wp gradUh_quad, Uh_quad, feval, grad_feval, Uh_val;
    TinyVector<VectComplex_wp, 1> Uh_dof;
    const ElementReference<Dimension2, 1>& Fb = var_helm.GetReferenceElementH1(i);
    var_helm.GetLocalUnknownVector(u_glob, i, Uh_dof);
    
    Uh_val = Uh_dof(0);
    // values of u (and/or grad u) on quadrature points are computed
    // Uh_quad contains u(xi_k) and gradUh_quad contains \nabla u(\xi_k)
    int Nquad = Fb.GetNbPointsQuadratureInside();
    Uh_quad.Reallocate(Nquad);
    if (stiff)
      gradUh_quad.Reallocate(Nquad*2);
    
    if (Fb.UseQuadraturePointsForRh())
      {
	Fb.ApplyChTranspose(Uh_dof(0), Uh_quad);
	if (stiff)
	  Fb.ApplyRhQuadratureTranspose(Uh_quad, gradUh_quad);
      }
    else
      {
	if (mass)
	  Fb.ApplyChTranspose(Uh_dof(0), Uh_quad);
	
	if (stiff)
	  Fb.ApplyRhTranspose(Uh_dof(0), gradUh_quad);
      }
    
    // we retrieve weight of integration
    VectReal_wp poids(Nquad);
    Matrix2_2 dfjm1;
    Vector<Matrix2_2sym> matB;
    if (stiff)
      matB.Reallocate(Nquad);
    
    bool affine = var_helm.mesh.IsElementAffine(i);
    if (affine)
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(0)*Fb.WeightsND(k)*real(var_helm.Glob_rtilde(i)(k));
	  if (stiff)
	    {
	      dfjm1 = var_helm.Glob_DFjm1(i)(0);
	      Mlt(1.0/var_helm.Glob_jacobian(i)(0), dfjm1);
	      MltTrans(dfjm1, dfjm1, matB(k));
	    }
	} 
    else
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(k)*real(var_helm.Glob_rtilde(i)(k));
	  if (stiff)
	    {
	      dfjm1 = var_helm.Glob_DFjm1(i)(k);
	      Mlt(Fb.WeightsND(k)/var_helm.Glob_jacobian(i)(k), dfjm1);
	      MltTrans(dfjm1, dfjm1, matB(k));
	    }
	}
    
    feval.Reallocate(Nquad);
    feval.Zero();
    
    if (stiff)
      {
	grad_feval.Reallocate(Nquad*2);
	grad_feval.Zero();
      }

    // we update feval and grad_feval by using the parameters
    int ref = var_helm.mesh.Element(i).GetReference();
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO:
	    {
	      Real_wp coef = alpha*omega2;
	      for (int k = 0; k < Nquad; k++)
		feval(k) += coef*poids(k)*c_quad(p)(offset+k)*Uh_quad(k);
	    }
	    break;
	  case PARAM_MU:
	    {
	      for (int k = 0; k < Nquad; k++)
		{
		  CopyVector(gradUh_quad, k, vec_u);
		  vec_u *= -alpha*poids(k)*c_quad(p)(offset+k);
		  Mlt(matB(k), vec_u, vec_v);
		  AddVector(vec_v, 2*k, grad_feval);
		}
	    }
	    break;
	  case PARAM_DENSITY:
	    {
	      Real_wp coef = -alpha*omega2;
	      Real_wp c = realpart(ref_c(ref)), rho;
	      for (int k = 0; k < Nquad; k++)
		{
		  if (pc >= 0)
		    c = this->eval_c_quad(pc)(offset+k);
		  
		  rho = this->eval_c_quad(p)(offset+k);
		  
		  feval(k) += coef*poids(k)/square(rho*c)*c_quad(p)(offset+k)*Uh_quad(k);
		  
		  CopyVector(gradUh_quad, k, vec_u);
		  vec_u *= alpha/(rho*rho)*poids(k)*c_quad(p)(offset+k);
		  Mlt(matB(k), vec_u, vec_v);
		  AddVector(vec_v, 2*k, grad_feval);
		}	      
	    }
	    break;
	  case PARAM_C:
	    {
	      Real_wp coef = -2*alpha*omega2;
	      Real_wp rho = realpart(ref_density(ref)), c;
	      for (int k = 0; k < Nquad; k++)
		{
		  if (pd >= 0)
		    rho = this->eval_c_quad(pd)(offset+k);
		  
		  c = this->eval_c_quad(p)(offset+k);
		  feval(k) += coef*poids(k)/(rho*c*c*c)*c_quad(p)(offset+k)*Uh_quad(k);
		}
	    }
	    break;
	  }
      }
	
    // y_col is updated with feval and grad_feval
    if (Fb.UseQuadraturePointsForRh())
      {
	if (stiff)
	  {
	    Fb.ApplyRhQuadrature(grad_feval, Uh_dof(0));
	    feval = feval + Uh_dof(0);
	  }
	
	Fb.ApplyCh(feval, Uh_dof(0));
      }
    else
      {	
	if (mass)
	  Fb.ApplyCh(feval, Uh_dof(0));
	else
	  Uh_dof(0).Zero();
	
	if (stiff)
	  {
	    Uh_quad.Reallocate(Uh_dof(0).GetM());
	    Fb.ApplyRh(grad_feval, Uh_quad);
	    Uh_dof(0) = Uh_dof(0) + Uh_quad;
	  }	
      }
    
    // contribution of absorbing boundary condition
    Complex_wp coef = Real_wp(0.5)*alpha*Iwp*sqrt(omega2);
    Complex_wp SigmaOm = Iwp*ref_sigma(ref)/sqrt(omega2);
    const MeshNumbering<Dimension2>& mesh_num_helm = var_helm.GetMeshNumbering(0);
    for (int num_loc = 0; num_loc < var_helm.mesh.Element(i).GetNbBoundary(); num_loc++)
      {
	int nf = var_helm.mesh.Element(i).numBoundary(num_loc);
	int ref_line = var_helm.mesh.Boundary(nf).GetReference();
	if (var_helm.mesh.GetBoundaryCondition(ref_line) == BoundaryConditionEnum::LINE_ABSORBING)
	  {
	    int nb_pts = mesh_num_helm.GetNbPointsQuadratureBoundary(nf);
	    poids.Reallocate(nb_pts);
	    
	    Uh_quad.Reallocate(nb_pts);
	    feval.Reallocate(nb_pts); feval.Zero();
	    Fb.ApplyShTranspose(num_loc, Uh_val, Uh_quad);
	    int offset2 = offset + Fb.GetNbPointsQuadratureInside();
	    for (int k = 0; k < num_loc; k++)
	      offset2 += mesh_num_helm.GetNbPointsQuadratureBoundary(var_helm.mesh.Element(i).numBoundary(k));
	    
	    for (int k = 0; k < nb_pts; k++)
	      poids(k) = var_helm.Glob_dsj(nf)(k)*Fb.WeightsQuadratureBoundary(k, num_loc)
		*real(var_helm.Glob_rtilde(i)(offset2-offset+k));

	    for (int p = 0; p < this->param_c_type.GetM(); p++)
	      {
		switch (this->param_c_type(p))
		  {
		  case PARAM_RHO:
		    {
		      Real_wp mu = realpart(ref_mu(ref)), rho;
		      for (int k = 0; k < nb_pts; k++)
			{
			  if (pmu >= 0)
			    mu = this->eval_c_quad(pmu)(offset2+k);
			  
			  rho = this->eval_c_quad(p)(offset2+k);
			  
			  feval(k) += coef*Uh_quad(k)*poids(k)*sqrt(mu/(rho+SigmaOm))*c_quad(p)(offset2+k);
			}
		    }
		    break;
		  case PARAM_MU:
		    {
		      Real_wp rho = realpart(ref_rho(ref)), mu;
		      for (int k = 0; k < nb_pts; k++)
			{
			  if (prho >= 0)
			    rho = this->eval_c_quad(prho)(offset2+k);
			  
			  mu = this->eval_c_quad(p)(offset2+k);
			  
			  feval(k) += coef*Uh_quad(k)*poids(k)*sqrt((rho+SigmaOm)/mu)*c_quad(p)(offset2+k);
			}
		    }
		    break;
		  case PARAM_DENSITY:
		    {
		      Real_wp c = realpart(ref_c(ref)), rho;
		      for (int k = 0; k < nb_pts; k++)
			{
			  if (pc >= 0)
			    c = this->eval_c_quad(pc)(offset2+k);
			  
			  rho = this->eval_c_quad(p)(offset2+k);
			  
			  feval(k) -= Real_wp(2)*coef*Uh_quad(k)*poids(k)/(rho*rho*c)*c_quad(p)(offset2+k);
			}
		    }
		    break;		    
		  case PARAM_C:
		    {
		      Real_wp rho = realpart(ref_density(ref)), c;
		      for (int k = 0; k < nb_pts; k++)
			{
			  if (pd >= 0)
			    rho = this->eval_c_quad(pd)(offset2+k);
			  
			  c = this->eval_c_quad(p)(offset2+k);
			  
			  feval(k) -= Real_wp(2)*coef*Uh_quad(k)*poids(k)/(rho*c*c)*c_quad(p)(offset2+k);
			}
		    }
		    break;		    
		  }
	      }
	    
	    Fb.ApplySh(Complex_wp(1, 0), num_loc, feval, Uh_dof(0), Fb.GetOrder());
	  }
      }
    
    var_helm.AddLocalUnknownVector(Complex_wp(1, 0), Uh_dof, i, y_col);
  }


  //! local product by transpose of derivative with respect to parameters   
  void BoundaryInverseProblem<HelmholtzEquationAxi>
  ::MltParamDerivTransLoc(int nk, const Real_wp& omega2, VectComplex_wp& u_glob,
			  int i, const Real_wp& alpha, int offset, 
			  const VectComplex_wp& y_col, Vector<VectReal_wp>& ProdC)
  {
    bool mass = false, stiff = false;
    int pd = -1, pc = -1, prho = -1, pmu = -1;
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO : mass = true; prho = p; break;
	  case PARAM_MU : stiff = true; pmu = p; break;
	  case PARAM_DENSITY : mass = true; stiff = true; pd = p; break;
	  case PARAM_C : mass = true; pc = p; break;
	  }
      }

    VectComplex_wp Uh_quad, Vh_quad, gradUh_quad, gradVh_quad;
    TinyVector<VectComplex_wp, 1> Uh_dof, Vh_dof;
    const ElementReference<Dimension2, 1>& Fb = var_helm.GetReferenceElementH1(i);
    var_helm.GetLocalUnknownVector(u_glob, i, Uh_dof);
    
    var_helm.GetLocalUnknownVector(y_col, i, Vh_dof);
    
    // values of u and v are computed on quadrature points (and gradient if needed)
    int Nquad = Fb.GetNbPointsQuadratureInside();
    Uh_quad.Reallocate(Nquad);
    Vh_quad.Reallocate(Nquad);
    
    if (stiff)
      {
	gradUh_quad.Reallocate(2*Nquad);
	gradVh_quad.Reallocate(2*Nquad);
      }
    
    if (Fb.UseQuadraturePointsForRh())
      {
	Fb.ApplyChTranspose(Uh_dof(0), Uh_quad);
	Fb.ApplyChTranspose(Vh_dof(0), Vh_quad);	
	if (stiff)
	  {
	    Fb.ApplyRhQuadratureTranspose(Uh_quad, gradUh_quad);
	    Fb.ApplyRhQuadratureTranspose(Vh_quad, gradVh_quad);
	  }
      }
    else
      {
	if (mass)
	  {
	    Fb.ApplyChTranspose(Uh_dof(0), Uh_quad);
	    Fb.ApplyChTranspose(Vh_dof(0), Vh_quad);
	  }
	
	if (stiff)
	  {
	    Fb.ApplyRhTranspose(Uh_dof(0), gradUh_quad);
	    Fb.ApplyRhTranspose(Vh_dof(0), gradVh_quad);
	  }
      }
    
    // we retrieve weight of integration
    VectReal_wp poids(Nquad);
    Matrix2_2 dfjm1;
    Vector<TinyMatrix<Real_wp, Symmetric, 2, 2> > matB;
    if (stiff)
      matB.Reallocate(Nquad);

    bool affine = var_helm.mesh.IsElementAffine(i);
    if (affine)
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(0)*Fb.WeightsND(k)*real(var_helm.Glob_rtilde(i)(k));
	  if (stiff)
	    {
	      dfjm1 = var_helm.Glob_DFjm1(i)(0);
	      Mlt(1.0/var_helm.Glob_jacobian(i)(0), dfjm1);
	      MltTrans(dfjm1, dfjm1, matB(k));
	    }
	} 
    else
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(k)*real(var_helm.Glob_rtilde(i)(k));
	  if (stiff)
	    {
	      dfjm1 = var_helm.Glob_DFjm1(i)(k);
	      Mlt(Fb.WeightsND(k)/var_helm.Glob_jacobian(i)(k), dfjm1);
	      MltTrans(dfjm1, dfjm1, matB(k));
	    }
	}
    
    // we update ProdC
    TinyVector<Real_wp, 2> vec_ur, vec_ui, vec_vr, vec_vi, vec_ur2, vec_ui2;
    int ref = var_helm.mesh.Element(i).GetReference();
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO:
	    {
	      Real_wp coef = alpha*omega2;
	      for (int k = 0; k < Nquad; k++)
		ProdC(p)(offset+k) += coef*poids(k)*(real(Uh_quad(k))*real(Vh_quad(k))
						     + imag(Uh_quad(k))*imag(Vh_quad(k)));
	    }
	    break;
	  case PARAM_MU:
	    {
	      Real_wp coef = -alpha;
	      for (int k = 0; k < Nquad; k++)
		{
		  CopyVector(gradUh_quad, k, vec_ur, vec_ui);
		  Mlt(matB(k), vec_ur, vec_ur2);
		  Mlt(matB(k), vec_ui, vec_ui2);
		  
		  CopyVector(gradVh_quad, k, vec_vr, vec_vi);
		  ProdC(p)(offset+k) += coef*poids(k)*(DotProd(vec_ur2, vec_vr) + DotProd(vec_ui2, vec_vi));
		}
	    }
	    break;
	  case PARAM_DENSITY:
	    {
	      Real_wp coef = -alpha*omega2;
	      Real_wp c = realpart(ref_c(ref)), rho;
	      for (int k = 0; k < Nquad; k++)
		{
		  if (pc >= 0)
		    c = this->eval_c_quad(pc)(offset+k);
		  
		  rho = this->eval_c_quad(p)(offset+k);
		  
		  ProdC(p)(offset+k) += coef*poids(k)/square(rho*c)*
		    (real(Uh_quad(k))*real(Vh_quad(k))
		     + imag(Uh_quad(k))*imag(Vh_quad(k)));
		  
		  CopyVector(gradUh_quad, k, vec_ur, vec_ui);
		  Mlt(matB(k), vec_ur, vec_ur2);
		  Mlt(matB(k), vec_ui, vec_ui2);
		  
		  CopyVector(gradVh_quad, k, vec_vr, vec_vi);
		  ProdC(p)(offset+k) += alpha*poids(k)/(rho*rho)*(DotProd(vec_ur2, vec_vr) + DotProd(vec_ui2, vec_vi));
		}
	    }
	    break;
	  case PARAM_C:
	    {
	      Real_wp coef = -2*alpha*omega2;
	      Real_wp rho = realpart(ref_density(ref)), c;
	      for (int k = 0; k < Nquad; k++)
		{
		  if (pd >= 0)
		    rho = this->eval_c_quad(pd)(offset+k);
		  
		  c = this->eval_c_quad(p)(offset+k);
		  
		  ProdC(p)(offset+k) += coef*poids(k)/(rho*c*c*c)*
		    (real(Uh_quad(k))*real(Vh_quad(k))
		     + imag(Uh_quad(k))*imag(Vh_quad(k)));
		}
	    }
	    break;
	  }
      }

    // contribution of absorbing boundary condition
    Complex_wp coef = Real_wp(0.5)*alpha*Iwp*sqrt(omega2);
    Complex_wp SigmaOm = Iwp*ref_sigma(ref)/sqrt(omega2);
    const MeshNumbering<Dimension2>& mesh_num_helm = var_helm.GetMeshNumbering(0);
    for (int num_loc = 0; num_loc < var_helm.mesh.Element(i).GetNbBoundary(); num_loc++)
      {
	int nf = var_helm.mesh.Element(i).numBoundary(num_loc);
	int ref_line = var_helm.mesh.Boundary(nf).GetReference();
	if (var_helm.mesh.GetBoundaryCondition(ref_line) == BoundaryConditionEnum::LINE_ABSORBING)
	  {
	    int nb_pts = mesh_num_helm.GetNbPointsQuadratureBoundary(nf);
	    poids.Reallocate(nb_pts);
	    
	    Uh_quad.Reallocate(nb_pts);
	    Vh_quad.Reallocate(nb_pts);
	    Fb.ApplyShTranspose(num_loc, Uh_dof(0), Uh_quad);
	    Fb.ApplyShTranspose(num_loc, Vh_dof(0), Vh_quad);
	    int offset2 = offset + Fb.GetNbPointsQuadratureInside();
	    for (int k = 0; k < num_loc; k++)
	      offset2 += mesh_num_helm.GetNbPointsQuadratureBoundary(var_helm.mesh.Element(i).numBoundary(k));

	    for (int k = 0; k < nb_pts; k++)
	      poids(k) = var_helm.Glob_dsj(nf)(k)*Fb.WeightsQuadratureBoundary(k, num_loc)
		*real(var_helm.Glob_rtilde(i)(offset2-offset+k));
	    
	    for (int p = 0; p < this->param_c_type.GetM(); p++)
	      {
		switch (this->param_c_type(p))
		  {
		  case PARAM_RHO:
		    {
		      Real_wp mu = realpart(ref_mu(ref)), rho;
		      for (int k = 0; k < nb_pts; k++)
			{
			  if (pmu >= 0)
			    mu = this->eval_c_quad(pmu)(offset2+k);
			  
			  rho = this->eval_c_quad(p)(offset2+k);
			  
			  Complex_wp vloc = coef*Uh_quad(k)*poids(k)*sqrt(mu/(rho+SigmaOm));
			  ProdC(p)(offset2+k) += real(Vh_quad(k))*real(vloc) + imag(Vh_quad(k))*imag(vloc);
			}
		    }
		    break;
		  case PARAM_MU:
		    {
		      Real_wp rho = realpart(ref_rho(ref)), mu;
		      for (int k = 0; k < nb_pts; k++)
			{
			  if (prho >= 0)
			    rho = this->eval_c_quad(prho)(offset2+k);
			  
			  mu = this->eval_c_quad(p)(offset2+k);
			  
			  Complex_wp vloc = coef*Uh_quad(k)*poids(k)*sqrt((rho+SigmaOm)/mu);
			  ProdC(p)(offset2+k) += real(Vh_quad(k))*real(vloc) + imag(Vh_quad(k))*imag(vloc);
			}
		    }
		    break;
		  case PARAM_DENSITY:
		    {
		      Real_wp c = realpart(ref_c(ref)), rho;
		      for (int k = 0; k < nb_pts; k++)
			{
			  if (pc >= 0)
			    c = this->eval_c_quad(pc)(offset2+k);
			  
			  rho = this->eval_c_quad(p)(offset2+k);
			  
			  Complex_wp vloc = Real_wp(2)*coef*Uh_quad(k)*poids(k)/(rho*rho*c);
			  ProdC(p)(offset2+k) -= real(Vh_quad(k))*real(vloc) + imag(Vh_quad(k))*imag(vloc);
			}
		    }
		    break;
		  case PARAM_C:
		    {
		      Real_wp rho = realpart(ref_density(ref)), c;
		      for (int k = 0; k < nb_pts; k++)
			{
			  if (pd >= 0)
			    rho = this->eval_c_quad(pd)(offset2+k);
			  
			  c = this->eval_c_quad(p)(offset2+k);
			  
			  Complex_wp vloc = Real_wp(2)*coef*Uh_quad(k)*poids(k)/(rho*c*c);
			  ProdC(p)(offset2+k) -= real(Vh_quad(k))*real(vloc) + imag(Vh_quad(k))*imag(vloc);
			}
		    }
		    break;		    
		  }
	      }
	  }
      }

  }

  
  /***************************************************************
   * BoundaryInverseProblem<HarmonicElasticEquation<Dimension> > *
   ***************************************************************/
  

  //! default constructor
  template<class Dimension>
  BoundaryInverseProblem<HarmonicElasticEquation<Dimension> >::BoundaryInverseProblem()
    : BoundaryInverseProblem_Base<Dimension>(var_helm, var_migration)
  {
  }


  //! computation of density and speed of sound prior to inversion
  template<class Dimension>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension> >::InitIndices()
  {
    ref_density.Reallocate(var_helm.ref_rho.GetM());
    ref_C.Reallocate(var_helm.ref_rho.GetM());
    for (int ref = 0; ref < var_helm.ref_rho.GetM(); ref++)
      {
	ref_density(ref) = var_helm.ref_rho(ref).GetConstant();
	ref_C(ref) = var_helm.ref_tensorC(ref).GetConstant();
      }
  }


#ifdef MONTJOIE_WITH_TWO_DIM
  // isotropic media (elasto)
  template<>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension2> >::
  SetInversionParameters(const VectString& param)
  {
    if (param(0) == "Isotrope")
      {
        this->ResizeNbParameters(4);
        this->SetExactParameter(0, "lambda", this->exact_fct_to_find[0]);
        this->SetExactParameter(1, "mu", this->exact_fct_to_find[1]);
        this->SetExactParameter(2, "lambda_i", this->exact_fct_to_find[2]);
        this->SetExactParameter(3, "mu_i", this->exact_fct_to_find[3]);
      }
    else if (param(0) == "Orthotrope")
      {
        // orthotropic media (elasto)
        this->ResizeNbParameters(8);
        this->SetExactParameter(0, "c11", this->exact_fct_to_find[0]);
        this->SetExactParameter(1, "c12", this->exact_fct_to_find[1]);
        this->SetExactParameter(2, "c22", this->exact_fct_to_find[2]);
        this->SetExactParameter(3, "c66", this->exact_fct_to_find[3]);
        this->SetExactParameter(4, "c11_i", this->exact_fct_to_find[4]);
        this->SetExactParameter(5, "c12_i", this->exact_fct_to_find[5]);
        this->SetExactParameter(6, "c22_i", this->exact_fct_to_find[6]);
        this->SetExactParameter(7, "c66_i", this->exact_fct_to_find[7]);
      }
    else if (param(0) == "Anisotrope")
      {
        // anisotropic media (elasto)
        this->ResizeNbParameters(12);
        this->SetExactParameter(0, "c11", this->exact_fct_to_find[0]);
        this->SetExactParameter(1, "c12", this->exact_fct_to_find[1]);
        this->SetExactParameter(2, "c16", this->exact_fct_to_find[2]);
        this->SetExactParameter(3, "c22", this->exact_fct_to_find[3]);
        this->SetExactParameter(4, "c26", this->exact_fct_to_find[4]);
        this->SetExactParameter(5, "c66", this->exact_fct_to_find[5]);
        this->SetExactParameter(6, "c11_i", this->exact_fct_to_find[6]);
        this->SetExactParameter(7, "c12_i", this->exact_fct_to_find[7]);
        this->SetExactParameter(8, "c16_i", this->exact_fct_to_find[8]);
        this->SetExactParameter(9, "c22_i", this->exact_fct_to_find[9]);
        this->SetExactParameter(10, "c26_i", this->exact_fct_to_find[10]);
        this->SetExactParameter(11, "c66_i", this->exact_fct_to_find[11]);
      }
    else
      {
        cout << "Parameter not known : " << param << endl;
        cout << "TypeInversion = Isotrope, Orthotrope or Anisotrope" << endl;
        abort();
      }

    if (param.GetM() > 1)
      {
        if (param(1) == "Force")
          this->nb_comp_source_invert = 2;
      }

    if (param.GetM() > 2)
      {
        if (param(2) == "MeasureAngle")
          {
            this->nb_angles_measure = 1;
            this->nb_angles_to_invert = 1;
          }
      }
  }

  //! linear combination of data
  template<>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension2> >
  ::CombineData(const Matrix<Complex_wp, General, ColMajor>& data0,
                Matrix<Complex_wp, General, ColMajor>& data)
  {
    if (this->nb_angles_measure == 1)
      {
        Real_wp c1 = coef_combination_measure(0);
        Real_wp c2 = coef_combination_measure(1);
        data.Reallocate(data0.GetM(), data0.GetN()/2);
        for (int i = 0; i < data.GetM(); i++)
          for (int j = 0; j < data.GetN(); j++)
            data(i, j) = data0(i, 2*j)*c1 + data0(i, 2*j+1)*c2;
        
        return;
      }
    
    // default choice : no combination
    data = data0;
  }


  //! linear combination of data (and derivative)
  template<>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension2> >
  ::CombineDeriveData(const VectReal_wp& h0, int k,
                      const Matrix<Complex_wp, General, ColMajor>& data_simu,
                      Matrix<Complex_wp, General, ColMajor>& data_comb)
  {
    if (this->nb_angles_measure == 1)
      {
        int offset = 0;
        for (int p = 0; p < this->c_discrete.GetM(); p++)
          offset += this->c_discrete(p)->GetNbDof();
        
        offset += this->nb_comp_source_invert;
        Real_wp teta_m(0);
        if (this->nb_angles_to_invert == 1)
          teta_m = h0(offset)*pi_wp/180;
        
        Real_wp c1 = this->coef_combination_measure(0);
        Real_wp c2 = this->coef_combination_measure(1);
        data_comb.Reallocate(data_simu.GetM(), data_simu.GetN()/2);
        for (int i = 0; i < data_comb.GetM(); i++)
          for (int j = 0; j < data_comb.GetN(); j++)
            {
              Complex_wp ux = u_observable(k)(i, 2*j), uy = u_observable(k)(i, 2*j+1);
              data_comb(i, j) = data_simu(i, 2*j)*c1 + data_simu(i, 2*j+1)*c2;
              data_comb(i, j) += (-c2*ux + c1*uy)*teta_m;
            }

        return;
      }

    // default choice : no combination
    data_comb = data_simu;
  }


  //! linear combination of data (and adjoint of derivative)
  template<>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension2> >
  ::CombineAdjointData(int k, const Matrix<Complex_wp, General, ColMajor>& data_simu,
                       Matrix<Complex_wp, General, ColMajor>& data_comb, VectReal_wp& h0)
  {
    if (this->nb_angles_measure == 1)
      {
        int offset = 0;
        for (int p = 0; p < this->c_discrete.GetM(); p++)
          offset += this->c_discrete(p)->GetNbDof();
        
        offset += this->nb_comp_source_invert;
        Real_wp c1 = this->coef_combination_measure(0);
        Real_wp c2 = this->coef_combination_measure(1);
        data_comb.Reallocate(data_simu.GetM(), data_simu.GetN()*2);
        for (int i = 0; i < data_simu.GetM(); i++)
          for (int j = 0; j < data_simu.GetN(); j++)
            {
              Complex_wp ux = u_observable(k)(i, 2*j), uy = u_observable(k)(i, 2*j+1);
              data_comb(i, 2*j) = data_simu(i, j)*c1;
              data_comb(i, 2*j+1) = data_simu(i, j)*c2;
              if (this->nb_angles_to_invert == 1)
                h0(offset) += (realpart(-c2*ux + c1*uy)*realpart(data_simu(i, j))
                               + imagpart(-c2*ux + c1*uy)*imagpart(data_simu(i, j)))*pi_wp/180;
            }
        
        return;
      }

    // default choice : no combination
    data_comb = data_simu;
  }
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
  // isotropic media (elasto)
  template<>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension3> >::
  SetInversionParameters(const VectString& param)
  {
    if (param(0) == "Isotrope")
      {
        this->ResizeNbParameters(4);
        this->SetExactParameter(0, "lambda", this->exact_fct_to_find[0]);
        this->SetExactParameter(1, "mu", this->exact_fct_to_find[1]);
        this->SetExactParameter(2, "lambda_i", this->exact_fct_to_find[2]);
        this->SetExactParameter(3, "mu_i", this->exact_fct_to_find[3]);
      }
    else if (param(0) == "Orthotrope")
      {
        // orthotropic media (elasto)
        this->ResizeNbParameters(18);
        this->SetExactParameter(0, "c11", this->exact_fct_to_find[0]);
        this->SetExactParameter(1, "c12", this->exact_fct_to_find[1]);
        this->SetExactParameter(2, "c13", this->exact_fct_to_find[2]);
        this->SetExactParameter(3, "c22", this->exact_fct_to_find[3]);
        this->SetExactParameter(4, "c23", this->exact_fct_to_find[4]);
        this->SetExactParameter(5, "c33", this->exact_fct_to_find[5]);
        this->SetExactParameter(6, "c44", this->exact_fct_to_find[6]);
        this->SetExactParameter(7, "c55", this->exact_fct_to_find[7]);
        this->SetExactParameter(8, "c66", this->exact_fct_to_find[8]);
        this->SetExactParameter(9, "c11_i", this->exact_fct_to_find[9]);
        this->SetExactParameter(10, "c12_i", this->exact_fct_to_find[10]);
        this->SetExactParameter(11, "c13_i", this->exact_fct_to_find[11]);
        this->SetExactParameter(12, "c22_i", this->exact_fct_to_find[12]);
        this->SetExactParameter(13, "c23_i", this->exact_fct_to_find[13]);
        this->SetExactParameter(14, "c33_i", this->exact_fct_to_find[14]);
        this->SetExactParameter(15, "c44_i", this->exact_fct_to_find[15]);
        this->SetExactParameter(16, "c55_i", this->exact_fct_to_find[16]);
        this->SetExactParameter(17, "c66_i", this->exact_fct_to_find[17]);
      }
    else if (param(0) == "Anisotrope")
      {
        // anisotropic media (elasto)
        this->ResizeNbParameters(42);
        this->SetExactParameter(0, "c11", this->exact_fct_to_find[0]);
        this->SetExactParameter(1, "c12", this->exact_fct_to_find[1]);
        this->SetExactParameter(2, "c13", this->exact_fct_to_find[2]);
        this->SetExactParameter(3, "c14", this->exact_fct_to_find[3]);
        this->SetExactParameter(4, "c15", this->exact_fct_to_find[4]);
        this->SetExactParameter(5, "c16", this->exact_fct_to_find[5]);
        this->SetExactParameter(6, "c22", this->exact_fct_to_find[6]);
        this->SetExactParameter(7, "c23", this->exact_fct_to_find[7]);
        this->SetExactParameter(8, "c24", this->exact_fct_to_find[8]);
        this->SetExactParameter(9, "c25", this->exact_fct_to_find[9]);
        this->SetExactParameter(10, "c26", this->exact_fct_to_find[10]);
        this->SetExactParameter(11, "c33", this->exact_fct_to_find[11]);
        this->SetExactParameter(12, "c34", this->exact_fct_to_find[12]);
        this->SetExactParameter(13, "c35", this->exact_fct_to_find[13]);
        this->SetExactParameter(14, "c36", this->exact_fct_to_find[14]);
        this->SetExactParameter(15, "c44", this->exact_fct_to_find[15]);
        this->SetExactParameter(16, "c45", this->exact_fct_to_find[16]);
        this->SetExactParameter(17, "c46", this->exact_fct_to_find[17]);
        this->SetExactParameter(18, "c55", this->exact_fct_to_find[18]);
        this->SetExactParameter(19, "c56", this->exact_fct_to_find[19]);
        this->SetExactParameter(20, "c66", this->exact_fct_to_find[20]);
        this->SetExactParameter(21, "c11_i", this->exact_fct_to_find[21]);
        this->SetExactParameter(22, "c12_i", this->exact_fct_to_find[22]);
        this->SetExactParameter(23, "c13_i", this->exact_fct_to_find[23]);
        this->SetExactParameter(24, "c14_i", this->exact_fct_to_find[24]);
        this->SetExactParameter(25, "c15_i", this->exact_fct_to_find[25]);
        this->SetExactParameter(26, "c16_i", this->exact_fct_to_find[26]);
        this->SetExactParameter(27, "c22_i", this->exact_fct_to_find[27]);
        this->SetExactParameter(28, "c23_i", this->exact_fct_to_find[28]);
        this->SetExactParameter(29, "c24_i", this->exact_fct_to_find[29]);
        this->SetExactParameter(30, "c25_i", this->exact_fct_to_find[30]);
        this->SetExactParameter(31, "c26_i", this->exact_fct_to_find[31]);
        this->SetExactParameter(32, "c33_i", this->exact_fct_to_find[32]);
        this->SetExactParameter(33, "c34_i", this->exact_fct_to_find[33]);
        this->SetExactParameter(34, "c35_i", this->exact_fct_to_find[34]);
        this->SetExactParameter(35, "c36_i", this->exact_fct_to_find[35]);
        this->SetExactParameter(36, "c44_i", this->exact_fct_to_find[36]);
        this->SetExactParameter(37, "c45_i", this->exact_fct_to_find[37]);
        this->SetExactParameter(38, "c46_i", this->exact_fct_to_find[38]);
        this->SetExactParameter(39, "c55_i", this->exact_fct_to_find[39]);
        this->SetExactParameter(40, "c56_i", this->exact_fct_to_find[40]);
        this->SetExactParameter(41, "c66_i", this->exact_fct_to_find[41]);
      }
    else
      {
        cout << "Parameter not known : " << param << endl;
        cout << "TypeInversion = Isotrope, Orthotrope or Anisotrope" << endl;
        abort();
      }
  }

  
  //! linear combination of data
  template<>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension3> >
  ::CombineData(const Matrix<Complex_wp, General, ColMajor>& data0,
                Matrix<Complex_wp, General, ColMajor>& data)
  {
    // default choice : no combination
    data = data0;
  }

  
  //! linear combination of data (and derivative)
  template<>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension3> >
  ::CombineDeriveData(const VectReal_wp& h0, int k,
                      const Matrix<Complex_wp, General, ColMajor>& data_simu,
                      Matrix<Complex_wp, General, ColMajor>& data_comb)
  {
    // default choice : no combination
    data_comb = data_simu;    
  }


  //! linear combination of data (and adjoint of derivative)
  template<>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension3> >
  ::CombineAdjointData(int k, const Matrix<Complex_wp, General, ColMajor>& data_simu,
                       Matrix<Complex_wp, General, ColMajor>& data_comb, VectReal_wp& h0)
  {
    // default choice : no combination
    data_comb = data_simu;    
  }
#endif
  
  
  //! sets the exact parameter c
  template<class Dimension>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension> >::
  SetExactParameter(int k, const string& name, type_ptr_fctC fct)
  {
    BoundaryInverseProblem_Base<Dimension>::SetExactParameter(k, name, fct);
    if (name == "rho")
      this->param_c_type(k) = PARAM_RHO;
    else if (name == "lambda")
      this->param_c_type(k) = PARAM_LAMBDA;
    else if (name == "mu")
      this->param_c_type(k) = PARAM_MU;
    else if (name == "lambda_i")
      this->param_c_type(k) = PARAM_LAMBDA_I;
    else if (name == "mu_i")
      this->param_c_type(k) = PARAM_MU_I;
    else if (name == "c11")
      this->param_c_type(k) = PARAM_C11;
    else if (name == "c12")
      this->param_c_type(k) = PARAM_C12;
    else if (name == "c13")
      this->param_c_type(k) = PARAM_C13;
    else if (name == "c14")
      this->param_c_type(k) = PARAM_C14;
    else if (name == "c15")
      this->param_c_type(k) = PARAM_C15;
    else if (name == "c16")
      this->param_c_type(k) = PARAM_C16;
    else if (name == "c22")
      this->param_c_type(k) = PARAM_C22;
    else if (name == "c23")
      this->param_c_type(k) = PARAM_C23;
    else if (name == "c24")
      this->param_c_type(k) = PARAM_C24;
    else if (name == "c25")
      this->param_c_type(k) = PARAM_C25;
    else if (name == "c26")
      this->param_c_type(k) = PARAM_C26;
    else if (name == "c33")
      this->param_c_type(k) = PARAM_C33;
    else if (name == "c34")
      this->param_c_type(k) = PARAM_C34;
    else if (name == "c35")
      this->param_c_type(k) = PARAM_C35;
    else if (name == "c36")
      this->param_c_type(k) = PARAM_C36;
    else if (name == "c44")
      this->param_c_type(k) = PARAM_C44;
    else if (name == "c45")
      this->param_c_type(k) = PARAM_C45;
    else if (name == "c46")
      this->param_c_type(k) = PARAM_C46;
    else if (name == "c55")
      this->param_c_type(k) = PARAM_C55;
    else if (name == "c56")
      this->param_c_type(k) = PARAM_C56;
    else if (name == "c66")
      this->param_c_type(k) = PARAM_C66;
    else if (name == "c11_i")
      this->param_c_type(k) = PARAM_C11_I;
    else if (name == "c12_i")
      this->param_c_type(k) = PARAM_C12_I;
    else if (name == "c13_i")
      this->param_c_type(k) = PARAM_C13_I;
    else if (name == "c14_i")
      this->param_c_type(k) = PARAM_C14_I;
    else if (name == "c15_i")
      this->param_c_type(k) = PARAM_C15_I;
    else if (name == "c16_i")
      this->param_c_type(k) = PARAM_C16_I;
    else if (name == "c22_i")
      this->param_c_type(k) = PARAM_C22_I;
    else if (name == "c23_i")
      this->param_c_type(k) = PARAM_C23_I;
    else if (name == "c24_i")
      this->param_c_type(k) = PARAM_C24_I;
    else if (name == "c25_i")
      this->param_c_type(k) = PARAM_C25_I;
    else if (name == "c26_i")
      this->param_c_type(k) = PARAM_C26_I;
    else if (name == "c33_i")
      this->param_c_type(k) = PARAM_C33_I;
    else if (name == "c34_i")
      this->param_c_type(k) = PARAM_C34_I;
    else if (name == "c35_i")
      this->param_c_type(k) = PARAM_C35_I;
    else if (name == "c36_i")
      this->param_c_type(k) = PARAM_C36_I;
    else if (name == "c44_i")
      this->param_c_type(k) = PARAM_C44_I;
    else if (name == "c45_i")
      this->param_c_type(k) = PARAM_C45_I;
    else if (name == "c46_i")
      this->param_c_type(k) = PARAM_C46_I;
    else if (name == "c55_i")
      this->param_c_type(k) = PARAM_C55_I;
    else if (name == "c56_i")
      this->param_c_type(k) = PARAM_C56_I;
    else if (name == "c66_i")
      this->param_c_type(k) = PARAM_C66_I;
    else
      {
	cout << "Unknown parameter" << endl;
	abort();
      }
  }

  
  //! Fills i, j, k, l from parameter number num
  template<class Dimension>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension> >
  ::GetIndexParam(int num, int& i, int& j, int& k, int& l, bool& is_imag)
  {
    is_imag = false;
    if (num >= PARAM_C11_I)
      is_imag = true;
    
    switch (num)
      {
      case PARAM_C11: case PARAM_C11_I: i = 0; j = 0; k = 0; l = 0; break;
      case PARAM_C12: case PARAM_C12_I: i = 0; j = 0; k = 1; l = 1; break;
      case PARAM_C13: case PARAM_C13_I: i = 0; j = 0; k = 2; l = 2; break;
      case PARAM_C14: case PARAM_C14_I: i = 0; j = 0; k = 1; l = 2; break;
      case PARAM_C15: case PARAM_C15_I: i = 0; j = 0; k = 0; l = 2; break;
      case PARAM_C16: case PARAM_C16_I: i = 0; j = 0; k = 0; l = 1; break;
      case PARAM_C22: case PARAM_C22_I: i = 1; j = 1; k = 1; l = 1; break;
      case PARAM_C23: case PARAM_C23_I: i = 1; j = 1; k = 2; l = 2; break;
      case PARAM_C24: case PARAM_C24_I: i = 1; j = 1; k = 1; l = 2; break;
      case PARAM_C25: case PARAM_C25_I: i = 1; j = 1; k = 0; l = 2; break;
      case PARAM_C26: case PARAM_C26_I: i = 1; j = 1; k = 0; l = 1; break;
      case PARAM_C33: case PARAM_C33_I: i = 2; j = 2; k = 2; l = 2; break;
      case PARAM_C34: case PARAM_C34_I: i = 2; j = 2; k = 1; l = 2; break;
      case PARAM_C35: case PARAM_C35_I: i = 2; j = 2; k = 0; l = 2; break;
      case PARAM_C36: case PARAM_C36_I: i = 2; j = 2; k = 0; l = 1; break;
      case PARAM_C44: case PARAM_C44_I: i = 1; j = 2; k = 1; l = 2; break;
      case PARAM_C45: case PARAM_C45_I: i = 1; j = 2; k = 0; l = 2; break;
      case PARAM_C46: case PARAM_C46_I: i = 1; j = 2; k = 0; l = 1; break;
      case PARAM_C55: case PARAM_C55_I: i = 0; j = 2; k = 0; l = 2; break;
      case PARAM_C56: case PARAM_C56_I: i = 0; j = 2; k = 0; l = 1; break;
      case PARAM_C66: case PARAM_C66_I: i = 0; j = 1; k = 0; l = 1; break;
      }
  }


  //! returns the variable index
  template<class Dimension>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension> >
  ::GetVariableIndex(int ref,
		     VarHarmonic_Base<Complex_wp, Dimension>& var,
		     Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
		     Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_cplx)
  {
    EllipticProblem<HarmonicElasticEquation<Dimension> >& var_
      = dynamic_cast<EllipticProblem<HarmonicElasticEquation<Dimension> >& >(var);
    
    TinySymmetricTensor<int, Dimension::dim_N> Ci;

    bool ortho = false, aniso = false;
    for (int num_index = 0; num_index < this->param_c_type.GetM(); num_index++)
      {
	if (this->param_c_type(num_index) == PARAM_RHO)
	  {
	    PhysicalVaryingMedia<Dimension, Complex_wp>& rho 
	      = var_.ref_rho(ref).GetVaryingIndex();
	    
	    var_.ref_rho(ref).SetMediaType(rho.USER);
	    rho_cplx.PushBack(&rho);
	  }
	else if ((this->param_c_type(num_index) == PARAM_LAMBDA)
		 || (this->param_c_type(num_index) == PARAM_MU)
		 || (this->param_c_type(num_index) == PARAM_LAMBDA_I)
		 || (this->param_c_type(num_index) == PARAM_MU_I))
	  {
	    if (var_.ref_tensorC(ref).GetAnisotropyType() != var_.ref_tensorC(ref).ISOTROPE)
	      {
		cout << "Material should be isotropic for lambda, mu parameters" << endl;
		DISP(ref); DISP(var_.ref_tensorC(ref).GetAnisotropyType());
		abort();
	      }

	    PhysicalVaryingMedia<Dimension, Complex_wp>& Lambda 
	      = var_.ref_tensorC(ref).GetVaryingIndex(0, 0, 0, 0);
	    
	    var_.ref_tensorC(ref).SetMediaType(Lambda.USER);
	    rho_cplx.PushBack(&Lambda);

	    PhysicalVaryingMedia<Dimension, Complex_wp>& Mu
	      = var_.ref_tensorC(ref).GetVaryingIndex(0, 1, 0, 1);
	    
	    var_.ref_tensorC(ref).SetMediaType(Mu.USER);
	    rho_cplx.PushBack(&Mu);
	  }
	else if (this->param_c_type(num_index) >= PARAM_C11)
	  {
	    if (!aniso)
	      ortho = true;
	    
	    int i(-1), j(-1), k(-1), l(-1); bool is_imag;
	    GetIndexParam(this->param_c_type(num_index), i, j, k, l, is_imag);
	    PhysicalVaryingMedia<Dimension, Complex_wp>& C
	      = var_.ref_tensorC(ref).GetVaryingIndex(i, j, k, l);
	    
	    Ci(i, j, k, l) = 1;
	    
	    var_.ref_tensorC(ref).SetMediaType(C.USER);
	    rho_cplx.PushBack(&C);

	    if ((i != j) || (k != l))
	      if ((i != k) || (j != l))
		{
		  aniso = true; ortho = false;
		}
	  }
      }

    if (ortho)
      if (var_.ref_tensorC(ref).GetAnisotropyType() != var_.ref_tensorC(ref).ORTHOTROPE)
	{
	  cout << "Material should be orthotropic for C_{i, j} parameters" << endl;
	  DISP(ref); DISP(var_.ref_tensorC(ref).GetAnisotropyType());
	  abort();
	}
    
    if (aniso)
      if (var_.ref_tensorC(ref).GetAnisotropyType() != var_.ref_tensorC(ref).ANISOTROPE)
	{
	  cout << "Material should be anisotropic for these parameters" << endl;
	  DISP(ref); DISP(var_.ref_tensorC(ref).GetAnisotropyType());
	  abort();
	}
    
    int nb_comp = Dimension::dim_N*(Dimension::dim_N+1)/2;
    int nb_param = 0;
    for (int i = 0; i < nb_comp; i++)
      for (int j = i; j < nb_comp; j++)
	nb_param += Ci(i, j);
    
    if (Dimension::dim_N == 2)
      {
	if (ortho && (nb_param != 4))
	  {
	    cout << "All components of an orthotropic media must be inverted" << endl;
	    abort();
	  } 
	
	if (aniso && (nb_param != 6))
	  {
	    cout << "All components of an anisotropic media must be inverted" << endl;
	    abort();
	  }
      }

    if (Dimension::dim_N == 3)
      {
	if (ortho && (nb_param != 9))
	  {
	    cout << "All components of an orthotropic media must be inverted" << endl;
	    abort();
	  } 
	
	if (aniso && (nb_param != 21))
	  {
	    cout << "All components of an anisotropic media must be inverted" << endl;
	    abort();
	  }
      }
  }


  //! Fills variables indexes on element i1, point j from values contained in eval_c
  template<class Dimension>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension> >
  ::SetVariableIndexPoint(int ref, int i, int i1, int j, int offset, const Vector<VectReal_wp>& eval_c,
			  Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
			  Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_cplx)
  {
    Complex_wp L = ref_C(ref)(0, 0, 0, 0), Mu = ref_C(ref)(0, 1, 0, 1), C;
    int num = 0, numI = 0;
    bool isotrope = false;
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO:
	    {
	      rho_cplx(num)->SetCoefficient(i1, j, eval_c(p)(offset+j));
	      num++;
	    }
	    break;
	  case PARAM_LAMBDA:
	    L = Complex_wp(eval_c(p)(offset+j), imag(L));
	    isotrope = true; numI = num; num += 2;
	    break;
	  case PARAM_MU:
	    Mu = Complex_wp(eval_c(p)(offset+j), imag(Mu));
	    isotrope = true; numI = num; num += 2;
	    break;
	  case PARAM_LAMBDA_I:
	    L = Complex_wp(real(L), eval_c(p)(offset+j));
	    isotrope = true; numI = num; num += 2;
	    break;
	  case PARAM_MU_I:
	    Mu = Complex_wp(real(Mu), eval_c(p)(offset+j));
	    isotrope = true; numI = num; num += 2;
	    break;
	  case PARAM_C11: case PARAM_C12: case PARAM_C13: case PARAM_C14: case PARAM_C15:
	  case PARAM_C16: case PARAM_C22: case PARAM_C23: case PARAM_C24: case PARAM_C25:
	  case PARAM_C26: case PARAM_C33: case PARAM_C34: case PARAM_C35: case PARAM_C36:
	  case PARAM_C44: case PARAM_C45: case PARAM_C46: case PARAM_C55: case PARAM_C56: case PARAM_C66:
	    C = rho_cplx(num)->GetCoefficient(i1, j);
	    rho_cplx(num)->SetCoefficient(i1, j, Complex_wp(eval_c(p)(offset+j), imag(C)));
	    num++;
	    break;
	  default:
	    C = rho_cplx(num)->GetCoefficient(i1, j);
	    rho_cplx(num)->SetCoefficient(i1, j, Complex_wp(real(C), eval_c(p)(offset+j)));
	    num++;
	    break;
	  }
      }
    
    if (isotrope)
      {
	rho_cplx(numI)->SetCoefficient(i1, j, L);
	rho_cplx(numI+1)->SetCoefficient(i1, j, Mu);
      }
  }
  

  //! local product by derivative with respect to parameters 
  template<class Dimension>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension> >
  ::MltParamDerivLocal(int nk, const Real_wp& omega2, VectComplex_wp& u_glob,
		       int i, const Real_wp& alpha, const Vector<VectReal_wp>& c_quad,
		       int offset, VectComplex_wp& y_col)
  {
    bool stiff = false, mass = false;
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO : mass = true; break;
	  default: stiff = true;
	  }
      }
    
    TinyVector<VectComplex_wp, Dimension::dim_N> gradUh_quad, grad_feval;
    TinyVector<VectComplex_wp, Dimension::dim_N> Uh_dof, Uh_quad, feval;
    const ElementReference<Dimension, 1>& Fb = var_helm.GetReferenceElementH1(i);
    var_helm.GetLocalUnknownVector(u_glob, i, Uh_dof);
    
    // values of u (and/or grad u) on quadrature points are computed
    // Uh_quad contains u(xi_k) and gradUh_quad contains \nabla u(\xi_k)
    int Nquad = Fb.GetNbPointsQuadratureInside();
    for (int k = 0; k < Dimension::dim_N; k++)
      {
	Uh_quad(k).Reallocate(Nquad);
	feval(k).Reallocate(Nquad);
	feval(k).Zero();
	grad_feval(k).Reallocate(Nquad*Dimension::dim_N);
    	grad_feval(k).Zero();
	
	if (stiff)
	  gradUh_quad(k).Reallocate(Nquad*Dimension::dim_N);
	
	if (Fb.UseQuadraturePointsForRh())
	  {
	    Fb.ApplyChTranspose(Uh_dof(k), Uh_quad(k));
	    
	    if (stiff)
	      Fb.ApplyRhQuadratureTranspose(Uh_quad(k), gradUh_quad(k));
	  }
	else
	  {
	    if (mass)
	      Fb.ApplyChTranspose(Uh_dof(k), Uh_quad(k));
	    
	    if (stiff)
	      Fb.ApplyRhTranspose(Uh_dof(k), gradUh_quad(k));
	  }
      }
    
    // we retrieve weight of integration
    VectReal_wp poids(Nquad);
    Vector<typename Dimension::MatrixN_N> dfjm1;
    Vector<TinySymmetricTensor<Complex_wp, Dimension::dim_N> > matC;
    if (stiff)
      {
	dfjm1.Reallocate(Nquad);
	matC.Reallocate(Nquad);
      }
    
    bool affine = var_helm.mesh.IsElementAffine(i);
    if (affine)
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(0)*Fb.WeightsND(k);
	  if (stiff)
	    {
	      dfjm1(k) = var_helm.Glob_DFjm1(i)(0);
	      Mlt(1.0/var_helm.Glob_jacobian(i)(0), dfjm1(k));
	    }
	} 
    else
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(k);
	  if (stiff)
	    {
	      dfjm1(k) = var_helm.Glob_DFjm1(i)(k);
	      Mlt(Fb.WeightsND(k)/var_helm.Glob_jacobian(i)(k), dfjm1(k));
	    }
	}
        
    // we update feval and grad_feval by using the parameters
    int ref = var_helm.mesh.Element(i).GetReference();
    bool ortho = true;
    if (var_helm.ref_tensorC(ref).GetAnisotropyType() == var_helm.ref_tensorC(ref).ANISOTROPE)
      ortho = false;
    
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO:
	    {
	      Real_wp coef = alpha*omega2;
	      for (int k = 0; k < Nquad; k++)
		{
		  Real_wp beta = coef*poids(k)*c_quad(p)(offset+k);
		  for (int s = 0; s < Dimension::dim_N; s++)
		    feval(s)(k) += beta*Uh_quad(s)(k);
		}
	    }
	    break;
	  case PARAM_LAMBDA:
	    {
	      for (int k = 0; k < Nquad; k++)
		{
		  Real_wp beta = -alpha*poids(k)*c_quad(p)(offset+k);
		  matC(k).AddIsotrope(Complex_wp(beta, 0), Complex_wp(0, 0));
		}
	    }
	    break;
	  case PARAM_MU:
	    {
	      for (int k = 0; k < Nquad; k++)
		{
		  Real_wp beta = -alpha*poids(k)*c_quad(p)(offset+k);
		  matC(k).AddIsotrope(Complex_wp(0, 0), Complex_wp(beta, 0));
		}
	    }
	    break;
	  case PARAM_LAMBDA_I:
	    {
	      for (int k = 0; k < Nquad; k++)
		{
		  Real_wp beta = -alpha*poids(k)*c_quad(p)(offset+k);
		  matC(k).AddIsotrope(Complex_wp(0, beta), Complex_wp(0, 0));
		}
	    }
	    break;
	  case PARAM_MU_I:
	    {
	      for (int k = 0; k < Nquad; k++)
		{
		  Real_wp beta = -alpha*poids(k)*c_quad(p)(offset+k);
		  matC(k).AddIsotrope(Complex_wp(0, 0), Complex_wp(0, beta));
		}
	    }
	    break;
	  default:
	    {
	      int i1(-1), j1(-1), k1(-1), l1(-1); bool is_imag(false);
	      GetIndexParam(this->param_c_type(p), i1, j1, k1, l1, is_imag);
	      if (!is_imag)
		for (int k = 0; k < Nquad; k++)
		  matC(k)(i1, j1, k1, l1) -= alpha*poids(k)*c_quad(p)(offset+k);
	      else
		for (int k = 0; k < Nquad; k++)
		  matC(k)(i1, j1, k1, l1) -= Complex_wp(0, alpha*poids(k)*c_quad(p)(offset+k));
	    }
	    break;
	  }
      }
    
    if (stiff)
      {
	TinyVector<Complex_wp, Dimension::dim_N> vec_u, vec_v;
	TinyVector<TinyVector<Complex_wp, Dimension::dim_N>, Dimension::dim_N> grad_u, grad_v;
	for (int k = 0; k < Nquad; k++)
	  {
	    for (int s = 0; s < Dimension::dim_N; s++)
	      {
		CopyVector(gradUh_quad(s), k, vec_u);
		MltTrans(dfjm1(k), vec_u, grad_u(s));
	      }
	    
	    if (ortho)
	      matC(k).MltOrthotrope(grad_u, grad_v);
	    else
	      matC(k).Mlt(grad_u, grad_v);
	    
	    for (int s = 0; s < Dimension::dim_N; s++)
	      {
		Mlt(dfjm1(k), grad_v(s), vec_v);
		CopyVector(vec_v, k, grad_feval(s));
	      }
	  }
      }
    
    // y_col is updated with feval and grad_feval
    for (int k = 0; k < Dimension::dim_N; k++)
      {
	if (Fb.UseQuadraturePointsForRh())
	  {
	    if (stiff)
	      {
		Fb.ApplyRhQuadrature(grad_feval(k), Uh_dof(k));
		feval(k) = feval(k) + Uh_dof(k);
	      }
	    
	    Fb.ApplyCh(feval(k), Uh_dof(k));
	  }
	else
	  {	
	    if (mass)
	      Fb.ApplyCh(feval(k), Uh_dof(k));
	    else
	      Uh_dof(k).Zero();
	    
	    if (stiff)
	      {
		Uh_quad(0).Reallocate(Uh_dof(k).GetM());
		Fb.ApplyRh(grad_feval(k), Uh_quad(0));
		Uh_dof(k) = Uh_dof(k) + Uh_quad(0);
	      }
	  }
      }
    
    var_helm.AddLocalUnknownVector(Complex_wp(1, 0), Uh_dof, i, y_col);
  }
  

  //! local product by transpose of derivative with respect to parameters   
  template<class Dimension>
  void BoundaryInverseProblem<HarmonicElasticEquation<Dimension> >
  ::MltParamDerivTransLoc(int nk, const Real_wp& omega2, VectComplex_wp& u_glob,
			  int i, const Real_wp& alpha, int offset, 
			  const VectComplex_wp& y_col, Vector<VectReal_wp>& ProdC)
  {
    bool stiff = false, mass = false;
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO : mass = true; break;
	  default : stiff = true; break;
	  }
      }
    
    TinyVector<VectComplex_wp, Dimension::dim_N> Uh_dof, Vh_dof, Uh_quad, Vh_quad;
    TinyVector<VectComplex_wp, Dimension::dim_N> gradUh_quad, gradVh_quad;
    const ElementReference<Dimension, 1>& Fb = var_helm.GetReferenceElementH1(i);
    var_helm.GetLocalUnknownVector(u_glob, i, Uh_dof);
    
    var_helm.GetLocalUnknownVector(y_col, i, Vh_dof);
    
    // values of u and v are computed on quadrature points (and gradient if needed)
    int Nquad = Fb.GetNbPointsQuadratureInside();
    for (int k = 0; k < Dimension::dim_N; k++)
      {
	Uh_quad(k).Reallocate(Nquad);
	Vh_quad(k).Reallocate(Nquad);
    
	if (stiff)
	  {
	    gradUh_quad(k).Reallocate(Dimension::dim_N*Nquad);
	    gradVh_quad(k).Reallocate(Dimension::dim_N*Nquad);
	  }
    
	if (Fb.UseQuadraturePointsForRh())
	  {
	    Fb.ApplyChTranspose(Uh_dof(k), Uh_quad(k));
	    Fb.ApplyChTranspose(Vh_dof(k), Vh_quad(k));	
	    if (stiff)
	      {
		Fb.ApplyRhQuadratureTranspose(Uh_quad(k), gradUh_quad(k));
		Fb.ApplyRhQuadratureTranspose(Vh_quad(k), gradVh_quad(k));
	      }
	  }
	else
	  {
	    if (mass)
	      {
		Fb.ApplyChTranspose(Uh_dof(k), Uh_quad(k));
		Fb.ApplyChTranspose(Vh_dof(k), Vh_quad(k));
	      }
	
	    if (stiff)
	      {
		Fb.ApplyRhTranspose(Uh_dof(k), gradUh_quad(k));
		Fb.ApplyRhTranspose(Vh_dof(k), gradVh_quad(k));
	      }
	  }
      }
    
    // we retrieve weight of integration
    VectReal_wp poids(Nquad);
    Vector<typename Dimension::MatrixN_N> dfjm1;
    if (stiff)
      dfjm1.Reallocate(Nquad);

    bool affine = var_helm.mesh.IsElementAffine(i);
    if (affine)
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(0)*Fb.WeightsND(k);
	  if (stiff)
	    {
	      dfjm1(k) = var_helm.Glob_DFjm1(i)(0);
	      Mlt(1.0/var_helm.Glob_jacobian(i)(0), dfjm1(k));
	    }
	} 
    else
      for (int k = 0; k < Nquad; k++)
	{
	  poids(k) = var_helm.Glob_jacobian(i)(k);
	  if (stiff)
	    {
	      dfjm1(k) = var_helm.Glob_DFjm1(i)(k);
	      Mlt(Fb.WeightsND(k)/var_helm.Glob_jacobian(i)(k), dfjm1(k));
	    }
	}
    
    // we update ProdC
    TinyVector<TinyVector<Real_wp, Dimension::dim_N>, Dimension::dim_N>
      vec_ur, vec_ui, vec_vr, vec_vi, vec_ur2, vec_ui2, vec_vr2, vec_vi2;
    //int ref = var_helm.mesh.Element(i).GetReference();
    for (int p = 0; p < this->param_c_type.GetM(); p++)
      {
	switch (this->param_c_type(p))
	  {
	  case PARAM_RHO:
	    {
	      Real_wp coef = alpha*omega2;
	      for (int k = 0; k < Nquad; k++)
		for (int s = 0; s < Dimension::dim_N; s++)
		  ProdC(p)(offset+k) += coef*poids(k)*(real(Uh_quad(s)(k))*real(Vh_quad(s)(k))
						       + imag(Uh_quad(s)(k))*imag(Vh_quad(s)(k)));
	    }
	    break;
	  case PARAM_LAMBDA:
	  case PARAM_LAMBDA_I:
	    {	      
	      for (int k = 0; k < Nquad; k++)
		{
		  Real_wp coef = -alpha*poids(k), div_ur(0), div_ui(0), div_vr(0), div_vi(0);
		  for (int s = 0; s < Dimension::dim_N; s++)
		    {
		      CopyVector(gradUh_quad(s), k, vec_ur(s), vec_ui(s));
		      CopyVector(gradVh_quad(s), k, vec_vr(s), vec_vi(s));
		      MltTrans(dfjm1(k), vec_ur(s), vec_ur2(s));
		      MltTrans(dfjm1(k), vec_ui(s), vec_ui2(s));
		      MltTrans(dfjm1(k), vec_vr(s), vec_vr2(s));
		      MltTrans(dfjm1(k), vec_vi(s), vec_vi2(s));
		      div_ur += vec_ur2(s)(s); div_ui += vec_ui2(s)(s);
		      div_vr += vec_vr2(s)(s); div_vi += vec_vi2(s)(s);
		    }
		  
		  if (this->param_c_type(p) == PARAM_LAMBDA_I)
		    ProdC(p)(offset+k) += coef*(-div_ui*div_vr + div_ur*div_vi);
		  else
		    ProdC(p)(offset+k) += coef*(div_ur*div_vr + div_ui*div_vi);
		}
	    }
	    break;
	  default:
	    {	      
	      for (int k = 0; k < Nquad; k++)
		{
		  Real_wp coef = -alpha*poids(k);
		  for (int s = 0; s < Dimension::dim_N; s++)
		    {
		      CopyVector(gradUh_quad(s), k, vec_ur(s), vec_ui(s));
		      CopyVector(gradVh_quad(s), k, vec_vr(s), vec_vi(s));
		      MltTrans(dfjm1(k), vec_ur(s), vec_ur2(s));
		      MltTrans(dfjm1(k), vec_ui(s), vec_ui2(s));
		      MltTrans(dfjm1(k), vec_vr(s), vec_vr2(s));
		      MltTrans(dfjm1(k), vec_vi(s), vec_vi2(s));
		    }
		  
		  if (this->param_c_type(p) == PARAM_MU_I)
		    for (int s = 0; s < Dimension::dim_N; s++)
		      {
			ProdC(p)(offset+k) += 2*coef*(-vec_ui2(s)(s)*vec_vr2(s)(s) + vec_ur2(s)(s)*vec_vi2(s)(s));
			for (int t = s+1; t < Dimension::dim_N; t++)
			  ProdC(p)(offset+k) += coef*((-vec_ui2(s)(t)-vec_ui2(t)(s))*(vec_vr2(s)(t)+vec_vr2(t)(s))
						      + (vec_ur2(s)(t)+vec_ur2(t)(s))*(vec_vi2(s)(t)+vec_vi2(t)(s)));
		      }
		  else if (this->param_c_type(p) == PARAM_MU)
		    for (int s = 0; s < Dimension::dim_N; s++)
		      {
			ProdC(p)(offset+k) += 2*coef*(vec_ur2(s)(s)*vec_vr2(s)(s) + vec_ui2(s)(s)*vec_vi2(s)(s));
			for (int t = s+1; t < Dimension::dim_N; t++)
			  ProdC(p)(offset+k) += coef*((vec_ur2(s)(t)+vec_ur2(t)(s))*(vec_vr2(s)(t)+vec_vr2(t)(s))
						      + (vec_ui2(s)(t)+vec_ui2(t)(s))*(vec_vi2(s)(t)+vec_vi2(t)(s)));
		      }
		  else
		    {
		      // symmetrization of grad u
		      vec_ur2(0)(1) += vec_ur2(1)(0); vec_ur2(1)(0) = vec_ur2(0)(1);
		      vec_ui2(0)(1) += vec_ui2(1)(0); vec_ui2(1)(0) = vec_ui2(0)(1);
		      vec_vr2(0)(1) += vec_vr2(1)(0); vec_vr2(1)(0) = vec_vr2(0)(1);
		      vec_vi2(0)(1) += vec_vi2(1)(0); vec_vi2(1)(0) = vec_vi2(0)(1);
		      if (Dimension::dim_N == 3)
			{
			  vec_ur2(0)(2) += vec_ur2(2)(0); vec_ur2(2)(0) = vec_ur2(0)(2);
			  vec_ur2(1)(2) += vec_ur2(2)(1); vec_ur2(2)(1) = vec_ur2(1)(2);
			  vec_ui2(0)(2) += vec_ui2(2)(0); vec_ui2(2)(0) = vec_ui2(0)(2);
			  vec_ui2(1)(2) += vec_ui2(2)(1); vec_ui2(2)(1) = vec_ui2(1)(2);
			  vec_vr2(0)(2) += vec_vr2(2)(0); vec_vr2(2)(0) = vec_vr2(0)(2);
			  vec_vr2(1)(2) += vec_vr2(2)(1); vec_vr2(2)(1) = vec_vr2(1)(2);
			  vec_vi2(0)(2) += vec_vi2(2)(0); vec_vi2(2)(0) = vec_vi2(0)(2);
			  vec_vi2(1)(2) += vec_vi2(2)(1); vec_vi2(2)(1) = vec_vi2(1)(2);
			}
		      
		      int i1(-1), j1(-1), k1(-1), l1(-1); bool is_imag(false);
		      GetIndexParam(this->param_c_type(p), i1, j1, k1, l1, is_imag);
		      if ((i1 == k1) && (j1 == l1))
			{
			  if (is_imag)
			    ProdC(p)(offset+k) += coef*(-vec_ui2(i1)(j1)*vec_vr2(k1)(l1) + vec_ur2(i1)(j1)*vec_vi2(k1)(l1));
			  else
			    ProdC(p)(offset+k) += coef*(vec_ur2(i1)(j1)*vec_vr2(k1)(l1) + vec_ui2(i1)(j1)*vec_vi2(k1)(l1));
			}
		      else
			{
			  if (is_imag)
			    ProdC(p)(offset+k) += coef*(-vec_ui2(i1)(j1)*vec_vr2(k1)(l1)
							- vec_ui2(k1)(l1)*vec_vr2(i1)(j1)
							+ vec_ur2(i1)(j1)*vec_vi2(k1)(l1)
							+ vec_ur2(k1)(l1)*vec_vi2(i1)(j1));
			  else
			    ProdC(p)(offset+k) += coef*(vec_ur2(i1)(j1)*vec_vr2(k1)(l1)
							+ vec_ur2(k1)(l1)*vec_vr2(i1)(j1)
							+ vec_ui2(i1)(j1)*vec_vi2(k1)(l1)
							+ vec_ui2(k1)(l1)*vec_vi2(i1)(j1));
			}
		    }
		}
	    }
	    break;
	  }
      }
  }

}

#define MONTJOIE_FILE_BOUNDARY_INVERSE_PROBLEM_CXX
#endif
