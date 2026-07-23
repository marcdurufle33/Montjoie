#ifndef MONTJOIE_FILE_HELMHOLTZ_POLAR_CXX

#include "HelmholtzPolarInline.cxx"

namespace Montjoie
{
  
  /**********************
   * VarHelmholtz_Polar *
   **********************/
  

  //! default constructor
  template<class Complexe>
  VarHelmholtz_Polar<Complexe>::VarHelmholtz_Polar()
  {
    Lmax = 0;
    number_mode_to_be_computed = false; 
    this->type_coordinate = this->POLAR;
    type_source = SRC_DIFFRACTED_FIELD;
    wave_vector.Init(0, 1);
    radius_gaussian_src = 0.0;
    radius_gaussian_cut_off = 0.0;
    display_du_dr = false;
  }  
  
  
  //! reads a line of the data file
  template<class Complexe>
  void VarHelmholtz_Polar<Complexe>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (description_field == "NumberModes")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHelmholtz_Polar" << endl;
	    cout << "NumberModes needs more parameters, for instance :" << endl;
	    cout << "NumberModes = n" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        if (parameters(0) == "AUTO")
          {
            if (parameters.GetM() <= 2)
              {
                cout << "In SetInputData of VarHelmholtz_Polar" << endl;
                cout << "NumberModes needs more parameters, for instance :" << endl;
                cout << "NumberModes = AUTO n threshold" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }
            
            number_mode_to_be_computed = true;
            Lmax = to_num<int>(parameters(1));
            threshold_mode = to_num<Real_wp>(parameters(2));
          }
        else
          {
            number_mode_to_be_computed = false;
            Lmax = to_num<int>(parameters(0));
          }
      }
    else if (description_field == "TypeSource")
      {
        if (parameters(0) == "SRC_DIFFRACTED_FIELD")
          type_source = SRC_DIFFRACTED_FIELD;
        else if (parameters(0) == "SRC_TOTAL_FIELD")
          type_source = SRC_TOTAL_FIELD;
        else if (parameters(0) == "SRC_DIRAC")
          {
            if (parameters.GetM() <= 2)
              {
                cout << "In SetInputData of VarHelmholtz_Polar" << endl;
                cout << "TypeSource needs more parameters, for instance :" << endl;
                cout << "TypeSource = SRC_DIRAC x0 y0" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();                
              }
            
            type_source = SRC_DIRAC;
            Real_wp x0 = to_num<Real_wp>(parameters(1));
            Real_wp y0 = to_num<Real_wp>(parameters(2));
            origin_dirac.Init(x0, y0);
          }
        else if ((parameters(0) == "SRC_GAUSSIAN") || (parameters(0) == "SRC_GAUSSIAN_R"))
          {
            if (parameters.GetM() <= 4)
              {
                cout << "In SetInputData of VarHelmholtz_Polar" << endl;
                cout << "TypeSource needs more parameters, for instance :" << endl;
                cout << "TypeSource = SRC_GAUSSIAN x0 y0 z0 r0" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();                
              }
            
            if (parameters(0) == "SRC_GAUSSIAN")
	      type_source = SRC_GAUSSIAN;
	    else
	      type_source = SRC_GAUSSIAN_R;
	    
            Real_wp x0 = to_num<Real_wp>(parameters(1));
            Real_wp y0 = to_num<Real_wp>(parameters(2));
            origin_dirac.Init(x0, y0);
            radius_gaussian_src = to_num<Real_wp>(parameters(3));
            if (parameters.GetM() > 4)
              {
                radius_gaussian_cut_off = to_num<Real_wp>(parameters(4));
              }
            else
              radius_gaussian_cut_off = radius_gaussian_src;
          }
      }
    else if (!description_field.compare("IncidentAngle"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarHelmholtz_Polar" << endl;
	    cout << "IncidentAngle needs more parameters, for instance :" << endl;
	    cout << "IncidentAngle = teta" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	// in 2-D, the wave vector is oriented along
	// (cos(teta), sin(teta))
	Real_wp teta;
        teta = to_num<Real_wp>(parameters(0))*pi_wp/180;
	wave_vector(0) = cos(teta); wave_vector(1) = sin(teta);    
      }
    else if (!description_field.compare("ElectricOrMagnetic"))
      {
	if (parameters(0) == "-1")
	  display_du_dr = true;
	else
	  display_du_dr = false;
      }
    else if ((description_field == "SismoPlane")
             || (description_field == "SismoPointsFile")
             || (description_field == "SismoLine")
             || (description_field == "SismoPoint")
             || (description_field == "SismoCircle"))
      {
	// 2-D outputs
	int i = var_grid.GetM();
	var_grid.Resize(var_grid.GetM()+1);
	var_grid(i).SetInputData(description_field, parameters);
      }
    else if ((!description_field.compare("FileOutputPlane")) ||
             (!description_field.compare("FileOutputLine")) 
	     || (!description_field.compare("FileOutputPoint"))  || 
             (!description_field.compare("FileOutputCircle")) ||
             (!description_field.compare("FileOutputPointsFile")))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarHelmholtzPolar" << endl;
	    cout << "FileOutputXXX needs more parameters, for instance :" << endl;
	    cout << "FileOutputPlane = file_diffrac file_total component coef" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	int type = GridInterpolationFull<Dimension1>::GetType(description_field);
	int nb = GetIndexOutputFiles(type);
	
	// Files of storage of solution on a regular grid. Diffracted and total field
	string ext = ".dat";
	
	this->output_grid_param(nb).SetDiffractedFieldFile(GetBaseString(parameters(0)) + ext);
	this->output_grid_param(nb).SetTotalFieldFile(GetBaseString(parameters(1)) + ext);
      }
    else
      VarHelmholtz_1D<Complexe>::SetInputData(description_field, parameters);
  }


  //! inits the computation before solving each mode
  template<class Complexe>
  void VarHelmholtz_Polar<Complexe>
  ::ConstructAll(const string& input_file, const string& name_elt, MontjoieTimer& var_chrono)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(this->comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    var_chrono.SetMessage("InputFile", "read the input file");
    var_chrono.Start("InputFile");

    ReadInputFile(input_file, *this);
    if ((this->print_level >= 0) && (rank_proc == 0))
      cout << "Input file has been successfully read" << endl;
    
    var_chrono.Stop("InputFile");
    
    var_chrono.SetMessage("ConstructMesh", "construct the mesh");
    var_chrono.Start("ConstructMesh");

    // computation of finite element and mesh
    this->ComputeMeshAndFiniteElement(name_elt);
    
    var_chrono.Stop("ConstructMesh");

    var_chrono.SetMessage("ConstructMass", "construct the mass matrix");
    var_chrono.Start("ConstructMass");

    this->ComputePhysicalCoefficients();

    // retrieving rho0 and mu0
    Complexe v0;
    int ne = this->mesh.GetNbElt()-1;
    int ref = this->mesh.Element(ne).GetReference();    
    int order = this->GetReferenceElement(ne).GetNbPointsQuadratureInside();
    this->GetPhysicalCoefficient(this->rho0, this->mu0, this->sigma0, v0,
				 ne, order+1, ref);        

    Complexe rho = this->rho0 + Iwp*this->sigma0/this->omega;
    Complexe kinf = this->omega*sqrt(rho/this->mu0);
    
    if (this->condition_right_side == BoundaryConditionEnum::LINE_TRANSPARENT)
      {
        // for a transparent condition, we precompute ratios h_n' / h_n
	Real_wp kr = realpart(kinf*this->xmax_);
        VectComplex_wp Jn, Hn, dJn, dHn;
	ComputeDeriveBesselAndHankel(Real_wp(0), Lmax+2, kr, Jn, Hn, dJn, dHn);
	
        // then forming the coefficient -mu r h_n' / h_n
        coef_impedance_transparent.Reallocate(Lmax+1);
        Complex_wp coef = this->mu0*this->xmax_*kinf;
        for (int l = 0; l <= Lmax; l++)
          coef_impedance_transparent(l) = -coef*dHn(l)/Hn(l);
      }
    
    // computing geometry quantities
    this->ComputeMassMatrix();
    
    int N = this->GetNbDof();
    int r = this->mesh.GetOrder();
    if ((this->print_level >= 0) && (rank_proc == 0))
      {
        cout << "Number of degrees of freedom = " << N << endl;
        cout << "Mesh size : " << this->mesh.GetMeshSize() << endl;
      }

    var_chrono.Stop("ConstructMass");

    var_chrono.SetMessage("ComputeFEM", "compute finite element matrices");
    var_chrono.Start("ComputeFEM");
    
    // computation of finite element matrices
    GlobalGenericMatrix<Complexe> nat_mat;

    if (this->mixed_formulation)
      {
	stiffness_matrix_sp.Reallocate(N, N);
	this->AddMatrixFEM(stiffness_matrix_sp, nat_mat);
	this->AddBoundaryTerms(stiffness_matrix_sp, nat_mat);

	stiffness_matrix_sp.WriteText("stiffness.dat");
      }
    else
      {
	stiffness_matrix.Reallocate(N, N, r, r);
	this->AddMatrixFEM(stiffness_matrix, nat_mat);
	this->AddBoundaryTerms(stiffness_matrix, nat_mat);

	if (this->print_level >= 4)
	  stiffness_matrix.WriteText("stiffness.dat");
      }
    
    // checking if mass lumped elements are used
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      if (!this->GetReferenceElement(i).LumpedMassMatrix())
        {
          cout << "In HelmholtzPolar, the mass matrix is "
               << "assumed to be diagonal" << endl;
          
          cout << "Use EDGE_LOBATTO or EDGE_RADAU element to achieve mass lumping" << endl;
          
          abort();
        }

    // diagonal part due to l^2 part
    mass_matrix.Reallocate(N);
    mass_matrix.Fill(0);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      if (this->mat_mass_mu.GetM() > 0)
	for (int j = 0; j < this->mat_mass_mu(i).GetM(); j++)
	  mass_matrix(this->mesh.GetNumberDof(i, j)) += this->mat_mass_mu(i)(j);

    SetComplexZero(mass_matrix(0));
    
    if (this->condition_right_side == BoundaryConditionEnum::LINE_DIRICHLET)
      SetComplexZero(mass_matrix(N-1));

    // intermediary arrays are cleared to save memory
    this->mat_mass_mu.Clear(); this->mat_mass.Clear();
    this->mat_massV.Clear(); this->mat_sigma.Clear(); this->mat_sigmaV.Clear();
    
    var_chrono.Stop("ComputeFEM");
  }
  
  
  //! computes the display grid
  template<class Complexe>
  void VarHelmholtz_Polar<Complexe>::ComputeVarGrid()
  {
    // boundaries of the grid are set to boundaries of the domain if required
    for (int i = 0; i < var_grid.GetM(); i++)
      if (var_grid(i).GridIntervalToBeComputed())
	{
	  var_grid(i).SetXmin0(this->xmin_);
	  var_grid(i).SetXmax0(this->xmax_);
        }
    
    // components of unknowns to store in output files 
    // are selected with the global variable choice_field_output
    // this->choice_field_output
    for (int i = 0; i < this->output_grid_param.GetM(); i++)
      output_grid_param(i).SetComponent(0);
    
    for (int i = 0; i < var_grid.GetM(); i++)
      var_grid(i).InitGrid(this->var_section, 2);

    // localization of the points of the predefined grids inside the mesh
    // implementation of this method in file GridInterpolation.cxx
    this->var_section.LocalizePoints(this->mesh);    
  }
  
  
  //! computes and factorizes finite element matrix for mode l
  template<class Complexe>
  void VarHelmholtz_Polar<Complexe>
  ::PerformFactorizationStep(Matrix<Complex_wp, General, BandedCol>& mat_lu,
			     All_MatrixLU<Complex_wp>& mat_lu_sp,
                             IVect& pivot, int l)
  {
    Real_wp coef = l*l;
    int N = this->mesh.GetNbDof();
    if (this->mixed_formulation)
      {
	Matrix<Complex_wp, General, ArrayRowSparse> mat_sp;
	mat_sp = stiffness_matrix_sp;
	if (l != 0)
	  for (int i = 0; i < this->mesh.GetNbDof(); i++)
	    mat_sp.AddInteraction(i, i, mass_matrix(i)*coef*Iwp/this->omega);
	
	// part due to transparent condition
	if (this->condition_right_side == BoundaryConditionEnum::LINE_TRANSPARENT)
	  mat_sp.AddInteraction(N-1, N-1, coef_impedance_transparent(l)*Iwp/this->omega);

	mat_sp.WriteText("mat.dat");
	mat_lu_sp.Factorize(mat_sp);

	int ierr = 0;
	int type = mat_lu_sp.GetInfoFactorization(ierr);
        PrintFactorizationFailed(type, ierr);	
      }
    else
      {
	mat_lu = stiffness_matrix;
	//mat_lu.WriteText("MatBefore.dat");
	
	for (int i = 0; i < this->mesh.GetNbDof(); i++)
	  mat_lu.AddInteraction(i, i, mass_matrix(i)*coef);

	// part due to transparent condition
	if (this->condition_right_side == BoundaryConditionEnum::LINE_TRANSPARENT)
	  mat_lu.AddInteraction(N-1, N-1, coef_impedance_transparent(l));
	
	// try with Dirichlet
	IVect last_row(1);
	last_row(0) = N-1;
	if (this->condition_right_side == BoundaryConditionEnum::LINE_DIRICHLET)
	  {
	    EraseCol(last_row, mat_lu);
	    mat_lu.Set(N-1, N-1, Complex_wp(1, 0));
	  }
	
	if (this->print_level >= 8)
	  mat_lu.WriteText("MatAfter.dat");
	
	GetLU(mat_lu, pivot);	
      }
    

  }
  
  
  //! solves a linear system with mat_lu
  template<class Complexe>
  void VarHelmholtz_Polar<Complexe>
  ::ComputeSolution(Matrix<Complex_wp, General, BandedCol>& mat_lu,
		    All_MatrixLU<Complex_wp>& mat_lu_sp,
                    IVect& pivot, Vector<Complex_wp>& x)
  {  
    if (this->mixed_formulation)
      mat_lu_sp.Solve(x);
    else
      SolveLU(mat_lu, pivot, x);
  }
  
  
  //! computes right hand side for all modes l
  template<class Complexe>
  void VarHelmholtz_Polar<Complexe>
  ::ComputeRightHandSide(Vector<Vector<Complexe> >& rhs, MontjoieTimer& var_chrono)
  {
    rhs.Reallocate(2*Lmax+1);
    rhs(0).Reallocate(this->GetNbDof());
    rhs(0).Zero();
    
    const ElementReference<Dimension1, 1>& Fb = this->GetReferenceElement(0);
    VectReal_wp phi;
    Fb.ComputeValuesPhiRef(Real_wp(0), phi);
    Real_wp coef = Real_wp(1) / (2.0*pi_wp);
    for (int i = 0; i < phi.GetM(); i++)
      rhs(0)(this->mesh.GetNumberDof(0, i)) = coef*phi(i);
  }


  //! write solution on required outputs
  template<class Complexe>
  void VarHelmholtz_Polar<Complexe>
  ::WriteDatas(Vector<Vector<Complexe> >& sol)
  {
    int L = 0;
    for (int l = 1; l < sol.GetM(); l++)
      if (sol(l).GetM() > 0)
        L = (l+1)/2;

    Vector<Complexe> u(2*L+1);
    VectReal_wp val_phi;
    Complex_wp phase, phase_elt;
    Vector<Vector<Real_wp> > P;
    for (int n = 0; n < var_grid.GetM(); n++)
      {
        const IVect& list_points = var_grid(n).GetPointNumber();
	Vector<Complexe> trace_u(list_points.GetM());
        for (int p = 0; p < list_points.GetM(); p++)
          {
            int num_point = list_points(p);
            int ne = this->var_section.GetElementNumber(num_point);
            Real_wp pt_loc = this->var_section.GetLocalCoordinate(num_point);
	    Complex_wp val_u(0, 0);
            if (ne >= 0)
              {
                u.Fill(0);
                const ElementReference<Dimension1, 1>& Fb = this->GetReferenceElement(ne);
                int nb_dof_loc = Fb.GetNbDof();
                val_phi.Reallocate(nb_dof_loc);
		if (display_du_dr)
		  {
		    TinyVector<Real_wp, 2> s;
		    Fb.ComputeGradientPhiRef(pt_loc, val_phi);
		    this->mesh.GetVerticesElement(ne, s);
		    Real_wp h = abs(s(1) - s(0));
		    Mlt(1.0/h, val_phi);
		  }
		else	  
		  Fb.ComputeValuesPhiRef(pt_loc, val_phi);
                
                for (int l = 0; l <= 2*L; l++)
		  if (sol(l).GetM() > 0)
		    for (int j = 0; j < nb_dof_loc; j++)
		      u(l) += sol(l)(this->mesh.GetNumberDof(ne, j)) * val_phi(j);
                
		Real_wp teta = this->var_section.TetaInterp(num_point);
                phase_elt = Complex_wp(cos(teta), sin(teta));
		phase = phase_elt;
		if (sol(0).GetM() > 0)
		  val_u += u(0);

                for (int l = 0; l < L; l++)
                  {
		    if (sol(2*l+1).GetM() > 0)
		      val_u += phase*u(2*l+1) + conjugate(phase)*u(2*l+2);
		    
		    phase *= phase_elt;
                  }     
              }
	    
	    trace_u(p) = val_u;
          }
        
        int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
        bool ascii = false;

        string file_name_diff = output_grid_param(n).GetDiffractedFieldFile();
        string file_name_tot = output_grid_param(n).GetTotalFieldFile();
                
        if ((type_source == SRC_DIFFRACTED_FIELD) ||
            (type_source == SRC_TOTAL_FIELD))
          {
            if (type_source == SRC_DIFFRACTED_FIELD)
              WriteMatlab(trace_u, var_grid(n), file_name_diff, double_prec, ascii);
            else
              WriteMatlab(trace_u, var_grid(n), file_name_tot, double_prec, ascii);
            
            // computing the total field or diffracted field by computing the incident field
            Real_wp kinf = this->omega*sqrt(abs(this->rho0/this->mu0));
            Real_wp kx = wave_vector(0), ky = wave_vector(1);
            Complex_wp u_inc;
            for (int p = 0; p < list_points.GetM(); p++)
              {
                int num_point = list_points(p);
                int ne = this->var_section.GetElementNumber(num_point);
                if (ne >= 0)
                  {
                    Real_wp r = this->var_section.GetGlobalCoordinate(num_point);
                    Real_wp teta = this->var_section.TetaInterp(num_point);
                    Real_wp x = r*cos(teta), y = r*sin(teta);
                    u_inc = Complex_wp(cos(kinf*(kx*x + ky*y)), sin(kinf*(kx*x + ky*y)));
                    if (type_source == SRC_DIFFRACTED_FIELD)
                      trace_u(p) += u_inc;
                    else
                      trace_u(p) -= u_inc;
                  }
              }
	    
            if (type_source == SRC_DIFFRACTED_FIELD)
              WriteMatlab(trace_u, var_grid(n), file_name_tot, double_prec, ascii);
            else
              WriteMatlab(trace_u, var_grid(n), file_name_diff, double_prec, ascii);
          }
        else
          WriteMatlab(trace_u, var_grid(n), file_name_tot, double_prec, ascii);
      }
  }
  
}

#define MONTJOIE_FILE_HELMHOLTZ_POLAR_CXX
#endif
