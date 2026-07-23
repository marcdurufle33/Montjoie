#ifndef MONTJOIE_FILE_OUTPUT_HARMONIC_CXX

namespace Montjoie
{
  
  //! default constructor
  void VarOutputProblem_Base::InitDefaultValues()
  {
    grid_to_be_computed = true;
    keep_grid_localization_array = false;
    display_solution = true;
    choice_field_output = -1;

    format_output_file = OutputTypeEnum::BINARY_MATLAB;
    precision_output_file = OutputTypeEnum::DOUBLE_PRECISION;

    output_lobatto_points = false;
    output_file_extension = ".dat";
    output_solution_file = false;

    this->move_points_on_output_mesh = false;
    this->component_used_to_translate = 0;
    this->coefficient_used_to_translate = 0.0;

    write_solution_on_quad_points = false;
    last_output_solution = false;
    DOSSIER_output = "./";
  }


  //! gets index for a type of output
  /*!
    \param[in] type type of output (sismo on points, lines, etc)
    \return the index array where the type is stored
  */
  int VarOutputProblem_Base::GetIndexOutputFiles(int type)
  {
    int index = -1;
    for (int i  = 0; i < this->output_grid_param.GetM(); i++)
      if (this->output_grid_param(i).GetNature() == type)
	index = i;
    
    if (index == -1)
      {
	int N = this->output_grid_param.GetM();
	this->output_grid_param.Resize(N+1);
	this->output_grid_param(N).SetNature(type);
	index = N;
      }
    
    return index;
  }
  
  
  //! writing of the solution on all types of outputs (plane, lines, points ...)
  template<class T>
  void VarOutputProblem_Base::WriteDatas(Vector<T> & U0)
  {
    Vector<Vector<T> > Uvec;

    Uvec.SetData(1, &U0);
    WriteDatas(Uvec);
    
    Uvec.Nullify();
  }


  //! writing of the solution on all types of outputs (plane, lines, points ...)
  void VarOutputProblem_Base::WriteDatas(Vector<VectReal_wp> & U0)
  {
    if (!grid_to_be_computed)
      {
	this->InitVarGrid();
	this->ComputeVarGrid();
	grid_to_be_computed = true;
      }
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    if ((rank_proc == 0)||(var_problem.print_level >= 4))
      if (var_problem.print_level >= 0)
        DISP(Norm2(U0(0)));
    
    this->WriteOutputFile(U0, 1);
  }


  //! writing of the solution on all types of outputs (plane, lines, points ...)
  void VarOutputProblem_Base::WriteDatas(Vector<VectComplex_wp> & U0)
  {
    if (!grid_to_be_computed)
      {
	this->InitVarGrid();
	this->ComputeVarGrid();
	grid_to_be_computed = true;
      }

#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    if ((rank_proc == 0)||(var_problem.print_level >= 2))
      if (var_problem.print_level >= 0)
        DISP(Norm2(U0(0)));

    Complex_wp one; SetComplexOne(one);
    if (var_source.GetSourceType(0) == var_source.SRC_DIFFRACTED_FIELD)
      {
	// diffracted field is written
	this->WriteOutputFile(U0, 0);
	
	var_source.AddIncidentWave(one, U0);
	
	// total field is written
        this->WriteOutputFile(U0, 1);
      }
    else
      {
	// total field is written
	this->WriteOutputFile(U0, 1);
	
	
	if (var_source.GetSourceType(0) == var_source.SRC_TOTAL_FIELD)
	  {
	    var_source.AddIncidentWave(-one, U0);
	    
	    // diffracted field is written
	    this->WriteOutputFile(U0, 0);
	  }
      }
  }
  
  
  //! internal function
  template<class T>
  void VarOutputProblem_Base
  ::GetModalOutput(int n, int nm, FftInterface<Complex_wp>& fft_interface,
		   Vector<T>& val_v, Vector<T>& val_grad_v, int off_u,                   
		   Vector<T>& val_u, int off_v, Vector<T>& val_grad_u, bool compute_grad) const
  {
    int nb_comp = var_problem.GetNbComponentsUnknown(nm);
    int nb_comp_grad = var_problem.GetNbComponentsGradient(nm);
    
    if (var_boundary.GetSymmetryType() == var_boundary.NO_SYMMETRY)
      {        
	// no modes are present, we take the first components
        for (int i = 0; i < nb_comp; i++)
          val_u(off_u + i) = val_v(i);
        
	if (compute_grad)
	  for (int i = 0; i < nb_comp_grad; i++)
            val_grad_u(off_v + i) = val_grad_v(i);
	
        return;
      }
    
    int ix, iy, iz; 
    bool sym_teta;
    var_boundary.GetPeriodicModes(n, ix, iy, iz, sym_teta);

    int nx, ny, nz; 
    var_boundary.GetPeriodicModes(var_boundary.GetCurrentModeNumber(), nx, ny, nz, sym_teta);
    
    if (iz < 0)
      {
        if (iy < 0)
          {
            if (ix < 0)
              {
                // no symmetry, we copy the first components
                for (int i = 0; i < nb_comp; i++)
                  val_u(off_u + i) = val_v(i);
                
		if (compute_grad)
                  for (int i = 0; i < nb_comp_grad; i++)
                    val_grad_u(off_v + i) = val_grad_v(i);
              }
            else
              {
                if (var_boundary.ModesNotStored())
                  {
                    T coef;
		    to_complex(fft_interface.GetCoefficient(ix, nx), coef);
                    for (int i = 0; i < nb_comp; i++)
                      val_u(off_u + i) = coef*val_v(i);

		    if (compute_grad)
                      for (int i = 0; i < nb_comp_grad; i++)
                        val_grad_u(off_v + i) = coef*val_grad_v(i);
                  }
                else
                  {
                    Vector<T> vec(nx);
                    for (int i = 0; i < nb_comp; i++)
                      {
                        for (int j = 0; j < nx; j++)
                          vec(j) = val_v(nb_comp*j + i);
                        
                        // 1-D discrete Fourier transform
                        this->ApplyInversePoint(fft_interface, ix, vec, val_u(off_u + i));
                      }
                    
                    if (compute_grad)
                      for (int i = 0; i < nb_comp_grad; i++)
                        {
                          for (int j = 0; j < nx; j++)
                            vec(j) = val_grad_v(nb_comp_grad*j + i);
                        
                          // 1-D discrete Fourier transform
                          this->ApplyInversePoint(fft_interface, ix, vec, val_grad_u(off_v + i));
                      }
                  }
              }
          }
        else
          {
            if (var_boundary.ModesNotStored())
              {
                T coef;
		to_complex(fft_interface.GetCoefficient(ix, iy, nx, ny), coef);
                for (int i = 0; i < nb_comp; i++)
                  val_u(off_u + i) = coef*val_v(i);
                
                if (compute_grad)
                  for (int i = 0; i < nb_comp_grad; i++)
                    val_grad_u(off_v + i) = coef*val_grad_v(i);
              }
            else
              {
                Vector<T> vec(nx*ny);
                for (int i = 0; i < nb_comp; i++)
                  {
                    for (int j = 0; j < nx*ny; j++)
                      vec(j) = val_v(nb_comp*j + i);
                    
                    // 2-D discrete Fourier transform
                    this->ApplyInversePoint(fft_interface, ix, iy, vec, val_u(off_u + i));
                  }
                
                if (compute_grad)
                  for (int i = 0; i < nb_comp_grad; i++)
                    {
                      for (int j = 0; j < nx*ny; j++)
                        vec(j) = val_grad_v(nb_comp_grad*j + i);
                      
                      // 2-D discrete Fourier transform
                      this->ApplyInversePoint(fft_interface, ix, iy, vec, val_grad_u(off_v + i));
                    }
              }
          }
      }
    else
      {
        if (var_boundary.ModesNotStored())
          {
            T coef;
	    to_complex(fft_interface.GetCoefficient(ix, iy, iz, nx, ny, nz), coef);
            for (int i = 0; i < nb_comp; i++)
              val_u(off_u + i) = coef*val_v(i);
            
            if (compute_grad)
              for (int i = 0; i < nb_comp_grad; i++)
                val_grad_u(off_v + i) = coef*val_grad_v(i);
          }
        else
          {
            Vector<T> vec(nx*ny*nz);
            for (int i = 0; i < nb_comp; i++)
              {
                for (int j = 0; j < nx*ny*nz; j++)
                  vec(j) = val_v(nb_comp*j + i);
                
                // 3-D discrete Fourier transform
                this->ApplyInversePoint(fft_interface, ix, iy, iz, vec, val_u(off_u + i));
              }
            
            if (compute_grad)
              for (int i = 0; i < nb_comp_grad; i++)
                {
                  for (int j = 0; j < nx*ny*nz; j++)
                    vec(j) = val_grad_v(nb_comp_grad*j + i);
                  
                  // 3-D discrete Fourier transform
                  this->ApplyInversePoint(fft_interface, ix, iy, iz, vec, val_grad_u(off_v + i));
                }
          }
      }
    
    if (sym_teta)
      {
        Real_wp cos_nalpha(1), sin_nalpha(0);
        fft_interface.GetCosSinAlpha(ix, cos_nalpha, sin_nalpha);
        
        if (nb_comp > 1)
          {
            T vx = cos_nalpha * val_u(off_u) - sin_nalpha * val_u(off_u+1);
            T vy = sin_nalpha * val_u(off_u) + cos_nalpha * val_u(off_u+1);
            val_u(off_u) = vx; val_u(off_u+1) = vy;
          }

        if ((nb_comp_grad > 1) && compute_grad)
          {
            T vx = cos_nalpha * val_grad_u(off_v) - sin_nalpha * val_grad_u(off_v+1);
            T vy = sin_nalpha * val_grad_u(off_v) + cos_nalpha * val_grad_u(off_v+1);
            val_grad_u(off_v) = vx; val_grad_u(off_v+1) = vy;
          }
      }
  }
  

  //! reading of parameters of the data file associated with outputs
  void VarOutputProblem_Base
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("OutputFormat"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarOutputProblem" << endl;
	    cout << "OutputFormat needs more parameters, for instance :" << endl;
	    cout << "OutputFormat = ASCII FLOAT" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	// OutputFormat = Format Precision
	// where Format can be BINARY, ASCII, VTK or ASCII_VTK
	// BINARY : output files are written in Matlab format (.dat) 
        //          and Medit format (.mesh / .bb) for outputs on meshes 
	// ASCII : the Matlab files are written in ascii, however the function
	//         loadND.m reads only binary files, so this format should be avoided
	// VTK : output files are written in binary vtk format (both for meshes and grids)
	//       they can be visualized with paraview
	// ASCII_VTK : output files are written in ascii vtk format
	//            it is better to use VTK since the files are smaller
	
	// Precision can be FLOAT (single precision) or DOUBLE (double precision)
	format_output_file = OutputTypeEnum::BINARY_MATLAB;
        output_file_extension = ".dat";
	if (!parameters(0).compare("ASCII"))
	  format_output_file = OutputTypeEnum::ASCII_MATLAB;
        else if (!parameters(0).compare("VTK"))
          {
            format_output_file = OutputTypeEnum::BINARY_VTK;
            output_file_extension = ".vtk";
          }
        else if (!parameters(0).compare("ASCII_VTK"))
          {
            format_output_file = OutputTypeEnum::ASCII_VTK;
            output_file_extension = ".vtk";
          }
	
	precision_output_file = OutputTypeEnum::SINGLE_PRECISION;
	if (parameters(1) == "DOUBLE")
	  precision_output_file = OutputTypeEnum::DOUBLE_PRECISION;
        else if (parameters(1) == "QUADRUPLE")
          precision_output_file = OutputTypeEnum::QUADRUPLE_PRECISION;
        else if (parameters(1) == "MULTIPLE")
          precision_output_file = OutputTypeEnum::MULTIPLE_PRECISION;
      }
    else if (description_field == "DirectoryOutput")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarOutputProblem" << endl;
	    cout << "DirectoryOutput needs one parameter, for instance :" << endl;
	    cout << "DirectoryOutput = chemin" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	if (parameters(0) == "AUTO")
          DOSSIER_output.clear();
        else
          DOSSIER_output = parameters(0);
      }
    else if (!description_field.compare("ElectricOrMagnetic"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarOutputProblem" << endl;
	    cout << "ElectricOrMagnetic needs more parameters, for instance :" << endl;
	    cout << "ElectricOrMagnetic = component" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	// ElectricOrMagnetic = numComponent
	// 0, 1, ..., n : only the i-th component of the solution will be written
	// -1 : all the components and gradients of the solution are written
	// -2 : all the components of the solution are written
	choice_field_output = to_num<int>(parameters(0));
      }
    else if (!description_field.compare("FineMeshLobatto"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarOutputProblem" << endl;
	    cout << "FineMeshLobatto needs more parameters, for instance :" << endl;
	    cout << "FineMeshLobatto = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	// if YES, subdivided meshes are split on Gauss-Lobatto points
	// instead of regular points (default choice)
	if (!parameters(0).compare("YES"))
	  output_lobatto_points = true;
      }
    else if (description_field == "MovePointsSurface")
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of VarOutputProblem" << endl;
	    cout << "MovePointsSurface needs more parameters, for instance :" << endl;
	    cout << "MovePointsSurface = YES component coef" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (parameters(0) == "YES")
          {
            this->move_points_on_output_mesh = true;
            this->component_used_to_translate = to_num<int>(parameters(1));
            this->coefficient_used_to_translate = to_num<Real_wp>(parameters(2));
          }
      }
    else if ((!description_field.compare("FileOutputMeshVolumetric"))
             ||(!description_field.compare("FileOutputMeshSurfacic")))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarOutputProblem" << endl;
	    cout << "FileOutputXXX needs more parameters, for instance :" << endl;
	    cout << "FileOutputMeshSurfacic = file_diffrac file_total component coef" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	ParamOutputClass output_par;
       
	for (int k = 0; k < parameters.GetM(); k++)
	  output_par.SetFileName(k, parameters(k));

      	this->output_mesh_param.PushBack(output_par);
      }
    else if (!description_field.compare("FileOutputSolution"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarOutputProblem" << endl;
	    cout << "FileOutputSolution needs more parameters, for instance :" << endl;
	    cout << "FileOutputSolution = file_diffrac file_total component coef" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        output_solution_file = true;
        for (int k = 0; k < parameters.GetM(); k++)
	  output_solution_param.SetFileName(k, parameters(k));
      }
    else if (description_field == "WriteSolutionQuadrature")
      {
        if (parameters(0) == "YES")
          {
            file_name_solution_quad = parameters(1);
            file_name_weights_quad = parameters(2);
            write_solution_on_quad_points = true;
          }
        else
          write_solution_on_quad_points = false;
      }
    else if (!description_field.compare("ParametersOutputSolution"))
      {        
        output_solution_param.SetInputData(parameters);
      }
    else if ((!description_field.compare("ParametersOutputMeshVolumetric"))
	     ||(!description_field.compare("ParametersOutputMeshSurfacic")))
      {
	int num_parameter_mesh = this->output_mesh_param.GetM()-1;
	this->output_mesh_param(num_parameter_mesh).SetInputData(parameters);
      }
    
  }


  //! writing the solution on file beginning with name_file
  /*!
    \param[in] U0 components of the solution
    \param[in] name_file prefix for each output file
   */
  template<class T>
  void VarOutputProblem_Base::
  WriteOutputFileReshaped(const Vector<Vector<T> >& U0, const string& name_file) const
  {
    int num_file = 0;
    for (int i = 0; i < this->output_grid_param.GetM(); i++)
      num_file = max(num_file, this->output_grid_param(i).GetNbFile());
    
    for (int i = 0; i < this->output_mesh_param.GetM(); i++)
      num_file = max(num_file, this->output_mesh_param(i).GetNbFile());
    
    // attention : on rajoute les noms de fichier dans les variables
    // output_grid_param et output_mesh_param : a modifier ulterieurement
    for (int i = 0; i < this->output_grid_param.GetM(); i++)
      const_cast<ParamOutputClass&>(this->output_grid_param(i)).SetFileName(num_file, name_file + to_str(i));
		
    for (int i = 0; i < this->output_mesh_param.GetM(); i++)
      const_cast<ParamOutputClass&>(this->output_mesh_param(i))
        .SetFileName(num_file, name_file + to_str(i+this->output_grid_param.GetM()));
    
    WriteOutputFileReshaped(U0, num_file);
  }
  
  
  //! reshapes U as a vector of vectors (each unknown is placed in a vector)
  template<class T>
  void VarOutputProblem_Base::ReshapeVector(const Vector<T>& U,
					    Vector<Vector<T> >& vecU, int nb_u) const
  {
    if (nb_u == -1)
      nb_u = var_problem.nb_unknowns;
    
    if (U.GetM() < var_problem.GetNbDof())
      nb_u = var_problem.nb_unknowns_scal;

    int m0 = 0, offset = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	m0 = var_problem.nb_unknowns_hdg;
	offset = var_problem.GetOffsetDofUnknown(m0);
      }
    
    vecU.Reallocate(nb_u);
    for (int i = 0; i < nb_u; i++)
      vecU(i).SetData(var_problem.GetNbDofUnknown(m0+i),
		      const_cast<T*>(&U(var_problem.GetOffsetDofUnknown(m0+i) - offset)));
  }
  

  //! writing the output indexed num_file in the matching output files
  /*!
    \param[in] U0 components of the solution on dofs
    \param[in] num_file output number
  */  
  template<class T>
  void VarOutputProblem_Base
  ::WriteOutputFile(const Vector<T>& U0, int num_file, int nb_u) const
  {
    if (nb_u == -1)
      nb_u = var_problem.nb_unknowns;
    
    Vector<int> num_unknown(nb_u);
    num_unknown.Fill();
    
    WriteOutputFile(U0, num_file, nb_u, num_unknown);
  }


  //! writing the output indexed num_file in the matching output files
  /*!
    \param[in] U0 components of the solution on dofs
    \param[in] num_file output number
  */  
  template<class T>
  void VarOutputProblem_Base
  ::WriteOutputFile(const Vector<T>& U0, int num_file, int nb_u,
                    const Vector<int>& num_unknown) const
  {
    // counting the number of unknowns stored in U0
    if (nb_u == -1)
      nb_u = var_problem.nb_unknowns;

    int m0 = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      m0 = var_problem.nb_unknowns_hdg;
    
    // U0 is reshaped as a vector of vectors
    Vector<Vector<T> > vecU(nb_u);
    for (int i = 0; i < nb_u; i++)
      {
        int offset = var_problem.GetOffsetDofUnknown(m0+num_unknown(i));
        //offset += var_problem.GetNbDof()*i/var_problem.nb_unknowns;
        vecU(i).SetData(var_problem.GetNbDofUnknown(m0+num_unknown(i)),
                        const_cast<T*>(&U0(offset)));
      }
    
    WriteOutputFileReshaped(vecU, num_file);
    
    for (int i = 0; i < nb_u; i++)
      vecU(i).Nullify();
  }


  //! writing the solution on file beginning with name_file
  /*!
    \param[in] U0 components of the solution
    \param[in] name_file prefix for each output file
  */
  template<class T>
  void VarOutputProblem_Base
  ::WriteOutputFile(const Vector<T>& U0, const string& name_file, int nb_u) const
  {
    if (nb_u == -1)
      nb_u = var_problem.nb_unknowns;
    
    Vector<int> num_unknown(nb_u);
    num_unknown.Fill();
    
    WriteOutputFile(U0, name_file, nb_u, num_unknown);
  }

  
  //! writing the solution on file beginning with name_file
  /*!
    \param[in] U0 components of the solution
    \param[in] name_file prefix for each output file
  */
  template<class T>
  void VarOutputProblem_Base
  ::WriteOutputFile(const Vector<T>& U0, const string& name_file, int nb_u, const Vector<int>& num_unknown) const
  {
    // counting the number of unknowns stored in U0
    if (nb_u == -1)
      nb_u = var_problem.nb_unknowns;

    int m0 = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      m0 = var_problem.nb_unknowns_hdg;

    // U0 is reshaped as a vector of vectors
    Vector<Vector<T> > vecU(nb_u);
    for (int i = 0; i < nb_u; i++)
      {
        int offset = var_problem.GetOffsetDofUnknown(m0+num_unknown(i));
        //offset += var_problem.GetNbDof()*i/var_problem.nb_unknowns;
        vecU(i).SetData(var_problem.GetNbDofUnknown(m0+num_unknown(i)),
                        const_cast<T*>(&U0(offset)));
      }
    
    WriteOutputFileReshaped(vecU, name_file);
    
    for (int i = 0; i < nb_u; i++)
      vecU(i).Nullify();
  }


  //! sets time in outputs
  void VarOutputProblem_Base
  ::ChangeTimeSnapshot(int nb_iter, const Real_wp& t)
  {
    for (int i = 0; i < output_grid_param.GetM(); i++)
      output_grid_param(i).ChangeTime(nb_iter, t);
    
    for (int i = 0; i < output_mesh_param.GetM(); i++)
      output_mesh_param(i).ChangeTime(nb_iter, t);
  }
  

  //! writing of the solution on all types of outputs (plane, lines, points ...)
  template<class T>
  void VarOutputProblem_Base::WriteOutputFile(const Vector<Vector<T> > & U0, int num_file, int nb_uloc) const
  {
    // counting the number of unknowns stored in U0
    int nb_u = 0;
    if (nb_uloc == -1)
      nb_uloc = var_problem.nb_unknowns;
    
    int nb_modes = var_boundary.GetNbModes();
    if (var_boundary.ModesNotStored())
      nb_modes = 1;
    
    int nb_vec = U0.GetM() / nb_modes;
    if (U0.GetM()%nb_modes != 0)
      {
	cout << "The size of U0 is not a multiple of the number of modes" << endl;
	abort();
      }

    nb_u = U0.GetM()*nb_uloc;

    int m0 = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      m0 = var_problem.nb_unknowns_hdg;

    // in vecU, modes are gathered together
    Vector<Vector<T> > vecU(nb_u);
    nb_u = 0; int offset = 0;
    for (int q = 0; q < nb_vec; q++)
      {
	for (int m = 0; m < nb_uloc; m++)
	  for (int p = 0; p < nb_modes; p++)
	    {
	      vecU(nb_u).SetData(var_problem.GetNbDofUnknown(m+m0),
				 const_cast<T*>(&U0(offset+p)(var_problem.GetOffsetDofUnknown(m0+m))));
	      nb_u++;
	    }
	
	offset += nb_modes;
      }
    
    WriteOutputFileReshaped(vecU, num_file);

    for (int i = 0; i < nb_u; i++)
      vecU(i).Nullify();
  }


  //! copies input parameters of another similar problem
  void VarOutputProblem_Base::CopyInputData(const VarOutputProblem_Base& var)
  {
    output_grid_param = var.output_grid_param;
    output_mesh_param = var.output_mesh_param;
    output_points_outside_param = var.output_points_outside_param;
    choice_field_output = var.choice_field_output;
    format_output_file = var.format_output_file;
    precision_output_file = var.precision_output_file;
    output_file_extension = var.output_file_extension;
    output_lobatto_points = var.output_lobatto_points;
    output_solution_param = var.output_solution_param;
    output_solution_file = var.output_solution_file;
    
    move_points_on_output_mesh = var.move_points_on_output_mesh;
    component_used_to_translate = var.component_used_to_translate;
    coefficient_used_to_translate = var.coefficient_used_to_translate;
    DOSSIER_output = var.DOSSIER_output;    
  }

  
  //! fills how the unknown u and gradient du can be written in a vtk file
  /*!
    By using scalar and vector outputs, we can write u and du
    on vtk files
  */
  void VarOutputProblem_Base
  ::GetVtkVolumeOutput(int nb_u, int nb_du, int dim_N, bool cplx_num,
		       IVect& type_u, Vector<IVect>& decomp_u, Vector<string>& name_u,
		       IVect& type_du, Vector<IVect>& decomp_du, Vector<string>& name_du) const
  {
    int cplx = 1;
    if (cplx_num)
      cplx = 2;
    
    type_u.Reallocate(nb_u*cplx);
    type_u.Fill(int(VtkDataTypeEnum::SCALAR));
    name_u.Reallocate(nb_u*cplx);
    decomp_u.Reallocate(nb_u*cplx);
    for (int i = 0; i < nb_u*cplx; i++)
      {
        decomp_u(i).Reallocate(1);
        decomp_u(i)(0) = i;
        name_u(i) = "u"+ to_str(i);
      }
    
    int nb_vec = nb_du;
    if (cplx == 1)
      {
	nb_vec = nb_du/dim_N;
	if (nb_du%dim_N > 0)
	  nb_vec++;
      }
    
    type_du.Reallocate(nb_vec*cplx);
    if (cplx == 1)
      type_du.Fill(int(VtkDataTypeEnum::VECTOR));
    else
      type_du.Fill(int(VtkDataTypeEnum::SCALAR));
    
    decomp_du.Reallocate(nb_vec*cplx);
    name_du.Reallocate(nb_vec*cplx);
    for (int i = 0; i < nb_vec*cplx; i++)
      {
	if (cplx == 1)
	  decomp_du(i).Reallocate(3);
	else
	  decomp_du(i).Reallocate(1);
	
	decomp_du(i).Fill(-1);
	name_du(i) = "du" + to_str(i);
        if (cplx == 2)
	  decomp_du(i)(0) = i;
        else
          for (int j = i*dim_N; j < min((i+1)*dim_N, nb_du); j++)
            decomp_du(i)(j-i*dim_N) = j;          
      }
  }
  

  //! fills how the unknown u and gradient du can be written in a vtk file
  /*!
    By using scalar and vector outputs, we can write u and du
    on vtk files
  */  
  void VarOutputProblem_Base
  ::GetVtkSurfaceOutput(int nb_u, int nb_du, int dim_N, bool cplx,
			IVect& type_u, Vector<IVect>& decomp_u, Vector<string>& name_u,
			IVect& type_du, Vector<IVect>& decomp_du, Vector<string>& name_du) const
  {
    GetVtkVolumeOutput(nb_u, nb_du, dim_N, cplx,
		       type_u, decomp_u, name_u, type_du, decomp_du, name_du);
  }
  
  
  /************************
   * VarOutputProblem_Dim *
   ************************/
  

#ifdef MONTJOIE_WITH_TWO_DIM    
  //! sets the number of modes for cyclic domain
  template<>
  void VarOutputProblem_Dim<Dimension2>::SetNbCyclicModes()
  {
    if (var_boundary.GetSymmetryType() == var_boundary.PERIODIC_THETA)    
      all_points_display.SetNbCyclicSections(var_boundary.GetNbPeriodicModesX());
    else
      all_points_display.SetNbCartesianSections(var_boundary.GetNbPeriodicModesX(),
						var_boundary.GetNbPeriodicModesY(), 1);    
  }
#endif
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! sets the number of modes for cyclic domain
  template<>
  void VarOutputProblem_Dim<Dimension3>::SetNbCyclicModes()
  {
    if (var_boundary.GetSymmetryType() == var_boundary.PERIODIC_THETA)    
      all_points_display.SetNbCyclicSections(var_boundary.GetNbPeriodicModesX());    
    else if (var_boundary.GetSymmetryType() == var_boundary.PERIODIC_ZTHETA)
      {
        all_points_display.SetNbCyclicSections(var_boundary.GetNbPeriodicModesX());
        all_points_display.SetNbCartesianSections(1, 1, var_boundary.GetNbPeriodicModesZ());
      }
    else
      all_points_display.SetNbCartesianSections(var_boundary.GetNbPeriodicModesX(),
						var_boundary.GetNbPeriodicModesY(),
						var_boundary.GetNbPeriodicModesZ());
    
  }
#endif
  

  //! adds memory used by the current object
  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>::GetMemoryUsed(map<string, size_t>& var) const
  {
    size_t taille = sizeof(*this);
    for (int i = 0; i < this->output_grid_param.GetM(); i++)
      taille += this->output_grid_param(i).GetMemorySize();

    for (int i = 0; i < this->output_mesh_param.GetM(); i++)
      taille += this->output_mesh_param(i).GetMemorySize();

    taille += this->output_solution_param.GetMemorySize();
    taille += this->all_points_display.GetMemorySize();
    for (int i = 0; i < var_grid.GetM(); i++)
      taille += var_grid(i).GetMemorySize();

    for (int i = 0; i < var_mesh_interp.GetM(); i++)
      taille += var_mesh_interp(i).GetMemorySize();
    
    var["Interpolation Grid"] = taille;
  }
  
  
  //! reading of parameters of the data file associated with outputs
  template<class Dimension> 
  void VarOutputProblem_Dim<Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarOutputProblem_Base::SetInputData(description_field, parameters);
    
    if (!description_field.compare("SismoGrille"))
      {
	// 3-D output on three planes
	int i = var_grid.GetM();
	var_grid.Resize(var_grid.GetM()+1);
	var_grid(i).SetInputData(description_field, parameters);
      }
    else if ((!description_field.compare("SismoPlane")) 
             || (!description_field.compare("SismoPlaneAxi")))
      {
	// output on a plane
	int i = var_grid.GetM();
	var_grid.Resize(var_grid.GetM()+1);
	var_grid(i).SetInputData(description_field, parameters);
      }
    else if ( (!description_field.compare("SismoPointsFile")) 
              || (!description_field.compare("SismoPointsFileAxi")))
      {
	// output on points given in a file
	int i = var_grid.GetM();
	var_grid.Resize(var_grid.GetM()+1);
	var_grid(i).SetInputData(description_field, parameters);
      }
    else if ( (!description_field.compare("SismoLine"))
              || (!description_field.compare("SismoLineAxi")))
      {
	// output on a plane
	int i = var_grid.GetM();
	var_grid.Resize(var_grid.GetM()+1);
	var_grid(i).SetInputData(description_field, parameters);
      }
    else if ( (!description_field.compare("SismoPoint"))
              || (!description_field.compare("SismoPointAxi")) )
      {
	// output on a point
	int i = var_grid.GetM();
	var_grid.Resize(var_grid.GetM()+1);
	var_grid(i).SetInputData(description_field, parameters);
      }
    else if ( (!description_field.compare("SismoCircle"))
              || (!description_field.compare("SismoCircleAxi")) )
      {
	// output on a circle
	int i = var_grid.GetM();
	var_grid.Resize(var_grid.GetM()+1);
	var_grid(i).SetInputData(description_field, parameters);
      }
    else if ( (!description_field.compare("SismoGrille3D")) )
      {
	// output on full 3-D grids
	int i = var_grid.GetM();
	var_grid.Resize(var_grid.GetM()+1);
	var_grid(i).SetInputData(description_field, parameters);
      }
    else if (!description_field.compare("SismoMeshVolumetric"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarOutputProblem" << endl;
	    cout << "SismoMeshVolumetric needs more parameters, for instance :" << endl;
	    cout << "SismoMeshVolumetric = AUTO" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	// output on a subdivided mesh
	MeshInterpolation<Dimension> var_interp;
	var_interp.SetType(var_interp.VOLUME_MESH);
	if (!parameters(0).compare("AUTO"))
	  var_interp.SetRegularSubdivisions(var_problem.GetDefaultOrder());
	else
	  var_interp.SetRegularSubdivisions(to_num<int>(parameters(0)));
	
	var_mesh_interp.PushBack(var_interp);
      }
    else if (!description_field.compare("SismoMeshSurfacic"))
      {
	// output on a subdivided mesh of the surface 
	MeshInterpolation<Dimension> var_interp;
	var_interp.SetType(var_interp.SURFACE_MESH);
        if (parameters.GetM() <= 2)
          {
	    cout << "For the keyword SismoMeshSurfacic" << endl;
            cout << "You have to specify the number of subdivisions"
                 << ", the type of surface and references, for instance: " << endl;
            cout << "SismoMeshSurfacic = AUTO REFERENCE 1 2 3" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
	if (!parameters(0).compare("AUTO"))
	  var_interp.SetRegularSubdivisions(var_problem.GetDefaultOrder());
	else
	  var_interp.SetRegularSubdivisions(to_num<int>(parameters(0)));
	
	if (!parameters(1).compare("BODY"))
	  {
	    var_interp.SetType(var_interp.BODY_MESH);
	    IVect ref_surf(parameters.GetM()-2);
	    for (int i = 2; i < parameters.GetM(); i++)
	      ref_surf(i-2) = to_num<int>(parameters(i));
            
            var_interp.SetSurfaceReference(ref_surf);
	  }
	else if (!parameters(1).compare("CONDITION"))
	  {
	    var_interp.SetType(var_interp.CONDITION_MESH);
	    IVect ref_surf(parameters.GetM()-2);
	    for (int i = 2; i < parameters.GetM(); i++)
	      ref_surf(i-2) = to_num<int>(parameters(i));

            var_interp.SetSurfaceReference(ref_surf);
	  }
	else if (!parameters(1).compare("REFERENCE"))
	  {
	    var_interp.SetType(var_interp.REFERENCE_MESH);
            IVect ref_surf;
	    if (!parameters(2).compare("ALL"))
	      {
                int nb_max = 100;
		if (parameters.GetM() > 3)
                  nb_max = to_num<int>(parameters(3));
                
                ref_surf.Reallocate(nb_max);
                ref_surf.Fill();
	      }
	    else
	      {
		ref_surf.Reallocate(parameters.GetM()-2);
		for (int i = 2; i < parameters.GetM(); i++)
		  ref_surf(i-2) = to_num<int>(parameters(i));
	      }
            
            var_interp.SetSurfaceReference(ref_surf);
	  }
	
	var_mesh_interp.PushBack(var_interp);
      }
    else if (!description_field.compare("SismoOutsidePoints"))
      {                
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of VarOutputProblem" << endl;
	    cout << "SismoOutsidePoints needs more parameters, for instance :" << endl;
	    cout << "SismoOutsidePoints = file_name file_diffrac file_total" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        // premier parametre : fichier ou les points sont stockes
        string file_points = parameters(0);
        int nb = 1; VectR_N PointsOutside;
        if (file_points == "CIRCLE")
          {
	    if (parameters.GetM() <= 4+Dimension::dim_N)
	      {
		cout << "In SetInputData of VarOutputProblem" << endl;
		cout << "SismoOutsidePoints needs more parameters, for instance :" << endl;
		cout << "SismoOutsidePoints = CIRCLE x0 y0 z0 r N file_diffrac file_total" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

            // centre du cercle (dans le plan Oz)
            R_N center;
	    for (int k = 1; k <= Dimension::dim_N; k++)
	      center(k-1) = to_num<Real_wp>(parameters(k));
	    
            // rayon
	    nb = Dimension::dim_N+1;
            Real_wp radius = to_num<Real_wp>(parameters(nb++));
            
            // nombre de points
            int nb_angles = to_num<int>(parameters(nb++));
	    
            PointsOutside.Reallocate(nb_angles);
            for (int i = 0; i < nb_angles; i++)
              {
		PointsOutside(i) = center;
                
                Real_wp teta = 2.0*pi_wp*Real_wp(i)/nb_angles;
		
                PointsOutside(i)(0) += radius*cos(teta);
                PointsOutside(i)(1) += radius*sin(teta);
              }
          }
        else if (file_points == "CIRCLE_STEREO")
          {
	    if (parameters.GetM() <= 5+Dimension::dim_N)
	      {
		cout << "In SetInputData of VarOutputProblem" << endl;
		cout << "SismoOutsidePoints needs more parameters, for instance :" << endl;
		cout << "SismoOutsidePoints "
		     << "= CIRCLE_STEREO x0 y0 z0 r N dist file_diffrac file_total" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

            // centre du cercle (dans le plan Oz)
	    R_N center;
	    for (int k = 1; k <= Dimension::dim_N; k++)
	      center(k-1) = to_num<Real_wp>(parameters(k));
            
            // rayon
	    nb = Dimension::dim_N+1;
            Real_wp radius = to_num<Real_wp>(parameters(nb++));
            
            // nombre de points
            int nb_angles = to_num<int>(parameters(nb++));
	    
	    // distance entre oreille gauche et droite
	    Real_wp dist = to_num<Real_wp>(parameters(nb++));
            Real_wp dteta = dist / radius;
	    
            PointsOutside.Reallocate(2*nb_angles);
            for (int i = 0; i < nb_angles; i++)
              {
		PointsOutside(2*i) = center;
		PointsOutside(2*i+1) = center;
		
                Real_wp teta = 2.0*pi_wp*Real_wp(i)/nb_angles;

                PointsOutside(2*i)(0) += radius*cos(teta);
		PointsOutside(2*i)(1) += radius*sin(teta);

                PointsOutside(2*i+1)(0) += radius*cos(teta+dteta);                
                PointsOutside(2*i+1)(1) += radius*sin(teta+dteta);
              }
          }
	else if (file_points == "LINE")
	  {
	    if (parameters.GetM() <= 3+2*Dimension::dim_N)
	      {
		cout << "In SetInputData of VarOutputProblem" << endl;
		cout << "SismoOutsidePoints needs more parameters, for instance :" << endl;
		cout << "SismoOutsidePoints "
		     << "= LINE x0 y0 z0 x1 y1 z1 N file_diffrac file_total" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    // extremites de la ligne
	    R_N pt0, pt1;
	    nb = 1;
	    for (int k = 0; k < Dimension::dim_N; k++)
	      pt0(k) = to_num<Real_wp>(parameters(nb++));
	    
	    for (int k = 0; k < Dimension::dim_N; k++)
	      pt1(k) = to_num<Real_wp>(parameters(nb++));

	    // nombre de points
            int nb_points = to_num<int>(parameters(nb++));

	    // creation des points sur la ligne
	    PointsOutside.Reallocate(nb_points);
            for (int i = 0; i < nb_points; i++)
              {
		Real_wp L(0);
		if (nb_points > 1)
		  L = Real_wp(i) / (nb_points-1);
		
		PointsOutside(i) = (1.0 - L)*pt0 + L*pt1;
	      }
	  }
	else if (file_points == "PLANE")
	  {
	    if (parameters.GetM() <= 4+3*Dimension::dim_N)
	      {
		cout << "In SetInputData of VarOutputProblem" << endl;
		cout << "SismoOutsidePoints needs more parameters, for instance :" << endl;
		cout << "SismoOutsidePoints "
		     << "= PLANE x0 y0 z0 x1 y1 z1 x2 y2 z2 Nx Ny file_diffrac file_total" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    // trois plans du parallelogramme
	    R_N pt0, pt1, pt2;
	    nb = 1;
	    for (int k = 0; k < Dimension::dim_N; k++)
	      pt0(k) = to_num<Real_wp>(parameters(nb++));
	    
	    for (int k = 0; k < Dimension::dim_N; k++)
	      pt1(k) = to_num<Real_wp>(parameters(nb++));

	    for (int k = 0; k < Dimension::dim_N; k++)
	      pt2(k) = to_num<Real_wp>(parameters(nb++));

	    // nombre de points en x et y
            int Nx = to_num<int>(parameters(nb++));
	    int Ny = to_num<int>(parameters(nb++));

	    R_N vec_u, vec_v;
	    vec_u = pt1 - pt0;
	    vec_v = pt2 - pt0;
	    
	    // creation des points sur la ligne
	    PointsOutside.Reallocate(Nx*Ny);
            for (int i = 0; i < Nx; i++)
	      for (int j = 0; j < Ny; j++)
		{
		  int num = j*Nx + i;
		  Real_wp Lx(0), Ly(0);
		  if (Nx > 1)
		    Lx = Real_wp(i) / (Nx-1);
		  
		  if (Ny > 1)
		    Ly = Real_wp(j) / (Ny-1);
		  
		  PointsOutside(num) = pt0 + Lx*vec_u + Ly*vec_v;
		}
	  }
        else if (file_points == "POINT")
          {
	    if (parameters.GetM() <= 2+Dimension::dim_N)
	      {
		cout << "In SetInputData of VarOutputProblem" << endl;
		cout << "SismoOutsidePoints needs more parameters, for instance :" << endl;
		cout << "SismoOutsidePoints "
		     << "= POINT x0 y0 z0 file_diffrac file_total" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

            // point
            typename Dimension::R_N point;
	    for (int k = 0; k < Dimension::dim_N; k++)
              point(k) = to_num<Real_wp>(parameters(k+1));
	    
            PointsOutside.Reallocate(1);
            PointsOutside(0) = point;
            nb = Dimension::dim_N+1;
          }
        else
          PointsOutside.ReadText(file_points);
        
        // deuxieme/troisieme parametre : fichiers de sortie
        string ext = this->output_file_extension;
        this->output_points_outside_param.SetDiffractedFieldFile(GetBaseString(parameters(nb++)) + ext);
        this->output_points_outside_param.SetTotalFieldFile(GetBaseString(parameters(nb++)) + ext);

        // dernier parametre : pas de temps
	Real_wp dt(0);
        if (parameters.GetM() > nb)
          dt = to_num<Real_wp>(parameters(nb++));
	
	this->SetComputationFarPoints(PointsOutside, dt);
      }
    else if ((!description_field.compare("FileOutputGrille")) 
             || (!description_field.compare("FileOutputGrille3D")) ||
	     (!description_field.compare("FileOutputPlane")) 
             || (!description_field.compare("FileOutputPlaneAxi")) ||
             (!description_field.compare("FileOutputLine")) 
             || (!description_field.compare("FileOutputLineAxi")) || 
	     (!description_field.compare("FileOutputPoint")) 
             || (!description_field.compare("FileOutputPointAxi")) || 
             (!description_field.compare("FileOutputCircle")) 
             || (!description_field.compare("FileOutputCircleAxi")) ||
             (!description_field.compare("FileOutputPointsFile")) 
             || (!description_field.compare("FileOutputPointsFileAxi")))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarOutputProblem" << endl;
	    cout << "FileOutputXXX needs more parameters, for instance :" << endl;
	    cout << "FileOutputPlane = file_diffrac file_total component coef" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	int type = GridInterpolationFull<Dimension>::GetType(description_field);
	int nb = this->GetIndexOutputFiles(type);
	
	// Files of storage of solution on a regular grid. Diffracted and total field
	string ext = this->output_file_extension;
	
	for (int k = 0; k < parameters.GetM(); k++)
	  this->output_grid_param(nb).SetFileName(k, GetBaseString(parameters(k)) + ext);

      }
    else if (!description_field.compare("ParametersOutputGrille"))
      {
	int nb = GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                     GetType("FileOutputGrille"));
	this->output_grid_param(nb).SetInputData(parameters);
      }
    else if (!description_field.compare("ParametersOutputGrille3D"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputGrille3D"));
	this->output_grid_param(nb).SetInputData(parameters);
      }
    else if (!description_field.compare("ParametersOutputPlane"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputPlane"));
	this->output_grid_param(nb).SetInputData(parameters);
      }
    else if (!description_field.compare("ParametersOutputPlaneAxi"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputPlaneAxi"));
	this->output_grid_param(nb).SetInputData(parameters);
      }
    else if (!description_field.compare("ParametersOutputLine"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputLine"));
	this->output_grid_param(nb).SetInputData(parameters);
      }
    else if (!description_field.compare("ParametersOutputLineAxi"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputLineAxi"));
	this->output_grid_param(nb).SetInputData(parameters);
      }
    else if (!description_field.compare("ParametersOutputPoint"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputPoint"));
        
	this->output_grid_param(nb).SetInputData(parameters, true);
      }
    else if (!description_field.compare("ParametersOutputPointAxi"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputPointAxi"));
	this->output_grid_param(nb).SetInputData(parameters, true);
      }
    else if (!description_field.compare("ParametersOutputCircle"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputCircle"));
	this->output_grid_param(nb).SetInputData(parameters);
      }
    else if (!description_field.compare("ParametersOutputCircleAxi"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputCircleAxi"));
	this->output_grid_param(nb).SetInputData(parameters);
      }
    else if (!description_field.compare("ParametersOutputPointsFile"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputPointsFile"));
	this->output_grid_param(nb).SetInputData(parameters);
      }
    else if (!description_field.compare("ParametersOutputPointsFileAxi"))
      {
	int nb = this->GetIndexOutputFiles(GridInterpolationFull<Dimension>::
                                           GetType("FileOutputPointsFileAxi"));
	this->output_grid_param(nb).SetInputData(parameters);
      }    
  }
  
  
  //! Initialization of the display grids
  /*!
    Points for predefined grids are computed and
    initialization before the localization in the mesh.
    Subdivided meshes for output on meshes, are written
  */
  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>::InitVarGrid()
  {
    // boundaries of the grid are set to boundaries of the domain if required
    for (int i = 0; i < var_grid.GetM(); i++)
      {
	var_grid(i).SetXmin0(var_problem.GetXmin());
	var_grid(i).SetXmax0(var_problem.GetXmax());
	var_grid(i).SetYmin0(var_problem.GetYmin());
	var_grid(i).SetYmax0(var_problem.GetYmax());
	var_grid(i).SetZmin0(var_problem.GetZmin());
	var_grid(i).SetZmax0(var_problem.GetZmax());
      }
    
    // components of unknowns to store in output files 
    // are selected with the global variable choice_field_output
    for (int i = 0; i < this->output_grid_param.GetM(); i++)
      this->output_grid_param(i).SetComponent(this->choice_field_output);
    
    for (int i = 0; i < this->output_mesh_param.GetM(); i++)
      this->output_mesh_param(i).SetComponent(this->choice_field_output);
    
    this->output_solution_param.SetComponent(this->choice_field_output);
    
    // For each line of the data file : SismoMeshVolumetric = ...
    // there should be an associated line : FileOutputMeshVolumetric = ...
    if (this->output_mesh_param.GetM() != var_mesh_interp.GetM())
      {
	cout<<"not same number of FileOutputMesh and SismoMesh in data file"<<endl;
	abort();
      }
    
    // loop over all the outputs on meshes
    for (int i = 0; i < var_mesh_interp.GetM(); i++)
      {
	string ext(".mesh");
	if ((this->format_output_file == OutputTypeEnum::ASCII_VTK)
	    || (this->format_output_file == OutputTypeEnum::BINARY_VTK))
	  ext = string(".vtk");
	
	// creation of a divided mesh 
	var_mesh_interp(i).SetFileName(GetBaseString(this->output_mesh_param(i).
                                                     GetTotalFieldFile())
                                       + ext); 
	
	if (var_mesh_interp(i).IsVolumetric())
	  InitVarMeshVolumetric(var_mesh_interp(i));
	else
	  InitVarMeshSurfacic(var_mesh_interp(i));
      }
    
    all_points_display.SetXmin(var_problem.mesh.GetXmin());
    all_points_display.SetXmax(var_problem.mesh.GetXmax());
    all_points_display.SetYmin(var_problem.mesh.GetYmin());
    all_points_display.SetYmax(var_problem.mesh.GetYmax());
    all_points_display.SetZmin(var_problem.mesh.GetZmin());
    all_points_display.SetZmax(var_problem.mesh.GetZmax());
    
#ifdef SELDON_WITH_MPI
    all_points_display.SetCommunicator(var_problem.comm_group_mode);
#endif

    this->SetNbCyclicModes();
    
    // predefined grids are generated
    // all_points_display is an union of all the predefined grids
    // so that the localization step is done once for the grid all_points_display
    for (int i = 0; i < var_grid.GetM(); i++)
      var_grid(i).InitGrid(all_points_display);

    // pre-localization of the elements of the mesh so that
    // each new point to search can be fastly found 
    // implementation of this method is in GridInterpolation.cxx
    all_points_display.InitInterpolationGrid(var_problem.mesh);
  }

  
  //! computation of the interpolation grid
  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>::ComputeVarGrid()
  {
    // signaling that this operation has already been performed
    this->display_solution = true;
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    if ((rank_proc == 0)||(var_problem.print_level >= 10))
      if (var_problem.print_level>=0)
	cout << rank_proc << " We compute grid now " << endl;
    
    // localization of the points of the predefined grids inside the mesh
    // implementation of this method in file GridInterpolation.cxx
    all_points_display.LocalizePoints(var_problem.mesh);

    // releasing memory used to prelocalize points
    if (!keep_grid_localization_array)
      all_points_display.ClearPrelocalizationArrays();
    
    // the grids are compressed
    all_points_display.CompressGrid(var_grid);
    
    if ((rank_proc == 0)||(var_problem.print_level >= 10))
      if (var_problem.print_level >= 1)
	cout<<rank_proc<<" Interpolation Grid OK"<<endl;
  }
  
  
  //! Construction of a subdivided mesh for outputs on medit
  /*!
    \param[in,out] var_interp definition of the subdivision to be done
  */
  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>
  ::InitVarMeshVolumetric(MeshInterpolation<Dimension>& var_interp)
  {
    // use of Gauss-Lobatto points if required
    if (this->output_lobatto_points)
      var_interp.SetLobattoSubdivisions(var_interp.GetNbSubdivisions());

#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif
    
    Mesh<Dimension> mesh_subdiv;
    mesh_subdiv = var_problem.mesh;

    // the mesh is subdivided
    const VectReal_wp& step_subdiv = var_interp.GetSubdivisionStep();
    
    Vector<IVect> NodleSurf;
    mesh_subdiv.SubdivideMesh(step_subdiv, var_interp.PointsReferenceSurface(),
                              var_interp.PointsReferenceVolume(),
                              var_interp.NumberingArray(), NodleSurf);
    
    // the mesh is written in a file
    string file_name = GetBaseString(var_interp.GetFileName())+to_str(rank_proc)
      +'.'+GetExtension(var_interp.GetFileName());
   
    if (nb_proc == 1)
      file_name = GetBaseString(var_interp.GetFileName())+'.'+GetExtension(var_interp.GetFileName());
     
    bool double_prec = true, ascii = true;
    if ((this->format_output_file == OutputTypeEnum::BINARY_MATLAB)
        || (this->format_output_file == OutputTypeEnum::BINARY_VTK))
      ascii = false;
    
    if (this->precision_output_file == OutputTypeEnum::SINGLE_PRECISION)
      double_prec = false;

    mesh_subdiv.Write(this->DOSSIER_output + file_name, double_prec, ascii);
    
    // definition of the projection operator
    // (from an original element to a subdivided element)
    var_interp.InitProjectionVolume(var_problem.mesh);
    
    // we wait that all the meshes are written
#ifdef SELDON_WITH_MPI
    MPI_Barrier(var_problem.comm_group_mode);
#endif
    
    if ((nb_proc > 1) &&(rank_proc == 0))
      {
	Mesh<Dimension> mesh_all;
	// merging all the meshes into a single big mesh
	for (int i = 0; i < nb_proc; i++)
	  {
	    file_name = GetBaseString(var_interp.GetFileName())+to_str(i)
	      +'.'+GetExtension(var_interp.GetFileName());
	    
            mesh_subdiv.Read(this->DOSSIER_output + file_name);
            // changing references to remove "neighbor" faces
            IVect ref_neighbor = var_problem.mesh.GetOriginalNeighborReference();
            for (int i = 0; i < mesh_subdiv.GetNbBoundaryRef(); i++)
              {
                int ref = mesh_subdiv.BoundaryRef(i).GetReference();
                if (ref > 0)
                  if (ref_neighbor(ref) == 0)
                    mesh_subdiv.BoundaryRef(i).SetReference(0);
              }
            
	    // the temporary file is removed
	    std::remove(file_name.data());
	    mesh_all.AppendMesh(mesh_subdiv, false);
	  }

	file_name = GetBaseString(var_interp.GetFileName())+'.'+GetExtension(var_interp.GetFileName());
	mesh_all.Write(this->DOSSIER_output + file_name, double_prec, ascii);        
      }
    
  }
  
  
  //! Extraction and numbering of a surfacic mesh from the volumetric mesh
  /*!
    \param[in,out] var_interp variables defining the surfacic mesh
    \param[in] dim dummy argument to select this method only in 3-D case
  */
  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>
  ::InitVarMeshSurfacic(MeshInterpolation<Dimension>& var_interp)
  {
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    // use of Gauss-Lobatto points if required
    if (this->output_lobatto_points)
      var_interp.SetLobattoSubdivisions(var_interp.GetNbSubdivisions());
    
    // construction of ref_cond, to regroup all the surfaces that need to be extracted
    IVect ref_cond(var_problem.mesh.GetNbReferences()+1); ref_cond.Fill(0);

    switch (var_interp.GetType())
      {
      case MeshInterpolation<Dimension>::SURFACE_MESH :
        {
          // in this case, all the referenced boundaries are extracted
          ref_cond.Fill(1);
          ref_cond(0) = 0;
        }
        break;
      case MeshInterpolation<Dimension>::REFERENCE_MESH :
        {
          // only referenced boundary, whose reference is contained in arry ref_surf
          // are extracted
          const IVect& ref_surf = var_interp.GetSurfaceReference();
	  for (int i = 0; i < ref_surf.GetM(); i++)
            if (ref_surf(i) < ref_cond.GetM())
              ref_cond(ref_surf(i)) = 1;
        }
        break;
      case MeshInterpolation<Dimension>::BODY_MESH :
        {
          // only boundaries, which belongs to bodies of reference ref_surf
          // are extracted
          const IVect& ref_surf = var_interp.GetSurfaceReference();
          for (int i = 0; i < ref_surf.GetM(); i++)
            for (int j = 0; j < var_problem.mesh.GetNbReferences(); j++)
              if (var_problem.mesh.GetBodyNumber(j) == ref_surf(i))
                ref_cond(j) = 1;
          
        }
        break;
      case MeshInterpolation<Dimension>::CONDITION_MESH :
        {
          // only boundaries, which belongs to boundary conditions of reference ref_surf
          // are extracted
          const IVect& ref_surf = var_interp.GetSurfaceReference();
          for (int i = 0; i < ref_surf.GetM(); i++)
            for (int j = 0; j < var_problem.mesh.GetNbReferences(); j++)
              if (var_problem.mesh.GetBoundaryCondition(j) == ref_surf(i))
                ref_cond(j) = 1;
        }
        break;
      }
    
    // removing neighbor references
    IVect ref_neighbor = var_problem.mesh.GetOriginalNeighborReference();
    for (int k = 0; k < ref_cond.GetM(); k++)
      {
        int ref_n = ref_neighbor(k);
        if (ref_n == 0)
          ref_cond(k) = 0;
        else if (ref_n > 0)
          if (ref_cond(ref_n) == 1)
            ref_cond(k) = 1;
      }
    
    Mesh<Dimension> mesh_subdiv;
    var_interp.ComputeSurfaceMesh(ref_cond, var_problem.mesh, mesh_subdiv);
    
    var_interp.InitProjectionSurface(var_problem.mesh);

    bool double_prec = true, ascii = true;
    if ((this->format_output_file == OutputTypeEnum::BINARY_MATLAB)
        || (this->format_output_file == OutputTypeEnum::BINARY_VTK))
      ascii = false;
    
    if (this->precision_output_file == OutputTypeEnum::SINGLE_PRECISION)
      double_prec = false;
    
    // writing of the mesh
    string file_name = GetBaseString(var_interp.GetFileName())+to_str(rank_proc)
      + '.' + GetExtension(var_interp.GetFileName());

    if (nb_proc == 1)
      file_name = GetBaseString(var_interp.GetFileName())
	+ '.' + GetExtension(var_interp.GetFileName());
    
    mesh_subdiv.Write(this->DOSSIER_output + file_name, double_prec, ascii); 
       
#ifdef SELDON_WITH_MPI 
    // we wait that all the meshes are written
    MPI_Barrier(var_problem.comm_group_mode);
#endif
    
    if ((nb_proc > 1) && (rank_proc == 0))
      {
	Mesh<Dimension> mesh_all;
	// appending all the meshes to a unique mesh
	for (int i = 0; i < nb_proc; i++)
	  {
	    file_name = GetBaseString(var_interp.GetFileName())+to_str(i)
              + '.' + GetExtension(var_interp.GetFileName());
	    
            mesh_subdiv.Read(this->DOSSIER_output + file_name);
	    // the temporary file is removed
	    std::remove(file_name.data());
	    mesh_all.AppendMesh(mesh_subdiv, false);
	  }

	file_name = GetBaseString(var_interp.GetFileName()) + '.' + GetExtension(var_interp.GetFileName());
	mesh_all.Write(this->DOSSIER_output + file_name, double_prec, ascii);
      }
  }
  
  
  //! interpolation of U0 over a grid, result is placed in trace_u
  /*!
    \param[in] U0 components of the solution on dofs
    \param[in] var_interp description of points where we want to know the solution
    \param[in] var_gr interpolation grid localizing all the points
    \param[out] trace_u solution on asked points
  */
  template<class Dimension> template<class T>
  void VarOutputProblem_Dim<Dimension>::
  ComputeInterpolationU_GradU(const Vector<Vector<T> > & U0,
			      const GridInterpolation<Dimension> & var_interp,
			      const GridInterpolationFull<Dimension>& var_gr,
			      Vector<Vector<T> >& trace_u, Vector<Vector<T> >& trace_grad_u,
			      int add_total_field, const Real_wp& t,
			      Vector<Vector<T> >& trace_u_diff, Vector<Vector<T> >& trace_grad_diff,
			      bool compute_grad) const
  {
    const IVect& list_points = var_gr.GetPointNumber();
    const Vector<IVect>& list_points_proc = var_gr.GetPointNumberAllProc();
    
    if (var_gr.GetNbGlobalPoints() <= 0)
      return;
    
    int nb_modes = var_boundary.GetNbModes();
    int nb_proc = var_problem.GetNbProcPerMode();
    if (var_boundary.ModesNotStored())
      nb_modes = 1;
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    if (var_problem.GetNbComponentsGradientAll() <= 0)
      compute_grad = false;
    
    int nnz = list_points.GetM();
    int nb_vec = U0.GetM()/nb_modes;
    int nb_comp_all = var_problem.GetNbComponentsAll(nb_vec);
    int nb_comp_grad = var_problem.GetNbComponentsGradientAll(nb_vec);
    trace_u.Reallocate(nb_comp_all);
    if (compute_grad)
      trace_grad_u.Reallocate(nb_comp_grad);
   
    Vector<T> trace_vec, trace_vec_diff, trace_grad_vec, trace_grad_vec_diff;
    
    // first step : we compute the values on U on each point of the grid for each processor
    // these values are stored in vector trace_vec
    this->ComputeInterpolationUloc(U0, var_interp, trace_vec, trace_grad_vec,
				   list_points, nnz, compute_grad);    

    // intermediary step, we compute the diffracted field/total field if needed
    if (add_total_field != 0)
      {
        trace_u_diff.Reallocate(nb_comp_all);
        trace_vec_diff.Reallocate(nnz*nb_comp_all);
        FillZero(trace_vec_diff);
	FillZero(trace_u_diff);

        if (compute_grad)
	  {
	    trace_grad_diff.Reallocate(nb_comp_grad);
	    trace_grad_vec_diff.Reallocate(nnz*nb_comp_grad);
	    FillZero(trace_grad_vec_diff);
	    FillZero(trace_grad_diff);
	  }
	
	IncidentWaveProjector<T, Dimension>& incident_wave
	  = *var_source.GetIncidentWaveProjector(0, T(0));

	incident_wave.SetTime(t, 0);
        
        int nb = 0;
        for (int i1 = 0; i1 < list_points.GetM(); i1++)
          {
            int i = list_points(i1); R_N pt_glob;
            Vector<T> u_inc(nb_comp_all), grad_u_inc(nb_comp_grad); T vloc;
	    u_inc.Fill(0);
            if (i >= 0)
              {
                pt_glob = var_interp.GetGlobalCoordinate(i);
		if (compute_grad)
		  incident_wave.EvaluateFunctionGradient(-1, -1, pt_glob, u_inc, grad_u_inc);
		else
		  incident_wave.EvaluateFunction(-1, -1, pt_glob, u_inc);

                if (add_total_field == 1)
                  {
                    for (int j = 0; j < nb_comp_all; j++)
                      {
                        vloc = trace_vec(nb+j*nnz);
                        trace_vec_diff(nb+j*nnz) = vloc - u_inc(j);
                      }

		    if (compute_grad)
		      for (int j = 0; j < nb_comp_grad; j++)
			{
			  vloc = trace_grad_vec(nb + j*nnz);
			  trace_grad_vec_diff(nb + j*nnz) = vloc - grad_u_inc(j);
			}
                  }
                else
                  {
                    for (int j = 0; j < nb_comp_all; j++)
                      {
                        vloc = trace_vec(nb + j*nnz);
                        trace_vec_diff(nb + j*nnz) = vloc;
                        trace_vec(nb + j*nnz) = vloc + u_inc(j);
                      }

		    if (compute_grad)
		      for (int j = 0; j < nb_comp_grad; j++)
			{
			  vloc = trace_grad_vec(nb + j*nnz);
			  trace_grad_vec_diff(nb + j*nnz) = vloc;
			  trace_grad_vec(nb + j*nnz) = vloc + grad_u_inc(j);
			}
                  }
		
                nb++;
              }
          }
      }

    // second step : each processor sends its own data to the root processor
#ifdef SELDON_WITH_MPI
    Vector<int64_t> xtmp;
    MPI_Status status;
    if ((nnz > 0)&&(rank_proc != 0))
      {
        MpiSsend(var_problem.comm_group_mode, trace_vec, xtmp, trace_vec.GetM(), 0, 31);
	if (compute_grad)
          MpiSsend(var_problem.comm_group_mode, trace_grad_vec, xtmp, trace_grad_vec.GetM(), 0, 36);
        
        if (add_total_field != 0)
          {
	    MpiSsend(var_problem.comm_group_mode, trace_vec_diff, xtmp, trace_vec_diff.GetM(), 0, 32);
	    if (compute_grad)
              MpiSsend(var_problem.comm_group_mode, trace_grad_vec_diff, xtmp, trace_grad_vec_diff.GetM(), 0, 37);
	  }
      }
#endif
    
    Vector<T> trace_sol, trace_sol_diff;
    Vector<T> trace_grad_sol, trace_grad_sol_diff;
    if (rank_proc == 0)
      {
	// third step :
	// root processor retrieves the values computed on each processor
	// in the variable trace_sol
	int nb_pts_all = var_gr.GetNbGlobalPoints();
	for (int p = 0; p < nb_comp_all; p++)
	  {
	    trace_u(p).Reallocate(nb_pts_all);
	    trace_u(p).Zero();
	  }
	
        if (add_total_field != 0)
          for (int p = 0; p < nb_comp_all; p++)
            {
	      trace_u_diff(p).Reallocate(nb_pts_all);
	      trace_u_diff(p).Zero();
	    }
	
	if (compute_grad)
	  for (int p = 0; p < nb_comp_grad; p++)
	    {
	      trace_grad_u(p).Reallocate(nb_pts_all);
	      trace_grad_u(p).Zero();
	    }
	
        if (compute_grad && (add_total_field != 0))
          for (int p = 0; p < nb_comp_grad; p++)
            {
	      trace_grad_diff(p).Reallocate(nb_pts_all);
	      trace_grad_diff(p).Zero();
	    }
	
	for (int i = 0; i < nb_proc; i++)
	  if (list_points_proc(i).GetM() > 0)
	    {
	      trace_sol.Reallocate(list_points_proc(i).GetM()*nb_comp_all);
	      if (compute_grad)
		trace_grad_sol.Reallocate(list_points_proc(i).GetM()*nb_comp_grad);
	      
#ifdef SELDON_WITH_MPI
	      if (i != rank_proc)
                {
                  MpiRecv(var_problem.comm_group_mode, trace_sol, xtmp,
                          trace_sol.GetM(), i, 31, status);
		  
		  if (compute_grad)
		    MpiRecv(var_problem.comm_group_mode, trace_grad_sol, xtmp,
			    trace_grad_sol.GetM(), i, 36, status);
		  
                  if (add_total_field != 0)
                    {
                      trace_sol_diff.Reallocate(list_points_proc(i).GetM()*nb_comp_all);
                      MpiRecv(var_problem.comm_group_mode, trace_sol_diff, xtmp,
                              trace_sol_diff.GetM(), i, 32, status);

		      if (compute_grad)
			{
			  trace_grad_sol_diff.Reallocate(list_points_proc(i).GetM()*nb_comp_grad);
			  MpiRecv(var_problem.comm_group_mode, trace_grad_sol_diff, xtmp,
				  trace_grad_sol_diff.GetM(), i, 37, status);
			}
                    }
                }
              else
#endif
		{
                  trace_sol = trace_vec;
		  if (compute_grad)
		    trace_grad_sol = trace_grad_vec;
		  
                  if (add_total_field != 0)
                    {
		      trace_sol_diff = trace_vec_diff;
		      if (compute_grad)
			trace_grad_sol_diff = trace_grad_vec_diff;
		    }
                }
	      
	      // now we convert trace_sol in trace_u
	      for (int j = 0; j < list_points_proc(i).GetM(); j++)
		{
		  int nnz2 = list_points_proc(i).GetM();
		  int num_pt = list_points_proc(i)(j);
		  for (int p = 0; p < nb_comp_all; p++)
		    trace_u(p)(num_pt) = trace_sol(j + nnz2*p);
		  
		  if (compute_grad)
		    for (int p = 0; p < nb_comp_grad; p++)
		      trace_grad_u(p)(num_pt) = trace_grad_sol(j + nnz2*p);
		  
		  if (add_total_field != 0)
		    {
		      for (int p = 0; p < nb_comp_all; p++)
			trace_u_diff(p)(num_pt) = trace_sol_diff(j + nnz2*p);
		      
		      if (compute_grad)
			for (int p = 0; p < nb_comp_grad; p++)
			  trace_grad_diff(p)(num_pt) = trace_grad_sol_diff(j + nnz2*p);
		    }
		}
	    }
      }
  }
  

  //! computes the interpolation of U0 on points of grid num_grid
  /*!
    \param[in] compute_grad if true the gradient is also computed
    \param[out] trace_u trace_u(m) is the value of the m-th component of u on selected grid
    \param[out] trace_grad_u trace_grad_u(m) is the value of the m-th component of grad(u) on selected grid
   */
  template<class Dimension> template<class T>
  void VarOutputProblem_Dim<Dimension>
  ::ComputeInterpolationU_GradU(const Vector<T> & U0, int num_grid, bool compute_grad,
                                Vector<Vector<T> >& trace_u, Vector<Vector<T> >& trace_grad_u,
                                int nb_u) const
  {
    Vector<Vector<T> > Uvec;
    if (nb_u == 1)
      Uvec.SetData(1, const_cast<Vector<T>* >(&U0));
    else
      {
        int N = U0.GetM() / nb_u;
        Uvec.Reallocate(nb_u);
        for (int i = 0; i < nb_u; i++)
          Uvec(i).SetData(N, const_cast<T*>(&U0(i*N)));
      }
    
    Vector<Vector<T> > trace_u_diff, trace_grad_diff;
    ComputeInterpolationU_GradU(Uvec, this->all_points_display,
                                var_grid(num_grid), trace_u, trace_grad_u, 0,
                                Real_wp(0), trace_u_diff, trace_grad_diff, compute_grad);
    
    if (nb_u == 1)
      Uvec.Nullify();
    else
      for (int i = 0; i < nb_u; i++)
        Uvec(i).Nullify();
  }
  

  //! interpolation of U0 over a grid, result is placed in trace_u
  /*!
    \param[in] Unodal components of the solution on nodal points
    \param[in] grad_nodal gradient of the solution on nodal points
    \param[in] var_gr description of points where we want to know the solution
    \param[out] trace_u solution on asked points
  */
  template<class Dimension> template<class T>
  void VarOutputProblem_Dim<Dimension>::
  ComputeInterpolationNodalU(const Vector<Vector<Vector<T> > > & Unodal,
			     const Vector<Vector<Vector<T> > >& GradNodal,
			     const GridInterpolationFull<Dimension>& var_gr,
			     int nb_comp_all, Vector<Vector<T> >& trace_u) const
  {
    const GridInterpolation<Dimension>& var_interp = this->all_points_display;
    const IVect& list_points = var_gr.GetPointNumber();
    const Vector<IVect>& list_points_proc = var_gr.GetPointNumberAllProc();
    
    if (var_gr.GetNbGlobalPoints() <= 0)
      return;
    
    int nb_proc = var_problem.GetNbProcPerMode();

#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    int nnz = list_points.GetM();
    trace_u.Reallocate(nb_comp_all);
    
    // first step : we compute the values on U on each point of the grid for each processor
    // these values are stored in vector trace_vec
    Vector<T> trace_vec;
    this->ComputeInterpolationNodalUloc(Unodal, GradNodal, var_interp, trace_vec,
					nb_comp_all, list_points, nnz);    
    
    // second step : each processor sends its own data to the root processor
#ifdef SELDON_WITH_MPI
    Vector<int64_t> xtmp;
    MPI_Status status;
    if ((nnz > 0)&&(rank_proc != 0))
      {
        MpiSsend(var_problem.comm_group_mode, trace_vec, xtmp, trace_vec.GetM(), 0, 31);
      }
#endif
    
    Vector<T> trace_sol;
    if (rank_proc == 0)
      {
	// third step :
	// root processor retrieves the values computed on each processor
	// in the variable trace_sol
	int nb_pts_all = var_gr.GetNbGlobalPoints();
	for (int p = 0; p < nb_comp_all; p++)
	  {
	    trace_u(p).Reallocate(nb_pts_all);
	    trace_u(p).Zero();
	  }
	
	for (int i = 0; i < nb_proc; i++)
	  if (list_points_proc(i).GetM() > 0)
	    {
	      trace_sol.Reallocate(list_points_proc(i).GetM()*nb_comp_all);
	      
#ifdef SELDON_WITH_MPI
	      if (i != rank_proc)
                {
                  MpiRecv(var_problem.comm_group_mode, trace_sol, xtmp,
                          trace_sol.GetM(), i, 31, status);
                }
              else
#endif
		{
                  trace_sol = trace_vec;
		}
	      
	      // now we convert trace_sol in trace_u
	      for (int j = 0; j < list_points_proc(i).GetM(); j++)
		{
		  int nnz2 = list_points_proc(i).GetM();
		  int num_pt = list_points_proc(i)(j);
		  for (int p = 0; p < nb_comp_all; p++)
		    trace_u(p)(num_pt) = trace_sol(j + nnz2*p);
		}
	    }
      }
  }


  //! interpolation of U0 over a grid, result is placed in trace_vec
  template<class Dimension> template<class T>
  void VarOutputProblem_Dim<Dimension>
  ::ComputeInterpolationNodalUloc(const Vector<Vector<Vector<T> > > & Unodal,
				  const Vector<Vector<Vector<T> > >& GradNodal,
				  const GridInterpolation<Dimension> & var_interp,
				  Vector<T>& trace_vec, int nb_comp_all,
				  const IVect& list_points, int nnz) const
  {
    int nb = 0;
    trace_vec.Reallocate(nnz*nb_comp_all);
    trace_vec.Zero();
    
    typename Dimension::R_N point_loc;
    VectReal_wp phi;
    for (int i1 = 0; i1 < list_points.GetM(); i1++)
      {
	int i = list_points(i1);
	if (i >= 0)
	  {
	    // iquad is the element where the point i is
	    int iquad = var_interp.GetElementNumber(i);
	    
	    // local coordinates in the element
	    point_loc = var_interp.GetLocalCoordinate(i);
	    
	    // if the point has been found in the grid
	    if ((iquad >= 0)&&(iquad < var_problem.mesh.GetNbElt()))
	      {
		const ElementGeomReference<Dimension>& Fb = var_problem.GetShapeElement(iquad);
		
		// we compute the values of shape functions at point_loc
		Fb.ComputeValuesPhiNodalRef(point_loc, phi);
                
		for (int j = 0; j < Fb.GetNbPointsNodalElt(); j++)
		  {
		    for (int k = 0; k < Unodal(iquad).GetM(); k++)
		      {
			// we add contribution of u_i phi_i
			trace_vec(nb + k*nnz) += Unodal(iquad)(k)(j)*phi(j);
		      }
		    
		    int nb_comp = Unodal(iquad).GetM();
		    if (GradNodal.GetM() > 0)
		      for (int k = 0; k < GradNodal(iquad).GetM(); k++)
			{
			  trace_vec(nb + (k+nb_comp)*nnz) += GradNodal(iquad)(k)(j)*phi(j);
			}
		  }
	      }
            
	    nb++;
	  }
      }    
  }

    
  //! copies input parameters of another similar problem
  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>::CopyInputData(const VarOutputProblem_Base& var_base)
  {
    VarOutputProblem_Base::CopyInputData(var_base);
    
    const VarOutputProblem_Dim<Dimension>& var
      = static_cast<const VarOutputProblem_Dim<Dimension>& >(var_base);
    
    var_grid = var.var_grid;
    var_mesh_interp = var.var_mesh_interp;
  }


  //! writing the solution on "mesh-type" output files
  /*!
    \param[in] U0 vector solution
    \param[in] var_interp definition of the meshes where the solution has to be computed
    \param[in] grid parameters of output (name of the file, binary or ascii, etc)
    \param[in] num_file number of the output file to be written
    if 0, files are written on grid.file_diffracted_field
    if 1, files are written on grid.file_total_field
    if >= 2, files are written on grid.file_other(num_file-2)
  */
  template<class Dimension> template<class T>
  void VarOutputProblem_Dim<Dimension>::
  WriteOutputFileGen(const Vector<Vector<T> >& U0,
		     const Vector<MeshInterpolation<Dimension> >& var_interp,
		     const Vector<ParamOutputClass>& grid, int num_file, int nb_u) const
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    //int append = OutputTypeEnum::APPEND_NO;
    string suffix_mode;
    if ((var_boundary.ModesNotStored()) && (var_boundary.GetNbModes() > 1))
      {
        //append = OutputTypeEnum::APPEND_YES;
        //if (var_boundary.GetCurrentModeNumber() == var_boundary.GetNbModes()-1)
	//append = OutputTypeEnum::APPEND_FINAL;
	
	// here, we store the result for each mode, the final combination is left to the 
	// attention of the user, append is therefore useless in this function
	suffix_mode = "_M"+to_str(var_boundary.GetCurrentModeNumber());
      }
    
    Vector<Vector<Vector<T> > > valU_nodal, gradU_nodal;
    
    // computation of U0 and its gradient on nodal values 
    this->ComputeNodalUgradU(U0, valU_nodal, gradU_nodal, true, true);
    
    // name of vectors to export in vtk files
    IVect type_vtk_u, type_vtk_du;
    Vector<IVect> decomp_vtk_u, decomp_vtk_du;
    Vector<string> name_vtk_u, name_vtk_du;
    bool cplx_number = IsComplexNumber(T(0));
    
    // type of vtk output
    bool double_prec = true, ascii = true;
    if ((format_output_file == OutputTypeEnum::BINARY_MATLAB)
        || (format_output_file == OutputTypeEnum::BINARY_VTK))
      ascii = false;
    
    if (precision_output_file == OutputTypeEnum::SINGLE_PRECISION)
      double_prec = false;

    if (nb_u == -1)
      nb_u = var_problem.nb_unknowns;
    
    int nb_comp_all = var_problem.GetNbComponentsAll(nb_u);
    int nb_comp_grad = var_problem.GetNbComponentsGradientAll(nb_u);
    
    // loop on all outputs on meshes
    Dimension dim;
    for (int num_mesh = 0; num_mesh < var_interp.GetM(); num_mesh++)
      {
	const Vector<IVect> & Nodle = var_interp(num_mesh).NumberingArray();
	int Nvol = 0;
	for (int i = 0; i < Nodle.GetM(); i++)
	  for (int j = 0; j < Nodle(i).GetM(); j++)
	    Nvol = max(Nvol, Nodle(i)(j)+1);
	
	Vector<T> Uglob, dUglob;
	if (var_interp(num_mesh).IsVolumetric())
	  {
	    // output on volumetric meshes
	    Vector<T> Uloc, Udiv;
	    
	    this->GetVtkVolumeOutput(nb_comp_all, nb_comp_grad, Dimension::dim_N, cplx_number,
				     type_vtk_u, decomp_vtk_u, name_vtk_u,
				     type_vtk_du, decomp_vtk_du, name_vtk_du);
	    
	    Uglob.Reallocate(Nvol*nb_comp_all);
	    dUglob.Reallocate(Nvol*nb_comp_grad);
	    // loop on elements
	    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
	      {
		int nb_points_nodal = valU_nodal(i)(0).GetM();
                int type_elt = var_problem.mesh.Element(i).GetHybridType();
		
		Uloc.Reallocate(nb_points_nodal);
		for (int m = 0; m < nb_comp_all; m++)
		  {
		    // we get nodal values of U
		    for (int j = 0; j < nb_points_nodal; j++)
		      Uloc(j) = valU_nodal(i)(m)(j);
		    
		    // projection on the subdivided element
		    var_interp(num_mesh).ProjectVolume(Uloc, Udiv, type_elt);
		    
		    // values are stored in Uglob
		    for (int j = 0; j < Udiv.GetM(); j++)
		      Uglob(m*Nvol + Nodle(i)(j)) = Udiv(j);
		  }
		
		for (int m = 0; m < nb_comp_grad; m++)
		  {
		    // we get nodal values of gradient of U
		    for (int j = 0; j < nb_points_nodal; j++)
		      Uloc(j) = gradU_nodal(i)(m)(j);
		    
		    // projection on the subdivided element
		    var_interp(num_mesh).ProjectVolume(Uloc, Udiv, type_elt);
		    
		    // values are stored in dUglob
		    for (int j = 0; j < Udiv.GetM(); j++)
		      dUglob(m*Nvol + Nodle(i)(j)) = Udiv(j);
		  }
		
	      }    
	  }
	else if (var_interp(num_mesh).IsSurfacic())
	  {
	    // case of a surfacic mesh
	    // initializations of some variables
	    Vector<T> Uloc, dUloc, Uboundary, dUboundary, GradBoundary;
	    Vector<Vector<T> > Udiv(nb_comp_all), GradDiv(nb_comp_grad), val_En, val_Hn;
	    
	    Uglob.Reallocate(Nvol*var_problem.nb_components_en);
	    dUglob.Reallocate(Nvol*var_problem.nb_components_hn);
            
	    // loop on all boundaries of the surfacic mesh
            for (int i = 0; i < var_interp(num_mesh).GetNbBoundary(); i++)
	      {
		int num_elem = var_interp(num_mesh).GetElementNumberOfSurface(i);
		int num_loc = var_interp(num_mesh).GetLocalPositionOfSurface(i);
		int num_face = var_problem.mesh.Element(num_elem).numBoundary(num_loc);
                int type_surf = var_problem.mesh.Boundary(num_face).GetHybridType();
		
		// we get values of U and gradient on nodal points
		int nb_points_nodal = valU_nodal(num_elem)(0).GetM();
		Uloc.Reallocate(nb_points_nodal); 
		dUloc.Reallocate(nb_points_nodal);
		for (int m = 0; m < nb_comp_all; m++)
		  {
		    for (int j = 0; j < nb_points_nodal; j++)
		      Uloc(j) = valU_nodal(num_elem)(m)(j);

		    // extraction of values of u in the surface
		    this->ComputeValueNodalBoundary(Uloc, Uboundary, num_elem, num_loc);

		    // now projection on a subdivided boundary
		    var_interp(num_mesh).ProjectSurface(Uboundary, Udiv(m), type_surf);
		  }

		for (int m = 0; m < nb_comp_grad; m++)
		  {
		    for (int j = 0; j < nb_points_nodal; j++)                  
		      dUloc(j) = gradU_nodal(num_elem)(m)(j);
		
		    this->ComputeValueNodalBoundary(dUloc, GradBoundary, num_elem, num_loc);
		    
		    var_interp(num_mesh).ProjectSurface(GradBoundary, GradDiv(m), type_surf);
		  }
		
		// we get u, du/dn for Helmholtz, E \times n, H\times n for Maxwell, etc
		var_problem.ComputeEnHnNodal(Udiv, GradDiv, num_elem,
					     var_interp(num_mesh).PointsQuadrature(i),
					     var_interp(num_mesh).GetNormale(i), val_En, val_Hn);
                
		// we store the computed values
		for (int m = 0; m < var_problem.nb_components_en; m++)
		  for (int j = 0; j < Udiv(0).GetM(); j++)
		    Uglob(m*Nvol + Nodle(i)(j)) = val_En(m)(j);
		
		for (int m = 0; m < var_problem.nb_components_hn; m++)
		  for (int j = 0; j < Udiv(0).GetM(); j++)
		    dUglob(m*Nvol + Nodle(i)(j)) = val_Hn(m)(j);		
	      }
	    
	    this->GetVtkSurfaceOutput(val_En.GetM(), val_Hn.GetM(), Dimension::dim_N, cplx_number,
				      type_vtk_u, decomp_vtk_u, name_vtk_u,
				      type_vtk_du, decomp_vtk_du, name_vtk_du);
	    
	    nb_comp_all = var_problem.nb_components_en;
	    nb_comp_grad = var_problem.nb_components_hn;
	  }
	
	// result is now stored on Uglob and dUglob
	// sending these vectors to root processor
	int n1 = Uglob.GetM();
	
#ifdef SELDON_WITH_MPI
        Vector<int64_t> Uglob_tmp, dUglob_tmp;
	if (rank_proc != 0)
	  {
	    MPI_Send(&n1, 1, MPI_INTEGER, 0, 7, var_problem.comm_group_mode);
            
	    if (n1 > 0)
              MpiSend(var_problem.comm_group_mode, Uglob, Uglob_tmp, Uglob.GetM(), 0, 8);
	    
	    n1 = dUglob.GetM();
	    MPI_Send(&n1, 1, MPI_INTEGER, 0, 9, var_problem.comm_group_mode);
            
	    if (n1 > 0)
              MpiSend(var_problem.comm_group_mode, dUglob, dUglob_tmp, dUglob.GetM(), 0, 12);
          }
#endif
	
	if (rank_proc == 0)
	  {
	    // storing Uglob and dUglob for each processor
	    Vector<Vector<T> > Uproc(nb_proc), dUproc(nb_proc);
	    IVect nb_dof_proc(nb_proc); Nvol = 0;
	    IVect offset_dof(nb_proc+1); offset_dof.Fill(0);
	    // receiving values
	    for (int i = 0; i < nb_proc; i++)
	      {
		if (i != rank_proc)
		  {
#ifdef SELDON_WITH_MPI
                    MPI_Status status;
                    MPI_Recv(&n1, 1, MPI_INTEGER, i, 7, var_problem.comm_group_mode, &status);
		    if (n1 > 0)
                      {
                        Uproc(i).Reallocate(n1);
                        MpiRecv(var_problem.comm_group_mode, Uproc(i), Uglob_tmp, n1, i, 8, status);
                      }
                    else
                      Uproc(i).Clear();
                    
		    nb_dof_proc(i) = n1/nb_comp_all;
		    Nvol += nb_dof_proc(i);
		    offset_dof(i+1) = Nvol;
		    
                    MPI_Recv(&n1, 1, MPI_INTEGER, i, 9, var_problem.comm_group_mode, &status);
		    
                    if (n1 > 0)
                      {
                        dUproc(i).Reallocate(n1);
                        MpiRecv(var_problem.comm_group_mode, dUproc(i), dUglob_tmp, n1, i, 12, status);
                      }
                    else
                      dUproc(i).Clear();
#endif
		  }
		else
		  {
		    Uproc(i) = Uglob; dUproc(i) = dUglob;
		    n1 = Uglob.GetM(); nb_dof_proc(i) = n1/nb_comp_all;
		    Nvol += nb_dof_proc(i);
                    offset_dof(i+1) = Nvol;
		    Uglob.Clear(); dUglob.Clear();
		  }
	      }
	    
	    // writing the result on bb files
	    // depending on the components of U and dU, the procedure is the following :
	    // 1- extracting each component on the vector output_vector
	    // 2- adapting the name of the output file
	    //      (_U0.bb, _U1.bb, _dU0.bb) according to the component
	    // 3- writing the file on the disk
	    // 4- make a symbolic link towards the mesh
	    //     (to avoid several copies of the mesh on the disk)
	    string file_name = grid(num_mesh).GetFileName(num_file);
	    string root = GetBaseString(file_name) + suffix_mode;
	    int component = grid(num_mesh).GetComponent();
            string name_mesh;	 
	    if ((format_output_file == OutputTypeEnum::ASCII_VTK)
                || (format_output_file == OutputTypeEnum::BINARY_VTK))
	      {
		name_mesh = GetBaseString(file_name) + suffix_mode + ".vtk";
                
		if (this->move_points_on_output_mesh)
                  {
#ifdef MONTJOIE_WITH_THREE_DIM
                    Mesh<Dimension> mesh_surf;
                    Mesh<Dimension3> mesh3d;
                    mesh_surf.Read(DOSSIER_output + var_interp(num_mesh).GetFileName());
                    
                    int m = this->component_used_to_translate;
                    if (Dimension::dim_N == 2)
                      {
                        mesh3d.ReallocateVertices(mesh_surf.GetNbVertices());
                        mesh3d.ReallocateBoundariesRef(mesh_surf.GetNbElt());
                        for (int i = 0; i < mesh_surf.GetNbVertices(); i++)
                          mesh3d.Vertex(i).Init(mesh_surf.Vertex(i)(0),
                                                mesh_surf.Vertex(i)(1), 0.0);
                        
                        for (int i = 0; i < mesh_surf.GetNbElt(); i++)
                          {
                            IVect num(mesh_surf.Element(i).GetNbVertices());
                            for (int j = 0; j < num.GetM(); j++)
                              num(j) = mesh_surf.Element(i).numVertex(j);
                            
                            int ref = 1;
                            mesh3d.BoundaryRef(i).Init(num, ref);
                          }
                        
                        for (int n = 0; n < nb_proc; n++)
                          for (int i = 0; i < nb_dof_proc(n); i++)
                            mesh3d.Vertex(i+offset_dof(n))(2)
                              = this->coefficient_used_to_translate
                              *realpart(Uproc(n)(i+m*nb_dof_proc(n)));
                        
                        mesh3d.FindConnectivity();
                        mesh3d.Write(DOSSIER_output + name_mesh, double_prec, ascii);
                      }
                    else
                      {
                        for (int n = 0; n < nb_proc; n++)
                          for (int i = 0; i < nb_dof_proc(n); i++)
                            mesh_surf.Vertex(i+offset_dof(n))(2)
                              += this->coefficient_used_to_translate
                              *realpart(Uproc(n)(i+m*nb_dof_proc(n)));
                        
                        mesh_surf.Write(DOSSIER_output + name_mesh, double_prec, ascii);
                      }
#endif
                  }
                else
                  CopyFile(var_interp(num_mesh).GetFileName(), name_mesh);
		
                name_mesh = DOSSIER_output + name_mesh;
                ofstream file_out(name_mesh.data(), ios::app);
		file_out << "POINT_DATA " << Nvol << '\n';
		
		Vector<Vector<T> > output_vector(nb_comp_all);
		
		if ((component == -1) || (component == -2))
		  {
		    // components of u are written
		    for (int m = 0; m < nb_comp_all; m++)
		      {
			output_vector(m).Reallocate(Nvol);
			for (int n = 0; n < nb_proc; n++)
			  for (int i = 0; i < nb_dof_proc(n); i++)
			    output_vector(m)(i+offset_dof(n)) = Uproc(n)(i+m*nb_dof_proc(n));
		      }
		    
		    WriteVtk(output_vector, type_vtk_u, decomp_vtk_u, name_vtk_u,
			     file_out, double_prec, ascii);
		  }
		
		output_vector.Reallocate(nb_comp_grad);
		if (component == -1)
		  {
		    for (int m = 0; m < nb_comp_grad; m++)
                      {
			output_vector(m).Reallocate(Nvol);
                        for (int n = 0; n < nb_proc; n++)
                          for (int i = 0; i < nb_dof_proc(n); i++)
                            output_vector(m)(i+offset_dof(n)) = dUproc(n)(i+m*nb_dof_proc(n));
		      }
		 		
		    WriteVtk(output_vector, type_vtk_du, decomp_vtk_du, name_vtk_du,
			     file_out, double_prec, ascii);
		  }

		if ((component >= 0) && (component < nb_comp_all))
		  {
		    Vector<T> output_vec(Nvol);
		    // one component of u is written
		    for (int n = 0; n < nb_proc; n++)
		      for (int i = 0; i < nb_dof_proc(n); i++)
			output_vec(i+offset_dof(n)) = Uproc(n)(i+component*nb_dof_proc(n));
		    
		    WriteVtk(output_vec, "u", file_out, double_prec, ascii);
		  }
		
		if ((component >= nb_comp_all) && (component < nb_comp_grad))
		  {
		    Vector<T> output_vec(Nvol);
		    // one component of dU is written
		    component -= nb_comp_all;
		    for (int n = 0; n < nb_proc; n++)
		      for (int i = 0; i < nb_dof_proc(n); i++)
			output_vec(i+offset_dof(n)) = dUproc(n)(i+component*nb_dof_proc(n));
		    
		    WriteVtk(output_vec, "u", file_out, double_prec, ascii);
		  }		
		
		file_out.close();
	      }
	    else
	      {
		Vector<T> output_vector(Nvol);
		name_mesh = GetBaseString(file_name) + suffix_mode + ".mesh";
                
		if (this->move_points_on_output_mesh)
                  {
#ifdef MONTJOIE_WITH_THREE_DIM
                    Mesh<Dimension> mesh_surf;
                    Mesh<Dimension3> mesh3d;
                    mesh_surf.Read(DOSSIER_output + var_interp(num_mesh).GetFileName());
                    
                    int m = this->component_used_to_translate;
                    if (Dimension::dim_N == 2)
                      {
                        mesh3d.ReallocateVertices(mesh_surf.GetNbVertices());
                        mesh3d.ReallocateBoundariesRef(mesh_surf.GetNbElt());
                        for (int i = 0; i < mesh_surf.GetNbVertices(); i++)
                          mesh3d.Vertex(i).Init(mesh_surf.Vertex(i)(0),
                                                mesh_surf.Vertex(i)(1), 0.0);
                        
                        for (int i = 0; i < mesh_surf.GetNbElt(); i++)
                          {
                            IVect num(mesh_surf.Element(i).GetNbVertices());
                            for (int j = 0; j < num.GetM(); j++)
                              num(j) = mesh_surf.Element(i).numVertex(j);
                            
                            int ref = 1;
                            mesh3d.BoundaryRef(i).Init(num, ref);
                          }
                        
                        for (int n = 0; n < nb_proc; n++)
                          for (int i = 0; i < nb_dof_proc(n); i++)
                            mesh3d.Vertex(i+offset_dof(n))(2)
                              = this->coefficient_used_to_translate
                              *realpart(Uproc(n)(i+m*nb_dof_proc(n)));
                        
                        mesh3d.FindConnectivity();
                        mesh3d.Write(DOSSIER_output+name_mesh, double_prec, ascii);
                      }
                    else
                      {
                        for (int n = 0; n < nb_proc; n++)
                          for (int i = 0; i < nb_dof_proc(n); i++)
                            mesh_surf.Vertex(i+offset_dof(n))(2)
                              += this->coefficient_used_to_translate
                              *realpart(Uproc(n)(i+m*nb_dof_proc(n)));
                        
                        mesh_surf.Write(DOSSIER_output+name_mesh, double_prec, ascii);
                      }
#endif
                  }
                else
                  MakeLink(DOSSIER_output + var_interp(num_mesh).GetFileName(),
			   DOSSIER_output + name_mesh);
		
		if ((component == -1) || (component == -2))
		  {
		    // components of u are written
		    for (int m = 0; m < nb_comp_all; m++)
		      {
			for (int n = 0; n < nb_proc; n++)
			  for (int i = 0; i < nb_dof_proc(n); i++)
			    output_vector(i+offset_dof(n)) = Uproc(n)(i+m*nb_dof_proc(n));
			
			string name = root + "_U" + to_str(m) + ".bb";
			if (root.find('.',0) != string::npos)
			  name = GetBaseString(root) + "_U" + to_str(m) +"."
			    + GetExtension(root) + ".bb";
			
			WriteMedit(output_vector, DOSSIER_output+name, dim, double_prec);
			name = GetBaseString(name) + ".mesh";
			MakeLink(DOSSIER_output+name_mesh, DOSSIER_output+name);
		      }
		  }
		
                if (component == -1)
                  {
                    // components of dU are written
                    for (int m = 0; m < nb_comp_grad; m++)
                      {
                        for (int n = 0; n < nb_proc; n++)
                          for (int i = 0; i < nb_dof_proc(n); i++)
                            output_vector(i+offset_dof(n)) = dUproc(n)(i+m*nb_dof_proc(n));
			
			string name = root + "_dU" + to_str(m) + ".bb";
			if (root.find('.',0) != string::npos)
			  name = GetBaseString(root) + "_dU" + to_str(m) +"."
			    + GetExtension(root) + ".bb";
			
			WriteMedit(output_vector, DOSSIER_output+name, dim, double_prec);
			name = GetBaseString(name) + ".mesh";
			MakeLink(DOSSIER_output+name_mesh, DOSSIER_output+name);
		      }
                  }
		
		if ((component >= 0) && (component < nb_comp_all))
		  {
		    // one component of u is written
                    for (int n = 0; n < nb_proc; n++)
		      for (int i = 0; i < nb_dof_proc(n); i++)
			output_vector(i+offset_dof(n)) = Uproc(n)
			  (i+component*nb_dof_proc(n));
		    
		    string name = root + ".bb";
		    WriteMedit(output_vector, DOSSIER_output+name, dim, double_prec);
		  }
		
		if ((component >= nb_comp_all) && (component < nb_comp_grad))
		  {
		    // one component of dU is written
		    component -= nb_comp_all;
		    for (int n = 0; n < nb_proc; n++)
		      for (int i = 0; i < nb_dof_proc(n); i++)
			output_vector(i+offset_dof(n))
			  = dUproc(n)(i+component*nb_dof_proc(n));
		    
		    string name = root + ".bb";
		    WriteMedit(output_vector, DOSSIER_output+name, dim, double_prec);
		  }
	      }
	  }
      }
    
    if (output_solution_file)
      {
        if (nb_proc > 1)
          {
            cout << "not implemented " << endl;
            abort();
          }
        
        int nb_elt = valU_nodal.GetM();
        Vector<Vector<T> > nu(nb_elt);
        for (int i = 0; i < nb_elt; i++)
          nu(i).Reallocate(valU_nodal(i)(0).GetM());

        int component = output_solution_param.GetComponent();
        if ((component == -1) || (component == -2))
          {            
            for (int m = 0; m < nb_comp_all; m++)
              {
                for (int i = 0; i < nb_elt; i++)
                  for (int j = 0; j < nu(i).GetM(); j++)
                    nu(i)(j) = valU_nodal(i)(m)(j);
                
                string root = output_solution_param.GetFileName(num_file);
                string name = root + "_U" + to_str(m) + ".elb";
                WriteMeshData(nu, DOSSIER_output+name);
              }
            
            if (component == -1)
              for (int m = 0; m < nb_comp_grad; m++)
                {
                  for (int i = 0; i < nb_elt; i++)
                    for (int j = 0; j < nu(i).GetM(); j++)
                      nu(i)(j) = gradU_nodal(i)(m)(j);
                  
                  string root = output_solution_param.GetFileName(num_file);
                  string name = root + "_dU" + to_str(m) + ".elb";
                  WriteMeshData(nu, DOSSIER_output+name);
                }
            
          }
        else if (component < nb_comp_all)
          {
            int m = component;
            for (int i = 0; i < nb_elt; i++)
              for (int j = 0; j < nu(i).GetM(); j++)
                nu(i)(j) = valU_nodal(i)(m)(j);
            
            string root = output_solution_param.GetFileName(num_file);
            string name = root + ".elb";
            WriteMeshData(nu, DOSSIER_output+name);
          }
        else
          {
            int m = component - nb_comp_all;
            for (int i = 0; i < nb_elt; i++)
              for (int j = 0; j < nu(i).GetM(); j++)
                nu(i)(j) = gradU_nodal(i)(m)(j);
            
            string root = output_solution_param.GetFileName(num_file);
            string name = root + ".elb";
            WriteMeshData(nu, DOSSIER_output+name);
          }        
      }
  }
  
  
  //! writing of U0 on output files specified by the user
  /*!
    \param[in] U0 solution vector
    \param[in] var_interp description of the interpolation points
    \param[in] num_file file number
    \param[in] file_name name of the output file
    \param[in] add_total_field if equal to 1, the incident field is added to obtain the total field
    \param[in] file_name_diff name of the output file for the diffracted field
    \param[in] append values to append on existing file ?
  */
  template<class Dimension> template<class T, class TypeInterpolation>
  void VarOutputProblem_Dim<Dimension>::
  WriteOutputFileGen(const Vector<Vector<T> > & U0, const TypeInterpolation & var_interp,
		     int num_file, const string & file_name, 
		     int add_total_field, const string& file_name_diff,
		     const ParamOutputClass& param, int append) const
  {  
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    int component = param.GetComponent();
    int type = param.GetNature();
    string ext = output_file_extension;
    int N = 0, ngrid = var_interp.GetM();
    for (int i = 0; i < ngrid; i++)
      if (var_interp(i).GetOutputType() == type)
	N += var_interp(i).GetNbGlobalPoints();
    
    bool vtk_output = false;
    if ((format_output_file == OutputTypeEnum::BINARY_VTK)
        || (format_output_file == OutputTypeEnum::ASCII_VTK))
      vtk_output = true;
    
    Vector<T> output_vector;
    if (rank_proc == 0)
      {
	output_vector.Reallocate(N);
	output_vector.Zero();
      }
    
    Vector<Vector<Vector<T> > > trace_u(ngrid), trace_u_diff(ngrid);
    Vector<Vector<Vector<T> > > trace_grad(ngrid), trace_grad_diff(ngrid);
    Vector<Vector<T> > output_u, output_du, output_u_diff, output_du_diff;

    int nb_modes = var_boundary.GetNbModes();
    //int nb_proc = var_problem.GetNbProcPerMode();
    if (var_boundary.ModesNotStored())
      nb_modes = 1;
    
    int nb_vec = U0.GetM()/nb_modes;
    int nb_comp_all = var_problem.GetNbComponentsAll(nb_vec);
    int nb_comp_grad = var_problem.GetNbComponentsGradientAll(nb_vec);
    
    bool compute_grad = false;
    if ((component == -1) || (component >= nb_comp_all))
      compute_grad = true;

    // on ecrit toutes les composantes de u et grad u
    for (int i = 0; i < ngrid; i++)
      if (var_interp(i).GetOutputType() == type)
	{
	  this->ComputeInterpolationU_GradU(U0, this->all_points_display,
					    var_interp(i), trace_u(i), trace_grad(i),
					    add_total_field, param.GetCurrentTime(),
					    trace_u_diff(i), trace_grad_diff(i), compute_grad);
	  
	}
    
    if ((component == -1) || (component == -2))
      {
        if (component == -1)
          param.InitResult(output_vector, nb_comp_all+nb_comp_grad);
        else
          param.InitResult(output_vector, nb_comp_all);
        
	if (rank_proc == 0)
	  {
	    output_u.Reallocate(nb_comp_all);
	    output_u_diff.Reallocate(nb_comp_all);
            
	    for (int m = 0; m < nb_comp_all; m++)
	      {
                string root = GetBaseString(file_name);
		int nb = 0;
		for (int i = 0; i < ngrid; i++)
		  if (var_interp(i).GetOutputType() == type)
		    for (int j = 0; j < var_interp(i).GetNbGlobalPoints(); j++)
		      output_vector(nb++) = trace_u(i)(m)(j);
		
		if (vtk_output)
		  output_u(m) = output_vector;
		else
		  {
		    string name = root + "_U" + to_str(m) + ext;
                    WriteOutputVector(var_interp, output_vector, output_u,
				      output_du, num_file, name, type, m, param, append);
		  }
		
                if (add_total_field != 0)
                  {
                    root = GetBaseString(file_name_diff);
                    nb = 0;
                    for (int i = 0; i < ngrid; i++)
                      if (var_interp(i).GetOutputType() == type)
                        for (int j = 0; j < var_interp(i).GetNbGlobalPoints(); j++)
                          output_vector(nb++) = trace_u_diff(i)(m)(j);
                    
                    if (vtk_output)
                      output_u_diff(m) = output_vector;
                    else
                      {
                        string name = root + "_U" + to_str(m) + ext;
                        WriteOutputVector(var_interp, output_vector, output_u_diff,
                                          output_du_diff, 0, name, type, m, param, append);
                      }                    
                  }
	      }   
	    
            int off_m = nb_comp_all;
            if (component == -1)
	      {
		output_du.Reallocate(nb_comp_grad);
		for (int m = 0; m < nb_comp_grad; m++)
		  {
                    string root = GetBaseString(file_name);
		    int nb = 0;
		    for (int i = 0; i < ngrid; i++)
		      if (var_interp(i).GetOutputType() == type)
			for (int j = 0; j < var_interp(i).GetNbGlobalPoints(); j++)
			  output_vector(nb++) = trace_grad(i)(m)(j);
		    
		    if (vtk_output)
		      output_du(m) = output_vector;
		    else
		      {
			string name = root + "_dU" + to_str(m) + ext;
			WriteOutputVector(var_interp, output_vector, output_u,
					  output_du, num_file, name, type,
					  m+off_m, param, append);
		      }

                    if (add_total_field != 0)
                      {
                        root = GetBaseString(file_name_diff);
                        nb = 0;
                        for (int i = 0; i < ngrid; i++)
                          if (var_interp(i).GetOutputType() == type)
                            for (int j = 0; j < var_interp(i).GetNbGlobalPoints(); j++)
                              output_vector(nb++) = trace_grad_diff(i)(m)(j);
                        
                        if (vtk_output)
                          output_du_diff(m) = output_vector;
                        else
                          {
                            string name = root + "_dU" + to_str(m) + ext;
                            WriteOutputVector(var_interp, output_vector, output_u_diff,
                                              output_du_diff, 0, name, type,
                                              m+off_m, param, append);
                          }
                      }
                  }
              }
	    
            if (vtk_output)
	      {
		string name = GetBaseString(file_name) + ext;
		WriteOutputVector(var_interp, output_vector, output_u,
				  output_du, num_file, name, type,
				  -1, param, append);
                
                if (add_total_field != 0)
                  {
                    name = GetBaseString(file_name_diff) + ext;
                    WriteOutputVector(var_interp, output_vector, output_u_diff,
                                      output_du_diff, 0, name, type,
                                      -1, param, append);                  
                  }
	      }
	  }
      }
    else if (component < nb_comp_all)
      {
        param.InitResult(output_vector, 1);
        
	if (rank_proc == 0)
	  {
	    int nb = 0;
	    for (int i = 0; i < ngrid; i++)
	      if (var_interp(i).GetOutputType() == type)
		for (int j = 0; j < var_interp(i).GetNbGlobalPoints(); j++)
		  output_vector(nb++) = trace_u(i)(component)(j);
	    
	    WriteOutputVector(var_interp, output_vector, output_u, output_du,
			      num_file, file_name, type, 0, param, append);
            
            if (add_total_field != 0)
              {
                nb = 0;
                for (int i = 0; i < ngrid; i++)
                  if (var_interp(i).GetOutputType() == type)
                    for (int j = 0; j < var_interp(i).GetNbGlobalPoints(); j++)
                      output_vector(nb++) = trace_u_diff(i)(component)(j);
                
                WriteOutputVector(var_interp, output_vector, output_u_diff, output_du_diff,
                                  0, file_name_diff, type, 0, param, append);
              }
          }
      }
    else
      {
        param.InitResult(output_vector, 1);
        
        if (rank_proc == 0)
          {
            component -= nb_comp_all;
            int nb = 0;
            for (int i = 0; i < ngrid; i++)
              if (var_interp(i).GetOutputType() == type)
                for (int j = 0; j < var_interp(i).GetNbGlobalPoints(); j++)
                  output_vector(nb++) = trace_grad(i)(component)(j);
            
            WriteOutputVector(var_interp, output_vector, output_u,
			      output_du, num_file, file_name, type, 0, param, append);
            
            if (add_total_field != 0)
              {
                nb = 0;
                for (int i = 0; i < ngrid; i++)
                  if (var_interp(i).GetOutputType() == type)
                    for (int j = 0; j < var_interp(i).GetNbGlobalPoints(); j++)
                      output_vector(nb++) = trace_grad_diff(i)(component)(j);
                
                WriteOutputVector(var_interp, output_vector, output_u_diff,
                                  output_du_diff, 0, file_name_diff, type, 0, param, append);
              }
          }
      }    
  }
  

  //! computes far field and writes it on output files
  template<class Dimension> template<class T>
  void VarOutputProblem_Dim<Dimension>::
  WriteOutputFileGen(const Vector<Vector<T> >& U0,
		     const VarComputationRCS_Base<Dimension>& var_rcs,
		     const ParamOutputClass& param, int num_file) const
  {
    int N = var_rcs.GetNbPointsOutside();
    if (N <= 0)
      return;
    
    Vector<T> trace_En, trace_Hn;
    var_problem.ComputeEnHnOnBoundary(var_rcs.GetInterpolationMesh(), U0(0), trace_En, trace_Hn, false);
    
    Vector<T> trace_p(N);
    Vector<T> val_u(1);
    for (int i = 0; i < N; i++)
      {
        var_rcs.ComputeIntegralRepresentation(trace_En, trace_Hn, var_rcs.GetInterpolationMesh(),
                                              var_rcs.GetOutsidePoint(i), val_u);
        
        trace_p(i) = val_u(0);
      }
    
    // assembling u
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    if (nb_proc > 1)
      {
        Vector<T> trace_ploc(trace_p);
        Vector<int64_t> trace_tmp;
        MpiReduce(var_problem.comm_group_mode, trace_ploc, trace_tmp,
                  trace_p, trace_p.GetM(), MPI_SUM, 0);                
      }
    
    if (rank_proc == 0)
      {
#endif
        GridInterpolationFull<Dimension> grid;
        grid.SetPointsOutput(var_rcs.GetOutsidePoint());
        
        string root = GetBaseString(param.GetFileName(num_file));
        string name = root + ".dat";
        bool ascii = true;
        if ((format_output_file == OutputTypeEnum::BINARY_MATLAB)
            || (format_output_file == OutputTypeEnum::BINARY_VTK))
          ascii = false;
        
        WriteMatlab(trace_p, grid, DOSSIER_output+name, precision_output_file, ascii);

#ifdef SELDON_WITH_MPI
      }
#endif
    
  }
  

  //! Writes the solution on output files by interpolation with nodal points
  template<class Dimension> template<class T, class TypeInterpolation>
  void VarOutputProblem_Dim<Dimension>::
  WriteOutputFileNodalGen(const Vector<Vector<T> > & U0, const TypeInterpolation & var_interp,
			  const Vector<ParamOutputClass>& grid, int num_file) const
  {
    if (grid.GetM() <= 0)
      return;

#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    int component = grid(0).GetComponent();
    int nb_u = var_problem.nb_unknowns_scal;
    int nb_comp = var_problem.GetNbComponentsAll(nb_u);
    int nb_comp_grad = var_problem.GetNbComponentsGradientAll(nb_u);
    int nb_comp_hess = var_problem.GetNbComponentsHessianAll(nb_u);
    bool u_comp = true, grad_comp = false, hess_comp = false;
    Vector<string> suffix(nb_comp);
    for (int k = 0; k < nb_comp; k++)
      suffix(k) = "_U" + to_str(k) + ".dat";

    if (component <= -4)
      {
	grad_comp = true;
	suffix.Resize(nb_comp + nb_comp_grad);
	for (int k = 0; k < nb_comp_grad; k++)
	  suffix(nb_comp+k) = "_dU" + to_str(k) + ".dat";
      }
    
    if (component <= -5)
      {
	hess_comp = true;
	suffix.Resize(nb_comp + nb_comp_grad + nb_comp_hess);
	for (int k = 0; k < nb_comp_hess; k++)
	  suffix(nb_comp+nb_comp_grad+k) = "_ddU" + to_str(k) + ".dat";
      }
    
    int nb_comp_all = suffix.GetM();
    
    // computation of U0 on nodal points
    Vector<Vector<Vector<T> > > valU_nodal, gradU_nodal;
    this->ComputeNodalUgradU(U0, valU_nodal, gradU_nodal, u_comp,
			     grad_comp, hess_comp, true, nb_u);
    
    int ngrid = var_interp.GetM();
    Vector<Vector<Vector<T> > > trace_u(ngrid);
    
    for (int n = 0; n < grid.GetM(); n++)
      {
	int type = grid(n).GetNature();
	int nb_true = 0, N = 0;
	for (int i = 0; i < ngrid; i++)
	  if (var_interp(i).GetOutputType() == type)
	    {
	      this->ComputeInterpolationNodalU(valU_nodal, gradU_nodal, var_interp(i),
					       nb_comp_all, trace_u(i));
	      nb_true++;
	      N += var_interp(i).GetNbGlobalPoints();
	    }
	
	if (rank_proc != 0)
	  continue;
	
	bool ascii = false;
	string file_name = grid(n).GetFileName(num_file);
	string root = GetBaseString(file_name);
	
	Vector<T> output_vector(N);
	for (int m = 0; m < nb_comp_all; m++)
	  {
	    int nb = 0;
	    for (int i = 0; i < ngrid; i++)
	      if (var_interp(i).GetOutputType() == type)
		for (int j = 0; j < var_interp(i).GetNbGlobalPoints(); j++)
		  output_vector(nb++) = trace_u(i)(m)(j);
	    
	    WriteMatlab(output_vector, nb_true, var_interp, type,
			DOSSIER_output+root+suffix(m), precision_output_file, ascii);
	  }
      }
  }
  

  //! writing of the interpolation of U0 over a grid in a file,
  //! it can be also a component of gradient of U0
  /*!
    \param[in] U0 solution vector
    \param[in] var_interp description of the interpolation points
    \param[in] grid output parameters
    \param[in] num_file identifying number for output file
  */
  template<class Dimension> template<class T, class TypeInterpolation>
  void VarOutputProblem_Dim<Dimension>::
  WriteOutputFileGen(const Vector<Vector<T> > & U0, const TypeInterpolation & var_interp,
		     const Vector<ParamOutputClass>& grid, int num_file) const
  {
    int append = OutputTypeEnum::APPEND_NO;
    if ((var_boundary.ModesNotStored()) && (var_boundary.GetNbModes() > 1))
      {
        append = OutputTypeEnum::APPEND_YES;
        if (this->last_output_solution)
          append = OutputTypeEnum::APPEND_FINAL;
      }
    
    if (grid.GetM() > 0)
      {
	if (grid(0).GetComponent() <= -3)
	  {
	    // in this case, U is interpolated on nodal points before interpolation
	    WriteOutputFileNodalGen(U0, var_interp, grid, num_file);
	    return;
	  }
      }

    for (int i = 0; i < grid.GetM(); i++)
      WriteOutputFileGen(U0, var_interp, num_file, grid(i).GetFileName(num_file),
                         0, grid(i).GetFileName(0), grid(i), append);
  }
  

  //! writing the output indexed num_file in the matching output files
  /*!
    \param[in] U0 components of the solution on dofs
    \param[in] num_file output number
  */
  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>::
  WriteOutputFileReshaped(const Vector<VectReal_wp>& U0, int num_file) const
  {
    if (this->output_grid_param.GetM() !=0)
      this->WriteOutputFileGen(U0, this->var_grid,
			       this->output_grid_param, num_file);
    
    if (this->output_mesh_param.GetM() !=0)
      this->WriteOutputFileGen(U0, this->var_mesh_interp,
			       this->output_mesh_param, num_file);
    
    if (this->output_rcs_param.GetNbPointsOutside() > 0)
      this->WriteOutputFileGen(U0, this->output_rcs_param, 
			       this->output_points_outside_param, num_file);

    if ((this->write_solution_on_quad_points) && (num_file == 1))
      this->WriteUquadrature(U0, this->file_name_solution_quad, this->file_name_weights_quad);    
  }


  //! writing the output indexed num_file in the matching output files
  /*!
    \param[in] U0 components of the solution on dofs
    \param[in] num_file output number
  */
  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>::
  WriteOutputFileReshaped(const Vector<VectComplex_wp>& U0, int num_file) const
  {
    if (this->output_grid_param.GetM() !=0)
      this->WriteOutputFileGen(U0, this->var_grid,
			       this->output_grid_param, num_file);
    
    if (this->output_mesh_param.GetM() !=0)
      this->WriteOutputFileGen(U0, this->var_mesh_interp,
			       this->output_mesh_param, num_file);

    // total field is not written for scattering problems
    // since integral representation discards the incident field
    bool compute_outside_field = true;
    if ((num_file == 1) && ((var_source_base.GetSourceType(0) == var_source_base.SRC_TOTAL_FIELD)
			    || (var_source_base.GetSourceType(0) == var_source_base.SRC_DIFFRACTED_FIELD)))
      compute_outside_field = false;
    
    if ((this->output_rcs_param.GetNbPointsOutside() > 0) && compute_outside_field)
      this->WriteOutputFileGen(U0, this->output_rcs_param, 
			       this->output_points_outside_param, num_file);

    if ((this->write_solution_on_quad_points) && (num_file == 1))
      this->WriteUquadrature(U0, this->file_name_solution_quad, this->file_name_weights_quad);    
  }

  
  //! writing result of interpolation on a file
  /*!
    \param[in] var_interp description of grid where the result is computed
    \param[in] output_vector result computed on the grid
    \param[in] file_name name of the ouput file
    \param[in] type type of the grid
    \param[in] append overwrite file or append at the end of the file ?
  */
  template<class Dimension> template<class T, class VarInterpolation>
  void VarOutputProblem_Dim<Dimension>::
  WriteOutputVector(const VarInterpolation & var_interp, Vector<T> & output_vector,
		    Vector<Vector<T> >& output_u, Vector<Vector<T> >& output_du,
		    int num_file, const string & file_name,
		    int type, int m, const ParamOutputClass& param, int append) const
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    if ((rank_proc == 0)||(var_problem.print_level >= 10))
      if (var_problem.print_level >= 4)
	cout << rank_proc << " Writing" << endl;
    
    int ngrid = var_interp.GetM(), nb_true = 0;
    for (int i = 0; i < ngrid; i++)
      if (var_interp(i).GetOutputType() == type)
	nb_true++;
    
    if (nb_true <= 0)
      return;
    
    bool ascii = true;
    if ((format_output_file == OutputTypeEnum::BINARY_VTK)
        || (format_output_file == OutputTypeEnum::BINARY_MATLAB))
      ascii = false;

    bool complex_number = IsComplexNumber(T(0));
    bool double_prec = false;
    if (precision_output_file == OutputTypeEnum::DOUBLE_PRECISION)
      double_prec = true;
    
    if ((append == OutputTypeEnum::APPEND_YES) || (append == OutputTypeEnum::APPEND_FINAL))
      {
	if (m == -1)
	  {
	    for (int m = 0; m < output_u.GetM(); m++)
	      param.UpdateResult(num_file, m, output_u(m));
	    
	    int off_m = output_u.GetM();
	    for (int m = 0; m < output_du.GetM(); m++)
	      param.UpdateResult(num_file, off_m+m, output_du(m));
	  }
	else
	  param.UpdateResult(num_file, m, output_vector);
      }

    if ((append == OutputTypeEnum::APPEND_NO) || (append == OutputTypeEnum::APPEND_FINAL))
      {        
        if (rank_proc != 0)
          return;
        
        if ((format_output_file == OutputTypeEnum::BINARY_MATLAB)
            || (format_output_file == OutputTypeEnum::ASCII_MATLAB))
	  {
	    WriteMatlab(output_vector, nb_true,
			var_interp, type, DOSSIER_output+file_name, precision_output_file, ascii);
	  }
        else if ((format_output_file == OutputTypeEnum::ASCII_VTK)
                 || (format_output_file == OutputTypeEnum::BINARY_VTK))
          {
            if (nb_true > 1)
              {
                cout << "not implemented " << endl;
                abort();
              }
            
            if (m == -1)
	      {
		int nb_u = output_u.GetM();
		int nb_du = output_du.GetM();
		IVect type_vtk_u, type_vtk_du;
		Vector<IVect> decomp_vtk_u, decomp_vtk_du;
		Vector<string> name_vtk_u, name_vtk_du;
		
		this->GetVtkVolumeOutput(nb_u, nb_du, Dimension::dim_N, complex_number,
					 type_vtk_u, decomp_vtk_u, name_vtk_u,
					 type_vtk_du, decomp_vtk_du, name_vtk_du);
		
                string file_name_tmp = DOSSIER_output + file_name;
		ofstream file_out(file_name_tmp.data());
		if (!file_out.is_open())
		  {
		    cout << "Unable to open " << file_name << endl;
		    abort();
		  }

		for (int i = 0; i < ngrid; i++)
		  if (var_interp(i).GetOutputType() == type)
		    {
		      var_interp(i).WriteVtk(file_out, double_prec, ascii);
		      WriteVtk(output_u, type_vtk_u, decomp_vtk_u, name_vtk_u,
			       file_out, double_prec, ascii);
		      
		      WriteVtk(output_du, type_vtk_du, decomp_vtk_du, name_vtk_du,
			       file_out, double_prec, ascii);
		    }
		
		file_out.close();
	      }
	    else
	      {
		for (int i = 0; i < ngrid; i++)
		  if (var_interp(i).GetOutputType() == type)
		    WriteVtk(output_vector, var_interp(i),
			     DOSSIER_output+file_name, double_prec, ascii);
	      }
          }
        else
          {
            cout << "not implemented " << endl;
            abort();
          }
      }
  }


  //! writing a snapshot on the disk
  /*!
    \param[in] nb_iter number of the iteration
    \param[in] t time associated to the current solution
    \param[in] Uh vector solution
    \param[in,out] var_interp definition of the output in different meshes
    \param[in] increment if true the snapshot number is incremented
  */
  template<class Dimension> template<class T>
  void VarOutputProblem_Dim<Dimension>::
  WriteSnapshot(int nb_iter, const Real_wp& t, Vector<T>& Uh,
		const Vector<MeshInterpolation<Dimension> >& var_interp,
		bool increment) const
  {
    bool snapshot = false;
    for (int i = 0; i < this->output_mesh_param.GetM(); i++)
      if (this->output_mesh_param(i).SnapshotToStore(t))
	{
	  string numero = to_str(this->output_mesh_param(i).GetSnapshotNumber());
	  if (increment)
	    this->output_mesh_param(i).IncrementSnapshot();
	  
	  // using the name of the diffracted to store the name of the file
	  // e.g. total.XXX.bb
	  string file_name = GetBaseString(this->output_mesh_param(i).GetTotalFieldFile()) 
	    + "." + numero + string(".bb");
	  
	  const_cast<ParamOutputClass&>(this->output_mesh_param(i)).SetDiffractedFieldFile(file_name);
	  snapshot = true;
	}
    
    Vector<Vector<T> > Uvec;
    this->ReshapeVector(Uh, Uvec);
    
    if (snapshot)
      this->WriteOutputFileGen(Uvec, var_interp, this->output_mesh_param, 0);
    
    Uvec.Nullify();
  }


  //! writing a snapshot on the disk
  /*!
    \param[in] nb_iter number of the iteration
    \param[in] t time associated to the current solution
    \param[in] Uh vector solution
    \param[in,out] var_interp definition of the output in different predefined grids
    \param[in] increment if true, the snapshot number is incremented
  */
  template<class Dimension> template<class T>
  void VarOutputProblem_Dim<Dimension>::
  WriteSnapshot(int nb_iter, const Real_wp& t, Vector<T>& Uh,
		const ParamOutputClass& output_param,
		const Vector<GridInterpolationFull<Dimension> >& var_interp,
		bool increment) const
  {
    int add_total_field = 0;
    if (increment)
      {
	if (var_source_base.GetSourceType(0) == var_source_base.SRC_DIFFRACTED_FIELD)
	  add_total_field = -1;
	else if (var_source_base.GetSourceType(0) == var_source_base.SRC_TOTAL_FIELD)
	  add_total_field = 1;
      }

#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    // for a non-sismogramm output
    // we insert the number of the snapshot in the file name
    string name_file, name_file_diff;
    if ( (!increment) || (!output_param.IsSeismogramOutput()))
      {
        string numero = NumberToString(output_param.GetSnapshotNumber());
        name_file = GetBaseString(output_param.GetFileName(1)) + numero + this->output_file_extension;
        if (add_total_field != 0)
          name_file_diff = GetBaseString(output_param.GetFileName(0)) + numero + this->output_file_extension;
      }
    else
      {
        name_file = GetBaseString(output_param.GetFileName(1)) + ".dat";
        if (add_total_field != 0)
          name_file_diff = GetBaseString(output_param.GetFileName(0)) + ".dat";
      }
    
    Vector<Vector<T> > Uvec;
    this->ReshapeVector(Uh, Uvec);
    
    if ( (!increment) || (!output_param.IsSeismogramOutput()))
      {	
        this->WriteOutputFileGen(Uvec, var_interp, 0, name_file,
				 add_total_field, name_file_diff, output_param, false);
        
	if (rank_proc == 0)
	  if (var_problem.print_level >= 1)
	    cout << "we write output on file at t= " << t << endl;
      }
    else
      {
	int ngrid = 0;
	for (int i = 0; i < var_interp.GetM(); i++)
	  if (var_interp(i).GetOutputType() == output_param.GetNature())
	    ngrid++;
	
	// in the case of a seismogramm, values are written on the same
	// file (values are appended along the time, to the same file)
	// and values are written in ascii format
	// each grid will contain only one point
	Vector<Vector<Vector<T> > > trace_u(ngrid), trace_u_diff(ngrid);	
        Vector<Vector<Vector<T> > > trace_grad(ngrid), trace_grad_diff(ngrid);
        int nb_comp_all = var_problem.GetNbComponentsAll();
        int nb_comp_grad_all = var_problem.GetNbComponentsGradientAll();
	int component = output_param.GetComponent();
	bool compute_grad = output_param.GradientToCompute(nb_comp_all);
	int nb_components_all = 1;
	if (component == -1)
	  nb_components_all = nb_comp_all + nb_comp_grad_all;
	else if (component == -2)
	  nb_components_all = nb_comp_all;
	
	// all components for u are computed
	int ng = 0;
	for (int i = 0; i < var_interp.GetM(); i++)
	  if (var_interp(i).GetOutputType() == output_param.GetNature())
	    {
	      this->ComputeInterpolationU_GradU(Uvec, this->all_points_display,
						var_interp(i), trace_u(ng), trace_grad(ng),
						add_total_field, t, trace_u_diff(ng),
						trace_grad_diff(ng), compute_grad);
	      ng++;
	    }
	
	// computed values are appended at the end of the file
	if (rank_proc == 0)
	  {
	    bool double_prec = (this->precision_output_file == OutputTypeEnum::DOUBLE_PRECISION);
	    const_cast<ParamOutputClass&>(output_param).InitSismo(this->DOSSIER_output, name_file, name_file_diff, 
                                                                  add_total_field, double_prec);
	    
            Vector<T> val(ngrid*nb_components_all + 1);
            Vector<T> val_diff;
            int nb_val_diff = 0;
            if (add_total_field != 0)
              {
                val_diff.Reallocate(ngrid*nb_components_all+1);
                val_diff(0) = t; nb_val_diff++;
              }
            
            val(0) = t; int nb_val = 1;
            
	    // storing the adequate components depending on the choice of the user
	    // we take into account that each grid consists of a single point
	    if ((component == -1) || (component == -2))
	      {
                for (int i = 0; i < ngrid; i++)
                  {
		    for (int j = 0; j < trace_u(i).GetM(); j++)
                      val(nb_val++) = trace_u(i)(j)(0);
                    
		    for (int j = 0; j < trace_u_diff(i).GetM(); j++)
                      val_diff(nb_val_diff++) = trace_u_diff(i)(j)(0);
                  }
		
		if (component == -1)
                  for (int i = 0; i < ngrid; i++)
                    {
		      for (int j = 0; j < trace_grad(i).GetM(); j++)
			val(nb_val++) = trace_grad(i)(j)(0);
                      
		      for (int j = 0; j < trace_grad_diff(i).GetM(); j++)
			val_diff(nb_val_diff++) = trace_grad_diff(i)(j)(0);
                    }
	      }
	    else if (component < nb_comp_all)
	      {
		for (int i = 0; i < ngrid; i++)
		  {
                    if (trace_u(i).GetM() > 0)
		      val(nb_val++) = trace_u(i)(component)(0);
		    
                    if (trace_u_diff(i).GetM() > 0)
		      val_diff(nb_val_diff++) = trace_u_diff(i)(component)(0);
                  }
	      }
	    else
	      {
		// asked component of the gradient is stored
		component -= nb_comp_all;
                for (int i = 0; i < ngrid; i++)
		  {
                    if (trace_grad(i).GetM() > 0)
		      val(nb_val++) = trace_grad(i)(component)(0);

                    if (trace_grad_diff(i).GetM() > 0)
		      val_diff(nb_val_diff++) = trace_grad_diff(i)(component)(0);
                  }
	      }
	    
	    const_cast<ParamOutputClass&>(output_param).AddVectorSismo(val, val_diff);
	  }
      }
    
    Uvec.Nullify();
    if (increment)
      output_param.IncrementSnapshot();
  }


  //! writing of outputs at time t
  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>
  ::WriteSnapshot(int& nb_iter, const Real_wp& t, VectReal_wp& Uh, int num_file) const
  {
    // snapshot
    bool increment = (num_file == 1);
    for (int i = 0; i < output_grid_param.GetM(); i++)
      if (output_grid_param(i).SnapshotToStore(t))
	this->WriteSnapshot(nb_iter, t, Uh, output_grid_param(i), var_grid, increment);
            
    if (output_mesh_param.GetM() != 0)
      this->WriteSnapshot(nb_iter, t, Uh, var_mesh_interp, increment);
  }


  //! writing of outputs at time t
  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>
  ::WriteSnapshot(int& nb_iter, const Real_wp& t, VectComplex_wp& Uh, int num_file) const
  {
    // snapshot
    bool increment = (num_file == 1);
    for (int i = 0; i < output_grid_param.GetM(); i++)
      if (output_grid_param(i).SnapshotToStore(t))
	this->WriteSnapshot(nb_iter, t, Uh, output_grid_param(i), var_grid, increment);
            
    if (output_mesh_param.GetM() != 0)
      this->WriteSnapshot(nb_iter, t, Uh, var_mesh_interp, increment);
  }
    

  /********************
   * VarOutputProblem *
   ********************/


  //! interpolation of U0 over a grid, result is placed in trace_vec
  template<class Dimension> template<class T>
  void VarOutputProblem<Dimension>
  ::ComputeInterpolationUlocGen(const Vector<Vector<T> > & U0,
				const GridInterpolation<Dimension> & var_interp,
				Vector<T>& trace_vec, Vector<T>& trace_grad_vec,
				const IVect& list_points, int nnz, bool compute_grad) const
  {
    int nb_modes = var_boundary.GetNbModes();
    if (var_boundary.ModesNotStored())
      nb_modes = 1;

    int nb_mesh_num = var_problem.GetNbMeshNumberings();
    Vector<int> type_u(nb_mesh_num);
    for (int n = 0; n < nb_mesh_num; n++)
      {
        if (n == 0)
          type_u(n) = var_problem.type_element;
        else
          type_u(n) = var_problem.other_type_element(n-1);
      }

    // counting the number of components (for the solution and its gradient)
    int nb_vec = U0.GetM()/nb_modes;
    Vector<int> vec_mesh_num(nb_vec);
    int nb_comp_all = 0; int nb_grad_all = 0;
    Vector<int> off_phi(nb_vec+1), off_dphi(nb_vec+1);
    off_phi(0) = 0; off_dphi(0) = 0;
    for (int i = 0; i < nb_vec; i++)
      {
        int n = var_problem.mesh_num_unknown(i%var_problem.nb_unknowns);
        vec_mesh_num(i) = n;
        
        int nb_u = var_problem.GetNbComponentsUnknown(n);
        int nb_du = var_problem.GetNbComponentsGradient(n);
        nb_comp_all += nb_u; off_phi(i+1) = off_phi(i) + nb_u;
        nb_grad_all += nb_du; off_dphi(i+1) = off_dphi(i) + nb_du;
      }

    Vector<IVect> Nodle(nb_mesh_num);
    Vector<T> val_u(nb_comp_all);
    Vector<T> val_grad_u(nb_grad_all);
    Vector<T> val_v(3*nb_modes);
    Vector<T> val_v_grad(3*nb_modes);

    Vector<Vector<TinyVector<Real_wp, 1> > > phi_H1(nb_mesh_num);
    Vector<Vector<TinyVector<Real_wp, Dimension::dim_N> > > grad_phi_H1(nb_mesh_num);

    Vector<Vector<TinyVector<Real_wp, Dimension::dim_N> > > phi_Hdiv(nb_mesh_num);
    Vector<Vector<TinyVector<Real_wp, 1> > > div_phi_Hdiv(nb_mesh_num);

    enum {nb_u_hcurl = FiniteElementHcurl<Dimension>::nb_components_u,
          nb_du_hcurl = FiniteElementHcurl<Dimension>::nb_components_grad};
    
    Vector<Vector<TinyVector<Real_wp, nb_u_hcurl> > > phi_Hcurl(nb_mesh_num);
    Vector<Vector<TinyVector<Real_wp, nb_du_hcurl> > > curl_phi_Hcurl(nb_mesh_num);
    
    Vector<R_N> s;
    R_N point_loc;
    
    FftInterface<Complex_wp> fft_interface;
    var_source.InitFftComputation(fft_interface);
    
    int nb = 0;
    trace_vec.Reallocate(nnz*nb_comp_all);
    trace_vec.Fill(0);
    if (compute_grad)
      {
	trace_grad_vec.Reallocate(nnz*nb_grad_all);
	trace_grad_vec.Fill(0);
      }
    
    for (int i1 = 0; i1 < list_points.GetM(); i1++)
      {
	int i = list_points(i1);
	FillZero(val_u);
	FillZero(val_grad_u);
	if (i >= 0)
	  {
	    // iquad is the element where the point i is
	    int iquad = var_interp.GetElementNumber(i);
	    // local coordinates in the element
	    point_loc = var_interp.GetLocalCoordinate(i);
	    int n = var_interp.GetSectionNumber(i);
	    
	    // if the point has been found in the grid
	    if ((iquad >= 0)&&(iquad < var_problem.mesh.GetNbElt()))
	      {
                // we get vertices of the element
		var_problem.mesh.GetVerticesElement(iquad, s);

                // basis functions are computed and dof numbers are retrieved
		for (int nm = 0; nm < nb_mesh_num; nm++)
                  {
                    Nodle(nm) = var_problem.GetDofNumberOnElement(iquad, nm);
                    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);

                    switch(type_u(nm))
                      {
                      case 1:
                        {
                          const ElementReference<Dimension, 1>& Fb = var_problem.GetReferenceElementH1(iquad, nm); 
                        
                        // we compute the values of basis functions at point_loc
                        Fb.ComputeValuesPhi(point_loc, phi_H1(nm), var_interp.GetDFjm1(i),
                                            mesh_num, iquad);
                        
                        if (compute_grad)
                          Fb.ComputeValuesGradientPhi(point_loc, grad_phi_H1(nm), var_interp.GetDFjm1(i),
                                                      mesh_num, iquad);
                        }
                        break;
                      case 2:
                        {
                          const ElementReference<Dimension, 2>& Fb = var_problem.GetReferenceElementHcurl(iquad, nm); 
                        
                        // we compute the values of basis functions at point_loc
                        Fb.ComputeValuesPhi(point_loc, phi_Hcurl(nm), var_interp.GetDFjm1(i),
                                            mesh_num, iquad);
                        
                        if (compute_grad)
                          Fb.ComputeValuesGradientPhi(point_loc, curl_phi_Hcurl(nm), var_interp.GetDFjm1(i),
                                                      mesh_num, iquad);
                        }
                        break;
                      case 3 :
                        {
                          const ElementReference<Dimension, 3>& Fb = var_problem.GetReferenceElementHdiv(iquad, nm); 
                        
                        // we compute the values of basis functions at point_loc
                        Fb.ComputeValuesPhi(point_loc, phi_Hdiv(nm), var_interp.GetDFjm1(i),
                                            mesh_num, iquad);
                        
                        if (compute_grad)
                          Fb.ComputeValuesGradientPhi(point_loc, div_phi_Hdiv(nm), var_interp.GetDFjm1(i),
                                                      mesh_num, iquad);
                        }
                        break;
                      }
                  }

                // loop over vectors contained in Uà
		for (int q = 0; q < nb_vec; q++)
		  {
		    FillZero(val_v); FillZero(val_v_grad);
                    int nm = vec_mesh_num(q);
		    for (int j = 0; j < Nodle(nm).GetM(); j++)
		      {
			int num_ddl = Nodle(nm)(j);
			switch (type_u(nm))
                          {
                          case 1:
                            {
                              if (num_ddl >= 0)
                                for (int p = 0; p < nb_modes; p++)
                                  {
                                    val_v(p) += U0(p + q*nb_modes)(num_ddl)*phi_H1(nm)(j)(0);
                                    if (compute_grad)
                                      this->AddContribInterpol(U0(p + q*nb_modes)(num_ddl),
                                                               grad_phi_H1(nm)(j), p, val_v_grad);
                                  }
                            }
                            break;
                          case 2:
                            {
                              if (num_ddl >= 0)
                                for (int p = 0; p < nb_modes; p++)
                                  {
                                    this->AddContribInterpol(U0(p + q*nb_modes)(num_ddl), phi_Hcurl(nm)(j), p, val_v);
                                    if (compute_grad)
                                      this->AddContribInterpol(U0(p + q*nb_modes)(num_ddl),
                                                               curl_phi_Hcurl(nm)(j), p, val_v_grad);
                                  }                              
                            }
                            break;
                          case 3:
                            {
                              if (num_ddl >= 0)
                                for (int p = 0; p < nb_modes; p++)
                                  {
                                    this->AddContribInterpol(U0(p + q*nb_modes)(num_ddl), phi_Hdiv(nm)(j), p, val_v);
                                    if (compute_grad)
                                      val_v_grad(p) += U0(p + q*nb_modes)(num_ddl)*div_phi_Hdiv(nm)(j)(0);
                                  }                              
                            }
                            break;
                          }
		      }
		    
		    this->GetModalOutput(n, nm, fft_interface, val_v, val_v_grad,
					 off_phi(q), val_u, off_dphi(q), val_grad_u, compute_grad);		    
                  }
                
		this->ModifyOutputUnknown(val_u, val_grad_u, i, var_interp, iquad, compute_grad);
	      }
            
	    for (int p = 0; p < nb_comp_all; p++)
              trace_vec(nb + p*nnz) = val_u(p);

            if (compute_grad)
              for (int p = 0; p < nb_grad_all; p++)
                trace_grad_vec(nb + p*nnz) = val_grad_u(p);
	    
	    nb++;
	  }
      }
  }


  //! Computation of the solution u on nodal points of the mesh
  /*!
    \param[in] U0 components of the solution on degrees of freedom
    \param[out] Unodal values of the solution on nodal points of the mesh
    Unodal(i,j) is the value of u on nodal point j of the element i
    \param[out] GradNodal values of the gradient on nodal points of the mesh
    \param[in] u_component  if true all components of u are computed
    \param[in] grad_component  if true all components of gradient of u are computed
  */
  template<class Dimension> template<class T>
  void VarOutputProblem<Dimension>
  ::ComputeNodalUgradU_Gen(const Vector<Vector<T> >& U0,
			   Vector<Vector<Vector<T> > >& Unodal,
			   Vector<Vector<Vector<T> > >& GradNodal,
			   bool u_component, bool grad_component,
			   bool hess_component, bool true_unknown, int nb_u) const
  {
    if (nb_u == -1)
      nb_u = var_problem.nb_unknowns;
    
    int nb_mesh_num = var_problem.GetNbMeshNumberings();
    Vector<int> type_u(nb_mesh_num);
    for (int n = 0; n < nb_mesh_num; n++)
      {
        if (n == 0)
          type_u(n) = var_problem.type_element;
        else
          type_u(n) = var_problem.other_type_element(n-1);
      }

    // counting the number of components (for the solution and its gradient)
    Vector<int> vec_mesh_num(nb_u);
    int nb_comp_all = 0; int nb_grad_all = 0;
    Vector<int> off_phi(nb_u+1), off_dphi(nb_u+1), off_hess(nb_u+1);
    off_phi(0) = 0; off_dphi(0) = 0; off_hess.Zero();
    Vector<int> IndexMesh(nb_u); Vector<IVect> IndexUnknown(nb_mesh_num);
    for (int i = 0; i < nb_u; i++)
      {
        int n = var_problem.mesh_num_unknown(i%var_problem.nb_unknowns);
        vec_mesh_num(i) = n;
        
        int nb_comp = var_problem.GetNbComponentsUnknown(n);
        int nb_grad = var_problem.GetNbComponentsGradient(n);
        nb_comp_all += nb_comp; off_phi(i+1) = off_phi(i) + nb_comp;
        nb_grad_all += nb_grad; off_dphi(i+1) = off_dphi(i) + nb_grad;

        IndexMesh(i) = IndexUnknown(n).GetM();
        IndexUnknown(n).PushBack(i);
        
        if (hess_component && (nb_comp > 1))
          {
            cout << "Not implemented" << endl;
            abort();
          }
      }

    if (hess_component)
      {
        off_hess(0) = nb_grad_all;
        int nb_comp_hess = Dimension::dim_N*(Dimension::dim_N+1)/2;
        for (int i = 0; i < nb_u; i++)
          {
            nb_grad_all += nb_comp_hess; 
            off_hess(i+1) = nb_grad_all;            
          }
      }
    
    if (nb_grad_all <= 0)
      grad_component = false;

    Vector<Vector<Vector<T> > > Uloc;
    Vector<Vector<T> > Uloc_node, dUloc_node, d2Uloc_node, curlUloc_node, zero_vec;
    Uloc.Reallocate(nb_mesh_num);
    for (int i = 0; i < nb_mesh_num; i++)
      Uloc(i).Reallocate(IndexUnknown(i).GetM());
    
    // allocation of output arrays
    if (u_component)
      {
	Unodal.Reallocate(var_problem.mesh.GetNbElt());
	for (int i = 0; i < Unodal.GetM(); i++)
	  Unodal(i).Reallocate(nb_comp_all);
      }
    
    if (grad_component)
      {
	GradNodal.Reallocate(var_problem.mesh.GetNbElt());
	for (int i = 0; i < GradNodal.GetM(); i++)
	  GradNodal(i).Reallocate(nb_grad_all);
      }
    
    // loop over elements
    VectR_N s;
    SetPoints<Dimension> PointsElem;
    SetMatrices<Dimension> MatricesElem;
    Vector<IVect> Nodle(nb_mesh_num);
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
	const ElementReference_Dim<Dimension>& Fb_g = var_problem.GetReferenceElement(i);
	
	// we get the vertices of the face
	var_problem.mesh.GetVerticesElement(i, s);
	
	// transformation on nodal points
	Fb_g.FjElemNodal(s, PointsElem, var_problem.mesh, i);
	Fb_g.DFjElemNodal(s, PointsElem, MatricesElem, var_problem.mesh, i);
	
	// we get values of U on the element i
	for (int nm = 0; nm < nb_mesh_num; nm++)
          Nodle(nm) = var_problem.GetDofNumberOnElement(i, nm);
        
	for (int m = 0; m < nb_u; m++)
	  {
            int nm = vec_mesh_num(m); int im = IndexMesh(m);
            const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
	    Uloc(nm)(im).Reallocate(Nodle(nm).GetM());
	    for (int j = 0; j < Nodle(nm).GetM(); j++)
	      Uloc(nm)(im)(j) = U0(m)(Nodle(nm)(j));
	    
	    mesh_num.number_map.
	      ModifyLocalComponentVector(mesh_num, Uloc(nm)(im), i);
	  }

	// projection on nodal points
        for (int n = 0; n < nb_mesh_num; n++)
          {
            int nb_comp_u = 1, nb_comp_g = 1;
            switch (type_u(n))
              {
              case 1:
                {
                  nb_comp_g = Dimension::dim_N;
                  const ElementReference<Dimension, 1>& Fb = var_problem.GetReferenceElementH1(i, n);
                  
                  Uloc_node.Reallocate(Uloc(n).GetM());
                  Fb.ComputeNodalValues(MatricesElem, Uloc(n), Uloc_node, var_problem.mesh, i);
                  if (grad_component)
                    {
                      dUloc_node.Reallocate(Dimension::dim_N*Uloc_node.GetM());
                      Fb.ComputeNodalGradient(MatricesElem, Uloc_node, dUloc_node);
                    }
                  
                  if (hess_component)
                    {
                      // if the hessian is required
                      // we store it in GradNodal
                      int d = Dimension::dim_N;
                      int d2 = Dimension::dim_N*Dimension::dim_N;
                      d2Uloc_node.Reallocate(Uloc_node.GetM()*d2);
                      Fb.ComputeNodalGradient(MatricesElem, dUloc_node, d2Uloc_node);
                      
                      for (int j = 0; j < Uloc_node.GetM(); j++)
                        {
                          int num = off_hess(IndexUnknown(n)(j));
                          for (int k = 0; k < Dimension::dim_N; k++)
                            for (int n = k; n < Dimension::dim_N; n++)
                              GradNodal(i)(num++) = d2Uloc_node(j*d2 + k*d + n);
                        }
                    }
                }
                break;
              case 2:
                {
                  nb_comp_u = Dimension::dim_N; nb_comp_g = FiniteElementHcurl<Dimension>::nb_components_grad;
                  const ElementReference<Dimension, 2>& Fb = var_problem.GetReferenceElementHcurl(i, n);
                  
                  Fb.ComputeNodalValues(MatricesElem, Uloc(n), Uloc_node, var_problem.mesh, i);                  
                  if (grad_component)
                    {
                      dUloc_node.Reallocate(Dimension::dim_N*Uloc_node.GetM());
                      Fb.ComputeNodalGradient(MatricesElem, Uloc_node, dUloc_node);
                      Fb.GetCurlFromGradient(dUloc_node, curlUloc_node);
                      dUloc_node = curlUloc_node;
                    }
                }
                break;
              case 3:
                {
                  nb_comp_u = Dimension::dim_N;
                  const ElementReference<Dimension, 3>& Fb = var_problem.GetReferenceElementHdiv(i, n);
                  
                  Fb.ComputeNodalValues(MatricesElem, Uloc(n), Uloc_node, var_problem.mesh, i);                  
                  if (grad_component)
                    {
                      dUloc_node.Reallocate(Dimension::dim_N*Uloc_node.GetM());
                      Fb.ComputeNodalGradient(MatricesElem, Uloc_node, dUloc_node);
                      Fb.GetCurlFromGradient(dUloc_node, curlUloc_node);
                      dUloc_node = curlUloc_node;
                    }
                }
                break;
              }

            if (u_component)
              for (int j = 0; j < Uloc_node.GetM(); j++)
                Unodal(i)(off_phi(IndexUnknown(n)(j/nb_comp_u)) + j%nb_comp_u) = Uloc_node(j);

            if (grad_component)
              for (int j = 0; j < dUloc_node.GetM(); j++)
                GradNodal(i)(off_dphi(IndexUnknown(n)(j/nb_comp_g)) + j%nb_comp_g) = dUloc_node(j);        
          }
        
	if (!true_unknown)
	  {
	    if (grad_component)
	      this->ModifyOutputUnknown(Unodal(i), GradNodal(i), i, u_component, grad_component);
	    else
	      this->ModifyOutputUnknown(Unodal(i), zero_vec, i, u_component, grad_component);
	  }
      }
  }
  

  //! Computation of the solution u on quadrature points of the mesh
  /*!
    \param[in] U0 components of the solution on degrees of freedom
    \param[out] Unodal values of the solution on quadrature points of the mesh
    Unodal(i,j) is the value of u on quadrature point j of the element i
    \param[out] GradNodal values of the gradient on quadrature points of the mesh
    \param[in] u_component if true u is computed on quadrature points
    \param[in] grad_component if true the gradient of u is computed on quadrature points
    \param[in] true_unknown if true the real value u is provided,
    otherwise ModifyOutputUnknown is called
  */
  template<class Dimension> template<class T>
  void VarOutputProblem<Dimension>
  ::ComputeQuadratureUgradU_Gen(const Vector<Vector<T> >& U0,
				Vector<Vector<Vector<T> > >& Unodal,
				Vector<Vector<Vector<T> > >& GradNodal,
				bool u_component, bool grad_component, bool true_unknown, int nb_u) const
  {
    if (nb_u == -1)
      nb_u = var_problem.nb_unknowns;
    
    int nb_comp_all = var_problem.GetNbComponentsAll(nb_u);
    int nb_comp_grad = var_problem.GetNbComponentsGradientAll(nb_u);
    if (nb_comp_grad <= 0)
      grad_component = false;

    // allocation of output arrays
    if (u_component)
      {
	Unodal.Reallocate(var_problem.mesh.GetNbElt());
	for (int i = 0; i < Unodal.GetM(); i++)
	  Unodal(i).Reallocate(nb_comp_all);
      }
    
    if (grad_component)
      {
	GradNodal.Reallocate(var_problem.mesh.GetNbElt());
	for (int i = 0; i < GradNodal.GetM(); i++)
	  GradNodal(i).Reallocate(nb_comp_grad);
      }
    
    // loop over elements
    VectR_N s;
    SetPoints<Dimension> PointsElem;
    SetMatrices<Dimension> MatricesElem;
    int nb_mesh_num = var_problem.GetNbMeshNumberings();
    Vector<IVect> Nodle(nb_mesh_num);
    Vector<Vector<Vector<T> > > Uloc, Uloc_node, dUloc_node;
    Uloc.Reallocate(nb_mesh_num); Uloc_node.Reallocate(nb_mesh_num);
    dUloc_node.Reallocate(nb_mesh_num);

    IVect nb_u_num(nb_mesh_num); nb_u_num.Zero();
    for (int m = 0; m < nb_u; m++)
      {
        int nm = var_problem.mesh_num_unknown(m%var_problem.nb_unknowns);
        nb_u_num(nm)++;
      }

    IVect type_u(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      {
        Uloc(nm).Reallocate(nb_u_num(nm));
        if (nm == 0)
          type_u(nm) = var_problem.type_element;
        else
          type_u(nm) = var_problem.other_type_element(nm-1);
      }

    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
	// we get the vertices of the face
	var_problem.mesh.GetVerticesElement(i, s);
	const ElementReference_Dim<Dimension>& Fb = var_problem.GetReferenceElement(i);
	int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    
	// transformation on quadrature points
	var_problem.FillQuadratureJacobian(i, nb_points_quad, Fb.GetGeometricElement(),
					   s, PointsElem, MatricesElem);
	
	// we get values of U on the element i
	for (int nm = 0; nm < nb_mesh_num; nm++)
          Nodle(nm) = var_problem.GetDofNumberOnElement(i, nm);
        
        nb_u_num.Zero();
	for (int m = 0; m < nb_u; m++)
	  {
            int nm = var_problem.mesh_num_unknown(m%var_problem.nb_unknowns);
            int p = nb_u_num(nm);
            const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
	    Uloc(nm)(p).Reallocate(Nodle(nm).GetM());
	    for (int j = 0; j < Nodle(nm).GetM(); j++)
	      Uloc(nm)(p)(j) = U0(m)(Nodle(nm)(j));
	    
            if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
              mesh_num.number_map.
                ModifyLocalComponentVector(mesh_num, Uloc(nm)(p), i);
            
            nb_u_num(nm)++;
	  }
        
	// projection on quadrature points
	for (int nm = 0; nm < nb_mesh_num; nm++)
          switch(type_u(nm))
            {
            case 1:
              var_problem.GetReferenceElementH1(i, nm).
                ComputeQuadratureValues(MatricesElem, Uloc(nm), Uloc_node(nm), dUloc_node(nm),
                                        u_component, grad_component, var_problem.mesh, i);
              break;
            case 2:
              var_problem.GetReferenceElementHcurl(i, nm).
                ComputeQuadratureValues(MatricesElem, Uloc(nm), Uloc_node(nm), dUloc_node(nm),
                                        u_component, grad_component, var_problem.mesh, i);
              break;
            case 3:
              var_problem.GetReferenceElementHdiv(i, nm).
                ComputeQuadratureValues(MatricesElem, Uloc(nm), Uloc_node(nm), dUloc_node(nm),
                                        u_component, grad_component, var_problem.mesh, i);
              break;
            }
        
        // assuming that unknowns are sorted
        if (u_component)
	  {
            int offset = 0;
            for (int nm = 0; nm < nb_mesh_num; nm++)
              {
                for (int j = 0; j < Uloc_node(nm).GetM(); j++)
                  Unodal(i)(offset++) = Uloc_node(nm)(j);
              }
          }
	
	if (grad_component)
	  {
            int offset = 0;
            for (int nm = 0; nm < nb_mesh_num; nm++)
              {
                for (int j = 0; j < dUloc_node(nm).GetM(); j++)
                  GradNodal(i)(offset++) = dUloc_node(nm)(j);
              }
          }
        
	// results are modified if necessary
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    if (!true_unknown)
	      {
		abort();
		// appeler ModifyOutputUnknown
	      }
	  }
      }
  }

  template<class Dimension> template<class T>
  void VarOutputProblem_Dim<Dimension>
  ::WriteUquadratureGen(const Vector<Vector<T> >& U0,
                        const string& file_name_u, const string& file_name_poids) const
  {
    Vector<Vector<Vector<T> > > Uquad, grad_Uquad;    
    this->ComputeQuadratureUgradU(U0, Uquad, grad_Uquad, true, false, true);

    // we compute the number of points
    int nb_pts = 0;
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      nb_pts += Uquad(i)(0).GetM();
    
    int nb_u = Uquad(0).GetM();

    // computes the integral of all components
    Vector<Vector<T> > all_val(nb_u); VectReal_wp poids(nb_pts);
    for (int i = 0; i < nb_u; i++)
      {
        all_val(i).Reallocate(nb_pts);
        all_val(i).Zero();
      }

    Vector<T> sum_u(nb_u); sum_u.Zero(); int count = 0;
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
        const VectReal_wp& omega = var_problem.WeightsND(i);
        bool affine = (var_problem.Glob_jacobian(i).GetM() == 1);

        Real_wp jacob;
        for (int j = 0; j < Uquad(i)(0).GetM(); j++)
          {
            if (affine)
              jacob = var_problem.Glob_jacobian(i)(0)*omega(j);
            else
              jacob = var_problem.Glob_jacobian(i)(j);
            
            for (int k = 0; k < nb_u; k++)
              {
                all_val(k)(count) = Uquad(i)(k)(j);
                sum_u(k) += jacob*Uquad(i)(k)(j);
              }

            poids(count) = jacob; 
            count++;            
          }
      }
    
    if (var_problem.print_level >= 2)
      cout << "Integral of u = " << sum_u << endl;
    
    string root = GetBaseString(file_name_u);
    for (int k = 0; k < nb_u; k++)
      all_val(k).Write(root + "_U" + to_str(k) + ".dat");
    
    poids.Write(file_name_poids);
  }

  template<class Dimension> 
  void VarOutputProblem_Dim<Dimension>
  ::WriteUquadrature(const Vector<VectReal_wp >& U0,
                     const string& file_name_poids, const string& file_name_u) const
  {
    WriteUquadratureGen(U0, file_name_poids, file_name_u);
  }


  template<class Dimension>
  void VarOutputProblem_Dim<Dimension>
  ::WriteUquadrature(const Vector<VectComplex_wp>& U0,
                     const string& file_name_poids, const string& file_name_u) const
  {
    WriteUquadratureGen(U0, file_name_poids, file_name_u);
  }

  template<class Dimension> template<class T>
  void VarOutputProblem_Dim<Dimension>
  :: ComputeValueNodalBoundary(const Vector<T>& Uloc, Vector<T>& Uboundary, int num_elem, int num_loc) const
  {
    switch (var_problem.type_element)
      {
      case 1:
        var_problem.GetReferenceElementH1(num_elem).ComputeValueNodalBoundary(Uloc, Uboundary, num_loc); break;
      case 2:
        var_problem.GetReferenceElementHcurl(num_elem).ComputeValueNodalBoundary(Uloc, Uboundary, num_loc); break;
      case 3:
        var_problem.GetReferenceElementHdiv(num_elem).ComputeValueNodalBoundary(Uloc, Uboundary, num_loc); break;
      }
  }
  
}

#define MONTJOIE_FILE_OUTPUT_HARMONIC_CXX
#endif
