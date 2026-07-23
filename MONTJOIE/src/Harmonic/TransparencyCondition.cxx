#ifndef MONTJOIE_FILE_TRANSPARENCY_CONDITION_CXX

namespace Montjoie
{

  /**************************
   * VarComputationRCS_Base *
   **************************/
  
  
  //! Sets attributes to default values
  template<class Dimension>
  void VarComputationRCS_Base<Dimension>::InitDefaultValues()
  {
    print_level = 0;
    rcs_to_be_computed = false;
    first_angle_RCS = 0;
    last_angle_RCS = pi_wp;
    type_rcs = BISTATIC_RCS;
    LINE_RCS = 1;
    nb_angles_RCS = 0;
    order_integration_rcs = 1;
    file_RCS = string("rcs.dat");    
    ref_domain_RCS = 1;
  }
  

  template<class Dimension>
  size_t VarComputationRCS_Base<Dimension>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille += var_mesh.GetMemorySize() + PointsOutside.GetMemorySize();
    return taille;
  }

  
  //! modification of parameters with a line of the data file
  template<class Dimension>
  void VarComputationRCS_Base<Dimension>
  ::SetInputData(const string& description_field,
		 const VectString& parameters)
  {
    if (!description_field.compare("ParametersRCS"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of VarComputationRCS_Base" << endl;
	    cout << "ParametersRCS needs more parameters, for instance :" << endl;
	    cout << "ParametersRCS = YES ref AUTO" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("YES"))
	  rcs_to_be_computed = true;
	else
	  rcs_to_be_computed = false;
	
        LINE_RCS = to_num<int>(parameters(1));
	if (!parameters(2).compare("AUTO"))
	  order_integration_rcs = var_problem.GetMeshNumbering(0).GetOrder();
	else
	  order_integration_rcs = to_num<int>(parameters(2));
	
	if (parameters.GetM() > 3)
	  {
	    if (!parameters(3).compare("MONOSTATIC"))
	      type_rcs = MONOSTATIC_RCS;
	    else
	      type_rcs = BISTATIC_RCS;	    
	  }
	else
	  type_rcs = BISTATIC_RCS;
        
        if (parameters.GetM() > 4)
          ref_domain_RCS = to_num<int>(parameters(4));
        else
          ref_domain_RCS = 1;
      }
    else if (!description_field.compare("AngleRCS"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of VarComputationRCS_Base" << endl;
	    cout << "AngleRCS needs more parameters, for instance :" << endl;
	    cout << "AngleRCS = teta0 teta_final N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	first_angle_RCS = atof(parameters(0).data())*pi_wp/180.0-pi_wp;
	last_angle_RCS = atof(parameters(1).data())*pi_wp/180.0-pi_wp;
	nb_angles_RCS = atoi(parameters(2).data());
      }
    else if (!description_field.compare("FileRCS"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarComputationRCS_Base" << endl;
	    cout << "FileRCS needs more parameters, for instance :" << endl;
	    cout << "FileRCS = file_name" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	file_RCS  = string(parameters(0));
	file_far_field.clear();
	if (parameters.GetM() > 1)
	  file_far_field = parameters(1);	
      }
  }

  
  //! initialisation before computing far field
  template<class Dimension>
  void VarComputationRCS_Base<Dimension>::InitComputationRCS(bool assemble)
  {
    if (!rcs_to_be_computed)
      return;
    
    // calcul des points et poids de quadrature
    int ref = this->LINE_RCS;
    
    this->var_mesh.SetGaussQuadrature(this->order_integration_rcs);
    IVect ref_cond(var_problem.mesh.GetNbReferences()+1); ref_cond.Fill(0);
    for (int i = 1; i <= var_problem.mesh.GetNbReferences(); i++)
      if (var_problem.mesh.GetBodyNumber(i) == ref)
        ref_cond(i) = 1;
    
    Mesh<Dimension> mesh_subdiv;
    this->var_mesh.InitProjectionSurface(var_problem.mesh);
    this->var_mesh.ComputeSurfaceMesh(ref_cond, var_problem.mesh, mesh_subdiv,
                                      var_problem, this->ref_domain_RCS);
    
#ifdef SELDON_WITH_MPI
    this->var_mesh.GatherQuadraturePoints(var_problem.comm_group_mode, assemble);
    int Nall = this->var_mesh.GetNbAllQuadraturePoints();
    int N(Nall);
    MPI_Allreduce(&N, &Nall, 1, MPI_INTEGER, MPI_MAX, var_problem.comm_group_mode);
#else
    this->var_mesh.GatherQuadraturePoints();
    int Nall = this->var_mesh.GetNbAllQuadraturePoints();
#endif
    
    if (Nall <= 0)
      {
        cout << "No quadrature points to compute radar cross section " << endl;
        abort();
      }
  }  
  

  /***************************
   * TransparencySolver_Base *
   ***************************/

  
  //! Sets values of attributes to default values
  void TransparencySolver_Base::InitDefaultValues()
  {
    order_integration_transparency_condition = 1;
    use_transparency_condition = false;
    ref_domain_transparent = 1;
    
    iterative_solver.SelectIterativeSolver(iterative_solver.GMRES);
    iterative_solver.SetMaxNumberIteration(100);
    iterative_solver.SetRestart(10);
    iterative_solver.SetStoppingCriterion(1e-6);
    iterative_solver.SetPrintLevel(4);

    display_messages = true;
    nb_prod = 0;
    first_point_to_treat = 0;
    nb_points_to_treat = 0;
  }


  //! modification of parameters with a line of the data file    
  void TransparencySolver_Base
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("TransparencyCondition"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of TransparencySolver_Base" << endl;
	    cout << "TransparencyCondition needs more parameters, for instance :" << endl;
	    cout << "TransparencyCondition = YES ref AUTO" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("YES"))
	  use_transparency_condition = true;
	else
	  use_transparency_condition = false;
	
	ref_transparent_body = to_num<int>(parameters(1));
	if (!parameters(2).compare("AUTO"))
	  this->order_integration_transparency_condition = var_problem.GetDefaultOrder();
	else
	  this->order_integration_transparency_condition = to_num<int>(parameters(2));
        
        if (parameters.GetM() > 3)
          ref_domain_transparent = to_num<int>(parameters(3));
        else
          ref_domain_transparent = 1;
      }
    else if (!description_field.compare("ParamResolutionTransparency"))
      {
	if (parameters.GetM() <= 3)
	  {
	    cout << "In SetInputData of TransparencySolver_Base" << endl;
	    cout << "ParamResolutionTransparency needs more parameters, for instance :" << endl;
	    cout << "ParamResolutionTransparency = type_resol nmax_iter tolerance restart" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	Vector<string> param1(2), param2(1), param_res(1);
	param_res(0) = parameters(0);
	param1(0) = parameters(1);
	param1(1) = parameters(3);
	param2(0) = parameters(2);
	string key_res("TypeResolution");
	string key1("NumberMaxIterations");
	string key2("Tolerance");
	iterative_solver.SetInputData(key_res, param_res);
	iterative_solver.SetInputData(key1, param1);
	iterative_solver.SetInputData(key2, param2);
      }
  }
  
  
  //! modification of the solution in order to take in account transparent condition
  /*!
    \param[in] b_source solution with first-order absorbing boundary condition on input,
    \param[in,out] x_solution solution with transparent condition (initial guess in input)
  */
  void TransparencySolver_Base
  ::Solve(VectComplex_wp & x_solution, VectComplex_wp& b_source)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    // use of an iterative method to solve
    // (I - A_c^{-1} A_p) U = G
    // where A_c is the sparse finite element matrix
    // and A_p is the dense matrix arising from integral representation
    // identity preconditioner
    if (nb_proc == 1)      
      iterative_solver.Solve(*this, x_solution, b_source);
    else
      {
#ifdef SELDON_WITH_MPI
        // using distributed vectors for parallel version
        DistributedVector<Complex_wp>* source
	  = var_comm.AllocateDistributedVector(b_source);
        
        DistributedVector<Complex_wp>* solution
	  = var_comm.AllocateDistributedVector(x_solution);
                
        iterative_solver.Solve(*this, *solution, *source);
	
	var_comm.NullifyDistributedVector(source);
	var_comm.NullifyDistributedVector(solution);
#endif
      }
  }
  
  
  //! initialization before iterations needed to obtain a transparent condition
  void TransparencySolver_Base::Init()
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int nb_proc(1);
#endif

    // size of the system to solve
    this->m_ = var_problem.GetNbDof();
    this->n_ = var_problem.GetNbDof();
    
    // reference of absorbing boundary condition (outside boundary Sigma)
    int ref_abc = BoundaryConditionEnum::LINE_ABSORBING;
    
    // reference of intern boundary (Gamma)
    int ref_gamma = ref_transparent_body;
    
    // computes surface mesh for Gamma and the absorbing boundary condition
    // Nall is the total number of quadrature points on the absorbing boundary condition
    Vector<int> offset_abc_proc;
    int Nall = this->ComputeSurfaceGammaAndAbsorbing(ref_abc, ref_gamma, offset_abc_proc);
    
    int nb_points_effective = Nall;
    int num_point0 = 0;
    first_point_to_treat = num_point0;
    nb_points_to_treat = nb_points_effective;
    
    if (nb_proc > 1)
      {
#ifdef SELDON_WITH_MPI
        // for points on absorbing boundary condition, we redistribute them
        // on all the processors        
        GetParallelDistributionPoints(nb_proc, rank_proc, Nall, nb_points_effective, num_point0);
        first_point_to_treat = num_point0;
        nb_points_to_treat = nb_points_effective;
        
        // which points to send ?
        points_to_send.Clear();
        int n0 = num_point0, n1 = n0 + nb_points_effective;
        TinyVector<int, 3> num;
        for (int i = 0; i < nb_proc; i++)
          {
            int p0 = offset_abc_proc(i);
            int p1 = offset_abc_proc(i+1);
            int na = 0, nb = -1;
            
            // intersection between [n0, n1] and [p0, p1] ?
            if (n0 <= p0)
              {
                if (n1 > p0)
                  {
                    na = p0;
                    nb = min(p1, n1);
                  }            
              }
            else
              {
                if (n0 < p1)
                  {
                    na = n0;
		    nb = min(p1, n1);
                  }
              }
            
            if (nb > na)
              {
                // we have to send points in [na, nb] to processor i
                num.Init(i, na-n0, nb-n0);
                points_to_send.PushBack(num);
              }
          }
        
        // which points to receive ?
        IVect nb_points_treat(nb_proc);
        nb_points_treat.Fill(0);
	MPI_Allgather(&nb_points_effective, 1, MPI_INTEGER,
                      nb_points_treat.GetData(), 1, MPI_INTEGER, var_comm.comm_group_mode);
	
        IVect offset_points(nb_proc+1);
        offset_points(0) = 0;
        for (int i = 0; i < nb_proc; i++)
          offset_points(i+1) = offset_points(i) + nb_points_treat(i);
        
        points_to_recv.Clear();        
        n0 = offset_abc_proc(rank_proc);
        n1 = offset_abc_proc(rank_proc+1);
        for (int i = 0; i < nb_proc; i++)
          {
            int p0 = offset_points(i), p1 = offset_points(i+1);
            int na = 0, nb = -1;
            
            // intersection between [n0, n1] and [p0, p1] ?
            if (n0 <= p0)
              {
                if (n1 > p0)
                  {
                    na = p0;
                    nb = min(p1, n1);
                  }            
              }
            else
              {
                if (n0 < p1)
                  {
                    na = n0;
                    nb = min(p1, n1);
                  }
              }
            
            if (nb > na)
              {
                // we have to send points in [na, nb] to processor i
                num.Init(i, na-n0, nb-n0);
                points_to_recv.PushBack(num);
              }
          }
        
#endif
        
      }
     
    iterative_solver.SetPrintLevel(4);
    if (var_problem.print_level < 6)
      {
        if (glob_solver.IsIterativeSolver())
          {
            if (var_problem.print_level < 3)
              glob_solver.SetPrintLevel(0);
            else
              glob_solver.SetPrintLevel(1);
          }
        else
          glob_solver.SetPrintLevel(0);
      }
  }
  
  
  //! Exchange values (EnSigma, HnSigma) and (EnQuad, HnQuad)
  void TransparencySolver_Base
  ::ExchangeValues(VectComplex_wp& EnSigma, VectComplex_wp& HnSigma,
		   VectComplex_wp& EnQuad, VectComplex_wp& HnQuad) const
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif

    if (nb_proc > 1)
      {
#ifdef SELDON_WITH_MPI
	int nb_en = var_problem.nb_components_en;
	int nb_hn = var_problem.nb_components_hn;
        int nsend = points_to_send.GetM();
        int nrecv = points_to_recv.GetM();
        Vector<MPI_Request> request_Esend(nsend), request_Hsend(nsend);
        Vector<MPI_Request> request_Erecv(nrecv), request_Hrecv(nrecv);
	Vector<int64_t> xtmp;
        for (int i = 0; i < nsend; i++)
          {
            int proc = points_to_send(i)(0);
            int n0 = points_to_send(i)(1);
            int n1 = points_to_send(i)(2);
            int size_send = nb_en*(n1-n0);
	    request_Esend(i) = 
	      MpiIsend(var_comm.comm_group_mode, &EnSigma(n0*nb_en), xtmp, size_send, proc, 77);
            
            size_send = nb_hn*(n1-n0);
	    request_Hsend(i) = 
	      MpiIsend(var_comm.comm_group_mode, &HnSigma(n0*nb_hn), xtmp, size_send, proc, 78);
          }
        
        for (int i = 0; i < nrecv; i++)
          {
            int proc = points_to_recv(i)(0);
            int n0 = points_to_recv(i)(1);
            int n1 = points_to_recv(i)(2);
            int size_recv = nb_en*(n1-n0);
	    request_Erecv(i) =
	      MpiIrecv(var_comm.comm_group_mode, &EnQuad(n0*nb_en), xtmp,
		       size_recv, proc, 77);

            size_recv = nb_hn*(n1-n0);
	    request_Hrecv(i) =
	      MpiIrecv(var_comm.comm_group_mode, &HnQuad(n0*nb_hn), xtmp,
		       size_recv, proc, 78);
          }
        
	MPI_Status status;
        for (int i = 0; i < nsend; i++)
          {
	    MPI_Wait(&request_Esend(i), &status);
	    MPI_Wait(&request_Hsend(i), &status);
          }

        for (int i = 0; i < nrecv; i++)
          {
	    MPI_Wait(&request_Erecv(i), &status);
	    MPI_Wait(&request_Hrecv(i), &status);
          }
        
#endif
      }
    else
      {
	Copy(EnSigma, EnQuad);
	Copy(HnSigma, HnQuad);
      }
  }
    
    
  //! Computes 3-D green's function (for Helmholtz), gradient and hessian
  /*!
    X, Y : two points for which phi(x,y) is evaluated
    phi  : green function phi(x,y)
    grad_phi : gradient of green function according y
    hessian_phi : Hessian matrix of green function (d^2  phi(x,y))/ dy_m dy_l
  */
  void TransparencySolver_Base::
  ComputeGreenKernel(const R3& X, const R3& Y, const Real_wp& omega, Complex_wp& phi,
		     R3_Complex_wp& grad_phi, Matrix3_3sym_Complex_wp& hessian_phi)
  {
    Real_wp omega2 = omega*omega;

    // in 3-D \phi(x,y) = exp(ik|x-y|) / (4 pi |x-y| )
    R3 xMinusy = X - Y;

    // T = |x-y|   T2 = |x-y|^2   arg = ik |x-y|
    Real_wp T = Norm2(xMinusy); Complex_wp arg = Iwp*omega*T;
    Real_wp inv_T = 1.0/T,inv_T2 = inv_T*inv_T;
    
    
    // evaluation of green function
    phi = exp(arg) * inv_T / (4*pi_wp); Complex_wp phi_div_T = phi * inv_T;
    
    
    // evaluation of its gradient
    grad_phi = xMinusy;    
    Complex_wp alpha = (-Iwp*omega + inv_T)*phi_div_T; Mlt( alpha , grad_phi);    
    
    // evaluation of hessian matrix
    // initialization to (x_m - y_m) (x_l - y_l)
    GetNormalProjector(xMinusy, hessian_phi);    

    alpha = (-omega2 - Real_wp(3)*Iwp*omega*inv_T + Real_wp(3) * inv_T2) * inv_T;
    Mlt(alpha, hessian_phi);
    alpha = Iwp*omega - inv_T;
    hessian_phi(0, 0) += alpha; hessian_phi(1, 1) += alpha; hessian_phi(2, 2) += alpha;
    
    Mlt(phi_div_T, hessian_phi);
  }
  
  
  //! Computes green's function (Helmholtz), gradient and hessian
  /*!
    X, Y : two points for which phi(x,y) is evaluated
    phi  : green function phi(x,y)
    grad_phi : gradient of green function according y
  */
  void TransparencySolver_Base
  ::ComputeGreenKernel(const R3& X, const R3& Y, const Real_wp& omega, Complex_wp& phi, R3_Complex_wp& grad_phi)
  {
    // in 3-D \phi(x,y) = exp(ik|x-y|) / (4 pi |x-y| )
    R3 xMinusy = X - Y;
    
    // T = |x-y|   T2 = |x-y|^2   arg = ik |x-y|
    Real_wp T = Norm2(xMinusy); Complex_wp arg = Iwp*omega*T;
    Real_wp inv_T = 1.0/T;
    
    // evaluation of green function
    phi = exp(arg) * inv_T / (4*pi_wp); Complex_wp phi_div_T = phi * inv_T;
    
    // evaluation of its gradient
    grad_phi = xMinusy;
    
    Complex_wp alpha = (-Iwp*omega + inv_T)*phi_div_T;
    Mlt(alpha, grad_phi);
  }


  void TransparencySolver_Base
  ::MltAddVector(const Real_wp& alpha, const VectReal_wp& B, const Real_wp& beta,
		 VectReal_wp& C) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }

  
  //! matrix-vector product C = beta C + alpha A B where A is the "transparent" iteration matrix
  void TransparencySolver_Base
  ::MltAddVector(const Complex_wp& alpha, const VectComplex_wp& B, const Complex_wp& beta,
		 VectComplex_wp& C) const
  {
    Complex_wp zero(0, 0), one(1, 0);
    if (beta == zero)
      C.Fill(zero);
    else
      Mlt(beta, C);
    
    int N = this->GetM();
    VectComplex_wp g_source(N), C_tmp(N);
    
    // g_source = Ap B
    this->ComputeRightHandSide(B, g_source);
    
    // we solve linear system Ac C = Ap B, so C = Ac^{-1} Ap B 
    C_tmp.Zero();
    this->ComputeSolution(g_source, C_tmp);
    
    // and finally C = alpha*(B - Ac^{-1} Ap B) 
    Add(-one, B, C_tmp);
    Add(-alpha, C_tmp, C);
    
    nb_prod++;
  }


  void TransparencySolver_Base
  ::MltAddVector(const Real_wp& alpha, const SeldonTranspose&,
                 const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }

  
  void TransparencySolver_Base
  ::MltAddVector(const Complex_wp& alpha, const SeldonTranspose&,
                 const VectComplex_wp& B, const Complex_wp& beta, VectComplex_wp& C) const
  {
    cout << "Matrix-vector product with transpose not implemented" << endl;
    abort();
  }
  

  void TransparencySolver_Base::MltVector(const VectReal_wp& B, VectReal_wp& C) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }

  
  //! matrix-vector product C = A B where A is the "transparent" iteration matrix
  void TransparencySolver_Base::MltVector(const VectComplex_wp& B, VectComplex_wp& C) const
  {
    Complex_wp zero(0, 0), one(1, 0);
    MltAddVector(one, B, zero, C);
  }


  void TransparencySolver_Base::MltVector(const SeldonTranspose&, const VectReal_wp& B, VectReal_wp& C) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }

  
  void TransparencySolver_Base::MltVector(const SeldonTranspose&, const VectComplex_wp& B, VectComplex_wp& C) const
  {
    cout << "Matrix-vector product with transpose not implemented" << endl;
    abort();

  }


  /**************************
   * TransparencySolver_Dim *
   **************************/
  

  template<class Dimension>
  size_t TransparencySolver_Dim<Dimension>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille += this->iterative_solver.GetMemorySize();
    taille += this->points_to_send.GetMemorySize() + this->points_to_recv.GetMemorySize();
    taille += interp_abc.GetMemorySize() + interp_gamma.GetMemorySize();
    return taille;
  }


  //! extracts meshes for Gamma and Absorbing boundary condition
  template<class Dimension>
  int TransparencySolver_Dim<Dimension>
  ::ComputeSurfaceGammaAndAbsorbing(int ref_abc, int ref_gamma, IVect& offset_abc_proc)
  {
    IVect ref_cond(var_problem.mesh.GetNbReferences()+1);
    ref_cond.Fill(0);
    for (int i = 1; i < ref_cond.GetM(); i++)
      if (var_problem.mesh.GetBodyNumber(i) == ref_gamma)
        ref_cond(i) = 1;
      
    Mesh<Dimension> mesh_subdiv;
    interp_gamma.SetGaussQuadrature(this->order_integration_transparency_condition);    
    interp_gamma.InitProjectionSurface(var_problem.mesh);
    interp_gamma.ComputeSurfaceMesh(ref_cond, var_problem.mesh, mesh_subdiv,
				    var_problem, ref_domain_transparent);
    
#ifdef SELDON_WITH_MPI
    interp_gamma.GatherQuadraturePoints(var_problem.comm_group_mode);
#else
    interp_gamma.GatherQuadraturePoints();
#endif
    
    ref_cond.Fill(0);
    for (int i = 1; i < ref_cond.GetM(); i++)
      if (var_problem.mesh.GetBoundaryCondition(i) == ref_abc)
        ref_cond(i) = 1;
    
    interp_abc.SetGaussQuadrature(var_problem.GetMeshNumbering(0).GetOrder());    
    interp_abc.InitProjectionSurface(var_problem.mesh);
    interp_abc.ComputeSurfaceMesh(ref_cond, var_problem.mesh, mesh_subdiv,
				  var_problem);
    
#ifdef SELDON_WITH_MPI
    interp_abc.GatherQuadraturePoints(var_problem.comm_group_mode);
#else
    interp_abc.GatherQuadraturePoints();
#endif
  
    int Nall = interp_abc.GetNbAllQuadraturePoints();
    offset_abc_proc = interp_abc.offset_quadrature_per_proc;
    return Nall;
  }


  //! right hand side for transparent condition
  template<class Dimension>
  void TransparencySolver_Dim<Dimension>
  ::ComputeRightHandSide(const VectComplex_wp& U0, VectComplex_wp& g_source) const
  {
    g_source.Zero();
    
    // loop on the edges/faces of sigma (absorbing boundary condition
    int nb_boundaries_ref = interp_abc.GetNbBoundary();
    Vector<Vector<Complex_wp> > contrib(var_problem.nb_unknowns);
    Vector<Vector<Complex_wp> > scal_g(var_problem.nb_unknowns);
    
    SetPoints<Dimension> PointsElem;
    SetMatrices<Dimension> MatricesElem;
    
    // computation of E \times n and H \times n
    // on all quadrature points of Gamma
    VectComplex_wp trace_En;
    VectComplex_wp trace_Hn;
    var_problem.ComputeEnHnOnBoundary(interp_gamma, U0, trace_En, trace_Hn);
    
    /* ofstream file_En("En.dat"), file_Hn("Hn.dat");
       file_En.precision(15); file_Hn.precision(15); 
       for (int i = 0; i < trace_En.GetM(); i++)
       {
       file_En << real(trace_En(i)) << " " << imag(trace_En(i)) << endl;
       file_Hn << real(trace_Hn(i)) << " " << imag(trace_Hn(i)) << endl;
       }
    
       file_En.close();
       file_Hn.close(); */
    
    //  computation of E^pot and H^pot on all quadrature points of Sigma
    // with representation formula
    R_N pointX, normaleX; 
    int k0 = first_point_to_treat;
    int k1 = k0 + nb_points_to_treat;
    VectComplex_wp EnSigma(var_problem.nb_components_en*nb_points_to_treat);
    VectComplex_wp HnSigma(var_problem.nb_components_hn*nb_points_to_treat);
    for (int k = k0; k < k1; k++)
      {
        pointX = interp_abc.GetQuadraturePoint(k);
        normaleX = interp_abc.GetQuadratureNormale(k);
	
	this->ComputeAndStoreEnPot(trace_En, trace_Hn, pointX, normaleX,
				   EnSigma, HnSigma, k-k0);	
      }
    
    // exchanging values between processors
    // final values are contained in EnQuad, HnQuad
    int Nquad = interp_abc.GetNbPointsQuadrature();
    VectComplex_wp EnQuad(Nquad*var_problem.nb_components_en);
    VectComplex_wp HnQuad(Nquad*var_problem.nb_components_hn);
    
    this->ExchangeValues(EnSigma, HnSigma, EnQuad, HnQuad);
    
    Real_wp omega = var_problem.GetOmega();
    Real_wp rho0(1), mu0(1);
    var_source.GetCoefAB_Infinity(mu0, rho0);
    Real_wp k_inf = omega*sqrt(rho0/mu0);

    int Nvol = var_problem.GetMeshNumbering(0).GetNbDof();
    int nb = 0;
    for (int i = 0; i < nb_boundaries_ref; i++)
      {        
	int num_elem = interp_abc.GetElementNumberOfSurface(i);

        // we retrieve number of dofs
	IVect Nodle = var_problem.GetDofNumberOnElement(num_elem);
        
        int num_loc = interp_abc.GetLocalPositionOfSurface(i);
	int nb_dof_loc = var_problem.GetMeshNumbering(0).GetNbLocalDof(num_elem);
	int nb_points_quadrature = interp_abc.GetNbPointsQuadrature(i);
        for (int m = 0; m < var_problem.nb_unknowns; m++)
	  {
            int nm = var_problem.mesh_num_unknown(m);
            scal_g(m).Reallocate(nb_points_quadrature*var_problem.GetNbComponentsUnknown(nm));
          }
        
	for (int j = 0; j < nb_points_quadrature; j++)
	  {
            pointX = interp_abc.PointsQuadrature(i, j);
	    normaleX = interp_abc.NormaleQuadrature(i, j);
            
	    this->GetSource(EnQuad, HnQuad, nb, k_inf, pointX, normaleX, scal_g, j);
	    nb++;	    
	  }

        const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
        int type_elt = var_problem.type_element;
	switch (type_elt)
          {
          case 1:
            var_problem.GetReferenceElementH1(num_elem).
              ComputeGaussIntegralSurface(interp_abc.Dfjm1Quadrature(i),
                                          interp_abc.WeightsQuadrature(i),
                                          interp_abc.DsQuadrature(i),
                                          scal_g, contrib, mesh_num, num_elem, num_loc);
            break;
          case 2:
            var_problem.GetReferenceElementHcurl(num_elem).
              ComputeGaussIntegralSurface(interp_abc.Dfjm1Quadrature(i),
                                          interp_abc.WeightsQuadrature(i),
                                          interp_abc.DsQuadrature(i),
                                          scal_g, contrib, mesh_num, num_elem, num_loc);
            break;
          case 3:
            var_problem.GetReferenceElementH1(num_elem).
              ComputeGaussIntegralSurface(interp_abc.Dfjm1Quadrature(i),
                                          interp_abc.WeightsQuadrature(i),
                                          interp_abc.DsQuadrature(i),
                                          scal_g, contrib, mesh_num, num_elem, num_loc);
            break;
          }
        
        for (int j = 0; j < nb_dof_loc; j++)
	  for (int m = 0; m < var_problem.nb_unknowns; m++)
	    {
              int num_dof_loc = Nodle(j);
              int num_dof = num_dof_loc + m*Nvol;
              if (num_dof_loc >= 0)
                g_source(num_dof) += contrib(m)(j);
            }
      }

    //int test_input; cout << "we wait " << endl; cin >> test_input;    
  } // end method ComputeRightHandSide


  //! computes E^pot and H^pot and store it in EnStore and HnStore
  template<class Dimension, int nb_en, int nb_hn>
  void TransparencySolver_Fem<Dimension, nb_en, nb_hn>
  ::ComputeAndStoreEnPot(const VectComplex_wp& trace_En, const VectComplex_wp& trace_Hn,
			 const R_N& pointX, const R_N& normaleX,
			 VectComplex_wp& EnStore, VectComplex_wp& HnStore, int k) const
  {
    TinyVector<Complex_wp, nb_en> val_En;
    TinyVector<Complex_wp, nb_hn> val_Hn;

    this->ComputeIntegralRepresentation(trace_En, trace_Hn, this->interp_gamma, 
					pointX, normaleX, val_En, val_Hn);
    
    CopyVector(val_En, k, EnStore);
    CopyVector(val_Hn, k, HnStore);
  }
  
}

#define MONTJOIE_FILE_TRANSPARENCY_CONDITION_CXX
#endif
