#ifndef MONTJOIE_FILE_TIME_REVERSAL_CXX

namespace Montjoie
{
  
  /**********************
   * TimeReversalSource *
   **********************/
  
  
  template<class Dimension>
  void TimeReversalSource<Dimension>::SetCurrentTime(const Real_wp& t)
  {
    t_current = t;
  }
  
  
  template<class Dimension>
  const Real_wp& TimeReversalSource<Dimension>::GetCurrentTime() const
  {
    return t_current;
  }
  
  
  //! inits interpolation (in time)
  template<class Dimension>
  void TimeReversalSource<Dimension>::InitInterpolation()
  {
    t_current = -1.0;
    n0_buffer = -1; n1_buffer = -1;
    t0_buffer = -1.0; t1_buffer = -1.0;
    
    // time intervals are regular
    VectReal_wp equi_points(order+1);
    for (int i = 0; i <= order; i++)
      equi_points(i) = Real_wp(i)/order;
    
    lob.AffectPoints(equi_points);
    lob.ComputeGradPhi(1e3*epsilon_machine);
  }


  //! computes an interpolation of v \cdot n at the time trev
  template<class Dimension>
  void TimeReversalSource<Dimension>::ComputeSource(const Real_wp& trev, const Real_wp& epsilon_time,
						    int nb_deriv)
  {
    if ((trev < t0_buffer-epsilon_time) || (trev > t1_buffer+epsilon_time))
      this->InitSourceSurface(trev, epsilon_time);
    
    t_current = trev;
    
    // computation of interpolation of v.n at the correct time
    int N = ValueVn(0).GetM();
    Real_wp lambda = (ValueTime(0) - trev) / (ValueTime(0) - ValueTime(order));
    
    Real_wp coef_scale = 1.0/(ValueTime(order) - ValueTime(0));
    if ((lambda <= -1e-5)||(lambda >=1.0001))
      {
	cout << "t outside [t0,t1]" << endl;
	DISP(trev); DISP(lambda);
	DISP(ValueTime);
        abort();
      }
    
    // computation of Lagrange function phi_i(lambda)
    VectReal_wp coef_phi(order+1), coef_dphi(order+1);
    for (int k = 0; k <= order; k++)
      coef_phi(k) = lob.EvaluatePhi(k, lambda);
    
    for (int i = 1; i < nb_deriv; i++)
      {
	// if derivative of source is needed
	// we compute the derivatives d^n phi_i(lambda)/dt^n
	coef_dphi.Fill(0);
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    coef_dphi(j) += lob.GradPhi(j, k)*coef_phi(k);
	
	Mlt(coef_scale, coef_dphi);
	coef_phi = coef_dphi;
      }
    
    // v \, = \, \sum v_i \phi_i (or derivatives)
    InterpolVn.Reallocate(N); InterpolVn.Fill(0);
    for (int j = 0; j <= order; j++)
      Add(coef_phi(j), ValueVn(j), InterpolVn);	    
  }


  //! initialization before computation of source for time reversal experiment
  template<class Dimension>
  void TimeReversalSource<Dimension>::InitSourceSurface(const Real_wp& t, const Real_wp& epsilon_time)
  {
    int N = time_reversal.mesh_interp.GetNbPointsQuadrature();
    typedef typename TimeAcousticReversal<Dimension>::value_type value_type;
    value_type t0, t1;
    
    // we extract the stored value of v dot n
    ifstream file_in(time_reversal.file_out_dtn.data());

    if (!file_in)
      {
	cout << "the file doesn't exist " << time_reversal.file_out_dtn << endl;
	abort();
      }    
    
    // size of an entry
    int size_entry = (N+1)*sizeof(value_type) + sizeof(int);
    Real_wp deltat = time_reversal.deltat;
    bool new_interval = false;
    if (n0_buffer == -1)
      {
	// first initialization
	n0_buffer = toInteger(round(t/deltat)) + 2;
	new_interval = true;
	
	// size of the file ?
	file_in.seekg(0, ios::end);
	int size_file = file_in.tellg();
	n0_buffer = min(n0_buffer, size_file/size_entry - 1);
	
	bool t_not_on_interval = true;
	// we decrement n0_buffer until we find t in interval [t1,t0]
	while (t_not_on_interval)
	  {
	    int new_position = size_entry*n0_buffer;
	    // DISP(N); DISP(new_position); DISP(sizeof(real));
	    file_in.seekg(new_position, ios::beg);
	    file_in.read(reinterpret_cast<char*>(&t0), sizeof(value_type)); // DISP(t0);
	    if (n0_buffer >= order)
	      {
		new_position = size_entry*(n0_buffer-order);
		file_in.seekg(new_position, ios::beg);
		file_in.read(reinterpret_cast<char*>(&t1), sizeof(value_type));
	      }
	    else
	      t1 = toDouble(time_reversal.tbegin);
	    
	    if (n0_buffer == 0)
	      t_not_on_interval = false;
	    else if ((t > Real_wp(t1)-epsilon_time)&&(t < Real_wp(t0)+epsilon_time))
	      t_not_on_interval = false;
	    else
	      n0_buffer--;
	  }	
	  
	if ((t <= Real_wp(t1)-epsilon_time)&&(t >= epsilon_time))
	  {
	    cout << "intervalle non trouve " << endl;
	    exit(0);
	  }
      }
    else
      {
	if ((t > t1_buffer-epsilon_time)&&(t < t0_buffer+epsilon_time))
	  {
	  }
	else
	  {
	    new_interval = true;
	    // decalage de n0_buffer
	    n0_buffer = n1_buffer;
	  }
      }
    
    if (new_interval)
      {
	n0_buffer = max(n0_buffer, order);
	n1_buffer = max(n0_buffer - order, 0);
	
	// we read all the snapshots between n0_buffer and n1_buffer
	ValueVn.Reallocate(order+1); ValueTime.Reallocate(order+1);
	Vector<value_type> v_dot_n;
	
	for (int i = 0; i <= order; i++)
	  {
	    int new_position = size_entry*(n0_buffer-i);
	    file_in.seekg(new_position, ios::beg);
	    file_in.read(reinterpret_cast<char*>(&t1), sizeof(value_type));
	    ValueTime(i) = t1;
	    v_dot_n.Read(file_in); ValueVn(i).Reallocate(N);
	    // DISP(t1); DISP(new_position); DISP(v_dot_n);
	    for (int j = 0; j < N; j++)
	      ValueVn(i)(j) = v_dot_n(j);
	    
	    // we can apply some filters on ValueVn
	    time_reversal.ApplyFilterVn(ValueVn(i), ValueTime(i));
	  }
        	
	t0_buffer = ValueTime(0); t1_buffer = ValueTime(order);
      }
    
    file_in.close();
  }

  
  //! returns true if there is a source here
  template<class Dimension>
  bool TimeReversalSource<Dimension>::IsNonNullSurfacicSource(int ref)
  {
    if (time_reversal.is_ref_on_mode(ref) == 1)
      return true;
    
    return false;
  }
  
  
  //! evaluation of surfacic source
  template<class Dimension>
  void TimeReversalSource<Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
                           const SetMatrices<Dimension>& MatricesElem, VectReal_wp& f)
  {
    f.Fill(0);
    int num_seg = time_reversal.IndexSurface(this->num_bound_ref_);
    if (num_seg != -1)
      f(0) = InterpolVn(time_reversal.OffsetQuad(num_seg) +  k);
    
  }
  
  
  /************************
   * TimeAcousticReversal *
   ************************/
  
  
  //! modification of the data file (ex : guide.ini)
  /*!
    \param[in] description_field : first keyword of a line of the data file
    \param[in] parameters : list of values on the line
    \param[in] nb_param : number of values
   */
  template<class Dimension>
  void TimeAcousticReversal<Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("TimeReversal"))
      {
	if (parameters.GetM() <= 6)
	  {
	    cout << "In SetInputData of TimeAcousticReversal" << endl;
	    cout << "TimeReversal needs more parameters, for instance :" << endl;
	    cout << "TimeReversal = DIRECT t0 tf dt file_mesh file_dtn ref" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	presence_time_reversal = true;
	if (!parameters(0).compare("DIRECT"))
	  {
	    // direct simulation (datas are stored in files)
	    type_simulation = DIRECT;
	  }
	else if (!parameters(0).compare("INVERSE"))
	  {
	    // inverse simulation
	    type_simulation = INVERSE;
	  }
        
	// common parameters to direct and inverse simulation
	tbegin = to_num<Real_wp>(parameters(1));
	tend = to_num<Real_wp>(parameters(2));
	deltat = to_num<Real_wp>(parameters(3));
	
	file_mesh = parameters(4);
	file_out_dtn = parameters(5);
	is_ref_on_mode.Reallocate(100);
        is_ref_on_mode.Fill(0);
        is_ref_on_mode(to_num<int>(parameters(6))) = 1;	
      }
    else if (!description_field.compare("TimeReversalParameters"))
      {
 	// if you want to read additional informations on the data file
	// you can put that in this section
	// related to a line of the data file : 
	// TimeReversalParameters = parameters(0) parameters(1) ...
	
	// if you want to define an another keyword, you add an else if group like
	// else if (!description_field.compare("MyKeyword"))
	// {
	// }	
      }
  }
  

  //! initialization of some variables
  template<class Dimension>
  void TimeAcousticReversal<Dimension>::InitComputation()
  {
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    
    // we compute quadrature points and normale for this boundary mesh
    this->mesh_interp.SetGaussQuadrature(mesh_num.GetOrder(),
					 Globatto<Real_wp>::QUADRATURE_LOBATTO);
    
    Mesh<Dimension> mesh_subdiv;
    this->mesh_interp.ComputeSurfaceMesh(this->is_ref_on_mode, var_problem.mesh,
					 mesh_subdiv, var_problem);
    
#ifdef SELDON_WITH_MPI
    this->mesh_interp.GatherQuadraturePoints(var_problem.comm_group_mode, false);
#else
    this->mesh_interp.GatherQuadraturePoints(false);
#endif
    
    if (type_simulation == DIRECT)
      {
	src_reverse.SetCurrentTime(tbegin);
	// old file is deleted if existing
	remove(this->file_out_dtn.data());
      }
    else
      {
	int N = this->mesh_interp.GetNbPointsQuadrature();
	OffsetQuad.Reallocate(N+1);
	OffsetQuad(0) = 0;
	for (int i = 0; i < this->mesh_interp.GetNbBoundary(); i++)
	  OffsetQuad(i+1) = OffsetQuad(i) + this->mesh_interp.GetNbPointsQuadrature(i);
	    
	IndexSurface.Reallocate(var_problem.mesh.GetNbBoundaryRef());
	IndexSurface.Fill(-1);
	for (int i = 0; i < this->mesh_interp.GetNbBoundary(); i++)
	  {
	    int num_elem = this->mesh_interp.GetElementNumberOfSurface(i);
	    int num_loc = this->mesh_interp.GetLocalPositionOfSurface(i);
	    int num_face = var_problem.mesh.Element(num_elem).numBoundary(num_loc);
	    IndexSurface(num_face) = i;
	  }
	
	src_reverse.InitInterpolation();
      }
  }
  
  
  //! write on files, values of du/dn, when a direct simulation is selected
  template<class Dimension>
  void TimeAcousticReversal<Dimension>
  ::WriteOutput(int nb_iter, const Real_wp& t, const VectReal_wp Vh)
  {
    if (type_simulation == DIRECT)
      if ((t >= src_reverse.GetCurrentTime() - var_time.epsilon_time)
          &&(t <= (tend + var_time.epsilon_time + var_time.GetTimeStep())))
	{
	  Vector<value_type> v_dot_n;          
	  // computation of v \cdot n on Gamma
	  this->GetVn(Vh, v_dot_n);
          
	  // appending values at the end of the output file
	  ofstream file_out(file_out_dtn.data(), ios::app);
	  
	  value_type tmp = t;
	  // we write the time
	  file_out.write(reinterpret_cast<char*>(&tmp), sizeof(value_type));
	  // and the vector
	  v_dot_n.Write(file_out);
	  
	  file_out.close();
	}
  }
  
  
  //! computation of v dot n
  template<class Dimension>
  void TimeAcousticReversal<Dimension>
  ::GetVn(const VectReal_wp& Vn, Vector<value_type>& v_dot_n) const
  {
    int N = this->mesh_interp.GetNbPointsQuadrature();
    v_dot_n.Reallocate(N); v_dot_n.Fill(0);
    // loop on faces of Gamma
    for (int i = 0; i < this->mesh_interp.GetNbBoundary(); i++)
      {
        int num_elem = this->mesh_interp.GetElementNumberOfSurface(i);
        int num_loc = this->mesh_interp.GetLocalPositionOfSurface(i);
        int nb_dof_face = this->mesh_interp.GetNbPointsQuadrature(i);
	const ElementReference<Dimension, 1>& Fb = var_problem.GetReferenceElementH1(num_elem);
                
	// we get values of v
	for (int j = 0; j < nb_dof_face; j++)
	  {
            int num_dof_loc = Fb.GetQuadNumber(num_loc, j);
	    
	    int offset = var_problem.GetOffsetDofV(num_elem) + Dimension::dim_N*num_dof_loc;
	    typename Dimension::R_N normale = this->mesh_interp.NormaleQuadrature(i)(j);
	    
	    // computation of scalar product v. n
            Real_wp sum; SetComplexZero(sum);
	    for (int k = 0; k < Dimension::dim_N; k++)
	      sum += normale(k)*Vn(offset+k);
	    
	    // we store the value in output vector
	    v_dot_n(i*nb_dof_face + j) = value_type(sum);
	  }
      }
  }
  

  //! application of filtering on v dot n
  /*!
    \param[in,out] dU_dn values of \f$ v \cdot n \f$ on the surface
    \param[in] t time
   */
  template<class Dimension>
  void TimeAcousticReversal<Dimension>
  ::ApplyFilterVn(VectReal_wp& dU_dn, const Real_wp& t)
  {
    // here, you can apply filter techniques on v \dot n
    // the time t is the time for the direct simulation 
    // (for a reverse simulation, you have to compute t1 -  t0 - t)
    
    // for a temporal filtering
    // Real_wp coef_filter = MY_FORMULE(t);
    // Mlt(coef_filter, dU_dn);    
  }
  
  
  //! add the source :  b_src <- b_src + alpha*df/dt^nb_deriv
  /*!
    \param[in] t time where we want to evaluate the source
    \param[in] alpha coefficient 
    \param[in] nb_deriv 0 for the primitive of the source f,
                        1 for the source, 2 for the derivative, etc
    \param[in,out] b_src vector to be modified
   */
  template<class Dimension>
  void TimeAcousticReversal<Dimension>
  ::AddScalarSourceAtTime(const Real_wp& alpha, const Real_wp& t, int nb_deriv, VectReal_wp& b_src)
  {
    // only called in the reverse case
    // REVERSE : source is given by values of the output file    
    if ((t >= -var_time.epsilon_time)
	&&(t <= tend-tbegin+var_time.epsilon_time))
      {
	Real_wp trev = tend - tbegin - t;
	// we read files, vectors at different times near the time t
	src_reverse.ComputeSource(trev, var_time.epsilon_time, nb_deriv);
	
	// computation of right-hand side
	Vector<VectReal_wp> bvec;
	bvec.SetData(1, &b_src);
	Vector<VirtualSourceFEM<Real_wp, Dimension>* > f_vec(1);
	f_vec(0) = &src_reverse;
	var_source.AddSurfaceSource(1.0, bvec, f_vec);
	bvec.Nullify();
      }
  }
  
  //! source depends on time ?
  template<class Dimension>
  bool TimeAcousticReversal<Dimension>::SourceDoesNotDependOnTime() const
  {
    if (presence_time_reversal)
      if (type_simulation == INVERSE)
	  return false;
    
    return true;
  }
    
}

#define MONTJOIE_FILE_TIME_REVERSAL_CXX
#endif
