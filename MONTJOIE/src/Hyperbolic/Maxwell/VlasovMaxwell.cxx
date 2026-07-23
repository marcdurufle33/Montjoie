#ifndef MONTJOIE_FILE_VLASOV_MAXWELL_CXX

#include "VlasovMaxwell.hxx"

namespace Montjoie
{
  
  /******************************************************
   * Beams of particles / space charge limited emission *
   ******************************************************/
  
    //! default constructor
  template<class Dimension>
  Vlasov_BeamCurrent<Dimension>::Vlasov_BeamCurrent() 
  {
    I0 = 0; charge = 1; ratio_q0_over_m0 = 1; ref = -1; type_beam = UNIFORM_BEAM; alpha = 0.0;
    gamv_normal_min = 0; gamv_normal_max = 0; height_min = 0; height_max = 0; delta_iter = 1;
    gamv_tangential_min = 0; gamv_tangential_max = 0; uniform_spatial = true;
  }
  
  //! reading of line of data file
  template<class Dimension>
  void Vlasov_BeamCurrent<Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!parameters(0).compare("UNIFORM"))
      type_beam = UNIFORM_BEAM;
    else if (!parameters(0).compare("GAUSSIAN"))
      type_beam = GAUSSIAN_BEAM;
    
    int nb = 1, dim_N = Dimension::dim_N;
    // two extremities of the beam [A,B]
    for (int k = 0; k < dim_N; k++)
      ptA(k) = to_num<Real_wp>(parameters(nb++));
    
    for (int k = 0; k < dim_N; k++)
      ptB(k) = to_num<Real_wp>(parameters(nb++));
    
    if (dim_N >= 3)
      {
	// in 3-D, two other extremites [C,D]
	for (int k = 0; k < dim_N; k++)
	  ptC(k) = to_num<Real_wp>(parameters(nb++));
	
	for (int k = 0; k < dim_N; k++)
	  ptD(k) = to_num<Real_wp>(parameters(nb++));
      }
    
    // relativistic factor multiplied by the speed
    for (int k = 0; k < dim_N; k++)
      v0(k) = to_num<Real_wp>(parameters(nb++));
    
    // emitting current of the beam
    I0 = to_num<Real_wp>(parameters(nb++));
    // definition of the particles emitted (charge and mass)
    SetInputSpecies(charge, ratio_q0_over_m0, parameters, nb);
    // number of particles to create by unit of length
    nb_part_per_unit = to_num<Real_wp>(parameters(nb++));
    
    // number of iterations between two creations of particles
    // if AUTO -> this number is computed depending the velocity and distance
    // between two consecutive particles , in order to have
    // |X_i^n - X_{i+1}^n| = |X_i^{n+1} - X_i^n| 
    if (!parameters(nb++).compare("AUTO"))
      delta_iter = 0;
    else
      delta_iter = to_num<int>(parameters(nb-1));
    
    // for gaussian beam computation of alpha and center
    Add(ptA, ptB, center); Mlt(0.5, center);
    Real_wp radius = ptA.Distance(center);
    alpha = -10.0/(radius*radius);
  }
  
  
  //! reading of a line of data file, related to breakdown current
  template<class Dimension>
  void Vlasov_BeamCurrent<Dimension>
  ::SetInputBreakdown(const string& description_field, const VectString& parameters)
  {
    int nb = 0;
    if (!parameters(nb++).compare("UNIFORM"))
      uniform_spatial = true;
    else
      uniform_spatial = false;
    
    // reference of the boundary
    ref = to_num<int>(parameters(nb++));
    // breakdown field
    breakdown_field = to_num<Real_wp>(parameters(nb++));
    // minimal and maximal normal gamv
    gamv_normal_min = to_num<Real_wp>(parameters(nb++));
    gamv_normal_max = to_num<Real_wp>(parameters(nb++));
    // minimal and maximal tangential gamv
    gamv_tangential_min = to_num<Real_wp>(parameters(nb++));
    gamv_tangential_max = to_num<Real_wp>(parameters(nb++));
    // minimal and maximal height
    height_min = to_num<Real_wp>(parameters(nb++));
    height_max = to_num<Real_wp>(parameters(nb++));
    // definition of the particles emitted (charge and mass)
    SetInputSpecies(charge, ratio_q0_over_m0, parameters, nb);
    // number of particles to create by unit of length
    nb_part_per_unit = to_num<Real_wp>(parameters(nb++));
    
    // number of time iterations before creation of particles
    if (!parameters(nb++).compare("AUTO"))
      delta_iter = 0;
    else
      delta_iter = to_num<int>(parameters(nb-1));
  }
  
  
  //! selecting species
  template<class Dimension>
  void Vlasov_BeamCurrent<Dimension>
  ::SetInputSpecies(Real_wp& charge, Real_wp& coef_q0_over_m0, const VectString& param, int& nb)
  {
    Real_wp mass;
    if (!param(nb).compare("ELECTRON"))
       {
	nb++;
	mass = PhysicalConstant::mass_electron;
	charge = -PhysicalConstant::elementary_charge;
      }
    else if (!param(nb).compare("PROTON"))
      {
	nb++;
	mass = PhysicalConstant::mass_proton;
	charge = PhysicalConstant::elementary_charge;
      }
    else
      {
	charge = to_num<Real_wp>(param(nb++));
	mass = to_num<Real_wp>(param(nb++));
      }
    
    coef_q0_over_m0 = charge/mass;
  }

  
  template<class Dimension>
  Real_wp Vlasov_BeamCurrent<Dimension>::GetFunction(const Real_wp& lambda)
  {
    if (type_beam == GAUSSIAN_BEAM)
      {
	Real_wp x = 2.0*lambda-1;
	return exp(-10.0*x*x);
      }
    
    if ((lambda < 0)||(lambda > 1))
      return 0;
    
    return 1.0;
  }
  
  template<class Dimension>
  Real_wp Vlasov_BeamCurrent<Dimension>::GetFunction(const Real_wp& lambda, const Real_wp& mu)
  {
    if (type_beam == GAUSSIAN_BEAM)
      {
	Real_wp x = (2.0*lambda-1)*(2.0*lambda-1) + (2.0*mu-1.0)*(2.0*mu-1.0);
	return exp(-10.0*x);
      }
    
    if ((lambda < 0)||(lambda > 1)||(mu < 0)||(mu > 1))
      return 0;
    
    return 1.0;
  }
  
  
  /************************
   * CurrentSource_Vlasov *
   ************************/
  
  
  //! initialization for the definition of particles shape functions 
  template<class Dimension> 
  CurrentSource_Vlasov<Dimension>
  ::CurrentSource_Vlasov(HyperbolicProblem<TypeEquation>& var) :
    var_time(var), mesh(var.var_harmonic.mesh), Glob_jacobian(var.Glob_jacobian)
  {
    InitDefaultValues();
  }


  template<class Dimension> 
  void CurrentSource_Vlasov<Dimension>::InitDefaultValues()
  {
    if (var_time.order_interpolation_current > 0)
      {
	evalJ.Reallocate(var_time.PointsInterpolation.GetM());
	evalRho.Reallocate(var_time.PointsInterpolation.GetM());
      }
    else
      {
	evalJ.Reallocate(var_time.PointsQuadrature.GetM());
	evalRho.Reallocate(var_time.PointsQuadrature.GetM());
      }
    
    FillZero(evalJ);     FillZero(evalRho);
    
    cut_off_Sshape = var_time.radius_influence_particle; Real_wp R = cut_off_Sshape;
    invRadius = 1.0/R; invR2 = invRadius*invRadius;
    
    // gaussian parameters
    // alpha = 10.0/(R*R); beta = alpha/pi_wp;
    
    // "Hesthaven" polynomial parameters
    alpha = var_time.alpha_influence_particle; beta = (alpha+1)/pi_wp*invR2;
    if (Dimension::dim_N == 3)
      {
	// computation of beta
	IVect coef_binom(toInteger(alpha)+1); coef_binom.Fill(1);
	for (int i = 2; i < coef_binom.GetM(); i++)
	  {
	    coef_binom(i) = 1;
	    for (int j = i-1; j > 0; j--)
	      coef_binom(j) += coef_binom(j-1);
	    
	  }
	
	Real_wp val_int(0);
	for (int i = 0; i < coef_binom.GetM(); i++)
	  {
	    if (i%2 == 0)
	      val_int += Real_wp(coef_binom(i))/Real_wp(3+2*i);
	    else
	      val_int -= Real_wp(coef_binom(i))/Real_wp(3+2*i);
	  }
	// DISP(coef_binom); DISP(val_int);
	beta = 1.0/(4.0*pi_wp*val_int*R*R*R);
	// DISP(beta);
      }
    
    xmin = var.var_harmonic.GetXmin(); xmax = var.var_harmonic.GetXmax();
    ymin = var.var_harmonic.GetYmin(); ymax = var.var_harmonic.GetYmax();
    zmin = var.var_harmonic.GetZmin(); zmax = var.var_harmonic.GetZmax();

    // PLATEAU H2 parameter
    // beta = 10.0/(3.0*pi_wp*R*R);
    
    // PLATEAU H1 parameter
    // beta = 3.0/(pi_wp*R*R);
    // alpha = beta*invRadius;
  }
  
  
  //! updating values of J for a particle
  /*!
    \param[in] poids_J weight of the particle multiplied by the charge
    \param[in] vk velocity of the particle
    \param[in] xk position of the particle
    \param[in] num_elem subdivision of the grid where the particle lies
    \param[in] pt_loc not used
    \param[in] ndiv not used
    \param[in] var considered problem
   */
  template<class Dimension>
  void CurrentSource_Vlasov<Dimension>::
  UpdateCurrent(const Real_wp& poids_J, const R_N& vk, const R_N& xk,
		int num_elem, const R_N& pt_loc)
  {
    int nb_points = var_time.CloudElem_QuadraturePoint(num_elem).GetM();
    CoefficientShape.Reallocate(nb_points);
    
    int dim_N = xk.GetM();
    //bool point_boundary = false;
    VectR_N vec_xk(1);
    // if the point if near of a periodic boundary, we generate "fictive" images of this point
    if (var_time.periodic_vector.GetM() > 0)
      if ((xk(0) < (xmin+cut_off_Sshape))||(xk(0) > (xmax-cut_off_Sshape))
	  ||(xk(1) < (ymin+cut_off_Sshape))||(xk(1) > (ymax-cut_off_Sshape))
	  ||((dim_N >= 3)&&(xk(2) < (zmin+cut_off_Sshape)))
          ||((dim_N>=3)&&(xk(2) > (zmax-cut_off_Sshape))))
	{
	  //point_boundary = true;
	  int Nvec = var_time.periodic_vector.GetM();
	  vec_xk.Reallocate(Nvec);
	  for (int m = 0; m < Nvec; m++)
	    Add(xk, var_time.periodic_vector(m), vec_xk(m));
	}
	  
    if (var_time.order_interpolation_current > 0)
      {
	// loop on all points which are supposed near
	Real_wp evalS, coefS; // DISP(nb_points);
	for (int k1 = 0; k1 < nb_points; k1++)
	  {
	    int nquad = var_time.CloudElem_QuadraturePoint(num_elem)(k1);
	    //Real_wp dist = xk.Distance(var_time.PointsInterpolation(nquad));
            Real_wp dist = 1;
            abort();
	    /*for (int m = 1; m < vec_xk.GetM(); m++)
	      dist = min(dist, vec_xk(m).Distance(var_time.PointsInterpolation(nquad)));*/
	    
	    if (dist < cut_off_Sshape)
	      {
		// new quadrature point
		evalS = GetShapeFunction(dist);
		CoefficientShape(k1) = evalS;
		
		coefS = evalS*poids_J;
	    
		// updating value of J and rho
		Add(coefS, vk, evalJ(nquad));
		evalRho(nquad) += coefS;
	      }
	    else
	      CoefficientShape(k1) = 0;
	  }
      }
    else
      {
		// loop on all points which are supposed near
	Real_wp evalS, coefS; // DISP(nb_points);
	for (int k1 = 0; k1 < nb_points; k1++)
	  {
	    int nquad = var_time.CloudElem_QuadraturePoint(num_elem)(k1);
            Real_wp dist = 1; abort();
	    /*Real_wp dist = xk.Distance(var_time.PointsQuadrature(nquad));
	    for (int m = 1; m < vec_xk.GetM(); m++)
	      dist = min(dist, vec_xk(m).Distance(var_time.PointsQuadrature(nquad)));
            */
	    if (dist < cut_off_Sshape)
	      {
		// new quadrature point
		evalS = GetShapeFunction(dist);
                abort();
		//CoefficientShape(k1) = evalS*Glob_jacobian(nquad);
		
		coefS = evalS*poids_J;
	    
		// updating value of J and rho
		Add(coefS, vk, evalJ(nquad));
		evalRho(nquad) += coefS;
	      }
	    else
	      CoefficientShape(k1) = 0;
	  }
      }
    
  }
  
  
  //! evaluation of the shape function \f$ \hat{S}(r) \f$
  /*!
    \param[in] dist value of r
   */
  template<class Dimension>
  inline Real_wp CurrentSource_Vlasov<Dimension>
  ::GetShapeFunction(const Real_wp& dist)
  {
    Real_wp Sval = 0;
    // gaussian
    // Sval = beta*exp(-alpha*dist*dist);
    
    // "Hesthaven" polynomial
    Sval = beta*pow((1-invR2*dist*dist), alpha);
    
    // PLATEAU_H2
    // Real_wp x = dist*invRadius;
    // Sval = 1.0 + x*x*(-3.0+2.0*x);
    // Sval *= beta;
    
    // PLATEAU_H1
    // Sval = beta - alpha*dist;
    
    return Sval;
  }


  /**********************
   * VlasovMaxwell_Base *
   **********************/
  
  
  //! constructor for base Vlasov-Maxwell problem
  template<class Dimension> template<class TypeEquation>
  VlasovMaxwell_Base<Dimension>
  ::VlasovMaxwell_Base(HyperbolicProblem<TypeEquation>& var)
    : laplace_solver(var_laplace), var_harmonic(var.var_harmonic), var_time(var), var_vlasov(var)
  {
    increment_particles = 100000; // size of blocks containg particle parameters
    name_file_position = "pos.dat"; // file where trajectory of first particle is stored
    relativist_particle = true; // relativistic computation
    type_correction_divergence = CORRECTION_NONE; // type of correction
    xsi_hyperbolic_correction = 1; sigma_hyperbolic = 1; // parameters for hyperbolic correction
    computation_energy = true; // computation of energy
    file_output_energy = "energy.dat"; // file where energy will be stored
    radius_influence_particle = 0.2; // parameter R for particles
    maximum_velocity = 0.9999; // used for non-relativistic simulation (for preload)
    type_space_charge_limited_emission = FULL_INTERVAL; 
    // algorithm used to create particle when E.n > Ebreakdown
    order_integration_current = 0; // order of integration
    regular_interpolation_current = true;
    nb_subdiv_interpolation = 1; 
    order_interpolation_current = 0;
    nb_points_quadrature_current = 0;
    nb_points_interpolation_current = 0;
    var_laplace.InitIndices(100);
  }
  
  
  //! Reading of an input file line
  template<class Dimension>
  void VlasovMaxwell_Base<Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    var_laplace.SetInputData(description_field, parameters);
    laplace_solver.SetInputData(description_field, parameters);
    
    if (!description_field.compare("RelativistParticle"))
      {
	// relativistic particles ?
	if (!parameters(0).compare("YES"))
	  relativist_particle = true;
	else
	  relativist_particle = false;
      }
    else if (!description_field.compare("InterpolationCurrent"))
      {
	// computation of current on interpolation points before intergation ?
	if (!parameters(0).compare("YES"))
	  {
	    order_interpolation_current = to_num<int>(parameters(1));
	    if (parameters.GetM() > 2)
	      {
		if (!parameters(2).compare("REGULAR"))
		  regular_interpolation_current = true;
		else
		  regular_interpolation_current = false;
		
		if (parameters.GetM() > 3)
		  nb_subdiv_interpolation = to_num<int>(parameters(3));
		
	      }
	  }
	else
	  order_interpolation_current = 0;
      }
    else if (!description_field.compare("MaximumVelocity"))
      {
	// only used for non-relativistic simulations
	// it is the maximum velocity allowed (to avoid velocities greater than 1)
	// put a float between 0 and 1
	maximum_velocity = to_num<Real_wp>(parameters(0));
	if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
	  maximum_velocity *= PhysicalConstant::speed_light;
      }
    else if (!description_field.compare("DivergenceCorrection"))
      {
	// technique to enforce charge conservation law
	if (!parameters(0).compare("BORIS"))
	  type_correction_divergence = CORRECTION_BORIS;
	else if (!parameters(0).compare("HYPERBOLIC"))
	  {
	    type_correction_divergence = CORRECTION_HYPERBOLIC;
	    xsi_hyperbolic_correction = to_num<Real_wp>(parameters(1));
	    sigma_hyperbolic = to_num<Real_wp>(parameters(2));
	  }
	else if (!parameters(0).compare("INITIAL"))
	  {
	    type_correction_divergence = CORRECTION_INITIAL;
	  }
	else
	  type_correction_divergence = CORRECTION_NONE;
      }
    else if (!description_field.compare("SpaceChargeLimitedEmission"))
      {
	// method to find intervals where |E.n| > E_{breakdown}
	if (!parameters(0).compare("SUBDIV"))
	  type_space_charge_limited_emission = SUBDIVISION_INTERVAL;
	else if (!parameters(0).compare("FULL"))
	  type_space_charge_limited_emission = FULL_INTERVAL;
	else
	  type_space_charge_limited_emission = FIND_INTERVAL;
      }
    else if (!description_field.compare("EnergyComputation"))
      {
	// computation of energy ?
	computation_energy = false;
	if (!parameters(0).compare("YES"))
	  {
	    computation_energy = true;
	    file_output_energy = parameters(1);
	  }
      }
    else if (!description_field.compare("RadiusParticle"))
      {
	// radius R of particles and \alpha
	radius_influence_particle = to_num<Real_wp>(parameters(0));
	alpha_influence_particle = to_num<Real_wp>(parameters(1));
      }
    else if (!description_field.compare("OrderIntegrationCurrent"))
      {
	// order of integration or interpolation for the computation of -\int J \varphi
	order_integration_current = to_num<int>(parameters(0));
      }
    else if (!description_field.compare("ParametersOutputRho"))
      {
	// density rho can be displayed
	output_rho_param.SetInputData(description_field, parameters);
	output_rho_param.SetTotalFieldFile(parameters(3));
	output_rho_param.SetDiffractedFieldFile(parameters(4));
      }
    else if (!description_field.compare("ParametersOutputParticle"))
      {
	// particles can be displayed
	output_particle_param.SetInputData(description_field, parameters);
	output_particle_param.SetComponent(to_num<int>(parameters(3)));
	output_particle_param.SetTotalFieldFile(parameters(4));
	output_particle_param.SetNature(to_num<int>(parameters(5)));
      }
    else if (!description_field.compare("InitialParticles"))
      {
	// preload of particles
	initial_particles_param.PushBack(parameters);
      }
    else if (!description_field.compare("BeamCurrent"))
      {
	// definition of a beam
	Vlasov_BeamCurrent<Dimension> Bsrc;
	Bsrc.SetInputData(description_field, parameters);
	Liste_beam_current.PushBack(Bsrc);
      }
    else if (!description_field.compare("BreakdownCurrent"))
      {
	// definition of a surface where emission of particles is allowed when E.n > Ebreakdown
	Vlasov_BeamCurrent<Dimension> Bsrc;
	Bsrc.SetInputBreakdown(description_field, parameters);
	breakdown_currents.PushBack(Bsrc);
      }

  }
  
  
  //! returns the maximal allowed number of particles without reallocation of arrays
  template<class Dimension>
  int VlasovMaxwell_Base<Dimension>::GetMaximumNumberParticles() const
  {
    int nb_max = 0;
    for (int i = 0; i < nb_max_particles.GetM(); i++)
      nb_max += nb_max_particles(i);
    
    return nb_max;
  }
  
  //! returns the number of particles in the system
  template<class Dimension>
  int VlasovMaxwell_Base<Dimension>::GetNumberParticles() const
  {
    int nb = 0;
    for (int i = 0; i < nb_particles.GetM(); i++)
      nb += nb_particles(i);
    
    return nb;
  }
  
  
  void TimesProd(const R2& u, const R2& v, R2& w)
  {
    abort();
  }
  
  
  //! initializations before time iterations
  template<class Dimension>
  void VlasovMaxwell_Base<Dimension>::InitTimeIterations()
  {
    Mesh<Dimension>& mesh = var_harmonic.mesh;
    MeshNumbering<Dimension>& mesh_num = var_harmonic.mesh_num;
    // initialization for divergence correction
    // Laplacian problem
    for (int i = 0; i < var_laplace.GetNbPhysicalIndices(); i++)
      {
	// var_laplace.ref_rho(i).Zero();
	// rho is not set to 0 to ensure well-posedness of the linear system
	var_laplace.ref_rho(i).SetConstant(1e-7);
	var_laplace.ref_mu(i).SetIdentity();
      }
    
    // var_laplace.mesh = var_harmonic.mesh;
    var_laplace.mesh_data = var_harmonic.mesh_data;
    
    // construction of the mesh and finite element
    var_laplace.ComputeMeshAndFiniteElement(string("TRIANGLE_LOBATTO"));
    var_laplace.PerformOtherInitializations();
    
    // computation and factorization of the Laplacian matrix 
    if (type_correction_divergence != CORRECTION_NONE)
      {  
	GlobalGenericMatrix<Real_wp> nat_mat;
	var_laplace.ComputeMassMatrix();
	laplace_solver.PerformFactorizationStep(nat_mat);
      }
    else
      var_laplace.ComputeMassMatrix();
    
    // changing boundary conditions
    if (type_correction_divergence == CORRECTION_HYPERBOLIC)
      for (int i = 0; i < var_laplace.mesh.GetNbReferences(); i++)
	if (var_laplace.mesh.GetBoundaryCondition(i) == BoundaryConditionEnum::LINE_DIRICHLET)
	  var_laplace.mesh.SetBoundaryCondition(i, BoundaryConditionEnum::LINE_ABSORBING);
    
    // modification of number of dofs to incorporate initial preloads
    int nb_part = GetInitialNumberParticles();
    //var_time.nodl += 4*nb_part;
    abort();
    std::remove(name_file_position.data());
        
    // initializations of outputs for rho and particles
    output_rho_param.InitTime();
    output_particle_param.InitTime();
    std::remove(file_output_energy.data());
        
    // incrementing dofs for hyperbolic correction
    if (type_correction_divergence == CORRECTION_HYPERBOLIC)
      {
	abort();
	//var_time.nodl += var_laplace.mesh_num.GetNbDof();
      }
    
    // computation of projectors
    VectR_N PtsQuadND, PtsInterpolND; VectReal_wp WeightsND;
    var_vlasov.ComputeProjectors(PtsInterpolND, PtsQuadND, WeightsND);
    
    // DISP(PtsInterpolND); DISP(PtsQuadND); DISP(WeightsND);
    WeightsQuadrature = WeightsND;
    nb_points_quadrature_current = PtsQuadND.GetM();
    nb_points_interpolation_current = PtsInterpolND.GetM();
    int nb_pts_quad = PtsQuadND.GetM();
    
    PointsQuadrature.Reallocate(mesh.GetNbElt()*nb_pts_quad);
    PointsInterpolation.Reallocate(mesh.GetNbElt()*PtsInterpolND.GetM());
    Glob_jacobian.Reallocate(mesh.GetNbElt());
    Glob_DFjm1.Reallocate(mesh.GetNbElt());
    VectR_N s;
    SetPoints<Dimension> PointsElem;
    SetMatrices<Dimension> MatricesElem;
    Glob_normale.Reallocate(mesh.GetNbBoundary());
    Glob_point_boundary.Reallocate(mesh.GetNbBoundary());
    Glob_DFjm1_Boundary.Reallocate(mesh.GetNbBoundary());
    
    inv_mass_matrix.Reallocate(var_laplace.mesh_num.GetNbDof());
    inv_mass_matrix.Fill(0);
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	const TypeElement& Fb = var_harmonic.GetReferenceElement(i);
	var_laplace.mesh.GetVerticesElement(i, s);
	Fb.FjElem(s, PointsElem, var_laplace.mesh, i);
	Fb.DFjElem(s, PointsElem, MatricesElem, var_laplace.mesh, i);
	
	MatrixN_N mat_dfj, dfjm1; Real_wp jacob; R_N pt_glob;
        Glob_jacobian(i).Reallocate(nb_pts_quad);
        Glob_DFjm1(i).Reallocate(nb_pts_quad);
	if (this->order_integration_current == 0)
	  for (int j = 0; j < nb_pts_quad; j++)
	    {
	      R_N x = PointsElem.GetPointQuadrature(j); 
	      PointsQuadrature(i*nb_pts_quad+j) = x;
	      jacob = Det(MatricesElem.GetPointQuadrature(j));
	      Glob_jacobian(i)(j) = jacob*Fb.WeightsND(j);
	    }
	
	for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
	  {
	    jacob = Det(MatricesElem.GetPointQuadrature(j));
	    inv_mass_matrix(var_laplace.mesh_num.Element(i).GetNumberDof(j))
              += jacob*Fb.WeightsND(j);
	  }
	  
	if (this->order_integration_current > 0)
	  for (int j = 0; j < nb_pts_quad; j++)
	    {
	      Fb.Fj(s, PointsElem, PtsQuadND(j), pt_glob, var_laplace.mesh, i);
	      Fb.DFj(s, PointsElem, PtsQuadND(j), mat_dfj, var_laplace.mesh, i);
	      PointsQuadrature(i*nb_pts_quad+j) = pt_glob;
	      jacob = Det(mat_dfj); GetInverse(mat_dfj, dfjm1);
	      Glob_jacobian(i)(j) = jacob*WeightsND(j);
	      Glob_DFjm1(i)(j) = dfjm1; Mlt(jacob, Glob_DFjm1(i)(j));
	    }
	
	if (this->order_interpolation_current > 0)
	  {
	    for (int j = 0; j < PtsInterpolND.GetM(); j++)
	      {
		Fb.Fj(s, PointsElem, PtsInterpolND(j), pt_glob, var_laplace.mesh, i);
		PointsInterpolation(i*PtsInterpolND.GetM()+j) = pt_glob;
	      }
	  }
 
	R_N normale; Real_wp dsj;
	for (int num_loc = 0; num_loc < var_laplace.mesh.Element(i).GetNbBoundary(); num_loc++)
	  {
	    int num_edge = var_laplace.mesh.Element(i).numBoundary(num_loc);
            int nb_pts_quad_bound = var_laplace.mesh_num.GetNbPointsQuadratureBoundary(num_edge);
            Glob_point_boundary(num_edge).Reallocate(nb_pts_quad_bound);
            Glob_normale(num_edge).Reallocate(nb_pts_quad_bound);
            Glob_DFjm1_Boundary(num_edge).Reallocate(nb_pts_quad_bound);
            Fb.FjSurfaceElem(s, PointsElem, var_laplace.mesh, i, num_loc);
	    Fb.DFjSurfaceElem(s, PointsElem, MatricesElem, var_laplace.mesh, i, num_loc);
	    if (var_laplace.mesh.Boundary(num_edge).numElement(0) == i)
	      {
		for (int k = 0; k < nb_pts_quad_bound; k++)
		  {
		    normale = MatricesElem.GetNormaleQuadratureBoundary(k);
		    dsj = MatricesElem.GetDsQuadratureBoundary(k);
		    GetInverse(MatricesElem.GetPointQuadratureBoundary(k), dfjm1);
		    
		    Glob_point_boundary(num_edge)(k) = PointsElem.GetPointQuadratureBoundary(k);
		    Glob_normale(num_edge)(k) = normale;
		    Glob_DFjm1_Boundary(num_edge)(k) = dfjm1;
		    Mlt(dsj, Glob_normale(num_edge)(k));
		  }
	      }
	  }
      }
    
    // DISP(PointsQuadrature);
    // mass matrix for nodal mesh
    for (int i = 0; i < inv_mass_matrix.GetM(); i++)
      inv_mass_matrix(i) = 1.0/inv_mass_matrix(i);
    
    // now precomputing clouds for each quadrature point (or interpolation point)
    abort();
    /*if (order_interpolation_current > 0)
      var_time.PreComputeCloudsGrid(PointsInterpolation, this->radius_influence_particle,
                                    CloudElem_QuadraturePoint, loc_particle);
    else
      var_time.PreComputeCloudsGrid(PointsQuadrature, this->radius_influence_particle,
                                    CloudElem_QuadraturePoint, loc_particle);
    */
    
    // initialization for beams
    for (int nbeam = 0; nbeam < Liste_beam_current.GetM(); nbeam++)
      {
	R_N v0 = Liste_beam_current(nbeam).v0;
	// true velocity
	R_N vk;
	this->GetVelocity_FromMomentum(v0, vk);
	// computing delta_iter
	if (Liste_beam_current(nbeam).delta_iter <= 0)
	  {
	    int n = toInteger(ceil(1.0/(Liste_beam_current(nbeam).nb_part_per_unit
                                        *Norm2(vk)*var_time.deltat)));
	    Liste_beam_current(nbeam).delta_iter = max(1,n); 
	  }
	// DISP(nbeam); DISP(Liste_beam_current(nbeam).delta_iter);
	
	// shifting the beam so that created particles are far from the domain of R
	// this is useful so that there is no need of a correction for the electric field
	R_N normale, t1, t2;
	Subtract(Liste_beam_current(nbeam).ptA, Liste_beam_current(nbeam).ptB, t1);
	Subtract(Liste_beam_current(nbeam).ptC, Liste_beam_current(nbeam).ptD, t2);
	if (normale.GetM() == 2)
	  {
	    normale(0) = t1(1);
	    normale(1) = -t1(0);
	  }
	else
	  TimesProd(t1, t2, normale);
	
	Mlt(1.0/Norm2(normale), normale);
	
	if (DotProd(normale, v0) < 0)
	  Mlt(-1.0, normale);
	
	Real_wp R = 0.999*this->radius_influence_particle;
	Add(-R, normale, Liste_beam_current(nbeam).ptA);
	Add(-R, normale, Liste_beam_current(nbeam).ptB);
	Add(-R, normale, Liste_beam_current(nbeam).ptC);
	Add(-R, normale, Liste_beam_current(nbeam).ptD);
	Add(Liste_beam_current(nbeam).ptA, Liste_beam_current(nbeam).ptB,
            Liste_beam_current(nbeam).center);
	Mlt(0.5, Liste_beam_current(nbeam).center);
	//DISP(normale); DISP(R); DISP(Liste_beam_current(nbeam).ptA);
        // DISP(Liste_beam_current(nbeam).ptB);
	// DISP( Liste_beam_current(nbeam).ptC); DISP( Liste_beam_current(nbeam).ptD);
      }
    
    for (int nbeam = 0; nbeam < breakdown_currents.GetM(); nbeam++)
      {
	// computing delta_iter
	// int n = toInteger(ceil(1.0/(breakdown_currents(nbeam).nb_part_per_unit
        // *GlobalVariables::speed_light*this->deltat)));
	// breakdown_currents(nbeam).delta_iter = max(1,n);
	if (breakdown_currents(nbeam).delta_iter <= 0)
	  breakdown_currents(nbeam).delta_iter = 1;
      } 

    for (int i = 0; i < var_harmonic.output_grid_param.GetM(); i++)
      {
	string file_name = var_harmonic.output_grid_param(i).GetTotalFieldFile();	
	var_harmonic.output_grid_param(i).SetFileName(2, GetBaseString(file_name) + "_H.dat");
      }
  }
  
  
  //! creation of a new species
  template<class Dimension>
  void VlasovMaxwell_Base<Dimension>
  ::CreateNewSpecies(const Real_wp& q0, const Real_wp& ratio_q0m0, VectReal_wp& Y)
  {
    nb_particles.PushBack(0);
    ratio_qm_species.PushBack(ratio_q0m0); charge_species.PushBack(q0);
    nb_max_particles.PushBack(increment_particles);
    VectReal_wp zero_vec(increment_particles); zero_vec.Fill(0);
    weight_particle.PushBack(zero_vec);
    
    int dim_N = Dimension::dim_N;
    Y.Resize(Y.GetM()+2*dim_N*increment_particles);
  }
  
  
  //! computing momentum from speed
  template<class Dimension>
  inline void VlasovMaxwell_Base<Dimension>
  ::GetMomentum_FromVelocity(const R_N& speed, R_N& momentum) const
  {
    momentum = speed;
    if (relativist_particle)
      {
	Real_wp gamma(1);
	if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
	  {
	    Real_wp vk = Norm2(speed);
	    gamma /= sqrt(1-vk*vk*PhysicalConstant::invSquare_speed_light);
	  }
	else
	  {
	    Real_wp vk = Norm2(speed);
	    gamma /= sqrt(1-vk*vk);
	  }
	momentum *= gamma;
      }
  }
  
  //! computing speed from momentum
  template<class Dimension>
  inline void VlasovMaxwell_Base<Dimension>
  ::GetVelocity_FromMomentum(const R_N& pk, R_N& vk) const
  {
    if (relativist_particle)
      {
	Real_wp pk2, gamma2, gamma;
	if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
	  {
	    pk2 = DotProd(pk, pk);
	    gamma2 = 1.0 + pk2*PhysicalConstant::invSquare_speed_light;
	    gamma = sqrt(gamma2);
	    vk = pk; vk *= 1.0/gamma;
	  }
	else
	  {
	    pk2 = DotProd(pk, pk);
	    gamma2 = 1.0 + pk2;
	    gamma = sqrt(gamma2);
	    vk = pk; vk *= 1.0/gamma;
	  }
      }
    else
      {
	vk = pk;
      }
  }
  
  //! computing momentum from speed
  template<class Dimension>
  inline void VlasovMaxwell_Base<Dimension>
  ::GetVelocity_FromMomentum(const R_N& pk, R_N& vk, Real_wp& gamma, Real_wp& gamma2) const
  {
    if (relativist_particle)
      {
	Real_wp pk2;
	if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
	  {
	    pk2 = DotProd(pk, pk);
	    gamma2 = 1.0 + pk2*PhysicalConstant::invSquare_speed_light;
	    gamma = sqrt(gamma2);
	    vk = pk; vk *= 1.0/gamma;
	  }
	else
	  {
	    pk2 = DotProd(pk, pk);
	    gamma2 = 1.0 + pk2;
	    gamma = sqrt(gamma2);
	    vk = pk; vk *= 1.0/gamma;
	  }
      }
    else
      {
	gamma = 1.0; gamma2 = 1.0;
	vk = pk;
      }
  }
  
  //! computing speed from momentum
  template<class Dimension>
  int VlasovMaxwell_Base<Dimension>::GetInitialNumberParticles() const
  {
    int nb_part = 0;
    for (int num_p = 0; num_p < initial_particles_param.GetM(); num_p++)
      {
	if (!initial_particles_param(num_p)(0).compare("ONE"))
	  nb_part += increment_particles;
	else if (!initial_particles_param(num_p)(0).compare("CIRCLE"))
	  nb_part += increment_particles;
	else
	  {
	    int Ntot = to_num<int>(initial_particles_param(num_p)(3));
	    int nb_div = (Ntot-1)/increment_particles + 2;
	    DISP(Ntot); DISP(nb_div);
	    nb_part += increment_particles*nb_div;
	  }
      }
    
    DISP(nb_part);
    return nb_part;
  }
  
}

#define MONTJOIE_FILE_VLASOV_MAXWELL_CXX
#endif
