#ifndef MONTJOIE_FILE_DEFINE_SOURCE_ELLIPTIC_CXX

namespace Montjoie
{
  
  /*************************
   * VarSourceProblem_Base *
   *************************/
  
  
  //! kind of default constructor
  void VarSourceProblem_Base::InitDefaultValues()
  {
    // by default : a null source
    type_source.Reallocate(1);
    incident_source.Reallocate(1);
    source_space_param.Reallocate(1);
    
    type_source(0) = SRC_NULL;
    incident_source(0) = INCIDENT_PLANE_WAVE;

    threshold_rhs = epsilon_machine;
  }
  
  
  //! modifies object with a line of the data file
  void VarSourceProblem_Base::SetInputData(const string& description_field,
					   const VectString& parameters)
  {
    if (!description_field.compare("TypeSource"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarSourceProblem" << endl;
	    cout << "TypeSource needs more parameters, for instance :" << endl;
	    cout << "TypeSource = SRC_DIRAC" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	this->UpdateSourceParam(parameters, this->type_source, this->incident_source);
      }
    else if (!description_field.compare("ThresholdRhs"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarSourceProblem" << endl;
	    cout << "ThresholdRhs needs more parameters, for instance :" << endl;
	    cout << "ThresholdRhs = epsilon" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	threshold_rhs = to_num<Real_wp>(parameters(0));
      }    
  }


  //! sets source parameters
  void VarSourceProblem_Base::SetParameterSource(const Vector<VectString>& param)
  {
    type_source.Reallocate(1);
    incident_source.Reallocate(1);
    source_space_param.Clear();
    source_space_param.Reallocate(1);
    
    for (int k = 0; k < param.GetM(); k++)
      this->UpdateSourceParam(param(k), this->type_source, this->incident_source);
  }
  
  
  //! returns the type of the incident field as given in name
  int VarSourceProblem_Base::GetIncidentFieldType(const string& name)
  {
    if (!name.compare("GAUSSIAN_BEAM"))
      return INCIDENT_GAUSSIAN_BEAM;
    else if (!name.compare("PLANE_WAVE"))
      return INCIDENT_PLANE_WAVE;
   else if (!name.compare("PLANE_WAVE_CPLX"))
      return INCIDENT_PLANE_WAVE_CPLX;
    else if (!name.compare("HANKEL_WAVE"))
      return INCIDENT_HANKEL;
    else if (!name.compare("LAYERED_PLANE_WAVE"))
      return INCIDENT_LAYERED_PLANE_WAVE;	      

    return INCIDENT_NONE;
  }


  //! updates the type of source from a line of the data file
  void VarSourceProblem_Base
  ::UpdateSourceParam(const VectString& parameters,
		      Vector<int>& type, Vector<int>& incident)
  {
    // we retrieve the source number
    int n = 0, nb = 0;
    if (isdigit(parameters(0)[0]))
      {
	n = to_num<int>(parameters(0));
	nb++;
      }
    
    // we resize the arrays if n is too large
    if (n >= type.GetM())
      {
	type.Resize(n+1);
	incident.Resize(n+1);
	source_space_param.Resize(n+1);
      }
    
    // parameters are stored in source_space_param
    // source number is not stored
    VectString param;
    param.Reallocate(parameters.GetM()-nb);
    for (int i = nb; i < parameters.GetM(); i++)
      param(i-nb) = parameters(i);
    
    source_space_param(n).PushBack(param);
    
    // if scattering = true, a scattering problem is solved
    // => need of defining the incident wave
    bool scattering = false;

    // SRC_SURFACE, SRC_VOLUME and SRC_MODE induce the same source SRC_VOLUME
    // since VolumetricSource object implements a generic source
    // where you can specify different volume and surface integrals
    if (!parameters(nb).compare("SRC_VOLUME")) 
      type(n) = SRC_VOLUME;
    else if (!parameters(nb).compare("SRC_DIRAC"))
      type(n) = SRC_DIRAC;
    else if (!parameters(nb).compare("SRC_SURFACE")) 
      type(n) = SRC_VOLUME;
    else if (!parameters(nb).compare("SRC_USER")) 
      type(n) = SRC_USER;
    else if (!parameters(nb).compare("SRC_MODE")) 
      type(n) = SRC_VOLUME;
    else if (!parameters(nb).compare("SRC_TOTAL_FIELD")) 
      {
	type(n) = SRC_TOTAL_FIELD;
	scattering = true;
      }
    else if (!parameters(nb).compare("SRC_DIFFRACTED_FIELD"))
      {
	type(n) = SRC_DIFFRACTED_FIELD;
	scattering = true;
      }
    else
      type(n) = SRC_NULL;
    
    if (scattering)
      {
	if (parameters.GetM() <= nb+1)
	  {
	    cout << "In SetInputData of VarSourceProblem" << endl;
	    cout << "TypeSource needs more parameters, for instance :" << endl;
	    cout << "TypeSource = SRC_DIFFRACTED_FIELD PLANE_WAVE" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	incident(n) = this->GetIncidentFieldType(parameters(nb+1));
      }
  }
  

  //! reads parameters of a gaussian
  template<class Dimension>
  void VarSourceProblem_Base
  ::InitGaussianParameter(GaussianSource<Dimension>& f, const VectString& param, int& nb) const
  {
    typename Dimension::R_N origin;
    Real_wp RadiusSource(0), RadiusSourceCutOff(0);
    
    if (param.GetM() <= nb+Dimension::dim_N)
      {
	cout << "Enter more parameters for the gaussian" << endl;
	cout << "Current parameters are " << endl << param << endl;
	abort();
      }
    
    for (int k = 0; k < Dimension::dim_N; k++)
      origin(k) = to_num<Real_wp>(param(nb++));
    
    RadiusSource = to_num<Real_wp>(param(nb++));
    
    if (param.GetM() == nb)
      RadiusSourceCutOff = RadiusSource;
    else
      {
	if (!isalpha(param(nb)[0]))
	  RadiusSourceCutOff = to_num<Real_wp>(param(nb++));
      }
    
    f.Init(origin, RadiusSource, RadiusSourceCutOff);
  }
	    

  //! reads parameters of a gaussian
  template<class Dimension>
  void VarSourceProblem_Base
  ::InitRandomGaussianParameter(GaussianSource<Dimension>& f, const VectString& param, int nb) const
  {
    typename Dimension::R_N origin;
    Real_wp RadiusSource(0), RadiusSourceCutOff(0);
    
    srand(Seed);
    for (int k = 0; k < Dimension::dim_N; k++)
      {
	Real_wp random_number = Real_wp(rand())/Real_wp(RAND_MAX);
	Real_wp minimal_value = to_num<Real_wp>(param(nb++));
	Real_wp maximal_value = to_num<Real_wp>(param(nb++));
	origin(k) = (maximal_value - minimal_value)*random_number
	  + minimal_value;
      }
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    // We write directly the origins in the file origins.txt
    if (rank_proc == 0)
      {
	ofstream origines ("origins.txt", ios::app);
	for (int k = 0; k < Dimension::dim_N; k++)
	  origines << origin(k) << " ";
	
	origines << endl;
	origines.close();
      }
    
    RadiusSource = to_num<Real_wp>(param(nb++));
    if (param.GetM() == nb)
      RadiusSourceCutOff = RadiusSource;
    else
      {
	if (!isalpha(param(nb)[0]))
	  RadiusSourceCutOff = to_num<Real_wp>(param(nb++));
      }
    
    f.Init(origin, RadiusSource, RadiusSourceCutOff);
  }


  //! initialisation of ffts for outputs of solutions, when we know their modes
  void VarSourceProblem_Base
  ::InitFftComputation(FftInterface<Complex_wp>& fft_interface) const
  {
    int nx, ny, nz; 
    bool sym_teta;
    var_boundary.GetPeriodicNumberModes(nx, ny, nz, sym_teta);
    
    if (nz <= 1)
      {
        if (ny <= 1)
          {
            if (nx <= 1)
              {
                // no fft to compute
              }
            else
              {
                // 1-D fft
                fft_interface.Init(nx);
              }
          }
        else
          {
            // 2-D fft
            fft_interface.Init(nx, ny);
          }
      }
    else
      {
        // 3-D fft
        fft_interface.Init(nx, ny, nz);
      }
  }


  //! copies input parameters of another similar problem
  void VarSourceProblem_Base::CopyInputData(const VarSourceProblem_Base& var)
  {
    type_source = var.type_source;
    incident_source = var.incident_source;
    threshold_rhs = var.threshold_rhs;
    source_space_param = var.source_space_param;
  }    
  
  
  /************************
   * VarSourceProblem_Dim *
   ************************/


  //! Retrieves parameters of an incident wave (wave vector, origin of phase, pulsation and waist)
  template<class Dimension>
  void VarSourceProblem_Dim<Dimension>
  ::ReadIncidentWaveParam(const VectString& parameters, R_N& kwave,
                          TinyVector<Complex_wp, Dimension::dim_N>& kc,
                          R_N& origin, Real_wp& omega, Real_wp& w) const
  {
    kwave = var_problem.GetWaveVector();
    origin = var_problem.GetPhaseOrigin();
    omega = var_problem.GetOmega();
    kc = kwave;

    // loop over parameters to find a different wave vector and/or origin
    Real_wp theta(0), phi(0); Complex_wp kx(0,0), ky(0,0);
    bool new_incidence = false;
    for (int k = 0; k < parameters.GetM(); k++)
      {
	if (parameters(k) == "Theta")
	  {
	    new_incidence = true;
	    theta = to_num<Real_wp>(parameters(k+1))*pi_wp/180;
	  }
	else if (parameters(k) == "Phi")
	  {
	    new_incidence = true;
	    phi = to_num<Real_wp>(parameters(k+1))*pi_wp/180;
	  }
        else if (parameters(k) == "Kx")
          kx = to_num<Complex_wp>(parameters(k+1));
        else if (parameters(k) == "Ky")
          ky = to_num<Complex_wp>(parameters(k+1));
	else if (parameters(k) == "Origin")
	  {
	    for (int p = 0; p < Dimension::dim_N; p++)
	      origin(p) = to_num<Real_wp>(parameters(k+p+1));
	  }
	else if (parameters(k) == "Waist")
	  w = to_num<Real_wp>(parameters(k+1));
	else if (parameters(k) == "Pulsation")
	  omega = to_num<Real_wp>(parameters(k+1));
      }
    
    kc(0) = kx; kc(1) = ky;
    if (Dimension::dim_N == 3)
      kc(2) = sqrt(omega*omega - kx*kx - ky*ky);
    else
      kc(1) = sqrt(omega*omega - kx*kx);

    if (new_incidence)
      SetIncidentAngle(omega, kwave, theta, phi);
  }
  
  
  //! allocation of incident field
  template<class Dimension>
  IncidentWaveField<Complex_wp, Dimension>* VarSourceProblem_Dim<Dimension>
  ::GetNewIncidentField(int n, const Vector<VectString>& param, const Complex_wp& val) const
  {
    // parameters of the incident wave
    R_N kwave, origin;
    TinyVector<Complex_wp, Dimension::dim_N> kc;
    Real_wp omega, w(0);
    
    ReadIncidentWaveParam(param(0), kwave, kc, origin, omega, w);
    
    // coef a, b at infinity for layered plane wave
    Real_wp a_infty = 1.0;
    Real_wp b_infty = 1.0;
    var_source.GetCoefAB_Infinity(a_infty, b_infty);

    // incident field is constructed and returned
    switch (n)
      {
      case VarSourceProblem_Base::INCIDENT_PLANE_WAVE:
	return new PlaneWaveIncidentField<Complex_wp, Dimension>(origin, kwave);
      case VarSourceProblem_Base::INCIDENT_PLANE_WAVE_CPLX:
	return new PlaneWaveIncidentFieldComplex<Complex_wp, Dimension>(origin, kc);
      case VarSourceProblem_Base::INCIDENT_GAUSSIAN_BEAM:
	return new GaussianBeamIncidentField<Complex_wp, Dimension>(origin, kwave, w);
      case VarSourceProblem_Base::INCIDENT_HANKEL:
	return new HankelIncidentField<Complex_wp, Dimension>(origin, kwave);
      case VarSourceProblem_Base::INCIDENT_LAYERED_PLANE_WAVE:
	return new LayeredPlaneWaveIncidentField<Complex_wp, Dimension>(origin, kwave, omega, 
									a_infty, b_infty, param(0));
      default :
	cout << "Unknown incident field" << endl;
	abort();
      }
  }

  
  //! allocation of incident field
  template<class Dimension>
  IncidentWaveField<Real_wp, Dimension>* VarSourceProblem_Dim<Dimension>
  ::GetNewIncidentField(int n, const Vector<VectString>& param, const Real_wp& val) const
  {
    R_N kwave, origin;
    TinyVector<Complex_wp, Dimension::dim_N> kc;
    Real_wp omega, w;
    
    ReadIncidentWaveParam(param(0), kwave, kc, origin, omega, w);
    
    // only plane waves are implemented in time-domain simulations
    return new PlaneWaveIncidentField<Real_wp, Dimension>(origin, kwave);
  }
  

  /*************************
   * VarSourceProblem_Cplx *
   *************************/


  //! Destructor
  template<class Complexe, class Dimension>
  VarSourceProblem_Cplx<Complexe, Dimension>::~VarSourceProblem_Cplx()
  {
    ClearIncidentArray();
  }
  
  
  //! Clears memory used by incident field 
  template<class Complexe, class Dimension>
  void VarSourceProblem_Cplx<Complexe, Dimension>::ClearIncidentArray()
  {
    for (int n = 0; n < incident_wave_proj.GetM(); n++)
      {
	if (incident_wave_proj(n) != NULL)
	  delete incident_wave_proj(n);

	if (incident_wave_field(n) != NULL)
	  delete incident_wave_field(n);
      }
    
    incident_wave_proj.Clear();
    incident_wave_field.Clear();
  }

  
  //! returns a pointer to the incident field as used to compute the source
  template<class Complexe, class Dimension>
  IncidentWaveField<Complexe, Dimension>* VarSourceProblem_Cplx<Complexe, Dimension>
  ::GetIncidentFieldGen(int n, const Complexe&) const
  {
    if (incident_wave_field(n) == NULL)
      {
	cout << "Unallocated incident field" << endl;
	abort();
      }
    
    return incident_wave_field(n);
  }


  //! returns a pointer to the incident field projector as used in OutputHarmonic.cxx
  template<class Complexe, class Dimension>
  IncidentWaveProjector<Complexe, Dimension>* VarSourceProblem_Cplx<Complexe, Dimension>
  ::GetIncidentWaveProjectorGen(int n, const Complexe&) const
  {
    if (incident_wave_proj(n) == NULL)
      {
	cout << "Unallocated projector" << endl;
	abort();
      }
    
    return incident_wave_proj(n);
  }
  

  //! returns the volume source f as given by the user in the data file
  template<class Complexe, class Dimension>
  void VarSourceProblem_Cplx<Complexe, Dimension>
  ::ConstructVolumeSourceFunctions(const Vector<VectString>& param,
                                   VolumetricSource_Base<Complexe, Dimension>& var) const
  {
    // initialisation with a null source
    VirtualSourceField<Complexe, Dimension>* fsrc;
    fsrc = NULL;

    Vector<Complexe> polar;
    
    // loop over lines defining the source
    for (int num = 0; num < param.GetM(); num++)
      {
	if ((param(num).GetM() > 1) && (param(num)(0) == "SRC_VOLUME"))
	  {
	    var_problem.GetPolarization(polar);
            int position_word = -1;
            // we look for the first parameter, which begins by a letter
            // This parameter is the description of the volume source
            for (int i = 1; i < param(num).GetM(); i++)
              if (isalpha(param(num)(i)[0]))
                {
                  position_word = i;
                  break;
                }
            
            IVect ref;
            if (position_word != -1)
              {
                // references
                ref.Reallocate(position_word-1);
                for (int i = 0; i < position_word-1; i++)
                  ref(i) = to_num<int>(param(num)(i+1));
              }
            
	    int nb = position_word;
            fsrc = GetNewVolumeSourceFunction(ref, param(num), nb, polar, var);
            if (fsrc != NULL)
              {
                if (ref.GetM() == 0)
                  var.SetVolumeSourceFunction(*fsrc);
                else
                  var.SetVolumeSource(ref, fsrc);
              }
	  }
      }
  }


  //! returns the volume source g as given by the user in the data file
  template<class Complexe, class Dimension>
  VirtualSourceField<Complexe, Dimension>* VarSourceProblem_Cplx<Complexe, Dimension>
  ::GetNewVolumeSourceFunction(const IVect& ref_vol, const VectString& param,
                               int& nb, Vector<Complexe>& polar,
                               VolumetricSource_Base<Complexe, Dimension>& var) const
  {
    for (int k = nb; k < param.GetM(); k++)
      if (param(k) == "Polarization")
        {
          for (int p = 0; p < polar.GetM(); p++)
            polar(p) = to_num<Complexe>(param(k+p+1));
          
          break;
        }
    
    if (param(nb) == "GAUSSIAN")
      {
        nb++;
        GaussianSourceField<Complexe, Dimension>* fgauss;
        fgauss = new GaussianSourceField<Complexe, Dimension>();
        var_source.InitGaussianParameter(*fgauss, param, nb);
        fgauss->SetPolarization(polar);
        return fgauss;
      }
    else if (param(nb) == "RANDOM_GAUSSIAN")
      {
        nb++;
        GaussianSourceField<Complexe, Dimension>* fgauss;
        fgauss = new GaussianSourceField<Complexe, Dimension>();
        var_source.InitRandomGaussianParameter(*fgauss, param, nb);
        fgauss->SetPolarization(polar);
        return fgauss;
      }
    else if (param(nb) == "UNIFORM")
      {
        nb++;
        UniformSourceField<Complexe, Dimension>* funif;
        funif = new UniformSourceField<Complexe, Dimension>();
        funif->SetPolarization(polar);
        return funif;
      }
    else if ((param(nb) == "VARIABLE_BINARY") || (param(nb) == "GRADIENT_BINARY"))
      {
	nb++;
	if (param.GetM() <= nb)
	  {
	    cout << "Waiting one file name after VARIABLE_BINARY " << endl;
	    cout << "Example : TypeSource = SRC_VOLUME 1 VARIABLE_BINARY value.dat" << endl;
	    abort();
	  }

        bool source_grad = false;        
        if (param(nb-1) == "GRADIENT_BINARY")
          source_grad = true;

        if (source_grad)
          var.SetVariableGradientSource(ref_vol);
        else
          var.SetVariableSource(ref_vol);

        string file_name = param(nb++);
        string root = GetBaseString(file_name);
        Vector<Complexe> val;
        int rank_proc = var_problem.GetRankProcMode();
        int nb_proc = var_problem.GetNbProcPerMode();
        if (nb_proc > 1)
          file_name = root + "_P" + to_str(rank_proc) + ".dat";
        
        val.Read(file_name);
        
        int nb_comp = var_problem.GetNbComponentsAll();
        if (source_grad)
          nb_comp = var_problem.GetNbComponentsGradientAll();

        if (source_grad)
          var.evalG.Reallocate(nb_comp);
        else
          var.evalS.Reallocate(nb_comp);
        
        int num = 0;
        for (int k = 0; k < nb_comp; k++)
          {
            int nb_elt = var_problem.mesh.GetNbElt();
            Vector<Vector<Complexe> >* eval_ptr;
            if (source_grad)
              eval_ptr = &var.evalG(k);
            else
              eval_ptr = &var.evalS(k);

            Vector<Vector<Complexe> >& eval = *eval_ptr;
            eval.Reallocate(nb_elt);
            for (int i = 0; i < nb_elt; i++)
              {
                int nb_quad = var_problem.GetNbPointsQuadratureInside(i);
                if (val.GetM() < num+nb_quad)
                  {
                    cout << "Invalid source for element " << i << " component " << k << endl;
                    cout << "Size equal to " << val.GetM() << " < " << num+nb_quad << endl;
                    abort();
                  }
                
                eval(i).Reallocate(nb_quad);
                for (int j = 0; j < nb_quad; j++)
                  eval(i)(j) = val(num++);
              }
          }
        
        return NULL;
      }
  
    return NULL;
  }
  
  
  //! sets the surface source
  template<class Complexe, class Dimension>
  void VarSourceProblem_Cplx<Complexe, Dimension>
  ::ConstructSurfaceSourceFunctions(const Vector<VectString>& param,
				    VolumetricSource_Base<Complexe, Dimension>& var) const
  {
    VirtualSourceField<Complexe, Dimension>* fsrc;
    
    Vector<Complexe> polar(var_problem.GetNbComponentsAll());
    
    // loop over lines beginning with TypeSource
    int num_mode = 0;
    for (int num = 0; num < param.GetM(); num++)
      {
	var_problem.GetPolarization(polar);
	fsrc = NULL;

	int position_word = -1;
	// we look for the first parameter, which begins by a letter
	// This parameter is the description of the surface source
	for (int i = 1; i < param(num).GetM(); i++)
	  if (isalpha(param(num)(i)[0]))
	    {
	      position_word = i;
	      break;
	    }
	
	IVect ref;
	if (position_word != -1)
	  {
	    // references
	    ref.Reallocate(position_word-1);
	    for (int i = 0; i < position_word-1; i++)
	      ref(i) = to_num<int>(param(num)(i+1));
	  }
	
	if ((param(num).GetM() > 1) && (param(num)(0) == "SRC_SURFACE"))
	  {
	    int nb = position_word;
	    fsrc = GetNewSurfaceSourceFunction(ref, param(num), nb, polar, var);

	    if (fsrc != NULL)
	      var.SetSurfaceSource(ref, fsrc);
	  }

	if ((param(num).GetM() > 1) && (param(num)(0) == "SRC_SURFACE_GRAD"))
	  {
	    int nb = position_word;
	    fsrc = GetNewSurfaceSourceFunction(ref, param(num), nb, polar, var);

	    if (fsrc != NULL)
	      var.SetSurfaceSourceGrad(ref, fsrc);
	  }

	if ((param(num).GetM() > 1) && (param(num)(0) == "SRC_MODE"))
	  {
	    if (position_word == -1)
	      {
		cout << "Usage: TypeSource = SRC_MODE ref KEYWORD" << endl;
		cout << "Parameters are " << endl << param << endl;
		abort();
	      }
	    
	    var.SetModalSource(ref);
	    
	    ModalSourceBoundary_Dim<Complexe, Dimension>* f_mode;
	    f_mode = this->GetNewModalSourceEquation();
	    VectString parameters(param(num).GetM()-position_word);
	    for (int k = position_word; k < param(num).GetM(); k++)
	      parameters(k-position_word) = param(num)(k);
	    
	    f_mode->ComputeMode(num_mode, ref, parameters, var.evalSurf);
	    num_mode++;
	    delete f_mode;
	  }
      }
  }
  

  template<class Complexe, class Dimension>
  VirtualSourceField<Complexe, Dimension>* VarSourceProblem_Cplx<Complexe, Dimension>
  ::GetNewSurfaceSourceFunction(const IVect& ref_surf, const VectString& param,
				int& nb, Vector<Complexe>& polar,
				VolumetricSource_Base<Complexe, Dimension>& var) const
  {
    for (int k = nb; k < param.GetM(); k++)
      if (param(k) == "Polarization")
	{
	  for (int p = 0; p < polar.GetM(); p++)
	    polar(p) = to_num<Complexe>(param(k+p+1));
		  
	  break;
	}
    
    if (param(nb) == "GAUSSIAN")
      {
	nb++;
	GaussianSourceField<Complexe, Dimension>* fgauss;
	fgauss = new GaussianSourceField<Complexe, Dimension>();
	var_source.InitGaussianParameter(*fgauss, param, nb);
	fgauss->SetPolarization(polar);
	return fgauss;
      }
    else if (param(nb) == "UNIFORM")
      {
	nb++;
	UniformSourceField<Complexe, Dimension>* f_unif;
	f_unif = new UniformSourceField<Complexe, Dimension>();
	f_unif->SetPolarization(polar);
	return f_unif;
      }
    else if (param(nb) == "INCIDENT_WAVE")
      {
	nb++;
	Vector<VectString> parameters(1);
	parameters(0) = param;
	
	int type = var_source.GetIncidentFieldType(param(nb++));
	IncidentWaveField<Complexe, Dimension>* u_inc;
	u_inc = this->GetNewIncidentField(type, parameters, Complexe(0));
	
	IncidentWaveSourceField<Complexe, Dimension>* f_inc;
	f_inc = new IncidentWaveSourceField<Complexe, Dimension>(u_inc);
	f_inc->SetPolarization(polar);
	
	return f_inc;
      }
    else if (param(nb) == "VARIABLE")
      {
	nb++;
	if (param.GetM() <= nb+1)
	  {
	    cout << "Waiting two file names after VARIABLE " << endl;
	    cout << "Example : TypeSource = SRC_SURFACE 1 VARIABLE Points.txt value.dat" << endl;
	    abort();
	  }
	
	typename Dimension::VectR_N Points;
	Points.ReadText(param(nb++));
		
	Matrix<Complexe> val;
	val.Read(param(nb++));
		
	var.SetModalSource(ref_surf);
	FillVariableSource(ref_surf, Points, val, var.evalSurf);
	
	return NULL;
      }
    
    return NULL;
  }	    

  
  //! Fills eval from values contained in val associated with Points
  template<class Complexe, class Dimension>
  void VarSourceProblem_Cplx<Complexe, Dimension>
  ::FillVariableSource(const IVect& ref, typename Dimension::VectR_N& Points,
		       Matrix<Complexe>& val, Vector<Vector<Vector<Complexe> > >& eval) const
  {
    const Mesh<Dimension>& mesh = var_problem.mesh;
    // ref_cond is constructed (ref_cond(ref) = 1)
    IVect ref_cond(mesh.GetNbReferences()+1);
    ref_cond.Zero();
    for (int i = 0; i < ref.GetM(); i++)
      ref_cond(ref(i)) = 1;

    // loop over edges/faces with the required references
    typename Dimension::DimensionBoundary::VectR_N LocalPoints;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      if (ref_cond(mesh.BoundaryRef(i).GetReference()) == 1)
	{
	  int num_face = i;
	  int num_elem = mesh.BoundaryRef(i).numElement(0);
	  int num_loc = mesh.Element(num_elem).GetPositionBoundary(num_face);
	  LocalPoints = var_problem.PointsQuadratureBoundary(num_elem, num_loc);
	  
	  // eval is allocated and set to 0
	  int nb_pts = LocalPoints.GetM();
	  for (int m = 0; m < eval.GetM(); m++)
	    {
	      eval(m)(i).Reallocate(nb_pts);
	      eval(m)(i).Zero();
	    }

	  FjInverseProblem<Dimension> inverseFj(mesh, num_elem);
	  
	  // we try to localize points on elements num_elem
	  typename Dimension::DimensionBoundary::VectR_N pts; IVect num;
	  var_problem.LocalizePointsBoundaryElement(Points, num_elem, num_loc,
						    inverseFj, pts, num);
	  
	  // the projector is computed between the two sets of points
	  Matrix<Real_wp> Proj;
	  var_problem.GetReferenceElement(num_elem).ComputeProjectionPointsSurf(num_loc, pts, LocalPoints, Proj);
          
	  // projecting val to quadrature points LocalPoints
	  Vector<Complexe> feval(pts.GetM()), contrib(LocalPoints.GetM());
	  contrib.Zero();
	  for (int m = 0; m < min(int(eval.GetM()), val.GetN()); m++)
	    {
	      for (int j = 0; j < pts.GetM(); j++)
		feval(j) = val(num(j), m);

	      Mlt(Proj, feval, contrib);
	      for (int j = 0; j < LocalPoints.GetM(); j++)
		eval(m)(i)(j) = contrib(j);
	    }
	}
  }
  

  /************************
   * VarSourceProblem_Fem *
   ************************/

  
  //! adding a volumetric source
  /*!
    \param[in,out] b_source right hand side to modify
    \param[in] f specification of the source 
    This method adds to b_source \f$ \int_\Omega f \varphi and \int g \nabla \varphi \f$
    where f is called by the method f.EvaluateVolumetricSource and g by f.EvaluateGradientSource
  */
  template<class Dimension> template<class T>
  void VarSourceProblem_Fem<Dimension>::
  AddVolumetricSource(const T& alpha, Vector<Vector<T> > & b_source,
		      Vector<VirtualSourceFEM<T, Dimension>* >& f) const
  {
    VectR_N s;

    // counting the number of source terms
    //IVect offset_source_num(f.GetM()); offset_source_num.Zero();
    int m0 = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      m0 = var_problem.nb_unknowns_hdg;

    int nb_source = 0;
    for (int k = 0; k < f.GetM(); k++)
      {
	//offset_source_num(k) = nb_source;
	if (f(k) != NULL)
	  nb_source += f(k)->GetNbUnknowns();
      }

    // no source => end of function
    if (nb_source <= 0)
      return;

    // number of components for each unknown (and gradients)
    int nb_mesh_num = var_problem.GetNbMeshNumberings();
    Vector<IVect> Nodle(nb_mesh_num);
    IVect type_u(nb_mesh_num);
    IVect nb_comp_mesh(nb_mesh_num), nb_grad_mesh(nb_mesh_num);
    Vector<bool> type_vec_f(nb_mesh_num), type_vec_df(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      {
        int nb_u = var_problem.GetNbComponentsUnknown(nm);
        int nb_du = var_problem.GetNbComponentsGradient(nm);
        nb_comp_mesh(nm) = nb_u; nb_grad_mesh(nm) = nb_du;
        if (nm == 0)
          type_u(nm) = var_problem.type_element;
        else
          type_u(nm) = var_problem.other_type_element(nm-1);

        type_vec_f(nm) = false;
        if (nb_u > 1)
          type_vec_f(nm) = true;

        type_vec_df(nm) = false;
        if (nb_du > 1)
          type_vec_df(nm) = true;
      }
    
    // number of components for each source
    IVect nb_comp_f(f.GetM()), nb_comp_df(f.GetM());
    IVect nb_f_rhs(nb_mesh_num), nb_df_rhs(nb_mesh_num);
    nb_f_rhs.Zero(); nb_df_rhs.Zero();
    for (int k = 0; k < f.GetM(); k++)
      {
        nb_comp_f(k) = 0; nb_comp_df(k) = 0;
	if (f(k) != NULL)
	  for (int m = m0; m < m0+f(k)->GetNbUnknowns(); m++)
	    {
              int nm = var_problem.mesh_num_unknown(m-m0);
              int nb_u = var_problem.GetNbComponentsUnknown(nm);
              int nb_du = var_problem.GetNbComponentsGradient(nm);
              nb_comp_f(k) += nb_u; nb_comp_df(k) += nb_du;
              nb_f_rhs(nm)++;
            }
      }

    Vector<IVect> offset_rhs(nb_mesh_num), offset_drhs(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      {
        offset_rhs(nm).Reallocate(nb_f_rhs(nm));
        offset_drhs(nm).Reallocate(nb_f_rhs(nm));
      }

    Vector<IVect> num_mesh_feval(f.GetM()), num_source_feval(f.GetM());
    Vector<Matrix<int> > pos_source_feval(f.GetM());
    Vector<IVect> num_mesh_geval(f.GetM()), num_source_geval(f.GetM());
    Vector<Matrix<int> > pos_source_geval(f.GetM());
    
    SetPoints<Dimension> PointsElem, PointsElemOrig;
    SetMatrices<Dimension> MatricesElem, MatricesElemOrig;
    int nb_modes = var_boundary.GetNbModesSource();
    Vector<Vector<Vector<Vector<T> > > > gradf_eval(nb_modes);
    Vector<Vector<Vector<Vector<T> > > > feval(nb_modes);
    Vector<Vector<T> > contrib(nb_source);
    Vector<bool> int_phi(nb_modes), int_grad(nb_modes);
    Vector<bool> source_phi(f.GetM()), source_grad(f.GetM());
    source_phi.Fill(false); source_grad.Fill(false);

    for (int i = 0; i < nb_modes; i++)
      {
        feval(i).Reallocate(nb_mesh_num);
        gradf_eval(i).Reallocate(nb_mesh_num);
      }
    
    // loop on elements of the mesh
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
	// we get the vertices of the face
	var_problem.mesh.GetVerticesElement(i, s);
	
	// finite element class associated with element i
	const ElementReference_Dim<Dimension>& Fb = var_problem.GetReferenceElement(i);
	
        int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();

        // we compute offsets to place each source in the right position of b_source
	bool source_present = false;
        int nb_dof_b = 0; nb_f_rhs.Zero(); nb_df_rhs.Zero();
	for (int ns = 0; ns < f.GetM(); ns++)
          {
            if (f(ns) != NULL)
              {
                f(ns)->InitElement(i, s);
                source_phi(ns) = f(ns)->IsNonNullVolumetricSource(s);
                source_grad(ns) = f(ns)->IsNonNullGradientSource(s);
                if (source_phi(ns) || source_grad(ns))
                  source_present = true;
                
                int nb_unknowns = f(ns)->GetNbUnknowns();
                if (source_phi(ns))                
                  {
                    num_mesh_feval(ns).Reallocate(nb_comp_f(ns));
                    num_source_feval(ns).Reallocate(nb_comp_f(ns));
                    pos_source_feval(ns).Reallocate(nb_comp_f(ns), nb_points_quadrature);
                    int offset_f_ = 0; 
                    for (int n = 0; n < nb_unknowns; n++)
                      {
                        int nm = var_problem.mesh_num_unknown(n);
                        int N = nb_comp_mesh(nm);
                        for (int k = 0; k < nb_comp_mesh(nm); k++)
                          {
                            num_mesh_feval(ns)(offset_f_ + k) = nm;
                            num_source_feval(ns)(offset_f_ + k) = nb_f_rhs(nm);
                            for (int j = 0; j < nb_points_quadrature; j++)
                              pos_source_feval(ns)(offset_f_+k, j) = N*j + k;
                          }
                        
                        offset_rhs(nm)(nb_f_rhs(nm)) = nb_dof_b + var_problem.GetOffsetDofUnknown(n);
                        offset_f_ += nb_comp_mesh(nm);
                        nb_f_rhs(nm)++;
                      }                      
                  }
                
                if (source_grad(ns))
                  {
                    int offset_df_ = 0;
                    num_mesh_geval(ns).Reallocate(nb_comp_df(ns));
                    num_source_geval(ns).Reallocate(nb_comp_df(ns));
                    pos_source_geval(ns).Reallocate(nb_comp_df(ns), nb_points_quadrature);
                    for (int n = 0; n < nb_unknowns; n++)
                      {
                        int nm = var_problem.mesh_num_unknown(n);
                        int N = nb_grad_mesh(nm);
                        for (int k = 0; k < N; k++)
                          {
                            num_mesh_geval(ns)(offset_df_ + k) = nm;
                            num_source_geval(ns)(offset_df_ + k) = nb_df_rhs(nm);
                            for (int j = 0; j < nb_points_quadrature; j++)
                              pos_source_geval(ns)(offset_df_+k, j) = N*j + k;
                          }
                        
                        offset_drhs(nm)(nb_df_rhs(nm)) = nb_dof_b + var_problem.GetOffsetDofUnknown(n);
                        offset_df_ += N;
                        nb_df_rhs(nm)++;
                      }                      
                  }
              }
            
            nb_dof_b += var_problem.GetNbDof();
          }
        
	// skipping to the next element if there is no source
	if (!source_present)
	  continue;
	
	// we retrieve number of dofs
	for (int nm = 0; nm < nb_mesh_num; nm++)
          Nodle(nm) = var_problem.GetDofNumberOnElement(i, nm);
    
	// we make the transformation Fi on quadrature points
	Fb.FjElemQuadrature(s, PointsElemOrig, var_problem.mesh, i);
	Fb.DFjElemQuadrature(s, PointsElemOrig, MatricesElemOrig, var_problem.mesh, i);

        // loop over modes
	int_phi.Fill(false); int_grad.Fill(false);
	for (int nsrc = 0; nsrc < nb_modes; nsrc++)
	  {
	    PointsElem = PointsElemOrig;
	    MatricesElem = MatricesElemOrig;
            // rotation/translation of points (because of modes computation)
	    if (f(0) != NULL)
	      f(0)->ModifyPoints(nsrc, s, PointsElem, MatricesElem);

            // loop over the sources contained in f
	    for (int ns = 0; ns < f.GetM(); ns++)
	      if (f(ns) != NULL)
		{
		  //int nb_unknowns = f(ns)->GetNbUnknowns();
		  Vector<T> f_eval(nb_comp_f(ns)), g_eval(nb_comp_df(ns));
		  f_eval.Zero(); g_eval.Zero();
		  if (source_phi(ns))
		    {                      
		      if (!int_phi(nsrc))
			{
                          // first time => we allocate the array feval
			  int_phi(nsrc) = true;
			  for (int nm = 0; nm < nb_mesh_num; nm++)
                            {
                              feval(nsrc)(nm).Reallocate(nb_f_rhs(nm));
                              for (int k = 0; k < nb_f_rhs(nm); k++)
                                {
                                  feval(nsrc)(nm)(k).Reallocate(nb_comp_mesh(nm)*nb_points_quadrature);
                                  feval(nsrc)(nm)(k).Zero();
                                }
                            }
                        }
                      
		      // \int_{K_i} f \phi_i = \int_\hat{K} J_i f(xi_i) \phi_i
                      // evaluation of f(xi_i) => they are stored in feval
		      for (int j = 0; j < nb_points_quadrature; j++)
			{
			  f(ns)->EvaluateVolumetricSource(i, j, PointsElem.
							  GetPointQuadrature(j), f_eval);

			  for (int k = 0; k < f_eval.GetM(); k++)
			    {
                              int nm = num_mesh_feval(ns)(k);
                              int ks = num_source_feval(ns)(k);
                              int pos = pos_source_feval(ns)(k, j);
                              feval(nsrc)(nm)(ks)(pos) = f_eval(k); 
			    }
			}
		    }
                  
		  if (source_grad(ns))
		    {
		      if (!int_grad((nsrc)))
			{
			  int_grad(nsrc) = true;
			  for (int nm = 0; nm < nb_mesh_num; nm++)
                            {
                              gradf_eval(nsrc)(nm).Reallocate(nb_df_rhs(nm));
                              for (int k = 0; k < nb_df_rhs(nm); k++)
                                {
                                  gradf_eval(nsrc)(nm)(k).Reallocate(nb_grad_mesh(nm)*nb_points_quadrature);
                                  gradf_eval(nsrc)(nm)(k).Zero();
                                }
                            }
			}
                      
                      // evaluation of g(xi_i) => they are stored in gradf_eval
                      for (int j = 0; j < nb_points_quadrature; j++)
                        {
                          f(ns)->EvaluateGradientSource(i, j, PointsElem.
                                                        GetPointQuadrature(j), g_eval);
                          
                          for (int k = 0; k < g_eval.GetM(); k++)
                            {
                              int nm = num_mesh_geval(ns)(k);
                              int ks = num_source_geval(ns)(k);
                              int pos = pos_source_geval(ns)(k, j);
                              gradf_eval(nsrc)(nm)(ks)(pos) = g_eval(k); 
                            }
                        }
		    }
		}
	  }
	
	if (f(0) != NULL)
	  f(0)->ModifyEvaluationVolume(int_phi, int_grad, feval, gradf_eval,
                                       type_vec_f, type_vec_df);
        
        // loop over modes
	for (int nsrc = 0; nsrc < b_source.GetM(); nsrc++)
	  {
	    if (int_phi(nsrc))
	      {
                for (int nm = 0; nm < nb_mesh_num; nm++)
                  {
                    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
                    
                    // we compute the integrals int( f phi_j dx ) j=0..nb_dof_elt-1 
                    switch (type_u(nm))
                      {
                      case 1:
                        var_problem.GetReferenceElementH1(i, nm).
                          ComputeIntegral(MatricesElemOrig, feval(nsrc)(nm), contrib,
                                       mesh_num, i);
                        break;
                      case 2:
                        var_problem.GetReferenceElementHcurl(i, nm).
                          ComputeIntegral(MatricesElemOrig, feval(nsrc)(nm), contrib,
                                       mesh_num, i);
                        break;
                      case 3:
                        var_problem.GetReferenceElementHdiv(i, nm).
                          ComputeIntegral(MatricesElemOrig, feval(nsrc)(nm), contrib,
                                       mesh_num, i);
                        break;
                      }
                    
                    // we add contributions to b_source
                    for (int m = 0; m < nb_f_rhs(nm); m++)
                      for (int j = 0; j < Nodle(nm).GetM(); j++)
                        {
                          int num_dof_loc = Nodle(nm)(j);
                          int num_dof = num_dof_loc + offset_rhs(nm)(m);
                          if (num_dof_loc >= 0)
                            b_source(nsrc)(num_dof) += alpha*contrib(m)(j);
                        }
                  }
              }
	    
	    if (int_grad(nsrc))
	      {
                for (int nm = 0; nm < nb_mesh_num; nm++)
                  {
                    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
                    
                    // we compute the integrals int( f \nabla phi_j dx ) j=0..nb_dof_elt-1 
                    switch (type_u(nm))
                      {
                      case 1:
                        var_problem.GetReferenceElementH1(i, nm).
                          ComputeIntegralGradient(MatricesElemOrig, gradf_eval(nsrc)(nm), contrib,
                                                  mesh_num, i);
                        break;
                      case 2:
                        var_problem.GetReferenceElementHcurl(i, nm).
                          ComputeIntegralGradient(MatricesElemOrig, gradf_eval(nsrc)(nm), contrib,
                                                  mesh_num, i);
                        break;
                      case 3:
                        var_problem.GetReferenceElementHdiv(i, nm).
                          ComputeIntegralGradient(MatricesElemOrig, gradf_eval(nsrc)(nm), contrib,
                                                  mesh_num, i);
                        break;
                      }
                    
                    for (int m = 0; m < nb_df_rhs(nm); m++)
                      for (int j = 0; j < Nodle(nm).GetM(); j++)
                        {
                          int num_dof_loc = Nodle(nm)(j);
                          int num_dof = num_dof_loc + offset_drhs(nm)(m);
                          if (num_dof_loc >= 0)
                            b_source(nsrc)(num_dof) += alpha*contrib(m)(j);
                        }
                  }
              }
	  }	
      }
  }
  
    
  //! we set dirichlet nodes at values fixed by f.EvaluateFunction
  //! (inhomogeneous dirichlet condition)
  /*!
    \param[in,out] b_source right hand side to modify
    \param[in] f specification of source
   */
  template<class Dimension> template<class T>
  void VarSourceProblem_Fem<Dimension>::
  SetDirichletSource(Vector<Vector<T> >& b_source, 
		     Vector<VirtualSourceFEM<T, Dimension>* >& f) const
  {
    // no Dirichlet for some variational formulations
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      return;

    if ((var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
	&& (var_problem.FirstOrderFormulationDG()))
      return;
    
    // we set dirichlet projection of f on surfacic degrees of freedom
    T one; SetComplexOne(one);
    SetSurfacicProjection(one, b_source, f);
  }
  
  
  //! Add the projection of f on all degrees of freedom
  /*!
    \param[in] alpha we add alpha*f
    \param[in,out] b_src vector whose we add components of f on basis functions
    \param[in] f function to project
   */
  template<class Dimension> template<class T>
  void VarSourceProblem_Fem<Dimension>
  ::AddVolumetricProjection(const T& alpha, Vector<Vector<T> >& b_src,
			    Vector<VirtualProjectorFEM<T, Dimension>* >& f) const
  {
    VectR_N s;
    int m0 = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      m0 = var_problem.nb_unknowns_hdg;

    int nb_source = 0;
    for (int k = 0; k < f.GetM(); k++)
      {
	//offset_source_num(k) = nb_source;
	if (f(k) != NULL)
	  nb_source += f(k)->GetNbUnknowns();
      }

    if (nb_source <= 0)
      return;
    
    // number of components for each unknown (and gradients)
    int nb_mesh_num = var_problem.GetNbMeshNumberings();
    Vector<IVect> Nodle(nb_mesh_num);
    IVect type_u(nb_mesh_num);
    IVect nb_comp_mesh(nb_mesh_num);
    Vector<bool> type_vec_f(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      {
        int nb_u = var_problem.GetNbComponentsUnknown(nm);
        nb_comp_mesh(nm) = nb_u;
        if (nm == 0)
          type_u(nm) = var_problem.type_element;
        else
          type_u(nm) = var_problem.other_type_element(nm-1);

        type_vec_f(nm) = false;
        if (nb_u > 1)
          type_vec_f(nm) = true;
      }

    // number of components for each source
    IVect nb_comp_f(f.GetM());
    IVect nb_f_rhs(nb_mesh_num);
    nb_comp_f.Zero(); nb_f_rhs.Zero();
    for (int k = 0; k < f.GetM(); k++)
      {
        nb_comp_f(k) = 0;
	if (f(k) != NULL)
	  for (int m = m0; m < m0+f(k)->GetNbUnknowns(); m++)
	    {
              int nm = var_problem.mesh_num_unknown(m-m0);
              int nb_u = var_problem.GetNbComponentsUnknown(nm);
              nb_comp_f(k) += nb_u;
              nb_f_rhs(nm)++;
            }
      }
    
    Vector<IVect> offset_rhs(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      offset_rhs(nm).Reallocate(nb_f_rhs(nm));
    
    Vector<IVect> num_mesh_feval(f.GetM()), num_source_feval(f.GetM());
    Vector<Matrix<int> > pos_source_feval(f.GetM());

    Vector<SetPoints<Dimension> > PointsElem(nb_mesh_num), PointsElemOrig(nb_mesh_num);
    Vector<SetMatrices<Dimension> > MatricesElem(nb_mesh_num), MatricesElemOrig(nb_mesh_num);
    
    int nb_modes = var_boundary.GetNbModesSource();
    Vector<Vector<Vector<Vector<T> > > > feval(nb_modes);
    Vector<Vector<T> > contrib(nb_source);

    Vector<VectBool> DofUsed(nb_mesh_num);
    for (int i = 0; i < nb_mesh_num; i++)
      {
        DofUsed(i).Reallocate(var_problem.GetMeshNumbering(i).GetNbDof());
        DofUsed(i).Fill(false);
      }
    
    for (int i = 0; i < nb_modes; i++)
      feval(i).Reallocate(nb_mesh_num);

    IVect nb_dof_elt(nb_mesh_num), nb_dof_loc(nb_mesh_num);
    int nb_dof_max = 0;
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
	// we get the vertices of the face
	var_problem.mesh.GetVerticesElement(i, s);
	for (int nm = 0; nm < nb_mesh_num; nm++)
          {
            const ElementReference_Dim<Dimension>& Fb = var_problem.GetReferenceElement(i, nm);
        
            // number of points where the function f is to be evaluated
            nb_dof_elt(nm) = Fb.GetNbPointsDofInside();
            
            // number of dofs in the face
            nb_dof_loc(nm) = Fb.GetNbDof();
            nb_dof_max = max(nb_dof_max, nb_dof_elt(nm));
          }
	
	// initialisation of the source for element i
        int nb_dof_b = 0; nb_f_rhs.Zero();
	for (int ns = 0; ns < f.GetM(); ns++)
	  {
            if (f(ns) != NULL)
              {
                f(ns)->InitElement(i, s);
                
                num_mesh_feval(ns).Reallocate(nb_comp_f(ns));
                num_source_feval(ns).Reallocate(nb_comp_f(ns));
                pos_source_feval(ns).Reallocate(nb_comp_f(ns), nb_dof_max);
                int offset_f_ = 0;
                int nb_unknowns = f(ns)->GetNbUnknowns();
                for (int n = 0; n < nb_unknowns; n++)
                  {
                    int nm = var_problem.mesh_num_unknown(n);
                    int N = nb_comp_mesh(nm);
                    for (int k = 0; k < nb_comp_mesh(nm); k++)
                      {
                        num_mesh_feval(ns)(offset_f_ + k) = nm;
                        num_source_feval(ns)(offset_f_ + k) = nb_f_rhs(nm);
                        for (int j = 0; j < nb_dof_elt(nm); j++)
                          pos_source_feval(ns)(offset_f_+k, j) = N*j + k;
                      }
                    
                    offset_rhs(nm)(nb_f_rhs(nm)) = nb_dof_b + var_problem.GetOffsetDofUnknown(n);
                    offset_f_ += nb_comp_mesh(nm);
                    nb_f_rhs(nm)++;
                  }                      
              }
            
            nb_dof_b += var_problem.GetNbDof();
          }
	
	// dof numbers
	for (int nm = 0; nm < nb_mesh_num; nm++)
          {
            Nodle(nm) = var_problem.GetDofNumberOnElement(i, nm);
            
            // we make the transformation Fi on dof points
            var_problem.GetReferenceElement(i, nm).FjElemDof(s, PointsElemOrig(nm), var_problem.mesh, i);
            
            // and DF_i
            var_problem.GetReferenceElement(i, nm).DFjElemDof(s, PointsElemOrig(nm), MatricesElemOrig(nm), var_problem.mesh, i);
          }
	
	for (int nsrc = 0; nsrc < nb_modes; nsrc++)
	  {
	    PointsElem = PointsElemOrig;
	    MatricesElem = MatricesElemOrig;
	    if (f(0) != NULL)
	      for (int nm = 0; nm < nb_mesh_num; nm++)
                f(0)->ModifyPoints(nsrc, s, PointsElem(nm), MatricesElem(nm));

            for (int nm = 0; nm < nb_mesh_num; nm++)
	      {
                feval(nsrc)(nm).Reallocate(nb_f_rhs(nm));
                for (int k = 0; k < nb_f_rhs(nm); k++)
                  {
                    feval(nsrc)(nm)(k).Reallocate(nb_comp_mesh(nm)*nb_dof_elt(nm));
                    feval(nsrc)(nm)(k).Zero();
                  }
              }

	    for (int ns = 0; ns < f.GetM(); ns++)
	      if (f(ns) != NULL)
		{
		  Vector<T> f_eval(nb_comp_f(ns));
		  f_eval.Zero();
		  
		  // evaluation of f on dof points
		  for (int nm = 0; nm < nb_mesh_num; nm++)
                    for (int j = 0; j < nb_dof_elt(nm); j++)
                      {
                        f(ns)->EvaluateFunction(i, j, PointsElem(nm).GetPointDof(j), f_eval);
                        
                        for (int k = 0; k < f_eval.GetM(); k++)
                          {
                            int nm0 = num_mesh_feval(ns)(k);
                            if (nm0 == nm)
                              {
                                int ks = num_source_feval(ns)(k);
                                int pos = pos_source_feval(ns)(k, j);
                                feval(nsrc)(nm)(ks)(pos) = f_eval(k); 
                              }
                          }
                      }	    
		}
	  }
	
	if (f(0) != NULL)
	  f(0)->ModifyEvaluationProjection(feval, type_vec_f);
	
	PointsElem = PointsElemOrig;
	MatricesElem = MatricesElemOrig;
	
	for (int nsrc = 0; nsrc < b_src.GetM(); nsrc++)
	  {
            for (int nm = 0; nm < nb_mesh_num; nm++)
              {
                const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
                
                // evaluation of the projection in the element i
                // for interpolatory element contrib = feval
                // for hierarchic basis, we need to solve a linear system
                switch(type_u(nm))
                  {
                  case 1:
                    var_problem.GetReferenceElementH1(i, nm)
                      .ComputeProjectionDof(MatricesElem(nm), feval(nsrc)(nm), contrib,
                                            mesh_num, i);
                    break;
                 case 2:
                    var_problem.GetReferenceElementHcurl(i, nm)
                      .ComputeProjectionDof(MatricesElem(nm), feval(nsrc)(nm), contrib,
                                            mesh_num, i);
                    break;
                 case 3:
                    var_problem.GetReferenceElementHdiv(i, nm)
                      .ComputeProjectionDof(MatricesElem(nm), feval(nsrc)(nm), contrib,
                                            mesh_num, i);
                    break;
                  }
                
                for (int j = 0; j < nb_dof_loc(nm); j++)
                  {
                    int num_dof_loc = Nodle(nm)(j);
                    if (num_dof_loc >= 0)
                      if (!DofUsed(nm)(num_dof_loc))
                        {
                          for (int m = 0; m < nb_f_rhs(nm); m++)
                            {
                              int num_dof = num_dof_loc + offset_rhs(nm)(m);
                              b_src(nsrc)(num_dof) += alpha*contrib(m)(j);
                            }
                          
                          if (nsrc == b_src.GetM()-1)
                            DofUsed(nm)(num_dof_loc) = true;
                        }
                  }
              }	
          }
      }
  }
  
  
  //! Sets the projection of f on surfacic degrees of freedom
  /*!
    \param[in,out] b_src vector whose we set components of f on basis functions
    \param[in] f function to project
   */
  template<class Dimension> template<class T>
  void VarSourceProblem_Fem<Dimension>
  ::SetSurfacicProjection(const T& alpha, Vector<Vector<T> >& b_src,
			  Vector<VirtualSourceFEM<T, Dimension>* >& f) const
  {
    // cancelling values on dirichlet dofs
    for (int nsrc = 0; nsrc < b_src.GetM(); nsrc++)
      for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
	b_src(nsrc)(var_boundary.GetDirichletDofNumber(i)) = 0;
    
    // then we can add terms due to f
    VectR_N s;
    
    int nb_source = 0;
    bool presence_dirichlet = false;
    for (int k = 0; k < f.GetM(); k++)
      {
	if (f(k) != NULL)
	  {
	    if (f(k)->PresenceDirichlet())
	      presence_dirichlet = true;
	    
	    nb_source += f(k)->GetNbUnknowns();
	  }
      }

    if (nb_source <= 0)
      return;

    if (!presence_dirichlet)
      return;

    // number of components for each unknown (and gradients)
    int nb_mesh_num = var_problem.GetNbMeshNumberings();
    Vector<IVect> Nodle(nb_mesh_num);
    IVect type_u(nb_mesh_num);
    IVect nb_comp_mesh(nb_mesh_num);
    Vector<bool> type_vec_f(nb_mesh_num), type_vec_df(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      {
        int nb_u = var_problem.GetNbComponentsUnknown(nm);
        nb_comp_mesh(nm) = nb_u;
        if (nm == 0)
          type_u(nm) = var_problem.type_element;
        else
          type_u(nm) = var_problem.other_type_element(nm-1);

        type_vec_f(nm) = false;
        if (nb_u > 1)
          type_vec_f(nm) = true;
      }

    int m0 = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      m0 = var_problem.nb_unknowns_hdg;

    IVect nb_comp_f(f.GetM());
    IVect nb_f_rhs(nb_mesh_num);    
    nb_comp_f.Zero(); nb_f_rhs.Zero();
    for (int k = 0; k < f.GetM(); k++)
      {
        nb_comp_f(k) = 0;
	if (f(k) != NULL)
	  for (int m = m0; m < m0+f(k)->GetNbUnknowns(); m++)
	    {
              int nm = var_problem.mesh_num_unknown(m-m0);
              int nb_u = var_problem.GetNbComponentsUnknown(nm);
              nb_comp_f(k) += nb_u;
              nb_f_rhs(nm)++;
            }
      }

    Vector<IVect> offset_rhs(nb_mesh_num), num_unknown_rhs(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      {
        offset_rhs(nm).Reallocate(nb_f_rhs(nm));
        num_unknown_rhs(nm).Reallocate(nb_f_rhs(nm));
      }
    
    Vector<IVect> num_mesh_feval(f.GetM()), num_source_feval(f.GetM());
    Vector<Matrix<int> > pos_source_feval(f.GetM());

    Vector<SetPoints<Dimension> > PointsElem(nb_mesh_num), PointsElemOrig(nb_mesh_num);
    Vector<SetMatrices<Dimension> > MatricesElem(nb_mesh_num), MatricesElemOrig(nb_mesh_num);
    
    int nb_modes = var_boundary.GetNbModesSource();
    Vector<Vector<Vector<Vector<T> > > > feval(nb_modes);
    Vector<Vector<T> > contrib(nb_source);
    
    // results of different processors and edges will be added
    // in NbDofAdded, we will just count the number of times each dof is added
    Vector<IVect> NbDofAdded(var_problem.nb_unknowns);
    for (int n = 0; n < var_problem.nb_unknowns; n++)
      {
        int i = var_problem.mesh_num_unknown(n);
        // to take into account Dirichlet dofs for PMLs dofs
        int Nmax = var_problem.GetMeshNumbering(i).GetNbDof();
        Nmax = max(Nmax, var_problem.GetNbDofUnknown(n));
        
        NbDofAdded(n).Reallocate(Nmax);
        NbDofAdded(n).Zero();
      }
    
    for (int i = 0; i < nb_modes; i++)
      feval(i).Reallocate(nb_mesh_num);
    
    // loop on all referenced boundaries
    for (int i = 0; i < var_problem.mesh.GetNbBoundaryRef(); i++)
      {
	// reference of the boundary
	int ref = var_problem.mesh.BoundaryRef(i).GetReference();

	bool source_present = false;
	for (int ns = 0; ns < f.GetM(); ns++)
	  if (f(ns) != NULL)
	    if (f(ns)->IsNonNullDirichletSource(ref))
	      source_present = true;
	
	// if there is no source, we skip to the next element
	if (!source_present)
	  continue;

	// we get global face number, element number, local face number ...
	int num_face = i;
	int num_elem = var_problem.mesh.Boundary(num_face).numElement(0);
	if (num_elem < 0)
	  continue;

        IVect nb_dof_surf(nb_mesh_num);
        int num_loc = var_problem.mesh.Element(num_elem).GetPositionBoundary(num_face);
        // we get the vertices of the element
	var_problem.mesh.GetVerticesElement(num_elem, s);

        int nb_dof_max = 0;
        for (int nm = 0; nm < nb_mesh_num; nm++)
          {
            const ElementReference_Dim<Dimension>& Fb = var_problem.GetReferenceElement(num_elem, nm);
            nb_dof_surf(nm) = Fb.GetNbPointsDofSurface(num_loc);
            nb_dof_max = max(nb_dof_max, nb_dof_surf(nm));
          }
        
	// initialization of some variables
        int nb_dof_b = 0; nb_f_rhs.Zero();
	for (int ns = 0; ns < f.GetM(); ns++)
          {
            if (f(ns) != NULL)
              if (f(ns)->IsNonNullDirichletSource(ref))
                {
                  f(ns)->InitSurface(i, num_face, num_elem, num_loc);
                  num_mesh_feval(ns).Reallocate(nb_comp_f(ns));
                  num_source_feval(ns).Reallocate(nb_comp_f(ns));
                  pos_source_feval(ns).Reallocate(nb_comp_f(ns), nb_dof_max);
                  int offset_f_ = 0;
                  int nb_unknowns = f(ns)->GetNbUnknowns();
                  for (int n = 0; n < nb_unknowns; n++)
                    {
                      int nm = var_problem.mesh_num_unknown(n);
                      int N = nb_comp_mesh(nm);
                      for (int k = 0; k < nb_comp_mesh(nm); k++)
                        {
                          num_mesh_feval(ns)(offset_f_ + k) = nm;
                          num_source_feval(ns)(offset_f_ + k) = nb_f_rhs(nm);
                          for (int j = 0; j < nb_dof_surf(nm); j++)
                            pos_source_feval(ns)(offset_f_+k, j) = N*j + k;
                        }
                      
                      num_unknown_rhs(nm)(nb_f_rhs(nm)) = n;
                      offset_rhs(nm)(nb_f_rhs(nm)) = nb_dof_b + var_problem.GetOffsetDofUnknown(n);
                      offset_f_ += nb_comp_mesh(nm);
                      nb_f_rhs(nm)++;
                    }
                }
            
            nb_dof_b += var_problem.GetNbDof();
          }
        
	// and dof numbers
	for (int nm = 0; nm < nb_mesh_num; nm++)
          {
            Nodle(nm) = var_problem.GetDofNumberOnElement(num_elem, nm);
            const ElementReference_Dim<Dimension>& Fb = var_problem.GetReferenceElement(num_elem, nm);
            
            // we apply the transformation Fi on dof points
            Fb.FjElemDof(s, PointsElemOrig(nm), var_problem.mesh, num_elem);
            Fb.DFjElemDof(s, PointsElemOrig(nm), MatricesElemOrig(nm), var_problem.mesh, num_elem);
            
            Fb.FjSurfaceElemDof(s, PointsElemOrig(nm), var_problem.mesh, num_elem, num_loc);		
            Fb.DFjSurfaceElemDof(s, PointsElemOrig(nm), MatricesElemOrig(nm), var_problem.mesh, num_elem, num_loc);
          }
        
	for (int nsrc = 0; nsrc < nb_modes; nsrc++)
	  {	    
	    PointsElem = PointsElemOrig;
	    MatricesElem = MatricesElemOrig;
	    
	    if (f(0) != NULL)
	      for (int nm = 0; nm < nb_mesh_num; nm++)
                f(0)->ModifyPoints(nsrc, s, PointsElem(nm), MatricesElem(nm));

            for (int nm = 0; nm < nb_mesh_num; nm++)
	      {
                feval(nsrc)(nm).Reallocate(nb_f_rhs(nm));
                for (int k = 0; k < nb_f_rhs(nm); k++)
                  {
                    feval(nsrc)(nm)(k).Reallocate(nb_comp_mesh(nm)*nb_dof_surf(nm));
                    feval(nsrc)(nm)(k).Zero();
                  }
              }

	    for (int ns = 0; ns < f.GetM(); ns++)
	      if (f(ns) != NULL)
		if (f(ns)->IsNonNullDirichletSource(ref))
		  {
		    //int nb_unknowns = f(ns)->GetNbUnknowns();
		    Vector<T> f_eval(nb_comp_f(ns));
		    f_eval.Zero();
                    
		    // evaluation of f on dof points
                    for (int nm = 0; nm < nb_mesh_num; nm++)
                      for (int j = 0; j < nb_dof_surf(nm); j++)
                        {
                          f(ns)->EvaluateFunction(num_elem, j, PointsElem(nm).GetPointDofBoundary(j), f_eval);
                          
                          for (int k = 0; k < f_eval.GetM(); k++)
                            {
                              int nm0 = num_mesh_feval(ns)(k);
                              if (nm0 == nm)
                                {
                                  int ks = num_source_feval(ns)(k);
                                  int pos = pos_source_feval(ns)(k, j);
                                  feval(nsrc)(nm)(ks)(pos) = f_eval(k);
                                }
                          }
                      }
                  }
	  }
        
	if (f(0) != NULL)
	  f(0)->ModifyEvaluationProjection(feval, type_vec_f);
	
	// evaluation of the projection in the element i
	// for interpolatory element contrib = feval
	// for hierarchic basis, we need to solve a linear system
	for (int nsrc = 0; nsrc < b_src.GetM(); nsrc++)
	  {
            for (int nm = 0; nm < nb_mesh_num; nm++)
              {
                const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
                
                switch (type_u(nm))
                  {
                  case 1:
                    var_problem.GetReferenceElementH1(num_elem, nm).
                      ComputeProjectionSurfaceDof(MatricesElemOrig(nm), feval(nsrc)(nm), contrib,
                                                  mesh_num, num_elem, num_loc);
                    break;
                  case 2:
                    var_problem.GetReferenceElementHcurl(num_elem, nm).
                      ComputeProjectionSurfaceDof(MatricesElemOrig(nm), feval(nsrc)(nm), contrib,
                                                  mesh_num, num_elem, num_loc);
                    break;
                  case 3:
                    var_problem.GetReferenceElementHdiv(num_elem, nm).
                      ComputeProjectionSurfaceDof(MatricesElemOrig(nm), feval(nsrc)(nm), contrib,
                                                  mesh_num, num_elem, num_loc);
                    break;
                  }

                const ElementReference_Dim<Dimension>& Fb = var_problem.GetReferenceElement(num_elem, nm);
                // addition only for surfacic dofs
                for (int k = 0; k < Fb.GetNbDofBoundary(num_loc); k++)
                  {
                    int j = Fb.GetLocalNumber(num_loc, k);
                    int num_dof_loc = Nodle(nm)(j);
                    
                    if (num_dof_loc >= 0)
                      {
                    
                        for (int m = 0; m < nb_f_rhs(nm); m++)
                          {
                            int num_dof = num_dof_loc + offset_rhs(nm)(m);
                            if (var_boundary.IsDofDirichlet(num_dof))
                              {
                                if (nsrc == 0)
                                  NbDofAdded(num_unknown_rhs(nm)(m))(num_dof_loc)++;

                                b_src(nsrc)(num_dof) += alpha*contrib(m)(k);
                              }
                          }
                      }
                  }
	      }
	  }
      }
   
    // finding for each dof the number of times where its contribution has been added
    for (int i = 0; i < var_problem.nb_unknowns; i++)
      var_problem.AssembleDirichlet(NbDofAdded(i), i, true);

    // then we set 1/NbDofAdded  b_src
    // such that b_src will contain the correct value on the final assembling
    for (int n = 0; n < var_problem.nb_unknowns; n++)
      {
        int nm = var_problem.mesh_num_unknown(n);
        const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
        for (int j = 0; j < mesh_num.GetNbDof(); j++)    
          if (NbDofAdded(n)(j) > 0)
            {
              Real_wp coef = Real_wp(1) / NbDofAdded(n)(j);
              for (int m = 0; m < f.GetM(); m++)
                {			
                  int num_dof = j + var_problem.offset_dof_unknown(n) + m*var_problem.GetNbDof();
                  if (var_boundary.IsDofDirichlet(num_dof))
                    for (int i = 0; i < b_src.GetM(); i++)
                      b_src(i)(num_dof) *= coef;
                }
            }
      }
  }
  
  
  //! adding a surfacic source
  /*!
    \param[in,out] b_source right hand side to modify
    \param[in] f specification of the source 
    This method adds to b_source \f$ \int_\Sigma g \varphi
    + \int_\Sigma h \cdot \nabla_\Sigma \varphi \f$
    where g is called by method f.EvaluateSurfacicSource
    and h by f.EvaluateSurfacicSourceGradient
  */
  template<class Dimension> template<class T>
  void VarSourceProblem_Fem<Dimension>
  ::AddSurfacicSource(const T& alpha, Vector<Vector<T> >& b_source,
		      Vector<VirtualSourceFEM<T, Dimension>* >& f) const
  {
    VectR_N s;

    int nb_source = 0;
    for (int k = 0; k < f.GetM(); k++)
      {
	if (f(k) != NULL)
	  nb_source += f(k)->GetNbUnknowns();
      }

    if (nb_source <= 0)
      return;

    // number of components for each unknown (and gradients)
    int nb_mesh_num = var_problem.GetNbMeshNumberings();
    IVect type_u(nb_mesh_num);
    IVect nb_comp_mesh(nb_mesh_num), nb_grad_mesh(nb_mesh_num);
    Vector<bool> type_vec_f(nb_mesh_num), type_vec_df(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      {
        int nb_u = var_problem.GetNbComponentsUnknown(nm);
        int nb_du = var_problem.GetNbComponentsGradient(nm);
        nb_comp_mesh(nm) = nb_u; nb_grad_mesh(nm) = nb_du;
        if (nm == 0)
          type_u(nm) = var_problem.type_element;
        else
          type_u(nm) = var_problem.other_type_element(nm-1);

        type_vec_f(nm) = false;
        if (nb_u > 1)
          type_vec_f(nm) = true;

        type_vec_df(nm) = false;
        if (nb_du > 1)
          type_vec_df(nm) = true;
      }

    int m0 = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      m0 = var_problem.nb_unknowns_hdg;

    IVect nb_comp_f(f.GetM()), nb_comp_df(f.GetM());
    IVect nb_f_rhs(nb_mesh_num), nb_df_rhs(nb_mesh_num), nb_f_hdg(nb_mesh_num);
    nb_f_rhs.Zero(); nb_df_rhs.Zero(); nb_f_hdg.Zero();
    for (int k = 0; k < f.GetM(); k++)
      {
        nb_comp_f(k) = 0; nb_comp_df(k) = 0;
	if (f(k) != NULL)
	  for (int m = 0; m < f(k)->GetNbUnknowns(); m++)
	    {
              int nm = var_problem.mesh_num_unknown(m);
              int nb_u = var_problem.GetNbComponentsUnknown(nm);
              int nb_du = var_problem.GetNbComponentsGradient(nm);
              nb_comp_f(k) += nb_u; nb_comp_df(k) += nb_du;
              if ((var_problem.FormulationDG() == ElementReference_Base::HDG) &&(m < m0))
                nb_f_hdg(nm)++;
              else
                nb_f_rhs(nm)++;
            }
      }

    Vector<IVect> offset_rhs(nb_mesh_num), offset_drhs(nb_mesh_num), offset_hdg(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      {
        offset_rhs(nm).Reallocate(nb_f_rhs(nm));
        offset_drhs(nm).Reallocate(nb_f_rhs(nm));
        offset_hdg(nm).Reallocate(nb_f_hdg(nm));
      }
    
    Vector<IVect> num_mesh_feval(f.GetM()), num_source_feval(f.GetM());
    Vector<Matrix<int> > pos_source_feval(f.GetM());
    Vector<IVect> num_mesh_geval(f.GetM()), num_source_geval(f.GetM());
    Vector<Matrix<int> > pos_source_geval(f.GetM());
    Vector<Vector<bool> > hdg_source_feval(f.GetM());
    
    int nb_modes = var_boundary.GetNbModesSource();
    Vector<Vector<Vector<Vector<T> > > > feval(nb_modes), feval_hdg(nb_modes);
    Vector<Vector<T> > contrib;
    Vector<Vector<Vector<Vector<T> > > > feval_diff(nb_modes);
    
    SetPoints<Dimension> PointsElemOrig, PointsElem;
    SetMatrices<Dimension> MatricesElemOrig, MatricesElem;

    Vector<Vector<T> > contrib_surf; Vector<T> contrib_dof;
    Vector<IVect> Nodle(nb_mesh_num);

    for (int i = 0; i < nb_modes; i++)
      {
        feval(i).Reallocate(nb_mesh_num);
        feval_hdg(i).Reallocate(nb_mesh_num);
        feval_diff(i).Reallocate(nb_mesh_num);
      }
    
    // inhomogeneous Neumann-like condition
    for (int i = 0; i < var_problem.mesh.GetNbBoundaryRef(); i++)
      {
	// reference of the boundary
	int ref = var_problem.mesh.BoundaryRef(i).GetReference();
	
	bool source_phi = false, source_grad = false;
	for (int k = 0; k < f.GetM(); k++)
	  if (f(k) != NULL)
	    {
	      if (f(k)->IsNonNullSurfacicSource(ref))
		source_phi = true;

	      if (f(k)->IsNonNullSurfacicSourceGradient(ref))
		source_grad = true;
	    }

	// this face is skipped if no source is present
	if (!source_phi && !source_grad)
	  continue;

	int num_face = i;
	int num_elem = var_problem.mesh.Boundary(num_face).numElement(0);
	if (num_elem < 0)
	  continue;

	int num_loc = var_problem.mesh.Element(num_elem).GetPositionBoundary(num_face);
		
	// we get the vertices of the element
	var_problem.mesh.GetVerticesElement(num_elem, s);
	const ElementReference_Dim<Dimension>& Fb = var_problem.GetReferenceElement(num_elem);
        
	int nb_points_quadrature_boundary = Fb.GetNbQuadBoundary(num_loc);
	
	// initialization of some variables
        int nb_dof_b = 0; nb_f_rhs.Zero(); nb_df_rhs.Zero(); nb_f_hdg.Zero();
        bool hdg_source = false;
	for (int ns = 0; ns < f.GetM(); ns++)
	  {
            if (f(ns) != NULL)
              {
                f(ns)->InitSurface(i, num_face, num_elem, num_loc);

                int nb_unknowns = f(ns)->GetNbUnknowns();
                if (f(ns)->IsNonNullSurfacicSource(ref))
                  {
                    hdg_source_feval(ns).Reallocate(nb_comp_f(ns));
                    num_mesh_feval(ns).Reallocate(nb_comp_f(ns));
                    num_source_feval(ns).Reallocate(nb_comp_f(ns));
                    pos_source_feval(ns).Reallocate(nb_comp_f(ns), nb_points_quadrature_boundary);
                    hdg_source_feval(ns).Fill(false);
                    int offset_f_ = 0; 
                    for (int n = 0; n < nb_unknowns; n++)
                      {
                        int nm = var_problem.mesh_num_unknown(n);
                        int N = nb_comp_mesh(nm);
                        if ((var_problem.FormulationDG() == ElementReference_Base::HDG) &&(n < m0))
                          {                            
                            for (int k = 0; k < nb_comp_mesh(nm); k++)
                              {
                                hdg_source_feval(ns)(offset_f_ + k) = true; hdg_source = true;
                                num_mesh_feval(ns)(offset_f_ + k) = nm;
                                num_source_feval(ns)(offset_f_ + k) = nb_f_hdg(nm);
                                for (int j = 0; j < nb_points_quadrature_boundary; j++)
                                  pos_source_feval(ns)(offset_f_+k, j) = N*j + k;
                              }
                            
                            offset_hdg(nm)(nb_f_hdg(nm)) = nb_dof_b + var_problem.GetOffsetDofUnknown(n);
                            nb_f_hdg(nm)++;
                          }
                        else
                          {
                            for (int k = 0; k < nb_comp_mesh(nm); k++)
                              {
                                num_mesh_feval(ns)(offset_f_ + k) = nm;
                                num_source_feval(ns)(offset_f_ + k) = nb_f_rhs(nm);
                                for (int j = 0; j < nb_points_quadrature_boundary; j++)
                                  pos_source_feval(ns)(offset_f_+k, j) = N*j + k;
                              }
                            
                            offset_rhs(nm)(nb_f_rhs(nm)) = nb_dof_b + var_problem.GetOffsetDofUnknown(n+m0);
                            nb_f_rhs(nm)++;
                          }
                        
                        offset_f_ += nb_comp_mesh(nm);
                      }                      
                  }
                
                if (f(ns)->IsNonNullSurfacicSourceGradient(ref))
                  {
                    int offset_df_ = 0;
                    num_mesh_geval(ns).Reallocate(nb_comp_df(ns));
                    num_source_geval(ns).Reallocate(nb_comp_df(ns));
                    pos_source_geval(ns).Reallocate(nb_comp_df(ns), nb_points_quadrature_boundary);
                    for (int n = 0; n < nb_unknowns; n++)
                      {
                        int nm = var_problem.mesh_num_unknown(n);
                        int N = nb_grad_mesh(nm);
                        for (int k = 0; k < N; k++)
                          {
                            num_mesh_geval(ns)(offset_df_ + k) = nm;
                            num_source_geval(ns)(offset_df_ + k) = nb_df_rhs(nm);
                            for (int j = 0; j < nb_points_quadrature_boundary; j++)
                              pos_source_geval(ns)(offset_df_+k, j) = N*j + k;
                          }
                        
                        offset_drhs(nm)(nb_df_rhs(nm)) = nb_dof_b + var_problem.GetOffsetDofUnknown(n);
                        offset_df_ += N;
                        nb_df_rhs(nm)++;
                      }                      
                  }
              }            
            nb_dof_b += var_problem.GetNbDof();
          }
        
	// and dof numbers
	for (int nm = 0; nm < nb_mesh_num; nm++)
          Nodle(nm) = var_problem.GetDofNumberOnElement(num_elem, nm);
	
	// we compute projection of quadrature points on real face, ds, and courbure
	Fb.FjElemQuadrature(s, PointsElemOrig, var_problem.mesh, num_elem);
	Fb.DFjElemQuadrature(s, PointsElemOrig, MatricesElemOrig, var_problem.mesh, num_elem);
	
	Fb.FjSurfaceElem(s, PointsElemOrig, var_problem.mesh, num_elem, num_loc);		
	Fb.DFjSurfaceElem(s, PointsElemOrig, MatricesElemOrig, var_problem.mesh, num_elem, num_loc);
	
	for (int nsrc = 0; nsrc < nb_modes; nsrc++)
	  {
	    PointsElem = PointsElemOrig;
	    MatricesElem = MatricesElemOrig;
	    
	    if (f(0) != NULL)
	      f(0)->ModifyPoints(nsrc, s, PointsElem, MatricesElem);

            if (source_phi)
              {
                for (int nm = 0; nm < nb_mesh_num; nm++)
                  {
                    feval(nsrc)(nm).Reallocate(nb_f_rhs(nm));
                    feval_hdg(nsrc)(nm).Reallocate(nb_f_hdg(nm));
                    for (int k = 0; k < nb_f_rhs(nm); k++)
                      {
                        feval(nsrc)(nm)(k).Reallocate(nb_comp_mesh(nm)*nb_points_quadrature_boundary);
                        feval(nsrc)(nm)(k).Zero();
                      }

                    for (int k = 0; k < nb_f_hdg(nm); k++)
                      {
                        feval_hdg(nsrc)(nm)(k).Reallocate(nb_comp_mesh(nm)*nb_points_quadrature_boundary);
                        feval_hdg(nsrc)(nm)(k).Zero();
                      }
                  }
              }

            if (source_grad)
              {
                for (int nm = 0; nm < nb_mesh_num; nm++)
                  {
                    feval_diff(nsrc)(nm).Reallocate(nb_df_rhs(nm));
                    for (int k = 0; k < nb_df_rhs(nm); k++)
                      {
                        feval_diff(nsrc)(nm)(k).Reallocate(nb_grad_mesh(nm)*nb_points_quadrature_boundary);
                        feval_diff(nsrc)(nm)(k).Zero();
                      }
                  }
              }
            
	    for (int ns = 0; ns < f.GetM(); ns++)
	      if (f(ns) != NULL)
		{
		  //int nb_unknowns = f(ns)->GetNbUnknowns();
		  Vector<T> f_eval(nb_comp_f(ns)), g_eval(nb_comp_df(ns));
		  f_eval.Zero(); g_eval.Zero();
                  
		  // computation of g at quadrature points
		  if (f(ns)->IsNonNullSurfacicSource(ref))
		    for (int j = 0; j < nb_points_quadrature_boundary; j++)
		      {
			f(ns)->EvaluateSurfacicSource(j, PointsElem, MatricesElem, f_eval);
                        
                        for (int k = 0; k < f_eval.GetM(); k++)
			    {
                              int nm = num_mesh_feval(ns)(k);
                              int ks = num_source_feval(ns)(k);
                              int pos = pos_source_feval(ns)(k, j);
                              bool hdg = hdg_source_feval(ns)(k);
                              if (hdg)
                                feval_hdg(nsrc)(nm)(ks)(pos) = f_eval(k);
                              else
                                feval(nsrc)(nm)(ks)(pos) = f_eval(k); 
			    }
		      }
		  
		  // computation of h at quadrature points
		  if (f(ns)->IsNonNullSurfacicSourceGradient(ref))
		    for (int j = 0; j < nb_points_quadrature_boundary; j++)
		      {
			f(ns)->EvaluateSurfacicSourceGradient(j, PointsElem,
							      MatricesElem, g_eval);
                        
                        for (int k = 0; k < g_eval.GetM(); k++)
                          {
                            int nm = num_mesh_geval(ns)(k);
                            int ks = num_source_geval(ns)(k);
                            int pos = pos_source_geval(ns)(k, j);
                            feval_diff(nsrc)(nm)(ks)(pos) = g_eval(k); 
                          }
		      }
		}
	  }
        
	if (f(0) != NULL)
	  f(0)->ModifyEvaluationSurface(source_phi, source_grad, feval, feval_diff,
                                        type_vec_f, type_vec_df);
        
	Vector<T> contrib_dof; //int rf(0), rot(0);
	
	for (int nsrc = 0; nsrc < b_source.GetM(); nsrc++)
	  {
	    if (source_phi)
	      {
		if (hdg_source)
		  {
                    for (int nm = 0; nm < nb_mesh_num; nm++)
                      {
                        const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
                        
                        const ElementReference_Dim<typename Dimension::DimensionBoundary>&
                          Fb_s = var_problem.GetSurfaceFiniteElement(num_face, nm);
                        
                        int rf = Fb_s.GetOrder();
                        if (rf != Fb.GetOrder())
                          {
                            cout << "Case not implemented" << endl;
                            cout << "there is not the same order for the surface and the element" << endl;
                            abort();
                          }
                        
                        int rot = var_problem.mesh.Element(num_elem).GetOrientationBoundary(num_loc);
                        
                        switch (type_u(nm))
                          {
                          case 1:
                            var_problem.GetReferenceElementH1(num_elem, nm).
                              ComputeIntegralSurfaceHDG(MatricesElemOrig, feval_hdg(nsrc)(nm),
                                                        contrib_surf, mesh_num,
                                                        num_elem, num_loc);
                            break;
                          case 2:
                            var_problem.GetReferenceElementHcurl(num_elem, nm).
                              ComputeIntegralSurfaceHDG(MatricesElemOrig, feval_hdg(nsrc)(nm),
                                                        contrib_surf, mesh_num,
                                                        num_elem, num_loc);
                            break;
                          case 3:
                            var_problem.GetReferenceElementHdiv(num_elem, nm).
                              ComputeIntegralSurfaceHDG(MatricesElemOrig, feval_hdg(nsrc)(nm),
                                                        contrib_surf, mesh_num,
                                                        num_elem, num_loc);
                            break;
                          }
                        
			int offset_loc = 0;
			for (int k = 0; k < num_loc; k++)
			  {
			    int nf = var_problem.mesh.Element(num_elem).numBoundary(k);
			    offset_loc += var_problem.GetSurfaceFiniteElement(nf).GetNbDof();
			  }
                        
                        contrib_dof.Reallocate(Fb_s.GetNbDof());
                        for (int m = 0; m < nb_f_hdg(nm); m++)
                          {
                            mesh_num.number_map.
                              ApplyOperatorFace(SeldonTrans, rf, rot,
                                                var_problem.mesh.Boundary(num_face),
                                                contrib_surf(m), contrib_dof);
                            
                            for (int k = 0; k < Fb.GetNbDofBoundary(num_loc); k++)
                              {
                                int num_dof_loc = mesh_num.Element(num_elem).GetNumberDof(offset_loc + k);
                                int num_dof = num_dof_loc + offset_hdg(nm)(m);
                                if (num_dof_loc >= 0)
                                  b_source(nsrc)(num_dof) += alpha*contrib_dof(k);
                              }
                          }
                      }
		  }
                
                for (int nm = 0; nm < nb_mesh_num; nm++)
                  {
                    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
                    int nb_dof_loc = var_problem.GetReferenceElement(num_elem, nm).GetNbDof();
                    
                    switch (type_u(nm))
                      {
                      case 1:
                        var_problem.GetReferenceElementH1(num_elem, nm).
                          ComputeIntegralSurface(MatricesElemOrig, feval(nsrc)(nm), contrib,
                                                 mesh_num, num_elem, num_loc);
                        break;
                      case 2:
                        var_problem.GetReferenceElementHcurl(num_elem, nm).
                          ComputeIntegralSurface(MatricesElemOrig, feval(nsrc)(nm), contrib,
                                                 mesh_num, num_elem, num_loc);
                        break;
                      case 3:
                        var_problem.GetReferenceElementHdiv(num_elem, nm).
                          ComputeIntegralSurface(MatricesElemOrig, feval(nsrc)(nm), contrib,
                                                 mesh_num, num_elem, num_loc);
                        break;
                      }

                    for (int m = 0; m < nb_f_rhs(nm); m++)
                      {
                        // computation on the unit face of \int_\hat{\Sigma} f \phi_i 
                        for (int k = 0; k < nb_dof_loc; k++)		  
                          {
                            int num_dof_loc = Nodle(nm)(k);
                            int num_dof = num_dof_loc + offset_rhs(nm)(m);
                            if (num_dof_loc >= 0)
                              b_source(nsrc)(num_dof) += alpha*contrib(m)(k);
                          }
                      }
		  }			  
	      }
	
	    if (source_grad)
	      {
                for (int nm = 0; nm < nb_mesh_num; nm++)
                  {
                    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
                    int nb_dof_loc = var_problem.GetReferenceElement(num_elem, nm).GetNbDof();
                    
                    // computation on the unit face of \int_\hat{\Sigma} h \nabla_\Sigma \phi_i 
                    switch (type_u(nm))
                      {
                      case 1:
                        var_problem.GetReferenceElementH1(num_elem, nm).
                          ComputeIntegralSurfaceGradient(MatricesElem, feval_diff(nsrc)(nm), contrib,
                                                         mesh_num, num_elem, num_loc);
                        break;
                      case 2:
                        var_problem.GetReferenceElementHcurl(num_elem, nm).
                          ComputeIntegralSurfaceGradient(MatricesElem, feval_diff(nsrc)(nm), contrib,
                                                         mesh_num, num_elem, num_loc);
                        break;
                      case 3:
                        var_problem.GetReferenceElementHdiv(num_elem, nm).
                          ComputeIntegralSurfaceGradient(MatricesElem, feval_diff(nsrc)(nm), contrib,
                                                         mesh_num, num_elem, num_loc);
                        break;
                      }

                    for (int m = 0; m < nb_f_rhs(nm); m++)
                      for (int k = 0; k < nb_dof_loc; k++)
                        {
                          int num_dof_loc = Nodle(nm)(k);
                          int num_dof = num_dof_loc + offset_rhs(nm)(m);
                          if (num_dof_loc >= 0)
                            b_source(nsrc)(num_dof) += alpha*contrib(m)(k);
                        }
                  }		
              }	    
          }    
      }
  }
  
  
  //! adding a Dirac source
  /*!
    \param[in,out] b_source right hand side to modify
    \param[in] f specification of the source 
    This method adds to b_source \f$ \alpha \varphi_i(x_0) \f$
  */
  template<class Dimension> template<class T>
  void VarSourceProblem_Fem<Dimension>
  ::AddDiracSourceGen(const T& alpha, Vector<Vector<T> > & b_source,
		      Vector<VirtualSourceFEM<T, Dimension>* >& f) const
  {
    VectR_N points_source;
    Vector<Vector<T> > polar;
    IVect point_number;
    for (int k = 0; k < f.GetM(); k++)
      if (f(k) != NULL)
	if (f(k)->IsDiracSource())
	  {	    
	    points_source.PushBack(f(k)->GetOrigin());
	    if (f(k)->IsGradientDirac())
              polar.PushBack(f(k)->GetPolarizationGrad());
            else
              polar.PushBack(f(k)->GetPolarization());
            
	    point_number.PushBack(k);
	  }
    
    if (point_number.GetM() > 0)
      AddDiracSourceGen(alpha, points_source, polar, point_number, b_source, f);
  }

  
  //! Dirac for axisymmetric equations
  template<class Dimension> template<class T>
  void VarSourceProblem_Fem<Dimension>::
  AddDiracSourceAxisym(const T& alpha, const Vector<bool>& vec_unknown, int m,
                       Vector<Vector<T> > & b_source, const VectR3& points,
		       Vector<VirtualSourceFEM<T, Dimension>* >& f) const
  {
    VectR2 points_source;    
    Vector<Vector<T> > polar_source;
    IVect point_number;

    for (int k = 0; k < f.GetM(); k++)
      if (f(k) != NULL)
	if (f(k)->IsDiracSource())
	  {	        
	    R3 pt = points(k);
	    
	    // conversion in cylindrical coordinates of the origin
	    Real_wp r0, z0(pt(2)), theta0, cos_theta0, sin_theta0;
	    CartesianToPolar(pt(0), pt(1), r0, theta0);
	    cos_theta0 = cos(theta0); sin_theta0 = sin(theta0);
            
	    R2 origin(r0, z0);
    
	    // checking if the point belongs to the axis Oz
	    bool point_on_axis = false;
	    if (r0 <= epsilon_machine)
	      point_on_axis = true;
	    
	    // beta = 1/(2 pi) if the point is on the axis, exp(im theta0) / (2pi) otherwise
	    T beta = 1.0/(2.0*pi_wp);
	    if (!point_on_axis)
	      to_complex(beta*exp(Iwp*Real_wp(m)*theta0), beta);    
	    
	    Vector<T> polar = f(k)->GetPolarization();    
	    Vector<T> coef(polar.GetM());
	    coef.Fill(0);
	    
	    int nb = 0;
	    //bool null_source = point_on_axis;
	    while (nb < coef.GetM())
	      {
		if (vec_unknown(nb))
		  {
		    if (nb == coef.GetM()-1)
		      {
			cout << "Give an even number of vectorial unknowns" << endl;
			abort();
		      }
		    
		    Real_wp ux = realpart(polar(nb)), uy = realpart(polar(nb+1));
		    if (point_on_axis)
		      {
			if (abs(m) == 1)
			  {
			    //null_source = false; 
			    if (m == -1)
			      {
				to_complex(0.5*beta*Complex_wp(ux, -uy), coef(nb));
				to_complex(0.5*beta*Complex_wp(uy, ux), coef(nb+1));
			      }
			    else
			      {
				to_complex(0.5*beta*Complex_wp(ux, uy), coef(nb));
				to_complex(0.5*beta*Complex_wp(uy, -ux), coef(nb+1));
			      }
			  }
		      }
		    else
		      {
			coef(nb) = beta*(ux*cos_theta0 + uy*sin_theta0);
			coef(nb+1) = beta*(-ux*sin_theta0 + uy*cos_theta0);
		      }
		    
		    nb += 2;
		  }
		else
		  {
		    Real_wp uz = realpart(polar(nb));
		    if (point_on_axis)
		      {
			if (m == 0)
			  {
			    //null_source = false; 
			    coef(nb) = beta*uz;
			  }
		      }
		    else
		      {
			coef(nb) = beta*uz;
		      }
		    
		    nb++;
		  }
	      }
	    
            // pour eviter les decalages, on calcule meme si la source est nulle
	    //if (!null_source)
            //{
            points_source.PushBack(origin);
            polar_source.PushBack(coef);
            point_number.PushBack(k);
            //}
	  }
    
    //DISP(points_source);
    //DISP(polar_source);
    //DISP(point_number);
    
    if (point_number.GetM() > 0)
      AddDiracSourceGen(alpha, points_source, polar_source, point_number, b_source, f);    
  }
  
  
  //! adding a Dirac source
  /*!
    \param[in,out] b_source right hand side to modify
    \param[in] f specification of the source 
    This method adds to b_source \f$ \alpha \varphi_i(x_0) \f$
  */
  template<class Dimension> template<class T>
  void VarSourceProblem_Fem<Dimension>
  ::AddDiracSourceGen(const T& alpha, const VectR_N& points_source,
		      const Vector<Vector<T> >& polar,
		      const IVect& point_number, Vector<Vector<T> > & b_source,
		      Vector<VirtualSourceFEM<T, Dimension>* >& f) const
  {
    if (points_source.GetM() <= 0)
      return;

    int nb_source = 0;
    IVect offset_source_num(f.GetM()); offset_source_num.Zero();
    for (int k = 0; k < f.GetM(); k++)
      {
	offset_source_num(k) = nb_source;
	if (f(k) != NULL)
	  nb_source += f(k)->GetNbUnknowns();
      }

    if (nb_source <= 0)
      return;

    int m0 = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      m0 = var_problem.nb_unknowns_hdg;

    int nb_mesh_num = var_problem.GetNbMeshNumberings();
    Vector<IVect> Nodle(nb_mesh_num);
    IVect type_u(nb_mesh_num);
    IVect nb_comp_mesh(nb_mesh_num);
    IVect nb_grad_mesh(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      {
        int nb_u = var_problem.GetNbComponentsUnknown(nm);
        int nb_du = var_problem.GetNbComponentsGradient(nm);
        nb_comp_mesh(nm) = nb_u; nb_grad_mesh(nm) = nb_du;
        if (nm == 0)
          type_u(nm) = var_problem.type_element;
        else
          type_u(nm) = var_problem.other_type_element(nm-1);
      }
    
    IVect offset_rhs(nb_source); offset_rhs.Zero();
    int nb_dof_b = 0; nb_source = 0;
    for (int k = 0; k < f.GetM(); k++)
      {
	if (f(k) != NULL)
	  for (int m = m0; m < f(k)->GetNbUnknowns()+m0; m++)
	    offset_rhs(nb_source++) = nb_dof_b + var_problem.GetOffsetDofUnknown(m);
	
	nb_dof_b += var_problem.GetNbDof();
      }

    Vector<Vector<T> > b_add(b_source.GetM());
    for (int i = 0; i < b_source.GetM(); i++)
      {
	b_add(i).Reallocate(b_source(i).GetM());
	b_add(i).Zero();
      }
    
    enum {nb_u_hcurl = FiniteElementHcurl<Dimension>::nb_components_u,
          nb_du_hcurl = FiniteElementHcurl<Dimension>::nb_components_grad };

    Vector<Vector<TinyVector<Real_wp, 1> > > val_H1(nb_mesh_num);    
    Vector<Vector<TinyVector<Real_wp, Dimension::dim_N> > > val_Hdiv(nb_mesh_num);
    Vector<Vector<TinyVector<Real_wp, nb_u_hcurl> > > val_Hcurl(nb_mesh_num);

    Vector<Vector<TinyVector<Real_wp, Dimension::dim_N> > > grad_H1(nb_mesh_num);    
    Vector<Vector<TinyVector<Real_wp, 1> > > div_Hdiv(nb_mesh_num);
    Vector<Vector<TinyVector<Real_wp, nb_du_hcurl> > > curl_Hcurl(nb_mesh_num);

    // then loop over  all elements
    SetMatrices<Dimension> MatricesElem;
    TinyVector<R_N, 2> enveloppe; R_N pt_loc; 
    IVect nb_eval(points_source.GetM());
    nb_eval.Zero();
    for (int num_elem = 0; num_elem < var_problem.mesh.GetNbElt(); num_elem++)
      {
	const ElementReference_Dim<Dimension>& Fb = var_problem.GetReferenceElement(num_elem);
	
	// initialization of F_i^-1
	FjInverseProblem<Dimension> inverseFj(var_problem.mesh, num_elem);

	SetPoints<Dimension>& PointsElem = inverseFj.GetSetPoints();
	VectR_N& s = inverseFj.GetVertices();
	
	var_problem.mesh.GetBoundingBox(num_elem, s, PointsElem, enveloppe);
	for (int nm = 0; nm < nb_mesh_num; nm++)
          Nodle(nm) = var_problem.GetDofNumberOnElement(num_elem, nm);
	
	// loop over points to localize
	for (int k = 0; k < points_source.GetM(); k++)
	  if (PointInsideBoundingBox(points_source(k), enveloppe))
	    {
	      R_N pt_source = points_source(k);
	      int ns = point_number(k);
	      int ks = offset_source_num(ns);
	      int nb_unknowns = f(ns)->GetNbUnknowns();
              bool grad_dirac = f(ns)->IsGradientDirac();
	      bool pt_inside = inverseFj.Solve(pt_source, pt_loc);
	      Real_wp distance_boundary_elt
		= var_problem.mesh.GetDistanceToBoundary(pt_loc, num_elem);
	      
	      if (distance_boundary_elt < 0)
		var_problem.mesh.ProjectPointOnBoundary(pt_loc, num_elem);
	      
	      typename Dimension::MatrixN_N dfj, dfjm1;
	      if (pt_inside)
		{
		  Fb.DFj(s, PointsElem, pt_loc, dfj, var_problem.mesh, num_elem);
		  GetInverse(dfj, dfjm1);

                  if (grad_dirac)
                    {
                      // we compute derivative of basis functions on the given point
                      for (int nm = 0; nm < nb_mesh_num; nm++)
                        switch (type_u(nm))
                          {
                          case 1:
                            var_problem.GetReferenceElementH1(num_elem, nm)
                              .ComputeValuesGradientPhi(pt_loc, grad_H1(nm), dfjm1, var_problem.GetMeshNumbering(nm), num_elem);
                            break;
                          case 2:
                            var_problem.GetReferenceElementHcurl(num_elem, nm)
                              .ComputeValuesGradientPhi(pt_loc, curl_Hcurl(nm), dfjm1, var_problem.GetMeshNumbering(nm), num_elem);
                            break;
                          case 3:
                            var_problem.GetReferenceElementHdiv(num_elem, nm)
                              .ComputeValuesGradientPhi(pt_loc, div_Hdiv(nm), dfjm1, var_problem.GetMeshNumbering(nm), num_elem);
                            break;
                          }
                      
                      T value; int offset = 0;
                      for (int m = 0; m < nb_unknowns; m++)
                        {
                          int nm = var_problem.mesh_num_unknown(m%var_problem.nb_unknowns);
                          int nb_dof = var_problem.GetReferenceElement(num_elem, nm).GetNbDof();
                          switch (type_u(nm))
                            {
                            case 1:
                              {
                                TinyVector<T, Dimension::dim_N> coef;
                                for (int p = 0; p < Dimension::dim_N; p++)
                                  coef(p) = polar(k)(offset+p);
                                
                                coef *= f(ns)->GetCoefficientVolume();
                                for (int i = 0; i < nb_dof; i++)
                                  {
                                    int num_dof = Nodle(nm)(i);
                                    if (num_dof >= 0)
                                      {
                                        value = DotProd(coef, grad_H1(nm)(i));
                                        b_add(0)(offset_rhs(m+ks) + num_dof) += alpha*value;
                                      }
                                  }
                              }
                              break;
                            case 2:
                              {
                                TinyVector<T, nb_du_hcurl> coef;
                                for (int p = 0; p < nb_du_hcurl; p++)
                                  coef(p) = polar(k)(offset+p);
                                
                                coef *= f(ns)->GetCoefficientVolume();
                                for (int i = 0; i < nb_dof; i++)
                                  {
                                    int num_dof = Nodle(nm)(i);
                                    if (num_dof >= 0)
                                      {
                                        value = DotProd(coef, curl_Hcurl(nm)(i));
                                        b_add(0)(offset_rhs(m+ks) + num_dof) += alpha*value;
                                      }
                                  }
                              }
                              break;
                            case 3:
                              {
                                T coef = polar(k)(offset)*f(ns)->GetCoefficientVolume();
                                for (int i = 0; i < nb_dof; i++)
                                  {
                                    int num_dof = Nodle(nm)(i);
                                    if (num_dof >= 0)
                                      {
                                        value = coef*div_Hdiv(nm)(i)(0);
                                        b_add(0)(offset_rhs(m+ks) + num_dof) += alpha*value;
                                      }
                                  }
                              }
                              break;
                            }
                          
                          offset += nb_grad_mesh(nm);
                        }
                    }
                  else
                    {
                      // we compute basis functions on the given point
                      for (int nm = 0; nm < nb_mesh_num; nm++)
                        switch (type_u(nm))
                          {
                          case 1:
                            var_problem.GetReferenceElementH1(num_elem, nm)
                              .ComputeValuesPhi(pt_loc, val_H1(nm), dfjm1, var_problem.GetMeshNumbering(nm), num_elem);
                            break;
                          case 2:
                            var_problem.GetReferenceElementHcurl(num_elem, nm)
                              .ComputeValuesPhi(pt_loc, val_Hcurl(nm), dfjm1, var_problem.GetMeshNumbering(nm), num_elem);
                            break;
                          case 3:
                            var_problem.GetReferenceElementHdiv(num_elem, nm)
                              .ComputeValuesPhi(pt_loc, val_Hdiv(nm), dfjm1, var_problem.GetMeshNumbering(nm), num_elem);
                            break;
                          }
                      
                      T value; int offset = 0;
                      for (int m = 0; m < nb_unknowns; m++)
                        {
                          int nm = var_problem.mesh_num_unknown(m%var_problem.nb_unknowns);
                          int nb_dof = var_problem.GetReferenceElement(num_elem, nm).GetNbDof();
                          switch (type_u(nm))
                            {
                            case 1:
                              {
                                T coef = polar(k)(offset)*f(ns)->GetCoefficientVolume();
                                for (int i = 0; i < nb_dof; i++)
                                  {
                                    int num_dof = Nodle(nm)(i);
                                    if (num_dof >= 0)
                                      {
                                        value = coef*val_H1(nm)(i)(0);
                                        b_add(0)(offset_rhs(m+ks) + num_dof) += alpha*value;
                                      }
                                  }
                              }
                              break;
                            case 2:
                              {
                                TinyVector<T, nb_u_hcurl> coef;
                                for (int p = 0; p < nb_u_hcurl; p++)
                                  coef(p) = polar(k)(offset+p);
                                
                                coef *= f(ns)->GetCoefficientVolume();
                                for (int i = 0; i < nb_dof; i++)
                                  {
                                    int num_dof = Nodle(nm)(i);
                                    if (num_dof >= 0)
                                      {
                                        value = DotProd(coef, val_Hcurl(nm)(i));
                                        b_add(0)(offset_rhs(m+ks) + num_dof) += alpha*value;
                                      }
                                  }
                              }
                              break;
                            case 3:
                              {
                                TinyVector<T, Dimension::dim_N> coef;
                                for (int p = 0; p < Dimension::dim_N; p++)
                                  coef(p) = polar(k)(offset+p);
                                
                                coef *= f(ns)->GetCoefficientVolume();
                                for (int i = 0; i < nb_dof; i++)
                                  {
                                    int num_dof = Nodle(nm)(i);
                                    if (num_dof >= 0)
                                      {
                                        value = DotProd(coef, val_Hdiv(nm)(i));
                                        b_add(0)(offset_rhs(m+ks) + num_dof) += alpha*value;
                                      }
                                  }
                              }
                              break;
                            }
                          
                          offset += nb_comp_mesh(nm);
                        }
                    }
                  
		  nb_eval(k)++;
		}
	    }
      }
    
    IVect nb_eval_all(nb_eval);
#ifdef SELDON_WITH_MPI
    MPI_Allreduce(nb_eval.GetData(), nb_eval_all.GetData(),
                  nb_eval.GetM(), MPI_INTEGER, MPI_SUM, var_problem.comm_group_mode);
#endif    
    
    for (int k = 0; k < nb_eval_all.GetM(); k++)
      if (nb_eval_all(k) > 1)
	{
	  Real_wp coef = Real_wp(1)/nb_eval_all(k);
	  int ks = point_number(k);
	  int Nvol = var_problem.offset_dof_unknown(var_problem.nb_unknowns);    
	  int offset = offset_rhs(offset_source_num(ks));
	  for (int i = 0; i < b_add.GetM(); i++)
	    for (int j = 0; j < Nvol; j++)
	      b_add(i)(offset + j) *= coef;
	}
    
    T one; SetComplexOne(one);
    for (int i = 0; i < b_source.GetM(); i++)
      Add(one, b_add(i), b_source(i));
  }
  
  
  //! computation of a generic (volumic, surfacic and inhomogeneous dirichlet)
  /*!
    \param[out] b_source right hand side to compute
    \param[in] f specification of the source
   */
  template<class Dimension> template<class T>
  void VarSourceProblem_Fem<Dimension>
  ::ComputeGenericSource(Vector<Vector<T> > & b_source,
			 Vector<VirtualSourceFEM<T, Dimension>* >& f, bool assemble) const
  {
    // allocation of source vectors
    int nb_rhs = f.GetM();
    for (int nsrc = 0; nsrc < b_source.GetM(); nsrc++)
      {
	b_source(nsrc).Reallocate(var_problem.GetNbDof()*nb_rhs);
	b_source(nsrc).Zero();
      }
    
    // different type of sources -> volumic source, surfacic sources, dirichlet sources
    T one; SetComplexOne(one);
    
    // volumic part of the right hand side
    var_source.AddVolumeSource(one, b_source, f);
    
    // for neumann boundaries and other Robin-like conditions
    var_source.AddSurfaceSource(one, b_source, f);
    
    // Dirac
    this->AddDiracSource(one, b_source, f);

    // then inhomogeneous Dirichlet condition
    this->SetDirichletSource(b_source, f);
    
    // assembling right hand side
    if (assemble)
      for (int nsrc = 0; nsrc < b_source.GetM(); nsrc++)
        var_problem.AddDomains(b_source(nsrc));    
  }


  /********************
   * VarSourceProblem *
   ********************/
  

  //! computation of the right hand side
  /*!
    \param[out] b_rhs right hand side computed
   */
  template<class Complexe, class Dimension>
  void VarSourceProblem<Complexe, Dimension>
  ::ComputeRightHandSideGen(VectComplexe & b_rhs, bool assemble)
  {
    Vector<VectComplexe> vec_rhs;
    vec_rhs.SetData(1, &b_rhs);
    
    ComputeRightHandSide(vec_rhs, assemble);
    
    vec_rhs.Nullify();
  }
  

  template<class Complexe, class Dimension>
  void VarSourceProblem<Complexe, Dimension>
  ::InitIncidentField()
  {
    int nb_rhs = this->GetNbRhs();
    
    // reallocation of incident arrays
    this->ClearIncidentArray();
    this->incident_wave_proj.Reallocate(nb_rhs);
    this->incident_wave_field.Reallocate(nb_rhs);
    this->incident_wave_proj.Fill((IncidentWaveProjector<Complexe, Dimension>*) NULL);
    this->incident_wave_field.Fill((IncidentWaveField<Complexe, Dimension>*) NULL);
    
    // incident fields and projectors are initialized
    for (int n = 0; n < nb_rhs; n++)
      {
	if ( (this->type_source(n) == VarSourceProblem_Base::SRC_DIFFRACTED_FIELD)
	     || (this->type_source(n) == VarSourceProblem_Base::SRC_TOTAL_FIELD) )
	  {
	    this->incident_wave_field(n)
	      = this->GetNewIncidentField(this->incident_source(n),
					  this->source_space_param(n), Complexe(0));
	    
	    this->incident_wave_proj(n)
	      = this->GetNewIncidentProjector(this->incident_source(n),
					      this->source_space_param(n),
					      *this->incident_wave_field(n));
	  }
      }
  }

  
  //! computation of the right hand side  
  template<class Complexe, class Dimension>
  void VarSourceProblem<Complexe, Dimension>
  ::ComputeRightHandSide(Vector<VectComplexe> & b_rhs, bool assemble)
  { 
    int nb_rhs = this->GetNbRhs();
    Vector<VirtualSourceFEM<Complexe, Dimension>* > f(nb_rhs);
    
    InitIncidentField();
    
    for (int n = 0; n < nb_rhs; n++)
      {
	// object describing the source n is created
	f(n) = this->GetNewSourceEquationObject(n);
      }
    
    // right hand side are computed together
    this->ComputeGenericSource(b_rhs, f, assemble);
    
    // then temporary objects stored in f are destructed
    for (int n = 0; n < nb_rhs; n++)
      {
	if (f(n) != NULL)
	  delete f(n);
      }
  }

  
  //! computes right hand sides
  template<class Complexe, class Dimension>
  void VarSourceProblem<Complexe, Dimension>::
  ComputeRightHandSide(Matrix<Complexe, General, ColMajor> & b_source, bool assemble)
  {
    Vector<VectComplexe> rhs(1);
    b_source.Reallocate(var_problem.GetNbDof(), this->GetNbRhs());
    rhs(0).SetData(b_source.GetDataSize(), b_source.GetData());
    
    ComputeRightHandSide(rhs, assemble);
    
    rhs.Nullify();
  }
  

  //! we add to U0 an incident wave alpha*uinc
  /*!
    \param[in] alpha 
    \param[in,out] U0
   */
  template<class Complexe, class Dimension>
  void VarSourceProblem<Complexe, Dimension>::
  AddIncidentWaveGen(const Complexe& alpha, Vector<Vector<Complexe> > & U0)
  {
    int N = this->incident_wave_proj.GetM();
    Vector<VirtualProjectorFEM<Complexe, Dimension>* > f(N);
    for (int i = 0; i < N; i++)
      f(i) = this->incident_wave_proj(i);
    
    this->AddVolumetricProjection(alpha, U0, f);
  }
  
  
  //! copies input parameters of another similar problem
  template<class Complexe, class Dimension>
  void VarSourceProblem<Complexe, Dimension>
  ::CopyInputData(const VarSourceProblem_Base& var_base)
  {
    VarSourceProblem_Fem<Dimension>::CopyInputData(var_base);
    
    const VarSourceProblem<Complexe, Dimension>& var
      = static_cast<const VarSourceProblem<Complexe, Dimension>& >(var_base);

    this->coefficient_volumic_source = var.coefficient_volumic_source;
  }
    
} // end namespace

#define MONTJOIE_FILE_DEFINE_SOURCE_ELLIPTIC_CXX
#endif  
