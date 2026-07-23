#ifndef MONTJOIE_FILE_CROSS_COVARIANCE_CXX

#include "Helio/CrossCovariance.hxx"

namespace Montjoie
{
  
  /************
   * MySource *
   ************/
  
  MySource::MySource(const R3& pt, const Real_wp& om, int m0, const Vector<Real_wp> dirac)
  {
    pt_source = pt;
    omega = om;
    m = m0;
    diracAmp = dirac;
    Lmax  = diracAmp.GetM();
    assocLegendre.Init(Lmax);
    Real_wp thetaSrc = atan2(pt_source(0),pt_source(2));
    assocLegendre.EvaluatePnm(Lmax,abs(m0),thetaSrc,Psrc);
  }

  
  void MySource::EvaluateFunction(const R3& x, Vector<Complex_wp>& f) const
  {
    //Real_wp alpha = 7.0 / square(radius);
    //f(0) = exp(-alpha*x.DistanceSquare(pt_source));
    int     am    = abs(m);
    Real_wp theta = atan2(x(0),x(2));
    Vector<VectReal_wp> P;
    assocLegendre.EvaluatePnm(Lmax,am,theta,P);

    f(0) = 0.e0;
    for (int l=am;l<Lmax;l++)
      f(0) += diracAmp(l) *P(l)(am) *Psrc(l)(am);

  } 
 

  RandomSource::RandomSource(const Real_wp& om, int m0, int lmax, string sourcePowerType, Real_wp ts, Vector<Real_wp> depthSources, Vector<Real_wp> deltaDepthSources, Vector<Real_wp> sourcePowerDepth)
  {
    omega = om;
    m = m0;
    Lmax  = lmax;
    assocLegendre.Init(lmax);
    source_power_type = sourcePowerType;
    Ts = ts;
    source_power_freq = GetSourcePower(omega);

    depth_sources = depthSources;
    delta_depth_sources = deltaDepthSources;
    source_power_depth = sourcePowerDepth;

    realisations.Reallocate(depth_sources.GetM());
    for (int i=0; i < depth_sources.GetM(); i++)
      {
	realisations(i).Reallocate(lmax);
	realisations(i).Zero();
	realisations(i).FillRand();
      }
  }

  Real_wp RandomSource::GetSourcePower(const Real_wp& omega)
  {
    if (source_power_type == "EXP")
      return exp(-(Ts*omega)*(Ts*omega));
    else
      {
	cout << "This type of damping is not implemented yet" << endl;
	abort();
      }
  }

  
  void RandomSource::EvaluateFunction(const R3& x, Vector<Complex_wp>& f) const
  {
    Real_wp r = sqrt(x(0)*x(0)+x(1)*x(1)+x(2)*x(2));
    f(0) = 0.e0;
    for (int i=0; i < depth_sources.GetM(); i++)
      {
	for (int l=abs(m);l<Lmax;l++)
	  f(0) += realisations(i)(l) * exp(-( (r-depth_sources(i))*(r-depth_sources(i)) ) / (2.* delta_depth_sources(i)* delta_depth_sources(i)) ) * source_power_depth(i);
      }
    f(0) *= source_power_freq;
  } 
  
  
  
  /****************************
   * CrossCovarianceParameter *
   ****************************/

  
  CrossCovarianceParameter::CrossCovarianceParameter()
  {
    nSrc = 1;
    filteredSrc = false;
    var = NULL;
    var_multi = NULL;
    display_solution = false;
    damping_type = "CONSTANT";
    gammaOld = -1;
    RSUN = 696.e6;
  }
  

  void CrossCovarianceParameter
  ::SetPulsation(int m, int mg, const Real_wp& omega,
                 VarComputationProblem& var_pb, MultiFrequencyProblem& pb)
  {
    var = &dynamic_cast<EllipticProblem<HelmholtzEquationAxi>& >(var_pb);
    var_multi = &pb;

    omega_stored = omega;
    var->UpdateWaveVector(omega);
    // Update the damping
    // var->ref_sigma(1).Mlt(alpha) to multiply sigma by alpha
    if (damping_type != "CONSTANT")
      {
	Real_wp gamma = GetDamping(omega);
	if (gammaOld < 0) // first frequency
	  {
	    var->ref_sigma(1).Mlt(gamma);
	  }
	else
	  {
	    var->ref_sigma(1).Mlt(gamma/gammaOld); 
	  }
	gammaOld = gamma;
      }

    var->ComputeMassMatrix(false);  // false : on recalcule pas les indices



    // the coefficients are computed for the covariance
    CrossCovarianceProblem& var_cross = dynamic_cast<CrossCovarianceProblem&>(*var_multi);
    var_cross.ComputeMassMatrix(*var);
    
#ifdef SELDON_WITH_MPI
    comm_freq = pb.GetFrequencyCommunicator();
#endif
  }

  Real_wp CrossCovarianceParameter
  ::GetDamping(const Real_wp& omega)
  {
    if (damping_type == "POWER_LAW")
      return gamma0 * pow(abs(omega/omega0), beta);
    else if (damping_type == "PROPORTIONAL")
      return beta * abs(omega);
    else if (damping_type == "CONSTANT")
      return gamma0;
    else
      {
	cout << "This type of damping is not implemented yet" << endl;
	abort();
      }
  }


  int CrossCovarianceParameter
  ::GetRealM(int mg)
  {
    if (var->include_flow_term)
      {
	if (mg % 2 == 0)
	  return mg/2;
	else
	  return -(mg+1)/2;
      }
    else
      {
	return mg;
      }
  }  

  void CrossCovarianceParameter
  ::SetParameters(int k, int kg, int m, int mg,
                  VarComputationProblem& var_pb, MultiFrequencyProblem&)
  {    
    
    IVect num(1); num(0) = GetRealM(mg);
    /*
    if (var->include_flow_term)
      {
	if (mg % 2 == 0)
	  num(0) = mg/2;
	else
	  num(0) = -(mg+1)/2;
        //cout << "Case with flow not handled" << endl;
        //abort();
      }
    else
      {
	// cas sans ecoulement : on factorise que les modes positifs
	num(0) = mg;
      }
    */

    var->SetModesToCompute(num);  
    var->SetCurrentModeNumber(num(0)); 
    var->UpdateDirichlet(num(0));
    
  }


  void CrossCovarianceParameter
  ::InitSource(Vector<string> lines_data_file)
  {
    // The source data file contains FILTERED on the first line if filtered then the list of sources (one per line)

    /* VectString words;
    VectString filterString;
    for (int i = 0; i < lines_data_file.GetM(); i++)
      {
	StringTokenize(lines_data_file(i), words, string(" \t"));
	if (words(0) == "Filter")
	  {
	    filterString.Reallocate(words.GetM()-2);
	    for (int j=0; j < words.GetM()-2; j++)
	      filterString(j) = words(j+2);
	    CreateFilter(filterString);
	  }
	  } */


    /* nSrc           = srcLines.GetM()-1; // one line is added due to the end of file
    //filteredSrc    = (srcLines(0)=="FILTERED"); 
    //nSrc -= filteredSrc;

    pt_source.Reallocate(nSrc);
    //diracAmps.Reallocate(nSrc);
    for (int ID_src = 0; ID_src < nSrc; ID_src++)
      {
        int num_src = ID_src;
        //if (filteredSrc)
        //  num_src = ID_src+1;

        // read source location
        istringstream flux(srcLines(num_src));
        R3 pt; flux >> pt(0) >> pt(1) >> pt(2);
        pt_source(ID_src) = pt;

        if (filteredSrc)
          {
	    // Read and store the filter coefficients (first line is the number of l coefficients)
	    string fileSrc;
	    flux >> ref_surface_observation;
	    flux >> fileSrc;
	    Vector<string> sourceAmpLines;
	    ReadLinesFile(fileSrc,sourceAmpLines);
	    int nL = to_num<int>(sourceAmpLines(0));
	    int nLine = 1;
	    diracAmps(ID_src).Reallocate(nL+1); diracAmps(ID_src).Fill(0.e0);
	    for (int l=0; l<=nL; l++) {
	      sourceAmpLines(nLine).append(" ");
	      diracAmps(ID_src)(l) = to_num<Real_wp>(sourceAmpLines(nLine));
	      nLine++;
	    }
	    }
      }
    
    srcLines.Clear(); */
  }

  void CrossCovarianceParameter
  ::InitDamping(VectString dampingString)
  {
    damping_type = dampingString(0);
    if (dampingString(0) == "POWER_LAW")
      {
	gamma0 = 2.*pi_wp*to_num<Real_wp>(dampingString(1))*2.*RSUN;	
	omega0 = 2.*pi_wp*to_num<Real_wp>(dampingString(2))*RSUN;
	beta = to_num<Real_wp>(dampingString(3));
      }
    else if (dampingString(0) == "PROPORTIONAL")
      {
	beta = to_num<Real_wp>(dampingString(1));
      }
    else if (dampingString(0) == "CONSTANT")
      {
	gamma0 = to_num<Real_wp>(dampingString(1));
      }

    
  }

  void CrossCovarianceParameter
  ::CreateFilter(VectString filterString)
  {
    filteredSrc = true;
    if (filterString(0) == "LFILTER")
      {
	int nL = to_num<int>(filterString(1));
	string typeOfFilter = filterString(2);
	diracAmps.Reallocate(1);
	diracAmps(0).Reallocate(nL+1);
	if (typeOfFilter == "TANH")
	  {
	    int threshold = to_num<int>(filterString(3));
	    Real_wp slope = to_num<Real_wp>(filterString(4));
	    for (int l=0; l<=nL; l++)
	      diracAmps(0)(l) = 0.5*(1.e0-tanh(slope*(l-threshold)));	  
	  }
	else if (typeOfFilter == "GAUSSIAN")
	  {
	    int threshold = to_num<int>(filterString(3));
	    Real_wp slope = to_num<Real_wp>(filterString(4));
	    for (int l=0; l<=nL; l++)
	      diracAmps(0)(l) = exp(-square(l - threshold) / (2*square(slope)));	  
	  }
	else if (typeOfFilter == "LINEAR")
	  {
	    int threshold = to_num<int>(filterString(3));
	    for (int l=0; l<=nL; l++)
	      diracAmps(0)(l) = max(1.e0-l/(threshold*1.e0),0.e0);	  
	  }
	else if (typeOfFilter == "SINGLE")
	  {
	    int lsingle = to_num<int>(filterString(3));
	    diracAmps(0).Fill(0.);
	    diracAmps(lsingle) = 1.;	  
	  }
	else if (typeOfFilter == "STEP")
	  {
	    int threshold = to_num<int>(filterString(3));
	    diracAmps(0).Fill(0.);
	    for (int l=0; l<=threshold; l++)
	      diracAmps(0)(l) = 1.;	  
	  }
	else
	  {
	    cout << "This type of filter (" << typeOfFilter << ") is not known" << endl;
	    cout << "Choose between LINEAR, TANH, GAUSSIAN, SINGLE or STEP" << endl;
	    abort();
	  }

	if (filterString(filterString.GetM()-1) == "SQRT")
	  {
	    for (int l=0; l<diracAmps(0).GetM(); l++)
	      diracAmps(0)(l) = sqrt(diracAmps(0)(l));
	  }
      }

  }

  void CrossCovarianceParameter
  ::InitSourcePower(VectString sourcePowerString)
  {
    source_power_type = sourcePowerString(0);
    if (source_power_type == "EXP")
      {
	Ts = to_num<Real_wp>(sourcePowerString(1));
      }
  }  
 

  void CrossCovarianceParameter
  ::SetSourceFrequency(int k, int kg, int m, int mg, const VarComputationProblem& var_pb,
                       Matrix<Complex_wp, General, ColMajor>& sol)
  {
    // cas ou on prend la source du fichier .ini
    /* 
    // source pour mode positif
    VectComplex_wp source_pos, source_neg;
    var->ComputeRightHandSide(source_pos);

    if (mg != 0)
      {
        var->SetCurrentModeNumber(-mg);
        var->ComputeRightHandSide(source_neg);
        
        var->SetCurrentModeNumber(mg);
        
        sol.Reallocate(source_pos.GetM(), 2);
        for (int i = 0; i < source_pos.GetM(); i++)
          {
            sol(i, 0) = source_pos(i);
            sol(i, 1) = source_neg(i);
          }
      }
    else
      {
        sol.Reallocate(source_pos.GetM(), 1);
        SetCol(source_pos, 0, sol);
      }
    */

    // cas ou on prend une source Dirac avec un fichier
    // contenant la liste des points source
    // param_source contient les parametres additionels du Dirac (polarisation pour une equation vectorielle)
    Vector<VectString> param_source(nSrc);
    int nb_comp = 1;
    if (mg != 0 and !var->include_flow_term)
      nb_comp = 2;
    
    sol.Reallocate(var->GetNbDof(), nSrc*nb_comp);

    VectComplex_wp b_source(sol.GetDataSize()/nb_comp, sol.GetData());
    CrossCovarianceProblem& var_cross = dynamic_cast<CrossCovarianceProblem&>(*var_multi);

    if (!filteredSrc and typeSource != "RANDOM")
      var->ComputeDiracSource(b_source, pt_source, param_source);
    else
      {
        Vector<VirtualSourceFEM<Complex_wp, Dimension2>* > vec_source(nSrc);
        VirtualSourceField<Complex_wp, Dimension3>* f_surf;
        for (int k = 0; k < nSrc; k++)
          {
	    VolumetricSource<HelmholtzEquationAxi>* my_source;
	    my_source = new VolumetricSource<HelmholtzEquationAxi>(*var, param_source);
	    // ici f_surf peut dependre de mg (numero de mode) et kg (numero de frequence)
	    if (typeSource == "RANDOM")
	      {
		//VolumetricSource<TypeElement, RandomSource, TypeEquation>  my_source;
	        f_surf = new RandomSource(omega_stored, GetRealM(mg), lMax, source_power_type, Ts, depth_sources, delta_depth_sources, source_power_depth);
		my_source->SetVolumeSourceFunction(f_surf);
	      }
	    else
	      {
		f_surf = new MySource(pt_source(k),omega_stored, GetRealM(mg), diracAmps(0));
		my_source->SetSurfaceSource(var_cross.ref_surface_observation, f_surf);
	      }
	    vec_source(k) = my_source;

          }
        
        Vector<VectComplex_wp> b_vec;
        b_vec.SetData(1, &b_source);
        var->ComputeGenericSource(b_vec, vec_source, false);

        for (int k = 0; k < nSrc; k++)
	  delete vec_source(k);
        
        b_vec.Nullify();
      }
    
    if (mg != 0 and !var->include_flow_term)
      {
        // nSrc premieres colonnes : mode positif mg (deja rempli)
        // nSrc dernieres colonnes : mode negatif -mg
        // pour Dirac, on a un simple conjugue
        for (int i = 0; i < nSrc; i++)
          {
            int offset = i*var->GetNbDof();
            for (int j = 0; j < sol.GetM(); j++)
	      sol(j, nSrc+i) = conjugate(b_source(offset + j));
          }
      }
    
    b_source.Nullify();
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(comm_freq, &rank_proc);
    if (rank_proc == 0)
      cout << "Size taken by b_source " << GetHumanReadableMemory(sol.GetMemorySize()) << endl;
#endif
    
  }
  

  void CrossCovarianceParameter
  ::GiveSolutionFrequency(int k, int kg, int m, int mg, const VarComputationProblem& var_exp,
                          Matrix<Complex_wp, General, ColMajor>& sol_exp)
  {
    const VarOutputProblem_Dim<Dimension2>& var_output = dynamic_cast<const VarOutputProblem_Dim<Dimension2>&>(var_exp);
    
    cout << "J'ecris la frequence " << kg << " mode " << GetRealM(mg) << endl;
    
    // cas basique, on ecrit la solution dans un fichier
    /* var->SetCurrentModeNumber(mg);
    VectComplex_wp x_sol(sol_exp.GetM());
    GetCol(sol_exp, 0, x_sol); DISP(x_sol(1000));
    var_output.WriteOutputFile(x_sol, "SolOm" + to_str(kg) + "P" + to_str(mg) + "_G", 1);

    if (mg != 0)
      {
        var->SetCurrentModeNumber(-mg);
        GetCol(sol_exp, 1, x_sol); DISP(x_sol(1000));
        var_output.WriteOutputFile(x_sol, "SolOm" + to_str(kg) + "M" + to_str(mg) + "_G", 1);    

        var->SetCurrentModeNumber(mg);
      }
    */
    
    // on ecrit toutes les solutions d'un coup si c'est demande
    if (display_solution)
      {
        int nb_comp = 1;
        if (mg != 0 and !var->include_flow_term)
          nb_comp = 2;
        
        var->SetCurrentModeNumber(GetRealM(mg));
        VectComplex_wp x_sol(sol_exp.GetDataSize()/nb_comp, sol_exp.GetData());
	if (mg == 0)
	  x_solFull = x_sol;
	else
	  x_solFull += x_sol;

	// Write only the sum over m.
	// !!! TO FIX if phi is different from 0.
        if (mg != 0 and !var->include_flow_term)
          {
            var->SetCurrentModeNumber(-mg);
	    x_sol.Nullify();
            x_sol.SetData(sol_exp.GetDataSize()/nb_comp, &sol_exp(0, this->nSrc));
	    x_solFull += x_sol;
	  }
	IVect num_unknown(sol_exp.GetN()/nb_comp); num_unknown.Zero();
	var_output.WriteOutputFile(x_solFull, "Sol_f" + to_str(kg) + "_G",
                                       sol_exp.GetN()/nb_comp, num_unknown);
	
	x_sol.Nullify();
	
	
	/*
        IVect num_unknown(sol_exp.GetN()/nb_comp); num_unknown.Zero();
        var_output.WriteOutputFile(x_sol, "SolOm" + to_str(kg) + "P" + to_str(mg) + "_G",
                                   sol_exp.GetN()/nb_comp, num_unknown);
        
        x_sol.Nullify();
        if (mg != 0)
          {
            var->SetCurrentModeNumber(-mg);
            x_sol.SetData(sol_exp.GetDataSize()/nb_comp, &sol_exp(0, this->nSrc));
            var_output.WriteOutputFile(x_sol, "SolOm" + to_str(kg) + "M" + to_str(mg) + "_G",
                                       sol_exp.GetN()/nb_comp, num_unknown);
            
            x_sol.Nullify();
          }
	*/
	
	
      }
    
    // compute cross-covariance or kernel if necessary
    CrossCovarianceProblem& var_cross = dynamic_cast<CrossCovarianceProblem&>(*var_multi);
    if (var_cross.compute_covariance)
      var_cross.computeXS(*var, k, GetRealM(mg), nSrc, sol_exp);    
    if (var_cross.compute_kernel)
      var_cross.computeKernel(*var, k, GetRealM(mg), nSrc, sol_exp);
  }


  /**************************
   * CrossCovarianceProblem *
   **************************/

  CrossCovarianceProblem::CrossCovarianceProblem(CrossCovarianceParameter& param)
    : pt_source(param.pt_source), all_param(param)
  {
    compute_covariance = false;
    compute_kernel = false;
    num_grid = 0;
    nb_modes = 1;
    filename = "mat_covariance.dat";
    ref_surface_observation = 1;
    rmin = 0; rmax = 2; 
    RSUN = 696.e6;
    refinement_vertex_source = true;
  }
  
  
  void CrossCovarianceProblem::SetInputData(const string& keyword, const Vector<string>& params)
  {
    MultiFrequencyProblem::SetInputData(keyword, params);
    
    if (keyword == "NumberModes")
      {
	string rangeMode = params(0);
        nb_modes = to_num<int>(params(1));
	if (rangeMode == "ALL")
	  this->nb_glob_parameters = 2*nb_modes+1;
	else if (rangeMode == "POSITIVE")
	  this->nb_glob_parameters = nb_modes+1;
	else
	  {
	    cout << "The first parameter of NumberModes should be POSITIVE or ALL" << endl;
	    abort();
	  }
      }
    else if (keyword == "Sources")
      {
	int nbSources;
        if (params.GetM() < 1) 
          {
            cout << "Sources requires a parameter :" << endl;
            cout << "Sources = LIST/CIRCLE/LINE followed by options" << endl;
            cout << "Current entry is : " << endl << params << endl;
            abort();
          }
	all_param.typeSource = params(0);
	if (all_param.typeSource == "LIST")
	  {
	    string coordSystem = params(1);
	    if (coordSystem != "CARTESIAN" and coordSystem != "SPHERICAL")
	      {
		cout << "Only CARTESIAN and SPHERICAL coordinate systems are implemented" << endl;
		abort();
	      }
	    nbSources = (params.GetM()-2)/3; // each source has 3 coordinates
	    pt_source.Reallocate(nbSources);
	    for (int i =0; i < nbSources; i++)
	      {
		for (int j = 0; j < 3; j++)
		  pt_source(i)(j) = to_num<Real_wp>(params(i*3+j+2));
	      }
	    if (coordSystem ==  "SPHERICAL")
	      {
		for (int i =0; i < nbSources; i++)
		  SphericalToCartesian(pt_source(i)(0), pt_source(i)(1)*pi_wp/180., pt_source(i)(2)*pi_wp/180.,pt_source(i)(0), pt_source(i)(1), pt_source(i)(2));
	      }

	  }
	else if (all_param.typeSource == "CIRCLE")
	  {
	    Real_wp R = to_num<Real_wp>(params(1));
	    Real_wp thMin = to_num<Real_wp>(params(2))*pi_wp/180.;
	    Real_wp thMax = to_num<Real_wp>(params(3))*pi_wp/180.;
	    nbSources = to_num<int>(params(4));
	    Real_wp phi;
	    if (params.GetM() > 5)
	      phi = to_num<Real_wp>(params(5));
	    else
	      phi = 0.;
	    //cout << nbSources << " sources located at R = " << R << ", phi = " << phi << " and theta between " << thMin << " and " << thMax << endl;
	    pt_source.Reallocate(nbSources);
	    for (int i = 0; i < nbSources; i++)
	      {
		SphericalToCartesian(R,thMin+i*(thMax-thMin)/nbSources,phi,pt_source(i)(0), pt_source(i)(1), pt_source(i)(2));
	      }
	  }
	else if (all_param.typeSource == "LINE")
	  {
	    Real_wp xMin = to_num<Real_wp>(params(1));
	    Real_wp xMax = to_num<Real_wp>(params(2));
	    Real_wp yMin = to_num<Real_wp>(params(3));
	    Real_wp yMax = to_num<Real_wp>(params(4));
	    Real_wp zMin = to_num<Real_wp>(params(5));
	    Real_wp zMax = to_num<Real_wp>(params(6));
	    nbSources = to_num<int>(params(7));
	    pt_source.Reallocate(nbSources);
	    for (int i = 0; i < nbSources; i++)
	      {
		pt_source(i)(0) = xMin+i*(xMax-xMin)/nbSources;
		pt_source(i)(1) = yMin+i*(yMax-yMin)/nbSources;
		pt_source(i)(2) = zMin+i*(zMax-zMin)/nbSources;
	      }
	  }
	else if (all_param.typeSource == "RANDOM")
	  {
	    nbSources = to_num<int>(params(1));
	    all_param.lMax = to_num<int>(params(2));
	  }
	else
	  {
	    cout << "This type of source is not implemented" << endl;
	    cout << "Only LINE, LIST and CIRCLE are implemented" << endl;
	    abort();
	  }
	all_param.nSrc = nbSources;
      }
    else if (keyword == "Filter")
      {
	all_param.CreateFilter(params);
      }
    else if (keyword == "Damping")
      {
	all_param.InitDamping(params);
      }
    else if (keyword == "SourcePowerFreq")
      {
	all_param.InitSourcePower(params);
      }
    else if (keyword == "SourcePowerDepth")
      {
	all_param.source_power_depth.ReadText(params(0));
      }
    else if (keyword == "DepthSources")
      {
	all_param.depth_sources.ReadText(params(0));
      }
    else if (keyword == "DeltaDepthSources")
      {
	all_param.delta_depth_sources.ReadText(params(0));
      }
    else if (keyword == "DisplaySolution")
      {
        if (params.GetM() < 1) 
          {
            cout << "DisplaySolution requires a parameter :" << endl;
            cout << "DisplaySolution = YES" << endl;
            cout << "Current entry is : " << endl << params << endl;
            abort();
          }

        if (params(0) == "YES")
          all_param.display_solution = true;
        else
          all_param.display_solution = false;
      }
    else if (keyword == "ReferenceObserv")
      {
        if (params.GetM() < 1)
          {
            cout << "ReferenceObserv requires an integer :" << endl;
            cout << "ReferenceObserv = 0" << endl;
            cout << "Current entry is : " << endl << params << endl;
            abort();
          }

        ref_surface_observation = to_num<int>(params(0));
      }
    else if (keyword == "RefinementVertexSource")
      {
	if (params(0) == "YES")
	  {
	    if (params.GetM() < 3)
	      {
		cout << "RefinementVertexSource Requires level and ratio :" << endl;
		cout << "RefinementVertexSource = YES 2 2.0" << endl;
		cout << "Current entry is : " << endl << params << endl;
		abort();
	      }
	    
	    refinement_vertex_source = true;
	    lvl_refinement = to_num<int>(params(1));
	    ratio_refinement = to_num<Real_wp>(params(2));
	  }
	else
	  refinement_vertex_source = false;
      }
    else if (!keyword.compare("KernelPairs"))
      {
        if (params.GetM() < 2)
          {
            cout << "KernelPairs Requires pairs of integers :" << endl;
            cout << "KernelPairs = 0 1" << endl;
            cout << "Current entry is : " << endl << params << endl;
            abort();
          }
        
        if (params(0) == "SingleRow")
          {
            int num_row = to_num<int>(params(1));
	    row_mat.Reallocate(pt_source.GetM());
	    col_mat.Reallocate(pt_source.GetM());
            for (int i = 0; i < pt_source.GetM(); i++)
              {
		row_mat(i) = num_row;
		col_mat(i) = i;
              }
          }
        else
          {
            int nb = 0;  
            int N = params.GetM()/2;
	    row_mat.Reallocate(N);
	    col_mat.Reallocate(N);	    
            for (int i = nb; i < params.GetM(); i += 2)
              {
                row_mat(i/2) = to_num<int>(params(i));
                col_mat(i/2) = to_num<int>(params(i+1));
              }
          }
      }
    else if (!keyword.compare("CrossCovariance"))
      {
        
        if (params.GetM() < 1)
          {
            cout << "CrossCovariance requires a filename" << endl;
            cout << "CrossCovariance = filename" << endl;
            cout << "Current entry is : " << endl << params << endl;
            abort();
          }
        
        filename = params(0);
	compute_covariance = true;
	if (params.GetM() > 1)
	  {
	    if (params(1) == "CUT")
	      {
		if (params.GetM() > 2)
		  rmin = to_num<Real_wp>(params(2));
		if (params.GetM() > 3)
		  rmax = to_num<Real_wp>(params(3));
	      }
	    else
	      {
		cout << "This type of method to compute the cross-covariance is not known" << endl;
		cout << "Only CUT rmin rmax is implemented." << endl;
	      }
	  }
      }

    else if (!keyword.compare("Kernel"))
      {
        if (params.GetM() < 1)
          {
            cout << "CrossCovariance requires a filename" << endl;
            cout << "CrossCovariance = filename" << endl;
            cout << "Current entry is : " << endl << params << endl;
            abort();
          }
	filename = params(0);
	compute_kernel = true;
	typeOfKernels.Reallocate(params.GetM()-1);
	for (int i=1; i<params.GetM(); i++)
	  typeOfKernels(i-1) = params(i);
      }

    else if (keyword == "FrequencyPonderation")
      {
        if (params.GetM() < 1)
          {
            cout << "FrequencyPonderation requires a string :" << endl;
            cout << "FrequencyPonderation = coef.dat" << endl;
            cout << "Current entry is : " << endl << params << endl;
            abort();
          }

        coef_freq.ReadText(params(0));
      }
    else if (keyword == "DirectoryOutput")
      {
        if (params.GetM() < 1)
          {
            cout << "DirectoryOutput requires a string :" << endl;
            cout << "DirectoryOutput = /here/" << endl;
            cout << "Current entry is : " << endl << params << endl;
            abort();
          }

        directoryOutput = params(0);
      }

  }  


  void CrossCovarianceProblem::SetProblem(DistributedProblem<Dimension2>& var)
  {
    if (refinement_vertex_source)
      {
	VectR2 pt_source2D(pt_source.GetM());
	for (int i = 0; i < pt_source.GetM(); i++)
	  pt_source2D(i).Init(sqrt(square(pt_source(i)(0)) + square(pt_source(i)(1))), pt_source(i)(2));

	// it is better to not add vertices on the mesh, because triangles are created
	// mass lumping is no longer true => computeXs false
	// var.mesh.SetVerticesToBeAdded(pt_source2D);

	var.mesh.SetVerticesToBeRefined(pt_source2D, lvl_refinement, ratio_refinement);
      }
  }
  
  
  void CrossCovarianceProblem::InitCovariance(VarOutputProblem_Dim<Dimension2>& var_output)
  {
    // on alloue mat_covariance
    if (compute_covariance)
      {
        mat_covariance.Reallocate(this->omega.GetM(), row_mat.GetM());
        mat_covariance.Zero();
      }
    else
      {
        int nb_pts_grid = var_output.var_grid(num_grid).GetNbGlobalPoints();
	mat_kernel.Reallocate(nb_pts_grid, row_mat.GetM()*typeOfKernels.GetM());
	mat_kernel.Zero();
	  
      }
  }

  // Fills the coefficients rho, c and gamma for the element defined by ID_elt and the point ID_pt.
  template<class TypeEquation> 
  void CrossCovarianceProblem::getPhysicalCoefficients(EllipticProblem<TypeEquation>& var, int ID_elt, int ID_pt, int ref, Complex_wp& rho, Complex_wp& c, Complex_wp& gamma)
  {
    // Montjoie coefficients
    Complex_wp rho_MJ    = var.ref_rho  (ref).GetCoefficient(var, ID_elt, ID_pt);
    Complex_wp mu_MJ     = var.ref_mu   (ref).GetCoefficient(var, ID_elt, ID_pt)(0,0);
    Complex_wp sigma_MJ  = var.ref_sigma(ref).GetCoefficient(var, ID_elt, ID_pt);
            
    // Physical quantities
    gamma = 0.5*sigma_MJ/rho_MJ/RSUN;
    rho   = 1.e0/mu_MJ;
    c     = sqrt(mu_MJ/rho_MJ);
  }


  template<class TypeEquation>
  void CrossCovarianceProblem::ComputeMassMatrix(EllipticProblem<TypeEquation>& var)
  {
    mat_mass.Reallocate(var.GetNbDof());
    mat_mass.Zero();

    // boucle elements (pour la partie volumique)
    Real_wp jacob;
    Complex_wp omega = var.GetOmega()/RSUN;
    Complex_wp rho, c, gamma, coeff;
    Real_wp r;
    int nb_elt = var.GetNbElt();
  
    for (int ID_elt = 0; ID_elt < nb_elt; ID_elt++)
      {	
        const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(ID_elt);
        int N   = Fb.GetNbPointsQuadratureInside();
        
        // The following matrices will be computed:
        // \int CoefMass  \phi_i\phi_j
        int  ref    = var.mesh.Element(ID_elt).GetReference();
        bool affine = var.mesh.IsElementAffine(ID_elt);
        
        // Loop on all volumetric quadrature points of the element (\xi_j)
        for (int ID_pt = 0; ID_pt < N; ID_pt++)
          {
            jacob = var.GetWeightedJacobian(ID_elt, ID_pt, affine, Fb.GetGeometricElement());
            r       = realpart(var.Glob_rtilde(ID_elt)(ID_pt));
 
	    getPhysicalCoefficients(var, ID_elt, ID_pt, ref, rho, c, gamma);
            coeff = gamma/(rho*square(c)); // coefficient in front of G_1* G_2 in the integral for the cross-covariance
            
            int cpt = var.mesh_num.Element(ID_elt).GetNumberDof(ID_pt);

            if (r > rmin and r < rmax)
              mat_mass(cpt) += realpart(coeff * jacob * r * RSUN)*2.*pi_wp; // 2pi for the phi integral
          }

        // loop on edges of the element
        for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
          {
            int ne = var.mesh.Element(ID_elt).numBoundary(num_loc);
            int refB = var.mesh.Boundary(ne).GetReference();
            if ((refB == ref_surface_observation) && var.IsNewFace(ID_elt)(num_loc))
              {
                int Nquad = Fb.GetNbQuadBoundary(num_loc);
                for (int j = 0; j < Nquad; j++)
                  {
                    jacob = var.GetSurfaceWeightedJacobian(ID_elt, num_loc, ne, j,
                                                           Fb.GetGeometricElement());
                    
                    int ID_pt = Fb.GetQuadNumber(num_loc, j);
                    int cpt = var.mesh_num.Element(ID_elt).GetNumberDof(ID_pt);
		    getPhysicalCoefficients(var, ID_elt, ID_pt, ref, rho, c, gamma);
		    Real_wp alpha = 6.6334629324389107e3/RSUN;
		    Real_wp BC_part = real(1-square(c*alpha/(2.*omega)));
		    //cout << "c = " << c << ", alpha = " << alpha << "omega = " << omega << endl;
		    if (BC_part > 0)
		      {
			BC_part = sqrt(BC_part);
			coeff = 1./(2.*rho*c*BC_part);
			r     = realpart(var.Glob_rtilde(ID_elt)(ID_pt));
                    
			mat_mass(cpt) += realpart(coeff*jacob)*2.*pi_wp;
		      }
                  }
              }
          }
      }
  }
                

  template<class TypeEquation>
  void CrossCovarianceProblem::computeXS(EllipticProblem<TypeEquation>& var, int nfreq,
					 int mg, int ns,
                                         Matrix<Complex_wp, General, ColMajor>& all_sol)
  {    
    for (int m = 0; m < row_mat.GetM(); m++)
      {
	int i = row_mat(m);
	int j = col_mat(m);
	Complex_wp vloc = 0;
	if (mg == 0 or var.include_flow_term)
	  for (int k = 0; k < var.GetNbDof(); k++)
	    vloc += mat_mass(k)*conjugate(all_sol(k, i)) * all_sol(k, j);
	else
	  for (int k = 0; k < var.GetNbDof(); k++)
	    vloc += mat_mass(k)*(conjugate(all_sol(k, i))*all_sol(k, j)
				 + conjugate(all_sol(k, i+ns))*all_sol(k, j+ns) );
	
	mat_covariance(nfreq, m) += vloc;
      }
  }
  
  
  template<class TypeEquation>
  void CrossCovarianceProblem
  ::computeKernel(EllipticProblem<TypeEquation>& var, int k, int mg, int ns,
                  Matrix<Complex_wp, General, ColMajor>& all_sol)
  {
    Real_wp RSUN = 696.e6;
    // 1ere etape : on interpole sur la grille de sortie
    // mettre true pour avoir le gradient
    VectComplex_wp U0; Vector<VectComplex_wp> trace_u, trace_grad_u; 
    U0.SetData(all_sol.GetDataSize(), all_sol.GetData());
    var.ComputeInterpolationU_GradU(U0, num_grid, true,
                                    trace_u, trace_grad_u, all_sol.GetN());
    
    U0.Nullify();

    Vector<VectComplex_wp> trace_rho, trace_grad_rho, trace_c, trace_grad_c, trace_gamma, trace_grad_gamma; 
    VectComplex_wp c0, rho0, gamma0;
    c0.Reallocate(all_sol.GetM());
    rho0.Reallocate(all_sol.GetM()); 
    gamma0.Reallocate(all_sol.GetM());
    int nb_elt = var.GetNbElt();
    for (int ID_elt = 0; ID_elt < nb_elt; ID_elt++)
      {	
        const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(ID_elt);
        int N   = Fb.GetNbPointsQuadratureInside();
        int  ref    = var.mesh.Element(ID_elt).GetReference();
	for (int ID_pt = 0; ID_pt < N; ID_pt++)
	  {
            
	    // Montjoie coefficients
	    Complex_wp rho_MJ    = var.ref_rho  (ref).GetCoefficient(var, ID_elt, ID_pt);
	    Complex_wp mu_MJ     = var.ref_mu   (ref).GetCoefficient(var, ID_elt, ID_pt)(0,0);
	    Complex_wp sigma_MJ  = var.ref_sigma(ref).GetCoefficient(var, ID_elt, ID_pt);
            
	    // Physical quantities            
	    int cpt = var.mesh_num.Element(ID_elt).GetNumberDof(ID_pt);
	    rho0(cpt) = 1.e0/mu_MJ;
	    c0(cpt) = sqrt(mu_MJ/rho_MJ);
	    gamma0(cpt) = 0.5*sigma_MJ/rho_MJ/RSUN;
	  }
      }
    DISP("U0 interpolated");
    //U0.SetData(rho0.GetDataSize(), rho0.GetData());
    var.ComputeInterpolationU_GradU(rho0, num_grid, true,
                                    trace_rho, trace_grad_rho, 1);
    DISP("rho0 interpolated");
    var.ComputeInterpolationU_GradU(c0, num_grid, true,
                                    trace_c, trace_grad_c, 1);
    DISP("c0 interpolated");
    var.ComputeInterpolationU_GradU(gamma0, num_grid, false,
                                    trace_gamma, trace_grad_gamma, 1);
    DISP("gamma0 interpolated");
    
    rho0.Nullify(); c0.Nullify(); gamma0.Nullify();
    
    // 2eme etape : on ajoute conjugate(Gi) Gj sur tous les points de la grille
    int nb_pts_grid = trace_u(0).GetM();
    Complex_wp rho,c,gamma,drho,dc,k2;
    Complex_wp               rho_MJ ,mu_MJ ,sigma_MJ; 
    TinyVector<Complex_wp,2> drho_MJ,dmu_MJ;
    Complex_wp a,b,F,drF;
    Real_wp r,theta,phi;
    Complex_wp omega = var.GetOmega()/RSUN;
    Complex_wp I_wp  = Complex_wp(0.e0,1.e0); 
    int cpt = 0;
    const IVect& list_points = var.var_grid(0).GetPointNumber();
    DISP(list_points.GetM());
    DISP(nb_pts_grid);
    for (int ID_pair = 0; ID_pair < row_mat.GetM(); ID_pair++)
      {
	int i = row_mat(ID_pair);
	int j = col_mat(ID_pair);
	for (int ID_type = 0; ID_type < typeOfKernels.GetM(); ID_type++)
	  {
	    cpt = 0;
	    Complex_wp vloc = 0;
	    for (int ID_pt = 0; ID_pt < list_points.GetM(); ID_pt++)
	      {
		int ip = list_points(ID_pt);
		if (ip != ID_pt)
		  {
		    cout << ip << ", " << ID_pt << endl;
		  }
		//if (Norm2(var.all_points_display.GetGlobalCoordinate(ip)) > 1)
		//  cout << "ip = " << ip << "pt(ip) = " << var.all_points_display.GetGlobalCoordinate(ip) << ", ID_pt = " << ID_pt << ", pt(ID_pt) = " << var.all_points_display.GetGlobalCoordinate(ID_pt) << endl;
		if (ip < 0)
		  {
		    DISP(var.all_points_display.GetGlobalCoordinate(ip));
		    cpt++;
		    continue; // remove points outside of the domain
		  }

 
		// Get physical indices at output points
		typename TypeEquation::Dimension::R_N point = var.all_points_display.GetGlobalCoordinate(ip);
		cout << ID_pt << ", " << point << endl;
		//Real_wp rp,thetap,phip;
		CartesianToSpherical(point(0),0.e0,point(1),r,theta,phi);

		/*var.ref_mu   (1).GetVaryingIndex(0,0).index_radial.ComputeIndexAtPoint(point,   mu_MJ, dmu_MJ);
		var.ref_rho  (1).GetVaryingIndex()   .index_radial.ComputeIndexAtPoint(point,  rho_MJ,drho_MJ);
		var.ref_sigma(1).GetVaryingIndex()   .index_radial.ComputeIndexAtPoint(point,sigma_MJ);

		dmu_MJ(0)  = sin(theta)*dmu_MJ(0) + cos(theta)*dmu_MJ(1);
		dmu_MJ(1)  = cos(theta)*dmu_MJ(0) - sin(theta)*dmu_MJ(1);
		drho_MJ(0) = sin(theta)*drho_MJ(0) + cos(theta)*drho_MJ(1);
		drho_MJ(1) = cos(theta)*drho_MJ(0) - sin(theta)*drho_MJ(1);

		// mu_MJ    is 1/rho
		// rho_MJ   is 1/rho*c*c
		// sigma_MJ is 2*gamma/rho*c*c
		gamma = 0.5*sigma_MJ/rho_MJ/RSUN;
		rho   = 1.e0/mu_MJ;
		drho  = -dmu_MJ(0)/square(mu_MJ);
		c     = sqrt(mu_MJ/rho_MJ);
		dc    = -0.5*(drho_MJ(0)/(square(rho_MJ)*rho*c) + drho*c/rho);
		k2    = -(square(omega)+2*I_wp*gamma*omega)/square(c);
		//scale the derivatives with respect to RSUN
		drho = drho/RSUN;
		dc   = dc  /RSUN; */
		// Note we scale by RSUN for the greens functions and RSUN squared for the derivatives
		// This is because of (1/R^3 from the Dirac, equation already multiplied by R^2

		rho = trace_rho(0)(ip);
		c = trace_c(0)(ip);
		gamma = trace_gamma(0)(ip);
		dc = (sin(theta)*trace_grad_c(0)(ip)+cos(theta)*trace_grad_c(0)(ip+nb_pts_grid))/RSUN;
		drho = (sin(theta)*trace_grad_rho(0)(ip)+cos(theta)*trace_grad_rho(0)(ip+nb_pts_grid))/RSUN;
		k2    = -(square(omega)+2*I_wp*gamma*omega)/square(c);

		// Coefficient to convert Green's function into covariance
		// C = (G - G*) / coeff
		Complex_wp coeffCov = (4.*I_wp*omega);
		if (omega == Complex_wp(0.,0.))
		  coeffCov = (4.*I_wp);
		

		// Convert the Cartesian gradient to spherical
		//Real_wp x(var.var_grid(0).GlobalPoints2D(ID_pt)(0));
		//Real_wp z(var.var_grid(0).GlobalPoints2D(ID_pt)(1));
		//CartesianToSpherical(x,0.e0,z,r,theta,phi);

		for (int signM = 0; signM < 2; signM++)
		  {
		    if (signM > 0 and mg == 0)
		      continue;

		    Complex_wp Ui = trace_u(i)(ip)/RSUN;
		    Complex_wp Uj = trace_u(j)(ip)/RSUN;
		    Complex_wp drUi = sin(theta)*trace_grad_u(i+signM*ns)(ip)+ cos(theta)*trace_grad_u(i+signM*ns)(ID_pt+nb_pts_grid)/square(RSUN);
		    Complex_wp dtUi = cos(theta)*trace_grad_u(i+signM*ns)(ip) - sin(theta)*trace_grad_u(i+signM*ns)(ID_pt+nb_pts_grid)/square(RSUN);
		    Complex_wp dpUi = trace_u(i+signM*ns)(ip)* I_wp*mg/(r*sin(theta))/square(RSUN);
		    Complex_wp drUj = sin(theta)*trace_grad_u(j+signM*ns)(ip)+ cos(theta)*trace_grad_u(j+signM*ns)(ID_pt+nb_pts_grid)/square(RSUN);
		    Complex_wp dtUj = cos(theta)*trace_grad_u(j+signM*ns)(ip) - sin(theta)*trace_grad_u(j+signM*ns)(ID_pt+nb_pts_grid)/square(RSUN);
		    Complex_wp dpUj = trace_u(j+signM*ns)(ip)* I_wp*mg/(r*sin(theta))/square(RSUN);

		    if (typeOfKernels(ID_type) == "UR")
		      {
			a = 2*I_wp*omega/c;
			F = 1/(rho*c);
			drF = drho*c + dc*rho;
			drF = -drF * square(F);		    
			vloc = a*Uj*(F*(drUi - conj(drUi))/coeffCov+drF*(Ui - conj(Ui))/coeffCov) + conj(a*Ui*(F*(drUj - conj(drUj))/coeffCov+drF*(Uj - conj(Uj))/coeffCov));
		      }
		    else if (typeOfKernels(ID_type) == "UTHETA")
		      {
			a = 2*I_wp*omega/c;
			F = 1/(rho*c);	    
			vloc = a*Uj*F*(dtUi - conj(dtUi))/coeffCov + conj(a*Ui*F*(dtUj - conj(dtUj))/coeffCov);
		      }
		    else if (typeOfKernels(ID_type) == "UPHI")
		      {
			a = 2*I_wp*omega/c;
			F = 1/(rho*c);	    
			vloc = a*Uj*F*(dpUi - conj(dpUi))/coeffCov + conj(a*Ui*F*(dpUj - conj(dpUj))/coeffCov);
		      }
		    else if (typeOfKernels(ID_type) == "C")
		      {
			a = 2*k2/(rho*c);	    
			vloc = a*Uj*(Ui - conj(Ui))/coeffCov + conj(a*Ui*(Uj - conj(Uj))/coeffCov);
		      }
		    else if (typeOfKernels(ID_type) == "RHO")
		      {
			a = 3.*k2/(rho*rho);
			vloc = a*Uj*(Ui - conj(Ui))/coeffCov + conj(a*Ui*(Uj - conj(Uj))/coeffCov);
			a = 1./(rho*rho);
			vloc += a*(drUj*(drUi - conj(drUi))+dtUj*(dtUi - conj(dtUi))+dpUj*(dpUi - conj(dpUi)))/coeffCov + conj(a*(drUi*(drUj - conj(drUj))+dtUi*(dtUj - conj(dtUj))+dpUi*(dpUj - conj(dpUj))) /coeffCov);
		      }
		    else if (typeOfKernels(ID_type) == "GAMMA")
		      {
			a = 2*I_wp*omega/(rho*square(c));
	                b = 1/(rho*square(c));
			vloc = a*Uj*(Ui - conj(Ui))/coeffCov + conj(a*Ui*(Uj - conj(Uj))/coeffCov) + b*conj(Ui)*Uj;
		      }
		    mat_kernel(ID_pt, ID_pair*typeOfKernels.GetM()+ID_type) += vloc;
		  } // +/- m
	      } // ID_pt
          } // ID_type
      } // ID_pair
    DISP(cpt);
  }
  
  
  void CrossCovarianceProblem::WriteCovariance()
  {
    
#ifdef SELDON_WITH_MPI
    // summing the contribution of modes (with comm_param)
    int nb_proc; MPI_Comm_size(comm_param, &nb_proc);
    Vector<int64_t> xtmp;
    if (nb_proc > 1)
      {
	Matrix<Complex_wp> mat_send(mat_covariance);
	MpiReduce(comm_param, mat_send.GetData(), xtmp,
		  mat_covariance.GetData(), mat_send.GetDataSize(), MPI_SUM, 0);
      }
    
    // gathering the contributions for different frequencies (with comm_first_param)
    if (this->offset_param == 0)
      {
	MPI_Comm_size(comm_first_param, &nb_proc);
	int rank_proc; MPI_Comm_rank(comm_first_param, &rank_proc);
	Vector<int> nb_freq_per_proc;
	if (rank_proc == 0)
	  nb_freq_per_proc.Reallocate(nb_proc);
	
	int nb_freq = this->omega.GetM();
	MPI_Gather(&nb_freq, 1, MPI_INTEGER, nb_freq_per_proc.GetData(), 1, MPI_INTEGER, 0, comm_first_param);
	MPI_Status status;
	if (rank_proc == 0)
	  {
	    int nb_all = 0;
	    for (int i = 0; i < nb_proc; i++)
	      nb_all += nb_freq_per_proc(i);

	    mat_covariance.Resize(nb_all, mat_covariance.GetN());
	    int offset = this->omega.GetM();
	    Matrix<Complex_wp> mat_tmp;
	    for (int i = 1; i < nb_proc; i++)
	      {
		mat_tmp.Reallocate(nb_freq_per_proc(i), mat_covariance.GetN());
		MpiRecv(comm_first_param, mat_tmp.GetData(), xtmp,
			mat_tmp.GetDataSize(), i, 19, status);

		for (int j = 0; j < nb_freq_per_proc(i); j++)
		  for (int k = 0; k < mat_covariance.GetN(); k++)
		    mat_covariance(offset + j, k) = mat_tmp(j, k);
		
		offset += nb_freq_per_proc(i);
	      }
	  }
	else
	  {
	    MpiSsend(comm_first_param, mat_covariance.GetData(), xtmp,
		     mat_covariance.GetDataSize(), 0, 19);

	    // only the proc 0 will write mat_covariance
	    return;
	  }
      }
    else
      return;
    
#endif
    
    mat_covariance.Write(directoryOutput + filename);
  }

  void CrossCovarianceProblem::WriteKernel()
  {
#ifdef SELDON_WITH_MPI
    // summing the contribution of modes and frequency (with comm_freq)
    int nb_proc; MPI_Comm_size(comm_freq, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm_freq, &rank_proc);
    Vector<int64_t> xtmp;
    if (nb_proc > 1)
      {
	Matrix<Complex_wp> mat_send(mat_kernel);
	MpiReduce(comm_freq, mat_send.GetData(), xtmp,
		  mat_kernel.GetData(), mat_send.GetDataSize(), MPI_SUM, 0);
      }
    
    // seul le maitre ecrit
    if (rank_proc > 0)
      return;
    
#endif
    
    mat_kernel.Write(directoryOutput + filename);
  }

}

#define MONTJOIE_FILE_CROSS_COVARIANCE_CXX
#endif
