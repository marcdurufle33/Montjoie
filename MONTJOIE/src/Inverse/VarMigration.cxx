#ifndef MONTJOIE_FILE_VAR_MIGRATION_CXX

#include "VarMigration.hxx"

namespace Montjoie
{
  
  //! default constructor
  template<class Dimension>
  VarMigration_Base<Dimension>::VarMigration_Base()
  {
    axisymmetric_problem = false;
    nb_points_global_measure = 0;
  }


  //! returns the number of quadrature points for the current processor
  template<class Dimension>
  int VarMigration_Base<Dimension>::GetNbPointsQuadrature() const
  {
    if (nb_points_global_measure > 0)
      return weights_measure.GetM();
    
    return mesh_interp.GetNbPointsQuadrature();
  }
  
  
  //! fills the weights (w_i ds_i for each quadrature point)
  template<class Dimension>
  void VarMigration_Base<Dimension>
  ::GetWeightsQuadrature(VectReal_wp& weights) const
  {
    if (nb_points_global_measure > 0)
      {
	weights = weights_measure;
	return;
      }
    
    weights.Reallocate(GetNbPointsQuadrature());
    int num = 0;
    for (int i = 0; i < mesh_interp.GetNbBoundary(); i++)
      for (int j = 0; j < mesh_interp.GetNbPointsQuadrature(i); j++)
	weights(num++) = mesh_interp.WeightsQuadrature(i, j);
  }
  
  
  //! a line of the data file is read
  template<class Dimension>
  void VarMigration_Base<Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (description_field == "DataFileExperiment")
      {
	input_file_exp = parameters(0);
	ref_measure_exp.Reallocate(parameters.GetM()-1);	
	for (int j = 1; j < parameters.GetM(); j++)
	  ref_measure_exp(j-1) = to_num<int>(parameters(j));
      }
    else if (description_field == "DataFileSimulation")
      {
	input_file_simu = parameters(0);
	ref_measure_simu.Reallocate(parameters.GetM()-1);	
	for (int j = 1; j < parameters.GetM(); j++)
	  ref_measure_simu(j-1) = to_num<int>(parameters(j));
      }
    else if (description_field == "PointsMeasure")
      {
	points_measure.ReadText(parameters(0));
	if (parameters.GetM() > 1)
	  weights_measure.ReadText(parameters(1));
	else
	  {
	    weights_measure.Reallocate(points_measure.GetM());
	    weights_measure.Fill(Real_wp(1));
	  }
        
        nb_points_global_measure = weights_measure.GetM();
      }
    else if (description_field == "TimeImpulsion")
      {
	nb_points_time = to_num<int>(parameters(0));
	if (parameters.GetM() <= 4)
	  {
	    cout << "In SetInputData of VarMigration" << endl;
	    cout << "TimeImpulsion needs at least six parameters, for instance :" << endl;
	    cout << "TimeImpulsion = N GaussianSinus tau Tmax t_sinus" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	  
	if (parameters(1) == "GaussianSinus")
	  {
	    Real_wp fwhm_impulse = to_num<Real_wp>(parameters(2));
	    //Relation between fwhm and tau
	    tau_impulse = fwhm_impulse /(Real_wp(2) * sqrt(Real_wp(2)*log(Real_wp(2)))); 	    
	    Tmax = to_num<Real_wp>(parameters(3));
	    Tsinus = to_num<Real_wp>(parameters(4));
	  }
      }
    else if (description_field == "TypeSourceLine")
      {
	int nb_source = to_num<int>(parameters(0));
	VectString param;
	int nb = 1, ref = ref_measure_simu(0);
	if (parameters(1) == "REF")
	  {
	    ref = to_num<int>(parameters(2));
	    nb = 3;
	  }
	
	int ns = param_source.GetM();
	param_source.Resize(nb_source + ns);
	int dim_N = Dimension::dim_N;
	if (axisymmetric_problem)
	  dim_N = 3;
	
	if (parameters(nb) == "Gaussian")
	  {
	    nb++;
	    // case of gaussian sources along a line
	    // the user puts the two extremities
	    R_N ptSourceA, ptSourceB, pt;
	    for (int k = 0; k < ptSourceA.GetM(); k++)
	      to_num(parameters(nb++), ptSourceA(k));
	      
	    for (int k = 0; k < ptSourceB.GetM(); k++)
	      to_num(parameters(nb++), ptSourceB(k));
	      
	    // other parameters are not changed
	    param.Reallocate(parameters.GetM() - nb + 4 + dim_N);
	    param(1) = "SRC_SURFACE";
	    param(2) = to_str(ref);	    
	    param(3) = "GAUSSIAN";
	    int offset = 4 + dim_N;
	    for (int k = nb; k < parameters.GetM(); k++)
	      param(offset + k - nb) = parameters(k);
	      
	    // loop over nb_source to generate parameters for each point
	    for (int n = 0; n < nb_source; n++)
	      {
		param(0) = to_str(ns + n);
		Real_wp L = Real_wp(n) / (nb_source - 1);
		pt = (Real_wp(1) - L)*ptSourceA + L*ptSourceB;
		if (axisymmetric_problem)
		  {
		    // axisymmetric case
		    param(4) = to_str(pt(0));
		    param(5) = "0.0";
		    param(6) = to_str(pt(1));
		  }
		else
		  for (int k = 0; k < ptSourceB.GetM(); k++)
		    param(4+k) = to_str(pt(k));
		  
		param_source(ns + n) = param;
	      }
	  }
	else if (parameters(nb) == "PlaneWave")
	  {
	    nb++;
	    // case of surface sources of the type exp(i k \cdot x)
	    // k = omega (cos teta, sin teta)
	    // where teta will vary between two angles given by the user
	    Real_wp tetaIncidenceA, tetaIncidenceB;
	    tetaIncidenceA = to_num<Real_wp>(parameters(nb++));
	    tetaIncidenceB = to_num<Real_wp>(parameters(nb++));
	    param.Reallocate(parameters.GetM() - nb + 7);
	    param(1) = "SRC_SURFACE";
	    param(2) = to_str(ref);
	    param(3) = "INCIDENT_WAVE";
	    param(4) = "PLANE_WAVE";
	    param(5) = "Theta";
	    int offset = 7;
	    for (int k = nb; k < parameters.GetM(); k++)
	      param(offset + k - nb) = parameters(k);
	      
	    for (int n = 0; n < nb_source; n++)
	      {
		param(0) = to_str(ns + n);
		Real_wp L = Real_wp(n) / (nb_source - 1);
		Real_wp teta = (Real_wp(1) - L)*tetaIncidenceA + L*tetaIncidenceB;
		param(6) = to_str(teta);
		  
		param_source(ns + n) = param;
	      }
	  }
	else if (parameters(nb) == "IncidentWave")
	  {
	    nb++;
	    // incident plane wave, the user gives two incident angles
	    Real_wp tetaIncidenceA, tetaIncidenceB;
	    tetaIncidenceA = to_num<Real_wp>(parameters(nb++));
	    tetaIncidenceB = to_num<Real_wp>(parameters(nb++));
	    param.Reallocate(5 + parameters.GetM() - nb);
	    param(1) = "SRC_TOTAL_FIELD";
	    param(2) = "PLANE_WAVE";
	    param(3) = "Theta";
	      
	    for (int k = nb; k < parameters.GetM(); k++)
	      param(5 + k - nb) = parameters(k);
	      
	    for (int n = 0; n < nb_source; n++)
	      {
		param(0) = to_str(ns + n);
		Real_wp L = Real_wp(n) / (nb_source - 1);
		Real_wp teta = (Real_wp(1) - L)*tetaIncidenceA + L*tetaIncidenceB;
		param(4) = to_str(teta);
		  
		param_source(ns + n) = param;
	      }
	  }
	else 
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }	
      }
    else if (description_field == "TypeSourceCircle")
      {
	int nb_source = to_num<int>(parameters(0));
	VectString param;
	int nb = 1, ref = ref_measure_simu(0);
	if (parameters(1) == "REF")
	  {
	    ref = to_num<int>(parameters(2));
	    nb = 3;
	  }
	
	int dim_N = Dimension::dim_N;
	if (axisymmetric_problem)
	  dim_N = 3;

	int ns = param_source.GetM();
	param_source.Resize(nb_source + ns);
	if (parameters(nb) == "Gaussian")
	  {
	    nb++;
	    // case of gaussian sources along a circle
	    // the user puts the two extremities
	    R_N ptSourceA;
	    for (int k = 0; k < ptSourceA.GetM(); k++)
	      to_num(parameters(nb++), ptSourceA(k));
	    
	    R_N pt(ptSourceA);
	    Real_wp radius = to_num<Real_wp>(parameters(nb++));
	    Real_wp tetaA = to_num<Real_wp>(parameters(nb++))*pi_wp/180;
	    Real_wp tetaB = to_num<Real_wp>(parameters(nb++))*pi_wp/180;
	    
	    // other parameters are not changed
	    param.Reallocate(parameters.GetM() - nb + 4 + dim_N);
	    param(1) = "SRC_SURFACE";
	    param(2) = to_str(ref);	    
	    param(3) = "GAUSSIAN";
	    int offset = 4 + dim_N;
	    for (int k = nb; k < parameters.GetM(); k++)
	      param(offset + k - nb) = parameters(k);
	      
	    // loop over nb_source to generate parameters for each point
	    Real_wp dteta = (tetaB - tetaA) / (nb_source-1);
	    if (nb_source == 1)
	      dteta = Real_wp(0);
	    
	    for (int n = 0; n < nb_source; n++)
	      {
		param(0) = to_str(ns + n);
		Real_wp teta = tetaA + n*dteta;
		pt(0) = ptSourceA(0) + radius*cos(teta);
		pt(1) = ptSourceA(1) + radius*sin(teta);
		if (axisymmetric_problem)
		  {
		    // axisymmetric case
		    param(4) = to_str(pt(0));
		    param(5) = "0.0";
		    param(6) = to_str(pt(1));
		  }
		else
		  for (int k = 0; k < ptSourceA.GetM(); k++)
		    param(4+k) = to_str(pt(k));
		  
		param_source(ns + n) = param;
	      }
	  }
	else if (parameters(nb) == "Dirac")
	  {
	    nb++;
	    // case of Dirac sources along a circle
	    // the user puts the two extremities
	    R_N ptSourceA;
	    for (int k = 0; k < ptSourceA.GetM(); k++)
	      to_num(parameters(nb++), ptSourceA(k));
	    
	    R_N pt(ptSourceA);
	    Real_wp radius = to_num<Real_wp>(parameters(nb++));
	    Real_wp tetaA = to_num<Real_wp>(parameters(nb++))*pi_wp/180;
	    Real_wp tetaB = to_num<Real_wp>(parameters(nb++))*pi_wp/180;
	    
	    // other parameters are not changed
	    param.Reallocate(parameters.GetM() - nb + 2 + Dimension::dim_N);
	    param(1) = "SRC_DIRAC";
	    int offset = 2 + Dimension::dim_N;
	    for (int k = nb; k < parameters.GetM(); k++)
	      param(offset + k - nb) = parameters(k);
	      
	    // loop over nb_source to generate parameters for each point
	    Real_wp dteta = (tetaB - tetaA) / (nb_source-1);
	    if (nb_source == 1)
	      dteta = Real_wp(0);
	    
	    for (int n = 0; n < nb_source; n++)
	      {
		param(0) = to_str(ns + n);
		Real_wp teta = tetaA + n*dteta;
		pt(0) = ptSourceA(0) + radius*cos(teta);
		pt(1) = ptSourceA(1) + radius*sin(teta);
		for (int k = 0; k < ptSourceA.GetM(); k++)
		  param(2+k) = to_str(pt(k));
		  
		param_source(ns + n) = param;
	      }
	  }
      }
    else if (description_field == "TypeSource")
      {
	// keyword used to specify a source as other targets
	param_source.PushBack(parameters);
      }
    else if (description_field == "WriteIntermediate")
      {
	if (parameters(0) == "YES")
	  write_intermediate_output = true;
	else
	  write_intermediate_output = false;
      }
  }

  
  //! assures uniqueness of each localized point (a point belongs to an unique processor)
  template<class Dimension>
  int VarMigration_Base<Dimension>::RemoveDuplicateLocalization(VarProblem<Dimension>& var, IVect& NumElement, IVect& ProcLocal)
  {
    ProcLocal.Reallocate(NumElement.GetM());
    
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc = 1, rank_proc(0);
#endif

    ProcLocal.Fill(nb_proc);
    for (int k = 0; k < NumElement.GetM(); k++)
      if (NumElement(k) >= 0)
        ProcLocal(k) = rank_proc;
    
#ifdef SELDON_WITH_MPI
    IVect Proc_tmp(ProcLocal);
    MPI_Allreduce(Proc_tmp.GetData(), ProcLocal.GetData(), ProcLocal.GetM(), MPI_INTEGER, MPI_MIN, comm);

    // we avoid duplicate points (localized on two processors)
    for (int k = 0; k < NumElement.GetM(); k++)
      if (ProcLocal(k) != rank_proc)
        NumElement(k) = -1;
#endif

    // we count the number of points
    int nb_points = 0;
    for (int k = 0; k < NumElement.GetM(); k++)
      if (NumElement(k) >= 0)
	nb_points++;
    
    return nb_points;
  }
  
  
  // the quadrature points are constructed and localized for the experimental data
  template<class Dimension>
  void VarMigration_Base<Dimension>::ComputeQuadrature(VarProblem<Dimension>& var,
						       VarProblem<Dimension>& var_exp)
  {
    // we check that references are the same for the simulation and the experience
    bool ref_diff = false;
    if (ref_measure_simu.GetM() != ref_measure_exp.GetM())
      ref_diff = true;
    else
      for (int i = 0; i < ref_measure_simu.GetM(); i++)
	if (ref_measure_simu(i) != ref_measure_exp(i))
	  ref_diff = true;
      
    if (ref_diff)
      {
	cout << "References for simulation and experiment must be the same" << endl;
	abort();
      }
      
    // quadrature points on the surface gamma are computed
    // and gathered for all processors    
    Mesh<Dimension> mesh_subdiv;
    
    Vector<int> ref_cond(var.mesh.GetNbReferences()+1);
    ref_cond.Zero();
    for (int i = 0; i < ref_measure_simu.GetM(); i++)
      ref_cond(ref_measure_simu(i)) = 1;

    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    if (nb_points_global_measure == 0)
      {
	mesh_interp.SetGaussQuadrature(mesh_num.GetOrder());
	mesh_interp.InitProjectionSurface(var.mesh);
	mesh_interp.ComputeSurfaceMesh(ref_cond, var.mesh, mesh_subdiv, var);
      }
      
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    if (nb_points_global_measure == 0)
      mesh_interp.GatherQuadraturePoints(var.comm_group_mode);
#else
    int nb_proc = 1, rank_proc(0);
    if (nb_points_global_measure == 0)
      mesh_interp.GatherQuadraturePoints();
#endif
      
    VectR_N* Points_ptr = NULL;
    if (nb_points_global_measure == 0)
      Points_ptr = const_cast<VectR_N*>(&mesh_interp.GetQuadraturePoints());
    else
      Points_ptr = &points_measure;
    
    const VectR_N& Points = *Points_ptr;
    /*
      VectReal_wp PointsX(Points.GetM());
      for (int i = 0; i < Points.GetM(); i++)
      PointsX(i) = Points(i)(0);
	
      PointsX.WriteText("PointsX.dat"); */
      
    // on each processor of experimental simulation we try to localize
    // these quadrature points
    // NumElement(i) : element number where the point i is
    // NumLocal(i) : local boundary on the element
    // LocalPos(i) : local coordinates on the boundary
    IVect NumElement(Points.GetM()), NumLocal(Points.GetM());
    typename Dimension::VectR_N CoorLocal;
    typename Dimension::DimensionBoundary::VectR_N LocalPos(Points.GetM());
    NumElement.Fill(-1);
    if ((var_exp.mesh.GetNbElt() > 0) && (nb_points_global_measure > 0))
      {
	VectReal_wp EpsilonPoint; Vector<bool> point_treated;
	Vector<typename Dimension::MatrixN_N> mat_dfjm1;
	var_exp.GetOutputProblem().all_points_display.
	  LocalizePoints(var_exp.mesh, Points, NumElement, CoorLocal,
			 EpsilonPoint, point_treated, mat_dfjm1);

      }
    else if (var_exp.mesh.GetNbElt() > 0)
      for (int i = 0; i < var_exp.mesh.GetNbBoundaryRef(); i++)
	if (ref_cond(var_exp.mesh.BoundaryRef(i).GetReference()) == 1)
	  {
	    int num_face = i;
	    int num_elem = var_exp.mesh.BoundaryRef(i).numElement(0);
	    int num_loc = var_exp.mesh.Element(num_elem).GetPositionBoundary(num_face);
	    
	    FjInverseProblem<Dimension> inverseFj(var_exp.mesh, num_elem);
	    
	    // we try to localize points on elements num_elem
	    typename Dimension::DimensionBoundary::VectR_N pts; IVect num;
	    var_exp.LocalizePointsBoundaryElement(Points, num_elem, num_loc,
						  inverseFj, pts, num);
	    
	    for (int k = 0; k < num.GetM(); k++)
	      {
		NumElement(num(k)) = num_elem;
		NumLocal(num(k)) = num_loc;
		LocalPos(num(k)) = pts(k);
	      }
	  }
    
    IVect ProcExp;
    int nb_points = RemoveDuplicateLocalization(var_exp, NumElement, ProcExp);
    
    //DISP(nb_points); DISP(NumElement); DISP(NumLocal); DISP(LocalPos);
    
    // then the projection operator is constructed
    // proj_exp achieves the interpolation of experimental data on quadrature points
    // of the simulated case (or measure points)
    IVect num_points_exp(nb_points);
    int nb = 0; VectReal_wp phi;
    const MeshNumbering<Dimension>& mesh_num_exp = var_exp.GetMeshNumbering(0);
    proj_exp.Reallocate(nb_points, mesh_num_exp.GetNbDof());
    if (var_exp.mesh.GetNbElt() > 0)
      for (int k = 0; k < NumElement.GetM(); k++)
        if (NumElement(k) >= 0)
          {
            int num_elem = NumElement(k);
            int num_loc = NumLocal(k);
            const ElementReference<Dimension, 1>& Fb = var_exp.GetReferenceElementH1(num_elem);
            
            R_N pt_loc;
            if (nb_points_global_measure == 0)
              {
                Fb.GetLocalCoordOnBoundary(num_loc, LocalPos(k), pt_loc);
                Fb.ComputeValuesPhiRef(pt_loc, phi);
                for (int j = 0; j < Fb.GetNbDofBoundary(num_loc); j++)
                  {
                    int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
                    if (abs(phi(num_dof_loc)) > epsilon_machine)
                      {
                        int num_dof = mesh_num_exp.Element(num_elem).GetNumberDof(num_dof_loc);
                        proj_exp.AddInteraction(nb, num_dof, phi(num_dof_loc));
                      }
                  }
              }
            else
              {
                Fb.ComputeValuesPhiRef(CoorLocal(k), phi);
                for (int j = 0; j < Fb.GetNbDof(); j++)
                  if (abs(phi(j)) > epsilon_machine)
                    {
                      int num_dof = mesh_num_exp.Element(num_elem).GetNumberDof(j);
                      proj_exp.AddInteraction(nb, num_dof, phi(j));
                    }
              }
            
            num_points_exp(nb) = k;
            nb++;
          }
    
    IVect ProcElementS, num_points_simu;
    if (nb_points_global_measure > 0)
      {
	// case where quadrature points are provided by the user
	// these points need to be localized also for the simulated case
	typename Dimension::VectR_N CoorLocalS;
	IVect NumElementS(Points.GetM());
        ProcElementS.Reallocate(Points.GetM());
	NumElementS.Fill(-1);
        ProcElementS.Fill(nb_proc);
	
	VectReal_wp EpsilonPoint; Vector<bool> point_treated;
	Vector<typename Dimension::MatrixN_N> mat_dfjm1;
	var.GetOutputProblem().all_points_display.
	  LocalizePoints(var.mesh, Points, NumElementS, CoorLocalS,
			 EpsilonPoint, point_treated, mat_dfjm1);
	
        int nb_points_S = RemoveDuplicateLocalization(var, NumElementS, ProcElementS);
        //DISP(Points); DISP(nb_points_S); DISP(NumElementS); DISP(NumLocalS); DISP(LocalPosS);
        
        num_points_simu.Reallocate(nb_points_S);
	proj_simu.Reallocate(nb_points_S, mesh_num.GetNbDof());
        weights_measure.Reallocate(nb_points_S);
	nb = 0;
        VectReal_wp old_weights(weights_measure);
	for (int k = 0; k < NumElementS.GetM(); k++)
	  if (NumElementS(k) >= 0)
	    {
	      int num_elem = NumElementS(k);
	      const ElementReference<Dimension, 1>& Fb = var.GetReferenceElementH1(num_elem);
	      
	      Fb.ComputeValuesPhiRef(CoorLocalS(k), phi);
	      for (int j = 0; j < Fb.GetNbDof(); j++)
		if (abs(phi(j)) > epsilon_machine)
		  {
		    int num_dof = mesh_num.Element(num_elem).GetNumberDof(j);
		    proj_simu.AddInteraction(nb, num_dof, phi(j));
		  }
              
              num_points_simu(nb) = k;
              weights_measure(nb) = old_weights(k);
	      nb++;
	    }	
	
        if (rank_proc == 0)
          {
            // we check that all points have been found
            for (int k = 0; k < ProcElementS.GetM(); k++)
              if (ProcElementS(k) == nb_proc)
                {
                  cout << "Point " << Points(k) << " not localized in the simulated problem " << endl;
                  abort();
                }
          }     
	//proj_simu.WriteText("proj_simu.dat");
      }
    else
      {
	// weights are multiplied by r for axisymmetric case
	if (axisymmetric_problem)
	  for (int i = 0; i < mesh_interp.GetNbBoundary(); i++)
	    for (int j = 0; j < mesh_interp.GetNbPointsQuadrature(i); j++)
	      mesh_interp.WeightsQuadrature(i, j) *=
		mesh_interp.PointsQuadrature(i, j)(0);	

        // ProcElementS is given in offset_quadrature_per_proc
        ProcElementS.Reallocate(Points.GetM());
        const IVect& offset_proc = mesh_interp.offset_quadrature_per_proc;
        for (int p = 0; p < nb_proc; p++)
          for (int k = offset_proc(p); k < offset_proc(p+1); k++)
            ProcElementS(k) = p;

        num_points_simu.Reallocate(offset_proc(rank_proc+1) - offset_proc(rank_proc));
        for (int k = offset_proc(rank_proc); k < offset_proc(rank_proc+1); k++)
          num_points_simu(k-offset_proc(rank_proc)) = k;
        
      }
    
    if (var_exp.mesh.GetNbElt() == 0)
      return;
    
    //DISP(num_points_exp);
    //proj_exp.WriteText("proj_exp.dat");

    if (rank_proc == 0)
      {
	// we check that all points have been found
	for (int k = 0; k < ProcExp.GetM(); k++)
 	  if (ProcExp(k) == nb_proc)
	    {
	      cout << "Point " << Points(k) << " not localized in the experimental problem " << endl;
	      abort();
	    }
      }
    
#ifdef SELDON_WITH_MPI
    
    // determining which points to send
    IVect nb_points_to_send(nb_proc);
    nb_points_to_send.Zero();
    for (int k = 0; k < num_points_exp.GetM(); k++)
      {
        int p = num_points_exp(k);
        nb_points_to_send(ProcElementS(p))++;
      }

    num_points_to_send.Reallocate(nb_proc);
    for (int p = 0; p < nb_proc; p++)
      num_points_to_send(p).Reallocate(nb_points_to_send(p));
    
    nb_points_to_send.Zero();
    for (int k = 0; k < num_points_exp.GetM(); k++)
      {
        int p = num_points_exp(k);
        int proc = ProcElementS(p);
        int n = nb_points_to_send(proc);
        num_points_to_send(proc)(n) = k;
        nb_points_to_send(ProcElementS(p))++;
      }
      
    // determining which points to receive
    IVect nb_points_to_recv(nb_proc);
    nb_points_to_recv.Zero();
    for (int k = 0; k < num_points_simu.GetM(); k++)
      {
        int p = num_points_simu(k);
        nb_points_to_recv(ProcExp(p))++;
      }

    num_points_to_recv.Reallocate(nb_proc);
    for (int p = 0; p < nb_proc; p++)
      num_points_to_recv(p).Reallocate(nb_points_to_recv(p));
    
    nb_points_to_recv.Zero();
    for (int k = 0; k < num_points_simu.GetM(); k++)
      {
        int p = num_points_simu(k);
        int proc = ProcExp(p);
        int n = nb_points_to_recv(proc);
        num_points_to_recv(proc)(n) = k;
        nb_points_to_recv(ProcExp(p))++;
      }      
#endif
      
      
  }
    
  // projects the experience solution on quadrature points of the simulated reference
  template<class Dimension> template<class T>
  void VarMigration_Base<Dimension>
  ::ProjectExperimentData(VarHarmonic_Base<T, Dimension>& var_exp,
			  const Matrix<T, General, ColMajor>& rhs,
			  Matrix<T, General, ColMajor>& data_exp)
  {
    int nb_rhs = rhs.GetN();
    int nb_u = var_exp.nb_unknowns_scal;
    Matrix<T, General, ColMajor> data_proj(proj_exp.GetM(), nb_rhs*nb_u);
    data_proj.Zero();
      
    // projection on quadrature points that are contained in the current processor
    const MeshNumbering<Dimension>& mesh_num_exp = var_exp.GetMeshNumbering(0);
    int Nvol = mesh_num_exp.GetNbDof();
    Vector<T> x_sol, y;
    if (proj_exp.GetM() > 0)
      for (int k = 0; k < rhs.GetN(); k++)
	for (int m = 0; m < nb_u; m++)
	  {
	    int offset_rhs = var_exp.GetOffsetDofUnknown(m);
	    x_sol.SetData(Nvol, const_cast<T*>(&rhs(offset_rhs, k)));
	    y.SetData(data_proj.GetM(), &data_proj(0, k*nb_u + m));
	    Mlt(proj_exp, x_sol, y);
	      
	    x_sol.Nullify(); y.Nullify();
	  }
    
    //DISP(data_proj);
    
    // then we have to send the results to the correct processor
    int nb_proc = 1;
    
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_exp.comm_group_mode;
    MPI_Comm_size(comm, &nb_proc);
    int my_proc; MPI_Comm_rank(comm, &my_proc);
#endif
      
    if (nb_proc <= 1)
      {
	// for sequential execution data_proj is already good
	// we exchange data_proj and data_exp to avoid a duplication
	data_exp.SetData(data_proj.GetM(), data_proj.GetN(), data_proj.GetData());
	data_proj.Nullify();
	return;
      }
      
#ifdef SELDON_WITH_MPI
      
    nb_rhs *= nb_u;
      
    Vector<MPI_Request> request_send(nb_proc), request_recv(nb_proc);
    Vector<int64_t> xtmp; Vector<Vector<T> > x_send(nb_proc);
    Vector<Vector<T> > x_recv(nb_proc);
    // values are exchanged between processors
    for (int p = 0; p < nb_proc; p++)
      if (num_points_to_send(p).GetM() > 0)
	{
	  int nb_points = num_points_to_send(p).GetM()*nb_rhs;
	  x_send(p).Reallocate(nb_points);
          int nb = 0;
	  for (int j = 0; j < num_points_to_send(p).GetM(); j++)
	    for (int k = 0; k < nb_rhs; k++)
	      x_send(p)(nb++) = data_proj(num_points_to_send(p)(j), k);
          
	  if (p != my_proc)
	    request_send(p) = MpiIsend(comm, x_send(p), xtmp, nb_points, p, 126);
	}
      
    for (int p = 0; p < nb_proc; p++)
      if (num_points_to_recv(p).GetM() > 0)
	{
	  int nb_points = num_points_to_recv(p).GetM()*nb_rhs;
	  x_recv(p).Reallocate(nb_points);
	  if (p != my_proc)
	    request_recv(p) = MpiIrecv(comm, x_recv(p), xtmp, nb_points, p, 126);
	  else
	    x_recv(p) = x_send(p);
	}
      
    MPI_Status status;
    for (int p = 0; p < nb_proc; p++)
      if (p != my_proc)
	{
          if (num_points_to_send(p).GetM() > 0)
	    MPI_Wait(&request_send(p), &status);
	    
	  if (num_points_to_recv(p).GetM() > 0)
	    {
	      MPI_Wait(&request_recv(p), &status);
	      MpiCompleteIrecv(x_recv(p), xtmp, x_recv(p).GetM());
	    }
	}
      
    // then we fill data_exp
    int nb_points = this->GetNbPointsQuadrature();
    data_exp.Reallocate(nb_points, nb_rhs);
    for (int p = 0; p < nb_proc; p++)
      {
	//DISP(p); DISP(x_recv(p));
	int nb = 0;
	for (int j = 0; j < num_points_to_recv(p).GetM(); j++)
	  for (int k = 0; k < nb_rhs; k++)
	    data_exp(num_points_to_recv(p)(j), k) = x_recv(p)(nb++);
      }
      
#endif
      
    //DISP(data_exp);
  }
    

  template<class Dimension> template<class T>
  void VarMigration_Base<Dimension>
  ::ProjectSimulationData(VarHarmonic_Base<T, Dimension>& var_simu,
			  const Matrix<T, General, ColMajor>& rhs,
			  Matrix<T, General, ColMajor>& data_simu)
  {
    int nb_u = var_simu.nb_unknowns_scal, nb_pts;
    if (nb_points_global_measure == 0)
      nb_pts = mesh_interp.GetNbPointsQuadrature();
    else
      nb_pts = proj_simu.GetM();

    const MeshNumbering<Dimension>& mesh_num_simu = var_simu.GetMeshNumbering(0);
    int Nvol = mesh_num_simu.GetNbDof();
    Vector<T> x_sol, y, trace_du_dn;
    data_simu.Reallocate(nb_pts, rhs.GetN()*nb_u);
    if (nb_pts > 0)
      for (int k = 0; k < rhs.GetN(); k++)
	{
	  if (nb_points_global_measure == 0)
	    {
	      x_sol.SetData(rhs.GetM(), const_cast<T*>(&rhs(0, k)));
	      mesh_interp.ComputeEnHnOnBoundary(var_simu, x_sol, y,
						trace_du_dn, false, false);
	      
	      for (int m = 0; m < nb_u; m++)
		{
		  int kcol = k*nb_u + m;
		  for (int j = 0; j < nb_pts; j++)
		    data_simu(j, kcol) = y(nb_u*j + m);
		}
	      
	      x_sol.Nullify();
	    }
	  else
	    for (int m = 0; m < nb_u; m++)
	      {
		int offset_rhs = var_simu.GetOffsetDofUnknown(m);
		x_sol.SetData(Nvol, const_cast<T*>(&rhs(offset_rhs, k)));
		y.SetData(proj_simu.GetM(), &data_simu(0, k*nb_u + m));
		Mlt(proj_simu, x_sol, y);
		
		x_sol.Nullify(); y.Nullify();
	      }
	}
    
    //DISP(data_simu);
    //data_simu.Write("data_simu.dat");
  }
    
    
  // completes production of experimental measurements
  // var_exp : object used to perform experiences
  // solver_exp : linear solver
  // data_exp : values of experience solutions on quadrature points of the simulated problem
  template<class Dimension>
  void VarMigration_Base<Dimension>
  ::ComputeExperiment(VarHarmonic_Base<Complex_wp, Dimension>& var_exp,
		      All_LinearSolver* solver_exp, Matrix<Complex_wp, General, ColMajor>& data_exp)
  {
    // source are modified with the parameters stored in the current class
    var_exp.SetParameterSource(param_source);
    
    // computation of all sources in a matrix (each column corresponds to a different source)
    Matrix<Complex_wp, General, ColMajor> rhs;
    var_exp.ComputeRightHandSide(rhs);
    
    // finite element matrix is computed and factorized
    GlobalGenericMatrix<Complex_wp> nat_mat;
    solver_exp->PerformFactorizationStep(nat_mat);
      
    // the solutions are computed (rhs is overwritten by solutions)
    solver_exp->ComputeSolution(rhs, nat_mat);
    solver_exp->ClearFactorization();
      
    if (write_intermediate_output)
      {
	VectComplex_wp x_sol(var_exp.GetNbDof());
	GetCol(rhs, 0, x_sol);    
	var_exp.GetOutputProblem().WriteOutputFile(x_sol, 1);
      }
    
    // the projection on quadrature points is computed and stored in data_exp
    ProjectExperimentData(var_exp, rhs, data_exp);
      
    //DISP(data_exp);
    //data_exp.Write("data_exp.dat");
  }
    
    
  template<class Dimension>
  void VarMigration_Base<Dimension>
  ::ComputeRightHandSide(VarHarmonic_Base<Complex_wp, Dimension>& var_simu,
			 Matrix<Complex_wp, General, ColMajor>& rhs)
  {
    // parameters describing sources are replaced by parameters stored in the current object
    var_simu.SetParameterSource(param_source);
      
    // the right hand sides are computed (as a matrix, since there can be several rhs)
    var_simu.ComputeRightHandSide(rhs);
  }

  
  // completes a simulation
  // var_simu : object containing simulation problem
  // solver_simu : linear solver used to factorize finite element matrix
  // data_simu : projection of the solutions on quadrature points 
  // rhs : solutions 
  template<class Dimension>
  void VarMigration_Base<Dimension>
  ::ComputeSimulation(VarHarmonic_Base<Complex_wp, Dimension>& var_simu,
		      All_LinearSolver* solver_simu,
		      Matrix<Complex_wp, General, ColMajor>& data_simu,
		      Matrix<Complex_wp, General, ColMajor>& rhs)
  {
    ComputeRightHandSide(var_simu, rhs);
      
    // the finite element matrix is computed and factorized
    GlobalGenericMatrix<Complex_wp> nat_mat;
    solver_simu->PerformFactorizationStep(nat_mat);
      
    // rhs is overwritten by solutions
    solver_simu->ComputeSolution(rhs, nat_mat);
      
    if (write_intermediate_output)
      {
	VectComplex_wp x_solb(var_simu.GetNbDof());
	GetCol(rhs, 0, x_solb);
	var_simu.GetOutputProblem().WriteOutputFile(x_solb, 0);
      }
    
    // the values on the quadrature points are computed
    ProjectSimulationData(var_simu, rhs, data_simu);
  }
    
  template<class Dimension>
  void VarMigration_Base<Dimension>
  ::ComputeIntegralResidu(VarProblem<Dimension>& var_simu,
			  Matrix<Complex_wp, General, ColMajor>& residu,
			  Matrix<Complex_wp, General, ColMajor>& sol_residu)
  {    
    VectComplex_wp r, x, rhs, r2(residu.GetM());
    const MeshNumbering<Dimension>& mesh_num_simu = var_simu.GetMeshNumbering(0);
    int Nvol = mesh_num_simu.GetNbDof();
    int nb_u = var_simu.nb_unknowns_scal;
    int nb_rhs = residu.GetN()/nb_u;
    sol_residu.Reallocate(var_simu.GetNbDof(), nb_rhs);
    sol_residu.Zero();
    for (int k = 0; k < nb_rhs; k++)
      {
	rhs.SetData(var_simu.GetNbDof(), &sol_residu(0, k));
	
	if (nb_points_global_measure > 0)
	  {
	    for (int m = 0; m < nb_u; m++)
	      {
		int offset_rhs = var_simu.GetOffsetDofUnknown(m);
		r.SetData(residu.GetM(), &residu(0, k*nb_u + m));
		x.SetData(Nvol, &sol_residu(offset_rhs, k)); 
		
		for (int p = 0; p < r2.GetM(); p++)
		  r2(p) = r(p)*weights_measure(p);
		
		Mlt(SeldonTrans, proj_simu, r2, x);
		r.Nullify(); x.Nullify();
	      }
	  }
	else
	  {
	    if (mesh_interp.GetNbPointsQuadrature() > 0)
	      for (int m = 0; m < nb_u; m++)
		{
		  int offset_rhs = var_simu.GetOffsetDofUnknown(m);
		  r.SetData(residu.GetM(), &residu(0, k*nb_u + m));
		  x.SetData(Nvol, &sol_residu(offset_rhs, k)); 
		  mesh_interp.AddSourceBoundary(var_simu, r, x);
		  r.Nullify(); x.Nullify();
		}
	    
	    var_simu.AddDomains(rhs, nb_u);
	  }
	
	rhs.Nullify();
      }
  }
  

  // computes adjoint solutions
  // var_simu : object containing datas about the solved problem
  // solver_simu : linear solver used to solve finite element matrix
  // residu : Neumann data on quadrature points (they are acting as sources for the adjoint problem)
  // sol_residu : solutions
  template<class Dimension>
  void VarMigration_Base<Dimension>
  ::ComputeAdjoint(VarProblem<Dimension>& var_simu, All_LinearSolver* solver_simu,
		   Matrix<Complex_wp, General, ColMajor>& residu,
		   Matrix<Complex_wp, General, ColMajor>& sol_residu)
  {
    // Neumann source is conjugated
    Conjugate(residu);
      
    // computing integral \int residu \varphi, the result is placed in sol_residu
    ComputeIntegralResidu(var_simu, residu, sol_residu);
    //DISP(MaxAbs(residu)); DISP(MaxAbs(sol_residu));
    
    // then the solutions are computed
    GlobalGenericMatrix<Complex_wp> nat_mat;
    solver_simu->ComputeSolution(sol_residu, nat_mat);

    if (write_intermediate_output)
      {
	VectComplex_wp x_solb(var_simu.GetNbDof());
	GetCol(sol_residu, 0, x_solb);
	var_simu.GetOutputProblem().WriteOutputFile(x_solb, 2);
      }
  }

}

#define MONTJOIE_FILE_VAR_MIGRATION_CXX
#endif

