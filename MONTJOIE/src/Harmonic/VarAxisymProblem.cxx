#ifndef MONTJOIE_FILE_VAR_AXISYM_PROBLEM_CXX

namespace Montjoie
{

#ifdef MONTJOIE_WITH_TWO_DIM
  /********************
   * VarAxisymProblem *
   ********************/


  //! Sets values of attributes to default values
  void VarAxisymProblem::InitDefaultValues()
  {
    number_mode_to_be_computed = true;
    threshold_mode = 1e-6;
    kwave3D.Fill(0);
    kwave3D(0) = var_problem.GetOmega();
    kwave_envelope = var_problem.GetOmega();
    force_diagonal_mass_invertible = false;
  }


  void VarAxisymProblem::GetMemoryUsed(map<string, size_t>& var) const
  {
    size_t taille = sizeof(*this);
    taille += Dof_On_Axe.GetMemorySize() + maximum_mode_for_each_incidence.GetMemorySize()
      + Seldon::GetMemorySize(Jn_precomputed) + Vertex_On_Axe.GetMemorySize()
      + Element_On_Axe.GetMemorySize();

    var["Axisymmetric Data"] = taille;
  }


  //! reads keywords related to axisymmetric problems
  void VarAxisymProblem
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    // 3-D incident plane wave, we need two angles teta phi (instead of one angle in 2-D)
    if (!description_field.compare("IncidentAngle"))
      {
        if (parameters.GetM() <= 1)
          {
            cout << "In SetInputData of VarHelmholtz_Axi" << endl;
            cout << "IncidentAngle needs 2 parameters, for instance :" << endl;
            cout << "IncidentAngle = teta phi" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        // in 3-D, the wave vector is oriented along 
        // (sin(teta)*cos(phi) ; sin(teta)*sin(phi) ; cos(teta))
        Real_wp teta = 0.0, phi = 0.0;
        teta = to_num<Real_wp>(parameters(0))*pi_wp/180;
        phi = to_num<Real_wp>(parameters(1))*pi_wp/180;
        
        Real_wp c0 = var_problem.GetVelocityOfInfinity();
        Real_wp omega = var_problem.GetOmega();
        kwave3D(0) = omega/c0*sin(teta)*cos(phi);
        kwave3D(1) = omega/c0*sin(teta)*sin(phi);
        kwave3D(2) = omega/c0*cos(teta);
        kwave_envelope = omega/c0;
      }
    else if (!description_field.compare("NumberModes"))
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of VarHelmholtz_Axi" << endl;
            cout << "NumberModes needs one parameter, for instance :" << endl;
            cout << "NumberModes = N" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        IVect list_number_mode(1);
        list_number_mode(0) = 0;

        // number of modes to take into account in the Fourier expansion of the solution over theta
        if (!parameters(0).compare("AUTO"))
          {
            if (parameters.GetM() <= 1)
              {
                cout << "In SetInputData of VarHelmholtz_Axi" << endl;
                cout << "NumberModes needs 2 parameters, for instance :" << endl;
                cout << "NumberModes = AUTO threshold" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }
            
            // in that case the computation starts from 0 until x is small enough
            number_mode_to_be_computed = true;
            threshold_mode = to_num<Real_wp>(parameters(1));
            if (parameters.GetM() > 2)
              {
                // you can specify a minimal number of modes to explore
                int number_maximum_mode = to_num<int>(parameters(2));
                list_number_mode.Reallocate(2*number_maximum_mode+1);
                list_number_mode(0) = 0;
                int nb = 1;
                for (int i = 1; i <= number_maximum_mode; i++)
                  {
                    list_number_mode(nb) = -i;
                    list_number_mode(nb+1) = i;
                    nb++;
                  }
              }
          }
        else if (!parameters(0).compare("SINGLE"))
          {
            if (parameters.GetM() <= 1)
              {
                cout << "In SetInputData of VarHelmholtz_Axi" << endl;
                cout << "NumberModes needs 2 parameters, for instance :" << endl;
                cout << "NumberModes = SINGLE mode_number" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }
            
            int number_maximum_mode = 0;
            to_num(parameters(1), number_maximum_mode);
            list_number_mode.Reallocate(1);
            list_number_mode(0) = number_maximum_mode;
            
            number_mode_to_be_computed = false;
          }
        else
          {
            int number_maximum_mode = 0;
            // to_num(parameters(0), number_maximum_mode);
            number_maximum_mode = to_num<int>(parameters(0)); //NATHAN
            list_number_mode.Reallocate(2*number_maximum_mode+1);
            list_number_mode(0) = 0;
            int nb = 1;
            for (int i = 1; i <= number_maximum_mode; i++)
              {
                list_number_mode(nb) = -i;
                list_number_mode(nb+1) = i;
                nb += 2;
              }
            
            number_mode_to_be_computed = false;
          }
        
        var_boundary.SetModesToCompute(list_number_mode);
      }
    else if (!description_field.compare("OriginePhase"))
      {
        if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of VarHelmholtz_Axi" << endl;
            cout << "OriginePhase needs 3 parameters, for instance :" << endl;
            cout << "OriginePhase = x0 y0 z0" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        Real_wp x0 = to_num<Real_wp>(parameters(0));
        Real_wp y0 = to_num<Real_wp>(parameters(1));
        Real_wp z0 = to_num<Real_wp>(parameters(2));
        origine_phase3D.Init(x0, y0, z0);        
      }
    else if (!description_field.compare("ForceDiagonalMass"))
      {
        if (parameters(0) == "YES")
          force_diagonal_mass_invertible = true;
        else
          force_diagonal_mass_invertible = false;
      }
  }
  

  //! Performs mesh treatment for axisymmetric problems
  void VarAxisymProblem::CheckSectionMeshAxi()
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    if (rank_proc == 0)
      {
        var_problem.mesh.SymmetrizeMeshToGetPositiveX();
        
        if (var_problem.mesh.GetNbElt() == 0)
          {
            cout<<"The mesh must have only negative x or only positive x"<<endl;
            if (nb_proc == 1)
              var_problem.mesh.Write("test.mesh");
            
            abort();
          }
      }
    
    // specific treatment for axisymmetric domains
    // we remove from the referenced edges, edges which are on the axis of revolution
    // and retrieve dofs on the axis as well
    IVect list_vertices_axis, list_edges_axis;
    var_problem.mesh.GetEdgesOnZaxis(list_vertices_axis, list_edges_axis, Vertex_On_Axe);

    // elements adjacent to the axis ?
    Element_On_Axe.Reallocate(var_problem.mesh.GetNbElt());
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
        Element_On_Axe(i) = false;
        for (int j = 0; j < var_problem.mesh.Element(i).GetNbVertices(); j++)
          if (Vertex_On_Axe(var_problem.mesh.Element(i).numVertex(j)))
            Element_On_Axe(i) = true;
      }
    
    if ((rank_proc == 0) && (this->force_diagonal_mass_invertible))
      {
        // reorienting elements close to the axis
        // such that x = 0 correspond to the axis for each quadrilateral
        for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
          if (Element_On_Axe(i))
            {
              int nv = var_problem.mesh.Element(i).GetNbVertices();
              if (nv != 4)
                {
                  cout << "triangles not allowed with mass lumping" << endl;
                  abort();
                }
              
              int n0 = var_problem.mesh.Element(i).numVertex(0);
              int n1 = var_problem.mesh.Element(i).numVertex(1);
              int n2 = var_problem.mesh.Element(i).numVertex(2);
              int n3 = var_problem.mesh.Element(i).numVertex(3);
              int ref = var_problem.mesh.Element(i).GetReference();
              
              if (Vertex_On_Axe(n0) && Vertex_On_Axe(n1))
                var_problem.mesh.Element(i).InitQuadrangular(n1, n2, n3, n0, ref);
              else if (Vertex_On_Axe(n1) && Vertex_On_Axe(n2))
                var_problem.mesh.Element(i).InitQuadrangular(n2, n3, n0, n1, ref);
              else if (Vertex_On_Axe(n2) && Vertex_On_Axe(n3))
                var_problem.mesh.Element(i).InitQuadrangular(n3, n0, n1, n2, ref);
              else if (Vertex_On_Axe(n3) && Vertex_On_Axe(n0))
                var_problem.mesh.Element(i).InitQuadrangular(n0, n1, n2, n3, ref);
              else
                {
                  cout << "Case not treated" << endl;
                  cout << "All the elements that touch the axis must have an edge completely on the axis" << endl;
                  var_problem.mesh.Write("test.mesh");
                  abort();
                }
            }
        
        var_problem.mesh.ReorientElements();
        var_problem.mesh.FindConnectivity();
        var_problem.mesh.ProjectPointsOnCurves();
      }
  }
  
  
  //! computes dofs on the axis Oz
  void VarAxisymProblem::ComputeDofOnAxe(const VarProblem<Dimension2>& var_pb)
  {
    // dofs located on the axis are retrieved
    // if IsDofOnAxe(i) is true, dof i is located on the z-axis
    Vector<bool> IsDofOnAxe(var_problem.GetNbDof());
    IsDofOnAxe.Fill(false);
    const MeshNumbering<Dimension2>& mesh_num = var_problem.GetMeshNumbering(0);
    if (!var_pb.FirstOrderFormulationDG())
      for (int i = 0; i < var_problem.mesh.GetNbBoundaryRef(); i++)
        {
          int n1 = var_problem.mesh.BoundaryRef(i).numVertex(0),
            n2 = var_problem.mesh.BoundaryRef(i).numVertex(1);

#ifdef MONTJOIE_WITH_TWO_DIM
          if (Vertex_On_Axe(n1)&&Vertex_On_Axe(n2))
            {
              // dofs on this edge
              int ne = i;
              int num_elem = var_problem.mesh.Boundary(ne).numElement(0);
              int num_loc = var_problem.mesh.Element(num_elem).GetPositionBoundary(ne);
              const ElementReference_Dim<Dimension2>& Fb = var_pb.GetReferenceElement(num_elem);
              for (int j = 0; j < Fb.GetNbDofBoundary(num_loc); j++)
                {
                  int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
                  int num_dof = mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
                  IsDofOnAxe(num_dof) = true;
                }
            }
#endif
        }
    
    const Vector<int>& DdlVol = var_gibc_base.GetLocalVolumeDofNumber();
    int Nvol = mesh_num.GetNbDof();
    for (int i = 0; i < DdlVol.GetM(); i++)
      if (IsDofOnAxe(DdlVol(i)))
        IsDofOnAxe(Nvol+i) = true;
    
    int nb_dof_axe = 0;
    for (int i = 0; i < IsDofOnAxe.GetM(); i++)
      if ((IsDofOnAxe(i)) && (!var_boundary.IsDofDirichlet(i)))
        nb_dof_axe++;
    
    Dof_On_Axe.Reallocate(nb_dof_axe);
    nb_dof_axe = 0;
    for (int i = 0; i < IsDofOnAxe.GetM(); i++)
      if ((IsDofOnAxe(i)) && (!var_boundary.IsDofDirichlet(i)))
        Dof_On_Axe(nb_dof_axe++) = i;    
  }
  
  
  //! Computes k_perp, kz and phase
  void VarAxisymProblem
  ::Get_KwavePerp_Kz_Phase(const Complex_wp& rho_tilde, const Complex_wp& mu, int m,
                           const R3& kwave, const Real_wp& omega,
                           Real_wp& k_perp, Real_wp& kz, bool& incidence_axial, Complex_wp& phase)
  {
    Real_wp invC0 = sqrt(abs(rho_tilde/mu));
    
    Real_wp kx = kwave(0)*invC0;
    Real_wp ky = kwave(1)*invC0;
    k_perp = sqrt(kx*kx + ky*ky);
    kz = kwave(2)*invC0;
    Real_wp teta0 = 0.0;
    incidence_axial = true;    
    if (k_perp/omega > 10.0*epsilon_machine)
      {
        incidence_axial = false; 
        teta0 = acos(kx/k_perp);
        if (ky < 0)
          teta0 = -teta0;
      }
    
    // computing phase = i^-m exp(im teta0)
    phase = Complex_wp(cos(m*teta0), sin(m*teta0));
    phase *= ComputePowerI(-m);
  }
  

  //! computation of the list of modes necessary to obtain an accurate solution
  void VarAxisymProblem
  ::ComputeListMode(VarComputationRCS_Axi& rcs_param)
  {
    Real_wp k_bot; int number_maximum_mode = 1;
    int nb_source = this->GetNbRightHandSide(rcs_param);
    IVect list_number_mode;
    if (nb_source > 1)
      {
        maximum_mode_for_each_incidence.Reallocate(nb_source);
        for (int n = 0; n < nb_source; n++)
          {
            R3 kwave, polar;
            rcs_param.GetWaveVectorMonostatic(n, var_problem.GetOmega(), kwave, polar);
            Real_wp k_bot = sqrt(kwave(0)*kwave(0) + kwave(1)*kwave(1));
            maximum_mode_for_each_incidence(n) = ComputeNbModes_Generic(k_bot);
            number_maximum_mode = max(number_maximum_mode, maximum_mode_for_each_incidence(n));
          }
        
        list_number_mode.Reallocate(2*number_maximum_mode+1);
        list_number_mode(0) = 0;
        for (int i = 1; i <= number_maximum_mode; i++)
          {
            list_number_mode(2*i-1) = -i;
            list_number_mode(2*i) = i;
          }
      }
    else
      {
        k_bot = sqrt(square(kwave3D(0)) + square(kwave3D(1))); // DISP(k_bot); 
        if (k_bot/var_problem.GetOmega() > 100.0*epsilon_machine)
          {
            number_maximum_mode = ComputeNbModes_Generic(k_bot);
            if (var_problem.print_level >= 1)
              cout<<"we take into account "<< number_maximum_mode<<" modes"<<endl;
            
            list_number_mode.Reallocate(2*number_maximum_mode+1);
            list_number_mode(0) = 0;
            for (int i = 1; i <= number_maximum_mode; i++)
              {
                list_number_mode(2*i-1) = -i;
                list_number_mode(2*i) = i;
              }
          }
        else
          {
            // for vectorial problems
            list_number_mode.Reallocate(2);
            list_number_mode(0) = -1; list_number_mode(1) = 1;
          }
      }    
    
    var_boundary.SetModesToCompute(list_number_mode);
  }
  
  
  //! returns an estimation of the number of modes involved in computation
  int VarAxisymProblem::ComputeNbModes_Generic(const Real_wp& k_bot)
  {
    int nb_modes = 0;
    // maximum of argument in bessel function :
    Real_wp kr = k_bot*var_problem.GetXmax();

    Real_wp lb = GetLambertW0(2.0/(3.0*pi_wp*square(var_source.GetThresholdSource())));
    nb_modes = toInteger(kr + 0.5*pow(1.5*lb, Real_wp(2)/3)*pow(kr,Real_wp(1)/3))+1;

    return nb_modes;
  }


  //! precomputes Bessel functions on quadrature points
  void VarAxisymProblem
  ::InitBesselArray(VarComputationRCS_Axi& rcs_param)
  {
    int number_maximum_mode = 0;
    for (int i = 0; i < var_boundary.GetNbModes(); i++)
      number_maximum_mode = max(number_maximum_mode, abs(var_boundary.GetModeNumber(i)));
    
    // computation of Bessel function on all quadrature points and for each incident direction
    Jn_precomputed.Reallocate(number_maximum_mode+2);
    
    const MeshInterpolationFEM<Dimension2>& var_mesh = rcs_param.GetInterpolationMesh();
    int nb_quadrature_points = var_mesh.GetNbAllQuadraturePoints();     
    for (int i = 0; i <= (number_maximum_mode+1); i++)
      Jn_precomputed(i).Reallocate(rcs_param.GetNbAngles()*nb_quadrature_points);
    
    Real_wp kr; R3 ur;
    VectReal_wp Jn(number_maximum_mode+2);
    
    for (int num_angle = 0; num_angle < rcs_param.GetNbAngles(); num_angle++)
      {
        ur = rcs_param.incident_direction(num_angle);
        Real_wp kbot = var_problem.GetOmega()*sqrt(ur(0)*ur(0)+ur(1)*ur(1));
        // Real_wp kz = this->GetOmega()*sin(teta);
        // loop on all quadrature points
        for (int ind = 0; ind < nb_quadrature_points; ind++)
          {
            kr = kbot*var_mesh.GetQuadraturePoint(ind)(0);
            ComputeBesselFunctions(0, number_maximum_mode+2, kr, Jn);
            
            for (int i = 0; i <= (number_maximum_mode+1); i++)
              {
                int num_point = num_angle*nb_quadrature_points+ind;
                Jn_precomputed(i)(num_point) = Jn(i);
              }
          }
      }
  }    
  
  
  //! constructs quadrature for radar cross sections
  void VarAxisymProblem
  ::InitRcs(VarComputationRCS_Axi& rcs_param)
  {
    rcs_param.InitComputationRCS(false);
    rcs_param.InitIncidentDirection();
  }
  
  
  //! only one source
  int VarAxisymProblem
  ::GetNbRightHandSide(VarComputationRCS_Axi& rcs_param) const
  {
    if (rcs_param.GetRcsType() == rcs_param.MONOSTATIC_RCS)
      return 2*rcs_param.GetNbAngles();
    
    return 1;
  }
  
  
  /*************************
   * VarComputationRCS_Axi *
   *************************/


  //! computes all incident vectors
  void VarComputationRCS_Axi
  ::InitIncidentDirection()
  {
    // direction is taken in the plane (Etilde, ktilde) where
    // Etilde is the unit vector oriented along the polarization of the incident plane wave
    // ktilde is the unit vector oriented along the wave vector
    R3 Etilde, ktilde; 
    var_problem.GetPolarization(Etilde);
    
    Mlt(1.0/Norm2(Etilde), Etilde);
    ktilde = var_axisym.GetWaveVector(); Mlt(1.0/Norm2(ktilde), ktilde);
    
    R3 ur;
    int nbPointsRCS = this->nb_angles_RCS;
    Real_wp step_angle = (this->last_angle_RCS-this->first_angle_RCS) / (nbPointsRCS-1);
    incident_direction.Reallocate(nbPointsRCS);
    for (int num_angle = 0; num_angle < nbPointsRCS; num_angle++)
      {
        Real_wp teta = this->first_angle_RCS + step_angle*num_angle;
        
        ur(0) = cos(teta)*ktilde(0) + sin(teta)*Etilde(0); 
        ur(1) = cos(teta)*ktilde(1) + sin(teta)*Etilde(1); 
        ur(2) = cos(teta)*ktilde(2) + sin(teta)*Etilde(2);
        
        incident_direction(num_angle) = ur;
      }
  }
  
  
  //! sets the wave vector and polarization of the k-th incident wave
  void VarComputationRCS_Axi
  ::GetWaveVectorMonostatic(int k, const Real_wp& omega, R3& kwave, R3& polar)
  {
    kwave = omega*incident_direction(k/2);
    if (abs(kwave(1)) > 1e-12)
      {
        cout << "Only implemented for wave vectors in the (x, z) plane" << endl;
        abort();
      }
    
    if (k%2 == 0)
      polar.Init(-incident_direction(k/2)(2), 0, incident_direction(k/2)(0));
    else
      polar.Init(0, Real_wp(1), 0);
  }
#endif

}

#define MONTJOIE_FILE_VAR_AXISYM_PROBLEM_CXX
#endif

