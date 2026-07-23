#ifndef MONTJOIE_FILE_VIBRO_ACOUSTIC_SESSION_CXX

namespace Montjoie
{
  
  //! we read a line of the data file
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  SetInputData(const string& keyword, const Vector<string>& param)
  {
    if(!keyword.compare("Density"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of VibroAcousticSession" << endl;
	    cout << "Density needs more parameters, for instance :" << endl;
	    cout << "Density = rho" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

        rho = to_num<Real_wp>(param(0));
      }
    else if (!keyword.compare("Velocity"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of VibroAcousticSession" << endl;
	    cout << "Velocity needs more parameters, for instance :" << endl;
	    cout << "Velocity = c" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

        c2 = square(to_num<Real_wp>(param(0)));
      }
    else if (!keyword.compare("ReferencePlate"))
      {
        ref_plate.Reallocate(param.GetM());
        for (int i = 0; i < ref_plate.GetM(); i++)
          ref_plate(i) = to_num<int>(param(i));
      }
    else if (keyword.find("Mecanic") == 0)
      {
        // for all keywords beginning with Mecanic
        // we call the method SetInputData of var_plate
        string keyword2 = keyword.substr(7);
        var_plate.SetInputData(keyword2, param);
      }
    else if (keyword.find("Acoustic") == 0)
      {
        // for all keywords beginning with Acoustic
        // we call the method SetInputData of var_volume
        string keyword2 = keyword.substr(8);
        var_volume.SetInputData(keyword2, param);
      }
    else
      {
        var_volume.SetInputData(keyword, param);
        var_plate.SetInputData(keyword, param);
      }
  }
  
  
  //! computation of mesh and finite element
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  ComputeMeshAndFiniteElement()
  {
    Mesh<Dimension3>& mesh = var_volume.var_harmonic.mesh;
    MeshNumbering<Dimension3>& mesh_num = var_volume.var_harmonic.GetMeshNumbering(0);
    Mesh<Dimension2>& mesh_plate = var_plate.var_harmonic.mesh;
    MeshNumbering<Dimension2>& mesh_plate_num = var_plate.var_harmonic.GetMeshNumbering(0);
    
    // copying order of approximation for 2-D mesh
    var_volume.var_harmonic.var_transmission.SetMeshPlate(mesh_plate);
    var_volume.var_harmonic.var_transmission.DuplicateExtremities(false);
    
    if (var_volume.Glob_mat_Dh == NULL)
      var_volume.Glob_mat_Dh = var_volume.GetNewMassMatrix();
    
    var_volume.Glob_mat_Dh->compute_Dh_default = true;
    var_volume.Glob_mat_Dh->compute_Sh_default = true;
    if (var_volume.FirstOrderScheme())
      var_volume.var_harmonic.SetFirstOrderFormulation(true);
    else
      var_volume.var_harmonic.SetFirstOrderFormulation(false);
    
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm_vol = var_volume.var_harmonic.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm_vol, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm_vol, &rank_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    // we construct the 3-D mesh
    if (rank_proc == 0)
      {
	cout << "Constructing the 3-D mesh ... " << endl;
        Vector<string> mesh_param = var_volume.var_harmonic.mesh_data(0);
        mesh.ConstructMesh(0, mesh_param);
        
        // reference of the plate
        var_volume.var_harmonic.var_transmission.SetReferences(ref_plate);
        
        // now we extract 2-D mesh of the plate, and 3-D vertices are duplicated
	cout << "Separating volume elements on the plate ..." << endl;
        var_volume.var_harmonic.var_transmission.PartMeshTransmission();
        
        // construction of 2-D mesh
        mesh_plate = var_volume.var_harmonic.var_transmission.GetMeshPlate();

	if (var_plate.var_harmonic.print_level >= 4)
	  mesh_plate.Write("plaque.msh");
	
        if (nb_proc == 1)
          {            
            // construction of 3-D finite element
            var_volume.var_harmonic.ConstructFiniteElement(string("TETRAHEDRON_LOBATTO"));
            
            // and 2-D finite element
            var_plate.var_harmonic.ConstructFiniteElement(string("TRIANGLE_LOBATTO"));
            
            var_volume.var_harmonic.FindElementsInsidePML();
            var_plate.var_harmonic.FindElementsInsidePML();
            
            mesh_num.NumberMesh();
            mesh_plate_num.NumberMesh();

            IVect MinimalProc;
            var_volume.var_harmonic.var_transmission.SetMeshPlate(mesh_plate);
            var_volume.var_harmonic.var_transmission.SetMeshPlateNumbering(mesh_plate_num);
            var_volume.var_harmonic.var_transmission.TreatTransmission(MinimalProc);

            var_plate.var_harmonic.ComputeNumberOfDofs();
            var_volume.var_harmonic.ComputeNumberOfDofs();            
          }        
      }
    
    if (nb_proc > 1)
      {
	if (rank_proc == 0)
	  cout << "Splitting 3-D mesh ... " << endl;
	
        // both 2-D and 3-D meshes are split into several subdomains
        var_volume.var_harmonic.SplitMeshForParallelComputation(string("TETRAHEDRON_LOBATTO"));
        
	if (rank_proc == 0)
	  cout << "Splitting 2-D mesh ..." << endl;
	else
          {
            var_volume.var_harmonic.mesh.print_level = 0;
            var_plate.var_harmonic.mesh.print_level = 0;
          }

        var_plate.var_harmonic.SplitMeshForParallelComputation(string("TRIANGLE_LOBATTO"));
      }
    
    if (rank_proc == 0)
      cout << "Exchanging volume and surface dofs between processors..." << endl;
    
    this->DistributeTransmissionDofs();
    
    // construction of projector from mecanic dofs to acoustic dofs
    VectReal_wp step_subdiv = mesh.GetQuadrilateralReferenceElement().PointsNodal1D();
    Vector<VectR2> points_div(2);
    points_div(0) = mesh.GetTriangleReferenceElement().PointsNodalND();
    points_div(1) = mesh.GetQuadrilateralReferenceElement().PointsNodalND();

    Vector<const ElementGeomReference<Dimension2>* > elt_plate;
    var_plate.var_harmonic.mesh.GetReferenceElementVolume(elt_plate);
    proj_meca.InitProjection(elt_plate, step_subdiv, points_div);
    
    if (rank_proc == 0)
      cout << "Computation of 2-D and 3-D interpolation grids..." << endl;

    var_plate.var_harmonic.InitVarGrid();
    var_plate.var_harmonic.ComputeVarGrid();
    
    // computation of interpolation grids
    var_volume.var_harmonic.InitVarGrid();
    var_volume.var_harmonic.ComputeVarGrid();
    
    if (rank_proc == 0)
      cout << "Computing number of degrees of freedom ..." << endl;
    
    var_plate.var_harmonic.PerformOtherInitializations();
    var_volume.var_harmonic.PerformOtherInitializations();
  }
  
  
  //! function used to link surface and volume dofs in parallel
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  DistributeTransmissionDofs()
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm_vol = var_volume.var_harmonic.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm_vol, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm_vol, &rank_proc);
#else
    int nb_proc(1);
#endif

    if (nb_proc == 1)
      {
	var_volume.var_harmonic.var_transmission.
	  CopyArray(ddl_vol_minus, ddl_vol_plus, NodleAcous, NodleMeca,
		    JacobianAcous, DirichletAcous);
	
	nodl_acous = ddl_vol_minus.GetM();
	
        return;
      }
    
#ifdef SELDON_WITH_MPI
    typedef typename TypeEqVol::TypeEquationStationary TypeEqVolS;
    VarTransmission<TypeEqVolS>& var_t = var_volume.var_harmonic.var_transmission;
    MeshNumbering<Dimension3>& mesh_num = var_volume.var_harmonic.GetMeshNumbering(0);
    Mesh<Dimension2>& mesh_plate = var_plate.var_harmonic.mesh;
    MeshNumbering<Dimension2>& mesh_plate_num = var_plate.var_harmonic.GetMeshNumbering(0);

    Vector<bool> proc_surf(nb_proc);
    proc_surf.Fill(true);
    
    IVect LocProcVolMinus, LocProcVolPlus, LocNumDofAcous, LocDirAcous, LocNbDofAcous;
    VectReal_wp LocJacobAcous;
    Vector<IVect> NumElem_Subdomain(nb_proc);
    int nb_ddl_max = 0;
    MPI_Status status;
    
    //mesh.Write("volume"+to_str(rank_proc)+".mesh");
    //mesh_plate.Write("surface"+to_str(rank_proc)+".mesh");
    Vector<IVect>& t_NodleAcous = var_t.GetNodleAcous();
    IVect& t_DirichletAcous = var_t.GetDirichletAcous();
    IVect& t_ProcVolMinus = var_t.GetProcVolMinus();
    IVect& t_ProcVolPlus = var_t.GetProcVolPlus();
    Vector<VectReal_wp>& t_JacobianAcous = var_t.GetJacobianAcous();
    Vector<IVect>& t_NodleVolMinus = var_t.GetNodleVolMinus();
    Vector<IVect>& t_NodleVolPlus = var_t.GetNodleVolPlus();
    if (rank_proc == 0)
      {
        int nb_faces;
        IVect NumGlobFace, ProcVolMinus, ProcVolPlus, NumDofAcous, DirAcous, NbDofAcous;
        VectReal_wp JacobAcous;
        for (int i = 0; i < t_NodleAcous.GetM(); i++)
          for (int j = 0; j < t_NodleAcous(i).GetM(); j++)
            nb_ddl_max = max(nb_ddl_max, t_NodleAcous(i)(j));
        
        nb_ddl_max++;
        Vector<bool> DofOnDirichlet(nb_ddl_max), DofUsed(nb_ddl_max);
        DofOnDirichlet.Fill(false); DirAcous.Reallocate(nb_ddl_max);
        for (int i = 0; i < t_DirichletAcous.GetM(); i++)
          DofOnDirichlet(t_DirichletAcous(i)) = true;
        
        for (int proc = 0; proc < nb_proc; proc++)
          if (proc_surf(proc))
            {
              if (proc == 0)
                {
                  NumGlobFace = mesh_plate.GlobElementNumber_Subdomain;
                  NumElem_Subdomain(proc) = NumGlobFace;
                  nb_faces = NumGlobFace.GetM();
                }
              else
                {
                  // on recoit les numeros globaux des faces
                  MPI_Recv(&nb_faces, 1, MPI_INTEGER, proc, 58, comm_vol, &status);
                  NumGlobFace.Reallocate(nb_faces);
                  MPI_Recv(NumGlobFace.GetData(), nb_faces, MPI_INTEGER, proc, 59, comm_vol, &status);
                }
              
              NumElem_Subdomain(proc) = NumGlobFace;

              // on renvoie les processeurs des elements qui s'appuient sur ces faces
              ProcVolMinus.Reallocate(nb_faces);
              ProcVolPlus.Reallocate(nb_faces);
              NbDofAcous.Reallocate(nb_faces);
              int nb_ddl_acous = 0;
              for (int i = 0; i < nb_faces; i++)
                {
                  int iglob = NumGlobFace(i);
                  ProcVolPlus(i) = t_ProcVolPlus(iglob);
                  ProcVolMinus(i) = t_ProcVolMinus(iglob);
                  NbDofAcous(i) = t_NodleAcous(iglob).GetM();
                  nb_ddl_acous += NbDofAcous(i);
                }
              
              if (proc == 0)
                {
                  LocProcVolMinus = ProcVolMinus;
                  LocProcVolPlus = ProcVolPlus;
                  LocNbDofAcous = NbDofAcous;
                }
              else
                {
                  MPI_Send(ProcVolMinus.GetData(), nb_faces, MPI_INTEGER, proc, 50, comm_vol);
                  MPI_Send(ProcVolPlus.GetData(), nb_faces, MPI_INTEGER, proc, 51, comm_vol);
                  MPI_Send(NbDofAcous.GetData(), nb_faces, MPI_INTEGER, proc, 52, comm_vol);
                }
              
              // ainsi que la numerotation NodleAcous
              // et les tableaux JacobianAcous et DirichletAcous
              DofUsed.Fill(false);
              NumDofAcous.Reallocate(nb_ddl_acous);
              JacobAcous.Reallocate(nb_ddl_acous);
              int nb = 0, nb_dir = 0;
              for (int i = 0; i < nb_faces; i++)
                {
                  int iglob = NumGlobFace(i);
                  int nb_dof = t_NodleAcous(iglob).GetM();
                  for (int j = 0; j < nb_dof; j++)
                    {
                      NumDofAcous(nb) = t_NodleAcous(iglob)(j);
                      JacobAcous(nb) = t_JacobianAcous(iglob)(j);
                      if (DofOnDirichlet(NumDofAcous(nb)))
                        {
                          if (!DofUsed(NumDofAcous(nb)))
                            {
                              DofUsed(NumDofAcous(nb)) = true;
                              DirAcous(nb_dir) = NumDofAcous(nb);
                              nb_dir++;
                            }
                        }
                      nb++;
                    }
                }
                  
              if (proc == 0)
                {
                  LocNumDofAcous = NumDofAcous;
                  LocJacobAcous = JacobAcous;
                  LocDirAcous.Reallocate(nb_dir);
                  for (int i = 0; i < nb_dir; i++)
                    LocDirAcous(i) = DirAcous(i);
                }
              else
                {
                  MPI_Send(NumDofAcous.GetData(), nb_ddl_acous, MPI_INTEGER, proc, 53, comm_vol);
                  Vector<int64_t> xtmp;
                  MpiSend(comm_vol, JacobAcous, xtmp, nb_ddl_acous, proc, 54);
                  MPI_Send(&nb_dir, 1, MPI_INTEGER, proc, 55, comm_vol);
                  MPI_Send(DirAcous.GetData(), nb_dir, MPI_INTEGER, proc, 56, comm_vol);
                }
            }
      }
    else
      {
        if (proc_surf(rank_proc))
          {
            int nb_elt = mesh_plate.GetNbElt();
            // we send to the master the global numbers of elements
            MPI_Send(&nb_elt, 1, MPI_INTEGER, 0, 58, comm_vol);
            MPI_Send(mesh_plate.GlobElementNumber_Subdomain.GetData(), nb_elt,
                      MPI_INTEGER, 0, 59, comm_vol);
            
            // we receive arrays related to plate variables
            LocProcVolMinus.Reallocate(nb_elt);
            LocProcVolPlus.Reallocate(nb_elt);
            LocNbDofAcous.Reallocate(nb_elt);
            MPI_Recv(LocProcVolMinus.GetData(), nb_elt, MPI_INTEGER,
                     0, 50, comm_vol, &status);
            MPI_Recv(LocProcVolPlus.GetData(), nb_elt, MPI_INTEGER,
                     0, 51, comm_vol, &status);
            MPI_Recv(LocNbDofAcous.GetData(), nb_elt, MPI_INTEGER,
                     0, 52, comm_vol, &status);
            
            int nb_ddl_acous = 0;
            for (int i = 0; i < nb_elt; i++)
              nb_ddl_acous += LocNbDofAcous(i);
            
            LocNumDofAcous.Reallocate(nb_ddl_acous);
            LocJacobAcous.Reallocate(nb_ddl_acous);
            int nb_dir = 0;
            MPI_Recv(LocNumDofAcous.GetData(), nb_ddl_acous, MPI_INTEGER,
                     0, 53, comm_vol, &status);
            
            Vector<int64_t> xtmp;
            MpiRecv(comm_vol, LocJacobAcous, xtmp, nb_ddl_acous, 0, 54, status);
            
            MPI_Recv(&nb_dir, 1, MPI_INTEGER, 0, 55, comm_vol, &status);

            LocDirAcous.Reallocate(nb_dir);
            MPI_Recv(LocDirAcous.GetData(), nb_dir, MPI_INTEGER, 0, 56, comm_vol, &status);
          }
      }

    if (proc_surf(rank_proc))
      {
        int nb_elt = mesh_plate.GetNbElt();
        
        // creating JacobianAcous, NodleAcous, NodleMeca and DirichletAcous
        int nb_ddl_acous = LocNumDofAcous.GetM();
        nb_ddl_max = 0;
        for (int i = 0; i < nb_ddl_acous; i++)
          nb_ddl_max = max(nb_ddl_max, LocNumDofAcous(i));
        
        nb_ddl_max++;
        IVect IndexDofAcous(nb_ddl_max);
        IndexDofAcous.Fill(-1);
        nb_ddl_acous = 0;
        for (int i = 0; i < LocNumDofAcous.GetM(); i++)
          {
            int j = LocNumDofAcous(i);
            if (IndexDofAcous(j) == -1)
              {
                IndexDofAcous(j) = nb_ddl_acous;
                nb_ddl_acous++;
              }
          }
        
        nodl_acous = nb_ddl_acous;
        NodleAcous.Reallocate(nb_elt);
        NodleMeca.Reallocate(nb_elt);
        JacobianAcous.Reallocate(nb_elt);
        DirichletAcous.Reallocate(LocDirAcous.GetM());
        for (int i = 0; i < LocDirAcous.GetM(); i++)
          DirichletAcous(i) = IndexDofAcous(LocDirAcous(i));
        
        int nb = 0;
        for (int i = 0; i < nb_elt; i++)
          {
            NodleMeca(i) = mesh_plate_num.Element(i).GetNodle();
            NodleAcous(i).Reallocate(LocNbDofAcous(i));
            JacobianAcous(i).Reallocate(LocNbDofAcous(i));
            for (int j = 0; j < NodleAcous(i).GetM(); j++)
              {
                NodleAcous(i)(j) = IndexDofAcous(LocNumDofAcous(nb));
                JacobianAcous(i)(j) = LocJacobAcous(nb);
                nb++;
              }
          }
        
        // then sorting dofs by processor number for exchange with volume variables
        Vector<bool> ProcUsed(nb_proc), DofUsed(nb_ddl_acous);
        ProcUsed.Fill(false);
        for (int i = 0; i < LocProcVolMinus.GetM(); i++)
          {
            ProcUsed(LocProcVolMinus(i)) = true;
            ProcUsed(LocProcVolPlus(i)) = true;
          }
        
        int nb_proc_used = 0;
        for (int i = 0; i < nb_proc; i++)
          if (ProcUsed(i))
            nb_proc_used++;
        
        NumberProc_DofSurf.Reallocate(nb_proc_used);
        DofSurfMinusOnProc.Reallocate(nb_proc_used);
        DofSurfPlusOnProc.Reallocate(nb_proc_used);
        int p = 0;
        for (int proc = 0; proc < nb_proc; proc++)
          if (ProcUsed(proc))
            {
              DofUsed.Fill(false);
              NumberProc_DofSurf(p) = proc;
              DofSurfMinusOnProc(p).Reallocate(nb_ddl_acous);
              nb = 0;
              for (int i = 0; i < LocProcVolMinus.GetM(); i++)
                if (LocProcVolMinus(i) == proc)
                  {
                    for (int j = 0; j < NodleAcous(i).GetM(); j++)
                      if (!DofUsed(NodleAcous(i)(j)))
                        {
                          DofUsed(NodleAcous(i)(j)) = true;
                          DofSurfMinusOnProc(p)(nb) = NodleAcous(i)(j);
                          nb++;
                        }
                  }
              
              if (nb == 0)
                DofSurfMinusOnProc(p).Clear();
              else
                DofSurfMinusOnProc(p).Resize(nb);
              
              DofUsed.Fill(false);
              nb = 0;
              DofSurfPlusOnProc(p).Reallocate(nb_ddl_acous);
              for (int i = 0; i < LocProcVolPlus.GetM(); i++)
                if (LocProcVolPlus(i) == proc)
                  {
                    for (int j = 0; j < NodleAcous(i).GetM(); j++)
                      if (!DofUsed(NodleAcous(i)(j)))
                        {
                          DofUsed(NodleAcous(i)(j)) = true;
                          DofSurfPlusOnProc(p)(nb) = NodleAcous(i)(j);
                          nb++;
                        }
                  }
              
              if (nb == 0)
                DofSurfPlusOnProc(p).Clear();
              else
                DofSurfPlusOnProc(p).Resize(nb);
              
              p++;
            }

      }

    // now treating volume dofs
    {
      IVect LocNumDofPlus, LocNumDofMinus;
      IVect LocNbDofMinus, LocNbDofPlus;
      if (rank_proc == 0)
        {
          // reconstitution of epart of the plate
          IVect Epart(t_NodleVolMinus.GetM());
          Epart.Fill(-1);
          
          for (int proc = 0; proc < nb_proc; proc++)
            for (int i = 0; i < NumElem_Subdomain(proc).GetM(); i++)
              Epart(NumElem_Subdomain(proc)(i)) = proc;
          
          IVect ProcVolMinus, ProcVolPlus, NumDofPlus, NumDofMinus, NbDofPlus, NbDofMinus;
          for (int proc = 0; proc < nb_proc; proc++)
            {            
              int nb_faces_minus = 0, nb_ddl_minus = 0;
              int nb_faces_plus = 0, nb_ddl_plus = 0;
              for (int i = 0; i < t_ProcVolMinus.GetM(); i++)
                if (t_ProcVolMinus(i) == proc)
                  {
                    nb_faces_minus++;
                    nb_ddl_minus += t_NodleVolMinus(i).GetM();
                  }
              
              for (int i = 0; i < t_ProcVolPlus.GetM(); i++)
                if (t_ProcVolPlus(i) == proc)
                  {
                    nb_faces_plus++;
                    nb_ddl_plus += t_NodleVolPlus(i).GetM();
                  }
              
              NbDofMinus.Reallocate(nb_faces_minus); NbDofPlus.Reallocate(nb_faces_plus);
              NumDofMinus.Reallocate(nb_ddl_minus); NumDofPlus.Reallocate(nb_ddl_plus);
              ProcVolMinus.Reallocate(nb_faces_minus); ProcVolPlus.Reallocate(nb_faces_plus);
              nb_faces_minus = 0; nb_faces_plus = 0;
              int nb = 0;
              for (int i = 0; i < t_ProcVolMinus.GetM(); i++)
                if (t_ProcVolMinus(i) == proc)
                  {
                    NbDofMinus(nb_faces_minus) = t_NodleVolMinus(i).GetM();
                    for (int j = 0; j < NbDofMinus(nb_faces_minus); j++)
                      {
                        NumDofMinus(nb) = t_NodleVolMinus(i)(j);
                        nb++;
                      }
                    
                    ProcVolMinus(nb_faces_minus) = Epart(i);
                    nb_faces_minus++;
                  }
              
              nb = 0;
              for (int i = 0; i < t_ProcVolPlus.GetM(); i++)
                if (t_ProcVolPlus(i) == proc)
                  {
                    NbDofPlus(nb_faces_plus) = t_NodleVolPlus(i).GetM();
                    for (int j = 0; j < NbDofPlus(nb_faces_plus); j++)
                      {
                        NumDofPlus(nb) = t_NodleVolPlus(i)(j);
                        nb++;
                      }
                    
                    ProcVolPlus(nb_faces_plus) = Epart(i);
                    nb_faces_plus++;
                  }
              
              // sending volume dofs and processor of each face
              if (proc == 0)
                {
                  LocNumDofMinus = NumDofMinus;
                  LocNumDofPlus = NumDofPlus;
                  LocNbDofPlus = NbDofPlus;
                  LocNbDofMinus = NbDofMinus;
                  LocProcVolMinus = ProcVolMinus;
                  LocProcVolPlus = ProcVolPlus;
                }
              else
                {
                  MPI_Send(&nb_faces_minus, 1, MPI_INTEGER, proc, 58, comm_vol);
                  MPI_Send(&nb_faces_plus, 1, MPI_INTEGER, proc, 59, comm_vol);
                  MPI_Send(NbDofMinus.GetData(), nb_faces_minus, MPI_INTEGER, proc, 50, comm_vol);
                  MPI_Send(NbDofPlus.GetData(), nb_faces_plus, MPI_INTEGER, proc, 51, comm_vol);
                  MPI_Send(ProcVolMinus.GetData(), nb_faces_minus, MPI_INTEGER, proc, 52, comm_vol);
                  MPI_Send(ProcVolPlus.GetData(), nb_faces_plus, MPI_INTEGER, proc, 53, comm_vol);
                  MPI_Send(NumDofMinus.GetData(), nb_ddl_minus, MPI_INTEGER, proc, 54, comm_vol);
                  MPI_Send(NumDofPlus.GetData(), nb_ddl_plus, MPI_INTEGER, proc, 55, comm_vol);
                }
            }
        }
      else
        {
          int nb_faces_minus = 0, nb_faces_plus = 0, nb_ddl_minus = 0, nb_ddl_plus = 0;
          MPI_Recv(&nb_faces_minus, 1, MPI_INTEGER, 0, 58, comm_vol, &status);
          MPI_Recv(&nb_faces_plus, 1, MPI_INTEGER, 0, 59, comm_vol, &status);
          LocNbDofMinus.Reallocate(nb_faces_minus); LocNbDofPlus.Reallocate(nb_faces_plus);
          MPI_Recv(LocNbDofMinus.GetData(), nb_faces_minus, MPI_INTEGER, 0, 50, comm_vol, &status);
          MPI_Recv(LocNbDofPlus.GetData(), nb_faces_plus, MPI_INTEGER, 0, 51, comm_vol, &status);
          for (int i = 0; i < LocNbDofMinus.GetM(); i++)
            nb_ddl_minus += LocNbDofMinus(i);
          
          for (int i = 0; i < LocNbDofPlus.GetM(); i++)
            nb_ddl_plus += LocNbDofPlus(i);
          
          //DISP(nb_faces_minus); DISP(nb_faces_plus); DISP(nb_ddl_minus); DISP(nb_ddl_plus);
          LocProcVolMinus.Reallocate(nb_faces_minus);
          LocProcVolPlus.Reallocate(nb_faces_plus);
          LocNumDofMinus.Reallocate(nb_ddl_minus);
          LocNumDofPlus.Reallocate(nb_ddl_plus);
          MPI_Recv(LocProcVolMinus.GetData(), nb_faces_minus, MPI_INTEGER,
                   0, 52, comm_vol, &status);
          MPI_Recv(LocProcVolPlus.GetData(), nb_faces_plus, MPI_INTEGER, 0, 53, comm_vol, &status);
          MPI_Recv(LocNumDofMinus.GetData(), nb_ddl_minus, MPI_INTEGER, 0, 54, comm_vol, &status);
          MPI_Recv(LocNumDofPlus.GetData(), nb_ddl_plus, MPI_INTEGER, 0, 55, comm_vol, &status);
        }
      
      //DISP(LocProcVolMinus.GetM());
      //DISP(LocNbDofMinus.GetM());
      
      IVect Glob_to_local(var_volume.var_harmonic.GetNbGlobalMeshDof());
      Glob_to_local.Fill(-1);
      for (int i = 0; i < mesh_num.GetNbDof(); i++)
        Glob_to_local(mesh_num.GlobDofNumber_Subdomain(i)) = i;
      
      Vector<bool> DofUsed(mesh_num.GetNbDof()); DofUsed.Fill(false);
      int nb_ddl_minus = 0, nb_ddl_plus = 0;
      for (int i = 0; i < LocNumDofMinus.GetM(); i++)
        {
          int n = Glob_to_local(LocNumDofMinus(i));
          if (!DofUsed(n))
            {
              DofUsed(n) = true;
              nb_ddl_minus++;
            }
        }
      
      ddl_vol_minus.Reallocate(nb_ddl_minus);
      nb_ddl_minus = 0; DofUsed.Fill(false);
      for (int i = 0; i < LocNumDofMinus.GetM(); i++)
        {
          int n = Glob_to_local(LocNumDofMinus(i));
          if (!DofUsed(n))
            {
              DofUsed(n) = true;
              ddl_vol_minus(nb_ddl_minus) = n;
              nb_ddl_minus++;
            }
        }
      
      DofUsed.Fill(false); nb_ddl_plus = 0;
      for (int i = 0; i < LocNumDofPlus.GetM(); i++)
        {
          int n = Glob_to_local(LocNumDofPlus(i));
          if (!DofUsed(n))
            {
              DofUsed(n) = true;
              nb_ddl_plus++;
            }
        }
      
      ddl_vol_plus.Reallocate(nb_ddl_plus);
      nb_ddl_plus = 0; DofUsed.Fill(false);
      for (int i = 0; i < LocNumDofPlus.GetM(); i++)
        {
          int n = Glob_to_local(LocNumDofPlus(i));
          if (!DofUsed(n))
            {
              DofUsed(n) = true;
              ddl_vol_plus(nb_ddl_plus) = n;
              nb_ddl_plus++;
            }
        }
      
      // index for dofs on Gamma+ and Gamma-
      IVect IndexDofPlus(mesh_num.GetNbDof()), IndexDofMinus(mesh_num.GetNbDof());
      IndexDofPlus.Fill(-1); IndexDofMinus.Fill(-1);
      for (int i = 0; i < ddl_vol_minus.GetM(); i++)
        IndexDofMinus(ddl_vol_minus(i)) = i;
      
      for (int i = 0; i < ddl_vol_plus.GetM(); i++)
        IndexDofPlus(ddl_vol_plus(i)) = i;
      
      // then sorting dofs by processor number for exchange with volume variables
      Vector<bool> ProcUsed(nb_proc);
      DofUsed.Reallocate(max(nb_ddl_minus, nb_ddl_plus));
      ProcUsed.Fill(false);
      for (int i = 0; i < LocProcVolMinus.GetM(); i++)
        ProcUsed(LocProcVolMinus(i)) = true;
      
      for (int i = 0; i < LocProcVolPlus.GetM(); i++)
        ProcUsed(LocProcVolPlus(i)) = true;
      
      int nb_proc_used = 0;
      for (int i = 0; i < nb_proc; i++)
        if (ProcUsed(i))
          nb_proc_used++;
      
      NumberProc_DofVol.Reallocate(nb_proc_used);
      DofVolMinusOnProc.Reallocate(nb_proc_used);
      DofVolPlusOnProc.Reallocate(nb_proc_used);
      int p = 0;
      for (int proc = 0; proc < nb_proc; proc++)
        if (ProcUsed(proc))
          {
            DofUsed.Fill(false);
            NumberProc_DofVol(p) = proc;
            DofVolMinusOnProc(p).Reallocate(nb_ddl_minus);
            int nb = 0; int nb_vol = 0;
            for (int i = 0; i < LocProcVolMinus.GetM(); i++)
              {
                if (LocProcVolMinus(i) == proc)
                  {
                    for (int j = 0; j < LocNbDofMinus(i); j++)
                      {
                        int num_dof = IndexDofMinus(Glob_to_local(LocNumDofMinus(nb_vol+j)));
                        if (!DofUsed(num_dof))
                          {
                            DofUsed(num_dof) = true;
                            DofVolMinusOnProc(p)(nb) = num_dof;
                            nb++;
                          }
                      }
                  }
                
                nb_vol += LocNbDofMinus(i);
              }
            
            if (nb == 0)
              DofVolMinusOnProc(p).Clear();
            else
              DofVolMinusOnProc(p).Resize(nb);
            
            DofUsed.Fill(false);
            DofVolPlusOnProc(p).Reallocate(nb_ddl_plus);
            nb = 0; nb_vol = 0;
            for (int i = 0; i < LocProcVolPlus.GetM(); i++)
              {
                if (LocProcVolPlus(i) == proc)
                  {
                    for (int j = 0; j < LocNbDofPlus(i); j++)
                      {
                        int num_dof = IndexDofPlus(Glob_to_local(LocNumDofPlus(nb_vol+j)));
                        if (!DofUsed(num_dof))
                          {
                            DofUsed(num_dof) = true;
                            DofVolPlusOnProc(p)(nb) = num_dof;
                            nb++;
                          }
                      }
                  }
                
                nb_vol += LocNbDofPlus(i);
              }
            
            if (nb == 0)
              DofVolPlusOnProc(p).Clear();
            else
              DofVolPlusOnProc(p).Resize(nb);
            
            p++;
          }
    }
    
    // creating MatchingDofAcous
    IVect& MatchingNumber_Plate = var_plate.var_harmonic.GetProcMatchingNeighbor();
    Vector<IVect>& MatchingDofOrig_Plate = var_plate.var_harmonic.GetOriginalMatchingDofNeighbor();
    
    if (mesh_plate_num.GetOrder() == mesh_num.GetOrder())
      {  
        IVect IndexDof(mesh_plate_num.GetNbDof());
        IndexDof.Fill(-1);
        for (int i = 0; i < NodleMeca.GetM(); i++)
          for (int j = 0; j < NodleMeca(i).GetM(); j++)
            IndexDof(NodleMeca(i)(j)) = NodleAcous(i)(j);
        
        int nb_domains =  MatchingNumber_Plate.GetM();
        NumberProcMatchingAcous = MatchingNumber_Plate;
        MatchingDofAcous.Reallocate(nb_domains);
        for (int i = 0; i < nb_domains; i++)
          {
            int nb_dof = MatchingDofOrig_Plate(i).GetM();
            MatchingDofAcous(i).Reallocate(nb_dof);
            for (int j = 0; j < nb_dof; j++)
              {
                int num_dof = MatchingDofOrig_Plate(i)(j);
                MatchingDofAcous(i)(j) = IndexDof(num_dof);
              }
          }
      }
    else
      {
        if (mesh_plate_num.GetOrder() == 1)
          {
            cout << "case not implemented" << endl;
            abort();
          }
        
        // creating a NodleMecaA adapted to this order
        int r = mesh_num.GetOrder();
        int rp = mesh_plate_num.GetOrder();
        int offset_int = mesh_plate.GetNbVertices() + (r-1)*mesh_plate.GetNbEdges();
        Vector<IVect> NodleMecaA(mesh_plate.GetNbElt());
        for (int i = 0; i < mesh_plate.GetNbElt(); i++)
          {
            int nb_dof = NodleAcous(i).GetM();
            NodleMecaA(i).Reallocate(nb_dof);
            // first vertices
            for (int j = 0; j < mesh_plate.Element(i).GetNbVertices(); j++)
              NodleMecaA(i)(j) = NodleMeca(i)(j);
            
            // then edges
            for (int j = 0; j < mesh_plate.Element(i).GetNbEdges(); j++)
              {
                int n1 = mesh_plate.Element(i).numVertex(j);
                int n2 = mesh_plate.Element(i).numVertex((j+1)%4);
                int n1_glob = mesh_plate.GlobVertexNumber_Subdomain(n1);
                int n2_glob = mesh_plate.GlobVertexNumber_Subdomain(n2);
                int ne = mesh_plate.Element(i).numEdge(j);
                int offset_loc = 4 + j*(r-1);
                int offset = mesh_plate.GetNbVertices() + ne*(r-1);
                
                if (n1_glob < n2_glob)
                  for (int k = 0; k < r-1; k++)
                    NodleMecaA(i)(offset_loc+k) = offset + k;
                else
                  for (int k = 0; k < r-1; k++)
                    NodleMecaA(i)(offset_loc+r-2-k) = offset+k;
              }
            
            // and interior
            for (int k = 4*r; k < nb_dof; k++)
              NodleMecaA(i)(k) = offset_int + k-4*r;
            
            offset_int += nb_dof-4*r;
          }
        
        // index array
        IVect IndexDof(offset_int);
        IndexDof.Fill(-1);
        for (int i = 0; i < NodleMecaA.GetM(); i++)
          for (int j = 0; j < NodleMecaA(i).GetM(); j++)
            IndexDof(NodleMecaA(i)(j)) = NodleAcous(i)(j);
        
        // then creating MatchingDofAcous
        int nb_domains =  MatchingNumber_Plate.GetM();
        NumberProcMatchingAcous = MatchingNumber_Plate;
        MatchingDofAcous.Reallocate(nb_domains);
        for (int i = 0; i < nb_domains; i++)
          {
            int nb_dof = MatchingDofOrig_Plate(i).GetM();
            int nb_vert = 0, nb_edges = 0;
            for (int j = 0; j < nb_dof; j++)
              {
                int num_dof = MatchingDofOrig_Plate(i)(j);
                if (num_dof < mesh_plate.GetNbVertices())
                  nb_vert++;
                else
                  nb_edges++;
              }
            
            nb_edges /= rp-1;
            
            int nb_dof_acous = nb_vert + (r-1)*nb_edges;
            MatchingDofAcous(i).Reallocate(nb_dof_acous);
            int j = 0, ja = 0;
            while (j < nb_dof)
              {
                int num_dof = MatchingDofOrig_Plate(i)(j);
                if (num_dof < mesh_plate.GetNbVertices())                  
                  {
                    MatchingDofAcous(i)(ja) = IndexDof(num_dof);
                    j++; ja++;
                  }
                else
                  {
                    int ne = (num_dof - mesh_plate.GetNbVertices())/(rp-1);
                    int offset = mesh_plate.GetNbVertices() + ne*(r-1);
                    for (int k = 0; k < r-1; k++)
                      MatchingDofAcous(i)(ja+k) = IndexDof(offset+k);
                    
                    ja += r-1; j += rp-1;
                  }
              }
          }
      }

    // displaying variables
    /*sleep(rank_proc);
    
    //DISP(mesh.GlobDofNumber_Subdomain);    
    DISP(ddl_vol_minus);
    DISP(ddl_vol_plus);
    //DISP(mesh_plate_num.GlobDofNumber_Subdomain);
    for (int i = 0; i < NodleAcous.GetM(); i++)
      {
        DISP(i);
        DISP(NodleAcous(i));
        DISP(NodleMeca(i));
        DISP(JacobianAcous(i));
      }
    
    DISP(DirichletAcous);
    
    for (int i = 0; i < NumberProc_DofVol.GetM(); i++)
      {
        DISP(i);
        DISP(NumberProc_DofVol(i));
        DISP(DofVolMinusOnProc(i));
        DISP(DofVolPlusOnProc(i));
      }
    
    for (int i = 0; i < NumberProc_DofSurf.GetM(); i++)
      {
        DISP(i);
        DISP(NumberProc_DofSurf(i));
        DISP(DofSurfMinusOnProc(i));
        DISP(DofSurfPlusOnProc(i));
      }    
    */
    // clearing variables of var_transmission
    var_t.Clear();
#endif
  }
  
  
  //! preparation of all variables needed for the simulation
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::ConstructAll(const string& data_file)
  {
    // several domains allowed for the plate
    var_plate.var_harmonic.InitIndices(100);
    var_plate.var_harmonic.SetTypeEquation("none");

    // only two domains for the volume
    var_volume.var_harmonic.InitIndices(100);
    var_volume.var_harmonic.SetTypeEquation("none");

    Vector<string> lines_data_file;

#ifdef SELDON_WITH_MPI
    MPI_Comm& comm_vol = var_volume.var_harmonic.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm_vol, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm_vol, &rank_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

#ifdef SELDON_WITH_MPI
    ReadLinesFile(data_file, lines_data_file, var_volume.var_harmonic.comm_group_mode);
#else
    ReadLinesFile(data_file, lines_data_file);
#endif

    // we read data file
    ReadInputFile(lines_data_file, *this);
    
    // if implicit-time scheme (currently not possible) is used for the plate, selection of solver
    var_plate.glob_solver = var_plate.var_harmonic.GetNewLinearSolver();
    ReadInputFile(lines_data_file, *var_plate.glob_solver);
    
    // for volume, we set rho = 1/(rho_f c^2) and mu = 1/rho_f
    TinyMatrix<Real_wp, Symmetric, 3, 3> Mu;
    Mu.SetIdentity();
    Mu *= 1.0/rho;
    for (int ref = 0; ref < var_volume.var_harmonic.ref_rho.GetM(); ref++)
      {
        var_volume.var_harmonic.ref_rho(ref).SetConstant(1.0/(rho*c2));
        var_volume.var_harmonic.ref_mu(ref).SetConstant(Mu);
      }
    
    // DISP(var_volume.var_harmonic.ref_rho(0));
    
    if (rank_proc == 0)
      cout << "Le fichier de donnees a ete lu " << endl;
    
    // mesh and finite elements are constructed
    this->ComputeMeshAndFiniteElement();
    
    if (nb_proc == 1)
      {
        var_plate.var_harmonic.mesh.Write("plate_vibro.mesh");
        var_volume.var_harmonic.mesh.Write("test_vibro.mesh");
      }
    
    /*--------Affichage des donnees physiques du probleme-----------------*/
    if (rank_proc == 0)
      {
        cout << "--Celerity: " << this->GetCelerity() << "m.s-1" << endl << endl;
        cout << "--Frequency: " << this->GetFrequency() << "Hz" <<endl << endl;
        cout << "--Time Step dt: "<< this->GetTimeStep() << endl << endl;
        cout << "Frequence supportee par le maillage :" 
             << (this->GetCelerity()*this->GetOrder())/(10*this->GetMeshSize())<<endl<<endl;
        cout << "---------------------------------------" << endl;
      }
    /*--------------------------------------------------------------------*/
    
    // we retrieve Dirichlet dofs if present
    var_volume.var_harmonic.TreatDirichletCondition();
    
    // right hand side is computed
    var_volume.ComputeRightHandSide();
    
    // computation of geometric quantities (F_i, DF_i, J_i, etc)
    var_volume.var_harmonic.ComputeMassMatrix();  
    
    // computation of mass and stiffness matrices
    var_volume.ComputeStiffnessMatrix();
    var_volume.ComputeMassMatrix();
    
    // 2-D stuff
    var_plate.var_harmonic.TreatDirichletCondition();
    var_plate.ComputeRightHandSide();
    var_plate.var_harmonic.ComputeMassMatrix();
    Vector<string> param(1); param(0) = string("LEAP_FROG");
    var_plate.var_time_scheme.SetTimeScheme(2, param);
    var_plate.ComputeStiffnessMatrix();
    var_plate.ComputeMassMatrix();
    var_plate.InitTimeIterations();
    this->m_ = GetNbMecanicDofSurface();
    this->n_ = GetNbMecanicDofSurface();

    // initializations before time iterations
    var_volume.InitTimeIterations();
    this->InitTimeIterations();
    
    /*
      computation of the global matrix
    if (nb_proc == 1)
      {
        int N = this->GetNbDof();
        Matrix<Real_wp, General, ArrayRowSparse> Ah(N, N);
        VectReal_wp Ones(N), AhOnes(N);
        AhOnes.Fill(0);
        for (int j = 0; j < N; j++)
          {
            Ones.Fill(0); Ones(j) = 1.0;
            this->EvaluateDerivativeFunction(0, 0, Ones, AhOnes, true, false);
            for (int i = 0; i < N; i++)
              if (abs(AhOnes(i)) > 1e-12)
                Ah.AddInteraction(i, j, AhOnes(i));
          }
          
        Ah.WriteText("AhOne.dat");
      }
    else
      {
#ifdef SELDON_WITH_MPI
        int N = this->GetNbDof();
        int Nvol_loc = this->var_volume.var_harmonic.GetNbDof();
        int Nplate_loc = this->var_plate.var_harmonic.GetNbDof();
        int Nvol = this->var_volume.var_harmonic.nodl_all;
        int Nplate = this->var_plate.var_harmonic.nodl_all;
        int Nall = this->var_volume.var_harmonic.nodl_all
        + 2*this->var_plate.var_harmonic.nodl_all;
        DISP(Nall);
        IVect Glob_to_local(Nall);
        Glob_to_local.Fill(-1);
        for (int i = 0; i < this->var_volume.var_harmonic
        .mesh_num.GlobDofNumber_Subdomain.GetM(); i++)
        Glob_to_local(this->var_volume.var_harmonic.mesh_num.GlobDofNumber_Subdomain(i)) = i;
        
        for (int i = 0; i < this->var_plate.var_harmonic
        .mesh_num.GlobDofNumber_Subdomain.GetM(); i++)
          {
            int num_dof = this->var_plate.var_harmonic.mesh_num.GlobDofNumber_Subdomain(i);
            Glob_to_local(Nvol + num_dof) = Nvol_loc + i;
            Glob_to_local(Nvol + Nplate + num_dof) = Nvol_loc + Nplate_loc + i;
          }
        
        Matrix<Real_wp, General, ArrayRowSparse> Ah(N, Nall);
        VectReal_wp Ones(N), AhOnes(N);
        AhOnes.Fill(0);
        for (int j = 0; j < Nall; j++)
          {
            int jloc = Glob_to_local(j);
            Ones.Fill(0);
            if (jloc >= 0)
              Ones(jloc) = 1.0;
            
            this->EvaluateDerivativeFunction(0, 0, Ones, AhOnes, true, false);
            for (int i = 0; i < N; i++)
              if (abs(AhOnes(i)) > 1e-12)
                Ah.AddInteraction(i, j, AhOnes(i));
            
          }
        
        DISP(this->var_volume.var_harmonic.mesh_num.GetNbDof());
        DISP(this->var_plate.var_harmonic.mesh_num.GetNbDof());
        this->var_volume.var_harmonic.mesh_num
        .GlobDofNumber_Subdomain.Write("num_vol"+to_str(rank_proc)+".dat");
        this->var_plate.var_harmonic.mesh_num
        .GlobDofNumber_Subdomain.Write("num_plate"+to_str(rank_proc)+".dat");
        Ah.WriteText("AhLoc" + to_str(rank_proc) + ".dat");
#endif
      }
    */
  }
  
  
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  ComputeDiagonalSchur(VectReal_wp& invMh, VectReal_wp& invMhMinusSurf,
		       VectReal_wp& invMhPlusSurf, bool assemble)
  {
    int nb_ddl_minus = ddl_vol_minus.GetM();
    int nb_ddl_plus = ddl_vol_plus.GetM();
    VectReal_wp MassVol = var_volume.Glob_mat_Dh->GetDiagonalDh();
    if (assemble)
      var_volume.Assemble(MassVol);
    
    int Nacous = this->GetNbAcousticDofSurface();
    invMhMinusSurf.Reallocate(Nacous);
    invMhPlusSurf.Reallocate(Nacous);
    VectReal_wp invMhMinus(nb_ddl_minus), invMhPlus(nb_ddl_plus);
    for (int i = 0; i < nb_ddl_minus; i++)
      invMhMinus(i) = 1.0 / MassVol(ddl_vol_minus(i));
    
    for (int i = 0; i < nb_ddl_plus; i++)
      invMhPlus(i) = 1.0 / MassVol(ddl_vol_plus(i));
    
    ExchangeVolumeToSurface(invMhMinus, invMhPlus, invMhMinusSurf, invMhPlusSurf);
      
    invMh.Reallocate(Nacous);
    for (int i = 0; i < Nacous; i++)
      invMh(i) = invMhMinusSurf(i) + invMhPlusSurf(i);    
  }
  

  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::InitTimeIterations()
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm_vol = var_volume.var_harmonic.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm_vol, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    // la matrice de masse de la mecanique
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(1.0);
    nat_mat.SetCoefStiffness(0.0);
    nat_mat.SetCoefDamping(0.0);
    
    FemMatrixFreeClass<Real_wp, typename TypeEqPlate::TypeEquationStationary>
      Kh(var_plate.var_harmonic);
    
    var_plate.var_harmonic.ComputeDiagonalMatrix(MassLambda, Kh, nat_mat);
    
    if (var_volume.GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_ORDER2)
      {
	VectReal_wp invMhMinusSurf, invMhPlusSurf;
	ComputeDiagonalSchur(invDiagonalSchur, invMhMinusSurf, invMhPlusSurf);
		
        if (nb_proc > 1)
          {
            OverlapDofsUz.Reallocate(var_plate.var_harmonic.GetNbOverlappedDof()/3);
            for (int i = 0; i < var_plate.var_harmonic.GetNbOverlappedDof()/3; i++)
              OverlapDofsUz(i) = var_plate.var_harmonic.GetOverlappedDofNumber(i);
          }
      }
  }
  
  
  //! returns the number of unknowns for the first-order evolution problem dU/dt = g(t, U)
  template<class TypeEqPlate, class TypeEqVol>
  int VibroAcousticSession<TypeEqPlate, TypeEqVol>::GetNbDof() const
  {
    int Nvol = this->GetNbDofVolume();
    int Nvec = this->var_volume.GetNbVectorialUnknowns();
    int Nplate = this->GetNbMecanicDofSurface();
    return Nvol + Nvec + 6*Nplate;
  }
  

  //! fills the volume variables Uminus Uplus from the surface variable u
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  ExchangeSurfaceToVolume(VectReal_wp& U, VectReal_wp& Uminus, VectReal_wp& Uplus)
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm_vol = var_volume.var_harmonic.comm_group_mode;
    int nb_process; MPI_Comm_size(comm_vol, &nb_process);
    int rank_proc; MPI_Comm_rank(comm_vol, &rank_proc);
#else
    int nb_process(1);
#endif

    if (nb_process == 1)
      {
        Uminus = U;
        Uplus = U;
        return;
      }

#ifdef SELDON_WITH_MPI
    VectReal_wp Uvol_minus_self, Uvol_plus_self;
    int nb_proc = NumberProc_DofSurf.GetM();
    Vector<MPI_Request> request_send(2*nb_proc);
    // sending u+ and u-
    Vector<bool> request_send_to_complete(2*nb_proc);
    request_send_to_complete.Fill(false);
    Vector<VectReal_wp> Vminus(nb_proc), Vplus(nb_proc);
    Vector<Vector<int64_t> > Vminus_tmp(nb_proc), Vplus_tmp(nb_proc);
    for (int p = 0; p < nb_proc; p++)
      {
        int proc = NumberProc_DofSurf(p);
        int nb_ddl_minus = DofSurfMinusOnProc(p).GetM();
        if (nb_ddl_minus > 0)
          {
            Vminus(p).Reallocate(nb_ddl_minus);
            for (int i = 0; i < nb_ddl_minus; i++)
              Vminus(p)(i) = U(DofSurfMinusOnProc(p)(i));
            
            if (proc == rank_proc)
              Uvol_minus_self = Vminus(p);
            else
              {
                request_send_to_complete(2*p) = true;                
                request_send(2*p) = MpiIsend(comm_vol, Vminus(p), Vminus_tmp(p),
                                             nb_ddl_minus, proc, 88);
              }
          }

        int nb_ddl_plus = DofSurfPlusOnProc(p).GetM();
        //cout << "Sending ... " << p << " " << proc << " " 
        // << nb_ddl_minus << " " << nb_ddl_plus << " " << rank_proc << endl;
        if (nb_ddl_plus > 0)
          {
            Vplus(p).Reallocate(nb_ddl_plus);
            for (int i = 0; i < nb_ddl_plus; i++)
              Vplus(p)(i) = U(DofSurfPlusOnProc(p)(i));
            
            if (proc == rank_proc)
              Uvol_plus_self = Vplus(p);
            else
              {
                //cout << "Sending to proc " << proc << " " << nb_ddl_plus << " dofs " << endl;
                request_send_to_complete(2*p+1) = true;
                request_send(2*p+1) = MpiIsend(comm_vol, Vplus(p), Vplus_tmp(p),
                                               nb_ddl_plus, proc, 89);
              }
          }
      }
    
    // then receiving u- and u+    
    nb_proc = NumberProc_DofVol.GetM();
    Vector<VectReal_wp> Uvol_minus(nb_proc), Uvol_plus(nb_proc);
    Vector<Vector<int64_t> > Uvol_minus_tmp(nb_proc), Uvol_plus_tmp(nb_proc);
    Vector<bool> request_recv_to_complete(2*nb_proc);
    request_recv_to_complete.Fill(false);
    Vector<MPI_Request> request_recv(2*nb_proc);    
    for (int p = 0; p < nb_proc; p++)
      {
        int proc = NumberProc_DofVol(p);
        int nb_ddl_minus = DofVolMinusOnProc(p).GetM();
        if (nb_ddl_minus > 0)
          {
            Uvol_minus(p).Reallocate(nb_ddl_minus);
            if (proc == rank_proc)
              Uvol_minus(p) = Uvol_minus_self;
            else
              {
                Uvol_minus(p).Fill(0);
                request_recv_to_complete(2*p) = true;
                request_recv(2*p) = MpiIrecv(comm_vol, Uvol_minus(p), Uvol_minus_tmp(p),
                                             nb_ddl_minus, proc, 88);
              }
          }

        int nb_ddl_plus = DofVolPlusOnProc(p).GetM();
        //cout << "Receiving ... " << p << " " << proc << " " << nb_ddl_minus 
        // << " " << nb_ddl_plus << " " << rank_proc << endl;
        if (nb_ddl_plus > 0)
          {
            Uvol_plus(p).Reallocate(nb_ddl_plus);
            
            if (proc == rank_proc)
              Uvol_plus(p) = Uvol_plus_self;
            else
              {
                Uvol_plus(p).Fill(0);
                //cout << "Receiving from proc " << proc << " " 
                // << nb_ddl_plus << " dofs " << endl;
                request_recv_to_complete(2*p+1) = true;
                request_recv(2*p+1) = MpiIrecv(comm_vol, Uvol_plus(p), Uvol_plus_tmp(p),
                                               nb_ddl_plus, proc, 89);
              }
          }
      }
    
    MPI_Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      if (request_send_to_complete(i))
        MPI_Wait(&request_send(i), &status);

    for (int i = 0; i < request_recv.GetM(); i++)
      if (request_recv_to_complete(i))
        MPI_Wait(&request_recv(i), &status);
    
    // then filling values
    Uplus.Fill(0); Uminus.Fill(0);
    for (int p = 0; p < nb_proc; p++)
      {
        int nb_ddl_minus = DofVolMinusOnProc(p).GetM();
        if ((nb_ddl_minus > 0) && (rank_proc != NumberProc_DofVol(p)))
          MpiCompleteIrecv(Uvol_minus(p), Uvol_minus_tmp(p), nb_ddl_minus);
        
        for (int i = 0; i < nb_ddl_minus; i++)
          Uminus(DofVolMinusOnProc(p)(i)) = Uvol_minus(p)(i);

        int nb_ddl_plus = DofVolPlusOnProc(p).GetM();
        if ((nb_ddl_plus > 0) && (rank_proc != NumberProc_DofVol(p)))
          MpiCompleteIrecv(Uvol_plus(p), Uvol_plus_tmp(p), nb_ddl_plus);
        
        for (int i = 0; i < nb_ddl_plus; i++)
          Uplus(DofVolPlusOnProc(p)(i)) = Uvol_plus(p)(i);
      }
#endif
  }
  
  
  //! fills surface variables Uminus_acous, Uplus_acous from volume variables Uminus and Uplus
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  ExchangeVolumeToSurface(VectReal_wp& Uminus, VectReal_wp& Uplus,
                          VectReal_wp& Uminus_acous, VectReal_wp& Uplus_acous)
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm_vol = var_volume.var_harmonic.comm_group_mode;
    int nb_process; MPI_Comm_size(comm_vol, &nb_process);
    int rank_proc; MPI_Comm_rank(comm_vol, &rank_proc);
#else
    int nb_process(1);
#endif

    if (nb_process == 1)
      {
        for (int i = 0; i < Uminus.GetM(); i++)
          Uminus_acous(i) = Uminus(i);

        for (int i = 0; i < Uplus.GetM(); i++)
          Uplus_acous(i) = Uplus(i);
        
        return;
      }
    
#ifdef SELDON_WITH_MPI
    VectReal_wp Uplus_self, Uminus_self;
    int nb_proc = NumberProc_DofVol.GetM();
    Vector<MPI_Request> request_send(2*nb_proc);
    // sending u+ and u-
    Vector<bool> request_send_to_complete(2*nb_proc);
    request_send_to_complete.Fill(false);
    Vector<VectReal_wp> Vplus(nb_proc), Vminus(nb_proc);
    Vector<Vector<int64_t> > Vplus_tmp(nb_proc), Vminus_tmp(nb_proc);
    for (int p = 0; p < nb_proc; p++)
      {
        int proc = NumberProc_DofVol(p);
        int nb_ddl_minus = DofVolMinusOnProc(p).GetM();
        if (nb_ddl_minus > 0)
          {
            Vminus(p).Reallocate(nb_ddl_minus);
            for (int i = 0; i < nb_ddl_minus; i++)
              Vminus(p)(i) = Uminus(DofVolMinusOnProc(p)(i));
            
            if (proc == rank_proc)
              Uminus_self = Vminus(p);
            else
              {
                request_send_to_complete(2*p) = true;
                request_send(2*p) = MpiIsend(comm_vol, Vminus(p), Vminus_tmp(p),
                                             nb_ddl_minus, proc, 80);
              }
          }

        int nb_ddl_plus = DofVolPlusOnProc(p).GetM();
        if (nb_ddl_plus > 0)
          {
            Vplus(p).Reallocate(nb_ddl_plus);
            for (int i = 0; i < nb_ddl_plus; i++)
              Vplus(p)(i) = Uplus(DofVolPlusOnProc(p)(i));
            
            if (proc == rank_proc)
              Uplus_self = Vplus(p);
            else
              {
                request_send_to_complete(2*p+1) = true;
                request_send(2*p+1) = MpiIsend(comm_vol, Vplus(p), Vplus_tmp(p),
                                               nb_ddl_plus, proc, 81);
              }
          }
      }
    
    // then receiving u- and u+    
    nb_proc = NumberProc_DofSurf.GetM();
    Vector<VectReal_wp> Uvol_minus(nb_proc), Uvol_plus(nb_proc);
    Vector<Vector<int64_t> > Uvol_minus_tmp(nb_proc), Uvol_plus_tmp(nb_proc);
    Vector<bool> request_recv_to_complete(2*nb_proc);
    request_recv_to_complete.Fill(false);
    Vector<MPI_Request> request_recv(2*nb_proc);    
    for (int p = 0; p < nb_proc; p++)
      {
        int proc = NumberProc_DofSurf(p);
        int nb_ddl_minus = DofSurfMinusOnProc(p).GetM();
        if (nb_ddl_minus > 0)
          {
            Uvol_minus(p).Reallocate(nb_ddl_minus);
            if (proc == rank_proc)
              Uvol_minus(p) = Uminus_self;
            else
              {
                Uvol_minus(p).Fill(0);
                request_recv_to_complete(2*p) = true;
                request_recv(2*p) = MpiIrecv(comm_vol, Uvol_minus(p), Uvol_minus_tmp(p),
                                             nb_ddl_minus, proc, 80);
              }
          }

        int nb_ddl_plus = DofSurfPlusOnProc(p).GetM();
        if (nb_ddl_plus > 0)
          {
            Uvol_plus(p).Reallocate(nb_ddl_plus);
            
            if (proc == rank_proc)
              Uvol_plus(p) = Uplus_self;
            else
              {
                Uvol_plus(p).Fill(0);
                request_recv_to_complete(2*p+1) = true;
                request_recv(2*p+1) = MpiIrecv(comm_vol, Uvol_plus(p), Uvol_minus_tmp(p),
                                               nb_ddl_plus, proc, 81);
              }
          }
      }
    
    MPI_Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      if (request_send_to_complete(i))
        MPI_Wait(&request_send(i), &status);

    for (int i = 0; i < request_recv.GetM(); i++)
      if (request_recv_to_complete(i))
        MPI_Wait(&request_recv(i), &status);
    
    // then filling Uplus_acous and Uminus_acous
    for (int p = 0; p < nb_proc; p++)
      {
        int nb_ddl_minus = DofSurfMinusOnProc(p).GetM();
        if ((nb_ddl_minus > 0) && (rank_proc != NumberProc_DofSurf(p)))
          MpiCompleteIrecv(Uvol_minus(p), Uvol_minus_tmp(p), nb_ddl_minus);
        
        for (int i = 0; i < nb_ddl_minus; i++)
          Uminus_acous(DofSurfMinusOnProc(p)(i)) = Uvol_minus(p)(i);
        
        int nb_ddl_plus = DofSurfPlusOnProc(p).GetM();
        if ((nb_ddl_plus > 0) && (rank_proc != NumberProc_DofSurf(p)))
          MpiCompleteIrecv(Uvol_plus(p), Uvol_plus_tmp(p), nb_ddl_plus);
        
        for (int i = 0; i < nb_ddl_plus; i++)
          Uplus_acous(DofSurfPlusOnProc(p)(i)) = Uvol_plus(p)(i);
      }
#endif
    
  }

  
  //! assembles the vector X (acoustic unknowns)
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  AssembleAcoustic(VectReal_wp& X) const
  {
#ifdef SELDON_WITH_MPI
    const MPI_Comm& comm_surf = var_plate.var_harmonic.comm_group_mode;
    int nb_process; MPI_Comm_size(comm_surf, &nb_process);
    
    glob_chrono.Start(VirtualTimer::COMM);
    
    // if only one processor, no communication needed, we leave method
    if (nb_process <= 1)
      return;
    
    int nb_dom = this->MatchingDofAcous.GetM();
    Vector<MPI_Request> request_send(nb_dom), request_recv(nb_dom);
    MPI_Status status;
    Vector<VectReal_wp > xsend(nb_dom);
    Vector<Vector<int64_t> > xsend_tmp(nb_dom);
    int tag = 85;
    
    // sending informations to other domains
    for (int i = 0; i < nb_dom; i++)
      {
	int j = NumberProcMatchingAcous(i);
	int nb = this->MatchingDofAcous(i).GetM();
	if (nb > 0)
          {
            xsend(i).Reallocate(nb);
            for (int k = 0; k < nb; k++)
              xsend(i)(k) = X(this->MatchingDofAcous(i)(k));
            
            // sending the value to the corresponding processor
            request_send(i) = MpiIsend(comm_surf, xsend(i), xsend_tmp(i), nb, j, tag);
          }
      }
    
    // receiving the informations
    Vector<VectReal_wp> xdom(nb_dom);
    Vector<Vector<int64_t> > xdom_tmp(nb_dom);
    for (int i = 0; i < nb_dom; i++)
      {
	int j = NumberProcMatchingAcous(i);
	int nb = this->MatchingDofAcous(i).GetM();
	if (nb > 0)
          {
            xdom(i).Reallocate(nb);
            xdom(i).Fill(0);
            
            // receiving the values of domain j
            request_recv(i) = MpiIrecv(comm_surf, xdom(i), xdom_tmp(i), nb, j, tag);
          }
      }
    
    // now waiting all communications are effective
    for (int i = 0; i < nb_dom; i++)
      if (this->MatchingDofAcous(i).GetM() > 0)
        MPI_Wait(&request_recv(i), &status);
    
    for (int i = 0; i < nb_dom; i++)
      if (this->MatchingDofAcous(i).GetM() > 0)
        MPI_Wait(&request_send(i), &status);
    
    for (int i = 0; i < nb_dom; i++)
      {
        if (this->MatchingDofAcous(i).GetM() > 0)
          MpiCompleteIrecv(xdom(i), xdom_tmp(i), xdom(i).GetM());
        
	// summing the values
        for (int k = 0; k < this->MatchingDofAcous(i).GetM(); k++)
          X(this->MatchingDofAcous(i)(k)) += xdom(i)(k);
      }
    
    glob_chrono.Stop(VirtualTimer::COMM);
#endif
  }

  
  //! Evaluation of the function g in the first-order formulation du/dt = g(t, u)
  /*!
    The following pde are considered :
    1/(rho_f c^2)  dp/dt + div v = 0
    rho_f  dv/dt + grad p = 0
    dtheta/dt = thetaPoint
    du/dt = uPoint
    rho delta^3/12 d thetaPoint/dt
    - delta^3/12 Div(C epsilon(theta)) + delta mu (grad u + theta) = 0
    rho delta du/dt - delta div(mu (grad u + theta))
   */
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                             const VectReal_wp& X, VectReal_wp& Y, bool invert, bool source)
  {
    int Nvol = this->GetNbDofVolume();
    int Nvec = this->var_volume.GetNbVectorialUnknowns();
    int Nplate = this->GetNbMecanicDofSurface();
    int Nacous = this->GetNbAcousticDofSurface();
    
    VectReal_wp P, V, Lambda, dLambda, Uz, dUz;
    VectReal_wp Prod_P, Prod_V, Prod_Lambda, Prod_dLambda, Prod_Uz, Prod_dUz;
    
    Real_wp* xptr = X.GetData();
    Real_wp* yptr = Y.GetData();
    
    P.SetData(Nvol, xptr); xptr += Nvol; 
    Prod_P.SetData(Nvol, yptr); yptr += Nvol; 

    V.SetData(Nvec, xptr); xptr += Nvec; 
    Prod_V.SetData(Nvec, yptr); yptr += Nvec; 
    
    Lambda.SetData(3*Nplate, xptr); xptr += 2*Nplate; 
    Prod_Lambda.SetData(3*Nplate, yptr); yptr += 2*Nplate; 

    Uz.SetData(Nplate, xptr); xptr += Nplate; 
    Prod_Uz.SetData(Nplate, yptr); yptr += Nplate;     
    
    dLambda.SetData(3*Nplate, xptr); xptr += 2*Nplate; 
    Prod_dLambda.SetData(3*Nplate, yptr); yptr += 2*Nplate; 
    
    dUz.SetData(Nplate, xptr); xptr += Nplate; 
    Prod_dUz.SetData(Nplate, yptr); yptr += Nplate; 
    
    // mu dP/dt + div V = f
    // => Mh dP/dt + RhS V + Sh P + Ch* dUz/dt = 0
    var_volume.ApplyOperatorRhScalar(-1.0, tn, V, 0, Prod_P);
    var_volume.ApplyOperatorSh(-1.0, tn, P, 1.0, Prod_P);
    if (source)
      this->var_volume.AddScalarSourceAtTime(1.0, tn, nb_deriv, Prod_P); 
    
    // rho dV/dt + grad P = 0
    // =>  Bh dV/dt + RhV P + ShV V = 0
    var_volume.ApplyOperatorRhVectorial(-1.0, tn, P, 0, Prod_V);
    var_volume.ApplyOperatorShVectorial(-1.0, tn, V, 1.0, Prod_V);    
    
    VectReal_wp Ch_dUz(Nacous);
    Ch_dUz.Fill(0);
    this->MltAddChTranspose(1.0, dUz, Ch_dUz);
    this->AssembleAcoustic(Ch_dUz);
    
    int nb_ddl_minus = this->ddl_vol_minus.GetM();
    int nb_ddl_plus = this->ddl_vol_plus.GetM();
    VectReal_wp Ch_dUz_minus(nb_ddl_minus), Ch_dUz_plus(nb_ddl_plus);
    this->ExchangeSurfaceToVolume(Ch_dUz, Ch_dUz_minus, Ch_dUz_plus);

    this->var_volume.Assemble(Prod_P);    
    for (int i = 0; i < nb_ddl_minus; i++)
      {
        int dof_m = this->ddl_vol_minus(i);
        Prod_P(dof_m) -= Ch_dUz_minus(i);
      }
    
    for (int i = 0; i < nb_ddl_plus; i++)
      {
        int dof_p = this->ddl_vol_plus(i);
        Prod_P(dof_p) += Ch_dUz_plus(i);
      }
    
    // dLambda/dt = LambdaPoint
    for (int i = 0; i < 3*Nplate; i++)
      Prod_Lambda(i) = dLambda(i);
    
    // Mlambda dLambdaPoint/dt - Kh Lambda - Ch P = 0
    this->var_plate.ApplyOperatorKh(1.0, tn, Lambda, 0.0, Prod_dLambda);
    if (source)
      this->var_plate.AddScalarSourceAtTime(1.0, tn, nb_deriv+1, Prod_dLambda); 

    VectReal_wp Pminus(nb_ddl_minus), Pplus(nb_ddl_plus),
      Pvol_minus(Nacous), Pvol_plus(Nacous), SautP(Nacous);
    for (int i = 0; i < nb_ddl_minus; i++)
      Pminus(i) = P(this->ddl_vol_minus(i));

    for (int i = 0; i < nb_ddl_plus; i++)
      Pplus(i) = P(this->ddl_vol_plus(i));
    
    this->ExchangeVolumeToSurface(Pminus, Pplus, Pvol_minus, Pvol_plus);
    for (int i = 0; i < Nacous; i++)
      SautP(i) = Pvol_minus(i) - Pvol_plus(i);
    
    this->MltAddCh(1.0, SautP, Prod_dUz);
    this->var_plate.Assemble(Prod_dLambda);
    
    const VectReal_wp& MassLambda = this->GetMecanicMass();
    
    // assembling before inversion
    if (invert)
      {
        this->var_volume.SolveOperatorDh(Prod_P);
        this->var_volume.SolveOperatorBh(Prod_V);
        
        for (int i = 0; i < 3*Nplate; i++)
          Prod_dLambda(i) /= MassLambda(i);
      }

    P.Nullify(); V.Nullify(); Lambda.Nullify(); dLambda.Nullify();
    Prod_P.Nullify(); Prod_V.Nullify(); Prod_Lambda.Nullify(); Prod_dLambda.Nullify();
    Uz.Nullify(); dUz.Nullify(); Prod_Uz.Nullify(); Prod_dUz.Nullify();
  }
  

  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  EvaluateFunction(const Real_wp& tn, const VectReal_wp& X, VectReal_wp& Y,
                   bool invert_mass, bool source)
  {
    EvaluateDerivativeFunction(tn, 0, X, Y, invert_mass, source);
  }


  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  GiveIterate(int n, const Real_wp& t, VectReal_wp& Y)
  {
    int Nvol = this->GetNbDofVolume();
    int Nplate = this->GetNbMecanicDofSurface();
    int Nvec = this->var_volume.GetNbVectorialUnknowns();
    
    VectReal_wp P, V, Lambda, dLambda, dY(Y);
    
    Real_wp* xptr = Y.GetData();
    
    P.SetData(Nvol, xptr); xptr += Nvol; 
    V.SetData(Nvec, xptr); xptr += Nvec;
    Lambda.SetData(3*Nplate, xptr); xptr += 3*Nplate; 
    dLambda.SetData(3*Nplate, xptr); xptr += 3*Nplate; 
    
    EvaluateFunction(t, Y, dY, true, false);

    VectReal_wp dP, dV, dLambda_bis, d2Lambda_dt2;
    
    Real_wp* yptr = dY.GetData();
    
    dP.SetData(Nvol, yptr); yptr += Nvol; 
    dV.SetData(Nvec, yptr); yptr += Nvec;
    dLambda_bis.SetData(3*Nplate, yptr); yptr += 3*Nplate; 
    d2Lambda_dt2.SetData(3*Nplate, yptr); yptr += 3*Nplate; 
    
    this->var_volume.WriteVectorialSnapshot(n, t, V);
    this->var_volume.WriteSnapshot(n, t, P);    
    this->var_plate.WriteSnapshot(n, t, dLambda);
    
    // on appelle le calcul du rayonnement
    if (var_volume.var_harmonic.output_rcs_param.GetNbPointsOutside() > 0)
      {
        // on appelle le calcul du champ en des points exterieurs au domaine
        var_volume.var_harmonic.output_rcs_param.
          WriteOutput(n, t, this->var_volume.GetTimeStep(), P, dP);
      }

    P.Nullify(); V.Nullify(); Lambda.Nullify(); dLambda.Nullify();
    dP.Nullify(); dV.Nullify(); dLambda_bis.Nullify(); d2Lambda_dt2.Nullify();
  }


  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  GiveNumberIterations(const Real_wp&, int n) {}
  
  
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  GiveFinalIterate(int n, const Real_wp& t, VectReal_wp& Y)
  {
    GiveIterate(n, t, Y);
  }
  

  template<class TypeEqPlate, class TypeEqVol>
  template<class MatrixSparse>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  ComputeCbarre(MatrixSparse& Ch)
  {
    int nFEM = var_plate.var_harmonic.GetMeshNumbering(0).GetNbDof();
    int nAcous = nodl_acous;
    Ch.Reallocate(nFEM, nAcous);
    
    // boucle sur les faces
    for (int n = 0; n < var_plate.var_harmonic.mesh.GetNbElt(); n++)
      {        
        int nb_dof_acous = NodleAcous(n).GetM();
        int nb_dof_meca = NodleMeca(n).GetM();
        VectReal_wp xloc(nb_dof_acous), yloc(nb_dof_meca);
        
        for (int i = 0; i < nb_dof_acous; i++)
          {
            xloc.Fill(0);
            xloc(i) = JacobianAcous(n)(i);
            
            int type_elt = var_plate.var_harmonic.mesh.GetTypeElement(n);
            
            proj_meca.TransposeProjectScalar(xloc, yloc, type_elt);
            
            for (int j = 0; j < nb_dof_meca; j++)
              {
                // boucle sur les ddl locaux
                int i_glob = NodleAcous(n)(i);
                int j_glob = NodleMeca(n)(j);
                Ch.AddInteraction(j_glob, i_glob, yloc(j));
	      }
          }	
      }
  }
    
  
  //! computation of Y = Y + alpha*Ch*X
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  MltAddCh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y) const
  {
    const Mesh<Dimension2>& mesh = var_plate.var_harmonic.mesh;

    VectReal_wp xloc, yloc;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        const IVect& NodleA = NodleAcous(i);
        int nb_dof_acous = NodleA.GetM();
        int type_elt = mesh.GetTypeElement(i);
        xloc.Reallocate(nb_dof_acous);
        for (int j = 0; j < nb_dof_acous; j++)
          xloc(j) = X(NodleA(j))*JacobianAcous(i)(j);
        
        const IVect& NodleM = NodleMeca(i);
        int nb_dof_meca = NodleM.GetM();
        yloc.Reallocate(nb_dof_meca); yloc.Fill(0);
        proj_meca.TransposeProjectScalar(xloc, yloc, type_elt);
        
        for (int j = 0; j < nb_dof_meca; j++)
          Y(NodleM(j)) += alpha*yloc(j);
      }
  }
  
    
  //! computation of Y = Y + alpha*Ch^T*X
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  MltAddChTranspose(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y) const
  {
    const Mesh<Dimension2>& mesh = var_plate.var_harmonic.mesh;

    VectReal_wp xloc, yloc;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        int type_elt = mesh.GetTypeElement(i);
        const IVect& NodleM = NodleMeca(i);
        int nb_dof_meca = NodleM.GetM();
        xloc.Reallocate(nb_dof_meca);
        for (int j = 0; j < nb_dof_meca; j++)
          xloc(j) = X(NodleM(j));
        
        const IVect& NodleA = NodleAcous(i);
        int nb_dof_acous = NodleA.GetM();
        yloc.Reallocate(nb_dof_meca); yloc.Fill(0);
        proj_meca.ProjectScalar(xloc, yloc, type_elt);
        
        for (int j = 0; j < nb_dof_acous; j++)
          Y(NodleA(j)) += alpha*yloc(j)*JacobianAcous(i)(j);
      }
    
    for (int i = 0; i < DirichletAcous.GetM(); i++)
      Y(DirichletAcous(i)) = 0;
    
  }
  
  
  //! computation of the Yvec = beta*Yvec + alpha*A*X 
  //! where A is the matrix to invert on the surface
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::MltAddVector(const Real_wp& alpha,
		 const VectReal_wp& X, const Real_wp& beta, VectReal_wp& Yvec) const
  {
    if (beta == Real_wp(0))
      Yvec.Fill(0);
    else
      Mlt(beta, Yvec);
    
    int Nplate = this->GetNbMecanicDofSurface();
    int Nacous = this->GetNbAcousticDofSurface();
    
    VectReal_wp Y(Nacous);
    VectReal_wp Y2(Nplate);
        
    Y.Zero();    
    this->MltAddChTranspose(1.0, X, Y);
    
    this->AssembleAcoustic(Y);
    
    for (int i = 0; i < Y.GetM(); i++)
      Y(i) *= this->invDiagonalSchur(i);
    
    Y2.Zero();
    this->MltAddCh(1.0, Y, Y2);
    this->var_plate.Assemble(Y2);
    
    Real_wp coef_stiff = 0.25*square(this->var_plate.GetTimeStep());
    
    for (int k = 0; k < Yvec.GetM(); k++)
      Yvec(k) += alpha*(this->MassLambda(2*Nplate + k)*X(k) + coef_stiff*Y2(k));
    
  }


  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::MltAddVector(const Complex_wp& alpha,
		 const VectComplex_wp& X, const Complex_wp& beta, VectComplex_wp& Yvec) const
  {
    cout << "not implemented" << endl;
    abort();
  }


  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::MltAddVector(const Real_wp& alpha, const SeldonTranspose&,
		 const VectReal_wp& X, const Real_wp& beta, VectReal_wp& Yvec) const
  {
    cout << "not implemented" << endl;
    abort();
  }


  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::MltAddVector(const Complex_wp& alpha, const SeldonTranspose&,
		 const VectComplex_wp& X, const Complex_wp& beta, VectComplex_wp& Yvec) const
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! computation of Yvec = A*X
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::MltVector(const VectReal_wp& X, VectReal_wp& Yvec) const
  {
    Yvec.Zero();
    this->MltAddVector(1.0, X, 0.0, Yvec);
  }


  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::MltVector(const VectComplex_wp& X, VectComplex_wp& Yvec) const
  {
    cout << "not implemented" << endl;
    abort();
  }


  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::MltVector(const SeldonTranspose&, const VectReal_wp& X, VectReal_wp& Yvec) const
  {
    cout << "not implemented" << endl;
    abort();
  }


  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::MltVector(const SeldonTranspose&, const VectComplex_wp& X, VectComplex_wp& Yvec) const
  {
    cout << "not implemented" << endl;
    abort();
  }  


  //! solving linear system A X = B
  template<class TypeEqPlate, class TypeEqVol>
  void VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  SolveSchur(VectReal_wp & X, VectReal_wp & B) 
  {
    int nb_iterations = 100;
    Real_wp stopping_criterion = 1e-12;
    Iteration<Real_wp> iter(nb_iterations, stopping_criterion);
    iter.SetInitGuess(false);
    Preconditioner_Base<Real_wp> pr;
    iter.HideMessages();
    //iter.ShowFullHistory();
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm_vol = var_volume.var_harmonic.comm_group_mode;
    int nb_process; MPI_Comm_size(comm_vol, &nb_process);

    if (nb_process > 1)
      {
        DistributedVector<Real_wp> source(OverlapDofsUz, var_plate.var_harmonic.comm_group_mode);
        DistributedVector<Real_wp> solution(OverlapDofsUz, var_plate.var_harmonic.comm_group_mode);
        
        source.SetData(B);
        solution.SetData(X);
        
        Cg(*this, solution, source, pr, iter);
        
        source.Nullify(); solution.Nullify();
      }
    else
#endif
      Cg(*this, X, B, pr, iter);
    
    nb_iterations = iter.GetNumberIteration();
  }
  
  
  VibroAcousticLeapFrogScheme::VibroAcousticLeapFrogScheme()
  {
    first_order = false;
    dt = 0.0;
  }
  
  
  //! initialization of arrays necessary for the leap frog scheme
  template<class GenericPb>
  void VibroAcousticLeapFrogScheme::
  SetInitialCondition(const Real_wp& t0, const Real_wp& dt_, GenericPb& var)
  {
    dt = dt_;
    
    // variables for the plate
    int Nacous = var.GetNbAcousticDofSurface();
    int Nplate = var.GetNbMecanicDofSurface();
    int Nvol = var.GetNbDofVolume();
    int Nvec = var.var_volume.GetNbVectorialUnknowns();
    
    if (first_order)
      {
        Vh_n.Reallocate(Nvec);
        Vh_np1.Reallocate(Nvec);
        Vh_n.Fill(0);
        Vh_np1.Fill(0);
      }
    
    Lambda_np1.Reallocate(3*Nplate); Lambda_np1.Fill(0);
    Lambda_n.Reallocate(3*Nplate); Lambda_n.Fill(0);
    Lambda_nm1.Reallocate(3*Nplate); Lambda_nm1.Fill(0);
    
    KlambdaLn.Reallocate(3*Nplate); KlambdaLn.Fill(0);
    Fuz.Reallocate(Nplate); Fuz.Fill(0);
    U_np1.Reallocate(Nplate); U_np1.Fill(0);
    
    Ph_n.Reallocate(Nvol); Ph_n.Fill(0);
    Ph_nm1.Reallocate(Nvol); Ph_nm1.Fill(0);
    
    KhPn.Reallocate(Nvol); KhPn.Fill(0);
    
    Gvol_minus.Reallocate(var.ddl_vol_minus.GetM()); Gvol_minus.Fill(0);
    Gvol_plus.Reallocate(var.ddl_vol_plus.GetM()); Gvol_plus.Fill(0);
    ChLambdaMinus_n.Reallocate(var.ddl_vol_minus.GetM()); ChLambdaMinus_n.Fill(0);
    ChLambdaPlus_n.Reallocate(var.ddl_vol_plus.GetM()); ChLambdaPlus_n.Fill(0);
    ChLambdaMinus_nm1.Reallocate(var.ddl_vol_minus.GetM()); ChLambdaMinus_nm1.Fill(0);
    ChLambdaPlus_nm1.Reallocate(var.ddl_vol_plus.GetM()); ChLambdaPlus_nm1.Fill(0);
    
    Fvol.Reallocate(Nacous); Fvol.Fill(0);
    Gsurf_minus.Reallocate(Nacous); Gsurf_minus.Fill(0);
    Gsurf_plus.Reallocate(Nacous); Gsurf_plus.Fill(0);
    
    /*
    
    Ph_n.FillRand(); Mlt(1e-9, Ph_n);
    Ph_nm1.FillRand(); Mlt(1e-9, Ph_nm1);
    Vh_n.FillRand(); Mlt(1e-9, Vh_n);
    Lambda_n.FillRand(); Mlt(1e-9, Lambda_n);
    Lambda_np1.FillRand(); Mlt(1e-9, Lambda_np1);
    
    Matrix<Real_wp, General, ArrayRowSparse> Ch;
    var.ComputeCbarre(Ch);
    Ch.WriteText("Ch.dat");
    
    var.var_plate.Glob_mat_Kh.mat_iterative_unsym.WriteText("Kh.dat");

    var.var_volume.Glob_mat_RhS.WriteText("RhS.dat");
    var.var_volume.Glob_mat_RhV.WriteText("RhV.dat");

    var.ddl_vol_minus.Write("DdlVolMinus.dat");
    var.ddl_vol_plus.Write("DdlVolPlus.dat");

    var.var_plate.var_harmonic.ImposeNullDirichletCondition(Lambda_n);
    var.var_plate.var_harmonic.ImposeNullDirichletCondition(Lambda_nm1);        
    */
  }

  
  //! we compute solution for the next time step
  template<class GenericPb>
  void VibroAcousticLeapFrogScheme::Advance(const Real_wp& t, int n, GenericPb& var)
  {
    // nombre de ddls
    int Nplate = var.GetNbMecanicDofSurface();
    int Nacous = var.GetNbAcousticDofSurface();
    Real_wp dt2 = dt*dt;
    
    const VectReal_wp& MassLambda = var.GetMecanicMass();
    const VectReal_wp& Mh = var.var_volume.Glob_mat_Dh->GetDiagonalDh();
    const VectReal_wp& Mh_minus_dt_Sh = var.var_volume.Glob_mat_Dh->GetDiagonalDhMinusdtSh();
    const VectReal_wp& invMh_plus_dt_Sh = var.var_volume.Glob_mat_Dh->GetInverseDiagonalDhPlusdtSh();
    
    if (first_order)
      {
        // schema completement explicite
        // on evalue d'abord Lambda^n+1
        // M^L (Lambda^n+1 - 2 Lambda^n + Lambda^n-1) + K^L Lambda^n - Ch P^n = 0
        
        // calcul de K^L Lambda^n - F^n
        var.var_plate.ApplyOperatorKh(-1.0, t, Lambda_n, 0.0, KlambdaLn);
        if (t < var.var_plate.GetFinalTimeSource())
          var.var_plate.AddScalarSourceAtTime(-1.0, t, 1, KlambdaLn);
        
        // on rajoute - Ch P^n
        for (int i = 0; i < var.ddl_vol_minus.GetM(); i++)
          Gvol_minus(i) = Ph_n(var.ddl_vol_minus(i));
        
        for (int i = 0; i < var.ddl_vol_plus.GetM(); i++)
          Gvol_plus(i) = Ph_n(var.ddl_vol_plus(i));
    
        var.ExchangeVolumeToSurface(Gvol_minus, Gvol_plus, Gsurf_minus, Gsurf_plus);
        for (int i = 0; i < Nacous; i++)
          Fvol(i) = Gsurf_minus(i) - Gsurf_plus(i);
        
        Fuz.Fill(0);
        var.MltAddCh(1.0, Fvol, Fuz);

        for (int i = 0; i < Nplate; i++)
          KlambdaLn(2*Nplate + i) -= Fuz(i);
        
        // on assemble ce vecteur
        var.var_plate.Assemble(KlambdaLn);
        
        // on obtient ainsi Lambda^n+1
        for (int i = 0; i < 3*Nplate; i++)
          Lambda_np1(i) = 2.0*Lambda_n(i) - Lambda_nm1(i) - dt2*KlambdaLn(i)/MassLambda(i);
        
        var.var_plate.var_harmonic.ImposeNullDirichletCondition(Lambda_np1);        
        
        // on calcule ensuite P^n+1
        // (Mh + dt/2 Sh) P^n+1 = (Mh - dt/2 Sh) P^n 
        // - dt ( - F + Rh^S V^{n+1/2} + Ch^T (Lambda^n+1 - Lambda^n-1)/(2 dt) )

        // d'abord KhPn = (Mh - dt/2 Sh) P^n - dt (Rh^S V^n+1/2 - F)
        var.var_volume.ApplyOperatorRhScalar(1.0, t+0.5*dt, Vh_n, 0.0, KhPn);
        
        if (t < var.var_volume.GetFinalTimeSource())
          var.var_volume.AddScalarSourceAtTime(-1.0, t+0.5*dt, 0, KhPn);
        
        Mlt(-dt, KhPn);
        var.var_volume.ApplyOperatorDhMinusdtSh(1.0, t, Ph_n, 1.0, KhPn);
        
        // on assemble ce vecteur
        var.var_volume.Assemble(KhPn);
        
        // on calcule d'une autre part Ch^T (Lambda^n+1 - Lambda^n-1)/ 2
        for (int i = 0; i < Nplate; i++)
          U_np1(i) = 0.5*(Lambda_np1(2*Nplate + i) - Lambda_nm1(2*Nplate+i));
        
        Fvol.Fill(0);
        var.MltAddChTranspose(1.0, U_np1, Fvol);
        var.AssembleAcoustic(Fvol);
        
        var.ExchangeSurfaceToVolume(Fvol, ChLambdaMinus_n, ChLambdaPlus_n);
        
        // et on met a jour P^n+1 avec ces valeurs
        for (int i = 0; i < var.ddl_vol_minus.GetM(); i++)
          KhPn(var.ddl_vol_minus(i)) -= ChLambdaMinus_n(i);
        
        for (int i = 0; i < var.ddl_vol_plus.GetM(); i++)
          KhPn(var.ddl_vol_plus(i)) += ChLambdaPlus_n(i);
        
        // on inverse par (Mh + dt/2) Ph^n+1
        var.var_volume.SolveOperatorDhPlusdtSh(KhPn);

        // on obtient ainsi P^n+1
        Copy(Ph_n, Ph_nm1);
        Copy(KhPn, Ph_n);
        
        // on met a jour Lambda^n
        Copy(Lambda_n, Lambda_nm1);
        Copy(Lambda_np1, Lambda_n);
        
        // on met a jour V^n+1/2
        const VectReal_wp& BhMinusDtSh = var.var_volume.Glob_mat_Dh->GetDiagonalBhMinusdtSh();
        const VectReal_wp& invBhPlusDtSh = var.var_volume.Glob_mat_Dh->GetInverseDiagonalBhPlusdtSh();
        var.var_volume.ApplyOperatorRhVectorial(-dt, t+dt, Ph_n, 0.0, Vh_np1);
        Real_wp tmp;
        for (int i = 0; i < Vh_np1.GetM(); i++)
          {
            tmp = invBhPlusDtSh(i)*(Vh_np1(i) + BhMinusDtSh(i)*Vh_n(i));
            Vh_n(i) = tmp;
          }

        //int test_input; cout << "we wait" << endl; cin >> test_input;
        return;
      }
    
    // partie rigidite : Calcul des produits matrices vecteurs K_h P^n et K^L Lambda^n
    var.var_volume.ApplyOperatorKh(-1.0, t, Ph_n, 0.0, KhPn);
    var.var_plate.ApplyOperatorKh(-1.0, t, Lambda_n, 0.0, KlambdaLn);
    
    // rajout de la source
    var.var_volume.AddScalarSourceAtTime(-1.0, t, 1, KhPn);
    var.var_plate.AddScalarSourceAtTime(-1.0, t, 1, KlambdaLn);
    
    // on calcule KhPn = 2 Mh Pn - (Mh - dt/2 Sh) Pn-1 - dt^2 (Kh Pn - F)
    for (int i = 0; i < KhPn.GetM(); i++)
      KhPn(i) = -dt2*KhPn(i) + 2.0*Mh(i)*Ph_n(i) - Mh_minus_dt_Sh(i)*Ph_nm1(i);
    
    // on assemble ce vecteur
    var.var_volume.Assemble(KhPn);
    
    // et on lui rajoute dt/2 Ch* Lambda^n-1
    for (int i = 0; i < var.ddl_vol_minus.GetM(); i++)
      KhPn(var.ddl_vol_minus(i)) += 0.5*dt*ChLambdaMinus_nm1(i);
    
    for (int i = 0; i < var.ddl_vol_plus.GetM(); i++)
      KhPn(var.ddl_vol_plus(i)) -= 0.5*dt*ChLambdaPlus_nm1(i);

    // et on inverse par (Mh + dt/2) Sh
    // au final on a KhPn qui est un P^n+1 sans la partie en Lambda^n+1 
    // KhPn = (Mh + dt/2 Sh)^-1 [ 2 Mh P^n  - (Mh - dt/2 Sh) P^(n-1) 
    // - dt^2 Kh P^n + dt/2 Ch* Lambda^(n-1) ]
    for (int i = 0; i < KhPn.GetM(); i++)
      KhPn(i) *= invMh_plus_dt_Sh(i);
    
    // on calcule Fuz = Ch (KhPn - P^(n-1))
    for (int i = 0; i < var.ddl_vol_minus.GetM(); i++)
      Gvol_minus(i) = KhPn(var.ddl_vol_minus(i)) - Ph_nm1(var.ddl_vol_minus(i));
    
    for (int i = 0; i < var.ddl_vol_plus.GetM(); i++)
      Gvol_plus(i) = KhPn(var.ddl_vol_plus(i)) - Ph_nm1(var.ddl_vol_plus(i));
    
    var.ExchangeVolumeToSurface(Gvol_minus, Gvol_plus, Gsurf_minus, Gsurf_plus);
    for (int i = 0; i < Nacous; i++)
      Fvol(i) = Gsurf_minus(i) - Gsurf_plus(i);
    
    Fuz.Fill(0);
    var.MltAddCh(1.0, Fvol, Fuz);
    
    // on calcule KlambdaLn = -dt^2 K^L Lambda^n  + dt/2 Ch (KhPn - P^(n-1))
    Mlt(-dt2, KlambdaLn);
    for (int i = 0; i < Nplate; i++)
      KlambdaLn(2*Nplate + i) += 0.5*dt*Fuz(i);
    
    // on assemble ce vecteur
    var.var_plate.Assemble(KlambdaLn);
    
    // on rajoute M^L (2 Lambda^n - Lambda^(n-1))
    for (int i = 0; i < 3*Nplate; i++)
      KlambdaLn(i) += MassLambda(i)*(2.0*Lambda_n(i) - Lambda_nm1(i));
    
    var.var_plate.var_harmonic.ImposeNullDirichletCondition(KlambdaLn);
    for (int i = 0; i < Nplate; i++)
      Fuz(i) = KlambdaLn(2*Nplate + i);
        
    // on resout le systeme :
    // (M^L + dt^2/4 Ch Mh^-1 Ch^*) Lambda^n+1 = Fuz   
    var.SolveSchur(U_np1, Fuz);
    
    // pour les ddls en teta, schema habituel
    for (int i = 0; i < Nplate; i++)
      {
        Lambda_np1(i) = KlambdaLn(i)/MassLambda(i);
        Lambda_np1(Nplate+i) = KlambdaLn(Nplate+i)/MassLambda(Nplate+i);
        // et pour u on prend la solution du complement de Schur
        Lambda_np1(2*Nplate+i) = U_np1(i);
      }
    
    var.var_plate.var_harmonic.ImposeNullDirichletCondition(Lambda_np1);
        
    // calcul de Ch* Lambda^(n+1)
    Fvol.Fill(0);
    var.MltAddChTranspose(1.0, U_np1, Fvol);
    var.AssembleAcoustic(Fvol);
    
    Copy(ChLambdaMinus_n, ChLambdaMinus_nm1);
    Copy(ChLambdaPlus_n, ChLambdaPlus_nm1);
    
    var.ExchangeSurfaceToVolume(Fvol, ChLambdaMinus_n, ChLambdaPlus_n);
    
    // et on met a jour P^n+1 avec ces valeurs
    for (int i = 0; i < var.ddl_vol_minus.GetM(); i++)
      KhPn(var.ddl_vol_minus(i)) -= 0.5*dt*ChLambdaMinus_n(i)
        *invMh_plus_dt_Sh(var.ddl_vol_minus(i));
    
    for (int i = 0; i < var.ddl_vol_plus.GetM(); i++)
      KhPn(var.ddl_vol_plus(i)) += 0.5*dt*ChLambdaPlus_n(i)*invMh_plus_dt_Sh(var.ddl_vol_plus(i));
    
    // on passe a l'itere suivant
    Copy(Ph_n, Ph_nm1);
    Copy(KhPn, Ph_n);
    
    Copy(Lambda_n, Lambda_nm1);
    Copy(Lambda_np1, Lambda_n);
  }
  
}

#define MONTJOIE_FILE_VIBRO_ACOUSTIC_SESSION_CXX
#endif
