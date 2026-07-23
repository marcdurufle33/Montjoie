#ifndef MONTJOIE_FILE_COUPLAGE_PIANO_CXX

namespace Montjoie
{
	
  //! constructeur par defaut
  CouplagePiano::CouplagePiano(int argc, char**argv)
  {
    CHRONO_NOTE = chrono_string.GetNumber();
    CHRONO_PLATE = chrono_string.GetNumber();
    CHRONO_COUPLAGE = chrono_string.GetNumber();
    CHRONO_ENERGY = chrono_string.GetNumber();
    CHRONO_FEM = chrono_string.GetNumber();
    display_max_time_single_proc = true;
    nb_iter_for_each_computation_time = 100;
    display_eigenvectors_plate = false;
    reprise_on_temporary_path = false;
    file_output_energy = "EnergyPiano.txt";
  }
  
  
  //! modification des parametres de la classe avec une ligne du fichier de donnees
  /*!
    \param[in] keyword mot-cle associe a la ligne du fichier de donnees
    \param[in] param liste des parametres associes au mot-cle
   */
  void CouplagePiano::SetInputData(const string & keyword, const VectString & param)
  {
#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX
    var_fluide.SetInputData(keyword,param);
#endif

    source_coupling_string.SetInputData(keyword, param, note.param_corde_data);
    
    if (keyword == "DisplayComputationTime")
      {
	if (param.GetM() <= 1)
	  {
	    cout << "In SetInputData of CouplagePiano" << endl;
	    cout << "DisplayComputationTime needs 2 parameters, for instance :" << endl;
	    cout << "DisplayComputationTime = GLOBAL 100" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }
	
        if (param(0) == "GLOBAL")
          display_max_time_single_proc = false;
        else
          display_max_time_single_proc = true;
        
        nb_iter_for_each_computation_time = to_num<int>(param(1));
      }
    else if (keyword == "DisplayEigenvectors")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of CouplagePiano" << endl;
	    cout << "DisplayEigenvectors needs one parameter, for instance :" << endl;
	    cout << "DisplayEigenvectors = YES" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

        if (param(0) == "YES")
          display_eigenvectors_plate = true;
        else
          display_eigenvectors_plate = false;
      }
    else if (keyword == "RegroupReprise")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of CouplagePiano" << endl;
	    cout << "RegroupReprise needs one parameter, for instance :" << endl;
	    cout << "RegroupReprise = N" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

        int nb_proc_per_node = to_num<int>(param(0));
        output_reprise.RegroupWritingOnMaster(nb_proc_per_node);
      }
    else if (keyword == "FileOutputEnergyPiano")
      file_output_energy = param(0);
  }
  
  
  //! affichage des parametres des cordes
  void CouplagePiano::Display()
  {
#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX

#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_fluide.var_vibro.var_volume.var_harmonic.comm_group_mode;
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int rank_proc(0);
#endif

#else
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
    int rank_proc(0);
#endif

#endif

    if (rank_proc == 0)
      {
        string file_name = note.DOSSIER + note.file_storage_param;
        ofstream file_out(file_name.data());
        note.Display(file_out, true);
        //~ soundboard.Display(file_out, true);
        file_out.close();
      }
  }
  
  
  //! construction des tableaux necessaires aux iterations en temps
  void CouplagePiano::ConstructAll(const string& data_file)
  {
#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX
    
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_fluide.var_vibro.var_volume.var_harmonic.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    chrono.SetCommunicator(var_fluide.var_vibro.var_volume.var_harmonic.comm_group_mode);
#else
    int nb_proc(1), rank_proc(0);
#endif
    
    bool maitre = (rank_proc == 0);
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    ReadLinesFile(data_file, lines_data_file, comm);
#else
    ReadLinesFile(data_file, lines_data_file);
#endif
    
    var_fluide.InitIndices(PhysicalConstant::nb_max_indices, PhysicalConstant::nb_max_indices);    
    var_fluide.var_vibro.var_volume.var_harmonic.SetTypeEquation("TIME_ACOUSTIC");
    var_fluide.var_vibro.var_plate.var_harmonic.SetTypeEquation("REISSNER_MINDLIN");
    ReadInputFile(lines_data_file, *this);

    string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
    string prefix_load = var_fluide.var_vibro.var_volume.path_reprise;
    this->prefix_save = prefix_load;
    if (this->reprise_on_temporary_path)
      {
        int pos = prefix_load.substr(0, prefix_load.size()-1).find_last_of("/");
        this->prefix_save = "/tmp/"+prefix_load.substr(pos+1, prefix_load.size()-pos-1);
        //string command_line = "mkdir " + prefix_save;
        //system(command_line.data());
      }
    
    IVect num_iter;
    //DISP(var_fluide.var_vibro.var_volume.load_reprise);
    if (var_fluide.var_vibro.var_volume.load_reprise)
      {
        // on lit les donnees de la reprise
        if (rank_proc == 0)
          {
            num_iter.Read(prefix_load+"IterationPianoSave"+suffix);
            if (num_iter.GetM() != 3)
              {
                // on pense alors que l'ecriture critique a ete interrompu,
                // on regarde l'autre fichier dans ce cas
                num_iter.Read(prefix_load+"IterationPianoSaveNew"+suffix);
                
                if (num_iter.GetM() != 3)
                  {
                    cout << "Echec de la reprise" << endl;
                    cout << "Veuillez corriger le probleme" << endl;
                    abort();
                  }
              }
          }
        else
          num_iter.Reallocate(3);
        
        // on broadcaste num_iter sur tous les procs
#ifdef SELDON_WITH_MPI
        MPI_Bcast(num_iter.GetData(), num_iter.GetM(), MPI_INTEGER, 0,
                  var_fluide.var_vibro.var_volume.var_harmonic.comm_group_mode);
#endif
                    
        var_fluide.var_vibro.var_volume.load_iter_reprise = num_iter(0);
        
        // on lit les infos sur la force du marteau
        if (maitre)
          {
            note.load_reprise = true;
            note.load_iter_reprise = num_iter(0);
            VectReal_wp tmp;
            tmp.Read(prefix_load + "ForceMarteauPianoSave"+suffix);
            note.beginning_time_force = tmp(0);
            note.last_time_force = tmp(1);
            note.final_time_shank = tmp(2);
            note.hammer.energy_free_hammer = tmp(3);
            note.shank.free_kinetic_energy = tmp(4);
            note.shank.free_potential_energy = tmp(5);
            note.shank.free_internal_energy = tmp(6);
            note.shank.free_shank_energy = tmp(7);
          }
      }
    
    //cout << "Marqueur 0" << endl;
    if (maitre)
      {
	// constructing arrays for the resolution of 1-D strings
        int num_dossier = -1;
        bool remove_file = true;
        if (var_fluide.var_vibro.var_volume.load_reprise)
          {
            num_dossier = num_iter(1);
            remove_file = false;
          }
        
	ReadInputFile(lines_data_file, note);
	note.name_data_file = data_file;
	if (note.piano_strings.GetM() > 0)
	  {
	    note.ConstructAll(num_dossier);
	    note.Init();
	    Display();
            cout << "Initialisation des cordes OK"<<endl;
            
            output_ener_all_but_note.Init(note.DOSSIER+file_output_energy,
                                          note.size_buffer_ener, remove_file);
          }
      }

#ifdef SELDON_WITH_MPI
    MPI_Bcast_string(note.DOSSIER, 0, comm);
#endif
    
    //cout << "Marqueur 1" << endl;
    // changing output directories for vibroacoustic 
    for(int i = 0; i < var_fluide.var_vibro.var_volume.var_harmonic.output_grid_param.GetM(); i++)
      {
	string tototemp = var_fluide.var_vibro.var_volume.var_harmonic
          .output_grid_param(i).GetTotalFieldFile();
	
        var_fluide.var_vibro.var_volume.var_harmonic
          .output_grid_param(i).SetTotalFieldFile(note.DOSSIER + tototemp);
      }
    
    for(int i = 0; i < var_fluide.var_vibro.var_plate.var_harmonic.output_grid_param.GetM(); i++)
      {
	string tototemp = var_fluide.var_vibro.var_plate.var_harmonic
          .output_grid_param(i).GetTotalFieldFile();
	
        var_fluide.var_vibro.var_plate.var_harmonic
          .output_grid_param(i).SetTotalFieldFile(note.DOSSIER + tototemp);
        
        tototemp = var_fluide.var_vibro.var_plate.var_harmonic
          .output_grid_param(i).GetDiffractedFieldFile();
	
        var_fluide.var_vibro.var_plate.var_harmonic
          .output_grid_param(i).SetDiffractedFieldFile(note.DOSSIER + tototemp);
      }
    
    for(int i = 0; i < var_fluide.var_vibro.var_volume.var_harmonic.output_mesh_param.GetM(); i++)
      {
	string tototemp = var_fluide.var_vibro.var_volume.var_harmonic
          .output_mesh_param(i).GetTotalFieldFile();
	
        var_fluide.var_vibro.var_volume.var_harmonic
          .output_mesh_param(i).SetTotalFieldFile(note.DOSSIER + tototemp);
      }
    
    for(int i = 0; i < var_fluide.var_vibro.var_plate.var_harmonic.output_mesh_param.GetM(); i++)
      {
	string tototemp = var_fluide.var_vibro.var_plate.var_harmonic
          .output_mesh_param(i).GetTotalFieldFile();
	
        var_fluide.var_vibro.var_plate.var_harmonic
          .output_mesh_param(i).SetTotalFieldFile(note.DOSSIER + tototemp);
      }    

    // on supprime les sismos mecaniques si deja existants
    string file_name; 
    for (int i = 0; i < var_fluide.var_vibro.var_plate.var_harmonic.output_grid_param.GetM(); i++)
      {
        if (var_fluide.var_vibro.var_plate.var_harmonic.output_grid_param(i).GetNature()
            == GridInterpolationFull<Dimension2>::POINT)
          file_name = var_fluide.var_vibro.var_plate.var_harmonic
            .output_grid_param(i).GetDiffractedFieldFile();
      }
    
    if ((file_name.size() > 0) && (!var_fluide.var_vibro.var_volume.load_reprise))
      remove(file_name.data());
    
    soundboard.InitOutput(var_fluide.var_vibro.var_plate.var_harmonic);

    {
      string tototemp = var_fluide.var_vibro.var_volume.var_harmonic
        .output_points_outside_param.GetTotalFieldFile();
      
      var_fluide.var_vibro.var_volume.var_harmonic
        .output_points_outside_param.SetTotalFieldFile(note.DOSSIER + tototemp);
      
      if (!var_fluide.var_vibro.var_volume.load_reprise)
        {
          file_name = note.DOSSIER + tototemp;
          remove(file_name.data());
        }
      
      tototemp = var_fluide.var_vibro.var_volume.var_harmonic.output_rcs_param.GetFileNameNormL2();
      var_fluide.var_vibro.var_volume.var_harmonic
        .output_rcs_param.SetFileNameNormL2(note.DOSSIER + tototemp);
      if (!var_fluide.var_vibro.var_volume.load_reprise)
        {
          file_name = note.DOSSIER + tototemp;
          remove(file_name.data());
        }
    }

    var_fluide.ConstructAcoustics();
    
    chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);    
    
    var_fluide.var_vibro.var_volume.var_harmonic.TreatDirichletCondition();
    var_fluide.var_vibro.var_plate.var_harmonic.TreatDirichletCondition();
    
    chrono.Stop(VirtualTimer::ALL);
    chrono.DisplayTime(VirtualTimer::ALL, "retrieve Dirichlet dofs");

    if (maitre)
      cout << "Acoustique OK"<<endl;
    
#ifdef SELDON_WITH_MPI
    soundboard.VarSoundboard.comm_group_mode = comm;	
#endif
    
    chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);    

    var_fluide.var_vibro.var_volume.ComputeRightHandSide();

    chrono.Stop(VirtualTimer::ALL);
    chrono.DisplayTime(VirtualTimer::ALL, "compute volume sources");
    
    chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);    
    
    soundboard.init_file = data_file;
    ReadInputFile(lines_data_file, soundboard);

    if( var_fluide.var_vibro.var_volume.save_reprise)
      {
        if (!soundboard.load_eigenvectors_asked)
          {
            soundboard.path_eigenvectors = prefix_load;
            soundboard.save_eigenvectors_asked = var_fluide.var_vibro.var_volume.save_reprise;
            soundboard.load_eigenvectors_asked = var_fluide.var_vibro.var_volume.load_reprise;
          }
      }    
    soundboard.PerformDiagonalization(var_fluide.var_vibro.var_plate.var_harmonic,
                                      soundboard.path_eigenvectors, soundboard.path_eigenvectors,
                                      lines_data_file, soundboard.save_eigenvectors_asked,
                                      soundboard.load_eigenvectors_asked,
                                      source_coupling_string, display_eigenvectors_plate);
    
    chrono.Stop(VirtualTimer::ALL);
    chrono.DisplayTime(VirtualTimer::ALL, "compute diagonalization of 2-D stiffness matrix");
    
    if (maitre)
      cout << "Diagonalization OK" << endl;
    
    if (nb_proc == 1)
      {
	var_fluide.var_vibro.var_plate.var_harmonic.mesh.Write("plate.mesh");
	var_fluide.var_vibro.var_volume.var_harmonic.mesh.Write("volume.mesh");
      }
    
    //var_fluide.var_vibro.var_plate.ComputeRightHandSide();
    
    chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);    
    
    soundboard.Deltat = var_fluide.var_vibro.var_volume.GetTimeStep();
    soundboard.ConstructAll();
    
    chrono.Stop(VirtualTimer::ALL);
    chrono.DisplayTime(VirtualTimer::ALL, "construct 2-D arrays for the soundboard");
    
    // soundboard.Rdelta
    // pour mettre une donnee initiale sur le soundboard : Initialise()
    //Real_wp rhodelta = var_fluide.var_vibro.var_plate.var_harmonic.ref_rho(1)
    // *var_fluide.var_vibro.var_plate.var_harmonic.ref_delta(1).GetConstant();
    //soundboard.Initialise(rhodelta);
            
    chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);    
    
    var_fluide.var_vibro.ComputeDiagonalSchur(MhCouplage, invMhSurfMinus, invMhSurfPlus,
                                              var_fluide.var_vibro.var_volume
                                              .var_harmonic.FirstOrderFormulation());
    
    soundboard.CalculeLU(MhCouplage, var_fluide.var_vibro);

    chrono.Stop(VirtualTimer::ALL);
    chrono.DisplayTime(VirtualTimer::ALL, "compute Cholesky factorisation of Schur complement");
    
    Real_wp scal_b = soundboard.ConstruitLambdaSchur();
    //~ note.bridge.SetScal_b(scal_b);
    if (maitre)
      note.bridge.SetVecScal_b(0,scal_b);

    ChLambdaNm1_Minus.Reallocate(var_fluide.var_vibro.ddl_vol_minus.GetM());
    ChLambdaNm1_Plus.Reallocate(var_fluide.var_vibro.ddl_vol_plus.GetM());
    ChLambdaNm1_Minus.Fill(0); ChLambdaNm1_Plus.Fill(0);
    
    if (maitre)
      {
	cout << "Initialisation of soundboard OK"<<endl;
	
	int Nvol = var_fluide.var_vibro.var_volume.var_harmonic.GetNbMeshDof();
        int nFEM = var_fluide.var_vibro.var_plate.var_harmonic.GetNbMeshDof();
        if (nb_proc > 1)
          {
            Nvol = var_fluide.var_vibro.var_volume.var_harmonic.GetNbGlobalMeshDof();
            nFEM = var_fluide.var_vibro.var_plate.var_harmonic.GetNbGlobalMeshDof();
          }
        
	cout << "Nombre de degres de liberte dans le volume : " << Nvol << endl;
	cout << "Nombre de degres de liberte sur la surface : " << nFEM << endl;
      }
    
    // on resette le chronometre
    chrono.Reset(VirtualTimer::ALL);
    
    if (var_fluide.var_vibro.var_volume.load_reprise)
      {        
        int nb_dof_u = var_fluide.Pn.GetM();
        int nb_dof_v = var_fluide.Vh_n.GetM();
        int nb_dof_lambda = soundboard.Lambda.GetM();
        int nb_dof_minus = ChLambdaNm1_Minus.GetM();
        int nb_dof_plus = ChLambdaNm1_Plus.GetM();
        int nb_dof_note = note.size;
        bool reprise_failed = false;
                  
        // l'utilisateur fait une reprise
        // on initialise les vecteurs avec ceux stockes en memoire
        suffix = to_str(num_iter(2)) + "_P" + suffix;
        output_reprise.StartReading(prefix_load + "SolPiano" + suffix);
        if (var_fluide.var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
          {
            output_reprise.Read(var_fluide.Pnm1);
            output_reprise.Read(var_fluide.Pn);
            output_reprise.Read(var_fluide.Vh_n);
            
            //var_fluide.Pnm1.Read(prefix_load+"Pnm1PianoSave"+suffix);
            //var_fluide.Pn.Read(prefix_load+"PnPianoSave"+suffix);
            //var_fluide.Vh_n.Read(prefix_load+"VnPianoSave"+suffix);
            
            if (var_fluide.Pnm1.GetM() != nb_dof_u)
              reprise_failed = true;

            if (var_fluide.Pn.GetM() != nb_dof_u)
              reprise_failed = true;
            
            if (var_fluide.Vh_n.GetM() != nb_dof_v)
              reprise_failed = true;
          }
        else
          {
            output_reprise.Read(var_fluide.Pnm1);
            output_reprise.Read(var_fluide.Pn);
            output_reprise.Read(var_fluide.KhPn);
            
            //var_fluide.Pnm1.Read(prefix_load+"Pnm1PianoSave"+suffix);
            //var_fluide.Pn.Read(prefix_load+"PnPianoSave"+suffix);
            //var_fluide.KhPn.Read(prefix_load+"KhPnPianoSave"+suffix);

            if (var_fluide.Pnm1.GetM() != nb_dof_u)
              reprise_failed = true;
            
            if (var_fluide.Pn.GetM() != nb_dof_u)
              reprise_failed = true;
            
            if (var_fluide.KhPn.GetM() != nb_dof_u)
              reprise_failed = true;
          }
        
        output_reprise.Read(soundboard.Lambda);
        output_reprise.Read(soundboard.LambdaPoint);
        
        //soundboard.Lambda.Read(prefix_load+"LambdaPianoSave"+suffix);
        //soundboard.LambdaPoint.Read(prefix_load+"LambdaPointPianoSave"+suffix);
        
        if (soundboard.Lambda.GetM() != nb_dof_lambda)
          reprise_failed = true;
        
        if (soundboard.LambdaPoint.GetM() != nb_dof_lambda)
          reprise_failed = true;
        
        output_reprise.Read(ChLambdaNm1_Minus);
        output_reprise.Read(ChLambdaNm1_Plus);
        
        //ChLambdaNm1_Minus.Read(prefix_load+"ChLambdaMinusPianoSave"+suffix);
        //ChLambdaNm1_Plus.Read(prefix_load+"ChLambdaPlusPianoSave"+suffix);        
        
        if (ChLambdaNm1_Minus.GetM() != nb_dof_minus)
          reprise_failed = true;
        
        if (ChLambdaNm1_Plus.GetM() != nb_dof_plus)
          reprise_failed = true;
        
        if (rank_proc == 0)
          {
            output_reprise.Read(note.val, true);
            
            //note.val.Read(prefix_load+"ValPianoSave"+suffix);
            if (note.val.GetM() != nb_dof_note)
              reprise_failed = true;
            
            //DISP(MaxAbs(note.val));
          }
        
        Real_wp temps = num_iter(0)*var_fluide.var_vibro.var_volume.GetTimeStep();
        for (int i = 0; i < var_fluide.var_vibro.var_volume
               .var_harmonic.output_grid_param.GetM(); i++)
          var_fluide.var_vibro.var_volume.var_harmonic
            .output_grid_param(i).ChangeTime(num_iter(0), temps);
        
        for(int i = 0; i < var_fluide.var_vibro.var_plate.var_harmonic
              .output_grid_param.GetM(); i++)
          var_fluide.var_vibro.var_plate.var_harmonic
            .output_grid_param(i).ChangeTime(num_iter(0), temps);
        
        for(int i = 0; i < var_fluide.var_vibro.var_volume
              .var_harmonic.output_mesh_param.GetM(); i++)
          var_fluide.var_vibro.var_volume.var_harmonic
            .output_mesh_param(i).ChangeTime(num_iter(0), temps);
        
        for(int i = 0; i < var_fluide.var_vibro.var_plate
              .var_harmonic.output_mesh_param.GetM(); i++)
          var_fluide.var_vibro.var_plate.var_harmonic
            .output_mesh_param(i).ChangeTime(num_iter(0), temps);
        
        var_fluide.var_vibro.var_volume.var_harmonic.output_rcs_param.LoadDatas(output_reprise);
        
        output_reprise.FinalizeReading();
        
        if (reprise_failed)
          {
            cout << "Echec de la reprise" << endl;
            cout << "Veuillez corriger le probleme" << endl;
            abort();
          }
      }
    
    return;
    
    /*
    //string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
    soundboard.Pu.Write("Pu"+suffix);
    var_fluide.var_vibro.ddl_vol_minus.Write("DdlVolMinus"+suffix);
    var_fluide.var_vibro.ddl_vol_plus.Write("DdlVolPlus"+suffix);
    
    var_fluide.Cbarre.WriteText("ChVol"+suffix);
    
    soundboard.MatriceResolExacte0.WriteText("S0dt"+suffix);
    soundboard.MatriceResolExacte1.WriteText("S1dt"+suffix);
    soundboard.MatriceResolExacteSecondMembre.WriteText("Rdt"+suffix);
    
    soundboard.MatriceResolExacte0Point.WriteText("S0dt_point"+suffix);
    soundboard.MatriceResolExacte1Point.WriteText("S1dt_point"+suffix);
    soundboard.MatriceResolExacteSecondMembrePoint.WriteText("Rdt_point"+suffix);
    soundboard.vect_masse.WriteText("LambdaMh"+suffix);
    soundboard.vect_lambda.WriteText("LambdaDh"+suffix);
    
    if (nb_proc == 1)
      {
        DISP(var_fluide.Pn.GetM());
        DISP(var_fluide.Vh_n.GetM());
        var_fluide.Pn.FillRand(); Mlt(1e-9, var_fluide.Pn);
        var_fluide.Pnm1.FillRand(); Mlt(1e-9, var_fluide.Pnm1);
        var_fluide.Vh_n.FillRand(); Mlt(1e-9, var_fluide.Vh_n);
        var_fluide.Vh_nm1.FillRand(); Mlt(1e-9, var_fluide.Vh_nm1);
      }
    else
      {
        int nb_dof_u = var_fluide.var_vibro.var_volume.var_harmonic.GetNbGlobalMeshDof();
        int nb_dof_v = var_fluide.var_vibro.var_volume.var_harmonic.GetNbGlobalDof() - nb_dof_u;
        DISP(nb_dof_u); DISP(nb_dof_v);
        VectReal_wp Pn(nb_dof_u), Pnm1(nb_dof_u), Vh_n(nb_dof_v), Vh_nm1(nb_dof_v);
        Pn.FillRand(); Mlt(1e-9, Pn);
        Pnm1.FillRand(); Mlt(1e-9, Pnm1);
        Vh_n.FillRand(); Mlt(1e-9, Vh_n);
        Vh_nm1.FillRand(); Mlt(1e-9, Vh_nm1);
        
        for (int i = 0; i < var_fluide.Pn.GetM(); i++)
          var_fluide.Pn(i) = Pn(num_vol(i));
        
        for (int i = 0; i < var_fluide.Pnm1.GetM(); i++)
          var_fluide.Pnm1(i) = Pnm1(num_vol(i));
        
        int Nvol = var_fluide.var_vibro.var_volume.var_harmonic.mesh_num.GetNbDof();
        for (int i = 0; i < var_fluide.Vh_n.GetM(); i++)
          var_fluide.Vh_n(i) = Vh_n(num_vol(Nvol + i)-nb_dof_u);
        
        for (int i = 0; i < var_fluide.Vh_nm1.GetM(); i++)
          var_fluide.Vh_nm1(i) = Vh_nm1(num_vol(Nvol + i)-nb_dof_u);
        
      }
    
    DISP(soundboard.Lambda.GetM());
    soundboard.Lambda.FillRand(); Mlt(1e-9, soundboard.Lambda);
    soundboard.LambdaPoint.FillRand(); Mlt(1e-9, soundboard.LambdaPoint);
    
    int nFEM = var_fluide.var_vibro.var_plate.var_harmonic.mesh_num.GetNbDof();
    VectReal_wp Force(nFEM); Force.Fill(0);
    Mlt(soundboard.Pu, soundboard.Lambda, Force);
    VectReal_wp FP(var_fluide.var_vibro.nodl_acous);
    FP.Fill(0);
    var_fluide.var_vibro.MltAddChTranspose(1.0, Force, FP);
    var_fluide.var_vibro.AssembleAcoustic(FP);
    
    // sending surface values to volume dofs
    var_fluide.var_vibro.
      ExchangeSurfaceToVolume(FP, ChLambdaNm1_Minus, ChLambdaNm1_Plus);

    ChLambdaNm1_Minus.Write("ChLambdaMinus"+suffix);
    ChLambdaNm1_Plus.Write("ChLambaPlus"+suffix);
    var_fluide.Kh.mat_iterative_sym.WriteText("Kh.dat");
    */

#else
    cout << "Recompile by including FluideExterieur.hxx" << endl;
    abort();
#endif
  }
  
  
  //! on sauve les instantanes en prevision d'une reprise
  void CouplagePiano::SaveDatas(int nt)
  {
#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX

#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_fluide.var_vibro.var_volume.var_harmonic.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif
    
    if (var_fluide.var_vibro.var_volume.save_reprise)
      if (nt % var_fluide.var_vibro.var_volume.save_iter_reprise == 0)
        {
          if (rank_proc == 0)
            cout << "Saving iterates... " << endl;
          
          chrono.Start(VirtualTimer::OUTPUT);
          
          // prefixe des fichiers de sauvegardes, et suffixe
          string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
          string prefix = this->prefix_save;
          
          // on cherche le numero de dossier dans STIFFOUT
          IVect num_iter(3);
          num_iter(0) = nt;
          int num_dossier = 0;
          if (rank_proc == 0)
            {
              size_t pos2 = note.DOSSIER.find_last_of("/");
              string chaine = note.DOSSIER.substr(0, pos2);
              size_t pos1 = chaine.find_last_of("/");
              num_dossier = to_num<int>(note.DOSSIER.substr(pos1+1, pos2-pos1-1));
              //DISP(num_dossier);
            }
          
          bool first_save = false;
          if (var_fluide.var_vibro.var_volume.bascule_reprise == -1)
            {
              // premiere fois qu'on imprime
              var_fluide.var_vibro.var_volume.bascule_reprise = 0;
              first_save = true;
            }
            
          // on imprime dans le fichier IterationPianoSave
          // le numero de dossier, la bascule de la reprise, le numero d'iteration
          num_iter(1) = num_dossier;
          num_iter(2) = var_fluide.var_vibro.var_volume.bascule_reprise;
          
          if (rank_proc == 0)
            {
              if (first_save)
                num_iter.Write(prefix+"IterationPianoSaveNew"+suffix);
              else
                {
                  // ecriture critique ou on affirme que les 
                  // anciens iteres sauves sont ceux qu'il faut considerer
                  // pour une reprise
                  IVect num_iter_old;
                  num_iter_old.Read(prefix+"IterationPianoSaveNew"+suffix);
                  num_iter_old.Write(prefix+"IterationPianoSave"+suffix);
                  
                  // fin de l'ecriture critique                  
                  num_iter.Write(prefix+"IterationPianoSaveNew"+suffix);
                }
              
              // on ecrit aussi les infos sur la force du marteau
              VectReal_wp tmp(4);
              tmp(0) = note.beginning_time_force;
              tmp(1) = note.last_time_force;
              tmp(2) = note.final_time_shank;
              tmp(3) = note.hammer.energy_free_hammer;
              
              tmp.Write(prefix + "ForceMarteauPianoSave"+suffix);
            }
          
          // on ecrit ici les nouveaux iteres (qui ne seront pas consideres pour la reprise)
          // le programme peut etre interrompu durant l'ecriture sans dommage collateral
          suffix = to_str(num_iter(2)) + "_P" + suffix;          
          output_reprise.StartWriting(prefix + "SolPiano" + suffix);
          
          if (var_fluide.var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
            {
              output_reprise.Write(var_fluide.Pnm1);
              output_reprise.Write(var_fluide.Pn);
              output_reprise.Write(var_fluide.Vh_n);
              
              //var_fluide.Pnm1.Write(prefix+"Pnm1PianoSave"+suffix);
              //var_fluide.Pn.Write(prefix+"PnPianoSave"+suffix);
              //var_fluide.Vh_n.Write(prefix+"VnPianoSave"+suffix);
            }
          else
            {
              output_reprise.Write(var_fluide.Pnm1);
              output_reprise.Write(var_fluide.Pn);
              output_reprise.Write(var_fluide.KhPn);

              //var_fluide.Pnm1.Write(prefix+"Pnm1PianoSave"+suffix);
              //var_fluide.Pn.Write(prefix+"PnPianoSave"+suffix);
              //var_fluide.KhPn.Write(prefix+"KhPnPianoSave"+suffix);
            }
          
          output_reprise.Write(soundboard.Lambda);
          output_reprise.Write(soundboard.LambdaPoint);
                    
          //soundboard.Lambda.Write(prefix+"LambdaPianoSave"+suffix);
          //soundboard.LambdaPoint.Write(prefix+"LambdaPointPianoSave"+suffix);
          
          output_reprise.Write(ChLambdaNm1_Minus);
          output_reprise.Write(ChLambdaNm1_Plus);
          
          //ChLambdaNm1_Minus.Write(prefix+"ChLambdaMinusPianoSave"+suffix);
          //ChLambdaNm1_Plus.Write(prefix+"ChLambdaPlusPianoSave"+suffix);        
          
          if (rank_proc == 0)
            {
              note.FlushBuffer();
              output_reprise.Write(note.val, true);
              //note.val.Write(prefix+"ValPianoSave"+suffix);
            }          
          
          var_fluide.var_vibro.var_volume.var_harmonic.output_rcs_param.SaveDatas(output_reprise);

          var_fluide.var_vibro.var_volume.bascule_reprise
            = 1 - var_fluide.var_vibro.var_volume.bascule_reprise;
                   
          output_reprise.FinalizeWriting();
          
#ifdef SELDON_WITH_MPI          
          MPI_Barrier(comm);
#endif

          chrono.Stop(VirtualTimer::OUTPUT);
          
          if (rank_proc == 0)
            cout << "Iterates saved" << endl;
        }
#else
    cout << "Recompile by including FluideExterieur.hxx" << endl;
    abort();
#endif
  }

  
  //! construction des tableaux necessaires aux iterations en temps (pas d'air ici)
  void CouplagePiano::ConstructAllImpedance(const string& data_file)
  {
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    ReadLinesFile(data_file, lines_data_file, MPI_COMM_WORLD);
#else
    ReadLinesFile(data_file, lines_data_file);
#endif

    ReadInputFile(lines_data_file, note);
    note.name_data_file = data_file;
    note.ConstructAll();
    note.Init();
    cout << "Corde : OK" << endl;
    
    output_ener_all_but_note.Init(note.DOSSIER+file_output_energy, note.size_buffer_ener);
    
    soundboard.init_file = data_file;
    soundboard.Deltat = note.Deltat;
    
    ReadInputFile(lines_data_file, soundboard);
    soundboard.ReadImpedance(note.DOSSIER);

    // soundboard.RunSoundboardDiagonalization();

    Real_wp scal_b = soundboard.ConstruitMatricesResolExacte();
    //~ note.bridge.SetScal_b(scal_b);
    note.bridge.SetVecScal_b(0,scal_b);
    
    soundboard.ConstruitLambda0etPoint(); // pour pouvoir utiliser le ActualiseLambda(coef);
    
    cout << "Soundboard : OK"<< endl;
  }
  

  //! construction des tableaux necessaires aux iterations en temps (pas d'air ici)  
  void CouplagePiano::ConstructAllStringSoundboard(const string& data_file)
  {
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
    ReadLinesFile(data_file, lines_data_file, MPI_COMM_WORLD);
#else
    int rank_proc(0);
    ReadLinesFile(data_file, lines_data_file);
#endif
    
    if (rank_proc == 0)
      {
        ReadInputFile(lines_data_file, note);
        note.name_data_file = data_file;
        note.ConstructAll();
        note.Init();
        
        cout << "Corde : OK" << endl;
        output_ener_all_but_note.Init(note.DOSSIER+file_output_energy, note.size_buffer_ener);
      }
    
    soundboard.init_file = data_file;
    ReadInputFile(lines_data_file, soundboard);
    soundboard.Deltat = note.Deltat;    
    
    soundboard.RunSoundboardDiagonalization(note.DOSSIER, source_coupling_string,
                                            lines_data_file, true);
    
    Real_wp scal_b = soundboard.ConstructAll();
    if(rank_proc == 0)
      note.bridge.SetVecScal_b(0,scal_b);
    
    for(int i = 0; i < soundboard.VarSoundboard.output_grid_param.GetM(); i++)
      {
	string tototemp = soundboard.VarSoundboard.output_grid_param(i).GetTotalFieldFile();
	soundboard.VarSoundboard.output_grid_param(i).SetTotalFieldFile(note.DOSSIER + tototemp);
        tototemp = soundboard.VarSoundboard.output_grid_param(i).GetDiffractedFieldFile();
	soundboard.VarSoundboard.output_grid_param(i)
          .SetDiffractedFieldFile(note.DOSSIER + tototemp);
      }
    

    
    for(int i = 0; i < soundboard.VarSoundboard.output_mesh_param.GetM(); i++)
      {
	string tototemp = soundboard.VarSoundboard.output_mesh_param(i).GetTotalFieldFile();
	soundboard.VarSoundboard.output_mesh_param(i)
          .SetTotalFieldFile(note.DOSSIER + tototemp);
      }    

    soundboard.InitOutput(soundboard.VarSoundboard);

    if(rank_proc == 0)
      cout << "Soundboard : OK"<< endl;        
  }
  
  
  //! calcul des energies des sous-systemes et ecriture
  Real_wp CouplagePiano::WriteEnergy()
  {
    chrono_string.Start(CHRONO_ENERGY);
    // reflechir au pb du decalage entre fluide et les autres
    Real_wp deltat_ener = note.deltat_ener;
    Real_wp t_courant = note.t_courant;
    Real_wp t_begin_ener = note.t_begin_ener;
    Real_wp t_end_ener = note.t_end_ener;
    Real_wp Deltat = note.Deltat;
    
    // initialisation
    Real_wp    ener_fluide =0.0;
    Real_wp ener_plaque_potential = 0.0;
    Real_wp ener_plaque_kinetic = 0.0;
        
    // energie des cordes
    Real_wp ener_cordes = 0.0;
    bool comes_from_elsewhere = true;
    //    VectReal_wp energy_sub;
    if (deltat_ener != 0)
      ener_cordes = note.WriteEnergy(comes_from_elsewhere);//le calcule et l ecrit sur un fichier
      //      ener_cordes = note.ComputeEnergy(energy_sub);
    
    Real_wp ener_tot = 0;
    Real_wp t_display = note.t_display;
    
    if (deltat_ener != 0)
      {
        int n = toInteger(round(t_courant/deltat_ener));
        TinyVector<Real_wp, 4> my_ener, previous_energy_tmp;
        bool write_energy = false;

#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX
        bool store_energy = false;
#endif
        
        if ((t_courant>= t_begin_ener) && (t_courant<= t_end_ener)
            && (abs(t_courant/deltat_ener - n) <(Deltat*0.5)))
          { 
#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX
            ener_fluide = var_fluide.GetEnergy();
#endif

            my_ener(2) = ener_fluide;
            
            ener_plaque_kinetic = soundboard.GetKineticEnergy();
            ener_plaque_potential = soundboard.GetPotentialEnergy();		
            
            my_ener(0) = ener_plaque_kinetic;
            my_ener(1) = ener_plaque_potential;
            
            write_energy = true;
            my_ener(3) = ener_fluide+ener_plaque_kinetic+ener_plaque_potential+ener_cordes;
            ener_tot = my_ener(3);
          }       
        
#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX
        if (!var_fluide.var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
          {
            if (write_energy)
              {
                previous_energy_tmp = my_ener;
                previous_energy_tmp(2) = ener_cordes;
                store_energy = true;
                write_energy = false;
              }
            
            t_courant -= var_fluide.var_vibro.var_volume.GetTimeStep();
            n = toInteger(round(t_courant/deltat_ener));
            if ((t_courant>= t_begin_ener) && (t_courant<= t_end_ener)
                && (abs(t_courant/deltat_ener - n) <(Deltat*0.5)))  
              {
                ener_fluide = var_fluide.GetEnergy();
                my_ener(2) = ener_fluide;
                
                my_ener(0) = previous_energy(0);
                my_ener(1) = previous_energy(1);
                
                my_ener(3) = my_ener(0) + my_ener(1) + previous_energy(2) + ener_fluide;
                ener_tot = my_ener(3);
                
                write_energy = true;
              }
            
            if (store_energy)
              previous_energy = previous_energy_tmp;
          }
#endif
        
        if (write_energy)
          {
            n = toInteger(round(t_courant/(t_display)));
            if (abs(t_courant/t_display - n) < Deltat/2)   
              {
                cout << "At time t = " << t_courant << " total energy = " << ener_tot << endl;
              }

            glob_chrono.Start(VirtualTimer::OUTPUT);
            output_ener_all_but_note.AddTinyVect(my_ener);
            glob_chrono.Stop(VirtualTimer::OUTPUT);
            
          }
      }
    else
      {
        int n = toInteger(round(t_courant/(t_display)));
        if (abs(t_courant/t_display - n) < Deltat/2)   
          {
            if (note.nb_dof_before.GetM() > 0)
              {
                Real_wp Lmax = 0.0, vmax = 0.0;
                int nb_corde = note.nb_corde;
                for (int i = 0; i < note.dimension*note.nb_dof_before(nb_corde); i++)
                  if (abs(note.val(i, 1)) > vmax)
                    vmax = abs(note.val(i, 1));
                
                for (int i = 0; i < soundboard.Lambda.GetM(); i++)
                  if (abs(soundboard.Lambda(i)) > Lmax)
                    Lmax = abs(soundboard.Lambda(i));
                
                cout << "At time t = " << t_courant << " || u,v,phi ||_inf = "
                     << vmax << ",  || Lambda||_inf = " << Lmax << endl;
              }
          }
      }
    
    chrono_string.Stop(CHRONO_ENERGY);
    return ener_tot;
  }
  

  //! Affichage des temps de calculs des differentes etapes de calcul
  void CouplagePiano::DisplayTimeChronos(int nt)
  {
#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX
    // on affiche les temps cumules (sommes des temps sur tous les processeurs)
    // et le temps de restitution (maximum des temps sur tous les processeurs)
    if (nt%nb_iter_for_each_computation_time != 0)
      return;
    
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_fluide.var_vibro.var_volume.var_harmonic.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif

    // Temps de calcul sur la rigidite 3-D
    double dt_loc, dt_sum, dt_min, dt_max, dt_loc2, dt_sum2, dt_max2, dt_min2;
    var_fluide.chrono.GetGlobalSeconds(VirtualTimer::STIFFNESS,
                                       dt_loc, dt_sum, dt_min, dt_max);
    
    if (rank_proc == 0)
      {
        cout << "Temps de calcul cumule pour la rigidite 3-D (RhS, Kh, RhV) : "
             << dt_sum << endl;
        
        if ((nb_proc > 1 ) && display_max_time_single_proc)
          cout << "Temps de calcul max sur un processeur pour la rigidite 3-D : "
               << dt_max << endl;
      }
    
    // temps de calcul sur la masse 3-D
    var_fluide.chrono.GetGlobalSeconds(VirtualTimer::MASS,
                                       dt_loc, dt_sum, dt_min, dt_max);
    
    if (rank_proc == 0)
      {
        cout << "Temps de calcul cumule pour la masse 3-D (Mh, Bh, Sh, ShV) : "
             << dt_sum << endl;
        
        if ((nb_proc > 1 ) && display_max_time_single_proc)
          cout << "Temps de calcul max sur un processeur pour la masse 3-D : "
               << dt_max << endl;
      }
    
    // temps de calcul pour la matrice Ch
    var_fluide.chrono.GetGlobalSeconds(VirtualTimer::FLUX,
                                       dt_loc, dt_sum, dt_min, dt_max);

    soundboard.chrono.GetGlobalSeconds(VirtualTimer::FLUX,
                                       dt_loc2, dt_sum2, dt_min2, dt_max2);
    
    if (rank_proc == 0)
      {
        cout << "Temps de calcul cumule pour le produit avec Ch : "
             << dt_sum+dt_sum2 << endl;
        
        if ((nb_proc > 1 ) && display_max_time_single_proc)
          cout << "Temps de calcul max sur un processeur pour Ch : "
               << dt_max+dt_max2 << endl;
      }
    
    // temps de calcul pour la matrice P
    chrono.GetGlobalSeconds(VirtualTimer::PROD,
                            dt_loc, dt_sum, dt_min, dt_max);

    soundboard.chrono.GetGlobalSeconds(VirtualTimer::PROD,
                                       dt_loc2, dt_sum2, dt_min2, dt_max2);
    
    if (rank_proc == 0)
      {
        cout << "Temps de calcul cumule pour le produit avec Pu : "
             << dt_sum+dt_sum2 << endl;
        
        if ((nb_proc > 1 ) && display_max_time_single_proc)
          cout << "Temps de calcul max sur un processeur pour Pu : "
               << dt_max+dt_max2 << endl;
      }
    
    // temps de calcul pour le schema 3-D
    var_fluide.chrono.GetGlobalSeconds(VirtualTimer::SCHEME,
                                       dt_loc, dt_sum, dt_min, dt_max);

    if (rank_proc == 0)
      {
        cout << "Temps de calcul cumule pour le schema 3-D : "
             << dt_sum << endl;
        
        if ((nb_proc > 1 ) && display_max_time_single_proc)
          cout << "Temps de calcul max sur un processeur pour le schema 3-D : "
               << dt_max << endl;
      }

    // temps de calcul pour le schema 2-D
    soundboard.chrono.GetGlobalSeconds(VirtualTimer::SCHEME,
                                       dt_loc, dt_sum, dt_min, dt_max);

    if (rank_proc == 0)
      {
        cout << "Temps de calcul cumule pour le schema 2-D (calcul duplique) : "
             << dt_sum << endl;
        
        if ((nb_proc > 1 ) && display_max_time_single_proc)
          cout << "Temps de calcul max sur un processeur pour le schema 2-D : "
               << dt_max << endl;
      }
    
    // temps de calcul pour l'energie
    var_fluide.chrono.GetGlobalSeconds(VirtualTimer::EXTRAPOL,
                                       dt_loc, dt_sum, dt_min, dt_max);

    dt_sum2 = chrono_string.GetSeconds(CHRONO_ENERGY);

    if (rank_proc == 0)
      cout << "Temps de calcul cumule pour l'energie : " << dt_sum + dt_sum2 << endl;
    
    // temps de calcul pour les cordes
    if (rank_proc == 0)
      {
        dt_sum = chrono.GetSeconds(VirtualTimer::FLUX);
        cout << "Temps de calcul sur le processeur maitre pour la corde : " << dt_sum << endl;
      }
    
    // temps de calcul pour les entrees sorties
    chrono.GetGlobalSeconds(VirtualTimer::OUTPUT, dt_loc, dt_sum, dt_min, dt_max);
    if (rank_proc == 0)
      cout << "Temps de calcul pour les entrees/sorties : " << dt_sum << endl;
    
    // temps de calcul pour les communications
    var_fluide.chrono.GetGlobalSeconds(VirtualTimer::COMM,
                                       dt_loc, dt_sum, dt_min, dt_max);

    soundboard.chrono.GetGlobalSeconds(VirtualTimer::COMM,
                                       dt_loc2, dt_sum2, dt_min2, dt_max2);
        
    if (rank_proc == 0)
      {
        cout << "Temps de calcul cumule pour les communications : "
             << dt_sum+dt_sum2 << endl;
        
        if (nb_proc > 1 )
          {
            //&& display_max_time_single_proc)
            cout << "Temps de calcul minimum pour les communications : "
                 << dt_min+dt_min2 << endl;

            cout << "Temps de calcul maximum pour les communications : "
                 << dt_max+dt_max2 << endl;
          }
      }
    
    dt_sum2 += dt_sum;
    
    // temps de calcul global et efficacite parallele
    chrono.GetGlobalSeconds(VirtualTimer::ALL,
                            dt_loc, dt_sum, dt_min, dt_max);
    
    if (rank_proc == 0)
      {
        cout << "Temps de calcul cumule global : " << dt_sum << endl;
        if (nb_proc > 1 )
          cout << "Efficacite parallele " << 100.0*(dt_sum - dt_sum2) / dt_sum << " %" << endl;
      }

    chrono.GetGlobalSeconds(VirtualTimer::FACTO, dt_loc, dt_sum, dt_min, dt_max);
    if (rank_proc == 0)
      cout << "Temps pour Gvol : " << dt_sum << endl;
    
    chrono.GetGlobalSeconds(VirtualTimer::SOLVE, dt_loc, dt_sum, dt_min, dt_max);
    if (rank_proc == 0)
      cout << "Temps pour 2-D : " << dt_sum << endl;

    chrono.GetGlobalSeconds(VirtualTimer::PML, dt_loc, dt_sum, dt_min, dt_max);
    if (rank_proc == 0)
      cout << "Temps pour Advance : " << dt_sum << endl;

#endif
  }
  
  
  //! Avancee du schema en temps (calcul de P^{n+1}, Lambda^{n+1}, etc)
  /*!
    \param[in] nt numero de l'iteration
    \param[in] temps temps t^n
    \param[in] save_snapshot si vrai, des instantanes sont eventuellement ecrits
   */
  void CouplagePiano::Advance(int nt, Real_wp temps, bool save_snapshot)
  {
#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX

#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_fluide.var_vibro.var_volume.var_harmonic.comm_group_mode;
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int rank_proc(0);
#endif

    chrono.Start(VirtualTimer::ALL);
    const MeshNumbering<Dimension2>& mesh_num_plate = var_fluide.var_vibro.var_plate.var_harmonic.GetMeshNumbering(0);
    int nFEM = mesh_num_plate.GetNbDof();
    Real_wp dt = var_fluide.var_vibro.var_volume.GetTimeStep(), invDt = 1.0/dt;
    
    // Sortir les snapshot (air et soundboard) si besoin
    if (save_snapshot)
    {
      chrono.Start(VirtualTimer::OUTPUT);
      soundboard.WriteOutputFEM(nt, temps, var_fluide.var_vibro.var_plate.var_harmonic);
        
      int Nvol = var_fluide.Pn.GetM();
      VectReal_wp Pn;
      if (var_fluide.var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
      {
        Pn.Reallocate(Nvol);
        for (int i = 0; i < Nvol; i++)
          Pn(i) = 0.5*(var_fluide.Pn(i) + var_fluide.Pnm1(i));

        var_fluide.var_vibro.var_volume.WriteSnapshot(nt, temps, Pn, false);
      }
      else
        var_fluide.var_vibro.var_volume.WriteSnapshot(nt, temps, var_fluide.Pn, false);
        
      if (var_fluide.var_vibro.var_volume.var_harmonic.output_rcs_param.GetNbPointsOutside() > 0)
      {
        // on evalue Pn et PnPoint au temps n+1/2
        VectReal_wp PnPoint(Nvol);
        for (int i = 0; i < Nvol; i++)
          PnPoint(i) = (var_fluide.Pn(i) - var_fluide.Pnm1(i))*invDt;
            
        // on appelle le calcul du champ en des points exterieurs au domaine
        var_fluide.var_vibro.var_volume.var_harmonic.output_rcs_param.
        WriteOutput(nt, temps, dt, Pn, PnPoint);
      }
        
      chrono.Stop(VirtualTimer::OUTPUT);
    }
    
    // suffixe utilise pour le debuggage
    //string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
    
    chrono.Start(VirtualTimer::FACTO);
    // Calculer Gvol = P^n-1/2 - dt/2 Mh^-1 (Rh^S V^n - F^n) depuis var_fluide
    VectReal_wp Gvol_minus, Gvol_plus;
    var_fluide.CalculeGvol(invMhSurfMinus, invMhSurfPlus,
                           ChLambdaNm1_Minus, ChLambdaNm1_Plus,
                           Gvol_minus, Gvol_plus, temps);
    
    //DISP(Norm2(Gvol_minus)); DISP(Norm2(Gvol_plus));
    chrono.Stop(VirtualTimer::FACTO);
    chrono.Start(VirtualTimer::SOLVE);
    // on calcule Lambda^n+1/2 en resolvant le complement de Schur associe
    // et on met a jour egalement LambdaPoint^n+1/2
    VectReal_wp rhs, Fvol;
    Real_wp scal_a = soundboard.ComputeScal_a(Gvol_minus, Gvol_plus, rhs, Fvol, var_fluide.var_vibro);
    //~ note.bridge.SetScal_a(scal_a);    

    chrono.Start(VirtualTimer::FLUX);
    if (rank_proc == 0)
    {
      note.bridge.SetVecScal_a(0, scal_a);
        
      // on resout la corde
      note.Advance(temps, nt,
                     save_snapshot, false);        
        
      //DISP(MaxAbs(note.val));
    }
    
    chrono.Stop(VirtualTimer::FLUX);
    
    // on met a jour Lambda et LambdaPoint
    soundboard.ActualiseLambdaSchur(note.bridge.Compute_coef_couplage_V(note,note.val), rhs, Fvol);
    //DISP(Norm2(soundboard.Lambda)); DISP(Norm2(soundboard.LambdaPoint));
    
    // on calcule alors P Lambda^n+1/2
    chrono.Start(VirtualTimer::PROD);
    VectReal_wp Force(nFEM); Force.Fill(0);
    Mlt(soundboard.Pu, soundboard.Lambda, Force);
    chrono.Stop(VirtualTimer::PROD);
    
    chrono.Stop(VirtualTimer::SOLVE); //DISP(Norm2(Force));
    chrono.Start(VirtualTimer::PML);
    // Calculer P^{n+1/2} et V^{n+1}, on stocke Ch Lambda^n+1/2
    var_fluide.Advance(Force, temps, ChLambdaNm1_Minus, ChLambdaNm1_Plus);
    
    //DISP(Norm2(ChLambdaNm1_Minus)); DISP(Norm2(ChLambdaNm1_Plus));
    chrono.Stop(VirtualTimer::ALL);     chrono.Stop(VirtualTimer::PML);
    this->DisplayTimeChronos(nt);
    //int test_input; cout << "we wait " << endl; cin >> test_input;
#endif
  }
  
  
  //! Avancee du schema en temps pour le fluide et la table d'harmonie (pas de corde)
  /*!
    \param[in] nt numero de l'iteration
    \param[in] temps temps t^n
    \param[in] save_snapshot si vrai, des instantanes sont eventuellement ecrits
  */
  void CouplagePiano::AdvancePlateFluid(int nt, Real_wp temps, bool save_snapshot)
  {
#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX
    // suffixe utilise pour le debuggage
    // string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
    
    //soundboard.Lambda.Write("Lambda"+suffix);
    //soundboard.LambdaPoint.Write("LambdaPoint"+suffix);
    //var_fluide.Pn.Write("Pn_"+suffix);
    //var_fluide.Pnm1.Write("Pnm1_"+suffix);
    //var_fluide.Vh_n.Write("Vh_"+suffix);
    
    chrono.Start(VirtualTimer::ALL);
    const MeshNumbering<Dimension2>& mesh_num_plate = var_fluide.var_vibro.var_plate.var_harmonic.GetMeshNumbering(0);
    int nFEM = mesh_num_plate.GetNbDof();
    Real_wp dt = var_fluide.var_vibro.var_volume.GetTimeStep(), invDt = 1.0/dt;
    note.t_courant = temps;
    
    // Sortir les snapshot (air et soundboard) si besoin
    if (save_snapshot)
      {
        chrono.Start(VirtualTimer::OUTPUT);
        soundboard.WriteOutputFEM(nt, temps, var_fluide.var_vibro.var_plate.var_harmonic);
        
        int Nvol = var_fluide.Pn.GetM();
        VectReal_wp Pn;
        if (var_fluide.var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
          {
            Pn.Reallocate(Nvol);
            for (int i = 0; i < Nvol; i++)
              Pn(i) = 0.5*(var_fluide.Pn(i) + var_fluide.Pnm1(i));

            var_fluide.var_vibro.var_volume.WriteSnapshot(nt, temps, Pn, false);
          }
        else
          var_fluide.var_vibro.var_volume.WriteSnapshot(nt, temps, var_fluide.Pn, false);
        
        if (var_fluide.var_vibro.var_volume.var_harmonic.output_rcs_param.GetNbPointsOutside() > 0)
          {
            // on evalue Pn et PnPoint au temps n+1/2
            VectReal_wp PnPoint(Nvol);
            for (int i = 0; i < Nvol; i++)
              PnPoint(i) = (var_fluide.Pn(i) - var_fluide.Pnm1(i))*invDt;
            
            // on appelle le calcul du champ en des points exterieurs au domaine
            var_fluide.var_vibro.var_volume.var_harmonic.output_rcs_param.
              WriteOutput(nt, temps, dt, Pn, PnPoint);
          }
        
        chrono.Stop(VirtualTimer::OUTPUT);               
      }
    
    // Calculer Gvol = P^n-1/2 - dt/2 Mh^-1 (Rh^S V^n - F^n) depuis var_fluide
    VectReal_wp Gvol_minus, Gvol_plus;
    var_fluide.CalculeGvol(invMhSurfMinus, invMhSurfPlus,
                           ChLambdaNm1_Minus, ChLambdaNm1_Plus,
                           Gvol_minus, Gvol_plus, temps);
    
    //Gvol_minus.Write("GvolMinus.dat"); Gvol_plus.Write("GvolPlus.dat");
    // Gvol_minus.Fill(0); Gvol_plus.Fill(0);
    // on calcule Lambda^n+1/2 en resolvant le complement de Schur associe
    // et on met a jour egalement LambdaPoint^n+1/2
    VectReal_wp LambdaNew;
    Real_wp coef = 0.0, pulse = 0.0;
    if (temps < var_fluide.var_vibro.var_plate.GetFinalTimeSource())
      var_fluide.var_vibro.var_plate.SourceOnlyTime(temps, 1, pulse);
    
    //DISP(pulse);
    soundboard.CalculeLambda(coef, pulse, Gvol_minus, Gvol_plus,
			     MhCouplage, var_fluide.var_vibro, LambdaNew);
    
    //DISP(Norm2(soundboard.Lambda));
    //  soundboard.Lambda.Fill(0); soundboard.LambdaPoint.Fill(0);
    // on calcule alors P Lambda^n+1/2
    VectReal_wp Force(nFEM); Force.Fill(0);
    Mlt(soundboard.Pu, soundboard.Lambda, Force);
    
    //Force.Write("Force.dat");
    // Force.Fill(0);
    // Calculer P^{n+1/2} et V^{n+1}, on stocke Ch Lambda^n+1/2
    var_fluide.Advance(Force, temps, ChLambdaNm1_Minus, ChLambdaNm1_Plus);
    
    //ChLambdaNm1_Minus.Write("ChLambdaMinus.dat");
    //ChLambdaNm1_Plus.Write("ChLambdaPlus.dat");
    
    chrono.Stop(VirtualTimer::ALL);    
    this->DisplayTimeChronos(nt);

    //soundboard.Lambda.Write("LambdaNext"+suffix);
    //soundboard.LambdaPoint.Write("LambdaPointNext"+suffix);
    //var_fluide.KhPn.Write("KhPn_"+suffix);
    //var_fluide.Pn.Write("PnNext_"+suffix);
    //var_fluide.Vh_n.Write("VhNext_"+suffix);

    //int test_input; cout << "we wait " << endl; cin >> test_input;
#endif
  }
  

  //! Avancee du schema en temps pour les cordes et la table d'harmonie (pas de fluide)
  /*!
    \param[in] nt numero de l'iteration
    \param[in] temps temps t^n
    \param[in] save_snapshot si vrai, des instantanes sont eventuellement ecrits
   */
  void CouplagePiano::AdvanceStringPlate(int nt, Real_wp temps, bool save_snapshot)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
    int rank_proc(0);
#endif

    chrono_string.Start(CHRONO_FEM);
    if (save_snapshot)
    {
      soundboard.WriteOutputFEM(nt, temps, soundboard.VarSoundboard);
      soundboard.WriteOutputLambda(nt, temps);
    }
    chrono_string.Stop(CHRONO_FEM);
    

    chrono_string.Start(CHRONO_COUPLAGE);
    // On calcule scal_a pour la corde
    Real_wp scal_a =  soundboard.CalculeScal_a();
    //~ note.bridge.SetScal_a(scal_a);
    chrono_string.Stop(CHRONO_COUPLAGE);
    

    // On résout la corde
    chrono_string.Start(CHRONO_NOTE);
    if (rank_proc == 0)
    {
	    note.bridge.SetVecScal_a(0,scal_a);
      note.Advance( temps, nt, save_snapshot);
    }
    
    chrono_string.Stop(CHRONO_NOTE);
    
    // On calcule Lambda et LambdaPoint
    chrono_string.Start(CHRONO_PLATE);
    soundboard.ActualiseLambda(note.bridge.Compute_coef_couplage_V(note,note.val), 0.0);
    chrono_string.Stop(CHRONO_PLATE);  
  }
  
}

#define MONTJOIE_FILE_COUPLAGE_PIANO_CXX
#endif
