#ifndef MONTJOIE_FILE_MULTI_STRING_CXX

namespace Montjoie
{
  
  //! default constructor
  template<int nb_base, int nb_quad, class TypeEquation>
  MultiString<nb_base, nb_quad, TypeEquation>::MultiString()
    : MultiString_Base(dimension), shank(hammer)
  {
    t_begin = t_end = t_courant = 0.0;
    t_display   = 1e-3;
    nb_max_iter   = 0;
    t_begin_interp  = t_end_interp  = deltat_interp = 0.0;
    t_begin_sismo = t_end_sismo = deltat_sismo  = 0.0;
    t_begin_ener  = t_end_ener  = deltat_ener = 0.0;
    t_begin_crush = t_end_crush = deltat_crush  = 0.0;
  
    nb_max_iter = -1;
    Deltat      = 1e-5;
    Deltat2     = square(Deltat);
    invDeltat   = 1.0/Deltat;
    invDeltat2    = square(invDeltat);
    invDeuxDeltat = 0.5*invDeltat;
    
    somme_interaction_hammer = 0.0;

    size      = 0;
    size_LM     = 0;
    nb_lagrange_mult = 1;
    nb_corde    = 0;
    nb_corde_struck = 0;
    nb_iter_courant = 0;
    
    size_buffer   = 100;
    size_buffer_ener= 100;
    
    nb_combinaisons = 0;
    permutations  = PERMUT_ALL;
    
    solver = THETA_NL_HALF;
    
    existence_donnee_ini  = false; 
    type_initial_condition  = NONE;
    amplitude_ini     = 0.0;
    nb_mode_ini       = 0;
    dimension_ini     = 0;
    
    // name_data_file = name;
    
    load_reprise = false; save_reprise = false;
    load_iter_reprise = 0; save_iter_reprise = 0;
    bascule_reprise = -1;
    path_reprise = "./";
    file_storage_param = "param.txt";

    shank.is_interacting_string = false;
    offset_hammer = 0;
    offset_hammer_aux = 0;
    offset_shank = 0;
    offset_force_hammer = 0;
    
    print_level = 0;
    beginning_time_force = -1.0; last_time_force = -2.0;
    final_time_shank = -1.0;
    
    drop_hammer = false;
    drop_hammer_when_possible = true;
    //Newton_solver.type_matrix_lu = Newton_solver.LU_TINY_ARROW;
    //Newton_solver.type_matrix_lu = Newton_solver.LU_ARROW;
    Newton_solver.type_matrix_lu = Newton_solver.LU_SPARSE;
    use_naive_non_linear_term = false;

    separate_energy_uv = false;
    shank_present_from_t0 = false;
    output_sismo_binary = false;
    output_sismo_double = false;
    
    number_of_duplex_string = 0;
    vect_duplex_string_ID = 0;    
    
    BridgeHasHeight = false;
    useQuaternion = false;
    
    file_output_LM = "Sismo_LM.txt";
    file_output_hammer = "Sismo_Hammer.txt";
    file_output_force = "ForceMarteau.txt";
    file_output_energy = "EnergyNote.txt";    
  }
  
  
  //! setting parameters with a line of the data file
  /*!
    \param[in] keyword keyword of the entry
    \param[in] param parameters associated with the keyword
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>
  ::SetInputData(const string & keyword, const Vector<string> &param)
  {
    bridge.SetInputData(keyword, param);
    
    if (!keyword.compare("TimeInterval"))
    {
      if (param.GetM() <= 1)
      {
        cout << "In SetInputData of MultiString" << endl;
        cout << "TimeInterval needs two parameters, for instance :" << endl;
        cout << "TimeInterval = t0 tf" << endl;
        cout << "Current parameters are : " << endl << param << endl;
        abort();
      }

        t_begin = to_num<Real_wp>(param(0));
        t_end = to_num<Real_wp>(param(1));
      }
    else if (!keyword.compare("DisplayFrequency"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "DisplayFrequency needs one parameter, for instance :" << endl;
      cout << "DisplayFrequency = dt" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }
  
        t_display = to_num<Real_wp>(param(0));
      }
    else if (!keyword.compare("DirectoryOutputString"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "DirectoryOutputString needs one parameter, for instance :" << endl;
      cout << "DirectoryOutputString = chemin" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }
  
  DOSSIER = param(0);
      }
    else if (!keyword.compare("StringInterpTimeGrid"))
      {
  if (param.GetM() <= 2)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "StringInterpTimeGrid needs three parameters, for instance :" << endl;
      cout << "StringInterpTimeGrid = t0 tf dt" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }
  
        t_begin_interp = to_num<Real_wp>(param(0));
        t_end_interp = to_num<Real_wp>(param(1));
        deltat_interp = to_num<Real_wp>(param(2));
  
        if (deltat_interp < Deltat)
          cout << "Warning : your time step for interpolation cannot be respected." <<endl;        
      }
    else if (!keyword.compare("StringSismoTimeGrid"))
      {
  if (param.GetM() <= 2)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "StringSismoTimeGrid needs three parameters, for instance :" << endl;
      cout << "StringSismoTimeGrid = t0 tf dt" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }
  
        t_begin_sismo = to_num<Real_wp>(param(0));
        t_end_sismo = to_num<Real_wp>(param(1));
        deltat_sismo = to_num<Real_wp>(param(2));
        if (param.GetM() > 3)
          size_buffer = to_num<int>(param(3));
        
        if (deltat_sismo < Deltat)
          cout << "Warning : your time step for sismo cannot be respected." <<endl;        
      }    
    else if (!keyword.compare("Energy"))
      {
  if (param.GetM() <= 2)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "Energy needs three parameters, for instance :" << endl;
      cout << "Energy = t0 tf dt" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        t_begin_ener = to_num<Real_wp>(param(0));
        t_end_ener = to_num<Real_wp>(param(1));
        deltat_ener = to_num<Real_wp>(param(2));
        if(param.GetM() > 3)
          size_buffer_ener = to_num<int>(param(3));
    
        if (param.GetM() > 4)
          {
            if (param(4) == "UV")
              separate_energy_uv = true;
            else
              separate_energy_uv = false;
          }
        
        if (deltat_ener < Deltat)
          cout << "Warning : your time step for energy cannot be respected." <<endl;
        
      }
    else if (!keyword.compare("SismoCrushTimeGrid"))
      {
  if (param.GetM() <= 2)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "SismoCrushTimeGrid needs three parameters, for instance :" << endl;
      cout << "SismoCrushTimeGrid = t0 tf dt" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        t_begin_crush = to_num<Real_wp>(param(0));
        t_end_crush = to_num<Real_wp>(param(1));
        deltat_crush = to_num<Real_wp>(param(2));
        if (param.GetM() > 3)
          size_buffer = to_num<int>(param(3));
        
        if (deltat_crush < Deltat)
          cout << "Warning : your time step for sismo crush cannot be respected." <<endl;        
      }   
    else if (!keyword.compare("TimeStep"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "TimeStep needs one parameter, for instance :" << endl;
      cout << "TimeStep = dt" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        Deltat = to_num<Real_wp>(param(0));
        Deltat2 = square(Deltat);
        invDeltat = 1.0/Deltat;
        invDeltat2 = square(invDeltat);
        invDeuxDeltat = 0.5*invDeltat;
      }
    else if (!keyword.compare("RemoveHammerAfterForce"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "RemoveHammerAfterForce needs one parameter, for instance :" << endl;
      cout << "RemoveHammerAfterForce = NO" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        if (param(0) == "NO")
          drop_hammer_when_possible = false;
      }
    else if (!keyword.compare("HammerDamping"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "HammerDamping needs one parameter, for instance :" << endl;
      cout << "HammerDamping = amo" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        // cout << "Hammer Damping" << param(0) << endl;
        param_corde_data.InitHammerDamping(param(0));
      }
    else if (!keyword.compare("PlanCordes"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "PlanCordes needs one parameter, for instance :" << endl;
      cout << "PlanCordes = fichier" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        param_corde_data.InitPlanCordes(param(0));
      }
    else if (!keyword.compare("StringDamping"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "StringDamping needs three parameters, for instance :" << endl;
      cout << "StringDamping = amo" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        //   cout << "String Damping" << param(0) << endl;
        param_corde_data.InitStringDamping(param(0));
      }
    else if (!keyword.compare("String"))
      {
        PianoString<nb_base, nb_quad, TypeEquation>
          piano_string_temp(param, hammer, param_corde_data);

        piano_strings.PushBack(piano_string_temp);
      }
    else if (!keyword.compare("DuplexString"))
      {
    if (param.GetM() <= 1)
      {
            cout << "In SetInputData of MultiString" << endl;
            cout << "DuplexString needs at least 2 parameters, for instance :" << endl;
            cout << "DuplexString = number_of_duplex_string duplex_string_ID" << endl;
            cout << "Current parameters are : " << endl << param << endl;
            abort();        
      }
      
      number_of_duplex_string = to_num<int>(param(0));
       cout << " we have set number of duplex string" << endl;
      if (param.GetM() <= number_of_duplex_string)
      {
            cout << "In SetInputData of MultiString" << endl;
            cout << "You specified " << number_of_duplex_string << " duplex string" << endl;
            cout << "but only gave " << param.GetM()-1 << " of them" << endl;
            cout << "Current parameters are : " << endl << param << endl;
            abort();        
      }   
      
    vect_duplex_string_ID.Reallocate(number_of_duplex_string);
        
    for (int j = 0; j < number_of_duplex_string; j++)
      {
      vect_duplex_string_ID(j) = to_num<int>(param(j+1));
        }
    
    }
    else if (!keyword.compare("LumpedBridgeHeight"))
      {
        if (param.GetM() != 3)
          {
            cout << "In SetInputData of MultiString" << endl;
            cout << "LumpedBridgeHeight needs 3 parameters, for instance :" << endl;
            cout << "LumpedBridgeHeight = height0 height1 height2" << endl;
            cout << "Current parameters are : " << endl << param << endl;
            abort();
          }
        
        BridgeHasHeight = true;
        
        for (int j = 0; j < param.GetM(); j++)
          {
            for (int i = 0; i < piano_strings.GetM(); i++)
              piano_strings(i).vect_bridge_height(j) = to_num<Real_wp>(param(j));  
          }       
      }
    else if (!keyword.compare("UseQuaternion"))
    {
    if (!param(0).compare("TRUE"))
      useQuaternion = true; 
  }
    else if (!keyword.compare("ExtraBridge"))
      {
        if(GetStringType() == 0 ||  GetStringType() == 2 || GetStringType() == 6 || GetStringType() == 7)
          {   
            if (!param(0).compare("TRUE"))
              {
                bridge.extra_bridge = true; 
    
              }  
          }
        else
          {
            cout << "ExtraBridge is only available for LIN, NL and StiffNL2T" << endl;
            cout << "There is no need to use ExtraBridge if you intend to solve" << endl;
            cout << "a string with stiffness. use MatchNLaMuToStringDOF = TRUE instead" << endl;
            cout << "to solve the rotation via the boundary conditions and LM." << endl;
            abort();
          } 
      }
    else if (!keyword.compare("MatchNLaMuToStringDOF"))
      {
        if(!param(0).compare("TRUE"))
          {
            bridge.matchNLaMuToStringDOF = true;
          } 
      }
    else if (!keyword.compare("Bridge"))
      {
        if (param.GetM() <= 0)
          {
            cout << "In SetInputData of MultiString" << endl;
            cout << "Bridge needs one parameter, for instance :" << endl;
            cout << "Bridge = ATTACHED" << endl;
            cout << "Current parameters are : " << endl << param << endl;
            abort();
          }
        
        if (!param(0).compare("ATTACHED") )
          {
            for (int i = 0; i < piano_strings.GetM(); i++)
              {
                piano_strings(i).is_attached = true;
                piano_strings(i).angle_chevalet = 0.0;
                cout << "String is attached. NO BRIDGE." << endl;
              }
          }
        else
          {           
            if (param.GetM() != piano_strings.GetM())
              {
                cout << "Please specify an angle for each string. "<< endl;
                DISP(param);
                abort();
              }
            
            for (int i = 0; i < param.GetM(); i++)
              {
                piano_strings(i).angle_chevalet = to_num<Real_wp>(param(i));
    piano_strings(i).is_attached = false;
              }            
          }
      }
    else if (!keyword.compare("BridgeStringPolarisationAngle"))
      {
    if (param.GetM() != piano_strings.GetM())
      {
      cout << "Please specify the angle offset of the polarisation for each string. "<< endl;
      DISP(param);
      abort();
      }
      
    for (int i = 0; i < param.GetM(); i++)
      {
      piano_strings(i).angleStringPolarisation = to_num<Real_wp>(param(i));
      }         
    }
    else if ( (!keyword.compare("Hammer")) || (!keyword.compare("HammerStiffness")))
      {
        if (piano_strings.GetM() <= 0)
          {
            cout << "The hammer requires at least one string" << endl;
            abort();
          }
        
        hammer.SetInputData(keyword, param, piano_strings(0), param_corde_data);
        

    if (hammer.strike_at_an_angle)
    {
      if(piano_strings(0).var_eq.type_equation == TypeEquation::WAVE_NL_2T_TIMO)  
        cout << "Hammer with angle is used with double polarised strings. OK!" << endl;
      else
      {
        cout << "Hammer with angle should be used with double polarised strings only" << endl;
        abort();
      }
    }
      }
    else if (keyword == "PrintLevel")
      {
        if (param.GetM() <= 0)
          {
            cout << "In SetInputData of MultiString" << endl;
            cout << "PrintLevel needs one parameter, for instance :" << endl;
            cout << "PrintLevel = 2" << endl;
            cout << "Current parameters are : " << endl << param << endl;
            abort();
          }
        
        print_level = to_num<int>(param(0));
      }
    else if (shank.AllowedKeyword(keyword))
      {
  shank.is_interacting_string = true;
  shank.SetInputData(keyword, param);
      }
    else if (!keyword.compare("NewtonParam"))
      {
        Newton_solver.SetInputData(keyword, param);
      }
    else if (!keyword.compare("OutputDiffMatrix"))
      {
        Newton_solver.SetInputData(keyword, param);
      }
    else if (!keyword.compare("LinearSolverNewton"))
      {
        if (param(0) == "Arrow")
          Newton_solver.type_matrix_lu = Newton_solver.LU_ARROW;
        else if (param(0) == "TinyArrow")
          Newton_solver.type_matrix_lu = Newton_solver.LU_TINY_ARROW;
        else
          Newton_solver.type_matrix_lu = Newton_solver.LU_SPARSE;
      }
    else if (!keyword.compare("ForceReevaluationJacobian"))
      {
        Newton_solver.SetInputData(keyword, param);
      }    
    else if (!keyword.compare("StringInterp"))
      {
        if (param.GetM() > piano_strings.GetM())
          {
            cout << "Too many StringInterp required ! "<< endl;
            abort();
    }
        
        for (int i=0; i<param.GetM(); i++)
          piano_strings(i).Nx_interp = to_num<int>(param(i));
    
      }
    else if (!keyword.compare("StringSismo"))
      {
        if (param.GetM() > piano_strings.GetM())
          {
            cout << "Too many StringSismo required ! "<< endl;
            abort();
    }
        
        for (int i = 0; i < param.GetM(); i++)
          {
            piano_strings(i).point_observation = to_num<Real_wp>(param(i));
    }
      }    
    else if (!keyword.compare("UseLobattoQuadrature"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "UseLobattoQuadrature needs one parameter, for instance :" << endl;
      cout << "UseLobattoQuadrature = YES" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        bool lob = false;
        if (param(0) == "YES")
          lob = true;
        
        for (int i = 0; i < param.GetM(); i++)
          piano_strings(i).use_lobatto_quadrature = lob;
        
      }    
    else if (!keyword.compare("StringSolver"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "StringSolver needs at least one parameter, for instance :" << endl;
      cout << "StringSolver = EXPLICIT" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        if (!param(0).compare("EXPLICIT"))
    {
            solver = EXPLICIT;
            cout << "solve explicit not implemented yet." << endl;
            abort();
    }
    
        if (!param(0).compare("THETA_NL_HALF"))
    {
            solver  = THETA_NL_HALF;
            theta.Reallocate(param.GetM()-1);
            for (int p = 0; p < theta.GetM(); p++)
              theta(p)  = to_num<Real_wp>(param(p+1));
            
    }
        if (!param(0).compare("THETA_NL_QUARTER"))
    {
            solver  = THETA_NL_QUARTER;
            theta.Reallocate(param.GetM()-1);
            for (int p = 0; p < theta.GetM(); p++)
              theta(p)  = to_num<Real_wp>(param(p+1));
            
    }
      }    
    else if (!keyword.compare("StringPermutations"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "StringPermutations needs one parameter, for instance :" << endl;
      cout << "StringPermutations = PERMUT_TWO_WAY" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        if (!param(0).compare("PERMUT_ONE"))
    {
            permutations = PERMUT_ONE;
    }
        if (!param(0).compare("PERMUT_TWO_WAY"))
    {
            permutations = PERMUT_TWO_WAY;
    }
        if (!param(0).compare("PERMUT_ALL"))
    {
            permutations = PERMUT_ALL;
    }
      }    
    else if (!keyword.compare("NaiveNonLinearScheme"))
      {
        if (param(0) == "YES")
          use_naive_non_linear_term = true;
        else
          use_naive_non_linear_term = false;
      }
    else if (!keyword.compare("StringInitialData"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "StringInitialData needs one parameter, for instance :" << endl;
      cout << "StringInitialData = NONE" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        if (!param(0).compare("NONE"))
    {
            existence_donnee_ini   = false;
            type_initial_condition = NONE;
    }
        if (!param(0).compare("SINUS") || !param(0).compare("DOUBLE_SINUS"))
    {
      if (param.GetM() <= 3)
        {
    cout << "In SetInputData of MultiString" << endl;
    cout << "StringInitialData needs 4 parameters, for instance :" << endl;
    cout << "StringInitialData = SINUS amplitude mode_num dimension" << endl;
    cout << "Current parameters are : " << endl << param << endl;
    abort();
        }

            existence_donnee_ini   = true;
      if(!param(0).compare("SINUS"))
        type_initial_condition = SINUS;
      else 
        type_initial_condition = DOUBLE_SINUS;
            amplitude_ini = to_num<Real_wp>(param(1));
            nb_mode_ini   = to_num<int>(param(2));
            dimension_ini = to_num<int>(param(3));
            
    }
        if (!param(0).compare("GAUSSIAN"))
    {
            existence_donnee_ini   = true;
            type_initial_condition = GAUSSIAN;
            cout << "Not implemented yet. Feel free to do it !" <<endl;
            abort();

    }
  if (!param(0).compare("SINUSMULTI"))
    {
      
        
      if ((param.GetM()-1) - 3*dimension != 0)
      {
        cout << "You missed out some data in SINUSMULTI" << endl;

        cout << "Right now, your total inputs are " << param.GetM()-1 << endl;
        cout << "So, you should have " << dimension*3 << " inputs. " << endl;
        
        cout << "Example: if your dimension == 2, you should have " << endl;
        cout << "SINUSMULTI = a0 n0 d0 a1 n1 d1 " << endl;
        abort();
      }
      
            existence_donnee_ini   = true;
            type_initial_condition = SINUSMULTI;
      vect_amplitude_ini.Reallocate(dimension);    
      vect_nb_mode_ini.Reallocate(dimension);  
      vect_dimension_ini.Reallocate(dimension);
      vect_amplitude_ini.Fill(0);    
      vect_nb_mode_ini.Fill(0);  
      vect_dimension_ini.Fill(0); 
          
      for (int i = 0; i < dimension; i++)
      {            
      vect_amplitude_ini(i) = to_num<Real_wp>(param(1+i*3));
            vect_nb_mode_ini(i)   = to_num<int>(param(2+i*3));
            vect_dimension_ini(i) = to_num<int>(param(3+i*3));
      }
      
      amplitude_ini = vect_amplitude_ini(0);
      nb_mode_ini = vect_nb_mode_ini(0);
      dimension_ini = vect_dimension_ini(0);
    }
  if (!param(0).compare("NONZERO"))
    {
            existence_donnee_ini   = true;
            type_initial_condition = NONZERO;
      
    }   
  if (!param(0).compare("PLUCK"))
    {
            existence_donnee_ini   = true;
            type_initial_condition = PLUCK;
            pluck_position = to_num<Real_wp>(param(1));
            amplitude_ini = to_num<Real_wp>(param(2));
      dimension_ini = to_num<Real_wp>(param(3));
      
    } 
  if (!param(0).compare("PLUCKANGLE"))
    {
            existence_donnee_ini   = true;
            type_initial_condition = PLUCKANGLE;
            pluck_position = to_num<Real_wp>(param(1));
            amplitude_ini = to_num<Real_wp>(param(2));
      dimension_ini = to_num<Real_wp>(param(3));
      pluckangle_ini = to_num<Real_wp>(param(4));
      
    }         
      }    
    else if (keyword == "OutputFormatSismo")
      {
  if (param.GetM() <= 1)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "OutputFormatSismo needs two parameters, for instance :" << endl;
      cout << "OutputFormatSismo = FLOAT BINARY" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

        if (param(0) == "DOUBLE")
          output_sismo_double = true;
        else
          output_sismo_double = false;

        if (param(1) == "BINARY")
          output_sismo_binary = true;
        else
          output_sismo_binary = false;
          
      }
    else if (!keyword.compare("LoadReprise"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "LoadReprise needs one parameter, for instance :" << endl;
      cout << "LoadReprise = YES" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

  if (param(0) == "YES")
          load_reprise = true;
        else
          load_reprise = false;
      }
    else if (!keyword.compare("SaveReprise"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "SaveReprise needs one parameter, for instance :" << endl;
      cout << "SaveReprise = YES" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }

  save_reprise = true;
  save_iter_reprise = to_num<int>(param(0));
      }
    else if (!keyword.compare("PathReprise"))
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "PathReprise needs one parameter, for instance :" << endl;
      cout << "PathReprise = path" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }
  
        path_reprise = param(0);
      }        
    else if (keyword == "FileOutputString")
      {
  if (param.GetM() <= 3)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "FileOutputString needs four parameters, for instance :" << endl;
      cout << "FileOutputString = file_LM file_hammer file_force file_energy" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }
        
        file_output_LM = param(0);
        file_output_hammer = param(1);
        file_output_force = param(2);
        file_output_energy = param(3);
      }
    else if (keyword == "FileOutputSismoString")
      {
  if (param.GetM() <= 3)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "FileOutputSismoString needs four parameters, for instance :" << endl;
      cout << "FileOutputSismoString = string0 string1 string2 string3" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }
        
        for (int i = 0; i < piano_strings.GetM(); i++)
          {
            piano_strings(i).file_output_sismo_string = param(0);
            piano_strings(i).file_output_last_point = param(1);
            piano_strings(i).file_output_crush = param(2);
            piano_strings(i).file_output_interp = param(3);
          }
      }
      else if (keyword == "HammerWithAuxiliaryVariable")
      {
  if (param.GetM() <= 0)
    {
      cout << "In SetInputData of MultiString" << endl;
      cout << "HammerWithAuxiliaryVariable needs one parameters, for instance :" << endl;
      cout << "HammerWithAuxiliaryVariable = YES" << endl;
      cout << "Current parameters are : " << endl << param << endl;
      abort();
    }
        if (param(0) == "YES")
        {
          hammer.with_auxiliary_variable = true;
        }
        else
        {
          hammer.with_auxiliary_variable = false;
        }
        
      }
  }
  
  
  //! construction of arrays needed for time iterations
  /*!
    \param[in] num numero du dossier dans le cas d'une reprise, 0 sinon
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>::ConstructAll(int num)
  {
    // construction du numero de dossier
    // si le dossier n'est pas connu
    if (DOSSIER.size() == 0)
      {
  DOSSIER = "[STIFFOUT]/"; 
  EcritDossier(DOSSIER, name_data_file, num);
      }
    
    // on initialise les sorties de type sismogrammes
    bool remove_file = true;
    if (num >= 1)
      remove_file = false;
    
    output_LM.SetDoublePrecision(output_sismo_double);
    output_hammer.SetDoublePrecision(output_sismo_double);
    output_force.SetDoublePrecision(output_sismo_double);
    output_ener.SetDoublePrecision(output_sismo_double);

    output_LM.SetBinary(output_sismo_binary);
    output_hammer.SetBinary(output_sismo_binary);
    output_force.SetBinary(output_sismo_binary);
    output_ener.SetBinary(output_sismo_binary);
    
    output_LM.    Init( DOSSIER + file_output_LM, size_buffer, remove_file);
    output_hammer.Init( DOSSIER + file_output_hammer, size_buffer, remove_file);
    output_force. Init( DOSSIER + file_output_force, size_buffer, remove_file);
    output_ener.  Init( DOSSIER + file_output_energy, size_buffer_ener, remove_file);
    
    // LES PERMUTATIONS
    nb_corde = piano_strings.GetM();
    GenerateCombinaisons(*this);
    // fin permutations
    
    // on regarde si on a reellement besoin du manche et du marteau
    shank_present_from_t0 = shank.is_interacting_string;
    if (load_reprise)
      {
        int n0 = load_iter_reprise;
        Real_wp t0 = t_begin + n0*Deltat;
        
        if (final_time_shank > 0)
          if (t0 > final_time_shank)
            {
              // on a depasse le moment de contact entre le marteau et la corde
              shank.is_interacting_string = false;
              hammer.is_interacting = false;
              drop_hammer = true;
            }
      }
    
    if (hammer.type_interaction == hammer.NONE)
      drop_hammer = true;
      
      

    bool struck_one = false;
    for (int i = 0; i < piano_strings.GetM(); i++)
      if (piano_strings(i).is_struck)
        struck_one = true;
    
    if (!struck_one)
      drop_hammer = true;

    // Conditions de transmission avec le chevalet
    for (int i = 0; i < nb_corde; i++)
      {
      
    piano_strings(i).transmission_chevalet.SetIdentity();
        Real_wp angle_chevalet = piano_strings(i).angle_chevalet;
        Real_wp angleStringPolarisation = piano_strings(i).angleStringPolarisation;
        
        if (piano_strings(i).var_eq.type_equation == TypeEquation::WAVE_LIN )
          {
            piano_strings(i).transmission_chevalet(0,0) = cos(angle_chevalet*pi_wp/180.0);

          }
        else if (piano_strings(i).var_eq.type_equation == TypeEquation::WAVE_LIN_TIMO)
          {
            piano_strings(i).transmission_chevalet(0,0) = cos(angle_chevalet*pi_wp/180.0);
            piano_strings(i).transmission_chevalet(0,1) = 0.0;
            piano_strings(i).transmission_chevalet(1,0) = 0.0;
            piano_strings(i).transmission_chevalet(1,1) = -sin(angle_chevalet*pi_wp/180.0);
            
      if (!piano_strings(i).is_attached)
              piano_strings(i).transmission_chevalet(0,1) = -GetBridgeHeight(2);
      
    }
        else if (piano_strings(i).var_eq.type_equation == TypeEquation::WAVE_NL)
          {
            piano_strings(i).transmission_chevalet(0,0) = cos(angle_chevalet*pi_wp/180.0);
            piano_strings(i).transmission_chevalet(0,1) = sin(angle_chevalet*pi_wp/180.0);
            
            piano_strings(i).transmission_chevalet(1,0) = -sin(angle_chevalet*pi_wp/180.0);
            piano_strings(i).transmission_chevalet(1,1) = cos(angle_chevalet*pi_wp/180.0);
    }
  else if (piano_strings(i).var_eq.type_equation == TypeEquation::WAVE_NL_2T)
          {
            piano_strings(i).transmission_chevalet(0,0) = cos(angle_chevalet*pi_wp/180.0);
            piano_strings(i).transmission_chevalet(0,1) = 0;
            piano_strings(i).transmission_chevalet(0,2) = sin(angle_chevalet*pi_wp/180.0);
            
      piano_strings(i).transmission_chevalet(1,0) = 0;
            piano_strings(i).transmission_chevalet(1,1) = 1;
            piano_strings(i).transmission_chevalet(1,2) = 0;
            
      piano_strings(i).transmission_chevalet(2,0) = -sin(angle_chevalet*pi_wp/180.0);
      piano_strings(i).transmission_chevalet(2,1) = 0;
            piano_strings(i).transmission_chevalet(2,2) = cos(angle_chevalet*pi_wp/180.0);
            
    }
        else if (piano_strings(i).var_eq.type_equation == TypeEquation::WAVE_NL_TIMO)
          {
            if (!piano_strings(i).is_attached)  
            {   
                piano_strings(i).transmission_chevalet(0,0) = cos(angle_chevalet*pi_wp/180.0);
                piano_strings(i).transmission_chevalet(0,1) = sin(angle_chevalet*pi_wp/180.0);
                piano_strings(i).transmission_chevalet(0,2) = 0; 
    
                piano_strings(i).transmission_chevalet(1,0) = -sin(angle_chevalet*pi_wp/180.0);
                piano_strings(i).transmission_chevalet(1,1) = cos(angle_chevalet*pi_wp/180.0);
                piano_strings(i).transmission_chevalet(1,2) = 0; 
                
        piano_strings(i).transmission_chevalet(2,0) = 0;
        piano_strings(i).transmission_chevalet(2,1) = 0;
        piano_strings(i).transmission_chevalet(2,2) = 0;                
                
            if (BridgeHasHeight)  
              { 
        piano_strings(i).transmission_chevalet(0,2) = -GetBridgeHeight(2);  
        
        piano_strings(i).transmission_chevalet(1,2) = -GetBridgeHeight(0);   
                      
        piano_strings(i).transmission_chevalet(2,0) = 0;
        piano_strings(i).transmission_chevalet(2,1) = 0;
        piano_strings(i).transmission_chevalet(2,2) = 1;
              }
              
        }
        else
        {
                piano_strings(i).transmission_chevalet(0,0) = 1;
                piano_strings(i).transmission_chevalet(0,1) = 0;
                piano_strings(i).transmission_chevalet(0,2) = 0; 
    
                piano_strings(i).transmission_chevalet(1,0) = 0;
                piano_strings(i).transmission_chevalet(1,1) = 1;
                piano_strings(i).transmission_chevalet(1,2) = 0; 
                
        piano_strings(i).transmission_chevalet(2,0) = 0;
        piano_strings(i).transmission_chevalet(2,1) = 0;
        piano_strings(i).transmission_chevalet(2,2) = 0;        
      }
         
    }
      else if (piano_strings(i).var_eq.type_equation == TypeEquation::WAVE_NL_2T_TIMO)
        {
            if (useQuaternion) //if you use angleStringPolarisation, angle_chevalet (and angle_chevalet_lateral)
                ComputeTransmissionChevaletQuaternion(i,angleStringPolarisation,angle_chevalet); 
            else
            {
        piano_strings(i).transmission_chevalet(0,0) = cos(angleStringPolarisation*pi_wp/180.0);
        piano_strings(i).transmission_chevalet(0,1) = -sin(angleStringPolarisation*pi_wp/180.0);
        piano_strings(i).transmission_chevalet(0,2) = 0;
        piano_strings(i).transmission_chevalet(0,3) = 0;
        piano_strings(i).transmission_chevalet(0,4) = 0;
        
        piano_strings(i).transmission_chevalet(1,0) = sin(angleStringPolarisation*pi_wp/180.0);
        piano_strings(i).transmission_chevalet(1,1) = cos(angleStringPolarisation*pi_wp/180.0);
        piano_strings(i).transmission_chevalet(1,2) = 0;
        piano_strings(i).transmission_chevalet(1,3) = 0;
        piano_strings(i).transmission_chevalet(1,4) = 0;
        
        piano_strings(i).transmission_chevalet(2,0) = 0;
        piano_strings(i).transmission_chevalet(2,1) = 0;
        piano_strings(i).transmission_chevalet(2,2) = 1;
        piano_strings(i).transmission_chevalet(2,3) = 0;  
        piano_strings(i).transmission_chevalet(2,4) = 0; 
         
        piano_strings(i).transmission_chevalet(3,0) = 0;
        piano_strings(i).transmission_chevalet(3,1) = 0;
        piano_strings(i).transmission_chevalet(3,2) = 0;
        piano_strings(i).transmission_chevalet(3,3) = 0;
        piano_strings(i).transmission_chevalet(3,4) = 1;
        
        piano_strings(i).transmission_chevalet(4,0) = 0;
        piano_strings(i).transmission_chevalet(4,1) = 0;
        piano_strings(i).transmission_chevalet(4,2) = 0;
        piano_strings(i).transmission_chevalet(4,3) = 1;
        piano_strings(i).transmission_chevalet(4,4) = 0;  
        
        if (BridgeHasHeight)
        {
          piano_strings(i).transmission_chevalet(0,3) = -GetBridgeHeight(2);
          piano_strings(i).transmission_chevalet(1,4) = GetBridgeHeight(2);
          piano_strings(i).transmission_chevalet(2,3) = GetBridgeHeight(0)*cos(angleStringPolarisation*pi_wp/180.0)
                               -GetBridgeHeight(1)*sin(angleStringPolarisation*pi_wp/180.0);  
          piano_strings(i).transmission_chevalet(2,4) =-GetBridgeHeight(1)*cos(angleStringPolarisation*pi_wp/180.0)
                               -GetBridgeHeight(0)*sin(angleStringPolarisation*pi_wp/180.0);        
          piano_strings(i).transmission_chevalet(3,3) = sin(angleStringPolarisation*pi_wp/180.0);
          piano_strings(i).transmission_chevalet(3,4) = cos(angleStringPolarisation*pi_wp/180.0);
          piano_strings(i).transmission_chevalet(4,3) = cos(angleStringPolarisation*pi_wp/180.0);
          piano_strings(i).transmission_chevalet(4,4) = -sin(angleStringPolarisation*pi_wp/180.0);  
        }
      }
      }
      else if (piano_strings(i).var_eq.type_equation == TypeEquation::WAVE_LIN_2T)
        {
            piano_strings(i).transmission_chevalet(0,0) = cos(angleStringPolarisation*pi_wp/180.0);
            piano_strings(i).transmission_chevalet(0,1) = -sin(angleStringPolarisation*pi_wp/180.0);
      
            piano_strings(i).transmission_chevalet(1,0) = sin(angleStringPolarisation*pi_wp/180.0);
            piano_strings(i).transmission_chevalet(1,1) = cos(angleStringPolarisation*pi_wp/180.0);

      }     
    else
      {
        cout << "Equation type not found." << endl;
        abort();
      }
  
  

  
  
  // setting transmission_chevalet from nu and nuorth
  //NU takes row 0, NUORTH row 1, NUORTH2 row 2

  //~ for (int k = 0; k < dimension; k++)
    //~ {
      //~ piano_strings(i).transmission_chevalet(0, k) = piano_strings(i).NU(k);
      //~ if (dimension > 1)
        //~ piano_strings(i).transmission_chevalet(1, k) = piano_strings(i).NUORTH(k);
      //~ if (dimension > 2)
        //~ piano_strings(i).transmission_chevalet(2, k) = piano_strings(i).NUORTH2(k);       
    //~ }
      }

          
    // we construct here all what is needed for the bridge
    if (piano_strings.GetM() <= 0)
      {
  cout << "No strings are present, the computation cannot continue" << endl;
  abort();
      }
 
    bridge.ConstructAll(*this, piano_strings(0).var_eq, DOSSIER);   
      
  // we set duplex strings here
  if (piano_strings.GetM() <= number_of_duplex_string)
  {
    cout << "We have only " << piano_strings.GetM() << " strings" << endl;
    cout << "We have specified " << number_of_duplex_string << " duplex strings" << endl;
    cout << "We need at least one string that is not duplex" << endl;
    abort();
  }
  
  for (int i = 0; i < nb_corde; i++)
    {
    for (int j = 0; j < vect_duplex_string_ID.GetM(); j++)
      {
      if (vect_duplex_string_ID(j)==i)    
        {
          piano_strings(i).is_duplex = true;  
          cout << " string " << i << " is duplex" << endl;
        }
      }
    
    }      


    // fin transmission 
    if (piano_strings(0).var_eq.only_one_LM || piano_strings(0).is_attached)
      {
        nb_lagrange_mult = 1;
        bridge.SetNbLagrangeMultipliers(nb_lagrange_mult);
      }
    else
      nb_lagrange_mult = bridge.GetNbLagrangeMultipliers();
    
    // on incremente au fur et a mesure size, en le stockant
    // pour avoir les offset pour chaque inconnue
    
    // Les cordes
    nb_dof_before.Reallocate(nb_corde+1);
    nb_dof_before(0) = 0;
    for(int i = 0; i < nb_corde; i++)
      {
        if (piano_strings(i).Nx <= 1)
          {
            cout << "The string must contain at least two points" << endl;
            abort();
          }

        piano_strings(i).ConstructMesh(piano_strings(i).L);
        nb_dof_before(i+1) = nb_dof_before(i) + piano_strings(i).mesh.GetNbDof();
        piano_strings(i).ConstruitImpact(hammer);
        piano_strings(i).ConstruitObserv();
        piano_strings(i).ConstruitInterp(piano_strings(i).L);
        piano_strings(i).numero_string = i;
        if (piano_strings(i).is_struck)
          nb_corde_struck += 1;
        
        piano_strings(i).InitBuffer(DOSSIER, size_buffer, i, remove_file,
                                    output_sismo_double, output_sismo_binary);
      }
    
    size = nb_dof_before(nb_corde)*dimension;
    
    // le manche
    if (shank.is_interacting_string)
      {
        offset_shank = size;
        shank.offset_shank = offset_shank;
        shank.SetTimeStep(Deltat);
        shank.t_begin = t_begin;
        shank.t_end = t_end;
        shank.t_display = t_display;
        
        shank.DOSSIER = DOSSIER;
        shank.ConstructAll(num, output_sismo_double, output_sismo_binary);
        size += shank.GetM();        
      }
    
    // the bridge
    bridge.offset_bridge = size;
    int size_bridge = bridge.GetSize();
    size += size_bridge;

    // multiplicateurs de Lagrange pour l'interaction avec la plaque
    size_LM = nb_corde*nb_lagrange_mult;
    offset_lagrange_mult.Reallocate(nb_lagrange_mult);
    for (int p = 0; p < nb_lagrange_mult; p++)
      offset_lagrange_mult(p) = size + nb_corde*p;
    
    size += size_LM;
    
    // ddl pour la position du marteau
    if (!drop_hammer)
      {
        offset_hammer = size;
        size += 1;
        if (hammer.with_auxiliary_variable)
        {
          offset_hammer_aux = size;
          size += nb_corde;
        }
      }
    else
      offset_hammer = -1;
    
    // ddls pour les forces de chaque corde sur le marteau
    if (hammer.is_interacting)
      {
        offset_force_hammer = size;
        size += nb_corde;
      }


    if (hammer.strike_at_an_angle)
    {
    if(piano_strings(0).var_eq.type_equation == TypeEquation::WAVE_NL_2T_TIMO)  
      cout << "Hammer = YES_ANGLE used with double polarised strings. OK!" << endl;
    else
    {
      cout << "Hammer = YES_ANGLE should be used with double polarised strings only" << endl;
      abort();
    }
  }
    // allocation de tableaux
    val.Reallocate(size, 3); // initialiser val
    val.Zero();
  }

  
  //! Affichage des parametres du probleme
  /*!
    \param[in] file_out flux de sortie pour ecrire les parametres
    \param[in] write_on_file si vrai, les parametres sont ecrits sur le flux de sortie
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>
  ::Display(ostream& file_out, bool write_on_file)
  {
    hammer.Display(file_out, write_on_file);
    shank.Display(file_out, write_on_file);
    cout << "We have " << nb_corde << " strings :" << endl;
    for(int i = 0; i < piano_strings.GetM(); i++)
      piano_strings(i).Display(file_out, write_on_file);
  
  
  cout << "We have " << number_of_duplex_string << " duplex strings" << endl;
  cout << "They are string" << vect_duplex_string_ID << endl;

  
    DISP(DOSSIER);
    DISP(theta);
    cout << "Simulation    from " <<t_begin<< " sec to " 
         << t_end<< " sec with time step "<<Deltat<< " sec." <<endl;
    cout << "Interpolation from " <<t_begin_interp<< " sec to " 
         << t_end_interp<< " sec with time step "<<deltat_interp<< " sec." <<endl;
    cout << "Sismo         from " <<t_begin_sismo<< " sec to " 
         << t_end_sismo<< " sec with time step "<<deltat_sismo<< " sec." <<endl;
    cout << "Sismo Crush   from " <<t_begin_crush<< " sec to " 
         << t_end_crush<< " sec with time step "<<deltat_crush<< " sec." <<endl;
   
 
    
    if (write_on_file)
      {
  file_out << "DOSSIER = " << DOSSIER << endl;
  file_out << "theta = " << theta << endl;
  file_out << "Simulation    from " 
                 <<t_begin<< " sec to " << t_end<< " sec with time step "<<Deltat<< " sec." <<endl;
  file_out << "Interpolation from " <<t_begin_interp<< " sec to " 
                 << t_end_interp<< " sec with time step "<<deltat_interp<< " sec." <<endl;
  file_out << "Sismo         from " <<t_begin_sismo<< " sec to " 
                 << t_end_sismo<< " sec with time step "<<deltat_sismo<< " sec." <<endl;
  file_out << "Sismo Crush   from " <<t_begin_crush<< " sec to " 
                 << t_end_crush<< " sec with time step "<<deltat_crush<< " sec." <<endl;
      }
  }
  
  
  //! Les buffers de sortie sont vides  
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>::FlushBuffer()
  {
    output_hammer.CloseBuffer();
    output_force.CloseBuffer();
    output_ener.CloseBuffer();
    output_LM.CloseBuffer();
    
    for(int i = 0; i < nb_corde; i++)
      {
  piano_strings(i).output_string.CloseBuffer();
  piano_strings(i).output_last_point.CloseBuffer();
  piano_strings(i).output_crush.CloseBuffer();
      }
    
    if (shank.is_interacting_string)
      {
        shank.output_Force_d.CloseBuffer();
        shank.output_shank.CloseBuffer();
      }
  }
  
  
  //! initialisation de la solution
  /*!
    \param[in] init_val si vrai, on initialise les cordes et marteau
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>::Init(bool init_val)
  {
    if (Deltat  != Real_wp(0))
      nb_max_iter = toInteger(ceil(abs(t_end - t_begin)/Deltat ) );
    else
      {
        cout << "Enter a time step different from 0 " << endl;
        abort();
      }
    
    if (!drop_hammer)
      hammer.ComputeInterpolatedForce(*this);
    
    if (init_val)
      {
        // initial conditions
        for (int i = 0; i < nb_corde; i++)
          piano_strings(i).SetInitialConditions(type_initial_condition, amplitude_ini,
            nb_mode_ini, dimension_ini, Deltat, *this);
        
        if (!drop_hammer)
          hammer.SetInitialConditions(*this, piano_strings);
      }
    
    if ((shank.is_interacting_string) && (init_val))
      {
        shank.Init(nb_max_iter, val);
        if (shank.type_function_jack == 3) // le enum ne fonctionne pas ici mais 3 = FROM_FILE
          shank.ComputeInterpolatedJackForce(*this); // shank.FROM_FILE doit marcher je pense
        if ((shank.tangential_component) && (shank.type_function_jack_tangential == 3))
          shank.ComputeInterpolatedTangentialJackForce(*this);
        if ((shank.torque_component) && (shank.type_function_jack_torque == 3))
          shank.ComputeInterpolatedTorqueJack(*this);

        // Cette initialisation écrase celle de hammer.SetInitialcondisions
        //val(this->GetHammerDofNumber(), 0) = val(shank.xi_dof_number+1, 0) - shank.Hprime_head;
        //val(this->GetHammerDofNumber(), 1) = val(shank.xi_dof_number+1, 1) - shank.Hprime_head;
        //val(this->GetHammerDofNumber(), 2) = val(shank.xi_dof_number+1, 2) - shank.Hprime_head;
        val(this->GetHammerDofNumber(), 0) = val(shank.xi_dof_number+1, 0) - shank.distance_corde_pivot;
        val(this->GetHammerDofNumber(), 1) = val(shank.xi_dof_number+1, 1) - shank.distance_corde_pivot;
        val(this->GetHammerDofNumber(), 2) = val(shank.xi_dof_number+1, 2) - shank.distance_corde_pivot;

        // Il faut aussi remettre à jour l'initialisation des variables aux avec les bonnes valeurs de xi
        if (hammer.with_auxiliary_variable)
        {
          for(int i = 0; i < piano_strings.GetM(); i++)
          {
            // on calcule le deplacement initial moyen de la corde au niveau de l'impact du marteau
            Real_wp ux0 = 0.0;
            Real_wp ux1 = 0.0;
            for (int j = 0; j < piano_strings(i).repartition_marteau_sparse.GetM(); j++)
            {
              int jj = piano_strings(i).repartition_marteau_sparse.Index(j);
              Real_wp repj = piano_strings(i).repartition_marteau_sparse.Value(j);
              ux0 += GetValue(i, 0, jj, 0, val)*repj;
              ux1 += GetValue(i, 0, jj, 1, val)*repj; 
            }

            SetHammerAux(i, 0, 0, val);
            SetHammerAux(i, 1, hammer.GAux(ux0 - GetHammer(0, val), i), val);
            SetHammerAux(i, 2, hammer.GAux(ux1 - GetHammer(1, val), i), val);
          }
        }
      }
    
    if ((Newton_solver.type_matrix_lu == Newton_solver.LU_ARROW)
        || (Newton_solver.type_matrix_lu == Newton_solver.LU_TINY_ARROW))
      {
        cout << "Resolution with arrow solver" << endl;        
        
        Newton_solver.size_band_lu = nb_base*dimension-1;
        if (drop_hammer)
          {
            if (nb_corde > 1)
              Newton_solver.size_last_row_lu = size_LM;
            else
              Newton_solver.size_last_row_lu = 0;
          }
        else
          {
            if (nb_corde > 1)
              Newton_solver.size_last_row_lu = 1 + size_LM;        
            else
              Newton_solver.size_last_row_lu = 1;
          }
        
        if (hammer.is_interacting)
          Newton_solver.size_last_row_lu += nb_corde;
        
        // Shank additional dofs : theta, lambda, Fcoupl and xsi 
        if (shank.is_interacting_string)
          Newton_solver.size_last_row_lu += 6;
  
  // bridge stuff
  // Newton_solver.size_last_row_lu += bridge.GetSize();
      }

    ConstructLinearMatrix();

    /* val.FillRand(); Mlt(1.0/RAND_MAX, val);
    //val.Zero();
    val.WriteText("x_ref.dat");
    //val.ReadText("x_example.dat");

    Matrix<Real_wp, General, ArrayRowSparse> DiffMatrix;
    DiffMatrix.Reallocate(val.GetM(), val.GetM());
    this->ComputeDiff(DiffMatrix, val);

    DiffMatrix.WriteText("mat_df.dat");

    VectReal_wp rhs(val.GetM()), scheme(val.GetM());
    rhs.Zero(); scheme.Zero();

    //for (int i = 0; i < val.GetM(); i++)
    //SetComplexZero(val(i, 2));
    
    this->ComputeRightHandSide(rhs, val, false);
    rhs.WriteText("rhs.dat");
    
    this->ComputeScheme(scheme, rhs, val);
    
    scheme.WriteText("val.dat");
    
    int test_input; cout << "waiting" << endl; cin >> test_input; */
    
    Real_wp celerity = sqrt(piano_strings(0).T0/(piano_strings(0).rho*piano_strings(0).A));
    cout << "The celerity of the transversal wave propagation is " << celerity <<" m/s."<<endl;
    Real_wp num_celerity = piano_strings(0).Deltax/Deltat;
    cout << "The num_celerity is "<<num_celerity<<endl;
    cout << "Delta x is " << piano_strings(0).Deltax << endl;
    
    // Initialisation of the Newton Solver
    Newton_solver.Init(*this, val);
  }

  
  //! Sauvegarde des donnees necessaires a la reprise
  /*!
    \param[in] nt numero d'iteration en temps
   */
  template<int nb_base, int nb_quad, class TypeEquation> void 
  MultiString<nb_base, nb_quad, TypeEquation>::SaveDatas(int nt)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
    int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif
    
    if (save_reprise)
      if (nt % save_iter_reprise == 0)
        {
          if (rank_proc == 0)
            cout << "Saving iterates... " << endl;
          
          // prefixe des fichiers de sauvegardes, et suffixe
          string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
          string prefix = path_reprise;
          
          // on cherche le numero de dossier dans STIFFOUT
          IVect num_iter(3);
          int num_dossier = 0;
          if (rank_proc == 0)
            {
              size_t pos2 = DOSSIER.find_last_of("/");
              string chaine = DOSSIER.substr(0, pos2);
              size_t pos1 = chaine.find_last_of("/");
              num_dossier = to_num<int>(DOSSIER.substr(pos1+1, pos2-pos1-1));
              //DISP(num_dossier);
            }
          
          bool first_save = false;
          if (bascule_reprise == -1)
            {
              // premiere fois qu'on imprime
              bascule_reprise = 0;
              first_save = true;
            }
            
          // on imprime dans le fichier IterationPianoSave
          // le numero de dossier, la bascule de la reprise, le numero d'iteration
          num_iter(0) = nt;
          num_iter(1) = num_dossier;
          num_iter(2) = bascule_reprise;
          
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
          
          // on ecrit ici les nouveaux iteres (qui ne seront pas consideres pour la reprise)
          // le programme peut etre interrompu durant l'ecriture sans dommage collateral
          suffix = to_str(num_iter(2)) + "_P" + suffix;
          
          FlushBuffer();

          // on ecrit val
          val.Write(prefix+"ValPianoSave"+suffix);

          // et les donnees sur le manche/marteau
          VectReal_wp tmp(8);
          tmp(0) = beginning_time_force;
          tmp(1) = last_time_force;
          tmp(2) = final_time_shank;
          tmp(3) = hammer.energy_free_hammer;
          tmp(4) = shank.free_kinetic_energy;
          tmp(5) = shank.free_potential_energy;
          tmp(6) = shank.free_internal_energy;
          tmp(7) = shank.free_shank_energy;
          
          tmp.Write(prefix + "ForceMarteauPianoSave"+suffix);
          
          bascule_reprise = 1 - bascule_reprise;
          
          if (rank_proc == 0)
            cout << "Iterates saved" << endl;
        }
  }
  
  
  //! Iterations en temps
  template<int nb_base, int nb_quad, class TypeEquation> void 
  MultiString<nb_base,nb_quad,TypeEquation>:: RunAll()
  {    
    Init(!load_reprise);
    
    int n0 = 0;
    if (load_reprise)
      n0 = load_iter_reprise;
    
    // boucle principale
    for (int n_time = n0; n_time < nb_max_iter; n_time++)
    {
      Real_wp t = t_begin + n_time*Deltat;
      t_courant = t;
      nb_iter_courant = n_time;
      
      // on sauvegarde les instantanes si besoin est
      this->SaveDatas(n_time);
        
      Advance(t, n_time);
  
      //cout << "Theta = " << shank.GetThetaShank(val, 2) << endl;
      //cout << "Fcoupl = " << val(shank.Fcoupl_dof_number+1, 2) << endl;
      //cout << "mg = " << shank.mass_head*shank.gravity << endl;
      //cout << "Fhead = " << shank.last_Fhead_stored << endl;

      WriteEnergy();        
    }
  }
  
  
  //! Avancee du schema en temps (calcul de u^n+1, puis u^n-1 = u^n et u^n = u^{n+1})
  /*!
    \param[in] temps t^n
    \param[in] n_time numero d'iteration n
    \param[in] save_snapshot si vrai, on ecrit les snapshots, si necessaire
    \param[in] compute_time si vrai, les temps de calcul sont affiches
   */
  template<int nb_base, int nb_quad, class TypeEquation> void
  MultiString<nb_base, nb_quad, TypeEquation>
  ::Advance( Real_wp temps, int n_time, bool save_snapshot, bool compute_time)
  {
    if (print_level >= 6)
      chrono_string.Start(VirtualTimer::ALL);
  
    t_courant = temps;
    nb_iter_courant = n_time;

    if (save_snapshot)
    {
      if (print_level >= 6)
        chrono_string.Start(VirtualTimer::OUTPUT);
        
      int n = toInteger(round(temps/(t_display)));
      if ((abs(temps/(t_display) - n) <Deltat/2))
      {
        if ((compute_time) && (print_level >= 6))
        {                
          DISP(chrono_string.GetSeconds(VirtualTimer::OUTPUT));
          DISP(chrono_string.GetSeconds(VirtualTimer::SCHEME));
          DISP(chrono_string.GetSeconds(VirtualTimer::FACTO));
          DISP(chrono_string.GetSeconds(VirtualTimer::JACOBIAN));
          DISP(chrono_string.GetSeconds(VirtualTimer::SOLVE));
          DISP(chrono_string.GetSeconds(VirtualTimer::ALL));
        }
      }
        
      // we write snapshots (interp et sismo) if necessary
      for (int i = 0; i < nb_corde; i++)
        piano_strings(i).WriteStringSnapshot(n_time, temps, *this, i);
        
      WriteSnapshots(n_time, temps);
        
      if (shank.is_interacting_string)
        shank.WriteSnapshots(n_time, temps, val);

      if (print_level >= 6)
        chrono_string.Stop(VirtualTimer::OUTPUT);
    }

    // We apply the numerical scheme !
    if ((TypeEquation::linear) && (drop_hammer))
    {
      VectReal_wp& b = Newton_solver.GetScheme();
                
      // une iteration de Newton
      for (int i = 0; i < val.GetM(); i++)
        val(i, 2) = 0.0;
        
      this->ComputeRightHandSide(b, val, false);
      Newton_solver.ApplyScaling(b);
        
      Newton_solver.SolveDifferential(b, print_level);
      for (int i = 0; i < val.GetM(); i++)
        val(i, 2) = -b(i);
        
      // on verifie que ca marche
      //this->ComputeScheme(b, val);
      //Newton_solver.ApplyScaling(b); DISP(Norm2(b));

      // on passe au temps suivant
      for (int i = 0; i < val.GetM(); i++)
      {
        val(i, 0) = val(i, 1);
        val(i, 1) = val(i, 2);
      }
    }
    else
      Newton_solver.Solve(*this, val);


    // checking if the hammer should be dropped
    if ((!drop_hammer) && (drop_hammer_when_possible))
    {
      bool remove_hammer = false;

      /* // Critère de départ de marteau basé sur le temps (ne marche pas bien à cause des rebonds sur la corde)
      if (beginning_time_force < 0)
      {
        // on est dans la phase d'approche du marteau
        // qui n'a pas encore touche de cordes
        if (abs(somme_interaction_hammer) > 0)
        {
          beginning_time_force = temps;
        }
      }
      else
      {
        if (last_time_force < 0)
        {
          // on est dans la phase de contact entre le marteau et la corde
          if (somme_interaction_hammer == 0)
          {
            // le marteau vient de quitter la corde
            last_time_force = temps;
          }
        }
        else
        {
          // on est dans la phase ou le marteau a quitte la corde
          if (somme_interaction_hammer != 0)
          {
            // le marteau vient de rebondir sur la corde
            // on repart sur une phase de contact
            // beginning_time_force = temps;
            last_time_force = -1.0;
          }
          else
          {
            // on regarde si ca fait longtemps que le marteau n'a pas touche la corde
            // dans ce cas on peut enlever le marteau
            if (abs(temps-last_time_force) > 2.05*(last_time_force - beginning_time_force))
            {
              final_time_shank = temps;
              remove_hammer = true;
            }
          }
        }
      }*/

      // Critère de départ de marteau basé sur son ordonnée
      if (beginning_time_force < 0)
      {
        // on est dans la phase d'approche du marteau
        // qui n'a pas encore touche de cordes
        if (abs(somme_interaction_hammer) > 0)
        {
          beginning_time_force = temps;
        }
      }
      else
      {
        // Le marteau est déjà entré en contact au moins une fois avec la corde
        // A partir de maintenant, si il redescend sous une certaine valeur
        // on remove_hammer
        if (GetHammer(1, val) < -2*hammer.delta)
        {
          final_time_shank = temps;
          remove_hammer = true;
          cout << "Hammer has been dropped" << endl;
        }
      }

        
      if (remove_hammer)
      {
        hammer.energy_free_hammer = hammer.ComputeKineticEnergy(*this);
            
        output_hammer.CloseBuffer();
        output_force.CloseBuffer();
    
        if (shank.is_interacting_string)
        {
          shank.free_shank_energy = shank.GetEnergy(val, shank.free_kinetic_energy,
                                                         shank.free_potential_energy,
                                                         shank.free_internal_energy);
          shank.output_Force_d.CloseBuffer();
          shank.output_shank.CloseBuffer();
        }
            
        drop_hammer = true;
        hammer.is_interacting = false;
        shank.is_interacting_string = false;
      
        // we recompute the new size of the linear system
        // and update the offset
        // first the strings
        size = nb_dof_before(nb_corde)*dimension;

        // then the bridge
        int old_offset_bridge = bridge.offset_bridge;
        bridge.offset_bridge = size;
        size += bridge.GetSize();
      
        // and finally the Lagrange multipliers
        Vector<int> old_offset_lagrange_mult = offset_lagrange_mult;
      
        for (int p = 0; p < nb_lagrange_mult; p++)
          offset_lagrange_mult(p) = size + nb_corde*p;
            
        size += size_LM;
      
        //  val is updated with the new numbering
        Matrix<Real_wp> old_val(val);            
        val.Resize(size, 3);
        for (int i = 0; i < bridge.GetSize(); i++)
        {
          val(bridge.offset_bridge+i, 0) = old_val(old_offset_bridge+i, 0);
          val(bridge.offset_bridge+i, 1) = old_val(old_offset_bridge+i, 1);
          val(bridge.offset_bridge+i, 2) = old_val(old_offset_bridge+i, 2);
        }
      
        for (int i = 0; i < nb_corde; i++)
        {
          for (int p = 0; p < nb_lagrange_mult; p++)
          {
            val(offset_lagrange_mult(p)+i, 0) = old_val(old_offset_lagrange_mult(p)+i, 0);
            val(offset_lagrange_mult(p)+i, 1) = old_val(old_offset_lagrange_mult(p)+i, 1);
            val(offset_lagrange_mult(p)+i, 2) = old_val(old_offset_lagrange_mult(p)+i, 2);
          }
        }
            
        // on reinitialise Newton
        if ( (Newton_solver.type_matrix_lu == Newton_solver.LU_ARROW)
                 || (Newton_solver.type_matrix_lu == Newton_solver.LU_TINY_ARROW) )
        {
          Newton_solver.size_band_lu = nb_base*dimension;
          Newton_solver.size_last_row_lu = size_LM;
          if (nb_corde == 1)
            Newton_solver.size_last_row_lu = 0;
        }
            
        Newton_solver.Init(*this, val);
      }
    }

    //coef_couplage = 0.0;
    //for (int p = 0; p < nb_corde; p++)
    //  coef_couplage -= GetLM(p, 2, val);
    
    if (print_level >= 6)
      chrono_string.Stop(VirtualTimer::ALL);
    
    // int test_input; cout << "Waiting..." << endl; cin >> test_input;
    //return (coef_couplage);
  }
  
  
  //! writes the energy if needed
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp MultiString<nb_base, nb_quad, TypeEquation>
  ::WriteEnergy(bool comes_from_elsewhere)
  {
    bool write_ener = false;
    if (deltat_ener != 0)
    {
      int n = toInteger(round(t_courant/deltat_ener));
      if((t_courant>= t_begin_ener) && (t_courant<= t_end_ener)
          && (abs(t_courant/deltat_ener - n) <(Deltat*0.5)))
        write_ener = true;
    }
    
    bool display_ener = false;
    
    if (t_display != 0)
    {
      int n = toInteger(round(t_courant/t_display));
      if ((abs(t_courant/(t_display) - n) <Deltat/2))
        display_ener = true;
    }
    
    Real_wp EnerTot = 0.0; VectReal_wp ENERGIE;
    if (display_ener || write_ener)
      EnerTot = ComputeEnergy(ENERGIE);
    
    if (write_ener)
      output_ener.AddVect(ENERGIE);
    
    // we display unless the function is called from another part of the piano
    if (display_ener && !comes_from_elsewhere)
    {
      cout << "At time : " << t_courant << ", ";
      cout << "Total energy =  " << EnerTot << endl;
    }
    
    return EnerTot;
  }
  

  //! computes the energy of the different sub-systems and returns the total energy
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp MultiString<nb_base,nb_quad,TypeEquation>
  ::ComputeEnergy(VectReal_wp& ENERGIE)
  {
    Real_wp EnerTot = 0.0, ener_u, ener_v;
    int nb_comp = 2;
    if (separate_energy_uv)
      nb_comp = 3;

    int add_for_bridge = 0;
    if (bridge.GetSize() > 0)    
      add_for_bridge = 1;

    if (shank_present_from_t0)
      ENERGIE.Reallocate(nb_comp*nb_corde + 4 + add_for_bridge);
    else
      ENERGIE.Reallocate(nb_comp*nb_corde + 2 + add_for_bridge);
    
    ENERGIE.Fill(0);
    for (int i = 0; i < nb_corde; i ++)
    {
      // energy due to u and v
      piano_strings(i).ComputeEnergy(*this, ener_u, ener_v);
      if (separate_energy_uv)
      {
        ENERGIE(2*i) = ener_u;
        ENERGIE(2*i+1) = ener_v;
      }
      else
        ENERGIE(i) = ener_u+ener_v;
        
      EnerTot += ener_u+ener_v;
        
      // energy due to the coupling string/hammer
      ENERGIE((nb_comp-1)*nb_corde+i) = ComputeCouplingEnergy(i);
      EnerTot += ENERGIE((nb_comp-1)*nb_corde+i);
    }
    
    int nb = nb_comp*nb_corde;
    if (shank.is_interacting_string)
    {
      // shank energy
      Real_wp kinetic_energy, potential_energy, internal_energy, shank_energy;
      shank_energy = shank.GetEnergy(val, kinetic_energy, potential_energy, internal_energy);
        
      ENERGIE(nb++) = kinetic_energy;
      ENERGIE(nb++) = potential_energy;
      ENERGIE(nb++) = internal_energy;
      EnerTot += shank_energy;
    }
    else if (shank_present_from_t0)
    {
      ENERGIE(nb++) = shank.free_kinetic_energy;
      ENERGIE(nb++) = shank.free_potential_energy;
      ENERGIE(nb++) = shank.free_internal_energy;
      EnerTot += shank.free_shank_energy;        
    }
    else
    {
      // hammer energy
      ENERGIE(nb) = hammer.ComputeKineticEnergy(*this);
      EnerTot += ENERGIE(nb); nb++;
    }
    // if present, bridge energy
    if (bridge.GetSize() > 0)
    {
      Real_wp ener_rot = bridge.ComputeEnergy(*this, val);
      ENERGIE(nb) = ener_rot; nb++;
      EnerTot += ener_rot;
    }
    
    // total energy
    ENERGIE(nb) = EnerTot;
    nb++;

    // difference with initial energy
    /*if (energy_init_system > -10.0)
      ENERGIE(nb) = EnerTot - energy_init_system;
    else
      {
        ENERGIE(nb) = 0.0;
        if (shank.is_interacting_string)
          {
            if (!shank.push_jack_until_theta)
              energy_init_system = EnerTot;
          }
        else
          energy_init_system = EnerTot;
          }*/
    return EnerTot;
  }

  
  //! calcul de l'energie de couplage (entre marteau et corde) de la corde i
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp MultiString<nb_base,nb_quad,TypeEquation>
  ::ComputeCouplingEnergy(int i)
  {
    Real_wp res = 0.0;  
    if (drop_hammer)
      return res;
    
    PianoString<nb_base, nb_quad, TypeEquation>& piano_string = piano_strings(i);
    if (hammer.type_interaction == hammer.COUPLING)
    {
      if(!piano_string.is_duplex)
      {
        Real_wp ux0p  = 0.0;
        Real_wp ux0m  = 0.0;
        
        for (int j = 0; j < piano_string.repartition_marteau_sparse.GetM(); j++)
        {
          int jj = piano_string.repartition_marteau_sparse.Index(j);
          Real_wp repj = piano_string.repartition_marteau_sparse.Value(j);
          if (hammer.strike_at_an_angle)
          {
            ux0p += GetValue(i, 0, jj, 1, val)*repj*cos(hammer.strike_angle*pi_wp/180.0) + GetValue(i, 1, jj, 1, val)*repj*sin(hammer.strike_angle*pi_wp/180.0);
            ux0m += GetValue(i, 0, jj, 0, val)*repj*cos(hammer.strike_angle*pi_wp/180.0) + GetValue(i, 1, jj, 0, val)*repj*sin(hammer.strike_angle*pi_wp/180.0); 
          } 
          else    
          {
            ux0p += GetValue(i, 0, jj, 1, val)*repj;
            ux0m += GetValue(i, 0, jj, 0, val)*repj; 
          }
        }

        res = hammer.ComputePotentialEnergy(*this, i, ux0p, ux0m);
        //res = hammer.Psi(ux0p - GetHammer(1, val)) + hammer.Psi(ux0m - GetHammer(0, val));
        //res *= hammer.Ki(i)*0.5;
      }
    }
    return res;
  }
  
  
  //! Ecriture des instantanes
  /*!
    \param[in] nb_iter numero d'iteration
    \param[in] t temps t^n
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>
  ::WriteSnapshots(int nb_iter, Real_wp t)
  {
    if (deltat_sismo > 0)
      {
        //string name_file = DOSSIER+"Sismo_Hammer.txt";
        int n = toInteger(round(t/deltat_sismo));
        if((t >= t_begin_sismo) && (t <= t_end_sismo) && (abs(t/deltat_sismo -n) < (Deltat*0.5)))
          {
            // sortie du deplacement du marteau
            if (!drop_hammer)
              {
                TinyVector<Real_wp, 2> temp;
                temp(0) = t;            
                temp(1) = GetHammer(0, val);
                
                output_hammer.AddTinyVect(temp);
                
                // sortie de la Force du marteau
                //  name_file = DOSSIER + "ForceMarteau.txt";
                TinyVector<Real_wp, 2> force_sismo;   
                force_sismo(0) = t_courant;
                force_sismo(1) = somme_interaction_hammer;
                output_force.AddTinyVect(force_sismo);
              }
            
      // sortie des multiplicateurs de Lagrange
            VectReal_wp temporaire(size_LM+1);
            
            temporaire(0) = t;
            for (int i = 0; i < nb_corde; i++)
        {
                for (int p = 0; p < nb_lagrange_mult; p++)
                  temporaire(i + nb_corde*p + 1) = GetLagrangeM(i, p, 0, val);
        }
            
      output_LM.AddVect(temporaire);
          }
      }
        
    if (deltat_crush != 0)
      {
        int n = toInteger(round(t_courant/deltat_crush));
        if ((t_courant >= t_begin_crush) && (t_courant <= t_end_crush)
            && (abs(t_courant/deltat_crush - n) < (Deltat*0.5)))
          {
            for ( int i = 0 ; i < nb_corde ; i++) 
              {
                if (piano_strings(i).is_struck)
                  {                    
                    if (!drop_hammer)
                      {
                        TinyVector<Real_wp, 3> vect_crush;
                        vect_crush (0) = t;
                        vect_crush(1) = piano_strings(i).ecrasement;
                        vect_crush(2) = piano_strings(i).force_marteau;
                        
                        piano_strings(i).output_crush.AddTinyVect(vect_crush);
                      }
                  }   
              }//end for i =1 : nb_corde
          }
      }// end deltat_crush != 0
    
  }
  
  
  //! Calcul de la force du marteau a appliquer sur la corde
  /*!
    \param[in] piano_string corde pour laquelle on veut connaitre la force du marteau
    \param[in] val2 deplacement de la corde val2(:, 0) = u^{n-1},
                    val2(:, 1) = u^{n}, val2(:, 2) = u^{n+1}
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp MultiString<nb_base, nb_quad, TypeEquation>
  ::ComputeHammerInteraction(PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                             Matrix<Real_wp> & val2)
  {
    if ((!piano_string.is_struck) || (drop_hammer))
    {
      // if the string is not struck we return 0
      return 0.0;
    }
  
    // if not, we compute the interaction.
    Real_wp result  = 0.0;
    int i = piano_string.numero_string;
    if (hammer.type_interaction == hammer.COUPLING)
    {
      if(!piano_string.is_duplex)   
      {
        Real_wp ux0p  = 0.0;
        Real_wp ux0   = 0.0;
        Real_wp ux0m  = 0.0;
        
        // on calcule le deplacement moyen de la corde au niveau de l'impact du marteau
        for (int j = 0; j < piano_string.repartition_marteau_sparse.GetM(); j++)
        {
          int jj = piano_string.repartition_marteau_sparse.Index(j);
          Real_wp repj = piano_string.repartition_marteau_sparse.Value(j);
            
          if (hammer.strike_at_an_angle)
          {
            ux0p += GetValue(i, 0, jj, 2, val2)*repj*cos(hammer.strike_angle*pi_wp/180.0)
                     + GetValue(i, 1, jj, 2, val2)*repj*sin(hammer.strike_angle*pi_wp/180.0);
            ux0  += GetValue(i, 0, jj, 1, val2)*repj*cos(hammer.strike_angle*pi_wp/180.0)
                     + GetValue(i, 1, jj, 1, val2)*repj*sin(hammer.strike_angle*pi_wp/180.0);
            ux0m += GetValue(i, 0, jj, 0, val2)*repj*cos(hammer.strike_angle*pi_wp/180.0)
                     + GetValue(i, 1, jj, 0, val2)*repj*sin(hammer.strike_angle*pi_wp/180.0); 
          }
          else
          {
            ux0p += GetValue(i, 0, jj, 2, val2)*repj;
            ux0  += GetValue(i, 0, jj, 1, val2)*repj;
            ux0m += GetValue(i, 0, jj, 0, val2)*repj; 
          }
        }
        
        if (hammer.with_auxiliary_variable)
        {
          result += 0.5 * hammer.DeriveeGAux(ux0-GetHammer(1, val2), i)
                         * (GetHammerAux(i, 2, val2) + GetHammerAux(i, 1, val2));
          result += invDeuxDeltat * hammer.Phi2Aux(ux0-GetHammer(1, val2), i)
                         * ((ux0p-GetHammer(2, val2)) - (ux0m-GetHammer(0, val2)));

          // on stocke cette valeur pour plus tard la sortir dans un output_crush
          piano_string.ecrasement = ux0 - GetHammer(1, val2);
          piano_string.force_marteau = result;
        }
        else
        {
          // la force vaut Ki phi(<u> - xi) - d/dt( Ri phi(<u> - xi))
          result += hammer.Ki(i)*hammer.DiffFiniePsi(ux0p-GetHammer(2, val2), ux0m-GetHammer(0, val2));
          result += -hammer.Ri(i) * (hammer.Phi(ux0p-GetHammer(2, val2))
                         - hammer.Phi(ux0m-GetHammer(0, val2))) * invDeuxDeltat;
          
          // on stocke cette valeur pour plus tard la sortir dans un output_crush
          piano_string.ecrasement = ux0 - GetHammer(1, val2);
          piano_string.force_marteau = result;
        }
      }
      else
        result = 1e-16; //(rand()%100)*1e-18;
    }

    else if (hammer.type_interaction == hammer.EXPERIMENT)
    {
      // force directement stockee dans la variable hammer
      result = hammer.ForceMarteau(nb_iter_courant)/nb_corde_struck;
    }

    else if (hammer.type_interaction == hammer.RICKER)
    {
      // todo
      result = 0.0;
      //    tmp3=(1/(sqrt(2.0*pi_wp)*20))*exp(-(t_courant/freq-1/(60.0*5.0))
        //             *(t_courant/freq-1/(60.0*5.0))/(2.0*20*20));
      //    gaussienne
    }

    else if (hammer.type_interaction == hammer.SINUS)
    {
      // todo
      result = 0.0;
      //tmp3=delta*sin(2.0*t_courant*PI);
      // sinus
      // cout << "do we get t_courant and delta here, t_courant = " << t_courant << " and delta = " << hammer.sinparam1 << endl; 
      result = hammer.sin_amplitude*sin(2.0*t_courant*pi_wp*hammer.sin_freq);
      // abort();
    }
    else if (hammer.type_interaction == hammer.CHIRP)
    {
      k_chirp = (hammer.chirp_fend - hammer.chirp_fstart)/t_end;
      result = hammer.chirp_amplitude*sin(hammer.chirp_phase + 2*pi_wp
        *(hammer.chirp_fstart*t_courant + 0.5*k_chirp*t_courant*t_courant)); 
    }
    
    //DISP(result);
    return result;  
  }
  

  //! Calcul de la derivee de la force du marteau a appliquer sur la corde
  /*!
    \param[in] piano_string corde pour laquelle on veut connaitre la force du marteau
    \param[in] val2 deplacement de la corde val2(:, 0) = u^{n-1},
                    val2(:, 1) = u^{n}, val2(:, 2) = u^{n+1}
   */  
  template<int nb_base, int nb_quad, class TypeEquation>
  Vector<Real_wp> MultiString<nb_base, nb_quad, TypeEquation>
  ::ComputeHammerInteractionDerivative(PianoString<nb_base, nb_quad,
                                       TypeEquation> & piano_string, Matrix<Real_wp> & val2)
  {
    if ((!piano_string.is_struck)|| (drop_hammer))
    {
      Vector<Real_wp> return_zeros;
      return_zeros.Reallocate(2);
      return_zeros(0) = 0.0;
      return_zeros(1) = 0.0;
      // if the string is not struck we return 0
      return return_zeros;
    }
    
    // if not, we compute the interaction.
    Real_wp resultHammer = 0.0;
    Real_wp resultAux = 0.0;
    Vector<Real_wp> result;
    int i = piano_string.numero_string;
    
    if (hammer.type_interaction == hammer.COUPLING)
    {
      if(!piano_string.is_duplex)
      {
        Real_wp ux0p = 0.0;
        Real_wp ux0 = 0.0;
        Real_wp ux0m = 0.0;
        
        for (int j = 0; j < piano_string.repartition_marteau_sparse.GetM(); j++)
        {
          int jj = piano_string.repartition_marteau_sparse.Index(j);
          Real_wp repj = piano_string.repartition_marteau_sparse.Value(j);
          if (hammer.strike_at_an_angle)
          {
            ux0p += GetValue(i, 0, jj, 2, val2)*repj*cos(hammer.strike_angle*pi_wp/180.0)
                     + GetValue(i, 1, jj, 2, val2)*repj*sin(hammer.strike_angle*pi_wp/180.0);
            ux0  += GetValue(i, 0, jj, 1, val2)*repj*cos(hammer.strike_angle*pi_wp/180.0)
                     + GetValue(i, 1, jj, 1, val2)*repj*sin(hammer.strike_angle*pi_wp/180.0);
            ux0m += GetValue(i, 0, jj, 0, val2)*repj*cos(hammer.strike_angle*pi_wp/180.0)
                     + GetValue(i, 1, jj, 0, val2)*repj*sin(hammer.strike_angle*pi_wp/180.0); 
          }
          else      
          {
            ux0p += GetValue(i,0,jj,2,val2)*repj;
            ux0  += GetValue(i,0,jj,1,val2)*repj;
            ux0m += GetValue(i,0,jj,0,val2)*repj; 
          }
        }
        
        if (hammer.with_auxiliary_variable)
        {
          resultHammer += invDeuxDeltat * hammer.Phi2Aux(ux0-GetHammer(1, val2), i);
          resultAux += -0.5 * hammer.DeriveeGAux(ux0-GetHammer(1, val2), i);
        }
        else
        {
          // Equation II.2.13 (p172) PhD Chabassier
          resultHammer += hammer.Ki(i)*hammer.DeriveeDiffFiniePsi(ux0p-GetHammer(2,val2), ux0m-GetHammer(0,val2));
          resultHammer += -hammer.Ri(i) * hammer.PhiPrime(ux0p-GetHammer(2,val2)) * invDeuxDeltat;
        }
      }
    }
    else if (hammer.type_interaction == hammer.EXPERIMENT)
    {
      // dans le cas d'une force imposee, la derivee par rapport a u vaut 0
      resultHammer = 0.0;
    }
    else if (hammer.type_interaction == hammer.RICKER)
    {
      // todo
      resultHammer = 0.0; 
    }
    else if (hammer.type_interaction == hammer.SINUS)
    {
      // todo
      resultHammer = 2.0*pi_wp*hammer.sin_freq*hammer.sin_amplitude*cos(2.0*t_courant*pi_wp*hammer.sin_freq);
    }
    else if (hammer.type_interaction == hammer.CHIRP)
    {
      resultHammer = (2*pi_wp*(hammer.chirp_fstart + k_chirp*t_courant))*cos(hammer.chirp_phase + 2*pi_wp
                  *(hammer.chirp_fstart*t_courant + 0.5*k_chirp*t_courant*t_courant));
    }
    result.Reallocate(2);
    result(0) = resultHammer;
    result(1) = resultAux;
    return result;  
  }
  
  
  //! calcul du terme non-lineaire du systeme non-lineaire F(x)
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>
  ::ComputeNablaApproxU(TinyMatrix<Real_wp, General, dimension, nb_quad>& Qplus,
      TinyMatrix<Real_wp, General, dimension, nb_quad>& Qminus,
      TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU,
      PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
  {
    // param[out] : nablaU = sum_sigma zeta(sigma) nabla_approx_U(sigma)
    
    nablaU.Zero();
    TinyVector<int, dimension-1> sigma;
    TinyVector<Real_wp, nb_quad> Qthetaplus,Qthetaminus;
    TinyMatrix<Real_wp, General, dimension-1, nb_quad> Qsigma;
    TinyVector<Real_wp, nb_quad> partialkU;
    TinyMatrix<Real_wp, General, dimension, nb_quad> nablaUtemp;
    
    for (int k = 0; k < dimension; k++)
      {
  // Loop over the permuts
  for (int s = 0; s < nb_combinaisons; s++)
    {
      Real_wp zeta = 0.0;
      Fillsigma(sigma,zeta,s,k);
      
      Qsigma.Zero();
      FillQsigma(Qplus, Qminus, sigma, Qsigma, k);
      GetRow(Qplus, k, Qthetaplus);
      GetRow(Qminus, k, Qthetaminus);

      piano_string.var_eq.ComputeDeltakU(Qthetaplus, Qthetaminus, Qsigma,
                                               partialkU, piano_string, k);
      
      // nablaU(k,:): zeta*partialkU(:);
      nablaUtemp.Zero();
      SetRow(zeta*partialkU, k, nablaUtemp);
      nablaU += nablaUtemp;
    } //end for permutations  
      }// end k
  }
  
  
  //! calcul de la jacobienne du terme non-lineaire du systeme non-lineaire F(x)
  template<int nb_base, int nb_quad, class TypeEquation> void 
  MultiString<nb_base, nb_quad, TypeEquation>
  ::ComputeJacobApproxU(TinyMatrix<Real_wp, General, dimension, nb_quad>& Qplus,
      TinyMatrix<Real_wp, General, dimension, nb_quad>& Qminus,
      TinyArray3D<Real_wp, dimension, dimension, nb_quad>& JacobU,
      PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
  {
    JacobU.Zero();
    TinyVector<int, dimension-1> sigma;
    TinyVector<Real_wp,nb_quad> Qthetaplus,Qthetaminus;
    TinyMatrix<Real_wp, General, dimension-1, nb_quad> Qsigma;
    TinyVector<Real_wp, nb_quad> DerivativeTemp;
    TinyMatrix<Real_wp, General, dimension, nb_quad> nablaUtemp;
    
    for (int k = 0; k < dimension; k++)
      {
  // Loop over the permuts
  for(int s = 0; s < nb_combinaisons; s++)
    {
      Real_wp zeta = 0.0;
      Fillsigma(sigma,zeta,s,k);
    
      Qsigma.Zero();
      FillQsigma(Qplus, Qminus, sigma, Qsigma, k);
      GetRow(Qplus, k, Qthetaplus);
      GetRow(Qminus, k, Qthetaminus);
      // We have always a non zero contribution for k,k
      piano_string.var_eq.ComputeDerivativeDeltakU(Qthetaplus, Qthetaminus, Qsigma,
                                                         DerivativeTemp, piano_string, k, k);
      // write a unlooped version of 
      for (int g = 0; g < nb_quad; g++)
        JacobU(k, k, g) += zeta*DerivativeTemp(g);
      
      for (int indice_sigma = 0; indice_sigma < dimension-1; indice_sigma++)
        if (sigma(indice_sigma) == 2)
    {
      // il faut faire quelque chose : 
      int k2;
      if(indice_sigma < k)
        k2 = indice_sigma;
      else
        k2 = indice_sigma+1;
      
      piano_string.var_eq.ComputeDerivativeDeltakU(Qthetaplus, Qthetaminus, Qsigma,
                                                               DerivativeTemp, piano_string,
                                                               k, k2);
      // write a unlooped version of 
      for (int g = 0;  g < nb_quad; g++)
        JacobU(k, k2, g) += zeta*DerivativeTemp(g);
    }
    }
      } 
  }
  
  
  //! computes Qsigma from Qminus and Qminus (obtained by a permutation sigma)
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>
  ::FillQsigma(TinyMatrix<Real_wp, General, dimension, nb_quad>& Qplus,
         TinyMatrix<Real_wp, General, dimension, nb_quad>& Qminus,
         TinyVector<int, dimension-1>& sigma,
               TinyMatrix<Real_wp, General, dimension-1, nb_quad>& Qsigma, int k)
  {
    TinyVector<Real_wp, nb_quad> ligne_c;
    for (int c = 0; c < k; c++)
      {
  ligne_c.Zero();
  if (sigma(c) == 2)
    GetRow(Qplus, c, ligne_c);
  else
    GetRow(Qminus, c, ligne_c);
  
  SetRow(ligne_c, c, Qsigma);
      }
    
    for (int c = k+1; c < dimension; c++)
      {
  ligne_c.Zero();
  if (sigma(c-1) == 2)
    GetRow(Qplus, c, ligne_c);
  else
    GetRow(Qminus, c, ligne_c);
  
  SetRow(ligne_c, c-1, Qsigma);
      }
  }
  
  
  //! Fills the permutation array sigma 
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>
  ::Fillsigma(TinyVector<int, dimension-1>& sigma, Real_wp& zeta, int s, int k)
  {
    if (permutations == PERMUT_ALL)
      {
  // all permutations
  zeta = Vect_zeta(s); // Theta(sigma)/N!  
  for (int c = 0; c < dimension-1; c++)
    sigma(c) = Sigma(s, c); 
      }
    else if(permutations == PERMUT_TWO_WAY)
      {
  // juste les permutations de Sonia
  if (s == 0)
    {
      for (int c = 0; c < k; c++)
        sigma(c) = 0;
      
      for (int c = k; c < dimension-1; c++)
        sigma(c) = 2;
    }
  else if (s == 1)
    {
      for (int c = 0; c < k; c++)
        sigma(c) = 2;
      
      for (int c = k; c < dimension-1; c++)
        sigma(c) = 0;
    }
  
  zeta = 0.5;
      }
    else if (permutations == PERMUT_ONE)
      {
  // juste la premiere permutation
  for (int c = 0; c < k; c++)
    sigma(c) = 2;
  
  for (int c = k; c < dimension-1; c++)
    sigma(c) = 0;
  
  zeta = 1.0;
      }
  }
  
  
  //! generation des combinaisons
  template<class GenericPb>
  void GenerateCombinaisons(GenericPb & var)
  {
    int dimension = var.dimension;
    if (var.permutations == var.PERMUT_ALL)
      {
  var.nb_combinaisons=1;
  for(int t = 0; t <dimension-1; t++)
          var.nb_combinaisons *= 2;
  
  var.Sigma.Reallocate(var.nb_combinaisons, var.dimension-1);
  var.Sigma.Fill(0);
      
  // matrice des combinaisons
  GenereCombinaisons(var.Sigma);
  CalculeZeta(var.Sigma, var.Vect_zeta);
      }
    else if(var.permutations==var.PERMUT_TWO_WAY)
      var.nb_combinaisons=2;
    else if(var.permutations==var.PERMUT_ONE)
      var.nb_combinaisons=1;  
  }
  
  
  //! generation des combinaisons
  void GenereCombinaisons(Matrix<int>& Sigma)
  {
    int taille = Sigma.GetN();
    if (taille == 1)
      {
  Sigma(0,0) = 0;
  Sigma(1,0) = 2;
      }
    
    if (taille > 1)
      {
  int Ptitnb = 1;
  for(int t = 0; t < taille-1; t++)
    Ptitnb *= 2; //2^(taille-1)
      
  Matrix<int> Sigma_temp(Ptitnb, taille-1);
  GenereCombinaisons(Sigma_temp); 
  for (int i = 0; i < Ptitnb; i++)
    {
      Sigma(i, 0) = 0;
      for(int j = 1; j < taille; j++)
        Sigma(i, j) = Sigma_temp(i, j-1);
    }
  
  for (int i = 0; i < Ptitnb; i++)
    { 
      Sigma(i+Ptitnb, 0) = 2;
      for(int j = 1; j < taille; j++)
        Sigma(i+Ptitnb, j) = Sigma_temp(i, j-1);      
    }
  
  Sigma_temp.Clear();
      }
  }
  
  
  //! factorielle de n
  Real_wp fact(int n)
  {
    if (n > 1)
      return n*fact(n-1);
    else
      return 1.0;
  }
  
  
  //! Calcul des poids pour les differentes permutations
  void CalculeZeta(Matrix<int>& sigma, Vector<Real_wp> & Vect_zeta)
  {    
    int taille = sigma.GetM();
    Vect_zeta.Reallocate(taille);
    for (int n = 0; n < taille; n++)
    {
      int compt_0 = 0;
      int compt_2 = 0;
      for(int i = 0; i < sigma.GetN(); i++)
      {
        if (sigma(n,i)==0) compt_0++;
        if (sigma(n,i)==2) compt_2++;
      }
  
      Vect_zeta(n) = fact(compt_0)*fact(compt_2)/fact(sigma.GetN()+1);  
    }    
  }
  

   // Compute transformation from string axes to bridge axes via quaternion.
  /*
    Rotation sequence: \gamma, \alpha, \beta (see Tan 2017 thesis page 98 for angle)
    \gamma = angle between two transverse displacements. angleStringPolarisation
    \alpha = angle on the plane of vertical transverse and longitudinal (used in Chabaassier's thesis). angle_chevalet
    \beta = angle on the plane of horizontal transverse and longitudinal. angle_chevalet_lateral
  */
  template<int nb_base, int nb_quad, class TypeEquation> void 
  MultiString<nb_base,nb_quad,TypeEquation>
  ::ComputeTransmissionChevaletQuaternion(int i_, Real_wp angleStringPolarisation, Real_wp angle_chevalet, Real_wp angle_chevalet_lateral)
  {
  cout << "We use quaternion to compute transmission_chevalet" << endl;
  angle_chevalet_lateral = 0;
  Real_wp C1 = cos(angle_chevalet_lateral/2/180*pi_wp); //angle rotated through x-axis, vertical transverse
  Real_wp C2 = cos(angle_chevalet/2/180*pi_wp); // angle rotated through y-axis, horizontal transverse
  Real_wp C3 = cos(angleStringPolarisation/2/180*pi_wp); //angle rotated through z-axis, longitudinal

  Real_wp S1 = sin(angle_chevalet_lateral/2/180*pi_wp);
  Real_wp S2 = sin(angle_chevalet/2/180*pi_wp);
  Real_wp S3 = sin(angleStringPolarisation/2/180*pi_wp);

  Real_wp quat0 = C1*C2*C3 - S1*S2*S3;
  Real_wp quat1 = S1*C2*C3 + C1*S2*S3;
  Real_wp quat2 = C1*S2*C3 - S1*C2*S3;
  Real_wp quat3 = C1*C2*S3 + S1*S2*C3;
  
  Real_wp RM11 = 2*quat0*quat0 - 1 + 2*quat1*quat1;
  Real_wp RM12 = 2*quat1*quat2 - 2*quat0*quat3;
  Real_wp RM13 = 2*quat1*quat3 + 2*quat0*quat2;
  Real_wp RM21 = 2*quat1*quat2 + 2*quat0*quat3;;
  Real_wp RM22 = 2*quat0*quat0 - 1 + 2*quat2*quat2;
  Real_wp RM23 = 2*quat2*quat3 - 2*quat0*quat1;
  Real_wp RM31 = 2*quat1*quat3 - 2*quat0*quat2;
  Real_wp RM32 = 2*quat2*quat3 + 2*quat0*quat1;
  Real_wp RM33 = 2*quat0*quat0 - 1 + 2*quat3*quat3;
  
  
  piano_strings(i_).transmission_chevalet(0,0) = RM11;
  piano_strings(i_).transmission_chevalet(0,1) = RM12;
  piano_strings(i_).transmission_chevalet(0,2) = RM13;
  piano_strings(i_).transmission_chevalet(0,3) = 0;
  piano_strings(i_).transmission_chevalet(0,4) = 0;
  
  piano_strings(i_).transmission_chevalet(1,0) = RM21;
  piano_strings(i_).transmission_chevalet(1,1) = RM22;
  piano_strings(i_).transmission_chevalet(1,2) = RM23;
  piano_strings(i_).transmission_chevalet(1,3) = 0;
  piano_strings(i_).transmission_chevalet(1,4) = 0;
  
  piano_strings(i_).transmission_chevalet(2,0) = RM31;
  piano_strings(i_).transmission_chevalet(2,1) = RM32;
  piano_strings(i_).transmission_chevalet(2,2) = RM33;
  piano_strings(i_).transmission_chevalet(2,3) = 0;  
  piano_strings(i_).transmission_chevalet(2,4) = 0; 
   
  piano_strings(i_).transmission_chevalet(3,0) = 0;
  piano_strings(i_).transmission_chevalet(3,1) = 0;
  piano_strings(i_).transmission_chevalet(3,2) = 0;
  piano_strings(i_).transmission_chevalet(3,3) = RM21;
  piano_strings(i_).transmission_chevalet(3,4) = RM22;
  
  piano_strings(i_).transmission_chevalet(4,0) = 0;
  piano_strings(i_).transmission_chevalet(4,1) = 0;
  piano_strings(i_).transmission_chevalet(4,2) = 0;
  piano_strings(i_).transmission_chevalet(4,3) = RM11;
  piano_strings(i_).transmission_chevalet(4,4) = RM12;  
  
  if (BridgeHasHeight)
  {
    piano_strings(i_).transmission_chevalet(0,3) = -GetBridgeHeight(2);
    piano_strings(i_).transmission_chevalet(1,4) = GetBridgeHeight(2);
    piano_strings(i_).transmission_chevalet(2,3) = GetBridgeHeight(0)*RM11
                         -GetBridgeHeight(1)*RM21;  
    piano_strings(i_).transmission_chevalet(2,4) =-GetBridgeHeight(1)*RM22
                         +GetBridgeHeight(0)*RM12;        
    //~ piano_strings(i_).transmission_chevalet(3,3) = sin(angleStringPolarisation*pi_wp/180.0);
    //~ piano_strings(i_).transmission_chevalet(3,4) = cos(angleStringPolarisation*pi_wp/180.0);
    //~ piano_strings(i_).transmission_chevalet(4,3) = cos(angleStringPolarisation*pi_wp/180.0);
    //~ piano_strings(i_).transmission_chevalet(4,4) = -sin(angleStringPolarisation*pi_wp/180.0); 
  }
    
  } 
}

#define MONTJOIE_FILE_MULTI_STRING_CXX
#endif
