#ifndef MONTJOIE_FILE_PIANO_STRING_CXX

namespace Montjoie
{
  
  //! default constructor
  template<int nb_base, int nb_quad, class TypeEquation> 
  PianoString<nb_base, nb_quad, TypeEquation>
  ::PianoString()
  {
    InitDefaultValues();
  }

  
  template<int nb_base, int nb_quad, class TypeEquation> 
  void PianoString<nb_base, nb_quad, TypeEquation>
  ::InitDefaultValues()
  {
    is_struck = true;
    is_attached = false;
    numero_string = 0;
    L = 1.05;
    E = 2.02e11;
    rho = 7850;
    A = 9.7993e-7;
    I = 7.6416e-14;
    G = 7.77e10;
    k_prime = 0.85;
    amo(0) = 7.69e-3;
    T0 = 880;
    angle_chevalet = 0.0;
    angleStringPolarisation = 0.0;
    angle_chevalet_lateral = 0.0;
    is_duplex = false;
  
    rho_detuning = 1.0;

    ecrasement = 0.0;
    force_marteau = 0.0;
    
    vect_bridge_height.Reallocate(3);
    vect_bridge_height.Fill(0.0);
  }
  
  
  //! setting parameters of the string with a line of the data file
  /*!
    \param[in] description_field keyword of the line (String)
    \param[in] parameters parameters associated with the keyword
    \param[in] param_corde_data object containing parameters of all the strings
    param_corde_data is used if the string is specified by its name (A0, C2, Dd5, etc)
  */
  template<int nb_base, int nb_quad, class TypeEquation>  
  void PianoString<nb_base, nb_quad, TypeEquation>
  ::SetInputData(const string & description_field, const Vector<string> &parameters,
                 ParameterMultistring& param_corde_data)
  {
    if (!description_field.compare("String"))
    {
      Hammer hammer;
      FillPianoString(parameters, hammer, param_corde_data);
    }
  }
  

  //! setting parameters of the string
  /*!
    \param[in] parameters parameters of the string (length, Young's modulus, etc)
    \param[inout] hammer hammer striking the string 
    \param[in] param_corde_data object containing parameters of all the strings
    param_corde_data is used if the string is specified by its name (A0, C2, Dd5, etc)
  */  
  template<int nb_base, int nb_quad, class TypeEquation> 
  PianoString<nb_base, nb_quad, TypeEquation>
  :: PianoString(const Vector<string> &parameters, Hammer& hammer,
                 ParameterMultistring& param_corde_data)
    : var_eq()
  {
    InitDefaultValues();
    FillPianoString(parameters, hammer, param_corde_data);
  }
  

  //! setting parameters of the string
  /*!
    \param[in] parameters parameters of the string (length, Young's modulus, etc)
    \param[inout] hammer hammer striking the string 
    \param[in] param_corde_data object containing parameters of all the strings
    param_corde_data is used if the string is specified by its name (A0, C2, Dd5, etc)
  */    
  template<int nb_base, int nb_quad, class TypeEquation> 
  void PianoString<nb_base, nb_quad, TypeEquation>
  ::FillPianoString(const Vector<string> &parameters, Hammer& hammer,
                    ParameterMultistring& param_corde_data)
  {
    if (parameters.GetM() <= 2)
      {
  cout << "In SetInputData of PianoString" << endl;
  cout << "String needs at least three parameters, for instance :" << endl;
  cout << "String = Note desaccord N" << endl;
  cout << "Current parameters are : " << endl << parameters << endl;
  abort();
      }

    int nb = 0;
    if (isalpha(parameters(0)[0]))
      {
  // cas ou on a String = Identifiant desaccord nb_points is_struck
  // exemple : String = C2 1e-3 100 0
  
  // on a un identifiant de corde (A0, C2, Dd5)
        // on recupere alors les parametres de la note a partir
        // du plan de cordes du piano considere
        param_corde_data.GetParamNote(parameters(0), L, A, rho, T0, E, I, G, k_prime,
              amo, amoB2);
        
  // modification de la tension pour le desaccord
        T0 *= 1.0 + to_num<Real_wp>(parameters(1));
        
        nb = 2;
        if (parameters(nb) == "DAMPING")
          {
            nb++;
      if (parameters.GetM() <= 10)
        {
    cout << "In SetInputData of PianoString" << endl;
    cout << "String needs six parameters for damping, for instance :" << endl;
    cout << "String = Note desaccord DAMPING " 
                     << "amo0 amo1 amo2 amoB2_0 amoB2_1 amoB2_2" << endl;
    cout << "Current parameters are : " << endl << parameters << endl;
    abort();
        }
      
            amo(0) = to_num<Real_wp>(parameters(nb++));
            amo(1) = to_num<Real_wp>(parameters(nb++));
            amo(2) = to_num<Real_wp>(parameters(nb++));
            
            amoB2(0) = to_num<Real_wp>(parameters(nb++));
            amoB2(1) = to_num<Real_wp>(parameters(nb++));
            amoB2(2) = to_num<Real_wp>(parameters(nb++));
          }
        // DISP(amo); DISP(amoB2);
      }
    else
      {
  if (parameters.GetM() <= 11)
    {
      cout << "In SetInputData of PianoString" << endl;
      cout << "String needs 12 parameters, for instance :" << endl;
      cout << "String = L A rho T0 E I G k amo0 amo1 N is_struck" << endl;
      cout << "Current parameters are : " << endl << parameters << endl;
      abort();
    }

        // length of the string
        L       = to_num<Real_wp>(parameters(0));
        // section of the string
        A       = to_num<Real_wp>(parameters(1));
        // density (e.g. 7850 for steel)
        rho     = to_num<Real_wp>(parameters(2));
        // tension of the string
        T0      = to_num<Real_wp>(parameters(3));
        // Young's modulus
        E       = to_num<Real_wp>(parameters(4));
        // moment of inertia (pi d^4/64 for a cylinder)
        I       = to_num<Real_wp>(parameters(5));
        // shear modulus
        G       = to_num<Real_wp>(parameters(6));
        // Timoshenko coefficient
        k_prime    = to_num<Real_wp>(parameters(7));
        nb = 8;
        
        if (parameters(nb) == "FREQ_DAMPING")
          {
      if (parameters.GetM() <= 14)
        {
    cout << "In SetInputData of PianoString" << endl;
    cout << "String needs 15 parameters, for instance :" << endl;
    cout << "String = L A rho T0 E I G k FREQ_DAMPING "
                     << "amo0 amo1 amoB2_0 amoB2_1 N is_struck" << endl;
    cout << "Current parameters are : " << endl << parameters << endl;
    abort();
        }

            nb++;
            amo(0)    = to_num<Real_wp>(parameters(nb++));
            amo(1)    = to_num<Real_wp>(parameters(nb++));
            amoB2(0)    = to_num<Real_wp>(parameters(nb++));
            amoB2(1)    = to_num<Real_wp>(parameters(nb++));
          }
        else if (parameters(nb) == "PHI_DAMPING")
          {
      if (parameters.GetM() <= 16)
        {
    cout << "In SetInputData of PianoString" << endl;
    cout << "String needs 17 parameter, for instance :" << endl;
    cout << "String = L A rho T0 E I G k PHI_DAMPING "
                     << "amo0 amo1 amo2 amoB2_0 amoB2_1 amoB2_2 N is_struck" << endl;
    cout << "Current parameters are : " << endl << parameters << endl;
    abort();
        }
      
            nb++;
            amo(0)    = to_num<Real_wp>(parameters(nb++));
            amo(1)    = to_num<Real_wp>(parameters(nb++));
            amo(2)    = to_num<Real_wp>(parameters(nb++));
            amoB2(0)    = to_num<Real_wp>(parameters(nb++));
            amoB2(1)    = to_num<Real_wp>(parameters(nb++));
            amoB2(2)    = to_num<Real_wp>(parameters(nb++));
          }
        else
          {
            amo(0)    = to_num<Real_wp>(parameters(nb++));
            amo(1)    = to_num<Real_wp>(parameters(nb++));
          }
      }
    
          
    if ( to_num<Real_wp>(parameters(nb)) <= 2 )  
      rho_detuning = to_num<Real_wp>(parameters(nb++)); //density detuning
    else
      {
        cout << "Rho_detuning obtained is: " << rho_detuning << ". It should be less than 2. Is there a mistake?" << endl;
        cout << "parameters are " << parameters << endl;
        abort();
      }
    
    this->Nx  = to_num<int>(parameters(nb++));
    is_struck = to_num<bool>(parameters(nb++));
    
    // on peut être frappee par une force exterieure
    // i.e. is_struck = 1 mais pas de Ki.
    if (parameters.GetM() > nb)
    {
      // parametres additionnels : raideur et amortissement du marteau
      Real_wp K_temp = to_num<Real_wp>(parameters(nb++));
      hammer.Ki.PushBack(K_temp);
      Real_wp R_temp = 0;
      if (parameters.GetM() > nb)
        R_temp = to_num<Real_wp>(parameters(nb++));
  
      hammer.Ri.PushBack(R_temp);
    }      
  }
  
  
  //! affichage des parametres de la corde
  /*!
    \param[inout] file_out flux de sortie ou on ecrit les parametres de la corde
    \param[in] write_on_file si vrai, les parametres sont ecrits dans le flux de sortie
  */
  template<int nb_base, int nb_quad, class TypeEquation> 
  void PianoString<nb_base, nb_quad, TypeEquation>
  ::Display(ostream& file_out, bool write_on_file)
  {
    cout << " ------- String "<< numero_string << " Parameters: ------- "<<endl;
    DISP(epsilon_machine);
    DISP(L);   
    DISP(A);
    DISP(rho);
    DISP(T0); 
    DISP(E);     
    DISP(I);     
    DISP(G);     
    DISP(k_prime);
    DISP(amo);
    DISP(amoB2);
    DISP(is_struck);
    DISP(is_attached);
    DISP(is_duplex);
    DISP(rho_detuning);
    DISP(this->Nx); 
    DISP(this->point_observation);
    DISP(angle_chevalet);
    DISP(angleStringPolarisation);
    DISP(vect_bridge_height);
    //~ DISP(NU);
    //~ DISP(NUORTH);
    cout << "long TC " << setprecision(30) << transmission_chevalet << endl;
    cout << "long stuff done" << setprecision(16) << endl;
    DISP(transmission_chevalet);
    DISP(var_eq.type_equation);
    Real_wp f0 = sqrt(T0/(rho*A))/(2*L);
    DISP(f0);
    cout << " ------------------------------------ "<<endl;
    
    if (write_on_file)
      {
        file_out << " ------- String "<< numero_string << " Parameters: ------- "<<endl;
        file_out << "L = " << L << endl;
        file_out << "A = " << A << endl;
        file_out << "rho = " << rho << endl;
        file_out << "T0 = " << T0 << endl; 
        file_out << "E = " << E << endl;
        file_out << "I = " << I << endl;
        file_out << "G = " << G << endl;
        file_out << "Kappa = " << k_prime << endl;
        file_out << "amo = " << amo << endl;
        file_out << "amoB2 = " << amoB2 << endl;
        if (is_struck)
          file_out << "The string is struck by the hammer " << endl;
        else
          file_out << "The string is not struck by the hammer " << endl;
        
        if (is_attached)
          file_out << "The string is attached " << endl;
        else
          file_out << "The string is not attached " << endl;
        
        file_out << "Nx = " << this->Nx << endl;
        file_out << "point_observation = " << this->point_observation << endl;
        file_out << "angle_chevalet = " << angle_chevalet << endl;
        //~ file_out << "NU = " << NU << endl;
        //~ file_out << "NUORTH = " << NUORTH << endl;
        Real_wp f0 = sqrt(T0/(rho*A))/(2*L);
        file_out << "f0 = " << f0 << endl;
        cout << " ------------------------------------ "<<endl;        
      }
  }
  
  
  //! sets initial condition on the current string
  /*!
    \param[in] type_initial_condition type of initial condition (NONE, SINUS or GAUSSIAN)
    \param[in] amplitude amplitude of the initial condition
    \param[in] nb_mode mode number (for a sinus, the initial condition is sin(k pi x/L)
    where k is the mode number
    \param[in] dimension_ini component of the solution to set
    \param[in] Deltat time step (not used actually)
    \param[in] var object multistring
  */
  template<int nb_base, int nb_quad, class TypeEquation> template<class GenericPb>
  void PianoString<nb_base, nb_quad, TypeEquation>
  ::SetInitialConditions(int type_initial_condition, Real_wp amplitude, int nb_mode,
       int dimension_ini, Real_wp Deltat, GenericPb & var )
  {
    switch (type_initial_condition)
      {
      case NONE:
        //on ne fait rien 
        break;
      case SINUS:
        {
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
      for (int p = 0; p < nb_base; p++)
        {
    int j = this->mesh.GetNumberDof(ne,p);
    Real_wp x = this->Deltax* ( ne + this->base.Points(p));// le point associe au ddl j
    Real_wp v = amplitude*sin(nb_mode*pi_wp*x/L); // valeur de la condition initiale
    
    var.SetValue(numero_string, dimension_ini, j, 0, v, var.val);
    var.SetValue(numero_string, dimension_ini, j, 1, v, var.val);
    var.SetValue(numero_string, dimension_ini, j, 2, v, var.val);
        }
  }

        break;
      case DOUBLE_SINUS:
        {
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
      for (int p = 0; p < nb_base; p++)
        {
    int j = this->mesh.GetNumberDof(ne,p);
    Real_wp x = this->Deltax* ( ne + this->base.Points(p));// le point associe au ddl j
    Real_wp v = amplitude*(sin(nb_mode*pi_wp*x/L)+sin((nb_mode+2)*pi_wp*x/L)); // valeur de la condition initiale
    
    var.SetValue(numero_string, dimension_ini, j, 0, v, var.val);
    var.SetValue(numero_string, dimension_ini, j, 1, v, var.val);
    var.SetValue(numero_string, dimension_ini, j, 2, v, var.val);
        }
  }

        break;

      case GAUSSIAN:
        {
    //A ecrire
    cout << "Gaussian ini cond not implemented yet.";
    abort();
  }
      case SINUSMULTI:
  {
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
      for (int p = 0; p < nb_base; p++)
        {
    int j = this->mesh.GetNumberDof(ne,p);
    Real_wp x = this->Deltax* ( ne + this->base.Points(p));// le point associe au ddl j
    
    for (int i = 0; i < dimension; i++)
      {
        Real_wp index_ini = var.vect_dimension_ini(i);
        Real_wp v = var.vect_amplitude_ini(i)*sin(var.vect_nb_mode_ini(i)*pi_wp*x/L); 
        
        var.SetValue(numero_string, index_ini, j, 0, v, var.val);
        var.SetValue(numero_string, index_ini, j, 1, v, var.val);
        var.SetValue(numero_string, index_ini, j, 2, v, var.val); 
      
      }

        }
  }
  break;
  case NONZERO:
  {
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
      {
      for (int p = 0; p < nb_base; p++)
      {
        int j = this->mesh.GetNumberDof(ne,p);
          for (int i = 0; i < dimension; i++)
            {
              var.SetValue(numero_string, i, j, 0, 1e-16, var.val);
              var.SetValue(numero_string, i, j, 1, 1e-16, var.val);
              var.SetValue(numero_string, i, j, 2, 1e-16, var.val); 
            }
      }
    } 
  }
  
        break;
        
        
  case PLUCK:
        {
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
      for (int p = 0; p < nb_base; p++)
        {
    int j = this->mesh.GetNumberDof(ne,p);
    Real_wp x = this->Deltax* ( ne + this->base.Points(p));// le point associe au ddl j
    Real_wp v = 0;
    
    if (x <= var.pluck_position)
      v = amplitude*x/var.pluck_position;
    else if ((x > var.pluck_position) && (x < 1))
      v = amplitude*(L-x)/(L-var.pluck_position);
    
    var.SetValue(numero_string, dimension_ini, j, 0, v, var.val);
    var.SetValue(numero_string, dimension_ini, j, 1, v, var.val);
    var.SetValue(numero_string, dimension_ini, j, 2, v, var.val);
        }
  }
  
        break;    
      case PLUCKANGLE:
        {
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
      for (int p = 0; p < nb_base; p++)
        {
    int j = this->mesh.GetNumberDof(ne,p);
    Real_wp x = this->Deltax* ( ne + this->base.Points(p));// le point associe au ddl j
    Real_wp v = 0;
    Real_wp v_cos = 0;
    Real_wp v_sin = 0;
    
    if (x <= var.pluck_position)
                  v = amplitude*x/var.pluck_position;
    else if ((x > var.pluck_position) && (x < 1))
                  v = amplitude*(L-x)/(L-var.pluck_position);
    
                if (abs(cos(var.pluckangle_ini*pi_wp/180)) < 1e-16)
                  v_cos = v*1e-16;
                else
                  v_cos = v*cos(var.pluckangle_ini*pi_wp/180);
    
                if (abs(sin(var.pluckangle_ini*pi_wp/180)) < 1e-16)
                  v_sin = v*1e-16;
                else
                  v_sin = v*sin(var.pluckangle_ini*pi_wp/180);
    
    
    var.SetValue(numero_string, dimension_ini, j, 0, v_cos, var.val);
    var.SetValue(numero_string, dimension_ini, j, 1, v_cos, var.val);
    var.SetValue(numero_string, dimension_ini, j, 2, v_cos, var.val);
    
    var.SetValue(numero_string, dimension_ini+1, j, 0, v_sin, var.val);
    var.SetValue(numero_string, dimension_ini+1, j, 1, v_sin, var.val);
    var.SetValue(numero_string, dimension_ini+1, j, 2, v_sin, var.val);
        }
  }
  
        break; 
      }     
  }
  
  
  //! Ecriture des instantanes de la corde (sismogramme et interpolee)
  /*!
    \param[in] nb snapshot number
    \param[in] t time
    \param[in] var object multistring
    \param[in] i numero de la solution a ecrire (0, 1 ou 2)
  */
  template<int nb_base, int nb_quad, class TypeEquation> template<class GenericPb>
  void PianoString<nb_base, nb_quad, TypeEquation>
  ::WriteStringSnapshot(int nb, Real_wp t, GenericPb & var, int i)
  {
    // ecriture de l'interpolation
    if (var.deltat_interp != 0)
      {
        int n = toInteger(round(t/var.deltat_interp));
        if ((t >= var.t_begin_interp) && (t <= var.t_end_interp)
            && (abs(t/var.deltat_interp - n) <(var.Deltat*0.5)))
    this->WriteInterp(n,var, i);
  
      }   
    
    // ecriture des sismos
    if(var.deltat_sismo!=0)
      {
  int n = toInteger(round(t/var.deltat_sismo));
  if ((t >= var.t_begin_sismo) && (t <= var.t_end_sismo)
            && (abs(t/var.deltat_sismo - n)<(var.Deltat*0.5)))
          this->WriteSismo(t,var, i);
      }
  }
  

  //! calcul du systeme non-lineaire a resoudre F(X)
  /*!
    \param[out] scheme vecteur F(X)
    \param[in] val variable X pour laquelle on veut evaluer F
    \param[in] var objet multistring
  */
  template<int nb_base, int nb_quad, class TypeEquation>  template<class GenericPb> 
  Real_wp PianoString<nb_base, nb_quad, TypeEquation> 
  ::ComputeScheme(VectReal_wp & scheme, Matrix<Real_wp> & val, GenericPb & var)
  {
    
    //Real_wp theta = var.theta(0);
    //Real_wp thetastiff = var.theta(1);
    //Real_wp invDeltat2 = var.invDeltat2;
    //Real_wp invDeuxDeltat = var.invDeuxDeltat;
    
    int i = numero_string;
    Real_wp interaction_hammer = 0.0;
    
    // ****************************
    // Treat Major Scheme Part (FV)
    // ****************************
    
    //TinyMatrix<Real_wp, General, dimension, nb_quad> TMplus, TM, TMminus;
    //TinyMatrix<Real_wp, General, dimension, nb_quad> dtt, delta, dx_delta,
    // thetastiffvect, dthetavect, dthetastiffvect;
    
    //TinyMatrix<Real_wp, General, dimension, nb_quad> Mdtt, Rdelta,
    //  Adthetavect, AdthetaStiffvect, Bthetavect,
    //  Ddthetavect, Cthetavect, Adamping_deltat;
    //TinyMatrix<Real_wp, General, dimension, nb_quad> vect_contre_phi, vect_contre_dphi;

    TinyVector<Real_wp, nb_quad> quad_values;

    TinyMatrix<Real_wp, General, dimension, nb_quad> dTMplus, dTM, dTMminus;
    TinyMatrix<Real_wp, General, dimension, nb_quad> nablaU;
    TinyMatrix<Real_wp, General, dimension, nb_base> res_int_phi, res_int_dphi;

    TinyVector<TinyVector<Real_wp, nb_base>, dimension> Unp1, res_Unp1, Un;
    Real_wp coefNL = (E*A - T0);
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
    {
      for (int j = 0; j < nb_base; j++)
      {
        int jglob = this->mesh.GetNumberDof(ne, j);
        for (int k = 0; k < dimension; k++)
        {
          int num_row = var.GetStringDofNumber(i, k, jglob);
          Unp1(k)(j) = val(num_row, 2);
          Un(k)(j) = 0.25*Unp1(k)(j) + 0.5*val(num_row, 1) + 0.25*val(num_row, 0);
        }
      }
        
      // linear part
      TypeEquation::ApplyMatrixEquation(mat_Unp1, Unp1, res_Unp1);
  
      // non-linear part
      if (!TypeEquation::linear)
      {
        for (int k = 0; k < dimension; k++)
        {
          MltTrans(this->valdPhi, Un(k), quad_values);
          quad_values *= this->invDeltax;
                
          SetRow(quad_values, k, dTM);
        }
            
        if (var.use_naive_non_linear_term)
          TypeEquation::ApplyNonLinearTerm(coefNL, dTM, nablaU);
        else
        {
          this->ComputedTM(i, ne, 2, val, var, dTMplus);
          this->ComputedTM(i, ne, 1, val, var, dTM);
          this->ComputedTM(i, ne, 0, val, var, dTMminus);
                
          if (var.solver == var.THETA_NL_HALF)
          {
            var.ComputeNablaApproxU (dTMplus , dTMminus , nablaU, *this);
          }   
          else if (var.solver ==var.THETA_NL_QUARTER)
          {
            TinyMatrix<Real_wp, General, dimension, nb_quad> dTMplushalf
                      = Real_wp(0.5)*(dTMplus+dTM);
                    
            TinyMatrix<Real_wp, General, dimension, nb_quad> dTMminushalf
                      = Real_wp(0.5)*(dTM+dTMminus);
                    
            var.ComputeNablaApproxU(dTMplushalf ,dTMminushalf , nablaU, *this);
          }           
        }
            
        this->CalculeIntGrad(nablaU, res_int_dphi);
            
        for (int j = 0; j < nb_base; j++)
        {
          int jglob = this->mesh.GetNumberDof(ne, j);
          for (int k = 0; k < dimension; k++)
            var.AddSchemeString(i, k, jglob, res_Unp1(k)(j) + res_int_dphi(k, j), scheme);
        }
      }
      else
      {
        for (int j = 0; j < nb_base; j++)
        {
          int jglob = this->mesh.GetNumberDof(ne, j);
          for (int k = 0; k < dimension; k++)
            var.AddSchemeString(i, k, jglob, res_Unp1(k)(j), scheme);
        }
      }
        
        /*
  // We fill the TM with values at quad points.
  ComputeTM (i, ne, 2, val, var, TMplus);
  ComputeTM (i, ne, 1, val, var, TM);
  ComputeTM (i, ne, 0, val, var, TMminus);
  
  ComputedTM(i, ne, 2, val, var, dTMplus);
  ComputedTM(i, ne, 1, val, var, dTM);
  ComputedTM(i, ne, 0, val, var, dTMminus);
      
  // We prep the tinymatrices
  dtt = (TMplus - Real_wp(2.0)*TM + TMminus)*invDeltat2;
  delta = (TMplus - TMminus) * invDeuxDeltat;
  dx_delta = (dTMplus - dTMminus) * invDeuxDeltat;
  thetastiffvect = thetastiff*TMplus  + (1.0-2.0*thetastiff)*TM  + thetastiff*TMminus;
  dthetavect = theta*dTMplus + (1.0-2.0*theta)*dTM + theta*dTMminus;
  dthetastiffvect = thetastiff*dTMplus + (1.0-2.0*thetastiff)*dTM + thetastiff*dTMminus;
  
  // We apply the linear terms
  var_eq.ApplyM(dtt, Mdtt, *this);
  var_eq.ApplyR(delta, Rdelta, *this);
  var_eq.ApplyAGrad(dthetavect, Adthetavect, *this);
  var_eq.ApplyAStiffGrad(dthetastiffvect, AdthetaStiffvect, *this);
  var_eq.ApplyB(thetastiffvect, Bthetavect, *this);
  var_eq.ApplyC(thetastiffvect, Cthetavect, *this);
  var_eq.ApplyD(dthetastiffvect, Ddthetavect, *this);
  var_eq.ApplyAdamping(dx_delta, Adamping_deltat, *this);
  
  // We compute the NL part
  if (TypeEquation::linear)
        {
        // linear equation, nablaU is equal to 0
        }
        else if   (var.solver == var.THETA_NL_HALF)
        {
        var.ComputeNablaApproxU (dTMplus , dTMminus , nablaU, *this);
        }   
  else if (var.solver ==var.THETA_NL_QUARTER)
        {
        TinyMatrix<Real_wp, General, dimension, nb_quad> dTMplushalf
        = Real_wp(0.5)*(dTMplus+dTM);
        TinyMatrix<Real_wp, General, dimension, nb_quad> dTMminushalf
        = Real_wp(0.5)*(dTM+dTMminus);
        var.ComputeNablaApproxU(dTMplushalf ,dTMminushalf , nablaU, *this);
        }
        
  // We prep vectors for numerical integration
  vect_contre_phi = Mdtt + Rdelta +  Cthetavect+ Ddthetavect;
  vect_contre_dphi = Adthetavect + AdthetaStiffvect  + Bthetavect + nablaU + Adamping_deltat;
  
  res_int_phi.Zero(); res_int_dphi.Zero();
  CalculeIntPhi(vect_contre_phi, res_int_phi );
  CalculeIntGrad(vect_contre_dphi, res_int_dphi);
    
  
  // We copy that in the scheme.
  for (int j = 0; j < nb_base; j++)
        {
        int jglob = this->mesh.GetNumberDof(ne, j);
        for (int k = 0; k < dimension; k++)
        {
        Real_wp value = res_int_phi(k, j) + res_int_dphi(k, j);
        var.AddSchemeString(i, k, jglob, value, scheme);
        }
        }
        */
    }
    
    
    // ****************************
    // Treat Hammer Part
    // ****************************
    
    // on recupere la force du marteau sur la corde
    interaction_hammer = var.ComputeHammerInteraction(*this, val);
    
    // on ecrit l'equation force = interaction_marteau
    // dans le cas ou le mouvement du marteau est regi par une equation
    Real_wp force_hammer(0);
    if (var.hammer.is_interacting)
    {
      // equation force_marteau - interaction_marteau = 0
      force_hammer = val(var.offset_force_hammer + i, 2);
      scheme(var.offset_force_hammer + i) = force_hammer - interaction_hammer;
    }
    else
      force_hammer = interaction_hammer;
    
    // si le marteau touche la corde
    if ((is_struck) && (!var.drop_hammer))
    {
      // on applique Fi sur la corde
      for (int j = 0; j < this->repartition_marteau_sparse.GetM(); j++)
      {
        int jj = this->repartition_marteau_sparse.Index(j);
        Real_wp repj = this->repartition_marteau_sparse.Value(j);
      
        if (var.hammer.strike_at_an_angle)
        {
          var.AddSchemeString(i, 0, jj, cos(var.hammer.strike_angle*pi_wp/180.0)*force_hammer*repj, scheme);
          var.AddSchemeString(i, 1, jj, sin(var.hammer.strike_angle*pi_wp/180.0)*force_hammer*repj, scheme);
        }
        else
          var.AddSchemeString(i, 0, jj, force_hammer*repj, scheme);
      }
    }
    
    // ****************************
    // Treat Boundary Conditions 
    // ****************************
    
    for (int k = 0; k < dimension; k++)
    {
      if (var_eq.Dirichlet_Agraffe(k))
      {
        // We set dirichlet condition on the first point
        var.SetSchemeString(i, k, 0, var.GetValue(i, k, 0, 2, val), scheme);
      }
      if (var_eq.Dirichlet_Bridge(k))
      {
        // the last point
        if (is_attached)
        {
          // We set Dirichlet condition where needed on the last point if is_attached
          var.SetSchemeString(i, k, this->mesh.GetNbDof()-1,
                                  var.GetValue(i, k, this->mesh.GetNbDof()-1, 2, val), scheme);
        }
        else
        {
          // term - lambda_k^n phi_i(L) in the string equation
          // since lambda_k is an unknown, it is taken at index 2 of val
          for (int p = 0; p < var.nb_lagrange_mult; p++)
          {
            var.AddSchemeString(i, k, this->mesh.GetNbDof()-1,
              -var.GetLagrangeM(i, p, 2, val)*transmission_chevalet(p, k), scheme);
          }
        }
      }
    }
    return force_hammer;
  }
  
  
  //! computation of right hand side for linear scheme
  template<int nb_base, int nb_quad, class TypeEquation>
  template<class GenericPb> 
  void PianoString<nb_base, nb_quad, TypeEquation> 
  ::ComputeRightHandSide(VectReal_wp & scheme, Matrix<Real_wp> & val,
                         GenericPb & var, bool newton_rhs)
  {
    
    int i = numero_string;
    
    // ****************************
    // Treat Major Scheme Part (FV)
    // ****************************
    
    /*Real_wp theta = var.theta(0);
      Real_wp thetastiff = var.theta(1);
      //Real_wp invDeltat = var.invDeltat;
      Real_wp invDeltat2 = var.invDeltat2;
      Real_wp invDeuxDeltat = var.invDeuxDeltat;
    
    
      TinyMatrix<Real_wp, General, dimension, nb_quad> TM, TMminus, dTM, dTMminus;
      TinyMatrix<Real_wp, General, dimension, nb_quad> dtt, delta, dx_delta,
      thetastiffvect, dthetavect, dthetastiffvect;
    
      TinyMatrix<Real_wp, General, dimension, nb_quad> Mdtt, Rdelta, Adthetavect,
      AdthetaStiffvect, Bthetavect,
      Ddthetavect, Cthetavect, Adamping_deltat;
    
      TinyMatrix<Real_wp, General, dimension, nb_quad> vect_contre_phi, vect_contre_dphi;
      TinyMatrix<Real_wp, General, dimension, nb_base> res_int_phi, res_int_dphi; */
    
    TinyVector<TinyVector<Real_wp, nb_base>, dimension> Un, Unm1, res_Un, res_Unm1;
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
    { 
      for (int j = 0; j < nb_base; j++)
      {
        int jglob = this->mesh.GetNumberDof(ne, j);
        for (int k = 0; k < dimension; k++)
        {
          int num_row = var.GetStringDofNumber(i, k, jglob);
          Unm1(k)(j) = val(num_row, 0);
          Un(k)(j) = val(num_row, 1);
        }
      }
        
      TypeEquation::ApplyMatrixEquation(mat_Unm1, Unm1, res_Unm1);

      TypeEquation::ApplyMatrixEquation(mat_Un, Un, res_Un);
        
      for (int j = 0; j < nb_base; j++)
      {
        int jglob = this->mesh.GetNumberDof(ne, j);
        for (int k = 0; k < dimension; k++)
        {
          //Real_wp value = res_Unm1(k*nb_base+j) + res_Un(k*nb_base+j);
          Real_wp value = res_Unm1(k)(j) + res_Un(k)(j);
          var.AddSchemeString(i, k, jglob, value, scheme);
        }
      }

        /*
  // We fill the TM with values at quad points.
        ComputeTM (i, ne, 1, val, var, TM);
  ComputeTM (i, ne, 0, val, var, TMminus);
  
        ComputedTM(i, ne, 1, val, var, dTM);
  ComputedTM(i, ne, 0, val, var, dTMminus);
      
  // We prep the tinymatrices
  dtt = (- Real_wp(2.0)*TM + TMminus)*invDeltat2;
  delta = (- TMminus) * invDeuxDeltat;
  dx_delta = - dTMminus * invDeuxDeltat;
  thetastiffvect = (1.0-2.0*thetastiff)*TM  + thetastiff*TMminus;
  dthetavect = (1.0-2.0*theta)*dTM + theta*dTMminus;
  dthetastiffvect = (1.0-2.0*thetastiff)*dTM + thetastiff*dTMminus;
  
  // We apply the linear terms
  var_eq.ApplyM(dtt, Mdtt, *this);
  var_eq.ApplyR(delta, Rdelta, *this);
  var_eq.ApplyAGrad(dthetavect, Adthetavect, *this);
  var_eq.ApplyAStiffGrad(dthetastiffvect, AdthetaStiffvect, *this);
  var_eq.ApplyB(thetastiffvect, Bthetavect, *this);
  var_eq.ApplyC(thetastiffvect, Cthetavect, *this);
  var_eq.ApplyD(dthetastiffvect, Ddthetavect, *this);
  var_eq.ApplyAdamping(dx_delta, Adamping_deltat, *this);
  
  // We prep vectors for numerical integration
  vect_contre_phi = Mdtt + Rdelta +  Cthetavect+ Ddthetavect;
  vect_contre_dphi = Adthetavect + AdthetaStiffvect  + Bthetavect + Adamping_deltat;
  
  res_int_phi.Zero(); res_int_dphi.Zero();
  CalculeIntPhi(vect_contre_phi, res_int_phi );
  CalculeIntGrad(vect_contre_dphi, res_int_dphi);
    
  
  // We copy that in the scheme.
  for (int j = 0; j < nb_base; j++)
        {
        int jglob = this->mesh.GetNumberDof(ne, j);
        for (int k = 0; k < dimension; k++)
        {
        Real_wp value = res_int_phi(k, j) + res_int_dphi(k, j);
        var.AddSchemeString(i, k, jglob, value, scheme);
        }
        }
        */
    }
        
    // Dirichlet condition
    for (int k = 0; k < dimension; k++)
    {
      if (var_eq.Dirichlet_Agraffe(k))
      {
        // We set dirichlet condition on the first point
        var.SetSchemeString(i, k, 0, 0.0, scheme);
      }
      if (var_eq.Dirichlet_Bridge(k))
      {
        if (is_attached)
        {
          // We set Dirichlet condition where needed on the last point if is_attached
          var.SetSchemeString(i, k, this->mesh.GetNbDof()-1, 0.0, scheme);
        }
      }
    }
  }

  
  //! computes u \cdot nu for u^{n-1} and u^{n+1}
  /*!
    \param[in] val solution u at times t^{n-1}, t^n and t^{n+1}
    \param[in] var multistring object
    \param[out] somme0  u^{n-1} \cdot nu
    \param[out] somme2  u^{n+1} \cdot nu
  */
  template<int nb_base, int nb_quad, class TypeEquation>
  template<class GenericPb> 
  void  PianoString<nb_base,nb_quad,TypeEquation>
  ::GetScalarProductAtBridge(Matrix<Real_wp> & val, GenericPb & var,
                             Real_wp& somme0, Real_wp& somme2)
  {
    int i = numero_string;
    somme2 = 0.0;
    somme0 = 0.0;
    
    for (int k = 0; k < dimension; k++)
    {
      //~ somme0 += var.GetValue(i, k, this->mesh.GetNbDof()-1, 0, val)*NU(k);
      //~ somme2 += var.GetValue(i, k, this->mesh.GetNbDof()-1, 2, val)*NU(k);   
      somme0 += var.GetValue(i, k, this->mesh.GetNbDof()-1, 0, val)*transmission_chevalet(0,k);
      somme2 += var.GetValue(i, k, this->mesh.GetNbDof()-1, 2, val)*transmission_chevalet(0,k); 
    } 
  }
  
  
  //! Calcul de la jacobienne DF(X)
  /*!
    \param[inout] DiffMatrix differential matrix
    \param[in] val solution X for which the differential matrix is computed
    \param[in] num_xsi dof number for the position of the shank (if a shank is present)
    \param[in] var multistring object
  */
  template<int nb_base, int nb_quad, class TypeEquation>
  template<class GenericPb> 
  Real_wp PianoString<nb_base, nb_quad, TypeEquation>
  ::ComputeDiff(VirtualMatrix<Real_wp>& DiffMatrix, Matrix<Real_wp> & val, int num_xsi, GenericPb & var)
  {
  
    int i = numero_string;
    Real_wp theta = var.theta(0);
    Real_wp thetastiff = var.theta(1);
    Real_wp invDeltat2 = var.invDeltat2;
    Real_wp invDeuxDeltat = var.invDeuxDeltat;
    
    // ****************************
    // Treat Major Scheme Part (FV)
    // ****************************
    
    TinyMatrix<Real_wp, General, dimension, nb_quad> TMBase, dTMBase, dTMplus, dTM, dTMminus;
    TinyMatrix<Real_wp, General, dimension, nb_quad> Mphi, Rphi, Adphi, AStiffdphi,
                                                    Adamping_dphi, Bphi, Cphi, Ddphi, NLpart;
    
    TinyMatrix<Real_wp, General, dimension, nb_quad> vect_contre_phi, vect_contre_dphi;
    TinyMatrix<Real_wp, General, dimension, nb_base> res_int_phi, res_int_dphi;
    TinyArray3D<Real_wp, dimension, dimension, nb_quad>  LocalJacobNL;
    
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
    {
      this->ComputedTM(i, ne, 2, val, var, dTMplus);
      this->ComputedTM(i, ne, 1, val, var, dTM);
      this->ComputedTM(i, ne, 0, val, var, dTMminus);
  
      // We compute the NL part one time per element
      Real_wp coef_NL = 0.0;
      if   (var.solver == var.THETA_NL_HALF)
      {
        var.ComputeJacobApproxU (dTMplus, dTMminus, LocalJacobNL, *this);
        coef_NL = 1.0;
      }   
      else if (var.solver == var.THETA_NL_QUARTER)
      {
        TinyMatrix<Real_wp, General, dimension, nb_quad> dTMplushalf
              = Real_wp(0.5)*(dTMplus+dTM);
        TinyMatrix<Real_wp, General, dimension, nb_quad> dTMminushalf
              = Real_wp(0.5)*(dTM+dTMminus);
        var.ComputeJacobApproxU (dTMplushalf, dTMminushalf, LocalJacobNL, *this);
        coef_NL = 0.5;
      }
  
      for (int k2 = 0; k2 < dimension; k2++)
      {
        for (int j2 = 0; j2 < nb_base; j2++)
        {
          // it will be the column indices
    
          // We need an elementary matrix
          TMBase.Zero(); dTMBase.Zero();
          for (int g = 0; g < nb_quad; g++)
          {
            TMBase(k2, g) = this->valPhi(j2, g);
            dTMBase(k2, g) = this->valdPhi(j2, g)*this->invDeltax;
          }
    
          // We apply every linear part
          var_eq.ApplyM (TMBase, Mphi, *this);
          var_eq.ApplyR (TMBase, Rphi, *this);
          var_eq.ApplyAGrad (dTMBase, Adphi, *this);
          var_eq.ApplyAStiffGrad (dTMBase, AStiffdphi, *this);
          var_eq.ApplyB (TMBase, Bphi, *this);
          var_eq.ApplyC (TMBase, Cphi, *this);
          var_eq.ApplyD (dTMBase, Ddphi, *this);
          var_eq.ApplyAdamping (dTMBase, Adamping_dphi, *this);
          
          // We multiply LocalJacobNL with dTMBase
          // Write an unlooped version :
          //MltVect(LocalJacobNL,dTMBase,NLpart);
          for (int k1 = 0; k1 < dimension; k1++)
          {
            for (int g = 0; g < nb_quad; g++)
            {
              //for (int j=0; j<dimension; j++) {
              NLpart(k1, g) = LocalJacobNL(k1, k2, g)*dTMBase(k2, g);
              //}
            }
          }
    
          vect_contre_phi  = invDeltat2*Mphi 
                  + invDeuxDeltat* Rphi + thetastiff * ( Cphi+Ddphi);// thetastiff * Bdphi +
    
          vect_contre_dphi = theta* Adphi + thetastiff* ( AStiffdphi + Bphi)
                  + coef_NL*NLpart + invDeuxDeltat*Adamping_dphi;
    
          res_int_phi.Zero();
          res_int_dphi.Zero();
          this->CalculeIntPhi  (vect_contre_phi  ,res_int_phi );
          this->CalculeIntGrad (vect_contre_dphi ,res_int_dphi);

          // We copy that in the DiffMatrix.
          int jglob2 = this->mesh.GetNumberDof(ne, j2);
          for (int j1 = 0; j1 < nb_base; j1++)
          {
            int jglob1 = this->mesh.GetNumberDof(ne, j1);
            for (int k1 = 0; k1 < dimension; k1++)
            {
              Real_wp value = res_int_phi(k1, j1) + res_int_dphi(k1, j1);
              var.AddDiffStringString(i, k1, jglob1, value, i, k2, jglob2, DiffMatrix);
            } // for k1 (component number of the basis function)
          }// for j1 (local ddl basis function )

        }// for k2 ( component number of up)
      }// for j2 (local ddl number of up)
    }// for ne

    // ****************************
    // Treat Boundary Conditions at the Bridge
    // ****************************
    
    int nodl = this->mesh.GetNbDof()-1;
    
    if (is_attached)
    {
      for (int k = 0; k < dimension; k++)
      {
        if (var_eq.Dirichlet_Bridge(k))
        {
          // We set dirichlet condition 
          var.SetDiffStringString(i, k, nodl, 1.0, i, k, nodl, DiffMatrix);
          for (int j = 1; j < nb_base; j++)
            var.SetDiffStringString(i, k, nodl, 0.0, i, k, nodl-j, DiffMatrix);
    
          // hence it does not depend on other components
          for(int k2 = 0; k2 < dimension; k2++)
            if(k2 != k)
              for (int j = 0; j < nb_base; j++)
                var.SetDiffStringString(i, k, nodl, 0.0, i, k2, nodl-j, DiffMatrix);
    
        }
      }
    }
    
    // ****************************
    // Treat LM Part
    // ****************************
    
    for (int k = 0; k < dimension; k++)
    {
      if (var_eq.Dirichlet_Bridge(k))
      {
        if (is_attached)
        {
          var.SetDiffLMLM(i, 1.0, i, DiffMatrix);
        }
        else
        {
          for (int p = 0; p < var.nb_lagrange_mult; p++)
          {
            int dof_string = var.GetStringDofNumber(i, k, nodl);
            int dof_LM = var.GetLagrangeM_DofNumber(i, p);
            DiffMatrix.SetEntry(dof_string, dof_LM, -transmission_chevalet(p, k));
          }
        }
      }
    }
    
    // ****************************
    // Treat Hammer Part 
    // ****************************
    
    // on derive la partie - dt \sum F_i par rapport a l'inconnue force_hammer
    if (!var.drop_hammer)
    {
      if ((var.hammer.is_interacting) && (!var.shank.is_interacting_string))
        DiffMatrix.AddInteraction(var.GetHammerDofNumber(),
                                    var.offset_force_hammer+i, -var.Deltat);
    }
    
    Vector<Real_wp> der_interaction_hammer;
    Real_wp interaction_hammer = 0.0;
    Real_wp interaction_hammer_aux = 0.0;
    // derivee de Fi par rapport a son argument y = <u>-xi
    der_interaction_hammer = var.ComputeHammerInteractionDerivative(*this, val);
    interaction_hammer = der_interaction_hammer(0);

    // derivee de l'equation force_marteau - interaction_marteau = 0
    // par rapport a force_marteau
    if ((var.hammer.is_interacting) && (!var.drop_hammer))
    {
      DiffMatrix.AddInteraction(var.offset_force_hammer+i, var.offset_force_hammer+i, 1.0);
      DiffMatrix.AddInteraction(var.offset_force_hammer+i,
                                  var.GetHammerDofNumber(), interaction_hammer);

      if (var.hammer.with_auxiliary_variable)
      {
        interaction_hammer_aux = der_interaction_hammer(1);
        DiffMatrix.AddInteraction(var.offset_force_hammer+i,
                                  var.GetHammerAuxDofNumber(i), interaction_hammer_aux);

        var.AddDiffHammerAuxHammer(i, -interaction_hammer_aux, DiffMatrix);
        var.AddDiffHammerAuxHammerAux(i, 1.0, i, DiffMatrix);
      }
    }
    
    // derivee par rapport a la corde
    if ((is_struck) && (!var.drop_hammer))
    {
      for (int jj = 0; jj < this->mesh.GetNbDof(); jj++)
      {
        Real_wp repj = this->repartition_marteau(jj);
        if (abs(repj) > 100.0*epsilon_machine)
        {
          if (var.hammer.is_interacting)
          {
            if (var.hammer.strike_at_an_angle)
            {
              DiffMatrix.AddInteraction(var.offset_force_hammer+i, var.GetStringDofNumber(i, 0, jj),
                      -interaction_hammer*repj*cos(var.hammer.strike_angle*pi_wp/180.0));
              DiffMatrix.AddInteraction(var.offset_force_hammer+i, var.GetStringDofNumber(i, 1, jj),
                      -interaction_hammer*repj*sin(var.hammer.strike_angle*pi_wp/180.0));
              DiffMatrix.AddInteraction(var.GetStringDofNumber(i, 0, jj), var.offset_force_hammer+i,
                                        repj*cos(var.hammer.strike_angle*pi_wp/180.0));
              DiffMatrix.AddInteraction(var.GetStringDofNumber(i, 1, jj), var.offset_force_hammer+i,
                                        repj*sin(var.hammer.strike_angle*pi_wp/180.0));
              
            }
            else
            {
              DiffMatrix.AddInteraction(var.offset_force_hammer+i, var.GetStringDofNumber(i, 0, jj),
                      -interaction_hammer*repj);
              DiffMatrix.AddInteraction(var.GetStringDofNumber(i, 0, jj), var.offset_force_hammer+i,
                                          repj);
            }
            if (var.hammer.with_auxiliary_variable)
            {
              if (var.hammer.strike_at_an_angle)
              {
                DiffMatrix.AddInteraction(var.GetHammerAuxDofNumber(i), var.GetStringDofNumber(i, 0, jj),
                      interaction_hammer_aux*repj*cos(var.hammer.strike_angle*pi_wp/180.0));
                DiffMatrix.AddInteraction(var.GetHammerAuxDofNumber(i), var.GetStringDofNumber(i, 1, jj),
                      interaction_hammer_aux*repj*sin(var.hammer.strike_angle*pi_wp/180.0));
              }
              else
              {
                DiffMatrix.AddInteraction(var.GetHammerAuxDofNumber(i), var.GetStringDofNumber(i, 0, jj),
                      interaction_hammer_aux*repj);
              }
            }
          }
        }
      }
    }
    
    // ****************************
    // Treat Boundary Conditions at the Agraffe
    // ****************************
    
    for (int k = 0; k < dimension; k++)
    {
      if (var_eq.Dirichlet_Agraffe(k))
      {
        // We set dirichlet condition on the first point
        DiffMatrix.ClearRow(var.GetStringDofNumber(i, k, 0));
        var.SetDiffStringString(i, k, 0, 1.0, i, k, 0, DiffMatrix);
      }
    }
    
    return interaction_hammer;
  }


  //! Calcul de la jacobienne DF(X)
  /*!
    \param[inout] DiffMatrix differential matrix
    \param[in] val solution X for which the differential matrix is computed
    \param[in] num_xsi dof number for the position of the shank (if a shank is present)
    \param[in] var multistring object
  */
  template<int nb_base, int nb_quad, class TypeEquation>
  template<class GenericPb> 
  void PianoString<nb_base, nb_quad, TypeEquation>
  ::UpdateDiff(VirtualMatrix<Real_wp> & DiffMatrix, Matrix<Real_wp> & val, GenericPb & var)
  {
    int i = numero_string;
    //Real_wp invDeltat2 = var.invDeltat2;
    //Real_wp invDeuxDeltat = var.invDeuxDeltat;

    // ****************************
    // Treat Hammer Part 
    // ****************************
    
    Vector<Real_wp> der_interaction_hammer;
    Real_wp interaction_hammer = 0.0;
    Real_wp interaction_hammer_aux = 0.0;
    // derivee de Fi par rapport a son argument xi
    der_interaction_hammer = var.ComputeHammerInteractionDerivative(*this, val);
    interaction_hammer = der_interaction_hammer(0);
    
    // derivee de l'equation force_marteau - interaction_marteau = 0
    // par rapport a force_marteau
    if ((var.hammer.is_interacting) && (!var.drop_hammer))
    {
      DiffMatrix.SetEntry(var.offset_force_hammer+i,
                                  var.GetHammerDofNumber(), interaction_hammer);

      if (var.hammer.with_auxiliary_variable)
      {
        interaction_hammer_aux = der_interaction_hammer(1);
        DiffMatrix.SetEntry(var.offset_force_hammer+i,
                                  var.GetHammerAuxDofNumber(i), interaction_hammer_aux);

        var.SetDiffHammerAuxHammer(i, -interaction_hammer_aux, DiffMatrix);
      }
    }
    
    // derivee par rapport a la corde
    if ((is_struck) && (!var.drop_hammer))
    {
      for (int jj = 0; jj < this->mesh.GetNbDof(); jj++)
      {
        Real_wp repj = this->repartition_marteau(jj);
        if (abs(repj) > 100.0*epsilon_machine)
        {
          if (var.hammer.is_interacting)
          {
            if (var.hammer.strike_at_an_angle)
            {
              DiffMatrix.SetEntry(var.offset_force_hammer+i, var.GetStringDofNumber(i, 0, jj),
                      -interaction_hammer*repj*cos(var.hammer.strike_angle*pi_wp/180.0));
              DiffMatrix.SetEntry(var.offset_force_hammer+i, var.GetStringDofNumber(i, 1, jj),
                      -interaction_hammer*repj*sin(var.hammer.strike_angle*pi_wp/180.0));
            }
            else
            {
              DiffMatrix.SetEntry(var.offset_force_hammer+i, var.GetStringDofNumber(i, 0, jj),
                      -interaction_hammer*repj);
            }
            if (var.hammer.with_auxiliary_variable)
            {
              if (var.hammer.strike_at_an_angle)
              {
                DiffMatrix.SetEntry(var.GetHammerAuxDofNumber(i), var.GetStringDofNumber(i, 0, jj),
                      interaction_hammer_aux*repj*cos(var.hammer.strike_angle*pi_wp/180.0));
                DiffMatrix.SetEntry(var.GetHammerAuxDofNumber(i), var.GetStringDofNumber(i, 1, jj),
                      interaction_hammer_aux*repj*sin(var.hammer.strike_angle*pi_wp/180.0));
              }
              else
              {
                DiffMatrix.SetEntry(var.GetHammerAuxDofNumber(i), var.GetStringDofNumber(i, 0, jj),
                      interaction_hammer_aux*repj);
              }
            }
          }
        }
      }
    }
  }
  
  
  //! computes linear matrix applied to Un and Unm1
  template<int nb_base, int nb_quad, class TypeEquation>
  template<class GenericPb>
  void PianoString<nb_base, nb_quad, TypeEquation>
  ::ComputeLinearMatrix(GenericPb& var)
  {
    Real_wp theta = var.theta(0);
    Real_wp thetastiff = var.theta(1);
    Real_wp invDeltat2 = var.invDeltat2;
    Real_wp invDeuxDeltat = var.invDeuxDeltat;
    
    // ****************************
    // Treat Major Scheme Part (FV)
    // ****************************
    
    TinyMatrix<Real_wp, General, dimension, nb_quad> TMBase, dTMBase;
    TinyMatrix<Real_wp, General, dimension, nb_quad> Mphi, Rphi, Adphi, AStiffdphi,
      Adamping_dphi, Bphi, Cphi, Ddphi;
    
    TinyMatrix<Real_wp, General, dimension, nb_quad> vect_contre_phi_un, vect_contre_dphi_un,
      vect_contre_phi_unm1, vect_contre_dphi_unm1, vect_contre_phi_unp1, vect_contre_dphi_unp1;
    
    TinyMatrix<Real_wp, General, dimension, nb_base>
      res_int_phi_un, res_int_dphi_un, res_int_phi_unm1,
      res_int_dphi_unm1, res_int_phi_unp1, res_int_dphi_unp1;

    for (int k2 = 0; k2 < dimension; k2++)
    {
      for (int j2 = 0; j2 < nb_base; j2++)
      {
        // it will be the column indices
            
        // We need an elementary matrix
        TMBase.Zero();
        dTMBase.Zero();
        for (int g = 0; g < nb_quad; g++)
        {
          TMBase(k2, g) = this->valPhi(j2, g);
          dTMBase(k2, g) = this->valdPhi(j2, g)*this->invDeltax;
        }
        
        // We apply every linear part
        var_eq.ApplyM (TMBase, Mphi, *this);
        var_eq.ApplyR (TMBase, Rphi, *this);
        var_eq.ApplyAGrad (dTMBase, Adphi, *this);
        var_eq.ApplyAStiffGrad (dTMBase, AStiffdphi, *this);

        var_eq.ApplyB (TMBase, Bphi, *this);
        var_eq.ApplyC (TMBase, Cphi, *this);
        var_eq.ApplyD (dTMBase, Ddphi, *this);
        var_eq.ApplyAdamping (dTMBase, Adamping_dphi, *this);
        
        vect_contre_phi_unm1  = invDeltat2*Mphi - invDeuxDeltat*Rphi + thetastiff*(Cphi + Ddphi);
        
        vect_contre_dphi_unm1 = theta*Adphi + thetastiff*(AStiffdphi + Bphi) - invDeuxDeltat*Adamping_dphi;
        
        vect_contre_phi_un  = -2.0*invDeltat2*Mphi + (1.0-2.0*thetastiff) * (Cphi + Ddphi);
        vect_contre_dphi_un = (1.0-2.0*theta)*Adphi + (1.0-2.0*thetastiff)*(AStiffdphi + Bphi);

        vect_contre_phi_unp1  = invDeltat2*Mphi + invDeuxDeltat*Rphi + thetastiff*(Cphi + Ddphi);
        
        vect_contre_dphi_unp1 = theta*Adphi + thetastiff*(AStiffdphi + Bphi) + invDeuxDeltat*Adamping_dphi;
        
        res_int_phi_unm1.Zero();
        res_int_dphi_unm1.Zero();
        this->CalculeIntPhi  (vect_contre_phi_unm1, res_int_phi_unm1);
        this->CalculeIntGrad (vect_contre_dphi_unm1, res_int_dphi_unm1);
        
        res_int_phi_un.Zero();
        res_int_dphi_un.Zero();
        this->CalculeIntPhi  (vect_contre_phi_un, res_int_phi_un);
        this->CalculeIntGrad (vect_contre_dphi_un, res_int_dphi_un);

        res_int_phi_unp1.Zero();
        res_int_dphi_unp1.Zero();
        this->CalculeIntPhi  (vect_contre_phi_unp1, res_int_phi_unp1);
        this->CalculeIntGrad (vect_contre_dphi_unp1, res_int_dphi_unp1);

        // We copy that in the DiffMatrix.
        for (int j1 = 0; j1 < nb_base; j1++)
        {
          for (int k1 = 0; k1 < dimension; k1++)
          {
            //int num_col = k2*nb_base + j2;
            //int num_row = k1*nb_base + j1;
            
            Real_wp value = res_int_phi_unm1(k1, j1) + res_int_dphi_unm1(k1, j1);
            mat_Unm1(k1, k2)(j1, j2) = value;
            
            //var.AddDiffStringString(i, k1, jglob1, value, i, k2, jglob2, A);

            value = res_int_phi_un(k1, j1) + res_int_dphi_un(k1, j1);
            mat_Un(k1, k2)(j1, j2) = value;

            value = res_int_phi_unp1(k1, j1) + res_int_dphi_unp1(k1, j1);
            mat_Unp1(k1, k2)(j1, j2) = value;
            
            //var.AddDiffStringString(i, k1, jglob1, value, i, k2, jglob2, B);
          } // for k1 (component number of the basis function)
        }// for j1 (local ddl basis function )
            
      }// for k2 ( component number of up)
    }// for j2 (local ddl number of up)    
  }  

  
  //! Calcule l'energie associee a la corde
  /*!
    \param[in] var objet multistring
    Retourne l'energie de la corde
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  template<class GenericPb> 
  Real_wp PianoString<nb_base, nb_quad, TypeEquation>
  ::ComputeEnergy(GenericPb & var, Real_wp& ener_u, Real_wp& ener_v)
  {
    Real_wp energy = 0.0; ener_u = 0.0; ener_v = 0.0;
    
    Real_wp theta = var.theta(0);
    Real_wp thetastiff = var.theta(1);
    Real_wp Deltat = var.Deltat;
    Real_wp Deltat2 = Deltat*Deltat;
    Real_wp invDeltat = var.invDeltat;
    
    int i = numero_string;
    
   
    // ****************************
    // Treat Major Scheme Part (FV)
    // ****************************
    
    TinyMatrix<Real_wp, General, dimension, nb_quad> TMplus, TM, dTMplus, dTM;
    TinyMatrix<Real_wp, General, dimension, nb_quad> TimeDerivative, TimeAverage,
      dTimeDerivative, dTimeAverage;
    
    TinyMatrix<Real_wp, General, dimension, nb_quad> MTimeDerivative, ATimeDerivative,
      AStiffTimeDerivative, BTimeDerivative;
    
    TinyMatrix<Real_wp, General, dimension, nb_quad>
      CTimeDerivative, DTimeDerivative, AStiffTimeAverage;
    TinyMatrix<Real_wp, General, dimension, nb_quad> BTimeAverage, CTimeAverage, DTimeAverage;
    TinyMatrix<Real_wp, General, dimension, nb_quad> ATimeAveragePlus, ATimeAverageCurrent,
      ATimeAverageHalf;
    TinyMatrix<Real_wp, General, dimension, nb_quad> vect_contre_TA, vect_contre_dTA,
      vect_contre_TD, vect_contre_dTD;
    TinyVector<Real_wp, nb_base> res_dot_prod_u, res_dot_prod_v, UPlus, UCurrent, Uhalf;
    
    // boucle sur les elements
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
      { 
  res_dot_prod_u.Zero(); res_dot_prod_v.Zero();
  
  // We fill the TM with values at quad points.
  this->ComputeTM (i, ne, 1, var.val, var, TMplus);
  this->ComputeTM (i, ne, 0, var.val, var, TM);
  
  this->ComputedTM (i, ne, 1, var.val, var, dTMplus);
  this->ComputedTM (i, ne, 0, var.val, var, dTM);
  // We prep the tinymatrices
  TimeDerivative = (TMplus - TM)* invDeltat;
  TimeAverage = (TMplus + TM)*Real_wp(0.5);
  dTimeDerivative = (dTMplus - dTM)* invDeltat;
  dTimeAverage = (dTMplus + dTM)*Real_wp(0.5);
  
  // We apply the linear terms
  var_eq.ApplyM (TimeDerivative, MTimeDerivative, *this);
  
  var_eq.ApplyAGrad (dTimeDerivative, ATimeDerivative, *this);
  
  var_eq.ApplyAStiffGrad (dTimeDerivative, AStiffTimeDerivative, *this);
  //  var_eq.ApplyB (dTimeDerivative, BTimeDerivative, *this);
  var_eq.ApplyB (TimeDerivative, BTimeDerivative, *this);
  var_eq.ApplyC (TimeDerivative, CTimeDerivative, *this);
  var_eq.ApplyD (dTimeDerivative, DTimeDerivative, *this);
  
  var_eq.ApplyAStiffGrad (dTimeAverage, AStiffTimeAverage, *this);
  //  var_eq.ApplyB (dTimeAverage, BTimeAverage, *this);
  var_eq.ApplyB (TimeAverage, BTimeAverage, *this);
  var_eq.ApplyC (TimeAverage, CTimeAverage, *this);
  var_eq.ApplyD (dTimeAverage, DTimeAverage, *this);
  
  // Depending on the NL part
  Real_wp gamma = 1.0;
  if   (var.solver == var.THETA_NL_HALF)
    {
      gamma = 1.0;
      var_eq.ApplyAGrad(dTMplus, ATimeAveragePlus, *this);
      var_eq.ApplyAGrad(dTM, ATimeAverageCurrent, *this);
      
      TypeEquation::DotProdColUV(ATimeAveragePlus, dTMplus, res_dot_prod_u, res_dot_prod_v);
      TypeEquation::DotProdColUV(ATimeAverageCurrent, dTM, res_dot_prod_u, res_dot_prod_v);
      res_dot_prod_u *= 0.5; res_dot_prod_v *= 0.5;
      
      var_eq.ComputeU (dTMplus, UPlus, *this);
      var_eq.ComputeU (dTM, UCurrent, *this);
      Uhalf = (UPlus + UCurrent)*Real_wp(0.5);
    } 
  else if (var.solver == var.THETA_NL_QUARTER)
    {
      gamma = 0.5;
      var_eq.ApplyAGrad (dTimeAverage, ATimeAverageHalf, *this); 
      TypeEquation::DotProdColUV(ATimeAverageHalf, dTimeAverage, res_dot_prod_u, res_dot_prod_v);
      var_eq.ComputeU (dTimeAverage, Uhalf, *this);
    }
  
  res_dot_prod_u += Real_wp(2)*Uhalf;
  
  Real_wp coef = Deltat2*(theta-0.5*gamma);
  Real_wp coefstiff = Deltat2*(thetastiff-0.25);
  
  vect_contre_TA = (CTimeAverage+ DTimeAverage);
  vect_contre_dTA = BTimeAverage +  AStiffTimeAverage ;
  
  vect_contre_TD = (MTimeDerivative + coefstiff* ( CTimeDerivative + DTimeDerivative));
        // BTimeDerivative +
  vect_contre_dTD = coef* ATimeDerivative
          + coefstiff*(BTimeDerivative + AStiffTimeDerivative);
  
  TypeEquation::DotProdColUV(vect_contre_TA, TimeAverage, res_dot_prod_u, res_dot_prod_v);
  
  TypeEquation::DotProdColUV(vect_contre_dTA, dTimeAverage, res_dot_prod_u, res_dot_prod_v);
  
  TypeEquation::DotProdColUV(vect_contre_TD, TimeDerivative, res_dot_prod_u, res_dot_prod_v);
  
  TypeEquation::DotProdColUV(vect_contre_dTD, dTimeDerivative, res_dot_prod_u, res_dot_prod_v);
  
  ener_u += DotProd(res_dot_prod_u, this->Weights);
        ener_v += DotProd(res_dot_prod_v, this->Weights);
      }
    
    ener_u *= this->Deltax*0.5;
    ener_v *= this->Deltax*0.5;
    energy = ener_u + ener_v;
    
    return energy;  
  }
  
}

#define MONTJOIE_FILE_PIANO_STRING_CXX
#endif
