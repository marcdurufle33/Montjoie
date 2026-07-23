#ifndef MONTJOIE_FILE_HAMMER_CXX

namespace Montjoie
{
  
  //! default constructor
  Hammer::Hammer()
  {
    is_interacting = true;
    type_interaction = COUPLING;
    delta = 0.1;
    masse_marteau = 20e-3;
    point_impact_marteau = 0.2;
    vitesse_marteau = 0.0;
    eloignement_marteau = -10.0;
    epaisseur_marteau = 0.02;
    pente_marteau = 2000;
    exposant = 2.1;
    DeltatForce = 0.0;
    epsilon = 1e-6;
    nb_base_time = 10;
    energy_free_hammer = 0.0;
    sin_amplitude = 1;
    sin_freq = 50;
    strike_angle = 0.0;
    strike_at_an_angle = false;
    with_auxiliary_variable = true;

    contact_model_with_string = ORIGINAL;
    regularised_contact_model = false;
  }
  
  
  //! modification des parametres du marteau avec une ligne du fichier de donnees
  /*!
    \param[in] description_field mot-cle de la ligne (par exemple Hammer = ...)
    \param[in] parameters liste des parametres associes au mot-cle
  */
  void Hammer::SetInputData(const string & description_field, const Vector<string> &parameters)
  {
    
    if (!description_field.compare("Hammer"))
    {
      if (parameters.GetM() <= 0)
      {
        cout << "In SetInputData of Hammer" << endl;
        cout << "Hammer needs at least one parameter, for instance :" << endl;
        cout << "Hammer = NO" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      if (!parameters(0).compare("YES"))
      {
        if (parameters.GetM() <= 9)
        {
          cout << "In SetInputData of Hammer" << endl;
          cout << "Hammer needs 9 parameters, for instance :" << endl;
          cout << "Hammer = YES masse point_impact vitesse exposant "
               << "epaisseur pente delta eloignement epsilon" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
      
        is_interacting = true;
        type_interaction = COUPLING;

        masse_marteau = to_num<Real_wp>(parameters(1));
        point_impact_marteau = to_num<Real_wp>(parameters(2));
        vitesse_marteau = to_num<Real_wp>(parameters(3));
        exposant = to_num<Real_wp>(parameters(4));
        epaisseur_marteau = to_num<Real_wp>(parameters(5));
        pente_marteau = to_num<Real_wp>(parameters(6));
        delta = to_num<Real_wp>(parameters(7));
        eloignement_marteau = to_num<Real_wp>(parameters(8));
        epsilon = to_num<Real_wp>(parameters(9));
      
      
        if(parameters.GetM() > 10)
        {
          if (parameters.GetM() <= 11)
          {
            cout << "In SetInputData of Hammer" << endl;
            cout << "Hammer needs 11 parameters, for instance :" << endl;
            cout << "Hammer = YES masse point_impact vitesse exposant "
                 << "epaisseur pente delta eloignement epsilon ANGLE angle" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }     
      
          if (!parameters(10).compare("ANGLE"))
          {
            strike_angle = to_num<Real_wp>(parameters(11));
            strike_at_an_angle = true;
          }
        }
      }

      else if (!parameters(0).compare("NO"))
      {
        is_interacting = false;
        type_interaction = NONE;
      }
      else if (!parameters(0).compare("FORCE"))
      {
        if (parameters.GetM() <= 1)
        {
          cout << "In SetInputData of Hammer" << endl;
          cout << "Hammer needs at least two parameters, for instance :" << endl;
          cout << "Hammer = FORCE RICKER" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
        
        is_interacting = false;
        if (!parameters(1).compare("EXPERIMENT"))
        {
          if (parameters.GetM() <= 6)
          {
            cout << "In SetInputData of Hammer" << endl;
            cout << "Hammer needs 7 parameters, for instance :" << endl;
            cout << "Hammer = FORCE EXPERIMENT file point_impact epaisseur "
              << "pente DeltaForce" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
      
          type_interaction = EXPERIMENT;
          string name_force_file(parameters(2));
          ForceMarteauExt.ReadText(name_force_file);
          cout << "Force coming from experiment"<<endl;
          point_impact_marteau  = to_num<Real_wp>(parameters(3));
          epaisseur_marteau = to_num<Real_wp>(parameters(4));
          pente_marteau = to_num<Real_wp>(parameters(5));
          DeltatForce = to_num<Real_wp>(parameters(6));
          nb_base_time = 2;
        }
        else if (!parameters(1).compare("RICKER"))
        {
          type_interaction = RICKER;
          cout<<"Ricker force not implemented yet. Feel free to do it ! " <<endl;
          abort();
          // a ecrire
        }
        else if (!parameters(1).compare("SINUS"))
        {
          type_interaction = SINUS;
          cout<<"Sinus force at: " << parameters(3) << "Hz" <<endl;
          sin_amplitude = to_num<Real_wp>(parameters(2)); //amplitude
          sin_freq = to_num<Real_wp>(parameters(3)); //frequency
          point_impact_marteau = to_num<Real_wp>(parameters(4));
          // delta = 1;
          // abort();
          // a ecrire
        }
        else if (!parameters(1).compare("CHIRP"))
        {
          type_interaction = CHIRP;
          cout << "Chirping at" << endl;
          chirp_amplitude = to_num<Real_wp>(parameters(2)); //amplitude
          chirp_fstart = to_num<Real_wp>(parameters(3)); //f_start
          chirp_fend = to_num<Real_wp>(parameters(4)); //f_end
          chirp_phase = to_num<Real_wp>(parameters(5)); //phase
          point_impact_marteau = to_num<Real_wp>(parameters(6));
        }       
      }
      else
      { 
      }
    }
    else if (!description_field.compare("HammerStiffness"))
    {
      if (!parameters(0).compare("ORIGINAL"))
      {
        if (parameters.GetM() <= 2)
        {
          cout << "In SetInputData of Hammer" << endl;
          cout << "HammerStiffness needs 2 parameters, for instance :" << endl;
          cout << "HammerStiffness = ORIGINAL K R" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
    
        Ki.PushBack(to_num<Real_wp>(parameters(1)));
        Ri.PushBack(to_num<Real_wp>(parameters(2)));
      }
      else if (!parameters(0).compare("ROBOTRAN"))
      {
        if (parameters.GetM() <= 3)
        {
          cout << "In SetInputData of Hammer" << endl;
          cout << "HammerStiffness needs 2 parameters, for instance :" << endl;
          cout << "HammerStiffness = ROBOTRAN A B C D" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
    
        Ai.PushBack(to_num<Real_wp>(parameters(1)));
        Bi.PushBack(to_num<Real_wp>(parameters(2)));
        Ci.PushBack(to_num<Real_wp>(parameters(3)));
        Di.PushBack(to_num<Real_wp>(parameters(4)));
        contact_model_with_string = ROBOTRAN;
      }
    }
    else if (!description_field.compare("RegulariseHammerStringContact"))
    {
      if (parameters.GetM() <= 0)
      {
        cout << "In SetInputData of Hammer" << endl;
        cout << "RegulariseHammerStringContact needs 1 parameters, for instance :" << endl;
        cout << "RegulariseHammerStringContact = YES" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
      if (parameters(0) == "YES")
      {
        regularised_contact_model = true;
      }
    }
  }
  
  
  //! setting hammer parameters with a line of the data file
  /*!
    \param[in] description_field keyword of the line
    \param[in] parameters parameters associated with the keyword
    \param[in] string_p string struck by the current hammer
    \param[in] param_corde parameters of all the strings of the piano
  */
  template<int p, int q, class TypeEq>
  void Hammer::SetInputData(const string & description_field, const Vector<string> &parameters,
     PianoString<p, q, TypeEq>& string_p, ParameterMultistring& param_corde)
  {
    
    if (!description_field.compare("HammerStiffness")) 
      SetInputData(description_field, parameters);
    
    if (!description_field.compare("Hammer")) 
    {
      if ((!parameters(0).compare("YES")) || (!parameters(0).compare("NO"))
       || (!parameters(0).compare("FORCE")) )
      {
        SetInputData(description_field, parameters);
      }
      else
      {
        if (parameters.GetM() <= 1)
        {
          cout << "In SetInputData of Hammer" << endl;
          cout << "Hammer needs two parameters, for instance :" << endl;
          cout << "Hammer = Note vitesse" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
      
        // on considere que c'est une note (e.g. C2)
        Real_wp Ks, Rs;
        param_corde.GetParamHammer(parameters(0), masse_marteau, exposant, Ks, Rs);
        point_impact_marteau = 0.12*string_p.L;
        vitesse_marteau = to_num<Real_wp>(parameters(1));              
        if (parameters.GetM() > 6)
        {
          epaisseur_marteau = to_num<Real_wp>(parameters(2));
          pente_marteau = to_num<Real_wp>(parameters(3));
          delta = to_num<Real_wp>(parameters(4));
          eloignement_marteau = to_num<Real_wp>(parameters(5));
          epsilon = to_num<Real_wp>(parameters(6));
        }
        else
        {
          epaisseur_marteau = 0.02;
          pente_marteau = 2011;
          delta = 0.02;
          eloignement_marteau = -0.02001;
          epsilon = 1e-10;
        }
      
        Ki.Fill(Ks);
        Ri.Fill(Rs);
      }
    } 
  }
  
  
  //! displaying the parameters of the hammer
    /*!
      \param[inout] file_out output stream where the parameters are written
      \param[in] write_on_file if true, the parameters are written in the output stream
    */
  void Hammer::Display(ostream& file_out, bool write_on_file)
  {
    cout << " -------- Hammer Parameters: -------- "<<endl;
    DISP(type_interaction);
    DISP(masse_marteau);
    DISP(point_impact_marteau);
    DISP(vitesse_marteau);
    DISP(exposant);
    DISP(epaisseur_marteau);
    DISP(pente_marteau);
    if (contact_model_with_string == ORIGINAL)
    {
      DISP(Ki);
      DISP(Ri);
    }
    else
    {
      DISP(Ai);
      DISP(Bi);
      DISP(Ci);
      DISP(Di);
    }
    DISP(delta);
    DISP(eloignement_marteau);
    DISP(epsilon); 
    cout << " ----------------------------------- "<<endl;
    
    if (write_on_file)
    {
      file_out << " -------- Hammer Parameters: -------- "<<endl;
      if (type_interaction == COUPLING)
        file_out << "Hammer is coupled with the string" << endl;
      else if (type_interaction == EXPERIMENT)
        file_out << "The force is given by datas of a file" << endl;
      else
        file_out << "Unknown hammer" << endl;
  
      file_out << "masse_marteau = " << masse_marteau << endl;
      file_out << "point_impact_marteau = " << point_impact_marteau << endl;
      file_out << "vitesse_marteau = " << vitesse_marteau << endl;
      file_out << "exposant = " << exposant << endl;
      file_out << "epaisseur_marteau = " << epaisseur_marteau << endl;
      file_out << "pente_marteau = " << pente_marteau << endl;
      file_out << "stiffness_Kmarteau = " << Ki << endl;
      file_out << "damping_Rmarteau = " << Ri << endl;
      file_out << "d0_marteau = " << delta << endl;
      file_out << "eloignement_marteau_corde" << eloignement_marteau << endl;
      file_out << " ----------------------------------- " << endl << endl;
    }
  }
  
  
  //! setting initial condition of the hammer (from initial position and velocity)
  template<class GenericPb, int p, int q, class TypeEq>
  void Hammer::SetInitialConditions(GenericPb & var, Vector<PianoString<p, q, TypeEq> > & piano_strings)
  {
    var.SetHammer(0, eloignement_marteau, var.val); //position initiale marteau
    var.SetHammer(1, eloignement_marteau + var.Deltat*vitesse_marteau, var.val);
    var.SetHammer(2, eloignement_marteau + 2.0*var.Deltat*vitesse_marteau, var.val);

    if (with_auxiliary_variable)
    {
      for(int i = 0; i < piano_strings.GetM(); i++) // on initialise les variables aux pour chaque corde
      {
        // on calcule le deplacement initial moyen de la corde au niveau de l'impact du marteau
        Real_wp ux0 = 0.0;
        Real_wp ux1 = 0.0;
        for (int j = 0; j < piano_strings(i).repartition_marteau_sparse.GetM(); j++)
        {
          int jj = piano_strings(i).repartition_marteau_sparse.Index(j);
          Real_wp repj = piano_strings(i).repartition_marteau_sparse.Value(j);
          ux0 += var.GetValue(i, 0, jj, 0, var.val)*repj;
          ux1 += var.GetValue(i, 0, jj, 1, var.val)*repj; 
        }

        var.SetHammerAux(i, 0, 0, var.val);
        var.SetHammerAux(i, 1, GAux(ux0 - var.GetHammer(0, var.val), i), var.val);
        var.SetHammerAux(i, 2, GAux(ux1 - var.GetHammer(1, var.val), i), var.val);
      }
    }
  }
  
  
  //! computation of the energy of the hammer
  /*!
    \param[in] var object multistring
    \return the kinetic energy of the hammer
  */
  template<class GenericPb>
  Real_wp Hammer::ComputeKineticEnergy(GenericPb & var)
  {
    if (var.drop_hammer)
      return energy_free_hammer;
    
    // m/2*[ (xi^{n+1}-xi{n})/Deltat ]^2
    Real_wp res = 0.0;
    res = (var.GetHammer(1, var.val)-var.GetHammer(0, var.val))
      *(var.GetHammer(1, var.val)-var.GetHammer(0, var.val));
    
    res *= masse_marteau*0.5*var.invDeltat2;
    return res;
  }


  //! computation of the potential energy of the hammer
  /*!
    \param[in] var object multistring
    \param[in] i string number
    \param[in] ux0p string vertical displacement (after)
    \param[in] ux0m string vertical displacement (before)
    \return the potential energy of the hammer
  */
  template<class GenericPb>
  Real_wp Hammer::ComputePotentialEnergy(GenericPb & var, int i, Real_wp ux0p, Real_wp ux0m)
  {
    Real_wp res = 0.0;
    if (with_auxiliary_variable)
    {
      Real_wp var_aux = var.GetHammerAux(i, 1, var.val);
      res = 0.5 * var_aux * var_aux;
    }
    else
    {
      // K_i/2 * [ \Psi(U_i^{n+1} - xi^{n+1}) + \Psi(U_i^{n} - xi^{n}) ]
      res = Psi(ux0p - var.GetHammer(1, var.val)) + Psi(ux0m - var.GetHammer(0, var.val));
      res *= Ki(i)*0.5;
    }
    return res;
  }
  
  
  //! computation of the force on all time steps in the case where
  //! the force is read from a file
  template<class GenericPb>
  void Hammer::ComputeInterpolatedForce(GenericPb & var)
  {
    if (var.drop_hammer)
      return;
    
    if (type_interaction == EXPERIMENT)
    {
      Real_wp t_begin = var.t_begin;
      Real_wp t_end = var.t_end;
  
      int nb_force = toInteger(ceil(abs(t_end - t_begin)/DeltatForce ));
      int taille = ForceMarteauExt.GetM();
      ForceMarteauExt.Resize(nb_force+1);
      for (int p = taille; p < nb_force+1; p++)
        ForceMarteauExt(p) = 0.0;
      
      string name_file = var.DOSSIER + "ForceMarteauExt.txt";
      remove(name_file.data());
      ForceMarteauExt.WriteText(name_file);
      
      mesh_time.CreateRegularMesh(t_begin,t_end,nb_force+1,1);
      mesh_time.SetOrder(nb_base_time-1);
      mesh_time.NumberMesh(false);
      base_time.ConstructQuadrature(nb_base_time-1,Globatto<Real_wp>::QUADRATURE_LOBATTO); 
      
      grid_interp_time.Init(t_begin, t_end, var.nb_max_iter+1);
      // en utilisant le grid_interp de montjoie
      grid_interp_time.LocalizePoints(mesh_time);
      
      ForceMarteau.Reallocate(var.nb_max_iter+1);
      ForceMarteau.Fill(Real_wp(0));
      for (int l = 0; l < grid_interp_time.GetNbPointsGrid(); l++)
      {
        int ne = grid_interp_time.GetElementNumber(l);
        Real_wp x = grid_interp_time.GetLocalCoordinate(l);
          
        for (int jj = 0; jj < nb_base_time; jj++)
        {
          int j = mesh_time.GetNumberDof(ne,jj);
          // minus sign otherwise the force is upside down
          ForceMarteau(l) -= ForceMarteauExt(j)*base_time.EvaluatePhi(jj, x);
        }       
      }
      
      name_file = var.DOSSIER + "ForceMarteauInterp.txt";
      remove(name_file.data());
      ForceMarteau.WriteText(name_file);
    }
  }
  
  
  //! Delta function used to distribute the hammer on the string
  Real_wp Hammer::Fct_Delta(Real_wp x)
  {
    // spatial repartition of the hammer on the string
    return ((1.0/(1.0+exp(-pente_marteau*(x+epaisseur_marteau*0.5))) 
       - 1.0/(1.0+exp(-pente_marteau*(x-epaisseur_marteau*0.5))))/epaisseur_marteau);
  }
  
  
  //! Function phi(d) = Positive part(d - u)^p where p is the exponent
  //! delta = \bar\xi dans la thèse de Juliette Chabassier
  Real_wp Hammer::Phi(Real_wp u)
  {
    if((delta-u)>0)
    {
      //DISP(" Phi ne confirme PAS que Hammer parti");
      return pow(delta-u, exposant);
    }
    else
    {
      //DISP(" Phi confirme, Hammer parti");
      return 0.0;  
    }
  }
  
  
  //! Derivative of phi
  Real_wp Hammer::PhiPrime(Real_wp u)
  {
    if((delta-u) > 0)
    {
      //cout << " PhiPrime ne confirme PAS que Hammer parti" << endl;
      return -exposant*pow(delta-u, exposant-1);
    }
    else
    {
      //cout << " PhiPrime confirme, Hammer parti" << endl;
      return 0.0;
    }      
  }
  
  
  //! Function Psi (primitive of function Phi)
  Real_wp Hammer::Psi(Real_wp u)
  {
    if((delta-u) > 0)
    {
      return pow(delta-u,exposant+1.0)/(exposant+1.0);
    }
    else
    {
      //cout << "Hammer parti" << endl;
      return 0.0;
    }
  }
  
  
  //! Evaluation of the numerical derivate of function Psi
  Real_wp Hammer::DiffFiniePsi(Real_wp u, Real_wp v)
  {
    Real_wp deno=(u-v);
    if(abs(deno)>epsilon)
    {
      Real_wp tmp1 = Psi(u);
      Real_wp tmp2 = Psi(v);
      return((tmp1-tmp2)/deno);
    }
    else
    {
      return(-Phi(u));
    }
  }
  
  
  //! Evaluation of the numerical derivate of function Psi    
  Real_wp Hammer::DiffFiniePsiOPT(Real_wp u, Real_wp v)
  {
    Real_wp deno = (u-v);   
    Real_wp tmp1 = Psi(u);
    Real_wp tmp2 = Psi(v);
    Real_wp QU = (tmp1-tmp2);
    // il faut que Xi=abs(deno), puis DE=-LoiPhi(exposant,delta,u); , QU=QU/Xi;
    // Alors sous epsilon/2 c'est DE, au dessus de epsilon c'est QU
    // et entre les deux c'est 2*(1-Xi/epsilon)*DE + 2(Xi/epsilon-1/2)*QU
    if (abs(deno) > 1)
    {
      return QU/deno;
    }
    else
    {
      Real_wp DE = (-Phi(u));
      Real_wp signe = 1.0;
      if (deno < 0)
        signe = -1.0;
  
      abort(); //bugue
      return ((1-abs(deno))*DE + signe* QU); // (1-Xi)f1 + Xi f2
    }
  } 
  
  
  //! Evaluation of the derivative of numerical derivate of function Psi
  Real_wp Hammer::DeriveeDiffFiniePsi(Real_wp u, Real_wp v)
  {
    Real_wp deno = (u-v);
    if(abs(deno) > epsilon)
    {
      Real_wp tmp1 = Psi(u);
      Real_wp tmp2 = Psi(v);
      Real_wp tmp3 = Phi(u);
      return (-tmp3/deno - (tmp1-tmp2)/(deno*deno));
    }
    else
    {
      return (-PhiPrime(u)*0.5);
    }
  }
  
  
  //! Evaluation of the derivative of numerical derivate of function Psi
  Real_wp Hammer::DeriveeDiffFiniePsiOPT(Real_wp u, Real_wp v)
  {
    Real_wp deno = (u-v);
    Real_wp deno2 = deno*deno;
    Real_wp tmp1 = Psi(u);
    Real_wp tmp2 = Psi(v);
    Real_wp tmp3 = Phi(u);
    Real_wp QU = -tmp3*deno - (tmp1-tmp2);
    abort(); //bugue
    if(abs(deno) > 1)
    {
      return QU/deno2;
    }
    else
    {
      Real_wp DE = -PhiPrime(u)*0.5;
      return ((1-deno2)*DE + QU);
    }
  }


  //! Function phi_1(d) = K * Positive part(d - u)^p where p is the exponent
  Real_wp Hammer::Phi1Aux(Real_wp u, int i)
  {
    if (regularised_contact_model)
    {
      if (contact_model_with_string == ORIGINAL)
      {
        return -DeriveeGAux(u, i) * GAux(u, i);
      }
      else if (contact_model_with_string == ROBOTRAN)
      {
        abort();
      }
    }
    else
    {
      if (contact_model_with_string == ORIGINAL)
      {
        if((delta-u) > 0)
        {
          return Ki(i) * pow(delta-u, exposant);
        }
        else
        {
          return 0.0;  
        }
      }
      else if (contact_model_with_string == ROBOTRAN)
      {
        if((delta-u) > 0)
        {
          return Ai(i) * pow(delta-u, 4) + Bi(i) * pow(delta-u, 3) + Ci(i) * pow(delta-u, 2) + Di(i) * (delta-u);
        }
        else
        {
          return 0.0;  
        }
      }
    }
    return 0; // Juste pour éviter un warning
  }


  //! Function phi_2(d) = R * p * Positive part(d - u)^(p-1) where p is the exponent
  Real_wp Hammer::Phi2Aux(Real_wp u, int i)
  {
    if (regularised_contact_model)
    {
      if (contact_model_with_string == ORIGINAL)
      {
        return Ri(i) * exposant * pow((delta-u)*Echelon(delta-u), exposant-1.0);
      }
      else if (contact_model_with_string == ROBOTRAN)
      {
        abort();
      }
    }
    else
    {
      if (contact_model_with_string == ORIGINAL)
      {
        if((delta-u) > 0)
        {
          return Ri(i) * exposant * pow(delta-u, exposant-1);
        }
        else
        {
          return 0.0; 
        }
      }
      else if (contact_model_with_string == ROBOTRAN)
      {
        if((delta-u) > 0)
        {
          return 0.8 * Phi1Aux(u, i);   // A modifier coef D !!
        }
        else
        {
          return 0.0;
        }
      }
    }
    return 0; // Juste pour éviter un warning
  }


  //! Function GAux(d) = sqrt(2*Psi1Aux) (auxiliary variable function)
  Real_wp Hammer::GAux(Real_wp u, int i)
  {
    if (regularised_contact_model)
    {
      int n = 10000;
      if (delta-u < 0)
          return 0.0;
      else
      {
        Real_wp res = 0;
        for (int j = 1; j < n; j++)
        {
          res += DeriveeGAux(u + j*(delta-u)/n, i);
        }
        return (delta-u)/n * ((DeriveeGAux(u, i)+DeriveeGAux(delta, i))/2 + res);
      }
    }
    else
    {
      if (contact_model_with_string == ORIGINAL)
      {
        if((delta-u) > 0)
        {
          return sqrt(2.0*Ki(i)/(exposant + 1)) * pow(delta-u, 0.5*(exposant+1.0));
        }
        else
        {
          return 0.0;
        }
      }
      else if (contact_model_with_string == ROBOTRAN)
      {
        if((delta-u) > 0)
        {
          return sqrt(2) * sqrt(Ai(i) * pow(delta-u, 5)/5 + Bi(i) * pow(delta-u, 4)/4 + Ci(i) * pow(delta-u, 3)/3 + Di(i) * pow(delta-u, 2)/2);
        }
        else
        {
          return 0.0;
        }
      }
    }
    return 0; // Juste pour éviter un warning
  }


  //! Derivative Function of GAux(d)
  Real_wp Hammer::DeriveeGAux(Real_wp u, int i)
  {
    if (regularised_contact_model)
    {
      return -0.5*sqrt(2.0*Ki(i)*(exposant + 1)) * pow((delta-u)*Echelon(delta-u), 0.5*(exposant-1.0));
    }
    else
    {
      if (contact_model_with_string == ORIGINAL)
      {
        if((delta-u) > 0)
        {
          //DISP("CONTACT avec la corde");
          return -0.5*sqrt(2.0*Ki(i)*(exposant + 1)) * pow(delta-u, 0.5*(exposant-1.0));
        }
        else
        {
          //DISP("PAS DE CONTACT avec la corde");
          return 0.0;
        }
      }
      else if (contact_model_with_string == ROBOTRAN)
      {
        if((delta-u) > 0)
        {
          return -Phi1Aux(u, i) / GAux(u, i);
        }
        else
        {
          return 0.0;
        }
      }
    }
    return 0; // Juste pour éviter un warning
  }


  //! Fonction échelon régularisé C2
  Real_wp Hammer::Echelon(Real_wp x)
  {
    Real_wp dx = 1e-4;

    if(x < 0)
    {
      return 0;
    }
    else if (x >= dx)
    {
      return 1;
    }
    else
    {
      return (x/dx)*(x/dx)*(x/dx)*(10 - 15*(x/dx) + 6*(x/dx)*(x/dx));
    }
  }


}

#define MONTJOIE_FILE_HAMMER_CXX
#endif
