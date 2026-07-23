#ifndef MONTJOIE_FILE_SHANK_CXX

namespace Montjoie
{
  
  //! default constructor
  template<int nb_base, int nb_quad, class TypeEquation>
  Shank<nb_base, nb_quad, TypeEquation>::Shank(Hammer& h) : hammer(h)
  {
    is_interacting_string = false;
    E = 1.0;
    rho = 1.0;
    A = 1.0;
    I = 1.0;
    G = 1.0;
    k_prime = 0.85;
    T0 = 0.0;
    L = 1.0;
    gravity = 9.81;
    
    t_begin = 0.0;
    t_end = 1.0;
    t_courant = 0.0;
    t_display = 1.0;
    Deltat = 0.0;
    Deltat2 = 0.0;
    invDeltat = 0.0;
    invDeltat2 = 0.0;
    invDeuxDeltat = 0.0;
    
    alpha = 0.25;
    DOSSIER = "./";
    
    t_begin_interp = 0.0;
    t_end_interp = 0.0;
    deltat_interp = 0.0;
    
    t_begin_sismo = 0.0;
    t_end_sismo = 0.0;
    deltat_sismo = 0.0;
    
    point_impact_jack = 0.01;
    thickness_jack = 0.0;
    slope_jack = 2012;
    
    t_begin_jack = 0.0;
    t_end_jack = 0.0;
    DeltatForceJack = 0.0;
    amplitude_jack = 1.0;
    push_jack_until_theta = false;
    theta_end_jack = 0.0;
    type_condition_jack = JACK_THETA;
    distance_min_jack = 1.8e-3;
    type_function_jack = JACK_RECTANGLE;

    tangential_component = false;
    torque_component = false;

    amo_theta = 0.0;
    time_attrape = 1000000;    
    amo_theta_after_attrape = 0.0;
    
    initial_teta = 0.0;
    H_head = 0.0;
    mass_head = 0.0;
    Hprime_head = 0.0;
    d0 = 0.0;
    y_wall = 1e40;
    distance_corde_pivot = 0.0;
    
    last_d_stored = 0.0;
    last_Fhead_stored = 0.0;
    last_Fjack_stored = 0.0;
    last_tangential_Fjack_stored = 0.0;
    last_torque_jack_stored = 0.0;
    
    offset_shank = 0;
    formulation_with_multiplier = true;
    Fcoupl_dof_number = 0;
    xi_dof_number = 0;
    neglect_w_square = false;
    print_level = 0;
    
    size_sismo_buffer = 100;
    
    free_kinetic_energy = 0.0;
    free_potential_energy = 0.0;
    free_internal_energy = 0.0;
    free_shank_energy = 0.0;
    
    has_first_contact_happened = false;
    
    file_output_sismo = "Sismo_Shank.txt";
    file_output_force = "Sismo_Force.txt";
    file_output_jack = "Sismo_Fjack.txt";
  }
  
  
  //! returns the number of unknowns of the system
  template<int nb_base, int nb_quad, class TypeEquation>
  int Shank<nb_base, nb_quad, TypeEquation>::GetM() const
  {
    // N unknowns for w, phi, 1 for theta, 1 to enforce mean value of w equal to 0
    if (formulation_with_multiplier)
      return TypeEquation::dimension*this->mesh.GetNbDof() + 6;
    
    return TypeEquation::dimension*this->mesh.GetNbDof() + 2;
  }
  
  
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::ComputeRightHandSide(VectReal_wp & scheme, Matrix<Real_wp> & val2,
       bool newton_rhs)
  {
    cout << "not implemented" << endl;
    abort();
  }
  
  
  //! construction de tous les tableaux, matrices necessaires aux iterations en temps
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::ConstructAll(int num, bool double_prec, bool binary)
  {
    
    output_shank.SetDoublePrecision(double_prec);
    output_Force_d.SetDoublePrecision(double_prec);
    output_Fjack.SetDoublePrecision(double_prec);
    
    output_shank.SetBinary(binary);
    output_Force_d.SetBinary(binary);
    output_Fjack.SetBinary(binary);
    
    bool remove_file = true;
    output_shank.Init(DOSSIER +file_output_sismo,
          size_sismo_buffer, remove_file);
    
    output_Force_d.Init(DOSSIER +file_output_force,
      size_sismo_buffer, remove_file);
    output_Fjack.Init(DOSSIER+file_output_jack,
          size_sismo_buffer, remove_file);
    
    this->ConstructMesh(L);
    this->ConstruitInterp(L);
    
    lambda_dof_number = TypeEquation::dimension*this->mesh.GetNbDof();
    teta_dof_number = lambda_dof_number+1;
      
    lambda_dof_number += offset_shank;
    teta_dof_number += offset_shank;
    
    if (formulation_with_multiplier)
    {
      Fcoupl_dof_number =  teta_dof_number + 1;
      xi_dof_number = Fcoupl_dof_number + 2;
    }
    
    ConstructLinearMatrix();
    
    // on calcule la repartition de Fjack
    VectReal_wp repartition;
    this->ComputeRepartitionVector(repartition, repartition_sparse_jack, thickness_jack,
           point_impact_jack, *this);
    // En fait juste ici on triche on utilise la fonction dans SpaceStringScheme mais pas
    // comme il faudrait. On prend pas un genericPb hammer mais un genericPb shank, donc on
    // utilise la Fct_delta de shank
    
    repartition_sparse_jack.Get(0) = 0.0;
    
    // position des ddls
    PositionDofs.Reallocate(this->mesh.GetNbDof());
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
      for (int g = 0; g < nb_base; g++)
      {
        Real_wp x = (ne + this->base.Points(g))*this->Deltax;
        PositionDofs(this->mesh.GetNumberDof(ne, g)) = x;
      }
    
    // sismo sur le manche
    this->ConstruitObserv();
    if (this->elt_obs >= this->mesh.GetNbElt())
      this->elt_obs = this->mesh.GetNbElt()-1;
    
    // DISP(PositionDofs);
    /* Matrix<Real_wp> Fs(this->GetM(), 1);
       Fs.Fill(0); Real_wp sum = 0;
       for (int i = 0; i < repartition.GetM(); i++)
       {
       sum += repartition(i);
       Fs(2*i, 0) = repartition(i);
       }
       
       WriteInterp(0, Fs); */
  }
  
  
  //! calcul de la partie lineaire de la differentielle
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>::ConstructLinearMatrix()
  {
    TinyMatrix<Real_wp, General, dimension, nb_quad>    TMBase, dTMBase;
    TinyMatrix<Real_wp, General, dimension, nb_quad>    Mphi, Rphi, Adphi, AStiffdphi,
      Adamping_dphi, Bphi, Cphi, Ddphi;
    TinyMatrix<Real_wp, General, dimension, nb_quad>    vect_contre_phi, vect_contre_dphi;
    TinyMatrix<Real_wp, General, dimension, nb_base>    res_int_phi, res_int_dphi;
    
    int N = this->GetM();
    linear_matrix.Reallocate(N, N);
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
    {
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
          var_eq.ApplyM       (TMBase ,Mphi   ,*this);
          var_eq.ApplyR       (TMBase ,Rphi   ,*this);
          var_eq.ApplyAGrad     (dTMBase  ,Adphi    ,*this);
          var_eq.ApplyAStiffGrad  (dTMBase  ,AStiffdphi ,*this);
          var_eq.ApplyB       (TMBase ,Bphi   ,*this);
          //  var_eq.ApplyB     (TMBase,  Bphi,    *this);
          var_eq.ApplyC       (TMBase ,Cphi   ,*this);
          var_eq.ApplyD       (dTMBase  ,Ddphi    ,*this);
          var_eq.ApplyAdamping      (dTMBase  ,Adamping_dphi  ,*this);
          
          
          vect_contre_phi  = invDeltat2*Mphi + invDeuxDeltat* Rphi
                              + alpha * ( Cphi+Ddphi);// thetastiff * Bdphi +
          
          vect_contre_dphi = alpha* Adphi + alpha* ( AStiffdphi + Bphi) 
                              + invDeuxDeltat*Adamping_dphi;
          
          
          res_int_phi.Zero();   res_int_dphi.Zero();
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
              linear_matrix.AddInteraction(jglob1*dimension + k1,
                     jglob2*dimension+k2, value);
            } // for k1 (component number of the basis function)
          }// for j1 (local ddl basis function )
              
        }// for j2 (local ddl number of up)
      }// for k2 ( component number of up)
    }// for ne
    
    ValueI1_Base.Zero();
    MassMatrixBase.Zero();
    
    // Petite mass matrix
    for (int i = 0; i < nb_base; i++)
      for (int k = 0; k < nb_quad; k++)
      {
        ValueI1_Base(i) += this->Deltax*this->quad.Weights(k)*this->valPhi(i, k);
          
        for (int j = i; j < nb_base; j++)
          MassMatrixBase(i,j) += this->Deltax*this->quad.Weights(k)
            *this->valPhi(i,k)*this->valPhi(j,k);
      }
    
    // int x \varphi sur chaque element
    ValueIx_Mesh.Reallocate(this->mesh.GetNbElt());
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
    {
      for (int k = 0; k < nb_quad; k++)
      {
        Real_wp x = (ne + this->quad.Points(k))*this->Deltax;    
        for (int i = 0; i < nb_base; i++)
          ValueIx_Mesh(ne)(i) += this->Deltax*x*this->quad.Weights(k)*this->valPhi(i, k);
      }
    }
  }
  
  
  //! Les buffers de sortie sont vides
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>::FlushBuffer()    
  {
    output_shank.CloseBuffer();
    output_Fjack.CloseBuffer();
  }
  
  
  //! Le deplacement du manche est calcule sur les points d'interpolation et ecrit sur le disque
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>::WriteInterp(int n, const Matrix<Real_wp>& val)
  {
    int dim = TypeEquation::dimension;
    
    Matrix<Real_wp> tmp;
    tmp.Reallocate(this->Nx_interp, dim);
    tmp.Fill(Real_wp(0));
    
    for (int k = 0; k < dim; k++)
      for (int l = 0; l < this->Nx_interp; l++)
  {
    for (int p = 0; p < this->MatInterp.GetRowSize(l); p++)
      {
        int indice_ddl = this->MatInterp.Index(l, p);
        tmp(l, k) += val(offset_shank+indice_ddl*dimension+k, 0)
    *this->MatInterp.Value(l,p);
      }
  }
    
    string name_file = DOSSIER + "Val" + to_str(n) + "shank.txt";
    
    tmp.WriteText(name_file);
  }
  
  
  //! L'ecriture des instantanes est effectuee
  /*!
    \param[in] nb_iter numero d'iteration
    \param[in] t temps t^n
    \param[in] val solution u au temps t^{n-1}, t^n, t^{n+1}
  */
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::WriteSnapshots(int nb_iter, Real_wp t, const Matrix<Real_wp>& val)
  {
    //DISP(t_begin_interp); DISP(t_end_interp); DISP(deltat_interp);
    if (deltat_sismo > 0)
    {
      int n = toInteger(round(t/deltat_sismo));
      if ((t >= t_begin_sismo) && (t <= t_end_sismo)
              && (abs(t/deltat_sismo -n) < (Deltat*0.5)))
      {
        int last_dof_w = offset_shank + (this->mesh.GetNbDof()-1)*dimension;;
        Vector<Real_wp> temp(10); temp.Fill(0);
        // on ecrit t, w(L), phi(L) et theta
        temp(0) = t;
        temp(1) = val(last_dof_w, 1);
        temp(2) = val(last_dof_w+1, 1);
        temp(3) = val(teta_dof_number, 1);
              
        // un point du manche, xi de la tete de marteau, xi au bout du marteau
        Real_wp cos_teta = cos(temp(3)), sin_teta = sin(temp(3));
        Real_wp w = 0, x = this->point_observation, wL = temp(1);
        for(int p = 0; p < nb_base; p++)
        {
          int pp = this->mesh.GetNumberDof(this->elt_obs, p);
          w += val(offset_shank + dimension*pp, 1)*this->observ(p);
        }
      
        temp(4) = x*cos_teta + w*sin_teta;
        temp(5) = x*sin_teta - w*cos_teta;
        temp(6) = val(xi_dof_number, 1);
        temp(7) = val(xi_dof_number+1, 1);
        temp(8) = L*cos_teta + (wL -Hprime_head)*sin_teta;
        temp(9) = L*sin_teta - (wL -Hprime_head)*cos_teta;
              
        output_shank.AddVect(temp);
      }
    }
    
    if (deltat_interp != 0)
    {
      int n = toInteger(round(t / deltat_interp));
      //DISP(n);
      if ((t >= t_begin_interp) && (t <= t_end_interp)
              && (abs(t / deltat_interp - n) < (Deltat*0.5)))
      {     
        //cout << "j'appelle WriteInterp"<< endl;
        WriteInterp(n, val);            
      }
    }
    
    Real_wp d_plus = last_d_stored;
    Real_wp Fhead = last_Fhead_stored;
    TinyVector<Real_wp, 3> tmp;
    tmp(0) = t_courant;
    tmp(1) = d_plus;
    tmp(2) = Fhead;
    output_Force_d.AddTinyVect(tmp);
    
    TinyVector<Real_wp,2> tempJack;
    tempJack(0) = t_courant;
    tempJack(1) = last_Fjack_stored;
    output_Fjack.AddTinyVect(tempJack);

  }
  
  
  //! Changement du pas de temps
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>::SetTimeStep(const Real_wp& dt)
  {
    Deltat      = dt;
    Deltat2     = square(Deltat);
    invDeltat     = 1.0/Deltat;
    invDeltat2    = square(invDeltat);
    invDeuxDeltat   = 0.5*invDeltat;
  }
  
  
  //! returns true is the keyword is specific to the shank
  template<int nb_base, int nb_quad, class TypeEquation>
  bool Shank<nb_base, nb_quad, TypeEquation>::AllowedKeyword(const string& keyword)
  {
    if (keyword == "Shank")
      return true;
    
    if (keyword == "InitialConditionShank")
      return true;
    
    if (keyword == "ShankDamping")
      return true;
    
    if (keyword == "NeglectWsquare")
      return true;
    
    if (keyword == "ThetaDampingAfterTime")
      return true;
    
    if (keyword == "ShankInterpTimeGrid")
      return true;
    
    if (keyword == "ShankSismoTimeGrid")
      return true;
    
    if (keyword == "ShankInterp")
      return true;

    if (keyword == "ForceJack")
      return true;
    
    if (keyword == "ForceTimeJack")
      return true;

    if (keyword == "TangentialForceTimeJack")
      return true;

    if (keyword == "TorqueTimeJack")
      return true;
    
    if (keyword == "HammerHead")
      return true;
    
    if (keyword == "FormulationMultiplier")
      return true;
    
    if (keyword == "FileOutputShank")
      return true;
    
    return false;
  }
  
  
  //! setting shank parameters with a line of the data file
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::SetInputData(const string & keyword, const Vector<string>& parameters)
  {   
    if (!keyword.compare("Shank"))
    {
      if (parameters.GetM() <= 1)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "Shank needs at least two parameters, for instance :" << endl;
        cout << "Shank = AUTO N" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      if (!parameters(0).compare("YES"))
      { 
        if (parameters.GetM() <= 9)
        {
          cout << "In SetInputData of Shank" << endl;
          cout << "Shank needs nine parameters, for instance :" << endl;
          cout << "Shank = YES L E rho A I G Kappa N d" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
      
        L = to_num<Real_wp>(parameters(1));
        E = to_num<Real_wp>(parameters(2));
        rho = to_num<Real_wp>(parameters(3));
        A = to_num<Real_wp>(parameters(4));
        I = to_num<Real_wp>(parameters(5));
        G = to_num<Real_wp>(parameters(6));
        k_prime = to_num<Real_wp>(parameters(7));
        this->Nx = to_num<int>(parameters(8));
        distance_corde_pivot = to_num<Real_wp>(parameters(9));
      }
      else
      {
        // parametres du bouleau
        E = 16.0e9;
        rho = 610.0;
        G = 1e9;
            
        // parametres dans Forss p 246
        Real_wp d = 3.9e-3;
        A = pi_wp*d*d/4.0;
        I = pi_wp*pow(d, 4.0)/64.0;
        L = 0.13;
            
        this->Nx = to_num<int>(parameters(1));
      }
    }
    else if (!keyword.compare("PrintLevel"))
    {
      if (parameters.GetM() <= 0)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "PrintLevel needs one parameter, for instance :" << endl;
        cout << "PrintLevel = 2" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      print_level = to_num<int>(parameters(0));
    }
    else if (!keyword.compare("InitialConditionShank"))
    {
      if (parameters.GetM() <= 0)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "InitialConditionShank needs one parameter, for instance :" << endl;
        cout << "InitialConditionShank = teta0" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      initial_teta = to_num<Real_wp>(parameters(0))*pi_wp/180.0;
    }
    else if (!keyword.compare("Hammer"))
    {
      if (!is_interacting_string)
        hammer.SetInputData(keyword, parameters);
    }
    else if (!keyword.compare("HammerStiffness"))
    {
      if (parameters.GetM() <= 1)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "HammerStiffness needs two parameters, for instance :" << endl;
        cout << "HammerStiffness = K R" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      hammer.Ki.Reallocate(1);
      hammer.Ri.Reallocate(1);
      hammer.Ki(0) = to_num<Real_wp>(parameters(0));
      hammer.Ri(0) = to_num<Real_wp>(parameters(1));
    }
    else if (!keyword.compare("ShankDamping"))
    {
      if (parameters.GetM() <= 6)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "ShankDamping needs seven parameters, for instance :" << endl;
        cout << "ShankDamping = amo0 amo1 amo2 amoB2_0 amoB2_1 amoB2_2 amo_teta" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      amo(0) = to_num<Real_wp>(parameters(0));
      amo(1) = to_num<Real_wp>(parameters(1));
      amo(2) = to_num<Real_wp>(parameters(2));
            
      amoB2(0) = to_num<Real_wp>(parameters(3));
      amoB2(1) = to_num<Real_wp>(parameters(4));
      amoB2(2) = to_num<Real_wp>(parameters(5));
            
      amo_theta = to_num<Real_wp>(parameters(6));
    }
    else if (!keyword.compare("NeglectWsquare"))
    {
      if (parameters.GetM() <= 0)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "NeglectWsquare needs one parameter, for instance :" << endl;
        cout << "NeglectWsquare = YES" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      if (parameters(0)=="YES")
        neglect_w_square = true;
      else
        neglect_w_square = false;
    }
    else if (!keyword.compare("ThetaDampingAfterTime"))
    {
      if (parameters.GetM() <= 1)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "ThetaDampingAfterTime needs two parameters, for instance :" << endl;
        cout << "ThetaDampingAfterTime = t0 amo_teta" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      time_attrape = to_num<Real_wp>(parameters(0));
      amo_theta_after_attrape = to_num<Real_wp>(parameters(1));
    }
    else if (!keyword.compare("TimeStep"))
    {
      if (parameters.GetM() <= 0)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "TimeStep needs one parameter, for instance :" << endl;
        cout << "TimeStep = dt" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      SetTimeStep(to_num<Real_wp>(parameters(0)));
    } 
    else if (!keyword.compare("NewtonParam"))
    {
      Newton_solver.SetInputData(keyword,parameters);
    }
    else if (!keyword.compare("ForceReevaluationJacobian"))
    {
      Newton_solver.SetInputData(keyword,parameters);
    }
    else if (!keyword.compare("TimeInterval"))
    {
      if (parameters.GetM() <= 1)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "TimeInterval needs two parameters, for instance :" << endl;
        cout << "TimeInterval = t0 tf" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      t_begin = to_num<Real_wp>(parameters(0));
      t_end   = to_num<Real_wp>(parameters(1));
    }
    else if (!keyword.compare("DisplayFrequency"))
    {
      if (parameters.GetM() <= 0)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "DisplayFrequency needs one parameter, for instance :" << endl;
        cout << "DisplayFrequency = dt" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      t_display     = to_num<Real_wp>(parameters(0));
    }      
    else if (!keyword.compare("TimeSchemeShank"))
    {
      if (parameters.GetM() <= 1)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "TimeSchemeShank needs two parameters, for instance :" << endl;
        cout << "TimeSchemeShank = THETA_SCHEME 0.25" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      alpha = to_num<Real_wp>(parameters(1));
    }
    else if (!keyword.compare("ShankInterpTimeGrid"))
    {
      if (parameters.GetM() <= 2)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "ShankInterpTimeGrid needs three parameters, for instance :" << endl;
        cout << "ShankInterpTimeGrid = t0 tf dt" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      t_begin_interp  = to_num<Real_wp>(parameters(0));
      t_end_interp    = to_num<Real_wp>(parameters(1));
      deltat_interp   = to_num<Real_wp>(parameters(2));
            
      if (deltat_interp < Deltat)
        cout << "Warning : your time step for interpolation cannot be respected." <<endl;
  
    }
    else if (!keyword.compare("ShankSismoTimeGrid"))
    {
      if (parameters.GetM() <= 3)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "ShankSismoTimeGrid needs four parameters, for instance :" << endl;
        cout << "ShankSismoTimeGrid = t0 tf dt size_buffer" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      t_begin_sismo   = to_num<Real_wp>(parameters(0));
      t_end_sismo   = to_num<Real_wp>(parameters(1));
      deltat_sismo    = to_num<Real_wp>(parameters(2));
      size_sismo_buffer = to_num<int>(parameters(3));
            
      if (deltat_sismo < Deltat)
        cout << "Warning : your time step for sismo cannot be respected." <<endl;        
  
      if (parameters.GetM() >= 5)
      {
        this->point_observation = to_num<Real_wp>(parameters(4));
        if (this->point_observation > L)
        {
          cout << "The observation point of the shank is too far" << endl;
          abort();
        }
      }
      else
        this->point_observation = 0.85*L;
    }
    else if (!keyword.compare("ShankInterp"))
    {
      if (parameters.GetM() <= 0)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "ShankInterp needs one parameter, for instance :" << endl;
        cout << "ShankInterp = N" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      this->Nx_interp = to_num<int>(parameters(0));   
    }
    else if (!keyword.compare("ForceJack"))
    {
      if (parameters.GetM() <= 2)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "ForceJack needs three parameters, for instance :" << endl;
        cout << "ForceJack = point_impact thickness slope" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      point_impact_jack = to_num<Real_wp>(parameters(0));
      thickness_jack = to_num<Real_wp>(parameters(1));
      slope_jack = to_num<Real_wp>(parameters(2));
    }
    else if (!keyword.compare("ForceTimeJack"))
    {
      if (parameters.GetM() <= 2)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "ForceTimeJack needs three parameters, for instance :" << endl;
        cout << "ForceTimeJack = THETA_END theta_end amplitude" << endl;
        cout << "ForceTimeJack = t0 tf amplitude" << endl;
        cout << "ForceTimeJack = FROMFILE file DeltatForce" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
      // Pour utiliser une force de jack d'un fichier
      if (parameters(0) == "FROMFILE") 
      {
        type_function_jack = FROM_FILE;
        string name_force_file(parameters(1));
        ForceJackExt.ReadText(name_force_file);
        DeltatForceJack = to_num<Real_wp>(parameters(2)); // Pas de temps utilisé pour générer le fichier de force
      }
      // Pour parametrer une force qui devient nulle
      // lorsque le marteau atteint un angle theta_end d'échappement
      else if (parameters(0) == "THETA_END")
      {
        type_condition_jack = JACK_THETA;
        theta_end_jack = to_num<Real_wp>(parameters(1))*pi_wp/180.0; // angle du marteau auquel se produit l'échappement
        amplitude_jack = to_num<Real_wp>(parameters(2));
        push_jack_until_theta = true;
        t_begin_jack = 0.0; t_end_jack = 10.0;
        if (parameters.GetM() > 3)
        {
          type_function_jack = JACK_RECTANGLE;
          if (parameters(3) == "COS_SQUARE")
            type_function_jack = JACK_COS_SQUARE;
          else if (parameters(3) == "X_SQRT_COS")
            type_function_jack = JACK_X_SQRT_COS;
          
          t_begin_jack = 0.0;
          t_end_jack = to_num<Real_wp>(parameters(4));
        }
      }
      // Pour parametrer une force qui devient nulle
      // lorsque la distance marteau/corde atteint une
      // distance dist_min d'échappement 
      else if (parameters(0) == "DIST_MIN")
      {
        type_condition_jack = JACK_DISTANCE;
        distance_min_jack = to_num<Real_wp>(parameters(1)); // = distance marteau/corde à laquelle se produit l'échappement
        amplitude_jack = to_num<Real_wp>(parameters(2));
        push_jack_until_theta = true;
        t_begin_jack = 0.0; t_end_jack = 10.0;
        if (parameters.GetM() > 3)
        {
          type_function_jack = JACK_RECTANGLE;
          if (parameters(3) == "COS_SQUARE")
            type_function_jack = JACK_COS_SQUARE;
          else if (parameters(3) == "X_SQRT_COS")
            type_function_jack = JACK_X_SQRT_COS;
          
          t_begin_jack = 0.0;
          t_end_jack = to_num<Real_wp>(parameters(4));
        }
      }
      else
      {
        t_begin_jack = to_num<Real_wp>(parameters(0));
        t_end_jack = to_num<Real_wp>(parameters(1));
        amplitude_jack = to_num<Real_wp>(parameters(2));
        if (parameters.GetM() > 3)
        {
          type_function_jack = JACK_RECTANGLE;
          if (parameters(3) == "COS_SQUARE")
            type_function_jack = JACK_COS_SQUARE;
          else if (parameters(3) == "X_SQRT_COS")
            type_function_jack = JACK_X_SQRT_COS;
        }
      }
    }
    else if (!keyword.compare("TangentialForceTimeJack"))
    {
      if (parameters.GetM() <= 2)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "TangentialForceTimeJack needs three parameters, for instance :" << endl;
        cout << "TangentialForceTimeJack = THETA_END theta_end amplitude" << endl;
        cout << "TangentialForceTimeJack = t0 tf amplitude" << endl;
        cout << "TangentialForceTimeJack = FROMFILE file DeltatForce" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
      // Pour utiliser une force de jack d'un fichier
      if (parameters(0) == "FROMFILE") 
      {
        tangential_component = true;
        type_function_jack_tangential = FROM_FILE;
        string name_force_file(parameters(1));
        TangentialForceJackExt.ReadText(name_force_file);
        DeltatTangentialForceJack = to_num<Real_wp>(parameters(2)); // Pas de temps utilisé pour générer le fichier de force
      }
      else
      {
        abort();  // A coder pour prendre en compte les autres types
      }
    }
    else if (!keyword.compare("TorqueTimeJack"))
    {
      if (parameters.GetM() <= 2)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "TorqueTimeJack needs three parameters, for instance :" << endl;
        cout << "TorqueTimeJack = THETA_END theta_end amplitude" << endl;
        cout << "TorqueTimeJack = t0 tf amplitude" << endl;
        cout << "TorqueTimeJack = FROMFILE file DeltatTorque" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
      // Pour utiliser une force de jack d'un fichier
      if (parameters(0) == "FROMFILE") 
      {
        torque_component = true;
        type_function_jack_torque = FROM_FILE;
        string name_force_file(parameters(1));
        TorqueJackExt.ReadText(name_force_file);
        DeltatTorqueJack = to_num<Real_wp>(parameters(2)); // Pas de temps utilisé pour générer le fichier de force
      }
      else
      {
        abort();  // A coder pour prendre en compte les autres types
      }
    }
    else if (!keyword.compare("HammerHead"))
    {
      if (parameters.GetM() <= 1)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "HammerHead needs at least two parameters, for instance :" << endl;
        cout << "HammerHead = mass H" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      mass_head = to_num<Real_wp>(parameters(0));
      H_head = to_num<Real_wp>(parameters(1));
      if (parameters.GetM() > 2)
      {
        d0 = to_num<Real_wp>(parameters(2));
        if (parameters.GetM() > 3)
          y_wall = to_num<Real_wp>(parameters(3));
        else
          y_wall = 1e100;
        
        Hprime_head = H_head + d0;
      }
    }
    else if (!keyword.compare("FormulationMultiplier"))
    {
      if (parameters.GetM() <= 0)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "FormulationMultiplier needs one parameter, for instance :" << endl;
        cout << "FormulationMultiplier = YES" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      if (parameters(0) == "YES")
        formulation_with_multiplier = true;
      else
        formulation_with_multiplier = false;
    }
    else if (keyword == "FileOutputShank")
    {
      if (parameters.GetM() <= 2)
      {
        cout << "In SetInputData of Shank" << endl;
        cout << "FileOutputShank needs three parameters, for instance :" << endl;
        cout << "FileOutputShank = file_sismo file_force file_jack" << endl;
        cout << "Current parameters are : " << endl << parameters << endl;
        abort();
      }
  
      file_output_sismo = parameters(0);
      file_output_force = parameters(1);
      file_output_jack = parameters(2);
    }
  }
  
  
  //! initial conditions for theta, phi, w
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>::Init(int nb_iterations, Matrix<Real_wp>& val)
  {
    if (!is_interacting_string)
      val.Fill(0);
    
    if (formulation_with_multiplier)
    {
      val(xi_dof_number, 0) = L*cos(initial_teta) - H_head*sin(initial_teta);
      val(xi_dof_number+1, 0) = L*sin(initial_teta) + H_head*cos(initial_teta);
      
      val(xi_dof_number, 1) = val(xi_dof_number, 0);
      val(xi_dof_number+1, 1) = val(xi_dof_number+1, 0);
      val(xi_dof_number, 2) = val(xi_dof_number, 0);
      val(xi_dof_number+1, 2) = val(xi_dof_number+1, 0);
    }
    
    val(teta_dof_number, 0) = initial_teta;
    val(teta_dof_number, 1) = initial_teta;
    val(teta_dof_number, 2) = initial_teta;
    
    //val.FillRand(); Mlt(1e-9, val);
    
    if (!is_interacting_string)
    {
      Newton_solver.type_matrix_lu = Newton_solver.LU_ARROW;
      Newton_solver.size_band_lu = 2*nb_base;
      Newton_solver.size_last_row_lu = 6;
            
      Newton_solver.Init(*this, val);
    }
  }
  
  
  //! computation of w^n+1, phi^n+1, theta^n+1 from previous values
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::Advance(int nt, const Real_wp& t, Matrix<Real_wp>& val)
  {
    Newton_solver.Solve(*this, val);
  }
  
  
  //! displaying the parameters of the shank
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::Display(ostream& file_out, bool write_on_file)
  {
    cout << " -------- Shank Parameters: -------- "<<endl;
    DISP(E);
    DISP(rho);
    DISP(A);
    DISP(I);
    DISP(G);
    DISP(L);
    DISP(amo);
    DISP(amoB2);
    DISP(amo_theta);
    DISP(t_begin_jack);
    DISP(t_end_jack);
    DISP(amplitude_jack);
    DISP(point_impact_jack);
    DISP(initial_teta);
    DISP(H_head);
    DISP(d0);
    cout << " ----------------------------------- "<<endl;
      
    if (write_on_file)
      {
  file_out << " -------- Shank Parameters: -------- "<<endl;
  file_out << "E = " << E << endl;
  file_out << "rho = " << rho << endl;
  file_out << "A = " << A << endl;
  file_out << "I = " << I << endl;
  file_out << "G = " << G << endl;
  file_out << "L = " << L << endl;
  file_out << "amo = " << amo << endl;
  file_out << "amoB2 = " << amoB2 << endl;
  file_out << "amo_theta = " << amo_theta << endl;
  file_out << "t_begin_jack = " << t_begin_jack << endl;
  file_out << "t_end_jack = " << t_end_jack << endl;
  file_out << "amplitude_jack = " << amplitude_jack << endl;
  file_out << "point_impact_jack = " << point_impact_jack << endl; 
  file_out << "initial_teta = " << initial_teta << endl;
  file_out << "H_head = " << H_head << endl;
  file_out << "d0 = " << d0 << endl;
  file_out << " ----------------------------------- "<<endl << endl;
      }
  }
  
  
  //! returns the value of theta at instant n
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp Shank<nb_base, nb_quad, TypeEquation>
  ::GetThetaShank(const Matrix<Real_wp>& val, int n)
  {
    return val(teta_dof_number, n);
  }    
  
  
  //! sets the scheme for the row corresponding to the equation governing theta
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::SetThetaShank(Vector<Real_wp>& scheme, const Real_wp& x)
  {
    scheme(teta_dof_number) = x;
  }    
  
  
  //! d/dx (sin(x)/x)
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp Shank<nb_base, nb_quad, TypeEquation>::deriv_sinc(const Real_wp& x)
  {
    // DL : -x/3 + x^3/30
    if (abs(x) < 1e-5)
      return -x/Real_wp(3);
    
    return (x*cos(x) - sin(x))/(x*x);
  }
  
  
  //! (cos(a) - cos(b)) / (a-b)
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp Shank<nb_base, nb_quad, TypeEquation>::delta_cos(const Real_wp& a, const Real_wp& b)
  {
    Real_wp sin_cardinal = 1.0;
    if (abs(a-b) > 2e-9)
      sin_cardinal = sin(0.5*(a-b))/(0.5*(a-b));
    
    Real_wp result = -sin(0.5*(a+b))*sin_cardinal;
    return result;
  }
  
  
  //! d/da (cos(a) - cos(b)) / (a-b)
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp Shank<nb_base, nb_quad, TypeEquation>::deriv_delta_cos(const Real_wp& a, const Real_wp& b)
  {
    Real_wp sin_cardinal = 1.0;
    if (abs(a-b) > 2e-9)
      sin_cardinal = sin(0.5*(a-b))/(0.5*(a-b));
    
    Real_wp deriv_sin_cardinal = deriv_sinc(0.5*(a-b));
    Real_wp result = -0.5*cos(0.5*(a+b))*sin_cardinal -0.5*sin(0.5*(a+b))*deriv_sin_cardinal;
    return result;
  } 
  
  
  //! (sin(a) - sin(b)) / (a-b)
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp Shank<nb_base, nb_quad, TypeEquation>::delta_sin(const Real_wp& a, const Real_wp& b)
  {
    Real_wp sin_cardinal = 1.0;
    if (abs(a-b) > 2e-9)
      sin_cardinal = sin(0.5*(a-b))/(0.5*(a-b));
    
    Real_wp result = cos(0.5*(a+b))*sin_cardinal;
    return result;
  } 
  
  
  //! d/da (sin(a) - sin(b)) / (a-b)
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp Shank<nb_base, nb_quad, TypeEquation>::deriv_delta_sin(const Real_wp& a, const Real_wp& b)
  {
    Real_wp sin_cardinal = 1.0;
    if (abs(a-b) > 2e-9)
      sin_cardinal = sin(0.5*(a-b))/(0.5*(a-b));
    
    Real_wp deriv_sin_cardinal = deriv_sinc(0.5*(a-b));
    Real_wp result = -0.5*sin(0.5*(a+b))*sin_cardinal + 0.5*cos(0.5*(a+b))*deriv_sin_cardinal;
    return result;
  } 
  
  
  //! computing values of w, phi on quadrature points of element ne and values of gradient
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::ComputeTM(int ne, int instant, const Matrix<Real_wp> & val,
        TinyMatrix<Real_wp, General, dimension, nb_quad> & TM,
        TinyMatrix<Real_wp, General, dimension, nb_quad> & dTM)
  {
    TinyVector<Real_wp, nb_quad> quad_values_ligne;
    TinyVector<Real_wp, nb_quad> quad_gradient_ligne;
    // loop over unknowns w, phi
    for (int k = 0; k < dimension; k++) 
    {
      // we get the nodal values of this component
      TinyVector<Real_wp, nb_base> nodal_values;
      for (int j = 0; j < nb_base; j++) 
      {
        int jj = this->mesh.GetNumberDof(ne, j);
        int num = offset_shank + jj*dimension + k;
        nodal_values(j) = val(num, instant);
      }
  
      // we compute the values at the quad points
      MltTrans(this->valPhi, nodal_values, quad_values_ligne);
      MltTrans(this->valdPhi, nodal_values, quad_gradient_ligne);
            
      // we store it in TM
      SetRow(quad_values_ligne, k, TM);
      SetRow(this->invDeltax*quad_gradient_ligne, k, dTM);
    }
  }
  
  
  //! computing values of w, phi on quadrature points of element ne
  //! and values on degrees of freedom
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::ComputeWbase(int ne, int instant, const Matrix<Real_wp> & val,
     TinyVector<Real_wp, nb_base> & ubase,
     TinyVector<Real_wp, nb_quad> & TM)
  {
    // we get the nodal values of w
    for (int j = 0; j < nb_base; j++) 
    {
      int jj = this->mesh.GetNumberDof(ne, j);
      int num = offset_shank + jj*dimension;
      ubase(j) = val(num, instant);
    }
    
    // we compute the values at the quad points
    MltTrans(this->valPhi, ubase, TM);
  }
  
  
  //! computation of F(val) (we try to solve F(U^n+1, U^n, U^n-1) = 0 where U^n+1 is the unknown)
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::ComputeScheme(VectReal_wp & scheme, Matrix<Real_wp> & val)
  {
    scheme.Fill(0);
    AddScheme(scheme, val, 0, 0.0);
  }
  
  
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::ComputeScheme(VectReal_wp & scheme, VectReal_wp& rhs, Matrix<Real_wp> & val)
  {
    ComputeScheme(scheme, val);
  }
  
  
  //! on rajoute F(x) a scheme
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::AddScheme(VectReal_wp & scheme, Matrix<Real_wp> & val,
        int hammer_dof, const Real_wp& interaction_hammer, 
        const string& file_name)
  {
    TinyMatrix<Real_wp, General, dimension, nb_quad> TMplus, TM, TMminus, dTMplus, dTM, dTMminus;
    TinyMatrix<Real_wp, General, dimension, nb_quad> dtt, delta, dx_delta, dthetavect, thetavect;
    
    TinyMatrix<Real_wp, General, dimension, nb_quad> Mdtt, Rdelta, Adthetavect,
      AdthetaStiffvect, Bthetavect,
      Ddthetavect, Cthetavect, nablaU, Adamping_deltat;
    
    TinyMatrix<Real_wp, General, dimension, nb_quad> vect_contre_phi, vect_contre_dphi;
    TinyMatrix<Real_wp, General, dimension, nb_base> res_int_phi, res_int_dphi;
    
    Real_wp teta_plus, teta, teta_minus, dteta_dt, dteta_dt2, teta_pondere;
    
    teta_plus = this->GetThetaShank(val, 2);
    teta = this->GetThetaShank(val, 1);
    teta_minus = this->GetThetaShank(val, 0);
    
    dteta_dt2 = (teta_plus - 2.0*teta + teta_minus)*invDeltat2;
    dteta_dt = (teta_plus - teta_minus)*invDeuxDeltat;
    teta_pondere = this->alpha*teta_plus  + (1.0-2.0*this->alpha)*teta  + this->alpha*teta_minus;
    
    Real_wp lambda_plus, lambda, lambda_minus, lambda_pondere, mean_value_w = 0.0;
    //DISP(this->alpha); DISP(lambda_dof_number);
    lambda_plus = val(lambda_dof_number, 2);
    lambda = val(lambda_dof_number, 1);
    lambda_minus = val(lambda_dof_number, 0);      
    lambda_pondere = this->alpha*lambda_plus  + (1.0-2.0*this->alpha)*lambda 
      + this->alpha*lambda_minus;
    //DISP(lambda_pondere);
    
    Real_wp force_grav_cos = rho*A*gravity*(this->alpha*cos(teta_plus)
              + (1.0-2.0*this->alpha)*cos(teta)
              + this->alpha*cos(teta_minus));
    
    Real_wp force_grav_dsin = rho*A*gravity*delta_sin(teta_plus, teta_minus);
    Real_wp force_grav_dcos = rho*A*gravity*delta_cos(teta_plus, teta_minus);
    
    Real_wp coef_tetapp = 0.0, coef_teta_plus(0), coef_teta_minus(0), coef_const = 0.0, poids;
    
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++) 
    {
      // We fill the TM with values at quad points.
      ComputeTM (ne, 2, val, TMplus, dTMplus);
      ComputeTM (ne, 1, val, TM, dTM);
      ComputeTM (ne, 0, val, TMminus, dTMminus);
      
      // We prep the tinymatrices
      dtt = (TMplus-Real_wp(2.0)*TM+TMminus)*invDeltat2;
      delta = (TMplus-TMminus)*invDeuxDeltat;
            
      dx_delta = (dTMplus-dTMminus)*invDeuxDeltat;
      thetavect = this->alpha*TMplus  + (1.0-2.0*this->alpha)*TM  + this->alpha*TMminus;
      dthetavect = this->alpha*dTMplus + (1.0-2.0*this->alpha)*dTM + this->alpha*dTMminus;
      
      // We apply the linear terms
      this->var_eq.ApplyM (dtt, Mdtt, *this);
      this->var_eq.ApplyR (delta, Rdelta, *this);
      this->var_eq.ApplyAGrad (dthetavect, Adthetavect, *this);
      this->var_eq.ApplyAStiffGrad (dthetavect, AdthetaStiffvect ,*this);
      this->var_eq.ApplyB (thetavect, Bthetavect, *this);
      this->var_eq.ApplyC (thetavect, Cthetavect, *this);
      this->var_eq.ApplyD (dthetavect, Ddthetavect, *this);
      this->var_eq.ApplyAdamping (dx_delta, Adamping_deltat, *this);
            
      // We prep vectors for numerical integration
      vect_contre_phi = Mdtt + Rdelta + Cthetavect+ Ddthetavect;
      vect_contre_dphi = Adthetavect + AdthetaStiffvect  + Bthetavect + Adamping_deltat;
            
      for (int g = 0; g < nb_quad; g++)
      {
        // part + lambda \int phi
        vect_contre_phi(0, g) += lambda_pondere;
              
        // part -rho A w (d theta/dt)^2
        if (!neglect_w_square)
        {
          //vect_contre_phi(0, g) -= rho*A*TM(0, g)*(teta_plus-teta)*(teta-teta_minus)
          //  *invDeltat2;
                      
          vect_contre_phi(0, g) -= rho*A*(Real_wp(0.25)*TMplus(0, g) +Real_wp(0.5)*TM(0, g)
                  + Real_wp(0.25)*TMminus(0, g))*(teta_plus-teta)*(teta-teta_minus)*invDeltat2;
        }
      
        // part -rho A d^2 theta / dt^2 x
        Real_wp x = (ne + this->quad.Points(g))*this->Deltax;
        vect_contre_phi(0, g) -= rho*A*dteta_dt2*x;
              
        // part -rho A g cos theta
        vect_contre_phi(0, g) -= force_grav_cos;
              
        // part -rho I d^2 theta / dt^2
        vect_contre_phi(1, g) -= rho*I*dteta_dt2;
              
        // part for theta_pp : rho A (w^2 + x^2) + rho I 
        poids = this->Deltax*this->quad.Weights(g);
        //DISP(x); DISP(poids);
        
        // part for teta_pp 
        coef_tetapp += (rho*A*x*x + rho*I)*poids;
              
        // part for theta_p : rho A w wp and rho A w wm
        if (!neglect_w_square)
        {
          //coef_teta_minus += rho*A*TM(0, g)*TMminus(0, g)*poids;
          //coef_teta_plus += rho*A*TM(0, g)*TMplus(0, g)*poids;
                      
          coef_teta_minus += rho*A*square(Real_wp(0.5)*(TMminus(0, g) + TM(0, g)))*poids;
          coef_teta_plus += rho*A*square(Real_wp(0.5)*(TMplus(0, g) + TM(0, g)))*poids;
        }
      
        // constant part : -rho A x wpp - rho I phi_pp + rho A g (x cos teta + w sin teta)
        coef_const += (-rho*A*x*dtt(0, g) - rho*I*dtt(1, g)
           + x*force_grav_dsin - thetavect(0, g)*force_grav_dcos)*poids;
        //coef_const += (x*force_grav_cos - rho*I*dtt(1, g)
        //  + thetavect(0, g)*force_grav_sin)*poids;
        //+ thetavect(0, g)*poids;
        
        // updating mean value of w
        mean_value_w += thetavect(0, g)*poids;
      }
  
      res_int_phi.Zero();
      res_int_dphi.Zero();
      this->CalculeIntPhi(vect_contre_phi, res_int_phi);
      this->CalculeIntGrad(vect_contre_dphi ,res_int_dphi);
      
      
      // We copy that in the scheme.
      for (int j = 0; j < nb_base; j++) 
      {
        int jglob = this->mesh.GetNumberDof(ne, j);
        for (int k = 0; k < dimension; k++)
        {
          int num = offset_shank + jglob*dimension + k;
          Real_wp value = res_int_phi(k,j) + res_int_dphi(k,j);
          scheme(num) += value;
        }
      }
    }
    
    scheme(lambda_dof_number) = mean_value_w;
    
    Real_wp amo_theta_true = amo_theta;
    if (t_courant > time_attrape)
      amo_theta_true = amo_theta_after_attrape;
    
    scheme(teta_dof_number) += dteta_dt2*coef_tetapp
      + dteta_dt*(rho*A*pow(L, 3.0)*amo_theta_true) + coef_const;
    
    if (!neglect_w_square)
      scheme(teta_dof_number) += (coef_teta_plus*(teta_plus-teta)
          - coef_teta_minus*(teta-teta_minus))*invDeltat2;
    
    // contribution due a M_H
    int last_dof_w = offset_shank + (this->mesh.GetNbDof()-1)*dimension;
    Real_wp w_pp = (val(last_dof_w, 2) - 2.0*val(last_dof_w, 1) + val(last_dof_w, 0))*invDeltat2;
    Real_wp w_p = (val(last_dof_w, 2) - val(last_dof_w, 0))*invDeuxDeltat;
    Real_wp w = this->alpha*(val(last_dof_w, 2) + val(last_dof_w, 0))
      + (1.0 - 2.0*this->alpha)*val(last_dof_w, 1);
    
    Real_wp wn = val(last_dof_w, 1), w_minus = val(last_dof_w, 0), w_plus = val(last_dof_w, 2);
    
    if (formulation_with_multiplier)
    {
      // Added terms for Fcoupl in shank and theta lines
      // -(Fcoupl.utheta) w*(L) (on the lhs of the w line of the shank)
      // thetaline rhs :  (w(L)-H) (Fcoupl . ur) -L(Fcoupl.utheta)
      Real_wp F_uteta = val(Fcoupl_dof_number, 2)*sin(teta)
        - val(Fcoupl_dof_number+1, 2)*cos(teta);
      
      Real_wp F_ur = val(Fcoupl_dof_number, 2)*cos(teta)
        + val(Fcoupl_dof_number+1, 2)*sin(teta);
      
      scheme(last_dof_w) -= F_uteta;
      scheme(teta_dof_number) -= (w - H_head)*F_ur - L*F_uteta;
      
      // mH d2xi_dt2 = -mH g uy - Fcoupl (written in the cartesian basis ux, uy) 
      Real_wp xi_x_pp = (val(xi_dof_number, 2) - 2.0*val(xi_dof_number, 1)
             + val(xi_dof_number, 0))*invDeltat2;
      
      Real_wp xi_y_pp = (val(xi_dof_number+1, 2) - 2.0*val(xi_dof_number+1, 1)
             + val(xi_dof_number+1, 0))*invDeltat2;
      
      scheme(Fcoupl_dof_number) += mass_head*xi_x_pp + val(Fcoupl_dof_number, 2);
      scheme(Fcoupl_dof_number+1) += mass_head*xi_y_pp + mass_head*gravity
        + val(Fcoupl_dof_number+1, 2);
      
      // We impose the continuity on the discrete time derivative 
      // dxi_dt = dtheta_dt (w(L)-H) ur + (dw(L)_dt - L dtheta_dt ) utheta
      Real_wp dxi_x_dt = (val(xi_dof_number, 2) - val(xi_dof_number, 0))*invDeuxDeltat;
      Real_wp dxi_y_dt = (val(xi_dof_number+1, 2) - val(xi_dof_number+1, 0))*invDeuxDeltat;
      scheme(xi_dof_number) += dxi_x_dt - dteta_dt*(w - H_head)*cos(teta) 
        - (w_p - L*dteta_dt)*sin(teta);
      
      scheme(xi_dof_number+1) += dxi_y_dt - dteta_dt*(w - H_head)*sin(teta)
        + (w_p - L*dteta_dt)*cos(teta);
    }
    else
    {
      force_grav_cos /= rho*A; 
      force_grav_dcos /= rho*A;
      force_grav_dsin /= rho*A;
      
      // part m_H ( \ddot{w(L)} - L \ddot{theta} - (w(L) - H) \dot{theta}^2 - g cos theta)
      scheme(last_dof_w) += mass_head*(w_pp - L*dteta_dt2 - force_grav_cos);
      if (!neglect_w_square)
        scheme(last_dof_w) -= mass_head*(w - H_head)
          *(teta_plus-teta)*(teta-teta_minus)*invDeltat2;
      
      // part m_H ( \ddot{theta} (w(L) - H)^2 + 2 \dot{theta} \dot{w(L)} (w(L) - H) 
      //     - L \ddot{w(L)} + L^2 \ddot{theta}  + L g cos theta + g (w(L) - H) sin theta )
      scheme(teta_dof_number) += mass_head*( - L*w_pp + L*L*dteta_dt2 + L*force_grav_dsin
                     - (w-H_head)*force_grav_dcos);
      if (!neglect_w_square)
        scheme(teta_dof_number) += mass_head*(square(Real_wp(0.5)*(wn + w_plus)-H_head)
                *(teta_plus-teta)
                - square(Real_wp(0.5)*(wn + w_minus)-H_head)
                *(teta-teta_minus))*invDeltat2;
    }
    
    // application de Fjack
    Real_wp Fjack = 0.0;
    Real_wp force_jack = 0.0;
    if (type_function_jack == FROM_FILE)
    {
      force_jack = NewForceJack(toInteger(round(t_courant*invDeltat)));
    }
    else if ((t_courant >= t_begin_jack) && (t_courant <= t_end_jack))
    {
      if (type_function_jack == JACK_COS_SQUARE)
      {
        Real_wp x = (t_courant - t_begin_jack) / (t_end_jack - t_begin_jack);
        force_jack = amplitude_jack*square(cos(0.5*pi_wp*(2.0*x-1.0)));
      }
      else if (type_function_jack == JACK_X_SQRT_COS)
      {
        Real_wp x = (t_courant - t_begin_jack) / (t_end_jack - t_begin_jack);
        force_jack = amplitude_jack*x*sqrt(cos(0.5*pi_wp*x));
      }
      else
      {
        force_jack = amplitude_jack;
      }
    }

    // application de Fjack composante tengentielle
    Real_wp Fjack_tangential = 0.0;
    Real_wp tangential_force_jack = 0.0;
    if (tangential_component)
    {
      if (type_function_jack_tangential == FROM_FILE)
      {
        tangential_force_jack = NewTangentialForceJack(toInteger(round(t_courant*invDeltat)));
      }
      else
      {
        abort();
      }
    }

    // application de Fjack composante couple pur
    Real_wp Fjack_torque = 0.0;
    Real_wp torque_jack = 0.0;
    if (torque_component)
    {
      if (type_function_jack_torque == FROM_FILE)
      {
        torque_jack = NewTorqueJack(toInteger(round(t_courant*invDeltat)));
      }
      else
      {
        abort();
      }
    }
    
    if (push_jack_until_theta)
    {
      if (type_condition_jack == JACK_THETA)
      {
        if (teta <= theta_end_jack)
          Fjack = force_jack;
        else
          push_jack_until_theta = false;
      }
      else if (type_condition_jack == JACK_DISTANCE)
      {
        if ( (val(xi_dof_number+1, 1) - H_head) < -distance_min_jack)
          Fjack = force_jack;
        else
          push_jack_until_theta = false;
      }
  
      if (!push_jack_until_theta)          
      {
        t_end_jack = t_courant;
        Real_wp vx = (val(xi_dof_number, 1) - val(xi_dof_number, 0))/Deltat;
        Real_wp vy = (val(xi_dof_number+1, 1) - val(xi_dof_number+1, 0))/Deltat;
        Real_wp velocity_let_off = sqrt(vx*vx + vy*vy);
        cout << "At t = " << t_courant << endl;
        cout << "Velocity at let-off position = " << velocity_let_off << endl;
              
        // printing the velocity in param.txt
        if (file_name.size() > 2)
        {
          ofstream file_out(file_name.data(), ios::app);
          file_out << "---------------------------------\n";
          file_out << "-------- Hammer simulation --------\n";
          file_out << "Velocity_let_off = " << velocity_let_off << '\n';
          file_out.close();
        }
      }
    }
    else if (type_function_jack == FROM_FILE)
      Fjack = force_jack;
    else if ((t_courant >= t_begin_jack) && (t_courant <= t_end_jack))
      Fjack = force_jack;

    if ((tangential_component) && (type_function_jack_tangential == FROM_FILE))
    {
      Fjack_tangential = tangential_force_jack;
      last_tangential_Fjack_stored = Fjack_tangential;
    }
    if ((torque_component) && (type_function_jack_torque == FROM_FILE))
    {
      Fjack_torque = torque_jack;
      last_torque_jack_stored = Fjack_torque;
    }
    
    last_Fjack_stored = Fjack;
    

    // Répartition des forces sur la surface de contact
    for (int i = 0; i < repartition_sparse_jack.GetM(); i++)
    {
      int dof = repartition_sparse_jack.Index(i);
      Real_wp value = repartition_sparse_jack.Value(i);
      scheme(offset_shank + dimension*dof) += Fjack * value;
      scheme(teta_dof_number) -= Fjack * value * PositionDofs(dof);
      if (tangential_component)
      {
        // +1 parce que c'est la deuxième equation dans l'ordre de rangement
        scheme(offset_shank + dimension*dof + 1) += Fjack_tangential * 0.003 * value;  //Fjack_tangential * (0.003+0.004684) * value
        scheme(teta_dof_number) -= Fjack_tangential * value * (TM(0, i) + 0.003);
      }
    }
    /*if (tangential_component)
    {
      scheme(teta_dof_number) -= 0.004684 * Fjack_tangential;
    }*/
    if (torque_component)
    {
      scheme(teta_dof_number) -= Fjack_torque;
    }
    
    
    if ((abs(interaction_hammer) > 0) && (!has_first_contact_happened))
    {
      has_first_contact_happened = true;
      Real_wp vx = (val(xi_dof_number, 1) - val(xi_dof_number, 0))/Deltat;
      Real_wp vy = (val(xi_dof_number+1, 1) - val(xi_dof_number+1, 0))/Deltat;
      Real_wp velocity_contact = sqrt(vx*vx + vy*vy);
      cout << "At t = " << t_courant << endl;
      cout << "Velocity before contact = " << velocity_contact << endl;
        
      // printing the velocity in param.txt
      if (file_name.size() > 2)
      {
        ofstream file_out(file_name.data(), ios::app);
        file_out << "Velocity_contact = " << velocity_contact << '\n';
        file_out.close();
      }
    }
    
    // application de Fhead
    Real_wp Fhead(0);
    if (is_interacting_string)
    {
      // on impose que xsi_hammer = d0 + L sin teta - (w(L) - H') cos teta - H'
      if (formulation_with_multiplier)
      {
        //scheme(hammer_dof) = val(hammer_dof, 2) - val(xi_dof_number+1, 2) + (H_head + d0);
        scheme(hammer_dof) = val(hammer_dof, 2) - val(xi_dof_number+1, 2) + distance_corde_pivot;
      }
      else
        scheme(hammer_dof) = val(hammer_dof, 2)
                              - ( hammer.delta + (sin(teta_plus)*L - (val(last_dof_w, 2)
                                    - Hprime_head)*cos(teta_plus)) - Hprime_head );
      
      Fhead = -interaction_hammer;
      last_Fhead_stored = Fhead;
    }
    else if (y_wall < 1e30)
    {
      // cas d'un mur
      Real_wp d_minus;
      if (formulation_with_multiplier)
      {
        Real_wp d_plus = y_wall - val(xi_dof_number+1, 2);
        d_minus = y_wall - val(xi_dof_number+1, 0);
        Fhead = hammer.Ki(0)*hammer.DiffFiniePsi(d_plus, d_minus)
          - hammer.Ri(0)*(hammer.Phi(d_plus)-hammer.Phi(d_minus))*invDeuxDeltat;
        
        Fhead *= -1.0;              
      }
      else
      {
        Real_wp d_plus = hammer.delta + y_wall
          - (sin(teta_plus)*L - (val(last_dof_w, 2) - Hprime_head)*cos(teta_plus));
        
        d_minus = hammer.delta + y_wall -
          (sin(teta_minus)*L - (val(last_dof_w, 0) - Hprime_head)*cos(teta_minus));
        
        Fhead = hammer.Ki(0)*hammer.DiffFiniePsi(d_plus, d_minus)
          - hammer.Ri(0)*(hammer.Phi(d_plus)-hammer.Phi(d_minus))*invDeuxDeltat;
        
        Fhead *= -1.0;
      }
  
      // on stocke d et F
      last_d_stored = d_minus;
      last_Fhead_stored = Fhead;
    }
    
    //DISP(cos(teta_minus)); DISP(Hprime_head*cos(teta_minus));
    //DISP(last_dof_w); DISP(teta_dof_number); DISP(Fhead);
    if (formulation_with_multiplier)
    {
      scheme(Fcoupl_dof_number+1) += Fhead; 
    }
    else
    {
      scheme(last_dof_w) -= Fhead*cos(teta_pondere);
      scheme(teta_dof_number) -= -sin(teta_pondere)*Fhead*(w-Hprime_head)
                                  - Fhead*cos(teta_pondere)*L;
    }
    
    // condition de Dirichlet pour x = 0 (pour w)
    scheme(offset_shank) = val(offset_shank, 2);
    
    /*if (t_courant > 0.9995e-4)
      {
      scheme.WriteText("scheme_after.dat");
      val.WriteText("value.dat");
      int test_input; cout << "we wait" << endl; cin >> test_input;
      }*/
  }
  
  
  //! forcing constraints to be satisfied by the new solution
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::ForceConstraint(Matrix<Real_wp>& x, int hammer_dof)
  {
    Real_wp teta = this->GetThetaShank(x, 1);
    Real_wp teta_plus = this->GetThetaShank(x, 2);
    Real_wp teta_minus = this->GetThetaShank(x, 0);
    Real_wp dteta_dt = (teta_plus - teta_minus)*invDeuxDeltat;
    
    int last_dof_w = offset_shank + (this->mesh.GetNbDof()-1)*dimension;
    Real_wp w = this->alpha*(x(last_dof_w, 2) + x(last_dof_w, 0))
      + (1.0 - 2.0*this->alpha)*x(last_dof_w, 1);
    
    Real_wp w_p = (x(last_dof_w, 2) - x(last_dof_w, 0))*invDeuxDeltat;
    
    if (formulation_with_multiplier)
    {
      Real_wp vx, vy;
      vx = dteta_dt*(w - H_head)*cos(teta) + (w_p - L*dteta_dt)*sin(teta);
      vy = dteta_dt*(w - H_head)*sin(teta) - (w_p - L*dteta_dt)*cos(teta);
      x(xi_dof_number, 2) = x(xi_dof_number, 0) + 2.0*Deltat*vx;
      x(xi_dof_number+1, 2) = x(xi_dof_number+1, 0) + 2.0*Deltat*vy;
    }
    
    if ((is_interacting_string) && (formulation_with_multiplier))
    {
      if (hammer_dof >= 0)
      {
        //x(hammer_dof, 2) = x(xi_dof_number+1, 2) - (H_head+d0);
        x(hammer_dof, 2) = x(xi_dof_number+1, 2) - distance_corde_pivot;
      }
    } 
  }
  
  
  //! computation of jacobian matrix related to the equation F(U^n+1, U^n, U^n-1) = 0
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp Shank<nb_base, nb_quad, TypeEquation>
  ::ComputeDiff(VirtualMatrix<Real_wp>& DiffMatrix, Matrix<Real_wp> & val)
  {
    DiffMatrix.Zero();
    AddDiff(DiffMatrix, val, 0, 0);
    return 0.0;
  }
    
  
  //! Adds differential matrix to DiffMatrix
  template<int nb_base, int nb_quad, class TypeEquation>
  void Shank<nb_base, nb_quad, TypeEquation>
  ::AddDiff(VirtualMatrix<Real_wp>& DiffMatrix, Matrix<Real_wp> & val,
      int hammer_dof, int force_hammer_dof)
  {
    for (int i = 0; i < linear_matrix.GetM(); i++)
    {
      int nb_interac = linear_matrix.GetRowSize(i);
      IVect col_interac(nb_interac);
      VectReal_wp val_interac(nb_interac);
      for (int j = 0; j < nb_interac; j++)
      {
        col_interac(j) = offset_shank + linear_matrix.Index(i, j);
        val_interac(j) = linear_matrix.Value(i, j);
      }
  
      DiffMatrix.AddInteractionRow(offset_shank + i, nb_interac, col_interac, val_interac);
    }
    
    // loop over elements of the mesh
    Real_wp teta_plus, teta, teta_minus, dteta_dt, teta_pondere;
    
    teta_plus = this->GetThetaShank(val, 2);
    teta = this->GetThetaShank(val, 1);
    teta_minus = this->GetThetaShank(val, 0);
    
    dteta_dt = (teta_plus - teta_minus)*invDeuxDeltat;
    teta_pondere = this->alpha*teta_plus  + (1.0-2.0*this->alpha)*teta  + this->alpha*teta_minus;
    
    Real_wp force_grav_sin = rho*A*gravity*sin(teta_plus);
    Real_wp force_grav_dcos = rho*A*gravity*delta_cos(teta_plus, teta_minus);
    Real_wp force_grav_ddcos = -rho*A*gravity*deriv_delta_cos(teta_plus, teta_minus);
    Real_wp force_grav_ddsin = rho*A*gravity*deriv_delta_sin(teta_plus, teta_minus);
    
    TinyVector<Real_wp, nb_base> wbase_plus, wbase, wbase_minus,
      walpha, M_walpha, wpoint, M_wpoint, M_wplus, wplus;
    
    TinyVector<Real_wp, nb_quad> wquad, wquad_plus, wquad_minus, w_point, w_alpha;
    
    // part -alpha rho A \dot{theta}^2 M to add in rows of W
    Real_wp coef_ww = 0.0;
    if (!neglect_w_square)
      coef_ww = -Real_wp(0.25)*rho*A*(teta_plus-teta)*(teta-teta_minus)*invDeltat2;
    
    // part -2 (theta^(n+1) - theta^(n-1)) / (2 Delta t^2) rho A M w^alpha
    Real_wp coef_wtetaM =0.0;
    if (!neglect_w_square)
      coef_wtetaM = -(teta-teta_minus)*invDeltat2*rho*A;
    
    // part -rho A Ix / dt^2
    Real_wp coef_wtetaI = -rho*A/this->Deltat2;
    
    // part rho A g alpha sin teta I1
    Real_wp coef_wtetaS = force_grav_sin*alpha;
    
    // part -rho I I1 / dt^2
    Real_wp coef_phiteta = -rho*I*invDeltat2;
    
    Real_wp coef_tetawM_plus = 0.0;
    if (!neglect_w_square)
      coef_tetawM_plus = rho*A*(teta_plus-teta)*invDeltat2;
    
    // part -rho A Ix / dt^2
    Real_wp coef_tetawIx = -rho*A*invDeltat2;
    
    // part +alpha rho A g sin teta I1
    Real_wp coef_tetawI1 = -alpha*force_grav_dcos;
    
    // part -rho I / dt^2 I1
    Real_wp coef_tetaphi = -rho*I*invDeltat2;
    
    // part 1/dt^2 rho I L - rho A g L^2/2 alpha sin teta + 1/dt^2 rho A L^3/3
    Real_wp coef_tetateta = rho*I*L*invDeltat2 + invDeltat2*rho*A*L*L*L/3.0
      + force_grav_ddsin*L*L/2;
    
    // integral of w^2, w wpoint and w
    Real_wp integral_wplus(0), integral_w(0), poids(0);
    
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++) 
    {
      // retrieving values of w on degrees of freedom and on quadrature points
      ComputeWbase(ne, 2, val, wbase_plus, wquad_plus);
      ComputeWbase(ne, 1, val, wbase, wquad);
      ComputeWbase(ne, 0, val, wbase_minus, wquad_minus);
            
      w_point = (wquad_plus - wquad_minus)*invDeuxDeltat;
      w_alpha = this->alpha*wquad_plus  + (1.0-2.0*this->alpha)*wquad
        + this->alpha*wquad_minus;
      
      wpoint = (wbase_plus - wbase_minus)*invDeuxDeltat;
      walpha = this->alpha*wbase_plus  + (1.0-2.0*this->alpha)*wbase + this->alpha*wbase_minus;
      wplus = Real_wp(0.5)*(wbase + wbase_plus);
            
      Mlt(MassMatrixBase, walpha, M_walpha);
      Mlt(MassMatrixBase, wplus, M_wplus);
      Mlt(MassMatrixBase, wpoint, M_wpoint);
            
      for (int i = 0; i < nb_quad; i++)
      {
        poids = this->Deltax*this->quad.Weights(i);
        integral_w += w_alpha(i)*poids;
        //integral_wminus += square(Real_wp(0.5)*(wquad_minus(i) + wquad(i)))*poids;
        integral_wplus += square(Real_wp(0.5)*(wquad_plus(i) + wquad(i)))*poids;
      }
  
      for (int i = 0; i < nb_base; i++)
      {
        int irow = this->mesh.GetNumberDof(ne, i);              
        for (int j = 0; j < nb_base; j++)
        {
          int jcol = this->mesh.GetNumberDof(ne, j);
          // interaction w <-> w
          DiffMatrix
            .AddInteraction(offset_shank+irow*dimension,
                offset_shank+jcol*dimension, coef_ww*MassMatrixBase(i, j));
        }
      
        // interaction w, theta
        DiffMatrix
          .AddInteraction(offset_shank+irow*dimension, teta_dof_number,
              coef_wtetaM*M_walpha(i) + coef_wtetaS*ValueI1_Base(i)
              + coef_wtetaI*ValueIx_Mesh(ne)(i));
        
        // interaction w <-> lambda
        DiffMatrix
          .AddInteraction(offset_shank+irow*dimension,
              lambda_dof_number, alpha*ValueI1_Base(i));
        
        // interaction phi <-> teta
        DiffMatrix
          .AddInteraction(offset_shank+irow*dimension+1, teta_dof_number,
              coef_phiteta*ValueI1_Base(i));
        
        // interaction lambda <-> w
        DiffMatrix
          .AddInteraction(lambda_dof_number,
              offset_shank+irow*dimension, alpha*ValueI1_Base(i));
        
        // interaction teta <-> w (equation sur theta 13.c dérivée par rapport à la variable w)
        DiffMatrix
          .AddInteraction(teta_dof_number, offset_shank+irow*dimension,
              coef_tetawM_plus*M_wplus(i)
              + coef_tetawI1*ValueI1_Base(i) 
              + coef_tetawIx*ValueIx_Mesh(ne)(i));
        
        // interaction teta <-> phi
        DiffMatrix
          .AddInteraction(teta_dof_number, offset_shank+irow*dimension+1,
              coef_tetaphi*ValueI1_Base(i));
      }
    }      
    
    Real_wp amo_theta_true = amo_theta;
    if (t_courant > time_attrape)
      amo_theta_true = amo_theta_after_attrape;
    
    // part 1/dt^2 rho I L - rho A g L^2/2 alpha sin teta + 1/dt^2 rho A L^3/3
    //   + rho A / dt^2 \int w^2  + rho A / dt \int w wpoint + rho A g cos teta alpha \int w
    Real_wp value = coef_tetateta + force_grav_ddcos*integral_w
      + rho*A*pow(L, 3.0)*amo_theta_true*invDeuxDeltat;
    
    if (!neglect_w_square)
      value += rho*A*invDeltat2*integral_wplus;
    
    DiffMatrix.AddInteraction(teta_dof_number, teta_dof_number, value);
    
    // contribution due a M_H
    int last_dof_w = offset_shank+(this->mesh.GetNbDof()-1)*dimension;
    //Real_wp w_pp = (val(last_dof_w, 2) - 2.0*val(last_dof_w, 1) 
    //                + val(last_dof_w, 0))*invDeltat2;
    Real_wp w = this->alpha*(val(last_dof_w, 2) + val(last_dof_w, 0))
      + (1.0 - 2.0*this->alpha)*val(last_dof_w, 1);
    
    Real_wp wn = val(last_dof_w, 1), w_plus = val(last_dof_w, 2);
    
    //DISP(last_dof_w); DISP(teta_dof_number); DISP(Fcoupl_dof_number); DISP(xi_dof_number);
    if (formulation_with_multiplier)
    {
      //Real_wp F_uteta = val(Fcoupl_dof_number, 2)*sin(teta)
      //           - val(Fcoupl_dof_number+1, 2)*cos(teta);
      Real_wp F_ur = val(Fcoupl_dof_number, 2)*cos(teta)
        + val(Fcoupl_dof_number+1, 2)*sin(teta);
      
      DiffMatrix.AddInteraction(last_dof_w, Fcoupl_dof_number, -sin(teta));
      DiffMatrix.AddInteraction(last_dof_w, Fcoupl_dof_number+1, cos(teta));
      
      DiffMatrix.AddInteraction(teta_dof_number, last_dof_w, -this->alpha*F_ur);
      DiffMatrix.AddInteraction(teta_dof_number, Fcoupl_dof_number,
              -(w - H_head)*cos(teta) + L*sin(teta));
      DiffMatrix.AddInteraction(teta_dof_number, Fcoupl_dof_number+1,
              -(w - H_head)*sin(teta) - L*cos(teta));
      
      DiffMatrix.AddInteraction(Fcoupl_dof_number, xi_dof_number, mass_head*invDeltat2);
      DiffMatrix.AddInteraction(Fcoupl_dof_number, Fcoupl_dof_number, 1.0);
      
      DiffMatrix.AddInteraction(Fcoupl_dof_number+1, xi_dof_number+1, mass_head*invDeltat2);
      DiffMatrix.AddInteraction(Fcoupl_dof_number+1, Fcoupl_dof_number+1, 1.0);
      
      DiffMatrix.AddInteraction(xi_dof_number, last_dof_w,
              -dteta_dt*this->alpha*cos(teta) - invDeuxDeltat*sin(teta));
      DiffMatrix.AddInteraction(xi_dof_number, teta_dof_number,
              -invDeuxDeltat*(w - H_head)*cos(teta)
              + L *invDeuxDeltat*sin(teta));
      DiffMatrix.AddInteraction(xi_dof_number, xi_dof_number, invDeuxDeltat);
      
      DiffMatrix.AddInteraction(xi_dof_number+1, last_dof_w,
              -dteta_dt*this->alpha*sin(teta) + invDeuxDeltat*cos(teta));
      DiffMatrix.AddInteraction(xi_dof_number+1, teta_dof_number,
              -invDeuxDeltat*(w - H_head)*sin(teta)
              - L *invDeuxDeltat*cos(teta));
      DiffMatrix.AddInteraction(xi_dof_number+1, xi_dof_number+1, invDeuxDeltat);
    }
    else
    {
      force_grav_dcos /= rho*A;
      force_grav_ddsin /= rho*A;
      force_grav_ddcos /= rho*A;
      
      // interaction w(L), w(L)
      value = mass_head*invDeltat2;
      if (!neglect_w_square)
        value -= mass_head*this->alpha*(teta_plus-teta)*(teta-teta_minus)*invDeltat2;
      
      DiffMatrix.AddInteraction(last_dof_w, last_dof_w, value);
      
      // interaction w(L), theta
      value = mass_head*(-L*invDeltat2 + gravity*this->alpha*sin(teta_plus));
      if (!neglect_w_square)
        value -= mass_head*(w-H_head)*(teta-teta_minus)*invDeltat2;
      
      DiffMatrix.AddInteraction(last_dof_w, teta_dof_number, value);
      
      // interaction theta, w(L)
      value = mass_head*( -L*invDeltat2 - this->alpha*force_grav_dcos);
      if (!neglect_w_square)
        value += mass_head*Real_wp(0.5)*(wn + w_plus-2.0*H_head)*(teta_plus-teta)*invDeltat2;
      
      DiffMatrix.AddInteraction(teta_dof_number, last_dof_w, value);
      
      // interaction theta, theta
      value = mass_head*(L*L*invDeltat2 + L*force_grav_ddsin + (w-H_head)*force_grav_ddcos);
      if (!neglect_w_square)
        value += mass_head*square(Real_wp(0.5)*(wn + w_plus)-H_head)*invDeltat2;
      
      DiffMatrix.AddInteraction(teta_dof_number, teta_dof_number, value);
    }
    
    // part due to Fhead      
    Real_wp Fhead(0), dFhead_dw(0), dFhead_dteta(0), dFhead_dxi(0);
    if (is_interacting_string)
    {
      // on impose que xsi_hammer = d0 + L sin teta - (w(L) - H') cos teta + H'
      if (formulation_with_multiplier)
      {
        DiffMatrix.AddInteraction(hammer_dof, hammer_dof, 1.0);
        DiffMatrix.AddInteraction(hammer_dof, xi_dof_number+1, -1.0);
      }
      else
      {
        DiffMatrix.AddInteraction(hammer_dof, hammer_dof, 1.0);
        DiffMatrix.AddInteraction(hammer_dof, last_dof_w, cos(teta_plus));
        DiffMatrix.AddInteraction(hammer_dof, teta_dof_number, -L*cos(teta_plus) 
                - (val(last_dof_w, 2) - Hprime_head)*sin(teta_plus));
      }
    }
    else if (y_wall < 1e30)
    {
      // cas d'un mur
      if (formulation_with_multiplier)
      {
        Real_wp d_plus = y_wall - val(xi_dof_number+1, 2);
        Real_wp d_minus = y_wall - val(xi_dof_number+1, 0);
        Fhead = hammer.Ki(0)*hammer.DiffFiniePsi(d_plus, d_minus)
          - hammer.Ri(0)*(hammer.Phi(d_plus)-hammer.Phi(d_minus))*invDeuxDeltat;
        
        Fhead *= -1.0;
              
        Real_wp dPsi_dplus = hammer.DeriveeDiffFiniePsi(d_plus, d_minus);
        Real_wp dPhi_dplus = hammer.PhiPrime(d_plus);
              
        dFhead_dxi = (hammer.Ki(0)*dPsi_dplus - hammer.Ri(0)*dPhi_dplus*invDeuxDeltat);
      }
      else
      {
        Real_wp d_plus = hammer.delta + y_wall
          - (sin(teta_plus)*L - (val(last_dof_w, 2) - Hprime_head)*cos(teta_plus));
        
        Real_wp d_plus_dteta = -cos(teta_plus)*L
          - (val(last_dof_w, 2)- Hprime_head)*sin(teta_plus);
        
        Real_wp d_plus_dw = cos(teta_plus);
        Real_wp d_minus = hammer.delta + y_wall
          - (sin(teta_minus)*L - (val(last_dof_w, 0) - Hprime_head)*cos(teta_minus));
        
        Fhead = hammer.Ki(0)*hammer.DiffFiniePsi(d_plus, d_minus)
          - hammer.Ri(0)*(hammer.Phi(d_plus)-hammer.Phi(d_minus))*invDeuxDeltat;
        
        Real_wp dPsi_dplus = hammer.DeriveeDiffFiniePsi(d_plus, d_minus);
        Real_wp dPhi_dplus = hammer.PhiPrime(d_plus);
              
        dFhead_dteta = d_plus_dteta*(hammer.Ki(0)*dPsi_dplus
             - hammer.Ri(0)*dPhi_dplus*invDeuxDeltat);
        
        dFhead_dw = d_plus_dw*(hammer.Ki(0)*dPsi_dplus
             - hammer.Ri(0)*dPhi_dplus*invDeuxDeltat);
        
        //DISP(last_dof_w); DISP(teta_dof_number);
        Fhead *= -1.0; dFhead_dteta *= -1.0; dFhead_dw *= -1.0;
      }
    }
    
    if (formulation_with_multiplier)
    {
      if (is_interacting_string)
      {
        DiffMatrix.AddInteraction(Fcoupl_dof_number+1, force_hammer_dof, -1.0);
      }
      else
        DiffMatrix.AddInteraction(Fcoupl_dof_number+1, xi_dof_number+1, dFhead_dxi);
    }
    else
    {
      DiffMatrix.AddInteraction(last_dof_w, last_dof_w, -dFhead_dw*cos(teta_pondere));
      DiffMatrix.AddInteraction(last_dof_w, teta_dof_number, -dFhead_dteta*cos(teta_pondere)
              + this->alpha*Fhead*sin(teta_pondere));
      
      DiffMatrix.AddInteraction(teta_dof_number, last_dof_w, sin(teta_pondere)
              *(dFhead_dw*(w-Hprime_head) +  Fhead*this->alpha)
              + dFhead_dw*cos(teta_pondere)*L);
      
      DiffMatrix.AddInteraction(teta_dof_number, teta_dof_number,
              (w-Hprime_head)*(sin(teta_pondere)*dFhead_dteta 
                   + cos(teta_pondere)*this->alpha*Fhead)
              + dFhead_dteta*cos(teta_pondere)*L
              - Fhead*sin(teta_pondere)*L*this->alpha);
    }
    
    // Dirichlet on first component of w
    DiffMatrix.ClearRow(offset_shank);
    DiffMatrix.AddInteraction(offset_shank, offset_shank, 1.0);
    
    return;
    
    /*Real_wp h = pow(epsilon_machine, Real_wp(1.0)/Real_wp(3.0));
    VectReal_wp scheme_minus(val.GetM()), scheme_plus(val.GetM());
    //      Real_wp a(0), b(0);
    Matrix<Real_wp> DF_num(val.GetM(), val.GetM());
    for(int p = 0; p < val.GetM(); p++)
      {
  //boucle sur les variables
  //variations et calcul du scheme
  Real_wp souvenir = val(p, 2);
  val(p, 2) += h;
  this->ComputeScheme(scheme_plus, val);
        
  val(p, 2) -= 2.0*h;
  this->ComputeScheme(scheme_minus, val);
        
  for(int q = 0; q < val.GetM(); q++)
    {
      Real_wp DF = (scheme_plus(q) - scheme_minus(q)) / (2.0*h);
      DF_num(q, p) = DF;
      if (abs(DF - DiffMatrix(q, p)) > 1e-30)
        {
    DISP(p); DISP(q); DISP(DF);
    DISP(DiffMatrix(q, p));
    DISP(DF-DiffMatrix(q, p));
    abort();
        }
    }
  
  val(p, 2) = souvenir;
      }//for p
    */
    //DiffMatrix.WriteText("df.dat");      
    //DF_num.WriteText("df_num.dat");
    
    cout << "computation successful" << endl;
    exit(0);
    
    
    return;
  }
  
  
  //! Delta function used to distribute the jack on the shank
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp Shank<nb_base, nb_quad, TypeEquation>::Fct_Delta(Real_wp x)
  {
    // gaussian for pente == -1.0
    if (slope_jack == -1.0)
    {
      if (abs(x/thickness_jack) > 1.0)
        return 0.0;

      Real_wp distrib = 1.4927/thickness_jack*exp(-7.0*square(x/thickness_jack));
      return distrib;
    }
    
    //  spatial repartition of the jack on the shank
    return ((1.0/(1.0+exp(-slope_jack*(x+thickness_jack*0.5)))
       -1.0/(1.0+exp(-slope_jack*(x-thickness_jack*0.5))))/thickness_jack);
  }


  //! computation of the force on all time steps in the case where
  //! the force is read from a file
  template<int nb_base, int nb_quad, class TypeEquation>
  template<class GenericPb>
  void Shank<nb_base, nb_quad, TypeEquation>::ComputeInterpolatedJackForce(GenericPb & var)
  {
    Real_wp t_begin = var.t_begin;
    Real_wp t_end = var.t_end;

    string name_file = var.DOSSIER + "ForceJackExt.txt";
    remove(name_file.data());
    ForceJackExt.WriteText(name_file);

    int nb_force = toInteger(ceil(abs(t_end - t_begin)/DeltatForceJack ));
    int taille = ForceJackExt.GetM();
    ForceJackExt.Resize(nb_force+1);
    for (int p = taille; p < nb_force+1; p++)
      ForceJackExt(p) = 0.0;
    /*
    string name_file = var.DOSSIER + "ForceJackExt.txt";
    remove(name_file.data());
    ForceJackExt.WriteText(name_file);
    */
    int nb_base_time = 2;
    mesh_time.CreateRegularMesh(t_begin, t_end, nb_force+1, 1);
    mesh_time.SetOrder(nb_base_time-1);
    mesh_time.NumberMesh(false);
    base_time.ConstructQuadrature(nb_base_time-1, Globatto<Real_wp>::QUADRATURE_LOBATTO); 
    
    grid_interp_time.Init(t_begin, t_end, var.nb_max_iter+1);
    // en utilisant le grid_interp de montjoie
    grid_interp_time.LocalizePoints(mesh_time);
    
    NewForceJack.Reallocate(var.nb_max_iter+1);
    NewForceJack.Fill(Real_wp(0));
    for (int l = 0; l < grid_interp_time.GetNbPointsGrid(); l++)
    {
      int ne = grid_interp_time.GetElementNumber(l);
      Real_wp x = grid_interp_time.GetLocalCoordinate(l);
        
      for (int jj = 0; jj < nb_base_time; jj++)
      {
        int j = mesh_time.GetNumberDof(ne, jj);
        NewForceJack(l) += ForceJackExt(j)*base_time.EvaluatePhi(jj, x);
      }       
    }
    
    name_file = var.DOSSIER + "ForceJackInterp.txt";
    remove(name_file.data());
    NewForceJack.WriteText(name_file);
  }


  //! computation of the force on all time steps in the case where
  //! the force is read from a file
  template<int nb_base, int nb_quad, class TypeEquation>
  template<class GenericPb>
  void Shank<nb_base, nb_quad, TypeEquation>::ComputeInterpolatedTangentialJackForce(GenericPb & var)
  {
    Real_wp t_begin = var.t_begin;
    Real_wp t_end = var.t_end;

    string name_file = var.DOSSIER + "TangentialForceJackExt.txt";
    remove(name_file.data());
    TangentialForceJackExt.WriteText(name_file);

    int nb_force = toInteger(ceil(abs(t_end - t_begin)/DeltatTangentialForceJack ));
    int taille = TangentialForceJackExt.GetM();
    TangentialForceJackExt.Resize(nb_force+1);
    for (int p = taille; p < nb_force+1; p++)
      TangentialForceJackExt(p) = 0.0;
    /*
    string name_file = var.DOSSIER + "ForceJackExt.txt";
    remove(name_file.data());
    ForceJackExt.WriteText(name_file);
    */
    int nb_base_time = 2;
    mesh_time.CreateRegularMesh(t_begin, t_end, nb_force+1, 1);
    mesh_time.SetOrder(nb_base_time-1);
    mesh_time.NumberMesh(false);
    base_time.ConstructQuadrature(nb_base_time-1, Globatto<Real_wp>::QUADRATURE_LOBATTO); 
    
    grid_interp_time.Init(t_begin, t_end, var.nb_max_iter+1);
    // en utilisant le grid_interp de montjoie
    grid_interp_time.LocalizePoints(mesh_time);
    
    NewTangentialForceJack.Reallocate(var.nb_max_iter+1);
    NewTangentialForceJack.Fill(Real_wp(0));
    for (int l = 0; l < grid_interp_time.GetNbPointsGrid(); l++)
    {
      int ne = grid_interp_time.GetElementNumber(l);
      Real_wp x = grid_interp_time.GetLocalCoordinate(l);
        
      for (int jj = 0; jj < nb_base_time; jj++)
      {
        int j = mesh_time.GetNumberDof(ne, jj);
        NewTangentialForceJack(l) += TangentialForceJackExt(j)*base_time.EvaluatePhi(jj, x);
      }       
    }
    
    name_file = var.DOSSIER + "TangentialForceJackInterp.txt";
    remove(name_file.data());
    NewTangentialForceJack.WriteText(name_file);
  }


  //! computation of the force on all time steps in the case where
  //! the force is read from a file
  template<int nb_base, int nb_quad, class TypeEquation>
  template<class GenericPb>
  void Shank<nb_base, nb_quad, TypeEquation>::ComputeInterpolatedTorqueJack(GenericPb & var)
  {
    Real_wp t_begin = var.t_begin;
    Real_wp t_end = var.t_end;

    string name_file = var.DOSSIER + "TorqueJackExt.txt";
    remove(name_file.data());
    TorqueJackExt.WriteText(name_file);

    int nb_force = toInteger(ceil(abs(t_end - t_begin)/DeltatTorqueJack ));
    int taille = TorqueJackExt.GetM();
    TorqueJackExt.Resize(nb_force+1);
    for (int p = taille; p < nb_force+1; p++)
      TorqueJackExt(p) = 0.0;
    /*
    string name_file = var.DOSSIER + "ForceJackExt.txt";
    remove(name_file.data());
    ForceJackExt.WriteText(name_file);
    */
    int nb_base_time = 2;
    mesh_time.CreateRegularMesh(t_begin, t_end, nb_force+1, 1);
    mesh_time.SetOrder(nb_base_time-1);
    mesh_time.NumberMesh(false);
    base_time.ConstructQuadrature(nb_base_time-1, Globatto<Real_wp>::QUADRATURE_LOBATTO); 
    
    grid_interp_time.Init(t_begin, t_end, var.nb_max_iter+1);
    // en utilisant le grid_interp de montjoie
    grid_interp_time.LocalizePoints(mesh_time);
    
    NewTorqueJack.Reallocate(var.nb_max_iter+1);
    NewTorqueJack.Fill(Real_wp(0));
    for (int l = 0; l < grid_interp_time.GetNbPointsGrid(); l++)
    {
      int ne = grid_interp_time.GetElementNumber(l);
      Real_wp x = grid_interp_time.GetLocalCoordinate(l);
        
      for (int jj = 0; jj < nb_base_time; jj++)
      {
        int j = mesh_time.GetNumberDof(ne, jj);
        NewTorqueJack(l) += TorqueJackExt(j)*base_time.EvaluatePhi(jj, x);
      }       
    }
    
    name_file = var.DOSSIER + "TorqueJackInterp.txt";
    remove(name_file.data());
    NewTorqueJack.WriteText(name_file);
  }
  
    
  //! computation of the energy of the shank
  template<int nb_base, int nb_quad, class TypeEquation>
  Real_wp Shank<nb_base, nb_quad, TypeEquation>
  ::GetEnergy(const Matrix<Real_wp>& val, Real_wp& kinetic_energy,
        Real_wp& potential_energy, Real_wp& internal_energy)
  {
    Real_wp res = 0.0;
    TinyMatrix<Real_wp, General, dimension, nb_quad> TMplus, TM, dTMplus, dTM;
    TinyMatrix<Real_wp, General, dimension, nb_quad> TimeDerivative, TimeAverage, dTimeAverage;
    
    Real_wp teta_plus, teta, teta_point, poids;
    
    teta_plus = this->GetThetaShank(val, 1);
    teta = this->GetThetaShank(val, 0);
    
    teta_point = (teta_plus - teta)*invDeltat;
    
    Real_wp force_grav_cos = 0.5*rho*A*(cos(teta_plus)+cos(teta))*gravity;
    // 1 is added to ensure positivity of the potential energy
    Real_wp force_grav_sin = rho*A*(1.0+0.5*(sin(teta_plus)+sin(teta)))*gravity;
    
    kinetic_energy = 0; potential_energy = 0; internal_energy = 0;
    for (int ne = 0; ne < this->mesh.GetNbElt(); ne++)
    {
      // We fill the TM with values at quad points.
      ComputeTM (ne, 1, val, TMplus, dTMplus);
      ComputeTM (ne, 0, val, TM, dTM);
      
      // We prep the tinymatrices
      TimeDerivative  = (TMplus - TM)* invDeltat;
      TimeAverage   = (TMplus + TM)*Real_wp(0.5);
      dTimeAverage    = (dTMplus + dTM)*Real_wp(0.5);
            
      // loop over quadrature points
      for (int g = 0; g < nb_quad; g++)
      {
        Real_wp x = (ne + this->quad.Points(g))*this->Deltax;
        Real_wp w = TimeAverage(0, g);
        Real_wp phi = TimeAverage(1, g);
        Real_wp w_x = dTimeAverage(0, g);
        Real_wp phi_x = dTimeAverage(1, g);
        Real_wp w_point = TimeDerivative(0, g);
        Real_wp phi_point = TimeDerivative(1, g);
        
        poids = this->Deltax*this->quad.Weights(g);
              
        kinetic_energy +=  0.5*(rho*A*square(w_point - x*teta_point)
              + rho*I*square(teta_point - phi_point))*poids;
        if (!neglect_w_square)
        {
          //kinetic_energy += 0.5*rho*A*square(teta_point)*TM(0, g)*TMplus(0, g)*poids;
          kinetic_energy += 0.5*rho*A*square(teta_point*w)*poids;
        }
      
        internal_energy += 0.5*(A*G*k_prime*square(w_x-phi) + E*I*square(phi_x))*poids;
        potential_energy += (-force_grav_cos*w + force_grav_sin*x)*poids;
      }
    }
    
    int last_dof_w = offset_shank + (this->mesh.GetNbDof()-1)*dimension;
    Real_wp w_point = (val(last_dof_w, 1) - val(last_dof_w, 0))*invDeltat;
    Real_wp w = 0.5*(val(last_dof_w, 1) + val(last_dof_w, 0));
    
    if (formulation_with_multiplier)
    {
      Real_wp xsi_x_plus = val(xi_dof_number, 1);
      Real_wp xsi_x = val(xi_dof_number, 0);
      Real_wp xsi_y_plus = val(xi_dof_number+1, 1);
      Real_wp xsi_y = val(xi_dof_number+1, 0);
      kinetic_energy += 0.5*mass_head*(square((xsi_x_plus - xsi_x)*invDeltat)
               + square((xsi_y_plus - xsi_y)*invDeltat));
      
      potential_energy += mass_head*gravity*(xsi_y_plus + xsi_y)*0.5;
    }
    else
    {
      kinetic_energy += 0.5*mass_head*square(w_point - L*teta_point);
      if (!neglect_w_square)
        kinetic_energy += 0.5*mass_head*square(teta_point*(w-H_head));
    
      force_grav_cos /= rho*A; force_grav_sin /= rho*A;
      potential_energy += mass_head*(-force_grav_cos*(w-H_head) + force_grav_sin*L);
    }
    
    res = kinetic_energy + internal_energy + potential_energy;
    return res;
  }
  
}

#define MONTJOIE_FILE_SHANK_CXX
#endif
