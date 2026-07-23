#ifndef MONTJOIE_FILE_BRIDGE_CXX

namespace Montjoie
{
  
  Bridge::Bridge()
  {
    nb_lagrange_mult = 2;
    scal_a = 0.0;
    scal_b = 0.0;
    scal_c = 0.0;
    scal_d = 0.0;
    coef_couplage_V = 0.0;
    coef_couplage_H = 0.0;
    coef_couplage_any = 0.0;    
    offset_bridge = -1;
    size_bridge = 0;
    phi_control = 0;
    height_of_bridge = 0;
    nev = 0;
    
    vect_inertia = 1;
    vect_damp = 0;
    vect_stiff = 1;
    
    extra_bridge = false;
    matchNLaMuToStringDOF = false;
  }
  

  //! modifies the object with a line of the data file
  void Bridge::SetInputData(const string & keyword, const Vector<string> &param)
  {
    if (keyword == "BridgeSolver")
    {
  	  if (param(0) == "NUMERICAL")
  	  {
  	    size_bridge = 1;
     
  	    vect_inertia = to_num<Real_wp>(param(1));
  	    vect_damp = to_num<Real_wp>(param(2));
  	    vect_stiff = to_num<Real_wp>(param(3));
  	    
  	    phi_control = to_num<Real_wp>(param(4));
  	    height_of_bridge = to_num<Real_wp>(param(5));
  	  }
    }
  }


  //! reads impedance directly on a file ?
  void Bridge::ReadBridgeImpedance(const string& DOSSIER)
  {
    //  vect_stiff.ReadText(stiff_file_name);
    //  nev = vect_stiff.GetSize();
	
    //  vect_damp.ReadText(damp_file_name);
    //  if (vect_damp.GetSize() != nev)
    //    {
    //      cout << "Probleme taille " << damp_file_name << endl;
    //      abort();
    //    }
    
    //  vect_inertia.ReadText(inertia_file_name);
    //  if(vect_inertia.GetSize() != nev)
    //    {
    //      cout << "Probleme taille " << inertia_file_name << endl;
    //      abort();
    //    }
  }
  

  //! constructs everything needed for ComputeScheme
  template<class TypeEquation>
  void Bridge::ConstructAll(MultiString_Base& multi,TypeEquation & var_eq, const string& DOSSIER)
  {
    ReadBridgeImpedance(DOSSIER);
    int dimension = multi.GetNbUnknownsPerString();
    
    
    //todo - change to match user's preference
    // to decide how many LM to use
        
    if(matchNLaMuToStringDOF)
    {
        
      nb_lagrange_mult = dimension;
	
      if (var_eq.type_equation == TypeEquation::WAVE_LIN_TIMO)
      {
        var_eq.only_one_LM = false;
        var_eq.Dirichlet_Bridge(1) = true;
      }
      else if (var_eq.type_equation == TypeEquation::WAVE_NL_TIMO)
      {
        var_eq.Dirichlet_Bridge(2) = true;  
      }
      else if (var_eq.type_equation == TypeEquation::WAVE_NL_2T_TIMO)
      {
        var_eq.Dirichlet_Bridge(3) = true;
        var_eq.Dirichlet_Bridge(4) = true;
      }
    } 
    if (size_bridge > 0)
    {
      nb_lagrange_mult = dimension;
	    if (dimension > 2)
	      var_eq.Dirichlet_Bridge(2) = true;
    }
    
    vect_couplage_plaque.Reallocate(nb_lagrange_mult);    
    vect_scal_a.Reallocate(nb_lagrange_mult);  
    vect_scal_b.Reallocate(nb_lagrange_mult);
    vect_scal_a.Fill(0.0);
    vect_scal_b.Fill(0.0);      
  }
  
  int Bridge::GetSize()
  {
    return size_bridge;
  }

  int Bridge::GetNbLagrangeMultipliers() const
  {
    return nb_lagrange_mult;
  }


  void Bridge::SetNbLagrangeMultipliers(int n)
  {
    nb_lagrange_mult = n;
  }
  
  
  //! updates the equations to satisfy for the bridge part 
  /*!
    Warning : scheme contains at the beginning the values as computed
    by ComputeRightHandSide
   */
  void Bridge::AddScheme(MultiString_Base& multi,
              const VectReal_wp& somme_corde_0, const VectReal_wp& somme_corde_2,
              Matrix<Real_wp>& val2, VectReal_wp& scheme)
  {
    int nb_corde = multi.GetNbStrings();
    int dimension = multi.GetNbUnknownsPerString();
    Real_wp couplage_plaque_theta = 0;
    
    vect_couplage_plaque.Fill(0.0);
    
    for (int g2 = 0; g2 < nb_lagrange_mult; g2++)
    {
      for (int g = 0; g < nb_corde; g++)
      {
        vect_couplage_plaque(g2) += multi.GetLaMu(g2, g, 2, val2);
      }
	    vect_couplage_plaque(g2) = vect_scal_a(g2) - vect_couplage_plaque(g2)*vect_scal_b(g2);
    }
    
    if(extra_bridge)
    {
      Real_wp BridgeHeight0 = multi.GetBridgeHeight(0);		
      Real_wp BridgeHeight1 = multi.GetBridgeHeight(1);		
      Real_wp BridgeHeight2 = multi.GetBridgeHeight(2);
    	
      for (int g = 0; g < nb_corde; g++)
	    { 
        if (multi.GetStringType() == 0) //0 = WAVE_LIN
        {
		      couplage_plaque_theta += BridgeHeight2 * multi.GetLaMu(0,g,2,val2);
        }
        else if (multi.GetStringType() == 2) //2 = WAVE_NL
        {
		      couplage_plaque_theta += (BridgeHeight2 * multi.GetLaMu(0,g,2,val2)
                                          -BridgeHeight0 * multi.GetLaMu(1,g,2,val2));
        }
        else if (multi.GetStringType() == 6) //6 = WAVE_NL_2T_TIMO
        {       
		      couplage_plaque_theta += (BridgeHeight0 * multi.GetLaMu(1,g,2,val2) 
                                          - BridgeHeight1 * multi.GetLaMu(0,g,2,val2));
        }
        else if (multi.GetStringType() == 7) //7 = WAVE_LIN_2T
        {     
		      couplage_plaque_theta += (BridgeHeight0 * multi.GetLaMu(1,g,2,val2) 
                                          - BridgeHeight1 * multi.GetLaMu(0,g,2,val2));
        }		
        else
        {
          cout << "extra_bridge currently works for StringBridgeLIN or StringBridgeStiffNL2T" << endl;
          abort();
        }		
	    }
    }
    
    for (int g = 0; g < nb_corde; g++)
    {
	    // orthogonal LM part (nu part)
      if (multi.IsAttachedString(g))
        multi.SetSchemeLaMu(0, g, 0.0, scheme);
      else
        multi.SetSchemeLaMu(0, g, ((somme_corde_2(g)-somme_corde_0(g))*0.5 - vect_couplage_plaque(0)), scheme);
        
        
	    // tangential LM part (nu_ortho part)      
      for (int g2 = 1; g2 < nb_lagrange_mult; g2++)
      {
        Real_wp diffe2 = 0.0;
        Real_wp diffe0 = 0.0;            
        int last_dof_string = multi.GetNbStringDofs(g) - 1;
            
        for (int k2 = 0; k2 < dimension; k2 ++)
        {	
          Real_wp nu_k2 = multi.GetTransmissionChevaletString(g, g2, k2);
          diffe0 += multi.GetValue(g, k2, last_dof_string, 0, val2)*nu_k2;
          diffe2 += multi.GetValue(g, k2, last_dof_string, 2, val2)*nu_k2;
        }

        multi.SetSchemeLaMu(g2, g, (diffe2-diffe0)*0.5 - vect_couplage_plaque(g2), scheme);    
      }	
	
	
      // if extra bridge is present, needs to include the additional part
      if(extra_bridge)
      {
        Real_wp BridgeHeight2 = multi.GetBridgeHeight(2);
        Real_wp BridgeHeight1 = multi.GetBridgeHeight(1);			
        Real_wp BridgeHeight0 = multi.GetBridgeHeight(0);	
        if (multi.GetStringType() == 0)	
        {
          int dof_LM = multi.GetLaMuDofNumber(0,g);
          scheme(dof_LM) -=  BridgeHeight2*(scal_c - scal_d*couplage_plaque_theta);
        }
        else if (multi.GetStringType() == 2)
        {
          int dof_LM1 = multi.GetLaMuDofNumber(0,g);
          scheme(dof_LM1) -=  BridgeHeight2*(scal_c - scal_d*couplage_plaque_theta);
	
          int dof_LM2 = multi.GetLaMuDofNumber(1,g);
          scheme(dof_LM2) +=  BridgeHeight0*(scal_c - scal_d*couplage_plaque_theta);
        }		  
        else if (multi.GetStringType() == 6)
        { 
          int dof_LM1 = multi.GetLaMuDofNumber(1,g);
          scheme(dof_LM1) -=  BridgeHeight0*(scal_c - scal_d*couplage_plaque_theta);
	
          int dof_LM2 = multi.GetLaMuDofNumber(0,g);
          scheme(dof_LM2) +=  BridgeHeight1*(scal_c - scal_d*couplage_plaque_theta);
        }
        else if (multi.GetStringType() == 7) //7 = WAVE_LIN_2T
        { 
          int dof_LM1 = multi.GetLaMuDofNumber(1,g);
          scheme(dof_LM1) -=  BridgeHeight0*(scal_c - scal_d*couplage_plaque_theta);
	
          int dof_LM2 = multi.GetLaMuDofNumber(0,g);
          scheme(dof_LM2) +=  BridgeHeight1*(scal_c - scal_d*couplage_plaque_theta);
        }		  
        else
        {
          cout << "extra_bridge currently works for system with dimension = 1 or 5" << endl;
          cout << "Please do not include the keyword 'ExtraBridge' in your .ini file " << endl; 
          abort();
        }
      } 
      // other parts are unknowns
      if(size_bridge > 0)
      {
  	    //~ Real_wp theta_np1 = val2(offset_bridge,2);
  	    //~ Real_wp theta_nm1 = val2(offset_bridge,0);
  	    //~ int last_dof_string = multi.GetNbStringDofs(g) - 1;
              //~ 
  	    //~ for (int k = 2; k < nb_lagrange_mult; k++)
              //~ {
              //~ int dof_LM = multi.GetLagrangeM_DofNumber(g, k);
              //~ Real_wp diffe0 = multi.GetValue(g, k, last_dof_string, 0, val2);
              //~ Real_wp diffe2 = multi.GetValue(g, k, last_dof_string, 2, val2);
              //~ 
              //~ 
              //~ scheme(dof_LM) = (diffe2-diffe0)*0.5 - (theta_np1-theta_nm1)*0.5 ;
              //~ }
      }
      
      
    	// equations for the bridge
    	if (size_bridge > 0)
      {	
  	    // we retrieve values of theta at different time steps
  	    Real_wp theta_n = val2(offset_bridge,1);
  	    Real_wp theta_np1 = val2(offset_bridge,2);
  	    Real_wp theta_nm1 = val2(offset_bridge,0);  
          
  	    // part I theta_ddot
  	    scheme(offset_bridge) = (theta_np1 - 2.0*theta_n + theta_nm1)*multi.invDeltat2*vect_inertia;
              
  	    // part B dtheta/dt
  	    scheme(offset_bridge) += vect_damp*(theta_np1 - theta_nm1)*multi.invDeuxDeltat;
              
  	    // we retrieve values of the force (it is an unknown at time n)
  	    // 0 : force in nu direction, 1 : force in nu_ortho, 2 : etc
  	    // need more work in describing the forces with different components
  	    int dof_LM = multi.GetLagrangeM_DofNumber(g, 0);
  	    Real_wp Force_nu = val2(dof_LM, 2);
          Real_wp Force_5 = 0;
  	    if (nb_lagrange_mult > 2)
            Force_5 = val2(multi.GetLagrangeM_DofNumber(g,2),2);
          
  	    //~ cout << "Force_nu and force_5 are " << Force_nu << " and " << Force_5 << endl;
  	    // part fsdlkfjsdlkfjsd
  	    scheme(offset_bridge) += vect_stiff*(phi_control*(theta_np1 + theta_nm1) + (1-2*phi_control)*theta_n) + height_of_bridge*Force_nu + Force_5 ;
          
          //~ 
          //~ DISP(offset_bridge);
          //~ DISP(dof_LM);
          //~ cout << "index of val2 for Force_5 " << multi.GetLagrangeM_DofNumber(g,2) << endl;
          //~ getchar();
      }
    }
  }
  
  Real_wp Bridge::ComputeEnergy(MultiString_Base & multi, Matrix<Real_wp>& val2)
  {
    Real_wp ener_bridge = 0.0;
    if(size_bridge > 0)
    {
      Real_wp osci_rot_nm1 = val2(offset_bridge,0);
      Real_wp osci_rot_n = val2(offset_bridge,1);
      ener_bridge = 0.5*(vect_inertia+(phi_control-0.25)*multi.Deltat2*vect_stiff)*square(multi.invDeltat*(osci_rot_nm1-osci_rot_n)) ;
      ener_bridge+= 0.5*vect_stiff*square(0.5*(osci_rot_nm1+osci_rot_n));
    }
    return(ener_bridge);
    
  }
  void Bridge
  ::AddRightHandSide(MultiString_Base& multi,
                     const VectReal_wp& somme_corde_0,
                     Matrix<Real_wp>& val2, VectReal_wp& scheme)
  {
	  
    int nb_corde = multi.GetNbStrings();
    int dimension = multi.GetNbUnknownsPerString();
    Real_wp couplage_plaque_theta = 0;
    vect_couplage_plaque.Fill(0.0);
    
    for (int g2 = 0; g2 <nb_lagrange_mult; g2++)
    {	
      for (int g = 0; g < nb_corde; g++)
      {
        vect_couplage_plaque(g2) += multi.GetLaMu(g2, g, 2, val2);
      }
      vect_couplage_plaque(g2) = vect_scal_a(g2) - vect_couplage_plaque(g2)*vect_scal_b(g2);
    }
    
    if(extra_bridge)
    {
      Real_wp BridgeHeight2 = multi.GetBridgeHeight(2);
      Real_wp BridgeHeight1 = multi.GetBridgeHeight(1);	
	    Real_wp BridgeHeight0 = multi.GetBridgeHeight(0);	
        
      for (int g = 0; g < nb_corde; g++)
	    {  
        if (multi.GetStringType() == 0)
        {
		      couplage_plaque_theta += BridgeHeight2 * multi.GetLaMu(0,g,2,val2);
        }
        else if (multi.GetStringType() == 2) //2 = WAVE_NL
        {
		      couplage_plaque_theta += (BridgeHeight2 * multi.GetLaMu(0,g,2,val2)
                                          -BridgeHeight0 * multi.GetLaMu(1,g,2,val2));
        }
        else if (multi.GetStringType() == 6)
        {
		      couplage_plaque_theta += (BridgeHeight0 * multi.GetLaMu(1,g,2,val2) 
                                          - BridgeHeight1 * multi.GetLaMu(0,g,2,val2));
        }
        else if (multi.GetStringType() == 7) //7 = WAVE_LIN_2T
        {
		      couplage_plaque_theta += (BridgeHeight0 * multi.GetLaMu(1,g,2,val2) 
                                          - BridgeHeight1 * multi.GetLaMu(0,g,2,val2));
        }		
        else
        {
          cout << "extra_bridge currently works for StringBridgeLIN or StringBridgeStiffNL2T" << endl;
          abort();
        }		
	    }
    }
    
    for (int g = 0; g < nb_corde; g++)
    {
      if (multi.IsAttachedString(g))
        multi.SetSchemeLaMu(0, g, 0.0, scheme);
      else
        multi.SetSchemeLaMu(0, g, -somme_corde_0(g)*0.5 - vect_couplage_plaque(0), scheme);
        
      for (int g2 = 1; g2 < nb_lagrange_mult; g2++)
      {
        Real_wp diffe0 = 0.0;
        int last_dof_string = multi.GetNbStringDofs(g) - 1;
        for (int k2 = 0; k2 < dimension; k2 ++)
        {
          Real_wp	nu_k2 = multi.GetTransmissionChevaletString(g, g2, k2);  
          diffe0 += multi.GetValue(g, k2, last_dof_string, 0, val2)*nu_k2;
        }
            
        // New line for the force orthogonal to NU
        multi.SetSchemeLaMu(g2, g, -diffe0*0.5 - vect_couplage_plaque(g2), scheme);
      }	
        
      // if extra bridge is present, we need to include the additional terms
      // for multi.GetStringType() == 0
      //  extra term is : h * du/dx (L,t)
      // for multi.GetStringType() == 6
      // extra term is h_2 * du_1/dt (L,t) + h_1 * du_2/dt (L,t)
      // where h = [h_0,h_1,h+2], u = [u_0,u_1,u_2].. i.e. index starts from 0
      if(extra_bridge)
      {
        Real_wp BridgeHeight2 = multi.GetBridgeHeight(2);
        Real_wp BridgeHeight1 = multi.GetBridgeHeight(1);
        Real_wp BridgeHeight0 = multi.GetBridgeHeight(0);	
            
        if (multi.GetStringType() == 0)	
        {
          int dof_LM = multi.GetLaMuDofNumber(0,g);
          scheme(dof_LM) -=  BridgeHeight2*(scal_c - scal_d*couplage_plaque_theta);
        }
        else if (multi.GetStringType() == 2)
        {    
          int dof_LM1 = multi.GetLaMuDofNumber(0,g);
          scheme(dof_LM1) -=  BridgeHeight2*(scal_c - scal_d*couplage_plaque_theta);
		
          int dof_LM2 = multi.GetLaMuDofNumber(1,g);
          scheme(dof_LM2) +=  BridgeHeight0*(scal_c - scal_d*couplage_plaque_theta);
		
        }			  
        else if (multi.GetStringType() == 6)
        {
          int dof_LM1 = multi.GetLaMuDofNumber(1,g);
          scheme(dof_LM1) -=  BridgeHeight0*(scal_c - scal_d*couplage_plaque_theta);
		
          int dof_LM2 = multi.GetLaMuDofNumber(0,g);
          scheme(dof_LM2) +=  BridgeHeight1*(scal_c - scal_d*couplage_plaque_theta);
		    }
        else if (multi.GetStringType() == 7) //7 = WAVE_LIN_2T
        {	
          int dof_LM1 = multi.GetLaMuDofNumber(1,g);
          scheme(dof_LM1) -=  BridgeHeight0*(scal_c - scal_d*couplage_plaque_theta);
		
          int dof_LM2 = multi.GetLaMuDofNumber(0,g);
          scheme(dof_LM2) +=  BridgeHeight1*(scal_c - scal_d*couplage_plaque_theta);
		    }		  
        else
        {
          cout << "extra_bridge currently works for system with dimension = 1 or 5" << endl;
          cout << "Please do not include the keyword 'ExtraBridge' in your .ini file " << endl; 
          abort();
        }
      } 

        
	    //equation for the bridge
	    if (size_bridge > 0)
	    {
  	    // we retrieve values of theta at different time steps
  	    Real_wp theta_n = val2(offset_bridge,1);
  	    Real_wp theta_nm1 = val2(offset_bridge,0);
              
              //~ int last_dof_string = multi.GetNbStringDofs(g) - 1;
              
  	    // lagrange multiplier part
  	    //~ for (int k = 2; k < nb_lagrange_mult; k++)
              //~ {
              //~ int dof_LM = multi.GetLagrangeM_DofNumber(g, k);
              //~ Real_wp diffe0 = multi.GetValue(g, k, last_dof_string, 0, val2);
              //~ 
              //~ 
              //~ scheme(dof_LM) = -diffe0*0.5 + theta_nm1*0.5 ;
              //~ }
  	    
  	    // part theta_n
  	    scheme(offset_bridge) = -theta_n*(2*vect_inertia*multi.invDeltat2 - vect_stiff*(1-2*phi_control));
              
  	    //part theta_n-1
  	    scheme(offset_bridge) -= theta_nm1*(-vect_inertia*multi.invDeltat2 + vect_damp*multi.invDeuxDeltat - vect_stiff*phi_control );
	    }
    }
  }

  
  void Bridge
  ::AddDiff(MultiString_Base& multi, VirtualMatrix<Real_wp>& DiffMatrix,
            Matrix<Real_wp>& val2)
  {
    
    int nb_corde = multi.GetNbStrings();
    int dimension = multi.GetNbUnknownsPerString();
    for (int i1 = 0; i1 < nb_corde; i1++)
    {
      if (!multi.IsAttachedString(i1))
      {
        for (int i2 = 0; i2 < nb_corde; i2++)
        {
          for (int g2 = 0; g2 < nb_lagrange_mult; g2++)
          {
            if(extra_bridge)
            {    
              Real_wp BridgeHeight2 = multi.GetBridgeHeight(2);
              Real_wp BridgeHeight1 = multi.GetBridgeHeight(1);	
              Real_wp BridgeHeight0 = multi.GetBridgeHeight(0);	
              if (multi.GetStringType() == 0)  // WAVE_LIN
              {
                multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2)+BridgeHeight2*BridgeHeight2*scal_d, i2, DiffMatrix);	
              }
              else if (multi.GetStringType() == 2) // WAVE_NL
              {                
                if (g2 == 0)
                {
                  multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2)+BridgeHeight2*BridgeHeight2*scal_d, i2, DiffMatrix);
                  int dof_LM = multi.GetLaMuDofNumber(g2, i1);
                  DiffMatrix.SetEntry(dof_LM , dof_LM+1, -BridgeHeight0*BridgeHeight2*scal_d);	
                }	
                else if (g2 == 1)
                {
                  multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2)+BridgeHeight0*BridgeHeight0*scal_d, i2, DiffMatrix);	
                  int dof_LM2= multi.GetLaMuDofNumber(g2, i1);
                  DiffMatrix.SetEntry(dof_LM2 , dof_LM2 -1, -BridgeHeight0*BridgeHeight2*scal_d);	
                }
                else 
                  multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2), i2, DiffMatrix);		
              }						
              else if (multi.GetStringType() == 6) // WAVE_NL_2T
              {              
                if (g2 == 0)
                {
                  multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2)+BridgeHeight1*BridgeHeight1*scal_d, i2, DiffMatrix);
                  int dof_LM = multi.GetLaMuDofNumber(g2, i1);
                  DiffMatrix.SetEntry(dof_LM , dof_LM+1, -BridgeHeight0*BridgeHeight1*scal_d);
                }	
                else if (g2 == 1)
                {
                  multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2)+BridgeHeight0*BridgeHeight0*scal_d, i2, DiffMatrix);	
                  int dof_LM2= multi.GetLaMuDofNumber(g2, i1);
                  DiffMatrix.SetEntry(dof_LM2 , dof_LM2 -1, -BridgeHeight0*BridgeHeight1*scal_d);	
                }
                else 
                  multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2), i2, DiffMatrix);		
              }
              else if (multi.GetStringType() == 7) //7 = WAVE_LIN_2T
              {              
                if (g2 == 0)
                {
                  multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2)+BridgeHeight1*BridgeHeight1*scal_d, i2, DiffMatrix);
                  int dof_LM = multi.GetLaMuDofNumber(g2, i1);
                  DiffMatrix.SetEntry(dof_LM , dof_LM+1, -BridgeHeight0*BridgeHeight1*scal_d);
                }	
                else if (g2 == 1)
                {
                  multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2)+BridgeHeight0*BridgeHeight0*scal_d, i2, DiffMatrix);	
                  int dof_LM2= multi.GetLaMuDofNumber(g2, i1);
                  DiffMatrix.SetEntry(dof_LM2 , dof_LM2 -1, -BridgeHeight0*BridgeHeight1*scal_d);	
                }
                else 
                  multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2), i2, DiffMatrix);		
              }						
              else 
              {
                cout << "extra_bridge currently works for system with dimension = 1 or 5" << endl;
                cout << "Please do not include the keyword 'ExtraBridge' in your .ini file " << endl; 
                abort();
              }						
            }
            else // WAVE_NL, WAVE_NL_TIMO, etc...
            {
              multi.SetDiffLaMuLaMu(g2,i1, vect_scal_b(g2), i2, DiffMatrix);	
            }
          }                  
        }
        //enum {WAVE_LIN, WAVE_LIN_TIMO, WAVE_NL, WAVE_NL_TIMO, WAVE_PIPO, WAVE_NL_2T,WAVE_NL_2T_TIMO,WAVE_LIN_2T};			
            
        int last_dof_string = multi.GetNbStringDofs(i1) - 1;
        Real_wp nu_k;
        Real_wp nu_ortho_k;
            
        for (int k = 0; k < dimension; k ++)
        {
          nu_k = multi.GetTransmissionChevaletString(i1, 0, k);
          multi.SetDiffLaMuString(0, i1, nu_k*0.5, i1, k, last_dof_string, DiffMatrix);
              
          for (int g2 = 1; g2 < nb_lagrange_mult; g2++)
          {	
            nu_ortho_k = multi.GetTransmissionChevaletString(i1, g2, k);        
            multi.SetDiffLaMuString(g2, i1, nu_ortho_k*0.5, i1, k, last_dof_string, DiffMatrix);	
          }      
        }  
      }
           
	    if(size_bridge > 0)
	    {    
	      int last_dof_string = multi.GetNbStringDofs(i1) - 1;
            
 	      DiffMatrix.SetEntry(offset_bridge, offset_bridge, vect_inertia*multi.invDeltat2+vect_damp*multi.invDeuxDeltat+vect_stiff*phi_control); 
	      DiffMatrix.SetEntry(offset_bridge, multi.GetLagrangeM_DofNumber(i1,0), height_of_bridge);
	    
        if (dimension > 2)
        {
          int dof_LM = multi.GetLagrangeM_DofNumber(i1, 2);
          int row_num = multi.GetStringDofNumber(i1, 2, last_dof_string);
          DiffMatrix.SetEntry(dof_LM, row_num, 0.5);
          DiffMatrix.SetEntry(dof_LM, offset_bridge, -0.5);
          DiffMatrix.SetEntry(offset_bridge, multi.GetLagrangeM_DofNumber(i1,2), 1);			
        }
	    }
    }
  }
  
  Real_wp Bridge::Compute_coef_couplage_V(MultiString_Base& multi, Matrix<Real_wp>& val2)
  {
    int nb_corde = multi.GetNbStrings();      
    coef_couplage_V = 0.0;
    for (int p = 0; p < nb_corde; p++)
      coef_couplage_V -= multi.GetLaMu(0,p, 2, val2);
    
    return(coef_couplage_V);
  }
  
  Real_wp Bridge::Compute_coef_couplage_H(MultiString_Base& multi, Matrix<Real_wp>& val2)
  {
    int nb_corde = multi.GetNbStrings();      
    coef_couplage_H = 0.0;
    if (nb_lagrange_mult > 1)
    {
      for (int p = 0; p < nb_corde; p++)
        coef_couplage_H -= multi.GetLaMu(1, p, 2, val2);
    }
    return(coef_couplage_H);
  }
  
  
  Real_wp Bridge::Compute_coef_couplage_any(int g2, MultiString_Base& multi, Matrix<Real_wp>& val2)
  {
    int nb_corde = multi.GetNbStrings();     
    
    coef_couplage_any = 0.0;
    
    for (int p = 0; p < nb_corde; p++)
      coef_couplage_any -= multi.GetLaMu(g2,p, 2, val2);
    
    return(coef_couplage_any);
  } 

  Real_wp Bridge::Compute_coef_couplage_Theta(MultiString_Base& multi, Matrix<Real_wp>& val2)
  {
    int nb_corde = multi.GetNbStrings();     
    coef_couplage_any = 0.0;
    
    for (int p = 0; p < nb_corde; p++)
    {
      if (multi.GetStringType() == 0)
      {  
        Real_wp BridgeHeight = multi.GetBridgeHeight(2);
        coef_couplage_any -= multi.GetLaMu(0,p, 2, val2)*BridgeHeight;
      }
      else if (multi.GetStringType() == 2)
      {
        Real_wp BridgeHeight0 = multi.GetBridgeHeight(0);
        Real_wp BridgeHeight2 = multi.GetBridgeHeight(2);	
        coef_couplage_any -= ( multi.GetLaMu(0,p, 2, val2)*BridgeHeight2 - multi.GetLaMu(1,p,2,val2)*BridgeHeight0 );
      }		
      else if (multi.GetStringType() == 6)
      {
        Real_wp BridgeHeight0 = multi.GetBridgeHeight(0);
        Real_wp BridgeHeight1 = multi.GetBridgeHeight(1);	
        coef_couplage_any -= ( multi.GetLaMu(1,p, 2, val2)*BridgeHeight0 - multi.GetLaMu(0,p,2,val2)*BridgeHeight1 );
      }
      else if (multi.GetStringType() == 7) //7 = WAVE_LIN_2T
      {
        Real_wp BridgeHeight0 = multi.GetBridgeHeight(0);
        Real_wp BridgeHeight1 = multi.GetBridgeHeight(1);	
        coef_couplage_any -= ( multi.GetLaMu(1,p, 2, val2)*BridgeHeight0 - multi.GetLaMu(0,p,2,val2)*BridgeHeight1 );
      }		
    }
    return(coef_couplage_any);
  }      
  
}

#define MONTJOIE_FILE_BRIDGE_CXX
#endif

