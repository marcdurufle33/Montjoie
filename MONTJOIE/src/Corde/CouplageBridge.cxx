#ifndef MONTJOIE_FILE_COUPLAGE_BRIDGE_CXX

namespace Montjoie
{

  //! default constructor
  CouplageBridge::CouplageBridge(const string& data_file)
  {
    file_output_energy = "EnergyCouplageBridge.txt";
  }


  void CouplageBridge::SetInputData(const string & description_field, const Vector<string> &param)
  {
		if (!description_field.compare("Impedance0"))
          {
			  if (note.nb_lagrange_mult > 0)
			  {
				oscillator[0].SetInputData("Impedance", param);
				cout << "oscillator0 is used" << endl;
			  }
          }
        else if(!description_field.compare("Impedance1"))
          {
			  if (note.nb_lagrange_mult > 1)
			  {
				oscillator[1].SetInputData("Impedance", param);
				cout << "oscillator1 is used" << endl;
			  }
          }
	    else if(!description_field.compare("Impedance2"))
          {
			  if (note.nb_lagrange_mult > 2)
			  {
				oscillator[2].SetInputData("Impedance", param);
				cout << "oscillator2 is used" << endl;
			  }
          }
	    else if(!description_field.compare("Impedance3"))
          {
			  if (note.nb_lagrange_mult > 3)
              {
				oscillator[3].SetInputData("Impedance", param);
				cout << "oscillator3 is used" << endl;
			  }
          }
		else if(!description_field.compare("Impedance4"))
          {
			  if (note.nb_lagrange_mult > 4)
              {
				oscillator[4].SetInputData("Impedance", param);
				cout << "oscillator4 is used" << endl;
			  }
          }
		else if(!description_field.compare("ImpedanceT"))
          {
			  if (note.bridge.extra_bridge)
              {
				oscillatorTheta.SetInputData("Impedance", param);
				
				cout << "oscillatorT is used" << endl;
			  }
          }
        else if(!description_field.compare("SismoLambdaExtra"))
          {
			  if (note.bridge.extra_bridge)
			  {
				oscillatorTheta.SetInputData("SismoLambda",param);
				cout << "oscillatorTheta will be output-ed" << endl;
			  }
		  }          
         
          

  }


  //! constructs everything needed to complete the simulation (such that AdvanceStringBridge will work)
  void CouplageBridge::ConstructAll(const string& data_file)
  {
    // the data file is read
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    ReadLinesFile(data_file, lines_data_file, MPI_COMM_WORLD);
#else
    ReadLinesFile(data_file, lines_data_file);
#endif

    ReadInputFile(lines_data_file, note);
    
    // constructing the strings
    note.name_data_file = data_file;
    note.ConstructAll();
    cout << "String ConstructAll : OK! " << endl;
    note.Init();
    cout << "Strings : OK" << endl;

    output_ener_all_but_note.Init(note.DOSSIER+file_output_energy, note.size_buffer_ener);
  
    
    ReadInputFile(lines_data_file, *this);
    Real_wp scal_b;

    for(int g2 = 0; g2 < note.nb_lagrange_mult; g2++)
	{
		oscillator[g2].Deltat=note.Deltat;
		oscillator[g2].ReadImpedance(note.DOSSIER);
		scal_b = oscillator[g2].ConstruitMatricesResolExacte();
		note.bridge.SetVecScal_b(g2,scal_b);
		oscillator[g2].ConstruitLambda0etPoint();
		
	}
	if (note.bridge.extra_bridge)
	{
		Real_wp scal_d;
		oscillatorTheta.Deltat=note.Deltat;
		oscillatorTheta.ReadImpedance(note.DOSSIER);
		scal_d = oscillatorTheta.ConstruitMatricesResolExacte();
		note.bridge.SetScal_d(scal_d);
		oscillatorTheta.ConstruitLambda0etPoint();
		cout << "Oscillator : OK"<< endl;
	}
  }


  //! Advancing the time scheme for the strings and the bridge
  /*!
    \param[in] nt iteration number
    \param[in] time current time t^n
    \param[in] save_snapshot if true, snapshots may be written
   */
  void CouplageBridge::AdvanceStringBridge(int nt, const Real_wp& temps, bool save_snapshot)
  {
    if (save_snapshot)
	{
	  oscillator[0].WriteOutputLambda(nt, temps);
  	  if (note.bridge.extra_bridge)
	    oscillatorTheta.WriteOutputLambda(nt, temps);
	}
	  
	  
	Real_wp scal_a;

    for(int g2 = 0; g2 < note.nb_lagrange_mult; g2++)
	{
		scal_a = oscillator[g2].CalculeScal_a();
		note.bridge.SetVecScal_a(g2,scal_a);
	}
	
	if (note.bridge.extra_bridge)
	{
	  Real_wp scal_c;
	  scal_c = oscillatorTheta.CalculeScal_a();
	  note.bridge.SetScal_c(scal_c);
	}
    // the string is solved 
    note.Advance(temps, nt, save_snapshot);
    
    // then we update the values of lambda and time derivative of lambda   
	for(int g2 = 0; g2 < note.nb_lagrange_mult; g2++)
	{	
		oscillator[g2].ActualiseLambda(note.bridge.Compute_coef_couplage_any(g2, note, note.val));		
	}
	
	if (note.bridge.extra_bridge)
	  oscillatorTheta.ActualiseLambda(note.bridge.Compute_coef_couplage_Theta(note, note.val));
  }
  
  
  //! Displays parameters of the class
  void CouplageBridge::Display()
  {

      {
        string file_name = note.DOSSIER + note.file_storage_param;
        ofstream file_out(file_name.data());
        note.Display(file_out, true);

        
        for(int g2 = 0; g2 < note.nb_lagrange_mult; g2++)
			oscillator[g2].Display(file_out,true);
			
		if(note.bridge.extra_bridge)
			oscillatorTheta.Display(file_out,true);
			
		DISP(note.nb_lagrange_mult);
		
        
        file_out.close();
      }
  }
  

  //! Computes energy of the different sub-systems
  Real_wp CouplageBridge::WriteEnergy()
  {
    Real_wp deltat_ener = note.deltat_ener;
    Real_wp t_courant = note.t_courant;
    Real_wp t_begin_ener = note.t_begin_ener;
    Real_wp t_end_ener = note.t_end_ener;
    Real_wp Deltat = note.Deltat;
    
    Real_wp ener_plaque_potential_V = 0.0;
    Real_wp ener_plaque_kinetic_V = 0.0;

    Real_wp ener_plaque_potential = 0.0;
    Real_wp ener_plaque_kinetic = 0.0;
        
    Real_wp ener_cordes = 0.0;
    VectReal_wp energy_sub;
    if (deltat_ener != 0)
      ener_cordes = note.ComputeEnergy(energy_sub);
    
    Real_wp ener_tot = 0;
    Real_wp t_display = note.t_display;
    
    if (deltat_ener != 0)
      {
        int n = toInteger(round(t_courant/deltat_ener));
        TinyVector<Real_wp, 6> my_ener;//, previous_energy_tmp;
        bool write_energy = false;
        
        if ((t_courant>= t_begin_ener) && (t_courant<= t_end_ener)
            && (abs(t_courant/deltat_ener - n) <(Deltat*0.5)))
          {
            for (int g2 = 0; g2 < note.nb_lagrange_mult; g2++)
            {
				ener_plaque_kinetic += oscillator[g2].GetKineticEnergy(); 
				ener_plaque_potential += oscillator[g2].GetPotentialEnergy();
			}
                     
            my_ener(0) = ener_plaque_kinetic;
            my_ener(1) = ener_plaque_potential;
            my_ener(2) = ener_cordes;
            write_energy = true;
            my_ener(3) = ener_plaque_kinetic+ener_plaque_potential+ener_cordes;
            
            my_ener(4) = ener_plaque_kinetic_V;
            my_ener(5) = ener_plaque_potential_V;
            ener_tot = my_ener(3);
          }       
        
        if (write_energy)
          {
            n = toInteger(round(t_courant/(t_display)));
            if (abs(t_courant/t_display - n) < Deltat/2)   
              {
                cout << "At time t = " << t_courant << " total energy (string & bridge) = " << ener_tot << endl;
              }

            output_ener_all_but_note.AddTinyVect(my_ener);            
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

                for (int i = 0; i < oscillator[0].Lambda.GetM(); i++)
                  if (abs(oscillator[0].Lambda(i)) > Lmax)
                    Lmax = abs(oscillator[0].Lambda(i));                                          
                
                cout << "At time t = " << t_courant << " || u,v,phi ||_inf = "
                     << vmax << ",  || Lambda||_inf = " << Lmax << endl;
              }
          }
      }
    
    return ener_tot;
  }
  
}

#define MONTJOIE_FILE_COUPLAGE_BRIDGE_CXX
#endif
