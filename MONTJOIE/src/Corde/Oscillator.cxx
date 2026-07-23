#ifndef MONTJOIE_FILE_OSCILLATOR_CXX

namespace Montjoie
{
  
  /*******************
   * Oscillator_Base *
   *******************/


  Oscillator_Base::Oscillator_Base()
  {
    t_begin_output_lambda = 0;
    t_end_output_lambda = 0;
    deltat_output_lambda = 0;
    size_buffer_lambda = 100;
    N_output_Lambda= -1;
    vect_second_membre_is_given = false;
    write_oscillator_param_bool= false;
  }

  void Oscillator_Base::Display(ostream& file_out, bool write_on_file){
    cout << " ---------- Oscillator Parameters : ----------- " <<endl;
    DISP(mass_file_name);
    DISP(vect_masse);
    DISP(lambda_file_name);
    DISP(vect_lambda);
    DISP(vect_second_membre_file_name);
    DISP(VectSecondMembreModal);
    cout << "Oscillator frequency is: " << sqrt(vect_lambda(0)/vect_masse(0))/2/pi_wp << endl;
    cout << " ---------------------------------------------- " <<endl;
  }
  
  
  void Oscillator_Base
  ::SetInputData(const string & description_field, const VectString & parameters)
  {
    if (!description_field.compare("Impedance"))
    {
  	  if (parameters.GetM() <= 2)
  	  {
  	    cout << "In SetInputData of Oscillator_Base" << endl;
  	    cout << "Impedance needs three parameters, for instance :" << endl;
  	    cout << "Impedance = mass lambda amo" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
  	  }

      mass_file_name = parameters(0);
      lambda_file_name = parameters(1);
      amo_file_name = parameters(2);
      if(parameters.GetM()>3)
      {
        vect_second_membre_file_name = parameters(3);
        vect_second_membre_is_given = true;
      }
    }
    else if(!description_field.compare("SismoLambda"))
    {
  	  if (parameters.GetM() <= 3)
  	  {
  	    cout << "In SetInputData of Soundboard" << endl;
  	    cout << "SismoLambda needs three parameters, for instance :" << endl;
  	    cout << "SismoLambda = t0 tf dt size_buffer" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
  	  }

      t_begin_output_lambda = to_num<Real_wp>(parameters(0));
      t_end_output_lambda = to_num<Real_wp>(parameters(1));
      deltat_output_lambda = to_num<Real_wp>(parameters(2));
      size_buffer_lambda = to_num<int>(parameters(3));
    }
    else if(!description_field.compare("SismoLambdaTruncate"))
    {
  	  if (parameters.GetM() <= 0)
  	  {
  	    cout << "In SetInputData of Soundboard" << endl;
  	    cout << "SismoLambdaTruncate needs one parameter, for instance :" << endl;
  	    cout << "SismoLambdaTruncate = N" << endl;
  	    cout << "Current parameters are : " << endl << parameters << endl;
  	    abort();
  	  }
	
      N_output_Lambda = to_num<int>(parameters(0));
    }
  }
  

  //! Lecture des coefficients de masse, d'amortissement et rigidite
  void Oscillator_Base::ReadImpedance(const string& DOSSIER)
  {
    vect_lambda.ReadText(lambda_file_name);
    nev = vect_lambda.GetSize();
	
    vect_amortissement.ReadText(amo_file_name);
    if (vect_amortissement.GetSize() != nev)
    {
      cout << "Probleme taille " << amo_file_name << endl;
      abort();
    }
    
    vect_masse.ReadText(mass_file_name);
    if(vect_masse.GetSize() != nev)
    {
      cout << "Probleme taille " << mass_file_name << endl;
      abort();
    }
    if(vect_second_membre_is_given)
    {
      VectSecondMembreModal.ReadText(vect_second_membre_file_name);
      if(VectSecondMembreModal.GetSize() != nev)
      {
        cout << "Probleme taille " << vect_second_membre_file_name << endl;
        abort();
      }
      cout << "***************************" << endl;
    }
    else
    {
      cout << "***************************" << endl;
           
      VectSecondMembreModal.Reallocate(nev);
      VectSecondMembreModal.Fill(1.0);
    }

    if((N_output_Lambda==-1)||(N_output_Lambda>nev))
      N_output_Lambda = nev;
    
    if (deltat_output_lambda != 0)
    {
      output_lambda.  Init(DOSSIER + "SismoLambda.bin", size_buffer_lambda, true);
      output_lambda.SetBinary(true);
      output_lambda.SetDoublePrecision(false);
      
      VectReal_wp freq_propres(nev), alpha(nev);
      for(int g = 0; g < nev; g++)
      {
        Real_wp invmasse = 1.0/vect_masse(g);
        freq_propres(g) = sqrt(vect_lambda(g)*invmasse)/(2.0*pi_wp);
        alpha(g) = vect_amortissement(g)*invmasse;
      }
        
      output_lambda.AddVect(freq_propres);
      output_lambda.AddVect(alpha);
    }
  }
  

  //! Computation of R_{\Delta t} for the analytical computation of the solution
  Real_wp Oscillator_Base::ConstruitMatricesResolExacte()
  {
    // avoir rempli les valeurs propres
    MatriceResolExacte0.Reallocate(nev); MatriceResolExacte0.Fill(Real_wp(0));
    MatriceResolExacte1.Reallocate(nev); MatriceResolExacte1.Fill(Real_wp(0));
    MatriceResolExacteSecondMembre.Reallocate(nev);
    MatriceResolExacteSecondMembre.Fill(Real_wp(0));
    MatriceResolExacte0Point.Reallocate(nev); MatriceResolExacte0Point.Fill(Real_wp(0));
    MatriceResolExacte1Point.Reallocate(nev); MatriceResolExacte1Point.Fill(Real_wp(0));
    MatriceResolExacteSecondMembrePoint.Reallocate(nev);
    MatriceResolExacteSecondMembrePoint.Fill(Real_wp(0));
    Real_wp scal_b = 0.0;


    for (int g = 0; g < nev; g++)
    {
      if (vect_masse(g)==0)
      {
        cout << "Probleme : masse nulle au "<< g << "ieme element"<<endl;
        abort();
      }
        
      Real_wp invmasse =1.0/vect_masse(g);
      Real_wp lambda = sqrt(vect_lambda(g)*invmasse);
      // sqrt(eigen_values(g));
      Real_wp alpha = vect_amortissement(g)*0.5*invmasse;
      
      Real_wp lambda2 = lambda*lambda;
      Real_wp invlambda2 = 1.0/lambda2;
      Real_wp lambdatilde,invlambdatilde;
      Real_wp co,si;
      Real_wp expo= exp(-alpha*Deltat);
      Real_wp coef;
        
      Real_wp Delta = abs(alpha)-lambda;
        
      if (Delta < 0)
      {
        lambdatilde = sqrt(lambda2 - alpha*alpha);
        invlambdatilde = 1.0/lambdatilde;
        co = cos(lambdatilde*Deltat);
        si = sin(lambdatilde*Deltat);		
        coef = -1;
      }
      else
      {
        lambdatilde = sqrt( alpha*alpha-lambda2);
        invlambdatilde = 1.0/lambdatilde;
        co = cosh(lambdatilde*Deltat);
        si = sinh(lambdatilde*Deltat);
        coef = 1;
      }
        
      MatriceResolExacte0(g) = expo*(co + alpha*invlambdatilde*si);
      MatriceResolExacte1(g) = expo*invlambdatilde*si;
      MatriceResolExacteSecondMembre(g) = invlambda2*invlambdatilde
          *(lambdatilde*(1-expo*co)-alpha*expo*si)*invmasse;
        
      MatriceResolExacte0Point(g) = -alpha*MatriceResolExacte0(g) 
          + expo*( coef*lambdatilde*si + alpha*co);
        
      MatriceResolExacte1Point(g) = - alpha*MatriceResolExacte1(g)+ expo*co;
      MatriceResolExacteSecondMembrePoint(g) = expo*invlambdatilde*si*invmasse;

      scal_b += MatriceResolExacteSecondMembre(g)*VectSecondMembreModal(g)
                    *VectSecondMembreModal(g); 
      // scal_b  sera ecrasé si on invoque ConstruitLambdaSchur
    }

    return(scal_b);
  }
  
  
  //! Construction of modal right hand sides
  void Oscillator_Base::ConstruitLambda0etPoint()
  {
    // Lambda0 =                R_∆t tr(P_h^U) J_h
    // Lambda0Point = \dot{R}_∆t tr(P_h^U) J_h
    // avoir déjà les matrices de résolution
    Lambda0.Reallocate(nev);
    Lambda0Point.Reallocate(nev);
    for (int g = 0; g < nev; g++)
    {
      Lambda0(g) = MatriceResolExacteSecondMembre(g) * VectSecondMembreModal(g);
      Lambda0Point(g) = MatriceResolExacteSecondMembrePoint(g) * VectSecondMembreModal(g);
    }
    
    Lambda.Reallocate(nev);		Lambda.Fill(Real_wp(0));
    LambdaPoint.Reallocate(nev);	LambdaPoint.Fill(Real_wp(0));
    
    // verify that N_output_Lambda is not to big
    if((N_output_Lambda==-1)||(N_output_Lambda>nev))
    {
      // Attention ca va faire un super gros vecteur !
      N_output_Lambda = nev;
    }  

 }


  //! pour le couplage avec la corde
  Real_wp Oscillator_Base::CalculeScal_a()
  {
    VectReal_wp V;
    V.Reallocate(nev);
    Real_wp scal_a = 0.0;
    for (int g = 0; g < nev; g++)
    {
      V(g) = MatriceResolExacte0(g) * Lambda(g)
          + MatriceResolExacte1(g) * LambdaPoint(g) - Lambda(g);
        
      scal_a += V(g)*VectSecondMembreModal(g);
    }
    
    return(scal_a);
  }
  
  
  //! Calcul de l'énergie
  Real_wp Oscillator_Base::GetKineticEnergy()
  {
    Real_wp res = 0.0;	
    for (int g = 0; g < nev; g++)
      res+= LambdaPoint(g)*LambdaPoint(g)*vect_masse(g);
    
    return(res*0.5);
  }
  
  
  //! energie potentielle
  Real_wp Oscillator_Base::GetPotentialEnergy()
  {
    Real_wp res = 0.0;	
    for (int g = 0; g < nev; g++)
      res+= Lambda(g)*Lambda(g)*vect_lambda(g);
    
    return(res*0.5);
  }
  

  //! Avance d'un pas de temps
  void Oscillator_Base::ActualiseLambda(Real_wp coef)
  {
    // avance d'un pas
    VectReal_wp Lambdaold(Lambda);
    VectReal_wp LambdaPointold(LambdaPoint);
    Lambda.Fill(Real_wp(0));
    LambdaPoint.Fill(Real_wp(0));
	
    for (int g = 0; g < nev; g++)
    {
      Lambda(g) = MatriceResolExacte0(g) * Lambdaold(g)
          + MatriceResolExacte1(g) * LambdaPointold(g) + coef*Lambda0(g);
        
      LambdaPoint(g) = MatriceResolExacte0Point(g) * Lambdaold(g)
          + MatriceResolExacte1Point(g) * LambdaPointold(g) + coef*Lambda0Point(g);
    }
  }

  
  //! sortie des instantanes, directement les composantes sur les modes propres
  void Oscillator_Base::WriteOutputLambda(int nt, Real_wp temps)
  {
    if(deltat_output_lambda != Real_wp(0.0))
    {
      int n = toInteger(round(temps / deltat_output_lambda));
      if((temps >= t_begin_output_lambda) && (temps <= t_end_output_lambda)
           && (abs(temps/deltat_output_lambda - n) <(Deltat*0.5)))
      {
        VectReal_wp vect_output(N_output_Lambda+1);
        //  VectReal_wp vect_output(Lambda.GetM()+1);
        vect_output(0) = temps;
        for (int i = 0; i < N_output_Lambda; i++)
          vect_output(i+1) = Lambda(i);
          
        output_lambda.AddVect(vect_output);
      }
    }
  }
  
}

#define MONTJOIE_FILE_OSCILLATOR_CXX
#endif

