#ifndef MONTJOIE_FILE_HELMHOLTZ_RADIAL_CXX

#include "HelmholtzRadialInline.cxx"

namespace Montjoie
{

  //! Default constructor
  AtmosphereImpedanceCondition::AtmosphereImpedanceCondition()
  {
    order = 0;
    alpha = Real_wp(0);
    exact_condition = false;
    abc_condition = false;
    ha_condition = false; // NATHAN
    threshold = 1e-20;
  }


  //! Reads parameters of the equivalent atmosphere condition
  int AtmosphereImpedanceCondition
    ::ReadParameters(int& nb, const VectString& parameters, int side)
    {
      nb++;
      if (parameters.GetM() < nb+2)
      {
        cout << "For ATMOSPHERE, give the order and alpha" << endl;
        abort();
      }

      exact_condition = false;
      abc_condition = false;    
      ha_condition = false; //NATHAN
      if (parameters(nb) == "EXACT")
      {
        exact_condition = true;	
        nb++;
        order = 0;
        if (parameters.GetM() > nb+1)
          order = to_num<int>(parameters(nb++));
      }
      else if (parameters(nb) == "ABC")
      {
        nb++;
        abc_condition = true;
        if(parameters(nb) == "HA"){ //NATHAN
          nb++;
          ha_condition = true;
        }
        order = to_num<int>(parameters(nb++));
        DISP(order);
      }
      else
        order = to_num<int>(parameters(nb++));

      alpha = to_num<Real_wp>(parameters(nb++));
      DISP(alpha);
      return BoundaryConditionEnum::LINE_IMPEDANCE;
    }


  //! Updates impedance corresponding to the atmosphere condition
  void AtmosphereImpedanceCondition::UpdateImpedance(const Real_wp& rt, const Complex_wp& k0,
      Complex_wp& impedance)
  {
    if (order <= 0 && !abc_condition)
      return;

    if (abc_condition && !ha_condition) //NATHAN
    {
      if(order==0)
      {
        impedance = -Iwp*k0; DISP(impedance);
        return;
      }
      else if (order==1)
      {
        Complex_wp kinf2 = k0*k0;
        //	    Complex_wp det = sqrt(Real_wp(4)*kinf2 - alpha*alpha);
        // Complex_wp km = Real_wp(0.5)*(Iwp*alpha - det);
        //	    Complex_wp kp = Real_wp(0.5)*(Iwp*alpha + det);
        //	    impedance = -Iwp*kp; DISP(impedance);
        impedance =-( -alpha*0.5+Iwp*k0*sqrt(Real_wp(1)-alpha*alpha*0.25/kinf2));
        //  DISP(impedance);
        return;
      }
      else if(order==2)
      {
        // abc condition with curvature
        //	    Real_wp alpha2 = alpha + Real_wp(2)/rt;
        Complex_wp kinf2 = k0*k0;
        //  Complex_wp det = sqrt(Real_wp(4)*kinf2 - alpha2*alpha2);
        // Complex_wp km = Real_wp(0.5)*(Iwp*alpha - det);
        //Complex_wp kp = Real_wp(0.5)*(Iwp*alpha2 + det);
        //impedance = -Iwp*kp; DISP(impedance);
        Real_wp kappa = Real_wp(1)/rt;
        impedance = -( -(kappa + alpha*0.5) + Iwp*k0*sqrt(Real_wp(1)-(alpha*kappa+alpha*alpha*0.25)/kinf2) ) ;
        DISP(impedance);
        return;
      }
      else if(order>=3)
      {
        // abc condition with curvature and surfacic laplace term
        //	    alpha = alpha + 2.0/rt;
        //Complex_wp kinf2 = k0*k0 - l*(l+1)/(rt*rt);
        //Complex_wp det = sqrt(Real_wp(4)*kinf2 - alpha*alpha);
        // // Complex_wp km = Real_wp(0.5)*(Iwp*alpha - det);
        // Complex_wp kp = Real_wp(0.5)*(Iwp*alpha + det);
        //impedance = -Iwp*kp; DISP(impedance);
        //return;

        // we will fill the coef inside GetLaplacianImpedance
        SetComplexZero(impedance);
        return;
      }
    }

    if (order == 3) 
    {
      SetComplexZero(impedance);
      return;
    }

    if(abc_condition && ha_condition) // NATHAN
    {
      if(order == 2)
      {

        // DISP("Impedance HA");

        Real_wp kappa = 1.0/rt;
        impedance = -( -(kappa+ 0.5*alpha) + Iwp*k0*sqrt(1-0.25*alpha*alpha/(k0*k0)*(1.0-(alpha*kappa)/(k0*k0-0.25*alpha*alpha))));
        // DISP(impedance);
//          std::ofstream file;
//          file.open("imped.txt",std::ios::app);
//          file << std::setprecision(12) << impedance.real() << " " << impedance.imag() <<std::endl;
        return;
      }
      else if(order>=3)
      {

        // we use GetLaplacianImpedance
        SetComplexZero(impedance);
        return;
      }
    }

    // conditions de Victor    
    Real_wp delta = Real_wp(1)/alpha;
    Complex_wp coef; SetComplexOne(coef);
    Real_wp one(1);
    switch (order)
    {
      case 1:
        coef = delta;
        break;
      case 2:
        coef = delta*(one-Real_wp(2)*delta/rt);
        break;
    }

    impedance = one/coef;
  }


  //! returns true if a Laplacian term is involved
  bool AtmosphereImpedanceCondition::LaplacianPresent() const
  {
    if (exact_condition)
      return false;

    if (order == 3)
      return true;

    if (order == 4)
      return true;

    if (order == 5)
      return true;

    if (order == 6)
      return true;

    return false;
  }


  //! returns true if the exact condition has to be used
  bool AtmosphereImpedanceCondition::ExactCondition() const
  {
    return exact_condition;
  }


  //! returns the impedance depending on l
  Complex_wp AtmosphereImpedanceCondition
    ::GetLaplacianImpedance(const Real_wp& factor_laplace, const Real_wp& rt, 
        const Complex_wp& kinf)
    {
      Complex_wp coef; SetComplexZero(coef);
      if (order < 3)
        return coef;

      if (abc_condition && !ha_condition) //NATHAN
      {
        if (order == 3) 
        {
          // localized for high frequencies (first order in 1/omega)
          //	    Real_wp alpha_loc = alpha + Real_wp(2)/rt;
          //Complex_wp det = kinf - (factor_laplace / (rt*rt) + alpha_loc*alpha_loc/Real_wp(4) ) /(Real_wp(2) * kinf);
          //Complex_wp kp = Real_wp(0.5)*Iwp*alpha_loc + det;
          //coef = -Iwp*kp; DISP(coef);	
          Real_wp kappa = Real_wp(1)/rt;
          coef = - ( -(kappa + Real_wp(0.5)*alpha) + Iwp*kinf  +( factor_laplace/(rt*rt) + alpha*alpha*Real_wp(0.25) + alpha*kappa ) / (Real_wp(2)*Iwp*kinf ) );
           DISP(coef);
          return coef;
        }
        else if (order == 4)
        {
          // non local
          //Real_wp alpha_loc = alpha + Real_wp(2)/rt;
          //Complex_wp kinf2 = kinf*kinf - factor_laplace/(rt*rt);
          //Complex_wp det = sqrt(Real_wp(4)*kinf2 - alpha_loc*alpha_loc);
          // Complex_wp km = Real_wp(0.5)*(Iwp*alpha_loc - det);
          //Complex_wp kp = Real_wp(0.5)*(Iwp*alpha_loc + det);
          //coef = -Iwp*kp; DISP(coef);	
          Real_wp kappa = Real_wp(1)/rt;
          coef = -( -(kappa +Real_wp(0.5)*alpha ) +Iwp*kinf*sqrt(Real_wp(1)-(factor_laplace/(rt*rt) + alpha*alpha*Real_wp(0.25) + alpha*kappa )/(kinf*kinf) ) );
            DISP(coef);
          return coef;
        }
        else if (order == 5)
        {
          // microlocal analysis
          //	    Real_wp alpha_loc = alpha + Real_wp(2)/rt;
          //Complex_wp det = kinf -( Real_wp(1)+Real_wp(1)/(rt*Iwp*kinf)) *(factor_laplace / (rt*rt) ) /(Real_wp(2) * kinf);
          //Complex_wp kp = Real_wp(0.5)*Iwp*alpha_loc + det - alpha*(Real_wp(1)/rt + alpha*Real_wp(0.25) - Real_wp(0.5)/(Iwp*kinf*rt*rt))/(Real_wp(2)*kinf);
          //coef = -Iwp*kp; DISP(coef);	
          Real_wp kappa = Real_wp(1)/rt;
          Complex_wp Beltra =  ( Real_wp(1)+ kappa/(Iwp*kinf) )*(factor_laplace / (rt*rt) ) /(Real_wp(2) *Iwp* kinf );
          Complex_wp correc = alpha*(kappa + alpha*Real_wp(0.25) +kappa*(Real_wp(1)- Real_wp(2.0)*kappa)/(Real_wp(2.0)*Iwp*kinf))/(Real_wp(2)*Iwp*kinf) ;
          coef = - ( -(kappa+alpha*Real_wp(0.5)) +Iwp*kinf + Beltra + correc) ;
          //  DISP(coef);
          return coef;	    
        }
        else if (order == 6)
        {
          // half localized analysis
          //Real_wp alpha_loc = alpha + Real_wp(2)/rt;
          //Complex_wp det = kinf*sqrt(Real_wp(1)-alpha_loc*alpha_loc/(Real_wp(4)*kinf*kinf)) - (factor_laplace / (rt*rt) ) /(Real_wp(2) * kinf);
          //Complex_wp kp = Real_wp(0.5)*Iwp*alpha_loc + det;
          //coef = -Iwp*kp; DISP(coef);	
          Real_wp kappa = Real_wp(1)/rt;
          Complex_wp term_sqrt = sqrt(Real_wp(1.0)-( alpha*alpha*Real_wp(0.25) + alpha*kappa  )/(kinf*kinf) );
          Complex_wp Beltra = factor_laplace/(rt*rt)/(Real_wp(2)*Iwp*kinf*term_sqrt) ;
          coef = -( -(kappa +Real_wp(0.5)*alpha ) +Iwp*kinf*term_sqrt + Beltra );
          //  DISP(coef);
          return coef;	    
        }
      }


      if(abc_condition && ha_condition) // NATHAN
      {
        if(order == 3)
        {
          // DISP("Impedance HA");DISP(order);
          Real_wp kappa = Real_wp(1)/rt;
          Real_wp alpha2 = Real_wp(0.5)*alpha;
          Complex_wp denom = Real_wp(2)*(kinf*kinf-alpha2*alpha2);

          Complex_wp racine = Complex_wp(1)-alpha2*alpha2/(kinf*kinf);
          Complex_wp fact1 = Complex_wp(1)-(alpha*kappa)/denom;
          Complex_wp fact2 = Complex_wp(1)-(Complex_wp(1)/denom)*(factor_laplace/(rt*rt));

          coef = -(-(kappa+alpha2)+Iwp*kinf*sqrt(racine)*fact1*fact2);

          DISP(coef);


          return coef;
        }

        if(order == 4)
        {
           DISP("Impedance HA");DISP(order);
          Real_wp kappa = Real_wp(1)/rt;
          Real_wp alpha2 = Real_wp(0.5)*alpha;
          Complex_wp denom = Real_wp(2)*(kinf*kinf-alpha2*alpha2);

          Complex_wp racine = Complex_wp(1)-alpha2*alpha2/(kinf*kinf);
          Complex_wp ord1 = Complex_wp(1)-(alpha*kappa)/denom;
          Complex_wp ord2 = Complex_wp(1)-(Complex_wp(1)/denom)*(factor_laplace/(rt*rt));
          Complex_wp ord3 = (alpha*alpha)/(rt*rt);
          ord3 += Real_wp(4)*alpha/rt*(factor_laplace/(rt*rt));
          ord3 += (factor_laplace*factor_laplace)/(rt*rt*rt*rt);
          ord3 /= (Real_wp(2)*denom*denom);

          coef = -(-(kappa+alpha2)+Iwp*kinf*sqrt(racine)*ord1*(ord2-ord3));

           DISP(coef);

          return coef;
        }


        if(order == 5)
        {
          //Non-local with Pade
          DISP("Impedance HA Pade");
          Real_wp kappa = Real_wp(1)/rt;
          Real_wp alpha2 = Real_wp(0.5)*alpha;
          Complex_wp kc = kinf*kinf-Real_wp(0.25)*alpha*alpha;

          Complex_wp up = Real_wp(2)*factor_laplace*kc + Real_wp(0.25)*alpha*alpha;
          Complex_wp down = Real_wp(4)*rt*rt*kc*kc - (Real_wp(4)*rt*alpha+factor_laplace)*kc + Real_wp(0.75)*alpha*alpha;

          coef = (kappa+alpha2) - Iwp*(sqrt(kc)-alpha2*kappa/sqrt(kc))*(-Real_wp(1) + up/down);
          DISP(coef);
          return coef;
        }
      }
      Real_wp delta = Real_wp(1)/alpha;
      Real_wp one(1);
      coef = (one + Real_wp(2)*delta/rt -
          square(delta)*(Real_wp(2.0)/(rt*rt) + square(kinf) - factor_laplace/(rt*rt))) / delta;

      return coef;
    }


  //! computes the exact impedances
  void AtmosphereImpedanceCondition::
    ComputeExactImpedance(int Lmax, const Real_wp& rt, const Complex_wp& k0,
        VectComplex_wp& impedance)
    {
      Globatto<Real_wp> lob;
      int p = 20;
      lob.ConstructQuadrature(20, lob.QUADRATURE_LOBATTO);
      lob.ComputeGradPhi(1e3*epsilon_machine);

      const Matrix<Real_wp>& dPhi = lob.GradPhi();

      impedance.Reallocate(Lmax+1);
      for (int l = 0; l <= Lmax; l++)
      {
        Complex_wp kinf2 = k0*k0 - l*(l+1) / (rt*rt);
        Complex_wp det = sqrt(Real_wp(4)*kinf2 - alpha*alpha);
        Complex_wp km = Real_wp(0.5)*(Iwp*alpha - det);
        Complex_wp kp = Real_wp(0.5)*(Iwp*alpha + det);

        Real_wp eta = realpart(-Iwp*(kp-km)); //DISP(kp); DISP(km); DISP(eta);
        //DISP(l); DISP(alpha); DISP(kinf2); DISP(det); DISP(km); DISP(kp); DISP(eta);

        Real_wp zmax = -log(threshold) / eta;
        //DISP(rt); DISP(zmax);

        Real_wp L = 2*pi_wp / abs(kp);
        int N = toInteger(ceil(Real_wp(0.5)*zmax / L));
        int nb_dof = N*p + 1;
        if (order > 0)
        {
          N = order;
          nb_dof = N*p + 1;
        }
        DISP(zmax); DISP(N);

        Real_wp dz = zmax / N;
        //Complex_wp beta(0, 0);
        Complex_wp beta = Iwp*kp;
        Matrix<Complex_wp, General, BandedCol> mat_lu;
        mat_lu.Reallocate(nb_dof, nb_dof, p, p);
        Matrix<Complex_wp> mat_elem(p+1, p+1);
        VectReal_wp Points(nb_dof);
        for (int i = 0; i < N; i++)
        {
          Real_wp r0 = rt + i*dz;
          Real_wp r1 = rt + (i+1)*dz;

          mat_elem.Zero();
          for (int k = 0; k <= p; k++)
          {
            Real_wp r = r0 + (r1-r0) * lob.Points(k);
            Points(i*p + k) = r;
            Real_wp poids = (r1-r0)*lob.Weights(k);

            // part - r^2 v'' + 2 r v'
            Real_wp coef_r = r*r*lob.Weights(k) / (r1-r0);
            for (int q = 0; q <= p; q++)
              for (int j = 0; j <= p; j++)
                mat_elem(q, j) += dPhi(q, k)*dPhi(j, k)*coef_r;

            // part - (2 beta r^2 + alpha r^2) v'
            Complex_wp coef_c = Real_wp(2)*beta*r*r + alpha*r*r;
            coef_c *= lob.Weights(k);
            for (int j = 0; j <= p; j++)
              mat_elem(k, j) -= dPhi(j, k)*coef_c;

            // part - [beta^2 r^2 + 2 beta r + \alpha \beta r^2 + k^2 r^2 - l (l+1) ] v
            coef_c = beta*(beta + alpha)*r*r + Real_wp(2)*beta*r + kinf2*r*r;
            coef_c *= poids;
            mat_elem(k, k) -= coef_c;
          }

          for (int q = 0; q <= p; q++)
            for (int j = 0; j <= p; j++)
              mat_lu.AddInteraction(i*p + q, i*p + j, mat_elem(q, j));

        }

        mat_lu.ClearRow(0);
        mat_lu.Set(0, 0, Complex_wp(1, 0));
        //mat_lu.WriteText("test_impedance.dat");

        VectComplex_wp rhs(nb_dof);
        rhs.Zero();
        rhs(0) = Complex_wp(1, 0);

        Vector<int> pivot;
        GetLU(mat_lu, pivot);
        SolveLU(mat_lu, pivot, rhs);

        VectReal_wp rhs_real(nb_dof), rhs_imag(nb_dof);
        for (int i = 0; i < nb_dof; i++)
        {
          rhs_real(i) = real(rhs(i));
          rhs_imag(i) = imag(rhs(i));
        }

        //Points.WriteText("points.dat");
        //rhs_real.WriteText("sol_real.dat");
        //rhs_imag.WriteText("sol_imag.dat");

        Complex_wp dphi_dn(0, 0);
        for (int i = 0; i <= p; i++)
          dphi_dn += dPhi(i, 0)*rhs(i);

        dphi_dn /= dz;
        dphi_dn += beta;
        //DISP(dphi_dn);
        impedance(l) = dphi_dn;
        // int test_input; cout << "waiting" << endl; cin >> test_input;
      }

      //DISP(impedance);
    }


  /***********************
   * VarHelmholtz_Radial *
   ***********************/


  //! default constructor
  template<class Complexe>
    VarHelmholtz_Radial<Complexe>::VarHelmholtz_Radial()
    {
      Lmax = 0;
      number_mode_to_be_computed = false; 
      this->type_coordinate = this->SPHERICAL;
      type_source = SRC_DIFFRACTED_FIELD;
      wave_vector.Init(0, 0, 1);
      diagonal_scaling = false;
      radius_gaussian_src = 0.0;
      radius_gaussian_cut_off = 0.0;
      display_du_dr = false;
    }  


  //! reads a line of the data file
  template<class Complexe>
    void VarHelmholtz_Radial<Complexe>
    ::SetInputData(const string& description_field, const VectString& parameters)
    {
      if (description_field == "NumberModes")
      {
        if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of VarHelmholtz_Radial" << endl;
          cout << "NumberModes needs more parameters, for instance :" << endl;
          cout << "NumberModes = n" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

        if (parameters(0) == "AUTO")
        {
          if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of VarHelmholtz_Radial" << endl;
            cout << "NumberModes needs more parameters, for instance :" << endl;
            cout << "NumberModes = AUTO n threshold" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

          number_mode_to_be_computed = true;
          Lmax = to_num<int>(parameters(1));
          threshold_mode = to_num<Real_wp>(parameters(2));
        }
        else
        {
          number_mode_to_be_computed = false;
          Lmax = to_num<int>(parameters(0));
        }
      }
      else if(description_field == "OrderQuadratureTheta")
      {
        int r = to_num<int>(parameters(0));
        gauss_teta.ConstructQuadrature(r);
      }
      else if ((description_field == "SismoGrille")
          || (description_field == "SismoPlaneAxi")
          || (description_field == "SismoPointsFileAxi")
          || (description_field == "SismoLineAxi")
          || (description_field == "SismoPointAxi")
          || (description_field == "SismoCircleAxi")
          || (description_field == "SismoGrille3D"))
      {
        // 3-D output on three planes
        int i = var_grid.GetM();
        var_grid.Resize(var_grid.GetM()+1);
        var_grid(i).SetInputData(description_field, parameters);
      }
      else if ((!description_field.compare("FileOutputGrille")) 
          || (!description_field.compare("FileOutputGrille3D")) ||
          (!description_field.compare("FileOutputPlaneAxi")) ||
          (!description_field.compare("FileOutputLine")) 
          || (!description_field.compare("FileOutputLineAxi")) || 
          (!description_field.compare("FileOutputPoint")) 
          || (!description_field.compare("FileOutputPointAxi")) || 
          (!description_field.compare("FileOutputCircleAxi")) ||
          (!description_field.compare("FileOutputPointsFileAxi")))
      {
        if (parameters.GetM() <= 1)
        {
          cout << "In SetInputData of VarHelmholtzRadial" << endl;
          cout << "FileOutputXXX needs more parameters, for instance :" << endl;
          cout << "FileOutputPlane = file_diffrac file_total component coef" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

        int type = GridInterpolationFull<Dimension1>::GetType(description_field);
        int nb = GetIndexOutputFiles(type);

        // Files of storage of solution on a regular grid. Diffracted and total field
        string ext = ".dat";

        this->output_grid_param(nb).SetDiffractedFieldFile(GetBaseString(parameters(0)) + ext);
        this->output_grid_param(nb).SetTotalFieldFile(GetBaseString(parameters(1)) + ext);
      }
      else if (description_field == "TypeSource")
      {
        if (parameters(0) == "SRC_DIFFRACTED_FIELD")
          type_source = SRC_DIFFRACTED_FIELD;
        else if (parameters(0) == "SRC_TOTAL_FIELD")
          type_source = SRC_TOTAL_FIELD;
        else if (parameters(0) == "SRC_DIRAC")
        {
          if (parameters.GetM() <= 3)
          {
            cout << "In SetInputData of VarHelmholtz_Radial" << endl;
            cout << "TypeSource needs more parameters, for instance :" << endl;
            cout << "TypeSource = SRC_DIRAC x0 y0 z0" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();                
          }

          type_source = SRC_DIRAC;
          Real_wp x0 = to_num<Real_wp>(parameters(1));
          Real_wp y0 = to_num<Real_wp>(parameters(2));
          Real_wp z0 = to_num<Real_wp>(parameters(3));
          origin_dirac.Init(x0, y0, z0);
        }
        else if ((parameters(0) == "SRC_GAUSSIAN") || (parameters(0) == "SRC_GAUSSIAN_R"))
        {
          if (parameters.GetM() <= 4)
          {
            cout << "In SetInputData of VarHelmholtz_Radial" << endl;
            cout << "TypeSource needs more parameters, for instance :" << endl;
            cout << "TypeSource = SRC_GAUSSIAN x0 y0 z0 r0" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();                
          }

          if (parameters(0) == "SRC_GAUSSIAN")
            type_source = SRC_GAUSSIAN;
          else
            type_source = SRC_GAUSSIAN_R;

          Real_wp x0 = to_num<Real_wp>(parameters(1));
          Real_wp y0 = to_num<Real_wp>(parameters(2));
          Real_wp z0 = to_num<Real_wp>(parameters(3));
          origin_dirac.Init(x0, y0, z0);
          radius_gaussian_src = to_num<Real_wp>(parameters(4));
          if (parameters.GetM() > 6)
          {
            radius_gaussian_cut_off = to_num<Real_wp>(parameters(6));
          }
          else
            radius_gaussian_cut_off = radius_gaussian_src;
        }
      }
      else if (!description_field.compare("IncidentAngle"))
      {
        if (parameters.GetM() <= 1)
        {
          cout << "In SetInputData of VarHelmholtz_Radial" << endl;
          cout << "IncidentAngle needs more parameters, for instance :" << endl;
          cout << "IncidentAngle = teta" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

        // in 3-D, the wave vector is oriented along
        // (sin(teta)*cos(phi) ; sin(teta)*sin(phi) ; cos(teta))
        Real_wp teta = 0.0, phi = 0.0;
        teta = to_num<Real_wp>(parameters(0))*pi_wp/180;
        phi = to_num<Real_wp>(parameters(1))*pi_wp/180;

        SetIncidentAngle(Real_wp(1), wave_vector, teta, phi);
      }
      else if (!description_field.compare("ScalingMatrix"))
      {
        if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of VarHelmholtz_Radial" << endl;
          cout << "ScalingMatrix needs more parameters, for instance :" << endl;
          cout << "ScalingMatrix = YES" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

        if (!parameters(0).compare("YES"))
          diagonal_scaling = true;
        else
          diagonal_scaling = false;
      }
      else if (!description_field.compare("ElectricOrMagnetic"))
      {
        if (parameters(0) == "-1")
          display_du_dr = true;
        else
          display_du_dr = false;
      }
      else
        VarHelmholtz_1D<Complexe>::SetInputData(description_field, parameters);
    }


  //! reads a boundary condition
  template<class Complexe>
    int VarHelmholtz_Radial<Complexe>
    ::ReadBoundaryCondition(int& nb, const VectString& parameters, int side)
    {
      if (parameters(nb) == "ATMOSPHERE")
        return atmo_condition.ReadParameters(nb, parameters, side);

      return VarHelmholtz_1D<Complexe>::ReadBoundaryCondition(nb, parameters, side);
    }


  //! inits the computation before solving each mode
  template<class Complexe>
    void VarHelmholtz_Radial<Complexe>
    ::ConstructAll(const string& input_file, const string& name_elt, MontjoieTimer& var_chrono)
    {
#ifdef SELDON_WITH_MPI
      int rank_proc; MPI_Comm_rank(this->comm_group_mode, &rank_proc);
#else
      int rank_proc(0);
#endif

      var_chrono.SetMessage("InputFile", "read the input file");
      var_chrono.Start("InputFile");

      ReadInputFile(input_file, *this);
      if ((this->print_level >= 0) && (rank_proc == 0))
        cout << "Input file has been successfully read" << endl;

      var_chrono.Stop("InputFile");

      var_chrono.SetMessage("ConstructMesh", "construct the mesh");
      var_chrono.Start("ConstructMesh");

      // computation of finite element and mesh
      this->ComputeMeshAndFiniteElement(name_elt);

      var_chrono.Stop("ConstructMesh");

      var_chrono.SetMessage("ConstructMass", "construct the mass matrix");
      var_chrono.Start("ConstructMass");

      this->ComputePhysicalCoefficients();

      // retrieving rho0 and mu0
      Complexe v0;
      int ne = this->mesh.GetNbElt()-1;
      int ref = this->mesh.Element(ne).GetReference();    
      int order = this->GetReferenceElement(ne).GetNbPointsQuadratureInside();
      this->GetPhysicalCoefficient(this->rho0, this->mu0, this->sigma0, v0,
          ne, order+1, ref);        

      Real_wp rt = this->xmax_;
      Complexe rho = this->rho0 + Iwp*this->sigma0/this->omega;
      Complexe kinf = this->omega*sqrt(rho/this->mu0);

      if (this->condition_right_side == BoundaryConditionEnum::LINE_TRANSPARENT)
      {
        // for a transparent condition, we precompute ratios h_n' / h_n
        Real_wp kr = realpart(kinf*this->xmax_);
        VectComplex_wp Hn, dHn;
        ComputeDeriveSphericalHankel(Lmax+2, kr, Hn, dHn);

        // then forming the coefficient -mu r^2 h_n' / h_n
        coef_impedance_transparent.Reallocate(Lmax+1);
        Complex_wp coef = this->mu0*square(this->xmax_)*kinf;
        for (int l = 0; l <= Lmax; l++)
          coef_impedance_transparent(l) = -coef*dHn(l)/Hn(l);
      }

      atmo_condition.UpdateImpedance(rt, kinf, this->beta_impedance_right);
      if (atmo_condition.ExactCondition())
      {	
        atmo_condition.ComputeExactImpedance(Lmax, rt, kinf,
            coef_impedance_transparent);

        DISP(coef_impedance_transparent);
        Complex_wp coef = this->mu0*square(this->xmax_);
        coef_impedance_transparent *= -coef;
      }

      // computing geometry quantities
      this->ComputeMassMatrix();

      int N = this->GetNbDof();
      int r = this->mesh.GetOrder();
      if ((this->print_level >= 0) && (rank_proc == 0))
      {
        cout << "Number of degrees of freedom = " << N << endl;
        cout << "Mesh size : " << this->mesh.GetMeshSize() << endl;
      }

      var_chrono.Stop("ConstructMass");

      var_chrono.SetMessage("ComputeFEM", "compute finite element matrices");
      var_chrono.Start("ComputeFEM");

      // computation of finite element matrices
      GlobalGenericMatrix<Complexe> nat_mat;

      stiffness_matrix.Reallocate(N, N, r, r);
      this->AddMatrixFEM(stiffness_matrix, nat_mat);
      this->AddBoundaryTerms(stiffness_matrix, nat_mat);

      if (this->print_level >= 4)
        stiffness_matrix.WriteText("stiffness.dat");

      // checking if mass lumped elements are used
      for (int i = 0; i < this->mesh.GetNbElt(); i++)
        if (!this->GetReferenceElement(i).LumpedMassMatrix())
        {
          cout << "In HelmholtzRadial, the mass matrix is "
            << "assumed to be diagonal" << endl;

          cout << "Use EDGE_LOBATTO element to achieve mass lumping" << endl;

          abort();
        }

      // diagonal part due to l(l+1) part
      mass_matrix.Reallocate(N);
      mass_matrix.Fill(0);
      for (int i = 0; i < this->mesh.GetNbElt(); i++)
        if (this->mat_mass_mu.GetM() > 0)
          for (int j = 0; j < this->mat_mass_mu(i).GetM(); j++)
            mass_matrix(this->mesh.GetNumberDof(i, j)) += this->mat_mass_mu(i)(j);

      if (this->condition_left_side == BoundaryConditionEnum::LINE_DIRICHLET)
        SetComplexZero(mass_matrix(0));

      if (this->condition_right_side == BoundaryConditionEnum::LINE_DIRICHLET)
        SetComplexZero(mass_matrix(N-1));

      // intermediary arrays are cleared to save memory
      this->mat_mass_mu.Clear(); this->mat_mass.Clear();
      this->mat_massV.Clear(); this->mat_sigma.Clear(); this->mat_sigmaV.Clear();

      var_chrono.Stop("ComputeFEM");

      // spherical harmonics
      fct_leg.Init(Lmax);
    }


  //! computes the display grid
  template<class Complexe>
    void VarHelmholtz_Radial<Complexe>::ComputeVarGrid()
    {
      // boundaries of the grid are set to boundaries of the domain if required
      for (int i = 0; i < var_grid.GetM(); i++)
        if (var_grid(i).GridIntervalToBeComputed())
        {
          var_grid(i).SetXmin0(this->xmin_);
          var_grid(i).SetXmax0(this->xmax_);
        }

      // components of unknowns to store in output files 
      // are selected with the global variable choice_field_output
      // this->choice_field_output
      for (int i = 0; i < this->output_grid_param.GetM(); i++)
        output_grid_param(i).SetComponent(0);

      for (int i = 0; i < var_grid.GetM(); i++)
        var_grid(i).InitGrid(this->var_section);

      // localization of the points of the predefined grids inside the mesh
      // implementation of this method in file GridInterpolation.cxx
      this->var_section.LocalizePoints(this->mesh);    
    }


  //! computes and factorizes finite element matrix for mode l
  template<class Complexe>
    void VarHelmholtz_Radial<Complexe>
    ::PerformFactorizationStep(Matrix<Complex_wp, General, BandedCol>& mat_lu,
        IVect& pivot, int l)
    {
      mat_lu = stiffness_matrix;
      //mat_lu.WriteText("MatBefore.dat");

      Real_wp coef = l*(l+1);
      for (int i = 0; i < this->mesh.GetNbDof(); i++)
        mat_lu.AddInteraction(i, i, mass_matrix(i)*coef);

      // part due to atmosphere condition
      if (atmo_condition.LaplacianPresent())
      {
        int irow = this->mesh.GetNbDof()-1;
        Real_wp rt = this->xmax_;
        Complexe rho = this->rho0 + Iwp*this->sigma0/this->omega;
        Complexe kinf = this->omega*sqrt(rho/this->mu0);

        Complex_wp mu, sigma, v0;
        int ref = this->mesh.Element(this->mesh.GetNbElt()-1).GetReference();
        int order = this->GetReferenceElement(this->mesh.GetNbElt()-1).GetNbPointsQuadratureInside();
        this->GetPhysicalCoefficient(rho, mu, sigma, v0, this->mesh.GetNbElt()-1, order+1, ref);

        Complexe impedance = atmo_condition.GetLaplacianImpedance(coef, this->xmax_, kinf);
        Complexe vloc = impedance*mu*rt*rt;
        mat_lu.AddInteraction(irow, irow, vloc);	
      }

      //DISP(coef_impedance_transparent);
      // part due to transparent condition
      int N = this->mesh.GetNbDof();
      if ((this->condition_right_side == BoundaryConditionEnum::LINE_TRANSPARENT)
          || (atmo_condition.ExactCondition()))
        mat_lu.AddInteraction(N-1, N-1, coef_impedance_transparent(l));

      // try with Dirichlet
      IVect last_row(1);
      last_row(0) = N-1;
      if (this->condition_right_side == BoundaryConditionEnum::LINE_DIRICHLET)
      {
        EraseCol(last_row, mat_lu);
        mat_lu.Set(N-1, N-1, Complex_wp(1, 0));
      }

      // scaling
      if (diagonal_scaling)
      {
        GetRowColSum(diagonal_scale_left, diagonal_scale_right, mat_lu);
        diagonal_scale_right.Fill(1.0);// = diagonal_scale_left;
        int n = mat_lu.GetM();
        for (int i = 0; i < n; i++)
        {
          if (diagonal_scale_left(i) != 0)
            diagonal_scale_left(i) = abs(1.0/diagonal_scale_left(i));
          else
            diagonal_scale_left(i) = 1.0;

          if (diagonal_scale_right(i) != 0)
            diagonal_scale_right(i) = sqrt(1.0/diagonal_scale_right(i));
          else
            diagonal_scale_right(i) = 1.0;
        }

        ScaleMatrix(mat_lu, diagonal_scale_left, diagonal_scale_right);
      }

      if (this->print_level >= 8)
        mat_lu.WriteText("MatAfter.dat");

      GetLU(mat_lu, pivot);
    }


  //! solves a linear system with mat_lu
  template<class Complexe>
    void VarHelmholtz_Radial<Complexe>
    ::ComputeSolution(Matrix<Complex_wp, General, BandedCol>& mat_lu,
        IVect& pivot, Vector<Complex_wp>& x)
    {  
      if (diagonal_scaling)
        for (int i = 0; i < x.GetM(); i++)
          x(i) *= diagonal_scale_left(i);

      SolveLU(mat_lu, pivot, x);
      //x.Write("sol.dat");

      if (diagonal_scaling)
        for (int i = 0; i < x.GetM(); i++)
          x(i) *= diagonal_scale_right(i);    
    }


  //! computes right hand side for all modes l, m
  template<class Complexe>
    void VarHelmholtz_Radial<Complexe>
    ::ComputeRightHandSide(Vector<Vector<Vector<Complexe> > >& rhs, MontjoieTimer& var_chrono)
    {
      rhs.Reallocate(Lmax+1);
      for (int l = 0; l <= Lmax; l++)
        rhs(l).Reallocate(2*l + 1);

      if (this->type_source == SRC_DIFFRACTED_FIELD)
        ComputeRhsDiffractedField(rhs);
      else if (this->type_source == SRC_TOTAL_FIELD)
        ComputeRhsTotalField(rhs);
      else if (this->type_source == SRC_GAUSSIAN_R)
      {
        // allocating the right hand side
        for (int l = 0; l <= Lmax; l++)
          for (int m = 0; m < rhs(l).GetM(); m++)
          {
            rhs(l)(m).Reallocate(this->GetNbDof());
            rhs(l)(m).Fill(0);
          }

        VectReal_wp s(2), points;
        VectComplex_wp feval;
        VectComplex_wp contrib;

        Real_wp z0 = origin_dirac(2);
        Real_wp RadiusSource2 = square(radius_gaussian_src);
        Real_wp alpha = -log(1.e-6)/RadiusSource2;
        Real_wp beta = sqrt(alpha/pi_wp);

        for (int i = 0; i < this->mesh.GetNbElt(); i++)
        {
          // we get vertices of the element i
          this->mesh.GetVerticesElement(i, s);
          const ElementReference<Dimension1, 1>& Fb = this->GetReferenceElementH1(i);
          // number of integration points
          int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
          int nb_dof_elt = Fb.GetNbDof();
          Real_wp h = abs(s(1) - s(0));
          int ref = this->mesh.Element(i).GetReference();

          Fb.FjElem(s, points);

          feval.Reallocate(nb_points_quadrature);
          feval.Fill(0);

          contrib.Reallocate(nb_dof_elt);

          for (int j = 0; j < nb_points_quadrature; j++)
          {
            Real_wp r_ = points(j);
            Complex_wp density(1, 0);
            if (this->helio_scaling_model)
            {
              if (this->log_rho_read(ref))
                density = exp(this->ref_rho(ref).GetCoefficient(*this, i, j+1));
              else
                density = this->ref_true_rho(ref).GetCoefficient(*this, i, j+1);
            }

            if (abs(r_ - z0) <= radius_gaussian_cut_off)
            {
              Complex_wp coef = density;
              coef *= h*r_*r_;

              Real_wp f = beta*exp(-alpha*(square(r_-z0)));
              feval(j) = coef*f;                    
            }
          }

          // we compute \int f \phi_i
          Fb.ComputeIntegralRef(feval, contrib);

          // we add contributions to b_source
          for (int j = 0; j < nb_dof_elt; j++)
          {
            int num_dof = this->mesh.GetNumberDof(i, j);
            rhs(0)(0)(num_dof) += contrib(j);
          }
        }

      }
      else if (this->type_source == SRC_GAUSSIAN)
      {
        if ((abs(origin_dirac(0)) <= 10.0*epsilon_machine) && (abs(origin_dirac(1)) <= 10.0*epsilon_machine))
        {
          // point on the axis Oz, only m = 0 is involved
          for (int l = 0; l <= Lmax; l++)
            rhs(l).Reallocate(1);
        }
        else
        {
          cout << "Not implemented" << endl;
          abort();
        }

        // allocating the right hand side
        for (int l = 0; l <= Lmax; l++)
          for (int m = 0; m < rhs(l).GetM(); m++)
          {
            rhs(l)(m).Reallocate(this->GetNbDof());
            rhs(l)(m).Fill(0);
          }

        VectReal_wp s(2), points;
        Vector<VectComplex_wp> feval(Lmax+1);
        VectComplex_wp contrib;

        Real_wp z0 = origin_dirac(2);
        Real_wp RadiusSource2 = square(radius_gaussian_src);
        Real_wp alpha = -log(1.e-6)/RadiusSource2;
        Real_wp beta = sqrt(alpha/pi_wp);

        var_chrono.SetMessage("ComputeRule", "compute quadrature rule");
        var_chrono.Start("ComputeRule");

        DISP(gauss_teta.GetOrder());
        if (gauss_teta.GetOrder() == 0)
          gauss_teta.ConstructQuadrature(2*Lmax+3);        

        var_chrono.Stop("ComputeRule");

        //var_chrono.SetMessage("QuadratureRhs", "compute rhs (quadrature)");
        //var_chrono.Start("QuadratureRhs");

        Vector<VectReal_wp> P;
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
        {
          // we get vertices of the element i
          this->mesh.GetVerticesElement(i, s);
          const ElementReference<Dimension1, 1>& Fb = this->GetReferenceElementH1(i);
          // number of integration points
          int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
          int nb_dof_elt = Fb.GetNbDof();
          Real_wp h = abs(s(1) - s(0));
          int ref = this->mesh.Element(i).GetReference();

          Fb.FjElem(s, points);

          for (int l = 0; l <= Lmax; l++)
          {
            feval(l).Reallocate(nb_points_quadrature);
            feval(l).Fill(0);
          }

          contrib.Reallocate(nb_dof_elt);

          for (int j = 0; j < nb_points_quadrature; j++)
          {
            Real_wp r_ = points(j);
            Complex_wp density(1, 0);
            if (this->helio_scaling_model)
            {
              if (this->log_rho_read(ref))
                density = exp(this->ref_rho(ref).GetCoefficient(*this, i, j+1));
              else
                density = this->ref_true_rho(ref).GetCoefficient(*this, i, j+1);
            }

            // DISP(density);		
            for (int k = 0; k <= gauss_teta.GetOrder(); k++)
            {
              Real_wp teta = pi_wp*gauss_teta.Points(k);
              Real_wp cos_teta = cos(teta);
              Real_wp sin_teta = sin(teta);
              Complex_wp coef = pi_wp*gauss_teta.Weights(k)*density;
              coef *= 2.0*pi_wp*sin_teta*h*r_*r_;

              Real_wp r2 = square(r_*sin_teta);
              Real_wp z = r_*cos_teta;
              Real_wp dist2 = square(z-z0) + r2;
              if (dist2 <= square(radius_gaussian_cut_off))
              {
                Complex_wp f = beta*exp(-alpha*dist2);
                fct_leg.EvaluatePnm(Lmax, 0, teta, P);
                for (int l = 0; l <= Lmax; l++)
                  feval(l)(j) += coef*P(l)(0)*f;   
              }                 		    
            }
          }

          // we compute \int f \phi_i
          for (int l = 0; l <= Lmax; l++)
          {
            Fb.ComputeIntegralRef(feval(l), contrib);

            // we add contributions to b_source
            for (int j = 0; j < nb_dof_elt; j++)
            {
              int num_dof = this->mesh.GetNumberDof(i, j);
              rhs(l)(0)(num_dof) += contrib(j);
            }
          }  
        }
        //var_chrono.Stop("QuadratureRhs");
      }
      else if (this->type_source == SRC_DIRAC)
      {
        Real_wp r0, theta0, phi0;
        CartesianToSpherical(origin_dirac(0), origin_dirac(1), origin_dirac(2),
            r0, theta0, phi0);

        // r0 is localized in the 1-D mesh
        int ne = 0; Real_wp pos_loc(0);
        while ((ne < this->mesh.GetNbVertices()) && (this->mesh.Vertex(ne) < r0))
          ne++;

        if (ne == 0)
        {
          if (abs(this->mesh.Vertex(0) - r0) <= 10*epsilon_machine)
            pos_loc = 0.0;
          else
          {
            cout << "Point " << r0 << " not found in the mesh " << endl;
            abort();
          }
        }
        else if (ne == this->mesh.GetNbVertices())
        {
          if (abs(this->mesh.Vertex(ne-1) - r0) <= 10*epsilon_machine)
          {
            ne = this->mesh.GetNbElt()-1;
            pos_loc = 1.0;
          }
          else
          {
            cout << "Point " << r0 << " not found in the mesh " << endl;
            abort();
          }
        }
        else
        {
          ne--;
          pos_loc = (r0 - this->mesh.Vertex(ne)) / (this->mesh.Vertex(ne+1) - this->mesh.Vertex(ne));
        }

        DISP(ne); DISP(pos_loc);

        int ref = this->mesh.Element(ne).GetReference();
        Complex_wp density(1, 0);
        if (this->helio_scaling_model)
        {
          if (this->log_rho_read(ref))
            density = exp(this->ref_rho(ref).EvaluateCoefficient(r0));
          else
            density = this->ref_true_rho(ref).EvaluateCoefficient(r0);

          DISP(density);
        }

        int Mmax = Lmax;
        if (abs(theta0) <= 10.0*epsilon_machine)
        {
          // point on the axis Oz, only m = 0 is involved
          Mmax = 0;
          for (int l = 0; l <= Lmax; l++)
            rhs(l).Reallocate(1);
        }

        // allocating the right hand side
        for (int l = 0; l <= Lmax; l++)
          for (int m = 0; m < rhs(l).GetM(); m++)
          {
            rhs(l)(m).Reallocate(this->GetNbDof());
            rhs(l)(m).Fill(0);
          }

        if (r0 == 0)
        {
          // point at the center only l = 0, m=0 is involved
          SetComplexReal(sqrt(1.0/(4.0*pi_wp)), rhs(0)(0)(0));
          rhs(0)(0)(0) *= density;
          return;
        }

        // then evaluating basis functions at r = r0
        const ElementReference<Dimension1, 1>& Fb = this->GetReferenceElementH1(ne);

        VectReal_wp phi;
        Fb.ComputeValuesPhiRef(pos_loc, phi);

        IVect num_dof(phi.GetM());
        for (int j = 0; j < num_dof.GetM(); j++)
          num_dof(j) = this->mesh.GetNumberDof(ne, j);

        // and spherical harmonics at theta0, phi0
        Vector<VectReal_wp> P;
        fct_leg.EvaluatePnm(Lmax, Mmax, theta0, P);

        Complex_wp phase_elt = Complex_wp(cos(phi0), sin(phi0)), phase;
        SetComplexOne(phase);
        for (int m = 0; m <= 2*Mmax; m++)
        {
          for (int l = (m+1)/2; l <= Lmax; l++)
          {
            for (int k = 0; k < num_dof.GetM(); k++)
              rhs(l)(m)(num_dof(k)) = density*phi(k)*GetYnm(l, m, P, phase);
          }

          if (m%2== 0)
            phase *= phase_elt;
        }        
      }

      if (this->condition_right_side == BoundaryConditionEnum::LINE_DIRICHLET)
        for (int l = 0; l < rhs.GetM(); l++)
          for (int m = 0; m < rhs(l).GetM(); m++)
            SetComplexZero(rhs(l)(m)(this->nodl-1));
    }


  //! computes decomposition of incident field on spherical harmonics
  template<class Complexe>
    void VarHelmholtz_Radial<Complexe>
    ::ComputeIncidentField(const Real_wp& kinf, const Real_wp& kr, Vector<VectComplex_wp>& u_inc,
        Vector<VectComplex_wp>& du_inc_dr, bool compute_grad)
    {
      u_inc.Reallocate(Lmax+1);
      if (compute_grad)
        du_inc_dr.Reallocate(Lmax+1);

      bool incidence_axial = true;
      Real_wp kx = wave_vector(0), ky = wave_vector(1), kz = wave_vector(2);
      Real_wp k_perp = sqrt(kx*kx + ky*ky), phi0(0);
      if (k_perp > 10.0*epsilon_machine)
      {
        incidence_axial = false; 
        phi0 = acos(kx/k_perp);
        if (ky < 0)
          phi0 = -phi0;
      }    

      VectReal_wp Jn, dJn;
      if (incidence_axial)
      {
        if (compute_grad)
          ComputeDeriveSphericalBessel(Lmax+1, kr, Jn, dJn);
        else
          ComputeSphericalBessel(Lmax+1, kr, Jn);

        for (int l = 0; l <= Lmax; l++)
        {
          u_inc(l).Reallocate(1);
          Complex_wp coef = sqrt(4.0*pi_wp*(2*l+1))*conj(ComputePowerI(l));
          if (compute_grad)
          {
            du_inc_dr(l).Reallocate(1);
            du_inc_dr(l)(0) = coef*dJn(l)*kinf;
          }

          u_inc(l)(0) = coef*Jn(l);
        }
      }
      else
      {
        if (gauss_teta.GetOrder() == 0)
          gauss_teta.ConstructQuadrature(2*Lmax+1);        

        Vector<Vector<Real_wp> > P;
        for (int l = 0; l <= Lmax; l++)
        {
          u_inc(l).Reallocate(2*l+1);
          du_inc_dr(l).Reallocate(2*l+1);            
          u_inc(l).Fill(0);
          du_inc_dr(l).Fill(0);
        }

        /* Vector<VectComplex_wp> u_test(Lmax+1);
           for (int l = 0; l <= Lmax; l++)
           {
           u_test(l).Reallocate(2*l+1);
           u_test(l).Fill(0);
           }

           DISP(kx); DISP(ky); DISP(kz); DISP(kr);
           for (int k = 0; k <= gauss_teta.GetOrder(); k++)
           for (int p = 0; p <= gauss_teta.GetOrder(); p++)
           {
           Real_wp teta = pi_wp*gauss_teta.Points(k);
           Real_wp phi = 2.0*pi_wp*gauss_teta.Points(p);
           Real_wp arg = kx*kr*sin(teta)*cos(phi) + ky*kr*sin(teta)*sin(phi) + kz*kr*cos(teta);
           Complex_wp phase = exp(Iwp*arg);
           Complex_wp poids = sin(teta)*pi_wp*gauss_teta.Weights(k)*2.0*pi_wp*gauss_teta.Weights(p);
           poids *= phase;
           fct_leg.EvaluatePnm(Lmax, Lmax, teta, P);

           Complex_wp phase_elt;
           phase_elt = Complex_wp(cos(phi), sin(phi));
           SetComplexOne(phase);

           for (int m = 0; m <= 2*Lmax; m++)
           {
           for (int l = (m+1)/2; l <= Lmax; l++)
           u_test(l)(m) += conj(GetYnm(l, m, P, phase))*poids;

           if (m%2 == 0)
           phase *= phase_elt;
           }
           } 

           for (int l = 0; l <= Lmax; l++)
           for (int m = 0; m <= 2*l; m++)
           {
           DISP(l); DISP(m); DISP(u_test(l)(m));
           }
           */

        Complex_wp phase_elt = Iwp*Complex_wp(cos(phi0), -sin(phi0)), phase;
        SetComplexOne(phase);
        for (int k = 0; k <= gauss_teta.GetOrder(); k++)
        {
          Real_wp teta = pi_wp*gauss_teta.Points(k);
          Real_wp cos_teta = cos(teta);
          Real_wp sin_teta = sin(teta);
          Real_wp arg = k_perp*kr*sin_teta;
          if (compute_grad)
            ComputeDeriveBesselFunctions(0, Lmax+1, arg, Jn, dJn);
          else
            ComputeBesselFunctions(0, Lmax+1, arg, Jn);

          SetComplexOne(phase);
          arg = kz*kr*cos_teta;
          Complex_wp coef = Complex_wp(cos(arg), sin(arg))*pi_wp*gauss_teta.Weights(k);
          coef *= 2.0*pi_wp*sin_teta;

          fct_leg.EvaluatePnm(Lmax, Lmax, teta, P);
          for (int m = 0; m <= 2*Lmax; m++)
          {
            int m2 = (m+1)/2;
            Complex_wp val = coef*Jn(m2);
            Complex_wp dval_dr = (Iwp*kz*kinf*cos_teta*Jn(m2) + k_perp*kinf*sin_teta*dJn(m2))*coef;
            for (int l = m2; l <= Lmax; l++)              
            {
              if (m%2 == 1)
              {
                u_inc(l)(m) += val*P(l)(m2)*conj(phase);
                du_inc_dr(l)(m) += dval_dr*P(l)(m2)*conj(phase);
              }
              else
              {
                u_inc(l)(m) += val*P(l)(m2)*phase;
                du_inc_dr(l)(m) += dval_dr*P(l)(m2)*phase;
              }
            }

            if (m%2== 0)
              phase *= phase_elt;
          }
        }

        /*for (int l = 0; l <= Lmax; l++)
          for (int m = 0; m <= 2*l; m++)
          {
          DISP(l); DISP(m); DISP(u_inc(l)(m));
          if (abs(u_inc(l)(m) - u_test(l)(m)) > 1e-12)
          {
          cout << "Difference ici" << endl;
          }
          }

          exit(0); */
      }
    }


  //! computes right hand side for diffracted field source
  template<class Complexe>
    void VarHelmholtz_Radial<Complexe>
    ::ComputeRhsDiffractedField(Vector<Vector<Vector<Complexe> > >& rhs)
    {
      Real_wp kinf = this->omega*sqrt(abs(this->rho0/this->mu0));
      Real_wp kr = kinf*this->xmin_;
      int N = this->mesh.GetNbDof();
      Vector<VectComplex_wp> u_inc, du_inc_dr;
      ComputeIncidentField(kinf, kr, u_inc, du_inc_dr, true);

      // for axial incidence rhs is resized
      for (int l = 0; l <= Lmax; l++)
        if (u_inc(l).GetM() == 1)
          rhs(l).Reallocate(1);

      for (int l = 0; l <= Lmax; l++)
        for (int m = 0; m < rhs(l).GetM(); m++)
        {
          rhs(l)(m).Reallocate(N);
          rhs(l)(m).Fill(0);
        }

      // part due to the boundary condition on the left side x = x0
      switch (this->condition_left_side)
      {
        case BoundaryConditionEnum::LINE_DIRICHLET :
          {
            for (int l = 0; l <= Lmax; l++)
              for (int m = 0; m < rhs(l).GetM(); m++)
                rhs(l)(m)(0) = -u_inc(l)(m);
          }
          break;
        case BoundaryConditionEnum::LINE_NEUMANN :
          {
            Real_wp r2 = square(this->xmin_);
            for (int l = 0; l <= Lmax; l++)
              for (int m = 0; m < rhs(l).GetM(); m++)
                rhs(l)(m)(0) = r2*this->mu0*du_inc_dr(l)(m);
          }
          break;
        case BoundaryConditionEnum::LINE_IMPEDANCE :
          {
            Real_wp r2 = square(this->xmin_);
            int ref = this->mesh.Element(0).GetReference();
            Complex_wp mu = this->ref_mu(ref).GetConstant();
            for (int l = 0; l <= Lmax; l++)
              for (int m = 0; m < rhs(l).GetM(); m++)
                rhs(l)(m)(0) = r2*this->mu0*du_inc_dr(l)(m) - r2*mu*this->beta_impedance_left*u_inc(l)(m);
          }
          break;
      }  

      // part due to dieletric layers
      Vector<Vector<VectComplex_wp> > feval(Lmax+1), feval_grad(Lmax+1);
      for (int l = 0; l <= Lmax; l++)
      {
        feval(l).Reallocate(rhs(l).GetM());
        feval_grad(l).Reallocate(rhs(l).GetM());
      }

      for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        int ref = this->mesh.Element(i).GetReference();
        Complex_wp rho = this->ref_rho(ref).GetConstant() + Iwp*this->ref_sigma(ref).GetConstant()/this->omega;
        Complex_wp mu = this->ref_mu(ref).GetConstant();
        VectReal_wp s(2);
        if ((rho != this->rho0) || (mu != this->mu0))
        {
          const ElementReference<Dimension1, 1>& Fb = this->GetReferenceElementH1(i);
          int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
          for (int l = 0; l <= Lmax; l++)
            for (int m = 0; m < rhs(l).GetM(); m++)
            {
              feval(l)(m).Reallocate(nb_points_quadrature);
              if (mu != this->mu0)
                feval_grad(l)(m).Reallocate(nb_points_quadrature);
            }

          // We get the two extremities of the edge
          VectComplex_wp contrib(Fb.GetNbDof());
          this->mesh.GetVerticesElement(i, s);

          // length of the edge
          Real_wp Ji = abs(s(1)-s(0));
          for (int k = 0; k < nb_points_quadrature; k++)
          {
            // quadrature point
            Real_wp x = Fb.Points(k)*s(1) + (1.0-Fb.Points(k))*s(0);
            Real_wp r2 = x*x;
            kr = kinf*x;
            Complex_wp coef_rho = this->omega2*r2*(rho-this->rho0)*Ji;
            Complex_wp coef_mu = (this->mu0 - mu);
            // computes incident field and its gradient
            ComputeIncidentField(kinf, kr, u_inc, du_inc_dr, true);

            for (int l = 0; l <= Lmax; l++)
              for (int m = 0; m < rhs(l).GetM(); m++)
              {
                if (rho != this->rho0)
                  feval(l)(m)(k) = coef_rho*u_inc(l)(m);
                else
                  SetComplexZero(feval(l)(m)(k));

                if (mu != this->mu0)
                {
                  Complex_wp coef_mul = coef_mu*Real_wp(l*(l+1))*Ji;
                  feval(l)(m)(k) += coef_mul*u_inc(l)(m);
                  feval_grad(l)(m)(k) = coef_mu*r2*du_inc_dr(l)(m);
                }
              }
          }


          for (int l = 0; l <= Lmax; l++)
            for (int m = 0; m < rhs(l).GetM(); m++)
            {
              Fb.ComputeIntegralRef(feval(l)(m), contrib);
              for (int j = 0; j < Fb.GetNbDof(); j++)
                rhs(l)(m)(this->mesh.GetNumberDof(i, j)) += contrib(j);

              if (mu != this->mu0)
              {
                Fb.ComputeIntegralGradientRef(feval_grad(l)(m), contrib);
                for (int j = 0; j < Fb.GetNbDof(); j++)
                  rhs(l)(m)(this->mesh.GetNumberDof(i, j)) += contrib(j);                    
              }
            }
        }
      }

      switch (this->condition_right_side)
      {
        case BoundaryConditionEnum::LINE_DIRICHLET :
        case BoundaryConditionEnum::LINE_NEUMANN :
        case BoundaryConditionEnum::LINE_IMPEDANCE :
          {
            cout << "not implemented" << endl;
            abort();
          }
      }  
    }


  //! computes right hand side for total field source
  template<class Complexe>
    void VarHelmholtz_Radial<Complexe>
    ::ComputeRhsTotalField(Vector<Vector<Vector<Complexe> > >& rhs)
    {
      Real_wp kinf = this->omega*sqrt(abs(this->rho0/this->mu0));
      Real_wp kr = kinf*this->xmax_;
      int N = this->mesh.GetNbDof();
      Vector<VectComplex_wp> u_inc, du_inc_dr;
      ComputeIncidentField(kinf, kr, u_inc, du_inc_dr, true);

      // for axial incidence rhs is resized
      for (int l = 0; l <= Lmax; l++)
        if (u_inc(l).GetM() == 1)
          rhs(l).Reallocate(1);

      for (int l = 0; l <= Lmax; l++)
        for (int m = 0; m < rhs(l).GetM(); m++)
        {
          rhs(l)(m).Reallocate(N);
          rhs(l)(m).Fill(0);
        }

      switch (this->condition_right_side)
      {
        case BoundaryConditionEnum::LINE_ABSORBING :
          {
            Real_wp r2 = square(this->xmax_);
            for (int l = 0; l <= Lmax; l++)
              for (int m = 0; m < rhs(l).GetM(); m++)
                rhs(l)(m)(N-1) = this->mu0*r2*(du_inc_dr(l)(m) -Iwp*kinf*u_inc(l)(m));          
          }
          break;
      }    

      if (this->condition_left_side == BoundaryConditionEnum::LINE_ABSORBING)
      {
        cout << "not implemented" << endl;
        abort();
      }
    }


  //! write solution on required outputs
  template<class Complexe>
    void VarHelmholtz_Radial<Complexe>
    ::WriteDatas(Vector<Vector<Vector<Complexe> > >& sol)
    {
      int L = 0;
      for (int l = 0; l < sol.GetM(); l++)
        if (sol(l).GetM() > 0)
          L = l;

      // counting M
      int M = 0;
      for (int l = 0; l <= L; l++)
        for (int m = 0; m < sol(l).GetM(); m++)
          if (sol(l)(m).GetM() > 0)
            M = max(m, M);

      Matrix<Complexe> u(L+1, M+1);
      VectReal_wp val_phi;
      Complex_wp phase, phase_elt;
      Vector<Vector<Real_wp> > P;
      for (int n = 0; n < var_grid.GetM(); n++)
      {
        const IVect& list_points = var_grid(n).GetPointNumber();
        Vector<Complexe> trace_u(list_points.GetM());
        for (int p = 0; p < list_points.GetM(); p++)
        {
          int num_point = list_points(p);
          int ne = this->var_section.GetElementNumber(num_point);
          Real_wp pt_loc = this->var_section.GetLocalCoordinate(num_point);
          Complex_wp val_u(0, 0);
          if (ne >= 0)
          {
            u.Fill(0);
            const ElementReference<Dimension1, 1>& Fb = this->GetReferenceElementH1(ne);
            int nb_dof_loc = Fb.GetNbDof();
            val_phi.Reallocate(nb_dof_loc);
            if (display_du_dr)
            {
              TinyVector<Real_wp, 2> s;
              Fb.ComputeGradientPhiRef(pt_loc, val_phi);
              this->mesh.GetVerticesElement(ne, s);
              Real_wp h = abs(s(1) - s(0));
              Mlt(1.0/h, val_phi);
            }
            else	  
              Fb.ComputeValuesPhiRef(pt_loc, val_phi);

            for (int l = 0; l <= L; l++)
              for (int m = 0; m < sol(l).GetM(); m++)
                if (sol(l)(m).GetM() > 0)
                  for (int j = 0; j < nb_dof_loc; j++)
                    u(l, m) += sol(l)(m)(this->mesh.GetNumberDof(ne, j)) * val_phi(j);

            Real_wp teta = this->var_section.TetaInterp(num_point);
            Real_wp phi = this->var_section.PhiInterp(num_point);
            fct_leg.EvaluatePnm(L, (M+1)/2, teta, P);
            phase_elt = Complex_wp(cos(phi), sin(phi));
            SetComplexOne(phase);
            for (int m = 0; m <= M; m++)
            {
              for (int l = (m+1)/2; l <= L; l++)
                if ((sol(l).GetM() > 0) && (sol(l)(m).GetM() > 0))
                  val_u += GetYnm(l, m, P, phase)*u(l, m);

              if (m%2 == 0)
                phase *= phase_elt;
            }     
          }            
          trace_u(p) = val_u;
        }

        int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
        bool ascii = false;

        string file_name_diff = output_grid_param(n).GetDiffractedFieldFile();
        string file_name_tot = output_grid_param(n).GetTotalFieldFile();

        if ((type_source == SRC_DIFFRACTED_FIELD) ||
            (type_source == SRC_TOTAL_FIELD))
        {
          if (type_source == SRC_DIFFRACTED_FIELD)
            WriteMatlab(trace_u, var_grid(n), file_name_diff, double_prec, ascii);
          else
            WriteMatlab(trace_u, var_grid(n), file_name_tot, double_prec, ascii);

          // computing the total field or diffracted field by computing the incident field
          Real_wp kinf = this->omega*sqrt(abs(this->rho0/this->mu0));
          Real_wp kx = wave_vector(0), ky = wave_vector(1), kz = wave_vector(2);
          Complex_wp u_inc;
          for (int p = 0; p < list_points.GetM(); p++)
          {
            int num_point = list_points(p);
            int ne = this->var_section.GetElementNumber(num_point);
            if (ne >= 0)
            {
              Real_wp r = this->var_section.GetGlobalCoordinate(num_point);
              Real_wp teta = this->var_section.TetaInterp(num_point);
              Real_wp phi = this->var_section.PhiInterp(num_point);
              Real_wp x, y, z;
              SphericalToCartesian(r, teta, phi, x, y, z);
              u_inc = Complex_wp(cos(kinf*(kx*x + ky*y + kz*z)), sin(kinf*(kx*x + ky*y +kz*z)));
              if (type_source == SRC_DIFFRACTED_FIELD)
                trace_u(p) += u_inc;
              else
                trace_u(p) -= u_inc;
            }
          }

          if (type_source == SRC_DIFFRACTED_FIELD)
            WriteMatlab(trace_u, var_grid(n), file_name_tot, double_prec, ascii);
          else
            WriteMatlab(trace_u, var_grid(n), file_name_diff, double_prec, ascii);
        }
        else
          WriteMatlab(trace_u, var_grid(n), file_name_tot, double_prec, ascii);
      }
    }

}

#define MONTJOIE_FILE_HELMHOLTZ_RADIAL_CXX
#endif
