#ifndef MONTJOIE_FILE_IMPEDANCE_HELMHOLTZ_CXX

#include "ImpedanceHelmholtzInline.cxx"
#include "PhysicalConstants.hxx"

namespace Montjoie
{

  //! computes coefficients alpha, beta, delta, eta from parameters gamma, theta and zeta
  void ImpedanceABC_HelmholtzBase
    ::ComputeCoefficientsABC(Complex_wp myk, Real_wp courbure, int order, bool grazing_abc,
        Complex_wp& alpha_cla, Complex_wp& beta_cla,
        Complex_wp& delta_cla, Complex_wp& eta_cla,
        Real_wp gamma_cla_coef, Real_wp theta_cla_coef,
        Real_wp zeta_cla_coef, Real_wp h, Real_wp g)
    {
    //ik contains -i omega/c
    // alpha_cla is pre-filled with 1
    // beta_cla is pre-filled with ik (= -i omega/c)
    Complex_wp ik = Iwp * myk;
    Complex_wp courburecomp = courbure;
    //      Complex_wp SixKappaUnTiers = pow(6.0*courbure, 1.0/3.0); attention nan quand courbure est negative
    Complex_wp SixKappaUnTiers = pow(6.0*abs(courbure), Real_wp(1)/3)*exp(Iwp*arg(courburecomp)/3.0);
    // racine troisieme de l'unite : celle qui a la partie imaginaire négative.
    Complex_wp coeff= pow(myk*myk, Real_wp(1)/3)*exp(-Iwp*pi_wp/3.0)
      *SixKappaUnTiers*tgamma(2.0/3.0)/tgamma(1.0/3.0);
    Real_wp theta_cla = courbure*theta_cla_coef;
    Real_wp gamma_cla = courbure*gamma_cla_coef;
    Real_wp zeta_cla  = courbure*zeta_cla_coef ;
    // C-ABC
    //    DISP(order); DISP(gamma_cla_coef);
    switch(order)

      /*    
	    ::ComputeCoefficientsABC(Complex_wp myk, Real_wp courbure, int order, bool grazing_abc,
	    Complex_wp& alpha_cla, Complex_wp& beta_cla,
	    Complex_wp& delta_cla, Complex_wp& eta_cla,
	    Real_wp gamma_cla_coef, Real_wp theta_cla_coef,
	    Real_wp zeta_cla_coef, Real_wp h, Real_wp g)
	    {
	    Complex_wp ik = Iwp * myk;
	    Complex_wp courburecomp = courbure;
	    //      Complex_wp SixKappaUnTiers = pow(6.0*courbure, 1.0/3.0); attention nan quand courbure est negative
	    Complex_wp SixKappaUnTiers = pow(6.0*abs(courbure), Real_wp(1)/3)*exp(Iwp*arg(courburecomp)/3.0);
	    // racine troisieme de l'unite : celle qui a la partie imaginaire négative.
	    Complex_wp coeff= pow(myk*myk, Real_wp(1)/3)*exp(-Iwp*pi_wp/3.0)
	    *SixKappaUnTiers*tgamma(2.0/3.0)/tgamma(1.0/3.0);
	    Real_wp theta_cla = courbure*theta_cla_coef;
	    Real_wp gamma_cla = courbure*gamma_cla_coef;
	    Real_wp zeta_cla  = courbure*zeta_cla_coef ;
	    
	    // C-ABC
	    //DISP(order); DISP(1.0/gamma_cla_coef); DISP(courbure);
	    switch(order) 
      */
      {
        case CLA_1://if (order == 1) // equation 40 (Km1 omega 0) C-ABC
          beta_cla += courbure/2.0;
          //DISP("CLA_1");
          break;

        case CLA_2: // equation 43 (Km1 omega1 ) C-ABC + Laplace
          beta_cla += courbure/2.0;
          delta_cla = 1.0/(2.0*ik);
          break;

        case KM2_DELTA:// (Km2 delta )
          alpha_cla += (gamma_cla+theta_cla-courbure/4.0)/ik +
            (courbure*courbure/8.0 - (theta_cla+gamma_cla)*courbure/4.0+(theta_cla*gamma_cla))/(ik*ik) -
            ((gamma_cla-courbure/2.0)*(zeta_cla+courbure/4.0)*courbure/4.0)/(ik*ik*ik);

          beta_cla+= ( gamma_cla+theta_cla+courbure/4.0 ) +
            ( ( theta_cla+gamma_cla)*courbure/4.0 + gamma_cla*theta_cla - courbure*courbure/8.0 )/(ik)-
            (courbure/2.0-gamma_cla)*(zeta_cla-courbure/4.0)*courbure/(4.0*ik*ik);

          break;

        case KM2_OMEGA:// km2 omega 1


          alpha_cla += (gamma_cla+theta_cla-courbure/4.0)/ik;

          beta_cla+= ( gamma_cla+theta_cla+courbure/4.0 ) +
            ( ( theta_cla+gamma_cla)*courbure/4.0 + gamma_cla*theta_cla - courbure*courbure/8.0 )/(ik);

          delta_cla = 1.0/(2.0*ik);

          break;

        case KM2_OMEGA2:// km2 omega 2 symétrique :  a utiliser avec gamma = kappa/4


          alpha_cla += (gamma_cla+theta_cla-courbure/4.0)/ik +
            (courbure*courbure/8.0 - (theta_cla+gamma_cla)*courbure/4.0+(theta_cla*gamma_cla))/(ik*ik);

          beta_cla+= ( gamma_cla+theta_cla+courbure/4.0 ) +
            ( ( theta_cla+gamma_cla)*courbure/4.0 + gamma_cla*theta_cla - courbure*courbure/8.0 )/(ik)-
            (courbure/2.0-gamma_cla)*(zeta_cla-courbure/4.0)*courbure/(4.0*ik*ik);

          delta_cla = (1.0-courbure/(2.0*ik))/(2.0*ik);

          break;

        case KM2_DELTA2:// Km2 delta2 (non symétrique mais pour cercle ok)
          alpha_cla += (gamma_cla+theta_cla-courbure/4.0)/ik +
            (courbure*courbure/8.0 - (theta_cla+gamma_cla)*courbure/4.0+(theta_cla*gamma_cla))/(ik*ik)+
            -(gamma_cla-courbure/2.0)*(zeta_cla+courbure/4.0)*courbure/(4.0*ik*ik*ik);

          beta_cla+= ( gamma_cla+theta_cla+courbure/4.0 ) +
            ( ( theta_cla+gamma_cla)*courbure/4.0 + gamma_cla*theta_cla - courbure*courbure/8.0 )/(ik)+
            (gamma_cla-courbure/2.0)*(zeta_cla-courbure/4.0)*courbure/(4.0*ik*ik);

          delta_cla = (1.0-courbure/(2.0*ik))/(2.0*ik)
            -(3*courbure*(gamma_cla-3.0*courbure/2.0)/4.0 +theta_cla*(gamma_cla+3.0*courbure/4.0))/(2*ik*ik*ik)
            + (3.0*courbure*courbure*(gamma_cla-courbure)/8.0 - courbure*zeta_cla*(gamma_cla-5.0*courbure/4.0))/(2.0*ik*ik*ik*ik);
          // a verifier le eta
          eta_cla = -(theta_cla+gamma_cla-3.0*courbure/4.0)/(2.0*ik*ik*ik)
            + (courbure*(gamma_cla-5.0*courbure/4.0)/2.0-theta_cla*(gamma_cla-courbure/2.0))/(ik*ik*ik*ik)
            +(courbure*(courbure*(7.0*gamma_cla-13.0*courbure/2.0)/4.0+zeta_cla*(5.0*gamma_cla-11.0*courbure/2.0)))/(8.0*ik*ik*ik*ik*ik) ;

          break;
        case PARAM_KM1_DELTA://if (order == 11) // Km1 delta avec parametre
          beta_cla += (gamma_cla+courbure/4.0);
          alpha_cla += (gamma_cla-courbure/4.0)/ik;
          break;

        case PARAM_KM1_OMEGA: // Km1 omega avec parametre
          beta_cla += (gamma_cla+courbure/4.0);
          alpha_cla +=(gamma_cla-courbure/4.0)/ik;
          delta_cla = 1.0/(2.0*ik);

          break;

        case PARAM_KM0_OMEGA://if (order == 13) // Km0 omega avec parametre
          beta_cla += (gamma_cla+courbure/4.0);
          break;

        case ATMO_1:// condition d'atmosphere d'ordre 1
          courbure = courbure/2.0;
          alpha_cla = 1.0;
          beta_cla = 1/gamma_cla_coef;
          //	DISP(beta_cla);
          break;

        case ATMO_2: // condition d'atmosphere d'ordre 2
          courbure = courbure/2.0;
          alpha_cla = gamma_cla_coef*(1-2*gamma_cla_coef*courbure);
          DISP(alpha_cla);
          beta_cla = 1.0;
          break;

        case ATMO_3: // condition d'atmosphere d'ordre 3
          courbure = courbure/2.0;

          // Neumann to Dirichlet expression
          //alpha_cla = gamma_cla_coef*(1-2*gamma_cla_coef*courbure+square(gamma_cla_coef)*(6*square(courbure)+myk*myk));
          //beta_cla = 1.0;
          //eta_cla = pow(gamma_cla_coef,3);

          // Dirichlet to Neumann expression
          alpha_cla = 1.0;
          beta_cla = ((1+2*gamma_cla_coef*courbure) - 2.0*square(gamma_cla_coef*courbure)
              -square(gamma_cla_coef*myk))/gamma_cla_coef;

          delta_cla = -gamma_cla_coef;
          break;

        case ATMO_RBC_1: // condition AtmoABC // cas 1 du radial
          { // gamma_cla_coef was filled as 1/ given coef in .ini 
            // see AsiSymHelmholtz.cxx 
            // ik a le signe oppose de celui du radial
            //	  Complex_wp alpha_loc = Real_wp(0.5)/gamma_cla_coef;
            // Complex_wp kinf2 = myk*myk;
            //Complex_wp det = sqrt(kinf2 - alpha_loc*alpha_loc);
            //Complex_wp kp = Iwp*alpha_loc + det;
            //beta_cla = -Iwp*kp;
            //DISP(beta_cla);
            Complex_wp k0 = -myk;
            Complex_wp kinf2 = k0*k0;
            Complex_wp alpha = Real_wp(1)/gamma_cla_coef;
            alpha_cla = 1.0;
            beta_cla =  -( - alpha*Real_wp(0.5) + Iwp*k0*sqrt( Real_wp(1.0)-alpha*alpha*Real_wp(0.25)/kinf2 ) );


          }
          break;

        case ATMO_SAI_0: // condition AtmoABC avec courbure // cas 2 du radial
          { // gamma_cla_coef was filled as 1/ given coef in .ini 
            // see AsiSymHelmholtz.cxx 
            //	  Complex_wp alpha = Real_wp(0.5)/gamma_cla_coef + h; // h is the mean curvature
            //Complex_wp kinf2 = myk*myk;
            //Complex_wp det = sqrt(kinf2 - alpha*alpha);
            //Complex_wp kp = Iwp*alpha + det;

            //	  beta_cla = -Iwp*kp;
            Complex_wp k0 = -myk;
            Complex_wp kinf2 = k0*k0;
            Complex_wp alpha = Real_wp(1)/gamma_cla_coef;
            Real_wp kappa = h;
            alpha_cla = Real_wp(1.0);
            beta_cla = -( -(kappa + alpha*0.5) + Iwp*k0*sqrt(Real_wp(1)-(alpha*kappa+alpha*alpha*0.25)/kinf2) ) ;
          }
          break;

        case ATMO_SAI_1: // condition AtmoABC avec courbure et laplace beltrami  // cas 6 du radial
          { // gamma_cla_coef was filled as 1/ given coef in .ini 
            // see AsiSymHelmholtz.cxx 
            //	  Complex_wp alpha = Real_wp(0.5)/gamma_cla_coef + h; // h is the mean curvature
            //Complex_wp kinf2 = myk*myk;
            //Complex_wp det = sqrt(kinf2 - alpha*alpha);
            //Complex_wp kp = Iwp*alpha + det;
            //alpha_cla = Real_wp(1.0);
            //beta_cla = -Iwp*kp;
            //delta_cla = Real_wp(1.0)/(Real_wp(2.0)*ik);

            Complex_wp k0 = -myk;
            Complex_wp alpha = Real_wp(1)/gamma_cla_coef;
            Real_wp kappa = h;
            Complex_wp term_sqrt = sqrt(Real_wp(1)-( alpha*alpha*Real_wp(0.25) + alpha*kappa  )/(k0*k0) );
            Complex_wp Beltra = Real_wp(1)/(Real_wp(2)*Iwp*k0*term_sqrt) ;
            alpha_cla = Real_wp(1);
            beta_cla = -( -(kappa +Real_wp(0.5)*alpha ) +Iwp*k0*term_sqrt );
            delta_cla = - Beltra ; 
          }
          break;
        case ATMO_HF_1: // cas 3 du radial
          {
            Complex_wp k0 = -myk;
            Complex_wp alpha = Real_wp(1)/gamma_cla_coef;
            Real_wp kappa = h;
            alpha_cla = Real_wp(1.0);
            beta_cla = - ( -(kappa + Real_wp(0.5)*alpha) + Iwp*k0  +(  alpha*alpha*Real_wp(0.25) + alpha*kappa ) / (Real_wp(2)*Iwp*k0 ) );
            delta_cla = - Real_wp(1)/(Real_wp(2)*Iwp*k0);

          }
          break;
        case ATMO_ML_SPHERE:// microlocal works only for sphere // cas 5 du radial
          {
            Complex_wp k0 = -myk;
            Complex_wp alpha = Real_wp(1)/gamma_cla_coef;
            Real_wp kappa = h;
            Complex_wp Beltra =  ( Real_wp(1)+ kappa/(Iwp*k0) ) /(Real_wp(2) *Iwp* k0 );
            Complex_wp correc = alpha*(kappa + alpha*Real_wp(0.25) +kappa*(Real_wp(1.0)-Real_wp(2.0)*kappa)/(Real_wp(2.0)*Iwp*k0))/(Real_wp(2)*Iwp*k0) ;
            alpha_cla = Real_wp(1.0);
            beta_cla = - ( -(kappa+alpha*Real_wp(0.5)) +Iwp*k0 + correc) ;
            delta_cla = - Beltra ;
          }
          break;
        case ATMO_ML_ANY: // microlocal general for any shape // equiv cas 5 du radial // A DEBUGGUER
          {
            Complex_wp k0 = -myk;
            Complex_wp alpha = Real_wp(1.0)/gamma_cla_coef;
            Complex_wp correc = ( (g-h*h)*(Real_wp(1.0)+Real_wp(2.0)*h/(Iwp*k0)) + alpha*h + alpha*alpha*Real_wp(0.25) + alpha*h*(Real_wp(0.5)-h)/(Iwp*k0)) /(Real_wp(2.0)*Iwp*k0) ;
            alpha_cla = Real_wp(1.0);
            beta_cla = - (  - (h + alpha*Real_wp(0.5)  ) + Iwp*k0 + correc) ;
            delta_cla =- Real_wp(1.0)/(Real_wp(2.0)*Iwp*k0); // to be multiplied by the right value in AxiSymHelmholtz.cxx L1251

          }
          break;

        case WHIT_1: // NATHAN -- Condition from Ha's ABC radiation condition
          {
            //
            // DISP("CLA HA");abort();
            Real_wp alpha = Real_wp(1.0)/gamma_cla_coef;
            Complex_wp kh = sqrt(myk*myk-Real_wp(0.25)*alpha*alpha);
            courbure = Real_wp(0.5)*courbure;
            // Complex_wp sq_root = sqrt(Real_wp(1.0)-Real_wp(0.25)*alpha*alpha/(myk*myk));
            // Complex_wp factor = Real_wp(1.0)-Real_wp(0.5)*(alpha*courbure)/(myk*myk-Real_wp(0.25)*alpha*alpha);
            // beta_cla = courbure + Real_wp(0.5)*alpha - Iwp*myk*sq_root*factor;

            beta_cla = courbure + Real_wp(0.5)*alpha - Iwp*(kh-Real_wp(0.5)*courbure*alpha/kh);
            alpha_cla = Real_wp(1.0);
            eta_cla = Real_wp(0.0);
            delta_cla = Real_wp(0.0);
            // DISP(alpha_cla); DISP(beta_cla);
          }
          break;


        case WHIT_2: // NATHAN -- Condition from Ha's ABC radiation condition
          {
            Real_wp alpha = Real_wp(1.0)/gamma_cla_coef;
            Complex_wp kh = sqrt(myk*myk-Real_wp(0.25)*alpha*alpha);
            courbure = Real_wp(0.5)*courbure;
            // Complex_wp sq_root = sqrt(Real_wp(1.0)-Real_wp(0.25)*alpha*alpha/(myk*myk));
            // Complex_wp factor = Real_wp(1.0)-Real_wp(0.5)*(alpha*courbure)/(myk*myk-Real_wp(0.25)*alpha*alpha);
            // beta_cla = courbure + Real_wp(0.5)*alpha - Iwp*myk*sq_root*factor;

            Complex_wp cplx_part = kh-Real_wp(0.5)*courbure*alpha/kh;
            beta_cla = courbure + Real_wp(0.5)*alpha - Iwp*cplx_part;
            alpha_cla = Real_wp(1.0);
            eta_cla = Real_wp(0.0);
            delta_cla = Iwp*cplx_part*Real_wp(0.5)/(kh*kh);
            // DISP(alpha_cla); DISP(beta_cla);
          }
          break;
        case WHIT_3: // NATHAN -- Condition from Ha's ABC radiation condition
          {
            cout << "Order 3 not implemented yet" << endl;
            abort();
            Real_wp alpha = Real_wp(1.0)/gamma_cla_coef;
            Complex_wp kh = sqrt(myk*myk-Real_wp(0.25)*alpha*alpha);
            courbure = Real_wp(0.5)*courbure;
            // Complex_wp sq_root = sqrt(Real_wp(1.0)-Real_wp(0.25)*alpha*alpha/(myk*myk));
            // Complex_wp factor = Real_wp(1.0)-Real_wp(0.5)*(alpha*courbure)/(myk*myk-Real_wp(0.25)*alpha*alpha);
            // beta_cla = courbure + Real_wp(0.5)*alpha - Iwp*myk*sq_root*factor;

            Complex_wp cplx_part = kh-Real_wp(0.5)*courbure*alpha/kh;
            beta_cla = courbure + Real_wp(0.5)*alpha - Iwp*cplx_part;
            alpha_cla = Real_wp(1.0);
            eta_cla = Real_wp(0.0);
            delta_cla = Iwp*cplx_part*Real_wp(0.5)/(kh*kh);
            // DISP(alpha_cla); DISP(beta_cla);
          }
          break;

        case VISCOTHERMAL: // condition of Berggren 2018 for viscothermal losses at a pipe boundary
          { // gamma_cla_coef was filled as the temperature coef in .ini 
            // see AxiSymHelmholtz.cxx 
	          // quand je fais des calculs j'ai l'impression que : 
	          // alpha \partial_r p + beta u - \delta \Delta_T u = 0
	          // alpha \partial_r p + beta u + \delta \xi^2 u = 0
	          Real_wp temperature = gamma_cla_coef;
	          Physics physics(temperature);
	          //	    physics.Display();
            Complex_wp k0 = -myk;
	          Complex_wp omega = k0*physics.c;
            Complex_wp kinf2 = k0*k0;
	          Complex_wp sqrtminusi = (1-Iwp)/sqrt(2);
            alpha_cla = Real_wp(1.0);
	          Complex_wp thermo = (physics.gamma-1)*sqrt(physics.kappa / (omega * physics.rho*physics.Cp));
	          Complex_wp visco = sqrt(physics.mu/(physics.rho*omega));
	          //	    beta_cla = -kinf2*(thermo)*sqrtminusi;
	          //	    delta_cla = -visco*sqrtminusi;
	          beta_cla = (Iwp-1)*0.5*kinf2*(physics.gamma-1)*sqrt(2*physics.kappa/(omega*physics.rho*physics.Cp));
	          delta_cla = (Iwp-1)*0.5*sqrt(2*physics.mu/(omega*physics.rho));
	          

          }
          break;

        case CREMER: // condition of Cremer for viscothermal losses at a pipe boundary
          { // gamma_cla_coef was filled as the temperature coef in .ini 
            // see AsiSymHelmholtz.cxx 
	          Real_wp temperature = gamma_cla_coef;
	          Physics physics(temperature);
	          //	    physics.Display();
            Complex_wp k0 = -myk;
	          Complex_wp omega = k0*physics.c;
            Complex_wp kinf2 = k0*k0;
	          Complex_wp sqrtminusi = (1-Iwp)/sqrt(2);
            alpha_cla = Real_wp(1.0);
	          Complex_wp thermo = (physics.gamma-1)*sqrt(physics.kappa / (omega * physics.rho*physics.Cp));
	          Complex_wp visco = sqrt(physics.mu/(physics.rho*omega));
	          //beta_cla = -kinf2*(thermo + visco)*sqrtminusi;
	          beta_cla = (Iwp-1)*0.5*kinf2*((physics.gamma-1)*sqrt(2*physics.kappa/(omega*physics.rho*physics.Cp)) +sqrt(2*physics.mu/(omega*physics.rho)));
	          delta_cla = 0.0;
	    
          }
          break;

        default :
          {
            cout << "Order " << order << " not implemented" << endl;
            abort();
          }

      }

      if (grazing_abc)
      {
        // Ondes Rampantes Vero thèse p 136
        if (order==0)
        {
          beta_cla = coeff;
        }
        else
        {
          // (beta_cla_orig ) \p_r u + (\alpha_cla_orig + \delta_cla_orig \xi^2 ) u = 0
          Complex_wp alpha_cla_orig = alpha_cla;
          Complex_wp beta_cla_orig = beta_cla;
          Complex_wp delta_cla_orig = delta_cla;

          // (beta_cla + eta_cla \xi ^2) \p_r u + (\alpha_cla + \delta_cla \xi^2) u = 0
          alpha_cla = beta_cla_orig+alpha_cla_orig*(coeff-courbure);
          beta_cla = beta_cla_orig*coeff+alpha_cla_orig*ik*ik;
          delta_cla = alpha_cla_orig+coeff*delta_cla_orig;
          eta_cla = delta_cla_orig;
        }
      }

    }

  // is called when the keyword "ConditionReference" is read in the input file.
  // fills periodic, order, gamma and take_curve with read values.
  int ImpedanceABC_HelmholtzBase
    ::GetBoundaryConditionId(const IVect& ref, int pos, const VectString& parameters,
        bool& periodic, int& order, Real_wp& gamma, bool& take_curve)
    {
      if(parameters(pos) == "VISCOTHERMAL" || parameters(pos) == "CREMER")
	{
	  periodic = false;
	  if(parameters.GetM() <=pos +1)
	    {
	      cout << "You need to give the temperature " << endl;
	      cout << "Current parameters are " << parameters << endl;
	      abort();
	    }
	  gamma = to_num<Real_wp>(parameters(pos+1));
	  if( parameters(pos) == "VISCOTHERMAL")
	    {
	      order = VISCOTHERMAL; // arbitrary , could be used to implement condition with and without Laplace Beltrami term
	    }	 
	  if(parameters(pos) == "CREMER")
	    {
	      order = CREMER;
	    }
	  take_curve = true;
	  return BoundaryConditionEnum::LINE_ABSORBING;

	}

      if (parameters(pos) == "ATMOSPHERE")
      {
        periodic = false;
        if (parameters.GetM() <= pos + 2)
        {
          cout << "You need to give the order of the atmosphere boundary condition and alpha parameter" << endl;
          cout << "Current parameters are " << parameters << endl;
          abort();
        }

        gamma = 1.0/to_num<Real_wp>(parameters(pos+2));
        if (parameters(pos+1) == "RBC_1")
          order = ATMO_RBC_1;
        else if (parameters(pos+1) == "SAI_0")
          order = ATMO_SAI_0;
        else if (parameters(pos+1) == "SAI_1")
          order = ATMO_SAI_1;
        else if (parameters(pos+1) == "HF_1")
          order = ATMO_HF_1;
        else if (parameters(pos+1) == "ML_SPHERE")
          order = ATMO_ML_SPHERE;
        else if (parameters(pos+1) == "ML_ANY")
          order = ATMO_ML_ANY;
        else if (parameters(pos+1) == "ABC")
        {
          if (parameters.GetM() <= pos + 3)
          {
            cout << "You need to give the order of the atmosphere boundary condition and alpha parameter" << endl;
            cout << "Current parameters are " << parameters << endl;
            abort();
          }

          int r = to_num<int>(parameters(pos+2));
          if ((parameters(pos+2).size() != 1) || (!isdigit(parameters(pos+2)[0]))
              ||  (r < 1) || (r > 6) )
          {
            cout << "You should given an order between 1 and 6 " << endl;
            cout << "Given order = " << parameters(pos+2) << endl;
            abort();
          }

          order = 30 + r;
          gamma = 1.0/to_num<Real_wp>(parameters(pos+3));
        }
        else if (parameters(pos+1) == "WHIT") // NATHAN
        {
          if (parameters.GetM() <= pos + 3)
          {
            cout << "You need to give the order of the atmosphere boundary condition and alpha parameter" << endl;
            cout << "Current parameters are " << parameters << endl;
            abort();
          }

          int r = to_num<int>(parameters(pos+2));
          if ((parameters(pos+2).size() != 1) || (!isdigit(parameters(pos+2)[0]))
              ||  (r < 1) || (r > 3) )
          {
            cout << "You should given an order between 1 and 3 " << endl;
            cout << "Given order = " << parameters(pos+2) << endl;
            abort();
          }

          order = 40 +r;
          gamma = Real_wp(1.0)/to_num<Real_wp>(parameters(pos+3));

        }
        else // ATMOSPHERE followed with a number (Victor conditions)
        {
          int r = to_num<int>(parameters(pos+1));
          if ((parameters(pos+1).size() != 1) || (!isdigit(parameters(pos+1)[0]))
              ||  (r < 1) || (r > 3) )
          {
            cout << "You should given an order between 1 and 3 " << endl;
            cout << "Given order = " << parameters(pos+1) << endl;
            abort();
          }

          order = 20 + r;
        }

        take_curve = true;
        return BoundaryConditionEnum::LINE_ABSORBING;
      }// end if ATMOSPHERE

      return 0;
    }


  //! Sets values of attributes to default values
  template<class Dimension>
    void ImpedanceABC_Helm<Dimension>::InitDefaultValues()
    {
      num_row_phi = -1;
      num_col_u = -1;
      add_gibc_term = false; 
      epsilon_pade = 0;
      k_infty = 0;
      coef_abc = 0; coef_abc_grad = 0;
      this->read_param_condition = false;
    }


  template<class Dimension>
  Complex_wp ImpedanceABC_Helm<Dimension>
  ::GetCoefficient(int i, int num_elem, int num_loc, int k, int ref_d, int ref,
                   const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat)
  {
    GlobalGenericMatrix<Complex_wp> nat_mat;
    this->EvaluateImpedancePhi(i, num_elem, 0, num_loc, k, nat_mat,
                               ref_d, Pts, Mat);
    
    Vector<Complex_wp> f_phi(1); TinyVector<Real_wp, 1> phi;
    TinyVector<Real_wp, Dimension::dim_N> grad_phi;
    phi(0) = 1.0;
    this->ApplyImpedancePhi_H1(0, k, 0, phi, grad_phi, f_phi);
    return f_phi(0);
  }
  
  
  //! computation of impedance
  template<class Dimension>
  void ImpedanceABC_Helm<Dimension>
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                         const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
                         const SetPoints<Dimension>& Pts,
                         const SetMatrices<Dimension>& Mat)
  {
    if (gibc_condition && add_gibc_term && var_gibc.UsePadeApproximationSquareRoot())
      { 
        this->coef_phi = this->coef_abc;
        return;
      }
    
    // variable impedance => allocating arrays
    if (k == 0)
      {
        int nb_points = Mat.GetNbPointsQuadratureBoundary();
        this->stored_coef_phi(0).Reallocate(nb_points);
        this->stored_coef_grad(0).Reallocate(nb_points);
        this->stored_coef_phi(0).Fill(0);
        this->stored_coef_grad(0).Fill(0);
        
        // storing the normale
        this->vec_normale.Reallocate(nb_points);
        for (int kp = 0; kp < nb_points; kp++)
          this->vec_normale(kp) = Mat.GetNormaleQuadratureBoundary(kp);
      }
    
    int num_point = var_problem.GetNbPointsQuadratureInside(num_elem);
    for (int pos = 0; pos < num_loc; pos++)
      num_point += var_problem.GetMeshNumbering(0).
        GetNbPointsQuadratureBoundary(var_problem.mesh.Element(num_elem).numBoundary(pos));
    
    // dealing with anisotropy c = sqrt( n \cdot mu n / rho)
    typename Dimension::R_N normale = Mat.GetNormaleQuadratureBoundary(k);
    Complex_wp mu_n = var_helm.GetMuNormale(ref, num_elem, num_point+k, normale);
    Complex_wp rho = var_helm.GetRhoTilde(ref, num_elem, num_point+k);
    Complex_wp c = sqrt(mu_n/rho);
    
    if (!var_problem.FirstOrderFormulation())
      {
        Complex_wp feval, coef_laplace;
        SetComplexZero(coef_laplace);
        SetComplexZero(feval);

        // -i omega mu/c = mu ik
        Complex_wp ik = -Iwp*var_problem.GetOmega()/c;
        Complex_wp myk = -var_problem.GetOmega()/c;
        
        Complex_wp alpha_cla(1, 0);
        Complex_wp beta_cla = ik;
        Complex_wp delta_cla(0, 0), eta_cla(0, 0);
        
        Real_wp courbure = 0.0;
        Real_wp h =	0.0;
        Real_wp g = 0.0;
        if (var_boundary.take_into_account_curvature_for_abc)
          {
            Real_wp k1 = Mat.GetK1QuadratureBoundary(k);
            Real_wp k2 = Mat.GetK2QuadratureBoundary(k);
            h= 0.5*(k1+k2); // h = 1/R in 3-D, 1/(2R) in 2-D, mean curvature
            g= k1*k2; // Gauss curvature
            courbure = 2.0*h;
            // h is multiplied by two to coincide with the definition of the 2-D curvature
            // absorbing boundary condition is however correct in 3-D as well
          }
        
        ComputeCoefficientsABC(myk, courbure, var_boundary.GetOrderAbsorbingCondition(),
                               var_boundary.grazing_abc, alpha_cla, beta_cla, delta_cla, eta_cla,
                               var_boundary.gamma_cla_coef, var_boundary.theta_cla_coef,
                               var_boundary.zeta_cla_coef, h, g);
        
        if (gibc_condition)
          {
            Complex_wp gamma(0, 0);
            bool unsym = var_gibc.UseUnsymmetricImplementation();
            if (!unsym)
              {
                // gamma = beta - delta alpha / eta
                if (eta_cla == Complex_wp(0,0))
                  {
                    cout << "Symmetrization does not work if eta is equal to 0" << endl;
                    abort();
                  }
                
                gamma = beta_cla - delta_cla*alpha_cla/eta_cla;
                gamma *= mu_n;
              }
            
            if (add_gibc_term)
              {
                if (num_row_phi >= 0)
                  {
                    if (num_col_u >= 0)
                      {
                        if (unsym)
                          {
                            // unsymmetric version :
                            // term psi \tilde{psi} => alpha_cla
                            feval = alpha_cla;
                            coef_laplace = eta_cla;
                          }
                        else
                          {
                            // symmetric version
                            // term -\alpha \gamma psi psi
                            feval = -alpha_cla*gamma;
                            coef_laplace = -eta_cla*gamma;
                          }
                      }
                    else
                      {
                        if (unsym)
                          {
                            // unsymmetric version
                            // term u \tilde{psi}  => beta_cla
                            feval = beta_cla;
                            coef_laplace = delta_cla;
                          }
                        else
                          {
                            // symmetric version
                            // term gamma u \tilde{psi}
                            feval = gamma;
                            SetComplexZero(coef_laplace);
                          }
                      }
                  }
                else
                  {
                    if (num_col_u >= 0)
                      {
                        if (unsym)
                          {
                            // term psi phi => -mu
                            feval = -mu_n;
                          }
                        else
                          {
                            // symmetric version
                            // term gamma phi psi
                            feval = gamma;
                          }
                      }
                  }
              }
            else
              {
                if (unsym)
                  {
                    // unsymmetric version
                    // term u phi => 0
                    SetComplexZero(feval);
                  }
                else
                  {
                    // symmetric version
                    // term mu delta/eta u phi
                    feval = mu_n*delta_cla/eta_cla;
                  }
              }
          }
        else
          {
            feval = mu_n*beta_cla/alpha_cla;
            coef_laplace = mu_n*delta_cla/alpha_cla;
          }
        
        this->stored_coef_phi(0)(k) = feval*nat_mat.GetCoefDamping();
        this->stored_coef_grad(0)(k) = coef_laplace*nat_mat.GetCoefDamping();
      }
    else
      {
        // first order formulation
        this->stored_coef_phi(0)(k) = rho*c*nat_mat.GetCoefDamping();
        
        if ((var_boundary.GetOrderAbsorbingCondition() >= 2) || (var_boundary.grazing_abc))
          {
            cout << "not implemented" << endl;
            abort();
          }
      }
  }
  
  
  //! computation of impedance for gradient
  template<class Dimension>
  void ImpedanceABC_Helm<Dimension>
  ::EvaluateImpedanceGrad(int i, int num_elem, int num_edge, int num_loc, int k,
                          const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
                          const SetPoints<Dimension>& Pts,
                          const SetMatrices<Dimension>& Mat)
  {
    if (gibc_condition && add_gibc_term && var_gibc.UsePadeApproximationSquareRoot())
      { 
        Real_wp eps = epsilon_pade;
        if (eps < 0)
          {
            if (k == 0)
              {
                int nb_points = Mat.GetNbPointsQuadratureBoundary();
                this->stored_coef_grad(0).Reallocate(nb_points);
                this->stored_coef_grad(0).Fill(0);
              }
            
            Real_wp k1 = Mat.GetK1QuadratureBoundary(k);
            Real_wp k2 = Mat.GetK2QuadratureBoundary(k);
            Real_wp Kappa = k1+k2;
            
            eps = 0.4*pow(k_infty, Real_wp(1)/3)*pow(abs(Kappa), Real_wp(2)/3);
          }
        
        Complex_wp keps = k_infty + Iwp*eps;
        Complex_wp mu_eps = 1.0/square(keps);
        
        if (eps < 0)
          this->stored_coef_grad(0)(k) = this->coef_abc_grad*mu_eps;
        else
          this->coef_grad = this->coef_abc_grad*mu_eps;
        
        return;
      }
  }  
  
  
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
  template<class Dimension>
  Complex_wp ImpedanceHighConductivity<Complex_wp, HelmholtzEquation<Dimension> >
  ::alpha(sqrt(Real_wp(0.5)), -sqrt(Real_wp(0.5)));
  
  
  template<class Dimension>
  ImpedanceHighConductivity<Complex_wp, HelmholtzEquation<Dimension> >
  ::ImpedanceHighConductivity(const EllipticProblem<HelmholtzEquation<Dimension> >& var)
    : ImpedanceFunction_Base<Complex_wp, Dimension>(var),
      var_helm(var), var_boundary(var), var_problem(var)
  {
  }
  
  
  //! impedance for high-conductivity boundary condition
  template<class Dimension>
  void ImpedanceHighConductivity<Complex_wp, HelmholtzEquation<Dimension> >::
  EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_d,
                       const SetPoints<Dimension>& Pts,
                       const SetMatrices<Dimension>& MatricesElem)
  {
    Complex_wp feval(0, 0);
    int ref = var_problem.mesh.BoundaryRef(i).GetReference();
    Complex_wp eps = var_boundary.GetParamCondition(ref, 0);
    
    int num_point = var_problem.GetNbPointsQuadratureInside(num_elem);
    for (int pos = 0; pos < num_loc; pos++)
      num_point += var_problem.GetMeshNumbering(0).
        GetNbPointsQuadratureBoundary(var_problem.mesh.Element(num_elem).numBoundary(pos));
    
    typename Dimension::R_N normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    Complex_wp mu_n = var_helm.GetMuNormale(ref_d, num_elem, num_point+k, normale);
    //Complex_wp rho = var_helm.GetRhoTilde(ref_d, num_elem, num_point+k);
    
    Real_wp k1 = MatricesElem.GetK1QuadratureBoundary(k);
    Real_wp k2 = MatricesElem.GetK2QuadratureBoundary(k);
    Real_wp h = 0.5*(k1+k2);
    Real_wp g = k1*k2;
    Real_wp omega2 = var_problem.GetSquareOmega();
    // DISP(h); DISP(MatricesElem.K1_curve(k)); DISP(MatricesElem.K2_curve(k));
    switch (var_boundary.GetHighConductivityOrder())
      {
      case 1 :
        feval = alpha/eps * mu_n; break;
      case 2 :
        feval = (alpha/eps + h) * mu_n; break;
      case 3 :
        feval = (alpha/eps + h - Real_wp(0.5)*eps/alpha*(h*h - g + omega2) ) * mu_n; break;
      }
    
    if (k == 0)
      {
        int nb_points = MatricesElem.GetNbPointsQuadratureBoundary();
        this->read_param_condition = false;
        this->stored_coef_phi(0).Reallocate(nb_points);
        this->stored_coef_grad(0).Reallocate(nb_points);
        this->stored_coef_phi(0).Fill(0);
        this->stored_coef_grad(0).Fill(0);
        
        // storing the normale
        this->vec_normale.Reallocate(MatricesElem.GetNbPointsQuadratureBoundary());
        for (int kp = 0; kp < MatricesElem.GetNbPointsQuadratureBoundary(); kp++)
          this->vec_normale(kp) = MatricesElem.GetNormaleQuadratureBoundary(kp);
      }
    
    this->stored_coef_phi(0)(k) = feval*nat_mat.GetCoefStiffness();
  }
  
  
  //! impedance related to dphi/ds for high-conductivity boundary condition
  template<class Dimension>
  void ImpedanceHighConductivity<Complex_wp, HelmholtzEquation<Dimension> >::
  EvaluateImpedanceGrad(int i, int num_elem, int num_edge, int num_loc, int k,
                        const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_d,
                        const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat)
  {
    
    int ref = var_problem.mesh.BoundaryRef(i).GetReference();
    Complex_wp eps = var_boundary.GetParamCondition(ref, 0);
    
    int num_point = var_problem.GetNbPointsQuadratureInside(num_elem);
    for (int pos = 0; pos < num_loc; pos++)
      num_point += var_problem.GetMeshNumbering(0).
        GetNbPointsQuadratureBoundary(var_problem.mesh.Element(num_elem).numBoundary(pos));
    
    typename Dimension::R_N normale = Mat.GetNormaleQuadratureBoundary(k);
    Complex_wp mu_n = var_helm.GetMuNormale(ref_d, num_elem, num_point+k, normale);
    //Complex_wp rho = var_helm.GetRhoTilde(ref_d, num_elem, num_point+k);

    Complex_wp feval(0, 0);
    if (var_boundary.GetHighConductivityOrder() == 3)
      feval = Real_wp(0.5)*eps/alpha*mu_n;
    
    this->stored_coef_grad(0)(k) = feval*nat_mat.GetCoefStiffness();
  }
#endif  
  
  
  //! initialization of GIBC
  template<class Dimension>
  void VarGeneralizedImpedance_Helm<Dimension>::TreatGibc(const IVect& Epart)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    if (var_helm.UseAdditionalUnknownForABC() || (this->pade_approx_sqrt))
      {
        ImpedanceABC_Helm<Dimension>::gibc_condition = true;
        
        int n = 0;
        for (int i = 1; i <= var_problem.mesh.GetNbReferences(); i++)
          if (var_problem.mesh.GetBoundaryCondition(i) == BoundaryConditionEnum::LINE_ABSORBING)
            n++;
        
        this->gibc_references.Reallocate(n);
        n = 0;
        for (int i = 1; i <= var_problem.mesh.GetNbReferences(); i++)
          if (var_problem.mesh.GetBoundaryCondition(i) == BoundaryConditionEnum::LINE_ABSORBING)
            this->gibc_references(n++) = i;
        
        if (nb_proc > 1)
          {
            cout << "not yet implemented in parallel" << endl;
            abort();
          }
      }
    
    VarGeneralizedImpedance_Base::FindDofsGibc(var_problem, Epart);
  }
  
  
  //! initialization of GIBC
  template<class Dimension>
  void VarGeneralizedImpedance_Helm<Dimension>::InitGIBC()
  {
    if (this->gibc_references.GetM() > 0)
      {
        int nb_dof_surf = this->GetNbVolumeDofOnSurface();
        
        if (this->pade_approx_sqrt)
          var_boundary.ResizeNbDof(var_problem.GetNbDof() + nb_dof_surf*this->CoefAl.GetM());
        else
          var_boundary.ResizeNbDof(var_problem.GetNbDof() + nb_dof_surf);
        
        var_boundary.NewRowNumbers_Impedance.Reallocate(var_problem.GetNbDof());
        var_boundary.NewColumnNumbers_Impedance.Reallocate(var_problem.GetNbDof());
        
        var_boundary.NewRowNumbers_Impedance.Fill(-1);
        var_boundary.NewColumnNumbers_Impedance.Fill(-1);
        
        var_boundary.ProcRowNumbers_Impedance.Reallocate(var_problem.GetNbDof());
        var_boundary.ProcColumnNumbers_Impedance.Reallocate(var_problem.GetNbDof());
        
        var_boundary.ProcRowNumbers_Impedance.Fill(0);
        var_boundary.ProcColumnNumbers_Impedance.Fill(0);
      }
  }
  
  
  //! adding terms of GIBC boundary condition
  template<class Dimension>
  void VarGeneralizedImpedance_Helm<Dimension>
  ::AddGibcTerms(const Complex_wp& alpha, const GlobalGenericMatrix<Complex_wp>& nat_mat,
                 VirtualMatrix<Complex_wp>& mat_sp, int offset_row, int offset_col)
  {
    ImpedanceFunction_Base<Complex_wp, Dimension>* fimped_ptr
      = var_boundary.GetNewImpedanceABC(Complex_wp(0));
    
    ImpedanceABC_Helm<Dimension>& fimped
      = dynamic_cast<ImpedanceABC_Helm<Dimension>& >(*fimped_ptr);
    
    fimped.add_gibc_term = true;
    
    if (this->pade_approx_sqrt)
      {        
        Real_wp kinf = var_problem.GetOmega()*sqrt(abs(var_helm.rho0/var_helm.mu0));
        fimped.epsilon_pade = this->epsilon_pade;
        fimped.k_infty = kinf;
        
        //DISP(this->CoefC0); DISP(this->CoefAl); DISP(this->CoefBl); DISP(keps);
        // term -ik (C0 + \sum A_l / B_l) mu_0 \int u phi
        Complex_wp C0_tilde = this->CoefC0;
        for (int l = 0; l < this->CoefAl.GetM(); l++)
          C0_tilde += this->CoefAl(l) / this->CoefBl(l);
        
        fimped.coef_abc = -Iwp*kinf*C0_tilde*var_helm.mu0; fimped.coef_abc_grad = 0;
        var_boundary
          .AddMatrixImpedanceBoundary(alpha, ref_cond, 1, nat_mat, mat_sp,
                                      offset_row, offset_col, fimped, false, false);
        
        for (int l = 0; l < this->CoefAl.GetM(); l++)
          {
            // changing row and column numbers to contain dof numbers of psi_l
            this->SetModifiedRowNumbers(l);
            this->SetModifiedColNumbers(l);
            
            // term ik mu0 A_l/B_l \int psi_l \phi
            fimped.coef_abc = Iwp*kinf*this->CoefAl(l)/this->CoefBl(l)*var_helm.mu0; fimped.coef_abc_grad = 0;
            var_boundary.
              AddMatrixImpedanceBoundary(alpha, ref_cond, 1, nat_mat, mat_sp,
                                         offset_row, offset_col, fimped, true, false);
            
            // term ik mu0 A_l/B_l \int u \phi_l
            fimped.coef_abc = Iwp*kinf*this->CoefAl(l)/this->CoefBl(l)*var_helm.mu0; fimped.coef_abc_grad = 0;
            var_boundary.
              AddMatrixImpedanceBoundary(alpha, ref_cond, 1, nat_mat, mat_sp,
                                         offset_row, offset_col, fimped, false, true);
            
            // term ik mu0 A_l/k_eps^2 \int \nabla \psi_l \nabla \phi_l
            // - ik mu0 A_l/B_l \int psi_l phi_l
            fimped.coef_abc = -Iwp*kinf*this->CoefAl(l)/this->CoefBl(l)*var_helm.mu0;
            fimped.coef_abc_grad = Iwp*kinf*var_helm.mu0*this->CoefAl(l);
            var_boundary.
              AddMatrixImpedanceBoundary(alpha, ref_cond, 1, nat_mat, mat_sp,
                                         offset_row, offset_col, fimped, true, true);
          }
        
        return;
      }
    
    // changing row and column numbers to contain dof numbers of psi
    this->SetModifiedRowNumbers(0);
    this->SetModifiedColNumbers(0);
    
    // term - \int mu psi phi ds 
    fimped.num_row_phi = -1; fimped.num_col_u = 0;
    var_boundary.
      AddMatrixImpedanceBoundary(alpha, ref_cond, 1, nat_mat, mat_sp,
                                 offset_row, offset_col, fimped, true, false);
    
    // term \int (alpha - \eta \Delta_\perp) psi psi^T
    fimped.num_row_phi = 0; fimped.num_col_u = 0;
    var_boundary.
      AddMatrixImpedanceBoundary(alpha, ref_cond, 1, nat_mat, mat_sp,
                                 offset_row, offset_col, fimped, true, true);
    
    // term \int (beta - \delta \Delta_\perp) u psi^T
    fimped.num_row_phi = 0; fimped.num_col_u = -1;
    var_boundary.
      AddMatrixImpedanceBoundary(alpha, ref_cond, 1, nat_mat, mat_sp,
                                 offset_row, offset_col, fimped, false, true);

    delete fimped_ptr;
  }
  
}

#define MONTJOIE_FILE_IMPEDANCE_HELMHOLTZ_CXX
#endif

