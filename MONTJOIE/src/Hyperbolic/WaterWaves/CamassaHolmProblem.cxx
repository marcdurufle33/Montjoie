#ifndef MONTJOIE_FILE_CAMASSA_HOLM_PROBLEM_CXX

namespace Montjoie
{
    
  CamassaHolmEquation::CamassaHolmEquation() : GenericWaterWaveEquation()
  {
  }
  
  
  template<class GenericPb>
  void CamassaHolmEquation::ComputePhysicalProperty(const GenericPb& var, double xmin, double xmax)
  {      
    // we approximate bottom with higher degree in order to have
    // a correct approximation of Cxxx
    InitOrder(var.mesh.GetOrder() + 3, var);
    int Nlob = lob.GetNbPointsQuad();
    VectReal_wp Xhigh(Nlob*var.mesh.GetNbElt());
    int nb = 0;
    
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      {
	Real_wp h = var.mesh.Vertex(i+1) - var.mesh.Vertex(i);
	for (int j = 0; j < Nlob; j++)
	  Xhigh(nb++) = var.mesh.Vertex(i) + h*lob.Points(j);
      }
      
    VectReal_wp Bhigh, Chigh;
    ComputeBottom(xmin, xmax, Xhigh, Bhigh, Chigh);
    
    // computation of Cx, Cxx and Cxxx
    VectReal_wp CxDof;
    GetGradient(var, Chigh, CxDof);
        
    GetUquadrature(var, Chigh, Cquad);
    GetUquadrature(var, CxDof, CxQuad);
        
    int Nquad = Cquad.GetM();
    invCquad.Reallocate(Nquad);
    invC2quad.Reallocate(Nquad);
    for (int i = 0; i < Nquad; i++)
      {
        invCquad(i) = 1.0/Cquad(i);
        invC2quad(i) = invCquad(i)*invCquad(i);
      }
    
    if (type_model == CamassaHolmEquation::STRONG)
      {
        VectReal_wp  Atwelve, Btwelve;
        Atwelve = Chigh;
        Btwelve = Chigh;
        for (int i = 0; i < Chigh.GetM(); i++)
          {
            Atwelve(i) = sqrt( pow(Chigh(i), 5.0)/6.0 );
            Btwelve(i) = sqrt( Chigh(i)/12.0 );
          }
        
        VectReal_wp AxDof, AxxDof, AxxxDof, BxDof, BxxDof, BxxxDof;
        GetGradient(var, Atwelve, AxDof);
        GetGradient(var, AxDof, AxxDof);
        GetGradient(var, AxxDof, AxxxDof);
        
        GetGradient(var, Btwelve, BxDof);
        GetGradient(var, BxDof, BxxDof);
        GetGradient(var, BxxDof, BxxxDof);
    
        VectReal_wp AQuad, AxQuad, AxxQuad, AxxxQuad, BQuad, BxQuad, BxxQuad, BxxxQuad;
        GetUquadrature(var, Atwelve, AQuad);
        GetUquadrature(var, AxDof, AxQuad);
        GetUquadrature(var, AxxDof, AxxQuad);
        GetUquadrature(var, AxxxDof, AxxxQuad);
        
        GetUquadrature(var, Btwelve, BQuad);
        GetUquadrature(var, BxDof, BxQuad);
        GetUquadrature(var, BxxDof, BxxQuad);
        GetUquadrature(var, BxxxDof, BxxxQuad);

        CoefZeta.Reallocate(Nquad);
        CoefDxZeta.Reallocate(Nquad);
        CoefDxxZeta.Reallocate(Nquad);
        CoefDxxxZeta.Reallocate(Nquad);
        for (int i = 0; i < Nquad; i++)
          {
            CoefDxxxZeta(i) = mu*(AQuad(i)*AQuad(i) - BQuad(i)*BQuad(i));
            CoefDxxZeta(i) = 3.0*mu*(AQuad(i)*AxQuad(i) - BQuad(i)*BxQuad(i));
            CoefDxZeta(i) = 3.0*mu*(AQuad(i)*AxxQuad(i) - BQuad(i)*BxxQuad(i));
            CoefZeta(i) = mu*(AQuad(i)*AxxxQuad(i) - BQuad(i)*BxxxQuad(i));
          }
      }
    else if (type_model == CamassaHolmEquation::ORIGINAL)
      {
        int Nh = Chigh.GetM();
        VectReal_wp Atilde(Nh), Etilde(Nh), Ftilde(Nh);
        for (int i = 0; i < Nh; i++)
          {
            double c = Chigh(i);
            Atilde(i) = mu/12*(2.0*pow(c, 5.0) - c);
            Etilde(i) = epsilon*mu/24*(-4.0*pow(c, 3.0) - 3.0/c);
            Ftilde(i) = epsilon*mu/24*(-5.0*pow(c, 3.0) - 9.0/c);
          }
        
        VectReal_wp FxTilde, FxxTilde;
        GetGradient(var, Ftilde, FxTilde);
        GetGradient(var, FxTilde, FxxTilde);
        
        GetUquadrature(var, Atilde, CoefAtilde);
        GetUquadrature(var, Etilde, CoefEtilde);
        GetUquadrature(var, Ftilde, CoefFtilde);
        GetUquadrature(var, FxTilde, CoefFxTilde);
        GetUquadrature(var, FxxTilde, CoefFxxTilde);
        //DISP(CoefAtilde); DISP(CoefEtilde);
        //DISP(CoefFtilde); DISP(CoefFxTilde); DISP(CoefFxxTilde);
      }
  }

  
  CamassaHolmMassOperatorDG::CamassaHolmMassOperatorDG(CamassaHolmEquationDG& var_kdv)
    : var(var_kdv)
  {
  }
    
  
  template<class Vector1, class Vector2>
  void CamassaHolmMassOperatorDG
  ::EvaluateNumericalFlux(int i, int k, const Vector1& Um,
                          const Vector1& Up, Vector2& flux1, Vector2& flux2) const
  {
    flux1.Zero(); flux2.Zero();
    if (k == 0)
      {
        flux1(0, 1) = var.mu/24;
        flux2(0, 1) = var.mu/24;
      }
    else
      {
        flux1(0, 1) = -var.mu/24;
        flux2(0, 1) = -var.mu/24;
      }
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void CamassaHolmMassOperatorDG
  ::EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
                                    const Vector1& Up, Vector2& flux1, Vector2& flux2,
                                    Matrix1& df1_dum, Matrix1& df1_dup,
                                    Matrix1& df2_dum, Matrix1& df2_dup) const
  {
    flux1.Zero(); flux2.Zero();
    
    if (k == 0)
      {
        flux1(0, 1) = var.mu/24;
        flux2(0, 1) = var.mu/24;
      }
    else
      {
        flux1(0, 1) = -var.mu/24;
        flux2(0, 1) = -var.mu/24;
      }
  }
  
  
  template<class Vector1, class Vector2>
  void CamassaHolmMassOperatorDG
  ::EvaluateFunction(int i, const Vector1& u, const Vector1& du, 
                     Vector2& f, Vector2& g, Vector2& fx, Vector2& gx) const
  {
    f(0, 1) = var.mu/12;
    g(0, 0) = 1.0;
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void CamassaHolmMassOperatorDG
  ::EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                       Vector2& f, Vector2& g, Vector2& fx, Vector2& gx,
                       Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                       Matrix1& dfx_du, Matrix1& dgx_du, Matrix1& dfx_ddu, Matrix1& dgx_ddu) const
  {
    f(0, 1) = var.mu/12;
    g(0, 0) = 1.0;
  }

  
  CamassaHolmStiffnessOperatorDG
  ::CamassaHolmStiffnessOperatorDG(CamassaHolmEquationDG& var_kdv) : var(var_kdv)
  {
  }
    
  
  template<class Vector1, class Vector2>
  void CamassaHolmStiffnessOperatorDG
  ::EvaluateFunctionFlux(int i, const Vector1& u, Vector2& f) const
  {
    f(1) = -u(0);
    f(2) = -u(1);
    double coef = var.epsilon*var.mu*7.0/24;
    switch (var.type_model)
      {
      case CamassaHolmEquation::GENTLE :
        {
          f(0) = u(0)*(var.Cquad(i) + 0.75*var.epsilon*u(0)
                       *(1.0 - 0.5*var.epsilon*u(0)*(1.0 - 0.375*var.epsilon*u(0))))
            + var.mu/12*u(2) + coef*(0.5*u(1)*u(1) + u(0)*u(2));
        }
        break;
      case CamassaHolmEquation::STRONG :
        {
          f(0) = u(0)*(var.Cquad(i) + 0.75*var.epsilon*var.invCquad(i)*u(0)
                       *(1.0 - 0.5*var.epsilon*var.invC2quad(i)*u(0)
                         *(1.0 - 0.375*var.invC2quad(i)*var.epsilon*u(0))))
            + var.CoefDxxxZeta(i)*u(2) + coef*(0.5*u(1)*u(1) + u(0)*u(2));
        }
        break;
      case CamassaHolmEquation::ORIGINAL :
        {
          f(0) = u(0)*(var.Cquad(i) + 0.75*var.epsilon*var.invCquad(i)*u(0)
                       *(1.0 - 0.5*var.epsilon*var.invC2quad(i)*u(0)
                         *(1.0 - 0.375*var.invC2quad(i)*var.epsilon*u(0))))
            + var.CoefAtilde(i)*u(2) - var.CoefEtilde(i)*(0.5*u(1)*u(1) + u(0)*u(2));
        }
        break;
      }
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void CamassaHolmStiffnessOperatorDG
  ::EvaluateDerivativeFlux(int i, const Vector1& u, Vector2& f, Matrix1& df_du) const
  {
    f(1) = -u(0);
    f(2) = -u(1);
    df_du(1, 0) = -1.0;
    df_du(2, 1) = -1.0;
    double coef = var.epsilon*var.mu*7.0/24;
    switch (var.type_model)
      {
      case CamassaHolmEquation::GENTLE :
        {
          f(0) = u(0)*(var.Cquad(i) + 0.75*var.epsilon*u(0)
                       *(1.0 - 0.5*var.epsilon*u(0)*(1.0 - 0.375*var.epsilon*u(0))))
            + var.mu/12*u(2) + coef*(0.5*u(1)*u(1) + u(0)*u(2));
          
          df_du(0, 0) = var.Cquad(i) + 1.5*var.epsilon*u(0)
            *(1.0 - 0.25*var.epsilon*u(0)*(1.0 - 0.5*var.epsilon*u(0))) + coef*u(2);
          df_du(0, 1) = coef*u(1);
          df_du(0, 2) = var.mu/12 + coef*u(0);
        }
        break;
      case CamassaHolmEquation::STRONG :
        {
          f(0) = u(0)*(var.Cquad(i) + 0.75*var.epsilon*var.invCquad(i)*u(0)
                       *(1.0 - 0.5*var.epsilon*var.invC2quad(i)*u(0)
                         *(1.0 - 0.375*var.invC2quad(i)*var.epsilon*u(0))))
            + var.CoefDxxxZeta(i)*u(2) + coef*(0.5*u(1)*u(1) + u(0)*u(2));

          df_du(0, 0) = var.Cquad(i) + 1.5*var.epsilon*var.invCquad(i)
            *u(0)*(1.0 - 0.25*var.invC2quad(i)*var.epsilon*u(0)
                   *(1.0 - 0.5*var.invC2quad(i)*var.epsilon*u(0))) + coef*u(2);
          df_du(0, 1) = coef*u(1);
          df_du(0, 2) = var.CoefDxxxZeta(i) + coef*u(0);
        }
        break;
      case CamassaHolmEquation::ORIGINAL :
        {
          f(0) = u(0)*(var.Cquad(i) + 0.75*var.epsilon*var.invCquad(i)*u(0)
                       *(1.0 - 0.5*var.epsilon*var.invC2quad(i)*u(0)
                         *(1.0 - 0.375*var.invC2quad(i)*var.epsilon*u(0))))
            + var.CoefAtilde(i)*u(2) - var.CoefEtilde(i)*(0.5*u(1)*u(1) + u(0)*u(2));
          
          df_du(0, 0) = var.Cquad(i) + 1.5*var.epsilon*var.invCquad(i)
            *u(0)*(1.0 - 0.25*var.invC2quad(i)*var.epsilon*u(0)
                   *(1.0 - 0.5*var.invC2quad(i)*var.epsilon*u(0)))
            - var.CoefEtilde(i)*u(2);
          
          df_du(0, 1) = -var.CoefEtilde(i)*u(1);
          df_du(0, 2) = var.CoefAtilde(i) - var.CoefEtilde(i)*u(0);
        }
        break;
      }
  }
  
  
  template<class Vector1, class Vector2>
  void CamassaHolmStiffnessOperatorDG::
  EvaluateNumericalFlux(int i, int k, const Vector1& Um,
			const Vector1& Up, Vector2& flux) const
  {
    Vector2 fm, fp;
    EvaluateFunctionFlux(i, Um, fm);
    EvaluateFunctionFlux(i, Up, fp);
    
    if (k == 0)
      {
        flux = 0.5*(fp - fm);
        flux(0) -= 0.5*var.alpha_lax*(Um(0) - Up(0));
      }
    else
      {
        flux = 0.5*(fp - fm);
        flux(0) += 0.5*var.alpha_lax*(Um(0) - Up(0));
      }
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void CamassaHolmStiffnessOperatorDG::
  EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um, const Vector1& Up,
                                  Vector2& flux, Matrix1& dflux_dum, Matrix1& dflux_dup) const
  {
    Vector2 fm, fp; 
    Matrix1 dfm_du, dfp_du;
    EvaluateDerivativeFlux(i, Um, fm, dfm_du);
    EvaluateDerivativeFlux(i, Up, fp, dfp_du);
    
    if (k == 0)
      {
        flux = 0.5*(fp - fm);
        dflux_dum = -0.5*dfm_du;
        dflux_dup = 0.5*dfp_du;
	flux(0) -= 0.5*var.alpha_lax*(Um(0) - Up(0));
	dflux_dum(0, 0) -= 0.5*var.alpha_lax;
	dflux_dup(0, 0) += 0.5*var.alpha_lax;
      }
    else
      {
        flux = 0.5*(fp - fm);
        dflux_dum = -0.5*dfm_du;
        dflux_dup = 0.5*dfp_du;
	flux(0) += 0.5*var.alpha_lax*(Um(0) - Up(0));
	dflux_dum(0, 0) += 0.5*var.alpha_lax;
	dflux_dup(0, 0) -= 0.5*var.alpha_lax;
      }
  }
  
  
  template<class Vector1, class Vector2>
  void CamassaHolmStiffnessOperatorDG::EvaluateFunction(int i, const Vector1& u, const Vector1& du,
							Vector2& g, Vector2& h) const
  {
    double eps = var.epsilon, mu = var.mu;
    double coef = 7.0/24*var.epsilon*var.mu;
    
    h(1) = u(1) - du(0);
    h(2) = u(2) - du(1);
    
    switch (var.type_model)
      {
      case CamassaHolmEquation::GENTLE :
        {
          h(0) = 0.5*var.CxQuad(i)*u(0) + du(0)*(var.Cquad(i) + 1.5*eps*u(0)
                                                 *(1.0 - 0.25*eps*u(0)*(1.0-0.5*eps*u(0)))); 
          h(0) += 2.0*coef*u(1)*du(1) + (mu/12 + coef*u(0))*du(2);
        }
        break;
      case CamassaHolmEquation::STRONG :
        {
          double coef_dz = 1.5*eps*u(0)*var.invCquad(i)*
            (1.0 - 0.25*u(0)*eps*var.invC2quad(i)*(1.0 - 0.5*var.invC2quad(i)*eps*u(0)));
          
          h(0) = 0.5*u(0)*var.CxQuad(i)
            *(1.0 - eps*var.invC2quad(i)*u(0)*
              (1.0 - 0.75*var.invC2quad(i)*eps*u(0)*(0.75 - 0.5*var.invC2quad(i)*eps*u(0))))
            + var.CoefZeta(i)*u(0) + var.CoefDxZeta(i)*u(1) + var.CoefDxxZeta(i)*u(2);
          
          h(0) += du(0)*(var.Cquad(i) + coef_dz) + du(1)*2.0*coef*u(1)
            + du(2)*(var.CoefDxxxZeta(i) + coef*u(0));
        }
        break;
      case CamassaHolmEquation::ORIGINAL :
        {
          double coef_dz = 1.5*eps*u(0)*var.invCquad(i)*
            (1.0 - 0.25*u(0)*eps*var.invC2quad(i)*(1.0 - 0.5*var.invC2quad(i)*eps*u(0)));
          
          h(0) = 0.5*u(0)*var.CxQuad(i) - (var.CoefFtilde(i) - 2.0*var.CoefEtilde(i))*u(2)*u(1)
            - 0.5*var.CoefFxTilde(i)*u(0)*u(2) - var.CoefFxTilde(i)*u(1)*u(1)
            - 0.5*var.CoefFxxTilde(i)*u(0)*u(1);
          
          h(0) += du(0)*(var.Cquad(i) + coef_dz) - 2.0*var.CoefEtilde(i)*u(1)*du(1)
            + du(2)*(var.CoefAtilde(i) - var.CoefEtilde(i)*u(0));
        }
        break;
      }
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void CamassaHolmStiffnessOperatorDG
  ::EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                       Vector2& g, Vector2& h,
                       Matrix1& dg_du, Matrix1& dh_du, Matrix1& dg_ddu, Matrix1& dh_ddu) const
  {
    double eps = var.epsilon, mu = var.mu;
    double coef = 7.0/24*var.epsilon*var.mu;
    
    h(1) = u(1) - du(0);
    h(2) = u(2) - du(1);
    
    dh_du(1, 1) = 1.0;    
    dh_du(2, 2) = 1.0;
    
    dh_ddu(1, 0) = -1.0;
    dh_ddu(2, 1) = -1.0;
    
    switch (var.type_model)
      {
      case CamassaHolmEquation::GENTLE :
        {
          Real_wp val = var.Cquad(i) + 1.5*eps*u(0)*(1.0 - 0.25*eps*u(0)*(1.0-0.5*eps*u(0)));
          h(0) = 0.5*var.CxQuad(i)*u(0) + du(0)*val; 
          h(0) += 2.0*coef*u(1)*du(1) + (mu/12 + coef*u(0))*du(2);
          
          dh_du(0, 0) = 0.5*var.CxQuad(i) + du(0)*1.5*eps
            *(1.0 - 0.5*eps*u(0)*(1.0-0.75*eps*u(0))) + coef*du(2);
          dh_du(0, 1) = 2.0*coef*du(1);
          
          dh_ddu(0, 0) = val;
          dh_ddu(0, 1) = 2.0*coef*u(1);
          dh_ddu(0, 2) = mu/12 + coef*u(0);
        }
        break;
      case CamassaHolmEquation::STRONG :
        {
          double coef_dz = 1.5*eps*u(0)*var.invCquad(i)*
            (1.0 - 0.25*u(0)*eps*var.invC2quad(i)*(1.0 - 0.5*var.invC2quad(i)*eps*u(0)));
                    
          h(0) = 0.5*u(0)*var.CxQuad(i)
            *(1.0 - eps*var.invC2quad(i)*u(0)*
              (1.0 - 0.75*var.invC2quad(i)*eps*u(0)*(0.75 - 0.5*var.invC2quad(i)*eps*u(0))))
            + var.CoefZeta(i)*u(0) + var.CoefDxZeta(i)*u(1) + var.CoefDxxZeta(i)*u(2);
          
          h(0) += du(0)*(var.Cquad(i) + coef_dz) 
            + du(1)*2.0*coef*u(1) + du(2)*(var.CoefDxxxZeta(i) + coef*u(0));
          
          dh_du(0, 0) = 0.5*var.CxQuad(i)
            - 0.5*eps*var.invC2quad(i)*u(0)*var.CxQuad(i)
            *(2.0 - 0.75*var.invC2quad(i)*eps*u(0)*(2.25 - 2.0*var.invC2quad(i)*eps*u(0)))
            + var.CoefZeta(i);

          dh_du(0, 0) += du(0)*(1.5*eps*var.invCquad(i)
                                *(1.0 - 0.5*eps*var.invC2quad(i)*u(0)
                                  *(1.0 - 0.75*var.invC2quad(i)*eps*u(0))));
          
          dh_du(0, 1) = var.CoefDxZeta(i);
          dh_du(0, 2) = var.CoefDxxZeta(i);
          
          dh_ddu(0, 0) = var.Cquad(i) + coef_dz;
          dh_ddu(0, 1) = 2.0*coef*u(1);
          dh_ddu(0, 2) = var.CoefDxxxZeta(i) + coef*u(0);
        }
        break;
      case CamassaHolmEquation::ORIGINAL :
        {
          double coef_dz = 1.5*eps*u(0)*var.invCquad(i)*
            (1.0 - 0.25*u(0)*eps*var.invC2quad(i)*(1.0 - 0.5*var.invC2quad(i)*eps*u(0)));
          
          h(0) = 0.5*u(0)*var.CxQuad(i) - (var.CoefFtilde(i) - 2.0*var.CoefEtilde(i))*u(2)*u(1)
            - 0.5*var.CoefFxTilde(i)*u(0)*u(2) - var.CoefFxTilde(i)*u(1)*u(1)
            - 0.5*var.CoefFxxTilde(i)*u(0)*u(1);
          
          h(0) += du(0)*(var.Cquad(i) + coef_dz) - 2.0*var.CoefEtilde(i)*u(1)*du(1)
            + du(2)*(var.CoefAtilde(i) - var.CoefEtilde(i)*u(0));
          
          dh_du(0, 0) = 0.5*var.CxQuad(i) - 0.5*var.CoefFxTilde(i)*u(2)
            - 0.5*var.CoefFxxTilde(i)*u(1);
          
          dh_du(0, 0) = 1.5*eps*var.invCquad(i)*du(0)
            *(1.0 - 0.5*eps*var.invC2quad(i)*u(0)*(1.0 - 0.75*var.invC2quad(i)*eps*u(0)))
            - var.CoefEtilde(i)*du(2);                    

          dh_du(0, 1) = -(var.CoefFtilde(i) - 2.0*var.CoefEtilde(i))*u(2)
            - 2.0*var.CoefFxTilde(i)*u(1) - 0.5*var.CoefFxxTilde(i)*u(0)
            - 2.0*var.CoefEtilde(i)*du(1);
          
          dh_du(0, 2) = -(var.CoefFtilde(i) - 2.0*var.CoefEtilde(i))*u(1)
            - 0.5*var.CoefFxTilde(i)*u(0);
          
          dh_ddu(0, 0) = var.Cquad(i) + coef_dz;
          dh_ddu(0, 1) = - 2.0*var.CoefEtilde(i)*u(1);
          dh_ddu(0, 2) = var.CoefAtilde(i) - var.CoefEtilde(i)*u(0);
        }
        break;
      }
  }

  
  CamassaHolmEquationDG::CamassaHolmEquationDG() : var_mass(*this), var_stiff(*this)
  {
    alpha_lax = 5.0;
  }
  
  
  CamassaHolmEquationDG::CamassaHolmEquationDG(const CamassaHolmEquationDG& var)
    : CamassaHolmEquation(var), var_mass(*this), var_stiff(*this)
  {
    alpha_lax = var.alpha_lax;
  }
  
  
  CamassaHolmMassOperatorFem::CamassaHolmMassOperatorFem(const CamassaHolmEquationFem& var_eq)
    : var(var_eq)
  {    
  }
  
  
  template<class Vector1, class Vector2>
  void CamassaHolmMassOperatorFem::
  EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                   Vector2& f, Vector2& g, Vector2& fx, Vector2& gx) const
  {
    g(0, 0) = 1.0;
    fx(0, 0) = var.mu/12;
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void CamassaHolmMassOperatorFem::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du, Vector2& f,
                     Vector2& g, Vector2& fx, Vector2& gx,
                     Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                     Matrix1& dfx_du, Matrix1& dgx_du, Matrix1& dfx_ddu, Matrix1& dgx_ddu) const
  {
    g(0, 0) = 1.0;
    fx(0, 0) = var.mu/12;
  }
  
  
  CamassaHolmStiffnessOperatorFem
  ::CamassaHolmStiffnessOperatorFem(const CamassaHolmEquationFem& var_eq) : var(var_eq)
  {
  }
    
  
  template<class Vector1, class Vector2>
  void CamassaHolmStiffnessOperatorFem::
  EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                   Vector2& g, Vector2& h) const
  {
    double eps = var.epsilon, mu = var.mu;
    
    h(1) = u(1) - du(0);
    h(2) = u(2) - du(1);
    
    switch (var.type_model)
      {
      case CamassaHolmEquation::GENTLE :
        {
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0)
            + 1.5*eps*u(0)*du(0)*(1.0 - 0.25*u(0)*eps*(1.0 - 0.5*eps*u(0)))
            + mu/12*du(2) + 7.0*eps*mu/24*(u(0)*du(2) + 2.0*du(0)*du(1));
        }
        break;
      case CamassaHolmEquation::STRONG :
        {
          double coef = 7.0*eps*mu/24;
          double coef_dz = 1.5*eps*u(0)*var.invCquad(i)*
            (1.0 - 0.25*u(0)*eps*var.invC2quad(i)*(1.0 - 0.5*var.invC2quad(i)*eps*u(0)));
          
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0) + du(0)*coef_dz
            - 0.5*eps*var.invC2quad(i)*u(0)*u(0)*var.CxQuad(i)
            *(1.0 - 0.75*var.invC2quad(i)*eps*u(0)*(0.75 - 0.5*var.invC2quad(i)*eps*u(0)))
            + var.CoefZeta(i)*u(0) + var.CoefDxZeta(i)*du(0) 
            + var.CoefDxxZeta(i)*du(1) + var.CoefDxxxZeta(i)*du(2) 
            + coef*(u(0)*du(2) + 2.0*du(0)*du(1));
        }
        break;
      case CamassaHolmEquation::ORIGINAL :
        {
          double coef_dz = 1.5*eps*u(0)*var.invCquad(i)*
            (1.0 - 0.25*u(0)*eps*var.invC2quad(i)*(1.0 - 0.5*var.invC2quad(i)*eps*u(0)));
          
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0) + du(0)*coef_dz
            + var.CoefAtilde(i)*du(2) - var.CoefEtilde(i)*u(0)*du(2) 
            - var.CoefFtilde(i)*du(0)*du(1) - var.CoefFxTilde(i)*(0.5*u(0)*du(1) + du(0)*du(0))
            - 0.5*var.CoefFxxTilde(i)*u(0)*du(0);
          
        }
        break;
      }    
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void CamassaHolmStiffnessOperatorFem::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                     Vector2& g, Vector2& h, Matrix1& dg_du,
                     Matrix1& dh_du, Matrix1& dg_ddu, Matrix1& dh_ddu) const
  {
    double eps = var.epsilon, mu = var.mu;
    double coef = 7.0*eps*mu/24;
    
    h(1) = u(1) - du(0);
    h(2) = u(2) - du(1);
    dh_du(1, 1) = 1.0;
    dh_du(2, 2) = 1.0;
    dh_ddu(1, 0) = -1.0;
    dh_ddu(2, 1) = -1.0;
    
    switch (var.type_model)
      {
      case CamassaHolmEquation::GENTLE :
        {
          double coef_dz = 1.5*eps*u(0)*(1.0 - 0.25*u(0)*eps*(1.0 - 0.5*eps*u(0)));
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0)
            + du(0)*coef_dz + mu/12*du(2) + coef*(u(0)*du(2) + 2.0*du(0)*du(1));
          
          dh_du(0, 0) = 0.5*var.CxQuad(i) + 1.5*eps*du(0)
            *(1.0 - 0.5*eps*u(0)*(1.0 - 0.75*eps*u(0))) + coef*du(2);
          dh_ddu(0, 0) = var.Cquad(i) + coef_dz + 2.0*coef*du(1);
          dh_ddu(0, 1) = 2.0*coef*du(0);
          dh_ddu(0, 2) = mu/12 + coef*u(0);
        }
        break;
      case CamassaHolmEquation::STRONG :
        {
          double coef_dz = 1.5*eps*u(0)*var.invCquad(i)*
            (1.0 - 0.25*u(0)*eps*var.invC2quad(i)*(1.0 - 0.5*var.invC2quad(i)*eps*u(0)));
          
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0) + du(0)*coef_dz
            - 0.5*eps*var.invC2quad(i)*u(0)*u(0)*var.CxQuad(i)
            *(1.0 - 0.75*var.invC2quad(i)*eps*u(0)*(0.75 - 0.5*var.invC2quad(i)*eps*u(0)))
            + var.CoefZeta(i)*u(0) + var.CoefDxZeta(i)*du(0) 
            + var.CoefDxxZeta(i)*du(1) + var.CoefDxxxZeta(i)*du(2) 
            + coef*(u(0)*du(2) + 2.0*du(0)*du(1));

          dh_du(0, 0) = 0.5*var.CxQuad(i) + 1.5*eps*var.invCquad(i)*du(0)
            *(1.0 - 0.5*eps*var.invC2quad(i)*u(0)*(1.0 - 0.75*var.invC2quad(i)*eps*u(0)))
            - 0.5*eps*var.invC2quad(i)*u(0)*var.CxQuad(i)
            *(2.0 - 0.75*var.invC2quad(i)*eps*u(0)*(2.25 - 2.0*var.invC2quad(i)*eps*u(0)))
            + var.CoefZeta(i) + coef*du(2);
          
          dh_ddu(0, 0) = var.Cquad(i) + coef_dz + var.CoefDxZeta(i) + 2.0*coef*du(1);
          dh_ddu(0, 1) = var.CoefDxxZeta(i) + 2.0*coef*du(0);
          dh_ddu(0, 2) = var.CoefDxxxZeta(i) + coef*u(0);
        }
        break;
      case CamassaHolmEquation::ORIGINAL :
        {
          double coef_dz = 1.5*eps*u(0)*var.invCquad(i)*
            (1.0 - 0.25*u(0)*eps*var.invC2quad(i)*(1.0 - 0.5*var.invC2quad(i)*eps*u(0)));
          
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0) + du(0)*coef_dz
            + var.CoefAtilde(i)*du(2) - var.CoefEtilde(i)*u(0)*du(2) 
            - var.CoefFtilde(i)*du(0)*du(1) - var.CoefFxTilde(i)*(0.5*u(0)*du(1) + du(0)*du(0))
            - 0.5*var.CoefFxxTilde(i)*u(0)*du(0);
          
          dh_du(0, 0) = 0.5*var.CxQuad(i) + 1.5*eps*var.invCquad(i)*du(0)
            *(1.0 - 0.5*eps*var.invC2quad(i)*u(0)*(1.0 - 0.75*var.invC2quad(i)*eps*u(0)))
            - var.CoefEtilde(i)*du(2) - 0.5*var.CoefFxTilde(i)*du(1)
            - 0.5*var.CoefFxxTilde(i)*du(0);
          
          dh_ddu(0, 0) = var.Cquad(i) + coef_dz - var.CoefFtilde(i)*du(1)
            - 2.0*var.CoefFxTilde(i)*du(0) - 0.5*var.CoefFxxTilde(i)*u(0);
          
          dh_ddu(0, 1) = -var.CoefFtilde(i)*du(0) - 0.5*var.CoefFxTilde(i)*u(0);
          dh_ddu(0, 2) = var.CoefAtilde(i) -var.CoefEtilde(i)*u(0);
        }
        break;
      }
   
  }
  
 
  CamassaHolmEquationFem::CamassaHolmEquationFem() : var_mass(*this), var_stiff(*this)
  {
  }
  
  
  CamassaHolmEquationFem::CamassaHolmEquationFem(const CamassaHolmEquationFem& A)
    :  CamassaHolmEquation(A), var_mass(*this), var_stiff(*this)
  {
  }

  
  bool CamassaHolmEquationFem::LinearMassTerm()
  {
    return true;
  }
  
  
  /******************************************
   * Resolution with pseudo-spectral method *
   ******************************************/


  CamassaHolmEquationSpectral::CamassaHolmEquationSpectral() : CamassaHolmEquation()
  {
  }
  
  //! computation of c and cx 
  void CamassaHolmEquationSpectral
  ::ComputePhysicalProperty(FftRealInterface& fft, const VectReal_wp& nu,
			    double xmin, double xmax, int N)
  {      
    VectReal_wp Xdof(N);
    Real_wp dx = (xmax-xmin)/N;
    for (int i = 0; i < N; i++)
      Xdof(i) = xmin + i*dx;
    
    VectReal_wp Bdof;
    ComputeBottom(xmin, xmax, Xdof, Bdof, Cquad);
    
    Cquad.WriteText("bottom_sp.dat");
    invCquad.Reallocate(N);
    invC2quad.Reallocate(N);
    VectReal_wp Aquad(N), Bquad(N);
    for (int i = 0; i < N; i++)
      {
        invCquad(i) = 1.0/Cquad(i);
        invC2quad(i) = invCquad(i)*invCquad(i);
        Aquad(i) = sqrt( pow(Cquad(i), 5.0)/6.0 );
        Bquad(i) = sqrt( Cquad(i)/12.0 );
      }
    
    // evaluating derivatives of c, a_{1/12} and b_{1/12}
    VectComplex_wp ChatDx(N/2+1), AhatDx(N/2+1), AhatDxx(N/2+1), AhatDxxx(N/2+1);
    VectComplex_wp BhatDx(N/2+1), BhatDxx(N/2+1), BhatDxxx(N/2+1);
    
    fft.ApplyForward(Cquad, ChatDx);
    fft.ApplyForward(Aquad, AhatDx);
    fft.ApplyForward(Bquad, BhatDx);
    
    for (int i = 0; i < nu.GetM(); i++)
      {
        ChatDx(i) *= Iwp*nu(i);
        AhatDx(i) *= Iwp*nu(i);
        AhatDxx(i) = Iwp*nu(i)*AhatDx(i);
        AhatDxxx(i) = Iwp*nu(i)*AhatDxx(i);
        BhatDx(i) *= Iwp*nu(i);
        BhatDxx(i) = Iwp*nu(i)*BhatDx(i);
        BhatDxxx(i) = Iwp*nu(i)*BhatDxx(i);
      }
    
    VectReal_wp AquadDx(N), AquadDxx(N), AquadDxxx(N);
    VectReal_wp BquadDx(N), BquadDxx(N), BquadDxxx(N);
    CxQuad.Reallocate(N);
    fft.ApplyInverse(ChatDx, CxQuad);
    fft.ApplyInverse(AhatDx, AquadDx);
    fft.ApplyInverse(AhatDxx, AquadDxx);
    fft.ApplyInverse(AhatDxxx, AquadDxxx);
    fft.ApplyInverse(BhatDx, BquadDx);
    fft.ApplyInverse(BhatDxx, BquadDxx);
    fft.ApplyInverse(BhatDxxx, BquadDxxx);
    
    if (type_model == CamassaHolmEquation::STRONG)
      {
        CoefZeta.Reallocate(N);
        CoefDxZeta.Reallocate(N);
        CoefDxxZeta.Reallocate(N);
        CoefDxxxZeta.Reallocate(N);
        for (int i = 0; i < N; i++)
          {
            CoefDxxxZeta(i) = mu*(Aquad(i)*Aquad(i) - Bquad(i)*Bquad(i));
            CoefDxxZeta(i) = 3.0*mu*(Aquad(i)*AquadDx(i) - Bquad(i)*BquadDx(i));
            CoefDxZeta(i) = 3.0*mu*(Aquad(i)*AquadDxx(i) - Bquad(i)*BquadDxx(i));
            CoefZeta(i) = mu*(Aquad(i)*AquadDxxx(i) - Bquad(i)*BquadDxxx(i));
          }
      }
    else if (type_model == CamassaHolmEquation::ORIGINAL)
      {
        VectReal_wp Atilde(N), Etilde(N), Ftilde(N);
        for (int i = 0; i < N; i++)
          {
            double c = Cquad(i);
            Atilde(i) = mu/12*(2.0*pow(c, 5.0) - c);
            Etilde(i) = epsilon*mu/24*(-4.0*pow(c, 3.0) - 3.0/c);
            Ftilde(i) = epsilon*mu/24*(-5.0*pow(c, 3.0) - 9.0/c);
          }
        
        VectComplex_wp FhatDx(N/2+1), FhatDxx(N/2+1);
        fft.ApplyForward(Ftilde, FhatDx);

        for (int i = 0; i < nu.GetM(); i++)
          {
            FhatDx(i) *= Iwp*nu(i);
            FhatDxx(i) = Iwp*nu(i)*FhatDx(i);
          }

        CoefAtilde = Atilde;
        CoefEtilde = Etilde;
        CoefFtilde = Ftilde;
        CoefFxTilde.Reallocate(N);
        CoefFxxTilde.Reallocate(N);
        
        fft.ApplyInverse(FhatDx, CoefFxTilde);
        fft.ApplyInverse(FhatDxx, CoefFxxTilde);
      }    
    
    CoefMass.Reallocate(nu.GetM());
    for (int i = 0; i < nu.GetM(); i++)
      CoefMass(i) = 1.0/(1.0 + mu/12.0*nu(i)*nu(i));

    ProdReal.Reallocate(N);
    XnReel.Reallocate(N);
    dXreal.Reallocate(N);
    d2Xreal.Reallocate(N);
    d3Xreal.Reallocate(N);
    
    ProdHat.Reallocate(N/2+1);
    XnHat.Reallocate(N/2+1);
    dXn.Reallocate(N/2+1);
    d2Xn.Reallocate(N/2+1);
    d3Xn.Reallocate(N/2+1);

    coef_Gamma1.Reallocate(nu.GetM());
    coef_Gamma2.Reallocate(nu.GetM());
    coef_Gamma3.Reallocate(nu.GetM());
    for (int i = 0; i < nu.GetM(); i++)
      {
        coef_Gamma1(i) = Iwp*nu(i);
        coef_Gamma2(i) = -nu(i)*nu(i);
        coef_Gamma3(i) = -Iwp*pow(nu(i), 3.0);
      }
  }

  
  void CamassaHolmEquationSpectral
  ::InitSinus(const Real_wp& dt, const VectReal_wp& nu)
  {
    coef_Gamma1.Reallocate(nu.GetM());
    coef_Gamma2.Reallocate(nu.GetM());
    coef_Gamma3.Reallocate(nu.GetM());
    for (int i = 0; i < nu.GetM(); i++)
      {
        coef_Gamma1(i) = Iwp/dt*sin(dt*nu(i));
        coef_Gamma2(i) = -1.0/dt*sin(dt*nu(i)*nu(i));
        coef_Gamma3(i) = -Iwp/dt*sin(dt*pow(nu(i), 3.0));
      }
  }


  //! applying mass operator
  void CamassaHolmEquationSpectral
  ::ApplyMass(const VectReal_wp& X, FftRealInterface& fft, VectReal_wp& nu, VectReal_wp& Y)
  {
    fft.ApplyForward(X, XnHat);
    
    for (int i = 0; i < XnHat.GetM(); i++)
      XnHat(i) *= CoefMass(i);
    
    fft.ApplyInverse(XnHat, Y);
  }
  
  void CamassaHolmEquationSpectral::
  AddFunction(const Real_wp& alpha, const VectReal_wp& X,
              FftRealInterface& fft, VectReal_wp& nu,
              VectReal_wp& Y, bool solve_mass)
  {
    VectReal_wp Prod(Y);
    
    fft.ApplyForward(X, XnHat);
    
    GetFunction(XnHat, X, fft, nu, ProdHat, solve_mass);
    
    fft.ApplyInverse(ProdHat, Prod);
    
    Add(alpha, Prod, Y);
  }
  
  
  void CamassaHolmEquationSpectral
  ::GetFunction(const VectComplex_wp& Xn, const VectReal_wp& Xreal,
                FftRealInterface& fft, VectReal_wp& nu, VectComplex_wp& ProdXn,
                bool solve_mass)
  {
    for (int i = 0; i < nu.GetM(); i++)
      {
        dXn(i) = Iwp*nu(i)*Xn(i);
        d2Xn(i) = coef_Gamma2(i)*Xn(i);
        d3Xn(i) = coef_Gamma3(i)*Xn(i);
      }    
    
    fft.ApplyInverse(dXn, dXreal);
    fft.ApplyInverse(d2Xn, d2Xreal);
    fft.ApplyInverse(d3Xn, d3Xreal);
    
    Real_wp coef = 7.0/24.0*epsilon*mu;
    if (type_model == GENTLE)
      {        
        for (int i = 0; i < Xreal.GetM(); i++)
          {
            ProdReal(i) = - Cquad(i)*dXreal(i) - 0.5*CxQuad(i)*Xreal(i) 
              - 1.5*epsilon*Xreal(i)*dXreal(i)
              *(1.0 - 0.25*epsilon*Xreal(i)*(1.0 - 0.5*epsilon*Xreal(i)))
              - mu/12.0*d3Xreal(i) - coef*(Xreal(i)*d3Xreal(i) + 2.0*dXreal(i)*d2Xreal(i));
          }
      }
    else if (type_model == STRONG)
      {
        for (int i = 0; i < Xreal.GetM(); i++)
          {
            ProdReal(i) = - Cquad(i)*dXreal(i) - 0.5*CxQuad(i)*Xreal(i) 
              - 1.5*epsilon*Xreal(i)*dXreal(i)*invCquad(i)
              *(1.0 - 0.25*invC2quad(i)*epsilon*Xreal(i)*(1.0 - 0.5*invC2quad(i)*epsilon*Xreal(i)))
              + 1.5*epsilon*Xreal(i)*Xreal(i)*CxQuad(i)*invC2quad(i)
              *(1.0/3 - 0.25*invC2quad(i)*epsilon*Xreal(i)
                *(3.0/4 - 0.5*invC2quad(i)*epsilon*Xreal(i)))
              - CoefZeta(i)*Xreal(i) - CoefDxZeta(i)*dXreal(i)
              - CoefDxxZeta(i)*d2Xreal(i) - CoefDxxxZeta(i)*d3Xreal(i)
              - coef*(Xreal(i)*d3Xreal(i) + 2.0*dXreal(i)*d2Xreal(i));
          }
      }
    else if (type_model == ORIGINAL)
      {
        for (int i = 0; i < Xreal.GetM(); i++)
          {
            ProdReal(i) = - Cquad(i)*dXreal(i) - 0.5*CxQuad(i)*Xreal(i) 
              - 1.5*epsilon*invCquad(i)*Xreal(i)*dXreal(i)
              *(1.0 - 0.25*epsilon*invC2quad(i)*Xreal(i)*(1.0 - 0.5*invC2quad(i)*epsilon*Xreal(i)))
              - CoefAtilde(i)*d3Xreal(i) + CoefEtilde(i)*Xreal(i)*d3Xreal(i)
              + d2Xreal(i)*(CoefFtilde(i)*dXreal(i) + 0.5*CoefFxTilde(i)*Xreal(i))
              + dXreal(i)*(CoefFxTilde(i)*dXreal(i) + 0.5*CoefFxxTilde(i)*Xreal(i));
          }
      }
    
    fft.ApplyForward(ProdReal, ProdXn);
    
    if (solve_mass)
      for (int i = 0; i < Xn.GetM(); i++)
        ProdXn(i) *= CoefMass(i);
  }
  

  //! evaluation of the operator that evolves in time
  void CamassaHolmEquationSpectral
  ::EvaluateFunction(const Real_wp& t, const VectComplex_wp& X,
		     FftRealInterface& fft, VectReal_wp& nu,
		     VectComplex_wp& Y, bool solve_mass)
  {
    fft.ApplyInverse(X, XnReel);
    
    GetFunction(X, XnReel, fft, nu, Y, solve_mass);    
  }


  //! evaluation of the operator that evolves in time
  void CamassaHolmEquationSpectral
  ::EvaluateFunction(const Real_wp& t, const VectReal_wp& X,
		     FftRealInterface& fft, VectReal_wp& nu,
		     VectReal_wp& Y, bool solve_mass)
  {
    fft.ApplyForward(X, XnHat);
    
    GetFunction(XnHat, X, fft, nu, ProdHat, solve_mass);
    
    fft.ApplyInverse(ProdHat, Y);
  }

}

#define MONTJOIE_FILE_CAMASSA_HOLM_PROBLEM_CXX
#endif
