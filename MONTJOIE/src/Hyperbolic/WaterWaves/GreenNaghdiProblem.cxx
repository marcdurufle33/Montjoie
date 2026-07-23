#ifndef MONTJOIE_FILE_GREEN_NAGHDI_PROBLEM_CXX

namespace Montjoie
{
  
  GreenNaghdiEquation::GreenNaghdiEquation() : CamassaHolmEquation()
  {
  }
    
  
  void GreenNaghdiEquation
  ::ComputeInitialCondition(double xmin, double xmax, const VectReal_wp& Xdof, VectReal_wp& Y0)
  {
    int nodl = Xdof.GetM();
    VectReal_wp Bdof, Cdof;
    ComputeBottom(xmin, xmax, Xdof, Bdof, Cdof);
    
    //if (type_initial == SOLITON_GN)
    if (false)
      {
        double alpha_s = freq_init;
        Y0.Reallocate(2*nodl);
        double k = sqrt( (3.0*alpha_s*epsilon)/(4.0*mu*(1.0 + alpha_s*epsilon) ) );
        double c = sqrt(1+alpha_s*epsilon);
        // DISP(c); DISP(k); DISP(mu); DISP(epsilon);
        for (int i = 0; i < nodl; i++)
          {
            double x = Xdof(i);
            Y0(i) = alpha_s/square(cosh(k*(x - center_init)));
            Y0(nodl+i) = c*Y0(i)/(1.0 + epsilon*Y0(i));
          }
      }
    else
      {
        // computing initial condition
        VectReal_wp Z0, Z0_x(nodl), Z0_xx(nodl);
        CamassaHolmEquation::ComputeInitialCondition(xmin, xmax, Xdof, Z0, Z0_x, Z0_xx);
        Y0.Reallocate(2*nodl); Y0.Fill(0);
        for (int i = 0; i < nodl; i++)
          Y0(i) = Z0(i);
        
        for (int i = 0; i < nodl; i++)
          {
            Real_wp c = Cdof(i), z = Z0(i), z_xx = Z0_xx(i), z_x = Z0_x(i);
            Y0(nodl+i) = (z + c*c/(c*c + epsilon*z)
                          *( -0.25*epsilon/(c*c)*z*z - pow(epsilon, 2.0)/(pow(c, 4.0)*8.0)
                             *pow(z, 3.0)
                             + 3.0*pow(epsilon, 3.0)/(64.0*pow(c, 6.0))*pow(z, 4.0)
                             + mu/6*pow(c, 4.0)*z_xx
                             + epsilon*mu*c*c*(5.0/12*z*z_xx + 13.0/48*z_x*z_x) ) ) / c;
          }
      }
  }

  
  template<class GenericPb>
  void GreenNaghdiEquation::ComputePhysicalProperty(const GenericPb& var, double xmin, double xmax)
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
    
    VectReal_wp Bhigh, Chigh, Cs;
    ComputeBottom(xmin, xmax, Xhigh, Bhigh, Chigh);
    Cs = Chigh;
    for (int i = 0; i < Cs.GetM(); i++)
      Cs(i) = square(Chigh(i));
    
    // computation of c^2 and (c^2)_x on quadrature points
    VectReal_wp CxDof;
    GetGradient(var, Cs, CxDof);
    
    GetUquadrature(var, Cs, Csquare);
    GetUquadrature(var, CxDof, CsquareDx);
  }

    
  GreenNaghdiMassOperatorDG
  ::GreenNaghdiMassOperatorDG(GreenNaghdiEquationDG& var_kdv) : var(var_kdv)
  {
  }
    
  
  template<class Vector1, class Vector2>
  void GreenNaghdiMassOperatorDG
  ::EvaluateNumericalFlux(int i, int k, const Vector1& Um,
                          const Vector1& Up, Vector2& flux1, Vector2&flux2) const
  {
    flux1.Zero(); flux2.Zero();
    
    Real_wp hm = var.Csquare(i) + var.epsilon*Um(0);
    Real_wp hp = var.Csquare(i) + var.epsilon*Up(0);
    
    flux1(1, 2) = var.mu/6*hm*hm;
    flux2(1, 2) = var.mu/6*hp*hp;
    
    if (k == 1)
      {
        flux1(1, 2) = -flux1(1, 2);
        flux2(1, 2) = -flux2(1, 2);
      }
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void GreenNaghdiMassOperatorDG::
  EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
                                  const Vector1& Up, Vector2& flux1, Vector2& flux2,
                                  Matrix1& df1_dum, Matrix1& df1_dup,
                                  Matrix1& df2_dum, Matrix1& df2_dup) const
  {
    flux1.Zero(); flux2.Zero();
    df1_dum.Zero(); df1_dup.Zero(); df2_dum.Zero(); df2_dup.Zero();
    
    Real_wp hm = var.Csquare(i) + var.epsilon*Um(0);
    Real_wp hp = var.Csquare(i) + var.epsilon*Up(0);
    
    flux1(1, 2) = var.mu/6*hm*hm;
    flux2(1, 2) = var.mu/6*hp*hp;
    
    df1_dum(1, 2, 0) = var.mu*var.epsilon/3*hm;
    df2_dup(1, 2, 0) = var.mu*var.epsilon/3*hp;
    
    if (k == 1)
      {
        flux1(1, 2) = -flux1(1, 2);
        flux2(1, 2) = -flux2(1, 2);
        
        df1_dum(1, 2, 0) = -df1_dum(1, 2, 0);
        df2_dup(1, 2, 0) = -df2_dup(1, 2, 0);
      }
    
  }
  
  
  template<class Vector1, class Vector2>
  void GreenNaghdiMassOperatorDG
  ::EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                     Vector2& f, Vector2& g, Vector2& fx, Vector2& gx) const
  {
    Real_wp h = var.Csquare(i) + var.epsilon*u(0);
    Real_wp hx = var.CsquareDx(i) + var.epsilon*du(0);
    
    g(0, 0) = 1.0;
    g(1, 1) = 1.0;
    g(1, 2) = -var.mu/3*h*hx;
    
    f(1, 2) = var.mu/3*h*h;
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void GreenNaghdiMassOperatorDG::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                     Vector2& f, Vector2& g, Vector2& fx, Vector2& gx,
                     Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                     Matrix1& dfx_du, Matrix1& dgx_du, Matrix1& dfx_ddu, Matrix1& dgx_ddu) const
  {
    Real_wp h = var.Csquare(i) + var.epsilon*u(0);
    Real_wp hx = var.CsquareDx(i) + var.epsilon*du(0);

    g(0, 0) = 1.0;
    g(1, 1) = 1.0;
    g(1, 2) = -var.mu/3*h*hx;
    
    f(1, 2) = var.mu/3*h*h;
    
    df_du(1, 2, 0) = 2.0*var.mu*var.epsilon/3*h;
    dg_du(1, 2, 0) = -var.mu/3*var.epsilon*hx;
    dg_ddu(1, 2, 0) = -var.mu/3*var.epsilon*h;
  }

  
  GreenNaghdiStiffnessOperatorDG
  ::GreenNaghdiStiffnessOperatorDG(GreenNaghdiEquationDG& var_kdv) : var(var_kdv)
  {    
  }
    
  
  
  template<class Vector1, class Vector2>
  void GreenNaghdiStiffnessOperatorDG
  ::EvaluateNumericalFlux(int i, int k, const Vector1& Um,
                          const Vector1& Up, Vector2& flux) const
  {
    flux.Zero();
    Real_wp hm = var.Csquare(i) + var.epsilon*Um(0);
    Real_wp hp = var.Csquare(i) + var.epsilon*Up(0);
    flux(0) = 0.5*(hm*Um(1) + hp*Up(1));
    flux(1) = 0.5*(Um(0) + 0.5*var.epsilon*Um(1)*Um(1)
                   - var.mu*var.epsilon/3*hm*hm*(Um(1)*Um(3) - Um(2)*Um(2))
                   + Up(0) + 0.5*var.epsilon*Up(1)*Up(1)
                   - var.mu*var.epsilon/3*hp*hp*(Up(1)*Up(3) - Up(2)*Up(2)) );
    
    flux(0) += 0.5*var.alpha_lax*(Um(0)-Up(0));
    flux(1) += 0.5*var.alpha_lax*(Um(1)-Up(1));    
    
    flux(2) = -0.5*(Um(1) + Up(1));
    flux(3) = -0.5*(Um(2) + Up(2));
    
    if (k == 0)
      Mlt(-1.0, flux);
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void GreenNaghdiStiffnessOperatorDG
  ::EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
                                    const Vector1& Up, Vector2& flux,
                                    Matrix1& dflux_dum, Matrix1& dflux_dup) const
  {
    flux.Zero(); dflux_dum.Zero(); dflux_dup.Zero();
    Real_wp hm = var.Csquare(i) + var.epsilon*Um(0);
    Real_wp hp = var.Csquare(i) + var.epsilon*Up(0);
    flux(0) = 0.5*(hm*Um(1) + hp*Up(1));
    flux(1) = 0.5*(Um(0) + 0.5*var.epsilon*Um(1)*Um(1) 
                   - var.mu*var.epsilon/3*hm*hm*(Um(1)*Um(3) - Um(2)*Um(2))
                   + Up(0) + 0.5*var.epsilon*Up(1)*Up(1) 
                   - var.mu*var.epsilon/3*hp*hp*(Up(1)*Up(3) - Up(2)*Up(2)) );
    
    dflux_dum(0, 0) = 0.5*var.epsilon*Um(1);
    dflux_dum(0, 1) = 0.5*hm;
    
    dflux_dup(0, 0) = 0.5*var.epsilon*Up(1);
    dflux_dup(0, 1) = 0.5*hp;

    dflux_dum(1, 0) = 0.5 - var.mu*var.epsilon*var.epsilon/3*hm*(Um(1)*Um(3) - Um(2)*Um(2));
    dflux_dum(1, 1) = 0.5*var.epsilon*(Um(1) - var.mu/3*hm*hm*Um(3));
    dflux_dum(1, 2) = var.mu*var.epsilon/3*hm*hm*Um(2);
    dflux_dum(1, 3) = -0.5*var.mu*var.epsilon/3*hm*hm*Um(1);
    
    dflux_dup(1, 0) = 0.5 - var.mu*var.epsilon*var.epsilon/3*hp*(Up(1)*Up(3) - Up(2)*Up(2));
    dflux_dup(1, 1) = 0.5*var.epsilon*(Up(1) - var.mu/3*hp*hp*Up(3));
    dflux_dum(1, 2) = var.mu*var.epsilon/3*hp*hp*Up(2);
    dflux_dup(1, 3) = -0.5*var.mu*var.epsilon/3*hp*hp*Up(1);
    
    flux(2) = -0.5*(Um(1) + Up(1));
    flux(3) = -0.5*(Um(2) + Up(2));
    
    dflux_dum(2, 1) = -0.5;
    dflux_dum(3, 2) = -0.5;
    dflux_dup(2, 1) = -0.5;
    dflux_dup(3, 2) = -0.5;

    flux(0) += 0.5*var.alpha_lax*(Um(0)-Up(0));
    flux(1) += 0.5*var.alpha_lax*(Um(1)-Up(1));
    
    dflux_dum(0, 0) += 0.5*var.alpha_lax;
    dflux_dum(1, 1) += 0.5*var.alpha_lax;
    
    dflux_dup(0, 0) -= 0.5*var.alpha_lax;
    dflux_dup(1, 1) -= 0.5*var.alpha_lax;
    
    if (k == 0)
      {
        Mlt(-1.0, flux);
        Mlt(-1.0, dflux_dum); Mlt(-1.0, dflux_dup);
      }
  }
  
  
  template<class Vector1, class Vector2>
  void GreenNaghdiStiffnessOperatorDG
  ::EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                     Vector2& g, Vector2& hs) const
  {
    Real_wp h = var.Csquare(i) + var.epsilon*u(0);
    Real_wp hx = var.CsquareDx(i) + var.epsilon*du(0);
    
    g(0) = -h*u(1);
    g(1) = -u(0) + var.epsilon*( -0.5*u(1)*u(1) + var.mu*h*h/3*(u(1)*u(3) - u(2)*u(2)));
    g(2) = u(1);
    g(3) = u(2);
    
    hs(1) = -var.mu*var.epsilon/3*h*hx*(u(1)*u(3) - u(2)*u(2));
    hs(2) = u(2);
    hs(3) = u(3);
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void GreenNaghdiStiffnessOperatorDG::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                     Vector2& g, Vector2& hs,
                     Matrix1& dg_du, Matrix1& dhs_du, Matrix1& dg_ddu, Matrix1& dhs_ddu) const
  {
    Real_wp h = var.Csquare(i) + var.epsilon*u(0);
    Real_wp hx = var.CsquareDx(i) + var.epsilon*du(0);
    Real_wp val = u(1)*u(3) - u(2)*u(2);
    
    g(0) = -h*u(1);
    g(1) = -u(0) + var.epsilon*( -0.5*u(1)*u(1) + var.mu*h*h/3*val );
    g(2) = u(1);
    g(3) = u(2);
    
    dg_du(0, 0) = -var.epsilon*u(1);
    dg_du(0, 1) = -h;
    dg_du(1, 0) = -1.0 + 2.0*square(var.epsilon)*var.mu*h/3*val;
    dg_du(1, 1) = var.epsilon*( -u(1) + var.mu*h*h/3*u(3));
    dg_du(1, 2) = -2.0*var.epsilon*var.mu*h*h/3*u(2);
    dg_du(1, 3) = var.epsilon*var.mu/3*h*h*u(1);
    dg_du(2, 1) = 1.0;
    dg_du(3, 2) = 1.0;

    hs(1) = -var.mu*var.epsilon/3*h*hx*(u(1)*u(3) - u(2)*u(2));
    hs(2) = u(2);
    hs(3) = u(3);

    dhs_du(1, 0) = -var.mu*square(var.epsilon)*hx/3*val;
    dhs_ddu(1, 0) = -var.mu*square(var.epsilon)*h/3*val;
    
    dhs_du(1, 1) = -var.mu*var.epsilon/3*h*hx*u(3);
    dhs_du(1, 2) = 2.0*var.mu*var.epsilon/3*h*hx*u(2);
    dhs_du(1, 3) = -var.mu*var.epsilon/3*h*hx*u(1);
    
    dhs_du(2, 2) = 1.0;
    dhs_du(3, 3) = 1.0;
  }

  
  GreenNaghdiEquationDG::GreenNaghdiEquationDG() : var_mass(*this), var_stiff(*this)
  {
    alpha_lax = 3.0;
  }
  
  
  GreenNaghdiEquationDG::GreenNaghdiEquationDG(const GreenNaghdiEquationDG& var)
    : GreenNaghdiEquation(var), var_mass(*this), var_stiff(*this)
  {
    alpha_lax = var.alpha_lax;
  }
  
  
  GreenNaghdiMassOperatorFem
  ::GreenNaghdiMassOperatorFem(const GreenNaghdiEquationFem& var_eq) : var(var_eq)
  {
  }

  
  template<class Vector1, class Vector2>
  void GreenNaghdiMassOperatorFem::
  EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                   Vector2& f, Vector2& g, Vector2& fx, Vector2& gx) const
  {
    Real_wp h = var.Csquare(i) + var.epsilon*u(0);
    Real_wp hx = var.CsquareDx(i) + var.epsilon*du(0);
    
    g(0, 0) = 1.0;
    g(1, 1) = 1.0;
    
    fx(1, 1) = var.mu/3*h*h;
    gx(1, 1) = -var.mu/3*h*hx;
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void GreenNaghdiMassOperatorFem::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du, Vector2& f,
                     Vector2& g, Vector2& fx, Vector2& gx,
                     Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                     Matrix1& dfx_du, Matrix1& dgx_du, Matrix1& dfx_ddu, Matrix1& dgx_ddu) const
  {
    Real_wp h = var.Csquare(i) + var.epsilon*u(0);
    Real_wp hx = var.CsquareDx(i) + var.epsilon*du(0);
    
    g(0, 0) = 1.0;
    g(1, 1) = 1.0;
    
    fx(1, 1) = var.mu/3*h*h;
    gx(1, 1) = -var.mu/3*h*hx;
    
    dfx_du(1, 1, 0) = 2.0*var.mu*var.epsilon/3*h;
    
    dgx_du(1, 1, 0) = -var.mu*var.epsilon/3*hx;
    dgx_ddu(1, 1, 0) = -var.mu*var.epsilon/3*h;
  }

  GreenNaghdiStiffnessOperatorFem
  ::GreenNaghdiStiffnessOperatorFem(const GreenNaghdiEquationFem& var_eq) : var(var_eq)
  {
  }
  
  
  template<class Vector1, class Vector2>
  void GreenNaghdiStiffnessOperatorFem::
  EvaluateFunction(int i, const Vector1& u, const Vector1& du, Vector2& f, Vector2& g) const
  {
    Real_wp h = var.Csquare(i) + var.epsilon*u(0);
    Real_wp hx = var.CsquareDx(i) + var.epsilon*du(0);
    Real_wp coef = var.epsilon*var.mu/3;
    
    f(0) = -h*u(1);
    f(1) = -(u(0) + 0.5*var.epsilon*u(1)*u(1)) + coef*h*h*(u(1)*du(2) - du(1)*du(1));
    
    g(1) = -coef*hx*h*(u(1)*du(2) - du(1)*du(1));
    g(2) = u(2) - du(1);
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void GreenNaghdiStiffnessOperatorFem::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du, Vector2& f,
                     Vector2& g, Matrix1& df_du, Matrix1& dg_du,
                     Matrix1& df_ddu, Matrix1& dg_ddu) const
  {
    Real_wp h = var.Csquare(i) + var.epsilon*u(0);
    Real_wp hx = var.CsquareDx(i) + var.epsilon*du(0);
    Real_wp coef = var.epsilon*var.mu/3;
    Real_wp val = u(1)*du(2) - du(1)*du(1);
    
    f(0) = -h*u(1);
    f(1) = -(u(0) + 0.5*var.epsilon*u(1)*u(1)) + coef*h*h*val;
    
    df_du(0, 0) = -var.epsilon*u(1);
    df_du(0, 1) = -h;
    
    df_du(1, 0) = -1.0 + 2.0*coef*var.epsilon*h*val;
    df_du(1, 1) = -var.epsilon*u(1) + coef*h*h*du(2);
    
    df_ddu(1, 1) = -2.0*coef*h*h*du(1);
    df_ddu(1, 2) = coef*h*h*u(1);
    
    g(1) = -coef*hx*h*val;
    g(2) = u(2) - du(1);
    
    dg_du(1, 0) = -coef*hx*var.epsilon*val;
    dg_du(1, 1) = -coef*h*hx*du(2);
    
    dg_ddu(1, 0) = -coef*var.epsilon*h*val;
    dg_ddu(1, 1) = 2.0*coef*hx*h*du(1);
    dg_ddu(1, 2) = -coef*h*hx*u(1);
    
    dg_du(2, 2) = 1.0;
    dg_ddu(2, 1) = -1.0;
  }
  
    
  GreenNaghdiEquationFem::GreenNaghdiEquationFem() : var_mass(*this), var_stiff(*this)
  {
  }
  
  
  GreenNaghdiEquationFem::GreenNaghdiEquationFem(const GreenNaghdiEquationFem& A)
    : GreenNaghdiEquation(A), var_mass(*this), var_stiff(*this)
  {
  }
  
  
  inline bool GreenNaghdiEquationFem::LinearMassTerm()
  {
    return true;
  }

}

#define MONTJOIE_FILE_GREEN_NAGHDI_PROBLEM_CXX
#endif
