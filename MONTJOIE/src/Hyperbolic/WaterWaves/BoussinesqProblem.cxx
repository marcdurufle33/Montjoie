#ifndef MONTJOIE_FILE_BOUSSINESQ_PROBLEM_CXX

namespace Montjoie
{
  //! default constructor
  BoussinesqEquation::BoussinesqEquation() : KdvEquation()
  {
  }
    
  
  //! computation of the initial condition
  void BoussinesqEquation
  ::ComputeInitialCondition(double xmin, double xmax, const VectReal_wp& Xdof, VectReal_wp& Y0)
  {
    // computation of b and c
    int nodl = Xdof.GetM();
    VectReal_wp Bdof, Cdof;
    ComputeBottom(xmin, xmax, Xdof, Bdof, Cdof);
    
    // computing initial condition
    VectReal_wp Z0, Z0_x(nodl), Z0_xx(nodl);
    KdvEquation::ComputeInitialCondition(xmin, xmax, Xdof, Z0, Z0_x, Z0_xx);
    Y0.Reallocate(2*nodl); Y0.Fill(0);
    for (int i = 0; i < nodl; i++)
      Y0(i) = Z0(i);
    
    // compatible initial condition with KdV
    for (int i = 0; i < nodl; i++)
      {
	Real_wp c = Cdof(i), z = Z0(i);
        Y0(nodl+i) = z/c;
	//Y0(nodl+i) = (z - 0.25*epsilon/(c*c)*z*z + mu/6*pow(c, 4.0)*Z0_xx(i)) / c;
      }
  }
  
  
  //! computation of physical coefficients needed to evaluate the scheme
  template<class GenericPb>
  void BoussinesqEquation::ComputePhysicalProperty(const GenericPb& var, double xmin, double xmax)
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
    
    // computation of C^2 and C^4 on nodal points
    VectReal_wp Bhigh, Chigh, Cs, C4, C4xDof;
    ComputeBottom(xmin, xmax, Xhigh, Bhigh, Chigh);
    Cs = Chigh; C4 = Chigh;
    for (int i = 0; i < Cs.GetM(); i++)
      {
	Cs(i) = square(Chigh(i));
	C4(i) = square(Cs(i));
      }
    
    // DISP(Cs); DISP(Chigh);
    // computation of c^2 and (c^2)_x on quadrature points
    VectReal_wp CxDof;
    GetGradient(var, Cs, CxDof);
    GetGradient(var, C4, C4xDof);
    
    GetUquadrature(var, Cs, Csquare);
    GetUquadrature(var, CxDof, CsquareDx);
    GetUquadrature(var, C4, C4Quad);
    GetUquadrature(var, C4xDof, C4xQuad);
  }


  /***********************
   * Boussinesq with LDG *
   ***********************/
  
  
  //! default constructor
  BoussinesqMassOperatorDG::BoussinesqMassOperatorDG(BoussinesqEquationDG& var_kdv) : var(var_kdv)
  {
  }
    
  
  //! numerical fluxes for mass operator
  template<class Vector1, class Vector2>
  void BoussinesqMassOperatorDG
  ::EvaluateNumericalFlux(int i, int k, const Vector1& Um,
                          const Vector1& Up, Vector2& flux1, Vector2& flux2) const
  {
    flux1.Zero(); flux2.Zero();
    if (k == 0)
      flux1(1, 2) = var.mu/3;
    else 
      flux1(1, 2) = -var.mu/3;
  }
  
  
  //! numerical fluxes for mass operator
  template<class Vector1, class Vector2, class Matrix1>
  void BoussinesqMassOperatorDG
  ::EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
                                    const Vector1& Up, Vector2& flux1, Vector2& flux2,
                                    Matrix1& df1_dum, Matrix1& df1_dup,
                                    Matrix1& df2_dum, Matrix1& df2_dup) const
  {
    flux1.Zero(); flux2.Zero();
    if (k == 0)
      flux1(1, 2) = var.mu/3;
    else 
      flux1(1, 2) = -var.mu/3;
  }
  
  
  //! volume integrals for mass operator
  template<class Vector1, class Vector2>
  void BoussinesqMassOperatorDG
  ::EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                     Vector2& f, Vector2& g, Vector2& fx, Vector2& gx) const
  {
    f(1, 2) = var.mu/3;
    g(0, 0) = 1.0;
    g(1, 1) = 1.0;
  }
  
  
  //! volume integrals for mass operator
  template<class Vector1, class Vector2, class Matrix1>
  void BoussinesqMassOperatorDG
  ::EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                       Vector2& f, Vector2& g, Vector2& fx, Vector2& gx,
                       Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                       Matrix1& dfx_du, Matrix1& dgx_du, Matrix1& dfx_ddu, Matrix1& dgx_ddu) const
  {
    f(1, 2) = var.mu/3;
    g(0, 0) = 1.0;
    g(1, 1) = 1.0;
  }
  
  
  BoussinesqStiffnessOperatorDG
  ::BoussinesqStiffnessOperatorDG(BoussinesqEquationDG& var_kdv) : var(var_kdv)
  {
  }
    
  
  //! numerical fluxes for stiffness operator
  template<class Vector1, class Vector2>
  void BoussinesqStiffnessOperatorDG
  ::EvaluateNumericalFlux(int i, int k, const Vector1& Um, const Vector1& Up, Vector2& flux) const
  {
    flux.Zero();
    if (k == 1)
      flux(2) = -var.C4Quad(i)*(Up(1) - Um(1));
    
    R2 flux_m, flux_p;
    flux_m(0) = (var.Csquare(i) + var.epsilon*Um(0))*Um(1);
    flux_m(1) = Um(0) + 0.5*var.epsilon*Um(1)*Um(1);

    flux_p(0) = (var.Csquare(i) + var.epsilon*Up(0))*Up(1);
    flux_p(1) = Up(0) + 0.5*var.epsilon*Up(1)*Up(1);
    
    flux(0) = 0.5*(flux_p(0) + flux_m(0));
    flux(1) = 0.5*(flux_p(1) + flux_m(1));
    if (k == 0)
      {
	flux(0) = -flux(0) - 0.5*var.alpha_lax*(Um(0)-Up(0));
	flux(1) = -flux(1) - 0.5*var.alpha_lax*(Um(1)-Up(1));
      }
    else
      {
	flux(0) += 0.5*var.alpha_lax*(Um(0)-Up(0));
	flux(1) += 0.5*var.alpha_lax*(Um(1)-Up(1));
      }
  }
  
  
  //! numerical fluxes for stiffness operator
  template<class Vector1, class Vector2, class Matrix1>
  void BoussinesqStiffnessOperatorDG
  ::EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
                                    const Vector1& Up, Vector2& flux,
                                    Matrix1& df_dum, Matrix1& df_dup) const
  {
    flux.Zero(); df_dum.Zero(); df_dup.Zero();
    if (k == 1)
      {
        flux(2) = -var.C4Quad(i)*(Up(1) - Um(1));
        df_dum(2, 1) = var.C4Quad(i);
        df_dup(2, 1) = -var.C4Quad(i);
      }
    
    R2 flux_m, flux_p;
    flux_m(0) = (var.Csquare(i) + var.epsilon*Um(0))*Um(1);
    flux_m(1) = Um(0) + 0.5*var.epsilon*Um(1)*Um(1);

    flux_p(0) = (var.Csquare(i) + var.epsilon*Up(0))*Up(1);
    flux_p(1) = Up(0) + 0.5*var.epsilon*Up(1)*Up(1);
    
    flux(0) = 0.5*(flux_p(0) + flux_m(0));
    flux(1) = 0.5*(flux_p(1) + flux_m(1));
    df_dum(0, 0) = 0.5*var.epsilon*Um(1);
    df_dum(0, 1) = 0.5*(var.Csquare(i) + var.epsilon*Um(0));
    df_dum(1, 0) = 0.5;
    df_dum(1, 1) = 0.5*var.epsilon*Um(1);
    
    df_dup(0, 0) = 0.5*var.epsilon*Up(1);
    df_dup(0, 1) = 0.5*(var.Csquare(i) + var.epsilon*Up(0));
    df_dup(1, 0) = 0.5;
    df_dup(1, 1) = 0.5*var.epsilon*Up(1);
    if (k == 0)
      {
	flux(0) = -flux(0) - 0.5*var.alpha_lax*(Um(0)-Up(0));
	flux(1) = -flux(1) - 0.5*var.alpha_lax*(Um(1)-Up(1));
        df_dum(0, 0) = -df_dum(0, 0) - 0.5*var.alpha_lax;
        df_dum(0, 1) = -df_dum(0, 1);
        df_dum(1, 0) = -df_dum(1, 0);
        df_dum(1, 1) = -df_dum(1, 1) - 0.5*var.alpha_lax;
        
        df_dup(0, 0) = -df_dup(0, 0) + 0.5*var.alpha_lax;
        df_dup(0, 1) = -df_dup(0, 1);
        df_dup(1, 0) = -df_dup(1, 0);
        df_dup(1, 1) = -df_dup(1, 1) + 0.5*var.alpha_lax;
      }
    else
      {
	flux(0) += 0.5*var.alpha_lax*(Um(0)-Up(0));
	flux(1) += 0.5*var.alpha_lax*(Um(1)-Up(1));

        df_dum(0, 0) += 0.5*var.alpha_lax;
        df_dum(1, 1) += 0.5*var.alpha_lax;

        df_dup(0, 0) -= 0.5*var.alpha_lax;
        df_dup(1, 1) -= 0.5*var.alpha_lax;
      }
  }
  

  //! volume integrals for stiffness operator
  template<class Vector1, class Vector2>
  void BoussinesqStiffnessOperatorDG
  ::EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                     Vector2& g, Vector2& h) const
  {
    g(0) = -(var.Csquare(i) + var.epsilon*u(0))*u(1);
    g(1) = -u(0) - 0.5*var.epsilon*u(1)*u(1);

    // equation p - u_x = 0
    h(2) = u(2) -var.C4Quad(i)*du(1);
  }
  
  
  //! volume integrals for stiffness operator
  template<class Vector1, class Vector2, class Matrix1>
  void BoussinesqStiffnessOperatorDG::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                     Vector2& g, Vector2& h,
                     Matrix1& dg_du, Matrix1& dh_du, Matrix1& dg_ddu, Matrix1& dh_ddu) const
  {
    g(0) = -(var.Csquare(i) + var.epsilon*u(0))*u(1);
    g(1) = -u(0) - 0.5*var.epsilon*u(1)*u(1);
    
    dg_du(0, 0) = -var.epsilon*u(1);
    dg_du(0, 1) = -(var.Csquare(i) + var.epsilon*u(0));
    dg_du(1, 0) = -1.0;
    dg_du(1, 1) = -var.epsilon*u(1);
    
    // equation p - u_x = 0
    h(2) = u(2) -var.C4Quad(i)*du(1);
    
    dh_du(2, 2) = 1.0;
    dh_ddu(2, 1) = -var.C4Quad(i);
  }

  
  //! default constructor
  BoussinesqEquationDG::BoussinesqEquationDG() : var_mass(*this), var_stiff(*this)
  {
    alpha_lax = 3.0;
  }
  
  
  //! copy constructor
  BoussinesqEquationDG::BoussinesqEquationDG(const BoussinesqEquationDG& var)
    : BoussinesqEquation(var), var_mass(*this), var_stiff(*this)
  {
    alpha_lax = var.alpha_lax;
  }
  

  /**********************************
   * Boussinesq with finite element *
   **********************************/

  
  BoussinesqMassOperatorFem::BoussinesqMassOperatorFem(const BoussinesqEquationFem& var_eq)
    : var(var_eq)
  {
  }

  
  //! volume integrals for mass operator
  template<class Vector1, class Vector2>
  void BoussinesqMassOperatorFem::
  EvaluateFunction(int i, const Vector1& u, const Vector1& du, Vector2& f,
                   Vector2& g, Vector2& fx, Vector2& gx) const
  {
    g(0, 0) = 1.0;
    g(1, 1) = 1.0;
    fx(1, 1) = var.mu/3*var.C4Quad(i);
  }
  
  
  //! volume integrals for mass operator
  template<class Vector1, class Vector2, class Matrix1>
  void BoussinesqMassOperatorFem::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du, Vector2& f,
                     Vector2& g, Vector2& fx, Vector2& gx,
                     Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                     Matrix1& dfx_du, Matrix1& dgx_du, Matrix1& dfx_ddu, Matrix1& dgx_ddu) const
  {
    g(0, 0) = 1.0;
    g(1, 1) = 1.0;
    fx(1, 1) = var.mu/3*var.C4Quad(i);
  }
  
  
  //! default constructor
  BoussinesqStiffnessOperatorFem
  ::BoussinesqStiffnessOperatorFem(const BoussinesqEquationFem& var_eq) : var(var_eq)
  {
  }
  
  
  //! volume integrals for stiffness operator
  template<class Vector1, class Vector2>
  void BoussinesqStiffnessOperatorFem::
  EvaluateFunction(int i, const Vector1& u, const Vector1& du, Vector2& f, Vector2& g) const
  {
    f(0) = -(var.Csquare(i) + var.epsilon*u(0))*u(1);
    f(1) = -(u(0) + 0.5*var.epsilon*u(1)*u(1));
  }
  
  
  //! volume integrals for stiffness operator
  template<class Vector1, class Vector2, class Matrix1>
  void BoussinesqStiffnessOperatorFem::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du, Vector2& f,
                     Vector2& g, Matrix1& df_du, Matrix1& dg_du,
                     Matrix1& df_ddu, Matrix1& dg_ddu) const
  {
    f(0) = -(var.Csquare(i) + var.epsilon*u(0))*u(1);
    f(1) = -(u(0) + 0.5*var.epsilon*u(1)*u(1));
    
    df_du(0, 0) = -var.epsilon*u(1);
    df_du(0, 1) = -(var.Csquare(i) + var.epsilon*u(0));
    
    df_du(1, 0) = -1.0;
    df_du(1, 1) = -var.epsilon*u(1);
  }
  
  
  //! default constructor
  BoussinesqEquationFem::BoussinesqEquationFem() : var_mass(*this), var_stiff(*this)
  {
  }
  
  
  //! copy constructor
  BoussinesqEquationFem::BoussinesqEquationFem(const BoussinesqEquationFem& A)
    : BoussinesqEquation(A), var_mass(*this), var_stiff(*this)
  {
  }

  
  //! mass operator is linear
  inline bool BoussinesqEquationFem::LinearMassTerm()
  {
    return true;
  }
  
}

#define MONTJOIE_FILE_BOUSSINESQ_PROBLEM_CXX
#endif
