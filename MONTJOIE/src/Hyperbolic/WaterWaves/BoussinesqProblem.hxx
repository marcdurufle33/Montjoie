#ifndef MONTJOIE_FILE_BOUSSINESQ_PROBLEM_HXX

namespace Montjoie
{
  
  //! base class for Boussinesq equation
  /*!
    We are considering the following Boussinesq system :
    zeta_t + [h u]_x = 0
    u_t + zeta_x + epsilon u u_x = mu/3 (c^4 u_{xt})_x
    where u_x is the partial derivative of u with respect to x
    and h = 1 + epsilon zeta - beta b^alpha
    alpha, beta, epsilon and mu are parameters of the model
   */
  class BoussinesqEquation : public KdvEquation
  {
  public :
    //! evaluation of c^2, d/dx(c^2), c^4 and d/dx(c^4)
    VectReal_wp Csquare, CsquareDx, C4Quad, C4xQuad;
    
    BoussinesqEquation();
    
    void ComputeInitialCondition(double, double, const VectReal_wp& x, VectReal_wp& z0);
    
    template<class GenericPb>
    void ComputePhysicalProperty(const GenericPb& var, double xmin, double xmax);
      
  };
    
  class BoussinesqEquationDG;

  
  //! class for stiffness operator for Boussinesq equation solved with LDG
  class BoussinesqStiffnessOperatorDG
  {
  protected :
    BoussinesqEquationDG& var;
    
  public :    
    BoussinesqStiffnessOperatorDG(BoussinesqEquationDG& var_kdv);
    
    template<class Vector1, class Vector2>
    void EvaluateNumericalFlux(int i, int k, const Vector1& Um,
			       const Vector1& Up, Vector2& flux) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
					 const Vector1& Up, Vector2& flux,
					 Matrix1& df_dum, Matrix1& df_dup) const;
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
			  Vector2& g, Vector2& h) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                            Vector2& g, Vector2& h, Matrix1& dg_du,
                            Matrix1& dh_du, Matrix1& dg_ddu, Matrix1& dh_ddu) const;
    
  };
  
  
  //! class for mass operator for Boussinesq equation solved with LDG
  class BoussinesqMassOperatorDG
  {
  protected :
    BoussinesqEquationDG& var;
    
  public :    
    BoussinesqMassOperatorDG(BoussinesqEquationDG& var_kdv);
    
    template<class Vector1, class Vector2>
    void EvaluateNumericalFlux(int i, int k, const Vector1& Um,
			       const Vector1& Up, Vector2& flux1, Vector2&flux2) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
					 const Vector1& Up, Vector2& flux1, Vector2& flux2,
					 Matrix1& df1_dum, Matrix1& df1_dup,
                                         Matrix1& df2_dum, Matrix1& df2_dup) const;
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                          Vector2& f, Vector2& g, Vector2& fx, Vector2& gx) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du, Vector2& f,
                            Vector2& g, Vector2& fx, Vector2& gx,
                            Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                            Matrix1& dfx_du, Matrix1& dgx_du,
                            Matrix1& dfx_ddu, Matrix1& dgx_ddu) const;
    
  };
  
  
  //! class for the resolution of Boussinesq equation with LDG method
  class BoussinesqEquationDG : public BoussinesqEquation
  {
  public :
    Real_wp alpha_lax;
    BoussinesqMassOperatorDG var_mass;
    BoussinesqStiffnessOperatorDG var_stiff;
    enum { nb_unknowns = 3, nb_time_unknowns = 2};

    static const bool FormulationDG = true;
    
    BoussinesqEquationDG();
    BoussinesqEquationDG(const BoussinesqEquationDG& var);
    
    static bool LinearMassTerm() { return true; }
    
  };

  
  class BoussinesqEquationFem;
  
  
  //! class for the mass operator for Boussinesq equation and continuous finite element
  class BoussinesqMassOperatorFem : public IdentityOperatorDg1D<BoussinesqEquationFem>
  {    
  protected :
    const BoussinesqEquationFem& var;
    
  public :    
    BoussinesqMassOperatorFem(const BoussinesqEquationFem& var_eq);
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                          Vector2& f, Vector2& g, Vector2& fx, Vector2& gx) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du, Vector2& f,
                            Vector2& g, Vector2& fx, Vector2& gx,
                            Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                            Matrix1& dfx_du, Matrix1& dgx_du,
                            Matrix1& dfx_ddu, Matrix1& dgx_ddu) const;
    
  };
  

  //! class for the stiffness operator for Boussinesq equation and continuous finite element  
  class BoussinesqStiffnessOperatorFem : public IdentityOperatorDg1D<BoussinesqEquationFem>
  {    
  protected :
    const BoussinesqEquationFem& var;
    
  public :
    
    BoussinesqStiffnessOperatorFem(const BoussinesqEquationFem& var_eq);
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                          Vector2& f, Vector2& g) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du, Vector2& f,
                            Vector2& g, Matrix1& df_du, Matrix1& dg_du,
                            Matrix1& df_ddu, Matrix1& dg_ddu) const;
    
  };

  
  //! class for the resolution of Boussinesq equation with continuous finite element
  class BoussinesqEquationFem : public BoussinesqEquation
  {
  public :
    BoussinesqMassOperatorFem var_mass;
    BoussinesqStiffnessOperatorFem var_stiff;
    enum { nb_unknowns = 2, nb_time_unknowns = 2};
    
    static const bool FormulationDG = false;
    
    BoussinesqEquationFem();

    BoussinesqEquationFem(const BoussinesqEquationFem& A);

    static bool LinearMassTerm();
    
  };
  
}

#define MONTJOIE_FILE_BOUSSINESQ_PROBLEM_HXX
#endif
