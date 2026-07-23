#ifndef MONTJOIE_FILE_GREEN_NAGHDI_PROBLEM_HXX

namespace Montjoie
{
    
  class GreenNaghdiEquation : public CamassaHolmEquation
  {
  public :
    VectReal_wp Csquare, CsquareDx;
    
    GreenNaghdiEquation();
    
    void ComputeInitialCondition(double, double, const VectReal_wp& x, VectReal_wp& z0);
    
    template<class GenericPb>
    void ComputePhysicalProperty(const GenericPb& var, double, double);
    
  };
  
  
  class GreenNaghdiEquationDG;
  
  class GreenNaghdiStiffnessOperatorDG
  {
  protected :
    GreenNaghdiEquationDG& var;
    
  public :    
    GreenNaghdiStiffnessOperatorDG(GreenNaghdiEquationDG& var_kdv);
        
    template<class Vector1, class Vector2>
    void EvaluateNumericalFlux(int i, int k, const Vector1& Um,
			       const Vector1& Up, Vector2& flux) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
					 const Vector1& Up, Vector2& flux,
					 Matrix1& df_dum, Matrix1& df_dup) const;
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
			  Vector2& f, Vector2& g) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
			    Vector2& f, Vector2& g,
			    Matrix1& df_du, Matrix1& dg_du,
                            Matrix1& df_ddu, Matrix1& dg_ddu) const;
    
  };
  
  
  class GreenNaghdiMassOperatorDG
  {
  protected :
    GreenNaghdiEquationDG& var;
    
  public :    
    GreenNaghdiMassOperatorDG(GreenNaghdiEquationDG& var_kdv);
    
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
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                            Vector2& f, Vector2& g, Vector2& fx, Vector2& gx,
                            Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                            Matrix1& dfx_du, Matrix1& dgx_du,
                            Matrix1& dfx_ddu, Matrix1& dgx_ddu) const;
    
  };
  
  
  class GreenNaghdiEquationDG : public GreenNaghdiEquation
  {
  public :
    Real_wp alpha_lax;
    GreenNaghdiMassOperatorDG var_mass;
    GreenNaghdiStiffnessOperatorDG var_stiff;
    enum { nb_unknowns = 4, nb_time_unknowns = 2};

    static const bool FormulationDG = true;
    
    GreenNaghdiEquationDG();
    GreenNaghdiEquationDG(const GreenNaghdiEquationDG& var);
    
    static bool LinearMassTerm() { return true; }
    
  };
  
  class GreenNaghdiEquationFem;
  
  class GreenNaghdiMassOperatorFem : public IdentityOperatorDg1D<GreenNaghdiEquationFem>
  {    
  protected :
    const GreenNaghdiEquationFem& var;
    
  public :
    GreenNaghdiMassOperatorFem(const GreenNaghdiEquationFem& var_eq);
    
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
  
  
  class GreenNaghdiStiffnessOperatorFem : public IdentityOperatorDg1D<GreenNaghdiEquationFem>
  {    
  protected :
    const GreenNaghdiEquationFem& var;
    
  public :
    GreenNaghdiStiffnessOperatorFem(const GreenNaghdiEquationFem& var_eq);
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                          Vector2& f, Vector2& g) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du, Vector2& f,
                            Vector2& g, Matrix1& df_du, Matrix1& dg_du,
                            Matrix1& df_ddu, Matrix1& dg_ddu) const;
  };


  class GreenNaghdiEquationFem : public GreenNaghdiEquation
  {
  public :
    GreenNaghdiMassOperatorFem var_mass;
    GreenNaghdiStiffnessOperatorFem var_stiff;
    enum { nb_unknowns = 3, nb_time_unknowns = 2};
    
    static const bool FormulationDG = false;
    
    GreenNaghdiEquationFem();

    GreenNaghdiEquationFem(const GreenNaghdiEquationFem& A);

    static bool LinearMassTerm();
    
  };

}

#define MONTJOIE_FILE_GREEN_NAGHDI_PROBLEM_HXX
#endif
