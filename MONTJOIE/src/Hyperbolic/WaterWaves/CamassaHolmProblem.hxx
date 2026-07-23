#ifndef MONTJOIE_FILE_CAMASSA_HOLM_PROBLEM_HXX

namespace Montjoie
{
    
  class CamassaHolmEquation : public GenericWaterWaveEquation
  {
  public :    
    VectReal_wp Cquad, CxQuad, invCquad, invC2quad;
    // for violent model
    VectReal_wp CoefZeta, CoefDxZeta, CoefDxxZeta, CoefDxxxZeta;
    // for original model
    VectReal_wp CoefEtilde, CoefAtilde, CoefFtilde, CoefFxTilde, CoefFxxTilde;
    
    CamassaHolmEquation();
    
    template<class GenericPb>
    void ComputePhysicalProperty(const GenericPb& var, double xmin, double xmax);
    
  };
  
  
  class CamassaHolmEquationDG;
  
  class CamassaHolmMassOperatorDG
  {
  protected :
    CamassaHolmEquationDG& var;
    
  public :    
    CamassaHolmMassOperatorDG(CamassaHolmEquationDG& var_kdv);
    
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
  
  
  class CamassaHolmStiffnessOperatorDG
  {
  protected :
    CamassaHolmEquationDG& var;
    
  public :    
    CamassaHolmStiffnessOperatorDG(CamassaHolmEquationDG& var_kdv);
    
    template<class Vector1, class Vector2>
    void EvaluateFunctionFlux(int i, const Vector1& Um, Vector2& flux) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivativeFlux(int i, const Vector1& Um, Vector2& flux, Matrix1& dflux_du) const;
    
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
			    Vector2& g, Vector2& h, Matrix1& df_du,
                            Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu) const;
    
  };
  
  
  class CamassaHolmEquationDG : public CamassaHolmEquation
  {
  public :  
    Real_wp alpha_lax;
    CamassaHolmMassOperatorDG var_mass;
    CamassaHolmStiffnessOperatorDG var_stiff;
    enum { nb_unknowns = 3, nb_time_unknowns = 1};
    
    static const bool FormulationDG = true;
    
    CamassaHolmEquationDG();
    CamassaHolmEquationDG(const CamassaHolmEquationDG& var);
    
    static bool LinearMassTerm() { return true; }
    
  };
  
  class CamassaHolmEquationFem;
  
  class CamassaHolmMassOperatorFem : public IdentityOperatorDg1D<CamassaHolmEquationFem>
  {    
  protected :
    const CamassaHolmEquationFem& var;
    
  public :
    CamassaHolmMassOperatorFem(const CamassaHolmEquationFem& var_eq);
    
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
  
  class CamassaHolmStiffnessOperatorFem : public IdentityOperatorDg1D<CamassaHolmEquationFem>
  {    
  protected :
    const CamassaHolmEquationFem& var;
    
  public :
        
    CamassaHolmStiffnessOperatorFem(const CamassaHolmEquationFem& var_eq);
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
			  Vector2& g, Vector2& h) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
			    Vector2& g, Vector2& h, Matrix1& df_du,
                            Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu) const;
    
  };


  class CamassaHolmEquationFem : public CamassaHolmEquation
  {
  public :
    CamassaHolmMassOperatorFem var_mass;
    CamassaHolmStiffnessOperatorFem var_stiff;
    enum { nb_unknowns = 3, nb_time_unknowns = 1};
    
    static const bool FormulationDG = false;
    
    CamassaHolmEquationFem();

    CamassaHolmEquationFem(const CamassaHolmEquationFem& A);

    static bool LinearMassTerm();
    
  };


  class CamassaHolmEquationSpectral : public CamassaHolmEquation
  {
  public :
    enum { nb_unknowns = 1, nb_time_unknowns = 1};
    typedef CamassaHolmEquation TypeEquationDiff;
    
  protected :
    VectReal_wp XnReel, dXreal, ProdReal, d2Xreal, d3Xreal;
    VectComplex_wp XnHat, dXn, d2Xn, d3Xn, ProdHat;
    VectReal_wp CoefMass;
    VectComplex_wp coef_Gamma1, coef_Gamma2, coef_Gamma3;
    
  public :
    CamassaHolmEquationSpectral();
    
    void ComputePhysicalProperty(FftRealInterface& fft, const VectReal_wp& nu,
				 double xmin, double xmax, int N);

    void InitSinus(const Real_wp& dt, const VectReal_wp& nu);

    void ApplyMass(const VectReal_wp&,
                   FftRealInterface& fft, VectReal_wp& nu, VectReal_wp&);
    
    void AddFunction(const Real_wp&, const VectReal_wp&,
                     FftRealInterface& fft, VectReal_wp& nu,
                     VectReal_wp&, bool solve_mass = true);
        
    void GetFunction(const VectComplex_wp& Xn, const VectReal_wp& Xreal,
                     FftRealInterface& fft, VectReal_wp& nu,
                     VectComplex_wp& ProdXn, bool solve_mass = true);
    
    void EvaluateFunction(const Real_wp& t, const VectComplex_wp& X,
			  FftRealInterface& fft, VectReal_wp& nu,
			  VectComplex_wp& Y, bool solve_mass = true);

    void EvaluateFunction(const Real_wp& t, const VectReal_wp& X,
			  FftRealInterface& fft, VectReal_wp& nu,
			  VectReal_wp& Y, bool solve_mass = true);
    
  };

}

#define MONTJOIE_FILE_CAMASSA_HOLM_PROBLEM_HXX
#endif
