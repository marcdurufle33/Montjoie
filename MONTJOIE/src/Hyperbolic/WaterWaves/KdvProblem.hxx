#ifndef MONTJOIE_FILE_KDV_PROBLEM_HXX

namespace Montjoie
{
    
  class KdvEquation : public GenericWaterWaveEquation
  {
  public :
    // evaluation of c, c_x, c^5, etc on quadrature points
    VectReal_wp Cquad, CxQuad, invCquad;
    VectReal_wp C5Quad, C5xQuad, C5xxQuad, C5xxxQuad;

    
    KdvEquation();
    
    template<class GenericPb>
    void ComputePhysicalProperty(const GenericPb& var, double xmin, double xmax);

  };
  
  
  class KdvEquationDG;
  
  class KdvStiffnessOperatorDG
  {
  protected :
    KdvEquationDG& var;
    
  public :    
    KdvStiffnessOperatorDG(KdvEquationDG& var_kdv);
    
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
  
  
  class KdvEquationDG : public KdvEquation
  {
  public :
    Real_wp alpha_lax;
    IdentityOperatorDg1D<KdvEquationDG> var_mass;
    KdvStiffnessOperatorDG var_stiff;
    enum { nb_unknowns = 3, nb_time_unknowns = 1};
    
    KdvEquationDG();
    KdvEquationDG(const KdvEquationDG& var);
    
    static const bool FormulationDG = true;
    static bool LinearMassTerm() { return false; }
        
  };

  class KdvEquationFem;
  
  class KdvStiffnessOperatorFem : public IdentityOperatorDg1D<KdvEquationFem>
  {    
  protected :
    const KdvEquationFem& var;
    
  public :
    enum{nb_terms = 1};
    
    KdvStiffnessOperatorFem(const KdvEquationFem& var_eq);
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
			  Vector2& g, Vector2& h) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                            Vector2& g, Vector2& h, Matrix1& df_du,
                            Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu) const;
    
  };


  class KdvEquationFem : public KdvEquation
  {
  public :
    IdentityOperatorDg1D<KdvEquationFem> var_mass;
    KdvStiffnessOperatorFem var_stiff;
    enum { nb_unknowns = 3, nb_time_unknowns = 1};
    
    static const bool FormulationDG = false;
    
    KdvEquationFem();

    KdvEquationFem(const KdvEquationFem& A);

    static bool LinearMassTerm();
    
  };
  
  class KdvEquationSpectral : public KdvEquation
  {
  public :
    enum { nb_unknowns = 1, nb_time_unknowns = 1};
    typedef KdvEquation TypeEquationDiff;
    
  protected :    
    bool constant_coef;
    VectReal_wp XnReel, dXreal, ProdReal, d2Xreal, d3Xreal;
    VectComplex_wp dXn, d2Xn, d3Xn, XnHat;
    VectComplex_wp coef_Gamma1, coef_Gamma2, coef_Gamma3;
    
  public :
    KdvEquationSpectral();
    
    void ComputePhysicalProperty(FftRealInterface& fft, const VectReal_wp& nu,
				 double xmin, double xmax, int N);
    
    void InitSinus(const Real_wp& dt, const VectReal_wp& nu);
    
    void ApplyMass(const VectReal_wp&,
                   FftRealInterface& fft, VectReal_wp& nu, VectReal_wp&);
    
    void AddFunction(const Real_wp&, const VectReal_wp&,
                     FftRealInterface& fft, VectReal_wp& nu,
                     VectReal_wp&, bool solve_mass = true);
    
    void EvaluateFunction(const Real_wp& t, const VectComplex_wp& X,
			  FftRealInterface& fft, VectReal_wp& nu,
			  VectComplex_wp& Y, bool solve_mass = true);

    void EvaluateFunction(const Real_wp& t, const VectReal_wp& X,
			  FftRealInterface& fft, VectReal_wp& nu,
			  VectReal_wp& Y, bool solve_mass = true);
    
    void GetFunction(const VectComplex_wp&, const VectReal_wp&,
                     FftRealInterface&, VectReal_wp&, VectReal_wp& Y);
    
  };
  
}

#define MONTJOIE_FILE_KDV_PROBLEM_HXX
#endif
