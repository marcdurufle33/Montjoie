#ifndef MONTJOIE_FILE_LOCAL_DISCONTINUOUS_GALERKIN_1D_HXX

namespace Montjoie
{
  
  //! base class for LDG method in 1-D
  template<class TypeEquation>
  class LdgScheme1D : public VirtualOdeSystem<Real_wp>, public InputDataProblem_Base
  {
  public :
    //! 1-D mesh
    Mesh<Dimension1> mesh;
    enum {CENTERED, LEFT_UPWIND, RIGHT_UPWIND, LAX_FRIEDRICHS};
    
    //! object containing datas of the solved equation
    TypeEquation var_eq;
    //! use of Gauss-Lobatto points for degrees of freedom
    Globatto<Real_wp> lob;
    //! integration with Gauss points
    Globatto<Real_wp> gauss;
    //! values of phi and dphi/dx on quadrature points
    Matrix<Real_wp> GradPhi, ValPhi, ValPhiWeight, GradPhiWeight, G_GL;
    //! h_i for each element, and 1/h_i
    VectReal_wp h_subdiv, invH_subdiv;
    //! position of degrees of freedom
    VectReal_wp Xdof;
    //! display grid
    GridInterpolation<Dimension1> grid;
    //! true if we use mass lumping for mass matrix
    bool mass_lumping;
    //! diagonal mass matrix with mass lumping
    VectReal_wp invDh, Dh;
    //! dense mass matrix without mass lumping
    Matrix<Real_wp> invMass, mat_mass, mat_rigid;
    //! 
    bool periodic_condition;
    VectReal_wp invGlobalMass;
    
  public :
    LdgScheme1D();

    void SetInputData(const string& keyword, const Vector<string>& param);
    
    int GetQuadNumber(int i, int j) const;
    void SetOrder(int order);
    void SetInterval(const Real_wp& x0, const Real_wp& xN, int n, int r);
    
    void InitGrid(const Real_wp& xmin, const Real_wp& xmax, int Nd);
      
    template<class T>
    void SolveMass(int i, Vector<T>& Zn) const;
    
    template<class T0, class T1>
    void MltMass(int i, const T0& alpha, const Vector<T1>& Un, Vector<T1>& Zn) const;
    
    template<class T>
    void GetFunction(int i, const Vector<T>& Zn, Vector<T>& ZnQuad) const;
    
    template<class T0, class T1>
    void IntegrateAgainstFunction(int i, const T0& alpha,
                                  const Vector<T1>& ZnQuad, Vector<T1>& Zn) const;
    
    template<class T0>
    void GetDerivative(int i, const Vector<T0>& Zn, Vector<T0>& dZn) const;
    
    template<class T0, class T1>
    void IntegrateAgainstDerivative(int i, const T0& alpha,
                                    const Vector<T1>& dZn, Vector<T1>& Zn) const;
    
    template<class T>
    void ProjectToDof(const Vector<T>&, Vector<T>& ) const;
    
    template<class T>
    void GetUplusMinus(int i, int offset, const Vector<T>& Un, const IVect&, T& Um, T& Up) const;
    
    template<int m, class T>
    void GetUplusMinus(int i, int offset, const Vector<T>& Un, const IVect&,
		       TinyVector<T, m>& Um, TinyVector<T, m>& Up) const;
    
    template<class OperatorEq, class T0, class T1>
    void MltAddStiffness(const T0& alpha, const Vector<T1>& Un,
                         int offset, const OperatorEq& var, const IVect&,
			 Vector<T1>& Zn, bool solve_mass = false) const;
    
    template<class OperatorEq, class T0, class T1>
    void MltAddMass(const T0& alpha, const Vector<T1>& Un, const Vector<T1>& dU_dt,
                    int offset, const OperatorEq& var, const IVect&, Vector<T1>& Zn) const;
    
    template<class T0, class T1>
    void AddSystem(const Vector<T1>& Un_mass, const Vector<T1>& Un, int offset,
		   const T0& alpha_mass, const T0& beta_stiff, const IVect&, Vector<T1>& Fn) const;
    
    template<class OperatorEq, class MatrixSparse, class T0, class T1>
    void AddStiffnessDifferential(const Vector<T0>& alpha, const Vector<T1>& Un,
                                  int offset, const OperatorEq& var,
                                  const IVect&, MatrixSparse& DF) const;

    template<class OperatorEq, class MatrixSparse, class T0, class T1>
    void AddMassDifferential(const Vector<T0>& alpha, const Vector<T1>& Un,
                             const Vector<T1>& dU_dt,
			     int offset, const OperatorEq& var,
                             const IVect&, MatrixSparse& DF) const;
    
    template<class MatrixSparse, class T0, class T1>
    void GetDifferential(const Vector<T1>& Un_mass, const Vector<T1>& Un, int offset,
			 const Vector<T0>& beta_stiff, const IVect&, MatrixSparse& DF) const;
    
    Real_wp GetEnergy(const VectReal_wp& Zn, const IVect&) const;
    Real_wp GetEnergy(const VectComplex_wp& Zn, const IVect&) const;
    
    template<class T>
    void GetInterpolateUn(int m, const Vector<T>& Zn, const IVect&, Vector<T>& ZnInterp) const;

    void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp&, Real_wp alpha = 1.0) {}
    
  };
 

  template<class TypeEquation>
  class TimeLdgScheme1D : public LdgScheme1D<TypeEquation>
  {
  public :
    LowStorageRK_Iterator<Real_wp> RK_scheme;
    int type_scheme, order_time_scheme;
    Real_wp deltat; VectReal_wp Un_time;
    Real_wp stopping_criterion;
    Matrix<Real_wp> AgaussRK; VectReal_wp BgaussRK;
    Matrix<Real_wp, General, ArrayRowSparse> mat_diff;
    bool compute_differential;
    int nb_eval_diff;
    IVect NewRowNumber, GaussRowNumber;
    //TinyArrowMatrix<double, 15, 5> mat_lu;
    //TinyArrowMatrix<double, 15, 5> mat_diff;
    
#ifdef SELDON_WITH_PASTIX
    typedef MatrixPastix<Real_wp> MatrixLU;
#else
#ifdef SELDON_WITH_MUMPS
    typedef MatrixMumps<Real_wp> MatrixLU;
#else
#ifdef SELDON_WITH_UMFPACK
    typedef MatrixUmfPack<double> MatrixLU;
#else
#ifdef SELDON_WITH_SUPERLU
    typedef MatrixSuperLU<double> MatrixLU;
#endif
#endif
#endif
#endif

    MatrixLU lu_mass_operator;
    MatrixLU mat_lu;
    
    TimeLdgScheme1D();

    int GetNbEvaluationDifferential() const
    {
      return nb_eval_diff;
    }

    void InitTimeScheme(int type_scheme, int, const Real_wp& dt, const Real_wp& Tf);
    
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                                    const VectReal_wp& Xn, VectReal_wp& ProdXn,
                                    bool invert_mass = true, bool source = true);
    
    void EvaluateFunction(const Real_wp& tn, const VectReal_wp& Xn, VectReal_wp& ProdXn,
                          bool invert_mass = true, bool source = true);

    void CheckDifferential(VectReal_wp& Un_sol, VectReal_wp& Un_half);    
    void Advance(const Real_wp& t, int n);
    Real_wp GetEnergy() const;

    const VectReal_wp& GetIterate() const;
    const VectReal_wp& GetIterateReal() const;
    
    void GiveIterate(int n, const Real_wp& t, VectReal_wp& Y);
    
    void GetInterpolateUn(int m, VectReal_wp& ZnInterp) const;
    
  };
  
  
  template<class TypeEquation>
  class IdentityOperatorDg1D
  {    
  public :
    IdentityOperatorDg1D() {}
    
    IdentityOperatorDg1D(const TypeEquation& var) {}
    
    template<class Vector1, class Vector2>
    void EvaluateNumericalFlux(int i, int k, const Vector1& Um,
			       const Vector1& Up, Vector2& flux1) const
    {
      flux1.Zero();
    }
    
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
					 const Vector1& Up, Vector2& flux1,
					 Matrix1& df1_dum, Matrix1& df1_dup) const
    {
      flux1.Zero();
      df1_dum.Zero(); df1_dup.Zero();
    }
    

    template<class Vector1, class Vector2>
    void EvaluateNumericalFlux(int i, int k, const Vector1& Um,
			       const Vector1& Up, Vector2& flux1, Vector2&flux2) const
    {
      flux1.Zero(); flux2.Zero();
    }
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
					 const Vector1& Up, Vector2& flux1, Vector2& flux2,
					 Matrix1& df1_dum, Matrix1& df1_dup,
                                         Matrix1& df2_dum, Matrix1& df2_dup) const
    {
      flux1.Zero(); flux2.Zero();
      df1_dum.Zero(); df1_dup.Zero(); df2_dum.Zero(); df2_dup.Zero();
    }
    
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                          Vector2& f, Vector2& g, Vector2& fx, Vector2& gx) const
    {
      fx.Zero(); gx.Zero();
      g.SetIdentity(); f.Zero();
      for (int i = TypeEquation::nb_time_unknowns; i < TypeEquation::nb_unknowns; i++)
        g(i, i) = 0;
    }
    
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
			    Vector2& f, Vector2& g, Vector2& fx, Vector2& gx,
                            Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                            Matrix1& dfx_du, Matrix1& dgx_du,
                            Matrix1& dfx_ddu, Matrix1& dgx_ddu) const
    {
      fx.Zero(); gx.Zero();
      g.SetIdentity(); f.Zero();
      df_du.Zero(); dg_du.Zero(); df_ddu.Zero(); dg_ddu.Zero();
      for (int i = TypeEquation::nb_time_unknowns; i < TypeEquation::nb_unknowns; i++)
        g(i, i) = 0;
    }

  };
  
  
  template<class TypeEquation>
  class IdentityOmegaOperatorDg1D
  {    
  public :
    Complex_wp iomega;
    
    IdentityOmegaOperatorDg1D() {}
    
    IdentityOmegaOperatorDg1D(const TypeEquation& var) {}
    
    template<class Vector1, class Vector2>
    void EvaluateNumericalFlux(int i, int k, const Vector1& Um,
			       const Vector1& Up, Vector2& flux1) const
    {
      flux1.Zero();
    }
    
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
					 const Vector1& Up, Vector2& flux1,
					 Matrix1& df1_dum, Matrix1& df1_dup) const
    {
      flux1.Zero();
      df1_dum.Zero(); df1_dup.Zero();
    }
    

    template<class Vector1, class Vector2>
    void EvaluateNumericalFlux(int i, int k, const Vector1& Um,
			       const Vector1& Up, Vector2& flux1, Vector2&flux2) const
    {
      flux1.Zero(); flux2.Zero();
    }
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
					 const Vector1& Up, Vector2& flux1, Vector2& flux2,
					 Matrix1& df1_dum, Matrix1& df1_dup,
                                         Matrix1& df2_dum, Matrix1& df2_dup) const
    {
      flux1.Zero(); flux2.Zero();
      df1_dum.Zero(); df1_dup.Zero(); df2_dum.Zero(); df2_dup.Zero();
    }
    
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                          Vector2& f, Vector2& g, Vector2& fx, Vector2& gx) const
    {
      fx.Zero(); gx.Zero();
      g.SetIdentity();
      g *= iomega;
      f.Zero();
      for (int i = TypeEquation::nb_time_unknowns; i < TypeEquation::nb_unknowns; i++)
        g(i, i) = 0;
    }
    
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
			    Vector2& f, Vector2& g, Vector2& fx, Vector2& gx,
                            Matrix1& df_du, Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu,
                            Matrix1& dfx_du, Matrix1& dgx_du,
                            Matrix1& dfx_ddu, Matrix1& dgx_ddu) const
    {
      fx.Zero(); gx.Zero();
      g.SetIdentity();
      g *= iomega;
      f.Zero();
      df_du.Zero(); dg_du.Zero();
      df_ddu.Zero(); dg_ddu.Zero();
      for (int i = TypeEquation::nb_time_unknowns; i < TypeEquation::nb_unknowns; i++)
        g(i, i) = 0;
    }

  };
  
}

#define MONTJOIE_FILE_LOCAL_DISCONTINUOUS_GALERKIN_1D_HXX
#endif
