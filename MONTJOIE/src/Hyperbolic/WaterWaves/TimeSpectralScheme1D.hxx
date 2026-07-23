#ifndef MONTJOIE_FILE_TIME_SPECTRAL_SCHEME1D_HXX

namespace Montjoie
{
  
  template<class TypeEquation>
  class TimeSpectralScheme1D : public VirtualOdeSystem<Real_wp>, public VirtualOdeSystem<Complex_wp>,
			       public InputDataProblem_Base
  {
  public :
    //! object for the definition of the equation
    TypeEquation var_eq;
    
  protected :
    //! talezer iterator
    Talezer_Iterator<Complex_wp> talezer_scheme;
    Talezer_Iterator<Real_wp> talezer_scheme_real;
    
    //! Runge-Kutta iterator
    RungeKutta_Iterator<Complex_wp> RK_scheme;
    RungeKutta_Iterator<Real_wp> RK_scheme_real;
    
    //! Low-storage Runge-Kutta iterator
    LowStorageRK_Iterator<Complex_wp> LowRK_scheme;
    LowStorageRK_Iterator<Real_wp> LowRK_scheme_real;
    
    //! Adams-Bashforth iterator
    AdamsBashforth_Moulton_Iterator<Complex_wp> AB_scheme;
    AdamsBashforth_Moulton_Iterator<Real_wp> AB_scheme_real;
    
    //! implicit Runge-Kutta (in real space only)
    GaussRungeKutta_Iterator<Real_wp> GaussRK_scheme;
    
    //! time-scheme used
    int type_scheme, order_time_scheme;    
    
    //! time step
    Real_wp deltat;
    
    //! iterates of time scheme in spectral domain
    VectComplex_wp Un_time, Un_prev, Un_next;
    //! iterates of time scheme in real domain
    VectReal_wp Ur_time, Ur_prev, Ur_next;
    VectReal_wp rhs, Ur_half, AhUn;
    //! object used to perform fft transforms
    FftRealInterface fft;
    //! number of points used for the discretization
    int nb_points_fft;
    //! computational interval 
    Real_wp xmin, xmax;
    //! variable nu used to compute derivatives of the solution
    //! i.e du/dx = F^-1 ( i nu * F(u))
    VectReal_wp nu_fft;
    //! position of points used for the discretization
    VectReal_wp Xdof;
    //! if true the unknown in the real domain is advanced in time, otherwise
    //! it is advanced in the spectral domain
    bool advance_in_real_space;
    //! display grid
    GridInterpolation<Dimension1> grid;
    
    Matrix<Real_wp> diff_matrix; IVect pivot;
    //! finite difference scheme 
    FiniteDifferenceScheme1D<typename TypeEquation::TypeEquationDiff> fd_scheme;
    
  public :
    TimeSpectralScheme1D();

    int GetNbEvaluationDifferential() const
    {
      return 0;
    }
    
    void SetInputData(const string& keyword, const Vector<string>& param);
    
    void SetInterval(const Real_wp&, const Real_wp&, int N, int r);
    void InitGrid(const Real_wp&, const Real_wp&, int N);
    
    void SetOrder(int r);
    
    void InitTimeScheme(int type_scheme, int, const Real_wp& dt, const Real_wp& Tf);

    void FactorizeOperatorDhPlusAijKh(const Matrix<Real_wp>& A, const Real_wp& dt,
                                      const VectReal_wp& Y, const VectReal_wp& ki);
    
    void SolveOperatorDhPlusAijKh(const VectReal_wp& tsub,
                                  const Matrix<Real_wp>& A, const Real_wp& dt,
                                  const VectReal_wp& Y, VectReal_wp& ki);

    void FactorizeOperatorDhPlusAijKh(const Matrix<Real_wp>& A, const Real_wp& dt,
                                      const VectComplex_wp& Y, const VectComplex_wp& ki) { abort(); }

    void SolveOperatorDhPlusAijKh(const VectReal_wp& tsub,
                                  const Matrix<Real_wp>& A, const Real_wp& dt,
                                  const VectComplex_wp& Y, VectComplex_wp& ki) { abort(); }

    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                                    const VectComplex_wp& Xn, VectComplex_wp& ProdXn,
                                    bool invert_mass = true, bool source = true);
    
    void EvaluateFunction(const Real_wp& tn, const VectComplex_wp& Xn, VectComplex_wp& ProdXn,
			  bool invert_mass = true, bool source = true);

    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                                    const VectReal_wp& Xn, VectReal_wp& ProdXn,
                                    bool invert_mass = true, bool source = true);
    
    void EvaluateFunction(const Real_wp& tn, const VectReal_wp& Xn, VectReal_wp& ProdXn,
			  bool invert_mass = true, bool source = true);
    
    const VectComplex_wp& GetIterate() const;
    const VectReal_wp& GetIterateReal() const;
    
    void Advance(const Real_wp& t, int n);
    
    void GetInterpolateUn(int m, VectReal_wp& ZnInterp);
    
    Real_wp GetEnergy()
    {
      if (advance_in_real_space)
	return GetEnergy(this->GetIterateReal());
      
      return 0.0;
    }
    
    Real_wp GetEnergy(const VectReal_wp& Zn)
    {
      return Norm2(Zn);
    }

    void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp&, Real_wp alpha = 1.0)
    {}
    
    void SetDirichletCondition(const Real_wp& t, int n, VectComplex_wp&, Real_wp alpha = 1.0)    
    {}
    
  };

}

#define MONTJOIE_FILE_TIME_SPECTRAL_SCHEME1D_HXX
#endif

