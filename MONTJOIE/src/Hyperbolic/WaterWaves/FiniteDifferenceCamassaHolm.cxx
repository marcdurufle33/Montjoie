#ifndef MONTJOIE_FILE_FINITE_DIFFERENCE_CAMASSA_HOLM_CXX

#include "Algebra/TinyBandMatrix.cxx"

namespace Montjoie
{
  
  template<>
  class FiniteDifferenceScheme1D<CamassaHolmEquation> 
    : public FiniteDifferenceScheme1D_Base<CamassaHolmEquation>
  {
  protected :
    Vector<Real_wp> Zn_mass, Zn_stiff;
    TinyArrowMatrix<Real_wp, 2, 2> linear_part, mat_sparse;
    TinyArrowMatrix<Real_wp, 2, 2> mat_lu;
    TinyArrowMatrix<Real_wp, 1, 1> mass_operator, mass_lu;;
    Vector<Real_wp> coefVar_uux, coefVar_u2ux, coefVar_u3ux, chalf, coef_ones;
    Vector<Real_wp> coefA_uxxx, coefB_uxxx, coefVar_Etilde, coefVar_Ftilde, coefVar_FtildeDx;

    TinyArrowMatrix<Real_wp, 5, 5> mat_lu22, mat_sparse22;
    TinyArrowMatrix<Real_wp, 8, 8> mat_lu33, mat_sparse33;
    TinyArrowMatrix<Real_wp, 11, 11> mat_lu44, mat_sparse44;
    
  public :
    
    //! default constructor
    FiniteDifferenceScheme1D() : FiniteDifferenceScheme1D_Base<CamassaHolmEquation>()
    {
    }
    
    
    //! initialisation of physical coefficients needed to evaluate the scheme
    void SetOrder(int order)
    {
      int N = nb_points;

      // subdivision x_i in space
      VectReal_wp x;
      GetSubdivisionSpace(x);

      // computation of b and c
      VectReal_wp b, c;
      var_eq.ComputeBottom(x0, xN, x, b, c);
      
      //c.WriteText("bottom.dat");
      
      // on calcule c_{i+1/2} = 1/2(c_i + c_i+1)
      chalf.Reallocate(N);
      for (int i = 0; i < N-1; i++)
        chalf(i) = 0.5*(c(i) + c(i+1));
      
      if (type_condition_right == PERIODIC)
        chalf(N-1) = 0.5*(c(N-1) + c(0));
      else
        chalf(N-1) = 1.5*c(N-1) - 0.5*c(N-2);
      
      coef_ones.Reallocate(N); coef_ones.Fill(1.0);
      if (var_eq.type_model == var_eq.GENTLE)
        {
          coefVar_uux.Reallocate(N);
          coefVar_u2ux.Reallocate(N);
          coefVar_u3ux.Reallocate(N);
          coefVar_uux.Fill(1); coefVar_u2ux.Fill(1); coefVar_u3ux.Fill(1);
          
          coefA_uxxx.Reallocate(N);
          coefB_uxxx.Reallocate(N);
          coefA_uxxx.Fill(sqrt(var_eq.mu/6));
          coefB_uxxx.Fill(sqrt(var_eq.mu/12));
        }
      else if (var_eq.type_model == var_eq.ORIGINAL)
        {
          coefVar_uux.Reallocate(N);
          coefVar_u2ux.Reallocate(N);
          coefVar_u3ux.Reallocate(N);
          coefA_uxxx.Reallocate(N);
          coefB_uxxx.Reallocate(N);
          coefVar_Etilde.Reallocate(N);
          coefVar_Ftilde.Reallocate(N);
          coefVar_FtildeDx.Reallocate(N);
          for (int i = 0; i < N; i++)
            {
              coefVar_uux(i) = 1.0/c(i);
              coefVar_u2ux(i) = pow(1.0/c(i), 3.0);
              coefVar_u3ux(i) = pow(1.0/c(i), 5.0);
              
              coefA_uxxx(i) = var_eq.mu*(pow(c(i), 5.0)/6.0 - c(i)/12.0);
              coefB_uxxx(i) = 0.0;
              
              Real_wp dc = 0;
              if (i == 0)
                dc = (chalf(1) - chalf(N-1))/(2.0*dx);
              else if (i == N-1)
                dc = (chalf(0) - chalf(N-2))/(2.0*dx);
              else
                dc = (chalf(i+1) - chalf(i-1))/(2.0*dx);
              
              coefVar_Etilde(i) = -pow(c(i), 3.0)/6.0 - 1.0/(8.0*c(i));
              coefVar_Ftilde(i) = -5.0/24*pow(chalf(i), 3.0) - 9.0/(24.0*chalf(i));
              coefVar_FtildeDx(i) = -5.0/24*3.0*dc*pow(chalf(i), 2.0)
                + 9.0*dc/(24.0*square(chalf(i)));
            }
        }
      else
        {
          coefVar_uux.Reallocate(N);
          coefVar_u2ux.Reallocate(N);
          coefVar_u3ux.Reallocate(N);
          coefA_uxxx.Reallocate(N);
          coefB_uxxx.Reallocate(N);
          for (int i = 0; i < N; i++)
            {
              coefVar_uux(i) = pow(1.0/c(i), 1.0/3.0);
              coefVar_u2ux(i) = pow(1.0/c(i), 0.75);
              coefVar_u3ux(i) = 1.0/c(i);
              
              coefA_uxxx(i) = sqrt(var_eq.mu*pow(c(i), 5.0)/6);
              coefB_uxxx(i) = sqrt(var_eq.mu*c(i)/12);
            }
        }
      
      // on calcule la partie lineaire 
      linear_part.Clear(); linear_part.Reallocate(N, N);
      
      if (var_eq.type_model == var_eq.ORIGINAL)
        {
          AddMatrixD3(1.0, coefA_uxxx, coef_ones, linear_part);
        }
      else
        {
          Vector<Real_wp> coef_mu(N);
          Copy(coefA_uxxx, coef_mu);
          AddMatrixD3(1.0, coef_mu, coef_mu, linear_part);
          
          Copy(coefB_uxxx, coef_mu);
          AddMatrixD3(-1.0, coef_mu, coef_mu, linear_part);
        }
      
      AddMatrixD1var(1.0, chalf, linear_part);

      // operateur de masse
      mass_operator.Clear(); mass_operator.Reallocate(N, N);
      AddMatrixD2(-var_eq.mu/12.0, coef_ones, coef_ones, mass_operator, 0, 0);
      for (int i = 0; i < N; i++)
	mass_operator.AddInteraction(i, i, 1.0);      
      
      GetLU(mass_operator, mass_lu, true);
      
      mat_sparse.Reallocate(N, N);
    }

    
    //! computes ProdXn = g(Xn) where the evolution system is equal to dXn/dt = g(Xn)
    void GetFunction(const VectReal_wp& Xn, VectReal_wp& ProdXn)
    {
      ProdXn.Fill(0);
      AddFunction(1.0, Xn, ProdXn);
    }
    

    //! computes Yn = Yn + alpha g(Xn)
    //! where the evolution system is equal to dXn/dt = g(Xn)
    void AddFunction(const Real_wp& alpha, const VectReal_wp& Xn,
                     VectReal_wp& Yn, bool invert_mass = true)
    {
      VectReal_wp ProdXn(nb_points); ProdXn.Fill(0);
      
      // linear part (c zeta_x + c_x/2 zeta + mu/12 zeta_xxx) for gentle model
      MltAdd(-alpha, linear_part, Xn, Real_wp(1), ProdXn);
      
      if (var_eq.type_model == var_eq.ORIGINAL)
        {
          // part 3/(2c) epsilon zeta zeta_x
          AddOperatorUpUx(-1.5*alpha*var_eq.epsilon, 1,
                          coefVar_uux, coef_ones, Xn, Xn, ProdXn);
          
          // part - 3/(8 c^3) epsilon^2 zeta^2 zeta_x
          AddOperatorUpUx(0.375*alpha*square(var_eq.epsilon), 2,
                          coefVar_u2ux, coef_ones, Xn, Xn, ProdXn);
          
          // part + 3/(16 c^5) epsilon^3 zeta^3 zeta_x
          AddOperatorUpUx(-3.0/16*alpha*pow(var_eq.epsilon, 3.0), 3,
                          coefVar_u3ux, coef_ones, Xn, Xn, ProdXn);
          
          // part - Etilde epsilon mu zeta zeta_xxx
          AddOperatorUxxxU(alpha*var_eq.epsilon*var_eq.mu,
                           coefVar_Etilde, coef_ones, Xn, Xn, ProdXn);
          
          // part - Ftilde epsilon mu zeta_xx (F zeta_x + F_x/2 zeta)
          AddOperatorUxxUx(alpha*var_eq.epsilon*var_eq.mu,
                           coef_ones, coef_ones, coefVar_Ftilde, Xn, Xn, ProdXn);

          // part - Ftilde epsilon mu zeta_x (F_x zeta_x + F_xx/2 zeta)
          AddOperatorUxUx(alpha*var_eq.epsilon*var_eq.mu,
                          coef_ones, coef_ones, coefVar_FtildeDx, Xn, Xn, ProdXn);
        }
      else
        {
          // part 3/2 epsilon zeta zeta_x
          AddOperatorUpUx(-1.5*alpha*var_eq.epsilon, 1,
                          coefVar_uux, coefVar_uux, Xn, Xn, ProdXn);
          
          // part - 3/8 epsilon^2 zeta^2 zeta_x
          AddOperatorUpUx(0.375*alpha*square(var_eq.epsilon), 2,
                          coefVar_u2ux, coefVar_u2ux, Xn, Xn, ProdXn);
          
          // part + 3/16 epsilon^3 zeta^3 zeta_x
          AddOperatorUpUx(-3.0/16*alpha*pow(var_eq.epsilon, 3.0), 3,
                          coefVar_u3ux, coefVar_u3ux, Xn, Xn, ProdXn);
          
          // part 7/24 epsilon mu (zeta zeta_xxx + 2 zeta_x zeta_xx)
          AddOperatorUxxUx_UxxxU(-7.0*var_eq.epsilon*var_eq.mu/24*alpha,
                                 coef_ones, coef_ones, Xn, Xn, ProdXn);
        }
      
      if (invert_mass)
        {
          // resolution by mass matrix
          SolveLU(mass_lu, ProdXn);
        }
      
      Add(1.0, ProdXn, Yn);
    }
    
    
    //! computation of A = A + alpha * M
    template<class MatrixSparse>
    void AddMass(MatrixSparse& A,
                 int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      int N = nb_points;
      for (int i = 0; i < N; i++)
        for (int j2 = i-1; j2 <= i+1; j2++)
          {
            int j = j2;
            if (j2 < 0)
              j = j2 + N;
            
            if (j2 >= N)
              j = j2 - N;
            
            A.AddInteraction(m + incx*i, n + incy*j, mass_operator(i, j));
          }
    }
    
    
    //! computation of A = A + alpha * M
    void AddMass(int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      switch (incx)
        {
        case 1 :
          AddMass(mat_sparse, m, n, incx, incy);
          break;
        case 2 :
          AddMass(mat_sparse22, m, n, incx, incy);
          break;
        case 3 :
          AddMass(mat_sparse33, m, n, incx, incy);
          break;
        case 4 :
          AddMass(mat_sparse44, m, n, incx, incy);
          break;
        default :
          {
            cout << "case not handled" << endl;
            abort();
          }
        }
    }

    //! computation of A = A + alpha * DF
    template<class MatrixSparse>
    void AddStiffJacobian(const Real_wp& alpha,
                          const VectReal_wp& Un_demi,
                          MatrixSparse& A,
                          int m = 0, int n = 0, int incx = 1, int incy = 1)
    {    
      int N = nb_points;
      for (int i = 0; i < N; i++)
        for (int j2 = i-2; j2 <= i+2; j2++)
          {
            int j = j2;
            if (j2 < 0)
              j = j2 + N;
            
            if (j2 >= N)
              j = j2 - N;
            
            A.AddInteraction(m + incx*i, n + incy*j, -alpha*linear_part(i, j));
          }
      
      if (var_eq.type_model == var_eq.ORIGINAL)
        {
          // part 3/(2c) epsilon zeta zeta_x
          AddDifferentialUpUx(-1.5*alpha*var_eq.epsilon, 1,
                              coefVar_uux, coef_ones, Un_demi, A, m, n, incx, incy);
          
          // part - 3/(8 c^3) epsilon^2 zeta^2 zeta_x
          AddDifferentialUpUx(0.375*alpha*square(var_eq.epsilon), 2,
                              coefVar_u2ux, coef_ones, Un_demi, A, m, n, incx, incy);
          
          // part + 3/(16 c^5) epsilon^3 zeta^3 zeta_x
          AddDifferentialUpUx(-3.0/16*alpha*pow(var_eq.epsilon, 3.0), 3,
                              coefVar_u3ux, coef_ones, Un_demi, A, m, n, incx, incy);
          
          // part - Etilde epsilon mu zeta zeta_xxx
          AddDifferentialUxxxU(alpha*var_eq.epsilon*var_eq.mu,
                               coefVar_Etilde, coef_ones, Un_demi, A, m, n, incx, incy);
          
          // part - Ftilde epsilon mu zeta_xx (F zeta_x + F_x/2 zeta)
          AddDifferentialUxxUx(alpha*var_eq.epsilon*var_eq.mu,
                               coef_ones, coef_ones, coefVar_Ftilde, Un_demi, A, m, n, incx, incy);

          // part - Ftilde epsilon mu zeta_x (F_x zeta_x + F_xx/2 zeta)
          AddDifferentialUxUx(alpha*var_eq.epsilon*var_eq.mu,
                              coef_ones, coef_ones, coefVar_FtildeDx, Un_demi,
                              A, m, n, incx, incy);
          
        }
      else
        {
          // part 3/2 epsilon zeta zeta_x
          AddDifferentialUpUx(-1.5*alpha*var_eq.epsilon, 1,
                              coefVar_uux, coefVar_uux, Un_demi, A, m, n, incx, incy);
          
          // part - 3/8 epsilon^2 zeta^2 zeta_x
          AddDifferentialUpUx(0.375*alpha*square(var_eq.epsilon), 2,
                              coefVar_u2ux, coefVar_u2ux, Un_demi, A, m, n, incx, incy);
          
          // part + 3/16 epsilon^3 zeta^3 zeta_x
          AddDifferentialUpUx(-3.0/16*alpha*pow(var_eq.epsilon, 3.0), 3,
                              coefVar_u3ux, coefVar_u3ux, Un_demi, A, m, n, incx, incy);
          
          // part 7/24 epsilon mu (zeta zeta_xxx + 2 zeta_x zeta_xx)
          AddDifferentialUxxUx_UxxxU(-7.0*var_eq.epsilon*var_eq.mu/24*alpha,
                                     coef_ones, coef_ones, Un_demi, A, m, n, incx, incy);
        }
    }
    
    
    //! computation of A = A + alpha * DF
    void AddStiffJacobian(const Real_wp& alpha,
                          const VectReal_wp& Un_demi,
                          int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      switch (incx)
        {
        case 1 :
          AddStiffJacobian(alpha, Un_demi, mat_sparse, m, n, incx, incy);
          break; 
        case 2 :
          AddStiffJacobian(alpha, Un_demi, mat_sparse22, m, n, incx, incy);
          break;
        case 3 :
          AddStiffJacobian(alpha, Un_demi, mat_sparse33, m, n, incx, incy);
          break;
        case 4 :
          AddStiffJacobian(alpha, Un_demi, mat_sparse44, m, n, incx, incy);
          break;
        default :
          {
            cout << "case not handled" << endl;
            abort();
          }
        }
    }
    
    
    //! allocation of jacobian matrix
    void InitJacobian(int incx)
    {
      int N = nb_points;
      switch (incx)
        {
        case 1 :
          mat_sparse.Reallocate(N, N);
          mat_sparse.Zero();
          break;
        case 2 :
          mat_sparse22.Reallocate(2*N, 2*N);
          mat_sparse22.Zero();
          break;
        case 3 :
          mat_sparse33.Reallocate(3*N, 3*N);
          mat_sparse33.Zero();
          break;
        case 4 :
          mat_sparse44.Reallocate(4*N, 4*N);
          mat_sparse44.Zero();
          break;
        default :
          {
            cout << "case forbidden" << endl;
            abort();
          }
        }
    }
    
    
    //! factorisation of I + alpha DF
    void FactoriseJacobian(int incx)
    {
      switch (incx)
        {
        case 1 :
          mat_lu.HideMessages();
          GetLU(mat_sparse, mat_lu, true);
          break;
        case 2 :
          mat_lu22.HideMessages();
          GetLU(mat_sparse22, mat_lu22, true);
          break;
        case 3 :
          mat_lu33.HideMessages();
          GetLU(mat_sparse33, mat_lu33, true);
          break;
        case 4 :
          mat_lu44.HideMessages();
          GetLU(mat_sparse44, mat_lu44, true);
          break;
        default :
          {
            cout << "Forbidden case " << endl;
            abort();
          }
        }
    }
    
    
    //! Computes matrix M + alpha DF and factorises this matrix
    void ComputeAndFactoriseJacobian(const Real_wp& alpha,
                                     const VectReal_wp& Un_demi)
    {
      mat_sparse = mass_operator;

      AddStiffJacobian(alpha, Un_demi, mat_sparse);
      
      mat_lu.HideMessages();
      GetLU(mat_sparse, mat_lu, true);
    }
    
    
    //! solve the system (I + alpha DF) x = rhs
    /*!
      \param[inout] rhs on input the right hand side, on output the solution
     */
    void SolveGaussJacobian(VectReal_wp& rhs, int incx)
    {
      switch (incx)
        {
        case 1 :
          SolveLU(mat_lu, rhs);
          break;
        case 2 :
          SolveLU(mat_lu22, rhs);
          break;
        case 3 :
          SolveLU(mat_lu33, rhs);
          break;
        case 4 :
          SolveLU(mat_lu44, rhs);
          break;
        default :
          {
            cout << "Not implemented" << endl;
            abort();
          }
        }
    }

    
    //! solve the system (M + alpha DF) x = rhs
    /*!
      \param[inout] rhs on input the right hand side, on output the solution
     */
    void SolveJacobian(VectReal_wp& rhs)
    {
      SolveLU(mat_lu, rhs);
    }
    
    
    //! application of mass matrix to a vector
    void ApplyMass(const VectReal_wp& U, VectReal_wp& V)
    {      
      MltAdd(1.0, mass_operator, U, 0.0, V);
    }

    
    //! advancing semi-implicit scheme 
    void AdvancePredictorScheme(Real_wp t, int iter_num)
    {
      glob_chrono.Start(VirtualTimer::MASS);
      
      // we consider linear part
      mat_sparse.Zero();
      Add(0.5*dt, linear_part, mat_sparse);
      
      // and we add non-linear part        
      if (var_eq.type_model == var_eq.ORIGINAL)
        {
          // part 3/(2c) epsilon zeta zeta_x
          AddDifferentialHalfUpUx(0.5*dt*1.5*var_eq.epsilon, 1,
                                  coefVar_uux, coef_ones, Un_demi, mat_sparse);

          // part - 3/(8 c^3) epsilon^2 zeta^2 zeta_x
          AddDifferentialHalfUpUx(-0.5*dt*0.375*square(var_eq.epsilon), 2,
                                  coefVar_u2ux, coef_ones, Un_demi, mat_sparse);
          
          // part + 3/(16 c^5) epsilon^3 zeta^3 zeta_x
          AddDifferentialHalfUpUx(0.5*dt*3.0/16*pow(var_eq.epsilon, 3.0), 3,
                                  coefVar_u3ux, coef_ones, Un_demi, mat_sparse);
          
          // part - Etilde epsilon mu zeta zeta_xxx
          AddDifferentialHalfUxxxU(-0.5*dt*var_eq.epsilon*var_eq.mu,
                                   coefVar_Etilde, coef_ones, Un_demi, mat_sparse);
          
          // part - Ftilde epsilon mu zeta_xx (F zeta_x + F_x/2 zeta)
          AddDifferentialHalfUxxUx(-0.5*dt*var_eq.epsilon*var_eq.mu,
                                   coef_ones, coef_ones, coefVar_Ftilde, Un_demi, mat_sparse);

          // part - Ftilde epsilon mu zeta_x (F_x zeta_x + F_xx/2 zeta)
          AddDifferentialHalfUxUx(-0.5*dt*var_eq.epsilon*var_eq.mu,
                                  coef_ones, coef_ones, coefVar_FtildeDx, Un_demi, mat_sparse);
        }
      else
        {
          // part 3/2 epsilon zeta zeta_x
          AddDifferentialHalfUpUx(0.5*dt*1.5*var_eq.epsilon, 1,
                                  coefVar_uux, coefVar_uux, Un_demi, mat_sparse);

          // part - 3/8 epsilon^2 zeta^2 zeta_x
          AddDifferentialHalfUpUx(-0.5*dt*0.375*square(var_eq.epsilon), 2,
                                  coefVar_u2ux, coefVar_u2ux, Un_demi, mat_sparse);
          
          // part + 3/16 epsilon^3 zeta^3 zeta_x
          AddDifferentialHalfUpUx(0.5*dt*3.0/16*pow(var_eq.epsilon, 3.0), 3,
                                  coefVar_u3ux, coefVar_u3ux, Un_demi, mat_sparse);
          
          // part 7/24 epsilon mu (zeta zeta_xxx + 2 zeta_x zeta_xx)
          AddDifferentialHalfUxxUx_UxxxU(0.5*dt*7.0*var_eq.epsilon*var_eq.mu/24,
                                         coef_ones, coef_ones, Un_demi, mat_sparse);
        }
      
      // on calcule (M-K) U^n
      ApplyMass(Un, AhUn);
      
      MltAdd(-1.0, mat_sparse, Un, 1.0, AhUn);
      
      // on factorise (M+K)
      Add(1.0, mass_operator, mat_sparse);
      
      glob_chrono.Stop(VirtualTimer::MASS);
      glob_chrono.Start(VirtualTimer::STIFFNESS);
      
      mat_lu.HideMessages();
      
      GetLU(mat_sparse, mat_lu, true);
      
      SolveLU(mat_lu, AhUn);  
      Copy(AhUn, Un);
	  
      // on calcule U_n+3/2 = 2 Un+1 - Un+1/2
      Mlt(-1.0, Un_demi); Add(2.0, Un, Un_demi);
      
      glob_chrono.Stop(VirtualTimer::STIFFNESS);
    }

    
  };
    
}

#define MONTJOIE_FILE_FINITE_DIFFERENCE_CAMASSA_HOLM_CXX
#endif
