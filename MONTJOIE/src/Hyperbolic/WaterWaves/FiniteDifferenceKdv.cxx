#ifndef MONTJOIE_FILE_FINITE_DIFFERENCE_KDV_CXX

#include "Algebra/TinyBandMatrix.cxx"

namespace Montjoie
{
  
  //! Class for the resolution of Kdv equation with finite-difference scheme
  template<>
  class FiniteDifferenceScheme1D<KdvEquation>
    : public FiniteDifferenceScheme1D_Base<KdvEquation>
  {
  public :
    //! sparse matrix storing the linear part of the operator
    TinyArrowMatrix<Real_wp, 2, 2> linear_part, mat_sparse;   
    //! factorisation of a pentadiagonal matrix with periodic conditions
    TinyArrowMatrix<Real_wp, 2, 2> mat_lu;
    
    TinyArrowMatrix<Real_wp, 5, 5> mat_lu22, mat_sparse22;
    TinyArrowMatrix<Real_wp, 8, 8> mat_lu33, mat_sparse33;
    TinyArrowMatrix<Real_wp, 11, 11> mat_lu44, mat_sparse44;

    //! coefD = (1/c)^(1/3) and invVelocity = 1/c
    Vector<Real_wp> coefD, invVelocity;
    
    //! default constructor
    FiniteDifferenceScheme1D() : FiniteDifferenceScheme1D_Base<KdvEquation>()
    {
    }
    
    
    //! computation of physical coefficients and other intermediary matrices
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
         
      // we compute c_{i+1/2} = 1/2(c_i + c_i+1)
      Vector<Real_wp> chalf(N);
      for (int i = 0; i < N-1; i++)
        chalf(i) = 0.5*(c(i) + c(i+1));
      
      if (type_condition_right == PERIODIC)
        chalf(N-1) = 0.5*(c(N-1) + c(0));
      else
        chalf(N-1) = 1.5*c(N-1) - 0.5*c(N-2);
      
      if (var_eq.type_model == var_eq.GENTLE)
	{	  
	  // we compute linear part -> D1^v + mu/6 D^3
	  linear_part.Clear(); linear_part.Reallocate(N, N);
	  Vector<Real_wp> coef_mu(N);
          coef_mu.Fill(var_eq.mu/6.0);
	  AddMatrixD3var(1.0, coef_mu, linear_part);
	  
	  AddMatrixD1var(1.0, chalf, linear_part);
	  coefD.Reallocate(N); coefD.Fill(1.0);
	}
      else
	{
          // we compute linear part of the stiffness operator
          linear_part.Clear(); linear_part.Reallocate(N, N);
	  Vector<Real_wp> coef_mu(N), ones(N);
          ones.Fill(1.0);
          coef_mu.Fill(var_eq.mu/6.0);
	  for (int i = 0; i < N; i++)
	    coef_mu(i) *= pow(c(i), 5.0);
	  
	  if (var_eq.type_model == var_eq.ORIGINAL)
            AddMatrixD3(1.0, coef_mu, ones, linear_part);
          else
            AddMatrixD3var(1.0, coef_mu, linear_part);
	  
	  AddMatrixD1var(1.0, chalf, linear_part);
	  
          // coefD = (1/c)^1/3 for violent model and 1 for original model
          coefD.Reallocate(N);
          if (var_eq.type_model == var_eq.ORIGINAL)
            {
              coefD.Fill(1.0);
              invVelocity.Reallocate(N);
              for (int i = 0; i < N; i++)
                invVelocity(i) = 1.0/c(i);
            }
          else
            for (int i = 0; i < N; i++)
              coefD(i) = pow(1.0/c(i), 1.0/3.0);	  
	}
      
      // matrix used for implicit schemes
      mat_sparse.Reallocate(N, N);	   
    }
    
    
    //! computes ProdXn = g(Xn) where the evolution system is equal to dXn/dt = g(Xn)
    void GetFunction(const VectReal_wp& Xn, VectReal_wp& ProdXn)
    {
      ProdXn.Fill(0);
      AddFunction(1.0, Xn, ProdXn);
    }
    

    //! computes ProdXn = ProdXn + alpha g(Xn)
    //! where the evolution system is equal to dXn/dt = g(Xn)
    void AddFunction(const Real_wp& alpha, const VectReal_wp& Xn,
                     VectReal_wp& ProdXn, bool invert_mass = true)
    {
      // we add the part -alpha ( c zeta_x + c_x/2 zeta + mu/6 Gamma_3 zeta)
      // where Gamma_3 zeta = c^5 zeta_xxx for original model
      MltAdd(-alpha, linear_part, Xn, Real_wp(1), ProdXn);
      
      if (var_eq.type_model == var_eq.ORIGINAL)
	{
          AddOperatorUpUx(-1.5*alpha*var_eq.epsilon, 1,
                          invVelocity, coefD, Xn, Xn, ProdXn);
        }
      else
        {
          // we add the part -alpha ( 3/2 epsilon (1/c)^1/3  (1/c)^1/3 zeta  ( (1/c)^1/3 zeta )_x
          // with the scheme (u_{i+1} - u_{i-1})/(2 dx) * (u_{i-1} + u_i + u_{i+1})/3
          AddOperatorUpUx(-1.5*alpha*var_eq.epsilon, 1,
                          coefD, coefD, Xn, Xn, ProdXn);
        }
    }
    
    
    //! computation of A = A + alpha * M
    template<class MatrixSparse>
    void AddMass(MatrixSparse& A,
                 int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      for (int i = 0; i < nb_points; i++)
        A.AddInteraction(m + incx*i, n + incy*i, 1.0);      
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
          AddDifferentialUpUx(-1.5*alpha*var_eq.epsilon, 1,
                              invVelocity, coefD, Un_demi, A,
                              m, n, incx, incy);
        }
      else
        AddDifferentialUpUx(-1.5*alpha*var_eq.epsilon, 1,
                            coefD, coefD, Un_demi, A,
                            m, n, incx, incy);
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
    
    
    //! Computes matrix I + alpha DF and factorises this matrix
    void ComputeAndFactoriseJacobian(const Real_wp& alpha,
                                     const VectReal_wp& Un_demi)
    {
      mat_sparse.SetIdentity();
      
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

    
    //! solve the system (I + alpha DF) x = rhs
    /*!
      \param[inout] rhs on input the right hand side, on output the solution
     */
    void SolveJacobian(VectReal_wp& rhs)
    {
      SolveLU(mat_lu, rhs);
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
        AddDifferentialHalfUpUx(0.75*var_eq.epsilon*dt, 1,
                                invVelocity, coefD, Un_demi, mat_sparse);
      else
        AddDifferentialHalfUpUx(0.75*var_eq.epsilon*dt, 1, coefD, coefD, Un_demi, mat_sparse);
      
      // on calcule (I-M) U^n
      Copy(Un, AhUn);
      
      MltAdd(-1.0, mat_sparse, Un, 1.0, AhUn);
      
      // on factorise (I+M)
      for (int i = 0; i < nb_points; i++)
        mat_sparse.AddInteraction(i, i, 1.0);
      
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

#define MONTJOIE_FILE_FINITE_DIFFERENCE_KDV_CXX
#endif
