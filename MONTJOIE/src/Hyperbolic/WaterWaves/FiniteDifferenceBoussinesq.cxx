#ifndef MONTJOIE_FILE_FINITE_DIFFERENCE_BOUSSINESQ_CXX

namespace Montjoie
{
  
  /* template<>
  class FiniteDifferenceScheme1D<BoussinesqEquation>
    : public FiniteDifferenceScheme1D_Base<BoussinesqEquation>
  {
  public :
    
    Vector<Real_wp> Xn, Xn_mass, Xn_stiff, Xn_next, rhs;
    
    Vector<Real_wp> C2, C4Half, Ones;
    TinyArrowMatrix<Real_wp, 3, 2> mat_sparse;
    TinyArrowMatrix<Real_wp, 3, 2> mat_lu;
    Real_wp dt, alpha, dx;
    bool compute_differential;
    int nb_eval_diff;
    
    FiniteDifferenceScheme1D() : FiniteDifferenceScheme1D_Base<BoussinesqEquation>()
    {
      dt = 0;
      alpha = 2.0/3;
      dx = 0;
      compute_differential = true;
      nb_eval_diff = 0;
    }
    
    VectReal_wp GetLastSolution(int n) const
    {
      VectReal_wp Zn(nb_points);
      for (int i = 0; i < nb_points; i++)
        Zn(i) = Xn(2*i+n);
      
      return Zn;
    }
    
    void InitTimeScheme(Real_wp dt_, Real_wp Tfinal)
    {
      dt = dt_;
      int N = nb_points;
      Xn.Reallocate(2*N); Xn.Fill(0);
      Xn_next.Reallocate(2*N); Xn_next.Fill(0);
      Xn_mass.Reallocate(2*N); Xn_mass.Fill(0);
      Xn_stiff.Reallocate(2*N); Xn_stiff.Fill(0);
      rhs.Reallocate(2*N); rhs.Fill(0);
      
      dx = (xN - x0)/(nb_points-1);
      VectReal_wp x(nb_points);
      for (int i = 0; i < nb_points; i++)
	x(i) = x0 + i*dx;
      
      VectReal_wp b, c;
      var_eq.ComputeBottom(x0, xN, x, b, c);
      
      // on calcule c_{i+1/2}
      Vector<Real_wp> xhalf(N), bhalf(N), chalf(N);
      for (int i = 0; i < N-1; i++)
	xhalf(i) = x(i) + 0.5*dx; 
      
      var_eq.ComputeBottom(x0, xN, xhalf, bhalf, chalf);
      
      C2.Reallocate(N);
      for (int i = 0; i < N; i++)
	C2(i) = c(i)*c(i);
      
      C4Half.Reallocate(N);
      for (int i = 0; i < N; i++)
	C4Half(i) = pow(chalf(i), 4.0);
      
      Ones.Reallocate(N); Ones.Fill(1.0);
      
      // initial condition
      VectReal_wp Y0(2*N);
      var_eq.ComputeInitialCondition(x0, xN, x, Y0);
      for (int i = 0; i < N; i++)
	{
	  Xn(2*i) = Y0(i);
	  Xn(2*i+1) = Y0(i+N);
	}
      
      mat_sparse.Reallocate(2*N, 2*N);
      compute_differential = true;
      nb_eval_diff = 0;
    }
    
    
    void ApplyMassOperator(Real_wp alpha, const Vector<Real_wp>& v, Vector<Real_wp>& y)
    {
      // part with identity
      for (int i = 0; i < 2*nb_points; i++)
        y(i) += alpha*v(i);
      
      // part - mu/3 (C4 u_x)_x
      Real_wp beta = alpha*var_eq.mu/(3.0*dx*dx);
      for (int i = 1; i < nb_points-1; i++)
        y(2*i+1) += beta*(-C4Half(i)*v(2*i+3) + (C4Half(i) + C4Half(i-1))*v(2*i+1)
                          - C4Half(i-1)*v(2*i-1));
      
      y(1) += beta*(-C4Half(0)*v(3) + (C4Half(0) + C4Half(nb_points-1))*v(1)
                    - C4Half(nb_points-1)*v(2*nb_points-1));
      
      y(2*nb_points-1) += beta*(-C4Half(nb_points-1)*v(1)
                                + (C4Half(nb_points-1) + C4Half(nb_points-2))*v(2*nb_points-1)
                                - C4Half(nb_points-2)*v(2*nb_points-3));
    }
    
    
    template<class MatrixSparse>
    void AddMassDifferential(Real_wp alpha, const Vector<Real_wp>& v, MatrixSparse& A)
    {
      // part with identity
      for (int i = 0; i < 2*nb_points; i++)
        A.AddInteraction(i, i, alpha);
      
      // part - mu/3 (C4 u_x)_x
      Real_wp beta = alpha*var_eq.mu/(3.0*dx*dx);
      for (int i = 1; i < nb_points-1; i++)
        {
          int irow = 2*i+1;
          A.AddInteraction(irow, irow, beta*(C4Half(i) + C4Half(i-1)));
          A.AddInteraction(irow, irow+2, -beta*C4Half(i));
          A.AddInteraction(irow, irow-2, -beta*C4Half(i-1));
        }
      
      int N = 2*nb_points-1;
      A.AddInteraction(1, 1, beta*(C4Half(0) + C4Half(nb_points-1)));
      A.AddInteraction(1, 3, -beta*C4Half(0));
      A.AddInteraction(1, N, -beta*C4Half(nb_points-1));

      A.AddInteraction(N, N, beta*(C4Half(nb_points-1) + C4Half(nb_points-2)));
      A.AddInteraction(N, 1, -beta*C4Half(nb_points-1));
      A.AddInteraction(N, N-2, -beta*C4Half(nb_points-2));      
    }
    
    
    void ApplyStiffnessOperator(Real_wp alpha, const Vector<Real_wp>& v, Vector<Real_wp>& y)
    {
      // part [h u]_x
      Real_wp beta = alpha/(2.0*dx);
      for (int i = 1; i < nb_points-1; i++)
        {
          int j = 2*i;
          y(j) += beta*( (C2(i+1) + var_eq.epsilon*v(j+2))*v(j+3)
                         -(C2(i-1) + var_eq.epsilon*v(j-2))*v(j-1) );
        }
      
      int N = 2*nb_points-1;
      y(0) += beta*( (C2(1) + var_eq.epsilon*v(2))*v(3)
                     -(C2(nb_points-1) + var_eq.epsilon*v(N-1))*v(N) );

      y(N-1) += beta*( (C2(0) + var_eq.epsilon*v(0))*v(1)
                       -(C2(nb_points-2) + var_eq.epsilon*v(N-3))*v(N-2) );
      
      // part zeta_x + epsilon u u_x
      Real_wp coef = var_eq.epsilon/3;
      for (int i = 1; i < nb_points-1; i++)
        {
          int j = 2*i;
          y(j+1) += beta*( (v(j+2) - v(j-2)) + coef*(v(j+3) - v(j-1))
                           *(v(j-1) + v(j+1) + v(j+3)) ); 
        }
      
      y(1) += beta*( (v(2) - v(N-1)) + coef*(v(3) - v(N))*(v(N) + v(1) + v(3)) ); 
      y(N) += beta*( (v(0) - v(N-3)) + coef*(v(1) - v(N-2))*(v(N-2) + v(N) + v(1)) );
      
    }
    

    template<class MatrixSparse>
    void AddStiffnessDifferential(Real_wp alpha, const Vector<Real_wp>& v, MatrixSparse& A)
    {
      // part [h u]_x
      Real_wp beta = alpha/(2.0*dx);
      for (int i = 1; i < nb_points-1; i++)
        {
          int j = 2*i;
          A.AddInteraction(j, j+2, beta*var_eq.epsilon*v(j+3));
          A.AddInteraction(j, j+3, beta*(C2(i+1) + var_eq.epsilon*v(j+2)));
          A.AddInteraction(j, j-2, -beta*var_eq.epsilon*v(j-1));
          A.AddInteraction(j, j-1, -beta*(C2(i-1) + var_eq.epsilon*v(j-2)));
        }
      
      int N = 2*nb_points-1;
      A.AddInteraction(0, 2, beta*var_eq.epsilon*v(3));
      A.AddInteraction(0, 3, beta*(C2(1) + var_eq.epsilon*v(2)));
      A.AddInteraction(0, N-1, -beta*var_eq.epsilon*v(N));
      A.AddInteraction(0, N, -beta*(C2(nb_points-1) + var_eq.epsilon*v(N-1)));
      
      A.AddInteraction(N-1, 0, beta*var_eq.epsilon*v(1));
      A.AddInteraction(N-1, 1, beta*(C2(0) + var_eq.epsilon*v(0)));
      A.AddInteraction(N-1, N-3, -beta*var_eq.epsilon*v(N-2));
      A.AddInteraction(N-1, N-2, -beta*(C2(nb_points-2) + var_eq.epsilon*v(N-3)));
      
      // part zeta_x + epsilon u u_x
      Real_wp coef = var_eq.epsilon/3;
      for (int i = 1; i < nb_points-1; i++)
        {
          int j = 2*i;
          A.AddInteraction(j+1, j+2, beta);
          A.AddInteraction(j+1, j-2, -beta);
          A.AddInteraction(j+1, j+3, beta*coef*(2.0*v(j+3) + v(j+1)));
          A.AddInteraction(j+1, j+1, beta*coef*(2.0*v(j+3) - v(j-1)));
          A.AddInteraction(j+1, j-1, beta*coef*(-2.0*v(j-1) - v(j+1)));
        }
      
      A.AddInteraction(1, 2, beta);
      A.AddInteraction(1, N-1, -beta);
      A.AddInteraction(1, 3, beta*coef*(2.0*v(3) + v(1)));
      A.AddInteraction(1, 1, beta*coef*(2.0*v(3) - v(N)));
      A.AddInteraction(1, N, beta*coef*(-2.0*v(N) - v(1)));
      
      A.AddInteraction(N, 0, beta);
      A.AddInteraction(N, N-3, -beta);
      A.AddInteraction(N, 1, beta*coef*(2.0*v(1) + v(N)));
      A.AddInteraction(N, N, beta*coef*(2.0*v(1) - v(N-2)));
      A.AddInteraction(N, N-2, beta*coef*(-2.0*v(N-2) - v(N)));
      
    }
    
    
    template<class Matrix1, class Matrix2>
    void CheckMatrix(const Matrix1& A, const Matrix2& B)
    {
      for (int j = 0; j < Xn.GetM(); j++)
        for (int i = 0; i < Xn.GetM(); i++)
          {
            Real_wp val_num = A(i, j);
            Real_wp val_exact = B(i, j);
            if (abs(val_num-val_exact) > 1e-4)
              {
                DISP(i); DISP(j);
                DISP(val_num); DISP(val_exact);
                // DF(i,j) = val_num;
                abort();
              }
          }
    }
    
    
    void Advance(Real_wp t, int iter_num)
    {
      Real_wp test = 1e30, test_prec = 2e30;
      Copy(Xn, Xn_next);
      int nb_iter = 0;
      Real_wp threshold = 10/dt*epsilon_machine*Norm2(Xn);
      while ((test > threshold)&&(test < test_prec))
        {
          test_prec = test;
          // calcul fonctionnelle g = M((z - z^n)/dt)  K( (z+z^n)/2)
          for (int i = 0; i < 2*nb_points; i++)
            {
              Xn_stiff(i) = 0.5*(Xn(i) + Xn_next(i));
              Xn_mass(i) = (Xn_next(i) - Xn(i))/dt;
            }
          
          rhs.Fill(0);
          ApplyMassOperator(1.0, Xn_mass, rhs);
          ApplyStiffnessOperator(1.0, Xn_stiff, rhs);
          test = Norm2(rhs);

          if (compute_differential)
            {
              mat_sparse.Zero();
              AddMassDifferential(1.0/dt, Xn_mass, mat_sparse);
              AddStiffnessDifferential(0.5, Xn_stiff, mat_sparse);
  */
              /*Matrix<Real_wp> DF_num(Xn.GetM(), Xn.GetM()); DF_num.Fill(0);
              VectReal_wp ZnTmp = Xn;
              VectReal_wp ProdZnp = Xn, ProdZnm = Xn;
              ProdZnp.Fill(0); ProdZnm.Fill(0);
              DF_num.Fill(0);
              Real_wp h = 1e-7;
              for (int j = 0; j < Xn.GetM(); j++)
                {
                  ZnTmp(j) += h;
                  for (int i = 0; i < 2*nb_points; i++)
                    {
                      Xn_stiff(i) = 0.5*(Xn(i) + ZnTmp(i));
                      Xn_mass(i) = (ZnTmp(i) - Xn(i))/dt;
                    }
                  
                  ProdZnp.Fill(0);
                  ApplyMassOperator(1.0, Xn_mass, ProdZnp);
                  ApplyStiffnessOperator(1.0, Xn_stiff, ProdZnp);
                  
                  ZnTmp(j) -= 2.0*h;
                  for (int i = 0; i < 2*nb_points; i++)
                    {
                      Xn_stiff(i) = 0.5*(Xn(i) + ZnTmp(i));
                      Xn_mass(i) = (ZnTmp(i) - Xn(i))/dt;
                    }
                  
                  ProdZnm.Fill(0);
                  ApplyMassOperator(1.0, Xn_mass, ProdZnm);
                  ApplyStiffnessOperator(1.0, Xn_stiff, ProdZnm);
                  
                  for (int i = 0; i < Xn.GetM(); i++)
                    DF_num(i, j) = (ProdZnp(i) - ProdZnm(i))/(2.0*h);
                  
                  ZnTmp(j) += h;
                }
                  
              //DF_num.WriteText("df_num.dat");
              //mat_sparse.WriteText("df_exact.dat");
              
              CheckMatrix(DF_num, mat_sparse);
              */
  /* mat_lu.HideMessages();
              GetLU(mat_sparse, mat_lu, true);
              compute_differential = false;
              nb_eval_diff++;
            }
          
          SolveLU(mat_lu, rhs);
          
          for (int i = 0; i < 2*nb_points; i++)
            Xn_next(i) -= rhs(i);
          
          nb_iter++;
        }
      // DISP(test);
      
      if (nb_iter > 20)
        compute_differential = true;
      
      Copy(Xn_next, Xn);
      
    }
    
    Real_wp GetEnergy()
    {
      Real_wp energy = Norm2(Xn);
      
      return sqrt(energy);
    }
    
  }; */
  
}

#define MONTJOIE_FILE_FINITE_DIFFERENCE_BOUSSINESQ_CXX
#endif
