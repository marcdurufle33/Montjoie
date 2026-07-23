#ifndef MONTJOIE_FILE_FINITE_DIFFERENCE_GREEN_NAGHDI_CXX

namespace Montjoie
{
  
  /* template<>
  class FiniteDifferenceScheme1D<GreenNaghdiEquation>
    : public FiniteDifferenceScheme1D_Base<GreenNaghdiEquation>
  {
  public :
    
    Vector<Real_wp> Un_demi, Un, Zn, AhUn, Xn, AhXn;
    Vector<Real_wp> oneMinusbetab, C4, Ones;
#ifdef SELDON_WITH_MUMPS
    MatrixMumps<Real_wp> mat_lu;
#endif
#ifdef SELDON_WITH_SUPERLU
    MatrixSuperLU<Real_wp> mat_lu;
#endif
#ifdef SELDON_WITH_UMFPACK
    MatrixUmfPack<Real_wp> mat_lu;
#endif
#ifdef SELDON_WITH_PASTIX
    MatrixPastix<Real_wp> mat_lu;
#endif
    Real_wp dt, alpha;
    
    FiniteDifferenceScheme1D() : FiniteDifferenceScheme1D_Base<GreenNaghdiEquation>()
    {
      dt = 0;
      alpha = 2.0/3;
    }
    
    const VectReal_wp& GetLastSolution(int n) const
    {
      if (n == 0)
	return Zn;
      
      return Un;
    }
    
    void InitTimeScheme(Real_wp dt_, Real_wp Tfinal)
    {
      dt = dt_;
      int N = nb_points;
      Un_demi.Reallocate(N); Un_demi.Fill(0);
      Un.Reallocate(N); Un.Fill(0);
      AhUn.Reallocate(N); AhUn.Fill(0);
      Zn.Reallocate(N); Zn.Fill(0);
      Xn.Reallocate(2*N); Xn.Fill(0);
      AhXn.Reallocate(2*N); AhXn.Fill(0);
      
      Real_wp dx = (xN - x0)/(nb_points-1);
      VectReal_wp x(nb_points);
      for (int i = 0; i < nb_points; i++)
	x(i) = x0 + i*dx;
      
      VectReal_wp b, c;
      var_eq.ComputeBottom(x0, xN, x, b, c);
      
      // on calcule c_{i+1/2} = 1/2(c_i + c_i+1)
      Vector<Real_wp> chalf(N);
      for (int i = 0; i < N-1; i++)
	chalf(i) = 0.5*(c(i) + c(i+1));
      
      chalf(N-1) = 0.5*(c(N-1) + c(0));
      
      oneMinusbetab.Reallocate(N);
      for (int i = 0; i < N; i++)
	oneMinusbetab(i) = 1.0 - var_eq.beta*b(i);
      
      C4.Reallocate(N);
      for (int i = 0; i < N; i++)
	C4(i) = pow(chalf(i), 4.0);
      
      Ones.Reallocate(N); Ones.Fill(1.0);
      
      // initial condition
      var_eq.ComputeInitialCondition(x0, xN, x, Xn);
      for (int i = 0; i < N; i++)
	{
	  Zn(i) = Xn(i);
	  Un(i) = Xn(i+N);
	}
      
      Copy(Un, Un_demi);
    }
    
    
    void Advance(Real_wp t, int iter_num)
    {
      int N = nb_points;
      // on calcule U_n+1/2 = 2 Un - Un-1/2
      Mlt(-1.0, Un_demi); Add(2.0, Un, Un_demi);
      
      // we consider linear part
      Matrix<Real_wp, General, ArrayRowSparse> M(2*N, 2*N);
      
      // terms coming from first equation
      Vector<Real_wp> ones(N); ones.Fill(1.0);
      AddMatrixD1(0.5*dt, ones, oneMinusbetab, M, 0, N);
      AddMatrixD1(0.5*dt*var_eq.epsilon, Un_demi, ones, M, 0, 0);
      MltMatrixD1(0.5*dt*var_eq.epsilon, Un_demi, AhUn);
      for (int i = 0; i < N; i++)
	M.AddInteraction(i, i, AhUn(i));
      
      // terms coming from second equation
      AddMatrixD1(0.5*dt, ones, ones, M, N, 0);
      AddMatrixD1(0.5*var_eq.epsilon*dt*alpha, Un_demi, ones, M, N, N);
      for (int i = 0; i < N; i++)
	M.AddInteraction(N+i, N+i, (1.0-alpha)*AhUn(i));
      
      // on calcule (I-M+mu/3 D2) X^n
      MltMatrixD2var(-var_eq.mu/3.0, C4, Un, AhUn);
      for (int i = 0; i < N; i++)
	{
	  Xn(i) = Zn(i);
	  Xn(i+N) = Un(i);
	  AhXn(i) = Zn(i);
	  AhXn(i+N) = AhUn(i) + Un(i);
	}
      MltAdd(-1.0, M, Xn, 1.0, AhXn);    
      
      // on factorise (I+M+mu/3 D2)
      AddMatrixD2var(-var_eq.mu/3.0, C4, M, N, N);
      for (int i = 0; i < 2*N; i++)
	M.AddInteraction(i, i, 1.0);
      
      mat_lu.HideMessages();
      GetLU(M, mat_lu);
      
      // et on resout
      SolveLU(mat_lu, AhXn);
      for (int i = 0; i < N; i++)
	{
	  Zn(i) = AhXn(i);
	  Un(i) = AhXn(i+N);
	}
      
    }
    
    Real_wp GetEnergy()
    {
      Real_wp energy = 0;
      for (int i = 0; i < nb_points; i++)
	energy += Un(i)*Un(i) + Zn(i)*Zn(i);
      
      return sqrt(energy);
    }
    
  };
  */
}

#define MONTJOIE_FILE_FINITE_DIFFERENCE_GREEN_NAGHDI_CXX
#endif
