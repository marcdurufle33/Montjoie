#ifndef MONTJOIE_FILE_FINITE_DIFFERENCE_1D_CXX

namespace Montjoie
{
  template<class TypeEquation>
  class FiniteDifferenceScheme1D
  {
  };
  
  
  //! base class for the resolution of non-linear equations with finite-difference methods
  template<class TypeEquation>
  class FiniteDifferenceScheme1D_Base  : public VirtualOdeSystem<Real_wp>, public InputDataProblem_Base
  {
  public :
    //! available boundary conditions
    enum {DIRICHLET, PERIODIC};
    
    //! equation considered
    TypeEquation var_eq;
        
  protected :
    //! time step and space step
    Real_wp dt, dx;
    
    //! type of time scheme
    int time_scheme;
    
    //! order of time scheme
    int order_time_scheme;
    
    //! talezer iterator
    Talezer_Iterator<Real_wp> talezer_scheme;
    
    //! Runge-Kutta iterator
    RungeKutta_Iterator<Real_wp> RK_scheme;
    
    //! Implicit Runge-Kutta iterator
    GaussRungeKutta_Iterator<Real_wp> GaussRK_scheme;
    
    //! Low-storage Runge-Kutta iterator
    LowStorageRK_Iterator<Real_wp> LowRK_scheme;
    
    //! Adams-Bashforth iterator
    AdamsBashforth_Moulton_Iterator<Real_wp> AB_scheme;
    
    //! if true the differential has to be computed and factorized
    bool compute_differential;
    
    //! number of times the differential matrix is evaluated
    int nb_eval_diff;
    
    //! iterates of time scheme in real domain
    VectReal_wp Ur_time, Ur_prev, Ur_next;

    //! temporary variables used by time schemes
    Vector<Real_wp> Un_demi, Un, AhUn, Un_next, rhs;
    VectReal_wp Un_demi_tmp, Un_tmp, Vn_tmp, Prod_tmp;
    
    //! interval of computation
    Real_wp x0, xN;
    
    //! number of points in the interval of computation
    /*!
      - If a periodic condition is set, the last extremity is not discretized
      - If Dirichlet conditions are set, the two extremities are not discretized      
     */
    int nb_points;
    
    //! conditions at the two extremities of the interval
    int type_condition_left, type_condition_right;
    
    //! interpolation grid (points on which the solution is computed
    GridInterpolation<Dimension1> grid;
    
  public :
    //! default constructor
    FiniteDifferenceScheme1D_Base()
    {
      dt = 0;
      dx = 0;
      time_scheme = TimeSchemeEnum::PREDICTOR;
      nb_eval_diff = 0;
      compute_differential = true;

      type_condition_left = PERIODIC;
      type_condition_right = PERIODIC;
      nb_points = 0;
    }    
    
    
    //! returns instance of the leaf class
    const FiniteDifferenceScheme1D<TypeEquation>& GetLeafClass() const
    {
      return static_cast<const FiniteDifferenceScheme1D<TypeEquation>& >(*this);
    }


    //! returns instance of the leaf class
    FiniteDifferenceScheme1D<TypeEquation>& GetLeafClass()
    {
      return static_cast<FiniteDifferenceScheme1D<TypeEquation>& >(*this);
    }
    
    
    //! returns the space step \Delta x
    const Real_wp& GetSpaceStep() const
    {
      return dx;
    }
    
    
    int GetNbEvaluationDifferential() const
    {
      return nb_eval_diff;
    }

    
    bool DifferentialToBeComputed() const
    {
      return compute_differential;
    }


    void SetDifferentialToCompute(bool to_compute)
    {
      compute_differential = to_compute;
    }
    
    
    void IncrementDifferential()
    {
      nb_eval_diff++;
    }
    
    
    //! modification of parameters because of a line of the data file
    void SetInputData(const string& keyword, const Vector<string>& param)
    {
    }
    

    //! sets the interval of computation
    void SetInterval(const Real_wp& x0_, const Real_wp& xN_, int n)
    {
      x0 = x0_;
      xN = xN_;
      nb_points = n;
      
      // computation of dx
      if (type_condition_left == PERIODIC)
        dx = (xN - x0)/nb_points;
      else
        {
          if ((type_condition_left == DIRICHLET) && (type_condition_right == DIRICHLET))
            dx = (xN - x0)/(nb_points+1);   
          else
            {
              cout << "not implemented" << endl;
              abort();
            }
        }
      
      // allocation of intermediary vectors
      Un_tmp.Reallocate(n);
      Un_demi_tmp.Reallocate(n);
      Vn_tmp.Reallocate(n);
      Prod_tmp.Reallocate(n);
      Un_tmp.Fill(0);
      Un_demi_tmp.Fill(0);
      Vn_tmp.Fill(0);
      Prod_tmp.Fill(0);
    }
    
    
    //! sets the interval of computation
    void SetInterval(const Real_wp& x0_, const Real_wp& xN_, int n, int r)
    {
      SetInterval(x0_, xN_, n*(r+1));
    }
    
    
    //! sets the boundary condition at the left and right extremity
    void SetBoundaryCondition(int cond_left, int cond_right)
    {
      if ( ((cond_left == PERIODIC) && (cond_right != PERIODIC)) ||
           ((cond_right == PERIODIC) && (cond_left != PERIODIC)) )
        {
          cout << "Case Forbidden, periodic conditions should be enforced "
               << "on the two extremities" << endl;
          abort();
        }
      
      type_condition_left = cond_left;
      type_condition_right = cond_right;
    }
    
    
    //! computes the interpolation of the current solution on the interpolation grid
    /*!
      \param[in] n current iterate number
      \param[out] Un values of the solution on points of the interpolation grid
     */
    void GetInterpolateUn(int n, VectReal_wp& Un) const
    {
      GetInterpolateUn(n, Un, this->GetIterate());
    }
    
    
    //! returns the current solution
    const VectReal_wp& GetIterate() const
    {
      switch (time_scheme)
	{
	case TimeSchemeEnum::LEAP_FROG :
	case TimeSchemeEnum::LEAP_FROG_EXACT :
	  return Ur_prev;
	case TimeSchemeEnum::RUNGE_KUTTA :
	  return RK_scheme.GetIterate();
	case TimeSchemeEnum::ADAMS_BASHFORTH :
	  return AB_scheme.GetIterate();
	case TimeSchemeEnum::TALEZER :
	  return talezer_scheme.GetIterate();
	case TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE :
	  return LowRK_scheme.GetIterate();
        case TimeSchemeEnum::CRANK_NICOLSON :
          return Un;
        case TimeSchemeEnum::PREDICTOR :
          return Un;
	case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
	  return GaussRK_scheme.GetIterate();
	}
      
      return Ur_prev;
    }
    
    
    //! returns the current solution
    const VectReal_wp& GetIterateReal() const
    {
      return GetIterate();
    }
    
    
    //! Computes the positions where unknowns are discretized
    /*!
      This routine also initializes the space step dx
     */
    void GetSubdivisionSpace(VectReal_wp& x)
    {
      if (type_condition_left == PERIODIC)
        {
          x.Reallocate(nb_points);
          for (int i = 0; i < nb_points; i++)
            x(i) = x0 + i*dx;
        }
      else
        {
          if ((type_condition_left == DIRICHLET) && (type_condition_right == DIRICHLET))
            {
              x.Reallocate(nb_points);
              for (int i = 0; i < nb_points; i++)
                x(i) = x0 + (i+1)*dx;
            }
          
          //dx = (xN - x0)/(nb_points-1);
          //x.Reallocate(nb_points);
          //for (int i = 0; i < nb_points; i++)
          //x(i) = x0 + i*dx;
        }
    }
    
    
    //! initialisation of the spatial scheme
    /*!
      \param[in] order dummy parameter required by analogy with LDG class
      In leaf classes, this method is overloaded
      in order to compute physical coefficients
    */
    void SetOrder(int order)
    {
      // subdivision x_i in space
      VectReal_wp x;
      GetSubdivisionSpace(x);
    }
    
    
    //! initialization of the interpolation grid
    /*!
      \param[in] xmin first extremity of the interval
      \param[in] xmax second extremity of the interval
     */
    void InitGrid(const Real_wp& xmin, const Real_wp& xmax, int Nd)
    {
      Mesh<Dimension1> mesh;
      if (type_condition_left == PERIODIC)
        mesh.CreateRegularMesh(x0, xN, nb_points+1, 1);
      else
        {
          if ((type_condition_left == DIRICHLET) && (type_condition_right == DIRICHLET))
            mesh.CreateRegularMesh(x0, xN, nb_points+2, 1);
          else
            {
              cout << "case not implemented" << endl;
              abort();
            }
        }
      
      grid.Init(xmin, xmax, Nd);
      grid.LocalizePoints(mesh);
    }
    
    
    //! Computes the interpolation of xsol on the interpolation grid
    /*!
      \param[in] n iterate number
      \param[in] xsol vector to interpolate
      \param[out] Un interpolation of xsol on the interpolation grid
     */
    void GetInterpolateUn(int n, VectReal_wp& Un, const VectReal_wp& xsol) const
    {
      Un.Reallocate(grid.GetNbPointsGrid());
      for (int i = 0; i < Un.GetM(); i++)
	{
	  int ne = grid.GetElementNumber(i);
	  if (ne >= 0)
	    {
	      Real_wp x = grid.GetLocalCoordinate(i);
              if ((type_condition_left == DIRICHLET) && (type_condition_right == DIRICHLET))
                {
                  if (ne == 0)
                    Un(i) = x*xsol(0);
                  else if (ne == nb_points)
                    Un(i) = (1.0-x)*xsol(nb_points-1);
                  else
                    Un(i) = (1.0-x)*xsol(ne-1) + x*xsol(ne);
                }
              else
                {
                  if (ne == nb_points-1)
                    Un(i) = (1.0-x)*xsol(ne) + x*xsol(0);
                  else
                    Un(i) = (1.0-x)*xsol(ne) + x*xsol(ne+1);
                }
	    }
	  else
	    Un(i) = 0.0;
	}
    }
    
    
    //! Computes A = A + alpha D1
    /*!
      \param[in] alpha coefficient used in the combination
      \param[in] coef_row coefficient applied on rows of the matrix D1
      \param[in] coef_col coefficient applied on rows of the matrix D1
      \param[inout] A sparse matrix on which alpha D1 is added
      \param[in] m offset applied to row numbers
      \param[in] n offset applied to column numbers
      D1 is here the matrix of the discretization of du/dx
      D1 U = (U_{i+1} - U_{i-1})/(2 dx)
      The coefficient 1.0 should be provided in alpha if you wish to
      have the discretization of u_x. coef_col and coef_row, are such that you can consider
      the discretization of coef_row (coef_col U)_x
     */
    template<class MatrixSparse>
    void AddMatrixD1(Real_wp alpha, const VectReal_wp& coef_row,
                     const VectReal_wp& coef_col, MatrixSparse& mat,
                     int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      alpha *= 0.5/dx;
      int N = nb_points;
      for (int i = 1; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i-1), -alpha*coef_row(i)*coef_col(i-1));
      
      for (int i = 0; i < N-1; i++)
	mat.AddInteraction(m + incx*i, n+incy*(i+1), alpha*coef_row(i)*coef_col(i+1));
      
      if (type_condition_left == DIRICHLET)
	{
	}
      else if (type_condition_left == PERIODIC)
	{
	  mat.AddInteraction(m, n + incy*(N-1), -alpha*coef_row(0)*coef_col(N-1));
	  mat.AddInteraction(m + incx*(N-1), n, alpha*coef_row(N-1)*coef_col(0));
	}
    }
    
    
    //! Matrix vector product Vn = alpha D1 Un
    /*!
      \param[in] alpha coefficient
      \param[in] Un vector
      \param[inout] Vn vector
      D1 is here the matrix of the discretization of du/dx
      D1 U = (U_{i+1} - U_{i-1})/(2 dx)
      The coefficient 1 should be provided in alpha if you wish to
      have the discretization of u_x.
     */
    void MltMatrixD1(Real_wp alpha, const VectReal_wp& Un, VectReal_wp& Vn)
    {
      alpha *= 0.5/dx;
      int N = nb_points;
      for (int i = 1; i < N-1; i++)
	Vn(i) = alpha*(Un(i+1) - Un(i-1));
      
      if (type_condition_left == DIRICHLET)
	{
	  Vn(0) = alpha*Un(1);	
	}
      else if (type_condition_left == PERIODIC)
	{
	  Vn(0) = alpha*(Un(1) - Un(N-1));
	  Vn(N-1) = alpha*(Un(0) - Un(N-2));
	}
      
      if (type_condition_right == DIRICHLET)
	Vn(N-1) = -alpha*Un(N-2);
    }
    
    
    //! Adds matrix coming from scheme c_{i+1/2} u_{i+1} - c_{i-1/2} u_{i-1}
    /*!
      \param[in] alpha coefficient
      \param[in] chalf evaluation of c_{i+1/2}
      \param[inout] mat replaced by mat + alpha * D1
      \param[in] m offset applied to row numbers
      \param[in] n offset applied to column numbers
      where D1 is the matrix such that
      (D1 U)_i = (c_{i+1/2} u_{i+1} - c_{i-1/2} u_{i-1})/(2 dx)
      it is the discretisation of c u_x + c_x/2 u
    */
    template<class MatrixSparse>
    void AddMatrixD1var(Real_wp alpha, const VectReal_wp& chalf,
                        MatrixSparse& mat, int m = 0, int n = 0,
                        int incx = 1, int incy = 1)
    {
      alpha *= 0.5/dx;
      int N = nb_points;
      for (int i = 1; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i-1), -alpha*chalf(i-1));
      
      for (int i = 0; i < N-1; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i+1), alpha*chalf(i));
      
      if (type_condition_left == DIRICHLET)
	{
	}
      else if (type_condition_left == PERIODIC)
	{
	  mat.AddInteraction(m, n + incy*(N-1), -alpha*chalf(N-1));
	  mat.AddInteraction(m + incx*(N-1), n, alpha*chalf(N-1));
	}
    }


    //! Adds matrix coming from scheme c_{i+1/2} u_{i+1} - c_{i-1/2} u_{i-1}
    /*!
      \param[in] alpha coefficient
      \param[in] coef_row coefficient applied to rows
      \param[in] coef_col coefficient applied to columns
      \param[in] chalf evaluation of c_{i+1/2}
      \param[inout] mat replaced by mat + alpha * D1
      \param[in] m offset applied to row numbers
      \param[in] n offset applied to column numbers
      where D1 is the matrix such that
      (D1 U)_i = (c_{i+1/2} u_{i+1} - c_{i-1/2} u_{i-1})/(2 dx)
      it is the discretisation of c u_x + c_x/2 u
    */
    template<class MatrixSparse>
    void AddMatrixD1var(Real_wp alpha, const VectReal_wp& coef_row,
                        const VectReal_wp& coef_col, const VectReal_wp& chalf,
                        MatrixSparse& mat, int m = 0, int n = 0,
                        int incx = 1, int incy = 1)
    {
      alpha *= 0.5/dx;
      int N = nb_points;
      for (int i = 1; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i-1),
                           -alpha*coef_row(i)*coef_col(i-1)*chalf(i-1));
      
      for (int i = 0; i < N-1; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i+1), alpha*coef_row(i)*coef_col(i+1)*chalf(i));
      
      if (type_condition_left == DIRICHLET)
	{
	}
      else if (type_condition_left == PERIODIC)
	{
	  mat.AddInteraction(m, n + incy*(N-1), -alpha*coef_row(0)*coef_col(N-1)*chalf(N-1));
	  mat.AddInteraction(m + incx*(N-1), n, alpha*coef_row(N-1)*coef_col(0)*chalf(N-1));
	}
    }
    

    //! Matrix vector product Vn = alpha D1 Un
    /*!
      \param[in] alpha coefficient
      \param[in] Un vector
      \param[inout] Vn vector
      D1 is here the matrix of the discretization of c du/dx + c_x/2 u
      D1 U = (c_{i+1/2} U_{i+1} - c_{i-1/2} U_{i-1})/(2 dx)
      The coefficient 1 should be provided in alpha if you wish to
      have the discretization of c u_x + c_x/2 u.
     */
    void MltMatrixD1var(Real_wp alpha, const VectReal_wp& chalf,
                        const VectReal_wp& Un, VectReal_wp& Vn)
    {
      alpha *= 0.5/dx;
      int N = nb_points;
      for (int i = 1; i < N-1; i++)
	Vn(i) = alpha*(chalf(i)*Un(i+1) - chalf(i-1)*Un(i-1));
      
      if (type_condition_left == DIRICHLET)
	{
	  Vn(0) = alpha*chalf(0)*Un(1);	
	}
      else if (type_condition_left == PERIODIC)
	{
	  Vn(0) = alpha*(chalf(0)*Un(1) - chalf(N-1)*Un(N-1));
	  Vn(N-1) = alpha*(chalf(N-1)*Un(0) - chalf(N-2)*Un(N-2));
	}
      
      if (type_condition_right == DIRICHLET)
	Vn(N-1) = -alpha*chalf(N-2)*Un(N-2);
    }

    
    //! Adds to mat alpha D2
    /*!
      \param[in] alpha coefficient
      \param[in] coef_row rows of D2 are multiplied by these coefficients
      \param[in] coef_col columns of D2 are multiplied by these coefficients
      \param[inout] mat replaced by mat + alpha D2
      \param[in] m offset applied to row numbers
      \param[in] n offset applied to column numbers
      D2 is a matrix such that (D2 U)_i = coef_row_i (U_{i+1} - 2 U_i + U_{i-1})/dx^2
      In order to have a discretization of coef_row u_xx, you have to set alpha = 1
     */
    template<class MatrixSparse>
    void AddMatrixD2(Real_wp alpha, const VectReal_wp& coef_row,
                     const VectReal_wp& coef_col, MatrixSparse& mat, int m = 0, int n = 0,
                     int incx = 1, int incy = 1)
    {
      alpha *= 1.0/(dx*dx);
      int N = nb_points;
      for (int i = 0; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*i, -2.0*alpha*coef_row(i)*coef_col(i));
      
      for (int i = 1; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i-1), alpha*coef_row(i)*coef_col(i-1));
      
      for (int i = 0; i < N-1; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i+1), alpha*coef_row(i)*coef_col(i+1));
      
      if (type_condition_left == DIRICHLET)
	{
	}
      else if (type_condition_left == PERIODIC)
	{
	  mat.AddInteraction(m, n + incy*(N-1), alpha*coef_row(0)*coef_col(N-1));
	  mat.AddInteraction(m + incx*(N-1), n, alpha*coef_row(N-1)*coef_col(0));
	}
    }
    
    
    //! Adds to mat alpha D2
    /*!
      \param[in] alpha coefficient
      \param[in] c coefficients used in the expression of D2
      \param[inout] mat replaced by mat + alpha D2
      \param[in] m offset applied to row numbers
      \param[in] n offset applied to column numbers
      c contains values of c(x_i + dx/2)      
      D2 is a matrix used for the discretisation of d/dx( c d/dx)
      D2 is associated with the scheme c_{i+1/2} (u_i+1 - u_i) - c_{i-1/2}(u_i - u_i-1) / dx^2
      In order to have a discretization of this operator, you have to set alpha = 1
     */
    template<class MatrixSparse>
    void AddMatrixD2var(Real_wp alpha, const VectReal_wp& c, MatrixSparse& mat,
                        int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      alpha *= 1.0/(dx*dx);
      int N = nb_points;
      for (int i = 1; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*i, -alpha*(c(i-1)+c(i)));
      
      for (int i = 1; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i-1), alpha*c(i-1));
      
      for (int i = 0; i < N-1; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i+1), alpha*c(i));
      
      if (type_condition_left == DIRICHLET)
	{
	  mat.AddInteraction(m, n, -2.0*alpha*c(0));
	}
      else if (type_condition_left == PERIODIC)
	{
	  mat.AddInteraction(m, n, -alpha*(c(0)+c(N-1)));
	  mat.AddInteraction(m, n + incy*(N-1), alpha*c(N-1));
	  mat.AddInteraction(m + incx*(N-1), n, alpha*c(N-1));
	}
    }
    
    //! Matrix vector product Vn = alpha D2 Un
    /*!
      \param[in] alpha coefficient
      \param[in] Un vector
      \param[inout] Vn vector
      D2 is here the matrix of the discretization of d^2 u/dx^2
      D2 U = (U_{i+1} - 2 U_i + U_{i-1}) / dx^2
      The coefficient 1 should be provided in alpha if you wish to
      have the discretization of u_xx.
     */    
    void MltMatrixD2(Real_wp alpha, const VectReal_wp& Un, VectReal_wp& Vn)
    {
      alpha *= 1.0/(dx*dx);
      int N = nb_points;
      for (int i = 1; i < N-1; i++)
	Vn(i) = alpha*(-2.0*Un(i) + Un(i+1) + Un(i-1));
      
      if (type_condition_left == DIRICHLET)
	{
	  Vn(0) = alpha*(-2.0*Un(0) + Un(1));
	}
      else if (type_condition_left == PERIODIC)
	{
	  Vn(0) = alpha*(-2.0*Un(0) + Un(1) + Un(N-1));
	  Vn(N-1) = alpha*(-2.0*Un(N-1) + Un(0) + Un(N-2));
	}
      
      if (type_condition_right == DIRICHLET)
	Vn(N-1) = alpha*(-2.0*Un(N-1) + Un(N-2));
      
    }
    
    
    //! Matrix vector product Vn = alpha D2 Un
    /*!
      \param[in] alpha coefficient
      \param[in] c physical coefficients
      \param[in] Un vector
      \param[inout] Vn vector
      D2 is here the matrix of the discretization of d/dx ( c du/dx)
      D2 U = (c_{i+1/2} U_{i+1} - (c_{i+1/2} + c_{i-1/2}) U_i + c_{i-1/2} U_{i-1}) / dx^2
      The coefficient 1 should be provided in alpha if you wish to
      have the discretization of d/dx ( c du/dx).
     */    
    void MltMatrixD2var(Real_wp alpha, const VectReal_wp& c,
                        const VectReal_wp& Un, VectReal_wp& Vn)
    {
      alpha *= 1.0/(dx*dx);
      int N = nb_points;
      for (int i = 1; i < N-1; i++)
	Vn(i) = alpha*(-(c(i) + c(i-1))*Un(i) + c(i)*Un(i+1) + c(i-1)*Un(i-1));
      
      if (type_condition_left == DIRICHLET)
	{
	  Vn(0) = alpha*(-2.0*c(0)*Un(0) + c(0)*Un(1));
	}
      else if (type_condition_left == PERIODIC)
	{
	  Vn(0) = alpha*(-(c(0)+c(N-1))*Un(0) + c(0)*Un(1) + c(N-1)*Un(N-1));
	  Vn(N-1) = alpha*(-(c(N-1)+c(N-2))*Un(N-1) + c(N-1)*Un(0) + c(N-2)*Un(N-2));
	}
      
      if (type_condition_right == DIRICHLET)
	Vn(N-1) = alpha*(-(c(N-1)+c(N-2))*Un(N-1) + c(N-2)*Un(N-2));
      
    }
    
    
    //! Adds to mat alpha D3
    /*!
      \param[in] alpha coefficient
      \param[in] c coefficients used in the expression of D3
      \param[inout] mat replaced by mat + alpha D3
      \param[in] m offset applied to row numbers
      \param[in] n offset applied to column numbers
      D3 is a matrix used for the discretisation of d^3/dx^3  (if c = 1)
      For c different from 1, D3 is the discretization of
      c u_xxx + 3/2 c_x u_xx + 3/4 c_xx u_x + 1/8 c_xxx u
      this operator is conservative whatever the values of c
      In order to have a discretization of this operator, you have to set alpha = 1.0
     */
    template<class MatrixSparse>
    void AddMatrixD3var(Real_wp alpha, const VectReal_wp& c,
                        MatrixSparse& mat, int m = 0, int n = 0,
                        int incx = 1, int incy = 1)
    {
      alpha *= 0.5/(dx*dx*dx);
      int N = nb_points;
      for (int i = 2; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i-2), -alpha*c(i-1));
      
      for (int i = 1; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i-1), alpha*(c(i) + c(i-1)));
      
      for (int i = 0; i < N-1; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i+1), -alpha*(c(i) + c(i+1)));
      
      for (int i = 0; i < N-2; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i+2), alpha*c(i+1));
      
      if (type_condition_left == PERIODIC)
	{
	  mat.AddInteraction(m, n + incy*(N-2), -alpha*c(N-1));
	  mat.AddInteraction(m + incx, n + incy*(N-1), -alpha*c(0));
	  
	  mat.AddInteraction(m, n + incy*(N-1), alpha*(c(0) + c(N-1)));
	  mat.AddInteraction(m + incx*(N-1), n, -alpha*(c(N-1) + c(0)));
	  
	  mat.AddInteraction(m + incx*(N-1), n + incy, alpha*c(0));
	  mat.AddInteraction(m + incx*(N-2), n, alpha*c(N-1));
	}
      else
        {
          cout << "not implemented" << endl;
          abort();
        }
    }


    //! Matrix vector product Vn = alpha D3 Un
    /*!
      \param[in] alpha coefficient
      \param[in] c physical coefficients
      \param[in] Un vector
      \param[inout] Vn vector
      D3 is a matrix used for the discretisation of d^3/dx^3  (if c = 1)
      For c different from 1, D3 is the discretization of 
      c u_xxx + 3/2 c_x u_xx + 3/4 c_xx u_x + 1/8 c_xxx u
      this operator is conservative whatever the values of c
      In order to have a discretization of this operator, you have to set alpha = 1.0
     */    
    void MltMatrixD3var(Real_wp alpha, const VectReal_wp& c,
                        const VectReal_wp& Un, VectReal_wp& Vn)
    {
      alpha *= 0.5/(dx*dx*dx);
      int N = nb_points;
      for (int i = 2; i < N-2; i++)
	Vn(i) = alpha*(-c(i-1)*Un(i-2) + (c(i)+c(i-1))*Un(i-1)
                       - (c(i)+c(i+1))*Un(i+1) + c(i+1)*Un(i+2));
      
      if (type_condition_left == PERIODIC)
	{
          Vn(N-2) = alpha*(-c(N-3)*Un(N-4) + (c(N-2)+c(N-3))*Un(N-3) 
                           - (c(N-2)+c(N-1))*Un(N-1) + c(N-1)*Un(0));
          Vn(N-1) = alpha*(-c(N-2)*Un(N-3) + (c(N-1)+c(N-2))*Un(N-2) 
                           - (c(N-1)+c(0))*Un(0) + c(0)*Un(1));
          Vn(0) = alpha*(-c(N-1)*Un(N-2) + (c(0)+c(N-1))*Un(N-1) - (c(0)+c(1))*Un(1) + c(1)*Un(2));
          Vn(1) = alpha*(-c(0)*Un(N-1) + (c(1)+c(0))*Un(0) - (c(1)+c(2))*Un(2) + c(2)*Un(3));
	}
      else
        {
          cout << "not implemented" << endl;
          abort();
        }
    }

    
    //! Adds to mat alpha D3
    /*!
      \param[in] alpha coefficient
      \param[in] coef_row coefficients applied to rows of the matrix
      \param[in] coef_col coefficients applieds to columns of the matrix
      \param[inout] mat replaced by mat + alpha D3
      \param[in] m offset applied to row numbers
      \param[in] n offset applied to column numbers
      D3 is a matrix used for the discretisation of d^3/dx^3  (if coef_row = 1 and coef_col = 1)
      For coefficients different from one, it should discretize coef_row d^3 ( coef_col u) / dx^3
      In order to have a discretization of this operator, you have to set alpha = 1
     */    
    template<class MatrixSparse>
    void AddMatrixD3(Real_wp alpha, const VectReal_wp& coef_row,
                     const VectReal_wp& coef_col,
		     MatrixSparse& mat, int m = 0, int n = 0,
                     int incx = 1, int incy = 1)
    {
      alpha *= 0.5/(dx*dx*dx);
      int N = nb_points;
      for (int i = 2; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i-2), -alpha*coef_row(i)*coef_col(i-2));
      
      for (int i = 1; i < N; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i-1), 2.0*alpha*coef_row(i)*coef_col(i-1));
      
      for (int i = 0; i < N-1; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i+1), -2.0*alpha*coef_row(i)*coef_col(i+1));
      
      for (int i = 0; i < N-2; i++)
	mat.AddInteraction(m + incx*i, n + incy*(i+2), alpha*coef_row(i)*coef_col(i+2));
      
      if (type_condition_left == PERIODIC)
	{
	  mat.AddInteraction(m, n + incy*(N-2), -alpha*coef_row(0)*coef_col(N-2));
	  mat.AddInteraction(m + incx, n + incy*(N-1), -alpha*coef_row(1)*coef_col(N-1));
	  
	  mat.AddInteraction(m, n + incy*(N-1), 2.0*alpha*coef_row(0)*coef_col(N-1));
	  mat.AddInteraction(m + incx*(N-1), n, -2.0*alpha*coef_row(N-1)*coef_col(0));
	  
	  mat.AddInteraction(m + incx*(N-1), n + incy, alpha*coef_row(N-1)*coef_col(1));
	  mat.AddInteraction(m + incx*(N-2), n, alpha*coef_row(N-2)*coef_col(0));
	}
      else
        {
          cout << "not implemented" << endl;
          abort();
        }
    }


    //! Matrix vector product Vn = alpha D3 Un
    /*!
      \param[in] alpha coefficient
      \param[in] Un vector
      \param[inout] Vn vector
      D3 is here the matrix of the discretization of d^3 u/dx^3
      D3 U = (U_{i+2} - 2 U_{i+1} + 2 U_{i-1} - U_{i-2}) / (2 dx^3)
      The coefficient 1 should be provided in alpha if you wish to
      have the discretization of u_xxx.
     */    
    void MltMatrixD3(Real_wp alpha, const VectReal_wp& Un, VectReal_wp& Vn)
    {
      alpha *= 0.5/(dx*dx*dx);
      int N = nb_points;
      for (int i = 2; i < N-2; i++)
	Vn(i) = alpha*(Un(i+2) - 2.0*(Un(i+1) - Un(i-1)) - Un(i-2));
      
      if (type_condition_left == PERIODIC)
	{
	  Vn(1) = alpha*(Un(3) - 2.0*(Un(2) - Un(0)) - Un(N-1));
          Vn(0) = alpha*(Un(2) - 2.0*(Un(1) - Un(N-1)) - Un(N-2));
	  Vn(N-1) = alpha*(Un(1) - 2.0*(Un(0) - Un(N-2)) - Un(N-3));
	  Vn(N-2) = alpha*(Un(0) - 2.0*(Un(N-1) - Un(N-3)) - Un(N-4));
	}
      else
        {
          cout << "not implemented" << endl;
          abort();
        }
    }

    
    //! computes Prod = Prod + alpha * g(Un, Un_demi)
    /*!
      where g(Un, Un_demi) = 1/(p+2) [ D1 Un . Un_demi^p + D1( Un Un_demi^p) ]
      g is a conservative discretization of the operator u^p u_x
     */
    void AddOperatorUpUx(Real_wp alpha, int p, 
                         const VectReal_wp& coef_row, const VectReal_wp& coef_col,
                         const VectReal_wp& Un, const VectReal_wp& Un_demi, VectReal_wp& Prod)
    {
      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = Un(i)*coef_col(i);
      
      if (p < 0)
        {
          cout << "Invalid value of p" << endl;
          abort();
        }
      else if (p == 0)
        {          
          MltMatrixD1(alpha, Un_tmp, Prod_tmp);
          
          for (int i = 0; i < nb_points; i++)
            Prod(i) += coef_row(i)*Prod_tmp(i);
          
          return;
        }
      
      // alpha is replaced by alpha (p+1) / (p+2)
      alpha *= Real_wp(1) / (p + 2);

      // Un_demi_tmp = (coef_col Un_demi)^p
      if (p == 1)
        for (int i = 0; i < nb_points; i++)
          Un_demi_tmp(i) = Un_demi(i)*coef_col(i);
      else if (p == 2)
        for (int i = 0; i < nb_points; i++)
          Un_demi_tmp(i) = square(Un_demi(i)*coef_col(i));
      else
        for (int i = 0; i < nb_points; i++)
          Un_demi_tmp(i) = pow(Un_demi(i)*coef_col(i), p);
      
      // part D1 Un . Un_demi^p
      MltMatrixD1(alpha, Un_tmp, Prod_tmp);
      
      for (int i = 0; i < nb_points; i++)
        Prod(i) += coef_row(i)*Prod_tmp(i)*Un_demi_tmp(i);
      
      // part D1(Un^p Un_demi)
      for (int i = 0; i < nb_points; i++)
        Vn_tmp(i) = Un_demi_tmp(i)*Un_tmp(i);
      
      MltMatrixD1(alpha, Vn_tmp, Prod_tmp);
      
      for (int i = 0; i < nb_points; i++)
        Prod(i) += coef_row(i)*Prod_tmp(i);
      
    }

    
    //! computes A = A + alpha * Dg(Un, Un_demi)
    /*!
      where g(Un, Un_demi) = 1/(p+2) [ D1 Un . Un_demi^p + D1( Un Un_demi^p) ]
      g is a conservative discretization of the operator u^p u_x
      Dg is the differential matrix of g with respect to Un
     */
    template<class MatrixSparse>
    void AddDifferentialHalfUpUx(Real_wp alpha, int p, 
                                 const VectReal_wp& coef_row, const VectReal_wp& coef_col,
                                 const VectReal_wp& Un_demi, MatrixSparse& mat,
                                 int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      if (p < 0)
        {
          cout << "Invalid value of p" << endl;
          abort();
        }
      else if (p == 0)
        {          
          AddMatrixD1(alpha, coef_row, coef_col, mat, m, n, incx, incy);
          
          return;
        }
      
      // alpha is replaced by alpha / (p+2)
      alpha *= Real_wp(1) / (p + 2);

      if (p == 1)
        for (int i = 0; i < nb_points; i++)
          Un_demi_tmp(i) = Un_demi(i)*coef_col(i);
      else if (p == 2)
        for (int i = 0; i < nb_points; i++)
          Un_demi_tmp(i) = square(Un_demi(i)*coef_col(i));
      else
        for (int i = 0; i < nb_points; i++)
          Un_demi_tmp(i) = pow(Un_demi(i)*coef_col(i), p);
      
      // part Un_demi^p D1
      for (int i = 0; i < nb_points; i++)
        Vn_tmp(i) = coef_row(i) * Un_demi_tmp(i);
      
      AddMatrixD1(alpha, Vn_tmp, coef_col, mat, m, n, incx, incy);
      
      // part D1( Un_demi^p )
      for (int i = 0; i < nb_points; i++)
        Vn_tmp(i) = coef_col(i) * Un_demi_tmp(i);
      
      AddMatrixD1(alpha, coef_row, Vn_tmp, mat, m, n, incx, incy);
    }
    

    //! computes A = A + alpha * Dg(Un)
    /*!
      where g(Un) = 1/(p+2) [ D1 Un . Un^p + D1( Un^{p+1} ) ]
      g is a conservative discretization of the operator u^p u_x
      Dg is the differential matrix of g with respect to Un
     */
    template<class MatrixSparse>
    void AddDifferentialUpUx(Real_wp alpha, int p, 
                             const VectReal_wp& coef_row, const VectReal_wp& coef_col,
                             const VectReal_wp& Un, MatrixSparse& mat, int m = 0, int n = 0,
                             int incx = 1, int incy = 1)
    {
      if (p < 0)
        {
          cout << "Invalid value of p" << endl;
          abort();
        }
      else if (p == 0)
        {          
          AddMatrixD1(alpha, coef_row, coef_col, mat, m, n, incx, incy);
          
          return;
        }
      
      // alpha is replaced by alpha / (p+2)
      alpha *= Real_wp(1) / (p + 2);

      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = Un(i)*coef_col(i);
      
      // Un_demi_tmp = Un^(p-1)
      if (p == 1)
        for (int i = 0; i < nb_points; i++)
          Un_demi_tmp(i) = 1.0;
      else if (p == 2)
        for (int i = 0; i < nb_points; i++)
          Un_demi_tmp(i) = Un_tmp(i);
      else
        for (int i = 0; i < nb_points; i++)
          Un_demi_tmp(i) = pow(Un_tmp(i), p-1);
      
      // part D1 Un  d/dun (Un^p) I
      MltMatrixD1(alpha, Un_tmp, Vn_tmp);
      
      Real_wp val;
      for (int i = 0; i < nb_points; i++)
        {
          val = p*coef_row(i)*Vn_tmp(i)*Un_demi_tmp(i)*coef_col(i);
          mat.AddInteraction(incx*i + m, incy*i + n, val);
        }
      
      // part Un^p D1 and Un_demi_tmp = Un^p
      for (int i = 0; i < nb_points; i++)
        {
          Un_demi_tmp(i) *= Un_tmp(i);
          Vn_tmp(i) = coef_row(i)*Un_demi_tmp(i);
        }
      
      AddMatrixD1(alpha, Vn_tmp, coef_col, mat, m, n, incx, incy);
      
      // part D1( d/dun( Un^p+1) )
      for (int i = 0; i < nb_points; i++)
        Vn_tmp(i) = (p+1)*coef_col(i) * Un_demi_tmp(i);
      
      AddMatrixD1(alpha, coef_row, Vn_tmp, mat, m, n, incx, incy);
    }

    
    //! computes Prod = Prod + alpha * g(Un, Un_demi)
    /*!
      where g(Un, Un_demi) = D1 Un . D2 Un_demi + D1( D2 Un_demi Un)
      g is a conservative discretization of the operator 2 u_xx u_x + u_xxx u
     */
    void AddOperatorUxxUx_UxxxU(Real_wp alpha,
                                const VectReal_wp& coef_row, const VectReal_wp& coef_col,
                                const VectReal_wp& Un, const VectReal_wp& Un_demi,
                                VectReal_wp& Prod)
    {
      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = Un(i)*coef_col(i);

      for (int i = 0; i < nb_points; i++)
        Un_demi_tmp(i) = Un_demi(i)*coef_col(i);
      
      // part D1 Un  D2 Un_demi
      MltMatrixD2(1.0, Un_demi_tmp, Vn_tmp);
      
      MltMatrixD1(alpha, Un_tmp, Prod_tmp);
      
      for (int i = 0; i < nb_points; i++)
        Prod(i) += coef_row(i)*Vn_tmp(i)*Prod_tmp(i);
      
      // part D1( D2 Un_demi Un)
      for (int i = 0; i < nb_points; i++)
        Vn_tmp(i) *= Un_tmp(i);
      
      MltMatrixD1(alpha, Vn_tmp, Prod_tmp);
      
      for (int i = 0; i < nb_points; i++)
        Prod(i) += coef_row(i)*Prod_tmp(i);      
    }
    

    //! computes A = A + alpha * Dg(Un, Un_demi)
    /*!
      where g(Un, Un_demi) = D1 Un . D2 Un_demi + D1( Un D2 Un_demi)
      g is a conservative discretization of the operator
      2 u_xx u_x + u_xxx u
      Dg is the differential matrix of g with respect to Un
     */
    template<class MatrixSparse>
    void AddDifferentialHalfUxxUx_UxxxU(Real_wp alpha,
                                        const VectReal_wp& coef_row, const VectReal_wp& coef_col,
                                        const VectReal_wp& Un_demi, MatrixSparse& mat,
                                        int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      for (int i = 0; i < nb_points; i++)
        Un_demi_tmp(i) = Un_demi(i)*coef_col(i);

      MltMatrixD2(1.0, Un_demi_tmp, Vn_tmp);
      
      // part   D2 Un_demi  D1
      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = coef_row(i)*Vn_tmp(i);
      
      AddMatrixD1(alpha, Un_tmp, coef_col, mat, m, n, incx, incy);
      
      // part D1( D2 Un_demi . )
      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = coef_col(i)*Vn_tmp(i);
      
      AddMatrixD1(alpha, coef_row, Un_tmp, mat, m, n, incx, incy);
    }
    
    
    //! computes A = A + alpha * Dg(Un)
    /*!
      where g(Un) = D1 Un . D2 Un + D1( Un D2 Un)
      g is a conservative discretization of the operator
      2 u_xx u_x + u_xxx u
      Dg is the differential matrix of g with respect to Un
     */
    template<class MatrixSparse>
    void AddDifferentialUxxUx_UxxxU(Real_wp alpha,
                                    const VectReal_wp& coef_row, const VectReal_wp& coef_col,
                                    const VectReal_wp& Un, MatrixSparse& mat,
                                    int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      int N = nb_points;
      Real_wp beta = alpha / (2.0*dx*dx*dx);
      for (int i = 0; i < N; i++)
        Un_tmp(i) = coef_col(i)*Un(i);
      
      for (int i = 2; i < N-2; i++)
        {
          mat.AddInteraction(m + incx*i, n + incy*(i+2),
                             beta*coef_row(i)*coef_col(i+2)*Un_tmp(i+1));
          
          mat.AddInteraction(m + incx*i, n + incy*(i+1),
                             beta*coef_row(i)*coef_col(i+1)
                             *(Un_tmp(i+2) -2.0*Un_tmp(i+1) - Un_tmp(i)) );
          
          mat.AddInteraction(m + incx*i, n + incy*i,
                             beta*coef_row(i)*coef_col(i)*(-Un_tmp(i+1) + Un_tmp(i-1))); 
          
          mat.AddInteraction(m + incx*i, n + incy*(i-1),
                             beta*coef_row(i)*coef_col(i-1)
                             *(Un_tmp(i) + 2.0*Un_tmp(i-1) - Un_tmp(i-2)) );
          
          mat.AddInteraction(m + incx*i, n + incy*(i-2),
                             -beta*coef_row(i)*coef_col(i-2)*Un_tmp(i-1));
        }
      
      mat.AddInteraction(m + incx*(N-2), n, beta*coef_row(N-2)*coef_col(0)*Un_tmp(N-1));
      mat.AddInteraction(m + incx*(N-2), n + incy*(N-1),
                         beta*coef_row(N-2)*coef_col(N-1)
                         *( Un_tmp(0) -2.0*Un_tmp(N-1) - Un_tmp(N-2)) );
      
      mat.AddInteraction(m + incx*(N-2), n + incy*(N-2),
                         beta*coef_row(N-2)*coef_col(N-2)*(-Un_tmp(N-1) + Un_tmp(N-3))); 
      
      mat.AddInteraction(m + incx*(N-2), n + incy*(N-3),
                         beta*coef_row(N-2)*coef_col(N-3)
                         *(Un_tmp(N-2) + 2.0*Un_tmp(N-3) - Un_tmp(N-4)) );
      
      mat.AddInteraction(m + incx*(N-2), n + incy*(N-4),
                         -beta*coef_row(N-2)*coef_col(N-4)*Un_tmp(N-3));

      mat.AddInteraction(m + incx*(N-1), n + incy, beta*coef_row(N-1)*coef_col(1)*Un_tmp(0));
      mat.AddInteraction(m + incx*(N-1), n, beta*coef_row(N-1)*coef_col(0)
                         *( Un_tmp(1) -2.0*Un_tmp(0) - Un_tmp(N-1)) );
      
      mat.AddInteraction(m + incx*(N-1), n + incy*(N-1), beta*coef_row(N-1)*coef_col(N-1)
                         *(-Un_tmp(0) + Un_tmp(N-2))); 
      
      mat.AddInteraction(m + incx*(N-1), n + incy*(N-2), beta*coef_row(N-1)*coef_col(N-2)
                         *(Un_tmp(N-1) + 2.0*Un_tmp(N-2) - Un_tmp(N-3)) );
      mat.AddInteraction(m + incx*(N-1), n + incy*(N-3), -beta*coef_row(N-1)*coef_col(N-3)
                         *Un_tmp(N-2));

      mat.AddInteraction(m, n + incy*2, beta*coef_row(0)*coef_col(2)*Un_tmp(1));
      mat.AddInteraction(m, n + incy, beta*coef_row(0)*coef_col(1)
                         *( Un_tmp(2) -2.0*Un_tmp(1) - Un_tmp(0)) );
      mat.AddInteraction(m, n, beta*coef_row(0)*coef_col(0)*(-Un_tmp(1) + Un_tmp(N-1))); 
      mat.AddInteraction(m, n + incy*(N-1), beta*coef_row(0)*coef_col(N-1)
                         *(Un_tmp(0) + 2.0*Un_tmp(N-1) - Un_tmp(N-2)) );
      mat.AddInteraction(m, n + incy*(N-2), -beta*coef_row(0)*coef_col(N-2)*Un_tmp(N-1));

      mat.AddInteraction(m + incx, n + incy*3, beta*coef_row(1)*coef_col(3)*Un_tmp(2));
      mat.AddInteraction(m + incx, n + incy*2, beta*coef_row(1)*coef_col(2)
                         *( Un_tmp(3) -2.0*Un_tmp(2) - Un_tmp(1)) );
      mat.AddInteraction(m + incx, n + incy, beta*coef_row(1)*coef_col(1)
                         *(-Un_tmp(2) + Un_tmp(0))); 
      mat.AddInteraction(m + incx, n, beta*coef_row(1)*coef_col(0)
                         *(Un_tmp(1) + 2.0*Un_tmp(0) - Un_tmp(N-1)) );
      mat.AddInteraction(m + incx, n + incy*(N-1), -beta*coef_row(1)*coef_col(N-1)*Un_tmp(0));
    }


    //! computes Prod = Prod + alpha * g(Un, Un_demi)
    /*!
      where g(Un, Un_demi) = Un_demi D3 Un
      g is a discretization of the operator u_xxx u
     */
    void AddOperatorUxxxU(Real_wp alpha,
                          const VectReal_wp& coef_row, const VectReal_wp& coef_col,
                          const VectReal_wp& Un, const VectReal_wp& Un_demi, VectReal_wp& Prod)
    {
      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = Un(i)*coef_col(i);

      for (int i = 0; i < nb_points; i++)
        Un_demi_tmp(i) = Un_demi(i)*coef_col(i);
      
      MltMatrixD3(1.0, Un_tmp, Vn_tmp);
      
      for (int i = 0; i < nb_points; i++)
        Prod(i) += alpha*coef_row(i)*Un_demi_tmp(i)*Vn_tmp(i);
    }

    
    //! computes A = A + alpha * Dg(Un, Un_demi)
    /*!
      where g(Un, Un_demi) = Un_demi . D1 Un
      g is a discretization of the operator u_xxx u
      Dg is the differential matrix of g with respect to Un
     */
    template<class MatrixSparse>
    void AddDifferentialHalfUxxxU(Real_wp alpha,
                                  const VectReal_wp& coef_row, const VectReal_wp& coef_col,
                                  const VectReal_wp& Un_demi, MatrixSparse& mat,
                                  int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      for (int i = 0; i < nb_points; i++)
        Un_demi_tmp(i) = Un_demi(i)*coef_col(i);
      
      for (int i = 0; i < nb_points; i++)
        Vn_tmp(i) = alpha*coef_row(i)*Un_demi_tmp(i);
      
      AddMatrixD3(1.0, Vn_tmp, coef_col, mat, m, n, incx, incy);
    }
    

    //! computes A = A + alpha * Dg(Un, Un)
    /*!
      where g(Un) = Un . D3 Un
      g is a discretization of the operator u_xxx u
      Dg is the differential matrix of g with respect to Un
     */
    template<class MatrixSparse>
    void AddDifferentialUxxxU(Real_wp alpha,
                              const VectReal_wp& coef_row, const VectReal_wp& coef_col,
                              const VectReal_wp& Un, MatrixSparse& mat, int m = 0, int n = 0,
                              int incx = 1, int incy = 1)
    {
      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = Un(i)*coef_col(i);
      
      MltMatrixD3(alpha, Un_tmp, Prod_tmp);
      
      for (int i = 0; i < nb_points; i++)
        {
          mat.AddInteraction(m+i, n+i, Prod_tmp(i)*coef_col(i)*coef_row(i));
          Vn_tmp(i) = alpha*coef_row(i)*Un_tmp(i);
        }
      
      AddMatrixD3(1.0, Vn_tmp, coef_col, mat, m, n, incx, incy);
    }

    
    //! computes Prod = Prod + alpha * g(Un, Un_demi)
    /*!
      where g(Un, Un_demi) = D1^v Un_demi D2 Un
      g is a discretization of the operator u_xx (c u_x + c_x/2 u)
     */
    void AddOperatorUxxUx(Real_wp alpha, const VectReal_wp& coef_row,
                          const VectReal_wp& coef_col, const VectReal_wp& c_half,
                          const VectReal_wp& Un, const VectReal_wp& Un_demi, VectReal_wp& Prod)
    {
      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = Un(i)*coef_col(i);

      for (int i = 0; i < nb_points; i++)
        Un_demi_tmp(i) = Un_demi(i)*coef_col(i);
      
      MltMatrixD1var(1.0, c_half, Un_demi_tmp, Vn_tmp);
      MltMatrixD2(1.0, Un_tmp, Prod_tmp);
      
      for (int i = 0; i < nb_points; i++)
        Prod(i) += alpha*coef_row(i)*Vn_tmp(i)*Prod_tmp(i);
    }

    
    //! computes A = A + alpha * Dg(Un, Un_demi)
    /*!
      where g(Un, Un_demi) = D1^v Un_demi . D2 Un
      g is a discretization of the operator u_xx (c u_x + c_x/2 u)
      Dg is the differential matrix of g with respect to Un
     */
    template<class MatrixSparse>
    void AddDifferentialHalfUxxUx(Real_wp alpha, const VectReal_wp& coef_row,
                                  const VectReal_wp& coef_col, const VectReal_wp& c_half,
                                  const VectReal_wp& Un_demi, MatrixSparse& mat,
                                  int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      for (int i = 0; i < nb_points; i++)
        Un_demi_tmp(i) = Un_demi(i)*coef_col(i);
      
      MltMatrixD1var(1.0, c_half, Un_demi_tmp, Vn_tmp);
      
      for (int i = 0; i < nb_points; i++)
        Vn_tmp(i) *= alpha*coef_row(i);
      
      AddMatrixD2(1.0, Vn_tmp, coef_col, mat, m, n, incx, incy);
    }
    

    //! computes A = A + alpha * Dg(Un, Un)
    /*!
      where g(Un) = D1^v Un . D2 Un
      g is a discretization of the operator u_xx (c u_x + c_x/2 u)
      Dg is the differential matrix of g with respect to Un
     */
    template<class MatrixSparse>
    void AddDifferentialUxxUx(Real_wp alpha, const VectReal_wp& coef_row,
                              const VectReal_wp& coef_col, const VectReal_wp& c_half,
                              const VectReal_wp& Un, MatrixSparse& mat,
                              int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = Un(i)*coef_col(i);
      
      MltMatrixD1var(1.0, c_half, Un_tmp, Vn_tmp);
      MltMatrixD2(1.0, Un_tmp, Prod_tmp);
      
      for (int i = 0; i < nb_points; i++)
        {
          Vn_tmp(i) *= coef_row(i);
          Prod_tmp(i) *= coef_row(i);
        }
      
      AddMatrixD2(alpha, Vn_tmp, coef_col, mat, m, n, incx, incy);
      AddMatrixD1var(alpha, Prod_tmp, coef_col, c_half, mat, m, n, incx, incy);
    }

    
    //! computes Prod = Prod + alpha * g(Un, Un_demi)
    /*!
      where g(Un, Un_demi) = D1^v Un_demi D1 Un
      g is a discretization of the operator u_x (c u_x + c_x/2 u)
     */
    void AddOperatorUxUx(Real_wp alpha, const VectReal_wp& coef_row,
                         const VectReal_wp& coef_col, const VectReal_wp& c_half,
                         const VectReal_wp& Un, const VectReal_wp& Un_demi, VectReal_wp& Prod)
    {
      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = Un(i)*coef_col(i);

      for (int i = 0; i < nb_points; i++)
        Un_demi_tmp(i) = Un_demi(i)*coef_col(i);
      
      MltMatrixD1var(1.0, c_half, Un_demi_tmp, Vn_tmp);
      MltMatrixD1(1.0, Un_tmp, Prod_tmp);
      
      for (int i = 0; i < nb_points; i++)
        Prod(i) += alpha*coef_row(i)*Vn_tmp(i)*Prod_tmp(i);
    }

    
    //! computes A = A + alpha * Dg(Un, Un_demi)
    /*!
      where g(Un, Un_demi) = D1^v Un_demi . D1 Un
      g is a discretization of the operator u_x (c u_x + c_x/2 u)
      Dg is the differential matrix of g with respect to Un
     */
    template<class MatrixSparse>
    void AddDifferentialHalfUxUx(Real_wp alpha, const VectReal_wp& coef_row,
                                 const VectReal_wp& coef_col, const VectReal_wp& c_half,
                                 const VectReal_wp& Un_demi, MatrixSparse& mat,
                                 int m = 0, int n = 0, int incx = 1, int incy = 1)
    {
      for (int i = 0; i < nb_points; i++)
        Un_demi_tmp(i) = Un_demi(i)*coef_col(i);
      
      MltMatrixD1var(1.0, c_half, Un_demi_tmp, Vn_tmp);
      
      for (int i = 0; i < nb_points; i++)
        Vn_tmp(i) *= alpha*coef_row(i);
      
      AddMatrixD1(1.0, Vn_tmp, coef_col, mat, m, n, incx, incy);
    }
    

    //! computes A = A + alpha * Dg(Un, Un)
    /*!
      where g(Un) = D1^v Un . D1 Un
      g is a discretization of the operator u_x (c u_x + c_x/2 u)
      Dg is the differential matrix of g with respect to Un
     */
    template<class MatrixSparse>
    void AddDifferentialUxUx(Real_wp alpha, const VectReal_wp& coef_row,
                             const VectReal_wp& coef_col, const VectReal_wp& c_half,
                             const VectReal_wp& Un, MatrixSparse& mat, int m = 0, int n = 0,
                             int incx = 1, int incy = 1)
    {
      for (int i = 0; i < nb_points; i++)
        Un_tmp(i) = Un(i)*coef_col(i);
      
      MltMatrixD1var(1.0, c_half, Un_tmp, Vn_tmp);
      MltMatrixD1(1.0, Un_tmp, Prod_tmp);
      
      for (int i = 0; i < nb_points; i++)
        {
          Vn_tmp(i) *= coef_row(i);
          Prod_tmp(i) *= coef_row(i);
        }
      
      AddMatrixD1(alpha, Vn_tmp, coef_col, mat, m, n, incx, incy);
      AddMatrixD1var(alpha, Prod_tmp, coef_col, c_half, mat, m, n, incx, incy);
    }

    
    //! initialisation of time scheme 
    /*!
      This method computes the initial condition (by calling var_eq.ComputeInitialCondition)
      and initializes the asked time scheme with this initial condition
    */
    void InitTimeScheme(int type_scheme_, int order_, const Real_wp& dt_, const Real_wp& Tf)
    {      
      time_scheme = type_scheme_;
      order_time_scheme = order_;
      int Nvol = nb_points;
      int nodl = Nvol;
      VectReal_wp Y0(nodl), x;
      GetSubdivisionSpace(x);
      Y0.Fill(0);
      this->var_eq.ComputeInitialCondition(this->x0, this->xN, x, Y0);
    
      dt = dt_;
      int N = nodl;
      
      switch (time_scheme)
	{
	case TimeSchemeEnum::LEAP_FROG :
	case TimeSchemeEnum::LEAP_FROG_EXACT :
	  {
	    // basic RK-scheme to get U^1 from U^0
	    Real_wp t = 0.0;
	    Ur_time.Reallocate(nodl);
	    Ur_next.Reallocate(nodl);
	    Ur_prev.Reallocate(nodl);
	    EvaluateFunction(t, Y0, Ur_time);
	    Copy(Ur_time, Ur_next);
	    Copy(Y0, Ur_prev);
	    Add(0.5*dt, Ur_time, Ur_prev);
	    
	    EvaluateFunction(t+0.5*dt, Ur_prev, Ur_time);
	    Add(Real_wp(2.0), Ur_time, Ur_next);
	    Copy(Y0, Ur_prev);
	    Add(0.5*dt, Ur_time, Ur_prev);
	    
	    EvaluateFunction(t+0.5*dt, Ur_prev, Ur_time);
	    Add(Real_wp(2.0), Ur_time, Ur_next);
	    Copy(Y0, Ur_prev);
	    Add(dt, Ur_time, Ur_prev);
	    
	    EvaluateFunction(t+dt, Ur_prev, Ur_time);
	    Add(Real_wp(1.0), Ur_time, Ur_next);
	    
	    for (int i = 0; i < nodl; i++)
	      {
		Ur_prev(i) = Y0(i);
		Ur_time(i) = Y0(i) + dt/6.0*Ur_next(i);
	      }
	  }
	  break;	    	  
	case TimeSchemeEnum::RUNGE_KUTTA :
	  {
	    RK_scheme.SetOrder(order_time_scheme);
	    RK_scheme.SetInitialCondition(0.0, dt, Y0, this->GetLeafClass());
	  }
	  break;
	case TimeSchemeEnum::ADAMS_BASHFORTH :
	  {
	    AB_scheme.SetOrder(order_time_scheme);
	    AB_scheme.SetInitialCondition(0.0, dt, Y0, this->GetLeafClass());
	  }
	  break;
	case TimeSchemeEnum::TALEZER :
	  {
	    talezer_scheme.SetOrder(order_time_scheme, dt);
	    talezer_scheme.SetInitialCondition(0.0, dt, Y0, this->GetLeafClass());
	  }
	  break;
	case TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE :
	  {
	    LowRK_scheme.SetOrder(order_time_scheme);
	    LowRK_scheme.SetInitialCondition(0.0, dt, Y0, this->GetLeafClass());
	  }
	  break;
	case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
	  {
            Un_demi.Reallocate(N); Un_demi.Fill(0);
	    rhs.Reallocate(N); rhs.Fill(0);
            Un.Reallocate(N); Un.Fill(0);
	    AhUn.Reallocate(order_time_scheme/2*N); AhUn.Fill(0);
	    
	    GaussRK_scheme.SetOrder(order_time_scheme);
	    GaussRK_scheme.SetInitialCondition(0.0, dt, Y0, this->GetLeafClass());
	  }
	  break;
	case TimeSchemeEnum::PREDICTOR :
	case TimeSchemeEnum::CRANK_NICOLSON :
 	  {
	    Un_demi.Reallocate(N); Un_demi.Fill(0);
	    
	    Un_next.Reallocate(N); Un_next.Fill(0);
	    rhs.Reallocate(N); rhs.Fill(0);
	    
	    Un.Reallocate(N); Un.Fill(0);
	    AhUn.Reallocate(N); AhUn.Fill(0);
	    
	    Copy(Y0, Un);
	    
	    if (time_scheme == TimeSchemeEnum::PREDICTOR)
	      {
		// for the semi-implicit time scheme, we compute U^n+1/2
		// with a step of Newton's method
		Copy(Un, Un_demi);
		
		Vector<Real_wp> U0 = Un;
		dt *= 0.5;
		this->GetLeafClass().AdvanceNewton(0.0, 0);
		Copy(Un, Un_demi);
		Copy(U0, Un);
		dt *= 2.0;
	      }
	    
	    nb_eval_diff = 0;
	    compute_differential = true;
	  }
	  break;
	}
    }


    //! Performs a step of Crank-Nicolson scheme
    void AdvanceNewton(const Real_wp& t, int iter_num)
    {
      // non-linear system to solve => Newton algorithm
      Real_wp test = 1e30, test_prec = 2e30;
      Copy(Un, Un_next);
      int nb_iter = 0;
      Real_wp threshold = 100.0*epsilon_machine, norme_init = 1.0;      
      // on continue tant que || g(Un_next) || > threshold
      // et tant que la suite des residus est bien decroissante
      while ((test > threshold)&&(test < 0.9*test_prec))
        {
          test_prec = test;
          // calcul fonctionnelle g = M (z - z^n) - dt f( (z+z^n)/2)
          // M : matrice de masse
          for (int i = 0; i < nb_points; i++)
            {
              Un_demi(i) = 0.5*(Un(i) + Un_next(i));
              Un_tmp(i) = Un_next(i) - Un(i);
            }
          
          this->GetLeafClass().ApplyMass(Un_tmp, rhs);          
          this->GetLeafClass().AddFunction(-dt, Un_demi, rhs, false);
          
          if (nb_iter == 0)
            {
              norme_init = Norm2(rhs);
              test = 1.0;
            }
          else
            test = Norm2(rhs)/norme_init;

          cout << "residu at iteration " << nb_iter << " = " << test << endl;
          
          if (compute_differential)
            {
              // computation of the differential : M - dt/2 DF  and factorisation
              this->GetLeafClass().ComputeAndFactoriseJacobian(-0.5*dt, Un_demi);
              
              compute_differential = false;
              nb_eval_diff++;
            }
          
          this->GetLeafClass().SolveJacobian(rhs);
          
          // Newton iterate : X^n+1 = X^n - DG^-1 G(X^n)
          for (int i = 0; i < nb_points; i++)
            Un_next(i) -= rhs(i);
          
          nb_iter++;
        }
      
      // if the number of iterations is large, we consider that the 
      // jacobian should be recomputed
      if (nb_iter > 10)
        compute_differential = true;
      
      // overwriting previous iterate with the computed iterate
      Copy(Un_next, Un);
    }

    
    //! computation of U^n+1 from U^n
    void Advance(const Real_wp& t, int n)
    {
      switch (time_scheme)
	{
        case TimeSchemeEnum::CRANK_NICOLSON:
          AdvanceNewton(t, n);
          break;
        case TimeSchemeEnum::PREDICTOR :
          this->GetLeafClass().AdvancePredictorScheme(t, n);
          break;
	case TimeSchemeEnum::LEAP_FROG :
	case TimeSchemeEnum::LEAP_FROG_EXACT :
	  {
	    EvaluateFunction(t, Ur_time, Ur_next);
	    Real_wp coef = 2.0*dt;
	    for (int i = 0; i < Ur_next.GetM(); i++)
	      Ur_next(i) = Ur_prev(i) + coef*Ur_next(i);
	    
	    Copy(Ur_time, Ur_prev);
	    Copy(Ur_next, Ur_time);
	  }
	  break;
	case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
	  GaussRK_scheme.Advance(t, n, *this);
	  break;
	case TimeSchemeEnum::RUNGE_KUTTA :
	  RK_scheme.Advance(t, n, *this);
	  break;
	case TimeSchemeEnum::ADAMS_BASHFORTH :
	  AB_scheme.Advance(t, n, *this);
	  break;
	case TimeSchemeEnum::TALEZER :
	  talezer_scheme.Advance(t, n, *this);
	  break;
	case TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE :
	  LowRK_scheme.Advance(t, n, *this);
	  break;
	}
    }
    
    
    //! by default, mass matrix M = Identity    
    void ApplyMass(const VectReal_wp& U, VectReal_wp& V)
    {
      Copy(U, V);
    }
    
    
    //! evaluation of g(t, X) where dX/dt = g(t, X) is the evolution system
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                                    const VectReal_wp& Xn, VectReal_wp& ProdXn,
                                    bool invert_mass = true, bool source = true)
    {
      this->GetLeafClass().GetFunction(Xn, ProdXn); 
    }
    
    
    //! evaluation of g(t, X) where dX/dt = g(t, X) is the evolution system
    void EvaluateFunction(const Real_wp& tn, const VectReal_wp& Xn, VectReal_wp& ProdXn,
			  bool invert_mass = true, bool source = true)
    {
      this->EvaluateDerivativeFunction(tn, 0, Xn, ProdXn);
    }

    
    //! factorisation of M - dt a_ij df/dx_j(y^n + dt \sum a_ij k_j)
    void FactorizeOperatorDhPlusAijKh(const Matrix<Real_wp>& A, const Real_wp& dt,
                                      const VectReal_wp& Y, const VectReal_wp& ki)
    {
      this->GetLeafClass().InitJacobian(A.GetM());
      
      for (int i = 0; i < A.GetM(); i++)
        {
          for (int k = 0; k < Y.GetM(); k++)
            {
              Un_demi(k) = Y(k);
              Un(k) = ki(k*A.GetM() + i);
              for (int j = 0; j < A.GetM(); j++)
                Un_demi(k) += dt*A(i, j)*ki(k*A.GetM() + j);
            }
          
          this->GetLeafClass().AddMass(i, i, A.GetM(), A.GetM());
          for (int j = 0; j < A.GetM(); j++)
            this->GetLeafClass().AddStiffJacobian(-dt*A(i, j), Un_demi,
                                                  i, j, A.GetM(), A.GetM());          
        }

      /*VectReal_wp kx(ki), fplus(ki), fminus(ki);
      Real_wp h = 1e-6;
      for (int jcol = 0; jcol < ki.GetM(); jcol++)
        {
          kx(jcol) += h;
          
          for (int i = 0; i < A.GetM(); i++)
            {
              for (int k = 0; k < Y.GetM(); k++)
                {
                  Un_demi(k) = Y(k);
                  Un(k) = kx(k*A.GetM() + i);
                }
              
              for (int j = 0; j < A.GetM(); j++)
                for (int k = 0; k < Y.GetM(); k++)
                  Un_demi(k) += dt*A(i, j)*kx(k*A.GetM() + j);
              
              this->GetLeafClass().ApplyMass(Un, rhs);
              this->GetLeafClass().AddFunction(-1.0, Un_demi, rhs, false);

              for (int k = 0; k < Y.GetM(); k++)
                fplus(k*A.GetM() + i) = rhs(k);
            }

          kx(jcol) -= 2.0*h;

          for (int i = 0; i < A.GetM(); i++)
            {
              for (int k = 0; k < Y.GetM(); k++)
                {
                  Un_demi(k) = Y(k);
                  Un(k) = kx(k*A.GetM() + i);
                }
              
              for (int j = 0; j < A.GetM(); j++)
                for (int k = 0; k < Y.GetM(); k++)
                  Un_demi(k) += dt*A(i, j)*kx(k*A.GetM() + j);
              
              this->GetLeafClass().ApplyMass(Un, rhs);
              this->GetLeafClass().AddFunction(-1.0, Un_demi, rhs, false);

              for (int k = 0; k < Y.GetM(); k++)
                fminus(k*A.GetM() + i) = rhs(k);
            }
          
          for (int irow = 0; irow < ki.GetM(); irow++)
            {
              Real_wp df = (fplus(irow) - fminus(irow)) / (2.0*h);
              if (abs(this->GetLeafClass().mat_sparse22(irow, jcol) - df) > 1e-6)
                {
                  DISP(irow); DISP(jcol); DISP(df);
                  DISP(this->GetLeafClass().mat_sparse(irow, jcol));
                  abort();
                }
            }

          kx(jcol) += h;
        }
      */
      this->GetLeafClass().FactoriseJacobian(A.GetM());
    }
    
    
    //! resolution of M ki - f(y^n + dt \sum a_ij k_j) = 0
    void SolveOperatorDhPlusAijKh(const VectReal_wp& tsub,
                                  const Matrix<Real_wp>& A, const Real_wp& dt,
                                  const VectReal_wp& Y, VectReal_wp& ki)
    {
      // using Newton's method to solve non-linear system
      Real_wp test = 1e30, test_prec = 2e30;
      int nb_iter = 0;
      Real_wp threshold = 100.0*epsilon_machine, norme_init(1);
      ki.Fill(0);
      // on continue tant que || g(Un_next) || > threshold
      // et tant que la suite des residus est bien decroissante
      while ((test > threshold)&&(test < 0.9*test_prec))
        {
          // calcul fonctionelle
          for (int i = 0; i < A.GetM(); i++)
            {
              for (int k = 0; k < Y.GetM(); k++)
                {
                  Un_demi(k) = Y(k);
                  Un(k) = ki(k*A.GetM() + i);
                }

              for (int j = 0; j < A.GetM(); j++)
                for (int k = 0; k < Y.GetM(); k++)
                  Un_demi(k) += dt*A(i, j)*ki(k*A.GetM() + j);
              
              this->GetLeafClass().ApplyMass(Un, rhs);
              this->GetLeafClass().AddFunction(-1.0, Un_demi, rhs, false);
              
              for (int k = 0; k < Y.GetM(); k++)
                AhUn(k*A.GetM() + i) = rhs(k);
            }
          
          test_prec = test;
          if (nb_iter == 0)
            {
              norme_init = Norm2(AhUn);
              test = 1.0;
            }
          else
            test = Norm2(AhUn)/norme_init;
          
          cout << "residu at iteration " << nb_iter << " = " << test << endl;
          //cout << "ratio = " << test_prec/test << endl;
          if (compute_differential)
            {
              FactorizeOperatorDhPlusAijKh(A, dt, Y, ki);
              compute_differential = false;
              nb_eval_diff++;
            }
          
          this->GetLeafClass().SolveGaussJacobian(AhUn, A.GetM());
          
          // Newton iterate : X^n+1 = X^n - DG^-1 G(X^n)
          for (int i = 0; i < ki.GetM(); i++)
            ki(i) -= AhUn(i);
          
          nb_iter++;
        }
      
      // if the number of iterations is large, we consider that the 
      // jacobian should be recomputed
      if (nb_iter > 10)
        compute_differential = true;      
    }
    
    
    //! returns L2 norm of the solution
    Real_wp GetEnergy()
    {
      return GetEnergy(this->GetIterate());
    }
    
    
    //! returns L2 norm of Zn
    Real_wp GetEnergy(const VectReal_wp& Zn)
    {
      Real_wp energy = 0;
      this->GetLeafClass().ApplyMass(Zn, Un_tmp);
      energy = DotProd(Zn, Un_tmp);
      
      return sqrt(energy);
    }

    void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp&, Real_wp alpha = 1.0)
    {
    }
    
  };
  
}

#define MONTJOIE_FILE_FINITE_DIFFERENCE_1D_CXX
#endif
