// Copyright (C) 2003-2009 Marc Duruflé
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#ifndef SELDON_FILE_ITERATIVE_GMRES_CXX

namespace Seldon
{

  //! Solves a linear system by using Generalized Minimum Residual (GMRES)
  /*!
    Solves the unsymmetric linear system Ax = b using restarted GMRES.

    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    See: Y. Saad and M. Schulter. GMRES: A generalized minimum residual
    algorithm for solving nonsysmmetric linear systems, SIAM
    J. Sci. Statist. Comp.  7(1986), pp, 856-869

    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] outer Iteration parameters
  */
#ifdef SELDON_WITH_VIRTUAL
  template<class T, class Vector1>
  int Gmres(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	    Preconditioner_Base<T>& M,
	    Iteration<typename ClassComplexType<T>::Treal>& outer)
#else
  template <class Titer, class MatrixSparse, class Vector1, class Preconditioner>
  int Gmres(const MatrixSparse& A, Vector1& x, const Vector1& b,
	    Preconditioner& M, Iteration<Titer> & outer)
#endif
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    Complexe zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);

    int m = outer.GetRestart();
    // V is the array of orthogonal basis contructed
    // from the Krylov subspace (v0,A*v0,A^2*v0,...,A^m*v0)
    std::vector<Vector1> V(m+1, b);

    // Upper triangular hessenberg matrix
    // we don't store the sub-diagonal
    // we apply rotations to eliminate this sub-diagonal
    Matrix<Complexe, General, RowUpTriang> H(m+1,m+1);
    H.Fill(zero);

    // s is the vector of residual norm for each inner iteration
    // w is used in the Arnoldi algorithm
    // u is a temporary vector which contains the product A*v(i)
    // r is the residual
    Vector1 w(b), r(b), u(b);
    Vector<Complexe> s(m+1);
    s.Fill(zero); w.Fill(zero); r.Fill(zero); u.Fill(zero);

    for (int i = 0; i < m+1; i++)
      V[i].Fill(zero);

    typedef typename ClassComplexType<Complexe>::Treal Treal;
    Vector<Complexe> rotations_sin(m+1);
    rotations_sin.Fill(zero);
    Vector<Treal> rotations_cos(m+1);
    rotations_cos.Fill(Treal(0));

    // we compute residual
    Copy(b, w);
    if (!outer.IsInitGuess_Null())
      outer.MltAdd(-one, A, x, one, w);
    else
      x.Fill(zero);

    // preconditioning
    M.Solve(A, w, r);
    Treal beta = Norm2(r);

    // we initialize outer
    int success_init = outer.Init(r);
    if (success_init != 0)
      return outer.ErrorCode();

    // the coefficient H(m+1,m)
    Complexe hi_ip1;

    outer.SetNumberIteration(0);
    // Loop until the stopping criteria are reached
    while (! outer.Finished(beta))
      {
	// we normalize V(0) and we init s
	Copy(r, V[0]);
	Mlt(one/beta, V[0]);
	s.Fill(zero);
	SetComplexReal(beta, s(0));

	int i = 0, k;

	// we initialize the iter iteration
	// m is the maximum number of inner iterations
	Iteration<Treal> inner(outer);
	inner.SetNumberIteration(outer.GetNumberIteration());
	inner.SetMaxNumberIteration(outer.GetNumberIteration()+m);
        H.Reallocate(m+1, m+1);
        H.Fill(zero);
        
	do
	  {
	    // product matrix vector u=A*V(i)
	    outer.Mlt(A, V[i], u);

	    // preconditioning
	    M.Solve(A, u, w);

	    // Arnoldi algorithm
	    for (k = 0; k <= i; k++)
	      {
		// h_{k,i} = \bar{v(k)} w
		H.Val(k, i) = DotProdConj(V[k], w);
		Add(-H(k,i), V[k], w);
	      }

	    // we compute h(i+1,i)
	    SetComplexReal(Norm2(w), hi_ip1);
	    Copy(w, V[i+1]);

	    // we normalize V(i+1)
	    if (hi_ip1 != zero)
	      Mlt(one/hi_ip1, V[i+1]);

	    // we apply precedent generated rotations
	    // to the last column we computed.
	    for (k = 0; k < i; k++)
	      ApplyRot(H.Val(k,i), H.Val(k+1,i),
		       rotations_cos(k), rotations_sin(k));

	    // we generate a new rotation Omega=[c s;-conj(s) c] in order to
	    // cancel h(i+1,i) and we store this rotation
	    if (hi_ip1 != zero)
	      {
		GenRot(H.Val(i,i), hi_ip1,
		       rotations_cos(i), rotations_sin(i));
	        // After this call we must have hi_ip1=0
		// GenRot must modify the entries H(i,i) hi_ip1
		// we apply the rotation to the right hand side s
		ApplyRot(s(i), s(i+1), rotations_cos(i), rotations_sin(i));
	      }

	    ++inner, ++outer, ++i;

	  } while (! inner.Finished(abs(s(i))));

	// Now we solve the triangular system H
	H.Resize(i, i); s.Resize(i);
	Solve(H, s); s.Resize(m+1);

	// new iterate x = x + sum_0^{i-1} s(k)*V(k)
	for (k = 0; k < i; k++)
	  Add(s(k), V[k], x);

	// we compute the new residual
	Copy(b, w);
	outer.MltAdd(-one, A, x, one, w);
	M.Solve(A, w, r);

	// residual norm
	beta = Norm2(r);
      }

    return outer.ErrorCode();

  }


  //! Solves a linear system by using DQGMRES
  /*!
    Solves the unsymmetric linear system Ax = b using DQGMRES.

    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] outer Iteration parameters
  */
#ifdef SELDON_WITH_VIRTUAL
  template<class T, class Vector1>
  int DQGmres(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
              Preconditioner_Base<T>& M,
              Iteration<typename ClassComplexType<T>::Treal>& outer)
#else
  template <class Titer, class MatrixSparse, class Vector1, class Preconditioner>
  int DQGmres(const MatrixSparse& A, Vector1& x, const Vector1& b,
              Preconditioner& M, Iteration<Titer> & outer)
#endif
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    Complexe zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);

    int k = outer.GetRestart();
    
    // V is the array of orthogonal basis contructed
    // from the Krylov subspace (v0,A*v0,A^2*v0,...,A^m*v0)
    std::vector<Vector1> V(k+1, b);
    std::vector<Vector1> P(k, b);
    
    // last column of upper triangular hessenberg matrix
    Vector<Complexe> Hcol(k+2);
    Hcol.Fill(zero);

    // s is the vector of residual norm for each inner iteration
    // w is used in the Arnoldi algorithm
    // u is a temporary vector which contains the product A*v(i)
    // r is the residual
    Vector1 w(b), r(b), u(b), pj(b);
    Vector<Complexe> sm(2);
    sm.Fill(zero); w.Fill(zero); r.Fill(zero); u.Fill(zero);

    for (int i = 0; i <= k; i++)
      V[i].Fill(zero);

    // rotations stored in c and s 
    typedef typename ClassComplexType<Complexe>::Treal Treal;
    Vector<Complexe> s(k+1); s.Fill(zero);
    Vector<Complexe> c(k+1); c.Fill(zero);

    // we compute residual
    Copy(b, w);
    if (!outer.IsInitGuess_Null())
      outer.MltAdd(-one, A, x, one, w);
    else
      x.Fill(zero);

    // preconditioning
    M.Solve(A, w, r);
    Treal beta = Norm2(r);

    // we initialize outer
    int success_init = outer.Init(r);
    if (success_init != 0)
      return outer.ErrorCode();
    
    outer.SetNumberIteration(0);

    // we normalize V(0) and we init s
    Copy(r, V[k-1]);
    Mlt(one/beta, V[k-1]);
    sm(0) = beta;
    
    // Loop until the stopping criteria are reached
    int j = 0;
    while (! outer.Finished(beta))
      {
        // product matrix vector u = A*V(i)
        outer.Mlt(A, V[k-1], u);
        
        // preconditioning
        M.Solve(A, u, w);

        // we compute the new column of H
        Hcol(0) = 0;
        for (int i = max(0, j-k+1); i <= j; i++)
          {
            int im = i-j+k-1;
            Hcol(im+1) = DotProdConj(V[im], w);
            Add(-Hcol(im+1), V[im], w);
          }

        // we compute h(i+1,i)
        SetComplexReal(Norm2(w), Hcol(k+1));
        Copy(w, V[k]);

        // we normalize the new column of V
        if (Hcol(k+1) != zero)
          Mlt(one/Hcol(k+1), V[k]);

        // we apply precedent generated rotations
        // to the last column we computed.
        for (int km = max(0, j-k); km < j; km++)
          {
            int im = km -j+k-1;
            Complexe L1 = Hcol(im+1), L2 = Hcol(im+2);
            Hcol(im+1) = conjugate(c(im+1))*L1 - conjugate(s(im+1))*L2;
            Hcol(im+2) = s(im+1)*L1 + c(im+1)*L2;            
          }

        // we compute the new rotation
        Treal rm = sqrt(absSquare(Hcol(k)) + absSquare(Hcol(k+1)));
        c(k) = Hcol(k)/rm; s(k) = -Hcol(k+1)/rm;
        Hcol(k) = rm; Hcol(k+1) = 0;

        // we compute the new column of Pm
        pj = V[k-1];
        for (int i = max(0, j-k); i < j; i++)
          {
            int im = i-j+k-1;
            Add(-Hcol(im+1), P[im+1], pj);
          }

        Mlt(one/Hcol(k), pj);
                
        // we apply new rotation on the right hand side
        sm(1) = s(k)*sm(0);
        sm(0) = conjugate(c(k))*sm(0);
        beta = abs(sm(1));
        
        // solution is updated
        Add(sm(0), pj, x);

        // next iteration
        j++;
        sm(0) = sm(1); sm(1) = 0;
        for (int i = 0; i < k; i++) 
          {
            V[i] = V[i+1];
            c(i) = c(i+1);
            s(i) = s(i+1);
            if (i < k-1)
              P[i] = P[i+1];
            else
              P[i] = pj;
          }

        ++outer;
      }
    
    return outer.ErrorCode();
  }

} // end namespace

#define SELDON_FILE_ITERATIVE_GMRES_CXX
#endif
