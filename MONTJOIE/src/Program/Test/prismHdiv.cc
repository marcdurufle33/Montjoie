#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_HDIV
#define MONTJOIE_WITH_NODAL_H1

#include "FiniteElement/MontjoieFiniteElement.hxx"

using namespace Montjoie;

namespace Montjoie
{
  //! edge element for hexahedron
  class WedgeHdivFirstFamilyTest : public WedgeReference<3>
  {
    public :
    // we store basis functions in P
    Vector< TinyVector< MultivariatePolynomial< Real_wp >, 3 > > P;
    // and divergence of basis functions
    Vector< MultivariatePolynomial< Real_wp > > Pdiv;
    Matrix<R3> Value_Phi;
    Matrix<Real_wp> Div_Phi;
    Matrix<Real_wp> mass_matrix_div; VectReal_wp tau_div;
    int type_basis;
    int num_fct_removed;
    
    enum {CLASSICAL, MONOMIAL, OPTIMAL};
    
    WedgeHdivFirstFamilyTest() : WedgeReference<3>()
    {
      this->elt_geom.dof_equal_nodal = false;
      type_basis = MONOMIAL;
      num_fct_removed = -1;
    }
    
    // construction of basis functions
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1)
    {
      WedgeReference<3>::ConstructFiniteElement(r+2);
      
      MultivariatePolynomial< Real_wp > x, y, z, one, zero, Q, Qx, Qy, Qz;
      x.SetOrder(3, 1); x(1, 0, 0) = 1;
      y.SetOrder(3, 1); y(0, 1, 0) = 1;
      z.SetOrder(3, 1); z(0, 0, 1) = 1;
      one.SetOrder(3, 0); one(0,0,0) = 1;
      zero.SetOrder(3, 0); zero(0,0,0) = 0;
      
      // nombre maximal de fonctions
      nb_dof_loc = 3*(order+3)*(order+1)*(order+1);
      P.Reallocate(nb_dof_loc);
      nb_dof_loc = 0;
      
      // Wr-1^3
      /*
      for (int i = 0; i <= r; i++)
	for (int j = 0; j <= r-i; j++)
	  for (int k = 0; k <= r; k++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
	      P(nb_dof_loc++).Init(Q, zero, zero);
	      P(nb_dof_loc++).Init(zero, Q, zero);	      
	    }
      
      for (int i = 0; i <= r; i++)
	for (int k = 0; k <= r; k++)
	  {
	    Q = Pow(x, i)*Pow(y, r-i)*Pow(z, k);
	    P(nb_dof_loc++).Init(x*Q, y*Q, zero);
	  }
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1-i; j++)
	  for (int k = 0; k <= r+2; k++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
	      P(nb_dof_loc++).Init(zero, zero, Q);
	    }
      */
      /*
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1-i; j++)
	  for (int k = 0; k <= r-1; k++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
	      P(nb_dof_loc++).Init(Q, zero, zero);
	      P(nb_dof_loc++).Init(zero, Q, zero);
	      P(nb_dof_loc++).Init(zero, zero, Q);
	    }
      
      for (int i = 0; i <= r-1; i++)
	for (int k = 0; k <= r-1; k++)
	  {
	    Q = Pow(x, i)*Pow(y, r-1-i)*Pow(z, k);      
	    P(nb_dof_loc++).Init(y*Q, zero, zero);
	    P(nb_dof_loc++).Init(zero, x*Q, zero);
	  }
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1-i; j++)
	  {		
	    Q = Pow(x, i)*Pow(y, j)*Pow(z, r-1);
	    P(nb_dof_loc++).Init(z*Q, zero, zero);
	    P(nb_dof_loc++).Init(zero, z*Q, zero);
	  }
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1-i; j++)
	  for (int k = 0; k <= r-1; k++)
	    if ( (i+j == r-1) || (k==r-1))
	      if ( (j == 0) || (k==r-1))
		{
		  Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
		  P(nb_dof_loc++).Init(x*Q, zero, -z*Q);
		}
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1-i; j++)
	  for (int k = 0; k <= r-1; k++)
	    if ( (i+j == r-1) || (k==r-1))
	      if ( (i == 0) || (k==r-1))
		if ( (i == r-1) || (j == r-1))
		  {
		    Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
		    P(nb_dof_loc++).Init(zero, y*Q, -z*Q);
		  }
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1-i; j++)
	  for (int k = 0; k <= r-1; k++)
	    if ( (i+j == r-1) || (k == r-1))
	      if ((i == r-1) || (k == r-1))
		if ((j == r-1) || (k == r-1))
		  {
		    Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
		    P(nb_dof_loc++).Init(-x*z*Q, -y*z*Q, z*z*Q);
		  }
      */
      
      /*if (r == 1)
	{
	  P(nb_dof_loc++).Init(one, zero, zero);
	  P(nb_dof_loc++).Init(zero, one, zero);
	  P(nb_dof_loc++).Init(zero, zero, one);
	  
	  P(nb_dof_loc++).Init(y, zero, zero);
	  P(nb_dof_loc++).Init(z, zero, zero);
	  P(nb_dof_loc++).Init(zero, x, zero);
	  P(nb_dof_loc++).Init(zero, z, zero);
	  
	  P(nb_dof_loc++).Init(x, zero, -z);
	  P(nb_dof_loc++).Init(zero, y, -z);
	  P(nb_dof_loc++).Init(-x*z, -y*z, z*z);
	}
      else
	{
      */
	  for (int i = 0; i <= r-1; i++)
	    for (int j = 0; j <= r-1-i; j++)
	      for (int k = 0; k <= r; k++)
		{
		  Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
		  P(nb_dof_loc++).Init(Q, zero, zero);
		  P(nb_dof_loc++).Init(zero, Q, zero);
		}
	  
	  for (int i = 0; i <= r; i++)
	    for (int k = 0; k <= r-1; k++)
	      {
		Q = Pow(x, i)*Pow(y, r-i)*Pow(z, k);
		P(nb_dof_loc++).Init(Q, zero, zero);
		P(nb_dof_loc++).Init(zero, Q, zero);
	      }

	  for (int i = 0; i <= r-2; i++)
	    for (int j = 0; j <= r-2-i; j++)
	      for (int k = 0; k <= r+1; k++)
		{
		  Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
		  P(nb_dof_loc++).Init(zero, zero, Q);
		}
	  
	  for (int i = 0; i <= r-1; i++)
	    for (int k = 0; k <= r; k++)
	      {
		Q = Pow(x, i)*Pow(y, r-1-i)*Pow(z, k);
		P(nb_dof_loc++).Init(zero, zero, Q);
	      }	  
	  
	  for (int i = 0; i <= r-1; i++)
	    {
	      Q = Pow(x, i)*Pow(y, r-1-i)*Pow(z, r);
	      P(nb_dof_loc++).Init(-x*Q, -y*Q, z*Q);
	    }	  
	  //}
      
	/*
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1-i; j++)
	  for (int k = 0; k <= r; k++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
	      P(nb_dof_loc++).Init(Q, zero, zero);
	      P(nb_dof_loc++).Init(zero, Q, zero);
	    }
	  
      for (int i = 0; i <= r; i++)
	for (int k = 0; k <= r-1; k++)
	  {
	    Q = Pow(x, i)*Pow(y, r-i)*Pow(z, k);
	    P(nb_dof_loc++).Init(Q, zero, zero);
	    P(nb_dof_loc++).Init(zero, Q, zero);
	  }
      
      for (int i = 0; i <= r; i++)
	for (int k = 0; k <= r; k++)
	  {
	    Q = Pow(x, i)*Pow(y, r-i)*Pow(z, k);
	    P(nb_dof_loc++).Init(x*Q, y*Q, zero);
	  }      
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1-i; j++)
	  for (int k = 0; k <= r+2; k++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
	      P(nb_dof_loc++).Init(zero, zero, Q);
	    }
      
      for (int i = 0; i <= r-1; i++)
	{
	  Q = Pow(x, i)*Pow(y, r-1-i)*Pow(z, r);
	  //P(nb_dof_loc++).Init(-x*Q, -y*Q, z*Q);
	  P(nb_dof_loc++).Init(x*Q, y*Q, zero);
	}	  
	*/
      // on rajoute toutes les fonctions de base necessaires pour avoir Pr-1 inclus dans div PrF 
      //AddFunctionsDivergence(r);

      if (num_fct_removed >= 0)
	{
	  // test quand on enleve une seule fonction de Pchapeau_r
	  // pour verifier que cette fonction est necessaire
	  for (int i = num_fct_removed; i < nb_dof_loc; i++)
	    P(i) = P(i+1);
	  
	  P.Resize(nb_dof_loc-1);
	}
      else
	P.Resize(nb_dof_loc);
      
      DISP(P);
      DISP(nb_dof_loc);
      
      // for projection on dofs, we use quadrature points
      VectR3 points_dof3d;
      points_dof3d.Reallocate(nb_points_quadrature_inside);
      for (int i = 0; i < nb_points_quadrature_inside; i++)
        points_dof3d(i) = this->PointsND(i);
      
      this->SetPointsDofND(points_dof3d);
      
      // computation of phi_i(\xi_k) for all i
      Value_Phi.Reallocate(nb_dof_loc, this->GetNbPointsQuadrature());
      for (int k = 0; k < this->GetNbPointsQuadrature(); k++)
        {
          VectR3 phi;
          ComputeValuesPhiRef(this->PointsND(k), phi);
          for (int i = 0; i < nb_dof_loc; i++)
            Value_Phi(i, k) = phi(i);
        }

      // computation of mass matrix
      nb_dof_loc = P.GetM();
      mass_matrix.Reallocate(nb_dof_loc, nb_dof_loc);
      mass_matrix.Fill(0);
      for (int k = 0; k < nb_points_quadrature_inside; k++)
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = i; j < nb_dof_loc; j++)
            mass_matrix(i, j) += this->WeightsND(k)*DotProd(Value_Phi(i, k), Value_Phi(j, k));
      
      
#ifndef MONTJOIE_WITH_MPFR
      // test pour verifier le rang de la matrice de masse
      Matrix<Real_wp, Symmetric, RowSymPacked> Mh = mass_matrix;
      VectReal_wp lambda; Matrix<Real_wp, General, RowMajor> V;
      GetEigenvaluesEigenvectors(Mh, lambda, V);
      DISP(lambda);
      int rank_mass = nb_dof_loc;
      for (int i = 0; i < lambda.GetM(); i++)
	if (abs(lambda(i)) < 10.0*epsilon_machine)
	  {
	    DISP(i); DISP(lambda(i));
	    Qx = zero; Qy = zero; Qz = zero;
	    for (int j = 0; j < lambda.GetM(); j++)
	      if (abs(V(j, i)) > 1e-5)
		{
		  DISP(j); DISP(V(j, i)); DISP(P(j));
		  Qx += V(j, i)*P(j)(0);
		  Qy += V(j, i)*P(j)(0);
		  Qz += V(j, i)*P(j)(0);
		}
	    
	    DISP(Qx); DISP(Qy); DISP(Qz);
	    rank_mass--;
	  }
      
      DISP(nb_dof_loc);
      DISP(rank_mass);
#endif
      
      //mass_matrix.WriteText("mat_mass.dat");
      
      // inversion de la matrice de masse
      mass_matrix_chol = mass_matrix;
      GetCholesky(mass_matrix_chol);
      //DISP(3*(r+3)*(r+1)*(r+1));
      //DISP(nb_dof_loc);

      elt_geom.ComputeCoefficientTransformation();
      
      nb_dof_loc = P.GetM();
      // on calcule la divergence des fonctions de Er
      Pdiv.Reallocate(nb_dof_loc);
      for (int i = 0; i < nb_dof_loc; i++)
	GetDivPolynomial(P(i), Pdiv(i));
      
      // computation of phi_i(\xi_k) for all i
      Div_Phi.Reallocate(nb_dof_loc, this->GetNbPointsQuadrature());
      for (int k = 0; k < this->GetNbPointsQuadrature(); k++)
        {
          VectReal_wp phi;
          ComputeDivPhiRef(this->PointsND(k), phi);
          for (int i = 0; i < nb_dof_loc; i++)
            Div_Phi(i, k) = phi(i);
        }
      
      ConstructMassMatrixDiv();
      //ConstructDivergenceSpace(r);
      
    }

    void GetDhatR(int r, Vector<MultivariatePolynomial<Real_wp> >& Wr, bool diff_wr = false)
    {
    }


    void GetDhatR(int r, Vector<TinyVector<int, 3> >& Wr, bool diff_wr = false)
    {
    }
		  
    void AddFunctionsDivergence(int r)
    {
      MultivariatePolynomial< Real_wp > x, y, z, zero, Q;
      x.SetOrder(3, 1); x(1, 0, 0) = 1;
      y.SetOrder(3, 1); y(0, 1, 0) = 1;
      z.SetOrder(3, 1); z(0, 0, 1) = 1;
      zero.SetOrder(3, 0); zero(0,0,0) = 0;

      if (r == 1)
	{
	  P(nb_dof_loc++).Init(x, y, z);
	  P(nb_dof_loc++).Init(x*x, x*y, zero);
	  P(nb_dof_loc++).Init(x*y, y*y, zero);
	  P(nb_dof_loc++).Init(x*x*z, x*y*z, zero);
	  P(nb_dof_loc++).Init(x*y*z, y*y*z, zero);
	  P(nb_dof_loc++).Init(zero, zero, z*z);
	  P(nb_dof_loc++).Init(zero, zero, z*z*z);
	}
      else
	{
	  for (int i = 0; i <= r-1; i++)
	    {
	      Q = Pow(x, i)*Pow(y, r-1-i)*Pow(z, r+1);
	      P(nb_dof_loc++).Init(zero, zero, Q);
	    }
	  
	  for (int i = 0; i <= r; i++)
	    for (int k = 0; k <= r; k++)
	      {
		Q = Pow(x, i)*Pow(y, r-i)*Pow(z, k);
		P(nb_dof_loc++).Init(x*Q, y*Q, zero);
	      }
	  
	  for (int i = 0; i <= r-1; i++)
	    for (int j = 0; j <= r-1-i; j++)
	      {
		Q = Pow(x, i)*Pow(y, j)*Pow(z, r+2);
		P(nb_dof_loc++).Init(zero, zero, Q);
	      }
	}
    }
    

    void ConstructMassMatrixDiv()
    {
      // on construit la matrice \int div(phi_i) div(phi_j)
      mass_matrix_div.Reallocate(nb_dof_loc, nb_dof_loc);
      mass_matrix_div.Fill(0);
      
      for (int k = 0; k < nb_points_quadrature_inside; k++)
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = 0; j < nb_dof_loc; j++)
            mass_matrix_div(i, j) += this->WeightsND(k)*Div_Phi(i, k)*Div_Phi(j, k);
      
      GetQR(mass_matrix_div, tau_div);
    }
		  
    // projection on reference element
    void ComputeProjectionDivRef(const VectReal_wp& feval, VectReal_wp& res)
    {
      // integration against basis functions
      // f = \int p phi_i(x) dx   (where p is stored in feval)
      VectReal_wp F(nb_dof_loc);
      res.Reallocate(nb_dof_loc); res.Fill(0); F.Fill(0);
      for (int k = 0; k < nb_points_quadrature_inside; k++)
        for (int i = 0; i < nb_dof_loc; i++)
          res(i) += this->WeightsND(k)*feval(k)*Div_Phi(i, k);
      
      // then solving by mass matrix
      MltQ_FromQR(SeldonTrans, mass_matrix_div, tau_div, res);
      
      for (int j = nb_dof_loc-1; j >= 0; j--)
	{
	  if (abs(mass_matrix_div(j, j)) > 1e3*epsilon_machine)
	    {
	      for (int k = j+1; k < nb_dof_loc; k++)
		res(j) -= mass_matrix_div(j, k)*res(k);
	      
	      res(j) /= mass_matrix_div(j, j);
	    }
	}
    }

    void ConstructDivergenceSpace(int r)
    {      
    }

        
    void SubtractFunction(int n)
    {
      num_fct_removed = n;
    }
    
    
    // projection on reference element
    void ComputeProjectionDofRef(const VectR3& feval, VectReal_wp& res)
    {
      // integration against basis functions
      // f = \int p phi_i(x) dx   (where p is stored in feval)
      res.Reallocate(nb_dof_loc); res.Fill(0);
      for (int k = 0; k < nb_points_quadrature_inside; k++)
        for (int i = 0; i < nb_dof_loc; i++)
          res(i) += this->WeightsND(k)*DotProd(feval(k), Value_Phi(i, k));
      
      // then solving by mass matrix
      SolveMassMatrix(res);
      
    }
    
    
    // computation of basis functions
    void ComputeValuesPhiRef(const R3& pt_loc, VectR3& phi) const
    {
      phi.Reallocate(nb_dof_loc);
      Real_wp xt = pt_loc(0);
      Real_wp yt = pt_loc(1);
      Real_wp zt = pt_loc(2);
      R3 pt(xt, yt, zt), vec_u;
      for (int i = 0; i < nb_dof_loc; i++)
	{
	  vec_u(0) = P(i)(0).Evaluate(pt);
	  vec_u(1) = P(i)(1).Evaluate(pt);
	  vec_u(2) = P(i)(2).Evaluate(pt);
          
	  phi(i) = vec_u;
	}         
      
    }
    
    // computation of divergence of basis functions
    void ComputeDivPhiRef(const R3& pt_loc, VectReal_wp& phi) const
    {
      phi.Reallocate(nb_dof_loc);
      FillZero(phi);
      for (int i = 0; i < nb_dof_loc; i++)
	phi(i) = Pdiv(i).Evaluate(pt_loc);      
    }
    
    
    void GetOptimalSpace(int r, Vector<TinyVector<MultivariatePolynomial<Real_wp>, 3> >& Dr)
    {
      MultivariatePolynomial< Real_wp > x, y, z, zero, Q;
      x.SetOrder(3, 1); x(1, 0, 0) = 1;
      y.SetOrder(3, 1); y(0, 1, 0) = 1;
      z.SetOrder(3, 1); z(0, 0, 1) = 1;
      zero.SetOrder(3, 0); zero(0,0,0) = 0;

      int nb = 0;
      for (int i = 0; i <= r; i++)
	for (int j = 0; j <= r-i; j++)
	  for (int k = 0; k <= r; k++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
	      Dr(nb++).Init(Q, zero, zero);
	      Dr(nb++).Init(zero, Q, zero);	      
	    }
      
      for (int i = 0; i <= r; i++)
	for (int k = 0; k <= r; k++)
	  {
	    Q = Pow(x, i)*Pow(y, r-i)*Pow(z, k);
	    Dr(nb++).Init(x*Q, y*Q, zero);
	  }
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1-i; j++)
	  for (int k = 0; k <= r+2; k++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
	      Dr(nb++).Init(zero, zero, Q);
	    }
    }
    
    
    void FindMonomialBasis(int r)
    {
      Vector<TinyVector<MultivariatePolynomial<Real_wp>, 3> > Dr;
      
      int nb_fcts_dr = (r+3)*(r+1)*(3*r+2)/2;
      Dr.Reallocate(nb_fcts_dr);
      GetOptimalSpace(r, Dr);
      
      int nb_pts_quad = nb_points_quadrature_inside;
      VectR3 feval(nb_pts_quad); VectReal_wp contrib(nb_dof_loc);
      Matrix<R3> ValPhiDr(nb_dof_loc, nb_pts_quad);
      ValPhiDr.Fill(0);
      nb_fcts_dr = 0;
      cout << endl << "The following functions of Dr are not belonging to Pr" <<endl;
      //cout << endl << "The following functions of Dr are belonging to Pr" <<endl;
      for (int i = 0; i < Dr.GetM(); i++)
	{
	  for (int n = 0; n < nb_pts_quad; n++)
	    {
	      feval(n)(0) = Dr(i)(0).Evaluate(PointsND(n));
	      feval(n)(1) = Dr(i)(1).Evaluate(PointsND(n));
	      feval(n)(2) = Dr(i)(2).Evaluate(PointsND(n));
	    }
	  
	  ComputeProjectionDofRef(feval, contrib);
	  
	  bool fct_in_Pr = true;
	  Real_wp err_max = 0;
	  for (int n = 0; n < nb_pts_quad; n++)
	    {
	      R3 val_ex = feval(n);
	      R3 val_proj;
	      for (int m = 0; m < nb_dof_loc; m++)
		Add(contrib(m), Value_Phi(m, n), val_proj);
	      
	      Real_wp dist = val_proj.Distance(val_ex);
	      if (dist > 1e9*epsilon_machine)
		{
		  fct_in_Pr = false;
		  if (dist > err_max)
		    err_max = dist;
		}
	    }
		
	  if (!fct_in_Pr)
	    {
	      cout << Dr(i) << endl;
	      DISP(err_max);
	    }
	  else
	    {
	      /*cout << Dr(i) << endl;
	      for (int n = 0; n < nb_dof_loc; n++)
		if (abs(contrib(n)) > 1e9*epsilon_machine)
		  {
		    DISP(contrib(n));
		    DISP(P(n));
		  }
	      */
	      for (int n = 0; n < nb_pts_quad; n++)
		ValPhiDr(nb_fcts_dr, n) = feval(n);
	      
	      nb_fcts_dr++;
	    }
	}
      
      
      DISP(nb_fcts_dr);
      // matrice de masse du sous-espace de Dr
      Matrix<Real_wp, Symmetric, RowSymPacked> mat_mass_Dr(nb_fcts_dr, nb_fcts_dr);
      mat_mass_Dr.Fill(0);
      for (int i = 0; i < nb_fcts_dr; i++)
	for (int j = i; j < nb_fcts_dr; j++)
	  for (int k = 0; k < nb_pts_quad; k++)
	    mat_mass_Dr(i, j) += WeightsND(k)*DotProd(ValPhiDr(i, k), ValPhiDr(j, k));
      
      GetCholesky(mat_mass_Dr);
      
      // on cherche les fonctions de Pr qui ne sont pas dans le sous-espace de Dr
      contrib.Reallocate(nb_fcts_dr);
      cout << endl << "The following functions of Pr are linked " << endl;
      for (int i = 0; i < nb_dof_loc; i++)
	{
	  for (int n = 0; n < nb_pts_quad; n++)
	    feval(n) = Value_Phi(i, n);
	  
	  contrib.Fill(0);
	  for (int n = 0; n < nb_pts_quad; n++)
	    for (int k = 0; k < nb_fcts_dr; k++)
	      contrib(k) += WeightsND(n)*DotProd(ValPhiDr(k, n), feval(n));
	  
	  Seldon::SolveCholesky(SeldonNoTrans, mat_mass_Dr, contrib);
	  Seldon::SolveCholesky(SeldonTrans, mat_mass_Dr, contrib);
	  
	  bool fct_in_Dr = true;
	  Real_wp err_max = 0;
	  for (int n = 0; n < nb_pts_quad; n++)
	    {
	      R3 val_ex = feval(n);
	      R3 val_proj;
	      for (int m = 0; m < nb_fcts_dr; m++)
		Add(contrib(m), ValPhiDr(m, n), val_proj);
	      
	      Real_wp dist = val_proj.Distance(val_ex);
	      if (dist > 1e9*epsilon_machine)
		{
		  fct_in_Dr = false;
		  if (dist > err_max)
		    err_max = dist;
		}
	    }
	  
	  if (!fct_in_Dr)
	    {
	      cout << "i = " << i << " phi = " << P(i) << endl;
	      //DISP(err_max);
	    }
	}
      
      exit(0);

    }     
    
  };

}


void GetDkBasis(int r, Vector<TinyVector<MultivariatePolynomial<Real_wp>, 3> >& psi)
{
  int nb_fcts = r*(r+1)*(r+2)/2 + r*(r+1)/2;

  MultivariatePolynomial< Real_wp > x, y, z, one, zero;
  x.SetOrder(3, 1); x(1, 0, 0) = 1;
  y.SetOrder(3, 1); y(0, 1, 0) = 1;
  z.SetOrder(3, 1); z(0, 0, 1) = 1;
  one.SetOrder(3, 0); one(0, 0, 0) = 1;
  zero.SetOrder(3, 0); zero(0, 0, 0) = 0;
      
  // we put P_{r-1}^3
  psi.Reallocate(nb_fcts);
  int node = 0;
  for (int i = 0; i < r; i++)
    for (int j = 0; j < r-i; j++)
      for (int k = 0; k < r-i-j; k++)
        {
          psi(node)(0) = Pow(x, i)*Pow(y, j)*Pow(z, k);
          psi(node)(1) = zero;
          psi(node)(2) = zero;
          node++;
          
          psi(node)(0) = zero;
          psi(node)(1) = Pow(x, i)*Pow(y, j)*Pow(z, k);
          psi(node)(2) = zero;
          node++;
          
          psi(node)(0) = zero;
          psi(node)(1) = zero;
          psi(node)(2) = Pow(x, i)*Pow(y, j)*Pow(z, k);
          node++;
        }
  
  // then P_r-1 (x, y, z)
  for (int i = 0; i < r; i++)
    for (int j = 0; j < r-i; j++)
      {
        psi(node)(0) = Pow(x, i+1)*Pow(y, j)*Pow(z, r-1-i-j);
        psi(node)(1) = Pow(x, i)*Pow(y, j+1)*Pow(z, r-1-i-j);
        psi(node)(2) = Pow(x, i)*Pow(y, j)*Pow(z, r-i-j);
        node++;
      }
}


void GetPkBasis(int r, Vector<TinyVector<MultivariatePolynomial<Real_wp>, 3> >& psi)
{
  int nb_fcts = r*(r+1)*(r+2)/2;

  MultivariatePolynomial< Real_wp > x, y, z, one, zero;
  x.SetOrder(3, 1); x(1, 0, 0) = 1;
  y.SetOrder(3, 1); y(0, 1, 0) = 1;
  z.SetOrder(3, 1); z(0, 0, 1) = 1;
  one.SetOrder(3, 0); one(0, 0, 0) = 1;
  zero.SetOrder(3, 0); zero(0, 0, 0) = 0;
      
  // we put P_{r-1}^3
  psi.Reallocate(nb_fcts);
  int node = 0;
  for (int i = 0; i < r; i++)
    for (int j = 0; j < r-i; j++)
      for (int k = 0; k < r-i-j; k++)
        {
          psi(node)(0) = Pow(x, i)*Pow(y, j)*Pow(z, k);
          psi(node)(1) = zero;
          psi(node)(2) = zero;
          node++;
          
          psi(node)(0) = zero;
          psi(node)(1) = Pow(x, i)*Pow(y, j)*Pow(z, k);
          psi(node)(2) = zero;
          node++;
          
          psi(node)(0) = zero;
          psi(node)(1) = zero;
          psi(node)(2) = Pow(x, i)*Pow(y, j)*Pow(z, k);
          node++;
        }
  
}

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc != 2)
    {
      cout << "Usage : ./pyramid.x order" <<endl;
      abort();
    }
  
  int r = atoi(argv[1]); int r_test = r;
  
  WedgeHdivFirstFamilyTest hex;
  hex.ConstructFiniteElement(r);
  //hex.FindMonomialBasis(r);
  //exit(0);
  Mesh<Dimension3> mesh;
  
  mesh.ReallocateVertices(6);
  mesh.ReallocateElements(1);
  mesh.Element(0).InitWedge(0, 1, 2, 3, 4, 5, 1);
  mesh.SetGeometryOrder(r);
  
  mesh.ReorientElements();

  VectR3 s; SetPoints<Dimension3> PointsElem;
  SetMatrices<Dimension3> MatricesElem;
  int nb_pts_quad = hex.GetNbPointsQuadratureInside();
  int nb_pts_dof = hex.GetNbPointsDof();
  int nb_dof_loc = hex.GetNbDof();
  VectR3 feval(nb_pts_dof);
  VectReal_wp feval_div(nb_pts_dof);
  //nb_dof_loc--;
  //for (int nt = 0; nt <= nb_dof_loc; nt++)
  //{
  //hex.SubtractFunction(nt);
  //  hex.ConstructFiniteElement(r);
      
  VectReal_wp res(nb_dof_loc), res_div(nb_dof_loc);
  R3 vec_u, u_approche, u_exact;
  //Real_wp div_exact, div_approche;
  VectR3 Points(nb_pts_quad);
  for (int i = 0; i < nb_pts_quad; i++)
    Points(i) = hex.PointsND(i);
  
  Matrix3_3 dfjm1; 
  Real_wp jacob;
  Vector<TinyVector<MultivariatePolynomial<Real_wp>, 3> > Pk;
  GetPkBasis(r_test, Pk);
  // GetDkBasis(r_test, Pk);
  R3 pt; Real_wp h = 1.0;
  
  // first test : we check if we can generate Pk with basis functions
  //bool fct_needed = false;
  for (int num = 1; num <= Pk.GetM(); num++)
    {
      mesh.Vertex(0).Init(-0.9*h, -1.2*h, -0.4*h);
      mesh.Vertex(1).Init(1.2*h, -0.95*h, -0.2*h);
      mesh.Vertex(2).Init(-0.87*h, 1.05*h, -0.3*h);
      mesh.Vertex(3).Init(-0.95*h, -1.08*h, 0.97*h);
      mesh.Vertex(4).Init(1.13*h, -0.96*h, 1.21*h);
      mesh.Vertex(5).Init(-0.91*h, 1.02*h, 1.24*h);
	  
      // boucle sur les elements
      for (int i = 0; i < mesh.GetNbElt(); i++)
	{
	  // on calcule F_i(x), DF_i(x)
	  mesh.GetVerticesElement(i, s);
	  hex.FjElem(s, PointsElem, mesh, i);
	  hex.DFjElem(s, PointsElem, MatricesElem, mesh, i);
	  
	  // on calcule P = J DF^-1 p avec p polynome de Pk-1^3
	  // et Pdiv = J q avec q polynome de Pk-1
	  for (int k = 0; k < nb_pts_dof; k++)
	    {
	      pt = PointsElem.GetPointDof(k);
	      vec_u(0) = Pk(num-1)(0).Evaluate(pt);
	      vec_u(1) = Pk(num-1)(1).Evaluate(pt);
	      vec_u(2) = Pk(num-1)(2).Evaluate(pt);
	      jacob = Det(MatricesElem.GetPointDof(k));
	      GetInverse(MatricesElem.GetPointDof(k), dfjm1);
	      Mlt(dfjm1, vec_u, feval(k));
	      feval(k) *= jacob;
	      feval_div(k) = vec_u(0)*jacob;
	    }
	  
	  // et on projette sur les fonctions de base de l'element de reference
	  hex.ComputeProjectionDofRef(feval, res);
	  
	  // on projette vis-a-vis de la divergence
	  if (num%3 == 1)
	    hex.ComputeProjectionDivRef(feval_div, res_div);
	  
	  // on calcule la norme L2 entre P = J DF^-1 p et son projete
	  Real_wp erreur = 0, norme_L2 = 0; R3 evalP, eval_proj;
	  Real_wp divP, div_proj, erreur_div = 0;
	  for (int k = 0; k < nb_pts_quad; k++)
	    {
	      eval_proj.Zero();
	      for (int i = 0; i < nb_dof_loc; i++)
		eval_proj += res(i)*hex.Value_Phi(i, k);
	      
	      // evaluation de P
	      pt = PointsElem.GetPointQuadrature(k);
	      vec_u(0) = Pk(num-1)(0).Evaluate(pt);
	      vec_u(1) = Pk(num-1)(1).Evaluate(pt);
	      vec_u(2) = Pk(num-1)(2).Evaluate(pt);
	      jacob = Det(MatricesElem.GetPointDof(k));
	      GetInverse(MatricesElem.GetPointDof(k), dfjm1);
	      Mlt(dfjm1, vec_u, evalP);
	      evalP *= jacob;
	      
	      erreur += hex.WeightsND(k)*DotProd(evalP - eval_proj, evalP - eval_proj);
	      norme_L2 += hex.WeightsND(k)*DotProd(evalP, evalP);
	      
	      if (num%3 == 1)
		{
		  div_proj = 0;
		  for (int i = 0; i < nb_dof_loc; i++)
		    div_proj += res_div(i)*hex.Div_Phi(i, k);
		  
		  divP = vec_u(0)*jacob;
		  erreur_div += hex.WeightsND(k)*square(divP - div_proj);
		}
	    }
	  
	  erreur = sqrt(erreur);
	  norme_L2 = sqrt(norme_L2);
	  //if (erreur/norme_L2 > 1e-30)
	  //fct_needed = true;
	  
	  DISP(num); DISP(Pk(num-1)); DISP(erreur/norme_L2);	      
	  if (num%3 == 1)
	    {
	      erreur_div = sqrt(erreur_div);
	      DISP(erreur_div/norme_L2);
	    }
	}	  	  
    }
  
  //if (fct_needed)
  // cout << "Function " << nt << " is necessary" << endl;
  //else
  //cout << "Function " << nt << " is not necessary" << endl;
  //}
  
  return FinalizeMontjoie();
}
