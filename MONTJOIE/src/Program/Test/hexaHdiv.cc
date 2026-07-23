#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_HDIV
#define MONTJOIE_WITH_NODAL_H1

#include "FiniteElement/MontjoieFiniteElement.hxx"

using namespace Montjoie;

// 3-D fonction, we want to interpolate
R3 f(const R3& x)
{
  R3 u;
  u(0) = 1.0 + 2.0*x(0) -3.0*x(1) + 4.0*x(2);
  u(1) = 2 + 2.5*x(0) + 2.0*x(1) + 3.3*x(2);
  u(2) = 3  - 0.5*x(0) -1.2*x(1) + 2.8*x(2);
  
  return u; 
}

// divergence of this function (for error with divergence)
Real_wp fdiv(const R3& x)
{
  Real_wp u = 0;
  
  return u; 
}


namespace Montjoie
{
  //! edge element for hexahedron
  class HexahedronHdivFirstFamilyBis : public HexahedronReference<3>
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
    //Matrix<Real_wp> Qmass_matrix, Rmass_matrix;
    
    enum {CLASSICAL, MONOMIAL, OPTIMAL};
    
    HexahedronHdivFirstFamilyBis() : HexahedronReference<3>()
    {
      this->elt_geom.dof_equal_nodal = false;
      // type_basis = OPTIMAL;
      // type_basis = CLASSICAL;
      type_basis = MONOMIAL;
      num_fct_removed = -1;
    }
    
    // construction of basis functions
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1)
    {
      HexahedronReference<3>::ConstructFiniteElement(r+2);
      
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

      /*
      for (int i = 0; i <= r+2; i++)
	for (int j = 0; j <= r; j++)
	  for (int k = 0; k <= r; k++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
	      P(nb_dof_loc).Init(Q, zero, zero);
	      Q = Pow(x, j)*Pow(y, i)*Pow(z, k);
	      P(nb_dof_loc+1).Init(zero, Q, zero);
	      Q = Pow(x, j)*Pow(y, k)*Pow(z, i);
	      P(nb_dof_loc+2).Init(zero, zero, Q);
              nb_dof_loc += 3;
	    }
      */      
      /*
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1; j++)
	  for (int k = 0; k <= r-1; k++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
	      P(nb_dof_loc).Init(Q, zero, zero);
	      P(nb_dof_loc+1).Init(zero, Q, zero);
	      P(nb_dof_loc+2).Init(zero, zero, Q);
	      nb_dof_loc += 3;
	    }
      
      for (int j = 0; j <= r-1; j++)
	for (int k = 0; k <= r-1; k++)
	  {
	    Q = Pow(x, r)*Pow(y, j)*Pow(z, k);
	    P(nb_dof_loc).Init(Q, zero, zero);
	    P(nb_dof_loc+1).Init(zero, Q, zero);
	    P(nb_dof_loc+2).Init(zero, zero, Q);
	    nb_dof_loc += 3;
	    
	    Q = Pow(y, r)*Pow(x, j)*Pow(z, k);
	    P(nb_dof_loc).Init(Q, zero, zero);
	    P(nb_dof_loc+1).Init(zero, Q, zero);
	    P(nb_dof_loc+2).Init(zero, zero, Q);
	    nb_dof_loc += 3;
	    
	    Q = Pow(z, r)*Pow(x, j)*Pow(y, k);
	    P(nb_dof_loc).Init(Q, zero, zero);
	    P(nb_dof_loc+1).Init(zero, Q, zero);
	    P(nb_dof_loc+2).Init(zero, zero, Q);
	    nb_dof_loc += 3;
	  }
      
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1; j++)
	  P(nb_dof_loc++).Init(Pow(x, i+2)*Pow(y, j)*Pow(z, r), zero, -Pow(x, i+1)*Pow(y, j)*Pow(z, r+1));
      
      for (int i = 0; i <= r-2; i++)
	{
	  P(nb_dof_loc++).Init(Pow(x, r)*Pow(y, r)*Pow(z, i), zero, -Pow(x, r+1)*Pow(y, r)*Pow(z, i+1));
	  P(nb_dof_loc++).Init(Pow(x, i+1)*Pow(y, r)*Pow(z, r-1), zero, -Pow(x, i)*Pow(y, r)*Pow(z, r));
	  
	  P(nb_dof_loc++).Init(Pow(x, r)*Pow(y, i)*Pow(z, r), -Pow(x, r-1)*Pow(y, i+1)*Pow(z, r),zero);
	  P(nb_dof_loc++).Init(Pow(x, i+1)*Pow(y, r-1)*Pow(z, r), -Pow(x, i)*Pow(y, r)*Pow(z, r),zero);
	  
	  P(nb_dof_loc++).Init(zero,Pow(x, r)*Pow(y, r)*Pow(z, i), -Pow(x, r)*Pow(y, r-1)*Pow(z, i+1));
	  P(nb_dof_loc++).Init(zero,Pow(x, r)*Pow(y, i+1)*Pow(z, r-1), -Pow(x, r)*Pow(y, i)*Pow(z, r));
	  
	  P(nb_dof_loc++).Init(Pow(x, r+1)*Pow(y, r)*Pow(z, i), -Pow(x, r)*Pow(y, r+1)*Pow(z, i),zero);
	  P(nb_dof_loc++).Init(Pow(x, r+1)*Pow(y, i+1), -Pow(x, r)*Pow(y, i+2)*Pow(z, r),zero);
	  
	  P(nb_dof_loc++).Init(zero,Pow(x, i)*Pow(y, r+1)*Pow(z, r), -Pow(x, i)*Pow(y, r)*Pow(z, r+1));
	  P(nb_dof_loc++).Init(zero,Pow(y, r+1)*Pow(z, i+1), -Pow(y, r)*Pow(z, i+2));
	}
      
      if (r>1)
	{
	  for (int i = 0; i <= r-1; i++)
	    {
	      P(nb_dof_loc++).Init(Pow(x, r+1)*Pow(y, i)*Pow(z, r), zero, -Pow(x, r)*Pow(y, i)*Pow(z, r+1));
	      P(nb_dof_loc++).Init(Pow(x, r)*Pow(y, i-1)*Pow(z, r-1), zero, -Pow(x, r-1)*Pow(y, i)*Pow(z, r));
	    }
	  
	  P(nb_dof_loc++).Init(Pow(x, r+1)*Pow(y, r-1), -Pow(x, r)*Pow(y, r),zero);
	  P(nb_dof_loc++).Init(zero,Pow(x, r+1)*Pow(y, r-1)*Pow(y, r-1), -Pow(y, r)*Pow(z, r));
	}
      
      P(nb_dof_loc++).Init(Pow(x, r)*Pow(y, r)*Pow(z, r-1), zero, -Pow(x, r-1)*Pow(y, r)*Pow(z, r));
      P(nb_dof_loc++).Init(Pow(x, r)*Pow(y, r-1)*Pow(z, r), -Pow(x, r-1)*Pow(y, r)*Pow(z, r),zero);
      P(nb_dof_loc++).Init(zero,Pow(x, r)*Pow(y, r)*Pow(z, r-1), -Pow(x, r)*Pow(y, r-1)*Pow(z, r));
      
      P(nb_dof_loc++).Init(Pow(x, r+1)*Pow(y, r)*Pow(z, r-1), -Pow(x, r)*Pow(y, r+1)*Pow(z, r-1), zero);
      P(nb_dof_loc++).Init(zero,Pow(y, r+1)*Pow(z, r), -Pow(y, r)*Pow(z, r+1));
      
      P(nb_dof_loc++).Init(Pow(x, r+1)*Pow(y, r-1)*Pow(z, r-1), -Pow(x, r)*Pow(y, r)*Pow(z, r-1), -Pow(x, r)*Pow(y, r-1)*Pow(z, r));
      P(nb_dof_loc++).Init(-Pow(x, r)*Pow(y, r)*Pow(z, r-1), Pow(x, r-1)*Pow(y, r+1)*Pow(z, r-1),-Pow(x, r-1)*Pow(y, r)*Pow(z, r));
      P(nb_dof_loc++).Init(-Pow(x, r)*Pow(y, r-1)*Pow(z, r), -Pow(x, r-1)*Pow(y, r)*Pow(z, r), Pow(x, r-1)*Pow(y, r-1)*Pow(z, r+1));	
      */
      
      // on calcule tous les monomes de Er
      /*
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1; j++)
	  for (int k = 0; k <= r-1; k++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
	      P(nb_dof_loc).Init(Q, zero, zero);
	      P(nb_dof_loc+1).Init(zero, Q, zero);
	      P(nb_dof_loc+2).Init(zero, zero, Q);
              nb_dof_loc += 3;
	    }

      for (int j = 0; j <= r-1; j++)
	for (int k = 0; k <= r-1; k++)
	  {
	    Q = Pow(x, r)*Pow(y, j)*Pow(z, k);
	    P(nb_dof_loc).Init(Q, zero, zero);
	    P(nb_dof_loc+1).Init(zero, Q, zero);
	    P(nb_dof_loc+2).Init(zero, zero, Q);
	    nb_dof_loc += 3;
	  }      

      for (int j = 0; j <= r-1; j++)
	for (int k = 0; k <= r-1; k++)
	  {
	    Q = Pow(y, r)*Pow(x, j)*Pow(z, k);
	    P(nb_dof_loc).Init(Q, zero, zero);
	    P(nb_dof_loc+1).Init(zero, Q, zero);
	    P(nb_dof_loc+2).Init(zero, zero, Q);
	    nb_dof_loc += 3;
	  }      
      
      for (int j = 0; j <= r-1; j++)
	for (int k = 0; k <= r-1; k++)
	  {
	    Q = Pow(z, r)*Pow(x, j)*Pow(y, k);
	    P(nb_dof_loc).Init(Q, zero, zero);
	    P(nb_dof_loc+1).Init(zero, Q, zero);
	    P(nb_dof_loc+2).Init(zero, zero, Q);
	    nb_dof_loc += 3;
	  }      
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1; j++)
	  for (int k = 0; k <= r-1; k++)
	    if ( ((i>=r-1)&&(j >= r-1)) || ((k>=r-1)&&(j >= r-1)) )
	      P(nb_dof_loc++).Init(Pow(x, i+1)*Pow(y, j+1)*Pow(z, k), zero, -Pow(x, i)*Pow(y, j+1)*Pow(z, k+1));

      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1; j++)
	  for (int k = 0; k <= r-1; k++)
	    if ( ((i>=r-1)&&(k >= r-1)) || ((k>=r-1)&&(j >= r-1)) )
	      P(nb_dof_loc++).Init(Pow(x, i+1)*Pow(y, j)*Pow(z, k+1), -Pow(x, i)*Pow(y, j+1)*Pow(z, k+1), zero);	

      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1; j++)
	  for (int k = 0; k <= r-1; k++)
	    if ( ((i>=r-1)&&(j >= r-1)) || ((k>=r-1)&&(i >= r-1)) )
	      if ((i == r-1) || (j == 0) || (k == 0))
		P(nb_dof_loc++).Init(zero, Pow(x, i+1)*Pow(y, j+1)*Pow(z, k), -Pow(x, i+1)*Pow(y, j)*Pow(z, k+1));
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1; j++)
	  for (int k = 0; k <= r-1; k++)
	    if ( (i>=r-1) || (k>=r-1) || ((i+2>= r) &&(k+1>=r)) || ((i+1>= r) &&(k+2>=r)) )
	      {
		P(nb_dof_loc++).Init(Pow(x, i+2)*Pow(y, j)*Pow(z, k+1), zero, -Pow(x, i+1)*Pow(y, j)*Pow(z, k+2));
		
		if ( ((k == r-1) || (j == 0)) && ((i == r-1) || (k == r-1)) )
		  P(nb_dof_loc++).Init(Pow(x, i+2)*Pow(z, j)*Pow(y, k+1), -Pow(x, i+1)*Pow(z, j)*Pow(y, k+2), zero);
		
		if ( ((i == r-1) || (k == r-1)) && ((j==0) || (i==r-1)) && ((k == r-1) || (j == 0)))
		  P(nb_dof_loc++).Init(zero, Pow(y, i+2)*Pow(x, j)*Pow(z, k+1), -Pow(y, i+1)*Pow(x, j)*Pow(z, k+2));
	      }
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1; j++)
	  for (int k = 0; k <= r-1; k++)
	    if ( (i>=r-1) || ((i+1>= r) &&(j+1>=r)) || ((i+1>= r) &&(k+1>=r)) )
	      if ( ((j == r-1) || (k==0)) && ((j == 0) || (k==r-1)))
		{
		  //DISP(i); DISP(j); DISP(k); DISP(nb_dof_loc);
		  P(nb_dof_loc++).Init(Pow(x, i+2)*Pow(y, j)*Pow(z, k), -Pow(x, i+1)*Pow(y, j+1)*Pow(z, k),
				       -Pow(x, i+1)*Pow(y, j)*Pow(z, k+1));
		}
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1; j++)
	  for (int k = 0; k <= r-1; k++)
	    if ( (j>=r-1) || ((i+1>= r) &&(j+1>=r)) || ((j+1>= r) &&(k+1>=r)) )
	      if ( ((i == r-1) || (k==0)) && ((i == 0) || (k==r-1)))
		if ( ((i == 0) || (j==r-1) || (k==r-1)) && ((i==0) || (j==r-1) || (k==0)) )
		  {
		    P(nb_dof_loc++).Init(-Pow(x, i+1)*Pow(y, j+1)*Pow(z, k), Pow(x, i)*Pow(y, j+2)*Pow(z, k),
					 -Pow(x, i)*Pow(y, j+1)*Pow(z, k+1));
		  }
      
      for (int i = 0; i <= r-1; i++)
	for (int j = 0; j <= r-1; j++)
	  for (int k = 0; k <= r-1; k++)
	    if ( (k>=r-1) || ((i+1>= r) &&(k+1>=r)) || ((j+1>= r) &&(k+1>=r)) )
	      if ( ((i == r-1) || (j==0)) && ((i == 0) || (j==r-1)))
		if ( ((i == 0) || (j==0) || (k==r-1)) && ((j==0) || (k==r-1)) && ((i == 0) || (k==r-1)) )
		  {
		    P(nb_dof_loc++).Init(-Pow(x, i+1)*Pow(y, j)*Pow(z, k+1), -Pow(x, i)*Pow(y, j+1)*Pow(z, k+1),
					 Pow(x, i)*Pow(y, j)*Pow(z, k+2));
		  }
      */
      /*
      if (r == 1)
	{
	  P(nb_dof_loc++).Init(one, zero, zero);
	  P(nb_dof_loc++).Init(x, zero, zero);
	  P(nb_dof_loc++).Init(y, zero, zero);
	  P(nb_dof_loc++).Init(z, zero, zero);
	  
	  P(nb_dof_loc++).Init(zero, one, zero);
	  P(nb_dof_loc++).Init(zero, x, zero);
	  P(nb_dof_loc++).Init(zero, y, zero);
	  P(nb_dof_loc++).Init(zero, z, zero);
	  
	  P(nb_dof_loc++).Init(zero, zero, one);
	  P(nb_dof_loc++).Init(zero, zero, x);
	  P(nb_dof_loc++).Init(zero, zero, y);
	  P(nb_dof_loc++).Init(zero, zero, z);
	  
	  P(nb_dof_loc++).Init(x*y, zero, -y*z);
	  P(nb_dof_loc++).Init(x*z, -y*z, zero);
	  P(nb_dof_loc++).Init(zero, x*y, -x*z);
	  
	  P(nb_dof_loc++).Init(x*x, -x*y, -x*z);
	  P(nb_dof_loc++).Init(-x*y, y*y, -y*z);
	  P(nb_dof_loc++).Init(-x*z, -y*z, z*z);

	  P(nb_dof_loc++).Init(x*x*z, zero, -x*z*z);
	  P(nb_dof_loc++).Init(x*x*y, -x*y*y, zero);
	  P(nb_dof_loc++).Init(zero, y*y*z, -y*z*z);
	}
      else
	{
	  for (int i = 0; i <= r+1; i++)
	    for (int j = 0; j <= r-1; j++)
	      for (int k = 0; k <= r-1; k++)
		{
		  Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
		  P(nb_dof_loc++).Init(Q, zero, zero);
		  Q = Pow(y, i)*Pow(x, j)*Pow(z, k);
		  P(nb_dof_loc++).Init(zero, Q, zero);
		  Q = Pow(z, i)*Pow(x, j)*Pow(y, k);
		  P(nb_dof_loc++).Init(zero, zero, Q);
		}
	  
	  for (int i = 0; i <= r; i++)
	    for (int k = 0; k <= r-1; k++)
	      {
		Q = Pow(x, i)*Pow(y, r)*Pow(z, k);
		P(nb_dof_loc++).Init(Q, zero, zero);
		Q = Pow(x, i)*Pow(y, k)*Pow(z, r);
		P(nb_dof_loc++).Init(Q, zero, zero);
		
		Q = Pow(x, r)*Pow(y, i)*Pow(z, k);
		P(nb_dof_loc++).Init(zero, Q, zero);
		Q = Pow(x, k)*Pow(y, i)*Pow(z, r);
		P(nb_dof_loc++).Init(zero, Q, zero);
		
		Q = Pow(x, r)*Pow(y, k)*Pow(z, i);
		P(nb_dof_loc++).Init(zero, zero, Q);
		Q = Pow(x, k)*Pow(y, r)*Pow(z, i);
		P(nb_dof_loc++).Init(zero, zero, Q);
	      }
	  
	  for (int i = 0; i <= r-1; i++)
	    {
	      P(nb_dof_loc++).Init(Pow(x, r+1)*Pow(y, r)*Pow(z, i), -Pow(x, r)*Pow(y, r+1)*Pow(z, i), zero);
	      P(nb_dof_loc++).Init(Pow(x, r+1)*Pow(y, i)*Pow(z, r), zero, -Pow(x, r)*Pow(y, i)*Pow(z, r+1));
	      P(nb_dof_loc++).Init(zero, Pow(x, i)*Pow(y, r+1)*Pow(z, r), -Pow(x, i)*Pow(y, r)*Pow(z, r+1));
	    }
	}
      */

      for (int i = 0; i <= r; i++)
        for (int j = 0; j <= r; j++)
          if ((i != r) || (j != r))
            {
              P(nb_dof_loc++).Init((one - x)*Pow(y, i)*Pow(z, j), zero, zero);
              P(nb_dof_loc++).Init(x*Pow(y, i)*Pow(z, j), zero, zero);
              P(nb_dof_loc++).Init(zero, (one - y)*Pow(x, i)*Pow(z, j), zero);
              P(nb_dof_loc++).Init(zero, y*Pow(x, i)*Pow(z, j), zero);
              P(nb_dof_loc++).Init(zero, zero, (one - z)*Pow(x, i)*Pow(y, j));
              P(nb_dof_loc++).Init(zero, zero, z*Pow(x, i)*Pow(y, j));
            }
      
      for (int i = 0; i <= r; i++)
        for (int j = 0; j <= r; j++)
          for (int k = 0; k <= r; k++)
            if ( (i == r) || (j != r) || (k != r))
              {
                P(nb_dof_loc++).Init(x*(one-x)*Pow(x, i)*Pow(y, j)*Pow(z, k), zero, zero);
                P(nb_dof_loc++).Init(zero, y*(one-y)*Pow(y, i)*Pow(x, j)*Pow(z, k), zero);
                P(nb_dof_loc++).Init(zero, zero, z*(one-z)*Pow(z, i)*Pow(x, j)*Pow(y, k));
              }
      
      P(nb_dof_loc++).Init((one-x)*Pow(x, r)*Pow(y, r)*Pow(z, r), (one-y)*Pow(x, r)*Pow(y, r)*Pow(z, r), (one-z)*Pow(x, r)*Pow(y, r)*Pow(z, r));
      
      DISP(nb_dof_loc);
      for (int i = 0; i < nb_dof_loc; i++)
	{
	  DISP(i);
	  DISP(P(i));
	}

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
      
      // DISP(nb_dof_loc);
      //DISP(P);      
      
      // for projection on dofs, we use quadrature points
      VectR3 points_dof3d;
      points_dof3d.Reallocate(nb_points_quadrature_inside);
      for (int i = 0; i < nb_points_quadrature_inside; i++)
        points_dof3d(i) = this->PointsND(i);
      
      this->SetPointsDofND(points_dof3d);
      // DISP(points_dof3d);
      
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
	    //DISP(i); DISP(lambda(i));
	    Qx = zero; Qy = zero; Qz = zero;
	    for (int j = 0; j < lambda.GetM(); j++)
	      if (abs(V(j, i)) > 1e-5)
		{
		  DISP(j); DISP(V(j, i)); DISP(P(j));
		  Qx += V(j, i)*P(j)(0);
		  Qy += V(j, i)*P(j)(0);
		  Qz += V(j, i)*P(j)(0);
		}
	    
	    //DISP(Qx); DISP(Qy); DISP(Qz);
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

      this->elt_geom.ComputeCoefficientTransformation();
      
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
      int nb_fcts_wr = (r+1)*(r+1)*(r+1) + 3*(r+1)*(r+1);
      MultivariatePolynomial<Real_wp> x, y, z, Q;
      x.SetOrder(3, 1); x(1, 0, 0) = 1;
      y.SetOrder(3, 1); y(0, 1, 0) = 1;
      z.SetOrder(3, 1); z(0, 0, 1) = 1;
      
      Wr.Reallocate(nb_fcts_wr);
      nb_fcts_wr = 0;
	
      for (int i = 0; i <= r; i++)
	for (int j = 0; j <= r; j++)
	  for (int k = 0; k <= r; k++)
	    {
	      bool fct_taken = true;
	      if (diff_wr)
		{
		  if ( ((i < r) && (j < r)) || ((i < r) && (k < r)) || ((j < r) && (k < r)) )
		    fct_taken = false;
		  
		  if ((r==1)&& ((i!=0) || (j!=0) || (k!=0)) )
		    fct_taken = true;
		}
	      
	      if (fct_taken)
		{
		  DISP(i); DISP(j); DISP(k);
		  Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
		  Wr(nb_fcts_wr++) = Q;
		}
	    }
      
      for (int i = 0; i <= r; i++)
	for (int j = 0; j <= r; j++)
	  {
	    Q = Pow(x, r+1)*Pow(y, i)*Pow(z, j);
	    Wr(nb_fcts_wr++) = Q;
	    Q = Pow(y, r+1)*Pow(x, i)*Pow(z, j);
	    Wr(nb_fcts_wr++) = Q;
	    Q = Pow(z, r+1)*Pow(x, i)*Pow(y, j);
	    Wr(nb_fcts_wr++) = Q;
	  }
      
      if (diff_wr)
	Wr.Resize(nb_fcts_wr);
    }


    void GetDhatR(int r, Vector<TinyVector<int, 3> >& Wr, bool diff_wr = false)
    {
      int nb_fcts_wr = (r+1)*(r+1)*(r+1) + 3*(r+1)*(r+1);
      Wr.Reallocate(nb_fcts_wr);
      nb_fcts_wr = 0;
	
      for (int i = 0; i <= r; i++)
	for (int j = 0; j <= r; j++)
	  for (int k = 0; k <= r; k++)
	    {
	      bool fct_taken = true;
	      if (diff_wr)
		{
		  if ( ((i < r) && (j < r)) || ((i < r) && (k < r)) || ((j < r) && (k < r)) )
		    fct_taken = false;
		  
		  if ((r==1)&& ((i!=0) || (j!=0) || (k!=0)) )
		    fct_taken = true;
		}
	      
	      if (fct_taken)
		{
		  Wr(nb_fcts_wr).Init(i, j, k);
		  nb_fcts_wr++;
		}
	    }
      
      for (int i = 0; i <= r; i++)
	for (int j = 0; j <= r; j++)
	  {
	    Wr(nb_fcts_wr++).Init(r+1, i, j);
	    Wr(nb_fcts_wr++).Init(i, r+1, j);
	    Wr(nb_fcts_wr++).Init(i, j, r+1);
	  }
      
      if (diff_wr)
	Wr.Resize(nb_fcts_wr);
    }
		  
    void AddFunctionsDivergence(int r)
    {
      MultivariatePolynomial<Real_wp> x, y, z, Q, Q2, zero;
      x.SetOrder(3, 1); x(1, 0, 0) = 1;
      y.SetOrder(3, 1); y(0, 1, 0) = 1;
      z.SetOrder(3, 1); z(0, 0, 1) = 1;
      zero.SetOrder(3, 0); zero(0, 0, 0) = 0;

      // on calcule la difference Wr \ Wr-1
      Vector<TinyVector<int, 3> > Wr;
      GetDhatR(r, Wr, true);
      
      //DISP(Wr.GetM());
      // on "primitive" l'espace Wr pour completer Er et obtenir Pchapeau_r
      for (int n = 0; n < Wr.GetM(); n++)
	{
	  int i = Wr(n)(0);
	  int j = Wr(n)(1);
	  int k = Wr(n)(2);
	  if ((i > j) && (i > k))
	    {
	      Q = Pow(x, i+1)*Pow(y, j)*Pow(z, k);
	      P(nb_dof_loc++).Init(Q, zero, zero);
	    }
	  else if ((j > i) && (j > k))
	    {
	      Q = Pow(x, i)*Pow(y, j+1)*Pow(z, k);
	      P(nb_dof_loc++).Init(zero, Q, zero);
	    }
	  else if ((k > i) && (k > j))
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(z, k+1);
	      P(nb_dof_loc++).Init(zero, zero, Q);
	    }
	  else if ((i == j) && (i > k))
	    {
	      Q = Pow(x, i+1)*Pow(y, j)*Pow(z, k);
	      Q2 = Pow(x, i)*Pow(y, j+1)*Pow(z, k);
	      P(nb_dof_loc++).Init(Q, Q2, zero);
	    }
	  else if ((i == k) && (i > j))
	    {
	      Q = Pow(x, i+1)*Pow(y, j)*Pow(z, k);
	      Q2 = Pow(x, i)*Pow(y, j)*Pow(z, k+1);
	      P(nb_dof_loc++).Init(Q, zero, Q2);
	    }
	  else if ((j == k) && (j > i))
	    {
	      Q = Pow(x, i)*Pow(y, j+1)*Pow(z, k);
	      Q2 = Pow(x, i)*Pow(y, j)*Pow(z, k+1);
	      P(nb_dof_loc++).Init(zero, Q, Q2);
	    }
	  else
	    {
	      P(nb_dof_loc++).Init(Pow(x, i+1)*Pow(y, j)*Pow(z, k), Pow(x, i)*Pow(y, j+1)*Pow(z, k), Pow(x, i)*Pow(y, j)*Pow(z, k+1));
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
      ConstructMassMatrixDiv();
      
      // on calcule Dr_chapeau qu'on note Wr (espace que doit inclure div PrF)
      Vector<MultivariatePolynomial<Real_wp> > Wr, Wr_final;
      GetDhatR(r, Wr);
      int nb_fcts_wr = (r+1)*(r+1)*(r+1) + 3*(r+1)*(r+1);      
      int nb_pts_quad = this->nb_points_quadrature_inside;
      VectReal_wp feval(nb_pts_quad), contrib(nb_dof_loc);
      feval.Fill(0); contrib.Fill(0);
      Wr_final.Reallocate(nb_fcts_wr);
      int nb_final = 0;
      for (int i = 0; i < Pdiv.GetM(); i++)
	{
	  //DISP(i); DISP(Pdiv(i));
	}
      
      // on determine les fonctions de Wr qui ne sont pas dans div Er
      for (int i = 0; i < nb_fcts_wr; i++)
	{
	  for (int j = 0; j < nb_pts_quad; j++)
	    feval(j) = Wr(i).Evaluate(PointsND(j));
	  
	  ComputeProjectionDivRef(feval, contrib);
	  
	  bool fct_in_divEr = true;
	  Real_wp err_max = 0;
	  for (int j = 0; j < nb_pts_quad; j++)
	    {
	      Real_wp val_ex = feval(j);
	      Real_wp val_proj = 0;
	      for (int k = 0; k < nb_dof_loc; k++)
		val_proj += Div_Phi(k, j)*contrib(k);
	      
	      if (abs(val_proj-val_ex) > 1e-8)
		{
		  fct_in_divEr = false;
		  if (abs(val_proj-val_ex) > err_max)
		    err_max = abs(val_proj-val_ex);
		}
	    }
	  
	  if (!fct_in_divEr)
	    {
	      DISP(i); DISP(Wr(i)); DISP(err_max);
	      Wr_final(nb_final) = Wr(i);
	      nb_final++;
	    }
	}     
                  
      Wr_final.Resize(nb_final);
      // on affiche le nombre de ddls de l'espace optimal : nb_dof_loc + nb_final
      //DISP(nb_final); DISP(nb_dof_loc+nb_final); DISP(3*(r+3)*(r+1)*(r+1));
      //DISP(3*(r+3)*(r+1)*(r+1) - (nb_dof_loc+nb_final));
      
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
    
    void FindMonomialBasis(int r)
    {
      MultivariatePolynomial< Real_wp > x, y, z, zero, Q;
      TinyVector<MultivariatePolynomial<Real_wp>, 3> Dr;
      x.SetOrder(3, 1); x(1, 0, 0) = 1;
      y.SetOrder(3, 1); y(0, 1, 0) = 1;
      z.SetOrder(3, 1); z(0, 0, 1) = 1;
      zero.SetOrder(3, 0); zero(0,0,0) = 0;
      
      int nb_pts_quad = nb_points_quadrature_inside;
      VectR3 feval(nb_pts_quad); VectReal_wp contrib(nb_dof_loc);
      Real_wp xp, yp, zp;
      Matrix<R3> ValPhiDr(nb_dof_loc, nb_pts_quad);
      ValPhiDr.Fill(0);
      int nb_fcts_dr = 0;
      //cout << endl << "The following functions of Dr are not belonging to Pr" <<endl;
      cout << endl << "The following functions of Dr are belonging to Pr" <<endl;
      for (int i = 0; i <= r+2; i++)
	for (int j = 0; j <= r; j++)
	  for (int k = 0; k <= r; k++)
	    for (int np = 0; np < 3; np++)
	      {
		if (np == 0)
		  {
		    Q = Pow(x, i)*Pow(y, j)*Pow(z, k);
		    Dr.Init(Q, zero, zero);
		  }
		else if (np == 1)
		  {
		    Q = Pow(x, j)*Pow(y, i)*Pow(z, k);
		    Dr.Init(zero, Q, zero);
		  }
		else
		  {
		    Q = Pow(x, j)*Pow(y, k)*Pow(z, i);
		    Dr.Init(zero, zero, Q);
		  }
		
		for (int n = 0; n < nb_pts_quad; n++)
		  {
		    xp = PointsND(n)(0);
		    yp = PointsND(n)(1);
		    zp = PointsND(n)(2);
		    feval(n).Fill(0);
		    if (np == 0)			
		      feval(n)(0) = pow(xp, i)*pow(yp, j)*pow(zp, k);
		    else if (np == 1)
		      feval(n)(1) = pow(xp, j)*pow(yp, i)*pow(zp, k);
		    else if (np == 2)
		      feval(n)(2) = pow(xp, j)*pow(yp, k)*pow(zp, i);
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
		    //cout << Dr << endl;
		    //DISP(err_max);
		  }
		else
		  {
		    cout << Dr << endl;
		    for (int n = 0; n < nb_dof_loc; n++)
		      if (abs(contrib(n)) > 1e9*epsilon_machine)
			{
			  DISP(contrib(n));
			  DISP(P(n));
			}
		    
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
	    mat_mass_Dr(i, j) += this->WeightsND(k)*DotProd(ValPhiDr(i, k), ValPhiDr(j, k));
      
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
	      contrib(k) += this->WeightsND(n)*DotProd(ValPhiDr(k, n), feval(n));
	  
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
      
      //exit(0);
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
  
  HexahedronHdivFirstFamilyBis hex;
  hex.ConstructFiniteElement(r);
  //hex.FindMonomialBasis(r);
  
  Mesh<Dimension3> mesh;
  
  mesh.ReallocateVertices(8);
  mesh.ReallocateElements(1);
  mesh.Element(0).InitHexahedral(0, 1, 2, 3, 4, 5, 6, 7, 1);
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
      R3 pt; Real_wp h = 1.0;
      
      // first test : we check if we can generate Pk with basis functions
      //bool fct_needed = false;
      for (int num = 1; num <= Pk.GetM(); num++)
	{
	  mesh.Vertex(0).Init(-0.9*h, -1.2*h, -0.4*h);
	  mesh.Vertex(1).Init(1.2*h, -0.95*h, -0.2*h);
	  mesh.Vertex(2).Init(1.1*h, 0.97*h, 0.1*h); 
	  mesh.Vertex(3).Init(-0.87*h, 1.05*h, -0.3*h);
	  mesh.Vertex(4).Init(-0.95*h, -1.08*h, 0.97*h);
	  mesh.Vertex(5).Init(1.13*h, -0.96*h, 1.21*h);
	  mesh.Vertex(6).Init(1.07*h, 0.99*h, 1.1*h); 
	  mesh.Vertex(7).Init(-0.91*h, 1.02*h, 1.24*h);
	  
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
