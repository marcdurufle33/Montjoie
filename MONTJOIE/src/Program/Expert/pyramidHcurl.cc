#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_HCURL

#include "Montjoie.hxx"

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

// curl of this function (for error with curl)
R3 fcurl(const R3& x)
{
  R3 u;
  u.Zero();
  
  return u; 
}


namespace Montjoie
{
  //! edge element for pyramid
  class PyramidHcurlFirst : public PyramidReference<2>
  {
    public :
    // we storage basis functions in P
    Vector< TinyVector< MultivariatePolynomial< Real_wp >, 3 > > P;
    // and curl of basis functions
    Vector< TinyVector< MultivariatePolynomial< Real_wp >, 3 > > Pcurl;
    Matrix<R3> Value_Phi, Curl_Phi;
    int type_basis;
    enum {CLASSICAL, MONOMIAL, OPTIMAL};
    
    PyramidHcurlFirst() : PyramidReference()
    {
      this->elt_geom.dof_equal_nodal = false;
      // type_basis = OPTIMAL;
      // type_basis = CLASSICAL;
	  type_basis = MONOMIAL;
    }
    
    // construction of basis functions
    void ConstructFiniteElement(int r)
    {
      PyramidReference::ConstructFiniteElement(r+1);
      
      MultivariatePolynomial< Real_wp > x, y, z, xp, yp, one, zero;
      x.SetOrder(3,1); x(1,0,0) = 1;
      y.SetOrder(3,1); y(0,1,0) = 1;
      z.SetOrder(3,1); z(0,0,1) = 1;
      one.SetOrder(3,0); one(0,0,0) = 1;
      zero.SetOrder(3,0); zero(0,0,0) = 0;
      xp = 2.0*x - one;
      yp = 2.0*y - one;

      if (type_basis == CLASSICAL)
        {
          int n = r*r*r + 7*r*(r+1)/2;
          P.Reallocate(n);
                    
          int ind = 0;
          for (int k = 0; k <= r ; k++)
            for (int i = 0; i <= k ; i++)
              for (int j = 0; j <= k ; j++)
                {
                  if (i > 0)
                    {
                      P(ind)(0) = Pow(x,i-1)*Pow(y,j)*Pow(one-z,k);
                      P(ind)(1) = zero;
                      P(ind)(2) = zero;
                      ind++;
                    }
                  
                  if (j > 0)
                    {
                      P(ind)(0) = zero;
                      P(ind)(1) = Pow(x,i)*Pow(y,j-1)*Pow(one-z,k);
                      P(ind)(2) = zero;
                      ind++;
                    }
                  
                  if (k > 0)
                    {
                      P(ind)(0) = zero;
                      P(ind)(1) = zero;
                      P(ind)(2) = Pow(x,i)*Pow(y,j)*Pow(one-z,k-1);
                      ind++;
                    }
                }
        }
	if (type_basis == MONOMIAL)
        {
          int n = r*(r+3)*(2*r+3)/2;
          //DISP(n);
          MultivariatePolynomial<Real_wp> R;
          P.Reallocate(n); n = 0;
          for (int k = 0; k <= r-1; k++)
            for (int i = 0; i <= k; i++)
              for (int j = 0; j <= k; j++)
                {
                  R = Pow(xp, i)*Pow(yp, j)*Pow(one-z, k);
                  P(n).Init(R, zero, zero); n++;
                  P(n).Init(zero, R, zero); n++;
                  P(n).Init(zero, zero, R); n++;
                }
          
          for (int k = 0; k <= r-2; k++)
            for (int j = 0; j <= k; j++)
              {
                int i = k+2;
                P(n)(0) = Pow(xp, j)*Pow(yp, i)*Pow(one-z, k+1);
                P(n)(1) = zero;
                P(n)(2) = Pow(xp, j+1)*Pow(yp, i)*Pow(one-z, k+1);
                n++;
                
                P(n)(0) = zero;
                P(n)(1) = Pow(xp, i)*Pow(yp, j)*Pow(one-z, k+1);
                P(n)(2) = Pow(xp, i)*Pow(yp, j+1)*Pow(one-z, k+1);
                n++;
              }
          
          int k = r-1;
          for (int j = 0; j <= k; j++)
            for (int i = 0; i <= k+2; i++)
              {
                P(n)(0) = Pow(xp, j)*Pow(yp, i)*Pow(one-z, k+1);
                P(n)(1) = zero;
                P(n)(2) = Pow(xp, j+1)*Pow(yp, i)*Pow(one-z, k+1);
                n++;
                
                P(n)(0) = zero;
                P(n)(1) = Pow(xp, i)*Pow(yp, j)*Pow(one-z, k+1);
                P(n)(2) = Pow(xp, i)*Pow(yp, j+1)*Pow(one-z, k+1);
                n++;
              }
          
         for (int k = 1; k <= r; k++)
            {              
              P(n).Init(Pow(xp, k-1)*Pow(yp, k)*Pow(one-z, k-1),
                        Pow(xp, k)*Pow(yp, k-1)*Pow(one-z, k-1), Pow(xp*yp, k)*Pow(one-z, k-1));
              n++;
            }
        }
	else
		{
				
		}

      nb_dof_loc = P.GetM();
      Pcurl.Reallocate(nb_dof_loc);
      for (int i = 0; i < nb_dof_loc; i++)
        {
          DISP(i);
          DISP(P(i));
          GetCurlPolynomial(P(i), Pcurl(i));
        }         
      
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

      // computation of phi_i(\xi_k) for all i
      Curl_Phi.Reallocate(nb_dof_loc, this->GetNbPointsQuadrature());
      for (int k = 0; k < this->GetNbPointsQuadrature(); k++)
        {
          VectR3 phi;
          ComputeCurlPhiRef(this->PointsND(k), phi);
          for (int i = 0; i < nb_dof_loc; i++)
            Curl_Phi(i, k) = phi(i);
        }
      
      // computation of mass matrix
      nb_dof_loc = P.GetM();
      mass_matrix.Reallocate(nb_dof_loc, nb_dof_loc);
      mass_matrix.Fill(0);
      for (int k = 0; k < nb_points_quadrature_inside; k++)
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = i; j < nb_dof_loc; j++)
            mass_matrix(i, j) += this->WeightsND(k)*DotProd(Value_Phi(i, k), Value_Phi(j, k));
      
      mass_matrix.WriteText("mat_mass.dat");
      
      // inversion
      mass_matrix_chol = mass_matrix;
      GetCholesky(mass_matrix_chol);
      
      this->elt_geom.ComputeCoefficientTransformation();
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
      // we are using xtilde coordinates (on the cube)
      Real_wp xt = 0.5 + 0.5*pt_loc(0)/(1.0-pt_loc(2));
      Real_wp yt = 0.5 + 0.5*pt_loc(1)/(1.0-pt_loc(2));
      Real_wp zt = pt_loc(2);
      R3 pt(xt, yt, zt);
      
      R3 vec_u, vec_v; phi.Reallocate(nb_dof_loc);
      if (type_basis == CLASSICAL)
        {
          // DF of transformation (x = (2 xt -1 ) (1-zt), y = (2 yt - 1) (1-zt), z = zt)
          // is equal to | 2 (1-zt)   0        (1 - 2 xt) |
          //              |   0     2 (1-zt)   (1 - 2 yt) |
          //              |   0       0            1      |
          Matrix3_3 DF, DFm1;
          DF(0, 0) = 2.0*(1.0-zt); DF(0, 2) = 1.0 - 2.0*xt;
          DF(1, 1) = 2.0*(1.0-zt); DF(1, 2) = 1.0 - 2.0*yt;
          DF(2, 2) = 1.0;
          
          GetInverse(DF, DFm1);
          
          for (int i = 0; i < nb_dof_loc; i++)
            {
              vec_u(0) = P(i)(0).Evaluate(pt);
              vec_u(1) = P(i)(1).Evaluate(pt);
              vec_u(2) = P(i)(2).Evaluate(pt);
              
              MltTrans(DFm1, vec_u, vec_v);
              
              phi(i) = vec_v;
            }         
        }
      else
        for (int i = 0; i < nb_dof_loc; i++)
          {
            vec_u(0) = P(i)(0).Evaluate(pt);
            vec_u(1) = P(i)(1).Evaluate(pt);
            vec_u(2) = P(i)(2).Evaluate(pt);
            
            phi(i) = vec_u;
          }         
      
    }
    
    // computation of curl of basis functions
    void ComputeCurlPhiRef(const R3& pt_loc, VectR3& phi) const
    {
      phi.Reallocate(nb_dof_loc);
      FillZero(phi);
    }
    
  };

}


void GetRkBasis(int r, Vector<TinyVector<MultivariatePolynomial<Real_wp>, 3> >& psi)
{
  int nb_fcts = r*(r+2)*(r+3)/2;

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
  
  // then S_r
  for (int i = 0; i < r; i++)
    for (int j = 0; j < r-i; j++)
      {
        psi(node)(0) = Pow(x, i)*Pow(y, j)*Pow(z, r-i-j);
        psi(node)(1) = zero;
        psi(node)(2) = -Pow(x, i+1)*Pow(y, j)*Pow(z, r-i-j-1);
        node++;
      }
  
  
  for (int i = 0; i < r; i++)
    for (int j = 0; j < r; j++)
      if (i+j == r-1)
        {
          psi(node)(0) = Pow(x, i)*Pow(y, j+1);
          psi(node)(1) = -Pow(x, i+1)*Pow(y, j);
          psi(node)(2) = zero;
          node++;
        }
  
  for (int i = 0; i < r; i++)
    for (int j = 0; j < r-i; j++)
      {
        psi(node)(0) = zero;
        psi(node)(1) = Pow(x, i)*Pow(y, j)*Pow(z, r-i-j);
        psi(node)(2) = -Pow(x, i)*Pow(y, j+1)*Pow(z, r-i-j-1);
        node++;
      }  
}

void GetHrotBasis(int r, Vector<TinyVector<MultivariatePolynomial<Real_wp>, 3> >& psi)
{
	MultivariatePolynomial< Real_wp > x, y, z, xp, yp, one, zero;
	x.SetOrder(3,1); x(1,0,0) = 1;
	y.SetOrder(3,1); y(0,1,0) = 1;
	z.SetOrder(3,1); z(0,0,1) = 1;
	one.SetOrder(3,0); one(0,0,0) = 1;
	zero.SetOrder(3,0); zero(0,0,0) = 0;
	xp = 2.0*x - one;
	yp = 2.0*y - one;

	int n = r*(r+3)*(2*r+3)/2;
	//DISP(n);
	MultivariatePolynomial<Real_wp> R;
	psi.Reallocate(n); n = 0;
	for (int k = 0; k <= r-1; k++)
		for (int i = 0; i <= k; i++)
			for (int j = 0; j <= k; j++)
			{
				R = Pow(xp, i)*Pow(yp, j)*Pow(one-z, k);
				psi(n).Init(R, zero, zero); n++;
				psi(n).Init(zero, R, zero); n++;
				psi(n).Init(zero, zero, R); n++;
			}
	
	for (int k = 0; k <= r-2; k++)
		for (int j = 0; j <= k; j++)
		{
			int i = k+2;
			psi(n)(0) = Pow(xp, j)*Pow(yp, i)*Pow(one-z, k+1);
			psi(n)(1) = zero;
			psi(n)(2) = Pow(xp, j+1)*Pow(yp, i)*Pow(one-z, k+1);
			n++;
			
			psi(n)(0) = zero;
			psi(n)(1) = Pow(xp, i)*Pow(yp, j)*Pow(one-z, k+1);
			psi(n)(2) = Pow(xp, i)*Pow(yp, j+1)*Pow(one-z, k+1);
			n++;
		}
	
	int k = r-1;
	for (int j = 0; j <= k; j++)
		for (int i = 0; i <= k+2; i++)
		{
			psi(n)(0) = Pow(xp, j)*Pow(yp, i)*Pow(one-z, k+1);
			psi(n)(1) = zero;
			psi(n)(2) = Pow(xp, j+1)*Pow(yp, i)*Pow(one-z, k+1);
			n++;
			
			psi(n)(0) = zero;
			psi(n)(1) = Pow(xp, i)*Pow(yp, j)*Pow(one-z, k+1);
			psi(n)(2) = Pow(xp, i)*Pow(yp, j+1)*Pow(one-z, k+1);
			n++;
		}
	
	for (int k = 1; k <= r; k++)
	{              
		psi(n).Init(Pow(xp, k-1)*Pow(yp, k)*Pow(one-z, k-1),
				  Pow(xp, k)*Pow(yp, k-1)*Pow(one-z, k-1), Pow(xp*yp, k)*Pow(one-z, k-1));
		n++;
	}
}


int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  int r = 1; int r_test = r;
  
  PyramidHcurlFirst pyr;
  pyr.ConstructFiniteElement(r);

  Mesh<Dimension3> mesh;
  
  mesh.ReallocateVertices(5);
  mesh.ReallocateElements(1);
  mesh.Element(0).InitPyramidal(0, 1, 2, 3, 4, 1);
    
  mesh.ReorientElements();

  VectR3 s; SetPoints<Dimension3> PointsElem;
  SetMatrices<Dimension3> MatricesElem;
  int nb_pts_quad = pyr.GetNbPointsQuadratureInside();
  int nb_pts_dof = pyr.GetNbPointsDof();
  int nb_dof_loc = pyr.GetNbDof();
  VectR3 feval(nb_pts_dof);
  VectReal_wp res(nb_dof_loc);
  R3 vec_u, u_approche, u_exact, curl_exact, curl_approche;
  VectR3 Points(nb_pts_quad);
  for (int i = 0; i < nb_pts_quad; i++)
    Points(i) = pyr.PointsND(i);
  
  Matrix3_3 dfjm1; 
  // Real_wp jacob;
  Vector<TinyVector<MultivariatePolynomial<Real_wp>, 3> > Rk;
  GetRkBasis(r_test, Rk);
  R3 pt; Real_wp h = 1.0;
  // first test : we check if we can generate Rk with basis functions
  for (int num = 1; num <= Rk.GetM(); num++)
    {
		mesh.Vertex(0).Init(-0.9*h, -1.2*h, -0.4*h);
		mesh.Vertex(1).Init(1.2*h, -0.95*h, -0.2*h);
		mesh.Vertex(2).Init(1.1*h, 0.97*h, 0.1*h); 
		mesh.Vertex(3).Init(-0.87*h, 1.05*h, -0.3*h);
		mesh.Vertex(4).Init(0, 0, h);
      
      // boucle sur les elements
      for (int i = 0; i < mesh.GetNbElt(); i++)
        {
          // on calcule F_i(x), DF_i(x)
          mesh.GetVerticesElement(i, s);
          pyr.FjElem(s, PointsElem, mesh, i);
          pyr.DFjElem(s, PointsElem, MatricesElem, mesh, i);
          
          // on calcule P = DF* p avec p polynome de Rk
          for (int k = 0; k < nb_pts_dof; k++)
            {
              pt = PointsElem.GetPointDof(k);
              vec_u(0) = Rk(num-1)(0).Evaluate(pt);
              vec_u(1) = Rk(num-1)(1).Evaluate(pt);
              vec_u(2) = Rk(num-1)(2).Evaluate(pt);
              MltTrans(MatricesElem.GetPointDof(k), vec_u, feval(k));
            }
          
          // et on projette sur les fonctions de base de l'element de reference
          pyr.ComputeProjectionDofRef(feval, res);
          
          // on calcule la norme L2 entre P = DF* p et son projete
          Real_wp erreur = 0, norme_L2 = 0; R3 evalP, eval_proj;
          for (int k = 0; k < nb_pts_quad; k++)
            {
              eval_proj.Zero();
              for (int i = 0; i < nb_dof_loc; i++)
                eval_proj += res(i)*pyr.Value_Phi(i, k);
              
              // evaluation de P
              pt = PointsElem.GetPointQuadrature(k);
              vec_u(0) = Rk(num-1)(0).Evaluate(pt);
              vec_u(1) = Rk(num-1)(1).Evaluate(pt);
              vec_u(2) = Rk(num-1)(2).Evaluate(pt);
              MltTrans(MatricesElem.GetPointQuadrature(k), vec_u, evalP);
              
              erreur += pyr.WeightsND(k)*DotProd(evalP - eval_proj, evalP - eval_proj);
              norme_L2 += pyr.WeightsND(k)*DotProd(evalP, evalP);
            }
          
          erreur = sqrt(erreur);
          norme_L2 = sqrt(norme_L2);
          DISP(num); DISP(Rk(num-1)); DISP(erreur/norme_L2);
          
        }
    }
    
  return FinalizeMontjoie();
}
