#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_HCURL

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

// curl of this function (for error with curl)
R3 fcurl(const R3& x)
{
  R3 u;
  u.Zero();
  
  return u; 
}


namespace Montjoie
{
  //! edge element for wedamid
  class WedgeHcurlFirstFamilyBis : public WedgeReference<2>
  {
    public :
    // we storage basis functions in P
    Vector< TinyVector< MultivariatePolynomial< Real_wp >, 3 > > P;
    // and curl of basis functions
    Vector< TinyVector< MultivariatePolynomial< Real_wp >, 3 > > Pcurl;
    Matrix<R3> Value_Phi, Curl_Phi;
    int type_basis;
    enum {CLASSICAL, OPTIMAL};
    
    WedgeHcurlFirstFamilyBis() : WedgeReference()
    {
      this->elt_geom.dof_equal_nodal = false;
      type_basis = OPTIMAL;
      // type_basis = CLASSICAL;
    }
    
    // construction of basis functions
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1)
    {
      WedgeReference::ConstructFiniteElement(r+1);
      
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
          int n = 5*r*(r+1)*(r+2)/2;
          P.Reallocate(n);
                    
          int ind = 0;
          for (int k = 0; k <= r ; k++)
            {
              for (int i = 0; i <= r-1 ; i++)
                for (int j = 0; j <= r-1 ; j++)
                  {
                    if (i+j <= r-1)
                      {
                        P(ind)(0) = Pow(x,i)*Pow(y,j);
                        P(ind)(1) = zero;
                        P(ind)(2) = Pow(z,k);
                        ind++;
                      
                        P(ind)(0) = zero;
                        P(ind)(1) = Pow(x,i)*Pow(y,j);
                        P(ind)(2) = Pow(z,k);
                        ind++;
                      }
                    if (i+j == r-1)
                      {
                        P(ind)(0) = Pow(x,i)*Pow(y,j+1);
                        P(ind)(1) = -Pow(x,i+1)*Pow(y,j);
                        P(ind)(2) = Pow(z,k);
                        ind++;
                      } 
                  }
                  
            }

          for (int k = 0; k <= r-1 ; k++)
            {
              for (int i = 0; i <= r ; i++)
                for (int j = 0; j <= r ; j++)
                  if (i+j <= r)
                    {
                      P(ind)(0) = Pow(x,i)*Pow(y,j);
                      P(ind)(1) = zero;
                      P(ind)(2) = Pow(z,k);
                      ind++;
                      
                      P(ind)(0) = zero;
                      P(ind)(1) = Pow(x,i)*Pow(y,j);
                      P(ind)(2) = Pow(z,k);
                      ind++;
                    }   
                 
              }
        }
      else
        {
          int n = r*(r+2)*(3*r+7)/2;
		  P.Reallocate(n);
	  
          int ind = 0;
          for (int k = 0; k <= r-1 ; k++)
            for (int i = 0; i <= r+1 ; i++)
              for (int j = 0; j <= r+1 ; j++)
                if (i+j <= r+1)
                  {
                    P(ind)(0) = zero;
                    P(ind)(1) = zero;
                    P(ind)(2) = Pow(x,i)*Pow(y,j)*Pow(z,k);
                    ind++;
                  }
                
          for (int k = 0; k <= r+1 ; k++)
            for (int i = 0; i <= r-1 ; i++)
              for (int j = 0; j <= r-1 ; j++)
                if (i+j <= r-1)
                  {
                    P(ind)(0) = Pow(x,i)*Pow(y,j)*Pow(z,k);
                    P(ind)(1) = zero;
                    P(ind)(2) = zero;
                    ind++;

                    P(ind)(0) = zero;
                    P(ind)(1) = Pow(x,i)*Pow(y,j)*Pow(z,k);
                    P(ind)(2) = zero;
                    ind++;
                  }

          for (int k = 0; k <= r+1 ; k++)
            for (int i = 0; i <= r-1 ; i++)
              for (int j = 0; j <= r-1 ; j++)
                if (i+j == r-1)
                  {
                    P(ind)(0) = Pow(x,j)*Pow(y,i+1)*Pow(z,k);
                    P(ind)(1) = -Pow(x,j+1)*Pow(y,i)*Pow(z,k);
                    P(ind)(2) = zero;
                    ind++;
                  }
               
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
      R3 vec_u, vec_v; phi.Reallocate(nb_dof_loc);
     
      for (int i = 0; i < nb_dof_loc; i++)
        {
          vec_u(0) = P(i)(0).Evaluate(pt_loc);
          vec_u(1) = P(i)(1).Evaluate(pt_loc);
          vec_u(2) = P(i)(2).Evaluate(pt_loc);
            
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


int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  int r = 1; int r_test = r;
  
  WedgeHcurlFirstFamilyBis wed;
  wed.ConstructFiniteElement(r);

  Mesh<Dimension3> mesh;
  
  mesh.ReallocateVertices(6);
  mesh.ReallocateElements(1);
  mesh.Element(0).InitWedge(0, 1, 2, 3, 4, 5, 1);
  
  mesh.ReorientElements();

  VectR3 s; SetPoints<Dimension3> PointsElem;
  SetMatrices<Dimension3> MatricesElem;
  int nb_pts_quad = wed.GetNbPointsQuadratureInside();
  int nb_pts_dof = wed.GetNbPointsDof();
  int nb_dof_loc = wed.GetNbDof();
  VectR3 feval(nb_pts_dof);
  VectReal_wp res(nb_dof_loc);
  R3 vec_u, u_approche, u_exact, curl_exact, curl_approche;
  VectR3 Points(nb_pts_quad);
  for (int i = 0; i < nb_pts_quad; i++)
    Points(i) = wed.PointsND(i);
  
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
      
      mesh.Vertex(3).Init(-0.9*h, -1.2*h, 1.4*h);
      mesh.Vertex(4).Init(1.2*h, -0.95*h, 1.2*h);
      mesh.Vertex(5).Init(1.1*h, 0.97*h, 0.9*h);
      
      // boucle sur les elements
      for (int i = 0; i < mesh.GetNbElt(); i++)
        {
          // on calcule F_i(x), DF_i(x)
          mesh.GetVerticesElement(i, s);
          wed.FjElem(s, PointsElem, mesh, i);
          wed.DFjElem(s, PointsElem, MatricesElem, mesh, i);
          
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
          wed.ComputeProjectionDofRef(feval, res);
          
          // on calcule la norme L2 entre P = DF* p et son projete
          Real_wp erreur = 0, norme_L2 = 0; R3 evalP, eval_proj;
          for (int k = 0; k < nb_pts_quad; k++)
            {
              eval_proj.Zero();
              for (int i = 0; i < nb_dof_loc; i++)
                eval_proj += res(i)*wed.Value_Phi(i, k);
              
              // evaluation de P
              pt = PointsElem.GetPointQuadrature(k);
              vec_u(0) = Rk(num-1)(0).Evaluate(pt);
              vec_u(1) = Rk(num-1)(1).Evaluate(pt);
              vec_u(2) = Rk(num-1)(2).Evaluate(pt);
              MltTrans(MatricesElem.GetPointQuadrature(k), vec_u, evalP);
              
              erreur += wed.WeightsND(k)*DotProd(evalP - eval_proj, evalP - eval_proj);
              norme_L2 += wed.WeightsND(k)*DotProd(evalP, evalP);
            }
          
          erreur = sqrt(erreur);
          norme_L2 = sqrt(norme_L2);
          DISP(num); DISP(Rk(num-1)); DISP(erreur/norme_L2);
          
        }
    }
    
  return FinalizeMontjoie();
}
