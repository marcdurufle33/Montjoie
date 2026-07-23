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
  //! edge element for pyramid
  class PyramidHdivFirstFamilyTest : public PyramidReference<3>
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
    enum {CLASSICAL, MONOMIAL, OPTIMAL};
    
    PyramidHdivFirstFamilyTest() : PyramidReference<3>()
    {
      this->elt_geom.dof_equal_nodal = false;
      // type_basis = OPTIMAL;
      // type_basis = CLASSICAL;
	  type_basis = MONOMIAL;
    }
    
    // construction of basis functions
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1)
    {
      PyramidReference<3>::ConstructFiniteElement(r+1);
      
      MultivariatePolynomial< Real_wp > x, y, z, one, zero, Q;
      x.SetOrder(3, 1); x(1, 0, 0) = 1;
      y.SetOrder(3, 1); y(0, 1, 0) = 1;
      z.SetOrder(3, 1); z(0, 0, 1) = 1;
      one.SetOrder(3, 0); one(0,0,0) = 1;
      zero.SetOrder(3, 0); zero(0,0,0) = 0;
      
      nb_dof_loc = (order+1)*(order+2)*(2*order+3);
      P.Reallocate(nb_dof_loc);
      nb_dof_loc = 0;
      for (int k = 0; k < r; k++)
	for (int i = 0; i <= k; i++)
	  for (int j = 0; j <= k; j++)
	    {
	      Q = Pow(x, i)*Pow(y, j)*Pow(one-z, k);
	      P(nb_dof_loc).Init(Q, zero, zero);
	      P(nb_dof_loc+1).Init(zero, Q, zero);
	      P(nb_dof_loc+2).Init(zero, zero, Q);
              nb_dof_loc += 3;
	    }
      
      for (int k = 0; k <= r-1; k++)
	for (int j = 0; j <= k; j++)
	  {
	    P(nb_dof_loc).Init(Pow(x, k+1)*Pow(y, j)*Pow(one-z, k), zero, zero);
	    P(nb_dof_loc+1).Init(zero, Pow(y, k+1)*Pow(x, j)*Pow(one-z, k), zero);
            nb_dof_loc += 2;
	  }
      
      for (int k = 0; k <= r-1; k++)
	for (int j = 0; j <= k; j++)
	  {
	    P(nb_dof_loc).Init(Pow(x, j+1)*Pow(y, k+1)*Pow(one-z, k), zero, -Pow(x, j)*Pow(y, k+1)*Pow(one-z, k));
	    P(nb_dof_loc+1).Init(zero, Pow(y, j+1)*Pow(x, k+1)*Pow(one-z, k), -Pow(y, j)*Pow(x, k+1)*Pow(one-z, k));
            nb_dof_loc += 2;
	  }
      
      for (int i = 0; i <= r; i++)
        for (int j = 0; j <= r; j++)
          {
            Q = Pow(x, i)*Pow(y, j)*Pow(one-z, r);
            P(nb_dof_loc).Init(Q*x, Q*y, -Q);
            nb_dof_loc++;
          }
      
      DISP(nb_dof_loc);
      P.Resize(nb_dof_loc);
      //DISP(P);
      
      nb_dof_loc = P.GetM();
      Pdiv.Reallocate(nb_dof_loc);
      MultivariatePolynomial<Real_wp> dPx_dx, dPy_dy, dPz_dx, dPz_dy, dPz_dz;
      for (int i = 0; i < nb_dof_loc; i++)
	{
	  DerivatePolynomial(P(i)(0), dPx_dx, 0);
	  DerivatePolynomial(P(i)(1), dPy_dy, 1);
	  DerivatePolynomial(P(i)(2), dPz_dx, 0);
	  DerivatePolynomial(P(i)(2), dPz_dy, 1);
	  DerivatePolynomial(P(i)(2), dPz_dz, 2);
	  Pdiv(i) = dPx_dx + dPy_dy + x*dPz_dx + y*dPz_dy + (one-z)*dPz_dz;
	  DISP(i); DISP(Pdiv(i));
	}
      
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

      // computation of phi_i(\xi_k) for all i
      Div_Phi.Reallocate(nb_dof_loc, this->GetNbPointsQuadrature());
      for (int k = 0; k < this->GetNbPointsQuadrature(); k++)
        {
          VectReal_wp phi;
          ComputeDivPhiRef(this->PointsND(k), phi);
          for (int i = 0; i < nb_dof_loc; i++)
            Div_Phi(i, k) = phi(i);
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
      
      ConstructMassMatrixDiv(r);
      
      elt_geom.ComputeCoefficientTransformation();
    }
    
    
    void ConstructMassMatrixDiv(int r)
    {
      // on construit la matrice \int div(phi_i) div(phi_j)
      mass_matrix_div.Reallocate(nb_dof_loc, nb_dof_loc);
      mass_matrix_div.Fill(0);
      
      for (int k = 0; k < nb_points_quadrature_inside; k++)
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = 0; j < nb_dof_loc; j++)
            mass_matrix_div(i, j) += this->WeightsND(k)*Div_Phi(i, k)*Div_Phi(j, k);
      
      GetQR(mass_matrix_div, tau_div);
      int rank_div = 0;
      for (int k = 0; k < nb_dof_loc; k++)
	if (abs(mass_matrix_div(k, k)) > 1e6*epsilon_machine)
	  rank_div++;
      
      DISP(rank_div);
      int dim_Fr = 0, dim_Frm1 = 0;
      
      for (int k = 0; k <= r-1; k++)
	for (int i = 0; i <= k+1; i++)
	  for (int j = 0; j <= k+1; j++)
	    dim_Fr++;

      for (int k = 0; k <= r-2; k++)
	for (int i = 0; i <= k+1; i++)
	  for (int j = 0; j <= k+1; j++)
	    dim_Frm1++;
      
      DISP(dim_Fr); DISP(dim_Frm1); DISP((r+1)*(r+1)); DISP(dim_Fr-dim_Frm1);
    }

    // projection on reference element
    void ComputeProjectionDofRef(const VectR3& feval, VectReal_wp& res) const
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

    
    // computation of basis functions
    void ComputeValuesPhiRef(const R3& pt_loc, VectR3& phi) const
    {
      phi.Reallocate(nb_dof_loc);
      // we are using xtilde coordinates (on the cube [-1, 1])
      Real_wp xt = pt_loc(0)/(1.0-pt_loc(2));
      Real_wp yt = pt_loc(1)/(1.0-pt_loc(2));
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
      Real_wp xt = pt_loc(0)/(1.0-pt_loc(2));
      Real_wp yt = pt_loc(1)/(1.0-pt_loc(2));
      Real_wp zt = pt_loc(2);
      R3 pt(xt, yt, zt);      
      for (int i = 0; i < nb_dof_loc; i++)
	{
	  phi(i) = Pdiv(i).Evaluate(pt_loc);      
	  phi(i) /= (1.0-zt);
	}
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

  {
    Real_wp h = 1e-8, teta = 0.7751073211093346;
    Real_wp y = (cos(teta-h)-cos(teta+h))/(2.0*h);
    Real_wp y2 = sin(teta)*sin(h)/h;
    DISP(y); DISP(y2);
  }
  
  
  PyramidHdivFirstFamilyTest pyr;
  pyr.ConstructFiniteElement(r);

  Mesh<Dimension3> mesh;
  
  mesh.ReallocateVertices(5);
  mesh.ReallocateElements(1);
  mesh.Element(0).InitPyramidal(0, 1, 2, 3, 4, 1);
  mesh.SetGeometryOrder(r);
  
  mesh.ReorientElements();

  VectR3 s; SetPoints<Dimension3> PointsElem;
  SetMatrices<Dimension3> MatricesElem;
  int nb_pts_quad = pyr.GetNbPointsQuadratureInside();
  int nb_pts_dof = pyr.GetNbPointsDof();
  int nb_dof_loc = pyr.GetNbDof();
  VectR3 feval(nb_pts_dof); VectReal_wp feval_div(nb_pts_dof);
  VectReal_wp res(nb_dof_loc), res_div(nb_dof_loc);
  R3 vec_u, u_approche, u_exact;
  //Real_wp div_exact, div_approche;
  VectR3 Points(nb_pts_quad);
  for (int i = 0; i < nb_pts_quad; i++)
    Points(i) = pyr.PointsND(i);
  
  Matrix3_3 dfjm1; 
  Real_wp jacob;
  Vector<TinyVector<MultivariatePolynomial<Real_wp>, 3> > Dk;
  GetPkBasis(r_test, Dk);
  R3 pt; Real_wp h = 1.0;
  // first test : we check if we can generate Rk with basis functions
  for (int num = 1; num <= Dk.GetM(); num++)
    {
      //mesh.Vertex(0).Init(-h, -h, 0);
      //mesh.Vertex(1).Init(h, -h, 0);
      //mesh.Vertex(2).Init(h, h, 0);
      //mesh.Vertex(3).Init(-h, h, 0);
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
          
          // on calcule P = J DF^-1 p avec p polynome de Rk
          for (int k = 0; k < nb_pts_dof; k++)
            {
              pt = PointsElem.GetPointDof(k);
              vec_u(0) = Dk(num-1)(0).Evaluate(pt);
              vec_u(1) = Dk(num-1)(1).Evaluate(pt);
              vec_u(2) = Dk(num-1)(2).Evaluate(pt);
	      jacob = Det(MatricesElem.GetPointDof(k));
	      GetInverse(MatricesElem.GetPointDof(k), dfjm1);
              Mlt(dfjm1, vec_u, feval(k));
	      feval(k) *= jacob;
	      feval_div(k) = vec_u(0)*jacob;
            }
          
          // et on projette sur les fonctions de base de l'element de reference
          pyr.ComputeProjectionDofRef(feval, res);
          
	  // on projette vis-a-vis de la divergence
	  if (num%3 == 1)
	    {
	      pyr.ComputeProjectionDivRef(feval_div, res_div);
	      //DISP(Dk(num-1)); DISP(res_div);
	    }
	  
          // on calcule la norme L2 entre P = J DF^-1 p et son projete
          Real_wp erreur = 0, norme_L2 = 0; R3 evalP, eval_proj;
          Real_wp divP, div_proj, erreur_div = 0;
	  for (int k = 0; k < nb_pts_quad; k++)
            {
              eval_proj.Zero();
              for (int i = 0; i < nb_dof_loc; i++)
                eval_proj += res(i)*pyr.Value_Phi(i, k);
              
              // evaluation de P
              pt = PointsElem.GetPointQuadrature(k);
              vec_u(0) = Dk(num-1)(0).Evaluate(pt);
              vec_u(1) = Dk(num-1)(1).Evaluate(pt);
              vec_u(2) = Dk(num-1)(2).Evaluate(pt);
	      jacob = Det(MatricesElem.GetPointDof(k));
	      GetInverse(MatricesElem.GetPointDof(k), dfjm1);
              Mlt(dfjm1, vec_u, evalP);
	      evalP *= jacob;
	      
              erreur += pyr.WeightsND(k)*DotProd(evalP - eval_proj, evalP - eval_proj);
              norme_L2 += pyr.WeightsND(k)*DotProd(evalP, evalP);

	      if (num%3 == 1)
		{
		  div_proj = 0;
		  for (int i = 0; i < nb_dof_loc; i++)
		    div_proj += res_div(i)*pyr.Div_Phi(i, k);
		  
		  divP = vec_u(0)*jacob;
		  erreur_div += pyr.WeightsND(k)*square(divP - div_proj);
		}
            }
          
          erreur = sqrt(erreur);
          norme_L2 = sqrt(norme_L2);
          DISP(num); DISP(erreur/norme_L2);
	  if (num%3 == 1)
	    {
	      erreur_div = sqrt(erreur_div);
	      DISP(erreur_div/norme_L2);
	    }	  
        }
    }
    
  return FinalizeMontjoie();
}
