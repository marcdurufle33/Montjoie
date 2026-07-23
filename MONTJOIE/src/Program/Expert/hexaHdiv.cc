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
  u(0) = 1.0 + cos(2.0*x(0) -3.0*x(1) + 4.0*x(2));
  u(1) = 2 + sin(2.5*x(0) + 2.3*x(1) + 3.3*x(2));
  u(2) = 3  - sin(0.5*x(0) -1.2*x(1) + 2.8*x(2));
  
  return u; 
}

// divergence of this function (for error with curl)
Real_wp fdiv(const R3& x)
{
  Real_wp u = -2.0*sin(2.0*x(0)) + 2.3*cos(2.3*x(1)) - 2.8*cos(2.8*x(2));
  
  return u; 
}


namespace Montjoie
{
  //! edge element for hexamid
  class HexahedronHdivFirstFamilyBis : public HexahedronReference<3>
  {
    public :
    // we storage basis functions in P
    Vector< TinyVector< MultivariatePolynomial< Real_wp >, 3 > > P;
    // and divergence of basis functions
    Vector< MultivariatePolynomial< Real_wp > > Pdiv;
    Matrix<R3> Value_Phi;
    Matrix<Real_wp> Div_Phi;
    int type_basis;
    enum {CLASSICAL, OPTIMAL};
    
    HexahedronHdivFirstFamilyBis() : HexahedronReference<3>()
    {
      this->elt_geom.dof_equal_nodal = false;
      type_basis = OPTIMAL;
      //type_basis = CLASSICAL;
    }
    
    // construction of basis functions
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1)
    {
      HexahedronReference::ConstructFiniteElement(r+2);
      
      MultivariatePolynomial< Real_wp > x, y, z, one, zero;
      x.SetOrder(3,1); x(1,0,0) = 1;
      y.SetOrder(3,1); y(0,1,0) = 1;
      z.SetOrder(3,1); z(0,0,1) = 1;
      one.SetOrder(3,0); one(0,0,0) = 1;
      zero.SetOrder(3,0); zero(0,0,0) = 0;
      
      /* nb_dof_loc = 55;
      P.Reallocate(nb_dof_loc);
      P(0).Init(one, zero, zero);
      P(1).Init(zero, one, zero);
      P(2).Init(zero, zero, one);
      
      P(3).Init(x, zero, zero);
      P(4).Init(y, zero, zero);
      P(5).Init(z, zero, zero);
      
      P(6).Init(zero, x, zero);
      P(7).Init(zero, y, zero);
      P(8).Init(zero, z, zero);

      P(9).Init(zero, zero, x);
      P(10).Init(zero, zero, y);
      P(11).Init(zero, zero, z);
      
      P(12).Init(x*y, zero, -y*z);
      P(13).Init(x*z, -y*z, zero);
      P(14).Init(zero, x*y, -x*z);

      P(15).Init(x*x, -x*y, -x*z);
      P(16).Init(-x*y, y*y, -y*z);
      P(17).Init(-x*z, -y*z, z*z);

      P(18).Init(x*x*y, -x*y*y, zero);
      P(19).Init(x*x*z, zero, -x*z*z);
      P(20).Init(zero, y*y*z, -y*z*z);
      
      P(21).Init(zero, zero, x*y);
      P(22).Init(x*x, x*y, zero);
      P(23).Init(x*y, y*y, zero);
      P(24).Init(x*x*y, x*y*y, -x*y*z);
      
      P(25).Init(-x*x, zero, -x*z);
      P(26).Init(zero, x*z, zero);
      P(27).Init(-x*z, zero, -z*z);
      P(28).Init(x*x*z, -x*y*z, x*z*z);

      P(29).Init(zero, y*y, y*z);
      P(30).Init(zero, y*z, z*z);
      P(31).Init(y*z, zero, zero);
      P(32).Init(x*y*z, -y*y*z, -y*z*z);

      P(33).Init(x*x*x, zero, zero);
      P(34).Init(x*x*y, zero, x*y*z);
      P(35).Init(-x*x*z, -x*y*z, zero);
      P(36).Init(zero, -x*y*y, -x*y*z);
      P(37).Init(zero, y*y*y, zero);
      P(38).Init(x*y*z, y*y*z, zero);
      P(39).Init(zero, x*y*z, x*z*z);      
      P(40).Init(x*y*z, zero, y*z*z);
      P(41).Init(zero, zero, z*z*z);
      
      P(42).Init(x*x*x*y, zero, zero);
      P(43).Init(x*x*x*z, zero, zero);
      P(44).Init(zero, x*y*y*z, -x*y*z*z);
      P(45).Init(zero, x*y*y*y, zero);
      P(46).Init(-x*x*y*z, zero, x*y*z*z);
      P(47).Init(zero, y*y*y*z, zero);
      P(48).Init(x*x*y*z, -x*y*y*z, zero);
      P(49).Init(zero, zero, x*z*z*z);
      P(50).Init(zero, zero, y*z*z*z);
      P(51).Init(x*x*y*z, x*y*y*z, x*y*z*z);
      P(52).Init(x*x*x*y*z, zero, zero);
      P(53).Init(zero, y*y*y*x*z, zero, zero);
      P(54).Init(zero, zero, z*z*z*x*y);*/
      
      /*nb_dof_loc = 3*(r+1)*(r+1)*(r+3);
      P.Reallocate(nb_dof_loc);
      int nb = 0;
      for (int i = 0; i <= r+2; i++)
        for (int j = 0; j <= r; j++)
          for (int k = 0; k <= r; k++)
            {
              P(nb).Init(Pow(x, i)*Pow(y, j)*Pow(z, k), zero, zero);
              P(nb+1).Init(zero, Pow(x, j)*Pow(y, i)*Pow(z, k), zero);
              P(nb+2).Init(zero, zero, Pow(x, j)*Pow(y, k)*Pow(z, i));
              nb += 3;
              } */
      
      nb_dof_loc = 40;
      nb_dof_loc = 21;
      P.Reallocate(nb_dof_loc);
      P(0).Init(one, zero, zero);
      P(1).Init(zero, one, zero);
      P(2).Init(zero, zero, one);
      P(3).Init(x, zero, zero);
      P(4).Init(y, zero, zero);
      P(5).Init(z, zero, zero);
      P(6).Init(zero, x, zero);
      P(7).Init(zero, y, zero);
      P(8).Init(zero, z, zero);
      P(9).Init(zero, zero, x);
      P(10).Init(zero, zero, y);
      P(11).Init(zero, zero, z);
      
      P(12).Init(x*y, zero, -y*z);
      P(13).Init(x*z, -y*z, zero);
      P(14).Init(zero, x*y, -x*z);
      P(15).Init(x*x, -x*y, -x*z);
      P(16).Init(-x*y, y*y, -y*z);
      P(17).Init(-x*z, -y*z, z*z);
      P(18).Init(x*x*y, -x*y*y, zero);
      P(19).Init(x*x*z, zero, -x*z*z);
      P(20).Init(zero, y*y*z, -y*z*z);
      
      /*P(12).Init(x*y, zero, zero);
      P(13).Init(x*z, zero, zero);
      P(14).Init(zero, x*y, zero);
      P(15).Init(zero, y*z, zero);
      P(16).Init(zero, zero, x*z);
      P(17).Init(zero, zero, y*z);
      
      P(18).Init(x*(one-x), zero, zero);
      P(19).Init(x*(one-x)*x, zero, zero);
      P(20).Init(x*(one-x)*y, zero, zero);
      P(21).Init(x*(one-x)*z, zero, zero);
      P(22).Init(x*(one-x)*x*y, zero, zero);
      P(23).Init(x*(one-x)*x*z, zero, zero);
      P(24).Init(x*(one-x)*x*y*z, zero, zero);

      P(25).Init(zero, y*(one-y), zero);
      P(26).Init(zero, y*(one-y)*x, zero);
      P(27).Init(zero, y*(one-y)*y, zero);
      P(28).Init(zero, y*(one-y)*z, zero);
      P(29).Init(zero, y*(one-y)*y*x, zero);
      P(30).Init(zero, y*(one-y)*y*z, zero);
      P(31).Init(zero, y*(one-y)*y*x*z, zero);
      
      P(32).Init(zero, zero, z*(one-z));
      P(33).Init(zero, zero, z*(one-z)*x);
      P(34).Init(zero, zero, z*(one-z)*y);
      P(35).Init(zero, zero, z*(one-z)*z);
      P(36).Init(zero, zero, z*(one-z)*z*x);
      P(37).Init(zero, zero, z*(one-z)*z*y);
      P(38).Init(zero, zero, z*(one-z)*z*y*x);
      
      P(39).Init(x*(one-x)*y*z, y*(one-y)*x*z, z*(one-z)*x*y);
      */
      
      Pdiv.Reallocate(nb_dof_loc);
      for (int i = 0; i < nb_dof_loc; i++)
        {
          DISP(i);
          DISP(P(i));
          GetDivPolynomial(P(i), Pdiv(i));
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
    void ComputeDivPhiRef(const R3& pt_loc, VectReal_wp& phi) const
    {
      phi.Reallocate(nb_dof_loc);
      for (int i = 0; i < nb_dof_loc; i++)
        phi(i) = Pdiv(i).Evaluate(pt_loc);
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
        psi(node)(0) = Pow(x, i+1)*Pow(y, j)*Pow(z, r-i-j-1);
        psi(node)(1) = Pow(x, i)*Pow(y, j+1)*Pow(z, r-i-j-1);
        psi(node)(2) = Pow(x, i)*Pow(y, j)*Pow(z, r-i-j);
        node++;
      }
  
  
}


int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  int r = 1; int r_test = r;
  
  HexahedronHdivFirstFamilyBis hex;
  hex.ConstructFiniteElement(r);

  Mesh<Dimension3> mesh;
  
  mesh.ReallocateVertices(8);
  mesh.ReallocateElements(1);
  mesh.Element(0).InitHexahedral(0, 1, 2, 3, 4, 5, 6, 7, 1);
  
  mesh.ReorientElements();

  VectR3 s; SetPoints<Dimension3> PointsElem;
  SetMatrices<Dimension3> MatricesElem;
  int nb_pts_quad = hex.GetNbPointsQuadratureInside();
  int nb_pts_dof = hex.GetNbPointsDof();
  int nb_dof_loc = hex.GetNbDof();
  VectR3 feval(nb_pts_dof);
  VectReal_wp res(nb_dof_loc);
  R3 vec_u, u_approche, u_exact;
  Real_wp div_exact, div_approche;
  VectR3 Points(nb_pts_quad);
  for (int i = 0; i < nb_pts_quad; i++)
    Points(i) = hex.PointsND(i);
  
  Matrix3_3 dfjm1; 
  // Real_wp jacob;
  Vector<TinyVector<MultivariatePolynomial<Real_wp>, 3> > Dk;
  Vector<MultivariatePolynomial<Real_wp> > divDk;
  GetDkBasis(r_test, Dk);
  
  divDk.Reallocate(Dk.GetM());
  for (int i = 0; i < Dk.GetM(); i++)
    GetDivPolynomial(Dk(i), divDk(i));
  
  R3 pt; Real_wp h = 1.0;
  // first test : we check if we can generate Dk with basis functions
  for (int num = 1; num <= Dk.GetM(); num++)
    {
      DISP(num);
      mesh.Vertex(0).Init(-0.9*h, -1.2*h, -0.4*h);
      mesh.Vertex(1).Init(1.2*h, -0.95*h, -0.2*h);
      mesh.Vertex(2).Init(1.1*h, 0.97*h, 0.1*h); 
      mesh.Vertex(3).Init(-0.87*h, 1.05*h, -0.3*h);
      
      mesh.Vertex(4).Init(-0.9*h, -1.2*h, 1.4*h);
      mesh.Vertex(5).Init(1.2*h, -0.95*h, 1.2*h);
      mesh.Vertex(6).Init(1.1*h, 0.97*h, 0.9*h); 
      mesh.Vertex(7).Init(-0.87*h, 1.05*h, 1.3*h);
      
      // boucle sur les elements
      for (int i = 0; i < mesh.GetNbElt(); i++)
        {
          // on calcule F_i(x), DF_i(x)
          mesh.GetVerticesElement(i, s);
          hex.FjElem(s, PointsElem, mesh, i);
          hex.DFjElem(s, PointsElem, MatricesElem, mesh, i);
          
          // on calcule P = J DF-1 p avec p polynome de Rk
          for (int k = 0; k < nb_pts_dof; k++)
            {
              pt = PointsElem.GetPointDof(k);
              vec_u(0) = Dk(num-1)(0).Evaluate(pt);
              vec_u(1) = Dk(num-1)(1).Evaluate(pt);
              vec_u(2) = Dk(num-1)(2).Evaluate(pt);
              GetInverse(MatricesElem.GetPointDof(k), dfjm1);
              Mlt(dfjm1, vec_u, feval(k));
              Mlt(Det(MatricesElem.GetPointDof(k)), feval(k));
            }
          
          // et on projette sur les fonctions de base de l'element de reference
          hex.ComputeProjectionDofRef(feval, res);
          
          // on calcule la norme L2 entre P = J DF-1 p et son projete
          // et la norme Hdiv aussi
          Real_wp err_Hdiv = 0, norme_Hdiv = 0;
          Real_wp erreur = 0, norme_L2 = 0; R3 evalP, eval_proj;
          Real_wp divP, div_proj;
          for (int k = 0; k < nb_pts_quad; k++)
            {
              eval_proj.Zero();
              for (int i = 0; i < nb_dof_loc; i++)
                eval_proj += res(i)*hex.Value_Phi(i, k);

              div_proj = 0;
              for (int i = 0; i < nb_dof_loc; i++)
                div_proj += res(i)*hex.Div_Phi(i, k);
              
              // evaluation de P
              pt = PointsElem.GetPointQuadrature(k);
              vec_u(0) = Dk(num-1)(0).Evaluate(pt);
              vec_u(1) = Dk(num-1)(1).Evaluate(pt);
              vec_u(2) = Dk(num-1)(2).Evaluate(pt);
              Real_wp jacob = Det(MatricesElem.GetPointDof(k));
              GetInverse(MatricesElem.GetPointDof(k), dfjm1);
              Mlt(dfjm1, vec_u, evalP);
              Mlt(jacob, evalP);
              
              // et de div P
              divP = divDk(num-1).Evaluate(pt);
              divP *= jacob;
              
              erreur += hex.WeightsND(k)*DotProd(evalP - eval_proj, evalP - eval_proj);
              err_Hdiv += hex.WeightsND(k)*square(divP - div_proj);
              norme_L2 += hex.WeightsND(k)*DotProd(evalP, evalP);
              norme_Hdiv += hex.WeightsND(k)*divP*divP;
            }
          
          erreur = sqrt(erreur);
          norme_L2 = sqrt(norme_L2);

          err_Hdiv = sqrt(err_Hdiv);
          norme_Hdiv = sqrt(norme_Hdiv);
          
          DISP(num); DISP(Dk(num-1)); DISP(erreur/norme_L2);
          if (abs(norme_Hdiv) > epsilon_machine)
            DISP(err_Hdiv/norme_Hdiv);
          else
            DISP(err_Hdiv);
          
        }
    }
  
  /* VectReal_wp errL2(2), errHdiv(2);
  for (int nh = 0; nh < 2; nh++)
    {
      h = 0.1/(nh+1);
      
      mesh.Vertex(0).Init(-0.9*h, -1.2*h, -0.4*h);
      mesh.Vertex(1).Init(1.2*h, -0.95*h, -0.2*h);
      mesh.Vertex(2).Init(1.1*h, 0.97*h, 0.1*h); 
      mesh.Vertex(3).Init(-0.87*h, 1.05*h, -0.3*h);
      
      mesh.Vertex(4).Init(-0.9*h, -1.2*h, 1.4*h);
      mesh.Vertex(5).Init(1.2*h, -0.95*h, 1.2*h);
      mesh.Vertex(6).Init(1.1*h, 0.97*h, 0.9*h); 
      mesh.Vertex(7).Init(-0.87*h, 1.05*h, 1.3*h);
      
      // boucle sur les elements
      for (int i = 0; i < mesh.GetNbElt(); i++)
        {
          // on calcule F_i(x), DF_i(x)
          mesh.GetVerticesElement(i, s);
          hex.FjElem(s, PointsElem, mesh, i);
          hex.DFjElem(s, PointsElem, MatricesElem, mesh, i);
          
          // on calcule P = J DF-1 f
          for (int k = 0; k < nb_pts_dof; k++)
            {
              pt = PointsElem.GetPointDof(k);
              vec_u = f(pt);
              GetInverse(MatricesElem.GetPointDof(k), dfjm1);
              Mlt(dfjm1, vec_u, feval(k));
              Mlt(Det(MatricesElem.GetPointDof(k)), feval(k));
            }
          
          // et on projette sur les fonctions de base de l'element de reference
          hex.ComputeProjectionDofRef(feval, res);
          
          // on calcule la norme L2 entre P = J DF-1 f et son projete
          // et la norme Hdiv aussi
          Real_wp err_Hdiv = 0, norme_Hdiv = 0;
          Real_wp erreur = 0, norme_L2 = 0; R3 evalP, eval_proj;
          Real_wp divP, div_proj;
          for (int k = 0; k < nb_pts_quad; k++)
            {
              eval_proj.Zero();
              for (int i = 0; i < nb_dof_loc; i++)
                eval_proj += res(i)*hex.Value_Phi(i, k);
              
              div_proj = 0;
              for (int i = 0; i < nb_dof_loc; i++)
                div_proj += res(i)*hex.Div_Phi(i, k);
              
              // evaluation de P
              pt = PointsElem.GetPointQuadrature(k);
              vec_u = f(pt);
              Real_wp jacob = Det(MatricesElem.GetPointDof(k));
              GetInverse(MatricesElem.GetPointDof(k), dfjm1);
              Mlt(dfjm1, vec_u, evalP);
              Mlt(jacob, evalP);
              
              // et de div P
              divP = fdiv(pt);
              divP *= jacob;
              
              erreur += hex.WeightsND(k)*DotProd(evalP - eval_proj, evalP - eval_proj);
              err_Hdiv += hex.WeightsND(k)*square(divP - div_proj);
              norme_L2 += hex.WeightsND(k)*DotProd(evalP, evalP);
              norme_Hdiv += hex.WeightsND(k)*divP*divP;
            }
          
          err_Hdiv = sqrt(err_Hdiv + erreur);
          norme_Hdiv = sqrt(norme_Hdiv + norme_L2);
 
          erreur = sqrt(erreur);
          norme_L2 = sqrt(norme_L2);
          
          errL2(nh) = abs(erreur/norme_L2);
          errHdiv(nh) = abs(err_Hdiv/norme_Hdiv);
          
        }
    }
  
  DISP(errL2); DISP(errHdiv);
  Real_wp rL2 = round((log(errL2(0)) - log(errL2(1)))/log(2.0));
  Real_wp rDiv = round((log(errHdiv(0)) - log(errHdiv(1)))/log(2.0));
  DISP(rL2);
  DISP(rDiv);
  */
  
  return FinalizeMontjoie();
}
