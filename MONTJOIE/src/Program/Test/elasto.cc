#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1
#define SELDON_WITH_PRECONDITIONING

#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"
#include "Hyperbolic/Elastic/MontjoieTimeElastic.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
#include "Elliptic/Helmholtz/Helmholtz1D.cxx"

using namespace Montjoie;

// definition de la source pour phi^s
template<class T>
class GaussianSourceP : public VirtualSourceField<T, Dimension2>
{
public :
  R2 origin;
  
  template<class GenericPb>
  GaussianSourceP(const GenericPb& var)
  {
    // centre de la source
    origin.Init(0.5, 0.8);
  }

  const R2& GetOrigin() const { return origin; }
  
  void EvaluateFunction(const R2& x, Vector<T>& f) const
  {
    Real_wp r = x.Distance(origin);
    // une bete gaussienne
    f(0) = exp(-10*r*r);
    
    // derivee en y d'une gaussienne pour eviter le mode 0
    //f(0) = -(x(1)-origin(1))*exp(-10.0*r*r);
  }
  
};

// definition de la source pour l'elastodynamique
template<class T>
class GaussianSourceElas : public VirtualSourceField<T, Dimension2>
{
public :
  R2 origin;

  template<class GenericPb>
  GaussianSourceElas(const GenericPb& var)
  {
    // centre de la source
    origin.Init(0.5, 0.8);
  }
  
  const R2& GetOrigin() const { return origin; }

  void EvaluateFunction(const R2& x, Vector<T>& f) const
  {
    Real_wp r = x.Distance(origin);
    Real_wp gaussian = exp(-10.0*r*r);
    // cas de la gaussienne
    f(1) = 20.0*(x(0)-origin(0))*gaussian;
    f(0) = -20.0*(x(1)-origin(1))*gaussian;
    
    // cas de la derivee en y d'une gaussienne pour eviter le mode 0
    //f(1) = -20.0*(x(0)-origin(0))*(x(1)-origin(1))*gaussian;
    //f(0) = (20.0*square(x(1)-origin(1)) - 1.0)*gaussian;
  }
  
};


// calcul du second membre pour phi et l'elasto 
// resultat dans rhs_s et rhs_elas
template<class TypeEqAcous, class TypeEqElas>
void ComputeRightHandSide(HyperbolicProblem<TypeEqAcous>& var_s, HyperbolicProblem<TypeEqAcous>& var_p,
			  HyperbolicProblem<TypeEqElas>& var_elas,
                          Vector<VectReal_wp>& rhs_pot, Vector<VectReal_wp>& rhs_elas, int Ndof)
{
  rhs_pot.Reallocate(1);
  rhs_elas.Reallocate(1);
  
  VolumetricSource<LaplaceEquation<Dimension2> > fS(var_s.var_harmonic);
  GaussianSourceP<Real_wp> fsrc(var_s.var_harmonic);
  fS.SetVolumeSourceFunction(fsrc);
  
  Vector<VectReal_wp> rhs_s(1);
  Vector<VirtualSourceFEM<Real_wp, Dimension2>* > f_vec(1);
  f_vec(0) = &fS;
  var_s.var_harmonic.ComputeGenericSource(rhs_s, f_vec, true);
  
  VolumetricSource<ElasticEquation<Dimension2> > fElas(var_elas.var_harmonic);
  GaussianSourceElas<Real_wp> fsrc_elas(var_elas.var_harmonic);
  fElas.SetVolumeSourceFunction(fsrc_elas);
  
  f_vec(0) = &fElas;
  var_elas.var_harmonic.ComputeGenericSource(rhs_elas, f_vec, true);
  
  // removing small values
  Real_wp threshold = epsilon_machine*Norm2(rhs_s(0));
  int Np = var_p.var_harmonic.GetNbDof();
  int Ns = var_s.var_harmonic.GetNbDof();
  rhs_pot(0).Reallocate(Ndof);
  rhs_pot(0).Fill(0);
  for (int i = 0; i < Ns; i++)
    {
      if (abs(rhs_s(0)(i)) <= threshold)
        rhs_pot(0)(Np+i) = 0;
      else
	rhs_pot(0)(Np+i) = rhs_s(0)(i);
    }
  
  // conversion to a sparse vector
  Copy(rhs_pot(0), var_s.GetSpaceSource());

  for (int i = 0; i < rhs_elas(0).GetM(); i++)
    if (abs(rhs_elas(0)(i)) <= threshold)
      rhs_elas(0)(i) = 0;
  
  Copy(rhs_elas(0), var_elas.GetSpaceSource());
}

// calcul de la rigidite pour Phi (que les termes interieurs)
// K = [K_{phi_p}, 0; 0, K_{phi_s}]
template<class TypeEqAcous, class MatSparse>
void GetStiffnessMatrix(EllipticProblem<TypeEqAcous>& var_p,
			EllipticProblem<TypeEqAcous>& var_s, MatSparse& mat_sp)
{
  int N = var_p.GetNbDof();
  GlobalGenericMatrix<Real_wp> nat_mat;
  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefDamping(0.0);
  nat_mat.SetCoefStiffness(1.0);
  
  // on calcule le laplacien pour phi^p
  var_p.AddMatrixFEM(mat_sp, nat_mat);
  
  // on rajoute le laplacien pour phi^s
  var_s.AddMatrixFEM(mat_sp, nat_mat, N, N);
}


// calcul de Ngamma = nombre de ddls sur la fontiere Neumann
// IndexDof(i) : numero local du ddl sur Gamma_neumann
template<class TypeEqAcous>
void GetDofOnGammaNeumann(EllipticProblem<TypeEqAcous>& var,
                          IVect& ListeDof, IVect& IndexDof, int& Ngamma)
{
  int N = var.mesh_num.GetNbDof();
  Ngamma = 0;
  IndexDof.Reallocate(N);
  IndexDof.Fill(-1);
  for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
    {
      int num_edge = i;
      int ref = var.mesh.BoundaryRef(i).GetReference();
      int cond = var.mesh.GetBoundaryCondition(ref);
      if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      	{
	  int num_elem = var.mesh.BoundaryRef(i).numElement(0);
	  int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(num_edge);
	  const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(num_elem);	  
	  for (int j = 0; j < Fb.GetNbDofBoundary(num_loc); j++)
	    {
	      int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
              int num_dof = var.mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
              if (IndexDof(num_dof) == -1)
                {
                  IndexDof(num_dof) = Ngamma;
                  Ngamma++;
                }
            }
        }
    }
  
  ListeDof.Reallocate(Ngamma);
  for (int i = 0; i < N; i++)
    if (IndexDof(i) >= 0)
      ListeDof(IndexDof(i)) = i;
}


// on recupere les ddls independants tels que
// \int_\Gamma phi = 0
// ListeCouple(i) = (a, b) tel que phi_a - coef phi_b est de moyenne nulle
// CoefCouple(i) coef used in the linear combination of phi_a, phi_b
template<class TypeEqAcous>
void GetIndependentDofs(EllipticProblem<TypeEqAcous>& var,
			Vector<TinyVector<int, 2> >& ListeCouple, 
			VectReal_wp& CoefCouple, IVect& removed_dof, int condition_gamma)
{
  int nb_edges_neu = 0;
  int r = var.mesh_num.GetOrder();
  for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
    {
      int ref = var.mesh.BoundaryRef(i).GetReference();
      int cond = var.mesh.GetBoundaryCondition(ref);
      if ((cond == condition_gamma))
	nb_edges_neu++;
    }
  
  if (nb_edges_neu <= 0)
    return;
  
  IVect ListeDofNeu(nb_edges_neu*(r+1));
  VectReal_wp MeanValueDof(nb_edges_neu*(r+1));
  int nb_dof_neu = 0;

  VectR2 s;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;
  
  // values of basis functions and gradient on the boundary of the elemnent
  Matrix<TinyVector<Real_wp, 1> > ValuePhi_Boundary;
  
  for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
    {
      int num_edge = i;
      int ref = var.mesh.BoundaryRef(i).GetReference();
      int cond = var.mesh.GetBoundaryCondition(ref);
      if ((cond == condition_gamma))
	{
	  int num_elem = var.mesh.BoundaryRef(i).numElement(0);
	  int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(num_edge);
	  const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(num_elem);

	  // we get the vertices of the element iquad
	  var.mesh.GetVerticesElement(num_elem, s);
    
	  // transformation F_i and DF_i are applied
	  // dans PointsElem il  n y a rien
	  // dans MatriceElem il y a DFj aux points de quadratures et aux points nodaux 
	  Fb.FjElem(s, PointsElem, var.mesh, num_elem);
	  Fb.DFjElem(s, PointsElem, MatricesElem, var.mesh, num_elem);
    
	  // restriction of DF_i on surface
	  Fb.FjSurfaceElem(s, PointsElem, var.mesh, num_elem, num_loc);
    
	  Fb.DFjSurfaceElem(s, PointsElem, MatricesElem,
			    var.mesh, num_elem, num_loc);

	  Fb.ComputeValuesPhiBoundary(ValuePhi_Boundary,
				      MatricesElem, var.mesh_num, num_elem, num_loc);
          
	  int nb_points_quadrature_edge = Fb.GetNbQuadBoundary(num_loc);
	  
	  for (int j = 0; j < Fb.GetNbDofBoundary(num_loc); j++)
	    {
	      int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
              int num_dof = var.mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
	      ListeDofNeu(nb_dof_neu) = num_dof;
	      MeanValueDof(nb_dof_neu) = 0;
	      for (int k = 0; k < nb_points_quadrature_edge; k++)
		{
		  Real_wp dsj = MatricesElem.GetDsQuadratureBoundary(k);
		  MeanValueDof(nb_dof_neu) += ValuePhi_Boundary(num_dof_loc, k)(0)*dsj*Fb.Weights1D(k);
		}
	      
	      nb_dof_neu++;
	    }
	}
    }
  
  // on assemble les dofs en ajoutant les moyennes
  Assemble(nb_dof_neu, ListeDofNeu, MeanValueDof);
  
  ListeDofNeu.Resize(nb_dof_neu);
  MeanValueDof.Resize(nb_dof_neu);
  
  //DISP(ListeDofNeu); DISP(MeanValueDof);
  
  ListeCouple.Reallocate(nb_dof_neu-1);
  CoefCouple.Reallocate(nb_dof_neu-1);
  removed_dof.Reallocate(1);
  removed_dof(0) = ListeDofNeu(0);
  for (int i = 1; i < nb_dof_neu; i++)
    {
      ListeCouple(i-1)(0) = ListeDofNeu(i);
      ListeCouple(i-1)(1) = ListeDofNeu(0);
      CoefCouple(i-1) = MeanValueDof(i)/MeanValueDof(0);      
    }
  
  //DISP(ListeCouple); DISP(CoefCouple); DISP(removed_dof);
  
  /*
  typedef typename TypeElement::TypeElement2 TypeElement2;
  bool test_loop = true;
  Vector<bool> EdgeFree(var.mesh.GetNbBoundaryRef());
  EdgeFree.Fill(true);
  IVect NbEdgesAroundVertex(var.mesh.GetNbVertices());
  IVect EdgeNumber(var.mesh.GetNbVertices());
  IVect ref_neu(var.mesh.GetNbReferences()+1);
  ref_neu.Fill(0);
  for (int i = 0; i <= var.mesh.GetNbReferences(); i++)
    if (var.mesh.GetBoundaryCondition(i) == condition_gamma)
      ref_neu(i) = 1;
  
  int r = var.mesh_num.GetOrder();
  const TypeElement2& Fb = var.GetFaceBasis2(r);
  int nb_edges_neu = 0;
  for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
    {
      int ref = var.mesh.BoundaryRef(i).GetReference();
      int cond = var.mesh.GetBoundaryCondition(ref);
      if ((cond == condition_gamma))
	nb_edges_neu++;
    }
  
  ListeCouple.Reallocate(nb_edges_neu*(r+1));
  CoefCouple.Reallocate(nb_edges_neu*(r+1));
  int nb_dof_neu = 0, nb_removed_dof = 0;
  removed_dof.Reallocate(nb_edges_neu*(r+1));
  while (test_loop)
    {
      NbEdgesAroundVertex.Fill(0);
      for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
	{
	  int num_edge = i;
	  int ref = var.mesh.BoundaryRef(i).GetReference();
	  int cond = var.mesh.GetBoundaryCondition(ref);
	  if ((cond == condition_gamma) && EdgeFree(i))
	    {
	      int n1 = var.mesh.BoundaryRef(i).numVertex(0);
	      int n2 = var.mesh.BoundaryRef(i).numVertex(1);
	      NbEdgesAroundVertex(n1)++;
	      NbEdgesAroundVertex(n2)++;
	      EdgeNumber(n1) = num_edge;
	      EdgeNumber(n2) = num_edge;
	    }
	}
      
      int n_init = -1;
      for (int i = 0; i < var.mesh.GetNbVertices(); i++)
	{
	  if (NbEdgesAroundVertex(i) == 1)
	    n_init = i;
	  else if ((NbEdgesAroundVertex(i) == 2) && (n_init == -1))
	    {
	      n_init = i;
	    }
	}
      
      if (n_init >= 0)
	{
	  int ne = EdgeNumber(n_init);
	  int n1 = n_init, n2 = -1;
	  if (var.mesh.Boundary(ne).numVertex(0) == n_init)
	    n2 = var.mesh.Boundary(ne).numVertex(1);
	  else
	    n2 = var.mesh.Boundary(ne).numVertex(0);
	  
	  EdgeFree(ne) = false;
          int last_dof_vertex = -1;
	  while (n2 != n_init)
	    {
	      int ne_next = -1;
	      int n3 = var.mesh.FindFollowingVertex(n1, n2, 1, ref_neu, ne_next);
	      //DISP(n1); DISP(n2); DISP(n3);
	      if (n3 != -1)
		{
		  int num_elem = var.mesh.Boundary(ne).numElement(0);
		  int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(ne);

		  int num_elem2 = var.mesh.Boundary(ne_next).numElement(0);
		  int num_loc2 = var.mesh.Element(num_elem2).GetPositionBoundary(ne_next);

                  int nv1 = -1, nv2 = -1;
                  for (int k = 0; k < 4; k++)
                    {
                      if (var.mesh.Element(num_elem).numVertex(k) == n1)
                        nv1 = k;
                      
                      if (var.mesh.Element(num_elem).numVertex(k) == n2)
                        nv2 = k;                      
                    }
                  
		  ListeCouple(nb_dof_neu)(0) = var.mesh_num.Element(num_elem).GetNumberDof(nv1);
		  ListeCouple(nb_dof_neu)(1) = var.mesh_num.Element(num_elem).GetNumberDof(nv2);
                  last_dof_vertex = ListeCouple(nb_dof_neu)(1);
		  nb_dof_neu++;
		  		  
		  bool same_orientation = true;
		  if ( (var.mesh.Element(num_elem).numVertex(num_loc) == var.mesh.Element(num_elem2).numVertex(num_loc2))
		       || (var.mesh.Element(num_elem).numVertex((num_loc+1)%4) == var.mesh.Element(num_elem2).numVertex((num_loc2+1)%4)) )
		    same_orientation = false;
		  
		  for (int k = 1; k < r; k++)
		    {
		      int local_dof = Fb.GetLocalNumber(num_loc, k);
		      int num_dof = var.mesh_num.Element(num_elem).GetNumberDof(local_dof);
		      ListeCouple(nb_dof_neu)(0) = num_dof;
		      
		      int local_dof2 = Fb.GetLocalNumber(num_loc2, k);
		      if (!same_orientation)
			local_dof2 = Fb.GetLocalNumber(num_loc2, r-k);
		      
		      int num_dof2 = var.mesh_num.Element(num_elem2).GetNumberDof(local_dof2);
		      ListeCouple(nb_dof_neu)(1) = num_dof2;
		      nb_dof_neu++;
		    }
		  
		  EdgeFree(ne_next) = false;

		  n1 = n2;
		  n2 = n3;
		  ne = ne_next;
		}
	      else
		n2 = n_init;
	    }
	  
	  if (last_dof_vertex >= 0)
            {
              removed_dof(nb_removed_dof) = last_dof_vertex;
              nb_removed_dof++;
            }
	}
      else
	test_loop = false;
    }
  
  ListeCouple.Resize(nb_dof_neu);
  removed_dof.Resize(nb_removed_dof);
  */
}


template<class TypeEqAcous>
void GetCoordinateDofs(EllipticProblem<TypeEqAcous>& var,
		       VectR2& CoorDofs)
{
  CoorDofs.Reallocate(var.mesh_num.GetNbDof());
  VectR2 points_elt;
  for (int i = 0; i < var.mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(i);
      var.mesh.GetDofsElement(i, points_elt, Fb.GetGeometricElement());
      for (int j = 0; j < points_elt.GetM(); j++)
	CoorDofs(var.mesh_num.Element(i).GetNumberDof(j)) = points_elt(j);
    }
}


// calcul de la matrice de masse : Mh d^2 Phi/dt^2
template<class TypeEqAcous, class MatSparse>
void GetMassMatrix(EllipticProblem<TypeEqAcous>& var_p,
		   EllipticProblem<TypeEqAcous>& var_s,
		   MatSparse& Mh)
{
  int N = var_p.GetNbDof();
  GlobalGenericMatrix<Real_wp> nat_mat;
  nat_mat.SetCoefMass(1.0);
  nat_mat.SetCoefDamping(0.0);
  nat_mat.SetCoefStiffness(0.0);
  
  // on calcule la matrice de masse pour phi^p
  var_p.AddMatrixFEM(Mh, nat_mat);
  
  // et pour phi^s
  var_s.AddMatrixFEM(Mh, nat_mat, N, N);
}


// rajout de alpha \int n \times grad phi_i n \times grad phi_j
template<class TypeEqAcous, class TypeEqElas, class MatSparse>
void AddMatrixLaplacian(const Real_wp& alpha, const Real_wp& beta, const Real_wp& gamma,
			EllipticProblem<TypeEqAcous>& var,
			EllipticProblem<TypeEqElas>& var_e,
			MatSparse& mat_sp, const IVect& IndexDofRow, const IVect& IndexDofCol,
			const Vector<bool>& condition_taken,
			int offset_row, int offset_col, bool mlt_mu_sur_rho)
{
  VectReal_wp Ones, grad_phi, grad_feval, contrib;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;
  VectR2 s;
  for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
    {
      int num_edge = i;
      int ref = var.mesh.BoundaryRef(i).GetReference();
      int cond = var_e.mesh.GetBoundaryCondition(ref);
          
      if (condition_taken(cond))
	{
	  int num_elem = var.mesh.BoundaryRef(i).numElement(0);
	  int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(num_edge);
	  const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(num_elem);	  
	  int nb_quad = Fb.GetNbQuadBoundary(num_loc);
	  R2 normale = Fb.NormaleLoc(num_loc);

          // calcul de mu/rho
          int ref = var.mesh.Element(num_elem).GetReference();
          Real_wp mu_rho = 1.0/var.ref_rho(ref).GetConstant();

          // restriction of DF_i on surface
          var.mesh.GetVerticesElement(num_elem, s);
          Fb.FjElem(s, PointsElem, var.mesh, num_elem);
          Fb.DFjElem(s, PointsElem, MatricesElem, var.mesh, num_elem);
	  
          Fb.FjSurfaceElem(s, PointsElem, var.mesh, num_elem, num_loc);
          
          Fb.DFjSurfaceElem(s, PointsElem, MatricesElem,
                            var.mesh, num_elem, num_loc);

	  Ones.Reallocate(Fb.GetNbDof());
	  grad_phi.Reallocate(2*nb_quad);
	  grad_feval.Reallocate(2*nb_quad);
	  contrib.Reallocate(Fb.GetNbDof());
	  for (int j = 0; j < Fb.GetNbDofBoundary(num_loc); j++)
	    {
	      int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
	      Ones.Fill(0); Ones(num_dof_loc) = 1.0;
	      Fb.ApplyNablaShTranspose(num_loc, Ones, grad_phi);
	      for (int k = 0; k < nb_quad; k++)
		{
		  Real_wp dsj = MatricesElem.GetDsQuadratureBoundary(k);
		  R2 normale_fj = MatricesElem.GetNormaleQuadratureBoundary(k);
		  Real_wp n_grad_phi = normale(0)*grad_phi(2*k+1) - normale(1)*grad_phi(2*k);
		  Real_wp coef = alpha + beta*normale_fj(0) + gamma*normale_fj(1);
		  grad_feval(2*k) = -n_grad_phi*normale(1)*Fb.WeightsQuadratureBoundary(k, num_loc)*coef/dsj;
		  grad_feval(2*k+1) = n_grad_phi*normale(0)*Fb.WeightsQuadratureBoundary(k, num_loc)*coef/dsj;
		}
	      
	      Fb.ComputeIntegralSurfaceGradientRef(grad_feval, contrib, num_loc);
	      for (int k = 0; k < Fb.GetNbDof(); k++)
		if (abs(contrib(k)) > 1e-12)
		  {
		    int num_dof_p = var.mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
		    int num_dof_s = var.mesh_num.Element(num_elem).GetNumberDof(k);
		    Real_wp val = contrib(k);
		    if (mlt_mu_sur_rho)
		      val *= mu_rho;
		    
		    mat_sp.AddInteraction(offset_row + IndexDofRow(num_dof_s),
					  offset_col + IndexDofCol(num_dof_p), val);
		  }
	    }
	}
    }
}

// rajout de alpha \int n \times grad phi_i phi_j (ou la transposee)
template<class TypeEqAcous, class TypeEqElas, class MatSparse>
void AddMatrixRpsi(const Real_wp& alpha, const Real_wp& beta, const Real_wp& gamma,
		   EllipticProblem<TypeEqAcous>& var_p,
		   EllipticProblem<TypeEqAcous>& var_s,
                   EllipticProblem<TypeEqElas>& var_e,
		   MatSparse& mat_sp, const IVect& IndexDofRow, const IVect& IndexDofCol,
                   const Vector<bool>& condition_taken,
		   int offset_row, int offset_col, bool transpose, bool mlt_mu_sur_rho)
{
  VectReal_wp Ones, val_phi, grad_feval, contrib;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;
  VectR2 s;
  for (int i = 0; i < var_p.mesh.GetNbBoundaryRef(); i++)
    {
      int num_edge = i;
      int ref = var_p.mesh.BoundaryRef(i).GetReference();
      int cond = var_e.mesh.GetBoundaryCondition(ref);
          
      if (condition_taken(cond))
	{
	  int num_elem = var_p.mesh.BoundaryRef(i).numElement(0);
	  int num_loc = var_p.mesh.Element(num_elem).GetPositionBoundary(num_edge);
	  const ElementReference<Dimension2, 1>& Fb_p = var_p.GetReferenceElement(num_elem);	  
	  const ElementReference<Dimension2, 1>& Fb_s = var_s.GetReferenceElement(num_elem); 
	  int nb_quad = Fb_s.GetNbQuadBoundary(num_loc);
	  R2 normale = Fb_s.NormaleLoc(num_loc);

          // calcul de mu/rho
          int ref = var_s.mesh.Element(num_elem).GetReference();
          Real_wp mu_rho = 1.0/var_s.ref_rho(ref).GetConstant();
	  
          // restriction of DF_i on surface
          var_s.mesh.GetVerticesElement(num_elem, s);
          Fb_s.FjElem(s, PointsElem, var_s.mesh, num_elem);
          Fb_s.DFjElem(s, PointsElem, MatricesElem, var_s.mesh, num_elem);
	  
          Fb_s.FjSurfaceElem(s, PointsElem, var_s.mesh, num_elem, num_loc);
          
          Fb_s.DFjSurfaceElem(s, PointsElem, MatricesElem,
			      var_s.mesh, num_elem, num_loc);
	  
	  Ones.Reallocate(Fb_p.GetNbDof());
	  val_phi.Reallocate(nb_quad);
	  grad_feval.Reallocate(2*nb_quad);
	  contrib.Reallocate(Fb_s.GetNbDof());
	  
	  // boucle sur les colonnes phi_j^p
	  for (int j = 0; j < Fb_p.GetNbDofBoundary(num_loc); j++)
	    {
	      int num_dof_loc = Fb_p.GetLocalNumber(num_loc, j);
	      Ones.Fill(0); Ones(num_dof_loc) = 1.0;
	      
	      //DISP(j); DISP(Ones);
	      // calcul sur les points de quadrature de phi^s
	      Fb_p.ApplyShTranspose(num_loc, Ones, val_phi, Fb_s.GetOrder());
	      //DISP(val_phi);

	      // boucle sur ces points de quadrature
	      for (int k = 0; k < nb_quad; k++)
		{
		  R2 normale_fj = MatricesElem.GetNormaleQuadratureBoundary(k);
		  Real_wp coef = alpha + beta*normale_fj(0) + gamma*normale_fj(1);
		  grad_feval(2*k) = -val_phi(k)*normale(1)*Fb_s.WeightsQuadratureBoundary(k, num_loc)*coef;
		  grad_feval(2*k+1) = val_phi(k)*normale(0)*Fb_s.WeightsQuadratureBoundary(k, num_loc)*coef;
		}
	      
	      //DISP(grad_feval);
	      
	      Fb_s.ComputeIntegralSurfaceGradientRef(grad_feval, contrib, num_loc);
	      //DISP(num_loc);
	      for (int k = 0; k < Fb_s.GetNbDof(); k++)
		if (abs(contrib(k)) > 1e-12)
		  {
		    int num_dof_p = var_p.mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
		    int num_dof_s = var_s.mesh_num.Element(num_elem).GetNumberDof(k);
		    Real_wp val = contrib(k);
		    //DISP(num_dof_p); DISP(num_dof_s); DISP(val);
		    if (mlt_mu_sur_rho)
		      val *= mu_rho;
		    
		    if (transpose)
		      mat_sp.AddInteraction(offset_row + IndexDofRow(num_dof_p),
					    offset_col + IndexDofCol(num_dof_s), val);
		    else
		      mat_sp.AddInteraction(offset_row + IndexDofRow(num_dof_s),
					    offset_col + IndexDofCol(num_dof_p), val);
		  }
	    }
	}
    }
}


// rajout de alpha \int (alpha + beta n1 + gamma n2) phi_i phi_j
template<class TypeEqAcous, class TypeEqElas, class MatSparse>
void AddMatrixCh(const Real_wp& alpha, const Real_wp& beta, const Real_wp& gamma,
		 EllipticProblem<TypeEqAcous>& var,
                 EllipticProblem<TypeEqElas>& var_e,
                 MatSparse& mat_sp, const IVect& IndexDofRow,
		 const IVect& IndexDofCol, const Vector<bool>& condition_taken,
		 int offset_row, int offset_col)
{
  VectReal_wp Ones, val_phi, feval, contrib;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;
  VectR2 s;
  for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
    {
      int num_edge = i;
      int ref = var.mesh.BoundaryRef(i).GetReference();
      int cond = var_e.mesh.GetBoundaryCondition(ref);
      if (condition_taken(cond))
      	{
	  int num_elem = var.mesh.BoundaryRef(i).numElement(0);
	  int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(num_edge);
	  const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(num_elem);
	  int nb_quad = Fb.GetNbQuadBoundary(num_loc);
          var.mesh.GetVerticesElement(num_elem, s);
          Fb.FjElem(s, PointsElem, var.mesh, num_elem);
          Fb.DFjElem(s, PointsElem, MatricesElem, var.mesh, num_elem);
          
          // restriction of DF_i on surface
          Fb.FjSurfaceElem(s, PointsElem, var.mesh, num_elem, num_loc);
          
          Fb.DFjSurfaceElem(s, PointsElem, MatricesElem,
                            var.mesh, num_elem, num_loc);
	  
          Ones.Reallocate(Fb.GetNbDof());
	  val_phi.Reallocate(nb_quad);
	  feval.Reallocate(nb_quad);
	  contrib.Reallocate(Fb.GetNbDof());
	  for (int j = 0; j < Fb.GetNbDofBoundary(num_loc); j++)
	    {
	      int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
	      Ones.Fill(0); Ones(num_dof_loc) = 1.0;
	      Fb.ApplyShTranspose(num_loc, Ones, val_phi);
	      for (int k = 0; k < nb_quad; k++)
                {
                  R2 normale = MatricesElem.GetNormaleQuadratureBoundary(k);
                  Real_wp dsj = MatricesElem.GetDsQuadratureBoundary(k);
		  feval(k) = val_phi(k)*Fb.WeightsQuadratureBoundary(k, num_loc)*dsj*
		    (alpha + beta*normale(0) + gamma*normale(1));
                }
              
	      Fb.ComputeIntegralSurfaceRef(feval, contrib, num_loc);
	      for (int k = 0; k < Fb.GetNbDof(); k++)
		if ( abs(contrib(k)) > 1e-12)
		  {
		    int icol = var.mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
		    int irow = var.mesh_num.Element(num_elem).GetNumberDof(k);
		    
		    mat_sp.AddInteraction(offset_row + IndexDofRow(irow),
					  offset_col + IndexDofCol(icol), contrib(k));
		  }
	    }
	}
    }
}


// rajout de \int (alpha + beta n1 + gamma n2) phi_j (ou transposee)
template<class TypeEqAcous, class TypeEqElas, class MatSparse>
void AddMatrixPhi(const Real_wp& alpha, const Real_wp& beta, const Real_wp& gamma,
		 EllipticProblem<TypeEqAcous>& var,
                 EllipticProblem<TypeEqElas>& var_e,
                 MatSparse& mat_sp, const IVect& IndexDof, const Vector<bool>& condition_taken,
		  int offset_row, int offset_col, bool transpose)
{
  VectReal_wp feval, contrib;
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;
  VectR2 s;
  for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
    {
      int num_edge = i;
      int ref = var.mesh.BoundaryRef(i).GetReference();
      int cond = var_e.mesh.GetBoundaryCondition(ref);
      if (condition_taken(cond))
      	{
	  int num_elem = var.mesh.BoundaryRef(i).numElement(0);
	  int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(num_edge);
	  const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(num_elem);	  
	  int nb_quad = Fb.GetNbQuadBoundary(num_loc);
          var.mesh.GetVerticesElement(num_elem, s);
          Fb.FjElem(s, PointsElem, var.mesh, num_elem);
          Fb.DFjElem(s, PointsElem, MatricesElem, var.mesh, num_elem);
          
          // restriction of DF_i on surface
          Fb.FjSurfaceElem(s, PointsElem, var.mesh, num_elem, num_loc);
          
          Fb.DFjSurfaceElem(s, PointsElem, MatricesElem,
                            var.mesh, num_elem, num_loc);
	  
	  feval.Reallocate(nb_quad);
	  contrib.Reallocate(Fb.GetNbDof());
	  for (int k = 0; k < nb_quad; k++)
	    {
	      R2 normale = MatricesElem.GetNormaleQuadratureBoundary(k);
	      Real_wp dsj = MatricesElem.GetDsQuadratureBoundary(k);
	      feval(k) = Fb.WeightsQuadratureBoundary(k, num_loc)*dsj*
		(alpha + beta*normale(0) + gamma*normale(1));
	    }
              
	  Fb.ComputeIntegralSurfaceRef(feval, contrib, num_loc);
	  for (int k = 0; k < Fb.GetNbDof(); k++)
	    if ( abs(contrib(k)) > 1e-12)
	      {
		int icol = var.mesh_num.Element(num_elem).GetNumberDof(k);
		
		if (transpose)
		  mat_sp.AddInteraction(offset_row + IndexDof(icol),
					offset_col, contrib(k));
		else
		  mat_sp.AddInteraction(offset_row,
					offset_col + IndexDof(icol), contrib(k));
	      }
	}
    }
}


// rajout de alpha \int_0^1 \int_0^y phi_j(\xi) d\xi phi_i(y) dy (ou transposee)
template<class TypeEqAcous, class TypeEqElas, class MatSparse>
void AddMatrixIntPhi(const Real_wp& alpha,
		     EllipticProblem<TypeEqAcous>& var,
		     EllipticProblem<TypeEqElas>& var_e,
		     MatSparse& mat_sp, const Vector<bool>& condition_taken,
		     int offset_row, int offset_col)
{
  // d'abord on calcule omega ds_j pour tous les points de quadrature
  int r = var.mesh_num.GetOrder();
  const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(0);
  Matrix<Real_wp> Poids(var.mesh.GetNbBoundaryRef(), r+1);
  Matrix<R2> Points(var.mesh.GetNbBoundaryRef(), r+1);
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;
  VectR2 s; VectReal_wp pos_edge(var.mesh.GetNbBoundaryRef());
  for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
    {
      int num_edge = i;
      int ref = var.mesh.BoundaryRef(i).GetReference();
      int cond = var_e.mesh.GetBoundaryCondition(ref);
      int n0 = var.mesh.BoundaryRef(i).numVertex(0);
      int n1 = var.mesh.BoundaryRef(i).numVertex(1);
      Real_wp y0 = var.mesh.Vertex(n0)(1);
      Real_wp y1 = var.mesh.Vertex(n1)(1);
      pos_edge(i) = 0.5*(y0 + y1);
      if (condition_taken(cond))
      	{
	  int num_elem = var.mesh.BoundaryRef(i).numElement(0);
	  int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(num_edge);
          var.mesh.GetVerticesElement(num_elem, s);
          Fb.FjElem(s, PointsElem, var.mesh, num_elem);
          Fb.DFjElem(s, PointsElem, MatricesElem, var.mesh, num_elem);
          
          // restriction of DF_i on surface
          Fb.FjSurfaceElem(s, PointsElem, var.mesh, num_elem, num_loc);
          
          Fb.DFjSurfaceElem(s, PointsElem, MatricesElem,
                            var.mesh, num_elem, num_loc);
	  
	  for (int k = 0; k <= r; k++)
	    {
	      Real_wp dsj = MatricesElem.GetDsQuadratureBoundary(k);
	      Points(i, k) = PointsElem.GetPointQuadratureBoundary(k);
	      Poids(i, k) = dsj*Fb.WeightsQuadratureBoundary(k, num_loc);
	    }
	}
    }

  IVect permut_edge(var.mesh.GetNbBoundaryRef());
  permut_edge.Fill();
  Sort(pos_edge, permut_edge);
  //DISP(pos_edge); DISP(permut_edge);

  // valeurs de phi aux points de quadrature
  const ElementReference<Dimension1, 1>& edge = Fb.GetSurfaceFiniteElement(0);
  Matrix<Real_wp> ValPhi(r+1, r+1);
  Vector<Matrix<Real_wp> > ValPhiInt(r+1);
  VectReal_wp phi(r+1); //DISP(edge.Points);
  for (int k = 0; k <= r; k++)
    {
      edge.ComputeValuesPhiRef(edge.Points(k), phi);
      for (int j = 0; j <= r; j++)
	ValPhi(j, k) = phi(j);
    }
  
  for (int m = 0; m <= r; m++)
    {
      ValPhiInt(m).Reallocate(r+1, r+1);
      for (int k = 0; k <= r; k++)
	{
	  edge.ComputeValuesPhiRef(edge.Points(m)*edge.Points(k), phi);
	  for (int j = 0; j <= r; j++)
	    ValPhiInt(m)(j, k) = phi(j);
	}
      //DISP(m); DISP(ValPhiInt(m));
    }
  
  // integrale des fonctions de base (incrémentées au fur et a mesure)
  VectReal_wp int_basis(var.mesh_num.GetNbDof());
  VectReal_wp int_basis_current(var.mesh_num.GetNbDof());
  int_basis.Fill(0);
  int_basis_current.Fill(0);
  
  // boucle pour l'integrale exterieure
  for (int ip = 0; ip < var.mesh.GetNbBoundaryRef(); ip++)
    {
      int i = permut_edge(ip);
      int num_edge = i;
      int ref = var.mesh.BoundaryRef(i).GetReference();
      int cond = var_e.mesh.GetBoundaryCondition(ref);
      if (condition_taken(cond))
	{
	  int num_elem = var.mesh.BoundaryRef(i).numElement(0);
	  int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(num_edge);
	  int n0 = var.mesh.Element(num_elem).numVertex(num_loc);
	  int n1 = var.mesh.Element(num_elem).numVertex((num_loc+1)%4);
	  Real_wp y0 = var.mesh.Vertex(n0)(1);
	  Real_wp y1 = var.mesh.Vertex(n1)(1);
	  //DISP(y0); DISP(y1);
	  if (y0 > y1)
	    {
	      cout << " arete mal oriente" << endl;
	      abort();
	    }

	  // on boucle sur les points de quadrature pour y
	  for (int k = 0; k <= r; k++)
	    {
	      //Real_wp y = Points(i, k)(1);
	      
	      // premiere etape on calcule la primitive de toutes
	      // les fonctions de base de 0 a y
	      // on part de la valeur jusqu'a l'intervalle precedent
	      for (int j = 0; j < var.mesh_num.GetNbDof(); j++)
		int_basis_current(j) = int_basis(j);
	      
	      // ensuite on modifie que les fonctions de base de l'intervalle
	      // qui auront une contribution non-nulle
	      for (int jloc = 0; jloc <= r; jloc++)
		{
		  int j2 = Fb.GetLocalNumber(num_loc, jloc);
		  int j = var.mesh_num.Element(num_elem).GetNumberDof(j2);
		  Real_wp vloc = 0.0;
		  for (int k2 = 0; k2 <= r; k2++)
		    vloc += ValPhiInt(k)(jloc, k2)*edge.Weights(k2);
		  
		  int_basis_current(j) += vloc*edge.Points(k)*abs(y1-y0);
		}
	      
	      //DISP(y); DISP(int_basis_current);
	      
	      // ensuite on integre contre la fonction phi_i
	      for (int irow = 0; irow <= r; irow++)
		{
		  int num_row = Fb.GetLocalNumber(num_loc, irow);
		  int num_row_glob = var.mesh_num.Element(num_elem).GetNumberDof(num_row);
		  Real_wp vloc = alpha*ValPhi(irow, k)*Poids(i, k);
		  if (abs(vloc) > 1e-14)
		    for (int j = 0; j < var.mesh_num.GetNbDof(); j++)
		      if (abs(int_basis_current(j)) > 1e-14)
			mat_sp
			  .AddInteraction(offset_row + num_row_glob,
					  offset_col + j,
					  vloc*int_basis_current(j));
		}	      
	    }
	  
	  // ensuite on incremente les primitives pour l'intervalle d'apres
	  for (int jloc = 0; jloc <= r; jloc++)
	    {
	      int j2 = Fb.GetLocalNumber(num_loc, jloc);
	      int j = var.mesh_num.Element(num_elem).GetNumberDof(j2);
	      Real_wp vloc = 0.0;
	      for (int k2 = 0; k2 <= r; k2++)
		vloc += ValPhi(jloc, k2)*edge.Weights(k2);
	      
	      int_basis(j) += vloc*abs(y1-y0);
	    }	    
	}
    }
}


template<class TypeEqAcous, class TypeEqElas, class MatSparse>
void AddMatrixDtN(EllipticProblem<TypeEqAcous>& var,
                  EllipticProblem<TypeEqElas>& var_e,
                  MatSparse& Sh, MatSparse& Kh, IVect& ListeDofNeu,
                  IVect& IndexDofNeumann, int Ngamma)
{
  // on recupere les ddls avec Dirichlet avant calcul du DtN
  IVect DirichletElas = var_e.GetDirichletDofNumber();
  int N = DirichletElas.GetM();
  
  // on impose Dirichlet sur Gamma_neumann pour calculer le DtN
  IVect NewDirichletElas(DirichletElas.GetM() + 2*Ngamma);
  for (int i = 0; i < N; i++)
    NewDirichletElas(i) = DirichletElas(i);
  
  for (int i = 0; i < Ngamma; i++)
    {
      NewDirichletElas(N+i) = ListeDofNeu(i);
      NewDirichletElas(N+Ngamma+i) = var_e.mesh_num.GetNbDof() + ListeDofNeu(i);
    }
  
  var_e.SetDirichletDofs(NewDirichletElas.GetM(), NewDirichletElas);
  
  // on calcule que la partie rigidite
  GlobalGenericMatrix<Real_wp> nat_mat;
  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefStiffness(1.0);
  
  All_LinearSolver solver_elas(var_e);
  solver_elas.SetDirectSolver();
  
  // bord Neumann
  IVect ref_cond(var_e.mesh.GetNbReferences()+1);
  ref_cond.Fill(0);
  for (int i = 1; i <= var_e.mesh.GetNbReferences(); i++)
    if (var_e.mesh.GetBoundaryCondition(i) == BoundaryConditionEnum::LINE_NEUMANN)
      ref_cond(i) = 1;

  int offset_gamma = var.mesh_num.GetNbDof()*2;
  
  // objet pour calculer C epsilon(u) n sur le bord Gamma
  MeshInterpolationFEM<Dimension2> interp_gamma;
  int r = var_e.mesh_num.GetOrder();
  Mesh<Dimension2> mesh_subdiv;
  interp_gamma.SetGaussQuadrature(r);
  interp_gamma.InitProjectionSurface(var_e.mesh);
  interp_gamma.ComputeSurfaceMesh(ref_cond, var_e.mesh, mesh_subdiv, var_e);
  
#ifdef SELDON_WITH_MPI
  interp_gamma.GatherQuadraturePoints(MPI_COMM_WORLD);
#else
  interp_gamma.GatherQuadraturePoints();
#endif
  
  // on calcule et factorise la matrice elements finis
  solver_elas.PerformFactorizationStep(nat_mat);
  
  // element fini
  const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(0);
  
  VectReal_wp fevalX(r+1), fevalY(r+1), contribX(Fb.GetNbDof()), contribY(Fb.GetNbDof());
  contribX.Fill(0); contribY.Fill(0);
  
  VectReal_wp feval_grad(2*Fb.GetNbPointsQuadratureInside()), feval_rot(2*Fb.GetNbPointsQuadratureInside());
  
  VectR2 trace_u, trace_du_dn;
  Vector<VectR2> Uquad;
  Vector<Vector<TinyVector<Real_wp, 4> > > gradU_quad;
  VectReal_wp rhs_elas(var_e.GetNbDof()), sol_elas(var_e.GetNbDof());
  Vector<VectReal_wp> mode_dtn(2*Ngamma);
  sol_elas.Fill(0);
  for (int i = 0; i < 2*Ngamma; i++)
    {
      rhs_elas.Fill(0);
      rhs_elas(NewDirichletElas(N+i)) = 1.0;
      
      // calcul de la solution
      sol_elas = rhs_elas;
      solver_elas.ComputeSolution(sol_elas, nat_mat);
      mode_dtn(i) = sol_elas;
      
      // calcul de C epsilon(u) n sur Gamma
      //interp_gamma.ComputeEnHnOnBoundary(var_e, sol_elas, trace_u, trace_du_dn);
      abort();
      
      // calcul de la solution sur les points de quadrature
      //var_e.ComputeNodalUgradU(sol_elas, Uquad, gradU_quad, true, false);
      abort();
      
      // on calcule la partie \int T u_r v 
      int offset = 0; Real_wp poids;
      for (int n = 0; n < interp_gamma.GetNbBoundary(); n++)
        {
          for (int k = 0; k < interp_gamma.GetNbPointsQuadrature(n); k++)
            {
              //DISP(n); DISP(k);
              //DISP(interp_gamma.PointsQuadrature(n)(k));
              //DISP(trace_u(offset+k)); 
              //DISP(trace_du_dn(offset+k));
              poids = interp_gamma.WeightsQuadrature(n)(k); //DISP(poids);
              fevalX(k) = trace_du_dn(offset+k)(0)*poids;
              fevalY(k) = trace_du_dn(offset+k)(1)*poids;
            }
          
          int num_elem = interp_gamma.GetElementNumberOfSurface(n);
          int num_loc = interp_gamma.GetLocalPositionOfSurface(n);
          
          Fb.ComputeGaussIntegralSurfaceRef(fevalX, contribX, num_loc);
          Fb.ComputeGaussIntegralSurfaceRef(fevalY, contribY, num_loc);
          
          for (int k = 0; k < Fb.GetNbDofBoundary(num_loc); k++)
            {
              int num_dof_loc = Fb.GetLocalNumber(num_loc, k);
              int num_dof = var.mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
              int num_row = IndexDofNeumann(num_dof);
              Kh.AddInteraction(offset_gamma+num_row, offset_gamma+i, contribX(num_dof_loc));
              Kh.AddInteraction(offset_gamma+num_row+Ngamma, offset_gamma+i, contribY(num_dof_loc));
            }
          
          offset += interp_gamma.GetNbPointsQuadrature(n);
        }
      
      // on calcule la partie (grad phi, Pv) et (rot psi, Pv)
      Matrix2_2 dfjm1; R2 vec_u, vec_v;
      for (int n = 0; n < var.mesh.GetNbElt(); n++)
        {
          bool affine = var.mesh.IsElementAffine(n);
          for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
            {
              // on recupere J(xi_k) DFj^-1(xi_k)
              poids = Fb.WeightsND(k);
              if (affine)
                dfjm1 = var.Glob_DFjm1(n)(0);
              else
                dfjm1 = var.Glob_DFjm1(n)(k);
              
              vec_u(0) = Uquad(n)(k)(0)*poids;
              vec_u(1) = Uquad(n)(k)(1)*poids;
              
              Mlt(dfjm1, vec_u, vec_v);
              
              feval_grad(2*k) = vec_v(0);
              feval_grad(2*k+1) = vec_v(1);
              
              vec_u(0) = Uquad(n)(k)(1)*poids;
              vec_u(1) = -Uquad(n)(k)(0)*poids;

              Mlt(dfjm1, vec_u, vec_v);
              
              feval_rot(2*k) = vec_v(0);
              feval_rot(2*k+1) = vec_v(1);
            }
          
          // on integre contre grad phi
          Fb.ApplyRh(feval_grad, contribX);
          Fb.ApplyRh(feval_rot, contribY);
	  //DISP(feval_grad); DISP(feval_rot);
	  //DISP(contribX); DISP(contribY);
          
          // on rajoute la contribution sur la matrice Sh
          for (int k = 0; k < Fb.GetNbDof(); k++)
            {
              int num_dof = var.mesh_num.Element(n).GetNumberDof(k);
              Sh.AddInteraction(num_dof, offset_gamma+i, -contribX(k));
              Sh.AddInteraction(offset_gamma+i, num_dof, contribX(k));
              
              num_dof += var.mesh_num.GetNbDof();
              Sh.AddInteraction(num_dof, offset_gamma+i, -contribY(k));
              Sh.AddInteraction(offset_gamma+i, num_dof, contribY(k));
            }
        }
      
      // on ecrit la solution sur le disque
      //var_e.WriteOutputFile(sol_elas, string("UnElas")+to_str(i));

      //cout << "we wait " << endl;
      //int test_input; cin >> test_input;
    }

  // on revient sur le Dirichlet initial
  var_e.SetDirichletDofs(DirichletElas.GetM(), DirichletElas);

  // expression equivalente du dtn
  DistributedMatrix<Real_wp, General, ArrayRowSparse> Ah_elas;
  var_e.AddMatrixWithBC(Ah_elas, nat_mat);
  
  Matrix<Real_wp> Kp(2*Ngamma, 2*Ngamma);
  for (int i = 0; i < 2*Ngamma; i++)
    {
      Mlt(Ah_elas, mode_dtn(i), sol_elas);
      for (int j = 0; j < 2*Ngamma; j++)
	{
	  Kp(i, j) = DotProd(sol_elas, mode_dtn(j));
	  Kh.Get(offset_gamma+i, offset_gamma+j) = Kp(i, j);
	}
    }
  
  //Kp.WriteText("Kp.dat");
  
  //exit(0);
  
}

template<class TypeEquation>
void ComputeNormale(const EllipticProblem<TypeEquation>& var,
		    const Mesh<Dimension2>& mesh, Vector<Vector<Vector<R2> > >& NormaleMesh)
{
  NormaleMesh.Reallocate(mesh.GetNbElt());
  SetPoints<Dimension2> PointsElem;
  SetMatrices<Dimension2> MatricesElem;
  VectR2 s;
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int num_edge = i;
      int num_elem = mesh.BoundaryRef(i).numElement(0);
      int num_loc = mesh.Element(num_elem).GetPositionBoundary(num_edge);
      const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElement(num_elem);	  
      int nb_quad = Fb.GetNbQuadBoundary(num_loc);
      mesh.GetVerticesElement(num_elem, s);
      Fb.FjElem(s, PointsElem, mesh, num_elem);
      Fb.DFjElem(s, PointsElem, MatricesElem, mesh, num_elem);
      
      // restriction of DF_i on surface
      Fb.FjSurfaceElem(s, PointsElem, mesh, num_elem, num_loc);
      
      Fb.DFjSurfaceElem(s, PointsElem, MatricesElem,
			mesh, num_elem, num_loc);
      
      NormaleMesh(num_elem).Resize(4);
      NormaleMesh(num_elem)(num_loc).Reallocate(nb_quad);
      for (int k = 0; k < nb_quad; k++)
	NormaleMesh(num_elem)(num_loc)(k) = MatricesElem.GetNormaleQuadratureBoundary(k);
    }
}

// solver using fft in y-direction
class FourierElastoSolver
{
public :
  FftInterface<Complex_wp> fft;
  //! the solution is equal to \sum u_m exp(-i nu(m) y)
  //! nu_fft contains the values nu(m)
  VectReal_wp nu_fft; 
  //! number of points used in y-direction
  int nb_points_y;
  //! object used to compute 1-D finite element matrix
  EllipticProblem<LaplaceEquation1D> var_laplace;
  //! stiffness matrix \int dphi_i/dx dphi_j/dx 
  Matrix<Real_wp, General, ArrayRowSparse> mat_stiff;
  //! diagonal mass matrix \int phi_i phi_j 
  VectReal_wp mat_mass;
  //! physical coefficients
  Real_wp rho, lambda, mu;
  //! positions y
  VectReal_wp interval_y;
  //! computational domain
  Real_wp xmin, xmax, ymin, ymax;

  FourierElastoSolver()
  {
    nb_points_y = 2;
    rho = 1.0; lambda = 1.0; mu = 1.0;
    xmin = 0.0; xmax = 1.0; ymin = 0.0; ymax = 1.0;
  }
  
  template<class TypeEquation>
  void ConstructAll(EllipticProblem<TypeEquation>& var)
  {
    xmin = var.GetXmin(); xmax = var.GetXmax();
    var_laplace.mesh.SetOrder(var.mesh_num.GetOrder());
    
    var_laplace.SetBoundaryCondition(BoundaryConditionEnum::LINE_NEUMANN,
				     BoundaryConditionEnum::LINE_NEUMANN);
    
    int nbx = var.var_grid(0).GetNbPointsX();
    var_laplace.GetSectionGrid().Init(xmin, xmax, nbx);
    
    // 1-D regular mesh
    nbx = to_num<int>(var.mesh_data(0)(1));
    Vector<VectString> mesh_data(1);
    mesh_data.Reallocate(1);
    mesh_data(0).PushBack("REGULAR");
    mesh_data(0).PushBack(to_str(nbx));
    mesh_data(0).PushBack(to_str(xmin));
    mesh_data(0).PushBack(to_str(xmax));
    var_laplace.SetMeshData(mesh_data);
    
    var_laplace.SetOmega(var.GetOmega());
    
    var_laplace.ComputeMeshAndFiniteElement(string("EDGE_LOBATTO"));
    var_laplace.InitIndices(10);
    
    // mass matrix
    int Nvol = var_laplace.GetNbDof();
    var_laplace.ComputeMassMatrix();    

    var_laplace.GetMassMatrix(mat_mass);    
    
    // stiffness matrix
    mat_stiff.Reallocate(Nvol, Nvol);
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(0.0);
    var_laplace.AddMatrixFEM(mat_stiff, nat_mat);
    
    cout << "Combien de points en y (multiple de 2, 3, 5 de preference) ?" << endl;
    cin >> nb_points_y;

    // FFT parameters
    fft.Init(nb_points_y);
    nu_fft.Reallocate(nb_points_y);
    ymin = var.GetYmin(); ymax = var.GetYmax();
    for (int i = 0; i < nb_points_y/2; i++)
      {
        nu_fft(i) = 2.0*pi_wp*Real_wp(i) / (ymax - ymin);
        nu_fft(nb_points_y-1-i) = -2.0*pi_wp*Real_wp(i+1) / (ymax - ymin);
      }
    
    Linspace(ymin, ymax, nb_points_y+1, interval_y);
    interval_y.Resize(nb_points_y);
    
    // DISP(nu_fft); DISP(interval_y); DISP(var_laplace.CoorDofs);
  }
  
  // factorisation of finite element matrices and computation of the solution
  // rhs is the right hand side on input, the solution on output
  template<class TypeEquation>  
  void ComputeSolution(const EllipticProblem<TypeEquation>& var,
                       Vector<VectComplex_wp>& rhs)
  {
    cout << "Computing solution in Fourier space..." << endl;
    int Nvol = var_laplace.GetNbDof();
    Complex_wp val;
    All_MatrixLU<Complex_wp> mat_lu;
    for (int m = 0; m < nb_points_y; m++)
      {
        Matrix<Complex_wp, General, ArrayRowSparse> mat_glob;
        
        // stiffness terms
        mat_glob.Reallocate(2*Nvol, 2*Nvol);
        for (int i = 0; i < Nvol; i++)
          {
            for (int j = 0; j < mat_stiff.GetRowSize(i); j++)
              {
                int jcol = mat_stiff.Index(i, j);
                val = mat_stiff.Value(i, j);
                mat_glob.AddInteraction(i, jcol, val);
                mat_glob.AddInteraction(i+Nvol, jcol+Nvol, val);
              }
          }
        
        // mass term
        Real_wp m_omega2 = -var.GetSquareOmega();
        for (int i = 0; i < Nvol; i++)
          {
            val = rho/(lambda+2.0*mu)*m_omega2 + square(nu_fft(m));
            val *= mat_mass(i);
            mat_glob.AddInteraction(i, i, val);
            
            val = rho/mu*m_omega2 + square(nu_fft(m));
            val *= mat_mass(i);
            mat_glob.AddInteraction(Nvol+i, Nvol+i, val);                        
          }
        
        // Dirichlet terms for both sides
        mat_glob.AddInteraction(0, Nvol, Iwp*nu_fft(m));
        mat_glob.AddInteraction(Nvol-1, 2*Nvol-1, -Iwp*nu_fft(m));

        mat_glob.AddInteraction(Nvol, 0, -Iwp*nu_fft(m));
        mat_glob.AddInteraction(2*Nvol-1, Nvol-1, Iwp*nu_fft(m));
        
        // Neumann terms on the right side only
        // mode 0 is not treated to avoid division by 0
        if (m != 0)
          {
            val = -Iwp*rho/(2.0*nu_fft(m)*mu)*m_omega2;
            mat_glob.AddInteraction(Nvol-1, 2*Nvol-1, val);
            mat_glob.AddInteraction(2*Nvol-1, Nvol-1, -val);
          }
        else
          rhs(m).Fill(0);
        
        // resolution of linear system
        mat_lu.Factorize(mat_glob);
        mat_lu.Solve(rhs(m));
      }
  }

  
  // computation of eigenvalues
  template<class TypeEquation>  
  void ComputeEigenvalues(const EllipticProblem<TypeEquation>& var)
  {
    cout << "Computing eigenvalues in Fourier space..." << endl;
    int Nvol = var_laplace.GetNbDof();
    Complex_wp val;
    VectComplex_wp eigenval(2*Nvol*nb_points_y);
    Matrix<Complex_wp> eigen_vec;
    int num = 0;
    for (int m = 0; m < nb_points_y; m++)
      {
        // on utilise des matrices pleines hermitiennes
        Matrix<Complex_wp, General, RowMajor> Mh, Kh, Ah;
        
        // stiffness terms
        Mh.Reallocate(2*Nvol, 2*Nvol);
        Kh.Reallocate(2*Nvol, 2*Nvol);
        Mh.Fill(0); Kh.Fill(0);
        for (int i = 0; i < Nvol; i++)
          {
            for (int j = 0; j < mat_stiff.GetRowSize(i); j++)
              {
                int jcol = mat_stiff.Index(i, j);
                if (jcol >= i)
                  {
                    val = mat_stiff.Value(i, j);
                    Kh.Get(i, jcol) += val;
                    Kh.Get(i+Nvol, jcol+Nvol) += val;
                  }
              }
          }
        
        // mass term
        for (int i = 0; i < Nvol; i++)
          {
            val = rho/(lambda+2.0*mu)*mat_mass(i);
            Mh.Get(i, i) += val;
            
            val = square(nu_fft(m))*mat_mass(i);
            Kh.Get(i, i) += val;
            
            val = rho/mu*mat_mass(i);
            Mh.Get(Nvol+i, Nvol+i) += val;  
            
            val = square(nu_fft(m))*mat_mass(i);
            Kh.Get(Nvol+i, Nvol+i) += val;  
          }
        
        // Dirichlet terms for both sides
        Kh.Get(0, Nvol) += Iwp*nu_fft(m);
        Kh.Get(Nvol-1, 2*Nvol-1) += -Iwp*nu_fft(m);
        
        // Neumann terms on the right side only
        // mode 0 is not treated to avoid division by 0
        if (m != 0)
          {
            val = -Iwp*rho/(2.0*nu_fft(m)*mu);
            Mh.Get(Nvol-1, 2*Nvol-1) += val;
          }
        
	// for general matrices, lower part is enforced due to hermitian matrix
	for (int i = 0; i < 2*Nvol; i++)
	  for (int j = 0; j < i; j++)
	    {
	      Mh.Set(i, j, conj(Mh(j, i)));
	      Kh.Set(i, j, conj(Kh(j, i)));
	    }
	
	// Mh.Write("M.dat"); Kh.Write("K.dat");
	
        VectReal_wp w(2*Nvol);
        //GetEigenvalues(Kh, Mh, w); DISP(w);
	
	Vector<Complex_wp> wc;
	Ah.Reallocate(2*Nvol, 2*Nvol);
	GetInverse(Mh); Mlt(Mh, Kh, Ah);
	//GetEigenvalues(Ah, wc);
	GetEigenvaluesEigenvectors(Ah, wc, eigen_vec);
	
	for (int j = 0; j < 2*Nvol; j++)
	  w(j) = real(wc(j));
	
        for (int j = 0; j < 2*Nvol; j++)
          {
            if (w(j) < -1e-12)
              {
                cout << "Negative eigenvalue for m = " << m << endl;
                cout << "Value of this eigenvalue = " << w(j) << endl;
		cout << "nu_m = " << nu_fft(m) << endl;
		
		int Nx = var_laplace.GetSectionGrid().GetNbPointsGrid();
		VectComplex_wp trace_u(Nx);
		VectComplex_wp mode(Nvol);
		for (int i = 0; i < Nvol; i++)
		  mode(i) = eigen_vec(i, j);
		
		var_laplace.ComputeInterpolationU(mode, var_laplace.GetSectionGrid(), trace_u);
		
		int double_prec = OutputTypeEnum::DOUBLE_PRECISION;;
		bool ascii = false;
		
		GridInterpolationFull<Dimension1> grid;
		int type = GridInterpolationFull<Dimension1>::LINE;
		grid.Init(type, xmin, xmax, Nx);
		
		string name_file = "ModeP_M" + to_str(m) + ".dat";
		WriteMatlab(trace_u, grid, name_file, double_prec, ascii);

		for (int i = 0; i < Nvol; i++)
		  mode(i) = eigen_vec(Nvol+i, j);
		
		var_laplace.ComputeInterpolationU(mode, var_laplace.GetSectionGrid(), trace_u);
		
		name_file = "ModeS_M" + to_str(m) + ".dat";
		WriteMatlab(trace_u, grid, name_file, double_prec, ascii);
              }
            
            eigenval(num) = wc(j);
            num++;
          }
      }
    
    eigenval.Write("eigenval_mode.dat");
  }
  
  // computation of source
  void ComputeSource(GaussianSourceP<Real_wp>& fsrc, Vector<VectComplex_wp>& rhs)
  {
    int Nvol = var_laplace.GetNbDof();
    // allocation of arrays
    rhs.Reallocate(nb_points_y);
    for (int i = 0; i < nb_points_y; i++)
      {
        rhs(i).Reallocate(2*Nvol);
        rhs(i).Fill(0);
      }
    
    Vector<Complex_wp> x(nb_points_y);
    // loop over dofs in x-direction
    for (int i = 0; i < Nvol; i++)
      {
        Vector<Real_wp> f(1);
        // computation of f(x_i, y) for all points y
        for (int m = 0; m < nb_points_y; m++)
          {
            R2 pt(var_laplace.GetCoordinateDof(i), interval_y(m));
            fsrc.EvaluateFunction(pt, f);
            x(m) = f(0);
          }
        
        // then Fourier transform of the source is computed
        fft.ApplyInverse(x);
        
        // source is placed in phi_s equation
        for (int m = 0; m < nb_points_y; m++)
          rhs(m)(Nvol+i) = mat_mass(i)*x(m);
      }
  }
 
  // returns solution in real space
  void TransformSpace(VectComplex_wp& u, int Nx, VectComplex_wp& phi)
  {
    VectComplex_wp utmp;
    for (int i = 0; i < Nx; i++)
      {
        utmp.SetData(nb_points_y, &u(i*nb_points_y));
        
        // Fourier transform to retrieve solution on points y
        fft.ApplyForward(utmp);
        
        // the last point is duplicated and the array transposed
        for (int m = 0; m < nb_points_y; m++)
          phi(m*Nx + i) = utmp(m);
        
        phi(nb_points_y*Nx + i) = utmp(0);
        
        utmp.Nullify();
      }
  }
 
  // writing solution on disk
  void WriteSolution(Vector<VectComplex_wp>& sol)
  {
    int Nx = var_laplace.GetSectionGrid().GetNbPointsGrid();
    Vector<Complex_wp> trace_sol(Nx*nb_points_y);
    Vector<Complex_wp> sol_phi(Nx*(nb_points_y+1));
    int Nvol = var_laplace.GetNbDof();
    Vector<Complex_wp> x(Nvol);
    Vector<Complex_wp> trace_u(Nx);
    // 1-D interpolation for phi_p
    for (int m = 0; m < nb_points_y; m++)
      {
        for (int i = 0; i < Nvol; i++)
          x(i) = sol(m)(i);
        
        var_laplace.ComputeInterpolationU(x, var_laplace.GetSectionGrid(), trace_u);
        
        for(int i = 0; i < Nx; i++)
          trace_sol(m + i*nb_points_y) = trace_u(i);
      }
    
    // writing phi_p
    TransformSpace(trace_sol, Nx, sol_phi);

    int double_prec = OutputTypeEnum::DOUBLE_PRECISION;;
    bool ascii = false;
    
    GridInterpolationFull<Dimension2> grid;
    grid.SetPlaneOutput(xmin, xmax, ymin, ymax, Nx, nb_points_y+1);
    
    string name_file("totalCarrePhiP.dat");
    WriteMatlab(sol_phi, grid, name_file, double_prec, ascii);    
    
    // 1-D interpolation for phi_s
    for (int m = 0; m < nb_points_y; m++)
      {
        for (int i = 0; i < Nvol; i++)
          x(i) = sol(m)(Nvol + i);
        
        var_laplace.ComputeInterpolationU(x, var_laplace.GetSectionGrid(), trace_u);
        
        for(int i = 0; i < Nx; i++)
          trace_sol(m + i*nb_points_y) = trace_u(i);
      }
    
    // writing phi_s
    TransformSpace(trace_sol, Nx, sol_phi);
 
    name_file = string("totalCarrePhiS.dat");
    WriteMatlab(sol_phi, grid, name_file, double_prec, ascii);

    // gradient of PhiP
    for (int m = 0; m < nb_points_y; m++)
      {
        for (int i = 0; i < Nvol; i++)
          x(i) = sol(m)(i);
        
        var_laplace.ComputeInterpolationU(x, var_laplace.GetSectionGrid(), trace_u, false, true);
        
        for(int i = 0; i < Nx; i++)
          trace_sol(m + i*nb_points_y) = trace_u(i);
      }
    
    TransformSpace(trace_sol, Nx, sol_phi);
    
    name_file = string("totalCarrePhiP_dU0.dat");
    WriteMatlab(sol_phi, grid, name_file, double_prec, ascii);    

    for (int m = 0; m < nb_points_y; m++)
      {
        for (int i = 0; i < Nvol; i++)
          x(i) = sol(m)(i);
        
        var_laplace.ComputeInterpolationU(x, var_laplace.GetSectionGrid(), trace_u);
        Mlt(-Iwp*nu_fft(m), trace_u);
	
        for(int i = 0; i < Nx; i++)
          trace_sol(m + i*nb_points_y) = trace_u(i);
      }
    
    TransformSpace(trace_sol, Nx, sol_phi);
    
    name_file = string("totalCarrePhiP_dU1.dat");
    WriteMatlab(sol_phi, grid, name_file, double_prec, ascii);    

    // gradient of PhiS
    for (int m = 0; m < nb_points_y; m++)
      {
        for (int i = 0; i < Nvol; i++)
          x(i) = sol(m)(Nvol + i);
        
        var_laplace.ComputeInterpolationU(x, var_laplace.GetSectionGrid(), trace_u, false, true);
        
        for(int i = 0; i < Nx; i++)
          trace_sol(m + i*nb_points_y) = trace_u(i);
      }
    
    TransformSpace(trace_sol, Nx, sol_phi);
    
    name_file = string("totalCarrePhiS_dU0.dat");
    WriteMatlab(sol_phi, grid, name_file, double_prec, ascii);    

    for (int m = 0; m < nb_points_y; m++)
      {
        for (int i = 0; i < Nvol; i++)
          x(i) = sol(m)(Nvol + i);
        
        var_laplace.ComputeInterpolationU(x, var_laplace.GetSectionGrid(), trace_u);
        Mlt(-Iwp*nu_fft(m), trace_u);
	
        for(int i = 0; i < Nx; i++)
          trace_sol(m + i*nb_points_y) = trace_u(i);
      }
    
    TransformSpace(trace_sol, Nx, sol_phi);
    
    name_file = string("totalCarrePhiS_dU1.dat");
    WriteMatlab(sol_phi, grid, name_file, double_prec, ascii);    
  }
  
};

class ElastoParam
{
public :
  
  // HARMONIC_RESOL : resolution en harmonique
  // TIME_RESOL : resolution en temporel
  // MODE_RESOL : calcul des modes
  enum{HARMONIC_RESOL = 0, TIME_RESOL = 1, MODE_RESOL = 2};
  
  // UU_FORMUL : formulation en introduisant u
  // VV_FORMUL : formulation en introduisant v
  // PSI_FORMUL : formulation en introduisant u et psi = dphi/dn
  // TAU_FORMUL : formulation en introduisant u, et en derivant par tau l'equation de u
  // TAUPHI_FORMUL : formulation en phi uniquement en remplacant u matriciellement (en inversant le laplacien)
  // PHI_FORMUL : formulation en phi uniquement en remplacant u par la primitive de phi en espace
  // DTN_FORMUL : formulation en phi et u avec un DtN sur u
  // FOURIER_FORMUL : Fourier en y pour avoir des pbs 1-D ne dependant que de x
  enum{UU_FORMUL = 0, VV_FORMUL = 1, PSI_FORMUL = 2, TAU_FORMUL = 3, TAUPHI_FORMUL = 4, PHI_FORMUL = 5,
       DTN_FORMUL = 6, FOURIER_FORMUL = 7};

};


void ModifyMatrixMean(Matrix<Real_wp, General, ArrayRowSparse>& A,
		      Vector<TinyVector<int, 2> >& couple, VectReal_wp& coef_couple,
		      IVect& removed_row, const Real_wp& coef_diag)
{
  // replacing rows
  for (int i = 0; i < couple.GetM(); i++)
    {
      int irow = couple(i)(0);
      int irow2 = couple(i)(1);
      int nb_val = A.GetRowSize(irow2);
      IVect col_num(nb_val); VectReal_wp val(nb_val);
      for (int j = 0; j < nb_val; j++)
	{
	  col_num(j) = A.Index(irow2, j);
	  val(j) = -coef_couple(i)*A.Value(irow2, j);
	}
      
      A.AddInteractionRow(irow, nb_val, col_num, val);
    }
  
  // replacing columns
  Vector<Vector<Real_wp, VectSparse>, VectSparse> cols;
  GetCol(A, removed_row, cols);
  
  for (int i = 0; i < couple.GetM(); i++)
    {
      int icol = couple(i)(0);
      int icol2 = couple(i)(1);
      for (int j = 0; j < cols(icol2).GetM(); j++)
	A.AddInteraction(cols(icol2).Index(j), icol, -coef_couple(i)*cols(icol2).Value(j));
    }
  
  for (int i = 0; i < removed_row.GetM(); i++)
    A.ClearRow(removed_row(i));
  
  EraseCol(removed_row, A);
  
  // diagonal coefficient ?
  if (coef_diag != 0.0)
    {
      for (int i = 0; i < removed_row.GetM(); i++)
	{
	  int irow = removed_row(i);
	  A.ReallocateRow(irow, 1);
	  A.Index(irow, 0) = irow;
	  A.Value(irow, 0) = coef_diag;
	}
    }
}

void ModifyVectorMean(Vector<Real_wp>& x, Vector<TinyVector<int, 2> >& couple,
		      VectReal_wp& coef_couple, IVect& removed_row)
{
  for (int i = 0; i < couple.GetM(); i++)
    {
      int irow = couple(i)(0);
      int irow2 = couple(i)(1);
      if (irow < x.GetM())
        x(irow) -= coef_couple(i)*x(irow2);
    }
  
  for (int i = 0; i < removed_row.GetM(); i++)
    if (removed_row(i) < x.GetM())
      x(removed_row(i)) = 0.0;

}

void RecomposeVectorMean(Vector<Real_wp>& x, Vector<TinyVector<int, 2> >& couple,
			 VectReal_wp& coef_couple, IVect& removed_row)
{
  VectReal_wp x0(couple.GetM());
  for (int i = 0; i < couple.GetM(); i++)
    if (couple(i)(0) < x.GetM())
      x0(i) = x(couple(i)(0));

  for (int i = 0; i < removed_row.GetM(); i++)
    if (removed_row(i) < x.GetM())
      x(removed_row(i)) = 0.0;
  
  for (int i = 0; i < couple.GetM(); i++)
    if (couple(i)(1) < x.GetM())
      x(couple(i)(1)) -= coef_couple(i)*x0(i);
}

template<class TypeEqAcous, class TypeEqElas>
void RunAll(HyperbolicProblem<TypeEqAcous>& var_s,
	    HyperbolicProblem<TypeEqElas>& var_elas, string input_file, string name_element,
	    int type_resolution, int type_formulation)
{
  //typedef typename TypeEqElas::TypeEquationStationary TypeEqElasStatic;
  HyperbolicProblem<TypeEqAcous> var_p;
  
  // si true, les ddls de phi_p et phi_s sur le bord sont combines de telle
  // sorte a avoir la moyenne de phi_p et phi_s egale a zero
  bool mean_value_phi_zero = false;
  cout << "Tapez 1 si vous voulez des fonctions de base a moyenne nulle" << endl;
  cin >> mean_value_phi_zero;
  DISP(mean_value_phi_zero);
  
  // on construit les matrices du laplacien pour les les ondes P et les ondes S
  var_p.var_harmonic.InitIndices(20);
  var_s.var_harmonic.InitIndices(20);
  var_elas.var_harmonic.InitIndices(20);
  
  ReadInputFile(input_file, var_s);
  ReadInputFile(input_file, var_p);
  ReadInputFile(input_file, var_elas);
  
  // on force Neumann pour avoir la matrice du laplacien sans termes de bord
  bool presence_neumann = false;
  for (int ref = 1; ref < var_s.var_harmonic.mesh.GetNbReferences(); ref++)
    {
      int cond = var_s.var_harmonic.mesh.GetBoundaryCondition(ref);
      if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
        {
          var_s.var_harmonic.mesh.SetBoundaryCondition(ref, BoundaryConditionEnum::LINE_NEUMANN);
          var_p.var_harmonic.mesh.SetBoundaryCondition(ref, BoundaryConditionEnum::LINE_NEUMANN);
        }
      else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
	presence_neumann = true;
    }
  
  // on met les coefs rho/(lambda + 2 mu) et rho/mu sur les 
  // matrices de masse pour phi_p et phi_s
  Real_wp lambda(1), mu(1), rho(1);
  for (int ref = 0; ref < var_elas.var_harmonic.ref_rho.GetM(); ref++)
    {
      lambda = var_elas.var_harmonic.ref_tensorC(ref)(0, 0, 1, 1);
      mu = var_elas.var_harmonic.ref_tensorC(ref)(0, 1, 0, 1);
      rho = var_elas.var_harmonic.ref_rho(ref).GetConstant();
      var_s.var_harmonic.ref_mu(ref).SetIdentity();
      var_p.var_harmonic.ref_mu(ref).SetIdentity();
      var_p.var_harmonic.ref_rho(ref).SetConstant(rho/(lambda+2.0*mu));
      var_s.var_harmonic.ref_rho(ref).SetConstant(rho/mu);
      var_s.var_harmonic.ref_sigma(ref).SetConstant(0.0);
      var_p.var_harmonic.ref_sigma(ref).SetConstant(0.0);
    }
  
  // on peut utiliser var_p et var_s pour sortir U1 et U2
  for (int i = 0; i < var_p.var_harmonic.output_grid_param.GetM(); i++)
    {
      string name = var_p.var_harmonic.output_grid_param(i).GetTotalFieldFile();
      name = GetBaseString(name) + "P";
      var_p.var_harmonic.output_grid_param(i).SetTotalFieldFile(name);
      var_p.var_harmonic.output_grid_param(i).SetDiffractedFieldFile(name+"U1");
    }

  for (int i = 0; i < var_s.var_harmonic.output_grid_param.GetM(); i++)
    {
      string name = var_s.var_harmonic.output_grid_param(i).GetTotalFieldFile();
      name = GetBaseString(name) + "S";
      var_s.var_harmonic.output_grid_param(i).SetTotalFieldFile(name);
      var_s.var_harmonic.output_grid_param(i).SetDiffractedFieldFile(name+"U2");
    }

  // dans le cas Dirichlet pur, on peut utiliser un ordre different
  if (!presence_neumann)
    {
      int rp = var_p.var_harmonic.mesh_num.GetOrder(); DISP(rp);
      int rs = 0;
      cout << "Tapez 0 pour utiliser le meme ordre pour phi_p et phi_s " << endl;
      cout << "Sinon tapez l'ordre a utiliser pour phi_s" << endl;
      cout << "L'ordre pour phi_p est celui du fichier de donnees" << endl;
      cin >> rs;
      
      if (rs > 0)
	{
	  var_s.var_harmonic.mesh_num.SetOrder(rs);	  
	  if (type_formulation > 1)
	    {
	      cout << "Different orders are not allowed for this formulation" << endl;
	      abort();
	    }
	}
    }
  
  // condition periodique imposee de maniere forte
  var_p.var_harmonic.mesh_num.SetSameNumberPeriodicDofs();
  var_s.var_harmonic.mesh_num.SetSameNumberPeriodicDofs();
  var_elas.var_harmonic.mesh_num.SetSameNumberPeriodicDofs();
  
  // on construit les elements finis et le maillage
  var_p.var_harmonic.ComputeMeshAndFiniteElement(name_element);
  var_s.var_harmonic.ComputeMeshAndFiniteElement(name_element);
  var_elas.var_harmonic.ComputeMeshAndFiniteElement(name_element);

  // dans le cas d'ordres differents, il faut calculer les operateurs d'interpolation
  if (var_s.var_harmonic.mesh_num.GetOrder() != var_p.var_harmonic.mesh_num.GetOrder())
    {
      int rp = var_p.var_harmonic.mesh_num.GetOrder();
      int rs = var_s.var_harmonic.mesh_num.GetOrder();
      Vector<bool> order_present(rs+1);
      order_present.Fill(false); order_present(rs) = true;
      if (var_p.var_harmonic.GetReferenceElement(0).GetOrder() == rp)
	var_p.var_harmonic.GetReferenceElement(0)
	  .ComputeVariableInterpolation(rs, order_present);
    }
  
  // autres initialisations
  var_p.var_harmonic.PerformOtherInitializations();
  var_s.var_harmonic.PerformOtherInitializations();
  var_elas.var_harmonic.PerformOtherInitializations();

  // initialisation resolution de Fourier
  FourierElastoSolver var_fft;
  if (type_formulation == ElastoParam::FOURIER_FORMUL)
    {
      var_fft.rho = rho; var_fft.lambda = lambda; var_fft.mu = mu;
      var_fft.ConstructAll(var_p.var_harmonic);
    }
  
  // calcul des transfos geometriques Fi, DFi, etc
  var_p.var_harmonic.ComputeMassMatrix();
  var_s.var_harmonic.ComputeMassMatrix();
  var_elas.var_harmonic.ComputeMassMatrix();
  Vector<Vector<Vector<R2> > > NormaleMesh;
  ComputeNormale(var_elas.var_harmonic, var_elas.var_harmonic.mesh, NormaleMesh);
  
  // pour la periodicite
  var_p.var_harmonic.ComputeQuasiPeriodicPhase();
  var_s.var_harmonic.ComputeQuasiPeriodicPhase();
  var_elas.var_harmonic.ComputeQuasiPeriodicPhase();

  // on construit ce qui est necessaire pour faire avancer l'elasto classique
  var_elas.ComputeStiffnessMatrix();
  var_elas.ComputeMassMatrix();
  
  // initialisation des schemas en temps
  var_p.InitTimeIterations();
  var_s.InitTimeIterations();
  var_elas.InitTimeIterations();

  // Ngamma : nombre de ddls sur la frontiere Neumann (Gamma_N
  // IndexDofNeumann(i) : numero local du ddl sur Gamma_N
  IVect IndexDofNeumann, ListeDofNeu, IndexDofDirichlet; int Ngamma;
  GetDofOnGammaNeumann(var_elas.var_harmonic, ListeDofNeu, IndexDofNeumann, Ngamma);
  
  // IndexDofDirichlet : identite
  IndexDofDirichlet.Reallocate(2*var_s.var_harmonic.mesh_num.GetNbDof());
  IndexDofDirichlet.Fill();

  var_p.var_harmonic.mesh.Write("test.mesh");
  
  // nombre de ddls
  int Np = var_p.var_harmonic.GetNbDof();
  int Ns = var_s.var_harmonic.GetNbDof();
  cout << endl;
  cout << "Number of dofs for phiP = " << Np << endl;
  cout << "Number of dofs for phiS = " << Ns << endl;  
  cout << endl;
  int Ndof = Np + Ns + 2*Ngamma;
  switch (type_formulation)
    {
    case ElastoParam::PSI_FORMUL :
      Ndof += 2*Ngamma;
      break;
    case ElastoParam::TAUPHI_FORMUL :
    case ElastoParam::PHI_FORMUL :
      Ndof -= 2*Ngamma;
      break;
    }

  Vector<TinyVector<int, 2> > ListeCouple; VectReal_wp CoefCouple; IVect removed_dof;
  VectR2 CoorDofs;
  GetCoordinateDofs(var_s.var_harmonic, CoorDofs);
  
  if (mean_value_phi_zero)
    {
      GetIndependentDofs(var_elas.var_harmonic, ListeCouple, CoefCouple, removed_dof, BoundaryConditionEnum::LINE_NEUMANN);
      GetCoordinateDofs(var_elas.var_harmonic, CoorDofs);
      //DISP(ListeCouple); DISP(CoefCouple);
      
      //for (int i = 0; i < CoorDofs.GetM(); i++)
      //cout << "Ddl numero " << i << " : " << CoorDofs(i) << endl;
      
      ofstream file_out("CoorDofs.dat");
      for (int i = 0; i < CoorDofs.GetM(); i++)
	file_out << "Dof " << i << " x = " << CoorDofs(i)(0) << " y = " << CoorDofs(i)(1) << endl;
      
      file_out.close();
      
      int nb_remove_dof = removed_dof.GetM();
      int nb_couple = ListeCouple.GetM();
      removed_dof.Resize(2*nb_remove_dof);
      ListeCouple.Resize(2*nb_couple);
      CoefCouple.Resize(2*nb_couple);
      for (int i = 0; i < nb_couple; i++)
        {
          ListeCouple(i+nb_couple)(0) = ListeCouple(i)(0) + Np;
          ListeCouple(i+nb_couple)(1) = ListeCouple(i)(1) + Np;
	  CoefCouple(i+nb_couple) = CoefCouple(i);
        }
      
      for (int i = 0; i < nb_remove_dof; i++)
        removed_dof(i+nb_remove_dof) = removed_dof(i) + Np;
    }
  
  // initialisation des matrices
  // on met le systeme sous la forme Mh d^2 U/dt^2 + Sh dU/dt + Kh U = 0
  DistributedMatrix<Real_wp, General, ArrayRowSparse> Kh(Ndof, Ndof), Mh(Ndof, Ndof), Sh(Ndof, Ndof);
  
  int nb_ref = 20;
  Vector<bool> dirichlet_cond(nb_ref+1), neumann_cond(nb_ref+1);
  dirichlet_cond.Fill(false); neumann_cond.Fill(false);
  dirichlet_cond(BoundaryConditionEnum::LINE_DIRICHLET) = true;
  neumann_cond(BoundaryConditionEnum::LINE_NEUMANN) = true;
  
  // on assemble la partie commune a Dirichlet/Neumann
  GetMassMatrix(var_p.var_harmonic, var_s.var_harmonic, Mh);
  GetStiffnessMatrix(var_p.var_harmonic, var_s.var_harmonic, Kh);
    
  AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
		var_s.var_harmonic, var_elas.var_harmonic,
                Kh, IndexDofDirichlet, IndexDofDirichlet, dirichlet_cond,
		Np, 0, false, false);

  AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
		var_s.var_harmonic, var_elas.var_harmonic,
                Kh, IndexDofDirichlet, IndexDofDirichlet, dirichlet_cond,
		0, Np, true, false);
  
  // est-ce qu'on veut partir d'une solution "lisse" ?
  // ca sert si on veut partir d'une solution jolie afin de toucher le bord
  // et de voir l'instabilite generee par le bord
  bool start_smooth_solution = false;
  //cout << " Voulez vous projeter la solution de l'elasto sur phi_P, phi_S ? " << endl;
  //cout << " Tapez 0 pour non" << endl;
  //cin >> start_smooth_solution; DISP(start_smooth_solution);
  
  //Real_wp t_start = 3.0;
  All_MatrixLU<Real_wp> mat_lu_laplace;
  if (start_smooth_solution)
    {
      // factorisation de la partie "Dirichlet" de la matrice
      Matrix<Real_wp, General, ArrayRowSparse> Kh_extract(Np+Ns, Np+Ns);
      for (int i = 0; i < Np+Ns; i++)
	{
	  int size_row = Kh.GetRowSize(i);
	  Kh_extract.ReallocateRow(i, size_row);
	  for (int j = 0; j < size_row; j++)
	    {
	      int jcol = Kh.Index(i, j);
	      if (jcol < (Np+Ns))
		{
		  Kh_extract.Index(i, j) = jcol;
		  Kh_extract.Value(i, j) = Kh.Value(i, j);
		}
	      else
		{
		  abort();
		}
	    }
	}
      
      // on rajoute le Rpsi "Neumann"
      AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
		    var_s.var_harmonic, var_elas.var_harmonic,
		    Kh_extract, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
		    Np, 0, false, false);
      
      AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
		    var_s.var_harmonic, var_elas.var_harmonic,
		    Kh_extract, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
		    0, Np, true, false);
      
      /* Kh_extract.WriteText("KhDir.dat");
      
         if (mean_value_phi_zero)
         ModifyMatrixMean(Kh_extract, ListeCouple, CoefCouple, removed_dof, 1.0);
      
         Kh_extract.WriteText("KhModif.dat"); */
      
      mat_lu_laplace.Factorize(Kh_extract);
    }
  

  if (Ngamma > 0)
    {
      // partie due a Neumann
      switch (type_formulation)
	{
	case ElastoParam::UU_FORMUL :
	  {
	    // Dirichlet part
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  Np, 0, false, false);
	    
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  0, Np, true, false);

	    // part u.n phi_p and n x u phi_s
	    AddMatrixCh(0.0, -1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Sh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, 0, Np+Ns);

	    AddMatrixCh(0.0, 0.0, -1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Sh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, 0, Np+Ns+Ngamma);

	    AddMatrixCh(0.0, 0.0, -1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Sh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, Np, Np+Ns);

	    AddMatrixCh(0.0, 1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Sh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, Np, Np+Ns+Ngamma);
	    
	    // part n1 dphi_p/dt + n2 dphi_s/dt
	    //      n2 dphi_p/dt - n1 dphi_s/dt
	    AddMatrixCh(0.0, 1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Sh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns, 0);

	    AddMatrixCh(0.0, 0.0, 1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Sh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns, Np);

	    AddMatrixCh(0.0, 0.0, 1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Sh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns+Ngamma, 0);

	    AddMatrixCh(0.0, -1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Sh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns+Ngamma, Np);
	    
	    // part 2 mu/rho( -n x grad u2  and n x grad u1 )
	    AddMatrixRpsi(-2.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofNeumann, IndexDofNeumann, neumann_cond,
			  Np+Ns, Np+Ns+Ngamma, true, true);

	    AddMatrixRpsi(2.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofNeumann, IndexDofNeumann, neumann_cond,
			  Np+Ns+Ngamma, Np+Ns, true, true);	    
	  }
	  break;
	case ElastoParam::VV_FORMUL :
	  {
	    // Dirichlet part
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  Np, 0, false, false);
	    
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  0, Np, true, false);

	    // part u.n phi_p and n x u phi_s
	    AddMatrixCh(0.0, -1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, 0, Np+Ns);

	    AddMatrixCh(0.0, 0.0, -1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, 0, Np+Ns+Ngamma);

	    AddMatrixCh(0.0, 0.0, -1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, Np, Np+Ns);

	    AddMatrixCh(0.0, 1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, Np, Np+Ns+Ngamma);
	    
	    // part n1 d^2 phi_p/dt + n2 d^2 phi_s/dt
	    //      n2 d^2 phi_p/dt - n1 d^2 phi_s/dt
	    AddMatrixCh(0.0, 1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Mh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns, 0);

	    AddMatrixCh(0.0, 0.0, 1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Mh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns, Np);

	    AddMatrixCh(0.0, 0.0, 1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Mh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns+Ngamma, 0);

	    AddMatrixCh(0.0, -1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Mh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns+Ngamma, Np);
	    
	    // part 2 mu/rho( -n x grad u2  and n x grad u1 )
	    AddMatrixRpsi(-2.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofNeumann, IndexDofNeumann, neumann_cond,
			  Np+Ns, Np+Ns+Ngamma, true, true);

	    AddMatrixRpsi(2.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofNeumann, IndexDofNeumann, neumann_cond,
			  Np+Ns+Ngamma, Np+Ns, true, true);	    
	  }
	  break;
	case ElastoParam::PSI_FORMUL :
	  {
	    // part - psi^S tilde{phi}_s  and -psi_p tilde{phi}_p
	    AddMatrixCh(-1.0, 0.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, 0, Np+Ns+2*Ngamma);

	    AddMatrixCh(-1.0, 0.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, Np, Np+Ns+3*Ngamma);
	    
	    // part n1 v1 + n2 v2 - psi_p - n1 dphi_S/dx_2 + n2 dphi_S/dx_1 = 0
	    AddMatrixCh(0.0, 1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofNeumann, IndexDofNeumann, neumann_cond, Np+Ns+2*Ngamma, Np+Ns);
	    
	    AddMatrixCh(0.0, 0.0, 1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofNeumann, IndexDofNeumann, neumann_cond, Np+Ns+2*Ngamma, Np+Ns+Ngamma);
	    
	    AddMatrixCh(-1.0, 0.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofNeumann, IndexDofNeumann, neumann_cond, Np+Ns+2*Ngamma, Np+Ns+2*Ngamma);
	    
	    AddMatrixRpsi(-1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofNeumann, IndexDofDirichlet, neumann_cond,
			  Np+Ns+2*Ngamma, Np, true, false);

	    // part v1 n2 - n1 v2 - psi_s - n2 dphi_p/dx_1 + n1 dphi_p/dx_2
	    AddMatrixCh(0.0, 0.0, 1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofNeumann, IndexDofNeumann, neumann_cond, Np+Ns+3*Ngamma, Np+Ns);
	    
	    AddMatrixCh(0.0, -1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofNeumann, IndexDofNeumann, neumann_cond, Np+Ns+3*Ngamma, Np+Ns+Ngamma);
	    
	    AddMatrixCh(-1.0, 0.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofNeumann, IndexDofNeumann, neumann_cond, Np+Ns+3*Ngamma, Np+Ns+3*Ngamma);
	    
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofNeumann, IndexDofDirichlet, neumann_cond,
			  Np+Ns+3*Ngamma, 0, true, false);	    
	    
	    // part n1 d^2 phi_p/dt + n2 d^2 phi_s/dt
	    //      n2 d^2 phi_p/dt - n1 d^2 phi_s/dt
	    AddMatrixCh(0.0, 1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Mh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns, 0);

	    AddMatrixCh(0.0, 0.0, 1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Mh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns, Np);

	    AddMatrixCh(0.0, 0.0, 1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Mh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns+Ngamma, 0);

	    AddMatrixCh(0.0, -1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Mh,
			IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns+Ngamma, Np);
	    
	    // part 2 mu/rho( -n x grad u2  and n x grad u1 )
	    AddMatrixRpsi(-2.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofNeumann, IndexDofNeumann, neumann_cond,
			  Np+Ns, Np+Ns+Ngamma, true, true);

	    AddMatrixRpsi(2.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofNeumann, IndexDofNeumann, neumann_cond,
			  Np+Ns+Ngamma, Np+Ns, true, true);	    	    
	  }
	  break;
	case ElastoParam::TAU_FORMUL :
	  {
	    // Dirichlet part
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  Np, 0, false, false);
	    
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  0, Np, true, false);

	    // part u.n phi_p and n x u phi_s
	    AddMatrixCh(0.0, -1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, 0, Np+Ns);

	    AddMatrixCh(0.0, 0.0, -1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, 0, Np+Ns+Ngamma);

	    AddMatrixCh(0.0, 0.0, -1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, Np, Np+Ns);

	    AddMatrixCh(0.0, 1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, Kh,
			IndexDofDirichlet, IndexDofNeumann, neumann_cond, Np, Np+Ns+Ngamma);
	    
	    // part n1 d^2 phi_p/dt + n2 d^2 phi_s/dt
	    //      n2 d^2 phi_p/dt - n1 d^2 phi_s/dt
	    AddMatrixRpsi(0.0, 1.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic, Mh,
			  IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns, 0, false, false);

	    AddMatrixRpsi(0.0, 0.0, 1.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic, Mh,
			  IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns, Np, false, false);

	    AddMatrixRpsi(0.0, 0.0, 1.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic, Mh,
			  IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns+Ngamma, 0, false, false);

	    AddMatrixRpsi(0.0, -1.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic, Mh,
			  IndexDofNeumann, IndexDofDirichlet, neumann_cond, Np+Ns+Ngamma, Np, false, false);
	    
	    // part 2 mu/rho( -n x grad u2  and n x grad u1 )
	    AddMatrixLaplacian(-2.0, 0.0, 0.0, var_s.var_harmonic, var_elas.var_harmonic,
			       Kh, IndexDofNeumann, IndexDofNeumann, neumann_cond,
			       Np+Ns, Np+Ns+Ngamma, true);
	    
	    AddMatrixLaplacian(2.0, 0.0, 0.0, var_s.var_harmonic, var_elas.var_harmonic,
			       Kh, IndexDofNeumann, IndexDofNeumann, neumann_cond,
			       Np+Ns+Ngamma, Np+Ns, true);

	    // mean value of phi_p and phi_s is equal to 0
	    //AddMatrixPhi(1.0, 0.0, 0.0,
	    //var_s.var_harmonic, var_elas.var_harmonic, Mh,
	    //IndexDofDirichlet, neumann_cond, Np+Ns+2*Ngamma, 0, false);

	    //AddMatrixPhi(1.0, 0.0, 0.0,
	    //var_s.var_harmonic, var_elas.var_harmonic, Mh,
	    //IndexDofNeumann, neumann_cond, Np+Ns, Np+Ns+2*Ngamma, true);
	    
	    //AddMatrixPhi(1.0, 0.0, 0.0,
	    //var_s.var_harmonic, var_elas.var_harmonic, Mh,
	    //IndexDofDirichlet, neumann_cond, Np+Ns+2*Ngamma+1, Np, false);
	    
	    //AddMatrixPhi(1.0, 0.0, 0.0,
	    //var_s.var_harmonic, var_elas.var_harmonic, Mh,
	    //IndexDofNeumann, neumann_cond, Np+Ns+Ngamma, Np+Ns+2*Ngamma+1, true);
	  }
	  break;
	case ElastoParam::TAUPHI_FORMUL :
	  {
	    // Dirichlet part
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  Np, 0, false, false);
	    
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  0, Np, true, false);
	    
	    DistributedMatrix<Real_wp, General, ArrayRowSparse> A12(2*Ngamma, 2*Ngamma), A21(2*Ngamma, 2*Ngamma), A22(2*Ngamma, 2*Ngamma);
	    // part u.n phi_p and n x u phi_s
	    AddMatrixCh(0.0, -1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, A12,
			IndexDofNeumann, IndexDofNeumann, neumann_cond, 0, 0);

	    AddMatrixCh(0.0, 0.0, -1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, A12,
			IndexDofNeumann, IndexDofNeumann, neumann_cond, 0, Ngamma);

	    AddMatrixCh(0.0, 0.0, -1.0, 
			var_s.var_harmonic, var_elas.var_harmonic, A12,
			IndexDofNeumann, IndexDofNeumann, neumann_cond, Ngamma, 0);

	    AddMatrixCh(0.0, 1.0, 0.0, 
			var_s.var_harmonic, var_elas.var_harmonic, A12,
			IndexDofNeumann, IndexDofNeumann, neumann_cond, Ngamma, Ngamma);
	    
	    // part n1 d^2 phi_p/dt + n2 d^2 phi_s/dt
	    //      n2 d^2 phi_p/dt - n1 d^2 phi_s/dt
	    AddMatrixRpsi(0.0, 1.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic, A21,
			  IndexDofNeumann, IndexDofNeumann, neumann_cond, 0, 0, false, false);

	    AddMatrixRpsi(0.0, 0.0, 1.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic, A21,
			  IndexDofNeumann, IndexDofNeumann, neumann_cond, 0, Ngamma, false, false);

	    AddMatrixRpsi(0.0, 0.0, 1.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic, A21,
			  IndexDofNeumann, IndexDofNeumann, neumann_cond, Ngamma, 0, false, false);

	    AddMatrixRpsi(0.0, -1.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic, A21,
			  IndexDofNeumann, IndexDofNeumann, neumann_cond, Ngamma, Ngamma, false, false);
	    
	    // part 2 mu/rho( -n x grad u2  and n x grad u1 )
	    AddMatrixLaplacian(-2.0, 0.0, 0.0, var_s.var_harmonic, var_elas.var_harmonic,
			       A22, IndexDofNeumann, IndexDofNeumann, neumann_cond,
			       0, Ngamma, true);
	    
	    AddMatrixLaplacian(2.0, 0.0, 0.0, var_s.var_harmonic, var_elas.var_harmonic,
			       A22, IndexDofNeumann, IndexDofNeumann, neumann_cond,
			       Ngamma, 0, true);	    
	    
	    // forming Schur complement
	    Matrix<Real_wp> a12, a21, a22, a_tmp(2*Ngamma, 2*Ngamma), schur(2*Ngamma, 2*Ngamma);
	    Copy(A12, a12); Copy(A21, a21); Copy(A22, a22);
	    	    
	    GetPseudoInverse(a22, 1e-10);
	    a_tmp.Fill(0); Mlt(a22, a21, a_tmp);
	    schur.Fill(0); Mlt(a12, a_tmp, schur);
	    
	    for (int i = 0; i < Ngamma; i++)
	      for (int j = 0; j < Ngamma; j++)
		{
		  Mh.AddInteraction(ListeDofNeu(i), ListeDofNeu(j), -schur(i, j));
		  Mh.AddInteraction(Np+ListeDofNeu(i), ListeDofNeu(j), -schur(Ngamma+i, j));
		  Mh.AddInteraction(ListeDofNeu(i), Np+ListeDofNeu(j), -schur(i, Ngamma+j));
		  Mh.AddInteraction(Np+ListeDofNeu(i), Np+ListeDofNeu(j), -schur(Ngamma+i, Ngamma+j));
		}
	  }
	  break;
	case ElastoParam::PHI_FORMUL :
	  {
	    // Dirichlet part
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  Np, 0, false, false);
	    
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  0, Np, true, false);
	    
	    // primitive part
	    int ref = var_s.var_harmonic.mesh.Element(0).GetReference();
	    Real_wp mu_rho = 1.0/var_s.var_harmonic.ref_rho(ref).GetConstant();
	    Real_wp rho_div_2mu = mu_rho/2.0;
	    AddMatrixIntPhi(-rho_div_2mu, var_s.var_harmonic, var_elas.var_harmonic,
			    Mh, neumann_cond, 0, Np);
	    
	    AddMatrixIntPhi(rho_div_2mu, var_s.var_harmonic, var_elas.var_harmonic,
			    Mh, neumann_cond, Np, 0);
	  }
	  break;
        case ElastoParam::DTN_FORMUL :
          {
            // Dirichlet part
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  Np, 0, false, false);
	    
	    AddMatrixRpsi(1.0, 0.0, 0.0, var_p.var_harmonic,
			  var_s.var_harmonic, var_elas.var_harmonic,
			  Kh, IndexDofDirichlet, IndexDofDirichlet, neumann_cond,
			  0, Np, true, false);
            
            // non-local part
            AddMatrixDtN(var_s.var_harmonic, var_elas.var_harmonic,
                         Sh, Kh, ListeDofNeu, IndexDofNeumann, Ngamma);
          }
          break;
	}
      
    }
  
  // calcul du second membre
  Vector<VectReal_wp> rhs_acous, rhs_elas;
  ComputeRightHandSide(var_s, var_p, var_elas, rhs_acous, rhs_elas, Ndof);
  
  DISP(Ndof);
  if (Ndof < 5000)
    {
      Mh.WriteText("Mh0.dat");
      Sh.WriteText("Sh0.dat");
      Kh.WriteText("Kh0.dat");
    }
  
  if (mean_value_phi_zero)
    {
      ModifyMatrixMean(Mh, ListeCouple, CoefCouple, removed_dof, 1.0);
      ModifyMatrixMean(Kh, ListeCouple, CoefCouple, removed_dof, 0.0);
      ModifyMatrixMean(Sh, ListeCouple, CoefCouple, removed_dof, 0.0);
      ModifyVectorMean(rhs_acous(0), ListeCouple, CoefCouple, removed_dof);
    }
  
  // on ecrit les matrices
  if (Ndof < 5000)
    {
      Mh.WriteText("Mh.dat");
      Sh.WriteText("Sh.dat");
      Kh.WriteText("Kh.dat");
    }
  
  VectReal_wp Phi_p(Np), Phi_s(Ns);
  
  if (type_resolution == ElastoParam::MODE_RESOL)
    {
      // we search eigenvalues and eigenvectors
      // of the first order system
      // dU/dt - V = 0
      // Mh dV/dt + Sh V + Kh U = 0
      // Mass dX/dt + rigid X = 0
      // with Mass = [I 0; 0 Mh]
      // and Rigid = [0 I; -Kh -Sh]
      
      if (type_formulation == ElastoParam::FOURIER_FORMUL)
        {
          var_fft.ComputeEigenvalues(var_p.var_harmonic);
          return;
        }

      int type_storage = 0;
      cout << "Voulez vous une resolution en plein ou en creux ?" << endl;
      cout << "1- Matrice pleine" << endl;
      cout << "2- Matrice creuse" << endl;
      cin >> type_storage;
      
      if (type_storage == 1)
        {
          Matrix<Real_wp> Mass(2*Ndof, 2*Ndof), Rigid(2*Ndof, 2*Ndof);
          Mass.SetIdentity(); Rigid.Fill(0);
          for (int i = 0; i < Ndof; i++)
            for (int j = 0; j < Mh.GetRowSize(i); j++)
              Mass(Ndof+i, Ndof+Mh.Index(i, j)) = Mh.Value(i, j);
          
          for (int i = 0; i < Ndof; i++)
            {
              Rigid(i, Ndof+i) = 1.0;
              
              for (int j = 0; j < Kh.GetRowSize(i); j++)
                Rigid(Ndof+i, Kh.Index(i, j)) = -Kh.Value(i, j);
              
              for (int j = 0; j < Sh.GetRowSize(i); j++)
                Rigid(Ndof+i, Ndof+Sh.Index(i, j)) = -Sh.Value(i, j);
            }
          
          //Mass.Write("Mass.dat");
          //Rigid.Write("Rigid.dat");
          
          // calcul de tous les valeurs propres et vecteurs propres
          VectReal_wp alpha_real, alpha_imag, beta;
          Matrix<Real_wp> eigen_vectors;
          GetEigenvaluesEigenvectors(Rigid, Mass, alpha_real, alpha_imag, beta, eigen_vectors);
          
          // on regroupe les valeurs propres dans un fichier
          // premiere colonne : partie reelle
          // seconde colonne : partie imaginaire
          //DISP(alpha_real); DISP(alpha_imag); DISP(beta);
          ofstream file_out("eigenvalue.dat"); file_out.precision(15);
          for (int i = 0; i < alpha_real.GetM(); i++)
            file_out << alpha_real(i)/beta(i) << " " << alpha_imag(i)/beta(i) << '\n';
          
          file_out.close();
          
          // on n'ecrit que les vecteurs propres "instables"
          Vector<string> nameP(10), nameS(10);
          for (int i = 0; i < var_p.var_harmonic.output_grid_param.GetM(); i++)
            nameP(i) = var_p.var_harmonic.output_grid_param(i).GetTotalFieldFile();
          
          for (int i = 0; i < var_s.var_harmonic.output_grid_param.GetM(); i++)
            nameS(i) = var_s.var_harmonic.output_grid_param(i).GetTotalFieldFile();
          
          int num = 0;
          for (int n = 0; n < alpha_real.GetM(); n++)
            {
              if (beta(n) != 0)
		if (abs(alpha_real(n)/beta(n)) > 1e-4)
		  // if (abs(alpha_real(n)) + abs(alpha_imag(n)) < 1e-5)
                  {
                    DISP(num);
                    DISP(alpha_real(n)/beta(n));
                    DISP(alpha_imag(n)/beta(n));
                    
                    for (int i = 0; i < Np; i++)
                      Phi_p(i) = eigen_vectors(i, n);
                    
                    for (int i = 0; i < Ns; i++)
                      Phi_s(i) = eigen_vectors(Np+i, n);
                    
                    if (mean_value_phi_zero)
                      {
                        RecomposeVectorMean(Phi_p, ListeCouple, CoefCouple, removed_dof);
                        RecomposeVectorMean(Phi_s, ListeCouple, CoefCouple, removed_dof);
                      }
                    
                    for (int i = 0; i < var_p.var_harmonic.output_grid_param.GetM(); i++)
                      {
                        string name = GetBaseString(nameP(i)) + "P" + to_str(num);
                        var_p.var_harmonic.output_grid_param(i).SetTotalFieldFile(name);
                      }
                    
                    for (int i = 0; i < var_s.var_harmonic.output_grid_param.GetM(); i++)
                      {
                        string name = GetBaseString(nameS(i)) + "S" + to_str(num);
                        var_s.var_harmonic.output_grid_param(i).SetTotalFieldFile(name);
                      }
                    
                    var_p.var_harmonic.WriteOutputFile(Phi_p, 1);
                    var_s.var_harmonic.WriteOutputFile(Phi_s, 1);
                    
                    num++;
                  }
            }
        }
      else
        {
          if (Sh.GetNonZeros() > 0)
            {
              cout << "Cas qui n'est pas pris en charge" << endl;
              cout << "Choisissez une formulation qui donne Sh = 0" << endl;
              cout << "Nombre d'elements non nuls : " << Sh.GetNonZeros() << endl;
              abort();
            }
          
          int nb_eigenval = 0;
          cout << "Nombre de valeurs propres a calculer ?" << endl;
          cin >> nb_eigenval;
          
          double shift_eigen;
          cout << "Les valeurs propres sont cherchees autour de shift " << endl;
          cout << "Shift a utiliser : " << endl;
          cin >> shift_eigen;

          SparseEigenProblem<double, Matrix<double, General, ArrayRowSparse>,
            Matrix<double, General, ArrayRowSparse> > var_eig;
          
          var_eig.InitMatrix(Kh, Mh);
          var_eig.SetStoppingCriterion(1e-12);
          var_eig.SetNbAskedEigenvalues(nb_eigenval);
          var_eig.SetComputationalMode(var_eig.INVERT_MODE);
          if (shift_eigen > 0)
            {
              //var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_REAL);
              var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, shift_eigen, var_eig.SORTED_REAL);
            }
          else
            var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, shift_eigen, var_eig.SORTED_REAL);
          
          Vector<double> lambda, lambda_imag;
          Matrix<double> eigen_vectors;
          GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vectors);
          
          // on ecrit les valeurs propres d'ordre 1 (on prend la racine carre)
          DISP(lambda); DISP(lambda_imag);
          ofstream file_out("eigenvalue.dat"); file_out.precision(15);
          for (int i = 0; i < lambda.GetM(); i++)
            {
              Real_wp Lr(0), Li(0);
              if (lambda(i) < 0)
                Lr = sqrt(-lambda(i));
              else
                Li = sqrt(lambda(i));
              
              file_out << Lr << " " << Li << '\n';
            }
          
          file_out.close();
          
          // on ecrit tous les vecteurs propres
          Vector<string> nameP(10), nameS(10);
          for (int i = 0; i < var_p.var_harmonic.output_grid_param.GetM(); i++)
            nameP(i) = var_p.var_harmonic.output_grid_param(i).GetTotalFieldFile();
          
          for (int i = 0; i < var_s.var_harmonic.output_grid_param.GetM(); i++)
            nameS(i) = var_s.var_harmonic.output_grid_param(i).GetTotalFieldFile();
          
          int num = 0;
          for (int n = 0; n < lambda.GetM(); n++)
            {
              DISP(num);
              Real_wp Lr(0), Li(0);
              if (lambda(n) < 0)
                Lr = sqrt(-lambda(n));
              else
                Li = sqrt(lambda(n));

              DISP(Lr); DISP(Li);
              
              for (int i = 0; i < Np; i++)
                Phi_p(i) = eigen_vectors(i, n);
              
              for (int i = 0; i < Ns; i++)
                Phi_s(i) = eigen_vectors(Np+i, n);
              
              if (mean_value_phi_zero)
                {
                  RecomposeVectorMean(Phi_p, ListeCouple, CoefCouple, removed_dof);
                        RecomposeVectorMean(Phi_s, ListeCouple, CoefCouple, removed_dof);
                }
              
              for (int i = 0; i < var_p.var_harmonic.output_grid_param.GetM(); i++)
                {
                  string entier = NumberToString(num);
                  string name = GetBaseString(nameP(i)) + "P" + entier;
                  var_p.var_harmonic.output_grid_param(i).SetTotalFieldFile(name);
                }
              
              for (int i = 0; i < var_s.var_harmonic.output_grid_param.GetM(); i++)
                {
                  string entier = NumberToString(num);
                  string name = GetBaseString(nameS(i)) + "S" + entier;
                  var_s.var_harmonic.output_grid_param(i).SetTotalFieldFile(name);
                }
              
              var_p.var_harmonic.WriteOutputFile(Phi_p, 1);
              var_s.var_harmonic.WriteOutputFile(Phi_s, 1);
              
              num++;
            }
        }
          
      return;
    }
  
  if (type_resolution == ElastoParam::HARMONIC_RESOL)
    {
      
      // on calcule Ah = -omega^2 Mh - i omega Sh + Kh
      DistributedMatrix<Complex_wp, General, ArrayRowSparse> Ah(Ndof, Ndof);
      
      Real_wp omega = var_s.var_harmonic.GetOmega();
      Real_wp omega2 = omega*omega;
      Vector<VectComplex_wp> sol_fft;      
      All_MatrixLU<Complex_wp> mat_lu;
      mat_lu.ShowMessages();
      VectComplex_wp xsol;
      
      if (type_formulation == ElastoParam::FOURIER_FORMUL)
        {
          GaussianSourceP<Real_wp> fsrc(var_p.var_harmonic);
          var_fft.ComputeSource(fsrc, sol_fft);          
          var_fft.ComputeSolution(var_p.var_harmonic, sol_fft);
        }
      else
        {
          IVect col_interac; int nb_interac; VectComplex_wp val_interac;
          for (int i = 0; i < Ndof; i++)
            {
              nb_interac = Kh.GetRowSize(i) + Mh.GetRowSize(i) + Sh.GetRowSize(i);
              col_interac.Reallocate(nb_interac);
              val_interac.Reallocate(nb_interac);
              nb_interac = 0;
              for (int j = 0; j < Kh.GetRowSize(i); j++)
                {
                  col_interac(nb_interac) = Kh.Index(i, j);
                  val_interac(nb_interac) = Kh.Value(i, j);
                  nb_interac++;
                }
              
              for (int j = 0; j < Sh.GetRowSize(i); j++)
                {
                  col_interac(nb_interac) = Sh.Index(i, j);
                  val_interac(nb_interac) = -Iwp*omega*Sh.Value(i, j);
                  nb_interac++;
                }
              
              for (int j = 0; j < Mh.GetRowSize(i); j++)
                {
                  col_interac(nb_interac) = Mh.Index(i, j);
                  val_interac(nb_interac) = -omega2*Mh.Value(i, j);
                  nb_interac++;
                }
              
              Ah.AddInteractionRow(i, nb_interac, col_interac, val_interac);
            }
      
          //Ah.WriteText("Ah.dat");
          
          // on calcule la solution xsol          
          mat_lu.Factorize(Ah);
          
          // printing error messages
          int ierr = 0;
          int type = mat_lu.GetInfoFactorization(ierr);
          PrintFactorizationFailed(type, ierr);
          
          xsol.Reallocate(Ndof);
          for (int i = 0; i < rhs_acous(0).GetM(); i++)
            xsol(i) = rhs_acous(0)(i);
        }
      
      // on calcule la matrice de l'elasto
      DistributedMatrix<Real_wp, General, ArrayRowSparse> Ah_elas;
      GlobalGenericMatrix<Real_wp> nat_mat;
      nat_mat.SetCoefMass(-var_s.var_harmonic.GetSquareOmega());
      nat_mat.SetCoefStiffness(1.0);
      
      var_elas.var_harmonic.AddMatrixWithBC(Ah_elas, nat_mat);
      
      All_MatrixLU<Real_wp> mat_lu_real; mat_lu_real.ShowMessages();      
      mat_lu_real.Factorize(Ah_elas);
      
      mat_lu_real.Solve(rhs_elas(0));      
      var_elas.var_harmonic.WriteDatas(rhs_elas(0));
      
      if (start_smooth_solution)
	{
	  // on resout le laplacien couple pour phi_p et phi_s avec les sources suivantes :
	  // - \int_\Omega div u varphi^p + \int_\Gamma v.n \varphi^p  
	  // \int_\Omega rot u varphi^s - \int_\Gamma n \times v \varphi^s  
	  abort();
	  /*Vector<Vector<VectReal_wp> > Uquad;
	  Vector<Vector<VectReal_wp> > gradU_quad;
	  var_elas.var_harmonic.ComputeQuadratureUgradU(rhs_elas, Uquad, gradU_quad, true, true, true);
	  
	  VolumetricSource<TypeEqElasStatic> source_elas(var_elas.var_harmonic);
	  source_elas.variable_source = true;
	  source_elas.evalS.Reallocate(2);
	  source_elas.evalSurf.Reallocate(2);
	  
	  int r = var_elas.var_harmonic.mesh_num.GetOrder();
	  int nb_elt = var_elas.var_harmonic.mesh.GetNbElt();
	  int nb_edges_ref = var_elas.var_harmonic.mesh.GetNbBoundaryRef();
	  source_elas.evalS(0).Reallocate(nb_elt, (r+1)*(r+1));
	  source_elas.evalS(1).Reallocate(nb_elt, (r+1)*(r+1));
	  
	  source_elas.evalSurf(0).Reallocate(nb_edges_ref, r+1);
	  source_elas.evalSurf(1).Reallocate(nb_edges_ref, r+1);
	  source_elas.evalSurf(0).Fill(0);
	  source_elas.evalSurf(1).Fill(0);
	  
	  for (int i = 0; i < nb_elt; i++)
	    for (int j = 0; j < (r+1)*(r+1); j++)
	      {
		source_elas.evalS(0)(i, j) = -(gradU_quad(i)(j)(0) + gradU_quad(i)(j)(3));
		source_elas.evalS(1)(i, j) = gradU_quad(i)(j)(2) - gradU_quad(i)(j)(1);
		//DISP(i); DISP(j); DISP(source_elas.evalS(1)(i, j));
	      }
	  
	  for (int i = 0; i < nb_elt; i++)
	    if (!var_elas.var_harmonic.GetFaceBasis(i).LumpedMassMatrix())
	      {
		cout << "Use Gauss-Lobatto elements for this functionnality" << endl;
		abort();
	      }
	  
	  for (int i = 0; i < nb_edges_ref; i++)
	    {
	      int num_elem = var_elas.var_harmonic.mesh.BoundaryRef(i).numElement(0);
	      int num_loc = var_elas.var_harmonic.mesh.Element(num_elem).GetPositionBoundary(i);
	      const ElementReference<Dimension2, 1>& Fb = var_elas.var_harmonic.GetReferenceElement(num_elem);
	      for (int j = 0; j <= r; j++)
		{
		  int npoint = Fb.GetLocalNumber(num_loc, j);
		  Real_wp Ux, Uy, Nx, Ny;
		  Ux = Uquad(num_elem)(npoint)(0);
		  Uy = Uquad(num_elem)(npoint)(1);
		  Nx = NormaleMesh(num_elem)(num_loc)(j)(0);
		  Ny = NormaleMesh(num_elem)(num_loc)(j)(1);
		  
		  source_elas.evalSurf(0)(i, j) = (Ux*Nx + Uy*Ny);
		  source_elas.evalSurf(1)(i, j) = -(Uy*Nx - Ux*Ny);
		}	      
	    }
	  
	  var_elas.var_harmonic.ComputeGenericSource(rhs_elas, source_elas);

	  if (mean_value_phi_zero)
	    ModifyVectorMean(rhs_elas(0), ListeCouple, CoefCouple, removed_dof);
	  
	  mat_lu_laplace.Solve(rhs_elas(0));
	  
	  for (int i = 0; i < Np; i++)
	    Phi_p(i) = real(rhs_elas(0)(i));
	  
	  for (int i = 0; i < Ns; i++)
	    Phi_s(i) = real(rhs_elas(0)(Np+i));

	  if (mean_value_phi_zero)
	    {
	      RecomposeVectorMean(Phi_p, ListeCouple, CoefCouple, removed_dof);
	      RecomposeVectorMean(Phi_s, ListeCouple, CoefCouple, removed_dof);
	      }*/
	}
      else
	{
	  if (type_formulation != ElastoParam::FOURIER_FORMUL)
            {            
              mat_lu.Solve(xsol);
              
              // on sort la solution
              for (int i = 0; i < Np; i++)
                Phi_p(i) = real(xsol(i));
              
              for (int i = 0; i < Ns; i++)
                Phi_s(i) = real(xsol(Np+i));
              
              if (mean_value_phi_zero)
                {
                  RecomposeVectorMean(Phi_p, ListeCouple, CoefCouple, removed_dof);
                  RecomposeVectorMean(Phi_s, ListeCouple, CoefCouple, removed_dof);
                }
            }
	}
      
      if (type_formulation == ElastoParam::FOURIER_FORMUL)
        {
          var_fft.WriteSolution(sol_fft);
        }
      else
        {
          var_s.var_harmonic.WriteDatas(Phi_s);
          var_p.var_harmonic.WriteDatas(Phi_p);
      
          if (Ndof >= Np+Ns+2*Ngamma)
            {
              Phi_p.Fill(0); Phi_s.Fill(0);
              for (int i = 0; i < Np; i++)
                {
                  int num = IndexDofNeumann(i);
                  if (num >= 0)
                    {
                      Phi_p(i) = real(xsol(Np+Ns+num));
                      Phi_s(i) = real(xsol(Np+Ns+Ngamma+num));
                    }
                }
              
              var_p.var_harmonic.WriteOutputFile(Phi_p, 0);
              var_s.var_harmonic.WriteOutputFile(Phi_s, 0);
            }
        }

      return;
    }
    
  // derniere possibilite : on calcule l'evolution en temps de u
  Real_wp dt = var_elas.GetTimeStep(), dt2 = dt*dt;
  Real_wp t_begin = var_elas.GetInitialTime();
  Real_wp t_end = var_elas.GetFinalTime();
  int nb_max_iter = toInteger(ceil(abs(t_end - t_begin)/dt ) );
  
  int Nelas = var_elas.var_harmonic.GetNbDof();
  VectReal_wp U0(Nelas), dU0(Nelas);
  ModifiedEquationIterator<Real_wp> scheme;
  U0.Fill(0); dU0.Fill(0);
  
  scheme.SetOrder(2);
  scheme.SetInitialConditionS(t_begin, dt, U0, dU0, var_elas);

  Real_wp theta = 0.25;
  cout << "Si vous prenez Delta t > 1, la CFL du schema semi-implicite est calculee " << endl;
  cout << "Entrez la valeur de theta (negatif pour avoir 1/4) " << endl;
  cin >> theta;
  
  if (theta < 0)
    theta = 0.25;
      
  DistributedMatrix<Real_wp, General, ArrayRowSparse> Ah;

  if (dt > Real_wp(1))
    {
      // we extract Bh part of Kh
      Ah.Reallocate(Ndof, Ndof);
      int nb_interac = 0;
      IVect col_interac; VectReal_wp val_interac;
      for (int i = 0; i < Np; i++)
	{
	  nb_interac = 0;
	  for (int j = 0; j < Kh.GetRowSize(i); j++)
	    if (Kh.Index(i, j) < Np)
	      nb_interac++;
	  
	  if (nb_interac > 0)
	    {
	      col_interac.Reallocate(nb_interac);
	      val_interac.Reallocate(nb_interac);
	      nb_interac = 0;
	      for (int j = 0; j < Kh.GetRowSize(i); j++)
		if (Kh.Index(i, j) < Np)
		  {
		    col_interac(nb_interac) = Kh.Index(i, j);
		    val_interac(nb_interac) = Kh.Value(i, j);
		    nb_interac++;
		  }
	      
	      Ah.AddInteractionRow(i, nb_interac, col_interac, val_interac);
	    }
	}

      for (int i = 0; i < Ns; i++)
	{
	  nb_interac = 0;
	  for (int j = 0; j < Kh.GetRowSize(Np+i); j++)
	    if (Kh.Index(Np+i, j) >= Np)
	      nb_interac++;
	  
	  if (nb_interac > 0)
	    {
	      col_interac.Reallocate(nb_interac);
	      val_interac.Reallocate(nb_interac);
	      nb_interac = 0;
	      for (int j = 0; j < Kh.GetRowSize(Np+i); j++)
		if (Kh.Index(Np+i, j) >= Np)
		  {
		    col_interac(nb_interac) = Kh.Index(Np+i, j);
		    val_interac(nb_interac) = Kh.Value(Np+i, j);
		    nb_interac++;
		  }
	      
	      Ah.AddInteractionRow(Np+i, nb_interac, col_interac, val_interac);
	    }
	}            
      
      // Ah.WriteText("Ah.dat");
      // Mh.WriteText("Mh.dat");
      
      // on calcule le rayon spectral de M^{-1/2} A M^{-1/2}
      Matrix<Real_wp, Symmetric, RowSymSparse> K;
      IVect Ptr(Ndof+1), Ind((Ah.GetNonZeros()+Ndof)/2);
      VectReal_wp Val(Ind.GetM()), scaleMh(Ndof);
      for (int i = 0; i < Ndof; i++)
	scaleMh(i) = 1.0/sqrt(Mh(i, i));
      
      Ptr(0) = 0; int nb = 0;
      for (int i = 0; i < Ndof; i++)
	{
	  for (int j = 0; j < Ah.GetRowSize(i); j++)
	    if (Ah.Index(i, j) >= i)
	      {
		Ind(nb) = Ah.Index(i, j);
		Val(nb) = Ah.Value(i, j)*scaleMh(i)*scaleMh(j);
		nb++;
	      }
	  
	  Ptr(i+1) = nb;
	}
      
      K.SetData(Ndof, Ndof, Val, Ptr, Ind);
      //K.WriteText("Ah_sym.dat");
      
      Real_wp rho = GetSpectralRadius(K);
      cout << "Spectral radius of M^{-1} A = " << rho << endl;
      
      Real_wp cfl = 1.0/sqrt(theta*rho);
      cout << "CFL semi-implicit = 1 / sqrt( theta rho) = " << cfl << endl;
      
      // on calcule le rayon spectral de M^{-1/2} K M^{-1/2}
      Ptr.Reallocate(Ndof+1); Ind.Reallocate((Kh.GetNonZeros()+Ndof)/2);
      Val.Reallocate(Ind.GetM());
      Ptr(0) = 0; nb = 0;
      for (int i = 0; i < Ndof; i++)
	{
	  for (int j = 0; j < Kh.GetRowSize(i); j++)
	    if (Kh.Index(i, j) >= i)
	      {
		Ind(nb) = Kh.Index(i, j);
		Val(nb) = Kh.Value(i, j)*scaleMh(i)*scaleMh(j);
		nb++;
	      }
	  
	  Ptr(i+1) = nb;
	}
      
      K.SetData(Ndof, Ndof, Val, Ptr, Ind);
      
      rho = GetSpectralRadius(K);
      cout << "Spectral radius of M^{-1} K = " << rho << endl;
      
      cfl = 2.0/sqrt(rho);
      cout << "CFL explicit = 2 / sqrt( rho) = " << cfl << endl;

      return;
    }

  bool fully_implicit_scheme = true;
  cout << "Si vous voulez un schema completement implicite, tapez 1" << endl;
  cin >> fully_implicit_scheme;
  
  if (presence_neumann)
    {
      if (!fully_implicit_scheme)
	{
	  cout << "With neumann condition, only fully-implicit scheme is implemented" << endl;
	  abort();
	}
    }
  
  if (var_elas.var_harmonic.GetFrequency() == Real_wp(0))
    {
      cout << "Entrez une frequence non nulle pour avoir une source differente de NaN" << endl;
      abort();
    }
  
  if (!fully_implicit_scheme)
    cout << "On utilise un schema implicite seulement sur Bh" << endl;
  else
    cout << "On utilise un schema totalement implicite" << endl;

  bool advance_elasto = true;
  cout << "Tapez un pour avancer l'elasto en meme temps " << endl;
  cin >> advance_elasto;

  Copy(Mh, Ah);
  if (fully_implicit_scheme)
    Add(theta*dt2, Kh, Ah);
  else
    {
      // we add only Bh part of Kh
      int nb_interac = 0;
      IVect col_interac; VectReal_wp val_interac;
      for (int i = 0; i < Np; i++)
	{
	  nb_interac = 0;
	  for (int j = 0; j < Kh.GetRowSize(i); j++)
	    if (Kh.Index(i, j) >= Np)
	      nb_interac++;
	  
	  if (nb_interac > 0)
	    {
	      col_interac.Reallocate(nb_interac);
	      val_interac.Reallocate(nb_interac);
	      nb_interac = 0;
	      for (int j = 0; j < Kh.GetRowSize(i); j++)
		if (Kh.Index(i, j) >= Np)
		  {
		    col_interac(nb_interac) = Kh.Index(i, j);
		    val_interac(nb_interac) = Kh.Value(i, j);
		    val_interac(nb_interac) *= theta*dt2;
		    nb_interac++;
		  }
	      
	      Ah.AddInteractionRow(i, nb_interac, col_interac, val_interac);
	    }
	}

      for (int i = 0; i < Ns; i++)
	{
	  nb_interac = 0;
	  for (int j = 0; j < Kh.GetRowSize(Np+i); j++)
	    if (Kh.Index(Np+i, j) < Np)
	      nb_interac++;
	  
	  if (nb_interac > 0)
	    {
	      col_interac.Reallocate(nb_interac);
	      val_interac.Reallocate(nb_interac);
	      nb_interac = 0;
	      for (int j = 0; j < Kh.GetRowSize(Np+i); j++)
		if (Kh.Index(Np+i, j) < Np)
		  {
		    col_interac(nb_interac) = Kh.Index(Np+i, j);
		    val_interac(nb_interac) = Kh.Value(Np+i, j);
		    val_interac(nb_interac) *= theta*dt2;
		    nb_interac++;
		  }
	      
	      Ah.AddInteractionRow(Np+i, nb_interac, col_interac, val_interac);
	    }
	}            
    }
  
  Add(0.5*dt, Sh, Ah);
  
  //Ah.WriteText("Ah.dat");
  All_MatrixLU<Real_wp> mat_lu; mat_lu.ShowMessages();
  mat_lu.Factorize(Ah);
  mat_lu.HideMessages();
    
  VectReal_wp Phi_next(Ndof);
  VectReal_wp Phi_n(Ndof), Phi_nm1(Ndof), DiffPhi(Ndof);
  Phi_next.Fill(0); Phi_n.Fill(0); Phi_nm1.Fill(0); DiffPhi.Fill(0);
  for (int nt = 0; nt < nb_max_iter; nt++)
    {
      Real_wp t = t_begin + (nt+1)*dt;
      if (advance_elasto)
	var_elas.GiveIterate(nt, t, scheme.GetIterate());
      
      var_p.WriteSnapshot(nt, t, Phi_p);
      var_s.WriteSnapshot(nt, t, Phi_s);
      
      Phi_next.Fill(0);
      var_s.AddPrimitiveSourceAtTime(dt2, t, 1, Phi_next);
      
      // subtracting dt^2 Kh U^n
      MltAdd(-dt2, Kh, Phi_n, 1.0, Phi_next);
      
      // subtracting dt Sh (Un - Un-1)
      for (int i = 0; i < Ndof; i++)
	DiffPhi(i) = Phi_n(i) - Phi_nm1(i);
      
      MltAdd(-dt, Sh, DiffPhi, 1.0, Phi_next);
      
      // then inverting by (Mh + dt Sh/2 + theta dt^2 Bh)
      mat_lu.Solve(Phi_next);
      
      // then adding 2 U^n - U^{n-1}
      for (int i = 0; i < Ndof; i++)
	Phi_next(i) += 2.0*Phi_n(i) - Phi_nm1(i);
      
      // next iterate
      Copy(Phi_n, Phi_nm1);
      Copy(Phi_next, Phi_n);
      
      // extracting phi_p and phi_s
      for (int i = 0; i < Np; i++)
        Phi_p(i) = Phi_n(i);

      for (int i = 0; i < Ns; i++)
        Phi_s(i) = Phi_n(Np+i);
      
      if (nt%10 == 1)
        DISP(Norm2(Phi_n));
      
      if (advance_elasto)
	scheme.Advance(t, nt, var_elas);
      
      //int test_input; cout << "we wait" << endl; cin >> test_input;
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 2)
    {
      cout << "Entrez le fichier de donnees a lire" << endl;
      abort();
    }

  string input_file(argv[1]);
  
  string type_element, type_equation;
  getElement_Equation(input_file, type_element, type_equation);

  cout << "0- Resolution en regime harmonique" << endl;
  cout << "1- Resolution en regime temporel" << endl;
  cout << "2- Calcul des modes" << endl;
  
  int type_resol(0);
  cin >> type_resol;
  
  //cout << "On utilise les points de Gauss-Legendre" << endl;
  cout << "Attention : on utilise les points de Gauss-Lobatto" << endl;
  cout << endl;
  
  cout << "0- Formulation avec (u, phi)" << endl;
  cout << "1- Formulation avec (v, phi)" << endl;
  cout << "2- Formulation avec (u, psi, phi)" << endl;
  cout << "3- Formulation avec (u, phi) et avec derivee seconde" << endl;
  cout << "4- Formulation avec (phi) en eliminant matriciellement u" << endl;
  cout << "5- Formulation avec (phi) avec la primitive de phi" << endl;
  cout << "6- Formulation avec l'operateur DtN" << endl;
  cout << "7- Formulation en faisant du Fourier en y" << endl;
  
  int type_formul(0);
  cin >> type_formul;
  
  HyperbolicProblem<TimeElasticEquation<Dimension2> > var_elas;
  HyperbolicProblem<AcousticEquation<Dimension2> > var_acous;

  RunAll(var_acous, var_elas, input_file, type_element,
	 type_resol, type_formul);
  
  return FinalizeMontjoie();
}
