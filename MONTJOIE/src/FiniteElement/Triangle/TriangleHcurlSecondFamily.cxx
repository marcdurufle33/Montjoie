#ifndef MONTJOIE_FILE_TRIANGLE_HCURL_SECOND_FAMILY_CXX

namespace Montjoie
{
  //! default constructor
  TriangleHcurlSecondFamily::TriangleHcurlSecondFamily() : TriangleReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  
  
  //! how to number a mesh
  void TriangleHcurlSecondFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TriangleReference<2>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order+1);
    nmap.SetNbDofQuadrangle(order, 0);
    nmap.SetNbDofTriangle(order, nb_dof_loc - 3*(order+1));
    
    nmap.SetOppositeEdgesDofSymmetry(order, order+1);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
  }
  
  
  //! constructing finite element
  void TriangleHcurlSecondFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf, int type_surf)
  {
    TriangleReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    ConstructFunctions();

    this->Fb_geom.ComputeCoefficientTransformation();
    ConstructHcurlElementaryMatrix();
  }


  size_t TriangleHcurlSecondFamily::GetMemorySize() const
  {
    size_t taille = TriangleReference<2>::GetMemorySize();
    taille += NodalDof.GetMemorySize();
    taille += DirectionDof.GetMemorySize();
    taille += normale_dof.GetMemorySize();
    taille += ListeDof_Node.GetMemorySize();
    taille += NumDofs2D.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    taille += AhDof.GetMemorySize();
    taille += InverseBasisVDM.GetMemorySize();
    taille += InvWeightFct.GetMemorySize();
    taille += NumFct2D.GetMemorySize();
    return taille;
  }

  
  //! construction of basis functions
  void TriangleHcurlSecondFamily::ConstructFunctions()
  {
    nb_dof_boundaries = 3*(order+1);
    
    // two dofs by nodal point
    nb_dof_loc = (order+2)*(order+1);
    
    // construction of nodes numbering
    int r = order;
    MeshNumbering<Dimension2>::ConstructTriangularNumbering(r, NumDofs2D, CoordinateDofs);
    
    VectReal_wp points_dof1d; VectR2 points_dof2d;
    Fb_geom.ConstructLobattoPoints(order, 0, points_dof1d, points_dof2d);
    
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefficientLegendre();
    const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefficientOddJacobi();
    
    InvWeightFct.Reallocate((r+1)*(r+2)/2);
    NumFct2D.Reallocate(r+1, r+1);
    NumFct2D.Fill(-1);
    int nb = 0;
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        {
          NumFct2D(i, j) = nb;
          InvWeightFct(nb) = CoefLegendre(i)*CoefOddJacobi(i, j);
          nb++;
        }
    
    InverseBasisVDM.Reallocate(points_dof2d.GetM(), points_dof2d.GetM());
    InverseBasisVDM.Fill(Real_wp(0));
    
    // evaluation of orthogonal polynomials
    VectReal_wp phi;
    for (int k = 0; k < points_dof2d.GetM(); k++)
      {
        Fb_geom.ComputeValuesPhiOrthoRef(r, NumFct2D, InvWeightFct, points_dof2d(k), phi);
        for (int j = 0; j < nb_dof_loc/2; j++)
          InverseBasisVDM(j, k) = phi(j);
      }
    
    GetInverse(InverseBasisVDM);
    
    EdgesDof.Reallocate(order+1, 3);
    NodalDof.Reallocate(nb_dof_loc);
    DirectionDof.Reallocate(nb_dof_loc);
    normale_dof.Reallocate(nb_dof_loc);
    
    // First Edge
    for (int i = 0; i <= order; i++)
      {
	EdgesDof(i, 0) = i;
	if (i == order)
	  {
	    DirectionDof(EdgesDof(i, 0)).Init(1.0,1.0);
	    normale_dof(EdgesDof(i, 0)).Init(1.0, -1.0);
	  }
	else
	  {
	    DirectionDof(EdgesDof(i, 0)).Init(1.0,0.0);
	    normale_dof(EdgesDof(i, 0)).Init(0.0, -1.0);
	  }
	
	NodalDof(EdgesDof(i, 0)) = NumDofs2D(i,0);
      }
    
    // Second Edge
    Real_wp coef = 0.5;
    for (int i = 0; i <= order; i++)
      {
	EdgesDof(i, 1) = i+order+1;
	if (i == 0)
	  {
	    DirectionDof(EdgesDof(i, 1)).Init(0.0,1.0);
	    normale_dof(EdgesDof(i, 1)).Init(1.0, 0.0);
	  }
	else if (i == order)
	  {
	    DirectionDof(EdgesDof(i, 1)).Init(-1.0,0.0);
	    normale_dof(EdgesDof(i, 1)).Init(0.0, 1.0);
	  }
	else
	  {
	    DirectionDof(EdgesDof(i, 1)).Init(-coef,coef);
	    normale_dof(EdgesDof(i, 1)).Init(coef, coef);
	  }
	
	NodalDof(EdgesDof(i, 1)) = NumDofs2D(order-i,i);
      }
    
    // Third Edge
    for (int i = 0; i <= order; i++)
      {
	EdgesDof(i, 2) = i+2*(order+1);
	if (i == 0)
	  {
	    DirectionDof(EdgesDof(i, 2)).Init(-1.0,-1.0);
	    normale_dof(EdgesDof(i, 2)).Init(-1.0, 1.0);
	  }
	else if (i == order)
	  {
	    DirectionDof(EdgesDof(i, 2)).Init(0.0,-1.0);
	    normale_dof(EdgesDof(i, 2)).Init(-1.0, 0.0);
	  }
	else
	  {
	    DirectionDof(EdgesDof(i, 2)).Init(0.0,-1.0);
	    normale_dof(EdgesDof(i, 2)).Init(-1.0, 0.0);
	  }
	
	NodalDof(EdgesDof(i, 2)) = NumDofs2D(0,order-i);
      }
    
    int numero = nb_dof_boundaries;
    // edge-based internal dofs
    for (int i = 1; i < order; i++)
      {
	DirectionDof(numero).Init(0.0,1.0);
	NodalDof(numero) = NumDofs2D(i,0);
	normale_dof(numero).Init(Real_wp(1), Real_wp(0));
	numero++;
      }
    
    for (int i = 1; i < order; i++)
      {
	DirectionDof(numero).Init(-coef,-coef);
	NodalDof(numero) = NumDofs2D(order-i,i);
	normale_dof(numero).Init(-coef, coef);
	numero++;
      }
    
    for (int i = 1; i < order; i++)
      {
	DirectionDof(numero).Init(1.0,0.0);
	NodalDof(numero) = NumDofs2D(0,order-i);
	normale_dof(numero).Init(0.0, -1.0);
	numero++;
      }
    
    // inside dofs
    int nb_points_inside = (order-2)*(order-1)/2;
    for (int i = 0; i < nb_points_inside; i++)
      {
	DirectionDof(numero).Init(1.0,0.0);
	NodalDof(numero) = i + points_dof2d.GetM() - nb_points_inside;
	normale_dof(numero).Init(0.0, -1.0);
	numero++;
	DirectionDof(numero).Init(0.0,1.0);
	NodalDof(numero) = i + points_dof2d.GetM() - nb_points_inside;
	normale_dof(numero).Init(1.0, 0.0);
	numero++;
      }

    this->nb_points_dof_inside = points_dof2d.GetM();
    this->num_dof_points_surf.Reallocate(3);
    for (int n = 0; n < 3; n++)
      this->num_dof_points_surf(n).Reallocate(order+1);
    
    this->num_dof_points_surf(0)(0) = 0; this->num_dof_points_surf(0)(order) = 1;
    this->num_dof_points_surf(1)(0) = 1; this->num_dof_points_surf(1)(order) = 2;
    this->num_dof_points_surf(2)(0) = 2; this->num_dof_points_surf(2)(order) = 0;
    for (int i = 1; i < order; i++)
      {
	this->num_dof_points_surf(0)(i) = 2+i;
	this->num_dof_points_surf(1)(i) = 1+order+i;
	this->num_dof_points_surf(2)(i) = 2*order+i;
      }
    
    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);        
    ListeDof_Node.Reallocate(points_dof2d.GetM(),2);
    IVect nb_dof_node(points_dof2d.GetM()); nb_dof_node.Zero();
    for (int i = 0; i < nb_dof_loc; i++)
      {
	int node = NodalDof(i);
	ListeDof_Node(node, nb_dof_node(node)) = i;
	nb_dof_node(node)++;
      }
    
    AhDof.Reallocate(points_dof2d.GetM());
    Matrix2_2 A_tmp;
    for (int i = 0; i < points_dof2d.GetM(); i++)
      {
	A_tmp(0,0) = DirectionDof(ListeDof_Node(i,0))(0);
	A_tmp(1,0) = DirectionDof(ListeDof_Node(i,0))(1);
	A_tmp(0,1) = DirectionDof(ListeDof_Node(i,1))(0);
	A_tmp(1,1) = DirectionDof(ListeDof_Node(i,1))(1);
	GetInverse(A_tmp, AhDof(i));
      }    

    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);
  }
  

  //! projection of feval on dofs
  /*!
    \param[in] feval evaluation of a function on dof points
    \param[out] contrib dof components
   */
  template<class Vector1, class Vector2>
  void TriangleHcurlSecondFamily::ComputeProjectionDofGen(const Vector1& feval,
                                                          Vector2& contrib) const
  {
    typedef typename Vector2::value_type Complexe;
    TinyVector<Complexe,2> Eloc, Enodal;
    for (int node = 0; node < this->PointsDofND().GetM(); node++)
      {
	Eloc(0) = feval(2*node);
        Eloc(1) = feval(2*node+1);
	Mlt(AhDof(node), Eloc, Enodal);
	contrib(ListeDof_Node(node,0)) = Enodal(0);
	contrib(ListeDof_Node(node,1)) = Enodal(1);
      }
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void TriangleHcurlSecondFamily::
  ComputeValuesPhiRef(const R2& point_loc, VectR2& res) const
  {
    res.Reallocate(nb_dof_loc);
    VectReal_wp psi, phi(this->PointsDofND().GetM());
    Fb_geom.ComputeValuesPhiOrthoRef(order, NumFct2D, InvWeightFct, point_loc, psi);
    
    Mlt(InverseBasisVDM, psi, phi);
    
    for (int num_dof = 0; num_dof < nb_dof_loc; num_dof++)
      {
	res(num_dof) = DirectionDof(num_dof);
	Mlt(phi(NodalDof(num_dof)), res(num_dof));
      }
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void TriangleHcurlSecondFamily::
  ComputeCurlPhiRef(const R2& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    
    VectR2 grad_phi;
    VectReal_wp psi, phi; VectR2 grad_psi; 
    psi.Reallocate(nb_dof_loc/2); psi.Fill(0);
    phi.Reallocate(nb_dof_loc/2); phi.Fill(0);
    grad_phi.Reallocate(nb_dof_loc/2);
    Fb_geom.ComputeGradientPhiOrthoRef(order, NumFct2D, InvWeightFct, point_loc, grad_psi);
    for (int i = 0; i < psi.GetM(); i++)
      psi(i) = grad_psi(i)(0);
    
    Mlt(InverseBasisVDM, psi, phi);
    
    for (int i = 0; i < psi.GetM(); i++)
      {
        grad_phi(i)(0) = phi(i);
        psi(i) = grad_psi(i)(1);
      }
    
    Mlt(InverseBasisVDM, psi, phi);
    for (int i = 0; i < psi.GetM(); i++)
      grad_phi(i)(1) = phi(i);
    
    for (int num_dof = 0; num_dof < nb_dof_loc; num_dof++)
      {
	R2 grad = grad_phi(NodalDof(num_dof));
	res(num_dof) = DirectionDof(num_dof)(1)*grad(0)-DirectionDof(num_dof)(0)*grad(1);
      }
  }
  
  
  //! displays informations about class TriangleHcurlSecondFamily
  ostream& operator <<(ostream& out, const TriangleHcurlSecondFamily& e)
  {
    out<<static_cast<const TriangleReference<2>&>(e);
    out<<"Triangle of Nedelec's first family "<<endl;
    out<<"Number of local degrees of freedom "<<e.GetNbDof()<<endl;
    return out;
  }

}
  
#define MONTJOIE_FILE_TRIANGLE_HCURL_SECOND_FAMILY_CXX
#endif
