#ifndef MONTJOIE_FILE_TRIANGLE_CLASSICAL_CXX

namespace Montjoie
{
  
#ifdef MONTJOIE_WITH_META_PROGRAMMING
  template<int r> 
  void OptTriangleClassical<r>::ConstructFiniteElement(const TriangleClassical& Fb)
  {
    for (int i = 0; i < Fb.GetNbDof(); i++)
      for (int j = i; j < Fb.GetNbDof(); j++)
	mat_mass(i, j) = Fb.mass_matrix(i, j);
    
    for (int i = 0; i < Fb.GetNbDof(); i++)
      for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
	{
	  mat_ch(i, j) = Fb.Value_Phi(i, j);
	  for (int k = 0; k < 2; k++)
	    mat_rigid(i, 2*j+k) = Fb.Gradient_Phi(i, j)(k);
	}
  }
#endif
  
  //! default constructor
  TriangleClassical::TriangleClassical() : TriangleReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }


  size_t TriangleClassical::GetMemorySize() const
  {
    size_t taille = TriangleReference<1>::GetMemorySize();
    taille += GL_G.GetMemorySize();
    taille += InverseBasisVDM.GetMemorySize();
    taille += InvWeightFct.GetMemorySize();
    taille += NumFct2D.GetMemorySize();
    taille += basis_phi1D.GetMemorySize();
    for (int i = 0; i < GradientPhi_Boundary.GetM(); i++)
      taille += GradientPhi_Boundary(i).GetMemorySize();
    
    return taille;
  }

  
  //! constructing finite element
  void TriangleClassical::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                                                 int rsurf, int type_surf)
  {
    TriangleReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad, rsurf, type_surf);
    
    ConstructFunctions();

    for (int num_loc = 0; num_loc < 3; num_loc++)
      {
        GradientPhi_Boundary(num_loc).Reallocate(nb_dof_loc, order+1);
        R2 pt_loc; VectR2 grad_phi;
        for (int k = 0; k <= order; k++)
          {
            this->GetLocalCoordOnBoundary(num_loc, this->PointsDof1D(k), pt_loc);
            this->ComputeGradientPhiRef(pt_loc, grad_phi);
            for (int i = 0; i < nb_dof_loc; i++)
              GradientPhi_Boundary(num_loc)(i, k) = grad_phi(i);            
          }
      }        
    
    GL_G.Reallocate(r+1, this->Points1D().GetM());
    for (int i = 0; i <= r; i++)
      for (int j = 0; j < this->Points1D().GetM(); j++)
        GL_G(i,j) = basis_phi1D.EvaluatePhi(i, this->Points1D(j));
    
    ConstructElementaryMatrix(*this);
    
#ifdef MONTJOIE_WITH_META_PROGRAMMING
    switch (r)
      {
      case 1:
	elt_r1.ConstructFiniteElement(*this);
	break;
      case 2:
	elt_r2.ConstructFiniteElement(*this);
	break;
      case 3:
	elt_r3.ConstructFiniteElement(*this);
	break;
      case 4:
	elt_r4.ConstructFiniteElement(*this);
	break;
      case 5:
	elt_r5.ConstructFiniteElement(*this);
	break;
      }
#endif

  }
  
  
  //! constructions of basis functions
  void TriangleClassical::ConstructFunctions()
  {
    int r = order;
    nb_dof_boundaries = 3*order;

    if (order == this->GetGeometryOrder())
      {
	this->elt_geom.dof_equal_nodal = true;
        EdgesDof = this->GetNodalNumber();
        
        nb_dof_loc = this->GetNbPointsNodalElt();
        this->SetPointsDof1D(this->PointsNodal1D());
	this->SetPointsDofND(this->PointsNodalND()); 

        InverseBasisVDM = Fb_geom.GetInverseVandermonde();
        NumFct2D = Fb_geom.GetNumOrtho2D();
        InvWeightFct = Fb_geom.GetInverseWeightPolynomial();
        
        basis_phi1D.AffectPoints(this->PointsDof1D());
      }
    else
      {
	this->elt_geom.dof_equal_nodal = false;

        Matrix<int> NumDofs2D, CoordinateDofs;
        // construction of nodes numbering
        MeshNumbering<Dimension2>::ConstructTriangularNumbering(r, NumDofs2D, CoordinateDofs);
    
        // EdgesNode(i, num_loc) : node number of i-th node of edge num_loc
        EdgesDof.Reallocate(r+1, 3);
        for (int i = 0; i <= r; i++)
          EdgesDof(i,0) = NumDofs2D(i,0);
        
        for (int i = 0; i <= r; i++)
          EdgesDof(i,1) = NumDofs2D(r-i, i);
        
        for (int i = 0; i <= r; i++)
          EdgesDof(i,2) = NumDofs2D(0, r-i);	
        
	VectReal_wp points_dof1d;
	VectR2 points_dof2d;
        switch (type_interpolation)
          {
          case TriangleGeomReference::REGULAR_BASIS :
            TriangleGeomReference::ConstructRegularPoints(r, points_dof1d, points_dof2d, NumDofs2D);
            break;
          case TriangleGeomReference::LOBATTO_BASIS :
            TriangleGeomReference::ConstructLobattoPoints(r, type_interpolation, points_dof1d, points_dof2d);
            break;
          case TriangleGeomReference::MASS_LUMPED_BASIS :
            {
              VectReal_wp omega;
              ComputeGaussLobatto(points_dof1d, omega, r);
              TriangleQuadrature::ConstructQuadrature(r, points_dof2d, omega,
                                                      TriangleQuadrature::QUADRATURE_MASS_LUMPED);
            }
            break;
          }    
        
	this->SetPointsDof1D(points_dof1d);
	this->SetPointsDofND(points_dof2d);
	
        basis_phi1D.AffectPoints(points_dof1d);
        
        nb_dof_loc = points_dof2d.GetM();
        if (nb_dof_loc == (r+1)*(r+2)/2)
          {
            InvWeightFct.Reallocate((r+1)*(r+2)/2);
            NumFct2D.Reallocate(r+1, r+1); NumFct2D.Fill(-1);
            int nb = 0;
	    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefficientLegendre();
	    const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefficientOddJacobi();
            for (int sum = 0; sum <= r; sum++)
              for (int i = 0; i <= sum; i++)
                {
                  int j = sum-i;
                  NumFct2D(i, j) = nb;
                  InvWeightFct(nb) = CoefLegendre(i)*CoefOddJacobi(i, j);
                  nb++;
                }
            
            InverseBasisVDM.Reallocate(nb_dof_loc, nb_dof_loc);
            InverseBasisVDM.Fill(Real_wp(0));
            
            // evaluation of orthogonal polynomials
            VectReal_wp phi;
            for (int k = 0; k < nb_dof_loc; k++)
              {
                Fb_geom.ComputeValuesPhiOrthoRef(r, NumFct2D, InvWeightFct,
                                                 points_dof2d(k), phi);
                
                for (int j = 0; j < nb_dof_loc; j++)
                  InverseBasisVDM(j, k) = phi(j);
              }
            
            GetInverse(InverseBasisVDM);
          }
	else
	  {
	    cout << "Impossible case" << endl;
	    abort();
	  }
      }
    
    this->nb_points_dof_inside = this->nb_dof_loc;
    this->num_dof_points_surf.Reallocate(3);
    for (int n = 0; n < 3; n++)
      {
	this->num_dof_points_surf(n).Reallocate(r+1);
	for (int i = 0; i <= r; i++)
	  this->num_dof_points_surf(n)(i) = EdgesDof(i, n);
      }

    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);        
    Fb_geom.ComputeCoefficientTransformation();
  }
  
  
  //! changing interpolation points
  void TriangleClassical::SetInterpolationPoints(const VectR2& points)
  {
    int r = order;
    this->SetPointsDofND(points);
    if (points.GetM() != (r+1)*(r+2)/2)
      {
	cout << "Case not treated" << endl;
	abort();
      }
    else
      {
	InvWeightFct.Reallocate((r+1)*(r+2)/2);
	NumFct2D.Reallocate(r+1, r+1); NumFct2D.Fill(-1);
	int nb = 0;
	const VectReal_wp& CoefLegendre = Fb_geom.GetCoefficientLegendre();
	const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefficientOddJacobi();
	for (int sum = 0; sum <= r; sum++)
	  for (int i = 0; i <= sum; i++)
	    {
              int j = sum-i;
	      NumFct2D(i, j) = nb;
	      InvWeightFct(nb) = CoefLegendre(i)*CoefOddJacobi(i, j);
	      nb++;
	    }
	
	InverseBasisVDM.Reallocate(nb_dof_loc, nb_dof_loc);
	InverseBasisVDM.Fill(Real_wp(0));
	
	// evaluation of orthogonal polynomials
	VectReal_wp phi;
	for (int k = 0; k < nb_dof_loc; k++)
	  {
	    Fb_geom.ComputeValuesPhiOrthoRef(r, NumFct2D, InvWeightFct,
					  points(k), phi);
	    
	    for (int j = 0; j < nb_dof_loc; j++)
	      InverseBasisVDM(j, k) = phi(j);
	  }
	
	GetInverse(InverseBasisVDM);
      }
	
    for (int num_loc = 0; num_loc < 3; num_loc++)
      {
        GradientPhi_Boundary(num_loc).Reallocate(nb_dof_loc, order+1);
        R2 pt_loc; VectR2 grad_phi;
        for (int k = 0; k <= order; k++)
          {
            this->GetLocalCoordOnBoundary(num_loc, this->PointsDof1D(k), pt_loc);
            this->ComputeGradientPhiRef(pt_loc, grad_phi);
            for (int i = 0; i < nb_dof_loc; i++)
              GradientPhi_Boundary(num_loc)(i, k) = grad_phi(i);            
          }
      }        
    
    GL_G.Reallocate(r+1, this->Points1D().GetM());
    for (int i = 0; i <= r; i++)
      for (int j = 0; j < this->Points1D().GetM(); j++)
        GL_G(i,j) = basis_phi1D.EvaluatePhi(i, this->Points1D(j));
    
    ConstructElementaryMatrix(*this);
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of basis functions on pointloc
  */
  void TriangleClassical::ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& phi) const
  {
    VectReal_wp psi;
    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    Fb_geom.ComputeValuesPhiOrthoRef(order, NumFct2D, InvWeightFct, pointloc, psi);
    Mlt(InverseBasisVDM, psi, phi);
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on pointloc
  */
  void TriangleClassical::ComputeGradientPhiRef(const R2& pointloc, VectR2& grad_phi) const
  {
    VectReal_wp psi, phi; VectR2 grad_psi; 
    psi.Reallocate(nb_dof_loc); psi.Fill(0);
    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    grad_phi.Reallocate(nb_dof_loc);
    Fb_geom.ComputeGradientPhiOrthoRef(order, NumFct2D, InvWeightFct, pointloc, grad_psi);
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
  }
  
      
  //! select dofs near the current dof
  void TriangleClassical::PickNearDofs(int pos, const VectBool& DofUsed,
				       IVect& ListeDof, int nb_dof) const
  {
    if (nb_dof <= 0)
      return;
    
    const Matrix<int>& CoordinateNodes = this->GetCoordinateNodes2D();
    int i0 = CoordinateNodes(pos,0); int j0 = CoordinateNodes(pos,1);
    ListeDof.Reallocate(nb_dof); ListeDof.Fill(-1);
    int node, k = 1, nb = 0;
    if (!DofUsed(pos))
      ListeDof(nb++) = pos;
    
    while (nb < nb_dof)
      {
	// loop on concentric square at distance k
	if ((j0-k) >= 0)
	  for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	    if (nb < nb_dof)
	      {
		node = NumFct2D(m,j0-k);
		if (!DofUsed(node))
		  ListeDof(nb++) = node;
	      }
	
	if (((i0+k) <= order)&&(nb < nb_dof))
	  for (int m = max(0,j0-k+1); m < min(order,j0+k); m++)
	    if (nb < nb_dof)
	      {
		node = NumFct2D(i0+k,m);
		if (!DofUsed(node))
		  ListeDof(nb++) = node;
	      }
	
	if (((j0+k) <= order)&&(nb < nb_dof))
	  for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	    if (nb < nb_dof)
	      {
		node = NumFct2D(m,j0+k);
		if (!DofUsed(node))
		  ListeDof(nb++) = node;
	      }
	
	if (((i0-k) >= 0)&&(nb < nb_dof))
	  for (int m = max(0,j0-k+1); m < min(order,j0+k); m++)
	    if (nb < nb_dof)
	      {
		node = NumFct2D(i0-k,m);
		if (!DofUsed(node))
		  ListeDof(nb++) = node;
	      }
	
	k++;
	if (k > order)
	  nb = nb_dof;
	
      }
	
  }

  
  //! computation of u on quadrature points of an edge
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the edge
    \param[in] r order of quadrature rule of the edge
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
  */
  template<class Vector1, class Vector2>
  void TriangleClassical
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    int num_dof; Vh.Fill(0);
    if ((r == 0) || (r == order_quad))
      {
	for (int i = 0; i <= order; i++)
	  {
	    num_dof = EdgesDof(i, num_loc);
	    for (int j = 0; j <= order_quad; j++)
	      Vh(j) += GL_G(i, j)*Uh(num_dof);
	  }
      }
    else
      {
	for (int i = 0; i <= order; i++)
	  {
	    num_dof = EdgesDof(i, num_loc);
	    for (int j = 0; j <= r; j++)
	      Vh(j) += ProjOperatorOrder(r)(i, j)*Uh(num_dof);
	  }
      }
  }
  
  
  //! integration against basis functions on a face
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc face number
    \param[in] Uquad vector containing \omega_k f(\xi_k) 
    \param[out] Vh res_i = res_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the edge
  */
  template<class T0, class Vector1, class Vector2>
  void TriangleClassical::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh,
				     Vector2& Vh, int r) const
  {
    int num_dof;
    typename Vector2::value_type vloc;
    if ((r == 0) || (r == order_quad))
      {
	for (int i = 0; i <= order; i++)
	  {
	    num_dof = EdgesDof(i, num_loc);
	    vloc = 0;
	    for (int j = 0; j <= order_quad; j++)
	      vloc += GL_G(i, j)*Uh(j);

	    Vh(num_dof) += alpha*vloc;
	  }
      }
    else
      {
	for (int i = 0; i <= order; i++)
	  {
	    num_dof = EdgesDof(i, num_loc);
	    vloc = 0;
	    for (int j = 0; j <= r; j++)
	      vloc += ProjOperatorOrder(r)(i, j)*Uh(j);
	    
	    Vh(num_dof) += alpha*vloc;
	  }
      }
  }
  
  
  //! computation of gradient of u on quadrature points of an edge
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points of the edge
    \param[in] r order of quadrature rule of the edge
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
  */
  template<class Vector1, class Vector2>
  void TriangleClassical::ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh,
						   Vector2& Vh, int r) const
  {
    Vh.Fill(0);
    if ((r == 0) || (r == order_quad))
      ElementReference<Dimension2, 1>::ApplyNablaShTranspose(num_loc, Uh, Vh, r);
    else
      {
        typedef typename Vector2::value_type T;
        TinyVector<T, 2> vloc; Real_wp val;
	T zero; SetComplexZero(zero);
	for (int i = 0; i <= order; i++)
	  {
            vloc.Fill(zero);
            for (int j = 0; j < nb_dof_loc; j++)
              Add(Uh(j), GradientPhi_Boundary(num_loc)(j, i), vloc);
            
	    for (int j = 0; j <= r; j++)
	      {
                val = ProjOperatorOrder(r)(i, j);
                Vh(2*j) += val*vloc(0);
                Vh(2*j+1) += val*vloc(1);
              }
	  }
      }
  }
  
  
  //! integration against gradient of basis functions on an edge
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc edge number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is the quadrature rule of edge num_loc
   */
  template<class T0, class Vector1, class Vector2>
  void TriangleClassical::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh,
					  Vector2& Vh, int r) const
  {
    if ((r == 0) || (r == order_quad))
      ElementReference<Dimension2, 1>::ApplyNablaSh(alpha, num_loc, Uh, Vh, r);
    else
      {
        typedef typename Vector2::value_type T;
	typename Vector2::value_type zero; SetComplexZero(zero);
        TinyVector<T, 2> vloc; Real_wp val;
	for (int i = 0; i <= order; i++)
	  {
            vloc.Fill(zero);
            for (int j = 0; j <= r; j++)
	      {
                val = ProjOperatorOrder(r)(i, j);
                vloc(0) += val*Uh(2*j);
                vloc(1) += val*Uh(2*j+1);
              }
            
            Mlt(alpha, vloc);
            for (int j = 0; j < nb_dof_loc; j++)
              Vh(j) += DotProd(vloc, GradientPhi_Boundary(num_loc)(j, i));
            	    
	  }
      }
  }


  //! displays informations about class TriangleReference
  ostream& operator <<(ostream& out, const TriangleClassical& e)
  {
    out<<static_cast<const TriangleReference<1>&>(e);
    out<<"Number of degrees of freedom on a face "<<e.nb_dof_loc<<endl;
    return out;
  }
  
  
  //! default constructor
  TriangleLobatto::TriangleLobatto() : TriangleClassical()
  {
    type_integration_edge = Globatto<Real_wp>::QUADRATURE_LOBATTO;
  }
  

  size_t TriangleLobatto::GetMemorySize() const
  {
    size_t taille = TriangleClassical::GetMemorySize();
    taille += GL_Gquad.GetMemorySize();
    return taille;
  }
  
  
  //! constructing finite element
  void TriangleLobatto::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					       int rsurf, int type_surf)
  {    
    // we use Gauss-Lobatto points for edges
    TriangleClassical::ConstructFiniteElement(r, rgeom, rquad, type_quad,
                                              rsurf, Globatto<Real_wp>::QUADRATURE_LOBATTO);

    if (rquad == 0)
      rquad = r;

    if (rsurf == 0)
      rsurf = rquad;
    
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(rsurf);
    GL_Gquad.Reallocate(r+1, rsurf+1);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= rsurf; j++)
        GL_Gquad(i, j) = basis_phi1D.EvaluatePhi(i, gauss.Points(j));
    
  }
  
  
  //! Integration against basis functions on an edge
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc edge number
  */  
  template<class Vector1,class Vector2>
  void TriangleLobatto::
  ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
                                 Vector2& res, int num_loc) const
  {
    typename Vector2::value_type zero; SetComplexZero(zero);
    res.Fill(zero);
    for (int i = 0; i <= order; i++)
      {
        int num_dof = EdgesDof(i, num_loc);
        typename Vector1::value_type vloc(zero);
        for (int j = 0; j < GL_Gquad.GetN(); j++)
          vloc += GL_Gquad(i, j)*feval(j);
        
        res(num_dof) = vloc;
      }
  }
  
}
  
#define MONTJOIE_FILE_TRIANGLE_CLASSICAL_CXX
#endif
