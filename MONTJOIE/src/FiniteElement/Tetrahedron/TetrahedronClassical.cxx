#ifndef MONTJOIE_FILE_TETRAHEDRON_CLASSICAL_CXX

#ifdef MONTJOIE_WITH_META_PROGRAMMING
//#include "OptTetrahedronClassical.cxx"
#endif

namespace Montjoie
{

#ifdef MONTJOIE_WITH_META_PROGRAMMING
  template<int r> 
  void OptTetrahedronClassical<r>::ConstructFiniteElement(const TetrahedronClassical& Fb)
  {
    for (int i = 0; i < Fb.GetNbDof(); i++)
      for (int j = i; j < Fb.GetNbDof(); j++)
	mat_mass(i, j) = Fb.mass_matrix(i, j);
    
    for (int i = 0; i < Fb.GetNbDof(); i++)
      for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
	{
	  mat_ch(i, j) = Fb.Value_Phi(i, j);
	  for (int k = 0; k < 3; k++)
	    mat_rigid(i, 3*j+k) = Fb.Gradient_Phi(i, j)(k);
	}
  }
#endif
  
  //! default constructor
  TetrahedronClassical::TetrahedronClassical() : TetrahedronReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }


  size_t TetrahedronClassical::GetMemorySize() const
  {
    size_t taille = TetrahedronReference<1>::GetMemorySize();
    taille += NumFct3D.GetMemorySize();
    taille += InvWeightFct.GetMemorySize();
    taille += InverseBasisVDM.GetMemorySize();
    for (int k = 0; k < GradientPhi_Boundary.GetM(); k++)
      taille += GradientPhi_Boundary(k).GetMemorySize();
    
    return taille;
  }
  
  
  //! constructing finite element
  void TetrahedronClassical::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                         int rsurf_tri, int rsurf_quad,
                         int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    TetrahedronReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    ConstructFunctions();
    this->ConstructElementaryMatrix(*this);

    this->Fb_geom.ComputeCoefficientTransformation();

    const VectR2& points_dof2d_tri = this->PointsDof2D_tri();
    for (int num_loc = 0; num_loc < 4; num_loc++)
      {
        GradientPhi_Boundary(num_loc).Reallocate(nb_dof_loc, function_basis_tri->GetNbDof());
        R3 pt_loc; VectR3 grad_phi;
        for (int k = 0; k < function_basis_tri->GetNbDof(); k++)
          {
            this->GetLocalCoordOnBoundary(num_loc, points_dof2d_tri(k), pt_loc);
            this->ComputeGradientPhiRef(pt_loc, grad_phi);
            for (int i = 0; i < nb_dof_loc; i++)
              GradientPhi_Boundary(num_loc)(i, k) = grad_phi(i);            
          }
      }        
    
#ifdef MONTJOIE_WITH_META_PROGRAMMING
    /*switch (r)
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
      }*/
#endif

  }


  //! construction of basis functions
  void TetrahedronClassical::ConstructFunctions()
  {
    int order_geom = this->GetGeometryOrder();
    if (order == order_geom)
      {
	this->elt_geom.dof_equal_nodal = true;
        FacesDof = this->Fb_geom.GetNodalNumber();

        this->SetPointsDof2D_tri(this->PointsNodal2D_tri());
        this->SetPointsDof2D_quad(this->PointsNodal2D_quad());
        this->SetPointsDofND(this->Fb_geom.PointsNodalND());
	
        nb_dof_loc = this->GetNbPointsNodalElt();
        NumFct3D = Fb_geom.GetNumOrtho3D();
        InvWeightFct = Fb_geom.GetInvWeightPolynomial();        
      }
    else
      {
	this->elt_geom.dof_equal_nodal = false;
	
        int r = order;
        Array3D<int> NumDofs3D;
        Matrix<int> CoordinateDofs, NumDofs2D_tri;
        // construction of nodes numbering
        MeshNumbering<Dimension2>::
          ConstructTriangularNumbering(r, NumDofs2D_tri, CoordinateDofs);
        
        MeshNumbering<Dimension3>::ConstructTetrahedralNumbering(r, NumDofs3D, CoordinateDofs);

        FacesDof.Reallocate((r+1)*(r+2)/2, 4);
        // nodes on faces
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            {
              int node = NumDofs2D_tri(i,j);
              FacesDof(node,0) = NumDofs3D(i, j, 0);
              FacesDof(node,1) = NumDofs3D(i, 0, j);
              FacesDof(node,2) = NumDofs3D(0, i, j);
              FacesDof(node,3) = NumDofs3D(r-i-j, i, j);
            }
        
	VectReal_wp points_dof1d; VectR2 points_dof2d_tri;
	VectR3 points_dof3d;
        switch (type_interpolation)
          {
          case TetrahedronGeomReference::REGULAR_BASIS :
            TetrahedronGeomReference::
	      ConstructRegularPoints(r, points_dof1d, points_dof2d_tri, NumDofs2D_tri,
				     points_dof3d, NumDofs3D);
            break;
          case TetrahedronGeomReference::LOBATTO_BASIS :
	    TetrahedronGeomReference::
	      ConstructLobattoPoints(r, points_dof1d, points_dof2d_tri, points_dof3d);
            break;
          }    
        
	this->SetPointsDof1D(points_dof1d);
	this->SetPointsDof2D_tri(points_dof2d_tri);
	this->SetPointsDofND(points_dof3d);
	
	const VectReal_wp& CoefLegendre = Fb_geom.GetCoefLegendre();
	const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefOddJacobi();
	const Matrix<Real_wp>& CoefEvenJacobi = Fb_geom.GetCoefEvenJacobi();
        nb_dof_loc = points_dof3d.GetM();
        if (nb_dof_loc == (r+1)*(r+2)*(r+3)/6)
          {
            InvWeightFct.Reallocate((r+1)*(r+2)*(r+3)/6);
            NumFct3D.Reallocate(r+1, r+1, r+1); NumFct3D.Fill(-1);
            int nb = 0;
            for (int i = 0; i <= r; i++)
              for (int j = 0; j <= r-i; j++)
                for (int k = 0; k <= r-i-j; k++)
                  {
                    NumFct3D(i, j, k) = nb;
                    InvWeightFct(nb) = CoefLegendre(i)*CoefOddJacobi(i, j)*CoefEvenJacobi(i+j, k);
                    nb++;
                  }
            
            InverseBasisVDM.Reallocate(nb_dof_loc, nb_dof_loc);
            InverseBasisVDM.Fill(Real_wp(0));
            
            // evaluation of orthogonal polynomials
            VectReal_wp phi;
            for (int k = 0; k < nb_dof_loc; k++)
              {
                Fb_geom.ComputeValuesPhiOrthoRef(r, NumFct3D, InvWeightFct, points_dof3d(k), phi);
                for (int j = 0; j < nb_dof_loc; j++)
                  InverseBasisVDM(j, k) = phi(j);
              }
            
            GetInverse(InverseBasisVDM);
          }
      }

    this->nb_points_dof_inside = this->nb_dof_loc;
    this->num_dof_points_surf.Reallocate(4);
    for (int n = 0; n < 4; n++)
      {
	this->num_dof_points_surf(n).Reallocate(FacesDof.GetM());
	for (int i = 0; i < FacesDof.GetM(); i++)
	  this->num_dof_points_surf(n)(i) = FacesDof(i, n);
      }
    
    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);
    
    TriangleClassical* tri = new TriangleClassical();
    tri->ConstructFiniteElement(order, order_geom, order_quad);
    function_basis_tri = tri;
    element_tri_surf = tri;
    
    this->nb_dof_tri = function_basis_tri->GetNbDof();
    this->nb_dof_quad = 0;
    
    nb_dof_boundaries = 2*order*order+2;
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] x local point where functions are evaluated
    \param[out] phi values of basis functions for x
  */
  void TetrahedronClassical::ComputeValuesPhiRef(const R3& x, VectReal_wp& phi) const
  {
    if (this->GetGeometryOrder() == order)
      ComputeValuesPhiNodalRef(x, phi);
    else
      {
	VectReal_wp psi; phi.Reallocate(nb_dof_loc); phi.Fill(0);
        Fb_geom.ComputeValuesPhiOrthoRef(order, NumFct3D, InvWeightFct, x, psi);
        Mlt(InverseBasisVDM, psi, phi);
      }
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] x local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions for x
  */
  void TetrahedronClassical::ComputeGradientPhiRef(const R3& x, VectR3& grad_phi) const
  {
    if (this->GetGeometryOrder() == order)
      ComputeGradientPhiNodalRef(x, grad_phi);
    else
      {
        VectReal_wp psi, phi; VectR3 grad_psi; 
        psi.Reallocate(nb_dof_loc); psi.Fill(0);
        phi.Reallocate(nb_dof_loc); phi.Fill(0);
        grad_phi.Reallocate(nb_dof_loc);
        Fb_geom.ComputeGradientPhiOrthoRef(order, NumFct3D, InvWeightFct, x, grad_psi);
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
          {
            grad_phi(i)(1) = phi(i);
            psi(i) = grad_psi(i)(2);
          }
        
        Mlt(InverseBasisVDM, psi, phi);
        for (int i = 0; i < psi.GetM(); i++)
          grad_phi(i)(2) = phi(i);
        
      }
  }
  
  
  //! selects dofs near the current dof
  /*!
    \param[in] pos current dof
    \param[in] DofUsed some dofs are already in use and can't be selected
    \param[out] ListeDof selected dofs
    \param[in] nb_dof number of dofs to select
   */
  void TetrahedronClassical::PickNearDofs(int pos, const VectBool& DofUsed,
					 IVect& ListeDof, int nb_dof) const
  {
    if (nb_dof <= 0)
      return;
    
    const Array3D<int>& NumNodes3D = this->Fb_geom.GetNumNodes3D();
    const Matrix<int>& CoordinateNodes = this->Fb_geom.GetCoordinateNodes3D();
    int i0 = CoordinateNodes(pos,0); int j0 = CoordinateNodes(pos,1), k0 = CoordinateNodes(pos,2);
    ListeDof.Reallocate(nb_dof); ListeDof.Fill(-1);
    int k = 1, nb = 0, node;
    if (!DofUsed(pos))
      ListeDof(nb++) = pos;
    
    while (nb < nb_dof)
      {
	// loop on concentric cube at distance k
	if ((j0-k) >= 0)
	  for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	    for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumNodes3D(m,j0-k,n);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	if (((j0+k) <= order)&&(nb < nb_dof))
	  for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	    for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumNodes3D(m,j0+k,n);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	if (((i0-k) >= 0)&&(nb < nb_dof))
	  for (int m = max(0,j0-k+1); m < min(order,j0+k); m++)
	    for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumNodes3D(i0-k,m,n);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	if (((i0+k) <= order)&&(nb < nb_dof))
	  for (int m = max(0,j0-k+1); m < min(order,j0+k); m++)
	    for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumNodes3D(i0+k,m,n);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	if (((k0-k) >= 0)&&(nb < nb_dof))
	  for (int m = max(0,i0-k+1); m < min(order,i0+k); m++)
	    for (int n = max(0,j0-k+1); n < min(order,j0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumNodes3D(m,n,k0-k);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	if (((k0+k) <= order)&&(nb < nb_dof))
	  for (int m = max(0,i0-k+1); m < min(order,i0+k); m++)
	    for (int n = max(0,j0-k+1); n < min(order,j0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumNodes3D(m,n,k0+k);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	
	k++;
	if (k > order)
	  nb = nb_dof;
	
      }
	
  }
  
  
  //! computation of u on quadrature points of a face
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the face
    \param[in] r order of quadrature rule of the face
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void TetrahedronClassical
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0)||(r == order_quad))
      {
	Vector1 Udof(this->nb_dof_tri);
	for (int i = 0; i < this->nb_dof_tri; i++)
	  Udof(i) = Uh(FacesDof(i, num_loc));
	
	function_basis_tri->ApplyChTranspose(Udof, Vh);
      }
    else
      {
	Vh.Fill(0);
	const Matrix<Real_wp>& Ah = ProjOperatorTriOrder(r);
	typename Vector1::value_type vloc;
	for (int i = 0; i < this->nb_dof_tri; i++)
	  {
	    vloc = Uh(FacesDof(i, num_loc));
	    for (int j = 0; j < Ah.GetN(); j++)
	      Vh(j) += Ah(i, j)*vloc;
	  }
      }
  }
   
    
  //! integration against basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void TetrahedronClassical
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0)||(r == order_quad))
      {
	Vector1 Udof(this->nb_dof_tri);
	function_basis_tri->ApplyCh(Uh, Udof);
	
	for (int i = 0; i < this->nb_dof_tri; i++)
	  Vh(FacesDof(i, num_loc)) += alpha*Udof(i);
      }
    else
      {
	const Matrix<Real_wp>& Ah = ProjOperatorTriOrder(r);
	typename Vector1::value_type vloc;
	for (int i = 0; i < this->nb_dof_tri; i++)
	  {
	    vloc = 0;
	    for (int j = 0; j < Ah.GetN(); j++)
	      vloc += Ah(i, j)*Uh(j);
	    
	    Vh(FacesDof(i, num_loc)) += alpha*vloc;
	  }
      }	
  }
  

  //! computation of gradient of u on quadrature points of a boundary
  /*!
    \param[in] num_loc boundary number
    \param[in] res components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* res
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void TetrahedronClassical
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0)||(r == order_quad))
      ElementReference<Dimension3, 1>::ApplyNablaShTranspose(num_loc, Uh, Vh, r);
    else
      {
	Vh.Fill(0);
        typedef typename Vector2::value_type T;
        TinyVector<T, 3> vloc; Real_wp val;
	const Matrix<Real_wp>& Ah = ProjOperatorTriOrder(r);
	for (int i = 0; i < function_basis_tri->GetNbDof(); i++)
	  {
            vloc.Fill(0);
            for (int j = 0; j < nb_dof_loc; j++)
              Add(Uh(j), GradientPhi_Boundary(num_loc)(j, i), vloc);
            
	    for (int j = 0; j < Ah.GetN(); j++)
	      {
                val = Ah(i, j);
                Vh(3*j) += val*vloc(0);
                Vh(3*j+1) += val*vloc(1);
                Vh(3*j+2) += val*vloc(2);
              }
	  }
      }
  }
   
    
  //! integration against gradient of basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh  vector containing \omega_k f(\xi_k) 
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void TetrahedronClassical
  ::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0)||(r == order_quad))
      ElementReference<Dimension3, 1>::ApplyNablaSh(alpha, num_loc, Uh, Vh, r);
    else
      {
        typedef typename Vector2::value_type T;
        TinyVector<T, 3> vloc; Real_wp val;
	const Matrix<Real_wp>& Ah = ProjOperatorTriOrder(r);
        for (int i = 0; i < function_basis_tri->GetNbDof(); i++)
	  {
            vloc.Fill(0);
            for (int j = 0; j < Ah.GetN(); j++)
	      {
                val = Ah(i, j);
                vloc(0) += val*Uh(3*j);
                vloc(1) += val*Uh(3*j+1);
                vloc(2) += val*Uh(3*j+2);
              }
            
            Mlt(alpha, vloc);
            for (int j = 0; j < nb_dof_loc; j++)
              Vh(j) += DotProd(vloc, GradientPhi_Boundary(num_loc)(j, i));
            	    
	  }
      }	
  }

  
  //! displays details about class TetrahedronClassical
  ostream& operator <<(ostream& out, const TetrahedronClassical& e)
  {
    out<<static_cast<const TetrahedronReference<1>&>(e);
    out<<"Number dof over the tetrahedron "<<e.nb_dof_loc<<endl;
    return out;
  }
  
} // end namespace
  
#define MONTJOIE_FILE_TETRAHEDRON_CLASSICAL_CXX
#endif
