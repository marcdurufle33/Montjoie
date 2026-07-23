#ifndef MONTJOIE_FILE_ELEMENT_REFERENCE_CXX

#include "ElementGeomReference.cxx"

#include "ElementReference.hxx"
#include "ElementReferenceInline.cxx"
namespace Montjoie
{
  
  /****************
   * H^1 elements *
   ****************/
  
#ifdef MONTJOIE_WITH_NODAL_H1  

  //! projection from quadrature points to dof components
  /*!
    \param[in] Fb leaf finite element class
    \param[in] Equad values of E on quadrature points
    \param[out] Edof components of E on degrees of freedom
  */
  template<class Dimension> template<class T1>
  void ElementReferenceType<Dimension, 1>::
  ProjectQuadratureToDofGen(const ElementReference<Dimension, 1>& Fb,
			    const Vector<T1>& Equad, Vector<T1>& Edof) const
  {
    Vector<T1> Eweight(Fb.GetNbPointsQuadratureInside());
    for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
      Eweight(j) = Fb.WeightsND(j)*Equad(j);
    
    Fb.ApplyCh(Eweight, Edof);
    Fb.SolveMassMatrix(Edof);
  }


  //! returns the size taken by the object
  template<class Dimension>
  size_t ElementReferenceType<Dimension, 1>::GetMemorySize() const
  {
    size_t size = sizeof(Real_wp)*Dimension::dim_N*GradPhi_Nodal.GetDataSize();
    size += sizeof(Real_wp)*ValuePhi_Nodal.GetDataSize();
    size += sizeof(Real_wp)*Value_Phi.GetDataSize();
    size += sizeof(Real_wp)*Dimension::dim_N*Gradient_Phi.GetDataSize();
    size += sizeof(Real_wp)*const_grad_matrix.GetDataSize();
    size += sizeof(Real_wp)*const_stiff_matrix.GetDataSize();
    return size;
  }
  

  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */
  template<class Dimension>
  void ElementReferenceType<Dimension, 1>::
  GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const
  {
    int nb_dof_loc = Value_Phi.GetM();
    phi.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      phi(i) = Value_Phi(i, k);
  }  


  //! Retrieving gradient of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] grad_phi gradient of basis functions on quadrature point k
   */
  template<class Dimension>
  void ElementReferenceType<Dimension, 1>::
  GetGradientPhiOnQuadraturePoint(int k, VectR_N& grad_phi) const
  {
    int nb_dof_loc = Gradient_Phi.GetM();
    grad_phi.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      grad_phi(i) = Gradient_Phi(i, k);
  }

  
  //! Fills the value of a single basis function on quadrature points
  template<class Dimension>
  void ElementReferenceType<Dimension, 1>
  ::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    for (int k = 0; k < phi.GetM(); k++)
      phi(k) = Value_Phi(i, k);    
  }
  
  
  //! Fills the value and gradient of a single basis function on quadrature points
  template<class Dimension>
  void ElementReferenceType<Dimension, 1>
  ::GetGradientSinglePhiQuadrature(int i, VectReal_wp& phi, VectR_N& grad_phi) const
  {
    for (int k = 0; k < phi.GetM(); k++)
      {
	phi(k) = Value_Phi(i, k);    
	grad_phi(k) = Gradient_Phi(i, k);
      }
  }
  

  //! computation of mass matrix, stiffness matrix, etc (for continuous elements)
  /*!
    \param[in] Fb leaf finite element class
   */
  template<class Dimension>
  void ElementReferenceType<Dimension, 1>::ConstructElementaryMatrix(ElementReference<Dimension, 1>& Fb)
  {
    int nb_dof_loc = Fb.GetNbDof();
    int nb_points_quadrature = Fb.GetNbPointsQuadrature();
    VectReal_wp phi;
    Value_Phi.Reallocate(nb_dof_loc, nb_points_quadrature);        
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	Fb.ComputeValuesPhiRef(Fb.PointsND(i), phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  Value_Phi(j, i) = phi(j);
      }

    ValuePhi_Nodal.Reallocate(nb_dof_loc, Fb.GetNbPointsNodalElt());
    for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
      {
	Fb.ComputeValuesPhiRef(Fb.PointsNodalND(i), phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  ValuePhi_Nodal(j, i) = phi(j);
      }
    
    Vector<R_N> grad_phi;
    GradPhi_Nodal.Reallocate(nb_dof_loc, Fb.GetNbPointsNodalElt());
    for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
      {
	Fb.ComputeGradientPhiRef(Fb.PointsNodalND(i), grad_phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  GradPhi_Nodal(j, i) = grad_phi(j);
      }
    
    Fb.mass_matrix.Reallocate(nb_dof_loc, nb_dof_loc);
    Real_wp vloc;
    
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = i; j < nb_dof_loc; j++)
	{
	  vloc = 0.0;
	  for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
            vloc += Fb.WeightsND(k)*Value_Phi(i, k)*Value_Phi(j, k);
          
	  Fb.mass_matrix(i,j) = vloc;
	}

    Fb.mass_matrix_chol = Fb.mass_matrix;
    GetCholesky(Fb.mass_matrix_chol);
    
    
    Gradient_Phi.Reallocate(nb_dof_loc, nb_points_quadrature);
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	Fb.ComputeGradientPhiRef(Fb.PointsND(i), grad_phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  Gradient_Phi(j, i) = grad_phi(j);
      }

    Matrix<Real_wp, General, ArrayRowSparse> Sh;
    int offset = Fb.GetNbPointsQuadratureInside();
    for (int n = 0; n < Fb.GetNbBoundaries(); n++)
      {
        int Nquad = Fb.GetNbQuadBoundary(n);
        Sh.Clear();
        Sh.Reallocate(nb_dof_loc, Nquad);
        for (int i = 0; i < nb_dof_loc; i++)
          for (int k = 0; k < Nquad; k++)
            {
              int kloc = Fb.GetQuadNumber(n, k);
              if (abs(Value_Phi(i, kloc)) > 20.0*epsilon_machine)
                Sh.AddInteraction(i, k, Value_Phi(i, kloc));
            }

        Copy(Sh, Fb.sparse_const_sh(n));

        Fb.const_nabla_sh(n).Reallocate(nb_dof_loc, Dimension::dim_N*Nquad);
        for (int i = 0; i < nb_dof_loc; i++)
          for (int k = 0; k < Nquad; k++)
            {
              int kloc = Fb.GetQuadNumber(n, k);
              for (int m = 0; m < Dimension::dim_N; m++)
                Fb.const_nabla_sh(n)(i, Dimension::dim_N*k + m) = Gradient_Phi(i, kloc)(m);
            }
        
        offset += Nquad;
      }
    
    const_stiff_matrix.Reallocate(Dimension::dim_N*nb_dof_loc, Dimension::dim_N*nb_dof_loc);
    const_grad_matrix.Reallocate(nb_dof_loc, Dimension::dim_N*nb_dof_loc);
        
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	{
          for (int m = 0; m < Dimension::dim_N; m++)
            {
              vloc = 0;
              for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
                vloc += Fb.WeightsND(k)*Gradient_Phi(i, k)(m)*Value_Phi(j, k);
              
              const_grad_matrix(i, Dimension::dim_N*j+m) = vloc;
	    }
          
          for (int m = 0; m < Dimension::dim_N; m++)
            for (int n = 0; n < Dimension::dim_N; n++)
              {
                vloc = 0;
                for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
                  vloc += Fb.WeightsND(k)*Gradient_Phi(i, k)(m)*Gradient_Phi(j, k)(n);
                
                const_stiff_matrix(Dimension::dim_N*i + m, Dimension::dim_N*j + n) = vloc;
              }
        }
    
  }   
 
#endif
  

  /******************
   * H(curl) in 2-D *
   ******************/

#ifdef MONTJOIE_WITH_NODAL_HCURL
      
  //! projection from quadrature points to dof components
  /*!
    \param[in] Fb leaf finite element class
    \param[in] Equad values of E on quadrature points
    \param[out] Edof components of E on degrees of freedom
   */
  template<class T1>
  void ElementReferenceType<Dimension2, 2>::
  ProjectQuadratureToDofGen(const ElementReference<Dimension2, 2>& Fb,
			    const Vector<T1>& Equad, Vector<T1>& Edof) const
  {
    int nb_points_quadrature_inside = Fb.GetNbPointsQuadratureInside();
    Vector<T1> Eweight(2*nb_points_quadrature_inside);
    int p = 0;
    for (int j = 0; j < nb_points_quadrature_inside; j++)
      {
        Eweight(p) = Fb.WeightsND(j)*Equad(p);
        Eweight(p+1) = Fb.WeightsND(j)*Equad(p+1);
        p += 2;
      }
    
    Fb.ApplyCh(Eweight, Edof);
    Fb.SolveMassMatrix(Edof);
  }


  //! return the size taken by the object
  size_t ElementReferenceType<Dimension2, 2>::GetMemorySize() const
  {
    size_t size = sizeof(Real_wp)*Curl_Phi2D.GetDataSize();
    size += sizeof(Real_wp)*ValuePhi_Nodal.GetDataSize();
    size += sizeof(Real_wp)*CurlPhi_Nodal.GetDataSize();
    size += const_mass_matrix(0, 0).GetMemorySize();
    size += const_mass_matrix(1, 0).GetMemorySize();
    size += const_mass_matrix(0, 1).GetMemorySize();
    size += const_mass_matrix(1, 1).GetMemorySize();
    size += const_curl_matrix.GetMemorySize();
    size += Value_PhiVec.GetMemorySize();
    size += tangente_dof.GetMemorySize();
    return size;
  }
  

  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */
  void ElementReferenceType<Dimension2, 2>::
  GetValuePhiOnQuadraturePoint(int k, VectR2& phi) const
  {
    int nb_dof_loc = Value_PhiVec.GetM();
    phi.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      phi(i) = Value_PhiVec(i, k);
  }
  

  //! Retrieving curl of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] grad_phi gradient of basis functions on quadrature point k
   */
  void ElementReferenceType<Dimension2, 2>::
  GetCurlPhiOnQuadraturePoint(int k, VectReal_wp& grad_phi) const
  {
    int nb_dof_loc = Curl_Phi2D.GetM();
    grad_phi.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      grad_phi(i) = Curl_Phi2D(i, k);
  }
  
  
  /******************
   * H(curl) in 3-D *
   ******************/
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! projection from quadrature points to dof components
  /*!
    \param[in] Fb leaf finite element class
    \param[in] Equad values of E on quadrature points
    \param[out] Edof components of E on degrees of freedom
  */
  template<class T1>
  void ElementReferenceType<Dimension3, 2>::
  ProjectQuadratureToDofGen(const ElementReference<Dimension3, 2>& Fb,
			    const Vector<T1>& Equad, Vector<T1>& Edof) const
  {
    int nb_points_quadrature_inside = Fb.GetNbPointsQuadratureInside();
    Vector<T1> Eweight(3*nb_points_quadrature_inside);
    int p = 0;
    for (int j = 0; j < nb_points_quadrature_inside; j++)
      {
        Eweight(p) = Fb.WeightsND(j)*Equad(p);
        Eweight(p+1) = Fb.WeightsND(j)*Equad(p+1);
        Eweight(p+2) = Fb.WeightsND(j)*Equad(p+2);
        p += 3;
      }
    
    Fb.ComputeIntegralRef(Eweight, Edof);
    Fb.SolveMassMatrix(Edof);
  }


  //! returns the size taken by the object
  size_t ElementReferenceType<Dimension3, 2>::GetMemorySize() const
  {
    size_t size = 3*sizeof(Real_wp)*Curl_Phi.GetDataSize();
    size += sizeof(Real_wp)*ValuePhi_Nodal.GetDataSize();
    size += sizeof(Real_wp)*CurlPhi_Nodal.GetDataSize();
    size += 3*sizeof(Real_wp)*Value_PhiVec.GetDataSize();
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	{
	  size += sizeof(Real_wp)*const_mass_matrix(i, j).GetDataSize();
	  size += sizeof(Real_wp)*const_curl_matrix(i, j).GetDataSize();
	}
    
    size += tangente_dof.GetMemorySize();
    return size;
  }
  

  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */
  void ElementReferenceType<Dimension3, 2>::
  GetValuePhiOnQuadraturePoint(int k, VectR3& phi) const
  {
    int nb_dof_loc = Value_PhiVec.GetM();
    phi.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      phi(i) = Value_PhiVec(i, k);
  }


  //! Retrieving curl of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] grad_phi gradient of basis functions on quadrature point k
   */  
  void ElementReferenceType<Dimension3, 2>::
  GetCurlPhiOnQuadraturePoint(int k, VectR3& grad_phi) const
  {
    int nb_dof_loc = Curl_Phi.GetM();
    grad_phi.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      grad_phi(i) = Curl_Phi(i, k);
  }
#endif
  
#endif

  /**********
   * H(div) *
   **********/

#ifdef MONTJOIE_WITH_NODAL_HDIV

  //! projection from quadrature points to dof components
  /*!
    \param[in] Fb leaf finite element class
    \param[in] Equad values of E on quadrature points
    \param[out] Edof components of E on degrees of freedom
   */
  template<> template<class T1>
  void ElementReferenceType<Dimension2, 3>::
  ProjectQuadratureToDofGen(const ElementReference<Dimension2, 3>& Fb,
			    const Vector<T1>& Equad, Vector<T1>& Edof) const
  {
    int nb_points_quadrature_inside = Fb.GetNbPointsQuadratureInside();
    Vector<T1> Eweight(2*nb_points_quadrature_inside);
    int p = 0;
    for (int j = 0; j < nb_points_quadrature_inside; j++)
      {
        Eweight(p) = Fb.WeightsND(j)*Equad(p);
        Eweight(p+1) = Fb.WeightsND(j)*Equad(p+1);
	p += 2;
      }
    
    Fb.ApplyCh(Eweight, Edof);
    Fb.SolveMassMatrix(Edof);
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! projection from quadrature points to dof components
  /*!
    \param[in] Fb leaf finite element class
    \param[in] Equad values of E on quadrature points
    \param[out] Edof components of E on degrees of freedom
   */
  template<> template<class T1>
  void ElementReferenceType<Dimension3, 3>::
  ProjectQuadratureToDofGen(const ElementReference<Dimension3, 3>& Fb,
			    const Vector<T1>& Equad, Vector<T1>& Edof) const
  {
    int nb_points_quadrature_inside = Fb.GetNbPointsQuadratureInside();
    Vector<T1> Eweight(3*nb_points_quadrature_inside);
    int p = 0;
    for (int j = 0; j < nb_points_quadrature_inside; j++)
      {
        Eweight(p) = Fb.WeightsND(j)*Equad(p);
        Eweight(p+1) = Fb.WeightsND(j)*Equad(p+1);
        Eweight(p+2) = Fb.WeightsND(j)*Equad(p+2);
        p += 3;
      }
    
    Fb.ComputeIntegralRef(Eweight, Edof);
    Fb.SolveMassMatrix(Edof);
  }
#endif


  //! returns the size take by the object
  template<class Dimension>
  size_t ElementReferenceType<Dimension, 3>::GetMemorySize() const
  {
    size_t size = sizeof(Real_wp)*DivPhi_Nodal.GetDataSize();
    size += sizeof(Real_wp)*ValuePhi_Nodal.GetDataSize();
    size += sizeof(Real_wp)*Div_Phi.GetDataSize();
    size += sizeof(Real_wp)*const_div_matrix.GetDataSize();
    size += Dimension::dim_N*sizeof(Real_wp)*Value_PhiVec.GetDataSize();
    size += tangente_dof.GetMemorySize();
    for (int i = 0; i < Dimension::dim_N; i++)
      for (int j = 0; j < Dimension::dim_N; j++)
        size += const_mass_matrix(i, j).GetMemorySize();
    
    return size;
  }

  
  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */
  template<class Dimension>
  void ElementReferenceType<Dimension, 3>::
  GetValuePhiOnQuadraturePoint(int k, VectR_N& phi) const
  {
    int nb_dof_loc = Value_PhiVec.GetM();
    phi.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      phi(i) = Value_PhiVec(i, k);
  }


  //! Retrieving divergence of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] div_phi gradient of basis functions on quadrature point k
   */
  template<class Dimension>
  void ElementReferenceType<Dimension, 3>::
  GetDivPhiOnQuadraturePoint(int k, VectReal_wp& div_phi) const
  {
    int nb_dof_loc = Div_Phi.GetM();
    div_phi.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      div_phi(i) = Div_Phi(i, k);
  }  
#endif

  
  /************************
   * ElementReference_Dim *
   ************************/
  

  //! default constructor
  template<class Dimension>
  ElementReference_Dim<Dimension>::ElementReference_Dim(ElementGeomReference<Dimension>& elt)
    : ElementReference_Base(elt.WeightsND()), ElementGeomReferenceContainer<Dimension>(elt)
  {
  }
  

  //! returns quadrature points associated with 3-D integrals only
  template<class Dimension>
  typename Dimension::VectR_N ElementReference_Dim<Dimension>::PointsQuadInsideND() const
  {
    typename Dimension::VectR_N pts_quad(this->nb_points_quadrature_inside);
    for (int i = 0; i < this->nb_points_quadrature_inside; i++)
      pts_quad(i) = this->PointsND(i);
    
    return pts_quad;
  }


  //! computation of functions associated with quadrature points on a point
  /*!
    \param[in] point local point where functions are evaluated
    \param[out] phi values of functions on point
    From quadrature points, it is usually possible to construct 
    Lagrange interpolation functions satisfying phi_i(\xi_j) = \detla_{i, j}
    These functions are associated with quadrature points
   */
  template<class Dimension>
  void ElementReference_Dim<Dimension>
  ::ComputeValuesPhiQuadratureRef(const R_N& point, VectReal_wp& phi) const
  {
    abort();
  }
  
  
  //! computes res = Fi(point) for points located on the edges
  /*!
    \param[in] s list of vertices of the triangle
    \param res references points after transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad number of the triangle in the mesh
    \param[in] ne local edge number
  */
  template<class Dimension>
  void ElementReference_Dim<Dimension>::
  FjSurfaceElem(const VectR_N& s, SetPoints<Dimension>& res,
		const Mesh<Dimension>& mesh, int nquad, int ne) const
  {
    int nb_points = this->GetNbQuadBoundary(ne);
    res.ReallocatePointsQuadratureBoundary(nb_points);
    for (int m = 0; m < nb_points; m++)
      {
	int num_point = this->num_quad_points_surf(ne)(m);
	res.SetPointQuadratureBoundary(m, res.GetPointQuadrature(num_point));
      }
  }
  
  
  //! computes res = DFi(point) for points located on the edges
  /*!
    \param[in] s list of vertices of the triangle
    \param[in] PTReel references points after transformation Fi
    \param res jacobian matrices
    \param[in] mesh given mesh
    \param[in] nquad triangle number in the mesh
    \param[in] ne local edge number
    this method computes the normal and surfacic element (ds) on each quadrature point
  */
  template<class Dimension>
  void ElementReference_Dim<Dimension>::
  DFjSurfaceElem(const VectR_N& s, const SetPoints<Dimension>& PTReel,
		 SetMatrices<Dimension>& res,
		 const Mesh<Dimension>& mesh, int nquad, int ne) const
  {
    
    R_N normale_fj;
    Real_wp dsj(0); 
    int nb_points = this->GetNbQuadBoundary(ne);    
    res.ReallocatePointsQuadratureBoundary(nb_points);
    MatrixN_N dfjm1, mat_dfj;
    R_N point_loc;
    for (int k = 0; k < nb_points; k++)
      {
	int num_point = this->num_quad_points_surf(ne)(k);
	res.SetPointQuadratureBoundary(k, res.GetPointQuadrature(num_point));
	mat_dfj = res.GetPointQuadrature(num_point);
	Real_wp deter = Det(mat_dfj);
	if (deter != Real_wp(0))
	  {
	    GetInverse(mat_dfj, dfjm1);
	    MltTrans(dfjm1, this->NormaleLoc(ne), normale_fj);
	    Mlt(deter, normale_fj);
	    dsj = Norm2(normale_fj);
	    Mlt(Real_wp(Real_wp(1)/dsj), normale_fj);
	  }
	else
	  {
	    dsj = 0.0;
	    normale_fj.Zero();
	  }
	
	res.SetNormaleQuadratureBoundary(k, normale_fj);
	res.SetDsQuadratureBoundary(k, dsj);
      }
    
    int num_edge = mesh.Element(nquad).numBoundary(ne);
    this->elt_geom.GetCurvatureOnSurface(s, PTReel, res, mesh, nquad, num_edge, ne);
  }


  //! computes res = Fi(point) for points located on the edges
  /*!
    \param[in] s list of vertices of the triangle
    \param res references points after transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad number of the triangle in the mesh
    \param[in] ne local edge number
  */
  template<class Dimension>
  void ElementReference_Dim<Dimension>::
  FjSurfaceElemDof(const VectR_N& s, SetPoints<Dimension>& res,
		   const Mesh<Dimension>& mesh, int nquad, int ne) const
  {
    int nb_points = this->GetNbPointsDofSurface(ne);
    res.ReallocatePointsDofBoundary(nb_points);
    for (int m = 0; m < nb_points; m++)
      {
	int num_point = this->num_dof_points_surf(ne)(m);
	res.SetPointDofBoundary(m, res.GetPointDof(num_point));
      }
  }
  
  
  //! computes res = DFi(point) for points located on the edges
  /*!
    \param[in] s list of vertices of the triangle
    \param[in] PTReel references points after transformation Fi
    \param res jacobian matrices
    \param[in] mesh given mesh
    \param[in] nquad triangle number in the mesh
    \param[in] ne local edge number
    this method computes the normal and surfacic element (ds) on each quadrature point
  */
  template<class Dimension>
  void ElementReference_Dim<Dimension>::
  DFjSurfaceElemDof(const VectR_N& s, const SetPoints<Dimension>& PTReel,
		    SetMatrices<Dimension>& res,
		    const Mesh<Dimension>& mesh, int nquad, int ne) const
  {
    int nb_points = this->GetNbPointsDofSurface(ne);
    res.ReallocatePointsDofBoundary(nb_points);
    for (int k = 0; k < nb_points; k++)
      {
	int num_point = this->num_dof_points_surf(ne)(k);
	res.SetPointDofBoundary(k, res.GetPointDof(num_point));
      }
  }
  

  //! computes the local prolongation operator
  template<class Dimension>
  void ElementReference_Dim<Dimension>::
  ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
                           const ElementReference_Dim<Dimension>& FaceCoarse,
                           const ElementReference_Dim<Dimension>& FaceFine) const
  {
    cout << "Not implemented in the general case" << endl;
    abort();
  }


  //! computes the local prolongation operator (with first-order basis functions)
  template<class Dimension>
  void ElementReference_Dim<Dimension>::
  ComputeLocalProlongationLowOrder(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
                                   const ElementReference_Dim<Dimension>& FaceCoarse) const
  {
    cout << "not implemented" << endl;
    abort();
  }
  

  //! Computes projection from ptsA to ptsB
  template<>
  void ElementReference_Dim<Dimension2>::
  ComputeProjectionPointsSurf(int num_loc, const VectReal_wp& ptsA,
                              const VectReal_wp& ptsB, Matrix<Real_wp>& proj) const
  {
    Globatto<Real_wp> Fb;
    Fb.AffectPoints(ptsB);
    if (ptsA.GetM() < ptsB.GetM())
      {
	cout << "Case not implemented in FillVariableSource" << endl;
	cout << "Add more points" << endl;
	abort();
      }
    
    Matrix<Real_wp> VDM(ptsA.GetM(), ptsB.GetM());
    VectReal_wp phi, tau;
    for (int i = 0; i < ptsA.GetM(); i++)
      {
	Fb.ComputeValuesPhiRef(ptsA(i), phi);
	SetRow(phi, i, VDM); 
      }

    GetQR(VDM, tau);
    
    proj.Reallocate(ptsB.GetM(), ptsA.GetM());
    for (int j = 0; j < ptsA.GetM(); j++)
      {
	phi.Reallocate(ptsA.GetM());
	phi.Zero();
	phi(j) = Real_wp(1);
	SolveQR(VDM, tau, phi);
	for (int i = 0; i < ptsB.GetM(); i++)
	  proj(i, j) = phi(i);
      }
  }
  

  //! Computes projection from ptsA to ptsB  
  template<>
  void ElementReference_Dim<Dimension3>
  ::ComputeProjectionPointsSurf(int num_loc, const VectR2& ptsA,
				const VectR2& ptsB, Matrix<Real_wp>& proj) const
  {
    if (this->elt_geom.IsLocalFaceQuadrilateral(num_loc))
      {
	VectR2 ptsB_s(ptsB);
	Vector<int> permut(ptsB.GetM());
	permut.Fill();
	Sort(ptsB_s, permut);

	int ny = ptsB_s.GetM();
	for (int i = 0; i < ptsB_s.GetM(); i++)
	  if (abs(ptsB_s(i)(0) - ptsB_s(0)(0)) > 1e-12)
	    {
	      ny = i;
	      break;
	    }

	int nx = ptsB_s.GetM() / ny;
	if (nx*ny != ptsB_s.GetM())
	  {
	    cout << "Non-tensorized quadrature points not treated" << endl;
	    abort();
	  }
	
	// tensorized points
	VectReal_wp x_div(nx), y_div(ny);
	for (int i = 0; i < nx; i++)
	  x_div(i) = ptsB_s(i*ny)(0);

	for (int i = 0; i < ny; i++)
	  y_div(i) = ptsB_s(i)(1);

	for (int i = 0; i < nx; i++)
	  for (int j = 0; j < ny; j++)
	    if (ptsB_s(j + i*ny) != R2(x_div(i), y_div(j)))	      
	      {
		cout << "Non-tensorized quadrature points not treated" << endl;
		abort();
	      }
	
	Globatto<Real_wp> Fb_x, Fb_y;
	Fb_x.AffectPoints(x_div);
	Fb_y.AffectPoints(y_div);
	
	Matrix<Real_wp> VDM(ptsA.GetM(), ptsB.GetM());
	VectReal_wp phi, phi_x, phi_y, tau;
	for (int i = 0; i < ptsA.GetM(); i++)
	  {
	    Fb_x.ComputeValuesPhiRef(ptsA(i)(0), phi_x);
	    Fb_y.ComputeValuesPhiRef(ptsA(i)(1), phi_y);
	    for (int j1 = 0; j1 < nx; j1++)
	      for (int j2 = 0; j2 < ny; j2++)
		{
		  int n = j2 + j1*ny;
		  VDM(i, permut(n)) = phi_x(j1)*phi_y(j2);
		}
	  }
	
	GetQR(VDM, tau);
    
	proj.Reallocate(ptsB.GetM(), ptsA.GetM());
	for (int j = 0; j < ptsA.GetM(); j++)
	  {
	    phi.Reallocate(ptsA.GetM());
	    phi.Zero();
	    phi(j) = Real_wp(1);
	    SolveQR(VDM, tau, phi);
	    for (int i = 0; i < ptsB.GetM(); i++)
	      proj(i, j) = phi(i);
	  }
      }
    else
      {
	cout << "not implemented" << endl;
	abort();
      }
  }

  
  //! Fills pos with integers such that pos(i) & power_two(num_loc) is
  //! true if the dof i belongs to boundary num_loc
  void FillPositionDofBoundaries(const Matrix<int>& FacesDof,
				 Vector<int>& power_two, Vector<int>& pos)
  {
    int nb_faces = FacesDof.GetN();
    int nb_dof = 0;
    for (int i = 0; i < FacesDof.GetM(); i++)
      for (int j = 0; j < FacesDof.GetN(); j++)
	nb_dof = max(nb_dof, FacesDof(i, j));
    
    power_two.Reallocate(nb_faces);
    pos.Reallocate(nb_dof+1);
    pos.Fill(0);
    power_two(0) = 2;
    for (int num_loc = 0; num_loc < nb_faces; num_loc++)
      {	
	if (num_loc > 0)
	  power_two(num_loc) = power_two(num_loc-1)*2;
	
	for (int j = 0; j < FacesDof.GetM(); j++)
	  if (FacesDof(j, num_loc) >= 0)
	    {
	      int i = FacesDof(j, num_loc);
	      pos(i) = pos(i) | power_two(num_loc);
	    }
      }
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_ELEMENT_REFERENCE_CXX
#endif
