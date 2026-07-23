#ifndef MONTJOIE_FILE_ELEMENT_GEOM_REFERENCE_CXX

namespace Montjoie
{

  //! computation of coefficients useful to interpolate
  //! from the mesh order to the finite element order
  template<class Dimension> template<class DimensionB>
  void ElementGeomReference_Base<Dimension>::SetMesh(Mesh<DimensionB>& mesh)
  {
    // we consider that the order of approximation of curves in the mesh
    // can be different from order of finite element
    int r = mesh.GetGeometryOrder();
    
    // so we construct interpolation
    CoefCurve1D.Reallocate(order_geom-1, r+1);
    for (int i = 1; i < order_geom; i++)
      for (int j = 0; j <= r; j++)
	CoefCurve1D(i-1, j) = mesh.GetInterpolate1D(j, points_nodal1d(i));
  }

    
  //! returns size of memory used by the object
  template<class Dimension>
  size_t ElementGeomReference_Base<Dimension>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille += sizeof(Real_wp)*(points_nodal1d.GetM()+weights1d.GetM()+points1d.GetM());
    taille += sizeof(Real_wp)*(points_dof1d.GetM()+CoefCurve1D.GetDataSize());
    taille += sizeof(Real_wp)*Dimension::dim_N*(normale.GetM()+GradientPhi_Nodal.GetDataSize());
    return taille;
  }
  
  
  //! computes normale on a point from the jacobian matrix \f$ DF_i^{-1} \f$
  //! boundary element of integration is computed as well
  /*!
    \param[in] dfjm1 the jacobian matrix \f$ DF_i^{-1} \f$
    \param[in] num_edge local number of the edge/face in the element
    \param[out] normale_fj unit normale (outward to the element)
    \param[out] dsj surface element of integration
   */
  template<class Dimension>
  void ElementGeomReference_Base<Dimension>::
  GetNormale(const MatrixN_N& dfjm1, R_N& normale_fj,
             Real_wp& dsj, int num_edge) const
  {
    MltTrans(dfjm1, normale(num_edge), normale_fj); 
    Real_wp length = Norm2(normale_fj); 
    Real_wp det = Det(dfjm1);
    if (det != Real_wp(0))
      dsj = length/det;
    else
      dsj = 0;
    
    if (length != Real_wp(0))
      length = Real_wp(1)/length;
    
    Mlt(length, normale_fj);
  }
  

  //! transformation Fi in the case of curved element
  template<class Dimension>
  void ElementGeomReference_Base<Dimension>
  ::FjCurvePhi(const SetPoints<Dimension>& PTReel,
	       const R_N& point, R_N& res, VectReal_wp& phi) const
  {
    // values of shape functions are stored in phi
    // such that it can be reused by DFjCurvePhi
    res.Zero();
    ComputeValuesPhiNodalRef(point, phi);
    for (int node = 0; node < phi.GetM(); node++)
      {
	Real_wp coef = phi(node);
	Add(coef, PTReel.GetPointNodal(node), res);
      }
  }


  //! transformation DFi in the case of curved element
  template<class Dimension>
  void ElementGeomReference_Base<Dimension>
  ::DFjCurvePhi(const SetMatrices<Dimension>& Mat, const R_N& point, MatrixN_N& res, 
		const VectReal_wp& phi) const
  {
    res.Zero();
    for (int node = 0; node < phi.GetM(); node++)
      {
	Real_wp coef = phi(node);
	Add(coef, Mat.GetPointNodal(node), res);
      }
  }

    
  //! computation of gradient of u on nodal points
  /*!
    \param[in] Un values of u on nodal points
    \param[out] Unode gradient of u on nodal points
   */
  template<class Dimension>
  void ElementGeomReference_Base<Dimension>
  ::ComputeNodalGradientRef(const Vector<Real_wp>& Un, VectR_N& Unode) const
  {
    for (int i = 0; i < GradientPhi_Nodal.GetN(); i++)
      {
	Unode(i).Zero();
	for (int j = 0; j < GradientPhi_Nodal.GetM(); j++)
	  Add(Un(j), GradientPhi_Nodal(j, i), Unode(i));
      }
  }
  

  //! computes gradient of u on nodal points
  template<class Dimension>
  void ElementGeomReference_Base<Dimension>
  ::ComputeNodalGradientRef(const Vector<Complex_wp>& Un,
			    Vector<TinyVector<Complex_wp, Dimension::dim_N> >& Unode) const
  {
    for (int i = 0; i < GradientPhi_Nodal.GetN(); i++)
      {
	Unode(i).Zero();
	for (int j = 0; j < GradientPhi_Nodal.GetM(); j++)
	  Add(Un(j), GradientPhi_Nodal(j, i), Unode(i));
      }
  }


  //! computation of gradient of U on nodal points from nodal values of U
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] Uloc nodal values of U
    \param[out] grad_Uloc values of grad(U) on nodal points
  */
  template<class Dimension> template<class T1>
  void ElementGeomReference_Base<Dimension>::
  ComputeNodalGradient(const SetMatrices<Dimension>& MatricesElem,
                       const Vector<T1>& Uloc, Vector<TinyVector<T1, Dimension::dim_N> >& grad_Uloc) const
  {
    int nb_points_nodal = static_cast<const ElementGeomReference<Dimension>& >(*this).GetNbPointsNodalElt();
    grad_Uloc.Reallocate(nb_points_nodal);
    TinyVector<T1, Dimension::dim_N> grad;
    Vector<TinyVector<T1, Dimension::dim_N> > Unode(nb_points_nodal);
    
    // projection on nodal points
    this->ComputeNodalGradientRef(Uloc, Unode);

    // loop on each nodal point
    MatrixN_N dfjm1;
    for (int i = 0; i < nb_points_nodal; i++)
      {
	// we get the transformation DF_i^{*-1}
	GetInverse(MatricesElem.GetPointNodal(i), dfjm1);
	// applying DF_i^{*-1}
	MltTrans(dfjm1, Unode(i), grad);
	// storing value
	grad_Uloc(i) = grad;
      }
  }


  //! computes the gradient of different components of U, result is given in gradU
  template<class Dimension> template<class T>
  void ElementGeomReference_Base<Dimension>
  ::ComputeNodalGradient(const SetMatrices<Dimension>& MatricesElem,
			 const Vector<Vector<T> >& U, Vector<Vector<T> >& gradU) const
  {
    int nb_points_nodal = static_cast<const ElementGeomReference<Dimension>& >(*this).GetNbPointsNodalElt();
    Vector<TinyVector<T, Dimension::dim_N> > Unode(nb_points_nodal);
    TinyVector<T, Dimension::dim_N> grad, grad_chap;
    int nb_unknowns = U.GetM();
    if (gradU.GetM() < nb_unknowns*Dimension::dim_N)
      gradU.Reallocate(nb_unknowns*Dimension::dim_N);

    for (int k = 0; k < gradU.GetM(); k++)
      gradU(k).Reallocate(nb_points_nodal);

    int offset = 0;
    for (int i = 0; i < nb_unknowns; i++)
      {
	this->ComputeNodalGradientRef(U(i), Unode);
	for (int j = 0; j < nb_points_nodal; j++)
	  for (int k = 0; k < Dimension::dim_N; k++)
	    gradU(offset+k)(j) = Unode(j)(k);
	
	offset += Dimension::dim_N;
      }
    
    // loop on each nodal point
    MatrixN_N dfjm1;
    for (int j = 0; j < nb_points_nodal; j++)
      {
	// we get the transformation DF_i^{*-1}
	GetInverse(MatricesElem.GetPointNodal(j), dfjm1);
	
	// applying DF_i^{*-1}
	offset = 0;
	for (int i = 0; i < nb_unknowns; i++)
	  {
	    for (int k = 0; k < Dimension::dim_N; k++)
	      grad_chap(k) = gradU(offset+k)(j);
	    
	    MltTrans(dfjm1, grad_chap, grad);
	    
	    for (int k = 0; k < Dimension::dim_N; k++)
	      gradU(offset+k)(j) = grad(k);
	    
	    offset += Dimension::dim_N;
	  }
      }
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_ELEMENT_GEOM_REFERENCE_CXX
#endif
