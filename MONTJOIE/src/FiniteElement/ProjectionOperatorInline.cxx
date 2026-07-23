#ifndef MONTJOIE_FILE_PROJECTION_OPERATOR_INLINE_CXX

namespace Montjoie
{

  //! Destructor
  inline FiniteElementProjector::~FiniteElementProjector()
  {
  }

  
  /******************
   * DenseProjector *
   ******************/

  
  template<class Dimension> template<class Dim>
  inline void DenseProjector<Dimension>::
  InitProjector(const ElementGeomReference<Dim>&, const typename Dim::VectR_N& Pts)
  {
    abort();
  }

  
  template<class Dimension> template<class Dim>
  inline void DenseProjector<Dimension>
  ::InitProjector(const ElementGeomReference<Dim>&, const VectReal_wp& subdiv,
                  const typename Dim::VectR_N& Pts)
  {
    abort();
  }

  
  //! implemented in derived classes
  template<class Dimension> 
  inline void DenseProjector<Dimension>::Init(const ElementGeomReference<Dimension1>& FaceBasis,
                                              const VectReal_wp& points_div)
  {
    InitProjector(FaceBasis, points_div);
  }
    

  //! implemented in derived classes
  template<class Dimension> 
  inline void DenseProjector<Dimension>::Init(const ElementGeomReference<Dimension1>& FaceBasis,
                                              const VectReal_wp& x,
                                              const VectReal_wp& points_div)
  {
    InitProjector(FaceBasis, x, points_div);
  }
  
  
  //! implemented in derived classes
  template<class Dimension> 
  inline void DenseProjector<Dimension>::Init(const ElementGeomReference<Dimension2>& FaceBasis,
                                              const VectR2& points_div) 
  {
    InitProjector(FaceBasis, points_div);
  }
    
  
  //! implemented in derived classes
  template<class Dimension> 
  inline void DenseProjector<Dimension>::Init(const ElementGeomReference<Dimension2>& FaceBasis,
                                              const VectReal_wp& x,
                                              const VectR2& points_div)
  {
    InitProjector(FaceBasis, x, points_div);
  }
  
  
  //! implemented in derived classes
  template<class Dimension> 
  inline void DenseProjector<Dimension>::Init(const ElementGeomReference<Dimension3>& FaceBasis,
                                              const VectR3& points_div)
  {
    InitProjector(FaceBasis, points_div);
  }
  
  
  //! implemented in derived classes
  template<class Dimension> 
  inline void DenseProjector<Dimension>::Init(const ElementGeomReference<Dimension3>& FaceBasis,
                                              const VectReal_wp& x,
                                              const VectR3& points_div)
  {
    InitProjector(FaceBasis, x, points_div);
  }


  /***********************
   * TensorizedProjector *
   ***********************/


  //! implemented in derived classes
  inline void TensorizedProjector<Dimension2>::Init(const ElementGeomReference<Dimension1>& FaceBasis,
                                                    const VectReal_wp& points_div)
  {
    abort();
  }
    

  //! implemented in derived classes
  inline void TensorizedProjector<Dimension2>::Init(const ElementGeomReference<Dimension1>&, const VectReal_wp&,
                                                    const VectReal_wp& points_div)
  {
    abort();
  }
  
  
  //! implemented in derived classes
  inline void TensorizedProjector<Dimension2>::Init(const ElementGeomReference<Dimension3>& FaceBasis,
                                                      const VectR3& points_div)
  {
    abort();
  }
  
  
  //! implemented in derived classes
  inline void TensorizedProjector<Dimension2>::Init(const ElementGeomReference<Dimension3>&, const VectReal_wp&,
                                                    const VectR3& points_div)
  {
    abort();
  }  
  

#ifdef MONTJOIE_WITH_THREE_DIM
  //! implemented in derived classes
  inline void TensorizedProjector<Dimension3>::Init(const ElementGeomReference<Dimension1>& FaceBasis,
                                                    const VectReal_wp& points_div)
  {
    abort();
  }
    

  //! implemented in derived classes
  inline void TensorizedProjector<Dimension3>::Init(const ElementGeomReference<Dimension1>&, const VectReal_wp&,
                                                    const VectReal_wp& points_div)
  {
    abort();
  }
  
  
  //! implemented in derived classes
  inline void TensorizedProjector<Dimension3>::Init(const ElementGeomReference<Dimension2>& FaceBasis,
                                                    const VectR2& points_div) 
  {
    abort();
  }
    
  
  //! implemented in derived classes
  inline void TensorizedProjector<Dimension3>::Init(const ElementGeomReference<Dimension2>&, const VectReal_wp&,
                                                      const VectR2& points_div)
  {
    abort();
  }
#endif
  

  /*****************************
   * FiniteElementInterpolator *
   *****************************/
  
  
  //! interpolation from nodal points of the finite element to specified points
  //! during initialisation
  /*!
    \param[in] U values of u on nodal points of the finite element
    \param[in] V values of u on points_div provided in InitProjection
    \param[in] type_elt hybrid type of finite element (triangle, quad, tet, etc)
    By denoting \xi_i points provided by the user when calling InitProjection,
    and phi_i the nodal shape functions of the finite element, this function computes
    V_j = \sum_i phi_i(xi_j) U_i
   */
  template<class T1, class T2, int nb_unknowns>
  inline void FiniteElementInterpolator
  ::Project(const Vector<TinyVector<T1, nb_unknowns> >& U,
	    Vector<TinyVector<T2, nb_unknowns> >& V, int type_elt) const
  {
    stored_proj(type_elt)->Project(U, V);
  }
  
  
  //! transpose operation of interpolation from nodal points of the finite element 
  //! to specified points during initialisation
  /*!
    \param[in] U values of u on points_div provided in InitProjection
    \param[in] V projection of u on nodal points of the finite element
    \param[in] type_elt hybrid type of finite element (triangle, quad, tet, etc)
    By denoting \xi_i points provided by the user when calling InitProjection,
    and phi_i the nodal shape functions of the finite element, this function computes
    V_i = \sum_j phi_i(xi_j) U_j
    This operation is the transpose of method Project
   */
  template<class T1, class T2, int nb_unknowns>
  inline void FiniteElementInterpolator
  ::TransposeProject(const Vector<TinyVector<T1, nb_unknowns> >& U,
		     Vector<TinyVector<T2, nb_unknowns> >& V, int type_elt) const
  {
    stored_proj(type_elt)->Project(U, V);
  }


  //! interpolation from nodal points of the finite element to specified points
  //! during initialisation
  /*!
    \param[in] U values of u on nodal points of the finite element
    \param[in] V values of u on points_div provided in InitProjection
    \param[in] type_elt hybrid type of finite element (triangle, quad, tet, etc)
    By denoting \xi_i points provided by the user when calling InitProjection,
    and phi_i the nodal shape functions of the finite element, this function computes
    V_j = \sum_i phi_i(xi_j) U_i
   */
  template<class T1, class T2, int p, int q>
  inline void FiniteElementInterpolator
  ::Project(const Vector<TinyMatrix<T1, General, p, q> >& U,
	    Vector<TinyMatrix<T2, General, p, q> >& V, int type_elt) const
  {
    stored_proj(type_elt)->Project(U, V);
  }

  
  //! interpolation from nodal points of the finite element 
  //! to specified points during initialisation
  /*!
    \param[in] U values of u on nodal points of the finite element
    \param[in] V values of u on points_div provided in InitProjection
    \param[in] type_elt hybrid type of finite element (triangle, quad, tet, etc)
    By denoting \xi_i points provided by the user when calling InitProjection,
    and phi_i the nodal shape functions of the finite element, this function computes
    V_j = \sum_i phi_i(xi_j) U_i
   */
  template<class Vector1, class Vector2>
  inline void FiniteElementInterpolator
  ::ProjectScalar(const Vector1& U, Vector2& V, int type_elt) const
  {
    stored_proj(type_elt)->Project(U, V);
  }
  
  
  //! transpose operation of interpolation from nodal points of the finite element 
  //! to specified points during initialisation
  /*!
    \param[in] U values of u on points_div provided in InitProjection
    \param[in] V projection of u on nodal points of the finite element
    \param[in] type_elt hybrid type of finite element (triangle, quad, tet, etc)
    By denoting \xi_i points provided by the user when calling InitProjection,
    and phi_i the nodal shape functions of the finite element, this function computes
    V_i = \sum_j phi_i(xi_j) U_j
    This operation is the transpose of method ProjectScalar
   */
  template<class Vector1, class Vector2>
  inline void FiniteElementInterpolator
  ::TransposeProjectScalar(const Vector1& U, Vector2& V, int type_elt) const
  {
    stored_proj(type_elt)->TransposeProject(U, V);
  }

  
  template<class T1, class T2>
  inline void FiniteElementInterpolator::Project(const Vector<T1>& U, Vector<T2>& V, int type_elt) const
  {
    ProjectScalar(U, V, type_elt);
  }
  
  
  template<class T1, class T2>
  inline void FiniteElementInterpolator
  ::TransposeProject(const Vector<T1>& U, Vector<T2>& V, int type_elt) const
  {
    TransposeProjectScalar(U, V, type_elt);
  }
  
}
  
#define MONTJOIE_FILE_PROJECTION_OPERATOR_INLINE_CXX
#endif
