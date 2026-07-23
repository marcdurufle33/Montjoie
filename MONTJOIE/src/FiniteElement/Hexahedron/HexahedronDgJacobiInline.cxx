#ifndef MONTJOIE_FILE_HEXAHEDRON_DG_JACOBI_INLINE_CXX

namespace Montjoie
{

  inline Real_wp HexahedronDgJacobi::GetValuePhi(int num_dof, int num_point) const
  {
    return Value_Phi(num_dof, num_point);
  }
  
  
  //! returns \f$ \varphi^{2D}_k(\xi_j) \f$
  inline Real_wp HexahedronDgJacobi::GetValuePhiOnBoundary(int k, int j, int num_loc) const
  {
    return Value_Phi(node, this->num_quad_points_surf(num_loc)(num_point));
  }

  
  //! computes \f$ \nabla \varphi_{numdof} (\xi_{numpoint}) \f$
  inline const R3& HexahedronDgJacobi::
  GetGradientPhi(int num_dof, int num_point) const
  {
    return Gradient_Phi(num_dof, num_point);
  }
    
  
  //! returns \f$ \varphi_k(\xi_j) \f$
  /*!
    where \f$ \varphi \mbox{ and } \xi \f$
    are restriction of functions on boundary and quadrature point on the boundary
  */
  inline const R3& HexahedronDgJacobi::
  GetGradientPhiOnBoundary(int node, int num_point, int num_loc) const
  {
    return GetGradientPhi(node, this->num_quad_points_surf(num_loc)(num_point));
  }
  
  
  //! computation of u on quadrature points of a face
  /*!
    \param[in] pts Points after transformation F_i
    \param[in] mat jacobian matrices DF_i
    \param[in] feval dof components of u
    \param[out] res values of u on quadrature points of the face
    \param[in] vars given problem
    \param[in] FaceBasis finite element
    \param[in] n element number
    \param[in] num_loc local position of the face in the element
   */
  template<class TypeElement, class Vector1, class Vector2>
  inline void HexahedronDgJacobi::
  ComputeValueBoundary(const SetPoints<Dimension3>& pts,
		       const SetMatrices<Dimension3>& mat,
		       const Vector1& feval, Vector2& res,
		       TypeElement& FaceBasis,
                       const Mesh<Dimension3>& mesh, int n, int num_loc) const
  {
    // to do
    abort();
  }
  
}

#define MONTJOIE_FILE_HEXAHEDRON_DG_JACOBI_INLINE_CXX  
#endif

