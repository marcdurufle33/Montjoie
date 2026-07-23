#ifndef MONTJOIE_FILE_ELEMENTARY_MATRIX_H1_HXX

namespace Montjoie
{

  template<class Complexe, class TypeEquation,
	   class Vector1, class Vector2, class Vector3>
  void ComputeVariableElementaryMatrix(int iquad, VirtualMatrix<Complexe> & mat,
				       Vector1& Amass,
				       Vector2& Cgrad_grad, Vector3& Dgrad_phi, Vector3& Ephi_grad,
				       const GlobalGenericMatrix<Complexe>& nat_mat,
				       const EllipticProblem<TypeEquation>& vars,
				       const ElementReference<typename TypeEquation::Dimension, 1>& Fb);

  template<class Complexe, class TypeEquation>
  void ComputeElementaryMatrix(int iquad, IVect& num_dof,
			       VirtualMatrix<Complexe>& mat_interac,
			       const GlobalGenericMatrix<Complexe>& nat_mat,
			       const EllipticProblem<TypeEquation>& vars,
			       const ElementReference<typename TypeEquation::Dimension, 1>& Fb);
  
  template<class Complexe, class TypeEquation>
  void AddElementaryFluxDG(VirtualMatrix<Complexe>& mat_sp, const GlobalGenericMatrix<Complexe>& nat_mat,
			   const EllipticProblem<TypeEquation>& vars,
			   int iquad, int num_pos1_face, int num_elem2,
			   int num_face, int ref, bool new_face, int& nb_neighbor,
                           int offset_row, int offset_col,
			   const ElementReference<typename TypeEquation::Dimension, 1>& Fb);
  
  template<class Complexe, class TypeEquation>
  void AddElementaryFluxesDG(VirtualMatrix<Complexe>& mat_sp,
			     const GlobalGenericMatrix<Complexe>& nat_mat,
			     const EllipticProblem<TypeEquation>& vars,
                             int offset_row, int offset_col);
  
} // namespace Montjoie

#define MONTJOIE_FILE_ELEMENTARY_MATRIX_H1_HXX
#endif
