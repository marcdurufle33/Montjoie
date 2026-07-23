#ifndef MONTJOIE_FILE_NUMBER_MAP_INLINE_CXX

namespace Montjoie
{
  
  //! default constructor
  inline NumberMap::NumberMap()
  {
    only_hexahedral = false;
    dg_element = 0;
    nb_dof_vertex = 0;
  }
  
  
  //! returns the number of dofs per vertex
  inline int NumberMap::GetNbDofVertex(int order) const
  {
    return nb_dof_vertex;
  }
  
  
  //! returns the number of dofs per edge for a given order
  inline int NumberMap::GetNbDofEdge(int order) const
  {
    if (nb_dof_edge.GetM() > order)
      return nb_dof_edge(order);
    
    return 0;
  }
  
  
  //! returns the number of dofs per triangle for a given order
  inline int NumberMap::GetNbDofTriangle(int order) const
  {
    if (nb_dof_tri.GetM() > order)
      return nb_dof_tri(order);
    
    return 0;
  }
  
  
  //! returns the number of dofs per quadrangle for a given order
  inline int NumberMap::GetNbDofQuadrangle(int order) const
  {
    if (nb_dof_quad.GetM() > order)
      return nb_dof_quad(order);
    
    return 0;
  }
  
  
  //! returns the number of dofs per tetrahedron for a given order
  inline int NumberMap::GetNbDofTetrahedron(int order) const
  {
    if (nb_dof_tetra.GetM() > order)
      return nb_dof_tetra(order);
    
    return 0;
  }
  
  
  //! returns the number of dofs per pyramid for a given order
  inline int NumberMap::GetNbDofPyramid(int order) const
  {
    if (nb_dof_pyramid.GetM() > order)
      return nb_dof_pyramid(order);
    
    return 0;
  }
  
  
  //! returns the number of dofs per prism for a given order
  inline int NumberMap::GetNbDofWedge(int order) const
  {
    if (nb_dof_wedge.GetM() > order)
      return nb_dof_wedge(order);
    
    return 0;
  }
  
  
  //! returns the number of dofs per hexahedron for a given order
  inline int NumberMap::GetNbDofHexahedron(int order) const
  {
    if (nb_dof_hexa.GetM() > order)
      return nb_dof_hexa(order);
    
    return 0;
  }
  

  //! returns the dof number after a symmetry of the edge
  inline int NumberMap::GetSymmetricEdgeDof(int order, int k) const
  {
    if (EdgesDof_Symmetry.GetM() > order)
      return EdgesDof_Symmetry(order)(k);
    
    return 0;
  }
  
  
  //! returns true if the dof is skew-symmetric 
  //! (i.e. the sign depends on the orientation of the edge)
  inline bool NumberMap::IsSkewSymmetricEdgeDof(int r, int k) const
  {
    return EdgesDof_SkewSymmetry(r)(k);
  }
  

  //! returns the number of dof points on an edge
  inline int NumberMap::GetNbDofBoundary(int order, const Edge<Dimension2>& elt) const
  {
    return nb_dof_edge(order);
  }
  
  
  //! returns the number of dof points on a face
  inline int NumberMap::GetNbDofBoundary(int order, const Face<Dimension3>& elt) const
  {
    int nv = elt.GetNbVertices();
    if (nv == 4)
      return nb_dof_quad(order);
    
    return nb_dof_tri(order);
  }


  //! returns the number of dofs strictly inside the element
  template<class Dimension>
  inline int NumberMap::GetNbDofElement(int order, const Edge<Dimension>& elt) const
  {
    return nb_dof_edge(order);
  }

  
  //! returns the number of quadrature points on an edge
  inline int NumberMap::GetNbPointsQuadBoundary(int order, const Edge<Dimension2>& elt) const
  {
    return points_quadrature_edge(order).GetM();
  }
  
  
  //! returns the number of quadrature points on a face
  inline int NumberMap::GetNbPointsQuadBoundary(int order, const Face<Dimension3>& elt) const
  {
    int nv = elt.GetNbVertices();
    if (nv == 4)
      return points_quadrature_quad(order).GetM();
    
    return points_quadrature_tri(order).GetM();
  }

  
  //! returns local position of a quadrature point of an edge
  inline Real_wp NumberMap::GetQuadraturePoint(int rf, int k, const Edge<Dimension2>&) const
  {
    return points_quadrature_edge(rf)(k);
  }
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! returns local position of a quadrature point of a face
  inline const R2& NumberMap
  ::GetQuadraturePoint(int rf, int k, const Face<Dimension3>& face) const
  {
    if (face.GetNbVertices() == 3)
      return points_quadrature_tri(rf)(k);
    
    return points_quadrature_quad(rf)(k);
  }
#endif
  
  
  //! returns quadrature points used for the unit interval
  inline const VectReal_wp& NumberMap::GetEdgeQuadrature(int order) const
  {
    return points_quadrature_edge(order);
  }
  
  
  //! returns quadrature points for the unit triangle
  inline const VectR2& NumberMap::GetTriangleQuadrature(int order) const
  {
    return points_quadrature_tri(order);
  }
  
  
  //! returns quadrature points for the unit square
  inline const VectR2& NumberMap::GetQuadrangleQuadrature(int order) const
  {
    return points_quadrature_quad(order);
  }

  
  //! sets the number of dofs per vertex
  inline void NumberMap::SetNbDofVertex(int order, int nb_dof)
  {
    nb_dof_vertex = nb_dof;
  }
  
  
  //! returns the half of quadrature weights to use for an edge
  inline const VectReal_wp& NumberMap
  ::GetFluxWeight(int order, const Edge<Dimension2>& elt) const
  {
    return half_weight_edge(order);
  }
  
  
  //! returns the half of quadrature weights to use for a face
  inline const VectReal_wp& NumberMap
  ::GetFluxWeight(int order, const Face<Dimension3>& elt) const
  {
    if (elt.GetNbVertices() == 3)
      return half_weight_tri(order);
    
    return half_weight_quad(order);
  }
  
  
  //! returns true if these numbering rules are intended to Discontinuous Galerkin
  inline int NumberMap::FormulationDG() const
  {
    return dg_element;
  }
  
  
  //! if dg_elt is true, these numbering rules are intended to Discontinuous Galerkin
  inline void NumberMap::SetFormulationDG(int dg_elt)
  {
    dg_element = dg_elt;
  }
  
  
  //! sets the new dof number (and its sign) after symmetry of dof i
  inline void NumberMap::SetEdgesDofSymmetry(int r, int i, int j, bool skew_sym)
  {
    EdgesDof_Symmetry(r)(i) = j;
    EdgesDof_SkewSymmetry(r)(i) = skew_sym;
  }
  
  
  //! All dofs are skew-symmetric (signs are changed after symmetry)
  inline void NumberMap::SetAllEdgesDofToSkewSymmetric(int r)
  {
    EdgesDof_SkewSymmetry(r).Fill(true); 
  }
  
  
  //! Odd dofs are skew-symmetric
  inline void NumberMap::SetOddEdgesDofToSkewSymmetric(int r)
  {
    for (int i = 0; i < EdgesDof_SkewSymmetry(r).GetM(); i++)
      EdgesDof_SkewSymmetry(r)(i) = (i%2 == 1); 
  }

  
  //! Even dofs are skew-symmetric
  inline void NumberMap::SetEvenEdgesDofToSkewSymmetric(int r)
  {
    for (int i = 0; i < EdgesDof_SkewSymmetry(r).GetM(); i++)
      EdgesDof_SkewSymmetry(r)(i) = (i%2 == 0); 
  } 
  
  
  //! Provides direclty signs of dofs after symmetry
  inline void NumberMap::SetSkewSymmetryEdgesDof(int r, const Vector<bool>& skew)
  {
    EdgesDof_SkewSymmetry(r) = skew;
  }
  
  
  //! returns true if the dofs of the face are invariant by rotation
  inline bool NumberMap::DofInvariantByRotation(int order, const Face<Dimension3>& f) const
  {
    if (f.GetNbVertices() == 3)
      {
        if (general_combination_tri.GetM() <= 0)
          return true;

        return !general_combination_tri(order);
      }
    else
      {
        if (general_combination_quad.GetM() <= 0)
          return true;

        return !general_combination_quad(order);
      }
  }
  
  
  //! returns linear operator to apply to dofs of a triangular face after rotation
  inline const Vector<Matrix<Real_wp, General, RowSparse> >& NumberMap
  ::GetFacesDofRotationTri(int order) const
  {
    return coef_combination_tri(order);
  }
  
  
  //! provides dof changes of sign after rotation of a triangular face  
  inline void NumberMap::SetSignDofRotationTri(int order, const Matrix<bool>& is_neg)
  {
    SignDof_Rotation_Tri(order) = is_neg;
  }


  //! provides dof changes of sign after rotation of a quadrilateral face  
  inline void NumberMap::SetSignDofRotationQuad(int order, const Matrix<bool>& is_neg)
  {
    SignDof_Rotation_Quad(order) = is_neg;
  }


  //! returns the dof number after rotation of the face
  inline int NumberMap::GetRotationFaceDof(int rot, int r, int k, const Face<Dimension3>& f) const
  {
    if (f.GetNbVertices() == 3)
      return FacesDof_Rotation_Tri(r)(rot, k);
    
    return FacesDof_Rotation_Quad(r)(rot, k);
  }
  

  //! returns the dof number after rotation of the edge
  inline int NumberMap::GetRotationFaceDof(int rot, int r, int k, const Edge<Dimension2>& f) const
  {
    if (rot == 0)
      return k;
    
    return EdgesDof_Symmetry(r)(k);
  }
  
  
  //! returns quadrature points numbers after rotation of the face
  inline const Matrix<int>& NumberMap::
  GetRotationQuadraturePoints(int r, const Face<Dimension3>& f) const
  {
    if (f.GetNbVertices() == 3)
      return FacesQuad_Rotation_Tri(r);
    
    return FacesQuad_Rotation_Quad(r);
  }


  //! returns quadrature points numbers after rotation of the face
  inline const Matrix<int>& NumberMap::
  GetRotationQuadraturePoints(int r, int nb_vertices) const
  {
    if (nb_vertices == 3)
      return FacesQuad_Rotation_Tri(r);
    
    return FacesQuad_Rotation_Quad(r);
  }
  

  //! returns quadrature points numbers after rotation of the edge
  inline const Matrix<int>& NumberMap::
  GetRotationQuadraturePoints(int r, const Edge<Dimension2>& f) const
  {
    return EdgesQuad_Rotation(r);
  }
  
  
  //! returns true if the sign is the same after rotation of the face
  inline bool NumberMap::
  IsNegativeDof(int rot, int r, int k, const Face<Dimension3>& f) const
  {
    if (f.GetNbVertices() == 3)
      return SignDof_Rotation_Tri(r)(rot, k);
    
    return SignDof_Rotation_Quad(r)(rot, k);
  }
  
  
  //! projection from local dofs to global dofs by applying the linear operator of the face
  template<class Vector1>
  inline void NumberMap
  ::ProjectToGlobalDofs(Vector1& U, int rot, int r,
			int offset, const Face<Dimension3>& f) const
  {
    if (f.GetNbVertices() == 3)
      ApplyOperator(coef_combination_tri(r)(rot), U, offset, nb_dof_tri(r));
    else
      ApplyOperator(coef_combination_quad(r)(rot), U, offset, nb_dof_quad(r));
  }
  

  //! transpose projection from local dofs to global dofs
  //! by applying the linear operator of the face
  template<class Vector1>
  inline void NumberMap
  ::TransposeProjectToGlobalDofs(Vector1& U, int rot, int r,
				 int offset, const Face<Dimension3>& f) const
  {
    if (f.GetNbVertices() == 3)
      ApplyTransposeOperator(coef_combination_tri(r)(rot), U, offset, nb_dof_tri(r));
    else
      ApplyTransposeOperator(coef_combination_quad(r)(rot), U, offset, nb_dof_quad(r));
  }
  
  
  //! projection from global to local dofs bt applying the linear operator of the face
  template<class Vector1>
  inline void NumberMap::ProjectToLocalDofs(Vector1& U, int rot, int r,
					    int offset, const Face<Dimension3>& f) const
  {
    if (f.GetNbVertices() == 3)
      ApplyTransposeOperator(invCoef_combination_tri(r)(rot), U, offset, nb_dof_tri(r));
    else
      ApplyTransposeOperator(invCoef_combination_quad(r)(rot), U, offset, nb_dof_quad(r));
  }
  
  
  //! modification of the rows of the elementary matrix because of the linear operator
  //! of the considered face
  template<class Matrix1>
  inline void NumberMap::ModifyRowToGlobalDofs(Matrix1& A, int rot, int r, int offset,
					       const Face<Dimension3>& f) const
  {
    if (f.GetNbVertices() == 3)
      ApplyOperatorRow(coef_combination_tri(r)(rot), A, offset, nb_dof_tri(r));
    else
      ApplyOperatorRow(coef_combination_quad(r)(rot), A, offset, nb_dof_quad(r));
  }


  //! modification of the columns of the elementary matrix because of the linear operator
  //! of the considered face
  template<class Matrix1>
  inline void NumberMap::ModifyColumnToGlobalDofs(Matrix1& A, int rot, int r,
						  int offset, const Face<Dimension3>& f) const
  {
    if (f.GetNbVertices() == 3)
      ApplyOperatorColumn(coef_combination_tri(r)(rot), A, offset, nb_dof_tri(r));
    else
      ApplyOperatorColumn(coef_combination_quad(r)(rot), A, offset, nb_dof_quad(r));
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_NUMBER_MAP_INLINE_CXX
#endif
