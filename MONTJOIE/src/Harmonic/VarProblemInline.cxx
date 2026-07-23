#ifndef MONTJOIE_FILE_VAR_PROBLEM_INLINE_CXX

namespace Montjoie
{
    
  /**************
   * VarProblem *
   **************/


  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline VarProblem<Dimension>::VarProblem(EllipticProblem<TypeEquation>& var)
    : VarComputationProblem(var), DistributedProblem<Dimension>(var), VarFiniteElement<Dimension>(),
      var_boundary(var)
  {
  }

  //! computes arrays needed for static condensation 
  template<class Dimension>
  inline void VarProblem<Dimension>
  ::GetStaticCondensedRows(IVect& IndexCondensedRows, IVect& global_row, IVect& overlap_row,
			   IVect& overlap_proc, int& nb_scalar_dof, int& nb_global_dof,
			   IVect& sharing_proc, Vector<IVect>& sharing_rows) const
  {
    DistributedProblem<Dimension>::
      GetStaticCondensedRows(IndexCondensedRows, global_row, overlap_row, overlap_proc,
			     nb_scalar_dof, nb_global_dof, sharing_proc, sharing_rows);
  }
  
  
  //! computes velocity on elements of the mesh
  template<class Dimension>
  inline void VarProblem<Dimension>::GetVelocityOnElements(VectReal_wp& velocity)
  {
    VarProblem_Base::GetVelocityOnElements(velocity, this->mesh);
  }
  

  //! returns coefficient involved in stiffness matrix
  template<class Dimension>
  inline Real_wp VarProblem<Dimension>::GetCoefficientPenaltyStiffness(int ref) const
  {
    return Real_wp(1);
  }

  //! copies finite elements stored in another problem
  template<class Dimension>
  inline void VarProblem<Dimension>::CopyFiniteElement(const DistributedProblem<Dimension>& var)
  {
    const VarProblem<Dimension>& var_p = static_cast<const VarProblem<Dimension>& >(var);
    VarFiniteElement<Dimension>
      ::CopyFiniteElement(var_p, this->mesh, this->mesh_num,
                          this->other_mesh_num, this->type_element,
                          this->other_type_element);
  }
  
  //! the mesh is split for parallel computation
  template<class Dimension>
  inline void VarProblem<Dimension>
  ::SplitMeshForParallelComputation(const string& name_elt)
  {
    DistributedProblem<Dimension>::SplitMeshForParallelComputation(name_elt);
  }
  
  
  //! returns the number of quadrature points inside element i
  template<class Dimension>
  inline int VarProblem<Dimension>::GetNbPointsQuadratureInside(int i) const
  {
    return this->GetReferenceElement(i).GetNbPointsQuadratureInside();
  }
  

  //! returns the number of dof points on all the boundaries of element i
  template<class Dimension>
  inline int VarProblem<Dimension>::GetNbDofBoundaries(int i, int nm) const
  {
    return this->GetReferenceElement(i, nm).GetNbDofBoundaries();
  }


  //! returns the number of dofs on element i
  template<class Dimension>
  inline int VarProblem<Dimension>::GetNbLocalDof(int i, int nm) const
  {
    return this->GetReferenceElement(i, nm).GetNbDof();
  }


  //! returns the number of dofs on surface i
  template<class Dimension>
  inline int VarProblem<Dimension>::GetNbSurfaceDof(int i, int nm) const
  {
    return this->GetSurfaceFiniteElement(i, nm).GetNbDof();
  }

  
  //! returns the number of elements
  template<class Dimension>
  inline int VarProblem<Dimension>::GetNbElt() const
  {
    return this->mesh.GetNbElt();
  }


  //! returns the print level
  template<class Dimension>
  inline int VarProblem<Dimension>::GetPrintLevel() const
  {
    return VarComputationProblem::GetPrintLevel();
  }
  
  
  //! returns weights of integration of element i
  template<class Dimension>
  inline const VectReal_wp& VarProblem<Dimension>::WeightsND(int i) const
  {
    return this->GetReferenceElement(i).WeightsND();
  }
  

  //! returns quadrature points on reference element (on face num_loc)
  template<class Dimension>
  inline const typename Dimension::DimensionBoundary::VectR_N& 
  VarProblem<Dimension>::PointsQuadratureBoundary(int i, int num_loc) const
  {
    return this->GetReferenceElement(i).PointsQuadratureBoundary(num_loc);
  }


  //! returns dofs points on reference element (on face num_loc)
  template<class Dimension>
  inline const typename Dimension::DimensionBoundary::VectR_N& 
  VarProblem<Dimension>::PointsDofBoundary(int i, int num_loc) const
  {
    return this->GetReferenceElement(i).PointsDofBoundary(num_loc);
  }
  
  
  //! returns reference to the geometric reference element of element i
  /*!
    The returned object contains the definition of "shape functions"
    that are used for transformation F_i
   */
  template<class Dimension>
  inline const ElementGeomReference<Dimension>& VarProblem<Dimension>::GetShapeElement(int i) const
  {
    return this->GetReferenceElement(i).GetGeometricElement();
  }


  //! returns a reference to the reference element associated with element i
  template<class Dimension>
  inline const ElementReference_Base& VarProblem<Dimension>::GetReferenceElementBase(int i, int n) const
  {
    return this->GetReferenceElement(i, n);
  }


  //! returns a reference to the reference element associated with boundary i
  template<class Dimension>
  inline const ElementReference_Base& VarProblem<Dimension>::GetSurfaceElementBase(int i, int n) const
  {
    return this->GetSurfaceFiniteElement(i, n);
  }


  //! returns true if the elementary mass matrix of element i is diagonal
  template<class Dimension>
  inline bool VarProblem<Dimension>
  ::DiagonalScalarMassMatrix(const ElementReference_Dim<Dimension>& Fb, int i) const
  {
    return Fb.DiagonalMassMatrix();
  }


  //! returns true if the elementary mass matrix of element i is block-diagonal
  template<class Dimension>
  inline bool VarProblem<Dimension>
  ::BlockDiagonalScalarMassMatrix(const ElementReference_Dim<Dimension>& Fb, int i) const
  {
    return Fb.LumpedMassMatrix();
  }
  
  
  //! returns quadrature points inside element i (on the reference element)
  template<class Dimension>
  inline typename Dimension::VectR_N VarProblem<Dimension>::PointsQuadInsideND(int i) const
  {
    return this->GetReferenceElement(i).PointsQuadInsideND();
  }

  
  //! computes E \times n and H \times n (for Maxwell's equations) on quadrature points of a surface
  template<class Dimension>
  inline void VarProblem<Dimension>
  ::ComputeEnHnOnBoundary(const MeshInterpolationFEM<Dimension>& interp, const Vector<Real_wp>& U0,
			  VectReal_wp& trace_En, VectReal_wp& trace_Hn,
			  bool assemble, bool compute_H) const
  {
    interp.ComputeEnHnOnBoundary(*this, U0, trace_En, trace_Hn, assemble, compute_H);
  }


  //! computes E \times n and H \times n (for Maxwell's equations) on quadrature points of a surface
  template<class Dimension>
  inline void VarProblem<Dimension>
  ::ComputeEnHnOnBoundary(const MeshInterpolationFEM<Dimension>& interp, const VectComplex_wp& U0,
			  VectComplex_wp& trace_En, VectComplex_wp& trace_Hn,
			  bool assemble, bool compute_H) const
  {
    interp.ComputeEnHnOnBoundary(*this, U0, trace_En, trace_Hn, assemble, compute_H);
  }
    
}
  
#define MONTJOIE_FILE_VAR_PROBLEM_INLINE_CXX
#endif
