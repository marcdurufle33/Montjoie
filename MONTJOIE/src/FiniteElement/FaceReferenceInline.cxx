#ifndef MONTJOIE_FILE_FACE_REFERENCE_INLINE_CXX

#include "FaceGeomReferenceInline.cxx"

namespace Montjoie
{
  //! default constructor
  template<int type>
  inline FaceReference<type>::FaceReference(ElementGeomReference<Dimension2>& elt)
    : ElementReference_Dim<Dimension2>(elt)
  {
    element_surface = NULL;
  }
  
  
  /**********************
   * Convenient methods *
   **********************/
  
    
  //! returns number of the dof k of the edge j
  /*!
    \param[in] k local dof number on the edge
    \param[in] j edge number
   */
  template<int type>
  inline int FaceReference<type>::GetLocalNumber(int j, int k) const
  {
    return EdgesDof(k, j);
  }
  
  
  //! returns the number of dofs on the boundary j
  template<int type>
  inline int FaceReference<type>::GetNbDofBoundary(int j) const
  {
    return EdgesDof.GetM();
  }


  template<int type>
  inline void FaceReference<type>
  ::ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const
  {
    const ElementReference<Dimension2, type>& elt
      = static_cast<const ElementReference<Dimension2, type>&>(*this);
    
    this->ProjectQuadratureToDofGen(elt, Equad, Edof);
  }

  
  template<int type>
  inline void FaceReference<type>
  ::ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const
  {
    const ElementReference<Dimension2, type>& elt
      = static_cast<const ElementReference<Dimension2, type>&>(*this);
    
    this->ProjectQuadratureToDofGen(elt, Equad, Edof);
  }


  template<int type>
  inline void FaceReference<type>
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    this->ComputeNodalValuesGen(Un, Unode);
  }

  template<int type>
  inline void FaceReference<type>
  ::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    this->ComputeNodalValuesGen(Un, Unode);
  }

  
  /****************
   * H^1 elements *
   ****************/


#ifdef MONTJOIE_WITH_NODAL_H1
  inline ElementReference<Dimension2, 1>
  ::ElementReference(ElementGeomReference<Dimension2>& elt) : FaceReference<1>(elt)
  {}

  inline void ElementReference<Dimension2, 1>
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    this->ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
    
  inline void ElementReference<Dimension2, 1>
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    this->ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
    
  inline void ElementReference<Dimension2, 1>
  ::ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    this->ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    this->ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }
  
  //! projection from values on dof points to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  inline void ElementReference<Dimension2, 1>
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    contrib = feval;
  }
  

  //! projection from values on dof points to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */  
  inline void ElementReference<Dimension2, 1>
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    contrib = feval;
  }

  inline void ElementReference<Dimension2, 1>
  ::ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int num_loc) const
  {
    this->element_surface->ComputeProjectionDofRef(feval, contrib);
  }

  inline void ElementReference<Dimension2, 1>
  ::ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int num_loc) const
  {
    this->element_surface->ComputeProjectionDofRef(feval, contrib);
  }
  
  inline void ElementReference<Dimension2, 1>::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    this->ApplyRhGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension2, 1>::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    this->ApplyRhGen(Vh, Uh);
  }
    
  inline void ElementReference<Dimension2, 1>
  ::ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy) const
  {
    this->ApplyRhSplitGen(Uh, Vx, Vy);
  }
    
  inline void ElementReference<Dimension2, 1>
  ::ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy) const
  {
    this->ApplyRhSplitGen(Uh, Vx, Vy);
  }

  inline void ElementReference<Dimension2, 1>
  ::ApplyRhQuadratureSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy) const
  {
    this->ApplyRhQuadratureSplitGen(Uh, Vx, Vy);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::ApplyRhQuadratureSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy) const
  {
    this->ApplyRhQuadratureSplitGen(Uh, Vx, Vy);
  }

  inline void ElementReference<Dimension2, 1>
  ::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    this->ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    this->ApplyRhTransposeGen(Uh, Vh);
  }

  inline void ElementReference<Dimension2, 1>
  ::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    this->ApplyChGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    this->ApplyChGen(Vh, Uh);
  }

  inline void ElementReference<Dimension2, 1>
  ::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    this->ApplyChTransposeGen(Vh, Uh);
  }
    
  inline void ElementReference<Dimension2, 1>
  ::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    this->ApplyChTransposeGen(Vh, Uh);
  }


  inline void ElementReference<Dimension2, 1>::ApplyConstantRh(const VectReal_wp& U, VectReal_wp& V) const
  {
    this->ApplyConstantRhGen(U, V);
  }
  
  inline void ElementReference<Dimension2, 1>::ApplyConstantRh(const VectComplex_wp& U, VectComplex_wp& V) const
  {
    this->ApplyConstantRhGen(U, V);
  }
    
  inline void ElementReference<Dimension2, 1>
  ::ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy) const
  {
    this->ApplyConstantRhSplitGen(Uh, Vx, Vy);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy) const
  {
    this->ApplyConstantRhSplitGen(Uh, Vx, Vy);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::ApplyConstantRhTranspose(const VectReal_wp& U, VectReal_wp& V) const
  {
    this->ApplyConstantRhTransposeGen(U, V);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::ApplyConstantRhTranspose(const VectComplex_wp& U, VectComplex_wp& V) const
  {
    this->ApplyConstantRhTransposeGen(U, V);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::AddConstantMassMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::AddConstantMassMatrix(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, mass, A);
  }
    
  inline void ElementReference<Dimension2, 1>
  ::AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
			  const TinyMatrix<Complex_wp, General, 2, 2>& C,
			  const R2_Complex_wp& D, const R2_Complex_wp& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }

  inline void ElementReference<Dimension2, 1>
  ::AddConstantElemMatrix(int m, int n, const Real_wp& mass,
			  const TinyMatrix<Real_wp, General, 2, 2>& C,
			  const R2& D, const R2& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			  const Vector<TinyMatrix<Real_wp, General, 2, 2> >& C,
			  const Vector<R2>& D, const Vector<R2>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    this->AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			  const Vector<TinyMatrix<Complex_wp, General, 2, 2> >& C,
			  const Vector<R2_Complex_wp>& D, const Vector<R2_Complex_wp>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
    
  inline void ElementReference<Dimension2, 1>
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Real_wp>& A, VirtualMatrix<Real_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
  
  inline void ElementReference<Dimension2, 1>
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Complex_wp>& A, VirtualMatrix<Complex_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }

  inline void ElementReference<Dimension2, 1>
  ::ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
                             const ElementReference_Dim<Dimension2>& FaceCoarse,
                             const ElementReference_Dim<Dimension2>& FaceFine) const
  {
    FiniteElementH1<Dimension2>::ComputeLocalProlongation(proj, LocalProlongation, FaceCoarse, FaceFine);
  }
#endif


#ifdef MONTJOIE_WITH_NODAL_HCURL
  //! returns orientations associated with each dof
  inline const Vector<R2>& ElementReferenceType<Dimension2, 2>::GetOrientationDofs() const
  {
    return tangente_dof;
  }

  
  inline ElementReference<Dimension2, 2>
  ::ElementReference(ElementGeomReference<Dimension2>& elt) : FaceReference<2>(elt)
  { nb_dof_H = 0; }

  inline int ElementReference<Dimension2, 2>::GetNbDofH() const
  {
    return nb_dof_H; 
  }
  
  inline void ElementReference<Dimension2, 2>
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    this->ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension2, 2>
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    this->ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
    
  inline void ElementReference<Dimension2, 2>
  ::ComputeCurlBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    this->ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension2, 2>
  ::ComputeCurlBoundaryRef(const VectComplex_wp& Un,
			   VectComplex_wp& Unode, int num_loc) const
  {
    this->ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension2, 2>::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    this->ApplyRhGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension2, 2>::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    this->ApplyRhGen(Vh, Uh);
  }

  inline void ElementReference<Dimension2, 2>
  ::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    this->ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void ElementReference<Dimension2, 2>
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    this->ApplyRhTransposeGen(Uh, Vh);
  }

  inline void ElementReference<Dimension2, 2>::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    this->ApplyChGen(Vh, Uh);
  }

  inline void ElementReference<Dimension2, 2>::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    this->ApplyChGen(Vh, Uh);
  }

  inline void ElementReference<Dimension2, 2>::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    this->ApplyChTransposeGen(Vh, Uh);
  }

  inline void ElementReference<Dimension2, 2>
  ::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    this->ApplyChTransposeGen(Vh, Uh);
  }

  inline void ElementReference<Dimension2, 2>
  ::ApplyChH(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    Mlt(Ch_Hloc, Vh, Uh);
  }

  
  inline void ElementReference<Dimension2, 2>
  ::ApplyChH(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    Mlt(Ch_Hloc, Vh, Uh);
  }


  inline void ElementReference<Dimension2, 2>
  ::ApplyChTransposeH(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    Mlt(SeldonTrans, Ch_Hloc, Vh, Uh);
  }

  
  inline void ElementReference<Dimension2, 2>
  ::ApplyChTransposeH(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    Mlt(SeldonTrans, Ch_Hloc, Vh, Uh);
  }

  inline void ElementReference<Dimension2, 2>
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void ElementReference<Dimension2, 2>
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh,
		     VectComplex_wp& Vh, int r) const
  {
    this->ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void ElementReference<Dimension2, 2>
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
	    VectReal_wp& Vh, int r) const
  {
    this->ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void ElementReference<Dimension2, 2>
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    this->ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void ElementReference<Dimension2, 2>
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh,
			  VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void ElementReference<Dimension2, 2>
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh,
			  VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void ElementReference<Dimension2, 2>
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void ElementReference<Dimension2, 2>
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void ElementReference<Dimension2, 2>
  ::AddConstantStiffnessMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void ElementReference<Dimension2, 2>
  ::AddConstantStiffnessMatrix(int m, int n,
			       const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void ElementReference<Dimension2, 2>
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 2, 2>& C,
			  VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, C, A);
  }
    
  inline void ElementReference<Dimension2, 2>
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 2, 2>& C,
			  VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void ElementReference<Dimension2, 2>
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<TinyMatrix<Real_wp, General, 2, 2> >& C,
			  VirtualMatrix<Real_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, C, mat);
  }

  inline void ElementReference<Dimension2, 2>
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<TinyMatrix<Complex_wp, General, 2, 2> >& C,
			  VirtualMatrix<Complex_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, C, mat);
  }
    
  inline void ElementReference<Dimension2, 2>
  ::AddVariableStiffnessMatrix(int off_row, int off_col,
			       const Vector<Real_wp>& A, VirtualMatrix<Real_wp>& mat) const
  {
    this->AddVariableStiffnessMatrixGen(off_row, off_col, A, mat);
  }

  inline void ElementReference<Dimension2, 2>
  ::AddVariableStiffnessMatrix(int off_row, int off_col,
			       const Vector<Complex_wp>& A, VirtualMatrix<Complex_wp>& mat) const
  {
    this->AddVariableStiffnessMatrixGen(off_row, off_col, A, mat);
  }
#endif


#ifdef MONTJOIE_WITH_NODAL_HDIV
  inline ElementReference<Dimension2, 3>
  ::ElementReference(ElementGeomReference<Dimension2>& elt) : FaceReference<3>(elt)
  {}

  inline void ElementReference<Dimension2, 3>
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension2, 3>
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension2, 3>
  ::ComputeDivBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeDivBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension2, 3>
  ::ComputeDivBoundaryRef(const VectComplex_wp& Un,
			  VectComplex_wp& Unode, int num_loc) const
  {
    ComputeDivBoundaryGen(Un, Unode, num_loc);
  }

  inline void ElementReference<Dimension2, 3>
  ::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension2, 3>
  ::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void ElementReference<Dimension2, 3>
  ::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
    
  inline void ElementReference<Dimension2, 3>
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void ElementReference<Dimension2, 3>
  ::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
    
  inline void ElementReference<Dimension2, 3>
  ::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension2, 3>
  ::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension2, 3>
  ::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void ElementReference<Dimension2, 3>
  ::AddConstantStiffnessMatrix(int m, int n, const Real_wp& mass,
			       VirtualMatrix<Real_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
    
  inline void ElementReference<Dimension2, 3>
  ::AddConstantStiffnessMatrix(int m, int n, const Complex_wp& mass,
			       VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void ElementReference<Dimension2, 3>
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 2, 2>& C,
			  VirtualMatrix<Real_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void ElementReference<Dimension2, 3>
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 2, 2>& C,
			  VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void ElementReference<Dimension2, 3>
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
	    VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void ElementReference<Dimension2, 3>
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void ElementReference<Dimension2, 3>
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void ElementReference<Dimension2, 3>
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
#endif

}
  
#define MONTJOIE_FILE_FACE_REFERENCE_INLINE_CXX
#endif
