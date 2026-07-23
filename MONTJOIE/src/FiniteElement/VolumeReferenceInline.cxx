#ifndef MONTJOIE_FILE_VOLUME_REFERENCE_INLINE_CXX

namespace Montjoie
{
  //! default constructor
  template<int type>
  inline VolumeReference<type>::VolumeReference(ElementGeomReference<Dimension3>& elt)
    : ElementReference_Dim<Dimension3>(elt)
  {
    nb_dof_tri = 0;
    nb_dof_quad = 0;
  }
  

  //! not used
  template<int type>
  inline void VolumeReference<type>
  ::SetInputData(const string& keyword, const Vector<string>& param)
  {
    // nothing (sigh)
  }
  
  
  /************************
   * Convenient functions *
   ************************/
  
  
  //! returns dof number on face
  template<int type>
  inline int VolumeReference<type>::GetLocalNumber(int num_face_loc, int num_dof) const
  {
    return FacesDof(num_dof, num_face_loc);
  }

  
  //! returns the number of dofs points on a face
  template<int type> 
  inline int VolumeReference<type>::GetNbDofBoundary(int i) const
  {
    if (this->elt_geom.IsLocalFaceQuadrilateral(i))
      return nb_dof_quad;
    
    return nb_dof_tri;
  }


  template<int type>
  inline void VolumeReference<type>
  ::ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const
  {
    const ElementReference<Dimension3, type>& elt
      = static_cast<const ElementReference<Dimension3, type>&>(*this);
    
    this->ProjectQuadratureToDofGen(elt, Equad, Edof);
  }

  
  template<int type>
  inline void VolumeReference<type>
  ::ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const
  {
    const ElementReference<Dimension3, type>& elt
      = static_cast<const ElementReference<Dimension3, type>&>(*this);
    
    this->ProjectQuadratureToDofGen(elt, Equad, Edof);
  }


  template<int type>
  inline void VolumeReference<type>
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    this->ComputeNodalValuesGen(Un, Unode);
  }

  template<int type>
  inline void VolumeReference<type>
  ::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    this->ComputeNodalValuesGen(Un, Unode);
  }


  /****************
   * H^1 elements *
   ****************/
  
  
#ifdef MONTJOIE_WITH_NODAL_H1
  inline ElementReference<Dimension3, 1>
  ::ElementReference(ElementGeomReference<Dimension3>& elt) : VolumeReference<1>(elt)
  {     element_tri_surf = NULL; element_quad_surf = NULL; }

  inline void ElementReference<Dimension3, 1>
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    this->ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    this->ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    this->ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    this->ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }
  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  inline void ElementReference<Dimension3, 1>
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    contrib = feval;
  }
  
  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  inline void ElementReference<Dimension3, 1>
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    contrib = feval;
  }

  inline void ElementReference<Dimension3, 1>
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
  }

  inline void ElementReference<Dimension3, 1>::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    this->ApplyRhGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 1>::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    this->ApplyRhGen(Vh, Uh);
  }
    
  inline void ElementReference<Dimension3, 1>
  ::ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp& Vz) const
  {
    this->ApplyRhSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy, VectComplex_wp& Vz) const
  {
    this->ApplyRhSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyRhQuadratureSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp& Vz) const
  {
    this->ApplyRhQuadratureSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyRhQuadratureSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy, VectComplex_wp& Vz) const
  {
    this->ApplyRhQuadratureSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    this->ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    this->ApplyRhTransposeGen(Uh, Vh);
  }

  inline void ElementReference<Dimension3, 1>
  ::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    this->ApplyChGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    this->ApplyChGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    this->ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    this->ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 1>::ApplyConstantRh(const VectReal_wp& U, VectReal_wp& V) const
  {
    this->ApplyConstantRhGen(U, V);
  }
  
  inline void ElementReference<Dimension3, 1>::ApplyConstantRh(const VectComplex_wp& U, VectComplex_wp& V) const
  {
    this->ApplyConstantRhGen(U, V);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp& Vz) const
  {
    this->ApplyConstantRhSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy, VectComplex_wp& Vz) const
  {
    this->ApplyConstantRhSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyConstantRhTranspose(const VectReal_wp& U, VectReal_wp& V) const
  {
    this->ApplyConstantRhTransposeGen(U, V);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::ApplyConstantRhTranspose(const VectComplex_wp& U, VectComplex_wp& V) const
  {
    this->ApplyConstantRhTransposeGen(U, V);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::AddConstantMassMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::AddConstantMassMatrix(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, mass, A);
  }

  inline void ElementReference<Dimension3, 1>
  ::AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
			  const TinyMatrix<Complex_wp, General, 3, 3>& C,
			  const R3_Complex_wp& D, const R3_Complex_wp& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }

  inline void ElementReference<Dimension3, 1>
  ::AddConstantElemMatrix(int m, int n, const Real_wp& mass,
			  const TinyMatrix<Real_wp, General, 3, 3>& C,
			  const R3& D, const R3& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			  const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
			  const Vector<R3>& D, const Vector<R3>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    this->AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			  const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
			  const Vector<R3_Complex_wp>& D, const Vector<R3_Complex_wp>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
    
  inline void ElementReference<Dimension3, 1>
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Real_wp>& A, VirtualMatrix<Real_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
  
  inline void ElementReference<Dimension3, 1>
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Complex_wp>& A, VirtualMatrix<Complex_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }

  inline void ElementReference<Dimension3, 1>
  ::ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
                             const ElementReference_Dim<Dimension3>& FaceCoarse,
                             const ElementReference_Dim<Dimension3>& FaceFine) const
  {
    FiniteElementH1<Dimension3>::ComputeLocalProlongation(proj, LocalProlongation, FaceCoarse, FaceFine);
  }
#endif

  
  /********************
   * H(curl) elements *
   ********************/
  
  
#ifdef MONTJOIE_WITH_NODAL_HCURL
  inline void ElementReference<Dimension3, 2>
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::ComputeCurlBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::ComputeCurlBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void ElementReference<Dimension3, 2>
  ::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void ElementReference<Dimension3, 2>
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void ElementReference<Dimension3, 2>::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 2>::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void ElementReference<Dimension3, 2>
  ::AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& mass,
			       VirtualMatrix<Real_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& mass,
			       VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }

  inline void ElementReference<Dimension3, 2>
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& C,
			  VirtualMatrix<Real_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& C,
			  VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
			  VirtualMatrix<Real_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, C, mat);
  }

  inline void ElementReference<Dimension3, 2>
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
			  VirtualMatrix<Complex_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, C, mat);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::AddVariableStiffnessMatrix(int off_row, int off_col,
			       const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
			       VirtualMatrix<Real_wp>& mat) const
  {
    AddVariableStiffnessMatrixGen(off_row, off_col, C, mat);
  }
  
  inline void ElementReference<Dimension3, 2>
  ::AddVariableStiffnessMatrix(int off_row, int off_col,
			       const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
			       VirtualMatrix<Complex_wp>& mat) const
  {
    AddVariableStiffnessMatrixGen(off_row, off_col, C, mat);
  }

  inline VolumeHcurlReference::VolumeHcurlReference(ElementReference<Dimension3, 1>& elt,
						    ElementReference<Dimension3, 2>& elt_hc)
    : elt_H1(elt), elt_Hcurl(elt_hc)
  {
  }
#endif
  
  
  /********************
   * H(div) elements *
   ********************/

  
#ifdef MONTJOIE_WITH_NODAL_HDIV
  inline ElementReference<Dimension3, 3>
  ::ElementReference(ElementGeomReference<Dimension3>& elt) : VolumeReference<3>(elt)
  { element_tri_surf = NULL; element_quad_surf = NULL; }

  inline void ElementReference<Dimension3, 3>
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension3, 3>
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }

  inline void ElementReference<Dimension3, 3>
  ::ComputeDivBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeDivBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension3, 3>
  ::ComputeDivBoundaryRef(const VectComplex_wp& Un,
			  VectComplex_wp& Unode, int num_loc) const
  {
    ComputeDivBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void ElementReference<Dimension3, 3>
  ::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 3>
  ::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void ElementReference<Dimension3, 3>
  ::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
    
  inline void ElementReference<Dimension3, 3>
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void ElementReference<Dimension3, 3>
  ::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
    
  inline void ElementReference<Dimension3, 3>
  ::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 3>
  ::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void ElementReference<Dimension3, 3>
  ::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void ElementReference<Dimension3, 3>
  ::AddConstantStiffnessMatrix(int m, int n, const Real_wp& mass,
			       VirtualMatrix<Real_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
    
  inline void ElementReference<Dimension3, 3>
  ::AddConstantStiffnessMatrix(int m, int n, const Complex_wp& mass,
			       VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void ElementReference<Dimension3, 3>
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& C,
			  VirtualMatrix<Real_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void ElementReference<Dimension3, 3>
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& C,
			  VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void ElementReference<Dimension3, 3>
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
	    VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void ElementReference<Dimension3, 3>
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void ElementReference<Dimension3, 3>
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void ElementReference<Dimension3, 3>
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
#endif
    
} // end namespace

#define MONTJOIE_FILE_VOLUME_REFERENCE_INLINE_CXX
#endif
