#ifndef MONTJOIE_FILE_AXISYM_HCURL_MAXWELL_HXX

namespace Montjoie
{
  //! class to solve time-harmonic Maxwell equations in axisymmetric domain and edge finite element
  template<class TypeEquation>
  class HarmonicMaxwellAxi_Hcurl : public HarmonicMaxwellAxi_Eq<TypeEquation>
  {
  public :
    typedef typename TypeEquation::Complexe Complexe;
    typedef Vector<Complexe> VectComplexe;
    
    //! offset for dofs for unknowns H
    IVect OffsetDofH;
    bool modified_formulation; //!< if true, Lacoste formulation is used
    //! coefficients used to compute elementary matrices
    Vector<Vector<TinyMatrix<Complexe, Symmetric, 2, 2> > > matmass_Hrz, matmass_Erz,
      matmass_Erz_m2, matmass_Erz_zz, matmass_Erz_r3, matmass_Hrz_r3;
    
    //! coefficients used to compute elementary matrices
    Vector<Vector<TinyVector<Complexe, 2> > > vecstiff_Erz, vecstiff_Eteta, vecstiff_Eteta_r2;
    //! coefficients used to compute elementary matrices
    Vector<Vector<Complexe> > matmass_Eteta, matmass_Hteta, matstiff_radius,
      matstiff_a10, matstiff_a11, matmass_Eteta_m0, matmass_Eteta_rr, matmass_Hteta_r3;
    
  public:
    HarmonicMaxwellAxi_Hcurl();
    
    void SetInputData(const string&, const Vector<string>&);
    void SetTypeEquation(const string&);
    void SetTypeElement(string& name_elt);

    void ComputeNumberOfDofs(); 
    void AllocateMassMatrices();
    void PerformOtherInitializations();
    void UpdateDirichlet(int n);
    
    void ComputeHcurlMassMatrix(int i, const ElementReference_Dim<Dimension2>& Fb);

    void GetLocalUnknownVectorE(const Vector<Complexe>& U, int iquad, const ElementReference<Dimension2, 2>& Fb_hcurl,
                                const ElementReference<Dimension2, 1>& Fb_h1, Vector<Complexe>& E, Vector<Complexe>& Eteta) const;

    void GetLocalUnknownVectorE(const Vector<Vector<Complexe> >& U, int iquad, const ElementReference<Dimension2, 2>& Fb_hcurl,
                                const ElementReference<Dimension2, 1>& Fb_h1, Vector<Complexe>& E, Vector<Complexe>& Eteta) const;
    
    void ComputeEnergy(const Vector<Complexe>& U0, VectReal_wp& normE, VectReal_wp& normH);
    
    void ComputeEpolar(const Vector<VectComplexe>& U0, int num_elem, const R2& point_loc, const R2&,
		       const Matrix2_2& dfjm1, TinyVector<Complexe, 3>& val_E,
                       const ElementReference<Dimension2, 2>& Fb_hcurl,
                       const ElementReference<Dimension2, 1>& Fb_h1) const;
    
    void ComputeHpolar(const Vector<VectComplexe>& U0, int num_elem, const R2& point_loc, const R2& point_glob,
		       const Matrix2_2& dfjm1, TinyVector<Complexe, 3>& val_H,
                       const ElementReference<Dimension2, 2>& Fb_hcurl,
                       const ElementReference<Dimension2, 1>& Fb_h1) const;

  protected:
    template<class T>
    void ComputeInterpolationUlocGen(const Vector<Vector<T> > & U0,
                                     const GridInterpolation<Dimension2> & var_interp,
                                     Vector<T>& trace_vec, Vector<T>& trace_grad_vec,
                                     const IVect& list_points, int nnz, bool compute_grad) const;

    void ComputeInterpolationUlocGen(const Vector<VectComplexe> & U0,
                                     const GridInterpolation<Dimension2> & var_interp,
                                     Vector<Complexe>& trace_vec, Vector<Complexe>& trace_grad_vec,
                                     const IVect& list_points, int nnz, bool compute_grad) const;

    void ComputeInterpolationUloc(const Vector<VectReal_wp> & U0,
				  const GridInterpolation<Dimension2> & var_interp,
				  Vector<Real_wp>& trace_vec, Vector<Real_wp>& trace_grad_vec,
				  const IVect& list_points, int nnz, bool compute_grad) const;
    
    void ComputeInterpolationUloc(const Vector<VectComplex_wp> & U0,
				  const GridInterpolation<Dimension2> & var_interp,
				  Vector<Complex_wp>& trace_vec, Vector<Complex_wp>& trace_grad_vec,
				  const IVect& list_points, int nnz, bool compute_grad) const;
    
  public:    
    template<class T>
    void AddVolumeSourceGen(const T& alpha, Vector<Vector<T> > & b_source,
                            Vector<VirtualSourceFEM<T, Dimension2>* >& f) const;
    
    void AddVolumeSourceGen(const Complexe& alpha, Vector<VectComplexe> & b_source,
                            Vector<VirtualSourceFEM<Complexe, Dimension2>* >& f) const;

    void AddVolumeSource(const Real_wp& alpha, Vector<VectReal_wp> & b_source,
			 Vector<VirtualSourceFEM<Real_wp, Dimension2>* >& f) const;

    void AddVolumeSource(const Complex_wp& alpha, Vector<VectComplex_wp> & b_source,
			 Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >& f) const;
    
    void AddDiracSource(const Real_wp& alpha, Vector<Vector<Real_wp> > & b_source,
                        Vector<VirtualSourceFEM<Real_wp, Dimension2>* >& f) const;

    void AddDiracSource(const Complex_wp& alpha, Vector<Vector<Complex_wp> > & b_source,
                        Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >& f) const;
    
    void ComputeSourceH_ComplementSchur(Vector<Complexe>& contrib_Hx, Vector<Complexe>& contrib_Hy,
                                        Vector<Complexe>& contrib_Hteta,
					Vector<Complexe>& contrib1, Vector<Complexe>& contrib2,
                                        const ElementReference<Dimension2, 2>& Fb1,
                                        const ElementReference<Dimension2, 1>& Fb2, int num_elem) const;
    
    void ComputeLocalMatrices_MixedFormulation(int iquad, Matrix<Complexe>& Ch, Matrix<Complexe>& Rh1,
					       Matrix<Complexe>& Rh2, Matrix<Complexe>& Dh1,
					       Matrix<Complexe>& Dh2, Matrix<Complexe>& Bh1,
					       Matrix<Complexe>& Bh2,
					       const ElementReference<Dimension2, 2>& Fb2,
                                               const ElementReference<Dimension2, 1>& Fb1) const;
    
    void EvaluateH_MixedFormulation(const Vector<Vector<Complexe> >& U0,
                                    Vector<Vector<Complexe> >& EvalH_Nodal) const;
    
  };

  
  // function implemented in file ElementaryMatrixMaxwellHcurlAxi.cxx
  template<class Complexe, class Prop, class Storage, class TypeEquation>
  void ComputeElementaryMatrix(int iquad, IVect& num_dof, Matrix<Complexe, Prop, Storage>& mat,
			       const GlobalGenericMatrix<Complexe>& nat_mat,
			       const EllipticProblem<TypeEquation>& vars,
                               const ElementReference<Dimension2, 2>& Fb_hcurl,
                               const ElementReference<Dimension2, 1>& Fb);
  
  //! class to solve time-harmonic Maxwell equations
  //! in axisymmetric domain and first family on quadrilaterals
  template<>
  class EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>
    : public HarmonicMaxwellAxi_Hcurl<HarmonicMaxwellEquation_HcurlAxi>
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
    IncidentWaveProjector<Complex_wp, Dimension2>*
    GetNewIncidentProjector(int n, const Vector<VectString>&,
			    IncidentWaveField<Complex_wp, Dimension2>& u_inc) const;
    
  };

}

#define MONTJOIE_FILE_AXISYM_HCURL_MAXWELL_HXX
#endif
