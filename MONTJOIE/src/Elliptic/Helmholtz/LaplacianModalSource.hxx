#ifndef MONTJOIE_FILE_LAPLACIAN_MODAL_SOURCE_HXX

namespace Montjoie
{
  
  template<class Dimension>
  class ModalSourceLaplace_Dim
  {
    typedef typename Dimension::DimensionBoundary DimensionB;
    typedef typename Dimension::DimensionBoundary::VectR_N VectR_Nm1;
    
  private:
    const VarProblem<Dimension>& var_problem;
    const ModalSourceBoundary_Base& modal_source;
    
  public:
    template<class TypeEquation>
    ModalSourceLaplace_Dim(const EllipticProblem<TypeEquation>& var,
			   const ModalSourceBoundary_Base& mode);
    
    void EvaluateModeSection(const Mesh<Dimension>& mesh, int proc_mode,
			     const Vector<int>& IndexBoundary, const Vector<int>&,
			     const Vector<int>&, const Vector<VectR_Nm1>&,
			     Vector<Vector<VectReal_wp> >&, VectReal_wp&);
    
    void GetCoefficientBoundaryCondition(const VectReal_wp& eigenval, VectReal_wp& coef_eigenval,
					 const Real_wp& imped, int boundary_condition);
    
    void GetCoefficientBoundaryCondition(const VectReal_wp& eigenval, VectComplex_wp& coef_eigenval,
					 const Complex_wp& imped, int boundary_condition);
    
    virtual void CopyIndices(VarHelmholtz_Base<Real_wp, DimensionB>& var_eig) = 0;
    
  };

  //! source on a section of the mesh (modes of the section are involved)
  template<class Complexe, class Dimension>
  class ModalSourceLaplace : public ModalSourceBoundary_Dim<Complexe, Dimension>,
			     public ModalSourceLaplace_Dim<Dimension>
  {
  private :
    typedef typename Dimension::DimensionBoundary DimensionB;
    
    const VarProblem<Dimension>& var_problem;
    const VarHelmholtz_Base<Complexe, Dimension>& var_laplace;

  public :
    template<class TypeEquation>
    ModalSourceLaplace(const EllipticProblem<TypeEquation>& var);
    
    void EvaluateMode(const Mesh<Dimension>& mesh, int proc_mode, const Vector<int>& IndexBoundary, 
		      const Vector<int>&, const Vector<int>& RotElement,
		      const Vector<typename DimensionB::VectR_N>& LocalQuadPoints,
		      Vector<Vector<Vector<Complexe> > >& eval);
    
    void CopyIndices(VarHelmholtz_Base<Real_wp, DimensionB>& var_eig);
    
  };

}

#define MONTJOIE_FILE_LAPLACIAN_MODAL_SOURCE_HXX
#endif

