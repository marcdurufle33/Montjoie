#ifndef MONTJOIE_FILE_FINITE_ELEMENT_HDIV_HXX

namespace Montjoie
{
  
  //! class to take into account transformation F_i on nodal finite element
  /*!
    This class defines effect of transform F_i
    for the integrands or expression of basis functions,
    in the case of nodal finite element :
    \f$ \varphi \circ F_i \, = \, \hat{\varphi} \f$
    and
    \f$ \nabla \varphi \, = \, DF_i^{*-1} \hat{\nabla} \hat{\varphi} \f$
   */
  template<class Dimension>
  class FiniteElementHdiv
  {
  public :
    enum{nb_components_u = Dimension::dim_N, nb_components_grad = 1};
    
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    typedef typename Dimension::MatrixN_N MatrixN_N; //!< 2x2 or 3x3 matrix
    
    typedef ElementReference<typename Dimension::DimensionBoundary, 1> ElementReferenceTrace;

    inline static int GetNatureElement() { return 3; }
    
    // computes res(i) = \phi_i(point_loc) for all basis function \phi_i
    void ComputeValuesPhi(const R_N& point_loc, Vector<R_N>& res, const MatrixN_N& dfjm1,
			  const MeshNumbering<Dimension>& mesh_num, int nquad) const;
    
    // computes res(i) = \nabla \cdot \phi_i (point_loc)
    void ComputeValuesGradientPhi(const R_N& point_loc,
                                  Vector<TinyVector<Real_wp, 1> >& res, const MatrixN_N& dfjm1,
				  const MeshNumbering<Dimension>& mesh_num, int nquad) const;
    
    template<class T>
    static void GetCurlFromGradient(const Vector<Vector<T> >& dU,
				    Vector<Vector<T> >& gradU);
    
    // computes ValuePhi_Boundary(i,k) = \phi_i(\xi_k)
    // where \xi_k is a quadrature point on the edge num_loc
    void ComputeValuesPhiBoundary(Matrix<R_N>& ValuePhi_Boundary,
				  const SetMatrices<Dimension>& MatricesElem,
				  const MeshNumbering<Dimension>& mesh_num,
				  int nquad, int num_loc) const;
    
    // computes Gradient_Boundary(i,k) = d\phi_i / dx (\xi_k),
    // where \xi_k is a quadrature point on the edge num_loc
    void ComputeValuesGradientPhiBoundary(Matrix<TinyVector<Real_wp, 1> >& GradientPhi_Boundary,
					  const SetMatrices<Dimension>& MatricesElem,
					  const MeshNumbering<Dimension>& mesh_num, 
					  int nquad, int num_loc) const;
    
    template<class T>
    void ComputeQuadratureValues(const SetMatrices<Dimension>& MatricesElem,
				 const Vector<Vector<T> >& Uloc,
                                 Vector<Vector<T> >& Uloc_node,
				 Vector<Vector<T> >& dUloc_node, bool compute_u, bool compute_grad,
				 const Mesh<Dimension>& mesh, int iquad) const;
    
    template<class T>
    void ComputeNodalValues(const SetMatrices<Dimension>& MatricesElem,
			    const Vector<Vector<T> >& Uloc, Vector<Vector<T> >& Uloc_node,
			    const Mesh<Dimension>& mesh, int nquad) const;
    
    template<class T>
    void ComputeProjectionDof(const SetMatrices<Dimension>& MatricesElem,
			      Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
			      const MeshNumbering<Dimension>& mesh_num, int n) const;

    template<class T>
    void ComputeProjectionSurfaceDof(const SetMatrices<Dimension>& MatricesElem,
				     Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
				     const MeshNumbering<Dimension>& mesh_num, int n, int num_loc) const;
    
    template<class T>
    void ComputeIntegral(const SetMatrices<Dimension>& mat,
			 Vector<Vector<T> >& feval, Vector<Vector<T> > & res,
			 const MeshNumbering<Dimension>& mesh_num, int nquad) const;
    
    template<class T>
    void ComputeIntegralGradient(const SetMatrices<Dimension>& mat,
				 Vector<Vector<T> >& feval,
				 Vector<Vector<T> >& res,
				 const MeshNumbering<Dimension>& mesh_num, int nquad) const;
    
    template<class T>
    void ComputeIntegralSurface(const SetMatrices<Dimension>& mat,
				Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
				const MeshNumbering<Dimension>& mesh_num,
				int n, int num_loc) const;

    template<class T>
    void ComputeIntegralSurfaceHDG(const SetMatrices<Dimension>& mat,
				   Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
				   const MeshNumbering<Dimension>& mesh_num,
				   int n, int num_loc) const {}

    template<class T>
    void ComputeGaussIntegralSurface(const Vector<typename Dimension::MatrixN_N>& dfjm1,
				     const VectReal_wp& weights, const VectReal_wp& ds,
				     Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
				     const MeshNumbering<Dimension>& mesh_num,
				     int n, int num_loc) const;
    
    template<class T>
    void ComputeIntegralSurfaceGradient(const SetMatrices<Dimension>& mat,
					Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
					const MeshNumbering<Dimension>& mesh_num,
					int n, int num_loc) const;
    
    template<class T>
    void ComputeValueBoundary(const SetMatrices<Dimension>& mat,
			      const Vector<Vector<T> >& feval,
                              Vector<Vector<TinyVector<T, Dimension::dim_N> > >& res,
			      const Mesh<Dimension>&, int n, int num_loc) const;
    
    template<class T>
    void ComputeGradientBoundary(const SetMatrices<Dimension>& mat,
				 const Vector<Vector<T> >& feval,
				 Vector<Vector<TinyVector<T, 1> > >& res,
                                 const Mesh<Dimension>&,
                                 int n, int num_loc) const;
    
    template<class T>
    void ComputeValueNodalBoundary(const Vector<T>& Uloc, Vector<T>& Uboundary,
                                   int num_loc) const;
        
    template<class T>
    static void ApplyRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
				    TinyVector<T, Dimension::dim_N>& val_u);
    
    template<class T>
    static void ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
					   TinyVector<T, Dimension::dim_N>& val_u);

    template<class T>
    static void ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
					   Vector<T>& val_u);
    
  };
  
}

#define MONTJOIE_FILE_FINITE_ELEMENT_HDIV_HXX
#endif
