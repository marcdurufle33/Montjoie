#ifndef MONTJOIE_FILE_FINITE_ELEMENT_HCURL3D_HXX

namespace Montjoie
{
  
  //! class for integration on 3D edge finite element
  /*!
    This class defines effect of transform F_i
    for the integrands or expression of basis functions,
    in the case of 3-D edge finite element :
    \f$ DF_i^* varphi \circ F_i \, = \, \hat{\varphi} \f$
   */
  template<>
  class FiniteElementHcurl<Dimension3>
  {
  public :
    enum{nb_components_u = 3, nb_components_grad = 3};

    typedef ElementReference<Dimension2, 2> ElementReferenceTrace;

    inline static int GetNatureElement() { return 2; }
    
    // computes res(i) = \phi_i(point_loc) for all basis function \phi_i
    void ComputeValuesPhi(const R3& point_loc, Vector<R3>& res, const Matrix3_3& dfjm1,
			  const MeshNumbering<Dimension3>& mesh_num, int nquad) const;
    
    // computes res(i) = \nabla \phi_i (point_loc)
    void ComputeValuesGradientPhi(const R3& point_loc, Vector<R3>& res, const Matrix3_3& dfjm1,
				  const MeshNumbering<Dimension3>& mesh_num, int nquad) const;
    
    template<class T>
    static void GetCurlFromGradient(const Vector<Vector<T> >& dU,
                                    Vector<Vector<T> >& curlU);
    
    // computes ValuePhi_Boundary(i,k) = \phi_i(\xi_k)
    // where \xi_k is a quadrature point on the edge num_loc
    void ComputeValuesPhiBoundary(Matrix<R3>& ValuePhi_Boundary,
                                  const SetMatrices<Dimension3>& MatricesElem,
                                  const MeshNumbering<Dimension3>& mesh_num,
				  int nquad, int num_loc) const;
    
    // computes GradientPhi_Boundary(i,k) = d\phi_i / ds (\xi_k),
    // where \xi_k is a quadrature point on the edge num_loc
    void ComputeValuesGradientPhiBoundary(Matrix<R3>& GradientPhi_Boundary,
                                          const SetMatrices<Dimension3>& MatricesElem,
                                          const MeshNumbering<Dimension3>& mesh_num, int nquad,
                                          int num_loc) const;
    
    template<class T>
    void ComputeNodalValues(const SetMatrices<Dimension3>& MatricesElem,
			    const Vector<Vector<T> >& Uloc, Vector<Vector<T> >& Uloc_node,
                            const Mesh<Dimension3>& mesh, int i) const;
    
    template<class T>
    void ComputeQuadratureValues(const SetMatrices<Dimension3>& MatricesElem,
				 const Vector<Vector<T> >& Uloc,
                                 Vector<Vector<T> >& Uloc_node,
				 Vector<Vector<T> >& dUloc_node, bool compute_u, bool compute_curl,
				 const Mesh<Dimension3>& mesh, int iquad) const;

    template<class T>
    void ComputeProjectionDof(const SetMatrices<Dimension3>& MatricesElem,
			      Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
			      const MeshNumbering<Dimension3>& mesh_num, int n) const;

    template<class T>
    void ComputeProjectionSurfaceDof(const SetMatrices<Dimension3>& MatricesElem,
				     Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
				     const MeshNumbering<Dimension3>& mesh_num, int n, int nloc) const;
    
    // computes integral \int f \phi_i  for all test functions phi_i of the element
    template<class T>
    void ComputeIntegral(const SetMatrices<Dimension3>& mat,
			 Vector<Vector<T> >& feval, Vector<Vector<T> > & res,
			 const MeshNumbering<Dimension3>& mesh_num, int nquad) const;
    
    // computes integral \int f \nabla \times \phi_i  for all test functions phi_i of the element
    template<class T>
    void ComputeIntegralGradient(const SetMatrices<Dimension3>& mat,
				 Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
                                 const MeshNumbering<Dimension3>& mesh_num, int nquad) const;
    
    // computes integral \int_\Gamma f \phi_i dx
    //  (surfacic integrals) for all test functions phi_i of the element
    template<class T>
    void ComputeGaussIntegralSurface(const Vector<Matrix3_3>& dfjm1,
                                     const VectReal_wp& weights, const VectReal_wp& ds,
                                     Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
                                     const MeshNumbering<Dimension3>& mesh_num,
                                     int n, int num_loc) const;
    
    // computes integral \int_\Gamma f \phi_i dx 
    // (surfacic integrals) for all test functions phi_i of the element
    template<class T>
    void ComputeIntegralSurface(const SetMatrices<Dimension3>& mat,
				Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
                                const MeshNumbering<Dimension3>& mesh_num,
                                int n, int num_loc) const;

    template<class T>
    void ComputeIntegralSurfaceHDG(const SetMatrices<Dimension3>& mat,
				   Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
                                   const MeshNumbering<Dimension3>& mesh_num,
				   int n, int num_loc) const;
    
    template<class T>
    void ComputeIntegralSurfaceGradient(const SetMatrices<Dimension3>& mat,
					Vector<Vector<T> >& feval,
                                        Vector<Vector<T> >& res,
					const MeshNumbering<Dimension3>& mesh_num,
                                        int n, int num_loc) const;
    
    template<class T>
    void ComputeValueBoundary(const SetMatrices<Dimension3>& mat,
			      const Vector<Vector<T> >& feval, Vector<Vector<TinyVector<T, 3> > >& res,
			      const Mesh<Dimension3>& mesh, int n, int num_loc) const;
    
    template<class T>
    void ComputeGradientBoundary(const SetMatrices<Dimension3>& mat,
				 const Vector<Vector<T> >& feval, Vector<Vector<TinyVector<T, 3> > >& res,
				 const Mesh<Dimension3>& mesh, int n, int num_loc) const;
    
    template<class T>
    void ComputeValueNodalBoundary(const Vector<T>& Uloc, Vector<T>& Uboundary,
                                   int num_loc) const;
    

    template<class T>
    static void ApplyRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
                                    TinyVector<T, 3>& val_u);

    template<class T>
    static void ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
                                           TinyVector<T, 3>& val_u);
    
    template<class T>
    static void ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
                                           Vector<T>& val_u);
    
  };

}

#define MONTJOIE_FILE_FINITE_ELEMENT_HCURL3D_HXX
#endif
