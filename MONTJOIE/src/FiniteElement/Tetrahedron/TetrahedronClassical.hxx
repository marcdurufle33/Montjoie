#ifndef MONTJOIE_FILE_TETRAHEDRON_CLASSICAL_HXX

namespace Montjoie
{
  
  //! nodal finite element for tetrahedra
  /*!
    Interpolation points are Hesthaven's points
    Basis functions are generating the polynomial space P_r
   */
  class TetrahedronClassical : public TetrahedronReference<1>
  {
  protected :
    
    Array3D<int> NumFct3D; VectReal_wp InvWeightFct;
    Matrix<Real_wp> InverseBasisVDM;
    TriangleClassical* function_basis_tri;
    TinyVector<Matrix<R3>, 4> GradientPhi_Boundary;
    
  public :
    TetrahedronClassical();
    size_t GetMemorySize() const;
  
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
  protected :
    void ConstructFunctions();
    
    template<class Vector1, class Vector2>
    void ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;

    template<class Vector1, class Vector2>
    void ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh,
			 Vector2& Vh, int r = 0) const;
    
  public :
    void ComputeValuesPhiRef(const R3& x, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R3& x, VectR3& grad_phi) const;
    
    void PickNearDofs(int pos, const VectBool& DofUsed, IVect& ListeDof, int nb_dof) const;

    void ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r = 0) const;
    
    void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r = 0) const;
    
    void ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r = 0) const;
    
    void ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		      VectReal_wp& Vh, int r = 0) const;
    
    void ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		      VectComplex_wp& Vh, int r = 0) const;
        
    friend ostream& operator <<(ostream& out, const TetrahedronClassical& e);
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_CLASSICAL_HXX
#endif
