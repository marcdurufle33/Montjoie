#ifndef MONTJOIE_FILE_TETRAHEDRON_HCURL_SECOND_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's second family on tetrahedra
  /*!
    not stabilized, do not use
   */
  class TetrahedronHcurlSecondFamily : public TetrahedronReference<2>
  {
  protected :
    
    Matrix<R3> Value_Phi; //!< \f$ \varphi_i(\xi_j) \f$
    // MatrixFullR2 Value_Phi2D,Curl_Phi2D;
    Matrix<R3> Curl_Phi; //!< \f$ \nabla \times \varphi_i(\xi_j) \f$
    //! mass_matrixij = \f$ \int u_i v_j  (u_i and v_j \f$
    //! are components of vectors u and v) 
    Matrix<Real_wp> mass_matrix00, mass_matrix11, mass_matrix01,
      mass_matrix02, mass_matrix12, mass_matrix22;
    //! stiffness_matrixij = \f$ \int curl(u)_i curl(v)_j  
    //! (curl(u)_i \mbox{ and } curl(v)_j \f$ 
    //! are components of vectors curl(u) and curl(v)) 
    Matrix<Real_wp> stiffness_matrix00, stiffness_matrix11, stiffness_matrix01,
      stiffness_matrix02, stiffness_matrix12, stiffness_matrix22;
    // VectR3 normale_dof;
    // VectBool IsDependantDof; int nb_dof_dependant, nb_dof_independant;
    // Matrix<Real_wp> CoefDependancy;
    
    Matrix<bool> SignFaceDof; //!< signs of dofs placed on face
    IVect NodalDof; //!< for each dof, related nodal point 
    VectR3 DirectionDof; //!< for each dof, orientation on unit element
    Matrix<int> ListeDof_Node; //!< for each node, three related dofs
    Vector<Matrix3_3> AhDof; //!< matrix used to project on dof components
    Matrix<int> FacesDof_Rotation_Tri;
    
  public :
    TetrahedronHcurlSecondFamily();

    void ConstructFunctions();
    void ConstructNumberMap(NumberMap& map, int dg) const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();
    
    Real_wp ComputeMassInteraction(int i, int j, const Matrix3_3sym& coef) const;
    Real_wp ComputeStiffnessInteraction(int i, int j, const Matrix3_3sym& coef) const;
    
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeCurlPhiRef(const R3& ptloc, VectR3& dphi) const;
    
    const R3& GetValuePhi(int node, int num_point) const;
    R3 GetValuePhiOnBoundary(int, int, int) const;
    
    const R3& GetCurlPhi(int node, int num_point) const;
    R3 GetCurlPhiOnBoundary(int, int, int) const;
    
    template<class Vector1, class Vector2>
    void ComputeProjectionDofRef(const Vector1& feval, Vector2& res) const;
    
    template<class Vector1, class Vector2>
    void ComputeIntegralRef(const Vector1 & feval, Vector2& res) const;
    
    template<class Vector1, class Vector2>
    void ComputeIntegralCurlRef(const Vector1& feval, Vector2& res) const;
    
    template<class Vector1, class Vector2>
    void ComputeIntegralSurfaceRef(const Vector1 & feval, Vector2& res,
				    int num_loc) const;

    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceRef(const Vector1 & feval, Vector2& res,
                                        int num_loc) const;

    template<class Vector1, class Vector2>
    void ComputeIntegralSurfaceCurlRef(const Vector1 & feval, Vector2& res,
                                       int num_loc) const;
    
    template<class Vector1,class Vector2>
    void ComputeValueBoundaryRef(const Vector1& u_loc, Vector2& u_boundary,
				  int num_loc) const;
    
    //! not implemented
    template<class Vector1,class Vector2>
    void ComputeCurlBoundaryRef(const Vector1& u_loc, Vector2& u_boundary,
				 int num_loc) const;
    
    //! not implemented
    template<class Vector1,class Vector2>
    void ComputeNodalValuesRef(const Vector1& u_loc, Vector2& u_boundary) const;
    
    //! not implemented
    template<class Vector1,class Vector2>
    void ComputeNodalCurlRef(const Vector1& u_loc, Vector2& u_boundary) const;
    
  };
  
  ostream& operator <<(ostream& out, const TetrahedronHcurlSecondFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_HCURL_SECOND_FAMILY_HXX
#endif
