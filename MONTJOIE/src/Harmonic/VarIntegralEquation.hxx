#ifndef MONTJOIE_FILE_VAR_INTEGRAL_EQUATION_HXX

namespace Montjoie
{
  //! base class for integral equation parameters
  template<class TypeElement, class TypeEquation>
  class VarIntegralEquation_Base
  {
  protected :
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeElement::Dimension Dimension;
    
    // WARNING 
    // integral equation is supposed to be placed on the absorbing boundary condition
    
    //! order of integration for regular integrals, singular (for x and y)
    int order_integration_regular, order_integrationX, order_integrationY;
    //! quadrature method used for singular integrals (GAUSS, DUFFY, POLAR ...)
    int singular_integration;
    //! number of dofs of integral equation
    int nodl_integral; 
    //! integral formulation used (EFIE, MFIE or CFIE)
    int type_integral_formulation;
    Complexe alpha_efie; //!< coefficient for CFIE
    Complexe beta_mfie; //!< coefficient for CFIE
    
  public :
    //! type of coupling between BIE (Boundary Integral Equation) and FEM (Finite Element Method)
    int coupling_integral_equation; 

    //! different types of coupling
    enum {FEM_ONLY, BIE_FEM, BIE_ONLY};
    
    //! integral equation formulations
    /*!
      EFIE = Electric Field Integral Equation
      MFIE = Magnetic Field Integral Equation 
      CFIE = Combined Field Integral Equation
    */
    enum {EFIE_FORMULATION, MFIE_FORMULATION, CFIE_FORMULATION};

    SurfacicMesh<Dimension> mesh_ie; //!< surfacic mesh associated with integral equation 
    string name_file_mesh3D_ie; //!< file name where currents will be stored
    
  protected :

    typename TypeEquation::Matrix_Schur schur_matrix; //!< schur matrix for IE-FEM coupling
    int condition_integral_equation;
    //!< boundary condition of integral equation (Dirichlet, Neumann ...)
    Complexe beta_imped_integral;
    
    //! arrays with face numbers and dofs numbers, related to integral equation
    IVect ListeFace_IE, ListeDof_IE, Inverse_IndexFace_IE, Inverse_IndexDof_IE;

  public :
    VarIntegralEquation_Base()
    {
      order_integration_regular = 0;
      order_integrationX = 0;
      order_integrationY = 0;
      singular_integration = 0;
      nodl_integral = 0;
      type_integral_formulation = CFIE_FORMULATION;
      alpha_efie = 0.5; beta_mfie = 0.5;
      coupling_integral_equation = FEM_ONLY;
      name_file_mesh3D_ie = "current3D";
      condition_integral_equation = 1;
      beta_imped_integral = 0;
      
    }
    
    //! overloaded
    void GetDofsBoundary(int ref, IVect& ListeFace, IVect& ListeDof,
			  IVect& Inverse_IndexFace, IVect& Inverse_IndexDof) {}
    
    //! overloaded
    template<class MatrixIntegral>
    void ComputeMatrixIntegralEquation(MatrixIntegral& mat_integral) {}
    
    //! overloaded
    void ComputeRightHandSideIntegral(Vector<Complexe>& b) {}
    //! overloaded
    void ComputeOutputMeshIE(const string& file_name) {}

  };

} // namespace Montjoie

#define MONTJOIE_FILE_VAR_INTEGRAL_EQUATION_HXX
#endif

