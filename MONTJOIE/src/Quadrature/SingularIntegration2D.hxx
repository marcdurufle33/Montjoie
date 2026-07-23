#ifndef MONTJOIE_FILE_SINGULAR_INTEGRAL_EQUATION_2D_HXX

namespace Montjoie
{
  template<class T, class Dim>
  class SingularDoubleQuadratureGalerkin_Base
  {
  };

  template<class Dimension>
  class EuclidianDistanceClass_Base
  {
    typedef typename Dimension::R_N R_N;
    
  public:
    virtual Real_wp GetDistance(const R_N& x, const R_N& y) const;
    virtual R_N GetDiff(const R_N& u, const R_N& v) const;
    
  };

  
  //! class for 2-D double integrals with a singularity |s-t|^\beta
  template<class T>
  class SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  {
  protected :
    //! quadrature rules used for close elements
    Globatto<Real_wp> quadrature_proche;
    //! eta parameter used to distinguish close elements from far elements.
    Real_wp param_eta_quadrature;
    //! pointer to finite element class
    ElementReference<Dimension1, 1>* fe_class;
    //! type of operator
    int type_operator;
    //! exponent of the singularity
    Real_wp beta;
    //! type of integration
    int type_integration;
    //! quadrature rules for 1-D integrals involved in Lenoir-Salles formulas
    Globatto<Real_wp> quadrature_joint;
    //! coefficients used in Lenoir-Salles formulas
    VectReal_wp InvBetaPlusN;
    //! quadrature rules for outer integral in joint elements
    Globatto<Real_wp> quadrature_outer_joint;
    //! quadrature rules for inner integral in joint elements
    Globatto<Real_wp> quadrature_inner_joint;
    
    //! values of basis functions on quadrature points
    Matrix<Real_wp> ValPhiFar, ValPhiClose, ValPhiJoint, ValPhiJointExt;
    //! for interpolation
    Matrix<Real_wp> QuadPhiClose, QuadPhiJoint, QuadPhiJointExt;
    //! gradients of basis functions on quadrature points
    Matrix<Real_wp> GradPhiFar, GradPhiClose, GradPhiJoint, GradPhiJointExt;

    //! for identical panels (and curved)
    Globatto<Real_wp> gauss_ext, gauss_int, quad_basis;
    Matrix<Real_wp> QuadPhiExt, ValPhiExt, GradPhiExt;
    Vector<Matrix<Real_wp> > GradPhiInt, ValPhiInt, QuadPhiInt;
    Vector<Vector<Real_wp> > PointsInt, WeightsInt;
    
    //! elementary matrix on the unit square [0, 1]^2
    Matrix<Real_wp> mat_elem_reference_identical;

    //! basis functions in canonical basis
    Vector<UnivariatePolynomial<T> > Phi;
    Vector<UnivariatePolynomial<T> > dPhi;
    
    void ComputeFarPanelIntegral(bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
				 const VectR2& PointsQuadI, const VectR2& PointsQuadJ,
				 const VectR2& NormaleQuadI, const VectR2& NormaleQuadJ,
				 const EuclidianDistanceClass_Base<Dimension2>& var, Matrix<Real_wp>& mat_elem,
				 Matrix<Real_wp>&, Matrix<Real_wp>&);

    void ComputeClosePanelIntegral(bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
				   const VectR2& PointsQuadI, const VectR2& PointsQuadJ,
				   const VectR2& NormaleQuadI, const VectR2& NormaleQuadJ,
				   const EuclidianDistanceClass_Base<Dimension2>& var, Matrix<Real_wp>& mat_elem,
				   Matrix<Real_wp>&, Matrix<Real_wp>&);

    void ComputeIdenticalPanelIntegral(bool curved, const VectReal_wp& Dsi, const VectR2& PointsQuadI,
				       const VectR2& NormaleQuadI,
				       const EuclidianDistanceClass_Base<Dimension2>& var,
				       Matrix<Real_wp>& mat_elem);

    void ComputeJointPanelIntegral(bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
				   const VectR2& PointsQuadI, const VectR2& PointsQuadJ,
				   const VectR2& NormaleQuadI, const VectR2& NormaleQuadJ,
				   const R2& u, const R2& v,
				   const EuclidianDistanceClass_Base<Dimension2>& var, Matrix<Real_wp>& mat_elem,
				   Matrix<Real_wp>&, Matrix<Real_wp>&);
    
  public :    
    //! available methods of integration
    /*!
      Lenoir-Salles : reduction strategy such that double integrals are reduced to simple integrals without singularity
      Gauss-Squared : change of variables to remove a singularity in 1/sqrt(x) for the inner integral. for the outer integral, a singularity in sqrt(x) or sqrt(1-x) is removed
     */
    enum{LENOIR_SALLES, GAUSS_SQUARED};

    //! available operators
    enum {OPERATOR_GRAD, OPERATOR_DS, OPERATOR_DIFF};

    SingularDoubleQuadratureGalerkin_Base();
    
    void SetOperator(int type, const Real_wp& b);    
    void SetAnalyticalIntegration(const Real_wp& eta, int rp, int rj);
    void SetNumericalIntegration(const Real_wp& eta, int rp, int rj, int rj_ext, int n_int, int n_ext);
    
    void SetFiniteElement(const ElementReference<Dimension1, 1>& elt);

    template<class MatrixJacob>
    void ComputeElemMatrix(const VectR2& PtsAi, const VectR2& PtsAj, 
			   bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			   const VectR2& PointsNodalI, const VectR2& PointsNodalJ,
			   const VectR2& PointsQuadI, const VectR2& PointsQuadJ,
			   const VectR2& NormaleNodalI, const VectR2& NormaleNodalJ,
			   const VectR2& NormaleQuadI, const VectR2& NormaleQuadJ,
			   const MatrixJacob& MatJacobNodalI, const MatrixJacob& MatJacobNodalJ,
			   const MatrixJacob& MatJacobI, const MatrixJacob& MatJacobJ,
			   const EuclidianDistanceClass_Base<Dimension2>& var, Matrix<Real_wp>& mat_elem,
			   Matrix<Real_wp>&, Matrix<Real_wp>&);
    
    void ComputeMassMatrix(const VectR2& PtsAi, bool curved,
			   const VectReal_wp& Dsi, Matrix<Real_wp>& mat_elem);

    template<class MatrixJacob>
    void ComputeStiffnessMatrix(const VectR2& PtsAi, bool curved,
				const VectReal_wp& Dsi, const VectR2& NormaleQuad,
				const MatrixJacob&, const Complex_wp& kinf, const Real_wp& epsilon,
				Matrix<Complex_wp>& mat_elem);
    
  };
  
}

#define MONTJOIE_FILE_SINGULAR_INTEGRAL_EQUATION_2D_HXX
#endif
