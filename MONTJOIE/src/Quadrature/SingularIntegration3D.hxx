#ifndef MONTJOIE_FILE_SINGULAR_INTEGRAL_EQUATION_3D_HXX

#include "SingularIntegration2D.hxx"

namespace Montjoie
{
  //! class for 3-D double integrals with a singularity |s-t|^\beta
  template<class T>
  class SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  {
  protected :
    //! eta parameter used to distinguish close elements from far elements.
    Real_wp param_eta_quadrature;
    //! pointer to triangular finite element class
    ElementReference<Dimension2, 1>* fe_tri;
    //! pointer to quadrangular finite element class
    ElementReference<Dimension2, 1>* fe_quad;
    //! type of operator
    int type_operator;
    //! exponent of the singularity
    Real_wp beta;
    //! type of integration
    int type_integration;

    //! quadrature rule for triangles
    VectR2 PointsTri; VectReal_wp WeightsTri;
    
    //! quadrature rule for edges
    Globatto<Real_wp> lob;

    Matrix<R2> GradPhiTriFar, GradPhiQuadFar;
    Matrix<Real_wp> ValPhiTriFar, ValPhiQuadFar;
    Vector<Real_wp> WeightsTriFar, WeightsQuadFar;

    TinyVector<VectReal_wp, 3> WeightsVertexTri;
    TinyVector<VectR2, 3> PointsVertexTri;
    TinyVector<Matrix<Real_wp>, 3> NodalPhiTriVertex, ValPhiTriVertex;
    TinyVector<Matrix<R2>, 3> GradPhiTriVertex;

    TinyVector<VectReal_wp, 4> WeightsVertexQuad;
    TinyVector<VectR2, 4> PointsVertexQuad;
    TinyVector<Matrix<Real_wp>, 4> NodalPhiQuadVertex, ValPhiQuadVertex;
    TinyVector<Matrix<R2>, 4> GradPhiQuadVertex;

    TinyVector<VectReal_wp, 3> WeightsEdgeTri;
    TinyVector<VectR2, 3> PointsEdgeTri;
    TinyVector<Matrix<Real_wp>, 3> NodalPhiTriEdge, ValPhiTriEdge;
    TinyVector<Matrix<R2>, 3> GradPhiTriEdge;    

    TinyVector<VectReal_wp, 4> WeightsEdgeQuad;
    TinyVector<VectR2, 4> PointsEdgeQuad;
    TinyVector<Matrix<Real_wp>, 4> NodalPhiQuadEdge, ValPhiQuadEdge;
    TinyVector<Matrix<R2>, 4> GradPhiQuadEdge;
    
    VectReal_wp WeightsCommonExtTri;
    VectR2 PointsCommonExtTri;
    Matrix<Real_wp> NodalPhiTriCommonExt, ValPhiTriCommonExt;
    Matrix<R2> GradPhiTriCommonExt;
    Vector<VectReal_wp> WeightsCommonTri;
    Vector<VectR2> PointsCommonTri;
    Vector<Matrix<Real_wp> > NodalPhiTriCommon, ValPhiTriCommon;
    Vector<Matrix<R2> > GradPhiTriCommon;

    VectReal_wp WeightsCommonExtQuad;
    VectR2 PointsCommonExtQuad;
    Matrix<Real_wp> NodalPhiQuadCommonExt, ValPhiQuadCommonExt;
    Matrix<R2> GradPhiQuadCommonExt;
    Vector<VectReal_wp> WeightsCommonQuad;
    Vector<VectR2> PointsCommonQuad;
    Vector<Matrix<Real_wp> > NodalPhiQuadCommon, ValPhiQuadCommon;
    Vector<Matrix<R2> > GradPhiQuadCommon;

    void GetOppositeEdge(int i, int& j, int& k);
    int GetOppositeVertex(int i, int j);
    void FillDofEdge(int n0, int n1, int order, int nb_vert, IVect& numI);
    
    template<class MatrixJacob>
    void ComputePanelIntegral(const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			      const VectReal_wp& WeightsI, const VectReal_wp& WeightsJ,
			      const VectR3& PointsQuadI, const VectR3& PointsQuadJ,
			      const VectR3& NormaleQuadI, const VectR3& NormaleQuadJ,
			      const MatrixJacob& MatJacobI, const MatrixJacob& MatJacobJ,
			      const Matrix<R2>& GradPhiI, const Matrix<R2>& GradPhiJ,
			      const Matrix<Real_wp>& ValPhiI, const Matrix<Real_wp>& ValPhiJ,
			      const EuclidianDistanceClass_Base<Dimension3>& var, Matrix<Real_wp>& mat_elem,
			      Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ,
			      const IVect& numI, const IVect& numJ);
    
    template<class MatrixJacob> void 
    ComputeFarPanelIntegral(int nb_i, int nb_j, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			    const VectR3& PointsQuadI, const VectR3& PointsQuadJ,
			    const VectR3& NormaleQuadI, const VectR3& NormaleQuadJ,
			    const MatrixJacob& MatJacobI, const MatrixJacob& MatJacobJ,
			    const EuclidianDistanceClass_Base<Dimension3>& var, Matrix<Real_wp>& mat_elem,
			    Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ);

    template<class MatrixJacob> void 
    ComputeVertexPanelIntegral(int pt_common_i, int pt_common_j, int nb_i, int nb_j,
			       bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			       const VectR3& PointsNodalI, const VectR3& PointsNodalJ,
			       const VectR3& NormaleNodalI, const VectR3& NormaleNodalJ,
			       const MatrixJacob& MatJacobNodalI, const MatrixJacob& MatJacobNodalJ,
			       const EuclidianDistanceClass_Base<Dimension3>& var, Matrix<Real_wp>& mat_elem,
			       Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ);
    
    template<class MatrixJacob>
    void ComputeEdgePanelIntegral(int pt_i0, int pt_i1, int pt_j0, int pt_j1, int nb_i, int nb_j,
				  bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
				  const VectR3& PointsNodalI, const VectR3& PointsNodalJ,
				  const VectR3& NormaleNodalI, const VectR3& NormaleNodalJ,
				  const MatrixJacob& MatJacobNodalI, const MatrixJacob& MatJacobNodalJ,
				  const EuclidianDistanceClass_Base<Dimension3>& var, Matrix<Real_wp>& mat_elem,
				  Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ);

    template<class MatrixJacob>
    void ComputeSamePanelIntegral(int nb_i, bool curved, const VectReal_wp& Dsi,
				  const VectR3& PointsNodalI, const VectR3& NormaleNodalI,
				  const MatrixJacob& MatJacobNodalI,
				  const EuclidianDistanceClass_Base<Dimension3>& var,
				  Matrix<Real_wp>& mat_elem,
				  Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ);
    
  public :    
    //! available methods of integration
    /*!
      Lenoir-Salles : reduction strategy such that integrals are reduced to integrals of lesser dimension without singularity
      
     */
    enum{LENOIR_SALLES, DUFFY};

    //! available operators
    enum {OPERATOR_GRAD, OPERATOR_DS, OPERATOR_DIFF};

    SingularDoubleQuadratureGalerkin_Base();
    
    void SetOperator(int type, const Real_wp& b);    
    void SetAnalyticalIntegration(const Real_wp& eta, int rp, int rj);
    void SetNumericalIntegration(const Real_wp& eta, int rp, int rj, int rj_ext, int n_int, int n_ext);
    
    void SetTriangularFiniteElement(const ElementReference<Dimension2, 1>& elt,
				    const VectR2&, const VectReal_wp&);
    
    void SetQuadrilateralFiniteElement(const ElementReference<Dimension2, 1>& elt,
				       const VectR2&, const VectReal_wp&);
    
    template<class MatrixJacob>
    void ComputeElemMatrix(const VectR3& PtsAi, const VectR3& PtsAj, 
			   bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			   const VectR3& PointsNodalI, const VectR3& PointsNodalJ,
			   const VectR3& PointsQuadI, const VectR3& PointsQuadJ,
			   const VectR3& NormaleNodalI, const VectR3& NormaleNodalJ,
			   const VectR3& NormaleQuadI, const VectR3& NormaleQuadJ,
			   const MatrixJacob& MatJacobNodalI, const MatrixJacob& MatJacobNodalJ,
			   const MatrixJacob& MatJacobI, const MatrixJacob& MatJacobJ,
			   const EuclidianDistanceClass_Base<Dimension3>& var, Matrix<Real_wp>& mat_elem,
			   Matrix<Real_wp>&, Matrix<Real_wp>&);
    
    void ComputeMassMatrix(const VectR3& PtsAi, bool curved,
			   const VectReal_wp& Dsi, Matrix<Real_wp>& mat_elem);

    template<class MatrixJacob>
    void ComputeStiffnessMatrix(const VectR3& PtsAi, bool curved,
				const VectReal_wp& Dsi, const VectR3& NormaleQuad,
				const MatrixJacob& MatJacobI,
				const Complex_wp& kinf, const Real_wp& epsilon,
				Matrix<Complex_wp>& mat_elem);
    
  };

}

#define MONTJOIE_FILE_SINGULAR_INTEGRAL_EQUATION_3D_HXX
#endif
