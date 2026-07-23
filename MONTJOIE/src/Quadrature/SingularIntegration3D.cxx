#ifndef MONTJOIE_FILE_SINGULAR_INTEGRAL_EQUATION_3D_CXX

#include "SingularIntegration3D.hxx"

namespace Montjoie
{
  
  //! (j, k) is the opposite edge of vertex i (triangle)
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::GetOppositeEdge(int i, int& j, int& k)
    {
      if (i == 0)
	{
	  j = 1; k = 2;
	}
      else if (i == 1)
	{
	  j = 2; k = 0;
	}
      else
	{
	  j = 0; k = 1;
	}
    }
    
  
  //! returns the opposite vertex of edge (i, j)
  template<class T>
  int SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::GetOppositeVertex(int i, int j)
  {
    if (i == 0)
      {
	if (j == 1)
	  return 2;
	else
	  return 1;
      }
    else if (i == 1)
      {
	if (j == 2)
	  return 0;
	else
	  return 2;
      }
    else
      {
	if (j == 0)
	  return 1;
	else
	  return 0;
      }
  }
  
  
  //! fills dof numbers located on the edge [n0, n1]
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::FillDofEdge(int n0, int n1, int order, int nb_vert, IVect& numI)
  {
    numI.Reallocate(1+order);
    numI(0) = n0; numI(order) = n1;
    if (nb_vert == 3)
      {
	if (n0 == 0)
	  {
	    if (n1 == 1)
	      for (int i = 1; i < order; i++)
		numI(i) = 2+i;
	    else
	      for (int i = 1; i < order; i++)
		numI(i) = 3*order-i;
	  }
	else if (n0 == 1)
	  {
	    if (n1 == 2)
	      for (int i = 1; i < order; i++)
		numI(i) = 1+order+i;
	    else
	      for (int i = 1; i < order; i++)
		numI(i) = 2+order-i;
	  }
	else
	  {
	    if (n1 == 0)
	      for (int i = 1; i < order; i++)
		numI(i) = 2*order+i;
	    else
	      for (int i = 1; i < order; i++)
		numI(i) = 1+2*order-i;
	  }
      }
    else
      {
	if (n0 == 0)
	  {
	    if (n1 == 1)
	      for (int i = 1; i < order; i++)
		numI(i) = 3+i;
	    else
	      for (int i = 1; i < order; i++)
		numI(i) = 4*order-i;
	  }
	else if (n0 == 1)
	  {
	    if (n1 == 2)
	      for (int i = 1; i < order; i++)
		numI(i) = 2+order+i;
	    else
	      for (int i = 1; i < order; i++)
		numI(i) = 3+order-i;
	  }
	else if (n0 == 2)
	  {
	    if (n1 == 3)
	      for (int i = 1; i < order; i++)
		numI(i) = 1+2*order+i;
	    else
	      for (int i = 1; i < order; i++)
		numI(i) = 2+2*order-i;	      
	  }
	else
	  {
	    if (n1 == 0)
	      for (int i = 1; i < order; i++)
		numI(i) = 3*order+i;
	    else
	      for (int i = 1; i < order; i++)
		numI(i) = 1+3*order-i;	      
	  }
      }    
  }
  
  
  //! Evaluates the double integral for panels which are far from each other
  template<class T> template<class MatrixJacob>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::ComputeFarPanelIntegral(int nb_i, int nb_j, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			    const VectR3& PointsQuadI, const VectR3& PointsQuadJ,
			    const VectR3& NormaleQuadI, const VectR3& NormaleQuadJ,
			    const MatrixJacob& MatJacobI, const MatrixJacob& MatJacobJ,
			    const EuclidianDistanceClass_Base<Dimension3>& var, Matrix<Real_wp>& mat_elem,
			    Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ)
  {
    Matrix<R2>* GradPhiI, *GradPhiJ;
    Matrix<Real_wp>* ValPhiI, *ValPhiJ;
    Vector<Real_wp>* WeightsI, *WeightsJ;
    if (nb_i == 3)
      {
	GradPhiI = &GradPhiTriFar;
	ValPhiI = &ValPhiTriFar;
	WeightsI = &WeightsTriFar;
      }
    else
      {
	GradPhiI = &GradPhiQuadFar;
	ValPhiI = &ValPhiQuadFar;
	WeightsI = &WeightsQuadFar;
      }

    if (nb_j == 3)
      {
	GradPhiJ = &GradPhiTriFar;
	ValPhiJ = &ValPhiTriFar;
	WeightsJ = &WeightsTriFar;
      }
    else
      {
	GradPhiJ = &GradPhiQuadFar;
	ValPhiJ = &ValPhiQuadFar;
	WeightsJ = &WeightsQuadFar;
      }
    
    IVect numI, numJ;
    ComputePanelIntegral(Dsi, Dsj, *WeightsI, *WeightsJ,
			 PointsQuadI, PointsQuadJ, NormaleQuadI, NormaleQuadJ,
			 MatJacobI, MatJacobJ, *GradPhiI, *GradPhiJ, *ValPhiI, *ValPhiJ,
			 var, mat_elem, mat_elem_diagI, mat_elem_diagJ, numI, numJ);
  }    
  

  //! Evaluates the double integral for panels which share a common vertex
  template<class T> template<class MatrixJacob>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::ComputeVertexPanelIntegral(int pt_common_i, int pt_common_j, int nb_i, int nb_j,
			       bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			       const VectR3& PointsNodalI, const VectR3& PointsNodalJ,
			       const VectR3& NormaleNodalI, const VectR3& NormaleNodalJ,
			       const MatrixJacob& MatJacobNodalI, const MatrixJacob& MatJacobNodalJ,
			       const EuclidianDistanceClass_Base<Dimension3>& var, Matrix<Real_wp>& mat_elem,
			       Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ)
  {
    if (type_integration == LENOIR_SALLES)
      {
	if ((curved) || (nb_i == 4) || (nb_j == 4))
	  {
	    cout << "Not implemented" << endl;
	    cout << "Select a numerical integration" << endl;
	    abort();
	  }
	
	// integral of |x-y|^\beta
	
	// points of the opposite edge to the common vertex in element K
	int nk1, nk2;
	GetOppositeEdge(pt_common_i, nk1, nk2);	
		
	// same stuff for element L
	int nl1, nl2;
	GetOppositeEdge(pt_common_j, nl1, nl2);	
		
	// evaluating integral h_K \int_{\partial K x L} |x-y|^\beta
	R3 vec_u = PointsNodalJ(1) - PointsNodalJ(0), vec_v = PointsNodalJ(2) - PointsNodalJ(0);
	Real_wp IntL(0);
	for (int i = 0; i < PointsTri.GetM(); i++)
	  {
	    R3 pt_x = PointsNodalJ(0) + PointsTri(i)(0)*vec_u + PointsTri(i)(1)*vec_v;
	    Real_wp poids_x = WeightsTri(i)*Dsi(0)*Dsj(0);
	    for (int j = 0; j <= lob.GetOrder(); j++)
	      {
		R3 pt_y = (1.0-lob.Points(j))*PointsNodalI(nk1) + lob.Points(j)*PointsNodalI(nk2);
		Real_wp poids_y = lob.Weights(j);
		IntL += poids_x*poids_y*pow(var.GetDistance(pt_x, pt_y), beta);
	      }
	  }
	
	// adding integral h_L \int_{\partial L x K} |x-y|^\beta
	vec_u = PointsNodalI(1) - PointsNodalI(0); vec_v = PointsNodalI(2) - PointsNodalI(0);
	for (int i = 0; i < PointsTri.GetM(); i++)
	  {
	    R3 pt_x = PointsNodalI(0) + PointsTri(i)(0)*vec_u + PointsTri(i)(1)*vec_v;
	    Real_wp poids_x = WeightsTri(i)*Dsi(0)*Dsj(0);
	    for (int j = 0; j <= lob.GetOrder(); j++)
	      {
		R3 pt_y = (1.0-lob.Points(j))*PointsNodalJ(nl1) + lob.Points(j)*PointsNodalJ(nl2);
		Real_wp poids_y = lob.Weights(j);
		IntL += poids_x*poids_y*pow(var.GetDistance(pt_x, pt_y), beta);
	      }
	  }
	
	Real_wp val = IntL / (beta + 4);

	const ElementReference<Dimension2, 1>& Fb = *fe_tri;
	int nb_dof = Fb.GetNbDof();
	if ((Fb.GetOrder() != 1) || (type_operator == OPERATOR_DIFF))
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
	
	R3 grad_phi_k, grad_phi_n;	
	switch (type_operator)
	  {
	  case OPERATOR_GRAD :
	    {
              for (int k = 0; k < nb_dof; k++)
                for (int n = 0; n < nb_dof; n++)
                  {
                    Mlt(MatJacobNodalI(0), GradPhiTriFar(k, 0), grad_phi_k);
                    Mlt(MatJacobNodalJ(0), GradPhiTriFar(n, 0), grad_phi_n);
                    mat_elem(k, n) += val*DotProd(grad_phi_k, grad_phi_n);
                  }	    
	      break;
            }
          case OPERATOR_DS :
            {
              for (int k = 0; k < nb_dof; k++)
		for (int n = 0; n < nb_dof; n++)
		  {
		    Mlt(MatJacobNodalI(0), GradPhiTriFar(k, 0), grad_phi_k);
		    Mlt(MatJacobNodalJ(0), GradPhiTriFar(n, 0), grad_phi_n);
		    TimesProd(NormaleNodalI(0), grad_phi_k, vec_u);
		    TimesProd(NormaleNodalJ(0), grad_phi_n, vec_v);
		    mat_elem(k, n) += val*DotProd(vec_u, vec_v);
		  }
            }
            break;	      
	  }
      }
    else
      {
	VectReal_wp* WeightsX_, *WeightsY_;
	Matrix<Real_wp>* NodalPhiX_, *ValPhiX_, *NodalPhiY_, *ValPhiY_;
	Matrix<R2>* GradPhiX_, *GradPhiY_;
	if (nb_i == 3)
	  {
	    WeightsX_ = &WeightsVertexTri(pt_common_i);
	    NodalPhiX_ = &NodalPhiTriVertex(pt_common_i);
	    ValPhiX_ = &ValPhiTriVertex(pt_common_i);
	    GradPhiX_ = &GradPhiTriVertex(pt_common_i);
	  }
	else
	  {
	    WeightsX_ = &WeightsVertexQuad(pt_common_i);
	    NodalPhiX_ = &NodalPhiQuadVertex(pt_common_i);
	    ValPhiX_ = &ValPhiQuadVertex(pt_common_i);
	    GradPhiX_ = &GradPhiQuadVertex(pt_common_i);
	  }
	
	if (nb_j == 3)
	  {
	    WeightsY_ = &WeightsVertexTri(pt_common_j);
	    NodalPhiY_ = &NodalPhiTriVertex(pt_common_j);
	    ValPhiY_ = &ValPhiTriVertex(pt_common_j);
	    GradPhiY_ = &GradPhiTriVertex(pt_common_j);
	  }
	else
	  {
	    WeightsY_ = &WeightsVertexQuad(pt_common_j);
	    NodalPhiY_ = &NodalPhiQuadVertex(pt_common_j);
	    ValPhiY_ = &ValPhiQuadVertex(pt_common_j);
	    GradPhiY_ = &GradPhiQuadVertex(pt_common_j);
	  }
	
	VectReal_wp& WeightsX = *WeightsX_;
	VectReal_wp& WeightsY = *WeightsY_;
	Matrix<Real_wp>& NodalPhiX = *NodalPhiX_;
	Matrix<Real_wp>& NodalPhiY = *NodalPhiY_;
	Matrix<Real_wp>& ValPhiX = *ValPhiX_;
	Matrix<Real_wp>& ValPhiY = *ValPhiY_;
	Matrix<R2>& GradPhiX = *GradPhiX_;
	Matrix<R2>& GradPhiY = *GradPhiY_;
	
	int Nx = WeightsX.GetM();
	VectReal_wp DsQuadX(Nx); DsQuadX.Fill(0);
	VectR3 PointsQuadX(Nx), NormaleQuadX(Nx);
	MatrixJacob MatJacobX(Nx); Real_wp val;
	for (int i = 0; i < NodalPhiX.GetM(); i++)
	  for (int j = 0; j < WeightsX.GetM(); j++)
	    {
	      val = NodalPhiX(i, j);
	      DsQuadX(j) += val*Dsi(i);
	      Add(val, PointsNodalI(i), PointsQuadX(j));
	      Add(val, NormaleNodalI(i), NormaleQuadX(j));
	      Add(val, MatJacobNodalI(i), MatJacobX(j));
	    }
	
	int Ny = WeightsY.GetM();
	VectReal_wp DsQuadY(Ny); DsQuadY.Fill(0);
	VectR3 PointsQuadY(Ny), NormaleQuadY(Ny);
	MatrixJacob MatJacobY(Ny);
	for (int i = 0; i < NodalPhiY.GetM(); i++)
	  for (int j = 0; j < WeightsY.GetM(); j++)
	    {
	      val = NodalPhiY(i, j);
	      DsQuadY(j) += val*Dsj(i);
	      Add(val, PointsNodalJ(i), PointsQuadY(j));
	      Add(val, NormaleNodalJ(i), NormaleQuadY(j));
	      Add(val, MatJacobNodalJ(i), MatJacobY(j));
	    }
	
	IVect numI(1), numJ(1);
	numI(0) = pt_common_i; numJ(0) = pt_common_j;
	ComputePanelIntegral(DsQuadX, DsQuadY, WeightsX, WeightsY, PointsQuadX, PointsQuadY,
			     NormaleQuadX, NormaleQuadY, MatJacobX, MatJacobY,
			     GradPhiX, GradPhiY, ValPhiX, ValPhiY,
			     var, mat_elem, mat_elem_diagI, mat_elem_diagJ, numI, numJ);
      }
  }
  

  //! Evaluates the double integral for panels which share a common edge
  template<class T> template<class MatrixJacob>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::ComputeEdgePanelIntegral(int pt_i0, int pt_i1, int pt_j0, int pt_j1, int nb_i, int nb_j,
			     bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			     const VectR3& PointsNodalI, const VectR3& PointsNodalJ,
			     const VectR3& NormaleNodalI, const VectR3& NormaleNodalJ,
			     const MatrixJacob& MatJacobNodalI, const MatrixJacob& MatJacobNodalJ,
			     const EuclidianDistanceClass_Base<Dimension3>& var, Matrix<Real_wp>& mat_elem,
			     Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ)
  {
    
    if (type_integration == LENOIR_SALLES)
      {
	if ((curved) || (nb_i == 4) || (nb_j == 4))
	  {
	    cout << "Not implemented" << endl;
	    cout << "Select a numerical integration" << endl;
	    abort();
	  }
	
	// integral of |x-y|^\beta 
	
	// opposite edges 
	int nk1, nk2;
	GetOppositeEdge(pt_i0, nk1, nk2);	
	int nl1, nl2;
	GetOppositeEdge(pt_j1, nl1, nl2);	
	
	// integral with these two edges
	Real_wp IntL(0);
	for (int i = 0; i <= lob.GetOrder(); i++)
	  {
	    R3 pt_x = (1.0-lob.Points(i))*PointsNodalI(nk1) + lob.Points(i)*PointsNodalI(nk2);
	    Real_wp poids_x = lob.Weights(i)*Dsi(0)*Dsj(0);
	    for (int j = 0; j <= lob.GetOrder(); j++)
	      {
		R3 pt_y = (1.0-lob.Points(j))*PointsNodalJ(nl1) + lob.Points(j)*PointsNodalJ(nl2);
		Real_wp poids_y = lob.Weights(j);
		IntL += poids_x*poids_y*pow(var.GetDistance(pt_x, pt_y), beta);
	      }
	  }
	
	// same stuff with the other edges
	GetOppositeEdge(pt_i1, nk1, nk2);	
	GetOppositeEdge(pt_j0, nl1, nl2);	
	
	// integral with these two edges
	for (int i = 0; i <= lob.GetOrder(); i++)
	  {
	    R3 pt_x = (1.0-lob.Points(i))*PointsNodalI(nk1) + lob.Points(i)*PointsNodalI(nk2);
	    Real_wp poids_x = lob.Weights(i)*Dsi(0)*Dsj(0);
	    for (int j = 0; j <= lob.GetOrder(); j++)
	      {
		R3 pt_y = (1.0-lob.Points(j))*PointsNodalJ(nl1) + lob.Points(j)*PointsNodalJ(nl2);
		Real_wp poids_y = lob.Weights(j);
		IntL += poids_x*poids_y*pow(var.GetDistance(pt_x, pt_y), beta);
	      }
	  }
	
	// integral with triangles
	nk1 = GetOppositeVertex(pt_i0, pt_i1);
	R3 vec_u = PointsNodalJ(1) - PointsNodalJ(0), vec_v = PointsNodalJ(2) - PointsNodalJ(0);
	R3 pt_y = PointsNodalI(nk1);
	for (int i = 0; i < PointsTri.GetM(); i++)
	  {
	    R3 pt_x = PointsNodalJ(0) + PointsTri(i)(0)*vec_u + PointsTri(i)(1)*vec_v;
	    Real_wp poids_x = WeightsTri(i)*Dsi(0)*Dsj(0);
	    IntL += poids_x*pow(var.GetDistance(pt_x, pt_y), beta);
	  }
	
	nl1 = GetOppositeVertex(pt_j0, pt_j1);
	vec_u = PointsNodalI(1) - PointsNodalI(0); vec_v = PointsNodalI(2) - PointsNodalI(0);
	pt_y = PointsNodalJ(nl1);
	for (int i = 0; i < PointsTri.GetM(); i++)
	  {
	    R3 pt_x = PointsNodalI(0) + PointsTri(i)(0)*vec_u + PointsTri(i)(1)*vec_v;
	    Real_wp poids_x = WeightsTri(i)*Dsi(0)*Dsj(0);
	    IntL += poids_x*pow(var.GetDistance(pt_x, pt_y), beta);
	  }
	
	Real_wp val = IntL / ((beta+4)*(beta+3));

	const ElementReference<Dimension2, 1>& Fb = *fe_tri;
	int nb_dof = Fb.GetNbDof();
	if ((Fb.GetOrder() != 1) || (type_operator == OPERATOR_DIFF))
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
	
	R3 grad_phi_k, grad_phi_n;	
	switch (type_operator)
	  {
	  case OPERATOR_GRAD :
	    {
              for (int k = 0; k < nb_dof; k++)
                for (int n = 0; n < nb_dof; n++)
                  {
                    Mlt(MatJacobNodalI(0), GradPhiTriFar(k, 0), grad_phi_k);
                    Mlt(MatJacobNodalJ(0), GradPhiTriFar(n, 0), grad_phi_n);
                    mat_elem(k, n) += val*DotProd(grad_phi_k, grad_phi_n);
                  }	
            }    
            break;
          case OPERATOR_DS :
            {
              for (int k = 0; k < nb_dof; k++)
		for (int n = 0; n < nb_dof; n++)
		  {
		    Mlt(MatJacobNodalI(0), GradPhiTriFar(k, 0), grad_phi_k);
		    Mlt(MatJacobNodalJ(0), GradPhiTriFar(n, 0), grad_phi_n);
		    TimesProd(NormaleNodalI(0), grad_phi_k, vec_u);
		    TimesProd(NormaleNodalJ(0), grad_phi_n, vec_v);
		    mat_elem(k, n) += val*DotProd(vec_u, vec_v);
		  }
            }
            break;	      
	  }
      }
    else
      {
	int ni = GetOppositeVertex(pt_i0, pt_i1);
	int nj = GetOppositeVertex(pt_j0, pt_j1);
	VectReal_wp* WeightsX_, *WeightsY_;
	Matrix<Real_wp>* NodalPhiX_, *ValPhiX_, *NodalPhiY_, *ValPhiY_;
	Matrix<R2>* GradPhiX_, *GradPhiY_;
	IVect numI, numJ;
	
	if (nb_i == 3)
	  {
	    WeightsX_ = &WeightsEdgeTri(ni);
	    NodalPhiX_ = &NodalPhiTriEdge(ni);
	    ValPhiX_ = &ValPhiTriEdge(ni);
	    GradPhiX_ = &GradPhiTriEdge(ni);
	    FillDofEdge(pt_i0, pt_i1, fe_tri->GetOrder(), nb_i, numI);
	  }
	else
	  {
	    WeightsX_ = &WeightsEdgeQuad(ni);
	    NodalPhiX_ = &NodalPhiQuadEdge(ni);
	    ValPhiX_ = &ValPhiQuadEdge(ni);
	    GradPhiX_ = &GradPhiQuadEdge(ni);
	    FillDofEdge(pt_i0, pt_i1, fe_quad->GetOrder(), nb_i, numI);
	  }
	
	if (nb_j == 3)
	  {
	    WeightsY_ = &WeightsEdgeTri(nj);
	    NodalPhiY_ = &NodalPhiTriEdge(nj);
	    ValPhiY_ = &ValPhiTriEdge(nj);
	    GradPhiY_ = &GradPhiTriEdge(nj);
	    FillDofEdge(pt_j0, pt_j1, fe_tri->GetOrder(), nb_j, numJ);
	  }
	else
	  {
	    WeightsY_ = &WeightsEdgeQuad(nj);
	    NodalPhiY_ = &NodalPhiQuadEdge(nj);
	    ValPhiY_ = &ValPhiQuadEdge(nj);
	    GradPhiY_ = &GradPhiQuadEdge(nj);
	    FillDofEdge(pt_j0, pt_j1, fe_quad->GetOrder(), nb_j, numJ);
	  }
	
	VectReal_wp& WeightsX = *WeightsX_;
	VectReal_wp& WeightsY = *WeightsY_;
	Matrix<Real_wp>& NodalPhiX = *NodalPhiX_;
	Matrix<Real_wp>& NodalPhiY = *NodalPhiY_;
	Matrix<Real_wp>& ValPhiX = *ValPhiX_;
	Matrix<Real_wp>& ValPhiY = *ValPhiY_;
	Matrix<R2>& GradPhiX = *GradPhiX_;
	Matrix<R2>& GradPhiY = *GradPhiY_;
	
	int Nx = WeightsX.GetM();
	VectReal_wp DsQuadX(Nx); DsQuadX.Fill(0);
	VectR3 PointsQuadX(Nx), NormaleQuadX(Nx);
	MatrixJacob MatJacobX(Nx); Real_wp val;
	for (int i = 0; i < NodalPhiX.GetM(); i++)
	  for (int j = 0; j < WeightsX.GetM(); j++)
	    {
	      val = NodalPhiX(i, j);
	      DsQuadX(j) += val*Dsi(i);
	      Add(val, PointsNodalI(i), PointsQuadX(j));
	      Add(val, NormaleNodalI(i), NormaleQuadX(j));
	      Add(val, MatJacobNodalI(i), MatJacobX(j));
	    }
	
	int Ny = WeightsY.GetM();
	VectReal_wp DsQuadY(Ny); DsQuadY.Fill(0);
	VectR3 PointsQuadY(Ny), NormaleQuadY(Ny);
	MatrixJacob MatJacobY(Ny);
	for (int i = 0; i < NodalPhiY.GetM(); i++)
	  for (int j = 0; j < WeightsY.GetM(); j++)
	    {
	      val = NodalPhiY(i, j);
	      DsQuadY(j) += val*Dsj(i);
	      Add(val, PointsNodalJ(i), PointsQuadY(j));
	      Add(val, NormaleNodalJ(i), NormaleQuadY(j));
	      Add(val, MatJacobNodalJ(i), MatJacobY(j));
	    }

	ComputePanelIntegral(DsQuadX, DsQuadY, WeightsX, WeightsY, PointsQuadX, PointsQuadY,
			     NormaleQuadX, NormaleQuadY, MatJacobX, MatJacobY,
			     GradPhiX, GradPhiY, ValPhiX, ValPhiY,
			     var, mat_elem, mat_elem_diagI, mat_elem_diagJ, numI, numJ);
      }
  }
  

  //! Evaluates the double integral for identical panels
  template<class T> template<class MatrixJacob>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::ComputeSamePanelIntegral(int nb_i, bool curved, const VectReal_wp& Dsi,
			     const VectR3& PointsNodalI, const VectR3& NormaleNodalI,
			     const MatrixJacob& MatJacobNodalI, const EuclidianDistanceClass_Base<Dimension3>& var,
			     Matrix<Real_wp>& mat_elem,
			     Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ)
  {
    if (type_integration == LENOIR_SALLES)
      {
	if ((curved) || (nb_i == 4))
	  {
	    cout << "Not implemented" << endl;
	    cout << "Select a numerical integration" << endl;
	    abort();
	  }
	
	// integral of |x-y|^\beta 
	// three integrals
	Real_wp IntL = 0;
	R3 pt_y = PointsNodalI(2);
	for (int i = 0; i <= lob.GetOrder(); i++)
	  {
	    R3 pt_x = (1.0-lob.Points(i))*PointsNodalI(0) + lob.Points(i)*PointsNodalI(1);
	    Real_wp poids_x = lob.Weights(i)*Dsi(0)*Dsi(0);
	    IntL += poids_x*pow(var.GetDistance(pt_x, pt_y), beta);
	  }

	pt_y = PointsNodalI(1);
	for (int i = 0; i <= lob.GetOrder(); i++)
	  {
	    R3 pt_x = (1.0-lob.Points(i))*PointsNodalI(0) + lob.Points(i)*PointsNodalI(2);
	    Real_wp poids_x = lob.Weights(i)*Dsi(0)*Dsi(0);
	    IntL += poids_x*pow(var.GetDistance(pt_x, pt_y), beta);
	  }

	pt_y = PointsNodalI(0);
	for (int i = 0; i <= lob.GetOrder(); i++)
	  {
	    R3 pt_x = (1.0-lob.Points(i))*PointsNodalI(1) + lob.Points(i)*PointsNodalI(2);
	    Real_wp poids_x = lob.Weights(i)*Dsi(0)*Dsi(0);
	    IntL += poids_x*pow(var.GetDistance(pt_x, pt_y), beta);
	  }
	
	Real_wp val = 2.0*IntL / ((beta+4)*(beta+3)*(beta+2));
	
	const ElementReference<Dimension2, 1>& Fb = *fe_tri;
	int nb_dof = Fb.GetNbDof();
	if ((Fb.GetOrder() != 1) || (type_operator == OPERATOR_DIFF))
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
	
	R3 grad_phi_k, grad_phi_n, vec_u, vec_v;	
	switch (type_operator)
	  {
	  case OPERATOR_GRAD :
	    {
              for (int k = 0; k < nb_dof; k++)
                for (int n = 0; n < nb_dof; n++)
                  {
                    Mlt(MatJacobNodalI(0), GradPhiTriFar(k, 0), grad_phi_k);
                    Mlt(MatJacobNodalI(0), GradPhiTriFar(n, 0), grad_phi_n);
                    mat_elem(k, n) += val*DotProd(grad_phi_k, grad_phi_n);
                  }	
            }    
            break;
          case OPERATOR_DS :
            {
              for (int k = 0; k < nb_dof; k++)
		for (int n = 0; n < nb_dof; n++)
		  {
		    Mlt(MatJacobNodalI(0), GradPhiTriFar(k, 0), grad_phi_k);
		    Mlt(MatJacobNodalI(0), GradPhiTriFar(n, 0), grad_phi_n);
		    TimesProd(NormaleNodalI(0), grad_phi_k, vec_u);
		    TimesProd(NormaleNodalI(0), grad_phi_n, vec_v);
		    mat_elem(k, n) += val*DotProd(vec_u, vec_v);
		  }
            }
            break;	      
	  }
      }
    else
      {
	VectReal_wp* WeightsX_;
	Matrix<Real_wp>* NodalPhiX_, *ValPhiX_;
	Matrix<R2>* GradPhiX_;
	
	Vector<VectReal_wp>* WeightsY_;
	Vector<Matrix<Real_wp> >* NodalPhiY_, *ValPhiY_;
	Vector<Matrix<R2> >* GradPhiY_;
	
	if (nb_i == 3)
	  {
	    WeightsX_ = &WeightsCommonExtTri;
	    NodalPhiX_ = &NodalPhiTriCommonExt;
	    ValPhiX_ = &ValPhiTriCommonExt;
	    GradPhiX_ = &GradPhiTriCommonExt;
	    
	    WeightsY_ = &WeightsCommonTri;
	    NodalPhiY_ = &NodalPhiTriCommon;
	    ValPhiY_ = &ValPhiTriCommon;
	    GradPhiY_ = &GradPhiTriCommon;
	  }
	else
	  {
	    WeightsX_ = &WeightsCommonExtQuad;
	    NodalPhiX_ = &NodalPhiQuadCommonExt;
	    ValPhiX_ = &ValPhiQuadCommonExt;
	    GradPhiX_ = &GradPhiQuadCommonExt;

	    WeightsY_ = &WeightsCommonQuad;
	    NodalPhiY_ = &NodalPhiQuadCommon;
	    ValPhiY_ = &ValPhiQuadCommon;
	    GradPhiY_ = &GradPhiQuadCommon;
	  }
	
	VectReal_wp& WeightsX = *WeightsX_;
	Vector<VectReal_wp>& WeightsY = *WeightsY_;
	Matrix<Real_wp>& NodalPhiX = *NodalPhiX_;
	Vector<Matrix<Real_wp> >& NodalPhiY = *NodalPhiY_;
	Matrix<Real_wp>& ValPhiX = *ValPhiX_;
	Vector<Matrix<Real_wp> >& ValPhiY = *ValPhiY_;
	Matrix<R2>& GradPhiX = *GradPhiX_;
	Vector<Matrix<R2> >& GradPhiY = *GradPhiY_;
	
	int Nx = WeightsX.GetM();
	VectReal_wp DsQuadX(Nx); DsQuadX.Fill(0);
	VectR3 PointsQuadX(Nx), NormaleQuadX(Nx);
	MatrixJacob MatJacobX(Nx); Real_wp val;
	for (int i = 0; i < NodalPhiX.GetM(); i++)
	  for (int j = 0; j < WeightsX.GetM(); j++)
	    {
	      val = NodalPhiX(i, j);
	      DsQuadX(j) += val*Dsi(i);
	      Add(val, PointsNodalI(i), PointsQuadX(j));
	      Add(val, NormaleNodalI(i), NormaleQuadX(j));
	      Add(val, MatJacobNodalI(i), MatJacobX(j));
	    }

	Vector<VectReal_wp> DsQuadY(Nx);
	Vector<VectR3> PointsQuadY(Nx), NormaleQuadY(Nx);
	Vector<MatrixJacob> MatJacobY(Nx);
	for (int i0 = 0; i0 < Nx; i0++)
	  {
	    int Ny = WeightsY(i0).GetM();
	    DsQuadY(i0).Reallocate(Ny); DsQuadY(i0).Fill(0);
	    PointsQuadY(i0).Reallocate(Ny);
	    NormaleQuadY(i0).Reallocate(Ny);
	    MatJacobY(i0).Reallocate(Ny);
	    for (int i = 0; i < NodalPhiY(i0).GetM(); i++)
	      for (int j = 0; j < WeightsY(i0).GetM(); j++)
		{
		  val = NodalPhiY(i0)(i, j);
		  DsQuadY(i0)(j) += val*Dsi(i);
		  Add(val, PointsNodalI(i), PointsQuadY(i0)(j));
		  Add(val, NormaleNodalI(i), NormaleQuadY(i0)(j));
		  Add(val, MatJacobNodalI(i), MatJacobY(i0)(j));
		}	    
	  }
	
	R3 grad_phi_k, grad_phi_n, vec_u, vec_v;
	for (int i = 0; i < PointsQuadX.GetM(); i++)
	  for (int j = 0; j < PointsQuadY(i).GetM(); j++)
	    {
	      val = WeightsX(i)*WeightsY(i)(j);
	      val *= pow(var.GetDistance(PointsQuadX(i), PointsQuadY(i)(j)), beta);
	      val *= DsQuadX(i)*DsQuadY(i)(j);
	      switch (type_operator)
		{
		case OPERATOR_GRAD :
		  for (int k = 0; k < GradPhiX.GetM(); k++)
		    for (int n = 0; n < GradPhiY(i).GetM(); n++)
		      {
			Mlt(MatJacobX(i), GradPhiX(k, i), grad_phi_k);
			Mlt(MatJacobY(i)(j), GradPhiY(i)(n, j), grad_phi_n);
			mat_elem(k, n) += val*DotProd(grad_phi_k, grad_phi_n);
		      }	      
		  break;
		case OPERATOR_DS :
		  for (int k = 0; k < GradPhiX.GetM(); k++)
		    for (int n = 0; n < GradPhiY(i).GetM(); n++)
		      {
			Mlt(MatJacobX(i), GradPhiX(k, i), grad_phi_k);
			Mlt(MatJacobY(i)(j), GradPhiY(i)(n, j), grad_phi_n);
			TimesProd(NormaleQuadX(i), grad_phi_k, vec_u);
			TimesProd(NormaleQuadY(i)(j), grad_phi_n, vec_v);
			mat_elem(k, n) += val*DotProd(vec_u, vec_v);
		      }
		  break;
		case OPERATOR_DIFF :
		  for (int k = 0; k < ValPhiX.GetM(); k++)
		    for (int n = 0; n < ValPhiY(i).GetM(); n++)
		      {
			mat_elem(k, n) += val*(ValPhiX(k, i)-ValPhiY(i)(k, j))
			  *(ValPhiX(n, i) - ValPhiY(i)(n, j));
		      }
		}
	    }
      }
  }
  
  
  //! internal function
  template<class T> template<class MatrixJacob>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::ComputePanelIntegral(const VectReal_wp& DsQuadX, const VectReal_wp& DsQuadY,
			 const VectReal_wp& WeightsX, const VectReal_wp& WeightsY,
			 const VectR3& PointsQuadX, const VectR3& PointsQuadY,
			 const VectR3& NormaleQuadX, const VectR3& NormaleQuadY,
			 const MatrixJacob& MatJacobX, const MatrixJacob& MatJacobY,
			 const Matrix<R2>& GradPhiX, const Matrix<R2>& GradPhiY,
			 const Matrix<Real_wp>& ValPhiX, const Matrix<Real_wp>& ValPhiY,
			 const EuclidianDistanceClass_Base<Dimension3>& var, Matrix<Real_wp>& mat_elem,
			 Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ,
			 const IVect& numI, const IVect& numJ)
  {
    Real_wp val; R3 grad_phi_k, grad_phi_n, vec_u, vec_v;
    for (int i = 0; i < PointsQuadX.GetM(); i++)
      for (int j = 0; j < PointsQuadY.GetM(); j++)
	{
	  val = WeightsX(i)*WeightsY(j);
	  val *= pow(var.GetDistance(PointsQuadX(i), PointsQuadY(j)), beta);
	  //DISP(i); DISP(j); DISP(PointsQuadX(i)); DISP(PointsQuadY(j));
	  val *= DsQuadX(i)*DsQuadY(j); //DISP(val);
	  switch (type_operator)
	    {
	    case OPERATOR_GRAD :
	      for (int k = 0; k < GradPhiX.GetM(); k++)
		for (int n = 0; n < GradPhiY.GetM(); n++)
		  {
		    Mlt(MatJacobX(i), GradPhiX(k, i), grad_phi_k);
		    Mlt(MatJacobY(j), GradPhiY(n, j), grad_phi_n);
		    mat_elem(k, n) += val*DotProd(grad_phi_k, grad_phi_n);
		  }	      
	      break;
	    case OPERATOR_DS :
	      for (int k = 0; k < GradPhiX.GetM(); k++)
		for (int n = 0; n < GradPhiY.GetM(); n++)
		  {
		    Mlt(MatJacobX(i), GradPhiX(k, i), grad_phi_k);
		    Mlt(MatJacobY(j), GradPhiY(n, j), grad_phi_n);
		    TimesProd(NormaleQuadX(i), grad_phi_k, vec_u);
		    TimesProd(NormaleQuadY(j), grad_phi_n, vec_v);
		    mat_elem(k, n) += val*DotProd(vec_u, vec_v);
		  }
	      break;
	    case OPERATOR_DIFF :
	      for (int k = 0; k < ValPhiX.GetM(); k++)
		for (int n = 0; n < ValPhiY.GetM(); n++)
		  {
		    //DISP(k); DISP(n); DISP(ValPhiX(k, i)); DISP(ValPhiY(n, j));
		    mat_elem(k, n) -= val*ValPhiX(k, i)*ValPhiY(n, j);
		    mat_elem_diagI(k, n) += val*ValPhiX(k, i)*ValPhiX(n, i);
		    mat_elem_diagJ(k, n) += val*ValPhiY(k, j)*ValPhiY(n, j);
		  }
	      
	      //DISP(numI); DISP(numJ);
	      for (int k = 0; k < numI.GetM(); k++)
		{
		  int n = k;
		  //DISP(numI(k)); DISP(numJ(n)); DISP(ValPhiX(numI(n), i)); DISP(ValPhiY(numJ(k), j));
		  mat_elem(numI(k), numJ(n)) -= val*ValPhiX(numI(n), i)*ValPhiY(numJ(k), j);
		}
	    }
	}
    
    //DISP(mat_elem(1, 2));
  }

  
  //! default constructor
  template<class T>
  SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::SingularDoubleQuadratureGalerkin_Base()
  {
    fe_tri = NULL; fe_quad = NULL;
    param_eta_quadrature = 0.0;
    type_operator = OPERATOR_DS;
    beta = 0.5;
    type_integration = LENOIR_SALLES;
  }
  
  
  //! sets the operator to be computed with the exponent b
  /*!
    the operator is either :
    - \int \int \grad phi(x) \cdot \grad psi(y) |x-y|^b dx dy
    - \int \int n \times \grad phi(x) \cdot n \times \grad psi(y) |x-y|^b dx dy
    - \int \int (phi(x)-phi(y)) (psi(x) - psi(y)) |x-y|^b dx dy
   */
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::SetOperator(int type, const Real_wp& b)
  {
    type_operator = type;
    beta = b;
  }
  
  
  //! uses analytical integration whenever possible
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::SetAnalyticalIntegration(const Real_wp& eta, int rp, int rj)
  {    
    type_integration = LENOIR_SALLES;
    param_eta_quadrature = eta;

    lob.ConstructQuadrature(rj);

    if (2*rj <= 20)
      TriangleQuadrature::ConstructQuadrature(2*rj, PointsTri, WeightsTri,
					      TriangleQuadrature::QUADRATURE_GAUSS);
    else
      TriangleQuadrature::ConstructQuadrature(2*rj, PointsTri, WeightsTri,
					      TriangleQuadrature::QUADRATURE_TENSOR);      
    
  }
  
  
  //! uses numerical integration
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::SetNumericalIntegration(const Real_wp& eta, int rp, int rj, int rj_ext, int n_int, int n_ext)
  {
    type_integration = DUFFY;
    param_eta_quadrature = eta;
    
    // integration rules for elements sharing a vertex
    TriangleQuadrature::ConstructQuadrature(2*rj, PointsTri, WeightsTri,
					    TriangleQuadrature::QUADRATURE_TENSOR);
    
    WeightsVertexTri(0) = WeightsTri;
    PointsVertexTri(0).Reallocate(PointsTri.GetM());
    WeightsVertexTri(1) = WeightsTri;
    PointsVertexTri(1).Reallocate(PointsTri.GetM());
    WeightsVertexTri(2) = WeightsTri;
    PointsVertexTri(2) = PointsTri;
    Real_wp one(1);
    for (int i = 0; i < PointsTri.GetM(); i++)
      {
	Real_wp x = PointsTri(i)(0), y = PointsTri(i)(1);
	PointsVertexTri(0)(i).Init(one-x-y, x);
	PointsVertexTri(1)(i).Init(y, x);
      }
    
    int N = PointsTri.GetM();
    for (int k = 0; k < 4; k++)
      {
	WeightsVertexQuad(k).Reallocate(2*N);
	PointsVertexQuad(k).Reallocate(2*N);
      }
    
    for (int i = 0; i < N; i++)
      {
	Real_wp x = PointsTri(i)(0), y = PointsTri(i)(1);
	WeightsVertexQuad(0)(i) = WeightsTri(i);
	PointsVertexQuad(0)(i).Init(one-y, x);
	WeightsVertexQuad(0)(N+i) = WeightsTri(i);
	PointsVertexQuad(0)(N+i).Init(x, one-y);

	WeightsVertexQuad(1)(i) = WeightsTri(i);
	PointsVertexQuad(1)(i).Init(y, x);
	WeightsVertexQuad(1)(N+i) = WeightsTri(i);
	PointsVertexQuad(1)(N+i).Init(one-x, one-y);
	
	WeightsVertexQuad(2)(i) = WeightsTri(i);
	PointsVertexQuad(2)(i).Init(x+y, y);
	WeightsVertexQuad(2)(N+i) = WeightsTri(i);
	PointsVertexQuad(2)(N+i).Init(y, x+y);

	WeightsVertexQuad(3)(i) = WeightsTri(i);
	PointsVertexQuad(3)(i).Init(x, y);
	WeightsVertexQuad(3)(N+i) = WeightsTri(i);
	PointsVertexQuad(3)(N+i).Init(one-y, one-x);
      }
    
    // integration rules for elements sharing an edge
    TriangleQuadrature::ConstructQuadrature(2*rj_ext, PointsTri, WeightsTri,
					    TriangleQuadrature::QUADRATURE_GAUSS);

    for (int k = 0; k < 3; k++)
      {
	WeightsEdgeTri(k) = WeightsTri;
	PointsEdgeTri(k) = PointsTri;
      }
    
    QuadrangleQuadrature::ConstructQuadrature(2*rj_ext, PointsTri, WeightsTri);
    
    for (int k = 0; k < 4; k++)
      {
	WeightsEdgeQuad(k) = WeightsTri;
	PointsEdgeQuad(k) = PointsTri;
      }
    
    // integration rules for same elements
    TriangleQuadrature::ConstructQuadrature(2*n_ext, PointsTri, WeightsTri,
					    TriangleQuadrature::QUADRATURE_GAUSS);
    
    WeightsCommonExtTri = WeightsTri;
    PointsCommonExtTri = PointsTri;

    TriangleQuadrature::ConstructQuadrature(2*n_int, PointsTri, WeightsTri,
					    TriangleQuadrature::QUADRATURE_TENSOR);
    
    N = PointsTri.GetM();
    int Ni = PointsCommonExtTri.GetM();
    WeightsCommonTri.Reallocate(Ni);
    PointsCommonTri.Reallocate(Ni);
    for (int i = 0; i < Ni; i++)
      {
	Real_wp xi = PointsCommonExtTri(i)(0), yi = PointsCommonExtTri(i)(1);
	WeightsCommonTri(i).Reallocate(3*N);
	PointsCommonTri(i).Reallocate(3*N);
	for (int j = 0; j < N; j++)
	  {
	    Real_wp x = PointsTri(j)(0), y = PointsTri(j)(1);
	    WeightsCommonTri(i)(j) = WeightsTri(j)*yi;
	    PointsCommonTri(i)(j).Init(x + xi*y, yi*y);
	    
	    WeightsCommonTri(i)(N+j) = WeightsTri(j)*xi;
	    PointsCommonTri(i)(N+j).Init(xi*y, x + yi*y);

	    WeightsCommonTri(i)(2*N+j) = WeightsTri(j)*(one-xi-yi);
	    PointsCommonTri(i)(2*N+j).Init(one-x-y + xi*y, x + yi*y);
	  }
      }
    
    QuadrangleQuadrature::ConstructQuadrature(2*n_ext, PointsCommonExtQuad, WeightsCommonExtQuad);
    
    Ni = PointsCommonExtQuad.GetM();
    WeightsCommonQuad.Reallocate(Ni);
    PointsCommonQuad.Reallocate(Ni);
    for (int i = 0; i < Ni; i++)
      {
	Real_wp xi = PointsCommonExtQuad(i)(0), yi = PointsCommonExtQuad(i)(1);
	WeightsCommonQuad(i).Reallocate(4*N);
	PointsCommonQuad(i).Reallocate(4*N);
	for (int j = 0; j < N; j++)
	  {
	    Real_wp x = PointsTri(j)(0), y = PointsTri(j)(1);
	    WeightsCommonQuad(i)(j) = WeightsTri(j)*yi;
	    PointsCommonQuad(i)(j).Init(x + xi*y, yi*y);
	    
	    WeightsCommonQuad(i)(N+j) = WeightsTri(j)*xi;
	    PointsCommonQuad(i)(N+j).Init(xi*y, x + yi*y);
	    
	    WeightsCommonQuad(i)(2*N+j) = WeightsTri(j)*(one-xi);
	    PointsCommonQuad(i)(2*N+j).Init(one-y+y*xi, one-x-y+yi*y);

	    WeightsCommonQuad(i)(3*N+j) = WeightsTri(j)*(one-yi);
	    PointsCommonQuad(i)(3*N+j).Init(one-x-y+y*xi, one-y+yi*y);
	  }
      }
    
    for (int k = 0; k < 4; k++)
      {
	string name = "VertexQuad" + to_str(k) + ".dat";
	ofstream file_out(name.data()); file_out.precision(15);
	for (int i = 0; i < PointsCommonQuad(k).GetM(); i++)
	  {
	    PrintNoBrackets(file_out, PointsCommonQuad(k)(i));
	    file_out << '\n';
	  }
	
	file_out.close();
      }

    for (int k = 0; k < 3; k++)
      {
	string name = "VertexTri" + to_str(k) + ".dat";
	ofstream file_out(name.data()); file_out.precision(15);
	for (int i = 0; i < PointsCommonTri(k).GetM(); i++)
	  {
	    PrintNoBrackets(file_out, PointsCommonTri(k)(i));
	    file_out << '\n';
	  }
	
	file_out.close();
      }
  }
   
  
  //! sets the finite element and pre-computes arrays needed to compute integrals
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::SetTriangularFiniteElement(const ElementReference<Dimension2, 1>& elt,
			       const VectR2& points, const VectReal_wp& weights)
  {
    fe_tri = const_cast<ElementReference<Dimension2, 1>* >(&elt);
    
    int nb_nodes = elt.GetNbPointsNodalElt();
    int nb_dof = elt.GetNbDof();
    int nb_quad = points.GetM();
    
    // integration for elements far from each other
    WeightsTriFar = weights;
    ValPhiTriFar.Reallocate(nb_dof, nb_quad);
    GradPhiTriFar.Reallocate(nb_dof, nb_quad);
    VectReal_wp phi; VectR2 grad_phi;
    for (int j = 0; j < nb_quad; j++)
      {
	elt.ComputeValuesPhiRef(points(j), phi);
	elt.ComputeGradientPhiRef(points(j), grad_phi);
	for (int i = 0; i < nb_dof; i++)
	  {
	    ValPhiTriFar(i, j) = phi(i);
	    GradPhiTriFar(i, j) = grad_phi(i);
	  }
      }
    
    // integration for elements sharing a vertex
    if (type_integration == DUFFY)
      {
	nb_quad = WeightsVertexTri(0).GetM();
	for (int k = 0; k < 3; k++)
	  {
	    NodalPhiTriVertex(k).Reallocate(nb_nodes, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiNodalRef(PointsVertexTri(k)(j), phi);
		for (int i = 0; i < nb_nodes; i++)
		  NodalPhiTriVertex(k)(i, j) = phi(i);
	      }
	    
	    ValPhiTriVertex(k).Reallocate(nb_dof, nb_quad);
	    GradPhiTriVertex(k).Reallocate(nb_dof, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiRef(PointsVertexTri(k)(j), phi);
		elt.ComputeGradientPhiRef(PointsVertexTri(k)(j), grad_phi);
		for (int i = 0; i < nb_dof; i++)
		  {
		    ValPhiTriVertex(k)(i, j) = phi(i);
		    GradPhiTriVertex(k)(i, j) = grad_phi(i);
		  }
	      }
	  }	
      }

    // integration for elements sharing an edge
    if (type_integration == DUFFY)
      {
	nb_quad = WeightsEdgeTri(0).GetM();
	for (int k = 0; k < 3; k++)
	  {
	    NodalPhiTriEdge(k).Reallocate(nb_nodes, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiNodalRef(PointsEdgeTri(k)(j), phi);
		for (int i = 0; i < nb_nodes; i++)
		  NodalPhiTriEdge(k)(i, j) = phi(i);
	      }
	    
	    ValPhiTriEdge(k).Reallocate(nb_dof, nb_quad);
	    GradPhiTriEdge(k).Reallocate(nb_dof, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiRef(PointsEdgeTri(k)(j), phi);
		elt.ComputeGradientPhiRef(PointsEdgeTri(k)(j), grad_phi);
		for (int i = 0; i < nb_dof; i++)
		  {
		    ValPhiTriEdge(k)(i, j) = phi(i);
		    GradPhiTriEdge(k)(i, j) = grad_phi(i);
		  }
	      }
	  }	
      }
    
    // integration for the same triangles
    if (type_integration == DUFFY)
      {
	nb_quad = WeightsCommonExtTri.GetM();
	NodalPhiTriCommonExt.Reallocate(nb_nodes, nb_quad);
	for (int j = 0; j < nb_quad; j++)
	  {
	    elt.ComputeValuesPhiNodalRef(PointsCommonExtTri(j), phi);
	    for (int i = 0; i < nb_nodes; i++)
	      NodalPhiTriCommonExt(i, j) = phi(i);
	  }
	
	ValPhiTriCommonExt.Reallocate(nb_dof, nb_quad);
	GradPhiTriCommonExt.Reallocate(nb_dof, nb_quad);
	for (int j = 0; j < nb_quad; j++)
	  {
	    elt.ComputeValuesPhiRef(PointsCommonExtTri(j), phi);
	    elt.ComputeGradientPhiRef(PointsCommonExtTri(j), grad_phi);
	    for (int i = 0; i < nb_dof; i++)
	      {
		ValPhiTriCommonExt(i, j) = phi(i);
		GradPhiTriCommonExt(i, j) = grad_phi(i);
	      }
	  }
	
	NodalPhiTriCommon.Reallocate(WeightsCommonTri.GetM());
	ValPhiTriCommon.Reallocate(WeightsCommonTri.GetM());
	GradPhiTriCommon.Reallocate(WeightsCommonTri.GetM());
	for (int k = 0; k < WeightsCommonTri.GetM(); k++)
	  {
	    nb_quad = WeightsCommonTri(k).GetM();
	    NodalPhiTriCommon(k).Reallocate(nb_nodes, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiNodalRef(PointsCommonTri(k)(j), phi);
		for (int i = 0; i < nb_nodes; i++)
		  NodalPhiTriCommon(k)(i, j) = phi(i);
	      }
	    
	    ValPhiTriCommon(k).Reallocate(nb_dof, nb_quad);
	    GradPhiTriCommon(k).Reallocate(nb_dof, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiRef(PointsCommonTri(k)(j), phi);
		elt.ComputeGradientPhiRef(PointsCommonTri(k)(j), grad_phi);
		for (int i = 0; i < nb_dof; i++)
		  {
		    ValPhiTriCommon(k)(i, j) = phi(i);
		    GradPhiTriCommon(k)(i, j) = grad_phi(i);
		  }
	      }
	  }	
      }
  }


  //! sets the finite element and pre-computes arrays needed to compute integrals
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::SetQuadrilateralFiniteElement(const ElementReference<Dimension2, 1>& elt,
				  const VectR2& points, const VectReal_wp& weights)
  {
    fe_quad = const_cast<ElementReference<Dimension2, 1>* >(&elt);
    
    int nb_nodes = elt.GetNbPointsNodalElt();
    int nb_dof = elt.GetNbDof();
    int nb_quad = points.GetM();
    
    WeightsQuadFar = weights;
    ValPhiQuadFar.Reallocate(nb_dof, nb_quad);
    GradPhiQuadFar.Reallocate(nb_dof, nb_quad);
    VectReal_wp phi; VectR2 grad_phi;
    for (int j = 0; j < nb_quad; j++)
      {
	elt.ComputeValuesPhiRef(points(j), phi);
	elt.ComputeGradientPhiRef(points(j), grad_phi);
	for (int i = 0; i < nb_dof; i++)
	  {
	    ValPhiQuadFar(i, j) = phi(i);
	    GradPhiQuadFar(i, j) = grad_phi(i);
	  }
      }

    // integration for elements sharing a vertex
    if (type_integration == DUFFY)
      {
	nb_quad = WeightsVertexQuad(0).GetM();
	for (int k = 0; k < 4; k++)
	  {
	    NodalPhiQuadVertex(k).Reallocate(nb_nodes, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiNodalRef(PointsVertexQuad(k)(j), phi);
		for (int i = 0; i < nb_nodes; i++)
		  NodalPhiQuadVertex(k)(i, j) = phi(i);
	      }
	    
	    ValPhiQuadVertex(k).Reallocate(nb_dof, nb_quad);
	    GradPhiQuadVertex(k).Reallocate(nb_dof, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiRef(PointsVertexQuad(k)(j), phi);
		elt.ComputeGradientPhiRef(PointsVertexQuad(k)(j), grad_phi);
		for (int i = 0; i < nb_dof; i++)
		  {
		    ValPhiQuadVertex(k)(i, j) = phi(i);
		    GradPhiQuadVertex(k)(i, j) = grad_phi(i);
		  }
	      }
	  }	
      }

    // integration for elements sharing an edge
    if (type_integration == DUFFY)
      {
	nb_quad = WeightsEdgeQuad(0).GetM();
	for (int k = 0; k < 4; k++)
	  {
	    NodalPhiQuadEdge(k).Reallocate(nb_nodes, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiNodalRef(PointsEdgeQuad(k)(j), phi);
		for (int i = 0; i < nb_nodes; i++)
		  NodalPhiQuadEdge(k)(i, j) = phi(i);
	      }
	    
	    ValPhiQuadEdge(k).Reallocate(nb_dof, nb_quad);
	    GradPhiQuadEdge(k).Reallocate(nb_dof, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiRef(PointsEdgeQuad(k)(j), phi);
		elt.ComputeGradientPhiRef(PointsEdgeQuad(k)(j), grad_phi);
		for (int i = 0; i < nb_dof; i++)
		  {
		    ValPhiQuadEdge(k)(i, j) = phi(i);
		    GradPhiQuadEdge(k)(i, j) = grad_phi(i);
		  }
	      }
	  }	
      }

    // integration for the same quadrilaterals
    if (type_integration == DUFFY)
      {
	nb_quad = WeightsCommonExtQuad.GetM();
	NodalPhiQuadCommonExt.Reallocate(nb_nodes, nb_quad);
	for (int j = 0; j < nb_quad; j++)
	  {
	    elt.ComputeValuesPhiNodalRef(PointsCommonExtQuad(j), phi);
	    for (int i = 0; i < nb_nodes; i++)
	      NodalPhiQuadCommonExt(i, j) = phi(i);
	  }
	
	ValPhiQuadCommonExt.Reallocate(nb_dof, nb_quad);
	GradPhiQuadCommonExt.Reallocate(nb_dof, nb_quad);
	for (int j = 0; j < nb_quad; j++)
	  {
	    elt.ComputeValuesPhiRef(PointsCommonExtQuad(j), phi);
	    elt.ComputeGradientPhiRef(PointsCommonExtQuad(j), grad_phi);
	    for (int i = 0; i < nb_dof; i++)
	      {
		ValPhiQuadCommonExt(i, j) = phi(i);
		GradPhiQuadCommonExt(i, j) = grad_phi(i);
	      }
	  }
	
	NodalPhiQuadCommon.Reallocate(WeightsCommonQuad.GetM());
	ValPhiQuadCommon.Reallocate(WeightsCommonQuad.GetM());
	GradPhiQuadCommon.Reallocate(WeightsCommonQuad.GetM());
	for (int k = 0; k < WeightsCommonQuad.GetM(); k++)
	  {
	    nb_quad = WeightsCommonQuad(k).GetM();
	    NodalPhiQuadCommon(k).Reallocate(nb_nodes, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiNodalRef(PointsCommonQuad(k)(j), phi);
		for (int i = 0; i < nb_nodes; i++)
		  NodalPhiQuadCommon(k)(i, j) = phi(i);
	      }
	    
	    ValPhiQuadCommon(k).Reallocate(nb_dof, nb_quad);
	    GradPhiQuadCommon(k).Reallocate(nb_dof, nb_quad);
	    for (int j = 0; j < nb_quad; j++)
	      {
		elt.ComputeValuesPhiRef(PointsCommonQuad(k)(j), phi);
		elt.ComputeGradientPhiRef(PointsCommonQuad(k)(j), grad_phi);
		for (int i = 0; i < nb_dof; i++)
		  {
		    ValPhiQuadCommon(k)(i, j) = phi(i);
		    GradPhiQuadCommon(k)(i, j) = grad_phi(i);
		  }
	      }
	  }	
      }
  }
  
  
  //! computes the elementary matrix associated with the integral operator considered
  template<class T> template<class MatrixJacob>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::ComputeElemMatrix(const VectR3& PtsAi, const VectR3& PtsAj, 
		      bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
		      const VectR3& PointsNodalI, const VectR3& PointsNodalJ,
		      const VectR3& PointsQuadI, const VectR3& PointsQuadJ,
		      const VectR3& NormaleNodalI, const VectR3& NormaleNodalJ,
		      const VectR3& NormaleQuadI, const VectR3& NormaleQuadJ,
		      const MatrixJacob& MatJacobNodalI, const MatrixJacob& MatJacobNodalJ,
		      const MatrixJacob& MatJacobI, const MatrixJacob& MatJacobJ,
		      const EuclidianDistanceClass_Base<Dimension3>& var, Matrix<Real_wp>& mat_elem,
		      Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ)
  {    
    mat_elem.Fill(0); mat_elem_diagI.Fill(0); mat_elem_diagJ.Fill(0);
    
    if (param_eta_quadrature == Real_wp(0))
      {
	// all elements are treated as far elements
	ComputeFarPanelIntegral(PtsAi.GetM(), PtsAj.GetM(), Dsi, Dsj, PointsQuadI, PointsQuadJ,
				NormaleQuadI, NormaleQuadJ,
				MatJacobI, MatJacobJ, var, mat_elem,
				mat_elem_diagI, mat_elem_diagJ);
	
	return;
      }
    
    // searching common points
    int nb_pts_common = 0;
    TinyVector<int, 4> pt_common_i, pt_common_j;
    for (int i = 0; i < PtsAi.GetM(); i++)
      for (int j = 0; j < PtsAj.GetM(); j++)
	{
	  if (var.GetDistance(PtsAi(i), PtsAj(j)) <= R3::threshold)
	    {
	      pt_common_i(nb_pts_common) = i;
	      pt_common_j(nb_pts_common) = j;
	      nb_pts_common++;
	      break;
	    }
	}
    
    //if (nb_pts_common != 2)
    //return;
    
    if (nb_pts_common == 0)
      {
	// far elements
	ComputeFarPanelIntegral(PtsAi.GetM(), PtsAj.GetM(), Dsi, Dsj, PointsQuadI, PointsQuadJ,
				NormaleQuadI, NormaleQuadJ,
				MatJacobI, MatJacobJ, var, mat_elem,
				mat_elem_diagI, mat_elem_diagJ);
      }
    else if (nb_pts_common == 1)
      {
	// elements with a common vertex
	ComputeVertexPanelIntegral(pt_common_i(0), pt_common_j(0), PtsAi.GetM(), PtsAj.GetM(),
				   curved, Dsi, Dsj, PointsNodalI, PointsNodalJ,
				   NormaleNodalI, NormaleNodalJ,
				   MatJacobNodalI, MatJacobNodalJ, var, mat_elem,
				   mat_elem_diagI, mat_elem_diagJ);
      }
    else if (nb_pts_common == 2)
      {
	// elements sharing an edge
	ComputeEdgePanelIntegral(pt_common_i(0), pt_common_i(1), pt_common_j(0), pt_common_j(1),
				 PtsAi.GetM(), PtsAj.GetM(), curved, Dsi, Dsj,
				 PointsNodalI, PointsNodalJ, NormaleNodalI, NormaleNodalJ,
				 MatJacobNodalI, MatJacobNodalJ, var, mat_elem,
				 mat_elem_diagI, mat_elem_diagJ);
      }
    else if (nb_pts_common > 2)
      {
	// same element
	ComputeSamePanelIntegral(PtsAi.GetM(), curved, Dsi, PointsNodalI, NormaleNodalI,
				 MatJacobNodalI, var, mat_elem, mat_elem_diagI, mat_elem_diagJ);
      }
  }
  
  
  //! computes elementary of mass matrix \int u(x) v(x) dx
  template<class T> 
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::ComputeMassMatrix(const VectR3& PtsAi, bool curved,
		      const VectReal_wp& Dsi, Matrix<Real_wp>& mat_elem)
  {
    mat_elem.Fill(0);    
    Real_wp val;

    if (PtsAi.GetM() == 3)
      {
	const ElementReference<Dimension2, 1>& Fb = *fe_tri;
	for (int i = 0; i < WeightsTriFar.GetM(); i++)
	  {
	    val = WeightsTriFar(i)*Dsi(i);
	    for (int k = 0; k < Fb.GetNbDof(); k++)
	      for (int n = 0; n < Fb.GetNbDof(); n++)
		mat_elem(k, n) += val*ValPhiTriFar(k, i)*ValPhiTriFar(n, i);
	  }
      }
    else
      {
	const ElementReference<Dimension2, 1>& Fb = *fe_quad;
	for (int i = 0; i < WeightsQuadFar.GetM(); i++)
	  {
	    val = WeightsQuadFar(i)*Dsi(i);
	    for (int k = 0; k < Fb.GetNbDof(); k++)
	      for (int n = 0; n < Fb.GetNbDof(); n++)
		mat_elem(k, n) += val*ValPhiQuadFar(k, i)*ValPhiQuadFar(n, i);
	  }
      }      
  }
  
  
  //! computes elementary of stiffness matrix \int 1/k_\eps^2 \grad u(x) \cdot \grad v(x) dx
  template<class T> template<class MatrixJacob>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension2>
  ::ComputeStiffnessMatrix(const VectR3& PtsAi, bool curved,
			   const VectReal_wp& Dsi, const VectR3& NormaleQuad,
			   const MatrixJacob& MatJacobI,
			   const Complex_wp& kinf, const Real_wp& epsilon,
			   Matrix<Complex_wp>& mat_elem)
  {
    mat_elem.Fill(0);
    
    Complex_wp val; R3 grad_phi_k, grad_phi_n;

    if (PtsAi.GetM() == 3)
      {
	const ElementReference<Dimension2, 1>& Fb = *fe_tri;
	for (int i = 0; i < WeightsTriFar.GetM(); i++)
	  {
	    Real_wp eps = 0;
	    
	    Complex_wp keps = kinf + Iwp*eps;
	    val = WeightsTriFar(i)*Dsi(i)/(keps*keps);
	    for (int k = 0; k < Fb.GetNbDof(); k++)
	      for (int n = 0; n < Fb.GetNbDof(); n++)
		{
		  Mlt(MatJacobI(i), GradPhiTriFar(k, i), grad_phi_k);
		  Mlt(MatJacobI(i), GradPhiTriFar(n, i), grad_phi_n);
		  mat_elem(k, n) += val*DotProd(grad_phi_k, grad_phi_n);
		}
	  }
      }
    else
      {
	const ElementReference<Dimension2, 1>& Fb = *fe_quad;
	for (int i = 0; i < WeightsQuadFar.GetM(); i++)
	  {
	    Real_wp eps = 0;
	    
	    Complex_wp keps = kinf + Iwp*eps;
	    val = WeightsQuadFar(i)*Dsi(i)/(keps*keps);
	    for (int k = 0; k < Fb.GetNbDof(); k++)
	      for (int n = 0; n < Fb.GetNbDof(); n++)
		{
		  Mlt(MatJacobI(i), GradPhiQuadFar(k, i), grad_phi_k);
		  Mlt(MatJacobI(i), GradPhiQuadFar(n, i), grad_phi_n);
		  mat_elem(k, n) += val*DotProd(grad_phi_k, grad_phi_n);
		}
	  }
      }      
  }
  
}

#define MONTJOIE_FILE_SINGULAR_INTEGRAL_EQUATION_3D_CXX
#endif
