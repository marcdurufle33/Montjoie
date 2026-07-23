#ifndef MONTJOIE_FILE_SINGULAR_INTEGRAL_EQUATION_2D_CXX

#include "SingularIntegration2D.hxx"

namespace Montjoie
{

  template<class Dimension>
  Real_wp EuclidianDistanceClass_Base<Dimension>::GetDistance(const R_N& x, const R_N& y) const
  {
    return x.Distance(y);
  }


  template<class Dimension>
  typename Dimension::R_N EuclidianDistanceClass_Base<Dimension>
  ::GetDiff(const R_N& x, const R_N& y) const
  {
    return y-x;
  }
  

  //! Evaluates the double integral for panels which are far from each other
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::ComputeFarPanelIntegral(bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			    const VectR2& PointsQuadI, const VectR2& PointsQuadJ,
			    const VectR2& NormaleQuadI, const VectR2& NormaleQuadJ,
			    const EuclidianDistanceClass_Base<Dimension2>& var, Matrix<Real_wp>& mat_elem,
			    Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ)
  {
    const ElementReference<Dimension1, 1>& Fb = *fe_class;
    mat_elem.Fill(0); mat_elem_diagI.Fill(0); mat_elem_diagJ.Fill(0);
    int p = Fb.GetOrder();
    int r = Fb.GetNbPointsQuadratureInside();
    
    Real_wp val;
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	{
	  val = Fb.Weights(i)*Fb.Weights(j);
	  val *= pow(var.GetDistance(PointsQuadI(i), PointsQuadJ(j)), beta);
          if (type_operator == OPERATOR_GRAD)
	      val *= DotProd(NormaleQuadI(i), NormaleQuadJ(j));
          
	  switch (type_operator)
	    {
	    case OPERATOR_GRAD :
	    case OPERATOR_DS :
	      for (int k = 0; k <= p; k++)
		for (int n = 0; n <= p; n++)
		  mat_elem(k, n) += val*GradPhiFar(k, i)*GradPhiFar(n, j);
	      
	      break;
	    case OPERATOR_DIFF :
	      val *= Dsi(i)*Dsj(j);
	      for (int k = 0; k <= p; k++)
		for (int n = 0; n <= p; n++)
		  {
		    mat_elem(k, n) -= val*ValPhiFar(k, i)*ValPhiFar(n, j);
		    mat_elem_diagI(k, n) += val*ValPhiFar(k, i)*ValPhiFar(n, i);
		    mat_elem_diagJ(k, n) += val*ValPhiFar(k, j)*ValPhiFar(n, j);
		  }
	      
	      break;
	    }		    
	}
  }
  

  //! Evaluates the double integral for panels which are close from each other
  //! These panels are not joints neither identical
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::ComputeClosePanelIntegral(bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			      const VectR2& PointsQuadI, const VectR2& PointsQuadJ,
			      const VectR2& NormaleQuadI, const VectR2& NormaleQuadJ,
			      const EuclidianDistanceClass_Base<Dimension2>& var, Matrix<Real_wp>& mat_elem,
			      Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ)
  {
    const ElementReference<Dimension1, 1>& Fb = *fe_class;
    mat_elem.Fill(0); mat_elem_diagI.Fill(0); mat_elem_diagJ.Fill(0);
    int p = Fb.GetOrder();
    int r = quadrature_proche.GetOrder();
    Real_wp val;
    
    // interpolation
    VectR2 PointsQuadInterp(r+1), PointsQuadJinterp(r+1);
    VectR2 NormaleQuadInterp(r+1), NormaleQuadJinterp(r+1);
    VectReal_wp DsiInterp(r+1), DsjInterp(r+1);
    DsiInterp.Fill(0); DsjInterp.Fill(0);
    for (int i = 0; i <= p; i++)
      for (int j = 0; j <= r; j++)
	{
	  val = QuadPhiClose(i, j);
	  DsiInterp(j) += val*Dsi(i); DsjInterp(j) += val*Dsj(i);
	  Add(val, PointsQuadI(i), PointsQuadInterp(j));
	  Add(val, PointsQuadJ(i), PointsQuadJinterp(j));
	  Add(val, NormaleQuadI(i), NormaleQuadInterp(j));
	  Add(val, NormaleQuadJ(i), NormaleQuadJinterp(j));
	}
    
    // integration with r+1 points
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  val = quadrature_proche.Weights(i)*quadrature_proche.Weights(j);
	  val *= pow(var.GetDistance(PointsQuadInterp(i), PointsQuadJinterp(j)), beta);
	  if (type_operator == OPERATOR_GRAD)
	      val *= DotProd(NormaleQuadInterp(i), NormaleQuadJinterp(j));
          
	  switch (type_operator)
	    {
	    case OPERATOR_GRAD :
	    case OPERATOR_DS :
	      for (int k = 0; k <= p; k++)
		for (int n = 0; n <= p; n++)
		  mat_elem(k, n) += val*GradPhiClose(k, i)*GradPhiClose(n, j);
	      
	      break;
	    case OPERATOR_DIFF :
	      val *= DsiInterp(i)*DsjInterp(j);
	      for (int k = 0; k <= p; k++)
		for (int n = 0; n <= p; n++)
		  {
		    mat_elem(k, n) -= val*ValPhiClose(k, i)*ValPhiClose(n, j);
		    mat_elem_diagI(k, n) += val*ValPhiClose(k, i)*ValPhiClose(n, i);
		    mat_elem_diagJ(k, n) += val*ValPhiClose(k, j)*ValPhiClose(n, j);
		  }
	      
	      break;
	    }		    
	}
  }
  
  
  //! Evaluates the double integral for identical panels
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::ComputeIdenticalPanelIntegral(bool curved, const VectReal_wp& Dsi, const VectR2& PointsQuadI,
				  const VectR2& NormaleQuadI, const EuclidianDistanceClass_Base<Dimension2>& var,
				  Matrix<Real_wp>& mat_elem)
  {
    if (curved)
      {
	if (type_integration == LENOIR_SALLES)
	  {
	    cout << "Not implemented" << endl;
	    cout << "Select a numerical integration" << endl;
	    abort();
	  }
	
        int p = fe_class->GetOrder();
        R2 pt_k, normale_k, pt_n, normale_n; Real_wp val, ds_k, ds_n, dist;
	mat_elem.Fill(0);
        // loop over quadrature points of the outer integral
        for (int k = 0; k <= gauss_ext.GetOrder(); k++)
	  {
            // interpolation
            ds_k = 0; pt_k.Zero(); normale_k.Zero();
            for (int i = 0; i <= p; i++)
              {
                val = QuadPhiExt(i, k);
                ds_k += val*Dsi(i);
                Add(val, PointsQuadI(i), pt_k);
                Add(val, NormaleQuadI(i), normale_k);
              }
            
            // loop over quadrature points of the inner integral
	    Real_wp val_int(0);
	    for (int n = 0; n < PointsInt(k).GetM(); n++)
	      {
                // interpolation
                ds_n = 0; pt_n.Zero(); normale_n.Zero();
                for (int i = 0; i <= p; i++)
                  {
                    val = QuadPhiInt(k)(i, n);
                    ds_n += val*Dsi(i);
                    Add(val, PointsQuadI(i), pt_n);
                    Add(val, NormaleQuadI(i), normale_n);
                  }
                
                dist = var.GetDistance(pt_k, pt_n);
		val_int = gauss_ext.Weights(k)*WeightsInt(k)(n)*pow(dist, beta);
		
		if (type_operator != OPERATOR_DIFF)
		  {
                    if (type_operator == OPERATOR_GRAD)
                      val_int *= DotProd(normale_k, normale_n);
                    
		    for (int i = 0; i <= p; i++)
		      for (int j = 0; j <= p; j++)
			mat_elem(i, j) += val_int*GradPhiInt(k)(i, n)*GradPhiExt(j, k);
		  }
		else
		  {
                    val_int *= ds_k*ds_n;
		    for (int i = 0; i <= p; i++)
		      for (int j = 0; j <= p; j++)
			mat_elem(i, j) += val_int*(ValPhiInt(k)(i, n) - ValPhiExt(i, k))
			  *(ValPhiInt(k)(j, n) - ValPhiExt(j, k));
		  }
	      }
	  }
      }
    else
      {
	mat_elem = mat_elem_reference_identical;
	if (type_operator != OPERATOR_DIFF)
	  Mlt(pow(Dsi(0), beta), mat_elem);
	else
	  Mlt(pow(Dsi(0), beta+2), mat_elem);
      }
  }
  
  
  //! Evaluates the double integral for joint panels
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::ComputeJointPanelIntegral(bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
			      const VectR2& PointsQuadI, const VectR2& PointsQuadJ,
			      const VectR2& NormaleQuadI, const VectR2& NormaleQuadJ,
			      const R2& u, const R2& v,
			      const EuclidianDistanceClass_Base<Dimension2>& var, Matrix<Real_wp>& mat_elem,
			      Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ)
  {
    int r = fe_class->GetOrder();
    Real_wp one(1), s, t; mat_elem.Fill(0); mat_elem_diagI.Fill(0); mat_elem_diagJ.Fill(0);
    if (type_integration == LENOIR_SALLES)
      {
	if (curved)
	  {
	    cout << "Not implemented" << endl;
	    cout << "Select a numerical integration" << endl;
	    abort();
	  }
	
	// methode semi-analytique avec n_ext+1 points
	int n_ext = quadrature_joint.GetOrder();
	if (type_operator != OPERATOR_DIFF)
	  {
	    VectReal_wp IntElemS(r), IntElemT(r);
	    IntElemS.Fill(0); IntElemT.Fill(0);
	    Real_wp val_int(0), Si;
	    for (int k = 0; k <= n_ext; k++)
	      {
		s = quadrature_joint.Points(k);
		val_int = pow(sqrt(square(s*u(0)-v(0)) + square(s*u(1)-v(1))), beta)
		  *quadrature_joint.Weights(k);
		
		Si = one;
		for (int i = 0; i <= r-1; i++)
		  {
		    IntElemS(i) += val_int*Si;
		    Si *= s;
		  }
		
		t = s;
		val_int = pow(sqrt(square(u(0)-t*v(0)) + square(u(1)-t*v(1))), beta)
		  *quadrature_joint.Weights(k);
		
		Si = one;
		for (int i = 0; i <= r-1; i++)
		  {
		    IntElemT(i) += val_int*Si;
		    Si *= t;
		  }
	      }
	    
	    for (int i = 0; i <= r; i++)
	      for (int j = 0; j <= r; j++)
		{
		  Real_wp val(0);
		  for (int k = 0; k < r; k++)
		    for (int n = 0; n < r; n++)
		      val += dPhi(i)(k)*dPhi(j)(n)
			*InvBetaPlusN(k+n+1)*(IntElemS(k)+IntElemT(n));
		  
		  if (type_operator == OPERATOR_GRAD)
		    val *= DotProd(NormaleQuadI(0), NormaleQuadJ(0));
		  
		  mat_elem(i, j) = -val;
		}
	  }
	else
	  {
	    VectReal_wp IntElemS(2*r+1), IntElemT(2*r+1);
	    IntElemS.Fill(0); IntElemT.Fill(0);
	    Real_wp val_int(0), Si;
	    for (int k = 0; k <= n_ext; k++)
	      {
		s = quadrature_joint.Points(k);
		val_int = pow(sqrt(square(s*u(0)-v(0)) + square(s*u(1)-v(1))), beta)
		  *quadrature_joint.Weights(k);
		
		Si = one;
		for (int i = 0; i <= 2*r; i++)
		  {
		    IntElemS(i) += val_int*Si;
		    Si *= s;
		  }
		
		t = s;
		val_int = pow(sqrt(square(u(0)-t*v(0)) + square(u(1)-t*v(1))), beta)
		  *quadrature_joint.Weights(k);
		
		Si = one;
		for (int i = 0; i <= 2*r; i++)
		  {
		    IntElemT(i) += val_int*Si;
		    Si *= t;
		  }
	      }
	    
	    MultivariatePolynomial<T> PolI, PolJ, PolI2, PolJ2, Q;
	    for (int i = 0; i <= r; i++)
	      for (int j = 0; j <= r; j++)
		{
		  Real_wp val(0);
		  PolI.SetOrder(2, r); PolI2.SetOrder(2, r);
		  PolJ.SetOrder(2, r); PolJ2.SetOrder(2, r);
		  PolI.Fill(0); PolJ.Fill(0); PolI2.Fill(0); PolJ2.Fill(0);
		  if (i == 0)
		    {
		      for (int n = 1; n <= r; n++)
			{
			  PolI(n, 0) = Phi(i)(n);
			  PolI(0, n) = -Phi(i)(n);
			  PolI2(n, 0) = Phi(i)(n);
			  PolI2(0, n) = -Phi(i)(n);
			}
		    }
		  else
		    {
		      for (int n = 0; n <= r; n++)
			{
			  PolI(n, 0) = Phi(i)(n);
			  PolI2(0, n) = -Phi(i)(n);
			}
		    }
		  
		  if (j == 0)
		    {
		      for (int n = 1; n <= r; n++)
			{
			  PolJ(n, 0) = Phi(j)(n);
			  PolJ(0, n) = -Phi(j)(n);
			  PolJ2(n, 0) = Phi(j)(n);
			  PolJ2(0, n) = -Phi(j)(n);
			}
		    }
		  else
		    {
		      for (int n = 0; n <= r; n++)
			{
			  PolJ(0, n) = -Phi(j)(n);
			  PolJ2(n, 0) = Phi(j)(n);
			}
		    }
		  
		  Q = PolI*PolJ;
		  for (int k = 0; k <= Q.GetOrder(); k++)
		    for (int n = 0; n <= Q.GetOrder()-k; n++)
		      if ((k != 0) || (n != 0))
			val += Q(k, n)*InvBetaPlusN(k+n+1)*(IntElemS(k)+IntElemT(n));
		  
		  mat_elem(i, j) = val*Dsi(0)*Dsj(0);

		  if ((i != 0) && (j != 0))
		    {
		      Q = PolI*PolJ2;
		      val = Real_wp(0);
		      for (int k = 0; k <= Q.GetOrder(); k++)
			for (int n = 0; n <= Q.GetOrder()-k; n++)
			  if ((k != 0) || (n != 0))
			    val += Q(k, n)*InvBetaPlusN(k+n+1)*(IntElemS(k)+IntElemT(n));
		      
		      mat_elem_diagI(i, j) = val*Dsi(0)*Dsj(0);
		      
		      Q = PolI2*PolJ;
		      val = Real_wp(0);
		      for (int k = 0; k <= Q.GetOrder(); k++)
			for (int n = 0; n <= Q.GetOrder()-k; n++)
			  if ((k != 0) || (n != 0))
			    val += Q(k, n)*InvBetaPlusN(k+n+1)*(IntElemS(k)+IntElemT(n));
		      
		      mat_elem_diagJ(i, j) = val*Dsi(0)*Dsj(0);
		    }
		}
	  }
      }
    else
      {
	// methode purement numerique
	
	// interpolation
	int r_int = quadrature_inner_joint.GetOrder();
	int r_ext = quadrature_outer_joint.GetOrder();
	VectR2 NormaleQuadInterp(r_int+1), NormaleQuadJinterp(r_ext+1);
	VectReal_wp DsiInterp(r_int+1), DsjInterp(r_ext+1);
	DsiInterp.Fill(0); DsjInterp.Fill(0);
	Real_wp val;
	for (int i = 0; i <= r; i++)
	  {
	    for (int j = 0; j <= r_int; j++)
	      {
		val = QuadPhiJoint(i, j);
		DsiInterp(j) += val*Dsi(i);
		Add(val, NormaleQuadI(i), NormaleQuadInterp(j));
	      }
	    
	    for (int j = 0; j <= r_ext; j++)
	      {
		val = QuadPhiJointExt(i, j);
		DsjInterp(j) += val*Dsj(i);
		Add(val, NormaleQuadJ(i), NormaleQuadJinterp(j));
	      }
	  }
	
	// quadrature
	for (int k = 0; k <= quadrature_outer_joint.GetOrder(); k++)
	  { 
	    t = quadrature_outer_joint.Points(k);
	    Real_wp val_int(0);
	    for (int n = 0; n <= quadrature_inner_joint.GetOrder(); n++)
	      {
		s = quadrature_inner_joint.Points(n);
		val_int = quadrature_outer_joint.Weights(k)*quadrature_inner_joint.Weights(n)
		  *pow(sqrt(square(s*u(0)-t*v(0)) + square(s*u(1)-t*v(1))), beta);
		
                if (type_operator == OPERATOR_GRAD)
		    val_int *= DotProd(NormaleQuadJinterp(k), NormaleQuadInterp(n));
                
		switch (type_operator)
		  {
		  case OPERATOR_GRAD :
		  case OPERATOR_DS :
		    for (int i = 0; i <= r; i++)
		      for (int j = 0; j <= r; j++)
			mat_elem(i, j) -= val_int*GradPhiJoint(i, n)*GradPhiJointExt(j, k);
		    break;
		  case OPERATOR_DIFF :
		    val_int *= DsiInterp(n)*DsjInterp(k);
		    mat_elem(0, 0) += val_int*square(ValPhiJoint(0, n)-ValPhiJointExt(0, k));
		    
		    for (int i = 1; i <= r; i++)
		      {
			mat_elem(i, 0) += val_int*(ValPhiJoint(0, n)-ValPhiJointExt(0, k))*ValPhiJoint(i, n);
			mat_elem(0, i) -= val_int*(ValPhiJoint(0, n)-ValPhiJointExt(0, k))*ValPhiJointExt(i, k);
		      }
		    
		    for (int i = 1; i <= r; i++)
		      for (int j = 1; j <= r; j++)
			{
			  mat_elem(i, j) -= val_int*ValPhiJoint(i, n)*ValPhiJointExt(j, k);
			  mat_elem_diagI(i, j) += val_int*ValPhiJoint(i, n)*ValPhiJoint(j, n);
			  mat_elem_diagJ(i, j) += val_int*ValPhiJointExt(i, k)*ValPhiJointExt(j, k);
			}
		  }
	      }
	  }
      }
  }
    
  
  //! default constructor
  template<class T>
  SingularDoubleQuadratureGalerkin_Base<T, Dimension1>::SingularDoubleQuadratureGalerkin_Base()
  {
    fe_class = NULL;
    param_eta_quadrature = 0.0;
    type_operator = OPERATOR_DS;
    beta = 0.5;
    type_integration = LENOIR_SALLES;      
  }
  
  
  //! selects the operator to be integrated and coefficient beta
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::SetOperator(int type, const Real_wp& b)
  {
    type_operator = type;
    beta = b;
  }
  
  
  //! selects an analytical integration
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::SetAnalyticalIntegration(const Real_wp& eta, int rp, int rj)
  {
    type_integration = LENOIR_SALLES;
    param_eta_quadrature = eta;
    
    quadrature_proche.ConstructQuadrature(rp);
    quadrature_joint.ConstructQuadrature(rj);
  }
  
  
  //! selects a numerical integration
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::SetNumericalIntegration(const Real_wp& eta, int rp, int rj, int rj_ext, int n_int, int n_ext)
  {
    type_integration = GAUSS_SQUARED;
    param_eta_quadrature = eta;
    
    Globatto<Real_wp> gauss;
    quadrature_proche.ConstructQuadrature(rp);

    // for joint panels
    quadrature_inner_joint.ConstructQuadrature(rj, gauss.QUADRATURE_GAUSS_SQUARED);
    
    gauss.ConstructQuadrature(rj_ext, gauss.QUADRATURE_GAUSS_SQUARED);
    VectReal_wp Points(2*rj_ext+2), Weights(2*rj_ext+2);
    for (int i = 0; i <= rj_ext; i++)
      {
	T xi_i = gauss.Points(i), one(1);
	Points(i) = xi_i/2;
	Weights(i) = gauss.Weights(i)/2;
	
	Points(i+rj_ext+1) = (one-xi_i/2);
	Weights(i+rj_ext+1) = gauss.Weights(i)/2;      
      } 
    
    quadrature_outer_joint.AffectPoints(Points);
    quadrature_outer_joint.AffectWeights(Weights);      

    // for identical panels
    gauss_int.ConstructQuadrature(n_int, gauss_int.QUADRATURE_GAUSS_SQUARED);
    
    gauss.ConstructQuadrature(n_ext, gauss.QUADRATURE_GAUSS_SQUARED);
    Points.Reallocate(2*n_ext+2); Weights.Reallocate(2*n_ext+2);
    Real_wp one(1);
    for (int i = 0; i <= n_ext; i++)
      {
	Real_wp xi_i = gauss.Points(i);
	Points(i) = xi_i/2;
	Weights(i) = gauss.Weights(i)/2;
	
	Points(i+n_ext+1) = (one-xi_i/2);
	Weights(i+n_ext+1) = gauss.Weights(i)/2;      
      } 
    
    gauss_ext.AffectPoints(Points); gauss_ext.AffectWeights(Weights);
  }
  
  
  //! initialisation of arrays needed to compute elementary matrices
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::SetFiniteElement(const ElementReference<Dimension1, 1>& elt)
  {
    fe_class = const_cast<ElementReference<Dimension1, 1>*>(&elt);  
    int r = elt.GetOrder();
    
    quad_basis.AffectPoints(elt.Points());
    quad_basis.ComputeGradPhi();
    
    // evaluation of reference element matrix for identical panels
    // computation of the singular integral :
    // \int_0^1 \int_0^1 |s-t|^\beta f(s, t) ds dt
    
    Matrix<Real_wp>& mat_elem = mat_elem_reference_identical;
    Real_wp xi_i, one(1);
    mat_elem.Reallocate(r+1, r+1);
    mat_elem.Fill(0);
    if (type_integration == GAUSS_SQUARED)
      {
	int n_int = gauss_int.GetOrder(), n_ext = gauss_ext.GetOrder();
	
	// values of basis functions on outer quadrature points
	GradPhiExt.Reallocate(r+1, n_ext+1);
        ValPhiExt.Reallocate(r+1, n_ext+1);
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= n_ext; j++)
	    {
	      ValPhiExt(i, j) = elt.GetValuePhi1D(i, gauss_ext.Points(j));
	      GradPhiExt(i, j) = elt.GetGradientPhi1D(i, gauss_ext.Points(j));
	    }
	
	// values of basis functions on inner quadrature points
	GradPhiInt.Reallocate(n_ext+1); ValPhiInt.Reallocate(n_ext+1);
	PointsInt.Reallocate(n_ext+1);
        WeightsInt.Reallocate(n_ext+1);
	
	for (int i = 0; i <= n_ext; i++)
	  {
	    xi_i = gauss_ext.Points(i);
	    // cas avec juste deux intervalles [0, xi_i] et [xi_i, 1]
	    PointsInt(i).Reallocate(2*(n_int+1));
	    WeightsInt(i).Reallocate(2*(n_int+1));
	    for (int j = 0; j <= n_int; j++)
	      {
		PointsInt(i)(j) = (one-gauss_int.Points(j))*xi_i;
		PointsInt(i)(j+n_int+1) = gauss_int.Points(j)*(one-xi_i) + xi_i;
		WeightsInt(i)(j) = gauss_int.Weights(j)*xi_i;
		WeightsInt(i)(j+n_int+1) = gauss_int.Weights(j)*(one-xi_i);
	      }
	    
	    GradPhiInt(i).Reallocate(r+1, PointsInt(i).GetM());
	    ValPhiInt(i).Reallocate(r+1, PointsInt(i).GetM());
	    for (int k = 0; k <= r; k++)
	      for (int j = 0; j < PointsInt(i).GetM(); j++)
		{
		  GradPhiInt(i)(k, j) = elt.GetGradientPhi1D(k, PointsInt(i)(j));
		  ValPhiInt(i)(k, j) = elt.GetValuePhi1D(k, PointsInt(i)(j));
		}	  	  
	  }
	
	// numerical integration
	for (int k = 0; k <= n_ext; k++)
	  {
	    xi_i = gauss_ext.Points(k);
	    Real_wp val_int(0);
	    for (int n = 0; n < PointsInt(k).GetM(); n++)
	      {
		val_int = gauss_ext.Weights(k)*WeightsInt(k)(n)*pow(abs(PointsInt(k)(n)-xi_i), beta);
		
		if (type_operator != OPERATOR_DIFF)
		  {
		    for (int i = 0; i <= r; i++)
		      for (int j = 0; j <= r; j++)
			mat_elem(i, j) += val_int*GradPhiInt(k)(i, n)*GradPhiExt(j, k);
		  }
		else
		  {
		    for (int i = 0; i <= r; i++)
		      for (int j = 0; j <= r; j++)
			mat_elem(i, j) += val_int*(ValPhiInt(k)(i, n) - ValPhiExt(i, k))
			  *(ValPhiInt(k)(j, n) - ValPhiExt(j, k));
		  }
	      }
	  }
	
        QuadPhiExt.Reallocate(r+1, gauss_ext.GetOrder()+1);
        for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= gauss_ext.GetOrder(); j++)
            QuadPhiExt(i, j) = quad_basis.EvaluatePhi(i, gauss_ext.Points(j));
	
        QuadPhiInt.Reallocate(n_ext+1);
	for (int k = 0; k <= n_ext; k++)
          {
            QuadPhiInt(k).Reallocate(r+1, PointsInt(k).GetM());
            for (int i = 0; i <= r; i++)
              for (int j = 0; j < PointsInt(k).GetM(); j++)
                QuadPhiInt(k)(i, j) = quad_basis.EvaluatePhi(i, PointsInt(k)(j));
          }
        
	// quadrature rules for joint integrals
	ValPhiJoint.Reallocate(r+1, quadrature_inner_joint.GetOrder()+1);
	QuadPhiJoint.Reallocate(r+1, quadrature_inner_joint.GetOrder()+1);
	GradPhiJoint.Reallocate(r+1, quadrature_inner_joint.GetOrder()+1);
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= quadrature_inner_joint.GetOrder(); j++)
	    {
	      ValPhiJoint(i, j) = elt.GetValuePhi1D(i, quadrature_inner_joint.Points(j));
	      QuadPhiJoint(i, j) = quad_basis.EvaluatePhi(i, quadrature_inner_joint.Points(j));
	      GradPhiJoint(i, j) = elt.GetGradientPhi1D(i, quadrature_inner_joint.Points(j));
	    }
	
	ValPhiJointExt.Reallocate(r+1, quadrature_outer_joint.GetOrder()+1);
	QuadPhiJointExt.Reallocate(r+1, quadrature_outer_joint.GetOrder()+1);
	GradPhiJointExt.Reallocate(r+1, quadrature_outer_joint.GetOrder()+1);
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= quadrature_outer_joint.GetOrder(); j++)
	    {
	      ValPhiJointExt(i, j) = elt.GetValuePhi1D(i, quadrature_outer_joint.Points(j));
	      QuadPhiJointExt(i, j) = quad_basis.EvaluatePhi(i, quadrature_outer_joint.Points(j));
	      GradPhiJointExt(i, j) = elt.GetGradientPhi1D(i, quadrature_outer_joint.Points(j));
	    }
	
      }
    else
      {
	elt.GetPolynomialFunctions(Phi, dPhi);
	
	// analytical value of the integral
	// \int_0^1 (1-t)^\beta t^j dt = j! / ( (beta+1+j) (beta+j) ... (beta+1) )
	Vector<Real_wp> IntElem(2*r+1);
	InvBetaPlusN.Reallocate(2*r+4);
	for (int i = 0; i < InvBetaPlusN.GetM(); i++)
	  InvBetaPlusN(i) = one/(beta+i+1);
	
	if (type_operator != OPERATOR_DIFF)
	  {
	    IntElem(0) = one/(beta+one);
	    for (int i = 0; i < IntElem.GetM()-1; i++)
	      IntElem(i+1) = IntElem(i)*(i+1)*InvBetaPlusN(i+1);
	  }
	else
	  {
	    IntElem(0) = one/(beta+3);
	    for (int i = 0; i < IntElem.GetM()-1; i++)
	      IntElem(i+1) = IntElem(i)*(i+1)*InvBetaPlusN(i+3);    
	  }
	
	MultivariatePolynomial<Real_wp> PolI, PolJ, Q;
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= r; j++)
	    {
	      Real_wp val(0);
	      if (type_operator != OPERATOR_DIFF)
		{
		  for (int k = 0; k < r; k++)
		    for (int n = 0; n < r; n++)
		      val += dPhi(i)(k)*dPhi(j)(n)*InvBetaPlusN(k+n+1)*(IntElem(k)+IntElem(n));
		}
	      else
		{
		  PolI.SetOrder(2, r-1);
		  PolJ.SetOrder(2, r-1);
		  for (int n = 1; n <= r; n++)
		    for (int k = 0; k < n; k++)
		      {
			PolI(n-1-k, k) = Phi(i)(n);
			PolJ(n-1-k, k) = Phi(j)(n);
		      }
		  
		  Q = PolI*PolJ;
		  for (int k = 0; k <= Q.GetOrder(); k++)
		    for (int n = 0; n <= Q.GetOrder()-k; n++)
		      if (Q(k, n) != Real_wp(0))
			val += Q(k, n)*InvBetaPlusN(k+n+3)*(IntElem(k)+IntElem(n));
		}
	      
	      mat_elem(i, j) = val;
	    }
      }
    
    ValPhiClose.Reallocate(r+1, quadrature_proche.GetOrder()+1);
    QuadPhiClose.Reallocate(r+1, quadrature_proche.GetOrder()+1);
    GradPhiClose.Reallocate(r+1, quadrature_proche.GetOrder()+1);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= quadrature_proche.GetOrder(); j++)
	{
	  ValPhiClose(i, j) = elt.GetValuePhi1D(i, quadrature_proche.Points(j));
	  QuadPhiClose(i, j) = quad_basis.EvaluatePhi(i, quadrature_proche.Points(j));
	  GradPhiClose(i, j) = elt.GetGradientPhi1D(i, quadrature_proche.Points(j));
	}
    
    ValPhiFar.Reallocate(r+1, elt.GetOrder()+1);
    GradPhiFar.Reallocate(r+1, elt.GetOrder()+1);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= elt.GetOrder(); j++)
	{
	  ValPhiFar(i, j) = elt.GetValuePhi1D(i, elt.Points(j));
	  GradPhiFar(i, j) = elt.GetGradientPhi1D(i, elt.Points(j));
	}
    
  }
  
  
  //! computes the elementary matrix \int_Ai^Bi \int_Aj^Bj |X(s)-X(t)|^\beta f(s, t) ds dt  
  template<class T> template<class MatrixJacob>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::ComputeElemMatrix(const VectR2& PtsAi, const VectR2& PtsAj, 
		      bool curved, const VectReal_wp& Dsi, const VectReal_wp& Dsj,
		      const VectR2& PointsNodalI, const VectR2& PointsNodalJ,
		      const VectR2& PointsQuadI, const VectR2& PointsQuadJ,
		      const VectR2& NormaleNodalI, const VectR2& NormaleNodalJ,
		      const VectR2& NormaleQuadI, const VectR2& NormaleQuadJ,
		      const MatrixJacob& MatJacobNodalI, const MatrixJacob& MatJacobNodalJ,
		      const MatrixJacob& MatJacobI, const MatrixJacob& MatJacobJ,
		      const EuclidianDistanceClass_Base<Dimension2>& var, Matrix<Real_wp>& mat_elem,
		      Matrix<Real_wp>& mat_elem_diagI, Matrix<Real_wp>& mat_elem_diagJ)
  {
    R2 Ai(PtsAi(0)), Bi(PtsAi(1));
    R2 Aj(PtsAj(0)), Bj(PtsAj(1));
    int r = fe_class->GetOrder();
    if (param_eta_quadrature == Real_wp(0))
      {
	// all elements are treated as far elements
	ComputeFarPanelIntegral(curved, Dsi, Dsj, PointsQuadI, PointsQuadJ, 
				NormaleQuadI, NormaleQuadJ, var, mat_elem,
				mat_elem_diagI, mat_elem_diagJ);
	return;
      }
    
    if (var.GetDistance(Ai, Aj) <= R2::threshold)
      {
	if (var.GetDistance(Bi, Bj) <= R2::threshold)
	  {
	    // identical elements
	    ComputeIdenticalPanelIntegral(curved, Dsi, PointsQuadI, NormaleQuadI,
					  var, mat_elem);	      
	    
	    return;
	  }
	
	// joint panels
	R2 vec_u = var.GetDiff(Ai, Bi); R2 vec_v = var.GetDiff(Aj, Bj);
	ComputeJointPanelIntegral(curved, Dsi, Dsj, PointsQuadI, PointsQuadJ,
				  NormaleQuadI, NormaleQuadJ, vec_u, vec_v, var, mat_elem,
				  mat_elem_diagI, mat_elem_diagJ);
	
	return;
      }
    
    if (var.GetDistance(Ai, Bj) <= R2::threshold )
      {
	// joints panels
	R2 vec_u = var.GetDiff(Ai, Bi); R2 vec_v = var.GetDiff(Bj, Aj);
	VectR2 PointsQuadJR(r+1), NormaleQuadJR(r+1);
	VectReal_wp DsjR(r+1);
	for (int i = 0; i <= r; i++)
	  {
	    DsjR(i) = Dsj(r-i);
	    PointsQuadJR(i) = PointsQuadJ(r-i);
	    NormaleQuadJR(i) = NormaleQuadJ(r-i);
	  }
	
	Matrix<Real_wp> mat_elemR(r+1, r+1), mat_elem_diagJR(r+1, r+1);
	ComputeJointPanelIntegral(curved, Dsi, DsjR, PointsQuadI, PointsQuadJR,
				  NormaleQuadI, NormaleQuadJR, vec_u, vec_v, var, mat_elemR,
				  mat_elem_diagI, mat_elem_diagJR);
	
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= r; j++)
	    {
	      mat_elem(i, j) = mat_elemR(i, r-j);
	      mat_elem_diagJ(i, j) = mat_elem_diagJR(r-i, r-j);
	    }
	
	return;
      }
    
    if (var.GetDistance(Bi, Aj) <= R2::threshold )
      {
	// joint panels
	R2 vec_u = var.GetDiff(Bi, Ai); R2 vec_v = var.GetDiff(Aj, Bj);
	VectR2 PointsQuadIR(r+1), NormaleQuadIR(r+1);
	VectReal_wp DsiR(r+1);
	for (int i = 0; i <= r; i++)
	  {
	    DsiR(i) = Dsi(r-i);
	    PointsQuadIR(i) = PointsQuadI(r-i);
	    NormaleQuadIR(i) = NormaleQuadI(r-i);
	  }
	
	Matrix<Real_wp> mat_elemR(r+1, r+1), mat_elem_diagIR(r+1, r+1);
	ComputeJointPanelIntegral(curved, DsiR, Dsj, PointsQuadIR, PointsQuadJ,
				  NormaleQuadIR, NormaleQuadJ, vec_u, vec_v, var, mat_elemR,
				  mat_elem_diagIR, mat_elem_diagJ);
	
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= r; j++)
	    {
	      mat_elem(i, j) = mat_elemR(r-i, j);
	      mat_elem_diagI(i, j) = mat_elem_diagIR(r-i, r-j);
	    }
	
	return;	  
      }
    
    if (var.GetDistance(Bj, Aj) <= R2::threshold )
      {
	// joint panels
	R2 vec_u = var.GetDiff(Bi, Ai); R2 vec_v = var.GetDiff(Bj, Aj);
	VectR2 PointsQuadIR(r+1), NormaleQuadIR(r+1), PointsQuadJR(r+1), NormaleQuadJR(r+1);
	VectReal_wp DsiR(r+1), DsjR(r+1);
	for (int i = 0; i <= r; i++)
	  {
	    DsiR(i) = Dsi(r-i); DsjR(i) = Dsj(r-i);
	    PointsQuadIR(i) = PointsQuadI(r-i); PointsQuadJR(i) = PointsQuadJ(r-i);
	    NormaleQuadIR(i) = NormaleQuadI(r-i); NormaleQuadJR(i) = NormaleQuadJ(r-i);
	  }
	
	Matrix<Real_wp> mat_elemR(r+1, r+1), mat_elem_diagIR(r+1, r+1), mat_elem_diagJR(r+1, r+1);
	ComputeJointPanelIntegral(curved, DsiR, DsjR, PointsQuadIR, PointsQuadJR,
				  NormaleQuadIR, NormaleQuadJR, vec_u, vec_v, var, mat_elemR,
				  mat_elem_diagI, mat_elem_diagJ);
	
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= r; j++)
	    {
	      mat_elem(i, j) = mat_elemR(r-i, r-j);
	      mat_elem_diagI(i, j) = mat_elem_diagIR(r-i, r-j);
	      mat_elem_diagJ(i, j) = mat_elem_diagJR(r-i, r-j);
	    }
	
	return;
      }
    
    // evaluating if the two elements are close or far
    R2 center1 = Real_wp(0.5)*(Ai+Bi);
    R2 center2 = Real_wp(0.5)*(Aj+Bj);
    Real_wp dist = var.GetDistance(center1, center2);
    Real_wp radius = max(var.GetDistance(Ai, Bi), var.GetDistance(Aj, Bj));
    
    if (dist < 2.0*param_eta_quadrature*radius)
      {
	// close panels
	ComputeClosePanelIntegral(curved, Dsi, Dsj, PointsQuadI, PointsQuadJ, 
				  NormaleQuadI, NormaleQuadJ, var, mat_elem,
				  mat_elem_diagI, mat_elem_diagJ);	  
	
	return;
      }
    
    // far panels
    ComputeFarPanelIntegral(curved, Dsi, Dsj, PointsQuadI, PointsQuadJ, 
			    NormaleQuadI, NormaleQuadJ, var, mat_elem,
			    mat_elem_diagI, mat_elem_diagJ);      
  }
  
  
  //! compute the elementary mass matrix \int_Ai^Bi \varphi_i \varphi_j ds
  template<class T>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::ComputeMassMatrix(const VectR2& PtsAi, bool curved,
		      const VectReal_wp& Dsi, Matrix<Real_wp>& mat_elem)
  {
    //R2 Ai(PtsAi(0)), Bi(PtsAi(1));
    const ElementReference<Dimension1, 1>& Fb = *fe_class;
    mat_elem.Fill(0);
    int p = Fb.GetOrder();
    int r = Fb.GetNbPointsQuadratureInside();
    
    Real_wp val;
    for (int i = 0; i < r; i++)
      {
	val = Fb.Weights(i)*Dsi(i);
	for (int k = 0; k <= p; k++)
	  for (int n = 0; n <= p; n++)
	    mat_elem(k, n) += val*ValPhiFar(k, i)*ValPhiFar(n, i);
      }
  }  

  
  //! computes the elementary stiffness matrix \int_Ai^Bi 1/keps^2 d/ds(\varphi_i) d/ds(\varphi_j) ds
  template<class T> template<class MatrixJacob>
  void SingularDoubleQuadratureGalerkin_Base<T, Dimension1>
  ::ComputeStiffnessMatrix(const VectR2& PtsAi, bool curved,
			   const VectReal_wp& Dsi, const VectR2& NormaleQuad,
			   const MatrixJacob&, const Complex_wp& kinf, const Real_wp& epsilon,
			   Matrix<Complex_wp>& mat_elem)
  {
    //R2 Ai(PtsAi(0)), Bi(PtsAi(1));
    const ElementReference<Dimension1, 1>& Fb = *fe_class;
    mat_elem.Fill(0);
    int p = Fb.GetOrder();
    int r = Fb.GetNbPointsQuadratureInside();

    VectReal_wp KappaQuad(r);    
    if (epsilon < 0)
      ElementGeomReference<Dimension2>::GetCurvatureNodal(quad_basis, NormaleQuad, Dsi, KappaQuad);
    
    Complex_wp val;
    for (int i = 0; i < r; i++)
      {
	Real_wp eps = epsilon;
	if (epsilon < 0)
	  eps = 0.4*pow(abs(kinf), Real_wp(1)/3)*pow(abs(KappaQuad(i)), Real_wp(2)/3);
	
	Complex_wp keps = kinf + Iwp*eps;
	val = Fb.Weights(i)/(Dsi(i)*keps*keps);
	for (int k = 0; k <= p; k++)
	  for (int n = 0; n <= p; n++)
	    mat_elem(k, n) += val*GradPhiFar(k, i)*GradPhiFar(n, i);
      }    
  }
  
}

#define MONTJOIE_FILE_SINGULAR_INTEGRAL_EQUATION_2D_CXX
#endif
