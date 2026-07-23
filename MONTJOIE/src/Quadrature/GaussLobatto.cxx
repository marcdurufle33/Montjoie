#ifndef MONTJOIE_FILE_GAUSS_LOBATTO_CXX

namespace Montjoie
{
  
  /************
   * Globatto *
   ************/
  
  
  //! returns size of memory used by the object
  template<class T>
  size_t Globatto<T>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille += sizeof(T)*(points_.GetM()+weights_.GetM()+cte_phi.GetM());
    taille += sizeof(T)*val_grad.GetDataSize();
    return taille;
  }
  
  
  //! computation of 1-D integration formula
  /*!
    \param[in] ORDER order of the quadrature formula
    (in 1-D, there will be (ORDER+1) points)
    \param[in] type_quadrature Gauss or Gauss-Lobatto rule
   */
  template<class T>
  void Globatto<T>::ConstructQuadrature(int ORDER, int type_quadrature,
					T alpha, T beta)
  {
    order = ORDER;
    ndloc = (order+1);
    nb_points_quadrature = ndloc;
    
    switch (type_quadrature)
      {
      case QUADRATURE_GAUSS :
        ComputeGaussLegendre(points_, weights_, order);
        break;
      case QUADRATURE_LOBATTO :
        ComputeGaussLobatto(points_, weights_, order);
        break;
      case QUADRATURE_JACOBI :
        ComputeGaussJacobi(points_, weights_, order, alpha, beta);
        break;
      case QUADRATURE_RADAU :
        ComputeGaussRadauJacobi(points_, weights_, order, T(0), T(0), false);
        break;
      case QUADRATURE_RADAU_RIGHT :
        ComputeGaussRadauJacobi(points_, weights_, order, T(0), T(0), true);
        break;
      case QUADRATURE_RADAU_JACOBI :
        ComputeGaussRadauJacobi(points_, weights_, order, alpha, beta, false);
        break;
      case QUADRATURE_LOBATTO_JACOBI :
        ComputeGaussLobattoJacobi(points_, weights_, order, alpha, beta);
        break;
      case QUADRATURE_GAUSS_SQUARED :
        {
          // we replace (xi_i, omega_i) by (xi_i^2, 2 xi_i omega_i)
          // so that \int 1/sqrt(x) becomes a regular integral
          ComputeGaussLegendre(points_, weights_, order);
          for (int i = 0; i < points_.GetM(); i++)
            {
              T xi = points_(i);
              points_(i) = xi*xi;
              weights_(i) *= 2.0*xi;
            }
        }
	break;
      case QUADRATURE_GAUSS_BLENDED :
	if (blending_default == T(-100))
	  ComputeGaussBlendedFormulas(order+1, points_, weights_, T(order)/T(order+1));
	else
	  ComputeGaussBlendedFormulas(order+1, points_, weights_, blending_default);
	
	break;
      case TCHEBYSHEV :
        {
          points_.Reallocate(order+1);
          for (int i = 0; i <= order; i++)
            points_(order-i) = 0.5 + 0.5*cos((2.0*i+1.0)/(2.0*(order+1))*pi_wp);          
        } 
        break;
      default:
        cout << "Unknown quadrature rule " << type_quadrature << endl;
        abort();
      }
    
    ComputeFactorPhi();
  }
  
  
  // method to compute cte_phi from points
  template<class T>
  void Globatto<T>::ComputeFactorPhi()
  {
    cte_phi.Reallocate(ndloc);
    for (int i = 0; i <= order; i++)
      {
	cte_phi(i) = T(1);
	for (int j = 0; j <= order; j++)
	  if (j != i) 
	    cte_phi(i) *= points_(i) - points_(j);
	
	cte_phi(i) = T(1)/cte_phi(i);
      }
  }
  
  
  //! it is possible for the user to provide his own set of points
  /*!
    \param[in] Points1D provided set of points
    This method has to be used to provide 1-D basis functions
    integration weights are not modified
   */
  template<class T> template<class Vector1>
  void Globatto<T>::AffectPoints(const Vector1& points)
  {
    points_ = points;
    order = points.GetM() - 1;
    ndloc = order + 1;
    
    ComputeFactorPhi();
  }
  
  
  //! changes weights
  template<class T> template<class Vector1>
  void Globatto<T>::AffectWeights(const Vector1& poids)
  {
    weights_ = poids;
  }
  

  //! computes basis functions at position x
  template<class T>
  void Globatto<T>::ComputeValuesPhiRef(const T& x, Vector<T>& phi) const
  {
    phi.Reallocate(order+1);
    if (order <= 6)
      {
	for (int i = 0; i <= order; i++)
	  phi(i) = EvaluatePhi(i, x);

	return;
      }

    T denom(0), term;
    for (int i = 0; i <= order; i++)
      {
	if (abs(x-points_(i)) <= threshold)
	  {
	    phi.Zero();
	    phi(i) = 1.0;
            return;
	  }
	else
          {
            term = cte_phi(i)/(x - points_(i));
            phi(i) = term;
            denom += term;
          }
      }
    
    denom = 1.0/denom;
    for (int i = 0; i <= order; i++)
      phi(i) *= denom;
  }
  
  
  //! computes derivatives of basis functions at nodal points
  template<class T>
  void Globatto<T>::ComputeGradPhi(T eps)
  {
    T zero(0);
    val_grad.Reallocate(ndloc, ndloc);
    for (int i = 0; i <= order; i++)
      {
        T sum = 0;
        for (int j = 0; j <= order; j++)
          if (j != i)
            {
              val_grad(j, i) = cte_phi(j)/(cte_phi(i)*(points_(i)-points_(j)));  
              sum += val_grad(j, i);
            }
        
        val_grad(i, i) = -sum;
	if (abs(val_grad(i, i)) <= eps)
	  val_grad(i, i) = zero;
      }
  }
  
  
  //! displays informations about object Globatto
  template<class T>
  ostream& operator <<(ostream& out, const Globatto<T>& e)
  {
    out<<"Order of integration "<<e.GetOrder()<<endl;
    out<<"Number of integration points "<<e.GetNbPointsQuad()<<endl;
    
    out<<"Points of integration on edge "<<endl;
    for (int i = 0; i < e.ndloc; i++)
      out<<e.points_(i)<<" ";
    
    out<<endl<<"Weights of integration on edge"<<endl;
    for (int i = 0; i < e.ndloc; i++)
      out<<e.weights_(i)<<" ";
    
    return out;
  }
  

#ifdef MONTJOIE_WITH_MPFR
  template<class T>
  T Globatto<T>::blending_default(-100, MONTJOIE_DEFAULT_PRECISION);
#else
  template<class T>
  T Globatto<T>::blending_default(-100);
#endif
  
  
  /******************
   * SubdivGlobatto *
   ******************/
  
  
  //! returns interpolation point i
  const Real_wp& SubdivGlobatto::Points(int i) const
  {
    int offset = 0;
    for (int k = 0; k < points.GetM(); k++)
      {
	if (i < offset+points(k).GetM())
	  return points(k)(i-offset);
	
	offset += points(k).GetM()-1;
      }
    
    return points(0)(0);
  }


  //! returns all interpolation points
  const Vector<Real_wp>& SubdivGlobatto::Points() const
  {
    return points1d;
  }


  //! initialisation of shape functions
  /*!
    \param[in] regular if true regular points are used instead of Gauss-Lobatto
    \param[in] nb_subdiv number of subdivisions n
    \param[in] r degree of Lagrange polynomials on each interval
   */
  void SubdivGlobatto::Init(bool regular, int nb_subdiv, int r, int type_quad)
  {
    VectReal_wp x(r+1);
    if (!regular)
      {
        Globatto<Real_wp> lob;
        lob.ConstructQuadrature(r, type_quad);
        x = lob.Points();
      }
    else
      {
	x.Fill();
	Mlt(Real_wp(1)/r, x);
      }
    
    InitPoints(nb_subdiv, x);
  }

  
  //! initializations of shape functions
  void SubdivGlobatto::InitPoints(int nb_subdiv, const VectReal_wp& x)
  {
    int r = x.GetM()-1;
    VectReal_wp cte;
    xpos.Reallocate(nb_subdiv+1);
    xpos.Fill(); Mlt(Real_wp(1)/nb_subdiv, xpos);
    xpos(0) = 0.0; xpos(nb_subdiv) = 1.0;
        
    cte.Reallocate(r+1); cte.Fill(1);
    for (int i = 0; i <= r; i++)
      {
	for (int j = 0; j <= r; j++)
	  if (i != j)
	    cte(i) *= x(i) - x(j);
	
	cte(i) = 1.0/cte(i)*pow(Real_wp(nb_subdiv), Real_wp(r));
      }
    
    CteG.Reallocate(nb_subdiv);
    points.Reallocate(nb_subdiv);
    points1d.Reallocate(nb_subdiv*r + 1);
    for (int i = 0; i < nb_subdiv; i++)
      {
	CteG(i) = cte;
	points(i).Reallocate(r+1);
	for (int j = 0; j <= r; j++)
          points(i)(j) = xpos(i)*(1.0-x(j)) + xpos(i+1)*x(j);
        
        for (int j = 0; j <= r; j++)
          points1d(i*r + j) = points(i)(j);
      }
  }
  
  
  //! returns the number of interpolation points - 1
  int SubdivGlobatto::GetOrder() const
  {
    int order = 0;
    for (int i = 0; i < points.GetM(); i++)
      order += points(i).GetM()-1;
    
    return order;
  }
  

  //! returns the number of interpolation points - 1
  int SubdivGlobatto::GetGeometryOrder() const
  {
    return GetOrder();
  }
  
  
  //! evalues shape function at a given point of [0, 1]
  Real_wp SubdivGlobatto::EvaluatePhi(int node, const Real_wp& x) const
  {
    if (x < 0)
      return (node == 0);
    
    int offset = 0;
    for (int k = 0; k < points.GetM(); k++)
      {
	if (x <= xpos(k+1))
	  {
	    int i = node - offset;
	    if ((i < 0)||(i >= points(k).GetM()))
	      return Real_wp(0);
	    else
	      {
		Real_wp phi = CteG(k)(i);
		for (int j = 0; j < points(k).GetM(); j++)
		  if (j != i)
		    phi *= x - points(k)(j);
	    
		return phi;
	      }
	  }
	
	offset += points(k).GetM() - 1;
      }
    
    return (node==offset);
  }  
    
} // end namespace

#define MONTJOIE_FILE_GAUSS_LOBATTO_CXX
#endif
