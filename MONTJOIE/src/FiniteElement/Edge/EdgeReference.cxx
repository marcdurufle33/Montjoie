#ifndef MONTJOIE_FILE_EDGE_REFERENCE_CXX

namespace Montjoie
{
  /*****************
   * EdgeReference *
   *****************/
  
  
  //! default constructor
  ElementGeomReference<Dimension1>::ElementGeomReference()
  {
    print_level = 0;
    order = 0;
    order_geom = 0;
    order_quad = 0;
    nb_points = 0;
    nb_dof_loc = 0;
  }
  
  
  //! returns the size of memory used by the object
  size_t ElementGeomReference<Dimension1>::GetMemorySize() const
  {
    size_t taille = sizeof(Real_wp)*(points_dof.GetM()+Points.GetM()+Weights.GetM());
    taille += stiffness_matrix.GetMemorySize() + gradient_matrix.GetMemorySize();
    taille += mass_matrix.GetMemorySize() + sizeof(*this);
    return taille;
  }
  
  
  //! constructs quadrature formulas in 1-D
  /*!
    \param order_ : order of approximation
  */
  void ElementGeomReference<Dimension1>::ConstructQuadrature(int r, int type_quadrature)
  {
    order_quad = r;
    nb_points = order_quad + 1;
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(order_quad, type_quadrature);
    
    Weights = gauss.Weights();
    Points = gauss.Points();
  }
  
  
  void ElementGeomReference<Dimension1>
  ::GetPolynomialFunctions(Vector<UnivariatePolynomial<Real_wp> >& Phi,
			   Vector<UnivariatePolynomial<Real_wp> >& dPhi) const
  {
    cout << "This function is not implemented for this finite element" << endl;
    abort(); 
  }
  

  //! computation of F_i
  /*!
    \param s : list of vertices of the edge
    \param points : list of quadrature points
   */
  void ElementGeomReference<Dimension1>::FjElem(const VectReal_wp& s, VectReal_wp& points) const
  {
    int nb_quad = Points.GetM();
    points.Reallocate(nb_quad);
    for (int i = 0; i < nb_quad; i++)
      points(i) = Points(i) * s(1) + (1.0-Points(i) ) * s(0) ;
    
  }


  //! computation of F_i
  /*!
    \param s : list of vertices of the edge
    \param points : list of dof points
   */
  void ElementGeomReference<Dimension1>::FjElemDof(const VectReal_wp& s, VectReal_wp& points) const
  {
    int nb_dof = points_dof.GetM();
    points.Reallocate(nb_dof);
    for (int i = 0; i < nb_dof; i++)
      points(i) = points_dof(i) * s(1) + (1.0-points_dof(i) ) * s(0) ;
    
  }

  
  //! displays information about the object ElementGeomReference<Dimension1>
  ostream& operator<<(ostream& out, const ElementGeomReference<Dimension1>& e)
  {
    out<<"Points of integration "<<endl;
    out<<e.Points<<endl;
    out<<"Weights of integration "<<endl;
    out<<e.Weights<<endl;
    return out;
  }


  /*************
   * EdgeGauss *
   *************/
  
  
  //! default constructor
  EdgeGauss::EdgeGauss() : ElementGeomReference<Dimension1>()
  {
    mass_lumping = false;
  }
  

  //! returns the size of memory used by the object
  size_t EdgeGauss::GetMemorySize() const
  {
    size_t taille = ElementGeomReference<Dimension1>::GetMemorySize();
    taille += sizeof(Real_wp)*(Value_Phi.GetDataSize()+Gradient_Phi.GetDataSize());
    taille += lob_geom.GetMemorySize()+lob_basis.GetMemorySize();
    taille += sizeof(*this) - sizeof(ElementGeomReference<Dimension1>);
    return taille;
  }
  
  
  //! how to number the 1D-mesh
  void EdgeGauss::ConstructNumberMap(NumberMap& nmap, int dg_form) const
  {
    if (dg_form == 1)
      {
	nmap.SetNbDofVertex(order, 0);
	nmap.SetNbDofEdge(order, order+1);
	return;
      }

    nmap.SetNbDofVertex(order, 1);
    nmap.SetNbDofEdge(order, order-1);
  }
  
  
  //! construction of finite element for order equal to r
  void EdgeGauss::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					 int type_func)
  {
    if (rgeom == 0)
      rgeom = r;
    
    if (rquad == 0)
      rquad = r;
    
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_GAUSS;
    
    ConstructQuadrature(rquad, type_quad);
    
    order_geom = rgeom;
    order = r;
    lob_geom.ConstructQuadrature(rgeom, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    nb_dof_loc = order+1;
    ConstructFunctions(type_func);
    ConstructStiffnessMatrix();
  }


  //! changes points for dofs
  void EdgeGauss::SetDofPoints(const VectReal_wp& pts)
  {
    lob_basis.AffectPoints(pts);
    
    order = pts.GetM()-1;
    nb_dof_loc = order+1;
    points_dof = pts;
    ConstructStiffnessMatrix();
  }
  
  
  //! construction of basis functions
  void EdgeGauss::ConstructFunctions(int type_func)
  {
    if (type_func == -1)
      type_func = LOBATTO;

    if (type_func == LOBATTO)
      lob_basis.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    else if (type_func == RADAU)
      lob_basis.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_RADAU_RIGHT);
    else if (type_func == GAUSS)
      lob_basis.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_GAUSS);
    else if (type_func == LOBATTO_INT)
      {
	VectReal_wp points_dof1d, points, weights;
	ComputeGaussLobatto(points, weights, order+2);
	points_dof1d.Reallocate(order+1);
	for (int i = 0; i <= order; i++)
	  points_dof1d(i) = points(i+1);
	
	lob_basis.AffectPoints(points_dof1d);    
      }
    
    points_dof = lob_basis.Points();
  }
  
  
  //! computation of stiffness matrix
  void EdgeGauss::ConstructStiffnessMatrix()
  {
    mass_matrix.Reallocate(nb_dof_loc, nb_dof_loc);
    stiffness_matrix.Reallocate(nb_dof_loc, nb_dof_loc);
    gradient_matrix.Reallocate(nb_dof_loc, nb_dof_loc);
    Gradient_Phi.Reallocate(nb_dof_loc, nb_points);
    Value_Phi.Reallocate(nb_dof_loc, nb_points);
    for (int i = 0; i < nb_dof_loc; i++)
      {
	for (int j = 0; j < nb_points; j++)
	  {
	    Value_Phi(i, j) = lob_basis.EvaluatePhi(i, Points(j));
	    Gradient_Phi(i, j) = lob_basis.EvaluatePhiGrad(i, Points(j));
	  }
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	{
	  Real_wp val = 0.0;
	  for (int k = 0; k < nb_points; k++)
	    val += Weights(k)*Value_Phi(i, k)*Value_Phi(j, k);
	  
	  mass_matrix(i, j) = val;

	  val = 0.0;
	  for (int k = 0; k < nb_points; k++)
	    val += Weights(k)*Gradient_Phi(i, k)*Gradient_Phi(j, k);
	  
	  stiffness_matrix(i, j) = val;
	  
	  val = 0.0;
	  for (int k = 0; k < nb_points; k++)
	    val += Weights(k)*Gradient_Phi(j, k)*Value_Phi(i, k);
	  
	  gradient_matrix(i, j) = val;
	}
  }
  
  
  //! compute res(i) = \f$ \hat{\varphi}_i(\mbox{pointloc}) \f$
  void EdgeGauss::ComputeValuesPhiRef(const Real_wp& point_loc, VectReal_wp& res) const
  {
    if (this->order >= 8)
      lob_basis.ComputeValuesPhiRef(point_loc, res);
    else
      {
        res.Reallocate(nb_dof_loc);
        for (int i = 0; i < nb_dof_loc; i++)
          res(i) = lob_basis.EvaluatePhi(i, point_loc);
      }
  }
  

  //! compute res(i) = \f$ \hat{\varphi}'_i(\mbox{pointloc}) \f$
  void EdgeGauss::ComputeGradientPhiRef(const Real_wp& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      res(i) = lob_basis.EvaluatePhiGrad(i, point_loc);
  }

  
  //! compute res(i) = \f$ \hat{\varphi}_i(\mbox{pointloc}) \f$
  void EdgeGauss::ComputeValuesPhiNodalRef(const Real_wp& point_loc, VectReal_wp& res) const
  {
    if (this->order >= 8)
      lob_geom.ComputeValuesPhiRef(point_loc, res);
    else
      {
        res.Reallocate(lob_geom.GetOrder()+1);
        for (int i = 0; i <= lob_geom.GetOrder(); i++)
          res(i) = lob_geom.EvaluatePhi(i, point_loc);
      }
  }

  
  void EdgeGauss::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_points);
    for (int j = 0; j < nb_points; j++)
      phi(j) = Value_Phi(i, j);
  }

  
  //! returns varphi_i(point_loc)
  Real_wp EdgeGauss::GetValuePhi1D(int i, const Real_wp& point_loc) const
  {
    return lob_basis.EvaluatePhi(i, point_loc);
  }
  

  //! returns varphi_i'(point_loc)
  Real_wp EdgeGauss::GetGradientPhi1D(int i, const Real_wp& point_loc) const
  {
    return lob_basis.EvaluatePhiGrad(i, point_loc);
  }

  
  //! fills polynomials associated with basis functions
  void EdgeGauss::GetPolynomialFunctions(Vector<UnivariatePolynomial<Real_wp> >& Phi,
					 Vector<UnivariatePolynomial<Real_wp> >& dPhi) const
  {
    Real_wp one(1), coef;
    
    // Lagrange polynomials with points of lob_basis
    int r = lob_basis.GetOrder();
    Phi.Reallocate(r+1);
    UnivariatePolynomial<Real_wp> monome;
    monome.SetOrder(1);
    for (int i = 0; i <= r; i++)
      {
	Phi(i).SetOrder(0); Phi(i)(0) = one;
	for (int j = 0; j <= r; j++)
	  if (i != j)
	    {
	      coef = one/(lob_basis.Points(i) - lob_basis.Points(j));
	      monome(0) = -lob_basis.Points(j)*coef;
	      monome(1) = coef;
	      Phi(i) = Phi(i)*monome;
	    }
      }
    
    // derivatives of these functions
    dPhi.Reallocate(r+1);
    for (int i = 0; i <= r; i++)
      DerivatePolynomial(Phi(i), dPhi(i));
  }

  
  //! compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i(x) dx \f$
  /*!
    \param[in] feval feval(j) is the evaluation of f on the quadrature point j
    \param[out] contrib result vector
  */
  template<class Vector1>
  void EdgeGauss::ComputeIntegralGen(const Vector1& feval, Vector1& contrib) const
  {
    contrib.Reallocate(this->nb_dof_loc);
    contrib.Fill(0);
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = 0; j < nb_points; j++)
	contrib(i) += feval(j) * Value_Phi(i, j) * Weights(j);
  }


  template<class Vector1>
  void EdgeGauss::ApplyChGen(const Vector1& feval, Vector1& contrib) const
  {
    contrib.Reallocate(this->nb_dof_loc);
    contrib.Zero();
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = 0; j < nb_points; j++)
	contrib(i) += feval(j) * Value_Phi(i, j);
  }


  template<class Vector1>
  void EdgeGauss::ApplyChTransposeGen(const Vector1& feval, Vector1& contrib) const
  {
    contrib.Reallocate(nb_points);
    contrib.Zero();
    for (int j = 0; j < nb_points; j++)
      for (int i = 0; i < this->nb_dof_loc; i++)
	contrib(j) += feval(i) * Value_Phi(i, j);
  }

  
  //! compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i'(x) dx \f$
  /*! 
    \param[in] feval feval(j) is the evalution of f on the quadrature point j
    \param[out] contrib result vector
  */
  template<class Vector1>
  void EdgeGauss::ComputeIntegralGradientGen(const Vector1& feval,Vector1& contrib) const
  {
    contrib.Reallocate(this->nb_dof_loc);
    contrib.Fill(0);
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = 0; j < nb_points; j++)
	contrib(i) += feval(j) * Gradient_Phi(i, j) * Weights(j);
  }


  void EdgeGauss::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    contrib = feval;
  }
  
  
  void EdgeGauss::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    contrib = feval;
  }
  
  
  //! computes A = A + mass \int \varphi_i \varphi_j dx + C \int d/dx(\varphi_i) d/dx(\varphi_j) dx
  //! + D \int d/dx(\varphi_i) \varphi_j dx + E \int d/dx(\varphi_j) \varphi_i dx
  /*!
    \param[in] m offset for rows
    \param[in] n offset for columns
    \param[in] mass, C, D, E coefficients
    \param[inout] A matrix to which elementary matrix is added
    This function performs the operation :
    A(m:m+N, n:n+N) += mass * M + C * S + D * R + E * Rt
    where N is the size of the stiffness matrix S
    M is the usual mass matrix M_{i, j} = \int \varphi_i \varphi_j
    S is the usual stiffness matrix S_{i, j} = \int d/dx(\varphi_i) d/dx(\varphi_j) dx
    R the gradient matrix R_{i, j} = \int d/dx(varphi_j) \varphi_i dx
    Rt is the transpose of R
   */
  template<class T0, class Matrix1>
  void EdgeGauss::AddConstantElemMatrixGen(int m, int n, const T0& mass, const T0& C,
					   const T0& D, const T0& E,
					   const TinyVector<bool, 4>& null_term, Matrix1& A) const
  {
    T0 vloc;
    bool sym = A.IsSymmetric();
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	if ((!sym) || (m+i <= n+j))
	  {
	    vloc = mass_matrix(i, j)*mass + stiffness_matrix(i, j)*C
	      + gradient_matrix(j, i)*D + gradient_matrix(i, j)*E;
	    
	    A.AddInteraction(m+i, n+j, vloc);
	  }
  }
  

  //! computes A = A + mass \int \varphi_i \varphi_j dx + C \int d/dx(\varphi_i) d/dx(\varphi_j) dx
  //! + D \int d/dx(\varphi_i) \varphi_j dx + E \int d/dx(\varphi_j) \varphi_i dx
  /*!
    \param[in] m offset for rows
    \param[in] n offset for columns
    \param[in] mass, C, D, E coefficients/h*omega_i on quadrature points
    \param[inout] A matrix to which elementary matrix is added
    This function performs the operation :
    A(m:m+N, n:n+N) += M + S + R + Rt
    where N is the size of the stiffness matrix S
    M is the usual mass matrix M_{i, j} = \int mass \varphi_i \varphi_j
    S is the usual stiffness matrix S_{i, j} = \int S d/dx(\varphi_i) d/dx(\varphi_j) dx
    R the gradient matrix R_{i, j} = \int E d/dx(varphi_j) \varphi_i dx
    Rt is given as Rt_{i, j} = \int D d/dx(varphi_i) \varphi_j dx
   */
  template<class T, class Matrix1>
  void EdgeGauss::AddVariableElemMatrixGen(int m, int n, const Vector<T>& mass,
					   const Vector<T>& C, const Vector<T>& D, const Vector<T>& E,
					   const TinyVector<bool, 4>& null_term, Matrix1& A) const
  {
    T val;
    bool sym = A.IsSymmetric();
    for (int i = 0; i < nb_dof_loc; i++)
      {
        for (int j = 0; j < nb_dof_loc; j++)
          {
	    if ((!sym) || (m+i <= n+j))
	      {
		SetComplexZero(val);
		if (!null_term(0))
		  for (int k = 0; k <= order_quad; k++)
		    val += mass(k)*Value_Phi(i, k)*Value_Phi(j, k);
		
		if (!null_term(1))
		  for (int k = 0; k <= order_quad; k++)
		    val += C(k)*Gradient_Phi(i, k)*Gradient_Phi(j, k);

		if (!null_term(2))
		  for (int k = 0; k <= order_quad; k++)
		    val += D(k)*Value_Phi(j, k)*Gradient_Phi(i, k);
		
		if (!null_term(3))
		  for (int k = 0; k <= order_quad; k++)
		    val += E(k)*Value_Phi(i, k)*Gradient_Phi(j, k);
				
		A.AddInteraction(m+i, n+j, val);
	      }
          }
      }
  }
  

  void EdgeGauss::ComputeIntegralRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeIntegralGen(feval, contrib);
  }

  void EdgeGauss::ComputeIntegralRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeIntegralGen(feval, contrib);
  }

  void EdgeGauss::ApplyCh(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ApplyChGen(feval, contrib);
  }

  void EdgeGauss::ApplyCh(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ApplyChGen(feval, contrib);
  }

  void EdgeGauss::ApplyChTranspose(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ApplyChTransposeGen(feval, contrib);
  }

  void EdgeGauss::ApplyChTranspose(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ApplyChTransposeGen(feval, contrib);
  }

  void EdgeGauss::ComputeIntegralGradientRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeIntegralGradientGen(feval, contrib);
  }

  void EdgeGauss::ComputeIntegralGradientRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeIntegralGradientGen(feval, contrib);
  }

  void EdgeGauss::AddConstantElemMatrix(int m, int n, const Real_wp& mass, const Real_wp& C,
					const Real_wp& D, const Real_wp& E,
					const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }
  
  void EdgeGauss::AddConstantElemMatrix(int m, int n, const Complex_wp& mass, const Complex_wp& C,
					const Complex_wp& D, const Complex_wp& E,
					const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }
        
  void EdgeGauss::AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
					const VectReal_wp& C, const VectReal_wp& D, const VectReal_wp& E,
					const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& mat) const
  {
    AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, mat);
  }
  
  void EdgeGauss::AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
					const VectComplex_wp& C, const VectComplex_wp& D, const VectComplex_wp& E,
					const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& mat) const
  {
    AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, mat);
  }
  
  
  //! displays information about object EdgeGauss
  ostream& operator<<(ostream& out, const EdgeGauss& e)
  {
    out<<"Stiffness Matrix "<<endl;
    out<<e.stiffness_matrix<<endl;
    return out;
  }
  
  
  /***************
   * EdgeLobatto *
   ***************/
  
  
  //! default constructor
  EdgeLobatto::EdgeLobatto() : EdgeGauss()
  {
    mass_lumping = true;
  }
  
  
  //! construction of finite element for order equal to r
  void EdgeLobatto::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad, int type_func)
  {
    rgeom = r;
    if (rquad == 0)
      rquad = r;
    
    if (type_quad == -1)
      type_quad = lob_geom.QUADRATURE_LOBATTO;
    
    ConstructQuadrature(rquad, type_quad);

    type_func = EdgeGauss::LOBATTO;
    if (type_quad == Globatto<Real_wp>::QUADRATURE_RADAU_RIGHT)
      type_func = RADAU;
    else if (type_quad != Globatto<Real_wp>::QUADRATURE_LOBATTO)
      {
	cout << "this element does not authorize specific quadrature rule " << endl;
        abort();
      }
    
    order_geom = rgeom;
    order = r;
    nb_dof_loc = order+1;
    ConstructFunctions(type_func);
    ConstructStiffnessMatrix();
  }
  
  
  //! compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i(x) dx \f$
  /*!
    \param[in] feval feval(j) is the evaluation of f on the quadrature point j
    \param[out] contrib result vector
  */
  template<class Vector1>
  void EdgeLobatto::ComputeIntegralGen(const Vector1& feval, Vector1& contrib) const
  {
    if (order != order_quad)
      return EdgeGauss::ComputeIntegralGen(feval, contrib);

    contrib.Reallocate(nb_dof_loc);
    contrib.Fill(0);
    for (int i = 0; i < nb_dof_loc; i++)
      contrib(i) = feval(i) * Weights(i);
  }


  template<class Vector1>
  void EdgeLobatto::ApplyChGen(const Vector1& feval, Vector1& contrib) const
  {
    if (order != order_quad)
      return EdgeGauss::ApplyChGen(feval, contrib);

    contrib = feval;
  }


  void EdgeLobatto::ComputeIntegralRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeIntegralGen(feval, contrib);
  }

  void EdgeLobatto::ComputeIntegralRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeIntegralGen(feval, contrib);
  }

  void EdgeLobatto::ApplyCh(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ApplyChGen(feval, contrib);
  }

  void EdgeLobatto::ApplyCh(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ApplyChGen(feval, contrib);
  }

  void EdgeLobatto::ApplyChTranspose(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ApplyChGen(feval, contrib);
  }

  void EdgeLobatto::ApplyChTranspose(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ApplyChGen(feval, contrib);
  }

  void EdgeLobatto::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    if (order != order_quad)
      return EdgeGauss::GetValueSinglePhiQuadrature(i, phi);

    phi.Reallocate(nb_points);
    phi.Zero(); phi(i) = Real_wp(1);
  }


  EdgeRadau::EdgeRadau()
  {
  }

  
  void EdgeRadau::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					 int type_func)
  {
    if (type_quad == -1)
      type_quad = lob_geom.QUADRATURE_RADAU_RIGHT;

    EdgeLobatto::ConstructFiniteElement(r, rgeom, rquad, type_quad, type_func);
  }


  /******************
   * EdgeHierarchic *
   ******************/
  
  
  //! default constructor
  EdgeHierarchic::EdgeHierarchic() : ElementGeomReference<Dimension1>(), Leg_pol(1)
  {
    this->mass_lumping = false;
    type_function = JACOBI_11;
  }
  
  
  //! how to number the 1D-mesh
  void EdgeHierarchic::ConstructNumberMap(NumberMap& nmap, int dg_form) const
  {
    if (dg_form == 1)
      {
	nmap.SetNbDofVertex(order, 0);
	nmap.SetNbDofEdge(order, order+1);
	return;
      }
    
    nmap.SetNbDofVertex(order, 1);
    nmap.SetNbDofEdge(order, order-1);
  }
  

  //! returns the size of memory used by the object
  size_t EdgeHierarchic::GetMemorySize() const
  {
    size_t taille = ElementGeomReference<Dimension1>::GetMemorySize();
    taille += sizeof(Real_wp)*(Value_Phi.GetDataSize()+Gradient_Phi.GetDataSize());
    taille += jacobi_11_pol.GetMemorySize();
    taille += CoefLeg11.GetMemorySize() + legendre_pol.GetMemorySize();
    taille += sizeof(*this) - sizeof(ElementGeomReference<Dimension1>);
    return taille;
  }


  //! construction of finite element for order equal to r  
  void EdgeHierarchic::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad, int type_func)
  {
    if (rgeom == 0)
      rgeom = r;
    
    if (rquad == 0)
      rquad = r;
    
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_GAUSS;
    
    if (type_func == -1)
      type_func = JACOBI_11;

    ConstructQuadrature(rquad, type_quad);
    
    order_geom = rgeom;
    order = r;
    nb_dof_loc = order+1;
    ConstructFunctions(type_func);
  }


  void EdgeHierarchic::SetDofPoints(const VectReal_wp& pts)
  {
    cout << "Degrees of freedom are not associated with points" << endl;
    abort();
  }
  
  
  //! computation of stiffness matrix
  void EdgeHierarchic::ConstructFunctions(int type_func)
  {
    jacobi_11_pol.Clear();
    legendre_pol.Clear();
    type_function = type_func;
    if ((type_func == LEGENDRE) || (type_func == LEGENDRE_COMBINED))
      {
	GetJacobiPolynomial(legendre_pol, order+1, Real_wp(0), Real_wp(0));	
      }
    else if (type_func == SHEN)
      {
	// basis similar to P_m^{1,1} with a different normalization
	Leg_pol = LegendrePolynomial<Real_wp>(order+1);
	CoefLeg11.Reallocate(order+1);
	CoefLeg11.Fill(1.0);
	//for (int k = 1; k < order; k++)
	//CoefLeg11(k+1) = 1.0/sqrt(2.0*(4.0*k+2));

	//DISP(CoefLeg11);
      }
    else
      {
	// computing P_m^{1,1}
	GetJacobiPolynomial(jacobi_11_pol, order, Real_wp(1), Real_wp(1));
	
	// orthonormalisation coefficients
	CoefLeg11.Reallocate(order-1); CoefLeg11.Fill(0);
	VectReal_wp Pn;
	for (int i = 0; i <= order_quad; i++)
	  {
	    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*Points(i) - 1.0, Pn);
	    for (int j = 1; j < order; j++)
	      CoefLeg11(j-1) += Weights(i)*square(Points(i)*(1.0-Points(i))*Pn(j-1));
	  }
	
	for (int j = 1; j < order; j++)
	  CoefLeg11(j-1) = 1.0/sqrt(CoefLeg11(j-1));
      }
    
    // stiffness and mass matrix
    stiffness_matrix.Reallocate(nb_dof_loc, nb_dof_loc);
    gradient_matrix.Reallocate(nb_dof_loc, nb_dof_loc);
    mass_matrix.Reallocate(nb_dof_loc, nb_dof_loc);
    Gradient_Phi.Reallocate(nb_dof_loc, nb_points);
    Value_Phi.Reallocate(nb_dof_loc, nb_points);
    VectReal_wp phi, grad_phi;
    for (int j = 0; j < nb_points; j++)
      {
        ComputeValuesPhiRef(Points(j), phi);
        ComputeGradientPhiRef(Points(j), grad_phi);
        for (int i = 0; i < nb_dof_loc; i++)
          {
	    Value_Phi(i, j) = phi(i);
	    Gradient_Phi(i, j) = grad_phi(i);
	  }
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	{
	  Real_wp val = 0.0;
	  for (int k = 0; k < nb_points; k++)
	    val += Weights(k)*Gradient_Phi(i, k)*Gradient_Phi(j, k);
	  
	  stiffness_matrix(i, j) = val;

	  val = 0.0;
	  for (int k = 0; k < nb_points; k++)
	    val += Weights(k)*Gradient_Phi(j, k)*Value_Phi(i, k);
	  
	  gradient_matrix(i, j) = val;

	  val = 0.0;
	  for (int k = 0; k < nb_points; k++)
	    val += Weights(k)*Value_Phi(i, k)*Value_Phi(j, k);
	  
	  mass_matrix(i, j) = val;
	}

    //DISP(stiffness_matrix);
    //DISP(gradient_matrix);
    //DISP(mass_matrix);
    
    GetCholesky(mass_matrix);
    
    // quadrature points are used to project on basis functions
    this->points_dof = this->Points;
  }
  
  
  //! compute res(i) = \f$ \hat{\varphi}_i(\mbox{pointloc}) \f$
  void EdgeHierarchic::ComputeValuesPhiRef(const Real_wp& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    if (type_function == SHEN)
      {
	VectReal_wp Pn;
	Leg_pol.EvaluatePn(order+1, 2.0*point_loc-1.0, Pn);
	
	res(0) = 1.0 - point_loc;
	res(1) = point_loc;
	
	for (int j = 1; j < order; j++)
	  res(j+1) = CoefLeg11(j+1)*(Pn(j-1) - Pn(j+1));
      }    
    else if (legendre_pol.GetM() > 0)
      {
	res.Fill(0);
	EvaluateJacobiPolynomial(legendre_pol, order, 2.0*point_loc-1.0, res);
	
	if (type_function == LEGENDRE_COMBINED)
	  {
	    VectReal_wp Pn(res);
	    res(0) = 1.0;
	    if (nb_dof_loc >= 2)
	      res(1) = 2.0*point_loc-1.0;
	    
	    Real_wp coef_phi0, coef_phi1;
	    for (int k = 2; k <= order; k++)
	      {
		coef_phi0 = -(k-1)*Pn(k-2)/k ;
		coef_phi1 = (2*k-1)*Pn(k-1)/k;
		res(k) = coef_phi0 + coef_phi1*res(1);
	      }
	  }
      }
    else
      {
	res(0) = 1.0 - point_loc;
	res(order) = point_loc;
	
	if (order >= 2)
	  {
	    VectReal_wp Pn;
	    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*point_loc - 1.0, Pn);
	    Real_wp vloc = (1.0-point_loc)*point_loc;
	    for (int j = 1; j < order; j++)
	      res(j) = vloc*Pn(j-1)*CoefLeg11(j-1);
	  }
      }
  }
  
  
  //! compute res(i) = \f$ \hat{\varphi}_i(\mbox{pointloc}) \f$
  void EdgeHierarchic::ComputeValuesPhiNodalRef(const Real_wp& point_loc, VectReal_wp& res) const
  {
    abort();
  }
  

  void EdgeHierarchic::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_points);
    for (int j = 0; j < nb_points; j++)
      phi(j) = Value_Phi(i, j);
  }

  
  //! computation of derivative of basis functions on point_loc
  void EdgeHierarchic::ComputeGradientPhiRef(const Real_wp& point_loc, VectReal_wp& res) const
  {  
    res.Reallocate(nb_dof_loc);
    
    if 	(type_function == SHEN)
      {
	VectReal_wp Pn;
	Leg_pol.EvaluatePn(order+1, 2.0*point_loc-1.0, Pn);
	
	res(0) = -1.0;
	res(1) = 1.0;
	
	for (int j = 1; j < order; j++)
	  res(j+1) = -2.0*CoefLeg11(j+1)*(2.0*j+1)*Pn(j);	
      }
    else if (legendre_pol.GetM() > 0)
      {
	VectReal_wp Pn;

	res.Fill(0);
	EvaluateJacobiPolynomial(legendre_pol, order, 2.0*point_loc-1.0, Pn, res);
	Mlt(Real_wp(2), res);

	if (type_function == LEGENDRE_COMBINED)
	  {
	    // To be done ...
	  }
      }
    else
      {
	res(0) = -1.0;
	res(order) = 1.0;
	
	if (order >= 2)
	  {
	    VectReal_wp Pn, dPn;
	    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*point_loc - 1.0, Pn, dPn);
	    Real_wp vloc = (1.0-point_loc)*point_loc;
	    Real_wp dvloc = 1.0 - 2.0*point_loc;
	    for (int j = 1; j < order; j++)
	      res(j) = (2.0*vloc*dPn(j-1) + dvloc*Pn(j-1))*CoefLeg11(j-1);
	  }
      }
  }
  
  
  //! compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i(x) dx \f$
  /*!
    \param[in] feval feval(j) is the evaluation of f on the quadrature point j
    \param[out] contrib result vector
  */
  template<class Vector1>
  void EdgeHierarchic::ComputeIntegralGen(const Vector1& feval, Vector1& contrib) const
  {
    contrib.Reallocate(this->nb_dof_loc);
    contrib.Fill(0);
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = 0; j < nb_points; j++)
	contrib(i) += feval(j) * Value_Phi(i, j) * Weights(j);
  }


  template<class Vector1>
  void EdgeHierarchic::ApplyChGen(const Vector1& feval, Vector1& contrib) const
  {
    contrib.Reallocate(this->nb_dof_loc);
    contrib.Zero();
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = 0; j < nb_points; j++)
	contrib(i) += feval(j) * Value_Phi(i, j);
  }


  template<class Vector1>
  void EdgeHierarchic::ApplyChTransposeGen(const Vector1& feval, Vector1& contrib) const
  {
    contrib.Reallocate(nb_points);
    contrib.Zero();
    for (int j = 0; j < nb_points; j++)
      for (int i = 0; i < this->nb_dof_loc; i++)
	contrib(j) += feval(i) * Value_Phi(i, j);
  }

  
  //! compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i'(x) dx \f$
  /*! 
    \param[in] feval feval(j) is the evalution of f on the quadrature point j
    \param[out] contrib result vector
  */
  template<class Vector1>
  void EdgeHierarchic::ComputeIntegralGradientGen(const Vector1& feval,Vector1& contrib) const
  {
    contrib.Reallocate(this->nb_dof_loc);
    contrib.Fill(0);
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = 0; j < nb_points; j++)
	contrib(i) += feval(j) * Gradient_Phi(i, j) * Weights(j);
  }
  
  
  //! x is overwritten by M^{-1} x where M is the mass matrix
  template<class Vector1>
  void EdgeHierarchic::SolveMassMatrix(Vector1& x) const
  {
    Seldon::SolveCholesky(SeldonNoTrans, mass_matrix, x);
    Seldon::SolveCholesky(SeldonTrans, mass_matrix, x);
  }
  
  
  void EdgeHierarchic::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeIntegralRef(feval, contrib);
    SolveMassMatrix(contrib);
  }
  
  
  void EdgeHierarchic::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeIntegralRef(feval, contrib);
    SolveMassMatrix(contrib);
  }
  

  void EdgeHierarchic::ComputeIntegralRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeIntegralGen(feval, contrib);
  }

  void EdgeHierarchic::ComputeIntegralRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeIntegralGen(feval, contrib);
  }


  void EdgeHierarchic::ApplyCh(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ApplyChGen(feval, contrib);
  }

  void EdgeHierarchic::ApplyCh(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ApplyChGen(feval, contrib);
  }

  void EdgeHierarchic::ApplyChTranspose(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ApplyChTransposeGen(feval, contrib);
  }

  void EdgeHierarchic::ApplyChTranspose(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ApplyChTransposeGen(feval, contrib);
  }

  void EdgeHierarchic::ComputeIntegralGradientRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeIntegralGradientGen(feval, contrib);
  }

  void EdgeHierarchic::ComputeIntegralGradientRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeIntegralGradientGen(feval, contrib);
  }

  
  //! displays information about object EdgeGauss
  ostream& operator<<(ostream& out, const EdgeHierarchic& e)
  {
    out<<static_cast<const ElementGeomReference<Dimension1>& >(e); return out;
    out<<"Stiffness Matrix "<<endl;
    out<<e.stiffness_matrix<<endl;
    return out;
  }

  
  /***********************************
   * ElementReference_Dim<Dimension1> *
   ***********************************/

  
  void ElementReference_Dim<Dimension1>::SolveMassMatrix(VectReal_wp& x) const
  {
    if (elt_geom.LumpedMassMatrix())
      {
	for (int i = 0; i < elt_geom.Weights.GetM(); i++)
	  x(i) /= elt_geom.Weights(i);

	return;
      }

    cout << "Not implemented" << endl;
    abort();
  }

  
  void ElementReference_Dim<Dimension1>::SolveMassMatrix(VectComplex_wp& x) const
  {
    if (elt_geom.LumpedMassMatrix())
      {
	for (int i = 0; i < elt_geom.Weights.GetM(); i++)
	  x(i) /= elt_geom.Weights(i);

	return;
      }

    cout << "Not implemented" << endl;
    abort();
  }


  void ElementReference_Dim<Dimension1>::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp&) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  void ElementReference_Dim<Dimension1>::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp&) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  void ElementReference_Dim<Dimension1>::MltMassMatrix(VectReal_wp& x) const
  {
    if (elt_geom.LumpedMassMatrix())
      {
	for (int i = 0; i < elt_geom.Weights.GetM(); i++)
	  x(i) *= elt_geom.Weights(i);

	return;
      }

    cout << "Not implemented" << endl;
    abort();
  }

  
  void ElementReference_Dim<Dimension1>::MltMassMatrix(VectComplex_wp& x) const
  {
    if (elt_geom.LumpedMassMatrix())
      {
	for (int i = 0; i < elt_geom.Weights.GetM(); i++)
	  x(i) *= elt_geom.Weights(i);
	
	return;
      }

    cout << "Not implemented" << endl;
    abort();
  }

  
  void ElementReference_Dim<Dimension1>
  ::ComputeMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
		      const VectReal_wp & coef) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  void ElementReference_Dim<Dimension1>
  ::IntegrateMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
			const VectReal_wp & coef) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  //! allocates a new projector for interpolation
  FiniteElementProjector* ElementReference_Dim<Dimension1>::GetNewNodalInterpolation() const
  {
    return elt_geom.GetNewNodalInterpolation();
  }
  
    
  //! computes local prolongation from a finite element to another one
  void ElementReference_Dim<Dimension1>
  ::ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
			     const ElementReference_Dim<Dimension1>& FaceCoarse,
			     const ElementReference_Dim<Dimension1>& FaceFine) const
  {
    // computing prolongation operation as a matrix
    LocalProlongation.Reallocate(FaceFine.GetNbDof(), FaceCoarse.GetNbDof());
    VectReal_wp phi;
    for (int j = 0; j < FaceFine.GetNbDof(); j++)
      {
	FaceCoarse.ComputeValuesPhiRef(FaceFine.PointsDof()(j), phi);
	for (int k = 0; k < FaceCoarse.GetNbDof(); k++)
	  LocalProlongation(j, k) = phi(k);
      }

    // we initialize proj as well
    proj.Init(FaceCoarse.GetGeometricElement(), FaceFine.PointsDof());
  }

}

#define MONTJOIE_FILE_EDGE_REFERENCE_CXX
#endif
