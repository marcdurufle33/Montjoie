#ifndef MONTJOIE_FILE_QUADRANGLE_GEOM_REFERENCE_CXX

namespace Montjoie
{
  
  /****************************
   * Initialization functions *
   ****************************/
  

  //! default constructor
  QuadrangleGeomReference::QuadrangleGeomReference() : ElementGeomReference<Dimension2>()
  {
    // four normales of the square
    this->normale.Reallocate(4);
    this->normale(0).Init(0,-1);
    this->normale(1).Init(1,0);
    this->normale(2).Init(0,1);
    this->normale(3).Init(-1,0); 
    
    this->hybrid_type_elt = 1;
    this->nb_vertices_elt = 4;
    this->nb_boundaries_elt = 4;
    this->nb_edges_elt = 4;
  }
  
    
  //! constructing finite element
  void QuadrangleGeomReference::
  ConstructFiniteElement(int rgeom)
  {
    this->order_geom = rgeom;
    
    ConstructNodalShapeFunctions(rgeom);
    
    ComputeCurvedTransformation();
    ComputeCoefficientTransformation();
  }
  

  size_t QuadrangleGeomReference::GetMemorySize() const
  {
    size_t taille = ElementGeomReference<Dimension2>::GetMemorySize();
    taille += SymmetricNodes.GetMemorySize();
    taille += coefFi.GetMemorySize();
    taille += Seldon::GetMemorySize(coefFi_curve);
    taille += Seldon::GetMemorySize(numberFi_curve);
    return taille;
  }

  
  //! changing nodal points of the quadrangle
  void QuadrangleGeomReference::SetNodalPoints(int r, const Vector<R2>& pts)
  {
    this->order_geom = r;
    this->points_nodal1d.Reallocate(r+1);
    this->points_nodal_nd = pts;
    // extracting 1-d nodal points
    this->points_nodal1d.Fill(0);
    int nb = 0;
    for (int i = 0; i < pts.GetM(); i++)
      if (abs(pts(i)(1)) < 1e4*epsilon_machine)
	this->points_nodal1d(nb++) = pts(i)(0);

    Sort(this->points_nodal1d.GetM(), this->points_nodal1d);

    this->NumNodes2D.Reallocate(r+1, r+1);
    this->CoordinateNodes.Reallocate((r+1)*(r+1), 2);
    // constructing NumNodes2D and CoordinateNodes
    for (int i = 0; i < pts.GetM(); i++)
      {
	int i1 = -1, i2 = -1;
	for (int j = 0; j <= this->order_geom; j++)
	  if (abs(this->points_nodal1d(j) - pts(i)(0)) < 1e4*epsilon_machine)
	    i1 = j;
	
	for (int j = 0; j <= this->order_geom; j++)
	  if (abs(this->points_nodal1d(j) - pts(i)(1)) < 1e4*epsilon_machine)
	    i2 = j;
	
	if ((i1 == -1)||(i2 == -1))
	  {
	    cout << "This method only accepts tensorized points "<<endl;
	    abort();
	  }
	
	this->NumNodes2D(i1, i2) = i;
	this->CoordinateNodes(i, 0) = i1;
	this->CoordinateNodes(i, 1) = i2;
      }
    
    // deducing EdgesNode
    this->EdgesNode.Reallocate(r+1, 4);
    for (int i = 0; i <= this->order_geom; i++)
      {
	this->EdgesNode(i, 0) = this->NumNodes2D(i, 0);
	this->EdgesNode(i, 1) = this->NumNodes2D(this->order_geom, i);
	this->EdgesNode(i, 2) = this->NumNodes2D(this->order_geom-i,this->order_geom);
	this->EdgesNode(i, 3) = this->NumNodes2D(0, this->order_geom-i);    
      }

    this->SymmetricNodes.Reallocate((r+1)*(r+1));
    for (int i = 0; i <= this->order_geom; i++)
      for (int j = 0; j <= this->order_geom; j++)
	SymmetricNodes(this->NumNodes2D(i,j)) = this->NumNodes2D(j,i);
    
    this->quadrature_equal_nodal = false;
    this->dof_equal_nodal = false;

    this->edge_geom.ConstructFiniteElement(r);
    lob_geom.AffectPoints(this->points_nodal1d);
    
    ComputeCurvedTransformation();
    ComputeCoefficientTransformation();
  }
  
  
  //! construction of nodal points on the square
  /*!
    \param[in] order_ order of approximation
   */  
  void QuadrangleGeomReference::ConstructNodalShapeFunctions(int r)
  {
    this->edge_geom.ConstructFiniteElement(r);
    lob_geom.ConstructQuadrature(r, lob_geom.QUADRATURE_LOBATTO);
    lob_geom.ComputeGradPhi(1e3*epsilon_machine);
    
    this->points_nodal1d = lob_geom.Points();
    
    MeshNumbering<Dimension2>::
      ConstructQuadrilateralNumbering(r, this->NumNodes2D, this->CoordinateNodes);
    
    this->EdgesNode.Reallocate(r+1, 4);
    for (int i = 0; i <= r; i++)
      {
	this->EdgesNode(i, 0) = this->NumNodes2D(i, 0);
	this->EdgesNode(i, 1) = this->NumNodes2D(r, i);
	this->EdgesNode(i, 2) = this->NumNodes2D(r-i, r);
	this->EdgesNode(i, 3) = this->NumNodes2D(0, r-i);    
      }
    
    this->points_nodal_nd.Reallocate((r+1)*(r+1));
    this->weights_nodal_nd.Reallocate((r+1)*(r+1));
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  this->points_nodal_nd(this->NumNodes2D(i,j))
	    = R2(lob_geom.Points(i), lob_geom.Points(j));

	  this->weights_nodal_nd(this->NumNodes2D(i,j))
	    = lob_geom.Weights(i) * lob_geom.Weights(j);
	}
    
    SymmetricNodes.Reallocate(this->points_nodal_nd.GetM());
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	SymmetricNodes(this->NumNodes2D(i,j)) = this->NumNodes2D(j,i);
  }
  
  
  //! computation of curved transformation (Gordon-Hall)  
  void QuadrangleGeomReference::ComputeCurvedTransformation()
  {    
    numberFi_curve.Reallocate(this->points_nodal_nd.GetM());
    coefFi_curve.Reallocate(this->points_nodal_nd.GetM());
    
    // use of Gordon-Hall transformation to get nodal points inside the quadrilateron
    int r = this->order_geom; Real_wp x, y;
    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)      
	{
	  int node = this->NumNodes2D(i, j);
	  x = this->points_nodal_nd(node)(0);
	  y = this->points_nodal_nd(node)(1);
	  
	  // Fi = (1-y) e0(x) + y e2(1-x) + x e1(y) + (1-x) e3(y)
	  //      - (1-x)(1-y) S0 - x(1-y) S1 - xy S2 - (1-x)y S3
	  numberFi_curve(node).Reallocate(8);
	  coefFi_curve(node).Reallocate(8);
	  numberFi_curve(node)(0) = this->NumNodes2D(i, 0);
	  coefFi_curve(node)(0) = 1.0-y;
	  
	  numberFi_curve(node)(1) = this->NumNodes2D(i, r);
	  coefFi_curve(node)(1) = y;

	  numberFi_curve(node)(2) = this->NumNodes2D(r, j);
	  coefFi_curve(node)(2) = x;

	  numberFi_curve(node)(3) = this->NumNodes2D(0, j);
	  coefFi_curve(node)(3) = 1.0-x;

	  numberFi_curve(node)(4) = this->NumNodes2D(0, 0);
	  coefFi_curve(node)(4) = -(1.0-y)*(1.0-x);
	  
	  numberFi_curve(node)(5) = this->NumNodes2D(r, 0);
	  coefFi_curve(node)(5) = -x*(1.0-y);

	  numberFi_curve(node)(6) = this->NumNodes2D(r, r);
	  coefFi_curve(node)(6) = -x*y;

	  numberFi_curve(node)(7) = this->NumNodes2D(0, r);
	  coefFi_curve(node)(7) = -(1.0-x)*y;
	}
    
  }
  
  
  //! computation of coefficients in order to obtain faster evaluations of Fi and DFi  
  void QuadrangleGeomReference::ComputeCoefficientTransformation()
  {
    // for linear case
    // F_i = (1-x)*(1-y)*A0 + x*(1-y)*A1 + x*y*A2 + (1-x)*y*A3
    int nb_points_all = this->points_nodal_nd.GetM() 
      + this->points_dof2d.GetM() + this->points2d.GetM();
    
    coefFi.Reallocate(nb_points_all, 4);
    Real_wp x, y;
    
    VectR2 all_points = this->points_nodal_nd;
    all_points.PushBack(this->points2d); all_points.PushBack(this->points_dof2d);
    for (int i = 0; i < nb_points_all; i++)
      {
	x = all_points(i)(0);
	y = all_points(i)(1);
	
	coefFi(i,0) = (1-x)*(1-y);
	coefFi(i,1) = x*(1-y);
	coefFi(i,2) = x*y;
	coefFi(i,3) = (1-x)*y;
      }
    
    numberFi_curve.Resize(all_points.GetM());
    coefFi_curve.Resize(all_points.GetM());
    Real_wp value; VectReal_wp phi;
    for (int i = this->points_nodal_nd.GetM(); i < nb_points_all; i++)
      {
	// first test to count number of non-zeros
	int nnz = 0;
	ComputeValuesPhiNodalRef(all_points(i), phi);
	for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
	  {
	    value = phi(j);
	    if (abs(value) > 1e3*epsilon_machine)
	      nnz++;
	  }
	
	coefFi_curve(i).Reallocate(nnz);
	numberFi_curve(i).Reallocate(nnz);
	nnz = 0;
	for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
	  {
	    value = phi(j);
	    if (abs(value) > 1e3*epsilon_machine)
	      {
		coefFi_curve(i)(nnz) = value;
		numberFi_curve(i)(nnz) = j;
		nnz++;
	      }
	  }
      }
  }


  //! returns the minimal length of the element
  Real_wp QuadrangleGeomReference::GetMinimalSize(const VectR2& s) const
  {
    Real_wp h = Distance(s(0), s(1));
    h = min(h, Distance(s(1), s(2)));
    h = min(h, Distance(s(2), s(3)));
    h = min(h, Distance(s(0), s(3)));
    return h;
  }
  

  //! returns true if the point is outside the unit square
  /*!
    \param[in] Xn local coordinates of the point
    \param[in] epsilon threshold 
    \return true if the point is outside
   */
  bool QuadrangleGeomReference
  ::OutsideReferenceElement(const VectR2& s, const R2& Xn, const Real_wp& epsilon) const
  {
    Real_wp Rmax = max(abs(s(0)(0)), abs(s(0)(1)));
    Rmax = max(Rmax, abs(s(1)(0))); Rmax = max(Rmax, abs(s(1)(1)));
    Rmax = max(Rmax, abs(s(2)(0))); Rmax = max(Rmax, abs(s(2)(1)));
    Rmax = max(Rmax, abs(s(3)(0))); Rmax = max(Rmax, abs(s(3)(1)));

    Real_wp dx = s(0).Distance(s(1)), dy = s(0).Distance(s(3));
    if ((Rmax == Real_wp(0)) || (dx == Real_wp(0)) || (dy == Real_wp(0)))
      return true;
    
    if (dx*Xn(0) < -epsilon*Rmax)
      return true;

    if (dy*Xn(1) < -epsilon*Rmax)
      return true;

    if (dx*(Xn(0) - Real_wp(1)) > epsilon*Rmax)
      return true;

    if (dy*(Xn(1) - Real_wp(1)) > epsilon*Rmax)
      return true;

    return false;
  }
  
  
  //! returns distance of the local point to the boundary of the element
  /*!
    \param[in] pointloc local coordinates of the point
    \return distance to the boundary (negative if the point is outside)
   */
  Real_wp QuadrangleGeomReference::GetDistanceToBoundary(const R2& pointloc) const
  {
    Real_wp distance = pointloc(0);
    distance = min(distance,Real_wp(Real_wp(1) - pointloc(0)));
    distance = min(distance,pointloc(1));
    distance = min(distance,Real_wp(Real_wp(1) - pointloc(1)));
    
    return distance;
  }
  
  
  //! returns distance of the local point to a given boundary of the element
  Real_wp QuadrangleGeomReference::GetDistanceToBoundary(const R2& pointloc, int num_loc) const
  {
    Real_wp distance(0);
    switch (num_loc)
      {
      case 0:
	{
	  distance = abs(pointloc(1));
	  distance = max(distance, -pointloc(0));
	  distance = max(distance, pointloc(0)-Real_wp(1));
	}
	break;
      case 1:
	{
	  distance = abs(pointloc(0)-Real_wp(1));
	  distance = max(distance, -pointloc(1));
	  distance = max(distance, pointloc(1)-Real_wp(1));
	}
	break;
      case 2:
	{
	  distance = abs(pointloc(1)-Real_wp(1));
	  distance = max(distance, -pointloc(0));
	  distance = max(distance, pointloc(0)-Real_wp(1));
	}
	break;
      case 3:
	{
	  distance = abs(pointloc(0));
	  distance = max(distance, -pointloc(1));
	  distance = max(distance, pointloc(1)-Real_wp(1));
	}
	break;	
      }
    
    return distance;
  }


  //! projects the local point on the boundary, if outside the unit square
  int QuadrangleGeomReference::ProjectPointOnBoundary(R2& pointloc) const
  {
    pointloc(0) = std::max(Real_wp(0), pointloc(0));
    pointloc(0) = std::min(Real_wp(1), pointloc(0)); 
    pointloc(1) = std::max(Real_wp(0), pointloc(1));
    pointloc(1) = std::min(Real_wp(1), pointloc(1));
    
    if (abs(pointloc(1)) <= R2::threshold)
      return 0;

    if (abs(pointloc(0)-1.0) <= R2::threshold)
      return 1;

    if (abs(pointloc(1)-1.0) <= R2::threshold)
      return 2;

    if (abs(pointloc(0)) <= R2::threshold)
      return 3;
    
    return -1;
  }
  
    
  /****************
   * Fj transform *
   ****************/
  
  
  //! transformation Fi in the case of curved quadrilateral
  void QuadrangleGeomReference::FjCurve(const SetPoints<Dimension2>& PTReel,
                                        const R2& pointloc, R2& res) const
  {
    res.Fill(0); VectReal_wp phi;
    ComputeValuesPhiNodalRef(pointloc, phi);
    for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
      {
	Real_wp coef = phi(node);
	Add(coef, PTReel.GetPointNodal(node), res);
      }
  }
  
  
  //! transformation DFi in the case of curved quadrilateral
  void QuadrangleGeomReference::DFjCurve(const SetPoints<Dimension2>& PTReel,
                                         const R2& pointloc, Matrix2_2& res) const
  {
    int r = this->order_geom;
    VectReal_wp lambda_x(r+1),lambda_y(r+1), dlambda_x(r+1), dlambda_y(r+1);
    
    for (int i = 0; i <= r; i++)
      {
	lambda_x(i) = lob_geom.EvaluatePhi(i, pointloc(0));
	lambda_y(i) = lob_geom.EvaluatePhi(i, pointloc(1));
	dlambda_x(i) = lob_geom.EvaluatePhiGrad(i, pointloc(0));
	dlambda_y(i) = lob_geom.EvaluatePhiGrad(i, pointloc(1));
      }
    
    res.Fill(0);
    R2 temp;
    int node = 0;
    for (int i = 0; i <= r; i++)
      {
	temp.Fill(0);
	for (int j = 0; j <= r; j++)
	  {
	    node = this->NumNodes2D(i,j);
	    Add(lambda_y(j), PTReel.GetPointNodal(node), temp);
	  }
	res(0,0) += temp(0)*dlambda_x(i);
	res(1,0) += temp(1)*dlambda_x(i);
      }
    
    for (int i = 0; i <= r; i++)
      {
	temp.Fill(0);
	for (int j = 0; j <= r; j++)
	  {
	    node = this->NumNodes2D(i,j);
	    Add(dlambda_y(j), PTReel.GetPointNodal(node), temp);
	  }
	res(0,1) += temp(0)*lambda_x(i);
	res(1,1) += temp(1)*lambda_x(i);
      }
    
  }


  /**********************
   * FjElem and DFjElem *
   **********************/
  
  
  //! computes res = Fi(point) for all reference points
  /*!
    \param[in] s list of vertices of the quad
    \param[out] res "reference points" after transformation Fi
    \param[in] mesh quadrilateral mesh considered
    \param[in] nquad quadrilateral element number
   */
  void QuadrangleGeomReference
  ::FjElem(const VectR2& s, SetPoints<Dimension2>& res,
           const Mesh<Dimension2>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
	if (!this->quadrature_equal_nodal)
	  FjElemQuadratureCurve(s, res, mesh, nquad);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  FjElemDofCurve(s, res, mesh, nquad);
      }
    else
      {
	FjElemNodalLinear(s, res);
	if (!this->quadrature_equal_nodal)
	  FjElemQuadratureLinear(s, res);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  FjElemDofLinear(s, res);
      }
  }
  
  
  //! computes res = Fi(point) for all nodal points
  void QuadrangleGeomReference
  ::FjElemNodal(const VectR2& s, SetPoints<Dimension2>& res,
                const Mesh<Dimension2>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
    else
      FjElemNodalLinear(s, res);
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! computes res = Fi(point) for all nodal points
  void QuadrangleGeomReference::
  FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
	      const Mesh<Dimension3>& mesh, int i, const Face<Dimension3>& f) const
  {
    if (mesh.IsBoundaryCurved(i))
      FjElemNodalCurve(s, res, mesh, i, f);
    else
      FjElemNodalLinear(s, res);
  }

  void QuadrangleGeomReference
  ::DFjElemNodal(const VectR3& pts_nodal, bool is_curved,
		 VectR3& dF_dx, VectR3& dF_dy) const
  {
    dF_dx.Reallocate(this->points_nodal_nd.GetM());
    dF_dy.Reallocate(this->points_nodal_nd.GetM());

    if (is_curved)
      {
	for (int k = 0; k <= this->order_geom; k++)
	  for (int l = 0; l <= this->order_geom; l++)
	    {
	      int node = this->NumNodes2D(k,l);
	      dF_dx(node).Fill(0);
	      dF_dy(node).Fill(0);
	      for (int i = 0; i <= this->order_geom; i++)
		{
		  int node2 = this->NumNodes2D(i,l);
		  Add(lob_geom.GradPhi(i, k), pts_nodal(node2), dF_dx(node)); 
		}
	      
	      for (int j = 0; j <= this->order_geom; j++)
		{
		  int node2 = this->NumNodes2D(k,j);
		  Add(lob_geom.GradPhi(j, l), pts_nodal(node2), dF_dy(node));
		}
	    }
      }
    else
      {
	Real_wp x, y, one(1);
	for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
	  {
	    x = this->points_nodal_nd(i)(0);
	    y = this->points_nodal_nd(i)(1);
	    dF_dx(i).Zero();
	    Add(one-y, pts_nodal(1), dF_dx(i));
	    Add(-(one-y), pts_nodal(0), dF_dx(i));
	    Add(y, pts_nodal(2), dF_dx(i));
	    Add(-y, pts_nodal(3), dF_dx(i));

	    dF_dy(i).Zero();
	    Add(one-x, pts_nodal(3), dF_dy(i));
	    Add(-(one-x), pts_nodal(0), dF_dy(i));
	    Add(x, pts_nodal(2), dF_dy(i));
	    Add(-x, pts_nodal(1), dF_dy(i));	    
	  }
      }
  }

  void QuadrangleGeomReference
  ::DFjElemQuadrature(const VectR3& pts_nodal, bool is_curved,
		      const VectR3& dF_dx_node, const VectR3& dF_dy_node,
		      VectR3& dF_dx, VectR3& dF_dy) const
  {
    dF_dx.Reallocate(this->points2d.GetM());
    dF_dy.Reallocate(this->points2d.GetM());

    if (is_curved)
      {
	for (int i = 0; i < this->points2d.GetM(); i++)
	  {
	    dF_dx(i).Fill(0);
	    dF_dy(i).Fill(0);
	    int N = this->points_nodal_nd.GetM();
	    for (int node = 0; node < coefFi_curve(N+i).GetM(); node++)
	      {
		Real_wp coef = coefFi_curve(N+i)(node);
		Add(coef, dF_dx_node(numberFi_curve(N+i)(node)), dF_dx(i));
		Add(coef, dF_dy_node(numberFi_curve(N+i)(node)), dF_dy(i));
	      }
	  }
      }
    else
      {
	Real_wp x, y, one(1);
	for (int i = 0; i < this->points2d.GetM(); i++)
	  {
	    x = this->points2d(i)(0);
	    y = this->points2d(i)(1);
	    dF_dx(i).Zero();
	    Add(one-y, pts_nodal(1), dF_dx(i));
	    Add(-(one-y), pts_nodal(0), dF_dx(i));
	    Add(y, pts_nodal(2), dF_dx(i));
	    Add(-y, pts_nodal(3), dF_dx(i));

	    dF_dy(i).Zero();
	    Add(one-x, pts_nodal(3), dF_dy(i));
	    Add(-(one-x), pts_nodal(0), dF_dy(i));
	    Add(x, pts_nodal(2), dF_dy(i));
	    Add(-x, pts_nodal(1), dF_dy(i));	    
	  }
      }
  }

  //! computes res = Fi(point) for all quadrature points
  void QuadrangleGeomReference::
  FjElemQuadrature(const VectR3& s, SetPoints<Dimension3>& res,
		   const Mesh<Dimension3>& mesh, int nquad, const Face<Dimension3>& f) const
  {
    if (mesh.IsBoundaryCurved(nquad))
      {
	FjElemNodalCurve(s, res, mesh, nquad, f);
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  FjElemQuadratureCurve(s, res, mesh, nquad);
      }
    else
      {
	FjElemNodalLinear(s, res);
	FjElemQuadratureLinear(s, res);
      }
  }
#endif
  
  
  //! computes res = Fi(point) for all quadrature points
  void QuadrangleGeomReference::
  FjElemQuadrature(const VectR2& s, SetPoints<Dimension2>& res,
		   const Mesh<Dimension2>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  FjElemQuadratureCurve(s,res,mesh,nquad);
      }
    else
      FjElemQuadratureLinear(s,res);
  }
  
  
  //! computes res = Fi(point) for all dof points
  void QuadrangleGeomReference::
  FjElemDof(const VectR2& s, SetPoints<Dimension2>& res,
	    const Mesh<Dimension2>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else
	  FjElemDofCurve(s, res, mesh, nquad);
      }
    else
      FjElemDofLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all reference points
  /*!
    \param[in] s list of vertices of the quad
    \param[in] PTReel "reference points" after transformation Fi
    \param[out] res jacobian matrices DFi
    \param[in] mesh quadrilateral mesh considered
    \param[in] nquad quadrilateral element number
   */
  void QuadrangleGeomReference::
  DFjElem(const VectR2& s, const SetPoints<Dimension2>& PTReel,
	  SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (!this->quadrature_equal_nodal)
	  DFjElemQuadratureCurve(s, PTReel, res, mesh, nquad);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  DFjElemDofCurve(s, PTReel, res, mesh, nquad);
      }
    else
      {
	DFjElemNodalLinear(s, res);
	if (!this->quadrature_equal_nodal)
	  DFjElemQuadratureLinear(s, res);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  DFjElemDofLinear(s, res);
      }
  }
  
  
  //! computes res = DFi(point) for all nodal points
  void QuadrangleGeomReference::
  DFjElemNodal(const VectR2& s, const SetPoints<Dimension2>& PTReel,
	       SetMatrices<Dimension2>& res,
	       const Mesh<Dimension2>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      DFjElemNodalCurve(s,PTReel,res,mesh,nquad);
    else
      DFjElemNodalLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all quadrature points
  void QuadrangleGeomReference::
  DFjElemQuadrature(const VectR2& s, const SetPoints<Dimension2>& PTReel, 
		    SetMatrices<Dimension2>& res,
		    const Mesh<Dimension2>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s,PTReel,res,mesh,nquad);
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  DFjElemQuadratureCurve(s, PTReel, res, mesh, nquad);
      }
    else
      DFjElemQuadratureLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all dof points
  void QuadrangleGeomReference::
  DFjElemDof(const VectR2& s, const SetPoints<Dimension2>& PTReel,
	     SetMatrices<Dimension2>& res,
	     const Mesh<Dimension2>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else
	  DFjElemDofCurve(s, PTReel, res, mesh, nquad);
      }
    else
      DFjElemDofLinear(s,res);
  }
  
  
  //! nodal points in the case of straight quadrilateral 
  template<class DimensionB>
  void QuadrangleGeomReference
  ::FjElemNodalLinear(const typename DimensionB::VectR_N& s,
                      SetPoints<DimensionB>& res) const
  {
    typename DimensionB::R_N res_n;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
      {
	res_n.Fill(0);
	Add(coefFi(i,0), s(0), res_n);
	Add(coefFi(i,1), s(1), res_n);
	Add(coefFi(i,2), s(2), res_n);
	Add(coefFi(i,3), s(3), res_n);
	
	res.SetPointNodal(i,res_n);
      }
  }
  
  
  //! quadrature points in the case of straight quadrilateral 
  void QuadrangleGeomReference
  ::FjElemQuadratureLinear(const VectR2& s,
                           SetPoints<Dimension2>& res) const
  {
    R2 res_n; int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	res_n.Fill(0);
	Add(coefFi(nb_points_nodal_elt+i, 0), s(0), res_n);
	Add(coefFi(nb_points_nodal_elt+i, 1), s(1), res_n);
	Add(coefFi(nb_points_nodal_elt+i, 2), s(2), res_n);
	Add(coefFi(nb_points_nodal_elt+i, 3), s(3), res_n);
	
	res.SetPointQuadrature(i,res_n);
      }
  }


  //! quadrature points in the case of straight quadrilateral 
  void QuadrangleGeomReference
  ::FjElemQuadratureLinear(const VectR3& s,
                           SetPoints<Dimension3>& res) const
  {
    R3 res_n; int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	res_n.Fill(0);
	Add(coefFi(nb_points_nodal_elt+i, 0), s(0), res_n);
	Add(coefFi(nb_points_nodal_elt+i, 1), s(1), res_n);
	Add(coefFi(nb_points_nodal_elt+i, 2), s(2), res_n);
	Add(coefFi(nb_points_nodal_elt+i, 3), s(3), res_n);
	
	res.SetPointQuadrature(i,res_n);
      }
  }

  
  //! dof points in the case of straight quadrilateral 
  void QuadrangleGeomReference
  ::FjElemDofLinear(const VectR2& s, SetPoints<Dimension2>& res) const
  {
    R2 res_n;
    int offset = this->points_nodal_nd.GetM() + this->points2d.GetM();
    res.ReallocatePointsDof(this->points_dof2d.GetM());
    for (int i = 0; i < this->points_dof2d.GetM(); i++)
      {
	res_n.Fill(0);
	Add(coefFi(offset+i,0), s(0), res_n);
	Add(coefFi(offset+i,1), s(1), res_n);
	Add(coefFi(offset+i,2), s(2), res_n);
	Add(coefFi(offset+i,3), s(3), res_n);
	
	res.SetPointDof(i,res_n);
      }
  }
  
  
  //! transformation DFi for nodal points in the case of straight quadrilateral
  void QuadrangleGeomReference
  ::DFjElemNodalLinear(const VectR2& s, SetMatrices<Dimension2>& res) const
  {
    Matrix2_2 mat; Real_wp x,y,one(1);
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
      {
	x = this->points_nodal_nd(i)(0); y = this->points_nodal_nd(i)(1);
	mat(0,0) = (one-y)*(s(1)(0)-s(0)(0))+y*(s(2)(0)-s(3)(0));
	mat(1,0) = (one-y)*(s(1)(1)-s(0)(1))+y*(s(2)(1)-s(3)(1));
	mat(0,1) = (one-x)*(s(3)(0)-s(0)(0))+x*(s(2)(0)-s(1)(0));
	mat(1,1) = (one-x)*(s(3)(1)-s(0)(1))+x*(s(2)(1)-s(1)(1));
	res.SetPointNodal(i,mat);
      }
  }
  
  
  //! transformation DFi for quadrature points in the case of straight quadrilateral  
  void QuadrangleGeomReference
  ::DFjElemQuadratureLinear(const VectR2& s,
                            SetMatrices<Dimension2>& res) const
  {
    Matrix2_2 mat; Real_wp x,y,one(1);
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	x = this->points2d(i)(0); y = this->points2d(i)(1);
	mat(0,0)=(one-y)*(s(1)(0)-s(0)(0))+y*(s(2)(0)-s(3)(0));
	mat(1,0)=(one-y)*(s(1)(1)-s(0)(1))+y*(s(2)(1)-s(3)(1));
	mat(0,1)=(one-x)*(s(3)(0)-s(0)(0))+x*(s(2)(0)-s(1)(0));
	mat(1,1)=(one-x)*(s(3)(1)-s(0)(1))+x*(s(2)(1)-s(1)(1));
        res.SetPointQuadrature(i,mat);
      }
  }

  
  //! transformation DFi for dof points in the case of straight quadrilateral  
  void QuadrangleGeomReference
  ::DFjElemDofLinear(const VectR2& s, SetMatrices<Dimension2>& res) const
  {
    Matrix2_2 mat; Real_wp x,y,one(1);
    res.ReallocatePointsDof(this->points_dof2d.GetM());
    for (int i = 0; i < this->points_dof2d.GetM(); i++)
      {
	x = this->points_dof2d(i)(0); y = this->points_dof2d(i)(1);
	mat(0,0)=(one-y)*(s(1)(0)-s(0)(0))+y*(s(2)(0)-s(3)(0));
	mat(1,0)=(one-y)*(s(1)(1)-s(0)(1))+y*(s(2)(1)-s(3)(1));
	mat(0,1)=(one-x)*(s(3)(0)-s(0)(0))+x*(s(2)(0)-s(1)(0));
	mat(1,1)=(one-x)*(s(3)(1)-s(0)(1))+x*(s(2)(1)-s(1)(1));
	res.SetPointDof(i,mat);
      }
  }
  
  
  //! transformation DFi for boundary points
  void QuadrangleGeomReference::
  FjSurfaceElemLinear(const VectR2& s, SetPoints<Dimension2>& res, int ne) const
  {
    // we get the two extremities of the edge
    R2 ptA = s(ne);  R2 ptB = s((ne+1)%4), ptC;
    // we project quadrature points
    res.ReallocatePointsQuadratureBoundary(this->points1d.GetM());
    for (int k = 0; k < this->points1d.GetM(); k++)
      {
	ptC.Fill(0); Add(Real_wp(1)-this->points1d(k), ptA, ptC);
	Add(this->points1d(k), ptB, ptC);
	res.SetPointQuadratureBoundary(k, ptC);
      }
  }
  
  
  //! transformation Fi for nodal points in the case of curved quadrilateral
  /*!
    \param[in] s four vertices of the quad
    \param[out] res references points after the transformation Fi
    \param[in] mesh considered mesh
    \param[in] nquad quadrilateral number in the mesh
    This method uses Gordon-Hall transformation on quads
   */
  template<class DimensionB, class TypeEntity>
  void QuadrangleGeomReference::
  FjElemNodalCurve(const typename DimensionB::VectR_N& s, SetPoints<DimensionB>& res,
		   const Mesh<DimensionB>& mesh, int nquad, const TypeEntity& elt) const
  {
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    
    // we consider that the order of approximation of curves in the mesh
    // can be different from order of finite element
    int r = mesh.GetGeometryOrder();
          
    // vertices
    res.SetPointNodal(0, s(0));
    res.SetPointNodal(1, s(1));
    res.SetPointNodal(2, s(2));
    res.SetPointNodal(3, s(3));
    typename DimensionB::R_N res_n;
    
    if ((this->CoefCurve1D.GetM() != this->order_geom-1)||(this->CoefCurve1D.GetN() != r+1))
      {
	cout<<"Did you call SetMesh after ConstructFiniteElement ?"<<endl;
	abort();
      }
    
    // edges
    for (int j = 0; j < 4; j++)
      {
	// we get the edge number and the reference
	int num_edge = elt.numEdge(j);
        int num_edge_ref = -1;
	if ((num_edge >= 0) && (num_edge < mesh.GetNbEdgesRef()))
          num_edge_ref = num_edge;        
	
	if (num_edge_ref == -1)
	  {
	    // case of a straight edge
	    // two extremities of the edge
	    typename DimensionB::R_N ptA = s(j);
	    typename DimensionB::R_N ptB = s((j+1)%4);	    
	    for (int k = 1; k < this->order_geom; k++)
	      {
		res_n.Zero();
		// linear interpolation
		Add(1.0-this->points_nodal1d(k), ptA, res_n); 
		Add(this->points_nodal1d(k), ptB, res_n);
		res.SetPointNodal(this->EdgesNode(k, j), res_n); 
	      }
	  }
	else
	  {
	    // first step : we get the points on the curved edge
	    Vector<typename DimensionB::R_N> PtsEdge(r+1);
	    PtsEdge(0) = s(j);
	    PtsEdge(r) = s((j+1)%4);
	    if (elt.GetOrientationEdge(j))
	      for (int k = 1; k < r; k++)
		PtsEdge(k) = mesh.GetPointInsideEdge(num_edge_ref, k-1);
	    else
	      for (int k = 1; k < r; k++)
		PtsEdge(k) = mesh.GetPointInsideEdge(num_edge_ref, r-k-1);
	    
	    // second step : interpolation to nodal points of the finite element
	    for (int k = 1; k < this->order_geom; k++)
	      {
		res_n.Zero();
		for (int p = 0; p <= r; p++)
		  Add(this->CoefCurve1D(k-1, p), PtsEdge(p), res_n);
		
		res.SetPointNodal(this->EdgesNode(k, j), res_n);
	      }
	  }
      }
    
    // for other nodal points, we use coefFi_curve
    for (int node = 4*this->order_geom; node < this->points_nodal_nd.GetM(); node++)
      {
	res_n.Zero();
	for (int k = 0; k < numberFi_curve(node).GetM(); k++)
	  Add(coefFi_curve(node)(k), res.GetPointNodal(numberFi_curve(node)(k)), res_n);
	
	res.SetPointNodal(node, res_n);
      }

  }

  
  //! transformation Fi for quadrature points in the case of curved quadrilateral
  void QuadrangleGeomReference::
  FjElemQuadratureCurve(const VectR2& s, SetPoints<Dimension2>& res,
			const Mesh<Dimension2>& mesh,int nquad) const
  {
    R2 res_n; int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < coefFi_curve(nb_points_nodal_elt+i).GetM(); node++)
	  {
	    Real_wp coef = coefFi_curve(nb_points_nodal_elt+i)(node);
	    Add(coef, res.GetPointNodal(numberFi_curve(nb_points_nodal_elt+i)(node)),
		res_n);
	  }
	res.SetPointQuadrature(i, res_n);
      }
  }


  //! transformation Fi for quadrature points in the case of curved quadrilateral
  void QuadrangleGeomReference::
  FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
			const Mesh<Dimension3>& mesh,int nquad) const
  {
    R3 res_n; int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	res_n.Zero();
	for (int node = 0; node < coefFi_curve(nb_points_nodal_elt+i).GetM(); node++)
	  {
	    Real_wp coef = coefFi_curve(nb_points_nodal_elt+i)(node);
	    Add(coef, res.GetPointNodal(numberFi_curve(nb_points_nodal_elt+i)(node)),
		res_n);
	  }
	res.SetPointQuadrature(i, res_n);
      }
  }

  
  //! transformation Fi for dof points in the case of curved quadrilateral  
  void QuadrangleGeomReference
  ::FjElemDofCurve(const VectR2& s, SetPoints<Dimension2>& res,
                   const Mesh<Dimension2>& mesh, int nquad) const
  {
    int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    int nb_points_quadrature = this->points2d.GetM();
    R2 res_n;  int offset = nb_points_nodal_elt+nb_points_quadrature;
    res.ReallocatePointsDof(this->points_dof2d.GetM());
    for (int i = 0; i < this->points_dof2d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < numberFi_curve(offset+i).GetM(); node++)
	  {
	    Real_wp coef = coefFi_curve(offset+i)(node);
	    Add(coef,res.GetPointNodal(numberFi_curve(offset+i)(node)),res_n);
	  }
	res.SetPointDof(i, res_n);
      }
  }
  
  
  //! transformation DFi for nodal points in the case of curved quadrilateral
  void QuadrangleGeomReference::
  DFjElemNodalCurve(const VectR2& s, const SetPoints<Dimension2>& PTReel, 
		    SetMatrices<Dimension2>& res,
		    const Mesh<Dimension2>& mesh, int nquad) const
  {
    Matrix2_2 res_n;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    // derivation of F_i = \sum A_i \varphi_i
    // where A_i are the nodal points already computed
    // and \varphi_i the basis functions
    for (int k = 0; k <= this->order_geom; k++)
      {
	for (int l = 0; l <= this->order_geom; l++)
	  {
	    int node = this->NumNodes2D(k,l);
	    res_n.Fill(0);
	    for (int i = 0; i <= this->order_geom; i++)
	      {
		int node2 = this->NumNodes2D(i,l);
		res_n(0,0) += lob_geom.GradPhi(i,k)*PTReel.GetPointNodal(node2)(0);
		res_n(1,0) += lob_geom.GradPhi(i,k)*PTReel.GetPointNodal(node2)(1);
	      }
            
	    for (int j = 0; j <= this->order_geom; j++)
	      {
		int node2 = this->NumNodes2D(k,j);
		res_n(0,1) += lob_geom.GradPhi(j,l)*PTReel.GetPointNodal(node2)(0);
		res_n(1,1) += lob_geom.GradPhi(j,l)*PTReel.GetPointNodal(node2)(1);
	      }
            
	    res.SetPointNodal(node,res_n);
	  }
      }
  }
  
  
  //! transformation DFi for quadrature points in the case of curved quadrilateral
  void QuadrangleGeomReference::
  DFjElemQuadratureCurve(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			 SetMatrices<Dimension2>& res,
			 const Mesh<Dimension2>& mesh, int nquad) const
  {
    Matrix2_2 res_n; int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < coefFi_curve(nb_points_nodal_elt+i).GetM(); node++)
	  {
	    Real_wp coef = coefFi_curve(nb_points_nodal_elt+i)(node);
	    Add(coef, res.GetPointNodal(numberFi_curve(nb_points_nodal_elt+i)(node)),
		res_n);
	  }
	res.SetPointQuadrature(i,res_n);
      }
  }

  
  //! transformation DFi for dof points in the case of curved quadrilateral
  void QuadrangleGeomReference::
  DFjElemDofCurve(const VectR2& s, const SetPoints<Dimension2>& PTReel,
		  SetMatrices<Dimension2>& res,
		  const Mesh<Dimension2>& mesh,int nquad) const
  {
    int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    int nb_points_quadrature = this->points2d.GetM();
    Matrix2_2 res_n; int offset = nb_points_nodal_elt + nb_points_quadrature;
    res.ReallocatePointsDof(this->points_dof2d.GetM());
    for (int i = 0; i < this->points_dof2d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < numberFi_curve(offset+i).GetM(); node++)
	  {
	    Real_wp coef = coefFi_curve(offset+i)(node);
	    Add(coef, res.GetPointNodal(numberFi_curve(offset+i)(node)), res_n);
	  }
	
	res.SetPointDof(i,res_n);
      }
  }
    

  //! transformation Fi in the case of straight quadrilateral
  void QuadrangleGeomReference
  ::FjLinear(const VectR2& s, const R2& point, R2& res) const
  {
    Real_wp one(1);
    res.Init(s(0)(0)*(one-point(0))*(one-point(1))+s(1)(0)*point(0)*(one-point(1))+
	     s(2)(0)*point(0)*point(1)+s(3)(0)*(one-point(0))*point(1),
	     s(0)(1)*(one-point(0))*(one-point(1))+s(1)(1)*point(0)*(one-point(1))+
	     s(2)(1)*point(0)*point(1)+s(3)(1)*(one-point(0))*point(1));
  }
  
  
  //! transformation DFi in the case of straight quadrilateral
  void QuadrangleGeomReference::DFjLinear(const VectR2& s,
                                          const R2& point, Matrix2_2& res) const
  {
    Real_wp one(1);
    res(0,0) = (one-point(1))*(s(1)(0)-s(0)(0))+point(1)*(s(2)(0)-s(3)(0));
    res(1,0) = (one-point(1))*(s(1)(1)-s(0)(1))+point(1)*(s(2)(1)-s(3)(1));
    res(0,1) = (one-point(0))*(s(3)(0)-s(0)(0))+point(0)*(s(2)(0)-s(1)(0));
    res(1,1) = (one-point(0))*(s(3)(1)-s(0)(1))+point(0)*(s(2)(1)-s(1)(1));
  }

  
  //! returns 2-D points of a point on a edge
  /*!
    \param[in] num_loc local edge number in the quad
    \param[in] t_loc abscisse on the edge
    \param[out] res local coordinates of the point in the unit square
  */
  void QuadrangleGeomReference::
  GetLocalCoordOnBoundary(int num_loc, const Real_wp& t_loc, R2& res) const
  {
    if (num_loc == 0)
      res.Init(t_loc,0);
    else if (num_loc == 1)
      res.Init(1,t_loc);
    else if (num_loc == 2)
      res.Init(Real_wp(1)-t_loc,Real_wp(1));
    else
      res.Init(0,Real_wp(1)-t_loc);
  }
  

  //! computes abscissa of a 2-D point
  /*!
    \param[in] num_loc local edge number in the quad
    \param[out] t_loc abscisse on the edge
    \param[in] pt local coordinates of the point in the unit square
  */
  void QuadrangleGeomReference::
  GetLocalCoordOnBoundary(int num_loc, const R2& pt, Real_wp& t_loc) const
  {
    if (num_loc == 0)
      t_loc = pt(0);
    else if (num_loc == 1)
      t_loc = pt(1);
    else if (num_loc == 2)
      t_loc = Real_wp(1)-pt(0);
    else
      t_loc = Real_wp(1)-pt(1);
  }

  
  //! Evaluating shape functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of shape functions on pointloc
  */
  void QuadrangleGeomReference
  ::ComputeValuesPhiNodalRef(const R2& pointloc, VectReal_wp& phi) const
  {
    int r = this->order_geom;
    phi.Reallocate((r+1)*(r+1));
    VectReal_wp phix(r+1), phiy(r+1);
    lob_geom.ComputeValuesPhiRef(pointloc(0), phix);
    lob_geom.ComputeValuesPhiRef(pointloc(1), phiy);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	phi(this->NumNodes2D(i, j)) = phix(i)*phiy(j);
    
  }
  
    
  //! Evaluating gradient of shape functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of shape functions on pointloc
  */
  void QuadrangleGeomReference
  ::ComputeGradientPhiNodalRef(const R2& pointloc, VectR2& grad_phi) const
  {
    int r = this->order_geom;
    grad_phi.Reallocate((r+1)*(r+1));
    VectReal_wp phix(r+1), phiy(r+1);
    VectReal_wp dphix(r+1), dphiy(r+1);
    for (int i = 0; i <= r; i++)
      {
	phix(i) = lob_geom.EvaluatePhi(i, pointloc(0));
	phiy(i) = lob_geom.EvaluatePhi(i, pointloc(1));
	dphix(i) = lob_geom.EvaluatePhiGrad(i, pointloc(0));
	dphiy(i) = lob_geom.EvaluatePhiGrad(i, pointloc(1));
      }
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	grad_phi(this->NumNodes2D(i, j)).Init(dphix(i)*phiy(j), phix(i)*dphiy(j));
    
  }
  
  
  //! computation of gradient of u, on nodal points
  template <class T0>
  void QuadrangleGeomReference::
  ComputeNodalGradientRefT(const Vector<T0>& Uloc, Vector<TinyVector<T0, 2> >& gradU) const
  {
    T0 zero; SetComplexZero(zero);
    int r = this->order_geom;
    gradU.Reallocate((r+1)*(r+1));
    for (int i1 = 0; i1 <= r; i1++)
      for (int i2 = 0; i2 <= r; i2++)
	{
	  int i = this->NumNodes2D(i1,i2);
	  gradU(i).Fill(zero);
	  for (int j = 0; j <= r; j++)
	    {
	      gradU(i)(0) += Uloc( this->NumNodes2D(j,i2))*lob_geom.GradPhi(j,i1);
	      gradU(i)(1) += Uloc(this->NumNodes2D(i1,j))*lob_geom.GradPhi(j,i2);
	    }
	}
  }


  //! displays informations class QuadrangleGeomReference
  ostream& operator <<(ostream& out, const QuadrangleGeomReference& e)
  {
    out << static_cast<const ElementGeomReference<Dimension2>& >(e);
    return out;
  }
}
  
#define MONTJOIE_FILE_QUADRANGLE_GEOM_REFERENCE_CXX
#endif
