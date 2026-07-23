#ifndef MONTJOIE_FILE_TRIANGLE_GEOM_REFERENCE_CXX

namespace Montjoie
{

#ifndef SELDON_WITH_COMPILED_LIBRARY
  int TriangleGeomReference::type_interpolation_nodal(0);
#endif

  //! default constructor
  TriangleGeomReference::TriangleGeomReference() : ElementGeomReference<Dimension2>()
  {
    // normales of the three edges
    this->normale.Reallocate(3);
    this->normale(0).Init(0,-1);
    this->normale(1).Init(1,1);
    this->normale(2).Init(-1,0);
    
    this->hybrid_type_elt = 0;
    this->nb_vertices_elt = 3;
    this->nb_boundaries_elt = 3;
    this->nb_edges_elt = 3;
  }
  

  size_t TriangleGeomReference::GetMemorySize() const
  {
    size_t taille = ElementGeomReference<Dimension2>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += InvWeightPolynomial.GetMemorySize();
    taille += CoefLegendre.GetMemorySize();
    taille += NumOrtho2D.GetMemorySize();
    taille += Seldon::GetMemorySize(OddJacobiPolynom);
    taille += CoefOddJacobi.GetMemorySize();
    taille += coefFi_curve.GetMemorySize();
    taille += InverseVDM.GetMemorySize();
    return taille;
  }

  
  //! construction of finite element
  void TriangleGeomReference::ConstructFiniteElement(int rgeom)
  {
    this->order_geom = rgeom;

    // orthogonal functions
    ComputeOrthogonalFunctions(rgeom);
    
    // nodal shape functions (used for determination of Fi)
    ConstructNodalShapeFunctions(rgeom);
  } 
  
  
  //! construction of nodal points on the triangle
  void TriangleGeomReference::ConstructNodalShapeFunctions(int r)
  {    
    this->edge_geom.ConstructFiniteElement(r);
    this->points_nodal_nd.Reallocate((r+1)*(r+2)/2);    
    
    // construction of nodes numbering
    MeshNumbering<Dimension2>
      ::ConstructTriangularNumbering(r, this->NumNodes2D, this->CoordinateNodes);
    
    // EdgesNode(i, num_loc) : node number of i-th node of edge num_loc
    this->EdgesNode.Reallocate(r+1, 3);
    for (int i = 0; i <= r; i++)
      this->EdgesNode(i,0) = this->NumNodes2D(i,0);
    
    for (int i = 0; i <= r; i++)
      this->EdgesNode(i,1) = this->NumNodes2D(r-i, i);
    
    for (int i = 0; i <= r; i++)
      this->EdgesNode(i,2) = this->NumNodes2D(0, r-i);	
    
    // Lobatto points as nodal points
    if (type_interpolation_nodal == REGULAR_BASIS)
      {
	this->points_nodal1d.Reallocate(r+1);
	this->points_nodal1d.Fill();
	Mlt(Real_wp(1)/r, this->points_nodal1d);
	
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= r-i; j++)
	    this->points_nodal_nd(this->NumNodes2D(i, j))
              .Init(this->points_nodal1d(i), this->points_nodal1d(j));
      }
    else
      ConstructLobattoPoints(r, type_interpolation_nodal, this->points_nodal1d,
                             this->points_nodal_nd);

    // computation of lagrangian functions (Vdm)
    ComputeLagrangianFunctions(r);
    
    VectR2 grad_phi;
    this->GradientPhi_Nodal.Reallocate(this->points_nodal_nd.GetM(),
                                       this->points_nodal_nd.GetM());
    for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
      {
	ComputeGradientPhiNodalRef(this->points_nodal_nd(j), grad_phi);
	for (int k = 0; k < this->points_nodal_nd.GetM(); k++)
	  this->GradientPhi_Nodal(k, j) = grad_phi(k);
      }

    ComputeCurvedTransformation();
    ComputeCoefficientTransformation();
  }
  
  
  //! construction of nodal points on the triangle (pascal's triangle)
  void TriangleGeomReference
  ::ConstructRegularPoints(int r, VectReal_wp& points_nodal1d_,
                           VectR2& points_nodal2d_, const Matrix<int>& NumNodes2D_)
  {
    // we construct equally spaced points on the triangle
    points_nodal1d_.Reallocate(r+1);
    for (int i = 0; i <= r; i++)
      points_nodal1d_(i) = Real_wp(i)/r;
    
    int nb_points_nodal_elt = (r+1)*(r+2)/2;
    points_nodal2d_.Reallocate(nb_points_nodal_elt);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  points_nodal2d_(NumNodes2D_(i, j))(0) = points_nodal1d_(i);
	  points_nodal2d_(NumNodes2D_(i, j))(1) = points_nodal1d_(j);
	}    
  }
  
  
  //! construction of nodal points on the triangle 
  void TriangleGeomReference
  ::ConstructLobattoPoints(int r, int type_interpolation, VectReal_wp& points_nodal1d_,
                           VectR2& points_nodal2d_)
  {
    Globatto<Real_wp> lob;
    lob.ConstructQuadrature(r, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    points_nodal1d_ = lob.Points();
	
    points_nodal2d_.Reallocate((r+1)*(r+2)/2);
    Real_wp zero(0), one(1);
    R2 center(one/Real_wp(3), one/Real_wp(3));
    
    // first, we put vertices of the triangle
    points_nodal2d_(0).Init(zero, zero);
    points_nodal2d_(1).Init(one, zero);
    points_nodal2d_(2).Init(zero, one);
    
    // then, Lobatto points on the edges
    for (int i = 0; i < r-1; i++)
      {
	points_nodal2d_(3+i).Init(points_nodal1d_(i+1), zero);
	points_nodal2d_(3+i+r-1).Init(one-points_nodal1d_(i+1), points_nodal1d_(i+1));
	points_nodal2d_(3+i+2*(r-1)).Init(zero, one-points_nodal1d_(i+1));
      }
    
    // then inside nodes
    int ind = 3*r;
    int n1 = 0, n3 = 0, n6 = 0;
    switch (r)
      {
      case 1:
      case 2:
	{
	  // no inside nodes
	}
	break;
      case 3:
	{
	  n1 = 1;
	  points_nodal2d_(ind) = center;
	}
	break;
      case 4 :
	{
	  n3 = 1;
	  points_nodal2d_(ind).Init(0.2371200168, 0.2371200168);
	}
	break;
      case 5 :
	{
	  n3 = 2;
	  points_nodal2d_(ind).Init(0.1575181512, 0.1575181512);
	  points_nodal2d_(ind+3).Init(0.410515151, 0.410515151);
	}
	break;
      case 6 :
	{
	  n1 = 1; n3 = 1; n6 = 1;
	  points_nodal2d_(ind) = center;
	  points_nodal2d_(ind+1).Init(0.1061169285, 0.1061169285);
	  points_nodal2d_(ind+4).Init(0.1332918645, 0.3097982151);
	  
	  if (type_interpolation == FEKETE_BASIS)
	    {
	      n1 = 1; n3 = 1; n6 = 1;
	      points_nodal2d_(ind) = center;
	      points_nodal2d_(ind+1).Init(0.1063354684, 0.1063354684);
	      points_nodal2d_(ind+4).Init(0.1171809171, 0.3162697959);
	    }
	  
	}
	break;
      case 7 :
	{
	  n1 = 0; n3 = 3; n6 = 1;
	  points_nodal2d_(ind).Init(0.0660520784, 0.0660520784);
	  points_nodal2d_(ind+3).Init(0.4477725053, 0.4477725053);
	  points_nodal2d_(ind+6).Init(0.2604038024, 0.2604038024);
	  points_nodal2d_(ind+9).Init(0.0914849272, 0.2325524777);
	}
	break;
      case 8 :
	{
	  n1 = 0; n3 = 3; n6 = 2;
	  points_nodal2d_(ind).Init(0.0469685351, 0.0469685351);
	  points_nodal2d_(ind+3).Init(0.2033467796, 0.2033467796);
	  points_nodal2d_(ind+6).Init(0.3905496216, 0.3905496216);
	  points_nodal2d_(ind+9).Init(0.0840385433, 0.3617970895);
	  points_nodal2d_(ind+15).Init(0.0679538347, 0.1801396087);
	}
	break;
      case 9 :
	{
	  n1 = 1; n3 = 3; n6 = 3;
	  points_nodal2d_(ind) = center;
	  points_nodal2d_(ind+1).Init(0.0355775717, 0.0355775717);
	  points_nodal2d_(ind+4).Init(0.4640303025, 0.4640303025);
	  points_nodal2d_(ind+7).Init(0.1633923069, 0.1633923069);
	  points_nodal2d_(ind+10).Init(0.0684032457, 0.2966333890);
	  points_nodal2d_(ind+16).Init(0.0529419344, 0.1439089974);
	  points_nodal2d_(ind+22).Init(0.1806100558, 0.3225890045);
	  
	  if (type_interpolation == FEKETE_BASIS)
	    {
	      points_nodal2d_(ind) = center;
	      points_nodal2d_(ind+1).Init(0.0489345696, 0.0489345696);
	      points_nodal2d_(ind+4).Init(0.4699587644, 0.4699587644);
	      points_nodal2d_(ind+7).Init(0.1704318201, 0.1704318201);
	      points_nodal2d_(ind+10).Init(0.0588564879, 0.3010242110);
	      points_nodal2d_(ind+16).Init(0.0551758079, 0.1543901944);
	      points_nodal2d_(ind+22).Init(0.1784337588, 0.3252434900);
	    }
	}
	break;
      case 10 :
	{
	  n1 = 0; n3 = 4; n6 = 4;
	  points_nodal2d_(ind).Init(0.0265250690, 0.0265250690);
	  points_nodal2d_(ind+3).Init(0.1330857076, 0.1330857076);
	  points_nodal2d_(ind+6).Init(0.4232062312, 0.4232062312);
	  points_nodal2d_(ind+9).Init(0.2833924371, 0.2833924371);
	  points_nodal2d_(ind+12).Init(0.0592706549, 0.3934913008);
	  points_nodal2d_(ind+18).Init(0.0545659823, 0.2462883939);
	  points_nodal2d_(ind+24).Init(0.0409265838, 0.1163195333);
	  points_nodal2d_(ind+30).Init(0.1481684519, 0.2707097521);
	}
	break;
      case 11 :
	{
	  n1 = 0; n3 = 5; n6 = 5;
	  points_nodal2d_(ind).Init(0.0204278105, 0.0204278105);
	  points_nodal2d_(ind+3).Init(0.4746683133, 0.4746683133);
	  points_nodal2d_(ind+6).Init(0.1104863580, 0.1104863580);
	  points_nodal2d_(ind+9).Init(0.2422268680, 0.2422268680);
	  points_nodal2d_(ind+12).Init(0.3764778430, 0.3764778430);
	  points_nodal2d_(ind+15).Init(0.0491488966, 0.3361965758);
	  points_nodal2d_(ind+21).Init(0.0442774816, 0.2073828199);
	  points_nodal2d_(ind+27).Init(0.0324707491, 0.0959885136);
	  points_nodal2d_(ind+33).Init(0.1289906484, 0.3649298091);
	  points_nodal2d_(ind+39).Init(0.1234642640, 0.2299163838);
	}
	break;
      case 12 :
	{
	  n1 = 1; n3 = 4; n6 = 7;
	  points_nodal2d_(ind) = center;
	  points_nodal2d_(ind+1).Init(0.0156346170, 0.0156346170);
	  points_nodal2d_(ind+4).Init(0.0926054449, 0.0926054449);
	  points_nodal2d_(ind+7).Init(0.4452623516, 0.4452623516);
	  points_nodal2d_(ind+10).Init(0.2091994115, 0.2091994115);
	  points_nodal2d_(ind+13).Init(0.0424612551, 0.4145270405);
	  points_nodal2d_(ind+19).Init(0.0401529071, 0.2890538328);
	  points_nodal2d_(ind+25).Init(0.0357329057, 0.1765067641);
	  points_nodal2d_(ind+31).Init(0.0257050155, 0.0799959921);
	  points_nodal2d_(ind+37).Init(0.1069206829, 0.3177076680);
	  points_nodal2d_(ind+43).Init(0.1010201598, 0.1969473706);
	  points_nodal2d_(ind+49).Init(0.2139555316, 0.3292603428);
	  
	  if (type_interpolation == FEKETE_BASIS)
	    {
	      points_nodal2d_(ind) = center;
	      points_nodal2d_(ind+1).Init(0.4005558262, 0.4005558262);
	      points_nodal2d_(ind+4).Init(0.2618405201, 0.2618405201); 
	      points_nodal2d_(ind+7).Init(0.0807386775, 0.0807386775);
	      points_nodal2d_(ind+10).Init(0.0336975736, 0.0336975736); 
	      points_nodal2d_(ind+13).Init(0.1089969290, 0.3837518758);
	      points_nodal2d_(ind+19).Init(0.1590834479, 0.2454317980);
	      points_nodal2d_(ind+25).Init(0.0887037176, 0.1697134458);
	      points_nodal2d_(ind+31).Init(0.0302317829, 0.4071849276);
	      points_nodal2d_(ind+37).Init(0.0748751152, 0.2874821712);
	      points_nodal2d_(ind+43).Init(0.0250122615, 0.2489279690);
	      points_nodal2d_(ind+49).Init(0.0262645218, 0.1206826354);
	    }
	}
	break;
      default :
	{
	  cout << "Interpolation points for r = " << r
	       << " are not implemented for triangles" << endl;
	  
	  abort();
	}
      }
    
    int nb_couples = n1 + n3 + n6;
    IVect NbPoints_couple(nb_couples);
    for (int i = 0; i < n1; i++)
      NbPoints_couple(i) = 1;
    
    for (int i = 0; i < n3; i++)
      NbPoints_couple(n1+i) = 3;

    for (int i = 0; i < n6; i++)
      NbPoints_couple(n1+n3+i) = 6;
    
    ind = 3*r;    
    // now we generate the other points by symmetry
    for (int i = 0; i < nb_couples; i++)
      {
	int np = NbPoints_couple(i);
	
	if (np==3)
	  {
	    Real_wp x = points_nodal2d_(ind)(0), y = one-2.0*x;
	    points_nodal2d_(ind).Init(x, x);
	    points_nodal2d_(ind+1).Init(x, y);
	    points_nodal2d_(ind+2).Init(y, x);
	  }
	else if (np==6)
	  {
	    Real_wp x = points_nodal2d_(ind)(0), y = points_nodal2d_(ind)(1); 
	    Real_wp t = one-x-y;
	    points_nodal2d_(ind).Init(x, y);
	    points_nodal2d_(ind+1).Init(x, t);
	    points_nodal2d_(ind+2).Init(y, x);
	    points_nodal2d_(ind+3).Init(y, t);
	    points_nodal2d_(ind+4).Init(t, x);
	    points_nodal2d_(ind+5).Init(t, y);
	  }
	
	ind += np;
      }
  }
  
  
  //! comstructing orthogonal functions generating the polynomial space P_r
  void TriangleGeomReference
  ::ComputeOrthogonalFunctions(int r)
  {
    // computation or orthogonal polynomials
    GetJacobiPolynomial(LegendrePolynom, r+1, Real_wp(0), Real_wp(0));
    OddJacobiPolynom.Reallocate(r+2);
    for (int i = 0; i <= r+1; i++)
      GetJacobiPolynomial(OddJacobiPolynom(i), r+1-i, Real_wp(2*i+1), Real_wp(0));
    
    // we compute weight of polynomials in order to obtain orthonormal expansion
    CoefLegendre.Reallocate(r+2); CoefLegendre.Fill(0);
    CoefOddJacobi.Reallocate(r+2, r+2); CoefOddJacobi.Fill(0);
    VectReal_wp xi, omega;
    ComputeGaussLegendre(xi, omega, r+1);
    for (int q = 0; q <= r+1; q++)
      {
	VectReal_wp Pn;
	EvaluateJacobiPolynomial(LegendrePolynom, r+1, 2.0*xi(q) - 1.0, Pn);
	for (int i = 0; i <= r+1; i++)
	  CoefLegendre(i) += omega(q)*Pn(i)*Pn(i);
	
	for (int i = 0; i <= r+1; i++)
	  {
	    EvaluateJacobiPolynomial(OddJacobiPolynom(i), r+1-i, 2.0*xi(q) - 1.0, Pn);
	    for (int j = 0; j <= r+1-i; j++)
	      CoefOddJacobi(i, j) += omega(q)*Pn(j)*Pn(j)*pow(1.0-xi(q), 2*i+1); 
	  }
      }
    
    // we keep 1/sqrt(omega)
    for (int i = 0; i <= r; i++)
      {
	CoefLegendre(i) = 1.0/sqrt(CoefLegendre(i));
	for (int j = 0; j <= r-i; j++)
	  CoefOddJacobi(i, j) = 1.0/sqrt(CoefOddJacobi(i, j));
      }
    
    // orthonormal functions
    InvWeightPolynomial.Reallocate((r+1)*(r+2)/2);
    NumOrtho2D.Reallocate(r+1, r+1); NumOrtho2D.Fill(-1);
    // we sort by increasing sum = i+j
    int nb = 0;
    for (int sum = 0; sum <= r; sum++)
      for (int i = 0; i <= sum; i++)
        {
          int j = sum-i;
	  NumOrtho2D(i, j) = nb;
	  InvWeightPolynomial(nb) = CoefLegendre(i)*CoefOddJacobi(i, j);
	  nb++;
	}

    //DISP(NumOrtho2D);
  }
  
  
  //! Computation of Lagrangian functions
  void TriangleGeomReference::ComputeLagrangianFunctions(int r)
  {
    // 1-D nodal functions
    lob_geom.AffectPoints(this->points_nodal1d);
    
    // 2-D nodal functions
    int nb_points_nodal_elt = this->points_nodal_nd.GetM();
    InverseVDM.Reallocate(nb_points_nodal_elt, nb_points_nodal_elt);
    InverseVDM.Fill(Real_wp(0));
    
    // evaluation of orthogonal polynomials
    VectReal_wp phi;
    for (int k = 0; k < nb_points_nodal_elt; k++)
      {
	ComputeValuesPhiOrthoRef(r, NumOrtho2D, InvWeightPolynomial,
                                 this->points_nodal_nd(k), phi);
        
	for (int j = 0; j < nb_points_nodal_elt; j++)
	  InverseVDM(j, k) = phi(j);
      }

    //InverseVDM.WriteText("VDM.dat");
    GetInverse(InverseVDM);
  }
  
  
  //! computation of precomputed coefficients for curved elements
  void TriangleGeomReference::ComputeCurvedTransformation()
  {
    // interpolation using Gauss-Lobatto points for curved edges
    coefFi_curve.Reallocate(this->points_nodal_nd.GetM(), this->points_nodal_nd.GetM());
    coefFi_curve.Fill(0);
    // loop on each intern nodal point (strictly inside the triangle)
    for (int i = 3*this->order_geom; i < this->points_nodal_nd.GetM(); i++)
      {
	// for the definition of transformation, see the book of P. Solin, 2003
	Real_wp x = this->points_nodal_nd(i)(0);
	Real_wp y = this->points_nodal_nd(i)(1);
	
	// part due to vertices
	coefFi_curve(i, 0) = 1.0 - x - y;
	coefFi_curve(i, 1) = x;
	coefFi_curve(i, 2) = y;
	
	// part due to first edge
	Real_wp zeta = x + 0.5*y;
	Real_wp coef = (1.0-x-y)*x/(zeta*(1.0-zeta));
	// contribution of points of the edge
	for (int k = 0; k <= this->order_geom; k++)
	  coefFi_curve(i, this->NumNodes2D(k, 0)) += lob_geom.EvaluatePhi(k, zeta)*coef;
	
	// contribution of the extremities of the edge
	coefFi_curve(i, 0) -= (1.0-zeta)*coef;
	coefFi_curve(i, 1) -= zeta*coef;
	
	// part due to second edge
	zeta = 0.5*y - 0.5*x + 0.5;
	coef = y*x/(zeta*(1.0-zeta));
	// contribution of points of the edge
	for (int k = 0; k <= this->order_geom; k++)
	  coefFi_curve(i, this->NumNodes2D(this->order_geom-k, k))
            += lob_geom.EvaluatePhi(k, zeta)*coef;
	
	// contribution of the extremities of the edge
	coefFi_curve(i, 1) -= (1.0-zeta)*coef;
	coefFi_curve(i, 2) -= zeta*coef;
	
	// part due to third edge
	zeta = 1.0 - 0.5*x - y;
	coef = y*(1.0-x-y)/(zeta*(1.0-zeta));
	// contribution of points of the edge
	for (int k = 0; k <= this->order_geom; k++)
	  coefFi_curve(i, this->NumNodes2D(0, this->order_geom-k))
            += lob_geom.EvaluatePhi(k, zeta)*coef;
	
	// contribution of the extremities of the edge
	coefFi_curve(i, 2) -= (1.0-zeta)*coef;
	coefFi_curve(i, 0) -= zeta*coef;

      }

  }
  
  
  //! computation of coefficients for a fast evaluation of Fi
  void TriangleGeomReference::ComputeCoefficientTransformation()
  {
    VectR2 all_points = this->points_nodal_nd;
    all_points.PushBack(this->points2d);
    all_points.PushBack(this->points_dof2d);
    coefFi_curve.Resize(all_points.GetM(), this->points_nodal_nd.GetM());
    
    // for the other points, we are using nodal basis functions
    VectReal_wp phi;
    for (int i = this->points_nodal_nd.GetM(); i < all_points.GetM(); i++)
      {
	ComputeValuesPhiNodalRef(all_points(i), phi);
	for (int k = 0; k < this->points_nodal_nd.GetM(); k++)
	  coefFi_curve(i, k) = phi(k);
      } 
    
  }
  

  //! sets nodal points
  void TriangleGeomReference::SetNodalPoints(int r, const VectR2& nodal_points)
  {
    this->order_geom = r;
    int nb_nodes = (r+1)*(r+2)/2;
    if (nodal_points.GetM() != nb_nodes)
      {
        cout << "Incorrect number of points " << nodal_points.GetM() << " instead of " << nb_nodes << endl;
        abort();
      }

    if (NumOrtho2D.GetM() <= r)
      {
        cout << "Not enough orthogonal polynomials : " << NumOrtho2D.GetM() << " <= " << r << endl;
        abort();
      }

    // points on vertices ?
    int nb = 0; Vector<bool> node_used(nodal_points.GetM());
    node_used.Fill(false);
    this->points_nodal_nd.Reallocate(nodal_points.GetM());
    for (int i = 0; i < nodal_points.GetM(); i++)
      {
        if (nodal_points(i) == R2(0, 0))
          {
            node_used(i) = true;
            this->points_nodal_nd(0) = R2(0, 0);
            nb++;
          }

        if (nodal_points(i) == R2(1, 0))
          {
            node_used(i) = true;
            this->points_nodal_nd(1) = R2(1, 0);
            nb++;
          }

        if (nodal_points(i) == R2(0, 1))
          {
            node_used(i) = true;
            this->points_nodal_nd(2) = R2(0, 1);
            nb++;
          }
      }

    if ((nb != 3) && (nb != 0))
      {
        cout << "You have to put all the vertices or none, nb = " << nb << endl;
        abort();
      }
    
    // we check the points on the edges
    VectReal_wp edge1, edge2, edge3;
    for (int i = 0; i < nodal_points.GetM(); i++)
      if (!node_used(i))
        {
          if (abs(nodal_points(i)(1)) <= R2::threshold)
            {
              node_used(i) = true;
              edge1.PushBack(nodal_points(i)(0));
            }
          else if (abs(nodal_points(i)(0)) <= R2::threshold)
            {
              node_used(i) = true;
              edge3.PushBack(nodal_points(i)(1));
            }
          else if (abs(Real_wp(1)-nodal_points(i)(0) - nodal_points(i)(1)) <= R2::threshold)
            {
              node_used(i) = true;
              edge2.PushBack(nodal_points(i)(1));
            }
        }
    
    Sort(edge1); Sort(edge2); Sort(edge3);
    if ((edge1.GetM() != edge2.GetM()) || (edge1.GetM() != edge3.GetM()))
      {
        cout << "Points on edge not invariant by rotation " << endl;
        abort();
      }

    bool pt_sym_edge = true;
    for (int i = 0; i < edge1.GetM(); i++)
      {
        if (abs(1.0 - edge1(i) - edge1(edge1.GetM()-1-i)) > R2::threshold)
          pt_sym_edge = false;

        if (abs(edge1(i)-edge2(i)) > R2::threshold)
          pt_sym_edge = false;
        
        if (abs(edge1(i)-edge3(i)) > R2::threshold)
          pt_sym_edge = false;        
      }

    if (!pt_sym_edge)
      {
        cout << "Points on edge not invariant by rotation " << endl;
        cout << edge1 << ", " << edge2 << ", " << edge3 << endl;
        abort();
      }

    //  we do not take into account points edges for 1-D finite element
    // to be modified
    this->edge_geom.ConstructFiniteElement(r);

    // constructing EdgesNode and NumNodes2D (if vertices are included)
    bool conforming = true;
    if (nb == 3)
      {
        if (edge1.GetM() != r-1)
          {
            cout << "Incorrect number of points on the edge" << endl;
            abort();
          }

        this->NumNodes2D.Reallocate(r+1, r+1); this->NumNodes2D.Fill(-1);
        this->NumNodes2D(0, 0) = 0; this->NumNodes2D(r, 0) = 1; this->NumNodes2D(0, r) = 2;
        this->points_nodal1d.Reallocate(r+1);
        this->points_nodal1d(0) = Real_wp(0);
        this->points_nodal1d(r) = Real_wp(1);
        for (int i = 1; i < r; i++)
          {
            this->points_nodal1d(0) = edge1(i-1);
            
            this->points_nodal_nd(2+i) = R2(edge1(i-1), 0);
            this->NumNodes2D(i, 0) = 2+i;

            this->points_nodal_nd(r+i+1) = R2(Real_wp(1)-edge1(i-1), edge1(i-1));
            this->NumNodes2D(r-i, i) = r+i+1;

            this->points_nodal_nd(2*r+i) = R2(Real_wp(0), Real_wp(1)-edge1(i-1));
            this->NumNodes2D(0, r-i) = 2*r+i;
          }
        
        this->EdgesNode.Reallocate(r+1, 3);
        for (int i = 0; i <= r; i++)
          this->EdgesNode(i,0) = this->NumNodes2D(i,0);
    
        for (int i = 0; i <= r; i++)
          this->EdgesNode(i,1) = this->NumNodes2D(r-i, i);
        
        for (int i = 0; i <= r; i++)
          this->EdgesNode(i,2) = this->NumNodes2D(0, r-i);	
    
        nb += 3*(r-1);
      }
    else
      {
        // nonconforming shape functions
        conforming = false;
        this->points_nodal1d = edge1;
        this->NumNodes2D.Clear(); this->EdgesNode.Clear();

        if (edge1.GetM() != 0)
          {
            if (edge1.GetM() != r)
              {
                cout << "Incorrect number of points on the edge : = " << edge1.GetM() << " different from " << r << endl;
                abort();
              }
            
            for (int i = 0; i < edge1.GetM(); i++)
              {
                this->points_nodal_nd(i) = R2(edge1(i), 0);
                this->points_nodal_nd(r+i) = R2(Real_wp(1)-edge1(i), edge1(i));
                this->points_nodal_nd(2*r+i) = R2(Real_wp(0), Real_wp(1)-edge1(i));
              }
            
            nb += 3*r;
          }
      }

    // and finally interior points
    for (int i = 0; i < nodal_points.GetM(); i++)
      if (!node_used(i))
        {
          node_used(i) = true;
          this->points_nodal_nd(nb) = nodal_points(i);
          nb++;
        }

    if (nb != nodal_points.GetM())
      {
        cout << "Not possible" << endl;
        abort();
      }

    // we recompute shape functions
    ComputeLagrangianFunctions(r);
    
    VectR2 grad_phi;
    this->GradientPhi_Nodal.Reallocate(this->points_nodal_nd.GetM(),
                                       this->points_nodal_nd.GetM());
    for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
      {
	ComputeGradientPhiNodalRef(this->points_nodal_nd(j), grad_phi);
	for (int k = 0; k < this->points_nodal_nd.GetM(); k++)
	  this->GradientPhi_Nodal(k, j) = grad_phi(k);
      }
    
    if (conforming)
      {
        ComputeCurvedTransformation();
        ComputeCoefficientTransformation();
      }
  }
  
  
  /************************
   * Convenient functions *
   ************************/
  
  
  //! evaluating orthogonal functions phi at pointloc
  void TriangleGeomReference
  ::ComputeValuesPhiOrthoRef(int r, const Matrix<int>& NumOrtho, const VectReal_wp& InvWeight,
                             const R2& pointloc, VectReal_wp& phi) const
  {
    int nb_fct = (r+2)*(r+1)/2;
    
    // orthogonal polynomials :
    // L_i(2x/(1-y) - 1) P_j^{2i+1}(2y-1) (1-y)^i
    // = L_i(xtilde) P_j^(2i+1)(ytilde) ((1-ytilde)/2)^i
    // where (xtilde, ytilde) is on the square [-1,1]^2 and (x,y) on the unit triangle
    Real_wp xtilde(0), ytilde(1);    
    if (abs(pointloc(1) - 1.0) > epsilon_machine)
      {
	xtilde = 2.0*pointloc(0)/(1.0-pointloc(1)) - 1.0;
	ytilde = 2.0*pointloc(1) - 1.0;
      }
    
    VectReal_wp Li, Pj; phi.Reallocate(nb_fct);
    EvaluateJacobiPolynomial(LegendrePolynom, r, xtilde, Li); 
    Real_wp pow_oneMinusytilde = 1.0, oneMinusy = 0.5*(1.0-ytilde), vali;
    for (int i = 0; i <= r; i++)
      {
        if (i == 0)
          vali = 1.0;
        else if (i == 1)
          vali = 2.0*pointloc(0) + pointloc(1) -1.0;
        else
          vali = pow_oneMinusytilde*Li(i);
        
	EvaluateJacobiPolynomial(OddJacobiPolynom(i), r-i, ytilde, Pj);
	for (int j = 0; j <= r-i; j++)
	  phi(NumOrtho(i, j)) = Pj(j)*vali;
	
	pow_oneMinusytilde *= oneMinusy;
      }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < phi.GetM(); i++)
      phi(i) *= InvWeight(i);
    
  }
  

  //! evaluating gradient of orthogonal functions phi at pointloc  
  void TriangleGeomReference
  ::ComputeGradientPhiOrthoRef(int r, const Matrix<int>& NumOrtho,
                               const VectReal_wp& InvWeight,
                               const R2& pointloc, VectR2& grad_phi) const
  {
    int nb_fct = (r+2)*(r+1)/2;
    
    // orthogonal polynomials :
    // L_i(2x/(1-y) - 1) P_j^{2i+1}(2y-1) (1-y)^i
    // = L_i(xtilde) P_j^(2i+1)(ytilde) ((1-ytilde)/2)^i
    // where (xtilde, ytilde) is on the square [-1,1]^2 and (x,y) on the unit triangle
    
    // derivative with respect to xtilde, ytilde is equal to :
    // (L'_i(xtilde) P_j^(2i+1)(ytilde) ((1-ytilde)/2)^i; L_i(xtilde)
    // ((1-ytilde)/2)^{i-1} (P'_j^(2i+1)(ytilde) ((1-ytilde)/2) - 1/2 i P_j^(2i+1)(ytilde) )
    Real_wp xtilde(0), ytilde(1), dxtilde_dx(0), dxtilde_dy(0);
    if (abs(pointloc(1) - 1.0) > epsilon_machine)
      {
        xtilde = 2.0*pointloc(0)/(1.0-pointloc(1)) - 1.0;
        dxtilde_dx = 2.0/(1.0-pointloc(1));
        dxtilde_dy = dxtilde_dx*pointloc(0)/(1.0-pointloc(1));
      }
	
    ytilde = 2.0*pointloc(1) - 1.0;
        
    VectReal_wp Li, Pj, dLi, dPj;
    grad_phi.Reallocate(nb_fct);
    EvaluateJacobiPolynomial(LegendrePolynom, r, xtilde, Li, dLi);
    Real_wp powm1_oneMinusytilde = 0.0, pow_oneMinusytilde = 1.0, oneMinusy = 1.0-pointloc(1);
    Real_wp vali; R2 gradi;
    for (int i = 0; i <= r; i++)
      {
        // derivating L_i(xtilde) (1-y)
        if (i == 0)
          {
            // derivative of 1 is 0
            vali = 1.0;
            gradi.Fill(0);
          }
        else if (i == 1)
          {
            // derivative of xtilde (1-y) = 2x -1 + y  is (2, 1)
            vali = 2.0*pointloc(0) + pointloc(1) - 1.0;
            gradi.Init(2.0, 1.0);
          }
        else
          {
            // in that case, we can use dxtilde_dx and dxtilde_dy
            vali = pow_oneMinusytilde*Li(i);
            gradi(0) = pow_oneMinusytilde*dLi(i)*dxtilde_dx;
            gradi(1) = -i*powm1_oneMinusytilde*Li(i) + pow_oneMinusytilde*dLi(i)*dxtilde_dy;
          }
        
        // taking into account P_j^(2i+1)(ytilde)
	EvaluateJacobiPolynomial(OddJacobiPolynom(i), r-i, ytilde, Pj, dPj);
	for (int j = 0; j <= r-i; j++)
          grad_phi(NumOrtho(i, j)).Init(Pj(j)*gradi(0), 2.0*vali*dPj(j) + gradi(1)*Pj(j));
          	
	powm1_oneMinusytilde = pow_oneMinusytilde;
	pow_oneMinusytilde *= oneMinusy;
      }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < grad_phi.GetM(); i++)
      grad_phi(i) *= InvWeight(i);

  }
    
  
  //! Evaluating shape functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of shape functions on pointloc
  */
  void TriangleGeomReference
  ::ComputeValuesPhiNodalRef(const R2& pointloc, VectReal_wp& phi) const
  {
    VectReal_wp psi;
    phi.Reallocate(this->points_nodal_nd.GetM());
    phi.Fill(0);
    ComputeValuesPhiOrthoRef(this->order_geom, NumOrtho2D,
                             InvWeightPolynomial, pointloc, psi);
    Mlt(InverseVDM, psi, phi);
  }
  
  
  //! Evaluating gradient of shape functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of shape functions on pointloc
  */
  void TriangleGeomReference
  ::ComputeGradientPhiNodalRef(const R2& pointloc, VectR2& grad_phi) const
  {
    VectReal_wp psi, phi; VectR2 grad_psi; 
    psi.Reallocate(this->points_nodal_nd.GetM()); psi.Fill(0);
    phi.Reallocate(this->points_nodal_nd.GetM()); phi.Fill(0);
    grad_phi.Reallocate(this->points_nodal_nd.GetM());
    ComputeGradientPhiOrthoRef(this->order_geom, NumOrtho2D,
                               InvWeightPolynomial, pointloc, grad_psi);
    
    for (int i = 0; i < psi.GetM(); i++)
      psi(i) = grad_psi(i)(0);
    
    Mlt(InverseVDM, psi, phi);
    
    for (int i = 0; i < psi.GetM(); i++)
      {
        grad_phi(i)(0) = phi(i);
        psi(i) = grad_psi(i)(1);
      }
    
    Mlt(InverseVDM, psi, phi);
    for (int i = 0; i < psi.GetM(); i++)
      grad_phi(i)(1) = phi(i);
  }

  //! returns the minimal length of the element
  Real_wp TriangleGeomReference
  ::GetMinimalSize(const VectR2& s) const
  {
    Real_wp h = Distance(s(0), s(1));
    h = min(h, Distance(s(1), s(2)));
    h = min(h, Distance(s(0), s(2)));
    return h;
  }
  
  
  //! returns true if the point is outside the unit triangle
  /*!
    \param[in] Xn local coordinates of the point
    \param[in] epsilon threshold 
    \return true if the point is outside
  */
  bool TriangleGeomReference
  ::OutsideReferenceElement(const VectR2& s, const R2& Xn, const Real_wp& epsilon) const
  {
    Real_wp Rmax = max(abs(s(0)(0)), abs(s(1)(0))); Rmax = max(Rmax, abs(s(2)(0)));
    Rmax = max(Rmax, abs(s(0)(1))); Rmax = max(Rmax, abs(s(1)(1))); Rmax = max(Rmax, abs(s(2)(1)));
    Real_wp dx = s(0).Distance(s(1)), dy = s(0).Distance(s(2));
    if ((Rmax == Real_wp(0)) || (dx == Real_wp(0)) || (dy == Real_wp(0)))
      return true;
    
    if (dx*Xn(0) < -epsilon*Rmax)
      return true;

    if (dy*Xn(1) < -epsilon*Rmax)
      return true;

    if (min(dx, dy)*(Xn(0) + Xn(1) - Real_wp(1)) > epsilon*Rmax)
      return true;
    
    return false;
  }
  
  
  //! returns distance of the local point to the boundary of the element
  /*!
    \param[in] pointloc local coordinates of the point
    \return distance to the boundary (negative if the point is outside)
  */
  Real_wp TriangleGeomReference::GetDistanceToBoundary(const R2& pointloc) const
  {
    // we take distance = min(\lambda_1,\lambda_2,\lambda_3)
    // where (\lambda_i)  are barycentric coordinates
    Real_wp distance = pointloc(0);
    distance = min(distance,pointloc(1));
    Real_wp dist2 = Real_wp(1)-pointloc(0)-pointloc(1);
    distance = min(distance, dist2);
    return distance;
  }
  
  
  //! returns distance of the local point to a given boundary of the element
  Real_wp TriangleGeomReference::GetDistanceToBoundary(const R2& pointloc, int num_loc) const
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
	  distance = abs(pointloc(0)+pointloc(1)-Real_wp(1));
	  distance = max(distance, -pointloc(0));
	  distance = max(distance, pointloc(0)-Real_wp(1));
	}
	break;
      case 2:
	{
	  distance = abs(pointloc(0));
	  distance = max(distance, -pointloc(1));
	  distance = max(distance, pointloc(1)-Real_wp(1));
	}
	break;	
      }
    
    return distance;
  }

  
  //! project the local point on the boundary, if outside the unit triangle
  int TriangleGeomReference
  ::ProjectPointOnBoundary(R2& pointloc) const
  {
    pointloc(0) = max(Real_wp(0), pointloc(0));
    pointloc(1) = max(Real_wp(0), pointloc(1));
    
    if (pointloc(0) < pointloc(1))
      {
	Real_wp dist2 = Real_wp(1)-pointloc(1);
	pointloc(0) = min(dist2, pointloc(0));
      }
    else
      {
	Real_wp dist2 = Real_wp(1)-pointloc(0);
	pointloc(1) = min(dist2, pointloc(1));
      }
    
    if (abs(pointloc(1)) <= R2::threshold)
      return 0;

    if (abs(pointloc(0)) <= R2::threshold)
      return 2;

    if (abs(pointloc(0)+pointloc(1)-1.0) <= R2::threshold)
      return 1;
    
    return -1;
  }
    
  
  /****************
   * Fj transform *
   ****************/
  
  
  //! transformation Fi in the case of curved triangle
  void TriangleGeomReference
  ::FjCurve(const SetPoints<Dimension2>& PTReel, const R2& pointloc, R2& res) const
  {
    res.Fill(0);
    VectReal_wp phi; ComputeValuesPhiNodalRef(pointloc, phi);
    for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
      Add(phi(node), PTReel.GetPointNodal(node), res);
  }
  
  
  //! transformation DFi in the case of curved triangle
  void TriangleGeomReference
  ::DFjCurve(const SetPoints<Dimension2>& PTReel, const R2& pointloc, Matrix2_2& res) const
  {
    res.Fill(0); R2 grad;
    VectR2 grad_phi; ComputeGradientPhiNodalRef(pointloc, grad_phi);
    for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
      {
	grad = grad_phi(node);
	res(0,0) += PTReel.GetPointNodal(node)(0)*grad(0);
	res(1,0) += PTReel.GetPointNodal(node)(1)*grad(0);
	res(0,1) += PTReel.GetPointNodal(node)(0)*grad(1);
	res(1,1) += PTReel.GetPointNodal(node)(1)*grad(1);
      }
  }


  /**********************
   * FjElem and DFjElem *
   **********************/
  
  
  //! computes res = Fi(point) for all reference points
  /*!
    \param[in] s list of vertices of the triangle
    \param[out] res "reference points" after transformation Fi
    \param[in] mesh triangular mesh considered
    \param[in] nquad triangular element number
  */
  void TriangleGeomReference::
  FjElem(const VectR2& s, SetPoints<Dimension2>& res,
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
  void TriangleGeomReference::
  FjElemNodal(const VectR2& s, SetPoints<Dimension2>& res,
	      const Mesh<Dimension2>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
    else
      FjElemNodalLinear(s, res);
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! computes res = Fi(point) for all nodal points
  void TriangleGeomReference::
  FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
	      const Mesh<Dimension3>& mesh, int i, const Face<Dimension3>& f) const
  {
    if (mesh.IsBoundaryCurved(i) > 0)
      FjElemNodalCurve(s, res, mesh, i, f);
    else
      FjElemNodalLinear(s, res);
  }


  //! computes res = Fi(point) for all quadrature points
  void TriangleGeomReference::
  FjElemQuadrature(const VectR3& s, SetPoints<Dimension3>& res,
		   const Mesh<Dimension3>& mesh, int i, const Face<Dimension3>& f) const
  {
    if (mesh.IsBoundaryCurved(i) > 0)
      {
	FjElemNodalCurve(s, res, mesh, i, f);
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  FjElemQuadratureCurve(s, res, mesh, i);
      }
    else
      {
	FjElemNodalLinear(s, res);
	FjElemQuadratureLinear(s, res);
      }
  }

  
  //! computes dF/dx and dF/dy on nodal points
  void TriangleGeomReference
  ::DFjElemNodal(const VectR3& pts_nodal, bool is_curved,
		 VectR3& dF_dx, VectR3& dF_dy) const
  {
    dF_dx.Reallocate(this->points_nodal_nd.GetM());
    dF_dy.Reallocate(this->points_nodal_nd.GetM());
    
    if (is_curved)
      {
	R2 grad;
	for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
	  {
	    dF_dx(i).Fill(0);
	    dF_dy(i).Fill(0);
	    for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	      {
		grad = this->GradientPhi_Nodal(node, i);
		Add(grad(0), pts_nodal(i), dF_dx(i));
		Add(grad(1), pts_nodal(i), dF_dy(i));
	      }
	  }
      } 
    else
      {
	R3 vec_u, vec_v;
	vec_u = pts_nodal(1) - pts_nodal(0);
	vec_v = pts_nodal(2) - pts_nodal(0);
	
	dF_dx.Fill(vec_u);
	dF_dy.Fill(vec_v);
      }
  }


  //! computes dF/dx and dF/dy on quadrature points
  void TriangleGeomReference
  ::DFjElemQuadrature(const VectR3& pts_nodal, bool is_curved,
		      const VectR3& dF_dx_node, const VectR3& dF_dy_node,
		      VectR3& dF_dx, VectR3& dF_dy) const
  {
    dF_dx.Reallocate(this->points2d.GetM());
    dF_dy.Reallocate(this->points2d.GetM());
    
    if (is_curved)
      {
	Real_wp coef;
	for (int i = 0; i < this->points2d.GetM(); i++)
	  {
	    dF_dx(i).Fill(0);
	    dF_dy(i).Fill(0);
	    int N = this->points_nodal_nd.GetM();
	    for (int node = 0; node < N; node++)
	      {
		coef = coefFi_curve(N+i, node);
		Add(coef, dF_dx_node(node), dF_dx(i));
		Add(coef, dF_dy_node(node), dF_dy(i));
	      }
	  }
      } 
    else
      {
	R3 vec_u, vec_v;
	vec_u = pts_nodal(1) - pts_nodal(0);
	vec_v = pts_nodal(2) - pts_nodal(0);
	
	dF_dx.Fill(vec_u);
	dF_dy.Fill(vec_v);
      }
  }
#endif
  
  
  //! computes res = Fi(point) for all quadrature points
  void TriangleGeomReference::
  FjElemQuadrature(const VectR2& s, SetPoints<Dimension2>& res,
		   const Mesh<Dimension2>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  FjElemQuadratureCurve(s, res, mesh, nquad);
      }
    else
      FjElemQuadratureLinear(s,res);
  }
  
  
  //! computes res = Fi(point) for all dof points
  void TriangleGeomReference::
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
    \param[in] s list of vertices of the triangle
    \param[in] PTReel "reference points" after transformation Fi
    \param[out] res jacobian matrices DFi
    \param[in] mesh triangular mesh considered
    \param[in] nquad triangular element number
  */
  void TriangleGeomReference::
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
  void TriangleGeomReference::
  DFjElemNodal(const VectR2& s, const SetPoints<Dimension2>& PTReel,
	       SetMatrices<Dimension2>& res,
	       const Mesh<Dimension2>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
    else
      DFjElemNodalLinear(s, res);
  }
  
  
  //! computes res = DFi(point) for all quadrature points
  void TriangleGeomReference::
  DFjElemQuadrature(const VectR2& s, const SetPoints<Dimension2>& PTReel,
		    SetMatrices<Dimension2>& res,
		    const Mesh<Dimension2>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  DFjElemQuadratureCurve(s, PTReel, res, mesh, nquad);
      }
    else
      DFjElemQuadratureLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all dof points
  void TriangleGeomReference::
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
  
  
  //! nodal points in the case of straight triangle  
  template<class DimensionB> void TriangleGeomReference::
  FjElemNodalLinear(const typename DimensionB::VectR_N& s, SetPoints<DimensionB>& res) const
  {
    Real_wp x,y,l1;
    typename DimensionB::R_N res_n;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
      {
	x = this->points_nodal_nd(i)(0);
	y = this->points_nodal_nd(i)(1);
	l1 = Real_wp(1)-x-y;
	res_n.Zero();
	Add(l1, s(0), res_n);
	Add(x, s(1), res_n);
	Add(y, s(2), res_n);
	res.SetPointNodal(i,res_n);
      }
  }
  
  
  //! quadrature points in the case of straight triangle
  void TriangleGeomReference::
  FjElemQuadratureLinear(const VectR2& s, SetPoints<Dimension2>& res) const
  {
    Real_wp x,y,l1;
    R2 res_n;
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	x = this->points2d(i)(0);
	y = this->points2d(i)(1);
	l1 = Real_wp(1)-x-y;
	res_n(0) = s(0)(0)*l1 + s(1)(0)*x + s(2)(0)*y;
	res_n(1) = s(0)(1)*l1 + s(1)(1)*x + s(2)(1)*y;
	res.SetPointQuadrature(i, res_n);
      }
  }


  //! quadrature points in the case of straight triangle
  void TriangleGeomReference::
  FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    Real_wp x,y,l1;
    R3 res_n;
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	x = this->points2d(i)(0);
	y = this->points2d(i)(1);
	l1 = Real_wp(1)-x-y;
	res_n(0) = s(0)(0)*l1 + s(1)(0)*x + s(2)(0)*y;
	res_n(1) = s(0)(1)*l1 + s(1)(1)*x + s(2)(1)*y;
	res_n(2) = s(0)(2)*l1 + s(1)(2)*x + s(2)(2)*y;
	res.SetPointQuadrature(i, res_n);
      }
  }

  
  //! dof points in the case of straight triangle
  void TriangleGeomReference::
  FjElemDofLinear(const VectR2& s, SetPoints<Dimension2>& res) const
  {
    Real_wp x,y,l1;
    R2 res_n;
    res.ReallocatePointsDof(this->points_dof2d.GetM());
    for (int i = 0; i < this->points_dof2d.GetM();i++)
      {
	x = this->points_dof2d(i)(0);
	y = this->points_dof2d(i)(1);
	l1 = Real_wp(1) - x - y;
	res_n(0) = s(0)(0)*l1+s(1)(0)*x+s(2)(0)*y;
	res_n(1) = s(0)(1)*l1+s(1)(1)*x+s(2)(1)*y;
	res.SetPointDof(i,res_n);
      }
  }
  
  
  //! transformation DFi for nodal points in the case of straight triangle
  void TriangleGeomReference::
  DFjElemNodalLinear(const VectR2& s, SetMatrices<Dimension2>& res) const
  {
    Matrix2_2 res_n;
    res_n(0, 0) = s(1)(0)-s(0)(0);
    res_n(0, 1) = s(2)(0)-s(0)(0);
    res_n(1, 0) = s(1)(1)-s(0)(1);
    res_n(1, 1) = s(2)(1)-s(0)(1);
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    res.FillNodal(res_n);
  }
  
  
  //! transformation DFi for quadrature points in the case of straight triangle
  void TriangleGeomReference::
  DFjElemQuadratureLinear(const VectR2& s, SetMatrices<Dimension2>& res) const
  {
    Matrix2_2 res_n;
    res_n(0, 0) = s(1)(0)-s(0)(0);
    res_n(0, 1) = s(2)(0)-s(0)(0);
    res_n(1, 0) = s(1)(1)-s(0)(1);
    res_n(1, 1) = s(2)(1)-s(0)(1);
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    res.FillQuadrature(res_n);
  }
  
  
  //! transformation DFi for dof points in the case of straight triangle
  void TriangleGeomReference::
  DFjElemDofLinear(const VectR2& s, SetMatrices<Dimension2>& res) const
  {
    Matrix2_2 res_n;
    res_n(0, 0) = s(1)(0)-s(0)(0);
    res_n(0, 1) = s(2)(0)-s(0)(0);
    res_n(1, 0) = s(1)(1)-s(0)(1);
    res_n(1, 1) = s(2)(1)-s(0)(1);
    res.ReallocatePointsDof(this->points_dof2d.GetM());
    res.FillDof(res_n);
  }

  
  //! transformation DFi for boundary points
  void TriangleGeomReference::
  FjSurfaceElemLinear(const VectR2& s, SetPoints<Dimension2>& res, int ne) const
  {
    R2 ptA = s(ne);  R2 ptB = s((ne+1)%3);
    R2 resn;
    res.ReallocatePointsQuadratureBoundary(this->points1d.GetM());
    for (int k = 0; k < this->points1d.GetM(); k++)
      {
	resn = ptA; Mlt((Real_wp(1)-this->points1d(k)), resn);
	Add(this->points1d(k), ptB, resn);
	res.SetPointQuadratureBoundary(k, resn);
      }
  }
  
  
  //! transformation Fi for nodal points in the case of curved triangle
  /*!
    \param[in] s four vertices of the triangle
    \param[out] res references points after the transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad triangle number in the mesh
    This method uses curved triangles as described in Solin
  */
  template<class DimensionB, class TypeEntity>
  void TriangleGeomReference::
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
    typename DimensionB::R_N res_n;
    
    if ((this->CoefCurve1D.GetM() != this->order_geom-1)||(this->CoefCurve1D.GetN() != r+1))
      {
        DISP(this->CoefCurve1D.GetM()); DISP(this->order_geom);
	cout<<"Did you call SetMesh after ConstructFiniteElement ?"<<endl;
	abort();
      }
    
    // edges
    for (int j = 0; j < 3; j++)
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
	    typename DimensionB::R_N ptB = s((j+1)%3);	    
	    for (int k = 1; k < this->order_geom; k++)
	      {
		res_n.Zero();
		// linear interpolation
		Add(1.0-this->PointsNodal1D(k), ptA, res_n); 
		Add(this->PointsNodal1D(k), ptB, res_n);
		res.SetPointNodal(this->EdgesNode(k, j), res_n); 
	      }
	  }
	else
	  {
	    // first step : we get the points on the curved edge
	    Vector<typename DimensionB::R_N> PtsEdge(r+1);
	    PtsEdge(0) = s(j);
	    PtsEdge(r) = s((j+1)%3);
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
    for (int node = 3*this->order_geom; node < this->points_nodal_nd.GetM(); node++)
      {
	res_n.Zero();
	for (int k = 0; k < 3*this->order_geom; k++)
	  Add(coefFi_curve(node, k), res.GetPointNodal(k), res_n);
	
	res.SetPointNodal(node, res_n);
      }
    
    // DISP(s); DISP(nquad); DISP(res);
  }
  
  
  //! transformation Fi for quadrature points in the case of curved triangle
  void TriangleGeomReference::
  FjElemQuadratureCurve(const VectR2& s, SetPoints<Dimension2>& res,
			const Mesh<Dimension2>& mesh, int nquad) const
  {
    // isoparametric element
    // other points are computed by interpolation based on nodal points
    R2 res_n; int N = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < N; node++)
	  Add(coefFi_curve(N+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointQuadrature(i, res_n);
      }
  }
  
  
  //! transformation Fi for quadrature points in the case of curved triangle
  void TriangleGeomReference::
  FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
			const Mesh<Dimension3>& mesh, int nquad) const
  {
    // isoparametric element
    // other points are computed by interpolation based on nodal points
    R3 res_n; int N = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < N; node++)
	  Add(coefFi_curve(N+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointQuadrature(i, res_n);
      }
  }
  
  
  //! transformation Fi for dof points in the case of curved triangle
  void TriangleGeomReference::
  FjElemDofCurve(const VectR2& s, SetPoints<Dimension2>& res,
		 const Mesh<Dimension2>& mesh, int nquad) const
  {
    // isoparametric element
    // other points are computed by interpolation based on nodal points
    R2 res_n; int N = this->points_nodal_nd.GetM();
    int offset = N + this->points2d.GetM();
    res.ReallocatePointsDof(this->points_dof2d.GetM());
    for (int i = 0; i < this->points_dof2d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < N; node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointDof(i,res_n);
      }
  }
  
  
  //! transformation DFi for nodal points in the case of curved triangle
  void TriangleGeomReference::
  DFjElemNodalCurve(const VectR2& s, const SetPoints<Dimension2>& PTReel,
		    SetMatrices<Dimension2>& res,
		    const Mesh<Dimension2>& mesh, int nquad) const
  {
    Matrix2_2 tmp;
    R2 grad; grad.Fill(0);
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
      {
	tmp.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  {
	    grad = this->GradientPhi_Nodal(node, i);
	    tmp(0,0) += PTReel.GetPointNodal(node)(0)*grad(0);
	    tmp(1,0) += PTReel.GetPointNodal(node)(1)*grad(0);
	    tmp(0,1) += PTReel.GetPointNodal(node)(0)*grad(1);
	    tmp(1,1) += PTReel.GetPointNodal(node)(1)*grad(1);
	  }
	res.SetPointNodal(i,tmp);
      } 
  }
  
  
  //! transformation DFi for quadrature points in the case of curved triangle
  void TriangleGeomReference::
  DFjElemQuadratureCurve(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			 SetMatrices<Dimension2>& res,
			 const Mesh<Dimension2>& mesh, int nquad) const
  {
    Matrix2_2 tmp;
    int N = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points2d.GetM());
    for (int i = 0; i < this->points2d.GetM(); i++)
      {
	tmp.Fill(0);
	for (int node = 0; node < N; node++)
	  Add(coefFi_curve(N+i, node), res.GetPointNodal(node), tmp);
	
	res.SetPointQuadrature(i, tmp);
      } 
  }
  
  
  //! transformation DFi for dof points in the case of curved triangle
  void TriangleGeomReference::
  DFjElemDofCurve(const VectR2& s, const SetPoints<Dimension2>& PTReel,
		  SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh, int nquad) const
  {
    Matrix2_2 tmp;
    int N = this->points_nodal_nd.GetM(); int offset = N + this->points2d.GetM();
    res.ReallocatePointsDof(this->points_dof2d.GetM());
    for (int i = 0; i < this->points_dof2d.GetM(); i++)
      {
	tmp.Fill(0);
	for (int node = 0; node < N; node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), tmp);
	
	res.SetPointDof(i,tmp);
      } 
  }
  
  
  //! displays informations about class TriangleGeomReference<type>  
  ostream& operator <<(ostream& out, const TriangleGeomReference& e)
  {
    out << static_cast<const ElementGeomReference<Dimension2>& >(e);
    return out;
  }
  
}


#define MONTJOIE_FILE_TRIANGLE_GEOM_REFERENCE_CXX
#endif
