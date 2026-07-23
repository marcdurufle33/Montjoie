#ifndef MONTJOIE_FILE_WEDGE_HIERARCHIC_CXX

namespace Montjoie
{
  
  //! default constructor
  WedgeHierarchic::WedgeHierarchic() : WedgeReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    use_quadrature_for_rh = true;
  }
  

  //! how to number mesh
  void WedgeHierarchic::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return WedgeReference<1>::ConstructNumberMap(nmap, dg);

    nmap.SetFormulationDG(dg);
    if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, (this->order+2)*(this->order+1)/2);
	nmap.SetNbDofQuadrangle(this->order, (this->order+1)*(this->order+1));
	nmap.SetNbDofWedge(this->order, 0);
        
        // rotation on triangular faces
        ElementReference<Dimension2, 1>::
        FindH1RotationTri(this->order, this->Points2D_tri(), this->Weights2D_tri(),
                          function_basis_tri->GetValuePhi(), nmap);
        
        const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
        Matrix<int> FacesDof_Rotation_Quad;
        MeshNumbering<Dimension3>::
          GetRotationQuadrilateralFace(NumQuad2D, FacesDof_Rotation_Quad);
    
        nmap.SetFacesDofRotationQuad(this->order, FacesDof_Rotation_Quad);
        
        return;
      }

    nmap.SetNbDofVertex(this->order, 1);
    nmap.SetNbDofEdge(this->order, this->order-1);
    nmap.SetNbDofTriangle(this->order, (this->order-1)*(this->order-2)/2);
    nmap.SetNbDofQuadrangle(this->order, (this->order-1)*(this->order-1));
    nmap.SetNbDofWedge(order, nb_dof_loc - nb_dof_boundaries);
    
    nmap.SetEqualEdgesDofSymmetry(order, order-1);
    nmap.SetOddEdgesDofToSkewSymmetric(order);
    
    FindH1LinearCombinationRotation(nmap, *function_basis_tri, *function_basis_quad);
  }
  

  size_t WedgeHierarchic::GetMemorySize() const
  {
    size_t taille = WedgeReference<1>::GetMemorySize();
    taille += jacobi_11_pol.GetMemorySize();
    taille += Seldon::GetMemorySize(jacobi_2ip1_pol);
    taille += CoefLeg11.GetMemorySize();
    taille += CoefJacobi.GetMemorySize();
    taille += NumDofs3D.GetMemorySize();    
    taille += NumDofs2D_tri.GetMemorySize() + NumDofs2D_quad.GetMemorySize() + CoordinateDofs3D.GetMemorySize();
    taille += lob_quad.GetMemorySize() + lob_y.GetMemorySize();
    taille += points1d_y.GetMemorySize() + weights1d_y.GetMemorySize();
    taille += ShTri.GetMemorySize() + ShLoc.GetMemorySize();
    taille += sh1_loc.GetMemorySize() + sh2_loc.GetMemorySize();
    taille += ch1_loc.GetMemorySize() + ch2_loc.GetMemorySize() + ch3_loc.GetMemorySize();
    taille += rh_loc.GetMemorySize() + const_rh.GetMemorySize();
    taille += MhLoc.GetMemorySize();
    taille += Seldon::GetMemorySize(NumProjOperator);
    taille += DerivDxtildeDx.GetMemorySize() + DerivDxtildeDy.GetMemorySize();
    
    return taille;
  }

  
  //! constructing finite element  
  void WedgeHierarchic::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					       int rsurf_tri, int rsurf_quad,
					       int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    gauss_z = Globatto<Real_wp>::QUADRATURE_GAUSS;
    
    WedgeReference<1>::
      ConstructFiniteElement(r, rgeom, rquad, TriangleQuadrature::QUADRATURE_TENSOR, 
			     0, 0, -1, -1, gauss_z);

    ConstructFunctions();
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructElementaryMatrix(*this);

    ConstructMassMatrix();
    ConstructStiffnessMatrix();
  }
  
  
  //! construction of basis functions
  void WedgeHierarchic::ConstructFunctions()
  {
    nb_dof_loc = (order+1)*(order+2)*(order+1)/2;
    
    nb_dof_quad = (order+1)*(order+1);
    nb_dof_tri = (order+1)*(order+2)/2;
    nb_dof_boundaries = 4*order*order+2;
    
    int r = order;
    
    function_basis_tri = new TriangleHierarchic();
    function_basis_tri->SetBasisType(TriangleHierarchic::INVARIANT_BASIS);
    function_basis_tri->ConstructFiniteElement(order);
    element_tri_surf = function_basis_tri;
    
    NumDofs2D_tri = function_basis_tri->GetNumDofs2D();
    
    int nb = 3 + 2*(order-1);
    for (int i = 1; i < order; i++)
      NumDofs2D_tri(0, order-i) = nb++;
    
    function_basis_quad = new QuadrangleHierarchic();
    function_basis_quad->ConstructFiniteElement(order);
    element_quad_surf = function_basis_quad;

    NumDofs2D_quad = function_basis_quad->GetNumDofs2D();
    
    MeshNumbering<Dimension3>::ConstructPrismaticNumbering(r, NumDofs3D, CoordinateDofs3D);
    
    // changing dofs on faces
    int offset = 6 + 9*(r-1) - 3*r;
    for (int i = 1; i < r; i++)
      for (int j = 1; j < r-i; j++)
	NumDofs3D(i, j, 0) = offset + NumDofs2D_tri(i, j);
    
    offset += (r-1)*(r-2)/2 - r;
    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)
	{
	  NumDofs3D(i, 0, j) = offset + NumDofs2D_quad(i, j);
	  NumDofs3D(r-i, i, j) = offset + (r-1)*(r-1) + NumDofs2D_quad(i, j);
	  NumDofs3D(0, i, j) = offset + 2*(r-1)*(r-1) + NumDofs2D_quad(i, j);
	}
    
    offset += r + 3*(r-1)*(r-1);
    for (int i = 1; i < r; i++)
      for (int j = 1; j < r-i; j++)
	NumDofs3D(i, j, r) = offset + NumDofs2D_tri(i, j);
    
    offset += 3*r + (r-1)*(r-2)/2;
    for (int diag = 1; diag < r; diag++)
      for (int i = 1; i <= diag; i++)
	for (int j = 1; j <= diag-i; j++)
	  for (int k = 1; k <= diag; k++)
	    if (max(i+j, k) == diag)
	      NumDofs3D(i, j, k) = offset++;
    
    for (int k = 0; k <= order; k++)
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order-i; j++)
	  {
	    CoordinateDofs3D(NumDofs3D(i, j, k), 0) = i;
	    CoordinateDofs3D(NumDofs3D(i, j, k), 1) = j;
	    CoordinateDofs3D(NumDofs3D(i, j, k), 2) = k;
	  }
    
    FacesDof.Reallocate((r+1)*(r+1), 5);
    FacesDof.Fill(-1);
    // quadrangular faces
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  int node = NumDofs2D_quad(i,j);
	  FacesDof(node, 1) = NumDofs3D(i, 0, j);
	  FacesDof(node, 2) = NumDofs3D(r-i, i, j);
	  FacesDof(node, 3) = NumDofs3D(0, i, j);
	}
    
    // triangular faces
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  int node = NumDofs2D_tri(i,j);
	  FacesDof(node,0) = NumDofs3D(i, j, 0);
	  FacesDof(node,4) = NumDofs3D(i, j, r);
	}

    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);
    
    VectReal_wp points_dof1d; VectR3 points_dof3d;
    VectR2 points_dof2d_tri, points_dof2d_quad;
    
    points_dof1d = this->Points1D();
    points_dof2d_tri = this->Points2D_tri();
    points_dof2d_quad = this->Points2D_quad();
    points_dof3d = this->PointsND();
    
    this->SetPointsDof1D(points_dof1d);

    this->SetPointsDof2D_tri(points_dof2d_tri);
    this->SetPointsDof2D_quad(points_dof2d_quad);
    
    this->SetPointsDofND(points_dof3d);
    this->elt_geom.dof_equal_nodal = false;

    this->nb_points_dof_inside = nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;    

    // computing P_m^{1,1}
    GetJacobiPolynomial(jacobi_11_pol, order, Real_wp(1), Real_wp(1));
    jacobi_2ip1_pol.Reallocate(order);
    for (int i = 0; i < order; i++)
      GetJacobiPolynomial(jacobi_2ip1_pol(i), order, Real_wp(2*(i+1)+1), Real_wp(1.0));
    
    CoefLeg11.Reallocate(order-1); CoefLeg11.Fill(0);
    CoefJacobi.Reallocate(order-1, order-1); CoefJacobi.Fill(0);
    VectReal_wp Pn;
    const VectReal_wp& points1d = this->Points1D();
    const VectReal_wp& weights1d = this->Weights1D();
    for (int i = 0; i <= order_quad; i++)
      {
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*points1d(i) - 1.0, Pn);
        for (int j = 1; j < order; j++)
          CoefLeg11(j-1) += weights1d(i)*square(points1d(i)*(1.0-points1d(i))*Pn(j-1));

        for (int k = 0; k < order-1; k++)
          {
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(k), order-2, 2*points1d(i) - 1.0, Pn);
            for (int j = 1; j < order; j++)
              CoefJacobi(k, j-1)
                += weights1d(i)*square(points1d(i)*pow(1.0-points1d(i), k+1)*Pn(j-1));
          }
      }
    
    for (int j = 1; j < order; j++)
      {
        CoefLeg11(j-1) = 1.0/sqrt(CoefLeg11(j-1));
        for (int k = 0; k < order - 1; k++)
          if (abs(CoefJacobi(k, j-1)) > 0)
            CoefJacobi(k, j-1) = 1.0/sqrt(CoefJacobi(k, j-1));
      }

    ComputeGaussJacobi(points1d_y, weights1d_y, order, Real_wp(1), Real_wp(0));
    
    lob_quad.AffectPoints(points1d);
    lob_y.AffectPoints(points1d_y);
  }
  
  
  //! construction of mass matrix
  void WedgeHierarchic::ConstructMassMatrix()
  {
    Array3D<Real_wp> ValPy(order+1, order+1, order_quad+1);
    ValPy.Fill(0);
    ShLoc.Reallocate(order+1, order_quad+1);
    Matrix<Real_wp> ShY(order+1, order_quad+1);
    VectReal_wp Pn;
    const VectReal_wp& points1d = this->Points1D();
    for (int i = 0; i <= order_quad; i++)
      {
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*points1d(i) - 1.0, Pn);
        for (int j = 0; j < order-1; j++)
          Pn(j) *= CoefLeg11(j);
        
        ShLoc(0, i) = 1.0 - points1d(i);
        ShLoc(order, i) = points1d(i);
        for (int j = 1; j < order; j++)
          ShLoc(j, i) = (1.0-points1d(i))*points1d(i)*Pn(j-1);        

        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*points1d_y(i) - 1.0, Pn);
        for (int j = 0; j < order-1; j++)
          Pn(j) *= CoefLeg11(j);
        
        ShY(0, i) = 1.0 - points1d_y(i);
        ShY(order, i) = points1d_y(i);
        for (int j = 1; j < order; j++)
          ShY(j, i) = (1.0-points1d_y(i))*points1d_y(i)*Pn(j-1);        

        for (int j = 1; j < order; j++)
          {
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(j-1),
                                     order-2-j, 2.0*points1d_y(i) - 1.0, Pn);
            for (int k = 0; k < order-1-j; k++)
              ValPy(j-1, k, i) = Pn(k)*CoefJacobi(j-1, k);
          }
        
      }

    const VectR2& points2d_tri = this->Points2D_tri();
    ShTri.Reallocate(nb_dof_tri, points2d_tri.GetM());
    for (int i = 0; i < nb_dof_tri; i++)
      {
        int num_dof = FacesDof(i, 0);
        for (int j = 0; j < points2d_tri.GetM(); j++)
          {
            int node = this->num_quad_points_surf(0)(j);
            ShTri(i, j) = Value_Phi(num_dof, node);
          }
      }
    
    EdgeHierarchic fe_edge;
    fe_edge.ConstructFiniteElement(order);
    
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2, Ch3;
    Ch1.Reallocate(nb_dof_loc, (order+1)*(order_quad+1)*(order+1));
    Ch2.Reallocate((order+1)*(order+1)*(order_quad+1), (order+1)*(order_quad+1)*(order_quad+1));
    Ch3.Reallocate((order+1)*(order_quad+1)*(order_quad+1), nb_points_quadrature_inside);
    
    Matrix<Real_wp, General, ArrayRowSparse> Sh1, Sh2;
    Sh1.Reallocate((order+1)*(order+1), (order_quad+1)*(order+1));
    Sh2.Reallocate((order_quad+1)*(order+1), (order_quad+1)*(order_quad+1));
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    for (int i = 0; i <= order_quad; i++)
      {
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order; k++)
            Sh1.AddInteraction(NumDofs2D_quad(j, k), (order_quad+1)*j + i, ShLoc(k, i));
        
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order_quad; k++)
            Sh2.AddInteraction((order_quad+1)*j + k, NumQuad2D(i, k), ShLoc(j, i));
        
        /* First interpolation along y */
        /* we compute u(x, \xi_j, z), where xi_j are quadrature points for y-coordinate */
        
        // first vertices
        Ch1.AddInteraction(0, (order_quad+1)*i, 1.0-points1d_y(i));
        Ch1.AddInteraction(1, (order_quad+1)*(order+1)*order+(order_quad+1)*i, 1.0-points1d_y(i));
        Ch1.AddInteraction(2, (order_quad+1)*i, points1d_y(i));
        Ch1.AddInteraction(2, (order_quad+1)*(order+1)*order+(order_quad+1)*i, points1d_y(i));

        Ch1.AddInteraction(3, (order_quad+1)*i+order_quad, 1.0-points1d_y(i));
        Ch1.AddInteraction(4, (order_quad+1)*(order+1)*order+(order_quad+1)*i+order_quad,
                           1.0-points1d_y(i));
        Ch1.AddInteraction(5, (order_quad+1)*i+order_quad, points1d_y(i));
        Ch1.AddInteraction(5, (order_quad+1)*(order+1)*order+(order_quad+1)*i+order_quad,
                           points1d_y(i));
        
        VectReal_wp pow_OneMinusY(order+1);
        pow_OneMinusY(0) = 1.0;
        for (int j = 0; j < order; j++)
          pow_OneMinusY(j+1) = (1.0-points1d_y(i))*pow_OneMinusY(j);
        
        VectReal_wp valLi(order_quad+1);
        VectReal_wp projLi(order+1);
        // horizontal edges
        int node = 0;
        for (int j = 1; j < order; j++)
          {
            node = 5 + j;
            for (int k = 0; k <= order_quad; k++)
              valLi(k) = Value_Phi(node, (order_quad+1)*((order_quad+1)*k + i))
                / (1.0-points1d(0));
            
            fe_edge.ComputeProjectionDofRef(valLi, projLi);
            for (int k = 0; k <= order; k++)
              if (abs(projLi(k)) > 100.0*epsilon_machine)
                {
                  Ch1.AddInteraction(node, (order_quad+1)*((order+1)*k + i), projLi(k));
                  Ch1.AddInteraction(node+6*(order-1), (order_quad+1)*((order+1)*k + i) + order,
                                     projLi(k));
                }
            
            node = 5 + (order-1) + j;
            for (int k = 0; k <= order_quad; k++)
              valLi(k) = Value_Phi(node, (order_quad+1)*((order_quad+1)*k + i))
                / (1.0-points1d(0));
            
            fe_edge.ComputeProjectionDofRef(valLi, projLi);
            for (int k = 0; k <= order; k++)
              if (abs(projLi(k)) > 100.0*epsilon_machine)
                {
                  Ch1.AddInteraction(node, (order_quad+1)*((order+1)*k + i), projLi(k));
                  Ch1.AddInteraction(node+6*(order-1), (order_quad+1)*((order+1)*k + i)+order,
                                     projLi(k));
                }

            node = 5 + 2*(order-1) + j;
            for (int k = 0; k <= order_quad; k++)
              valLi(k) = Value_Phi(node, (order_quad+1)*((order_quad+1)*k + i))
                / (1.0-points1d(0));
            
            fe_edge.ComputeProjectionDofRef(valLi, projLi);
            for (int k = 0; k <= order; k++)
              if (abs(projLi(k)) > 100.0*epsilon_machine)
                {
                  Ch1.AddInteraction(node, (order_quad+1)*((order+1)*k + i), projLi(k));
                  Ch1.AddInteraction(node+6*(order-1), (order_quad+1)*((order+1)*k + i)+order,
                                     projLi(k));
                }
          }
        
        // vertical edges
        for (int j = 1; j < order; j++)
          {
            node = 5 + 3*(order-1) + j;
            Ch1.AddInteraction(node, (order_quad+1)*i + j, 1.0-points1d_y(i));
            
            node += order-1;
            Ch1.AddInteraction(node, (order_quad+1)*(order+1)*order+(order_quad+1)*i+j,
                               1.0-points1d_y(i));

            node += order-1;
            Ch1.AddInteraction(node, (order_quad+1)*i + j, points1d_y(i));
            Ch1.AddInteraction(node, (order_quad+1)*(order+1)*order+(order_quad+1)*i + j,
                               points1d_y(i));
          }
        
        // face z = 0
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order-j; k++)
            {
              node = NumDofs3D(j, k, 0);
              Ch1.AddInteraction(node, (order_quad+1)*((order+1)*j + i),
                                 points1d_y(i)*pow_OneMinusY(j+1)*ValPy(j-1, k-1, i));
            }

        // vertical faces
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order; k++)
            {
	      node = NumDofs3D(j, 0, k);
              Ch1.AddInteraction(node, (order_quad+1)*((order+1)*j + i) + k, pow_OneMinusY(j+1));
            }

        for (int j = 1; j < order; j++)
          for (int k = 1; k < order; k++)
            {
	      node = NumDofs3D(order-j, j, k);
              Ch1.AddInteraction(node, (order_quad+1)*((order+1)*order + i) + k, ShY(j, i));
            }

        for (int j = 1; j < order; j++)
          for (int k = 1; k < order; k++)
            {
	      node = NumDofs3D(0, j, k);
              Ch1.AddInteraction(node, (order_quad+1)*i + k, ShY(j, i));
            }
            
        // face z = 1
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order-j; k++)
            {
	      node = NumDofs3D(j, k, order);
              Ch1.AddInteraction(node, (order_quad+1)*((order+1)*j + i) + order,
                                 points1d_y(i)*pow_OneMinusY(j+1)*ValPy(j-1, k-1, i));
            }
        
        // interior
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order-j; k++)
            for (int m = 1; m < order; m++)
              {
		node = NumDofs3D(j, k, m);
                Ch1.AddInteraction(node, (order_quad+1)*((order+1)*j + i) + m,
                                   points1d_y(i)*pow_OneMinusY(j+1)*ValPy(j-1, k-1, i));
              }
        
        /* Second interpolation along z */
        /* we compute u(x, \xi_i, \xi_j), where xi_i are quadrature points in y-coordinate */
        
        // contribution of 1-D basis functions on quadrature points
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order; k++)
            for (int m = 0; m <= order_quad; m++)
              Ch2.AddInteraction((order_quad+1)*((order+1)*j + m) + k,
                                 (order_quad+1)*((order_quad+1)*j+m)+i, ShLoc(k, i));
        
        /* Third interpolation along x */
        /* we compute u(\xi_i, \xi_j, \xi_k), where xi_i are quadrature points in x-coordinate */
        
        // contribution of 1-D basis functions on quadrature points
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order_quad; k++)
            for (int m = 0; m <= order_quad; m++)
              Ch3.AddInteraction((order_quad+1)*((order_quad+1)*j + k) + m,
                                 (order_quad+1)*((order_quad+1)*i + k) + m, ShLoc(j, i));
        
      }
    
    Copy(Ch1, ch1_loc);
    Copy(Ch2, ch2_loc);
    Copy(Ch3, ch3_loc);

    Copy(Sh1, sh1_loc);
    Copy(Sh2, sh2_loc);
  }
  
  
  //! construction of stiffness matrix
  void WedgeHierarchic::ConstructStiffnessMatrix()
  {
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    lob_y.ComputeGradPhi(1e3*epsilon_machine);

    const VectR3& points3d = this->PointsND();
    DerivDxtildeDx.Reallocate(points3d.GetM());
    DerivDxtildeDy.Reallocate(points3d.GetM());
        
    Matrix<Real_wp, General, ArrayRowSparse> Rh;
    Rh.Reallocate(nb_points_quadrature_inside, 3*nb_points_quadrature_inside);
    for (int i = 0; i <= order_quad; i++)
      for (int j = 0; j <= order_quad; j++)
        for (int k = 0; k <= order_quad; k++)
          {
            int node = (order_quad+1)*((order_quad+1)*i + j) + k;
            for (int m = 0; m <= order_quad; m++)
              {
                int nx = (order_quad+1)*((order_quad+1)*m + j) + k;
                Rh.AddInteraction(node, 3*nx, lob_quad.GradPhi(i, m));
                
                int ny = (order_quad+1)*((order_quad+1)*i + m) + k;
                Rh.AddInteraction(node, 3*ny+1, lob_y.GradPhi(j, m));
                
                int nz = (order_quad+1)*((order_quad+1)*i + j) + m;
                Rh.AddInteraction(node, 3*nz+2, lob_quad.GradPhi(k, m));
              }

	    Real_wp x = points3d(node)(0), y = points3d(node)(1);
            DerivDxtildeDx(node) = 1.0/(1.0-y);
	    DerivDxtildeDy(node) = x/square(1.0-y);
          }
    
    Copy(Rh, rh_loc);
    ConvertToSparse(const_grad_matrix, const_rh, 1e4*epsilon_machine);
    ConvertToSparse(mass_matrix, MhLoc, 1e4*epsilon_machine);

  }
  
  
  //! Vh is overwritten by M Vh where M is the mass matrix
  template<class Vector1>
  void WedgeHierarchic::MltMassMatrixGen(Vector1& Vh) const
  {
    Vector1 Uh(Vh);
    Mlt(MhLoc, Uh, Vh);
  }
   
  
  //! Integration against basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Ch Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void WedgeHierarchic::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Ux(ch3_loc.GetM()), Uy(ch2_loc.GetM());
    Mlt(ch3_loc, Uh, Ux);
    Mlt(ch2_loc, Ux, Uy);
    Mlt(ch1_loc, Uy, Vh);
  }
  
  
  //! computation of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Ch* Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void WedgeHierarchic::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Ux(ch3_loc.GetM()), Uy(ch2_loc.GetM());
    Mlt(SeldonTrans, ch1_loc, Uh, Uy);
    Mlt(SeldonTrans, ch2_loc, Uy, Ux);
    Mlt(SeldonTrans, ch3_loc, Ux, Vh);
  }
  
  
  //! Integration against gradient of basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \nabla \varphi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Rh Uh
    where (Rh)_{i,j} = \nabla \varphi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void WedgeHierarchic::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(nb_points_quadrature_inside);
    Vector1 Uh2(Uh);
    ApplyRhQuadratureGen(Uh2, Uquad);
    ApplyChGen(Uquad, Vh);    
  }
  
  
  //! computation of gradient of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j grad phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = grad phi_i(xi_j)
   */
  template<class Vector1, class Vector2>
  void WedgeHierarchic::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(nb_points_quadrature_inside);
    ApplyChTransposeGen(Uh, Uquad);
    ApplyRhQuadratureTransposeGen(Uquad, Vh);
  }
  

  //! Integration against gradient of basis functions associated with quadrature points
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \nabla \psi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Rh Uh
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void WedgeHierarchic::ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const
  {
    // applying chaine rule to get derivatives on symmetric pyramid
    Vector1 dUh(Uh.GetM());
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	dUh(3*i) = DerivDxtildeDx(i)*Uh(3*i) + DerivDxtildeDy(i)*Uh(3*i+1);
	dUh(3*i+1) = Uh(3*i+1);
	dUh(3*i+2) = Uh(3*i+2);
      }
    
    Mlt(rh_loc, dUh, Vh);
  }
  
  
  //! Computation of gradient on quadrature points from values on quadrature points
  /*!
    \param[in] Uh values of u on quadrature points
    \param[out] Vh gradient of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Rh* Uh
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void WedgeHierarchic::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);

    // applying chaine rule to get derivatives on unit prism
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	Vh(3*i+1) += Vh(3*i)*DerivDxtildeDy(i);
	Vh(3*i) *= DerivDxtildeDx(i);
      }
  }

  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void WedgeHierarchic::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(const_rh, Uh, Vh);
  }
  
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void WedgeHierarchic::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, const_rh, Uh, Vh);
  }
  
  
  //! computation of u on quadrature points of a face
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the face
    \param[in] r order of quadrature rule of the face
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void WedgeHierarchic::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vh.Fill(0);
    if ((r == 0) || (r == order_quad))
      {
        if (num_loc%4 != 0)
          {
            Vector1 Uquad(this->nb_dof_quad), Ux(sh1_loc.GetM());
            for (int i = 0; i < this->nb_dof_quad; i++)
              Uquad(i) = Uh(FacesDof(i, num_loc));
            
            Mlt(SeldonTrans, sh1_loc, Uquad, Ux);
            Mlt(SeldonTrans, sh2_loc, Ux, Vh);
          }
        else
          {
            Vector1 Utri(this->nb_dof_tri);
            for (int i = 0; i < this->nb_dof_tri; i++)
              Utri(i) = Uh(FacesDof(i, num_loc));
            
            Mlt(SeldonTrans, ShTri, Utri, Vh);
          }
      }
    else
      {
        if (num_loc%4 != 0)
          {
	    int Nquad = this->Points2D_quad().GetM();
            Vector2 Vlob(Nquad);
            for (int i = 0; i < Nquad; i++)
              Vlob(i) = Uh(FacesDof(i, num_loc));
            
            Vh.Fill(0);
            typename Vector1::value_type vloc;
            const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
            const Matrix<int>& Num = NumProjOperator(r);
            for (int i = 0; i <= order; i++)
              for (int m = 0; m <= r; m++)
                {
                  vloc = 0;
                  for (int j = 0; j <= order; j++)
                    vloc += Ah(j, m)*Vlob(NumDofs2D_quad(i, j));
                  
                  for (int n = 0; n <= r; n++)
                    Vh(Num(n, m)) += Ah(i, n)*vloc;
                }
          }
        else
          {
            Vector1 Utri(this->nb_dof_tri);
            for (int i = 0; i < this->nb_dof_tri; i++)
              Utri(i) = Uh(FacesDof(i, num_loc));
            
            Mlt(SeldonTrans, ProjOperatorTriOrder(r), Utri, Vh);
          }
      }
  }
  
  
  //! integration against basis functions on a face
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void WedgeHierarchic::
  ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0) || (r == order_quad))
      {
        if (num_loc%4 != 0)
          {
            Vector1 Uquad(this->nb_dof_quad), Ux(sh1_loc.GetM());
            Mlt(sh2_loc, Uh, Ux);
            Mlt(sh1_loc, Ux, Uquad);
            
            for (int i = 0; i < this->nb_dof_quad; i++)
              Vh(FacesDof(i, num_loc)) += alpha*Uquad(i);
          }
        else
          {
            Vector1 Utri(this->nb_dof_tri);
            Mlt(ShTri, Uh, Utri);
            for (int i = 0; i < this->nb_dof_tri; i++)
              Vh(FacesDof(i, num_loc)) += alpha*Utri(i);
          }
      }
    else
      {
	if (num_loc%4 != 0)
          {
            Vector1 Vlob((order+1)*(order+1));
            Vlob.Fill(0);
            typename Vector1::value_type vloc;
            const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
            const Matrix<int>& Num = NumProjOperator(r);
            for (int i = 0; i <= order; i++)
              for (int m = 0; m <= r; m++)
                {
                  vloc = 0;
                  for (int n = 0; n <= r; n++)
                    vloc += Ah(i, n)*Uh(Num(n, m));
                  
                  for (int j = 0; j <= order; j++)
                    Vlob(NumDofs2D_quad(i, j)) += Ah(j, m)*vloc;
                }
            
	    int Nquad = this->Points2D_quad().GetM();
            int num_dof;
            for (int i = 0; i < Nquad; i++)
              {
                num_dof = FacesDof(i, num_loc);
                Vh(num_dof) += alpha*Vlob(i);
              }
          }
        else
          {
            Vector1 Utri(this->nb_dof_tri);
            Mlt(ProjOperatorTriOrder(r), Uh, Utri);
            for (int i = 0; i < this->nb_dof_tri; i++)
              Vh(FacesDof(i, num_loc)) += alpha*Utri(i);
          }
      }    
  }

  
  //! computation of prolongation operator inside an element
  /*!
    \param[in,out] proj prolongation operator
    \param[in,out] LocalProlongation prolongation operator
    \param[in] FaceCoarse coarse finite element
    \param[in] FaceFine fine finite element
   */
  void WedgeHierarchic::
  ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
                           const ElementReference_Dim<Dimension3>& FaceCoarse,
			   const ElementReference_Dim<Dimension3>& FaceFine) const
  {
    int rc = FaceCoarse.GetOrder();
    int rf = FaceFine.GetOrder();
    IVect num_coarse(FaceCoarse.GetNbDof()), num_fine(FaceCoarse.GetNbDof());
    num_coarse.Fill();
    num_fine.Fill(-1);
    
    // vertices
    for (int i = 0; i < FaceCoarse.GetNbVertices(); i++)
      num_fine(i) = i;
    
    int nc = FaceCoarse.GetNbVertices();
    int nf = nc;
    
    // edges
    for (int e = 0; e < FaceCoarse.GetNbEdges(); e++)
      for (int i = 0; i < rc-1; i++)
        num_fine(nc + (rc-1)*e + i) = nf + (rf-1)*e + i;
    
    nc += (rc-1)*FaceCoarse.GetNbEdges();
    nf += (rf-1)*FaceCoarse.GetNbEdges();
    
    // faces
    int type = FaceCoarse.GetHybridType();
    for (int f = 0; f < FaceCoarse.GetNbBoundaries(); f++)
      {
        if (MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(type, f))
          {
            for (int i = 0; i < rc-1; i++)
              for (int j = 0; j < rc-1; j++)
                num_fine(nc + (rc-1)*i + j) = nf + (rf-1)*i + j;
            
            nc += (rc-1)*(rc-1);
            nf += (rf-1)*(rf-1);
          }
        else
          {
            for (int i = 1; i < rf; i++)
              for (int j = 1; j < rf-i; j++)
                {
                  if ((i+j) < rc)
                    num_fine(nc++) = nf;
                  
                  nf++;
                }
          }
      }
    
    // interior
    for (int i = 1; i < rf; i++)
      for (int j = 1; j < rf-i; j++)
        for (int k = 1; k < rf; k++)
          {
            if ((i+j < rc) && (k < rc))
              num_fine(nc++) = nf;
            
            nf++;
          }
    
    proj.SetIdentity(num_coarse, num_fine, FaceFine.GetNbDof());
  }


  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void WedgeHierarchic::ComputeValuesPhiRef(const R3& point_loc, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);
    Real_wp L1 = 1.0-x-y, L2 = x, L3 = y;
    // dof on vertices
    phi(0) = (1.0-x-y)*(1.0-z);
    phi(1) = x*(1.0-z);
    phi(2) = y*(1.0-z);
    phi(3) = (1.0-x-y)*z;
    phi(4) = x*z;
    phi(5) = y*z;
    
    if (order <= 1)
      return;
    
    // dof on edges
    int node = 6;
    Real_wp a = 0, b = 2.0*y - 1.0, c = 2.0*z - 1.0;
    if (abs(1.0-y) > epsilon_machine)
      {
        a = 2.0*x/(1.0-y) - 1.0;
      }
    
    VectReal_wp P1, P2, P3, Px, Py2, Pz; Vector<VectReal_wp> Py(order-1);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L2-L1, P1);    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L3-L2, P2);    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L3-L1, P3);
    for (int k = 0; k < order-1; k++)
      {
        P1(k) *= CoefLeg11(k);
        P2(k) *= CoefLeg11(k);
        P3(k) *= CoefLeg11(k);
      }

    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, a, Px);    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, b, Py2);    
    for (int k = 0; k < order-1; k++)
      {
        Px(k) *= CoefLeg11(k);
        Py2(k) *= CoefLeg11(k);
      }
    
    VectReal_wp pow_OneMinusY(order);
    pow_OneMinusY(0) = 1.0;
    for (int i = 1; i < order; i++)
      {
        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, b, Py(i-1));    
        for (int k = 0; k < order-1-i; k++)
          Py(i-1)(k) *= CoefJacobi(i-1, k);
    
        pow_OneMinusY(i) = pow_OneMinusY(i-1)*(1.0-y);
      }
    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, c, Pz);
    for (int k = 0; k < order-1; k++)
      Pz(k) *= CoefLeg11(k);
            
    Real_wp vloc = L1*L2*(1.0-z);
    for (int i = 1; i < order; i++)
      phi(node++) = vloc*P1(i-1);
    
    vloc = L2*L3*(1.0-z);
    for (int i = 1; i < order; i++)
      phi(node++) = vloc*P2(i-1);
    
    vloc = L1*L3*(1.0-z);
    for (int i = 1; i < order; i++)
      phi(node++) = vloc*P3(i-1);
    
    vloc = L1*z*(1.0-z);
    for (int i = 1; i < order; i++)
      phi(node++) = vloc*Pz(i-1);
    
    vloc = L2*z*(1.0-z);
    for (int i = 1; i < order; i++)
      phi(node++) = vloc*Pz(i-1);

    vloc = L3*z*(1.0-z);
    for (int i = 1; i < order; i++)
      phi(node++) = vloc*Pz(i-1);

    vloc = L1*L2*z;
    for (int i = 1; i < order; i++)
      phi(node++) = vloc*P1(i-1);

    vloc = L2*L3*z;
    for (int i = 1; i < order; i++)
      phi(node++) = vloc*P2(i-1);

    vloc = L1*L3*z;
    for (int i = 1; i < order; i++)
      phi(node++) = vloc*P3(i-1);
    
    // dofs on faces
    vloc = L1*L2*L3*(1.0-z); Real_wp vali;
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          vali = CoefLeg11(i-1);
        else if (i == 2)
          vali = (L2-L1)*CoefLeg11(i-1);
        else
          vali = pow_OneMinusY(i-1)*Px(i-1);
        
        for (int j = 1; j < order-i; j++)
          {
	    node = NumDofs3D(i, j, 0);
	    phi(node) = vloc*vali*Py(i-1)(j-1);
	  }
      }
    
    vloc = L1*L2*z*(1.0-z);
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          vali = CoefLeg11(i-1);
        else if (i == 2)
          vali = (L2-L1)*CoefLeg11(i-1);
        else
          vali = pow_OneMinusY(i-1)*Px(i-1);
        
        for (int j = 1; j < order; j++)
          phi(NumDofs3D(i, 0, j)) = vloc*vali*Pz(j-1);
      }
    
    vloc = L2*L3*z*(1.0-z);
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
        phi(NumDofs3D(order-i, i, j)) = vloc*Py2(i-1)*Pz(j-1);


    vloc = L1*L3*z*(1.0-z);
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
        phi(NumDofs3D(0, i, j)) = vloc*Py2(i-1)*Pz(j-1);
    
    vloc = L1*L2*L3*z;
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          vali = CoefLeg11(i-1);
        else if (i == 2)
          vali = (L2-L1)*CoefLeg11(i-1);
        else
          vali = pow_OneMinusY(i-1)*Px(i-1);
        
        for (int j = 1; j < order-i; j++)
          phi(NumDofs3D(i, j, order)) = vloc*vali*Py(i-1)(j-1);
      }
    
    // dofs inside
    vloc = L1*L2*L3*z*(1.0-z);
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          vali = CoefLeg11(i-1);
        else if (i == 2)
          vali = (L2-L1)*CoefLeg11(i-1);
        else
          vali = pow_OneMinusY(i-1)*Px(i-1);
        
        for (int j = 1; j < order-i; j++)
          for (int k = 1; k < order; k++)
            phi(NumDofs3D(i, j, k)) = vloc*vali*Py(i-1)(j-1)*Pz(k-1);
      }
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res gradient of basis functions on point_loc
  */
  void WedgeHierarchic::ComputeGradientPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);
    Real_wp L1 = 1.0-x-y, L2 = x, L3 = y;
    // dofs on vertices
    res(0).Init(-(1.0-z), -(1.0-z), -(1.0-x-y));
    res(1).Init((1.0-z), 0, -x);
    res(2).Init(0, 1.0-z, -y);
    res(3).Init(-z, -z, 1.0-x-y);
    res(4).Init(z, 0, x);
    res(5).Init(0, z, y);
    
    if (order <= 1)
      return;
    
    // dofs on edges
    int node = 6;
    Real_wp a = 0, b = 2.0*y - 1.0, c = 2.0*z - 1.0;
    Real_wp da_dx = 0, da_dy = 0;
    if (abs(1.0-y) > epsilon_machine)
      {
        a = 2.0*x/(1.0-y) - 1.0;
        da_dx = 2.0/(1.0-y);
        da_dy = 2.0*x/square(1.0-y);
      }
    
    VectReal_wp P1, dP1, P2, dP2, P3, dP3, Px, dPx, Py2, dPy2, Pz, dPz;
    Vector<VectReal_wp> Py(order-1), dPy(order-1);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L2-L1, P1, dP1);    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L3-L2, P2, dP2);    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L3-L1, P3, dP3);
    for (int k = 0; k < order-1; k++)
      {
        P1(k) *= CoefLeg11(k); dP1(k) *= CoefLeg11(k);
        P2(k) *= CoefLeg11(k); dP2(k) *= CoefLeg11(k);
        P3(k) *= CoefLeg11(k); dP3(k) *= CoefLeg11(k);
      }

    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, a, Px, dPx);    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, b, Py2, dPy2);    
    for (int k = 0; k < order-1; k++)
      {
        Px(k) *= CoefLeg11(k); dPx(k) *= CoefLeg11(k);
        Py2(k) *= CoefLeg11(k); dPy2(k) *= CoefLeg11(k);
      }

    VectReal_wp pow_OneMinusY(order), dpow_OneMinusY(order);
    pow_OneMinusY(0) = 1.0;
    dpow_OneMinusY(0) = 0.0;
    for (int i = 1; i < order; i++)
      {
        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, b, Py(i-1), dPy(i-1));    
        for (int k = 0; k < order-1-i; k++)
          {
            Py(i-1)(k) *= CoefJacobi(i-1, k);
            dPy(i-1)(k) *= CoefJacobi(i-1, k);
          }
        
        pow_OneMinusY(i) = pow_OneMinusY(i-1)*(1.0-y);
        dpow_OneMinusY(i) = -Real_wp(i)*pow_OneMinusY(i-1);
      }
    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, c, Pz, dPz);
    for (int k = 0; k < order-1; k++)
      {
        Pz(k) *= CoefLeg11(k);
        dPz(k) *= CoefLeg11(k);
      }

    Real_wp vloc; R3 grad_vloc;
    vloc = L1*L2*(1.0-z);
    grad_vloc.Init((L1-L2)*(1.0-z), -L2*(1.0-z), -L1*L2);
    for (int i = 1; i < order; i++)
      res(node++).Init(2.0*vloc*dP1(i-1) + P1(i-1)*grad_vloc(0),
                       vloc*dP1(i-1) + P1(i-1)*grad_vloc(1), P1(i-1)*grad_vloc(2));

    vloc = L2*L3*(1.0-z);
    grad_vloc.Init(L3*(1.0-z), L2*(1.0-z), -L2*L3);
    for (int i = 1; i < order; i++)
      res(node++).Init(-vloc*dP2(i-1) + P2(i-1)*grad_vloc(0),
                       vloc*dP2(i-1) + P2(i-1)*grad_vloc(1), P2(i-1)*grad_vloc(2));
    
    vloc = L1*L3*(1.0-z);
    grad_vloc.Init(-L3*(1.0-z), (L1-L3)*(1.0-z), -L1*L3);
    for (int i = 1; i < order; i++)
      res(node++).Init(vloc*dP3(i-1) + P3(i-1)*grad_vloc(0),
                       2.0*vloc*dP3(i-1) + P3(i-1)*grad_vloc(1), P3(i-1)*grad_vloc(2));
    
    vloc = L1*z*(1.0-z);
    grad_vloc.Init(-z*(1.0-z), -z*(1.0-z), L1*(1.0-2.0*z));
    for (int i = 1; i < order; i++)
      res(node++).Init(Pz(i-1)*grad_vloc(0), Pz(i-1)*grad_vloc(1),
                       Pz(i-1)*grad_vloc(2) + 2.0*vloc*dPz(i-1));

    vloc = L2*z*(1.0-z);
    grad_vloc.Init(z*(1.0-z), 0, L2*(1.0-2.0*z));
    for (int i = 1; i < order; i++)
      res(node++).Init(Pz(i-1)*grad_vloc(0), Pz(i-1)*grad_vloc(1),
                       Pz(i-1)*grad_vloc(2) + 2.0*vloc*dPz(i-1));

    vloc = L3*z*(1.0-z);
    grad_vloc.Init(0, z*(1.0-z), L3*(1.0-2.0*z));
    for (int i = 1; i < order; i++)
      res(node++).Init(Pz(i-1)*grad_vloc(0), Pz(i-1)*grad_vloc(1),
                       Pz(i-1)*grad_vloc(2) + 2.0*vloc*dPz(i-1));

    vloc = L1*L2*z;
    grad_vloc.Init((L1-L2)*z, -L2*z, L1*L2);
    for (int i = 1; i < order; i++)
      res(node++).Init(2.0*vloc*dP1(i-1) + P1(i-1)*grad_vloc(0),
                       vloc*dP1(i-1) + P1(i-1)*grad_vloc(1), P1(i-1)*grad_vloc(2));

    vloc = L2*L3*z;
    grad_vloc.Init(L3*z, L2*z, L2*L3);
    for (int i = 1; i < order; i++)
      res(node++).Init(-vloc*dP2(i-1) + P2(i-1)*grad_vloc(0),
                       vloc*dP2(i-1) + P2(i-1)*grad_vloc(1), P2(i-1)*grad_vloc(2));
    
    vloc = L1*L3*z;
    grad_vloc.Init(-L3*z, (L1-L3)*z, L1*L3);
    for (int i = 1; i < order; i++)
      res(node++).Init(vloc*dP3(i-1) + P3(i-1)*grad_vloc(0),
                       2.0*vloc*dP3(i-1) + P3(i-1)*grad_vloc(1), P3(i-1)*grad_vloc(2));
    
    // dofs on faces
    vloc = L1*L2*L3*(1.0-z); Real_wp vali; R2 gradi;
    grad_vloc.Init(L3*(L1 - L2)*(1.0-z), L2*(L1 - L3)*(1.0-z), -L1*L2*L3);
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          {
            vali = CoefLeg11(i-1);
            gradi.Fill(0);
          }
        else if (i == 2)
          {
            vali = (L2 - L1)*CoefLeg11(i-1);
            gradi.Init(2.0, 1.0);
            gradi *= CoefLeg11(i-1);
          }
        else
          {
            vali = pow_OneMinusY(i-1)*Px(i-1);
            gradi.Init(da_dx*pow_OneMinusY(i-1)*dPx(i-1), da_dy*pow_OneMinusY(i-1)*dPx(i-1)
                       + dpow_OneMinusY(i-1)*Px(i-1));
          }
        
        for (int j = 1; j < order-i; j++)
	  {
	    node = NumDofs3D(i, j, 0);
	    res(node).Init(vloc*gradi(0)*Py(i-1)(j-1) + grad_vloc(0)*vali*Py(i-1)(j-1),
                           vloc*gradi(1)*Py(i-1)(j-1) + grad_vloc(1)*vali*Py(i-1)(j-1)
                           + 2.0*vloc*vali*dPy(i-1)(j-1),
                           grad_vloc(2)*vali*Py(i-1)(j-1));
	  }
      }
    
    vloc = L1*L2*z*(1.0-z);
    grad_vloc.Init((L1-L2)*z*(1.0-z), -L2*z*(1.0-z), L1*L2*(1.0-2.0*z));    
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          {
            vali = CoefLeg11(i-1);
            gradi.Fill(0);
          }
        else if (i == 2)
          {
            vali = (L2 - L1)*CoefLeg11(i-1);
            gradi.Init(2.0, 1.0);
            gradi *= CoefLeg11(i-1);
          }
        else
          {
            vali = pow_OneMinusY(i-1)*Px(i-1);
            gradi.Init(da_dx*pow_OneMinusY(i-1)*dPx(i-1), da_dy*pow_OneMinusY(i-1)*dPx(i-1)
                       + dpow_OneMinusY(i-1)*Px(i-1));
          }
        
        for (int j = 1; j < order; j++)
          res(NumDofs3D(i, 0, j))
	    .Init(Pz(j-1)*(grad_vloc(0)*vali + vloc*gradi(0)),
		  Pz(j-1)*(grad_vloc(1)*vali + vloc*gradi(1)), 
		  Pz(j-1)*grad_vloc(2)*vali + 2.0*vloc*vali*dPz(j-1));
      }
    
    vloc = L2*L3*z*(1.0-z);
    grad_vloc.Init(L3*z*(1.0-z), L2*z*(1.0-z), L2*L3*(1.0-2.0*z));
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
        {
          res(NumDofs3D(order-i, i, j))
	    .Init(grad_vloc(0)*Py2(i-1)*Pz(j-1),
		  grad_vloc(1)*Py2(i-1)*Pz(j-1) + 2.0*vloc*dPy2(i-1)*Pz(j-1),
		  grad_vloc(2)*Py2(i-1)*Pz(j-1) + 2.0*vloc*Py2(i-1)*dPz(j-1));
        }

    vloc = L1*L3*z*(1.0-z);
    grad_vloc.Init(-L3*z*(1.0-z), (L1-L3)*z*(1.0-z), L1*L3*(1.0-2.0*z));
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
        res(NumDofs3D(0, i, j))
	  .Init(grad_vloc(0)*Py2(i-1)*Pz(j-1),
		grad_vloc(1)*Py2(i-1)*Pz(j-1) + 2.0*vloc*dPy2(i-1)*Pz(j-1),
		grad_vloc(2)*Py2(i-1)*Pz(j-1) + 2.0*vloc*Py2(i-1)*dPz(j-1));
    
    vloc = L1*L2*L3*z;
    grad_vloc.Init(L3*(L1 - L2)*z, L2*(L1 - L3)*z, L1*L2*L3);
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          {
            vali = CoefLeg11(i-1);
            gradi.Fill(0);
          }
        else if (i == 2)
          {
            vali = (L2 - L1)*CoefLeg11(i-1);
            gradi.Init(2.0, 1.0);
            gradi *= CoefLeg11(i-1);
          }
        else
          {
            vali = pow_OneMinusY(i-1)*Px(i-1);
            gradi.Init(da_dx*pow_OneMinusY(i-1)*dPx(i-1), da_dy*pow_OneMinusY(i-1)*dPx(i-1)
                       + dpow_OneMinusY(i-1)*Px(i-1));
          }
        
        for (int j = 1; j < order-i; j++)
          res(NumDofs3D(i, j, order))
	    .Init(vloc*gradi(0)*Py(i-1)(j-1) + grad_vloc(0)*vali*Py(i-1)(j-1),
		  vloc*gradi(1)*Py(i-1)(j-1) + grad_vloc(1)*vali*Py(i-1)(j-1)
		  + 2.0*vloc*vali*dPy(i-1)(j-1),
		  grad_vloc(2)*vali*Py(i-1)(j-1));
      }
    
    // dofs inside
    vloc = L1*L2*L3*z*(1.0-z);
    grad_vloc.Init(L3*(L1-L2)*z*(1.0-z), L2*(L1-L3)*z*(1.0-z), L1*L2*L3*(1.0-2.0*z));
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          {
            vali = CoefLeg11(i-1);
            gradi.Fill(0);
          }
        else if (i == 2)
          {
            vali = (L2 - L1)*CoefLeg11(i-1);
            gradi.Init(2.0, 1.0);
            gradi *= CoefLeg11(i-1);
          }
        else
          {
            vali = pow_OneMinusY(i-1)*Px(i-1);
            gradi.Init(da_dx*pow_OneMinusY(i-1)*dPx(i-1), da_dy*pow_OneMinusY(i-1)*dPx(i-1)
                       + dpow_OneMinusY(i-1)*Px(i-1));
          }
        
        for (int j = 1; j < order-i; j++)
          for (int k = 1; k < order; k++)
            res(NumDofs3D(i, j, k))
	      .Init(Py(i-1)(j-1)*Pz(k-1)*(grad_vloc(0)*vali + vloc*gradi(0)),
		    Py(i-1)(j-1)*Pz(k-1)*(grad_vloc(1)*vali + vloc*gradi(1))
		    + 2.0*vloc*vali*dPy(i-1)(j-1)*Pz(k-1),
		    Py(i-1)(j-1)*Pz(k-1)*grad_vloc(2)*vali
		    + 2.0*vloc*vali*Py(i-1)(j-1)*dPz(k-1));
      }
  }
  
  

  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void WedgeHierarchic
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    VectComplex_wp feval_weight(feval.GetM());
    for (int i = 0; i < feval.GetM(); i++)
      feval_weight(i) = feval(i)*this->WeightsND(i);
    
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }
  
  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void WedgeHierarchic::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    VectReal_wp feval_weight(feval.GetM());
    for (int i = 0; i < feval.GetM(); i++)
      feval_weight(i) = feval(i)*this->WeightsND(i);
    
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }


  void WedgeHierarchic::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {
    VectReal_wp contrib_orig(contrib);
    if (num_loc%4 == 0)
      {
	for (int i = 1; i < this->order; i++)
	  {
	    if (i%2 == order%2)
	      contrib(2*this->order+i) = -contrib_orig(3*this->order-i);
	    else
	      contrib(2*this->order+i) = contrib_orig(3*this->order-i);
	  }
      }
    else
      {
	for (int i = 1; i < this->order; i++)
	  {
	    if (i%2 == this->order%2)
	      {
		contrib(2*this->order+i+1) = -contrib_orig(3*this->order-i+1);
		contrib(3*this->order+i) = -contrib_orig(4*this->order-i);
	      }
	    else
	      {
		contrib(2*this->order+i+1) = contrib_orig(3*this->order-i+1);
		contrib(3*this->order+i) = contrib_orig(4*this->order-i);
	      }	
	  }
      }
  }


  void WedgeHierarchic::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    VectComplex_wp contrib_orig(contrib);
    if (num_loc%4 == 0)
      {
	for (int i = 1; i < this->order; i++)
	  {
	    if (i%2 == order%2)
	      contrib(2*this->order+i) = -contrib_orig(3*this->order-i);
	    else
	      contrib(2*this->order+i) = contrib_orig(3*this->order-i);
	  }
      }
    else
      {
	for (int i = 1; i < this->order; i++)
	  {
	    if (i%2 == this->order%2)
	      {
		contrib(2*this->order+i+1) = -contrib_orig(3*this->order-i+1);
		contrib(3*this->order+i) = -contrib_orig(4*this->order-i);
	      }
	    else
	      {
		contrib(2*this->order+i+1) = contrib_orig(3*this->order-i+1);
		contrib(3*this->order+i) = contrib_orig(4*this->order-i);
	      }	
	  }
      }
  }

  
  //! computation of projection between finite element of different orders  
  void WedgeHierarchic::
  ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                            const ElementReference<Dimension2, 1>& Fb,
                                            const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi; VectR2 Points; IVect perm;
    ProjOperatorQuadOrder.Reallocate(rmax+1);
    NumProjOperator.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorQuadOrder(r).Reallocate(order+1, r+1);
	NumProjOperator(r).Reallocate(r+1, r+1);
	Points = Pts(i); perm.Reallocate(Points.GetM()); perm.Fill();
	Sort(Points, perm);
	int nb = 0;
	for (int j1 = 0; j1 <= r; j1++)
	  for (int j2 = 0; j2 <= r; j2++)
	    NumProjOperator(r)(j1, j2) = perm(nb++);
	
        VectReal_wp Pn;
	for (int i2 = 0; i2 <= r; i2++)
	  {
            ProjOperatorQuadOrder(r)(0, i2) = 1.0 - Points(i2)(1);
            ProjOperatorQuadOrder(r)(order, i2) = Points(i2)(1);
            
            EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*Points(i2)(1) - 1.0, Pn);
            for (int i1 = 1; i1 < order; i1++)
              ProjOperatorQuadOrder(r)(i1, i2) = Pn(i1-1)*CoefLeg11(i1-1);
          }
      }
  }


  //! displays details about class WedgeHierarchic
  ostream& operator <<(ostream& out, const WedgeHierarchic& e)
  {
    out<<static_cast<const WedgeReference<1>&>(e);
    out<<"Number dof over the Wedge "<<e.nb_dof_loc<<endl;
    return out;
  }
  
} // end namespace
  
#define MONTJOIE_FILE_WEDGE_HIERARCHIC_CXX
#endif
  
