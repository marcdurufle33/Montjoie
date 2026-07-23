#ifndef MONTJOIE_FILE_PYRAMID_HIERARCHIC_CXX

namespace Montjoie
{
  
  //! default constructor
  PyramidHierarchic::PyramidHierarchic() : PyramidReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    use_quadrature_for_rh = true;
  }
  

  //! how to number mesh
  void PyramidHierarchic::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return PyramidReference<1>::ConstructNumberMap(nmap, dg);

    nmap.SetFormulationDG(dg);    
    if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, (this->order+2)*(this->order+1)/2);
	nmap.SetNbDofQuadrangle(this->order, (this->order+1)*(this->order+1));
	nmap.SetNbDofPyramid(this->order, 0);
        
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
    nmap.SetNbDofPyramid(this->order, this->nb_dof_loc - 3*this->order*this->order - 2);
    
    nmap.SetEqualEdgesDofSymmetry(order, order-1);
    nmap.SetOddEdgesDofToSkewSymmetric(order);
    
    FindH1LinearCombinationRotation(nmap, *function_basis_tri, *function_basis_quad);    
  }
  

  size_t PyramidHierarchic::GetMemorySize() const
  {
    size_t taille = PyramidReference<1>::GetMemorySize();
    taille += jacobi_11_pol.GetMemorySize();
    taille += Seldon::GetMemorySize(jacobi_2ip1_pol) + Seldon::GetMemorySize(jacobi_2ip2_pol);
    taille += CoefLeg11.GetMemorySize();
    taille += CoefJacobi.GetMemorySize();
    taille += lob_quad.GetMemorySize() + lob_z.GetMemorySize();
    taille += NumDofs3D.GetMemorySize();    
    taille += NumDofs2D_tri.GetMemorySize() + NumDofs2D_quad.GetMemorySize() + CoordinateDofs.GetMemorySize();
    taille += ShTri.GetMemorySize() + ShLoc.GetMemorySize();
    taille += sh1_loc.GetMemorySize() + sh2_loc.GetMemorySize();
    taille += ch1_loc.GetMemorySize() + ch2_loc.GetMemorySize() + ch3_loc.GetMemorySize();
    taille += rh_loc.GetMemorySize() + const_rh.GetMemorySize();
    taille += MhLoc.GetMemorySize();
    taille += Seldon::GetMemorySize(NumProjOperator);
    taille += DerivDxtildeDx.GetMemorySize() + DerivDxtildeDz.GetMemorySize()
      + DerivDytildeDy.GetMemorySize() + DerivDytildeDz.GetMemorySize();
        
    return taille;
  }

  
  //! constructing finite element
  void PyramidHierarchic::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						 int rsurf_tri, int rsurf_quad,
						 int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    PyramidReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);

    ConstructFunctions();
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructElementaryMatrix(*this);

    ConstructMassMatrix();
    ConstructStiffnessMatrix();
  }

  
  //! construction of basis functions
  void PyramidHierarchic::ConstructFunctions()
  {
    nb_dof_loc = (order+1)*(order+2)*(2*order+3)/6;
    nb_dof_boundaries = 3*order*order + 2;
    nb_dof_tri = (order+1)*(order+2)/2;
    nb_dof_quad = (order+1)*(order+1);
    
    int r = order;
    // construction of nodes numbering
    function_basis_tri = new TriangleHierarchic();
    function_basis_tri->SetBasisType(TriangleHierarchic::INVARIANT_BASIS);
    function_basis_tri->ConstructFiniteElement(order);
    element_tri_surf = function_basis_tri;
    
    NumDofs2D_tri = function_basis_tri->GetNumDofs2D();
    
    function_basis_quad = new QuadrangleHierarchic();
    function_basis_quad->ConstructFiniteElement(order);
    element_quad_surf = function_basis_quad;

    NumDofs2D_quad = function_basis_quad->GetNumDofs2D();
    
    MeshNumbering<Dimension3>::ConstructPyramidalNumbering(r, NumDofs3D, CoordinateDofs);
    
    // changing numbering for faces
    int offset = 5 + 8*(order-1) - 4*order;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumDofs3D(i, j, 0) = offset + NumDofs2D_quad(i, j);
        
    offset += order + (order-1)*(order-1);
    int ntri = (order-1)*(order-2)/2;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	{
	  NumDofs3D(i, 0, j) = offset + NumDofs2D_tri(i, j);
	  NumDofs3D(r-j, i, j) = offset + ntri + NumDofs2D_tri(i, j);
	  NumDofs3D(i, r-j, j) = offset + 2*ntri + NumDofs2D_tri(i, j);
	  NumDofs3D(0, i, j) = offset + 3*ntri + NumDofs2D_tri(i, j);
	}
    
    offset += 3*order + 4*ntri;
    for (int diag = 1; diag < order; diag++)
      {
	for (int k = 1; k < order; k++)
	  for (int i = 1; i <= order-k; i++)
	    for (int j = 1; j <= order-k; j++)
	      if (max(i, j) + k == diag)
		NumDofs3D(i, j, k) = offset++;
      }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order-max(i, j); k++)
          {
            CoordinateDofs(NumDofs3D(i, j, k), 0) = i;
            CoordinateDofs(NumDofs3D(i, j, k), 1) = j;
            CoordinateDofs(NumDofs3D(i, j, k), 2) = k;
          }
    
    FacesDof.Reallocate((r+1)*(r+1), 5);
    FacesDof.Fill(-1);
    // dofs on faces
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	FacesDof(NumDofs2D_quad(i, j), 0) = NumDofs3D(i, j, 0);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        {
          int node = NumDofs2D_tri(i,j);
          FacesDof(node,1) = NumDofs3D(i, 0, j);
          FacesDof(node,2) = NumDofs3D(r-j, i, j);
          FacesDof(node,3) = NumDofs3D(i, r-j, j);
          FacesDof(node,4) = NumDofs3D(0, i, j);
        }

    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);
    
    // computing P_m^{1,1}
    GetJacobiPolynomial(jacobi_11_pol, order, Real_wp(1), Real_wp(1));
    jacobi_2ip1_pol.Reallocate(order);
    for (int i = 0; i < order; i++)
      GetJacobiPolynomial(jacobi_2ip1_pol(i), order, Real_wp(2*(i+1)+1), Real_wp(1));

    jacobi_2ip2_pol.Reallocate(order);
    for (int i = 0; i < order; i++)
      GetJacobiPolynomial(jacobi_2ip2_pol(i), order, Real_wp(2*(i+1)+2), Real_wp(1));
    
    CoefLeg11.Reallocate(order-1); CoefLeg11.Fill(0);
    CoefJacobi.Reallocate(order-1, order-1); CoefJacobi.Fill(0);
    VectReal_wp Pn;
    const VectReal_wp& points1d = this->Points1D();
    const VectReal_wp& weights1d = this->Weights1D();
    for (int i = 0; i <= order_quad; i++)
      {
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*points1d(i) - 1.0, Pn);
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
    
    lob_quad.AffectPoints(points1d);
    lob_z.AffectPoints(points1d_z);

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
  }
  
  
  //! construction of mass matrix
  void PyramidHierarchic::ConstructMassMatrix()
  {
    Array3D<Real_wp> ValPz(order+1, order+1, order_quad+1);
    Array3D<Real_wp> ValPz2(order+1, order+1, order_quad+1);
    ValPz.Fill(0);     ValPz2.Fill(0);
    ShLoc.Reallocate(order+1, order_quad+1);
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

        for (int j = 1; j < order; j++)
          {
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(j-1),
                                     order-2-j, 2.0*points1d_z(i) - 1.0, Pn);
            for (int k = 0; k < order-1-j; k++)
              ValPz(j-1, k, i) = Pn(k)*CoefJacobi(j-1, k);

            EvaluateJacobiPolynomial(jacobi_2ip2_pol(j-1),
                                     order-2-j, 2.0*points1d_z(i) - 1.0, Pn);
            for (int k = 0; k < order-1-j; k++)
              ValPz2(j-1, k, i) = Pn(k)*CoefJacobi(j-1, k);
          }
        
      }

    const VectR2& points2d_tri = this->Points2D_tri();
    ShTri.Reallocate(nb_dof_tri, points2d_tri.GetM());
    for (int i = 0; i < nb_dof_tri; i++)
      {
        int num_dof = FacesDof(i, 1);
        for (int j = 0; j < points2d_tri.GetM(); j++)
          {
	    int node = this->num_quad_points_surf(1)(j);
            ShTri(i, j) = Value_Phi(num_dof, node);
          }
      }
    
    QuadrangleHierarchic fe_quad;
    fe_quad.ConstructFiniteElement(order, this->GetGeometryOrder(), order_quad);
    
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2, Ch3;
    Ch1.Reallocate(nb_dof_loc, (order+1)*(order_quad+1)*(order+1));
    Ch2.Reallocate((order+1)*(order+1)*(order_quad+1), (order+1)*(order_quad+1)*(order_quad+1));
    Ch3.Reallocate((order+1)*(order_quad+1)*(order_quad+1), nb_points_quadrature_inside);

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();    
    Matrix<Real_wp, General, ArrayRowSparse> Sh1, Sh2;
    Sh1.Reallocate((order+1)*(order+1), (order_quad+1)*(order+1));
    Sh2.Reallocate((order_quad+1)*(order+1), (order_quad+1)*(order_quad+1));
    for (int i = 0; i <= order_quad; i++)
      {
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order; k++)
            Sh1.AddInteraction(NumDofs2D_quad(j, k), (order_quad+1)*j + i, ShLoc(k, i));
        
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order_quad; k++)
            Sh2.AddInteraction((order_quad+1)*j + k, NumQuad2D(i, k), ShLoc(j, i));
        
        /* First interpolation along z */
        /* we compute u(x, y, \xi_j), where xi_j are quadrature points for z-coordinate */
        
        // first vertices
        Ch1.AddInteraction(0, i, 1.0-points1d_z(i));
        Ch1.AddInteraction(1, (order_quad+1)*(order+1)*order+i, 1.0-points1d_z(i));
        Ch1.AddInteraction(2, (order_quad+1)*(order+2)*order+i, 1.0-points1d_z(i));
        Ch1.AddInteraction(3, (order_quad+1)*order + i, 1.0-points1d_z(i));
        Ch1.AddInteraction(4, i, points1d_z(i));
        Ch1.AddInteraction(4, (order_quad+1)*(order+1)*order+i, points1d_z(i));
        Ch1.AddInteraction(4, (order_quad+1)*(order+2)*order+i, points1d_z(i));
        Ch1.AddInteraction(4, (order_quad+1)*order + i, points1d_z(i));
        
        VectReal_wp pow_OneMinusZ(order+1);
        pow_OneMinusZ(0) = 1.0;
        for (int j = 0; j < order; j++)
          pow_OneMinusZ(j+1) = (1.0-points1d_z(i))*pow_OneMinusZ(j);
        
        // horizontal edges
        for (int j = 1; j < order; j++)
          {
            Ch1.AddInteraction(5+j-1, (order_quad+1)*(order+1)*j + i, pow_OneMinusZ(j+1));
            Ch1.AddInteraction(5+(order-1)+j-1,
                               (order_quad+1)*(order+1)*order + (order_quad+1)*j + i,
                               pow_OneMinusZ(j+1));
            Ch1.AddInteraction(5+2*(order-1)+j-1, (order_quad+1)*(order+1)*j
                               + (order_quad+1)*order + i, pow_OneMinusZ(j+1));
            Ch1.AddInteraction(5+3*(order-1)+j-1, (order_quad+1)*j + i, pow_OneMinusZ(j+1));
          }
        
        // vertical edges
        VectReal_wp valPi((order_quad+1)*(order_quad+1));
        VectReal_wp projPi((order+1)*(order+1));
        for (int num_loc = 4; num_loc < 8; num_loc++)
          for (int j = 1; j < order; j++)
            {
              int node = 5 + num_loc*(order-1) + j - 1;
              valPi.Fill(0);
              for (int k = 0; k <= order_quad; k++)
                for (int m = 0; m <= order_quad; m++)
                  {
                    int np = (order_quad+1)*((order_quad+1)*k + m) + i;
                    valPi(NumQuad2D(k, m)) = Value_Phi(node, np);
                  }
              
              fe_quad.ComputeProjectionDofRef(valPi, projPi);              
              for (int k = 0; k <= order; k++)
                for (int m = 0; m <= order; m++)
                  {
                    int np = (order_quad+1)*((order+1)*k + m) + i;
                    int ndof = NumDofs2D_quad(k, m);
                    Real_wp vloc = projPi(ndof);
                    if (k == 0)
                      if ((m > 0) && (m < order))
                        {
                          vloc = projPi(NumDofs2D_quad(k, order-m));
                          if (m%2 == 0)
                            vloc = -vloc;
                        }

                    if (m == order)
                      if ((k > 0) && (k < order))
                        {
                          vloc = projPi(NumDofs2D_quad(order-k, m));
                          if (k%2 == 0)
                            vloc = -vloc;
                        }
                    
                    if (abs(vloc) > 100.0*epsilon_machine)
                      Ch1.AddInteraction(node, np, vloc);
                  }              
            }
        
        // quadrilateral base
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order; k++)
            {
              int node = NumDofs3D(j, k, 0);
              Ch1.AddInteraction(node, (order_quad+1)*((order+1)*j + k) + i,
                                 pow_OneMinusZ(max(j, k)+1));
            }
        
        // triangular faces
        for (int num_loc = 1; num_loc < 5; num_loc++)
          for (int j = 1; j < order; j++)
            for (int k = 1; k < order-j; k++)
              {
                int np = i, node = -1;
                switch (num_loc)
                  {
                  case 1 :
		    node = NumDofs3D(j, 0, k);
                    np += (order_quad+1)*(order+1)*j;
                    break;
                  case 2 :
		    node = NumDofs3D(order-k, j, k);
                    np += (order_quad+1)*(order+1)*order + (order_quad+1)*j;
                    break;
                  case 3 :
		    node = NumDofs3D(j, order-k, k);
                    np += (order_quad+1)*(order+1)*j + (order_quad+1)*order;
                    break;
                  case  4 :
		    node = NumDofs3D(0, j, k);
                    np += (order_quad+1)*j;
                    break;
                  }
                
                Ch1.AddInteraction(node, np, ValPz(j-1, k-1, i)*points1d_z(i)*pow_OneMinusZ(j+1));
              }
            
        // interior nodes
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order; k++)
            for (int m = 1; m < order-max(j,k); m++)
              {
                int np = (order_quad+1)*((order+1)*j + k) + i;
		int node = NumDofs3D(j, k, m);
                Ch1.AddInteraction(node, np, ValPz2(max(j,k)-1, m-1, i)
                                   *points1d_z(i)*pow_OneMinusZ(max(j,k)+1));
              }

        /* Second interpolation along y */
        /* we compute u(x, \xi_i, \xi_j), where xi_i are quadrature points in y-coordinate */
        
        // contribution of 1-D basis functions on quadrature points
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order; k++)
            for (int m = 0; m <= order_quad; m++)
              Ch2.AddInteraction((order_quad+1)*((order+1)*j + k) + m,
                                 (order_quad+1)*((order_quad+1)*j+i)+m, ShLoc(k, i));
        
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
  void PyramidHierarchic::ConstructStiffnessMatrix()
  {
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    lob_z.ComputeGradPhi(1e3*epsilon_machine);

    const VectR3& points3d = this->PointsND();
    DerivDxtildeDx.Reallocate(points3d.GetM());
    DerivDxtildeDz.Reallocate(points3d.GetM());
    DerivDytildeDy.Reallocate(points3d.GetM());
    DerivDytildeDz.Reallocate(points3d.GetM());
        
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
                Rh.AddInteraction(node, 3*ny+1, lob_quad.GradPhi(j, m));
                
                int nz = (order_quad+1)*((order_quad+1)*i + j) + m;
                Rh.AddInteraction(node, 3*nz+2, lob_z.GradPhi(k, m));
              }

	    Real_wp x = points3d(node)(0), y = points3d(node)(1), z = points3d(node)(2);
	    DerivDxtildeDx(node) = 0.5/(1.0-z);
	    DerivDxtildeDz(node) = 0.5*x/square(1.0-z);
	    DerivDytildeDy(node) = 0.5/(1.0-z);
	    DerivDytildeDz(node) = 0.5*y/square(1.0-z);
          }
    
    Copy(Rh, rh_loc);
    ConvertToSparse(const_grad_matrix, const_rh, 1e3*epsilon_machine);
    ConvertToSparse(mass_matrix, MhLoc, 1e3*epsilon_machine);
  }
  
  
  //! Vh is overwritten by M Vh where M is the mass matrix
  template<class Vector1>
  void PyramidHierarchic::MltMassMatrixGen(Vector1& Vh) const
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
  void PyramidHierarchic::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    ElementReference<Dimension3, 1>::ApplyChGen(Uh, Vh);
    /*Vector2 Ux(ch3_loc.GetM()), Uy(ch2_loc.GetM());
    Mlt(ch3_loc, Uh, Ux);
    Mlt(ch2_loc, Ux, Uy);
    Mlt(ch1_loc, Uy, Vh);*/
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
  void PyramidHierarchic::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    ElementReference<Dimension3, 1>::ApplyChTransposeGen(Uh, Vh);
    /* Vector2 Ux(ch3_loc.GetM()), Uy(ch2_loc.GetM());
    Mlt(SeldonTrans, ch1_loc, Uh, Uy);
    Mlt(SeldonTrans, ch2_loc, Uy, Ux);
    Mlt(SeldonTrans, ch3_loc, Ux, Vh); */
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
  void PyramidHierarchic::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
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
  void PyramidHierarchic::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  void PyramidHierarchic::ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 dUh(Uh.GetM());
    // applying chaine rule to get derivatives on symmetric pyramid
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	dUh(3*i) = DerivDxtildeDx(i)*Uh(3*i) + DerivDxtildeDz(i)*Uh(3*i+2);
	dUh(3*i+1) = Uh(3*i+1)*DerivDytildeDy(i) + Uh(3*i+2)*DerivDytildeDz(i);
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
  void PyramidHierarchic::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);

    // applying chaine rule to get derivatives on symmetric pyramid
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	Vh(3*i+2) += Vh(3*i+1)*DerivDytildeDz(i) + Vh(3*i)*DerivDxtildeDz(i);
	Vh(3*i+1) *= DerivDytildeDy(i);
	Vh(3*i) *= DerivDxtildeDx(i);
      }
  }

  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void PyramidHierarchic::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(const_rh, Uh, Vh);
  }
  
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void PyramidHierarchic::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  void PyramidHierarchic
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vh.Fill(0);
    if ((r == 0) || (r == order_quad))
      {
        if (num_loc == 0)
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
        if (num_loc == 0)
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
  void PyramidHierarchic::
  ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0) || (r == order_quad))
      {
        if (num_loc == 0)
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
	if (num_loc == 0)
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
  void PyramidHierarchic::
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
      for (int j = 1; j < rf; j++)
        for (int k = 1; k < rf-max(i,j); k++)
          {
            if ((i+k < rc) && (j+k < rc))
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
  void PyramidHierarchic::ComputeValuesPhiRef(const R3& point_loc, VectReal_wp& phi) const
  {
    // dof on a vertex ?
    Real_wp gamma1(-1), gamma2(-1), L1, L2, L3, L4, N1, N2, N3, N4, c;
    Real_wp x = point_loc(0);
    Real_wp y = point_loc(1);
    Real_wp z = point_loc(2);
    
    c = 2.0*z - 1.0;
    L1 = 0.5*(1.0 - x - z);
    L2 = 0.5*(1.0 + y - z);
    L3 = 0.5*(1.0 + x - z);
    L4 = 0.5*(1.0 - y - z);
    
    if (point_loc(2) != 1.0)
      {
	gamma1 = x/(1.0 - z);
	gamma2 = y/(1.0 - z);
        N1 = L1*L4/(1.0 - z);
        N2 = L3*L4/(1.0 - z);
        N3 = L2*L3/(1.0 - z);
        N4 = L1*L2/(1.0 - z);
      }
    else
      {
        N1 = 0;
        N2 = 0;
        N3 = 0;
        N4 = 0;
      }
    
    phi.Reallocate(nb_dof_loc);
    // dofs on vertices
    phi(0) = N1;
    phi(1) = N2;
    phi(2) = N3;
    phi(3) = N4;
    phi(4) = z;
    
    if (order <= 1)
      return;
    
    // dofs on edges
    int node = 5;
    Real_wp vloc, vali;
    VectReal_wp Px, Py, pow_OneMinusZ(order), Px2, Py2;
    Vector<VectReal_wp> Pz(order-1);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma1, Px);    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma2, Py);

    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, x, Px2);    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, y, Py2);
    for (int k = 0; k < order-1; k++)
      {
        Px(k) *= CoefLeg11(k); Px2(k) *= CoefLeg11(k);
        Py(k) *= CoefLeg11(k); Py2(k) *= CoefLeg11(k);
      }
    
    pow_OneMinusZ(0) = 1.0;
    for (int i = 1; i < order; i++)
      {
        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2, c, Pz(i-1));
        for (int k = 0; k < order-1-i; k++)
          Pz(i-1)(k) *= CoefJacobi(i-1, k);
        
        pow_OneMinusZ(i) = pow_OneMinusZ(i-1)*(1.0-point_loc(2));
      }
    
    // horizontal edges
    vloc = N1*L3;
    for (int i = 1; i < order; i++)
      {
        phi(node) = vloc*Px2(i-1);
        node++;
      }
    
    vloc = L2*N2;
    for (int i = 1; i < order; i++)
      {
        phi(node) = vloc*Py2(i-1);
        node++;
      }
    
    vloc = L1*N3;
    for (int i = 1; i < order; i++)
      {
        phi(node) = vloc*Px2(i-1);
        node++;
      }
    
    vloc = N4*L4;
    for (int i = 1; i < order; i++)
      {
        phi(node) = vloc*Py2(i-1);
        node++;
      }
    
    // vertical edges
    vloc = N1*z;
    VectReal_wp P1;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, z + 0.5*(x+y), P1);
    for (int k = 0; k < order-1; k++)
      P1(k) *= CoefLeg11(k);

    for (int i = 1; i < order; i++)
      {
        vali = vloc*P1(i-1);
        phi(node) = vali; node++;
      }
    
    vloc = N2*z;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, z + 0.5*(-x+y), P1);
    for (int k = 0; k < order-1; k++)
      P1(k) *= CoefLeg11(k);
    
    for (int i = 1; i < order; i++)
      {
        vali = vloc*P1(i-1);
        phi(node) = vali; node++;
      }

    vloc = N3*z;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, z + 0.5*(-x-y), P1);
    for (int k = 0; k < order-1; k++)
      P1(k) *= CoefLeg11(k);
    
    for (int i = 1; i < order; i++)
      {
        vali = vloc*P1(i-1);
        phi(node) = vali; node++;
      }

    vloc = N4*z;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, z + 0.5*(x-y), P1);
    for (int k = 0; k < order-1; k++)
      P1(k) *= CoefLeg11(k);
    
    for (int i = 1; i < order; i++)
      {
        vali = vloc*P1(i-1);
        phi(node) = vali; node++;
      }

    // dofs on quadrilateral base
    vloc = N1*N3;
    int offset = 5 + 8*(order-1) - 4*order;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	{
	  node = offset + NumDofs2D_quad(i, j);
	  vali = pow_OneMinusZ(max(i, j) - 1)*Px(i-1)*Py(j-1);
	  phi(node) = vloc*vali;
	}
    
    // dofs on triangular faces
    vloc = N1*L3*z;
    offset += order + (order-1)*(order-1);
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
        {
	  node = offset + NumDofs2D_tri(i, j);
          vali = pow_OneMinusZ(i-1)*Px(i-1)*Pz(i-1)(j-1);
          phi(node) = vloc*vali;
        }

    vloc = N2*L2*z;
    offset += (order-1)*(order-2)/2;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
        {
	  node = offset + NumDofs2D_tri(i, j);
          vali = pow_OneMinusZ(i-1)*Py(i-1)*Pz(i-1)(j-1);
          phi(node) = vloc*vali;
        }

    vloc = N3*L1*z;
    offset += (order-1)*(order-2)/2;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
        {
	  node = offset + NumDofs2D_tri(i, j);
          vali = pow_OneMinusZ(i-1)*Px(i-1)*Pz(i-1)(j-1);
          phi(node) = vloc*vali;
        }

    vloc = N4*L4*z;
    offset += (order-1)*(order-2)/2;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
        {
	  node = offset + NumDofs2D_tri(i, j);
          vali = pow_OneMinusZ(i-1)*Py(i-1)*Pz(i-1)(j-1);
          phi(node) = vloc*vali;
        }
    
    // dofs inside
    vloc = N1*N3*z;
    for (int i = 1; i < order; i++)
      {
        EvaluateJacobiPolynomial(jacobi_2ip2_pol(i-1), order-2, c, Pz(i-1));
        for (int k = 0; k < order-1-i; k++)
          Pz(i-1)(k) *= CoefJacobi(i-1, k);
      }
    
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
        for (int k = 1; k < order-max(i,j); k++)
          {
	    node = NumDofs3D(i, j, k);
            vali = pow_OneMinusZ(max(i, j)-1)*Px(i-1)*Py(j-1)*Pz(max(i, j)-1)(k-1);
            phi(node) = vloc*vali;
          }
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on point_loc
  */
  void PyramidHierarchic::ComputeGradientPhiRef(const R3& point_loc, VectR3& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc);
    Real_wp gamma1(-1), gamma2(-1), L1, L2, L3, L4, N1, N2, N3, N4, c;
    R3 dg1, dg2, dN1, dN2, dN3, dN4;
    Real_wp x = point_loc(0);
    Real_wp y = point_loc(1);
    Real_wp z = point_loc(2);
    
    c = 2.0*z - 1.0;
    L1 = 1.0 - x - z;
    L2 = 1.0 + y - z;
    L3 = 1.0 + x - z;
    L4 = 1.0 - y - z;
    
    if (point_loc(2) != 1.0)
      {
	gamma1 = x/(1.0 - z);
	gamma2 = y/(1.0 - z);
        N1 = 0.25*L1*L4/(1.0 - z);
        N2 = 0.25*L3*L4/(1.0 - z);
        N3 = 0.25*L2*L3/(1.0 - z);
        N4 = 0.25*L1*L2/(1.0 - z);
        
        dg1.Init(1.0/(1.0-z), 0, x/square(1.0-z));
        dg2.Init(0, 1.0/(1.0-z), y/square(1.0-z));
        dN1.Init(-L4/(1.0-z), -L1/(1.0-z), -(L1+L4)/(1.0-z) + L1*L4/square(1.0-z));
        dN2.Init(L4/(1.0-z), -L3/(1.0-z), -(L3+L4)/(1.0-z) + L3*L4/square(1.0-z));
        dN3.Init(L2/(1.0-z), L3/(1.0-z), -(L2+L3)/(1.0-z) + L2*L3/square(1.0-z));
        dN4.Init(-L2/(1.0-z), L1/(1.0-z), -(L1+L2)/(1.0-z) + L1*L2/square(1.0-z));
        dN1 *= 0.25; dN2 *= 0.25; dN3 *= 0.25; dN4 *= 0.25;
      }
    else
      {
        N1 = 0;
        N2 = 0;
        N3 = 0;
        N4 = 0;
      }
    
    L1 *= 0.5; L2 *= 0.5; L3 *= 0.5; L4 *= 0.5;
    
    // dofs on vertices
    grad_phi(0) = dN1;
    grad_phi(1) = dN2;
    grad_phi(2) = dN3;
    grad_phi(3) = dN4;
    grad_phi(4).Init(0, 0, 1.0);
    
    if (order <= 1)
      return;
    
    // dofs on edges
    int node = 5;
    Real_wp vloc, vali; R3 grad_vloc, gradi;
    VectReal_wp Px, Py, dPx, dPy, pow_OneMinusZ(order), dpow_OneMinusZ(order);
    VectReal_wp Px2, Py2, dPx2, dPy2;
    Vector<VectReal_wp> Pz(order-1), dPz(order-1);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma1, Px, dPx);    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma2, Py, dPy);    

    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, x, Px2, dPx2);    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, y, Py2, dPy2);    
    for (int k = 0; k < order-1; k++)
      {
        Px(k) *= CoefLeg11(k); Px2(k) *= CoefLeg11(k);
        dPx(k) *= CoefLeg11(k); dPx2(k) *= CoefLeg11(k);
        Py(k) *= CoefLeg11(k); Py2(k) *= CoefLeg11(k);
        dPy(k) *= CoefLeg11(k); dPy2(k) *= CoefLeg11(k);
      }
    
    pow_OneMinusZ(0) = 1.0;
    dpow_OneMinusZ(0) = 0.0;
    for (int i = 1; i < order; i++)
      {
        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2, c, Pz(i-1), dPz(i-1));
        for (int k = 0; k < order-1-i; k++)
          {
            Pz(i-1)(k) *= CoefJacobi(i-1, k);
            dPz(i-1)(k) *= CoefJacobi(i-1, k);
          }

        pow_OneMinusZ(i) = pow_OneMinusZ(i-1)*(1.0-point_loc(2));
        dpow_OneMinusZ(i) = -Real_wp(i)*pow_OneMinusZ(i-1);
      }
    
    // horizontal edges
    vloc = N1*L3;
    grad_vloc = dN1; Mlt(L3, grad_vloc);
    grad_vloc(0) += 0.5*N1;
    grad_vloc(2) -= 0.5*N1;
    for (int i = 1; i < order; i++)
      {
        vali = Px2(i-1);
        grad_phi(node).Init(grad_vloc(0)*vali+vloc*dPx2(i-1),
                            grad_vloc(1)*vali, grad_vloc(2)*vali);
        node++;
      }
    
    vloc = L2*N2;
    grad_vloc = dN2; Mlt(L2, grad_vloc);
    grad_vloc(1) += 0.5*N2;
    grad_vloc(2) -= 0.5*N2;
    for (int i = 1; i < order; i++)
      {
        vali = Py2(i-1);
        grad_phi(node).Init(grad_vloc(0)*vali,
                            grad_vloc(1)*vali+vloc*dPy2(i-1), grad_vloc(2)*vali);
        node++;
      }
    
    vloc = L1*N3;
    grad_vloc = dN3; Mlt(L1, grad_vloc);
    grad_vloc(0) -= 0.5*N3;
    grad_vloc(2) -= 0.5*N3;
    for (int i = 1; i < order; i++)
      {
        vali = Px2(i-1);
        grad_phi(node).Init(grad_vloc(0)*vali+vloc*dPx2(i-1),
                            grad_vloc(1)*vali, grad_vloc(2)*vali);
        node++;
      }
    
    vloc = N4*L4;
    grad_vloc = dN4; Mlt(L4, grad_vloc);
    grad_vloc(1) -= 0.5*N4;
    grad_vloc(2) -= 0.5*N4;
    for (int i = 1; i < order; i++)
      {
        vali = Py2(i-1);
        grad_phi(node).Init(grad_vloc(0)*vali,
                            grad_vloc(1)*vali+vloc*dPy2(i-1), grad_vloc(2)*vali);
        node++;
      }
    
    // vertical edges
    vloc = N1*z;
    grad_vloc = dN1; Mlt(z, grad_vloc);
    grad_vloc(2) += N1;
    VectReal_wp P1, dP1;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, z + 0.5*(x+y), P1, dP1);
    for (int k = 0; k < order-1; k++)
      {
        P1(k) *= CoefLeg11(k);
        dP1(k) *= CoefLeg11(k);
      }

    for (int i = 1; i < order; i++)
      {
        vali = P1(i-1);
        gradi.Init(0.5*dP1(i-1), 0.5*dP1(i-1), dP1(i-1));
        grad_phi(node).Init(grad_vloc(0)*vali+vloc*gradi(0),
                            grad_vloc(1)*vali+vloc*gradi(1), grad_vloc(2)*vali+vloc*gradi(2));
        node++;
      }
    
    vloc = N2*z;
    grad_vloc = dN2; Mlt(z, grad_vloc);
    grad_vloc(2) += N2;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, z + 0.5*(-x+y), P1, dP1);
    for (int k = 0; k < order-1; k++)
      {
        P1(k) *= CoefLeg11(k);
        dP1(k) *= CoefLeg11(k);
      }
    
    for (int i = 1; i < order; i++)
      {
        vali = P1(i-1);
        gradi.Init(-0.5*dP1(i-1), 0.5*dP1(i-1), dP1(i-1));
        grad_phi(node).Init(grad_vloc(0)*vali+vloc*gradi(0),
                            grad_vloc(1)*vali+vloc*gradi(1), grad_vloc(2)*vali+vloc*gradi(2));
        node++;
      }

    vloc = N3*z;
    grad_vloc = dN3; Mlt(z, grad_vloc);
    grad_vloc(2) += N3;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, z + 0.5*(-x-y), P1, dP1);
    for (int k = 0; k < order-1; k++)
      {
        P1(k) *= CoefLeg11(k);
        dP1(k) *= CoefLeg11(k);
      }
    
    for (int i = 1; i < order; i++)
      {
        vali = P1(i-1);
        gradi.Init(-0.5*dP1(i-1), -0.5*dP1(i-1), dP1(i-1));
        grad_phi(node).Init(grad_vloc(0)*vali+vloc*gradi(0),
                            grad_vloc(1)*vali+vloc*gradi(1), grad_vloc(2)*vali+vloc*gradi(2));
        node++;
      }

    vloc = N4*z;
    grad_vloc = dN4; Mlt(z, grad_vloc);
    grad_vloc(2) += N4;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, z + 0.5*(x-y), P1, dP1);
    for (int k = 0; k < order-1; k++)
      {
        P1(k) *= CoefLeg11(k);
        dP1(k) *= CoefLeg11(k);
      }
        
    for (int i = 1; i < order; i++)
      {
        vali = P1(i-1);
        gradi.Init(0.5*dP1(i-1), -0.5*dP1(i-1), dP1(i-1));
        grad_phi(node).Init(grad_vloc(0)*vali+vloc*gradi(0),
                            grad_vloc(1)*vali+vloc*gradi(1), grad_vloc(2)*vali+vloc*gradi(2));
        node++;
      }
    
    // dofs on quadrilateral base
    vloc = N1*N3;
    grad_vloc = dN1*N3 + N1*dN3;
    Real_wp omz, domz;
    int offset = 5 + 8*(order-1) - 4*order;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	{
	  node = offset + NumDofs2D_quad(i, j);
          omz = pow_OneMinusZ(max(i, j) - 1);
          domz = dpow_OneMinusZ(max(i, j) - 1);
	  vali = omz*Px(i-1)*Py(j-1);
          gradi.Init(dg1(0)*dPx(i-1)*omz*Py(j-1), dg2(1)*dPy(j-1)*Px(i-1)*omz,
                     omz*(dg1(2)*dPx(i-1)*Py(j-1)
                          + dg2(2)*Px(i-1)*dPy(j-1)) + Px(i-1)*Py(j-1)*domz);
	  
          grad_phi(node).Init(grad_vloc(0)*vali+vloc*gradi(0),
                              grad_vloc(1)*vali+vloc*gradi(1), grad_vloc(2)*vali+vloc*gradi(2));
	}
    
    // dofs on triangular faces
    vloc = N1*L3*z;
    grad_vloc = dN1; Mlt(L3*z, grad_vloc);
    grad_vloc(0) += 0.5*N1*z;
    grad_vloc(2) += N1*L3 - 0.5*N1*z;
    offset += order + (order-1)*(order-1);
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
        {
	  node = offset + NumDofs2D_tri(i, j);
          omz = pow_OneMinusZ(i-1);
          domz = dpow_OneMinusZ(i-1);
          vali = omz*Px(i-1)*Pz(i-1)(j-1);
          gradi.Init(dg1(0)*dPx(i-1)*Pz(i-1)(j-1)*omz, 0,
                     dg1(2)*dPx(i-1)*Pz(i-1)(j-1)*omz + 2.0*omz*Px(i-1)*dPz(i-1)(j-1)
                     + domz*Px(i-1)*Pz(i-1)(j-1));          
          grad_phi(node).Init(grad_vloc(0)*vali+vloc*gradi(0), grad_vloc(1)*vali+vloc*gradi(1),
                              grad_vloc(2)*vali+vloc*gradi(2));
        }

    vloc = N2*L2*z;
    grad_vloc = dN2; Mlt(L2*z, grad_vloc);
    grad_vloc(1) += 0.5*N2*z;
    grad_vloc(2) += N2*L2 - 0.5*N2*z;
    offset += (order-1)*(order-2)/2;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
        {
	  node = offset + NumDofs2D_tri(i, j);
          omz = pow_OneMinusZ(i-1);
          domz = dpow_OneMinusZ(i-1);
          vali = omz*Py(i-1)*Pz(i-1)(j-1);
          gradi.Init(0, dg2(1)*dPy(i-1)*Pz(i-1)(j-1)*omz, dg2(2)*dPy(i-1)*Pz(i-1)(j-1)*omz
                     + 2.0*omz*Py(i-1)*dPz(i-1)(j-1) + domz*Py(i-1)*Pz(i-1)(j-1));          
          grad_phi(node).Init(grad_vloc(0)*vali+vloc*gradi(0), grad_vloc(1)*vali+vloc*gradi(1),
                              grad_vloc(2)*vali+vloc*gradi(2));
        }

    vloc = N3*L1*z;
    grad_vloc = dN3; Mlt(L1*z, grad_vloc);
    grad_vloc(0) -= 0.5*N3*z;
    grad_vloc(2) += N3*L1 - 0.5*N3*z;
    offset += (order-1)*(order-2)/2;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
        {
	  node = offset + NumDofs2D_tri(i, j);
          omz = pow_OneMinusZ(i-1);
          domz = dpow_OneMinusZ(i-1);
          vali = omz*Px(i-1)*Pz(i-1)(j-1);
          gradi.Init(dg1(0)*dPx(i-1)*Pz(i-1)(j-1)*omz, 0, dg1(2)*dPx(i-1)*Pz(i-1)(j-1)*omz
                     + 2.0*omz*Px(i-1)*dPz(i-1)(j-1) + domz*Px(i-1)*Pz(i-1)(j-1));          
          grad_phi(node).Init(grad_vloc(0)*vali+vloc*gradi(0), grad_vloc(1)*vali+vloc*gradi(1),
                              grad_vloc(2)*vali+vloc*gradi(2));
        }

    vloc = N4*L4*z;
    grad_vloc = dN4; Mlt(L4*z, grad_vloc);
    grad_vloc(1) -= 0.5*N4*z;
    grad_vloc(2) += N4*L4 - 0.5*N4*z;
    offset += (order-1)*(order-2)/2;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
        {
	  node = offset + NumDofs2D_tri(i, j);
          omz = pow_OneMinusZ(i-1);
          domz = dpow_OneMinusZ(i-1);
          vali = omz*Py(i-1)*Pz(i-1)(j-1);
          gradi.Init(0, dg2(1)*dPy(i-1)*Pz(i-1)(j-1)*omz, dg2(2)*dPy(i-1)*Pz(i-1)(j-1)*omz
                     + 2.0*omz*Py(i-1)*dPz(i-1)(j-1) + domz*Py(i-1)*Pz(i-1)(j-1));          
          grad_phi(node).Init(grad_vloc(0)*vali+vloc*gradi(0), grad_vloc(1)*vali+vloc*gradi(1),
                              grad_vloc(2)*vali+vloc*gradi(2));
        }
    
    // dofs inside
    vloc = N1*N3*z;
    grad_vloc = z*N1*dN3 + z*N3*dN1;
    grad_vloc(2) += N1*N3;
    for (int i = 1; i < order; i++)
      {
        EvaluateJacobiPolynomial(jacobi_2ip2_pol(i-1), order-2, c, Pz(i-1), dPz(i-1));
        for (int k = 0; k < order-1-i; k++)
          {
            Pz(i-1)(k) *= CoefJacobi(i-1, k);
            dPz(i-1)(k) *= CoefJacobi(i-1, k);
          }
      }
    
    Real_wp vz, dvz;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
        for (int k = 1; k < order-max(i,j); k++)
          {
            node = NumDofs3D(i, j, k);
	    omz = pow_OneMinusZ(max(i, j)-1);
            domz = dpow_OneMinusZ(max(i, j)-1);
            vz = Pz(max(i, j)-1)(k-1);
            dvz = 2.0*dPz(max(i, j)-1)(k-1);
            vali = omz*Px(i-1)*Py(j-1)*vz;
            gradi(0) = dg1(0)*omz*dPx(i-1)*Py(j-1)*vz;
            gradi(1) = dg2(1)*omz*Px(i-1)*dPy(j-1)*vz;
            gradi(2) = dg1(2)*omz*dPx(i-1)*Py(j-1)*vz + dg2(2)*omz*Px(i-1)*dPy(j-1)*vz
              + Px(i-1)*Py(j-1)*(domz*vz + omz*dvz);
            grad_phi(node).Init(grad_vloc(0)*vali+vloc*gradi(0), grad_vloc(1)*vali+vloc*gradi(1),
                                grad_vloc(2)*vali+vloc*gradi(2));
	  }
  }
  
    
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void PyramidHierarchic
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
  void PyramidHierarchic::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    VectReal_wp feval_weight(feval.GetM());
    for (int i = 0; i < feval.GetM(); i++)
      feval_weight(i) = feval(i)*this->WeightsND(i);
    
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }


  //! computation of projection between finite element of different orders  
  void PyramidHierarchic::
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

  
  void PyramidHierarchic::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {
    VectReal_wp contrib_orig(contrib);
    if (num_loc == 0)
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
    else
      for (int i = 1; i < this->order; i++)
	if (i%2 == 0)
	  contrib(2*this->order+i) = -contrib(2*this->order+i);
  }


  void PyramidHierarchic::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    VectComplex_wp contrib_orig(contrib);
    if (num_loc == 0)
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
    else
      for (int i = 1; i < this->order; i++)
	if (i%2 == 0)
	  contrib(2*this->order+i) = -contrib(2*this->order+i);
  }


  //! displays details about class PyramidHierarchic
  ostream& operator <<(ostream& out, const PyramidHierarchic& e)
  {
    out << static_cast<const PyramidReference<1>&>(e);
    out<<"Number dof over the pyramid "<<e.nb_dof_loc<<endl;
    return out;
  }
  
} // end namespace
  
#define MONTJOIE_FILE_PYRAMID_HIERARCHIC_CXX
#endif
