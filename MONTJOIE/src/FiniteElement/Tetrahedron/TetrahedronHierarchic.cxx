#ifndef MONTJOIE_FILE_TETRAHEDRON_HIERARCHIC_CXX

namespace Montjoie
{
  
  //! default constructor
  TetrahedronHierarchic::TetrahedronHierarchic() : TetrahedronReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    use_quadrature_for_rh = true;
  }
  

  //! how to number mesh
  void TetrahedronHierarchic::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TetrahedronReference<1>::ConstructNumberMap(nmap, dg);
    
    nmap.SetFormulationDG(dg);
    if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, (this->order+2)*(this->order+1)/2);
	nmap.SetNbDofQuadrangle(this->order, (this->order+1)*(this->order+1));
	nmap.SetNbDofTetrahedron(this->order, 0);
        
        // rotation on triangular faces
        ElementReference<Dimension2, 1>::
          FindH1RotationTri(this->order, this->Points2D_tri(), this->Weights2D_tri(),
                            function_basis_tri->GetValuePhi(), nmap);
        
        return;
      }

    nmap.SetNbDofVertex(order, 1);
    nmap.SetNbDofEdge(order, order-1);
    nmap.SetNbDofTriangle(order, (order-2)*(order-1)/2);
    nmap.SetNbDofTetrahedron(order, (order-3)*(order-2)*(order-1)/6);
    
    nmap.SetEqualEdgesDofSymmetry(order, order-1);
    nmap.SetOddEdgesDofToSkewSymmetric(order);
    
    //FindH1SignEdge(*this, nmap);
    
    FindH1LinearCombinationRotation(nmap, *function_basis_tri, *function_basis_tri);    
  }


  size_t TetrahedronHierarchic::GetMemorySize() const
  {
    size_t taille = TetrahedronReference<1>::GetMemorySize();
    taille += jacobi_11_pol.GetMemorySize();
    taille += Seldon::GetMemorySize(jacobi_2ip1_pol);
    taille += CoefLeg11.GetMemorySize();
    taille += CoefJacobi.GetMemorySize();
    taille += ch1_loc.GetMemorySize() + ch2_loc.GetMemorySize() + ch3_loc.GetMemorySize();
    taille += rh_loc.GetMemorySize() + const_rh.GetMemorySize();
    taille += MhLoc.GetMemorySize();
    taille += points1d_y.GetMemorySize() + points1d_z.GetMemorySize();
    taille += gauss_x.GetMemorySize() + gauss_y.GetMemorySize() + gauss_z.GetMemorySize();
    taille += ShTri.GetMemorySize();
    taille += DerivDxtildeDx.GetMemorySize() + DerivDxtildeDy.GetMemorySize()
      + DerivDytildeDy.GetMemorySize() + DerivDytildeDz.GetMemorySize();
    taille += NumDofs2D_tri.GetMemorySize() + CoordinateDofs.GetMemorySize();
    taille += NumDofs3D.GetMemorySize();    
    return taille;
  }

  
  //! constructing finite element
  void TetrahedronHierarchic::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						     int rsurf_tri, int rsurf_quad,
						     int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    TetrahedronReference<1>::
      ConstructFiniteElement(r, rgeom, rquad, TetrahedronQuadrature::QUADRATURE_TENSOR);
    
    ConstructFunctions();
    this->Fb_geom.ComputeCoefficientTransformation();
    
    this->ConstructElementaryMatrix(*this);
    
    ConstructMassMatrix();
    ConstructStiffnessMatrix();
    
  }


  //! construction of basis functions
  void TetrahedronHierarchic::ConstructFunctions()
  {
    nb_dof_loc = (order+3)*(order+2)*(order+1)/6;

    int nb_dof_inside = (order-3)*(order-2)*(order-1)/6;
    nb_dof_boundaries = nb_dof_loc-nb_dof_inside;
        
    int r = order;
    Matrix<int> NumDofs2D_quad;
    // construction of nodes numbering
    function_basis_tri = new TriangleHierarchic();
    function_basis_tri->SetBasisType(TriangleHierarchic::INVARIANT_BASIS);
    function_basis_tri->ConstructFiniteElement(order);
    element_tri_surf = function_basis_tri;

    this->nb_dof_tri = function_basis_tri->GetNbDof();
    this->nb_dof_quad = 0;
    NumDofs2D_tri = function_basis_tri->GetNumDofs2D();
    
    MeshNumbering<Dimension3>::
      ConstructTetrahedralNumbering(order, NumDofs3D, CoordinateDofs);
    
    // changing numbering for faces
    int ntri = (order-1)*(order-2)/2;
    int offset = 4 + 6*(order-1) - 3*order;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
        {
          NumDofs3D(i, j, 0) = offset + NumDofs2D_tri(i, j);
          NumDofs3D(i, 0, j) = offset + ntri + NumDofs2D_tri(i, j);
          NumDofs3D(0, i, j) = offset + 2*ntri + NumDofs2D_tri(i, j);
          NumDofs3D(order-i-j, i, j) = offset + 3*ntri + NumDofs2D_tri(i, j);
        }
    
    // and internal nodes
    offset += 4*ntri + 3*order;
    for (int diag = 3; diag < order; diag++)
      for (int i = 1; i < order; i++)
        for (int j = 1; j < order-i; j++)
          {
            int k = diag - i - j;
            if (k >= 1)
              {
                NumDofs3D(i, j, k) = offset++;
              }
          }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        for (int k = 0; k <= order-i-j; k++)
          {
            CoordinateDofs(NumDofs3D(i, j, k), 0) = i;
            CoordinateDofs(NumDofs3D(i, j, k), 1) = j;
            CoordinateDofs(NumDofs3D(i, j, k), 2) = k;
          }
    
    FacesDof.Reallocate((r+1)*(r+2)/2, 4);
    // nodes on faces
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        {
          int node = NumDofs2D_tri(i,j);
          FacesDof(node,0) = NumDofs3D(i, j, 0);
          FacesDof(node,1) = NumDofs3D(i, 0, j);
          FacesDof(node,2) = NumDofs3D(0, i, j);
          FacesDof(node,3) = NumDofs3D(r-i-j, i, j);
        }
        
    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);

    // computing P_m^{1,1}
    GetJacobiPolynomial(jacobi_11_pol, order, Real_wp(1), Real_wp(1));
    jacobi_2ip1_pol.Reallocate(order);
    for (int i = 0; i < order; i++)
      GetJacobiPolynomial(jacobi_2ip1_pol(i), order, Real_wp(2*(i+1)+1), Real_wp(1));
    
    CoefLeg11.Reallocate(order-1); CoefLeg11.Fill(0);
    CoefJacobi.Reallocate(order-1, order-1); CoefJacobi.Fill(0);
    VectReal_wp Pn;
    const VectReal_wp& points1d = this->Points1D();
    const VectReal_wp& weights1d = this->Weights1D();
    const VectR3& points3d = this->PointsND();
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
    
    // retrieving points1d_z and points1d_y
    points1d_z.Reallocate(order+1);
    points1d_y.Reallocate(order+1);
    for (int i = 0; i <= order; i++)
      {
        points1d_z(i) = points3d(i)(2);
        points1d_y(i) = points3d(i*(order+1))(1)/(1.0-points1d_z(0));
      }

    gauss_x.AffectPoints(points1d);
    gauss_y.AffectPoints(points1d_y);
    gauss_z.AffectPoints(points1d_z);

    VectR2 points_dof2d_tri, points_dof2d_quad;
    VectR3 points_dof3d;
    VectReal_wp points_dof1d = this->Points1D();
    
    points_dof2d_tri = this->Points2D_tri();
    points_dof3d = this->PointsND();

    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDof2D_tri(points_dof2d_tri);
    
    this->elt_geom.dof_equal_nodal = false;
    
    this->SetPointsDofND(points_dof3d);
  }
  
  
  //! construction of mass matrix
  void TetrahedronHierarchic::ConstructMassMatrix()
  {
    Array3D<Real_wp> ValPy(order+1, order+1, order_quad+1);
    Array3D<Real_wp> ValPz(order+1, order+1, order_quad+1);
    ValPy.Fill(0);     ValPz.Fill(0);
    Matrix<Real_wp> ShX(order+1, order_quad+1);
    Matrix<Real_wp> ShY(order+1, order_quad+1);
    Matrix<Real_wp> ShZ(order+1, order_quad+1);
    VectReal_wp Pn;
    const VectReal_wp& points1d = this->Points1D();
    const VectR2& points2d_tri = this->Points2D_tri();
    for (int i = 0; i <= order_quad; i++)
      {
        // evaluating ShX
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*points1d(i) - 1.0, Pn);
        for (int j = 0; j < order-1; j++)
          Pn(j) *= CoefLeg11(j);
        
        ShX(0, i) = 1.0 - points1d(i);
        ShX(order, i) = points1d(i);
        for (int j = 1; j < order; j++)
          ShX(j, i) = (1.0-points1d(i))*points1d(i)*Pn(j-1);        
        
        // evaluating ShY
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*points1d_y(i) - 1.0, Pn);
        for (int j = 0; j < order-1; j++)
          Pn(j) *= CoefLeg11(j);
        
        ShY(0, i) = 1.0 - points1d_y(i);
        ShY(order, i) = points1d_y(i);
        for (int j = 1; j < order; j++)
          ShY(j, i) = (1.0-points1d_y(i))*points1d_y(i)*Pn(j-1);        

        // evaluating ShZ
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*points1d_z(i) - 1.0, Pn);
        for (int j = 0; j < order-1; j++)
          Pn(j) *= CoefLeg11(j);
        
        ShZ(0, i) = 1.0 - points1d_z(i);
        ShZ(order, i) = points1d_z(i);
        for (int j = 1; j < order; j++)
          ShZ(j, i) = (1.0-points1d_z(i))*points1d_z(i)*Pn(j-1);        

        for (int j = 1; j < order; j++)
          {
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(j-1), order-2-j, 2.0*points1d_y(i) - 1.0, Pn);
            for (int k = 0; k < order-1-j; k++)
              ValPy(j-1, k, i) = Pn(k)*CoefJacobi(j-1, k);

            EvaluateJacobiPolynomial(jacobi_2ip1_pol(j-1), order-2-j, 2.0*points1d_z(i) - 1.0, Pn);
            for (int k = 0; k < order-1-j; k++)
              ValPz(j-1, k, i) = Pn(k)*CoefJacobi(j-1, k);
          }
        
      }

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
    
    TriangleHierarchic fe_tri;
    fe_tri.SetBasisType(TriangleHierarchic::TENSOR_BASIS);
    fe_tri.ConstructFiniteElement(order, this->GetGeometryOrder(), order_quad,
                                  TriangleQuadrature::QUADRATURE_TENSOR);
    
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2, Ch3;
    Ch1.Reallocate(nb_dof_loc, (order+1)*(order_quad+1)*(order+2)/2);
    Ch2.Reallocate((order+1)*(order+2)*(order_quad+1)/2, (order+1)*(order_quad+1)*(order_quad+1));
    Ch3.Reallocate((order+1)*(order_quad+1)*(order_quad+1), nb_points_quadrature_inside);
    
    Array3D<int> NumPrism(order+1, order+1, order_quad+1);
    int node = 0;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        for (int k = 0; k <= order_quad; k++)
          NumPrism(i, j, k) = node++;
    
    for (int i = 0; i <= order_quad; i++)
      {
        
        /* First interpolation along z */
        /* we compute u(x, y, \xi_j), where xi_j are quadrature points for z-coordinate */
        
        // first vertices
        Ch1.AddInteraction(0, NumPrism(0, 0, i), 1.0-points1d_z(i));
        Ch1.AddInteraction(1, NumPrism(order, 0, i), 1.0-points1d_z(i));
        Ch1.AddInteraction(2, NumPrism(0, order, i), 1.0-points1d_z(i));
        Ch1.AddInteraction(3, NumPrism(0, 0, i), points1d_z(i));
        Ch1.AddInteraction(3, NumPrism(order, 0, i), points1d_z(i));
        Ch1.AddInteraction(3, NumPrism(0, order, i), points1d_z(i));
        
        VectReal_wp pow_OneMinusZ(order+1);
        pow_OneMinusZ(0) = 1.0;
        for (int j = 0; j < order; j++)
          pow_OneMinusZ(j+1) = (1.0-points1d_z(i))*pow_OneMinusZ(j);

        VectReal_wp pow_OneMinusY(order+1);
        pow_OneMinusY(0) = 1.0;
        for (int j = 0; j < order; j++)
          pow_OneMinusY(j+1) = (1.0-points1d_y(i))*pow_OneMinusY(j);
        
        // edges
        VectReal_wp valPi((order_quad+1)*(order_quad+1));
        VectReal_wp projPi((order+1)*(order+2)/2);
        for (int num_loc = 0; num_loc < 6; num_loc++)
          for (int j = 1; j < order; j++)
            {
              int node = 4 + num_loc*(order-1) + j - 1;
              valPi.Fill(0);
              for (int k = 0; k <= order_quad; k++)
                for (int m = 0; m <= order_quad; m++)
                  {
                    int np = (order_quad+1)*((order_quad+1)*k + m) + i;
                    valPi((order_quad+1)*k+m) = Value_Phi(node, np);
                  }
              
              fe_tri.ComputeProjectionDofRef(valPi, projPi);
              for (int k = 0; k <= order; k++)
                for (int m = 0; m <= order-k; m++)
                  {
                    int np = NumPrism(k, m, i);
                    int ndof = NumDofs2D_tri(k, m);
                    Real_wp vloc = projPi(ndof);
                    if (k == 0)
                      if ((m > 0) && (m < order))
                        {
                          vloc = projPi(NumDofs2D_tri(k, m));
                          if (m%2 == 0)
                            vloc = -vloc;
                        }
                    
                    if (abs(vloc) > 100.0*epsilon_machine)
                      Ch1.AddInteraction(node, np, vloc);
                  }
            }
        
        // faces
        int offset = 4 + 6*(order-1) - 3*order;
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order-j; k++)
            {
              node = offset + NumDofs2D_tri(j, k);
              Ch1.AddInteraction(node, NumPrism(j, k, i), pow_OneMinusZ(j+k+1));
            }
        
        offset += (order-1)*(order-2)/2;
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order-j; k++)
            {
              node = offset + NumDofs2D_tri(j, k);
              Ch1.AddInteraction(node, NumPrism(j, 0, i),
                                 points1d_z(i)*pow_OneMinusZ(j+1)*ValPz(j-1, k-1, i));
            }
        
        offset += (order-1)*(order-2)/2;
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order-j; k++)
            {
              node = offset + NumDofs2D_tri(j, k);
              Ch1.AddInteraction(node, NumPrism(0, j, i),
                                 points1d_z(i)*pow_OneMinusZ(j+1)*ValPz(j-1, k-1, i));
            }

        offset += (order-1)*(order-2)/2;
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order-j; k++)
            {
              node = offset + NumDofs2D_tri(j, k);
              Ch1.AddInteraction(node, NumPrism(order-j, j, i),
                                 points1d_z(i)*pow_OneMinusZ(j+1)*ValPz(j-1, k-1, i));
            }
        
        // interior of tetrahedron
        node = offset + (order-1)*(order-2)/2 + 3*order;
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order-j; k++)
            for (int m = 1; m < order-k-j; m++)
              {
                Ch1.AddInteraction(node, NumPrism(j, k, i),
                                   points1d_z(i)*pow_OneMinusZ(j+k+1)*ValPz(j+k-1, m-1, i));
                node++;
              }
        
        /* Second interpolation along y */
        /* we compute u(x, \xi_i, \xi_j), where xi_i are quadrature points in y-coordinate */
        
        // contribution of 1-D basis functions on quadrature points
        for (int m = 0; m <= order_quad; m++)
          {
            // vertices of the triangle
            Ch2.AddInteraction(NumPrism(0, 0, m), (order_quad+1)*i + m, 1.0-points1d_y(i));
            Ch2.AddInteraction(NumPrism(order, 0, m),
                               (order_quad+1)*((order+1)*order + i) + m, 1.0-points1d_y(i));
            Ch2.AddInteraction(NumPrism(0, order, m),
                               (order_quad+1)*((order+1)*order + i) + m, points1d_y(i));
            Ch2.AddInteraction(NumPrism(0, order, m), (order_quad+1)*i + m, points1d_y(i));
            
            // edges of the triangle
            for (int k = 1; k < order; k++)
              {
                Ch2.AddInteraction(NumPrism(k, 0, m),
                                   (order_quad+1)*((order_quad+1)*k + i) + m, pow_OneMinusY(k+1));
                Ch2.AddInteraction(NumPrism(order-k, k, m),
                                   (order_quad+1)*((order_quad+1)*order + i) + m, ShY(k, i));
                Ch2.AddInteraction(NumPrism(0, k, m), (order_quad+1)*i + m, ShY(k, i));
              }
            
            // interior of the triangle
            for (int j = 1; j < order; j++)
              for (int k = 1; k < order-j; k++)
                Ch2.AddInteraction(NumPrism(j, k, m), (order_quad+1)*((order_quad+1)*j + i) + m,
                                   points1d_y(i)*pow_OneMinusY(j+1)*ValPy(j-1, k-1, i));
            
            //for (int j = 0; j <= order; j++)
            //for (int k = 0; k <= order-j; k++)
            //  Ch2.AddInteraction(NumPrism(j, k, m),
            // (order_quad+1)*((order_quad+1)*j+i)+m, ShY(k, i));
          }
        
        /* Third interpolation along x */
        /* we compute u(\xi_i, \xi_j, \xi_k), where xi_i are quadrature points in x-coordinate */
        
        // contribution of 1-D basis functions on quadrature points
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order_quad; k++)
            for (int m = 0; m <= order_quad; m++)
              Ch3.AddInteraction((order_quad+1)*((order_quad+1)*j + k) + m,
                                 (order_quad+1)*((order_quad+1)*i + k) + m, ShX(j, i));
        
      }
    
    Copy(Ch1, ch1_loc);
    Copy(Ch2, ch2_loc);
    Copy(Ch3, ch3_loc);
  }
  
  
  //! construction of stiffness matrix
  void TetrahedronHierarchic::ConstructStiffnessMatrix()
  {    
    gauss_x.ComputeGradPhi(1e3*epsilon_machine);
    gauss_y.ComputeGradPhi(1e3*epsilon_machine);
    gauss_z.ComputeGradPhi(1e3*epsilon_machine);
    
    const VectR3& points3d = this->PointsND();
    DerivDxtildeDx.Reallocate(points3d.GetM());
    DerivDxtildeDy.Reallocate(points3d.GetM());
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
                Rh.AddInteraction(node, 3*nx, gauss_x.GradPhi(i, m));
                
                int ny = (order_quad+1)*((order_quad+1)*i + m) + k;
                Rh.AddInteraction(node, 3*ny+1, gauss_y.GradPhi(j, m));
                
                int nz = (order_quad+1)*((order_quad+1)*i + j) + m;
                Rh.AddInteraction(node, 3*nz+2, gauss_z.GradPhi(k, m));
              }

	    Real_wp x = points3d(node)(0), y = points3d(node)(1), z = points3d(node)(2);
            DerivDxtildeDx(node) = 1.0/(1.0-y-z);
	    DerivDxtildeDy(node) = x/square(1.0-y-z);
	    DerivDytildeDy(node) = 1.0/(1.0-z);
	    DerivDytildeDz(node) = y/square(1.0-z);
          }
    
    Copy(Rh, rh_loc);
    ConvertToSparse(const_grad_matrix, const_rh, 1e4*epsilon_machine);
    ConvertToSparse(mass_matrix, MhLoc, 1e3*epsilon_machine);

  }
  
    
  //! Vh is overwritten by M Vh where M is the mass matrix  
  template<class Vector1>
  void TetrahedronHierarchic::MltMassMatrixGen(Vector1& Vh) const
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
  void TetrahedronHierarchic::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
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
  void TetrahedronHierarchic::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  void TetrahedronHierarchic::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
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
  void TetrahedronHierarchic::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  void TetrahedronHierarchic::ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const
  {
    // applying chaine rule to get derivatives on symmetric pyramid
    Vector1 dUh(Uh.GetM());
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	dUh(3*i) = DerivDxtildeDx(i)*Uh(3*i) + DerivDxtildeDy(i)*Uh(3*i+1)
          + DerivDxtildeDy(i)*Uh(3*i+2);
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
  void TetrahedronHierarchic::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);
    
    // applying chaine rule to get derivatives on unit tetrahedron
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	Vh(3*i+2) += Vh(3*i+1)*DerivDytildeDz(i) + Vh(3*i)*DerivDxtildeDy(i);
	Vh(3*i+1) = Vh(3*i+1)*DerivDytildeDy(i) + Vh(3*i)*DerivDxtildeDy(i);
	Vh(3*i) *= DerivDxtildeDx(i);
      }
  }

  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void TetrahedronHierarchic::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(const_rh, Uh, Vh);
  }
  
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void TetrahedronHierarchic::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  void TetrahedronHierarchic
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vh.Fill(0);
    if ((r == 0) || (r == order_quad))
      {
        Vector1 Utri(this->nb_dof_tri);
        for (int i = 0; i < this->nb_dof_tri; i++)
          Utri(i) = Uh(FacesDof(i, num_loc));
        
        Mlt(SeldonTrans, ShTri, Utri, Vh);
      }
    else
      {
        Vector1 Utri(this->nb_dof_tri);
        for (int i = 0; i < this->nb_dof_tri; i++)
          Utri(i) = Uh(FacesDof(i, num_loc));
        
        Mlt(SeldonTrans, ProjOperatorTriOrder(r), Utri, Vh);
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
  void TetrahedronHierarchic::
  ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0) || (r == order_quad))
      {
        Vector1 Utri(this->nb_dof_tri);
        Mlt(ShTri, Uh, Utri);
        for (int i = 0; i < this->nb_dof_tri; i++)
          Vh(FacesDof(i, num_loc)) += alpha*Utri(i);
      }
    else
      {
        Vector1 Utri(this->nb_dof_tri);
        Mlt(ProjOperatorTriOrder(r), Uh, Utri);
        for (int i = 0; i < this->nb_dof_tri; i++)
          Vh(FacesDof(i, num_loc)) += alpha*Utri(i);
      }    
  }

  
  //! computation of prolongation operator inside an element
  /*!
    \param[in,out] proj prolongation operator
    \param[in,out] LocalProlongation prolongation operator
    \param[in] FaceCoarse coarse finite element
    \param[in] FaceFine fine finite element
   */
  void TetrahedronHierarchic::
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
        for (int k = 1; k < rf-i-j; k++)
          {
            if (i+j+k < rc)
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
  void TetrahedronHierarchic::ComputeValuesPhiRef(const R3& point_loc, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    
    // barycentric coordinates
    Real_wp L0 = 1.0 - point_loc(0) - point_loc(1) - point_loc(2);
    Real_wp L1 = point_loc(0);
    Real_wp L2 = point_loc(1);
    Real_wp L3 = point_loc(2);
    
    // 4 vertices
    phi(0) = L0;
    phi(1) = L1;
    phi(2) = L2;
    phi(3) = L3;
    
    // dofs on edges
    VectReal_wp P0, P1, P2, P3, P4, P5;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L1-L0, P0);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L2-L0, P1);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L3-L0, P2);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L2-L1, P3);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L3-L1, P4);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L3-L2, P5);
    
    for (int k = 0; k < order-1; k++)
      {
        P0(k) *= CoefLeg11(k);
        P1(k) *= CoefLeg11(k);
        P2(k) *= CoefLeg11(k);
        P3(k) *= CoefLeg11(k);
        P4(k) *= CoefLeg11(k);
        P5(k) *= CoefLeg11(k);
      }
    
    int num_dof = 4;
    for (int i = 1; i < order; i++)
      phi(num_dof++) = L0*L1*P0(i-1);

    for (int i = 1; i < order; i++)
      phi(num_dof++) = L0*L2*P1(i-1);

    for (int i = 1; i < order; i++)
      phi(num_dof++) = L0*L3*P2(i-1);
    
    for (int i = 1; i < order; i++)
      phi(num_dof++) = L1*L2*P3(i-1);
    
    for (int i = 1; i < order; i++)
      phi(num_dof++) = L1*L3*P4(i-1);
    
    for (int i = 1; i < order; i++)
      phi(num_dof++) = L2*L3*P5(i-1);

    // dofs on faces    
    Real_wp gamma1 = 0, gamma2 = 0;
    if (abs(L1+L0) > epsilon_machine)
      gamma1 = (L1-L0)/(L1+L0);
    
    if (abs(L0+L1+L2) > epsilon_machine)
      gamma2 = (L2-L1-L0)/(L2+L1+L0);
    
    Real_wp gamma3 = 2.0*point_loc(2) - 1.0;
    
    // face z = 0
    Real_wp vloc = L0*L1*L2;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma1, P0);
    for (int k = 0; k < order-1; k++)
      P0(k) *= CoefLeg11(k);
    
    Real_wp vali = 0, valj = 0, pow_L0pL1 = 1, pow_OneMinusZ = 1;
    int offset = 4 + 6*(order-1) - 3*order;
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          vali = CoefLeg11(i-1);
        else if (i == 2)
          vali = (L1 - L0)*CoefLeg11(i-1);
        else
          vali = P0(i-1)*pow_L0pL1;
        
        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, gamma2, P1);
        for (int k = 0; k < order-1-i; k++)
          P1(k) *= CoefJacobi(i-1, k);
        
        pow_OneMinusZ = 1.0;
        for (int j = 1; j < order-i; j++)
          {
            if (j == 1)
              valj = CoefJacobi(i-1, j-1);
            else if (j == 2)
              valj = (L2-L1-L0 - jacobi_2ip1_pol(i-1)(0,0)*(1.0-point_loc(2)))
                *CoefJacobi(i-1, j-1);
            else
              valj = P1(j-1)*pow_OneMinusZ;
            
            num_dof = offset + NumDofs2D_tri(i, j);
            phi(num_dof) = vloc*vali*valj;
            pow_OneMinusZ *= 1.0 - point_loc(2);
          }
        
        pow_L0pL1 *= L0+L1;
      }
    
    // face y = 0
    vloc = L0*L1*L3;
    pow_L0pL1 = 1.0;
    offset += (order-1)*(order-2)/2;
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          vali = CoefLeg11(i-1);
        else if (i == 2)
          vali = (L1 - L0)*CoefLeg11(i-1);
        else
          vali = P0(i-1)*pow_L0pL1;
        
        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, gamma3, P1);
        for (int k = 0; k < order-1-i; k++)
          P1(k) *= CoefJacobi(i-1, k);
        
        for (int j = 1; j < order-i; j++)
          {
            num_dof = offset + NumDofs2D_tri(i, j);
            phi(num_dof) = vloc*vali*P1(j-1);
          }
        
        pow_L0pL1 *= L0+L1;
      }
    
    // face x = 0
    offset += (order-1)*(order-2)/2;
    vloc = L0*L2*L3;
    pow_OneMinusZ = 1;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma2, P0);
    for (int k = 0; k < order-1; k++)
      P0(k) *= CoefLeg11(k);
 
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          vali = CoefLeg11(i-1);
        else if (i == 2)
          vali = (L2-L1-L0)*CoefLeg11(i-1);
        else
          vali = P0(i-1)*pow_OneMinusZ;
        
        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, gamma3, P1);
        for (int k = 0; k < order-1-i; k++)
          P1(k) *= CoefJacobi(i-1, k);
        
        for (int j = 1; j < order-i; j++)
          {
            num_dof = offset + NumDofs2D_tri(i, j);
            phi(num_dof) = vloc*vali*P1(j-1);
          }
        
        pow_OneMinusZ *= 1.0 - point_loc(2);
      }
    
    // face x + y + z = 1
    offset += (order-1)*(order-2)/2;
    vloc = L1*L2*L3;
    pow_OneMinusZ = 1;
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          vali = CoefLeg11(i-1);
        else if (i == 2)
          vali = (L2-L1-L0)*CoefLeg11(i-1);
        else
          vali = P0(i-1)*pow_OneMinusZ;
        
        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, gamma3, P1);
        for (int k = 0; k < order-1-i; k++)
          P1(k) *= CoefJacobi(i-1, k);
        
        for (int j = 1; j < order-i; j++)
          {
            num_dof = offset + NumDofs2D_tri(i, j);
            phi(num_dof) = vloc*vali*P1(j-1);
          }
        
        pow_OneMinusZ *= 1.0 - point_loc(2);
      }
    
    // dofs inside the tetrahedron
    offset += (order-1)*(order-2)/2 + 3*order;
    vloc = L0*L1*L2*L3;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma1, P0);
    for (int k = 0; k < order-1; k++)
      P0(k) *= CoefLeg11(k);
    
    pow_L0pL1 = 1.0;
    num_dof = offset;
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          vali = CoefLeg11(i-1);
        else if (i == 2)
          vali = (L1-L0)*CoefLeg11(i-1);
        else
          vali = pow_L0pL1*P0(i-1);
        
	EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, gamma2, P1);
        for (int k = 0; k < order-1-i; k++)
          P1(k) *= CoefJacobi(i-1, k);
 
        pow_OneMinusZ = 1;
	for (int j = 1; j < order-i; j++)
	  {
            if (j == 1)
              valj = CoefJacobi(i-1, j-1);
            else if (j == 2)
              valj = (L2-L1-L0 - jacobi_2ip1_pol(i-1)(0,0)*(1.0-point_loc(2)))
                *CoefJacobi(i-1, j-1);
            else
              valj = pow_OneMinusZ*P1(j-1);
            
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(i+j-1), order-2-i-j, gamma3, P2);
            for (int k = 0; k < order-1-i-j; k++)
              P2(k) *= CoefJacobi(i+j-1, k);
 
	    for (int k = 1; k < order-i-j; k++)
              {
                num_dof = NumDofs3D(i, j, k);
                phi(num_dof) = vloc*vali*valj*P2(k-1);
              }
            
            pow_OneMinusZ *= 1.0 - point_loc(2);
	  }
        
        pow_L0pL1 *= L0+L1;
      }
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on point_loc
  */
  void TetrahedronHierarchic::ComputeGradientPhiRef(const R3& point_loc, VectR3& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc);
    
    // barycentric coordinates
    Real_wp L0 = 1.0 - point_loc(0) - point_loc(1) - point_loc(2);
    Real_wp L1 = point_loc(0);
    Real_wp L2 = point_loc(1);
    Real_wp L3 = point_loc(2);
    
    // dof on vertices
    grad_phi(0).Init(-1.0, -1.0, -1.0);
    grad_phi(1).Init(1.0, 0, 0);
    grad_phi(2).Init(0, 1.0, 0);
    grad_phi(3).Init(0, 0, 1.0);
    
    // dofs on edges
    VectReal_wp P0, P1, P2, P3, P4, P5, dP0, dP1, dP2, dP3, dP4, dP5, dP6;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L1-L0, P0, dP0);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L2-L0, P1, dP1);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L3-L0, P2, dP2);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L2-L1, P3, dP3);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L3-L1, P4, dP4);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L3-L2, P5, dP5);
    for (int k = 0; k < order-1; k++)
      {
        P0(k) *= CoefLeg11(k); dP0(k) *= CoefLeg11(k);
        P1(k) *= CoefLeg11(k); dP1(k) *= CoefLeg11(k);
        P2(k) *= CoefLeg11(k); dP2(k) *= CoefLeg11(k);
        P3(k) *= CoefLeg11(k); dP3(k) *= CoefLeg11(k);
        P4(k) *= CoefLeg11(k); dP4(k) *= CoefLeg11(k);
        P5(k) *= CoefLeg11(k); dP5(k) *= CoefLeg11(k);
      }
    
    int num_dof = 4;
    Real_wp du, dv;
    for (int i = 1; i < order; i++)
      {
        du = L0*L1*dP0(i-1);
        dv = -L1*P0(i-1) + du;
        grad_phi(num_dof).Init((L0-L1)*P0(i-1) + 2.0*du, dv, dv);
        num_dof++;
      }

    for (int i = 1; i < order; i++)
      {
        du = L0*L2*dP1(i-1);
        dv = -L2*P1(i-1) + du;
        grad_phi(num_dof).Init(dv, (L0-L2)*P1(i-1) + 2.0*du, dv);
        num_dof++;
      }

    for (int i = 1; i < order; i++)
      {
        du = L0*L3*dP2(i-1);
        dv = -L3*P2(i-1) + du;
        grad_phi(num_dof).Init(dv, dv, (L0-L3)*P2(i-1) + 2.0*du);
        num_dof++;
      }
    
    for (int i = 1; i < order; i++)
      grad_phi(num_dof++).Init(L2*P3(i-1) - L1*L2*dP3(i-1), L1*P3(i-1) + L1*L2*dP3(i-1), 0);
    
    for (int i = 1; i < order; i++)
      grad_phi(num_dof++).Init(L3*P4(i-1) - L1*L3*dP4(i-1), 0, L1*P4(i-1) + L1*L3*dP4(i-1));
    
    for (int i = 1; i < order; i++)
      grad_phi(num_dof++).Init(0, L3*P5(i-1) - L2*L3*dP5(i-1), L2*P5(i-1) + L2*L3*dP5(i-1));

    // dofs on faces
    Real_wp gamma1 = 0, gamma2 = 0;
    R3 dgamma1, dgamma2;
    if (abs(L1+L0) > epsilon_machine)
      {
        // gamma1 = 2x/(1-y-z) - 1
        gamma1 = (L1-L0)/(L1+L0);
        dgamma1(0) = 2.0/(L1+L0);
        dgamma1(1) = 2.0*L1/square(L1+L0);
        dgamma1(2) = dgamma1(1);
      }
    
    if (abs(L0+L1+L2) > epsilon_machine)
      {
        // gamma2 = 2y/(1-z) - 1
        gamma2 = (L2-L1-L0)/(L2+L1+L0);
        dgamma2(0) = 0;
        dgamma2(1) = 2.0/(L2+L1+L0);
        dgamma2(2) = 2.0*L2/square(L0+L1+L2);
      }
    
    Real_wp gamma3 = 2.0*point_loc(2) - 1.0;
    Real_wp vali = 0, valj = 0, pow_L0pL1 = 1, pow_OneMinusZ = 1;
    Real_wp powM1_L0pL1 = 0, powM1_OneMinusZ = 0;
    R3 dvali, dvalj;
    
    // face z = 0
    Real_wp vloc = L0*L1*L2;
    Real_wp dv_dx = (L0-L1)*L2;
    Real_wp dv_dy = (L0-L2)*L1;
    Real_wp dv_dz = -L1*L2;
    Real_wp prod; R3 grad_prod;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma1, P0, dP0);
    for (int k = 0; k < order-1; k++)
      {
        P0(k) *= CoefLeg11(k);
        dP0(k) *= CoefLeg11(k);
      }

    int offset = 4 + 6*(order-1) - 3*order; 
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          {
            vali = CoefLeg11(i-1);
            dvali.Fill(0);
          }
        else if (i == 2)
          {
            vali = (L1 - L0)*CoefLeg11(i-1);
            dvali.Init(2.0, 1.0, 1.0);
            dvali *= CoefLeg11(i-1);
          }
        else
          {
            vali = pow_L0pL1*P0(i-1);
            dv = -Real_wp(i-1)*powM1_L0pL1*P0(i-1);
            du = pow_L0pL1*dP0(i-1);
            dvali.Init(dgamma1(0)*du, dgamma1(1)*du + dv, dgamma1(2)*du + dv);
          }
        
        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, gamma2, P1, dP1);
        for (int k = 0; k < order-1-i; k++)
          {
            P1(k) *= CoefJacobi(i-1, k);
            dP1(k) *= CoefJacobi(i-1, k);
          }
        
        pow_OneMinusZ = 1.0;
        for (int j = 1; j < order-i; j++)
          {
            if (j == 1)
              {
                valj = CoefJacobi(i-1, j-1);
                dvalj.Fill(0);
              }
            else if (j == 2)
              {
                valj = L2-L1-L0 - jacobi_2ip1_pol(i-1)(0,0)*(1.0-point_loc(2));
                valj *= CoefJacobi(i-1, j-1);
                dvalj.Init(0, 2.0, 1.0 + jacobi_2ip1_pol(i-1)(0,0));
                dvalj *= CoefJacobi(i-1, j-1);
              }
            else
              {
                valj = P1(j-1)*pow_OneMinusZ;
                dv = -Real_wp(j-1)*powM1_OneMinusZ*P1(j-1);
                du = pow_OneMinusZ*dP1(j-1);
                dvalj.Init(dgamma2(0)*du, dgamma2(1)*du, dgamma2(2)*du + dv);
              }
            
            num_dof = offset + NumDofs2D_tri(i, j);
            prod = vali*valj;
            grad_prod(0) = dvali(0)*valj + vali*dvalj(0);
            grad_prod(1) = dvali(1)*valj + vali*dvalj(1);
            grad_prod(2) = dvali(2)*valj + vali*dvalj(2);
            grad_phi(num_dof).Init(dv_dx*prod + vloc*grad_prod(0),
                                   dv_dy*prod + vloc*grad_prod(1), dv_dz*prod + vloc*grad_prod(2));
            powM1_OneMinusZ = pow_OneMinusZ;
            pow_OneMinusZ *= 1.0 - point_loc(2);
          }
        
        powM1_L0pL1 = pow_L0pL1;
        pow_L0pL1 *= L0+L1;
      }
    
    // face y = 0
    offset += (order-1)*(order-2)/2;
    vloc = L0*L1*L3;
    dv_dx = (L0-L1)*L3;
    dv_dy = -L1*L3;
    dv_dz = (L0-L3)*L1;
    pow_L0pL1 = 1.0;
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          {
            vali = CoefLeg11(i-1);
            dvali.Fill(0);
          }
        else if (i == 2)
          {
            vali = (L1 - L0)*CoefLeg11(i-1);
            dvali.Init(2.0, 1.0, 1.0);
            dvali *= CoefLeg11(i-1);
          }
        else
          {
            vali = pow_L0pL1*P0(i-1);
            dv = -Real_wp(i-1)*powM1_L0pL1*P0(i-1);
            du = pow_L0pL1*dP0(i-1);
            dvali.Init(dgamma1(0)*du, dgamma1(1)*du + dv, dgamma1(2)*du + dv);
          }
        
        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, gamma3, P1, dP1);
        for (int k = 0; k < order-1-i; k++)
          {
            P1(k) *= CoefJacobi(i-1, k);
            dP1(k) *= CoefJacobi(i-1, k);
          }
        
        for (int j = 1; j < order-i; j++)
          {
            num_dof = offset + NumDofs2D_tri(i, j);
            prod = vali*P1(j-1);
            grad_prod(0) = dvali(0)*P1(j-1);
            grad_prod(1) = dvali(1)*P1(j-1);
            grad_prod(2) = dvali(2)*P1(j-1) + 2.0*vali*dP1(j-1);
            grad_phi(num_dof).Init(dv_dx*prod + vloc*grad_prod(0),
                                   dv_dy*prod + vloc*grad_prod(1), dv_dz*prod + vloc*grad_prod(2));
          }
        
        powM1_L0pL1 = pow_L0pL1;
        pow_L0pL1 *= L0+L1;
      }
    
    // face x = 0
    offset += (order-1)*(order-2)/2;
    vloc = L0*L2*L3;
    dv_dx = -L2*L3;
    dv_dy = (L0-L2)*L3;
    dv_dz = (L0-L3)*L2;
    pow_OneMinusZ = 1.0;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma2, P0, dP0);
    for (int k = 0; k < order-1; k++)
      {
        P0(k) *= CoefLeg11(k);
        dP0(k) *= CoefLeg11(k);
      }
    
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          {
            vali = CoefLeg11(i-1);
            dvali.Fill(0);
          }
        else if (i == 2)
          {
            vali = (L2-L1-L0)*CoefLeg11(i-1);
            dvali.Init(0, 2.0, 1.0);
            dvali *= CoefLeg11(i-1);
          }
        else
          {
            vali = P0(i-1)*pow_OneMinusZ;
            dv = -Real_wp(i-1)*powM1_OneMinusZ*P0(i-1);
            du = pow_OneMinusZ*dP0(i-1);
            dvali.Init(dgamma2(0)*du, dgamma2(1)*du, dgamma2(2)*du + dv);
          }

        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, gamma3, P1, dP1);
        for (int k = 0; k < order-1-i; k++)
          {
            P1(k) *= CoefJacobi(i-1, k);
            dP1(k) *= CoefJacobi(i-1, k);
          }
        
        for (int j = 1; j < order-i; j++)
          {
            num_dof = offset + NumDofs2D_tri(i, j);
            prod = vali*P1(j-1);
            grad_prod(0) = dvali(0)*P1(j-1);
            grad_prod(1) = dvali(1)*P1(j-1);
            grad_prod(2) = dvali(2)*P1(j-1) + 2.0*vali*dP1(j-1);
            grad_phi(num_dof).Init(dv_dx*prod + vloc*grad_prod(0),
                                   dv_dy*prod + vloc*grad_prod(1), dv_dz*prod + vloc*grad_prod(2));
          }
        
        powM1_OneMinusZ = pow_OneMinusZ;
        pow_OneMinusZ *= 1.0 - point_loc(2);
      }
    
    // face x + y + z = 1
    offset += (order-1)*(order-2)/2;
    vloc = L1*L2*L3;
    dv_dx = L2*L3;
    dv_dy = L1*L3;
    dv_dz = L1*L2;
    pow_OneMinusZ = 1.0;
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          {
            vali = CoefLeg11(i-1);
            dvali.Fill(0);
          }
        else if (i == 2)
          {
            vali = (L2-L1-L0)*CoefLeg11(i-1);
            dvali.Init(0, 2.0, 1.0);
            dvali *= CoefLeg11(i-1);
          }
        else
          {
            vali = P0(i-1)*pow_OneMinusZ;
            dv = -Real_wp(i-1)*powM1_OneMinusZ*P0(i-1);
            du = pow_OneMinusZ*dP0(i-1);
            dvali.Init(dgamma2(0)*du, dgamma2(1)*du, dgamma2(2)*du + dv);
          }

        EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, gamma3, P1, dP1);
        for (int k = 0; k < order-1-i; k++)
          {
            P1(k) *= CoefJacobi(i-1, k);
            dP1(k) *= CoefJacobi(i-1, k);
          }
        
        for (int j = 1; j < order-i; j++)
          {
            num_dof = offset + NumDofs2D_tri(i, j);
            prod = vali*P1(j-1);
            grad_prod(0) = dvali(0)*P1(j-1);
            grad_prod(1) = dvali(1)*P1(j-1);
            grad_prod(2) = dvali(2)*P1(j-1) + 2.0*vali*dP1(j-1);
            grad_phi(num_dof).Init(dv_dx*prod + vloc*grad_prod(0),
                                   dv_dy*prod + vloc*grad_prod(1), dv_dz*prod + vloc*grad_prod(2));
          }

        powM1_OneMinusZ = pow_OneMinusZ;
        pow_OneMinusZ *= 1.0 - point_loc(2);
      }
        
    // dofs inside the tetrahedron
    offset += (order-1)*(order-2)/2 + 3*order;
    vloc = L0*L1*L2*L3;
    dv_dx = (L0-L1)*L2*L3;
    dv_dy = (L0-L2)*L1*L3;
    dv_dz = (L0-L3)*L1*L2;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma1, P0, dP0);
    for (int k = 0; k < order-1; k++)
      {
        P0(k) *= CoefLeg11(k);
        dP0(k) *= CoefLeg11(k);
      }
    
    pow_L0pL1 = 1.0;
    num_dof = offset;
    for (int i = 1; i < order; i++)
      {
        if (i == 1)
          {
            vali = CoefLeg11(i-1);
            dvali.Fill(0);
          }
        else if (i == 2)
          {
            vali = (L1-L0)*CoefLeg11(i-1);
            dvali.Init(2.0, 1.0, 1.0);
            dvali *= CoefLeg11(i-1);
          }
        else
          {
            vali = pow_L0pL1*P0(i-1);
            dv = -Real_wp(i-1)*powM1_L0pL1*P0(i-1);
            du = pow_L0pL1*dP0(i-1);
            dvali.Init(dgamma1(0)*du, dgamma1(1)*du + dv, dgamma1(2)*du + dv);
          }
        
	EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, gamma2, P1, dP1);
        for (int k = 0; k < order-1-i; k++)
          {
            P1(k) *= CoefJacobi(i-1, k);
            dP1(k) *= CoefJacobi(i-1, k);
          }
        
        pow_OneMinusZ = 1; powM1_OneMinusZ = 0;
	for (int j = 1; j < order-i; j++)
	  {
            if (j == 1)
              {
                valj = CoefJacobi(i-1, j-1);
                dvalj.Fill(0);
              }
            else if (j == 2)
              {
                valj = L2-L1-L0 - jacobi_2ip1_pol(i-1)(0,0)*(1.0-point_loc(2));
                dvalj.Init(0, 2.0, 1.0 + jacobi_2ip1_pol(i-1)(0,0));
                valj *= CoefJacobi(i-1, j-1);
                dvalj *= CoefJacobi(i-1, j-1);
              }
            else
              {
                valj = pow_OneMinusZ*P1(j-1);
                dv = -Real_wp(j-1)*powM1_OneMinusZ*P1(j-1);
                du = pow_OneMinusZ*dP1(j-1);
                dvalj.Init(dgamma2(0)*du, dgamma2(1)*du, dgamma2(2)*du + dv);
              }
            
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(i+j-1), order-2-i-j, gamma3, P2, dP2);
	    for (int k = 0; k < order-1-i-j; k++)
              {
                P2(k) *= CoefJacobi(i+j-1, k);
                dP2(k) *= CoefJacobi(i+j-1, k);
              }
            
            for (int k = 1; k < order-i-j; k++)
              {
                num_dof = NumDofs3D(i, j, k);
                prod = vali*valj*P2(k-1);
                grad_prod(0) = (dvali(0)*valj + vali*dvalj(0))*P2(k-1);
                grad_prod(1) = (dvali(1)*valj + vali*dvalj(1))*P2(k-1);
                grad_prod(2) = (dvali(2)*valj + vali*dvalj(2))*P2(k-1) + 2.0*vali*valj*dP2(k-1);
                grad_phi(num_dof).Init(dv_dx*prod + vloc*grad_prod(0),
                                       dv_dy*prod + vloc*grad_prod(1),
                                       dv_dz*prod + vloc*grad_prod(2));
              }
            
            powM1_OneMinusZ = pow_OneMinusZ;
            pow_OneMinusZ *= 1.0 - point_loc(2);
	  }
        
        powM1_L0pL1 = pow_L0pL1;
        pow_L0pL1 *= L0+L1;
      }
  }


  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void TetrahedronHierarchic
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
  void TetrahedronHierarchic
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    VectReal_wp feval_weight(feval.GetM());
    for (int i = 0; i < feval.GetM(); i++)
      feval_weight(i) = feval(i)*this->WeightsND(i);
    
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }


  void TetrahedronHierarchic::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {
    for (int i = 1; i < this->order; i++)
      if (i%2 == 0)
	contrib(2*this->order+i) = -contrib(2*this->order+i);
  }


  void TetrahedronHierarchic::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    for (int i = 1; i < this->order; i++)
      if (i%2 == 0)
	contrib(2*this->order+i) = -contrib(2*this->order+i);
  }
  

  //! displays details about class TetrahedronHierarchic
  ostream& operator <<(ostream& out, const TetrahedronHierarchic& e)
  {
    out<< static_cast<const TetrahedronReference<1>&>(e);
    out<<"Number dof over the tetrahedron "<<e.nb_dof_loc<<endl;
    return out;
  }
  
} // end namespace
  
#define MONTJOIE_FILE_TETRAHEDRON_HIERARCHIC_CXX
#endif
