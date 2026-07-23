#ifndef MONTJOIE_FILE_TETRAHEDRON_HCURL_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  TetrahedronHcurlFirstFamily::TetrahedronHcurlFirstFamily() : TetrahedronReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    type_basis = NODAL_GAUSS;
    //type_basis = NODAL_REGULAR;
  }
  
  
  //! how to number mesh
  void TetrahedronHcurlFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TetrahedronReference<2>::ConstructNumberMap(nmap, dg);
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, this->order*(this->order+2));
	nmap.SetNbDofTetrahedron(this->order, 0);

	// rotation of dofs on faces
	element_tri_surf->FindHcurlRotationTri(nmap, 0);
	
	return;
      }
    
    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    // nothing on quadrangular face !
    nmap.SetNbDofQuadrangle(order, 0);
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order-1));
    
    // dofs inside the tetrahedron
    nmap.SetNbDofTetrahedron(order, nb_dof_loc-nb_dof_boundaries);
    nmap.SetNbDofHexahedron(order, 0);

    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
    
    // sign changes on edges
    //FindHcurlSignEdge(nmap);
    
    // rotation of dofs on faces
    FindHcurlLinearCombinationRotation(nmap, *element_tri_surf, *element_tri_surf);
    
  }


  size_t TetrahedronHcurlFirstFamily::GetMemorySize() const
  {
    size_t taille = TetrahedronReference<2>::GetMemorySize();
    taille += InverseVDM_Edge.GetMemorySize();
    taille += InvWeightBasisRr.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    return taille;
  }


  //! construction of orthogonal functions
  void TetrahedronHcurlFirstFamily
  ::ConstructOrthogonalBasis(int r, int rgeom, int rquad, int type_quad)
  {
    TetrahedronReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    int nb_dof_inside_tri = order*(order-1);
    int nb_dof_inside_elt = order*(order-1)*(order-2)/2;    
    nb_dof_boundaries = 4*nb_dof_inside_tri + 6*order;
    nb_dof_loc = nb_dof_inside_elt + nb_dof_boundaries;
    
    const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefOddJacobi();
    const Matrix<Real_wp>& CoefEvenJacobi = Fb_geom.GetCoefEvenJacobi();
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefLegendre();
    
    // weight for orthonormalization
    InvWeightBasisRr.Reallocate(nb_dof_loc);
    int num = 0;
    for (int i = 0; i < order; i++)
      {
        for (int j = 0; j <= order-1-i; j++)
          {
            for (int k = 0; k < order-i-j; k++)
              {
                Real_wp coef = CoefLegendre(i)*CoefOddJacobi(i, j)*CoefEvenJacobi(i+j, k);
                InvWeightBasisRr(num) = coef;
                InvWeightBasisRr(num+1) = coef;
                InvWeightBasisRr(num+2) = coef;
                num += 3;
              }
            
            InvWeightBasisRr(num) = CoefLegendre(i)*CoefOddJacobi(i, j)
              *CoefEvenJacobi(i+j, order-1-i-j);
            
            InvWeightBasisRr(num+1) = InvWeightBasisRr(num);
            num += 2;
          }
        
        InvWeightBasisRr(num) = CoefLegendre(i)*CoefOddJacobi(i, order-1-i);
        num++;
      }
    
    InvWeightBasisRr.Fill(1.0);
    VectReal_wp Scale(nb_dof_loc);
    VectR3 phi(nb_dof_loc); Scale.Fill(0);
    for (int k = 0; k < nb_points_quadrature_inside; k++)
      {
        ComputeValuesPhiOrthoRef(this->PointsND(k), phi);
        for (int i = 0; i < nb_dof_loc; i++)
          Scale(i) += DotProd(phi(i), phi(i))*this->WeightsND(k);
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      InvWeightBasisRr(i) = 1.0/sqrt(Scale(i));

    // DISP(InvWeightBasisRr);

  }
  
  
  //! construction of finite element
  void TetrahedronHcurlFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf_tri, int rsurf_quad,
			   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    ConstructOrthogonalBasis(r, rgeom, rquad, type_quad);
    
    // computation of basis functions
    ConstructFunctions();
    this->nb_points_dof_inside = this->nb_dof_loc;
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    // construction of Value_PhiVec, Curl_Phi and elementary matrices
    ConstructHcurlElementaryMatrix();
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace(true);    
    
    TriangleHcurlFirstFamily* Fb_tri = new TriangleHcurlFirstFamily();
    if (type_basis == NODAL_GAUSS)
      Fb_tri->type_nodal_basis = Fb_tri->NODAL_GAUSS;
    else if (type_basis == NODAL_LOBATTO)
      Fb_tri->type_nodal_basis = Fb_tri->NODAL_LOBATTO;
    else if (type_basis == NODAL_REGULAR)
      Fb_tri->type_nodal_basis = Fb_tri->NODAL_REGULAR;
    
    Fb_tri->ConstructFiniteElement(order);
    element_tri_surf = Fb_tri;

    VectR2 points_dof2d = Fb_tri->PointsDofND();
    this->SetPointsDof2D_tri(points_dof2d);    

    this->num_dof_points_surf.Reallocate(4);
    for (int n = 0; n < 4; n++)
      {
	this->num_dof_points_surf(n).Reallocate(this->GetNbDofBoundary(n));
	for (int i = 0; i < this->GetNbDofBoundary(n); i++)
	  this->num_dof_points_surf(n)(i) = FacesDof(i, n);
      }    
  }
  
  
  //! construction of basis functions
  void TetrahedronHcurlFirstFamily::ConstructFunctions()
  {
    nb_dof_tri = order*3 + order*(order-1);
    nb_dof_quad = 0;
    
    VectReal_wp points_dof1d; VectR3 points_dof3d;
    if (type_basis != NODAL_REGULAR)
      {
        // points where dofs are defined, and tangent on these points
        VectR3 phi, points_inside;
	tangente_dof.Reallocate(nb_dof_loc);
        VectR2 points_tri, points_quad;
        VectReal_wp points_lob, weights_lob;
        Matrix<int> NumTri, NumQuad; Array3D<int> NumNodes;        
        ComputeGaussLobatto(points_lob, weights_lob, order+1);
	
        if (type_basis == NODAL_GAUSS)
          ComputeGaussLegendre(points_dof1d, weights_lob, order-1);
        else
          {
            points_dof1d.Reallocate(order);
            for (int i = 0; i < order; i++)
              points_dof1d(i) = points_lob(i+1);
          }
        
        Fb_geom.ConstructLobattoPoints(order+1, points_lob, points_tri, points_inside);
        
        points_dof3d.Reallocate(nb_dof_loc);
        // dofs on edges
        for (int i = 0; i < order; i++)
          {
            points_dof3d(i).Init(points_dof1d(i), 0, 0);
            tangente_dof(i).Init(1, 0, 0);

            points_dof3d(order+i).Init(0, points_dof1d(i), 0);
            tangente_dof(order+i).Init(0, 1, 0);
            
            points_dof3d(2*order+i).Init(0, 0, points_dof1d(i));
            tangente_dof(2*order+i).Init(0, 0, 1);
            
            points_dof3d(3*order+i).Init(1.0-points_dof1d(i), points_dof1d(i), 0);
            tangente_dof(3*order+i).Init(-1, 1, 0);
            
            points_dof3d(4*order+i).Init(1.0-points_dof1d(i), 0, points_dof1d(i));
            tangente_dof(4*order+i).Init(-1, 0, 1);
            
            points_dof3d(5*order+i).Init(0, 1.0-points_dof1d(i), points_dof1d(i));
            tangente_dof(5*order+i).Init(0, -1, 1);
          }
        
        // dofs on faces
        int offset = 6*order; int Ntri = (order-1)*order;
        for (int i = 0; i < (order-1)*order/2; i++)
          {
            Real_wp x = points_tri(3*(order+1) + i)(0);
            Real_wp y = points_tri(3*(order+1) + i)(1);
            points_dof3d(offset).Init(x, y, 0);
            points_dof3d(offset + 1).Init(x, y, 0);
            
            tangente_dof(offset).Init(1, 0, 0);
            tangente_dof(offset + 1).Init(0, 1, 0);

            points_dof3d(Ntri+offset).Init(x, 0, y);
            points_dof3d(Ntri+offset + 1).Init(x, 0, y);
            
            tangente_dof(Ntri+offset).Init(1, 0, 0);
            tangente_dof(Ntri+offset + 1).Init(0, 0, 1);
            
            points_dof3d(2*Ntri+offset).Init(0, x, y);
            points_dof3d(2*Ntri+offset + 1).Init(0, x, y);
            
            tangente_dof(2*Ntri+offset).Init(0, 1, 0);
            tangente_dof(2*Ntri+offset + 1).Init(0, 0, 1);

            points_dof3d(3*Ntri+offset).Init(1.0-x-y, x, y);
            points_dof3d(3*Ntri+offset + 1).Init(1.0-x-y, x, y);
            
            tangente_dof(3*Ntri+offset).Init(-1, 1, 0);
            tangente_dof(3*Ntri+offset + 1).Init(-1, 0, 1);
            
            offset += 2;
          }
        
        // dofs inside
        offset += 3*Ntri;
        for (int i = 0; i < (order-2)*(order-1)*order/6; i++)
          {
            R3 pt = points_inside(2*(order+1)*(order+1)+2+i);
            points_dof3d(offset) = pt;
            points_dof3d(offset+1) = pt;
            points_dof3d(offset+2) = pt;
            
            tangente_dof(offset).Init(1, 0, 0);
            tangente_dof(offset+1).Init(0, 1, 0);
            tangente_dof(offset+2).Init(0, 0, 1);
            offset += 3;
          }

        Matrix<Real_wp> Mh(nb_dof_loc, nb_dof_loc);
        Mh.Fill(0);
        for (int k = 0; k < nb_points_quadrature_inside; k++)
          {
            ComputeValuesPhiOrthoRef(this->PointsND(k), phi);
            for (int i = 0; i < nb_dof_loc; i++)
              for (int j = 0; j < nb_dof_loc; j++)
                Mh(i, j) += this->WeightsND(k)*DotProd(phi(i), phi(j));
          }
        
        //Mh.WriteText("MassTetOpt"+to_str(order)+".dat");

        //DISP(points_dof3d); DISP(tangente_dof);        
        // computation of VDM = psi_j(xi_k) \cdot t_k
        InverseVDM_Edge.Reallocate(nb_dof_loc, nb_dof_loc);
        for (int k = 0; k < nb_dof_loc; k++)
          {
            ComputeValuesPhiOrthoRef(points_dof3d(k), phi);
            for (int j = 0; j < nb_dof_loc; j++)
              InverseVDM_Edge(j, k) = DotProd(phi(j), tangente_dof(k));
          }
        
        //InverseVDM_Edge.WriteText("VdmTet"+to_str(order)+".dat");
        GetInverse(InverseVDM_Edge);
        
        // then checking phi_i(xi_k) \cdot t_k = delta_{i, k}
        for (int k = 0; k < nb_dof_loc; k++)
          {
            ComputeValuesPhiRef(points_dof3d(k), phi);            
            for (int j = 0; j < nb_dof_loc; j++)
              {
                Real_wp val = DotProd(phi(j), tangente_dof(k));
                Real_wp val_exact = 0;
                if (j == k)
                  val_exact = 1.0;
                
                if (abs(val - val_exact) > 1e4*epsilon_machine)
                  {
                    DISP(val); DISP(val_exact);
                    abort();
                  }
              }
          }
      }
    else
      {
        FacesDof.Reallocate(nb_dof_tri, 4);
        CoordinateDofs.Reallocate(nb_dof_loc, 5);
        points_dof3d.Reallocate(nb_dof_loc);
        VectReal_wp sm(order+2);
        for (int i = 0; i <= (order+1); i++)
          sm(i) = Real_wp(i)/Real_wp(order+1);
        
        // dofs on the first edge
        nb_dof_loc = 0;
        for (int i = 0; i < order; i++)
          {
            FacesDof(i, 0) = nb_dof_loc;
            FacesDof(i, 1) = nb_dof_loc;
            CoordinateDofs(nb_dof_loc, 0) = order-i;
            CoordinateDofs(nb_dof_loc, 1) = i+1;
            CoordinateDofs(nb_dof_loc, 2) = 0;
            CoordinateDofs(nb_dof_loc, 3) = 0;
            CoordinateDofs(nb_dof_loc, 4) = 0;
            points_dof3d(nb_dof_loc).Init(sm(i+1), 0.0, 0.0);
            nb_dof_loc++;
          }
        
        // second edge
        for (int i = 0; i < order; i++)
          {
            FacesDof(2*order+i,0) = nb_dof_loc;
            FacesDof(i,2) = nb_dof_loc;
            CoordinateDofs(nb_dof_loc,0) = order-i;
            CoordinateDofs(nb_dof_loc,1) = 0;
            CoordinateDofs(nb_dof_loc,2) = i+1;
            CoordinateDofs(nb_dof_loc,3) = 0;
            CoordinateDofs(nb_dof_loc,4) = 1;
            points_dof3d(nb_dof_loc).Init(0.0, sm(i+1), 0.0);
            nb_dof_loc++;
          }
        
        // third edge
        for (int i = 0; i < order; i++)
          {
            FacesDof(2*order+i,1) = nb_dof_loc;
            FacesDof(2*order+i,2) = nb_dof_loc;
            CoordinateDofs(nb_dof_loc,0) = order-i;
            CoordinateDofs(nb_dof_loc,1) = 0;
            CoordinateDofs(nb_dof_loc,2) = 0;
            CoordinateDofs(nb_dof_loc,3) = i+1;
            CoordinateDofs(nb_dof_loc,4) = 2;
            points_dof3d(nb_dof_loc).Init(0.0, 0.0, sm(i+1));
            nb_dof_loc++;
          }
        
        // fourth edge
        for (int i = 0; i < order; i++)
          {
            FacesDof(order+i,0) = nb_dof_loc;
            FacesDof(i,3) = nb_dof_loc;
            CoordinateDofs(nb_dof_loc,0) = 0;
            CoordinateDofs(nb_dof_loc,1) = order-i;
            CoordinateDofs(nb_dof_loc,2) = i+1;
            CoordinateDofs(nb_dof_loc,3) = 0;
            CoordinateDofs(nb_dof_loc,4) = 3;
            points_dof3d(nb_dof_loc).Init(sm(order-i), sm(i+1), 0.0);
            nb_dof_loc++;
          }
        
        // fifth edge
        for (int i = 0; i < order; i++)
          {
            FacesDof(order+i,1) = nb_dof_loc;
            FacesDof(2*order+i,3) = nb_dof_loc;
            CoordinateDofs(nb_dof_loc,0) = 0;
            CoordinateDofs(nb_dof_loc,1) = order-i;
            CoordinateDofs(nb_dof_loc,2) = 0;
            CoordinateDofs(nb_dof_loc,3) = i+1;
            CoordinateDofs(nb_dof_loc,4) = 4;
            points_dof3d(nb_dof_loc).Init(sm(order-i), 0.0, sm(i+1));
            nb_dof_loc++;
          }
        
        // sixth edge
        for (int i = 0; i < order; i++)
          {
            FacesDof(order+i,2) = nb_dof_loc;
            FacesDof(order+i,3) = nb_dof_loc;
            CoordinateDofs(nb_dof_loc,0) = 0;
            CoordinateDofs(nb_dof_loc,1) = 0;
            CoordinateDofs(nb_dof_loc,2) = order-i;
            CoordinateDofs(nb_dof_loc,3) = i+1;
            CoordinateDofs(nb_dof_loc,4) = 5;
            points_dof3d(nb_dof_loc).Init(0.0, sm(order-i), sm(i+1));
            nb_dof_loc++;
          }
        
        // first face
        int num = 3*order;
        for (int i = 1; i < order; i++)
          for (int j = 1; j < order; j++)
            if ((i+j) < (order+1))
              {
                CoordinateDofs(nb_dof_loc,0) = order+1-i-j;
                CoordinateDofs(nb_dof_loc,1) = i;
                CoordinateDofs(nb_dof_loc,2) = j;
                CoordinateDofs(nb_dof_loc,3) = 0;
                CoordinateDofs(nb_dof_loc,4) = 0;
                FacesDof(num++,0) = nb_dof_loc;
                points_dof3d(nb_dof_loc).Init(sm(i), sm(j), 0.0);
                nb_dof_loc++;
                CoordinateDofs(nb_dof_loc,0) = order+1-i-j;
                CoordinateDofs(nb_dof_loc,1) = i;
                CoordinateDofs(nb_dof_loc,2) = j;
                CoordinateDofs(nb_dof_loc,3) = 0;
                CoordinateDofs(nb_dof_loc,4) = 3;
                FacesDof(num++,0) = nb_dof_loc;
                points_dof3d(nb_dof_loc).Init(sm(i), sm(j), 0.0);
                nb_dof_loc++;
              }
        
        // second face
        num = 3*order;
        for (int i = 1; i < order; i++)
          for (int j = 1; j < order; j++)
            if ((i+j) < (order+1))
              {
                CoordinateDofs(nb_dof_loc,0) = order+1-i-j;
                CoordinateDofs(nb_dof_loc,1) = i;
                CoordinateDofs(nb_dof_loc,2) = 0;
                CoordinateDofs(nb_dof_loc,3) = j;
                CoordinateDofs(nb_dof_loc,4) = 0;
                FacesDof(num++,1) = nb_dof_loc;
                points_dof3d(nb_dof_loc).Init(sm(i), 0.0, sm(j));
                nb_dof_loc++;
                CoordinateDofs(nb_dof_loc,0) = order+1-i-j;
                CoordinateDofs(nb_dof_loc,1) = i;
                CoordinateDofs(nb_dof_loc,2) = 0;
                CoordinateDofs(nb_dof_loc,3) = j;
                CoordinateDofs(nb_dof_loc,4) = 4;
                FacesDof(num++,1) = nb_dof_loc;
                points_dof3d(nb_dof_loc).Init(sm(i), 0.0, sm(j));
                nb_dof_loc++;
              }
        
        // third face
        num = 3*order;
        for (int i = 1; i < order; i++)
          for (int j = 1; j < order; j++)
            if ((i+j) < (order+1))
              {
                CoordinateDofs(nb_dof_loc,0) = order+1-i-j;
                CoordinateDofs(nb_dof_loc,1) = 0;
                CoordinateDofs(nb_dof_loc,2) = i;
                CoordinateDofs(nb_dof_loc,3) = j;
                CoordinateDofs(nb_dof_loc,4) = 1;
                FacesDof(num++,2) = nb_dof_loc;
                points_dof3d(nb_dof_loc).Init(0.0, sm(i), sm(j));
                nb_dof_loc++;
                CoordinateDofs(nb_dof_loc,0) = order+1-i-j;
                CoordinateDofs(nb_dof_loc,1) = 0;
                CoordinateDofs(nb_dof_loc,2) = i;
                CoordinateDofs(nb_dof_loc,3) = j;
                CoordinateDofs(nb_dof_loc,4) = 5;
                FacesDof(num++,2) = nb_dof_loc;
                points_dof3d(nb_dof_loc).Init(0.0, sm(i), sm(j));
                nb_dof_loc++;
              }
        
        // fourth face
        num = 3*order;
        for (int i = 1; i < order; i++)
          for (int j = 1; j < order; j++)
            if ((i+j) < (order+1))
              {
                CoordinateDofs(nb_dof_loc,0) = 0;
                CoordinateDofs(nb_dof_loc,1) = order+1-i-j;
                CoordinateDofs(nb_dof_loc,2) = i;
                CoordinateDofs(nb_dof_loc,3) = j;
                CoordinateDofs(nb_dof_loc,4) = 3;
                FacesDof(num++,3) = nb_dof_loc;
                points_dof3d(nb_dof_loc).Init(sm(order+1-i-j), sm(i), sm(j));
                nb_dof_loc++;
                CoordinateDofs(nb_dof_loc,0) = 0;
                CoordinateDofs(nb_dof_loc,1) = order+1-i-j;
                CoordinateDofs(nb_dof_loc,2) = i;
                CoordinateDofs(nb_dof_loc,3) = j;
                CoordinateDofs(nb_dof_loc,4) = 5;
                FacesDof(num++,3) = nb_dof_loc;
                points_dof3d(nb_dof_loc).Init(sm(order+1-i-j), sm(i), sm(j));
                nb_dof_loc++;
              }
        
        // dofs inside the tetrahedron
        for (int i = 1; i < order; i++)
          for (int j = 1; j < order; j++)
            for (int k = 1; k < order; k++)
              if ((i+j+k) < (order+1))
                {
                  CoordinateDofs(nb_dof_loc,0) = order+1-i-j-k;
                  CoordinateDofs(nb_dof_loc,1) = i;
                  CoordinateDofs(nb_dof_loc,2) = j;
                  CoordinateDofs(nb_dof_loc,3) = k;
                  CoordinateDofs(nb_dof_loc,4) = 0;
                  points_dof3d(nb_dof_loc).Init(sm(i), sm(j), sm(k));
                  nb_dof_loc++;
                  CoordinateDofs(nb_dof_loc,0) = order+1-i-j-k;
                  CoordinateDofs(nb_dof_loc,1) = i;
                  CoordinateDofs(nb_dof_loc,2) = j;
                  CoordinateDofs(nb_dof_loc,3) = k;
                  CoordinateDofs(nb_dof_loc,4) = 1;
                  points_dof3d(nb_dof_loc).Init(sm(i), sm(j), sm(k));
                  nb_dof_loc++;
                  CoordinateDofs(nb_dof_loc,0) = order+1-i-j-k;
                  CoordinateDofs(nb_dof_loc,1) = i;
                  CoordinateDofs(nb_dof_loc,2) = j;
                  CoordinateDofs(nb_dof_loc,3) = k;
                  CoordinateDofs(nb_dof_loc,4) = 2;
                  points_dof3d(nb_dof_loc).Init(sm(i), sm(j), sm(k));
                  nb_dof_loc++;
                }
        
        points_dof3d = this->PointsNodalND();
        
        // fonctions on polynoms
        Real_wp one(1), zero(0), coef;
        Real_wp rp1 = order+1;
        MultivariatePolynomial<Real_wp> L0(3,1), L1(3,1), L2(3,1),
          L3(3,1), P, Q, R, Omega_X, Omega_Y, Omega_Z, Zero;
        L0(0) = Real_wp(1); L0(1) = Real_wp(-1); L0(2) = Real_wp(-1); L0(3) = Real_wp(-1);
        L1(0) = Real_wp(0); L1(1) = Real_wp(0); L1(2) = Real_wp(0); L1(3) = Real_wp(1);
        L2(0) = Real_wp(0); L2(1) = Real_wp(0); L2(2) = Real_wp(1); L2(3) = Real_wp(0);
        L3(0) = Real_wp(0); L3(1) = Real_wp(1); L3(2) = Real_wp(0); L3(3) = Real_wp(0);
        Zero.SetOrder(3,0); Zero(0) = zero;
        
        Poly_Phi_X.Reallocate(nb_dof_loc);
        Poly_Phi_Y.Reallocate(nb_dof_loc);
        Poly_Phi_Z.Reallocate(nb_dof_loc);
        for (int node = 0; node < nb_dof_loc; node++)
          {    
            int  i = CoordinateDofs(node,0);
            int  j = CoordinateDofs(node,1);
            int  k = CoordinateDofs(node,2);
            int  l = CoordinateDofs(node,3);
            int num_edge = CoordinateDofs(node,4);
            
            if (num_edge == 0)
              {
                Omega_X = L0; Add(one, L1, Omega_X);
                Omega_Y = L1; Omega_Z = L1;
                R.SetOrder(3,0);
                R(0) = one;
                for (int m = 1; m < i; m++)
                  {
                    P = L0;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 1; m < j; m++)
                  {
                    P = L1;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 0; m < k; m++)
                  {
                    P = L2;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 0; m < l; m++)
                  {
                    P = L3;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
              }
            else if (num_edge == 1)
              {
                Omega_Y = L0; Add(one, L2, Omega_Y);
                Omega_X = L2; Omega_Z = L2;
                R.SetOrder(3,0);
                R(0) = one;
                for (int m = 1; m < i; m++)
                  {
                    P = L0;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 0; m < j; m++)
                  {
                    P = L1;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 1; m < k; m++)
                  {
                    P = L2;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 0; m < l; m++)
                  {
                    P = L3;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
              }
            else if (num_edge == 2)
              {
                Omega_Z = L0; Add(one, L3, Omega_Z);
                Omega_X = L3; Omega_Y = L3;
                R.SetOrder(3,0);
                R(0) = one;
                for (int m = 1; m < i; m++)
                  {
                    P = L0;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 0; m < j; m++)
                  {
                    P = L1;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 0; m < k; m++)
                  {
                    P = L2;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 1; m < l; m++)
                  {
                    P = L3;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
              }
            else if (num_edge == 3)
              {
                Omega_X = L2; Mlt(-one, Omega_X);
                Omega_Y = L1; Omega_Z = Zero;
                R.SetOrder(3,0);
                R(0) = one;
                for (int m = 0; m < i; m++)
                  {
                    P = L0;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 1; m < j; m++)
                  {
                    P = L1;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 1; m < k; m++)
                  {
                    P = L2;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 0; m < l; m++)
                  {
                    P = L3;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
              }
            else if (num_edge == 4)
              {
                Omega_Z = L1;
                Omega_X = L3; Omega_Y = Zero; Mlt(Real_wp(-1),Omega_X);
                R.SetOrder(3,0);
                R(0) = one;
                for (int m = 0; m < i; m++)
                  {
                    P = L0;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 1; m < j; m++)
                  {
                    P = L1;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 0; m < k; m++)
                  {
                    P = L2;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 1; m < l; m++)
                  {
                    P = L3;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
              }
            else if (num_edge == 5)
              {
                Omega_Y = L3; Mlt(-one, Omega_Y);
                Omega_X = Zero; Omega_Z = L2;
                R.SetOrder(3,0);
                R(0) = one;
                for (int m = 0; m < i; m++)
                  {
                    P = L0;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 0; m < j; m++)
                  {
                    P = L1;
                    coef = one/Real_wp(m+1);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 1; m < k; m++)
                  {
                    P = L2;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
                
                for (int m = 1; m < l; m++)
                  {
                    P = L3;
                    coef = one/Real_wp(m);
                    Mlt(Real_wp(rp1), P); P(0) -= Real_wp(m);
                    Q = R;
                    MltAdd(coef, P, Q, zero, R);
                  }
              }
            
            // DISP(node); DISP(Omega_X); DISP(Omega_Y); DISP(Omega_Z); 
            Poly_Phi_X(node).SetOrder(3,0);
            MltAdd(one, Omega_X, R, zero, Poly_Phi_X(node));
            Poly_Phi_Y(node).SetOrder(3,0);
            MltAdd(one, Omega_Y, R, zero, Poly_Phi_Y(node));
            Poly_Phi_Z(node).SetOrder(3,0);
            MltAdd(one, Omega_Z, R, zero, Poly_Phi_Z(node));
            // DISP(node); DISP(Poly_Phi_X(node));
            // DISP(Poly_Phi_Y(node)); DISP(Poly_Phi_Z(node));
          }
        
        Poly_CurlPhi_X.Reallocate(nb_dof_loc);
        Poly_CurlPhi_Y.Reallocate(nb_dof_loc);
        Poly_CurlPhi_Z.Reallocate(nb_dof_loc);
        for (int node = 0; node < nb_dof_loc; node++)
          {
            DerivatePolynomial(Poly_Phi_Z(node), P, 1);
            DerivatePolynomial(Poly_Phi_Y(node), Q, 2);
            Add(Real_wp(-one), Q, P);
            Poly_CurlPhi_X(node) = P;
            
            DerivatePolynomial(Poly_Phi_X(node), P, 2);
            DerivatePolynomial(Poly_Phi_Z(node), Q, 0);
            Add(Real_wp(-one), Q, P);
            Poly_CurlPhi_Y(node) = P;
            
            DerivatePolynomial(Poly_Phi_Y(node), P, 0);
            DerivatePolynomial(Poly_Phi_X(node), Q, 1);	
            Add(Real_wp(-one), Q, P);
            Poly_CurlPhi_Z(node) = P;
          }
        
      }

    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof3d);
  }
   
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void TetrahedronHcurlFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    if (type_basis != NODAL_REGULAR)
      {
        VectR3 psi;
        ComputeValuesPhiOrthoRef(point_loc, psi);
        FillZero(res);
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = 0; j < nb_dof_loc; j++)
            Add(InverseVDM_Edge(i, j), psi(j), res(i));
        
      }
    else
      {
        for (int node = 0; node < nb_dof_loc; node++)
          {
            res(node)(0) = Poly_Phi_X(node).Evaluate(point_loc);
            res(node)(1) = Poly_Phi_Y(node).Evaluate(point_loc);
            res(node)(2) = Poly_Phi_Z(node).Evaluate(point_loc);
          }
      }
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void TetrahedronHcurlFirstFamily::ComputeCurlPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    if (type_basis != NODAL_REGULAR)
      {
        VectR3 curl_psi;
        ComputeCurlPhiOrthoRef(point_loc, curl_psi);
        FillZero(res);
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = 0; j < nb_dof_loc; j++)
            res(i) += InverseVDM_Edge(i, j)*curl_psi(j);
        
      }
    else
      {
        for (int node = 0; node < nb_dof_loc; node++)
          {
            res(node)(0) = Poly_CurlPhi_X(node).Evaluate(point_loc);
            res(node)(1) = Poly_CurlPhi_Y(node).Evaluate(point_loc);
            res(node)(2) = Poly_CurlPhi_Z(node).Evaluate(point_loc);
          }
      }
  }


  //! Evaluating nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of nearly orthogonal functions on pointloc
  */
  void TetrahedronHcurlFirstFamily::
  ComputeValuesPhiOrthoRef(const R3& pointloc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    int r = order;
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    // coordinates on the cube [-1, 1]
    Real_wp a(-1), b(-1), c(-1);
    if (abs(1.0-y-z) > epsilon_machine)
      a = 2.0*x/(1.0-y-z) - 1.0;
    
    if (abs(1.0-z) > epsilon_machine)
      b = 2.0*y/(1.0-z) - 1.0;
    
    c = 2.0*z - 1.0;

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& OddJacobiPolynom = Fb_geom.GetOddJacobiPolynomial();
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();

    VectReal_wp Px;
    EvaluateJacobiPolynomial(LegendrePolynom, r, a, Px);
    Vector<VectReal_wp> Py(r+1);
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(OddJacobiPolynom(i), r-i, b, Py(i));
    
    Vector<VectReal_wp> Pz(r+1);
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, c, Pz(i));
    
    Real_wp valx(0), valy(0), vloc;
    Real_wp pow_oneMinusYZ = 1.0, pow_oneMinusZ = 1.0;
    int num = 0;
    for (int i = 0; i < r; i++)
      {
        // computation of valx = L_i(a) (1-y-z)^i
        if (i == 0)
          valx = 1.0;
        else if (i == 1)
          valx = 2.0*x + y + z - 1.0;
        else
          valx = pow_oneMinusYZ*Px(i);
        
        pow_oneMinusZ = 1.0;
        for (int j = 0; j < r-i; j++)
          {
            // computation of valy = P_j^(2i+1)(a) (1-z)^j
            if (j == 0)
              valy = 1.0;
            else if (j == 1)
              valy = 2.0*y + z - 1.0 - OddJacobiPolynom(i)(0,0)*(1.0-z);
            else
              valy = pow_oneMinusZ*Py(i)(j);
            
            for (int k = 0; k < r-i-j; k++)
              {
                vloc = valx*valy*Pz(i+j)(k);
                // part due to P_r-1^3
                phi(num).Init(vloc, 0, 0);
                phi(num+1).Init(0, vloc, 0);
                phi(num+2).Init(0, 0, vloc);
                num += 3;
              }
            
            // when vloc = x^i y^j z^r-1-i-j
            // part due to (0, vloc z, -vloc y) and (vloc z, 0, -vloc x)
            vloc = valx*valy*Pz(i+j)(r-1-i-j);
            phi(num).Init(0, vloc*z, -vloc*y);
            phi(num+1).Init(vloc*z, 0, -vloc*x);
            num += 2;
            
            pow_oneMinusZ *= 1.0-z;
          }
        
        // part due to (x^i y^r-i, -x^i+1 y^r-1-i, 0)
        vloc = valx*valy;
        phi(num).Init(vloc*y, -vloc*x, 0);
        num++;
            
        pow_oneMinusYZ *= 1.0-y-z;
      }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < phi.GetM(); i++)
      phi(i) *= InvWeightBasisRr(i);
    
  }

  
  //! Evaluating curl of nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] curl_phi curl of nearly orthogonal functions on pointloc
  */
  void TetrahedronHcurlFirstFamily::
  ComputeCurlPhiOrthoRef(const R3& pointloc, VectR3& curl_phi) const
  {
    // orthogonal polynomials :
    // L_i(2x/(1-y-z) - 1) P_j^{2i+1}(2y/(1-z)-1) (1-y)^i P_k^{2i+2j+2}(2z-1) (1-z)^(i+j)
    // = L_i(a) P_j^(2i+1)(b) ((1-b)/2)^i P_k^(2i+2j+2)(c) ((1-c)/2)^(i+j)
    // where (a, b, c) is on the cube [-1,1]^3 and (x,y,z) on the unit tetrahedron
    
    // derivative with respect to a, b, c is equal to :
    // L'_i(a) P_j^(2i+1)(b) ((1-b)/2)^i P_k^(2i+2j+2)(c) ((1-c)/2)^(i+j)
    // L_i(a) ((1-b)/2)^(i-1) ( (1-b)/2 P'_j^(2i+1)(b) - i/2 P_j^(2i+1)(b) )
    //      P_k^(2i+2j+2)(c) ((1-c)/2)^(i+j)
    // L_i(a) P_j^(2i+1)(b) ((1-b)/2)^i ((1-c)/2)^(i+j-1) 
    //    ( (1-c)/2 P'_k^(2i+2j+2)(c) - (i+j)/2 P_k^(2i+2j+2)(c)) 
    curl_phi.Reallocate(nb_dof_loc);
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    // coordinates on the cube [-1, 1]
    Real_wp a(-1), b(-1), c(-1);
    Real_wp da_dx(0), da_dy(0), db_dy(0), db_dz(0);
    if (abs(1.0-y-z) > epsilon_machine)
      {
	a = 2.0*x/(1.0-y-z) - 1.0;
	da_dx = 2.0/(1.0-y-z);
	da_dy = 2.0*x/square(1.0-y-z);
      }
    
    if (abs(1.0-z) > epsilon_machine)
      {
	b = 2.0*y/(1.0-z) - 1.0;
	db_dy = 2.0/(1.0-z);
	db_dz = 2.0*y/square(1.0-z);
      }
    
    c = 2.0*z - 1.0;

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& OddJacobiPolynom = Fb_geom.GetOddJacobiPolynomial();
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    
    VectReal_wp Px, dPx;
    int r = order;
    EvaluateJacobiPolynomial(LegendrePolynom, r, a, Px, dPx);
    Vector<VectReal_wp> Py(r+1), dPy(r+1);
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(OddJacobiPolynom(i), r-i, b, Py(i), dPy(i));
    
    Vector<VectReal_wp> Pz(r+1), dPz(r+1);
    for (int i = 0; i <= order; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, c, Pz(i), dPz(i));
    
    Real_wp valx(0), valy(0), valz(0), phi(0); R3 grad_phi;
    Real_wp dvalx_dx(0), dvalx_dy(0), dvaly_dy(0), dvaly_dz(0), dvalz_dz(0);
    Real_wp pow_oneMinusYZ = 1.0, pow_oneMinusZ = 1.0;    
    Real_wp powm1_oneMinusYZ = 0.0, powm1_oneMinusZ = 0.0;    
    int num = 0;
    for (int i = 0; i < r; i++)
      {
        if (i == 0)
          {
            valx = 1.0;
            dvalx_dx = 0;
            dvalx_dy = 0;
          }
        else if (i == 1)
          {
            valx = 2.0*x + y + z - 1.0;
            dvalx_dx = 2.0;
            dvalx_dy = 1.0;
          }
        else
          {
            valx = pow_oneMinusYZ*Px(i);
            dvalx_dx = pow_oneMinusYZ*dPx(i)*da_dx;
            dvalx_dy = pow_oneMinusYZ*dPx(i)*da_dy - i*powm1_oneMinusYZ*Px(i);
          }
        
        powm1_oneMinusZ = 0; pow_oneMinusZ = 1.0;
        for (int j = 0; j < r-i; j++)
          {
            if (j == 0)
              {
                valy = 1.0;
                dvaly_dy = 0;
                dvaly_dz = 0;
              }
            else if (j == 1)
              {
                valy = 2.0*y + z - 1.0 - OddJacobiPolynom(i)(0,0)*(1.0-z);
                dvaly_dy = 2.0;
                dvaly_dz = 1.0 + OddJacobiPolynom(i)(0,0);
              }
            else
              {
                valy = pow_oneMinusZ*Py(i)(j);
                dvaly_dy = pow_oneMinusZ*dPy(i)(j)*db_dy;
                dvaly_dz = pow_oneMinusZ*dPy(i)(j)*db_dz - j*powm1_oneMinusZ*Py(i)(j);
              }
                   
            for (int k = 0; k < r-i-j; k++)
              {
                valz = Pz(i+j)(k);
                dvalz_dz = 2.0*dPz(i+j)(k);
                phi = valx*valy*valz;
                grad_phi.Init(dvalx_dx*valy*valz, (dvalx_dy*valy + valx*dvaly_dy)*valz, 
                              dvalx_dy*valy*valz + valx*dvaly_dz*valz + valx*valy*dvalz_dz);
                
                // part due to P_r-1^3
                curl_phi(num).Init(0, grad_phi(2), -grad_phi(1));
                curl_phi(num+1).Init(-grad_phi(2), 0, grad_phi(0));
                curl_phi(num+2).Init(grad_phi(1), -grad_phi(0), 0);
                num += 3;
              }
            
            // when vloc = x^i y^j z^r-1-i-j
            // part due to (0, vloc z, -vloc y) and (vloc z, 0, -vloc x)
            curl_phi(num).Init(-grad_phi(1)*y -grad_phi(2)*z - 2.0*phi,
                               grad_phi(0)*y, grad_phi(0)*z);
            curl_phi(num+1).Init(-grad_phi(1)*x, grad_phi(2)*z + 2.0*phi + grad_phi(0)*x,
                                 -grad_phi(1)*z);
            num += 2;
            
            powm1_oneMinusZ = pow_oneMinusZ;
            pow_oneMinusZ *= 1.0-z;
          }
        
        curl_phi(num).Init((valx*dvaly_dz+dvalx_dy*valy)*x, (valx*dvaly_dz+dvalx_dy*valy)*y,
                           -dvalx_dx*valy*x - (dvalx_dy*valy + valx*dvaly_dy)*y - 2.0*valx*valy);
        num++;
        
        powm1_oneMinusYZ = pow_oneMinusYZ;
	pow_oneMinusYZ *= 1.0-y-z;
      }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < nb_dof_loc; i++)
      curl_phi(i) *= InvWeightBasisRr(i);

  }
  
  
  //! displays details of class TetrahedronHcurlFirstFamily
  ostream& operator <<(ostream& out, const TetrahedronHcurlFirstFamily& e)
  {
    out<<static_cast<const TetrahedronReference<2>&>(e);
    return  out;
  }

  
  //! default constructor
  TetrahedronHcurlOptimalFirstFamily::TetrahedronHcurlOptimalFirstFamily()
    : TetrahedronHcurlFirstFamily()
  {
    type_basis = NODAL_LOBATTO;
  }


  //! displays details of class TetrahedronHcurlOptimalFirstFamily
  ostream& operator <<(ostream& out, const TetrahedronHcurlOptimalFirstFamily& e)
  {
    out<<static_cast<const TetrahedronHcurlFirstFamily&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_TETRAHEDRON_HCURL_FIRST_FAMILY_CXX
#endif
