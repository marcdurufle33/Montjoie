#ifndef MONTJOIE_FILE_PYRAMID_HCURL_OTHER_CXX

namespace Montjoie
{
  // File with all pyramidal functions described in the litterature
  
  //! default constructor
  PyramidHcurlOther::PyramidHcurlOther() : PyramidReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    // type_basis = GRADINARU;
    // type_basis = DOUCET;
    //type_basis = GRAGLIA;
    type_basis = ZGAINSKI;
    //type_basis = NIGAM_PHILLIPS_2;
    //type_basis = ZAGLMAYR;
    //type_basis = NIGAM_PHILLIPS_1;
  }
  
    
  //! how to number mesh
  void PyramidHcurlOther::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return PyramidReference<2>::ConstructNumberMap(nmap, dg);

    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    nmap.SetNbDofQuadrangle(order, 2*order*(order-1));
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order-1));
    
    // dofs inside the pyramid
    nmap.SetNbDofPyramid(order, nb_dof_loc-nb_dof_boundaries);

    if (type_basis == GRAGLIA)
      {
        nmap.SetOppositeEdgesDofSymmetry(order, order);
        nmap.SetAllEdgesDofToSkewSymmetric(order);
      }
    else
      nmap.SetEqualEdgesDofSymmetry(order, order+1);
    
    // sign changes on edges
    FindHcurlSignEdge(nmap);
    
    // rotation of dofs on faces
    FindHcurlLinearCombinationRotation(nmap, *element_tri_surf, *element_quad_surf);
  }

    
  //! construction of finite element
  void PyramidHcurlOther::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						 int rsurf_tri, int rsurf_quad,
						 int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    PyramidReference<2>::ConstructFiniteElement(r, rgeom, r+1, type_quad);
    
    // computation of basis functions
    ConstructFunctions();
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    // construction of Value_PhiVec, Curl_Phi and elementary matrices
    ConstructHcurlElementaryMatrix();
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace(false);
    
    TriangleHcurlOther* Fb_tri = new TriangleHcurlOther();
    QuadrangleHcurlOther* Fb_quad = new QuadrangleHcurlOther();
    element_tri_surf = Fb_tri;
    element_quad_surf = Fb_quad;
    
    Fb_tri->type_basis = type_basis;
    Fb_quad->type_basis = type_basis;
    Fb_tri->ConstructFiniteElement(order);
    Fb_quad->ConstructFiniteElement(order);
    
  }
  
  
  //! construction of basis functions
  void PyramidHcurlOther::ConstructFunctions()
  {
    int r = order;
    int nb_dof_inside_tri = order*(order-1);
    int nb_dof_inside_quad = 2*order*(order-1);
    nb_dof_tri = order*3 + nb_dof_inside_tri;
    nb_dof_quad = order*4 + nb_dof_inside_quad;
    nb_dof_boundaries = 4*nb_dof_inside_tri + nb_dof_inside_quad + 8*order;
    nb_dof_loc = r*(r+1)*(2*r+1)/2 + r + (r-2)*(r-1) + 2*r*(r+1);
    
    if (type_basis == NIGAM_PHILLIPS_2)
      nb_dof_loc = order*(7 + 7*order + 2*order*order)/2;
    else if (type_basis == NIGAM_PHILLIPS_1)
      {
        //nb_dof_loc = order*(7 + 7*order + 2*order*order)/2;
        nb_dof_loc = 3*r*r*r + 5*r;
      }
    else if (type_basis == ZAGLMAYR)
      {
        nb_dof_loc = 0;
        for (int k = 0; k <= r-1; k++)
          for (int i = 0; i <= k; i++)
            for (int j = 0; j <= k; j++)
              nb_dof_loc += 3;
        
        for (int k = 1; k <= r; k++)
          for (int i = 1; i <= k; i++)
            for (int j = 1; j <= k; j++)
              nb_dof_loc++;

        for (int k = 0; k <= r-1; k++)
          for (int i = 0; i <= k; i++)
            for (int j = 0; j < k; j++)
              nb_dof_loc += 2;
      }
    else if (type_basis == GRAGLIA)
      {
        if (order == 2)
          nb_dof_loc++;
      }
    
    VectR3 points_dof3d;
    points_dof3d.Reallocate(nb_points_quadrature_inside);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      points_dof3d(i) = this->PointsND(i);

    this->SetPointsDofND(points_dof3d);
  }
  
  
  //! compute \f$ \varphi_{node} (point\_loc) \f$
  void PyramidHcurlOther::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    // on recupere x, y, z de la pyramide symetrique
    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);
    int r = order;
    switch (type_basis)
      {
      case ZGAINSKI :
        {
          Real_wp beta1 = 0.5*(1.0-x-z);
          Real_wp beta2 = 0.5*(1.0-y-z);
          Real_wp beta3 = 0.5*(1.0+x-z);
          Real_wp beta4 = 0.5*(1.0+y-z);
          
          R3 grad_beta1(-0.5, 0, -0.5);
          R3 grad_beta2(0, -0.5, -0.5);
          R3 grad_beta3(0.5, 0, -0.5);
          R3 grad_beta4(0, 0.5, -0.5);
          
          Real_wp N1 = beta1*beta2/(1.0-z);
          Real_wp N2 = beta2*beta3/(1.0-z);
          Real_wp N3 = beta3*beta4/(1.0-z);
          Real_wp N4 = beta4*beta1/(1.0-z);
          Real_wp N5 = z;
          
          R3 grad_N1 = grad_beta1*beta2 + beta1*grad_beta2;
          grad_N1 *= 1.0/(1.0-z); grad_N1(2) += beta1*beta2/square(1.0-z);
          R3 grad_N2 = grad_beta2*beta3 + beta2*grad_beta3;
          grad_N2 *= 1.0/(1.0-z); grad_N2(2) += beta2*beta3/square(1.0-z);
          R3 grad_N3 = grad_beta3*beta4 + beta3*grad_beta4;
          grad_N3 *= 1.0/(1.0-z); grad_N3(2) += beta3*beta4/square(1.0-z);
          R3 grad_N4 = grad_beta4*beta1 + beta4*grad_beta1;
          grad_N4 *= 1.0/(1.0-z); grad_N4(2) += beta4*beta1/square(1.0-z);
          R3 grad_N5(0, 0, 1);
          
          res(0) = N1*grad_beta3 - N2*grad_beta1;
          res(r) = N2*grad_beta4 - N3*grad_beta2;
          res(2*r) = N3*grad_beta1 - N4*grad_beta3; Mlt(-1.0, res(2*r));
          res(3*r) = N4*grad_beta2 - N1*grad_beta4; Mlt(-1.0, res(3*r));
          
          res(4*r) = N1*grad_N5 - N5*grad_N1;
          res(5*r) = N2*grad_N5 - N5*grad_N2;
          res(6*r) = N3*grad_N5 - N5*grad_N3;
          res(7*r) = N4*grad_N5 - N5*grad_N4;
          
          if (r >= 2)
            {
              res(1) = N1*grad_beta3 + N2*grad_beta1;
              res(r+1) = N2*grad_beta4 + N3*grad_beta2;
              res(2*r+1) = N3*grad_beta1 + N4*grad_beta3;
              res(3*r+1) = N4*grad_beta2 + N1*grad_beta4;

              res(4*r+1) = N1*grad_N5 + N5*grad_N1;
              res(5*r+1) = N2*grad_N5 + N5*grad_N2;
              res(6*r+1) = N3*grad_N5 + N5*grad_N3;
              res(7*r+1) = N4*grad_N5 + N5*grad_N4;
              
              res(16) = N1*N2*grad_beta4;
              res(17) = N2*N3*grad_beta1;
              res(18) = N3*N4*grad_beta2;
              res(19) = N4*N1*grad_beta3;
              
              res(20) = N1*N5*grad_beta3;
              res(21) = N2*N5*grad_beta1;
              
              res(22) = N2*N5*grad_beta4;
              res(23) = N3*N5*grad_beta2;

              res(25) = N3*N5*grad_beta1;
              res(24) = N4*N5*grad_beta3;
              
              res(27) = N4*N5*grad_beta2;
              res(26) = N1*N5*grad_beta4;
              
              res(28) = N1*N3*grad_N5;
            }
        }
        break;
      case DOUCET :
        {          
          // les expressions de l'article
          VectR3 phi(8);
          res(0).Init(1.0 , 0, 0);
          res(1).Init(0, 1.0, 0);
          res(2).Init(0, 0, 1.0);
          res(3).Init(y, 0, x*y/(1.0-z));
          res(4).Init(0, x, x*y/(1.0-z));
          res(5).Init(z, 0, -x);
          res(6).Init(0, z, -y);
          res(7).Init(y*z/(1.0-z), x*z/(1.0-z), x*y*z/square(1.0-z)-x*y/(1.0-z));
        }
        break;
      case GRAGLIA :
        {
          // on calcule xu, yu, zu de la pyramide unite
          Real_wp xu = 0.5*(1.0 + x - z);
          Real_wp yu = 0.5*(1.0 + y - z);
          Real_wp zu = z;
          
          // on calcule DF-1
          Matrix3_3 dfjm1;
          dfjm1(0, 0) = 0.5; dfjm1(1, 1) = 0.5; dfjm1(2, 2) = 1.0;
          dfjm1(0, 2) = -0.5; dfjm1(1, 2) = -0.5;
          
          Real_wp N1 = xu;
          Real_wp N2 = yu;
          Real_wp N3 = 1.0-zu-xu;
          Real_wp N4 = 1.0-zu-yu;
          Real_wp N5 = zu;
          Real_wp invOmz = 1.0/(1.0 - N5);
          Real_wp invOmz2 = invOmz*invOmz;
          
          R3 grad_N1(1.0, 0, 0);
          R3 grad_N2(0, 1.0, 0);
          R3 grad_N3(-1.0, 0, -1.0);
          R3 grad_N4(0, -1.0, -1.0);
          R3 grad_N5(0, 0, 1.0);
          
          VectR3 phi(nb_dof_loc);
          R3 v1, v2, v3, v4;
          R3 v5, v6, v7, v8;
          v5 = N3*N4*grad_N5*invOmz - N4*N5*grad_N3*invOmz
            - N3*N5*grad_N4*invOmz - N3*N4*N5*grad_N5*invOmz2;
          
          v6 = N1*N4*grad_N5*invOmz - N1*N5*grad_N4*invOmz
            - N4*N5*grad_N1*invOmz - N1*N4*N5*grad_N5*invOmz2;
          
          v7 = N1*N2*grad_N5*invOmz - N2*N5*grad_N1*invOmz
            - N1*N5*grad_N2*invOmz - N1*N2*N5*grad_N5*invOmz2;
          
          v8 = N2*N3*grad_N5*invOmz - N3*N5*grad_N2*invOmz
            - N2*N5*grad_N3*invOmz - N2*N3*N5*grad_N5*invOmz2;
          
          v4 = N2*N3*grad_N4*invOmz - N3*N4*grad_N2*invOmz; Mlt(-1.0, v4);
          v1 = N4*N3*grad_N1*invOmz - N1*N4*grad_N3*invOmz;
          v2 = N4*N1*grad_N2*invOmz - N1*N2*grad_N4*invOmz; 
          v3 = N2*N1*grad_N3*invOmz - N3*N2*grad_N1*invOmz; Mlt(-1.0, v3);
                    
          if (r > 1)
            {
              // dofs on edges
              phi(0) = v1; Mlt((3.0*N3 - 1.0)*(2.0*N4 + 2.0*N5 - 1.0), phi(0));
              phi(1) = v1; Mlt((3.0*N1 - 1.0)*(2.0*N4 + 2.0*N5 - 1.0), phi(1));
              
              phi(2) = v2; Mlt((3.0*N4 - 1.0)*(2.0*N1 + 2.0*N5 - 1.0), phi(2));
              phi(3) = v2; Mlt((3.0*N2 - 1.0)*(2.0*N1 + 2.0*N5 - 1.0), phi(3));

              phi(4) = v3; Mlt((3.0*N3 - 1.0)*(2.0*N2 + 2.0*N5 - 1.0), phi(4));
              phi(5) = v3; Mlt((3.0*N1 - 1.0)*(2.0*N2 + 2.0*N5 - 1.0), phi(5));

              phi(6) = v4; Mlt((3.0*N4 - 1.0)*(2.0*N3 + 2.0*N5 - 1.0), phi(6));
              phi(7) = v4; Mlt((3.0*N2 - 1.0)*(2.0*N3 + 2.0*N5 - 1.0), phi(7));

              phi(9) = v5; Mlt((3.0*N5 - 1.0), phi(9));
              phi(8) = v5; Mlt((3.0*N3 - 3.0*N2 - 1.0), phi(8));

              phi(11) = v6; Mlt((3.0*N5 - 1.0), phi(11));
              phi(10) = v6; Mlt((3.0*N4 - 3.0*N3 - 1.0), phi(10));

              phi(13) = v7; Mlt((3.0*N5 - 1.0), phi(13));
              phi(12) = v7; Mlt((3.0*N1 - 3.0*N4 - 1.0), phi(12));

              phi(15) = v8; Mlt((3.0*N5 - 1.0), phi(15));
              phi(14) = v8; Mlt((3.0*N2 - 3.0*N1 - 1.0), phi(14));
              
              // dofs on quadrilateral face
              phi(16) = v1; Mlt(N2*(3.0*N1 - 1.0), phi(16));
              phi(17) = v1; Mlt(N2*(3.0*N3 - 1.0), phi(17));
              
              phi(18) = v4; Mlt(N1*(3.0*N4 - 1.0), phi(18));
              phi(19) = v4; Mlt(N1*(3.0*N2 - 1.0), phi(19));
              
              // dofs on triangular faces
              phi(20) = v1; Mlt(3.0*N5*(1.0-3.0*N2), phi(20));
              phi(21) = v5; Mlt(3.0*N1, phi(21));
              
              phi(22) = v2; Mlt(3.0*N5*(1.0-3.0*N3), phi(22));
              phi(23) = v6; Mlt(3.0*N2, phi(23)); 
              
              phi(24) = v3; Mlt(3.0*N5*(1.0-3.0*N4), phi(24));
              phi(25) = v7; Mlt(3.0*N3, phi(25)); phi(25) += phi(24);
              
              phi(26) = v4; Mlt(3.0*N5*(1.0-3.0*N1), phi(26));
              phi(27) = v8; Mlt(3.0*N4, phi(27)); phi(27) += phi(26);
              
              // dofs inside
              phi(28) = v1; Mlt(9.0*N2*N5, phi(28));
              phi(29) = v4; Mlt(9.0*N1*N5, phi(29));
            }
          else
            {
              phi(0) = v1;
              phi(1) = v2;
              phi(2) = v3;
              phi(3) = v4;
              
              phi(4) = v5;
              phi(5) = v6;              
              phi(6) = v7;
              phi(7) = v8;
            }
          
          // on applique DFi*-1
          for (int k = 0; k < nb_dof_loc; k++)
            MltTrans(dfjm1, phi(k), res(k));
        }
        break;
      case GRADINARU :
        {          
          // on calcule xu, yu, zu de la pyramide unite
          Real_wp xu = 0.5*(1.0 + x - z);
          Real_wp yu = 0.5*(1.0 + y - z);
          Real_wp zu = z;
          
          // on calcule DF-1
          Matrix3_3 dfjm1;
          dfjm1(0, 0) = 0.5; dfjm1(1, 1) = 0.5; dfjm1(2, 2) = 1.0;
          dfjm1(0, 2) = -0.5; dfjm1(1, 2) = -0.5;
          
          // puis on met les expressions de l'article
          VectR3 phi(8);
          phi(0).Init(1.0 - zu - yu, 0, xu - xu*yu/(1.0-zu));
          phi(1).Init(0, xu, xu*yu/(1.0-zu));
          phi(2).Init(yu, 0, xu*yu/(1.0-zu));
          phi(3).Init(0, 1.0 - zu - xu, yu - xu*yu/(1.0-zu));
          phi(4).Init(zu - yu*zu/(1.0-zu), zu - xu*zu/(1.0-zu),
                      1.0 - xu - yu + xu*yu/(1.0-zu) - xu*yu*zu/square(1.0-zu) );
          phi(5).Init(-zu + yu*zu/(1.0-zu), xu*zu/(1.0-zu),
                      xu - xu*yu/(1.0-zu) + xu*yu*zu/square(1.0-zu) );
          phi(6).Init(yu*zu/(1.0-zu), -zu + xu*zu/(1.0-zu),
                      yu - xu*yu/(1.0-zu) + xu*yu*zu/square(1.0-zu) );
          phi(7).Init(-yu*zu/(1.0-zu), -xu*zu/(1.0-zu), xu*yu/(1.0-zu) - xu*yu*zu/square(1.0-zu));
          
          // on applique DFi*-1
          for (int k = 0; k < 8; k++)
            MltTrans(dfjm1, phi(k), res(k));
        }
        break;
      case NIGAM_PHILLIPS_1 :
        {
          // x, y, z de la pyramide infinie
          Real_wp xi = 0.5*(x/(1.0-z) + 1.0);
          Real_wp yi = 0.5*(y/(1.0-z) + 1.0);
          Real_wp zi = z/(1.0-z);
			
          Real_wp q1, q2, q3, q4;
          
          // DF^-1
          Matrix3_3 dfjm1;
          dfjm1(0, 0) = 0.5/(1.0-z);
          dfjm1(1, 1) = 0.5/(1.0-z);
          dfjm1(2, 2) = 1.0/square(1.0-z);
          dfjm1(0, 2) = 0.5*x/square(1.0-z);
          dfjm1(1, 2) = 0.5*y/square(1.0-z);
          
          VectR3 phi(nb_dof_loc);
          int node = 0;
          // edges
          for (int i = 0; i <= r-1; i++)
            {
              phi(node).Init(-1.0/pow(1+zi,r+1-i)*zi*(yi-1),
                             -1.0/pow(1+zi,r+1-i)*zi*(xi-1),
                             1.0/pow(1+zi,r+1)*(xi-1)*(yi-1)*pow(1+zi,i));
              phi(node+1).Init(-1.0/pow(1+zi,r+1-i)*zi*(yi-1),-1.0/pow(1+zi,r+1-i)*zi*xi,
                               1.0/pow(1+zi,r+1)*xi*(yi-1)*pow(1+zi,i));
              phi(node+2).Init(-1.0/pow(1+zi,r+1-i)*zi*yi,-1.0/pow(1+zi,r+1-i)*zi*xi,
                               1.0/pow(1+zi,r+1)*xi*yi*pow(1+zi,i));
              phi(node+3).Init(-1.0/pow(1+zi,r+1-i)*zi*yi,-1.0/pow(1+zi,r+1-i)*zi*(xi-1),
                               1.0/pow(1+zi,r+1)*yi*(xi-1)*pow(1+zi,i));
              node += 4;			  
            }
          
          for (int i = 0; i <= r-1; i++)
            {
              q1 = 1.0/pow(1+zi,r+1)*pow(xi,i);
              q2 = 1.0/pow(1+zi,r+1)*pow(yi,i);
              
              phi(node).Init((1-yi)*q1,0,0);
              phi(node+1).Init(0,xi*q2,0);
              phi(node+2).Init(yi*q1,0,0);
              phi(node+3).Init(0,(1-xi)*q2,0);
              node += 4;
            }
          
          // faces
          for (int i = 0; i <= r-2; i++)
            for (int j = 0; j <= r-2; j++)
              if (i+j <= r-2)
                {
                  phi(node).Init(zi*(1-yi)/pow(1+zi,r+1)*pow(xi,i)*pow(1+zi,j),0,0);
                  phi(node+1).Init(0,zi*xi/pow(1+zi,r+1)*pow(yi,i)*pow(1+zi,j),0);
                  phi(node+2).Init(zi*yi/pow(1+zi,r+1)*pow(xi,i)*pow(1+zi,j),0,0);
                  phi(node+3).Init(0,zi*(1-xi)/pow(1+zi,r+1)*pow(yi,i)*pow(1+zi,j),0);
                  node += 4;
                }
          
          for (int i = 0; i <= r-3; i++)
            for (int j = 0; j <= r-3; j++)
              if (i+j <= r-3)
                {
                  phi(node).Init(0,0,xi*(1-xi)*(1-yi)/pow(1+zi,r+1)*pow(xi,i)*pow(zi,j));
                  phi(node+1).Init(0,0,yi*(1-yi)*xi/pow(1+zi,r+1)*pow(yi,i)*pow(zi,j));
                  phi(node+2).Init(0,0,xi*(1-xi)*yi/pow(1+zi,r+1)*pow(xi,i)*pow(zi,j));
                  phi(node+3).Init(0,0,yi*(1-yi)*(1-xi)/pow(1+zi,r+1)*pow(yi,i)*pow(zi,j));
                  node += 4;
                }
          
          for (int i = 0; i <= r-2; i++)
            {
              q1 = (1-xi)*(1-yi)*pow(xi,i)*pow(1+zi,r-i-2)/pow(1+zi,r+1);
              q2 = xi*(1-yi)*pow(yi,i)*pow(1+zi,r-i-2)/pow(1+zi,r+1);
              q3 = xi*yi*pow(xi,i)*pow(1+zi,r-i-2)/pow(1+zi,r+1);
              q4 = (1-xi)*yi*pow(yi,i)*pow(1+zi,r-i-2)/pow(1+zi,r+1);
              
              phi(node).Init(zi*q1,0,-xi*q1);
              phi(node+1).Init(0,zi*q2,-yi*q2);
              phi(node+2).Init(zi*q3,0,-xi*q3);
              phi(node+3).Init(0,zi*q4,-yi*q4);
              node += 4;
            }
          
          for (int i = 0; i <= r-1; i++)
            for (int j = 0; j <= r-2; j++)
              {
                phi(node).Init(1.0/pow(1+zi,r+1)*yi*(1-yi)*pow(xi,i)*pow(yi,j),0,0);
                phi(node+1).Init(0,1.0/pow(1+zi,r+1)*xi*(1-xi)*pow(xi,j)*pow(yi,i),0);
                node += 2;
              }
          
          // interior
          for (int j = 0; j <= r-2; j++)
            for (int k = 0; k <= r-2; k++)
              {
                for (int i = 0; i <= r-1; i++)
                  {
                    phi(node).Init(yi*(1-yi)*zi/pow(1+zi,r+1)*pow(xi,i)*pow(yi,j)*pow(zi,k),0,0);
                    phi(node+1).Init(0,xi*(1-xi)*zi/pow(1+zi,r+1)*pow(xi,k)*pow(yi,i)*pow(zi,j),0);
                    node += 2;
                  }	
                
                for (int i = 0; i <= r-2; i++)
                  {
                    phi(node).Init(0,0,xi*(1-xi)*yi*(1-yi)/pow(1+zi,r+1)
                                   *pow(xi,i)*pow(yi,j)*pow(zi,k));
                    node++;
                  }
                
                q1 = -pow(zi,r)/pow(1+zi,r+1)*(-1+yi)*pow(yi,1+k)
                  *(pow(xi,j)-2*pow(xi,1+j)+pow(xi,j)*j-pow(xi,1+j)*j);
                q2 = pow(zi,r)/pow(1+zi,r+1)*pow(xi,1+j)*(-1+xi)
                  *(-pow(yi,k)+2*pow(yi,1+k)-pow(yi,k)*k+pow(yi,1+k)*k);
                q3 = pow(zi,r-1)/pow(1+zi,r+1)*xi*(1-xi)*yi*(1-yi)*pow(xi,j)*pow(yi,k);
		
                phi(node).Init(q1,q2,-q3);
                node++; 				   
              }
          
          /*
          for (int k = 0; k <= r-2; k++)
            for (int j = 0; j <= k; j++)
              {
                for (int i = 0; i <= k+1; i++)
                  {
                    phi(node).Init(yi*(1-yi)*zi/pow(1+zi,r+1)*pow(xi,i)*pow(yi,j)*pow(1+zi,k),0,0);
                    phi(node+1).Init(0,xi*(1-xi)*zi/pow(1+zi,r+1)*pow(xi,j)
                    *pow(yi,i)*pow(1+zi,k),0);
                    node += 2;
                  }
                  
                for (int i = 0; i <= k; i++)
                  {
                    phi(node).Init(0,0,xi*(1-xi)*yi*(1-yi)/pow(1+zi,r+1)
                    *pow(xi,i)*pow(yi,j)*pow(1+zi,k));
                    node++;
                  }
              }          
          
          for (int k = 0; k <= r-2; k++)
            for (int j = 0; j <= r-2; j++)
              {
                q1 = -pow(zi,r)/pow(1+zi,r+1)*(-1+yi)*pow(yi,1+k)
                *(pow(xi,j)-2*pow(xi,1+j)+pow(xi,j)*j-pow(xi,1+j)*j);
                
                q2 = pow(zi,r)/pow(1+zi,r+1)*pow(xi,1+j)*(-1+xi)
                *(-pow(yi,k)+2*pow(yi,1+k)-pow(yi,k)*k+pow(yi,1+k)*k);
                
                q3 = pow(zi,r-1)/pow(1+zi,r+1)*xi*(1-xi)*yi*(1-yi)*pow(xi,j)*pow(yi,k);
                
                phi(node).Init(q1,q2,-q3);
                node++; 				   
              }
          */
          
          // on applique DFi*-1
          for (int k = 0; k < nb_dof_loc; k++)
            MltTrans(dfjm1, phi(k), res(k));
        }
        break;
      case NIGAM_PHILLIPS_2 :
        {
          // x, y, z de la pyramide infinie
          Real_wp xi = 0.5*(x/(1.0-z) + 1.0);
          Real_wp yi = 0.5*(y/(1.0-z) + 1.0);
          Real_wp zi = z/(1.0-z);
          
          // DF^-1
          Matrix3_3 dfjm1;
          dfjm1(0, 0) = 0.5/(1.0-z);
          dfjm1(1, 1) = 0.5/(1.0-z);
          dfjm1(2, 2) = 1.0/square(1.0-z);
          dfjm1(0, 2) = 0.5*x/square(1.0-z);
          dfjm1(1, 2) = 0.5*y/square(1.0-z);
          
          VectR3 phi(nb_dof_loc);
          int node = 0;
          for (int c = 0; c <= r+1; c++)
            {
              for (int a = 0; a <= c-2; a++)
                for (int b = 0; b <= c-1; b++)
                  {
                    phi(node).Init(pow(xi, a)*pow(yi, b)/pow(1.0+zi, c), 0, 0);
                    phi(node+1).Init(0, pow(xi, b)*pow(yi, a)/pow(1.0+zi, c), 0);
                    node += 2;
                  }
            }
          
          for (int c = 0; c <= r; c++)
            for (int a = 0; a <= c; a++)
              for (int b = 0; b <= c; b++)
                {
                  if ((a != 0) || (b != 0) || (c != 0))
                    {
                      Real_wp vx = 0;
                      if (a != 0)
                        vx = Real_wp(a)*pow(xi, a-1)*pow(yi, b)/pow(1.0+zi, c);
                      
                      Real_wp vy = 0;
                      if (b != 0)
                        vy = Real_wp(b)*pow(xi, a)*pow(yi, b-1)/pow(1.0+zi, c);
                      
                      Real_wp vz = 0;
                      if (c != 0)
                        vz = -Real_wp(c)*pow(xi, a)*pow(yi, b)/pow(1.0+zi, c+1);
                      
                      phi(node).Init(vx, vy, vz);
                      node++;
                    }
                }
          
          // on applique DFi*-1
          for (int k = 0; k < nb_dof_loc; k++)
            MltTrans(dfjm1, phi(k), res(k));
        }
        break;
      case ZAGLMAYR :
        {
          // on calcule xu, yu, zu de la pyramide unite
          Real_wp xu = 0.5*(1.0 + x - z);
          Real_wp yu = 0.5*(1.0 + y - z);
          Real_wp zu = z;
          
          // on calcule DF-1
          Matrix3_3 dfjm1;
          dfjm1(0, 0) = 0.5; dfjm1(1, 1) = 0.5; dfjm1(2, 2) = 1.0;
          dfjm1(0, 2) = -0.5; dfjm1(1, 2) = -0.5;
          
          int node = 0;
          VectR3 phi(nb_dof_loc);
          for (int k = 0; k <= r-1; k++)
            for (int i = 0; i <= k; i++)
              for (int j = 0; j <= k; j++)
                {
                  Real_wp coef = pow(xu, i)*pow(yu, j)*pow(1.0-zu, k-i-j);
                  phi(node).Init(coef, 0, 0);
                  phi(node+1).Init(0, coef, 0);
                  phi(node+2).Init(0, 0, coef);
                  node += 3;
                }
          
          for (int k = 1; k <= r; k++)
            for (int i = 1; i <= k; i++)
              for (int j = 1; j <= k; j++)
                {
                  DISP(node);
                  phi(node++).Init(pow(xu, i-1)*pow(yu, k)*pow(1.0-zu, -i),
                                   pow(xu, k)*pow(yu, j-1)*pow(1.0-zu, -j),
                                   pow(xu, k)*pow(yu, k)*pow(1.0-zu, -k-1));
                }
          
          for (int k = 0; k < r; k++)
            for (int i = 0; i <= k; i++)
              for (int j = 0; j < k; j++)
                {
                  phi(node++).Init(pow(xu, i)*pow(yu, k+1)*pow(1.0-zu, -i-1),
                                   0, pow(xu, j+1)*pow(yu, k+1)*pow(1.0-zu, -j-2));
                  
                  phi(node++).Init(0, pow(xu, k+1)*pow(yu, i)*pow(1.0-zu, -i-1),
                                   pow(xu, k+1)*pow(yu, j+1)*pow(1.0-zu, -j-2));
                }
          
          // on applique DFi*-1
          for (int k = 0; k < nb_dof_loc; k++)
            MltTrans(dfjm1, phi(k), res(k));

        }
        break;
      }
  }
  
  
  //! computes \f$ curl(\varphi_{node})(point\_loc) \f$
  void PyramidHcurlOther::ComputeCurlPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    FillZero(res);

    int r = order;
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    if (type_basis)
      {
        // on calcule xu, yu, zu de la pyramide unite
        Real_wp xu = 0.5*(1.0 + x - z);
        Real_wp yu = 0.5*(1.0 + y - z);
        Real_wp zu = z;
        
        // on calcule 1/Ji DFi
        Matrix3_3 invJi_DFi;
        invJi_DFi(0, 0) = 0.5; invJi_DFi(1, 1) = 0.5; invJi_DFi(2, 2) = 0.25;
        invJi_DFi(0, 2) = 0.25; invJi_DFi(1, 2) = 0.25;
          
        Real_wp N1 = xu;
        Real_wp N2 = yu;
        Real_wp N3 = 1.0-zu-xu;
        Real_wp N4 = 1.0-zu-yu;
        Real_wp N5 = zu;
        Real_wp invOmz = 1.0/(1.0 - N5);
        Real_wp invOmz2 = invOmz*invOmz;
        
        R3 grad_N1(1.0, 0, 0);
        R3 grad_N2(0, 1.0, 0);
        R3 grad_N3(-1.0, 0, -1.0);
        R3 grad_N4(0, -1.0, -1.0);
        R3 grad_N5(0, 0, 1.0);
        
        VectR3 curl_phi(nb_dof_loc);
        R3 v1, v2, v3, v4, curl_v1, curl_v2, curl_v3, curl_v4;
        R3 v5, v6, v7, v8, curl_v5, curl_v6, curl_v7, curl_v8;
        v5 = N3*N4*grad_N5*invOmz - N4*N5*grad_N3*invOmz
          - N3*N5*grad_N4*invOmz - N3*N4*N5*grad_N5*invOmz2;
        
        v6 = N1*N4*grad_N5*invOmz - N1*N5*grad_N4*invOmz
          - N4*N5*grad_N1*invOmz - N1*N4*N5*grad_N5*invOmz2;
        
        v7 = N1*N2*grad_N5*invOmz - N2*N5*grad_N1*invOmz
          - N1*N5*grad_N2*invOmz - N1*N2*N5*grad_N5*invOmz2;
        
        v8 = N2*N3*grad_N5*invOmz - N3*N5*grad_N2*invOmz
          - N2*N5*grad_N3*invOmz - N2*N3*N5*grad_N5*invOmz2;
        
        v4 = N2*N3*grad_N4*invOmz - N3*N4*grad_N2*invOmz; Mlt(-1.0, v4);
        v1 = N4*N3*grad_N1*invOmz - N1*N4*grad_N3*invOmz;
        v2 = N4*N1*grad_N2*invOmz - N1*N2*grad_N4*invOmz; 
        v3 = N2*N1*grad_N3*invOmz - N3*N2*grad_N1*invOmz; Mlt(-1.0, v3);
        
        curl_v5 = 2.0*invOmz*(-N3*grad_N1 + N4*grad_N2);
        curl_v6 = 2.0*invOmz*(-N4*grad_N2 - N1*grad_N1);
        curl_v7 = 2.0*invOmz*(N1*grad_N1 - N2*grad_N2);
        curl_v8 = 2.0*invOmz*(N2*grad_N2 + N3*grad_N1);
        
        curl_v1 = 0.5*curl_v6; curl_v1(1) += -1.0; curl_v1(2) += 1.0;
        curl_v2 = 0.5*curl_v7; curl_v2(2) += 1.0;
        curl_v3 = 0.5*curl_v8; curl_v3(0) -= 1.0; curl_v3(2) += 1.0; Mlt(-1.0, curl_v3);
        curl_v4 = 0.5*curl_v5; curl_v4(0) -= 1.0; curl_v4(1) += -1.0;
        curl_v4(2) += 1.0; Mlt(-1.0, curl_v4);
        
        if (r > 1)
          {
            // dofs on edges
            Real_wp alpha = (3.0*N3 - 1.0)*(2.0*N4 + 2.0*N5 - 1.0);
            R3 grad_alpha = 3.0*grad_N3*(2.0*N4 + 2.0*N5 - 1.0)
              + (3.0*N3-1.0)*2.0*(grad_N4+grad_N5);            
            
            TimesProd(grad_alpha, v1, curl_phi(0));
            Add(alpha, curl_v1, curl_phi(0));
            
            alpha = (3.0*N1 - 1.0)*(2.0*N4 + 2.0*N5 - 1.0);
            grad_alpha = 3.0*grad_N1*(2.0*N4 + 2.0*N5 - 1.0) + (3.0*N1-1.0)*2.0*(grad_N4+grad_N5);
            TimesProd(grad_alpha, v1, curl_phi(1));
            Add(alpha, curl_v1, curl_phi(1));
            
            alpha = (3.0*N4 - 1.0)*(2.0*N1 + 2.0*N5 - 1.0);
            grad_alpha = 3.0*grad_N4*(2.0*N1 + 2.0*N5 - 1.0)
              + (3.0*N4 - 1.0)*2.0*(grad_N1+grad_N5);
            
            TimesProd(grad_alpha, v2, curl_phi(2));
            Add(alpha, curl_v2, curl_phi(2));
            
            alpha = (3.0*N2 - 1.0)*(2.0*N1 + 2.0*N5 - 1.0);
            grad_alpha = 3.0*grad_N2*(2.0*N1 + 2.0*N5 - 1.0)
              + (3.0*N2 - 1.0)*2.0*(grad_N1+grad_N5);
            
            TimesProd(grad_alpha, v2, curl_phi(3));
            Add(alpha, curl_v2, curl_phi(3));
            
            alpha = (3.0*N3 - 1.0)*(2.0*N2 + 2.0*N5 - 1.0);
            grad_alpha = 3.0*grad_N3*(2.0*N2 + 2.0*N5 - 1.0)
              + (3.0*N3 - 1.0)*2.0*(grad_N2+grad_N5);
            
            TimesProd(grad_alpha, v3, curl_phi(4));
            Add(alpha, curl_v3, curl_phi(4));
            
            alpha = (3.0*N1 - 1.0)*(2.0*N2 + 2.0*N5 - 1.0);
            grad_alpha = 3.0*grad_N1*(2.0*N2 + 2.0*N5 - 1.0)
              + (3.0*N1 - 1.0)*2.0*(grad_N2+grad_N5);
            
            TimesProd(grad_alpha, v3, curl_phi(5));
            Add(alpha, curl_v3, curl_phi(5));
            
            alpha = (3.0*N4 - 1.0)*(2.0*N3 + 2.0*N5 - 1.0);
            grad_alpha = 3.0*grad_N4*(2.0*N3 + 2.0*N5 - 1.0)
              + (3.0*N4 - 1.0)*2.0*(grad_N3+grad_N5);
            
            TimesProd(grad_alpha, v4, curl_phi(6));
            Add(alpha, curl_v4, curl_phi(6));
            
            alpha = (3.0*N2 - 1.0)*(2.0*N3 + 2.0*N5 - 1.0);
            grad_alpha = 3.0*grad_N2*(2.0*N3 + 2.0*N5 - 1.0)
              + (3.0*N2 - 1.0)*2.0*(grad_N3+grad_N5);
            
            TimesProd(grad_alpha, v4, curl_phi(7));
            Add(alpha, curl_v4, curl_phi(7));

            alpha = (3.0*N5 - 1.0);
            grad_alpha = 3.0*grad_N5;
            TimesProd(grad_alpha, v5, curl_phi(9));
            Add(alpha, curl_v5, curl_phi(9));
            
            alpha = (3.0*N3 - 3.0*N2 - 1.0);
            grad_alpha = 3.0*grad_N3 - 3.0*grad_N2;
            TimesProd(grad_alpha, v5, curl_phi(8));
            Add(alpha, curl_v5, curl_phi(8));
            
            alpha = (3.0*N5 - 1.0);
            grad_alpha = 3.0*grad_N5;
            TimesProd(grad_alpha, v6, curl_phi(11));
            Add(alpha, curl_v6, curl_phi(11));
            
            alpha = (3.0*N4 - 3.0*N3 - 1.0);
            grad_alpha = 3.0*grad_N4 - 3.0*grad_N3;
            TimesProd(grad_alpha, v6, curl_phi(10));
            Add(alpha, curl_v6, curl_phi(10));
            
            alpha = (3.0*N5 - 1.0);
            grad_alpha = 3.0*grad_N5;
            TimesProd(grad_alpha, v7, curl_phi(13));
            Add(alpha, curl_v7, curl_phi(13));
            
            alpha = (3.0*N1 - 3.0*N4 - 1.0);
            grad_alpha = 3.0*grad_N1 - 3.0*grad_N4;
            TimesProd(grad_alpha, v7, curl_phi(12));
            Add(alpha, curl_v7, curl_phi(12));
            
            alpha = 3.0*N5-1.0;
            grad_alpha = 3.0*grad_N5;
            TimesProd(grad_alpha, v8, curl_phi(15));
            Add(alpha, curl_v8, curl_phi(15));
            
            alpha = (3.0*N2 - 3.0*N1 - 1.0);
            grad_alpha = 3.0*grad_N2 - 3.0*grad_N1;
            TimesProd(grad_alpha, v8, curl_phi(14));
            Add(alpha, curl_v8, curl_phi(14));
            
            // dofs on quadrilateral face
            alpha = N2*(3.0*N1 - 1.0);
            grad_alpha = grad_N2*(3.0*N1 - 1.0) + 3.0*N2*grad_N1;
            TimesProd(grad_alpha, v1, curl_phi(16));
            Add(alpha, curl_v1, curl_phi(16));
            
            alpha = N2*(3.0*N3 - 1.0);
            grad_alpha = grad_N2*(3.0*N3 - 1.0) + 3.0*N2*grad_N3;
            TimesProd(grad_alpha, v1, curl_phi(17));
            Add(alpha, curl_v1, curl_phi(17));
            
            alpha = N1*(3.0*N4 - 1.0);
            grad_alpha = grad_N1*(3.0*N4 - 1.0) + 3.0*N1*grad_N4;
            TimesProd(grad_alpha, v4, curl_phi(18));
            Add(alpha, curl_v4, curl_phi(18));
            
            alpha = N1*(3.0*N2 - 1.0);
            grad_alpha = grad_N1*(3.0*N2 - 1.0) + 3.0*N1*grad_N2;
            TimesProd(grad_alpha, v4, curl_phi(19));
            Add(alpha, curl_v4, curl_phi(19));
            
            // dofs on triangular faces
            alpha = 3.0*N5*(1.0-3.0*N2);
            grad_alpha = 3.0*grad_N5*(1.0-3.0*N2) - 9.0*N5*grad_N2;
            TimesProd(grad_alpha, v1, curl_phi(20));
            Add(alpha, curl_v1, curl_phi(20));
            
            alpha = 3.0*N1;
            grad_alpha = 3.0*grad_N1;
            TimesProd(grad_alpha, v5, curl_phi(21));
            Add(alpha, curl_v5, curl_phi(21));

            alpha = 3.0*N5*(1.0-3.0*N3);
            grad_alpha = 3.0*grad_N5*(1.0-3.0*N3) - 9.0*N5*grad_N3;
            TimesProd(grad_alpha, v2, curl_phi(22));
            Add(alpha, curl_v2, curl_phi(22));
            
            alpha = 3.0*N2;
            grad_alpha = 3.0*grad_N2;
            TimesProd(grad_alpha, v6, curl_phi(23));
            Add(alpha, curl_v6, curl_phi(23));
            
            alpha = 3.0*N5*(1.0-3.0*N4);
            grad_alpha = 3.0*grad_N5*(1.0-3.0*N4) - 9.0*N5*grad_N4;
            TimesProd(grad_alpha, v3, curl_phi(24));
            Add(alpha, curl_v3, curl_phi(24));
            
            alpha = 3.0*N3;
            grad_alpha = 3.0*grad_N3;
            TimesProd(grad_alpha, v7, curl_phi(25));
            Add(alpha, curl_v7, curl_phi(25)); curl_phi(25) += curl_phi(24);

            alpha = 3.0*N5*(1.0-3.0*N1);
            grad_alpha = 3.0*grad_N5*(1.0-3.0*N1) - 9.0*N5*grad_N1;
            TimesProd(grad_alpha, v4, curl_phi(26));
            Add(alpha, curl_v4, curl_phi(26));
            
            alpha = 3.0*N4;
            grad_alpha = 3.0*grad_N4;
            TimesProd(grad_alpha, v8, curl_phi(27));
            Add(alpha, curl_v8, curl_phi(27)); curl_phi(27) += curl_phi(26);
            
            // dofs inside
            alpha = 9.0*N2*N5;
            grad_alpha = 9.0*grad_N2*N5 + 9.0*N2*grad_N5;
            TimesProd(grad_alpha, v1, curl_phi(28));
            Add(alpha, curl_v1, curl_phi(28));

            alpha = 9.0*N1*N5;
            grad_alpha = 9.0*grad_N1*N5 + 9.0*N1*grad_N5;
            TimesProd(grad_alpha, v4, curl_phi(29));
            Add(alpha, curl_v4, curl_phi(29));
          }
        else
          {
            curl_phi(0) = curl_v1;
            curl_phi(1) = curl_v2;
            curl_phi(2) = curl_v3;
            curl_phi(3) = curl_v4;
              
            curl_phi(4) = curl_v5;
            curl_phi(5) = curl_v6;              
            curl_phi(6) = curl_v7;
            curl_phi(7) = curl_v8;
          }
          
        // on applique 1/Ji DFi
        for (int k = 0; k < nb_dof_loc; k++)
          Mlt(invJi_DFi, curl_phi(k), res(k));
        
        return;
      }
    
    if (type_basis != ZGAINSKI)
        return;
    
    Real_wp omz, invOmz(0);
    omz = 1-z;
	  
    // parameters for the pyramid
    Real_wp beta1 = 0.5*(1.0-x-z);
    Real_wp beta2 = 0.5*(1.0-y-z);
    Real_wp beta3 = 0.5*(1.0+x-z);
    Real_wp beta4 = 0.5*(1.0+y-z);
        
    // coordinates on the cube
    //Real_wp a(0), b(0), c;
    //da_dx(0), da_dz(0), db_dy(0), db_dz(0);
    Real_wp lambda1(0), lambda2(0), lambda3(0), lambda4(0);
    if (abs(omz) > epsilon_machine)
      {
        invOmz = 1.0/(1.0-z);
        //a = x*invOmz;
        //b = y*invOmz;
        //da_dx = invOmz;
        //db_dy = invOmz;
        //da_dz = a*invOmz;
        //db_dz = b*invOmz;
        
        lambda1 = beta1*beta2*invOmz;
        lambda2 = beta2*beta3*invOmz;
        lambda3 = beta3*beta4*invOmz;
        lambda4 = beta4*beta1*invOmz;
      }      
    
    //c = 2.0*z - 1.0;
    
    // elementary functions for vertical edges
    R3 grad_L1(-0.5*beta2, -0.5*beta1, -0.5*(beta2+beta1));
    grad_L1 *= invOmz; grad_L1(2) += beta1*beta2*invOmz*invOmz;
    
    R3 grad_L2(0.5*beta2, -0.5*beta3, -0.5*(beta3+beta2));
    grad_L2 *= invOmz; grad_L2(2) += beta2*beta3*invOmz*invOmz;
    
    R3 grad_L3(0.5*beta4, 0.5*beta3, -0.5*(beta4+beta3));
    grad_L3 *= invOmz; grad_L3(2) += beta3*beta4*invOmz*invOmz;
    
    R3 grad_L4(-0.5*beta4, 0.5*beta1, -0.5*(beta1+beta4));
    grad_L4 *= invOmz; grad_L4(2) += beta4*beta1*invOmz*invOmz;
    
    // curl of elementary functions
    res(0).Init(0.5*(-grad_L1(1)+grad_L2(1)), 0.5*(grad_L1(2)+grad_L2(2)+grad_L1(0)-grad_L2(0)),
                -0.5*(grad_L1(1)+grad_L2(1)));
    
    res(r).Init(0.5*(-grad_L2(1)+grad_L3(1)-grad_L2(2)-grad_L3(2)), 0.5*(grad_L2(0)-grad_L3(0)),
                0.5*(grad_L2(0)+grad_L3(0)));
    
    res(2*r).Init(0.5*(grad_L3(1)-grad_L4(1)), 0.5*(grad_L3(2)+grad_L4(2)-grad_L3(0)+grad_L4(0)),
                  -0.5*(grad_L3(1)+grad_L4(1)));
    
    res(3*r).Init(0.5*(grad_L4(1)-grad_L1(1)-grad_L4(2)-grad_L1(2)), 0.5*(-grad_L4(0)+grad_L1(0)),
                  0.5*(grad_L4(0)+grad_L1(0)));
    
    res(4*r)(0) = 2.0*grad_L1(1); res(4*r)(1) = -2.0*grad_L1(0);
    res(5*r)(0) = 2.0*grad_L2(1); res(5*r)(1) = -2.0*grad_L2(0);
    res(6*r)(0) = 2.0*grad_L3(1); res(6*r)(1) = -2.0*grad_L3(0);
    res(7*r)(0) = 2.0*grad_L4(1); res(7*r)(1) = -2.0*grad_L4(0);
    
    if (r > 1)
      {
        res(1).Init(0.5*(-grad_L1(1)-grad_L2(1)), 0.5*(grad_L1(2)-grad_L2(2)+grad_L1(0)
                                                       +grad_L2(0)), -0.5*(grad_L1(1)-grad_L2(1)));
        
        res(r+1).Init(0.5*(-grad_L2(1)-grad_L3(1)-grad_L2(2)+grad_L3(2)),
                      0.5*(grad_L2(0)+grad_L3(0)), 0.5*(grad_L2(0)-grad_L3(0)));
        
        res(2*r+1).Init(-0.5*(grad_L3(1)+grad_L4(1)),
                        -0.5*(grad_L3(2)-grad_L4(2)-grad_L3(0)-grad_L4(0)),
                        0.5*(grad_L3(1)-grad_L4(1)));
        
        res(3*r+1).Init(-0.5*(grad_L4(1)+grad_L1(1)-grad_L4(2)+grad_L1(2)),
                        -0.5*(-grad_L4(0)-grad_L1(0)), -0.5*(grad_L4(0)-grad_L1(0)));
        
        R3 grad_beta1(-0.5, 0, -0.5);
        R3 grad_beta2(0, -0.5, -0.5);
        R3 grad_beta3(0.5, 0, -0.5);
        R3 grad_beta4(0, 0.5, -0.5);
        R3 grad_L5(0, 0, 1.0);
        
        R3 vec_u;
        vec_u = grad_L1*lambda2 + grad_L2*lambda1;
        TimesProd(vec_u, grad_beta4, res(16));
        
        vec_u = grad_L2*lambda3 + grad_L3*lambda2;
        TimesProd(vec_u, grad_beta1, res(17));
        
        vec_u = grad_L3*lambda4 + grad_L4*lambda3;
        TimesProd(vec_u, grad_beta2, res(18));
        
        vec_u = grad_L4*lambda1 + grad_L1*lambda4;
        TimesProd(vec_u, grad_beta3, res(19));
        
        vec_u = grad_L1*z + grad_L5*lambda1;
        TimesProd(vec_u, grad_beta3, res(20));
        
        vec_u = grad_L2*z + grad_L5*lambda2;
        TimesProd(vec_u, grad_beta1, res(21));
        
        vec_u = grad_L2*z + grad_L5*lambda2;
        TimesProd(vec_u, grad_beta4, res(22));
        
        vec_u = grad_L3*z + grad_L5*lambda3;
        TimesProd(vec_u, grad_beta2, res(23));
        
        vec_u = grad_L3*z + grad_L5*lambda3;
        TimesProd(vec_u, grad_beta1, res(25));
        
        vec_u = grad_L4*z + grad_L5*lambda4;
        TimesProd(vec_u, grad_beta3, res(24));
        
        vec_u = grad_L4*z + grad_L5*lambda4;
        TimesProd(vec_u, grad_beta2, res(27));
        
        vec_u = grad_L1*z + grad_L5*lambda1;
        TimesProd(vec_u, grad_beta4, res(26));

        vec_u = grad_L1*lambda3 + grad_L3*lambda1;
        TimesProd(vec_u, grad_L5, res(28));
      }
  }

  
  void TriangleHcurlOther::
  ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad, int rsurf, int type_surf)
  {
    TriangleReference<2>::ConstructFiniteElement(r, rgeom, r+1, type_quad);
    
    nb_dof_loc = order*(order+2);    
    this->Value_PhiVec.Reallocate(nb_dof_loc, this->PointsND().GetM());
    VectR2 phi;
    for (int i = 0; i < this->PointsND().GetM(); i++)
      {
	this->ComputeValuesPhiRef(this->PointsND(i), phi);
	for (int j = 0; j < nb_dof_loc; j++)
          this->Value_PhiVec(j, i) = phi(j);
      }
    
  }
  
  
  void TriangleHcurlOther::ComputeValuesPhiRef(const R2& pt, VectR2& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    int r = order;
    Real_wp x = pt(0), y = pt(1);
    phi(0).Init(1.0-y, x);
    phi(r).Init(-y, x);
    phi(2*r).Init(-y, x-1.0);
    
    if (r > 1)
      {
        phi(1).Init(1.0-2.0*x-y, -x);
        phi(r+1).Init(y, x);
        phi(2*r+1).Init(-y, 1.0-x-2.0*y);
        
        if (type_basis == PyramidHcurlOther::ZGAINSKI)
          {
            phi(3*r).Init((1.0-x-y)*y, 0);
            phi(3*r+1).Init(-x*y, -x*y);
          }
        else
          {
            phi(3*r).Init(3.0*(1.0-y)*y, 3.0*x*y);
            phi(3*r+1).Init(3.0*y*x, 3.0*(1.0-x)*x);
          }
      }
  }


  void TriangleHcurlOther
  ::ComputeCurlPhiRef(const R2& pt, VectReal_wp& curl_phi) const
  {
    cout << "Not implemented for this finite element" << endl;
    abort(); 
  }  
  

  void QuadrangleHcurlOther::
  ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad, int rsurf, int type_surf)
  {    
    QuadrangleReference<2>::ConstructFiniteElement(r, rgeom, r+1,
                                                   Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    nb_dof_loc = 2*order*(order+1);    
    this->Value_PhiVec.Reallocate(nb_dof_loc, this->PointsND().GetM());
    VectR2 phi;
    for (int i = 0; i < this->PointsND().GetM(); i++)
      {
	this->ComputeValuesPhiRef(this->PointsND(i), phi);
	for (int j = 0; j < nb_dof_loc; j++)
          this->Value_PhiVec(j, i) = phi(j);
      }

  }
  
  
  void QuadrangleHcurlOther::ComputeValuesPhiRef(const R2& pt, VectR2& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    int r = order;
    Real_wp x = pt(0), y = pt(1);
    phi(0).Init(1.0-y, 0);
    phi(r).Init(0, x);
    phi(2*r).Init(-y, 0);
    phi(3*r).Init(0, x-1.0);
    
    if (r > 1)
      {
        phi(1).Init((2.0*x-1.0)*(1.0-y), 0);
        phi(r+1).Init(0, x*(2.0*y-1.0));
        phi(2*r+1).Init(-(2.0*x-1.0)*y, 0);
        phi(3*r+1).Init(0, (x-1.0)*(2.0*y-1.0));
        
        if (type_basis == PyramidHcurlOther::ZGAINSKI)
          {
            phi(4*r).Init(0, x*(1.0-x)*(1.0-y)*(1.0-y));
            phi(4*r+1).Init(-y*(1.0-y)*x*x, 0);
            phi(4*r+2).Init(0, -x*(1.0-x)*y*y);
            phi(4*r+3).Init(y*(1.0-y)*(1.0-x)*(1.0-x), 0);
          }
        else
          {
            phi(4*r+1).Init((1.0-y)*y*(3.0*(1.0-x) - 1.0), 0);
            phi(4*r).Init((1.0-y)*y*(3.0*x - 1.0), 0);
            phi(4*r+2).Init(0, (1.0-x)*x*(3.0*(1.0-y) - 1.0));
            phi(4*r+3).Init(0, (1.0-x)*x*(3.0*y - 1.0));
          }
      }
  }
  

  void QuadrangleHcurlOther
  ::ComputeCurlPhiRef(const R2& pt, VectReal_wp& curl_phi) const
  {
    cout << "Not implemented for this finite element" << endl;
    abort(); 
  }  
  
  
  //! displays details of class PyramidHcurlFirstFamily
  ostream& operator <<(ostream& out, const PyramidHcurlOther& e)
  {
    out<<static_cast<const PyramidReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_HCURL_OTHER_CXX
#endif
