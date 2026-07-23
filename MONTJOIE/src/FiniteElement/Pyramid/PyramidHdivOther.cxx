#ifndef MONTJOIE_FILE_PYRAMID_HDIV_OTHER_CXX

#include "PyramidHdivOther.hxx"

namespace Montjoie
{
  // File with all pyramidal functions described in the litterature
  
  //! default constructor
  PyramidHdivOther::PyramidHdivOther() : PyramidReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    type_basis = NIGAM_PHILLIPS_2;
    //type_basis = OPTIMAL;
    //type_basis = GRAGLIA;
    //type_basis = GRADINARU;
  }
  
    
  //! how to number mesh
  void PyramidHdivOther::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return PyramidReference<3>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, 0);
    nmap.SetNbDofQuadrangle(order, nb_dof_quad);
    
    // triangular face
    nmap.SetNbDofTriangle(order, nb_dof_tri);
    
    // dofs inside the Pyramid
    nmap.SetNbDofPyramid(order, nb_dof_loc-nb_dof_boundaries);
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    ElementReference<Dimension2, 3>::FindHdivRotationTri(order, this->Points2D_tri(),
							 this->Weights2D_tri(),
							 ValuePhiTri2D, nmap);
    
    ElementReference<Dimension2, 3>::FindHdivRotationQuad(order, NumQuad2D,
							  ValuePhiQuad2D, nmap);
  }

    
  //! construction of finite element
  void PyramidHdivOther::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						int rsurf_tri, int rsurf_quad,
						int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (rquad == 0)
      rquad = r+1;
    
    type_quad = PyramidQuadrature::QUADRATURE_GAUSS;
    //type_quad = PyramidQuadrature::QUADRATURE_JACOBI1;
    PyramidReference<3>::ConstructFiniteElement(r, rgeom, r+1, type_quad, r, r);
    
    // computation of basis functions
    ConstructFunctions();
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    int nb_points_quadrature = this->PointsND().GetM();
    Value_PhiVec.Reallocate(nb_dof_loc, nb_points_quadrature);
    VectR3 phi;
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	this->ComputeValuesPhiRef(this->PointsND(i), phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  Value_PhiVec(j, i) = phi(j);
      }
    
    /*for (int num_loc = 0; num_loc < 5; num_loc++)
      {
	for (int j = 0; j < nb_dof_loc; j++)
	  {
	    bool dof_on_face = false;
	    for (int k = 0; k < this->GetNbQuadBoundary(num_loc); k++)
	      {
		int kglob = this->GetQuadNumber(num_loc, k);
		Real_wp val = DotProd(Value_PhiVec(j, kglob), normale(num_loc));
		if (abs(val) > 1e4*epsilon_machine)
		    dof_on_face = true;
	      }
		
	    if (dof_on_face)
	      cout << " Dof " << j << " associated with face " << num_loc << endl;
	  }
      }
    */
    
    // construction of Value_PhiVec, Curl_Phi and elementary matrices
    ConstructHdivElementaryMatrix();
    
    //this->mass_matrix.Write("MhJex_Other"+to_str(this->order)+".dat");
    //this->const_div_matrix.Write("KhJex_Other"+to_str(this->order)+".dat");
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace();
    
    ValuePhiQuad2D.Reallocate(nb_dof_quad, this->Points2D_quad().GetM());
    for (int i = 0; i < this->Points2D_quad().GetM(); i++)
      {
        Real_wp x = 2.0*this->Points2D_quad()(i)(0) - 1.0;
        Real_wp y = 2.0*this->Points2D_quad()(i)(1) - 1.0;
        R3 pt(x, y, 0.0); VectR3 phi;
        ComputeValuesPhiRef(pt, phi);
        for (int k = 0; k < nb_dof_quad; k++)
          ValuePhiQuad2D(k, i) = DotProd(phi(k), this->NormaleLoc(0));
      }    
    
    ValuePhiTri2D.Reallocate(nb_dof_tri, this->Points2D_tri().GetM());
    for (int i = 0; i < this->Points2D_tri().GetM(); i++)
      {
        Real_wp x = 2.0*this->Points2D_tri()(i)(0) + this->Points2D_tri()(i)(1) - 1.0;
        Real_wp y = this->Points2D_tri()(i)(1) - 1.0;
        Real_wp z = this->Points2D_tri()(i)(1);
        R3 pt(x, y, z); VectR3 phi;
        ComputeValuesPhiRef(pt, phi);
        for (int k = nb_dof_quad; k < nb_dof_quad+nb_dof_tri; k++)
          {
            ValuePhiTri2D(k-nb_dof_quad, i) = DotProd(phi(k), this->NormaleLoc(1));
          }
      }    
    
    /*
    VectR2 Pts2D(3);
    Pts2D(0).Init(0, 0);
    Pts2D(1).Init(1, 0);
    Pts2D(2).Init(0, 0.8);
    for (int i = 0; i < Pts2D.GetM(); i++)
      {
        Real_wp x = 2.0*Pts2D(i)(0) + Pts2D(i)(1) - 1.0;
        Real_wp y = Pts2D(i)(1) - 1.0;
        Real_wp z = Pts2D(i)(1);
        R3 pt(x, y, z); VectR3 phi;
        ComputeValuesPhiRef(pt, phi);
        for (int k = nb_dof_quad; k < nb_dof_quad+nb_dof_tri; k++)
          {
            ValuePhiTri2D(k-nb_dof_quad, i) = DotProd(phi(k), normale(1));
          }
      }
    
    cout << "face 1 " << endl;
    DISP(ValuePhiTri2D);

    for (int i = 0; i < Pts2D.GetM(); i++)
      {
        Real_wp x = 1.0 - Pts2D(i)(1);
        Real_wp y = 2.0*Pts2D(i)(0) + Pts2D(i)(1) - 1.0;
        Real_wp z = Pts2D(i)(1);
        R3 pt(x, y, z); VectR3 phi;
        ComputeValuesPhiRef(pt, phi);
        for (int k = nb_dof_quad+nb_dof_tri; k < nb_dof_quad+2*nb_dof_tri; k++)
          {
            ValuePhiTri2D(k-nb_dof_quad-nb_dof_tri, i) = DotProd(phi(k), normale(2));
          }
      }
    
    cout << "face 2 " << endl;
    DISP(ValuePhiTri2D);

    for (int i = 0; i < Pts2D.GetM(); i++)
      {
        Real_wp x = 2.0*Pts2D(i)(0) + Pts2D(i)(1) - 1.0;
        Real_wp y = 1.0 - Pts2D(i)(1);
        Real_wp z = Pts2D(i)(1);
        R3 pt(x, y, z); VectR3 phi;
        ComputeValuesPhiRef(pt, phi);
        for (int k = nb_dof_quad+2*nb_dof_tri; k < nb_dof_quad+3*nb_dof_tri; k++)
          {
            ValuePhiTri2D(k-nb_dof_quad-2*nb_dof_tri, i) = DotProd(phi(k), normale(3));
          }
      }

    cout << "face 3 " << endl;
    DISP(ValuePhiTri2D);
    
    for (int i = 0; i < Pts2D.GetM(); i++)
      {
        Real_wp x = Pts2D(i)(1) - 1.0;
        Real_wp y = 2.0*Pts2D(i)(0) + Pts2D(i)(1) - 1.0;
        Real_wp z = Pts2D(i)(1);
        R3 pt(x, y, z); VectR3 phi;
        ComputeValuesPhiRef(pt, phi);
        for (int k = nb_dof_quad+3*nb_dof_tri; k < nb_dof_quad+4*nb_dof_tri; k++)
          {
            ValuePhiTri2D(k-nb_dof_quad-3*nb_dof_tri, i) = DotProd(phi(k), normale(4));
          }
      }

    cout << "face 4 " << endl;
    DISP(ValuePhiTri2D);
    */
  }
  

  //! construction of basis functions
  void PyramidHdivOther::ConstructOrthogonalBasis(int r)
  {
    PyramidReference<3>::ConstructFiniteElement(r);
    
    if (type_basis == NIGAM_PHILLIPS_2)
      {
        order = r;
        nb_dof_loc = 0;
        for (int k = 0; k <= order+2; k++)
          for (int i = 0; i <= k-3; i++)
            for (int j = 0; j <= k-3; j++)
              nb_dof_loc++;
        
        for (int k = 0; k <= order+1; k++)
          for (int i = 0; i <= k-2; i++)
            for (int j = 0; j <= k-1; j++)
              nb_dof_loc += 2;    
      }
    else
      {
        order = r;
        nb_dof_loc = 2*(order+1)*(order-1)*order + order*order*order;
        nb_dof_loc += 2*order*(order+1);
      }
  }
  
  
  //! construction of basis functions
  void PyramidHdivOther::ConstructFunctions()
  {
    nb_dof_tri = order*(order+1)/2;
    nb_dof_quad = order*order;
    nb_dof_boundaries = 4*nb_dof_tri + nb_dof_quad;
    nb_dof_loc = 5;
    
    if ((type_basis == OPTIMAL) || (type_basis == OPTIMAL_HP))
      {
        nb_dof_quad = (order+1)*(order+1);
        nb_dof_loc = (order+1)*(2*order*order+7*order+2)/2;
      }
    else
      {
        if (type_basis == GRAGLIA)
          {
            if (order == 2)
              nb_dof_loc = 19;
          }
        else if (type_basis == NIGAM_PHILLIPS_1)
          {
            nb_dof_loc = 2.0*order + 3.0*order*order*order;
          }
        else if (type_basis == NIGAM_PHILLIPS_2)
          {
            nb_dof_loc = 0;
            for (int k = 0; k <= order+2; k++)
              for (int i = 0; i <= k-3; i++)
                for (int j = 0; j <= k-3; j++)
                  nb_dof_loc++;
            
            for (int k = 0; k <= order+1; k++)
              for (int i = 0; i <= k-2; i++)
                for (int j = 0; j <= k-1; j++)
                  nb_dof_loc += 2;
          }
      }
    
    VectR3 points_dof3d;    
    points_dof3d = this->PointsND();
    
    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    
    this->SetPointsDofND(points_dof3d);
  }
  
  
  //! compute \f$ \varphi_{node} (point\_loc) \f$
  void PyramidHdivOther::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();

    // on recupere x, y, z de la pyramide symetrique
    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);
    int r = order;
    res.Reallocate(nb_dof_loc);
    switch (type_basis)
      {
      case GRAGLIA :
        {
          // on calcule xu, yu, zu de la pyramide unite
	  // x = 2 xu + zu - 1
	  // y = 2 yu + zu - 1
	  // z = zu
          Real_wp xu = 0.5*(1.0 + x - z);
          Real_wp yu = 0.5*(1.0 + y - z);
          Real_wp zu = z;
          
          // on calcule 1/J DF
          Matrix3_3 DF;
          DF(0, 0) = 2.0; DF(1, 1) = 2.0;
	  DF(0, 2) = 1.0; DF(1, 2) = 1.0; DF(2, 2) = 1.0;
	  Real_wp invJacob = 1.0/Det(DF);
	  Mlt(invJacob, DF);
          
	  // calcul de phi
	  VectR3 phi(nb_dof_loc);
          
          if (order == 1)
            {
              phi(4).Init((-2.0*(1.0-xu-zu) - xu*zu)/(1.0-zu), -yu*zu/(1.0-zu), zu);
              phi(1).Init(-xu*zu/(1.0-zu), (-2.0*(1.0-yu-zu) -yu*zu)/(1.0-zu), zu);
              phi(2).Init(xu*(2.0-zu)/(1.0-zu), -yu*zu/(1.0-zu), zu);
              phi(3).Init(-xu*zu/(1.0-zu), yu*(2.0-zu)/(1.0-zu), zu);
              phi(0).Init(xu, yu, -(1.0-zu));
              
              Mlt(-1.0, phi(3));
              Mlt(-1.0, phi(4));
            }
          else
            {
              R3 v0, v1, v2, v3, v4;
              v4.Init((-2.0*(1.0-xu-zu) - xu*zu)/(1.0-zu), -yu*zu/(1.0-zu), zu);
              v1.Init(-xu*zu/(1.0-zu), (-2.0*(1.0-yu-zu) -yu*zu)/(1.0-zu), zu);
              v2.Init(xu*(2.0-zu)/(1.0-zu), -yu*zu/(1.0-zu), zu);
              v3.Init(-xu*zu/(1.0-zu), yu*(2.0-zu)/(1.0-zu), zu);
              v0.Init(xu, yu, -(1.0-zu));
              
              v3 *= -1.0; v4 *= -1.0;
              
              phi(0) = v0; phi(1) = v0; phi(2) = v0; phi(3) = v0;
              phi(0) *= (3.0*(1.0-xu) - 1.0)*(3.0*(1.0-yu) - 1.0) - 3.0*zu;
              phi(1) *= (3.0*yu - 1.0)*(3.0*(1.0-xu) - 1.0);
              phi(2) *= (3.0*xu - 1.0)*(3.0*(1.0-yu) - 1.0);
              phi(3) *= (3.0*xu - 1.0)*(3.0*yu - 1.0);
              
              phi(4) = v1; phi(5) = v1; phi(6) = v1;
              phi(4) *= 4.0*(1.0-xu-zu) - yu - 1.0;
              phi(5) *= 4.0*xu - yu - 1.0;
              phi(6) *= 4.0*zu - yu - 1.0; 
              
              phi(7) = v2; phi(8) = v2; phi(9) = v2;
              phi(7) *= 4.0*(1.0-yu-zu) - (1.0-xu-zu) - 1.0;
              phi(8) *= 4.0*yu - (1.0-xu-zu) - 1.0;
              phi(9) *= 4.0*zu - (1.0-xu-zu) - 1.0;
              
              phi(10) = v3; phi(11) = v3; phi(12) = v3;
              phi(10) *= 4.0*(1.0-xu-zu) - (1.0-yu-zu) - 1.0;
              phi(11) *= 4.0*xu - (1.0-yu-zu) - 1.0;
              phi(12) *= 4.0*zu - (1.0-yu-zu) - 1.0;
              
              phi(13) = v4; phi(14) = v4; phi(15) = v4;
              phi(13) *= 4.0*(1.0-yu-zu)-xu-1.0;
              phi(14) *= 4.0*yu - xu - 1.0;
              phi(15) *= 4.0*zu - xu - 1.0;              
              
              phi(16) = v4; phi(16) *= 3.0*xu;
              phi(17) = v1; phi(17) *= 3.0*yu;
              phi(18) = v0; phi(18) *= 3.0*zu;              
            }
          
          // on applique 1/J DF
          for (int k = 0; k < nb_dof_loc; k++)
            Mlt(DF, phi(k), res(k));
        }
        break;
      case GRADINARU :
        {          
          // on calcule xu, yu, zu de la pyramide unite
	  // x = 2 xu + zu - 1
	  // y = 2 yu + zu - 1
	  // z = zu
          Real_wp xu = 0.5*(1.0 + x - z);
          Real_wp yu = 0.5*(1.0 + y - z);
          Real_wp zu = z;
          
          // on calcule 1/J DF
          Matrix3_3 DF;
          DF(0, 0) = 2.0; DF(1, 1) = 2.0;
	  DF(0, 2) = 1.0; DF(1, 2) = 1.0; DF(2, 2) = 1.0;
	  Real_wp invJacob = 1.0/Det(DF);
	  Mlt(invJacob, DF);
          
	  //DISP(DF);
	  // calcul de phi (ordre 1 uniquement)
	  VectR3 phi(nb_dof_loc);
	  phi(1).Init(-xu*zu/(1.0-zu), yu - 2.0 + yu/(1.0-zu), zu);
	  phi(4).Init(xu - 2.0 + xu/(1.0-zu), -yu*zu/(1.0-zu), zu);
	  phi(2).Init(xu + xu/(1.0-zu), -yu*zu/(1.0-zu), zu);
	  phi(3).Init(-xu*zu/(1.0-zu), yu + yu/(1.0-zu), zu);
	  phi(0).Init(xu, yu, zu-1.0);
          
          Mlt(-1.0, phi(3));
          Mlt(-1.0, phi(4));
	  
          // on applique 1/J DF
          for (int k = 0; k < nb_dof_loc; k++)
            Mlt(DF, phi(k), res(k));
        }
        break;
      case NIGAM_PHILLIPS_1 :
        {
          // x, y, z de la pyramide infinie
          Real_wp xi = 0.5*(x/(1.0-z) + 1.0);
          Real_wp yi = 0.5*(y/(1.0-z) + 1.0);
          Real_wp zi = z/(1.0-z);
			
          // DF^-1
          Matrix3_3 dfjm1, DF;
          dfjm1(0, 0) = 0.5/(1.0-z);
          dfjm1(1, 1) = 0.5/(1.0-z);
          dfjm1(2, 2) = 1.0/square(1.0-z);
          dfjm1(0, 2) = 0.5*x/square(1.0-z);
          dfjm1(1, 2) = 0.5*y/square(1.0-z);
	  GetInverse(dfjm1, DF);
	  Real_wp invJacob = 1.0/Det(DF);
	  Mlt(invJacob, DF);

          VectReal_wp Px, Py, Pz(order+1), dPx, dPy;
	  EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*xi-1.0, Px, dPx);
	  EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*yi-1.0, Py, dPy);
          for (int i = 0; i <= order; i++)
            Pz(i) = pow(zi, i);
          
	  // calcul de phi
          VectR3 phi(nb_dof_loc);
          Real_wp den = 1.0/pow(1.0+zi, order+2), coef(0);
          int num = 0;
	  for (int i = 0; i < order; i++)
            for (int j = 0; j < order; j++)
              {
                coef = Px(i)*Py(j)*den;
                phi(num).Init(0, 0, coef);
                num++;
              }
          
          for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                coef = Px(i)*Pz(j)*den;
                phi(num).Init(0.0, 2.0*(1.0-yi)*coef, -zi*coef);
                num++;
              }
          
          for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                coef = Py(i)*Pz(j)*den;
                phi(num).Init(-2.0*xi*coef, 0, -zi*coef);
                num++;
              }
          
          for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                coef = Px(i)*Pz(j)*den;
                phi(num).Init(0.0, 2.0*yi*coef, zi*coef);
                num++;
              }
          
          for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                coef = Py(i)*Pz(j)*den;
                phi(num).Init(-2.0*(1.0-xi)*coef, 0, zi*coef);
                num++;
              }
          
          // part x(1-x) \khi_1
          for (int i = 0; i <= order-2; i++)
            for (int j = 0; j <= order-1; j++)
              for (int k = 0; k <= order-2; k++)
                {
                  coef = Px(i)*Py(j)*Pz(k)*den;
                  phi(num).Init(xi*(1.0-xi)*coef, 0, 0);
                  num++;
                }
          
          // part y(1-y) \khi_2
          for (int i = 0; i <= order-1; i++)
            for (int j = 0; j <= order-2; j++)
              for (int k = 0; k <= order-2; k++)
                {
                  coef = Px(i)*Py(j)*Pz(k)*den;
                  phi(num).Init(0, yi*(1.0-yi)*coef, 0);
                  num++;
                }
          
          // part z \khi_3
          for (int i = 0; i <= order-1; i++)
            for (int j = 0; j <= order-1; j++)
              for (int k = 0; k <= order-2; k++)
                {
                  coef = Px(i)*Py(j)*Pz(k)*den;
                  phi(num).Init(0, 0, zi*coef);
                  num++;
                }
          
          // part z^(k-1)/(1+z)^{k+2} (2t, 0, (1+z) t_x)
          Real_wp t, s, tx, sy;
          for (int i = 0; i <= order-2; i++)
            for (int j = 0; j <= order-1; j++)
              {
                // t et t_x
                t = xi*(1.0-xi)*Px(i)*Py(j);
                tx = (1.0 - 2.0*xi)*Px(i)*Py(j) + 2.0*xi*(1.0-xi)*dPx(i)*Py(j);
                coef = pow(zi, order-1)*den;
                
                phi(num).Init(2.0*t*coef, 0, (1.0+zi)*tx*coef);
                num++;
              }
          
          // part z^(k-1)/(1+z)^{k+2} (0, 2s, (1+z) s_y)
          for (int i = 0; i <= order-1; i++)
            for (int j = 0; j <= order-2; j++)
              {
                // s et s_y
                s = yi*(1.0-yi)*Px(i)*Py(j);
                sy = (1.0 - 2.0*yi)*Px(i)*Py(j) + 2.0*yi*(1.0-yi)*Px(i)*dPy(j);
                coef = pow(zi, order-1)*den;
                
                phi(num).Init(0.0, 2.0*s*coef, (1.0+zi)*sy*coef);
                num++;
              }
          
          // on applique 1/J DF
          for (int k = 0; k < nb_dof_loc; k++)
            Mlt(DF, phi(k), res(k));
        }
        break;
      case NIGAM_PHILLIPS_2 :
        {
          // x, y, z de la pyramide infinie
          Real_wp xi = 0.5*(x/(1.0-z) + 1.0);
          Real_wp yi = 0.5*(y/(1.0-z) + 1.0);
          Real_wp zi = z/(1.0-z);
			
          // DF^-1
          Matrix3_3 dfjm1, DF;
          dfjm1(0, 0) = 0.5/(1.0-z);
          dfjm1(1, 1) = 0.5/(1.0-z);
          dfjm1(2, 2) = 1.0/square(1.0-z);
          dfjm1(0, 2) = 0.5*x/square(1.0-z);
          dfjm1(1, 2) = 0.5*y/square(1.0-z);
	  GetInverse(dfjm1, DF);
	  Real_wp invJacob = 1.0/Det(DF);
	  Mlt(invJacob, DF);
          
	  // calcul de phi sur la pyramide infinie
          VectR3 phi(nb_dof_loc);
          
          VectReal_wp Px, Py, Pz(order+4), dPx, dPy;
	  EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*xi-1.0, Px, dPx);
	  EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*yi-1.0, Py, dPy);
          for (int i = 0; i <= order+3; i++)
            Pz(i) = pow(1.0 + zi, i);

          int num = 0;
          Real_wp coef;
          // dofs associated with the quadrilateral base
          Real_wp den = 1.0/Pz(order+2);
          for (int i = 0; i <= order-1; i++)
            for (int j = 0; j <= order-1; j++)
              {
                coef = Px(i)*Py(j)*den;
                phi(num).Init(0.0, 0.0, coef);
                num++;
              }
          
          // first triangular face
	  for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                int k = order+1-j;
                phi(num).Init(0.0, Real_wp(k)*(1.0-yi)*Px(i)*Pz(j)*den,
                              -Px(i)*(Pz(j+1)-1.0)*den);
                num++;
              }
          
          // second triangular face
          for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                int k = order+1-j;
                phi(num).Init(-Real_wp(k)*xi*Py(i)*Pz(j)*den,
                              0.0, -Py(i)*(Pz(j+1)-1.0)*den);
                num++;
              }
          
          // third triangular face
          for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                int k = order+1-j;
                phi(num).Init(0.0, Real_wp(k)*yi*Px(i)*Pz(j)*den,
                              Px(i)*(Pz(j+1)-1.0)*den);
                num++;
              }
          
          // last triangular face
          for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                int k = order+1-j;
                phi(num).Init(-Real_wp(k)*(1.0-xi)*Py(i)*Pz(j)*den,
                              0.0, Py(i)*(Pz(j+1)-1.0)*den);
                num++;
              }
          
          // interior dofs          
          for (int k = 0; k <= order+1; k++)
            for (int i = 0; i <= k-3; i++)
              for (int j = 0; j <= k-3; j++)
                {
                  coef = Px(i)*Py(j)*(1.0/Pz(k) - den);
                  phi(num).Init(0.0, 0.0, coef);
                  num++;
                }
          
          
          // part rot u avec u in Q_{k+1}^{[k-1, k]} x 0 x 0
          for (int k = 0; k <= order+1; k++)
            for (int i = 0; i <= k-2; i++)
              for (int j = 0; j <= k-3; j++)
                {
                  // coef = Px(i)*yi*(1.0-yi)*Py(j)/Pz(k);
                  // phi = rot coef + correction to cancel z-component for z = 0
                  phi(num).Init(0.0, Real_wp(k)*Px(i)*yi*(1.0-yi)*Py(j)/Pz(k+1),
                                ((1.0 - 2.0*yi)*Py(j)
                                 +2.0*yi*(1.0-yi)*dPy(j))*Px(i)*(-den+1.0/Pz(k)));
                  num++;
                }

          // part rot u avec u in 0 x Q_{k+1}^{[k, k-1]} x 0
          for (int k = 0; k <= order+1; k++)
            for (int i = 0; i <= k-3; i++)
              for (int j = 0; j <= k-2; j++)
                {
                  // coef = xi*(1.0-xi)*Px(i)*Py(j)/Pz(k);
                  // phi = rot coef + correction to cancel z-component for z = 0
                  phi(num).Init(Real_wp(k)*Py(j)*xi*(1.0-xi)*Px(i)/Pz(k+1), 0.0,
                                ((1.0 - 2.0*xi)*Px(i)
                                 +2.0*xi*(1.0-xi)*dPx(i))*Py(j)*(-den+1.0/Pz(k)));
                  num++;
                }
          
          // on applique 1/J DF
          for (int k = 0; k < nb_dof_loc; k++)
            Mlt(DF, phi(k), res(k));
        }
        break;
      case OPTIMAL :
	{
	  Real_wp xt = x/(1.0-z);
	  Real_wp yt = y/(1.0-z);
	  Real_wp zt = z;
	  
	  int num = 0;
	  for (int k = 0; k < r; k++)
	    for (int i = 0; i <= k; i++)
	      for (int j = 0; j <= k; j++)
		{
		  Real_wp val = pow(xt, i)*pow(yt, j)*pow(1.0-zt, k);
		  res(num).Init(val, 0.0, 0.0);
		  res(num+1).Init(0.0, val, 0.0);
		  res(num+2).Init(0.0, 0.0, val);
		  num += 3;
		}
	  
	  for (int k = 0; k <= r-1; k++)
	    for (int j = 0; j <= k; j++)
	      {
		Real_wp val = pow(xt, k+1)*pow(yt, j)*pow(1.0-zt, k);
		res(num).Init(val, 0.0, 0.0);
		val = pow(yt, k+1)*pow(xt, j)*pow(1.0-zt, k);
		res(num+1).Init(0.0, val, 0.0);
		num += 2;
	      }
	  
	  for (int k = 0; k <= r-1; k++)
	    for (int j = 0; j <= k; j++)
	      {
		res(num).Init(pow(xt, j+1)*pow(yt, k+1)*pow(1.0-zt, k),
			      0.0, -pow(xt, j)*pow(yt, k+1)*pow(1.0-zt, k));
		res(num+1).Init(0.0, pow(yt, j+1)*pow(xt, k+1)*pow(1.0-zt, k),
                                -pow(yt, j)*pow(xt, k+1)*pow(1.0-zt, k));		
		num += 2;
	      }
      
	  for (int i = 0; i <= r; i++)
	    for (int j = 0; j <= r; j++)
	      {
		Real_wp val = pow(xt, i)*pow(yt, j)*pow(1.0-z, r);
		res(num).Init(val*xt, val*yt, -val);
		num++;
	      }
	  
	}
	break;
      case OPTIMAL_HP :
	{
	  Real_wp omz, invOmz(0);
	  omz = 1.0 - z;
	  
	  int num = 0;
	  
	  // coordinates on the cube
	  Real_wp a(0), b(0), c;
	  if (abs(omz) > epsilon_machine)
	    {
	      invOmz = 1.0/(1.0-z);
	      a = x*invOmz;
	      b = y*invOmz;	      
	    }      
	  
	  c = 2.0*z - 1.0;
          
	  VectReal_wp powOneMinusZ(order+1);
	  powOneMinusZ(0) = 1.0;
	  for (int i = 0; i <= order-1; i++)
	    powOneMinusZ(i+1) = powOneMinusZ(i)*omz;
	  
	  // values of Legendre polynomials for edges
	  VectReal_wp Px, Py, Legz, Px2, Py2;
	  EvaluateJacobiPolynomial(LegendrePolynom, order, x, Px2);
	  EvaluateJacobiPolynomial(LegendrePolynom, order, y, Py2);
	  
	  EvaluateJacobiPolynomial(LegendrePolynom, order, a, Px);
	  EvaluateJacobiPolynomial(LegendrePolynom, order, b, Py);
	  
	  Vector<VectReal_wp> Pz(order);
	  for (int i = 0; i <= order-1; i++)
	    EvaluateJacobiPolynomial(EvenJacobiPolynom(i), order-1-i, c, Pz(i));
	  
	  EvaluateJacobiPolynomial(LegendrePolynom, order, c, Legz);
	  
	  R3 v1((2.0*(1.0+x-z) - x*z)/(1.0-z), -y*z/(1.0-z), z);
	  R3 v2((-2.0*(1.0-x-z) - x*z)/(1.0-z), -y*z/(1.0-z), z);
	  R3 v3(-x*z/(1.0-z), (2.0*(1.0+y-z) - y*z)/(1.0-z), z);
	  R3 v4(-x*z/(1.0-z), (-2.0*(1.0-y-z) - y*z)/(1.0-z), z);
	  R3 v5(-x, -y, (1.0-z));
		
	  R3 v6(x*z*y/(1.0-z), (1.0-z-y*y), -y*z);
	  R3 v7(1.0-z-x*x, x*z*y/(1.0-z), -x*z);
	  
	  for (int i = 0; i <= order; i++)
	    for (int j = 0; j <= order; j++)
	      {
		Real_wp val = Px(i)*Py(j)*powOneMinusZ(max(i, j))*invOmz;
		res(num) = val*v5;
		num++;
	      }
	  
	  for (int i = 0; i <= order-1; i++)
	    for (int j = 0; j <= order-1-i; j++)
	      {
		Real_wp val = Px2(i)*Legz(j);
		res(num++) = val*v4;
	      }
	  
	  for (int i = 0; i <= order-1; i++)
	    for (int j = 0; j <= order-1-i; j++)
	      {	
		Real_wp val = Py2(i)*Legz(j);
		res(num++) = val*v1;
	      }
	  
	  for (int i = 0; i <= order-1; i++)
	    for (int j = 0; j <= order-1-i; j++)
	      {	
		Real_wp val = Px2(i)*Legz(j);
		res(num++) = val*v3;
	      }
	  
	  for (int i = 0; i <= order-1; i++)
	    for (int j = 0; j <= order-1-i; j++)
	      {	
		Real_wp val = Py2(i)*Legz(j);
		res(num++) = val*v2;
	      }
	  
	  for (int i = 0; i <= order-1; i++)
	    for (int j = 0; j <= order-1; j++)
	      for (int k = 0; k <= order -1 - max(i, j); k++)
		{
		  int m = max(i, j);
		  Real_wp val = Px(i)*Py(j)*powOneMinusZ(m)*Legz(k)*invOmz;		  
		  res(num) = val*z*v5;
		  res(num+1) = val*v6;
		  res(num+2) = val*v7;
		  
		  num += 3;
		}
	}
      case ZAGLMAYR :
        {
	}
	break;
      }
  }
  
  
  //! computes \f$ curl(\varphi_{node})(point\_loc) \f$
  void PyramidHdivOther::ComputeDivPhiRef(const R3& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    FillZero(res);

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();

    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);
    
    switch (type_basis)
      {
      case GRAGLIA :
        {
          // on calcule xu, yu, zu de la pyramide unite
	  // x = 2 xu + zu - 1
	  // y = 2 yu + zu - 1
	  // z = zu
          Real_wp xu = 0.5*(1.0 + x - z);
          Real_wp yu = 0.5*(1.0 + y - z);
          Real_wp zu = z;
          
          // on calcule 1/J DF
          Matrix3_3 DF;
          DF(0, 0) = 2.0; DF(1, 1) = 2.0;
	  DF(0, 2) = 1.0; DF(1, 2) = 1.0; DF(2, 2) = 1.0;
	  Real_wp invJacob = 1.0/Det(DF);
          
          if (order == 1)
            {
              res(0) = 3.0;
              res(1) = 3.0;
              res(2) = 3.0;
              res(3) = -3.0;
              res(4) = -3.0;
            }
          else
            {
              R3 v0, v1, v2, v3, v4;
              v4.Init((-2.0*(1.0-xu-zu) - xu*zu)/(1.0-zu), -yu*zu/(1.0-zu), zu);
              v1.Init(-xu*zu/(1.0-zu), (-2.0*(1.0-yu-zu) -yu*zu)/(1.0-zu), zu);
              v2.Init(xu*(2.0-zu)/(1.0-zu), -yu*zu/(1.0-zu), zu);
              v3.Init(-xu*zu/(1.0-zu), yu*(2.0-zu)/(1.0-zu), zu);
              v0.Init(xu, yu, -(1.0-zu));
              
              v3 *= -1.0; v4 *= -1.0;
              
              res(0) = 3.0*((3.0*(1.0-xu) - 1.0)*(3.0*(1.0-yu) - 1.0) - 3.0*zu) 
                + (-3.0*(3.0*(1.0-yu)-1.0)*v0(0) - 3.0*(3.0*(1.0-xu) - 1.0)*v0(1) - 3.0*v0(2));
              
              res(1) = 3.0*(3.0*yu - 1.0)*(3.0*(1.0-xu) - 1.0) 
                + (-3.0*(3.0*yu - 1.0)*v0(0) + 3.0*(3.0*(1.0-xu)-1.0)*v0(1));
              
              res(2) = 3.0*(3.0*xu - 1.0)*(3.0*(1.0-yu) - 1.0)
                + ( 3.0*(3.0*(1.0-yu) - 1.0)*v0(0) - 3.0*(3.0*xu-1.0)*v0(1));
              
              res(3) = 3.0*(3.0*xu - 1.0)*(3.0*yu - 1.0)
                + ( 3.0*(3.0*yu-1.0)*v0(0) + 3.0*(3.0*xu-1.0)*v0(1));
              
              res(4) = 3.0*(4.0*(1.0-xu-zu) - yu - 1.0) - 4.0*v1(0) - v1(1) - 4.0*v1(2);
              res(5) = 3.0*(4.0*xu - yu - 1.0) + 4.0*v1(0) - v1(1);
              res(6) = 3.0*(4.0*zu - yu - 1.0) - v1(1) + 4.0*v1(2);
              
              res(7) = 3.0*(4.0*(1.0-yu-zu) - (1.0-xu-zu) - 1.0) + v2(0) - 4.0*v2(1) -3.0*v2(2);
              res(8) = 3.0*(4.0*yu - (1.0-xu-zu) - 1.0) + v2(0) + 4.0*v2(1) + v2(2);
              res(9) = 3.0*(4.0*zu - (1.0-xu-zu) - 1.0) + v2(0) + 5.0*v2(2);
              
              res(10) = -3.0*(4.0*(1.0-xu-zu) - (1.0-yu-zu) - 1.0) - 4.0*v3(0) + v3(1) - 3.0*v3(2);
              res(11) = -3.0*(4.0*xu - (1.0-yu-zu) - 1.0) + 4.0*v3(0) + v3(1) + v3(2);
              res(12) = -3.0*(4.0*zu - (1.0-yu-zu) - 1.0) + v3(1) + 5.0*v3(2);
              
              res(13) = -3.0*(4.0*(1.0-yu-zu)-xu-1.0) - v4(0) - 4.0*v4(1) - 4.0*v4(2);
              res(14) = -3.0*(4.0*yu - xu - 1.0) - v4(0) + 4.0*v4(1);
              res(15) = -3.0*(4.0*zu - xu - 1.0) - v4(0) + 4.0*v4(2);
              
              res(16) = -3.0*(3.0*xu) + 3.0*v4(0);
              res(17) = 3.0*(3.0*yu) + 3.0*v1(1);
              res(18) = 3.0*(3.0*zu) + 3.0*v0(2);
            }
          
          // on applique 1/J
          for (int k = 0; k < nb_dof_loc; k++)
            res(k) *= invJacob;
        }
        break;
      case GRADINARU :
        {          
          // on calcule xu, yu, zu de la pyramide unite
	  // x = 2 xu + zu - 1
	  // y = 2 yu + zu - 1
	  // z = zu
          //Real_wp xu = 0.5*(1.0 + x - z);
          //Real_wp yu = 0.5*(1.0 + y - z);
          //Real_wp zu = z;
          
          // on calcule 1/J DF
          Matrix3_3 DF;
          DF(0, 0) = 2.0; DF(1, 1) = 2.0;
	  DF(0, 2) = 1.0; DF(1, 2) = 1.0; DF(2, 2) = 1.0;
	  Real_wp invJacob = 1.0/Det(DF);
          
          res(0) = 3.0;
          res(1) = 3.0;
          res(2) = 3.0;
          res(3) = -3.0;
          res(4) = -3.0;
          
          // on applique 1/J
          for (int k = 0; k < nb_dof_loc; k++)
            res(k) *= invJacob;
        }
        break;
      case NIGAM_PHILLIPS_1 :
        {
          // x, y, z de la pyramide infinie
          Real_wp xi = 0.5*(x/(1.0-z) + 1.0);
          Real_wp yi = 0.5*(y/(1.0-z) + 1.0);
          Real_wp zi = z/(1.0-z);
			
          // DF^-1
          Matrix3_3 dfjm1, DF;
          dfjm1(0, 0) = 0.5/(1.0-z);
          dfjm1(1, 1) = 0.5/(1.0-z);
          dfjm1(2, 2) = 1.0/square(1.0-z);
          dfjm1(0, 2) = 0.5*x/square(1.0-z);
          dfjm1(1, 2) = 0.5*y/square(1.0-z);
	  GetInverse(dfjm1, DF);
	  Real_wp invJacob = 1.0/Det(DF);
          
          VectReal_wp Px, Py, Pz(order+1), dPx, dPy, dPz(order+2);
	  EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*xi-1.0, Px, dPx);
	  EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*yi-1.0, Py, dPy);
          dPz.Fill(0);
          for (int i = 0; i <= order; i++)
            {
              Pz(i) = pow(zi, i);
              dPz(i+1) = Real_wp(i+1)*Pz(i);
            }
          
          Real_wp den2 = 1.0/pow(1.0+zi, order+2);
          Real_wp den = Real_wp(order+2)/pow(1.0+zi, order+3);
          Real_wp coef(0);
          int num = 0;
	  for (int i = 0; i < order; i++)
            for (int j = 0; j < order; j++)
              {
                res(num) = -Px(i)*Py(j)*den;
                num++;
              }
          
	  for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                res(num) = Px(i)*Pz(j)*(-3.0*den2 + den*zi) - zi*den2*Px(i)*dPz(j);
                num++;
              }
          
          for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                res(num) = Py(i)*Pz(j)*(-3.0*den2 + den*zi) - zi*den2*Py(i)*dPz(j);
                num++;
              }
          
          for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                res(num) = Px(i)*Pz(j)*(3.0*den2 - den*zi)  + zi*den2*Px(i)*dPz(j);
                num++;
              }
          
          for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                res(num) = Py(i)*Pz(j)*(3.0*den2 - zi*den)  + zi*den2*Py(i)*dPz(j);
                num++;
              }
          
          // part x(1-x) \khi_1
          for (int i = 0; i <= order-2; i++)
            for (int j = 0; j <= order-1; j++)
              for (int k = 0; k <= order-2; k++)
                {
                  res(num) = Py(j)*Pz(k)*den2*((1.0-2.0*xi)*Px(i) + 2.0*dPx(i)*xi*(1.0-xi));
                  num++;
                }
          
          // part y(1-y) \khi_2
          for (int i = 0; i <= order-1; i++)
            for (int j = 0; j <= order-2; j++)
              for (int k = 0; k <= order-2; k++)
                {
                  res(num) = Px(i)*Pz(k)*den2*((1.0-2.0*yi)*Py(j) + 2.0*dPy(j)*yi*(1.0-yi));
                  num++;
                }
          
          // part z \khi_3
          for (int i = 0; i <= order-1; i++)
            for (int j = 0; j <= order-1; j++)
              for (int k = 0; k <= order-2; k++)
                {
                  res(num) = Px(i)*Py(j)*(-den*Pz(k+1) + Real_wp(k+1)*den2*Pz(k)); 
                  num++;
                }
          
          // part z^(k-1)/(1+z)^{k+2} (2t, 0, (1+z) t_x)
          Real_wp tx, sy, dcoef;
          for (int i = 0; i <= order-2; i++)
            for (int j = 0; j <= order-1; j++)
              {
                tx = (1.0 - 2.0*xi)*Px(i)*Py(j) + 2.0*xi*(1.0-xi)*dPx(i)*Py(j);
                coef = pow(zi, order-1)*den2;
                dcoef = Real_wp(order-1)*pow(zi, order-2)*den2 - pow(zi, order-1)*den;
                
                res(num) = tx*(3.0*coef + (1.0+zi)*dcoef);
                num++;
              }
          
          // part z^(k-1)/(1+z)^{k+2} (0, 2s, (1+z) s_y)
          for (int i = 0; i <= order-1; i++)
            for (int j = 0; j <= order-2; j++)
              {
                // s et s_y
                sy = (1.0 - 2.0*yi)*Px(i)*Py(j) + 2.0*yi*(1.0-yi)*Px(i)*dPy(j);
                coef = pow(zi, order-1)*den2;
                dcoef = Real_wp(order-1)*pow(zi, order-2)*den2 - pow(zi, order-1)*den;
                                
                res(num) = sy*(3.0*coef + (1.0+zi)*dcoef);
                num++;
              }
          
          // on applique 1/J
          for (int k = 0; k < nb_dof_loc; k++)
            res(k) *= invJacob;
        }
        break;
      case NIGAM_PHILLIPS_2 :
        {
          // x, y, z de la pyramide infinie
          Real_wp xi = 0.5*(x/(1.0-z) + 1.0);
          Real_wp yi = 0.5*(y/(1.0-z) + 1.0);
          Real_wp zi = z/(1.0-z);
			
          // DF^-1
          Matrix3_3 dfjm1, DF;
          dfjm1(0, 0) = 0.5/(1.0-z);
          dfjm1(1, 1) = 0.5/(1.0-z);
          dfjm1(2, 2) = 1.0/square(1.0-z);
          dfjm1(0, 2) = 0.5*x/square(1.0-z);
          dfjm1(1, 2) = 0.5*y/square(1.0-z);
	  GetInverse(dfjm1, DF);
	  Real_wp invJacob = 1.0/Det(DF);

          VectReal_wp Px, Py, Pz(order+4), dPx, dPy;
	  EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*xi-1.0, Px, dPx);
	  EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*yi-1.0, Py, dPy);
          for (int i = 0; i <= order+3; i++)
            Pz(i) = pow(1.0 + zi, i);
          
          // dofs associated with the quadrilateral base
          Real_wp den = Real_wp(order+2)/Pz(order+3);
          int num = 0;
          for (int i = 0; i <= order-1; i++)
            for (int j = 0; j <= order-1; j++)
              {
                res(num) = -Px(i)*Py(j)*den;
                num++;
              }
          
          // first triangular face
	  for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                res(num) = -Px(i)*den;
                num++;
              }

          // second triangular face
	  for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                res(num) = -Py(i)*den;
                num++;
              }

          // third triangular face
	  for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                res(num) = Px(i)*den;
                num++;
              }

          // last triangular face
	  for (int i = 0; i < order; i++)
            for (int j = 0; j < order-i; j++)
              {
                res(num) = Py(i)*den;
                num++;
              }

          // interior dofs
          for (int k = 0; k <= order+1; k++)
            for (int i = 0; i <= k-3; i++)
              for (int j = 0; j <= k-3; j++)
                {
                  res(num) = -Px(i)*Py(j)*(Real_wp(k)/Pz(k+1) - den);
                  num++;
                }
          
          for (int k = 0; k <= order+1; k++)
            for (int i = 0; i <= k-2; i++)
              for (int j = 0; j <= k-3; j++)
                {
                  res(num) = ((1.0 - 2.0*yi)*Py(j)+2.0*yi*(1.0-yi)*dPy(j))*Px(i)*den;
                  num++;
                }
          
          for (int k = 0; k <= order+1; k++)
            for (int i = 0; i <= k-3; i++)
              for (int j = 0; j <= k-2; j++)
                {
                  res(num) = ((1.0 - 2.0*xi)*Px(i)+2.0*xi*(1.0-xi)*dPx(i))*Py(j)*den;
                  num++;
                }
          
          // on applique 1/J
          for (int k = 0; k < nb_dof_loc; k++)
            res(k) *= invJacob;
        }
        break;
      case OPTIMAL :
	{
	}
	break;
      case OPTIMAL_HP :
	{
	}
      case ZAGLMAYR :
        {
	}
	break;
      }

  }

  
  void PyramidHdivOther::ComputeValuesPhiOrthoRef(const R3& pointloc, VectR3& res) const
  {
    VectR3 phi;
    phi.Reallocate(nb_dof_loc);
    res.Reallocate(nb_dof_loc);

    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    
    Real_wp xi = 0.5*(x/(1.0-z) + 1.0);
    Real_wp yi = 0.5*(y/(1.0-z) + 1.0);
    Real_wp zi = z/(1.0-z);

    // DF^-1
    Matrix3_3 dfjm1, DF;
    dfjm1(0, 0) = 0.5/(1.0-z);
    dfjm1(1, 1) = 0.5/(1.0-z);
    dfjm1(2, 2) = 1.0/square(1.0-z);
    dfjm1(0, 2) = 0.5*x/square(1.0-z);
    dfjm1(1, 2) = 0.5*y/square(1.0-z);
    GetInverse(dfjm1, DF);
    Real_wp invJacob = 1.0/Det(DF);
    Mlt(invJacob, DF);

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    
    VectReal_wp Px, Py, Pz(order+4), dPx, dPy;
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*xi-1.0, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*yi-1.0, Py, dPy);
    for (int i = 0; i <= order+3; i++)
      Pz(i) = pow(1.0 + zi, i);
    
    // part (0, 0, Q_{k+2}^{[k-1, k-1]}
    int num = 0;
    Real_wp coef, dcoef_dx, dcoef_dy, dcoef_dz;
    
    if (type_basis == NIGAM_PHILLIPS_1)
      {
        Real_wp den = 1.0/Pz(order+2);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order-1; j++)
            for (int k = 0; k <= order-2; k++)
              {
                phi(num).Init(Px(i)*Py(j)*pow(zi, k)*den, 0, 0);
                num++;
              }
        
        for (int i = 0; i <= order-1; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order-2; k++)
              {
                phi(num).Init(0, Px(i)*Py(j)*pow(zi, k)*den, 0);
                num++;
              }
        
        for (int i = 0; i <= order-1; i++)
          for (int j = 0; j <= order-1; j++)
            for (int k = 0; k <= order-1; k++)
              {
                phi(num).Init(0, 0, Px(i)*Py(j)*pow(zi, k)*den);
                num++;
              }
        
        for (int i = 0; i <= order-1; i++)
          for (int j = 0; j <= order; j++)
            {
              coef = pow(zi, order-1)*den;
              phi(num).Init(0, 2.0*Px(i)*Py(j)*coef, 2.0*Px(i)*dPy(j)*(1.0+zi)*coef);
              num++;
            }
        
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order-1; j++)
            {
              coef = pow(zi, order-1)*den;
              phi(num).Init(2.0*Px(i)*Py(j)*coef, 0.0, 2.0*dPx(i)*Py(j)*(1.0+zi)*coef);
              num++;
            }
      }
    else
      {
        // interior dofs
        for (int k = 0; k <= order+2; k++)
          for (int i = 0; i <= k-3; i++)
            for (int j = 0; j <= k-3; j++)
              {
                coef = Px(i)*Py(j)*(1.0/Pz(k));
                phi(num).Init(0.0, 0.0, coef);
                num++;
              }
        
        // part rot u avec u in Q_{k+1}^{[k-1, k]} x 0 x 0
        for (int k = 0; k <= order+1; k++)
          for (int i = 0; i <= k-2; i++)
            for (int j = 0; j <= k-1; j++)
              {
                coef = Px(i)*Py(j)/Pz(k);
                dcoef_dy = 2.0*dPy(j)*Px(i)/Pz(k);
                dcoef_dz = -Real_wp(k)*Px(i)*Py(j)/Pz(k+1);
                phi(num).Init(0.0, dcoef_dz, -dcoef_dy);
                num++;
              }
        
        // part rot u avec u in 0 x Q_{k+1}^{[k, k-1]} x 0
        for (int k = 0; k <= order+1; k++)
          for (int i = 0; i <= k-1; i++)
            for (int j = 0; j <= k-2; j++)
              {
                coef = Px(i)*Py(j)/Pz(k);
                dcoef_dx = 2.0*dPx(i)*Py(j)/Pz(k);
                dcoef_dz = -Real_wp(k)*Px(i)*Py(j)/Pz(k+1);
                phi(num).Init(-dcoef_dz, 0.0, dcoef_dx);
                num++;
              }
      }
    
    // on applique 1/J DF
    for (int k = 0; k < nb_dof_loc; k++)
      Mlt(DF, phi(k), res(k));

  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \cdot \varphi_i \f$ for all i
  /*!
    \param[in] feval values of f on quadrature points of the face
    \param[out] res result vector (dof components)
    \param[in] num_loc local face number
   */
  template<class Vector1, class Vector2>
  void PyramidHdivOther::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    ComputeIntegralSurfaceRef(feval, res, num_loc);
  }

  void PyramidHdivOther
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  void PyramidHdivOther
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }


  void PyramidHdivOther
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
    // to be corrected
    cout << "Projection on surface dofs is not working for this element" << endl;
    abort();
  }


  void PyramidHdivOther
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {    
    // to be corrected
    cout << "Projection on surface dofs is not working for this element" << endl;
    abort();
  }

  
  ostream& operator <<(ostream& out, const PyramidHdivOther& e)
  {
    out<<static_cast<const PyramidReference<3>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_HDIV_OTHER_CXX
#endif
