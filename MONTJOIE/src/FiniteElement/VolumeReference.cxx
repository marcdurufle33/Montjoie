#ifndef MONTJOIE_FILE_VOLUME_REFERENCE_CXX

namespace Montjoie
{
  //! returns size of memory used by the object
  template<int type>
  size_t VolumeReference<type>::GetMemorySize() const
  {
    size_t taille = ElementReference_Dim<Dimension3>::GetMemorySize();
    taille += ElementReferenceType<Dimension3, type>::GetMemorySize();
    taille += sizeof(int)*(FacesDof.GetDataSize()+FacesNode.GetDataSize());
    for (int i = 0; i < ProjOperatorTriOrder.GetM(); i++)
      taille += sizeof(Real_wp)*ProjOperatorTriOrder(i).GetDataSize();

    for (int i = 0; i < ProjOperatorQuadOrder.GetM(); i++)
      taille += sizeof(Real_wp)*ProjOperatorQuadOrder(i).GetDataSize();
    
    taille += sizeof(*this) - sizeof(ElementReference_Dim<Dimension3>);
    return taille;
  }
  
  
  /*****************
   * Other methods *
   *****************/

  
  //! computation of projection between finite element of different orders
  template<int type>
  void VolumeReference<type>::
  SetVariableOrder(const Mesh<Dimension3>& mesh,
		   const MeshNumbering<Dimension3>& mesh_num)
  {
    if (!mesh_num.IsOrderVariable())
      return;
    
    int rmax = 0;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      rmax = max(rmax, mesh_num.GetOrderElement(i));

    for (int i = 0; i < mesh.GetNbBoundary(); i++)
      rmax = max(rmax, mesh_num.GetOrderQuadrature(i));
    
    Vector<bool> order_quad(rmax + 1), order_tri(rmax+1);
    order_quad.Fill(false); order_tri.Fill(false);
    int type_elt = this->GetHybridType();
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	int r = mesh_num.GetOrderElement(i);
	if ((r == this->order) && (type_elt == mesh.Element(i).GetHybridType()))
	  {
	    // we are on an element with the order of this finite element
	    for (int j = 0; j < mesh.Element(i).GetNbFaces(); j++)
	      {
		int ne = mesh.Element(i).numFace(j);
		int re = mesh_num.GetOrderQuadrature(ne);
		if (re != r)
		  {
		    if (mesh.Boundary(ne).GetNbVertices() == 3)
		      order_tri(re) = true;
		    else
		      order_quad(re) = true;
		  }
	      }
	  }
      }

    // elements of other processors
    for (int i = 0; i < mesh_num.GetNbNeighborElt(); i++)
      {
        int r = mesh_num.GetOrderNeighborElement(i);
        if ( (r == this->order) && (type_elt = mesh_num.GetTypeNeighborElement(i)))
          {
            int ne = mesh_num.GetEdgeNeighborElement(i);
            int re = mesh_num.GetOrderQuadrature(ne);
            
            if (re != r)
              {
                if (mesh.Boundary(ne).GetNbVertices() == 3)
                  order_tri(re) = true;
                else
                  order_quad(re) = true;
              }
          }
      }
    
    int nb_order_tri = 0, nb_order_quad = 0;
    for (int r = 1; r <= rmax; r++)
      {
	if (order_quad(r))
	  nb_order_quad++;

	if (order_tri(r))
	  nb_order_tri++;
      }
    
    if ((nb_order_tri > 0) && (this->GetHybridType() != 3))
      {
	IVect order_elt(nb_order_tri);
	Vector<VectR2> points_quad(nb_order_tri); 
	VectR2 points; VectReal_wp poids;
	int type_quad = this->GetTypeIntegrationTriangle();
	nb_order_tri = 0;
	for (int r = 1; r <= rmax; r++)
	  if (order_tri(r))
	    {
	      TriangleQuadrature::ConstructQuadrature(2*r, points, poids, type_quad);
	      points_quad(nb_order_tri) = points;
	      order_elt(nb_order_tri) = r;
	      nb_order_tri++;
	    }
	
	if (type == 1)
          this->ComputeTriangularInterpolationProjectorOrder(order_elt, points_quad);
      }

    if ((nb_order_quad > 0) && (this->GetHybridType() != 0))
      {
	IVect order_elt(nb_order_quad);
	Vector<VectR2> points_quad(nb_order_quad); 
	VectR2 points, poids;
	int type_quad = this->GetTypeIntegrationQuadrangle();
	nb_order_quad = 0;
	for (int r = 1; r <= rmax; r++)
	  if (order_quad(r))
	    {
	      VectReal_wp points1d_, weights1d_;
	      Globatto<Real_wp> gauss;
	      gauss.ConstructQuadrature(r, type_quad);
	      points1d_ = gauss.Points(); weights1d_ = gauss.Weights();
	      Matrix<int> NumQuad, coor;
	      MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(r, NumQuad, coor);
	      poids.Reallocate((r+1)*(r+1));
	      points.Reallocate((r+1)*(r+1));
	      for (int i = 0; i <= r; i++)
		for (int j = 0; j <= r; j++)
		  {
		    points(NumQuad(i, j)).Init(points1d_(i), points1d_(j));
		    poids(NumQuad(i, j)) = weights1d_(i)*weights1d_(j);
		  }

	      points_quad(nb_order_quad) = points;
	      order_elt(nb_order_quad) = r;
	      nb_order_quad++;
	    }
	
	if (type == 1)
          this->
            ComputeQuadrangularInterpolationProjectorOrder(order_elt, points_quad);
      }
  }


  //! computation of projection between finite element of different orders
  template<int type>
  void VolumeReference<type>::
  ComputeTriangularInterpolationProjectorOrder(const IVect& order_elt,
					       const Vector<VectR2>& Pts)
  {
    cout << "Not implemented " << endl;
    abort();
  }


  //! computation of projection between finite element of different orders
  template<int type>
  void VolumeReference<type>::
  ComputeQuadrangularInterpolationProjectorOrder(const IVect& order_elt,
						 const Vector<VectR2>& Pts)
  {
    cout << "Not implemented " << endl;
    abort();
  }
    

  //! computation of u on nodal points of the element
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points
   */
  template<int type> template<class Vector1, class Vector2>
  void VolumeReference<type>
  ::ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const
  {
    Unode.Reallocate(this->ValuePhi_Nodal.GetN());
    Mlt(SeldonTrans, this->ValuePhi_Nodal, Un, Unode);
  }
  
  
  //! displays information about class VolumeReference<type>
  template<int type>
  ostream& operator <<(ostream& out, const VolumeReference<type>& e)
  {
    return out;
  }

  
  /****************
   * H^1 elements *
   ****************/
  
  
#ifdef MONTJOIE_WITH_NODAL_H1
  //! Destructor
  ElementReference<Dimension3, 1>::~ElementReference()
  {
    if (element_tri_surf != NULL)
      {
	delete element_tri_surf;
	element_tri_surf = NULL;
      }

    if (element_quad_surf != NULL)
      {
	delete element_quad_surf;
	element_quad_surf = NULL;
      }
  }
  


  //! returns 2-D triangular finite element (traces on triangular faces)
  const ElementReference<Dimension2, 1>& 
  ElementReference<Dimension3, 1>::GetTriangularSurfaceFiniteElement() const
  {
    if (element_tri_surf == NULL)
      {
	cout << "Null pointer for GetTriangularSurfaceFiniteElement" << endl;
	abort();
      }

    return *element_tri_surf;
  }


  //! returns 2-D quadrangular finite element (traces on quadrangular faces)
  const ElementReference<Dimension2, 1>& 
  ElementReference<Dimension3, 1>::GetQuadrangularSurfaceFiniteElement() const
  {
    if (element_quad_surf == NULL)
      {
	cout << "Null pointer for GetQuadrangularSurfaceFiniteElement" << endl;
	abort();
      }

    return *element_quad_surf;
  }


  const ElementReference<Dimension2, 1>& 
  ElementReference<Dimension3, 1>::GetSurfaceFiniteElement(int n) const
  {
    if (this->elt_geom.IsLocalFaceQuadrilateral(n))
      return GetQuadrangularSurfaceFiniteElement();
    else
      return GetTriangularSurfaceFiniteElement();
  }
  
  
  size_t ElementReference<Dimension3, 1>::GetMemorySize() const
  {
    size_t taille = VolumeReference<1>::GetMemorySize() + sizeof(FiniteElementH1<Dimension3>);
    if (element_tri_surf != NULL)
      taille += element_tri_surf->GetMemorySize();

    if (element_quad_surf != NULL)
      taille += element_quad_surf->GetMemorySize();
    
    return taille;
  }
  

  //! computation of projection between finite element of different orders
  void ElementReference<Dimension3, 1>::
  ComputeTriangularInterpolationProjectorOrder(const IVect& order_elt,
					       const Vector<VectR2>& Pts)
  {
    ComputeTriangularInterpolationProjector(order_elt,
                                            this->GetTriangularSurfaceFiniteElement(), Pts);
  }


  //! computation of projection between finite element of different orders
  void ElementReference<Dimension3, 1>::
  ComputeQuadrangularInterpolationProjectorOrder(const IVect& order_elt,
						 const Vector<VectR2>& Pts)
  {
    ComputeQuadrangularInterpolationProjector(order_elt,
                                              this->GetQuadrangularSurfaceFiniteElement(), Pts);
  }

    
  //! computation of projection between finite element of different orders
  void ElementReference<Dimension3, 1>::
  ComputeTriangularInterpolationProjector(const IVect& order_elt,
                                          const ElementReference<Dimension2, 1>& Fb,
                                          const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi;
    ProjOperatorTriOrder.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorTriOrder(r).Reallocate(Fb.GetNbDof(), Pts(i).GetM());
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
	    Fb.ComputeValuesPhiRef(Pts(i)(j), phi);
	    for (int k = 0; k < Fb.GetNbDof(); k++)
	      ProjOperatorTriOrder(r)(k, j) = phi(k);
	  }
      }
  }
  
  
  //! computation of projection between finite element of different orders
  void ElementReference<Dimension3, 1>::
  ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                            const ElementReference<Dimension2, 1>& Fb,
                                            const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi;
    ProjOperatorQuadOrder.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorQuadOrder(r).Reallocate(Fb.GetNbDof(), Pts(i).GetM());
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
	    Fb.ComputeValuesPhiRef(Pts(i)(j), phi);
	    for (int k = 0; k < Fb.GetNbDof(); k++)
	      ProjOperatorQuadOrder(r)(k, j) = phi(k);
	  }
      }
  }

  
  //! computation of values of u on nodal points of a face
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points of face num_loc
    \param[in] num_loc face number
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 1>
  ::ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    for (int j = 0; j < Unode.GetM(); j++)
      {
        Unode(j) = 0;
        int node = FacesNode(j, num_loc);
        for (int k = 0; k < this->nb_dof_loc; k++)
          Unode(j) += Un(k)*this->ValuePhi_Nodal(k, node);
      }
  }


  //! computation of gradient of u on nodal points of a face
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode gradient of u on nodal points of face num_loc
    \param[in] num_loc face number
  */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 1>
  ::ComputeGradientBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    for (int j = 0; j < Unode.GetM()/3; j++)
      {
        Unode(3*j) = 0; Unode(3*j+1) = 0; Unode(3*j+2) = 0;
        int node = FacesNode(j, num_loc);
        for (int k = 0; k < this->nb_dof_loc; k++)
          {
            Unode(3*j) += Un(k)*this->GradPhi_Nodal(k, node)(0);
            Unode(3*j+1) += Un(k)*this->GradPhi_Nodal(k, node)(1);
            Unode(3*j+2) += Un(k)*this->GradPhi_Nodal(k, node)(2);
          }
      }
  }


  //! computing operator to apply when the faces are rotated (for continuous elements)
  void ElementReference<Dimension3, 1>::
  FindH1LinearCombinationRotation(NumberMap& nmap,
                                  const ElementReference<Dimension2, 1>& Fb_tri,
				  const ElementReference<Dimension2, 1>& Fb_quad,
                                  bool display_message) const
  {
    // vertex numbers of triangular face after rotation
    TinyVector<TinyVector<int, 3>, 6> rot_tri_num;
    rot_tri_num(0).Init(0, 1, 2);
    rot_tri_num(1).Init(1, 2, 0);
    rot_tri_num(2).Init(2, 0, 1);
    rot_tri_num(3).Init(0, 2, 1);
    rot_tri_num(4).Init(1, 0, 2);
    rot_tri_num(5).Init(2, 1, 0);
    
    // vertex numbers of quadrangular face after rotation
    TinyVector<TinyVector<int, 4>, 8> rot_quad_num;
    rot_quad_num(0).Init(0, 1, 2, 3);
    rot_quad_num(1).Init(1, 2, 3, 0);
    rot_quad_num(2).Init(2, 3, 0, 1);
    rot_quad_num(3).Init(3, 0, 1, 2);
    rot_quad_num(4).Init(0, 3, 2, 1);
    rot_quad_num(5).Init(1, 0, 3, 2);
    rot_quad_num(6).Init(2, 1, 0, 3);
    rot_quad_num(7).Init(3, 2, 1, 0);
    
    Matrix<int> IndexDofFace(this->GetNbDof(), this->GetNbBoundaries());
    IndexDofFace.Fill(-1);
    bool first_tri = true, first_quad = true;
    int type_elt = this->GetHybridType();
    Real_wp threshold = 1e4*epsilon_machine;
    for (int num_loc = 0; num_loc < this->GetNbBoundaries(); num_loc++)
      {
        bool quad = MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(type_elt, num_loc);
        bool check_face = false;
        if (quad)
          {
            if (first_tri)
              {
                first_tri = false;
                check_face = true;
              }
          }
        else
          {
            if (first_quad)
              {
                first_quad = false;
                check_face = true;
              }
          }
        
        // check_face = true;
        if (check_face)
          {
            Matrix<int> FaceRot;
            Matrix<bool> OppositeSigneDof;
            bool linear_combination = false;
            
            // dof_on_face(i) will be true for each dof strictly inside face num_loc
            Vector<bool> dof_on_face(this->GetNbDof()); dof_on_face.Fill(true);
            for (int i = this->GetNbDofBoundaries(); i < this->GetNbDof(); i++)
              dof_on_face(i) = false;
            
            int N = this->GetNbQuadBoundary(num_loc);
            VectReal_wp phi, phi1, phi2;
            for (int num_loc2 = 0; num_loc2 < this->GetNbBoundaries(); num_loc2++)
              if (num_loc2 != num_loc)
                for (int i = 0; i < this->GetNbQuadBoundary(num_loc2); i++)
                  {
                    int k = this->GetQuadNumber(num_loc2, i);
                    this->GetValuePhiOnQuadraturePoint(k, phi);
                    
                    for (int j = 0; j < this->GetNbDof(); j++)
                      { 
                        if (abs(phi(j)) > R3::threshold)
                          {
                            dof_on_face(j) = false;
                          }
                      }
                  }
            
            IVect IndexDof(this->GetNbDof()); IndexDof.Fill(-1);
            int nb_inside = 0;
            for (int j = 0; j < this->GetNbDof(); j++)
              if (dof_on_face(j))
                {
                  IndexDof(j) = nb_inside;
                  IndexDofFace(j, num_loc) = nb_inside;
                  nb_inside++;
                }
            
            int nb_orientations = 6; int nb_vert_on_face = 3;
            if (quad)
              {
                nb_orientations = 8;
                nb_vert_on_face = 4;
              }
            
            if (nb_inside > 0)
              {    
                FaceRot.Reallocate(nb_orientations, nb_inside); FaceRot.Fill(-1); 
                OppositeSigneDof.Reallocate(nb_orientations, nb_inside);
                OppositeSigneDof.Fill(false);
                Vector<Matrix<Real_wp> > CoefCombination(nb_orientations);
                
                // small mesh with one element
                Mesh<Dimension3> mesh;
                mesh.ReallocateVertices(this->GetNbVertices());
                
                switch (type_elt)
                  {
                  case 0 :
                    {
                      mesh.Vertex(0).Init(0, 0, 0);
                      mesh.Vertex(1).Init(1, 0, 0);
                      mesh.Vertex(2).Init(0, 1, 0);
                      mesh.Vertex(3).Init(0, 0, 1);
                    }
                    break;
                  case 1 :
                    {
                      mesh.Vertex(0).Init(-1, -1, 0);
                      mesh.Vertex(1).Init(1, -1, 0);
                      mesh.Vertex(2).Init(1, 1, 0);
                      mesh.Vertex(3).Init(-1, 1, 0);
                      mesh.Vertex(4).Init(0, 0, 1);
                    }
                    break;
                  case 2 :
                    {
                      mesh.Vertex(0).Init(0, 0, 0);
                      mesh.Vertex(1).Init(1, 0, 0);
                      mesh.Vertex(2).Init(0, 1, 0);
                      mesh.Vertex(3).Init(0, 0, 1);
                      mesh.Vertex(4).Init(1, 0, 1);
                      mesh.Vertex(5).Init(0, 1, 1);
                    }
                    break;
                  case 3 :
                    {
                      mesh.Vertex(0).Init(0, 0, 0);
                      mesh.Vertex(1).Init(1, 0, 0);
                      mesh.Vertex(2).Init(1, 1, 0);
                      mesh.Vertex(3).Init(0, 1, 0);
                      mesh.Vertex(4).Init(0, 0, 1);
                      mesh.Vertex(5).Init(1, 0, 1);
                      mesh.Vertex(6).Init(1, 1, 1);
                      mesh.Vertex(7).Init(0, 1, 1);
                    }
                  }
                
                mesh.ReallocateElements(1);
                int ref = 1; IVect num(this->GetNbVertices());
                num.Fill();
                mesh.Element(0).Init(num, ref);
                mesh.FindConnectivity();
                
                VectR3 s, sFace;
                R3 pt, ptGlob1, ptGlob2;
                VectReal_wp phi_2d;
                VectReal_wp phi_glob1(this->GetNbDof()), phi_glob2;
                // we test all orientations
                for (int rot = 0; rot < nb_orientations; rot++)
                  {
                    IVect num(3);
                    if (quad)
                      num.Reallocate(4);
                    
                    sFace.Reallocate(num.GetM());
                    for (int j = 0; j < nb_vert_on_face; j++)
                      {
                        int jrot = -1;
                        if (quad)
                          jrot = rot_quad_num(rot)(j);
                        else
                          jrot = rot_tri_num(rot)(j);
                        
                        int nv = MeshNumbering<Dimension3>::
                          GetVertexNumberOfFace(type_elt, num_loc, j);
                        
                        num(jrot) = nv;
                        sFace(jrot) = mesh.Vertex(nv);
                      }

                    mesh.GetVerticesElement(0, s);
                    
                    // on boucle sur les points de quadrature de la face            
                    Matrix<Real_wp> Mh(nb_inside, nb_inside), Bh(nb_inside, nb_inside);
                    Mh.Fill(0); Bh.Fill(0); Real_wp x, y;
                    for (int ks = 0; ks < N; ks++)
                      {
                        Real_wp poids = this->WeightsQuadratureBoundary(ks, num_loc);
                        
                        // treating 3-D face of the element
                        int k = this->GetQuadNumber(num_loc, ks);
                        pt = this->PointsND(k);
                        this->FjLinear(s, pt, ptGlob1);
                        
                        this->GetValuePhiOnQuadraturePoint(k, phi1);
                        
                        // treating 3-D face from a reference triangle/quadrangle
                        int k2 = -1;
                        if (quad)
                          {
                            // searching point on 3-D face
                            for (int i2 = 0; i2 < Fb_quad.GetNbPointsQuadratureInside(); i2++)
                              {
                                x = Fb_quad.PointsND(i2)(0);
                                y = Fb_quad.PointsND(i2)(1);
                                ptGlob2 = (1.0-x)*(1.0-y)*sFace(0) + x*(1.0-y)*sFace(1)
                                  + x*y*sFace(2) + (1.0-x)*y*sFace(3);
                                
                                if (ptGlob2.Distance(ptGlob1) < threshold)
                                  k2 = i2;                            
                              }
                            
                            Fb_quad.GetValuePhiOnQuadraturePoint(k2, phi2);
                          }
                        else
                          {
                            // searching point on 3-D face
                            for (int i2 = 0; i2 < Fb_tri.GetNbPointsQuadratureInside(); i2++)
                              {
                                x = Fb_tri.PointsND(i2)(0);
                                y = Fb_tri.PointsND(i2)(1);
                                ptGlob2 = (1.0-x-y)*sFace(0) + x*sFace(1) + y*sFace(2);
                                if (ptGlob2.Distance(ptGlob1) < threshold)
                                  k2 = i2;                            
                              }
                            
                            Fb_tri.GetValuePhiOnQuadraturePoint(k2, phi2);
                          }
                        
                        // mise a jour matrice M et B
                        for (int iv = 0; iv < this->GetNbDof(); iv++)
                          for (int jv = 0; jv < this->GetNbDof(); jv++)
                            {
                              int i = IndexDof(iv);
                              int j = IndexDof(jv);
                              if ((i>= 0) && (j>= 0))
                                Mh(i, j) += poids*phi1(iv)*phi1(jv);
                            }
                        
                        for (int iv = 0; iv < this->GetNbDof(); iv++)
                          for (int jv = phi2.GetM()-nb_inside; jv < phi2.GetM(); jv++)
                            {
                              int i = IndexDof(iv);
                              int j = jv - phi2.GetM()+nb_inside;
                              if ((i>= 0) && (j>= 0))
                                Bh(i, j) += poids*phi1(iv)*phi2(jv);
                            }
                      }
                    
                    GetInverse(Mh);
                    Matrix<Real_wp> Ah(nb_inside, nb_inside);            
                    Ah.Fill(0);
                    Mlt(Mh, Bh, Ah);
                    
                    // DISP(num_loc); DISP(rot); DISP(Ah);
                    if (rot == 0)
                      {
                        bool test_identity = true;
                        for (int i = 0; i < nb_inside; i++)
                          for (int j = 0; j < nb_inside; j++)
                            {
                              if (j == i)
                                {
                                  if (abs(Ah(i, j) - 1.0) > threshold)
                                    test_identity = false;
                                }
                              else
                                {
                                  if (abs(Ah(i, j)) > threshold)
                                    test_identity = false;
                                }
                            }
                        
                        if (display_message)
                          if (!test_identity)
                            {
                              cout << "Warning : Operator not equal to identity for face "
                                 << num_loc << endl;
                              DISP(Ah);
                            }
                        
                      }
                    
                    Real_wp threshold = 1000*epsilon_machine;
                    // pour tous les ddls de la face
                    Vector<bool> DofUsed(nb_inside); DofUsed.Fill(false);            
                    for (int iv = 0; iv < this->GetNbDof(); iv++)
                      {
                        int i = IndexDof(iv);
                        if (i>= 0)
                          {
                            int nb_coef = 0;
                            for (int j = 0; j < nb_inside; j++)
                              if (abs(Ah(j, i)) > threshold)
                                nb_coef++;
                            
                            if (nb_coef == 0)
                              {
                                cout << "Cas impossible " << endl;
                                abort();
                              }
                            else if (nb_coef == 1)
                              {
                                // on est tombe sur un autre ddl, avec eventuellement un signe
                                for (int j = 0; j < nb_inside; j++)
                                  if (abs(Ah(j, i)) > threshold)
                                    {
                                      FaceRot(rot, i) = j;
                                      if ( abs(Ah(j, i) -1.0) <= threshold)
                                        OppositeSigneDof(rot, i) = false;
                                      else if ( abs(Ah(j, i) + 1.0) <= threshold)
                                        OppositeSigneDof(rot, i) = true;
                                      else
                                        {
                                          // this dof can't be handled with a sign
                                          linear_combination = true;
                                        }
                                      
                                      DofUsed(j) = true;
                                    }                       
                              }
                            else
                              {
                                // no direct relation between dofs
                                linear_combination = true;
                              }
                          }
                      }
                    
                    CoefCombination(rot) = Ah;
                    
                  }
                
                // DISP(CoefCombination);
                //DISP(linear_combination);
                //DISP(FaceRot);
                //DISP(OppositeSigneDof);
                
                if (linear_combination)
                  {
                    if (quad)
                      nmap.SetFacesDofRotationQuad(this->GetOrder(), CoefCombination);
                    else
                      nmap.SetFacesDofRotationTri(this->GetOrder(), CoefCombination);
                  }
                else
                  {
                    if (quad)
                      {
                        nmap.SetFacesDofRotationQuad(this->GetOrder(), FaceRot);
                        nmap.SetSignDofRotationQuad(this->GetOrder(), OppositeSigneDof);
                      }
                    else
                      {
                        nmap.SetFacesDofRotationTri(this->GetOrder(), FaceRot);
                        nmap.SetSignDofRotationTri(this->GetOrder(), OppositeSigneDof);
                      }
                  }
                
              }
          }
      }
  }
  
  
  //! computes the sign of dofs after symmetry of edges
  /*!
    This function is only a function to check that restriction of functions
    on edges are the same between the edges, and that a symmetry
    of an edge may change the sign only.
    This function does not modify the finite element class, neither nmap
   */
  void ElementReference<Dimension3, 1>::FindH1SignEdge(NumberMap& nmap) const
  {
    VectR3 s;
    switch (this->GetNbVertices())
      {
      case 4 :
        {
          s.Reallocate(4);
          s(0).Init(0, 0, 0);
          s(1).Init(1, 0, 0);
          s(2).Init(0, 1, 0);
          s(3).Init(0, 0, 1);
        }
        break;
      case 5 :
        {
          s.Reallocate(5);
          s(0).Init(-1, -1, 0);
          s(1).Init(1, -1, 0);
          s(2).Init(1, 1, 0);
          s(3).Init(-1, 1, 0);
          s(4).Init(0, 0, 1);
        }
        break;
      case 6 :
        {
          s.Reallocate(6);
          s(0).Init(0, 0, 0);
          s(1).Init(1, 0, 0);
          s(2).Init(0, 1, 0);
          s(3).Init(0, 0, 1);
          s(4).Init(1, 0, 1);
          s(5).Init(0, 1, 1);
        }
        break;
      case 8 :
        {
          s.Reallocate(8);
          s(0).Init(0, 0, 0);
          s(1).Init(1, 0, 0);
          s(2).Init(1, 1, 0);
          s(3).Init(0, 1, 0);
          s(4).Init(0, 0, 1);
          s(5).Init(1, 0, 1);
          s(6).Init(1, 1, 1);
          s(7).Init(0, 1, 1);
        }
        break;
      }

    int nb_quad = this->Points1D().GetM();
    int nb_dof_edge = nmap.GetNbDofEdge(this->order);
    if (nb_dof_edge == 0)
      return;
    
    Array3D<Real_wp> ValPhiEdge(this->GetNbEdges(), nb_dof_edge, nb_quad);
    ValPhiEdge.Fill(0);
    Real_wp threshold = 1e4*epsilon_machine;
    
    // computing tangent for each edge, and values of basis functions
    VectReal_wp phi;
    int type_elt = this->GetHybridType();
    int nb_vertices = this->GetNbVertices();
    for (int num_loc = 0; num_loc < this->GetNbEdges(); num_loc++)
      {
        int n1 = MeshNumbering<Dimension3>::FirstExtremityEdge(type_elt, num_loc);
        int n2 = MeshNumbering<Dimension3>::SecondExtremityEdge(type_elt, num_loc);
        for (int k = 0; k < nb_quad; k++)
          {
            R3 ptA = (1.0-this->Points1D(k))*s(n1) + this->Points1D(k)*s(n2);
            this->ComputeValuesPhiRef(ptA, phi);
            for (int j = 0; j < nb_dof_edge; j++)
              {
                int num_dof = num_loc*nb_dof_edge + j + nb_vertices;
                ValPhiEdge(num_loc, j, k) = phi(num_dof);
                for (int num_loc2 = 0; num_loc2 < this->GetNbEdges(); num_loc2++)
                  if (num_loc2 != num_loc)
                    {
                      int num_dof = num_loc2*nb_dof_edge + j + nb_vertices;
                      Real_wp val = phi(num_dof);
                      
                      if (abs(val) > threshold)
                        {
                          cout << "Basis function "<<num_dof << " is not null on edge "
                               << num_loc2 << endl;
                          
                          abort();
                        }
                    }
              }
          }
      }
    
    //DISP(ValPhiEdge);
    // we check tangential traces on each edge
    for (int num_loc = 0; num_loc < this->GetNbEdges(); num_loc++)
      {
        Matrix<Real_wp> A(nb_dof_edge, nb_dof_edge), B = A, M = A;
        A.Fill(0); B.Fill(0); M.Fill(0);
        for (int j = 0; j < nb_dof_edge; j++)
          {
            // we check that it is the same as first edge
            for (int k = 0; k < nb_quad; k++)
              if (abs(ValPhiEdge(num_loc, j, k) - ValPhiEdge(0, j, k)) > threshold)
                {
                  cout << "values of tangential trace of function " << j 
                       << " on edge " <<  num_loc 
                       << " does not coincide with function of first edge " << endl;
                  
                  DISP(ValPhiEdge(num_loc, j, k));
                  DISP(ValPhiEdge(0, j, k));
                  abort();
                }
            
            for (int k = 0; k < nb_quad; k++)
              {
                Real_wp val2 = ValPhiEdge(num_loc, j, nb_quad-1-k);
                Real_wp val = ValPhiEdge(num_loc, j, k);
                for (int j2 = 0; j2 < nb_dof_edge; j2++)
                  {
                    Real_wp val1 = ValPhiEdge(num_loc, j2, k);
                    B(j, j2) += val2*this->Weights1D(k)*val1;
                    M(j, j2) += val*this->Weights1D(k)*val1;
                  }
              }
          }
        
        GetInverse(M);
        Mlt(B, M, A);
        // DISP(num_loc); DISP(A);
        
        // we check that matrix A is a permutation matrix with sign
        for (int j = 0; j < nb_dof_edge; j++)
          {
            int nb = 0;
            for (int jb = 0; jb < nb_dof_edge; jb++)
              if (abs(A(j, jb)) > threshold)
                nb++;
            
            if (nb != 1)
              {
                abort();
              }
            
            //int j2 = -1; int signe = 0;
            for (int jb = 0; jb < nb_dof_edge; jb++)
              {
                Real_wp val = A(j, jb);
                if (abs(val - 1.0) <= threshold)
                  {
                    //j2 = jb;
                    //signe = 1;
                  }
                else if (abs(val + 1.0) <= threshold)
                  {
                    //j2 = jb;
                    //signe = -1;
                  }
                else if (abs(A(j, jb)) > threshold)
                  {
                    abort();
                  }
              }
            // DISP(num_loc); DISP(j); DISP(j2); DISP(signe);
          }
      }
  }


  void ElementReference<Dimension3, 1>
  ::ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int num_loc) const
  {
    if (this->elt_geom.IsLocalFaceQuadrilateral(num_loc))
      this->element_quad_surf->ComputeProjectionDofRef(feval, contrib);
    else
      this->element_tri_surf->ComputeProjectionDofRef(feval, contrib);

    this->ModifySignProjectionSurface(contrib, num_loc);
  }


  void ElementReference<Dimension3, 1>
  ::ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int num_loc) const
  {
    if (this->elt_geom.IsLocalFaceQuadrilateral(num_loc))
      this->element_quad_surf->ComputeProjectionDofRef(feval, contrib);
    else
      this->element_tri_surf->ComputeProjectionDofRef(feval, contrib);

    this->ModifySignProjectionSurface(contrib, num_loc);
  }
  
  
  //! Integration against gradient of basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
             where omega_k is the weight of integration
             and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \nabla \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = \nabla \varphi_i(\xi_j)
  */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 1>::ApplyRhGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < nb_dof_loc; i++)
      {
	val = Real_wp(0);
	for (int j = 0; j < nb_points_quadrature_inside; j++)
	  val += Gradient_Phi(i, j)(0)*Vh(3*j)
	    + Gradient_Phi(i, j)(1)*Vh(3*j+1) + Gradient_Phi(i, j)(2)*Vh(3*j+2);
	
	Uh(i) = val;
      }
  }
  
  
  //! Integration against derivatives of basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vx Vx_i = \int_K f d/dx(\varphi_i)  dx
    \param[out] Vy Vy_i = \int_K f d/dy(\varphi_i)  dx
    \param[out] Vz Vz_i = \int_K f d/dz(\varphi_i)  dx
    This operation is equivalent to a matrix vector product
    Vx = Rh^1 Uh, Vy = Rh^2 Uh, Vz = Rh^3 Uh
    where (Rh)^1_{i,j} = d/dx(\varphi_i(\xi_j) )
    where (Rh)^2_{i,j} = d/dy(\varphi_i(\xi_j) )
    where (Rh)^3_{i,j} = d/dz(\varphi_i(\xi_j) )
   */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 1>
  ::ApplyRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy, Vector2& Vz) const
  {
    typename Vector2::value_type val_x, val_y, val_z;
    for (int i = 0; i < nb_dof_loc; i++)
      {
        val_x = Real_wp(0);
        val_y = Real_wp(0);
        val_z = Real_wp(0);
        for (int j = 0; j < nb_points_quadrature_inside; j++)
          {
	    val_x += Gradient_Phi(i, j)(0)*Uh(3*j);
	    val_y += Gradient_Phi(i, j)(1)*Uh(3*j+1);
	    val_z += Gradient_Phi(i, j)(2)*Uh(3*j+2);
	  }
	
        Vx(i) = val_x;
        Vy(i) = val_y;
        Vz(i) = val_z;
      }
  }


  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 1>
  ::ApplyRhQuadratureSplitGen(const Vector1& Uh, Vector2& Vh, Vector2&, Vector2&) const
  {
    cout << "not implemented for this finite element" << endl;
    abort();
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
  void ElementReference<Dimension3, 1>::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    TinyVector<typename Vector2::value_type, 3> grad;
   
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        grad.Fill(0);
       
        for (int j = 0; j < nb_dof_loc; j++)
          Add(Uh(j), Gradient_Phi(j, i), grad);
        
        Vh(3*i) = grad(0); Vh(3*i+1) = grad(1); Vh(3*i+2) = grad(2);
      }
  }
  

  //! Integration against basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Ch Vh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 1>::ApplyChGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < nb_dof_loc; i++)
      {
	val = 0;
	for (int j = 0; j < nb_points_quadrature_inside; j++)
	  val += Value_Phi(i, j)*Vh(j);
        
        Uh(i) = val;
      }
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
  void ElementReference<Dimension3, 1>
  ::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	val = 0;
	for (int j = 0; j < nb_dof_loc; j++)
	  val += Value_Phi(j, i)*Uh(j);
	
	Vh(i) = val;
      }
  }


  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 1>
  ::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(const_grad_matrix, Uh, Vh);
  }
  

  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 1>
  ::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, const_grad_matrix, Uh, Vh);
  }


  //! Multiplication by derivative matrices
  /*!
    If we denote (Rh^1)_{i, j} = \int_K \varphi_j d/dx(\varphi_i) dx
    (Rh^2)_{i, j} = \int_K \varphi_j d/dy(\varphi_i) dx
    (Rh^3)_{i, j} = \int_K \varphi_j d/dz(\varphi_i) dx
    then Vx = Rh^1 Uh, Vy = Rh^2 Uh, Vz = Rh^3 Uh
  */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 1>::
  ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy, Vector2& Vz) const
  {
    int m = const_grad_matrix.GetM();
    int n = const_grad_matrix.GetN();
    typename Vector2::value_type val;
    for (int i = 0; i < m; i++)
      {
	val = 0;
	for (int j = 0; j < n; j += 3)
	  val += const_grad_matrix(i, j)*Uh(j);
	
	Vx(i) = val;

	val = 0;
	for (int j = 1; j < n; j += 3)
	  val += const_grad_matrix(i, j)*Uh(j);
	
	Vy(i) = val;

	val = 0;
	for (int j = 2; j < n; j += 3)
	  val += const_grad_matrix(i, j)*Uh(j);
	
	Vz(i) = val;
      }
  }    


  //! we add constant mass matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass coefficient
    \param[inout] A matrix modified
    A(m:, n:) = A(m:, n:) + mass M
    where M is the mass matrix, M_ij = \int \varphi_j \varphi_i dx
   */
  template<class T>
  void ElementReference<Dimension3, 1>::
  AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const
  {
    T zero; SetComplexZero(zero);
    if (mass != zero)
      {
	Vector<int> num(this->nb_dof_loc);
	Vector<T> val(this->nb_dof_loc);
	for (int i = 0; i < this->nb_dof_loc; i++)
	  num(i) = n + i;
	
	for (int i = 0; i < this->nb_dof_loc; i++)
	  {
	    for (int j = 0; j < this->nb_dof_loc; j++)
	      val(j) = mass*mass_matrix(i, j);
	    
	    A.AddInteractionRow(m + i, this->nb_dof_loc, num, val);
	  }
      }
  }


  //! we add constant elementary matrix to the matrix A
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass mass coefficient
    \param[in] C stiffness coefficient
    \param[in] D gradient coefficients
    \param[in] E gradient coefficients
    \param[in] null_term if null_term(i) is true, A, C, D or E is null (i=0, 1, 2, 3)
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + mass M + \sum_p  E(p) (R^p)* + \sum_p  D(p) R^p + \sum_p  C(p, q) S^{p,q}
    where M is the mass matrix equal to
    M_{i, j} = \int_K \phi_j \phi_i dx
    R^p is a gradient matrix equal to :
    (R^p)_{i, j} = \int_K \phi_j d(\phi_i)/dx_p  dx
    S^{p,q} is a stiffness matrix equal to :
    (S^{p,q })_{i, j} = \int_K d(\phi_j)/dx_q d(\phi_i)/dx_p  dx
   */
  template<class T, class Prop>
  void ElementReference<Dimension3, 1>::
  AddConstantElemMatrixGen(int m, int n, const T& mass,
			   const TinyMatrix<T, Prop, 3, 3>& C,
			   const TinyVector<T, 3>& D, 
			   const TinyVector<T, 3>& E,
			   const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;

    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	val.Zero();
	if (!null_term(0))
	  for (int j = 0; j < this->nb_dof_loc; j++)
	    val(j) += mass*mass_matrix(i, j);
	
	if (!null_term(1))
	  for (int j = 0; j < this->nb_dof_loc; j++)
	    val(j) += C(0, 0)*const_stiff_matrix(3*j, 3*i)
	      + C(0, 1)*const_stiff_matrix(3*j+1, 3*i)
	      + C(0, 2)*const_stiff_matrix(3*j+2, 3*i)
	      + C(1, 0)*const_stiff_matrix(3*j, 3*i+1)
	      + C(1, 1)*const_stiff_matrix(3*j+1, 3*i+1)
	      + C(1, 2)*const_stiff_matrix(3*j+2, 3*i+1)
	      + C(2, 0)*const_stiff_matrix(3*j, 3*i+2)
	      + C(2, 1)*const_stiff_matrix(3*j+1, 3*i+2)
	      + C(2, 2)*const_stiff_matrix(3*j+2, 3*i+2);
	
	if (!null_term(2))
	  for (int j = 0; j < this->nb_dof_loc; j++)
	    val(j) += D(0)*const_grad_matrix(i, j*3)
	      + D(1)*const_grad_matrix(i, j*3+1)
	      + D(2)*const_grad_matrix(i, j*3+2);
	
	if (!null_term(3))
	  for (int j = 0; j < this->nb_dof_loc; j++)
	    val(j) += E(0)*const_grad_matrix(j, i*3)
	      + E(1)*const_grad_matrix(j, i*3+1)
	      + E(2)*const_grad_matrix(j, i*3+2);
	
	A.AddInteractionRow(m + i, this->nb_dof_loc, num, val);
      }
  }
   
  
  //! adds \int C \nabla varphi_j \nabla varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] mass mass coefficient
    \param[in] C stiffness coefficient
    \param[in] D gradient coefficients
    \param[in] E gradient coefficients
    \param[in] null_term if null_term(i) is true, A, C, D or E is null (i=0, 1, 2, 3)
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + M + R + Rt + S
    M_ij = \int mass \varphi_j \varphi_i
    R_ij = \int D \cdot \nabla \varphi_i \varphi_j
    Rt_ij = \int E \cdot \nabla \varphi_j \varphi_i
    where S_ij = \int C \nabla varphi_j \nabla varphi_i
   */
  template<class T, class Prop>
  void ElementReference<Dimension3, 1>
  ::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& mass,
			     const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
			     const Vector<TinyVector<T, 3> >& D,
			     const Vector<TinyVector<T, 3> >& E, 
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    T vloc;
    TinyVector<T, 3> vec_u;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
	    SetComplexZero(vloc);
	    if (!null_term(0))
	      for (int k = 0; k < this->nb_points_quadrature_inside; k++)
		vloc += this->Value_Phi(i, k)*this->Value_Phi(j, k)*mass(k);
	    
	    if (!null_term(1))
	      for (int k = 0; k < this->nb_points_quadrature_inside; k++)
		{
		  Mlt(C(k), this->Gradient_Phi(j, k), vec_u);
		  vloc += DotProd(vec_u, this->Gradient_Phi(i, k));
		}
	    
	    if (!null_term(2))
	      for (int k = 0; k < this->nb_points_quadrature_inside; k++)
		vloc += this->Value_Phi(j, k)*DotProd(D(k), this->Gradient_Phi(i, k));

	    if (!null_term(3))
	      for (int k = 0; k < this->nb_points_quadrature_inside; k++)
		vloc += this->Value_Phi(i, k)*DotProd(E(k), this->Gradient_Phi(j, k));
	    
	    val(j) = vloc;
	  }
	
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }


  //! adds \int A varphi_j \varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A coefficient A on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int A varphi_j  varphi_i
   */
  template<class T>
  void ElementReference<Dimension3, 1>
  ::AddVariableMassMatrixGen(int off_row, int off_col,
			     const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	val.Zero();
	for (int j = 0; j < this->nb_dof_loc; j++)
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	    val(j) += this->Value_Phi(i, k)*this->Value_Phi(j, k)*A(k);
	
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }


  //! adds \int C \nabla varphi_j \nabla varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A mass coefficients
    \param[in] C tensor C on each quadrature point
    \param[in] D gradient coefficients
    \param[in] E gradient coefficients
    \param[in] null_term if null_term(i) is true, A, C, D or E is null (i=0, 1, 2, 3)
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + M + R + Rt + S
    M_ij = \int mass \varphi_j \varphi_i
    R_ij = \int D \cdot \nabla \varphi_i \varphi_j
    Rt_ij = \int E \cdot \nabla \varphi_j \varphi_i
    where S_ij = \int C \nabla varphi_j \nabla varphi_i
    This function uses function ApplyCh, ApplyRh, considering that
    they are optimized
   */
  template<class T, class Prop>
  void ElementReference<Dimension3, 1>::
  AddVariableElemMatrixOpt(int off_row, int off_col, const Vector<T>& A,
			   const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
			   const Vector<TinyVector<T, 3> >& D,
			   const Vector<TinyVector<T, 3> >& E, 
			   const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    int Nquad = this->nb_points_quadrature_inside;
    Vector<Real_wp> Ones(this->nb_dof_loc), grad_phi(3*Nquad), val_phi(Nquad);
    Vector<T> contrib(this->nb_dof_loc), feval_grad(3*Nquad), feval(Nquad);
    R3 vec_u; TinyVector<T, 3> vec_v;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
        Ones.Zero(); SetComplexOne(Ones(i));
	this->ApplyChTranspose(Ones, val_phi);
        this->ApplyRhTranspose(Ones, grad_phi);

	feval.Zero();
        for (int k = 0; k < Nquad; k++)
          {
            vec_u(0) = grad_phi(3*k);
            vec_u(1) = grad_phi(3*k+1);
            vec_u(2) = grad_phi(3*k+2);
	    
	    if (!null_term(0))
	      feval(k) = A(k)*val_phi(k);
	    
	    if (!null_term(1))
	      MltTrans(C(k), vec_u, vec_v);
	    else
	      vec_v.Zero();
	    
	    if (!null_term(2))
	      feval(k) += DotProd(D(k), vec_u);
	    
	    if (!null_term(3))
	      {
		vec_v(0) += E(k)(0)*val_phi(k);
		vec_v(1) += E(k)(1)*val_phi(k);
		vec_v(2) += E(k)(2)*val_phi(k);
	      }

            feval_grad(3*k) = vec_v(0);
            feval_grad(3*k+1) = vec_v(1);
            feval_grad(3*k+2) = vec_v(2);
          }
	
	this->ApplyCh(feval, val);	
        this->ApplyRh(feval_grad, contrib);
        for (int j = 0; j < this->nb_dof_loc; j++)
	  val(j) += contrib(j);
	
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }
        

  //! adds \int A varphi_j \varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A coefficient A on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int A varphi_j  varphi_i
    This function uses function ApplyCh, ApplyRh, considering that
    they are optimized
   */
  template<class T>
  void ElementReference<Dimension3, 1>::
  AddVariableMassMatrixOpt(int off_row, int off_col,
                           const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    int Nquad = this->nb_points_quadrature_inside;
    Vector<Real_wp> val_phi(Nquad);
    Vector<T> feval(Nquad);
    Vector<Real_wp> Ones(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
        Ones.Zero(); Ones(i) = 1.0;
        this->ApplyChTranspose(Ones, val_phi);
        for (int k = 0; k < Nquad; k++)
          feval(k) = val_phi(k)*A(k);
        
        this->ApplyCh(feval, val);
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }


  //! computes integrals against gradient of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension3, 1>
  ::ComputeIntegralGradientRef(const Vector1 & feval, Vector2& res) const
  {
    this->ApplyRh(feval, res);
  }
  
  
  //! computes surface integrals against gradient of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension3, 1>
  ::ComputeIntegralSurfaceGradientRef(const Vector1 & feval, Vector2& res, int n) const
  {
    typename Vector2::value_type one; SetComplexOne(one);
    res.Fill(0);
    this->ApplyNablaSh(one, n, feval, res);
  }
#endif

    
  /********************
   * H(curl) elements *
   ********************/


#ifdef MONTJOIE_WITH_NODAL_HCURL
  //! Constructor
  ElementReference<Dimension3, 2>::ElementReference(ElementGeomReference<Dimension3>& elt)
    : VolumeReference<2>(elt)
  {
    element_tri_surf = NULL; element_quad_surf = NULL;
    element_scal_vol = NULL;
  }
  

  //! Destructor
  ElementReference<Dimension3, 2>::~ElementReference()
  {
    if (element_tri_surf != NULL)
      {
	delete element_tri_surf;
	element_tri_surf = NULL;
      }

    if (element_quad_surf != NULL)
      {
	delete element_quad_surf;
	element_quad_surf = NULL;
      }
  }
  


  //! returns 2-D triangular finite element (tangential traces on triangular faces)
  const ElementReference<Dimension2, 2>& 
  ElementReference<Dimension3, 2>::GetTriangularSurfaceFiniteElement() const
  {
    if (element_tri_surf == NULL)
      {
	cout << "Null pointer for GetTriangularSurfaceFiniteElement" << endl;
	abort();
      }

    return *element_tri_surf;
  }
  
  
  //! returns 2-D quadrangular finite element (tangential traces on quadrangular faces)
  const ElementReference<Dimension2, 2>& 
  ElementReference<Dimension3, 2>::GetQuadrangularSurfaceFiniteElement() const
  {
    if (element_quad_surf == NULL)
      {
	cout << "Null pointer for GetQuadrangularSurfaceFiniteElement" << endl;
	abort();
      }

    return *element_quad_surf;
  }


  const ElementReference<Dimension2, 2>& 
  ElementReference<Dimension3, 2>::GetSurfaceFiniteElement(int n) const
  {
    if (this->elt_geom.IsLocalFaceQuadrilateral(n))
      return GetQuadrangularSurfaceFiniteElement();
    else
      return GetTriangularSurfaceFiniteElement();
  }


  const ElementReference<Dimension3, 1>& 
  ElementReference<Dimension3, 2>::GetScalarElement() const
  {
    if (element_scal_vol == NULL)
      {
	cout << "Null pointer for GetScalarElement" << endl;
	abort();
      }

    return *element_scal_vol;
  }


  size_t ElementReference<Dimension3, 2>::GetMemorySize() const
  {
    size_t taille = VolumeReference<2>::GetMemorySize() + sizeof(FiniteElementHcurl<Dimension3>);
    if (element_tri_surf != NULL)
      taille += element_tri_surf->GetMemorySize();

    if (element_quad_surf != NULL)
      taille += element_quad_surf->GetMemorySize();
    
    return taille;
  }
  

  void ElementReference<Dimension3, 2>
  ::ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int num_loc) const
  {
    int Npts = feval.GetM()/3;
    VectReal_wp feval_surf(2*Npts);
    TinyVector<Real_wp, 3> vec_u;
    for (int i = 0; i < Npts; i++)
      {
	CopyVector(feval, i, vec_u);
	
	feval_surf(2*i) = DotProd(vec_u, this->TangenteLocX(num_loc));
	feval_surf(2*i+1) = DotProd(vec_u, this->TangenteLocY(num_loc));
      }

    contrib.Zero();
    if (this->elt_geom.IsLocalFaceQuadrilateral(num_loc))
      this->element_quad_surf->ComputeProjectionDofRef(feval_surf, contrib);
    else
      this->element_tri_surf->ComputeProjectionDofRef(feval_surf, contrib);
    
    this->ModifySignProjectionSurface(contrib, num_loc);
  }

  
  void ElementReference<Dimension3, 2>
  ::ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int num_loc) const
  {
    int Npts = feval.GetM()/3;
    VectComplex_wp feval_surf(2*Npts);
    TinyVector<Complex_wp, 3> vec_u;
    for (int i = 0; i < Npts; i++)
      {
	CopyVector(feval, i, vec_u);
	
	feval_surf(2*i) = DotProd(vec_u, this->TangenteLocX(num_loc));
	feval_surf(2*i+1) = DotProd(vec_u, this->TangenteLocY(num_loc));
      }

    if (this->elt_geom.IsLocalFaceQuadrilateral(num_loc))
      this->element_quad_surf->ComputeProjectionDofRef(feval_surf, contrib);
    else
      this->element_tri_surf->ComputeProjectionDofRef(feval_surf, contrib);

    this->ModifySignProjectionSurface(contrib, num_loc);
  }


  void ElementReference<Dimension3, 2>
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {
    if (this->elt_geom.IsLocalFaceQuadrilateral(num_loc))
      {
	for (int k = 2*this->order; k < 4*this->order; k++)
	  contrib(k) = -contrib(k);
      }
    else
      {
	for (int k = 2*this->order; k < 3*this->order; k++)
	  contrib(k) = -contrib(k);
      }
  }


  void ElementReference<Dimension3, 2>
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    if (this->elt_geom.IsLocalFaceQuadrilateral(num_loc))
      {
	for (int k = 2*this->order; k < 4*this->order; k++)
	  contrib(k) = -contrib(k);
      }
    else
      {
	for (int k = 2*this->order; k < 3*this->order; k++)
	  contrib(k) = -contrib(k);
      }
  }


  //! computation of mass matrix for edge elements
  void ElementReference<Dimension3, 2>::ConstructHcurlElementaryMatrix()
  {
    int nb_points_quadrature = this->GetNbPointsQuadrature();
    VectR3 phi, curl_phi;
    // evaluation of basis functions on quadrature points
    this->Value_PhiVec.Reallocate(this->nb_dof_loc, nb_points_quadrature);        
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	this->ComputeValuesPhiRef(this->PointsND(i), phi);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  this->Value_PhiVec(j, i) = phi(j);
      }
    
    // mass matrix on reference element
    this->mass_matrix.Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    Real_wp vloc;
    
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = i; j < this->nb_dof_loc; j++)
	{
	  vloc = 0.0;
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*DotProd(this->Value_PhiVec(i, k), this->Value_PhiVec(j, k));
          
	  this->mass_matrix(i,j) = vloc;
	}
    
    //VectReal_wp lambda; 
    //GetEigenvalues(this->mass_matrix, lambda);
    //DISP(lambda);
	  
    this->mass_matrix_chol = this->mass_matrix;
    GetCholesky(this->mass_matrix_chol);
    
    // projection on dofs
    if (this->tangente_dof.GetM() > 0)
      {
	this->MatProjectionDof.Reallocate(this->GetNbPointsDofInside(), this->nb_dof_loc);
	for (int i = 0; i < this->GetNbPointsDofInside(); i++)
	  {
	    this->ComputeValuesPhiRef(this->PointsDofND(i), phi);
	    for (int j = 0; j < this->nb_dof_loc; j++)
	      this->MatProjectionDof(i, j) = DotProd(phi(j), this->tangente_dof(i));
	  }
	
	GetQR(this->MatProjectionDof, this->tauProjectionDof);
      }
    
    // values of basis functions on nodal points
    this->ValuePhi_Nodal.Reallocate(this->nb_dof_loc, 3*this->GetNbPointsNodalElt());
    for (int i = 0; i < this->GetNbPointsNodalElt(); i++)
      {
	this->ComputeValuesPhiRef(this->PointsNodalND(i), phi);
	
        for (int j = 0; j < this->nb_dof_loc; j++)
          {
            this->ValuePhi_Nodal(j, 3*i) = phi(j)(0);
            this->ValuePhi_Nodal(j, 3*i+1) = phi(j)(1);
            this->ValuePhi_Nodal(j, 3*i+2) = phi(j)(2);
          }
      }
    
    // curl of basis functions on nodal points
    this->CurlPhi_Nodal.Reallocate(this->nb_dof_loc, 3*this->GetNbPointsNodalElt());
    for (int i = 0; i < this->GetNbPointsNodalElt(); i++)
      {
        this->ComputeCurlPhiRef(this->PointsNodalND(i), phi);
        for (int j = 0; j < this->nb_dof_loc; j++)
          {
            this->CurlPhi_Nodal(j, 3*i) = phi(j)(0);
            this->CurlPhi_Nodal(j, 3*i+1) = phi(j)(1);
            this->CurlPhi_Nodal(j, 3*i+2) = phi(j)(2);
          }
      }
    
    // curl of basis on quadrature points
    this->Curl_Phi.Reallocate(this->nb_dof_loc, nb_points_quadrature);
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	this->ComputeCurlPhiRef(this->PointsND(i), curl_phi);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  this->Curl_Phi(j, i) = curl_phi(j);
      }    
    
    // mass and stiffness matrix for constant coefficients
    // \int_K C phi_i phi_j   and   \int_K C curl(phi_i)  curl(phi_j)
    // where C is a 3x3 symmetric matrix
    for (int p = 0; p < 3; p++)
      for (int q = 0; q < 3; q++)
        {
          this->const_mass_matrix(p, q).Reallocate(this->nb_dof_loc, this->nb_dof_loc);
          this->const_curl_matrix(p, q).Reallocate(this->nb_dof_loc, this->nb_dof_loc);
        }
    
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = 0; j < this->nb_dof_loc; j++)
        {
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i, k)(0)*this->Value_PhiVec(j, k)(0);
          
          this->const_mass_matrix(0, 0)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(0)*this->Value_PhiVec(j,k)(1);
          
          this->const_mass_matrix(0, 1)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(1)*this->Value_PhiVec(j,k)(0);
          
          this->const_mass_matrix(1, 0)(i, j) = vloc;

          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(0)*this->Value_PhiVec(j,k)(2);
             
              this->const_mass_matrix(0, 2)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc +=  this->WeightsND(k)*this->Value_PhiVec(i,k)(2)*this->Value_PhiVec(j,k)(0);
          
          this->const_mass_matrix(2, 0)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(1)*this->Value_PhiVec(j,k)(1);
          
          this->const_mass_matrix(1, 1)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(1)*this->Value_PhiVec(j,k)(2);
              
          this->const_mass_matrix(1, 2)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(2)*this->Value_PhiVec(j,k)(1);
          
          this->const_mass_matrix(2, 1)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(2)*this->Value_PhiVec(j,k)(2);
          
          this->const_mass_matrix(2, 2)(i, j) = vloc;

          // curl curl matrix
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Curl_Phi(i,k)(0)*this->Curl_Phi(j,k)(0);
          
          this->const_curl_matrix(0, 0)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Curl_Phi(i,k)(0)*this->Curl_Phi(j,k)(1);
          
          this->const_curl_matrix(0, 1)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Curl_Phi(i,k)(1)*this->Curl_Phi(j,k)(0);
          
          this->const_curl_matrix(1, 0)(i, j) = vloc;

          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Curl_Phi(i,k)(0)*this->Curl_Phi(j,k)(2);
          
          this->const_curl_matrix(0, 2)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Curl_Phi(i,k)(2)*this->Curl_Phi(j,k)(0);
          
          this->const_curl_matrix(2, 0)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Curl_Phi(i,k)(1)*this->Curl_Phi(j,k)(1);
          
          this->const_curl_matrix(1, 1)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Curl_Phi(i,k)(1)*this->Curl_Phi(j,k)(2);
          
          this->const_curl_matrix(1, 2)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Curl_Phi(i,k)(2)*this->Curl_Phi(j,k)(1);
          
          this->const_curl_matrix(2, 1)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Curl_Phi(i,k)(2)*this->Curl_Phi(j,k)(2);
          
          this->const_curl_matrix(2, 2)(i, j) = vloc;
        }

    Matrix<Real_wp, General, ArrayRowSparse> Sh;
    int offset = this->nb_points_quadrature_inside;
    for (int n = 0; n < this->GetNbBoundaries(); n++)
      {
        int Nquad = this->GetNbQuadBoundary(n);
        Sh.Clear();
        Sh.Reallocate(this->nb_dof_loc, 3*Nquad);
        for (int i = 0; i < this->nb_dof_loc; i++)
          for (int k = 0; k < Nquad; k++)
            {
              if (abs(this->Value_PhiVec(i, offset + k)(0)) > 20.0*epsilon_machine)
                Sh.AddInteraction(i, 3*k, this->Value_PhiVec(i, offset + k)(0));

              if (abs(this->Value_PhiVec(i, offset + k)(1)) > 20.0*epsilon_machine)
                Sh.AddInteraction(i, 3*k+1, this->Value_PhiVec(i, offset + k)(1));

              if (abs(this->Value_PhiVec(i, offset + k)(2)) > 20.0*epsilon_machine)
                Sh.AddInteraction(i, 3*k+2, this->Value_PhiVec(i, offset + k)(2));
            }
        
        Copy(Sh, this->sparse_const_sh(n));

        Sh.Clear();
        Sh.Reallocate(this->nb_dof_loc, 3*Nquad);
        for (int i = 0; i < this->nb_dof_loc; i++)
          for (int k = 0; k < Nquad; k++)
            {
              if (abs(this->Curl_Phi(i, offset + k)(0)) > 20.0*epsilon_machine)
                Sh.AddInteraction(i, 3*k, this->Curl_Phi(i, offset + k)(0));

              if (abs(this->Curl_Phi(i, offset + k)(1)) > 20.0*epsilon_machine)
                Sh.AddInteraction(i, 3*k+1, this->Curl_Phi(i, offset + k)(1));

              if (abs(this->Curl_Phi(i, offset + k)(2)) > 20.0*epsilon_machine)
                Sh.AddInteraction(i, 3*k+2, this->Curl_Phi(i, offset + k)(2));
            }
        
        Copy(Sh, this->const_nabla_sh(n));

        offset += Nquad;
      }
  }
  
  
  //! filling FacesDof and checking that dofs are numbered correctly 
  //! (first dof on edges, then faces, then interior)
  void ElementReference<Dimension3, 2>::FindDofsOnFace(bool nodal_element)
  {
    int N = this->GetNbDof();
    if (N <= 0)
      {
        cout << "no degree of freedom "<<endl;
        abort();
      }
    
    Real_wp threshold = 1e4*epsilon_machine;
    // for each face we are computing phi \times n 
    Matrix<bool> tangential_proj_null(N, this->GetNbBoundaries());
    tangential_proj_null.Fill(true);
    VectR3 phi; R3 phi_n;
    for (int num_loc = 0; num_loc < this->GetNbBoundaries(); num_loc++)
      for (int i = 0; i < this->GetNbQuadBoundary(num_loc); i++)
        {
          int k = this->GetQuadNumber(num_loc, i);
          this->GetValuePhiOnQuadraturePoint(k, phi);
          for (int j = 0; j < N; j++)
            {
              TimesProd(phi(j), this->NormaleLoc(num_loc), phi_n);
              if (Norm2(phi_n) > threshold)
                {
                  tangential_proj_null(j, num_loc) = false;
                }
            }
        }
    
    // finding dofs inside element
    Vector<bool> DofInside(N); DofInside.Fill(true);
    int nb = 0;
    for (int i = 0; i < this->GetNbDof(); i++)
      {
        for (int j = 0; j < this->GetNbBoundaries(); j++)
          if (!tangential_proj_null(i, j))
            {
              DofInside(i) = false;
              // cout << "tangential component of dof " << i << " on face " 
              // << j << " is different from 0 " << endl;     
            }
        
        if (DofInside(i))
          nb++;
      }
    
    for (int i = N-nb; i < N; i++)
      if (!DofInside(i))
        {
          cout << "Internal dofs must be numbered at the end " << endl;
          abort();
        }
    
    if (this->nb_dof_boundaries <= 0)
      this->nb_dof_boundaries = N - nb;
    else
      {
        if (this->nb_dof_boundaries != N - nb)
          {
            cout << "Number of dofs on the boundary is noticed as " 
                 << this->nb_dof_boundaries << endl;
            cout << "But we found " << N-nb << " dofs on the boundary " << endl;
            abort();
          }
      }
    
    // finding dofs inside faces
    int nb_max = this->nb_dof_boundaries;
    IVect nb_dof_face(this->GetNbBoundaries()); nb_dof_face.Fill(0);
    for (int num_loc = this->GetNbBoundaries()-1; num_loc >= 0; num_loc--)
      {
        DofInside.Fill(true);
        nb = 0;
        for (int i = 0; i < nb_max; i++)
          {
            for (int num_loc2 = 0; num_loc2 < this->GetNbBoundaries(); num_loc2++)
              if (num_loc2 != num_loc)
                if (!tangential_proj_null(i, num_loc2))
                  DofInside(i) = false;
            
            if (DofInside(i))
              nb++;
          }
        
        for (int i = nb_max-nb; i < nb_max; i++)
          if (!DofInside(i))
            {
              cout << "dof number for this face is incorrect " << endl;
              abort();
            }       
        
        nb_dof_face(num_loc) = nb;
        nb_max -= nb;
      }
    
    // dofs on edges
    IVect nb_dof_edge(this->GetNbEdges()); nb_dof_edge.Fill(0);
    int type_elt = this->GetHybridType();
    for (int num_edge = this->GetNbEdges()-1; num_edge >= 0; num_edge--)
      {
        DofInside.Fill(true);
        nb = 0;
        for (int i = 0; i < nb_max; i++)
          {
            for (int j = 0; j < this->GetNbBoundaries(); j++)
              {
                bool face_adj = false;
                for (int k = 0; k < 4; k++)
                  if (MeshNumbering<Dimension3>::GetEdgeNumberOfFace(type_elt, j, k) == num_edge)
                    face_adj = true;
                
                if (!face_adj)
                  if (!tangential_proj_null(i, j))
                    DofInside(i) = false;
              }
            
            if (DofInside(i))
              nb++;
          }
        
        for (int i = nb_max-nb; i < nb_max; i++)
          if (!DofInside(i))
            {
              cout << "dof number for this edge is incorrect " << endl;
              abort();
            }       
        
        nb_dof_edge(num_edge) = nb;
        nb_max -= nb;
      }
    
    // filling FacesDof
    nb_max = 0;
    for (int j = 0; j < this->GetNbBoundaries(); j++)
      {
        nb = nb_dof_face(j);
        for (int k = 0; k < 4; k++)
          {
            int ne = MeshNumbering<Dimension3>::GetEdgeNumberOfFace(type_elt, j, k);
            if (ne >= 0)
              nb += nb_dof_edge(ne);
          }
        
        nb_max = max(nb_max, nb);
      }
    
    // eventually filling FacesDof
    FacesDof.Reallocate(nb_max, this->GetNbBoundaries());
    FacesDof.Fill(-1);
    int nb_dof = nb_dof_edge(0);
    int offset = nb_dof*this->GetNbEdges();
    for (int j = 0; j < this->GetNbBoundaries(); j++)
      {
        nb = 0;
        for (int k = 0; k < 4; k++)
          {
            int ne = MeshNumbering<Dimension3>::GetEdgeNumberOfFace(type_elt, j, k);
            if (ne >= 0)
              {
		int n1 = MeshNumbering<Dimension3>::GetVertexNumberOfFace(type_elt, j, k);
		int n1b = MeshNumbering<Dimension3>::FirstExtremityEdge(type_elt, ne);
		if (!nodal_element)
		  n1b = n1;
		
		if (n1 == n1b)
		  {
		    // edges in the same direction
		    for (int i = 0; i < nb_dof; i++)
		      FacesDof(nb++, j) = ne*nb_dof + i;
		  }
		else
		  {
		    // edges in the opposite direction
		    for (int i = 0; i < nb_dof; i++)
		      FacesDof(nb++, j) = ne*nb_dof + nb_dof-1-i;
		  }
              }
          }
        
        for (int i = 0; i < nb_dof_face(j); i++)
          FacesDof(nb++, j) = offset+i;
        
        offset += nb_dof_face(j);
      }

    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);
  }


  //! computing operator to apply when the faces are rotated (for edge elements)
  void ElementReference<Dimension3, 2>::
  FindHcurlLinearCombinationRotation(NumberMap& nmap,
                                     const ElementReference<Dimension2, 2>& Fb_tri,
				     const ElementReference<Dimension2, 2>& Fb_quad,
                                     bool take_all_dofs, bool display_message) const
  {
    // vertex numbers of triangular face after rotation
    TinyVector<TinyVector<int, 3>, 6> rot_tri_num;
    rot_tri_num(0).Init(0, 1, 2);
    rot_tri_num(1).Init(1, 2, 0);
    rot_tri_num(2).Init(2, 0, 1);
    rot_tri_num(3).Init(0, 2, 1);
    rot_tri_num(4).Init(1, 0, 2);
    rot_tri_num(5).Init(2, 1, 0);
    
    // vertex numbers of quadrangular face after rotation
    TinyVector<TinyVector<int, 4>, 8> rot_quad_num;
    rot_quad_num(0).Init(0, 1, 2, 3);
    rot_quad_num(1).Init(1, 2, 3, 0);
    rot_quad_num(2).Init(2, 3, 0, 1);
    rot_quad_num(3).Init(3, 0, 1, 2);
    rot_quad_num(4).Init(0, 3, 2, 1);
    rot_quad_num(5).Init(1, 0, 3, 2);
    rot_quad_num(6).Init(2, 1, 0, 3);
    rot_quad_num(7).Init(3, 2, 1, 0);
    
    Matrix<int> IndexDofFace(this->GetNbDof(), this->GetNbBoundaries());
    IndexDofFace.Fill(-1);
    bool first_tri = true, first_quad = true;
    int type_elt = this->GetHybridType();
    Real_wp threshold = 1e5*epsilon_machine;
    for (int num_loc = 0; num_loc < this->GetNbBoundaries(); num_loc++)
      {
        bool quad = MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(type_elt, num_loc);
        bool check_face = false;
        if (quad)
          {
            if (first_tri)
              {
                first_tri = false;
                check_face = true;
              }
          }
        else
          {
            if (first_quad)
              {
                first_quad = false;
                check_face = true;
              }
          }
        
        if (check_face)
          {
            Matrix<int> FaceRot;
            Matrix<bool> OppositeSigneDof;
            bool linear_combination = false;
            
            // finding internal dofs of the face
            Vector<bool> dof_on_face(this->GetNbDof());
	    int N = this->GetNbQuadBoundary(num_loc);
            VectR3 phi, phi1, phi2; R3 phi_n, phi_t;
	    if (take_all_dofs)
	      {
		dof_on_face.Fill(false);
		for (int i = 0; i < this->GetNbDofBoundaries(); i++)
		  if (this->IsTangentialDof(i, num_loc))
		    dof_on_face(i) = true;
	      }
	    else
	      {
		dof_on_face.Fill(true);
		for (int i = this->GetNbDofBoundaries(); i < this->GetNbDof(); i++)
		  dof_on_face(i) = false;
		
		for (int num_loc2 = 0; num_loc2 < this->GetNbBoundaries(); num_loc2++)
		  if (num_loc2 != num_loc)
		    for (int i = 0; i < this->GetNbQuadBoundary(num_loc2); i++)
		      {
			int k = this->GetQuadNumber(num_loc2, i);
			this->GetValuePhiOnQuadraturePoint(k, phi);
			
			for (int j = 0; j < this->GetNbDof(); j++)
			  { 
			    TimesProd(phi(j), this->NormaleLoc(num_loc2), phi_n);
			    if (Norm2(phi_n) > threshold)
			      {
				dof_on_face(j) = false;
			      }
			  }
		      }           
	      }
	    
	    //DISP(dof_on_face);
            
            IVect IndexDof(this->GetNbDof()); IndexDof.Fill(-1);
            int nb_inside = 0;
            for (int j = 0; j < this->GetNbDof(); j++)
              if (dof_on_face(j))
                {
                  IndexDof(j) = nb_inside;
                  IndexDofFace(j, num_loc) = nb_inside;
                  nb_inside++;
                }
	    //DISP(num_loc); DISP(nb_inside);
            
            int nb_orientations = 6; int nb_vert_on_face = 3;
            if (quad)
              {
                nb_orientations = 8;
                nb_vert_on_face = 4;
              }
            
            if (nb_inside > 0)
              {    
                FaceRot.Reallocate(nb_orientations, nb_inside); FaceRot.Fill(-1); 
                OppositeSigneDof.Reallocate(nb_orientations, nb_inside);
                OppositeSigneDof.Fill(false);
                Vector<Matrix<Real_wp> > CoefCombination(nb_orientations);
                
                // small mesh with one element
                Mesh<Dimension3> mesh;
                mesh.ReallocateVertices(this->GetNbVertices());
                
                switch (type_elt)
                  {
                  case 0 :
                    {
                      mesh.Vertex(0).Init(0, 0, 0);
                      mesh.Vertex(1).Init(1, 0, 0);
                      mesh.Vertex(2).Init(0, 1, 0);
                      mesh.Vertex(3).Init(0, 0, 1);
                    }
                    break;
                  case 1 :
                    {
                      mesh.Vertex(0).Init(-1, -1, 0);
                      mesh.Vertex(1).Init(1, -1, 0);
                      mesh.Vertex(2).Init(1, 1, 0);
                      mesh.Vertex(3).Init(-1, 1, 0);
                      mesh.Vertex(4).Init(0, 0, 1);
                    }
                    break;
                  case 2 :
                    {
                      mesh.Vertex(0).Init(0, 0, 0);
                      mesh.Vertex(1).Init(1, 0, 0);
                      mesh.Vertex(2).Init(0, 1, 0);
                      mesh.Vertex(3).Init(0, 0, 1);
                      mesh.Vertex(4).Init(1, 0, 1);
                      mesh.Vertex(5).Init(0, 1, 1);
                    }
                    break;
                  case 3 :
                    {
                      mesh.Vertex(0).Init(0, 0, 0);
                      mesh.Vertex(1).Init(1, 0, 0);
                      mesh.Vertex(2).Init(1, 1, 0);
                      mesh.Vertex(3).Init(0, 1, 0);
                      mesh.Vertex(4).Init(0, 0, 1);
                      mesh.Vertex(5).Init(1, 0, 1);
                      mesh.Vertex(6).Init(1, 1, 1);
                      mesh.Vertex(7).Init(0, 1, 1);
                    }
                  }
                
                mesh.ReallocateElements(1);
                int ref = 1; IVect num(this->GetNbVertices());
                num.Fill();
                mesh.Element(0).Init(num, ref);
                mesh.FindConnectivity();
                
                VectR3 s, sFace; Matrix3_3 DF1, invDF1, DF2, invDF2;
                R3 pt, ptGlob1, ptGlob2, normale_fj;
                TinyVector<R3, 3> vec_DF;
                VectR2 phi_2d;
                VectR3 phi_glob1(this->GetNbDof()), phi_glob2;
                // we test all orientations
                for (int rot = 0; rot < nb_orientations; rot++)
                  {
                    IVect num(3);
                    if (quad)
                      num.Reallocate(4);
                    
                    sFace.Reallocate(num.GetM());
                    for (int j = 0; j < nb_vert_on_face; j++)
                      {
                        int jrot = -1;
                        if (quad)
                          jrot = rot_quad_num(rot)(j);
                        else
                          jrot = rot_tri_num(rot)(j);
                        
                        int nv = MeshNumbering<Dimension3>::
                          GetVertexNumberOfFace(type_elt, num_loc, j);
                        
                        num(jrot) = nv;
                        sFace(jrot) = mesh.Vertex(nv);
                      }
                    
                    // on calcule les DF de l'element 1 et 2
                    mesh.GetVerticesElement(0, s);
                    
                    // on boucle sur les points de quadrature de la face            
                    Matrix<Real_wp> Mh(nb_inside, nb_inside), Bh(nb_inside, nb_inside);
                    Mh.Fill(0); Bh.Fill(0); Real_wp x, y;
                    for (int ks = 0; ks < N; ks++)
                      {
                        Real_wp poids = this->WeightsQuadratureBoundary(ks, num_loc);
                        
                        // treating 3-D face of the element
                        int k = this->GetQuadNumber(num_loc, ks);
                        pt = this->PointsND(k);
                        this->FjLinear(s, pt, ptGlob1);
                        
                        this->DFjLinear(s, pt, DF1);
                        GetInverse(DF1, invDF1);
                        
                        Real_wp dsj;
                        this->GetNormale(invDF1, normale_fj, dsj, num_loc);
                        this->GetValuePhiOnQuadraturePoint(k, phi1);
                        
                        // treating 3-D face from a reference triangle/quadrangle
                        int k2 = -1;
                        if (quad)
                          {
                            // computation of DF2^*-1
                            vec_DF(0) = sFace(1) - sFace(0);
                            vec_DF(1) = sFace(3) - sFace(0);
                            TimesProd(vec_DF(0), vec_DF(1), vec_DF(2));
                            DF2 = Matrix3_3(vec_DF);
                            GetInverse(DF2, invDF2);
                            
                            // searching point on 3-D face
                            for (int i2 = 0; i2 < Fb_quad.GetNbPointsQuadratureInside(); i2++)
                              {
                                x = Fb_quad.PointsND(i2)(0);
                                y = Fb_quad.PointsND(i2)(1);
                                ptGlob2 = (1.0-x)*(1.0-y)*sFace(0) + x*(1.0-y)*sFace(1)
                                  + x*y*sFace(2) + (1.0-x)*y*sFace(3);
                                
                                if (ptGlob2.Distance(ptGlob1) < threshold)
                                  k2 = i2;                            
                              }
                            
                            Fb_quad.GetValuePhiOnQuadraturePoint(k2, phi_2d);
                            phi2.Reallocate(Fb_quad.GetNbDof());
                            for (int j = 0; j < Fb_quad.GetNbDof(); j++)
                              phi2(j).Init(phi_2d(j)(0), phi_2d(j)(1), 0);
                            
                            phi_glob2.Reallocate(Fb_quad.GetNbDof());
                            for (int j = 0; j < Fb_quad.GetNbDof(); j++)
                              {
                                MltTrans(invDF2, phi2(j), phi_t);
                                TimesProd(normale_fj, phi_t, phi_glob2(j));
                              }
                          }
                        else
                          {
                            // computation of DF2^*-1
                            vec_DF(0) = sFace(1) - sFace(0);
                            vec_DF(1) = sFace(2) - sFace(0);
                            TimesProd(vec_DF(0), vec_DF(1), vec_DF(2));
                            DF2 = Matrix3_3(vec_DF);
                            GetInverse(DF2, invDF2);
                            
                            // searching point on 3-D face
                            for (int i2 = 0; i2 < Fb_tri.GetNbPointsQuadratureInside(); i2++)
                              {
                                x = Fb_tri.PointsND(i2)(0);
                                y = Fb_tri.PointsND(i2)(1);
                                ptGlob2 = (1.0-x-y)*sFace(0) + x*sFace(1) + y*sFace(2);
                                if (ptGlob2.Distance(ptGlob1) < threshold)
                                  k2 = i2;                            
                              }
                            
                            Fb_tri.GetValuePhiOnQuadraturePoint(k2, phi_2d);
                            phi2.Reallocate(Fb_tri.GetNbDof());
                            for (int j = 0; j < Fb_tri.GetNbDof(); j++)
                              phi2(j).Init(phi_2d(j)(0), phi_2d(j)(1), 0);
                            
                            phi_glob2.Reallocate(Fb_tri.GetNbDof());
                            for (int j = 0; j < Fb_tri.GetNbDof(); j++)
                              {
                                MltTrans(invDF2, phi2(j), phi_t);
                                TimesProd(normale_fj, phi_t, phi_glob2(j));
                              }
                          }
                        
                        // phi globaux
                        for (int j = 0; j < this->GetNbDof(); j++)
                          {
                            MltTrans(invDF1, phi1(j), phi_t);
                            TimesProd(normale_fj, phi_t, phi_glob1(j));                       
                          }
                        
                        // mise a jour matrice M et B
                        for (int iv = 0; iv < this->GetNbDof(); iv++)
                          for (int jv = 0; jv < this->GetNbDof(); jv++)
                            {
                              int i = IndexDof(iv);
                              int j = IndexDof(jv);
                              if ((i>= 0) && (j>= 0))
                                Mh(i, j) += poids*DotProd(phi_glob1(iv), phi_glob1(jv));
                            }
                        
                        for (int iv = 0; iv < this->GetNbDof(); iv++)
                          for (int jv = phi_glob2.GetM()-nb_inside; jv < phi_glob2.GetM(); jv++)
                            {
                              int i = IndexDof(iv);
                              int j = jv - phi_glob2.GetM()+nb_inside;
                              if ((i>= 0) && (j>= 0))
                                Bh(i, j) += poids*DotProd(phi_glob1(iv), phi_glob2(jv));
                            }
                      }
                    
                    GetInverse(Mh);
                    Matrix<Real_wp> Ah(nb_inside, nb_inside);            
                    Ah.Fill(0);
                    Mlt(Mh, Bh, Ah);
                    
                    //DISP(num_loc); DISP(rot); DISP(Ah);
                    if (rot == 0)
                      {
                        bool test_identity = true;
                        for (int i = 0; i < nb_inside; i++)
                          for (int j = 0; j < nb_inside; j++)
                            {
                              if (j == i)
                                {
                                  if (abs(Ah(i, j) - 1.0) > threshold)
                                    test_identity = false;
                                }
                              else
                                {
                                  if (abs(Ah(i, j)) > threshold)
                                    test_identity = false;
                                }
                            }
                        
                        if (display_message)
                          if (!test_identity)
                            {
                              cout << "Warning : Operator not equal to identity for face "
                                   << num_loc << endl;
                              DISP(Ah);
                            }
                        
                      }
                    
                    Real_wp threshold = 1000*epsilon_machine;
                    // pour tous les ddls de la face
                    Vector<bool> DofUsed(nb_inside); DofUsed.Fill(false);            
                    for (int iv = 0; iv < this->GetNbDof(); iv++)
                      {
                        int i = IndexDof(iv);
                        if (i>= 0)
                          {
                            int nb_coef = 0;
                            for (int j = 0; j < nb_inside; j++)
                              if (abs(Ah(j, i)) > threshold)
                                nb_coef++;
                            
                            if (nb_coef == 0)
                              {
                                cout << "Cas impossible " << endl;
                                abort();
                              }
                            else if (nb_coef == 1)
                              {
                                // on est tombe sur un autre ddl, avec eventuellement un signe
                                for (int j = 0; j < nb_inside; j++)
                                  if (abs(Ah(j, i)) > threshold)
                                    {
                                      FaceRot(rot, i) = j;
                                      if ( abs(Ah(j, i) -1.0) <= threshold)
                                        OppositeSigneDof(rot, i) = false;
                                      else if ( abs(Ah(j, i) + 1.0) <= threshold)
                                        OppositeSigneDof(rot, i) = true;
                                      else
                                        {
                                          // this dof can't be handled with a sign
                                          linear_combination = true;
                                        }
                                      
                                      DofUsed(j) = true;
                                    }                       
                              }
                            else
                              {
                                // no direct relation between dofs
                                linear_combination = true;
                              }
                          }
                      }
                    
                    CoefCombination(rot) = Ah;

		    //DISP(rot); DISP(CoefCombination(rot));
                  }
                

		//DISP(linear_combination);
		//DISP(FaceRot); DISP(OppositeSigneDof);
                
                if (linear_combination)
                  {
                    if (quad)
                      nmap.SetFacesDofRotationQuad(this->GetOrder(), CoefCombination);
                    else
                      nmap.SetFacesDofRotationTri(this->GetOrder(), CoefCombination);
                  }
                else
                  {
                    if (quad)
                      {
                        nmap.SetFacesDofRotationQuad(this->GetOrder(), FaceRot);
                        nmap.SetSignDofRotationQuad(this->GetOrder(), OppositeSigneDof);
                      }
                    else
                      {
                        nmap.SetFacesDofRotationTri(this->GetOrder(), FaceRot);
                        nmap.SetSignDofRotationTri(this->GetOrder(), OppositeSigneDof);
                      }
                  }
                
              }
          }
      }
  }
  

  //! computes the sign of dofs after symmetry of edges
  /*!
    This function is only a function to check that restriction of functions
    on edges are the same between the edges, and that a symmetry
    of an edge may change the sign only.
    This function does not modify the finite element class, neither nmap
   */
  void ElementReference<Dimension3, 2>::
  FindHcurlSignEdge(NumberMap& nmap) const
  {
    VectR3 s;
    switch (this->GetNbVertices())
      {
      case 4 :
        {
          s.Reallocate(4);
          s(0).Init(0, 0, 0);
          s(1).Init(1, 0, 0);
          s(2).Init(0, 1, 0);
          s(3).Init(0, 0, 1);
        }
        break;
      case 5 :
        {
          s.Reallocate(5);
          s(0).Init(-1, -1, 0);
          s(1).Init(1, -1, 0);
          s(2).Init(1, 1, 0);
          s(3).Init(-1, 1, 0);
          s(4).Init(0, 0, 1);
        }
        break;
      case 6 :
        {
          s.Reallocate(6);
          s(0).Init(0, 0, 0);
          s(1).Init(1, 0, 0);
          s(2).Init(0, 1, 0);
          s(3).Init(0, 0, 1);
          s(4).Init(1, 0, 1);
          s(5).Init(0, 1, 1);
        }
        break;
      case 8 :
        {
          s.Reallocate(8);
          s(0).Init(0, 0, 0);
          s(1).Init(1, 0, 0);
          s(2).Init(1, 1, 0);
          s(3).Init(0, 1, 0);
          s(4).Init(0, 0, 1);
          s(5).Init(1, 0, 1);
          s(6).Init(1, 1, 1);
          s(7).Init(0, 1, 1);
        }
        break;
      }

    int nb_quad = this->Points1D().GetM();
    int nb_dof_edge = nmap.GetNbDofEdge(this->order);
    if (nb_dof_edge == 0)
      {
        cout << "You must have at least one dof per edge " << endl;
        abort();
      }

    Real_wp threshold = 1e4*epsilon_machine;    
    Array3D<Real_wp> ValPhiEdge(this->GetNbEdges(), nb_dof_edge, nb_quad);
    ValPhiEdge.Fill(0);
    
    // computing tangent for each edge, and values of basis functions
    VectR3 tangente(this->GetNbEdges()), phi;
    int type_elt = this->GetHybridType();
    for (int num_loc = 0; num_loc < this->GetNbEdges(); num_loc++)
      {
        int n1 = MeshNumbering<Dimension3>::FirstExtremityEdge(type_elt, num_loc);
        int n2 = MeshNumbering<Dimension3>::SecondExtremityEdge(type_elt, num_loc);
        tangente(num_loc) = s(n2) - s(n1);
        //DISP(num_loc); DISP(tangente(num_loc));
        for (int k = 0; k < nb_quad; k++)
          {
            R3 ptA = (1.0-this->Points1D(k))*s(n1) + this->Points1D(k)*s(n2);
            this->ComputeValuesPhiRef(ptA, phi);
            for (int j = 0; j < nb_dof_edge; j++)
              {
                int num_dof = num_loc*nb_dof_edge + j;
                ValPhiEdge(num_loc, j, k) = DotProd(phi(num_dof), tangente(num_loc));
                for (int num_loc2 = 0; num_loc2 < this->GetNbEdges(); num_loc2++)
                  if (num_loc2 != num_loc)
                    {
                      int num_dof2 = num_loc2*nb_dof_edge + j;
                      Real_wp val = DotProd(phi(num_dof2), tangente(num_loc));
                      
                      if (abs(val) > threshold)
                        {
                          cout << "Basis function "<<num_dof2 << " is not null on edge "
                               << num_loc2 << endl;
                          
                          abort();
                        }
                    }
              }
          }
      }
    
    //DISP(ValPhiEdge);
    // we check tangential traces on each edge
    for (int num_loc = 0; num_loc < this->GetNbEdges(); num_loc++)
      {
        Matrix<Real_wp> A(nb_dof_edge, nb_dof_edge), B = A, M = A;
        A.Fill(0); B.Fill(0); M.Fill(0);
        for (int j = 0; j < nb_dof_edge; j++)
          {
            // we check that it is the same as first edge
            for (int k = 0; k < nb_quad; k++)
              if (abs(ValPhiEdge(num_loc, j, k) - ValPhiEdge(0, j, k)) > threshold)
                {
                  cout << "values of tangential trace of function " << j 
                       << " on edge " <<  num_loc 
                       << " does not match with function of first edge " << endl;
                  
                  DISP(ValPhiEdge(num_loc, j, k));
                  DISP(ValPhiEdge(0, j, k));
                  abort();
                }
            
            for (int k = 0; k < nb_quad; k++)
              {
                Real_wp val2 = ValPhiEdge(num_loc, j, nb_quad-1-k);
                Real_wp val = ValPhiEdge(num_loc, j, k);
                for (int j2 = 0; j2 < nb_dof_edge; j2++)
                  {
                    Real_wp val1 = ValPhiEdge(num_loc, j2, k);
                    B(j, j2) += val2*this->Weights1D(k)*val1;
                    M(j, j2) += val*this->Weights1D(k)*val1;
                  }
              }
          }
        
        GetInverse(M);
        Mlt(B, M, A);
        //DISP(num_loc); DISP(A);
        
        // we check that matrix A is a permutation matrix with sign
        for (int j = 0; j < nb_dof_edge; j++)
          {
            int nb = 0;
            for (int jb = 0; jb < nb_dof_edge; jb++)
              if (abs(A(j, jb)) > threshold)
                nb++;
            
            if (nb != 1)
              {
                abort();
              }
            
            //int signe = 0;
            for (int jb = 0; jb < nb_dof_edge; jb++)
              {
                Real_wp val = A(j, jb);
                if (abs(val - 1.0) <= threshold)
                  {
                    //signe = 1;
                  }
                else if (abs(val + 1.0) <= threshold)
                  {
                    //signe = -1;
                  }
                else if (abs(A(j, jb)) > threshold)
                  {
                    abort();
                  }
              }
            // DISP(num_loc); DISP(j); DISP(j2); DISP(signe);
          }
      }
  }

   
  //! computation of values of u on nodal points of a face
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points of face num_loc
    \param[in] num_loc face number
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 2>
  ::ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    typename Vector2::value_type vx, vy, vz;
    for (int j = 0; j < Unode.GetM()/3; j++)
      {
        int node = FacesNode(j, num_loc);
        vx = 0; vy = 0; vz = 0;
        for (int k = 0; k < this->nb_dof_loc; k++)
          {
            vx += Un(k)*this->ValuePhi_Nodal(k, 3*node);
            vy += Un(k)*this->ValuePhi_Nodal(k, 3*node+1);
            vz += Un(k)*this->ValuePhi_Nodal(k, 3*node+2);
          }
        
        Unode(3*j) = vx; Unode(3*j+1) = vy; Unode(3*j+2) = vz;
      }
  }


  //! computation of curl u on nodal points of a face
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode curl of u on nodal points of face num_loc
    \param[in] num_loc face number
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 2>
  ::ComputeCurlBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    typename Vector2::value_type vx, vy, vz;
    for (int j = 0; j < Unode.GetM()/3; j++)
      {
        int node = FacesNode(j, num_loc);
        vx = 0; vy = 0; vz = 0;
        for (int k = 0; k < this->nb_dof_loc; k++)
          {
            vx += Un(k)*this->CurlPhi_Nodal(k, 3*node);
            vy += Un(k)*this->CurlPhi_Nodal(k, 3*node+1);
            vz += Un(k)*this->CurlPhi_Nodal(k, 3*node+2);
          }
        
        Unode(3*j) = vx; Unode(3*j+1) = vy; Unode(3*j+2) = vz;
      }
  }


  //! Integration against curl of basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f curl(\varphi_i) dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = curl(\varphi_i)(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 2>::ApplyRhGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < nb_dof_loc; i++)
      {
        val = Real_wp(0);
        for (int j = 0; j < nb_points_quadrature_inside; j++)
          val += Curl_Phi(i, j)(0)*Vh(3*j)
            + Curl_Phi(i, j)(1)*Vh(3*j+1) + Curl_Phi(i, j)(2)*Vh(3*j+2);
	
        Uh(i) = val;
      }
  }

  
  //! computation of curl of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh curl of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j curl phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = curl phi_i(xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 2>::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    TinyVector<typename Vector2::value_type, 3> curl;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        curl.Fill(0);
        for (int j = 0; j < nb_dof_loc; j++)
          Add(Uh(j), Curl_Phi(j, i), curl);
        
        Vh(3*i) = curl(0); Vh(3*i+1) = curl(1); Vh(3*i+2) = curl(2);
      }
  }
  

  //! Integration against basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Ch Vh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 2>::ApplyChGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < nb_dof_loc; i++)
      {
	val = 0;
	for (int j = 0; j < nb_points_quadrature_inside; j++)
	  val += Value_PhiVec(i, j)(0)*Vh(3*j)
            + Value_PhiVec(i, j)(1)*Vh(3*j+1)
            + Value_PhiVec(i, j)(2)*Vh(3*j+2);
	
	Uh(i) = val;
      }
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
  void ElementReference<Dimension3, 2>
  ::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    typename Vector2::value_type vx, vy, vz;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	vx = 0; vy = 0; vz = 0;
	for (int j = 0; j < nb_dof_loc; j++)
	  {
            vx += Value_PhiVec(j, i)(0)*Uh(j);
            vy += Value_PhiVec(j, i)(1)*Uh(j);
            vz += Value_PhiVec(j, i)(2)*Uh(j);
          }
	
	Vh(3*i) = vx; Vh(3*i+1) = vy; Vh(3*i+2) = vz;
      }
  }


  //! we add constant mass matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass coefficient
    \param[inout] A matrix modified
    A(m:, n:) = A(m:, n:) + mass M
    where M is the mass matrix, M_ij = \int \varphi_j \varphi_i dx
   */
  template<class T, class Prop>
  void ElementReference<Dimension3, 2>::
  AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& mass, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;

    for (int i = 0; i < nb_dof_loc; i++)
      {
	for (int j = 0; j < nb_dof_loc; j++)
	  val(j) = const_mass_matrix(0, 0)(i, j)*mass(0, 0)
	    + const_mass_matrix(0, 1)(i, j)*mass(0, 1)
	    + const_mass_matrix(1, 0)(i, j)*mass(1, 0)
	    + const_mass_matrix(2, 0)(i, j)*mass(2, 0)
	    + const_mass_matrix(0, 2)(i, j)*mass(0, 2)
	    + const_mass_matrix(1, 1)(i, j)*mass(1, 1) 
	    + const_mass_matrix(1, 2)(i, j)*mass(1, 2)
	    + const_mass_matrix(2, 1)(i, j)*mass(2, 1)
	    + const_mass_matrix(2, 2)(i, j)*mass(2, 2);

	A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
      }                    
  }
  

  //! we add constant stiffness matrix C \int_K curl phi_j curl phi_i
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] C coefficients
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + \sum_p  C(p, q) S^{p,q}
    where S^{p,q} is a stiffness matrix equal to :
    (S^{p,q })_{i, j} = \int_K curl(\phi_j)_q  curl(\phi_i)_p  dx
   */
  template<class T, class Prop>
  void ElementReference<Dimension3, 2>::
  AddConstantStiffnessMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;

    for (int i = 0; i < nb_dof_loc; i++)
      {
	for (int j = 0; j < nb_dof_loc; j++)
	  val(j) = const_curl_matrix(0, 0)(i, j)*C(0, 0)
	    + const_curl_matrix(0, 1)(i, j)*C(0, 1)
	    + const_curl_matrix(1, 0)(i, j)*C(1, 0)
	    + const_curl_matrix(2, 0)(i, j)*C(2, 0)
	    + const_curl_matrix(0, 2)(i, j)*C(0, 2)
	    + const_curl_matrix(1, 1)(i, j)*C(1, 1)
	    + const_curl_matrix(1, 2)(i, j)*C(1, 2)
	    + const_curl_matrix(2, 1)(i, j)*C(2, 1)
	    + const_curl_matrix(2, 2)(i, j)*C(2, 2);

	A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
      }    
  }


  //! adds \int C \nabla \times varphi_j \nabla \times varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] C tensor C on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int C \nabla \times varphi_j \nabla \times varphi_i
   */
  template<class T, class Prop>
  void ElementReference<Dimension3, 2>
  ::AddVariableStiffnessMatrixGen(int off_row, int off_col,
				  const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
				  VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    T vloc; TinyVector<T, 3> vec_u;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
	    SetComplexZero(vloc);
	    for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	      {
		Mlt(C(k), this->Curl_Phi(i, k), vec_u);
		vloc += DotProd(vec_u, this->Curl_Phi(j, k));
	      }
	    
	    val(j) = vloc;
	  }
	
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }


  //! adds \int A varphi_j \varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A tensor A on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int A varphi_j  varphi_i
   */
  template<class T, class Prop>
  void ElementReference<Dimension3, 2>
  ::AddVariableMassMatrixGen(int off_row, int off_col,
			     const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
			     VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    T vloc; TinyVector<T, 3> vec_u;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
	    SetComplexZero(vloc);
	    for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	      {
		Mlt(C(k), this->Value_PhiVec(i, k), vec_u);
		vloc += DotProd(vec_u, this->Value_PhiVec(j, k));
	      }
	    
	    val(j) = vloc;
	  }
	
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }


  //! computes integrals against curl of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension3, 2>
  ::ComputeIntegralCurlRef(const Vector1 & feval, Vector2& res) const
  {
    this->ApplyRh(feval, res);
  }
  
  
  //! computes surface integrals against curl of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension3, 2>
  ::ComputeIntegralSurfaceCurlRef(const Vector1 & feval, Vector2& res, int n) const
  {
    typename Vector2::value_type one; SetComplexOne(one);
    res.Fill(0);
    this->ApplyNablaSh(one, n, feval, res);
  }
  
  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void ElementReference<Dimension3, 2>
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    contrib.Reallocate(feval.GetM()/3);
    for (int i = 0; i < contrib.GetM(); i++)
      contrib(i) = feval(3*i)*tangente_dof(i)(0) 
        + feval(3*i+1)*tangente_dof(i)(1) + feval(3*i+2)*tangente_dof(i)(2);
    
    SolveQR(MatProjectionDof, tauProjectionDof, contrib);
  }

  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void ElementReference<Dimension3, 2>
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    VectReal_wp xreal(feval.GetM()/3), ximag(feval.GetM()/3);
    for (int i = 0; i < xreal.GetM(); i++)
      {
        Complex_wp vloc = feval(3*i)*tangente_dof(i)(0)
          + feval(3*i+1)*tangente_dof(i)(1) + feval(3*i+2)*tangente_dof(i)(2);
        
        xreal(i) = real(vloc);
        ximag(i) = imag(vloc);
      }
    
    SolveQR(MatProjectionDof, tauProjectionDof, xreal);
    SolveQR(MatProjectionDof, tauProjectionDof, ximag);
    
    //DISP(MatProjectionDof.GetM());
    //DISP(MatProjectionDof.GetN());
    
    contrib.Reallocate(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      contrib(i) = Complex_wp(xreal(i), ximag(i));
  }


  /************************
   * VolumeHcurlReference *
   ************************/

  
  void VolumeHcurlReference::ConstructFiniteElement()
  {
    elt_Hcurl.order = elt_H1.GetOrder();
    elt_Hcurl.order_quad = elt_H1.GetQuadratureOrder();
    elt_Hcurl.nb_points_quadrature_inside = elt_H1.GetNbPointsQuadratureInside();
    elt_Hcurl.nb_points_quadrature_boundaries = elt_H1.GetNbPointsQuadBoundaries();

    elt_Hcurl.nb_dof_loc = 3*elt_H1.GetNbDof();
    elt_Hcurl.nb_dof_boundaries = 0;
    elt_Hcurl.nb_points_dof_inside = 0;
    
    elt_Hcurl.mass_lumping = elt_H1.LumpedMassMatrix();
    elt_Hcurl.mass_lumping_ortho = elt_H1.LumpedMassMatrix();
    elt_Hcurl.diagonal_mass = elt_H1.DiagonalMassMatrix();

    elt_Hcurl.discontinuous_element = true;
    
    elt_Hcurl.optimized_mass_matrix = elt_H1.OptimizedComputationMassMatrix();
    elt_Hcurl.optimized_elem_matrix = elt_H1.OptimizedComputationElementaryMatrix();
    elt_Hcurl.sparse_mass_matrix = elt_H1.SparseMassMatrix();
    elt_Hcurl.linear_sparse_mass_matrix = elt_H1.LinearSparseMassMatrix();

    elt_Hcurl.use_quadrature_for_rh = elt_H1.UseQuadraturePointsForRh();
    elt_Hcurl.use_quadrature_for_sh = elt_H1.UseQuadraturePointsForSh();
    elt_Hcurl.use_quadrature_free_sh = elt_H1.UseQuadratureFreeSh();

    elt_Hcurl.type_integration_tri = elt_H1.GetTypeIntegrationTriangle();
    elt_Hcurl.type_integration_quad = elt_H1.GetTypeIntegrationQuadrangle();
    elt_Hcurl.num_quad_points_surf = elt_H1.GetQuadNumbersBoundary();
  }

  
  void VolumeHcurlReference::ComputeValuesPhiRef(const R3& x, VectR3& phi) const
  {
    VectReal_wp phi_scal;
    elt_H1.ComputeValuesPhiRef(x, phi_scal);

    int Ns = elt_H1.GetNbDof();
    phi.Reallocate(3*Ns);
    for (int i = 0; i < Ns; i++)
      {
	phi(i).Init(phi_scal(i), Real_wp(0), Real_wp(0));
	phi(Ns+i).Init(Real_wp(0), phi_scal(i), Real_wp(0));
	phi(2*Ns+i).Init(Real_wp(0), Real_wp(0), phi_scal(i));
      }
  }


  void VolumeHcurlReference::GetValuePhiOnQuadraturePoint(int k, VectR3& phi) const
  {
    VectReal_wp phi_scal;
    elt_H1.GetValuePhiOnQuadraturePoint(k, phi_scal);

    int N = elt_H1.GetNbDof();
    phi.Reallocate(3*N);
    for (int i = 0; i < N; i++)
      {
	phi(i).Init(phi_scal(i), Real_wp(0), Real_wp(0));
	phi(N+i).Init(Real_wp(0), phi_scal(i), Real_wp(0));
	phi(2*N+i).Init(Real_wp(0), Real_wp(0), phi_scal(i));
      }
  }


  template<class T>
  void VolumeHcurlReference::ApplyChGen(const Vector<T>& U, Vector<T>& V) const
  {
    Vector<T> feval(elt_H1.GetNbPointsQuadratureInside());
    Vector<T> contrib(elt_H1.GetNbDof());
    int Ns = elt_H1.GetNbDof();

    for (int i = 0; i < feval.GetM(); i++)
      feval(i) = U(3*i);

    elt_H1.ApplyCh(feval, contrib);

    for (int i = 0; i < contrib.GetM(); i++)
      V(i) = contrib(i);

    for (int i = 0; i < feval.GetM(); i++)
      feval(i) = U(3*i+1);

    elt_H1.ApplyCh(feval, contrib);

    for (int i = 0; i < contrib.GetM(); i++)
      V(Ns+i) = contrib(i);

    for (int i = 0; i < feval.GetM(); i++)
      feval(i) = U(3*i+2);

    elt_H1.ApplyCh(feval, contrib);

    for (int i = 0; i < contrib.GetM(); i++)
      V(2*Ns+i) = contrib(i);    
  }
#endif


  /*******************
   * H(div) elements *
   *******************/


#ifdef MONTJOIE_WITH_NODAL_HDIV
  //! Destructor
  ElementReference<Dimension3, 3>::~ElementReference()
  {
    if (element_tri_surf != NULL)
      {
	delete element_tri_surf;
	element_tri_surf = NULL;
      }

    if (element_quad_surf != NULL)
      {
	delete element_quad_surf;
	element_quad_surf = NULL;
      }
  }


  //! returns 2-D triangular finite element (normal traces on triangular faces)
  const ElementReference<Dimension2, 1>& 
  ElementReference<Dimension3, 3>::GetTriangularSurfaceFiniteElement() const
  {
    if (element_tri_surf == NULL)
      {
	cout << "Null pointer for GetTriangularSurfaceFiniteElement" << endl;
	abort();
      }

    return *element_tri_surf;
  }
  
  
  //! returns 2-D quadrangular finite element (normal traces on quadrangular faces)
  const ElementReference<Dimension2, 1>& 
  ElementReference<Dimension3, 3>::GetQuadrangularSurfaceFiniteElement() const
  {
    if (element_quad_surf == NULL)
      {
	cout << "Null pointer for GetQuadrangularSurfaceFiniteElement" << endl;
	abort();
      }
    
    return *element_quad_surf;
  }


  const ElementReference<Dimension2, 1>& 
  ElementReference<Dimension3, 3>::GetSurfaceFiniteElement(int n) const
  {
    if (this->elt_geom.IsLocalFaceQuadrilateral(n))
      return GetQuadrangularSurfaceFiniteElement();
    else
      return GetTriangularSurfaceFiniteElement();
  }
  
  
  size_t ElementReference<Dimension3, 3>::GetMemorySize() const
  {
    size_t taille = VolumeReference<3>::GetMemorySize() + sizeof(FiniteElementHdiv<Dimension3>);
    if (element_tri_surf != NULL)
      taille += element_tri_surf->GetMemorySize();

    if (element_quad_surf != NULL)
      taille += element_quad_surf->GetMemorySize();
    
    return taille;
  }
  

  void ElementReference<Dimension3, 3>
  ::ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int num_loc) const
  {
    VectReal_wp feval_normal(feval.GetM()/3);
    R3 vec_u;
    for (int i = 0; i < feval_normal.GetM(); i++)
      {
	vec_u.Init(feval(3*i), feval(3*i+1), feval(3*i+2));
	feval_normal(i) = DotProd(this->NormaleLoc(num_loc), vec_u);
      }
    
    if (this->elt_geom.IsLocalFaceQuadrilateral(num_loc))
      this->element_quad_surf->ComputeProjectionDofRef(feval_normal, contrib);
    else
      this->element_tri_surf->ComputeProjectionDofRef(feval_normal, contrib);

    this->ModifySignProjectionSurface(contrib, num_loc);
  }


  void ElementReference<Dimension3, 3>
  ::ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int num_loc) const
  {
    VectComplex_wp feval_normal(feval.GetM()/3);
    R3_Complex_wp vec_u;
    for (int i = 0; i < feval_normal.GetM(); i++)
      {
	vec_u.Init(feval(3*i), feval(3*i+1), feval(3*i+2));
	feval_normal(i) = DotProd(this->NormaleLoc(num_loc), vec_u);
      }

    if (this->elt_geom.IsLocalFaceQuadrilateral(num_loc))
      this->element_quad_surf->ComputeProjectionDofRef(feval_normal, contrib);
    else
      this->element_tri_surf->ComputeProjectionDofRef(feval_normal, contrib);

    this->ModifySignProjectionSurface(contrib, num_loc);
  }  


  void ElementReference<Dimension3, 3>
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
  }
  

  void ElementReference<Dimension3, 3>
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
  }


  //! computation of mass matrix, stiffness matrix, etc (for facet elements)  
  /*!
    \param[in] Fb leaf finite element class
   */
  void ElementReference<Dimension3, 3>::ConstructHdivElementaryMatrix()
  {
    int nb_points_quadrature = this->GetNbPointsQuadrature();
    Vector<R3> phi;
    Value_PhiVec.Reallocate(nb_dof_loc, nb_points_quadrature);        
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	this->ComputeValuesPhiRef(this->PointsND(i), phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  Value_PhiVec(j, i) = phi(j);
      }
    
    // mass matrix on reference element
    this->mass_matrix.Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    Real_wp vloc;
    
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = i; j < this->nb_dof_loc; j++)
	{
	  vloc = 0.0;
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*DotProd(this->Value_PhiVec(i, k), this->Value_PhiVec(j, k));
          
	  this->mass_matrix(i,j) = vloc;
	}
    
    //VectReal_wp lambda; 
    //GetEigenvalues(this->mass_matrix, lambda);
    //DISP(lambda);
	  
    this->mass_matrix_chol = this->mass_matrix;
    GetCholesky(this->mass_matrix_chol);
    
    // projection on dofs
    this->MatProjectionDof.Reallocate(Dimension::dim_N*this->GetNbPointsDofInside(), this->nb_dof_loc);
    for (int i = 0; i < this->GetNbPointsDofInside(); i++)
      {
        this->ComputeValuesPhiRef(this->PointsDofND(i), phi);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
	    this->MatProjectionDof(3*i, j) = phi(j)(0);
	    this->MatProjectionDof(3*i+1, j) = phi(j)(1);
	    this->MatProjectionDof(3*i+2, j) = phi(j)(2);
	  }
      } 
    
    GetQR(this->MatProjectionDof, this->tauProjectionDof);
    
    ValuePhi_Nodal.Reallocate(nb_dof_loc, Dimension::dim_N*this->GetNbPointsNodalElt());
    for (int i = 0; i < this->GetNbPointsNodalElt(); i++)
      {
	this->ComputeValuesPhiRef(this->PointsNodalND(i), phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  {
	    ValuePhi_Nodal(j, 3*i) = phi(j)(0);
	    ValuePhi_Nodal(j, 3*i+1) = phi(j)(1);
	    ValuePhi_Nodal(j, 3*i+2) = phi(j)(2);
	  }
      }
    
    Vector<Real_wp> div_phi;
    DivPhi_Nodal.Reallocate(nb_dof_loc, this->GetNbPointsNodalElt());
    for (int i = 0; i < this->GetNbPointsNodalElt(); i++)
      {
	this->ComputeDivPhiRef(this->PointsNodalND(i), div_phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  DivPhi_Nodal(j, i) = div_phi(j);
      }
    
    Div_Phi.Reallocate(nb_dof_loc, nb_points_quadrature);
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	this->ComputeDivPhiRef(this->PointsND(i), div_phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  Div_Phi(j, i) = div_phi(j);
      }
    
    // mass and stiffness matrix for constant coefficients
    // \int_K C phi_i phi_j   and   \int_K D div(phi_i)  div(phi_j)
    // where C is a 3x3 symmetric matrix
    for (int p = 0; p < Dimension::dim_N; p++)
      for (int q = 0; q < Dimension::dim_N; q++)
	this->const_mass_matrix(p, q).Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    
    this->const_div_matrix.Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = 0; j < this->nb_dof_loc; j++)
        {
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i, k)(0)*this->Value_PhiVec(j, k)(0);
          
          this->const_mass_matrix(0, 0)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(0)*this->Value_PhiVec(j,k)(1);
          
          this->const_mass_matrix(0, 1)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(1)*this->Value_PhiVec(j,k)(0);
          
          this->const_mass_matrix(1, 0)(i, j) = vloc;

	  vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(1)*this->Value_PhiVec(j,k)(1);
          
          this->const_mass_matrix(1, 1)(i, j) = vloc;
          
	  vloc = 0;
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	    vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(0)*this->Value_PhiVec(j,k)(2);
	  
	  this->const_mass_matrix(0, 2)(i, j) = vloc;
	  
	  vloc = 0;
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	    vloc +=  this->WeightsND(k)*this->Value_PhiVec(i,k)(2)*this->Value_PhiVec(j,k)(0);
	  
	  this->const_mass_matrix(2, 0)(i, j) = vloc;
	  
	  vloc = 0;
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	    vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(1)*this->Value_PhiVec(j,k)(2);
	  
	  this->const_mass_matrix(1, 2)(i, j) = vloc;
	  
	  vloc = 0;
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	    vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(2)*this->Value_PhiVec(j,k)(1);
	  
	  this->const_mass_matrix(2, 1)(i, j) = vloc;
	  
	  vloc = 0;
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	    vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(2)*this->Value_PhiVec(j,k)(2);
	  
	  this->const_mass_matrix(2, 2)(i, j) = vloc;
	  
          // div div matrix
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Div_Phi(i, k)*this->Div_Phi(j, k);
          
          this->const_div_matrix(i, j) = vloc;
        }
  }


  //! filling FacesDof and checking that dofs are numbered correctly 
  //! (first dof on faces, then interior)
  void ElementReference<Dimension3, 3>::FindDofsOnFace()
  {
    int N = this->GetNbDof();
    if (N <= 0)
      {
        cout << "no degree of freedom "<<endl;
        abort();
      }
    
    Real_wp threshold = 1e4*epsilon_machine;
    // for each face we are computing phi \times n 
    Matrix<bool> normal_proj_null(N, this->GetNbBoundaries());
    normal_proj_null.Fill(true);
    VectR3 phi; Real_wp phi_n;
    //DISP(this->GetNbBoundaries());
    for (int num_loc = 0; num_loc < this->GetNbBoundaries(); num_loc++)
      for (int i = 0; i < this->GetNbQuadBoundary(num_loc); i++)
        {
          int k = this->GetQuadNumber(num_loc, i);
          this->GetValuePhiOnQuadraturePoint(k, phi);
	  for (int j = 0; j < N; j++)
            {
              phi_n = DotProd(phi(j), this->NormaleLoc(num_loc));
              //DISP(i); DISP(j); DISP(phi_n);
              if (abs(phi_n) > threshold)
                {
                  normal_proj_null(j, num_loc) = false;
                }
            }
        }
    
    // finding dofs inside element
    Vector<bool> DofInside(N); DofInside.Fill(true);
    int nb = 0;
    for (int i = 0; i < this->GetNbDof(); i++)
      {
        for (int j = 0; j < this->GetNbBoundaries(); j++)
          if (!normal_proj_null(i, j))
            {
              DofInside(i) = false;
              // cout << "tangential component of dof " << i << " on face " 
              // << j << " is different from 0 " << endl;     
            }
        
        if (DofInside(i))
          nb++;
      }
    
    for (int i = N-nb; i < N; i++)
      if (!DofInside(i))
        {
          cout << "Internal dofs must be numbered at the end " << endl;
          abort();
        }
    
    if (this->nb_dof_boundaries <= 0)
      this->nb_dof_boundaries = N - nb;
    else
      {
        if (this->nb_dof_boundaries != N - nb)
          {
            cout << "Number of dofs on the boundary is noticed as " 
                 << this->nb_dof_boundaries << endl;
            cout << "But we found " << N-nb << " dofs on the boundary " << endl;
            abort();
          }
      }
    
    // finding dofs inside faces
    int nb_max = this->nb_dof_boundaries;
    IVect nb_dof_face(this->GetNbBoundaries()); nb_dof_face.Fill(0);
    for (int num_loc = this->GetNbBoundaries()-1; num_loc >= 0; num_loc--)
      {
        DofInside.Fill(true);
        nb = 0;
        for (int i = 0; i < nb_max; i++)
          {
            for (int num_loc2 = 0; num_loc2 < this->GetNbBoundaries(); num_loc2++)
              if (num_loc2 != num_loc)
                if (!normal_proj_null(i, num_loc2))
                  DofInside(i) = false;
            
            if (DofInside(i))
              nb++;
          }
        
        for (int i = nb_max-nb; i < nb_max; i++)
          if (!DofInside(i))
            {
              cout << "dof number for this face is incorrect " << endl;
              abort();
            }       
        
        nb_dof_face(num_loc) = nb;
        nb_max -= nb;
      }
    
    // filling FacesDof
    nb_max = 0;
    for (int j = 0; j < this->GetNbBoundaries(); j++)
      nb_max = max(nb_dof_face(j), nb_max);
    
    // eventually filling FacesDof
    FacesDof.Reallocate(nb_max, this->GetNbBoundaries());
    FacesDof.Fill(-1);
    int offset = 0;
    for (int j = 0; j < this->GetNbBoundaries(); j++)
      {
        nb = 0;
        for (int i = 0; i < nb_dof_face(j); i++)
          FacesDof(nb++, j) = offset+i;
        
        offset += nb_dof_face(j);
      }

    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);
  }


  //! computing operator to apply when the faces are rotated (for edge elements)
  void ElementReference<Dimension3, 3>::
  FindHdivLinearCombinationRotation(NumberMap& nmap,
				    const VectR2& PointsTri, const VectR2& PointsQuad,
				    const Matrix<Real_wp>& ValuePhiTri,
				    const Matrix<Real_wp>& ValuePhiQuad,
				    bool display_message) const
  {
    // vertex numbers of triangular face after rotation
    TinyVector<TinyVector<int, 3>, 6> rot_tri_num;
    rot_tri_num(0).Init(0, 1, 2);
    rot_tri_num(1).Init(1, 2, 0);
    rot_tri_num(2).Init(2, 0, 1);
    rot_tri_num(3).Init(0, 2, 1);
    rot_tri_num(4).Init(1, 0, 2);
    rot_tri_num(5).Init(2, 1, 0);
    
    // vertex numbers of quadrangular face after rotation
    TinyVector<TinyVector<int, 4>, 8> rot_quad_num;
    rot_quad_num(0).Init(0, 1, 2, 3);
    rot_quad_num(1).Init(1, 2, 3, 0);
    rot_quad_num(2).Init(2, 3, 0, 1);
    rot_quad_num(3).Init(3, 0, 1, 2);
    rot_quad_num(4).Init(0, 3, 2, 1);
    rot_quad_num(5).Init(1, 0, 3, 2);
    rot_quad_num(6).Init(2, 1, 0, 3);
    rot_quad_num(7).Init(3, 2, 1, 0);
    
    Matrix<int> IndexDofFace(this->GetNbDof(), this->GetNbBoundaries());
    IndexDofFace.Fill(-1);
    bool first_tri = true, first_quad = true;
    int type_elt = this->GetHybridType();
    Real_wp threshold = 1e4*epsilon_machine;
    for (int num_loc = 0; num_loc < this->GetNbBoundaries(); num_loc++)
      {
        bool quad = MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(type_elt, num_loc);
        bool check_face = false;
        if (quad)
          {
            if (first_tri)
              {
                first_tri = false;
                check_face = true;
              }
          }
        else
          {
            if (first_quad)
              {
                first_quad = false;
                check_face = true;
              }
          }
        
	check_face = true;
        if (check_face)
          {
            Matrix<int> FaceRot;
            Matrix<bool> OppositeSigneDof;
            bool linear_combination = false;
            
            // tetrahedron
            Vector<bool> dof_on_face(this->GetNbDof()); dof_on_face.Fill(true);
            for (int i = this->GetNbDofBoundaries(); i < this->GetNbDof(); i++)
              dof_on_face(i) = false;
            
            int N = this->GetNbQuadBoundary(num_loc);
            VectR3 phi; VectReal_wp phi2; Real_wp phi_n;
            for (int num_loc2 = 0; num_loc2 < this->GetNbBoundaries(); num_loc2++)
              if (num_loc2 != num_loc)
                for (int i = 0; i < this->GetNbQuadBoundary(num_loc2); i++)
                  {
                    int k = this->GetQuadNumber(num_loc2, i);
                    this->GetValuePhiOnQuadraturePoint(k, phi);
                    
                    for (int j = 0; j < this->GetNbDof(); j++)
                      { 
                        phi_n = DotProd(phi(j), this->NormaleLoc(num_loc2));
                        if (abs(phi_n) > threshold)
                          {
                            dof_on_face(j) = false;
                          }
                      }
                  }           
            // DISP(dof_on_face);
            
            IVect IndexDof(this->GetNbDof()); IndexDof.Fill(-1);
            int nb_inside = 0;
            for (int j = 0; j < this->GetNbDof(); j++)
              if (dof_on_face(j))
                {
                  IndexDof(j) = nb_inside;
                  IndexDofFace(j, num_loc) = nb_inside;
                  nb_inside++;
                }
            // DISP(num_loc); DISP(nb_inside);
            
            int nb_orientations = 6; int nb_vert_on_face = 3;
            if (quad)
              {
                nb_orientations = 8;
                nb_vert_on_face = 4;
              }
            
            if (nb_inside > 0)
              {    
                FaceRot.Reallocate(nb_orientations, nb_inside); FaceRot.Fill(-1); 
                OppositeSigneDof.Reallocate(nb_orientations, nb_inside);
                OppositeSigneDof.Fill(false);
                Vector<Matrix<Real_wp> > CoefCombination(nb_orientations);
                
                // small mesh with one element
                Mesh<Dimension3> mesh;
                mesh.ReallocateVertices(this->GetNbVertices());
                
                switch (type_elt)
                  {
                  case 0 :
                    {
                      mesh.Vertex(0).Init(0, 0, 0);
                      mesh.Vertex(1).Init(1, 0, 0);
                      mesh.Vertex(2).Init(0, 1, 0);
                      mesh.Vertex(3).Init(0, 0, 1);
                    }
                    break;
                  case 1 :
                    {
                      mesh.Vertex(0).Init(-1, -1, 0);
                      mesh.Vertex(1).Init(1, -1, 0);
                      mesh.Vertex(2).Init(1, 1, 0);
                      mesh.Vertex(3).Init(-1, 1, 0);
                      mesh.Vertex(4).Init(0, 0, 1);
                    }
                    break;
                  case 2 :
                    {
                      mesh.Vertex(0).Init(0, 0, 0);
                      mesh.Vertex(1).Init(1, 0, 0);
                      mesh.Vertex(2).Init(0, 1, 0);
                      mesh.Vertex(3).Init(0, 0, 1);
                      mesh.Vertex(4).Init(1, 0, 1);
                      mesh.Vertex(5).Init(0, 1, 1);
                    }
                    break;
                  case 3 :
                    {
                      mesh.Vertex(0).Init(0, 0, 0);
                      mesh.Vertex(1).Init(1, 0, 0);
                      mesh.Vertex(2).Init(1, 1, 0);
                      mesh.Vertex(3).Init(0, 1, 0);
                      mesh.Vertex(4).Init(0, 0, 1);
                      mesh.Vertex(5).Init(1, 0, 1);
                      mesh.Vertex(6).Init(1, 1, 1);
                      mesh.Vertex(7).Init(0, 1, 1);
                    }
                  }
                
                mesh.ReallocateElements(1);
                int ref = 1; IVect num(this->GetNbVertices());
                num.Fill();
                mesh.Element(0).Init(num, ref);
                mesh.FindConnectivity();
                
                VectR3 s, sFace; Matrix3_3 DF1, DF2, invDF1, invDF2;
                R3 pt, ptGlob1, ptGlob2, normale_fj, phi_t;
		VectReal_wp phi1(this->GetNbDof());
                // we test all orientations
                for (int rot = 0; rot < nb_orientations; rot++)
                  {
                    IVect num(3);
                    if (quad)
                      num.Reallocate(4);
                    
                    sFace.Reallocate(num.GetM());
                    for (int j = 0; j < nb_vert_on_face; j++)
                      {
                        int jrot = -1;
                        if (quad)
                          jrot = rot_quad_num(rot)(j);
                        else
                          jrot = rot_tri_num(rot)(j);
                        
                        int nv = MeshNumbering<Dimension3>::
                          GetVertexNumberOfFace(type_elt, num_loc, j);
                        
                        num(jrot) = nv;
                        sFace(jrot) = mesh.Vertex(nv);
                      }
                    
                    // on calcule les DF de l'element 1 et 2
                    mesh.GetVerticesElement(0, s);
                    
                    // on boucle sur les points de quadrature de la face            
                    Matrix<Real_wp> Mh(nb_inside, nb_inside), Bh(nb_inside, nb_inside);
                    Mh.Fill(0); Bh.Fill(0); Real_wp x, y;
                    for (int ks = 0; ks < N; ks++)
                      {
                        Real_wp poids = this->WeightsQuadratureBoundary(ks, num_loc);
                        
                        // treating 3-D face of the element
                        int k = this->GetQuadNumber(num_loc, ks);
                        pt = this->PointsND(k);
                        this->FjLinear(s, pt, ptGlob1);
                        
                        this->DFjLinear(s, pt, DF1);
                        GetInverse(DF1, invDF1);
                        
                        Real_wp dsj;
                        this->GetNormale(invDF1, normale_fj, dsj, num_loc);
                        this->GetValuePhiOnQuadraturePoint(k, phi);
                        
                        // treating 3-D face from a reference triangle/quadrangle
                        int k2 = -1;
                        if (quad)
                          {
                            // searching point on 3-D face
                            for (int i2 = 0; i2 < ValuePhiQuad.GetN(); i2++)
                              {
                                x = PointsQuad(i2)(0);
                                y = PointsQuad(i2)(1);
                                ptGlob2 = (1.0-x)*(1.0-y)*sFace(0) + x*(1.0-y)*sFace(1)
                                  + x*y*sFace(2) + (1.0-x)*y*sFace(3);
                                
                                if (ptGlob2.Distance(ptGlob1) < threshold)
                                  k2 = i2;                            
                              }
                            
			    phi2.Reallocate(ValuePhiQuad.GetM());
                            for (int i2 = 0; i2 < ValuePhiQuad.GetM(); i2++)
			      phi2(i2) = ValuePhiQuad(i2, k2);
                          }
                        else
                          {
                            // searching point on 3-D face
                            for (int i2 = 0; i2 < ValuePhiTri.GetN(); i2++)
                              {
                                x = PointsTri(i2)(0);
                                y = PointsTri(i2)(1);
                                ptGlob2 = (1.0-x-y)*sFace(0) + x*sFace(1) + y*sFace(2);
                                if (ptGlob2.Distance(ptGlob1) < threshold)
                                  k2 = i2;                            
                              }
                            
			    phi2.Reallocate(ValuePhiTri.GetM());
                            for (int i2 = 0; i2 < ValuePhiTri.GetM(); i2++)
			      phi2(i2) = ValuePhiTri(i2, k2);
                          }
                        
                        // phi globaux
                        for (int j = 0; j < this->GetNbDof(); j++)
                          {
                            Mlt(DF1, phi(j), phi_t);
			    Mlt(1.0/Det(DF1), phi_t);
			    phi1(j) = DotProd(normale_fj, phi_t);                       
                          }
                        
                        // mise a jour matrice M et B
                        for (int iv = 0; iv < this->GetNbDof(); iv++)
                          for (int jv = 0; jv < this->GetNbDof(); jv++)
                            {
                              int i = IndexDof(iv);
                              int j = IndexDof(jv);
                              if ((i>= 0) && (j>= 0))
                                Mh(i, j) += poids*phi1(iv)*phi1(jv);
                            }
                        
                        for (int iv = 0; iv < this->GetNbDof(); iv++)
                          for (int jv = phi2.GetM()-nb_inside; jv < phi2.GetM(); jv++)
                            {
                              int i = IndexDof(iv);
                              int j = jv - phi2.GetM()+nb_inside;
                              if ((i>= 0) && (j>= 0))
                                Bh(i, j) += poids*phi1(iv)*phi2(jv);
                            }
                      }
                    
                    GetInverse(Mh);
                    Matrix<Real_wp> Ah(nb_inside, nb_inside);            
                    Ah.Fill(0);
                    Mlt(Mh, Bh, Ah);
                    
                    // DISP(num_loc); DISP(rot); DISP(Ah);
                    if (rot == 0)
                      {
                        bool test_identity = true;
                        for (int i = 0; i < nb_inside; i++)
                          for (int j = 0; j < nb_inside; j++)
                            {
                              if (j == i)
                                {
                                  if (abs(Ah(i, j) - 1.0) > threshold)
                                    test_identity = false;
                                }
                              else
                                {
                                  if (abs(Ah(i, j)) > threshold)
                                    test_identity = false;
                                }
                            }
                        
                        if (display_message)
                          if (!test_identity)
                            {
                              cout << "Warning : Operator not equal to identity for face "
                                   << num_loc << endl;
                              DISP(Ah);
                            }
                        
                      }
                    
                    Real_wp threshold = 1000*epsilon_machine;
                    // pour tous les ddls de la face
                    Vector<bool> DofUsed(nb_inside); DofUsed.Fill(false);            
                    for (int iv = 0; iv < this->GetNbDof(); iv++)
                      {
                        int i = IndexDof(iv);
                        if (i>= 0)
                          {
                            int nb_coef = 0;
                            for (int j = 0; j < nb_inside; j++)
                              if (abs(Ah(j, i)) > threshold)
                                nb_coef++;
                            
                            if (nb_coef == 0)
                              {
                                cout << "Cas impossible " << endl;
                                abort();
                              }
                            else if (nb_coef == 1)
                              {
                                // on est tombe sur un autre ddl, avec eventuellement un signe
                                for (int j = 0; j < nb_inside; j++)
                                  if (abs(Ah(j, i)) > threshold)
                                    {
                                      FaceRot(rot, i) = j;
                                      if ( abs(Ah(j, i) -1.0) <= threshold)
                                        OppositeSigneDof(rot, i) = false;
                                      else if ( abs(Ah(j, i) + 1.0) <= threshold)
                                        OppositeSigneDof(rot, i) = true;
                                      else
                                        {
                                          // this dof can't be handled with a sign
                                          linear_combination = true;
                                        }
                                      
                                      DofUsed(j) = true;
                                    }                       
                              }
                            else
                              {
                                // no direct relation between dofs
                                linear_combination = true;
                              }
                          }
                      }
                    
                    CoefCombination(rot) = Ah;
                    
                  }
                
                // DISP(CoefCombination);
                // DISP(linear_combination);
                // DISP(OppositeSigneDof);
                
                if (linear_combination)
                  {
                    if (quad)
                      nmap.SetFacesDofRotationQuad(this->GetOrder(), CoefCombination);
                    else
                      nmap.SetFacesDofRotationTri(this->GetOrder(), CoefCombination);
                  }
                else
                  {
                    if (quad)
                      {
                        nmap.SetFacesDofRotationQuad(this->GetOrder(), FaceRot);
                        nmap.SetSignDofRotationQuad(this->GetOrder(), OppositeSigneDof);
                      }
                    else
                      {
                        nmap.SetFacesDofRotationTri(this->GetOrder(), FaceRot);
                        nmap.SetSignDofRotationTri(this->GetOrder(), OppositeSigneDof);
                      }
                  }
                
              }
          }
      }
  }


  //! computation of values of u on nodal points of a face
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points of face num_loc
    \param[in] num_loc face number
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 3>
  ::ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    typename Vector2::value_type vx, vy, vz;
    for (int j = 0; j < Unode.GetM()/3; j++)
      {
        int node = FacesNode(j, num_loc);
        vx = 0; vy = 0; vz = 0;
        for (int k = 0; k < this->nb_dof_loc; k++)
          {
            vx += Un(k)*this->ValuePhi_Nodal(k, 3*node);
            vy += Un(k)*this->ValuePhi_Nodal(k, 3*node+1);
            vz += Un(k)*this->ValuePhi_Nodal(k, 3*node+2);
          }
        
        Unode(3*j) = vx; Unode(3*j+1) = vy; Unode(3*j+2) = vz;
      }
  }


  //! computation of divergence of u on nodal points of a face
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode divergence of u on nodal points of face num_loc
    \param[in] num_loc face number
  */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 3>::
  ComputeDivBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    for (int j = 0; j < Unode.GetM(); j++)
      {
        Unode(j) = 0;
        int node = FacesNode(j, num_loc);
        for (int k = 0; k < this->nb_dof_loc; k++)
          Unode(j) += Un(k)*this->DivPhi_Nodal(k, node);
      }
  }
  

  //! Integration against divergence of basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f div(\varphi_i) dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = div(\varphi_i)(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 3>::ApplyRhGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < nb_dof_loc; i++)
      {
        val = Real_wp(0);
        for (int j = 0; j < nb_points_quadrature_inside; j++)
          val += Div_Phi(i, j)*Vh(j);
	
        Uh(i) = val;
      }
  }


  //! computation of divergence of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh divergence of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j div phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = div phi_i(xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 3>::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    typename Vector2::value_type curl;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        curl = 0;
        for (int j = 0; j < nb_dof_loc; j++)
          curl += Uh(j)*Div_Phi(j, i);
        
        Vh(i) = curl;
      }
  }


  //! Integration against basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Ch Vh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension3, 3>::ApplyChGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < nb_dof_loc; i++)
      {
	val = 0;
	for (int j = 0; j < nb_points_quadrature_inside; j++)
	  val += Value_PhiVec(i, j)(0)*Vh(3*j)
            + Value_PhiVec(i, j)(1)*Vh(3*j+1)
            + Value_PhiVec(i, j)(2)*Vh(3*j+2);
	
	Uh(i) = val;
      }
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
  void ElementReference<Dimension3, 3>
  ::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    typename Vector2::value_type vx, vy, vz;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	vx = 0; vy = 0; vz = 0;
	for (int j = 0; j < nb_dof_loc; j++)
	  {
            vx += Value_PhiVec(j, i)(0)*Uh(j);
            vy += Value_PhiVec(j, i)(1)*Uh(j);
            vz += Value_PhiVec(j, i)(2)*Uh(j);
          }
	
	Vh(3*i) = vx; Vh(3*i+1) = vy; Vh(3*i+2) = vz;
      }
  }
  

  //! Integration against basis functions on boundaries
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] n boundary number
   */  
  template<class T>
  void ElementReference<Dimension3, 3>::
  ApplyShGen(const T& alpha, int n, const Vector<T>& feval, Vector<T>& res, int r) const
  {
    TinyVector<T, 3> f;
    for (int i = offset_faceSh(n); i < offset_faceSh(n+1); i++)
      {
        int p = i - offset_faceSh(n);
	int j = this->GetQuadNumber(n, p);
        CopyVector(feval, p, f);
        for (int k = 0; k < nb_dof_loc; k++)
          res(k) += alpha*DotProd(Value_PhiVec(k, j), f);
      }
  }


  //! Integration against divergence of basis functions on boundaries
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f div(\varphi_i) dx
    \param[in] n boundary number
   */  
  template<class T>
  void ElementReference<Dimension3, 3>::
  ApplyNablaShGen(const T& alpha, int n, const Vector<T>& feval, Vector<T>& res, int r) const
  {
    for (int i = offset_faceSh(n); i < offset_faceSh(n+1); i++)
      {
        int p = i - offset_faceSh(n);
	int j = this->GetQuadNumber(n, p);
        for (int k = 0; k < nb_dof_loc; k++)
          res(k) += alpha*Div_Phi(k, j)*feval(p);
      }
  }

  
  //! we add constant mass matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass coefficient
    \param[inout] A matrix modified
    A(m:, n:) = A(m:, n:) + mass M
    where M is the mass matrix, M_ij = \int \varphi_j \varphi_i dx
   */
  template<class T, class Prop>
  void ElementReference<Dimension3, 3>::
  AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& mass, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;

    for (int i = 0; i < nb_dof_loc; i++)
      {
	for (int j = 0; j < nb_dof_loc; j++)
	  val(j) = const_mass_matrix(0, 0)(i, j)*mass(0, 0)
	    + const_mass_matrix(0, 1)(i, j)*mass(0, 1)
	    + const_mass_matrix(1, 0)(i, j)*mass(1, 0)
	    + const_mass_matrix(2, 0)(i, j)*mass(2, 0)
	    + const_mass_matrix(0, 2)(i, j)*mass(0, 2)
	    + const_mass_matrix(1, 1)(i, j)*mass(1, 1) 
	    + const_mass_matrix(1, 2)(i, j)*mass(1, 2)
	    + const_mass_matrix(2, 1)(i, j)*mass(2, 1)
	    + const_mass_matrix(2, 2)(i, j)*mass(2, 2);
	
	A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
      }  
  }


  //! we add constant stiffness matrix C \int_K div phi_j div phi_i
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] C coefficients
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + C S
    where S is a stiffness matrix equal to :
    (S)_{i, j} = \int_K div(\phi_j) div(\phi_i)  dx
   */
  template<class T>
  void ElementReference<Dimension3, 3>::
  AddConstantStiffnessMatrixGen(int m, int n, const T& C, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;

    for (int i = 0; i < nb_dof_loc; i++)
      {
	for (int j = 0; j < nb_dof_loc; j++)
	  val(j) = const_div_matrix(i, j)*C;

	A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
      }
          
  }

   
  //! computes integrals against divergence of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension3, 3>
  ::ComputeIntegralDivRef(const Vector1 & feval, Vector2& res) const
  {
    this->ApplyRh(feval, res);
  }
  

  //! computes surface integrals against divergence of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension3, 3>
  ::ComputeIntegralSurfaceDivRef(const Vector1 & feval, Vector2& res, int n) const
  {
    typename Vector2::value_type one; SetComplexOne(one);
    res.Fill(0);
    this->ApplyNablaSh(one, n, feval, res);
  }

 
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void ElementReference<Dimension3, 3>
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    contrib = feval;
    SolveQR(MatProjectionDof, tauProjectionDof, contrib);
  }

  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void ElementReference<Dimension3, 3>
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    VectReal_wp xreal(feval.GetM()), ximag(feval.GetM());
    for (int i = 0; i < feval.GetM(); i++)
      {
        xreal(i) = real(feval(i));
        ximag(i) = imag(feval(i));
      }

    SolveQR(MatProjectionDof, tauProjectionDof, xreal);
    SolveQR(MatProjectionDof, tauProjectionDof, ximag);
    
    contrib.Reallocate(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      contrib(i) = Complex_wp(xreal(i), ximag(i));
  }
#endif
  
} // end namespace

#define MONTJOIE_FILE_VOLUME_REFERENCE_CXX
#endif
