#ifndef MONTJOIE_FILE_FUNCTIONS_MESH_CXX

namespace Montjoie
{

#ifdef MONTJOIE_WITH_TWO_DIM
  //! translates the mesh
  void TranslateMesh(Mesh<Dimension2>& mesh, const R2& vec_u)
  {
    R2 ptA;
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      mesh.Vertex(i) += vec_u;
    
    int r = mesh.GetGeometryOrder();
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      for (int k = 0; k < r-1; k++)
        {
          ptA = mesh.GetPointInsideEdge(i, k);
          ptA += vec_u;
          mesh.SetPointInsideEdge(i, k, ptA);
        }
  }

  
  //! rotates the mesh
  void RotateMesh(Mesh<Dimension2>& mesh, const R2& center, const Real_wp& teta)
  {
    R2 ptA, ptB;
    Real_wp cos_teta = cos(teta);
    Real_wp sin_teta = sin(teta);
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      {
        ptA = mesh.Vertex(i);
        ptA -= center;
        ptB(0) = cos_teta*ptA(0) - sin_teta*ptA(1);
        ptB(1) = sin_teta*ptA(0) + cos_teta*ptA(1);
        ptB += center;
        mesh.Vertex(i) = ptB;
      }
    
    int r = mesh.GetGeometryOrder();
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      for (int k = 0; k < r-1; k++)
        {
          ptA = mesh.GetPointInsideEdge(i, k);
          ptA -= center;
          ptB(0) = cos_teta*ptA(0) - sin_teta*ptA(1);
          ptB(1) = sin_teta*ptA(0) + cos_teta*ptA(1);
          ptB += center;
          mesh.SetPointInsideEdge(i, k, ptB);
        }
  }

  
  //! applies a scaling on the mesh
  /*!
    \param[inout] mesh mesh to modify
    \param[in] vec_u scaling factor on each coordinate
   */
  void ScaleMesh(Mesh<Dimension2>& mesh, const R2& vec_u)
  {
    R2 ptA;
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      {
        mesh.Vertex(i)(0) *= vec_u(0);
        mesh.Vertex(i)(1) *= vec_u(1);
      }
    
    int r = mesh.GetGeometryOrder();
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      for (int k = 0; k < r-1; k++)
        {
          ptA = mesh.GetPointInsideEdge(i, k);
          ptA(0) *= vec_u(0);
          ptA(1) *= vec_u(1);
          mesh.SetPointInsideEdge(i, k, ptA);
        }
  }
#endif
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! translates a mesh
  void TranslateMesh(Mesh<Dimension3>& mesh, const R3& vec_u)
  {
    R3 ptA;
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      mesh.Vertex(i) += vec_u;
    
    int r = mesh.GetGeometryOrder();
    for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
      for (int k = 0; k < r-1; k++)
        {
          ptA = mesh.GetPointInsideEdge(i, k);
          ptA += vec_u;
          mesh.SetPointInsideEdge(i, k, ptA);
        }
    
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      for (int k = 0; k < mesh.GetNbPointsInsideFace(i); k++)
        {
          ptA = mesh.GetPointInsideFace(i, k);
          ptA += vec_u;
          mesh.SetPointInsideFace(i, k, ptA);
        }
    
    for (int ref = 0; ref < mesh.GetNbReferences(); ref++)
      if ((mesh.GetCurveType(ref) > 0) && (mesh.GetCurveType(ref) != mesh.CURVE_FILE))
        mesh.SetCurveType(ref, mesh.CURVE_FILE);
  }

  
  //! rotates a mesh
  void RotateMesh(Mesh<Dimension3>& mesh, const R3& center, const R3& axis_, const Real_wp& teta)
  {
    TinyVector<R3, 3> axis;
    axis(0) = axis_; 
    Mlt(1.0/Norm2(axis(0)), axis(0));
    
    GetVectorPlane(axis(0), axis(1), axis(2));
    Matrix3_3 matQ(axis), matRot, A, matTeta;
    
    matRot.Fill(0);
    matRot(0, 0) = 1.0;
    matRot(1, 1) = cos(teta);     matRot(1, 2) = -sin(teta);
    matRot(2, 2) = cos(teta);     matRot(2, 1) = sin(teta);
    
    MltTrans(matRot, matQ, A);
    Mlt(matQ, A, matTeta);
    // DISP(center); DISP(axis_); DISP(axis); DISP(matQ); DISP(matRot); DISP(matTeta);
    R3 ptA, ptB;
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      {
        ptA = mesh.Vertex(i);
        ptA -= center;
        Mlt(matTeta, ptA, ptB);
        ptB += center;
        mesh.Vertex(i) = ptB;
      }
    
    int r = mesh.GetGeometryOrder();
    for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
      for (int k = 0; k < r-1; k++)
        {
          ptA = mesh.GetPointInsideEdge(i, k);
          ptA -= center;
          Mlt(matTeta, ptA, ptB);
          ptB += center;
          mesh.SetPointInsideEdge(i, k, ptB);
        }
    
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      for (int k = 0; k < mesh.GetNbPointsInsideFace(i); k++)
        {
          ptA = mesh.GetPointInsideFace(i, k);
          ptA -= center;
          Mlt(matTeta, ptA, ptB);
          ptB += center;
          mesh.SetPointInsideFace(i, k, ptB);
        }

    for (int ref = 0; ref < mesh.GetNbReferences(); ref++)
      if ((mesh.GetCurveType(ref) > 0) && (mesh.GetCurveType(ref) != mesh.CURVE_FILE))
        mesh.SetCurveType(ref, mesh.CURVE_FILE);
  }


  //! applies a scaling on the mesh
  /*!
    \param[inout] mesh mesh to modify
    \param[in] vec_u scaling factor on each coordinate
   */
  void ScaleMesh(Mesh<Dimension3>& mesh, const R3& vec_u)
  {
    R3 ptA;
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      {
        mesh.Vertex(i)(0) *= vec_u(0);
        mesh.Vertex(i)(1) *= vec_u(1);
        mesh.Vertex(i)(2) *= vec_u(2);
      }
    
    int r = mesh.GetGeometryOrder();
    for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
      for (int k = 0; k < r-1; k++)
        {
          ptA = mesh.GetPointInsideEdge(i, k);
          ptA(0) *= vec_u(0);
          ptA(1) *= vec_u(1);
          ptA(2) *= vec_u(2);
          mesh.SetPointInsideEdge(i, k, ptA);
        }
    
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      for (int k = 0; k < mesh.GetNbPointsInsideFace(i); k++)
        {
          ptA = mesh.GetPointInsideFace(i, k);
          ptA(0) *= vec_u(0);
          ptA(1) *= vec_u(1);
          ptA(2) *= vec_u(2);
          mesh.SetPointInsideFace(i, k, ptA);
        }

    for (int ref = 0; ref < mesh.GetNbReferences(); ref++)
      if ((mesh.GetCurveType(ref) > 0) && (mesh.GetCurveType(ref) != mesh.CURVE_FILE))
        mesh.SetCurveType(ref, mesh.CURVE_FILE);
  }
  
  
  //! applies a scaling for each z
  void ScaleVariableMesh(Mesh<Dimension3>& mesh, ScalingMeshVirtualFunction3D& fct)
  {
    Real_wp coef_x, coef_y; R3 ptA;
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      {
        fct.EvaluateCoefficient(mesh.Vertex(i)(2), coef_x, coef_y);        
        mesh.Vertex(i)(0) *= coef_x;
        mesh.Vertex(i)(1) *= coef_y;
      }
    
    int r = mesh.GetGeometryOrder();
    for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
      for (int k = 0; k < r-1; k++)
        {
          ptA = mesh.GetPointInsideEdge(i, k);
          fct.EvaluateCoefficient(ptA(2), coef_x, coef_y);
          ptA(0) *= coef_x;
          ptA(1) *= coef_y;
          mesh.SetPointInsideEdge(i, k, ptA);
        }
    
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      for (int k = 0; k < mesh.GetNbPointsInsideFace(i); k++)
        {
          ptA = mesh.GetPointInsideFace(i, k);
          fct.EvaluateCoefficient(ptA(2), coef_x, coef_y);
          ptA(0) *= coef_x;
          ptA(1) *= coef_y;
          mesh.SetPointInsideFace(i, k, ptA);
        }

    for (int ref = 0; ref < mesh.GetNbReferences(); ref++)
      if ((mesh.GetCurveType(ref) > 0) && (mesh.GetCurveType(ref) != mesh.CURVE_FILE))
        mesh.SetCurveType(ref, mesh.CURVE_FILE);    
    
  }
#endif
  
  //! writes an element of the mesh on a file (the element is subdivided)
  template<class T>
  void WriteElementMesh(const Mesh<Dimension2, T>& mesh,
			const ElementGeomReference<Dimension2>& Fb,
			const SetPoints<Dimension2>& PointsElem,
                        const string& file_name, int i)
  {
    Mesh<Dimension2, T> new_mesh;
    new_mesh.SetGeometryOrder(Fb.GetGeometryOrder());
    
    int nv = mesh.Element(i).GetNbVertices();
    new_mesh.ReallocateVertices(nv);
    for (int j = 0; j < nv; j++)
      new_mesh.Vertex(j) = mesh.Vertex(mesh.Element(i).numVertex(j));
    
    new_mesh.ReallocateElements(1); IVect num(nv); num.Fill();
    new_mesh.Element(0).Init(num, 1);

    Vector<T> step_subdiv(Fb.GetGeometryOrder()+1);
    for (int j = 0; j < step_subdiv.GetM(); j++)
      step_subdiv(j) = Fb.PointsNodal1D(j);
    
    new_mesh.ReorientElements();
    
    new_mesh.FindConnectivity();
    new_mesh.AddBoundaryEdges();
    
    new_mesh.ConstructCrossReferenceBoundaryRef();

    new_mesh.SubdivideMesh(step_subdiv);
    
    R2 point;
    VectR2 s;
    mesh.GetVerticesElement(i, s);
    for (int j = 0; j < new_mesh.GetNbVertices(); j++)
      {
	Fb.FjLinear(s, Fb.PointsNodalND(j), point);
	for (int k = 0; k < new_mesh.GetNbVertices(); k++)
	  if (point == new_mesh.Vertex(k))
	    new_mesh.Vertex(k) = PointsElem.GetPointNodal(j);
      }
    
    new_mesh.AddBoundaryEdges();
    new_mesh.Write(file_name);
  }


  //! writes an element of the mesh on a file (the element is subdivided)
  template<class T>
  void WriteElementMesh(const Mesh<Dimension3, T>& mesh,
			const ElementGeomReference<Dimension3>& Fb,
			const SetPoints<Dimension3>& PointsElem,
                        const string& file_name, int i)
  {
    Mesh<Dimension3, T> new_mesh;
    int r = Fb.GetGeometryOrder();
    new_mesh.SetGeometryOrder(Fb.GetGeometryOrder());
    
    int nv = mesh.Element(i).GetNbVertices();
    new_mesh.ReallocateVertices(nv);
    Vector<int> permut(nv), inv_permut(nv), list_vertex(nv);
    for (int j = 0; j < nv; j++)
      list_vertex(j) = mesh.Element(i).numVertex(j);
    
    permut.Fill();
    Sort(list_vertex, permut);

    inv_permut.Fill();
    Vector<int> tmp(permut);
    Sort(tmp, inv_permut);

    Vector<int> InverseGlob(mesh.GetNbVertices());
    InverseGlob.Fill(-1);
    for (int j = 0; j < nv; j++)
      {
        new_mesh.Vertex(j) = mesh.Vertex(list_vertex(j));
        InverseGlob(list_vertex(j)) = j;
      }

    new_mesh.ResizeNbReferences(mesh.Element(i).GetNbBoundary()+3);
    for (int p = 1; p < new_mesh.GetNbReferences(); p++)
      new_mesh.SetCurveType(p, new_mesh.CURVE_FILE);
    
    new_mesh.ReallocateElements(1);
    new_mesh.Element(0).Init(inv_permut, 1);

    int type = mesh.Element(i).GetHybridType();
    int nb_edges_ref = 0;
    for (int j = 0; j < mesh.Element(i).GetNbEdges(); j++)
      {
        int ne = mesh.Element(i).numEdge(j);
        if (ne < mesh.GetNbEdgesRef())
          nb_edges_ref++;
      }
    
    new_mesh.ReallocateEdgesRef(nb_edges_ref);
    nb_edges_ref = 0;
    for (int j = 0; j < mesh.Element(i).GetNbEdges(); j++)
      {
        int ne = mesh.Element(i).numEdge(j);
        if (ne < mesh.GetNbEdgesRef())
          {
            int n1 = inv_permut(MeshNumbering<Dimension3, T>::FirstExtremityEdge(type, j));
            int n2 = inv_permut(MeshNumbering<Dimension3, T>::SecondExtremityEdge(type, j));
            cout << j << " Curved edge " << n1+1 << " " << n2+1 << endl;
            DISP(mesh.EdgeRef(ne).GetReference());
            new_mesh.EdgeRef(nb_edges_ref).Init(n1, n2, 1);
            for (int k = 0; k < r-1; k++)
              {
                new_mesh.SetPointInsideEdge(nb_edges_ref, k, mesh.GetPointInsideEdge(ne, k));
                DISP(Norm2(mesh.GetPointInsideEdge(ne, k)));
              }
            
            cout << endl;
            nb_edges_ref++;            
          }
      }
    
    cout << endl << endl;
    
    new_mesh.ReallocateBoundariesRef(mesh.Element(i).GetNbBoundary());
    for (int j = 0; j < mesh.Element(i).GetNbBoundary(); j++)
      {
        int nf = mesh.Element(i).numBoundary(j);
        int n1 = InverseGlob(mesh.Boundary(nf).numVertex(0));
        int n2 = InverseGlob(mesh.Boundary(nf).numVertex(1));
        int n3 = InverseGlob(mesh.Boundary(nf).numVertex(2)), n4 = -1;
        if (MeshNumbering<Dimension3, T>::IsLocalFaceQuadrilateral(type, j))
          {
            n4 = InverseGlob(mesh.Boundary(nf).numVertex(3));
            new_mesh.BoundaryRef(j).InitQuadrangular(n1, n2, n3, n4, j+2);
          }
        else
          new_mesh.BoundaryRef(j).InitTriangular(n1, n2, n3, j+2);

        if (mesh.IsBoundaryCurved(nf))
          {
            cout << j << " Curved face " << n1+1 << " " << n2+1 << " " << n3+1 << " " << n4+1 << endl;
            DISP(mesh.Boundary(nf).GetReference()); 
            for (int k = 0; k < mesh.GetPointInsideFace(nf).GetM(); k++)
              DISP(Norm2(mesh.GetPointInsideFace(nf)(k)));
            
            new_mesh.SetPointInsideFace(j, mesh.GetPointInsideFace(nf));
            cout << endl;
          }
        else
          new_mesh.SetCurveType(j+2, new_mesh.NO_CURVE);
      }
    
    int N = 21;
    Vector<T> step_subdiv(N);
    step_subdiv.Fill(); Mlt(Real_wp(1) / (N-1), step_subdiv);
    
    new_mesh.ReorientElements();
    new_mesh.FindConnectivity();
    new_mesh.ProjectPointsOnCurves();
    
    new_mesh.SubdivideMesh(step_subdiv);
    
    new_mesh.Write(file_name);
  }

  
  //! returns the integer associated with the condition number nom
  int BoundaryConditionEnum::GetInteger(string nom)
  {
    if (nom == "DIRICHLET")
      return LINE_DIRICHLET;
    else if (nom == "NEUMANN")
      return LINE_NEUMANN;    
    
    return LINE_INSIDE;
  }

  //! composition of periodicity types
  int BoundaryConditionEnum::GetCompositionPeriodicity(int t1, int t2)
  {
    if (t2 < 0)
      return t1;
    else
      {
        switch (t1)
          {
          case BoundaryConditionEnum::PERIODIC_CTE :
            return t2;
          case BoundaryConditionEnum::PERIODIC_THETA :
            if (t2 == BoundaryConditionEnum::PERIODIC_Z)
              return BoundaryConditionEnum::PERIODIC_ZTHETA;
            
            return BoundaryConditionEnum::PERIODIC_THETA;
          case BoundaryConditionEnum::PERIODIC_X :
            {
              switch (t2)
                {
                case BoundaryConditionEnum::PERIODIC_Y :
                  return BoundaryConditionEnum::PERIODIC_XY;
                case BoundaryConditionEnum::PERIODIC_Z :
                  return BoundaryConditionEnum::PERIODIC_XZ;
                }
              
              return BoundaryConditionEnum::PERIODIC_X;
            }
          case BoundaryConditionEnum::PERIODIC_Y :
            {
              switch (t2)
                {
                case BoundaryConditionEnum::PERIODIC_X :
                  return BoundaryConditionEnum::PERIODIC_XY;
                case BoundaryConditionEnum::PERIODIC_Z :
                  return BoundaryConditionEnum::PERIODIC_YZ;
                }
              
              return BoundaryConditionEnum::PERIODIC_Y;
            }
          case BoundaryConditionEnum::PERIODIC_Z :
            {
              switch (t2)
                {
                case BoundaryConditionEnum::PERIODIC_X :
                  return BoundaryConditionEnum::PERIODIC_XZ;
                case BoundaryConditionEnum::PERIODIC_Y :
                  return BoundaryConditionEnum::PERIODIC_YZ;
                }
              
              return BoundaryConditionEnum::PERIODIC_Z;
            }
            
          }
      }
    
    return -1;
  }
  

  //! composition of periodicity types  
  int BoundaryConditionEnum::GetCompositionPeriodicity(int t1, int t2, int t3)
  {
    if (t2 < 0)
      return t1;
    else
      {
        if (t3 < 0)
          return GetCompositionPeriodicity(t1, t2);
        else
          return BoundaryConditionEnum::PERIODIC_XYZ;
      }
    
    return -1;
  }

}

#define MONTJOIE_FILE_FUNCTIONS_MESH_CXX
#endif
