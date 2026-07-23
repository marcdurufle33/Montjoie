#ifndef MONTJOIE_FILE_MESH_VERIF_CXX

template<class T>
bool EqualVector(const Vector<T>& x, const Vector<T>& y)
{
  if (x.GetM() != y.GetM())
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if (isnan(abs(x(i) - y(i))) || (abs(x(i) - y(i)) > R2::threshold))
      {
	DISP(i); DISP(x(i)); DISP(y(i));
	return false;
      }
  
  return true;
}

bool EqualEdge(const Edge<Dimension2>& x, const Edge<Dimension2>& y)
{
  if (x.numVertex(0) != y.numVertex(0))
    return false;

  if (x.numVertex(1) != y.numVertex(1))
    return false;

  if (x.GetReference() != y.GetReference())
    return false;

  if (x.GetNbElements() != y.GetNbElements())
    return false;

  for (int k = 0; k < x.GetNbElements(); k++)
    if (x.numElement(k) != y.numElement(k))
      return false;

  return true;
}

bool EqualEdge(const Edge<Dimension3>& x, const Edge<Dimension3>& y)
{
  if (x.numVertex(0) != y.numVertex(0))
    return false;

  if (x.numVertex(1) != y.numVertex(1))
    return false;

  if (x.GetReference() != y.GetReference())
    return false;

  if (x.GetNbElements() != y.GetNbElements())
    return false;

  for (int k = 0; k < x.GetNbElements(); k++)
    if (x.numElement(k) != y.numElement(k))
      return false;

  if (x.GetNbFaces() != y.GetNbFaces())
    return false;

  for (int k = 0; k < x.GetNbFaces(); k++)
    if (x.numFace(k) != y.numFace(k))
      return false;

  return true;
}

bool EqualFace(const Face<Dimension3>& x, const Face<Dimension3>& y)
{
  if (x.GetNbVertices() != y.GetNbVertices())
    return false;

  for (int k = 0; k < x.GetNbVertices(); k++)
    if (x.numVertex(k) != y.numVertex(k))
      return false;

  if (x.GetReference() != y.GetReference())
    return false;

  if (x.GetNbElements() != y.GetNbElements())
    return false;

  for (int k = 0; k < x.GetNbElements(); k++)
    if (x.numElement(k) != y.numElement(k))
      return false;

  if (x.GetNbEdges() != y.GetNbEdges())
    return false;

  for (int k = 0; k < x.GetNbEdges(); k++)
    if (x.numEdge(k) != y.numEdge(k))
      return false;

  return true;
}

R2 GetPointCurve(int ref, const Mesh<Dimension2>& mesh, const Real_wp& teta, const R2& pt)
{
  switch (mesh.GetCurveType(ref))
    {
    case Mesh<Dimension2>::CURVE_CIRCLE:
      {
	Real_wp xc = mesh.GetCurveParameter(ref)(0);
	Real_wp yc = mesh.GetCurveParameter(ref)(1);
	Real_wp rc = mesh.GetCurveParameter(ref)(2);
	return R2(xc + rc*cos(teta), yc + rc*sin(teta));
      }
      break;
    case Mesh<Dimension2>::CURVE_ELLIPSE:
      {
	Real_wp xc = mesh.GetCurveParameter(ref)(0);
	Real_wp yc = mesh.GetCurveParameter(ref)(1);
	Real_wp a = mesh.GetCurveParameter(ref)(2);
	Real_wp b = mesh.GetCurveParameter(ref)(3);
	return R2(xc + a*cos(teta), yc + b*sin(teta));
      }
      break;
    case Mesh<Dimension2>::CURVE_PEANUT:
      {
	Real_wp xc = mesh.GetCurveParameter(ref)(0);
	Real_wp yc = mesh.GetCurveParameter(ref)(1);
	Real_wp a = mesh.GetCurveParameter(ref)(2);
	Real_wp b = mesh.GetCurveParameter(ref)(3);
	Real_wp c = mesh.GetCurveParameter(ref)(4);
	return R2(xc + a*cos(teta)*sqrt(1.0+c*square(sin(teta))), yc + b*sin(teta));
      }
      break;
    default:
      return pt;
    }
}

bool IsMeshCorrect(const Mesh<Dimension2>& mesh, int order)
{
  const VectR2& val_teta = mesh.GetValueParameterOnBoundary();
  const Globatto<Real_wp>& lob = mesh.GetCurveFunctions1D();
  VectR2 points;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      Real_wp teta0 = val_teta(i)(0);
      Real_wp teta1 = val_teta(i)(1);
      int n0 = mesh.BoundaryRef(i).numVertex(0);
      int n1 = mesh.BoundaryRef(i).numVertex(1);
      R2 pt0 = mesh.Vertex(n0); R2 pt1 = mesh.Vertex(n1);

      // we check that value of parameters is correct
      if ( (GetPointCurve(ref, mesh, teta0, pt0) != pt0) || (GetPointCurve(ref, mesh, teta1, pt1) != pt1))
	{
	  DISP(i); DISP(ref);
	  DISP(GetPointCurve(ref, mesh, teta0, pt0)); DISP(pt0);
	  DISP(GetPointCurve(ref, mesh, teta1, pt1)); DISP(pt1);
	  return false;
	}

      // then we check that PointsEdgeRef is correct
      mesh.GetPointsOnCurve(i, ref, lob.Points(), points);
      for (int k = 1; k < lob.GetOrder(); k++)
	if (mesh.GetPointInsideEdge(i, k-1) != points(k))
	  {
	    DISP(mesh.GetPointInsideEdge(i, k-1)); DISP(points(k));
	    return false;
	  }
    }

  const EdgeLobatto& edge = mesh.GetSurfaceFiniteElement();
  if (!EqualVector(edge.Points, lob.Points()))
    {
      return false;
    }
  
  if (edge.GetOrder() != lob.GetOrder())
    {
      DISP(edge.GetOrder()); DISP(lob.GetOrder());
      return false;
    }

  if ((order != mesh.GetGeometryOrder()) || (lob.GetOrder() != order))
    {
      return false;
    }

  const TriangleGeomReference& tri = mesh.GetTriangleReferenceElement();
  if (tri.GetGeometryOrder() != order)
    {
      return false;
    }

  const QuadrangleGeomReference& quad = mesh.GetQuadrilateralReferenceElement();
  if (quad.GetGeometryOrder() != order)
    {
      return false;
    }

  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      if (!EqualEdge(mesh.BoundaryRef(i), mesh.Boundary(i)))
	{
	  return false;
	}
    }

  for (int i = 0; i < mesh.GetNbBoundary(); i++)
    {
      int ref = mesh.Boundary(i).GetReference();
      if (ref > 0)
	{
	  if (i >= mesh.GetNbBoundaryRef())
	    {
	      return false;
	    }
	}

      int n1 = mesh.Boundary(i).numVertex(0);
      int n2 = mesh.Boundary(i).numVertex(1);
      for (int k = 0; k < mesh.Boundary(i).GetNbElements(); k++)
	{
	  int ne = mesh.Boundary(i).numElement(k);
	  bool edge_found = false;
	  for (int j = 0; j < mesh.Element(ne).GetNbEdges(); j++)
	    {
	      int j0 = mesh.Element(ne).numVertex(j);
	      int j1 = mesh.Element(ne).numVertex((j+1)%mesh.Element(ne).GetNbEdges());
	      Sort(j0, j1);
	      if ((j0 == n1) && (j1 == n2))
		edge_found = true;
	    }

	  if (!edge_found)
	    {
	      return false;
	    }
	}
    }

  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      for (int j = 0; j < mesh.Element(i).GetNbEdges(); j++)
	{
	  int j0 = mesh.Element(i).numVertex(j);
	  int j1 = mesh.Element(i).numVertex((j+1)%mesh.Element(i).GetNbEdges());
	  Sort(j0, j1);
	  int ne = mesh.Element(i).numEdge(j);
	  if ((j0 != mesh.Boundary(ne).numVertex(0)) || (j1 != mesh.Boundary(ne).numVertex(1)))
	    {
	      return false;
	    }
	}

      if (mesh.Element(i).GetNbEdges() != mesh.GetReferenceElement(i).GetNbEdges())
	{
	  return false;
	}
    }
  
  return true;
}


bool IsMeshCorrect(const Mesh<Dimension3>& mesh, int order)
{
  const Globatto<Real_wp>& lob = mesh.GetCurveFunctions1D();
    
  /* const VectR2& val_teta = mesh.GetValueParameterOnBoundary();
  VectR2 points;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      Real_wp teta0 = val_teta(i)(0);
      Real_wp teta1 = val_teta(i)(1);
      int n0 = mesh.BoundaryRef(i).numVertex(0);
      int n1 = mesh.BoundaryRef(i).numVertex(1);
      R2 pt0 = mesh.Vertex(n0); R2 pt1 = mesh.Vertex(n1);

      // we check that value of parameters is correct
      if ( (GetPointCurve(ref, mesh, teta0, pt0) != pt0) || (GetPointCurve(ref, mesh, teta1, pt1) != pt1))
	{
	  DISP(i); DISP(ref);
	  DISP(GetPointCurve(ref, mesh, teta0, pt0)); DISP(pt0);
	  DISP(GetPointCurve(ref, mesh, teta1, pt1)); DISP(pt1);
	  return false;
	}

      // then we check that PointsEdgeRef is correct
      mesh.GetPointsOnCurve(i, ref, lob.Points(), points);
      for (int k = 1; k < lob.GetOrder(); k++)
	if (mesh.GetPointInsideEdge(i, k-1) != points(k))
	  {
	    DISP(mesh.GetPointInsideEdge(i, k-1)); DISP(points(k));
	    return false;
	  }
    }
  */

  if ((order != mesh.GetGeometryOrder()) || (lob.GetOrder() != order))
    {
      return false;
    }

  const TriangleGeomReference& tri = mesh.GetSurfaceFiniteElement();
  if (tri.GetGeometryOrder() != order)
    {
      return false;
    }

  const QuadrangleGeomReference& quad = mesh.GetSurfaceFiniteElement2();
  if (quad.GetGeometryOrder() != order)
    {
      return false;
    }

  for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
    {
      if (!EqualEdge(mesh.EdgeRef(i), mesh.GetEdge(i)))
	{
	  return false;
	}
    }

  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      if (!EqualFace(mesh.BoundaryRef(i), mesh.Boundary(i)))
	{
	  return false;
	}
    }

  for (int i = 0; i < mesh.GetNbEdges(); i++)
    {
      int ref = mesh.GetEdge(i).GetReference();
      if (ref > 0)
	{
	  if (i >= mesh.GetNbEdgesRef())
	    {
	      return false;
	    }
	}

      int n1 = mesh.GetEdge(i).numVertex(0);
      int n2 = mesh.GetEdge(i).numVertex(1);
      for (int k = 0; k < mesh.GetEdge(i).GetNbFaces(); k++)
	{
	  int ne = mesh.GetEdge(i).numFace(k);
	  bool edge_found = false;
	  for (int j = 0; j < mesh.Boundary(ne).GetNbEdges(); j++)
	    {
	      int j0 = mesh.Boundary(ne).numVertex(j);
	      int j1 = mesh.Boundary(ne).numVertex((j+1)%mesh.Boundary(ne).GetNbEdges());
	      Sort(j0, j1);
	      if ((j0 == n1) && (j1 == n2))
		edge_found = true;
	    }

	  if (!edge_found)
	    {
	      return false;
	    }
	}

      for (int k = 0; k < mesh.GetEdge(i).GetNbElements(); k++)
	{
	  int ne = mesh.Boundary(i).numElement(k);
	  int hyb = mesh.Element(ne).GetHybridType();
	  bool edge_found = false;
	  for (int j = 0; j < mesh.Element(ne).GetNbEdges(); j++)
	    {
	      int j0 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::FirstExtremityEdge(hyb, j));
	      int j1 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::SecondExtremityEdge(hyb, j));
	      Sort(j0, j1);
	      if ((j0 == n1) && (j1 == n2))
		edge_found = true;
	    }
  
	  if (!edge_found)
	    {
	      return false;
	    }
	}
    }

  for (int i = 0; i < mesh.GetNbBoundary(); i++)
    {
      int ref = mesh.Boundary(i).GetReference();
      if (ref > 0)
	{
	  if (i >= mesh.GetNbBoundaryRef())
	    {
	      return false;
	    }
	}

      bool quad = mesh.Boundary(i).IsQuadrangular();
      int n1 = mesh.Boundary(i).numVertex(0);
      int n2 = mesh.Boundary(i).numVertex(1);
      int n3 = mesh.Boundary(i).numVertex(2);
      int n4 = -1;
      if (quad)
	n4 = mesh.Boundary(i).numVertex(3);

      Sort(n1, n2, n3, n4);
      for (int j = 0; j < mesh.Boundary(i).GetNbEdges(); j++)
	{
	  int j0 = mesh.Boundary(i).numVertex(j);
	  int j1 = mesh.Boundary(i).numVertex((j+1)%mesh.Boundary(i).GetNbEdges());
	  Sort(j0, j1);
	  int ne = mesh.Boundary(i).numEdge(j);
	  if ((j0 != mesh.GetEdge(ne).numVertex(0)) || (j1 != mesh.GetEdge(ne).numVertex(1)))
	    {
	      return false;
	    }
	}
      
      for (int k = 0; k < mesh.Boundary(i).GetNbElements(); k++)
	{
	  int ne = mesh.Boundary(i).numElement(k);
	  int hyb = mesh.Element(ne).GetHybridType();
	  bool face_found = false;
	  for (int j = 0; j < mesh.Element(ne).GetNbFaces(); j++)
	    {
	      int j0 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::FirstExtremityFace(hyb, j));
	      int j1 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::SecondExtremityFace(hyb, j));
	      int j2 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::ThirdExtremityFace(hyb, j));
	      int j3 = -1;
	      if (quad)
		mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::FourthExtremityFace(hyb, j));

	      Sort(j0, j1, j2, j3);
	      if ((j0 == n1) && (j1 == n2) && (j2 == n3) && (j3 == n4))
		face_found = true;
	    }

	  if (!face_found)
	    {
	      return false;
	    }
	}
    }

  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      int hyb = mesh.Element(i).GetHybridType();
      for (int j = 0; j < mesh.Element(i).GetNbEdges(); j++)
	{
	  int j0 = mesh.Element(i).numVertex(MeshNumbering<Dimension3>::FirstExtremityEdge(hyb, j));
	  int j1 = mesh.Element(i).numVertex(MeshNumbering<Dimension3>::SecondExtremityEdge(hyb, j));
	  Sort(j0, j1);
	  int ne = mesh.Element(i).numEdge(j);
	  if ((j0 != mesh.GetEdge(ne).numVertex(0)) || (j1 != mesh.GetEdge(ne).numVertex(1)))
	    {
	      return false;
	    }
	}

      for (int j = 0; j < mesh.Element(i).GetNbFaces(); j++)
	{
	  int j0 = mesh.Element(i).numVertex(MeshNumbering<Dimension3>::FirstExtremityFace(hyb, j));
	  int j1 = mesh.Element(i).numVertex(MeshNumbering<Dimension3>::SecondExtremityFace(hyb, j));
	  int j2 = mesh.Element(i).numVertex(MeshNumbering<Dimension3>::ThirdExtremityFace(hyb, j));
	  bool quad = mesh.Boundary(i).IsQuadrangular();
	  int j3 = -1;
	  if (quad)
	    mesh.Element(i).numVertex(MeshNumbering<Dimension3>::FourthExtremityFace(hyb, j));

	  int nf = mesh.Element(i).numBoundary(j);
	  int n0 = mesh.Boundary(nf).numVertex(0), n1 = mesh.Boundary(nf).numVertex(1);
	  int n2 = mesh.Boundary(nf).numVertex(2), n3 = -1;
	  if (quad)
	    n3 = mesh.Boundary(nf).numVertex(3);
	  
	  Sort(n0, n1, n2, n3);
	  Sort(j0, j1, j2, j3);
	  if ((n0 != j0) || (n1 != j1) || (n2 != j2) || (n3 != j3))
	    {
	      return false;
	    }
	}
      
      if (mesh.Element(i).GetNbEdges() != mesh.GetReferenceElement(i).GetNbEdges())
	{
	  return false;
	}

      if (mesh.GetReferenceElement(i).GetGeometryOrder() != order)
	{
	  return false;
	}
    }
  
  return true;
}

template<class Dimension>
bool IsMeshEqual(const Mesh<Dimension>& mesh, const Mesh<Dimension>& mesh2)
{
  if (mesh.GetNbVertices() != mesh2.GetNbVertices())
    {
      cout << "Number of vertices different" << endl;
      return false;
    }
  
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    if (mesh.Vertex(i) != mesh2.Vertex(i))
      {
	DISP(i); DISP(mesh.Vertex(i)); DISP(mesh2.Vertex(i));
	return false;
      }
  
  if (mesh.GetNbBoundaryRef() != mesh2.GetNbBoundaryRef())
    {
      cout << "Number of faces different" << endl;
      return false;
    }
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    if ((mesh.BoundaryRef(i) != mesh2.BoundaryRef(i))
	|| (mesh.BoundaryRef(i).GetReference() != mesh2.BoundaryRef(i).GetReference()))
      {
	DISP(i); DISP(mesh.BoundaryRef(i)); DISP(mesh2.BoundaryRef(i));
	return false;
      }
  
  if (mesh.GetNbElt() != mesh2.GetNbElt())
    {
      cout << "Number of elements different" << endl;
      return false;
    }
  
  for (int i = 0; i < mesh.GetNbElt(); i++)
    if ((mesh.Element(i) != mesh2.Element(i))
	|| (mesh.Element(i).GetReference() != mesh2.Element(i).GetReference()))
      {
	DISP(mesh.Element(i)); DISP(mesh2.Element(i));
	return false;
      }
  
  return true;
}

#define MONTJOIE_FILE_MESH_VERIF_CXX
#endif
