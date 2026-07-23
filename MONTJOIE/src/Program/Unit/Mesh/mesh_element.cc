#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

void CheckEdge(Edge<Dimension1>&)
{
  Edge<Dimension1> edge;

  int n1 = 2, n2 = 10;
  edge.Init(n1, n2);

  if ((edge.numVertex(0) != n1) || (edge.numVertex(1) != n2))
    {
      cout << "Init or numVertex incorrect" << endl;
      abort();
    }

  int ref = 3;
  edge.Init(n1, n2, ref);

  if ((edge.numVertex(0) != n1) || (edge.numVertex(1) != n2)
      || edge.GetReference() != ref)
    {
      cout << "Init incorrect" << endl;
      abort();
    }

  ref++;
  edge.SetReference(ref);
  if (edge.GetReference() != ref)
    {
      cout << "SetReference/GetReference incorrect" << endl;
      abort();
    }

  if (edge.GetNbVertices() != 2)
    {
      cout << "GetNbVertices incorrect" << endl;
      abort();
    }

}


void CheckEdge(Edge<Dimension2>&)
{
  Edge<Dimension2> edge;

  int n1 = 2, n2 = 10;
  edge.Init(n1, n2);

  if ((edge.numVertex(0) != n1) || (edge.numVertex(1) != n2))
    {
      cout << "Init or numVertex incorrect" << endl;
      abort();
    }

  int ref = 3;
  edge.Init(n1, n2, ref);

  if ((edge.numVertex(0) != n1) || (edge.numVertex(1) != n2)
      || edge.GetReference() != ref)
    {
      cout << "Init incorrect" << endl;
      abort();
    }

  ref++;
  edge.SetReference(ref);
  if (edge.GetReference() != ref)
    {
      cout << "SetReference/GetReference incorrect" << endl;
      abort();
    }

  if (edge.GetNbVertices() != 2)
    {
      cout << "GetNbVertices incorrect" << endl;
      abort();
    }

  if (edge.GetNbElements() != 0)
    {
      cout << "GetNbElements incorrect" << endl;
      abort();
    }

  edge.AddElement(2);
  if ((edge.GetNbElements() != 1) || (edge.numElement(0) != 2))
    {
      cout << "GetNbElements/numElement incorrect" << endl;
      abort();
    }

  edge.AddElement(5);
  if ((edge.GetNbElements() != 2) || (edge.numElement(0) != 2) || (edge.numElement(1) != 5))
    {
      cout << "GetNbElements/numElement incorrect" << endl;
      abort();
    }

  Edge<Dimension2> edge2;
  edge2.Init(3, 9, 1);
  edge2.SetElement(3, 6);
  if ((edge2.GetNbElements() != 2) || (edge2.numElement(0) != 3) || (edge2.numElement(1) != 6))
    {
      cout << "SetElement incorrect" << endl;
      abort();
    }
  
  if (edge.GetHybridType() != 0)
    {
      cout << "GetHybridType incorrect" << endl;
      abort();
    }

  edge.SetGeometryReference(4);
  if (edge.GetGeometryReference() != 4)
    {
      cout << "SetGeometryReference incorrect" << endl;
      abort();
    }

  edge.ClearConnectivity();
  if (edge.GetNbElements() != 0)
    {
      cout << "ClearConnectivity incorrect" << endl;
      abort();
    }
}


void CheckEdge(Edge<Dimension3>&)
{
  Edge<Dimension3> edge;

  int n1 = 2, n2 = 10;
  edge.Init(n1, n2);

  if ((edge.numVertex(0) != n1) || (edge.numVertex(1) != n2))
    {
      cout << "Init or numVertex incorrect" << endl;
      abort();
    }

  int ref = 3;
  edge.Init(n1, n2, ref);

  if ((edge.numVertex(0) != n1) || (edge.numVertex(1) != n2)
      || edge.GetReference() != ref)
    {
      cout << "Init incorrect" << endl;
      abort();
    }

  ref++;
  edge.SetReference(ref);
  if (edge.GetReference() != ref)
    {
      cout << "SetReference/GetReference incorrect" << endl;
      abort();
    }

  if (edge.GetNbVertices() != 2)
    {
      cout << "GetNbVertices incorrect" << endl;
      abort();
    }

  if ((edge.GetNbFaces() != 0) || (edge.GetNbElements() != 0))
    {
      cout << "GetNbFaces/GetNbElements incorrect" << endl;
      abort();
    }

  int nf1 = 12, nf2 = 4;
  int ne1 = 5, ne2 = 13, ne3 = 20;
  edge.AddFace(nf1); edge.AddFace(nf2);
  if ((edge.GetNbFaces() != 2) || (edge.numFace(0) != nf1) || (edge.numFace(1) != nf2))
    {
      cout << "GetNbFaces/numFace incorrect" << endl;
      abort();
    }

  edge.AddElement(ne1);
  edge.AddElement(ne2);
  edge.AddElement(ne3);
  if ((edge.GetNbFaces() != 2) || (edge.numFace(0) != nf1) || (edge.numFace(1) != nf2)
       || (edge.GetNbElements() != 3) || (edge.numElement(0) != ne1) || (edge.numElement(1) != ne2) || (edge.numElement(2) != ne3) )
    {
      cout << "GetNbElements/numElement incorrect" << endl;
      abort();
    }

  int nf3 = 22, nf4 = 15;
  edge.AddFace(nf3); edge.AddFace(nf4);
  if ((edge.GetNbFaces() != 4) || (edge.numFace(0) != nf1) || (edge.numFace(1) != nf2)
      || (edge.numFace(2) != nf3) || (edge.numFace(3) != nf4)
       || (edge.GetNbElements() != 3) || (edge.numElement(0) != ne1) || (edge.numElement(1) != ne2) || (edge.numElement(2) != ne3) )
    {
      cout << "GetNbElements/numElement incorrect" << endl;
      abort();
    }

  edge.ClearConnectivity();
  if ((edge.GetNbFaces() != 0) || (edge.GetNbElements() != 0))
    {
      cout << "ClearConnectivity incorrect" << endl;
      abort();
    }
  
}


void CheckFace(Face<Dimension2>&)
{
  Face<Dimension2> face;

  int n1 = 2, n2 = 17, n3 = 12, ref = 5;
  face.InitTriangular(n1, n2, n3, ref);

  if (!face.IsTriangular() || face.IsQuadrangular())
    {
      cout << "IsTriangular/IsQuadrangular incorrect" << endl;
      abort();
    }

  if ((face.GetHybridType() != 0) || (face.GetNbVertices() != 3) || (face.GetNbEdges() != 3))
    {
      cout << "GetHybridType/GetNbVertices/GetNbEdges incorrect" << endl;
      abort();
    }

  if ((face.numVertex(0) != n1) || (face.numVertex(1) != n2) || (face.numVertex(2) != n3)
      || (face.GetReference() != ref))
    {
      cout << "numVertex or GetReference incorrect" << endl;
      abort();
    }

  ref++;
  face.SetReference(ref);
  if (face.GetReference() != ref)
    {
      cout << "GetReference/SetReference incorrect" << endl;
      abort();
    }
  
  int ne1 = 6, ne2 = 13, ne3 = 20;
  face.SetEdge(0, ne1); face.SetEdge(1, ne2); face.SetEdge(2, ne3);

  if ((face.numEdge(0) != ne1) || (face.numEdge(1) != ne2) || (face.numEdge(2) != ne3))
    {
      cout << "numEdge/SetEdge incorrect" << endl;
      abort();
    }

  if (!face.GetOrientationEdge(0) || face.GetOrientationEdge(1) || face.GetOrientationEdge(2))
    {
      cout << "GetOrientationEdge incorrect" << endl;
      abort();
    }

  if ((face.GetPositionBoundary(ne1) != 0) || (face.GetPositionBoundary(ne2) != 1) ||
      (face.GetPositionBoundary(ne3) != 2) )
    {
      cout << "GetPositionBoundary incorrect" << endl;
      abort();
    }

  int n4 = 21;
  face.InitQuadrangular(n1, n2, n3, n4, ref);
    
  if (face.IsTriangular() || !face.IsQuadrangular())
    {
      cout << "IsTriangular/IsQuadrangular incorrect" << endl;
      abort();
    }

  if ((face.GetHybridType() != 1) || (face.GetNbVertices() != 4) || (face.GetNbEdges() != 4))
    {
      cout << "GetHybridType/GetNbVertices/GetNbEdges incorrect" << endl;
      abort();
    }

  if ((face.numVertex(0) != n1) || (face.numVertex(1) != n2) || (face.numVertex(2) != n3)
      || (face.numVertex(3) != n4) || (face.GetReference() != ref))
    {
      cout << "numVertex or GetReference incorrect" << endl;
      abort();
    }

  int ne4 = 32;
  face.SetEdge(0, ne1); face.SetEdge(1, ne2); face.SetEdge(2, ne3); face.SetEdge(3, ne4);

  if ((face.numEdge(0) != ne1) || (face.numEdge(1) != ne2) || (face.numEdge(2) != ne3) || (face.numEdge(3) != ne4))
    {
      cout << "numEdge/SetEdge incorrect" << endl;
      abort();
    }
  
  if (!face.GetOrientationEdge(0) || face.GetOrientationEdge(1) || !face.GetOrientationEdge(2) || face.GetOrientationEdge(3))
    {
      cout << "GetOrientationEdge incorrect" << endl;
      abort();
    }

  if ((face.GetPositionBoundary(ne1) != 0) || (face.GetPositionBoundary(ne2) != 1) ||
      (face.GetPositionBoundary(ne3) != 2) || (face.GetPositionBoundary(ne4) != 3) )
    {
      cout << "GetPositionBoundary incorrect" << endl;
      abort();
    }

  // specific to 2-D
  face.SetPML(0);
  if (!face.IsPML())
    {
      cout << "IsPML incorrect" << endl;
      abort();
    }

  if (face.GetTypePML() != 0)
    {
      cout << "GetTypePML incorrect" << endl;
      abort();
    }

  face.UnsetPML();
  if (face.IsPML())
    {
      cout << "UnsetPML incorrect" << endl;
      abort();      
    }

  face.SetCurved();
  if (!face.IsCurved())
    {
      cout << "SetCurved/IsCurved incorrect" << endl;
      abort();
    }

  face.UnsetCurved();
  if (face.IsCurved())
    {
      cout << "UnsetCurved incorrect" << endl;
      abort();
    }
  
  if ((face.GetNbBoundary() != 4) || (face.numBoundary(0) != ne1) || (face.numBoundary(1) != ne2) || (face.numBoundary(2) != ne3) || (face.numBoundary(3) != ne4))
    {
      cout << "numBoundary incorrect" << endl;
      abort();
    }
  
}


void CheckFace(Face<Dimension3>&)
{
  Face<Dimension3> face;

  int n1 = 2, n2 = 17, n3 = 12, ref = 5;
  face.InitTriangular(n1, n2, n3, ref);

  if (!face.IsTriangular() || face.IsQuadrangular())
    {
      cout << "IsTriangular/IsQuadrangular incorrect" << endl;
      abort();
    }

  if ((face.GetHybridType() != 0) || (face.GetNbVertices() != 3) || (face.GetNbEdges() != 3))
    {
      cout << "GetHybridType/GetNbVertices/GetNbEdges incorrect" << endl;
      abort();
    }

  if ((face.numVertex(0) != n1) || (face.numVertex(1) != n2) || (face.numVertex(2) != n3)
      || (face.GetReference() != ref))
    {
      cout << "numVertex or GetReference incorrect" << endl;
      abort();
    }

  ref++;
  face.SetReference(ref);
  if (face.GetReference() != ref)
    {
      cout << "GetReference/SetReference incorrect" << endl;
      abort();
    }
  
  int ne1 = 6, ne2 = 13, ne3 = 20;
  face.SetEdge(0, ne1); face.SetEdge(1, ne2); face.SetEdge(2, ne3);

  if ((face.numEdge(0) != ne1) || (face.numEdge(1) != ne2) || (face.numEdge(2) != ne3))
    {
      cout << "numEdge/SetEdge incorrect" << endl;
      abort();
    }

  if (!face.GetOrientationEdge(0) || face.GetOrientationEdge(1) || face.GetOrientationEdge(2))
    {
      cout << "GetOrientationEdge incorrect" << endl;
      abort();
    }

  if ((face.GetPositionBoundary(ne1) != 0) || (face.GetPositionBoundary(ne2) != 1) ||
      (face.GetPositionBoundary(ne3) != 2) )
    {
      cout << "GetPositionBoundary incorrect" << endl;
      abort();
    }

  int n4 = 21;
  face.InitQuadrangular(n1, n2, n3, n4, ref);
    
  if (face.IsTriangular() || !face.IsQuadrangular())
    {
      cout << "IsTriangular/IsQuadrangular incorrect" << endl;
      abort();
    }

  if ((face.GetHybridType() != 1) || (face.GetNbVertices() != 4) || (face.GetNbEdges() != 4))
    {
      cout << "GetHybridType/GetNbVertices/GetNbEdges incorrect" << endl;
      abort();
    }

  if ((face.numVertex(0) != n1) || (face.numVertex(1) != n2) || (face.numVertex(2) != n3)
      || (face.numVertex(3) != n4) || (face.GetReference() != ref))
    {
      cout << "numVertex or GetReference incorrect" << endl;
      abort();
    }

  int ne4 = 32;
  face.SetEdge(0, ne1); face.SetEdge(1, ne2); face.SetEdge(2, ne3); face.SetEdge(3, ne4);

  if ((face.numEdge(0) != ne1) || (face.numEdge(1) != ne2) || (face.numEdge(2) != ne3) || (face.numEdge(3) != ne4))
    {
      cout << "numEdge/SetEdge incorrect" << endl;
      abort();
    }
  
  if (!face.GetOrientationEdge(0) || face.GetOrientationEdge(1) || !face.GetOrientationEdge(2) || face.GetOrientationEdge(3))
    {
      cout << "GetOrientationEdge incorrect" << endl;
      abort();
    }

  if ((face.GetPositionBoundary(ne1) != 0) || (face.GetPositionBoundary(ne2) != 1) ||
      (face.GetPositionBoundary(ne3) != 2) || (face.GetPositionBoundary(ne4) != 3) )
    {
      cout << "GetPositionBoundary incorrect" << endl;
      abort();
    }

  if (face.GetNbElements() != 0)
    {
      cout << "GetNbElements incorrect" << endl;
      abort();
    }

  face.AddElement(2);
  if ((face.GetNbElements() != 1) || (face.numElement(0) != 2))
    {
      cout << "AddElement/numElement incorrect" << endl;
      abort();
    }
  
  face.AddElement(5);
  if ((face.GetNbElements() != 2) || (face.numElement(0) != 2) || (face.numElement(1) != 5))
    {
      cout << "AddElement/numElement incorrect" << endl;
      abort();
    }

  face.ClearConnectivity();
  
  if (face.GetNbElements() != 0)
    {
      cout << "ClearConnectivity incorrect" << endl;
      abort();
    }

  face.SetElement(3, 7);
  if ((face.GetNbElements() != 2) || (face.numElement(0) != 3) || (face.numElement(1) != 7))
    {
      cout << "SetElement incorrect" << endl;
      abort();
    }

  face.SetGeometryReference(11);
  if (face.GetGeometryReference() != 11)
    {
      cout << "GetGeometryReference/SetGeometryReference incorrect" << endl;
      abort();
    }

}

void CheckVolume()
{
  Volume elt;

  int n0 = 4, n1 = 12, n2 = 9, n3 = 15; int ref = 5;
  elt.InitTetrahedral(n0, n1, n2, n3, ref);

  if ((elt.numVertex(0) != n0) || (elt.numVertex(1) != n1) || (elt.numVertex(2) != n2)
      || (elt.numVertex(3) != n3) || (elt.GetReference() != ref) || (elt.GetHybridType() != 0)
      || (elt.GetNbVertices() != 4) || (elt.GetNbEdges() != 6) || (elt.GetNbFaces() != 4) || (elt.GetNbBoundary() != 4))
    {
      cout << "InitTetrahedral incorrect" << endl;
      abort();
    }

  ref++;
  elt.SetReference(ref);

  if (elt.GetReference() != ref)
    {
      cout << "GetReference/SetReference incorrect" << endl;
      abort();
    }

  elt.SetPML(0);
  if (!elt.IsPML())
    {
      cout << "SetPML/IsPML incorrect" << endl;
      abort();
    }

  if (elt.GetTypePML() != 0)
    {
      cout << "GetTypePML incorrect" << endl;
      abort();
    }

  elt.UnsetPML();
  if (elt.IsPML())
    {
      cout << "UnsetPML incorrect" << endl;
      abort();
    }

  elt.SetCurved();
  if (!elt.IsCurved())
    {
      cout << "SetCurved/IsCurved incorrect" << endl;
      abort();
    }

  elt.UnsetCurved();
  if (elt.IsCurved())
    {
      cout << "UnsetCurved incorrect" << endl;
      abort();
    }

  int ne0 = 5, ne1 = 15, ne2 = 21, ne3 = 36, ne4 = 8, ne5 = 13;
  elt.SetEdge(0, ne0); elt.SetEdge(1, ne1); elt.SetEdge(2, ne2); elt.SetEdge(3, ne3);
  elt.SetEdge(4, ne4); elt.SetEdge(5, ne5);
  if ((elt.numEdge(0) != ne0) || (elt.numEdge(1) != ne1) || (elt.numEdge(2) != ne2)
      || (elt.numEdge(3) != ne3) || (elt.numEdge(4) != ne4) || (elt.numEdge(5) != ne5))
    {
      cout << "SetEdge/numEdge incorrect" << endl;
      abort();
    }

  if ((elt.GetPositionEdge(ne0) != 0) || (elt.GetPositionEdge(ne1) != 1) || (elt.GetPositionEdge(ne2) != 2)
      || (elt.GetPositionEdge(ne3) != 3) || (elt.GetPositionEdge(ne4) != 4) || (elt.GetPositionEdge(ne5) != 5))
    {
      cout << "GetPositionEdge incorrect" << endl;
      abort();
    }

  int nf0 = 2, nf1 = 5, nf2 = 8, nf3 = 7;
  elt.SetFace(0, nf0); elt.SetFace(1, nf1); elt.SetFace(2, nf2); elt.SetFace(3, nf3);
  if ((elt.numFace(0) != nf0) || (elt.numFace(1) != nf1) || (elt.numFace(2) != nf2)
      || (elt.numFace(3) != nf3))
    {
      cout << "SetFace/numFace incorrect" << endl;
      abort();
    }

  int rot0 = 1, rot1 = 3, rot2 = 5, rot3 = 2;
  elt.SetOrientationFace(0, rot0); elt.SetOrientationFace(1, rot1);
  elt.SetOrientationFace(2, rot2); elt.SetOrientationFace(3, rot3);
  if ((elt.GetOrientationFace(0) != rot0) || (elt.GetOrientationFace(1) != rot1)
      || (elt.GetOrientationFace(2) != rot2) || (elt.GetOrientationFace(3) != rot3)
      || (elt.GetOrientationBoundary(0) != rot0) || (elt.GetOrientationBoundary(1) != rot1)
      || (elt.GetOrientationBoundary(2) != rot2) || (elt.GetOrientationBoundary(3) != rot3))
    {
      cout << "GetOrientationFace/SetOrientationFace incorrect" << endl;
      abort();
    }

  nf0++; nf1++; nf2++; nf3++; rot0++; rot1++; rot2++; rot3++;
  elt.SetFace(0, nf0, rot0); elt.SetFace(1, nf1, rot1);
  elt.SetFace(2, nf2, rot2); elt.SetFace(3, nf3, rot3);
  if ((elt.numFace(0) != nf0) || (elt.numFace(1) != nf1) || (elt.numFace(2) != nf2)
      || (elt.numFace(3) != nf3))
    {
      cout << "SetFace/numFace incorrect" << endl;
      abort();
    }

  if ((elt.GetOrientationFace(0) != rot0) || (elt.GetOrientationFace(1) != rot1)
      || (elt.GetOrientationFace(2) != rot2) || (elt.GetOrientationFace(3) != rot3)
      || (elt.GetOrientationBoundary(0) != rot0) || (elt.GetOrientationBoundary(1) != rot1)
      || (elt.GetOrientationBoundary(2) != rot2) || (elt.GetOrientationBoundary(3) != rot3))
    {
      cout << "GetOrientationFace/SetOrientationFace incorrect" << endl;
      abort();
    }

  n0 = 13; n1 = 3; n2 = 7; n3 = 11; int n4 = 9; ref--;
  elt.InitPyramidal(n0, n1, n2, n3, n4, ref);

  if ((elt.numVertex(0) != n0) || (elt.numVertex(1) != n1) || (elt.numVertex(2) != n2)
      || (elt.numVertex(3) != n3) || (elt.numVertex(4) != n4) || (elt.GetReference() != ref) || (elt.GetHybridType() != 1)
      || (elt.GetNbVertices() != 5) || (elt.GetNbEdges() != 8) || (elt.GetNbFaces() != 5) || (elt.GetNbBoundary() != 5))
    {
      cout << "InitPyramidal incorrect" << endl;
      abort();
    }

  n0 = 5; n1 = 1; n2 = 3; n3 = 8; n4 = 14; int n5 = 7; ref--;
  elt.InitWedge(n0, n1, n2, n3, n4, n5, ref);

  if ((elt.numVertex(0) != n0) || (elt.numVertex(1) != n1) || (elt.numVertex(2) != n2)
      || (elt.numVertex(3) != n3) || (elt.numVertex(4) != n4) || (elt.numVertex(5) != n5) || (elt.GetReference() != ref) || (elt.GetHybridType() != 2)
      || (elt.GetNbVertices() != 6) || (elt.GetNbEdges() != 9) || (elt.GetNbFaces() != 5) || (elt.GetNbBoundary() != 5))
    {
      cout << "InitWedge incorrect" << endl;
      abort();
    }

  n0 = 7; n1 = 2; n2 = 5; n3 = 12; n4 = 18; n5 = 16; int n6 = 21, n7 = 10; ref--;
  elt.InitHexahedral(n0, n1, n2, n3, n4, n5, n6, n7, ref);

  if ((elt.numVertex(0) != n0) || (elt.numVertex(1) != n1) || (elt.numVertex(2) != n2)
      || (elt.numVertex(3) != n3) || (elt.numVertex(4) != n4) || (elt.numVertex(5) != n5)
      || (elt.numVertex(6) != n6) || (elt.numVertex(7) != n7) || (elt.GetReference() != ref) || (elt.GetHybridType() != 3)
      || (elt.GetNbVertices() != 8) || (elt.GetNbEdges() != 12) || (elt.GetNbFaces() != 6) || (elt.GetNbBoundary() != 6))
    {
      cout << "InitPyramidal incorrect" << endl;
      abort();
    }
  
}

void CheckElementNumbering()
{
  ElementNumbering elt;

  if (elt.GetNbDof() != 0)
    {
      cout << "GetNbDof incorrect" << endl;
      abort();
    }

  int nb_dof = 8;
  Vector<int> num(8);
  
  elt.ReallocateDof(8);
  for (int i = 0; i < num.GetM(); i++)
    {
      num(i) = rand()%1000;
      elt.SetNumberDof(i, num(i));
    }

  if (elt.GetNbDof() != 8)
    {
      cout << "ReallocateDof incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < num.GetM(); i++)
    if ((elt.GetNumberDof(i) != num(i)) || (elt.GetNodle()(i) != num(i)))
      {
	cout << "GetNumberDof/GetNodle incorrect" << endl;
	abort();
      }

  num.PushBack(rand()%1000);
  elt.PushDof(num(8));

  if (elt.GetNbDof() != 9)
    {
      cout << "PushDof incorrect" << endl;
      abort();
    }

  for (int i = 0; i < num.GetM(); i++)
    if (elt.GetNumberDof(i) != num(i))
      {
	cout << "PushDof incorrect" << endl;
	abort();
      }

  Vector<int> num2(4);
  elt.AddDof(4);
  for (int i = 0; i < num2.GetM(); i++)
    {
      num2(i) = rand()%1000;
      elt.SetNumberDof(num.GetM() + i, num2(i));
    }

  num.PushBack(num2);
  if (elt.GetNbDof() != num.GetM())
    {
      cout << "AddDof incorrect" << endl;
      abort();
    }

  for (int i = 0; i < num.GetM(); i++)
    if (elt.GetNumberDof(i) != num(i))
      {
	cout << "AddDof/SetNumberDof incorrect" << endl;
	abort();
      }  
  
  Vector<int> neg(3);
  for (int i = 0; i < 3; i++)
    neg(i) = rand()%num.GetM();

  elt.SetNegativeDofNumber(neg.GetM(), neg);

  for (int i = 0; i < 3; i++)
    if (elt.GetNegativeDofNumber()(i) != neg(i))
      {
	cout << "SetNegativeDofNumber incorrect" << endl;
	abort();
      }

  elt.ReallocateFaces(4);
  int r1 = 3, r2 = 5, r3 = 6, r4 = 4;
  int off1 = 0, off2 = 9, off3 = 16, off4 = 32;
  elt.SetOrderFace(0, r1, off1);
  elt.SetOrderFace(1, r2, off2);
  elt.SetOrderFace(2, r3, off3);
  elt.SetOrderFace(3, r4, off4);
  if ( (elt.GetOrderFace(0) != r1) || (elt.GetOffsetFace(0) != off1) ||
       (elt.GetOrderFace(1) != r2) || (elt.GetOffsetFace(1) != off2) ||
       (elt.GetOrderFace(2) != r3) || (elt.GetOffsetFace(2) != off3) ||
       (elt.GetOrderFace(3) != r4) || (elt.GetOffsetFace(3) != off4) )
    {
      cout << "GetOrderFace/GetOffsetFace incorrect" << endl;
      abort();
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  {
    Edge<Dimension1> edge;
    CheckEdge(edge);
  }

  {
    Edge<Dimension2> edge;
    CheckEdge(edge);
  }

  {
    Edge<Dimension3> edge;
    CheckEdge(edge);
  }

  {
    Face<Dimension2> face;
    CheckFace(face);
  }

  {
    Face<Dimension3> face;
    CheckFace(face);
  }

  {
    CheckVolume();
  }

  {
    CheckElementNumbering();
  }

  return FinalizeMontjoie();
}
