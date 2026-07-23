#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{

  InitMontjoie(argc, argv);
  
  {
    // test of functions in Mesh_Base<Dimension2>
    Mesh<Dimension2> mesh;

    mesh.ReallocateVertices(5);
    mesh.Vertex(0).Init(-1.5, -1.5); mesh.Vertex(1) = R2(0, -1.5);
    mesh.Vertex(2).Init(1.5, 0); mesh.Vertex(3) = R2(1.5, 1.5);
    mesh.Vertex(4) = R2(0, 1.5);
    mesh.ResizeVertices(6);
    mesh.Vertex(5).Init(0.2, 0.2);
    cout << "Number of vertices in the mesh " << mesh.GetNbVertices() << endl;
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      cout << "Vertex " << i << " : " << mesh.Vertex(i) << endl;
    
    mesh.ReallocateElements(2);
    int ref = 1;
    mesh.Element(0).InitTriangular(0, 1, 5, ref);
    mesh.Element(1).InitTriangular(1, 2, 5, ref);
    mesh.ResizeElements(4);
    IVect num(3); num(0) = 2; num(1) = 4; num(2) = 5;
    mesh.Element(2).Init(num, ref);
    num(0) = 2; num(1) = 4; num(2) = 3;
    mesh.Element(3).Init(num, ref);
    
    cout << "Number of elements " << mesh.GetNbElt() << endl;
    cout << "Type of second element " << mesh.GetTypeElement(2) << endl;
    cout << "Basic mesh " << endl << mesh << endl;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      cout << "Element " << i << " : " << mesh.Element(i) << endl;
    
    // we change orientation of elements
    // and compute if elements are affine
    mesh.ReorientElements();
    
    cout << "Mesh after permutation " << endl << mesh << endl;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      cout << "Element " << i << " : " << mesh.Element(i) << endl;
    
    // third element affine ?
    cout << "Element 3 affine ? " << mesh.IsElementAffine(2) << endl;
        
    mesh.FindConnectivity();
    cout << "Number of edges " << mesh.GetNbEdges() << endl;
    // we print all the edges
    for (int i = 0; i < mesh.GetNbEdges(); i++)
      cout << "Edge " << i << " : " << mesh.GetEdge(i) << endl;
    
    // we display boundaries
    cout << "xmin = " << mesh.GetXmin() << endl;
    cout << "ymin = " << mesh.GetYmin() << endl;
    cout << "xmax = " << mesh.GetXmax() << endl;
    cout << "ymax = " << mesh.GetYmax() << endl;
        
    // we retrieve vertices of each element
    VectR2 s;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	mesh.GetVerticesElement(i, s);
	cout << "Vertices of element " << i << endl << s << endl;
      }
    
    VectReal_wp subdiv(3); subdiv.Fill(); Mlt(0.5, subdiv);    
    mesh.SubdivideMesh(subdiv);
    
    // changing the path where mesh is searched when a call to ConstructMesh is done
    Vector<string> param(1); param(0) = string("MAILLAGES/");
    mesh.SetInputData("MeshPath", param);
    
    // checking if two points are forming an edge
    int ne = -1;
    cout << "Point 2 and 5 on same edge ? " << mesh.PointsOnSameEdge(2, 5, ne) << endl;
    cout << "Point 0 and 4 on same edge ? " << mesh.PointsOnSameEdge(0, 4, ne) << endl;
    
    mesh.GetNeighboringElementsAroundVertex(2, num);
    cout << "Elements near vertex 2 : " << num << endl;
    
    mesh.AddBoundaryEdges();
    cout << "Number of referenced edges "<< mesh.GetNbBoundaryRef() << endl;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      cout << "Referenced Edge " << i << " : " << mesh.BoundaryRef(i) << endl;
    
    IVect ref_vertex;
    mesh.CreateReferenceVertices(ref_vertex);
    cout << "reference of vertices " << endl << ref_vertex << endl;
    
    Real_wp mesh_size = mesh.GetMeshSize();
    cout << "Mesh size : " << mesh_size << endl;
    
    // you can clear edges
    mesh.ClearConnectivity();
    cout << "Number of edges " << mesh.GetNbEdges() << endl;
    
    // and clear all the mesh
    mesh.Clear();
    
    // reload a new mesh
    mesh.Read("MAILLAGES/DisqueCond_QUAD.mesh");
    
    // we split it with metis
    IVect weight_elt(mesh.GetNbElt());
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	mesh.GetVerticesElement(i, s);
	R2 center;
	for (int j = 0; j < s.GetM(); j++)
	  center += s(j);
	
	Mlt(1.0/s.GetM(), center);
	if (center(0) < 0)
	  weight_elt(i) = 4;
	else
	  weight_elt(i) = 1;
      }
    
    int nb_parts = 4;
    Vector<IVect> num_elt; Vector<Mesh<Dimension2> > sub_mesh;

#ifdef MONTJOIE_WITH_SCOTCH
    mesh.SplitScotch(nb_parts, weight_elt, num_elt, sub_mesh, 1);
#else
    mesh.SplitMetis(nb_parts, weight_elt, num_elt, sub_mesh, 1);
#endif
    
    for (int j = 0; j < nb_parts; j++)
      sub_mesh(j).Write(string("subdomain")+to_str(j)+".mesh");
    
    // splitting into boxes
    mesh.SplitIntoBoxes(2, 2, 2, nb_parts, num_elt, sub_mesh, 1);
    for (int j = 0; j < nb_parts; j++)
      sub_mesh(j).Write(string("box")+to_str(j)+".mesh");

  }
  
  {
    // test of functions in NumberedMesh_Base<Dimension2> and NumberedMesh<Dimension2>
    Mesh<Dimension2> mesh;
    mesh.Read("example/mesh2d/test1.msh");
    
    mesh.Read("example/mesh2d/test2.msh");
    
    MeshNumbering<Dimension2> mesh_num(mesh);
    Vector<string> param(3);
    param(0) = string("MAX_EDGE"); param(1) = string("AUTO"); param(2) = to_str(1.0);
    mesh_num.SetInputData(string("OrderDiscretization"), param);
    
    mesh_num.ComputeVariableOrder(false);    

    TinyVector<IVect, 4> order_;
    mesh_num.GetOrderQuadrature(order_);
    mesh_num.number_map.ConstructQuadrature2D(order_);
    
    mesh_num.NumberMesh();
        
    mesh.WriteOrder("result1.msh", mesh_num);
    
    mesh.Read("example/mesh2d/test3.msh");
    VectReal_wp step(20); step.Zero();
    for (int i = 2; i < 20; i++)
      step(i) = Real_wp(i-1)/10;
    
    VectReal_wp coef(mesh.GetNbElt()); coef.Fill(1.0);
    // specifying a different velocity for domain 2
    for (int i = 0; i < mesh.GetNbElt(); i++)
      if (mesh.Element(i).GetReference() == 2)
	coef(i) = 2.0;
    
    mesh_num.SetMeshSizeVariableOrder(step);
    mesh_num.SetCoefficientVariableOrder(coef);
    mesh_num.ComputeVariableOrder(false);
    
    cout << "Variable order " << mesh_num.IsOrderVariable() << endl;

    mesh_num.GetOrderQuadrature(order_);
    mesh_num.number_map.ConstructQuadrature2D(order_);
    
    mesh_num.NumberMesh();
    cout << "Number of dofs " << mesh_num.GetNbDof() << endl;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      cout<<" Dofs of element " << i << endl << mesh_num.Element(i).GetNodle()<<endl;
    
    // rotation ?
    for (int i = 0; i < mesh.GetNbBoundary(); i++)
      if (mesh.Boundary(i).GetNbElements() == 2)
	{
	  int e1 = mesh.Boundary(i).numElement(0);
	  int e2 = mesh.Boundary(i).numElement(1);
	  int pos1 = mesh.Element(e1).GetPositionBoundary(i);
	  int rot = -1, pos2 = -1;
	  mesh_num.GetBoundaryRotation(i, e1, e2, pos1, pos2, rot);
	  cout <<" edge " << i <<" elements " << e1 <<" "<<e2
	       <<" local positions " << pos1 << "  " << pos2 << " rotation " << rot<<endl;
	}
    
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	cout << "nominal number of dofs on element " << i <<" "<<mesh_num.GetNbDofElement(i, mesh.Element(i))<<endl;
	cout << "current number of dofs on element " << i <<" "<<mesh_num.GetNbLocalDof(i)<<endl; 
      }
    
    mesh.WriteOrder("result2.msh", mesh_num);
  }
  
  {
    Mesh<Dimension2> mesh;
    // test of functions in MeshBoundaries_Base<Dimension2> and MeshBoundaries<Dimension2>
    mesh.Clear();
    MeshNumbering<Dimension2> mesh_num(mesh);
    mesh_num.SetOrder(3); mesh.SetGeometryOrder(3);  
    mesh.Read("example/mesh2d/test4.msh");
    
    // defining number of dofs for all orders
    mesh_num.number_map.SetFormulationDG(ElementReference_Base::CONTINUOUS);
    for (int r = 20; r >= 1; r--)
      {
	mesh_num.number_map.SetNbDofVertex(r, 1);
	mesh_num.number_map.SetNbDofEdge(r, r-1);
	mesh_num.number_map.SetNbDofQuadrangle(r, (r-1)*(r-1));
	mesh_num.number_map.SetNbDofTriangle(r, (r-1)*(r-2)/2);
	mesh_num.number_map.SetOppositeEdgesDofSymmetry(r, r-1);
      }

    TinyVector<IVect, 4> order_;
    mesh_num.GetOrderQuadrature(order_);
    mesh_num.number_map.ConstructQuadrature2D(order_);
    
    mesh_num.NumberMesh();
    
    mesh.ResizeNbReferences(14);
    cout << "Number of references " << mesh.GetNbReferences() << endl;
    mesh.SetCurveType(1, mesh.CURVE_CIRCLE);
    cout << "Type of curve for ref 1 : " << mesh.GetCurveType(1) << endl;
    mesh.SetBoundaryCondition(3, 1);
    cout << "Boundary condition for ref 3 : " << mesh.GetBoundaryCondition(3) << endl;
    Mesh<Dimension2> mesh2;
    mesh2.SetBoundaryCondition(mesh);
    cout << "Boundary condition for ref 3 : " << mesh2.GetBoundaryCondition(3) << endl;
    mesh.SetBodyNumber(1, 5);
    cout << "Body of curve 1  : " << mesh.GetBodyNumber(1) << endl;
    VectReal_wp rparam(3); rparam(0) = 0; rparam(1) = 0; rparam(2) = 1.0;
    mesh.SetCurveParameter(1, rparam); rparam.Clear();
    mesh.GetCurveParameter(1, rparam);
    cout << "Parameters for curve 1 : " << rparam << endl;
    mesh.ClearCurves();
    
    int newref = mesh.GetNewReference(); cout << " New reference " << newref << endl;
    mesh.SetNewReference(newref, 0, 3);
    cout << "Boundary condition for ref " << newref << " : " << mesh.GetBoundaryCondition(newref) << endl;

    mesh.ReallocatePmlAreas(1);
    PmlRegionParameter<Dimension2>& pml = mesh.GetPmlArea(0);

    pml.SetThicknessPML(1.5);
    cout << "Thickness for PML " << pml.GetThicknessPML() << endl;
    pml.SetAdditionPML(pml.PML_BOTH_SIDES, pml.PML_POSITIVE_SIDE, pml.PML_NEGATIVE_SIDE, 2);
    cout << "Number of layers " << pml.GetNbLayersPML() << endl;
    cout << "pml to add for y " << pml.GetTypeAdditionPML(1) << endl;
    
    // initialization of periodic conditions
    mesh.ClearPeriodicCondition();
    TinyVector<int, 2> num; num(0) = 2; num(1) = 4;
    mesh.AddPeriodicCondition(num, BoundaryConditionEnum::PERIODIC_CTE);
    num(0) = 3; num(1) = 5; mesh.AddPeriodicCondition(num, BoundaryConditionEnum::PERIODIC_CTE);
    mesh_num.TreatPeriodicCondition();
    cout << "number of periodic dofs " << mesh_num.GetNbPeriodicDof() << endl;
    for (int i = 0; i < mesh_num.GetNbPeriodicDof(); i++)
      cout << "dof " << mesh_num.GetPeriodicDof(i) << " translated from " << mesh_num.GetOriginalPeriodicDof(i)
	   <<" with vector " << mesh_num.GetTranslationPeriodicDof(i) << endl;
    
    MeshNumbering<Dimension2> mesh_num2(mesh_num);
    mesh_num2.number_map = mesh_num.number_map;
    mesh_num2.SetOrder(3); mesh2.SetGeometryOrder(3);
    mesh2.Read("example/mesh2d/test5.msh"); mesh_num2.NumberMesh();
    num(0) = 1; num(1) = 2; mesh2.AddPeriodicCondition(num, BoundaryConditionEnum::PERIODIC_THETA);
    mesh_num2.TreatPeriodicCondition();
    for (int i = 0; i < mesh_num2.GetNbPeriodicDof(); i++)
      cout << "dof " << mesh_num2.GetPeriodicDof(i) << " translated from " << mesh_num2.GetOriginalPeriodicDof(i)
	   <<" with vector " << mesh_num2.GetTranslationPeriodicDof(i) << endl;
    
    mesh2.SetCurveType(3, mesh2.CURVE_CIRCLE);
    mesh2.ProjectPointsOnCurves();
    for (int i = 0; i < mesh2.GetNbElt(); i++)
      cout << "safety coef for element " << i << " : " << mesh2.GetSafetyCoef(i) << endl;
    
    cout << "curve 1 and 2 equal ? " << mesh2.AreCurveEqual(1, 2) << endl;
    cout << "order of geometry " << mesh2.GetGeometryOrder() << endl;
    for (int i = 1; i < mesh2.GetGeometryOrder(); i++)
      {
	cout << "Inside point " << i << " : " << mesh2.GetNodalPointInsideEdge(i-1) << endl;
	cout << "Phi(0.23) = " << mesh2.GetInterpolate1D(i, 0.23) << endl;
      }
    
    mesh2.FindConnectivity();
    cout << "Number of referenced edges " << mesh2.GetNbBoundaryRef() << endl;
    for (int i = 0; i < mesh2.GetNbBoundaryRef(); i++)
      {
	int num_edge = i;
	cout << "Global edge number of referenced edge " << i << " : " <<  num_edge << endl;
	int i2 = num_edge;
	cout << "Referenced edge number of global edge " << num_edge << " : " << i2 << endl;
	
	// nodal points inside the edge
	VectR2 InsidePoint(mesh2.GetGeometryOrder()-1);
	for (int k = 0; k < InsidePoint.GetM(); k++)
	  InsidePoint(k) = mesh2.GetPointInsideEdge(i, k); 
	
	cout << "Extremities : " << mesh2.Vertex(mesh2.BoundaryRef(i).numVertex(0)) << "  "<< mesh2.Vertex(mesh2.BoundaryRef(i).numVertex(1))<<endl;
	cout << "Intemediary points on this edge " << InsidePoint <<endl;
      }
    
    // looping on edges of reference 1
    // choosing a first edge
    int iedge = -1;
    for (int i = 0; i < mesh2.GetNbBoundaryRef(); i++)
      if (mesh2.BoundaryRef(i).GetReference() == 1)
	iedge = i;
    
    int n1 = mesh2.BoundaryRef(iedge).numVertex(0);
    int n2 = mesh2.BoundaryRef(iedge).numVertex(1);
    bool test_loop = true; IVect ref_cond(mesh2.GetNbReferences()); ref_cond.Fill();
    int nfirst = n1;
    cout<<"Reference 1 with vertices " << endl << n1 << " ";
    while (test_loop)
      {
	cout << n2 << " ";
	int inext = mesh2.FindFollowingVertex(n1, n2, 1, ref_cond, iedge);	
	if ((inext == -1)||(inext == nfirst))
	  test_loop = false;
	else
	  {
	    n1 = n2; n2 = inext;
	  }
      }
    cout << endl;
    
    iedge = mesh2.FindEdgeRefWithExtremities(n1, n2);
    cout << "Referenced edge with vertices " << n1 << " " << n2 << " number : " << iedge << endl;
    
    n1 = mesh2.FindVertexNumber(R2(0, 4));
    cout << "vertex (0,4) of number " << n1 << endl;
    
    R2 middle = 0.5*(mesh2.Vertex(n1) + mesh2.Vertex(n2));
    bool point_on_edge = mesh2.IsPointOnBoundaryRef(iedge, middle);
    cout << "point " << middle <<" on the edge ? " << point_on_edge << endl;
    IVect all_ref;
    mesh2.GetAllReferencesAroundPoint(n1, all_ref);
    cout << "Point " << n1 << " surrounded by references " << all_ref << endl;
    
    int other_ref; n2 = mesh2.GetPointWithOtherReference(n1, all_ref(0), middle, other_ref);
    cout << " Point " << n2 << " of coordinates " << middle << " close to " << n1 << " but with another reference " << endl;
    
    DISP(mesh2.GetNbBoundaryRef());
    mesh2.RemoveReference(1); DISP(mesh2.GetNbBoundaryRef());
    mesh2.AddBoundaryEdges(); DISP(mesh2.GetNbBoundaryRef());
    mesh2.SortBoundariesRef(); 
    DISP(mesh2.GetNbBoundaryRef());
    mesh2.Write("result3.msh");
    
    mesh.ReadCurveType(all_ref, string("CIRCLE"));
    for (int i = 0; i < all_ref.GetM(); i++)
      cout << "Type of curve " << mesh.GetCurveType(all_ref(i)) << endl;

    mesh.ClearCurves();
    mesh.SetGeometryOrder(4);
    mesh.SetCurveType(1, mesh.CURVE_CIRCLE);
    mesh.Read("example/mesh2d/test6.msh");
    mesh.FindParametersCurve(); VectReal_wp param;
    mesh.GetCurveParameter(1, param); DISP(param);
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      if (mesh.BoundaryRef(i).GetReference() == 1)
	{
	  n1 = mesh.BoundaryRef(i).numVertex(0);
	  mesh.ProjectToCurve(mesh.Vertex(n1), 1);
	}
    
    mesh.ProjectPointsOnCurves();
    // finding middle points
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      if (mesh.BoundaryRef(i).GetReference() == 1)
	{
	  VectReal_wp xi(1); VectR2 pts_interm(1); xi(0) = 0.5;
	  mesh.GetPointsOnCurve(i, 1, xi, pts_interm);
	  cout << "Middle of edge " << i << " :  " << middle(0) << endl;
	}
    
    VectReal_wp coef_curve;
    mesh.FindParametersPlane(2, coef_curve);
    cout << "Coefficient of line " << coef_curve <<endl;

    mesh.Write("my_mesh.mesh");
    int num_elem = 4;
    const ElementGeomReference<Dimension2>& Fb = mesh.GetReferenceElement(num_elem);
    SetPoints<Dimension2> PointsElem;
    VectR2 s; mesh.GetVerticesElement(num_elem, s);
    Fb.FjElemNodal(s, PointsElem, mesh, num_elem);
    WriteElementMesh(mesh, Fb, PointsElem, "elem.mesh", 4);
    
    VectReal_wp subdiv(3); subdiv.Fill(); Mlt(0.5, subdiv);    
    mesh.SubdivideMesh(subdiv);
    mesh.Write("result4.msh");
    
    mesh.ClearCurves();
    mesh.SetGeometryOrder(4);
    mesh.SetCurveType(1, mesh.CURVE_SPLINE);
    mesh.Read("example/mesh2d/test7.msh");

    mesh.SubdivideMesh(subdiv);
    mesh.Write("result5.msh");

    mesh.ClearCurves();
    mesh.SetGeometryOrder(4);
    mesh.SetCurveType(3, mesh.CURVE_ELLIPSE);
    coef_curve.Reallocate(4); coef_curve(0) = 0; coef_curve(1) = 0;
    coef_curve(2) = 1.5; coef_curve(3) = 3.0;
    mesh.SetCurveParameter(3, coef_curve);
    mesh.Read("example/mesh2d/test8.msh");
    
    mesh.SubdivideMesh(subdiv);
    mesh.Write("result6.msh");

    mesh.ClearCurves();
    mesh.SetGeometryOrder(4);
    coef_curve.Reallocate(5); coef_curve(0) = 0; coef_curve(1) = 0;
    coef_curve(2) = 5.0; coef_curve(3) = 10.0; coef_curve(4) = 8.0;
    mesh.SetCurveType(1, mesh.CURVE_PEANUT);
    mesh.SetCurveParameter(1, coef_curve);
    mesh.Read("example/mesh2d/test9.msh");
    
    mesh.SubdivideMesh(subdiv);
    mesh.Write("result7.msh");

    // adding PML layers
    pml.AddPML(3, mesh);
    mesh.SetGeometryOrder(2);
    mesh.Write("result8.msh");
    
    mesh.Read("example/mesh2d/test10.msh");
    // extracting surface
    SurfacicMesh<Dimension2> mesh_surf;
    mesh.GetBoundaryMesh(1, mesh_surf, ref_cond);
    mesh_surf.Write("result9.msh");
    
    // revolution to get 3-D surfacic mesh
    Mesh<Dimension3> mesh3d; VectReal_wp AngleVertex;
    IVect InverseBoundaries, InverseVertex;
    mesh.GenerateSurfaceOfRevolution(1, mesh3d, InverseBoundaries, InverseVertex, AngleVertex, ref_cond);
    mesh3d.Write("result10.msh");
    
  }
  
  {
    // test of functions in Mesh<Dimension2>
    Mesh<Dimension2> mesh; mesh.SetGeometryOrder(4);
    mesh.Read("example/mesh2d/test11.msh");
    
    cout << "Number of referenced edges " << mesh.GetNbBoundaryRef() << endl;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      cout << "Referenced edge " << i << " : " << mesh.BoundaryRef(i) << endl;
    
    cout << "Number of global edges " << mesh.GetNbBoundary() << endl;
    for (int i = 0; i < mesh.GetNbBoundary(); i++)
      cout << "Global edge " << i << " : " << mesh.Boundary(i) << endl;
    
    cout << "Mesh containing triangles and quadrilaterals ? " << mesh.HybridMesh() << endl;
    cout << "Mesh containing only triangles ? " << mesh.IsOnlyTriangular() << endl;
    cout << "Mesh containing only quadrilaterals ? " << mesh.IsOnlyQuadrilateral() << endl;
    cout << "Number of triangles " << mesh.GetNbTriangles() << endl;
    cout << "Number of quadrangles " << mesh.GetNbQuadrangles() << endl;

    //VectReal_wp subdiv(3); subdiv.Fill(); Mlt(0.5, subdiv);    
    // mesh.SubdivideMesh(subdiv);    
        
    // changing the order of approximation
    mesh.SetGeometryOrder(3);
    mesh.Write("result11.msh");
    
    // retrieving nodes
    Vector<IVect> Nodle; VectR2 PosNodes;
    mesh.GetNodesCurvedMesh(PosNodes, Nodle, 3);
    cout << "Position of nodes "<<endl;
    for (int i = 0; i < PosNodes.GetM(); i++)
      cout << PosNodes(i) << endl;
    
    for (int i = 0; i < Nodle.GetM(); i++)
      cout << "Nodes of element " << i << endl << Nodle(i) <<endl;
    
    mesh.Write("result12.msh");
    
    mesh.RemoveDuplicateVertices();
    mesh.ForceCoherenceMesh();
    
    mesh.ApplyLocalRefinement(4, 3, 3.0);
    mesh.Write("result13.msh");
    
    R2 normale;
    mesh.SymmetrizeMesh(0, 0.0, normale);
    mesh.Write("result14.msh");
    
    mesh.ClearCurves(); mesh.SetGeometryOrder(2);
    mesh.Read("example/mesh2d/test12.msh");
    mesh.PeriodizeMeshTeta(R2(0,0));
    
    mesh.Write("result15.msh");
    
    mesh.SplitIntoTriangles();
    mesh.SplitIntoQuadrilaterals();
    
    mesh.Write("result16.msh");
    
    VectReal_wp step(4); step.Fill(); Mlt(1.0/(step.GetM()-1), step);
    mesh.SubdivideMesh(step);
    
    mesh.Write("result17.msh");
    
    mesh.Clear(); mesh.ClearCurves(); mesh.SetGeometryOrder(2);
    TinyVector<int, 4> ref_boundary; ref_boundary.Fill(2);
    TinyVector<int, 2> nbPoints; nbPoints.Fill(11);
    mesh.CreateRegularMesh(R2(-1,-1), R2(0, 1), nbPoints, 1, ref_boundary, mesh.QUADRILATERAL_MESH);
    
    mesh.SymmetrizeMeshToGetPositiveX();
    mesh.ExtrudeCoordinate(0, 2, 1.0, 0.5);
    mesh.ExtrudeCoordinate(1, 2, 1.0, 0.5);
    
    mesh.Write("result18.msh");
    
    mesh.ClearCurves(); mesh.SetGeometryOrder(2);
    Vector<string> param(3);
    param(0) = "example/mesh2d/test13.msh"; param(1) = "REFINED"; param(2) = "2";
    mesh.ConstructMesh(mesh.QUADRILATERAL_MESH, param);
    
    mesh.Write("result19.msh");
    
    Mesh<Dimension2> mesh2;
    mesh2.SetGeometryOrder(4);
    mesh2.Read("example/mesh2d/test14.msh");
    mesh.AppendMesh(mesh2);
    mesh.Write("result20.msh");
    
    // creating a sub-mesh
    VectR2 s; int nb_elt = 0; int nb_vertices = 0;
    Vector<bool> ElementOnSubdomain(mesh.GetNbElt()); ElementOnSubdomain.Fill(false);
    Vector<bool> VertexOnSubdomain(mesh.GetNbVertices()); VertexOnSubdomain.Fill(false);
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	mesh.GetVerticesElement(i, s);
	R2 center; 
	int nb_vert = mesh.Element(i).GetNbVertices();
	for (int k = 0; k < nb_vert; k++)
	  center += s(k);
	
	Mlt(Real_wp(1)/nb_vert, center);
	bool element_taken = (center(1) > -2.0);
	if (element_taken)
	  {
	    ElementOnSubdomain(i) = true;
	    nb_elt++;
	    for (int k = 0; k < nb_vert; k++)
	      {
		int num = mesh.Element(i).numVertex(k);
		if (!VertexOnSubdomain(num))
		  {
		    VertexOnSubdomain(num) = true;
		    nb_vertices++;
		  }
	      }
	  }
      }
    
    mesh2.Clear(); mesh2.ClearCurves();
    mesh.CreateSubmesh(mesh2, nb_vertices, nb_elt, VertexOnSubdomain, ElementOnSubdomain);
    mesh2.Write("result21.msh");
    
    IVect list_vertices, list_edges; Vector<bool> VertexOnAxe;
    mesh.GetEdgesOnZaxis(list_vertices, list_edges, VertexOnAxe);
    cout<<list_vertices<<endl;
    cout<<list_edges<<endl;
  }

  {
    Mesh<Dimension2> mesh;
    mesh.CreateRegularMesh(R2(-3, -2), R2(3, 2), TinyVector<int, 2>(5, 4), 1,
                           TinyVector<int, 4>(1, 2, 3, 4), mesh.QUADRILATERAL_MESH);

    mesh.Write("test.mesh");

    Real_wp teta = pi_wp/6;
    RotateMesh(mesh, R2(0, 0), teta);
    mesh.Write("rotat.mesh");
    
    R2 axis(-cos(teta), -sin(teta));
    Real_wp scal = DotProd(mesh.Vertex(0), axis);
    DISP(mesh.Vertex(4)); DISP(scal);
    Vector<int> ref_cond(mesh.GetNbReferences()+1); ref_cond.Zero();
    //ref_cond(2) = 1;
    ExtrudePMLLayer(mesh, axis, scal, ref_cond, -1, 0.5, 2, false);
    mesh.Write("pml.mesh");
  }
  
  return FinalizeMontjoie();  
}
