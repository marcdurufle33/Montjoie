#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{

  InitMontjoie(argc, argv);

  if (true)
  {

    Mesh<Dimension3> mesh;

    // starting from an initial mesh
    mesh.Read("quart_sphere_hexa.mesh");
    
    // you have to create an extern boundary
    // such that you want to fill the space between the initial mesh and the boundary
    // the vertices of this extern boundary must match vertices of the initial mesh
    // one way is to extract the boundary of the initial mesh :
    int ref_boundary = 5; Vector<int> ref_cond(mesh.GetNbReferences()+1); ref_cond.Fill();
    SurfacicMesh<Dimension3> mesh_extern;
    mesh.GetBoundaryMesh(ref_boundary, mesh_extern, ref_cond);
    
    // we multiply coordinates by two to have a distinct surface
    for (int i = 0; i < mesh_extern.GetNbVertices(); i++)
      mesh_extern.Vertex(i) *= 2.0;
    
    // we want to add two layers
    IVect nb_intervals(2), ref_surf(3), ref_domain(2);
    nb_intervals(0) = 4; nb_intervals(1) = 3;
    ref_surf(0) = ref_boundary; ref_surf(1) = 6; ref_surf(2) = 7;
    ref_domain(0) = 3; ref_domain(1) = 4;
    
    // interface is located at mid-point => 0.5
    VectReal_wp step_subdiv(1); step_subdiv(0) = 0.5;
    mesh.ExtrudeOrigin(nb_intervals, mesh_extern, ref_surf, ref_domain, step_subdiv, false);
    
    /*
    int reference_sphere = 3;
    mesh.ReallocateVertices(7);
    mesh.Vertex(0).Init(0.0, 0.0, 0.0);
    mesh.Vertex(1).Init(1.0, 0.0, 0.0);
    mesh.Vertex(2).Init(-1.0, 0.0, 0.0);
    mesh.Vertex(3).Init(0.0, 1.0, 0.0);
    mesh.Vertex(4).Init(0.0, -1.0, 0.0);
    mesh.Vertex(5).Init(0.0, 0.0, 1.0);
    mesh.Vertex(6).Init(0.0, 0.0, -1.0);
    mesh.ReallocateElements(8); int ref_inside = 2;
    mesh.Element(0).InitTetrahedral(0, 1, 3, 5, ref_inside);
    mesh.Element(1).InitTetrahedral(0, 4, 1, 5, ref_inside);
    mesh.Element(2).InitTetrahedral(0, 4, 6, 1, ref_inside);
    mesh.Element(3).InitTetrahedral(0, 6, 3, 1, ref_inside);
    mesh.Element(4).InitTetrahedral(0, 3, 2, 5, ref_inside);
    mesh.Element(5).InitTetrahedral(0, 2, 4, 5, ref_inside);
    mesh.Element(6).InitTetrahedral(0, 6, 4, 2, ref_inside);
    mesh.Element(7).InitTetrahedral(0, 3, 6, 2, ref_inside);
    mesh.ReallocateBoundariesRef(8);
    mesh.BoundaryRef(0).InitTriangular(1, 3, 5, reference_sphere);
    mesh.BoundaryRef(1).InitTriangular(1, 6, 3, reference_sphere);
    mesh.BoundaryRef(2).InitTriangular(1, 4, 5, reference_sphere);
    mesh.BoundaryRef(3).InitTriangular(1, 6, 4, reference_sphere);
    mesh.BoundaryRef(4).InitTriangular(2, 3, 5, reference_sphere);
    mesh.BoundaryRef(5).InitTriangular(2, 5, 4, reference_sphere);
    mesh.BoundaryRef(6).InitTriangular(2, 6, 3, reference_sphere);
    mesh.BoundaryRef(7).InitTriangular(2, 4, 6, reference_sphere);
                
    mesh.ReorientElements();
    mesh.FindConnectivity();
    mesh.ProjectPointsOnCurves();
    
    mesh.CreateSphericalBall(1.0, 2.0, 4.0, reference_sphere,
                             true, mesh.HEXAHEDRAL_MESH);
*/
 mesh.Write("exemple.mesh"); 
    return 0;
  }
  
  if (false)
  {
    if (argc < 4)
      {
        abort();
      }

   
    string file_mesh(argv[1]);
    int nb_proc = to_num<int>(argv[2]);
    int ref = to_num<int>(argv[3]);
    
    Mesh<Dimension3> mesh;
    mesh.Read(file_mesh);
    
    IVect weight(mesh.GetNbElt());
    weight.Fill(1);
    Vector<IVect> NumElem;
    Vector<Mesh<Dimension3> > sub_mesh;
    
    mesh.SplitScotch(nb_proc, weight, NumElem, sub_mesh, 0);
    
    int sum = 0, nb_faces = 0;
    for (int i = 0; i < nb_proc; i++)
      {
        sub_mesh(i).Write(string("subdomain") + to_str(i) + ".mesh");
        nb_faces = 0;
        for (int j = 0; j < sub_mesh(i).GetNbBoundaryRef(); j++)
          if (sub_mesh(i).BoundaryRef(j).GetReference() == ref)
            nb_faces++;
        
        sum += nb_faces;
        DISP(i); DISP(nb_faces);
      }
    
    DISP(sum); 
    nb_faces = 0;
    for (int j = 0; j < mesh.GetNbBoundaryRef(); j++)
      if (mesh.BoundaryRef(j).GetReference() == ref)
        nb_faces++;
    
    DISP(nb_faces);
  }
  
  //return 0;

  {
    // testing cylinder
    Mesh<Dimension3> mesh ;
    mesh.SetGeometryOrder(4);
    mesh.print_level = 6;
    mesh.ClearCurves();    
    mesh.SetCurveType(1, mesh.CURVE_CYLINDER);
    VectReal_wp param_cyl(7);
    param_cyl.Fill(0);
    param_cyl(5) = 1.0; param_cyl(6.0) = 1.0;
    mesh.SetCurveParameter(1, param_cyl);
    
    mesh.Clear();
    mesh.ReallocateVertices(8);
    mesh.Vertex(0).Init(0.2, 0, 0);
    mesh.Vertex(1).Init(0, 0.2, 0);
    mesh.Vertex(2).Init(0.2, 0, 1.0);
    mesh.Vertex(3).Init(0, 0.2, 1.0);

    mesh.Vertex(4).Init(1.0, 0, 0);
    mesh.Vertex(5).Init(0, 1.0, 0);
    mesh.Vertex(6).Init(1.0, 0, 1.0);
    mesh.Vertex(7).Init(0, 1.0, 1.0);
    
    mesh.ReallocateElements(1);
    mesh.Element(0).InitHexahedral(0, 1, 3, 2, 4, 5, 7, 6, 1);
    
    mesh.ReallocateBoundariesRef(1);
    mesh.BoundaryRef(0).InitQuadrangular(4, 5, 7, 6, 1);
    
    mesh.ReorientElements();
    mesh.FindConnectivity();
    mesh.AddBoundaryEdges();
    mesh.FindConnectivity();
    mesh.ProjectPointsOnCurves();
    mesh.Write("init.mesh");

    VectReal_wp step(5); step.Fill(); Mlt(0.25, step);
    mesh.SubdivideMesh(step);    
    mesh.Write("cylinder.mesh");        

  }
  
  {    
    // test of functions in MeshBoundaries<Dimension3>
    Mesh<Dimension3> mesh;
    
    mesh.SetGeometryOrder(4);
    mesh.Read("test3d/test1.msh");
    
    mesh.RemoveReference(3);
    
    mesh.Write("result1.msh");
    
    Vector<Real_wp> step(5); step.Fill(); Mlt(1.0/4, step);
    mesh.SubdivideMesh(step);
    
    mesh.Write("result2.msh");
    
    mesh.SetGeometryOrder(1);
    mesh.Read("test3d/test2.msh");
    
    mesh.AddBoundaryFaces();
    
    mesh.Write("result3.mesh");    
    
    mesh.RedistributeReferences();
    
    mesh.Write("result4.mesh");
    
    mesh.Read("test3d/test3.msh");
    
    IVect num;
    mesh.GetNeighboringBoundariesAroundVertex(0, num);
    
    cout << "Faces around vertex 0 : " << num << endl;
    
    // parameters of a sphere ?
    mesh.ClearCurves();
    mesh.SetCurveType(1, mesh.CURVE_SPHERE);
    mesh.Read("test3d/test4.msh");
    
    VectReal_wp param;
    mesh.GetCurveParameter(1, param);
    cout << "Parameters of the sphere : " << param << endl;

    mesh.SubdivideMesh(step);    
    mesh.Write("result5.msh");
    
    // parameters of a cylinder
    mesh.ClearCurves();
    mesh.SetCurveType(1, mesh.CURVE_CYLINDER);
    mesh.Read("test3d/test5.msh");
    
    mesh.GetCurveParameter(1, param);
    cout << "Parameters of the cylinder : " << param << endl;
    
    mesh.SubdivideMesh(step);    
    mesh.Write("result6.msh");
    
    // parameters of a conic cylinder
    mesh.ClearCurves();
    mesh.SetCurveType(1, mesh.CURVE_CONIC);
    mesh.SetCurveType(2, mesh.CURVE_SPHERE);
    mesh.Read("test3d/test6.msh");
    
    mesh.GetCurveParameter(1, param);
    cout << "Parameters of the cone : " << param << endl;
    
    mesh.SubdivideMesh(step);
    mesh.Write("result7.msh");
    
    // adding PML layers to a tetrahedral mesh, then hexahedral
    mesh.ClearCurves(); mesh.SetGeometryOrder(2);
    mesh.Read("test3d/test7.msh");

    mesh.ReallocatePmlAreas(1);
    PmlRegionParameter<Dimension3>& pml = mesh.GetPmlArea(0);
    
    pml.SetThicknessPML(3.0);
    pml.SetAdditionPML(pml.PML_BOTH_SIDES, pml.PML_BOTH_SIDES, pml.PML_BOTH_SIDES, 1);
    pml.AddPML(3, mesh);
    
    step.Reallocate(4); step.Fill(); Mlt(1.0/3, step);
    mesh.SubdivideMesh(step);
    mesh.Write("result8.mesh");    
    
    mesh.ClearCurves(); mesh.SetGeometryOrder(2);
    mesh.Read("test3d/test8.msh");

    mesh.ReallocatePmlAreas(1);
    PmlRegionParameter<Dimension3>& pml2 = mesh.GetPmlArea(0);

    pml2.SetThicknessPML(1.0);
    pml2.SetAdditionPML(pml2.PML_BOTH_SIDES, pml2.PML_NO, pml2.PML_NO, 1);
    pml2.AddPML(2, mesh);
    
    mesh.SubdivideMesh(step);
    mesh.Write("result9.mesh");
    
    SurfacicMesh<Dimension3> mesh_surf;
    IVect ref_cond(10); ref_cond.Fill(1); ref_cond(0) = 0; ref_cond(3) = 0;
    mesh.GetBoundaryMesh(1, mesh_surf, ref_cond);
    
    mesh_surf.Write("result10.msh");
    
  }

  {    
    // test of functions in Mesh<Dimension3>
    Mesh<Dimension3> mesh;
    
    mesh.SetGeometryOrder(4);    
    mesh.Read("test3d/test10.msh");
    
    VectReal_wp step(3); step.Fill(); Mlt(1.0/2, step);
    
    // we change the order of approximation for the geometry
    mesh.SetGeometryOrder(3);
    mesh.Write("result11.msh");
    
    // you can apply Warburton algorithm, but call after FindConnectivity
    mesh.ChangeLocalNumberingMesh();
    mesh.FindConnectivity();
    mesh.Write("result12.msh");
    
    mesh.Read("test3d/test11.neu");
    // you can call ReorientElements to force positive jacobian
    // but this method is already called whenever it is necessary in other methods
    mesh.ReorientElements();
    mesh.FindConnectivity();
    
    mesh.SubdivideMesh(step);
    
    mesh.Write("result13.vtk");
    
    // you can reorient faces in the specific case of a planar surface
    mesh.Read("test3d/test12.msh");
    
    IVect ref_cond(mesh.GetNbReferences()+1); ref_cond.Fill(1);
    mesh.ReorientFaces(1, ref_cond, R3(0,0,1));
    mesh.Write("result14.msh");
    
    mesh.SetGeometryOrder(4);
    mesh.Read("test3d/test13.msh");
    
    // symmetry with plane x = 0
    mesh.SymmetrizeMesh(0, 0.0, R3(0, 0, 0));
    // followed by symmetry with plane y = 0
    mesh.SymmetrizeMesh(1, 0.0, R3(0, 0, 0));
    // and finally with plane z = 0
    mesh.SymmetrizeMesh(2, 0.0, R3(0, 0, 0));

    mesh.SubdivideMesh(step);    
    mesh.Write("result15.mesh");
    
    //mesh.Read("test3d/test14.msh");
    // you can remove duplicate vertices
    //mesh.RemoveDuplicateVertices();
    //mesh.Write("result16.msh");

    mesh.SetGeometryOrder(4);
    mesh.Read("test3d/test15.msh");
    // if your mesh contains vertices/edges/faces not connected to a volumic element, it is dropped
    mesh.ForceCoherenceMesh();

    mesh.Write("result17.msh");
    
    //VectReal_wp step;
    step.Reallocate(3); step.Fill(); Mlt(1.0/2, step);
    mesh.SplitIntoHexahedra();
    
    //mesh.SubdivideMesh(step);
    
    mesh.Write("result18.mesh");
    

    mesh.Read("test3d/test16.neu");
    mesh.SplitIntoTetrahedra();
    
    //mesh.SubdivideMesh(step);

    mesh.Write("result19.mesh");
    
    TinyVector<int, 3> nb_points; nb_points.Fill(3);
    TinyVector<int, 6> ref_boundary; ref_boundary.Fill(3);
    mesh.CreateRegularMesh(R3(-2, -2, -2), R3(2, 2, 2), nb_points,
                           1, ref_boundary, mesh.HYBRID_MESH);
    
    step.Reallocate(4); step.Fill(); Mlt(1.0/3, step);
    mesh.SubdivideMesh(step);    
    mesh.Write("result20.mesh");
    
    mesh.ClearCurves();
    mesh.SetCurveType(1, mesh.CURVE_SPHERE);
    VectReal_wp step_sphere(1); step_sphere(0) = 1.5;
    mesh.CreateSphericalCrown(1.0, 2.0, 2.0, true, false, step_sphere);
    DISP(mesh.GetCurveParameter(1));
    
    mesh.SubdivideMesh(step);    
    mesh.Write("result21.mesh");
    
    mesh.ClearCurves(); mesh.SetGeometryOrder(4);
    mesh.SetCurveType(1, mesh.CURVE_SPHERE);
    mesh.SetCurveType(3, mesh.CURVE_SPHERE);
    Vector<string> param(10);
    param(0) = "BALL"; param(1) = "2.0"; param(2) = "3.0";
    param(3) = "2.0"; param(4) = "1"; param(5) = "REGULAR";
    mesh.ConstructMesh(mesh.HYBRID_MESH, param);
    
    mesh.SubdivideMesh(step);    
    mesh.Write("result22.mesh");
    
    mesh.ClearCurves(); mesh.SetGeometryOrder(4);
    mesh.Read("test3d/test17.msh");
    mesh.ExtrudeCoordinate(2, 3, 4.0, 1.5);
    
    mesh.SubdivideMesh(step);    
    mesh.Write("result23.mesh");
    
    mesh.Read("test3d/test18.msh");
    Mesh<Dimension3> add_mesh; add_mesh.SetGeometryOrder(4);
    add_mesh.Read("test3d/test19.msh");
    mesh.AppendMesh(add_mesh);
    
    mesh.SubdivideMesh(step);    
    mesh.Write("result24.mesh");
    
    mesh.SetGeometryOrder(2); add_mesh.SetGeometryOrder(2);
    mesh.Read("test3d/test20.msh");
    add_mesh.Read("test3d/test21.msh");
    IVect nb_interv(1), ref_surf(2), ref_domain(1);
    nb_interv(0) = 1;
    ref_surf(0) = 1; ref_surf(1) = 3; ref_domain(0) = 1;
    VectReal_wp step_layer;
    mesh.ExtrudeOrigin(nb_interv, add_mesh, ref_surf, ref_domain, step_layer);
    
    mesh.SubdivideMesh(step);    
    mesh.Write("result25.msh");
    
    mesh.SetGeometryOrder(4); add_mesh.SetGeometryOrder(4);
    Mesh<Dimension2> mesh2d; mesh2d.SetGeometryOrder(5);
    mesh2d.Read("test3d/test22.msh");

    nb_interv.Reallocate(4);
    step_layer.Reallocate(5);
    nb_interv(0) = 2; nb_interv(1) = 1; nb_interv(2) = 4; nb_interv(3) = 2;
    step_layer(0) = 0.0; step_layer(1) = 0.1;
    step_layer(2) = 0.14; step_layer(3) = 0.5; step_layer(4) = 0.7;
    R3 vec_e1(1,0,0), vec_e2(0,1,0), vec_e3(0,0,1);
    Vector<IVect> ref_surfacic(8), ref_volume(4);
    for (int i = 0; i < 4; i++)
      {
        ref_surfacic(2*i).Reallocate(4); ref_surfacic(2*i).Fill(0);
        ref_surfacic(2*i+1).Reallocate(4); ref_surfacic(2*i+1).Fill(0);
        ref_volume(i).Reallocate(4); ref_volume(i).Fill(0);
      }
    ref_surfacic(0)(1) = 4;
    ref_surfacic(1)(2) = 5; ref_surfacic(1)(3) = 5; 
    ref_surfacic(2)(1) = 4; ref_surfacic(2)(2) = 6;
    ref_surfacic(3)(3) = 10; 
    ref_surfacic(4)(1) = 4; ref_surfacic(4)(2) = 6; ref_surfacic(4)(3) = 7;
    ref_surfacic(5)(2) = 8; 
    ref_surfacic(6)(1) = 4;
    ref_surfacic(7)(1) = 9; ref_surfacic(7)(2) = 9; ref_surfacic(7)(3) = 9; 
    
    ref_volume(0)(1) = 1; ref_volume(0)(2) = 1; ref_volume(0)(3) = 1;
    ref_volume(1)(1) = 1; ref_volume(1)(2) = 0; ref_volume(1)(3) = 0;
    ref_volume(2)(1) = 1; ref_volume(2)(2) = 0; ref_volume(2)(3) = 1;
    ref_volume(3)(1) = 1; ref_volume(3)(2) = 1; ref_volume(3)(3) = 1;
    
    mesh.ExtrudeSurfaceMesh(mesh2d, nb_interv, step_layer, ref_surfacic,
                            ref_volume, vec_e1, vec_e2, vec_e3);
    
    step.Reallocate(3); step.Fill(); Mlt(1.0/2, step);
    
    mesh.SubdivideMesh(step);    
    mesh.Write("result26.mesh");
        
    int nb_vert = 0, nb_elt = 0;
    Vector<bool> VertexOnSubdomain(mesh.GetNbVertices()), ElementOnSubdomain(mesh.GetNbElt());
    VertexOnSubdomain.Fill(false); ElementOnSubdomain.Fill(false);
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        int nv = mesh.Element(i).GetNbVertices();
        R3 center;
        for (int j = 0; j < nv; j++)
          center += mesh.Vertex(mesh.Element(i).numVertex(j));
        
        Mlt(1.0/nv, center);
        if (center(0) < 0)
          {
            ElementOnSubdomain(i) = true;
            for (int j = 0; j < nv; j++)
              if (!VertexOnSubdomain(mesh.Element(i).numVertex(j)))
                {
                  VertexOnSubdomain(mesh.Element(i).numVertex(j)) = true;
                  nb_vert++;
                }
            
            nb_elt++;
          }        
      }
    
    mesh.CreateSubmesh(add_mesh, nb_vert, nb_elt, VertexOnSubdomain, ElementOnSubdomain);
    
    add_mesh.SubdivideMesh(step);    
    add_mesh.Write("result27.mesh");
  }

  return FinalizeMontjoie();  
}
