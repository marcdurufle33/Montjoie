#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

#include "MeshVerif.cxx"

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  Mesh<Dimension3> mesh;
  int order_geom = 6;
  mesh.SetGeometryOrder(order_geom);
  
  // testing ProjectToCurve
  mesh.ResizeNbReferences(10);
  VectReal_wp param(4);
  param(0) = Real_wp(1); param(1) = Real_wp(2); param(2) = Real_wp(1.5); param(3) = Real_wp(1.8);
  
  VectReal_wp x, y, z; int Nx = 90, Ny = 80, Nz = 70;

  double xmin = -4.0, xmax = 4.0, ymin = -5.0, ymax = 5.0, zmin = -4.5, zmax = 4.5;
  GridInterpolationFull<Dimension3> grid;  
  grid.SetVolumeOutput(R3(xmin, ymin, zmin), R3(xmax, ymax, zmax), Nx, Ny, Nz);

  Linspace(xmin, xmax, Nx, x);
  Linspace(ymin, ymax, Ny, y);
  Linspace(zmin, zmax, Nz, z);

  int N = 20;
  VectReal_wp teta;
  Linspace(Real_wp(0), pi_wp, N+1, teta);

  // generating a spherical mesh
  // testing GenerateSurfaceOfRevolution
  Mesh<Dimension2> mesh_curve;
  mesh_curve.ResizeVertices(N+1);
  for (int i = 0; i <= N; i++)
    mesh_curve.Vertex(i).Init(param(3)*sin(teta(i)), param(3)*cos(teta(i)));
  
  mesh_curve.ResizeBoundariesRef(N);
  for (int i = 0; i < N; i++)
    {
      int n0 = i, n1 = i+1;
      int ref = 1;
      if (i >= N/2)
	ref = 2;
      
      mesh_curve.BoundaryRef(i).Init(n0, n1, ref);
    }

  {
    IVect ref_cond(10);
    ref_cond.Zero(); ref_cond(1) = 1; ref_cond(2) = 1;

    IVect IndexFace, IndexVertex; VectReal_wp AngleVertex;
    mesh_curve.GenerateSurfaceOfRevolution(1, mesh, IndexFace, IndexVertex, AngleVertex, ref_cond);
  }

  mesh.Write("sphere.mesh");  
  
  Mesh<Dimension3> mesh2;
  mesh2.Read("src/Program/Unit/Mesh/sphere_ref.mesh");

  if (!IsMeshEqual(mesh, mesh2))
    {
      cout << "GenerateSurfaceOfRevolution incorrect" << endl;
      abort();
    }

  TranslateMesh(mesh, R3(param(0), param(1), param(2)));
  mesh.SetCurveType(1, mesh.CURVE_SPHERE);
  mesh.SetCurveParameter(1, param);

  VectReal_wp u(Nx*Ny*Nz);
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      for (int k = 0; k < Nz; k++)
	{
	  R3 pt(x(i), y(j), z(k));
	  R3 proj(pt);
	  mesh.ProjectToCurve(proj, 1);
	  u(k*Nx*Ny + j*Nx + i) = pt.Distance(proj);
	  
	  Real_wp dist = pt.Distance(proj);
	  for (int p = 0; p < mesh.GetNbVertices(); p++)
	    if (pt.Distance(mesh.Vertex(p)) < dist-1e-12)
	      {
		DISP(pt); DISP(proj); DISP(dist); DISP(mesh.Vertex(k));
		DISP(pt.Distance(mesh.Vertex(k)));
		cout << "ProjectToCurve (sphere) incorrect" << endl;
		abort();
	      }	  
	}
  
  WriteMatlab(u, grid, "dist_sphere.dat", true, false);

  Real_wp h2 = 2.0, h1 = -1.9; Real_wp radius = 2.4;
  Real_wp dx = radius / (N/3);
  for (int i = 0; i <= N/3; i++)
    mesh_curve.Vertex(i).Init(i*dx, h2);

  int i1 = N/3, i2 = 2*N/3;
  Real_wp dy = (h1-h2)/(i2-i1);
  for (int i = i1+1; i <= i2; i++)
    mesh_curve.Vertex(i).Init(radius, h2+(i-i1)*dy);
  
  dx = -radius / (N-i2);
  for (int i = i2+1; i <= N; i++)
    mesh_curve.Vertex(i).Init(radius+(i-i2)*dx, h1);
  
  for (int i = 0; i < N; i++)
    {
      DISP(mesh_curve.Vertex(i));
      int n0 = i, n1 = i+1;
      int ref = 4;
      if ((i >= N/3) && (i < 2*N/3))
	ref = 2;
      
      mesh_curve.BoundaryRef(i).Init(n0, n1, ref);
    }

  {
    IVect ref_cond(10);
    ref_cond.Fill(1);

    IVect IndexFace, IndexVertex; VectReal_wp AngleVertex;
    mesh_curve.GenerateSurfaceOfRevolution(1, mesh2, IndexFace, IndexVertex, AngleVertex, ref_cond);
    mesh2.Write("cylinder.mesh");
  }

  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    mesh.BoundaryRef(i).SetReference(1);

  TranslateMesh(mesh2, R3(param(0), param(1), param(2)));  
  mesh.AppendMesh(mesh2);
  mesh.Write("sphere_cyl.mesh");  

  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
