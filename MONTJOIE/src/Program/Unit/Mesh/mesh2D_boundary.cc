#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

#include "MeshVerif.cxx"

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  Mesh<Dimension2> mesh;
  int order_geom = 6;
  mesh.SetGeometryOrder(order_geom);
  
  // testing ProjectToCurve
  mesh.ResizeNbReferences(10);
  mesh.SetCurveType(1, mesh.CURVE_CIRCLE);
  VectReal_wp param(3);
  param(0) = Real_wp(1); param(1) = Real_wp(2); param(2) = Real_wp(1.5);
  mesh.SetCurveParameter(1, param);

  VectReal_wp x, y; int Nx = 200, Ny = 300;

  double xmin = -4.0, xmax = 4.0, ymin = -5.0, ymax = 5.0;
  GridInterpolationFull<Dimension2> grid;  
  grid.SetPlaneOutput(xmin, xmax, ymin, ymax, Nx, Ny);

  Linspace(xmin, xmax, Nx, x);
  Linspace(ymin, ymax, Ny, y);

  int N = 100;
  VectReal_wp teta;
  Linspace(Real_wp(0), 2*pi_wp, N+1, teta); teta.Resize(N);
  
  VectR2 ptCurve(N);
  for (int i = 0; i < N; i++)
    ptCurve(i).Init(param(0) + param(2)*cos(teta(i)), param(1) + param(2)*sin(teta(i)));

  mesh.ResizeVertices(N);
  for (int i = 0; i < N; i++)
    mesh.Vertex(i) = ptCurve(i);

  mesh.ResizeBoundariesRef(N);
  for (int i = 0; i < N; i++)
    {
      int n0 = i, n1 = i+1;
      if (i == N-1)
	{
	  n0 = 0;
	  n1 = N-1;
	}
      
      mesh.BoundaryRef(i).Init(n0, n1, 1);
    }

  int offset = N;
  
  VectReal_wp u(Nx*Ny);
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      {
	R2 pt(x(i), y(j));
	R2 proj(pt);
	mesh.ProjectToCurve(proj, 1);
	u(j*Nx+i) = pt.Distance(proj);

	Real_wp dist = pt.Distance(proj);
	for (int k = 0; k < N; k++)
	  if (pt.Distance(ptCurve(k)) < dist-1e-12)
	    {
	      DISP(pt); DISP(proj); DISP(dist); DISP(ptCurve(k));
	      DISP(pt.Distance(ptCurve(k)));
	      cout << "ProjectToCurve (circle) incorrect" << endl;
	      abort();
	    }	  
      }
  
  WriteMatlab(u, grid, "dist_circle.dat", true, false);

  mesh.SetCurveType(2, mesh.CURVE_ELLIPSE);
  param.Reallocate(4);
  param(0) = Real_wp(1); param(1) = Real_wp(2);
  param(2) = Real_wp(2.5); param(3) = Real_wp(1.2);
  mesh.SetCurveParameter(2, param);

  for (int i = 0; i < N; i++)
    ptCurve(i).Init(param(0) + param(2)*cos(teta(i)), param(1) + param(3)*sin(teta(i)));

  mesh.ResizeVertices(offset+N);
  for (int i = 0; i < N; i++)
    mesh.Vertex(offset+i) = ptCurve(i);

  mesh.ResizeBoundariesRef(offset+N);
  for (int i = 0; i < N; i++)
    {
      int n0 = offset+i, n1 = offset+i+1;
      if (i == N-1)
	{
	  n0 = offset;
	  n1 = offset+N-1;
	}
      
      mesh.BoundaryRef(offset+i).Init(n0, n1, 2);
    }

  offset += N;
  
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      {
	R2 pt(x(i), y(j));
	R2 proj(pt);
	mesh.ProjectToCurve(proj, 2);
	u(j*Nx+i) = pt.Distance(proj);

	Real_wp dist = pt.Distance(proj);
	for (int k = 0; k < N; k++)
	  if (pt.Distance(ptCurve(k)) < dist-1e-4)
	    {
	      DISP(pt); DISP(proj); DISP(dist); DISP(ptCurve(k));
	      DISP(pt.Distance(ptCurve(k)));
	      cout << "ProjectToCurve (ellipse) incorrect" << endl;
	      abort();
	    }
      }
  
  WriteMatlab(u, grid, "dist_ellipse.dat", true, false);

  mesh.SetCurveType(3, mesh.CURVE_PEANUT);
  param.Reallocate(5);
  param(0) = Real_wp(1); param(1) = Real_wp(2);
  param(2) = Real_wp(1.2); param(3) = Real_wp(2.5); param(4) = Real_wp(4.0);
  mesh.SetCurveParameter(3, param);

  for (int i = 0; i < N; i++)
    ptCurve(i).Init(param(0) + param(2)*cos(teta(i))*sqrt(1.0+param(4)*square(sin(teta(i)))),
		    param(1) + param(3)*sin(teta(i)));
  
  mesh.ResizeVertices(offset+N);
  for (int i = 0; i < N; i++)
    mesh.Vertex(offset+i) = ptCurve(i);

  mesh.ResizeBoundariesRef(offset+N);
  for (int i = 0; i < N; i++)
    {
      int n0 = offset+i, n1 = offset+i+1;
      if (i == N-1)
	{
	  n0 = offset;
	  n1 = offset+N-1;
	}
      
      mesh.BoundaryRef(offset+i).Init(n0, n1, 3);
    }

  offset += N;
  
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      {
	//R2 pt(0.9, 2.4);
	R2 pt(x(i), y(j));
	R2 proj(pt);
	mesh.ProjectToCurve(proj, 3);
	u(j*Nx+i) = pt.Distance(proj);

	Real_wp dist = pt.Distance(proj);
	if (dist < Real_wp(0.5))
	  for (int k = 0; k < N; k++)
	    if (pt.Distance(ptCurve(k)) < dist-1e-4)
	      {
		DISP(pt); DISP(proj); DISP(dist); DISP(ptCurve(k));
		DISP(pt.Distance(ptCurve(k)));
		cout << "ProjectToCurve (peanut) incorrect" << endl;
		abort();
	      }
      }
  
  WriteMatlab(u, grid, "dist_peanut.dat", true, false);

  mesh.FindConnectivity();
  
  // testing ComputeParameterValueOnCurves
  mesh.ComputeParameterValueOnCurves();

  const VectR2& val_teta = mesh.GetValueParameterOnBoundary();
  Real_wp xc = mesh.GetCurveParameter(1)(0);
  Real_wp yc = mesh.GetCurveParameter(1)(1);
  Real_wp rc = mesh.GetCurveParameter(1)(2);
  
  Real_wp xe = mesh.GetCurveParameter(2)(0);
  Real_wp ye = mesh.GetCurveParameter(2)(1);
  Real_wp ae = mesh.GetCurveParameter(2)(2);
  Real_wp be = mesh.GetCurveParameter(2)(3);

  Real_wp xp = mesh.GetCurveParameter(3)(0);
  Real_wp yp = mesh.GetCurveParameter(3)(1);
  Real_wp ap = mesh.GetCurveParameter(3)(2);
  Real_wp bp = mesh.GetCurveParameter(3)(3);
  Real_wp cp = mesh.GetCurveParameter(3)(4);

  for (int i = 0; i < N; i++)
    {
      // circle
      Real_wp teta0 = val_teta(i)(0);
      Real_wp teta1 = val_teta(i)(1);
      int n0 = mesh.BoundaryRef(i).numVertex(0);
      int n1 = mesh.BoundaryRef(i).numVertex(1);
      R2 pt0 = mesh.Vertex(n0), pt1 = mesh.Vertex(n1);
      R2 pt0_b(xc + rc*cos(teta0), yc + rc*sin(teta0));
      R2 pt1_b(xc + rc*cos(teta1), yc + rc*sin(teta1));
      if ((pt0 != pt0_b) || (pt1 != pt1_b))
	{
	  cout << "ComputeParameterValueOnCurves incorrect" << endl;
	  abort();
	}
      
      // ellipse
      teta0 = val_teta(N+i)(0);
      teta1 = val_teta(N+i)(1);
      n0 = mesh.BoundaryRef(N+i).numVertex(0);
      n1 = mesh.BoundaryRef(N+i).numVertex(1);
      pt0 = mesh.Vertex(n0); pt1 = mesh.Vertex(n1);
      pt0_b.Init(xe + ae*cos(teta0), ye + be*sin(teta0));
      pt1_b.Init(xe + ae*cos(teta1), ye + be*sin(teta1));
      if ((pt0 != pt0_b) || (pt1 != pt1_b))
	{
	  cout << "ComputeParameterValueOnCurves incorrect" << endl;
	  abort();
	}

      // peanut
      teta0 = val_teta(2*N+i)(0);
      teta1 = val_teta(2*N+i)(1);
      n0 = mesh.BoundaryRef(2*N+i).numVertex(0);
      n1 = mesh.BoundaryRef(2*N+i).numVertex(1);
      pt0 = mesh.Vertex(n0); pt1 = mesh.Vertex(n1);
      pt0_b.Init(xp + ap*cos(teta0)*sqrt(1.0+cp*sin(teta0)*sin(teta0)), yp + bp*sin(teta0));
      pt1_b.Init(xp + ap*cos(teta1)*sqrt(1.0+cp*sin(teta1)*sin(teta1)), yp + bp*sin(teta1));      
      if ((pt0 != pt0_b) || (pt1 != pt1_b))
	{
	  DISP(i); DISP(pt0); DISP(pt0_b); DISP(pt1); DISP(pt1_b);
	  cout << "ComputeParameterValueOnCurves incorrect" << endl;
	  abort();
	}

      // spline et local_spline a faire
    }

  Globatto<Real_wp> lob;
  lob.ConstructQuadrature(10, lob.QUADRATURE_LOBATTO);
  
  // testing GetPointsOnCurve
  VectR2 points;
  for (int i = 0; i < N; i++)
    {
      // circle
      mesh.GetPointsOnCurve(i, 1, lob.Points(), points);
      Real_wp teta0 = val_teta(i)(0);
      Real_wp teta1 = val_teta(i)(1);
      for (int k = 0; k <= lob.GetOrder(); k++)
	{
	  Real_wp teta = (1.0-lob.Points(k))*teta0 + lob.Points(k)*teta1;
	  R2 pt(xc + rc*cos(teta), yc + rc*sin(teta));
	  if (pt != points(k))
	    {
	      cout << "GetPointsOnCurve incorrect" << endl;
	      abort();
	    }
	}

      // ellipse
      mesh.GetPointsOnCurve(N+i, 2, lob.Points(), points);
      teta0 = val_teta(N+i)(0);
      teta1 = val_teta(N+i)(1);
      int n0 = mesh.BoundaryRef(N+i).numVertex(0);
      int n1 = mesh.BoundaryRef(N+i).numVertex(1);
      R2 pt0 = mesh.Vertex(n0); R2 pt1 = mesh.Vertex(n1);
      if ((points(0) != pt0) || (points(lob.GetOrder()) != pt1))
	{
	  cout << "GetPointsOnCurve (ellipse0) incorrect" << endl;
	  abort();
	}

      VectReal_wp length(lob.GetOrder()+1);
      length.Zero();
      for (int k = 1; k <= lob.GetOrder(); k++)
	{
	  Real_wp teta_k = acos((points(k)(0) - xe)/ae);
	  if (points(k)(1) - ye < 0)
	    teta_k = -teta_k;

	  if (teta_k > teta1 + pi_wp)
	    teta_k -= 2.0*pi_wp;

	  if (teta_k < teta0 - pi_wp)
	    teta_k += 2.0*pi_wp;
	  
	  for (int q = 0; q <= lob.GetOrder(); q++)
	    {
	      Real_wp teta = (1.0-lob.Points(q))*teta0 + lob.Points(q)*teta_k;
	      length(k) += lob.Weights(q)*sqrt(square(ae*sin(teta)) + square(be*cos(teta)));
	    }
	  
	  length(k) *= (teta_k - teta0);
	}

      Mlt(1.0/length(lob.GetOrder()), length);
      for (int k = 1; k <= lob.GetOrder(); k++)
	if (abs(lob.Points(k) - length(k)) > 1e-6)
	  {
	    DISP(i); DISP(teta0); DISP(teta1); DISP(length(k)); DISP(lob.Points(k));
	    cout << "GetPointsOnCurve (ellipse) incorrect" << endl;
	    abort();
	  }

      // peanut
      mesh.GetPointsOnCurve(2*N+i, 3, lob.Points(), points);
      teta0 = val_teta(2*N+i)(0);
      teta1 = val_teta(2*N+i)(1);
      n0 = mesh.BoundaryRef(2*N+i).numVertex(0);
      n1 = mesh.BoundaryRef(2*N+i).numVertex(1);
      pt0 = mesh.Vertex(n0); pt1 = mesh.Vertex(n1);
      if ((points(0) != pt0) || (points(lob.GetOrder()) != pt1))
	{
	  cout << "GetPointsOnCurve (peanut0) incorrect" << endl;
	  abort();
	}

      length.Zero();
      for (int k = 1; k <= lob.GetOrder(); k++)
	{
	  Real_wp teta_k = asin((points(k)(1) - yp)/bp);
	  if (points(k)(0) - xp < 0)
	    teta_k = pi_wp-teta_k;

	  if (teta_k > teta1 + pi_wp)
	    teta_k -= 2.0*pi_wp;

	  if (teta_k < teta0 - pi_wp)
	    teta_k += 2.0*pi_wp;
	  
	  for (int q = 0; q <= lob.GetOrder(); q++)
	    {
	      Real_wp teta = (1.0-lob.Points(q))*teta0 + lob.Points(q)*teta_k;
	      Real_wp radical = sqrt(1.0+cp*square(sin(teta)));
	      Real_wp coef = -radical + cp*square(cos(teta))/radical;
	      length(k) += lob.Weights(q)*sqrt(square(ap*sin(teta)*coef) + square(bp*cos(teta)));
	    }
	  
	  length(k) *= (teta_k - teta0);
	}

      Mlt(1.0/length(lob.GetOrder()), length);
      for (int k = 1; k <= lob.GetOrder(); k++)
	if (abs(lob.Points(k) - length(k)) > 1e-6)
	  {
	    DISP(length); DISP(lob.Points());
	    DISP(length(k)); DISP(lob.Points(k));
	    cout << "GetPointsOnCurve (peanut) incorrect" << endl;
	    abort();
	  }

    }

  // testing ProjectPointsOnCurves()
  lob.ConstructQuadrature(mesh.GetGeometryOrder(), lob.QUADRATURE_LOBATTO);
  mesh.ProjectPointsOnCurves();

  if (!IsMeshCorrect(mesh, order_geom))
    {
      cout << "ProjectPointsOnCurves incorrect" << endl;
      abort();
    }

  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
