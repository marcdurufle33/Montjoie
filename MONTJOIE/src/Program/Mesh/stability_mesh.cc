#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

Real_wp coef_g = 10;
Real_wp coef_alpha = 1;
Real_wp coef_beta = 1;

void GenerateListPoint(const Mesh<Dimension2>& mesh, int num_init, int inext, IVect& Point_curve)
{
  IVect ref_cond(mesh.GetNbReferences()); ref_cond.Fill();
  int ref = 4;
  int ipoint = inext, iback = num_init, iedge = -1, nb_points = 1;
  Point_curve.Reallocate(mesh.GetNbVertices());
  Point_curve(0) = num_init;
  while ((nb_points == 0)||((ipoint != num_init)&&(ipoint != -1)))
    {
      Point_curve(nb_points) = ipoint;
      inext = mesh.FindFollowingVertex(iback, ipoint, ref, ref_cond, iedge);
      iback = ipoint; ipoint = inext;
      nb_points++;
    }
  
  Point_curve.Resize(nb_points);
}

void GenerateDamier(const Mesh<Dimension2>& mesh0, const Mesh<Dimension2>& mesh1,
		    int n0, int n0_next, int n1, int n1_next,
		    const Real_wp& x0, const Real_wp& x1, Mesh<Dimension3>& mesh_add)
{
  IVect num_curve0, num_curve1;
  GenerateListPoint(mesh0, n0, n0_next, num_curve0);

  /* R2 pt0 = mesh0.Vertex(num_curve0(0));
  Real_wp dist_min = pt0.Distance(mesh1.Vertex(0));
  for (int k = 1; k < mesh1.GetNbVertices(); k++)
    if (pt0.Distance(mesh1.Vertex(k)) < dist_min)
      {
	n1 = k;
	dist_min = pt0.Distance(mesh1.Vertex(k));
      }
  */
  
  GenerateListPoint(mesh1, n1, n1_next, num_curve1);

  int N0 = mesh0.GetNbVertices(), N1 = mesh1.GetNbVertices();
  mesh_add.Clear();
  mesh_add.ReallocateVertices(N0+N1);
  mesh_add.ReallocateBoundariesRef(2*(N0+N1));
  for (int i = 0; i < N0; i++)
    mesh_add.Vertex(i).Init(coef_g*x0, coef_alpha*mesh0.Vertex(i)(0), coef_beta*mesh0.Vertex(i)(1));

  for (int i = 0; i < N1; i++)
    mesh_add.Vertex(N0+i).Init(coef_g*x1, coef_alpha*mesh1.Vertex(i)(0), coef_beta*mesh1.Vertex(i)(1));

  int ne = 0;
  int j = 0; int iback = num_curve0(0);
  for (int i = 1; i <= N0; i++)
    {
      int i2 = num_curve0(i%N0);
      int j2 = num_curve1(j%N1);
      Real_wp L0 = Real_wp(i) / N0;
      Real_wp L1 = Real_wp(j) / N1;

      mesh_add.BoundaryRef(ne).InitTriangular(iback, i2, N0+j2, 1);
      ne++;
      
      while ((L1 < L0) && (j < N1))
	{
	  int jnext = num_curve1((j+1)%N1);
	  mesh_add.BoundaryRef(ne).InitTriangular(i2, N0+j2, N0+jnext, 1);
	  ne++; j++;
	  L1 = Real_wp(j) / N1;
	}

      iback = i2;
    }

  mesh_add.ResizeBoundariesRef(ne);

  /* Mesh<Dimension2> mesh2(mesh0);
  mesh2.AppendMesh(const_cast<Mesh<Dimension2>& >(mesh1));
  mesh2.Write("contour.mesh"); */
}


void AppendMesh(const Mesh<Dimension2>& mesh0, const Mesh<Dimension2>& mesh1,
		const Mesh<Dimension3>& mesh_add, Mesh<Dimension3>& mesh)
{
  int nb_old_vert = mesh.GetNbVertices();
  int nb_vert0 = mesh0.GetNbVertices();
  mesh.ResizeVertices(mesh1.GetNbVertices() + nb_old_vert);
  int offset = nb_old_vert - nb_vert0;
  for (int i = 0; i < mesh1.GetNbVertices(); i++)
    mesh.Vertex(nb_old_vert + i) = mesh_add.Vertex(nb_vert0 + i);

  int nb_old_tri = mesh.GetNbBoundaryRef();
  mesh.ResizeBoundariesRef(nb_old_tri + mesh_add.GetNbBoundaryRef());
  for (int i = 0; i < mesh_add.GetNbBoundaryRef(); i++)
    {
      int n0 = mesh_add.BoundaryRef(i).numVertex(0);
      int n1 = mesh_add.BoundaryRef(i).numVertex(1);
      int n2 = mesh_add.BoundaryRef(i).numVertex(2);
      mesh.BoundaryRef(nb_old_tri + i).InitTriangular(offset + n0, offset+n1, offset+n2, 1);
    }
}

//./test.x 7 3 0.1361 0.13615 1.0041 1.0041 -12.5532 -12.5488 501 1 501 1
// ./test.x 7 3 0.18 0.26 0 4 -35 -15 201 201 201 1
// ./test.x 9 3 0.14 0.24 0 1 -40 -20 251 251 501 1

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 2)
    {
      cout<<"Cette commande demande au moins un argument"<<endl;
      cout << "Usage : ./test.x append" << endl;
      return -1;
    }

  string file_name(argv[1]);
  if (file_name == "append")
    {
      // cas ou concatene plusieurs maillages 2-D pour faire un maillage 3-D
      if (argc < 6)
	{
	  cout<<"Cette commande demande cinq arguments"<<endl;
	  cout << "Usage : ./test.x append n0 n1 root list_pt" << endl;
	  return -1;
	}

      Mesh<Dimension2> mesh0, mesh1;
      Mesh<Dimension3> mesh, mesh_add;
      Real_wp gamma0, gamma1;

      int n0 = atoi(argv[2]);
      int n1 = atoi(argv[3]);
      string base_name(argv[4]);
      string pt_name(argv[5]);
      Real_wp gamma_b(-1), gamma_f(-1);
      if (argc == 8)
	{
	  gamma_b = atof(argv[6]);
	  gamma_f = atof(argv[7]);
	}
      
      IVect num_pt;
      num_pt.ReadText(pt_name);
      
      mesh0.Read(base_name + to_str(n0) + ".mesh");
      for (int n = n0; n < n1; n++)
	{
	  int k = n-n0;
	  mesh1.Read(base_name + to_str(n+1) + ".mesh");
	  gamma0 = Real_wp(n) / 1000;
	  gamma1 = Real_wp(n+1) / 1000;
	  if ((n == n0) && (gamma_b > 0))
	    gamma0 = gamma_b;

	  if ((n == n1-1) && (gamma_f > 0))
	    gamma1 = gamma_f;

	  //DISP(gamma0); DISP(gamma1);
	  GenerateDamier(mesh0, mesh1, num_pt(2*k)-1, num_pt(2*k+1)-1,
			 num_pt(2*k+2)-1, num_pt(2*k+3)-1, gamma0, gamma1, mesh_add);
	  
	  // mesh_add.Write("add.mesh");
	  
	  if (n == n0)
	    mesh = mesh_add;
	  else
	    {
	      //mesh.AppendMesh(mesh_add);
	      AppendMesh(mesh0, mesh1, mesh_add, mesh);
	    }
	  
	  mesh0 = mesh1;
	}

      mesh.Write("final.mesh");
      
      return FinalizeMontjoie();
      
    }

  if (argc != 11)
    {
      cout<<"Cette commande demande au quatre arguments"<<endl;
      cout << "Usage : ./test.x file.dat nx ny nz xmin xmax ymin ymax zmin zmax" << endl;
      return -1;
    }

  int Nx = atoi(argv[2]);
  int Ny = atoi(argv[3]);
  int Nz = atoi(argv[4]);
  Real_wp xmin = atof(argv[5]);
  Real_wp xmax = atof(argv[6]);
  Real_wp ymin = atof(argv[7]);
  Real_wp ymax = atof(argv[8]);
  Real_wp zmin = atof(argv[9]);
  Real_wp zmax = atof(argv[10]);

  VectReal_wp cfl;
  cfl.ReadText(file_name);

  if (Nz > 1)
    {
      // constructing a 3-D mesh
      Array3D<int> IndexVertex(Nx, Ny, Nz);
      IndexVertex.Fill(-1);
      
      int nb_vert = 0;
      for (int i = 0; i < Nx; i++)
	for (int j = 0; j < Ny; j++)
	  for (int k = 0; k < Nz; k++)
	    {
	      int num = Nz*(i*Ny + j) + k;
	      if (cfl(num) > 1e100)
		IndexVertex(i, j, k) = nb_vert++;
	    }
      
      VectReal_wp x_div, y_div, z_div;
      Linspace(xmin, xmax, Nx, x_div);
      Linspace(ymin, ymax, Ny, y_div);
      Linspace(zmin, zmax, Nz, z_div);
      
      Mesh<Dimension3> mesh;
      mesh.ReallocateVertices(nb_vert);
      for (int i = 0; i < Nx; i++)
	for (int j = 0; j < Ny; j++)
	  for (int k = 0; k < Nz; k++)
	    if (IndexVertex(i, j, k) >= 0)
	      {
		int num = IndexVertex(i, j, k);
		mesh.Vertex(num).Init(x_div(i), y_div(j), z_div(k));
	      }
      
      DISP(nb_vert);
      
      Array3D<int> IndexHexa(Nx, Ny, Nz);
      Array3D<int> IndexTetra(Nx, Ny, Nz);
      IndexHexa.Fill(-1); IndexTetra.Fill(-1);
      int nb_hexa = 0, nb_tetra = 0;
      for (int i = 0; i < Nx-1; i++)
	for (int j = 0; j < Ny-1; j++)
	  for (int k = 0; k < Nz-1; k++)
	    {
	      nb_vert = 0;
	      if (IndexVertex(i, j, k) >= 0)
		nb_vert++;
	  
	      if (IndexVertex(i+1, j, k) >= 0)
		nb_vert++;
	      
	      if (IndexVertex(i, j+1, k) >= 0)
		nb_vert++;
	      
	      if (IndexVertex(i, j, k+1) >= 0)
		nb_vert++;
	      
	      if (IndexVertex(i, j+1, k+1) >= 0)
		nb_vert++;
	      
	      if (IndexVertex(i+1, j, k+1) >= 0)
		nb_vert++;
	      
	      if (IndexVertex(i+1, j+1, k) >= 0)
		nb_vert++;
	      
	      if (IndexVertex(i+1, j+1, k+1) >= 0)
		nb_vert++;
	      
	      if (nb_vert == 8)
		IndexHexa(i, j, k) = nb_hexa++;
	      //else if (nb_vert == 4)
	      //IndexTetra(i, j, k) = nb_tetra++;
	    }
      
      DISP(nb_tetra);
      
      mesh.ReallocateElements(nb_hexa + nb_tetra);
      IVect num_pt(8), nump(8), permut(8), index_p(8);
      TinyVector<int, 4> num_face, num_face_opp;
      for (int i = 0; i < Nx-1; i++)
	for (int j = 0; j < Ny-1; j++)
	  for (int k = 0; k < Nz-1; k++)
	    {
	      num_pt(0) = IndexVertex(i, j, k);
	      num_pt(1) = IndexVertex(i+1, j, k);
	      num_pt(2) = IndexVertex(i+1, j+1, k);
	      num_pt(3) = IndexVertex(i, j+1, k);
	      num_pt(4) = IndexVertex(i, j, k+1);
	      num_pt(5) = IndexVertex(i+1, j, k+1);
	      num_pt(6) = IndexVertex(i+1, j+1, k+1);
	      num_pt(7) = IndexVertex(i, j+1, k+1);
	      
	      nb_vert = 0;
	      index_p.Zero();
	      for (int m = 0; m < 8; m++)
		if (num_pt(m) >= 0)
		  {
		    index_p(m) = 1;
		    nb_vert++;
		  }
	      
	      permut.Fill();
	      nump = num_pt;
	      Sort(nump, permut);
	      
	      if (IndexTetra(i, j, k) >= 0)
		{
		  int num = IndexTetra(i, j, k);	      
		  mesh.Element(num).InitTetrahedral(nump(4), nump(5), nump(6), nump(7), 1);
		}
	      else if (IndexHexa(i, j, k) >= 0)
		{
		  int num = nb_tetra + IndexHexa(i, j, k);
		  if (nb_vert == 8)
		    mesh.Element(num).InitHexahedral(num_pt, 1);
		  else
		    {
		      //FindOppositeFaces(index_p, num_pt, num_face, num_face_opp);
		      //mesh.Element(num).InitHexahedral(num_face(0), num_face(1), num_face(2), num_face(3),
		      //				       num_face_opp(0), num_face_opp(1),
		      //				       num_face_opp(2), num_face_opp(3), 1);
		      
		      //int test_input; cout << "waiting" << endl; cin >> test_input;
		    }
		}
	    }
  
      DISP(cfl.GetNormInf()); DISP(nb_hexa);
      
      mesh.Write("test.mesh");
    }
  else
    {
      // constructing a 2-D mesh

      Matrix<int> IndexVertex(Nx, Ny);
      IndexVertex.Fill(-1);
      int nb_vert = 0;
      for (int i = 0; i < Nx; i++)
	for (int j = 0; j < Ny; j++)
	  {
	    int num = i*Ny + j;
	    if (cfl(num) > 1e100)
	      IndexVertex(i, j) = nb_vert++;
	  }
      
      VectReal_wp x_div, y_div;
      Linspace(xmin, xmax, Nx, x_div);
      Linspace(ymin, ymax, Ny, y_div);
      
      Mesh<Dimension2> mesh;
      mesh.ReallocateVertices(nb_vert);
      for (int i = 0; i < Nx; i++)
	for (int j = 0; j < Ny; j++)
	  if (IndexVertex(i, j) >= 0)
	    {
	      int num = IndexVertex(i, j);
	      mesh.Vertex(num).Init(x_div(i), y_div(j));
	    }
      
      DISP(nb_vert);

      Matrix<int> IndexQuad(Nx, Ny);
      Matrix<int> IndexTri(Nx, Ny);
      IndexQuad.Fill(-1); IndexTri.Fill(-1);
      int nb_quad = 0, nb_tri = 0;
      for (int i = 0; i < Nx-1; i++)
	for (int j = 0; j < Ny-1; j++)
	  {
	    nb_vert = 0;
	    if (IndexVertex(i, j) >= 0)
	      nb_vert++;
	    
	    if (IndexVertex(i+1, j) >= 0)
	      nb_vert++;
	    
	    if (IndexVertex(i, j+1) >= 0)
	      nb_vert++;
	    
	    if (IndexVertex(i+1, j+1) >= 0)
	      nb_vert++;
	    
	    if (nb_vert == 4)
		IndexQuad(i, j) = nb_quad++;
	    else if (nb_vert == 3)
	      IndexTri(i, j) = nb_tri++;
	  }
      
      DISP(nb_tri); DISP(nb_quad);
      
      mesh.ReallocateElements(nb_quad + nb_tri);
      IVect num_pt(4), nump(4);
      for (int i = 0; i < Nx-1; i++)
	for (int j = 0; j < Ny-1; j++)
	  {
	    num_pt(0) = IndexVertex(i, j);
	    num_pt(1) = IndexVertex(i+1, j);
	    num_pt(2) = IndexVertex(i+1, j+1);
	    num_pt(3) = IndexVertex(i, j+1);

	    if (IndexTri(i, j) >= 0)
	      {
		nump = num_pt;
		Sort(nump);
		int num = IndexTri(i, j);
		mesh.Element(num).InitTriangular(nump(1), nump(2), nump(3), 1);
	      }
	    else if (IndexQuad(i, j) >= 0)
	      {
		int num = nb_tri + IndexQuad(i, j);
		mesh.Element(num).InitQuadrangular(num_pt(0), num_pt(1), num_pt(2), num_pt(3), 1);
	      }	
	  }

      mesh.Write("test.mesh");
    }
  
  return FinalizeMontjoie();
}
