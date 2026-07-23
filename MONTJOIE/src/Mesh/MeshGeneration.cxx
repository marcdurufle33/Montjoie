#ifndef MONTJOIE_FILE_MESH_GENERATION_CXX

#include "MeshGeneration.hxx"

namespace Montjoie
{
  void CreatePyramidalLayer(const Mesh<Dimension3>& mesh_surf, Mesh<Dimension3>& mesh_tri,
			    Mesh<Dimension3>& mesh_pyramid,
                            const Real_wp& height_pyramid, int new_ref)
  {
    // on decoupe le maillage de peau en triangles
    int nb_tri = mesh_surf.GetNbTrianglesRef();
    int nb_quad = mesh_surf.GetNbQuadranglesRef();
    mesh_tri = mesh_surf;
    Mesh<Dimension3> mesh_tet;
    
    mesh_tri.SplitIntoTetrahedra();
    
    // on mets des tetras dedans pour savoir ou est l'interieur
    CreateTetrahedralMesh(mesh_tri, mesh_tet);
    
    // puis sur chaque quadrangle on essaie de rajouter une pyramide
    int nb_vert_surf = mesh_surf.GetNbVertices();
    int nb_vertices = nb_vert_surf + nb_quad;
    int nb_faces = nb_tri + 5*nb_quad;
    mesh_tri.ReallocateVertices(nb_vertices);
    mesh_tri.ReallocateBoundariesRef(nb_tri + 4*nb_quad);

    mesh_pyramid.ReallocateVertices(nb_vertices);
    mesh_pyramid.ReallocateBoundariesRef(nb_faces);
    mesh_pyramid.ReallocateElements(nb_quad);
    for (int i = 0; i < nb_vert_surf; i++)
      {
	mesh_tri.Vertex(i) = mesh_surf.Vertex(i);
	mesh_pyramid.Vertex(i) = mesh_surf.Vertex(i);
      }
    
    R3 vec_u, vec_v, vec_w, center, ptTet;
    nb_vertices = nb_vert_surf;
    nb_faces = nb_tri + nb_quad;
    int nb = 0, nb_elt = 0;
    for (int i = 0; i < mesh_surf.GetNbBoundaryRef(); i++)
      {
	int nb_vert = mesh_surf.BoundaryRef(i).GetNbVertices();
	// on copie l'ancienne face
	mesh_pyramid.BoundaryRef(i) = mesh_surf.BoundaryRef(i);
	
	if (nb_vert == 3)
	  {
	    mesh_tri.BoundaryRef(nb) = mesh_surf.BoundaryRef(i);
	    nb++;
	  }
	else
	  {
	    int n1 = mesh_surf.BoundaryRef(i).numVertex(0);
	    int n2 = mesh_surf.BoundaryRef(i).numVertex(1);
	    int n3 = mesh_surf.BoundaryRef(i).numVertex(2);
	    int n4 = mesh_surf.BoundaryRef(i).numVertex(3);
	    	    
	    center = 0.25*(mesh_surf.Vertex(n1) + mesh_surf.Vertex(n2)
                           + mesh_surf.Vertex(n3) + mesh_surf.Vertex(n4));
	    vec_u = mesh_surf.Vertex(n2) - mesh_surf.Vertex(n1);
	    vec_v = mesh_surf.Vertex(n4) - mesh_surf.Vertex(n1);
	    
	    TimesProd(vec_u, vec_v, vec_w);
	    Mlt(1.0/sqrt(Norm2(vec_w)), vec_w);
	    
	    // on chope le centre du tetra
	    int ne = mesh_tet.BoundaryRef(i).numElement(0);
	    for (int k = 0; k < 4; k++)
	      {
		int nv = mesh_tet.Element(ne).numVertex(k);
		if ((nv != n1) && (nv != n2) && (nv != n3) && (nv != n4))
		  ptTet = mesh_tet.Vertex(nv);
	      }
	    
	    ptTet -= mesh_surf.Vertex(n1);
	    // on change le signe si besoin
	    if (DotProd(ptTet, vec_w) < 0)
	      vec_w *= -1.0;
	    
	    // on cree l'apex de la pyramide
	    mesh_pyramid.Vertex(nb_vertices) = center + height_pyramid*vec_w;
	    mesh_tri.Vertex(nb_vertices) = center + height_pyramid*vec_w;
	    
	    // et les 4 faces triangulaires
	    mesh_pyramid.BoundaryRef(nb_faces).InitTriangular(n1, n2, nb_vertices, 0);
	    mesh_pyramid.BoundaryRef(nb_faces+1).InitTriangular(n2, n3, nb_vertices, 0);
	    mesh_pyramid.BoundaryRef(nb_faces+2).InitTriangular(n3, n4, nb_vertices, 0);
	    mesh_pyramid.BoundaryRef(nb_faces+3).InitTriangular(n1, n4, nb_vertices, 0);

	    mesh_tri.BoundaryRef(nb).InitTriangular(n1, n2, nb_vertices, new_ref);
	    mesh_tri.BoundaryRef(nb+1).InitTriangular(n2, n3, nb_vertices, new_ref);
	    mesh_tri.BoundaryRef(nb+2).InitTriangular(n3, n4, nb_vertices, new_ref);
	    mesh_tri.BoundaryRef(nb+3).InitTriangular(n1, n4, nb_vertices, new_ref);
	    
	    // la pyramide
	    mesh_pyramid.Element(nb_elt).InitPyramidal(n1, n2, n3, n4, nb_vertices, new_ref);
	    
	    nb_faces += 4; nb += 4;
	    nb_vertices++; nb_elt++;
	  }
      }
    
    mesh_tri.FindConnectivity();
    mesh_pyramid.ReorientElements();
    mesh_pyramid.FindConnectivity();
  }
  
    
  // creation d'un maillage tetraedrique a partir d'un maillage triangulaire en utilisant Ghs3D
  void CreateTetrahedralMesh(Mesh<Dimension3>& mesh_surf, Mesh<Dimension3>& mesh_tet)
  {
    std::remove("tri.mesh");
    std::remove("tri.noboiteb"); std::remove("tri.bb");
    std::remove("tet.meshb"); std::remove("don_ghs");
    
    // si on a des quadrangles dans le maillage de peau, on rajoute
    // une pyramide sur chaque quadrangle pour n'avoir que des triangles
    bool pyramid_layer = false;
    Mesh<Dimension3> mesh_pyramid, mesh_tri;
    int new_ref(0);
    if (mesh_surf.GetNbQuadranglesRef() > 0)
      {
	new_ref = mesh_surf.GetNewReference();
	pyramid_layer = true;
	CreatePyramidalLayer(mesh_surf, mesh_tri, mesh_pyramid, 0.25, new_ref);
	mesh_tri.Write("tri.mesh");
      }
    else
      {
	mesh_surf.Write("tri.mesh");
      }

    // on appelle Ghs3D pour avoir les tetras
    ofstream file_ghs("don_ghs");
    file_ghs << "10" << endl;
    file_ghs << "tri" << endl;
    file_ghs << "0" << endl;
    
    string command_gmsh = string("ghs3d < don_ghs > sortie");
    system(command_gmsh.data());
    
    command_gmsh = string("noboiteb2meshb tri.noboiteb tet.meshb > sortie");
    system(command_gmsh.data());
    
    mesh_tet.Read("tet.meshb");    
    
    // on ecrase les premiers sommets parce que Ghs3d marche en simple precision
    // et donc on recupere les sommets double precision a partir du maillage de peau initial
    if (pyramid_layer)
      {
	for (int i = 0; i < mesh_tri.GetNbVertices(); i++)
	  mesh_tet.Vertex(i) = mesh_tri.Vertex(i);
	
	// on remet les bords
	mesh_tet.ReallocateBoundariesRef(mesh_tri.GetNbBoundaryRef());
	for (int i = 0; i < mesh_tri.GetNbBoundaryRef(); i++)
	  mesh_tet.BoundaryRef(i) = mesh_tri.BoundaryRef(i);
      }
    else
      {
	for (int i = 0; i < mesh_surf.GetNbVertices(); i++)
	  mesh_tet.Vertex(i) = mesh_surf.Vertex(i);
	
	// on remet les bords
	mesh_tet.ReallocateBoundariesRef(mesh_surf.GetNbBoundaryRef());
	for (int i = 0; i < mesh_surf.GetNbBoundaryRef(); i++)
	  mesh_tet.BoundaryRef(i) = mesh_surf.BoundaryRef(i);
      }
    
    mesh_tet.FindConnectivity();
    
    // on rajoute les pyramides si necessaire
    if (pyramid_layer)
      {
	mesh_tet.AppendMesh(mesh_pyramid);
	mesh_tet.RemoveReference(new_ref);
      }
    
    std::remove("tri.mesh");
    std::remove("tri.noboiteb"); std::remove("tri.bb");
    std::remove("tet.meshb"); std::remove("don_ghs");
  }
  
  
  // creation d'un maillage quadrangulaire
  void CreateSquareMesh(int nb_points_x, int nb_points_y, const R3& ptA, const R3& ptB,
                        const R3& ptC, const R3& ptD, Mesh<Dimension3>& mesh)
  {
    Mesh<Dimension2> mesh2d;
    int ref = 1;
    TinyVector<int, 4> ref_boundary;
    ref_boundary.Fill(1);
    mesh2d.CreateRegularMesh(R2(0,0), R2(1,1), TinyVector<int, 2>(nb_points_x, nb_points_y),
			     ref, ref_boundary, mesh2d.QUADRILATERAL_MESH);
    
    mesh.ReallocateVertices(mesh2d.GetNbVertices());
    for (int i = 0; i < mesh2d.GetNbVertices(); i++)
      {
	Real_wp x = mesh2d.Vertex(i)(0);
	Real_wp y = mesh2d.Vertex(i)(1);
	R3 point = (1.0-x)*(1.0-y)*ptA + x*(1.0-y)*ptB + x*y*ptC + (1.0-x)*y*ptD;
	mesh.Vertex(i) = point;
      }
    
    mesh.ReallocateBoundariesRef(mesh2d.GetNbElt());
    IVect num(4);
    for (int i = 0; i < mesh2d.GetNbElt(); i++)
      {
	for (int k = 0; k < 4; k++)
	  num(k) = mesh2d.Element(i).numVertex(k);
	
	mesh.BoundaryRef(i).Init(num, 1);
      }
    
    mesh.FindConnectivity();
  }
  
  
  // creation d'un maillage triangulaire
  void CreateTriangularMesh(int n1, int n2, int n3, int n4,
			    const R3& ptA, const R3& ptB, const R3& ptC,
                            const R3& ptD, Mesh<Dimension3>& mesh)
  {
    ofstream file_out("tmp_geo.geo");
    file_out.precision(15);
    Real_wp dx;
    dx = 0.25*(ptA.Distance(ptB) + ptB.Distance(ptC) + ptC.Distance(ptD) + ptA.Distance(ptD));
    file_out << "lc = " << dx << ";" << endl << endl;
    file_out << "Point(1) = {" << ptA(0) << ", " << ptA(1) << ", " << ptA(2) << ", lc};" << endl;
    file_out << "Point(2) = {" << ptB(0) << ", " << ptB(1) << ", " << ptB(2) << ", lc};" << endl;
    file_out << "Point(3) = {" << ptC(0) << ", " << ptC(1) << ", " << ptC(2) << ", lc};" << endl;
    file_out << "Point(4) = {" << ptD(0) << ", " << ptD(1) 
             << ", " << ptD(2) << ", lc};" << endl << endl;
    
    file_out << "Line(1) = {1, 2};" << endl;
    file_out << "Line(2) = {2, 3};" << endl;
    file_out << "Line(3) = {3, 4};" << endl;
    file_out << "Line(4) = {4, 1};" << endl << endl;
    
    file_out << "Transfinite Line{1} = " << n1 <<";" << endl;
    file_out << "Transfinite Line{2} = " << n2 <<";" << endl;
    file_out << "Transfinite Line{3} = " << n3 <<";" << endl;
    file_out << "Transfinite Line{4} = " << n4 <<";" << endl;
    
    file_out << "Line Loop(1) = {1, 2, 3, 4};" << endl;
    file_out << "Plane Surface(1) = {1};" << endl;
    file_out << "Physical Surface(1) = {1};" << endl;
    
    file_out.close();
    
    string command_gmsh("gmsh -o plan0.msh tmp_geo.geo -3 > sortie");
    system(command_gmsh.data());
    
    mesh.Read("plan0.msh");
    std::remove("tmp_geo.geo");
    std::remove("plan0.msh");
  }

  int GetEdgeNumber(int i, int j)
  {
    switch (i)
      {
      case 0:
	{
	  if (j == 1)
	    return 0;
	  else
	    return 3;
	}
	break;
      case 1:
	{
	  if (j == 2)
	    return 1;
	  else
	    return 0;
	}
	break;
      case 2:
	{
	  if (j == 3)
	    return 2;
	  else
	    return 1;
	}
	break;
      case 3:
	{
	  if (j == 0)
	    return 3;
	  else
	    return 2;
	}
	break;    
      }

    return -1;
  }
  
  void CreateTransitionLayer(Mesh<Dimension3>& mesh, Mesh<Dimension3>& mesh2, int ref)
  {
    // pas de subdivision => on le fait a la main pour bien gerer les numeros
    //VectReal_wp step(3); step(0) = 0.0; step(1) = 0.5; step(2) = 1.0;
    //mesh2.SubdivideMesh(step);

    int Ne = mesh2.GetNbEdges();
    
    // on part d'un element qui est dans un coin
    int n0_face = -1, n0_pt = -1;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int nb_voisins = 0;
	TinyVector<bool, 4> point_on_edge(false, false, false, false);
	for (int j = 0; j < mesh.BoundaryRef(i).GetNbEdges(); j++)
	  {
	    int ne = mesh.BoundaryRef(i).numEdge(j);
	    if (mesh.GetEdge(ne).GetNbFaces() == 2)
	      {
		point_on_edge(j) = true;
		point_on_edge((j+1)%4) = true;
		nb_voisins++;
	      }
	  }
	
	if (nb_voisins == 2)
	  {
	    n0_face = i;
	    for (int j = 0; j < 4; j++)
	      if (!point_on_edge(j))
		n0_pt = j;
	    
	    break;
	  }	    
      }
    
    Vector<bool> black_point(mesh.GetNbVertices());
    black_point.Fill(false);

    {
      int nb_points = 0;
      int n1_face = n0_face; int inc1 = 1;
      int i = 1, j = 1; bool test_loop = true;
      int num_loc = (n0_pt+1)%4;
      int num_loc2 = (n0_pt+2)%4;
      while (test_loop)
	{
	  int nv = mesh.BoundaryRef(n1_face).numVertex((num_loc+inc1)%4);
	  int ne = mesh.BoundaryRef(n1_face).numEdge(num_loc);
          
          if ((i%2 == 1) && (j%2 == 1))
	    black_point(nv) = true;

	  if (mesh.GetEdge(ne).GetNbFaces() == 2)
	    {
              int nold = n1_face;
	      if (mesh.GetEdge(ne).numFace(0) == nold)
		n1_face = mesh.GetEdge(ne).numFace(1);
	      else
		n1_face = mesh.GetEdge(ne).numFace(0);

              int pos = mesh.BoundaryRef(n1_face).GetPositionBoundary(ne);
              bool same_orientation = true;
              if (mesh.BoundaryRef(nold).GetOrientationEdge(num_loc) ==
                  mesh.BoundaryRef(n1_face).GetOrientationEdge(pos))
                same_orientation = false;
              
	      if (!same_orientation)
                inc1 = 1 - inc1;
              
              num_loc = (pos+2)%4;
	    }
	  else
	    {
	      int ne2 = mesh.BoundaryRef(n0_face).numEdge(num_loc2);
	      j++; i = 0;
	      if (mesh.GetEdge(ne2).GetNbFaces() == 2)
		{
                  int nold = n0_face;
		  if (mesh.GetEdge(ne2).numFace(0) == nold)
		    n0_face = mesh.GetEdge(ne2).numFace(1);
		  else
		    n0_face = mesh.GetEdge(ne2).numFace(0);
                  
                  int pos = mesh.BoundaryRef(n0_face).GetPositionBoundary(ne2);
                  num_loc2 = (pos+2)%4;
                  num_loc = (num_loc2+1)%4;
                  int ne3 = mesh.BoundaryRef(n0_face).numEdge(num_loc);
                  if (mesh.GetEdge(ne3).GetNbFaces() == 1)
                    num_loc = (num_loc+2)%4;
                  
                  int nv2 = mesh.BoundaryRef(n0_face).numVertex((num_loc+inc1)%4);
                  bool good_point = true;
                  for (int k = 0; k < 4; k++)
                    if (mesh.BoundaryRef(nold).numVertex(k) == nv2)
                      good_point = false;
                  
                  if (!good_point)
                    inc1 = 1-inc1;
                  
                  n1_face = n0_face;
                  
		  /*n1_face = n0_face;
		  for (int k = 0; k < 4; k++)
		    if (mesh.BoundaryRef(n1_face).numVertex(k) == nv)
		      n0_pt = k;

                  if (same_orientation)
                    num_loc2 = (n0_pt+1)%4;
                  else
                    num_loc2 = (n0_pt+3)%4;
                  
		  num_loc = (num_loc2+1)%4;
		  n1_pt = n0_pt;*/
		}
	      else
		test_loop = false;
	    }
	  
	  i++; nb_points++;
	}
    }
    
    mesh.ref_vertices.Reallocate(mesh.GetNbVertices());
    mesh.ref_vertices.Zero();
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      if (black_point(i))
	mesh.ref_vertices(i) = 1;

    mesh.Write("black_point.mesh");
    mesh.ref_vertices.Clear();
    
    // on rajoute les sommets
    int N = mesh.GetNbBoundaryRef();
    //int Ne = mesh.GetNbEdges();
    int Ns = mesh.GetNbVertices();
    mesh.ResizeVertices(2*Ns + Ne + 11*N);
    int off_int = 2*Ns + Ne + N;
    //mesh.AppendMesh(mesh2);
    mesh.ResizeBoundariesRef(5*N);
    mesh.ReallocateElements(7*N);
    VectBool VertexCreated(mesh.GetNbVertices());
    VertexCreated.Fill(false);
    for (int i = 0; i < N; i++)
      {
	VectR3 s(4);
	
	int n0 = mesh.BoundaryRef(i).numVertex(0);
	int n1 = mesh.BoundaryRef(i).numVertex(1);
	int n2 = mesh.BoundaryRef(i).numVertex(2);
	int n3 = mesh.BoundaryRef(i).numVertex(3);
	s(0) = mesh.Vertex(n0); s(1) = mesh.Vertex(n1);
	s(2) = mesh.Vertex(n2); s(3) = mesh.Vertex(n3);

	IVect permut(4); permut.Fill();
	Sort(s, permut); permut(1) = permut(2);
	if (permut(1) == (permut(0) + 2)%4)
	  permut(1) = permut(3);
	
	if (permut(1) == (permut(0) + 1)%4)
	  {
	    permut(2) = (permut(0)+2)%4;
	    permut(3) = (permut(0)+3)%4;
	  }
	else if (permut(1) == (permut(0) + 3)%4)
	  {
	    permut(2) = (permut(0)+2)%4;
	    permut(3) = (permut(0)+1)%4;
	  }
	else
	  {
	    cout << "Impossible" << endl;
	    abort();
	  }

	n0 = mesh.BoundaryRef(i).numVertex(permut(0));
	n1 = mesh.BoundaryRef(i).numVertex(permut(1));
	n2 = mesh.BoundaryRef(i).numVertex(permut(2));
	n3 = mesh.BoundaryRef(i).numVertex(permut(3));

	s(0) = mesh.Vertex(n0);
	s(1) = mesh.Vertex(n1);
	s(2) = mesh.Vertex(n2);
	s(3) = mesh.Vertex(n3);

	int ne0 = GetEdgeNumber(permut(0), permut(1));
	int ne1 = GetEdgeNumber(permut(1), permut(2));
	int ne2 = GetEdgeNumber(permut(2), permut(3));
	int ne3 = GetEdgeNumber(permut(3), permut(0));

	int p0 = n0 + Ns;
	int p1 = n1 + Ns;
	int p2 = n2 + Ns;
	int p3 = n3 + Ns;
	
	 int p4 = mesh2.BoundaryRef(i).numEdge(ne0) + 2*Ns;
	 int p5 = mesh2.BoundaryRef(i).numEdge(ne1) + 2*Ns;
	 int p6 = mesh2.BoundaryRef(i).numEdge(ne2) + 2*Ns;
	 int p7 = mesh2.BoundaryRef(i).numEdge(ne3) + 2*Ns;

	 int p8 = i + 2*Ns + Ne;
	 VectR3 A(9), translat(9);
	 if (!VertexCreated(p0))
	   {
	     mesh.Vertex(p0) = mesh2.Vertex(p0-Ns);	     
	     VertexCreated(p0) = true;
	   }

	 A(0) = mesh.Vertex(p0-Ns);
	 translat(0) = A(0) - mesh.Vertex(p0);
	 
	 if (!VertexCreated(p1))
	   {	     
	     mesh.Vertex(p1) = mesh2.Vertex(p1-Ns);
	     VertexCreated(p1) = true;
	   }

	 A(1) = mesh.Vertex(p1-Ns);
	 translat(1) = A(1) - mesh.Vertex(p1);
	 
	 if (!VertexCreated(p2))
	   {
	     mesh.Vertex(p2) = mesh2.Vertex(p2-Ns);	     
	     VertexCreated(p2) = true;
	   }

	 A(2) = mesh.Vertex(p2-Ns);
	 translat(2) = A(2) - mesh.Vertex(p2);
	 
	 if (!VertexCreated(p3))
	   {
	     mesh.Vertex(p3) = mesh2.Vertex(p3-Ns);
	     VertexCreated(p3) = true;
	   }

	 A(3) = mesh.Vertex(p3-Ns);
	 translat(3) = A(3) - mesh.Vertex(p3);
	 
	 if (!VertexCreated(p4))
	   {
	     mesh.Vertex(p4) = 0.5*(mesh.Vertex(p0) + mesh.Vertex(p1));	     
	     VertexCreated(p4) = true;
	   }

	 A(4) = 0.5*(mesh.Vertex(p0-Ns) + mesh.Vertex(p1-Ns));
	 translat(4) = A(4) - mesh.Vertex(p4);

	 if (!VertexCreated(p5))
	   {
	     mesh.Vertex(p5) = 0.5*(mesh.Vertex(p1) + mesh.Vertex(p2));
	     VertexCreated(p5) = true;
	   }
	 
	 A(5) = 0.5*(mesh.Vertex(p1-Ns) + mesh.Vertex(p2-Ns));
	 translat(5) = A(5) - mesh.Vertex(p5);

	 if (!VertexCreated(p6))
	   {
	     mesh.Vertex(p6) = 0.5*(mesh.Vertex(p2) + mesh.Vertex(p3));
	     VertexCreated(p6) = true;
	   }

	 A(6) = 0.5*(mesh.Vertex(p2-Ns) + mesh.Vertex(p3-Ns));
	 translat(6) = A(6) - mesh.Vertex(p6);

	 if (!VertexCreated(p7))
	   {
	     mesh.Vertex(p7) = 0.5*(mesh.Vertex(p3) + mesh.Vertex(p0));
	     VertexCreated(p7) = true;
	   }

	 A(7) = 0.5*(mesh.Vertex(p3-Ns) + mesh.Vertex(p0-Ns));
	 translat(7) = A(7) - mesh.Vertex(p7);

	 A(8) = 0.25*(mesh.Vertex(p0-Ns) + mesh.Vertex(p1-Ns) + mesh.Vertex(p2-Ns) + mesh.Vertex(p3-Ns));
	 mesh.Vertex(p8) = 0.25*(mesh.Vertex(p0) + mesh.Vertex(p1) + mesh.Vertex(p2) + mesh.Vertex(p3));
	 translat(8) = A(8) - mesh.Vertex(p8);

	 mesh.BoundaryRef(N+i).InitQuadrangular(p0, p4, p8, p7, 1);
	 mesh.BoundaryRef(2*N + i).InitQuadrangular(p4, p1, p5, p8, 1);
	 mesh.BoundaryRef(3*N + i).InitQuadrangular(p7, p8, p6, p3, 1);
	 mesh.BoundaryRef(4*N + i).InitQuadrangular(p8, p5, p2, p6, 1);

	 // 1er cas	
	 if (black_point(n2))
	   {
	     // points du plan -0.25
	     mesh.Vertex(off_int) = A(2) - 0.25*translat(2);
	     mesh.Vertex(off_int+1) = A(1) - 0.25*translat(1);
	     mesh.Vertex(off_int+2) = A(4) - 0.25*translat(4);
	     mesh.Vertex(off_int+3) = A(6) - 0.25*translat(6);
	    
	     // points du plan -0.75
	     mesh.Vertex(off_int+4) = A(2) - 0.75*translat(2);
	     mesh.Vertex(off_int+5) = A(5) - 0.75*translat(5);
	     mesh.Vertex(off_int+6) = A(8) - 0.75*translat(8);
	     mesh.Vertex(off_int+7) = A(6) - 0.75*translat(6);

	     // points du plan -0.5
	     mesh.Vertex(off_int+8) = A(3) - 0.5*translat(3);
	     mesh.Vertex(off_int+9) = A(7) - 0.5*translat(7);
	     
	     mesh.Element(i).InitHexahedral(off_int, off_int+1, off_int+2, off_int+3, n2, n1, n0, n3, ref);
	     mesh.Element(N+i).InitHexahedral(off_int+4, off_int+5, off_int+6, off_int+7,off_int,off_int+1, off_int+2, off_int+3, ref);
	     mesh.Element(2*N+i).InitHexahedral(off_int+4, off_int+5, off_int+6, off_int+7,p2,p5,p8,p6, ref);
	     mesh.Element(3*N+i).InitHexahedral(off_int+5, off_int+1, off_int+2, off_int+6,p5,p1,p4,p8, ref);
	     mesh.Element(4*N+i).InitHexahedral(off_int+7, off_int+6, off_int+9, off_int+8,p6,p8,p7,p3, ref);
	     mesh.Element(5*N+i).InitHexahedral(off_int+7, off_int+6, off_int+9, off_int+8,off_int+3,off_int+2,n0,n3, ref);
	     mesh.Element(6*N+i).InitHexahedral(off_int+6, off_int+2, n0, off_int+9,p8,p4,p0,p7, ref);
	     
	     off_int += 10;
	   }
	
	 // 2eme cas	
	 if (black_point(n1))
	   {
	     // points du plan -0.25
	     mesh.Vertex(off_int) = A(2) - 0.25*translat(2);
	     mesh.Vertex(off_int+1) = A(1) - 0.25*translat(1);
	     mesh.Vertex(off_int+2) = A(4) - 0.25*translat(4);
	     mesh.Vertex(off_int+3) = A(6) - 0.25*translat(6);
	     
	     // points du plan -0.75
	     mesh.Vertex(off_int+4) = A(5) - 0.75*translat(5);
	     mesh.Vertex(off_int+5) = A(1) - 0.75*translat(1);
	     mesh.Vertex(off_int+6) = A(4) - 0.75*translat(4);
	     mesh.Vertex(off_int+7) = A(8) - 0.75*translat(8);

	     // points du plan -0.5
	     mesh.Vertex(off_int+8) = A(0) - 0.5*translat(0);
	     mesh.Vertex(off_int+9) = A(7) - 0.5*translat(7);
	    
	     mesh.Element(i).InitHexahedral( off_int, off_int+1, off_int+2, off_int+3, n2, n1, n0, n3, ref);
	     mesh.Element(N+i).InitHexahedral( off_int, off_int+4, off_int+7, off_int+3, p2, p5, p8, p6, ref);
	     mesh.Element(2*N+i).InitHexahedral( off_int+4, off_int+5, off_int+6, off_int+7, p5, p1, p4, p8, ref);
	     mesh.Element(3*N+i).InitHexahedral( off_int, off_int+1, off_int+2, off_int+3, off_int+4, off_int+5, off_int+6, off_int+7, ref);
	     mesh.Element(4*N+i).InitHexahedral( off_int+2, off_int+3, n3, n0, off_int+6, off_int+7, off_int+9, off_int+8, ref);
	     mesh.Element(5*N+i).InitHexahedral( off_int+6, off_int+7, off_int+9, off_int+8, p4, p8, p7, p0, ref);
	     mesh.Element(6*N+i).InitHexahedral( off_int+7, off_int+9, n3, off_int+3, p8, p7, p3, p6, ref);
	     
	     off_int += 10;
	   }
	 
	 //3eme cas
	 if (black_point(n0))
	   {
	     // points du plan -0.25
	     mesh.Vertex(off_int) = A(0) - 0.25*translat(0);
	     mesh.Vertex(off_int+1) = A(3) - 0.25*translat(3);
	     mesh.Vertex(off_int+2) = A(6) - 0.25*translat(6);
	     mesh.Vertex(off_int+3) = A(4) - 0.25*translat(4);

	     // points du plan -0.75
	     mesh.Vertex(off_int+4) = A(7) - 0.75*translat(7);
	     mesh.Vertex(off_int+5) = A(0) - 0.75*translat(0);
	     mesh.Vertex(off_int+6) = A(4) - 0.75*translat(4);
	     mesh.Vertex(off_int+7) = A(8) - 0.75*translat(8);
	     
	     // points du plan -0.5
	     mesh.Vertex(off_int+8) = A(1) - 0.5*translat(1);
	     mesh.Vertex(off_int+9) = A(5) - 0.5*translat(5);
	     
	     mesh.Element(i).InitHexahedral( off_int,off_int+1, off_int+2, off_int+3,n0,n3,n2,n1, ref);
	     mesh.Element(N+i).InitHexahedral( off_int+1,off_int+2, off_int+7, off_int+4,p3,p6,p8,p7, ref);
	     mesh.Element(2*N+i).InitHexahedral( off_int+4,off_int+7, off_int+6, off_int+5,p7,p8,p4,p0, ref);
	     mesh.Element(3*N+i).InitHexahedral( off_int+4,off_int+7, off_int+6, off_int+5,off_int+1,off_int+2, off_int+3, off_int, ref);
	     mesh.Element(4*N+i).InitHexahedral( off_int+2,n2,off_int+9,off_int+7,p6,p2,p5,p8, ref);
	     mesh.Element(5*N+i).InitHexahedral( off_int+2,n2,n1, off_int+3,off_int+7,off_int+9,off_int+8,off_int+6, ref);
	     mesh.Element(6*N+i).InitHexahedral( off_int+9,off_int+8,off_int+6,off_int+7,p5,p1,p4,p8, ref);
	    
	     off_int += 10;
	   }
	 
	 //4eme cas
	 if (black_point(n3))
	   {
	     // points du plan -0.25
	     mesh.Vertex(off_int) = A(0) - 0.25*translat(0);
	     mesh.Vertex(off_int+1) = A(3) - 0.25*translat(3);
	     mesh.Vertex(off_int+2) = A(6) - 0.25*translat(6);
	     mesh.Vertex(off_int+3) = A(4) - 0.25*translat(4);

	     // points du plan -0.75
	     mesh.Vertex(off_int+4) = A(3) - 0.75*translat(3);
	     mesh.Vertex(off_int+5) = A(7) - 0.75*translat(7);
	     mesh.Vertex(off_int+6) = A(8) - 0.75*translat(8);
	     mesh.Vertex(off_int+7) = A(6) - 0.75*translat(6);

	     // points du plan -0.5
	     mesh.Vertex(off_int+8) = A(2) - 0.5*translat(2);
	     mesh.Vertex(off_int+9) = A(5) - 0.5*translat(5);
	     mesh.Element(i).InitHexahedral( off_int,off_int+1, off_int+2, off_int+3,n0,n3,n2,n1, ref);
	     mesh.Element(N+i).InitHexahedral( off_int+4,off_int+7, off_int+6, off_int+5,p3,p6,p8,p7, ref);
	     mesh.Element(2*N+i).InitHexahedral( off_int+5,off_int+6, off_int+3, off_int,p7,p8,p4,p0, ref);
	     mesh.Element(3*N+i).InitHexahedral( off_int+4,off_int+7, off_int+6, off_int+5,off_int+1,off_int+2, off_int+3, off_int, 1);
	     mesh.Element(4*N+i).InitHexahedral( off_int+7,off_int+8,off_int+9,off_int+6,p6,p2,p5,p8, ref);
	     mesh.Element(5*N+i).InitHexahedral( off_int+7,off_int+8,off_int+9,off_int+6,off_int+2,n2,n1,off_int+3, ref);
	     mesh.Element(6*N+i).InitHexahedral( off_int+9,n1,off_int+3,off_int+6,p5,p1,p4,p8, ref);
	     
	     off_int += 10;
	   }
      }

    mesh.SortBoundariesRef();
    mesh.FindConnectivity();
    mesh.RemoveDuplicateVertices();
    mesh.AddBoundaryFaces();
    mesh.FindConnectivity();
  }
  
}

#define MONTJOIE_FILE_MESH_GENERATION_CXX
#endif

