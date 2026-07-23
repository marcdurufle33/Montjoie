#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

#include "Mesh/MeshGeneration.cxx"

using namespace Montjoie;

void CreateHexahedralMeshHelice(const Real_wp& alpha, const Real_wp& delta1, const Real_wp& delta2, const Real_wp& delta3,
				int nb_subdiv_circle, int nb_subdiv_teta, const Real_wp& ymax, const Real_wp& zmax, Mesh<Dimension3>& mesh)
{
  int order = mesh.GetGeometryOrder();
  Vector<R2> PtsQuad = mesh.GetQuadrilateralReferenceElement().PointsNodalND();
  
  int nb_points = nb_subdiv_teta + 1 + 3*nb_subdiv_circle*(nb_subdiv_teta+1);
  int nb_faces = 6*nb_subdiv_circle + 3*nb_subdiv_circle*nb_subdiv_teta;
  int nb_elt = 3*nb_subdiv_circle*nb_subdiv_teta;
  int offset_point = 2 + 3*nb_subdiv_circle*(nb_subdiv_teta+1);
  mesh.ReallocateVertices(nb_points);
  mesh.ReallocateBoundariesRef(nb_faces);
  mesh.ReallocateElements(nb_elt);
  mesh.ResizeNbReferences(20);
  Vector<VectR3> PtsCurvedFace(nb_faces);

  R3 ex(1, 0, 0), ey(0, 1, 0), ez(0, 0, 1), ux, uy, uz;
  nb_points = 0; nb_faces = 0; nb_elt = 0;
  R3 ptA, ptB;
  IVect num(8);
  
  // axis on the section
  ux = ex;
  uz = cos(alpha)*ey + sin(alpha)*ez;
  TimesProd(uz, ux, uy);
  
  // loop over teta-subdivisions of a wire
  int nc = nb_subdiv_circle;
  for (int nt = 0; nt <= nb_subdiv_teta; nt++)
    {
      Real_wp teta = Real_wp(nt)*ymax/nb_subdiv_teta;
      Real_wp z = Real_wp(nt)*zmax/nb_subdiv_teta;
      // center of the small section
      ptA.Init(0, teta, z);
            
      // section pour z = 0
      if (nt == 0)
        {
	  mesh.Vertex(nb_points) = ptA;
          
          for (int i = 0; i < nb_subdiv_circle; i++)
	    {
	      int ip2 = i+2;
	      if (i == nb_subdiv_circle-1)
		ip2 = 1;
	      
	      mesh.BoundaryRef(nb_faces).InitTriangular(nb_points, nb_points + (i+1), nb_points + ip2, 1);
	      
	      mesh.BoundaryRef(nb_faces+1).InitQuadrangular(nb_points + (i+1), nb_points + ip2,
							  nb_points + nc + ip2, nb_points + nc + (i+1), 1);
	      
	      mesh.BoundaryRef(nb_faces+2).InitQuadrangular(nb_points + nc + (i+1), nb_points+nc+ ip2,
							  nb_points + 2*nc + ip2, nb_points + 2*nc + (i+1), 1);
	      nb_faces += 3;
	    }
          
	  nb_points++;
        }
      
      for (int i = 0; i < nb_subdiv_circle; i++)
        {
          Real_wp phi = 2.0*pi_wp*Real_wp(i)/nb_subdiv_circle;
	  Real_wp lambda = Real_wp(nt)/nb_subdiv_teta;
	  Real_wp z0 = delta1*sin(phi)*cos(alpha);
	  Real_wp z1 = zmax/ymax*(1.0-delta1*sin(phi)*sin(alpha));
	  z = (1.0-lambda)*z0 + lambda*z1;
	  teta = z*ymax/zmax;
	  ptA.Init(0, teta, z);
          ptB = ptA + cos(phi)*delta1*ux + sin(phi)*delta1*uy;
	  // DISP(uy); DISP(cos(alpha)); DISP(z0); DISP(nt); DISP(lambda); DISP(z); DISP(teta); DISP(ptB);
	  mesh.Vertex(nb_points + i) = ptB;

	  z0 = delta2*sin(phi)*cos(alpha);
	  z1 = zmax/ymax*(1.0-delta2*sin(phi)*sin(alpha));
	  z = (1.0-lambda)*z0 + lambda*z1;
	  teta = z*ymax/zmax;
	  ptA.Init(0, teta, z);
          ptB = ptA + cos(phi)*delta2*ux + sin(phi)*delta2*uy;
	  mesh.Vertex(nb_points + nc + i) = ptB;

	  z0 = delta3*sin(phi)*cos(alpha);
	  z1 = zmax/ymax*(1.0-delta3*sin(phi)*sin(alpha));
	  z = (1.0-lambda)*z0 + lambda*z1;
	  teta = z*ymax/zmax;
	  ptA.Init(0, teta, z);
          ptB = ptA + cos(phi)*delta3*ux + sin(phi)*delta3*uy;
	  mesh.Vertex(nb_points + 2*nc + i) = ptB;
	  
	  if (nt < nb_subdiv_teta)
	    {
	      int ip1 = i+1;
	      if (i == nb_subdiv_circle-1)
		ip1 = 0;
	      
	      mesh.BoundaryRef(nb_faces).InitQuadrangular(nb_points + i, nb_points + ip1,
							  nb_points + ip1 + 3*nc, nb_points + i + 3*nc, 2);

	      mesh.BoundaryRef(nb_faces+1).InitQuadrangular(nb_points + nc+i, nb_points + nc+ip1,
							    nb_points + ip1 + 4*nc, nb_points + i + 4*nc, 3);
	      
	      mesh.BoundaryRef(nb_faces+2).InitQuadrangular(nb_points + 2*nc + i, nb_points + 2*nc + ip1,
							    nb_points + ip1 + 5*nc, nb_points + i + 5*nc, 4);
	      
	      num(0) = nb_points + i; num(1) = nb_points + ip1;
	      num(2) = nb_points + ip1 + 3*nc; num(3) = nb_points + i + 3*nc;
	      num(4) = nb_points + nc+i; num(5) = nb_points + nc+ip1;
	      num(6) = nb_points + ip1 + 4*nc; num(7) = nb_points + i + 4*nc;	      
	      mesh.Element(nb_elt).InitHexahedral(num, 2);
	      
	      num(0) = num(4); num(1) = num(5);
	      num(2) = num(6); num(3) = num(7);
	      num(4) += nc;  num(5) += nc;  num(6) += nc;  num(7) += nc; 
	      mesh.Element(nb_elt+1).InitHexahedral(num, 3);
	      
	      if (nt > 0)
		{
		  if (i == 0)
		    {
		      z = lambda*zmax/ymax;
		      teta = z*ymax/zmax;
		      ptA.Init(0, teta, z);
		      mesh.Vertex(offset_point) = ptA;
		    }
		  
		  if (nt < nb_subdiv_teta - 1)
		    mesh.Element(nb_elt+2).InitWedge(offset_point, nb_points+i, nb_points + ip1,
						   offset_point+1, nb_points+3*nc+i, nb_points + 3*nc + ip1, 2);
		  else
		    mesh.Element(nb_elt+2).InitWedge(offset_point, nb_points+i, nb_points + ip1,
						   nb_points+6*nc, nb_points+3*nc+i, nb_points + 3*nc + ip1, 2);		  		  
		}
	      else
		{
		  mesh.Element(nb_elt+2).InitWedge(0, nb_points+i, nb_points + ip1,
						 offset_point, nb_points+3*nc+i, nb_points + 3*nc + ip1, 2);
		}
	      
	      if (order > 1)
		{
		  PtsCurvedFace(nb_faces+1).Reallocate(PtsQuad.GetM());
		  for (int j = 0; j < PtsQuad.GetM(); j++)
		    {
		      // points on the curved boundary
		      Real_wp xj = PtsQuad(j)(0);
		      Real_wp yj = PtsQuad(j)(1);
		      
		      phi = 2.0*pi_wp*((1.0-xj)*i + xj*(i+1))/nb_subdiv_circle;
		      lambda = ((1.0-yj)*nt + yj*(nt+1))/nb_subdiv_teta;
		      z0 = delta2*sin(phi)*cos(alpha);
		      z1 = zmax/ymax*(1.0-delta2*sin(phi)*sin(alpha));
		      z = (1.0-lambda)*z0 + lambda*z1;
		      teta = z*ymax/zmax;
		      ptA.Init(0, teta, z);
		      ptB = ptA + cos(phi)*delta2*ux + sin(phi)*delta2*uy;
		      PtsCurvedFace(nb_faces+1)(j) = ptB;
		    }
		}
	      
	      nb_elt += 3;
	      nb_faces += 3;
	    }	  
	}
      
      if (nt > 0)
	offset_point++;
      
      nb_points += 3*nc;
      
      teta = Real_wp(nt)*ymax/nb_subdiv_teta;
      z = Real_wp(nt)*zmax/nb_subdiv_teta;
      // center of the small section
      ptA.Init(0, teta, z);

      // section z = zmax
      if (nt == nb_subdiv_teta)
        {
	  mesh.Vertex(nb_points) = ptA;
          
          for (int i = 0; i < nb_subdiv_circle; i++)
	    {
	      int ip2 = i+2;
	      if (i == nb_subdiv_circle-1)
		ip2 = 1;
                
	      mesh.BoundaryRef(nb_faces).InitQuadrangular(nb_points - nc - (i+1), nb_points-nc- ip2,
							  nb_points - 2*nc - ip2, nb_points - 2*nc - (i+1), 1);
	      
	      mesh.BoundaryRef(nb_faces+1).InitQuadrangular(nb_points - (i+1), nb_points - ip2,
							  nb_points - nc - ip2, nb_points - nc - (i+1), 1);
	      
	      mesh.BoundaryRef(nb_faces+2).InitTriangular(nb_points, nb_points - 2*nc - (i+1), nb_points - 2*nc - ip2, 1);
	      nb_faces += 3;
	    }
          
          nb_points++;
        }
        

    }
  
  //DISP(nb_points); DISP(mesh.GetNbVertices());
  //DISP(nb_faces); DISP(mesh.GetNbBoundaryRef());
  mesh.FindConnectivity();

  // putting curved boundaries
  if (order > 1)
    {
      mesh.SetCurveType(3, mesh.CURVE_FILE);
      IVect IndexEdge(mesh.GetNbEdges());
      Vector<VectR3> PtsCurvedEdge(mesh.GetNbEdges());
      IndexEdge.Fill(-1); int nb_curved_edges = 0;
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	if (PtsCurvedFace(i).GetM() > 0)
	  for (int j = 0; j < mesh.BoundaryRef(i).GetNbEdges(); j++)
	    {
	      int ne = mesh.BoundaryRef(i).numEdge(j);
	      if (IndexEdge(ne) == -1)
		{
		  IndexEdge(ne) = nb_curved_edges;
		  PtsCurvedEdge(ne).Reallocate(order-1);
		  if (mesh.BoundaryRef(i).GetOrientationEdge(j))
		    {
		      for (int k = 0; k < order-1; k++)
			PtsCurvedEdge(ne)(k) = PtsCurvedFace(i)(4+j*(order-1)+k);
		    }
		  else
		    {
		      for (int k = 0; k < order-1; k++)
			PtsCurvedEdge(ne)(order-2-k) = PtsCurvedFace(i)(4+j*(order-1)+k);
		    }
		  
		  nb_curved_edges++;
		}
	    }
      
      mesh.ReallocateEdgesRef(nb_curved_edges);
      for (int i = 0; i < mesh.GetNbEdges(); i++)
	if (IndexEdge(i) >= 0)
	  {
	    int ne = IndexEdge(i);
	    mesh.EdgeRef(ne) = mesh.GetEdge(i);
	    mesh.EdgeRef(ne).SetReference(3);
	    for (int k = 0; k < order-1; k++)
	      mesh.SetPointInsideEdge(ne, k, PtsCurvedEdge(i)(k));
	  }

      VectR3 PtsInside(PtsQuad.GetM() - 4*order);
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	if (PtsCurvedFace(i).GetM() > 0)
	  {
	    for (int k = 0; k < PtsInside.GetM(); k++)
	      PtsInside(k) = PtsCurvedFace(i)(4*order+k);
	    
	    mesh.SetPointInsideFace(i, PtsInside);
	  }
    }
  
  mesh.SortBoundariesRef();
  mesh.ReorientElements();
  mesh.FindConnectivity();
  mesh.ProjectPointsOnCurves();
}


void CreateMeshAroundCircleZ0(const Mesh<Dimension3>& mesh, const Real_wp& h, int nb_points,
			      const Real_wp& xmax, const Real_wp& ymax,  const Real_wp& zmax,
			      int nb_subdiv_circle, int nb_subdiv_teta, Mesh<Dimension3>& mesh_z0)
{
  // on extrait les points de la section z = 0
  ofstream file_geo("tmp_geo.geo"); file_geo.precision(15);
  file_geo << "lc = " << h << ";" << endl;
  file_geo << endl;
  file_geo << "Point(1) = {-"<<xmax<<", -"<<ymax<<", 0, lc}; " << endl;
  file_geo << "Point(2) = {"<<xmax<<", -"<<ymax<<", 0, lc}; " << endl;
  file_geo << "Point(3) = {"<<xmax<<", "<<ymax<<", 0, lc}; " << endl;
  file_geo << "Point(4) = {-"<<xmax<<", "<<ymax<<", 0, lc}; " << endl;
  int nc = nb_subdiv_circle;
  for (int i = 0; i < nb_subdiv_circle; i++)
    file_geo << "Point(" << i+5 << ") = {" << mesh.Vertex(i+2*nc+1)(0) << ", " << mesh.Vertex(i+2*nc+1)(1) << ", 0, lc};" << endl;  
  
  file_geo << endl;
  file_geo << "Line(1) = {1, 2};" << endl;
  file_geo << "Line(2) = {2, 3};" << endl;
  file_geo << "Line(3) = {3, 4};" << endl;
  file_geo << "Line(4) = {4, 1};" << endl;
  file_geo << "Transfinite Line{1, 2, 3, 4} = "<<nb_points<<";" << endl;
  
  for (int i = 0; i < nb_subdiv_circle; i++)
    file_geo << "Line(" << i+5 << ") = {" << i+5 << ", " << (i+1)%nb_subdiv_circle + 5 << "};" << endl;
  
  file_geo << endl;
  file_geo << "Line Loop(1) = {1, 2, 3, 4};" << endl;
  file_geo << "Line Loop(2) = {5";
  for (int i = 1; i < nb_subdiv_circle; i++)
    file_geo << ", " << i+5;
  
  file_geo << "};" << endl << endl;
  file_geo << "Physical Line(1) = {5";
  for (int i = 1; i < nb_subdiv_circle; i++)
    file_geo << ", " << i+5;

  file_geo << "};" << endl;
  file_geo << "Physical Line(2) = {1, 2, 3, 4};" << endl << endl;
  
  file_geo << "Plane Surface(1) = {1, 2};" << endl;
  file_geo << "Physical Surface(1) = {1};" << endl;
  file_geo.close();
  
  string command_gmsh("gmsh -o plan0.msh tmp_geo.geo -3 -format msh2 > sortie");
  system(command_gmsh.data());
  
  mesh_z0.Read("plan0.msh");
  std::remove("tmp_geo.geo");
  std::remove("plan0.msh");
}


void CreateMeshAroundCircleTeta1(const Mesh<Dimension3>& mesh, const Real_wp& h, int nb_points, int nb_points_z,
				 const Real_wp& xmax, const Real_wp& ymax,  const Real_wp& zmax,
				 int nb_subdiv_circle, int nb_subdiv_teta, Mesh<Dimension3>& mesh_teta1)
{
  // on extrait les points de la section teta = 1
  ofstream file_geo;
  file_geo.precision(15);
  file_geo.open("tmp_geo.geo");
  file_geo << "lc = " << h << ";" << endl;
  file_geo << endl;
  file_geo << "Point(1) = {-"<<xmax<<", "<<ymax<<", 0, lc}; " << endl;
  file_geo << "Point(2) = {"<<xmax<<", "<<ymax<<", 0, lc}; " << endl;
  file_geo << "Point(3) = {"<<xmax<<", "<<ymax<<", "<<2*zmax<<", lc}; " << endl;
  file_geo << "Point(4) = {-"<<xmax<<", "<<ymax<<", "<<2*zmax<<", lc}; " << endl;
  int nc = nb_subdiv_circle;
  int offset = 3*nb_subdiv_teta*nb_subdiv_circle + 2*nc + 1;
  for (int i = 0; i < nb_subdiv_circle; i++)
    file_geo << "Point(" << i+5 << ") = {" << mesh.Vertex(offset+i)(0) << ", " << mesh.Vertex(offset+i)(1) << ", " << mesh.Vertex(offset+i)(2) << ", lc};" << endl;  
  
  file_geo << endl;
  file_geo << "Line(1) = {1, 2};" << endl;
  file_geo << "Line(2) = {2, 3};" << endl;
  file_geo << "Line(3) = {3, 4};" << endl;
  file_geo << "Line(4) = {4, 1};" << endl;
  file_geo << "Transfinite Line{1, 3} = "<<nb_points<<";" << endl;
  file_geo << "Transfinite Line{2, 4} = "<<nb_points_z<<";" << endl;
  
  for (int i = 0; i < nb_subdiv_circle; i++)
    file_geo << "Line(" << i+5 << ") = {" << i+5 << ", " << (i+1)%nb_subdiv_circle + 5 << "};" << endl;
  
  file_geo << endl;
  file_geo << "Line Loop(1) = {1, 2, 3, 4};" << endl;
  file_geo << "Line Loop(2) = {5";
  for (int i = 1; i < nb_subdiv_circle; i++)
    file_geo << ", " << i+5;
  
  file_geo << "};" << endl << endl;
  file_geo << "Physical Line(1) = {5";
  for (int i = 1; i < nb_subdiv_circle; i++)
    file_geo << ", " << i+5;

  file_geo << "};" << endl;
  file_geo << "Physical Line(2) = {1, 2, 3, 4};" << endl << endl;
  
  file_geo << "Plane Surface(1) = {1, 2};" << endl;
  file_geo << "Physical Surface(1) = {1};" << endl;
  file_geo.close();
  
  string command_gmsh = string("gmsh -o plan1.msh tmp_geo.geo -3 -format msh2 > sortie");  
  system(command_gmsh.data());

  mesh_teta1.Read("plan1.msh");
}


void ConstructPyramidalTransition(const Mesh<Dimension3>& mesh, int ref, int ref_domain, const Real_wp& ratio,
				  Mesh<Dimension3>& support_mesh, int new_ref, Mesh<Dimension3>& mesh_pyramid)
{
  IVect ref_cond(mesh.GetNbReferences()+1); ref_cond.Fill();
  IVect IndexFace, IndexVert, NumElem, NumLoc;
  mesh.GetBoundaryMesh(ref, mesh_pyramid, IndexFace, IndexVert, NumElem, NumLoc, ref_cond, ref_domain);
  
  // on construit les pyramides
  int nb_old = mesh_pyramid.GetNbVertices();
  int nb_old_face = mesh_pyramid.GetNbBoundaryRef();
  mesh_pyramid.ResizeVertices(nb_old + mesh_pyramid.GetNbBoundaryRef());  
  mesh_pyramid.ReallocateElements(mesh_pyramid.GetNbBoundaryRef());  
  mesh_pyramid.ResizeBoundariesRef(5*mesh_pyramid.GetNbBoundaryRef());  
  R3 ptB, ux, uy, uz, ptA, center;
  for (int i = 0; i < nb_old_face; i++)
    {
      int num_elem = NumElem(i);
      int nb_vert = mesh.Element(num_elem).GetNbVertices();
      // barycenter of the element 
      center.Fill(0);
      for (int k = 0; k < nb_vert; k++)
	center += mesh.Vertex(mesh.Element(num_elem).numVertex(k));
      
      Mlt(1.0/nb_vert, center);
				 
      IVect num_supp(4), num_pyr(5);
      int nb_vert_on_supp = 0;
      for (int k = 0; k < 4; k++)
	{
	  num_pyr(k) = mesh_pyramid.BoundaryRef(i).numVertex(k);
	  int nv = support_mesh.FindVertexNumber(mesh_pyramid.Vertex(num_pyr(k)));
	  if (nv >= 0)
	    {
	      num_supp(nb_vert_on_supp) = nv;
	      nb_vert_on_supp++;
	    }
	}
      
      bool new_point = true;
      num_pyr(4) = nb_old;
      
      if (nb_vert_on_supp >= 2)
	{
	  int num_edge = support_mesh.FindEdgeWithExtremities(num_supp(0), num_supp(1));
	  int num_face = support_mesh.GetEdge(num_edge).numFace(0);	  
	  int old_point = -1;
	  for (int k = 0; k < support_mesh.Boundary(num_face).GetNbVertices(); k++)
	    {
	      int nv = support_mesh.Boundary(num_face).numVertex(k);
	      if ((nv != num_supp(0)) && (nv != num_supp(1)))
		old_point = nv;
	    }
	  
	  new_point = false;
	  ptB = support_mesh.Vertex(old_point);
	}
      
      ux = mesh_pyramid.Vertex(num_pyr(1)) - mesh_pyramid.Vertex(num_pyr(0));
      uy = mesh_pyramid.Vertex(num_pyr(3)) - mesh_pyramid.Vertex(num_pyr(0));
      TimesProd(ux, uy, ptA);
      Mlt(1.0/sqrt(Norm2(ptA)), ptA);
      
      uz = 0.25*(mesh_pyramid.Vertex(num_pyr(0)) + mesh_pyramid.Vertex(num_pyr(1))
		 + mesh_pyramid.Vertex(num_pyr(2)) + mesh_pyramid.Vertex(num_pyr(3)));
      
      // checking if ptA is well oriented
      center -= uz;
      if (DotProd(center, ptA) > 0)
	ptA *= -1.0;
      
      if (new_point)
	mesh_pyramid.Vertex(nb_old) = uz + ratio*ptA;
      else
	mesh_pyramid.Vertex(nb_old) = ptB;
      
      mesh_pyramid.Element(i).InitPyramidal(num_pyr(0), num_pyr(1), num_pyr(2),
					    num_pyr(3), num_pyr(4), ref_domain);
      
      mesh_pyramid.BoundaryRef(nb_old_face+i).InitTriangular(num_pyr(0), num_pyr(1), num_pyr(4), new_ref);
      mesh_pyramid.BoundaryRef(2*nb_old_face+i).InitTriangular(num_pyr(1), num_pyr(2), num_pyr(4), new_ref);
      mesh_pyramid.BoundaryRef(3*nb_old_face+i).InitTriangular(num_pyr(2), num_pyr(3), num_pyr(4), new_ref);
      mesh_pyramid.BoundaryRef(4*nb_old_face+i).InitTriangular(num_pyr(3), num_pyr(0), num_pyr(4), new_ref);
      
      nb_old++;
    }
  
  mesh_pyramid.ReorientElements();
  mesh_pyramid.SortBoundariesRef();
  mesh_pyramid.FindConnectivity();
}

void CreatePyramidLayer(const Mesh<Dimension3>& mesh_square, const R3& center,
			Mesh<Dimension3>& mesh_tri, Mesh<Dimension3>& mesh_pyramid, int ref_domain)
{
  mesh_pyramid.Clear();
  mesh_tri.Clear();
  
  mesh_tri.ReallocateBoundariesRef(4*mesh_square.GetNbBoundaryRef());
  mesh_tri.ReallocateVertices(mesh_square.GetNbVertices() + mesh_square.GetNbBoundaryRef());
  
  mesh_pyramid.ReallocateVertices(mesh_square.GetNbVertices() + mesh_square.GetNbBoundaryRef());
  mesh_pyramid.ReallocateElements(mesh_square.GetNbBoundaryRef());
  mesh_pyramid.ReallocateBoundariesRef(4*mesh_square.GetNbBoundaryRef());
  int nb_old = mesh_square.GetNbVertices();
  for (int i = 0; i < nb_old; i++)
    {
      mesh_pyramid.Vertex(i) = mesh_square.Vertex(i);
      mesh_tri.Vertex(i) = mesh_square.Vertex(i);
    }
  
  for (int i = 0; i < mesh_square.GetNbBoundaryRef(); i++)
    {
      int n1 = mesh_square.BoundaryRef(i).numVertex(0);
      int n2 = mesh_square.BoundaryRef(i).numVertex(1);
      int n3 = mesh_square.BoundaryRef(i).numVertex(2);
      int n4 = mesh_square.BoundaryRef(i).numVertex(3);
      Real_wp dx = mesh_square.Vertex(n1).Distance(mesh_square.Vertex(n2));
      R3 ptA = 0.25*(mesh_square.Vertex(n1) + mesh_square.Vertex(n2) + mesh_square.Vertex(n3) + mesh_square.Vertex(n4));
      Real_wp diff = center(0) - ptA(0);
      if (diff < 0)
	ptA(0) -= 0.25*dx;
      else
	ptA(0) += 0.25*dx;
      
      mesh_pyramid.Vertex(nb_old + i) = ptA;
      mesh_tri.Vertex(nb_old + i) = ptA;
      mesh_pyramid.Element(i).InitPyramidal(n1, n2, n3, n4, nb_old+i, ref_domain);
      mesh_pyramid.BoundaryRef(4*i).InitTriangular(n1, n2, nb_old+i, 11);
      mesh_pyramid.BoundaryRef(4*i+1).InitTriangular(n2, n3, nb_old+i, 11);
      mesh_pyramid.BoundaryRef(4*i+2).InitTriangular(n3, n4, nb_old+i, 11);
      mesh_pyramid.BoundaryRef(4*i+3).InitTriangular(n1, n4, nb_old+i, 11);
      
      mesh_tri.BoundaryRef(4*i).InitTriangular(n1, n2, nb_old+i, 1);
      mesh_tri.BoundaryRef(4*i+1).InitTriangular(n2, n3, nb_old+i, 1);
      mesh_tri.BoundaryRef(4*i+2).InitTriangular(n3, n4, nb_old+i, 1);
      mesh_tri.BoundaryRef(4*i+3).InitTriangular(n1, n4, nb_old+i, 1);
    }
  
  mesh_tri.SortBoundariesRef();
  mesh_tri.FindConnectivity();
  mesh_pyramid.FindConnectivity();
}

int main(int argc, char **argv)
{
  InitMontjoie(argc, argv);
  
  R3::threshold = 1e-6;
  
  if (argc < 10)
    {
      cout<<"Cette commande demande au moins six arguments"<<endl;
      cout<<"cell_wire N R alpha delta1 delta2 delta3 toto.mesh Nc Nt"<<endl;
      cout<<"the file toto.mesh is created"<<endl;
      return -1;
    }
  
  double alpha = pi_wp/180.0*atof(argv[1]);
  double delta1 = atof(argv[2]);
  double delta2 = atof(argv[3]);
  double delta3 = atof(argv[4]);
  string file_output(argv[5]);
  int order = 4;

  double xmax = atof(argv[6]);
  double zmax = 1.0/tan(alpha), ymax = 1.0;
  int nb_subdiv_circle = atoi(argv[7]);
  int nb_subdiv_teta = atoi(argv[8]);
  int nb_points_square = atoi(argv[9]);
  
  Mesh<Dimension3> mesh, mesh2;
  mesh.SetGeometryOrder(order);
  mesh2.SetGeometryOrder(order);
  
  CreateHexahedralMeshHelice(alpha, delta1, delta2, delta3,
			     nb_subdiv_circle, nb_subdiv_teta, ymax, zmax, mesh);
  
  //  mesh.Write("helice.msh"); exit(0);
  
  Mesh<Dimension3> mesh_z0, mesh_teta1;
  
  double h = 0.5;
  int nb_points = 4;
  int nb_points_z = 5;
  CreateMeshAroundCircleZ0(mesh, h, nb_points, xmax, ymax, zmax, nb_subdiv_circle, nb_subdiv_teta, mesh_z0);
  
  CreateMeshAroundCircleTeta1(mesh, h, nb_points, nb_points_z, xmax, ymax, zmax, nb_subdiv_circle, nb_subdiv_teta, mesh_teta1);
       
  // on cree maillage du second fil par symetrie
  mesh2 = mesh;
  R3 ptA(0, 0, 0);
  R3 ex(1, 0, 0), ey(0, 1, 0), ez(0, 0, 1);
  RotateMesh(mesh2, ptA, ez, pi_wp);
  RotateMesh(mesh2, ptA, ey, pi_wp);
  TranslateMesh(mesh2, 2.0*zmax*ez);
  mesh.AppendMesh(mesh2);
  
  Mesh<Dimension3> support_mesh;
    
  // on rajoute les maillages z = 0 et teta = 1
  for (int i = 0; i < mesh_z0.GetNbBoundaryRef(); i++)
    mesh_z0.BoundaryRef(i).SetReference(5);
  
  for (int i = 0; i < mesh_teta1.GetNbBoundaryRef(); i++)
    mesh_teta1.BoundaryRef(i).SetReference(7);
 
  support_mesh.AppendMesh(mesh_z0);
  support_mesh.AppendMesh(mesh_teta1);

  mesh.AppendMesh(mesh_z0);
  mesh.AppendMesh(mesh_teta1);
  
  // on translate les maillages z = 0 et teta = 1 pour forcer la periodicite
  TranslateMesh(mesh_z0, 2.0*zmax*ez);
  TranslateMesh(mesh_teta1, -2.0*ymax*ey);
  
  for (int i = 0; i < mesh_z0.GetNbBoundaryRef(); i++)
    mesh_z0.BoundaryRef(i).SetReference(6);

  for (int i = 0; i < mesh_teta1.GetNbBoundaryRef(); i++)
    mesh_teta1.BoundaryRef(i).SetReference(8);
  
  support_mesh.AppendMesh(mesh_z0);
  support_mesh.AppendMesh(mesh_teta1);

  mesh.AppendMesh(mesh_z0);
  mesh.AppendMesh(mesh_teta1);
  
  // on recupere les faces de reference 4 pour construire des pyramides dessus
  Mesh<Dimension3> mesh_pyramid;
  ConstructPyramidalTransition(mesh, 4, 3, 0.25, support_mesh, 11, mesh_pyramid);
  
  mesh.AppendMesh(mesh_pyramid);
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      bool face_on_z0 = true;
      for (int k = 0; k < mesh.BoundaryRef(i).GetNbVertices(); k++)
	if (abs(mesh.Vertex(mesh.BoundaryRef(i).numVertex(k))(2)) > 1e-10)
	  face_on_z0 = false;
      
      if (face_on_z0)
	mesh.BoundaryRef(i).SetReference(5);
    }
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      bool face_on_teta1 = true;
      for (int k = 0; k < mesh.BoundaryRef(i).GetNbVertices(); k++)
	if (abs(mesh.Vertex(mesh.BoundaryRef(i).numVertex(k))(1) - 1.0) > 1e-10)
	  face_on_teta1 = false;
      
      if (face_on_teta1)
	mesh.BoundaryRef(i).SetReference(7);
    }
    
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      bool face_on_z0 = true;
      for (int k = 0; k < mesh.BoundaryRef(i).GetNbVertices(); k++)
	if (abs(mesh.Vertex(mesh.BoundaryRef(i).numVertex(k))(2) - 2.0*zmax) > 1e-10)
	  face_on_z0 = false;
      
      if (face_on_z0)
	mesh.BoundaryRef(i).SetReference(6);
    }
   
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      bool face_on_teta1 = true;
      for (int k = 0; k < mesh.BoundaryRef(i).GetNbVertices(); k++)
	if (abs(mesh.Vertex(mesh.BoundaryRef(i).numVertex(k))(1) + 1.0) > 1e-10)
	  face_on_teta1 = false;
      
      if (face_on_teta1)
	mesh.BoundaryRef(i).SetReference(8);
    }
  
  Mesh<Dimension3> mesh_x0;  
  CreateSquareMesh(nb_points, nb_points_z, R3(-xmax, -ymax, 0), R3(-xmax, ymax, 0),
		   R3(-xmax, ymax, 2.0*zmax), R3(-xmax, -ymax, 2.0*zmax), mesh_x0);
  
  for (int i = 0; i < mesh_x0.GetNbBoundaryRef(); i++)
    mesh_x0.BoundaryRef(i).SetReference(9);

  Mesh<Dimension3> mesh_square;  
  mesh_square = mesh_x0;
  mesh.AppendMesh(mesh_x0);

  TranslateMesh(mesh_x0, 2.0*xmax*ex);
  
  for (int i = 0; i < mesh_x0.GetNbBoundaryRef(); i++)
    mesh_x0.BoundaryRef(i).SetReference(10);
 
  mesh.AppendMesh(mesh_x0);  
  mesh_square.AppendMesh(mesh_x0);

  Mesh<Dimension3> mesh_tri;
  CreatePyramidLayer(mesh_square, R3(0, 0, 0), mesh_tri, mesh_pyramid, 3);
  
  mesh.AppendMesh(mesh_pyramid);
  
  int nb_faces_tri = 0, nb_vertices_tri = 0;
  Vector<int> IndexVertex(mesh.GetNbVertices());
  IndexVertex.Fill(-1);
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    if ((mesh.BoundaryRef(i).GetReference() == 11) || (mesh.BoundaryRef(i).GetNbElements() == 0))
      {
	for (int k = 0; k < mesh.BoundaryRef(i).GetNbVertices(); k++)
	  {
	    int nv = mesh.BoundaryRef(i).numVertex(k);
	    if (IndexVertex(nv) == -1)
	      {
		IndexVertex(nv) = nb_vertices_tri;
		nb_vertices_tri++;
	      }
	  }
	
	nb_faces_tri++;
      }
  
  mesh_tri.Clear();
  mesh_tri.ReallocateVertices(nb_vertices_tri);
  mesh_tri.ReallocateBoundariesRef(nb_faces_tri);
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    if (IndexVertex(i) >= 0)
      mesh_tri.Vertex(IndexVertex(i)) = mesh.Vertex(i);
  
  nb_faces_tri = 0;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    if ((mesh.BoundaryRef(i).GetReference() == 11) || (mesh.BoundaryRef(i).GetNbElements() == 0))
      {
	int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
	IVect num(nb_vert);
	for (int k = 0; k < nb_vert; k++)
	  {
	    int nv = mesh.BoundaryRef(i).numVertex(k);
	    num(k) = IndexVertex(nv);
	  }
	
	mesh_tri.BoundaryRef(nb_faces_tri).Init(num, 1);
	
	nb_faces_tri++;
      }  
  
  Mesh<Dimension3> mesh_tet;
  CreateTetrahedralMesh(mesh_tri, mesh_tet);

  mesh.AppendMesh(mesh_tet);
  
  // on change quelques references
  mesh.RemoveReference(11);
  for (int i = 0; i < mesh.GetNbElt(); i++)
    if (mesh.Element(i).GetReference() == 0)
      mesh.Element(i).SetReference(3);
  
  // maillage d'une face x = 2*xmax avec des carres  
  CreateSquareMesh(nb_points_square, nb_points_square, R3(2.0*xmax, -ymax, 0), R3(2.0*xmax, ymax, 0),
		   R3(2.0*xmax, ymax, 2.0*zmax), R3(2.0*xmax, -ymax, 2.0*zmax), mesh_square);
  
  mesh_square.ResizeNbReferences(15);
  for (int i = 0; i < mesh_square.GetNbBoundaryRef(); i++)
    mesh_square.BoundaryRef(i).SetReference(12);
  
  CreateTriangularMesh(nb_points, nb_points_square, nb_points_square, nb_points_square,
		       R3(xmax, -ymax, 0), R3(xmax, ymax, 0), R3(2*xmax, ymax, 0.0),
		       R3(2*xmax, -ymax, 0.0), mesh_z0);

  CreateTriangularMesh(nb_points_square, nb_points_square, nb_points_square, nb_points_z,
		       R3(xmax, ymax, 0), R3(2.0*xmax, ymax, 0), R3(2*xmax, ymax, 2.0*zmax),
		       R3(xmax, ymax, 2.0*zmax), mesh_teta1);
  
  mesh_square.AppendMesh(mesh_x0);
  CreatePyramidLayer(mesh_square, R3(1.5*xmax, 0, zmax), mesh_tri, mesh_pyramid, 1);

  mesh_tri.AppendMesh(mesh_z0);
  mesh_tri.AppendMesh(mesh_teta1);
  for (int i = 0; i < mesh_z0.GetNbBoundaryRef(); i++)
    mesh_z0.BoundaryRef(i).SetReference(5);
  
  for (int i = 0; i < mesh_teta1.GetNbBoundaryRef(); i++)
    mesh_teta1.BoundaryRef(i).SetReference(7);
  
  mesh.AppendMesh(mesh_z0);
  mesh.AppendMesh(mesh_teta1);
  
  TranslateMesh(mesh_z0, 2.0*zmax*ez);
  mesh_tri.AppendMesh(mesh_z0);
  
  TranslateMesh(mesh_teta1, -2.0*ymax*ey);
  mesh_tri.AppendMesh(mesh_teta1);
  
  for (int i = 0; i < mesh_z0.GetNbBoundaryRef(); i++)
    mesh_z0.BoundaryRef(i).SetReference(6);
  
  for (int i = 0; i < mesh_teta1.GetNbBoundaryRef(); i++)
    mesh_teta1.BoundaryRef(i).SetReference(8);
  
  mesh.AppendMesh(mesh_z0);
  mesh.AppendMesh(mesh_teta1);
  
  CreateTetrahedralMesh(mesh_tri, mesh_tet);
  for (int i = 0; i < mesh_tet.GetNbElt(); i++)
    mesh_tet.Element(i).SetReference(1);

  for (int i = 0; i < mesh_square.GetNbBoundaryRef(); i++)
    mesh_square.BoundaryRef(i).SetReference(12);
  
  mesh.AppendMesh(mesh_tet);
  mesh.AppendMesh(mesh_pyramid);
  mesh.AppendMesh(mesh_square);
  
  // on ecrit le maillage final
  mesh.Write(file_output);

  /* VectReal_wp subdiv(5);
  subdiv.Fill(); Mlt(0.25, subdiv); 
  mesh.SubdivideMesh(subdiv);
  mesh.Write("subdiv.mesh"); */
  
  // on enleve les fichiers temporaires
  std::remove("tmp_geo.geo"); std::remove("sortie");
  
  return FinalizeMontjoie();
}  

