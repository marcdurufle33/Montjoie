#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

#include "Output/ParameterOutput.hxx"
#include "Output/GridInterpolation.hxx"
#include "Output/ParameterOutputInline.cxx"
#include "Output/ParameterOutput.cxx"
#include "Output/GridInterpolationInline.cxx"
#include "Output/GridInterpolation.cxx"

using namespace Montjoie;


int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc < 2)
    {
      cout<<"This command requires at least one argument"<<endl;
      cout<<"geo3pack carre.geo"<<endl;
      cout<<"the file carre.mesh is created"<<endl;
      
      cout<<"You can also specify the type of mesh and algorithm to use "<<endl;
      cout<<"geo2pack carre.geo quad advance"<<endl;
      cout<<"this command will generate a quadrilateral mesh with front advance algorithm"<<endl<<endl;
      cout<<"geo2pack carre.geo hybrid poly"<<endl;
      cout<<"geo2pack carre.geo tri advance"<<endl;
      cout<<"are other possibilities"<<endl;
    }
  
  string file_name(argv[1]);
  if (file_name.find(".geo") != string::npos)
    file_name = file_name.substr(0, file_name.size()-4);
  
  string file_geo = file_name + ".geo";
  string file_mesh = file_name + ".mesh";
  
  string tmp_geo = "temp.geo";
  string tmp_msh = "temp.msh";
  
  ifstream file_in(file_geo.data());
  if (!file_in)
    {
      cout<<"File "<<file_geo<<" does not exist "<<endl;
      return -1;
    }
  
  // counting number of lines
  char line[1024]; Vector<string> lignes; int nb_lines = 0;
  while (!file_in.eof())
    {
      file_in.getline(line, 1024);
      string ligne(line);
      StringTrim(ligne);
      if (ligne.size() > 2)
	if (ligne.compare(0, 2, "//") != 0)
	  nb_lines++;
    }
  
  // storing lines of geo files
  lignes.Reallocate(nb_lines); nb_lines = 0;
  file_in.close(); file_in.open(file_geo.data());
  while (!file_in.eof())
    {
      file_in.getline(line, 1024);
      // spaces on the extremities are removed
      string ligne(line);
      DeleteSpaceAtExtremityOfString(ligne);
      if (ligne.size() > 2)
	if (ligne.compare(0, 2, "//") != 0)
	  {
	    // we don't store commented lines
	    lignes(nb_lines) = ligne;
	    nb_lines++;
	  }
    }
  
  file_in.close();
  
  // now, we look for surfaces and volumes
  Vector<IVect> surface_loops; IVect num_surface_loops;  
  Vector<IVect> volumes; IVect num_volumes;  
  int ref_max_loop = 0;
  IVect VolumeReferenced;
  for (int i = 0; i < nb_lines; i++)
    {
      string::size_type index = lignes(i).find("Surface Loop");
      if (index != string::npos)
	{
	  // we get the number between parenthesis
	  string::size_type index1 = lignes(i).find('(',index);
	  string::size_type index2 = lignes(i).find(')',index);
	  if ((index1 != string::npos)&&(index2 > (index1+1)))
	    {
	      int iloop = to_num<int>(lignes(i).substr(index1+1, index2-index1-1));	      
	      num_surface_loops.PushBack(iloop); int j = num_surface_loops.GetM()-1;
	      if (surface_loops.GetM() <= j)
		surface_loops.Resize(j+1);
	      
	      // we get numbers of all surfaces defining the loop
	      index1 = lignes(i).find('{',index);
	      index2 = lignes(i).find('}',index);
	      if ((index1 != string::npos)&&(index2 > (index1+1)))
		{
		  Vector<string> param;
		  StringTokenize(lignes(i).substr(index1+1,index2-index1-1), param, string(","));
		  surface_loops(j).Reallocate(param.GetM());		  
		  for (int k = 0; k < param.GetM(); k++)
		    {
		      surface_loops(j)(k) = abs(to_num<int>(param(k)));
		      ref_max_loop = max(ref_max_loop, surface_loops(j)(k));
		    }
		}
	      else
		{
		  cout<<"Is it a valid geo file ?"<<endl;
		  abort();
		}
	    }
	  else
	    {
	      cout<<"Is it a valid geo file ?"<<endl;
	      abort();
	    }
	}      
         
      index = lignes(i).find("Physical Volume");
      if (index == 0)
	{
	  // we get the number between parenthesis
	  string::size_type index1 = lignes(i).find('(',index);
	  string::size_type index2 = lignes(i).find(')',index);
	  if ((index1 != string::npos)&&(index2 > (index1+1)))
	    {
	      int iphys = to_num<int>(lignes(i).substr(index1+1, index2-index1-1));
	      
	      // we get numbers of all volumes defining the physical volume
	      index1 = lignes(i).find('{',index);
	      index2 = lignes(i).find('}',index);
	      if ((index1 != string::npos)&&(index2 > (index1+1)))
		{
		  Vector<string> param;
		  StringTokenize(lignes(i).substr(index1+1,index2-index1-1), param, string(","));
		  for (int k = 0; k < param.GetM(); k++)
		    {
		      int ivol = to_num<int>(param(k));		  
		      if (ivol >= VolumeReferenced.GetM())
			{
			  int nbef = VolumeReferenced.GetM();
			  VolumeReferenced.Resize(ivol+1);
			  for (int k = nbef; k <= ivol; k++)
			    VolumeReferenced(k) = 0;
			}
		      
		      VolumeReferenced(ivol) = iphys;
		    }
		}
	      else
		{
		  cout<<"Is it a valid geo file ?"<<endl;
		  abort();
		}
	    }
	  else
	    {
	      cout<<"Is it a valid geo file ?"<<endl;
	      abort();
	    }
	}
      
      index = lignes(i).find("Volume");
      if (index == 0)
	{
	  // we get the number between parenthesis
	  string::size_type index1 = lignes(i).find('(',index);
	  string::size_type index2 = lignes(i).find(')',index);
	  if ((index1 != string::npos)&&(index2 > (index1+1)))
	    {
	      int ivol = to_num<int>(lignes(i).substr(index1+1, index2-index1-1));
	      num_volumes.PushBack(ivol); int j = num_volumes.GetM()-1;
	      if (volumes.GetM() <= j)
		volumes.Resize(j+1);
	      
	      // we get numbers of all loops defining the surface
	      index1 = lignes(i).find('{',index);
	      index2 = lignes(i).find('}',index);
	      if ((index1 != string::npos)&&(index2 > (index1+1)))
		{
		  Vector<string> param;
		  StringTokenize(lignes(i).substr(index1+1,index2-index1-1), param, string(","));
		  volumes(j).Reallocate(param.GetM());		  
		  for (int k = 0; k < param.GetM(); k++)
		    volumes(j)(k) = to_num<int>(param(k));		  
		}
	      else
		{
		  cout<<"Is it a valid geo file ?"<<endl;
		  abort();
		}
	    }
	  else
	    {
	      cout<<"Is it a valid geo file ?"<<endl;
	      abort();
	    }
	}
    }
  DISP(VolumeReferenced);
  
  // we add a physical surface for boundary faces not referenced in the initial geo file
  int ref_unused = 0; // reference unused in Physical Surfaces of the geo file
  // SurfaceReferenced(i) = -1 -> surface not present in Surface Loop
  // SurfaceReferenced(i) = 0 -> surface present in Surface Loop, but not in physical surface
  // SurfaceReferenced(i) > 0 -> surface present in Surface Loop and in physical surface
  Vector<int> SurfaceReferenced(ref_max_loop+1); SurfaceReferenced.Fill(-1);
  for (int i = 0; i < surface_loops.GetM(); i++)
    for (int j = 0; j < surface_loops(i).GetM(); j++)
      SurfaceReferenced(surface_loops(i)(j)) = 0;
  
  for (int i = 0; i < nb_lines; i++)
    {
      string::size_type index = lignes(i).find("Physical Surface");
      if (index != string::npos)
	{
	  // we get the number between parenthesis
	  string::size_type index1 = lignes(i).find('(',index);
	  string::size_type index2 = lignes(i).find(')',index);
	  if ((index1 != string::npos)&&(index2 > (index1+1)))
	    {
	      int iref = to_num<int>(lignes(i).substr(index1+1, index2-index1-1));
	      ref_unused = max(ref_unused, iref);
	      
	      // we get numbers of all lines defining the loop
	      index1 = lignes(i).find('{',index);
	      index2 = lignes(i).find('}',index);
	      if ((index1 != string::npos)&&(index2 > (index1+1)))
		{
		  Vector<string> param;
		  StringTokenize(lignes(i).substr(index1+1,index2-index1-1), param, string(","));
		  for (int k = 0; k < param.GetM(); k++)
		    {
		      int num_line = to_num<int>(param(k));
		      if (abs(num_line) <= ref_max_loop)
			{
			  SurfaceReferenced(abs(num_line)) = iref;
			}
		    }		      
		}
	      else
		{
		  cout<<"Is it a valid geo file ?"<<endl;
		  abort();
		}
	    }
	  else
	    {
	      cout<<"Is it a valid geo file ?"<<endl;
	      abort();
	    }
	}
    }
  
  // we regroup all the non referenced surfaces
  IVect surface_to_add_physical(ref_max_loop); int nb = 0;
  ref_unused++;
  for (int i = 1; i < SurfaceReferenced.GetM(); i++)
    if (SurfaceReferenced(i) == 0)
      surface_to_add_physical(nb++) = i;
  
  if (nb == 0)
    surface_to_add_physical.Clear();
  else
    surface_to_add_physical.Resize(nb);  
  
  // generation of mesh with gmsh
  ofstream file_out(tmp_geo.data());
  DISP(lignes);
  for (int i = 0; i < nb_lines; i++)
    if (lignes(i).find("Physical Volume") == string::npos)
      file_out << lignes(i) << '\n';
  
  // adding a Physical Surface for unused surface
  if (surface_to_add_physical.GetM() > 0)
    {
      file_out << "Physical Surface("<<ref_unused<<") = {"<<surface_to_add_physical(0);
      for (int j = 1; j < surface_to_add_physical.GetM(); j++)
	file_out << ", " << surface_to_add_physical(j);
      
      file_out << "};\n";
    }
  
  // adding physical volume (for each different volume)
  for (int i = 0; i < num_volumes.GetM(); i++)
    file_out<<"Physical Volume("<<i+1<<") = {"<<num_volumes(i)<<"};\n";
  
  file_out.close();
  
  // calling gmsh
  string command_line = string("gmsh -o ") + tmp_msh + " " + tmp_geo + " -3";
  int ierr = system(command_line.data());
  if (ierr != 0)
    {
      cout<<"Error during the generation of msh file"<<endl;
      abort();
    }
  
  // we read the generated mesh
  Mesh<Dimension3> mesh;
  mesh.Read(tmp_msh);
  mesh.SetGeometryOrder(1);

  // we force coherence mesh to remove unused vertices
  mesh.ForceCoherenceMesh();

  // we change references of boundaries
  int geom_max = 0;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetGeometryReference();
      geom_max = max(geom_max, ref);
      mesh.BoundaryRef(i).SetReference(ref);
    }
  
  // finding all the edges and faces of the mesh
  mesh.FindConnectivity();
  
  // counting edges on the interface between two different geometric surfaces
  /* mesh.ReallocateEdgesRef(mesh.GetNbBoundaryRef()*4);
  int nb_edge_ref = 0;
  for (int i = 0; i < mesh.GetNbEdges(); i++)
    {
      int ref1 = 0, ref_current;
      bool interface_edge = false;
      for (int j = 0; j < mesh.GetEdge(i).GetNbFaces(); j++)
	{
	  int nf = mesh.GetEdge(i).numFace(j);
	  ref_current = mesh.Boundary(nf).GetGeometryReference();
	  if (ref_current > 0)
	    {
	      if (ref1 == 0)
		ref1 = ref_current;
	      else if (ref1 != ref_current)
		interface_edge = true;
	    }
	}
      
      if (interface_edge)
	{
	  mesh.GetEdge(i).SetReference(1);
	  mesh.EdgeRef(nb_edge_ref) = mesh.GetEdge(i);
	  nb_edge_ref++;
	}
    }
    mesh.ResizeEdgesRef(nb_edge_ref);*/
  
  // loops on all the geometrical surfaces
  Vector<Mesh<Dimension3> > surface_mesh(geom_max+1);
  for (int i = 1; i <= geom_max; i++)
    {
      // we write mhs file
      IVect IndexVert(mesh.GetNbVertices()); IndexVert.Fill(-1);
      int nb_vertices = 0, nb_faces = 0; IVect ListeVertex(mesh.GetNbVertices());
      for (int nf = 0; nf < mesh.GetNbBoundaryRef(); nf++)
	if (mesh.BoundaryRef(nf).GetReference() == i)
	  {
	    for (int j = 0; j < mesh.BoundaryRef(nf).GetNbVertices(); j++)
	      {
		int nv = mesh.BoundaryRef(nf).numVertex(j);
		if (IndexVert(nv) == -1)
		  {
		    IndexVert(nv) = nb_vertices;
		    ListeVertex(nb_vertices) = nv;
		    nb_vertices++;
		  }
	      }
	    nb_faces++;
	  }
      
      if (nb_faces > 0)
	{
	  
	  ofstream surface_file("old_surface.mhs");
	  surface_file.precision(15);
	  surface_file << nb_vertices << '\n';
	  for (int nv = 0; nv < nb_vertices; nv++)
	    surface_file << mesh.Vertex(ListeVertex(nv)) << " 1 \n";
	  
	  surface_file << "\n 0 \n \n";
	  
	  surface_file << "103 "<<nb_faces<<'\n'; 
	  for (int nf = 0; nf < mesh.GetNbBoundaryRef(); nf++)
	    if (mesh.BoundaryRef(nf).GetReference() == i)
	      {
		for (int j = 0; j < mesh.BoundaryRef(nf).GetNbVertices(); j++)
		  {
		    int nv = mesh.BoundaryRef(nf).numVertex(j);
		    surface_file << IndexVert(nv) + 1 <<" ";
		  }
	      }
	  
	  for (int nf = 0; nf < mesh.GetNbBoundaryRef(); nf++)
	    if (mesh.BoundaryRef(nf).GetReference() == i)
	      {
		surface_file << " 5 ";
		for (int j = 0; j < mesh.BoundaryRef(nf).GetNbVertices(); j++)
		  {
		    // int ne = mesh.BoundaryRef(nf).numEdge(j);
		    surface_file <<  "5 ";
		  }
		surface_file << '\n';
	      }
	  
	  surface_file.close();
	  
	  ofstream curve_file("curve_input");
	  curve_file<<"0 0"<<endl;
	  curve_file.close();
	  
	  ofstream main_file("main_input");
	  int cs_in_type(1), mesh_in_type(1), mesh_out_type(1);
	  
	  main_file<<"251 0 "<<cs_in_type<<" "<<mesh_in_type<<" 0 0 "<<mesh_out_type<<'\n';
	  main_file<<"curve_input"<<'\n';
	  main_file<<"old_surface.mhs"<<'\n';
	  main_file<<"new_surface.mhs"<<'\n';
	  main_file<<"1e-10 0 0"<<'\n';
	  int implmidn(0), nodelem(-4), mtype(3), nimpiter(0), nsmpas(0);
	  Real_wp quadmu(0.01), optmu(0.3);
	  main_file<<implmidn<<" "<<nodelem<<" "<<mtype<<" "<<quadmu<<" "<<nimpiter<<" "<<nsmpas<<" "<<optmu<<'\n';
	  main_file.close();
	  
	  // we generate a mesh file with geompack
	  command_line = string("zgp0803_linux main_input test.out");
	  ierr = system(command_line.data());
	  if (ierr != 0)
	    {
	      cout<<"Error during the generation of mh3 file"<<endl;
	      abort();
	    }
	  
	  surface_mesh(i).Read("new_surface.mhs");
	  surface_mesh(i).SplitIntoHexahedra();      
	  
	  // affecting the geometry reference
	  for (int j = 0; j < surface_mesh(i).GetNbBoundaryRef(); j++)
	    surface_mesh(i).BoundaryRef(j).SetReference(i);
	}
    }

  // mesh that will concatenate all the other meshes
  Mesh<Dimension3> all_mesh;
  
  GridInterpolation<Dimension3> grid;
  grid.InitInterpolationGrid(mesh);
  
  // loop on all volumes to generate
  for (int ivol = 0;  ivol < num_volumes.GetM(); ivol++)
    {
      // we group all the surfaces on the boundary of the volume
      Mesh<Dimension3> boundary_mesh;
      int nb_bound_ref = 0;
      for (int j = 0; j < volumes(ivol).GetM(); j++)
	{
	  int isurf = volumes(ivol)(j);
	  int num_surf = -1;
	  for (int iloop = 0; iloop < num_surface_loops.GetM(); iloop++)
	    if (num_surface_loops(iloop) == isurf)
	      num_surf = iloop;
	  
	  for (int is = 0; is < surface_loops(num_surf).GetM(); is++)
	    {
	      int igeom = surface_loops(num_surf)(is);
	      if (SurfaceReferenced(igeom) > 0)
		nb_bound_ref += surface_mesh(igeom).GetNbBoundaryRef();
	      
	      boundary_mesh.AppendMesh(surface_mesh(igeom));
	    }
	}
      
      
      // we check the orientation of faces (counterclockwise when viewing from outside the volume)
      VectR3& GlobalCoord = grid.GetGlobalCoordinate();
      GlobalCoord.Reallocate(2*boundary_mesh.GetNbBoundaryRef());
      for (int i = 0; i < boundary_mesh.GetNbBoundaryRef(); i++)
	{
	  int nf = boundary_mesh.BoundaryRef(i).GetNbVertices()-1;
	  int n0 = boundary_mesh.BoundaryRef(i).numVertex(0);
	  int n1 = boundary_mesh.BoundaryRef(i).numVertex(1);
	  int n2 = boundary_mesh.BoundaryRef(i).numVertex(nf);
	  R3 ptO = boundary_mesh.Vertex(n0), t1, t2, normale, vec_u, center;
	  t1 = boundary_mesh.Vertex(n1) - ptO;
	  t2 = boundary_mesh.Vertex(n2) - ptO;
	  TimesProd(t1, t2, normale);
	  Mlt(0.25/sqrt(Norm2(normale)), normale);
	  
	  Add(boundary_mesh.Vertex(n0), boundary_mesh.Vertex(n1), center);
	  center += boundary_mesh.Vertex(n2); Mlt(1.0/3, center);
	  
	  GlobalCoord(2*i) = center + normale;
	  GlobalCoord(2*i+1) = center - normale;
	}
      
      // we try to localize Vertices in the volume mesh
      grid.LocalizePoints(mesh);
      
      // DISP(i); DISP(ptO); DISP(t1); DISP(t2); DISP(vec_u); DISP(normale);
      // DISP(DotProd(vec_u, normale));
      for (int i = 0; i < boundary_mesh.GetNbBoundaryRef(); i++)
	{
	  // int num_elem = grid.GetElementNumber(2*i);
	  int num_elem2 = grid.GetElementNumber(2*i+1);
	  int ref = boundary_mesh.BoundaryRef(i).GetReference();
	  int nf = boundary_mesh.BoundaryRef(i).GetNbVertices()-1;
	  int n0 = boundary_mesh.BoundaryRef(i).numVertex(0);
	  int n1 = boundary_mesh.BoundaryRef(i).numVertex(1);
	  int n2 = boundary_mesh.BoundaryRef(i).numVertex(nf);	  
	  
	  if ((num_elem2 < 0)||(mesh.Element(num_elem2).GetReference() != ivol+1))
	    {
	      if (nf == 3)
		{
		  int ns = boundary_mesh.BoundaryRef(i).numVertex(2);
		  boundary_mesh.BoundaryRef(i).InitQuadrangular(n0, n2, ns, n1, ref);
		}
	      else
		{
		  boundary_mesh.BoundaryRef(i).InitTriangular(n0, n2, n1, ref);
		}
	    }
	}
      	  
      boundary_mesh.Write("new_surface.mhs");
      
      // generation of geompack files
      ofstream main_file("main_input");
      ofstream curve_file("curve_input");
      curve_file<<"0 0"<<endl;
      curve_file.close();
      
      main_file<<"351 0 1 1 0 0 1"<<endl;
      main_file<<"curve_input"<<endl;
      main_file<<"new_surface.mhs"<<endl;
      main_file<<"mesh_output.mh3"<<endl;
      main_file<<"1e-10 0 0"<<endl;
      
      int implmidn(0), nodelem(-8), mtype(3), nimpiter(0), nsmpas(0), regcod(0);
      Real_wp quadmu(0.1), optmu(0.3), mergmu(0.05);
      main_file << "3 " << implmidn << " " << regcod << " " << nodelem << " " << mtype
		<< " " << quadmu << " " << nimpiter << " " << nsmpas << " " << optmu << " " << mergmu << endl;
      
      main_file.close();
            
      // we generate a mesh file with geompack
      command_line = string("zgp0803_linux main_input test.out");
      ierr = system(command_line.data());
      if (ierr != 0)
	{
	  cout<<"Error during the generation of mh3 file"<<endl;
	  abort();
	}
      
      Mesh<Dimension3> output_mesh;
      output_mesh.Read(string("mesh_output.mh3"));
      
      // we change references
      output_mesh.ReallocateBoundariesRef(nb_bound_ref);
      DISP(nb_bound_ref);
      nb_bound_ref = 0;
      for (int i = 0; i < boundary_mesh.GetNbBoundaryRef(); i++)
	{
	  int ref = SurfaceReferenced(boundary_mesh.BoundaryRef(i).GetReference());
	  if (ref > 0)
	    {
	      output_mesh.BoundaryRef(nb_bound_ref) = boundary_mesh.BoundaryRef(i);
	      output_mesh.BoundaryRef(nb_bound_ref).SetReference(ref);
	      nb_bound_ref++;
	    }
	}
      DISP(nb_bound_ref);
      
      // volume references
      int ref = VolumeReferenced(num_volumes(ivol));
      for (int i = 0; i < output_mesh.GetNbElt(); i++)
	output_mesh.Element(i).SetReference(ref);
      
      // we append this mesh to the global mesh
      all_mesh.AppendMesh(output_mesh, true);
    }
  
  // checking if symmetry is needed
  for (int i = 0; i < nb_lines; i++)
    {
      if (lignes(i).find("SymmetryX") == 0)
	all_mesh.SymmetrizeMesh(0, 0.0, R3(0, 0, 0));
      
      if (lignes(i).find("SymmetryY") == 0)
	all_mesh.SymmetrizeMesh(1, 0.0, R3());

      if (lignes(i).find("SymmetryZ") == 0)
	all_mesh.SymmetrizeMesh(2, 0.0, R3());
    }
  
  // writing final mesh obtained
  all_mesh.Write(file_mesh);
  
  return FinalizeMontjoie();
}  
