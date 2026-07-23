
#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;


void FindAllLoops(Mesh<Dimension2>& mesh, Vector<IVect>& all_loops, IVect& type_loop)
{
  int ref_max = 0;
  Vector<bool> BoundaryUsed(mesh.GetNbBoundary());
  BoundaryUsed.Fill(false);
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    ref_max = max(ref_max, mesh.BoundaryRef(i).GetReference());
  
  IVect list_points(mesh.GetNbVertices()); list_points.Fill(-1);
  R2 tangente, normale;
  IVect ref_cond(mesh.GetNbReferences()+1); ref_cond.Fill();
  for (int ref = 1; ref <= ref_max; ref++)
    {
      bool test_loop = true;
      while (test_loop)
	{
	  int num_init = -1;
	  IVect NbBoundary(mesh.GetNbVertices()); NbBoundary.Fill(0);
	  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	    {
	      int num_edge = i;
	      if ((!BoundaryUsed(num_edge))&&(mesh.BoundaryRef(i).GetReference() == ref))
		{
		  num_init = i;
		  NbBoundary(mesh.BoundaryRef(i).numVertex(0))++;
		  NbBoundary(mesh.BoundaryRef(i).numVertex(1))++;
		}
	    }
	  
	  if (num_init == -1)
	    test_loop = false;
	  else
	    {
	      // checking that the loop is closed
	      int ipoint = -1;
	      for (int k = 0; k < mesh.GetNbVertices(); k++)
		if (NbBoundary(k) == 1)
		  {
		    cout<<"Only closed loops are allowed "<<endl;
		    mesh.Write(string("test.mesh"));
		    abort();
		  }
	      
	      ipoint = mesh.BoundaryRef(num_init).numVertex(0);
	      
	      num_init = ipoint;
	      // looking all other points on the loop
	      int iback = num_init, inext = -1, nb = 0, num_edge = -1;
	      while (inext != num_init)
		{
		  iback = ipoint; ipoint = inext;
		  if (inext == -1)
		    ipoint = num_init;
		  
		  inext = mesh.FindFollowingVertex(iback, ipoint, ref, ref_cond, BoundaryUsed, num_edge);
		  if (nb == 0)
		    {
		      // checking the orientation of the loop
		      int num_elem = mesh.Boundary(num_edge).numElement(0);
		      int iortho = -1;
		      for (int k = 0; k < mesh.Element(num_elem).GetNbVertices(); k++)
			if ((mesh.Element(num_elem).numVertex(k) != ipoint)
			    &&(mesh.Element(num_elem).numVertex(k) != inext))
			  iortho = mesh.Element(num_elem).numVertex(k);
		      
		      tangente = mesh.Vertex(inext) - mesh.Vertex(ipoint);
		      normale = mesh.Vertex(iortho) - mesh.Vertex(ipoint);
		      
		      // if t \vec n < 0 on the outer boundary (or t vec n > 0 on a hole boundary
		      // we have to change the orientation, otherwise it is okay
		      Real_wp prod = tangente(0)*normale(1) - tangente(1)*normale(0);
		      if (((prod < 0)&&(ref==1))
			  ||((prod > 0)&&(ref>1)))
			{
			  iback = inext;
			  inext = ipoint;
			  ipoint = iback;
			  num_init = iback;
			}
		    }
		  
		  if (num_edge != -1)
		    BoundaryUsed(num_edge) = true;
		  
		  
		  if (inext == -1)
		    {
		      // no more edges on the loop
		      inext = num_init;
		    }
		  else
		    {
		      list_points(nb) = ipoint; nb++;
		    }
		}
	      
	      // adding the loop found to the list all_loops
	      IVect num_points(nb);
	      for (int i = 0; i < nb; i++)
		num_points(i) = list_points(i);
	      
	      all_loops.PushBack(num_points);
	      if (ref == 1)
		type_loop.PushBack(0);
	      else
		type_loop.PushBack(1);
	    }
	}
    }
}

int main(int argc, char **argv) 
{
  
  InitMontjoie(argc, argv);
  
  if (argc < 2)
    {
      cout<<"This command requires at least one argument"<<endl;
      cout<<"geo2pack carre.geo"<<endl;
      cout<<"the file carre.mesh is created"<<endl;
      
      cout<<"You can also specify the type of mesh and algorithm to use "<<endl;
      cout<<"geo2pack carre.geo quad advance"<<endl;
      cout<<"this command will generate a quadrilateral mesh with front advance algorithm"<<endl<<endl;
      cout<<"geo2pack carre.geo hybrid poly"<<endl;
      cout<<"geo2pack carre.geo tri advance"<<endl;
      cout<<"are other possibilities"<<endl;
    }
  
  int type_operation = 207; // (202 : polygon, 203 : Delaunay, 207 : advancing front)
  int nodelem = 4; // (-4 : hybrid, 3 : triangle, 4 : quadrilateral
  if (argc == 4)
    {
      string type_mesh(argv[2]);
      string type_algo(argv[3]);
      if (!type_mesh.compare("tri"))
	nodelem = 3;
      else if (!type_mesh.compare("quad"))
	nodelem = 4;
      else
	nodelem = -4;
      
      if (!type_algo.compare("poly"))
	type_operation = 202;
      else
	type_operation = 207;
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
  
  // now, we look for all loops and surfaces
  IVect num_loops; int iloop = -1, isurf = -1;
  Vector<IVect> line_loops;
  Vector<IVect> surface_loops; IVect num_surfaces;  
  int ref_max_loop = 0; int ref_max_surf = 0;
  for (int i = 0; i < nb_lines; i++)
    {
      string::size_type index = lignes(i).find("Line Loop");
      if (index != string::npos)
	{
	  // we get the number between parenthesis
	  string::size_type index1 = lignes(i).find('(',index);
	  string::size_type index2 = lignes(i).find(')',index);
	  if ((index1 != string::npos)&&(index2 > (index1+1)))
	    {
	      iloop = to_num<int>(lignes(i).substr(index1+1, index2-index1-1));	      
	      num_loops.PushBack(iloop); int j = num_loops.GetM()-1;
	      if (line_loops.GetM() <= j)
		line_loops.Resize(j+1);
	      
	      // we get numbers of all lines defining the loop
	      index1 = lignes(i).find('{',index);
	      index2 = lignes(i).find('}',index);
	      if ((index1 != string::npos)&&(index2 > (index1+1)))
		{
		  Vector<string> param;
		  StringTokenize(lignes(i).substr(index1+1,index2-index1-1), param, string(","));
		  line_loops(j).Reallocate(param.GetM());		  
		  for (int k = 0; k < param.GetM(); k++)
		    {
		      line_loops(j)(k) = abs(to_num<int>(param(k)));
		      ref_max_loop = max(ref_max_loop, line_loops(j)(k));
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
            
      index = lignes(i).find("Plane Surface");
      if (index != string::npos)
	{
	  // we get the number between parenthesis
	  string::size_type index1 = lignes(i).find('(',index);
	  string::size_type index2 = lignes(i).find(')',index);
	  if ((index1 != string::npos)&&(index2 > (index1+1)))
	    {
	      isurf = to_num<int>(lignes(i).substr(index1+1, index2-index1-1));
	      ref_max_surf = max(ref_max_surf, isurf);
	      num_surfaces.PushBack(isurf); int j = num_surfaces.GetM()-1;
	      if (surface_loops.GetM() <= j)
		surface_loops.Resize(j+1);
	      
	      // we get numbers of all loops defining the surface
	      index1 = lignes(i).find('{',index);
	      index2 = lignes(i).find('}',index);
	      if ((index1 != string::npos)&&(index2 > (index1+1)))
		{
		  Vector<string> param;
		  StringTokenize(lignes(i).substr(index1+1,index2-index1-1), param, string(","));
		  surface_loops(j).Reallocate(param.GetM());		  
		  for (int k = 0; k < param.GetM(); k++)
		    surface_loops(j)(k) = to_num<int>(param(k));
		  
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
  
  // we add a physical line for boundary edges not referenced in the initial geo file
  int ref_unused = 0; // reference unused in Physical Lines of the geo file
  // LineReferenced(i) = -1 -> line not present in Line Loop
  // LineReferenced(i) = 0 -> line present in Line Loop, but not in physical line
  // LineReferenced(i) = 1 -> line present in Line Loop and in physical line
  Vector<int> LineReferenced(ref_max_loop+1); LineReferenced.Fill(-1);
  for (int i = 0; i < line_loops.GetM(); i++)
    for (int j = 0; j < line_loops(i).GetM(); j++)
      LineReferenced(line_loops(i)(j)) = 0;
  
  IVect ref_surfaces(ref_max_surf+1); ref_surfaces.Fill(0);
  for (int i = 0; i < nb_lines; i++)
    {
      string::size_type index = lignes(i).find("Physical Line");
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
			  LineReferenced(abs(num_line)) = 1;
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
      
      // we modify ref_surfaces
      index = lignes(i).find("Physical Surface");
      if (index != string::npos)
	{
	  // we get the number between parenthesis
	  string::size_type index1 = lignes(i).find('(',index);
	  string::size_type index2 = lignes(i).find(')',index);
	  if ((index1 != string::npos)&&(index2 > (index1+1)))
	    {
	      int ref = to_num<int>(lignes(i).substr(index1+1, index2-index1-1));
	      	      
	      // we get numbers of all lines defining the loop
	      index1 = lignes(i).find('{',index);
	      index2 = lignes(i).find('}',index);
	      if ((index1 != string::npos)&&(index2 > (index1+1)))
		{
		  Vector<string> param;
		  StringTokenize(lignes(i).substr(index1+1,index2-index1-1), param, string(","));
		  for (int k = 0; k < param.GetM(); k++)
		    {
		      int num_surf = to_num<int>(param(k));
		      if (num_surf <= ref_max_surf)
			{
			  ref_surfaces(num_surf) = ref;
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
  
  // we regroup all the non referenced lines
  IVect line_to_add_physical_line(ref_max_loop); int nb = 0;
  ref_unused++;
  for (int i = 1; i < LineReferenced.GetM(); i++)
    if (LineReferenced(i) == 0)
      line_to_add_physical_line(nb++) = i;
  
  if (nb == 0)
    line_to_add_physical_line.Clear();
  else
    line_to_add_physical_line.Resize(nb);  
  
  // mesh that will concatenate all the other meshes
  Mesh<Dimension2> all_mesh;
  
  // loop on all surfaces to generate
  for (int isurf = 0;  isurf < num_surfaces.GetM(); isurf++)
    {
      // we create a geo file with only the generation of the surface
      ofstream file_out(tmp_geo.data());
      for (int i = 0; i < nb_lines; i++)
	{
	  string::size_type index = lignes(i).find("Plane Surface");
	  string::size_type index2 = lignes(i).find("Physical Surface");
	  if ((index == string::npos)&&(index2 == string::npos))
	    {
	      // this line is not a line with Plane Surface or Physical Surface
	      file_out << lignes(i) << '\n';
	    }
	}
      
      // we add a physical line for unreferenced lines
      if (line_to_add_physical_line.GetM() > 0)
	{
	  file_out << "Physical Line("<<ref_unused<<") = {"<<line_to_add_physical_line(0);
	  for (int i = 1; i < line_to_add_physical_line.GetM(); i++)
	    file_out << ", " << line_to_add_physical_line(i);
	  
	  file_out<<"};\n";
	}
      
      // now we add the treated Surface
      file_out << "Plane Surface(1) = {";
      for (int j = 0; j < surface_loops(isurf).GetM(); j++)
	{
	  file_out<<surface_loops(isurf)(j);
	  if (j == surface_loops(isurf).GetM()-1)
	    file_out << "};"<<'\n';
	  else
	    file_out<<", ";
	}
      
      file_out << "Physical Surface(1) = {1}; "<<endl;      
      file_out.close();
      
      // we generate the mesh file with gmsh
      string command_line = string("gmsh -o ") + tmp_msh + " " + tmp_geo + " -2";
      int ierr = system(command_line.data());
      if (ierr != 0)
	{
	  cout<<"Error during the generation of msh file"<<endl;
	  abort();
	}
      
      // we read the generated mesh
      Mesh<Dimension2> mesh;
      mesh.Read(tmp_msh);
      
      // we force coherence mesh to remove unused vertices
      mesh.ForceCoherenceMesh();
      
      // we change references of boundaries
      IVect OldReferences(mesh.GetNbBoundaryRef());
      int geom_max = 0;
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	{
	  OldReferences(i) = mesh.BoundaryRef(i).GetReference();
	  geom_max = max(geom_max, mesh.BoundaryRef(i).GetGeometryReference());
	}
      
      // outer boundary is referenced 1, and hole boundaries referenced 2, 3, etc
      IVect new_ref_geom(geom_max+1); new_ref_geom.Fill(0);
      for (int i = 0; i < surface_loops(isurf).GetM(); i++)
	{
	  int iloop = -1;
	  for (int k = 0; k < num_loops.GetM(); k++)
	    if (num_loops(k) == surface_loops(isurf)(i))
	      iloop = k;
	  
	  for (int j = 0; j < line_loops(iloop).GetM(); j++)
	    new_ref_geom(abs(line_loops(iloop)(j))) = i+1;
	}
      
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	{
	  int ref = new_ref_geom(mesh.BoundaryRef(i).GetGeometryReference());
	  mesh.BoundaryRef(i).SetReference(ref);
	}
      
      // we find all the edges of the mesh, and link them with boundary edges
      mesh.FindConnectivity();
      
      // generation of geompack files
      ofstream main_file("main_input");
      ofstream region_file("region_input");
      ofstream curve_file("curve_input");
      
      main_file<<type_operation<<" 1 1 0 0 0 1"<<endl;
      main_file<<"region_input"<<endl;
      main_file<<"curve_input"<<endl;
      main_file<<"mesh_output.mh2"<<endl;
      main_file<<"1e-10 0 0"<<endl;
      Real_wp munif = 1.0; // (1.0 for uniform mesh, 0.25 for coarse mesh)      
      if (type_operation == 202)
	main_file<<"4 0 45.0 0.04 0.2 30.0 20.0 " << munif <<" 0.4 10 200 " << nodelem<<" 1 0.01 2 2 0.3"<<endl;
      else if (type_operation == 203)
	main_file<<"3 0 45.0 0.04 0.2 0.5 1.0 " << nodelem << " 3 2 2 0.3"<<endl;
      else if (type_operation == 207)
	main_file<<"3 0 45.0 0.04 0.2  1.0 " << nodelem << " 3 2 2 0.3"<<endl;
      
      main_file.close();
      
      // we retrieve all the vertices of the boundary
      Vector<int> VertexBoundary(mesh.GetNbVertices()); VertexBoundary.Fill(-1);
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	{
	  VertexBoundary(mesh.BoundaryRef(i).numVertex(0)) = 1;
	  VertexBoundary(mesh.BoundaryRef(i).numVertex(1)) = 1;
	}
      
      int nb_vertices_ref = 0;
      for (int i = 0; i < mesh.GetNbVertices(); i++)
	if (VertexBoundary(i) != -1)
	  {
	    VertexBoundary(i) = nb_vertices_ref;
	    nb_vertices_ref++;
	  }
      
      region_file<<nb_vertices_ref<<endl;
      
      // each vertex of the boundary is constrained, should be in the final mesh
      for (int i = 0; i < mesh.GetNbVertices(); i++)
	if (VertexBoundary(i) != -1)
	  region_file<<mesh.Vertex(i)(0) << " " << mesh.Vertex(i)(1) <<" 2\n";
      
      // no extra-vertex
      region_file<<"\n 0 \n\n";
      
      // we retrieve all the loops
      Vector<IVect> all_loops; IVect type_loop;
      FindAllLoops(mesh, all_loops, type_loop);
      region_file<<all_loops.GetM()<<"\n";
      for (int i = 0; i < all_loops.GetM(); i++)
	{
	  int regcode = 0;
	  region_file<<all_loops(i).GetM()<<" "<<regcode<<" "<<type_loop(i)<<'\n';
	  for (int j = 0; j < all_loops(i).GetM(); j++)
	    {
	      region_file<<VertexBoundary(all_loops(i)(j))+1<<" ";
	      if ((j+1)%20 == 0)
		region_file<<'\n';
	    }
	  
	  region_file<<'\n';
	  for (int j = 0; j < all_loops(i).GetM(); j++)
	    {
	      region_file<<"1 ";
	      if ((j+1)%20 == 0)
		region_file<<'\n';
	    }
	}
  
      region_file.close();
      curve_file<<"0"<<endl;
      curve_file.close();
      
      // we generate a mesh file with geompack
      command_line = string("zgp0803_linux main_input test.out");
      ierr = system(command_line.data());
      if (ierr != 0)
	{
	  cout<<"Error during the generation of mh2 file"<<endl;
	  cout<<"Do you have an even number of points on the boundary for all-quadrilateral meshing ?"<<endl;
	  cout<<"Number of points on the boundary "<<nb_vertices_ref<<endl;
	  abort();
	}
      
      Mesh<Dimension2> output_mesh;
      output_mesh.Read(string("mesh_output.mh2"));
      
      // we replace referenced boundaries with initial ones
      int nb_boundary_ref = 0;
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	if (OldReferences(i) != ref_unused)
	  nb_boundary_ref++;
      
      output_mesh.ReallocateBoundariesRef(nb_boundary_ref); nb = 0;
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	if (OldReferences(i) != ref_unused)
	  {
	    output_mesh.BoundaryRef(nb) = mesh.BoundaryRef(i);
	    output_mesh.BoundaryRef(nb).Init(VertexBoundary(mesh.BoundaryRef(i).numVertex(0)),
					     VertexBoundary(mesh.BoundaryRef(i).numVertex(1)));
	    
	    output_mesh.BoundaryRef(nb).SetReference(OldReferences(i));
	    nb++;
	  }
      
      // chaning references of elements
      int ref = ref_surfaces(num_surfaces(isurf));
      for (int i = 0; i < output_mesh.GetNbElt(); i++)
	output_mesh.Element(i).SetReference(ref);
      
      // we append this mesh to the global mesh
      all_mesh.AppendMesh(output_mesh, true);
      
    }
  
  // checking if symmetry is needed
  for (int i = 0; i < nb_lines; i++)
    {
      if (lignes(i).find("SymmetryX") == 0)
	all_mesh.SymmetrizeMesh(0, 0.0, R2(0, 0));
      
      if (lignes(i).find("SymmetryY") == 0)
	all_mesh.SymmetrizeMesh(1, 0.0, R2(0, 0));
    }
  
  // writing final mesh obtained
  all_mesh.Write(file_mesh);
  
  return FinalizeMontjoie();
}  
