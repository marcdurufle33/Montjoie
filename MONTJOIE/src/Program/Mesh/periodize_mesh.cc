#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

template<class Dimension>
void ChangeReference(Mesh<Dimension>& mesh, int ref_max, int ref_to_increment,
		     int& new_ref)
{
  if (ref_to_increment <= 0)
    return;
  
  int nb_ref = mesh.GetNbReferences();
  IVect ref_to_affect(nb_ref+1);
  ref_to_affect.Fill(-1);
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      if ((ref == ref_to_increment) || (ref >= ref_max))
	{
	  if (ref_to_affect(ref) == -1)
	    ref_to_affect(ref) = new_ref++;
	      
	  mesh.BoundaryRef(i).SetReference(ref_to_affect(ref));
	}
    }
  
  if (new_ref > nb_ref)
    mesh.ResizeNbReferences(new_ref+1);
}

template<class Dimension>
void TransformMesh(Mesh<Dimension>& mesh, bool cyclic,
                   int nx_min, int nx_max, int ny_min, int ny_max, int nz_min, int nz_max,
		   int ref_to_increment, const Real_wp& dx_ref, const Real_wp& dy_ref, const Real_wp& dz_ref)
{
  // we find all the references on the boundary
  Vector<bool> ref_on_boundary(mesh.GetNbReferences()+1);
  ref_on_boundary.Fill(false);
  int new_ref = 0;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      if (ref > 0)
        {
	  if (mesh.BoundaryRef(i).GetNbElements() == 1)
	    ref_on_boundary(ref) = true;
	  
	  new_ref = max(new_ref, ref);
	}
    }

  new_ref++; int ref_max = new_ref;
  if (ref_to_increment > 0)
    cout << "First incremented reference = " << new_ref << endl;
  
  typename Dimension::R_N vec_u;
  if (nx_min <= nx_max)
    {
      Real_wp dx = mesh.GetXmax() - mesh.GetXmin();
      if (dx_ref > 0)
	dx = dx_ref;
      
      Mesh<Dimension> mesh_cell, mesh_tmp;
      mesh_cell = mesh;
      
      // first translation
      vec_u(0) = nx_min*dx;
      TranslateMesh(mesh, vec_u);
      
      for (int n = nx_min+1; n <= nx_max; n++)
        {
          // then appending other meshes
          mesh_tmp = mesh_cell;
          vec_u(0) = n*dx;
          TranslateMesh(mesh_tmp, vec_u);
	  ChangeReference(mesh_tmp, ref_max, ref_to_increment, new_ref);
          mesh.AppendMesh(mesh_tmp);
        }
      
      vec_u(0) = 0;
    }

  if (ny_min <= ny_max)
    {
      Real_wp dy = mesh.GetYmax() - mesh.GetYmin();
      if (dy_ref > 0)
	dy = dy_ref;

      Mesh<Dimension> mesh_cell, mesh_tmp;
      mesh_cell = mesh;
      
      vec_u(1) = ny_min*dy;
      TranslateMesh(mesh, vec_u);
      for (int n = ny_min+1; n <= ny_max; n++)
        {
          mesh_tmp = mesh_cell;
          vec_u(1) = n*dy;
          TranslateMesh(mesh_tmp, vec_u);
	  ChangeReference(mesh_tmp, ref_max, ref_to_increment, new_ref);
          mesh.AppendMesh(mesh_tmp);
        }
      vec_u(1) = 0;
    }

  if (nz_min <= nz_max)
    {
      Real_wp dz = mesh.GetZmax() - mesh.GetZmin();
      if (dz_ref > 0)
	dz = dz_ref;

      Mesh<Dimension> mesh_cell, mesh_tmp;
      mesh_cell = mesh;
      
      vec_u(2) = nz_min*dz;
      TranslateMesh(mesh, vec_u);
      for (int n = nz_min+1; n <= nz_max; n++)
        {
          mesh_tmp = mesh_cell;
          vec_u(2) = n*dz;
          TranslateMesh(mesh_tmp, vec_u);
	  ChangeReference(mesh_tmp, ref_max, ref_to_increment, new_ref);
          mesh.AppendMesh(mesh_tmp);
        }      
    }
  
  vec_u.Zero();
  if (cyclic)
    mesh.PeriodizeMeshTeta(vec_u);
  
  // we are removing referenced edges shared by two elements
  // which were on the boundary before
  int nb_new = mesh.GetNbBoundaryRef();
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      if ((ref > 0) && (ref < ref_on_boundary.GetM()))
        {
          if (ref_on_boundary(ref))
            if (mesh.BoundaryRef(i).GetNbElements() == 2)
              {
                mesh.BoundaryRef(i).SetReference(0);
                nb_new--;
              }
        }      
    }
  
  if (nb_new != mesh.GetNbBoundaryRef())
    mesh.ClearZeroBoundaryRef();
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  if ((argc == 2) && ((string(argv[1]) == "--help")
		      || (string(argv[1]) == "-help")))
    {
      cout << "Usage :" << endl;
      cout << "./per_mesh.x input.mesh output.mesh 1 ..." << endl;
      cout << endl;
      cout << "The first argument (here input.mesh) is the name of the file containing the mesh of a cell " << endl;
      cout << "The second argument (here output.mesh) is the name of the output file" << endl;
      cout << "The third argument (here 1) is the order of the produced mesh" << endl;
      cout << "If equal to 1, the mesh is not curved" << endl << endl;
      cout << "The arguments to put after this integer are one the followings:" << endl;
      cout << "(They can be given with any order)" << endl << endl;
      cout << "x n0 n1    :  periodicity in x, the x-coordinate will range from xmin + n0*dx until xmax + n1*dx " << endl << endl;
      cout << "y n0 n1    :  periodicity in y, the y-coordinate will range from ymin + n0*dy until ymax + n1*dy " << endl << endl;
      cout << "z n0 n1    :  periodicity in z, the z-coordinate will range from zmin + n0*dz until zmax + n1*dz " << endl << endl;
      cout << "dx val    :  value for dx (if not mentioned, equal to xmax-xmin) " << endl << endl;
      cout << "dy val    :  value for dy (if not mentioned, equal to ymax-ymin) " << endl << endl;
      cout << "dz val    :  value for dz (if not mentioned, equal to zmax-zmin) " << endl << endl;
      cout << "r ref    :  if the surfaces of reference ref must be incremented " << endl << endl;
      cout << "cyclic    :  periodicity in theta " << endl << endl;

      cout << endl << "Examples : " <<endl;
      cout << "./per_mesh.x input.mesh output.mesh 1 z -2 2 cyclic" << endl
	   << "        will produce a mesh periodized in theta and z with 5 cells in z" << endl << endl;
      cout << "./per_mesh.x input.mesh output.mesh 1 y 1 3 x -4 -1 r 2" << endl
	   << "        will produce a mesh periodized in x and y with three cells in y and four cells in x," << endl
	   << "         the reference 2 will be incremented for each new cell " << endl << endl;

      cout << "./per_mesh.x input.mesh output.mesh 1 z -2 2 dz 1.0" << endl
	   << "        will produce a mesh periodized in z, the height for each cell is equal to 1.0 (it will not be computed from the mesh) " << endl << endl;
      
      
      exit(0);
    }
  
  if (argc < 5)
    {
      cout<<"Cette commande demande au moins quatre arguments"<<endl;
      cout<<"per_mesh source.mesh result.mesh order cyclic"<<endl;
      cout<<"the file result.mesh is created"<<endl;
      abort();
    }
  
  bool cyclic = false;
  int nx_min = 1, nx_max = 0, ny_min = 1, ny_max = 0, nz_min = 1, nz_max = 0;
  int nb = 4; int ref_to_increment = -1;
  Real_wp dx_ref(-1), dy_ref(-1), dz_ref(-1);
  while (argc > nb)
    {
      string type_per(argv[nb]);
      nb++;
      if (!type_per.compare("x"))
        {
          if (argc > nb+1)
            {
              nx_min = atoi(argv[nb]);
              nx_max = atoi(argv[nb+1]);
              nb += 2;
            }
          else
            {
              cout << "Specify the interval of periodicity " << endl;
              abort();
            }
        }

      if (!type_per.compare("dx"))
        {
	  dx_ref = atof(argv[nb]);
	  nb++;
	}

      if (!type_per.compare("dy"))
        {
	  dy_ref = atof(argv[nb]);
	  nb++;
	}

      if (!type_per.compare("dz"))
        {
	  dz_ref = atof(argv[nb]);
	  nb++;
	}

      if (!type_per.compare("y"))
        {
          if (argc > nb+1)
            {
              ny_min = atoi(argv[nb]);
              ny_max = atoi(argv[nb+1]);
              nb += 2;
            }
          else
            {
              cout << "Specify the interval of periodicity " << endl;
              abort();
            }
        }
      
      if (!type_per.compare("z"))
        {
          if (argc > nb+1)
            {
              nz_min = atoi(argv[nb]);
              nz_max = atoi(argv[nb+1]);
              nb += 2;
            }
          else
            {
              cout << "Specify the interval of periodicity " << endl;
              abort();
            }
        }
      
      if (type_per == "r")
	{
	  if (argc > nb)
	    {
	      ref_to_increment = atoi(argv[nb]);
	      nb++;
	    }
          else
            {
              cout << "Specify the reference to increment " << endl;
              abort();
            }
	}	  
	  
      if (!type_per.compare("cyclic"))
        cyclic = true;
    }

  string file_input(argv[1]);
  string file_output(argv[2]);
  int order = atoi(argv[3]);
  string extension_input = GetExtension(file_input);
  
  int dimension = GetDimensionMesh(file_input, extension_input);
  if (dimension == 2)
    {
      Mesh<Dimension2> mesh;
      mesh.SetGeometryOrder(order);
      
      // on lit le maillage
      mesh.Read(file_input);
      
      TransformMesh(mesh, cyclic, nx_min, nx_max, ny_min, ny_max, nz_min, nz_max, ref_to_increment,
		    dx_ref, dy_ref, dz_ref);
      
      // on ecrit maillage
      mesh.Write(file_output);
    }
  else
    {
      Mesh<Dimension3> mesh;
      mesh.SetGeometryOrder(order);
      
      // on lit le maillage
      mesh.Read(file_input);
      
      TransformMesh(mesh, cyclic, nx_min, nx_max, ny_min, ny_max, nz_min, nz_max, ref_to_increment,
		    dx_ref, dy_ref, dz_ref);
      
      // on ecrit maillage
      mesh.Write(file_output);
    }
  
  return FinalizeMontjoie();
}

