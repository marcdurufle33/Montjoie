#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int GetBoundaryVertices(const Face<Dimension3>& face, const IVect& Index, bool use_index,
			int& n1, int& n2, int& n3, int& n4)
{
  n4 = -1;
  if (face.GetNbVertices() == 3)
    {
      n1 = face.numVertex(0);
      n2 = face.numVertex(1);
      n3 = face.numVertex(2);
    }
  else
    {
      n1 = face.numVertex(0);
      n2 = face.numVertex(1);
      n3 = face.numVertex(2);
      n4 = face.numVertex(3); 
    }

  if (use_index)
    {
      n1 = Index(n1);
      n2 = Index(n2);
      n3 = Index(n3);

      if (face.GetNbVertices() == 4)
	n4 = Index(n4);
    }

  if (face.GetNbVertices() == 3)
    Sort(n1, n2, n3);
  else
    Sort(n1, n2, n3, n4);
  
  return face.GetReference();
}

void CopySurfaceMesh(const Mesh<Dimension3>& mesh_surf,
		     Mesh<Dimension3>& mesh_vol)
{
  VectR3 vert_surf(mesh_surf.Vertex());
  VectR3 vert_vol(mesh_vol.Vertex());

  Vector<int> permut_surf(vert_surf.GetM()), permut_vol(vert_vol.GetM());
  permut_surf.Fill(); permut_vol.Fill();

  Sort(vert_surf, permut_surf);
  Sort(vert_vol, permut_vol);

  Vector<int> num_vertex_surf(vert_surf.GetM());
  Vector<int> IndexVertex(vert_vol.GetM());
  IndexVertex.Fill(-1);
  int k = 0;
  for (int i = 0; i < permut_surf.GetM(); i++)
    {
      DISP(i); DISP(permut_surf(i)); DISP(vert_surf(i));
      DISP(vert_vol(k));
      while ((k < vert_vol.GetM()) && (vert_vol(k) < vert_surf(i)))
	k++;

      if ((k < vert_vol.GetM()) && (vert_vol(k) == vert_surf(i)))
	{
	  num_vertex_surf(permut_surf(i)) = permut_vol(k);
	  IndexVertex(permut_vol(k)) = permut_surf(i);
	}
      else
	{
	  cout << "Vertex not localized" << endl;
	  cout << i << " : " << vert_surf(i) << endl;
	  abort();
	}
    }

  // changing references of faces
  Vector<int> num, all_num;
  for (int i = 0; i < mesh_vol.GetNbBoundaryRef(); i++)
    {
      int nb_vert = mesh_vol.BoundaryRef(i).GetNbVertices();
      int n1, n2, n3, n4;
      GetBoundaryVertices(mesh_vol.BoundaryRef(i), IndexVertex, false,
			  n1, n2, n3, n4);
      
      all_num.Clear();
      for (int j = 0; j < nb_vert; j++)
	{
	  int nv = mesh_vol.BoundaryRef(i).numVertex(j);
	  nv = IndexVertex(nv);
	  if (nv == -1)
	    {
	      cout << "Vertex not localized" << endl;
	      cout << nv << " : " << mesh_vol.Vertex(nv) << endl;
	      abort();
	    }
	  
	  mesh_surf.GetNeighboringBoundariesAroundVertex(nv, num);
	  all_num.PushBack(num);
	}

      RemoveDuplicate(all_num);
      DISP(i); DISP(all_num);
      
      for (int j = 0; j < all_num.GetM(); j++)
	{
	  int nf = all_num(j);
	  int n1b, n2b, n3b, n4b;
	  int ref
	    = GetBoundaryVertices(mesh_surf.BoundaryRef(nf), IndexVertex, true,
				  n1b, n2b, n3b, n4b);

	  if ((n1 == n1b) && (n2 == n2b) && (n3 == n3b) && (n4 == n4b))
	    {
	      mesh_vol.BoundaryRef(i).SetReference(ref);
	      break;
	    }
	}
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  if (argc != 3)
    {
      cout << "Invalid number of arguments" << endl;
      cout << "Usage : ./copy_surface surface.mesh vol.mesh" << endl;
      abort();
    }
  
  string name_surf(argv[1]);
  string name_vol(argv[2]);

  Mesh<Dimension3> mesh_surf, mesh_vol;
  mesh_surf.Read(name_surf);
  mesh_vol.Read(name_vol);

  CopySurfaceMesh(mesh_surf, mesh_vol);

  mesh_vol.Read("toto.mesh");
  
  return FinalizeMontjoie();
}
