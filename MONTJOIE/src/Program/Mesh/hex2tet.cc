#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc != 3)
    {
      cout<<"Cette commande demande deux arguments"<<endl;
      cout<<"hex2tet hex.mesh tet.mesh"<<endl;
      cout<<"the file tet.mesh is created"<<endl;
      return -1;
    }
  
  string file_input(argv[1]);
  string file_output(argv[2]);
  
  Mesh<Dimension3> mesh, mesh_tet;
  
  // on lit le maillage
  mesh.Read(file_input);
  
  if (mesh.GetNbTetrahedra() > 0)
    {
      cout<<"Maillage partiellement compose de tetraedres/triangles "<<endl;
      cout<<"On ne peut pas revenir en arriere"<<endl;
      return -2;
    }
  
  // on copie tous les sommets
  mesh_tet.ReallocateVertices(mesh.GetNbVertices());
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    mesh_tet.Vertex(i) = mesh.Vertex(i);
  
  // on construit tous les triangles
  IVect num(12);
  mesh_tet.ReallocateBoundariesRef(mesh.GetNbBoundaryRef()/3);
  for (int i = 0; i < (mesh.GetNbBoundaryRef()/3); i++)
    {
      int nb = 0, ref = mesh.Boundary(3*i).GetReference();
      for (int j = 0; j < 3; j++)
	for (int k = 0; k < 4; k++)
	  num(nb++) = mesh.BoundaryRef(3*i+j).numVertex(k);
      
      // on trie
      Sort(nb, num);
      
      // les trois premiers indices sont les sommets du triangles
      mesh_tet.BoundaryRef(i).InitTriangular(num(0), num(1), num(2), ref);
    }
  
  // on construit tous les tetraedres
  num.Reallocate(32);
  mesh_tet.ReallocateElements(mesh.GetNbElt()/4);
  for (int i = 0; i < (mesh.GetNbElt()/4); i++)
    {
      int nb = 0, ref = mesh.Element(4*i).GetReference();
      for (int j = 0; j < 4; j++)
	for (int k = 0; k < 8; k++)
	  num(nb++) = mesh.Element(4*i+j).numVertex(k);
      
      // on trie
      Sort(nb, num);
      
      // les quatre premiers indices sont les sommets du tetraedre
      mesh_tet.Element(i).InitTetrahedral(num(0), num(1), num(2), num(3), ref);
    }
  
  // on elimine tous les sommets inutiles
  mesh_tet.ForceCoherenceMesh();
  
  // on ecrit le maillage
  mesh_tet.Write(file_output);
  
  return FinalizeMontjoie();
}
