#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  
  InitMontjoie(argc, argv);

  if (argc != 3)
    {
      cout<<"Cette commande demande trois arguments"<<endl;
      cout<<"locate_face toto.mesh 2"<<endl;
      cout<<"the file face_location.mesh is created with the face and neighbouring faces (or elements) "<<endl;
    }
  
  Mesh<Dimension3> mesh;
  mesh.print_level = 3;
  // lecture du maillage
  mesh.Read(string(argv[1]));
  int level = atoi(argv[2]);
  
  cout << "Voulez vous ? " << endl;
  cout << "1- Donner des numeros de points " << endl;
  cout << "2- Donner des numeros de faces " << endl;
  cout << "3- Donner des references " << endl;
  cout << "4- Donner des coordonnees et un rayon de distribution " << endl;
  
  int type_input; cin >> type_input;
  string ligne; Vector<string> param;
  IVect num, num_vertex; list<int> face_to_add;
  string delim(", ");
  getline(cin, ligne);
  
  switch (type_input)
    {
    case 1 :
      {
        cout << "Entrez la liste des sommets " << endl;
        getline(cin, ligne);
        StringTokenize(ligne, param, delim);
        num.Reallocate(param.GetM());
        for (int i = 0; i < param.GetM(); i++)
          num(i) = to_num<int>(param(i)) - 1;
        
        for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
          {
            int nb_vertices_loc = mesh.BoundaryRef(i).GetNbVertices();
            bool add_face = false;
            for (int j = 0; j < nb_vertices_loc; j++)
              {
                int nv = mesh.BoundaryRef(i).numVertex(j);
                for (int k = 0; k < num.GetM(); k++)
                  if (nv == num(k))
                    add_face = true;
              }
            
            if (add_face)
              face_to_add.push_back(i);
            
          }
      }
      break;

    case 2 :
      {
        cout << "Entrez la liste des faces " << endl;
        getline(cin, ligne);
        StringTokenize(ligne, param, delim);
        for (int i = 0; i < param.GetM(); i++)
          face_to_add.push_back(to_num<int>(param(i)) - 1);
      }
      break;

    case 3 :
      {
        cout << "Entrez la liste des references " << endl;
        getline(cin, ligne);
        StringTokenize(ligne, param, delim);
        num.Reallocate(param.GetM());
        for (int i = 0; i < param.GetM(); i++)
          num(i) = to_num<int>(param(i));
        
        for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
          {
            bool add_face = false;
            for (int k = 0; k < num.GetM(); k++)
              if (mesh.BoundaryRef(i).GetReference() == num(k))
                add_face = true;
            
            if (add_face)
              face_to_add.push_back(i);
          }
      }
      break;

    case 4 :
      {
        cout << "Entrez le centre et rayon : x y z r " << endl;
        R3 coor; Real_wp radius;
        cin >> coor(0) >> coor(1) >> coor(2) >> radius;
        
        for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
          {
            bool add_face = false;
            int nb_vertices_loc = mesh.BoundaryRef(i).GetNbVertices();
            for (int j = 0; j < nb_vertices_loc; j++)
              {
                int nv = mesh.BoundaryRef(i).numVertex(j);
                if (coor.Distance(mesh.Vertex(nv)) < radius)
                  add_face = true;
              }

            if (add_face)
              face_to_add.push_back(i);
          }
      }
      break;
      
    default :
      {
        cout << "Choix non autorise " << endl;
        return -1;
      }
      
    }
      
  IVect num_face;   Copy(face_to_add, num_face);
  IVect num_edge; IVect num_elt; int nb_elt = 0;
  // on itere level fois, la recuperation des faces voisines
  for (int num_level = 0; num_level < level; num_level++)
    {
      int nb_faces = num_face.GetM();
      // on boucle sur toutes les faces, pour recuperer toutes les aretes
      num_edge.Reallocate(4*nb_faces);
      int nb_edges = 0;
      for (int i = 0; i < nb_faces; i++)
	for (int j = 0; j < mesh.Boundary(num_face(i)).GetNbEdges(); j++)
	  {
	    num_edge(nb_edges) = mesh.Boundary(num_face(i)).numEdge(j);
	    nb_edges++;
	  }
      
      // on elimine les doublons
      Assemble(nb_edges, num_edge);
      
      // on recupere toutes les faces autour de ces aretes
      nb_faces = 0;
      num_face.Reallocate(12*nb_edges);
      for (int i = 0; i < nb_edges; i++)
	for (int j = 0; j < mesh.GetEdge(num_edge(i)).GetNbFaces(); j++)
	  {
            int nf = mesh.GetEdge(num_edge(i)).numFace(j);
            num_face(nb_faces) =  nf;
            nb_faces++;
          }
      
      // on elimine doublons
      Assemble(nb_faces, num_face); num_face.Resize(nb_faces);
      
      // on recupere tous les elements autour de ces faces
      nb_elt = 0;
      num_elt.Reallocate(2*nb_faces);
      for (int i = 0; i < nb_faces; i++)
        for (int j = 0; j < mesh.Boundary(num_face(i)).GetNbElements(); j++)
          num_elt(nb_elt++) = mesh.Boundary(num_face(i)).numElement(j);
      
      // on elimine doublons
      Assemble(nb_elt, num_elt); num_elt.Resize(nb_elt);
      
    }
  
  // sommets mis en jeu ?
  VectBool VertexTaken(mesh.GetNbVertices()); VertexTaken.Fill(false);
  for (int i = 0; i < num_face.GetM(); i++)
    for (int j = 0; j < mesh.Boundary(num_face(i)).GetNbVertices(); j++)
      VertexTaken(mesh.Boundary(num_face(i)).numVertex(j)) = true;

  for (int i = 0; i < num_elt.GetM(); i++)
    for (int j = 0; j < mesh.Element(num_elt(i)).GetNbVertices(); j++)
      VertexTaken(mesh.Element(num_elt(i)).numVertex(j)) = true;
  
  // on calcule xmin, ymin, zmin, xmax, ymax, zmax
  Real_wp xmin = 1e300, ymin = 1e300, zmin = 1e300, xmax = -1e300, ymax = -1e300, zmax = -1e300;
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    if (VertexTaken(i))
      {
	xmin = min(xmin, mesh.Vertex(i)(0)); xmax = max(xmax, mesh.Vertex(i)(0));
	ymin = min(ymin, mesh.Vertex(i)(1)); ymax = max(ymax, mesh.Vertex(i)(1));
	zmin = min(zmin, mesh.Vertex(i)(2)); zmax = max(zmax, mesh.Vertex(i)(2));
      }

  int nb_new_vertices = 0;
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    if (VertexTaken(i))
      nb_new_vertices++;

  Vector<bool> ElementTaken(mesh.GetNbElt());
  ElementTaken.Fill(false);
  for (int i = 0; i < num_elt.GetM(); i++)
    ElementTaken(num_elt(i)) = true;

  // creation d'un maillage extrait  
  Mesh<Dimension3> sub_mesh;
  mesh.CreateSubmesh(sub_mesh, nb_new_vertices, nb_elt, VertexTaken, ElementTaken);  

  sub_mesh.FindConnectivity();
  sub_mesh.AddBoundaryFaces();
  //sub_mesh.ClearElements();
  sub_mesh.Write("sub.mesh");
  
  // position des points "fantomes"
  R3 pt_accumulation; pt_accumulation(0) = xmin + 0.5*(xmin-xmax);
  pt_accumulation(1) = ymin + 0.5*(ymin-ymax); pt_accumulation(2) = zmin + 0.5*(zmin-zmax);
  
  int nb_vertices = mesh.GetNbVertices();
  mesh.ResizeVertices(nb_vertices+8);
  for (int i = nb_vertices; i < mesh.GetNbVertices(); i++)
    mesh.Vertex(i) = pt_accumulation;
  
  // on modifie le maillage

  // sommets du nouveau maillage
  for (int i = 0; i < nb_vertices; i++)
    {
      if (!VertexTaken(i))
        mesh.Vertex(i) = pt_accumulation;
    }
  
  // modification faces
  int p = nb_vertices;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      bool face_taken = true;
      for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
        if (!VertexTaken(mesh.BoundaryRef(i).numVertex(j)))
          face_taken = false;
      
      if (!face_taken)
        mesh.BoundaryRef(i).InitTriangular(p, p+1, p+2, 1);
      
    }

  // modification elements
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      bool face_taken = true;
      for (int j = 0; j < mesh.Element(i).GetNbVertices(); j++)
        if (!VertexTaken(mesh.Element(i).numVertex(j)))
          face_taken = false;
      
      if (!face_taken)
        mesh.Element(i).InitTetrahedral(p, p+1, p+2, p+3, 1);
      
    }
  
  mesh.Write(string("face_location.mesh"));
  
  return FinalizeMontjoie();
}
