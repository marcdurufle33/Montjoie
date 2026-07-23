#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  if (argc < 2)
    {
      cout << "Entrez un nom de fichier .geo" << endl;
      cout << "Exemple : ./create_geo test.geo" << endl;
      abort();
    }
  
  Mesh<Dimension3> mesh;
  string nom_mesh;
  
  int offset = 0, ns = 0;
  cout << "De quel fichier .mesh veut on extraire les points " << endl;
  getline(cin, nom_mesh);
  
  mesh.Read(nom_mesh);
  
  IVect IndexVertex(mesh.GetNbVertices());
  IndexVertex.Fill(-1);
  
  // reconstructing IndexVert from previous geo file
  {
    string ligne;
    
    ifstream file_geo(argv[1]);
    if (!file_geo.is_open())
      {
	cout << "File " << argv[1] << "does not exist ?" << endl;
	abort();
      }
    
    list<R3> all_pts;
    list<int> num;
    while (!file_geo.eof())
      {
	getline(file_geo, ligne);
	
	if (file_geo.good())
	  {
	    if (ligne.substr(0, 5) == string("Point"))
	      {
		int pos = ligne.find(')');
		int n = to_num<int>(ligne.substr(6, pos-6));
		Real_wp x, y, z;
		int posA = ligne.find('{');
		pos = ligne.find(',', posA);
		x = to_num<Real_wp>(ligne.substr(posA+1, pos-posA-1));
		posA = pos; pos = ligne.find(',', posA+1);
		y = to_num<Real_wp>(ligne.substr(posA+1, pos-posA-1));
		posA = pos; pos = ligne.find(',', posA+1);
		z = to_num<Real_wp>(ligne.substr(posA+1, pos-posA-1));
		all_pts.push_back(R3(x, y, z));
		num.push_back(n);
		
		if (n > offset)
		  offset = n;
	      }
	    else if (ligne.substr(0, 4) == string("Line"))
	      {
		int pos = ligne.find(')');
		int n = to_num<int>(ligne.substr(5, pos-5));
		if (n > ns)
		  ns = n;
	      }
	    else if (ligne.substr(0, 6) == string("Spline"))
	      {
		int pos = ligne.find(')');
		int n = to_num<int>(ligne.substr(7, pos-7));
		if (n > ns)
		  ns = n;
	      }
	  }
      }
    
    file_geo.close();
    
    VectR3 Vertices;
    IVect Num;
    Copy(all_pts, Vertices);
    Copy(num, Num);
    
    // sorting vertices so that they can be quickly compared
    VectR3 VertMesh(mesh.GetNbVertices());
    IVect permut(mesh.GetNbVertices());
    for (int i = 0; i < VertMesh.GetM(); i++)
      VertMesh(i) = mesh.Vertex(i);
    
    permut.Fill();
    Sort(VertMesh, permut);
    
    // sorting vertices of geo file
    Sort(Vertices, Num);
    
    int n = 0;
    for (int p = 0; p < mesh.GetNbVertices(); p++)
      {
	while ((n < Num.GetM()) && (Vertices(n) < VertMesh(p)))
	  n++;
	
	if (n < Num.GetM())
	  if (Vertices(n) == VertMesh(p))
	    {
	      IndexVertex(permut(p)) = Num(n)-1;
	    }
      }
    
    file_geo.close();
  }
  
  DISP(offset); DISP(ns);
  while (true)
    {
      string chaine_num;
      cout << "Entrez la liste des points a rajouter " << endl;
      getline(cin, chaine_num);
      
      IVect num;
      istringstream flux(chaine_num);
      num.ReadText(flux);
      
      ofstream file_geo(argv[1], ios::app);
      file_geo.precision(8);
      
      for (int i = 0; i < num.GetM(); i++)
	{
	  int p = num(i) - 1;
	  int n = offset;
	  if (IndexVertex(p) != -1)
	    n = IndexVertex(p);
	  else
	    {
	      file_geo << "Point(" << n+1 << ") = {"
		       << mesh.Vertex(p)(0) << ", "
		       << mesh.Vertex(p)(1) << ", "
		       << mesh.Vertex(p)(2) << ", lc};\n";
	      
	      IndexVertex(p) = offset;
	      offset++;
	    }
	  
	  num(i) = n+1;
	}
      
      file_geo << '\n';
      file_geo << "Spline(" << ns+1 << ") = {" << num(0);
      for (int i = 1; i < num.GetM(); i++)
	file_geo << ", " << num(i);
      
      file_geo << "};\n";
      file_geo << '\n';
      
      ns++;
      file_geo.close();
    }
  
  return 0;
}
