#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

void CheckPermutation(const IVect& num2)
{
  IVect num(num2);
  Sort(num);
  for (int i = 0; i < num.GetM(); i++)
    if (num(i) != i)
      {
	cout << " Not a permutation array" << endl;
	abort();
      }
}

template<class T>
void ReadParameters(istream& file_in, Vector<T>& param)
{
  string chaine;
  // cin>>flush;
  // getline(cin, chaine ,'\n');
  getline(file_in, chaine ,'\n');
  istringstream stream_data(chaine);
  
  param.ReadText(stream_data);
}

void ReadParameters(istream& file_in, string& param)
{
  getline(file_in, param ,'\n');
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  R2::threshold = 1e-8;
  R3::threshold = 1e-8;
  
  IVect num, num_gmsh;
  Matrix<int> coor, NumNodes2D;
  string command_line;
  Vector<double> coord;
  VectR2 Points2D;
  
  {
    cout << "Checking Triangles" << endl;
    
    // checking high order triangles
    for (int r = 2; r <= 10; r++)
      {
	Vector<double> step(r+1);
	step.Fill(); Mlt(1.0/r, step);
	
	ofstream file_out("triangle.geo");
	file_out << "lc = 10.0;\n";
	file_out << "Point(1) = {0, 0, 0, lc};\n";
	file_out << "Point(2) = {1, 0, 0, lc};\n";
	file_out << "Point(3) = {0, 1, 0, lc};\n";
	file_out << "Line(1) = {1, 2};\n Line(2) = {2, 3};\n Line(3) = {3, 1};\n";
	file_out << "Line Loop(1) = {1,2,3} ;\n";
	file_out << "Transfinite Line {1,2,3} = 2;\n";
	file_out << "Plane Surface(1) = {1};\n";
	file_out << "Transfinite Surface {1} = {1, 2, 3};\n";
	file_out << "Physical Line(1) = {1, 2, 3} ;\n";
	file_out << "Mesh.SecondOrderIncomplete = 0;\n";
	file_out << "Physical Surface(1) = {1} ;\n";
	file_out.close();
	
	NumberMap::GetGmshTriangularNumbering(r, num_gmsh);
	MeshNumbering<Dimension2>::ConstructTriangularNumbering(r, NumNodes2D, coor);
	CheckPermutation(num_gmsh);
	
	string command_line = string("geo2mesh triangle ") + to_str(r) + " > sortie";
	system(command_line.data());
	
	ifstream file_in("triangle.msh");
	string ligne;
	getline(file_in, ligne);
	getline(file_in, ligne);
	getline(file_in, ligne);
	getline(file_in, ligne);
	ReadParameters(file_in, num); int nb_vertices = num(0);
	Points2D.Reallocate(nb_vertices);
	for (int i = 0; i < nb_vertices; i++)
	  {
	    ReadParameters(file_in, coord);
	    Points2D(i).Init(coord(1), coord(2));
	  }
	
	getline(file_in, ligne);
	getline(file_in, ligne);
	ReadParameters(file_in, num); int nb_elt = num(0);
	for (int i = 0; i < nb_elt; i++)
	  {
	    ReadParameters(file_in, num);
	    if (num(0) == 4)
	      {
		for (int j = 0; j < num_gmsh.GetM(); j++)
		  {
		    int n = num(5+j);
		    int node = num_gmsh(j);
		    int i1 = coor(node, 1);
		    int i2 = coor(node, 2);
		    double x = step(i1);
		    double y = step(i2);
		    if (Points2D(n-1) != R2(x, y))
		      {
			DISP(r); DISP(j); DISP(node); DISP(i1); DISP(i2);
			DISP(x); DISP(y); DISP(Points2D(n-1));
			cout << "point not coincident" << endl;
		      abort();
		      }
		  }
	      }
	  }
      }
    
  }


  {
    cout << "Checking Quadrangles" << endl;
    
    // checking high order quadrangles
    for (int r = 2; r <= 10; r++)
      {
	Vector<double> step(r+1);
	step.Fill(); Mlt(1.0/r, step);
	
	ofstream file_out("quadrangle.geo");
	file_out << "lc = 10.0;\n";
	file_out << "Point(1) = {0, 0, 0, lc};\n";
	file_out << "Point(2) = {1, 0, 0, lc};\n";
	file_out << "Point(3) = {1, 1, 0, lc};\n";
	file_out << "Point(4) = {0, 1, 0, lc};\n";
	file_out << "Line(1) = {1, 2};\n Line(2) = {2, 3};\n Line(3) = {3, 4};\n Line(4) = {4, 1};\n";
	file_out << "Line Loop(1) = {1,2,3,4} ;\n";
	file_out << "Transfinite Line {1,2,3,4} = 2;\n";
	file_out << "Plane Surface(1) = {1};\n";
	file_out << "Transfinite Surface {1} = {1, 2, 3, 4};\n";
	file_out << "Recombine Surface {1} = 80;\n";
	file_out << "Physical Line(1) = {1, 2, 3, 4} ;\n";
	file_out << "Mesh.SecondOrderIncomplete = 0;\n";
	file_out << "Physical Surface(1) = {1} ;\n";
	file_out.close();
	
	NumberMap::GetGmshQuadrilateralNumbering(r, num_gmsh);
	MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(r, NumNodes2D, coor);
	CheckPermutation(num_gmsh);
	
	string command_line = string("geo2mesh quadrangle ") + to_str(r) + " > sortie";
	system(command_line.data());
	
	ifstream file_in("quadrangle.msh");
	string ligne;
	getline(file_in, ligne);
	getline(file_in, ligne);
	getline(file_in, ligne);
	getline(file_in, ligne);
	ReadParameters(file_in, num); int nb_vertices = num(0);
	Points2D.Reallocate(nb_vertices);
	for (int i = 0; i < nb_vertices; i++)
	  {
	    ReadParameters(file_in, coord);
	    Points2D(i).Init(coord(1), coord(2));
	  }
	
	getline(file_in, ligne);
	getline(file_in, ligne);
	ReadParameters(file_in, num); int nb_elt = num(0);
	for (int i = 0; i < nb_elt; i++)
	  {
	    ReadParameters(file_in, num);
	    if (num(0) == 5)
	      {
		for (int j = 0; j < num_gmsh.GetM(); j++)
		  {
		    int n = num(5+j);
		    int node = num_gmsh(j);
		    int i1 = coor(node, 0);
		    int i2 = coor(node, 1);
		    double x = step(i1);
		    double y = step(i2);
		    if (Points2D(n-1) != R2(x, y))
		      {
			DISP(r); DISP(j); DISP(node); DISP(i1); DISP(i2);
			DISP(x); DISP(y); DISP(Points2D(n-1));
			cout << "point not coincident" << endl;
			abort();
		      }
		  }
	      }
	  }
      }
    
  }
   
  
  Array3D<int> NumNodes3D;
  VectR3 Points3D;
  
  if (false)
  {
    cout << "Checking Tetrahedra" << endl;
    // checking high order tetrahedra
    for (int r = 2; r <= 10; r++)
      {
	Vector<double> step(r+1);
	step.Fill(); Mlt(1.0/r, step);
	
	ofstream file_out("tetrahedron.geo");
	file_out << "lc = 10.0;\n";
	file_out << "Point(1) = {0, 0, 0, lc};\n";
	file_out << "Point(2) = {1, 0, 0, lc};\n";
	file_out << "Point(3) = {0, 1, 0, lc};\n";
	file_out << "Point(4) = {0, 0, 1, lc};\n";
	file_out << "Line(1) = {1, 2};\n Line(2) = {2, 3};\n Line(3) = {3, 4};\n";
	file_out << "Line(4) = {4, 1};\n Line(5) = {2, 4};\n Line(6) = {1, 3};\n";
	file_out << "Transfinite Line {1,2,3,4, 5, 6} = 2;\n";
	file_out << "Line Loop(1) = {5, 4, 1} ;\n Plane Surface(1) = {1};\n";
	file_out << "Line Loop(2) = {3, 4, 6} ;\n Plane Surface(2) = {2};\n";
	file_out << "Line Loop(3) = {1, 2, -6} ;\n Plane Surface(3) = {3};\n";
	file_out << "Line Loop(4) = {3, -5, 2} ;\n Plane Surface(4) = {4};\n";
	file_out << "Surface Loop(1) = {4, 1, 2, 3};\n";
	file_out << "Volume(1) = {1};\n";
	file_out << "Physical Surface(1) = {1, 2, 3, 4} ;\n";
	file_out << "Mesh.SecondOrderIncomplete = 0;\n";
	file_out << "Physical Volume(1) = {1} ;\n";
	file_out.close();
	
	NumberMap::GetGmshTetrahedralNumbering(r, num_gmsh);
	MeshNumbering<Dimension3>::ConstructTetrahedralNumbering(r, NumNodes3D, coor);
	CheckPermutation(num_gmsh);
	
	string command_line = string("geo2mesh tetrahedron ") + to_str(r) + " > sortie";
	system(command_line.data());
	
	ifstream file_in("tetrahedron.msh");
	string ligne;
	getline(file_in, ligne);
	getline(file_in, ligne);
	getline(file_in, ligne);
	getline(file_in, ligne);
	ReadParameters(file_in, num); int nb_vertices = num(0);
	Points3D.Reallocate(nb_vertices);
	for (int i = 0; i < nb_vertices; i++)
	  {
	    ReadParameters(file_in, coord);
	    Points3D(i).Init(coord(1), coord(2), coord(3));
	  }
	
	getline(file_in, ligne);
	getline(file_in, ligne);
	ReadParameters(file_in, num); int nb_elt = num(0);
	R3 ptA, ptB, ptC, ptD;
	for (int i = 0; i < nb_elt; i++)
	  {
	    ReadParameters(file_in, num);
	    if (num(0) == 5)
	      {
		ptA = Points3D(num(5)-1);
		ptB = Points3D(num(6)-1);
		ptC = Points3D(num(7)-1);
		ptD = Points3D(num(8)-1);
		
		for (int j = 0; j < num_gmsh.GetM(); j++)
		  {
		    int n = num(5+j);
		    int node = num_gmsh(j);
		    int i0 = coor(node, 0);
		    int i1 = coor(node, 1);
		    int i2 = coor(node, 2);
		    int i3 = coor(node, 3);
		    double x = step(i1);
		    double y = step(i2);
		    double z = step(i3);
		    double l0 = step(i0);
		    R3 pt = l0 * ptA + x * ptB + y * ptC + z * ptD;
		    if (Points3D(n-1) != pt)
		      {
			DISP(r); DISP(j); DISP(node); DISP(i1); DISP(i2);
			DISP(i3); DISP(n); DISP(pt); DISP(Points3D(n-1));
			cout << "point not coincident" << endl;
			abort();
		      }
		  }
	      }
	  }
      }
    
  }   
  
  if (false)
  {
    cout << "Checking Hexahedra" << endl;
    // checking high order hexahedra
    for (int r = 2; r <= 9; r++)
      {
	Vector<double> step(r+1);
	step.Fill(); Mlt(1.0/r, step);
	
	ofstream file_out("hexahedron.geo");
	file_out << "lc = 10.0;\n";
	file_out << "Point(1) = {0, 0, 0, lc};\n";
	file_out << "Point(2) = {0, 0, 1, lc};\n";
	file_out << "Point(3) = {0, 1, 1, lc};\n";
	file_out << "Point(4) = {0, 1, 0, lc};\n";
	file_out << "Point(5) = {1, 0, 0, lc};\n";
	file_out << "Point(6) = {1, 0, 1, lc};\n";
	file_out << "Point(7) = {1, 1, 1, lc};\n";
	file_out << "Point(8) = {1, 1, 0, lc};\n";
	file_out << "Line(1) = {1, 2};\n Line(2) = {2, 3};\n Line(3) = {3, 4};\n";
	file_out << "Line(4) = {4, 1};\n Line(5) = {5, 6};\n Line(6) = {6, 7};\n";
	file_out << "Line(7) = {7, 8};\n Line(8) = {8, 5};\n Line(9) = {1, 5};\n";
	file_out << "Line(10) = {2, 6};\n Line(11) = {3, 7};\n Line(12) = {4, 8};\n";
	file_out << "Transfinite Line {1,2,3,4, 5, 6, 7, 8, 9, 10, 11, 12} = 2;\n";
	file_out << "Line Loop(1) = {1, 2, 3, 4} ;\n Plane Surface(1) = {1};\n";
	file_out << "Line Loop(2) = {5, 6, 7, 8} ;\n Plane Surface(2) = {2};\n";
	file_out << "Line Loop(3) = {1, 10, -5, -9} ;\n Plane Surface(3) = {3};\n";
	file_out << "Line Loop(4) = {2, 11, -6, -10} ;\n Plane Surface(4) = {4};\n";
	file_out << "Line Loop(5) = {3, 12, -7, -11} ;\n Plane Surface(5) = {5};\n";
	file_out << "Line Loop(6) = {4, 9, -8, -12} ;\n Plane Surface(6) = {6};\n";
	file_out << "Transfinite Surface {1} = {1, 2, 3, 4};\n";
	file_out << "Transfinite Surface {2} = {5, 6, 7, 8};\n";
	file_out << "Transfinite Surface {3} = {1, 2, 6, 5};\n";
	file_out << "Transfinite Surface {4} = {2, 3, 7, 6};\n";
	file_out << "Transfinite Surface {5} = {3, 4, 8, 7};\n";
	file_out << "Transfinite Surface {6} = {4, 1, 5, 8};\n";
	file_out << "Surface Loop(1) = {1, 2, 3, 4, 5, 6};\n";
	file_out << "Volume(1) = {1};\n";
	file_out << "Transfinite Volume {1} = {1, 2, 3, 4, 5, 6, 7, 8};\n";
	file_out << "Recombine Surface {1, 2, 3, 4, 5, 6} = 80;\n";
	file_out << "Physical Surface(1) = {1, 2, 3, 4, 5, 6} ;\n";
	file_out << "Mesh.SecondOrderIncomplete = 0;\n";
	file_out << "Physical Volume(1) = {1} ;\n";
	file_out.close();
	
	NumberMap::GetGmshHexahedralNumbering(r, num_gmsh);
	MeshNumbering<Dimension3>::ConstructHexahedralNumbering(r, NumNodes3D, coor);
	CheckPermutation(num_gmsh);
	
	string command_line = string("geo2mesh hexahedron ") + to_str(r) + " > sortie";
	system(command_line.data());
	
	ifstream file_in("hexahedron.msh");
	string ligne;
	getline(file_in, ligne);
	getline(file_in, ligne);
	getline(file_in, ligne);
	getline(file_in, ligne);
	ReadParameters(file_in, num); int nb_vertices = num(0);
	Points3D.Reallocate(nb_vertices);
	for (int i = 0; i < nb_vertices; i++)
	  {
	    ReadParameters(file_in, coord);
	    Points3D(i).Init(coord(1), coord(2), coord(3));
	  }
	
	getline(file_in, ligne);
	getline(file_in, ligne);
	ReadParameters(file_in, num); int nb_elt = num(0);
	R3 ptA, ptB, ptC, ptD, ptE, ptF, ptG, ptH;
	for (int i = 0; i < nb_elt; i++)
	  {
	    ReadParameters(file_in, num);
	    ptA = Points3D(num(5)-1); ptB = Points3D(num(6)-1);
	    ptC = Points3D(num(7)-1); ptD = Points3D(num(8)-1);
	    ptE = Points3D(num(9)-1); ptF = Points3D(num(10)-1);
	    ptG = Points3D(num(11)-1); ptH = Points3D(num(12)-1);
	    if (num(0) == 7)
	      {
		for (int j = 0; j < num_gmsh.GetM(); j++)
		  {
		    int n = num(5+j);
		    int node = num_gmsh(j);
		    int i1 = coor(node, 0);
		    int i2 = coor(node, 1);
		    int i3 = coor(node, 2);
		    double x = step(i1);
		    double y = step(i2);
		    double z = step(i3);
		    R3 pt = (1.0-z)*( (1.0-x)*(1.0-y)*ptA + x*(1.0-y)*ptB + x*y*ptC + (1.0-x)*y*ptD)
		      + z*( (1.0-x)*(1.0-y)*ptE + x*(1.0-y)*ptF + x*y*ptG + (1.0-x)*y*ptH);
		    if (Points3D(n-1) != pt)
		      {
			DISP(r); DISP(j); DISP(node); DISP(i1); DISP(i2);
			DISP(i3); DISP(n); DISP(pt); DISP(Points3D(n-1));
			cout << "point not coincident" << endl;
			abort();
		      }
		  }
	      }
	  }
      }
    
  }   

  //NumberMap::GetGmshTetrahedralNumbering(r, num);
  
  //CheckPermutation(num);
    
  if (argc > 1)
    {
      int r = atoi(argv[1]);
      Mesh<Dimension3> mesh;
      mesh.SetGeometryOrder(r);
      mesh.Read("tet_quad.msh");
      
      Vector<double> step(r+1);
      step.Fill(); Mlt(1.0/r, step);
      mesh.SubdivideMesh(step);
      
      mesh.Write("test.mesh"); 
    }
  
  return 0;
}
