#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  
  InitMontjoie(argc, argv);

  if ((argc == 2) && ((string(argv[1]) == "--help")
		      || (string(argv[1]) == "-help")))
    {
      cout << "Usage :" << endl;
      cout << "./sym_mesh.x name.mesh X" << endl;
      cout << "./sym_mesh.x name.mesh XY" << endl;
      cout << "./sym_mesh.x name.mesh XYZ" << endl;
      cout << "./sym_mesh.x name.mesh XYZ ordre" << endl;
      cout << "./sym_mesh.x name.mesh Plane ordre a b c d" << endl;

      cout << endl << "The first argument is the input and output mesh. The mesh is overwritten by its symmetrization. " << endl;
      cout << "The second argument describes the type of symmetry. It can be a combination of axis X, Y and Z (X, Y, Z, XZ, XY, YX or XYZ) or Plane" << endl;
      cout << "The optional third argument is the order of approximation for the mesh" << endl;
      cout << "If the second argument is Plane, we expect the constants a, b, c (and d in 3-D) " << endl;
      cout << "such that a x + b y + c = 0 (a x + b y + c z + d = 0 in 3-D) is the equation of the plane of symmetry" << endl;

      cout << endl << "Examples" << endl << endl;
      cout << "./sym_mesh.x toto.mesh Y" << endl
	   << "         produces the mesh obtained by symmetry with respect to plane y=0" << endl << endl;
      
      cout << "./sym_mesh.x toto.mesh XZ" << endl
	   << "         produces the mesh obtained by symmetry with respect to plane x=0, followed by symmetry with respect to plane z=0" << endl;

      cout << "./sym_mesh.x toto.mesh YZ 2" << endl
	   << "         produces the mesh obtained by symmetry with respect to plane y=0, followed by symmetry with respect to plane z=0." << endl
	   << "         The produced mesh is a curved mesh of degree 2" << endl << endl;
      
      cout << "./sym_mesh.x toto.mesh Plane 1 1.0 0.5 2.0 " << endl
	   << "         produces the mesh obtained by symmetry with respect to plane x + y/2 + 2 = 0 (2-D mesh) " << endl << endl;
      
      exit(0);
    }

  if (argc < 3)
    {
      cout<<"This command requires at least two arguments "<<endl;
      cout<<"sym_mesh.x input.mesh X"<<endl;
      abort();
    }
  
  int r = 1;
  if (argc >= 4)
    r = to_num<int>(argv[3]);
  
  string input_name(argv[1]), symmetry_info(argv[2]);
  int dim_N = GetDimensionMesh(input_name, GetExtension(input_name));
  
  if (dim_N == 2)
    {
      Mesh<Dimension2> mesh;
      R2 normale;
      
      mesh.SetGeometryOrder(r);
      mesh.Read(input_name);

      if (symmetry_info == "Plane")
        {
          if (argc < 7)
            {
              cout << "Enter parameters of the plane " << endl;
              abort();
            }
          
          normale(0) = to_num<Real_wp>(argv[4]);
          normale(1) = to_num<Real_wp>(argv[5]);
          Real_wp cte = to_num<Real_wp>(argv[6]);
          mesh.SymmetrizeMesh(-1, cte, normale);
        }
      else
        {      
          if (symmetry_info.find("X") != string::npos)
            mesh.SymmetrizeMesh(0, 0.0, normale);
          
          if (symmetry_info.find("Y") != string::npos)
            mesh.SymmetrizeMesh(1, 0.0, normale);
        }
      
      mesh.Write(input_name);
    }
  else if (dim_N == 3)
    {
      Mesh<Dimension3> mesh;
      R3 normale;
      
      mesh.SetGeometryOrder(r);
      mesh.Read(input_name);
      
      if (symmetry_info == "Plane")
        {
          if (argc < 8)
            {
              cout << "Enter parameters of the plane " << endl;
              abort();
            }
          
          normale(0) = to_num<Real_wp>(argv[4]);
          normale(1) = to_num<Real_wp>(argv[5]);
          normale(2) = to_num<Real_wp>(argv[6]);
          Real_wp cte = to_num<Real_wp>(argv[7]);
          mesh.SymmetrizeMesh(-1, cte, normale);
        }
      else
        {
          if (symmetry_info.find("X") != string::npos)
            mesh.SymmetrizeMesh(0, 0.0, normale);
          
          if (symmetry_info.find("Y") != string::npos)
            mesh.SymmetrizeMesh(1, 0.0, normale);
          
          if (symmetry_info.find("Z") != string::npos)
            mesh.SymmetrizeMesh(2, 0.0, normale);
        }
      
      mesh.Write(input_name);
    }
  
  return FinalizeMontjoie();
  
}
