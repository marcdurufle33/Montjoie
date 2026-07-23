#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc < 9)
    {
      cout<<"Cette commande demande au moins six arguments"<<endl;
      cout<<"create_wire N R alpha zmax delta toto.mesh Nc Nt"<<endl;
      cout<<"the file toto.mesh is created"<<endl;
      return -1;
    }
  
  int N = atoi(argv[1]);
  double R = atof(argv[2]);
  double alpha = pi_wp/180.0*atof(argv[3]);
  double zmax = atof(argv[4]);
  //double teta_max = pi_wp/180.0*atof(argv[4]);
  double delta = atof(argv[5]);
  string file_output(argv[6]);
  
  Mesh<Dimension3> mesh;
  
  // valeur maximale de z 
  //double zmax = teta_max*R/tan(alpha);
  double teta_max = zmax*tan(alpha)/R;
  
  R3 ex(1, 0, 0), ey(0, 1, 0), ez(0, 0, 1), center(0, 0, 0);
  //center.Init(2.24, -1.8110442, 3.0838817);
  ex.Init(0, 1, 0); ey.Init(0, 0, 1); ez.Init(1, 0, 0);
  int nb_subdiv_circle = atoi(argv[7]);
  int nb_subdiv_teta = atoi(argv[8]);
  int nb_points = 2*N + nb_subdiv_circle*(nb_subdiv_teta+1)*N;
  int nb_elt = 2*nb_subdiv_circle*N + nb_subdiv_circle*nb_subdiv_teta*N;
  int Nteta = nb_subdiv_circle*N;
  mesh.ReallocateVertices(nb_points);
  mesh.ReallocateBoundariesRef(nb_elt);
  // loop over teta-subdivisions of a wire
  R3 ptA, ptB, ux, uy, uplan, ptC; int ref = 2;
  nb_points = 0; nb_elt = 0;
  VectReal_wp cos_tet(N), sin_tet(N), tet(N);
  Vector<Matrix3_3> RotTet(N); 
  TinyVector<R3, 3> axis;
  axis(0) = ex; axis(1) = ey; axis(2) = ez;
  Matrix3_3 matQ(axis), matRot, A;  
  for (int i = 0; i < N; i++)
    {
      tet(i) = 2.0*pi_wp*Real_wp(i)/N;
      cos_tet(i) = cos(tet(i));
      sin_tet(i) = sin(tet(i));
      matRot.Fill(0);
      matRot(0, 0) = cos_tet(i); matRot(0, 1) = -sin_tet(i);
      matRot(1, 0) = sin_tet(i); matRot(1, 1) = cos_tet(i);
      matRot(2, 2) = 1.0;
      MltTrans(matRot, matQ, A);
      Mlt(matQ, A, RotTet(i));
    }
  
  for (int nt = 0; nt <= nb_subdiv_teta; nt++)
    {
      Real_wp teta = Real_wp(nt)*teta_max/nb_subdiv_teta;
      Real_wp z = Real_wp(nt)*zmax/nb_subdiv_teta;
      // center of the small section
      ptA = R*cos(teta)*ex + R*sin(teta)*ey + z*ez;
      
      // axis on this section
      ux = cos(teta)*ex + sin(teta)*ey;
      uplan = sin(teta)*ex - cos(teta)*ey;
      uy = cos(alpha)*uplan + sin(alpha)*ez;
      
      // section pour z = 0
      if (nt == 0)
        {
          for (int i = 0; i < N; i++)
            {
              Mlt(RotTet(i), ptA, ptB);
              mesh.Vertex(nb_points+i) = center + ptB;
            }
          
          for (int i = 0; i < nb_subdiv_circle; i++)
            for (int j = 0; j < N; j++)
              {
                int ip2 = i+2;
                if (i == nb_subdiv_circle-1)
                  ip2 = 1;
                
                mesh.BoundaryRef(nb_elt).InitTriangular(nb_points+j, nb_points + (i+1)*N + j, nb_points + ip2*N + j, ref);
                nb_elt++;
              }
          
          nb_points += N; 
        }
      
      for (int i = 0; i < nb_subdiv_circle; i++)
        {
          Real_wp phi = 2.0*pi_wp*Real_wp(i)/nb_subdiv_circle;
          ptB = ptA + cos(phi)*delta*ux + sin(phi)*delta*uy;
          for (int j = 0; j < N; j++)
            {
              Mlt(RotTet(j), ptB, ptC);
              mesh.Vertex(nb_points + i*N + j) = center + ptC;
              if (nt < nb_subdiv_teta)
                {
                  int ip1 = i+1;
                  if (i == nb_subdiv_circle-1)
                    ip1 = 0;
                  
                  mesh.BoundaryRef(nb_elt).InitQuadrangular(nb_points + i*N + j, nb_points + ip1*N + j,
                                                            nb_points + ip1*N + j + Nteta, nb_points + i*N + j + Nteta, ref);
                  
                  nb_elt++;
                }
            }
        }
      
      nb_points += N*nb_subdiv_circle;

      // section z = zmax
      if (nt == nb_subdiv_teta)
        {
          for (int i = 0; i < N; i++)
            {
              Mlt(RotTet(i), ptA, ptB);
              mesh.Vertex(nb_points+i) = center + ptB;
            }
          
          for (int i = 0; i < nb_subdiv_circle; i++)
            for (int j = 0; j < N; j++)
              {
                int ip2 = i+2;
                if (i == nb_subdiv_circle-1)
                  ip2 = 1;
                
                mesh.BoundaryRef(nb_elt).InitTriangular(nb_points+j, nb_points - (i+1)*N + j, nb_points - ip2*N + j, ref);
                nb_elt++;
              }
          
          nb_points += N; 
        }
        

    }
  
  //DISP(nb_points); DISP(mesh.GetNbVertices());
  //DISP(nb_elt); DISP(mesh.GetNbBoundaryRef());
  mesh.FindConnectivity();
  
  // on ecrit le maillage
  mesh.Write(file_output);

  return FinalizeMontjoie();
}  

