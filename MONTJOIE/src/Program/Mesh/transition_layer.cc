#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

void CreateTransitionLayer(Mesh<Dimension3>& mesh)
{
  Real_wp h = 0.3;
  Mesh<Dimension3> mesh2(mesh);
  VectReal_wp step(3); step(0) = 0.0; step(1) = 0.5; step(2) = 1.0;
  mesh2.SubdivideMesh(step);
  
  TranslateMesh(mesh2, R3(0, 0, -h));
  mesh2.Write("refined.mesh");

  // on trie tous les sommets du maillage grossier
  VectR3 vert_coarse(mesh.Vertex());
  IVect permut(vert_coarse.GetM()); permut.Fill();
  Sort(vert_coarse, permut);

  Vector<bool> black_point(vert_coarse.GetM());
  black_point.Fill(false);
  {
    Real_wp dist_prec = 0; int i = 1, j = 0;
    for (int n = 1; n < vert_coarse.GetM(); n++)
      {
	Real_wp dist = vert_coarse(i).Distance(vert_coarse(i-1));
	if (n > 1)
	  {
	    if (dist > 1.8*dist_prec)
	      {
		i = 0;
		j++;
	      }
	  }

	if ((i%2 == 1) && (j%2 == 1))
	  black_point(permut(n)) = true;
	
	dist_prec = dist; i++;
      }
  }

  DISP(black_point);
  
  // on rajoute les sommets
  int N = mesh.GetNbBoundaryRef();
  //int Ne = mesh.GetNbEdges();
  int Ns = mesh.GetNbVertices(), Ns2 = mesh2.GetNbVertices();
  mesh.ResizeVertices(Ns + Ns2 + 10*N);
  for (int i = 0; i < Ns2; i++)
    mesh.Vertex(Ns+i) = mesh2.Vertex(i);

  int off_int = Ns + Ns2;
  //mesh.AppendMesh(mesh2);
  mesh.ResizeBoundariesRef(5*N);
  mesh.ReallocateElements(7*N);
  std::cout<<N<<endl;
  for (int i = 0; i < N; i++)
    {
      int n0 = mesh.BoundaryRef(i).numVertex(0);
      int n1 = mesh.BoundaryRef(i).numVertex(1);
      int n2 = mesh.BoundaryRef(i).numVertex(2);
      int n3 = mesh.BoundaryRef(i).numVertex(3);

      VectR3 A(9);
      A(0) = mesh.Vertex(n0); A(1) = mesh.Vertex(n1);
      A(2) = mesh.Vertex(n2); A(3) = mesh.Vertex(n3);
      
      A(4) = 0.5*(A(0) + A(1)); A(5) = 0.5*(A(1) + A(2));
      A(6) = 0.5*(A(2) + A(3)); A(7) = 0.5*(A(0) + A(3));
      A(8) = 0.25*(A(0) + A(1) + A(2) + A(3));

      IVect num(9); num.Fill(-1);
      
      for (int j = 0; j < 4; j++)
	for (int k = 0; k < 4; k++)
	  {
	    int nB = mesh2.BoundaryRef(4*i+j).numVertex(k);
	    R3 B = mesh2.Vertex(nB); B(2) += h;
	   
	    for (int p = 0; p < A.GetM(); p++)
	      if (A(p) == B)
		num(p) = nB;
	        
	  }
      
      int p0 = num(0) + Ns, p1 = num(1) + Ns, p2 = num(2) + Ns;
      int p3 = num(3) + Ns, p4 = num(4) + Ns;
      int p5 = num(5) + Ns, p6 = num(6) + Ns, p7 = num(7) + Ns, p8 = num(8) + Ns;
      
      mesh.BoundaryRef(N+i).InitQuadrangular(p0, p4, p8, p7, 1);
      mesh.BoundaryRef(2*N+i).InitQuadrangular(p4, p1, p5, p8, 1);
      mesh.BoundaryRef(3*N+i).InitQuadrangular(p7, p8, p6, p3, 1);
      mesh.BoundaryRef(4*N+i).InitQuadrangular(p8, p5, p2, p6, 1);
      //int m=int(sqrt(N));
      
      //1er cas
      
      //if (i%8==0 or i%8==2)
      if (black_point(n2))
	{
	  // points du plan -0.25
	  VectR3 B(4);
      
	  B(0)=A(2); B(0)(2) -= 0.25*h;
	  B(1)=A(1); B(1)(2) -= 0.25*h;
	  B(2)=A(4); B(2)(2) -= 0.25*h;
	  B(3)=A(6); B(3)(2) -= 0.25*h;
	  
	  for (int l=0; l<4;l++)
	    mesh.Vertex(off_int+l) = B(l);
	  
	  // points du plan -0.75
	  VectR3 D(4);
	  
	  D(0)=A(2); D(0)(2) -= 0.75*h;
	  D(1)=A(5); D(1)(2) -= 0.75*h;
	  D(2)=A(8); D(2)(2) -= 0.75*h;
	  D(3)=A(6); D(3)(2) -= 0.75*h;
	  
	  for (int l=4; l<8;l++)
	    mesh.Vertex(off_int+l) = D(l-4);
	  
	  // points du plan -0.5
	  VectR3 C(2);
	  
	  C(0)=A(3); C(0)(2) -= 0.5*h;
	  C(1)=A(7); C(1)(2) -= 0.5*h;
	  
	  for (int l=8; l<10;l++)
	    mesh.Vertex(off_int+l) = C(l-8);
	  
	  
	  
	  mesh.Element(i).InitHexahedral( off_int,off_int+1, off_int+2, off_int+3,n2,n1,n0,n3, 1);
     
	  mesh.Element(N+i).InitHexahedral( off_int+4,off_int+5, off_int+6, off_int+7,off_int,off_int+1, off_int+2, off_int+3, 1);
   
	  mesh.Element(2*N+i).InitHexahedral( off_int+4,off_int+5, off_int+6, off_int+7,p2,p5,p8,p6, 1);
      
	  mesh.Element(3*N+i).InitHexahedral( off_int+5,off_int+1, off_int+2, off_int+6,p5,p1,p4,p8, 1);
        
	  mesh.Element(4*N+i).InitHexahedral( off_int+7,off_int+6, off_int+9, off_int+8,p6,p8,p7,p3, 1);
      
	  mesh.Element(5*N+i).InitHexahedral( off_int+7,off_int+6, off_int+9, off_int+8,off_int+3,off_int+2,n0,n3, 1);

	  mesh.Element(6*N+i).InitHexahedral( off_int+6,off_int+2, n0, off_int+9,p8,p4,p0,p7, 1);
      
	  off_int += 10;
	}
      
      //2eme cas
      //if (i%8==4 or i%8==6)
      if (black_point(n1))
	{
	  // points du plan -0.25
	  VectR3 B(4);
      
	  B(0)=A(2); B(0)(2) -= 0.25*h;
	  B(1)=A(1); B(1)(2) -= 0.25*h;
	  B(2)=A(4); B(2)(2) -= 0.25*h;
	  B(3)=A(6); B(3)(2) -= 0.25*h;
	
	  for (int l=0; l<4;l++)
	    mesh.Vertex(off_int+l) = B(l);
	
	  // points du plan -0.75
	  VectR3 D(4);
	
	  D(0)=A(5); D(0)(2) -= 0.75*h;
	  D(1)=A(1); D(1)(2) -= 0.75*h;
	  D(2)=A(4); D(2)(2) -= 0.75*h;
	  D(3)=A(8); D(3)(2) -= 0.75*h;
	
	  for (int l=4; l<8;l++)
	    mesh.Vertex(off_int+l) = D(l-4);
	
	  // points du plan -0.5
	  VectR3 C(2);
	
	  C(0)=A(0); C(0)(2) -= 0.5*h;
	  C(1)=A(7); C(1)(2) -= 0.5*h;
	
	  for (int l=8; l<10;l++)
	    mesh.Vertex(off_int+l) = C(l-8);
	
      
      
	  mesh.Element(i).InitHexahedral( off_int,off_int+1, off_int+2, off_int+3,n2,n1,n0,n3, 1);
     
	  mesh.Element(N+i).InitHexahedral( off_int,off_int+4, off_int+7, off_int+3,p2,p5,p8,p6, 1);
   
	  mesh.Element(2*N+i).InitHexahedral( off_int+4,off_int+5, off_int+6, off_int+7,p5,p1,p4,p8, 1);
      
	  mesh.Element(3*N+i).InitHexahedral( off_int,off_int+1, off_int+2, off_int+3,off_int+4,off_int+5, off_int+6, off_int+7, 1);
        
	  mesh.Element(4*N+i).InitHexahedral( off_int+2,off_int+3, n3,n0,off_int+6,off_int+7, off_int+9, off_int+8, 1);
      
	  mesh.Element(5*N+i).InitHexahedral( off_int+6,off_int+7, off_int+9, off_int+8,p4,p8,p7,p0, 1);

	  mesh.Element(6*N+i).InitHexahedral( off_int+7,off_int+9, n3, off_int+3,p8,p7,p3,p6, 1);
      
	  off_int += 10;
	}

      //3eme cas
      //if (i%8==5 or i%8==7)
      if (black_point(n0))
	{
	  // points du plan -0.25
	  VectR3 B(4);
      
	  B(0)=A(0); B(0)(2) -= 0.25*h;
	  B(1)=A(3); B(1)(2) -= 0.25*h;
	  B(2)=A(6); B(2)(2) -= 0.25*h;
	  B(3)=A(4); B(3)(2) -= 0.25*h;
	
	  for (int l=0; l<4;l++)
	    mesh.Vertex(off_int+l) = B(l);
	
	  // points du plan -0.75
	  VectR3 D(4);
	
	  D(0)=A(7); D(0)(2) -= 0.75*h;
	  D(1)=A(0); D(1)(2) -= 0.75*h;
	  D(2)=A(4); D(2)(2) -= 0.75*h;
	  D(3)=A(8); D(3)(2) -= 0.75*h;
	
	  for (int l=4; l<8;l++)
	    mesh.Vertex(off_int+l) = D(l-4);
	
	  // points du plan -0.5
	  VectR3 C(2);
	
	  C(0)=A(1); C(0)(2) -= 0.5*h;
	  C(1)=A(5); C(1)(2) -= 0.5*h;
	
	  for (int l=8; l<10;l++)
	    mesh.Vertex(off_int+l) = C(l-8);
	
      
      
	  mesh.Element(i).InitHexahedral( off_int,off_int+1, off_int+2, off_int+3,n0,n3,n2,n1, 1);
     
	  mesh.Element(N+i).InitHexahedral( off_int+1,off_int+2, off_int+7, off_int+4,p3,p6,p8,p7, 1);
   
	  mesh.Element(2*N+i).InitHexahedral( off_int+4,off_int+7, off_int+6, off_int+5,p7,p8,p4,p0, 1);
      
	  mesh.Element(3*N+i).InitHexahedral( off_int+4,off_int+7, off_int+6, off_int+5,off_int+1,off_int+2, off_int+3, off_int, 1);
        
	  mesh.Element(4*N+i).InitHexahedral( off_int+2,n2,off_int+9,off_int+7,p6,p2,p5,p8, 1);
      
	  mesh.Element(5*N+i).InitHexahedral( off_int+2,n2,n1, off_int+3,off_int+7,off_int+9,off_int+8,off_int+6, 1);

	  mesh.Element(6*N+i).InitHexahedral( off_int+9,off_int+8,off_int+6,off_int+7,p5,p1,p4,p8, 1);
      
	  off_int += 10;
	}

      //4eme cas
      //if (i%8==1 or i%8==3)
      if (black_point(n3))
	{
	  // points du plan -0.25
	  VectR3 B(4);
      
	  B(0)=A(0); B(0)(2) -= 0.25*h;
	  B(1)=A(3); B(1)(2) -= 0.25*h;
	  B(2)=A(6); B(2)(2) -= 0.25*h;
	  B(3)=A(4); B(3)(2) -= 0.25*h;
	
	  for (int l=0; l<4;l++)
	    mesh.Vertex(off_int+l) = B(l);
	
	  // points du plan -0.75
	  VectR3 D(4);
	
	  D(0)=A(3); D(0)(2) -= 0.75*h;
	  D(1)=A(7); D(1)(2) -= 0.75*h;
	  D(2)=A(8); D(2)(2) -= 0.75*h;
	  D(3)=A(6); D(3)(2) -= 0.75*h;
	
	  for (int l=4; l<8;l++)
	    mesh.Vertex(off_int+l) = D(l-4);
	
	  // points du plan -0.5
	  VectR3 C(2);
	
	  C(0)=A(2); C(0)(2) -= 0.5*h;
	  C(1)=A(5); C(1)(2) -= 0.5*h;
	
	  for (int l=8; l<10;l++)
	    mesh.Vertex(off_int+l) = C(l-8);
	
      
      
	  mesh.Element(i).InitHexahedral( off_int,off_int+1, off_int+2, off_int+3,n0,n3,n2,n1, 1);
     
	  mesh.Element(N+i).InitHexahedral( off_int+4,off_int+7, off_int+6, off_int+5,p3,p6,p8,p7, 1);
   
	  mesh.Element(2*N+i).InitHexahedral( off_int+5,off_int+6, off_int+3, off_int,p7,p8,p4,p0, 1);
      
	  mesh.Element(3*N+i).InitHexahedral( off_int+4,off_int+7, off_int+6, off_int+5,off_int+1,off_int+2, off_int+3, off_int, 1);
        
	  mesh.Element(4*N+i).InitHexahedral( off_int+7,off_int+8,off_int+9,off_int+6,p6,p2,p5,p8, 1);
      
	  mesh.Element(5*N+i).InitHexahedral( off_int+7,off_int+8,off_int+9,off_int+6,off_int+2,n2,n1,off_int+3, 1);

	  mesh.Element(6*N+i).InitHexahedral( off_int+9,n1,off_int+3,off_int+6,p5,p1,p4,p8, 1);
      
	  off_int += 10;
	}
    }
  
  mesh.FindConnectivity();
  mesh.RemoveDuplicateVertices();
  mesh.AddBoundaryFaces();
  mesh.FindConnectivity();
}

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc < 3)
    {
      cout<<"Cette commande demande au moins trois arguments" << endl;
      cout<<"transition_layer surface.mesh final.mesh" << endl;
      cout<<"the file final.mesh is created" << endl;
      return -1;
    }

  DISP(R3::threshold);
  string file_input(argv[1]);
  string file_output(argv[2]);

  Mesh<Dimension3> mesh;
  mesh.Read(file_input);

  CreateTransitionLayer(mesh);
  
  mesh.Write(file_output);
}  
  
