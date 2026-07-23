#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc != 5)
    {
      cout<<"Cette commande demande trois arguments"<<endl;
      cout<<"reg.x cube.geo nbx nby nbz"<<endl;
      cout<<"the file cube.geo is created"<<endl;
    }
    
  /* Mesh<Dimension3> mesh;

  int nbx(2), nby(2), nbz(2), ref_domain = 1;
  int ref = 1;
  TinyVector<int, 6> ref_boundary; ref_boundary = ref; 

  mesh.CreateRegularMesh(R3(-1,-1,-1), R3(1,1,1), TinyVector<int, 3>(nbx, nby, nbz),
			 ref_domain, ref_boundary, mesh.TETRAHEDRAL_MESH);
  
  mesh.Write("tet.mesh");
  
  DISP(mesh.GetNbEdges());
  DISP(mesh.GetNbEdgesRef());
  for (int i = 0; i < mesh.GetNbEdges(); i++)
    {
      DISP(i); DISP(mesh.GetEdge(i));
    }
  
  mesh.SplitIntoHexahedra();
  
  mesh.Write("hex.mesh");
  
  return 0;
  */

  Mesh<Dimension3> mesh;
  int nbx(2), nby(2), nbz(2), ref_domain = 1;
  int ref = 1;
  TinyVector<int, 6> ref_boundary; ref_boundary.Fill(ref); 

  string file_name(argv[1]);
  to_num(argv[2], nbx);
  to_num(argv[3], nby);
  to_num(argv[4], nbz);
  
  mesh.CreateRegularMesh(R3(-1,-1,-1), R3(1,1,1), TinyVector<int, 3>(nbx, nby, nbz),
			 ref_domain, ref_boundary, mesh.HEXAHEDRAL_MESH);
  
  ofstream file_out(file_name.data());
  file_out<<"lc = 0.1; "<<endl<<endl;
  
  file_out<<endl;
  for (int i = 0; i < nbx; i++)
    file_out << "x" << i << " = " << Real_wp(i)/(nbx-1) << ";\n";
  
  file_out << endl;
  for (int i = 0; i < nby; i++)
    file_out << "y" << i << " = " << Real_wp(i)/(nby-1) << ";\n";

  file_out << endl;
  for (int i = 0; i < nbz; i++)
    file_out << "z" << i << " = " << Real_wp(i)/(nbz-1) << ";\n";

  file_out << endl;
  for (int i = 0; i < nbx-1; i++)
    file_out << "Nx" << i+1 << " = Ceil((x" << i+1 << " - x" << i <<")/lc) + 1;\n";
  
  file_out << endl;
  for (int i = 0; i < nby-1; i++)
    file_out << "Ny" << i+1 << " = Ceil((y" << i+1 << " - y" << i <<")/lc) + 1;\n";
  
  file_out << endl;
  for (int i = 0; i < nbz-1; i++)
    file_out << "Nz" << i+1 << " = Ceil((z" << i+1 << " - z" << i <<")/lc) + 1;\n";
  
  file_out << endl;
  int numv = 1;
  for (int k = 0; k < nbz; k++)
    for (int j = 0; j < nby; j++)
      for (int i = 0; i < nbx; i++)
        file_out<<"Point("<<(numv++)<<") = {x"<<i<<", y"<<j<<", z"<<k<<", lc};"<<endl;

  file_out << endl;
  for (int i = 0; i < mesh.GetNbEdges(); i++)
    file_out<<"Line("<<(i+1)<<") = {"<<mesh.GetEdge(i).numVertex(0)+1<<","<<mesh.GetEdge(i).numVertex(1)+1<<"};"<<endl;
  
  Real_wp dx = (mesh.GetXmax()-mesh.GetXmin())/(nbx-1);
  Real_wp dy = (mesh.GetYmax()-mesh.GetYmin())/(nby-1);
  Real_wp dz = (mesh.GetZmax()-mesh.GetZmin())/(nbz-1);
  R3 ptA, ptB, vec_u;
  file_out<<endl;
  for (int i = 0; i < (nbx-1); i++)
    {
      file_out<<"Transfinite Line{";
      int nb = 0;
      for (int k = 0; k < mesh.GetNbEdges(); k++)
	{
	  int ne0 = mesh.GetEdge(k).numVertex(0);
	  int ne1 = mesh.GetEdge(k).numVertex(1);
	  ptA = mesh.Vertex(ne0);
	  ptB = mesh.Vertex(ne1);
	  vec_u = ptA - ptB;
	  if (abs(vec_u(0)) > 1e-5)
	    {
	      int j = toInteger(floor((0.5*ptA(0) + 0.5*ptB(0) - mesh.GetXmin())/dx));
	      if (j == i)
		{
		  if (nb == 0)
		    file_out<<k+1;
		  else
		    file_out<<","<<k+1;
		  
		  nb++;
		}
	    }
	}
      
      file_out<<"} = Nx"<<i+1<<";"<<endl;
    }

  file_out<<endl;  
  for (int i = 0; i < (nby-1); i++)
    {
      file_out<<"Transfinite Line{";
      int nb = 0;
      for (int k = 0; k < mesh.GetNbEdges(); k++)
	{
	  int ne0 = mesh.GetEdge(k).numVertex(0);
	  int ne1 = mesh.GetEdge(k).numVertex(1);
	  ptA = mesh.Vertex(ne0);
	  ptB = mesh.Vertex(ne1);
	  vec_u = ptA - ptB;
	  if (abs(vec_u(1)) > 1e-5)
	    {
	      int j = toInteger(floor((0.5*ptA(1) + 0.5*ptB(1) - mesh.GetYmin())/dy));
	      if (j == i)
		{
		  if (nb == 0)
		    file_out<<k+1;
		  else
		    file_out<<","<<k+1;
		  
		  nb++;
		}
	    }
	}
      
      file_out<<"} = Ny"<<i+1<<";"<<endl;
    }      
  
  file_out << endl;
  for (int i = 0; i < (nbz-1); i++)
    { 
      file_out<<"Transfinite Line{";
      int nb = 0;
      for (int k = 0; k < mesh.GetNbEdges(); k++)
	{
	  int ne0 = mesh.GetEdge(k).numVertex(0);
	  int ne1 = mesh.GetEdge(k).numVertex(1);
	  ptA = mesh.Vertex(ne0);
	  ptB = mesh.Vertex(ne1);
	  vec_u = ptA - ptB;
	  if (abs(vec_u(2)) > 1e-5)
	    {
	      int j = toInteger(floor((0.5*ptA(2) + 0.5*ptB(2) - mesh.GetZmin())/dz));
	      if (j == i)
		{
		  if (nb == 0)
		    file_out<<k+1;
		  else
		    file_out<<","<<k+1;
		  
		  nb++;
		}
	    }
	}
      
      file_out<<"} = Nz"<<i+1<<";"<<endl;
    }

  file_out<<endl;
  for (int i = 0; i < mesh.GetNbBoundary(); i++)
    {
      file_out<<"Line Loop("<<i+1<<") = {";
      for (int j = 0; j < 4; j++)
	{
	  int nv0 = mesh.Boundary(i).numVertex(j)+1;
	  // int nv1 = mesh.faces(i).numVertex((j+1)%4)+1;
	  int num_edge = mesh.Boundary(i).numEdge(j);
	  int ne0 = mesh.GetEdge(num_edge).numVertex(0)+1;
	  // int ne1 = mesh.(num_edge).numVertex(1)+1;
	  if (nv0 == ne0)
	    file_out<<num_edge+1;
	  else
	    file_out<<-num_edge-1;
	  
	  if (j < 3)
	    file_out<<",";
	}
      file_out<<"};"<<endl;
    }
  
  file_out<<endl;
  for (int i = 0; i < mesh.GetNbBoundary(); i++)
    file_out<<"Plane Surface("<<i+1<<") = {"<<i+1<<"};"<<endl;
  
  file_out<<endl;
  for (int i = 0; i < mesh.GetNbBoundary(); i++)
    {
      file_out<<"Transfinite Surface{"<<i+1<<"} = {"<<mesh.Boundary(i).numVertex(0)+1<<","<<mesh.Boundary(i).numVertex(1)+1;
      file_out<<","<<mesh.Boundary(i).numVertex(2)+1<<","<<mesh.Boundary(i).numVertex(3)+1<<"};"<<endl;
    }
  
  file_out<<endl;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      file_out<<"Surface Loop("<<i+1<<") = {";
      for (int j = 0; j < 6; j++)
	{
	  file_out<<mesh.Element(i).numBoundary(j)+1;
	  if (j < 5)
	    file_out<<",";
	}
      file_out<<"};"<<endl;
    }
  
  file_out<<endl;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    file_out<<"Volume("<<i+1<<") = {"<<i+1<<"};"<<endl;
  
  file_out<<endl;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      file_out<<"Transfinite Volume{"<<i+1<<"} = {";
      for (int j = 0; j < 8; j++)
	{
	  file_out<<mesh.Element(i).numVertex(j)+1;
	  if (j < 7)
	    file_out<<",";
	}
      file_out<<"};"<<endl;
    }
  
  file_out<<endl;
  file_out<<"Recombine Surface{1:"<<mesh.GetNbBoundary()<<"} = 200;"<<endl;
  
  file_out.close();
  
  return 0;
}  
