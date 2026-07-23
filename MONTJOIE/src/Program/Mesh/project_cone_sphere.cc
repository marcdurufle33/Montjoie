#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

class ProjectorConeSphere
{
public :
  double L, R, teta, ri, sin_teta, yi;
  R3 center, axis;
  
  void Init(double teta_, double R_, double ri_)
  {
    teta = teta_;
    R = R_;
    ri = ri_;
    L = R/sin(teta);
    sin_teta = sin(teta);
    yi = L - ri/sin_teta;    
    center.Init(0, L, 0);
    axis.Init(0, 1, 0);
  }
  
  void Project(R3& ptA)
  {
    double x = ptA(0);
    double y = ptA(1);
    double z = ptA(2);
    
    // on repere l'angle phi de rotation autour de l'axe Oy
    double radius = sqrt(x*x + z*z);
    if (radius <= epsilon_machine)
      {
	// on est sur l'axe de rotation Oy
	ptA(0) = 0;
	ptA(2) = 0;
	if (y < 0)
	  ptA(1) = -R;
	else
	  ptA(1) = yi+ri;
	
	return;
      }
    
    // sur la section, le point considere est (radius, y)
    // on regarde l'angle avec la ligne Ox
    double rp = sqrt(radius*radius + y*y);
    double alpha = 0;
    if (radius >= rp)
      alpha = 0;
    else
      alpha = acos(radius/rp);
    
    if (y < 0)
      alpha = -alpha;
    
    if (alpha < teta)
      {
	// on est sur la grosse sphere
	Mlt(R/Norm2(ptA), ptA);
	
	return;
      }
    
    // angle avec la ligne y = yi
    rp = sqrt(radius*radius + (y-yi)*(y-yi));
    alpha = 0;
    if (radius >= rp)
      alpha = 0;
    else
      alpha = acos(radius/rp);
    
    if (y < yi)
      alpha = -alpha;
    
    if (alpha > teta)
      {
	// on est sur la petite sphere
        ptA(1) -= yi;
	Mlt(ri/Norm2(ptA), ptA);
	ptA(1) += yi;
	
	return;
      }
    
    // on est sur le cone
    Mesh<Dimension3>::ProjectToConicCylinder(ptA, center, axis, teta);

  }
  
};

void ProjectOnConeSphere(Mesh<Dimension3>& mesh, const IVect& ref_cond, 
			 double teta, double R, double ri)
{
  ProjectorConeSphere proj;
  proj.Init(teta, R, ri);
  
  mesh.SetCurveType(1, mesh.CURVE_FILE);
  
  Vector<bool> EdgeUsed(mesh.GetNbEdges()), VertexUsed(mesh.GetNbVertices());
  EdgeUsed.Fill(false); VertexUsed.Fill(false);
  int nb_edges_ref = 0;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      if ((ref > 0) && (ref_cond(ref) == 1))
	{
	  for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
            VertexUsed(mesh.BoundaryRef(i).numVertex(j)) = true;
          
	  for (int j = 0; j < mesh.BoundaryRef(i).GetNbEdges(); j++)
	    {
	      int ne = mesh.BoundaryRef(i).numEdge(j);
	      if (!EdgeUsed(ne))
		{
		  EdgeUsed(ne) = true;
		  nb_edges_ref++;
		}
	    }
	}
    }

  for (int i = 0; i < mesh.GetNbVertices(); i++)
    if (VertexUsed(i))
      proj.Project(mesh.Vertex(i));
  
  mesh.ReallocateEdgesRef(nb_edges_ref);
  nb_edges_ref = 0;
  int r = mesh.GetGeometryOrder();
  R3 ptA;
  for (int i = 0; i < mesh.GetNbEdges(); i++)
    if (EdgeUsed(i))
      {
	mesh.EdgeRef(nb_edges_ref) = mesh.GetEdge(i);
	mesh.EdgeRef(nb_edges_ref).SetReference(1);
	mesh.SetStraightEdge(nb_edges_ref);
	for (int k = 0; k < r-1; k++)
	  {
	    ptA = mesh.GetPointInsideEdge(nb_edges_ref, k);
	    proj.Project(ptA);
	    mesh.SetPointInsideEdge(nb_edges_ref, k, ptA);
	  }
	
	nb_edges_ref++;
      }
  
  mesh.FindConnectivity();
  mesh.ProjectPointsOnCurves();
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    mesh.SetStraightFace(i);
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      if ((ref > 0) && (ref_cond(ref) == 1))
	for (int k = 0; k < mesh.GetNbPointsInsideFace(i); k++)
	  {
	    ptA = mesh.GetPointInsideFace(i, k);
	    proj.Project(ptA);
	    mesh.SetPointInsideFace(i, k, ptA);
	  }
    }
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc != 7)
    {
      cout << "Donnez l'ordre d'approximation, le maillage source, le maillage destination, teta, R et ri " << endl;
      cout << "Exemple : " <<endl;
      cout << "./project.x 4 coneSphere.mesh res.mesh 30  2.5 0.5 " << endl;
      abort();
    }
  
  int r = atoi(argv[1]);
  Mesh<Dimension3> mesh;
  mesh.SetGeometryOrder(r);
  
  string input_name(argv[2]), output_name(argv[3]);
  mesh.Read(input_name);
  
  double teta = pi_wp/180.0*atof(argv[4]);
  double R = atof(argv[5]);
  double ri = atof(argv[6]);
  IVect ref_cond(mesh.GetNbReferences()+1);
  ref_cond.Fill(0);
  ref_cond(1) = 1;
  
  ProjectOnConeSphere(mesh, ref_cond, teta, R, ri);
  
  mesh.Write(output_name);  
  
  return FinalizeMontjoie();
}
