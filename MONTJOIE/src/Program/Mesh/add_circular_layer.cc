#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

void AddCircularLayer(Mesh<Dimension2>& mesh, const IVect& old_num,
		      const Vector<R2>& new_points, int r, int ref_domain,
		      int ref_left, int ref_right, int ref_top, bool closed)
{
  if (old_num.GetM() <= 1)
    return;
  
  if (r*(old_num.GetM()-1) != new_points.GetM()-1)
    {
      cout << "points not compatible" << endl;
      abort();
    }
  
  if (r == 1)
    {
      int nb_old_vert = mesh.GetNbVertices();
      int nb_pts = new_points.GetM();
      int nb_edges = old_num.GetM()-1;
      if (closed)
        nb_edges++;
      
      mesh.ResizeVertices(nb_old_vert + nb_pts);
      
      // adding vertices
      int nb = nb_old_vert;
      for (int i = 0; i < nb_pts; i++)
	mesh.Vertex(nb++) = new_points(i);
      
      // then adding elements
      int nb_old_elt = mesh.GetNbElt();
      mesh.ResizeElements(nb_old_elt + nb_edges);
      int num = nb_old_elt; int n1, n2, n3, n4;
      for (int i = 0; i < nb_edges; i++)
	{
	  if ((closed) && (i == nb_edges-1))
            {
              n1 = old_num(i); n2 = old_num(0);
              n3 = nb_old_vert; n4 = nb_old_vert+i; 
            }
          else
            {
              n1 = old_num(i); n2 = old_num(i+1);
              n3 = nb_old_vert + i+1; n4 = nb_old_vert + i;
            }
	  
	  mesh.Element(num++).InitQuadrangular(n1, n2, n3, n4, ref_domain);
	}
      
      // adding edges
      int nb_new_edges = 2;
      if (closed)
        nb_new_edges = 0;
      
      if (ref_top > 0)
	nb_new_edges += nb_edges;
      
      int nb_old_edges = mesh.GetNbBoundaryRef();
      mesh.ResizeBoundariesRef(nb_old_edges + nb_new_edges);
      num = nb_old_edges;
      if (!closed)
        {
          mesh.BoundaryRef(num++).Init(old_num(0), nb_old_vert, ref_left);
          mesh.BoundaryRef(num++).Init(old_num(nb_edges), nb_old_vert+nb_edges, ref_right);
        }
      
      if (ref_top > 0)
	for (int i = 0; i < nb_edges; i++)
	  {
            if ((closed) && (i == nb_edges-1))
              mesh.BoundaryRef(num++).Init(nb_old_vert, nb_old_vert+i, ref_top); 
            else
              mesh.BoundaryRef(num++).Init(nb_old_vert+i, nb_old_vert+i+1, ref_top); 
          }
      
    }
  else if (r == 2)
    {
      int nb_old_vert = mesh.GetNbVertices();
      int nb_pts = new_points.GetM();
      int nb_edges = old_num.GetM()-1;
      int nb_pts_add = nb_edges + (nb_edges+2)/2 + nb_pts;      
      mesh.ResizeVertices(nb_old_vert + nb_pts_add);
      
      // adding vertices
      int nb = nb_old_vert;
      for (int i = 0; i < nb_pts; i++)
	mesh.Vertex(nb++) = new_points(i);
      
      Real_wp r1, r2, theta1, theta2(0);
      for (int i = 0; i < nb_edges; i++)
	{
	  int n = old_num(i);
	  CartesianToPolar(mesh.Vertex(n)(0), mesh.Vertex(n)(1), r1, theta1);
	  CartesianToPolar(new_points(2*i+2)(0), new_points(2*i+2)(1), r2, theta2);
	  
	  if (theta2 < theta1 - pi_wp)
	    theta2 += 2.0*pi_wp;
	  else if (theta2 > theta1 + pi_wp)
	    theta2 -= 2.0*pi_wp;
	  
	  Real_wp thetam = (theta1 + theta2)/2;
	  Real_wp rm = (r1+r2)/2, xm(0), ym(0);
	  
	  if (i%2 == 0)
	    {
	      // middle point on the first edge
	      xm = rm*cos(theta1), ym = rm*sin(theta1);
	      mesh.Vertex(nb++).Init(xm, ym);
	    }

	  // point in the center of the cell
	  xm = rm*cos(thetam); ym = rm*sin(thetam);
	  mesh.Vertex(nb++).Init(xm, ym);	  
	}

      if (nb_edges%2 == 0)
	{
	  // last point
	  int n_inf = old_num(nb_edges);
	  mesh.Vertex(nb++) =
	    Real_wp(0.5)*(mesh.Vertex(n_inf) + new_points(2*nb_edges));
	}
      
      // then adding elements
      nb = nb_old_vert + nb_pts;
      int nb_old_elt = mesh.GetNbElt();
      mesh.ResizeElements(nb_old_elt + 3*nb_edges);
      int num = nb_old_elt; int n1, n2, n3, n4;
      for (int i = 0; i < nb_edges; i++)
	{
	  if (i%2 == 0)
	    {
	      n1 = old_num(i); n2 = nb;
	      n3 = nb+1; n4 = old_num(i+1);
	    }
	  else
	    {
	      n1 = old_num(i); n2 = nb_old_vert + 2*i;
	      n3 = nb_old_vert + 2*i + 1; n4 = nb;
	    }
	  
	  mesh.Element(num++).InitQuadrangular(n1, n2, n3, n4, ref_domain);
	  
	  if (i%2 == 0)
	    {
	      n1 = nb; n2 = nb_old_vert + 2*i;
	      n3 = nb_old_vert + 2*i + 1; n4 = nb+1;
	    }
	  else
	    {
	      n1 = old_num(i); n2 = nb;
	      n3 = nb+1; n4 = old_num(i+1);
	    }

	  mesh.Element(num++).InitQuadrangular(n1, n2, n3, n4, ref_domain);
	  
	  if (i%2 == 0)
	    {
	      n1 = nb+1; n2 = nb_old_vert + 2*i+1;
	      n3 = nb_old_vert + 2*(i+1); n4 = old_num(i+1);
	      nb += 2;
	    }
	  else
	    {
	      n1 = nb; n2 = nb_old_vert + 2*i + 1;
	      n3 = nb_old_vert + 2*(i+1); n4 = nb+1;
	      nb++;
	    }
	  
	  mesh.Element(num++).InitQuadrangular(n1, n2, n3, n4, ref_domain);
	}
      
      // adding edges
      int nb_new_edges = 3;
      if (ref_top > 0)
	nb_new_edges += 2*nb_edges;
      
      if (nb_edges%2 == 0)
	nb_new_edges += 1;
      
      int nb_old_edges = mesh.GetNbBoundaryRef();
      mesh.ResizeBoundariesRef(nb_old_edges + nb_new_edges);
      num = nb_old_edges;
      mesh.BoundaryRef(num++).Init(old_num(0), nb_old_vert+nb_pts, ref_left);
      mesh.BoundaryRef(num++).Init(nb_old_vert, nb_old_vert+nb_pts, ref_left);
      
      if (nb_edges%2 == 0)
	{
	  mesh.BoundaryRef(num++).Init(old_num(nb_edges), mesh.GetNbVertices()-1, ref_right);
	  mesh.BoundaryRef(num++).Init(nb_old_vert+2*nb_edges, mesh.GetNbVertices()-1, ref_right);
	}
      else
	{
	  mesh.BoundaryRef(num++).Init(old_num(nb_edges), nb_old_vert+2*nb_edges, ref_right);
	}
      
      if (ref_top > 0)
	for (int i = 0; i < nb_edges; i++)
	  {
	    mesh.BoundaryRef(num++).Init(nb_old_vert+2*i, nb_old_vert+2*i+1, ref_top); 
	    mesh.BoundaryRef(num++).Init(nb_old_vert+2*i+1, nb_old_vert+2*i+2, ref_top); 
	  }
    }
  else if (r == 3)
    {
      int nb_old_vert = mesh.GetNbVertices();
      int nb_pts = new_points.GetM();
      int nb_edges = old_num.GetM()-1;
      int nb_pts_add = 2*nb_edges + nb_pts;
      mesh.ResizeVertices(nb_old_vert + nb_pts_add);
      
      // adding vertices
      int nb = nb_old_vert;
      for (int i = 0; i < nb_pts; i++)
	mesh.Vertex(nb++) = new_points(i);
      
      Real_wp r1, r2, theta1, theta2;
      for (int i = 0; i < nb_edges; i++)
	{
	  int n = old_num(i);
	  CartesianToPolar(mesh.Vertex(n)(0), mesh.Vertex(n)(1), r1, theta1);
	  CartesianToPolar(new_points(3*i+3)(0), new_points(3*i+3)(1), r2, theta2);
	  
	  if (theta2 < theta1 - pi_wp)
	    theta2 += 2.0*pi_wp;
	  else if (theta2 > theta1 + pi_wp)
	    theta2 -= 2.0*pi_wp;
	  
	  Real_wp thetam = (2.0*theta1 + theta2)/3;
	  Real_wp rm = (r1+r2)/2;
	  
	  Real_wp xm = rm*cos(thetam), ym = rm*sin(thetam);
	  mesh.Vertex(nb++).Init(xm, ym);

	  thetam = (theta1 + 2.0*theta2)/3;
	  xm = rm*cos(thetam); ym = rm*sin(thetam);
	  mesh.Vertex(nb++).Init(xm, ym);	  
	}
      
      // then adding elements
      nb = nb_old_vert + nb_pts;
      int nb_old_elt = mesh.GetNbElt();
      mesh.ResizeElements(nb_old_elt + 4*nb_edges);
      int num = nb_old_elt; int n1, n2, n3, n4;
      for (int i = 0; i < nb_edges; i++)
	{
	  n1 = old_num(i); n2 = nb;
	  n3 = nb_old_vert + 3*i + 1; n4 = nb_old_vert + 3*i;
	  	  
	  mesh.Element(num++).InitQuadrangular(n1, n2, n3, n4, ref_domain);
	  
	  n1 = old_num(i); n2 = old_num(i+1);
	  n3 = nb+1; n4 = nb;
	  
	  mesh.Element(num++).InitQuadrangular(n1, n2, n3, n4, ref_domain);
	  
	  n1 = nb; n2 = nb+1;
	  n3 = nb_old_vert + 3*i+2; n4 = nb_old_vert + 3*i + 1;
	  	  
	  mesh.Element(num++).InitQuadrangular(n1, n2, n3, n4, ref_domain);
	  
	  n1 = nb+1; n2 = old_num(i+1);
	  n3 = nb_old_vert + 3*i+3; n4 = nb_old_vert + 3*i+2;
	  
	  mesh.Element(num++).InitQuadrangular(n1, n2, n3, n4, ref_domain);
	  nb += 2;
	}
      
      // adding edges
      int nb_new_edges = 2;
      if (ref_top > 0)
	nb_new_edges += 3*nb_edges;
      
      int nb_old_edges = mesh.GetNbBoundaryRef();
      mesh.ResizeBoundariesRef(nb_old_edges + nb_new_edges);
      num = nb_old_edges;
      mesh.BoundaryRef(num++).Init(old_num(0), nb_old_vert, ref_left);
      mesh.BoundaryRef(num++).Init(old_num(nb_edges), nb_old_vert+3*nb_edges, ref_right);
      
      if (ref_top > 0)
	for (int i = 0; i < nb_edges; i++)
	  {
	    mesh.BoundaryRef(num++).Init(nb_old_vert+3*i, nb_old_vert+3*i+1, ref_top); 
	    mesh.BoundaryRef(num++).Init(nb_old_vert+3*i+1, nb_old_vert+3*i+2, ref_top); 
	    mesh.BoundaryRef(num++).Init(nb_old_vert+3*i+2, nb_old_vert+3*i+3, ref_top); 
	  }      
    }
  else
    {
      cout << "not implemented" << endl;
      abort();
    }
}

Real_wp ExtractCircularReference(const Mesh<Dimension2>& mesh, int ref, IVect& old_num,
				 VectReal_wp& theta_pts, int& ref_left, int& ref_right)
{
  Vector<bool> PtUsed(mesh.GetNbVertices());
  PtUsed.Fill(false);
  int nb_pts = 0;
  // extracting all the points on the edges of reference ref
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    if (mesh.BoundaryRef(i).GetReference() == ref)
      {
	int n1 = mesh.BoundaryRef(i).numVertex(0);
	int n2 = mesh.BoundaryRef(i).numVertex(1);
	if (!PtUsed(n1))
	  {
	    PtUsed(n1) = true;
	    nb_pts++;
	  }
	
	if (!PtUsed(n2))
	  {
	    PtUsed(n2) = true;
	    nb_pts++;
	  }
      }
  
  old_num.Reallocate(nb_pts);
  theta_pts.Reallocate(nb_pts);
  int nb = 0; Real_wp r0(0);
  for (int i = 0; i < PtUsed.GetM(); i++)
    if (PtUsed(i))
      {
	old_num(nb) = i;
	CartesianToPolar(mesh.Vertex(i)(0), mesh.Vertex(i)(1), r0, theta_pts(nb));
	nb++;
      }
  
  // sorting by theta
  Sort(theta_pts, old_num);
  
  // we retrieve the reference on the left and right extremity
  ref_left = 0; ref_right = 0;
  
  IVect all_ref;
  mesh.GetAllReferencesAroundPoint(old_num(0), all_ref);
  for (int k = 0; k < all_ref.GetM(); k++)
    if (all_ref(k) != ref)
      ref_left = all_ref(k);
  
  all_ref.Clear();
  mesh.GetAllReferencesAroundPoint(old_num(old_num.GetM()-1), all_ref);
  for (int k = 0; k < all_ref.GetM(); k++)
    if (all_ref(k) != ref)
      ref_right = all_ref(k);  
    
  return r0;
}

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc < 5)
    {
      cout<<"Cette commande demande au moins cinq arguments" << endl;
      cout<<"add_circular_layer initial.mesh final.mesh ref file_r_teta" << endl;
      cout<<"the file final.mesh is created" << endl;
      return -1;
    }
  
  string file_input(argv[1]);
  string file_output(argv[2]);
  int ref = atoi(argv[3]);
  string file_r(argv[4]);
  
  VectReal_wp step_r, step_teta;
  IVect ref_top;
  ifstream file_in(file_r.data());
  
  if (!file_in.is_open())
    {
      cout << "File " << file_r << " does not exist" << endl;
      abort();
    }

  string ligne; VectString mots;
  while (!file_in.eof())
    {
      getline(file_in, ligne);
      if (file_in.good())
	{
	  StringTokenize(ligne, mots, " \t");
	  if (mots.GetM() >= 2)
	    {
	      step_r.PushBack(to_num<Real_wp>(mots(0)));
	      step_teta.PushBack(to_num<Real_wp>(mots(1)));
	      if (mots.GetM() == 3)
		ref_top.PushBack(to_num<int>(mots(2)));
	      else
		ref_top.PushBack(0);
	    }
	}
    }

  file_in.close();
  
  Mesh<Dimension2> mesh;
  mesh.Read(file_input);
  
  // extracting points on the circular boundary
  IVect old_num; VectReal_wp theta_pts; Real_wp r1;
  int ref_left, ref_right;
  r1 = ExtractCircularReference(mesh, ref, old_num, theta_pts, ref_left, ref_right);
  cout << "Initial radius = " << r1 << endl;

  bool closed = false;
  int n = mesh.FindEdgeRefWithExtremities(old_num(0), old_num(old_num.GetM()-1));
  if (n >= 0)
    closed = true;
  
  int nb_layers = step_r.GetM();
  int raff_global = 1;
  
  // loops on layers to add
  for (int p = 0; p < nb_layers; p++)
    {
      Real_wp r2 = step_r(p);
      int raff = toInteger(ceil(step_teta(p)/raff_global));
      raff_global *= raff;
      cout << "At radius r = " << r2 << " theta-refinement = " << raff << endl;
      
      // points on the extern boundary of the new layer
      Vector<R2> new_pts(raff*(old_num.GetM()-1)+1);
      int nb = 0; Real_wp theta1, theta2(0), dtheta;
      VectReal_wp new_theta(new_pts.GetM());
      IVect new_num(new_pts.GetM());
      int nb_old = mesh.GetNbVertices();
      for (int i = 0; i < old_num.GetM()-1; i++)
	{
	  theta1 = theta_pts(i);
	  theta2 = theta_pts(i+1);
	  dtheta = (theta2 - theta1);
	  if (raff > 1)
	    dtheta /= raff;
	  
	  for (int j = 0; j < raff; j++)
	    {
	      Real_wp thetam = theta1 + j*dtheta;
	      new_num(nb) = nb_old + nb;
	      new_theta(nb) = thetam;
	      new_pts(nb++).Init(r2*cos(thetam), r2*sin(thetam));
	    }
	}
      
      new_num(nb) = nb_old + nb;
      new_theta(nb) = theta2;
      new_pts(nb++).Init(r2*cos(theta2), r2*sin(theta2));
      
      if (closed && (raff > 1))
        {
          cout << "Not implemented " << endl;
          abort();
        }

      // completing mesh until the extern boundary
      int ref_domain = 1;
      
      AddCircularLayer(mesh, old_num, new_pts, raff, ref_domain,
		       ref_left, ref_right, ref_top(p), closed);
      
      //DISP(new_num); DISP(new_theta);
      // updating old_num, theta_pts
      old_num = new_num;
      theta_pts = new_theta;
    }
  
  mesh.Write(file_output);
}  
  
