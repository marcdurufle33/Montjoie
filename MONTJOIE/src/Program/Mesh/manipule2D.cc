#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

void AddVertex(Mesh<Dimension2>& mesh, const R2& point)
{
  int nv = mesh.GetNbVertices();
  mesh.ResizeVertices(nv+1);
  mesh.Vertex(nv) = point;
  cout<<"Numero du nouveau point"<<nv+1<<endl;
}

void AddEdge(Mesh<Dimension2>& mesh, int n1, int n2, int ref)
{
  int nb_vert = mesh.GetNbVertices();
  int n = mesh.GetNbBoundaryRef();
  mesh.ResizeBoundariesRef(n+1);
  
  n1--;
  if ((n1 < 0)||(n1 >= nb_vert))
    {
      cout<<"Vertex number must be between 1 and "<<nb_vert<<endl;
      return;
    }
  
  n2--;
  if ((n2 < 0)||(n2 >= nb_vert))
    {
      cout<<"Vertex number must be between 1 and "<<nb_vert<<endl;
      return;
    }
  
  mesh.BoundaryRef(n).Init(n1, n2, ref);
  cout<<"Numero de la nouvelle arete "<<n+1<<endl;
}

void AddElement(Mesh<Dimension2>& mesh, IVect num, int ref)
{
  int n = mesh.GetNbElt();
  int nb_vert = mesh.GetNbVertices();
  mesh.ResizeElements(n+1);
  for (int i = 0; i < num.GetM(); i++)
    {
      num(i)--;
     if ((num(i) < 0)||(num(i) >= nb_vert))
       {
	 cout<<"Vertex number must be between 1 and "<<nb_vert<<endl;
	 return;
       }
    }
  
  mesh.Element(n).Init(num, ref);
  cout<<"Numero du nouvel element "<<n+1<<endl;
}

void RemoveVertex(Mesh<Dimension2>& mesh, const Vector<bool>& VertexUsed)
{
  int nb_vert = mesh.GetNbVertices();
  IVect new_number(nb_vert);
  nb_vert = 0;
  for (int i = 0; i < VertexUsed.GetM(); i++)
    {
      if (VertexUsed(i))
	{
	  mesh.Vertex(nb_vert) = mesh.Vertex(i);
	  new_number(i) = nb_vert++;
	}
      else
	new_number(i) = -1;
    }
  
  mesh.ResizeVertices(nb_vert);
    
  // on modifie les aretes de reference
  int nb_edge = 0;
  int r = mesh.GetGeometryOrder();
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int n1 = mesh.BoundaryRef(i).numVertex(0);
      int n2 = mesh.BoundaryRef(i).numVertex(1);
      if ( (VertexUsed(n1)) && (VertexUsed(n2)) )
	{
	  int ref = mesh.BoundaryRef(i).GetReference();
	  mesh.BoundaryRef(nb_edge).Init(new_number(n1), new_number(n2), ref);
	    
	  if (nb_edge != i)
	    {
	      for (int k = 0; k < r-1; k++)
		mesh.SetPointInsideEdge(nb_edge, k, mesh.GetPointInsideEdge(i, k));
	    }
	  
	  nb_edge++;
	}
    }
  
  mesh.ResizeBoundariesRef(nb_edge);
  
  // on modifie les elements
  int nb_elt = 0;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      bool keep_elt = true;
      nb_vert = mesh.Element(i).GetNbVertices();
      int ref = mesh.Element(i).GetReference();
      IVect numv(nb_vert);
      for (int j = 0; j < nb_vert; j++)
	{
	  int nv = mesh.Element(i).numVertex(j);
	  numv(j) = new_number(nv);
	  if (!VertexUsed(nv))
	    keep_elt = false;
	}
      
      if (keep_elt)
	mesh.Element(nb_elt++).Init(numv, ref);      
    }
  
  if (nb_elt != mesh.GetNbElt())
    mesh.ResizeElements(nb_elt);
  
  mesh.FindConnectivity();
}


void RemoveVertex(Mesh<Dimension2>& mesh, IVect num)
{
  int nb_vert = mesh.GetNbVertices();
  for (int i = 0; i < num.GetM(); i++)
    {
      num(i)--;
      if ((num(i) < 0)||(num(i) >= nb_vert))
	{
	  cout<<"Vertex number must be between 1 and "<<nb_vert<<endl;
	  return;
	}
    }
  
  Vector<bool> VertexUsed(nb_vert); VertexUsed.Fill(true);
  for (int i = 0; i < num.GetM(); i++)
    VertexUsed(num(i)) = false;
  
  RemoveVertex(mesh, VertexUsed);
}

void RemoveVertex(Mesh<Dimension2>& mesh, const VectReal_wp& delim)
{
  if (delim.GetM() < 4)
    {
      cout<<"Enter at least four real numbers : xmin, xmax, ymin, ymax"<<endl;
      return;
    }
  
  Real_wp xmin = delim(0), xmax = delim(1), ymin = delim(2), ymax = delim(3);
  int nb_vert = mesh.GetNbVertices();
  Vector<bool> VertexUsed(nb_vert);
  for (int i = 0; i < nb_vert; i++)
    {
      if ((mesh.Vertex(i)(0) > xmin)&&(mesh.Vertex(i)(0) < xmax)&&
	  (mesh.Vertex(i)(1) > ymin)&&(mesh.Vertex(i)(1) < ymax) )
	VertexUsed(i) = false;
      else
	VertexUsed(i) = true;
    }
  
  RemoveVertex(mesh, VertexUsed);
}

void RemoveEdge(Mesh<Dimension2>& mesh, const Vector<bool>& EdgeUsed)
{
  int nb_edge = 0;
  int r = mesh.GetGeometryOrder();
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      if (EdgeUsed(i))
	{
	  if (nb_edge != i)
	    {
              mesh.BoundaryRef(nb_edge) = mesh.BoundaryRef(i);
              for (int k = 0; k < r-1; k++)
		mesh.SetPointInsideEdge(nb_edge, k, mesh.GetPointInsideEdge(i, k));
	      
            }
          
	  nb_edge++;
	}
    }
  
  if (nb_edge != mesh.GetNbBoundaryRef())
    mesh.ResizeBoundariesRef(nb_edge);
  
}

void RemoveEdge(Mesh<Dimension2>& mesh, IVect num)
{
  int nb_edge = mesh.GetNbBoundaryRef();
  Vector<bool> EdgeUsed(nb_edge); EdgeUsed.Fill(true);
  for (int i = 0; i < num.GetM(); i++)
    {
      num(i)--;
      if ((num(i) < 0)||(num(i) >= nb_edge))
	{
	  cout<<"Edge number must be between 1 and "<<nb_edge<<endl;
	  return;
	}
      
      EdgeUsed(num(i)) = false;
    }
  
  RemoveEdge(mesh, EdgeUsed);
}

void RemoveEdgeRef(Mesh<Dimension2>& mesh, IVect num)
{
  int nb_edge = mesh.GetNbBoundaryRef();
  Vector<bool> EdgeUsed(nb_edge); EdgeUsed.Fill(true);
  int ref_max = 0;
  for (int i = 0; i < num.GetM(); i++)
    ref_max = max(ref_max, num(i));
  
  Vector<bool> RemoveRef(ref_max+1); RemoveRef.Fill(false);
  for (int i = 0; i < num.GetM(); i++)
    RemoveRef(num(i)) = true;
  
  for (int i = 0; i < nb_edge; i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      if (ref >= 0)
        if (ref <= ref_max)
          if (RemoveRef(ref))
            EdgeUsed(i) = false;
    }
  
  RemoveEdge(mesh, EdgeUsed);
}

void RemoveElement(Mesh<Dimension2>& mesh, const Vector<bool>& ElementUsed)
{
  int nb_elt = 0;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      if (ElementUsed(i))
	{
	  if (nb_elt != i)
	    mesh.Element(nb_elt) = mesh.Element(i);
	  
	  nb_elt++;
	}
    }
  
  if (nb_elt != mesh.GetNbElt())
    mesh.ResizeElements(nb_elt);

}

void RemoveElement(Mesh<Dimension2>& mesh, IVect num)
{
  int nb_elt = mesh.GetNbElt();
  Vector<bool> ElementUsed(nb_elt); ElementUsed.Fill(true);
  for (int i = 0; i < num.GetM(); i++)
    {
      num(i)--;
      if ((num(i) < 0)||(num(i) >= nb_elt))
	{
	  cout<<"Element number must be between 1 and "<<nb_elt<<endl;
	  return;
	}
      
      ElementUsed(num(i)) = false;
    }
  
  RemoveElement(mesh, ElementUsed);
}

void RemoveElementRef(Mesh<Dimension2>& mesh, IVect num)
{
  int nb_elt = mesh.GetNbElt();
  Vector<bool> ElementUsed(nb_elt); ElementUsed.Fill(true);
  int ref_max = 0;
  for (int i = 0; i < num.GetM(); i++)
    ref_max = max(ref_max, num(i));
  
  Vector<bool> RemoveRef(ref_max+1); RemoveRef.Fill(false);
  for (int i = 0; i < num.GetM(); i++)
    RemoveRef(num(i)) = true;
  
  for (int i = 0; i < nb_elt; i++)
    {
      int ref = mesh.Element(i).GetReference();
      if (ref >= 0)
        if (ref <= ref_max)
          if (RemoveRef(ref))
            ElementUsed(i) = false;
    }
  
  RemoveElement(mesh, ElementUsed);
}

void ModifyVertex(Mesh<Dimension2>& mesh, int n, const R2& pt)
{
  if ((n <= 0)||(n > mesh.GetNbVertices()))
    {
      cout<<"Vertex number must be between 1 and "<<mesh.GetNbVertices()<<endl;
      return;
    }
  
  mesh.Vertex(n-1) = pt;
}

void ModifyEdge(Mesh<Dimension2>& mesh, int n, int n1, int n2)
{
  if ((n <= 0)||(n > mesh.GetNbBoundaryRef()))
    {
      cout<<"Edge number must be between 1 and "<<mesh.GetNbBoundaryRef()<<endl;
      return;
    }
  
  n1--;
  if ((n1 < 0)||(n1 >= mesh.GetNbVertices()))
    {
      cout<<"Vertex number must be between 1 and "<<mesh.GetNbVertices()<<endl;
      return;
    }

  n2--;
  if ((n2 < 0)||(n2 >= mesh.GetNbVertices()))
    {
      cout<<"Vertex number must be between 1 and "<<mesh.GetNbVertices()<<endl;
      return;
    }
    
  mesh.BoundaryRef(n-1).Init(n1, n2, mesh.BoundaryRef(n-1).GetReference());
}

void ModifyEdgeRef(Mesh<Dimension2>& mesh, int n, int ref)
{
  if ((n <= 0)||(n > mesh.GetNbBoundaryRef()))
    {
      cout<<"Edge number must be between 1 and "<<mesh.GetNbBoundaryRef()<<endl;
      return;
    }
    
  mesh.BoundaryRef(n-1).SetReference(ref);
}

void ModifyElement(Mesh<Dimension2>& mesh, int n, IVect num)
{
  if ((n <= 0)||(n > mesh.GetNbElt()))
    {
      cout<<"Element number must be between 1 and "<<mesh.GetNbElt()<<endl;
      return;
    }
  
  for (int i = 0; i < num.GetM(); i++)
    {
      num(i)--;
      if ((num(i) < 0)||(num(i) >= mesh.GetNbVertices()))
	{
	  cout<<"Vertex number must be between 1 and "<<mesh.GetNbVertices()<<endl;
	  return;
	}
    }
    
  mesh.Element(n-1).Init(num, mesh.Element(n-1).GetReference());
}

void ModifyElementRef(Mesh<Dimension2>& mesh, int n, int ref)
{
  if ((n <= 0)||(n > mesh.GetNbElt()))
    {
      cout<<"Element number must be between 1 and "<<mesh.GetNbElt()<<endl;
      return;
    }
    
  mesh.Element(n-1).SetReference(ref);
}


void ModifySurfaceRef(Mesh<Dimension2>& mesh, int old_ref, int new_ref, int coor_mode)
{
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    if (mesh.BoundaryRef(i).GetReference() == old_ref)
      {
	if (coor_mode == -1)
	  mesh.BoundaryRef(i).SetReference(new_ref);  
	else
	  {
	    R2 ptA, ptB;
	    ptA = mesh.Vertex(mesh.BoundaryRef(i).numVertex(0));
	    ptB = mesh.Vertex(mesh.BoundaryRef(i).numVertex(1));
	    R2 center = 0.5*(ptA + ptB);
	    if (center(coor_mode) < 0)
	      mesh.BoundaryRef(i).SetReference(new_ref);  
	  }
      }
}

void ModifyVolumeRef(Mesh<Dimension2>& mesh, int old_ref, int new_ref, bool modify_all)
{
  for (int i = 0; i < mesh.GetNbElt(); i++)
    if (mesh.Element(i).GetReference() == old_ref)
      { 
	if(modify_all)
	  mesh.Element(i).SetReference(new_ref);  
	else
	  {
	    R2 ptA, ptB;
	    ptA = mesh.Vertex(mesh.Element(i).numVertex(0));
	    ptB = mesh.Vertex(mesh.Element(i).numVertex(1));
	    R2 center = 0.5*(ptA + ptB); 
	    DISP(center);
	    if (center(1) < 0)
	      mesh.Element(i).SetReference(new_ref);
	  }
      }
}

template<class T>
void ReadParameters(Vector<T>& param)
{
  string chaine;
  // cin>>flush;
  // getline(cin, chaine ,'\n');
  getline(cin, chaine ,'\n');
  istringstream stream_data(chaine);
  
  param.ReadText(stream_data);
}

void ReadParameters(string& param)
{
  getline(cin, param ,'\n');
}


void SaveMesh(Mesh<Dimension2>& mesh, const string& output_name)
{
  // we detect unused vertices
  int nb_vert = mesh.GetNbVertices();
  Vector<bool> VertexUsed(nb_vert); VertexUsed.Fill(false);
  for (int i = 0; i < mesh.GetNbElt(); i++)
    for (int j = 0; j < mesh.Element(i).GetNbVertices(); j++)
      VertexUsed(mesh.Element(i).numVertex(j)) = true;
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    for (int j = 0; j < 2; j++)
      VertexUsed(mesh.BoundaryRef(i).numVertex(j)) = true;
  
  // we remove those vertices
  RemoveVertex(mesh, VertexUsed);
  
  mesh.Write(output_name);  
}


template<class T>
void TransformPolarToCartesianMesh(Mesh<Dimension2>& mesh, const T& rmin, const T& rmax,
                                   const T& teta_min, const T& teta_max)
{
  Real_wp xmin = mesh.GetXmin();
  Real_wp xmax = mesh.GetXmax();
  Real_wp ymin = mesh.GetYmin();
  Real_wp ymax = mesh.GetYmax();
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    {
      Real_wp lambda = (mesh.Vertex(i)(1) - ymin)/(ymax-ymin);
      Real_wp teta = (1.0-lambda)*teta_min + lambda*teta_max;

      lambda = (mesh.Vertex(i)(0) - xmin)/(xmax-xmin);
      Real_wp r = (1.0-lambda)*rmin + lambda*rmax;
      
      Real_wp x = r*cos(teta);
      Real_wp y = r*sin(teta);
      mesh.Vertex(i).Init(x, y);
    }
  
  int r = mesh.GetGeometryOrder();
  R2 ptA, ptB;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    for (int k = 0; k < r-1; k++)
      {
        ptA = mesh.GetPointInsideEdge(i, k);
        Real_wp lambda = (ptA(1) - ymin)/(ymax-ymin);
        Real_wp teta = (1.0-lambda)*teta_min + lambda*teta_max;
        
        lambda = (ptA(0) - xmin)/(xmax-xmin);
        Real_wp r = (1.0-lambda)*rmin + lambda*rmax;
        
        Real_wp x = r*cos(teta);
        Real_wp y = r*sin(teta);
        ptB.Init(x, y);
        mesh.SetPointInsideEdge(i, k, ptB);
      }
}

void AddRefinedLayer(Mesh<Dimension2>& mesh, int ref1, int ref2,
                     int ref_edge1, int ref_edge2, bool add_vertex, int ref_domain)
{
  IVect num1, num2; int nb1 = 0, nb2 = 0;
  Vector<int> IndexPoint(mesh.GetNbVertices());
  IndexPoint.Fill(-1);
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int n1 = mesh.BoundaryRef(i).numVertex(0);
      int n2 = mesh.BoundaryRef(i).numVertex(1);
      int ref = mesh.BoundaryRef(i).GetReference();
      if (ref == ref1)
        {
          if (IndexPoint(n1) == -1)
            {
              IndexPoint(n1) = ref1;
              nb1++;
            }
          
          if (IndexPoint(n2) == -1)
            {
              IndexPoint(n2) = ref1;
              nb1++;
            }
        }
      else if (ref == ref2)
        {
          if (IndexPoint(n1) == -1)
            {
              IndexPoint(n1) = ref2;
              nb2++;
            }
          
          if (IndexPoint(n2) == -1)
            {
              IndexPoint(n2) = ref2;
              nb2++;
            }
        }
    }
  
  VectR2 ListPoint1(nb1), ListPoint2(nb2);
  num1.Reallocate(nb1); num2.Reallocate(nb2);
  nb1 = 0; nb2 = 0;
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    {
      if (IndexPoint(i) == ref1)
        {
          num1(nb1) = i;
          ListPoint1(nb1) = mesh.Vertex(i);
          nb1++;
        }

      if (IndexPoint(i) == ref2)
        {
          num2(nb2) = i;
          ListPoint2(nb2) = mesh.Vertex(i);
          nb2++;
        }
    }

  Sort(ListPoint1, num1);
  Sort(ListPoint2, num2);

  if (nb1 != 2*nb2-1)
    {
      cout << "Refined layer incorrect" << endl;
      abort();
    }

  // creation des sommets aux milieux des aretes
  int n0 = 1;
  if (add_vertex)
    n0 = 0;
  
  int off_v = mesh.GetNbVertices();
  int nb_new_vert = 0;
  for (int i = n0; i < nb2; i += 2)
    nb_new_vert++;
  
  mesh.ResizeVertices(off_v + nb_new_vert + nb2-1);
  nb_new_vert = 0; R2 ptA;
  for (int i = n0; i < nb2; i += 2)
    {
      ptA = Real_wp(0.5)*(ListPoint2(i) + ListPoint1(2*i));
      mesh.Vertex(off_v + nb_new_vert) = ptA;
      nb_new_vert++;
    }

  // creation des sommets aux milieux des elements
  // et les quadrilateres
  int off_elt = mesh.GetNbElt();
  mesh.ResizeElements(off_elt + 3*(nb2-1));
  int ne = 0;
  for (int i = 0; i < nb2-1; i++)
    {
      ptA = Real_wp(0.25)*(ListPoint2(i) + ListPoint2(i+1) + ListPoint1(2*i) + ListPoint1(2*i+2));
      int nv = off_v + nb_new_vert + i;
      mesh.Vertex(nv) = ptA;
      
      if (i%2 == n0)
        {
          mesh.Element(off_elt + 3*i).InitQuadrangular(num1(2*i), num1(2*i+1), nv, off_v + ne, ref_domain);
          mesh.Element(off_elt + 3*i+1).InitQuadrangular(off_v + ne, nv, num2(i+1), num2(i), ref_domain);
          mesh.Element(off_elt + 3*i+2).InitQuadrangular(num1(2*i+1), num1(2*i+2), num2(i+1), nv, ref_domain);
          ne++;
        }
      else
        {
          mesh.Element(off_elt + 3*i).InitQuadrangular(num1(2*i), num1(2*i+1), nv, num2(i), ref_domain);
          mesh.Element(off_elt + 3*i+1).InitQuadrangular(nv, off_v + ne, num2(i+1), num2(i), ref_domain);
          mesh.Element(off_elt + 3*i+2).InitQuadrangular(num1(2*i+1), num1(2*i+2), off_v + ne, nv, ref_domain);
        }
    }

  // creation des aretes de bord
  int nb_edges_ref = mesh.GetNbBoundaryRef();
  if ((ref_edge1 != 0) && (ref_edge2 != 0))
    {
      if (add_vertex)
        {
          mesh.ResizeBoundariesRef(nb_edges_ref+4);
          mesh.BoundaryRef(nb_edges_ref).Init(num2(0), off_v, ref_edge1);
          mesh.BoundaryRef(nb_edges_ref+1).Init(off_v, num1(0), ref_edge1);
          mesh.BoundaryRef(nb_edges_ref+2).Init(num2(nb2-1), off_v+nb_new_vert-1, ref_edge2);
          mesh.BoundaryRef(nb_edges_ref+3).Init(off_v+nb_new_vert-1, num1(2*nb2-2), ref_edge2);
        }
      else
        {
          mesh.ResizeBoundariesRef(nb_edges_ref+2);
          mesh.BoundaryRef(nb_edges_ref).Init(num2(0), num1(0), ref_edge1);
          mesh.BoundaryRef(nb_edges_ref+1).Init(num2(nb2-1), num1(2*nb2-2), ref_edge2);
        }
    }

  mesh.SortBoundariesRef();
  mesh.FindConnectivity();
}

void CreateLocallyRefinedRegularMesh(const R2& ptMin, const R2& ptMax,
                                     const R2& pt0, const R2& ptN, const Real_wp& h,
                                     int ref_domain, const TinyVector<int, 4>& ref_boundary,
                                     Mesh<Dimension2>& mesh)
{
  mesh.ResizeNbReferences(20);
  // on cree la zone centrale
  int nx = toInteger(round((ptN(0) - pt0(0))/h)) + 1;
  int ny = toInteger(round((ptN(1) - pt0(1))/h)) + 1;
  if (nx%2 != 1)
    nx++;

  if (ny%2 != 1)
    ny++;

  TinyVector<int, 4> ref_bound(1, 2, 3, 4);
  mesh.CreateRegularMesh(pt0, ptN, TinyVector<int, 2>(2*nx-1, 2*ny-1),
                         ref_domain, ref_bound, mesh.QUADRILATERAL_MESH);

  //mesh.Write("carre.mesh");
  
  // on cree un carre a gauche
  int nx1 = toInteger(round((pt0(0) - ptMin(0))/h)) + 1;
  Real_wp xm = ptMin(0);
  if (nx1 > 2)
    {
      Real_wp dx = (pt0(0) - ptMin(0)) / (nx1 - 1);
      xm = pt0(0) - dx;
      nx1--;
    }
  
  Mesh<Dimension2> reg_mesh; reg_mesh.ResizeNbReferences(20);
  ref_bound.Init(0, 5, 0, 14);
  reg_mesh.CreateRegularMesh(R2(ptMin(0), pt0(1)), R2(xm, ptN(1)), TinyVector<int, 2>(nx1, ny),
                             ref_domain, ref_bound, mesh.QUADRILATERAL_MESH);

  //reg_mesh.Write("carre_left.mesh");
  
  mesh.AppendMesh(reg_mesh);
  AddRefinedLayer(mesh, 4, 5, 1, 3, true, ref_domain);

  //mesh.Write("carre_merge.mesh");
  // on cree un carre a droite
  int nx2 = toInteger(round((ptMax(0) - ptN(0))/h)) + 1;
  Real_wp xp = ptN(0);
  if (nx2 > 2)
    {
      Real_wp dx = (ptMax(0) - ptN(0)) / (nx2 - 1);
      xp = ptN(0) + dx;
      nx2--;
    }

  ref_bound.Init(0, 12, 0, 6); 
  reg_mesh.CreateRegularMesh(R2(xp, pt0(1)), R2(ptMax(0), ptN(1)), TinyVector<int, 2>(nx2, ny),
                             ref_domain, ref_bound, mesh.QUADRILATERAL_MESH);
  
  mesh.AppendMesh(reg_mesh);
  AddRefinedLayer(mesh, 2, 6, 1, 3, true, ref_domain);
  
  // on rajoute les quatre coins
  int ny1 = toInteger(round((pt0(1) - ptMin(1))/h)) + 1;
  int ny2 = toInteger(round((ptMax(1) - ptN(1))/h)) + 1;
  ref_bound.Init(11, 0, 0, 14);
  reg_mesh.CreateRegularMesh(R2(ptMin(0), ptMin(1)), R2(xm, pt0(1)), TinyVector<int, 2>(nx1, ny1),
                             ref_domain, ref_bound, mesh.QUADRILATERAL_MESH);

  mesh.AppendMesh(reg_mesh);

  ref_bound.Init(0, 0, 13, 14);
  reg_mesh.CreateRegularMesh(R2(ptMin(0), ptN(1)), R2(xm, ptMax(1)), TinyVector<int, 2>(nx1, ny2),
                             ref_domain, ref_bound, mesh.QUADRILATERAL_MESH);
  
  mesh.AppendMesh(reg_mesh);

  ref_bound.Init(11, 12, 0, 0);
  reg_mesh.CreateRegularMesh(R2(xp, ptMin(1)), R2(ptMax(0), pt0(1)), TinyVector<int, 2>(nx2, ny1),
                             ref_domain, ref_bound, mesh.QUADRILATERAL_MESH);

  mesh.AppendMesh(reg_mesh);

  ref_bound.Init(0, 12, 13, 0);
  reg_mesh.CreateRegularMesh(R2(xp, ptN(1)), R2(ptMax(0), ptMax(1)), TinyVector<int, 2>(nx2, ny2),
                             ref_domain, ref_bound, mesh.QUADRILATERAL_MESH);
  
  mesh.AppendMesh(reg_mesh);
  
  // on cree les deux dernieres regions verticales
  Real_wp ym = ptMin(1);
  if (ny1 > 2)
    {
      Real_wp dy = (pt0(1) - ptMin(1))/(ny1-1);
      ym = pt0(1) - dy;
      ny1--;
    }
  
  ref_bound.Init(11, 0, 7, 0);  
  reg_mesh.CreateRegularMesh(R2(xm, ptMin(1)), R2(xp, ym), TinyVector<int, 2>(nx+2, ny1),
                             ref_domain, ref_bound, mesh.QUADRILATERAL_MESH);

  mesh.AppendMesh(reg_mesh);
  AddRefinedLayer(mesh, 1, 7, 0, 0, false, ref_domain);
  
  Real_wp yp = ptMax(1);
  if (ny2 > 2)
    {
      Real_wp dy = (ptMax(1) - ptN(1))/(ny2-1);
      yp = ptN(1) + dy;
      ny2--;
    }
  
  ref_bound.Init(8, 0, 13, 0);  
  reg_mesh.CreateRegularMesh(R2(xm, yp), R2(xp, ptMax(1)), TinyVector<int, 2>(nx+2, ny2),
                             ref_domain, ref_bound, mesh.QUADRILATERAL_MESH);
  
  mesh.AppendMesh(reg_mesh);
  AddRefinedLayer(mesh, 3, 8, 0, 0, false, ref_domain);

  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      if (ref < 10)
        mesh.BoundaryRef(i).SetReference(0);
      else
        mesh.BoundaryRef(i).SetReference(ref-10);
    }

  mesh.ClearZeroBoundaryRef();  
}

void RefineLocallyMesh(Mesh<Dimension2>& mesh, const Vector<bool>& vertex_to_refine,
                       const Vector<bool>& edge_to_refine, bool anisotrope)
{
  TinyVector<int, 4> num;
  int nb_new_vert = 0, nb_new_elt = 0;
  // on compte le nombre de sommets et d'aretes a rajouter
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      if (mesh.Element(i).GetNbVertices() != 4)
        {
          cout << "Triangular elements not allowed" << endl;
          abort();
        }

      num(0) = mesh.Element(i).numVertex(0);
      num(1) = mesh.Element(i).numVertex(1);
      num(2) = mesh.Element(i).numVertex(2);
      num(3) = mesh.Element(i).numVertex(3);
      int nb_vert_refined = 0;
      for (int j = 0; j < 4; j++)
        if (vertex_to_refine(num(j)))
          nb_vert_refined++;

      int nb_edge_refined = 0;
      for (int j = 0; j < 4; j++)
        if (edge_to_refine(mesh.Element(i).numEdge(j)))
          nb_edge_refined++;

      if (nb_vert_refined == 1)
        {
          nb_new_vert++;
          nb_new_elt += 2;
        }
      else if (nb_vert_refined >= 2)
        {
          if (anisotrope)
            {
              if (nb_vert_refined == 2)
                {
                  if (nb_edge_refined == 0)
                    nb_new_elt++;
                  else
                    {
                      nb_new_vert += 4;
                      nb_new_elt += 6;
                    }
                }
              else if (nb_vert_refined == 3)
                {
                  cout << "Impossible" << endl;
                  abort();
                }
              else
                {
                  if (nb_edge_refined == 2)
                      nb_new_elt += 2;
                  else
                    {
                      nb_new_vert += 4;
                      nb_new_elt += 8;
                    }
                }              
            }
          else
            {
              nb_new_vert += 4;
              nb_new_elt += 4 + nb_vert_refined;
            }
        }
    }

  int nb_new_edges = 0;
  for (int i = 0; i < mesh.GetNbEdges(); i++)
    {
      num(0) = mesh.GetEdge(i).numVertex(0);
      num(1) = mesh.GetEdge(i).numVertex(1);
      int nb_vert_refined = 0;
      if (vertex_to_refine(num(0)))
        nb_vert_refined++;

      if (vertex_to_refine(num(1)))
        nb_vert_refined++;
      
      if (nb_vert_refined == 1)
        {
          nb_new_vert++;
          if (i < mesh.GetNbBoundaryRef())
            nb_new_edges++;
        }

      if (nb_vert_refined == 2)
        {
          if (edge_to_refine(i) || !anisotrope)
            {
              nb_new_vert += 2;
              if (i < mesh.GetNbBoundaryRef())
                nb_new_edges += 2;
            }
        }
    }

  int nb_old_vert = mesh.GetNbVertices();
  int nb_old_edges = mesh.GetNbBoundaryRef();
  int nb_old_elt = mesh.GetNbElt();
  mesh.ResizeVertices(nb_old_vert + nb_new_vert);
  mesh.ResizeBoundariesRef(nb_old_edges + nb_new_edges);
  mesh.ResizeElements(nb_old_elt + nb_new_elt);  

  // interpolation operator for interior points
  Vector<Real_wp> step_subdiv(2); Vector<VectR2> points_div(2);
  step_subdiv(0) = Real_wp(1)/3; step_subdiv(1) = Real_wp(2)/3;
  points_div(1).Reallocate(4);
  points_div(1)(0).Init(step_subdiv(0), step_subdiv(0));
  points_div(1)(1).Init(step_subdiv(1), step_subdiv(0));
  points_div(1)(2).Init(step_subdiv(1), step_subdiv(1));
  points_div(1)(3).Init(step_subdiv(0), step_subdiv(1));
  
  FiniteElementInterpolator interp;    
  interp.InitProjection(mesh.GetReferenceElement(), step_subdiv, points_div);

  // on rajoute les points sur les aretes et les aretes de bord
  Vector<int> offset_edge(mesh.GetNbEdges()); offset_edge.Fill(-1);
  nb_new_vert = nb_old_vert; nb_new_edges = nb_old_edges;
  nb_new_elt = nb_old_elt;
  Real_wp coef1 = Real_wp(2)/3, coef2 = Real_wp(1)-coef1;
  SetPoints<Dimension2> PointsElem;
  for (int i = 0; i < mesh.GetNbEdges(); i++)
    {
      num(0) = mesh.GetEdge(i).numVertex(0);
      num(1) = mesh.GetEdge(i).numVertex(1);
      int nb_vert_refined = 0;
      if (vertex_to_refine(num(0)))
        nb_vert_refined++;

      if (vertex_to_refine(num(1)))
        nb_vert_refined++;

      if (nb_vert_refined == 1)
        {
          int ext0 = num(0), ext1 = num(1);
          if (vertex_to_refine(num(1)))
            {
              ext0 = num(1); ext1 = num(0); 
            }

          if (i < nb_old_edges)
            {
              int ref = mesh.BoundaryRef(i).GetReference();
              mesh.BoundaryRef(i).Init(ext0, nb_new_vert, ref);
              mesh.BoundaryRef(nb_new_edges).Init(ext1, nb_new_vert, ref);
              nb_new_edges++;
            }
          
          offset_edge(i) = nb_new_vert;
          mesh.Vertex(nb_new_vert) = coef1*mesh.Vertex(ext0) + coef2*mesh.Vertex(ext1);          

          nb_new_vert++;          
        }
      
      if (nb_vert_refined == 2)
        {
          offset_edge(i) = nb_new_vert;

          if (edge_to_refine(i) || !anisotrope)
            {
              
              if (i < nb_old_edges)
                {
                  int ref = mesh.BoundaryRef(i).GetReference();
                  mesh.BoundaryRef(i).Init(num(0), nb_new_vert, ref);
                  mesh.BoundaryRef(nb_new_edges).Init(num(1), nb_new_vert+1, ref);
                  mesh.BoundaryRef(nb_new_edges+1).Init(nb_new_vert, nb_new_vert+1, ref);
                  nb_new_edges += 2;
                }
              
              mesh.Vertex(nb_new_vert) = coef1*mesh.Vertex(num(0)) + coef2*mesh.Vertex(num(1));
              mesh.Vertex(nb_new_vert+1) = coef1*mesh.Vertex(num(1)) + coef2*mesh.Vertex(num(0));
              nb_new_vert += 2;
            }
        }      
    }

  // on rajoute les points des elements et les nouveaux elements
  TinyMatrix<int, General, 4, 4> nv, nv2; VectR2 s, InterpPoints;
  for (int i = 0; i < nb_old_elt; i++)
    {
      num(0) = mesh.Element(i).numVertex(0);
      num(1) = mesh.Element(i).numVertex(1);
      num(2) = mesh.Element(i).numVertex(2);
      num(3) = mesh.Element(i).numVertex(3);
      int nb_vert_refined = 0;
      for (int j = 0; j < 4; j++)
        if (vertex_to_refine(num(j)))
          nb_vert_refined++;

      int nb_edge_refined = 0;
      for (int j = 0; j < 4; j++)
        if (edge_to_refine(mesh.Element(i).numEdge(j)))
          nb_edge_refined++;
      
      if (nb_vert_refined == 0)
        continue;

      nv.Fill(-1);
      nv(0, 0) = num(0); nv(3, 0) = num(1); nv(3, 3) = num(2); nv(0, 3) = num(3);
      // numeros sur la premiere arete
      int ne = mesh.Element(i).numEdge(0);
      int ne_refined = 0;
      if (vertex_to_refine(num(0)))
        ne_refined++;
      
      if (vertex_to_refine(num(1)))
        ne_refined++;
      
      if (ne_refined == 2)
        {
          if (num(0) < num(1))
            {
              nv(1, 0) = offset_edge(ne);
              nv(2, 0) = offset_edge(ne)+1;
            }
          else
            {
              nv(1, 0) = offset_edge(ne)+1;
              nv(2, 0) = offset_edge(ne);
            }
        }
      else if (ne_refined == 1)
        {
          if (vertex_to_refine(num(0)))
            nv(1, 0) = offset_edge(ne);
          else
            nv(2, 0) = offset_edge(ne);
        }

      // numeros sur la seconde arete
      ne = mesh.Element(i).numEdge(1);
      ne_refined = 0;
      if (vertex_to_refine(num(1)))
        ne_refined++;
      
      if (vertex_to_refine(num(2)))
        ne_refined++;
      
      if (ne_refined == 2)
        {
          if (num(1) < num(2))
            {
              nv(3, 1) = offset_edge(ne);
              nv(3, 2) = offset_edge(ne)+1;
            }
          else
            {
              nv(3, 1) = offset_edge(ne)+1;
              nv(3, 2) = offset_edge(ne);
            }
        }
      else if (ne_refined == 1)
        {
          if (vertex_to_refine(num(1)))
            nv(3, 1) = offset_edge(ne);
          else
            nv(3, 2) = offset_edge(ne);
        }

      // numeros sur la troisieme arete
      ne = mesh.Element(i).numEdge(2);
      ne_refined = 0;
      if (vertex_to_refine(num(2)))
        ne_refined++;
      
      if (vertex_to_refine(num(3)))
        ne_refined++;
      
      if (ne_refined == 2)
        {
          if (num(3) < num(2))
            {
              nv(1, 3) = offset_edge(ne);
              nv(2, 3) = offset_edge(ne)+1;
            }
          else
            {
              nv(1, 3) = offset_edge(ne)+1;
              nv(2, 3) = offset_edge(ne);
            }
        }
      else if (ne_refined == 1)
        {
          if (vertex_to_refine(num(3)))
            nv(1, 3) = offset_edge(ne);
          else
            nv(2, 3) = offset_edge(ne);
        }

      // numeros sur la quatrieme arete
      ne = mesh.Element(i).numEdge(3);
      ne_refined = 0;
      if (vertex_to_refine(num(3)))
        ne_refined++;
      
      if (vertex_to_refine(num(0)))
        ne_refined++;
      
      if (ne_refined == 2)
        {
          if (num(0) < num(3))
            {
              nv(0, 1) = offset_edge(ne);
              nv(0, 2) = offset_edge(ne)+1;
            }
          else
            {
              nv(0, 1) = offset_edge(ne)+1;
              nv(0, 2) = offset_edge(ne);
            }
        }
      else if (ne_refined == 1)
        {
          if (vertex_to_refine(num(0)))
            nv(0, 1) = offset_edge(ne);
          else
            nv(0, 2) = offset_edge(ne);
        }

      // points interieurs
      int n_first = 0; char type = 'a';
      mesh.GetVerticesElement(i, s);
      mesh.FjElemNodal(s, PointsElem, mesh, i);

      interp.Project(PointsElem.GetPointNodal(), InterpPoints, mesh.GetTypeElement(i));

      if (nb_vert_refined == 1)
        {
          if (vertex_to_refine(num(0)))
            nv(1, 1) = nb_new_vert;
          
          if (vertex_to_refine(num(1)))
            {
              n_first = 1;
              nv(2, 1) = nb_new_vert;
            }
          
          if (vertex_to_refine(num(2)))
            {
              n_first = 2;
              nv(2, 2) = nb_new_vert;
            }
          
          if (vertex_to_refine(num(3)))
            {
              n_first = 3;
              nv(1, 2) = nb_new_vert;
            }
          
          mesh.Vertex(nb_new_vert) = InterpPoints(n_first);
          nb_new_vert++;
        }
      else
        {
          bool create_new_points = true;
          if (anisotrope)
            {
              if ((nb_vert_refined == 2) && (nb_edge_refined == 0))
                create_new_points = false;

              if ((nb_vert_refined == 4) && (nb_edge_refined == 2))
                create_new_points = false;
            }
          
          if (create_new_points)
            {
              mesh.Vertex(nb_new_vert) = InterpPoints(0);
              mesh.Vertex(nb_new_vert+1) = InterpPoints(1);
              mesh.Vertex(nb_new_vert+2) = InterpPoints(2);
              mesh.Vertex(nb_new_vert+3) = InterpPoints(3);
            }
          
          nv(1, 1) = nb_new_vert; nv(2, 1) = nb_new_vert+1;
          nv(2, 2) = nb_new_vert+2; nv(1, 2) = nb_new_vert+3;

          if (nb_vert_refined == 2)
            {
              if (vertex_to_refine(num(0)))
                {
                  if (vertex_to_refine(num(2)))
                    type = 'b';
                  else if (vertex_to_refine(num(3)))
                    n_first = 3;
                }
              else if (vertex_to_refine(num(1)))
                {
                  if (vertex_to_refine(num(2)))
                    n_first = 1;
                  else if (vertex_to_refine(num(3)))
                    {
                      n_first = 1;
                      type = 'b';
                    }
                }
              else if (vertex_to_refine(num(2)))
                {
                  n_first = 2;
                }
            }
          else if (nb_vert_refined == 3)
            {
              if (!vertex_to_refine(num(0)))
                n_first = 1;

              if (!vertex_to_refine(num(1)))
                n_first = 2;

              if (!vertex_to_refine(num(2)))
                n_first = 3;
            }
          else if (nb_vert_refined == 4)
            {
              if (anisotrope && (nb_edge_refined == 2))
                {
                  if (edge_to_refine(mesh.Element(i).numEdge(1)))
                    n_first = 1;
                }
            }
          
          if (create_new_points)
            nb_new_vert += 4;
        }

      // on tourne les numeros si besoin
      if (n_first > 0)
        {
          nv2 = nv;
          if (n_first == 1)
            {
              for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                  nv(j, k) = nv2(3-k, j);
            }
          else if (n_first == 2)
            {
              for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                  nv(j, k) = nv2(3-j, 3-k);
            }
          else
            {
              for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                  nv(j, k) = nv2(k, 3-j);
            }            
        }
      
      // on rajoute les nouveaux elements
      int ref = mesh.Element(i).GetReference();
      if (nb_vert_refined == 1)
        {
          mesh.Element(i).InitQuadrangular(nv(0, 0), nv(1, 0), nv(1, 1), nv(0, 1), ref);
          mesh.Element(nb_new_elt).InitQuadrangular(nv(1, 0), nv(3, 0), nv(3, 3), nv(1, 1), ref);
          mesh.Element(nb_new_elt+1).InitQuadrangular(nv(0, 1), nv(1, 1), nv(3, 3), nv(0, 3), ref);
          nb_new_elt += 2;
        }
      else if (nb_vert_refined == 2)
        {
          if (anisotrope && nb_edge_refined == 0)
            {
              mesh.Element(i).InitQuadrangular(nv(0, 0), nv(3, 0), nv(3, 1), nv(0, 1), ref);
              mesh.Element(nb_new_elt).InitQuadrangular(nv(0, 1), nv(3, 1), nv(3, 3), nv(0, 3), ref);
              nb_new_elt++;
              continue;
            }
          
          if (type == 'a')
            {
              mesh.Element(i).InitQuadrangular(nv(0, 0), nv(1, 0), nv(1, 1), nv(0, 1), ref);
              mesh.Element(nb_new_elt).InitQuadrangular(nv(1, 0), nv(2, 0), nv(2, 1), nv(1, 1), ref);
              mesh.Element(nb_new_elt+1).InitQuadrangular(nv(2, 0), nv(3, 0), nv(3, 1), nv(2, 1), ref);
              mesh.Element(nb_new_elt+2).InitQuadrangular(nv(0, 1), nv(1, 1), nv(1, 2), nv(0, 3), ref);
              mesh.Element(nb_new_elt+3).InitQuadrangular(nv(1, 1), nv(2, 1), nv(2, 2), nv(1, 2), ref);
              mesh.Element(nb_new_elt+4).InitQuadrangular(nv(2, 1), nv(3, 1), nv(3, 3), nv(2, 2), ref);
              mesh.Element(nb_new_elt+5).InitQuadrangular(nv(1, 2), nv(2, 2), nv(3, 3), nv(0, 3), ref);
            }
          else
            {
              mesh.Element(i).InitQuadrangular(nv(0, 0), nv(1, 0), nv(1, 1), nv(0, 1), ref);
              mesh.Element(nb_new_elt).InitQuadrangular(nv(1, 0), nv(3, 0), nv(2, 1), nv(1, 1), ref);
              mesh.Element(nb_new_elt+1).InitQuadrangular(nv(0, 1), nv(1, 1), nv(1, 2), nv(0, 3), ref);
              mesh.Element(nb_new_elt+2).InitQuadrangular(nv(1, 1), nv(2, 1), nv(2, 2), nv(1, 2), ref);
              mesh.Element(nb_new_elt+3).InitQuadrangular(nv(2, 1), nv(3, 0), nv(3, 2), nv(2, 2), ref);
              mesh.Element(nb_new_elt+4).InitQuadrangular(nv(1, 2), nv(2, 2), nv(2, 3), nv(0, 3), ref);
              mesh.Element(nb_new_elt+5).InitQuadrangular(nv(2, 2), nv(3, 2), nv(3, 3), nv(2, 3), ref);
            }

          nb_new_elt += 6;
        }
      else if (nb_vert_refined == 3)
        {
          mesh.Element(i).InitQuadrangular(nv(0, 0), nv(1, 0), nv(1, 1), nv(0, 1), ref);
          mesh.Element(nb_new_elt).InitQuadrangular(nv(1, 0), nv(2, 0), nv(2, 1), nv(1, 1), ref);
          mesh.Element(nb_new_elt+1).InitQuadrangular(nv(2, 0), nv(3, 0), nv(3, 1), nv(2, 1), ref);
          mesh.Element(nb_new_elt+2).InitQuadrangular(nv(0, 1), nv(1, 1), nv(1, 2), nv(0, 3), ref);
          mesh.Element(nb_new_elt+3).InitQuadrangular(nv(1, 1), nv(2, 1), nv(2, 2), nv(1, 2), ref);
          mesh.Element(nb_new_elt+4).InitQuadrangular(nv(2, 1), nv(3, 1), nv(3, 2), nv(2, 2), ref);
          mesh.Element(nb_new_elt+5).InitQuadrangular(nv(1, 2), nv(2, 2), nv(2, 3), nv(0, 3), ref);
          mesh.Element(nb_new_elt+6).InitQuadrangular(nv(2, 2), nv(3, 2), nv(3, 3), nv(2, 3), ref);
          nb_new_elt += 7;
        }
      else
        {
          if (anisotrope && nb_edge_refined == 2)
            {
              mesh.Element(i).InitQuadrangular(nv(0, 0), nv(1, 0), nv(1, 3), nv(0, 3), ref);
              mesh.Element(nb_new_elt).InitQuadrangular(nv(1, 0), nv(2, 0), nv(2, 3), nv(1, 3), ref);
              mesh.Element(nb_new_elt+1).InitQuadrangular(nv(2, 0), nv(3, 0), nv(3, 3), nv(2, 3), ref);
              nb_new_elt += 2;
              continue;
            }

          mesh.Element(i).InitQuadrangular(nv(0, 0), nv(1, 0), nv(1, 1), nv(0, 1), ref);
          mesh.Element(nb_new_elt).InitQuadrangular(nv(1, 0), nv(2, 0), nv(2, 1), nv(1, 1), ref);
          mesh.Element(nb_new_elt+1).InitQuadrangular(nv(2, 0), nv(3, 0), nv(3, 1), nv(2, 1), ref);
          mesh.Element(nb_new_elt+2).InitQuadrangular(nv(0, 1), nv(1, 1), nv(1, 2), nv(0, 2), ref);
          mesh.Element(nb_new_elt+3).InitQuadrangular(nv(1, 1), nv(2, 1), nv(2, 2), nv(1, 2), ref);
          mesh.Element(nb_new_elt+4).InitQuadrangular(nv(2, 1), nv(3, 1), nv(3, 2), nv(2, 2), ref);
          mesh.Element(nb_new_elt+5).InitQuadrangular(nv(0, 2), nv(1, 2), nv(1, 3), nv(0, 3), ref);
          mesh.Element(nb_new_elt+6).InitQuadrangular(nv(1, 2), nv(2, 2), nv(2, 3), nv(1, 3), ref);
          mesh.Element(nb_new_elt+7).InitQuadrangular(nv(2, 2), nv(3, 2), nv(3, 3), nv(2, 3), ref);
          nb_new_elt += 8;
        }
    }

  mesh.ReorientElements();
  mesh.FindConnectivity();
}
  
void RefineLocallyMesh(Mesh<Dimension2>& mesh, int ref, const Real_wp& h, bool anisotrope)
{
  bool test_loop = true;
  while (test_loop)
    {
      test_loop = false;
      Vector<bool> vertex_to_refine(mesh.GetNbVertices());
      Vector<bool> edge_to_refine(mesh.GetNbEdges());
      vertex_to_refine.Fill(false);
      edge_to_refine.Fill(false);
      for (int i = 0; i < mesh.GetNbEdges(); i++)
        {
          int e1 = mesh.GetEdge(i).numElement(0);
          int ref1 = mesh.Element(e1).GetReference(), ref2 = -1;
          if (mesh.GetEdge(i).GetNbElements() > 1)
            {
              int e2 = mesh.GetEdge(i).numElement(1);
              ref2 = mesh.Element(e2).GetReference();
            }
          
          if ((ref1 == ref) || (ref2 == ref))
            {
              // on calcule la longueur de l'arete et on raffine si besoin
              int n0 = mesh.GetEdge(i).numVertex(0), n1 = mesh.GetEdge(i).numVertex(1);
              Real_wp hi = mesh.Vertex(n0).Distance(mesh.Vertex(n1));
              if (hi > h)
                {
                  test_loop = true;
                  edge_to_refine(i) = true;
                  vertex_to_refine(n0) = true;
                  vertex_to_refine(n1) = true;
                }
            }
        }

      // si pas anisotrope, on raffine les aretes avec deux sommets raffines
      if (!anisotrope)
        for (int i = 0; i < mesh.GetNbEdges(); i++)
          {
            int n0 = mesh.GetEdge(i).numVertex(0), n1 = mesh.GetEdge(i).numVertex(1);
            if (vertex_to_refine(n0) && vertex_to_refine(n1))
              edge_to_refine(i) = true;
          }
      
      RefineLocallyMesh(mesh, vertex_to_refine, edge_to_refine, anisotrope);
      mesh.Write("refined.mesh");
    }
}

void ExtrudeOriginLayer(const R2& center, const Real_wp& radius, int ref,
                        int N, bool circle, const R2& ptMin, const R2& ptMax,
                        int ref_domain, Mesh<Dimension2>& mesh, int new_ref = -1)
{
  Vector<int> IndexVertex(mesh.GetNbVertices());
  IndexVertex.Fill(-1);
  int nb_new_edge = 0;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    if (mesh.BoundaryRef(i).GetReference() == ref)
      {
        int n1 = mesh.BoundaryRef(i).numVertex(0);
        int n2 = mesh.BoundaryRef(i).numVertex(1);
        IndexVertex(n1) = 0;
        IndexVertex(n2) = 0;
        nb_new_edge++;
      }

  int nb_new_vert = 0;
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    if (IndexVertex(i) == 0)
      IndexVertex(i) = nb_new_vert++;

  VectR2 vert_ext(nb_new_vert);
  VectReal_wp pos_left, pos_right, pos_top, pos_bottom;
  Vector<int> num_left, num_right, num_top, num_bottom;
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    if (IndexVertex(i) >= 0)
      {
        R2 pt = mesh.Vertex(i) - center;
        if (circle)
          {
            Mlt(radius / Norm2(pt), pt);
            vert_ext(IndexVertex(i)) = center + pt;
          }
        else
          {
            if ((pt(1) >= pt(0)-R2::threshold) && (pt(1) >= -pt(0)-R2::threshold))
              {
                pos_top.PushBack(pt(0));
                num_top.PushBack(i);
              }
            
            if ((pt(1) <= pt(0)+R2::threshold) && (pt(1) >= -pt(0)-R2::threshold))
              {
                pos_right.PushBack(pt(1));
                num_right.PushBack(i);
              }
            
            if ((pt(1) <= pt(0)+R2::threshold) && (pt(1) <= -pt(0)+R2::threshold))
              {
                pos_bottom.PushBack(pt(0));
                num_bottom.PushBack(i);
              }
            
            if ((pt(1) >= pt(0)-R2::threshold) && (pt(1) <= -pt(0)+R2::threshold))
              {
                pos_left.PushBack(pt(1));
                num_left.PushBack(i);
              }
          }
      }

  Sort(pos_top, num_top);
  Sort(pos_bottom, num_bottom);
  Sort(pos_left, num_left);
  Sort(pos_right, num_right);

  Real_wp dx = (ptMax(0) - ptMin(0)) / (num_top.GetM()-1);
  for (int i = 0; i < num_top.GetM(); i++)
    {
      int n = IndexVertex(num_top(i));
      if (i == 0)
        vert_ext(n).Init(ptMin(0), ptMax(1));
      else if (i == num_top.GetM()-1)
        vert_ext(n).Init(ptMax(0), ptMax(1));
      else
        vert_ext(n).Init(ptMin(0) + i*dx, ptMax(1));
    }

  for (int i = 0; i < num_bottom.GetM(); i++)
    {
      int n = IndexVertex(num_bottom(i));
      if (i == 0)
        vert_ext(n).Init(ptMin(0), ptMin(1));
      else if (i == num_bottom.GetM()-1)
        vert_ext(n).Init(ptMax(0), ptMin(1));
      else
        vert_ext(n).Init(ptMin(0) + i*dx, ptMin(1));
    }

  Real_wp dy = (ptMax(1) - ptMin(1))/(num_left.GetM()-1);
  for (int i = 0; i < num_left.GetM(); i++)
    {
      int n = IndexVertex(num_left(i));
      if (i == 0)
        vert_ext(n).Init(ptMin(0), ptMin(1));
      else if (i == num_left.GetM()-1)
        vert_ext(n).Init(ptMin(0), ptMax(1));
      else
        vert_ext(n).Init(ptMin(0), ptMin(1) + i*dy);
    }

  for (int i = 0; i < num_right.GetM(); i++)
    {
      int n = IndexVertex(num_right(i));
      if (i == 0)
        vert_ext(n).Init(ptMax(0), ptMin(1));
      else if (i == num_left.GetM()-1)
        vert_ext(n).Init(ptMax(0), ptMax(1));
      else
        vert_ext(n).Init(ptMax(0), ptMin(1) + i*dy);
    }
  
  int nb_old_vert = mesh.GetNbVertices();
  int nb_old_elt = mesh.GetNbElt();
  mesh.ResizeVertices(nb_old_vert + nb_new_vert*N);
  mesh.ResizeElements(nb_old_elt + nb_new_edge*N);
  for (int k = 0; k < N; k++)
    {
      Real_wp la(1);
      if (k < N-1)
        la = Real_wp(k+1) / Real_wp(N);
      
      for (int i = 0; i < nb_old_vert; i++)
        if (IndexVertex(i) >= 0)
          {
            int n = IndexVertex(i);
            R2 pt = (Real_wp(1) - la)*mesh.Vertex(i) + la*vert_ext(n);
            mesh.Vertex(nb_old_vert + k*nb_new_vert + n) = pt;            
          }

      int offset = 0;
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
        if (mesh.BoundaryRef(i).GetReference() == ref)
          {
            int p1 = mesh.BoundaryRef(i).numVertex(0);
            int p2 = mesh.BoundaryRef(i).numVertex(1);
            int n3 = nb_old_vert + k*nb_new_vert + IndexVertex(p2);
            int n4 = nb_old_vert + k*nb_new_vert + IndexVertex(p1);
            int n1 = p1, n2 = p2;
            if (k > 0)
              {
                n1 = n4 - nb_new_vert;
                n2 = n3 - nb_new_vert;
              }

            mesh.Element(nb_old_elt + k*nb_new_edge + offset).InitQuadrangular(n1, n2, n3, n4, ref_domain);
              offset++;
          }
    }
  
  if (new_ref > 0)
    {
      int nb_old_ref = mesh.GetNbBoundaryRef();
      mesh.ResizeBoundariesRef(nb_old_ref + nb_new_edge);
      int offset = nb_old_ref;
      for (int i = 0; i < nb_old_ref; i++)
        if (mesh.BoundaryRef(i).GetReference() == ref)
          {
            int p1 = mesh.BoundaryRef(i).numVertex(0);
            int p2 = mesh.BoundaryRef(i).numVertex(1);
            int n1 = nb_old_vert + (N-1)*nb_new_vert + IndexVertex(p1);
            int n2 = nb_old_vert + (N-1)*nb_new_vert + IndexVertex(p2);
            mesh.BoundaryRef(offset).Init(n1, n2, new_ref);
            offset++;
          }
    }
  else
    {
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
        if (mesh.BoundaryRef(i).GetReference() == ref)
          {
            int p1 = mesh.BoundaryRef(i).numVertex(0);
            int p2 = mesh.BoundaryRef(i).numVertex(1);
            int n1 = nb_old_vert + (N-1)*nb_new_vert + IndexVertex(p1);
            int n2 = nb_old_vert + (N-1)*nb_new_vert + IndexVertex(p2);
            mesh.BoundaryRef(i).Init(n1, n2, ref);
          }
    }

  mesh.ReorientElements();
  mesh.FindConnectivity();
  
}

void CreateCircularBoundary(const R2& center, const Real_wp& radius, int N,
                            int ref, Mesh<Dimension2>& mesh)
{
  mesh.ReallocateVertices(4*N);
  mesh.ReallocateBoundariesRef(4*N);
  for (int i = 0; i < N; i++)
    {
      Real_wp teta = pi_wp/4 + Real_wp(i)*pi_wp/(2*N);
      Real_wp x0 = radius*cos(teta), y0 = radius*sin(teta);
      mesh.Vertex(i).Init(center(0) + x0, center(1) + y0);
      mesh.Vertex(i+N).Init(center(0) - y0, center(1) + x0);
      mesh.Vertex(i+2*N).Init(center(0) - x0, center(1) - y0);
      mesh.Vertex(i+3*N).Init(center(0) + y0, center(1) - x0);

      mesh.BoundaryRef(i).Init(i, i+1, ref);
      mesh.BoundaryRef(i+N).Init(i+N, i+N+1, ref);
      mesh.BoundaryRef(i+2*N).Init(i+2*N, i+2*N+1, ref);
      if (i == N-1)
        mesh.BoundaryRef(i+3*N).Init(0, i+3*N, ref);
      else
        mesh.BoundaryRef(i+3*N).Init(i+3*N, i+3*N+1, ref);
    }  
}

void CreateCircularMesh(const R2& center, const Real_wp& radius, int N,
                        int ref, Mesh<Dimension2>& mesh)
{
  if (N%2 == 0)
    {
      mesh.ReallocateVertices(9);
      mesh.Vertex(0) = center;
      mesh.Vertex(1).Init(center(0) + 0.7*radius, center(1));
      mesh.Vertex(2).Init(center(0) + 0.5*radius, center(1) + 0.5*radius);
      mesh.Vertex(3).Init(center(0), center(1) + 0.7*radius);
      mesh.Vertex(4).Init(center(0) - 0.5*radius, center(1) + 0.5*radius);
      mesh.Vertex(5).Init(center(0) - 0.7*radius, center(1));
      mesh.Vertex(6).Init(center(0) - 0.5*radius, center(1)-0.5*radius);
      mesh.Vertex(7).Init(center(0), center(1)  - 0.7*radius);
      mesh.Vertex(8).Init(center(0) + 0.5*radius, center(1) - 0.5*radius);

      mesh.ReallocateElements(4);
      mesh.Element(0).InitQuadrangular(0, 1, 2, 3, 2);
      mesh.Element(1).InitQuadrangular(3, 4, 5, 0, 2);
      mesh.Element(2).InitQuadrangular(5, 6, 7, 0, 2);
      mesh.Element(3).InitQuadrangular(7, 8, 1, 0, 2);

      mesh.ReallocateBoundariesRef(8);
      mesh.BoundaryRef(0).Init(1, 2, ref);
      mesh.BoundaryRef(1).Init(2, 3, ref);
      mesh.BoundaryRef(2).Init(3, 4, ref);
      mesh.BoundaryRef(3).Init(4, 5, ref);
      mesh.BoundaryRef(4).Init(5, 6, ref);
      mesh.BoundaryRef(5).Init(6, 7, ref);
      mesh.BoundaryRef(6).Init(7, 8, ref);
      mesh.BoundaryRef(7).Init(1, 8, ref);

      mesh.ResizeNbReferences(ref+1);
      mesh.ReorientElements();
      mesh.FindConnectivity();
      
      if (N > 2)
        {
          VectReal_wp step_subdiv(N/2 + 1);
          step_subdiv.Fill(); Mlt(Real_wp(2)/Real_wp(N), step_subdiv);
          mesh.SubdivideMesh(step_subdiv);
        }

      int Nc = max(1, N/3);
      ExtrudeOriginLayer(center, radius, ref, Nc, true, R2(), R2(), 2, mesh);
    }
  else
    {
      mesh.ReallocateVertices(4);
      mesh.Vertex(0).Init(center(0) + 0.5*radius, center(1) + 0.5*radius);
      mesh.Vertex(1).Init(center(0) - 0.5*radius, center(1) + 0.5*radius);
      mesh.Vertex(2).Init(center(0) - 0.5*radius, center(1) - 0.5*radius);
      mesh.Vertex(3).Init(center(0) + 0.5*radius, center(1) - 0.5*radius);

      mesh.ReallocateElements(1);
      mesh.Element(0).InitQuadrangular(0, 1, 2, 3, 2);

      mesh.ReallocateBoundariesRef(4);
      mesh.BoundaryRef(0).Init(0, 1, ref);
      mesh.BoundaryRef(1).Init(1, 2, ref);
      mesh.BoundaryRef(2).Init(2, 3, ref);
      mesh.BoundaryRef(3).Init(0, 3, ref);

      mesh.ResizeNbReferences(ref+1);
      mesh.ReorientElements();
      mesh.FindConnectivity();
      
      if (N > 1)
        {
          VectReal_wp step_subdiv((N+1)/2 + 1);
          step_subdiv.Fill(); Mlt(Real_wp(2)/Real_wp(N+1), step_subdiv);
          mesh.SubdivideMesh(step_subdiv);
        }

      int Nc = max(1, (N+1)/4);
      ExtrudeOriginLayer(center, radius, ref, Nc, true, R2(), R2(), 2, mesh);
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  Mesh<Dimension2> mesh;
  int order_geom = 1;
  string output_name("toto.mesh");
  if (argc > 1)
    {      
      string nom_init(argv[1]);
      int dimension = GetDimensionMesh(nom_init, GetExtension(nom_init));
      if (dimension == 3)
        {
          cout << " Only two-dimensional mesh are handled " << endl;
          abort();
        }
      
      if (argc > 2)
        {
          order_geom = atoi(argv[2]);
          mesh.SetGeometryOrder(order_geom);
	  if (order_geom > 1)
	    output_name = string("toto.msh");
          
	  if (argc > 3)
            output_name = string(argv[3]);
        }
      
      // lecture du maillage
      mesh.Read(nom_init);
    }
  
  // boucle de manipulation
  bool test_loop = true;
  IVect num, num2; VectReal_wp coord;
  while (test_loop)
    {
      cout << "1 - AJOUTER sommet, arete, element, maillage"<<endl;
      cout << "2 - ENLEVER sommet, arete, element"<<endl;
      cout << "3 - MODIFIER sommet, arete, element, references "<<endl;
      cout << "4 - DECOUPER"<<endl;      
      cout << "5 - VERIFIER / INFORMER"<<endl;      
      cout << "6 - SAUVER"<<endl;
      cout << "7 - FINIR"<<endl;
      cout << "8 - METTRE A JOUR CONNECTIVITE"<<endl;
      cout << "9- EXTRAIRE sommet, arete, element" << endl;
      
      ReadParameters(num); int choix = num(0);
      
      switch (choix)
	{
	case 1 :
	  {
	    // menu pour ajouter
	    cout<<"1- Sommet"<<endl;
	    cout<<"2- Arete"<<endl;
	    cout<<"3- Element"<<endl;
	    cout<<"4- Maillage"<<endl;
	    cout<<"5- Zone PML" << endl;
	    cout<<"Autre- Retour au menu principal" << endl;
	    ReadParameters(num); int choix_ajout = num(0);
	    Real_wp x, y; int ref;
	    string nom_maillage; 
	    switch (choix_ajout)
	      {
	      case 1 :
		cout<<"Entrez les coordonnees x, y du sommet"<<endl;
		ReadParameters(coord);
		x = coord(0); y = coord(1);
		AddVertex(mesh, R2(x, y));
		break;
	      case 2 :
		cout<<"Entrez les numeros des extremites de l'arete"<<endl;
		ReadParameters(num);
		cout<<"Entrez la reference de cette arete "<<endl;
		ReadParameters(num2); ref = num2(0);
		AddEdge(mesh, num(0), num(1), ref);
		break;
	      case 3 :
		cout<<"Entrez les numeros des sommets de l'element"<<endl;
		ReadParameters(num);
		cout<<"Entrez la reference de cet element "<<endl;
		ReadParameters(num2); ref = num2(0);
		AddElement(mesh, num, ref);
		break;
	      case 4 :
                {
                  cout << "1- Rajouter un rectangle regulier " << endl;
                  cout << "2- Rajouter un maillage regulier d'un carre avec trou circulaire " << endl;
                  cout << "3- Rajouter un maillage regulier d'un carre avec un dielectrique circulaire" << endl;
                  cout << "4- Rajouter un maillage d'un fichier " << endl;
                  cout << "5- Rajouter un rectangle regulier avec progression geometrique" << endl;
                  cout << "6- Rajouter un maillage d'un carre avec raffinement au centre" << endl;
                  cout << "7- Rajouter un maillage hexagonal" << endl;
		  cout<<"Autre- Retour au menu principal" << endl;
                  ReadParameters(num); int choix_maillage = num(0);
                  Mesh<Dimension2> new_mesh;
                  new_mesh.SetGeometryOrder(order_geom);
                  switch (choix_maillage)
                    {
                    case 1 :
                      {
                        cout << "Donnez xmin, xmax, ymin, ymax " << endl;
                        Real_wp xmin, xmax, ymin, ymax;
			ReadParameters(coord);
			xmin = coord(0); xmax = coord(1); ymin = coord(2); ymax = coord(3);
			cout << "Donnez le nombre de points suivant x et y " << endl;
                        ReadParameters(num);
			int nbx = num(0), nby = num(1);
                        cout << "Reference du domaine " << endl;
                        ReadParameters(num); int ref_domain = num(0);
                        cout << "References des bords y = ymin, x = xmax, y = ymax, x = xmin " << endl;
                        TinyVector<int, 4> ref_boundary;
			ReadParameters(num);
			for (int i = 0; i < 4; i++)
			  ref_boundary(i) = num(i);
			
                        cout << "Type de maillage " << endl;
                        cout << "0- Triangles " << endl;
                        cout << "1- Quadrangles " << endl;
                        ReadParameters(num); int type_mesh = num(0);
                        new_mesh.CreateRegularMesh(R2(xmin, ymin), R2(xmax, ymax), TinyVector<int, 2>(nbx, nby),
                                                   ref_domain, ref_boundary, type_mesh);
			
                      }
                      break;
                    case 2 :
                      {
                        cout << "Donnez le centre du cercle et le rayon" << endl;
                        ReadParameters(coord);
                        Real_wp x0 = coord(0), y0 = coord(1), radius = coord(2);
                        cout << "Entrez le nombre de subdivisions pour un quart de cercle et sur la tranche et la reference du cercle" << endl;
                        ReadParameters(num);
                        int N = num(0), Nt = num(1), ref = num(2);
                        CreateCircularBoundary(R2(x0, y0), radius, N, ref, new_mesh);
                        cout << "Entrez xmin, xmax, ymin, ymax de la boite" << endl;
                        ReadParameters(coord);
                        R2 ptMin(coord(0), coord(2));
                        R2 ptMax(coord(1), coord(3));
                        cout << "Entrez la reference de la frontiere externe" << endl;
                        ReadParameters(num);
                        int new_ref = num(0);
                        ExtrudeOriginLayer(R2(x0, y0), radius, ref, Nt, false,
                                           ptMin, ptMax, 1, new_mesh, new_ref);
                      }
                      break;
                    case 3 :
                      {
                        cout << "Donnez le centre du cercle et le rayon" << endl;
                        ReadParameters(coord);
                        Real_wp x0 = coord(0), y0 = coord(1), radius = coord(2);
                        cout << "Entrez le nombre de subdivisions pour un quart de cercle et sur la tranche et la reference du cercle" << endl;
                        ReadParameters(num);
                        int N = num(0), Nt = num(1), ref = num(2);
			CreateCircularMesh(R2(x0, y0), radius, N, ref, new_mesh);
                        cout << "Entrez xmin, xmax, ymin, ymax de la boite" << endl;
                        ReadParameters(coord);
                        R2 ptMin(coord(0), coord(2));
                        R2 ptMax(coord(1), coord(3));
                        cout << "Entrez la reference de la frontiere externe" << endl;
                        ReadParameters(num);
                        int new_ref = num(0);
                        ExtrudeOriginLayer(R2(x0, y0), radius, ref, Nt, false,
                                           ptMin, ptMax, 1, new_mesh, new_ref);
                      }
                      break;
                    case 4 :
                      {
                        cout << "Donnez le nom du fichier de maillage " << endl;
                        ReadParameters(nom_maillage);
                        new_mesh.Read(nom_maillage);
                      }
                      break;
                    case 5 :
                      {
                        cout << "Donnez xmin, xmax, ymin, ymax " << endl;
			ReadParameters(coord);
                        Real_wp xmin = coord(0), xmax  = coord(1), ymin  = coord(2), ymax  = coord(3);
			cout << "Donnez le nombre de points suivant x et y " << endl;
                        ReadParameters(num); int nbx = num(0), nby = num(1);
                        cout << "Donnez le ratio pour x y " << endl;
                        R2 ratio;
                        ReadParameters(coord);
			ratio(0) = coord(0); ratio(1) = coord(1);
                        cout << "Reference du domaine " << endl;
                        ReadParameters(num); int ref_domain = num(0);
                        cout << "References du bord y = ymin, x = xmax, y = ymax et x = xmin " << endl;
                        TinyVector<int, 4> ref_boundary;
                        ReadParameters(num);
			for (int i = 0; i < 4; i++)
			  ref_boundary(i) = num(i);
                                                
                        cout << "Type de maillage " << endl;
                        cout << "0- Triangles " << endl;
                        cout << "1- Quadrangles " << endl;
                        ReadParameters(num); int type_mesh = num(0);
                        new_mesh.CreateRegularMesh(R2(xmin, ymin), R2(xmax, ymax), TinyVector<int, 2>(nbx, nby),
                                                   ref_domain, ref_boundary, type_mesh, ratio);
                        
                      }
                      break;
                    case 6:
                      {
                        cout << "Donnez xmin, xmax, ymin, ymax " << endl;
			ReadParameters(coord);
                        Real_wp xmin = coord(0), xmax  = coord(1), ymin  = coord(2), ymax  = coord(3);
			cout << "Donnez le pas de maillage " << endl;
                        ReadParameters(coord); Real_wp h = coord(0);
                        cout << "Donnez la zone a raffiner xmin, xmax, ymin, ymax" << endl;
                        ReadParameters(coord);
                        Real_wp x0 = coord(0), xN  = coord(1), y0  = coord(2), yN  = coord(3);
                        
                        cout << "Reference du domaine " << endl;
                        ReadParameters(num); int ref_domain = num(0);
                        cout << "References du bord y = ymin, x = xmax, y = ymax et x = xmin " << endl;
                        TinyVector<int, 4> ref_boundary;
                        ReadParameters(num);
			for (int i = 0; i < 4; i++)
			  ref_boundary(i) = num(i);

                        CreateLocallyRefinedRegularMesh(R2(xmin, ymin), R2(xmax, ymax),
                                                        R2(x0, y0), R2(xN, yN), h, ref_domain,
                                                        ref_boundary, new_mesh);
                      }
                      break;
                    case 7:
                      {
                        cout << "Donnez le centre de la structure, la distance entre chaque trou et le rayon de chaque trou" << endl;
                        ReadParameters(coord);
                        Real_wp x0 = coord(0), y0 = coord(1), dx = coord(2), r = coord(3);
                        Real_wp dy = sqrt(Real_wp(3))*dx/2;
                        Real_wp coef = sqrt(Real_wp(3));
                        cout << "Donnez le nombre de lignes, le nombre de subdivisions pour chaque cellule et sur la tranche" << endl;
                        ReadParameters(num);
                        int Nx = num(0), N = num(1), Nt = num(2);
                        cout << "Donnez le nombre de cellules en plus gauche/droite et haut/bas" << endl;
                        ReadParameters(num);
                        int N1 = num(0), N2 = num(1);
                        cout << "Pour la premiere couronne, donnez la reduction du rayon et decalage" << endl;
                        ReadParameters(coord);
                        Real_wp ratio = coord(0), dec = coord(1);
                        
                        // motif de base
                        Mesh<Dimension2> cell_mesh;
                        CreateCircularMesh(R2(0, 0), r, N, 2, cell_mesh);
                        ExtrudeOriginLayer(R2(0, 0), r, 2, Nt, false,
                                           R2(-dx/2, -dy/2), R2(dx/2, dy/2), 1, cell_mesh, 1);

                        // brique pour completer le maillage
                        Mesh<Dimension2> brique, demi_brique;
                        demi_brique.CreateRegularMesh(R2(0, -dy/2), R2(dx/2, dy/2), TinyVector<int, 2>(N/2+1, N+1), 1, TinyVector<int, 4>(1, 1, 1, 1), brique.QUADRILATERAL_MESH);
                        brique.CreateRegularMesh(R2(-dx/2, -dy/2), R2(dx/2, dy/2), TinyVector<int, 2>(N+1, N+1), 1, TinyVector<int, 4>(1, 1, 1, 1), brique.QUADRILATERAL_MESH);

                        int bascule = 0; int new_ref = 2;
                        new_mesh.ResizeNbReferences((2*Nx+1)*(2*Nx+1));
                        for (int j = -Nx-N2; j <= Nx+N2; j++)
                          {
                            Real_wp yj = y0 + j*dy;
                            Mesh<Dimension2> row_mesh;
                            row_mesh.ResizeNbReferences((2*Nx+1)*(2*Nx+1));
                            if (bascule == 0)
                              row_mesh.AppendMesh(demi_brique, true, R2(x0-(Nx+N1+1)*dx, yj));
                            
                            for (int i = -Nx-N1; i <= Nx+N1; i++)
                              {
                                Real_wp xi = x0 + i*dx;
                                if (bascule == 1)
                                  xi -= 0.5*dx;

                                bool trou = false;
                                if ((i < -Nx) || (i > Nx) || (j < -Nx) || (j > Nx))
                                  trou = true;

                                if ((i == 0) && (j == 0))
                                  trou = true;

                                if ((coef*xi+yj <= -coef*(Nx+0.6)*dx) || (-coef*xi+yj >= coef*(Nx+0.6)*dx) || (coef*xi+yj >= coef*(Nx+0.4)*dx) || (-coef*xi+yj <= -coef*(Nx+0.4)*dx))
                                  trou = true;

                                Real_wp xc = 0, yc = 0, ri = r;
                                bool trou_decale = false;
                                if ((i == -1) && (j == 0))
                                  {
                                    ri *= ratio; trou_decale = true;
                                    xc -= dec*dx;
                                  }
                                else if ((i == 1) && (j == 0))
                                  {
                                    ri *= ratio; trou_decale = true;
                                    xc += dec*dx;
                                  }
                                else if ((i == -1) && (j == 1))
                                  {
                                    ri *= ratio; trou_decale = true;
                                    xc -= dec*dx/2;
                                    yc += dec*dy;
                                  }
                                else if ((i == 0) && (j == 1))
                                  {
                                    ri *= ratio; trou_decale = true;
                                    xc += dec*dx/2;
                                    yc += dec*dy;
                                  }
                                else if ((i == 0) && (j == -1))
                                  {
                                    ri *= ratio; trou_decale = true;
                                    xc += dec*dx/2;
                                    yc -= dec*dy;
                                  }
                                else if ((i == -1) && (j == -1))
                                  {
                                    ri *= ratio; trou_decale = true;
                                    xc -= dec*dx/2;
                                    yc -= dec*dy;
                                  }

                                if (trou)
                                  row_mesh.AppendMesh(brique, true, R2(xi, yj));
                                else if (trou_decale)
                                  {
                                    Mesh<Dimension2> mesh_hole;
                                    CreateCircularMesh(R2(xc, yc), ri, N, new_ref, mesh_hole);
                                    ExtrudeOriginLayer(R2(xc, yc), ri, new_ref, Nt, false,
                                                       R2(-dx/2, -dy/2), R2(dx/2, dy/2), 1, mesh_hole, 1);
                                    row_mesh.AppendMesh(mesh_hole, true, R2(xi, yj));
                                    new_ref++;
                                  }
                                else
                                  {
                                    for (int i = 0; i < cell_mesh.GetNbBoundaryRef(); i++)
                                      if (cell_mesh.BoundaryRef(i).GetReference() == 2)
                                        cell_mesh.BoundaryRef(i).SetReference(new_ref);
                                    
                                    row_mesh.AppendMesh(cell_mesh, true, R2(xi, yj));

                                    for (int i = 0; i < cell_mesh.GetNbBoundaryRef(); i++)
                                      if (cell_mesh.BoundaryRef(i).GetReference() == new_ref)
                                        cell_mesh.BoundaryRef(i).SetReference(2);

                                    new_ref++;
                                  }
                              }

                            if (bascule == 1)
                              row_mesh.AppendMesh(demi_brique, true, R2(x0+(Nx+N1)*dx, yj));
                            
                            new_mesh.AppendMesh(row_mesh);
                            bascule = 1-bascule;
                          }

                        num.Reallocate(1); num(0) = 1;
                        RemoveEdgeRef(new_mesh, num);

                        new_mesh.SortBoundariesRef();
                        new_mesh.FindConnectivity();
                        new_mesh.ReorientElements();
                        new_mesh.AddBoundaryEdges();
                        new_mesh.FindConnectivity();

                        for (int i = 0; i < new_mesh.GetNbBoundaryRef(); i++)
                          if (new_mesh.BoundaryRef(i).GetReference() == new_ref)
                            new_mesh.BoundaryRef(i).SetReference(1);
                      }
                      break;
                    case 8:
                      {
                        cout << "Donnez le nombre de couches suivant x y et z " << endl;
                        ReadParameters(num);
			int nbx = num(0), nby = num(1);      
                        cout << "Donnez les positions x " << endl;
                        VectReal_wp pos_x;
			ReadParameters(pos_x);
                        cout << "Donnez les positions y " << endl;
                        VectReal_wp pos_y;
			ReadParameters(pos_y);
                        cout << "Donnez le nombre de sous-couches en x " << endl;
                        Vector<int> nb_layers_x; ReadParameters(nb_layers_x);
                        if ((nb_layers_x.GetM() != pos_x.GetM()-1) || (nb_layers_x.GetM() != nbx))
                          {
                            cout << "Nombre couches (N) = " << nbx << endl;
                            cout << "Nombre d'entiers = (egal a N) " << nb_layers_x.GetM() << endl;
                            cout << "Nombre de positions = (egal a N+1) " << pos_x.GetM() << endl;
                            cout << "Probleme de coherence" << endl;
                            return FinalizeMontjoie();
                          }
                        
                        cout << "Donnez le nombre de sous-couches en y " << endl;
                        Vector<int> nb_layers_y; ReadParameters(nb_layers_y);
                        if ((nb_layers_y.GetM() != pos_y.GetM()-1) || (nb_layers_y.GetM() != nby))
                          {
                            cout << "Nombre couches (N) = " << nby << endl;
                            cout << "Nombre d'entiers = (egal a N) " << nb_layers_y.GetM() << endl;
                            cout << "Nombre de positions = (egal a N+1) " << pos_y.GetM() << endl;
                            cout << "Probleme de coherence" << endl;
                            return FinalizeMontjoie();
                          }
                                                
                        cout << "Reference du domaine " << endl;
                        ReadParameters(num); int ref_domain = num(0);
                        cout << "References sous-domaine (liste de i, j, ref) " << endl;
                        Matrix<int> ref_domain_layer(nbx, nby);
                        ref_domain_layer.Fill(ref_domain);
                        ReadParameters(num);
                        for (int i = 0; i < num.GetM(); i += 3)
                          ref_domain_layer(num(i), num(i+1)) = num(i+2);
                        
                        // sub-divisions
                        cout << "Souhaitez-vous donner les subdivisions ? 0- Oui, 1-Non" << endl;
                        ReadParameters(num);
                        if (num(0) == 0)
                          {
                            cout << "Donnez les subdivisions pour x" << endl;
                            ReadParameters(pos_x);

                            cout << "Donnez les subdivisions pour y" << endl;
                            ReadParameters(pos_y);

                            Matrix<int> old_ref(ref_domain_layer);
                            IVect nx(nb_layers_x.GetM()+1), ny(nb_layers_y.GetM()+1);
                            nx(0) = 0; ny(0) = 0;
                            // on calcule les indices cumules dans nx/ny/nz
                            for (int i = 0; i < nb_layers_x.GetM(); i++)
                              nx(i+1) = nx(i) + nb_layers_x(i);

                            for (int i = 0; i < nb_layers_y.GetM(); i++)
                              ny(i+1) = ny(i) + nb_layers_y(i);
                            
                            if ((pos_x.GetM() != nx(nb_layers_x.GetM())+1)
                                || (pos_y.GetM() != ny(nb_layers_y.GetM())+1))
                              {
                                cout << "Probleme de coherence" << endl;
                                DISP(nx); DISP(ny); DISP(pos_x.GetM());
                                DISP(pos_y.GetM()); 
                                return FinalizeMontjoie();
                              }
                            
                            // on recalcule ref_domain_layer                            
                            ref_domain_layer.Reallocate(pos_x.GetM()-1, pos_y.GetM()-1);
                            ref_domain_layer.Fill(ref_domain);
                            for (int i = 0; i < nbx; i++)
                              for (int j = 0; j < nby; j++)
                                if (old_ref(i, j) != ref_domain)
                                  {
                                    for (int i2 = 0; i2 < nb_layers_x(i); i2++)
                                      for (int j2 = 0; j2 < nb_layers_y(j); j2++)
                                        {
                                          int in = nx(i) + i2;
                                          int jn = ny(j) + j2;
                                          ref_domain_layer(in, jn) = old_ref(i, j);
                                        }
                                  }
                            
                            nb_layers_x.Reallocate(pos_x.GetM()-1); nb_layers_x.Fill(1);
                            nb_layers_y.Reallocate(pos_y.GetM()-1); nb_layers_y.Fill(1);
                          }
                        
                        cout << "References du bord y = ymin, x = xmax, y = ymax, x = xmin" << endl;
                        TinyVector<int, 4> ref_boundary;
			ReadParameters(num);
			for (int i = 0; i < 4; i++)
			  ref_boundary(i) = num(i);
			
                        new_mesh.CreateStructuredMesh(pos_x, pos_y, nb_layers_x, nb_layers_y,
                                                      ref_domain_layer, ref_boundary);
                      }
                      break;
                    }
                  
		  // on rajoute le maillage
                  mesh.AppendMesh(new_mesh, true);
                }
		break;
	      case 5 :
		{
		  
		  Vector<PmlRegionParameter<Dimension2> > pml_areas;
		  bool continue_pml = true;
		  while ( continue_pml)
		    {
		      cout << "1- Indiquez les parametres ini de la pml desiree, par exemple YES PML_X 4.0 15" << endl;
		      string line;
		      //	     getline(cin, line ,'\n');		      
		      ReadParameters(line);
		      VectString parameters;
		      StringTokenize(line, parameters, string(" \t"));
		      PmlRegionParameter<Dimension2> new_pml;
		      new_pml.SetParameters(parameters);
		      pml_areas.PushBack(new_pml);
		      
		      cout << "2- Continuer l'ajout de Pml ? 1 pour oui, 0 pour non" << endl;
		      VectReal_wp reponse;
		      ReadParameters(reponse);
		      if(reponse(0)==0)
			continue_pml = false;
		    }
		  
		  for (int i = 0; i < pml_areas.GetM(); i++)
		    pml_areas(i).AddPML(0, mesh);

		  break;
		  cout << "1- Indiquer les bornes xmin, xmax, ymin, ymax" << endl;
		  ReadParameters(coord);
		  Real_wp xmin = coord(0), xmax  = coord(1), ymin  = coord(2), ymax  = coord(3);
		  cout << "2- Indiquer les pas de mailles dx et dy souhaites" << endl;
		  ReadParameters(coord);
		  Real_wp dx = coord(0), dy = coord(1);
		  
		  Real_wp xmin_m = mesh.GetXmin(), xmax_m = mesh.GetXmax();
		  Real_wp ymin_m = mesh.GetYmin(), ymax_m = mesh.GetYmax();
		  Real_wp delta_xmin(0), delta_ymin(0), delta_xmax(0), delta_ymax(0);
		  int nb_div_x0(0), nb_div_xN(0), nb_div_y0(0), nb_div_yN(0);
		  if (abs(xmin-xmin_m) > R2::threshold)
		    {
		      delta_xmin = xmin - xmin_m;
		      nb_div_x0 = toInteger(round(abs(delta_xmin)/dx));
		    }

		  if (abs(xmax-xmax_m) > R2::threshold)
		    {
		      delta_xmax = xmax - xmax_m;
		      nb_div_xN = toInteger(round(abs(delta_xmax)/dx));
		    }

		  if (abs(ymin-ymin_m) > R2::threshold)
		    {
		      delta_ymin = ymin - ymin_m;
		      nb_div_y0 = toInteger(round(abs(delta_ymin)/dy));
		    }

		  if (abs(xmin-xmin_m) > R2::threshold)
		    {
		      delta_ymax = ymax - ymax_m;
		      nb_div_yN = toInteger(round(abs(delta_ymax)/dy));
		    }
		  
		  if (nb_div_x0 > 0)
		    mesh.ExtrudeCoordinate(0, nb_div_x0, xmin_m, delta_xmin);
		  
		  if (nb_div_xN > 0)
		    mesh.ExtrudeCoordinate(0, nb_div_xN, xmax_m, delta_xmax);
		  
		  if (nb_div_y0 > 0)
		    mesh.ExtrudeCoordinate(1, nb_div_y0, ymin_m, delta_ymin);
		  
		  if (nb_div_yN > 0)
		    mesh.ExtrudeCoordinate(1, nb_div_yN, ymax_m, delta_ymax);
		}
		break;
	      }
	    break;
	  }
	case 2 :
	  {
	    // menu pour supprimer
	    cout<<"1- Sommet"<<endl;
	    cout<<"2- Sommets contenus dans un pave"<<endl;
	    cout<<"3- Arete"<<endl;
	    cout<<"4- Aretes de meme reference"<<endl;
	    cout<<"5- Element"<<endl;
	    cout<<"6- Elements de meme reference"<<endl;
	    cout<<"7- Tout sauf certains elements"<<endl;
	    cout<<"Autre- Retour au menu principal" << endl;
	    ReadParameters(num); int choix_retrait = num(0);
	    VectReal_wp delim;
	    switch (choix_retrait)
	      {
	      case 1 :
		cout<<"Entrez les numeros des sommets"<<endl;
		ReadParameters(num);
		RemoveVertex(mesh, num);
		break;
	      case 2 :
		cout<<"Entrez xmin xmax ymin ymax de la zone a supprimer"<<endl;
		ReadParameters(delim);
		RemoveVertex(mesh, delim);
		break;
	      case 3 :
		cout<<"Entrez les numeros des aretes"<<endl;
		ReadParameters(num);
		RemoveEdge(mesh, num);
		break;
	      case 4 :
		cout<<"Entrez les references des aretes a enlever"<<endl;
		ReadParameters(num);
		RemoveEdgeRef(mesh, num);
		break;
	      case 5 :
		cout<<"Entrez les numeros des elements a enlever"<<endl;
		ReadParameters(num);
		RemoveElement(mesh, num);
		break;
	      case 6 :
		cout<<"Entrez les references des elements a enlever"<<endl;
		ReadParameters(num);
		RemoveElementRef(mesh, num);
		break;
	      case 7 :
		{
		  cout<<"Entrez les numeros des elements a conserver"<<endl;
		  ReadParameters(num);
		  Vector<bool> ElementToRemove(mesh.GetNbElt());
		  ElementToRemove.Fill(true);
		  for (int i = 0; i < num.GetM(); i++)
		    ElementToRemove(num(i)) = false;
		  
		  IVect other_num(mesh.GetNbElt()-num.GetM());
		  int ind = 0;
		  for (int i = 0; i < mesh.GetNbElt(); i++)
		    if (ElementToRemove(i))
		      other_num(ind++) = i+1;
		
		  RemoveElement(mesh, other_num);
		}
		break;
	      }
	    break;
	  }
	case 3 :
	  {
	    // menu pour modifier
	    cout<<"1- Coordonnees d'un sommet"<<endl;
	    cout<<"2- Connectique d'une arete"<<endl;
	    cout<<"3- Reference d'une arete"<<endl;
	    cout<<"4- Connectique d'un element"<<endl;
	    cout<<"5- Reference d'un element"<<endl;
	    cout<<"6- Reference surfacique"<<endl;
	    cout<<"7- Reference volumique"<<endl;
	    cout<<"8- Rotation"<<endl;
	    cout<<"9- Translation"<<endl;
	    cout<<"10- Homothetie"<<endl;
	    cout<<"11- (r,teta) -> (x, y)"<<endl;
            cout<<"12- Ordre d'approximation de la geometrie"<<endl;
            cout<<"13- Forcer la coherence du maillage" << endl; 
	    cout<<"14- Spécifier le type de courbe d'une référence" << endl;
	    cout<<"Autre- Retour au menu principal" << endl;
	    ReadParameters(num); int choix_modif = num(0);
	    VectReal_wp delim; int n, ref, nref; Real_wp x, y, teta;
            Real_wp x0, y0;
	    switch (choix_modif)
	      {
	      case 1 :
		cout<<"Entrez le numero du sommet a modifier"<<endl;
		ReadParameters(num); n = num(0);
		cout<<"Entrez les nouvelles coordoonees de ce sommet :  x y "<<endl;
		ReadParameters(coord); x = coord(0); y = coord(1);
		ModifyVertex(mesh, n, R2(x, y));
		break;
	      case 2 :
		cout<<"Entrez le numero de l'arete a modifier"<<endl;
		ReadParameters(num); n = num(0); cout<<"Arete : "<<mesh.BoundaryRef(n-1)<<endl;
		cout<<"Entrez les numeros des sommets de cette arete"<<endl;
		ReadParameters(num);
		ModifyEdge(mesh, n, num(0), num(1));
		break;
	      case 3 :
		cout<<"Entrez le numero de l'arete a modifier"<<endl;
		ReadParameters(num); n = num(0);
		cout<<"Arete : "<<mesh.BoundaryRef(n-1)<<endl;
		cout<<"Entrez la nouvelle reference "<<endl;
		ReadParameters(num); ref = num(0);
		ModifyEdgeRef(mesh, n, ref);
		break;
	      case 4 :
		cout<<"Entrez le numero de l'element a modifier"<<endl;
		ReadParameters(num); n = num(0);
		cout<<"Element : "<<mesh.Element(n-1)<<endl;
		cout<<"Entrez les numeros des sommets de cet element"<<endl;
		ReadParameters(num);
		ModifyElement(mesh, n, num);
		break;
	      case 5 :
		cout<<"Entrez le numero de l'element a modifier"<<endl;
		ReadParameters(num); n = num(0);
		cout<<"Element : "<<mesh.Element(n-1)<<endl;
		cout<<"Entrez la nouvelle reference "<<endl;
		ReadParameters(num); ref = num(0);
		ModifyElementRef(mesh, n, ref);
		break;
	      case 6 :
		{
		  cout<<"Entrez la reference surfacique a changer"<<endl;
		  ReadParameters(num); ref = num(0);
		  cout<<"Entrez la nouvelle reference"<<endl;
		  ReadParameters(num); nref = num(0);
		  cout << "Voulez-vous changer cette referencer sur : " << endl;
		  cout << "1- Tout le maillage " << endl;
		  cout << "2- Pour tous les points y < 0" << endl;
                  cout << "3- Pour tous les points x < 0" << endl;
		  ReadParameters(num);
                  int coor_mode = -1;
		  if (num(0) == 2)
		    coor_mode = 1;
                  else if (num(0) == 3)
                    coor_mode = 0;
		  
		  ModifySurfaceRef(mesh, ref, nref, coor_mode);
		}
		break;
	      case 7 :
		{
		  cout<<"Entrez la reference volumique a changer"<<endl;
		  ReadParameters(num); ref = num(0);
		  cout<<"Entrez la nouvelle reference"<<endl;
		  ReadParameters(num); nref = num(0);
		  
		  cout << "Voulez-vous changer cette referencer sur : " << endl;
		  cout << "1- Tout le maillage " << endl;
		  cout << "2- Pour tous les points y < 0" << endl;
		  ReadParameters(num);
		  bool modify_all = true;
		  if (num(0) == 2)
		    modify_all = false;
		  
		  ModifyVolumeRef(mesh, ref, nref, modify_all);
		}
		break;
	      case 8 :
		cout<<"Entrez le centre de la rotation"<<endl;
		ReadParameters(coord); x0 = coord(0); y0 = coord(1);
		cout<<"Entrez l'angle de rotation"<<endl;
		ReadParameters(coord); teta = coord(0)*pi_wp/180;
		RotateMesh(mesh, R2(x0, y0), teta);
		break;
	      case 9 :
		cout<<"Entrez le vecteur de translation"<<endl;
		ReadParameters(coord); x = coord(0); y = coord(1);
		TranslateMesh(mesh, R2(x, y));
		break;
	      case 10 :
		cout<<"Entrez les facteurs d'echelle"<<endl;
		ReadParameters(coord); x = coord(0); y = coord(1);
		ScaleMesh(mesh, R2(x, y));
		break;
	      case 11 :
		cout<<"Entrez rmin, rmax, teta_min and teta_max in degrees"<<endl;
		ReadParameters(coord); 
		x0= coord(0); y0 = coord(1); x = coord(2); y = coord(3);
		x *= pi_wp/180;
                y *= pi_wp/180;
		TransformPolarToCartesianMesh(mesh, x0, y0, x, y);
		break;
              case 12 :
                cout<<"Entrez l'ordre d'approximation de la geometrie"<<endl;
                ReadParameters(num); order_geom = num(0);
                mesh.SetGeometryOrder(order_geom);
                if (order_geom > 1)
                  output_name = string("toto.msh");
		else
		  output_name = string("toto.mesh");
                
                break;
              case 13 :
                mesh.ForceCoherenceMesh();
                break;
	      case 14 : 
       		// ici on donne le TypeCurve pour rendre courbe un bord
		bool continue_type = true;
		while (continue_type)
		  {
		    cout<< "1- Indiquez les parametres de type de courbe " << endl;
		    string line;
		    VectString parameters;
		    string param;
		    ReadParameters(line);
		    StringTokenize(line, parameters, string(" \t"));
		    IVect ref; ref.Reallocate(1);
		    ref(0) = to_num<int>(parameters(0));
		    param = parameters(1);
		    mesh.ReadCurveType(ref, param);
		    mesh.FindParametersCurve();
		    cout<< "2- Continuer l'ajout d'éléments courbes ? 1 pour oui, 0 pour non"<< endl;
		    VectReal_wp reponse;
		    ReadParameters(reponse);
		    if(reponse(0)==0)
		      continue_type = false;
		  }
		break;
	      }
	    break;
	  }
	case 4 :
	  {
	    // menu pour decouper
	    cout<<"1- Decouper chaque element en triangles (un quad -> 2 triangles) "<<endl;
	    cout<<"2- Decouper chaque element en petits elements (chaque arete est subdivise en n morceaux) "<<endl;
	    cout<<"3- Decouper chaque triangle en trois quads, et chaque quad en quatre quads "<<endl;
	    cout<<"4- Decouper chaque element en petits elements (avec points de Gauss-Lobatto) "<<endl;
            cout << "5- Raffiner les elements d'une meme reference" << endl;
	    cout<<"Autre- Retour au menu principal" << endl;
	    ReadParameters(num); int choix_decoupe = num(0);
	    switch (choix_decoupe)
	      {
	      case 1 :
		mesh.SplitIntoTriangles();
		break;
              case 2 :
                {
		  cout<<"Entrez le nombre de subdivisions"<<endl;
		  ReadParameters(num); int nb_subdiv = num(0);
		  VectReal_wp step_x(nb_subdiv+1);
		  step_x.Fill();
		  Mlt(1.0/nb_subdiv, step_x);
		  mesh.SubdivideMesh(step_x);
		}
		break;
	      case 3 :
		mesh.SplitIntoQuadrilaterals();
		break;
              case 4 :
                {
		  cout<<"Entrez le nombre de subdivisions"<<endl;
		  ReadParameters(num); int nb_subdiv = num(0);
		  VectReal_wp step_x, poids;
		  ComputeGaussLobatto(step_x, poids, nb_subdiv);
		  mesh.SubdivideMesh(step_x);
		}                
		break;
              case 5 :
                {
                  cout << "Entrez la reference a raffiner" << endl;
                  ReadParameters(num); int ref = num(0);
                  cout << "Entrez le pas de maillage cible" << endl;
                  ReadParameters(coord); Real_wp h = coord(0);
                  cout << "Entrez 1 pour un raffinement anisotrope, 0 sinon" << endl;
                  ReadParameters(num);
                  bool anisotrope = false;
                  if (num(0) == 1)
                    anisotrope = true;
                  
                  RefineLocallyMesh(mesh, ref, h, anisotrope);
                }
                break;
	      }
	    break;
	  }
	case 5 :
	  {
	    cout << "1- Verifier conformite maillage volumique / maillage surfacique " << endl;
	    cout << "2- Afficher les statistiques du maillage " << endl;
            cout << "3- Afficher les details d'une face "<< endl;
	    cout << "5- Distance entre deux sommets" << endl;
	    cout<<"Autre- Retour au menu principal" << endl;
	    ReadParameters(num); int choix_verif = num(0);
	    switch (choix_verif)
	      {
	      case 1 :
		try 
		  {
		    mesh.FindConnectivity();
                    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
                      if (mesh.BoundaryRef(i).GetNbElements() == 0)
                        {
                          cout << "Arete de reference " << i+1 << " n appartient pas au maillage volumique " << endl;
                          cout << endl;
                        }                   
		  }
		catch (const WrongMesh& err)
		  {
		    cout<<"WARNING : Mesh not conform"<<endl;
		  }
		break;
	      case 2 :
		{
                  int ref_max = 0;
                  for (int i = 0; i < mesh.GetNbElt(); i++)
                    ref_max = max(ref_max, mesh.Element(i).GetReference());
                  
                  Vector<int> nb_elt_per_ref(ref_max+1);
                  nb_elt_per_ref.Fill(0);
                  for (int i = 0; i < mesh.GetNbElt(); i++)
                    nb_elt_per_ref(mesh.Element(i).GetReference())++;
                  
                  VectReal_wp length_edge(mesh.GetNbEdges());
                  Real_wp h_min(1e300), h_max(0), h_average(0);
                  for (int i = 0; i < mesh.GetNbEdges(); i++)
                    {
                      int n1 = mesh.GetEdge(i).numVertex(0);
                      int n2 = mesh.GetEdge(i).numVertex(1);
                      length_edge(i) = mesh.Vertex(n1).Distance(mesh.Vertex(n2));
                      h_min = min(h_min, length_edge(i));
                      h_max = max(h_max, length_edge(i));
                      h_average += length_edge(i);
                    }
                  
                  h_average /= mesh.GetNbEdges();
                  Sort(length_edge);
                  length_edge.WriteText("h.dat");
                  
                  cout << "Nombre global d'elements : " << mesh.GetNbElt() << endl;
                  cout << "Nombre global d'aretes : " << mesh.GetNbEdges() << endl;
                  cout << "Pas de maillage (longueur moyenne des aretes) : " << h_average << endl;
                  cout << "h_min, h_max = " << h_min << " " << h_max << endl;
                  
                  for (int ref = 0; ref < nb_elt_per_ref.GetM(); ref++)
                    if (nb_elt_per_ref(ref) > 0)
                      cout << "Nombre elements dans le domaine " << ref << " = " << nb_elt_per_ref(ref) << endl;
                  
		}
		break;
              case 3 :
                {
                  cout << "Entrez le numero de la face de reference " <<endl;
                  ReadParameters(num); int nf = num(0) - 1;
                  cout << mesh.BoundaryRef(nf) << endl;
                }
                break;
              case 5 :
		{
		  cout << "Entrez les numeros des deux sommets" << endl;
		  ReadParameters(num); 
		  int nv0 = num(0) - 1, nv1 = num(1) - 1;
		  cout << "Distance entre les deux : " << mesh.Vertex(nv0).Distance(mesh.Vertex(nv1)) << endl;
		}
		break;

	      }
	  }
	  break;
	case 6 :
	  SaveMesh(mesh, output_name);
	  break;
	case 7 :
	  test_loop = false;
	  break;
	case 8 :
	  mesh.FindConnectivity();
	  break;
	case 9:
	  {
	    cout << "1- Extraire les sommets sur un axe" << endl;
	    cout << "2- Extraire les sommets sur un cercle" << endl;
	    ReadParameters(num); int choix_extract = num(0);
	    switch (choix_extract)
	      {
	      case 1:
		{
		  cout << "Entrez le numero de l'axe (0 ou 1) " << endl;
		  ReadParameters(num); int choix_axis = num(0);
		  
		  int nb_vertices = 0;
		  for (int i = 0; i < mesh.GetNbVertices(); i++)
		    if (abs(mesh.Vertex(i)(1-choix_axis)) < R2::threshold)
		      nb_vertices++;

		  VectR2 points(nb_vertices);
		  nb_vertices = 0;
		  for (int i = 0; i < mesh.GetNbVertices(); i++)
		    if (abs(mesh.Vertex(i)(1-choix_axis)) < R2::threshold)
		      points(nb_vertices++) = mesh.Vertex(i);

		  Sort(points);
		  ofstream file_out("Points.txt"); file_out.precision(15);
		  for (int i = 0; i < nb_vertices; i++)
		    file_out << points(i)(0) << " 0 " << points(i)(1) << '\n';

		  file_out.close();
		}
		break;
	      case 2:
		{
		  cout << "Entrez le rayon du cercle " << endl;
		  ReadParameters(coord);
		  Real_wp radius = coord(0);
		  
		  int nb_vertices = 0;
		  for (int i = 0; i < mesh.GetNbVertices(); i++)
		    if (abs(Norm2(mesh.Vertex(i)) - radius) < R2::threshold)
		      nb_vertices++;
		  
		  VectR2 points(nb_vertices);
		  VectReal_wp teta(nb_vertices);
		  nb_vertices = 0; Real_wp r, theta;
		  for (int i = 0; i < mesh.GetNbVertices(); i++)
		    if (abs(Norm2(mesh.Vertex(i)) - radius) < R2::threshold)
		      {
			CartesianToPolar(mesh.Vertex(i)(0), mesh.Vertex(i)(1), r, theta);
			points(nb_vertices) = mesh.Vertex(i);
			teta(nb_vertices++) = theta;
		      }

		  Vector<int> permut(teta.GetM()); permut.Fill();
		  Sort(teta, permut);
		  ofstream file_out("Points.txt"); file_out.precision(15);
		  for (int i = 0; i < nb_vertices; i++)
		    file_out << points(permut(i))(0) << " 0 " << points(permut(i))(1) << '\n';
		  
		  file_out.close();
		}
		break;
	      }
	  }
	  break;
	}
    }
  
  return FinalizeMontjoie();
}
