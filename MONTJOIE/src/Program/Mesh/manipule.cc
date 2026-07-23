#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

#include "Mesh/MeshGeneration.cxx"

using namespace Montjoie;

void AddVertex(Mesh<Dimension3>& mesh, const R3& point)
{
  int nv = mesh.GetNbVertices();
  mesh.ResizeVertices(nv+1);
  mesh.Vertex(nv) = point;
  cout<<"Numero du nouveau point"<<nv+1<<endl;
}

void AddFace(Mesh<Dimension3>& mesh, IVect num, int ref)
{
  int nb_vert = mesh.GetNbVertices();
  int n = mesh.GetNbBoundaryRef();
  mesh.ResizeBoundariesRef(n+1);
  for (int i = 0; i < num.GetM(); i++)
    {
      num(i)--;
      if ((num(i) < 0)||(num(i) >= nb_vert))
	{
	  cout<<"Vertex number must be between 1 and "<<nb_vert<<endl;
	  return;
	}
    }
  
  mesh.BoundaryRef(n).Init(num, ref);
  cout<<"Numero de la nouvelle face "<<n+1<<endl;
}

void AddElement(Mesh<Dimension3>& mesh, IVect num, int ref)
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

void RemoveVertex(Mesh<Dimension3>& mesh, const Vector<bool>& VertexUsed)
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
  for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
    {
      int n1 = mesh.EdgeRef(i).numVertex(0);
      int n2 = mesh.EdgeRef(i).numVertex(1);
      if ( (VertexUsed(n1)) && (VertexUsed(n2)) )
	{
	  int ref = mesh.EdgeRef(i).GetReference();
	  mesh.EdgeRef(nb_edge).Init(new_number(n1), new_number(n2), ref);
	    
	  if (nb_edge != i)
	    {
	      for (int k = 0; k < r-1; k++)
		mesh.SetPointInsideEdge(nb_edge, k, mesh.GetPointInsideEdge(i, k));
	    }
	  
	  nb_edge++;
	}
    }
  
  mesh.ResizeEdgesRef(nb_edge);
  
  // on modifie les faces de reference
  int nb_face = 0;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      bool keep_face = true;
      nb_vert = mesh.BoundaryRef(i).GetNbVertices();
      int ref = mesh.BoundaryRef(i).GetReference();
      IVect numv(nb_vert);
      for (int j = 0; j < nb_vert; j++)
	{
	  int nv = mesh.BoundaryRef(i).numVertex(j);
	  numv(j) = new_number(nv);
	  if (!VertexUsed(nv))
	    keep_face = false;
	}
      
      if (keep_face)
	{
	  mesh.BoundaryRef(nb_face).Init(numv, ref);      

	  if (i != nb_face)
	    {
	      int nb_pts = mesh.GetNbPointsInsideFace(i);	  
	      VectR3 Pts;
	      if (nb_pts > 0)
		{
		  Pts.Reallocate(nb_pts);
		  for (int k = 0; k < nb_pts; k++)
		    Pts(k) = mesh.GetPointInsideFace(i, k);
		}
	      
	      mesh.SetPointInsideFace(nb_face, Pts);
	    }
	  
	  nb_face++;
	}
    }
  
  if (nb_face != mesh.GetNbBoundaryRef())
    mesh.ResizeBoundariesRef(nb_face);
  
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


void RemoveVertex(Mesh<Dimension3>& mesh, IVect num)
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

void RemoveVertex(Mesh<Dimension3>& mesh, const VectReal_wp& delim)
{
  if (delim.GetM() < 6)
    {
      cout<<"Enter at least six real numbers : xmin, xmax, ymin, ymax, zmin, zmax"<<endl;
      return;
    }
  
  Real_wp xmin = delim(0), xmax = delim(1), ymin = delim(2), ymax = delim(3), zmin = delim(4), zmax = delim(5);
  int nb_vert = mesh.GetNbVertices();
  Vector<bool> VertexUsed(nb_vert);
  for (int i = 0; i < nb_vert; i++)
    {
      if ((mesh.Vertex(i)(0) > xmin)&&(mesh.Vertex(i)(0) < xmax)&&
	  (mesh.Vertex(i)(1) > ymin)&&(mesh.Vertex(i)(1) < ymax)&&
	  (mesh.Vertex(i)(2) > zmin)&&(mesh.Vertex(i)(2) < zmax))
	VertexUsed(i) = false;
      else
	VertexUsed(i) = true;
    }
  
  RemoveVertex(mesh, VertexUsed);
}


void ExtractVertex(Mesh<Dimension3>& mesh, const VectReal_wp& delim)
{
  if (delim.GetM() < 6)
    {
      cout<<"Enter at least six real numbers : xmin, xmax, ymin, ymax, zmin, zmax"<<endl;
      return;
    }
  
  Real_wp xmin = delim(0), xmax = delim(1), ymin = delim(2), ymax = delim(3), zmin = delim(4), zmax = delim(5);
  int nb_vert = mesh.GetNbVertices();
  Vector<bool> VertexUsed(nb_vert);
  for (int i = 0; i < nb_vert; i++)
    {
      if ((mesh.Vertex(i)(0) > xmin)&&(mesh.Vertex(i)(0) < xmax)&&
	  (mesh.Vertex(i)(1) > ymin)&&(mesh.Vertex(i)(1) < ymax)&&
	  (mesh.Vertex(i)(2) > zmin)&&(mesh.Vertex(i)(2) < zmax))
	VertexUsed(i) = true;
      else
	VertexUsed(i) = false;
    }
  
  RemoveVertex(mesh, VertexUsed);
}


void RemoveFace(Mesh<Dimension3>& mesh, const Vector<bool>& FaceUsed)
{
  int nb_face = 0;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      if (FaceUsed(i))
	{
	  if (nb_face != i)
	    {
              mesh.BoundaryRef(nb_face) = mesh.BoundaryRef(i);
	      int nb_pts = mesh.GetNbPointsInsideFace(i);
	      VectR3 Pts(nb_pts);
              for (int k = 0; k < mesh.GetNbPointsInsideFace(i); k++)
		Pts(k) = mesh.GetPointInsideFace(i, k);
	      
	      mesh.SetPointInsideFace(nb_face, Pts); 
            }
          
	  nb_face++;
	}
    }
  
  if (nb_face != mesh.GetNbBoundaryRef())
    mesh.ResizeBoundariesRef(nb_face);
  
}

void RemoveFace(Mesh<Dimension3>& mesh, IVect num)
{
  int nb_face = mesh.GetNbBoundaryRef();
  Vector<bool> FaceUsed(nb_face); FaceUsed.Fill(true);
  for (int i = 0; i < num.GetM(); i++)
    {
      num(i)--;
      if ((num(i) < 0)||(num(i) >= nb_face))
	{
	  cout<<"Face number must be between 1 and "<<nb_face<<endl;
	  return;
	}
      
      FaceUsed(num(i)) = false;
    }
  
  RemoveFace(mesh, FaceUsed);
}

void RemoveFaceRef(Mesh<Dimension3>& mesh, IVect num)
{
  int nb_face = mesh.GetNbBoundaryRef();
  Vector<bool> FaceUsed(nb_face); FaceUsed.Fill(true);
  int ref_max = 0;
  for (int i = 0; i < num.GetM(); i++)
    ref_max = max(ref_max, num(i));
  
  Vector<bool> RemoveRef(ref_max+1); RemoveRef.Fill(false);
  for (int i = 0; i < num.GetM(); i++)
    RemoveRef(num(i)) = true;
  
  for (int i = 0; i < nb_face; i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      if (ref >= 0)
        if (ref <= ref_max)
          if (RemoveRef(ref))
            FaceUsed(i) = false;
    }
  
  RemoveFace(mesh, FaceUsed);
}

void RemoveElement(Mesh<Dimension3>& mesh, const Vector<bool>& ElementUsed)
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

void RemoveElement(Mesh<Dimension3>& mesh, IVect num)
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

void RemoveElementRef(Mesh<Dimension3>& mesh, IVect num)
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

void ModifyVertex(Mesh<Dimension3>& mesh, int n, const R3& pt)
{
  if ((n <= 0)||(n > mesh.GetNbVertices()))
    {
      cout<<"Vertex number must be between 1 and "<<mesh.GetNbVertices()<<endl;
      return;
    }
  
  mesh.Vertex(n-1) = pt;
}

void ModifyFace(Mesh<Dimension3>& mesh, int n, IVect num)
{
  if ((n <= 0)||(n > mesh.GetNbBoundaryRef()))
    {
      cout<<"Face number must be between 1 and "<<mesh.GetNbBoundaryRef()<<endl;
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
    
  mesh.BoundaryRef(n-1).Init(num, mesh.BoundaryRef(n-1).GetReference());
}

void ModifyFaceRef(Mesh<Dimension3>& mesh, int n, int ref)
{
  if ((n <= 0)||(n > mesh.GetNbBoundaryRef()))
    {
      cout<<"Face number must be between 1 and "<<mesh.GetNbBoundaryRef()<<endl;
      return;
    }
    
  mesh.BoundaryRef(n-1).SetReference(ref);
}

void ModifyFaceRef(Mesh<Dimension3>& mesh, Vector<bool>& FaceUsed, int ref)
{
  for (int i = 0; i < FaceUsed.GetM(); i++)
    if (FaceUsed(i))
      mesh.BoundaryRef(i).SetReference(ref);
}

void ModifyElement(Mesh<Dimension3>& mesh, int n, IVect num)
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

void ModifyElementRef(Mesh<Dimension3>& mesh, int n, int ref)
{
  if ((n <= 0)||(n > mesh.GetNbElt()))
    {
      cout<<"Element number must be between 1 and "<<mesh.GetNbElt()<<endl;
      return;
    }
    
  mesh.Element(n-1).SetReference(ref);
}


void ModifyElementRef(Mesh<Dimension3>& mesh, const VectReal_wp& delim, int ref)
{
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      R3 center;
      for (int j = 0; j < mesh.Element(i).GetNbVertices(); j++)
        center += mesh.Vertex(mesh.Element(i).numVertex(j));

      Mlt(Real_wp(1) / mesh.Element(i).GetNbVertices(), center);
      if ((center(0) > delim(0)) && (center(0) < delim(1)) &&
          (center(1) > delim(2)) && (center(1) < delim(3)) &&
          (center(2) > delim(4)) && (center(2) < delim(5)) )        
        mesh.Element(i).SetReference(ref);
    }
}


void ModifySurfaceRef(Mesh<Dimension3>& mesh, int old_ref, int new_ref)
{
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    if (mesh.BoundaryRef(i).GetReference() == old_ref)
      mesh.BoundaryRef(i).SetReference(new_ref);  
}

void ModifyVolumeRef(Mesh<Dimension3>& mesh, int old_ref, int new_ref)
{
  for (int i = 0; i < mesh.GetNbElt(); i++)
    if (mesh.Element(i).GetReference() == old_ref)
      mesh.Element(i).SetReference(new_ref);  
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

void ReadParameters(Vector<int>& param)
{
  string chaine;
  getline(cin, chaine ,'\n');
  Vector<string> list_entier;
  StringTokenize(chaine, list_entier, string(" "));
  
  list<int> num;
  for (int i = 0; i < list_entier.GetM(); i++)
    {
      size_t pos = list_entier(i).find("-");
      if (pos != string::npos)
	{
	  int n1 = to_num<int>(list_entier(i).substr(0, pos));
	  int n2 = to_num<int>(list_entier(i).substr(pos+1));
	  for (int n = n1; n <= n2; n++)
	    num.push_back(n);
	}
      else
	num.push_back(to_num<int>(list_entier(i)));
    }
  
  Copy(num, param);
}

void ReadParameters(string& param)
{
  getline(cin, param ,'\n');
}


void SaveMesh(Mesh<Dimension3>& mesh, const string& output_name)
{
  // we detect unused vertices
  int nb_vert = mesh.GetNbVertices();
  Vector<bool> VertexUsed(nb_vert); VertexUsed.Fill(false);
  for (int i = 0; i < mesh.GetNbElt(); i++)
    for (int j = 0; j < mesh.Element(i).GetNbVertices(); j++)
      VertexUsed(mesh.Element(i).numVertex(j)) = true;
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
      VertexUsed(mesh.BoundaryRef(i).numVertex(j)) = true;
  
  // we remove those vertices
  RemoveVertex(mesh, VertexUsed);
  
  mesh.Write(output_name);  
}


template<class T>
void TransformPolarToCartesianMesh(Mesh<Dimension3>& mesh, const T& rmin, const T& rmax,
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
      Real_wp z = mesh.Vertex(i)(2);
      mesh.Vertex(i).Init(x, y, z);
    }
  
  int r = mesh.GetGeometryOrder();
  R3 ptA, ptB;
  for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
    for (int k = 0; k < r-1; k++)
      {
        ptA = mesh.GetPointInsideEdge(i, k);
        Real_wp lambda = (ptA(1) - ymin)/(ymax-ymin);
        Real_wp teta = (1.0-lambda)*teta_min + lambda*teta_max;
        
        lambda = (ptA(0) - xmin)/(xmax-xmin);
        Real_wp r = (1.0-lambda)*rmin + lambda*rmax;
        
        Real_wp x = r*cos(teta);
        Real_wp y = r*sin(teta);
        Real_wp z = ptA(2);
        ptB.Init(x, y, z);
        mesh.SetPointInsideEdge(i, k, ptB);
      }
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    for (int k = 0; k < mesh.GetNbPointsInsideFace(i); k++)
      {
        ptA = mesh.GetPointInsideFace(i, k);
        Real_wp lambda = (ptA(1) - ymin)/(ymax-ymin);
        Real_wp teta = (1.0-lambda)*teta_min + lambda*teta_max;
        
        lambda = (ptA(0) - xmin)/(xmax-xmin);
        Real_wp r = (1.0-lambda)*rmin + lambda*rmax;
        
        Real_wp x = r*cos(teta);
        Real_wp y = r*sin(teta);
        Real_wp z = ptA(2);
        ptB.Init(x, y, z);
        mesh.SetPointInsideFace(i, k, ptB);
      }
  
}

void PickFaceZone(const Mesh<Dimension3>& mesh, const R3& ptA, const R3& ptB, Vector<bool>& FaceUsed)
{
  FaceUsed.Reallocate(mesh.GetNbBoundaryRef());
  FaceUsed.Fill(false);
  TinyVector<R3, 2> box;
  box(0) = ptA; box(1) = ptB;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      R3 center;
      for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
	center += mesh.Vertex(mesh.BoundaryRef(i).numVertex(j));

      Mlt(Real_wp(1) / mesh.BoundaryRef(i).GetNbVertices(), center);
      if (PointInsideBoundingBox(center, box))
	FaceUsed(i) = true;
    }
}

void CreateSphericalBall(const Real_wp& radius_sphere_inside, const Real_wp& radius_sphere_outside,
			 R3& vec_u, int nbx, int nby, int type_mesh, int ref_inside, int ref_outside,
			 int reference_outside_boundary, TinyVector<int, 6>& ref_boundary, Mesh<Dimension3>& new_mesh)
{
  if (type_mesh == Mesh<Dimension3>::HEXAHEDRAL_MESH)
    {
      new_mesh.ReallocateVertices(7);
      new_mesh.Vertex(0).Init(Real_wp(0), Real_wp(0), Real_wp(0));
      new_mesh.Vertex(1).Init(Real_wp(1), Real_wp(0), Real_wp(0));
      new_mesh.Vertex(2).Init(Real_wp(-1), Real_wp(0), Real_wp(0));
      new_mesh.Vertex(3).Init(Real_wp(0), Real_wp(1), Real_wp(0));
      new_mesh.Vertex(4).Init(Real_wp(0), Real_wp(-1), Real_wp(0));
      new_mesh.Vertex(5).Init(Real_wp(0), Real_wp(0), Real_wp(1));
      new_mesh.Vertex(6).Init(Real_wp(0), Real_wp(0), Real_wp(-1));
      new_mesh.ReallocateElements(8);
      new_mesh.Element(0).InitTetrahedral(0, 1, 3, 5, ref_inside);
      new_mesh.Element(1).InitTetrahedral(0, 4, 1, 5, ref_inside);
      new_mesh.Element(2).InitTetrahedral(0, 4, 6, 1, ref_inside);
      new_mesh.Element(3).InitTetrahedral(0, 6, 3, 1, ref_inside);
      new_mesh.Element(4).InitTetrahedral(0, 3, 2, 5, ref_inside);
      new_mesh.Element(5).InitTetrahedral(0, 2, 4, 5, ref_inside);
      new_mesh.Element(6).InitTetrahedral(0, 6, 4, 2, ref_inside);
      new_mesh.Element(7).InitTetrahedral(0, 3, 6, 2, ref_inside);
      new_mesh.ReallocateBoundariesRef(8);
      new_mesh.BoundaryRef(0).InitTriangular(1, 3, 5, reference_outside_boundary);
      new_mesh.BoundaryRef(1).InitTriangular(1, 6, 3, reference_outside_boundary);
      new_mesh.BoundaryRef(2).InitTriangular(1, 4, 5, reference_outside_boundary);
      new_mesh.BoundaryRef(3).InitTriangular(1, 6, 4, reference_outside_boundary);
      new_mesh.BoundaryRef(4).InitTriangular(2, 3, 5, reference_outside_boundary);
      new_mesh.BoundaryRef(5).InitTriangular(2, 5, 4, reference_outside_boundary);
      new_mesh.BoundaryRef(6).InitTriangular(2, 6, 3, reference_outside_boundary);
      new_mesh.BoundaryRef(7).InitTriangular(2, 4, 6, reference_outside_boundary);
      
      new_mesh.ReorientElements();
      new_mesh.FindConnectivity();
      new_mesh.ProjectPointsOnCurves();
    }
  else
    {
      TinyVector<int, 6> ref_bound;
      TinyVector<int, 3> nb_pts_; nb_pts_.Fill(3);
      ref_bound.Fill(reference_outside_boundary);
      new_mesh.CreateRegularMesh(R3(-1,-1,-1), R3(1,1,1), nb_pts_,
				 ref_inside, ref_bound, type_mesh);
    }
  
  R3 point; Real_wp norme_init, coef_max;
  for (int i = 0; i < new_mesh.GetNbVertices(); i++)
    {
      norme_init = Norm2(new_mesh.Vertex(i));
      if (norme_init > 1e-5)
	{
	  point = new_mesh.Vertex(i);
	  coef_max = max(abs(point(0)), abs(point(1)));
	  coef_max = max(abs(point(2)), coef_max);
	  point *= 1.0/coef_max;
	  new_mesh.Vertex(i) *= 1.0/Norm2(point);
	}
    }
  
  //new_mesh.Write("init.mesh");
  
  Real_wp step_h = new_mesh.GetMeshSize();
  Real_wp radius_mesh(0);
  for (int i = 0; i < new_mesh.GetNbVertices(); i++)
    radius_mesh = max(radius_mesh, Norm2(new_mesh.Vertex(i)));
  
  Real_wp coef = radius_sphere_inside/radius_mesh;
  step_h *= coef;
  int nb_points_by_wavelength = nbx/step_h;
  bool is_ext_sphere = false;
  int ref_extern_sphere = 3;
  if (reference_outside_boundary == 3)
    ref_extern_sphere = 4;
  

  new_mesh.CreateSphericalBall(radius_sphere_inside, radius_sphere_outside,
			       nb_points_by_wavelength, reference_outside_boundary,
			       is_ext_sphere, type_mesh, ref_outside, 
                               nbx, nby, ref_extern_sphere);
  
  new_mesh.SetCurveType(reference_outside_boundary, new_mesh.CURVE_SPHERE);
  new_mesh.ProjectPointsOnCurves();
  
  //new_mesh.Write("test.mesh");
  
  // changing points and reference on the boundary
  Real_wp a = radius_sphere_outside;
  
  // x = xmin
  int nb_vert_on_boundary = 0;                        
  for (int i = 0; i < new_mesh.GetNbVertices(); i++)
    if (abs(a + new_mesh.Vertex(i)(0)) < 1e-7)
      nb_vert_on_boundary++;
  
  nbx = toInteger(round(sqrt(Real_wp(nb_vert_on_boundary))));
  VectReal_wp step_div(nbx);
  Real_wp dx = 2.0*a/Real_wp(nbx-1);
  for (int i = 0; i < nbx; i++)
    step_div(i) = -a + i*dx;
  
  IVect IndexVert(nb_vert_on_boundary);
  VectReal_wp x_bound(nb_vert_on_boundary);
  VectReal_wp y_bound(nb_vert_on_boundary);
  nb_vert_on_boundary = 0;
  for (int i = 0; i < new_mesh.GetNbVertices(); i++)
    if (abs(a + new_mesh.Vertex(i)(0)) < 1e-7)
      {
	x_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(1);
	y_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(2);
	IndexVert(nb_vert_on_boundary) = i;
	nb_vert_on_boundary++;
      }
  
  IVect IndexVert2 = IndexVert;
  
  Sort(x_bound, IndexVert);
  Sort(y_bound, IndexVert2);
  
  int offset = 0;
  for (int i = 0; i < nbx; i++)
    {
      Real_wp x = step_div(i);
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert(offset+j);
	  new_mesh.Vertex(nv)(1) = x;
	}
      
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert2(offset+j);
	  new_mesh.Vertex(nv)(2) = x;
	}
      
      offset += nbx;
    }
  
  for (int i = 0; i < new_mesh.GetNbBoundaryRef(); i++)
    {
      R3 center;
      int nb_vert = new_mesh.BoundaryRef(i).GetNbVertices();
      for (int j = 0; j < nb_vert; j++)
	center += new_mesh.Vertex(new_mesh.BoundaryRef(i).numVertex(j));
      
      Mlt(1.0/nb_vert, center);
      if (abs(a + center(0)) < 1e-7)
	new_mesh.BoundaryRef(i).SetReference(ref_boundary(0));
    }                        

  //new_mesh.Write("test1.mesh");  
  // y = ymin
  nb_vert_on_boundary = 0;
  for (int i = 0; i < new_mesh.GetNbVertices(); i++)
    if (abs(a + new_mesh.Vertex(i)(1)) < 1e-7)
      {
	x_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(0);
	y_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(2);
	IndexVert(nb_vert_on_boundary) = i;
	nb_vert_on_boundary++;
      }
  
  IndexVert2 = IndexVert;
  
  Sort(x_bound, IndexVert);
  Sort(y_bound, IndexVert2);
  
  offset = 0;
  for (int i = 0; i < nbx; i++)
    {
      Real_wp x = step_div(i);
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert(offset+j);
	  new_mesh.Vertex(nv)(0) = x;
	}
      
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert2(offset+j);
	  new_mesh.Vertex(nv)(2) = x;
	}
      
      offset += nbx;
    }
  
  for (int i = 0; i < new_mesh.GetNbBoundaryRef(); i++)
    {
      R3 center;
      int nb_vert = new_mesh.BoundaryRef(i).GetNbVertices();
      for (int j = 0; j < nb_vert; j++)
	center += new_mesh.Vertex(new_mesh.BoundaryRef(i).numVertex(j));
      
      Mlt(1.0/nb_vert, center);
      if (abs(a + center(1)) < 1e-7)
	new_mesh.BoundaryRef(i).SetReference(ref_boundary(1));
    }                        
  
  // z = zmin
  nb_vert_on_boundary = 0;
  for (int i = 0; i < new_mesh.GetNbVertices(); i++)
    if (abs(a + new_mesh.Vertex(i)(2)) < 1e-7)
      {
	x_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(0);
	y_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(1);
	IndexVert(nb_vert_on_boundary) = i;
	nb_vert_on_boundary++;
      }
  
  IndexVert2 = IndexVert;
  
  Sort(x_bound, IndexVert);
  Sort(y_bound, IndexVert2);
  
  offset = 0;
  for (int i = 0; i < nbx; i++)
    {
      Real_wp x = step_div(i);
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert(offset+j);
	  new_mesh.Vertex(nv)(0) = x;
	}
      
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert2(offset+j);
	  new_mesh.Vertex(nv)(1) = x;
	}
      
      offset += nbx;
    }
  
  for (int i = 0; i < new_mesh.GetNbBoundaryRef(); i++)
    {
      R3 center;
      int nb_vert = new_mesh.BoundaryRef(i).GetNbVertices();
      for (int j = 0; j < nb_vert; j++)
	center += new_mesh.Vertex(new_mesh.BoundaryRef(i).numVertex(j));
      
      Mlt(1.0/nb_vert, center);
      if (abs(a + center(2)) < 1e-7)
	new_mesh.BoundaryRef(i).SetReference(ref_boundary(2));
    }                        
  
  // x = xmax
  nb_vert_on_boundary = 0;
  for (int i = 0; i < new_mesh.GetNbVertices(); i++)
    if (abs(a - new_mesh.Vertex(i)(0)) < 1e-7)
      {
	x_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(1);
	y_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(2);
	IndexVert(nb_vert_on_boundary) = i;
	nb_vert_on_boundary++;
      }
  
  IndexVert2 = IndexVert;
  
  Sort(x_bound, IndexVert);
  Sort(y_bound, IndexVert2);
  
  offset = 0;
  for (int i = 0; i < nbx; i++)
    {
      Real_wp x = step_div(i);
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert(offset+j);
	  new_mesh.Vertex(nv)(1) = x;
	}
      
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert2(offset+j);
	  new_mesh.Vertex(nv)(2) = x;
	}
      
      offset += nbx;
    }
  
  for (int i = 0; i < new_mesh.GetNbBoundaryRef(); i++)
    {
      R3 center;
      int nb_vert = new_mesh.BoundaryRef(i).GetNbVertices();
      for (int j = 0; j < nb_vert; j++)
	center += new_mesh.Vertex(new_mesh.BoundaryRef(i).numVertex(j));
      
      Mlt(1.0/nb_vert, center);
      if (abs(a - center(0)) < 1e-7)
	new_mesh.BoundaryRef(i).SetReference(ref_boundary(5));
    }                        
  
  // y = ymax
  nb_vert_on_boundary = 0;
  for (int i = 0; i < new_mesh.GetNbVertices(); i++)
    if (abs(a - new_mesh.Vertex(i)(1)) < 1e-7)
      {
	x_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(0);
	y_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(2);
	IndexVert(nb_vert_on_boundary) = i;
	nb_vert_on_boundary++;
      }
  
  IndexVert2 = IndexVert;
  
  Sort(x_bound, IndexVert);
  Sort(y_bound, IndexVert2);
  
  offset = 0;
  for (int i = 0; i < nbx; i++)
    {
      Real_wp x = step_div(i);
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert(offset+j);
	  new_mesh.Vertex(nv)(0) = x;
	}
      
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert2(offset+j);
	  new_mesh.Vertex(nv)(2) = x;
	}
      
      offset += nbx;
    }
  
  for (int i = 0; i < new_mesh.GetNbBoundaryRef(); i++)
    {
      R3 center;
      int nb_vert = new_mesh.BoundaryRef(i).GetNbVertices();
      for (int j = 0; j < nb_vert; j++)
	center += new_mesh.Vertex(new_mesh.BoundaryRef(i).numVertex(j));
      
      Mlt(1.0/nb_vert, center);
      if (abs(a - center(1)) < 1e-7)
	new_mesh.BoundaryRef(i).SetReference(ref_boundary(4));
    }                        
  
  // z = zmin
  nb_vert_on_boundary = 0;
  for (int i = 0; i < new_mesh.GetNbVertices(); i++)
    if (abs(a - new_mesh.Vertex(i)(2)) < 1e-7)
      {
	x_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(0);
	y_bound(nb_vert_on_boundary) = new_mesh.Vertex(i)(1);
	IndexVert(nb_vert_on_boundary) = i;
	nb_vert_on_boundary++;
      }
  
  IndexVert2 = IndexVert;
  
  Sort(x_bound, IndexVert);
  Sort(y_bound, IndexVert2);
  
  offset = 0;
  for (int i = 0; i < nbx; i++)
    {
      Real_wp x = step_div(i);
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert(offset+j);
	  new_mesh.Vertex(nv)(0) = x;
	}
      
      for (int j = 0; j < nbx; j++)
	{
	  int nv = IndexVert2(offset+j);
	  new_mesh.Vertex(nv)(1) = x;
	}
      
      offset += nbx;
    }
  
  for (int i = 0; i < new_mesh.GetNbBoundaryRef(); i++)
    {
      R3 center;
      int nb_vert = new_mesh.BoundaryRef(i).GetNbVertices();
      for (int j = 0; j < nb_vert; j++)
	center += new_mesh.Vertex(new_mesh.BoundaryRef(i).numVertex(j));
      
      Mlt(1.0/nb_vert, center);
      if (abs(a - center(2)) < 1e-7)
	new_mesh.BoundaryRef(i).SetReference(ref_boundary(3));
    }                        
  
  TranslateMesh(new_mesh, vec_u);
}

class FctScaling : public ScalingMeshVirtualFunction3D
{
public :
  Real_wp z0, z1, max_ratio;
  
  void EvaluateCoefficient(const Real_wp& z,
                           Real_wp& coef_x, Real_wp& coef_y)
  {
    Real_wp lambda = 2.0*(z - 0.5*(z0+z1))/(z1-z0);
    coef_x = Real_wp(1);
    coef_y = Real_wp(1);
    if (abs(lambda) <= 1)
      {
        coef_x = 1.0 + (max_ratio-1.0)*square(cos(0.5*pi_wp*lambda));
        coef_y = 1.0/coef_x;
      }
  }
};

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  Mesh<Dimension3> mesh;
  int order_geom = 1;
  string output_name("toto.mesh");
  if (argc > 1)
    {      
      string nom_init(argv[1]);
      int dimension = GetDimensionMesh(nom_init, GetExtension(nom_init));
      if (dimension == 2)
        {
          cout << " Only three-dimensional mesh are handled " << endl;
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

      Real_wp L = mesh.GetXmax() - mesh.GetXmin();
      L = max(L, mesh.GetYmax() - mesh.GetYmin());
      L = max(L, mesh.GetZmax() - mesh.GetZmin());
      R3::threshold *= L;
      DISP(R3::threshold);
      
      cout << "Taille prise par le maillage " << GetHumanReadableMemory(mesh.GetMemorySize()) << endl;
    }
  
  // boucle de manipulation
  bool test_loop = true;
  IVect num, num2; VectReal_wp coord;
  while (test_loop)
    {
      cout<<"1 - AJOUTER sommet, face, element, maillage"<<endl;
      cout<<"2 - ENLEVER sommet, face, element"<<endl;
      cout<<"3 - MODIFIER sommet, face, element, references "<<endl;
      cout<<"4 - DECOUPER"<<endl;      
      cout<<"5 - VERIFIER / INFORMER"<<endl;      
      cout<<"6 - SAUVER"<<endl;
      cout<<"7 - FINIR"<<endl;
      cout<<"8 - METTRE A JOUR CONNECTIVITE"<<endl;
      cout<<"9 - EXTRAIRE maillage" << endl;
      
      ReadParameters(num); int choix = num(0);
      
      switch (choix)
	{
	case 1 :
	  {
	    // menu pour ajouter
	    cout<<"1- Sommet"<<endl;
	    cout<<"2- Face"<<endl;
	    cout<<"3- Element"<<endl;
	    cout<<"4- Maillage"<<endl;
	    cout<<"5- Mailler l'interieur d'une surface fermee (contenant des trous eventuels) " << endl;
	    cout<<"Autre- Retour au menu principal" << endl;
	    ReadParameters(num); int choix_ajout = num(0);
	    Real_wp x, y, z; int ref;
	    string nom_maillage; 
	    switch (choix_ajout)
	      {
	      case 1 :
		cout<<"Entrez les coordonnees x, y, z du sommet"<<endl;
		ReadParameters(coord);
		x = coord(0); y = coord(1); z = coord(2);
		AddVertex(mesh, R3(x, y, z));
		break;
	      case 2 :
		cout<<"Entrez les numeros des sommets de la face"<<endl;
		ReadParameters(num);
		cout<<"Entrez la reference de cette face "<<endl;
		ReadParameters(num2); ref = num2(0);
		AddFace(mesh, num, ref);
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
                  cout << "1- Rajouter un parallelepipede regulier " << endl;
                  cout << "2- Rajouter un maillage regulier d'un cube avec un trou spherique " << endl;
                  cout << "3- Rajouter un maillage regulier d'un cube avec un dielectrique spherique " << endl;
                  cout << "4- Rajouter un maillage d'un fichier " << endl;
                  cout << "5- Rajouter un parallelepipede regulier avec progression geometrique" << endl;
		  cout << "6- Rajouter un maillage extrude a partir d'un maillage 2-D " << endl; 
		  cout << "7- Rajouter une surface quadrilaterale " << endl;
                  cout << "8- Rajouter un parallelepipede regulier avec points personalises" << endl;
		  cout<<"Autre- Retour au menu principal" << endl;
                  ReadParameters(num); int choix_maillage = num(0);
                  Mesh<Dimension3> new_mesh;
                  new_mesh.SetGeometryOrder(order_geom);
                  switch (choix_maillage)
                    {
                    case 1 :
                      {
                        cout << "Donnez xmin, xmax, ymin, ymax, zmin, zmax " << endl;
                        Real_wp xmin, xmax, ymin, ymax, zmin, zmax;
			ReadParameters(coord);
			xmin = coord(0); xmax = coord(1); ymin = coord(2);
			ymax = coord(3); zmin = coord(4); zmax = coord(5);
			cout << "Donnez le nombre de points suivant x y et z " << endl;
                        ReadParameters(num);
			int nbx = num(0), nby = num(1), nbz = num(2);                        
                        cout << "Reference du domaine " << endl;
                        ReadParameters(num); int ref_domain = num(0);
                        cout << "References du bord x = xmin, y = ymin, z = zmin, z = zmax, y = ymax, x = xmax " << endl;
                        TinyVector<int, 6> ref_boundary;
			ReadParameters(num);
			for (int i = 0; i < 6; i++)
			  ref_boundary(i) = num(i);
			
                        cout << "Type de maillage " << endl;
                        cout << "0- Hybride " << endl;
                        cout << "1- Tetrahedres " << endl;
                        cout << "2- Hexahedres " << endl;
                        cout << "3- Pyramides " << endl;
                        cout << "4- Prismes " << endl;
                        ReadParameters(num); int type_mesh = num(0);
                        new_mesh.CreateRegularMesh(R3(xmin, ymin, zmin), R3(xmax, ymax, zmax), TinyVector<int, 3>(nbx, nby, nbz),
                                                   ref_domain, ref_boundary, type_mesh);
                                                
                      }
                      break;
                    case 2 :
                      {
                        Real_wp radius_sphere_inside, radius_sphere_outside, nb_points_by_wavelength;
                        VectReal_wp step_sphere; int nbx; R3 vec_u;
                        cout << "Donnez le rayon de la sphere interieure " << endl;
                        ReadParameters(coord);
			radius_sphere_inside = coord(0);
                        cout << "Donnez le parametre a du cube [-a,a]^3 " << endl;
                        ReadParameters(coord);
			radius_sphere_outside = coord(0);
                        cout << "Donnez le vecteur de translation (centre de la sphere) : ux, uy, uz "<<endl;
                        ReadParameters(coord);
			vec_u(0) = coord(0); vec_u(1) = coord(1); vec_u(2) = coord(2);
                        cout << "Number of points on each side" << endl;
                        ReadParameters(num); nbx = num(0);
                        nb_points_by_wavelength = 6.0*nbx/(sqrt(4.0*pi_wp)*(radius_sphere_inside + radius_sphere_outside));
                        
                        new_mesh.CreateSphericalCrown(radius_sphere_inside, radius_sphere_outside,
                                                      nb_points_by_wavelength, true, false, step_sphere);
                        
                        new_mesh.SetCurveType(1, new_mesh.CURVE_SPHERE);
                        new_mesh.ProjectPointsOnCurves();
                        
                        TranslateMesh(new_mesh, vec_u);
                      }
                      break;
                    case 3 :
                      {                        
                        Real_wp radius_sphere_inside;
                        Real_wp radius_sphere_outside;
                        int reference_outside_boundary, nbx, nby; R3 vec_u;
                        
                        cout << "Donnez le rayon de la sphere interieure " << endl;
                        ReadParameters(coord); radius_sphere_inside = coord(0);
                        cout << "Donnez le parametre a du cube [-a,a]^3 " << endl;
                        ReadParameters(coord); radius_sphere_outside = coord(0);
                        cout << "Donnez le vecteur de translation (centre de la sphere) : ux, uy, uz "<<endl;
                        ReadParameters(coord);
			vec_u(0) = coord(0); vec_u(1) = coord(1); vec_u(2) = coord(2);                        
                        cout << "Factor refinements (surface and interior line)" << endl;
                        ReadParameters(num); nbx = num(0); nby = num(1);
                        
                        cout << "Type de maillage " << endl;
                        cout << "0- Hybride " << endl;
                        cout << "1- Tetrahedres " << endl;
                        cout << "2- Hexahedres " << endl;
                        cout << "3- Pyramides " << endl;
                        cout << "4- Prismes " << endl;
                        ReadParameters(num); int type_mesh = num(0);
                        int ref_outside, ref_inside;
                        cout << "Reference for the inside sphere and outside sphere " << endl;
                        ReadParameters(num); ref_inside = num(0); ref_outside = num(1);
                        cout << "Reference for internal surface (sphere) " << endl;
                        ReadParameters(num); reference_outside_boundary = num(0);
                        

                        cout << "References du bord x = xmin, y = ymin, z = zmin, z = zmax, y = ymax, x = xmax " << endl;
                        TinyVector<int, 6> ref_boundary;
			ReadParameters(num);
			for (int i = 0; i < 6; i++)
			  ref_boundary(i) = num(i);
                        
			CreateSphericalBall(radius_sphere_inside, radius_sphere_outside, vec_u, nbx, nby, type_mesh,
					    ref_inside, ref_outside, reference_outside_boundary, ref_boundary, new_mesh);
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
                        cout << "Donnez xmin, xmax, ymin, ymax, zmin, zmax " << endl;
			ReadParameters(coord);
                        Real_wp xmin = coord(0), xmax  = coord(1), ymin  = coord(2), ymax  = coord(3), zmin  = coord(4), zmax  = coord(5);
			cout << "Donnez le nombre de points suivant x y et z " << endl;
                        ReadParameters(num); int nbx = num(0), nby = num(1), nbz = num(2);
                        cout << "Donnez le ratio pour x y z " << endl;
                        R3 ratio;
                        ReadParameters(coord);
			ratio(0) = coord(0); ratio(1) = coord(1); ratio(2) = coord(2);          
                        cout << "Reference du domaine " << endl;
                        ReadParameters(num); int ref_domain = num(0);
                        cout << "References du bord x = xmin, y = ymin, z = zmin, z = zmax, y = ymax, x = xmax " << endl;
                        TinyVector<int, 6> ref_boundary;
                        ReadParameters(num);
			for (int i = 0; i < 6; i++)
			  ref_boundary(i) = num(i);
                                                
                        cout << "Type de maillage " << endl;
                        cout << "0- Hybride " << endl;
                        cout << "1- Tetrahedres " << endl;
                        cout << "2- Hexahedres " << endl;
                        cout << "3- Pyramides " << endl;
                        cout << "4- Prismes " << endl;
                        ReadParameters(num); int type_mesh = num(0);
                        new_mesh.CreateRegularMesh(R3(xmin, ymin, zmin), R3(xmax, ymax, zmax), TinyVector<int, 3>(nbx, nby, nbz),
                                                   ref_domain, ref_boundary, type_mesh, ratio);
                        
                      }
                      break;
		    case 6 :
		      {
			VectReal_wp step_layer;
			IVect nbCells_layer;
			string file_name;
			
			cout<<"Entrez les differentes positions z = z_i"<<endl;
			ReadParameters(step_layer);
			
			cout << "Entrez le nom du fichier de maillage 2-D " << endl;
			ReadParameters(file_name);
			Mesh<Dimension2> mesh_2d;
			mesh_2d.Read(file_name);
			
			cout << "Entrez le nombre de subdivisions pour chaque couche " << endl;
			ReadParameters(nbCells_layer);
			
			int ref_max_line = 1, ref_max_surf = 1;
			for (int i = 0; i < mesh_2d.GetNbBoundaryRef(); i++)
			  ref_max_line = max(ref_max_line, mesh_2d.BoundaryRef(i).GetReference());
			
			for (int i = 0; i < mesh_2d.GetNbElt(); i++)
			  ref_max_surf = max(ref_max_surf, mesh_2d.Element(i).GetReference());
			
			cout << "Numero maximum des references lineiques = " << ref_max_line << endl;
			cout << "Numero maximum des references surfaciques = " << ref_max_surf << endl;

			int nb_layers = step_layer.GetM() - 1;
			Vector<IVect> ref_surface(2*nb_layers);
			Vector<IVect> ref_volume(2*nb_layers);
			for (int i = 0; i < nb_layers; i++)
			  {
			    ref_surface(2*i).Reallocate(ref_max_line+1);
			    ref_surface(2*i).Fill(0);
			    cout << "Pour la couche numero " << i+1 << ", donnez les nouvelles references des surfaces laterales (" << ref_max_line << " refs a donner) " << endl;
			    ReadParameters(num);
			    for (int j = 0; j < min(int(num.GetM()), ref_max_line); j++)
			      ref_surface(2*i)(j+1) = num(j);

			    ref_surface(2*i+1).Reallocate(ref_max_surf+1);
			    ref_surface(2*i+1).Fill(0);
			    cout << "Pour la couche numero " << i+1 << ", donnez les nouvelles references des surfaces de la section d'au dessus (" << ref_max_surf << " refs a donner) " << endl;
			    ReadParameters(num);
			    for (int j = 0; j < min(int(num.GetM()), ref_max_surf); j++)
			      ref_surface(2*i+1)(j+1) = num(j);
			    
			    ref_volume(i).Reallocate(ref_max_surf+1);
			    ref_volume(i).Fill(0);
			    cout << "Pour la couche numero " << i+1 << ", donnez les nouvelles references des volumes extrudes (" << ref_max_surf << " refs a donner) " << endl;
			    ReadParameters(num);
			    for (int j = 0; j < min(int(num.GetM()), ref_max_surf); j++)
			      ref_volume(i)(j+1) = num(j);
			    
			  }
			
			// ez -> axis of extrusion
			R3 vec_e1, vec_e2, vec_e3;
			vec_e1(0) = 1.0; vec_e2(1) = 1.0; vec_e3(2) = 1.0;
			new_mesh.ExtrudeSurfaceMesh(mesh_2d, nbCells_layer, step_layer,
						    ref_surface, ref_volume, vec_e1, vec_e2, vec_e3);
			
		      }
		      break;
		    case 7 :
		      {
			R3 ptA, ptB, ptC, ptD;
			cout << "Entrez le premier point du quadrilatere " << endl;
			ReadParameters(coord);
			ptA.Init(coord(0), coord(1), coord(2));
			cout << "Entrez le second point du quadrilatere " << endl;
			ReadParameters(coord);
			ptB.Init(coord(0), coord(1), coord(2));
			cout << "Entrez le troisieme point du quadrilatere " << endl;
			ReadParameters(coord);
			ptC.Init(coord(0), coord(1), coord(2));
			cout << "Entrez le quatrieme point du quadrilatere " << endl;
			ReadParameters(coord);
			ptD.Init(coord(0), coord(1), coord(2));	
			cout << "Donnez la reference surfacique du maillage a creer " << endl;
			ReadParameters(num);
			int ref = num(0);
			cout << "Tapez 0 pour un maillage regulier, 1 pour un maillage triangulaire" << endl;
			ReadParameters(num);
			int choix = num(0);
			if (choix == 0)
			  {
			    cout << "Donnez le nombre de points pour le premier segment et le second segment " << endl;
			    ReadParameters(num);
			    CreateSquareMesh(num(0), num(1), ptA, ptB, ptC, ptD, new_mesh);
			    for (int i = 0; i < new_mesh.GetNbBoundaryRef(); i++)
			      new_mesh.BoundaryRef(i).SetReference(ref);
			  }
			else
			  {
			    cout << "Donnez le nombre de points pour les quatre segments " << endl;
			    ReadParameters(num);
			    CreateTriangularMesh(num(0), num(1), num(2), num(3), ptA, ptB, ptC, ptD, new_mesh);
			    for (int i = 0; i < new_mesh.GetNbBoundaryRef(); i++)
			      new_mesh.BoundaryRef(i).SetReference(ref);
			  }
		      }
		      break;
                    case 8 :
                      {
                        cout << "Donnez le nombre de couches suivant x y et z " << endl;
                        ReadParameters(num);
			int nbx = num(0), nby = num(1), nbz = num(2);                       
                        cout << "Donnez les positions x " << endl;
                        VectReal_wp pos_x;
			ReadParameters(pos_x);
                        cout << "Donnez les positions y " << endl;
                        VectReal_wp pos_y;
			ReadParameters(pos_y);
                        cout << "Donnez les positions z " << endl;
                        VectReal_wp pos_z;
			ReadParameters(pos_z);
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
                        
                        cout << "Donnez le nombre de sous-couches en z " << endl;
                        Vector<int> nb_layers_z; ReadParameters(nb_layers_z);
                        if ((nb_layers_z.GetM() != pos_z.GetM()-1) || (nb_layers_z.GetM() != nbz))
                          {
                            cout << "Nombre couches (N) = " << nbz << endl;
                            cout << "Nombre d'entiers = (egal a N) " << nb_layers_z.GetM() << endl;
                            cout << "Nombre de positions = (egal a N+1) " << pos_z.GetM() << endl;
                            cout << "Probleme de coherence" << endl;
                            return FinalizeMontjoie();
                          }
                        
                        cout << "Reference du domaine " << endl;
                        ReadParameters(num); int ref_domain = num(0);
                        cout << "References sous-domaine (liste de i, j, k, ref) " << endl;
                        Array3D<int> ref_domain_layer(nbx, nby, nbz);
                        ref_domain_layer.Fill(ref_domain);
                        ReadParameters(num);
                        for (int i = 0; i < num.GetM(); i += 4)
                          ref_domain_layer(num(i), num(i+1), num(i+2)) = num(i+3);

                        // sub-divisions
                        cout << "Souhaitez-vous donner les subdivisions ? 0- Oui, 1-Non" << endl;
                        ReadParameters(num);
                        if (num(0) == 0)
                          {
                            cout << "Donnez les subdivisions pour x" << endl;
                            ReadParameters(pos_x);

                            cout << "Donnez les subdivisions pour y" << endl;
                            ReadParameters(pos_y);

                            cout << "Donnez les subdivisions pour z" << endl;
                            ReadParameters(pos_z);

                            Array3D<int> old_ref(ref_domain_layer);
                            IVect nx(nb_layers_x.GetM()+1), ny(nb_layers_y.GetM()+1), nz(nb_layers_z.GetM()+1);
                            nx(0) = 0; ny(0) = 0; nz(0) = 0;
                            // on calcule les indices cumules dans nx/ny/nz
                            for (int i = 0; i < nb_layers_x.GetM(); i++)
                              nx(i+1) = nx(i) + nb_layers_x(i);

                            for (int i = 0; i < nb_layers_y.GetM(); i++)
                              ny(i+1) = ny(i) + nb_layers_y(i);
                            
                            for (int i = 0; i < nb_layers_z.GetM(); i++)
                              nz(i+1) = nz(i) + nb_layers_z(i);

                            if ((pos_x.GetM() != nx(nb_layers_x.GetM())+1)
                                || (pos_y.GetM() != ny(nb_layers_y.GetM())+1)
                                || (pos_z.GetM() != nz(nb_layers_z.GetM())+1))
                              {
                                cout << "Probleme de coherence" << endl;
                                DISP(nx); DISP(ny); DISP(nz); DISP(pos_x.GetM());
                                DISP(pos_y.GetM()); DISP(pos_z.GetM());
                                return FinalizeMontjoie();
                              }
                            
                            // on recalcule ref_domain_layer                            
                            ref_domain_layer.Reallocate(pos_x.GetM()-1, pos_y.GetM()-1, pos_z.GetM()-1);
                            ref_domain_layer.Fill(ref_domain);
                            for (int i = 0; i < nbx; i++)
                              for (int j = 0; j < nby; j++)
                                for (int k = 0; k < nbz; k++)
                                  if (old_ref(i, j, k) != ref_domain)
                                    {
                                      for (int i2 = 0; i2 < nb_layers_x(i); i2++)
                                        for (int j2 = 0; j2 < nb_layers_y(j); j2++)
                                          for (int k2 = 0; k2 < nb_layers_z(k); k2++)
                                            {
                                              int in = nx(i) + i2;
                                              int jn = ny(j) + j2;
                                              int kn = nz(k) + k2;
                                              ref_domain_layer(in, jn, kn) = old_ref(i, j, k);
                                            }
                                    }
                            
                            nb_layers_x.Reallocate(pos_x.GetM()-1); nb_layers_x.Fill(1);
                            nb_layers_y.Reallocate(pos_y.GetM()-1); nb_layers_y.Fill(1);
                            nb_layers_z.Reallocate(pos_z.GetM()-1); nb_layers_z.Fill(1);
                          }
                        
                        cout << "References du bord x = xmin, y = ymin, z = zmin, z = zmax, y = ymax, x = xmax " << endl;
                        TinyVector<int, 6> ref_boundary;
			ReadParameters(num);
			for (int i = 0; i < 6; i++)
			  ref_boundary(i) = num(i);
			
                        new_mesh.CreateStructuredMesh(pos_x, pos_y, pos_z, nb_layers_x, nb_layers_y, nb_layers_z,
                                                      ref_domain_layer, ref_boundary);
                                                
                      }
                      break;

                    }
                  mesh.AppendMesh(new_mesh, true);
                }
		break;
	      case 5 :
		{
		  Mesh<Dimension3> mesh_tri = mesh;
		  cout << "Donnez la reference des elements interieurs " << endl;
		  ReadParameters(num);
		  int ref = num(0);
		  CreateTetrahedralMesh(mesh_tri, mesh);		  
		  for (int i = 0; i < mesh.GetNbElt(); i++)
		    mesh.Element(i).SetReference(ref);
		  
		  cout << "On a rempli l'interieur du maillage surfacique avec des tetras + pyramides " << endl;
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
	    cout<<"3- Face"<<endl;
	    cout<<"4- Faces de meme reference"<<endl;
	    cout<<"5- Element"<<endl;
	    cout<<"6- Elements de meme reference"<<endl;
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
		cout<<"Entrez xmin xmax ymin ymax zmin zmax de la zone a supprimer"<<endl;
		ReadParameters(delim);
		RemoveVertex(mesh, delim);
		break;
	      case 3 :
		cout<<"Entrez les numeros des faces"<<endl;
		ReadParameters(num);
		RemoveFace(mesh, num);
		break;
	      case 4 :
		cout<<"Entrez les references des faces a enlever"<<endl;
		ReadParameters(num);
		RemoveFaceRef(mesh, num);
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
	      }
	    break;
	  }
	case 3 :
	  {
	    // menu pour modifier
	    cout<<"1- Coordonnees d'un sommet"<<endl;
	    cout<<"2- Connectique d'une face"<<endl;
	    cout<<"3- Reference d'une face"<<endl;
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
            cout<<"14- Projeter les points d'une surface sur un cone/sphere/cylindre" << endl;
            cout<<"15- Separer une surface en plusieurs parties" << endl;
            cout<<"16- Homothetie variable" << endl;
            cout<<"17- Symmetrie par rapport a un plan" << endl;
	    cout<<"18- reference des faces contenus dans un pave" << endl;
	    cout<<"Autre- Retour au menu principal" << endl;
	    ReadParameters(num); int choix_modif = num(0);
	    VectReal_wp delim; int n, ref, nref; Real_wp x, y, z, teta;
            Real_wp x0, y0, z0;
	    switch (choix_modif)
	      {
	      case 1 :
		cout<<"Entrez le numero du sommet a modifier"<<endl;
		ReadParameters(num); n = num(0);
		cout<<"Entrez les nouvelles coordoonees de ce sommet :  x y z "<<endl;
		ReadParameters(coord); x = coord(0); y = coord(1); z = coord(2);
		ModifyVertex(mesh, n, R3(x, y, z));
		break;
	      case 2 :
		cout<<"Entrez le numero de la face a modifier"<<endl;
		ReadParameters(num); n = num(0); cout<<"Face : "<<mesh.BoundaryRef(n-1)<<endl;
		cout<<"Entrez les numeros des sommets de cette face"<<endl;
		ReadParameters(num);
		ModifyFace(mesh, n, num);
		break;
	      case 3 :
		cout<<"Entrez le numero de la face a modifier"<<endl;
		ReadParameters(num); n = num(0);
                
		cout<<"Face : "<<mesh.BoundaryRef(n-1)<<endl;
		cout<<"Entrez la nouvelle reference "<<endl;
		ReadParameters(num); ref = num(0);
		ModifyFaceRef(mesh, n, ref);
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
                if (n <= 0)
                  {
                    cout<<"Entrez xmin xmax ymin ymax zmin zmax de la zone a changer"<<endl;
                    ReadParameters(delim);
                    
                    cout<<"Entrez la nouvelle reference "<<endl;
                    ReadParameters(num); ref = num(0);
                    
                    ModifyElementRef(mesh, delim, ref);
                  }
                else
                  {
                    cout<<"Element : "<<mesh.Element(n-1)<<endl;
                    cout<<"Entrez la nouvelle reference "<<endl;
                    ReadParameters(num); ref = num(0);
                    ModifyElementRef(mesh, n, ref);
                  }
                
		break;
	      case 6 :
		cout<<"Entrez la reference surfacique a changer"<<endl;
		ReadParameters(num); ref = num(0);
		cout<<"Entrez la nouvelle reference"<<endl;
		ReadParameters(num); nref = num(0);
		ModifySurfaceRef(mesh, ref, nref);
		break;
	      case 7 :
		cout<<"Entrez la reference volumique a changer"<<endl;
		ReadParameters(num); ref = num(0);
		cout<<"Entrez la nouvelle reference"<<endl;
		ReadParameters(num); nref = num(0);
		ModifyVolumeRef(mesh, ref, nref);
		break;
	      case 8 :
		cout<<"Entrez un point de l'axe de rotation"<<endl;
		ReadParameters(coord); x0 = coord(0); y0 = coord(1); z0 = coord(2);
		cout<<"Entrez l'axe de rotation"<<endl;
		ReadParameters(coord); x = coord(0); y = coord(1); z = coord(2);
		cout<<"Entrez l'angle de rotation"<<endl;
		ReadParameters(coord); teta = coord(0)*pi_wp/180;
		RotateMesh(mesh, R3(x0, y0, z0), R3(x, y, z), teta);
		break;
	      case 9 :
		cout<<"Entrez le vecteur de translation"<<endl;
		ReadParameters(coord); x = coord(0); y = coord(1); z = coord(2);
		TranslateMesh(mesh, R3(x, y, z));
		break;
	      case 10 :
		cout<<"Entrez les facteurs d'echelle"<<endl;
		ReadParameters(coord); x = coord(0); y = coord(1); z = coord(2);
		ScaleMesh(mesh, R3(x, y, z));
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
              case 14:
                {
                  cout << "Quelles sont les references surfaciques impliquees ?" << endl;
                  ReadParameters(num);
                  IVect ref_cond(mesh.GetNbReferences()+1);
                  ref_cond.Fill(0);
                  for (int i = 0; i < num.GetM(); i++)
                    if ((num(i) >= 0) && (num(i) <= mesh.GetNbReferences()))
                      ref_cond(num(i)) = 1;
                  
                  Vector<bool> VertexUsed(mesh.GetNbVertices()); VertexUsed.Fill(false);
                  
                  cout << "1- Projeter sur une sphere" << endl;
                  cout << "2- Projeter sur un cylindre" << endl;
                  cout << "3- Projeter sur un cone" << endl;
                  ReadParameters(num); int type_surf = num(0);
                  if (type_surf == 1)
                    {
                      cout << "Entrez le centre de la sphere suivi du rayon : x, y, z, r" << endl;
                      ReadParameters(coord);
                      R3 center(coord(0), coord(1), coord(2)); Real_wp radius = coord(3);
                      
                      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
                        if (ref_cond(mesh.BoundaryRef(i).GetReference()) == 1)
                          for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
                            {
                              int nv = mesh.BoundaryRef(i).numVertex(j);
                              if (!VertexUsed(nv))
                                {
                                  mesh.ProjectToSphere(mesh.Vertex(nv), center, radius);
                                  VertexUsed(nv) = true;
                                }
                            }
                    }
                  else if (type_surf == 2)
                    {
                      cout << "Entrez un point de l'axe, l'axe et le rayon : x, y, z, u, v, w, r" << endl;
                      ReadParameters(coord);
                      R3 ptA(coord(0), coord(1), coord(2));
                      R3 axis(coord(3), coord(4), coord(5));
                      Real_wp radius = coord(6);
                      
                      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
                        if (ref_cond(mesh.BoundaryRef(i).GetReference()) == 1)
                          for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
                            {
                              int nv = mesh.BoundaryRef(i).numVertex(j);
                              if (!VertexUsed(nv))
                                {
                                  mesh.ProjectToCylinder(mesh.Vertex(nv), ptA, axis, radius);
                                  VertexUsed(nv) = true;
                                }
                            }
                    }
                  else
                    {
                      cout << "Entrez le centre du cone-sphere, l'axe et l'angle : x, y, z, u, v, w, theta " << endl;
                      ReadParameters(coord);
                      R3 ptA(coord(0), coord(1), coord(2));
                      R3 axis(coord(3), coord(4), coord(5));
                      teta = coord(6);
                      
                      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
                        if (ref_cond(mesh.BoundaryRef(i).GetReference()) == 1)
                          for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
                            {
                              int nv = mesh.BoundaryRef(i).numVertex(j);
                              if (!VertexUsed(nv))
                                {
                                  mesh.ProjectToConicCylinder(mesh.Vertex(nv), ptA, axis, teta);
                                  VertexUsed(nv) = true;
                                }
                            }
                    }
                }
                break;
              case 15 :
                {
                  cout << "Entrez la reference de la surface a separer" << endl;
                  ReadParameters(num);
                  int ref_surf = num(0);
                  int new_ref = mesh.GetNewReference();
		  DISP(new_ref);
		  cout << "Entrez un seuil pour determiner si le point est sur une primitive" << endl;
                  ReadParameters(coord);
                  Real_wp threshold = coord(0);
                  Vector<bool> VertexUsed(mesh.GetNbVertices()); VertexUsed.Fill(false);
                  R3 pt_surf;
                  
                  cout << "1- Isoler la partie spherique de la surface" << endl;
                  cout << "2- Isoler la partie cylindrique de la surface" << endl;
                  cout << "3- Isoler la partie conique de la surface" << endl;
                  
                  ReadParameters(num);
                  int type_surf = num(0);
                  if (type_surf == 1)
                    {
                      ReadParameters(coord);
                      R3 center(coord(0), coord(1), coord(2)); Real_wp radius = coord(3);
                      
                      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
                        if (mesh.BoundaryRef(i).GetReference() == ref_surf)
                          {
                            bool face_on_surf = true;
                            for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
                              {
                                int nv = mesh.BoundaryRef(i).numVertex(j);
                                pt_surf = mesh.Vertex(nv);
                                mesh.ProjectToSphere(pt_surf, center, radius);
                                if (pt_surf.Distance(mesh.Vertex(nv)) > threshold)
                                  face_on_surf = false;
                              }
                            
                            if (face_on_surf)
                              {
                                mesh.BoundaryRef(i).SetReference(new_ref);
                                for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
                                  {
                                    int nv = mesh.BoundaryRef(i).numVertex(j);
                                    if (!VertexUsed(nv))
                                      {
                                        mesh.ProjectToSphere(mesh.Vertex(nv), center, radius);
                                        VertexUsed(nv) = true;
                                      }
                                  }
                              }
                          }

                    }
                  else if (type_surf == 2)
                    {
                      cout << "Entrez un point de l'axe, l'axe et le rayon : x, y, z, u, v, w, r" << endl;
                      ReadParameters(coord);
                      R3 ptA(coord(0), coord(1), coord(2));
                      R3 axis(coord(3), coord(4), coord(5));
                      Real_wp radius = coord(6);
                      
                      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
                        if (mesh.BoundaryRef(i).GetReference() == ref_surf)
                          {
                            bool face_on_surf = true;
                            for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
                              {
                                int nv = mesh.BoundaryRef(i).numVertex(j);
                                pt_surf = mesh.Vertex(nv);
                                mesh.ProjectToCylinder(pt_surf, ptA, axis, radius);
                                if (pt_surf.Distance(mesh.Vertex(nv)) > threshold)
                                  face_on_surf = false;
                              }
                            
                            if (face_on_surf)
                              {
                                mesh.BoundaryRef(i).SetReference(new_ref);
                                for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
                                  {
                                    int nv = mesh.BoundaryRef(i).numVertex(j);
                                    if (!VertexUsed(nv))
                                      {
                                        mesh.ProjectToCylinder(mesh.Vertex(nv), ptA, axis, radius);
                                        VertexUsed(nv) = true;
                                      }
                                  }
                              }
                          }
                    }
                  else
                    {
                      cout << "Entrez le centre du cone-sphere, l'axe et l'angle : x, y, z, u, v, w, theta " << endl;
                      ReadParameters(coord);
                      R3 ptA(coord(0), coord(1), coord(2));
                      R3 axis(coord(3), coord(4), coord(5));
                      teta = coord(6);
                      
                      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
                        if (mesh.BoundaryRef(i).GetReference() == ref_surf)
                          {
                            bool face_on_surf = true;
                            for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
                              {
                                int nv = mesh.BoundaryRef(i).numVertex(j);
                                pt_surf = mesh.Vertex(nv);
                                mesh.ProjectToConicCylinder(pt_surf, ptA, axis, teta);
                                if (pt_surf.Distance(mesh.Vertex(nv)) > threshold)
                                  face_on_surf = false;
                              }
                            
                            if (face_on_surf)
                              {
                                mesh.BoundaryRef(i).SetReference(new_ref);
                                for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
                                  {
                                    int nv = mesh.BoundaryRef(i).numVertex(j);
                                    if (!VertexUsed(nv))
                                      {
                                        mesh.ProjectToConicCylinder(mesh.Vertex(nv), ptA, axis, teta);
                                        VertexUsed(nv) = true;
                                      }
                                  }
                              }
                          }
                    }
                  
                }
                break;
              case 16 :
                {
                  FctScaling fct;
                  cout << "Entrez les extremites z0 et z1, la section variera sur l'intervalle [z0, z1] " << endl;
                  ReadParameters(coord);
                  fct.z0 = coord(0); fct.z1 = coord(1);
                  cout << "Entrez le ratio maximal b/a " << endl;
                  ReadParameters(coord);
                  fct.max_ratio = coord(0);
                  
                  ScaleVariableMesh(mesh, fct);
                }
                break;
              case 17 :
                {
                  cout << "Entrez les coefficients a b, c, d du plan de symetrie" << endl;
                  ReadParameters(coord);
                  mesh.SymmetrizeMesh(-1, coord(3), R3(coord(0), coord(1), coord(2)));
                }
                break;
	      case 18 :
		{
		  cout << "Entrez la zone ou sont les faces dont la reference va changer" << endl;
		  cout << "Entrez xmin xmax ymin ymax zmin zmax" << endl;
		  ReadParameters(coord); Vector<bool> FaceUsed;
		  PickFaceZone(mesh, R3(coord(0), coord(2), coord(4)), R3(coord(1), coord(3), coord(5)), FaceUsed);
		  cout<<"Entrez la nouvelle reference "<<endl;
		  ReadParameters(num); ref = num(0);
		  ModifyFaceRef(mesh, FaceUsed, ref);
		}
		break;
	      }
	    break;
	  }
	case 4 :
	  {
	    // menu pour decouper
	    cout << "1- Decouper chaque element en triangles/tetras (un hexa -> 6 tetras) " << endl;
	    cout << "2- Decouper chaque element en petits elements "
		 << "(chaque arete est subdivise en n morceaux) " << endl;
	    cout << "3- Decouper chaque tetra en quatre hexas " << endl;
	    cout << "4- Raffiner un maillage surfacique avec yams " << endl;
	    cout << "5- Raffiner une partie du maillage " << endl;
	    cout << "Autre- Retour au menu principal" << endl;
	    ReadParameters(num); int choix_decoupe = num(0);
	    switch (choix_decoupe)
	      {
	      case 1 :
		mesh.SplitIntoTetrahedra();
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
		mesh.SplitIntoHexahedra();
		break;
	      case 4 :
		{
                  cout << "Entrez le pas de maillage cible " << endl;
                  ReadParameters(coord);
		  Real_wp hmax = coord(0);

		  // on ecrit le fichier au format mesh
		  mesh.Write("totoYams.mesh");
		  
		  string tmps = "totoYams.bb";
		  ofstream input_mit(tmps.data());
		  input_mit << "3 1 1 " << mesh.GetNbVertices() << " 2" << endl;
		  for (int i = 0; i < mesh.GetNbVertices(); i++)
		    input_mit << hmax << endl;
		  
		  input_mit.close();

		  input_mit.open("totoYams.yams");
		  //input_mit << "MinSize " << endl << hmin << endl;
		  input_mit << "MaxSize " << endl << hmax << endl;
		  input_mit.close();
		  
		  string yams_command = "./yams -m 200";
		  int type_refinement = 1;
		  string command_line = yams_command + " -O " + to_str(type_refinement) + " totoYams outputYams";
		  int success_mesh = system(command_line.data());
		  if (success_mesh != 0)
		    {
		      cout << "Echec lors de l'adaptation de maillage avec YAMS" << endl;
		      abort();
		    }
		  
		  mesh.Read("outputYams.meshb");
		  
		  remove("totoYams.mesh");
		  remove("totoYams.bb");
		  remove("totoYams.yams");
		  remove("outputYams.meshb");
		}
		break;
	      case 5:
		{
		  cout << "Entrez la reference volumique qui doit etre raffinee" << endl;
		  ReadParameters(num);
		  int ref_refine = num(0); int ref_domain = 0;
		  DISP(ref_refine);
		  
		  // on recupere les references surfaciques (bord) des elements en lien avec des references qui ne sont pas raffinees
		  IVect ref_bord(mesh.GetNbReferences()+1);
		  ref_bord.Fill(0); int nb_ref_surf = 0;
		  for (int i = 0; i < mesh.GetNbElt(); i++)
		    {
		      int ref = mesh.Element(i).GetReference();
		      if (ref == ref_refine)
			{
			  for (int j = 0; j < mesh.Element(i).GetNbBoundary(); j++)
			    {
			      int nf = mesh.Element(i).numBoundary(j);
			      if (mesh.Boundary(nf).GetNbElements() == 2)
				{
				  int ne = mesh.Boundary(nf).numElement(0);
				  if (ne == i)
				    ne = mesh.Boundary(nf).numElement(1);
				  
				  int ref2 = mesh.Element(ne).GetReference();
				  if (ref2 != ref_refine)
				    {
				      int ref_surf = mesh.Boundary(nf).GetReference();
				      if (ref_surf == 0)
					{
					  cout << "Il faut creer une Physical Surface pour l'interface" << endl;
					  abort();
					}
				      
				      if (ref_bord(ref_surf) != 1)
					{
					  ref_domain = ref2;
					  ref_bord(ref_surf) = 1;
					  nb_ref_surf++;
					}
				    }
				}
			    }
			}
		    }

		  Vector<int> list_ref_surf(nb_ref_surf);
		  nb_ref_surf = 0;
		  for (int ref = 1; ref < ref_bord.GetM(); ref++)
		    if (ref_bord(ref) == 1)
		      list_ref_surf(nb_ref_surf++) = ref;

		  DISP(list_ref_surf); DISP(ref_domain);
		  SurfacicMesh<Dimension3> mesh_surf;
		  mesh.GetBoundaryMesh(1, mesh_surf, ref_bord, ref_domain);
		  
		  mesh_surf.Write("interface.mesh");

		  IVect Global_to_local(mesh.GetNbVertices());
		  Global_to_local.Fill(-1);
		  for (int i = 0; i < mesh_surf.ListeVertices.GetM(); i++)
		    Global_to_local(mesh_surf.ListeVertices(i)) = i;
		  
		  SurfacicMesh<Dimension3> mesh2(mesh_surf);
		  for (int i = 0; i < mesh_surf.NumElement.GetM(); i++)
		    {
		      int ne = mesh_surf.NumElement(i);
		      int num_loc = mesh_surf.NumLocalBoundary(i);
		      int num_loc2 = 5-num_loc; // pour avoir la surface opposee
		      for (int nv = 0; nv < 4; nv++)
			{
			  int n0 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::FirstExtremityFace(3, num_loc));
			  int n1 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::SecondExtremityFace(3, num_loc));
			  int n2 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::ThirdExtremityFace(3, num_loc));
			  int n3 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::FourthExtremityFace(3, num_loc));

			  int p0 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::FirstExtremityFace(3, num_loc2));
			  int p1 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::SecondExtremityFace(3, num_loc2));
			  int p2 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::ThirdExtremityFace(3, num_loc2));
			  int p3 = mesh.Element(ne).numVertex(MeshNumbering<Dimension3>::FourthExtremityFace(3, num_loc2));

			  int n0_f = Global_to_local(n0);
			  int n1_f = Global_to_local(n1);
			  int n2_f = Global_to_local(n2);
			  int n3_f = Global_to_local(n3);
			  
			  mesh2.Vertex(n0_f) = mesh.Vertex(p0);
			  mesh2.Vertex(n1_f) = mesh.Vertex(p1);
			  mesh2.Vertex(n2_f) = mesh.Vertex(p2);
			  mesh2.Vertex(n3_f) = mesh.Vertex(p3);
			}
		    }

		  mesh2.Write("interface2.mesh");
		  CreateTransitionLayer(mesh2, mesh_surf, ref_domain);
		  mesh2.Write("layer.mesh");
		  
		  VectBool ElementUsed(mesh.GetNbElt());
		  ElementUsed.Fill(false);
		  Mesh<Dimension3> mesh_vol(mesh);
		  for (int i = 0; i < mesh.GetNbElt(); i++)
		    {
		      int ref = mesh.Element(i).GetReference();
		      if (ref == ref_refine)
			ElementUsed(i) = true;
		    }

		  RemoveElement(mesh_vol, ElementUsed);
		  mesh_vol.ForceCoherenceMesh();
		  VectReal_wp step(3); step(0) = 0.0; step(1) = 0.5; step(2) = 1.0;
		  mesh_vol.SubdivideMesh(step);

		  mesh_vol.AppendMesh(mesh2);
		  mesh_vol.Write("refine.mesh");

		  ElementUsed.Fill(false);
		  for (int i = 0; i < mesh.GetNbElt(); i++)
		    {
		      int ref = mesh.Element(i).GetReference();
		      if (ref != ref_refine)
			ElementUsed(i) = true;
		    }
		  
		  for (int i = 0; i < mesh_surf.NumElement.GetM(); i++)
		    {
		      int ne = mesh_surf.NumElement(i);
		      ElementUsed(ne) = false;
		    }
		  
		  RemoveElement(mesh, ElementUsed);
		  mesh.Write("bottom.mesh");
		  mesh.ForceCoherenceMesh();

		  mesh.AppendMesh(mesh_vol);
		}
		break;
	      }
	    break;
	  }
	case 5 :
	  {
	    cout<<"1- Verifier conformite maillage volumique / maillage surfacique "<<endl;
	    cout<<"2- Verifier maillage surfacique (surfaces fermees) "<<endl;
            cout << "3- Afficher les details d'une face "<< endl;
            cout << "4- Trouver un numero de sommet" << endl;
            cout << "5- Distance entre deux sommets" << endl;
            cout << "6- Parametres d'un plan" << endl;
	    cout<<"Autre- Retour au menu principal" << endl;
	    ReadParameters(num); int choix_verif = num(0);
	    bool mesh_ok = true;
	    switch (choix_verif)
	      {
	      case 1 :
		try 
		  {
		    mesh.FindConnectivity();
                    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
                      if (mesh.BoundaryRef(i).GetNbElements() == 0)
                        {
                          cout << "Face de reference " << i+1 << " n appartient pas au maillage volumique " << endl;
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
		  mesh.FindConnectivity();
		  mesh_ok = true;
		  for (int i = 0; i < mesh.GetNbEdges(); i++)
		    if (mesh.GetEdge(i).GetNbFaces() != 2)
		      mesh_ok = false;
		  
		  if (mesh_ok)
		    {
		      cout<<"Surface mesh is correct"<<endl;
		    }
		  else
		    {
		      cout<<"Surface mesh is incorrect "<<endl;
		    }
		}
		break;
              case 3 :
                {
                  cout << "Entrez le numero de la face de reference " <<endl;
                  ReadParameters(num); int nf = num(0) - 1;
                  cout << mesh.BoundaryRef(nf) << endl;
                }
                break;
              case 4 :
                {
                  cout << "Entrez le centre de la boule qui contient le sommet cherche" << endl;
                  ReadParameters(coord);
                  R3 center(coord(0), coord(1), coord(2));
                  cout << "Entrez le rayon de la boule qui contient le sommet cherche" << endl;
                  ReadParameters(coord);
                  Real_wp radius = coord(0);
                  for (int i = 0; i < mesh.GetNbVertices(); i++)
                    if (mesh.Vertex(i).Distance(center) <= radius)
                      cout << "Sommet " << i+1 << " de coordonnees " << mesh.Vertex(i) << endl;
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
	      case 6 :
		{
		  cout << "Attention ca ne marche pas si x = qqchose ou z = qqchose" << endl;
		  cout << "Entrez les numeros des sommets du plan" << endl;
		  ReadParameters(num);
		  int N = num.GetM();
		  for (int i = 0; i<N ; i++)
                    num(i)-- ;
                  
		  Matrix<Real_wp> Vdm(N,3);
		  for (int i = 0; i < N; i++)
		    {
		      Vdm(i,0) = 1.0;
		      Vdm(i,1) = mesh.Vertex(num(i))(0);
		      Vdm(i,2) = mesh.Vertex(num(i))(2);
		    }
                  
		  Matrix<Real_wp> Vdm0 = Vdm;
		  VectReal_wp tau;
		  GetQR(Vdm,tau);
		  VectReal_wp rhs(N), x_sol; 
		  for (int i = 0; i < N; i++)
                    rhs(i) = -mesh.Vertex(num(i))(1);

		  x_sol = rhs;
		  SolveQR(Vdm,tau,x_sol);
                  
		  // rhs contient les parametres du plan
                  cout << "Le plan le plus proche des points donnees est le plan a x + b y + c z + d = 0, avec" << endl;
                  cout << "a = " << x_sol(1) << endl;
                  cout << "b = 1" << endl;
                  cout << "c = " << x_sol(2) << endl;
                  cout << "d = " << x_sol(0) << endl;
                  
		  MltAdd(-1.0, Vdm0, x_sol, 1.0, rhs);
		  cout << "Erreur sur chaque point : " << rhs << endl;
                  cout << "Erreur globale : " << Norm2(rhs) << endl;
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
	case 9 :
	  {
	    cout << "1- Extraire le maillage de peau " << endl;
	    cout << "2- Extraire les elements voisins d'un element central " << endl;
	    cout << "3- Extraire les elements d'une zone" << endl;
	    cout << "Autre- Retour au menu principal " << endl;
	    ReadParameters(num);
	    int choix_extrac = num(0);
	    Mesh<Dimension3> sub_mesh;
	    sub_mesh.SetGeometryOrder(mesh.GetGeometryOrder());
	    
	    switch (choix_extrac)
	      {
	      case 1 :
		{
		  IVect ListeFace, ListeVertex, NumElem, NumLocal;
		  IVect ref_cond(mesh.GetNbReferences()+1);
		  ref_cond.Fill(0);
		  cout << "Entrez les references des faces a extraire " << endl;
		  ReadParameters(num);
		  for (int i = 0; i < num.GetM(); i++)
		    ref_cond(num(i)) = 1;
		  
		  mesh.GetBoundaryMesh(1, sub_mesh, ListeFace, ListeVertex, 
				       NumElem, NumLocal, ref_cond, 0);
		}
		break;
	      case 2 :
		{
		  cout << "Entrez le numero de l'element autour duquel vous voulez extraire les voisins" << endl;
		  ReadParameters(num); int ne = num(0) - 1;
		  
		  if ((ne < 0) || (ne >= mesh.GetNbElt()))
		    {
		      cout << "Echec : Entrez un numero d'element compris entre 1 et " << mesh.GetNbElt() << endl;
		    }
		  else
		    {
		      int nb_vertex = 0;
		      int nb_elt = 0;
		      Vector<bool> VertexOnSubdomain(mesh.GetNbVertices()), ElementOnSubdomain(mesh.GetNbElt());
		      VertexOnSubdomain.Fill(false); ElementOnSubdomain.Fill(false);
		      for (int i = 0; i < mesh.Element(ne).GetNbFaces(); i++)
			{
			  int nf = mesh.Element(ne).numFace(i);
			  for (int k = 0; k < mesh.Boundary(nf).GetNbElements(); k++)
			    {
			      int num_elem = mesh.Boundary(nf).numElement(0);			
			      if (!ElementOnSubdomain(num_elem))
				{
				  ElementOnSubdomain(num_elem) = true;
				  nb_elt++;
				  for (int j = 0; j < mesh.Element(num_elem).GetNbVertices(); j++)
				    {
				      int nv = mesh.Element(num_elem).numVertex(j);
				      if (!VertexOnSubdomain(nv))
					{
					  VertexOnSubdomain(nv) = true;
					  nb_vertex++;
					}
				    }
				}
			    }
			}
		      
		      mesh.CreateSubmesh(sub_mesh, nb_vertex, nb_elt, VertexOnSubdomain, ElementOnSubdomain);
		      sub_mesh.AddBoundaryFaces();
		      sub_mesh.FindConnectivity();
		    }
		}
		break;
	      case 3 :
		{
		  VectReal_wp delim;
		  cout << "Entrez xmin, xmax, ymin, zmin et zmax" << endl;
		  ReadParameters(delim);
		  sub_mesh = mesh;
		  ExtractVertex(sub_mesh, delim);
		}
		break;
	      default:
		{
		}
		break;
	      }
	    
	    sub_mesh.Write("sub.mesh");
	    cout << "Maillage sauvegarde sous le nom sub.mesh" << endl;
	  }
	  break;
	}
    }
  
  return FinalizeMontjoie();
}
