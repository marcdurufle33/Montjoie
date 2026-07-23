#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

void CountLateralFaces(const SurfacicMesh<Dimension3>& surf_mesh0, int i, 
                       const Mesh<Dimension3>& mesh, Vector<int>& RefEdge, int& nb_lateral_edge)
{
  // finding lateral faces
  int ref = surf_mesh0.BoundaryRef(i).GetReference();
  for (int j = 0; j < 4; j++)
    {
      int nf = surf_mesh0.ListeBoundaries(i);
      int ne0 = surf_mesh0.BoundaryRef(i).numEdge(j);
      if (RefEdge(ne0) == 0)
        {
          int ne = mesh.BoundaryRef(nf).numEdge(j);
          int refb = 0;
          for (int k = 0; k < mesh.GetEdge(ne).GetNbFaces(); k++)
            {
              int nf2 = mesh.GetEdge(ne).numFace(k);
              int ref2 = mesh.Boundary(nf2).GetReference();
              if ( (ref2 > 0) && (ref2 != ref))
                {
                  refb = ref2;
                  break;
                }
            }
          
          if (refb > 0)
            {
              RefEdge(ne0) = refb;
              nb_lateral_edge++;
            }
          else
            RefEdge(ne0) = -1;
        }
    }
}


void AddSphericalLayer(Mesh<Dimension3>& mesh, const SurfacicMesh<Dimension3>& surf_mesh0,
		       const SurfacicMesh<Dimension3>& surf_mesh, const Vector<IVect>& NodleSurf,
		       const Real_wp& r0, const Real_wp& r1, int r, int ref_domain, bool keep_ref, int ref_top)
{
  for (int i = 0; i < surf_mesh.GetNbBoundaryRef(); i++)
    if (surf_mesh.BoundaryRef(i).GetNbVertices() != 4)
      {
	cout << "Implemented only for hexahedral meshes" << endl;
	abort();
      }

  if (r == 1)
    {
      int nb_old_vert = mesh.GetNbVertices();
      int nb_pts = surf_mesh.GetNbVertices();
      
      mesh.ResizeVertices(nb_old_vert + nb_pts);
      
      // adding vertices
      int nb = nb_old_vert;
      for (int i = 0; i < nb_pts; i++)
	mesh.Vertex(nb++) = surf_mesh.Vertex(i);
      
      // counting new faces
      int nb_lateral_edge = 0;
      Vector<int> RefEdge(surf_mesh0.GetNbEdges()); RefEdge.Zero();
      for (int i = 0; i < surf_mesh.GetNbBoundaryRef(); i++)
        CountLateralFaces(surf_mesh0, i, mesh, RefEdge, nb_lateral_edge);
      
      // then adding elements and faces
      int nb_old_face = mesh.GetNbBoundaryRef();
      int offset_face = nb_old_face;
      if (keep_ref)
	{
          mesh.ResizeBoundariesRef(nb_old_face + surf_mesh.GetNbBoundaryRef() + nb_lateral_edge);
          offset_face += surf_mesh.GetNbBoundaryRef();
        }
      else if (nb_lateral_edge > 0)
        mesh.ResizeBoundariesRef(nb_old_face + nb_lateral_edge);

      int nb_old_elt = mesh.GetNbElt();
      mesh.ResizeElements(nb_old_elt + surf_mesh.GetNbBoundaryRef());
      int num = nb_old_elt;
      for (int i = 0; i < surf_mesh.GetNbBoundaryRef(); i++)
	{
	  int n1 = nb_old_vert + surf_mesh.BoundaryRef(i).numVertex(0);
	  int n2 = nb_old_vert + surf_mesh.BoundaryRef(i).numVertex(1);
	  int n3 = nb_old_vert + surf_mesh.BoundaryRef(i).numVertex(2);
	  int n4 = nb_old_vert + surf_mesh.BoundaryRef(i).numVertex(3);
	  int nv1 = surf_mesh.ListeVertices(n1-nb_old_vert);
	  int nv2 = surf_mesh.ListeVertices(n2-nb_old_vert);
	  int nv3 = surf_mesh.ListeVertices(n3-nb_old_vert);
	  int nv4 = surf_mesh.ListeVertices(n4-nb_old_vert);
	  mesh.Element(num++).InitHexahedral(nv1, nv2, nv3, nv4, n1, n2, n3, n4, ref_domain);
          
	  int nf = surf_mesh.ListeBoundaries(i);
	  if (keep_ref)
	    nf = nb_old_face + i;
          
	  mesh.BoundaryRef(nf).InitQuadrangular(n1, n2, n3, n4, ref_top);
          
	  int ne = surf_mesh0.Boundary(i).numEdge(0);
	  if (RefEdge(ne) > 0)
            {
              int ref2 = RefEdge(ne);
              mesh.BoundaryRef(offset_face++).InitQuadrangular(nv1, nv2, n2, n1, ref2);
              RefEdge(ne) = 0;
            }

	  ne = surf_mesh0.Boundary(i).numEdge(1);
	  if (RefEdge(ne) > 0)
            {
              int ref2 = RefEdge(ne);
              mesh.BoundaryRef(offset_face++).InitQuadrangular(nv2, nv3, n3, n2, ref2);
              RefEdge(ne) = 0;
            }

	  ne = surf_mesh0.Boundary(i).numEdge(2);
	  if (RefEdge(ne) > 0)
            {
              int ref2 = RefEdge(ne);
              mesh.BoundaryRef(offset_face++).InitQuadrangular(nv3, nv4, n4, n3, ref2);
              RefEdge(ne) = 0;
            }

	  ne = surf_mesh0.Boundary(i).numEdge(3);
	  if (RefEdge(ne) > 0)
            {
              int ref2 = RefEdge(ne);
              mesh.BoundaryRef(offset_face++).InitQuadrangular(nv1, nv4, n4, n1, ref2);
              RefEdge(ne) = 0;
            }
	}	
    }
  else if (r == 2)
    {
      cout << "Not implemented" << endl;
      abort();
    }
  else if (r == 3)
    {
      Matrix<int> NumNodes2D, CoordinateNodes;
      MeshNumbering<Dimension2>::
	ConstructQuadrilateralNumbering(r, NumNodes2D, CoordinateNodes);

      int nb_old_vert = mesh.GetNbVertices();
      int nb_pts = surf_mesh.GetNbVertices();
      int nb_edges = surf_mesh0.GetNbEdges();
      int nb_face0 = surf_mesh0.GetNbBoundaryRef();
      
      mesh.ResizeVertices(nb_old_vert + nb_pts + 2*nb_edges + 4*nb_face0);
      
      // adding vertices on the refined surface
      int nb = nb_old_vert;
      for (int i = 0; i < nb_pts; i++)
	mesh.Vertex(nb++) = surf_mesh.Vertex(i);

      // on lateral faces
      Real_wp rm = Real_wp(2)*r0/3 + r1/3;
      int offset = nb;
      for (int i = 0; i < nb_edges; i++)
	{
	  int n0 = surf_mesh0.GetEdge(i).numVertex(0);
	  int n1 = surf_mesh0.GetEdge(i).numVertex(1);
	  R3 pt0 = surf_mesh0.Vertex(n0);
	  R3 pt1 = surf_mesh0.Vertex(n1);
	  R3 ptM = Real_wp(2)/3*pt0 + Real_wp(1)/3*pt1;
	  ptM *= rm / Norm2(ptM);
	  mesh.Vertex(nb++) = ptM;
	  ptM = Real_wp(1)/3*pt0 + Real_wp(2)/3*pt1;
	  ptM *= rm / Norm2(ptM);
	  mesh.Vertex(nb++) = ptM;	  
	}

      // in the volume
      int nb_lateral_edge = 0;
      Vector<int> RefEdge(surf_mesh0.GetNbEdges()); RefEdge.Zero();
      rm = Real_wp(2)/3 * r1 + r0/3;
      int offset_f = nb;
      for (int i = 0; i < nb_face0; i++)
	{
	  R3 pt1 = surf_mesh0.Vertex(surf_mesh0.BoundaryRef(i).numVertex(0));
	  R3 pt2 = surf_mesh0.Vertex(surf_mesh0.BoundaryRef(i).numVertex(1));
	  R3 pt3 = surf_mesh0.Vertex(surf_mesh0.BoundaryRef(i).numVertex(2));
	  R3 pt4 = surf_mesh0.Vertex(surf_mesh0.BoundaryRef(i).numVertex(3));
	  R3 ptM = Real_wp(4)/9*pt1 + Real_wp(2)/9*pt2 + Real_wp(1)/9*pt3 + Real_wp(2)/9*pt4;
	  ptM *= rm / Norm2(ptM);
	  mesh.Vertex(nb++) = ptM;

	  ptM = Real_wp(2)/9*pt1 + Real_wp(4)/9*pt2 + Real_wp(2)/9*pt3 + Real_wp(1)/9*pt4;
	  ptM *= rm / Norm2(ptM);
	  mesh.Vertex(nb++) = ptM;

	  ptM = Real_wp(1)/9*pt1 + Real_wp(2)/9*pt2 + Real_wp(4)/9*pt3 + Real_wp(2)/9*pt4;
	  ptM *= rm / Norm2(ptM);
	  mesh.Vertex(nb++) = ptM;

	  ptM = Real_wp(2)/9*pt1 + Real_wp(1)/9*pt2 + Real_wp(2)/9*pt3 + Real_wp(4)/9*pt4;
	  ptM *= rm / Norm2(ptM);
	  mesh.Vertex(nb++) = ptM;

          CountLateralFaces(surf_mesh0, i, mesh, RefEdge, nb_lateral_edge);
	}

      // then adding elements and faces
      int nb_old_elt = mesh.GetNbElt();
      int nb_new_elt = nb_face0*13;
      mesh.ResizeElements(nb_old_elt + nb_new_elt);

      int nb_old_face = mesh.GetNbBoundaryRef();
      int offset_face = nb_old_face + 9*nb_face0;
      if (keep_ref)
	mesh.ResizeBoundariesRef(nb_old_face + 9*nb_face0 + 4*nb_lateral_edge);
      else
	{
	  mesh.ResizeBoundariesRef(nb_old_face + 8*nb_face0 + 4*nb_lateral_edge);
	  offset_face = nb_old_face + 8*nb_face0;
	}
      
      int num = nb_old_elt;
      for (int i = 0; i < nb_face0; i++)
	{
	  int n1 = nb_old_vert + NodleSurf(i)(NumNodes2D(0, 0));
	  int n2 = nb_old_vert + NodleSurf(i)(NumNodes2D(1, 0));
	  int n3 = nb_old_vert + NodleSurf(i)(NumNodes2D(2, 0));
	  int n4 = nb_old_vert + NodleSurf(i)(NumNodes2D(3, 0));
	  int n5 = nb_old_vert + NodleSurf(i)(NumNodes2D(0, 1));
	  int n6 = nb_old_vert + NodleSurf(i)(NumNodes2D(1, 1));
	  int n7 = nb_old_vert + NodleSurf(i)(NumNodes2D(2, 1));
	  int n8 = nb_old_vert + NodleSurf(i)(NumNodes2D(3, 1));
	  int n9 = nb_old_vert + NodleSurf(i)(NumNodes2D(0, 2));
	  int n10 = nb_old_vert + NodleSurf(i)(NumNodes2D(1, 2));
	  int n11 = nb_old_vert + NodleSurf(i)(NumNodes2D(2, 2));
	  int n12 = nb_old_vert + NodleSurf(i)(NumNodes2D(3, 2));
	  int n13 = nb_old_vert + NodleSurf(i)(NumNodes2D(0, 3));
	  int n14 = nb_old_vert + NodleSurf(i)(NumNodes2D(1, 3));
	  int n15 = nb_old_vert + NodleSurf(i)(NumNodes2D(2, 3));
	  int n16 = nb_old_vert + NodleSurf(i)(NumNodes2D(3, 3));
	  
	  int n17, n18, n19, n20, n21, n22, n23, n24;
	  int ne = surf_mesh0.Boundary(i).numEdge(0);
	  if (surf_mesh0.Boundary(i).GetOrientationEdge(0))
	    {
	      n17 = offset + 2*ne;
	      n18 = offset + 2*ne+1;
	    }
	  else
	    {
	      n17 = offset + 2*ne+1;
	      n18 = offset + 2*ne;
	    }

	  ne = surf_mesh0.Boundary(i).numEdge(1);
	  if (surf_mesh0.Boundary(i).GetOrientationEdge(1))
	    {
	      n21 = offset + 2*ne;
	      n22 = offset + 2*ne+1;
	    }
	  else
	    {
	      n21 = offset + 2*ne+1;
	      n22 = offset + 2*ne;
	    }

	  ne = surf_mesh0.Boundary(i).numEdge(2);
	  if (surf_mesh0.Boundary(i).GetOrientationEdge(2))
	    {
	      n24 = offset + 2*ne;
	      n23 = offset + 2*ne+1;
	    }
	  else
	    {
	      n24 = offset + 2*ne+1;
	      n23 = offset + 2*ne;
	    }

	  ne = surf_mesh0.Boundary(i).numEdge(3);
	  if (surf_mesh0.Boundary(i).GetOrientationEdge(3))
	    {
	      n20 = offset + 2*ne;
	      n19 = offset + 2*ne+1;
	    }
	  else
	    {
	      n20 = offset + 2*ne+1;
	      n19 = offset + 2*ne;
	    }

	  int n25 = surf_mesh0.ListeVertices(surf_mesh0.BoundaryRef(i).numVertex(0));
	  int n26 = surf_mesh0.ListeVertices(surf_mesh0.BoundaryRef(i).numVertex(1));
	  int n27 = surf_mesh0.ListeVertices(surf_mesh0.BoundaryRef(i).numVertex(2));
	  int n28 = surf_mesh0.ListeVertices(surf_mesh0.BoundaryRef(i).numVertex(3));

	  int n29 = offset_f + 4*i;
	  int n30 = offset_f + 4*i+1;
	  int n31 = offset_f + 4*i+2;
	  int n32 = offset_f + 4*i+3;

	  /*DISP(n1); DISP(n2); DISP(n3); DISP(n4); DISP(n5); DISP(n6); DISP(n7); DISP(n8);
	  DISP(n9); DISP(n10); DISP(n11); DISP(n12); DISP(n13); DISP(n14); DISP(n15); DISP(n16);
	  DISP(n17); DISP(n18); DISP(n19); DISP(n20); DISP(n21); DISP(n22); DISP(n23); DISP(n24);
	  DISP(n25); DISP(n26); DISP(n27); DISP(n28); DISP(n29); DISP(n30); DISP(n31); DISP(n32);*/
	  
	  mesh.Element(num++).InitHexahedral(n1, n2, n6, n5, n25, n17, n29, n19, ref_domain);
	  mesh.Element(num++).InitHexahedral(n2, n3, n7, n6, n17, n18, n30, n29, ref_domain);
	  mesh.Element(num++).InitHexahedral(n3, n4, n8, n7, n18, n26, n21, n30, ref_domain);
	  mesh.Element(num++).InitHexahedral(n5, n6, n10, n9, n19, n29, n32, n20, ref_domain);
	  mesh.Element(num++).InitHexahedral(n6, n7, n11, n10, n29, n30, n31, n32, ref_domain);
	  mesh.Element(num++).InitHexahedral(n7, n8, n12, n11, n30, n21, n22, n31, ref_domain);
	  mesh.Element(num++).InitHexahedral(n9, n10, n14, n13, n20, n32, n23, n28, ref_domain);
	  mesh.Element(num++).InitHexahedral(n10, n11, n15, n14, n32, n31, n24, n23, ref_domain);
	  mesh.Element(num++).InitHexahedral(n11, n12, n16, n15, n31, n22, n27, n24, ref_domain);
	  mesh.Element(num++).InitHexahedral(n17, n18, n30, n29, n25, n26, n21, n19, ref_domain);
	  mesh.Element(num++).InitHexahedral(n29, n30, n31, n32, n19, n21, n22, n20, ref_domain);
	  mesh.Element(num++).InitHexahedral(n32, n31, n24, n23, n20, n22, n27, n28, ref_domain);
	  mesh.Element(num++).InitHexahedral(n19, n21, n22, n20, n25, n26, n27, n28, ref_domain);

	  int nf = surf_mesh0.ListeBoundaries(i);
	  if (keep_ref)
	    nf = nb_old_face + 9*i;
	  
	  mesh.BoundaryRef(nf).InitQuadrangular(n1, n2, n6, n5, ref_top);
	  if (keep_ref)
	    nf++;
	  else
	    nf = nb_old_face + 8*i;
	  
	  mesh.BoundaryRef(nf++).InitQuadrangular(n2, n3, n7, n6, ref_top);
	  mesh.BoundaryRef(nf++).InitQuadrangular(n3, n4, n8, n7, ref_top);
	  mesh.BoundaryRef(nf++).InitQuadrangular(n5, n6, n10, n9, ref_top);
	  mesh.BoundaryRef(nf++).InitQuadrangular(n6, n7, n11, n10, ref_top);
	  mesh.BoundaryRef(nf++).InitQuadrangular(n7, n8, n12, n11, ref_top);
	  mesh.BoundaryRef(nf++).InitQuadrangular(n9, n10, n14, n13, ref_top);
	  mesh.BoundaryRef(nf++).InitQuadrangular(n10, n11, n15, n14, ref_top);
	  mesh.BoundaryRef(nf++).InitQuadrangular(n11, n12, n16, n15, ref_top);

	  ne = surf_mesh0.Boundary(i).numEdge(0);
	  if (RefEdge(ne) > 0)
	    {
	      int ref2 = RefEdge(ne);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n1, n2, n17, n25, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n2, n3, n18, n17, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n3, n4, n26, n18, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n17, n18, n26, n25, ref2);
	      RefEdge(ne) = 0;
	    }

	  ne = surf_mesh0.Boundary(i).numEdge(1);
	  if (RefEdge(ne) > 0)
	    {
	      int ref2 = RefEdge(ne);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n4, n8, n21, n26, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n8, n12, n22, n21, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n12, n16, n27, n22, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n21, n22, n27, n26, ref2);
	      RefEdge(ne) = 0;
	    }

	  ne = surf_mesh0.Boundary(i).numEdge(2);
	  if (RefEdge(ne) > 0)
	    {
	      int ref2 = RefEdge(ne);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n13, n14, n23, n28, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n14, n15, n24, n23, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n15, n16, n27, n24, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n23, n24, n27, n28, ref2);
	      RefEdge(ne) = 0;
	    }

	  ne = surf_mesh0.Boundary(i).numEdge(3);
	  if (RefEdge(ne) > 0)
	    {
	      int ref2 = RefEdge(ne);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n1, n5, n19, n25, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n5, n9, n20, n19, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n9, n13, n28, n20, ref2);
	      mesh.BoundaryRef(offset_face++).InitQuadrangular(n19, n20, n28, n25, ref2);
	      RefEdge(ne) = 0;
	    }
	}
    }
  else
    {
      cout << "not implemented" << endl;
      abort();
    }

  mesh.SortBoundariesRef();
  mesh.ReorientElements();
  mesh.FindConnectivity();
  mesh.AddBoundaryFaces();
}

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc < 6)
    {
      cout<<"Cette commande demande au moins cinq arguments" << endl;
      cout<<"add_spherical_layer initial.mesh final.mesh ref file_r_teta radius" << endl;
      cout<<"the file final.mesh is created" << endl;
      return -1;
    }
  
  string file_input(argv[1]);
  string file_output(argv[2]);
  int ref = atoi(argv[3]);
  string file_r(argv[4]);
  Real_wp radius = to_num<Real_wp>(argv[5]);
  
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

  Mesh<Dimension3> mesh;
  mesh.Read(file_input);

  // finding an unused reference
  int nb_max = mesh.GetNbReferences()+1;
  for (int i = 0; i < ref_top.GetM(); i++)
    nb_max = max(nb_max, ref_top(i));
  
  Vector<bool> RefUsed(nb_max+2);
  RefUsed.Fill(false);

  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    RefUsed(mesh.BoundaryRef(i).GetReference()) = true;

  for (int i = 0; i < ref_top.GetM(); i++)
    RefUsed(ref_top(i)) = true;

  int ref_temp = -1;
  for (int i = 1; i < RefUsed.GetM(); i++)
    if (!RefUsed(i))
      {
	ref_temp = i;
	break;
      }

  mesh.print_level = 2;

  // enforcing a sphere on the given reference
  mesh.SetCurveType(ref, mesh.CURVE_SPHERE);
  VectReal_wp param(4);
  cout << "Radius is equal to " << radius << endl;
  param(0) = Real_wp(0); param(1) = Real_wp(0); param(2) = Real_wp(0);
  param(3) = radius;
  mesh.SetCurveParameter(ref, param);
  mesh.ProjectPointsOnCurves();

  SurfacicMesh<Dimension3> surf_mesh, surf_mesh_init;
  IVect ref_cond(mesh.GetNbReferences()+1); ref_cond.Fill();

  VectReal_wp step_subdiv; 
  Vector<IVect> Nodle; Vector<VectR3> points_div;
  Vector<IVect> NodleSurf; Vector<VectR2> points_surf;    
    
  // loops on layers to add
  int nb_layers = step_r.GetM();
  int raff_global = 1;
  int ref_domain = 1, ref_surf = ref;  
  for (int p = 0; p < nb_layers; p++)
    {
      mesh.GetBoundaryMesh(ref_surf, surf_mesh_init, ref_cond, 0);
      
      Real_wp r2 = step_r(p);
      // int raff = toInteger(ceil(step_teta(p)/raff_global));
      // we select between 1 and 3
      int raff = 1;
      if (step_teta(p)/raff_global > 1.7)
	raff = 3;
      
      raff_global *= raff;
      cout << "At radius r = " << r2 << " theta-refinement = " << raff << endl;
      
      // the surface mesh is refined
      surf_mesh = surf_mesh_init;
      if (raff > 1)
	{
	  Linspace(Real_wp(0), Real_wp(1), 4, step_subdiv);
	  surf_mesh.SubdivideMesh(step_subdiv, points_surf, points_div, Nodle, NodleSurf);
	}

      // then vertices are scaled to fit the new radius
      Real_wp ratio = r2 / radius;
      surf_mesh.Vertex() *= ratio;
      
      // surf_mesh.Write("surf.mesh");

      int new_ref = ref_top(p);
      if (new_ref <= 0)
	new_ref = ref_temp;

      bool keep_ref = true;
      if (ref_surf == ref_temp)
	keep_ref = false;
      
      AddSphericalLayer(mesh, surf_mesh_init, surf_mesh, NodleSurf,
			radius, r2, raff, ref_domain, keep_ref, new_ref);

      mesh.SetCurveType(new_ref, mesh.CURVE_SPHERE);
      param(3) = r2; mesh.SetCurveParameter(new_ref, param);
      mesh.ProjectPointsOnCurves();
      
      //mesh.Write("vol.mesh");
      
      //cout << "waiting" << endl; int test_input; cin >> test_input;

      ref_surf = new_ref;
      radius = r2;
    }
  
  mesh.Write(file_output);
}  
  
