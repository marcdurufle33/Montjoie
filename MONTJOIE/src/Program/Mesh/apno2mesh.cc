#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);

  R2::threshold = 1e-6;
  
  if (argc < 3)
    {
      cout<<"Cette commande demande deux arguments"<<endl;
      cout<<"apno2mesh test.msh toto.mesh"<<endl;
      cout<<"the file toto.mesh is created"<<endl;
      return -1;
    }
  
  string file_input(argv[1]);
  string file_output(argv[2]);
  int order = 1;
  if (argc > 3)
    order = atoi(argv[3]);
  
  Mesh<Dimension2> mesh, sub_mesh, res_mesh;
  
  // on lit le maillage
  mesh.SetGeometryOrder(order);
  mesh.Read(file_input);

  // sorting vertices
  mesh.RemoveDuplicateVertices();
  
  int ref_max = 0;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    ref_max = max(ref_max, mesh.Element(i).GetReference());

  // on enleve des references surfaciques
  IVect ref_surface_to_remove;
  
  // on change les references volumiques
  IVect new_ref_volume(ref_max+1);
  new_ref_volume.Fill();
  cout << "Nombre de references surfaciques a changer ? " << endl;
  int n = 0; cin >> n;
  for (int i = 0; i < n; i++)
    {
      cout << "Reference " << i << endl;
      cout << "Entrez l'ancien et le nouveau numero a la suite " << endl;
      int ref1, ref2;
      cin >> ref1 >> ref2;
      if ((ref1 > 0)&&(ref2 > 0))
        {
          new_ref_volume(ref1) = ref2;
        }
    }
  
  // boucle sur toutes les references
  Vector<bool> VertexOnSubdomain(mesh.GetNbVertices());
  Vector<bool> ElementOnSubdomain(mesh.GetNbElt());
  for (int ref = 1; ref <= ref_max; ref++)
    {
      int nb_elt = 0;
      int nb_vertices = 0;
      VertexOnSubdomain.Fill(false);
      ElementOnSubdomain.Fill(false);
      for (int i = 0; i < mesh.GetNbElt(); i++)
        if (mesh.Element(i).GetReference() == ref)
          {
            ElementOnSubdomain(i) = true;
            for (int j = 0; j < mesh.Element(i).GetNbVertices(); j++)
              {
                int nv = mesh.Element(i).numVertex(j);
                if (!VertexOnSubdomain(nv))
                  {
                    VertexOnSubdomain(nv) = true;
                    nb_vertices++;
                  }
              }
            nb_elt++;
          }
      
      if (nb_elt > 0)
        {
          mesh.CreateSubmesh(sub_mesh, nb_vertices, nb_elt,
                             VertexOnSubdomain, ElementOnSubdomain);
          
          sub_mesh.AddBoundaryEdges();
          sub_mesh.FindConnectivity();
          
          // initial points
          int n1 = -1, n2 = -1;
          for (int i = 0; i < sub_mesh.GetNbBoundaryRef(); i++)
            if (sub_mesh.BoundaryRef(i).GetNbElements() == 1)
              {
                n1 = sub_mesh.BoundaryRef(i).numVertex(0);
                n2 = sub_mesh.BoundaryRef(i).numVertex(1);
                int num_elem = sub_mesh.BoundaryRef(i).numElement(0);
                int num_loc = sub_mesh.Element(num_elem).GetPositionBoundary(i);
                if (!sub_mesh.Element(num_elem).GetOrientationEdge(num_loc))
                  {
                    int n3 = n2;
                    n2 = n1;
                    n1 = n3;
                  }
                
                break;
              }
          
          Vector<bool> edge_used(sub_mesh.GetNbBoundary());
          edge_used.Fill(true);
          for (int i = 0; i < sub_mesh.GetNbBoundary(); i++)
            if (sub_mesh.Boundary(i).GetNbElements() == 1)
              edge_used(i) = false;
          
          // finding points on this curve
          IVect ref_cond(sub_mesh.GetNbReferences()+1);
          ref_cond(0) = 0;
          ref_cond.Fill(1);
          int num_init = n1;
          IVect VertexBound(sub_mesh.GetNbVertices());
          int nb_vert = 1, ne;
          VertexBound(0) = n1;
          while (n2 != num_init)
            {
              VertexBound(nb_vert) = n2;
              nb_vert++;
              
              int inext = sub_mesh.FindFollowingVertex(n1, n2, 1, ref_cond, edge_used, ne);
              n1 = n2;
              n2 = inext;              
            }

          ofstream file_mesh("SurfTest.mesh"); file_mesh.precision(15);
          file_mesh << "MeshVersionFormatted 1"<<endl;
          file_mesh << "Dimension"<<endl;
          file_mesh << "2"<<endl;
          file_mesh << "Vertices"<<endl;
          file_mesh << nb_vert << endl;
          for (int i = 0; i < nb_vert; i++)
            {
              int nv = VertexBound(i);
              Real_wp x  = sub_mesh.Vertex(nv)(0);
              Real_wp y  = sub_mesh.Vertex(nv)(1);
              file_mesh << x << "  " << y << " 1 " << endl;
            }          
      
          file_mesh << "Edges"<<endl;
          file_mesh << nb_vert << endl;
          for (int i = 0; i < nb_vert-1; i++)
            file_mesh << i+1 << " " << i+2 << " 1 " <<endl;
          
          file_mesh << nb_vert << " 1 1 " <<endl;
          file_mesh << "End"<<endl;
          file_mesh.close();
          
          // creating data file for apnoxx
          ofstream file_apno("TestApno");
          file_apno.precision(14);
          file_apno << "'TestApno  '" << endl;
          file_apno << "'POIN    '"<<endl;
          file_apno << "    1 " << nb_vert << "  $ IMPRE NPOINT  $" << endl;
          file_apno << "$   NOP   NOREF(NOP)     X(NOP).            Y(NOP).  $" << endl;
          for (int i = 0; i < nb_vert; i++)
            {
              int nv = VertexBound(i);
              Real_wp x  = sub_mesh.Vertex(nv)(0);
              Real_wp y  = sub_mesh.Vertex(nv)(1);
              file_apno << i+1 << "  1 " << x << "  " << y << endl;
            }
          
          n1 = nb_vert/4;
          n2 = 2*nb_vert/4;
          int n3 = 3*nb_vert/4;
          
          file_apno << "'LIGN    '"<<endl;
          file_apno << "  1  4    $ IMPRE  NDLM   $ " << endl;
          file_apno << "$   NOLIG NOELIG NEXTR1 NEXTR2 NOREFL NFFRON       RAISON $" << endl;
          file_apno << "   1 "<< n1+1 << "  1  " << n1+1 << "  1   -1    1.0 " << endl;
          for (int i = 1; i < n1; i++)
            file_apno << "  " << i+1 << "  $NOEL " << endl;
          
          file_apno << "   2 "<< n2-n1+1 << "  " << n1+1 <<"  " << n2+1 << "  1   -1    1.0 " << endl;
          for (int i = n1+1; i < n2; i++)
            file_apno << "  " << i+1 << "  $NOEL " << endl;

          file_apno << "   3 "<< n3-n2+1 << "  " << n2+1 <<"  " << n3+1 << "  1   -1    1.0 " << endl;
          for (int i = n2+1; i < n3; i++)
            file_apno << "  " << i+1 << "  $NOEL " << endl;

          file_apno << "   4 "<< nb_vert-n3+1 << "  " << n3+1 <<"  " << " 1 " << "  1   -1    1.0 " << endl;
          for (int i = n3+1; i < nb_vert; i++)
            file_apno << "  " << i+1 << "  $NOEL " << endl;
          
          file_apno << "'TRIA   '" << endl;
          file_apno << "  1     1    1     4     1          $ IMPRE NIVEAU NUDSD NBRELI NS1L" << endl;
          file_apno << "$ LISTE DES LIGNES DU CONTOUR :" << endl;
          file_apno << "   1  2   3   4 " << endl;
          file_apno << "   1     0     21                      $ NCOMP NBRINT IOPT $" << endl;
          file_apno << "   4                                  $ COMPOSANTE" << endl;
          
          file_apno<< "'SAUV       '" << endl;
          file_apno << " 1     1     0                     $ IMPRE NINOPO NTNOPO" << endl;
          file_apno << "TestApno.nopo                        $ NOM FICHIER  " << endl;
          file_apno << "'FIN   '" << endl;
          
          file_apno.close();
          
          cout << "Executez apno2mesh - tapez 1 quand c'est termine " << endl;
          int test_input; cin >> test_input;
          
          sub_mesh.Read("TestApno.mesh");
          int refb = new_ref_volume(ref);
          for (int i = 0; i < sub_mesh.GetNbElt(); i++)
            sub_mesh.Element(i).SetReference(refb);
          
          res_mesh.AppendMesh(sub_mesh);
        }
    }

  // sorting vertices
  VectR2 VertInit;
  VectR2 VertRes;
  IVect original_to_new(mesh.GetNbVertices());
  original_to_new.Fill(-1);
  
  VertInit = mesh.Vertex();
  VertRes = res_mesh.Vertex();
  
  // searching coincident points
  int k = 0;
  for (int i = 0; i < VertInit.GetM(); i++)
    {
      while ( (k < VertRes.GetM()) && (VertRes(k) < VertInit(i)))
        k++;
      
      if (k < VertRes.GetM())
        if (VertRes(k) == VertInit(i))
          original_to_new(i) = k;
    }
  
  // then searching coincident edges
  res_mesh.SetGeometryOrder(order);
  res_mesh.ProjectPointsOnCurves();
  Vector<bool> ref_surf_to_remove(mesh.GetNbReferences()+1);
  ref_surf_to_remove.Fill(false);
  for (int i = 0; i < ref_surface_to_remove.GetM(); i++)
    {
      int ref = ref_surface_to_remove(i);
      if ((ref > 0) && (ref < ref_surf_to_remove.GetM()))
        ref_surf_to_remove(ref) = true;
    }
    
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      if (ref > 0)
        if (!ref_surf_to_remove(ref))
          {
            int n1 = mesh.BoundaryRef(i).numVertex(0);
            int n2 = mesh.BoundaryRef(i).numVertex(1);
            int k1 = original_to_new(n1);
            int k2 = original_to_new(n2);
            if ((k1 >= 0) && (k2 >= 0))
              {
                int ne = res_mesh.FindEdgeRefWithExtremities(k1, k2);
                // DISP(i); DISP(n1); DISP(n2); DISP(k1); DISP(k2); DISP(ne); DISP(ref);
                if (ne >= 0)
                  {
                    res_mesh.BoundaryRef(ne).SetReference(ref);
                    res_mesh.Boundary(ne).SetReference(ref);
                    res_mesh.SetCurveType(ref, res_mesh.CURVE_FILE);
                    for (int k = 0; k < order-1; k++)
                      res_mesh.SetPointInsideEdge(ne, k, mesh.GetPointInsideEdge(i, k));
                  }
              }
          }
    }
  
  res_mesh.ProjectPointsOnCurves();
  res_mesh.SplitIntoQuadrilaterals();
      
  // on ecrit maillage
  res_mesh.Write(file_output);
  
  return FinalizeMontjoie();
}  
