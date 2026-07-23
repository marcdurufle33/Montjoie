#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

void Permute(Mesh<Dimension2>& mesh)
{
  int ref = 1;
  
  // first step : we find an initial point 
  IVect NbEdgesVertices(mesh.GetNbVertices());
  NbEdgesVertices.Fill(0);
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    if (mesh.BoundaryRef(i).GetReference() == ref)
      {
        NbEdgesVertices(mesh.BoundaryRef(i).numVertex(0))++;
        NbEdgesVertices(mesh.BoundaryRef(i).numVertex(1))++;
      }
  
  int num_init = -1;
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    if (NbEdgesVertices(i) > 0)
      {
        if (num_init == -1)
          num_init = i;
        else if (NbEdgesVertices(i) == 1)
          num_init = i;
        else if (NbEdgesVertices(i) >= 3)
          {
            cout << "Curve with a branch"<<endl;
            abort();
          }
      }
  
  IVect Point_curve(mesh.GetNbBoundaryRef()+1), Edge_curve(mesh.GetNbBoundaryRef());
  IVect ref_cond(mesh.GetNbReferences()); ref_cond.Fill();
  int ipoint = num_init, iback = -1, inext = -1, iedge = -1, nb_points = 0;
  int nb_edges = 0;
  while ((nb_points == 0)||((ipoint != num_init)&&(ipoint != -1)))
    {
      if (iedge != -1)
        {
          Edge_curve(nb_edges) = iedge;
          nb_edges++;
        }
      
      Point_curve(nb_points) = ipoint;
      inext = mesh.FindFollowingVertex(iback, ipoint, ref, ref_cond, iedge);
      iback = ipoint; ipoint = inext;
      nb_points++;
    }
  
  VectR2 OldVertex(mesh.GetNbVertices());
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    OldVertex(i) = mesh.Vertex(i);
  
  IVect OldFirstExt(mesh.GetNbBoundaryRef());
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    OldFirstExt = mesh.BoundaryRef(i).numVertex(0);
  
  int order = mesh.GetGeometryOrder();
  Matrix<R2> OldPoints(mesh.GetNbBoundaryRef(), order-1);
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    for (int j = 0; j < order-1; j++)
      OldPoints(i, j) = mesh.GetPointInsideEdge(i, j);
  
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    mesh.Vertex(i) = OldVertex(Point_curve(i));
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ne = Edge_curve(i);
      int nv = OldFirstExt(ne);
      mesh.BoundaryRef(i).Init(i, i+1, ref);
      if (nv == Point_curve(i))
        for (int j = 0; j < order-1; j++)
          mesh.SetPointInsideEdge(i, j, OldPoints(ne, j));
      else
        for (int j = 0; j < order-1; j++)
          mesh.SetPointInsideEdge(i, j, OldPoints(ne, order-2-j));
    }
}

void AddPoint(Vector<int>& num_pt, R2& normale, VectR2& Points,
	      int n, const R2& pt)
{
  // checking if the point already exists
  for (int k = 0; k < num_pt.GetM(); k++)
    if (n == num_pt(k))
      return;
  
  // if it doesn't exist, we add it to the correct position
  Real_wp scal_pt = DotProd(pt-Points(0), normale);
  int N = num_pt.GetM(); int pos = N;
  for (int k = 1; k < N; k++)
    {
      R2 vec_u = Points(k) - Points(0);
      Real_wp scal_k = DotProd(vec_u, normale);
      if (scal_pt < scal_k)
	{
	  pos = k;
	  break;
	}
    }
  
  Vector<int> old_num(num_pt); VectR2 OldPoints(Points);
  num_pt.Reallocate(N+1); Points.Reallocate(N+1);
  for (int k = 0; k < pos; k++)
    {
      num_pt(k) = old_num(k);
      Points(k) = OldPoints(k);
    }
  
  num_pt(pos) = n;
  Points(pos) = pt;
  
  for (int k = pos; k < N; k++)
    {
      num_pt(k+1) = old_num(k);
      Points(k+1) = OldPoints(k);
    }
}

int ChangeNum(int n, int offset)
{
  if (n < 0)
    return -n;
  else
    return offset+n;
}

int main(int argc, char** argv)
{  
  InitMontjoie(argc, argv);
  
  if (argc < 4)
    {
      cout << "Entrez deux noms de fichiers et l'epaisseur " << endl;
      cout << "Exemple : ./intersec.x CourbeChevalet.mesh Raidisseur.mesh 0.5" << endl;
      abort();
    }
  
  string input_curve(argv[1]);
  string input_mesh(argv[2]);
  
  int ref = 1;
  bool add_mesh_rib = true;
  Real_wp thickness = atof(argv[3]);
  int order = 4;
  Mesh<Dimension2> mesh_curve, mesh;
  mesh_curve.SetGeometryOrder(order);
  
  mesh_curve.Read(input_curve);
  mesh.Read(input_mesh);
  
  //VectReal_wp step(5); step.Fill(); Mlt(0.25, step);
  //mesh_curve.SubdivideMesh(step);
  
  //mesh_curve.Write("refined.mesh");
  
  Permute(mesh_curve);

  // boucle sur tous les segments du chevalet
  Real_wp threshold = 1e-12;
  Globatto<Real_wp> lob = mesh_curve.GetCurveFunctions1D();
  lob.ComputeGradPhi();
  ofstream file_out("toto.geo");
  file_out.precision(15);
  file_out << "Mesh.RecombineAll = 1;" << endl;
  file_out << "lc = 1.0;\n" << endl;
  file_out << "scalex = 1.0;\n";
  file_out << "scaley = 1.0;\n" << endl;
  //file_out << "scalex = scalex/10.674;\n";
  //file_out << "scaley = scaley/10.674;\n";
  file_out << "nmpt = 0; \n nmline = 0; nmsurf = 0;\n"<< endl;  
  
  int nmpt = 1, nmline = 1, nmsurf = 1;
  int nb_intersection = 0;
  Vector<TinyVector<int, 2> > liste_point_joint(1000), list_edge_joint(1000), list_spline_joint(1000);
  VectReal_wp length_spline(1000); length_spline.Fill(0);
  R2 tangente, normale, vec_u;
  
  // liste des points pour chaque cote de raidisseur
  Vector<IVect> list_points_edge(mesh.GetNbBoundaryRef());
  Vector<VectR2> list_points_coor(mesh.GetNbBoundaryRef());
  Vector<R2> list_normale(mesh.GetNbBoundaryRef());
  for (int p = 0; p < 2; p++)
    {
      nb_intersection = 0;
      int first_ext = nmpt; length_spline(nb_intersection) = 0;
      for (int i = 0; i < mesh_curve.GetNbBoundaryRef(); i++)
        {      
          R2 ptA = mesh_curve.Vertex(mesh_curve.BoundaryRef(i).numVertex(0));
          R2 ptB = mesh_curve.Vertex(mesh_curve.BoundaryRef(i).numVertex(1));
          VectR2 PtsMo(order+1), PtsM(order+1), NormaleM(order+1);
          PtsMo(0) = ptA;
          for (int j = 0; j < order-1; j++)
            PtsMo(j+1) = mesh_curve.GetPointInsideEdge(i, j);
          
          PtsMo(order) = ptB;
          
          // on bouge les points suivant la normale
          Real_wp delta = thickness;
          if (p == 0)
            delta = -thickness;
          
          for (int j = 0; j <= order; j++)
            {
              if (i >0 && (j ==0 ))
                {
                }
              else
                {
                  tangente.Fill(0); normale.Fill(0);
                  for (int k = 0; k <= order; k++)
                    Add(lob.GradPhi(k, j), PtsMo(k), tangente);
                  
                  normale(0) = -tangente(1); normale(1) = tangente(0);
                }
	      
	      Mlt(1.0/Norm2(normale), normale);
	      NormaleM(j) = 0.5*delta*normale;
              PtsM(j) = PtsMo(j) + 0.5*delta*normale;
            }
          
          //DISP(PtsM);
          
          ptA = PtsM(0); ptB = PtsM(order);
          
          // on rajoute le point initial du segment
          file_out << "Point(nmpt+"<<nmpt<<") = {" << PtsM(0)(0) << "*scalex, " << PtsM(0)(1) << "*scaley, 0, lc};\n";
          
          // pour  i = 0 on veut raccorder avec le point oppose
          if (i == 0)
            {
              liste_point_joint(nb_intersection)(p) = nmpt;
	      list_edge_joint(nb_intersection)(p) = -1;
	      list_spline_joint(nb_intersection)(p) = -1;
              nb_intersection++;
            }
          
          nmpt++;
          // point suivant, par defaut 1 (il peut changer si on trouve une intersection)
          int j0 = 1;
          
          // on boucle sur les aretes de reference ref du maillage mesh
          for (int k = 0; k < mesh.GetNbBoundaryRef(); k++)
            if (mesh.BoundaryRef(k).GetReference() == ref)
              {
                R2 pt1 = mesh.Vertex(mesh.BoundaryRef(k).numVertex(0));
                R2 pt2 = mesh.Vertex(mesh.BoundaryRef(k).numVertex(1));
                R2 point;
                
                // on teste si on trouve une intersection
                int success = IntersectionEdges(ptA, ptB, pt1, pt2, point, threshold);
                //DISP(ptA); DISP(ptB); DISP(pt1); DISP(pt2); DISP(point); DISP(success);
                
                // on a trouve une intersection
                if (success == 1)
                  {
                    //DISP(mesh_curve.BoundaryRef(i).numVertex(0));
                    //DISP(mesh_curve.BoundaryRef(i).numVertex(1));
                    //DISP(point); 
                    // on calcule l'intersection plus precisement avec Newton
                    Matrix2_2 df, inv_df;
                    Real_wp err = 1; R2 st, df_dt, dx;    
                    df(0, 1) = pt1(0) - pt2(0);
                    df(1, 1) = pt1(1) - pt2(1);
                    while (err > threshold)
                      {
                        Real_wp t = st(0), s = st(1);
                        df_dt.Fill(0); point.Fill(0);
                        for (int j = 0; j <= order; j++)
                          {
                            Add(lob.EvaluatePhi(j, t), PtsM(j), point);
                            Add(lob.EvaluatePhiGrad(j, t), PtsM(j), df_dt);
                          }
                        
                        point -= (1.0-s)*pt1 + s*pt2;
                        
                        // on calcule dx = DF^-1 F(x)
                        df(0, 0) = df_dt(0);
                        df(1, 0) = df_dt(1);
                        GetInverse(df, inv_df);
                        Mlt(inv_df, point, dx);
                        
                        // itere suivant
                        st(0) -= dx(0);
                        st(1) -= dx(1);
                        
                        err = Norm2(point);
                      }
                      
                    Real_wp t = st(0);
                    point.Fill(0);
                    for (int j = 0; j <= order; j++)
                      Add(lob.EvaluatePhi(j, t), PtsM(j), point);
                    
		    //Real_wp s = st(1);
                    //R2 ptC = (1.0-s)*pt1 + s*pt2;
                    
                    // on repere a quel j faut s'arreter
                    int j1 = 0;
                    for (int j = 0; j <= order; j++)
                      if (lob.Points(j) >= t)
                        {
                          j1 = j;
                          break;
                        }
                    
                    // on ecrit les points entre j0 et j1
                    for (int j = j0; j < j1; j++)
                      {                    
                        file_out << "Point(nmpt+"<<nmpt<<") = {" << PtsM(j)(0) << "*scalex, " << PtsM(j)(1) << "*scaley, 0, lc};\n";
			length_spline(nb_intersection) += PtsM(j).Distance(PtsM(j-1));
                        nmpt++;
                      }
                    
                    j0 = j1;
                    
		    if (!add_mesh_rib)
		      cout << "Numero du point intersection : " << nmpt << endl;
		    
                    // on cree le point intersection
                    //DISP(nmpt); DISP(point);
                    file_out << "Point(nmpt+"<<nmpt << ") = {" << point(0) << "*scalex, " << point(1) << "*scaley, 0, lc};\n"<<endl;
		    length_spline(nb_intersection) += point.Distance(PtsM(j0-1));
                    liste_point_joint(nb_intersection)(p) = nmpt;
		    list_edge_joint(nb_intersection)(p) = k;
		    list_spline_joint(nb_intersection)(p) = nmline;
                    nb_intersection++;
                    
                    // on cree une spline entre first_ext et ce point
                    file_out << "Spline(nmline+"<<nmline<<") = {nmpt+"<<first_ext<<":nmpt+"<<nmpt<<"};\n"<<endl; 
                    //file_out << "Line(nmline+"<<nmline<<") = {nmpt+"<<first_ext<<",nmpt+"<<nmpt<<"};\n"<<endl; 
		    
		    // on stocke le point intersection dans list_points_edge
		    vec_u = point - pt1;
		    int num_pt_add(-1);
		    if (DotProd(vec_u, NormaleM(j1)) > 0)
		      num_pt_add = mesh.BoundaryRef(k).numVertex(0);
		    else
		      num_pt_add = mesh.BoundaryRef(k).numVertex(1);
		    
		    //DISP(k); DISP(mesh.BoundaryRef(k)); DISP(num_pt_add);
		    //DISP(p); DISP(nb_intersection); DISP(pt1); DISP(pt2); DISP(num_pt_add); DISP(mesh.Vertex(num_pt_add));
		    //DISP(list_points_edge(k));
		    if (list_points_edge(k).GetM() == 0)
		      {
			list_points_edge(k).Reallocate(2); list_points_coor(k).Reallocate(2);
			list_normale(k) = NormaleM(j1);			
			list_points_edge(k)(0) = num_pt_add;
			list_points_edge(k)(1) = -nmpt;
			list_points_coor(k)(0) = mesh.Vertex(list_points_edge(k)(0));
			list_points_coor(k)(1) = point;
		      }
		    else
		      {
			AddPoint(list_points_edge(k), list_normale(k), list_points_coor(k),
				 num_pt_add, mesh.Vertex(num_pt_add));
			
			AddPoint(list_points_edge(k), list_normale(k), list_points_coor(k),
				 -nmpt, point);
		      }

		    //DISP(list_points_edge(k));
		    //DISP(list_points_coor(k));
		    length_spline(nb_intersection) += point.Distance(PtsM(j0));		    
                    first_ext = nmpt;
                    nmpt++; nmline++;
                  }  
              }
          
          // fin des points
          for (int j = j0; j < order; j++)
            {                    
              file_out << "Point(nmpt+"<<nmpt<<") = {" << PtsM(j)(0) << "*scalex, " << PtsM(j)(1) << "*scaley, 0, lc};\n";
	      length_spline(nb_intersection) += PtsM(j+1).Distance(PtsM(j));
              nmpt++;
            }          

          // on rajoute le point final
          if (i == mesh_curve.GetNbBoundaryRef()-1)
            {
              //DISP(nmpt); DISP(PtsM(order));
              file_out << "Point(nmpt+"<<nmpt<<") = {" << PtsM(order)(0) << "*scalex, " << PtsM(order)(1) << "*scaley, 0, lc};\n";
              
              liste_point_joint(nb_intersection)(p) = nmpt;
	      list_edge_joint(nb_intersection)(p) = -1;
	      list_spline_joint(nb_intersection)(p) = nmline;
              nb_intersection++;
              
              // spline de fin
              file_out << "Spline(nmline+"<<nmline<<") = {nmpt+"<<first_ext<<":nmpt+"<<nmpt<<"};\n"<<endl; 
              //file_out << "Line(nmline+"<<nmline<<") = {nmpt+"<<first_ext<<",nmpt+"<<nmpt<<"};\n"<<endl; 
              first_ext = nmpt; nmline++;               nmpt++;
            }
        }            
    }
  //DISP(length_spline);
    
  if (add_mesh_rib)
    {
      file_out << endl;
      
      // tableaux pour trouver l'intersection a partir du numero de point
      IVect point_to_intersec(nmpt+mesh.GetNbVertices());
      IVect point_to_loc(nmpt+mesh.GetNbVertices());
      point_to_intersec.Fill(-1); point_to_loc.Fill(-1);
      for (int i = 0; i < nb_intersection; i++)
	{
	  int n1 = liste_point_joint(i)(0);
	  int n2 = liste_point_joint(i)(1);
	  point_to_intersec(n1) = i; point_to_loc(n1) = 0;
	  point_to_intersec(n2) = i; point_to_loc(n2) = 1;
	}
      
      int offset_point = nmpt;
      // on rajoute tous les points des raidisseurs
      Vector<bool> VertexUsed(mesh.GetNbVertices());
      VertexUsed.Fill(false);
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	if (mesh.BoundaryRef(i).GetReference() == ref)
	  {
	    int n1 = mesh.BoundaryRef(i).numVertex(0);
	    int n2 = mesh.BoundaryRef(i).numVertex(1);
	    VertexUsed(n1) = true;
	    VertexUsed(n2) = true;
	  }
      
      for (int i = 0; i < mesh.GetNbVertices(); i++)
	if (VertexUsed(i))
	  file_out << "Point(nmpt+" << offset_point+i << ") = {" << mesh.Vertex(i)(0) << "*scalex, " << mesh.Vertex(i)(1) << "*scaley, 0, lc};\n";
      
      file_out << endl;
      IVect offset_points_edge(mesh.GetNbBoundaryRef());
      offset_points_edge.Fill(-1);
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	{
	  offset_points_edge(i) = nmline;
	  for (int k = 0; k < list_points_edge(i).GetM(); k++)
	    {
	      int num = list_points_edge(i)(k);
	      if (k > 0)
		{
		  // on rajoute l'arete si besoin entre deux points consecutifs 
		  int num_prev = list_points_edge(i)(k-1);
		  int n1 = -num_prev, n2 = -num;
		  if (num_prev >= 0)
		    n1 = offset_point+num_prev;
		  
		  if (num >= 0)
		    n2 = offset_point+num;
		  
		  file_out << "Line(nmline+" << nmline << ") = {nmpt+" << n1 << ", nmpt+" << n2 << "};\n";
		  nmline++;
		}	    
	    }
	}
      
      // on stocke les numeros des aretes joignant les points
      IVect list_lateral_joint(nb_intersection);
      list_lateral_joint.Fill(-1);
      list_lateral_joint(0) = nmline;
      list_lateral_joint(nb_intersection-1) = nmline+1;
      for (int i = 1; i < nb_intersection-1; i++)
	{
	  int ne = list_edge_joint(i)(0);
	  //DISP(i); DISP(ne); DISP(offset_points_edge(ne)); DISP(list_points_edge(ne));
	  //DISP(list_points_coor(ne));
	  int n1 = liste_point_joint(i)(0);
	  int n2 = liste_point_joint(i)(1);
	  for (int k = 0; k < list_points_edge(ne).GetM()-1; k++)
	    {
	      if ((list_points_edge(ne)(k) == -n1) && (list_points_edge(ne)(k+1) == -n2))
		list_lateral_joint(i) = offset_points_edge(ne) + k;
	      
	      if ((list_points_edge(ne)(k) == -n2) && (list_points_edge(ne)(k+1) == -n1))
		list_lateral_joint(i) = -(offset_points_edge(ne) + k);
	    }
	}
      
      // on met des lignes pour joindre les points finaux
      list<int> edge_bridge, edge_rib;
      edge_bridge.push_back(nmline); edge_bridge.push_back(nmline+1);
      file_out << endl;
      for (int i = 0; i < nb_intersection; i+= nb_intersection-1)
	{
	  file_out << "Line(nmline+" << nmline << ") = {nmpt+"<<liste_point_joint(i)(0) <<", nmpt+" << liste_point_joint(i)(1)<<"};\n";
	  nmline++;
	}      
      
      // on construit les lignes des raidisseurs qui n'intersectent pas le chevalet
      file_out << endl;
      IVect NumLine(mesh.GetNbBoundaryRef()); NumLine.Fill(-1);
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	if (mesh.BoundaryRef(i).GetReference() == ref)
	  {
	    if (list_points_edge(i).GetM() == 0)
	      {
		int n1 = mesh.BoundaryRef(i).numVertex(0);
		int n2 = mesh.BoundaryRef(i).numVertex(1);
		NumLine(i) = nmline;
		file_out << "Line(nmline+" << nmline << ") = {nmpt+" << offset_point+n1 <<
		  ", nmpt+" << offset_point+n2 << "};\n";
		nmline++;
	      }
	  }
      
      // les lignes "laterales" sont definies transfinite en prenant l'epaisseur
      file_out << endl;
      file_out << "Transfinite Line{nmline+" << list_lateral_joint(0);
      for (int i = 1; i < nb_intersection; i++)
	file_out << ", nmline+" << abs(list_lateral_joint(i));
      
      file_out << "} = Ceil(" << thickness << "/lc)+1;" << endl << endl;
      
      // on construit les surfaces
      list<int> surface_bridge, surface_rib, surface_intersec;
      for (int i = 1; i < nb_intersection; i++)
	{
	  int ne = list_edge_joint(i)(0);
	  if (ne != list_edge_joint(i)(1))
	    {
	      cout << "Cas interdit, les intersections doivent appartenir a la meme arete" << endl;
	      abort();
	    }
	  
	  int ne_prev = list_edge_joint(i-1)(0);
	  int p1 = -1, p2 = -1;
	  if ((ne_prev >= 0) && (ne >= 0))
	    {
	      int n1 = list_points_edge(ne_prev)(0);
	      int n2 = list_points_edge(ne)(0);
	      p1 = mesh.FindEdgeRefWithExtremities(n1, n2);

	      int N1 = list_points_edge(ne_prev).GetM();
	      int N2 = list_points_edge(ne).GetM();
	      n1 = list_points_edge(ne_prev)(N1-1);
	      n2 = list_points_edge(ne)(N2-1);
	      p2 = mesh.FindEdgeRefWithExtremities(n1, n2);
	    }
	  
	  if ((p1 >= 0) && (p2 >= 0))
	    {
	      // raidisseur
	      int ne1, ne2, ne3, ne4;
	      ne4 = NumLine(p1);
	      if (mesh.BoundaryRef(p1).numVertex(0) != list_points_edge(ne_prev)(0))
		ne4 = -ne4;
	      
	      // boucle sur les aretes
	      int na = -1;
	      string transf_line = "Transfinite Line{nmline+" + to_str(abs(ne4));
	      for (int k = 0; k < list_points_edge(ne).GetM()-1; k++)
		{
		  // numero des 4 sommets du quadrilatere a rajouter
		  int nv1 = list_points_edge(ne_prev)(k);
		  int nv2 = list_points_edge(ne_prev)(k+1);
		  int nv3 = list_points_edge(ne)(k+1);
		  int nv4 = list_points_edge(ne)(k);

		  if (k == list_points_edge(ne).GetM()-2)
		    {
		      na = NumLine(p2);
		      if (mesh.BoundaryRef(p2).numVertex(0) != nv2)
			na = -na;		      
		    }
		  else
		    {
		      // on trouve l'intersection	  
		      int ni1 = point_to_intersec(-nv2);
		      int ni2 = point_to_intersec(-nv3);
		      int loc1 = point_to_loc(-nv2);
		      //int loc2 = point_to_loc(nv3);
		      
		      // numero de l'arete face a ne4		      
		      if (ni2 == ni1+1)
			na = list_spline_joint(ni2)(loc1);
		      else if (ni1 == ni2+1)
			na = -list_spline_joint(ni1)(loc1);
		    }
		  
		  nv1 = ChangeNum(nv1, offset_point); 
		  nv2 = ChangeNum(nv2, offset_point); 
		  nv3 = ChangeNum(nv3, offset_point); 
		  nv4 = ChangeNum(nv4, offset_point); 
		  
		  // numero des 4 aretes du quadrilatere a rajouter
		  ne1 = offset_points_edge(ne_prev) + k;
		  ne2 = abs(na);
		  ne3 = offset_points_edge(ne) + k;
		  //DISP(nv1); DISP(nv2); DISP(nv3); DISP(nv4);
		  //DISP(ne1); DISP(ne2); DISP(ne3); DISP(abs(ne4)); DISP(na);
		  
		  string ne1s = "nmline+" + to_str(ne1);
		  string ne3s = "-(nmline+" + to_str(ne3) + ")";
		  string ne2s = "nmline+" + to_str(ne2);
		  if (na < 0)
		    ne2s = "-(nmline+" + to_str(ne2) + ")";
		  
		  string ne4s = "nmline+" + to_str(abs(ne4));
		  if (ne4 > 0)
		    ne4s = "-(nmline+" + to_str(ne4) + ")";
		  
		  file_out << '\n';
		  
		  Real_wp L = list_points_coor(ne)(k).Distance(list_points_coor(ne)(k+1));
		  transf_line += ", nmline+" + to_str(abs(na));
		  if (k%2 == 0)
		    file_out << "Transfinite Line{nmline+" << abs(ne1) << ", nmline+" << abs(ne3)
			     << "} = Ceil(" << L << "/lc)+1;\n";
		  
		  file_out << "Line Loop(nmsurf+" << nmsurf << ") = {" << ne1s << ", " << ne2s 
			   << ", " << ne3s << ", " << ne4s << "};\n";
		  
		  file_out << "Plane Surface(nmsurf+" << nmsurf << ") = {nmsurf+" << nmsurf << "};\n";
		  file_out << "Transfinite Surface{nmsurf+" << nmsurf << "} = {nmpt+" << nv1 
			   << ", nmpt+" << nv2 << ", nmpt+" << nv3 << ", nmpt+" << nv4 << "};\n";
		  		  
		  if (k%2 == 0)
		    surface_rib.push_back(nmsurf);
		  else
		    surface_intersec.push_back(nmsurf);
		  
		  edge_rib.push_back(abs(ne1));
		  edge_rib.push_back(abs(ne3));
		  if (k == 0)
		    edge_rib.push_back(abs(ne4));
		  else
		    edge_bridge.push_back(abs(ne4));
		  
		  ne4 = na;
		  nmsurf++;
		}
	      
	      edge_rib.push_back(abs(na));
	      
	      Real_wp L = list_points_coor(ne)(0).Distance(list_points_coor(ne_prev)(0));
	      file_out << transf_line << "} = Ceil(" << L << "/lc)+1;\n";
	    }
	  else
	    {
	      // chevalet	      
	      // numero des 4 aretes du "quadrilatere" a rajouter
	      int ne1 = list_spline_joint(i)(0);
	      int ne3 = list_spline_joint(i)(1);
	      int ne2 = list_lateral_joint(i);
	      int ne4 = list_lateral_joint(i-1);
	      //DISP(ne1); DISP(ne2); DISP(ne3); DISP(ne4);
	      
	      // numero des 4 sommets du "quadrilatere" a rajouter
	      int nv1 = liste_point_joint(i-1)(0);
	      int nv2 = liste_point_joint(i-1)(1);
	      int nv3 = liste_point_joint(i)(1);
	      int nv4 = liste_point_joint(i)(0);

	      // chaines de caractere pour les 4 aretes
	      string ne1s = "nmline+" + to_str(ne1);
	      string ne3s = "-(nmline+"+to_str(ne3) + ")";
	      string ne2s = "nmline+" + to_str(ne2);
	      string ne4s = "-(nmline+" + to_str(ne4) + ")";
	      if (ne2 < 0)
		ne2s = "-(nmline+" + to_str(-ne2) + ")";
	      
	      if (ne4 < 0)
		ne4s = "nmline+" + to_str(-ne4);
	      	      
	      // on rajoute le quadrilatere dans le fichier geo avec des transfinite
	      file_out << '\n';
	      file_out << "Transfinite Line{nmline+" << ne1s << ", nmline+" << ne3 << "} = Ceil(" 
		       << length_spline(i) << "/lc)+1;\n"; 
	      
	      file_out << "Line Loop(nmsurf+" << nmsurf << ") = {" << ne1s << ", " 
		       << ne2s << ", " << ne3s << ", " << ne4s << "};\n";
	      
	      file_out << "Plane Surface(nmsurf+" << nmsurf << ") = {nmsurf+" << nmsurf << "};\n";
	      file_out << "Transfinite Surface{nmsurf+" << nmsurf << "} = {nmpt+" << nv1 
		       << ", nmpt+" << nv2 << ", nmpt+" << nv3 << ", nmpt+" << nv4 << "};\n";
	      
	      file_out << '\n';
	      surface_bridge.push_back(nmsurf);
	      edge_bridge.push_back(abs(ne1));
	      edge_bridge.push_back(abs(ne3));
	      nmsurf++;
	    }
	}
      
      file_out << endl;
      
      // aretes du chevalet
      IVect num_edge;
      Copy(edge_bridge, num_edge);
      file_out << "Physical Line(1) = {nmline+" << num_edge(0);
      for (int i = 1; i < num_edge.GetM(); i++)
	file_out << ", nmline+" << num_edge(i);
      
      file_out << "};" << endl;
      
      // aretes des raidisseurs
      Copy(edge_rib, num_edge);
      file_out << "Physical Line(2) = {nmline+" << num_edge(0);
      for (int i = 1; i < num_edge.GetM(); i++)
	file_out << ", nmline+" << num_edge(i);
      
      file_out << "};" << endl;

      // physical surface pour le chevalet uniquement
      IVect num_surf;
      Copy(surface_bridge, num_surf);
      file_out << "Physical Surface(1) = {nmsurf+" << num_surf(0);
      for (int i = 1; i < num_surf.GetM(); i++)
	file_out << ", nmsurf+" << num_surf(i);
      
      file_out << "};" << endl;

      // physical surface pour les raidisseurs uniquement
      Copy(surface_rib, num_surf);
      file_out << "Physical Surface(2) = {nmsurf+" << num_surf(0);
      for (int i = 1; i < num_surf.GetM(); i++)
	file_out << ", nmsurf+" << num_surf(i);
      
      file_out << "};" << endl;

      // physical surface pour les intersections
      Copy(surface_intersec, num_surf);
      file_out << "Physical Surface(3) = {nmsurf+" << num_surf(0);
      for (int i = 1; i < num_surf.GetM(); i++)
	file_out << ", nmsurf+" << num_surf(i);
      
      file_out << "};" << endl;      
    }
  else
    {
      // on met des lignes pour joindre les intersections
      file_out << endl;
      for (int i = 0; i < nb_intersection; i++)
	{
	  file_out << "Line(nmline+" << nmline << ") = {nmpt+"<<liste_point_joint(i)(0) <<", nmpt+" << liste_point_joint(i)(1)<<"};\n";
	  nmline++;
	}
    }

  file_out.close();


  return FinalizeMontjoie();
}
