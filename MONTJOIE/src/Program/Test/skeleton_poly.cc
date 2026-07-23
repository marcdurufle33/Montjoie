#include "Montjoie.hxx"

using namespace Montjoie;

class InputData : public InputDataProblem_Base
{
public :
  VectR2 poly_line;
  
  void SetInputData(const string& keyword, const Vector<string>& parameters)
  {    
    if (keyword == "PolyLine")
      {
        poly_line.Reallocate(parameters.GetM());
        for (int k = 0; k < poly_line.GetM(); k++)
          {
            Complex_wp z = to_num<Complex_wp>(parameters(k));
            poly_line(k).Init(real(z), imag(z));
          }
      }
  }
  
};

inline double isLeft(double P0x, double P0y, double P1x, double P1y, double P2x, double P2y)
{
  return ( (P1x - P0x) * (P2y - P0y)
           - (P2x -  P0x) * (P1y - P0y) );
}

bool InsidePolygon(const R2& pt, const VectR2& polygon)
{
  int wn = 0;    // the  winding number counter

  int n = polygon.GetM();
  
  // loop through all edges of the polygon
  for (int i = 0; i < n-1; i++)
    {
      // edge from V[i] to  V[i+1]
      if (polygon(i)(1) <= pt(1))
        {
          // start y <= P.y
          if (polygon(i+1)(1)  > pt(1))      // an upward crossing
            if (isLeft( polygon(i)(0), polygon(i)(1), polygon(i+1)(0), polygon(i+1)(1), pt(0), pt(1)) > 0)  // P left of  edge
              ++wn;            // have  a valid up intersect
        }
      else
        {
          // start y > P.y (no test needed)
          if (polygon(i+1)(1)  <= pt(1))     // a downward crossing
            if (isLeft(polygon(i)(0), polygon(i)(1), polygon(i+1)(0), polygon(i+1)(1), pt(0), pt(1)) < 0)
              --wn;            // have  a valid down intersect
        }
    }
  
  if (wn != 0)
    return true;
  
  return false;
}

int FindSolutionPol2(const Real_wp& a, const Real_wp& b, const Real_wp& c, Real_wp& sol1, Real_wp& sol2)
{
  if (abs(a) <= 1e-12)
    {
      if (abs(b) <= 1e-12)
        {
          if (abs(c) <= 1e-12)
            return 3;
          else
            return 0;
        }
      else
        {
          sol1 = -c / (2*b);
          return 1;
        }
    }
  else
    {
      Real_wp d = b*b - a*c;
      if (abs(d) <= 1e-12)
        {
          sol1 = -b/a;
          return 1;
        }
      else if (d < 0)
        return 0;
      else
        {
          d = sqrt(d);
          sol1 = (-b + d) / a;
          sol2 = (-b - d) / a;
          return 2;
        }
    }
}

class EdgeVoronoi
{
public:
  R2 point, v;
  int num_i, num_j;
  bool pt_i, pt_j;
  R2 ext;
  bool linear;

public :
  void ProjectToBisector(const VectR2& polygon, const VectR2& normale, R2& M)
  {
    R2 pt;
    if (pt_i)
      pt = polygon(num_i);

    if (pt_j)
      pt = polygon(num_j);

    R2 n, pt0;
    if (pt_i)
      {
        n = normale(num_j);
        pt0 = polygon(num_j);
      }

    if (pt_j)
      {
        n = normale(num_i);
        pt0 = polygon(num_i);
      }

    R2 A = pt;
    R2 v = M - pt;
    v *= 1.0/Norm2(v);
    
    R2 vec_u0 = A - pt0;
    Real_wp scal0 = DotProd(n, vec_u0);
    Real_wp vn = DotProd(n, v);

    // parametres parabole de la droite
    Real_wp a = vn*vn, b = vn*scal0, c = scal0*scal0;

    // on retranche la parabole du point
    a -= 1.0;

    Real_wp sol1, sol2;
    FindSolutionPol2(a, b, c, sol1, sol2);

    Real_wp x = sol2;
    if (sol1 > 0)
      x = sol1;

    M = pt + x*v;
  }
  
};

class ParaboleVoronoi
{
public:
  Real_wp a, b, c;
  Real_wp tau1, tau2;
  int num; bool edge;

public:

  Real_wp FindIntervalPositive(const Real_wp& a, const Real_wp& b, const Real_wp& c)
  {
    if (abs(a) <= 1e-12)
      {
        if (abs(b) <= 1e-12)
          {
            if (c >= 0)
              return Real_wp(1e30);
            else
              {
                if (tau1 > 0)
                  return tau1;
              }
          }
        else
          {
            Real_wp sol1 = -c / (2*b);
            if (b > 0)
              {
                if ((tau1 > 0) && (tau1 < sol1))
                  return tau1;
              }
            else
              {
                if (sol1 < 0)
                  {
                    if (tau1 > 0)
                      return tau1;
                  }
                else
                  {
                    if (tau1 > sol1)
                      return tau1;
                    else
                      {
                        if (tau2 > sol1)
                          return sol1;
                      }
                  }
              }
          }
      }
    else
      {
        Real_wp d = b*b - a*c;
        if (d <= 1e-12)
          {
            if (a < 0)
              {
                if (tau1 > 0)
                  return tau1;
              }
          }
        else
          {
            d = sqrt(d);
            Real_wp sol1 = (-b + d) / a;
            Real_wp sol2 = (-b - d) / a;
            if (sol1 > sol2)
              {
                Real_wp sol0 = sol1;
                sol1 = sol2;
                sol2 = sol0;
              }
            
            if (a > 0)
              {
                Real_wp ext1 = tau1, ext2 = tau2;
                if (sol1 > tau1)
                  {
                    if (sol1 < tau2)
                      ext1 = sol1;
                    else
                      return Real_wp(1e30);
                  }

                if (sol2 < tau2)
                  {
                    if (sol2 > tau1)
                      ext2 = sol2;
                    else
                      return Real_wp(1e30);
                  }

                if (ext1 > 0)
                  return ext1;
              }
            else
              {
                if (sol2 > tau1)
                  {
                    if (sol1 > tau1)
                      {
                        if (sol2 >= tau2)
                          {
                            if (tau1 > 0)
                              return tau1;
                          }
                        else
                          {
                            if (sol1*sol2 <= 0)
                              return sol2;
                            else
                              {
                                if (sol1 > 0)
                                  return tau1;
                              }
                          }
                      }
                    else
                      {
                        if (sol2 < tau2)
                          {
                            if (sol2 > 0)
                              return sol2;
                          }
                      }
                  }
                else
                  {
                    if (tau1 > 0)
                      return tau1;
                  }
              }
          }
      }
    
    return Real_wp(1e30);
  }

};

ostream& operator<<(ostream& out, const ParaboleVoronoi& p)
{
  out << "a, b, c: " << p.a << " " << p.b << " " << p.c << " num : " << p.num
      << " edge: " << p.edge << " tau1 " << p.tau1 <<" tau2 " << p.tau2;
  
  return out;                                                     
}

ostream& operator<<(ostream& out, const EdgeVoronoi& e)
{
  out << "Pt: " << e.point << " , Axis: " << e.v << " , num: " << e.num_i << " " << e.num_j
      << ", is_point " << e.pt_i << " " << e.pt_j <<", Linear : " << e.linear;
  
  return out;                                                     
}

template<class T>
ostream& operator<<(ostream& out, const vector<T>& v)
{
  for (int i = 0; i < v.size(); i++)
    out << v[i] << '\n';
  
  return out;
}


EdgeVoronoi CreateEdgeVoronoi(const VectR2& polygon, const Real_wp& teta1, const Real_wp& teta2,
                              int i, bool ie, int j, bool je, const R2& ptA)
{
  EdgeVoronoi e;
  e.point = ptA;

  // angle between teta1 and teta2
  Real_wp teta = 0.5*(teta1+teta2);
  if (teta1 > teta2)
    teta += pi_wp;

  // bisector
  e.v.Init(cos(teta), sin(teta));
  e.num_i = i; e.num_j = j;
  e.pt_i = !ie; e.pt_j = !je;      
  e.linear = true;
  int N = polygon.GetM()-1;
  if (ie && !je)
    {
      e.linear = false;
      if ((j == i) || (j == (i+1)%N))
        e.linear = true;
    }

  if (je && !ie)
    {
      e.linear = false;
      if ((i == j) || (i == (j+1)%N))
        e.linear = true;
    }
  
  return e;
}

void FillDistancePolygon(const VectR2& polygon, const VectR2& normale,
                         const R2& M, VectReal_wp& dist)
{
  dist.Reallocate(2*(polygon.GetM()-1));
  dist.Fill(1e30);
  for (int k = 0; k < polygon.GetM()-1; k++)
    {
      R2 x0 = polygon(k), x1 = polygon(k+1);
      
      // distance avec les aretes
      Real_wp alpha = DotProd(M-x0, normale(k));
      R2 ptP = M - alpha*normale(k);
      Real_wp longueur = AbsSquare(x1-x0);
      if ((AbsSquare(ptP-x0) <= longueur+1e-12)
          && (AbsSquare(ptP-x1) <= longueur+1e-12))
        dist(k) = abs(alpha);
              
      alpha = x0.Distance(M);
      dist(polygon.GetM()-1 + k) = alpha;
    }
}

vector<R2> GetVoronoiVertices(const VectR2& polygon)
{
  Real_wp diam = 0;
  for (int i = 0; i < polygon.GetM(); i++)
    for (int j = 0; j < polygon.GetM(); j++)
      if (i != j)
        diam += AbsSquare(polygon(i)-polygon(j));

  diam = sqrt(diam);
  
  VectR2 normale;

  vector<EdgeVoronoi> edge;
  
  // on initialise les aretes de Voronoi
  int N = polygon.GetM()-1;
  normale.Reallocate(N);
  for (int i = 0; i < polygon.GetM()-1; i++)
    {
      EdgeVoronoi e;
      R2 ptA, ptB(polygon(i)), ptC(polygon(i+1));
      int nj = i-1;
      if (i == 0)
        {
          ptA = polygon(N-1);
          nj = N-1;
        }
      else
        ptA = polygon(i-1);
      
      e.point = polygon(i);
      
      R2 vec_u = ptA - ptB; vec_u *= 1.0/Norm2(vec_u);
      R2 vec_v = ptC - ptB; vec_v *= 1.0/Norm2(vec_v);
      normale(i).Init(-vec_v(1), vec_v(0));
      
      R2 vec_w = 0.5*(vec_u + vec_v); vec_w *= 1.0/Norm2(vec_w);
      //DISP(i); DISP(ptA); DISP(ptB); DISP(ptC);
      //DISP(vec_u); DISP(vec_v); DISP(vec_w);
      R2 ptM = ptB + 1e-6*vec_w; //DISP(ptM);
      if (InsidePolygon(ptM, polygon))
        {
          //cout << "Inside polygon" << endl;
          e.v = vec_w;
          e.num_i = nj;
          e.num_j = i;
          e.pt_i = false;
          e.pt_j = false;
          e.linear = true;
          edge.push_back(e);
        }
      else
        {
          e.v.Init(vec_u(1), -vec_u(0));
          if (DotProd(e.v, vec_w) > 0)
            e.v = -e.v;
          
          e.num_i = nj;
          e.num_j = i;
          e.pt_i = false;
          e.pt_j = true;
          e.linear = true;
          edge.push_back(e);

          e.v.Init(-vec_v(1), vec_v(0));
          if (DotProd(e.v, vec_w) > 0)
            e.v = -e.v;
                    
          e.num_i = i;
          e.num_j = i;
          e.pt_i = false;
          e.pt_j = true;
          edge.push_back(e);
        }
    }

  DISP(normale);
  //
  cout << "Liste des aretes = " << edge << endl;
  
  // tant qu'on a des aretes pas traitees, on continue
  vector<EdgeVoronoi> all_edge;
  vector<R2> all_vertex;
  bool use_bisection_algo = false;
  while (edge.size() > 0)
    {
      EdgeVoronoi e = edge[edge.size()-1];
      R2 A = e.point;
      R2 v = e.v;
      cout << endl;
      cout << "Arete traitee " << endl;
      DISP(e);
      
      // on regarde si l'extremite de l'arete n'est pas presente dans edge
      bool stop_loop = false;
      for (int i = 0; i < edge.size()-1; i++)
        {
          bool edge_found = false;
          if ((e.num_i == edge[i].num_i) && (e.pt_i == edge[i].pt_i)
              && (e.num_j == edge[i].num_j) && (e.pt_j == edge[i].pt_j))
            edge_found = true;

          if ((e.num_i == edge[i].num_j) && (e.pt_i == edge[i].pt_j)
              && (e.num_j == edge[i].num_i) && (e.pt_j == edge[i].pt_i))
            edge_found = true;

          if (edge_found)
            {
              cout << "Arete trouvee " << i << endl;
              stop_loop = true;
              e.ext = edge[i].point;
              edge.erase(edge.begin()+i);
              all_edge.push_back(e);
              edge.pop_back();
              break;
            }
        }

      if (stop_loop)
        continue;
      
      // on calcule la parabole pour chaque arete
      vector<ParaboleVoronoi> liste_parabole;
      if (e.linear)
        {
          for (int i = 0; i < N; i++)
            {
              ParaboleVoronoi p;
              bool add_parabole = true;
              p.num = i; p.edge = true;
              
              R2 pt0(polygon(i)), pt1(polygon(i+1));
              Real_wp tau0(0), tau1(0), tau2(0);
              
              // on cherche les tau pour lesquels la distance a la droite est egale a la distance
              // au pt0 et pt1
              R2 vec_u0 = A - pt0;
              Real_wp scal0 = DotProd(normale(i), vec_u0);
              Real_wp vn = DotProd(normale(i), v);
              Real_wp dist0 = AbsSquare(vec_u0);
              
              R2 vec_u1 = A - pt1;
              Real_wp scal1 = DotProd(normale(i), vec_u1);
              Real_wp dist1 = AbsSquare(vec_u1);
              
              // parametres de la parabole
              p.a = vn*vn; p.b = vn*scal0; p.c = scal0*scal0;          
              Real_wp a = 1.0 - vn*vn;
              if (abs(a) <= 1e-15)
                {
                  Real_wp c0 = dist0 - scal0*scal0;
                  Real_wp c1 = dist1 - scal1*scal1;
                  Real_wp longueur = AbsSquare(pt0 - pt1);
                  if ((c0 <= longueur+1e-12) && (c1 <= longueur+1e-12))
                    {
                      // cas ou la droite passe entre les deux points en etant perpendiculaire
                      p.tau1 = -1e30; p.tau2 = 1e30;
                    }
                  else
                    add_parabole = false; 
                }
              else
                {
                  // cas d'une droite oblique : on cherche les intersections
                  Real_wp b = DotProd(v, vec_u0) - vn*scal0; 
                  Real_wp c = dist0 - scal0*scal0;
                  
                  int nb_sol = FindSolutionPol2(a, b, c, tau0, tau2);
                  
                  
                  b = DotProd(v, vec_u1) - vn*scal1; 
                  c = dist1 - scal1*scal1;
                  
                  nb_sol = FindSolutionPol2(a, b, c, tau1, tau2);
                  if (tau0 > tau1)
                    {
                      tau2 = tau0;
                      tau0 = tau1;
                      tau1 = tau2;
                    }
                  
                  p.tau1 = tau0 - 1e-14; p.tau2 = tau1+1e-14;
                }
              
              if (add_parabole)
                {
                  liste_parabole.push_back(p);
                  //DISP(p.tau1); DISP(p.tau2);
                }
            }

          // on calcule la parabole pour chaque point
          for (int i = 0; i < N; i++)
            {
              ParaboleVoronoi p;
              p.tau1 = -1e30;  p.tau2 = 1e30;
              p.num = i; p.edge = false;
              
              R2 pt0(polygon(i));
              R2 vec_u0 = A - pt0;
              
              p.a = 1.0; p.b = DotProd(vec_u0, v);
              p.c = AbsSquare(vec_u0);
              
              liste_parabole.push_back(p);
            }
        }
      
      //DISP(liste_parabole);
      Vector<Real_wp> all_dist(liste_parabole.size());

      // dichotomie pour trouver le point final
      if (use_bisection_algo || !e.linear)
        {
          Real_wp a = 0, b = diam;
          Real_wp coef = Real_wp(1) + 1e-14;
          while (abs(a-b) > 1e-12*diam)
            {
              Real_wp tau = (a+b)/2;
              
              // on evalue toutes les distances
              Real_wp disti = 0, distj = 0;
              bool outside_edge = false;
              if (e.linear)
                {
                  for (int i = 0; i < liste_parabole.size(); i++)
                    {
                      ParaboleVoronoi& p = liste_parabole[i];
                      if ((tau >= p.tau1) && (tau <= p.tau2))
                        all_dist(i) = p.c + tau*(2.0*p.b + tau*p.a);
                      else
                        all_dist(i) = 1e30;
                      
                      if ((p.edge != e.pt_i) && (p.num == e.num_i))
                        disti = all_dist(i);
                      
                      if ((p.edge != e.pt_j) && (p.num == e.num_j))
                        distj = all_dist(i);
                    }
                  
                  for (int i = 0; i < liste_parabole.size(); i++)
                    if (coef*all_dist(i) < disti)
                      outside_edge = true;
                }
              else
                {
                  R2 M = A + tau*v; DISP(M);
                  e.ProjectToBisector(polygon, normale, M);
                  DISP(M);

                  FillDistancePolygon(polygon, normale, M, all_dist);
                  DISP(all_dist);
                  int num_i = e.num_i;
                  if (e.pt_i)
                    num_i += polygon.GetM()-1;

                  int num_j = e.num_j;
                  if (e.pt_j)
                    num_j += polygon.GetM()-1;
                  
                  disti = all_dist(num_i); distj = all_dist(num_j);
                  DISP(num_i); DISP(num_j); DISP(disti); DISP(distj);
                  for (int i = 0; i < all_dist.GetM(); i++)
                    if ((i != num_i) && (i != num_j))
                      if (coef*all_dist(i) < disti)
                        outside_edge = true;

                  DISP(outside_edge);
                }
              
              if (outside_edge)
                b = tau;
              else
                a = tau;
              
            }
          
          //DISP(a); DISP(b);
          e.ext = A + a*v;
          if (!e.linear)
            {
              e.ProjectToBisector(polygon, normale, e.ext);
              FillDistancePolygon(polygon, normale, e.ext, all_dist);
              liste_parabole.resize(all_dist.GetM());
              for (int i = 0; i < liste_parabole.size(); i++)
                {
                  ParaboleVoronoi& p = liste_parabole[i];
                  if (i < polygon.GetM()-1)
                    {
                      p.num = i;
                      p.edge = true;
                    }
                  else
                    {
                      p.num = i-(polygon.GetM()-1);
                      p.edge = false;
                    }
                }
            }
        }
      else
        {
          int numi = -1, numj = -1;
          for (int i = 0; i < liste_parabole.size(); i++)
            {
              ParaboleVoronoi& p = liste_parabole[i];
              if ((p.edge != e.pt_i) && (p.num == e.num_i))
                numi = i;
              
              if ((p.edge != e.pt_j) && (p.num == e.num_j))
                numj = i;
            }
          
          //DISP(numi); DISP(numj);
          
          Real_wp a = liste_parabole[numi].a;
          Real_wp b = liste_parabole[numi].b;
          Real_wp c = liste_parabole[numi].c;
          Real_wp tau = min(liste_parabole[numi].tau2, liste_parabole[numj].tau2);
          
          for (int i = 0; i < liste_parabole.size(); i++)
            if ((i != numi) && (i != numj))
              {
                ParaboleVoronoi& p = liste_parabole[i];
                Real_wp tau1 = p.FindIntervalPositive(p.a - a, p.b - b, p.c - c);
                //DISP(i); DISP(tau1);
                tau = min(tau, tau1);
              }

          //DISP(tau);
          // on evalue toutes les distances
          for (int i = 0; i < liste_parabole.size(); i++)
            {
              ParaboleVoronoi& p = liste_parabole[i];
              if ((tau >= p.tau1) && (tau <= p.tau2))
                all_dist(i) = p.c + tau*(2.0*p.b + tau*p.a);
              else
                all_dist(i) = 1e30;
            }

          e.ext = A + tau*v; DISP(e.ext);
        }
      
      for (int i = 0; i < all_vertex.size(); i++)
        if (all_vertex[i] == e.ext)
          {
            DISP(i); DISP(all_vertex[i]);
            cout << "Impossible case " << endl;
            abort();
          }

      // on rajoute l'arete et le sommet cree
      all_edge.push_back(e);
      all_vertex.push_back(e.ext);
      edge.pop_back();      
      
      // on cherche toutes les aretes partant du point nouvellement cree
      Vector<int> permut(all_dist.GetM());
      permut.Fill();
      Sort(all_dist, permut);
      //DISP(numi); DISP(all_dist);
      int nb_sites = 0;
      while ((nb_sites < all_dist.GetM()) && (abs(all_dist(nb_sites) - all_dist(0)) <= 1e-8))
        nb_sites++;

      DISP(e.ext);
      DISP(nb_sites); permut.Resize(nb_sites);
      
      // directions for each site
      VectR2 dir(nb_sites);
      VectReal_wp theta(nb_sites);
      for (int k = 0; k < nb_sites; k++)
        {
          ParaboleVoronoi& p = liste_parabole[permut(k)];
          int i = p.num;
          if (p.edge)
            {
              R2 vec_u = e.ext - polygon(i);
              if (DotProd(normale(i), vec_u) > 0)
                dir(k) = -normale(i);
              else
                dir(k) = normale(i);
            }
          else
            {
              dir(k) = polygon(i) - e.ext;
              Mlt(1.0/Norm2(dir(k)), dir(k));
            }

          Real_wp r;
          CartesianToPolar(dir(k)(0), dir(k)(1), r, theta(k));
        }
      
      DISP(dir);

      // on trie les directions par l'angle
      Sort(theta, permut);

      DISP(theta); DISP(permut);
      for (int k = 0; k < nb_sites; k++)
        {
          ParaboleVoronoi& p = liste_parabole[permut(k)];
          ParaboleVoronoi& p2 = liste_parabole[permut((k+1)%nb_sites)];
          
          Real_wp teta1 = theta(k);
          Real_wp teta2 = theta((k+1)%nb_sites);
          
          int i = p.num; int j = p2.num;
          bool ie = p.edge; int je = p2.edge;

          bool new_edge = true;
          for (int q = 0; q < all_edge.size(); q++)
            {
              int i0 = all_edge[q].num_i;
              bool i0_e = !all_edge[q].pt_i;

              int j0 = all_edge[q].num_j;
              bool j0_e = !all_edge[q].pt_j;
              
              if ((i == i0) && (ie == i0_e) && (j == j0) && (je == j0_e))
                new_edge = false;

              if ((i == j0) && (ie == j0_e) && (j == i0) && (je == i0_e))
                new_edge = false;
            }

          if (new_edge)
            {
              EdgeVoronoi e2 = CreateEdgeVoronoi(polygon, teta1, teta2,
                                                 i, ie, j, je, e.ext);

              int num_i = e2.num_i;
              if (e2.pt_i)
                num_i += polygon.GetM()-1;

              int num_j = e2.num_j;
              if (e2.pt_j)
                num_j += polygon.GetM()-1;

              // on verifie que e2.v est bien oriente
              R2 M = e2.point + 1e-6*e2.v;
              FillDistancePolygon(polygon, normale, M, all_dist);

              Real_wp disti = all_dist(num_i), distj = all_dist(num_j);
              for (int k = 0; k < all_dist.GetM(); k++)
                if ((k != num_i) && (k != num_j))
                  {
                    if ((all_dist(k) < disti) || (all_dist(k) < distj))
                      {
                        e2.v = -e2.v;
                        break;
                      }
                  }
                                    
              if (!e2.linear)
                {
                  DISP(e2.point); DISP(e2.v);
                }
              
              // on insere apres les aretes paraboliques
              vector<EdgeVoronoi>::iterator it;
              for (it = edge.begin(); it != edge.end(); it++)
                {
                  if (it->linear)
                    break;
                }

              DISP(e2);
              edge.insert(it, e2);
            }
        }
    }

  ofstream file_out("aretes.dat");
  file_out.precision(15);
  for (int i = 0; i < all_edge.size(); i++)
    {
      EdgeVoronoi& e = all_edge[i];
      file_out << e.point(0) << " " << e.point(1) << " " << e.ext(0) << " " << e.ext(1) << '\n';
    }

  file_out.close();
  
  return all_vertex;
}

void DrawVoronoiDiagram(const VectR2& polygon)
{
  Real_wp xmin = 2.0, xmax = 7.0, ymin = -1.0, ymax = 12.0;
  int nb_x = 200, nb_y = 200;

  VectR2 normale(polygon.GetM()-1);
  for (int i = 0; i < polygon.GetM()-1; i++)
    {
      R2 ptB(polygon(i)), ptC(polygon(i+1));
      R2 vec_v = ptC - ptB; vec_v *= 1.0/Norm2(vec_v);
      normale(i).Init(-vec_v(1), vec_v(0));
    }

  VectReal_wp x_div, y_div;
  Linspace(xmin, xmax, nb_x, x_div);
  Linspace(ymin, ymax, nb_y, y_div);

  ofstream file_out("diag.dat");
  for (int i = 0; i < nb_x; i++)
    for (int j = 0; j < nb_y; j++)
      {
        R2 M(x_div(i), y_div(j));
        
        Real_wp dist_min = 1e300; int num = -1;
        if (InsidePolygon(M, polygon))
          for (int k = 0; k < polygon.GetM()-1; k++)
            {
              R2 x0 = polygon(k), x1 = polygon(k+1);
              
              // distance avec les aretes
              Real_wp alpha = DotProd(M-x0, normale(k));
              R2 ptP = M - alpha*normale(k);
              Real_wp longueur = AbsSquare(x1-x0);
              if ((AbsSquare(ptP-x0) <= longueur+1e-12)
                  && (AbsSquare(ptP-x1) <= longueur+1e-12))
                {
                  if (abs(alpha) < dist_min)
                    {
                      num = k;
                      dist_min = abs(alpha);
                    }
                }
              
              alpha = x0.Distance(M);
              if (alpha < dist_min)
                {
                  num = polygon.GetM()-1 + k;
                  dist_min = alpha;
                }
            }
        
        file_out << num << " ";
        if (j == nb_y-1)
          file_out << '\n';
      }

  file_out.close();        
}

int main(int argc, char**argv)
{
  InitMontjoie(argc, argv);

  if (argc != 2)
    {
      cout << "Donnez un fichier de donnees" << endl;
      abort();
    }

  string input_file(argv[1]);
  InputData var;
  
  ReadInputFile(input_file, var);
  
  GetVoronoiVertices(var.poly_line);
  //DrawVoronoiDiagram(var.poly_line);
  
  return FinalizeMontjoie();
}
