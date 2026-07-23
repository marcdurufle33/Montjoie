#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

typedef Complex_wp Complexe;

// parametres de l'ellipse
Real_wp a = 0.3, b = 2.1;
Real_wp sigma = 0.02;
Real_wp alpha(1);

Real_wp GetFunction(const Real_wp& x, const Real_wp& y, const Real_wp& teta)
{
  return sqrt(square(x - a*cos(teta)) + square(y - b*sin(teta)));
}

Real_wp GetDerivative(const Real_wp& x, const Real_wp& y, const Real_wp& teta)
{
  Real_wp cos_teta = cos(teta), sin_teta = sin(teta);
  Real_wp dx = x - a*cos_teta, dy = y-b*sin_teta;
  Real_wp f = sqrt(dx*dx + dy*dy);
  return (dx*a*sin_teta - dy*b*cos_teta)/f;
}

void GetSecondDerivative(const Real_wp& x, const Real_wp& y, const Real_wp& teta, Real_wp& df, Real_wp& d2f)
{
  Real_wp cos_teta = cos(teta), sin_teta = sin(teta);
  Real_wp dx = x - a*cos_teta, dy = y-b*sin_teta;
  Real_wp f = sqrt(dx*dx + dy*dy);
  df = (dx*a*sin_teta - dy*b*cos_teta)/f;
  d2f = -square(dx*a*sin_teta - dy*b*cos_teta)/(f*f*f)
    + (square(a*sin_teta) - square(a*cos_teta) + x*a*cos_teta
       + square(b*cos_teta) - square(b*sin_teta) + b*y*sin_teta) / f;  
}


Real_wp GetDistance(const Real_wp& x, const Real_wp& y)
{
  // petit balayage pour partir d'un point assez bon
  int N = 10; Real_wp teta0(0), dist_min(1e300);
  for (int i = 0; i < N; i++)
    {
      Real_wp teta = 2*i*pi_wp/N;
      Real_wp dist = GetFunction(x, y, teta);
      if (dist < dist_min)
	{
	  dist_min = dist;
	  teta0 = teta;
	}
    }

  R2 center;
  EllipseMinimizer fct(x, y, center, a, b, teta0);
  TinyVector<Real_wp, 1> xsol;
  Real_wp dist = SolveLeastSquaresLvm(fct, xsol, 1e-16, -1, 30, 0, 0.1);

  // algo de Newton basique
  /* int k = 0; Real_wp df, d2f;
  int nb_iter_max = 20; Real_wp teta_next, teta_n = xsol(0);
  GetSecondDerivative(x, y, teta_n, df, d2f);
  Real_wp err = abs(df), err_prev = 1e300;
  while ((err_prev > err) && (k < nb_iter_max))
    {
      teta_next = teta_n - df/d2f;
      GetSecondDerivative(x, y, teta_next, df, d2f);

      err_prev = err;
      err = abs(df);
      teta_n = teta_next;
      k++; //DISP(k); DISP(err); DISP(teta_n);
    }

  dist_min = GetFunction(x, y, teta_n);
  return dist_min; */

  //DISP(dist); DISP(xsol);
  return dist;
}

void ComputeIndex(const Real_wp& x, const Real_wp& y, int n, Complexe& coef)
{
  Real_wp rho, c;
  Real_wp d = (x*x)/(a*a) + (y*y)/(b*b);
  if (d <= 1)
    {
      rho = 2.0;
      c = 1.0;
    }
  else
    {
      Real_wp dist = GetDistance(x, y);
      rho = 2.0*exp(-alpha*dist);
      //DISP(x); DISP(y); DISP(dist); DISP(rho);
      c = 1.0;
    }
  
  if (n == 0)
    {
      coef = 1.0/(rho*c*c)*Complex_wp(1, sigma);
    }
  else
    coef = 1.0/rho;
  
  
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  if (argc < 6)
    {
      cout << "Usage : write_index_helm alpha a b source destination" << endl;
      abort();
    }
  
  /* VectReal_wp xs, ys;
  Linspace(Real_wp(-6), Real_wp(6), 400, xs);
  Linspace(Real_wp(-6), Real_wp(6), 400, ys);

  Matrix<Real_wp> rho(xs.GetM(), ys.GetM());  
  for (int i = 0; i < xs.GetM(); i++)
    for (int j = 0; j < ys.GetM(); j++)
      {
	Complex_wp coef;
	ComputeIndex(xs(i), ys(j), 1, coef);
	rho(i, j) = Real_wp(1) / realpart(coef);
      }

      rho.WriteText("rho.dat"); */

  alpha = to_num<Real_wp>(string(argv[1]));
  sigma = to_num<Real_wp>(string(argv[2]));
  a = to_num<Real_wp>(string(argv[3]));
  b = to_num<Real_wp>(string(argv[4]));
  DISP(alpha);
  DISP(a); DISP(b);DISP(sigma);
  string name_source(argv[5]);
  Vector<string> name_dest(argc-6);
  for (int i = 6; i < argc; i++)
    name_dest(i-6) = string(argv[i]);
  
  ifstream file_in(name_source.data());
  int nb_elt, nb_nodes;
  file_in >> nb_elt >> nb_nodes;
  // on lit tous les points
  Vector<R2> points(nb_nodes);
  for (int i = 0; i < nb_nodes; i++)
    file_in >> points(i)(0) >> points(i)(1);  

  // on evalue l'indice sur tous les points (cas d'un indice continu)
  Vector<Vector<Complexe> > nu_node(name_dest.GetM());
  for (int j = 0; j < name_dest.GetM(); j++)
    {
      nu_node(j).Reallocate(nb_nodes);  
      for (int i = 0; i < nb_nodes; i++)
	ComputeIndex(points(i)(0), points(i)(1), j, nu_node(j)(i));      
    }
  
  // ensuite on remplit l'indice par element
  Vector<Vector<Vector<Complexe> > > nu_elt(name_dest.GetM());
  for (int j = 0; j < name_dest.GetM(); j++)
    nu_elt(j).Reallocate(nb_elt);
  
  int nb_nodes_elt, num_node;
  for (int i = 0; i < nb_elt; i++)
    {
      file_in >> nb_nodes_elt;
      for (int j = 0; j < name_dest.GetM(); j++)
	nu_elt(j)(i).Reallocate(nb_nodes_elt);
      
      for (int j = 0; j < nb_nodes_elt; j++)
	{
          file_in >> num_node;
          for (int k = 0; k < name_dest.GetM(); k++)
	    {
              nu_elt(k)(i)(j) = nu_node(k)(num_node);              
            }
	}
    }
  
  file_in.close();
  
  for (int j = 0; j < name_dest.GetM(); j++)
    {
      cout << "Ecriture de l'indice " << j << " sur le fichier " << name_dest(j) << endl;
      WriteMeshData(nu_elt(j), name_dest(j));
    }
  
  return FinalizeMontjoie();
}
