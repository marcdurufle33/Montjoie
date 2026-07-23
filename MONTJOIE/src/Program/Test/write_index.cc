#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

typedef Real_wp Complexe;
SplineInterpolation<Real_wp> spline_c, spline_rho, spline_pressure;
Real_wp coef_adim_p = 1.0/10.0;// passer de dyn/cm^2 a N/m^2
Real_wp coef_adim_c = 1.0/100; // passer de cm/s a m/s
Real_wp coef_adim_rho = 1.0e3; // passer de g/cm^3 a kg/m^3
//Real_wp  R =1.0007126;

Real_wp  R =1.0; //CHANGEMENT POUR CARRE !!!!


void Init ( string file_name)
{
  Matrix<double> val;
  val.ReadText(file_name);
  int N=val.GetM();
  VectReal_wp radius(N), density(N), celerity(N), pressure(N);
  for (int i = 0; i < N; i++)
    {
      double c = val(N-1-i, 1)*coef_adim_c;
      double r = val(N-1-i, 0);
      double rho = val(N-1-i, 2)*coef_adim_rho;
      Real_wp pre = val(N-1-i, 3)*coef_adim_p;

      radius(i)=r;
      celerity(i)=c;
      density(i)=rho;
      pressure(i) = pre;
    }

  spline_c.Init(radius,celerity);
  spline_rho.Init(radius, density);
  spline_pressure.Init(radius, pressure);

}

// mettre un boolean qui choisit entre splines et polynomes
void ComputeIndex(const Real_wp& x, const Real_wp& y, int n, Complexe& coef, bool poly)
{

  // n=0 : rho0.don
  // n=1 : c0.don
  // n=2 : mx.don
  // n=3 : my.don
  // n=4 : p0.don
  // n=7 : fx.don
  // n=8 : fy.don

  Real_wp r = sqrt(x*x+y*y);
  Real_wp rho, c0, mx, my, p, p_prime, p_seconde, fx, fy;//AJOUT NATHAN
  coef = 0.0;
  
  if (poly)
    {
      rho = 153.9*coef_adim_rho/cosh(r/0.09);
      UnivariatePolynomial<Real_wp> VectPoly(9);
      VectPoly(9) = -3.1618e+10;
      VectPoly(8) =  1.4151e+11;
      VectPoly(7) = -2.6703e+11;
      VectPoly(6) =  2.7569e+11;
      VectPoly(5) = -1.6830e+11;
      VectPoly(4) =  6.0535e+10;
      VectPoly(3) = -1.1626e+10;
      VectPoly(2) =  8.0179e+08;
      VectPoly(1) = -8.8174e+06;
      VectPoly(0) =  5.0673e+07;
      c0 = VectPoly.Evaluate(r);
      
      //    c = -3.1618e+10*pow(r,9)+ 1.4151e+11*pow(r,8)-2.6703e+11*pow(r,7)+ 2.7569e+11*pow(r,6) -1.6830e+11*pow(r,5)
      // + 6.0535e+10*pow(r,4) -1.1626e+10*pow(r,3) +8.0179e+08*pow(r,2)-8.8174e+06*pow(r,1)+5.0673e+07
      c0=c0*coef_adim_c;
      p = 4.1029*1e10*pow(rho,1.11); // rho*c0*c0;
      mx = -1.0*y/R;
      my = 1.0*x/R;
      p_prime = 4.1029*1e10*1.11*pow(r, 0.11);


      //AJOUT NATHAN
      //Calcul du second membre (no gravity)
      //Expression analytique de la dérivée
      // fx=rho(m_x\partical_xm_X+m_y\partial_ym_x)
      // fy = \rho(m_x\partial_xm_y/m_y\partial_ym_y)
      Real_wp deriv = 1.0/R;
      fx = rho*my*deriv+p_prime*(x/r);
      fy = -rho*mx*deriv+p_prime*(y/r); 

    }
  else
    {
      rho = spline_rho.Evaluate(r);
      c0 = spline_c.Evaluate(r); 
      mx = -c0*y/R; // 
      my = c0*x/R;

      // si on utilise Galbrun H1 il faut verifier l'equation d'etat 
      spline_pressure.EvaluateDerivative(r, p, p_prime, p_seconde);

      //AJOUT NATHAN
      //Calcul du second membre (no gravity)
      //Expression analytique de la dérivée
      // fx=rho(m_x\partical_xm_X+m_y\partial_ym_x)
      // fy = \rho(m_x\partial_xm_y/m_y\partial_ym_y)
      Real_wp deriv = c0/R;
      fx = rho*my*deriv+p_prime*(x/r);
      fy = -rho*mx*deriv+p_prime*(y/r); 
    }
  
  switch (n)
    {
    case 0 : coef = rho; break;
    case 1 : coef = c0; break;
    case 2 : coef = mx; break;
    case 3 : coef = my; break;
    case 4 : coef = p; break;
//    case 5 : coef = p_prime*x; break;//!!!
//    case 6 : coef = p_prime*y; break;//!!!
    //AJOUT NATHAN
    case 5 : coef = fx; break;
    case 6 : coef = fy; break;
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  // mettre le fichier en dur (le mettre dans un repertoire de MONTJOIE)
  string file_name("cptrho_no_leg");
  bool use_poly = false;
  //init splines
  Init(file_name);   
  if (argc < 3)
    {
      cout << "Usage : write_index_helio source destination" << endl;
      abort();
    }
  
  string name_source(argv[1]);
  Vector<string> name_dest(argc-2);
  for (int i = 2; i < argc; i++)
    name_dest(i-2) = string(argv[i]);
  
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
	ComputeIndex(points(i)(0), points(i)(1), j, nu_node(j)(i), use_poly);
      
      Complexe val;
      ComputeIndex(R, 0.0, j, val, use_poly);
      cout << "Index " << j << " for r=R : " << val << endl;
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
    WriteMeshData(nu_elt(j), name_dest(j));
  
  return FinalizeMontjoie();
}
