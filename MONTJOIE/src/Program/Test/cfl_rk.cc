#include "Solver/MontjoieSolver.hxx"
#include "Output/MontjoieOutput.hxx"
#include "Instationary/MontjoieTime.hxx"

using namespace Montjoie;

// we compute Pmod = |P(i \sqrt{x})|^2
void GetSquareModulus(const UnivariatePolynomial<Real_wp>& P,
		      UnivariatePolynomial<Real_wp>& Pmod)
{
  UnivariatePolynomial<Real_wp> Pr, Pi, z;
  z.SetOrder(1); z(0) = Real_wp(0); z(1) = Real_wp(1);
  
  // computing Pr(i \sqrt{z}) (to avoid high degrees)
  Pr.SetOrder(P.GetOrder()/2);  
  for (int i = 0; i <= P.GetOrder(); i += 2)
    {
      Pr(i/2) = P(i);
      if (i%4 == 2)
	Pr(i/2) = -P(i);
    }

  // computing \sqrt{z} Pi(i \sqrt{z}) (to avoid high degrees)
  Pi.SetOrder((P.GetOrder()-1)/2);  
  for (int i = 1; i <= P.GetOrder(); i += 2)
    {
      Pi((i-1)/2) = P(i);
      if (i%4 == 3)
	Pi((i-1)/2) = -P(i);
    }

  // modulus |P(i sqrt(z))|^2
  Pmod = Pr*Pr + z*Pi*Pi;
  
}

// returns the CFL associated with the stability function R = P / Q
Real_wp GetCFL(const UnivariatePolynomial<Real_wp>& P,
               bool display = false)
{
  // finding the values z for which |P|^2 = |Q|^2
  UnivariatePolynomial<Real_wp> Pmod, one, Pdiff;

  GetSquareModulus(P, Pmod);
  one.SetOrder(0); one(0) = Real_wp(1);
  
  if (display) { DISP(Pmod); }

  // we search the roots of |P|^2 - 1
  UnivariatePolynomial<Real_wp> R;
  VectReal_wp Li, Lr;
  R = Pmod - one;
  
  for (int i = 0; i <= R.GetOrder(); i++)
    if (abs(R(i)) <= 1e4*epsilon_machine)
      R(i) = Real_wp(0);

  // checking the sign of P-1 at origin
  for (int i = 0; i <= R.GetOrder(); i++)
    if (R(i) != Real_wp(0))
      {
	if (R(i) > Real_wp(0))
	  return Real_wp(0);
	
	break;
      }
  
  //DISP(epsilon_machine);
  SolvePolynomialEquation(R, Lr, Li);

  Sort(Lr, Li);
  if (display) { DISP(R); DISP(Lr); DISP(Li); }
  
  Real_wp cfl(1e300);
  int num_root = 0;
  while (num_root < Lr.GetM())
    {
      // we select only real roots different from 0
      if ((abs(Li(num_root)) <= epsilon_machine) && (Lr(num_root) >= Real_wp(1e-12)))
	{
	  int multiplicity = 1;
	  while ((num_root < Lr.GetM()-1) && (abs(Lr(num_root+1) - Lr(num_root)) <= 1e4*epsilon_machine))
	    {
	      multiplicity++;
	      num_root++;
	    }
          
          // for odd multiplicity, the CFL is lower than this root
	  if (multiplicity%2 == 1)
	    cfl = min(cfl, sqrt(Lr(num_root)));
	}

      num_root++;
    }

  return cfl;
}

void GetProfileCabane(int N, VectR2& points)
{
  VectReal_wp step_t;
  Linspace(Real_wp(0), Real_wp(1), N, step_t);

  points.Reallocate(3*N);
  ofstream file_out("points.dat");
  for (int i = 0; i < N; i++)
    {
      points(i).Init(Real_wp(0), step_t(i));
      points(i+N).Init(-step_t(i), Real_wp(1));
      points(i+2*N).Init(step_t(i)-Real_wp(2), step_t(i)*(Real_wp(14)-4*step_t(i))/10);
    }

  for (int i = 0; i < points.GetM(); i++)
    file_out << points(i)(0) << " " << points(i)(1) << '\n';

  file_out.close();
}

bool StableSchemeForPoint(const UnivariatePolynomial<Real_wp>& P,
			  const Real_wp& h, const VectR2& points)
{
  for (int i = 0; i < points.GetM(); i++)
    {
      Complex_wp val = P.Evaluate(h*Complex_wp(points(i)(0), points(i)(1)));
      if (abs(val) > Real_wp(1) + 1e4*epsilon_machine)
	return false;
    }

  return true;
}


Real_wp GetCFL_Cabane(const UnivariatePolynomial<Real_wp>& P)
{
  int N = 11;
  VectR2 points;
  GetProfileCabane(N, points);

  // first evaluation of an interval [h0, h1] for the CFL
  Real_wp h0(0), h1(1);
  bool scheme_stable(true);
  while (scheme_stable)
    {
      if (StableSchemeForPoint(P, h1, points))
	{
	  h0 = h1;
	  h1 *= Real_wp(2);
	}
      else
	scheme_stable = false;
    }

  // DISP(h0); DISP(h1);

  // then bisection method
  while (abs(h1-h0) > 1e-12)
    {
      Real_wp hm = (h0+h1)/2;
      scheme_stable = StableSchemeForPoint(P, hm, points);
      if (scheme_stable)
	h0 = hm;
      else
	h1 = hm;
    }
  
  return h0;
}

void GenerateMesh2D_Stab(int Nx, int Ny,
			 const UnivariatePolynomial<Real_wp>& P, 
			 const Real_wp& xmin, const Real_wp& xmax,
			 const Real_wp& ymin, const Real_wp& ymax)
{
  VectReal_wp x_div, y_div;
  Linspace(xmin, xmax, Nx, x_div);
  Linspace(ymin, ymax, Ny, y_div);

  
  Matrix<int> IndexVertex(Nx, Ny);
  IndexVertex.Fill(-1);
  int nb_vert = 0;
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      {
	Complex_wp val = P.Evaluate(Complex_wp(x_div(i), y_div(j)));
	if (abs(val) <= Real_wp(1)+1e4*epsilon_machine)
	  IndexVertex(i, j) = nb_vert++;
      }
  
  Mesh<Dimension2> mesh;
  mesh.ReallocateVertices(nb_vert);
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      if (IndexVertex(i, j) >= 0)
	{
	  int num = IndexVertex(i, j);
	  mesh.Vertex(num).Init(x_div(i), y_div(j));
	}
  
  Matrix<int> IndexQuad(Nx, Ny);
  Matrix<int> IndexTri(Nx, Ny);
  IndexQuad.Fill(-1); IndexTri.Fill(-1);
  int nb_quad = 0, nb_tri = 0;
  for (int i = 0; i < Nx-1; i++)
    for (int j = 0; j < Ny-1; j++)
      {
	nb_vert = 0;
	if (IndexVertex(i, j) >= 0)
	  nb_vert++;
	
	if (IndexVertex(i+1, j) >= 0)
	  nb_vert++;
	
	if (IndexVertex(i, j+1) >= 0)
	  nb_vert++;
	
	if (IndexVertex(i+1, j+1) >= 0)
	  nb_vert++;
	
	if (nb_vert == 4)
	  IndexQuad(i, j) = nb_quad++;
	else if (nb_vert == 3)
	  IndexTri(i, j) = nb_tri++;
      }
  
  mesh.ReallocateElements(nb_quad + nb_tri);
  IVect num_pt(4), nump(4);
  for (int i = 0; i < Nx-1; i++)
    for (int j = 0; j < Ny-1; j++)
      {
	num_pt(0) = IndexVertex(i, j);
	num_pt(1) = IndexVertex(i+1, j);
	num_pt(2) = IndexVertex(i+1, j+1);
	num_pt(3) = IndexVertex(i, j+1);
	
	if (IndexTri(i, j) >= 0)
	  {
	    nump = num_pt;
	    Sort(nump);
	    int num = IndexTri(i, j);
	    mesh.Element(num).InitTriangular(nump(1), nump(2), nump(3), 1);
	  }
	else if (IndexQuad(i, j) >= 0)
	  {
	    int num = nb_tri + IndexQuad(i, j);
	    mesh.Element(num).InitQuadrangular(num_pt(0), num_pt(1), num_pt(2), num_pt(3), 1);
	  }	
      }
  
  mesh.Write("test.mesh");
}


bool ProjectPointStabilityRegion(const UnivariatePolynomial<Real_wp>& P, R2& pt,
				 const R2& normale)
{
  Complex_wp val = P.Evaluate(Complex_wp(pt(0), pt(1)));
  if (abs(val) > Real_wp(1) + 1e4*epsilon_machine)
    {
      pt -= Real_wp(1e-12)*normale;
      val = P.Evaluate(Complex_wp(pt(0), pt(1)));
    }
  
  if (abs(val) > Real_wp(1) + 1e4*epsilon_machine)
    {
      cout << "Impossible " << endl;
      cout << "the point " << pt << " is not in the stability domain" << endl;
      abort();
    }

  Real_wp coef(1); bool test_loop = true;
  R2 ptA = pt, ptB = pt + normale;
  while (test_loop)
    {
      ptB = pt + coef*normale;
      val = P.Evaluate(Complex_wp(ptB(0), ptB(1)));
      if (abs(val) <= Real_wp(1) + 1e4*epsilon_machine)
	coef *= Real_wp(2);
      else
	test_loop = false;
      
      if (coef > Real_wp(16))
	{
	  cout << "Point not on the boundary ? " << pt << endl;
	  DISP(pt); DISP(normale); DISP(coef);
	  return false;
	}
    }
  
  // bisection method between ptA and ptB
  R2 ptC;
  while (ptA.Distance(ptB) > 1e-12)
    {
      ptC = Real_wp(0.5)*(ptA + ptB);
      val = P.Evaluate(Complex_wp(ptC(0), ptC(1)));
      if (abs(val) > Real_wp(1) + 1e4*epsilon_machine)
	ptB = ptC;
      else
	ptA = ptC;
    }
  
  pt = ptA;
  return true;
}

void RefineBoundary(const Mesh<Dimension2>& mesh, const UnivariatePolynomial<Real_wp>& P,
		    Mesh<Dimension2>& mesh_dest)
{
  // main loop over points
  Vector<bool> VertexTreated(mesh.GetNbVertices());
  VertexTreated.Fill(false);
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int n0 = mesh.BoundaryRef(i).numVertex(0);
      int n1 = mesh.BoundaryRef(i).numVertex(1);
      R2 pt0 = mesh.Vertex(n0), pt1 = mesh.Vertex(n1);
      R2 normale(pt1(1)-pt0(1), pt0(0)-pt1(0));
      R2 milieu = Real_wp(0.5)*(pt0 + pt1);
      
      int num_elem = mesh.BoundaryRef(i).numElement(0);
      R2 center;
      int nb_vert = mesh.Element(num_elem).GetNbVertices();
      for (int j = 0; j < nb_vert; j++)
	center += mesh.Vertex(mesh.Element(num_elem).numVertex(j));
          
      Mlt(Real_wp(1)/nb_vert, center);
      R2 diff = milieu - center;
      if (DotProd(diff, normale) < 0)
	normale = -normale;

      if (!VertexTreated(n0))
	{
	  VertexTreated(n0) = true;
	  bool pt_found = ProjectPointStabilityRegion(P, pt0, normale);
	  if (!pt_found)
	    VertexTreated(n0) = false;
	  
	  mesh_dest.Vertex(n0) = pt0;                  
	}

      if (!VertexTreated(n1))
	{
	  VertexTreated(n1) = true;
	  bool pt_found = ProjectPointStabilityRegion(P, pt1, normale);
	  if (!pt_found)
	    VertexTreated(n1) = false;
	  
	  mesh_dest.Vertex(n1) = pt1;
	}
    }
  
  mesh_dest.ForceCoherenceMesh(true);
  
  //int test_input; cout << "waiting" << endl; cin >> test_input;
}


void ComputeStabilityDomain(const UnivariatePolynomial<Real_wp>& P,
			    Mesh<Dimension2>& mesh_dest)
{
  cout << "Entrez xmin xmax ymin ymax" << endl;
  Real_wp xmin, xmax, ymin, ymax;
  cin >> xmin >> xmax >> ymin >> ymax;
  int Nx, Ny;
  cout << "Entrez Nx Ny" << endl;
  cin >> Nx >> Ny;
  
  GenerateMesh2D_Stab(Nx, Ny, P, xmin, xmax, ymin, ymax);

  Mesh<Dimension2> mesh;
  mesh.Read("test.mesh");

  mesh_dest.Read("test.mesh");
  mesh_dest.ClearElements();
  
  RefineBoundary(mesh, P, mesh_dest);
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc != 3)
    {
      cout << "Entrez ordre et nombre d'etapes additionnelles" << endl;
      cout << "Usage ./cfl.x 4 2" << endl;
      abort();
    }

  int order = atoi(argv[1]);
  int l = atoi(argv[2]);
  
  LinearRungeKutta_Iterator<Real_wp> rk;
  rk.SetOrder(order, l);

  const UnivariatePolynomial<Real_wp>& P = rk.GetStabilityFunction();

  cout << "Que voulez-vous faire pour ERK" << order << "-" << l << " ? " << endl;
  cout << "1- Calculer la CFL sur l'axe imaginaire" << endl;
  cout << "2- Calculer la CFL de la cabane " << endl;
  cout << "3- Calculer le contour du domaine de stabilite " << endl;
  int reponse(0); cin >> reponse;
  
  if (reponse == 1)
    {
      Real_wp cfl = GetCFL(P);
      cout << "CFL on the imaginary axis = " << cfl  << endl;
      cout << "Efficiency = " << cfl / (order + l) << endl;
    }
  else if (reponse == 2)
    {
      Real_wp cfl = GetCFL_Cabane(P);
      cout << "CFL for the cabane profile = " << cfl << endl;
      cout << "Efficiency = " << cfl / (order + l) << endl;
    } 
  else if (reponse == 3)
    {
      Mesh<Dimension2> mesh;
      ComputeStabilityDomain(P, mesh);
      cout << "Writing the contour in contour.mesh" << endl;
      mesh.Write("contour.mesh");
    }
  
  return 0;
}
  
