#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

template<class T>
void GeneratePermutation(const Vector<T>& x, Vector<Vector<T> >& all_perm)
{
  // recursive algorithm
  if (x.GetM() == 1)
    {
      all_perm.Reallocate(1);
      all_perm(0) = x;
    }
  else if (x.GetM() == 2)
    {
      all_perm.Reallocate(2);
      all_perm(1).Reallocate(2);
      all_perm(0) = x;
      all_perm(1)(0) = x(1); all_perm(1)(1) = x(0);
    }
  else
    {
      int offset = 0;
      for (int i = 0; i < x.GetM(); i++)
	{
	  Vector<T> x_reduce(x.GetM()-1);
	  for (int j = 0; j < i; j++)
	    x_reduce(j) = x(j);

	  for (int j = i+1; j < x.GetM(); j++)
	    x_reduce(j-1) = x(j);

	  Vector<Vector<T> > reduce_perm;
	  GeneratePermutation(x_reduce, reduce_perm);

	  if (i == 0)
	    all_perm.Reallocate(x.GetM()*reduce_perm.GetM());

	  for (int j = 0; j < reduce_perm.GetM(); j++)
	    {
	      all_perm(offset+j).Reallocate(x.GetM());
	      all_perm(offset+j)(0) = x(i);
	      for (int k = 0; k < x.GetM()-1; k++)
		all_perm(offset+j)(k+1) = reduce_perm(j)(k);
	    }

	  offset += reduce_perm.GetM();
	}
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  // testing string stuff
  string nom("fileTest1.5toto.08.dat");
  if (GetExtension(nom) != "dat")
    {
      cout << "GetExtension incorrect" << endl;
      abort();
    }
  
  if (GetBaseString(nom) != "fileTest1.5toto.08")
    {
      cout << "GetBaseString incorrect" << endl;
      abort();
    }
  
  nom = string("SansNom");
  if (GetExtension(nom).size() > 0)
    {
      cout << "GetExtension incorrect" << endl;
      abort();
    }

  if (GetBaseString(nom) != "SansNom")
    {
      cout << "GetBaseString incorrect" << endl;
      abort();
    }
  
  nom = string("value1.0_@#.dat totoY -9");
  Vector<string> param;
  StringTokenize(nom, param, " ");
  if ((param.GetM() != 3) || (param(0) != "value1.0_@#.dat") || (param(1) != "totoY") || (param(2) != "-9"))
    {
      cout << "StringTokenize incorrect" << endl;
      abort();
    }

  nom = string(" value1 .0_@#.dat \t  totoY \t-9  ");
  StringTokenize(nom, param, " \t");
  if ((param.GetM() != 4) || (param(0) != "value1") || (param(1) != ".0_@#.dat") || (param(2) != "totoY") || (param(3) != "-9"))
    {
      cout << "StringTokenize incorrect" << endl;
      abort();
    }
  
  nom = string("OnlyOneParam");
  StringTokenize(nom, param, " \t");
  if ((param.GetM() != 1) && (param(0) != nom))
    {
      cout << "StringTokenize incorrect" << endl;
      abort();
    }

  nom = string("   \t \t   ");
  StringTokenize(nom, param, " \t");
  if (param.GetM() != 0)
    {
      cout << "StringTokenize incorrect" << endl;
      abort();
    }
  
  nom = string("Nom TOto \tprod \t");
  StringTrim(nom);
  if (nom != "NomTOtoprod")
    {
      cout << "StringTrim incorrect" << endl;
      abort();
    }

  StringTrim(nom);
  if (nom != "NomTOtoprod")
    {
      cout << "StringTrim incorrect" << endl;
      abort();
    }

  nom = string("\t\t   \t \t");
  StringTrim(nom);
  if (nom.size() > 0)
    {
      cout << "StringTrim incorrect" << endl;
      abort();
    }

  nom = string("  MyValue \t Prend  cette \t\t ");
  DeleteSpaceAtExtremityOfString(nom);
  if (nom != "MyValue \t Prend  cette")
    {
      cout << "DeleteSpaceAtExtremityOfString incorrect" << endl;
      abort();
    }
  
  nom = string("NoSpaceString");
  DeleteSpaceAtExtremityOfString(nom);
  if (nom != "NoSpaceString")
    {
      cout << "DeleteSpaceAtExtremityOfString incorrect" << endl;
      abort();
    }

  nom = string("   \t \t   ");
  DeleteSpaceAtExtremityOfString(nom);
  if (nom.size() > 0)
    {
      cout << "StringTokenize incorrect" << endl;
      abort();
    }
  
  if (NumberToString(2) != "0002")
    {
      cout << "NumberToString incorrect" << endl;
      abort();
    }

  if (NumberToString(34) != "0034")
    {
      cout << "NumberToString incorrect" << endl;
      abort();
    }

  if (NumberToString(328) != "0328")
    {
      cout << "NumberToString incorrect" << endl;
      abort();
    }

  if (NumberToString(1238) != "1238")
    {
      cout << "NumberToString incorrect" << endl;
      abort();
    }

  if (NumberToString(19234, 7) != "0019234")
    {
      cout << "NumberToString incorrect" << endl;
      abort();
    }
  
  // testing Copy
  Vector<string> xvec, xref; list<string> xl;
  xl.push_back(string("toto")); xref.PushBack(string("toto"));
  xl.push_back(string("other")); xref.PushBack(string("other"));
  xl.push_back(string("nom")); xref.PushBack(string("nom"));
  Copy(xl, xvec);
  
  if (xvec.GetM() != xref.GetM())
    {
      cout << "Copy incorrect" << endl;
      abort();
    }
  else
    {
      for (int j = 0; j < xvec.GetM(); j++)
	if (xvec(j) != xref(j))
	  {
	    cout << "Copy incorrect" << endl;
	    abort();
	  }
    }
  
  xl.clear(); xref.Clear();
  Copy(xl, xvec);
  
  if (xvec.GetM() != 0)
    {
      cout << "Copy incorrect" << endl;
      abort();
    }
  
  Vector<double, VectSparse> x_sparse(4), x_sp;
  Vector<double> x_dense(20);
  x_sparse.Index(0) = 1; x_sparse.Value(0) = 1.2;
  x_sparse.Index(1) = 4; x_sparse.Value(1) = -0.5;
  x_sparse.Index(2) = 9; x_sparse.Value(2) = 2.6;
  x_sparse.Index(3) = 13; x_sparse.Value(3) = -3.2;
  x_dense.Fill(0);
  for (int i = 0; i < x_sparse.GetM(); i++)
    x_dense(x_sparse.Index(i)) = x_sparse.Value(i);
  
  Copy(x_dense, x_sp);
  if (x_sparse.GetM() != x_sp.GetM())
    {
      cout << "Copy incorrect" << endl;
      abort();
    }
  else
    {
      for (int i = 0; i < x_sp.GetM(); i++)
	if ((x_sp.Index(i) != x_sparse.Index(i)) || (x_sp.Value(i) != x_sparse.Value(i)))
	  {
	    cout << "Copy incorrect" << endl;
	    abort();
	  }
    }
  
  DISP(GetHumanReadableMemory(int64_t(12300)));
  DISP(GetHumanReadableMemory(int64_t(123000)));
  DISP(GetHumanReadableMemory(int64_t(1230000)));
  DISP(GetHumanReadableMemory(int64_t(12300000)));
  DISP(GetHumanReadableMemory(int64_t(123000000)));
  DISP(GetHumanReadableMemory(int64_t(1230000000)));
  DISP(GetHumanReadableMemory(int64_t(12300000000)));
  DISP(GetHumanReadableMemory(int64_t(123000000000)));
  DISP(GetHumanReadableMemory(int64_t(1230000000000)));
  DISP(GetHumanReadableMemory(int64_t(12300000000000)));

  // testing sort
  Vector<int> num(2); num.FillRand();
  Vector<Vector<int> > all_perm;
  GeneratePermutation(num, all_perm);

  Vector<int> num_res(num);
  Sort(num_res);
  
  for (int i = 0; i < all_perm.GetM(); i++)
    {
      int a = all_perm(i)(0);
      int b = all_perm(i)(1);
      Sort(a, b);
      if ((a != num_res(0)) || (b != num_res(1)))
	{
	  cout << "Sort incorrect" << endl;
	  abort();
	}
    }

  num.Reallocate(3); num.FillRand();
  GeneratePermutation(num, all_perm);

  num_res = num;
  Sort(num_res);
  
  for (int i = 0; i < all_perm.GetM(); i++)
    {
      int a = all_perm(i)(0);
      int b = all_perm(i)(1);
      int c = all_perm(i)(2);
      Sort(a, b, c);
      if ((a != num_res(0)) || (b != num_res(1)) || (c != num_res(2)))
	{
	  cout << "Sort incorrect" << endl;
	  abort();
	}
    }

  num.Reallocate(4); num.FillRand();
  GeneratePermutation(num, all_perm);

  num_res = num;
  Sort(num_res);
  
  for (int i = 0; i < all_perm.GetM(); i++)
    {
      int a = all_perm(i)(0);
      int b = all_perm(i)(1);
      int c = all_perm(i)(2);
      int d = all_perm(i)(3);
      Sort(a, b, c, d);
      if ((a != num_res(0)) || (b != num_res(1)) || (c != num_res(2)) || (d != num_res(3)))
	{
	  cout << "Sort incorrect" << endl;
	  abort();
	}
    }

  // testing conversion of coordinates  
  for (int n = 0; n < 20; n++)
    {
      Real_wp x = 2.0*Real_wp(rand())/RAND_MAX-1.0;
      Real_wp y = 2.0*Real_wp(rand())/RAND_MAX-1.0;
      Real_wp z = 2.0*Real_wp(rand())/RAND_MAX-1.0;
      Real_wp r, theta;
      
      if (n == 0) { x = 0; y = 0; z = 0;}
      if (n == 1) { x = 0; y = 0;}
      if (n == 2) { x = 0; z = 0;}
      if (n == 3) { y = 0; z = 0;}
      if (n == 4) { x = 0; }
      if (n == 5) { y = 0; }
      if (n == 6) { z = 0; }
      
      CartesianToPolar(x, y, r, theta);
      
      if ((abs(x-r*cos(theta)) > 1e-12) || isnan(r) || isnan(theta))
	{
	  DISP(x); DISP(y); DISP(r); DISP(theta);
	  cout << "CartesianToPolar incorrect" << endl;
	  abort();
	}

      Real_wp phi, cos_theta, sin_theta;
      R2 pt;
      SetIncidentAngle(r, pt, theta, phi);

      if ((abs(x-pt(0)) > 1e-12) || isnan(pt(0)) || (abs(y-pt(1)) > 1e-12) || isnan(pt(1)))
	{
	  cout << "SetIncidentAngle incorrect" << endl;
	  abort();
	}      

      CartesianToSpherical(x, y, z, r, theta, phi, cos_theta, sin_theta);

      if ((abs(x-r*sin(theta)*cos(phi)) > 1e-12) || isnan(r) || isnan(theta) || isnan(phi) ||
	  (abs(y-r*sin(theta)*sin(phi)) > 1e-12) || (abs(cos_theta-cos(theta)) > 1e-12) ||
	  (abs(z-r*cos(theta)) > 1e-12) || (abs(sin_theta-sin(theta)) > 1e-12))
	{
	  cout << "CartesianToSpherical incorrect" << endl;
	  abort();
	}

      SetComplexZero(theta); SetComplexZero(phi);
      CartesianToSpherical(x, y, z, r, theta, phi);

      if ((abs(x-r*sin(theta)*cos(phi)) > 1e-12) || isnan(r) || isnan(theta) || isnan(phi) ||
	  (abs(y-r*sin(theta)*sin(phi)) > 1e-12) || (abs(z-r*cos(theta)) > 1e-12))
	{
	  cout << "CartesianToSpherical incorrect" << endl;
	  abort();
	}

      R3 pt3;
      SetIncidentAngle(r, pt3, theta, phi);

      if ((abs(x-pt3(0)) > 1e-12) || isnan(pt3(0)) || (abs(y-pt3(1)) > 1e-12) || isnan(pt3(1)) || (abs(z-pt3(2)) > 1e-12) || isnan(pt3(2)))
	{
	  cout << "SetIncidentAngle incorrect" << endl;
	  abort();
	}      

      SphericalToCartesian(r, theta, phi, x, y, z);

      if ((abs(x-pt3(0)) > 1e-12) || isnan(pt3(0)) || (abs(y-pt3(1)) > 1e-12) || isnan(pt3(1)) || (abs(z-pt3(2)) > 1e-12) || isnan(pt3(2)))
	{
	  cout << "SphericalToCartesian incorrect" << endl;
	  abort();
	}      

    }
  

  // testing linspace
  Real_wp xmin = Real_wp(rand()) / RAND_MAX;
  Real_wp xmax = xmin + Real_wp(rand()) / RAND_MAX;
  int N = rand()%500+2;
  VectReal_wp x_sub;
  Linspace(xmin, xmax, N, x_sub);
  
  if ((x_sub.GetM() != N) || (x_sub(0) != xmin) || (x_sub(N-1) != xmax))
    {
      cout << "Linspace incorrect" << endl;
      abort();
    }

  Real_wp dx = x_sub(1) - x_sub(0);
  for (int i = 2; i < x_sub.GetM(); i++)
    if (abs(dx - (x_sub(i) - x_sub(i-1))) > 1e-12)
      {
	cout << "Linspace incorrect" << endl;
	abort();      
      }

  xmin = 10.0*Real_wp(rand()) / RAND_MAX - 5.0;
  xmax = xmin + 10.0*Real_wp(rand()) / RAND_MAX;
  Real_wp ymin = 10.0*Real_wp(rand()) / RAND_MAX - 5.0;
  Real_wp ymax = ymin + 10.0*Real_wp(rand()) / RAND_MAX;
  Real_wp zmin = 10.0*Real_wp(rand()) / RAND_MAX - 5.0;
  Real_wp zmax = zmin + 10.0*Real_wp(rand()) / RAND_MAX;

  R2 kwave2D;
  kwave2D(0) = 4.0*Real_wp(rand()) / RAND_MAX - 2.0;
  kwave2D(1) = 4.0*Real_wp(rand()) / RAND_MAX - 2.0;

  Real_wp k0_dot_u = GetMinimumDotProd(kwave2D, xmin, xmax, ymin, ymax, zmin, zmax);
  Real_wp k1_dot_u = GetMaximumDotProd(kwave2D, xmin, xmax, ymin, ymax, zmin, zmax);
  VectReal_wp k0; k0.PushBack(DotProd(kwave2D, R2(xmin, ymin)));
  k0.PushBack(DotProd(kwave2D, R2(xmax, ymin)));
  k0.PushBack(DotProd(kwave2D, R2(xmax, ymax)));
  k0.PushBack(DotProd(kwave2D, R2(xmin, ymax)));
  Real_wp kmin(1e300), kmax(-1e300);
  for (int i = 0; i < k0.GetM(); i++)
    {
      kmin = min(kmin, k0(i));
      kmax = max(kmax, k0(i));
    }

  if (abs(k0_dot_u - kmin) > 1e-12)
    {
      cout << "GetMinimumDotProd incorrect" << endl;
      abort();
    }
  
  if (abs(k1_dot_u - kmax) > 1e-12)
    {
      cout << "GetMaximumDotProd incorrect" << endl;
      abort();
    }
  
  R3 kwave3D;
  kwave3D(0) = 4.0*Real_wp(rand()) / RAND_MAX - 2.0;
  kwave3D(1) = 4.0*Real_wp(rand()) / RAND_MAX - 2.0;
  kwave3D(2) = 4.0*Real_wp(rand()) / RAND_MAX - 2.0;
  
  k0_dot_u = GetMinimumDotProd(kwave3D, xmin, xmax, ymin, ymax, zmin, zmax);
  k1_dot_u = GetMaximumDotProd(kwave3D, xmin, xmax, ymin, ymax, zmin, zmax);
  k0.Clear();
  k0.PushBack(DotProd(kwave3D, R3(xmin, ymin, zmin)));
  k0.PushBack(DotProd(kwave3D, R3(xmax, ymin, zmin)));
  k0.PushBack(DotProd(kwave3D, R3(xmax, ymax, zmin)));
  k0.PushBack(DotProd(kwave3D, R3(xmin, ymax, zmin)));
  k0.PushBack(DotProd(kwave3D, R3(xmin, ymin, zmax)));
  k0.PushBack(DotProd(kwave3D, R3(xmax, ymin, zmax)));
  k0.PushBack(DotProd(kwave3D, R3(xmax, ymax, zmax)));
  k0.PushBack(DotProd(kwave3D, R3(xmin, ymax, zmax)));
  kmin= 1e300; kmax = -1e300;
  for (int i = 0; i < k0.GetM(); i++)
    {
      kmin = min(kmin, k0(i));
      kmax = max(kmax, k0(i));
    }

  if (abs(k0_dot_u - kmin) > 1e-12)
    {
      cout << "GetMinimumDotProd incorrect" << endl;
      abort();
    }

    if (abs(k1_dot_u - kmax) > 1e-12)
      {
      cout << "GetMaximumDotProd incorrect" << endl;
      abort();
    }

  cout << "All tests passed successfully" << endl;

  return FinalizeMontjoie();
}
