#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x)
{
  int type = rand()%3;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}

template<class T>
bool IsEqualData(const Vector<Vector<T> >& x, const Vector<Vector<T> >& y)
{
  if (x.GetM() != y.GetM())
    {
      cout << "Tailles de x et y different : " << x.GetM() << " et " << y.GetM() << endl;
      return false;
    }

  for (int i = 0; i < x.GetM(); i++)
    {
      if (x(i).GetM() != y(i).GetM())
	{
	  cout << "Tailles de x et y different : " << x(i).GetM() << " et " << y(i).GetM() << endl;
	  return false;
	}
      
      for (int j = 0; j < x(i).GetM(); j++)
	if (x(i)(j) != y(i)(j))
	  {
	    cout << "x et y different " << x(i)(j) << " et " << y(i)(j) << endl;
	    return false;
	  }
    }

  return true;
}

template<class T>
void CheckReadMeshData(Vector<Vector<T> >& nu)
{
  int N = 10;

  // on remplit aleatoirement nu
  nu.Reallocate(N);
  for (int i = 0; i < N; i++)
    {
      int nb_coef = rand()%20;
      nu(i).Reallocate(nb_coef);
      GenerateRandomVector(nu(i), nb_coef);
    }
  
  // on ecrit nu sur le disque
  WriteMeshData(nu, "data.don");
  
  // on lit nu dans un autre tableau
  Vector<Vector<T> > nu_copy;  
  ReadMeshData(nu_copy, "data.don");
  
  if (!IsEqualData(nu, nu_copy))
    {
      cout << "WriteMeshData ou ReadMeshData incorrect" << endl;
      abort();
    }
  
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  {
    Vector<Vector<double> > nu;
    CheckReadMeshData(nu);
  }

  {
    Vector<Vector<complex<double> > > nu;
    CheckReadMeshData(nu);
  }
  
  {
    Vector<Vector<float> > nu;
    CheckReadMeshData(nu);
  }

  {
    Vector<Vector<complex<float> > > nu;
    CheckReadMeshData(nu);
  }

  cout << "All tests passed successfully" << endl;

  return 0;
}
