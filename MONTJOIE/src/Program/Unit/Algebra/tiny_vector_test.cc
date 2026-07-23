#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

Real_wp threshold;

template<class T, int m>
bool EqualVector(const TinyVector<T, m>& x, const TinyVector<T, m>& y)
{
  for (int i = 0; i < m; i++)
    if (isnan(abs(x(i)-y(i))) || (abs(x(i)-y(i)) > threshold))
      return false;

  return true;
}

template<class T, int m>
void CheckComparisonTinyVector(TinyVector<T, m>& x)
{
  TinyVector<T, m> y = x;
  
  y(0) += 0.1*threshold;
  y(1) -= 0.1;
  y(2) += 0.2;
  
  if ((y < x) && ( y <= x))
    {
    }
  else
    {
      cout << "Operators <, <= incorrect" << endl;
      abort();
    }
  
  if ((y > x) || (y >= x))
    {
      cout << "Operators >, >= incorrect" << endl;
      abort();
    }
  
  y = x;
  y(0) += 0.1*threshold;
  y(1) -= 0.14*threshold;
  y(2) += 0.23*threshold;
  
  if ((y >= x) && ( y <= x))
    {
    }
  else
    {
      cout << "Operators <=, >= incorrect" << endl;
      abort();
    }
  
  if ((y > x) || (y < x))
    {
      cout << "Operators >, < incorrect" << endl;
      abort();
    }

  y = x;
  y(0) += 0.5;
  y(1) -= 0.14*threshold;
  y(2) -= 0.23;
  
  if ((y >= x) && ( y > x))
    {
    }
  else
    {
      cout << "Operators >=, > incorrect" << endl;
      abort();
    }
  
  if ((y < x) || (y <= x))
    {
      cout << "Operators <, <= incorrect" << endl;
      abort();
    }
}

template<class T, int m>
void CheckComparisonTinyVector(TinyVector<complex<T>, m>& x)
{
  // no comparison possible for complex numbers
}


template<class T, int m>
void CheckScalarProductTinyVector(TinyVector<complex<T>, m>& x, TinyVector<complex<T>, m>& y)
{
  
  for (int i = 0; i < m; i++)
    {
      x(i) = complex<T>(rand(), rand());
      y(i) = complex<T>(rand(), rand());
    }
  
  Mlt(Real_wp(1e-8), x); Mlt(Real_wp(1e-8), y);

  // checking DotProd
  complex<T> scal = DotProd(x, y);
  complex<T> scal_ex(0, 0);
  for (int i = 0; i < m; i++)
    scal_ex += x(i)*y(i);
  
  if ((abs(scal - scal_ex) > threshold) || isnan(abs(scal-scal_ex)))
    {
      cout << "DotProd incorrect" << endl;
      abort();
    }

  // checking DotProdConj
  scal = DotProdConj(x, y);
  scal_ex = complex<T>(0, 0);
  for (int i = 0; i < m; i++)
    scal_ex += conj(x(i))*y(i);
  
  if ((abs(scal - scal_ex) > threshold)  || isnan(abs(scal-scal_ex)))
    {
      cout << "DotProdConj incorrect" << endl;
      abort();
    }
  
  // checking Norm2, Distance and DistanceSquare
  T module = sqrt(abs(DotProdConj(x, x)));
  if ((abs(module - Norm2(x)) > threshold) || isnan(abs(module - Norm2(x))))
    {
      cout << "Norm2 incorrect " << endl;
      abort();
    }

  // checking AbsSquare
  module = abs(DotProdConj(x, x));
  if ((abs(module - AbsSquare(x)) > threshold) || isnan(abs(module - AbsSquare(x))))
    {
      cout << "AbsSquare incorrect" << endl;
      abort();
    }
  
  TinyVector<complex<T>, m> diff;
  diff = x-y;
  module = Norm2(diff);
  if ((abs(module - x.Distance(y)) > threshold) || isnan(abs(module - x.Distance(y))))
    {
      cout << "Distance incorrect " << endl;
      abort();
    }

  if ((abs(module - Distance(x, y)) > threshold) || isnan(abs(module - Distance(x, y))))
    {
      cout << "Distance incorrect " << endl;
      abort();
    }

  if ((abs(module*module - x.DistanceSquare(y)) > threshold) || isnan(abs(module*module - x.DistanceSquare(y))))
    {
      cout << "DistanceSquare incorrect " << endl;
      abort();
    }
  
  // checking Norm1
  scal_ex = complex<T>(0, 0);
  for (int i = 0; i < m; i++)
    scal_ex += abs(x(i));
  
  scal = complex<T>(Norm1(x), 0);
  if ((abs(scal - scal_ex) > threshold) || isnan(abs(scal - scal_ex)))
    {
      cout << "Norm1 incorrect" << endl;
      abort();
    }

  // checking NormInf
  scal_ex = complex<T>(0, 0);
  for (int i = 0; i < m; i++)
    scal_ex = complex<T>(max(abs(scal_ex), abs(x(i))), 0);
  
  scal = complex<T>(NormInf(x), 0);
  if ((abs(scal - scal_ex) > threshold) || isnan(abs(scal - scal_ex)))
    {
      cout << "NormInf incorrect" << endl;
      abort();
    }  
}

template<class T, int m>
void CheckScalarProductTinyVector(TinyVector<T, m>& x, TinyVector<T, m>& y)
{
  T scal = DotProd(x, y);
  T scal_ex = 0;
  for (int i = 0; i < m; i++)
    scal_ex += x(i)*y(i);
  
  if ((abs(scal - scal_ex) > threshold) || isnan(abs(scal - scal_ex)))
    {
      cout << "DotProd incorrect" << endl;
      abort();
    }
  
  // checking Norm2, Distance and DistanceSquare
  T module = sqrt(abs(DotProd(x, x)));
  if ((abs(module - Norm2(x)) > threshold) || isnan(abs(module - Norm2(x))))
    {
      cout << "Norm2 incorrect " << endl;
      abort();
    }

  module = abs(DotProd(x, x));
  if ((abs(module - AbsSquare(x)) > threshold) || isnan(abs(module - AbsSquare(x))))
    {
      cout << "AbsSquare incorrect " << endl;
      abort();
    }
  
  TinyVector<T, m> diff;
  diff = x-y;
  module = Norm2(diff);
  if ((abs(module - x.Distance(y)) > threshold) || isnan(abs(module - x.Distance(y))))
    {
      cout << "Distance incorrect " << endl;
      abort();
    }

  if ((abs(module - Distance(x, y)) > threshold) || isnan(abs(module - Distance(x, y))))
    {
      cout << "Distance incorrect " << endl;
      abort();
    }

  if ((abs(module*module - x.DistanceSquare(y)) > threshold) || isnan(abs(module*module - x.DistanceSquare(y))))
    {
      cout << "DistanceSquare incorrect " << endl;
      abort();
    }

  // checking Norm1
  scal_ex = 0;
  for (int i = 0; i < m; i++)
    scal_ex += abs(x(i));
  
  scal = Norm1(x);
  if ((abs(scal - scal_ex) > threshold) || isnan(abs(scal - scal_ex)))
    {
      cout << "Norm1 incorrect" << endl;
      abort();
    }

  // checking NormInf
  scal_ex = 0;
  for (int i = 0; i < m; i++)
    scal_ex = max(scal_ex, abs(x(i)));
  
  scal = NormInf(x);
  if ((abs(scal - scal_ex) > threshold) || isnan(abs(scal - scal_ex)))
    {
      cout << "NormInf incorrect" << endl;
      abort();
    }

  // checking UpdateMinimum/UpdateMaximum
  TinyVector<T, m> u, v, w, z, zmax, zmin;

  u.FillRand(); u *= 1e-9;
  v.FillRand(); v *= 1e-9;
  w.FillRand(); w *= 1e-9;
  z.FillRand(); z *= 1e-9;

  zmin = u;
  zmax = u;
  
  UpdateMinimum(v, zmin);
  UpdateMinimum(w, zmin);
  UpdateMinimum(z, zmin);

  UpdateMaximum(v, zmax);
  UpdateMaximum(w, zmax);
  UpdateMaximum(z, zmax);

  for (int i = 0; i < m; i++)
    {
      scal_ex = min(u(i), v(i));
      scal_ex = min(scal_ex, w(i));
      scal_ex = min(scal_ex, z(i));

      if (isnan(abs(scal_ex - zmin(i))) || (abs(zmin(i) - scal_ex) > threshold))
	{
	  cout << "UpdateMinimum incorrect" << endl;
	  abort();
	}

      scal_ex = max(u(i), v(i));
      scal_ex = max(scal_ex, w(i));
      scal_ex = max(scal_ex, z(i));

      if (isnan(abs(scal_ex - zmax(i))) || (abs(zmax(i) - scal_ex) > threshold))
	{
	  cout << "UpdateMaximum incorrect" << endl;
	  abort();
	}
    }
}

template<class T, int m>
void CheckTinyVector(TinyVector<T, m>& x)
{
  TinyVector<T, m>::threshold = threshold;

  // checking operator =, and copy constructor
  x.FillRand(); Mlt(Real_wp(1e-8), x);
  TinyVector<T, m> y(x), z = x, u, v, w, zero;
  for (int i = 0; i < m; i++)
    if ((abs(x(i) - y(i)) > threshold) || (abs(x(i) - z(i)) > threshold) ||
        isnan(abs(x(i) - y(i))) || isnan(abs(x(i) - z(i))))      
      {
        cout << "Operator = incorrect" << endl;
        abort();
      }
  
  // checking comparison operator
  for (int i = 0; i < m; i++)
    y(i) += 0.5*threshold*double(rand())/RAND_MAX;
  
  if (x == y)
    {
    }
  else
    {
      cout << "Operator == incorrect" << endl;
      abort();
    }
  
  if (x != y)
    {
      cout << "Operator != incorrect" << endl;
      abort();
    }
  
  y(1) += 0.1;
  if (x == y)
    {
      cout << "Operator == incorrect" << endl;
      abort();
    }

  if (x != y)
    {
    }
  else
    {
      cout << "Operator != incorrect" << endl;
      abort();
    }
  
  // checking other operators for real numbers only
  CheckComparisonTinyVector(x);
  
  // checking operator +=
  x.FillRand(); y.FillRand();
  x *= 1e-8; y *= 1e-8;
  z = x;
  z += y;
  for (int i = 0; i < m; i++)
    if ((abs(z(i) - x(i) - y(i)) > threshold) || isnan(abs(z(i) - x(i) - y(i))))
      {
        cout << "Operator += incorrect " << endl;
        abort();
      }
  
  // checking operator -=
  z = x;
  z -= y;
  for (int i = 0; i < m; i++)
    if ((abs(z(i) - x(i) + y(i)) > threshold) || isnan(abs(z(i) - x(i) + y(i))))
      {
        cout << "Operator -= incorrect " << endl;
        abort();
      }
  
  // checking operator *=
  z = x;
  z = z*y;
  for (int i = 0; i < m; i++)
    if ((abs(z(i) - x(i)*y(i)) > threshold) || isnan(abs(z(i) - x(i)*y(i))))
      {
        cout << "Operator *= incorrect " << endl;
        abort();
      }
  
  u = z;
  u *= 2.5;
  for (int i = 0; i < m; i++)
    if ((abs(u(i) - 2.5*z(i)) > threshold) || isnan(abs(u(i) - 2.5*z(i))))
      {
        cout << "Operator *= incorrect " << endl;
        abort();
      }

  // checking operator /=
  z = x;
  z = z/y;
  for (int i = 0; i < m; i++)
    if ((abs(z(i) - x(i)/y(i)) > threshold) || isnan(abs(z(i) - x(i)/y(i))))
      {
        cout << "Operator /= incorrect " << endl;
        abort();
      }

  u = z;
  u /= Real_wp(2.5);
  for (int i = 0; i < m; i++)
    if ((abs(u(i) - z(i)/Real_wp(2.5)) > threshold) || isnan(abs(u(i) - z(i)/Real_wp(2.5))))
      {
        DISP(i); DISP(u(i)); DISP(z(i)/Real_wp(2.5)); DISP(abs(u(i)-z(i)/Real_wp(2.5)));
        cout << "Operator /= incorrect " << endl;
        abort();
      }
  
  // checking operator -
  u = -y;
  u += y;
  if (u != zero)
    {
      cout << "Operator - incorrect " << endl;
      abort();
    }

  u = x-y;
  v = x; v -= y;
  if (u != v)
    {
      cout << "Operator - incorrect " << endl;
      abort();
    }
  
  // checking operator +
  u = x+y;
  v = x; v += y;
  if (u != v)
    {
      cout << "Operator + incorrect " << endl;
      abort();
    }

  // checking operator *
  u = x*y;
  v = x; v = v*y;
  if (u != v)
    {
      cout << "Operator * incorrect " << endl;
      abort();
    }
  
  v = 2.5*x;
  u = x; u *= 2.5;
  if (u != v)
    {
      cout << "Operator * incorrect " << endl;
      abort();
    }
  
  v = x*2.5;
  u = x; u *= 2.5;
  if (u != v)
    {
      cout << "Operator * incorrect " << endl;
      abort();
    }

  // checking operator /
  u = x/y;
  v = x; v = v/y;
  if (u != v)
    {
      cout << "Operator / incorrect " << endl;
      abort();
    }
  
  v = 2.5/x;
  for (int i = 0; i < m; i++)
    if ((abs(v(i) - 2.5/x(i)) > threshold) || isnan(abs(v(i) - 2.5/x(i))))
      {
        cout << "Operator / incorrect " << endl;
        abort();
      }
  
  v = x/2.5;
  u = x; u /= 2.5;
  if (u != v)
    {
      cout << "Operator / incorrect " << endl;
      abort();
    }
  
  // checking GetM, GetLength, GetSize
  if ((w.GetM() != m) || (w.GetLength() != m) || (w.GetSize() != m))
    {
      cout << "GetM incorrect" << endl;
      abort();
    }
  
  // checking Zero
  u.Zero();
  if (u != zero)
    {
      cout << "Zero incorrect" << endl;
      abort();
    }
  
  // checking Fill
  u.Fill();
  for (int i = 0; i < m; i++)
    if ((abs(u(i) - i) > threshold) || isnan(abs(u(i) - i)))
      {
        cout << "Fill incorrect" << endl;
        abort();
      }

  u.Fill(3.4);
  for (int i = 0; i < m; i++)
    if ((abs(u(i) - 3.4) > threshold) || isnan(abs(u(i) - 3.4)))
      {
        cout << "Fill incorrect" << endl;
        abort();
      }
  
  // checking DotProd
  CheckScalarProductTinyVector(x, y);
  
  // checking CopyVector
  TinyVector<Vector<T>, m> V;
  for (int j = 0; j < m; j++)
    {
      V(j).Reallocate(10);
      V(j).FillRand();
      Mlt(Real_wp(1e-8), V(j));
    }
  
  for (int j = 0; j < 10; j++)
    {
      CopyVector(V, j, u);
      for (int k = 0; k < m; k++)
        if ((abs(u(k) - V(k)(j)) > threshold) || isnan(abs(u(k) - V(k)(j))))
          {
            cout << "CopyVector incorrect" << endl;
            abort();
          }
      
      u *= 1.3;
      CopyVector(u, j, V);
      for (int k = 0; k < m; k++)
        if ((abs(u(k) - V(k)(j)) > threshold) || isnan(abs(u(k) - V(k)(j))))
          {
            cout << "CopyVector incorrect" << endl;
            abort();
          }
    }
  
  Vector<T> W(10*m);
  W.FillRand(); Mlt(Real_wp(1e-8), W);
  for (int j = 0; j < 10; j++)
    {
      CopyVector(W, j, u);
      for (int k = 0; k < m; k++)
        if ((abs(u(k) - W(j*m + k)) > threshold) || isnan(abs(u(k) - W(j*m + k))))
          {
            cout << "CopyVector incorrect" << endl;
            abort();
          }
      
      u *= 1.3;
      CopyVector(u, j, W);
      for (int k = 0; k < m; k++)
        if ((abs(u(k) - W(j*m + k)) > threshold) || isnan(abs(u(k) - W(j*m + k))))
          {
            cout << "CopyVector incorrect" << endl;
            abort();
          }
    }
  
  // checking Add
  u = x+y;
  if (u != x+y)
    {
      cout << "Add incorrect"  << endl;
      abort();
    }
  
  v = u; u += x;
  v = v+x;
  if (v != u)
    {
      cout << "Add incorrect"  << endl;
      abort();
    }
  
  v = u; v += 2.5*y;
  u += 2.5*y;
  if (u != v)
    {
      cout << "Add incorrect"  << endl;
      abort();
    }
  
  // checking Subtract
  u = x-y;
  v = x-y;
  if (u != v)
    {
      cout << "Subtract incorrect"  << endl;
      abort();
    }
  
  // checking Mlt
  u = x; v = x;
  u *= 2.2; 
  Mlt(Real_wp(2.2), v);
  if (u != v)
    {
      cout << "Mlt incorrect"  << endl;
      abort();
    }
  
  // checking MltVector
  for (int j = 0; j < 10; j++)
    {
      CopyVector(V, j, u);
      T coef;
      SetComplexReal(double(j+1)/20, coef);
      u *= coef;
      MltVector(coef, j, V);
      
      for (int k = 0; k < m; k++)
        if ((abs(u(k) - V(k)(j)) > threshold) || isnan(abs(u(k) - V(k)(j))))
          {
            cout << "MltVector incorrect" << endl;
            abort();
          }
    }

  // checking SymmetrizePointPlane
  TinyVector<T, m> normale;
  normale.FillRand(); Mlt(1.0/Norm2(normale), normale);
  T d = rand() / Real_wp(RAND_MAX);

  u.FillRand(); u *= 1e-9;

  SymmetrizePointPlane(u, normale, d, v);

  x = Real_wp(0.5)*(u+v);

  T vloc = DotProd(normale, x) + d;
  if (isnan(abs(vloc)) || (abs(vloc) > threshold))
    {
      cout << "SymmetrizePointPlane incorrect" << endl;
      abort();
    }

  TinyVector<T, m> u1, u2, v1, v2;
  u1 = u - x;
  u2 = v - x;

  T alpha1 = DotProd(u1, normale);
  T alpha2 = DotProd(u2, normale);

  if (abs(alpha1+alpha2) > threshold)
    {
      cout << "SymmetrizePointPlane incorrect" << endl;
      abort();
    }

  v1 = alpha1*normale;
  v2 = alpha2*normale;

  if ((!EqualVector(u1, v1)) || (!EqualVector(u2, v2)))
    {
      cout << "SymmetrizePointPlane incorrect" << endl;
      abort();
    }
    
}

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  bool overall_success = true;
  threshold = 1e4*epsilon_machine;
  
  {
    // testing contructors, init, and operator()
    TinyVector<Real_wp, 1> y(2.3);
    if (abs(y(0) - 2.3) > threshold)
      {
        cout << "Constructor incorrect " << endl;
        abort();
      }
    
    y.Init(0.5);
    if (abs(y(0) - 0.5) > threshold)
      {
        cout << "Constructor incorrect " << endl;
        abort();
      }
    
    TinyVector<Real_wp, 3> x(-0.7, 1.2, 3.7), x2(3.4);
    if ( (abs(x(0) + 0.7) > threshold) || (abs(x(1) - 1.2) > threshold)
         || (abs(x(2) - 3.7) > threshold))
      {
        cout << "Constructor incorrect " << endl;
        abort();
      }
    
    if ((abs(x2(0) - 3.4) > threshold) || (abs(x2(1)-3.4) > threshold) 
        || (abs(x2(2)-3.4) > threshold))
      {
        cout << "Constructor incorrect " << endl;
        abort();
      }
    
    // you can change only first components with Init    
    x.Init(0.3, 4.1);
    if ( (abs(x(0) - 0.3) > threshold) || (abs(x(1) - 4.1) > threshold) 
      || (abs(x(2) - 3.7) > threshold))
      {
        cout << "Constructor incorrect " << endl;
        abort();
      }
    
    // you can use constructor up to 5 arguments
    TinyVector<Real_wp, 5> z(0.6, -1.4, 0.9, 2.1, 4.5);
    if ( abs( (z(0) - 0.6)*(z(1) + 1.4)*(z(2) - 0.9)
               *(z(3) - 2.1)*(z(4) - 4.5)) > threshold)
      {
        cout << "Constructor incorrect " << endl;
        abort();
      }

    z.Init(0.8, 1.9, -1.2, 2.5, -3.2);
    if ( abs( (z(0) - 0.8)*(z(1) - 1.9)*(z(2) + 1.2)
               *(z(3) - 2.5)*(z(4) + 3.2)) > threshold)
      {
        cout << "Init incorrect " << endl;
        abort();
      }

  }

  {
    TinyVector<Real_wp, 3> x, u, v, w;
    CheckTinyVector(x);
    
    // functions specific to points in R3
    w.Init(0, Real_wp(3)/5, Real_wp(4)/5);
    GetVectorPlane(w, u, v);
    if ( (abs(Norm2(u) -1.0) <= threshold) && (abs(Norm2(v) -1.0) <= threshold)
         && (abs(DotProd(w, u)) <= threshold) && (abs(DotProd(w, v)) <= threshold)
         && (abs(DotProd(u, v)) <= threshold) )
      {
        // correct
      }
    else
      {
        cout << "GetVectorPlane incorrect" << endl;
        abort();
      }
    
    if (abs(Determinant(w, u, v) - 1.0) > threshold)
      {
        DISP(w); DISP(u); DISP(v); DISP(Determinant(w, u, v) - Real_wp(1));
        cout << "Determinant not correct" << endl;
        abort();
      }
    
    TimesProd(u, v, x);
    if (x != w)
      {
        cout << "TimesProd incorrect " << endl;
        abort();
      }

    TinyVector<Real_wp, 3> normale;
    u.FillRand(); u *= Real_wp(1e-9);
    normale.FillRand(); Mlt(1.0/Norm2(normale), normale);
    ForceZeroVdotN(normale, u);

    if (isnan(DotProd(normale, u)) || (abs(DotProd(normale, u)) > threshold))
      {
	cout << "ForceZeroVdotN incorrect" << endl;
	abort();
      }

    TinyVector<TinyVector<Real_wp, 3>, 2> box;
    box(0).FillRand(); box(0) *= 1e-9;
    box(1).FillRand(); box(1) *= 1e-9;
    u = box(0);
    UpdateMinimum(box(1), box(0));
    UpdateMaximum(u, box(1));

    DISP(box(0)); DISP(box(1));

    for (int i = 0; i < 3; i++)
      u(i) = box(0)(i) + (box(1)(i) - box(0)(i))*Real_wp(rand())/RAND_MAX;

    DISP(u);
    if (!PointInsideBoundingBox(u, box))
      {
	cout << "PointInsideBoundingBox incorrect" << endl;
	abort();
      }

    for (int i = 0; i < 3; i++)
      {
	v = u;
	v(i) += box(1)(i)-box(0)(i);
	if (PointInsideBoundingBox(v, box))
	  {
	    cout << "PointInsideBoundingBox incorrect" << endl;
	    abort();
	  }

	v = u;
	v(i) -= box(1)(i)-box(0)(i);
	if (PointInsideBoundingBox(v, box))
	  {
	    cout << "PointInsideBoundingBox incorrect" << endl;
	    abort();
	  }
      }
    
  }
  
  {
    TinyVector<Complex_wp, 3> x;
    CheckTinyVector(x);
  }

  {
    TinyVector<Real_wp, 5> x;
    CheckTinyVector(x);
  }

  {
    TinyVector<Complex_wp, 5> x;
    CheckTinyVector(x);
  }

  if (overall_success)
    cout << "All tests passed successfully " << endl;
  
  return FinalizeMontjoie();
}  

