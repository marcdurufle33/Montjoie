#ifndef MONTJOIE_FILE_COMMON_INLINE_CXX

namespace std 
{
  //! returns the square of x
  template<class T>
  inline T square(const T& x)
  {
    return x*x;
  }
}

namespace Seldon
{

  //! sorts a and b
  inline void Sort (int& a, int& b)
  {
    if (b < a) 
      {  
	int temp = a;
	a = b;
	b = temp;
      }
  }
  
  
  //! sorts a and b
  template<class T, int m>
  inline void Sort (TinyVector<T, m>& a, TinyVector<T, m>& b)
  {
    if (b < a) 
      {  
	TinyVector<T, m> temp = a;
	a = b;
	b = temp;
      }
  }
  
  //! sorts a, b and c
  inline void Sort (int& a, int& b, int& c)
  {
    if (b < a) 
      {  
	int temp = a;
	a = b;
	b = temp;
      }
    if (c < a) 
      {  
	int temp = a;
	a = c;
	int temp2 = b;
	b = temp;
	c = temp2;
      }
    else if (c < b)
      {
	int temp = c;
	c = b;
	b = temp;
      }
  }
  
  
  //! Sorts i, j, k and l
  inline void Sort(int& i, int& j, int& k, int& l)
  {
    int i0 = i,j0 = j,k0 = k,l0 = l;
    if (i > j)
      {
	i0 = j;
	j0 = i;
      }
    
    if (k > l)
      {
	k0 = l;
	l0 = k;
      }
    
    if (i0 < k0)
      {
	i = i0;
	if (j0 < l0) 
	  {
	    l = l0;
	    if (j0 < k0)
	      {
		j = j0;
		k = k0;
	      }
	    else
	      {
		j = k0;
		k = j0;
	      }
	  }
	else
	  {
	    j = k0;
	    k = l0;
	    l = j0;
	  }
      }
    else
      {
	i = k0;
	if (l0 < j0) 
	  {
	    l = j0;
	    if (l0 < i0)
	      {
		j = l0;
		k = i0;
	      }
	    else
	      {
		j = i0;
		k = l0;
	      }
	  }
	else
	  {
	    j = i0;
	    k = j0;
	    l = l0;
	  }
      }  
  }

  //! Checks if X is large enough and resize if necessary
  template<class T, class Storage1, class Allocator1>
  inline void Check_And_ReallocateIfNecessary(Vector<T, Storage1, Allocator1>& X, int n)
  {
    if (n >= X.GetM())
      X.Resize(int(ceil(1.5*n+2)));
  }
  
  
  //! calculates wave vector from frequency and incident angle
  template<class T>
  inline void SetIncidentAngle(const T& omega, TinyVector<T, 2>& kwave,
			       const T& teta, const T& phi)
  {
    kwave(0) = omega*cos(teta);
    kwave(1) = omega*sin(teta);
  }
  
  //! calculates wave vector from frequency and incident angles
  template<class T>
  inline void SetIncidentAngle(const T& omega, TinyVector<T, 3>& kwave,
                               const T& teta, const T& phi)
  {
    kwave(0) = omega*sin(teta)*cos(phi);
    kwave(1) = omega*sin(teta)*sin(phi);
    kwave(2) = omega*cos(teta);
  }


  //! prints the vector u without brackets for inner vectors
  template<class T>
  inline void PrintNoBrackets(ostream& out, const Vector<T>& u)
  {
    for (int i = 0; i < u.GetM()-1; i++)
      {
        PrintNoBrackets(out, u(i));
        out << " ";
      }
    
    if (u.GetM() > 0)
      PrintNoBrackets(out, u(u.GetM()-1));    
  }

  
  //! prints a scalar
  template<class T>
  inline void PrintNoBrackets(ostream& out, const T& u)
  {
    out << u;
  }


  //! prints a tiny vector without brackets (only spaces)
  template<class T, int m>
  inline void PrintNoBrackets(ostream& out, const TinyVector<T, m>& u)
  {
    u.Print(out);
  }


  template<class T>
  inline size_t GetMemorySize(const Vector<Vector<T> >& x)
  {
    size_t taille = sizeof(Vector<T>) + sizeof(Vector<T>*)*x.GetM();
    for (int i = 0; i < x.GetM(); i++)
      taille += x(i).GetMemorySize();
    
    return taille;
  }


  template<class T>
  inline size_t GetMemorySize(const Vector<Vector<T, VectSparse> >& x)
  {
    size_t taille = sizeof(Vector<T>) + sizeof(Vector<T>*)*x.GetM();
    for (int i = 0; i < x.GetM(); i++)
      taille += x(i).GetMemorySize();
    
    return taille;
  }


  template<class T, class Prop, class Storage>
  inline size_t GetMemorySize(const Vector<Matrix<T, Prop, Storage> >& x)
  {
    size_t taille = sizeof(Vector<T>) + sizeof(void*)*x.GetM();
    for (int i = 0; i < x.GetM(); i++)
      taille += x(i).GetMemorySize();
    
    return taille;
  }


  template<class T>
  inline size_t GetMemorySize(const Array2D<Matrix<T> >& x)
  {
    size_t taille = sizeof(Matrix<T>) + sizeof(Matrix<T>*)*x.GetM();
    for (int i = 0; i < x.GetM(); i++)
      for (int j = 0; j < x.GetN(); j++)
        taille += x(i, j).GetMemorySize();
    
    return taille;
  }
  
  
  template<class T>
  inline size_t GetMemorySize(const Vector<Vector<Vector<T> > >& x)
  {
    size_t taille = sizeof(Vector<T>) + sizeof(Vector<T>*)*x.GetM();
    for (int i = 0; i < x.GetM(); i++)
      {
        taille += sizeof(Vector<T>) + sizeof(Vector<T>*)*x(i).GetM();
        for (int j = 0; j < x(i).GetM(); j++)
          taille += x(i)(j).GetMemorySize();
      }
    
    return taille;
  }


  template<class T, class Prop, class Storage>
  inline size_t GetMemorySize(const Vector<Vector<Matrix<T, Prop, Storage> > >& x)
  {
    size_t taille = sizeof(Vector<T>) + sizeof(Vector<T>*)*x.GetM();
    for (int i = 0; i < x.GetM(); i++)
      {
        taille += sizeof(Vector<T>) + sizeof(Vector<T>*)*x(i).GetM();
        for (int j = 0; j < x(i).GetM(); j++)
          taille += x(i)(j).GetMemorySize();
      }
    
    return taille;
  }
  
  
  template<class T, int q>
  inline size_t GetMemorySize(const Vector<TinyVector<Vector<T>, q> >& x)
  {
    size_t taille = sizeof(Vector<T>) + sizeof(Vector<T>*)*x.GetM();
    for (int i = 0; i < x.GetM(); i++)
      {
        taille += sizeof(Vector<T>*)*x(i).GetM();
        for (int j = 0; j < x(i).GetM(); j++)
          taille += x(i)(j).GetMemorySize();
      }
    
    return taille;

  }

  template<class T, int q>
  inline size_t GetMemorySize(const TinyVector<Vector<T>, q>& x)
  {
    size_t taille = sizeof(Vector<T>*)*x.GetM();
    for (int j = 0; j < x.GetM(); j++)
      taille += x(j).GetMemorySize();
    
    return taille;
  }
  
}

#define MONTJOIE_FILE_COMMON_INLINE_CXX
#endif

