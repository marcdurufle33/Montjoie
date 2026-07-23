#ifndef MONTJOIE_FILE_OPT_HEXAHEDRA_LOBATTO_HXX

namespace Montjoie
{
  template<int r>
  class OptHexahedronLobatto
  {
  public :
    enum { nb_dof = (r+1)*(r+1)*(r+1), nb_quad = (r+1)*(r+1)*(r+1) };

    static TinyVector<Real_wp, (r+1)*(r+1)*(r+1)> WeightsND;
    static TinyMatrix<Real_wp, General, r+1, r+1> dphi;
    
    static void ConstructFiniteElement(const HexahedronLobatto& Fb);
    
    template<class Vector1, class Vector2>
    static void ApplyCh(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void ApplyChTranspose(const Vector1& Uh, Vector2& Vh);
    
    template<class T>
    static void ApplyRh(const T* Uh, T* Vh);

    template<class T>
    static void ApplyRh(const Vector<T>& Uh, Vector<T>& Vh);

    template<class Vector1, class Vector2>
    static void ApplyRhSplit(const Vector1& Uh, Vector2&, Vector2&, Vector2&);
    
    template<class T>
    static void ApplyRhTranspose(const T* Uh, T* Vh);

    template<class T>
    static void ApplyRhTranspose(const Vector<T>& Uh, Vector<T>& Vh);

    template<class Vector1, class Vector2>
    static void MltMassMatrix(const Vector1&, Vector2&);
    
  };

#ifndef SELDON_HEADER_ONLY
  template<int r>
  TinyVector<Real_wp, (r+1)*(r+1)*(r+1)> OptHexahedronLobatto<r>::WeightsND;
  template<int r>
  TinyMatrix<Real_wp, General, r+1, r+1> OptHexahedronLobatto<r>::dphi;  
#endif
  
  template<int i, int j, int k, int r>
  class CubeNumNodes3D
  {
  public :
    enum
      {
	num = (6*r*r+2+(i-1)*(r-1)*(r-1) + (j-1)*(r-1) +k-1)*(i>0)*(j>0)*(k>0)*(i<r)*(j<r)*(k<r)
	+ (j==0)*(k==0)*(7+i) + (i==r)*(k==0)*(7+(r-1)+j) + (j==r)*(k==0)*(7+2*(r-1)+i)
        + (i==0)*(k==0)*(7+3*(r-1)+j)
	+ (i==0)*(j==0)*(7+4*(r-1)+k) + (i==r)*(j==0)*(7+5*(r-1)+k)
        + (i==r)*(j==r)*(7+6*(r-1)+k) + (i==0)*(j==r)*(7+7*(r-1)+k)
	+ (j==0)*(k==r)*(7+8*(r-1)+i) + (i==r)*(k==r)*(7+9*(r-1)+j)
        + (j==r)*(k==r)*(7+10*(r-1)+i) + (i==0)*(k==r)*(7+11*(r-1)+j)
	+ (i==0)*(j>0)*(j<r)*(k>0)*(k<r)*(12*r-5+(j-1)*(r-1)+k)
        + (i==r)*(j>0)*(j<r)*(k>0)*(k<r)*(12*r-5+5*(r-1)*(r-1)+(j-1)*(r-1)+k)
	+ (j==0)*(i>0)*(i<r)*(k>0)*(k<r)*(12*r-5+(r-1)*(r-1)+(i-1)*(r-1)+k)
        + (j==r)*(i>0)*(i<r)*(k>0)*(k<r)*(12*r-5+4*(r-1)*(r-1)+(i-1)*(r-1)+k)
	+ (k==0)*(i>0)*(i<r)*(j>0)*(j<r)*(12*r-5+2*(r-1)*(r-1)+(i-1)*(r-1)+j)
	+ (k==r)*(i>0)*(i<r)*(j>0)*(j<r)*(12*r-5+3*(r-1)*(r-1)+(i-1)*(r-1)+j)
      };
  };
  
  template<int r>
  class CubeNumNodes3D<0, 0, 0, r>
  {
  public :
    enum { num = 0};
  };
  
  template<int r>
  class CubeNumNodes3D<r, 0, 0, r>
  {
  public :
    enum { num = 1};
  };
  
  template<int r>
  class CubeNumNodes3D<r, r, 0, r>
  {
  public :
    enum { num = 2};
  };
  
  template<int r>
  class CubeNumNodes3D<0, r, 0, r>
  {
  public :
    enum { num = 3};
  };
  
  template<int r>
  class CubeNumNodes3D<0, 0, r, r>
  {
  public :
    enum { num = 4};
  };
  
  template<int r>
  class CubeNumNodes3D<r, 0, r, r>
  {
  public :
    enum { num = 5};
  };
  
  template<int r>
  class CubeNumNodes3D<r, r, r, r>
  {
  public :
    enum { num = 6};
  };
  
  template<int r>
  class CubeNumNodes3D<0, r, r, r>
  {
  public :
    enum { num = 7};
  };

  template<int m, int i, int j, int k, int r>
  class HexLobattoLoop4
  {
  public :
    template<class T>
    static void EvalV(const T* Uh, T* Vh);
    
    template<class T>
    static void EvalU(const T* Vh, T* Uh);

    template<class Vector1, class Vector2>
    static void EvalUsplit(const Vector1&, Vector2&, Vector2&, Vector2&);
    
  };
  
  template<int i, int j, int k, int r>
  class HexLobattoLoop4<-1, i, j, k, r>
  {
  public :
    template<class T>
    static inline void EvalV(const T* Uh, T* Vh)
    {
    }
    
    template<class T>
    static inline void EvalU(const T* Vh, T* Uh)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void EvalUsplit(const Vector1&, Vector2&, Vector2&, Vector2&) {}
    
  };
  
  template<int i, int j, int k, int r>
  class HexLobattoLoop3
  {
  public :
    template<class T>
    static void EvalV(const T* Uh, T*);
    
    template<class T>
    static void EvalU(const T* Vh, T* Uh);
    
    template<class Vector1, class Vector2>
    static void EvalUsplit(const Vector1&, Vector2&, Vector2&, Vector2&);
    
  };
  
  template<int i, int j, int r>
  class HexLobattoLoop3<i, j, -1, r>
  {
  public :
    template<class T>
    static inline void EvalV(const T* Uh, T* Vh)
    {  }
    
    template<class T>
    static inline void EvalU(const T* Vh, T* Uh)
    {  }
    
    template<class Vector1, class Vector2>
    static inline void EvalUsplit(const Vector1&, Vector2&, Vector2&, Vector2&) {}
    
  };
  
  template<int i, int j, int r>
  class HexLobattoLoop2
  {
  public :
    template<class T>
    static void EvalV(const T* Uh, T*);
    
    template<class T>
    static void EvalU(const T* Vh, T* Uh);

    template<class Vector1, class Vector2>
    static void EvalUsplit(const Vector1&, Vector2&, Vector2&, Vector2&);
    
  };
  
  template<int i, int r>
  class HexLobattoLoop2<i, -1, r>
  {
  public :
    template<class T>
    static inline void EvalV(const T* Uh, T* Vh)
    {  }
    
    template<class T>
    static inline void EvalU(const T* Vh, T* Uh)
    {  }

    template<class Vector1, class Vector2>
    static inline void EvalUsplit(const Vector1&, Vector2&, Vector2&, Vector2&) {}
    
  };


  template<int i, int r>
  class HexLobattoLoop1
  {
  public :
    template<class T>
    static void EvalV(const T* Uh, T*);
    
    template<class T>
    static void EvalU(const T* Vh, T* Uh);
    
    template<class Vector1, class Vector2>
    static void EvalUsplit(const Vector1&, Vector2&, Vector2&, Vector2&);
    
  };
  
  template<int r>
  class HexLobattoLoop1<-1, r>
  {
  public :
    template<class T>
    static inline void EvalV(const T* Uh, T* Vh)
    {  }
    
    template<class T>
    static inline void EvalU(const T* Vh, T* Uh)
    {  }

    template<class Vector1, class Vector2>
    static inline void EvalUsplit(const Vector1&, Vector2&, Vector2&, Vector2&) {}
        
  };
    
}

#define MONTJOIE_FILE_OPT_HEXAHEDRA_LOBATTO_HXX
#endif

