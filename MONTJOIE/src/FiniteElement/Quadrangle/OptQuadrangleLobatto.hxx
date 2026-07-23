#ifndef MONTJOIE_FILE_OPT_SQUARE_LOBATTO_HXX

namespace Montjoie
{
  template<int r>
  class OptQuadrangleLobatto
  {
  public :
    enum { nb_dof = (r+1)*(r+1), nb_quad = (r+1)*(r+1) };

    static TinyVector<Real_wp, (r+1)*(r+1)> Weights_ND;
    static TinyMatrix<Real_wp, General, r+1, r+1> dphi;
    
    static void ConstructFiniteElement(const QuadrangleLobatto& Fb);
    
    template<class Vector1, class Vector2>
    static void ApplyCh(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void ApplyCh_Transpose(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void ApplyRh(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void ApplyRh_Transpose(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void MltMassMatrix(const Vector1&, Vector2&);
    
  };

#ifndef SELDON_HEADER_ONLY
  template<int r>
  TinyVector<Real_wp, (r+1)*(r+1)> OptQuadrangleLobatto<r>::Weights_ND;
  template<int r>
  TinyMatrix<Real_wp, General, r+1, r+1> OptQuadrangleLobatto<r>::dphi;  
#endif
  
  template<int i, int j, int r>
  class SquareNumNodes2D
  {
  public :
    enum
      {
	num = (i==r)*(j==0) + 2*(i==r)*(j==r) + 3*(i==0)*(j==r)
	+ (3+i)*(i >= 1)*(i < r)*(j==0) + (3+(r-1)+j)*(j >= 1)*(j < r)*(i==r)
	+ (3+2*(r-1)+r-i)*(i >= 1)*(i < r)*(j==r) + (3+3*(r-1)+r-j)*(j >= 1)*(j < r)*(i==0)
	+ (4*r + (i-1)*(r-1) + j-1)*(i >= 1)*(i < r)*(j >= 1)*(j < r)
      };
  };

  
  template<int i, int j, int k, int r>
  class QuadLobattoLoop3
  {
  public :
    template<class Vector1, class Vector2>
    static void EvalV(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void EvalU(const Vector1& Vh, Vector2& Uh);
    
  };
  
  template<int i, int j, int r>
  class QuadLobattoLoop3<-1, i, j, r>
  {
  public :
    template<class Vector1, class Vector2>
    static inline void EvalV(const Vector1& Uh, Vector2& Vh)
    {  }
    
    template<class Vector1, class Vector2>
    static inline void EvalU(const Vector1& Vh, Vector2& Uh)
    {  }
    
  };
  
  template<int i, int j, int r>
  class QuadLobattoLoop2
  {
  public :
    template<class Vector1, class Vector2>
    static void EvalV(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void EvalU(const Vector1& Vh, Vector2& Uh);
    
  };
  
  template<int i, int r>
  class QuadLobattoLoop2<i, -1, r>
  {
  public :
    template<class Vector1, class Vector2>
    static inline void EvalV(const Vector1& Uh, Vector2& Vh)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void EvalU(const Vector1& Vh, Vector2& Uh)
    {
    }
    
  };


  template<int i, int r>
  class QuadLobattoLoop1
  {
  public :
    template<class Vector1, class Vector2>
    static void EvalV(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void EvalU(const Vector1& Vh, Vector2& Uh);
    
  };
  
  template<int r>
  class QuadLobattoLoop1<-1, r>
  {
  public :
    template<class Vector1, class Vector2>
    static inline void EvalV(const Vector1& Uh, Vector2& Vh)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void EvalU(const Vector1& Vh, Vector2& Uh)
    {
    }
        
  };

  template<int r>
  class ImplementationOptimOrder<QuadrangleLobatto, r>
  {
  public : 
    static const bool implemented = (r <= 8);
  };
    
}

#define MONTJOIE_FILE_OPT_SQUARE_LOBATTO_HXX
#endif

