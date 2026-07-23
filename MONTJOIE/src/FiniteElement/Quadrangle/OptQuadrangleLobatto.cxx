#ifndef MONTJOIE_FILE_OPT_SQUARE_LOBATTO_CXX

namespace Montjoie
{
  
  template<int r> template<class Vector1, class Vector2>
  inline void OptQuadrangleLobatto<r>::ApplyCh(const Vector1& Uh, Vector2& Vh)
  {
    Vh = Uh;
  }
    
  
  template<int r> template<class Vector1, class Vector2>
  inline void OptQuadrangleLobatto<r>::ApplyCh_Transpose(const Vector1& Uh, Vector2& Vh)
  {
    Vh = Uh;
  }
  
  
  template<int r> template<class Vector1, class Vector2>
  inline void OptQuadrangleLobatto<r>::ApplyRh(const Vector1& Vh, Vector2& Uh)
  {
    QuadLobattoLoop1<r, r>::EvalU(Vh, Uh);
  }
    
  
  template<int r> template<class Vector1, class Vector2>
  inline void OptQuadrangleLobatto<r>::ApplyRh_Transpose(const Vector1& Uh, Vector2& Vh)
  {
    QuadLobattoLoop1<r, r>::EvalV(Uh, Vh);
  }
    
  
  template<int r>  template<class Vector1, class Vector2>
  inline void OptQuadrangleLobatto<r>::MltMassMatrix(const Vector1& Uh, Vector2& Vh)
  {
    MltDot(Uh, Weights_ND, Vh);
  }
  
  
  template<int m, int i, int j, int r> template<class Vector1, class Vector2>
  inline void QuadLobattoLoop3<m, i, j, r>::EvalV(const Vector1& Uh, Vector2& Vh)
  {
    Vh(2*SquareNumNodes2D<i, j, r>::num)
      += OptQuadrangleLobatto<r>::dphi(m, i)*Uh(SquareNumNodes2D<m, j, r>::num);
    
    Vh(2*SquareNumNodes2D<i, j, r>::num + 1)
      += OptQuadrangleLobatto<r>::dphi(m, j)*Uh(SquareNumNodes2D<i, m, r>::num);
    
    QuadLobattoLoop3<m-1, i, j, r>::EvalV(Uh, Vh);
  }
  
  template<int m, int i, int j, int r> template<class Vector1, class Vector2>
  inline void QuadLobattoLoop3<m, i, j, r>::EvalU(const Vector1& Vh, Vector2& Uh)
  {
    Uh(SquareNumNodes2D<i, j, r>::num)
      += OptQuadrangleLobatto<r>::dphi(i, m)*Vh(2*SquareNumNodes2D<m, j, r>::num)
      + OptQuadrangleLobatto<r>::dphi(j, m)*Vh(2*SquareNumNodes2D<i, m, r>::num +1);
    
    QuadLobattoLoop3<m-1, i, j, r>::EvalU(Vh, Uh);
  }
  
  template<int i, int j, int r> template<class Vector1, class Vector2>
  inline void QuadLobattoLoop2<i, j, r>::EvalV(const Vector1& Uh, Vector2& Vh)
  {
    Vh(2*SquareNumNodes2D<i, j, r>::num) = 0;
    Vh(2*SquareNumNodes2D<i, j, r>::num + 1) = 0;
    QuadLobattoLoop3<r, i, j, r>::EvalV(Uh, Vh);
    QuadLobattoLoop2<i, j-1, r>::EvalV(Uh, Vh);
  }
  
  template<int i, int j, int r> template<class Vector1, class Vector2>
  inline void QuadLobattoLoop2<i, j, r>::EvalU(const Vector1& Vh, Vector2& Uh)
  {
    Uh(SquareNumNodes2D<i, j, r>::num) = 0;
    QuadLobattoLoop3<r, i, j, r>::EvalU(Vh, Uh);
    QuadLobattoLoop2<i, j-1, r>::EvalU(Vh, Uh);
  }

  template<int i, int r> template<class Vector1, class Vector2>
  inline void QuadLobattoLoop1<i, r>::EvalV(const Vector1& Uh, Vector2& Vh)
  {
    QuadLobattoLoop2<i, r, r>::EvalV(Uh, Vh);
    QuadLobattoLoop1<i-1, r>::EvalV(Uh, Vh);
  }
  
  template<int i, int r> template<class Vector1, class Vector2>
  inline void QuadLobattoLoop1<i, r>::EvalU(const Vector1& Vh, Vector2& Uh)
  {
    QuadLobattoLoop2<i, r, r>::EvalU(Vh, Uh);
    QuadLobattoLoop1<i-1, r>::EvalU(Vh, Uh);
  }

}

#define MONTJOIE_FILE_OPT_SQUARE_LOBATTO_CXX
#endif
