#ifndef MONTJOIE_FILE_OPT_HEXAHEDRA_LOBATTO_CXX

namespace Montjoie
{
  
  template<int r> template<class Vector1, class Vector2>
  inline void OptHexahedronLobatto<r>::ApplyCh(const Vector1& Uh, Vector2& Vh)
  {
    Vh = Uh;
  }
    
  
  template<int r> template<class Vector1, class Vector2>
  inline void OptHexahedronLobatto<r>::ApplyChTranspose(const Vector1& Uh, Vector2& Vh)
  {
    Vh = Uh;
  }
  
  
  template<int r>
  void OptHexahedronLobatto<r>::ConstructFiniteElement(const HexahedronLobatto& hex)
  {
    Globatto<Real_wp> lob;
    lob.ConstructQuadrature(r, lob.QUADRATURE_LOBATTO);
    lob.ComputeGradPhi();
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        OptHexahedronLobatto<r>::dphi(i,j) = lob.GradPhi(i, j);

    const Array3D<int>& num = hex.GetNumDofs3D();
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        for (int k = 0; k <= r; k++)
          OptHexahedronLobatto<r>::WeightsND(num(i, j, k)) = lob.Weights(i)*lob.Weights(j)*lob.Weights(k);
    
  }

  template<int r> template<class T>
  void OptHexahedronLobatto<r>::ApplyRh(const Vector<T>& Vh, Vector<T>& Uh)
  {
    HexLobattoLoop1<r, r>::EvalU(Vh.GetData(), Uh.GetData());
  }

  template<int r> template<class T>
  void OptHexahedronLobatto<r>::ApplyRh(const T* Vh, T* Uh)
  {
    HexLobattoLoop1<r, r>::EvalU(Vh, Uh);
  }


  template<int r> template<class Vector1, class Vector2>
  void OptHexahedronLobatto<r>
  ::ApplyRhSplit(const Vector1& Vh, Vector2& Ux, Vector2& Uy, Vector2& Uz)
  {
    HexLobattoLoop1<r, r>::EvalUsplit(Vh, Ux, Uy, Uz);
  }
    
  
  template<int r> template<class T>
  void OptHexahedronLobatto<r>::ApplyRhTranspose(const Vector<T>& Uh, Vector<T>& Vh)
  {
    HexLobattoLoop1<r, r>::EvalV(Uh.GetData(), Vh.GetData());
  }

  
    template<int r> template<class T>
  void OptHexahedronLobatto<r>::ApplyRhTranspose(const T* Uh, T* Vh)
  {
    HexLobattoLoop1<r, r>::EvalV(Uh, Vh);
  }

  
  template<int r>  template<class Vector1, class Vector2>
  void OptHexahedronLobatto<r>::MltMassMatrix(const Vector1& Uh, Vector2& Vh)
  {
    MltDot(Uh, WeightsND, Vh);
  }
  
  
  template<int m, int i, int j, int k, int r> template<class T>
  inline void HexLobattoLoop4<m, i, j, k, r>::EvalV(const T* Uh, T* Vh)
  {
    Vh[3*CubeNumNodes3D<i, j, k, r>::num]
      += OptHexahedronLobatto<r>::dphi(m,i)*Uh[CubeNumNodes3D<m, j, k, r>::num];
    
    Vh[3*CubeNumNodes3D<i, j, k, r>::num + 1]
      += OptHexahedronLobatto<r>::dphi(m, j)*Uh[CubeNumNodes3D<i, m, k, r>::num];
    
    Vh[3*CubeNumNodes3D<i, j, k, r>::num + 2]
      += OptHexahedronLobatto<r>::dphi(m, k)*Uh[CubeNumNodes3D<i, j, m, r>::num];
    
    HexLobattoLoop4<m-1, i, j, k, r>::EvalV(Uh, Vh);
  }
  
  template<int m, int i, int j, int k, int r> template<class T>
  inline void HexLobattoLoop4<m, i, j, k, r>::EvalU(const T* Vh, T* Uh)
  {
    Uh[CubeNumNodes3D<i, j, k, r>::num]
      += OptHexahedronLobatto<r>::dphi(i,m)*Vh[3*CubeNumNodes3D<m, j, k, r>::num]
      + OptHexahedronLobatto<r>::dphi(j,m)*Vh[3*CubeNumNodes3D<i, m, k, r>::num +1]
      + OptHexahedronLobatto<r>::dphi(k,m)*Vh[3*CubeNumNodes3D<i, j, m, r>::num +2];
    
    HexLobattoLoop4<m-1, i, j, k, r>::EvalU(Vh, Uh);
  }

  template<int m, int i, int j, int k, int r> template<class Vector1, class Vector2>
  inline void HexLobattoLoop4<m, i, j, k, r>
  ::EvalUsplit(const Vector1& Vh, Vector2& Ux, Vector2& Uy, Vector2& Uz)
  {
    Ux(CubeNumNodes3D<i, j, k, r>::num)
      += OptHexahedronLobatto<r>::dphi(i,m)*Vh(3*CubeNumNodes3D<m, j, k, r>::num);
    
    Uy(CubeNumNodes3D<i, j, k, r>::num)
      += OptHexahedronLobatto<r>::dphi(j,m)*Vh(3*CubeNumNodes3D<i, m, k, r>::num +1);
    
    Uz(CubeNumNodes3D<i, j, k, r>::num)
      += OptHexahedronLobatto<r>::dphi(k,m)*Vh(3*CubeNumNodes3D<i, j, m, r>::num +2);
    
    HexLobattoLoop4<m-1, i, j, k, r>::EvalUsplit(Vh, Ux, Uy, Uz);
  }
  
  template<int i, int j, int k, int r> template<class T>
  inline void HexLobattoLoop3<i, j, k, r>::EvalV(const T* Uh, T* Vh)
  {
    Vh[3*CubeNumNodes3D<i, j, k, r>::num] = 0;
    Vh[3*CubeNumNodes3D<i, j, k, r>::num + 1] = 0;
    Vh[3*CubeNumNodes3D<i, j, k, r>::num + 2] = 0;
    HexLobattoLoop4<r, i, j, k, r>::EvalV(Uh, Vh);
    HexLobattoLoop3<i, j, k-1, r>::EvalV(Uh, Vh);
  }
  
  template<int i, int j, int k, int r> template<class T>
  inline void HexLobattoLoop3<i, j, k, r>::EvalU(const T* Vh, T* Uh)
  {
    Uh[CubeNumNodes3D<i, j, k, r>::num] = 0;
    HexLobattoLoop4<r, i, j, k, r>::EvalU(Vh, Uh);
    HexLobattoLoop3<i, j, k-1, r>::EvalU(Vh, Uh);
  }

  template<int i, int j, int k, int r> template<class Vector1, class Vector2>
  inline void HexLobattoLoop3<i, j, k, r>
  ::EvalUsplit(const Vector1& Vh, Vector2& Ux, Vector2& Uy, Vector2& Uz)
  {
    Ux(CubeNumNodes3D<i, j, k, r>::num) = 0;
    Uy(CubeNumNodes3D<i, j, k, r>::num) = 0;
    Uz(CubeNumNodes3D<i, j, k, r>::num) = 0;
    
    HexLobattoLoop4<r, i, j, k, r>::EvalUsplit(Vh, Ux, Uy, Uz);
    HexLobattoLoop3<i, j, k-1, r>::EvalUsplit(Vh, Ux, Uy, Uz);
  }

  template<int i, int j, int r> template<class T>
  inline void HexLobattoLoop2<i, j, r>::EvalV(const T* Uh, T* Vh)
  {
    HexLobattoLoop3<i, j, r, r>::EvalV(Uh, Vh);
    HexLobattoLoop2<i, j-1, r>::EvalV(Uh, Vh);
  }
  
  
  template<int i, int j, int r> template<class T>
  inline void HexLobattoLoop2<i, j, r>::EvalU(const T* Vh, T* Uh)
  {
    HexLobattoLoop3<i, j, r, r>::EvalU(Vh, Uh);
    HexLobattoLoop2<i, j-1, r>::EvalU(Vh, Uh);
  }

  template<int i, int j, int r> template<class Vector1, class Vector2>
  inline void HexLobattoLoop2<i, j, r>
  ::EvalUsplit(const Vector1& Vh, Vector2& Ux, Vector2& Uy, Vector2& Uz)
  {
    HexLobattoLoop3<i, j, r, r>::EvalUsplit(Vh, Ux, Uy, Uz);
    HexLobattoLoop2<i, j-1, r>::EvalUsplit(Vh, Ux, Uy, Uz);
  }
  
  template<int i, int r> template<class T>
  inline void HexLobattoLoop1<i, r>::EvalV(const T* Uh, T* Vh)
  {
    HexLobattoLoop2<i, r, r>::EvalV(Uh, Vh);
    HexLobattoLoop1<i-1, r>::EvalV(Uh, Vh);
  }
  
  template<int i, int r> template<class T>
  inline void HexLobattoLoop1<i, r>::EvalU(const T* Vh, T* Uh)
  {
    HexLobattoLoop2<i, r, r>::EvalU(Vh, Uh);
    HexLobattoLoop1<i-1, r>::EvalU(Vh, Uh);
  }

  template<int i, int r> template<class Vector1, class Vector2>
  inline void HexLobattoLoop1<i, r>
  ::EvalUsplit(const Vector1& Vh, Vector2& Ux, Vector2& Uy, Vector2& Uz)
  {
    HexLobattoLoop2<i, r, r>::EvalUsplit(Vh, Ux, Uy, Uz);
    HexLobattoLoop1<i-1, r>::EvalUsplit(Vh, Ux, Uy, Uz);
  }

}

#define MONTJOIE_FILE_OPT_HEXAHEDRA_LOBATTO_CXX
#endif
