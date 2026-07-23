#ifndef MONTJOIE_FILE_OPT_PYRAMID_DG_ORTHO_CXX

namespace Montjoie
{
  template<int r>
  class OptPyramidDgOrtho;
  
  template<int r, int i, int j, int k>
  class PyramidNode
  {
  public :
    enum { dof = k + (i > j)*(r+1-i)*j + (j >= i)*( (r+1-i)*i + (j-i)*(2*r-i-j+3)/2)
           + i*(6*(r+1)*(r+2) - (i-1)*(2*i+2))/12, quad = i*(r+1)*(r+1) + j*(r+1) + k};
    
    
  };

  template<int r, int i, int j>
  class PyramidNode<r, i, j, -1>
  {
  public :
    enum { dof = 0, quad = 0, num = (i==r)*(j==0) + 2*(i==r)*(j==r) + 3*(i==0)*(j==r)
           + (3+i)*(i >= 1)*(i < r)*(j==0) + (3+(r-1)+j)*(j >= 1)*(j < r)*(i==r)
           + (3+2*(r-1)+r-i)*(i >= 1)*(i < r)*(j==r) + (3+3*(r-1)+r-j)*(j >= 1)*(j < r)*(i==0)
           + (4*r + (i-1)*(r-1) + j-1)*(i >= 1)*(i < r)*(j >= 1)*(j < r),
           numz0 = num, numx0 = numz0+(r+1)*(r+1), numy1 = numx0 + (r+1)*(r+1),
           numx1 = numy1 + (r+1)*(r+1), numy0 = numx1 + (r+1)*(r+1)};
    
    
  };
    
  template<int r, int i, int j, int k, int m>
  class PyramidQuadLoop
  {
  public :

    template<class Vector1, class Vector2>
    static inline void ApplyChX(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::quad)
        += OptPyramidDgOrtho<r>::PvX(i, m)*U(PyramidNode<r, m, j, k>::quad);
      PyramidQuadLoop<r, i, j, k, m-1>::ApplyChX(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChY(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::quad)
        += OptPyramidDgOrtho<r>::PvX(j, m)*U(PyramidNode<r, i, m, k>::quad);
      PyramidQuadLoop<r, i, j, k, m-1>::ApplyChY(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChZ(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::dof)
        += OptPyramidDgOrtho<r>::PvZ(i > j ? i : j, k, m)*U(PyramidNode<r, i, j, m>::quad);
      PyramidQuadLoop<r, i, j, k, m-1>::ApplyChZ(U, V);
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChTransX(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::quad)
        += OptPyramidDgOrtho<r>::PvX(m, i)*U(PyramidNode<r, m, j, k>::quad);
      PyramidQuadLoop<r, i, j, k, m-1>::ApplyChTransX(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransY(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::quad)
        += OptPyramidDgOrtho<r>::PvX(m, j)*U(PyramidNode<r, i, m, k>::quad);
      PyramidQuadLoop<r, i, j, k, m-1>::ApplyChTransY(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransZ(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::quad)
        += OptPyramidDgOrtho<r>::PvZ(i > j ? i : j, m, k)*U(PyramidNode<r, i, j, m>::dof);
      PyramidQuadLoop<r, i, j, k, m-1>::ApplyChTransZ(U, V);
    }

  };
  

  template<int r, int i, int j, int k>
  class PyramidQuadLoop<r, i, j, k, -1>
  {
  public :

    template<class Vector1, class Vector2>
    static inline void ApplyChX(const Vector1& U, Vector2& V)
    { }


    template<class Vector1, class Vector2>
    static inline void ApplyChY(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChZ(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChTransX(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransY(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransZ(const Vector1& U, Vector2& V)
    {
    }
  };
  
  template<int r, int i, int j, int k>
  class PyramidTripleLoop
  {
  public :

    template<class Vector1, class Vector2>
    static inline void ApplyChX(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::quad) = 0;
      PyramidQuadLoop<r, i, j, k, r>::ApplyChX(U, V);
      PyramidTripleLoop<r, i, j, k-1>::ApplyChX(U, V);
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChY(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::quad) = 0;
      PyramidQuadLoop<r, i, j, k, r>::ApplyChY(U, V);
      PyramidTripleLoop<r, i, j, k-1>::ApplyChY(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChZ(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::dof) = 0;
      PyramidQuadLoop<r, i, j, k, r>::ApplyChZ(U, V);
      PyramidTripleLoop<r, i, j, k-1>::ApplyChZ(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransX(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::quad) = 0;
      PyramidQuadLoop<r, i, j, k, r>::ApplyChTransX(U, V);
      PyramidTripleLoop<r, i, j, k-1>::ApplyChTransX(U, V);
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChTransY(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::quad) = 0;
      PyramidQuadLoop<r, i, j, k, r>::ApplyChTransY(U, V);
      PyramidTripleLoop<r, i, j, k-1>::ApplyChTransY(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransZ(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, k>::quad) = 0;
      PyramidQuadLoop<r, i, j, k, r-(i > j ? i : j)>::ApplyChTransZ(U, V);
      PyramidTripleLoop<r, i, j, k-1>::ApplyChTransZ(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyShLocTrans(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, -1>::numz0)
        += OptPyramidDgOrtho<r>::ValPhiZ0(k)*U(PyramidNode<r, i, j, k>::quad);
      V(PyramidNode<r, i, j, -1>::numx0)
        += OptPyramidDgOrtho<r>::ValPhiX0(k)*U(PyramidNode<r, i, k, j>::quad);
      V(PyramidNode<r, j, i, -1>::numy1)
        += OptPyramidDgOrtho<r>::ValPhiX1(k)*U(PyramidNode<r, k, j, i>::quad);
      V(PyramidNode<r, i, j, -1>::numx1)
        += OptPyramidDgOrtho<r>::ValPhiX1(k)*U(PyramidNode<r, i, k, j>::quad);
      V(PyramidNode<r, j, i, -1>::numy0)
        += OptPyramidDgOrtho<r>::ValPhiX0(k)*U(PyramidNode<r, k, j, i>::quad);
      PyramidTripleLoop<r, i, j, k-1>::ApplyShLocTrans(U, V);    
    }
    
    template<class T0, class Vector1, class Vector2>
    static inline void ApplyShLoc(const T0& alpha, const Vector1& U, Vector2& V)
    {
      Real_wp val = OptPyramidDgOrtho<r>::ValPhiZ0(k)*U(PyramidNode<r, i, j, -1>::numz0);
      val += OptPyramidDgOrtho<r>::ValPhiX0(j)*U(PyramidNode<r, i, k, -1>::numx0);
      val += OptPyramidDgOrtho<r>::ValPhiX1(i)*U(PyramidNode<r, j, k, -1>::numy1);
      val += OptPyramidDgOrtho<r>::ValPhiX1(j)*U(PyramidNode<r, i, k, -1>::numx1);
      val += OptPyramidDgOrtho<r>::ValPhiX0(i)*U(PyramidNode<r, j, k, -1>::numy0);
      V(PyramidNode<r, i, j, k>::quad) += alpha*val;
      PyramidTripleLoop<r, i, j, k-1>::ApplyShLoc(alpha, U, V);    
    }
  };
  
  
  template<int r, int i, int j>
  class PyramidTripleLoop<r, i, j, -1>
  {
  public :

    template<class Vector1, class Vector2>
    static inline void ApplyChX(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChY(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChZ(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransX(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChTransY(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransZ(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyShLocTrans(const Vector1& U, Vector2& V)
    {
    }
    
    template<class T0, class Vector1, class Vector2>
    static inline void ApplyShLoc(const T0& alpha, const Vector1& U, Vector2& V)
    {
    }
    
  };
  
  
  template<int r, int i, int j>
  class PyramidDoubleLoop
  {
  public :
    template<class Vector1, class Vector2>
    static inline void ApplyChX(const Vector1& U, Vector2& V)
    {
      PyramidTripleLoop<r, i, j, r>::ApplyChX(U, V);
      PyramidDoubleLoop<r, i, j-1>::ApplyChX(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChY(const Vector1& U, Vector2& V)
    {
      PyramidTripleLoop<r, i, j, r>::ApplyChY(U, V);
      PyramidDoubleLoop<r, i, j-1>::ApplyChY(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChZ(const Vector1& U, Vector2& V)
    {
      PyramidTripleLoop<r, i, j, r-(i > j ? i : j)>::ApplyChZ(U, V);
      PyramidDoubleLoop<r, i, j-1>::ApplyChZ(U, V);
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChTransX(const Vector1& U, Vector2& V)
    {
      PyramidTripleLoop<r, i, j, r>::ApplyChTransX(U, V);
      PyramidDoubleLoop<r, i, j-1>::ApplyChTransX(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransY(const Vector1& U, Vector2& V)
    {
      PyramidTripleLoop<r, i, j, r>::ApplyChTransY(U, V);
      PyramidDoubleLoop<r, i, j-1>::ApplyChTransY(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransZ(const Vector1& U, Vector2& V)
    {
      PyramidTripleLoop<r, i, j, r>::ApplyChTransZ(U, V);
      PyramidDoubleLoop<r, i, j-1>::ApplyChTransZ(U, V);
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyShLocTrans(const Vector1& U, Vector2& V)
    {
      V(PyramidNode<r, i, j, -1>::numz0) = 0;
      V(PyramidNode<r, i, j, -1>::numx0) = 0;
      V(PyramidNode<r, j, i, -1>::numy1) = 0;
      V(PyramidNode<r, i, j, -1>::numx1) = 0;
      V(PyramidNode<r, j, i, -1>::numy0) = 0;
      PyramidTripleLoop<r, i, j, r>::ApplyShLocTrans(U, V);
      PyramidDoubleLoop<r, i, j-1>::ApplyShLocTrans(U, V);
    }
    
    template<class T0, class Vector1, class Vector2>
    static inline void ApplyShLoc(const T0& alpha, const Vector1& U, Vector2& V)
    {
      PyramidTripleLoop<r, i, j, r>::ApplyShLoc(alpha, U, V);
      PyramidDoubleLoop<r, i, j-1>::ApplyShLoc(alpha, U, V);
    }
  };
  
  template<int r, int i>
  class PyramidDoubleLoop<r, i, -1>
  {
  public :
    template<class Vector1, class Vector2>
    static inline void ApplyChX(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChY(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChZ(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChTransX(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransY(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransZ(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyShLocTrans(const Vector1& U, Vector2& V)
    {
    }

    template<class T0, class Vector1, class Vector2>
    static inline void ApplyShLoc(const T0& alpha, const Vector1& U, Vector2& V)
    {
    }
  };
  
  template<int r, int i>
  class PyramidLoop
  {
  public :
    template<class Vector1, class Vector2>
    static inline void ApplyChX(const Vector1& U, Vector2& V)
    {
      PyramidDoubleLoop<r, i, r>::ApplyChX(U, V);
      PyramidLoop<r, i-1>::ApplyChX(U, V);
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChY(const Vector1& U, Vector2& V)
    {
      PyramidDoubleLoop<r, i, r>::ApplyChY(U, V);
      PyramidLoop<r, i-1>::ApplyChY(U, V);
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChZ(const Vector1& U, Vector2& V)
    {
      PyramidDoubleLoop<r, i, r>::ApplyChZ(U, V);
      PyramidLoop<r, i-1>::ApplyChZ(U, V);
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransX(const Vector1& U, Vector2& V)
    {
      PyramidDoubleLoop<r, i, r>::ApplyChTransX(U, V);
      PyramidLoop<r, i-1>::ApplyChTransX(U, V);
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChTransY(const Vector1& U, Vector2& V)
    {
      PyramidDoubleLoop<r, i, r>::ApplyChTransY(U, V);
      PyramidLoop<r, i-1>::ApplyChTransY(U, V);
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChTransZ(const Vector1& U, Vector2& V)
    {
      PyramidDoubleLoop<r, i, r>::ApplyChTransZ(U, V);
      PyramidLoop<r, i-1>::ApplyChTransZ(U, V);
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyShLocTrans(const Vector1& U, Vector2& V)
    {
      PyramidDoubleLoop<r, i, r>::ApplyShLocTrans(U, V);
      PyramidLoop<r, i-1>::ApplyShLocTrans(U, V);
    }

    template<class T0, class Vector1, class Vector2>
    static inline void ApplyShLoc(const T0& alpha, const Vector1& U, Vector2& V)
    {
      PyramidDoubleLoop<r, i, r>::ApplyShLoc(alpha, U, V);
      PyramidLoop<r, i-1>::ApplyShLoc(alpha, U, V);
    }
  };
  
  template<int r>
  class PyramidLoop<r, -1>
  {
  public :
    template<class Vector1, class Vector2>
    static inline void ApplyChX(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChY(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChZ(const Vector1& U, Vector2& V)
    {
    }

    template<class Vector1, class Vector2>
    static inline void ApplyChTransX(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChTransY(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyChTransZ(const Vector1& U, Vector2& V)
    {
    }
    
    template<class Vector1, class Vector2>
    static inline void ApplyShLocTrans(const Vector1& U, Vector2& V)
    { }

    template<class T0, class Vector1, class Vector2>
    static inline void ApplyShLoc(const T0& alpha, const Vector1& U, Vector2& V)
    { }
    
  };
  
  template<int r>
  class OptPyramidDgOrtho
  {
  public :
    enum { nb_quad = (r+1)*(r+1)*(r+1), nb_dof = (r+1)*(r+2)*(2*r+3)/6,
           nb_quad2d = (r+1)*(r+1), nb_tri = (r==1)*3 + (r==2)*6 + (r==3)*12  + (r==4)*16
           + (r==5)*25 + (r==6)*33 + (r==7)*42 + (r==8)*52 };
    
    static TinyMatrix<Real_wp, General, r+1, r+1> PvX;
    static TinyArray3D<Real_wp, r+1, r+1, r+1> PvZ;
    static TinyVector<Real_wp, r+1> ValPhiX0, ValPhiX1, ValPhiZ0;
    static TinyVector<TinyMatrix<Real_wp, General, (r+1)*(r+1),
                                 (r==1)*3 + (r==2)*6 + (r==3)*12  + (r==4)*16
                                 + (r==5)*25 + (r==6)*33 + (r==7)*42 + (r==8)*52>, 4> ShSurf;
    
    template<class Vector1, class Vector2>
    static void ApplyCh(const Vector1& U, Vector2& V)
    {
      TinyVector<Real_wp, nb_quad> Ux;
      TinyVector<Real_wp, nb_quad> Uy;
      
      PyramidLoop<r, r>::ApplyChX(U, Ux);
      PyramidLoop<r, r>::ApplyChY(Ux, Uy);
      PyramidLoop<r, r>::ApplyChZ(Uy, V);
    }
    
    template<class Vector1, class Vector2>
    static void ApplyChTranspose(const Vector1& U, Vector2& V)
    {
      TinyVector<Real_wp, nb_quad> Ux;
      TinyVector<Real_wp, nb_quad> Uy;
      
      PyramidLoop<r, r>::ApplyChTransZ(U, Uy);
      PyramidLoop<r, r>::ApplyChTransY(Uy, Ux);
      PyramidLoop<r, r>::ApplyChTransX(Ux, V);
    }
    
    static void CopyPvX(const Array3D<Real_wp>& Pv)
    {
      for (int i = 0; i <= r; i++)
        for (int j = 0; j <= r; j++)
          PvX(i, j) = Pv(0, i, j);
      
    }
    
    static void CopyPvZ(const Array3D<Real_wp>& Pv)
    {
      for (int i = 0; i <= r; i++)
        for (int j = 0; j <= r; j++)
          for (int k = 0; k <= r; k++)
            PvZ(i, j, k) = Pv(i, j, k);      
    }
    
    static void CopyShLoc(const VectReal_wp& phix0, const VectReal_wp& phix1,
                          const VectReal_wp& phiz0)
    {
      for (int i = 0; i <= r; i++)
        {
          ValPhiX0(i) = phix0(i);
          ValPhiX1(i) = phix1(i);
          ValPhiZ0(i) = phiz0(i);
        }
    }
    
    template<class Matrix1>
    static void CopyShSurf(const Matrix1& ShSurf2)
    {
      for (int num = 0; num < 4; num++)
        {
          for (int i = 0; i < ShSurf(num).GetM(); i++)
            for (int j = 0; j < ShSurf(num).GetN(); j++)
              ShSurf(num)(i, j) = ShSurf2(num)(i, j);
        }
    }
    
    template<class Vector1, class Vector2>
    static void ApplyShQuadratureTranspose(const Vector1& U, Vector2& V)
    {
      TinyVector<Real_wp, 5*(r+1)*(r+1)> Usurf;
      PyramidLoop<r, r>::ApplyShLocTrans(U, Usurf);
      
      //DISP(Usurf);
      for (int i = 0; i < nb_quad2d; i++)
        V(i) = Usurf(i);
      
      TinyVector<Real_wp, nb_tri> Utri;
      TinyVector<Real_wp, nb_quad2d> Uquad;
      for (int num_loc = 1; num_loc <= 4; num_loc++)
        {
          int offset = nb_quad2d*num_loc;
          for (int i = 0; i < nb_quad2d; i++)
            Uquad(i) = Usurf(offset+i);
          
          MltTrans(ShSurf(num_loc-1), Uquad, Utri);
          //DISP(Uquad); DISP(Utri);
          offset = nb_quad2d + nb_tri*(num_loc-1);
          for (int i = 0; i < nb_tri; i++)
            V(offset+i) = Utri(i);
        }
    }
    
    template<class T0, class Vector1, class Vector2>
    static void ApplyShQuadrature(const T0& alpha, const Vector1& U, Vector2& V)
    {
      TinyVector<Real_wp, nb_tri> Utri;
      TinyVector<Real_wp, nb_quad2d> Uquad;
      TinyVector<Real_wp, 5*nb_quad2d> Usurf;
      for (int i = 0; i < nb_quad2d; i++)
        Usurf(i) = U(i);
      
      for (int num_loc = 1; num_loc <= 4; num_loc++)
        {
          int offset = nb_quad2d + nb_tri*(num_loc-1);
          for (int i = 0; i < nb_tri; i++)
            Utri(i) = U(offset+i);
          
          Mlt(ShSurf(num_loc-1), Utri, Uquad);
          offset = nb_quad2d*num_loc;
          for (int i = 0; i < nb_quad2d; i++)
            Usurf(offset+i) = Uquad(i);
        }
      
      PyramidLoop<r, r>::ApplyShLoc(alpha, Usurf, V);
    }
    
  };
  
  template<int r> TinyMatrix<Real_wp, General, r+1, r+1> OptPyramidDgOrtho<r>::PvX;
  template<int r> TinyArray3D<Real_wp, r+1, r+1, r+1> OptPyramidDgOrtho<r>::PvZ;
  template<int r> TinyVector<Real_wp, r+1> OptPyramidDgOrtho<r>::ValPhiX0;
  template<int r> TinyVector<Real_wp, r+1> OptPyramidDgOrtho<r>::ValPhiX1;
  template<int r> TinyVector<Real_wp, r+1> OptPyramidDgOrtho<r>::ValPhiZ0;
  
  template<int r>
  TinyVector<TinyMatrix<Real_wp, General, (r+1)*(r+1),
                        (r==1)*3 + (r==2)*6 + (r==3)*12  + (r==4)*16
                        + (r==5)*25 + (r==6)*33 + (r==7)*42 + (r==8)*52>, 4>
  OptPyramidDgOrtho<r>::ShSurf;
  
  void PyramidDgOrtho::CopyPvX(const Array3D<Real_wp>& Pv)
  {
    switch (order)
      {
      case 1 :
        OptPyramidDgOrtho<1>::CopyPvX(Pv);
        break;
      case 2 :
        OptPyramidDgOrtho<2>::CopyPvX(Pv);
        break;
      case 3 :
        OptPyramidDgOrtho<3>::CopyPvX(Pv);
        break;
        
      case 4 :
        OptPyramidDgOrtho<4>::CopyPvX(Pv);
        break;
      case 5 :
        OptPyramidDgOrtho<5>::CopyPvX(Pv);
        break;
      case 6 :
        OptPyramidDgOrtho<6>::CopyPvX(Pv);
        break;
        /*
      case 7 :
        OptPyramidDgOrtho<7>::CopyPvX(Pv);
        break;
      case 8 :
        OptPyramidDgOrtho<8>::CopyPvX(Pv);
        break;
        */
      }
  }
  
  
  void PyramidDgOrtho::CopyPvZ(const Array3D<Real_wp>& Pv)
  {
    switch (order)
      {
      case 1 :
        OptPyramidDgOrtho<1>::CopyPvZ(Pv);
        break;
      case 2 :
        OptPyramidDgOrtho<2>::CopyPvZ(Pv);
        break;
      case 3 :
        OptPyramidDgOrtho<3>::CopyPvZ(Pv);
        break;      
      case 4 :
        OptPyramidDgOrtho<4>::CopyPvZ(Pv);
        break;
      case 5 :
        OptPyramidDgOrtho<5>::CopyPvZ(Pv);
        break;
      case 6 :
        OptPyramidDgOrtho<6>::CopyPvZ(Pv);
        break;
        /*
      case 7 :
        OptPyramidDgOrtho<7>::CopyPvZ(Pv);
        break;
      case 8 :
        OptPyramidDgOrtho<8>::CopyPvZ(Pv);
        break;
        */
      }
  }
  
  template<class Matrix1>
  void PyramidDgOrtho::CopyShLoc(const Matrix1& A, const Matrix1& B, const Matrix1& C)
  {
    switch (order)
      {
      case 1 :
        OptPyramidDgOrtho<1>::CopyShLoc(A, B, C);
        break;
      case 2 :
        OptPyramidDgOrtho<2>::CopyShLoc(A, B, C);
        break;
      case 3 :
        OptPyramidDgOrtho<3>::CopyShLoc(A, B, C);
        break;        
      case 4 :
        OptPyramidDgOrtho<4>::CopyShLoc(A, B, C);
        break;
      case 5 :
        OptPyramidDgOrtho<5>::CopyShLoc(A, B, C);
        break;
      case 6 :
        OptPyramidDgOrtho<6>::CopyShLoc(A, B, C);
        break;
        /*
      case 7 :
        OptPyramidDgOrtho<7>::CopyShLoc(A, B, C);
        break;
      case 8 :
        OptPyramidDgOrtho<8>::CopyShLoc(A, B, C);
        break;
        */
      }
  }
  
  
  template<class Matrix1>
  void PyramidDgOrtho::CopyShSurf(const Matrix1& A)
  {
    switch (order)
      {
      case 1 :
        OptPyramidDgOrtho<1>::CopyShSurf(A);
        break;
      case 2 :
        OptPyramidDgOrtho<2>::CopyShSurf(A);
        break;
      case 3 :
        OptPyramidDgOrtho<3>::CopyShSurf(A);
        break;
      case 4 :
        OptPyramidDgOrtho<4>::CopyShSurf(A);
        break;
      case 5 :
        OptPyramidDgOrtho<5>::CopyShSurf(A);
        break;
      case 6 :
        OptPyramidDgOrtho<6>::CopyShSurf(A);
        break;
        /*
      case 7 :
        OptPyramidDgOrtho<7>::CopyShSurf(A);
        break;
      case 8 :
        OptPyramidDgOrtho<8>::CopyShSurf(A);
        break;
        */
      }
  }
  
  
  template<class Vector1, class Vector2>
  void PyramidDgOrtho::ApplyShQuadratureTranspose(const Vector1& U, Vector2& V) const
  {
    switch (order)
      {
      case 1 :
        OptPyramidDgOrtho<1>::ApplyShQuadratureTranspose(U, V);
        break;
      case 2 :
        OptPyramidDgOrtho<2>::ApplyShQuadratureTranspose(U, V);
        break;
      case 3 :
        OptPyramidDgOrtho<3>::ApplyShQuadratureTranspose(U, V);
        break;
      case 4 :
        OptPyramidDgOrtho<4>::ApplyShQuadratureTranspose(U, V);
        break;
      case 5 :
        OptPyramidDgOrtho<5>::ApplyShQuadratureTranspose(U, V);
        break;
      case 6 :
        OptPyramidDgOrtho<6>::ApplyShQuadratureTranspose(U, V);
        break;
      case 7 :
        break;
      case 8 :
        break;
      }
  }

  template<class T0, class Vector1, class Vector2>
  void PyramidDgOrtho::ApplyShQuadrature(const T0& alpha, const Vector1& U, Vector2& V) const
  {
    switch (order)
      {
      case 1 :
        OptPyramidDgOrtho<1>::ApplyShQuadrature(alpha, U, V);
        break;
      case 2 :
        OptPyramidDgOrtho<2>::ApplyShQuadrature(alpha, U, V);
        break;
      case 3 :
        OptPyramidDgOrtho<3>::ApplyShQuadrature(alpha, U, V);
        break;
      case 4 :
        OptPyramidDgOrtho<4>::ApplyShQuadrature(alpha, U, V);
        break;
      case 5 :
        OptPyramidDgOrtho<5>::ApplyShQuadrature(alpha, U, V);
        break;
      case 6 :
        OptPyramidDgOrtho<6>::ApplyShQuadrature(alpha, U, V);
        break;
      case 7 :
        break;
      case 8 :
        break;
      }
  }
  
  template<class Vector1, class Vector2>
  void PyramidDgOrtho::ApplyCh(const Vector1& U, Vector2& V) const
  {
    //Vector1 Ux(ChX.GetM()), Uy(ChY.GetM());
    //Ux.Fill(0); Uy.Fill(0);
    //Mlt(ChX, U, Ux);
    //Mlt(ChY, Ux, Uy);
    //Mlt(ChZ, Uy, V);
    
    //return;
    
    switch (order)
      {
      case 1 :
        OptPyramidDgOrtho<1>::ApplyCh(U, V);
        break;
      case 2 :
        OptPyramidDgOrtho<2>::ApplyCh(U, V);
        break;
      case 3 :
        OptPyramidDgOrtho<3>::ApplyCh(U, V);
        break;
      case 4 :
        OptPyramidDgOrtho<4>::ApplyCh(U, V);
        break;
      case 5 :
        OptPyramidDgOrtho<5>::ApplyCh(U, V);
        break;
      case 6 :
        OptPyramidDgOrtho<6>::ApplyCh(U, V);
        break;
      case 7 :
        break;
      case 8 :
        break;
      }
  }
   
  
  template<class Vector1, class Vector2>
  void PyramidDgOrtho::ApplyChTranspose(const Vector1& U, Vector2& V) const
  {
    //if (!UseQuadraturePointsForSh())
    //return PyramidReference::ApplyChTranspose(U, V);
    
    //Vector1 Uy(ChTransposeY.GetM()), Uz(ChTransposeZ.GetM());
    //Uy.Fill(0); Uz.Fill(0);
    //Mlt(ChTransposeZ, U, Uz);
    //Mlt(ChTransposeY, Uz, Uy);
    //Mlt(ChTransposeX, Uy, V);
    
    //return;
    
    switch (order)
      {
      case 1 :
        OptPyramidDgOrtho<1>::ApplyChTranspose(U, V);
        break;
      case 2 :
        OptPyramidDgOrtho<2>::ApplyChTranspose(U, V);
        break;
      case 3 :
        OptPyramidDgOrtho<3>::ApplyChTranspose(U, V);
        break;
      case 4 :
        OptPyramidDgOrtho<4>::ApplyChTranspose(U, V);
        break;
      case 5 :
        OptPyramidDgOrtho<5>::ApplyChTranspose(U, V);
        break;
      case 6 :
        OptPyramidDgOrtho<6>::ApplyChTranspose(U, V);
        break;
      case 7 :
        break;
      case 8 :
        break;
      }
  }
    
}

#define MONTJOIE_FILE_OPT_PYRAMID_DG_ORTHO_CXX
#endif
