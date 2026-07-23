#ifndef MONTJOIE_FILE_OPT_TRIANGLE_CLASSICAL_HXX

namespace Montjoie
{
  
  template<int r>
  class OptTriangleClassical
  {
  public :
    enum { nb_dof = (r+1)*(r+2)/2, nb_quad = 3*(r==1) + 6*(r==2) + 12*(r==3) + 16*(r==4)};
    
    TinyMatrix<Real_wp, Symmetric, nb_dof, nb_dof> mat_mass;
    TinyMatrix<Real_wp, General, nb_dof, nb_quad> mat_ch;
    TinyMatrix<Real_wp, General, nb_dof, 2*nb_quad> mat_rigid;
    
    void ConstructFiniteElement(const TriangleClassical&);
    
    template<class Vector1, class Vector2>
    void ApplyCh(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyCh_Transpose(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRh(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRh_Transpose(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void MltMassMatrix(const Vector1&, Vector2&) const;
    
  };
  
  template<>
  class OptTriangleClassical<5>
  {
  public :
    void ConstructFiniteElement(const TriangleClassical&) {}
  };
  
  template<>
  class OptTriangleClassical<6>
  {
  public :
    void ConstructFiniteElement(const TriangleClassical&) {}
  };
  
  template<>
  class OptTriangleClassical<7>
  {
  public :
    void ConstructFiniteElement(const TriangleClassical&) {}
  };
  
  template<>
  class OptTriangleClassical<8>
  {
  public :
    void ConstructFiniteElement(const TriangleClassical&) {}
  };
  
  template<>
  class ImplementationOptimOrder<TriangleClassical, 5>
  {
  public : 
    static const bool implemented = false;
  };
  
}

#define MONTJOIE_FILE_OPT_TRIANGLE_CLASSICAL_HXX
#endif
