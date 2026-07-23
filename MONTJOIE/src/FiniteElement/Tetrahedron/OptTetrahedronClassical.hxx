#ifndef MONTJOIE_FILE_OPT_TETRAHEDRON_CLASSICAL_HXX

namespace Montjoie
{
  
  template<int r>
  class OptTetrahedronClassical
  {
  public :
    enum { nb_dof = (r+1)*(r+2)*(r+3)/6, nb_quad = 4*(r==1) + 11*(r==2) + 24*(r==3) + 43*(r==4)};
    static TinyVector<Real_wp, (r+1)*(r+1)*(r+1)> WeightsND;    
    static TinyMatrix<Real_wp, Symmetric, nb_dof, nb_dof> mat_mass;
    static TinyMatrix<Real_wp, General, nb_dof, nb_quad> mat_ch;
    static TinyMatrix<Real_wp, General, nb_dof, 3*nb_quad> mat_rigid;
    
    void ConstructFiniteElement(const TetrahedronClassical&);
    
    template<class Vector1, class Vector2>
    void ApplyCh(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTranspose(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRh(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTranspose(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void MltMassMatrix(const Vector1&, Vector2&) const;
    
  };

  template<int r>
  TinyVector<Real_wp, (r+1)*(r+1)*(r+1)> OptTetrahedronClassical<r>::WeightsND;    
  
  template<int r>
  TinyMatrix<Real_wp, Symmetric, OptTetrahedronClassical<r>::nb_dof,
             OptTetrahedronClassical<r>::nb_dof> OptTetrahedronClassical<r>::mat_mass;
  
  template<int r>
  TinyMatrix<Real_wp, General, OptTetrahedronClassical<r>::nb_dof,
             OptTetrahedronClassical<r>::nb_quad> OptTetrahedronClassical<r>::mat_ch;
  
  template<int r>
  TinyMatrix<Real_wp, General, OptTetrahedronClassical<r>::nb_dof,
             3*OptTetrahedronClassical<r>::nb_quad> OptTetrahedronClassical<r>::mat_rigid;
    
  template<int r>
  class ImplementationOptimOrder<TetrahedronClassical, r>
  {
  public : 
    static const bool implemented = (r<=3);
  };
  
}

#define MONTJOIE_FILE_OPT_TETRAHEDRON_CLASSICAL_HXX
#endif
