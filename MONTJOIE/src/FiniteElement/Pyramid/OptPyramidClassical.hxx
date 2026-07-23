#ifndef MONTJOIE_FILE_OPT_PYRAMID_CLASSICAL_HXX

namespace Montjoie
{
  
  template<int r>
  class OptPyramidClassical
  {
  public :
    enum { nb_dof = (r+1)*(r+2)*(2*r+3)/6, nb_quad = (r+1)*(r+1)*(r+1)};
    static TinyVector<Real_wp, (r+1)*(r+1)*(r+1)> WeightsND;    
    static TinyMatrix<Real_wp, Symmetric, nb_dof, nb_dof> mat_mass;
    static TinyMatrix<Real_wp, General, nb_dof, nb_quad> mat_ch;
    static TinyMatrix<Real_wp, General, nb_dof, 3*nb_quad> mat_rigid;
    
    static void ConstructFiniteElement(const PyramidClassical&);
    
    template<class Vector1, class Vector2>
    static void ApplyCh(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void ApplyChTranspose(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void ApplyRh(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void ApplyRhTranspose(const Vector1& Uh, Vector2& Vh);
    
    template<class Vector1, class Vector2>
    static void MltMassMatrix(const Vector1&, Vector2&);
    
  };

  template<int r> 
  TinyVector<Real_wp, (r+1)*(r+1)*(r+1)> OptPyramidClassical<r>::WeightsND;    
  
  template<int r>
  TinyMatrix<Real_wp, Symmetric, OptPyramidClassical<r>::nb_dof,
             OptPyramidClassical<r>::nb_dof> OptPyramidClassical<r>::mat_mass;
  
  template<int r>
  TinyMatrix<Real_wp, General, OptPyramidClassical<r>::nb_dof,
             OptPyramidClassical<r>::nb_quad> OptPyramidClassical<r>::mat_ch;
  
  template<int r>
  TinyMatrix<Real_wp, General, OptPyramidClassical<r>::nb_dof,
             3*OptPyramidClassical<r>::nb_quad> OptPyramidClassical<r>::mat_rigid;
    
}

#define MONTJOIE_FILE_OPT_PYRAMID_CLASSICAL_HXX
#endif
