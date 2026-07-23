#ifndef MONTJOIE_FILE_OPT_WEDGE_CLASSICAL_HXX

namespace Montjoie
{
  
  template<int r>
  class OptWedgeClassical
  {
  public :
    enum { nb_dof = (r+1)*(r+1)*(r+2)/2, nb_quad = (r+1)*(r+1)*(r+1) };
    static TinyMatrix<Real_wp, Symmetric, nb_dof, nb_dof> mat_mass;
    static TinyMatrix<Real_wp, General, nb_dof, nb_quad> mat_ch;
    static TinyMatrix<Real_wp, General, nb_dof, 3*nb_quad> mat_rigid;
    static TinyVector<Real_wp, (r+1)*(r+1)*(r+1)> WeightsND;
    
    static void ConstructFiniteElement(const WedgeClassical&);
    
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
  TinyMatrix<Real_wp, Symmetric, OptWedgeClassical<r>::nb_dof,
             OptWedgeClassical<r>::nb_dof> OptWedgeClassical<r>::mat_mass;
  
  template<int r>
  TinyMatrix<Real_wp, General, OptWedgeClassical<r>::nb_dof,
             OptWedgeClassical<r>::nb_quad> OptWedgeClassical<r>::mat_ch;
  
  template<int r>
  TinyMatrix<Real_wp, General, OptWedgeClassical<r>::nb_dof,
             3*OptWedgeClassical<r>::nb_quad> OptWedgeClassical<r>::mat_rigid;
  
  template<int r>
  TinyVector<Real_wp, (r+1)*(r+1)*(r+1)> OptWedgeClassical<r>::WeightsND;
    
}

#define MONTJOIE_FILE_OPT_WEDGE_CLASSICAL_HXX
#endif
