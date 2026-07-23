#ifndef MONTJOIE_FILE_WEDGE_DG_CLASSICAL_HXX

namespace Montjoie
{
  
  //! class for nodal finite element on wedges and with DG formulation
  class WedgeDgClassical : public WedgeClassical
  {
  public :
    WedgeDgClassical();
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rtri = 0, int r_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
  public :
    void IntegrateMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
                             const VectReal_wp & coef) const;
    
  };
  
} // end namespace

#define MONTJOIE_FILE_WEDGE_DG_CLASSICAL_HXX
#endif
