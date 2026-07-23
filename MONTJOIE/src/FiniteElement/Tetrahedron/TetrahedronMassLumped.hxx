#ifndef MONTJOIE_FILE_TETRAHEDRON_MASSLUMPED_HXX

namespace Montjoie
{
  
  //! nodal finite element on tetrahedra with mass lumping
  class TetrahedronMassLumped : public TetrahedronReference<1>
  {
  protected :
    
    Vector<MultivariatePolynomial<Real_wp> > BasisFunction;
    Vector<MultivariatePolynomial<Real_wp> > BasisDx_Function, BasisDy_Function, BasisDz_Function;
    VectR2 PointsDof2D, PointsDof2D_Inside;
    TriangleClassical function_basis_tri;
    
  public :
    TetrahedronMassLumped();
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    
    void ConstructFunctions();
    void ComputeLagrangianFunctions();
    
    void ComputeValuesPhiRef(const R3& ptloc, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R3& ptloc, VectR3& dphi) const;
    
    friend ostream& operator <<(ostream& out, const TetrahedronMassLumped& e);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_MASSLUMPED_HXX
#endif
