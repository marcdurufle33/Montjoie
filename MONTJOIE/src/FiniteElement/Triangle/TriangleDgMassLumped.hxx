#ifndef MONTJOIE_FILE_TRIANGLE_DG_MASSLUMPED_HXX

namespace Montjoie
{
    
  //! class for nodal triangular finite element with mass-lumping and with DG formulation
  /*!
    The meshes can be hybrid (triangles + quadrilaterals)
   */
  class TriangleMassLumped : public TriangleClassical
  {
  protected :
    //! basis functions
    Vector<MultivariatePolynomial<Real_wp> >
    LagrangeFunction, DxLagrangeFunction, DyLagrangeFunction;

  public :
    TriangleMassLumped();
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
    void ComputeLagrangianBubbleFunctions(int r);
    
    void ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R2& pointloc, VectR2& grad_phi) const;
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_DG_MASSLUMPED_HXX
#endif

