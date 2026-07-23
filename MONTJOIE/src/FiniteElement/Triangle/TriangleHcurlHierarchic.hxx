#ifndef MONTJOIE_FILE_TRIANGLE_HCURL_HIERARCHIC_HXX

namespace Montjoie
{
  
  //! Nedelec's second family on triangles, but with hierarchic basis
  class TriangleHcurlHierarchic : public TriangleReference<2>
  {
  protected :
    //! number of dofs based on edges (but inside the element)
    int nb_dof_based_edges;
    //! Legendre polynoms
    Matrix<Real_wp> LegendrePolynom;
    //! scalar <-> tensorial
    Matrix<int> CoordinateDofs;
            
  public :
    TriangleHcurlHierarchic();

    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;

    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
  protected:
    void ConstructFunctions();
        
  public:
    void ComputeValuesPhiRef(const R2& pointloc, VectR2& ) const;
    void ComputeCurlPhiRef(const R2& pointloc, VectReal_wp& ) const;
    
  };
  
  ostream& operator <<(ostream& out,const TriangleHcurlHierarchic& e);
  
} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_HCURL_HIERARCHIC_HXX
#endif

