#ifndef MONTJOIE_FILE_TRIANGLE_HIERARCHIC_INLINE_CXX

namespace Montjoie
{
  
  inline void TriangleHierarchic::SetBasisType(int type)
  {
    type_basis = type;
  }
  
  inline const Matrix<int>& TriangleHierarchic::GetNumDofs2D() const
  {
    return NumDofsTri;
  }

  inline void TriangleHierarchic
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
    
  inline void TriangleHierarchic
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
}

#define MONTJOIE_FILE_TRIANGLE_HIERARCHIC_INLINE_CXX
#endif
