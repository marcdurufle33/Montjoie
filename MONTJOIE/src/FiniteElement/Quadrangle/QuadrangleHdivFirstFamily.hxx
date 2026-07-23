#ifndef MONTJOIE_FILE_QUADRANGLE_HDIV_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! facet finite element of Nedelec's first family on quadrilateral elements
  /*!
    Gauss and Gauss-Lobatto points are used for dof definitions
    Basis functions are written as
    phi_i^GL(x) phi_j^G(y)  e_x   i <= r, j <= r-1
    phi_j^G(x) phi_i^GL(y)  e_y   i <= r, j <= r-1
    where phi_i^GL are Lagrange functions associated with Gauss-Lobatto points
    phi_i^G are Lagrange functions associated with Gauss points
   */
  class QuadrangleHdivFirstFamily : public QuadrangleReference<3>
  {    
  protected :
    //! tensorial to scalar numbering for dofs orientated along X
    Matrix<int> NumDofs_X;
    //! tensorial to scalar numbering for dofs orientated along Y
    Matrix<int> NumDofs_Y;
    //! scalar to tensorial matching arry
    Matrix<int> CoordinateDofs;
    //! On each dof, tangente
    VectR2 tangente_dof;
    Matrix<Real_wp> GL_Gexact;
    Globatto<Real_wp> lob_GL;
    
  public :
    QuadrangleHdivFirstFamily();
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);

  protected:
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();
    
    void ConstructFunctions();

    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
					Vector2& res, int num_loc) const;
    
  public:
    void ComputeValuesPhiRef(const R2& x, VectR2& phi) const;
    void ComputeDivPhiRef(const R2& x, VectReal_wp& phi) const;

    virtual void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
						VectReal_wp& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
						VectComplex_wp& res, int num_loc) const;
        
  };
  
  ostream& operator <<(ostream& out, const QuadrangleHdivFirstFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_QUADRANGLE_HDIV_FIRST_FAMILY_HXX
#endif
