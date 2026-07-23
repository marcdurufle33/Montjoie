#ifndef MONTJOIE_FILE_PROD_MAT_VECT_HELMHOLTZ_HXX

namespace Montjoie
{
  
  template<class Dimension>
  class FemMatrixFreeClass<Complex_wp, HelmholtzEquation<Dimension> >
    : public FemMatrixFreeClass_Eq<Complex_wp, HelmholtzEquation<Dimension> >
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<HelmholtzEquation<Dimension> >& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };


  template<class Dimension>
  class FemMatrixFreeClass<Complex_wp, HelmholtzEquationDG<Dimension> >
    : public FemMatrixFreeClass_Eq<Complex_wp, HelmholtzEquationDG<Dimension> >
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<HelmholtzEquationDG<Dimension> >& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };

}

#define MONTJOIE_FILE_PROD_MAT_VECT_HELMHOLTZ_HXX
#endif
