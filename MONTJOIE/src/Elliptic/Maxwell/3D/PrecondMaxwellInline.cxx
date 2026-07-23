#ifndef MONTJOIE_FILE_PRECOND_MAXWELL_INLINE_CXX

namespace Montjoie
{

  template<class T> template<class TypeEquation>
  inline LowOrderPreconditioning_Maxwell<T>
  ::LowOrderPreconditioning_Maxwell(EllipticProblem<TypeEquation>& var)
    : LowOrderPreconditioning_Dim<T, Dimension3>(var)
  {    
  }
  

  template<class TypeEquation>
  inline LowOrderPreconditioning<Complex_wp, HarmonicMaxwellEquation_3D>
  ::LowOrderPreconditioning(EllipticProblem<TypeEquation>& var)
    : LowOrderPreconditioning_Maxwell<Complex_wp>(var)
  {    
  }


  template<class T> template<class TypeEquation>
  inline MultigridPreconditioning_Maxwell<T>
  ::MultigridPreconditioning_Maxwell(EllipticProblem<TypeEquation>& var)
    : MultigridPreconditioning_Dim<T, Dimension3>(var), var_fine(var)
  {
  }
  
  
  template<class TypeEquation>
  inline MultigridPreconditioning<Complex_wp, HarmonicMaxwellEquation_3D>
  ::MultigridPreconditioning(EllipticProblem<TypeEquation>& var)
    : MultigridPreconditioning_Maxwell<Complex_wp>(var)
  {
  }


  template<class T> template<class TypeEquation>    
  inline All_Preconditioner_Maxwell<T>
  ::All_Preconditioner_Maxwell(EllipticProblem<TypeEquation>& var)
    : All_Preconditioner_Base<T>(var)
  {
  }


  template<class T> template<class TypeEquation>    
  inline All_Preconditioner<T, HarmonicMaxwellEquation_3D>
  ::All_Preconditioner(EllipticProblem<TypeEquation>& var)
    : All_Preconditioner_Maxwell<T>(var)
  {
  }
  
}
  
#define MONTJOIE_FILE_PRECOND_MAXWELL_INLINE_CXX
#endif
