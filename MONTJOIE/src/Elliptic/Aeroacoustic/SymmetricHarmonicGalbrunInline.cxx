#ifndef MONTJOIE_FILE_SYMMETRIC_HARMONIC_GALBRUN_INLINE_CXX

namespace Montjoie
{

  template<class T, class Dimension> template<class TypeEquation>
  inline ImpedanceABC<T, HarmonicGalbrunEquationSipg<Dimension> >
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var) 
    : ImpedanceFunction_Base<T, Dimension>(var), var_problem(var), var_galbrun(var)
  {
  }
  

  template<class T> template<class TypeEquation>
  inline ImpedanceABC<T, HarmonicGalbrunEquationH1>
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceABC<T, HarmonicGalbrunEquationSipg<Dimension2> >(var)
  {
  }


  template<class T> template<class TypeEquation>
  inline VolumetricSource_GalbrunH1<T>
  ::VolumetricSource_GalbrunH1(const EllipticProblem<TypeEquation>& var,
			       const Vector<VectString>& param)
    : VirtualSourceFEM<T, Dimension2>(var), var_problem(var), var_galbrun(var), var_galbrun_base(var)
  {
    coef_vol = var.GetCoefficientVolumeSource();
    int nb = 2;
    var.InitGaussianParameter(fsrc, param(0), nb);
  }

  
  template<class TypeEquation>
  inline VolumetricSource<HarmonicGalbrunEquationH1>
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_GalbrunH1<Complex_wp>(var, param)
  {
  }

}

#define MONTJOIE_FILE_SYMMETRIC_HARMONIC_GALBRUN_INLINE_CXX
#endif
