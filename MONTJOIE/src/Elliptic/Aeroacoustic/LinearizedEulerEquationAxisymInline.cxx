#ifndef MONTJOIE_FILE_LINEARIZED_EULER_EQUATION_AXISYM_INLINE_CXX

namespace Montjoie
{

  template<class Complexe> template<class TypeEquation>
  inline VarAeroacoustic_Axi<Complexe>::VarAeroacoustic_Axi(EllipticProblem<TypeEquation>& var)
    : VarGalbrun_Axi<Complexe>(var), var_problem(var)
  {
    
    if (TypeEquation::nb_unknowns == 5)
      this->type_model = this->LEE_MODEL;
    else
      this->type_model = this->SIMPLIFIED_LEE;
  }


  template<class TypeEquation>
  inline void VarAeroacousticAxi_Eq<TypeEquation>
  ::ConstructAll(const string& name_file, const string& name_element,
		 Vector<string>& lines_data_file)
  {
    VarAeroacoustic_Axi<Complexe>::ConstructAll(name_file, name_element, lines_data_file);
  }
  

  //! default constructor with given problem
  template<class TypeEquation>
  inline VolumetricSource_AxiAero::
  VolumetricSource_AxiAero(const EllipticProblem<TypeEquation>& var,
			   const Vector<VectString>& param)
    : VirtualSourceFEM<Complex_wp, Dimension2>(var),
      var_aero(var), var_boundary(var)
  { 
    coef_vol = var.GetCoefficientVolumeSource();
    int nb = 2;
    var.InitGaussianParameter(fsrc, param(0), nb);
    var.GetPolarization(polar);
  }
  
  
  template<class TypeEquation>
  inline VolumetricSource<HarmonicAeroacousticAxiEquation>
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_AxiAero(var, param)
  {
  }


  template<class TypeEquation>
  inline VolumetricSource<HarmonicLinearizedEulerEquationAxi>
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_AxiAero(var, param)
  {
  }

}

#define MONTJOIE_FILE_LINEARIZED_EULER_EQUATION_AXISYM_INLINE_CXX
#endif
