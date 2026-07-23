#ifndef MONTJOIE_FILE_AXISYM_GALBRUN_INLINE_CXX

namespace Montjoie
{

  //! default constructor
  template<class TypeEquation>
  inline VarGalbrunIndex_Axi::VarGalbrunIndex_Axi(EllipticProblem<TypeEquation>& var)
    : VarAxisymProblem(var), var_problem(var), var_output(var)
  {
    apply_convective_derivate_source = false;
    
    this->type_model = GALBRUN_DIV;
    drop_unstable_terms = DROP_NONE;
    coef_convective_term = 1.0;
    compute_gravity = false; //AJOUT NATHAN
  }
  

  template<class Complexe> template<class TypeEquation>
  inline VarGalbrun_Axi<Complexe>::VarGalbrun_Axi(EllipticProblem<TypeEquation>& var)
    : VarGalbrunIndex_Axi(var), var_problem(var), var_boundary(var)
  {
  }

  template<class TypeEquation>
  inline void VarGalbrunAxi_Eq<TypeEquation>
  ::ConstructAll(const string& name_file, const string& name_element,
		 Vector<string>& lines_data_file)
  {
    VarGalbrun_Axi<Complexe>::ConstructAll(name_file, name_element, lines_data_file);
  }


  //! default constructor with given problem
  template<class TypeEquation>
  inline VolumetricSource_AxiGalbrun::
  VolumetricSource_AxiGalbrun(const EllipticProblem<TypeEquation>& var,
			      const Vector<VectString>& param)
    : VirtualSourceFEM<Complex_wp, Dimension2>(var), var_problem(var),
      var_galbrun(var), var_boundary(var)
  { 
    coef_vol = var.GetCoefficientVolumeSource();
    int nb = 2;
    var.InitGaussianParameter(fsrc, param(0), nb);
    var.GetPolarization(polar);
  }


  template<class TypeEquation>
  inline VolumetricSource<HarmonicGalbrunEquationAxi>
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_AxiGalbrun(var, param)
  {
  }

}

#define MONTJOIE_FILE_AXISYM_GALBRUN_INLINE_CXX
#endif
