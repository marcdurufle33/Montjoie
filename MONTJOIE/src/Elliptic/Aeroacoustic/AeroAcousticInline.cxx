#ifndef MONTJOIE_FILE_AERO_ACOUSTIC_INLINE_CXX

namespace Montjoie
{

  //! Default constructor for aero-acoustic equation
  template<class Dimension> template<class TypeEquation>
  inline AeroAcoustic_Base<Dimension>::AeroAcoustic_Base(EllipticProblem<TypeEquation>& var)
    : VarGalbrunIndex_Base<Dimension>(var), var_problem(var)
  {
    if (TypeEquation::nb_unknowns == 2+Dimension::dim_N)
      type_model = LEE_MODEL;
    else
      type_model = SIMPLIFIED_LEE;
  }


  template<class Dimension>
  inline FemMatrixFreeClass<Real_wp, StationaryLinearizedEulerEquation<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<StationaryLinearizedEulerEquation<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Real_wp, StationaryLinearizedEulerEquation<Dimension> >(var_)
  {
  }


  template<class Dimension>
  inline FemMatrixFreeClass<Real_wp, AeroStationaryEquation<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<AeroStationaryEquation<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Real_wp, AeroStationaryEquation<Dimension> >(var_)
  {
  }

}

#define MONTJOIE_FILE_AERO_ACOUSTIC_INLINE_CXX
#endif

  
