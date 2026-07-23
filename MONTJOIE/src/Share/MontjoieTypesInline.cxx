#ifndef MONTJOIE_FILE_MONTJOIE_TYPES_INLINE_CXX

namespace Montjoie
{

  template<class TypeEquation>
  inline DiffractedWaveSource<TypeEquation>
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& var,
			 IncidentWaveField<typename TypeEquation::Complexe,
			 typename TypeEquation::Dimension>& fsrc_)
    : VirtualSourceFEM<typename TypeEquation::Complexe,
		       typename TypeEquation::Dimension>(var)
  {
  }


  template<class TypeEquation>
  inline TotalWaveSource<TypeEquation>
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& var,
		    IncidentWaveField<typename TypeEquation::Complexe,
		    typename TypeEquation::Dimension>& fsrc_)
    : VirtualSourceFEM<typename TypeEquation::Complexe,
		       typename TypeEquation::Dimension>(var)	
  {
  }
  
}

#define MONTJOIE_FILE_MONTJOIE_TYPES_INLINE_CXX
#endif
