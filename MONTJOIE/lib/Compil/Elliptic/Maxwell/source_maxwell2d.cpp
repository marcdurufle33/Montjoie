#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwell2D_Header.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell2D_Inline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/2D/DefineSourceHarmonicMaxwell2D.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template DiffractedWaveSource_Maxwell2D::DiffractedWaveSource_Maxwell2D(const EllipticProblem<HarmonicMaxwellEquation_2D>&, IncidentWaveField<Complex_wp, Dimension2>&);
  SELDON_EXTERN template DiffractedWaveSource<HarmonicMaxwellEquation_2D>::DiffractedWaveSource(const EllipticProblem<HarmonicMaxwellEquation_2D>&, IncidentWaveField<Complex_wp, Dimension2>&);

  SELDON_EXTERN template TotalWaveSource_Maxwell2D::TotalWaveSource_Maxwell2D(const EllipticProblem<HarmonicMaxwellEquation_2D>&, IncidentWaveField<Complex_wp, Dimension2>&);
  SELDON_EXTERN template TotalWaveSource<HarmonicMaxwellEquation_2D>::TotalWaveSource(const EllipticProblem<HarmonicMaxwellEquation_2D>&, IncidentWaveField<Complex_wp, Dimension2>&);
  
  SELDON_EXTERN template IncidentWaveProjector_Maxwell2D::IncidentWaveProjector_Maxwell2D(const EllipticProblem<HarmonicMaxwellEquation_2D>&, IncidentWaveField<Complex_wp, Dimension2>&);
  
}
