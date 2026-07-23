#ifndef MONTJOIE_FILE_MONTJOIE_MAXWELL_3D_HEADER_HXX

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtzHeader.hxx"

namespace Montjoie
{
  template<class T>
  class HarmonicMaxwell_3D;
  
  class HarmonicMaxwellEquation_3D;
}

#ifdef MONTJOIE_WITH_TRANSMISSION
#include "Elliptic/Maxwell/3D/TransmissionModelMaxwell3D.hxx"
#endif

#include "Elliptic/Maxwell/3D/HarmonicMaxwell3D.hxx"
#include "Elliptic/Maxwell/3D/ComputationModeMaxwell.hxx"
#include "Elliptic/Maxwell/3D/DefineSourceHarmonicMaxwell3D.hxx"
#include "Elliptic/Maxwell/3D/StaticMaxwell3D.hxx"

#ifdef SELDON_WITH_PRECONDITIONING  
#include "Elliptic/Maxwell/3D/PrecondMaxwell.hxx"
#endif

#include "Computation/ElementaryMatrixHcurl2D.hxx"
#include "Computation/ProdMatVectHcurl2D.hxx"
#include "Computation/ProdMatVectHcurl3D.hxx"

#define MONTJOIE_FILE_MONTJOIE_MAXWELL_3D_HEADER_HXX
#endif
