#include "Montjoie.hxx"

namespace Montjoie
{
  class MultiString_Base;
  class ParameterMultistring;
}

#include "Corde/DonneesString.hxx"

#include "Corde/MegaNewtonSolver.hxx"
#include "Corde/SpaceStringScheme.hxx"
#include "Corde/PianoString.hxx"
#include "Corde/EquationString.hxx"
#include "Corde/Hammer.hxx"
#include "Corde/Shank.hxx"
#include "Corde/Bridge.hxx"
#include "Corde/MultiString.hxx"

#include "Corde/MultiStringInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Corde/DonneesString.cxx"
#else
#include "vector/Vector.cxx"
#endif

#include "Corde/MegaNewtonSolver.cxx"
#include "Corde/SpaceStringScheme.cxx"
#include "Corde/PianoString.cxx"
#include "Corde/Hammer.cxx"
#include "Corde/Shank.cxx"
#include "Corde/Bridge.cxx"
#include "Corde/MultiString.cxx"
#include "Corde/ComputeSchemeAndDiffLM.cxx"
