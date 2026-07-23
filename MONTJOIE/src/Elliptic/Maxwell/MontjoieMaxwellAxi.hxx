#ifndef MONTJOIE_FILE_MONTJOIE_MAXWELL_AXI_HXX

#include "MontjoieMaxwellAxiHeader.hxx"
#include "MontjoieMaxwellAxiInline.hxx"

#include "Harmonic/MontjoieHarmonic.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "Elliptic/Maxwell/Axi/DefineSourceMaxwellAxi.cxx"

#ifdef MONTJOIE_WITH_INTEGRAL_EQUATION
#include "Elliptic/Maxwell/Axi/IntegralEquationMaxwellAxi.cxx"
#endif

#include "Elliptic/Maxwell/Axi/MaxwellAxiSymHarmonic.cxx"
#include "Elliptic/Maxwell/Axi/AxiSymHcurlMaxwell.cxx"

#include "Elliptic/Maxwell/Axi/StaticMaxwellAxi.cxx"

#endif

#define MONTJOIE_FILE_MONTJOIE_MAXWELL_AXI_HXX
#endif
