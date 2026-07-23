#ifndef MONTJOIE_FILE_MONTJOIE_ELASTIC_HXX

#include "MontjoieElasticHeader.hxx"
#include "MontjoieElasticInline.hxx"

#include "Harmonic/MontjoieHarmonic.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "Elliptic/Elastic/VarElastic.cxx"

#ifdef MONTJOIE_WITH_TRANSMISSION
#include "Elliptic/Elastic/TransmissionModelElastic.cxx"
#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_ELASTIC_HXX
#endif

