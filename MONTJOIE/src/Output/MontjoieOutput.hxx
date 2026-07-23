#ifndef MONTJOIE_FILE_MONTJOIE_OUTPUT_HXX

#include "MontjoieOutputHeader.hxx"
#include "MontjoieOutputInline.hxx"

// common functions
#include "Mesh/MontjoieMesh.hxx"

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "Output/MontjoieOutputCompil.hxx"
#endif

#include "Output/CommonInputOutput.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY

#include "Output/GridInterpolation.cxx"
#include "Output/SplineInterpolation.cxx"
#include "Output/MeshInterpolation.cxx"
#include "Output/ParameterOutput.cxx"
#include "Output/ParameterOutputReprise.cxx"

#endif

#endif

#define MONTJOIE_FILE_MONTJOIE_OUTPUT_HXX
#endif
