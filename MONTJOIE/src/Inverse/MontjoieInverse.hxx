#ifndef MONTJOIE_FILE_MONTJOIE_INVERSE_HXX

#include "MontjoieInverseHeader.hxx"
#include "MontjoieInverseInline.hxx"

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"
#include "Elliptic/Elastic/MontjoieElastic.hxx"

#include "Regularization/Itreg.hxx"

#ifndef MONTJOIE_WITH_FAST_LIBRARY

#include "Elliptic/Helmholtz/AxiSymHelmholtz.cxx"

#include "Elliptic/MultiDimPhysicalIndex.cxx"
#include "Inverse/VarMigration.cxx"
#include "Inverse/BoundaryInverseProblem.cxx"

#endif

#define MONTJOIE_FILE_MONTJOIE_INVERSE_HXX
#endif

