#ifndef MONTJOIE_FILE_MONTJOIE_VLASOV_MAXWELL_HXX

#include "Elliptic/Helmholtz/MontjoieLaplace.hxx"
#include "Hyperbolic/Maxwell/VlasovMaxwell.hxx"

#ifdef MONTJOIE_WITH_TWO_DIM
#include "Hyperbolic/Maxwell/MontjoieTimeMaxwell2D.hxx"
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
#include "Hyperbolic/Maxwell/MontjoieTimeMaxwell3D.hxx"
#endif

#include "Hyperbolic/Maxwell/VlasovMaxwell.cxx"

#ifdef MONTJOIE_WITH_TWO_DIM
#include "Hyperbolic/Maxwell/2D/VlasovMaxwell2D.cxx"
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
#include "Hyperbolic/Maxwell/3D/VlasovMaxwell3D.cxx"
#endif

#define MONTJOIE_FILE_MONTJOIE_VLASOV_MAXWELL_HXX
#endif
