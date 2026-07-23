#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Hyperbolic/Elastic/MontjoieTimeElastic.hxx"
#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"

#include "Hyperbolic/Elastic/FluidStructureInteraction.hxx"
#include "Hyperbolic/Elastic/FluidStructureInteractionInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Elastic/FluidStructureInteraction.cxx"
#endif

using namespace Montjoie;

template<class TypeEqFluid, class TypeEqSolid>
void RunAll(FluidStructureInteraction<TypeEqFluid, TypeEqSolid>& var, const string& data_file,
	    const string& name_element)
{
  typedef Montjoie::Real_wp Real_wp;
  
  var.ConstructAll(data_file, name_element);
  
  // time parameters
  Real_wp t0 = var.GetInitialTime();
  Real_wp dt = var.GetTimeStep();
  Real_wp tf = var.GetFinalTime();  

  // then main loop in time
  int nb_max_iter;
  if (dt != Real_wp(0))
    nb_max_iter = toInteger(ceil(abs(tf - t0)/dt ) );
  else
    {
      cout << "Enter a time step different from 0 " << endl;
      abort();
    }
  
  int Np = var.var_fluid.GetNbScalarUnknowns();
  int Nv = var.var_fluid.GetNbVectorialUnknowns();
  int Nu = var.var_solid.GetNbScalarUnknowns();
  int Ns = var.var_solid.GetNbVectorialUnknowns();
  var.var_solid.var_harmonic.print_level = 2;
  
  if (var.var_fluid.GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_SYSTEM)
    {
      FluidStructureLeapFrogScheme scheme;
      scheme.SetInitialCondition(t0, dt, var);
      
      for (int n_time = 0; n_time < nb_max_iter; n_time++)
        {
          Real_wp t = t0 + n_time*dt;
          
          // we write snapshots if necessary
          var.var_fluid.WriteSnapshot(n_time, t, scheme.Pn);
          var.var_solid.WriteSnapshot(n_time, t, scheme.Un_half);
          
          // solution for (n+1) dt is computed
          scheme.Advance(t, n_time, var);      
        }
    }
  else
    {
      LowStorageRK_Iterator<Real_wp> scheme;
      scheme.SetOrder(4);
      
      VectReal_wp Y0(var.GetNbDof());
      Y0.Fill(0);
      scheme.SetInitialCondition(t0, dt, Y0, var);
      
      VectReal_wp Pn, Un, Vn, Sn;
      Pn.SetData(Np, &(scheme.GetIterate()(0)));
      Vn.SetData(Nv, &(scheme.GetIterate()(Np)));
      Un.SetData(Nu, &(scheme.GetIterate()(Np+Nv)));
      Sn.SetData(Ns, &(scheme.GetIterate()(Np+Nv+Nu)));
      
      for (int n_time = 0; n_time < nb_max_iter; n_time++)
        {
          Real_wp t = t0 + n_time*dt;
          
          // we write snapshots if necessary
          var.var_fluid.WriteVectorialSnapshot(n_time, t, Vn);
          var.var_fluid.WriteSnapshot(n_time, t, Pn);
          var.var_solid.WriteVectorialSnapshot(n_time, t, Sn);
          var.var_solid.WriteSnapshot(n_time, t, Un);
          
          // solution for (n+1) dt is computed
          scheme.Advance(t, n_time, var);      
        }
      
      Pn.Nullify();
      Un.Nullify();
      Vn.Nullify(); Sn.Nullify();
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  // premier argument -> fichier de donnees
  if (argc != 2)
    {
      cout << "Cette commande demande un argument " << endl;
      cout << "Par exemple, ./fluid_structure.x EXECUTION/data_file.ini " << endl;
      abort();
    }
  
  string name_data_file(argv[1]);
  
  // we get the type of element selected by the user, and type of equation
  string type_element, type_equation;
  getElement_Equation(name_data_file, type_element, type_equation);
  
  int dim_N = 3;
  if ((type_element.find("TRIANGLE") == 0) || (type_element.find("QUADRANGLE") == 0))
    dim_N = 2;
  
  // couplage entre les equations de l'elastodynamique et les equations de l'acoustique
  if (dim_N == 3)
    {
      // on lance la simulation
      if (type_equation == "HARMONIC_ELASTO_ACOUSTIC")
	{
	  FluidStructureInteractionHarmonic<HelmholtzEquation<Dimension3>, HarmonicElasticEquation<Dimension3> > var;
	  var.RunAll(name_data_file, type_element);
	}
      else
	{
	  FluidStructureInteraction<AcousticEquation<Dimension3>, TimeElasticEquation<Dimension3> > var;
	  RunAll(var, name_data_file, type_element);
	}
    }
  else
    {
      // on lance la simulation
      if (type_equation == "HARMONIC_ELASTO_ACOUSTIC")
	{
	  FluidStructureInteractionHarmonic<HelmholtzEquation<Dimension2>, HarmonicElasticEquation<Dimension2> > var;
	  var.RunAll(name_data_file, type_element);
	}
      else
	{
	  FluidStructureInteraction<AcousticEquation<Dimension2>, TimeElasticEquation<Dimension2> > var;
	  RunAll(var, name_data_file, type_element);
	}
    }
  
  return FinalizeMontjoie();
}

