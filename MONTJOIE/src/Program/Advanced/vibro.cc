#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_TRANSMISSION

#include "Hyperbolic/Elastic/MontjoieVibroacoustic.hxx"

using namespace Montjoie;

// main function
template<class TypeEqPlate, class TypeEqVol>
void RunAll(VibroAcousticSession<TypeEqPlate, TypeEqVol>& var, const string& data_file)
{
  HyperbolicProblem<TypeEqPlate>& var_plate = var.var_plate;
  HyperbolicProblem<TypeEqVol>& var_volume = var.var_volume;
  
  typedef Montjoie::Real_wp Real_wp;
  
  var.ConstructAll(data_file);
  
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
  
  // null initial condition
  if (var_volume.GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_ORDER2)
    {
      VibroAcousticLeapFrogScheme scheme;
      scheme.SetInitialCondition(t0, dt, var);
      if (var_volume.var_harmonic.GetNbProcPerMode() == 1)
        var_volume.var_harmonic.mesh.Write("test.mesh");
      
      for (int n_time = 0; n_time < nb_max_iter; n_time++)
        {
          Real_wp t = t0 + n_time*dt;
          
          // we write snapshots if necessary
          var_volume.WriteSnapshot(n_time, t, scheme.Ph_n);
          var_plate.WriteSnapshot(n_time, t, scheme.Lambda_n);
          
          // solution for (n+1) dt is computed
          scheme.Advance(t, n_time, var);      
        }
    }
  else if (var_volume.GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_SYSTEM)
    {
      VibroAcousticLeapFrogScheme scheme;
      scheme.first_order = true;
      scheme.SetInitialCondition(t0, dt, var);
      if (var_volume.var_harmonic.GetNbProcPerMode() == 1)
        var_volume.var_harmonic.mesh.Write("test.mesh");
      
      int Nvol = scheme.Ph_n.GetM();
      VectReal_wp PnPoint(Nvol);
      for (int n_time = 0; n_time < nb_max_iter; n_time++)
        {
          Real_wp t = t0 + n_time*dt;
                    
          /*Real_wp invDt = 1.0/dt;
          if (var_volume.var_harmonic.output_rcs_param.PointsOutside.GetM() > 0)
            for (int i = 0; i < Nvol; i++)
	    PnPoint(i) = var_fluide.Ph_nm1(i);*/
          
          // we write snapshots if necessary
          var_volume.WriteSnapshot(n_time, t, scheme.Ph_n);
          var_plate.WriteSnapshot(n_time, t, scheme.Lambda_n);
          
          // solution for (n+1) dt is computed
          scheme.Advance(t, n_time, var);      

          // computing solution on outside points if needed
          /*if (var_volume.var_harmonic.output_rcs_param.PointsOutside.GetM() > 0)
            {
              // on evalue PnPoint au temps n
              for (int i = 0; i < Nvol; i++)
                PnPoint(i) = (var_fluide.Ph_n(i) - PnPoint(i))*invDt;
              
              // on appelle le calcul du champ en des points exterieurs au domaine
              var_volume.var_harmonic.output_rcs_param.
                WriteOutput(n_time, t, dt, var_fluide.Ph_nm1, PnPoint);
		}*/
        }
    }
  else
    {
      RungeKutta_Iterator<Real_wp> RK;
      RK.SetOrder(4);
      
      int N = var.GetNbDof();
      VectReal_wp X0(N); X0.Fill(0);
      
      if (var_volume.var_harmonic.GetRankProcMode() == 0)
        cout << " starting time scheme" << endl;
      
      RunTimeScheme(t0, tf, dt, X0, var, RK);
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  // premier argument -> fichier de donnees
  if (argc != 2)
    {
      cout << "Cette commande demande un argument " << endl;
      cout << "Par exemple, ./vibro.x EXECUTION/data_file.ini " << endl;
      abort();
    }
  
  string name_data_file(argv[1]);

  // couplage entre les equations de Reissner-Mindlin et les equations de l'acoustique
  VibroAcousticSession<TimeReissnerMindlinEquation, AcousticEquation<Dimension3> > var;

  if (var.var_volume.var_harmonic.GetRankProcMode() == 0)  
    cout<<"******Application : VIBROACOUSTIQUE *******"<<endl<<endl;
    
  // on lance la simulation
  RunAll(var, name_data_file);
  
  return FinalizeMontjoie();
}
