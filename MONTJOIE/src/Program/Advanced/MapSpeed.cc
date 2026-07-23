#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "CrossCorrelation/MontjoieCrossCorrelation.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      cout << "Argument missing" << endl
           << "crosscorrelation.x file is a good syntax" << endl
           << "Aborting..." << endl;
      abort();
    }

  InitMontjoie(argc,argv);

  CrossCorrelation CrossCorrelationProblem;
  CrossCorrelationProblem.InitFromFile(string(argv[1]));
  CrossCorrelationProblem.ReadDataFile();
  CrossCorrelationProblem.InitFourierInterface();
  CrossCorrelationProblem.ApplyFourier();
  CrossCorrelationProblem.CreateMapMaximumTime();
  CrossCorrelationProblem.ExportMapTime();
  CrossCorrelationProblem.CreateMapSpeed();
  //  CrossCorrelationProblem.CreateMapAverageSpeed();
  CrossCorrelationProblem.ExportMapSpeed();
  FinalizeMontjoie();
  return(0);
	
}
