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
  if (CrossCorrelationProblem.GetLightConfiguration() == BOTH)
    {
      CrossCorrelationProblem.SetLightConfiguration(BACKLIGHT);
      CrossCorrelationProblem.ComputePictureFromSnapshots();
      CrossCorrelationProblem.ExportPicture(1);
      CrossCorrelationProblem.SetLightConfiguration(DAYLIGHT);
      CrossCorrelationProblem.ComputePictureFromSnapshots();
      CrossCorrelationProblem.ExportPicture(1);
    }
  else
    {
      CrossCorrelationProblem.ComputePictureFromSnapshots();
      CrossCorrelationProblem.ExportPicture(1);
    }

  FinalizeMontjoie();
  return(0);

}
