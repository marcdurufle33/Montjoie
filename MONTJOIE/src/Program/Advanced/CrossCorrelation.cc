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
  CrossCorrelationProblem.ComputeCrossCorrelations(1);
  CrossCorrelationProblem.ApplyInverseFourier();
  CrossCorrelationProblem.ExportCrossCorrelations(1);
  if (CrossCorrelationProblem.GetLightConfiguration() == BOTH)
    {
      for (int CurrentMode = 1; CurrentMode <= 1 ; CurrentMode++)
        {
          CrossCorrelationProblem.SetLightConfiguration(BACKLIGHT);
          CrossCorrelationProblem.ComputePicture(CurrentMode);
          CrossCorrelationProblem.ExportPicture(CurrentMode);
          CrossCorrelationProblem.SetLightConfiguration(DAYLIGHT);
          CrossCorrelationProblem.ComputePicture(CurrentMode);
          CrossCorrelationProblem.ExportPicture(CurrentMode);
        }
    }
  else
    {
      for (int CurrentMode = 1; CurrentMode <= 1 ; CurrentMode++)
        {
          CrossCorrelationProblem.ComputePicture(CurrentMode);
          CrossCorrelationProblem.ExportPicture(CurrentMode);
        }
    }

  FinalizeMontjoie();
  return(0);
	
}
