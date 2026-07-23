#ifndef MONTJOIE_FILE_CROSS_CORRELATION_HXX

enum {BACKLIGHT, DAYLIGHT, BOTH};

namespace Montjoie
{

  class CrossCorrelation : public InputDataProblem_Base
  {
  protected:
    int NumberModes;
    int Dimension;
    int CentralCrossCorrelation;
    string DataFileName;
    string DataFileNameHomogeneous;
    string OutputPictureFileName;
    string OutputCrossFileName;
    VectR3 CrossCorrelationPoints;
    Vector<Vector<Vector<VectComplex_wp> > > CrossCorrelationValues;
    Vector<Vector<Vector<VectComplex_wp> > > CrossCorrelationValuesHomogeneous;
    Vector<VectComplex_wp> DataValues;
    Vector<VectComplex_wp> DataValuesHomogeneous;
    Vector<Real_wp> Time;
    FftInterface<Complex_wp> Fourier;
    bool IsOnFourier;
    bool IsHomogeneousDefined;
    Vector<bool> IsComputedCrossCorrelation;
    int LightConfiguration;
    Vector<Real_wp> ComputationBox;
    Vector<int> DiscretizationBox;
    int NumberPictureModes;
    VectR3 PicturePoints;
    Vector<VectComplex_wp> PictureValues;
    Real_wp Velocity;
    int PrintLevel;
    Vector<Real_wp> MapMaximumTime;
    Vector<Real_wp> MapSpeed;
  public:

    CrossCorrelation();
    ~CrossCorrelation();
    int GetNumberModes() const;
    void InitFromFile(const string &);
    void SetInputData(const string &, const Vector<string> &);
    void ReadDataFile();
    void ReadDataFileForSecondOrder();
    int GetLightConfiguration() const;
    void SetLightConfiguration(int);
    void DispCrossCorrelationPoints() const;
    void InitFourierInterface();
    void ApplyFourier();
    void ApplyInverseFourier();
    void ComputeCrossCorrelations(int);
    void ExportCrossCorrelations(int);
    void ComputePicture(int);
    void ComputePictureFromSnapshots();
    void ComputePictureSecondOrderFromSnapshots();
    void ExportPicture(int);
    void CreateMapMaximumTime();
    void CreateMapMaximumTime(const int);
    void CreateMapMaximumTimeEnvelop();
    void CreateMapMaximumTimeEnvelop(const int);
    void ExportMapTime();
    void CreateMapSpeed();
    void CreateMapAverageSpeed();
    void ExportMapSpeed();
  };

}

#define MONTJOIE_FILE_CROSS_CORRELATION_HXX
#endif
