#ifndef MONTJOIE_FILE_CROSS_CORRELATION_CXX

#include "CrossCorrelation/CrossCorrelation.hxx"
#include "CrossCorrelation/CrossCorrelationFunction.hxx"

namespace Montjoie
{

  CrossCorrelation::CrossCorrelation()
  {
    this->NumberModes = 0;
    this->Dimension = 0;
    this->CentralCrossCorrelation = 0;
    this->CrossCorrelationPoints.Resize(0);
    this->CrossCorrelationValues.Resize(4);
    this->CrossCorrelationValuesHomogeneous.Resize(1);
    this->DataValues.Resize(0);
    this->DataValuesHomogeneous.Resize(0);
    this->Time.Resize(0);
    this->IsOnFourier = false;
    this->IsHomogeneousDefined = false;
    this->IsComputedCrossCorrelation.Resize(4);

    for (int i=0; i<4; i++)
      this->IsComputedCrossCorrelation(i) = false;

    this->LightConfiguration = BACKLIGHT;
    this->ComputationBox.Resize(0);
    this->DiscretizationBox.Resize(0);
    this->NumberPictureModes = 0;
    this->PicturePoints.Resize(0);
    this->PictureValues.Resize(4);
    this->Velocity = 1.0;
    this->PrintLevel = 0;
    this->MapMaximumTime.Resize(0);
  }

  CrossCorrelation::~CrossCorrelation()
  {
    this->ComputationBox.Clear();
    this->DiscretizationBox.Clear();
    this->PicturePoints.Clear();
    this->PictureValues.Clear();
    this->MapMaximumTime.Clear();
  }

  void CrossCorrelation::InitFromFile(const string & DataFile)
  {
    ReadInputFile(DataFile, *this);
  }

  void CrossCorrelation::SetInputData(const string & keywords,
                                      const Vector<string> & parameters)
  {
    if (!keywords.compare("SismoPoint"))
      {
        int dimension = parameters.GetM();
        if (this->Dimension == 0)
          {
            this->Dimension = dimension;
            if (this->PrintLevel >= 2)
              cout << "Dimension of the problem: " << dimension << endl;
          }
        else
          {
            if (this->Dimension != dimension)
              {
                cerr << "Error: expecting point of dimension "
                     << this->Dimension << ", found one of dimension "
                     << dimension << endl;
                abort();
              }
          }
        R3 NewPoint;
        for (int index = 0; index < dimension; index++)
          {
            NewPoint(index) = to_num<Real_wp>(parameters(index));
          }
        this->CrossCorrelationPoints.PushBack(NewPoint);
        this->NumberModes++;
      }
    else if (!keywords.compare("FileOutputPoint"))
      {
        this->DataFileName = parameters(1);
      }
    else if (!keywords.compare("FileOutputPointHomogeneous"))
      {
        this->DataFileNameHomogeneous = parameters(1);
        this->IsHomogeneousDefined = true;
      }
    else if (!keywords.compare("LightConfiguration"))
      {
        if (!parameters(0).compare("Daylight"))
          {
            this->LightConfiguration = DAYLIGHT;
          }
        else if (!parameters(0).compare("Backlight"))
          {
            this->LightConfiguration = BACKLIGHT;
          }
        else if (!parameters(0).compare("Both"))
          {
            this->LightConfiguration = BOTH;
          }
        else
          {
            cerr << "Error: option " << parameters(0) << " for keyword "
                 << keywords << " not allowed. Abort." << endl;
            abort();
          }
      }
    else if (!keywords.compare("ComputationBox"))
      {
        this->ComputationBox.Resize(2*this->Dimension);
        for (int index = 0; index < 2*this->Dimension; index++)
          {
            this->ComputationBox(index) = to_num<Real_wp>(parameters(index));
          }
      }
    else if (!keywords.compare("CentralCrossCorrelation"))
      {
        this->CentralCrossCorrelation = to_num<int>(parameters(0));
      }
    else if (!keywords.compare("DiscretizationBox"))
      {
        this->DiscretizationBox.Resize(this->Dimension);
        for (int index = 0; index < this->Dimension; index++)
          {
            this->DiscretizationBox(index) = to_num<int>(parameters(index));
            this->DiscretizationBox(index)++;
          }
      }
    else if (!keywords.compare("PictureOutputFilename"))
      {
        this->OutputPictureFileName = parameters(0);
      }
    else if (!keywords.compare("PictureNumberModes"))
      {
        this->NumberPictureModes = to_num<int>(parameters(0));
      }
    else if (!keywords.compare("CrossCorrelationOutputFilename"))
      {
        this->OutputCrossFileName = parameters(0);
      }
    else if (!keywords.compare("HomogeneousVelocity"))
      {
        this->Velocity  = to_num<Real_wp>(parameters(0));
      }
    else if (!keywords.compare("BackgroundVelocity"))
      {
        this->Velocity  = to_num<Real_wp>(parameters(0));
      }
    else if (!keywords.compare("PrintLevel"))
      {
        this->PrintLevel = to_num<int>(parameters(0));
      }
  }

  int CrossCorrelation::GetNumberModes() const
  {
    return (this->NumberModes);
  }

  int CrossCorrelation::GetLightConfiguration() const
  {
    return (this->LightConfiguration);
  }

  void CrossCorrelation::SetLightConfiguration(int CurrentLightConfiguration)
  {
    this->LightConfiguration = CurrentLightConfiguration;
  }

  void CrossCorrelation::ReadDataFile()
  {

    this->DataValues.Resize(this->CrossCorrelationPoints.GetM());

    string line;

    VectString LineRead;

    int NumberLines = 0;
    int NumberLineRead = 0;

    ifstream FileIn(this->DataFileName.data());
    if (!FileIn.is_open())
      {
        cerr << "File " << this->DataFileName << " cannot be read. Abort..."
             << endl;
        abort();
      }

    while (!FileIn.eof())
      {
        getline(FileIn, line, '\n');

        // We count only non-empty lines
        if ((line.data()!=NULL)&&(strlen(line.data())!=0))
          NumberLines++;
      }
    FileIn.close();

    if (this->PrintLevel >= 2)
      {
        cout << "File " << this->DataFileName << " contains " << NumberLines
             << " lines." << endl;
      }

    this->Time.Resize(2*NumberLines);
    for(int index=0; index < this->NumberModes; index++)
      {
        this->DataValues(index).Resize(2*NumberLines);
      }

    FileIn.open(this->DataFileName.data());
    if (!FileIn.is_open())
      {
        cerr << "File " << this->DataFileName << " cannot be read. Abort..."
             << endl;
        abort();
      }

    // We fullfill second part of the array
    while (!FileIn.eof())
      {
        getline(FileIn, line, '\n');
        if ((line.data()!=NULL)&&(strlen(line.data())!=0))
          {
            StringTokenize(line, LineRead, " ");
            if (LineRead.GetM() != (this->NumberModes+1))
              {
                cerr << "Number of values read should be equal to "
                     << this->NumberModes+1 << ", but it's equal to "
                     << LineRead.GetM() << endl;
                abort();
              }
            this->Time(NumberLines + NumberLineRead)
              = to_num<Real_wp>(LineRead(0));
            for (int index=0; index < this->NumberModes; index++)
              {
                this->DataValues(index)(NumberLines + NumberLineRead)
                  = to_num<Complex_wp>(LineRead(index+1));
              }
            NumberLineRead++;
            if (this->PrintLevel >= 10)
              cout << "Number of lines read: " << NumberLineRead << endl;
          }
      }
    FileIn.close();

    // We compute DeltaTime to fullfill first part of the array
    Real_wp DeltaTime = this->Time(NumberLines+1)
      - this->Time(NumberLines);

    for(int tvalue = 0; tvalue < NumberLines; tvalue++)
      {
        this->Time(tvalue) = this->Time(NumberLines)
          - Real_wp(NumberLines - tvalue) * DeltaTime;
        for(int index = 0; index < this->NumberModes; index++)
          {
            this->DataValues(index)(tvalue) = 0.;
          }
      }

    // Finally, we substract "initial time" to get Time(NumberLineRead)=0
    Real_wp ShiftTime = this->Time(NumberLines);
    for(int tvalue = 0; tvalue < 2 * NumberLines; tvalue++)
      {
        this->Time(tvalue) -= ShiftTime;
      }
  }

  void CrossCorrelation::ReadDataFileForSecondOrder()
  {

    this->DataValues.Resize(this->CrossCorrelationPoints.GetM());

    string line;

    VectString LineRead;

    int NumberLines = 0;
    int NumberLineRead = 0;

    ifstream FileIn(this->DataFileName.data());
    if (!FileIn.is_open())
      {
        cerr << "File " << this->DataFileName << " cannot be read. Abort..."
             << endl;
        abort();
      }

    while (!FileIn.eof())
      {
        getline(FileIn, line, '\n');

        // We count only non-empty lines
        if ((line.data()!=NULL)&&(strlen(line.data())!=0))
          NumberLines++;
      }
    FileIn.close();

    if (this->PrintLevel >= 2)
      {
        cout << "File " << this->DataFileName << " contains " << NumberLines
             << " lines." << endl;
      }

    this->Time.Resize(4*NumberLines);
    for(int index=0; index < this->NumberModes; index++)
      {
        this->DataValues(index).Resize(4*NumberLines);
      }

    FileIn.open(this->DataFileName.data());
    if (!FileIn.is_open())
      {
        cerr << "File " << this->DataFileName << " cannot be read. Abort..."
             << endl;
        abort();
      }

    // We fullfill second part of the array
    while (!FileIn.eof())
      {
        getline(FileIn, line, '\n');
        if ((line.data()!=NULL)&&(strlen(line.data())!=0))
          {
            StringTokenize(line, LineRead, " ");
            if (LineRead.GetM() != (this->NumberModes+1))
              {
                cerr << "Number of values read should be equal to "
                     << this->NumberModes+1 << ", but it's equal to "
                     << LineRead.GetM() << endl;
                abort();
              }
            this->Time(2*NumberLines + NumberLineRead)
              = to_num<Real_wp>(LineRead(0));
            for (int index=0; index < this->NumberModes; index++)
              {
                this->DataValues(index)(2*NumberLines + NumberLineRead)
                  = to_num<Complex_wp>(LineRead(index+1));
              }
            NumberLineRead++;
            if (this->PrintLevel >= 10)
              cout << "Number of lines read: " << NumberLineRead << endl;
          }
      }
    FileIn.close();

    // We compute DeltaTime to fullfill other parts of the array
    Real_wp DeltaTime = this->Time(2*NumberLines+1)
      - this->Time(2*NumberLines);

    for(int tvalue = 0; tvalue < (2*NumberLines) ; tvalue++)
      {
        this->Time(tvalue) = this->Time(2*NumberLines)
          - Real_wp(2*NumberLines - tvalue) * DeltaTime;
        for(int index = 0; index < this->NumberModes; index++)
          {
            this->DataValues(index)(tvalue) = 0.;
          }
      }

    for(int tvalue = 0; tvalue < NumberLines ; tvalue++)
      {
        this->Time(3*NumberLines+tvalue) = this->Time(2*NumberLines)
          + Real_wp(NumberLines + tvalue) * DeltaTime;
        for(int index = 0; index < this->NumberModes; index++)
          {
            this->DataValues(index)(3*NumberLines+tvalue) = 0.;
          }
      }

    // Finally, we substract "initial time" to get Time(2*NumberLines)=0
    Real_wp ShiftTime = this->Time(2*NumberLines);
    for(int tvalue = 0; tvalue < 4*NumberLines; tvalue++)
      {
        this->Time(tvalue) -= ShiftTime;
      }
  }

  void CrossCorrelation::DispCrossCorrelationPoints() const
  {
    for (int index = 0; index < CrossCorrelationPoints.GetM(); index++)
      {
        cout << "Cross correlation point " << index << " with coordinates ("
             << CrossCorrelationPoints(index)(0) << ","
             << CrossCorrelationPoints(index)(1) << ","
             << CrossCorrelationPoints(index)(2) << ")" << endl;
      }
  }

  void CrossCorrelation::InitFourierInterface()
  {
    this->Fourier.Init(this->Time.GetM());
  }

  void CrossCorrelation::ApplyFourier()
  {
    if (this->IsOnFourier == true)
      {
        cerr << "Error: already set on Fourier" << endl;
        abort();
      }
    this->IsOnFourier = true;
    for(int index = 0; index < this->NumberModes ; index++)
      {
        this->Fourier.ApplyForward(this->DataValues(index));
      }
    for(int i=0; i<4; i++)
      {
        if (this->IsComputedCrossCorrelation(i) == true)
          {
            for (int IndexI = 0; IndexI < this->NumberPictureModes ; IndexI++)
              for (int IndexJ = 0; IndexJ < this->NumberPictureModes ; IndexJ++)
                this->Fourier.ApplyForward(this->CrossCorrelationValues(i)(IndexI)(IndexJ));
          }
      }
  }

  void CrossCorrelation::ApplyInverseFourier()
  {
    if (this->IsOnFourier == false)
      {
        cerr << "Error: already unset on Fourier" << endl;
        abort();
      }
    this->IsOnFourier = false;
    for(int index = 0; index < this->NumberModes ; index++)
      {
        this->Fourier.ApplyInverse(this->DataValues(index));
      }
    for(int i=0; i<4; i++)
      {
        if (this->IsComputedCrossCorrelation(i) == true)
          {
            for (int IndexI = 0; IndexI < this->NumberPictureModes ; IndexI++)
              for (int IndexJ = 0; IndexJ < this->NumberPictureModes ; IndexJ++)
                this->Fourier.ApplyInverse(this->CrossCorrelationValues(i)(IndexI)(IndexJ));
          }
      }
  }

  void CrossCorrelation::ComputeCrossCorrelations(int mode)
  {
    if (this->NumberPictureModes == 0)
      this->NumberPictureModes = this->NumberModes;

    if (this->IsComputedCrossCorrelation(mode-1) ==true)
      return;

    if (this->PrintLevel >= 2)
      cout << "Computation of cross correlation signal mode " << mode << endl;

    this->CrossCorrelationValues(mode-1).Resize(this->NumberPictureModes);

    for(int index=0; index< this->NumberPictureModes; index++)
      this->CrossCorrelationValues(mode-1)(index).Resize(this->NumberPictureModes);

    if (this->PrintLevel >= 2)
      cout << "Vector resized with " << this->NumberPictureModes << " modes" << endl;

    int DataSize = this->Time.GetM();

    if (this->IsOnFourier == true)
      {
        if (this->PrintLevel >= 2)
          {
            cout << "Computation of the cross-correlated modes in frequency"
                 << endl;
          }
        for (int IndexI = 0; IndexI < this->NumberPictureModes; IndexI++)
          for (int IndexJ = 0; IndexJ < this->NumberPictureModes; IndexJ++)
            {
              if (this->PrintLevel >= 4)
                {
                  cout << "Computation of the cross-correlated mode ("
                       << IndexI << "," << IndexJ << ")" << endl;
                }
              this->CrossCorrelationValues(mode-1)(IndexI)(IndexJ).Resize(DataSize);
              this->CrossCorrelationValues(mode-1)(IndexI)(IndexJ).Zero();

              // Computations of frequencies (// version)

              // Second order cross correlation
              if (mode == 1)
                {
                  ComputeFirstOrderCrossCorrelation(this->DataValues(IndexI),
                                                    this->DataValues(IndexJ),
                                                    this->CrossCorrelationValues(mode-1)
                                                    (IndexI)(IndexJ));

                }

            }
        this->IsComputedCrossCorrelation(mode-1) = true;
      }
    else
      {

        // Time domain computation
        if (this->PrintLevel >= 2)
          {
            cout << "Computation of the cross-correlated modes in time"
                 << endl;
          }
        for (int IndexI = 0; IndexI < this->NumberPictureModes; IndexI++)
          for (int IndexJ = 0; IndexJ < this->NumberPictureModes; IndexJ++)
            {
              if (this->PrintLevel >= 4)
                {
                  cout << "Computation of the cross-correlated mode ("
                       << IndexI << "," << IndexJ << ")" << endl;
                }
              this->CrossCorrelationValues(mode-1)(IndexI)(IndexJ).Resize(DataSize);
              this->CrossCorrelationValues(mode-1)(IndexI)(IndexJ).Zero();
              if (mode == 1)
                {
#pragma omp parallel for
                  for (int index = 0; index < DataSize; index++)
                    {
                      int lowerbound = max(0, DataSize / 2 - index);
                      int upperbound = min(DataSize, 3 * DataSize / 2 - index);
                      this->CrossCorrelationValues(0)(IndexI)(IndexJ)(index) = 0;
                      for(int indexd = lowerbound; indexd < upperbound; indexd++)
                        {
                          this->CrossCorrelationValues(0)(IndexI)(IndexJ)(index)
                            +=this->DataValues(IndexI)(indexd)
                            * this->DataValues(IndexJ)(indexd + index - DataSize/2);
                        }
                      this->CrossCorrelationValues(0)(IndexI)(IndexJ)(index)
                        *= 2. / Real_wp(DataSize);
                    }
                }

              if (mode == 2)
                {
                  // Even if defined in time, this will use Fourier computation
                  VectComplex_wp F,G,LF,LG,LPFG;
                  F.Resize(DataSize); G.Resize(DataSize);
                  LF.Resize(2*DataSize); LG.Resize(2*DataSize);
                  LPFG.Resize(2*DataSize); LPFG.Zero();
                  FftInterface<Complex_wp> LFourier;
                  LFourier.Init(2*DataSize);

                  for (int IndexK = 0; IndexK < this->NumberPictureModes; IndexK++)
                    {
                      F = CrossCorrelationValues(0)(IndexK)(IndexI);
                      G = CrossCorrelationValues(0)(IndexK)(IndexJ);
                      LF.Zero(); LG.Zero();
#pragma omp parallel for
                      for (int index = 0; index < DataSize; index++)
                        {
                          LF(DataSize / 2 + index) = F(index);
                          LG(DataSize / 2 + index) = G(index);
                        }
                      LFourier.ApplyForward(LF);
                      LFourier.ApplyForward(LG);
#pragma omp parallel for
                      for (int index = 0; index < 2*DataSize; index++)
                        {
                          // We test if index is even by looking weakest bit
                          if ((index & 1) == 0)
                            {
                              // index is even
                              LPFG(index) += LF(index) * LG(index);
                            }
                          else
                            {
                              // index is odd
                              LPFG(index) -= LF(index) * LG(index);
                            }
                        }
                    }
                  LFourier.ApplyInverse(LPFG);
#pragma omp parallel for
                  for (int index = 0; index < DataSize; index++)
                    {
                      this->CrossCorrelationValues(1)(IndexI)(IndexJ)(index)
                        += LPFG(DataSize / 2 + index);
                    }
                  F.Clear(); G.Clear(); LF.Clear(); LG.Clear(); LPFG.Clear();
                }

            }
        this->IsComputedCrossCorrelation(mode-1) = true;
      }
  }

  void CrossCorrelation::ExportCrossCorrelations(int mode)
  {
    if (IsComputedCrossCorrelation(mode-1) == false)
      {
        cerr << "Error while exporting cross correlations - not yet computed"
             << endl;
        abort();
      }
    if (this->IsOnFourier == true)
      {
        if (this->PrintLevel >= 2)
          cout << "Exporting in time-harmonic mode" << endl;
      }
    else
      {
        if (this->PrintLevel >= 2)
          cout << "Exporting in time-domain mode" << endl;
      }

    string CrossFileName = OutputCrossFileName;
    stringstream Smode; Smode << mode;
    CrossFileName += Smode.str();
    CrossFileName += ".dat";

    if (PrintLevel >= 2)
      cout << "Exporting file " << CrossFileName << endl;

    ofstream Output(CrossFileName.data());
    if (!Output.is_open())
      {
        cerr << "Output file " << this->OutputCrossFileName
             << " cannot be opened." << endl;
        abort();
      }

    int NumberCOmputations = this->Time.GetM();

    for(int Index = 0; Index < NumberCOmputations; Index++)
      {
        if (this->PrintLevel >= 10)
          {
            cout << "Writing line " << Index+1 << " out of "
                 << NumberCOmputations << endl;
          }
        if (this->IsOnFourier == true)
          {
            Output << 2 * pi_wp * Real_wp(Index) / Real_wp(NumberCOmputations);
          }
        else
          {
            Output << this->Time(Index);
          }
        for (int IndexI = 0; IndexI < this->NumberPictureModes; IndexI++)
          for (int IndexJ = 0; IndexJ < this->NumberPictureModes; IndexJ++)
            {
              Output << " "
                     << this->CrossCorrelationValues(mode-1)(IndexI)(IndexJ)(Index).real()
                     << " ";

              if (this->IsOnFourier == true)
                {
                  Output << this->CrossCorrelationValues(mode-1)(IndexI)(IndexJ)(Index).imag()
                         << " ";
                }
            }
        Output << endl;
      }
    Output.close();
  }


  void CrossCorrelation::ComputePicture(int mode)
  {
    if (this->NumberPictureModes == 0)
      this->NumberPictureModes = this->NumberModes;

    if (this->PrintLevel >= 2)
      {
        cout << "Computation of picture using " << this->NumberPictureModes
             << " modes" << endl;
      }

    int GlobalDimension = 1;
    for (int index=0; index < this->Dimension; index++)
      GlobalDimension *= this->DiscretizationBox(index);
    this->PictureValues(mode-1).Resize(GlobalDimension);

    this->PicturePoints.Resize(GlobalDimension);

    if (this->Dimension == 2)
      {
        this->ComputationBox.Resize(6);
        this->ComputationBox(4) = 0.0;
        this->ComputationBox(5) = 0.0;
        this->DiscretizationBox.Resize(3);
        this->DiscretizationBox(2) = 2;
      }

    Real_wp Xmin = ComputationBox(0), Ymin = ComputationBox(2),
      Zmin = ComputationBox(4);
    Real_wp DeltaX = (ComputationBox(1)-ComputationBox(0))
      /(Real_wp(DiscretizationBox(0)-1));
    Real_wp DeltaY = (ComputationBox(3)-ComputationBox(2))
      /(Real_wp(DiscretizationBox(1)-1));
    Real_wp DeltaZ = (ComputationBox(5)-ComputationBox(4))
      /(Real_wp(DiscretizationBox(2)-1));

    if (this->PrintLevel >= 2)
      cout << "Computation of output points" << endl;

    if (this->LightConfiguration == BACKLIGHT)
      {
        if (this->PrintLevel >= 2)
          cout << "Computation under backlight configuration" << endl;
        if (this->IsOnFourier == true)
          {
            if (this->PrintLevel >= 2)
              {
                cout << "KM imaging functionnal is computed in time domain"
                     << endl;
              }
            this->ApplyInverseFourier();
          }
      }

    if (this->LightConfiguration == DAYLIGHT)
      {
        if (this->PrintLevel >= 2)
          cout << "Computation under daylight configuration" << endl;
        if (this->IsOnFourier == true)
          {
            if (this->PrintLevel >= 2)
              {
                cout << "KM imaging functionnal is computed in time domain"
                     << endl;
              }
            this->ApplyInverseFourier();
          }
      }

#pragma omp parallel for
    for(int Index = 0; Index < GlobalDimension; Index++)
      {
        int IndexX = Index % this->DiscretizationBox(0);
        int IndexY = (Index / this->DiscretizationBox(0))
          % this->DiscretizationBox(1);
        int IndexZ = (Index / this->DiscretizationBox(0)
                      / this->DiscretizationBox(1))
          % this->DiscretizationBox(2);
        R3 NewPoint;
        NewPoint(0) = Xmin + Real_wp(IndexX) * DeltaX;
        NewPoint(1) = Ymin + Real_wp(IndexY) * DeltaY;
        NewPoint(2) = Zmin + Real_wp(IndexZ) * DeltaZ;
        this->PicturePoints(Index) = NewPoint;
        this->PictureValues(mode-1)(Index) = 0;

        // Computation of picture - backlight configuration
        if (this->LightConfiguration == BACKLIGHT)
          {
            for (int IndexJ = 0; IndexJ < this->NumberPictureModes; IndexJ++)
              {
                for (int IndexL = 0; IndexL < this->NumberPictureModes; IndexL++)
                  {
                    // Computation of distance between current point and xJ
                    Real_wp Distance_Point_To_J = 0.;
                    for (int IndexPoint = 0; IndexPoint < 3; IndexPoint++)
                      {
                        Distance_Point_To_J
                          += (NewPoint(IndexPoint)
                              - CrossCorrelationPoints(IndexJ)(IndexPoint))
                          * (NewPoint(IndexPoint)
                             - CrossCorrelationPoints(IndexJ)(IndexPoint));
                      }
                    Distance_Point_To_J = sqrt(Distance_Point_To_J);

                    // Computation of distance between current point and xL
                    Real_wp Distance_Point_To_L = 0.;
                    for (int IndexPoint = 0; IndexPoint < 3; IndexPoint++)
                      {
                        Distance_Point_To_L
                          += (NewPoint(IndexPoint)
                              - CrossCorrelationPoints(IndexL)(IndexPoint))
                          * (NewPoint(IndexPoint)
                             - CrossCorrelationPoints(IndexL)(IndexPoint));
                      }
                    Distance_Point_To_L = sqrt(Distance_Point_To_L);

                    // Computation of index K
                    Real_wp IndexK = Real_wp(this->Time.GetM())/2
                      + (Distance_Point_To_L - Distance_Point_To_J)
                      / (this->Velocity*(this->Time(1)-this->Time(0)));

                    if ((toInteger(IndexK) >= 0) && (toInteger(IndexK) < (Time.GetM()-1)))
                      {
                        this->PictureValues(mode-1)(Index)
                          += this->CrossCorrelationValues(mode-1)(IndexJ)(IndexL)(toInteger(IndexK));
                      }


                  }
              }
          }

        if (this->LightConfiguration == DAYLIGHT)
          {
            for (int IndexJ = 0; IndexJ < this->NumberPictureModes; IndexJ++)
              {
                for (int IndexL = 0; IndexL < this->NumberPictureModes; IndexL++)
                  {
                    // Computation of distance between current point and xJ
                    Real_wp Distance_Point_To_J = 0.;
                    for (int IndexPoint = 0; IndexPoint < 3; IndexPoint++)
                      {
                        Distance_Point_To_J
                          += (NewPoint(IndexPoint)
                              - CrossCorrelationPoints(IndexJ)(IndexPoint))
                          * (NewPoint(IndexPoint)
                             - CrossCorrelationPoints(IndexJ)(IndexPoint));
                      }
                    Distance_Point_To_J = sqrt(Distance_Point_To_J);

                    // Computation of distance between current point and xL
                    Real_wp Distance_Point_To_L = 0.;
                    for (int IndexPoint = 0; IndexPoint < 3; IndexPoint++)
                      {
                        Distance_Point_To_L
                          += (NewPoint(IndexPoint)
                              - CrossCorrelationPoints(IndexL)(IndexPoint))
                          * (NewPoint(IndexPoint)
                             - CrossCorrelationPoints(IndexL)(IndexPoint));
                      }
                    Distance_Point_To_L = sqrt(Distance_Point_To_L);

                    // Computation of index K
                    Real_wp IndexK = this->Time.GetM()/2
                      + (Distance_Point_To_L + Distance_Point_To_J)
                      / (this->Velocity*(this->Time(1)-this->Time(0)));


                    if ((IndexK >= 0) && (IndexK <= (Time.GetM()-1)))
                      {
                        this->PictureValues(mode-1)(Index)
                          += 2. * this->CrossCorrelationValues(mode-1)
                          (IndexJ)(IndexL)(toInteger(IndexK));
                      }


                  }
              }
          }

      }

  }

  void CrossCorrelation::ComputePictureFromSnapshots()
  {
    if (this->NumberPictureModes == 0)
      this->NumberPictureModes = this->NumberModes;

    if (this->PrintLevel >= 2)
      {
        cout << "Computation of picture using " << this->NumberPictureModes
             << " modes" << endl;
      }
  
    VectComplex_wp CurrentCrossCorrelation;
    int DataSize = this->Time.GetM();
    CurrentCrossCorrelation.Resize(DataSize);

    int GlobalDimension = 1;
    for (int index=0; index < this->Dimension; index++)
      GlobalDimension *= this->DiscretizationBox(index);
    this->PictureValues(0).Resize(GlobalDimension);

    this->PicturePoints.Resize(GlobalDimension);

    if (this->Dimension == 2)
      {
        this->ComputationBox.Resize(6);
        this->ComputationBox(4) = 0.0;
        this->ComputationBox(5) = 0.0;
        this->DiscretizationBox.Resize(3);
        this->DiscretizationBox(2) = 2;
      }
  
    Real_wp Xmin = ComputationBox(0), Ymin = ComputationBox(2),
      Zmin = ComputationBox(4);
    Real_wp DeltaX = (ComputationBox(1)-ComputationBox(0))
      /(Real_wp(DiscretizationBox(0)-1));
    Real_wp DeltaY = (ComputationBox(3)-ComputationBox(2))
      /(Real_wp(DiscretizationBox(1)-1));
    Real_wp DeltaZ = (ComputationBox(5)-ComputationBox(4))
      /(Real_wp(DiscretizationBox(2)-1));
  
    if (this->PrintLevel >= 2)
      cout << "Computation of output points" << endl;
  
    if (this->LightConfiguration == BACKLIGHT)
      {
        if (this->PrintLevel >= 2)
          cout << "Computation under backlight configuration" << endl;
      }
  
    if (this->LightConfiguration == DAYLIGHT)
      {
        if (this->PrintLevel >= 2)
          cout << "Computation under daylight configuration" << endl;
      }

    if (this->IsOnFourier == false)
      {
        if (this->PrintLevel >= 2)
          {
            cout << "KM imaging functionnal is computed in time domain,"
                 << "but we'll use time frequency functions"
                 << endl;
          }
        this->ApplyFourier();
      }
  
    // Start to prepare image
#pragma omp parallel for
    for(int Index = 0; Index < GlobalDimension; Index++)
      {
        int IndexX = Index % this->DiscretizationBox(0);
        int IndexY = (Index / this->DiscretizationBox(0))
          % this->DiscretizationBox(1);
        int IndexZ = (Index / this->DiscretizationBox(0)
                      / this->DiscretizationBox(1))
          % this->DiscretizationBox(2);
        R3 NewPoint;
        NewPoint(0) = Xmin + Real_wp(IndexX) * DeltaX;
        NewPoint(1) = Ymin + Real_wp(IndexY) * DeltaY;
        NewPoint(2) = Zmin + Real_wp(IndexZ) * DeltaZ;
        this->PicturePoints(Index) = NewPoint;
        this->PictureValues(0)(Index) = 0;
      }

    for (int IndexJ = 0; IndexJ < this->NumberPictureModes; IndexJ++)
      {
        for (int IndexL = 0; IndexL < this->NumberPictureModes; IndexL++)
          {

            if (this->PrintLevel >= 2)
              {
                cout << "Adding contribution of cross correlation C(t,"
                     << IndexJ+1 << "," << IndexL+1 << ")" << endl;
              }

            // We compute Cross Correlation under its Fourier Form
            ComputeFirstOrderCrossCorrelation(this->DataValues(IndexJ),
                                              this->DataValues(IndexJ),
                                              CurrentCrossCorrelation);
            this->Fourier.ApplyInverse(CurrentCrossCorrelation);
            // We add contribution of the Fourier form
#pragma omp parallel for
            for(int Index = 0; Index < GlobalDimension; Index++)
              {
                Real_wp Distance_Point_To_J = 0.;
                for (int IndexPoint = 0; IndexPoint < 3; IndexPoint++)
                  {
                    Distance_Point_To_J
                      += (PicturePoints(Index)(IndexPoint)
                          - CrossCorrelationPoints(IndexJ)(IndexPoint))
                      * (PicturePoints(Index)(IndexPoint)
                         - CrossCorrelationPoints(IndexJ)(IndexPoint));
                  }
                Distance_Point_To_J = sqrt(Distance_Point_To_J);

                // Computation of distance between current point and xL
                Real_wp Distance_Point_To_L = 0.;
                for (int IndexPoint = 0; IndexPoint < 3; IndexPoint++)
                  {
                    Distance_Point_To_L
                      += (PicturePoints(Index)(IndexPoint)
                          - CrossCorrelationPoints(IndexL)(IndexPoint))
                      * (PicturePoints(Index)(IndexPoint)
                         - CrossCorrelationPoints(IndexL)(IndexPoint));
                  }
                Distance_Point_To_L = sqrt(Distance_Point_To_L);
                Real_wp IndexK=0., Ponderation=1.;

                if (this->LightConfiguration == BACKLIGHT)
                  {
                    IndexK = this->Time.GetM()/2
                      + (Distance_Point_To_L - Distance_Point_To_J)
                      / (this->Velocity*(this->Time(1)-this->Time(0)));
                    Ponderation = 1.;
                  }
                if (this->LightConfiguration == DAYLIGHT)
                  {
                    IndexK = this->Time.GetM()/2
                      + (Distance_Point_To_L + Distance_Point_To_J)
                      / (this->Velocity*(this->Time(1)-this->Time(0)));
                    Ponderation = 2.;
                  }
                if ((IndexK >= 0) && (IndexK <= (Time.GetM()-1)))
                  {
                    this->PictureValues(0)(Index)
                      += Ponderation
                      * CurrentCrossCorrelation(toInteger(IndexK));
                  }

                    

              }
          }
      }
  
  
  }

  void CrossCorrelation::ComputePictureSecondOrderFromSnapshots()
  {
    if (this->NumberPictureModes == 0)
      this->NumberPictureModes = this->NumberModes;

    if (this->PrintLevel >= 2)
      {
        cout << "Computation of picture using " << this->NumberPictureModes
             << " modes" << endl;
      }
  
    VectComplex_wp CurrentCrossCorrelationJK, CurrentCrossCorrelationKL;
    VectComplex_wp CurrentCrossCorrelation;
    int DataSize = this->Time.GetM();
    CurrentCrossCorrelationJK.Resize(DataSize);
    CurrentCrossCorrelationKL.Resize(DataSize);
    CurrentCrossCorrelation.Resize(DataSize);

    Real_wp MaxDistance = 0;
    Real_wp DeltaTime = this->Time(DataSize/2+1) - this->Time(DataSize);
    Real_wp DistanceJL = 0;

    for (int IndexJ = 0; IndexJ < this->NumberPictureModes; IndexJ++)
      {
        for (int IndexL = 0; IndexL < this->NumberPictureModes; IndexL++)
          {
            DistanceJL = 0;
            for (int IndexPoint = 0; IndexPoint < 3; IndexPoint++)
              {
                DistanceJL
                  += (CrossCorrelationPoints(IndexL)(IndexPoint))
                  - CrossCorrelationPoints(IndexJ)(IndexPoint)
                  * (CrossCorrelationPoints(IndexL)(IndexPoint)
                     - CrossCorrelationPoints(IndexJ)(IndexPoint));
              }
            DistanceJL = sqrt(DistanceJL);
            if (DistanceJL > MaxDistance)
              {
                MaxDistance = DistanceJL;
              }
          }
      }
    
    int IndexTruncate = toInteger(MaxDistance / (this->Velocity * DeltaTime));


    int GlobalDimension = 1;
    for (int index=0; index < this->Dimension; index++)
      GlobalDimension *= this->DiscretizationBox(index);
    this->PictureValues(1).Resize(GlobalDimension);

    this->PicturePoints.Resize(GlobalDimension);

    if (this->Dimension == 2)
      {
        this->ComputationBox.Resize(6);
        this->ComputationBox(4) = 0.0;
        this->ComputationBox(5) = 0.0;
        this->DiscretizationBox.Resize(3);
        this->DiscretizationBox(2) = 2;
      }
  
    Real_wp Xmin = ComputationBox(0), Ymin = ComputationBox(2),
      Zmin = ComputationBox(4);
    Real_wp DeltaX = (ComputationBox(1)-ComputationBox(0))
      /(Real_wp(DiscretizationBox(0)-1));
    Real_wp DeltaY = (ComputationBox(3)-ComputationBox(2))
      /(Real_wp(DiscretizationBox(1)-1));
    Real_wp DeltaZ = (ComputationBox(5)-ComputationBox(4))
      /(Real_wp(DiscretizationBox(2)-1));
  
    if (this->PrintLevel >= 2)
      cout << "Computation of output points" << endl;
  
    if (this->LightConfiguration == BACKLIGHT)
      {
        if (this->PrintLevel >= 2)
          cout << "Computation under backlight configuration" << endl;
      }
  
    if (this->LightConfiguration == DAYLIGHT)
      {
        if (this->PrintLevel >= 2)
          cout << "Computation under daylight configuration" << endl;
      }

    if (this->IsOnFourier == false)
      {
        if (this->PrintLevel >= 2)
          {
            cout << "KM imaging functionnal is computed in time domain,"
                 << "but we'll use time frequency functions"
                 << endl;
          }
        this->ApplyFourier();
      }
  
    // Start to prepare image
#pragma omp parallel for
    for(int Index = 0; Index < GlobalDimension; Index++)
      {
        int IndexX = Index % this->DiscretizationBox(0);
        int IndexY = (Index / this->DiscretizationBox(0))
          % this->DiscretizationBox(1);
        int IndexZ = (Index / this->DiscretizationBox(0)
                      / this->DiscretizationBox(1))
          % this->DiscretizationBox(2);
        R3 NewPoint;
        NewPoint(0) = Xmin + Real_wp(IndexX) * DeltaX;
        NewPoint(1) = Ymin + Real_wp(IndexY) * DeltaY;
        NewPoint(2) = Zmin + Real_wp(IndexZ) * DeltaZ;
        this->PicturePoints(Index) = NewPoint;
        this->PictureValues(1)(Index) = 0;
      }

    for (int IndexJ = 0; IndexJ < this->NumberPictureModes; IndexJ++)
      {
        for (int IndexL = 0; IndexL < this->NumberPictureModes; IndexL++)
          {

            if (this->PrintLevel >= 2)
              {
                cout << "Adding contribution of cross correlation C(t,"
                     << IndexJ+1 << "," << IndexL+1 << ")" << endl;
              }
#pragma omp parallel for
            for (int IndexOmega = 0; IndexOmega < DataSize; IndexOmega++)
              {
                CurrentCrossCorrelation(IndexOmega)=0;
              }
            for (int IndexK = 0; IndexK < this->NumberPictureModes; IndexK++)
              {
                if (this->PrintLevel >= 6)
                  {
                    cout << "Adding cross contribution of cross correlations "
                         << "C(" << IndexJ+1 << "," << IndexK+1 << ") and"
                         << "C(" << IndexK+1 << "," << IndexL+1 << ")"
                         << endl;
                  }
                // Computation of Cross Correlation J-K
                ComputeFirstOrderCrossCorrelation(this->DataValues(IndexJ),
                                                  this->DataValues(IndexK),
                                                  CurrentCrossCorrelationJK);
                // Computation of Cross Correlation K-L
                ComputeFirstOrderCrossCorrelation(this->DataValues(IndexK),
                                                  this->DataValues(IndexL),
                                                  CurrentCrossCorrelationKL);
                // Apply inverse Fourier transform
                this->Fourier.ApplyInverse(CurrentCrossCorrelationJK);
                this->Fourier.ApplyInverse(CurrentCrossCorrelationKL);
                // Truncate indexes
#pragma omp parallel for
                for (int IndexTime = DataSize/2 - IndexTruncate;
                     IndexTime < DataSize/2 + IndexTruncate;
                     IndexTime++)
                  {
                    CurrentCrossCorrelationJK(IndexTime)=0;
                    CurrentCrossCorrelationKL(IndexTime)=0;
                  }
                // Apply Forward Fourier transform
                this->Fourier.ApplyForward(CurrentCrossCorrelationJK);
                this->Fourier.ApplyForward(CurrentCrossCorrelationKL);

                // We compute Second order Cross Correlation
                ComputeFirstOrderCrossCorrelation(CurrentCrossCorrelationJK,
                                                  CurrentCrossCorrelationKL,
                                                  CurrentCrossCorrelation);
                this->Fourier.ApplyInverse(CurrentCrossCorrelation);
                
              }
            

            // We add contribution of the Fourier form
#pragma omp parallel for
            for(int Index = 0; Index < GlobalDimension; Index++)
              {
                Real_wp Distance_Point_To_J = 0.;
                for (int IndexPoint = 0; IndexPoint < 3; IndexPoint++)
                  {
                    Distance_Point_To_J
                      += (PicturePoints(Index)(IndexPoint)
                          - CrossCorrelationPoints(IndexJ)(IndexPoint))
                      * (PicturePoints(Index)(IndexPoint)
                         - CrossCorrelationPoints(IndexJ)(IndexPoint));
                  }
                Distance_Point_To_J = sqrt(Distance_Point_To_J);

                // Computation of distance between current point and xL
                Real_wp Distance_Point_To_L = 0.;
                for (int IndexPoint = 0; IndexPoint < 3; IndexPoint++)
                  {
                    Distance_Point_To_L
                      += (PicturePoints(Index)(IndexPoint)
                          - CrossCorrelationPoints(IndexL)(IndexPoint))
                      * (PicturePoints(Index)(IndexPoint)
                         - CrossCorrelationPoints(IndexL)(IndexPoint));
                  }
                Distance_Point_To_L = sqrt(Distance_Point_To_L);
                Real_wp IndexK=0., Ponderation=1.;

                if (this->LightConfiguration == BACKLIGHT)
                  {
                    IndexK = this->Time.GetM()/2
                      + (Distance_Point_To_L - Distance_Point_To_J)
                      / (this->Velocity*(this->Time(1)-this->Time(0)));
                    Ponderation = 1.;
                  }
                if (this->LightConfiguration == DAYLIGHT)
                  {
                    IndexK = this->Time.GetM()/2
                      + (Distance_Point_To_L + Distance_Point_To_J)
                      / (this->Velocity*(this->Time(1)-this->Time(0)));
                    Ponderation = 2.;
                  }
                if ((IndexK >= 0) && (IndexK <= (Time.GetM()-1)))
                  {
                    this->PictureValues(1)(Index)
                      += Ponderation
                      * CurrentCrossCorrelation(toInteger(IndexK));
                  }

                    

              }
          }
      }
  
  
  }

  void CrossCorrelation::ExportPicture(int mode)
  {

    string PictureFileName = this->OutputPictureFileName;
    stringstream Smode; Smode << mode;
    PictureFileName += Smode.str();

    if (this->LightConfiguration == BACKLIGHT)
      {
        PictureFileName += "B.dat";
      }

    if (this->LightConfiguration == DAYLIGHT)
      {
        PictureFileName += "D.dat";
      }

    if (this->PrintLevel >= 2)
      cout << "Exportation of picture " << PictureFileName << endl;


    ofstream Output(PictureFileName.data());
    if (!Output.is_open())
      {
        cerr << "Output file " << this->OutputPictureFileName
             << " cannot be opened." << endl;
        abort();
      }

    int GlobalDimension = 1;
    for (int index=0; index < this->Dimension; index++)
      GlobalDimension *= this->DiscretizationBox(index);

    for(int Index = 0; Index < GlobalDimension; Index++)
      {
        if (this->PrintLevel >= 10)
          {
            cout << "Writing line " << Index+1 << " out of "
                 << GlobalDimension << endl;
          }
        Output << this->PicturePoints(Index)(0) << " "
               << this->PicturePoints(Index)(1) << " "
               << this->PicturePoints(Index)(2) << " "
               << this->PictureValues(mode-1)(Index).real()
               << endl;
      }
    Output.close();
  }

  void CrossCorrelation::CreateMapMaximumTime()
  {
    this->CreateMapMaximumTime(this->CentralCrossCorrelation);
  }

  void CrossCorrelation::CreateMapMaximumTime(const int Central)
  {
    const int NumberCrossCorrelationPoints 
      = this->CrossCorrelationPoints.GetM();
    this->MapMaximumTime.Resize(NumberCrossCorrelationPoints);
#pragma omp parallel for
    for (int Index = 0; Index < NumberCrossCorrelationPoints; Index++)
      {
        VectComplex_wp CurrentCrossCorrelation;
        VectReal_wp Data;
        ComputeFirstOrderCrossCorrelation(this->DataValues(Index),
                                          this->DataValues(Central),
                                          CurrentCrossCorrelation);

        this->Fourier.ApplyInverse(CurrentCrossCorrelation);
        Data.Resize(CurrentCrossCorrelation.GetM());
#pragma omp parallel for
        for(int I = 0; I < CurrentCrossCorrelation.GetM(); I++)
          {
            Data(I) = real(CurrentCrossCorrelation(I));
          }

        this->MapMaximumTime(Index)
          = ReturnArgmin(this->Time, Data);
        CurrentCrossCorrelation.Clear();
        Data.Clear();
      }
  }

  void CrossCorrelation::CreateMapMaximumTimeEnvelop()
  {
    this->CreateMapMaximumTimeEnvelop(this->CentralCrossCorrelation);
  }

  void CrossCorrelation::CreateMapMaximumTimeEnvelop(const int Central)
  {
    const int NumberCrossCorrelationPoints 
      = this->CrossCorrelationPoints.GetM();
    this->MapMaximumTime.Resize(NumberCrossCorrelationPoints);
#pragma omp parallel for
    for (int Index = 0; Index < NumberCrossCorrelationPoints; Index++)
      {
        VectComplex_wp CurrentCrossCorrelation;
        VectReal_wp Data;
        ComputeFirstOrderCrossCorrelation(this->DataValues(Index),
                                          this->DataValues(Central),
                                          CurrentCrossCorrelation);

        ComputeEnvelopFourier(CurrentCrossCorrelation, Data);

        this->MapMaximumTime(Index)
          = ReturnArgmin(this->Time, Data);
        CurrentCrossCorrelation.Clear();
        Data.Clear();
      }
  }

  void CrossCorrelation::CreateMapSpeed()
  {
    int Size = this->CrossCorrelationPoints.GetM();
    int nx, ny, nz = Size;
    for (nx = 0; (nx < nz) &&
          (CrossCorrelationPoints(nx)(1) == CrossCorrelationPoints(0)(1)); nx++);
    for (ny = 0; (ny < nz) &&
           (CrossCorrelationPoints(nx)(2) == CrossCorrelationPoints(0)(2)); ny++);
    nz /= ny;
    ny /= nx;
    VectComplex_wp GradientX, GradientY, GradientZ;
    GradientX.Resize(Size);
    GradientY.Resize(Size);
    GradientZ.Resize(Size);

    this->MapSpeed.Resize(Size);

    // Calcul de la composante en X du gradient

#pragma omp parallel for
    for(int Index = 0; Index < Size; Index++)
      {
        int IdX = Index % nx;

        if ( (IdX != 0) && (IdX != (nx-1) ) )
          {
            GradientX(Index) = (MapMaximumTime(Index+1) - MapMaximumTime(Index-1)) 
              / ( CrossCorrelationPoints(Index+1)(0) - CrossCorrelationPoints(Index-1)(0) );
          }

        if ( IdX == 0 )
          {
            GradientX(Index) = (MapMaximumTime(Index+1) - MapMaximumTime(Index)) 
              / ( CrossCorrelationPoints(Index+1)(0) - CrossCorrelationPoints(Index)(0) );
          }

        if ( IdX == (nx-1) )
          {
            GradientX(Index) = (MapMaximumTime(Index) - MapMaximumTime(Index-1))
              / ( CrossCorrelationPoints(Index)(0) - CrossCorrelationPoints(Index-1)(0) );
          }

      }

    // Calcul de la composante en Y du gradient

    if (ny == 1)
      {
#pragma omp parallel for
        for(int Index = 0; Index < Size; Index++)
          {
            GradientY(Index) = 0;
          }
      }
    else
      {
#pragma omp parallel for
        for(int Index = 0; Index < Size; Index++)
          {
            int IdY = (Index / nx) % ny;

            if ( (IdY != 0) && (IdY != (ny-1) ) )
              {
                GradientY(Index) = (MapMaximumTime(Index+nx) - MapMaximumTime(Index-nx)) 
                  / ( CrossCorrelationPoints(Index+nx)(1) - CrossCorrelationPoints(Index-nx)(1) );
              }
            
            if ( IdY == 0 )
              {
                GradientY(Index) = (MapMaximumTime(Index+nx) - MapMaximumTime(Index))
                  / ( CrossCorrelationPoints(Index+nx)(1) - CrossCorrelationPoints(Index)(1) );
              }

            if ( IdY == (ny-1) )
              {
                GradientY(Index) = (MapMaximumTime(Index) - MapMaximumTime(Index-nx))
                  / ( CrossCorrelationPoints(Index)(1) - CrossCorrelationPoints(Index-nx)(1) );
              }
            
          }
      }

    // Calcul de la composante en Z du gradient

    if (nz == 1)
      {
#pragma omp parallel for
        for(int Index = 0; Index < Size; Index++)
          {
            GradientZ(Index) = 0;
          }
      }
    else
      {
#pragma omp parallel for
        for(int Index = 0; Index < Size; Index++)
          {
            int IdZ = Index / (nx * ny);

            if ( (IdZ != 0) && (IdZ != (nz-1) ) )
              {
                GradientZ(Index) = (MapMaximumTime(Index+(nx*ny)) - MapMaximumTime(Index-(nx*ny)))
                  / ( CrossCorrelationPoints(Index+(nx*ny))(2)
                      - CrossCorrelationPoints(Index-(nx*ny))(2) );
              }
            
            if ( IdZ == 0 )
              {
                GradientZ(Index) = (MapMaximumTime(Index+(nx*ny)) - MapMaximumTime(Index))
                  / ( CrossCorrelationPoints(Index+(nx*ny))(2)
                      - CrossCorrelationPoints(Index)(2) );
              }

            if ( IdZ == (nz-1) )
              {
                GradientZ(Index) = (MapMaximumTime(Index) - MapMaximumTime(Index-(nx*ny))) /
                  ( CrossCorrelationPoints(Index)(2) - CrossCorrelationPoints(Index-(nx*ny))(2) );
              }
            
          }
      }
    
    // Calcul de la norme du gradient
#pragma omp parallel for
    for(int Index = 0; Index < Size; Index++)
      {
        MapSpeed(Index) = sqrt(abs(GradientX(Index)) * abs(GradientX(Index))
                               + abs(GradientY(Index)) * abs(GradientY(Index))
                               + abs(GradientZ(Index)) * abs(GradientZ(Index)));
      }
    
  }

  void CrossCorrelation::CreateMapAverageSpeed()
  {
    int Size = this->CrossCorrelationPoints.GetM();
    Vector<Real_wp> AverageSpeed;
    AverageSpeed.Resize(Size);
    AverageSpeed.Zero();
    for(int Index=0; Index < Size; Index++)
      {
        if (this->PrintLevel >= 2)
          {
            cout << "Treating point " << Index+1 << " out of "
                 << Size << endl;
          }
        this->CreateMapMaximumTime(Index);
        this->CreateMapSpeed();
#pragma omp parallel for
        for(int IdPoint=0; IdPoint < Size; IdPoint++)
          {
            AverageSpeed(IdPoint) += this->MapSpeed(IdPoint);
          }
      }
#pragma omp parallel for
        for(int IdPoint=0; IdPoint < Size; IdPoint++)
          {
            this->MapSpeed(IdPoint) = AverageSpeed(IdPoint) / Real_wp(Size);
          }

    AverageSpeed.Clear();

  }

  void CrossCorrelation::ExportMapTime()
  {
    ofstream Output("Time.dat");
    if (!Output.is_open())
      {
        cerr << "Output file " << this->OutputPictureFileName
             << " cannot be opened." << endl;
        abort();
      }
    int Size = this->CrossCorrelationPoints.GetM();
    for(int Index = 0; Index < Size; Index++)
      {
        if (this->PrintLevel >= 6)
          {
            cout << "Writing line " << Index+1 << " out of "
                 << Size << endl;
          }
        Output << this->CrossCorrelationPoints(Index)(0) << " "
               << this->CrossCorrelationPoints(Index)(1) << " "
               << this->CrossCorrelationPoints(Index)(2) << " "
               << this->MapMaximumTime(Index)
               << endl;
      }
    Output.close();
  }

  void CrossCorrelation::ExportMapSpeed()
  {
    ofstream Output("Speed.dat");
    if (!Output.is_open())
      {
        cerr << "Output file " << this->OutputPictureFileName
             << " cannot be opened." << endl;
        abort();
      }
    int Size = this->CrossCorrelationPoints.GetM();
    for(int Index = 0; Index < Size; Index++)
      {
        if (this->PrintLevel >= 6)
          {
            cout << "Writing line " << Index+1 << " out of "
                 << Size << endl;
          }
        Output << this->CrossCorrelationPoints(Index)(0) << " "
               << this->CrossCorrelationPoints(Index)(1) << " "
               << this->CrossCorrelationPoints(Index)(2) << " "
               << this->MapSpeed(Index)
               << endl;
      }
    Output.close();
  }

}

#define MONTJOIE_FILE_CROSS_CORRELATION_CXX
#endif
