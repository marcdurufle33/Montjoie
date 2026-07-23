#define MONTJOIE_WITH_TWO_DIM

#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  if (argc < 3)
    {
      cout << "Argument missing" << endl
           << "AppendFile.x resultfile file1 file2 ... is a good syntax"
           << endl << "Aborting..." << endl;
      abort();
    }

  string line;
  int NumberLines = 0, NumberColumns = 0;

  ifstream FileIn;

  // Detect number of lines
  FileIn.open(argv[2]);
  if (!FileIn.is_open())
    {
      cerr << "File " << argv[2] << " cannot be read. Abort..." << endl;
      abort();
    }

  cout << "Reading once " << argv[2] << " to determine number of lines and"
       << " number of columns." << endl;

  while (!FileIn.eof())
    {
      getline(FileIn, line, '\n');
      // We consider that a non-empty line contains datas
      if ( (line.data() != NULL) && (strlen(line.data()) != 0) )
        {
          NumberLines++;
          // We count once number of columns
          if (NumberColumns == 0)
            {
              Vector<string> LineRead;
              StringTokenize(line, LineRead, " ");
              NumberColumns = LineRead.GetM();
            }
        }
    }
  
  FileIn.close();

  cout << "File contains " << NumberLines << " lines and " << NumberColumns
       << " columns." << endl;


  Vector<Vector<Real_wp> > Data;
  Data.Resize(NumberLines);
  for(int index=0; index < NumberLines; index++)
    {
      Data(index).Resize(NumberColumns);
      Data(index).Zero();
    }

  for(int currentfile=2; currentfile<argc; currentfile++)
    {
      FileIn.open(argv[currentfile]);
      if (!FileIn.is_open())
        {
          cerr << "File " << argv[currentfile] << " cannot be read. Abort..."
               << endl;
          abort();
        }
      
      cout << "Reading " << argv[currentfile] << "." << endl;
      
      int CurrentLine = 0;

      while (!FileIn.eof())
        {
          getline(FileIn, line, '\n');
          // We consider that a non-empty line contains datas
          if ( (line.data() != NULL) && (strlen(line.data()) != 0) )
            {
              Vector<string> LineRead;
              StringTokenize(line, LineRead, " ");
              if ( LineRead.GetM() != NumberColumns )
                {
                  cerr << "Number of values read should be equal to "
                       << NumberColumns << ", but it's equal to"
                       << LineRead.GetM() << endl;
                  abort();
                }
              Data(CurrentLine)(0) = to_num<Real_wp>(LineRead(0));
              for (int Column = 1; Column < NumberColumns; Column++)
                Data(CurrentLine)(Column) += 
                  to_num<Real_wp>(LineRead(Column));
            }
          CurrentLine++;
        }
  
      FileIn.close();
      
    }

  ofstream FileOut;

  // Detect number of lines
  FileOut.open(argv[1]);
  if (!FileOut.is_open())
    {
      cerr << "File " << argv[1] << " cannot be written. Abort..." << endl;
      abort();
    }
  
  for(int Line = 0; Line < NumberLines; Line++)
    {
      for(int Column = 0; Column < NumberColumns; Column++)
        {
          FileOut << Data(Line)(Column) << " ";
        }
      FileOut << endl;
    }
  FileOut.close();

  return(0);
	
}
