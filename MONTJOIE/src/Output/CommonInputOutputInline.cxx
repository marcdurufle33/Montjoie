#ifndef MONTJOIE_FILE_COMMON_INPUT_OUTPUT_INLINE_CXX

namespace Montjoie
{

  //! adds value at the end of the file data_file (ascii format)
  inline void WriteOnTheGo(const string & data_file, Real_wp value)
  {
    ofstream file(data_file.data(), ios::app);
    file.precision(16);
    file << value << '\n';
    file.close();
  }
  
  
  //! adds a tiny vector of values at the end of the file data_file
  template<class T, int nb>
  inline void WriteOnTheGo(const string & data_file, TinyVector<T, nb> & value)
  {
    //cout << "Writing on the file " << endl;
    ofstream file(data_file.data(), ios::app);
    file.precision(16);
    PrintNoBrackets(file, value);
    file << '\n';
    file.close();
  }
  
  
  //! adds an array of values at the end of the file data_file
  template<class T>
  inline void WriteOnTheGo(const string & data_file, Vector<T> & value)
  {
    ofstream file(data_file.data(), ios::app);
    file.precision(16);
    for (int p = 0; p < value.GetM(); p++)
      file << value(p) << " ";
    
    file << '\n';
    file.close();
  }
  
  
}

#define MONTJOIE_FILE_COMMON_INPUT_OUTPUT_INLINE_CXX
#endif  

