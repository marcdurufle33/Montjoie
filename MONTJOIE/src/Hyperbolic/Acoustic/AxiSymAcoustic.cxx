#ifndef MONTJOIE_FILE_AXISYM_ACOUSTIC_CXX

namespace Montjoie
{
  
  template<class TypeEquation>
  void TimeAcoustic_Axi<TypeEquation>
  ::SetInputData(const string& keyword, const Vector<string>& param)
  {
    VarInstationary<TypeEquation>::SetInputData(keyword, param);
  }

}

#define MONTJOIE_FILE_AXISYM_ACOUSTIC_CXX
#endif
