#ifndef MONTJOIE_FILE_HELMHOLTZ_RADIAL_INLINE_CXX

namespace Montjoie
{

  template<class Complexe>
  inline int VarHelmholtz_Radial<Complexe>::GetLmax() const
  {
    return Lmax;
  }
  

  template<class Complexe>
  inline int VarHelmholtz_Radial<Complexe>::GetIndexOutputFiles(int type)
  {
    int index = -1;
    for (int i  = 0; i < this->output_grid_param.GetM(); i++)
      if (this->output_grid_param(i).GetNature() == type)
        index = i;
    
    if (index == -1)
      {
        int N = this->output_grid_param.GetM();
        this->output_grid_param.Resize(N+1);
        this->output_grid_param(N).SetNature(type);
        index = N;
      }
    
    return index;
  }
  
  
  template<class Complexe>
  inline Complex_wp VarHelmholtz_Radial<Complexe>
  ::GetYnm(int l, int m2, const Vector<VectReal_wp>& P,
           const Complex_wp& phase)
  {
    if (m2 == 0)
      return Complex_wp(P(l)(0), 0);
    
    int m = (m2+1)/2;
    if (m2%2 == 1)
      {
        if (m%2 == 1)
          return -P(l)(m)*conj(phase);
        else
          return P(l)(m)*conj(phase);
      }
    else
      return P(l)(m)*phase;
  }
  
}

#define MONTJOIE_FILE_HELMHOLTZ_RADIAL_INLINE_CXX
#endif
