#ifndef MONTJOIE_FILE_HELMHOLTZ_POLAR_INLINE_CXX

namespace Montjoie
{

  template<class Complexe>
  inline int VarHelmholtz_Polar<Complexe>::GetLmax() const
  {
    return Lmax;
  }
  

  template<class Complexe>
  inline int VarHelmholtz_Polar<Complexe>::GetIndexOutputFiles(int type)
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
    
}

#define MONTJOIE_FILE_HELMHOLTZ_POLAR_INLINE_CXX
#endif
