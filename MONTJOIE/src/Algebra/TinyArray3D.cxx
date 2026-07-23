#ifndef SELDON_FILE_TINY_ARRAY3D_CXX

#include "TinyArray3D.hxx"

namespace Seldon
{  
  //! displays 3-D array
  template<class T, int m, int n, int p> inline
  ostream & operator << (ostream & out , const TinyArray3D<T, m, n, p> & A)
  {
    // we use normal loops (not optimal)
    for (int i = 0; i < m; i++)
      {
        for (int j = 0;  j < n; j++)
          {
            for (int k = 0; k < p; k++)
              out << A(i, j, k) << ' ';
            
            out << '\n';
          }
        
        out << endl << endl;
      }
    
    return out;
  }
      
} // end namespace


#define SELDON_FILE_TINY_ARRAY3D_CXX
#endif

