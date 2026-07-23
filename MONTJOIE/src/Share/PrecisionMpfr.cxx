#ifndef MONTJOIE_FILE_PRECISION_MPFR_CXX

namespace std
{
  
  //! reads a multiple precision complex from input stream
  istream& operator>>(istream& in, complex<mpfr::mpreal>& z)
  {
    string s;
    in >> s;
    if (s.size() > 0)
      {
        int pos_comma = -1, pos_left = -1, pos_right = -1;
        for (unsigned i = 0; i < s.size(); i++)
          {
            switch (s[i])
              {
              case ',' : pos_comma = i; break;
              case '(' : pos_left = i; break;
              case ')' : pos_right = i; break;
              }
          }
        
        if (pos_comma == -1)
          {
            // case where z is a real number
            z = complex<mpfr::mpreal>(mpfr::mpreal(s), 0);
          }
        else
          {
            // complex number
            mpfr::mpreal xr(s.substr(pos_left+1, pos_comma-pos_left-1));
            mpfr::mpreal xi(s.substr(pos_comma+1, pos_right-pos_comma-1));
            z = complex<mpfr::mpreal>(xr, xi);
          }
      }
    
    return in;
  }

}
 

#define MONTJOIE_FILE_PRECISION_MPFR_CXX
#endif
