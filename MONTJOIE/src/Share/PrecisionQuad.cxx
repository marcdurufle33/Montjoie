#ifndef MONTJOIE_FILE_PRECISION_QUAD_CXX

namespace std
{
  
  ostream& operator<<(ostream& out, const __float128& x)
  {
    string format = "%";
    ios_base::fmtflags f = out.flags();
    if(f & ios_base::showpos)
      format += "+";
    if(f & ios_base::showpoint)
      format += "#";
    
    int digits = min(35, int(out.precision()));  
    format += "." + Seldon::to_str(digits);
    
    format += "Q";
    if(f & ios_base::scientific)
      format += "e";
    else if(f & ios_base::fixed)
      format += "f";
    else
      format += "g";
    
    char* tab;
    tab = reinterpret_cast<char*>(malloc(digits+10));
    quadmath_snprintf(tab, digits+10, format.data(), digits, x);
    out << tab;
    free(tab);
    return out;
  }


  istream& operator>>(istream& in, __float128& x)
  {
    string s; 
    in >> s;
    x = strtoflt128(s.data(), NULL);
    return in;
  }


  istream& operator>>(istream& in, complex<__float128>& z)
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
            z = complex<__float128>(Seldon::to_num<__float128>(s), 0);
          }
        else
          {
            // complex number
            __float128 xr = Seldon::to_num<__float128>(s.substr(pos_left+1, pos_comma-pos_left-1));
            __float128 xi = Seldon::to_num<__float128>(s.substr(pos_comma+1, pos_right-pos_comma-1));
            z = complex<__float128>(xr, xi);
          }
      }

    return in;
  }

}
 

#define MONTJOIE_FILE_PRECISION_QUAD_CXX
#endif
