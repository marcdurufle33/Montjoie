#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

// #define COMPLEX_BESSEL

void WriteVal(const string& name, const Real_wp& x)
{
  ofstream file_out(name.data(), ios::app);
  file_out.precision(16);
  file_out << x << endl;
  file_out.close();
}

void WriteVal(const string& name, const Complex_wp& x)
{
  ofstream file_out(name.data(), ios::app);
  file_out.precision(16);
  file_out << real(x) << '\n' << imag(x) << '\n';
  file_out.close();
}

void WriteVal(const string& name, const VectReal_wp& x)
{
  ofstream file_out(name.data(), ios::app);
  file_out.precision(16);
  for (int i = 0; i < x.GetM(); i++)
    file_out << x(i) << '\n';
  
  file_out.close();
}

void WriteVal(const string& name, const VectComplex_wp& x)
{
  ofstream file_out(name.data(), ios::app);
  file_out.precision(16);
  for (int i = 0; i < x.GetM(); i++)
    file_out << real(x(i)) << '\n' << imag(x(i)) << '\n';
  
  file_out.close();
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  Real_wp z = GetWhittakerM(Real_wp(0.2), Real_wp(0.3), Real_wp(1.1));
  
  z = GetWhittakerW(Real_wp(0.2), Real_wp(0.3), Real_wp(1.1));
  
  if ((abs(GetLambertW0(Real_wp(0.7)) - 0.447470259269655) > 1e-12)
      || isnan(GetLambertW0(Real_wp(0.7))))
    {
      cout << "GetLambertW0 incorrect" << endl;
      abort();
    }
  
  if ((abs(GetLambertW0(Real_wp(2.3)) - 0.918223536799827) > 1e-12)
      || isnan(GetLambertW0(Real_wp(2.3))))
    {
      cout << "GetLambertW0 incorrect" << endl;
      abort();
    }

  if ((abs(GetLambertW0(Real_wp(3.9)) - 1.18838299651336) > 1e-12)
      || isnan(GetLambertW0(Real_wp(3.9))))
    {
      cout << "GetLambertW0 incorrect" << endl;
      abort();
    }

  if ((abs(GetLambertW0(Real_wp(15.2)) - 2.01879475195702) > 1e-12)
      || isnan(GetLambertW0(Real_wp(15.2))))
    {
      cout << "GetLambertW0 incorrect" << endl;
      abort();
    }
  
  remove("output.dat");
  string name("output.dat");
  
  cout.precision(14);
  cout.setf(ios::scientific);
 
  {
#ifdef COMPLEX_BESSEL
    DISP(GetKn(4, Real_wp(0.34)));
    DISP(GetKn(-4, Real_wp(0.34)));
#endif

    DISP(GetJn(3, 0.87)); WriteVal(name, GetJn(3, 0.87));
    DISP(GetJn(-3, 0.87)); WriteVal(name, GetJn(-3, 0.87));
    
    DISP(GetYn(2, 0.87)); WriteVal(name, GetYn(2, 0.87));
    DISP(GetYn(-2, 0.87)); WriteVal(name, GetYn(-2, 0.87));
    
    DISP(GetJn(0, 0.0)); WriteVal(name, GetJn(0, 0.0));
    DISP(GetJn(4, 0.0)); WriteVal(name, GetJn(4, 0.0));
  }

  {
    Real_wp jn, djn;
    GetDeriveJn(3, 0.87, jn, djn);
    DISP(jn); DISP(djn); WriteVal(name, jn); WriteVal(name, djn);
  }
  
  cout << endl << endl;
  
  cout << "Testing ComputeBesselFunctions..." << endl;
  
  {
    VectComplex_wp Jn, Yn;
#ifdef COMPLEX_BESSEL
    ComputeBesselFunctions(4, 10, Complex_wp(2.4, 0.2), Jn);
    DISP(Jn); Jn.Clear();
    ComputeBesselFunctions(4, 10, Complex_wp(2.4, 0.2), Jn, Yn);
    DISP(Jn); DISP(Yn);
#endif

    ComputeBesselFunctions(-4, 14, Complex_wp(2.2, 0.0), Jn);
    DISP(Jn); Jn.Clear(); WriteVal(name, Jn);
    ComputeBesselFunctions(-4, 14, Complex_wp(2.2, 0.0), Jn, Yn);
    DISP(Jn); DISP(Yn); WriteVal(name, Jn); WriteVal(name, Yn);
    
    VectReal_wp JnR;
    ComputeBesselFunctions(-2, 12, Real_wp(0.8), JnR);
    DISP(JnR); WriteVal(name, JnR);
    
    VectReal_wp dJnR; JnR.Clear();
    ComputeDeriveBesselFunctions(-2, 12, Real_wp(1.2), JnR, dJnR);
    DISP(JnR); DISP(dJnR); WriteVal(name, JnR); WriteVal(name, dJnR);
  }

  cout << endl << endl;
  cout << "Testing Hankel functions... " << endl;
  
  {
    VectComplex_wp Jn, Hn, dJn, dHn;
#ifdef COMPLEX_BESSEL
    ComputeBesselAndHankelFunctions(Real_wp(0), 11, Complex_wp(1.8, 0.1),
                                    Jn, Hn);
    DISP(Jn); DISP(Hn);

    ComputeBesselAndHankelFunctions(Real_wp(0.5), 11, Complex_wp(1.8, 0.1),
                                    Jn, Hn);
    DISP(Jn); DISP(Hn);
    
    Jn.Clear(); Hn.Clear();
    ComputeDeriveBesselAndHankel(Real_wp(0), 11, Complex_wp(1.8, 0.1),
                                 Jn, Hn, dJn, dHn); 
    DISP(Jn); DISP(Hn); DISP(dJn); DISP(dHn);
    
    ComputeDeriveBesselAndHankel(Real_wp(0.5), 11, Complex_wp(1.8, 0.1),
                                 Jn, Hn, dJn, dHn); 
    DISP(Jn); DISP(Hn); DISP(dJn); DISP(dHn);
#endif

    ComputeBesselAndHankelFunctions(Real_wp(0), 11, Complex_wp(1.4, 0.0),
                                    Jn, Hn);
    DISP(Jn); DISP(Hn); WriteVal(name, Jn); WriteVal(name, Hn);
    
#ifdef COMPLEX_BESSEL
    ComputeBesselAndHankelFunctions(Real_wp(0.5), 11, Complex_wp(1.4, 0.0),
                                    Jn, Hn);
    DISP(Jn); DISP(Hn);
#endif
    
    Jn.Clear(); Hn.Clear();
    ComputeDeriveBesselAndHankel(Real_wp(0), 11, Complex_wp(1.4, 0.0),
                                 Jn, Hn, dJn, dHn); 
    DISP(Jn); DISP(Hn); DISP(dJn); DISP(dHn);
    WriteVal(name, Jn); WriteVal(name, Hn);
    WriteVal(name, dJn); WriteVal(name, dHn);
    
#ifdef COMPLEX_BESSEL    
    ComputeDeriveBesselAndHankel(Real_wp(0.5), 11, Complex_wp(1.4, 0.0),
                                 Jn, Hn, dJn, dHn); 
    DISP(Jn); DISP(Hn); DISP(dJn); DISP(dHn);
#endif
  }
  
  cout << endl << endl;
  cout << "Testing spherical Bessel functions... " << endl;
  
  {
    VectReal_wp Jn, dJn; int n = 20;
    
    ComputeSphericalBessel(n, Real_wp(0.45), Jn);
    DISP(Jn); WriteVal(name, Jn);

    ComputeSphericalBessel(n, Real_wp(0.0), Jn);
    DISP(Jn); WriteVal(name, Jn);

    Jn.Clear();
    ComputeDeriveSphericalBessel(n, Real_wp(0.45), Jn, dJn);
    DISP(Jn); DISP(dJn); WriteVal(name, Jn); WriteVal(name, dJn);

    ComputeDeriveSphericalBessel(n, Real_wp(0.0), Jn, dJn);
    DISP(Jn); DISP(dJn); WriteVal(name, Jn); WriteVal(name, dJn);
  }
  
  cout << endl << endl;
  cout << "Testing spherical Hankel functions... " << endl;
  
  {
    VectComplex_wp Hn, dHn; int n = 18;
    
    ComputeSphericalHankel(n, Real_wp(0.58), Hn);
    DISP(Hn); WriteVal(name, Hn);

    ComputeSphericalHankel(n, Real_wp(0.0), Hn);
    DISP(Hn); WriteVal(name, Hn);

    ComputeDeriveSphericalHankel(n, Real_wp(0.58), Hn, dHn);
    DISP(Hn); DISP(dHn); WriteVal(name, Hn); WriteVal(name, dHn);

    ComputeDeriveSphericalHankel(n, Real_wp(0.0), Hn, dHn);
    DISP(Hn); DISP(dHn); WriteVal(name, Hn); WriteVal(name, dHn);
  }

  cout << endl << endl;
  cout << "Testing Ricatti functions ... " << endl;
  
  {
    VectComplex_wp Jn, Hn, dJn, dHn, Psi_n, Xi_n, dPsi, dXi;
#ifdef COMPLEX_BESSEL
    ComputeDeriveSphericalBesselHankel(Real_wp(0.5), 14, Complex_wp(1.26, 0.13),
                                       Jn, Hn, dJn, dHn);
    DISP(Jn); DISP(Hn); DISP(dJn); DISP(dHn);
    
    Jn.Clear(); Hn.Clear(); dJn.Clear(); dHn.Clear();
    ComputeDeriveRiccatiBessel(Real_wp(0.5), 16, Complex_wp(1.26, 0.13),
                               Jn, Hn, dJn, dHn, Psi_n, Xi_n, dPsi, dXi);
    
    DISP(Jn); DISP(Hn); DISP(dJn); DISP(dHn); DISP(Psi_n); DISP(Xi_n);
    DISP(dPsi); DISP(dXi);
#endif

    Jn.Clear(); Hn.Clear(); dJn.Clear(); dHn.Clear();
    ComputeDeriveSphericalBesselHankel(Real_wp(0.5), 8, Complex_wp(1.06, 0.0),
                                       Jn, Hn, dJn, dHn);
    DISP(Jn); DISP(Hn); DISP(dJn); DISP(dHn);
    WriteVal(name, Jn); WriteVal(name, dJn);
    WriteVal(name, Hn); WriteVal(name, dHn);
    
    Jn.Clear(); Hn.Clear(); dJn.Clear(); dHn.Clear();
    ComputeDeriveRiccatiBessel(Real_wp(0.5), 8, Complex_wp(1.06, 0.0),
                               Jn, Hn, dJn, dHn, Psi_n, Xi_n, dPsi, dXi);
    
    DISP(Jn); DISP(Hn); DISP(dJn); DISP(dHn); DISP(Psi_n); DISP(Xi_n);
    DISP(dPsi); DISP(dXi);
    WriteVal(name, Jn); WriteVal(name, dJn);
    WriteVal(name, Hn); WriteVal(name, dHn);
    WriteVal(name, Psi_n); WriteVal(name, Xi_n);
    WriteVal(name, dPsi); WriteVal(name, dXi);
   
  }

  {
    VectReal_wp val, val_ref;
    val.ReadText(name);
    val_ref.ReadText("src/Program/Unit/Share/output_bessel_ref.dat");
    DISP(val_ref.GetM()); DISP(val.GetM());
    Real_wp err = 0;
    for (int i = 0; i < val.GetM(); i++)
      if ((val(i) != 0) && (val_ref(i) != 0))
        {
          Real_wp err_i = abs(val(i) - val_ref(i)) / abs(val_ref(i));
          err = max(err, err_i);
        }
    
    DISP(err);
  }

  return FinalizeMontjoie();
}
