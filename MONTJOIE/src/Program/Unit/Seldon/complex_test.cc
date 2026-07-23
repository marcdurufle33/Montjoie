#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  cout.precision(34);
  
  Real_wp threshold = 1e-15;

  Complex_wp a;
  a = to_num<Complex_wp>("2.4309583409580349580239403");
  DISP(a);
  
  a = to_num<Complex_wp>("(2.4309583409580349580239403,-0.345)");
  DISP(a);
  
  //cin >> a;
  //DISP(a);

  //Real_wp arg(0.4);
  //double arg_d(0.4);

  //DISP(yn(10, arg_d));
  //DISP(_yn(10, arg));
  
  complex<double> x, y, z;
  Complex_wp xt, yt, zt;
  
  {
    // testing operator =
    // operator = cannot be overloaded outside class
    // therefore z = x doesn't work for complex and multiple precision
    // use SetComplexReal instead
    x = 0;
    SetComplexReal(0, xt);
    if ((abs(x-xt) > threshold) || isnan(xt))
      {
        abort();
      }
    
    x = -3;
    SetComplexReal(-3, xt);
    if ((abs(x-xt) > threshold) || isnan(xt))
      {
        abort();
      }
    
    x = 21;
    SetComplexReal(21, xt);
    if ((abs(x-xt) > threshold) || isnan(xt))
      {
        abort();
      }
        
    x = 0.45;
    SetComplexReal(0.45, xt);
    if ((abs(x-xt) > threshold) || isnan(xt))
      {
        abort();
      }
    
    x = 2.3e-5;
    SetComplexReal(2.3e-5, xt);
    if ((abs(x-xt) > threshold) || isnan(xt))
      {
        abort();
      }
    
    x = to_num<double>("0.38");
    SetComplexReal(to_num<Real_wp>("0.38"), xt);
    if ((abs(x-xt) > threshold) || isnan(xt))
      {
        abort();
      }
    
    /*xt = 0; DISP(xt);
      xt = -3; DISP(xt);
      xt = 21; DISP(xt);
      xt = 0.45; DISP(xt);
      xt = 2.3e-5; DISP(xt);
      xt = to_num<Real_wp>("0.38"); DISP(xt);*/
  }
  
  {
    // testing operator ==
    x = 0;
    SetComplexReal(0, xt);
    bool test = (x==0), test_bis = (xt == 0);
    if (test != test_bis)
      {
        abort();
      }

    x = -3;
    SetComplexReal(-3, xt);
    test = (x==-3); test_bis = (xt == -3);
    if (test != test_bis)
      {
        abort();
      }

    x = 21;
    SetComplexReal(21, xt);
    test = (x==5); test_bis = (xt == 5);
    if (test != test_bis)
      {
        abort();
      }    

    x = 0.45;
    SetComplexReal(0.45, xt);
    test = (x == 0.45); test_bis = (xt == 0.45);
    if (test != test_bis)
      {
        abort();
      }    
    
    x = 2.3e-5;
    SetComplexReal(2.3e-5, xt);
    test = (x == 1.3e4); test_bis = (xt == 1.3e4);
    if (test != test_bis)
      {
        abort();
      }    
    
    x = to_num<double>("0.38");
    SetComplexReal(to_num<Real_wp>("0.38"), xt);
    test = (x == to_num<double>("0.38")); test_bis = (xt == to_num<Real_wp>("0.38")); 
    if (test != test_bis)
      {
        abort();
      }    

    x = to_num<double>("0.25");
    SetComplexReal(to_num<Real_wp>("0.25"), xt);
    test = (x == to_num<double>("0.34")); test_bis = (xt == to_num<Real_wp>("0.34")); 
    if (test != test_bis)
      {
        abort();
      }        
    
    x = complex<double>(1.6, 0.8);
    xt = Complex_wp(1.6, 0.8);
    if (xt == x)
      {
      }
    else
      {
        abort();
      }

    x = complex<double>(1.2, 0.7);
    xt = Complex_wp(1.6, 0.8);
    if (xt == x)
      {
        abort();
      }

    // testing operator !=
    x = 0;
    SetComplexReal(0, xt);
    test = (x != 0); test_bis = (xt != 0);
    if (test != test_bis)
      {
        abort();
      }

    x = -3;
    SetComplexReal(-3, xt);
    test = (x != -3); test_bis = (xt != -3);
    if (test != test_bis)
      {
        abort();
      }

    x = 21;
    SetComplexReal(21, xt);
    test = (x != 5); test_bis = (xt != 5);
    if (test != test_bis)
      {
        abort();
      }    

    x = 0.45;
    SetComplexReal(0.45, xt);
    test = (x != 0.45); test_bis = (xt != 0.45);
    if (test != test_bis)
      {
        abort();
      }    
    
    x = 2.3e-5;
    SetComplexReal(2.3e-5, xt);
    test = (x != 1.3e4); test_bis = (xt != 1.3e4);
    if (test != test_bis)
      {
        abort();
      }    
    
    x = to_num<double>("0.38");
    SetComplexReal(to_num<Real_wp>("0.38"), xt);
    test = (x != to_num<double>("0.38")); test_bis = (xt != to_num<Real_wp>("0.38")); 
    if (test != test_bis)
      {
        abort();
      }    

    x = to_num<double>("0.25");
    SetComplexReal(to_num<Real_wp>("0.25"), xt);
    test = (x != to_num<double>("0.34")); test_bis = (xt != to_num<Real_wp>("0.34")); 
    if (test != test_bis)
      {
        abort();
      }        
    
    x = complex<double>(1.6, 0.8);
    xt = Complex_wp(1.6, 0.8);
    if (xt != x)
      {
        abort();
      }

    x = complex<double>(1.2, 0.7);
    xt = Complex_wp(1.6, 0.8);
    if (xt != x)
      {
      }
    else
      {
        abort();
      }
    
    xt = Complex_wp(5, 6); yt = Complex_wp(5, -2); zt = Complex_wp(-1, 6);
    if (xt == yt)
      {
        abort();
      }

    if (xt == zt)
      {
        abort();
      }
    
    yt = Complex_wp(5, 6);
    if (xt == yt)
      {
      }
    else
      {
        abort();
      }

    xt = Complex_wp(5, 6); yt = Complex_wp(5, -2); zt = Complex_wp(-1, 6);
    if (xt != yt)
      {
      }
    else
      {
        abort();
      }

    if (xt != zt)
      {
      }
    else
      {
        abort();
      }
    
    yt = Complex_wp(5, 6);
    if (xt != yt)
      {
        abort();
      }    
  }

  {
    // testing operator +
    x = complex<double>(5, 1);
    y = x + 3;

    xt = complex<double>(5, 1);
    yt = xt + 3;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x + 4;

    xt = complex<double>(0, -3);
    yt = xt + 4;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = 3 + x;

    xt = complex<double>(5, 1);
    yt = 3 + xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = 4 + x;

    xt = complex<double>(0, -3);
    yt = 4 + xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }


    x = complex<double>(5, 1);
    y = x + 2.2;

    xt = complex<double>(5, 1);
    yt = xt + 2.2;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x + 0.4;

    xt = complex<double>(0, -3);
    yt = xt + 0.4;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = 2.9 + x;

    xt = complex<double>(5, 1);
    yt = 2.9 + xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = 3.2 + x;

    xt = complex<double>(0, -3);
    yt = 3.2 + xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = x + complex<double>(2, 3);

    xt = complex<double>(5, 1);
    yt = xt + complex<double>(2, 3);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x + complex<double>(0, 2.4);

    xt = complex<double>(0, -3);
    yt = xt + complex<double>(0, 2.4);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = complex<double>(3,0.8) + x;

    xt = complex<double>(5, 1);
    yt = complex<double>(3,0.8) + xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = complex<double>(0, 9.1) + x;

    xt = complex<double>(0, -3);
    yt = complex<double>(0, 9.1) + xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = x + double(3.8);

    xt = complex<double>(5, 1);
    yt = xt + Real_wp(3.8);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x + double(0.8);

    xt = complex<double>(0, -3);
    yt = xt + Real_wp(0.8);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = double(2.3) + x;

    xt = complex<double>(5, 1);
    yt = Real_wp(2.3) + xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = double(1.3) + x;

    xt = complex<double>(0, -3);
    yt = Real_wp(1.3) + xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(2, 4);
    y = x + complex<double>(5, 1);
    
    xt = Complex_wp(2, 4);
    yt = xt + Complex_wp(5, 1);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(3.2, 0.4);
    y = x + complex<double>(0, 1.8);
    
    xt = Complex_wp(3.2, 0.4);
    yt = xt + Complex_wp(0, 1.8);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 0.4);
    y = x + complex<double>(2.3, 0);
    
    xt = Complex_wp(0, 0.4);
    yt = xt + Complex_wp(2.3, 0);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(2, 4);
    y = complex<double>(5, 1) + x;
    
    xt = Complex_wp(2, 4);
    yt = Complex_wp(5, 1) + xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(3.2, 0.4);
    y = complex<double>(0, 1.8) + x;
    
    xt = Complex_wp(3.2, 0.4);
    yt = Complex_wp(0, 1.8) + xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 0.4);
    y = complex<double>(2.3, 0) + x;
    
    xt = Complex_wp(0, 0.4);
    yt = Complex_wp(2.3, 0) + xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

  }
  
  {
    // testing Operator -
    z = -complex<double>(3.2, 5.1);
    zt = -Complex_wp(3.2, 5.1);
    
    if (abs(zt-z) > threshold || isnan(abs(zt-z)))
      {
        abort();
      }

    z = -complex<double>(2.2, 0);
    zt = -Complex_wp(2.2, 0);
    
    if (abs(zt-z) > threshold || isnan(abs(zt-z)))
      {
        abort();
      }

    z = -complex<double>(0, 0.2);
    zt = -Complex_wp(0, 0.2);
    
    if (abs(zt-z) > threshold || isnan(abs(zt-z)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = x - 3;

    xt = complex<double>(5, 1);
    yt = xt - 3;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x - 4;

    xt = complex<double>(0, -3);
    yt = xt - 4;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = 3 - x;

    xt = complex<double>(5, 1);
    yt = 3 - xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = 4 - x;

    xt = complex<double>(0, -3);
    yt = 4 - xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }


    x = complex<double>(5, 1);
    y = x - 2.2;

    xt = complex<double>(5, 1);
    yt = xt - 2.2;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x - 0.4;

    xt = complex<double>(0, -3);
    yt = xt - 0.4;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = 2.9 - x;

    xt = complex<double>(5, 1);
    yt = 2.9 - xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = 3.2 - x;

    xt = complex<double>(0, -3);
    yt = 3.2 - xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = x - complex<double>(2, 3);

    xt = complex<double>(5, 1);
    yt = xt - complex<double>(2, 3);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x - complex<double>(0, 2.4);

    xt = complex<double>(0, -3);
    yt = xt - complex<double>(0, 2.4);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = complex<double>(3,0.8) - x;

    xt = complex<double>(5, 1);
    yt = complex<double>(3,0.8) - xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = complex<double>(0, 9.1) - x;

    xt = complex<double>(0, -3);
    yt = complex<double>(0, 9.1) - xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = x - double(3.8);

    xt = complex<double>(5, 1);
    yt = xt - Real_wp(3.8);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x - double(0.8);

    xt = complex<double>(0, -3);
    yt = xt - Real_wp(0.8);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = double(2.3) - x;

    xt = complex<double>(5, 1);
    yt = Real_wp(2.3) - xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = double(1.3) - x;

    xt = complex<double>(0, -3);
    yt = Real_wp(1.3) - xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(2, 4);
    y = x - complex<double>(5, 1);
    
    xt = Complex_wp(2, 4);
    yt = xt - Complex_wp(5, 1);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(3.2, 0.4);
    y = x - complex<double>(0, 1.8);
    
    xt = Complex_wp(3.2, 0.4);
    yt = xt - Complex_wp(0, 1.8);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 0.4);
    y = x - complex<double>(2.3, 0);
    
    xt = Complex_wp(0, 0.4);
    yt = xt - Complex_wp(2.3, 0);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(2, 4);
    y = complex<double>(5, 1) - x;
    
    xt = Complex_wp(2, 4);
    yt = Complex_wp(5, 1) - xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(3.2, 0.4);
    y = complex<double>(0, 1.8) - x;
    
    xt = Complex_wp(3.2, 0.4);
    yt = Complex_wp(0, 1.8) - xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 0.4);
    y = complex<double>(2.3, 0) - x;
    
    xt = Complex_wp(0, 0.4);
    yt = Complex_wp(2.3, 0) - xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }    
  }

  {
    // testing operator *
    x = complex<double>(5, 1);
    y = x * 3;

    xt = complex<double>(5, 1);
    yt = xt * 3;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x * 4;

    xt = complex<double>(0, -3);
    yt = xt * 4;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = 3 * x;

    xt = complex<double>(5, 1);
    yt = 3 * xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = 4 * x;

    xt = complex<double>(0, -3);
    yt = 4 * xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }


    x = complex<double>(5, 1);
    y = x * 2.2;

    xt = complex<double>(5, 1);
    yt = xt * 2.2;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x * 0.4;

    xt = complex<double>(0, -3);
    yt = xt * 0.4;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = 2.9 * x;

    xt = complex<double>(5, 1);
    yt = 2.9 * xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = 3.2 * x;

    xt = complex<double>(0, -3);
    yt = 3.2 * xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = x * complex<double>(2, 3);

    xt = complex<double>(5, 1);
    yt = xt * complex<double>(2, 3);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x * complex<double>(0, 2.4);

    xt = complex<double>(0, -3);
    yt = xt * complex<double>(0, 2.4);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = complex<double>(3,0.8) * x;

    xt = complex<double>(5, 1);
    yt = complex<double>(3,0.8) * xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = complex<double>(0, 9.1) * x;

    xt = complex<double>(0, -3);
    yt = complex<double>(0, 9.1) * xt;
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = x * double(3.8);

    xt = complex<double>(5, 1);
    yt = xt * Real_wp(3.8);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x * double(0.8);

    xt = complex<double>(0, -3);
    yt = xt * Real_wp(0.8);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = double(2.3) * x;

    xt = complex<double>(5, 1);
    yt = Real_wp(2.3) * xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = double(1.3) * x;

    xt = complex<double>(0, -3);
    yt = Real_wp(1.3) * xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(2, 4);
    y = x * complex<double>(5, 1);
    
    xt = Complex_wp(2, 4);
    yt = xt * Complex_wp(5, 1);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(3.2, 0.4);
    y = x * complex<double>(0, 1.8);
    
    xt = Complex_wp(3.2, 0.4);
    yt = xt * Complex_wp(0, 1.8);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 0.4);
    y = x * complex<double>(2.3, 0);
    
    xt = Complex_wp(0, 0.4);
    yt = xt * Complex_wp(2.3, 0);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(2, 4);
    y = complex<double>(5, 1) * x;
    
    xt = Complex_wp(2, 4);
    yt = Complex_wp(5, 1) * xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(3.2, 0.4);
    y = complex<double>(0, 1.8) * x;
    
    xt = Complex_wp(3.2, 0.4);
    yt = Complex_wp(0, 1.8) * xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 0.4);
    y = complex<double>(2.3, 0) * x;
    
    xt = Complex_wp(0, 0.4);
    yt = Complex_wp(2.3, 0) * xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
  }

  {
    // testing operator /
    x = complex<double>(5, 1);
    y = x / 3;

    xt = complex<double>(5, 1);
    yt = xt / 3;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x / 4;

    xt = complex<double>(0, -3);
    yt = xt / 4;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = 3 / x;

    xt = complex<double>(5, 1);
    yt = 3 / xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = 4 / x;

    xt = complex<double>(0, -3);
    yt = 4 / xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }


    x = complex<double>(5, 1);
    y = x / 2.2;

    xt = complex<double>(5, 1);
    yt = xt / 2.2;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x / 0.4;

    xt = complex<double>(0, -3);
    yt = xt / 0.4;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = 2.9 / x;

    xt = complex<double>(5, 1);
    yt = 2.9 / xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = 3.2 / x;

    xt = complex<double>(0, -3);
    yt = 3.2 / xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = x / complex<double>(2, 3);

    xt = complex<double>(5, 1);
    yt = xt / complex<double>(2, 3);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x / complex<double>(0, 2.4);

    xt = complex<double>(0, -3);
    yt = xt / complex<double>(0, 2.4);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = complex<double>(3,0.8) / x;

    xt = complex<double>(5, 1);
    yt = complex<double>(3,0.8) / xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = complex<double>(0, 9.1) / x;

    xt = complex<double>(0, -3);
    yt = complex<double>(0, 9.1) / xt;
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = x / double(3.8);

    xt = complex<double>(5, 1);
    yt = xt / Real_wp(3.8);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = x / double(0.8);

    xt = complex<double>(0, -3);
    yt = xt / Real_wp(0.8);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(5, 1);
    y = double(2.3) / x;

    xt = complex<double>(5, 1);
    yt = Real_wp(2.3) / xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -3);
    y = double(1.3) / x;

    xt = complex<double>(0, -3);
    yt = Real_wp(1.3) / xt;
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(2, 4);
    y = x / complex<double>(5, 1);
    
    xt = Complex_wp(2, 4);
    yt = xt / Complex_wp(5, 1);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(3.2, 0.4);
    y = x / complex<double>(0, 1.8);
    
    xt = Complex_wp(3.2, 0.4);
    yt = xt / Complex_wp(0, 1.8);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 0.4);
    y = x / complex<double>(2.3, 0);
    
    xt = Complex_wp(0, 0.4);
    yt = xt / Complex_wp(2.3, 0);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(2, 4);
    y = complex<double>(5, 1) / x;
    
    xt = Complex_wp(2, 4);
    yt = Complex_wp(5, 1) / xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(3.2, 0.4);
    y = complex<double>(0, 1.8) / x;
    
    xt = Complex_wp(3.2, 0.4);
    yt = Complex_wp(0, 1.8) / xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 0.4);
    y = complex<double>(2.3, 0) / x;
    
    xt = Complex_wp(0, 0.4);
    yt = Complex_wp(2.3, 0) / xt;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

  }

  {
    // Testing operator +=
    y = complex<double>(2, 3);
    y += 1;
    
    yt = Complex_wp(2, 3);
    yt += 1;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y += 2;
    
    yt = Complex_wp(0, 4);
    yt += 2;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y += 5;
    
    yt = Complex_wp(-5, 0);
    yt += 5;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }


    y = complex<double>(2, 3);
    y += 1.2;
    
    yt = Complex_wp(2, 3);
    yt += 1.2;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y += 2.9;
    
    yt = Complex_wp(0, 4);
    yt += 2.9;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y += 5.01;
    
    yt = Complex_wp(-5, 0);
    yt += 5.01;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y += complex<double>(1.2, -3);
    
    yt = Complex_wp(2, 3);
    yt += complex<double>(1.2, -3);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y += complex<double>(2.9, 0.8);
    
    yt = Complex_wp(0, 4);
    yt += complex<double>(2.9, 0.8);;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y += complex<double>(0, 5.01);
    
    yt = Complex_wp(-5, 0);
    yt += complex<double>(0, 5.01);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y += 1.2;
    
    yt = Complex_wp(2, 3);
    yt += Real_wp(1.2);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y += double(2.9);
    
    yt = Complex_wp(0, 4);
    yt += Real_wp(2.9);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y += double(5.01);
    
    yt = Complex_wp(-5, 0);
    yt += Real_wp(5.01);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y += complex<double>(1.2, -3);
    
    yt = Complex_wp(2, 3);
    yt += Complex_wp(1.2, -3);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y += complex<double>(2.9, 0.8);
    
    yt = Complex_wp(0, 4);
    yt += Complex_wp(2.9, 0.8);;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y += complex<double>(0, 5.01);
    
    yt = Complex_wp(-5, 0);
    yt += Complex_wp(0, 5.01);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
  }

  {
    // Testing operator -=
    y = complex<double>(2, 3);
    y -= 1;
    
    yt = Complex_wp(2, 3);
    yt -= 1;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y -= 2;
    
    yt = Complex_wp(0, 4);
    yt -= 2;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y -= 5;
    
    yt = Complex_wp(-5, 0);
    yt -= 5;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }


    y = complex<double>(2, 3);
    y -= 1.2;
    
    yt = Complex_wp(2, 3);
    yt -= 1.2;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y -= 2.9;
    
    yt = Complex_wp(0, 4);
    yt -= 2.9;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y -= 5.01;
    
    yt = Complex_wp(-5, 0);
    yt -= 5.01;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y -= complex<double>(1.2, -3);
    
    yt = Complex_wp(2, 3);
    yt -= complex<double>(1.2, -3);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y -= complex<double>(2.9, 0.8);
    
    yt = Complex_wp(0, 4);
    yt -= complex<double>(2.9, 0.8);;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y -= complex<double>(0, 5.01);
    
    yt = Complex_wp(-5, 0);
    yt -= complex<double>(0, 5.01);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y -= 1.2;
    
    yt = Complex_wp(2, 3);
    yt -= Real_wp(1.2);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y -= 2.9;
    
    yt = Complex_wp(0, 4);
    yt -= Real_wp(2.9);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y -= 5.01;
    
    yt = Complex_wp(-5, 0);
    yt -= Real_wp(5.01);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y -= complex<double>(1.2, -3);
    
    yt = Complex_wp(2, 3);
    yt -= Complex_wp(1.2, -3);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y -= complex<double>(2.9, 0.8);
    
    yt = Complex_wp(0, 4);
    yt -= Complex_wp(2.9, 0.8);;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y -= complex<double>(0, 5.01);
    
    yt = Complex_wp(-5, 0);
    yt -= Complex_wp(0, 5.01);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
  }


  {
    // Testing operator *=
    y = complex<double>(2, 3);
    y *= 3;
    
    yt = Complex_wp(2, 3);
    yt *= 3;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y *= 2;
    
    yt = Complex_wp(0, 4);
    yt *= 2;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y *= 5;
    
    yt = Complex_wp(-5, 0);
    yt *= 5;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }


    y = complex<double>(2, 3);
    y *= 1.2;
    
    yt = Complex_wp(2, 3);
    yt *= 1.2;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y *= 2.9;
    
    yt = Complex_wp(0, 4);
    yt *= 2.9;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y *= 5.01;
    
    yt = Complex_wp(-5, 0);
    yt *= 5.01;

    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y *= complex<double>(1.2, -3);
    
    yt = Complex_wp(2, 3);
    yt *= complex<double>(1.2, -3);

    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y *= complex<double>(2.9, 0.8);
    
    yt = Complex_wp(0, 4);
    yt *= complex<double>(2.9, 0.8);;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y *= complex<double>(0, 5.01);
    
    yt = Complex_wp(-5, 0);
    yt *= complex<double>(0, 5.01);

    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y *= 1.2;
    
    yt = Complex_wp(2, 3);
    yt *= Real_wp(1.2);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y *= 2.9;
    
    yt = Complex_wp(0, 4);
    yt *= Real_wp(2.9);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y *= 5.01;
    
    yt = Complex_wp(-5, 0);
    yt *= Real_wp(5.01);

    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y *= complex<double>(1.2, -3);
    
    yt = Complex_wp(2, 3);
    yt *= Complex_wp(1.2, -3);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y *= complex<double>(2.9, 0.8);
    
    yt = Complex_wp(0, 4);
    yt *= Complex_wp(2.9, 0.8);;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y *= complex<double>(0, 5.01);
    
    yt = Complex_wp(-5, 0);
    yt *= Complex_wp(0, 5.01);

    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
  }


  {
    // Testing operator /=
    y = complex<double>(2, 3);
    y /= 7;
    
    yt = Complex_wp(2, 3);
    yt /= 7;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y /= 11;
    
    yt = Complex_wp(0, 4);
    yt /= 11;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y /= 3;
    
    yt = Complex_wp(-5, 0);
    yt /= 3;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }


    y = complex<double>(2, 3);
    y /= 1.2;
    
    yt = Complex_wp(2, 3);
    yt /= 1.2;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y /= 2.9;
    
    yt = Complex_wp(0, 4);
    yt /= 2.9;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y /= 5.01;
    
    yt = Complex_wp(-5, 0);
    yt /= 5.01;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y /= complex<double>(1.2, -3);
    
    yt = Complex_wp(2, 3);
    yt /= complex<double>(1.2, -3);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y /= complex<double>(2.9, 0.8);
    
    yt = Complex_wp(0, 4);
    yt /= complex<double>(2.9, 0.8);;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y /= complex<double>(0, 5.01);
    
    yt = Complex_wp(-5, 0);
    yt /= complex<double>(0, 5.01);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y /= 1.2;
    
    yt = Complex_wp(2, 3);
    yt /= Real_wp(1.2);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y /= 2.9;
    
    yt = Complex_wp(0, 4);
    yt /= Real_wp(2.9);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y /= 5.01;
    
    yt = Complex_wp(-5, 0);
    yt /= Real_wp(5.01);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(2, 3);
    y /= complex<double>(1.2, -3);
    
    yt = Complex_wp(2, 3);
    yt /= Complex_wp(1.2, -3);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    y = complex<double>(0, 4);
    y /= complex<double>(2.9, 0.8);
    
    yt = Complex_wp(0, 4);
    yt /= Complex_wp(2.9, 0.8);;

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    y = complex<double>(-5, 0);
    y /= complex<double>(0, 5.01);
    
    yt = Complex_wp(-5, 0);
    yt /= Complex_wp(0, 5.01);

    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
  }
  
  {
    // Testing conj
    x = complex<double>(0.36, 2.13);
    y = conj(x);

    xt = Complex_wp(0.36, 2.13);
    yt = conj(xt);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.54, 0);
    y = conj(x);

    xt = Complex_wp(-0.54, 0);
    yt = conj(xt);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -4.2);
    y = conj(x);

    xt = Complex_wp(0, -4.2);
    yt = conj(xt);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }    
  }
  
  double val; 
  Real_wp val_t;
  
  {
    // Testing abs
    x = complex<double>(0.36, 2.13);
    val = abs(x);

    xt = Complex_wp(0.36, 2.13);
    val_t = abs(xt);
    
    if (abs(val_t-val) > threshold || isnan(abs(val_t-val)))
      {
        abort();
      }

    x = complex<double>(-0.54, 0);
    val = abs(x);

    xt = Complex_wp(-0.54, 0);
    val_t = abs(xt);
    
    if (abs(val_t-val) > threshold || isnan(abs(val_t-val)))
      {
        abort();
      }

    x = complex<double>(0, -4.2);
    val = abs(x);

    xt = Complex_wp(0, -4.2);
    val_t = abs(xt);
    
    if (abs(val_t-val) > threshold || isnan(abs(val_t-val)))
      {
        abort();
      }
  }

  {
    // Testing sqrt
    x = complex<double>(0.36, 2.13);
    y = sqrt(x);

    xt = Complex_wp(0.36, 2.13);
    yt = sqrt(xt);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.54, 0);
    y = sqrt(x);

    xt = Complex_wp(-0.54, 0);
    yt = sqrt(xt);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(43.2, 0);
    y = sqrt(x);

    xt = Complex_wp(43.2, 0);
    yt = sqrt(xt);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -4.2);
    y = sqrt(x);

    xt = Complex_wp(0, -4.2);
    yt = sqrt(xt);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 1.8);
    y = sqrt(x);

    xt = Complex_wp(0, 1.8);
    yt = sqrt(xt);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0.46, -3.63);
    y = sqrt(x);

    xt = Complex_wp(0.46, -3.63);
    yt = sqrt(xt);
    
    if (abs(yt-y) > threshold || isnan(abs(yt-y)))
      {
        abort();
      }
  }


  {
    // Testing exp
    x = complex<double>(0.36, 2.13);
    y = exp(x);
    
    xt = Complex_wp(0.36, 2.13);
    yt = exp(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(0, 10.46);
    y = exp(x);
    
    xt = Complex_wp(0, 10.46);
    yt = exp(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.83, 0);
    y = exp(x);
    
    xt = Complex_wp(-0.83, 0);
    yt = exp(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
  }
  
  {
    // Testing log
    x = complex<double>(0.36, 2.13);
    y = log(x);
    
    xt = Complex_wp(0.36, 2.13);
    yt = log(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(0, 10.46);
    y = log(x);
    
    xt = Complex_wp(0, 10.46);
    yt = log(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.83, 0);
    y = log(x);
    
    xt = Complex_wp(-0.83, 0);
    yt = log(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }    
  }

  {
    // Testing log10
    x = complex<double>(0.36, 2.13);
    y = log10(x);
    
    xt = Complex_wp(0.36, 2.13);
    yt = log10(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(0, 10.46);
    y = log10(x);
    
    xt = Complex_wp(0, 10.46);
    yt = log10(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.83, 0);
    y = log10(x);
    
    xt = Complex_wp(-0.83, 0);
    yt = log10(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }    
  }
  
  {
    // Testing pow
    x = complex<double>(0.36, 2.13);
    y = pow(x, 3);

    xt = Complex_wp(0.36, 2.13);
    yt = pow(xt, 3);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.54, 0);
    y = pow(x, 4);

    xt = Complex_wp(-0.54, 0);
    yt = pow(xt, 4);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(43.2, 0);
    y = pow(x, 2);

    xt = Complex_wp(43.2, 0);
    yt = pow(xt, 2);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -4.2);
    y = pow(x, 3);

    xt = Complex_wp(0, -4.2);
    yt = pow(xt, 3);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 1.8);
    y = pow(x, -4);

    xt = Complex_wp(0, 1.8);
    yt = pow(xt, -4);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0.46, -3.63);
    y = pow(x, -3);

    xt = Complex_wp(0.46, -3.63);
    yt = pow(xt, -3);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    // with a double
    x = complex<double>(0.36, 2.13);
    y = pow(x, 2.8);

    xt = Complex_wp(0.36, 2.13);
    yt = pow(xt, 2.8);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.54, 0);
    y = pow(x, 3.9);

    xt = Complex_wp(-0.54, 0);
    yt = pow(xt, 3.9);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(43.2, 0);
    y = pow(x, 1.7);

    xt = Complex_wp(43.2, 0);
    yt = pow(xt, 1.7);
    
    if (abs(yt-y) > 200*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -4.2);
    y = pow(x, 3.4);

    xt = Complex_wp(0, -4.2);
    yt = pow(xt, 3.4);
    
    if (abs(yt-y) > 100*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 1.8);
    y = pow(x, -4.2);

    xt = Complex_wp(0, 1.8);
    yt = pow(xt, -4.2);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0.46, -3.63);
    y = pow(x, -3.5);

    xt = Complex_wp(0.46, -3.63);
    yt = pow(xt, -3.5);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    // with a real
    x = complex<double>(0.36, 2.13);
    y = pow(x, 2.8);

    xt = Complex_wp(0.36, 2.13);
    yt = pow(xt, Real_wp(2.8));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.54, 0);
    y = pow(x, 3.9);

    xt = Complex_wp(-0.54, 0);
    yt = pow(xt, Real_wp(3.9));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(43.2, 0);
    y = pow(x, 1.7);

    xt = Complex_wp(43.2, 0);
    yt = pow(xt, Real_wp(1.7));
    
    if (abs(yt-y) > 200*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -4.2);
    y = pow(x, 3.4);

    xt = Complex_wp(0, -4.2);
    yt = pow(xt, Real_wp(3.4));
    
    if (abs(yt-y) > 100*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 1.8);
    y = pow(x, -4.2);

    xt = Complex_wp(0, 1.8);
    yt = pow(xt, Real_wp(-4.2));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0.46, -3.63);
    y = pow(x, -3.5);

    xt = Complex_wp(0.46, -3.63);
    yt = pow(xt, Real_wp(-3.5));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    // with a double complex
    x = complex<double>(0.36, 2.13);
    y = pow(x, complex<double>(2.8, 0.7));

    xt = Complex_wp(0.36, 2.13);
    yt = pow(xt, complex<double>(2.8, 0.7));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.54, 0);
    y = pow(x, complex<double>(3.9, -0.6));

    xt = Complex_wp(-0.54, 0);
    yt = pow(xt, complex<double>(3.9, -0.6));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(43.2, 0);
    y = pow(x, complex<double>(2.1, -2.2));

    xt = Complex_wp(43.2, 0);
    yt = pow(xt, complex<double>(2.1, -2.2));
    
    if (abs(yt-y) > 1000*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -4.2);
    y = pow(x, complex<double>(-3, 0));

    xt = Complex_wp(0, -4.2);
    yt = pow(xt, complex<double>(-3, 0));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 1.8);
    y = pow(x, complex<double>(-4.2, 0.3));

    xt = Complex_wp(0, 1.8);
    yt = pow(xt, complex<double>(-4.2, 0.3));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0.46, -3.63);
    y = pow(x, complex<double>(0, 1.1));

    xt = Complex_wp(0.46, -3.63);
    yt = pow(xt, complex<double>(0, 1.1));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    // with a complex
    x = complex<double>(0.36, 2.13);
    y = pow(x, complex<double>(2.8, 0.7));

    xt = Complex_wp(0.36, 2.13);
    yt = pow(xt, Complex_wp(2.8, 0.7));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.54, 0);
    y = pow(x, complex<double>(3.9, -0.6));

    xt = Complex_wp(-0.54, 0);
    yt = pow(xt, Complex_wp(3.9, -0.6));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(43.2, 0);
    y = pow(x, complex<double>(2.1, -2.2));

    xt = Complex_wp(43.2, 0);
    yt = pow(xt, Complex_wp(2.1, -2.2));
    
    if (abs(yt-y) > 1000*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, -4.2);
    y = pow(x, complex<double>(-3, 0));

    xt = Complex_wp(0, -4.2);
    yt = pow(xt, Complex_wp(-3, 0));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0, 1.8);
    y = pow(x, complex<double>(-4.2, 0.3));

    xt = Complex_wp(0, 1.8);
    yt = pow(xt, Complex_wp(-4.2, 0.3));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(0.46, -3.63);
    y = pow(x, complex<double>(0, 1.1));

    xt = Complex_wp(0.46, -3.63);
    yt = pow(xt, Complex_wp(0, 1.1));
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

  }

  {
    // Testing cos
    x = complex<double>(0.36, 2.13);
    y = cos(x);
    
    xt = Complex_wp(0.36, 2.13);
    yt = cos(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(0, 1.46);
    y = cos(x);
    
    xt = Complex_wp(0, 1.46);
    yt = cos(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.83, 0);
    y = cos(x);
    
    xt = Complex_wp(-0.83, 0);
    yt = cos(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }    
  }


  {
    // Testing sin
    x = complex<double>(0.36, 2.13);
    y = sin(x);
    
    xt = Complex_wp(0.36, 2.13);
    yt = sin(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(0, 1.46);
    y = sin(x);
    
    xt = Complex_wp(0, 1.46);
    yt = sin(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.83, 0);
    y = sin(x);
    
    xt = Complex_wp(-0.83, 0);
    yt = sin(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }    
  }  

  {
    // Testing tan
    x = complex<double>(0.36, 2.13);
    y = tan(x);
    
    xt = Complex_wp(0.36, 2.13);
    yt = tan(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(0, 1.46);
    y = tan(x);
    
    xt = Complex_wp(0, 1.46);
    yt = tan(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.83, 0);
    y = tan(x);
    
    xt = Complex_wp(-0.83, 0);
    yt = tan(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }    
  }


  {
    // Testing cosh
    x = complex<double>(0.36, 2.13);
    y = cosh(x);
    
    xt = Complex_wp(0.36, 2.13);
    yt = cosh(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(0, 1.46);
    y = cosh(x);
    
    xt = Complex_wp(0, 1.46);
    yt = cosh(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.83, 0);
    y = cosh(x);
    
    xt = Complex_wp(-0.83, 0);
    yt = cosh(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }    
  }


  {
    // Testing sinh
    x = complex<double>(0.36, 2.13);
    y = sinh(x);
    
    xt = Complex_wp(0.36, 2.13);
    yt = sinh(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(0, 1.46);
    y = sinh(x);
    
    xt = Complex_wp(0, 1.46);
    yt = sinh(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.83, 0);
    y = sinh(x);
    
    xt = Complex_wp(-0.83, 0);
    yt = sinh(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }    
  }  

  {
    // Testing tanh
    x = complex<double>(0.36, 2.13);
    y = tanh(x);
    
    xt = Complex_wp(0.36, 2.13);
    yt = tanh(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }
    
    x = complex<double>(0, 1.46);
    y = tanh(x);
    
    xt = Complex_wp(0, 1.46);
    yt = tanh(xt);
    
    if (abs(yt-y) > 100*threshold || isnan(abs(yt-y)))
      {
        abort();
      }

    x = complex<double>(-0.83, 0);
    y = tanh(x);
    
    xt = Complex_wp(-0.83, 0);
    yt = tanh(xt);
    
    if (abs(yt-y) > 10*threshold || isnan(abs(yt-y)))
      {
        abort();
      }    
  }


  cout << "All tests passed successfully" << endl;

  return FinalizeMontjoie();
}
