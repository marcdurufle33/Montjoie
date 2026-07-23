#define MONTJOIE_WITH_ONE_DIM

#include "Hyperbolic/WaterWaves/MontjoieKdv.hxx"

using namespace Montjoie;

// file testing 1-D finite-difference schemes

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  {
    cout << "Testing periodic conditions" << endl;
    
    // testing periodic conditions
    FiniteDifferenceScheme1D<KdvEquation> scheme, scheme2;
    
    scheme.SetBoundaryCondition(scheme.PERIODIC, scheme.PERIODIC);
    
    Real_wp xmin = 0.0, xmax = 3.0;
    int N = 200;
    scheme.SetInterval(xmin, xmax, N);
    
    Vector<Real_wp> x, x2;
    scheme.GetSubdivisionSpace(x);
    
    int Nd = 500;
    scheme.InitGrid(xmin, xmax, Nd);
    VectReal_wp xgrid;
    Linspace(xmin, xmax, Nd, xgrid);
    
    Vector<Real_wp> Usol(N);
    Real_wp L = xmax-xmin;
    for (int i = 0; i < N; i++)
      Usol(i) = cos(2.0*pi_wp*x(i)/L);
    
    x.Write("x.dat");
    Usol.Write("Usol.dat");
    
    VectReal_wp Un;
    scheme.GetInterpolateUn(0, Un, Usol);
        
    Un.Write("Uinterp.dat");    

    // checking that the interpolation is of second order
    Real_wp err1, err2;
    err1 = 0;
    for (int i = 0; i < Nd; i++)
      {
        Real_wp u_exact = cos(2.0*pi_wp*xgrid(i)/L);
        err1 += square(u_exact - Un(i));
      }
    
    err1 = sqrt(err1/Nd);
    
    scheme2.SetBoundaryCondition(scheme.PERIODIC, scheme.PERIODIC);
    scheme2.SetInterval(xmin, xmax, 2*N);
    
    scheme2.GetSubdivisionSpace(x2);
    
    scheme2.InitGrid(xmin, xmax, Nd);
    
    Vector<Real_wp> Vsol(2*N);
    for (int i = 0; i < 2*N; i++)
      Vsol(i) = cos(2.0*pi_wp*x2(i)/L);
    
    scheme2.GetInterpolateUn(0, Un, Vsol);
    
    err2 = 0;
    for (int i = 0; i < Nd; i++)
      {
        Real_wp u_exact = cos(2.0*pi_wp*xgrid(i)/L);
        err2 += square(u_exact - Un(i));
      }
    
    err2 = sqrt(err2/Nd);
    cout << "Error on interpolation = " << err1 << endl;
    if (abs(err1/err2 - 4.0) > 0.1)
      {
        cout << "Not a second order interpolation" << endl;
        abort();
      }
    
    // checking the derivatives of u by using matrices D1, D2 and D3
    Real_wp dx = scheme.GetSpaceStep(), dx2 = scheme2.GetSpaceStep();
    Real_wp coef_dx = 1.0, coef2_dx = 1.0;
    Real_wp coef_dx2 = 1.0, coef2_dx2 = 1.0;
    Real_wp coef_dx3 = 1.0, coef2_dx3 = 1.0;
    
    TinyArrowMatrix<Real_wp, 1, 1> D1, D2, D1_r, D2_r;
    TinyArrowMatrix<Real_wp, 2, 2> D3, D3_r;
    D1.Reallocate(N, N);     D2.Reallocate(N, N);
    D1_r.Reallocate(2*N, 2*N);     D2_r.Reallocate(2*N, 2*N);
    D3.Reallocate(N, N); D3_r.Reallocate(2*N, 2*N);
    Vector<Real_wp> ones(N), ones2(2*N); ones.Fill(1.0); ones2.Fill(1.0);
    scheme.AddMatrixD1(coef_dx, ones, ones, D1);
    scheme.AddMatrixD2(coef_dx2, ones, ones, D2);
    scheme.AddMatrixD3(coef_dx3, ones, ones, D3);

    scheme2.AddMatrixD1(coef2_dx, ones2, ones2, D1_r);
    scheme2.AddMatrixD2(coef2_dx2, ones2, ones2, D2_r);
    scheme2.AddMatrixD3(coef2_dx3, ones2, ones2, D3_r);
    
    Vector<Real_wp> dUsol(N), dUsolb(N), dVsol(2*N), dVsolb(2*N);
    Vector<Real_wp> d2Usol(N), d2Usolb(N), d2Vsol(2*N), d2Vsolb(2*N);
    Vector<Real_wp> d3Usol(N), d3Usolb(N), d3Vsol(2*N), d3Vsolb(2*N);
    dUsol.Fill(0); dUsolb.Fill(0); dVsol.Fill(0); dVsolb.Fill(0);
    d2Usol.Fill(0); d2Usolb.Fill(0); d2Vsol.Fill(0); d2Vsolb.Fill(0);
    d3Usol.Fill(0); d3Usolb.Fill(0); d3Vsol.Fill(0); d3Vsolb.Fill(0);
    
    Mlt(D1, Usol, dUsol);     Mlt(D1_r, Vsol, dVsol);
    Mlt(D2, Usol, d2Usol);    Mlt(D2_r, Vsol, d2Vsol);
    Mlt(D3, Usol, d3Usol);    Mlt(D3_r, Vsol, d3Vsol);
    
    scheme.MltMatrixD1(coef_dx, Usol, dUsolb);     scheme2.MltMatrixD1(coef2_dx, Vsol, dVsolb);
    scheme.MltMatrixD2(coef_dx2, Usol, d2Usolb);    scheme2.MltMatrixD2(coef2_dx2, Vsol, d2Vsolb);
    scheme.MltMatrixD3(coef_dx3, Usol, d3Usolb);    scheme2.MltMatrixD3(coef2_dx3, Vsol, d3Vsolb);
    
    Real_wp err_du, err_d2u, err_d3u, err_dv, err_d2v, err_d3v;
    err_du = 0.0; err_d2u = 0.0; err_d3u = 0.0;
    for (int i = 0; i < N; i++)
      {
        Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
        Real_wp d2u_ex = -square(2.0*pi_wp/L)*cos(2.0*pi_wp*x(i)/L);
        Real_wp d3u_ex = pow(2.0*pi_wp/L, 3.0)*sin(2.0*pi_wp*x(i)/L);
        
        err_du += square(du_ex - dUsol(i)); 
        err_d2u += square(d2u_ex - d2Usol(i)); 
        err_d3u += square(d3u_ex - d3Usol(i)); 
      }
    
    err_du = sqrt(err_du*dx);
    err_d2u = sqrt(err_d2u*dx);
    err_d3u = sqrt(err_d3u*dx);
    
    err_dv = 0.0; err_d2v = 0.0; err_d3v = 0.0;
    for (int i = 0; i < 2*N; i++)
      {
        Real_wp dv_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
        Real_wp d2v_ex = -square(2.0*pi_wp/L)*cos(2.0*pi_wp*x2(i)/L);
        Real_wp d3v_ex = pow(2.0*pi_wp/L, 3.0)*sin(2.0*pi_wp*x2(i)/L);
        
        err_dv += square(dv_ex - dVsol(i)); 
        err_d2v += square(d2v_ex - d2Vsol(i)); 
        err_d3v += square(d3v_ex - d3Vsol(i)); 
        
        if (abs(dVsol(i) - dv_ex) > 1e-4)
          {
            DISP(i);
            DISP(dVsol(i)); DISP(dv_ex);
            abort();
          }
        
        if (abs(dVsol(i) - dVsolb(i)) > 1e-12)
          {
            DISP(i);
            DISP(dVsol(i)); DISP(dVsolb(i));
            abort();
          }

        if (abs(d2Vsol(i) - d2v_ex) > 1e-4)
          {
            DISP(i);
            DISP(d2Vsol(i)); DISP(d2v_ex);
            abort();
          }
        
        if (abs(d2Vsol(i) - d2Vsolb(i)) > 1e-10)
          {
            DISP(i);
            DISP(d2Vsol(i)); DISP(d2Vsolb(i));
            abort();
          }

        if (abs(d3Vsol(i) - d3v_ex) > 1e-2)
          {
            DISP(i);
            DISP(d3Vsol(i)); DISP(d3v_ex);
            abort();
          }
        
        if (abs(d3Vsol(i) - d3Vsolb(i)) > 1e-8)
          {
            DISP(d3Vsol(i)); DISP(d3Vsolb(i));
            DISP(d3Vsol(i)-d3Vsolb(i));
            abort();
          }
      }

    err_dv = sqrt(err_dv*dx2);
    err_d2v = sqrt(err_d2v*dx2);
    err_d3v = sqrt(err_d3v*dx2);
    
    cout << "Erreur sur du/dx = " << err_dv << endl;
    cout << "Erreur sur d^2 u/dx^2 = " << err_d2v << endl;
    cout << "Erreur sur d^3 u/dx^3 = " << err_d3v << endl;
    
    if (abs(err_du/err_dv - 4.0) > 0.1)
      {
        cout << "Not second order approximation of du/dx" << endl;
        DISP(err_du/err_dv);
        abort();
      }
    
    if (abs(err_d2u/err_d2v - 4.0) > 0.1)
      {
        cout << "Not second order approximation of d^2 u/dx^2" << endl;
        DISP(err_d2u/err_d2v);
        abort();
      }    
    
    if (abs(err_d3u/err_d3v - 4.0) > 0.1)
      {
        cout << "Not second order approximation of d^3 u/dx^3" << endl;
        DISP(err_d3u/err_d3v);
        abort();
      }
    
    // checking operator c du/dx + c_x/2 u
    Vector<Real_wp> c_half(N), c_half2(2*N), c_eval(N), c_eval2(2*N);
    for (int i = 0; i < 2*N; i++)
      {
        c_eval2(i) = 1.0 + 0.2*cos(2.0*pi_wp*x2(i)/L);
        c_half2(i) = 1.0 + 0.2*cos(2.0*pi_wp*(x2(i)+0.5*dx2)/L);
      }
    
    for (int i = 0; i < N; i++)
      {
        c_eval(i) = 1.0 + 0.2*cos(2.0*pi_wp*x(i)/L);
        c_half(i) = 1.0 + 0.2*cos(2.0*pi_wp*(x(i)+0.5*dx)/L);
      }
    
    D1.Zero(); D1_r.Zero();
    scheme.AddMatrixD1var(coef_dx, c_half, D1);
    scheme2.AddMatrixD1var(coef2_dx, c_half2, D1_r);
    
    Mlt(D1, Usol, dUsol);
    Mlt(D1_r, Vsol, dVsol);
    
    scheme.MltMatrixD1var(coef_dx, c_half, Usol, dUsolb);
    scheme2.MltMatrixD1var(coef2_dx, c_half2, Vsol, dVsolb);

    err_du = 0.0;
    for (int i = 0; i < N; i++)
      {
        Real_wp u_ex = cos(2.0*pi_wp*x(i)/L);
        Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
        Real_wp c = 1.0 + 0.2*cos(2.0*pi_wp*x(i)/L);
        Real_wp dc = -0.4*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
        Real_wp val_ex = c*du_ex + 0.5*dc*u_ex;
        
        err_du += square(val_ex - dUsol(i)); 
      }
    
    err_du = sqrt(err_du*dx);
    
    cout << "Erreur sur c du/dx + c_x/2 u = " << err_du << endl;
    
    err_dv = 0.0;
    for (int i = 0; i < 2*N; i++)
      {
        Real_wp u_ex = cos(2.0*pi_wp*x2(i)/L);
        Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
        Real_wp c = 1.0 + 0.2*cos(2.0*pi_wp*x2(i)/L);
        Real_wp dc = -0.4*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
        Real_wp val_ex = c*du_ex + 0.5*dc*u_ex;
        
        err_dv += square(val_ex - dVsol(i)); 
        
        if (abs(dVsol(i) - val_ex) > 1e-3)
          {
            DISP(i);
            DISP(dVsol(i)); DISP(val_ex);
            abort();
          }
        
        if (abs(dVsol(i) - dVsolb(i)) > 1e-12)
          {
            DISP(i);
            DISP(dVsol(i)); DISP(dVsolb(i));
            abort();
          }
      }

    err_dv = sqrt(err_dv*dx2);
    
    if (abs(err_du/err_dv - 4.0) > 0.1)
      {
        cout << "Not second order approximation of c u_x + c_x/2 u" << endl;
        DISP(err_du/err_dv);
        abort();
      }
    
    // checking operator d/dx( c du/dx)
    D2.Zero(); D2_r.Zero();
    scheme.AddMatrixD2var(coef_dx2, c_half, D2);
    scheme2.AddMatrixD2var(coef2_dx2, c_half2, D2_r);
    
    Mlt(D2, Usol, d2Usol);
    Mlt(D2_r, Vsol, d2Vsol);
    
    scheme.MltMatrixD2var(coef_dx2, c_half, Usol, d2Usolb);
    scheme2.MltMatrixD2var(coef2_dx2, c_half2, Vsol, d2Vsolb);

    err_d2u = 0.0;
    for (int i = 0; i < N; i++)
      {
        //Real_wp u_ex = cos(2.0*pi_wp*x(i)/L);
        Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
        Real_wp d2u_ex = -square(2.0*pi_wp/L)*cos(2.0*pi_wp*x(i)/L);
        Real_wp c = 1.0 + 0.2*cos(2.0*pi_wp*x(i)/L);
        Real_wp dc = -0.4*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
        Real_wp val_ex = c*d2u_ex + dc*du_ex;
        
        err_d2u += square(val_ex - d2Usol(i)); 
      }
    
    err_d2u = sqrt(err_d2u*dx);
    
    cout << "Erreur sur d/dx( c du/dx) = " << err_d2u << endl;
    
    err_d2v = 0.0;
    for (int i = 0; i < 2*N; i++)
      {
        //Real_wp u_ex = cos(2.0*pi_wp*x2(i)/L);
        Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
        Real_wp d2u_ex = -square(2.0*pi_wp/L)*cos(2.0*pi_wp*x2(i)/L);
        Real_wp c = 1.0 + 0.2*cos(2.0*pi_wp*x2(i)/L);
        Real_wp dc = -0.4*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
        Real_wp val_ex = c*d2u_ex + dc*du_ex;
        
        err_d2v += square(val_ex - d2Vsol(i)); 
        
        if (abs(d2Vsol(i) - val_ex) > 1e-3)
          {
            DISP(i);
            DISP(d2Vsol(i)); DISP(val_ex);
            abort();
          }
        
        if (abs(d2Vsol(i) - d2Vsolb(i)) > 1e-10)
          {
            DISP(i);
            DISP(d2Vsol(i)); DISP(d2Vsolb(i));
            abort();
          }
      }

    err_d2v = sqrt(err_d2v*dx2);
    
    if (abs(err_d2u/err_d2v - 4.0) > 0.1)
      {
        cout << "Not second order approximation of d/dx( c du/dx ) " << endl;
        DISP(err_d2u/err_d2v);
        abort();
      }
    
    // checking operator c u_xxx + 3/2 c_x u_xx + 3/4 c_xx u_x + 1/8 c_xxx u
    D3.Zero(); D3_r.Zero();
    scheme.AddMatrixD3var(coef_dx3, c_eval, D3);
    scheme2.AddMatrixD3var(coef2_dx3, c_eval2, D3_r);
    
    Mlt(D3, Usol, d3Usol);
    Mlt(D3_r, Vsol, d3Vsol);
    
    scheme.MltMatrixD3var(coef_dx3, c_eval, Usol, d3Usolb);
    scheme2.MltMatrixD3var(coef2_dx3, c_eval2, Vsol, d3Vsolb);

    err_d3u = 0.0;
    for (int i = 0; i < N; i++)
      {
        Real_wp u_ex = cos(2.0*pi_wp*x(i)/L);
        Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
        Real_wp d2u_ex = -square(2.0*pi_wp/L)*cos(2.0*pi_wp*x(i)/L);
        Real_wp d3u_ex = pow(2.0*pi_wp/L, 3.0)*sin(2.0*pi_wp*x(i)/L);
        Real_wp c = 1.0 + 0.2*cos(2.0*pi_wp*x(i)/L);
        Real_wp dc = -0.4*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
        Real_wp d2c = -0.2*square(2.0*pi_wp/L)*cos(2.0*pi_wp*x(i)/L);
        Real_wp d3c = 0.2*pow(2.0*pi_wp/L, 3.0)*sin(2.0*pi_wp*x(i)/L);
        Real_wp val_ex = c*d3u_ex + 1.5*dc*d2u_ex + 0.75*d2c*du_ex + 0.125*d3c*u_ex;
        
        err_d3u += square(val_ex - d3Usol(i)); 
      }
    
    err_d3u = sqrt(err_d3u*dx);
    
    cout << "Erreur sur c u_xxx + 3/2 c_x u_xx + 3/4 c_xx u_x + 1/8 c_xxx u = " << err_d3u << endl;
    
    err_d3v = 0.0;
    for (int i = 0; i < 2*N; i++)
      {
        Real_wp u_ex = cos(2.0*pi_wp*x2(i)/L);
        Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
        Real_wp d2u_ex = -square(2.0*pi_wp/L)*cos(2.0*pi_wp*x2(i)/L);
        Real_wp d3u_ex = pow(2.0*pi_wp/L, 3.0)*sin(2.0*pi_wp*x2(i)/L);
        Real_wp c = 1.0 + 0.2*cos(2.0*pi_wp*x2(i)/L);
        Real_wp dc = -0.4*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
        Real_wp d2c = -0.2*square(2.0*pi_wp/L)*cos(2.0*pi_wp*x2(i)/L);
        Real_wp d3c = 0.2*pow(2.0*pi_wp/L, 3.0)*sin(2.0*pi_wp*x2(i)/L);
        Real_wp val_ex = c*d3u_ex + 1.5*dc*d2u_ex + 0.75*d2c*du_ex + 0.125*d3c*u_ex;
        
        err_d3v += square(val_ex - d3Vsol(i)); 
        
        if (abs(d3Vsol(i) - val_ex) > 1.0)
          {
            DISP(i);
            DISP(d3Vsol(i)); DISP(val_ex);
            abort();
          }
        
        if (abs(d3Vsol(i) - d3Vsolb(i)) > 1e-8)
          {
            DISP(i);
            DISP(d3Vsol(i)); DISP(d3Vsolb(i));
            abort();
          }
      }
    
    err_d3v = sqrt(err_d3v*dx2);
    DISP(err_d3v);
    
    //if (abs(err_d3u/err_d3v - 4.0) > 0.1)
    //{
    //  cout << "Not second order approximation of Gamma_3 u" << endl;
    //  DISP(err_d3u/err_d3v);
    //  abort();
    //}
    
    // checking operator u^p u_x
    for (int p = 1; p <= 3; p++)
      {
        dUsol.Fill(0);
        scheme.AddOperatorUpUx(1.0, p, ones, ones, Usol, Usol, dUsol);
        
        dVsol.Fill(0);
        scheme2.AddOperatorUpUx(1.0, p, ones2, ones2, Vsol, Vsol, dVsol);
        
        err_du = 0;
        for (int i = 0; i < N; i++)
          {
            Real_wp u_ex = cos(2.0*pi_wp*x(i)/L);
            Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
            Real_wp val_ex = pow(u_ex, p)*du_ex;
            
            err_du += square(val_ex - dUsol(i));
          }
        
        err_du = sqrt(err_du*dx);
        cout << "Erreur sur u^"<<p<<" u_x = " << err_du << endl;
        
        err_dv = 0.0;
        for (int i = 0; i < 2*N; i++)
          {
            Real_wp u_ex = cos(2.0*pi_wp*x2(i)/L);
            Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
            Real_wp val_ex = pow(u_ex, p)*du_ex;
            
            err_dv += square(val_ex - dVsol(i));
            if (abs(val_ex - dVsol(i)) > 1e-3)
              {
                DISP(i); DISP(val_ex); DISP(dVsol(i));
                abort();
              }
          }
        
        err_dv = sqrt(err_dv*dx2);
        
        if (abs(err_du/err_dv - 4.0) > 0.1)
          {
            cout << "Not second order approximation of u^"<<p<<" u_x" << endl;
            DISP(err_du/err_dv);
            abort();
          }
        
        // checking conservation of energy
        VectReal_wp Us(N), Un_demi(N), Utmp(N), gU_plus(N), gU_minus(N);
        VectReal_wp coef_row(N), coef_col(N);
        coef_row.FillRand(); Mlt(1e-9, coef_row);
        coef_col.FillRand(); Mlt(1e-9, coef_col);
        Us.FillRand(); Mlt(1e-9, Us);        
        Un_demi.FillRand(); Mlt(1e-9, Un_demi);
        
        dUsol.Fill(0);
        scheme.AddOperatorUpUx(1.0, p, coef_col, coef_col, Us, Un_demi, dUsol);
        
        err_du = DotProd(dUsol, Us);
        DISP(err_du);
        if (abs(err_du) > 1e-12)
          {
            cout << "Energy not conserved" << endl;
            DISP(err_du);
            abort();
          }
        
        // checking differential matrix
        TinyArrowMatrix<Real_wp, 1, 1> DF;
        DF.Reallocate(N, N); 
        Real_wp h = 1e-6;
        scheme.AddDifferentialHalfUpUx(1.0, p, coef_row, coef_col, Un_demi, DF);
        for (int j = 0; j < N; j++)
          {
            Copy(Us, Utmp);
            Utmp(j) += h;
            
            gU_plus.Fill(0);
            scheme.AddOperatorUpUx(1.0, p, coef_row, coef_col, Utmp, Un_demi, gU_plus);
            
            Utmp(j) -= 2.0*h;
            gU_minus.Fill(0);
            scheme.AddOperatorUpUx(1.0, p, coef_row, coef_col, Utmp, Un_demi, gU_minus);
            
            for (int i = 0; i < N; i++)
              {
                Real_wp df_num = (gU_plus(i) - gU_minus(i))/(2.0*h);
                if (abs(df_num - DF(i, j)) > 1e-6)
                  {
                    DISP(i); DISP(j); DISP(DF(i, j)); DISP(df_num);
                    abort();
                  }
              }
          }
        
        DF.Zero();
        scheme.AddDifferentialUpUx(1.0, p, coef_row, coef_col, Us, DF);
        for (int j = 0; j < N; j++)
          {
            Copy(Us, Utmp);
            Utmp(j) += h;
            
            gU_plus.Fill(0);
            scheme.AddOperatorUpUx(1.0, p, coef_row, coef_col, Utmp, Utmp, gU_plus);
            
            Utmp(j) -= 2.0*h;
            gU_minus.Fill(0);
            scheme.AddOperatorUpUx(1.0, p, coef_row, coef_col, Utmp, Utmp, gU_minus);
            
            for (int i = 0; i < N; i++)
              {
                Real_wp df_num = (gU_plus(i) - gU_minus(i))/(2.0*h);
                if (abs(df_num - DF(i, j)) > 1e-6)
                  {
                    DISP(i); DISP(j); DISP(DF(i, j)); DISP(df_num);
                    abort();
                  }
              }
          }
      }
    
    // checking operator 2 u_xx u_x + u_xxx u
    {
      dUsol.Fill(0);
      scheme.AddOperatorUxxUx_UxxxU(1.0, ones, ones,
                                    Usol, Usol, dUsol);
        
      dVsol.Fill(0);
      scheme2.AddOperatorUxxUx_UxxxU(1.0, ones2, ones2,
                                     Vsol, Vsol, dVsol);
      
      err_du = 0;
      for (int i = 0; i < N; i++)
        {
          Real_wp u_ex = cos(2.0*pi_wp*x(i)/L);
          Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
          Real_wp d2u_ex = -square(2.0*pi_wp/L)*cos(2.0*pi_wp*x(i)/L);
          Real_wp d3u_ex = pow(2.0*pi_wp/L, 3.0)*sin(2.0*pi_wp*x(i)/L);
          Real_wp val_ex = 2.0*d2u_ex*du_ex + d3u_ex*u_ex;
          
          err_du += square(val_ex - dUsol(i));
        }
      
      err_du = sqrt(err_du*dx);
      cout << "Erreur sur 2 u_xx u_x + u_xxx u = " << err_du << endl;
      
      err_dv = 0.0;
      for (int i = 0; i < 2*N; i++)
        {
          Real_wp u_ex = cos(2.0*pi_wp*x2(i)/L);
          Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
          Real_wp d2u_ex = -square(2.0*pi_wp/L)*cos(2.0*pi_wp*x2(i)/L);
          Real_wp d3u_ex = pow(2.0*pi_wp/L, 3.0)*sin(2.0*pi_wp*x2(i)/L);
          Real_wp val_ex = 2.0*d2u_ex*du_ex + d3u_ex*u_ex;
          
          err_dv += square(val_ex - dVsol(i));
          if (abs(val_ex - dVsol(i)) > 1e-2)
            {
              DISP(i); DISP(val_ex); DISP(dVsol(i));
              abort();
            }
        }
        
      err_dv = sqrt(err_dv*dx2);
        
      if (abs(err_du/err_dv - 4.0) > 0.1)
        {
          cout << "Not second order approximation of 2 u_xx u_x + u_xxx u" << endl;
          DISP(err_du/err_dv);
          abort();
        }
        
      // checking conservation of energy
      VectReal_wp Us(N), Un_demi(N), Utmp(N), gU_plus(N), gU_minus(N);
      VectReal_wp coef_row(N), coef_col(N);
      coef_row.FillRand(); Mlt(1e-9, coef_row);
      coef_col.FillRand(); Mlt(1e-9, coef_col);
      Us.FillRand(); Mlt(1e-9, Us);        
      Un_demi.FillRand(); Mlt(1e-9, Un_demi);
      
      dUsol.Fill(0);
      scheme.AddOperatorUxxUx_UxxxU(1.0, coef_col, coef_col, Us, Un_demi, dUsol);
      
      err_du = DotProd(dUsol, Us);
      DISP(err_du);
      if (abs(err_du) > 1e-8)
        {
          cout << "Energy not conserved" << endl;
          DISP(err_du);
          abort();
        }

      // checking differential matrix
      TinyArrowMatrix<Real_wp, 2, 2> DF;
      DF.Reallocate(N, N); 
      Real_wp h = 1e-6;
      scheme.AddDifferentialHalfUxxUx_UxxxU(1.0, coef_row, coef_col, Un_demi, DF);
      for (int j = 0; j < N; j++)
        {
          Copy(Us, Utmp);
          Utmp(j) += h;
          
          gU_plus.Fill(0);
          scheme.AddOperatorUxxUx_UxxxU(1.0, coef_row, coef_col, Utmp, Un_demi, gU_plus);
          
          Utmp(j) -= 2.0*h;
          gU_minus.Fill(0);
          scheme.AddOperatorUxxUx_UxxxU(1.0, coef_row, coef_col, Utmp, Un_demi, gU_minus);
          
          for (int i = 0; i < N; i++)
            {
              Real_wp df_num = (gU_plus(i) - gU_minus(i))/(2.0*h);
              if (abs(df_num - DF(i, j)) > 1e-3)
                {
                  DISP(i); DISP(j); DISP(DF(i, j)); DISP(df_num);
                  abort();
                }
            }
        }
      
      DF.Zero();
      scheme.AddDifferentialUxxUx_UxxxU(1.0, coef_row, coef_col, Us, DF);
      for (int j = 0; j < N; j++)
        {
          Copy(Us, Utmp);
          Utmp(j) += h;
          
          gU_plus.Fill(0);
          scheme.AddOperatorUxxUx_UxxxU(1.0, coef_row, coef_col, Utmp, Utmp, gU_plus);
          
          Utmp(j) -= 2.0*h;
          gU_minus.Fill(0);
          scheme.AddOperatorUxxUx_UxxxU(1.0, coef_row, coef_col, Utmp, Utmp, gU_minus);
          
          for (int i = 0; i < N; i++)
            {
              Real_wp df_num = (gU_plus(i) - gU_minus(i))/(2.0*h);
              if (abs(df_num - DF(i, j)) > 1e-2)
                {
                  DISP(i); DISP(j); DISP(DF(i, j)); DISP(df_num);
                  abort();
                }
            }
        }

    }


    // checking operator u_xxx u
    {
      dUsol.Fill(0);
      scheme.AddOperatorUxxxU(1.0, ones, ones,
                              Usol, Usol, dUsol);
        
      dVsol.Fill(0);
      scheme2.AddOperatorUxxxU(1.0, ones2, ones2,
                               Vsol, Vsol, dVsol);
      
      err_du = 0;
      for (int i = 0; i < N; i++)
        {
          Real_wp u_ex = cos(2.0*pi_wp*x(i)/L);
          Real_wp d3u_ex = pow(2.0*pi_wp/L, 3.0)*sin(2.0*pi_wp*x(i)/L);
          Real_wp val_ex = d3u_ex*u_ex;
          
          err_du += square(val_ex - dUsol(i));
        }
      
      err_du = sqrt(err_du*dx);
      cout << "Erreur sur u_xxx u = " << err_du << endl;
      
      err_dv = 0.0;
      for (int i = 0; i < 2*N; i++)
        {
          Real_wp u_ex = cos(2.0*pi_wp*x2(i)/L);
          Real_wp d3u_ex = pow(2.0*pi_wp/L, 3.0)*sin(2.0*pi_wp*x2(i)/L);
          Real_wp val_ex = d3u_ex*u_ex;
          
          err_dv += square(val_ex - dVsol(i));
          if (abs(val_ex - dVsol(i)) > 1e-2)
            {
              DISP(i); DISP(val_ex); DISP(dVsol(i));
              abort();
            }
        }
        
      err_dv = sqrt(err_dv*dx2);
        
      if (abs(err_du/err_dv - 4.0) > 0.1)
        {
          cout << "Not second order approximation of u_xxx u" << endl;
          DISP(err_du/err_dv);
          abort();
        }
        
      VectReal_wp Us(N), Un_demi(N), Utmp(N), gU_plus(N), gU_minus(N);
      VectReal_wp coef_row(N), coef_col(N);
      coef_row.FillRand(); Mlt(1e-9, coef_row);
      coef_col.FillRand(); Mlt(1e-9, coef_col);
      Us.FillRand(); Mlt(1e-9, Us);        
      Un_demi.FillRand(); Mlt(1e-9, Un_demi);
      
      dUsol.Fill(0);
      scheme.AddOperatorUxxxU(1.0, coef_col, coef_col, Us, Un_demi, dUsol);

      // checking differential matrix
      TinyArrowMatrix<Real_wp, 2, 2> DF;
      DF.Reallocate(N, N); 
      Real_wp h = 1e-6;
      scheme.AddDifferentialHalfUxxxU(1.0, coef_row, coef_col, Un_demi, DF);
      for (int j = 0; j < N; j++)
        {
          Copy(Us, Utmp);
          Utmp(j) += h;
          
          gU_plus.Fill(0);
          scheme.AddOperatorUxxxU(1.0, coef_row, coef_col, Utmp, Un_demi, gU_plus);
          
          Utmp(j) -= 2.0*h;
          gU_minus.Fill(0);
          scheme.AddOperatorUxxxU(1.0, coef_row, coef_col, Utmp, Un_demi, gU_minus);
          
          for (int i = 0; i < N; i++)
            {
              Real_wp df_num = (gU_plus(i) - gU_minus(i))/(2.0*h);
              if (abs(df_num - DF(i, j)) > 1e-3)
                {
                  DISP(i); DISP(j); DISP(DF(i, j)); DISP(df_num);
                  abort();
                }
            }
        }
      
      DF.Zero();
      scheme.AddDifferentialUxxxU(1.0, coef_row, coef_col, Us, DF);
      for (int j = 0; j < N; j++)
        {
          Copy(Us, Utmp);
          Utmp(j) += h;
          
          gU_plus.Fill(0);
          scheme.AddOperatorUxxxU(1.0, coef_row, coef_col, Utmp, Utmp, gU_plus);
          
          Utmp(j) -= 2.0*h;
          gU_minus.Fill(0);
          scheme.AddOperatorUxxxU(1.0, coef_row, coef_col, Utmp, Utmp, gU_minus);
          
          for (int i = 0; i < N; i++)
            {
              Real_wp df_num = (gU_plus(i) - gU_minus(i))/(2.0*h);
              if (abs(df_num - DF(i, j)) > 1e-2)
                {
                  DISP(i); DISP(j); DISP(DF(i, j)); DISP(df_num);
                  abort();
                }
            }
        }

    }


    // checking operator u_xx (c u_x + c_x/2 u)
    {
      dUsol.Fill(0);
      scheme.AddOperatorUxxUx(1.0, ones, ones, c_half,
                              Usol, Usol, dUsol);
        
      dVsol.Fill(0);
      scheme2.AddOperatorUxxUx(1.0, ones2, ones2, c_half2,
                               Vsol, Vsol, dVsol);
      
      err_du = 0;
      for (int i = 0; i < N; i++)
        {
          Real_wp u_ex = cos(2.0*pi_wp*x(i)/L);
          Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
          Real_wp d2u_ex = -square(2.0*pi_wp/L)*cos(2.0*pi_wp*x(i)/L);
          Real_wp c = 1.0 + 0.2*cos(2.0*pi_wp*x(i)/L);
          Real_wp dc = -0.4*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
          Real_wp val_ex = d2u_ex*(c*du_ex + 0.5*dc*u_ex);
          
          err_du += square(val_ex - dUsol(i));
        }
      
      err_du = sqrt(err_du*dx);
      cout << "Erreur sur u_xx (c u_x + c_x/2 u) = " << err_du << endl;
      
      err_dv = 0.0;
      for (int i = 0; i < 2*N; i++)
        {
          Real_wp u_ex = cos(2.0*pi_wp*x2(i)/L);
          Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
          Real_wp d2u_ex = -square(2.0*pi_wp/L)*cos(2.0*pi_wp*x2(i)/L);
          Real_wp c = 1.0 + 0.2*cos(2.0*pi_wp*x2(i)/L);
          Real_wp dc = -0.4*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
          Real_wp val_ex = d2u_ex*(c*du_ex + 0.5*dc*u_ex);
          
          err_dv += square(val_ex - dVsol(i));
          if (abs(val_ex - dVsol(i)) > 1e-2)
            {
              DISP(i); DISP(val_ex); DISP(dVsol(i));
              abort();
            }
        }
        
      err_dv = sqrt(err_dv*dx2);
        
      if (abs(err_du/err_dv - 4.0) > 0.1)
        {
          cout << "Not second order approximation of u_xx u_x" << endl;
          DISP(err_du/err_dv);
          abort();
        }
        
      VectReal_wp Us(N), Un_demi(N), Utmp(N), gU_plus(N), gU_minus(N);
      VectReal_wp coef_row(N), coef_col(N);
      coef_row.FillRand(); Mlt(1e-9, coef_row);
      coef_col.FillRand(); Mlt(1e-9, coef_col);
      Us.FillRand(); Mlt(1e-9, Us);        
      Un_demi.FillRand(); Mlt(1e-9, Un_demi);
      
      dUsol.Fill(0);
      scheme.AddOperatorUxxUx(1.0, coef_col, coef_col, c_half, Us, Un_demi, dUsol);

      // checking differential matrix
      TinyArrowMatrix<Real_wp, 2, 2> DF;
      DF.Reallocate(N, N); 
      Real_wp h = 1e-6;
      scheme.AddDifferentialHalfUxxUx(1.0, coef_row, coef_col, c_half, Un_demi, DF);
      for (int j = 0; j < N; j++)
        {
          Copy(Us, Utmp);
          Utmp(j) += h;
          
          gU_plus.Fill(0);
          scheme.AddOperatorUxxUx(1.0, coef_row, coef_col, c_half, Utmp, Un_demi, gU_plus);
          
          Utmp(j) -= 2.0*h;
          gU_minus.Fill(0);
          scheme.AddOperatorUxxUx(1.0, coef_row, coef_col, c_half, Utmp, Un_demi, gU_minus);
          
          for (int i = 0; i < N; i++)
            {
              Real_wp df_num = (gU_plus(i) - gU_minus(i))/(2.0*h);
              if (abs(df_num - DF(i, j)) > 1e-3)
                {
                  DISP(i); DISP(j); DISP(DF(i, j)); DISP(df_num);
                  abort();
                }
            }
        }
      
      DF.Zero();
      scheme.AddDifferentialUxxUx(1.0, coef_row, coef_col, c_half, Us, DF);
      for (int j = 0; j < N; j++)
        {
          Copy(Us, Utmp);
          Utmp(j) += h;
          
          gU_plus.Fill(0);
          scheme.AddOperatorUxxUx(1.0, coef_row, coef_col, c_half, Utmp, Utmp, gU_plus);
          
          Utmp(j) -= 2.0*h;
          gU_minus.Fill(0);
          scheme.AddOperatorUxxUx(1.0, coef_row, coef_col, c_half, Utmp, Utmp, gU_minus);
          
          for (int i = 0; i < N; i++)
            {
              Real_wp df_num = (gU_plus(i) - gU_minus(i))/(2.0*h);
              if (abs(df_num - DF(i, j)) > 1e-2)
                {
                  DISP(i); DISP(j); DISP(DF(i, j)); DISP(df_num);
                  abort();
                }
            }
        }

    }


    // checking operator u_x (c u_x + c_x/2 u)
    {
      dUsol.Fill(0);
      scheme.AddOperatorUxUx(1.0, ones, ones, c_half,
                             Usol, Usol, dUsol);
      
      dVsol.Fill(0);
      scheme2.AddOperatorUxUx(1.0, ones2, ones2, c_half2,
                              Vsol, Vsol, dVsol);
      
      err_du = 0;
      for (int i = 0; i < N; i++)
        {
          Real_wp u_ex = cos(2.0*pi_wp*x(i)/L);
          Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
          Real_wp c = 1.0 + 0.2*cos(2.0*pi_wp*x(i)/L);
          Real_wp dc = -0.4*pi_wp/L*sin(2.0*pi_wp*x(i)/L);
          Real_wp val_ex = du_ex*(c*du_ex + 0.5*dc*u_ex);
          
          err_du += square(val_ex - dUsol(i));
        }
      
      err_du = sqrt(err_du*dx);
      cout << "Erreur sur u_x (c u_x + c_x/2 u) = " << err_du << endl;
      
      err_dv = 0.0;
      for (int i = 0; i < 2*N; i++)
        {
          Real_wp u_ex = cos(2.0*pi_wp*x2(i)/L);
          Real_wp du_ex = -2.0*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
          Real_wp c = 1.0 + 0.2*cos(2.0*pi_wp*x2(i)/L);
          Real_wp dc = -0.4*pi_wp/L*sin(2.0*pi_wp*x2(i)/L);
          Real_wp val_ex = du_ex*(c*du_ex + 0.5*dc*u_ex);
          
          err_dv += square(val_ex - dVsol(i));
          if (abs(val_ex - dVsol(i)) > 1e-2)
            {
              DISP(i); DISP(val_ex); DISP(dVsol(i));
              abort();
            }
        }
        
      err_dv = sqrt(err_dv*dx2);
        
      if (abs(err_du/err_dv - 4.0) > 0.1)
        {
          cout << "Not second order approximation of u_x u_x" << endl;
          DISP(err_du/err_dv);
          abort();
        }
        
      VectReal_wp Us(N), Un_demi(N), Utmp(N), gU_plus(N), gU_minus(N);
      VectReal_wp coef_row(N), coef_col(N);
      coef_row.FillRand(); Mlt(1e-9, coef_row);
      coef_col.FillRand(); Mlt(1e-9, coef_col);
      Us.FillRand(); Mlt(1e-9, Us);        
      Un_demi.FillRand(); Mlt(1e-9, Un_demi);
      
      dUsol.Fill(0);
      scheme.AddOperatorUxUx(1.0, coef_col, coef_col, c_half, Us, Un_demi, dUsol);

      // checking differential matrix
      TinyArrowMatrix<Real_wp, 1, 1> DF;
      DF.Reallocate(N, N); 
      Real_wp h = 1e-6;
      scheme.AddDifferentialHalfUxUx(1.0, coef_row, coef_col, c_half, Un_demi, DF);
      for (int j = 0; j < N; j++)
        {
          Copy(Us, Utmp);
          Utmp(j) += h;
          
          gU_plus.Fill(0);
          scheme.AddOperatorUxUx(1.0, coef_row, coef_col, c_half, Utmp, Un_demi, gU_plus);
          
          Utmp(j) -= 2.0*h;
          gU_minus.Fill(0);
          scheme.AddOperatorUxUx(1.0, coef_row, coef_col, c_half, Utmp, Un_demi, gU_minus);
          
          for (int i = 0; i < N; i++)
            {
              Real_wp df_num = (gU_plus(i) - gU_minus(i))/(2.0*h);
              if (abs(df_num - DF(i, j)) > 1e-3)
                {
                  DISP(i); DISP(j); DISP(DF(i, j)); DISP(df_num);
                  abort();
                }
            }
        }
      
      DF.Zero();
      scheme.AddDifferentialUxUx(1.0, coef_row, coef_col, c_half, Us, DF);
      for (int j = 0; j < N; j++)
        {
          Copy(Us, Utmp);
          Utmp(j) += h;
          
          gU_plus.Fill(0);
          scheme.AddOperatorUxUx(1.0, coef_row, coef_col, c_half, Utmp, Utmp, gU_plus);
          
          Utmp(j) -= 2.0*h;
          gU_minus.Fill(0);
          scheme.AddOperatorUxUx(1.0, coef_row, coef_col, c_half, Utmp, Utmp, gU_minus);
          
          for (int i = 0; i < N; i++)
            {
              Real_wp df_num = (gU_plus(i) - gU_minus(i))/(2.0*h);
              if (abs(df_num - DF(i, j)) > 1e-2)
                {
                  DISP(i); DISP(j); DISP(DF(i, j)); DISP(df_num);
                  abort();
                }
            }
        }

    }
    
  }
  
  {
    cout << "Testing Dirichlet conditions" << endl;
    
    // testing Dirichlet conditions
    FiniteDifferenceScheme1D<KdvEquation> scheme;
    
    scheme.SetBoundaryCondition(scheme.DIRICHLET, scheme.DIRICHLET);
    
    Real_wp xmin = 0.0, xmax = 3.0;
    int N = 200;
    scheme.SetInterval(xmin, xmax, N);
    
    Vector<Real_wp> x;
    scheme.GetSubdivisionSpace(x);
    
    int Nd = 250;
    scheme.InitGrid(xmin, xmax, Nd);
    
    Vector<Real_wp> Usol(N), xgrid;
    Linspace(xmin, xmax, Nd, xgrid);
    
    Real_wp L = xmax-xmin;
    for (int i = 0; i < N; i++)
      Usol(i) = sin(2.0*pi_wp*x(i)/L);
    
    x.Write("x_dir.dat");
    Usol.Write("Usol_dir.dat");
    
    VectReal_wp Un;
    scheme.GetInterpolateUn(0, Un, Usol);
    
    Un.Write("Uinterp_dir.dat");    
    // checking that the interpolation is of second order
    Real_wp err1, err2;
    err1 = 0;
    for (int i = 0; i < Nd; i++)
      {
        Real_wp u_exact = sin(2.0*pi_wp*xgrid(i)/L);
        err1 += square(u_exact - Un(i));
      }
    
    err1 = sqrt(err1/Nd);
    
    N *= 2;
    scheme.SetInterval(xmin, xmax, N);
    scheme.GetSubdivisionSpace(x);
    
    scheme.InitGrid(xmin, xmax, Nd);
    
    Usol.Reallocate(N);
    for (int i = 0; i < N; i++)
      Usol(i) = sin(2.0*pi_wp*x(i)/L);
    
    scheme.GetInterpolateUn(0, Un, Usol);
    
    err2 = 0;
    for (int i = 0; i < Nd; i++)
      {
        Real_wp u_exact = sin(2.0*pi_wp*xgrid(i)/L);
        err2 += square(u_exact - Un(i));
      }
    
    err2 = sqrt(err2/Nd);
    DISP(err2); DISP(err1);
    if (abs(err1/err2 - 4.0) > 0.1)
      {
        cout << "Not a second order interpolation" << endl;
        abort();
      }

    // checking the derivatives of u by using matrix D1
    //Real_wp dx = scheme.GetSpaceStep();
    Real_wp alpha = 1.0;
    
    TinyBandMatrix<Real_wp, 1> D1;
    D1.Reallocate(N, N);
    Vector<Real_wp> ones(N); ones.Fill(1.0);
    scheme.AddMatrixD1(alpha, ones, ones, D1);
    
    Vector<Real_wp> dUsol(N), dUsolb(N);
    dUsol.Fill(0); dUsolb.Fill(0);
    Mlt(D1, Usol, dUsol);
    
    scheme.MltMatrixD1(alpha, Usol, dUsolb);
    
    for (int i = 0; i < N; i++)
      {
        Real_wp du_ex = 2.0*pi_wp/L*cos(2.0*pi_wp*x(i)/L);
        if (abs(dUsol(i) - du_ex) > 1e-4)
          {
            DISP(i);
            DISP(dUsol(i)); DISP(du_ex);
            abort();
          }
        
        if (abs(dUsol(i) - dUsolb(i)) > 1e-12)
          {
            DISP(dUsol(i)); DISP(dUsolb(i));
            abort();
          }
      }
  }
}
