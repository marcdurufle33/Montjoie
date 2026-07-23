#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Harmonic/MontjoieHarmonic.hxx"

using namespace Montjoie;

Real_wp threshold = 1e-12;

template<class T>
void CheckRectangleCutOff(RectangleCutOff<Dimension2, T>& cut)
{
  // testing with no cut-off
  R2 xmin(-2, -1), xmax(3, 4), coef(0, 0);  
  TinyVector<R2, 2> enveloppe; enveloppe(0) = xmin; enveloppe(1) = xmax;
  cut.Init(xmin, xmax, coef);
  
  int N = 200;
  VectReal_wp x(N+1);
  for (int i = 0; i <= N; i++)
    x(i) = -5.0 + 10.0*Real_wp(i)/N;
  
  Matrix<Real_wp> nu(N+1, N+1);
  for (int i = 0; i <= N; i++)
    for (int j = 0; j <= N; j++)
      {
        R2 point(x(i), x(j));
        if (PointInsideBoundingBox(point, enveloppe))
          {
            nu(i, j) = 1.0;
            cut.ApplyCutOff(point, nu(i, j));
          }
        else
          nu(i, j) = 0.0;
      }
  
  nu.Write("no_cutoff.dat");
  
  coef.Init(0.01, 0.01);
  cut.Init(xmin, xmax, coef);
  for (int i = 0; i <= N; i++)
    for (int j = 0; j <= N; j++)
      {
        R2 point(x(i), x(j));
        if (PointInsideBoundingBox(point, enveloppe))
          {
            nu(i, j) = 1.0;
            cut.ApplyCutOff(point, nu(i, j));
          }
        else
          nu(i, j) = 0.0;
      }
  
  nu.Write("gaussian_cutoff.dat");
  
  // testing derivatives
  R2 grad, deriv, point_p;
  Real_wp h = 1e-6, val_m, val_p;
  for (int i = 0; i <= N; i++)
    for (int j = 0; j <= N; j++)
      {
        R2 point(x(i), x(j));
        if (PointInsideBoundingBox(point, enveloppe))
          {
            nu(i, j) = sin(pi_wp*x(i))*sin(2.0*pi_wp*x(j));
            grad.Init(pi_wp*cos(pi_wp*x(i))*sin(2.0*pi_wp*x(j)),
                      2.0*pi_wp*sin(pi_wp*x(i))*cos(2.0*pi_wp*x(j)));            
            
            cut.ApplyCutOff(point, nu(i, j), grad);
            
            point_p = point; point_p(0) += h;
            val_p = sin(pi_wp*(x(i)+h))*sin(2.0*pi_wp*x(j));
            cut.ApplyCutOff(point_p, val_p);
            
            point_p = point; point_p(0) -= h;
            val_m = sin(pi_wp*(x(i)-h))*sin(2.0*pi_wp*x(j));
            cut.ApplyCutOff(point_p, val_m);
            
            deriv(0) = (val_p - val_m) / (2.0*h);
            
            point_p = point; point_p(1) += h;
            val_p = sin(pi_wp*(x(i)))*sin(2.0*pi_wp*(x(j)+h));
            cut.ApplyCutOff(point_p, val_p);

            point_p = point; point_p(1) -= h;
            val_m = sin(pi_wp*(x(i)))*sin(2.0*pi_wp*(x(j)-h));
            cut.ApplyCutOff(point_p, val_m);
            deriv(1) = (val_p - val_m) / (2.0*h);
            
            if (grad.Distance(deriv) > 1e-6)
              {
                cout << "Derivatives not correctly computed" << endl;
                DISP(deriv); DISP(grad); abort();
              }
          }
        else
          nu(i, j) = 0;
      }

  nu.Write("sinus_cutoff.dat");
}


template<class T>
void CheckRectangleCutOff(RectangleCutOff<Dimension3, T>& cut)
{
  // to do
}

template<class T>
void CheckSinusIndex(PhysicalSinusoidalMedia<Dimension2, T>& index)
{
  R2 xmin(-2, -1), xmax(3, 4), tronc(0.01, 0.01);
  R2 kwave(pi_wp, 2.0*pi_wp), point;
  index.Init(xmin, xmax, tronc, 2.2, 0.5, kwave);
  
  if (index.IsComplex())
    {
      cout << "IsComplex incorrect" << endl;
      abort();
    }
  
  int N = 200;
  VectReal_wp x(N+1);
  for (int i = 0; i <= N; i++)
    x(i) = -5.0 + 10.0*Real_wp(i)/N;
  
  Matrix<Real_wp> nu;
  nu.Read("sinus_cutoff.dat");
  T val;
  for (int i = 0; i <= N; i ++)
    for (int j = 0; j <= N; j++)
      {
        point.Init(x(i), x(j));
        index.ComputeIndexAtPoint(point, val);
        if (abs(val - 2.2 - 0.5*nu(i,j)) > 1e-10)
          {
            cout << "ComputeIndexAtPoint incorrect" << endl;
            DISP(val); DISP(2.2+0.5*nu(i,j));
            abort();
          }
      }
  
  Complex_wp coef(1.5, 0.2);
  index.Mlt(coef);
  if (!index.IsComplex())
    {
      cout << "IsComplex incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i <= N; i ++)
    for (int j = 0; j <= N; j++)
      {
        point.Init(x(i), x(j));
        index.ComputeIndexAtPoint(point, val);
        if (abs(val - coef*(2.2 + 0.5*nu(i,j))) > 1e-10)
          {
            cout << "Mlt incorrect" << endl;
            abort();
          }
      }
  
  index.GetInverse();
  for (int i = 0; i <= N; i ++)
    for (int j = 0; j <= N; j++)
      {
        point.Init(x(i), x(j));
        index.ComputeIndexAtPoint(point, val);
        if (abs(val - 1.0/(coef*(2.2 + 0.5*nu(i,j)))) > 1e-10)
          {
            cout << "GetInverse incorrect" << endl;
            abort();
          }
      }
  
  index.GetInverse();
  TinyVector<T, 2> grad, deriv;
  R2 point_p;
  Real_wp h = 1e-6; T val_m, val_p;
  for (int i = 0; i <= N; i ++)
    for (int j = 0; j <= N; j++)
      {
        point.Init(x(i), x(j));
        index.ComputeIndexAtPoint(point, val, grad);
        
        point_p = point; point_p(0) += h;
        index.ComputeIndexAtPoint(point_p, val_p);

        point_p = point; point_p(0) -= h;
        index.ComputeIndexAtPoint(point_p, val_m);
        deriv(0) = (val_p - val_m) / (2.0*h);
        
        point_p = point; point_p(1) += h;
        index.ComputeIndexAtPoint(point_p, val_p);

        point_p = point; point_p(1) -= h;
        index.ComputeIndexAtPoint(point_p, val_m);
        deriv(1) = (val_p - val_m) / (2.0*h);
        if (abs(grad.Distance(deriv)) > 1e-6)
          {
            cout << "Derivatives not correctly computed" << endl;
            DISP(deriv); DISP(grad); abort();
          }
      }
  
}

template<class T>
void CheckSinusIndex(PhysicalSinusoidalMedia<Dimension3, T>& index)
{
}

template<class T>
void CheckRegularIndex(PhysicalRegularMedia<Dimension2, T>& index)
{  
  R2 xmin(-2.0, -1.0), xmax(3.0, 4.0), tronc(0.01, 0.01);
  int N = 150;
  Matrix<double> nu(N+1, N+1);
  for (int i = 0; i <= N; i++)
    for (int j = 0; j <= N; j++)
      {
        double x = xmin(0) + (xmax(0)-xmin(0))*double(i)/N;
        double y = xmin(1) + (xmax(1)-xmin(1))*double(j)/N;
        nu(i, j) = sin(pi_wp*x)*sin(2.0*pi_wp*y);
      }
  
  nu.Write("index.dat");
  
  index.Init(xmin, xmax, tronc, 2.2, 0.5, 2, 
             string("index.dat"), false, true);
  
  if (index.IsComplex())
    {
      cout << "IsComplex incorrect" << endl;
      abort();
    }
  
  int Ng = 200;
  Vector<double> x(Ng+1);
  for (int i = 0; i <= Ng; i++)
    x(i) = -5.0 + 10.0*double(i)/Ng;
  
  nu.Read("sinus_cutoff.dat");
  R2 point; T val, val_ex;
  for (int i = 0; i <= Ng; i++)
    for (int j = 0; j <= Ng; j++)
      {
        point.Init(x(i), x(j));
        index.ComputeIndexAtPoint(point, val);
        if (abs(val - (2.2 + 0.5*nu(i,j)) ) > 1e-5)
          {
            cout << "ComputeIndexAtPoint incorrect" << endl;
            DISP(2.2+0.5*nu(i, j)); DISP(val);
            abort();
          }
      }
  
  xmin.Init(-2.0, -4.0);
  xmax.Init(2.0, 4.0); tronc.Fill(0);
  Matrix<float> nu_f(N+1,N+1);
  for (int i = 0; i <= N; i++)
    for (int j = 0; j <= N; j++)
      {
        double x = xmin(0) + (xmax(0)-xmin(0))*double(i)/N;
        double y = xmin(1) + (xmax(1)-xmin(1))*double(j)/N;
        nu_f(i, j) = cos(2.0*pi_wp*x)*sin(pi_wp*y);
      }
  
  nu_f.Write("index.dat");

  index.Init(xmin, xmax, tronc, 1.8, 0.4, 2, 
             string("index.dat"), false, false);
  
  index.SetPeriodic();
  for (int i = 0; i <= Ng; i++)
    for (int j = 0; j <= Ng; j++)
      {
        point.Init(x(i), x(j));
        val_ex = cos(2.0*pi_wp*x(i))*sin(pi_wp*x(j));
        index.ComputeIndexAtPoint(point, val);
        if (abs(val - (1.8 + 0.4*val_ex) ) > 1e-5)
          {
            cout << "ComputeIndexAtPoint/SetPeriodic incorrect" << endl;
            DISP(1.8 + 0.4*val_ex); DISP(val);
            abort();
          }
      }
  
  Matrix<complex<double> > nu_z(N+1, N+1);
  for (int i = 0; i <= N; i++)
    for (int j = 0; j <= N; j++)
      {
        double x = xmin(0) + (xmax(0)-xmin(0))*double(i)/N;
        double y = xmin(1) + (xmax(1)-xmin(1))*double(j)/N;
        nu_z(i, j) = exp(-1.5*(x*x+y*y))*exp(Iwp*((2.0*pi_wp*x) + (pi_wp*y)));
      }
  
  nu_z.Write("index.dat");

  index.Init(xmin, xmax, tronc, 1.8, 0.4, 2, 
             string("index.dat"), true, true);

  index.SetQuasiPeriodic();
  index.GetInverse();
  if (!index.IsComplex())
    {
      cout << "IsComplex incorrect" << endl;
      abort();
    }
  
  nu_z.Reallocate(Ng+1, Ng+1);
  for (int i = 0; i <= Ng; i++)
    for (int j = 0; j <= Ng; j++)
      {
        point.Init(x(i), x(j));
        double xp = x(i);
        if (x(i) > 2.0)
          xp -= 4.0;
        else if (x(i) < - 2.0)
          xp += 4.0;
        
        double yp = x(j);
        if (x(j) > 4.0)
          yp -= 8.0;
        else if (x(j) < - 4.0)
          yp += 8.0;
                
        double gauss = exp(-1.5*(xp*xp+yp*yp));
        val_ex = 1.0/(1.8 + 0.4*gauss*exp(Iwp*((2.0*pi_wp*x(i)) + (pi_wp*x(j)))));
        if ( (x(i) < xmax(0)-1e-6) && (x(i) > xmin(0)+1e-6)
             && (x(j) < xmax(1)-1e-6) && (x(j) > xmin(1)+1e-6) )
          val_ex = 1.0/1.8;
        
        index.ComputeIndexAtPoint(point, val);
        nu_z(i, j) = val;
        if (abs(val - val_ex ) > 4e-4)
          {
            cout << "ComputeIndexAtPoint/SetPeriodic incorrect" << endl;
            DISP(point); DISP(val_ex); DISP(val);
            abort();
          }
      }
  
  nu_z.Write("quasi_periodic.dat");
  
  Matrix<complex<float> > nu_c(N+1, N+1);
  for (int i = 0; i <= N; i++)
    for (int j = 0; j <= N; j++)
      {
        double x = xmin(0) + (xmax(0)-xmin(0))*double(i)/N;
        double y = xmin(1) + (xmax(1)-xmin(1))*double(j)/N;
        nu_c(i, j) = exp(-2.0*(x*x+y*y))*exp(Iwp*((2.0*pi_wp*x) - (pi_wp*y)));
      }
  
  nu_c.Write("index.dat");
  //tronc.Init(0.01, 0.01);
  
  index.Init(xmin, xmax, tronc, 1.8, 0.4, 2, 
             string("index.dat"), true, false);

  index.SetNoPeriodic();
  index.Mlt(complex<double>(1.2, 0.8) );
  for (int i = 0; i <= Ng; i++)
    for (int j = 0; j <= Ng; j++)
      {
        point.Init(x(i), x(j));
        double xp = x(i);
        double yp = x(j);
        double gauss = exp(-2.0*(xp*xp+yp*yp));
        val_ex = complex<double>(1.2, 0.8)*(1.8 + 0.4*gauss*exp(Iwp*((2.0*pi_wp*x(i)) - (pi_wp*x(j)))));
        if ( (x(i) > xmax(0)+1e-10) && (x(i) < xmin(0)-1e-10)
             && (x(j) > xmax(1)+1e-10) && (x(j) < xmin(1)-1e-10) )
          val_ex = complex<double>(1.2, 0.8)*1.8;
        
        index.ComputeIndexAtPoint(point, val);
        if (abs(val - val_ex ) > 2e-4)
          {
            cout << "ComputeIndexAtPoint/SetNoPeriodic incorrect" << endl;
            DISP(point); DISP(val_ex); DISP(val);
            abort();
          }
      }

  tronc.Init(0.01,0.01);
  index.Init(xmin, xmax, tronc, 1.8, 0.4, 2, 
             string("index.dat"), true, false);
  
  // testing derivatives
  TinyVector<T, 2> grad, deriv;
  R2 point_p;
  Real_wp h = 1e-6; T val_m, val_p;
  for (int i = 0; i <= Ng; i ++)
    for (int j = 0; j <= Ng; j++)
      {
        point.Init(x(i), x(j));
        index.ComputeIndexAtPoint(point, val, grad);
        
        point_p = point; point_p(0) += h;
        index.ComputeIndexAtPoint(point_p, val_p);

        point_p = point; point_p(0) -= h;
        index.ComputeIndexAtPoint(point_p, val_m);
        deriv(0) = (val_p - val_m) / (2.0*h);
        
        point_p = point; point_p(1) += h;
        index.ComputeIndexAtPoint(point_p, val_p);

        point_p = point; point_p(1) -= h;
        index.ComputeIndexAtPoint(point_p, val_m);
        deriv(1) = (val_p - val_m) / (2.0*h);
        if (abs(grad.Distance(deriv)) > 1e-2)
          {
            cout << "Derivatives not correctly computed" << endl;
            DISP(point); DISP(deriv); DISP(grad); abort();
          }
      }
}

template<class T>
void CheckRegularIndex(PhysicalRegularMedia<Dimension3, T>& index)
{  
  // to do
}

template<class T>
void CheckRadialIndex(RadialVaryingMedia<Dimension2, T>& index)
{
  Real_wp rmax = 10.0; int nb_radius = 100;
  string data_file("toto.dat");
  ofstream file_out(data_file.data()); file_out.precision(15);
  for (int i = 0; i < nb_radius; i++)
    {
      Real_wp r = Real_wp(i)*rmax / (nb_radius-1);
      file_out << r << " " << sin(r) << '\n';
    }

  file_out.close();
  
  T offset(0.8), amplitude(0.3); T cte(1);
  index.InitSpline(data_file, offset, amplitude, cte);

  if (index.IsComplex())
    {
      cout << "IsComplex incorrect" << endl;
      abort();
    }

  Real_wp xmin = -6.0, xmax = 6.0, ymin = -6.0, ymax = 6.0;
  int nbx = 201, nby = 200; VectReal_wp x_div, y_div;
  file_out.open("interp.dat");
  Linspace(xmin, xmax, nbx, x_div);
  Linspace(ymin, ymax, nby, y_div);

  for (int i = 0; i < nbx; i++)
    for (int j = 0; j < nby; j++)
      {
        Real_wp x = x_div(i);
        Real_wp y = y_div(j);
        R2 pt(x, y); Real_wp coef;
        
        index.ComputeIndexAtPoint(pt, coef);
        file_out << coef << " ";
        if (j == nby-1)
          file_out << '\n';

        Real_wp r = Norm2(pt);
        Real_wp coef_ref = offset + amplitude*sin(r);
        if (abs(coef-coef_ref) > 1e-4)
          {
            DISP(x); DISP(y); DISP(coef); DISP(coef_ref); DISP(abs(coef-coef_ref));
            cout << "ComputeIndexAtPoint incorrect" << endl;
            abort();
          }

        R2 grad_coef; TinyMatrix<T, Symmetric, 2, 2> hess_coef;
        index.ComputeIndexAtPoint(pt, coef, grad_coef, hess_coef, false);

        Real_wp h = 1e-5;
        Real_wp coef_xp, coef_yp, coef_xm, coef_ym;
        Real_wp r_xp = sqrt((x+h)*(x+h) + y*y);
        Real_wp r_xm = sqrt((x-h)*(x-h) + y*y);
        Real_wp r_yp = sqrt((y+h)*(y+h) + x*x);
        Real_wp r_ym = sqrt((y-h)*(y-h) + x*x);
        
        coef_xp = offset + amplitude*sin(r_xp);
        coef_xm = offset + amplitude*sin(r_xm);
        coef_yp = offset + amplitude*sin(r_yp);
        coef_ym = offset + amplitude*sin(r_ym);

        R2 grad_ref((coef_xp - coef_xm)/(2.0*h), (coef_yp - coef_ym)/(2.0*h));
        //DISP(grad_ref); DISP(grad_coef);DISP(grad_ref.Distance(grad_coef));
        if (grad_ref.Distance(grad_coef) > 1e-4)
          {
            DISP(x); DISP(y); DISP(grad_coef); DISP(grad_ref); DISP(grad_ref.Distance(grad_coef));
            cout << "ComputeIndexAtPoint incorrect" << endl;
            abort();

          }

        grad_coef.Zero();
        index.ComputeIndexAtPoint(pt, coef, grad_coef, hess_coef, true);
        if (grad_ref.Distance(grad_coef) > 1e-4)
          {
            DISP(x); DISP(y); DISP(grad_coef); DISP(grad_ref); DISP(grad_ref.Distance(grad_coef));
            cout << "ComputeIndexAtPoint incorrect" << endl;
            abort();
          }

        Real_wp dcoef_dx2 = (coef_xm + coef_xp - 2.0*coef_ref) / (h*h);
        Real_wp dcoef_dy2 = (coef_ym + coef_yp - 2.0*coef_ref) / (h*h);

        r_xp = sqrt((x+h)*(x+h) + (y+h)*(y+h));
        r_xm = sqrt((x-h)*(x-h) + (y-h)*(y-h));
        r_yp = sqrt((y+h)*(y+h) + (x-h)*(x-h));
        r_ym = sqrt((y-h)*(y-h) + (x+h)*(x+h));
        
        coef_xp = offset + amplitude*sin(r_xp);
        coef_xm = offset + amplitude*sin(r_xm);
        coef_yp = offset + amplitude*sin(r_yp);
        coef_ym = offset + amplitude*sin(r_ym);

        Real_wp dcoef_dxdy = (coef_xp + coef_xm - coef_yp - coef_ym) / (4.0*h*h);
        TinyMatrix<T, Symmetric, 2, 2> hess_ref;
        hess_ref(0, 0) = dcoef_dx2; hess_ref(1, 1) = dcoef_dy2;
        hess_ref(0, 1) = dcoef_dxdy;

        if ( (abs(hess_ref(0, 0) - hess_coef(0, 0)) > 1e-3)
             || (abs(hess_ref(1, 1) - hess_coef(1, 1)) > 1e-3)
             || (abs(hess_ref(0, 1) - hess_coef(0, 1)) > 1e-3))
          {
            DISP(x); DISP(y); DISP(hess_ref); DISP(hess_coef); 
            cout << "ComputeIndexAtPoint incorrect" << endl;
            abort();
          }
        
      }

  file_out.close();

  // testing with inverse
  index.GetInverse();
  for (int i = 0; i < nbx; i++)
    for (int j = 0; j < nby; j++)
      {
        Real_wp x = x_div(i);
        Real_wp y = y_div(j);
        R2 pt(x, y); Real_wp coef;
        
        index.ComputeIndexAtPoint(pt, coef);
        
        Real_wp r = Norm2(pt);
        Real_wp coef_ref = Real_wp(1) / (offset + amplitude*sin(r));
        if (abs(coef-coef_ref) > 1e-4)
          {
            DISP(x); DISP(y); DISP(coef); DISP(coef_ref); DISP(abs(coef-coef_ref));
            cout << "ComputeIndexAtPoint incorrect" << endl;
            abort();
          }
        
        R2 grad_coef; TinyMatrix<T, Symmetric, 2, 2> hess_coef;
        index.ComputeIndexAtPoint(pt, coef, grad_coef, hess_coef, false);

        Real_wp h = 1e-5;
        Real_wp coef_xp, coef_yp, coef_xm, coef_ym;
        Real_wp r_xp = sqrt((x+h)*(x+h) + y*y);
        Real_wp r_xm = sqrt((x-h)*(x-h) + y*y);
        Real_wp r_yp = sqrt((y+h)*(y+h) + x*x);
        Real_wp r_ym = sqrt((y-h)*(y-h) + x*x);
        
        coef_xp = Real_wp(1) / (offset + amplitude*sin(r_xp));
        coef_xm = Real_wp(1) / (offset + amplitude*sin(r_xm));
        coef_yp = Real_wp(1) / (offset + amplitude*sin(r_yp));
        coef_ym = Real_wp(1) / (offset + amplitude*sin(r_ym));

        R2 grad_ref((coef_xp - coef_xm)/(2.0*h), (coef_yp - coef_ym)/(2.0*h));
        //DISP(grad_ref); DISP(grad_coef);DISP(grad_ref.Distance(grad_coef));
        if (grad_ref.Distance(grad_coef) > 1e-4)
          {
            DISP(x); DISP(y); DISP(grad_coef); DISP(grad_ref); DISP(grad_ref.Distance(grad_coef));
            cout << "ComputeIndexAtPoint incorrect" << endl;
            abort();

          }

        grad_coef.Zero();
        index.ComputeIndexAtPoint(pt, coef, grad_coef, hess_coef, true);
        if (grad_ref.Distance(grad_coef) > 1e-4)
          {
            DISP(x); DISP(y); DISP(grad_coef); DISP(grad_ref); DISP(grad_ref.Distance(grad_coef));
            cout << "ComputeIndexAtPoint incorrect" << endl;
            abort();
          }

        Real_wp dcoef_dx2 = (coef_xm + coef_xp - 2.0*coef_ref) / (h*h);
        Real_wp dcoef_dy2 = (coef_ym + coef_yp - 2.0*coef_ref) / (h*h);

        r_xp = sqrt((x+h)*(x+h) + (y+h)*(y+h));
        r_xm = sqrt((x-h)*(x-h) + (y-h)*(y-h));
        r_yp = sqrt((y+h)*(y+h) + (x-h)*(x-h));
        r_ym = sqrt((y-h)*(y-h) + (x+h)*(x+h));
        
        coef_xp = Real_wp(1) / (offset + amplitude*sin(r_xp));
        coef_xm = Real_wp(1) / (offset + amplitude*sin(r_xm));
        coef_yp = Real_wp(1) / (offset + amplitude*sin(r_yp));
        coef_ym = Real_wp(1) / (offset + amplitude*sin(r_ym));
        
        Real_wp dcoef_dxdy = (coef_xp + coef_xm - coef_yp - coef_ym) / (4.0*h*h);
        TinyMatrix<T, Symmetric, 2, 2> hess_ref;
        hess_ref(0, 0) = dcoef_dx2; hess_ref(1, 1) = dcoef_dy2;
        hess_ref(0, 1) = dcoef_dxdy;

        if ( (abs(hess_ref(0, 0) - hess_coef(0, 0)) > 1e-2)
             || (abs(hess_ref(1, 1) - hess_coef(1, 1)) > 1e-2)
             || (abs(hess_ref(0, 1) - hess_coef(0, 1)) > 1e-2))
          {
            DISP(x); DISP(y); DISP(hess_ref); DISP(hess_coef); 
            cout << "ComputeIndexAtPoint incorrect" << endl;
            abort();
          }
        
      }

  // testing ComputeIndex
  index.GetInverse();
  
  Mesh<Dimension2> mesh;
  mesh.CreateRegularMesh(R2(xmin, ymin), R2(xmax, ymax), TinyVector<int, 2>(4, 4),
                         1, TinyVector<int, 4>(1, 1, 1, 1), 0); 

  IVect ElementRho(mesh.GetNbElt()); ElementRho.Fill();
  Vector<VectR2> Points(mesh.GetNbElt());
  Vector<Vector<T> > rho(mesh.GetNbElt());
  Vector<Vector<TinyVector<T, 2> > > grad_rho(mesh.GetNbElt());
  Vector<Vector<TinyMatrix<T, Symmetric, 2, 2> > > hess_rho(mesh.GetNbElt());
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      int nb_quad = rand()%10+1;
      
      Points(i).Reallocate(nb_quad);
      rho(i).Reallocate(nb_quad);
      grad_rho(i).Reallocate(nb_quad);
      hess_rho(i).Reallocate(nb_quad);
      for (int j = 0; j < nb_quad; j++)
        {
          Real_wp x = xmin + Real_wp(rand()) / RAND_MAX * (xmax-xmin);
          Real_wp y = ymin + Real_wp(rand()) / RAND_MAX * (ymax-ymin);
          Points(i)(j).Init(x, y);
        }
    }
    
  index.ComputeIndex(mesh, ElementRho, 1, Points, true, true, rho, grad_rho, hess_rho);
  
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      DISP(i);
      DISP(Points(i));
      DISP(rho(i)); DISP(grad_rho(i));
      for (int j = 0; j < Points(i).GetM(); j++)
        DISP(hess_rho(i)(j));
    }
  
}


template<class T>
void CheckVaryingIndex(PhysicalVaryingMedia<Dimension2, T>& index)
{  
  // checking constant media
  index.InitConstant(Complex_wp(0, 0));
  if (index.IsComplex())
    {
      cout << "IsComplex incorrect" << endl;
      abort();
    }
  
  if (!index.IsZero())
    {
      cout << "IsZero incorrect" << endl;
      abort();
    }
  
  if (index.GetMediaType() != index.CONSTANT)
    {
      cout << "GetMediaType incorrect" << endl;
      abort();
    }
  
  index.InitConstant(Complex_wp(1.3, 0.1));
  if (index.IsVarying())
    {
      cout << "IsVarying incorrect" << endl;
      abort();
    }
  
  if (!index.IsComplex())
    {
      cout << "IsComplex incorrect" << endl;
      abort();
    }
  
  if (index.IsZero())
    {
      cout << "IsZero incorrect" << endl;
      abort();
    }
  
  if (index.GetConstant() != Complex_wp(1.3, 0.1))
    {
      cout << "GetConstant incorrect" << endl;
      abort();
    }  
  
  // testing sinusoidal media
  T offset(2.2), amplitude(0.5);
  R2 xmin(-2, -1), xmax(3, 5), coef(0.01, 0.01), kwave(pi_wp, 2.0*pi_wp);
  TinyVector<R2, 2> enveloppe; enveloppe(0) = xmin; enveloppe(1) = xmax;
  index.InitSinus(xmin, xmax, coef, offset, amplitude, kwave);
  
  if (index.GetMediaType() != index.SINUSOIDE)
    {
      cout << "GetMediaType incorrect" << endl;
      abort();
    }
  
  if (index.IsComplex())
    {
      cout << "IsComplex incorrect" << endl;
      abort();
    }
  
  // testing regular media
  int N = 200;
  Matrix<complex<double> > nu(N+1, N+1);
  for (int i = 0; i <= N; i++)
    {
      for (int j = 0; j <= N; j++)
        {
          double x = xmin(0) + (xmax(0)-xmin(0))*double(i)/N;
          double y = xmin(1) + (xmax(1)-xmin(1))*double(j)/N;
          nu(i, j) = exp(Iwp*((2.0*pi_wp*x) - (pi_wp*y)));
        }
    }
  
  nu.Write("regular.dat");
  
  coef.Fill(0);
  index.InitRandom(xmin, xmax, coef, offset, amplitude, 2,
                   string("regular.dat"), true, true);
  
  index.index_regular.SetPeriodic();
  
  if (!index.IsComplex())
    {
      cout << "IsComplex incorrect" << endl;
      abort();
    }
  
  if (index.IsZero())
    {
      cout << "IsZero incorrect" << endl;
      abort();
    }
  
  if (index.GetMediaType() != index.REGULAR_GRID)
    {
      cout << "GetMediaType incorrect" << endl;
      abort();
    }
  
  // quadrature points on quadrangles
  VectR2 PointsQuad; VectReal_wp Weights;
  int r = 4;
  QuadrangleQuadrature::ConstructQuadrature(r, PointsQuad, Weights);
  
  // creating a mesh
  Mesh<Dimension2> mesh;
  mesh.CreateRegularMesh(R2(-5, -5), R2(5, 5), TinyVector<int, 2>(30, 30),
                         1, TinyVector<int, 4>(1, 1, 1, 1), mesh.QUADRILATERAL_MESH);
  
  // and points
  VectR2 s; SetPoints<Dimension2> PtsElem;
  Vector<VectR2> Points(mesh.GetNbElt());
  R2 pt_glob;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      mesh.GetVerticesElement(i, s);
      mesh.FjElemNodal(s, PtsElem, mesh, i);
      Points(i).Reallocate(PointsQuad.GetM());
      for (int j = 0; j < PointsQuad.GetM(); j++)
        {
          mesh.Fj(s, PtsElem, PointsQuad(j), pt_glob, mesh, i);
          Points(i)(j) = pt_glob;
        }
    }
  
  // computing values of index on these points
  index.Reallocate(mesh.GetNbElt(), true);
  for (int i = 0; i < mesh.GetNbElt(); i++)
    index.ReallocateOnElement(i, Points(i).GetM(), true);
  
  IVect ElementRho(mesh.GetNbElt()); ElementRho.Fill();
  index.index_regular.ComputeIndex(mesh, ElementRho, 1, Points, true, false,
                                   index.GetValue(), index.GetGradient(), index.GetHessian());
  
  // checking values
  T val_ex, val_p, val_m; TinyVector<T, 2> grad;
  double h = 1e-7;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    for (int j = 0; j < PointsQuad.GetM(); j++)
      {
        double x = Points(i)(j)(0);
        double y = Points(i)(j)(1);
        bool inside = PointInsideBoundingBox(R2(x, y), enveloppe);
        if (inside)
          {
            val_ex = offset + amplitude*exp(Iwp*((2.0*pi_wp*x) - (pi_wp*y)));
            
            if (abs(val_ex - index.GetCoefficient(i, j)) > 1e-4)
              {
                cout << "GetCoefficient/GetValue incorrect" << endl;
                DISP(i); DISP(j); DISP(x); DISP(y); DISP(val_ex); DISP(index.GetCoefficient(i, j));
                abort();
              }
            
            val_m = offset + amplitude*exp(Iwp*((2.0*pi_wp*(x-h)) - (pi_wp*y)));
            val_p = offset + amplitude*exp(Iwp*((2.0*pi_wp*(x+h)) - (pi_wp*y)));
            grad(0) = (val_p - val_m) / (2.0*h);

            val_m = offset + amplitude*exp(Iwp*((2.0*pi_wp*x) - (pi_wp*(y-h))));
            val_p = offset + amplitude*exp(Iwp*((2.0*pi_wp*x) - (pi_wp*(y+h))));
            grad(1) = (val_p - val_m) / (2.0*h);
            
            if (abs(grad.Distance(index.GetCoefGradient(i, j))) > 5e-4)
              {
                DISP(i); DISP(j); DISP(x); DISP(y); 
                DISP(grad); DISP(index.GetCoefGradient(i, j));
                cout << "GetCoefGradient incorrect" << endl;
                abort();
              }
          }
      }
  
  // you can also give directly coefficients
  index.SetMediaType(index.USER);
  
  if (index.GetMediaType() != index.USER)
    {
      cout << "GetMediaType incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < mesh.GetNbElt(); i++)
    for (int j = 0; j < PointsQuad.GetM(); j++)
      {
        double x = Points(i)(j)(0);
        double y = Points(i)(j)(1);
        val_ex = x*x + 2.0 + 1.5*x*y + pow(y, 3.0);
        grad(0) = 2.0*x + 1.5*y;
        grad(1) = 1.5*x + 3.0*y*y;
        index.SetCoefficient(i, j, val_ex);
        index.SetCoefGradient(i, j, grad);
      }
  
  index.GetInverse();
  
  for (int i = 0; i < mesh.GetNbElt(); i++)
    for (int j = 0; j < PointsQuad.GetM(); j++)
      {
        double x = Points(i)(j)(0);
        double y = Points(i)(j)(1);
        val_ex = 1.0/(x*x + 2.0 + 1.5*x*y + pow(y, 3.0));
        if (abs(val_ex - index.GetCoefficient(i, j)) > 1e-4)
          {
            cout << "GetInverse incorrect" << endl;
            abort();
          }
        
        x -= h;
        val_m = 1.0/(x*x + 2.0 + 1.5*x*y + pow(y, 3.0));
        x += 2.0*h;
        val_p = 1.0/(x*x + 2.0 + 1.5*x*y + pow(y, 3.0));
        grad(0) = (val_p - val_m) / (2.0*h);
        
        x -= h; y -= h;
        val_m = 1.0/(x*x + 2.0 + 1.5*x*y + pow(y, 3.0));
        y += 2.0*h;
        val_p = 1.0/(x*x + 2.0 + 1.5*x*y + pow(y, 3.0));
        grad(1) = (val_p - val_m) / (2.0*h);
        
        //DISP(i); DISP(j); DISP(grad); DISP(index.GetCoefGradient(i, j));
        if (abs(grad.Distance(index.GetCoefGradient(i, j))) > 1e-4*abs(val_ex))
          {
            cout << "GetInverse incorrect" << endl;
            abort();
          }
      }
  
        
}

template<class T>
void CheckVaryingIndex(PhysicalVaryingMedia<Dimension3, T>& index)
{  
  // to do
}

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x)
{
  int type = rand()%3;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void CheckScalarIndex(ScalarPhysicalIndice<Dimension2, T>& index)
{
  T alpha, zero, one, beta;
  SetComplexZero(zero);
  SetComplexOne(one);
  
  GetRandNumber(alpha);
  GetRandNumber(beta);
  index.SetConstant(alpha);
  
  if (index.IsZero() || index.IsVarying() )
    {
      cout << "SetConstant incorrect" << endl;
      abort();
    }
  
  if (index.GetConstant() != alpha)
    {
      cout << "GetConstant incorrect" << endl;
      abort();
    }
  
  index.Mlt(beta);
  index.GetInverse();
  if (abs(index.GetConstant() - one/(alpha*beta)) > threshold )
    {
      cout << "Mlt/GetInverse incorrect" << endl;
      abort();
    }
    
  index.Zero();
  if (!index.IsZero() || index.IsVarying() )
    {
      cout << "Zero incorrect" << endl;
      abort();
    }
  
  if (index.GetConstant() != zero)
    {
      cout << "Zero incorrect" << endl;
      abort();
    }
  
  index.SetIdentity();
  if (index.IsZero() || index.IsVarying() )
    {
      cout << "SetIdentity incorrect" << endl;
      abort();
    }
  
  if (index.GetConstant() != one)
    {
      cout << "SetIdentity incorrect" << endl;
      abort();
    }
  
  // testing variable index
  index.SetMediaType(PhysicalVaryingMedia<Dimension2, T>::SINUSOIDE);
  if ( !index.IsVarying() || index.GetNbVaryingMedia() != 1)
    {
      cout << "SetMediaType incorrect" << endl;
      abort();
    }
  
  Vector<PhysicalVaryingMedia<Dimension2, T>* > rho(1);
  IVect num(1);
  int nb = 0;
  index.GetVaryingMedia(nb, rho, num);
  
  R2 xmin, xmax, coef, kwave;
  xmin.Init(-2.0, -3.0);
  xmax.Init(3.0, 4.0);
  kwave.Init(pi_wp, 2.0*pi_wp);
  rho(0)->index_sinus.Init(xmin, xmax, coef, alpha, beta, kwave);
  
  R2 point(0.8, 0.7);
  T val_inv, val; TinyVector<T, 2> grad, grad_inv;
  rho(0)->index_sinus.ComputeIndexAtPoint(point, val, grad);
  
  index.GetInverse();
  rho(0)->index_sinus.ComputeIndexAtPoint(point, val_inv, grad_inv);
  
  if ( (abs(val - one/val_inv) > threshold)
       || Norm2(grad_inv + square(val_inv)*grad) > threshold)
    {
      cout << "Inverse incorrect" << endl;
      abort();
    }
}  

template<class T>
void CheckScalarIndex(ScalarPhysicalIndice<Dimension3, T>& index)
{
}

class FakeObject
{
public :
  IVect ElementRho;
};

template<class T>
void CheckElasticIndex(ElasticPhysicalIndice<Dimension2, 2, T>& index)
{
  TinySymmetricTensor<T, 2> C, D;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      GetRandNumber(C(i, j));
  
  C(0, 0) = 2.0; C(0, 1) = 1.7; C(0, 2) = 5.1;
  C(1, 0) = -0.8; C(1, 1) = 3.2; C(1, 2) = 1.2;
  C(2, 0) = 2.5; C(2, 1) = 3.1; C(2, 2) = 4.8;
  
  index.SetConstant(C);
  
  if (index.GetAnisotropyType() != index.ANISOTROPE)
    {
      cout << "GetAnisotropyType/SetConstant incorrect" << endl;
      abort();
    }

  //DISP(index.GetConstant());
  
  TinyVector<T, 4> dU, dV, dU2;
  for (int i = 0; i < 4; i++)
    GetRandNumber(dU(i));  
  
  dU(2) = dU(1);

  VarPhysicalProblem var;
  index.MltVector(var, 0, 0, dU, dV);
  
  index.GetInverse();
  
  index.MltVector(var, 0, 0, dV, dU2);
  for (int i = 0; i < 4; i++)
    {
      if (abs(dU(i) - dU2(i)) > threshold)
        {
          cout << "MltVector/GetInverse incorrect" << endl;
          abort();
        }
    }      
  
  T zero; SetComplexZero(zero);
  D(0, 0) = 2.0; D(0, 1) = zero; D(0, 2) = 5.1;
  D(1, 0) = zero; D(1, 1) = 3.2; D(1, 2) = zero;
  D(2, 0) = 2.5; D(2, 1) = zero; D(2, 2) = 4.8;

  index.SetConstant(D);
  
  if (index.GetAnisotropyType() != index.ORTHOTROPE)
    {
      cout << "GetAnisotropyType/SetConstant incorrect" << endl;
      abort();
    }
  
  TinyVector<TinyVector<T, 2>, 2> du, dv, dv_ref;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      {
        GetRandNumber(du(i)(j));
        dU(i*2 + j) = du(i)(j);
      }
  
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        for (int l = 0; l < 2; l++)
          dv_ref(i)(j) += index(i, j, k, l)*du(k)(l);
  
  index.MltMatrix(var, 0, 0, du, dv);
  //D.MltOrthotrope(du, dv);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      if (abs(dv_ref(i)(j) - dv(i)(j)) > threshold)
        {
          cout << "MltMatrix incorrect" << endl;
          abort();
        }
  
  FillZero(dV);
  index.MltVector(var, 0, 0, dU, dV);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      if (abs(dv_ref(i)(j) - dV(i*2 + j)) > threshold)
        {
          cout << "MltVector incorrect" << endl;
          abort();
        }

  index.SetConstant(C);
  
  FillZero(dv_ref);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        for (int l = 0; l < 2; l++)
          dv_ref(i)(j) += index(i, j, k, l)*du(k)(l);
  
  index.MltMatrix(var, 0, 0, du, dv);
  //C.Mlt(du, dv);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      if (abs(dv_ref(i)(j) - dv(i)(j)) > threshold)
        {
          cout << "MltMatrix incorrect" << endl;
          abort();
        }
  
  FillZero(dV);
  index.MltVector(var, 0, 0, dU, dV);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      if (abs(dv_ref(i)(j) - dV(i*2 + j)) > threshold)
        {
          cout << "MltVector incorrect" << endl;
          abort();
        }
  
  C.FillIsotrope(2.0, 3.0);
  
  C.GetInverse();
}

template<class T>
void CheckElasticIndex(ElasticPhysicalIndice<Dimension3, 3, T>& index)
{
  TinySymmetricTensor<T, 3> C, D;
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      GetRandNumber(C(i, j));
  
  index.SetConstant(C);
  
  if (index.GetAnisotropyType() != index.ANISOTROPE)
    {
      cout << "GetAnisotropyType/SetConstant incorrect" << endl;
      abort();
    }

  //DISP(index.GetConstant());
  
  TinyVector<T, 9> dU, dV, dU2;
  for (int i = 0; i < 9; i++)
    GetRandNumber(dU(i));  
  
  dU(3) = dU(1);   dU(6) = dU(2);   dU(7) = dU(5);
  
  VarPhysicalProblem var;
  index.MltVector(var, 0, 0, dU, dV);
  
  index.GetInverse();
  
  index.MltVector(var, 0, 0, dV, dU2);
  for (int i = 0; i < 9; i++)
    {
      if (abs(dU(i) - dU2(i)) > threshold)
        {
          cout << "MltVector/GetInverse incorrect" << endl;
          abort();
        }
    }      
  
  T zero; SetComplexZero(zero);
  D = C;
  D(0, 1) = D(0, 2) = D(0, 4) = zero;
  D(3, 1) = D(3, 2) = D(3, 4) = zero;
  D(5, 1) = D(5, 2) = D(5, 4) = zero;
  D(1, 2) = zero; D(1, 4) = zero; D(2, 4) = zero;
  //DISP(D);
  
  index.SetConstant(D);
  
  if (index.GetAnisotropyType() != index.ORTHOTROPE)
    {
      cout << "GetAnisotropyType/SetConstant incorrect" << endl;
      abort();
    }
  
  TinyVector<TinyVector<T, 3>, 3> du, dv, dv_ref;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {
        GetRandNumber(du(i)(j));
        dU(i*3 + j) = du(i)(j);
      }
  
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          dv_ref(i)(j) += index(i, j, k, l)*du(k)(l);
  
  index.MltMatrix(var, 0, 0, du, dv);
  //D.MltOrthotrope(du, dv);
  //DISP(du); DISP(dv); DISP(dv_ref);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (abs(dv_ref(i)(j) - dv(i)(j)) > threshold)
        {
          cout << "MltMatrix incorrect" << endl;
          abort();
        }
  
  FillZero(dV);
  index.MltVector(var, 0, 0, dU, dV);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (abs(dv_ref(i)(j) - dV(i*3 + j)) > threshold)
        {
          cout << "MltVector incorrect" << endl;
          abort();
        }

  index.SetConstant(C);
  
  FillZero(dv_ref);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          dv_ref(i)(j) += index(i, j, k, l)*du(k)(l);
  
  index.MltMatrix(var, 0, 0, du, dv);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (abs(dv_ref(i)(j) - dv(i)(j)) > threshold)
        {
          cout << "MltMatrix incorrect" << endl;
          abort();
        }
  
  FillZero(dV);
  index.MltVector(var, 0, 0, dU, dV);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (abs(dv_ref(i)(j) - dV(i*3 + j)) > threshold)
        {
          cout << "MltVector incorrect" << endl;
          abort();
        }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  {
    RadialVaryingMedia<Dimension2, Real_wp> rho;
    CheckRadialIndex(rho);
  }
  
  // testing RectangleCutOff  
  //srand(time(NULL));
  {
    RectangleCutOff<Dimension2, Real_wp> cut;
    CheckRectangleCutOff(cut);
  }

  {
    RectangleCutOff<Dimension3, Real_wp> cut;
    CheckRectangleCutOff(cut);
  }
  
  // testing PhysicalSinusMedia
  {
    PhysicalSinusoidalMedia<Dimension2, Complex_wp> index;
    CheckSinusIndex(index);
  }

  {
    PhysicalSinusoidalMedia<Dimension3, Real_wp> index;
    CheckSinusIndex(index);
  }
  
  // testing PhysicalRegularMedia
  {
    PhysicalRegularMedia<Dimension2, Complex_wp> index;
    CheckRegularIndex(index);
  }
  
  {
    PhysicalRegularMedia<Dimension3, Complex_wp> index;
    CheckRegularIndex(index);
  }
  
  // testing PhysicalVaryingMedia
  {
    PhysicalVaryingMedia<Dimension2, Complex_wp> index;
    CheckVaryingIndex(index);
  }

  {
    PhysicalVaryingMedia<Dimension3, Complex_wp> index;
    CheckVaryingIndex(index);
  }
  
  // testing ScalarPhysicalIndice
  {
    ScalarPhysicalIndice<Dimension2, Real_wp> index;
    CheckScalarIndex(index);
  }

  {
    ScalarPhysicalIndice<Dimension2, Complex_wp> index;
    CheckScalarIndex(index);
  }

  {
    ScalarPhysicalIndice<Dimension3, Real_wp> index;
    CheckScalarIndex(index);
  }

  {
    ScalarPhysicalIndice<Dimension3, Complex_wp> index;
    CheckScalarIndex(index);
  }

  // testing ElasticPhysicalIndice
  {
    ElasticPhysicalIndice<Dimension2, 2, Real_wp> index;
    CheckElasticIndex(index);
  }

  {
    ElasticPhysicalIndice<Dimension3, 3, Real_wp> index;
    CheckElasticIndex(index);
  }


  cout << "All tests passed successfully" << endl;
  
  return 0;
}
