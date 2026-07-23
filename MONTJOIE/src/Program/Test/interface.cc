#include "Montjoie.hxx"

using namespace Montjoie;

class LeastSquaresSplineInterp : public VirtualLeastSquaresFunction<Real_wp>
{
  VectR3 all_points, points_init;
  
public:
  void Init(const VectR3& all_pts, const VectR3& pts_init)
  {
    all_points = all_pts;
    points_init = pts_init;
    m_ = all_pts.GetM();
    n_ = pts_init.GetM();
  }

  void FindInitGuess(VectReal_wp& param)
  {
    param.Reallocate(points_init.GetM());
    for (int i = 0; i < points_init.GetM(); i++)
      param(i) = points_init(i)(2);
  }

  void EvaluateF(const VectReal_wp& z, VectReal_wp& f)
  {
    SplineInterpolation<Real_wp> spline;
    VectReal_wp x(points_init.GetM());
    for (int i = 0; i < points_init.GetM(); i++)
      x(i) = points_init(i)(0);
    
    spline.Init(x, z);

    f.Reallocate(this->m_);
    for (int i = 0; i < this->m_; i++)
      f(i) = spline.Evaluate(all_points(i)(0)) - all_points(i)(2);
  }
  
  void EvaluateF0(const VectReal_wp& z, VectReal_wp& f)
  {
    SplineInterpolation<Real_wp> spline;
    VectReal_wp x(points_init.GetM());
    for (int i = 0; i < points_init.GetM(); i++)
      x(i) = points_init(i)(0);
    
    spline.Init(x, z);

    f.Reallocate(this->m_);
    for (int i = 0; i < this->m_; i++)
      f(i) = spline.Evaluate(all_points(i)(0));
  }
  
  void EvaluateJacobian(const VectReal_wp& sol, VectReal_wp& f, Matrix<Real_wp, General, ColMajor>& fjac)
  {
    VectReal_wp z(sol);
    VectReal_wp fp;
    Real_wp h = 1e-6;
    EvaluateF(z, f);
    fjac.Reallocate(this->m_, this->n_);
    for (int j = 0; j < this->n_; j++)
      {
        z(j) += h;
        EvaluateF(z, fp);
        for (int i = 0; i < this->m_; i++)
          fjac(i, j) = (fp(i) - f(i)) / h;

        z(j) -= h;
      }
  }
  
};

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  Real_wp x1 = -10, y1 = -10, x2 = 10, y2 = 10, z1 = 0, z2 = -20;
  int Nx = 676, Ny = 676, Nz = 210;

  VectReal_wp xinterp, zinterp;
  Linspace(x1, x2, Nx, xinterp);
  Linspace(z1, z2, Nz, zinterp);
  Real_wp dz = zinterp(1) - zinterp(0); DISP(dz);

  Vector<float> v(Nx*Ny*Nz);
  v.Read("vel.H@", false);

  if (false)
  { 
  int Nxy = Nx*Ny;
  ofstream file_out("sel.geo");
  int npt = 0;
  for (int i = 1; i < Nx-1; i++)
      for (int j = 1; j < Ny-1; j++)
        for (int k = 1; k < Nz-1; k++)
          {
            int ind = j*Nx + i; int num = k*Nxy + ind;
            if (v(num) >= 4400.0)
              {
                int nb_neigh = 0;
                if (v(num-1) >= 4400.0)
                  nb_neigh++;

                if (v(num+1) >= 4400.0)
                  nb_neigh++;
                
                if (v(num-Nx) >= 4400.0)
                  nb_neigh++;
                                
                if (v(num+Nx) >= 4400.0)
                  nb_neigh++;
                
                if (v(num+Nxy) >= 4400.0)
                  nb_neigh++;

                if (v(num-Nxy) >= 4400.0)
                  nb_neigh++;

                if (nb_neigh <= 5)
                  {
                    Real_wp x = xinterp(i), y = xinterp(j), z = zinterp(k);
                    if ((z > -5.0) && (z >= -4.96))
                      file_out<<"Point(" + to_str(npt+1) + ") = { " << x << ", " << y << ", " << z << "};\n";
                    
                    npt++;
                  }
              }
          }

  file_out.close();
  }
  
  //return FinalizeMontjoie();

  
  Matrix<int> level_water(676, 676);
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      {
        int ind = j*Nx + i;
        int k = 0;
        while (v(k*Nx*Ny + ind) <= 1501.0)
          k++;

        //DISP(i); DISP(j); DISP(v(k*Nx*Ny + ind)); DISP(v((k-1)*Nx*Ny + ind));
        level_water(i, j) = k;
      }

  level_water.WriteText("lvl.dat");
  v.Clear();

  ofstream file_out("contour.geo");
  file_out.precision(7);
  file_out << "lc = 200;\n" << endl;
  int np = 1;
  for (int j = 0; j < level_water.GetN(); j++)
    {
      int lp = level_water(0, j);
      Real_wp y = 20.0*j;
      for (int i = 1; i < level_water.GetM(); i++)
        {
          if (level_water(i, j) != lp)
            {
              Real_wp x = 20.0*i - 10.0;
              Real_wp z = -20.0*lp;
              if (level_water(i, j) < lp)
                z = 20.0*level_water(i, j);
              
              file_out << "Point(" << np << ") = {" << x << ", " << y << ", " << z << "};\n";
              np++;
            }

          lp = level_water(i, j);
        }
    }
    
  /*
  int nb_splines = 10, nb_pts = 13;
  VectR3 points(2*676); VectR3 pts_init;
  int num0 = 0;
  if (argc > 1)
    num0 = atoi(argv[1]);

  DISP(num0);
  //for (int k = num0; k <= nb_splines; k++)
    {
      int j = num0;
      //if (k < nb_splines)
      //j = (676*k)/nb_splines;

      DISP(j);
      Real_wp x = 20.0*j;
      for (int i = 0; i < 676; i++)
        {
          Real_wp y = 20.0*i;
          Real_wp z = -20.0*level_water(j, i);
          points(2*i).Init(x, y, z);
          points(2*i+1).Init(x, y, z+20.0);
        }

      pts_init.Reallocate(nb_pts);
      for (int i = 0; i < nb_pts; i++)
        {
          int iv = 675;
          if (i < nb_pts-1)
            iv = (676*i)/(nb_pts-1);

          Real_wp y = 20.0*iv;
          pts_init(i).Init(x, y, -20.0*level_water(j, iv)+10.0);
        }

      DISP(pts_init);

      ofstream file_out("points.dat"); file_out.precision(15);
      for (int i = 0; i < points.GetM(); i++)
        file_out << points(i)(1) << " " << points(i)(2) << endl;
      
      file_out.close(); file_out.open("pts_init.dat");
      for (int i = 0; i < pts_init.GetM(); i++)
        file_out << pts_init(i)(1) << " " << pts_init(i)(2) << endl;

      file_out.close();
      exit(0);
          points.Reallocate(101);
      for (int i = 0 ; i < 101; i++)
        points(i).Init(0.032*i, 0.0, sin(0.032*i));

      pts_init.Reallocate(7);
      pts_init(0).Init(0.0, 0.0, 0.0);
      pts_init(1).Init(0.48, 0.0, 0.42);
      pts_init(2).Init(0.96, 0.0, 0.83);
      pts_init(3).Init(1.44, 0.0, 1.05);
      pts_init(4).Init(1.92, 0.0, 0.88);
      pts_init(5).Init(2.4, 0.0, 0.61);
      pts_init(6).Init(3.2, 0.0, 0.16);
          LeastSquaresSplineInterp var;
      var.Init(points, pts_init);

      VectReal_wp sol;
      SolveLeastSquaresLvm(var, sol, 1e-12, -1, 5000, 2);
      DISP(sol);
      VectReal_wp f;
      var.EvaluateF(sol, f);  DISP(f);
      var.EvaluateF0(sol, f);
      f.Write("fapp.dat");

      pts_init(0).Init(0.0, 0.0, -118.0);
      pts_init(1).Init(850.0, 0.0, -120.0);
      pts_init(2).Init(3000.0, 0.0, -128.0);
      pts_init(3).Init(5070.0, 0.0, -120.0);
      pts_init(4).Init(6700.0, 0.0, -102.0);
      pts_init(5).Init(8250.0, 0.0, -120.0);
      pts_init(6).Init(9430.0, 0.0, -140.0);
      pts_init(7).Init(9850.0, 0.0, -160.0);
      pts_init(8).Init(10270.0, 0.0, -180.0);
      pts_init(9).Init(10710.0, 0.0, -200.0);
      pts_init(10).Init(11350.0, 0.0, -220.0);
      pts_init(11).Init(12690.0, 0.0, -240.0);
      pts_init(12).Init(13500.0, 0.0, -252.0);
      var.Init(points, pts_init);
      var.FindInitGuess(sol);
      
      var.EvaluateF0(sol, f);
      f.Write("fapp_init.dat");

      SolveLeastSquaresLvm(var, sol, 1e-12, -1, 5000, 2);
      DISP(sol);
      var.EvaluateF0(sol, f);
      f.Write("fapp_sol.dat");
      
      exit(0);
    }
  */
  //BSplineInterpolation<Real_wp> spline;
  //int order = 4; VectReal_wp knots; Linspace(x1, x2, 5, knots);
  //DISP(knots);
  /*
  for (int i = 0; i < knots.GetM(); i++)
    {
      int j = (Ny-1)*Real_wp(i)/(knots.GetM()-1);
      DISP(j);
      VectReal_wp xdiv(Nx), ydiv(Nx);
      for (int k = 0; k < Nx; k++)
        {
          xdiv(k) = xinterp(k);
          ydiv(k) = zinterp(level_water(k, j)) + 0.5*dz;
        }

      xdiv.Write("xd" + to_str(j) + ".dat");
      ydiv.Write("yd" + to_str(j) + ".dat");

      spline.Init(xdiv, ydiv, order, knots);
      VectReal_wp yinterp(Nx);
      for (int k = 0; k < Nx; k++)
        yinterp(k) = spline.Evaluate(xinterp(k));
      
      xinterp.Write("x" + to_str(j) + ".dat");
      yinterp.Write("y" + to_str(j) + ".dat");

    }
  */

  // filtering with FFT
  /*int Ny2 = Ny/2 + 1;
  Vector<double> vr(Nx*Ny); VectComplex_wp vhat(Nx*Ny2);
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      vr(i*Ny + j) = level_water(i, j);

  FftRealInterface fft;
  fft.Init(Nx, Ny);

  fft.ApplyForward(vr, vhat);

  vr.WriteText("vr.dat");
  vhat.Write("vhat.dat");

  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny2; j++)
      {
        Real_wp coef = 0.0;
        if ((i < 10) || (i >= Nx-10))
          coef = 1.0;
        else if (i < 20)
          coef = Real_wp(20-i) / 10.0;
        else if (i >= Nx-20)
          coef = Real_wp(i-(Nx-21))/10.0;

        Real_wp coefy = 0.0;
        if (j < 10)
          coefy = 1.0;
        else if (j < 20)
          coefy = Real_wp(20-j) / 10.0;

        coef *= coefy;
        vhat(i*Ny2 + j) *= coef;
      }

  vhat.Write("vhat2.dat");

  fft.ApplyInverse(vhat, vr);
  vr.WriteText("vr2.dat");
  */
  // 2-D polynoms
  /*int order = 4;
  int p = (order+1)*(order+2)/2;
  Matrix<Real_wp> VDM(Nx*Ny, p);
  VectReal_wp rhs(Nx*Ny);
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      {
        Real_wp x = xinterp(i)/10, y = xinterp(j)/10;
        int num = 0;
        for (int k = 0; k <= order; k++)
          for (int l = 0; l <= order-k; l++)
            VDM(i*Ny+j, num++) = pow(x, k)*pow(y, l);

        rhs(i*Ny+j) = level_water(i, j);
      }

  VectReal_wp tau;
  GetQR(VDM, tau);
  SolveQR(VDM, tau, rhs);

  DISP(rhs);

  Matrix<Real_wp> vr(Nx, Ny);
  Real_wp err_max = 0;
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      {
        Real_wp x = xinterp(i)/10, y = xinterp(j)/10;
        Real_wp vloc = 0; int num = 0;
        for (int k = 0; k <= order; k++)
          for (int l = 0; l <= order-k; l++)
            vloc += pow(x, k)*pow(y, l)*rhs(num++);

        vr(i, j) = vloc;
        err_max = max(err_max, abs(vloc - level_water(i, j)));
      }

  DISP(err_max);
  vr.WriteText("vr.dat");

  //
  ofstream file_out("water.geo"); file_out.precision(5);
  int nb_spline = 5; int nb_pts = 10;
  int npt = 0;
  for (int i = 0; i < nb_spline; i++)
    {
      Real_wp x = x1 + (x2-x1)*Real_wp(i) / (nb_spline-1); x = x/10;
      for (int j = 0; j < nb_pts; j++)
        {
          Real_wp y = y1 + (y2-y1)*Real_wp(j) / (nb_pts-1); y = y/10;
          Real_wp vloc = 0; int num = 0;
          for (int k = 0; k <= order; k++)
            for (int l = 0; l <= order-k; l++)
              vloc += pow(x, k)*pow(y, l)*rhs(num++);

          file_out << " Point(" << npt+1 << ") = {" << x*10 << ", " << y*10 << ", " << dz*vloc << ", lc};\n";
          npt++;
        }

      file_out << "Spline(" << i+1 << " " <<") = {" << npt-nb_pts+1 << ":" << npt << "}\n";
      file_out << "\n" << endl;
    }

  file_out.close();
  */
  return FinalizeMontjoie();
}
