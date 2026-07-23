#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  {
    // we test writing on vtk and Matlab formats for PLANE with vectors, tensors, scalars
    int nx = 120, ny = 90;
    double xmin = -2.0, xmax = 2.0, ymin = -2.0, ymax = 2.0;
    GridInterpolationFull<Dimension2> grid;  
    grid.SetPlaneOutput(xmin, xmax, ymin, ymax, nx, ny);
    
    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    VectComplex_wp valc(points2d.GetM());
    VectComplex_wp valc2(points2d.GetM());
    Complex_wp vloc;
    Vector<TinyVector<Real_wp, 3> > val_vec(points2d.GetM());
    Vector<TinyVector<Complex_wp, 3> > valc_vec(points2d.GetM());
    Vector<TinyMatrix<Real_wp, Symmetric, 3, 3> > val_tens(points2d.GetM());
    Vector<TinyMatrix<Complex_wp, Symmetric, 3, 3> > valc_tens(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp x = points2d(i)(0);
	Real_wp y = points2d(i)(1);
	val(i) = sin(2.0*pi_wp*x)*cos(pi_wp*y);
	valc(i) = exp(-2.0*(x*x+y*y))*exp(-2.0*Iwp*pi_wp*x);
	valc2(i) = exp(-1.0*(x*x+y*y))*exp(-3.0*Iwp*pi_wp*x);
	vloc = cos(2.0*pi_wp*x*x)*Complex_wp(1.0+x, 0.5-x*y);
	val_vec(i).Init(val(i), real(valc(i)), real(vloc));
	valc_vec(i).Init(valc(i), valc2(i), vloc);
	val_tens(i)(0, 0) = val(i); val_tens(i)(0, 1) = real(valc(i)); val_tens(i)(0, 2) = imag(valc(i));
	val_tens(i)(1, 1) = real(valc2(i)); val_tens(i)(1, 2) = imag(valc2(i)); val_tens(i)(2, 2) = imag(vloc);
	
	valc_tens(i)(0, 0) = valc(i); valc_tens(i)(1, 1) = valc2(i); valc_tens(i)(2, 2) = vloc;
      }
    
    int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
    bool ascii = false;
    WriteVtk(val, grid, "scal_real.vtk", double_prec, ascii);
    WriteVtk(valc, grid, "scal_complex.vtk", double_prec, ascii);
    
    ofstream file_out("two_scal.vtk");
    // writing header
    grid.WriteVtk(file_out, double_prec, ascii);
    
    // then datas
    WriteVtk(val, "u", file_out, double_prec, ascii);
    WriteVtk(valc, "v", file_out, double_prec, ascii);
    //  WriteVtk(valv, "w", file_out, double_prec, ascii);
    file_out.close();
    
    WriteVtk(val_vec, grid, "vec.vtk", double_prec, ascii);
    WriteVtk(valc_vec, grid, "vec_cplx.vtk", double_prec, ascii);
    
    WriteVtk(val_tens, grid, "tens.vtk", double_prec, ascii);
    WriteVtk(valc_tens, grid, "tens_cplx.vtk", double_prec, ascii);
    
    WriteMatlab(val, grid, "scal_real.dat", double_prec, ascii);
    WriteMatlab(valc, grid, "scal_cplx.dat", double_prec, ascii);  
  }

  {
    // we test reading of vtk and Matlab files
    bool ascii = false;
    int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
    VectReal_wp val;     VectComplex_wp valc;
    GridInterpolationFull<Dimension2> grid;
    
    ReadMatlab(val, grid, "scal_real.dat", ascii);
    WriteVtk(val, grid, "Rscal_real.vtk", double_prec, ascii);

    ReadMatlab(valc, grid, "scal_cplx.dat", ascii);    
    WriteVtk(valc, grid, "Rscal_cplx.vtk", double_prec, ascii);
    
    val.Clear(); valc.Clear();
    ReadVtk(val, grid, "scal_real.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "Rscal_real.dat", double_prec, ascii);
      
    ReadVtk(valc, grid, "scal_complex.vtk", double_prec, ascii);
    WriteMatlab(valc, grid, "Rscal_cplx.dat", double_prec, ascii);

    VectR3 val_vec;
    Vector<TinyVector<Complex_wp, 3> > valc_vec;
    ReadVtk(val_vec, grid, "vec.vtk", double_prec, ascii);
    for (int i = 0; i < val_vec.GetM(); i++)
      val(i) = val_vec(i)(0);
    
    WriteMatlab(val, grid, "Rvec_real0.dat", double_prec, ascii);

    for (int i = 0; i < val_vec.GetM(); i++)
      val(i) = val_vec(i)(1);
    
    WriteMatlab(val, grid, "Rvec_real1.dat", double_prec, ascii);

    for (int i = 0; i < val_vec.GetM(); i++)
      val(i) = val_vec(i)(2);
    
    WriteMatlab(val, grid, "Rvec_real2.dat", double_prec, ascii);
    
    ReadVtk(valc_vec, grid, "vec_cplx.vtk", double_prec, ascii);
    for (int i = 0; i < valc_vec.GetM(); i++)
      valc(i) = valc_vec(i)(0);
    
    WriteMatlab(valc, grid, "Rvec_cplx0.dat", double_prec, ascii);

    for (int i = 0; i < valc_vec.GetM(); i++)
      valc(i) = valc_vec(i)(1);
    
    WriteMatlab(valc, grid, "Rvec_cplx1.dat", double_prec, ascii);

    for (int i = 0; i < valc_vec.GetM(); i++)
      valc(i) = valc_vec(i)(2);
    
    WriteMatlab(valc, grid, "Rvec_cplx2.dat", double_prec, ascii);

    Vector<TinyMatrix<Real_wp, Symmetric, 3, 3> > val_tens;
    Vector<TinyMatrix<Complex_wp, Symmetric, 3, 3> > valc_tens;
    ReadVtk(val_tens, grid, "tens.vtk", double_prec, ascii);
    ReadVtk(valc_tens, grid, "tens_cplx.vtk", double_prec, ascii);
    
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++)
	{
	  for (int k = 0; k < val_tens.GetM(); k++)
	    {
	      val(k) = val_tens(k)(i, j);
	      valc(k) = valc_tens(k)(i, j);
	    }
	  
	  string entier = to_str(i) + to_str(j);
	  WriteMatlab(val, grid, "Rtens_real"+entier+".dat", double_prec, ascii);
	  WriteMatlab(valc, grid, "Rtens_cplx"+entier+".dat", double_prec, ascii);
	}
  }

  {
    int nx = 120, ny = 100;
    R3 ptO, ptA(1, 0.2, 0.1), ptB(-0.2, 1.0, 0.1);
    
    // testing PLANE_AXI
    GridInterpolationFull<Dimension2> grid;  
    grid.SetPlaneAxiOutput(ptO, ptA, ptB, nx, ny);

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp r = points2d(i)(0);
	Real_wp z = points2d(i)(1);
	Real_wp x = r*cos(teta(i));
	Real_wp y = r*sin(teta(i));
	//DISP(r); DISP(x); DISP(y); DISP(z);
	val(i) = sin(10.0*pi_wp*x)*cos(3.0*pi_wp*y)*(z*z + 2.0*z + 0.2);
      }

    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "plane_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "plane_axi.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "plane_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rplane_axi.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "plane_axi.dat", ascii);
    WriteVtk(val, grid_r, "Rplane_axi.vtk", double_prec, ascii);
  }


  {
    int nx = 200;;
    R2 ptA(-2.0, -1.0), ptB(2.0, 1.0);
    
    // testing LINE
    GridInterpolationFull<Dimension2> grid;  
    grid.SetLineOutput(ptA(0), ptB(0), ptA(1), ptB(1), nx);

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp x = points2d(i)(0);
	Real_wp y = points2d(i)(1);
	val(i) = sin(5.0*pi_wp*x)*cos(3.0*pi_wp*y);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "line.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "line.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "line.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rline.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "line.dat", ascii);
    WriteVtk(val, grid_r, "Rline.vtk", double_prec, ascii);
  }

  {
    int nx = 200;;
    R3 ptA(-2.0, -1.0, 1.5), ptB(2.0, 1.0, 3.5);
    
    // testing LINE_AXI
    GridInterpolationFull<Dimension2> grid;  
    grid.SetLineAxiOutput(ptA, ptB, nx);

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp r = points2d(i)(0);
	Real_wp z = points2d(i)(1);
	Real_wp x = r*cos(teta(i));
	Real_wp y = r*sin(teta(i));
	val(i) = sin(5.0*pi_wp*x)*cos(3.0*pi_wp*y)*(-1.0+z);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "line_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "line_axi.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "line_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rline_axi.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "line_axi.dat", ascii);
    WriteVtk(val, grid_r, "Rline_axi.vtk", double_prec, ascii);
  }


  {
    R3 ptA(-2.0, -1.0, 1.5);
    
    // testing POINT_AXI
    GridInterpolationFull<Dimension2> grid;  
    grid.SetPointAxiOutput(ptA);

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp r = points2d(i)(0);
	Real_wp z = points2d(i)(1);
	Real_wp x = r*cos(teta(i));
	Real_wp y = r*sin(teta(i));
	val(i) = sin(5.0*pi_wp*x)*cos(3.0*pi_wp*y)*(-1.0+z);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "point_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "point_axi.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "point_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rpoint_axi.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "point_axi.dat", ascii);
    WriteVtk(val, grid_r, "Rpoint_axi.vtk", double_prec, ascii);
  }

  {
    R2 ptA(-2.0, -1.0);
    
    // testing POINT
    GridInterpolationFull<Dimension2> grid;  
    grid.SetPointOutput(ptA);

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp x = points2d(i)(0);
	Real_wp y = points2d(i)(1);
	val(i) = sin(5.0*pi_wp*x)*cos(3.0*pi_wp*y);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "point.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "point.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "point.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rpoint.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "point.dat", ascii);
    WriteVtk(val, grid_r, "Rpoint.vtk", double_prec, ascii);
  }

  {
    int nx = 200;
    R2 center(-2.0, -1.0);
    Real_wp radius = 3.0;
    
    // testing CIRCLE
    GridInterpolationFull<Dimension2> grid;  
    grid.SetCircleOutput(center, radius, nx);

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp x = points2d(i)(0);
	Real_wp y = points2d(i)(1);
	val(i) = sin(4.5*pi_wp*x)*cos(2.8*pi_wp*y);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "circle.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "circle.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "circle.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rcircle.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "circle.dat", ascii);
    WriteVtk(val, grid_r, "Rcircle.vtk", double_prec, ascii);
  }

  {
    int nx = 200;
    R3 center(-2.0, -1.0, 1.5), normale(1.0, 2.0, 0.5);
    Real_wp radius_x = 2.0, radius_y = 3.0;
    
    // testing CIRCLE_AXI
    GridInterpolationFull<Dimension2> grid;  
    grid.SetCircleAxiOutput(center, normale, radius_x, radius_y, nx);

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp r = points2d(i)(0);
	Real_wp z = points2d(i)(1);
	Real_wp x = r*cos(teta(i));
	Real_wp y = r*sin(teta(i));
	val(i) = sin(4.2*pi_wp*x)*cos(2.2*pi_wp*y)*(-1.0+z+0.3*z*z);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "circle_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "circle_axi.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "circle_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rcircle_axi.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "circle_axi.dat", ascii);
    WriteVtk(val, grid_r, "Rcircle_axi.vtk", double_prec, ascii);
  }

  {
    int nx = 200;
    R2 center(-2.0, -1.0);
    Real_wp a = 3.0, b = 1.5;
    
    ofstream file_out("PointsFile.txt");
    file_out.precision(15);
    for (int i = 0; i < nx; i++)
      {
	Real_wp teta = Real_wp(i)*2.0*pi_wp/nx;
	Real_wp x = center(0) + a*cos(teta);
	Real_wp y = center(1) + b*sin(teta);
	file_out << x << " " << y << '\n';
      }
    
    file_out.close();

    // testing POINTS_FILE    
    GridInterpolationFull<Dimension2> grid;  
    grid.SetPointsFileOutput("PointsFile.txt");

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp x = points2d(i)(0);
	Real_wp y = points2d(i)(1);
	val(i) = sin(4.5*pi_wp*x)*cos(2.8*pi_wp*y);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "points_file.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "points_file.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "points_file.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rpoints_file.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "points_file.dat", ascii);
    WriteVtk(val, grid_r, "Rpoints_file.vtk", double_prec, ascii);
  }

  {
    int nx = 100, ny = 100;
    R3 center(-2.0, -1.0, 1.5);
    Real_wp a = 2.0, b = 1.5, c = 2.5;

    ofstream file_out("PointsFile.txt");
    file_out.precision(15);
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
	{
	  Real_wp teta = Real_wp(i)*pi_wp/nx;
	  Real_wp phi = Real_wp(j)*2.0*pi_wp/ny;
	  Real_wp x = center(0) + a*sin(teta)*cos(phi);
	  Real_wp y = center(1) + b*sin(teta)*sin(phi);
	  Real_wp z = center(2) + c*cos(teta);
	  file_out << x << " " << y << ' ' << z << '\n';
	}
    
    file_out.close();
    
    // testing POINTS_FILE_AXI
    GridInterpolationFull<Dimension2> grid;  
    grid.SetPointsFileAxiOutput("PointsFile.txt");

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp r = points2d(i)(0);
	Real_wp z = points2d(i)(1);
	Real_wp x = r*cos(teta(i));
	Real_wp y = r*sin(teta(i));
	val(i) = sin(4.2*pi_wp*x)*cos(2.2*pi_wp*y)*(-1.0+z+0.3*z*z);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "points_file_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "points_file_axi.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "points_file_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rpoints_file_axi.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "points_file_axi.dat", ascii);
    WriteVtk(val, grid_r, "Rpoints_file_axi.vtk", double_prec, ascii);
  }

  {
    int nx = 100, ny = 120, nz = 110;
    R3 ptO(0.23, 0.51, 0.4), ptMin(-1, -1, -1), ptMax(1.0, 1.0, 1.0);
    
    // testing THREE_PLANES_AXI
    GridInterpolationFull<Dimension2> grid;  
    grid.SetThreePlanesAxiOutput(ptO, ptMin, ptMax, nx, ny, nz);

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp r = points2d(i)(0);
	Real_wp z = points2d(i)(1);
	Real_wp x = r*cos(teta(i));
	Real_wp y = r*sin(teta(i));
	//DISP(r); DISP(x); DISP(y); DISP(z);
	val(i) = sin(9.23*pi_wp*x)*cos(2.8*pi_wp*y)*cos(3.3*pi_wp*z);
      }

    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "three_planes_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "three_planes_axi.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "three_planes_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rthree_planes_axi.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "three_planes_axi.dat", ascii);
    WriteVtk(val, grid_r, "Rthree_planes_axi.vtk", double_prec, ascii);
  }

  {
    int nx = 100, ny = 120, nz = 110;
    R3 ptMin(-1, -1, -1), ptMax(1.0, 1.0, 1.0);
    
    // testing VOLUME_AXI
    GridInterpolationFull<Dimension2> grid;  
    grid.SetVolumeAxiOutput(ptMin, ptMax, nx, ny, nz);

    VectReal_wp teta; VectR2 points2d;
    grid.GenerateGridPoints(points2d, teta);
    
    VectReal_wp val(points2d.GetM());
    for (int i = 0; i < points2d.GetM(); i++)
      {
	Real_wp r = points2d(i)(0);
	Real_wp z = points2d(i)(1);
	Real_wp x = r*cos(teta(i));
	Real_wp y = r*sin(teta(i));
	//DISP(r); DISP(x); DISP(y); DISP(z);
	val(i) = sin(9.23*pi_wp*x)*cos(2.8*pi_wp*y)*cos(3.3*pi_wp*z);
      }

    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "volume_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "volume_axi.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension2> grid_r;
    ReadVtk(val, grid_r, "volume_axi.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rvolume_axi.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "volume_axi.dat", ascii);
    WriteVtk(val, grid_r, "Rvolume_axi.vtk", double_prec, ascii);
  }

  /* 3-D interpolation grid */

  {
    int nx = 120, ny = 100;
    R3 ptO, ptA(1, 0.2, 0.1), ptB(-0.2, 1.0, 0.1);
    
    // testing PLANE for 3-D 
    GridInterpolationFull<Dimension3> grid;  
    grid.SetPlaneOutput(ptO, ptA, ptB, nx, ny);

    VectR3 points3d;
    grid.GenerateGridPoints(points3d);
    
    VectReal_wp val(points3d.GetM());
    for (int i = 0; i < points3d.GetM(); i++)
      {
	Real_wp x = points3d(i)(0);
	Real_wp y = points3d(i)(1);
	Real_wp z = points3d(i)(2);
	val(i) = sin(10.0*pi_wp*x)*cos(3.0*pi_wp*y)*(z*z + 2.0*z + 0.2);
      }

    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "plane_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "plane_3d.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension3> grid_r;
    ReadVtk(val, grid_r, "plane_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rplane_3d.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "plane_3d.dat", ascii);
    WriteVtk(val, grid_r, "Rplane_3d.vtk", double_prec, ascii);
  }

  {
    int nx = 200;;
    R3 ptA(-2.0, -1.0, 1.5), ptB(2.0, 1.0, 3.5);
    
    // testing LINE
    GridInterpolationFull<Dimension3> grid;  
    grid.SetLineOutput(ptA, ptB, nx);

    VectR3 points3d;
    grid.GenerateGridPoints(points3d);
    
    VectReal_wp val(points3d.GetM());
    for (int i = 0; i < points3d.GetM(); i++)
      {
	Real_wp x = points3d(i)(0);
	Real_wp y = points3d(i)(1);
	Real_wp z = points3d(i)(2);
	val(i) = sin(5.0*pi_wp*x)*cos(3.0*pi_wp*y)*(-1.0+z);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "line_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "line_3d.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension3> grid_r;
    ReadVtk(val, grid_r, "line_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rline_3d.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "line_3d.dat", ascii);
    WriteVtk(val, grid_r, "Rline_3d.vtk", double_prec, ascii);
  }


  {
    R3 ptA(-2.0, -1.0, 1.5);
    
    // testing POINT
    GridInterpolationFull<Dimension3> grid;  
    grid.SetPointOutput(ptA);

    VectR3 points3d;
    grid.GenerateGridPoints(points3d);
    
    VectReal_wp val(points3d.GetM());
    for (int i = 0; i < points3d.GetM(); i++)
      {
	Real_wp x = points3d(i)(0);
	Real_wp y = points3d(i)(1);
	Real_wp z = points3d(i)(2);
	val(i) = sin(5.0*pi_wp*x)*cos(3.0*pi_wp*y)*(-1.0+z);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "point_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "point_3d.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension3> grid_r;
    ReadVtk(val, grid_r, "point_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rpoint_3d.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "point_3d.dat", ascii);
    WriteVtk(val, grid_r, "Rpoint_3d.vtk", double_prec, ascii);
  }

  {
    int nx = 200;
    R3 center(-2.0, -1.0, 1.5), normale(1.0, 2.0, 0.5);
    Real_wp radius_x = 2.0, radius_y = 3.0;
    
    // testing CIRCLE_AXI
    GridInterpolationFull<Dimension3> grid;  
    grid.SetCircleOutput(center, normale, radius_x, radius_y, nx);

    VectR3 points3d;
    grid.GenerateGridPoints(points3d);
    
    VectReal_wp val(points3d.GetM());
    for (int i = 0; i < points3d.GetM(); i++)
      {
	Real_wp x = points3d(i)(0);
	Real_wp y = points3d(i)(1);
	Real_wp z = points3d(i)(2);
	val(i) = sin(4.2*pi_wp*x)*cos(2.2*pi_wp*y)*(-1.0+z+0.3*z*z);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "circle_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "circle_3d.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension3> grid_r;
    ReadVtk(val, grid_r, "circle_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rcircle_3d.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "circle_3d.dat", ascii);
    WriteVtk(val, grid_r, "Rcircle_3d.vtk", double_prec, ascii);
  }

  {
    int nx = 100, ny = 100;
    R3 center(-2.0, -1.0, 1.5);
    Real_wp a = 2.0, b = 1.5, c = 2.5;

    ofstream file_out("PointsFile.txt");
    file_out.precision(15);
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
	{
	  Real_wp teta = Real_wp(i)*pi_wp/nx;
	  Real_wp phi = Real_wp(j)*2.0*pi_wp/ny;
	  Real_wp x = center(0) + a*sin(teta)*cos(phi);
	  Real_wp y = center(1) + b*sin(teta)*sin(phi);
	  Real_wp z = center(2) + c*cos(teta);
	  file_out << x << " " << y << ' ' << z << '\n';
	}
    
    file_out.close();
    
    // testing POINTS_FILE
    GridInterpolationFull<Dimension3> grid;  
    grid.SetPointsFileOutput("PointsFile.txt");

    VectR3 points3d;
    grid.GenerateGridPoints(points3d);
    
    VectReal_wp val(points3d.GetM());
    for (int i = 0; i < points3d.GetM(); i++)
      {
	Real_wp x = points3d(i)(0);
	Real_wp y = points3d(i)(1);
	Real_wp z = points3d(i)(2);
	val(i) = sin(4.2*pi_wp*x)*cos(2.2*pi_wp*y)*(-1.0+z+0.3*z*z);
      }
    
    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "points_file_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "points_file_3d.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension3> grid_r;
    ReadVtk(val, grid_r, "points_file_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rpoints_file_3d.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "points_file_3d.dat", ascii);
    WriteVtk(val, grid_r, "Rpoints_file_3d.vtk", double_prec, ascii);
  }

  {
    int nx = 100, ny = 120, nz = 110;
    R3 ptO(0.23, 0.51, 0.4), ptMin(-1, -1, -1), ptMax(1.0, 1.0, 1.0);
    
    // testing THREE_PLANES
    GridInterpolationFull<Dimension3> grid;  
    grid.SetThreePlanesOutput(ptO, ptMin, ptMax, nx, ny, nz);

    VectR3 points3d;
    grid.GenerateGridPoints(points3d);
    
    VectReal_wp val(points3d.GetM());
    for (int i = 0; i < points3d.GetM(); i++)
      {
	Real_wp x = points3d(i)(0);
	Real_wp y = points3d(i)(1);
	Real_wp z = points3d(i)(2);
	val(i) = sin(9.23*pi_wp*x)*cos(2.8*pi_wp*y)*cos(3.3*pi_wp*z);
      }

    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "three_planes_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "three_planes_3d.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension3> grid_r;
    ReadVtk(val, grid_r, "three_planes_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rthree_planes_3d.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "three_planes_3d.dat", ascii);
    WriteVtk(val, grid_r, "Rthree_planes_3d.vtk", double_prec, ascii);
  }

  {
    int nx = 100, ny = 120, nz = 110;
    R3 ptMin(-1, -1, -1), ptMax(1.0, 1.0, 1.0);
    
    // testing VOLUME
    GridInterpolationFull<Dimension3> grid;  
    grid.SetVolumeOutput(ptMin, ptMax, nx, ny, nz);

    VectR3 points3d;
    grid.GenerateGridPoints(points3d);
    
    VectReal_wp val(points3d.GetM());
    for (int i = 0; i < points3d.GetM(); i++)
      {
	Real_wp x = points3d(i)(0);
	Real_wp y = points3d(i)(1);
	Real_wp z = points3d(i)(2);
	val(i) = sin(9.23*pi_wp*x)*cos(2.8*pi_wp*y)*cos(3.3*pi_wp*z);
      }

    int double_prec = OutputTypeEnum::SINGLE_PRECISION;
    bool ascii = false;    
    WriteVtk(val, grid, "volume_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid, "volume_3d.dat", double_prec, ascii);
    
    val.Clear();
    GridInterpolationFull<Dimension3> grid_r;
    ReadVtk(val, grid_r, "volume_3d.vtk", double_prec, ascii);
    WriteMatlab(val, grid_r, "Rvolume_3d.dat", double_prec, ascii);
    
    val.Clear();
    ReadMatlab(val, grid_r, "volume_3d.dat", ascii);
    WriteVtk(val, grid_r, "Rvolume_3d.vtk", double_prec, ascii);
  }
  
}
