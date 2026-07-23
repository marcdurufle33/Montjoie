#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

template<class T>
void GetRandNumber(T & x)
{
  x = T(rand())/RAND_MAX;
}

void GetRandNumber(int & x)
{
  x = rand()%1000;
}

template<class T>
void GetRandNumber(complex<T> & x)
{
  x = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}

class MyProblem : public InputDataProblem_Base
{
public:
  void SetInputData(const string& keyword, const VectString& param)
  {
    DISP(keyword); DISP(param);
    if (param.GetM() > 0)
      {
	DISP(param(0)+"test.mesh");
      }
  }
  
};

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  Real_wp xmin(-5), xmax(4); int Nx = 121;
  Real_wp ymin(-6), ymax(3); int Ny = 231;
  Real_wp zmin(-4.5), zmax(5.5); int Nz = 107;
  Real_wp ra(2.4), rb(0.8);
  
  VectReal_wp x, y, z;
  Linspace(xmin, xmax, Nx, x);
  Linspace(ymin, ymax, Ny, y);
  Linspace(zmin, zmax, Nz, z);

  // testing GetParallelDistributionPoints
  int nbPoints = 53; int nloc, num0;
  int sumP = 0;
  for (int rank = 0; rank < 5; rank++)
    {
      GetParallelDistributionPoints(5, rank, nbPoints, nloc, num0, 1);
      if (num0 != sumP)
	{
	  cout << "GetParallelDistributionPoints incorrect" << endl;
	  DISP(rank); DISP(nloc); DISP(num0);
	  abort();
	}
      
      sumP += nloc;
    }

  if (sumP != nbPoints)
    {
      cout << "GetParallelDistributionPoints incorrect" << endl;
      DISP(sumP); DISP(nbPoints);
      abort();      
    }
  
  // testing 1-D outputs
  {
    GridInterpolationFull<Dimension1> grid;
    grid.Init(grid.LINE, xmin, xmax, Nx);
    
    VectReal_wp val(Nx);
    for (int i = 0; i < Nx; i++)
      val(i) = sin(x(i));
    
    val.WriteText("output1D_ascii.dat");
    WriteMatlab(val, grid, "output1D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);
  }
  
  // testing 2-D outputs
  {
    // output on a plane
    GridInterpolationFull<Dimension2> grid;
    grid.SetPlaneOutput(xmin, xmax, ymin, ymax, Nx, Ny);
    
    Matrix<Real_wp> u(Nx, Ny); VectReal_wp v(Nx*Ny);
    for (int i = 0; i < Nx; i++)
      for (int j = 0; j < Ny; j++)
	{
	  u(i, j) = sin(x(i))*cos(y(j));
	  v(j*Nx + i) = u(i, j);
	}

    // here loadND loads the transpose
    u.WriteText("plane2D_ascii.dat");
    WriteMatlab(v, grid, "plane2D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);

    // output on a line
    grid.SetLineOutput(xmin, xmax, ymin, ymax, Nx);

    v.Reallocate(Nx);
    for (int i = 0; i < Nx; i++)
      v(i) = sin(x(i))*cos(y(i));	

    v.WriteText("line2D_ascii.dat");
    WriteMatlab(v, grid, "line2D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);

    // output on a point
    grid.SetPointOutput(R2(xmin, ymin));

    v.Reallocate(1);
    v(0) = sin(xmin)*cos(ymin);

    v.WriteText("point2D_ascii.dat");
    WriteMatlab(v, grid, "point2D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);

    // output on a circle
    grid.SetCircleOutput(R2(xmin, ymin), ra, Nx);

    v.Reallocate(Nx);
    for (int i = 0; i < Nx; i++)
      {
	Real_wp teta = 2*i*pi_wp/Nx;
	Real_wp x_ = xmin + ra*cos(teta);
	Real_wp y_ = ymin + ra*sin(teta);
	v(i) = sin(x_)*cos(y_);
      }

    v.WriteText("circle2D_ascii.dat");
    WriteMatlab(v, grid, "circle2D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);

    // output on a set of points
    grid.SetOutputType(grid.POINTS_FILE);

    grid.GlobalPoints2D.Reallocate(Nx);
    v.Reallocate(Nx);
    for (int i = 0; i < Nx; i++)
      {
	Real_wp x_ = 2*Real_wp(rand())/RAND_MAX - 1.0;
	Real_wp y_ = 2*Real_wp(rand())/RAND_MAX - 1.0;
	grid.GlobalPoints2D(i).Init(x_, y_);
	v(i) = sin(x_)*cos(y_);
      }
    
    v.WriteText("points2D_ascii.dat");
    WriteMatlab(v, grid, "points2D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);
  }

  // testing 3-D outputs
  {
    // output on a 3-D grid
    GridInterpolationFull<Dimension3> grid;
    grid.SetVolumeOutput(R3(xmin, ymin, zmin), R3(xmax, ymax, zmax), Nx, Ny, Nz);
    
    Array3D<Real_wp> u(Nx, Ny, Nz); VectReal_wp v(Nx*Ny*Nz);
    for (int i = 0; i < Nx; i++)
      for (int j = 0; j < Ny; j++)
	for (int k = 0; k < Nz; k++)
	  {
	    u(i, j, k) = sin(x(i))*cos(y(j))*exp(z(k)/3);
	    v(k*Nx*Ny + j*Nx + i) = u(i, j, k);
	  }
    
    // here loadND loads the transpose
    // you need to do [X, Y, Z, coor, V] = loadND('volume3D.dat')
    // Vfinal = transpose(V,[2,1,0]) to have the same result as loadArray3D_real
    u.Write("volume3D_binary.dat");
    WriteMatlab(v, grid, "volume3D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);

    // output on three planes
    R3 center(0.4, 0.8, 0.2);
    grid.SetThreePlanesOutput(center, R3(xmin, ymin, zmin), R3(xmax, ymax, zmax),
			      Nx, Ny, Nz);

    Matrix<Real_wp> u1(Ny, Nz), u2(Nx, Nz), u3(Nx, Ny);
    v.Reallocate(Nx*Ny + Nx*Nz + Ny*Nz);
    int offset = 0;
    for (int i = 0; i < Ny; i++)
      for (int j = 0; j < Nz; j++)
	{
	  Real_wp x_(center(0)), y_(y(i)), z_(z(j));
	  u1(i, j) = sin(x_)*cos(y_)*exp(z_/3);
	  v(offset + i*Nz + j) = u1(i, j);
	}

    offset += Ny*Nz;
    for (int i = 0; i < Nx; i++)
      for (int j = 0; j < Nz; j++)
	{
	  Real_wp x_(x(i)), y_(center(1)), z_(z(j));
	  u2(i, j) = sin(x_)*cos(y_)*exp(z_/3);
	  v(offset + i*Nz + j) = u2(i, j);
	}

    offset += Nx*Nz;
    for (int i = 0; i < Nx; i++)
      for (int j = 0; j < Ny; j++)
	{
	  Real_wp x_(x(i)), y_(y(j)), z_(center(2));
	  u3(i, j) = sin(x_)*cos(y_)*exp(z_/3);
	  v(offset + i*Ny + j) = u3(i, j);
	}
    
    // the transpose is performed on each matrix when using loadND
    u1.Write("three_planes3D_u1.dat");
    u2.Write("three_planes3D_u2.dat");
    u3.Write("three_planes3D_u3.dat");
    WriteMatlab(v, grid, "three_planes3D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);

    // output on a plane
    R3 vec_u(2.4, -1.4, 0.4);
    R3 vec_v(0.7, -0.2, 2.8);
    v.Reallocate(Nx*Ny);
    for (int i = 0; i < Nx; i++)
      for (int j = 0; j < Ny; j++)
	{
	  R3 pt = center + Real_wp(i)/(Nx-1)*vec_u + Real_wp(j)/(Ny-1)*vec_v;
	  Real_wp x_(pt(0)), y_(pt(1)), z_(pt(2));
	  u3(i, j) = sin(x_)*cos(y_)*exp(z_/3);
	  v(j*Nx + i) = u3(i, j);
	}

    grid.SetPlaneOutput(center, center+vec_u, center+vec_v, Nx, Ny);
    
    // de nouveau transposee (visu 2-D)
    u3.WriteText("plane3D_ascii.dat");
    WriteMatlab(v, grid, "plane3D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);

    // output on a line
    v.Reallocate(Nx);
    for (int i = 0; i < Nx; i++)
      {
	R3 pt = center + Real_wp(i)/(Nx-1)*vec_u;
	Real_wp x_(pt(0)), y_(pt(1)), z_(pt(2));
	v(i) = sin(x_)*cos(y_)*exp(z_/3);
      }

    grid.SetLineOutput(center, center+vec_u, Nx);
    
    v.WriteText("line3D_ascii.dat");
    WriteMatlab(v, grid, "line3D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);

    // output on a circle
    for (int i = 0; i < Nx; i++)
      {
	Real_wp teta = Real_wp(i)/(Nx-1)*pi_wp*2;
	R3 pt = center + R3(ra*cos(teta), 0, rb*sin(teta));
	Real_wp x_(pt(0)), y_(pt(1)), z_(pt(2));
	v(i) = sin(x_)*cos(y_)*exp(z_/3);
      }
    
    grid.SetCircleOutput(center, R3(0, 1, 0), ra, rb, Nx);
    
    v.WriteText("circle3D_ascii.dat");
    WriteMatlab(v, grid, "circle3D.dat", OutputTypeEnum::DOUBLE_PRECISION, false);
    
  }

  // testing loadSismo
  {    
    VectReal_wp t, tsub;
    Linspace(Real_wp(0), Real_wp(3.2), Nx, tsub); t = tsub;
    Linspace(Real_wp(3.0), Real_wp(6), Ny, tsub); t.PushBack(tsub);
    Linspace(Real_wp(6), Real_wp(9.3), Nz, tsub); t.PushBack(tsub);
    Linspace(Real_wp(7.5), Real_wp(8.6), Ny, tsub); t.PushBack(tsub);

    VectReal_wp u, v;
    u.Reallocate(t.GetM());
    v.Reallocate(t.GetM());
    ofstream file_out("Sismo.dat"); file_out.precision(16);
    for (int i = 0; i < t.GetM(); i++)
      {
	u(i) = sin(t(i));
	v(i) = cos(t(i));
	file_out << t(i) << " " << u(i) << " " << v(i) << '\n';
      }

    file_out.close();
  }

  // testing load1D
  {
    VectComplex_wp u(Nx);
    GenerateRandomVector(u, Nx);

    u.Write("vecCplx_binary.dat");

    ofstream file_out("vecCplx_ascii.dat"); file_out.precision(16);
    for (int i = 0; i < u.GetM(); i++)
      file_out << real(u(i)) << ' ' << imag(u(i)) << '\n';
    
    file_out.close();
  }

  // testing load1D_real
  {
    VectReal_wp u(Nx);
    GenerateRandomVector(u, Nx);

    u.Write("vecReal_binary.dat");
    u.WriteText("vecReal_ascii.dat");    
  }

  // testing load1D_int
  {
    Vector<int> u(Nx);
    GenerateRandomVector(u, Nx);

    u.Write("vecInt_binary.dat");
    u.WriteText("vecInt_ascii.dat");    
  }
  
  return FinalizeMontjoie();
}
