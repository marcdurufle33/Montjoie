#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Output/MontjoieOutput.hxx"

#include "Output/OutputOpenCV.cxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc < 5)
    {
      cout<<"Cette commande demande au moins quatre arguments"<<endl;
      cout<<"film2D toto _U0.jpg cmin cmax"<<endl;
      return -1;
    }
  
  string root_input(argv[1]), ext_input(argv[2]);
  Real_wp cmin = atof(argv[3]), cmax = atof(argv[4]);  
  string file_coef("ONE");

  string base_ext = GetBaseString(ext_input);
  string ext_img = GetExtension(ext_input);
  
  int n0_inst = -1, n1_inst = -1;
  int nb_pixel_x = 0, nb_pixel_y = 0;

  Real_wp xmin(0), xmax(1), ymin(0), ymax(1);
  
  enum {REAL_PART, IMAG_PART, MODULUS_CPLX, REAL_FIELD};
  int type_field = REAL_FIELD;
  bool nan_enabled = true;
  string file_background("NONE");
  for (int i = 1; i < argc; i++)
    {
      if (string(argv[i]) == "-real")
	type_field = REAL_PART;
      else if (string(argv[i]) == "-imag")
	type_field = IMAG_PART;
      else if (string(argv[i]) == "-abs")
	type_field = MODULUS_CPLX;
      else if (string(argv[i]) == "-a")
	{
	  if (argc <= i+2)
	    {
	      cout << "-a must be followed by n0 and n1" << endl;
	      cout << "argc = " << argc << " and i+2 = " << i+2 << endl;
	      abort();
	    }
	  
	  // cas de plusieurs instantanes
	  n0_inst = atoi(argv[i+1]);
	  n1_inst = atoi(argv[i+2]);
	}
      else if (string(argv[i]) == "-c")
	{
	  if (argc <= i+1)
	    {
	      cout << "-c must be followed by the file name" << endl;
	      cout << "argc = " << argc << " and i+1 = " << i+1 << endl;
	      abort();
	    }
	  
	  file_coef = argv[i+1];
	}
      else if (string(argv[i]) == "-blank")
	{
	  if (argc <= i+2)
	    {
	      cout << "-blank must be followed by x0 and y0" << endl;
	      cout << "argc = " << argc << " and i+2 = " << i+2 << endl;
	      abort();
	    }
	  
	  // pixels blancs en plus en x et y
	  nb_pixel_x = atoi(argv[i+1]);
	  nb_pixel_y = atoi(argv[i+2]);
	}
      else if (string(argv[i]) == "-background")
	{
	  if (argc <= i+1)
	    {
	      cout << "-background must be followed by the file name" << endl;
	      cout << "argc = " << argc << " and i+1 = " << i+1 << endl;
	      abort();
	    }

	  file_background = argv[i+1];
	}
      else if (string(argv[i]) == "-bb")
	{
	  if (argc <= i+4)
	    {
	      cout << "-bb must be followed by xmin, xmax, ymin, ymax" << endl;
	      cout << "argc = " << argc << " and i+4 = " << i+2 << endl;
	      abort();
	    }

	  xmin = Real_wp(atof(argv[i+1]));
	  xmax = Real_wp(atof(argv[i+2]));
	  ymin = Real_wp(atof(argv[i+3]));
	  ymax = Real_wp(atof(argv[i+4]));
	}
      else if (string(argv[i]) == "-disable-nan")
	nan_enabled = false;
    }

  string extension;

  bool write_contour = false;
  Matrix<Real_wp> coef;
  Mesh<Dimension2> mesh;
  if (file_coef != "ONE")
    {
      extension = GetExtension(file_coef);
      if (extension != "mesh")
	coef.Read(file_coef);
      else
	{
	  write_contour = true;
	  mesh.Read(file_coef);
	}
    }

  Vector<int> field_bg;
  int m_bg, n_bg;
  
  if (file_background != "NONE")
    {
      ReadPicture(file_background, field_bg, n_bg, m_bg);
      DISP(n_bg); DISP(m_bg);
    }
  
  // on lit le fichier matlab
  bool ascii = false;
  VectReal_wp val; VectComplex_wp val_cplx;
  GridInterpolationFull<Dimension2> grid;  
  
  string file_name = root_input + base_ext + ".dat";
  if (n0_inst >= 0)
    file_name = root_input + NumberToString(n1_inst-1) + base_ext + ".dat";
  
  DISP(file_name);
  
  if (type_field == REAL_FIELD)
    ReadMatlab(val, grid, file_name, ascii);
  else
    {
      ReadMatlab(val_cplx, grid, file_name, ascii);
      val.Reallocate(val_cplx.GetM());
      if (type_field == REAL_PART)
	for (int i = 0; i < val.GetM(); i++)
	  val(i) = realpart(val_cplx(i));
      else if (type_field == IMAG_PART)
	for (int i = 0; i < val.GetM(); i++)
	  val(i) = imagpart(val_cplx(i));
      else if (type_field == MODULUS_CPLX)
	for (int i = 0; i < val.GetM(); i++)
	  val(i) = abs(val_cplx(i));
    }
  
  int m = grid.GetNbPointsX();
  int n = grid.GetNbPointsY();

  int m1 = m + nb_pixel_x, n1 = n + nb_pixel_y;
  int m0 = m1; int n0 = n1; int i0 = 0; int j0 = 0;
  // on prend un multiple de 4 pour eviter un bug int / unsigned char
  int factor = 4;
  if (m1%factor != 0)
    {
      m0 = (m1/factor + 1)*factor;
      i0 = (m0-m)/2;
    }

  if (n1%factor != 0)
    {
      n0 = (n1/factor + 1)*factor;
      j0 = (n0 - n)/2;
    }

  // on place les nan si demande
  Vector<bool> place_nan(val.GetM());
  place_nan.Fill(false);
  if (nan_enabled)
    for (int k = 0; k < val.GetM(); k++)
      {
	if (val(k) == 0)
	  {
	    val(k) = NAN;
	    place_nan(k) = true;
	  }
	else
	  {
	    int i = k%m, j = k/m;
	    if (coef.GetM() > 0)
	      val(k) *= coef(i, j);
	  }
      }

  Vector<int> vgray;
  if (file_background != "NONE")
    vgray.Reallocate(val.GetM());
  
  // on ecrit les fichier jpeg
  string file_output;
  if (n0_inst >= 0)
    {
      for (int p = n0_inst; p < n1_inst; p++)
	{
	  file_name = root_input + NumberToString(p) + base_ext + ".dat";
	  val.Reallocate(n*m);
	  if (type_field == REAL_FIELD)
	    ReadMatlab(val, grid, file_name, ascii);
	  else
	    {
	      ReadMatlab(val_cplx, grid, file_name, ascii);
	      if (type_field == REAL_PART)
		for (int i = 0; i < val.GetM(); i++)
		  val(i) = realpart(val_cplx(i));
	      else if (type_field == IMAG_PART)
		for (int i = 0; i < val.GetM(); i++)
		  val(i) = imagpart(val_cplx(i));
	      else if (type_field == MODULUS_CPLX)
		for (int i = 0; i < val.GetM(); i++)
		  val(i) = abs(val_cplx(i));
	    }
	  
	  for (int k = 0; k < val.GetM(); k++)
	    {
	      if (place_nan(k))
		val(k) = NAN;
	      else
		{
		  int i = k%m, j = k/m;
		  if (coef.GetM() > 0)
		    val(k) *= coef(i, j);
		}
	    }

	  if (write_contour)
	    WriteContour(mesh, xmin, xmax, ymin, ymax, m, n, val);

	  if (file_background != "NONE")
	    {
	      Real_wp dx = (xmax - xmin) / (m-1);
	      Real_wp dy = (ymax - ymin) / (n-1);
	      for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		  {
		    Real_wp x = xmin + i*dx;
		    Real_wp y = ymin + j*dy;

		    int ix = toInteger(round(x));
		    int iy = toInteger(round(y));
		    if ((ix > 0) && (ix < m_bg) && (iy > 0) && (iy < n_bg))
		      if (field_bg((n_bg-iy)*m_bg + ix) <= 230)
			{
			  val(j*m + i) = INFINITY;
			  vgray(j*m + i) = max(0, field_bg((n_bg-iy)*m_bg + ix)-60);
			}
		  }
	    }
	  
	  file_output = root_input + NumberToString(p) + base_ext + "." + ext_img;
	  WriteJpeg(val, file_output, n, n0, m0,
		    i0, j0, vgray, cmin, cmax, ColorMapEnum::JET);

	}
    }
  else
    {
      file_output = root_input + base_ext + "." + ext_img;
      WriteJpeg(val, file_output, m, m0, n0,
		i0, j0, vgray, cmin, cmax, ColorMapEnum::JET);
    }
  
  return FinalizeMontjoie();
}  

