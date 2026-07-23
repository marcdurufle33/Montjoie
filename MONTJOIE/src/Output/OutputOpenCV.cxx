#ifndef MONTJOIE_FILE_OUTPUT_OPEN_CV_CXX

#include "OutputOpenCV.hxx"

namespace Montjoie
{
  
  //! fills coefficients for predefined color maps
  /*!
    \param[in] cmap type of color map (JET, HOT, etc)
    \param[out] coef RGB components of level i are contained in coef(i, :)
    \param[out] r Red component for nan
    \param[out] g Green component for nan
    \param[out] b Blue component for nan
   */
  void ColorMapEnum::GetColorMap(int cmap, Matrix<int>& coef, int& r, int&g, int&b)
  {
    switch (cmap)
      {
      case JET :
        {
          coef.Reallocate(4*255, 3);
          r = 255; g = 255; b = 255;
          int nb = 0;
          for (int i = 0; i < 128; i++)
            {
              coef(nb, 0) = 0;
              coef(nb, 1) = 0;
              coef(nb, 2) = 128+i;
              nb++;
            }
          
          for (int i = 1; i < 256; i++)
            {
              coef(nb, 0) = 0;
              coef(nb, 1) = i;
              coef(nb, 2) = 255;
              nb++;
            }
          
          for (int i = 1; i < 256; i++)
            {
              coef(nb, 0) = i;
              coef(nb, 1) = 255;
              coef(nb, 2) = 255-i;
              nb++;
            }
          
          for (int i = 1; i < 256; i++)
            {
              coef(nb, 0) = 255;
              coef(nb, 1) = 255-i;
              coef(nb, 2) = 0;
              nb++;
            }
          
          for (int i = 1; i < 128; i++)
            {
              coef(nb, 0) = 255-i;
              coef(nb, 1) = 0;
              coef(nb, 2) = 0;
              nb++;
            }
        }
        break;
      case HOT:
        {
          coef.Reallocate(639, 3);
          r = 0; g = 0; b = 255;
          int nb = 0;
          for (int i = 0; i < 256; i++)
            {
              coef(nb, 0) = i;
              coef(nb, 1) = 0;
              coef(nb, 2) = 0;
              nb++;
            }
          
          for (int i = 1; i < 256; i++)
            {
              coef(nb, 0) = 255;
              coef(nb, 1) = i;
              coef(nb, 2) = 0;
              nb++;
            }
          
          for (int i = 0; i < 128; i++)
            {
              coef(nb, 0) = 255;
              coef(nb, 1) = 255;
              coef(nb, 2) = 2*i+1;
              nb++;
            }          
        }
        break;
      }
  }


  // Lecture d'une image en niveaux de gris
  void ReadPicture(const string& file_name, Vector<int>& v, int& m, int& n)
  {
    cv::Mat_<int> A;
    A = cv::imread(file_name, 0);
    
    cv::Size s = A.size();
    m = s.height;
    n = s.width;
    v.Resize(m*n);

    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
	v(i*n + j) = A(i, j);
    
  }
  
  
  // Ecriture d'une image en niveaux de gris
  void WritePicture(const string& file_name, const Vector<int>& v, int m, int n)
  {
    cv::Mat_<int> A(m, n);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
	A(i, j) = v(i*n + j);
    
    cv::imwrite(file_name, A);
    
  }
  
  
  // Lecture d'une image en couleur
  void ReadColorPicture(const string& file_name, Vector<int>& vred,
			Vector<int>& vgreen, Vector<int>& vblue, int& m, int& n)
  {
    IplImage* image = cvLoadImage(file_name.data());
    m = image->height;
    n = image->width;
    //int nc = image->nChannels;
    vred.Resize(m*n); vgreen.Resize(m*n); vblue.Resize(m*n);
    
    unsigned char *data= reinterpret_cast<unsigned char *>(image->imageData);
    
    int nb = 0;
    for (int i = 0;  i < m; i++)
      for (int j = 0; j < n; j++)
	{
	  vblue(i*n + j) = int(data[nb]);
	  vgreen(i*n + j) = int(data[nb+1]);
	  vred(i*n + j) = int(data[nb+2]);
	  
	  nb += 3;
	}
    
    cvReleaseImage(&image);
  }
  
  
  //! writes a vector in a picture by using opencv routines
  /*!
    \param[in] file_name output file name
    \param[in] vred Red components
    \param[in] vgreen Green components
    \param[in] vblue Blue components
    \param[in] m number of rows
    \param[in] n number of columns
    the user provides directly RGB components of a_{i, j}
    where a_{i, j} = v(i*n + j)
  */
  void WriteColorPicture(const string& file_name, const Vector<int>& vred,
			 const Vector<int>& vgreen, const Vector<int>& vblue, int m, int n)
  {
    IplImage* image = cvCreateImage(cvSize(n, m), IPL_DEPTH_8U, 3);
    
    unsigned char *data= reinterpret_cast<unsigned char *>(image->imageData);
    
    int nb = 0;
    for (int i = 0;  i < m; i++)
      for (int j = 0; j < n; j++)      
	{
	  data[nb] = (unsigned char)(vblue(i*n + j));
	  data[nb+1] = (unsigned char)(vgreen(i*n + j));
	  data[nb+2] = (unsigned char)(vred(i*n + j));
	  
	  nb += 3;
	}
    
    cvSaveImage(file_name.data(), image);
    
    cvReleaseImage(&image);
  }
  

  //! writes a vector in a picture by using OpenCV routines
  /*!
    \param[in] val vector to be written
    \param[in] file_output output file name
    \param[in] Nx number of rows associated with the vector val
    \param[in] m0 number of rows of the output picture
    \param[in] n0 number of columns of the output picture
    \param[in] i0 offset for the rows
    \param[in] j0 offset for the columns
    \param[in] vgray possible background
    \param[in] cmin values below cmin have the same minimal color
    \param[in] cmax values below cmax have the same minimal color
    \param[in] cmap type of color map
    the vector val is a matrix a_{i, j} where a_{i, j} = val(j*Nx + i)
    this matrix is written in a picture of size m0 x n0 
    in the range i+i0:i+i0+m  x  j+j0:j+j0+n
    [cmin, cmax] is the color range
   */
  void WriteJpeg(VectReal_wp& val, const string& file_output,
		 int Nx, int m0, int n0, int i0, int j0, const Vector<int>& vgray,
		 const double& cmin, const double& cmax, int cmap)
  {
    Matrix<int> rgb_map;
    int red_nan, green_nan, blue_nan;
    ColorMapEnum::GetColorMap(cmap, rgb_map, red_nan, green_nan, blue_nan);

    int m = Nx;
    if (m <= 0)
      return;
    
    int n = val.GetM() / m;
    Vector<int> vred(m0*n0), vgreen(m0*n0), vblue(m0*n0);
    for (int i = 0; i < vred.GetM(); i++)
      {
	vred(i) = red_nan;
	vgreen(i) = green_nan;
	vblue(i) = blue_nan;
      }

    bool gray_scale = false;
    if (vgray.GetM() > 0)
      gray_scale = true;

    int nb_levels = rgb_map.GetM();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
	{
	  int k = (m-1-i)*n + j;	  
	  int k2 = (i+i0)*n0 + (j+j0);
	  if (isnan(val(k)))
	    {
	      vred(k2) = red_nan;
	      vgreen(k2) = green_nan;
	      vblue(k2) = blue_nan;
	    }
	  else if (isinf(val(k)))
	    {
	      if (gray_scale)
		{
		  vred(k2) = vgray(k);
		  vgreen(k2) = vgray(k);
		  vblue(k2) = vgray(k);
		}
	      else
		{
		  vred(k2) = 0;
		  vgreen(k2) = 0;
		  vblue(k2) = 0;
		}
	    }
	  else
	    {
	      int lvl = int(round((nb_levels-1)*(val(k) - cmin) / (cmax - cmin)));
	      if (lvl < 0)
		lvl = 0;
	      
	      if (lvl >= nb_levels)
		lvl = nb_levels - 1;

	      vred(k2) = rgb_map(lvl, 0);
	      vgreen(k2) = rgb_map(lvl, 1);
	      vblue(k2) = rgb_map(lvl, 2);
	    }
	}
    
    WriteColorPicture(file_output, vred, vgreen, vblue, m0, n0);
  }


  //! writes an edge between point (x1, y1) and (x2, y2)
  void WriteEdge(int x1, int y1, int x2, int y2, int m, int n, VectReal_wp& Val)
  {
    int x = x1, y = y1;
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int s1 = sign(x2 - x1);
    int s2 = sign(y2 - y1);
    int swap = 0;
    if (dy > dx)
      {
	int temp = dx;
	dx = dy;
	dy = temp;
	swap = 1;
      }
    
    int D = 2*dy - dx;
    for (int i = 0; i < dx; i++)
      {
	Val(m*y + x) = INFINITY; 
	while (D >= 0)
	  {
	    D = D - 2*dx;
	    if (swap)
	      x += s1;
	    else
	      y += s2;
	  }
	
	D = D + 2*dy;
	if (swap)
	  y += s2;
	else
	  x += s1;
      }
  }
  

  //! writes a contour in the picture
  void WriteContour(const Mesh<Dimension2>& mesh, const Real_wp& xmin, const Real_wp& xmax,
		    const Real_wp& ymin, const Real_wp& ymax, int m, int n, VectReal_wp& Val)
  {
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int n1 = mesh.BoundaryRef(i).numVertex(0);
	int n2 = mesh.BoundaryRef(i).numVertex(1);
	Real_wp x1 = mesh.Vertex(n1)(0);
	Real_wp y1 = mesh.Vertex(n1)(1);
	Real_wp x2 = mesh.Vertex(n2)(0);
	Real_wp y2 = mesh.Vertex(n2)(1);
	
	int i1 = toInteger(round(m*(x1 - xmin)/(xmax - xmin)));
	int j1 = toInteger(round(n*(y1 - ymin)/(ymax - ymin)));
	int i2 = toInteger(round(m*(x2 - xmin)/(xmax - xmin)));
	int j2 = toInteger(round(n*(y2 - ymin)/(ymax - ymin)));
	
	i1 = max(0, i1); i1 = min(m-1, i1);
	i2 = max(0, i2); i2 = min(m-1, i2);
	j1 = max(0, j1); j1 = min(n-1, j1);
	j2 = max(0, j2); j2 = min(n-1, j2);
	
	WriteEdge(i1, j1, i2, j2, m, n, Val);
      }
  }
  
}

#define MONTJOIE_FILE_OUTPUT_OPEN_CV_CXX
#endif
