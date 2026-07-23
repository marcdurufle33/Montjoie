#ifndef MONTJOIE_FILE_OUTPUT_OPEN_CV_HXX

#include "opencv.hpp"

namespace Montjoie
{
  class ColorMapEnum
  {
  public :
    enum {JET, HOT};
    
    static void GetColorMap(int cmap, Matrix<int>& coef, int& r, int&g, int&b);
    
  };

  // lecture d'une image noir et bleu dane une matrice A
  // on stocke la matrice dans une vecteur avec v(i*n+j) = A_(i,j)
  void ReadPicture(const string& file_name, Vector<int>& v, int& m, int& n);

  // ecriture d'une image noir et blanc
  void WritePicture(const string& file_name, const Vector<int>& v, int m, int n);
  
  // lecture d'une image couleur
  // trois matrices pour le vert rouge et bleu
  void ReadColorPicture(const string& file_name, Vector<int>& vred,
			Vector<int>& vgreen, Vector<int>& vblue, int& m, int& n);

  // ecriture d'une image couleur
  void WriteColorPicture(const string& file_name, const Vector<int>& vred,
			 const Vector<int>& vgreen, const Vector<int>& vblue, int m, int n);
  
  void WriteJpeg(VectReal_wp& val, const string& file_output,
		 int Nx, int m0, int n0, int i0, int j0, const Vector<int>& vgray,
		 const Real_wp& cmin, const Real_wp& cmax, int cmap);

  void WriteEdge(int x1, int y1, int x2, int y2, int m, int n, VectReal_wp& Val);
  
  void WriteContour(const Mesh<Dimension2>& mesh, const Real_wp& xmin, const Real_wp& xmax,
		    const Real_wp& ymin, const Real_wp& ymax, int m, int n, VectReal_wp& Val);

}

#define MONTJOIE_FILE_OUTPUT_OPEN_CV_HXX
#endif
