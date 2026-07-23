#include "utils.hpp"

using namespace std;

string Get_file_extension (string pouet) {

  return pouet.substr(pouet.find_last_of('.')+1);

}

string Get_file_root (string pouet) {

  string debut = pouet.substr(0,pouet.find_last_of('.'));
  return debut.substr(pouet.find_last_of('/')+1);

}

// Geometry class

geom::geom() {

  N = Index(0);
  h = Vertex(1.e0);
  O = Vertex(0.e0);

}

geom::geom(string filename) {

  string plop = Get_file_extension(filename);

  if (plop=="vtk")
    Read_VTK   (filename);
  else if (plop=="dat")
    Read_Montjoie(filename);
  else {
    cerr << "Could not open " << filename << " to get geometry. Abandon ship !" << endl;
    abort();
  }

}

void geom::Read_VTK (string filename) {

  ifstream f;
  f.open(filename.c_str());
  if (!f.is_open()) {
    cerr << "Could not open " << filename << " to get geometry. Abandon ship !" << endl;
    abort();
  }
  string machin;
  getline(f,machin);      // #vtk legacy blah blah
  getline(f,machin);      // comment line
  getline(f,machin);      // ASCII
  if (machin!="ASCII") {
    cerr << "Sorry, my programmer is too dumb or lazy to read binary files." << endl;
    abort();
  }
  getline(f,machin);      // DATASET structured grid
  if (machin!="DATASET STRUCTURED_GRID") {
    cerr << "Sorry, my programmer is too dumb or lazy to read anything else than structured grids." << endl;
    abort();
  }
  getline(f,machin,' ');  // DIMENSIONS
#if NDIM == 2
  f >> N(0) >> N(1);
#endif
#if NDIM == 3
  f >> N(0) >> N(1) >> N(2);
#endif
  getline(f,machin);      // "\n" or "1 \n" 
  getline(f,machin);      // POINTS n double/float
 
  double junk;
  Vertex3D O3D,coord_next;
#if NDIM == 2
  // !!!!!
  // PlaneAxi Outputs are in the xOz plane
    f >> O3D(0) >> O3D(1) >> O3D(2);
    f >> coord_next(0) >> coord_next(1) >> coord_next(2); 
    h(0) = norm2(coord_next-O3D);
    O(0) = O3D(0);
    O(1) = O3D(2);

    for (int i=0;i<N(0)-1;i++) getline(f,machin);
    f >> coord_next(0) >> coord_next(1) >> coord_next(2);
    h(1) = norm2(coord_next-O3D);
#endif
#if NDIM == 3
    f >> O(0) >> O(1) >> O(2);
    f >> coord_next(0) >> coord_next(1) >> coord_next(2); 
    h(0) = norm2(coord_next-O);
    for (int i=0;i<N(0)-1;i++) getline(f,machin);
    f >> coord_next(0) >> coord_next(1) >> coord_next(2);
    h(1) = norm2(coord_next-O);
    for (int i=0;i<N(0)*(N(1)-1);i++) getline(f,machin);
    f >> coord_next(0) >> coord_next(1) >> coord_next(2);
    h(2) = norm2(coord_next-O);
#endif
 
  f.close();

}

void geom::Read_Montjoie (string filename) {

  ifstream f;
  f.open(filename.c_str());
  if (!f.is_open()) {
    cerr << "Could not open " << filename << " to get geometry. Abandon ship !" << endl;
    abort();
  }

#if NDIM==2
  double xmax,ymax,junk; int junk2;
  f >> junk2 >> junk2 >> junk2 >> junk2; // 4 first lines
  f >> O(0) >> xmax >> O(2) >> junk >> O(1) >> ymax;
  f >> junk >> junk >> junk; // "center"
  f >> N(0) >> N(1);
  N(2) = 1;

  h(0) = (xmax-O(0))/(N(0)-1.e0);
  h(1) = (ymax-O(1))/(N(1)-1.e0);
  h(2) = 1.e0;
#endif
#if NDIM==3
  // A voir avec un exemple de Sismogrille
#endif

  f.close();

}

void geom::Display() {

  cout.precision(PRECISION);
  cout << "------------------------------------------" << endl;
  cout << "Geometry" << endl;
  cout << "Points  : " << N << endl;
  cout << "Spacing : " << h << endl;
  cout << "Origin  : " << O << endl;
  cout << "------------------------------------------" << endl;

}

// IO Fields

void Read_VTK (Field& U, const string filename, const geom& g) {

  ifstream f;
  f.open(filename.c_str());
  if (!f.is_open()) {
    cerr << "Could not open " << filename << " to read field data. Mayday !" << endl;
    abort();
  }

  string machin;
  while (machin.substr(0,machin.find_first_of(' ')) != "POINT_DATA") getline(f,machin);

  int Npoints_file = atoi((machin.substr(machin.find_first_of(' '))).c_str());
  int Npoints_geom = product(g.N);
  if (Npoints_geom != Npoints_file ) {
    cerr << "In file " << filename << ", the number of data points ("
         << Npoints_file<< ") does not fit geometry ("<< Npoints_geom << " points)." << endl
         << "Will abort if not enough points." << endl;
    if (Npoints_geom > Npoints_file) abort();
  }

  getline(f,machin); // SCALARS u0 double/float 1
  getline(f,machin); // LOOKUP_TABLE default


  // We need to transpose the data as the iterator read the other way vtk is written
  U.resize(g.N);
  #if NDIM==2
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f >> U(i,j);
    }}
  #endif
  #if NDIM==3
    for (int k=0;k<g.N(2);k++) {
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f >> U(i,j,k);
    }}}
  #endif

  f.close();

}

void Read_VTK (CField& U, const string filename, const geom& g) {

  ifstream f;
  f.open(filename.c_str());
  if (!f.is_open()) {
    cerr << "Could not open " << filename << " to read field data. Mayday !" << endl;
    abort();
  }

  string machin;
  while (machin.substr(0,machin.find_first_of(' ')) != "POINT_DATA") getline(f,machin);

  int Npoints_file = atoi((machin.substr(machin.find_first_of(' '))).c_str());
  int Npoints_geom = product(g.N);
  if (Npoints_geom != Npoints_file ) {
    cerr << "In file " << filename << ", the number of data points ("
         << Npoints_file<< ") does not fit geometry ("<< Npoints_geom << " points)." << endl
         << "Will abort if not enough points." << endl;
    if (Npoints_geom > Npoints_file) abort();
  }

  getline(f,machin); // SCALARS u0 double/float 1
  getline(f,machin); // LOOKUP_TABLE default


  // We need to transpose the data as the iterator read the other way vtk is written
  U.resize(g.N);
  #if NDIM==2
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f >> real(U(i,j));
    }}
  #endif
  #if NDIM==3
    for (int k=0;k<g.N(2);k++) {
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f >> real(U(i,j,k));
    }}}
  #endif

  // Imaginary part
  getline(f,machin); // \n
  getline(f,machin); // SCALARS u1 double/float 1
  getline(f,machin); // LOOKUP_TABLE default
  #if NDIM==2
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f >> imag(U(i,j));
    }}
  #endif
  #if NDIM==3
    for (int k=0;k<g.N(2);k++) {
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f >> imag(U(i,j,k));
    }}}
  #endif

  f.close();

}

void Write_VTK (const Field& U, const string filename, const string varname, const geom& g) {

  ofstream f;
  f.open(filename.c_str());
  if (!f.is_open()) {
    cerr << "Could not open " << filename << " to write field data. Help !" << endl;
    abort();
  }
  f.precision(PRECISION);

  f << "# vtk DataFile Version 3.0" << endl
    << "post-treatment of Montjoie output" << endl
    << "ASCII" << endl
    << "DATASET STRUCTURED_POINTS" << endl
#if NDIM==2
    << "DIMENSIONS " << g.N(0) << " " << g.N(1) << " 1" << endl
    << "ORIGIN " << g.O(0) << " " << g.O(1) << " 0" << endl
    << "SPACING " << g.h(0) << " " << g.h(1) << " 1" << endl
#endif
#if NDIM==3
    << "DIMENSIONS " << g.N(0) << " " << g.N(1) << " " << g.N(2) << endl
    << "ORIGIN " << g.O(0) << " " << g.O(1) << " " << g.O(2)  << endl
    << "SPACING " << g.h(0) << " " << g.h(1) << " " << g.h(2) << endl
#endif 
    << "POINT_DATA " << product(g.N) << endl
    << "SCALARS " << varname << " double 1" << endl
    << "LOOKUP_TABLE default" << endl;

  #if NDIM==2
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f << U(i,j) << endl;
    }}
  #endif
  #if NDIM==3
    for (int k=0;k<g.N(2);k++) {
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f << U(i,j,k) << endl;
    }}}
  #endif

  f.close();

}

void Write_VTK (const CField& U, const string filename, const string varname, const geom& g) {

  ofstream f;
  f.open(filename.c_str());
  if (!f.is_open()) {
    cerr << "Could not open " << filename << " to write field data. Help !" << endl;
    abort();
  }
  f.precision(PRECISION);

  f << "# vtk DataFile Version 3.0" << endl
    << "post-treatment of Montjoie output" << endl
    << "ASCII" << endl
    << "DATASET STRUCTURED_POINTS" << endl
    << "DIMENSIONS " << g.N(0) << " " << g.N(1)
    #if NDIM==2
    << " 1" << endl
    #endif
    #if NDIM==3
    << " " << g.N(2) << endl
    #endif 
    << "ORIGIN " << g.O(0) << " " << g.O(1)
    #if NDIM==2
    << " 0" << endl
    #endif
    #if NDIM==3
    << " " << g.O(2) << endl
    #endif 
    << "SPACING " << g.h(0) << " " << g.h(1)
    #if NDIM==2
    << " 1" << endl
    #endif
    #if NDIM==3
    << " " << g.h(2) << endl
    #endif 
    << "POINT_DATA " << product(g.N) << endl
    << "SCALARS " << varname << "_real double 1" << endl
    << "LOOKUP_TABLE default" << endl;

  #if NDIM==2
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f << real(U(i,j)) << endl;
    }}
  #endif
  #if NDIM==3
    for (int k=0;k<g.N(2);k++) {
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f << real(U(i,j,k)) << endl;
    }}}
  #endif

  f << "SCALARS " << varname << "_imag double 1" << endl
    << "LOOKUP_TABLE default" << endl; 
  #if NDIM==2
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f << imag(U(i,j)) << endl;
    }}
  #endif
  #if NDIM==3
    for (int k=0;k<g.N(2);k++) {
    for (int j=0;j<g.N(1);j++) {
    for (int i=0;i<g.N(0);i++) {
      f << imag(U(i,j,k)) << endl;
    }}}
  #endif
  f.close();

}

void Display (const Field& U,const string label) {

  cout << endl << "Field " << label << ": size " << U.extent() << endl;
    for (Field::const_iterator it=U.begin();it!=U.end();it++)
      cout << "Index " << it.position() << " value " << *it << endl;
  cout << endl;

}

void Display (const CField& U,const string label) {

  cout << endl << "Field " << label << ": size " << U.extent() << endl;
    for (CField::const_iterator it=U.begin();it!=U.end();it++)
      cout << "Index " << it.position() << " value " << *it << endl;
  cout << endl;

}


// Field computations. There should be an assert T double or complex
// Why the fuck Array<T,NDIM>::iterator does not work ?

double integral(const Field& U, const geom& g) {
  return sum(U)*product(g.h);
};

double L1_norm(const Field& U, const geom& g) {
  Field AU(U.shape()); AU = fabs(U);
  return integral(AU,g);
};

double L2_norm(const Field& U, const geom& g) {
  Field UU(U.shape()); UU = U*U;
  return sqrt(integral(UU,g));
};

Field dx(const Field& U, const geom& g) {
  return derive(U,g,0);
};

Field dy(const Field& U, const geom& g) {
  return derive(U,g,1);
};

Field dz(const Field& U, const geom& g) {
  return derive(U,g,2);
};

Field derive(const Field& U, const geom& g, const int& dir) {

  Field dU(U.shape());
  dU = 0.e0;

  for (Field::iterator it=dU.begin(); it!=dU.end(); ++it) {

    Index I   = it.position();
    Index Im2 = I; Im2(dir)-=2;
    Index Im1 = I; Im1(dir)--;
    Index Ip1 = I; Ip1(dir)++;
    Index Ip2 = I; Ip2(dir)+=2;

    int i = I(dir);
    int N = g.N(dir);

    if (i==0) { // O(h^2)
      dU(I) = -3.e0*U(I) + 4.e0*U(Ip1) - U(Ip2); 
    }
    if (i==1 || i==N-2) { // O(h^2)
      dU(I) = U(Ip1) - U(Im1);
    }
    if (i==N-1) {// sym of i==0
      dU(I) = 3.e0*U(I) - 4.e0*U(Im1) + U(Im2);
    }
    if (i>1 && i<N-2) { // O(h^4)
      dU(I) = (U(Im2)+8.e0*(U(Ip1)-U(Im1))-U(Ip2))/6.e0;
    }

    dU(I) /= 2.e0*g.h(dir);
  }
  return dU;

};

// The same for complex -_______________-

Complex integral(const CField& U, const geom& g) {
  return sum(U)*product(g.h);
};

double L2_norm(const CField& U, const geom& g) {
  Field  UU(U.shape()); UU = pow(abs(U),2); 
  return sqrt(integral(UU,g));
};

CField dx(const CField& U, const geom& g) {
  return derive(U,g,0);
};

CField dy(const CField& U, const geom& g) {
  return derive(U,g,1);
};

CField dz(const CField& U, const geom& g) {
  return derive(U,g,2);
};

CField derive(const CField& U, const geom& g, const int& dir) {

  CField dU(U.shape());
  dU = Complex(0.e0);

  for (CField::iterator it=dU.begin(); it!=dU.end(); ++it) {

    Index I   = it.position();
    Index Im2 = I; Im2(dir)-=2;
    Index Im1 = I; Im1(dir)--;
    Index Ip1 = I; Ip1(dir)++;
    Index Ip2 = I; Ip2(dir)+=2;

    int i = I(dir);
    int N = g.N(dir);

    if (i==0) { // O(h^2)
      dU(I) = -3.e0*U(I) + 4.e0*U(Ip1) - U(Ip2); 
    }
    if (i==1 || i==N-2) { // O(h^2)
      dU(I) = U(Ip1) - U(Im1);
    }
    if (i==N-1) { // sym of i==0
      dU(I) = 3.e0*U(I) - 4.e0*U(Im1) + U(Im2);
    }
    if (i>1 && i<N-2) { // inside O(h^4)
      dU(I) = (U(Im2)+8.e0*(U(Ip1)-U(Im1))-U(Ip2))/6.e0;
    }

    dU(I) /= 2.e0*g.h(dir);
  }
  return dU;

};

