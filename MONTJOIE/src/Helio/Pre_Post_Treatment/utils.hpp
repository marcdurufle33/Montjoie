#ifndef _UTILS_HXX_
#define _UTILS_HXX_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>

#include <complex>
#include <vector>
#include <string>
#include <blitz/array.h>
BZ_USING_NAMESPACE(blitz)

//#define NDIM 2
#define PRECISION 12

using namespace std;

// Typedefs

typedef complex<double> Complex;

typedef TinyVector<int    ,NDIM> Index;
typedef TinyVector<double ,NDIM> Vertex;
typedef TinyVector<Complex,NDIM> Certex;

typedef TinyVector<int    ,3> Index3D;
typedef TinyVector<double ,3> Vertex3D;
typedef TinyVector<Complex,3> Certex3D;

typedef Array<int    ,NDIM> IField;
typedef Array<double ,NDIM> Field;
typedef Array<Complex,NDIM> CField;
typedef Array<Index  ,NDIM> Index_Field;
typedef Array<Vertex ,NDIM> Vertex_Field;
typedef Array<Certex ,NDIM> Certex_Field;

// Geometry class

class geom {

  public :

    geom();
    geom(string);

    Index  N; // Npoints
    Vertex h; // grid spacing
    Vertex O; // origin

    void Read_VTK(string);
    void Read_Montjoie(string);

    void Display();

};

// IO 

string Get_file_extension (string pouet);
string Get_file_root      (string pouet);

void Read_VTK  ( Field&,const string,const geom&);
void Read_VTK  (CField&,const string,const geom&);
void Write_VTK (const  Field&,const string,const string,const geom&);
void Write_VTK (const CField&,const string,const string,const geom&);

void Display (const  Field&,const string="");
void Display (const CField&,const string="");

// Field computations.
// Thanks Blitz++ or my stupidity to not allow me to template those functions...

double integral(const Field&,const geom& g);
double L1_norm (const Field&,const geom& g);
double L2_norm (const Field&,const geom& g);
Field  dx(const Field&,const geom& g);
Field  dy(const Field&,const geom& g);
Field  dz(const Field&,const geom& g);
Field derive(const Field& U, const geom& g, const int& dir);

Complex integral(const CField&,const geom& g);
double  L2_norm (const CField&,const geom& g);
CField dx(const CField&,const geom& g);
CField dy(const CField&,const geom& g);
CField dz(const CField&,const geom& g);
CField derive(const CField& U, const geom& g, const int& dir);


// Tinyvector operations

template <typename T,int ndim>
TinyVector<T,ndim> operator+ (const TinyVector<T,ndim>& a, const TinyVector<T,ndim>& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a(i) + b(i);
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator+ (const TinyVector<T,ndim>& a, const T& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a(i) + b;
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator+ (const T a, const TinyVector<T,ndim>& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a + b(i);
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator- (const TinyVector<T,ndim>& a, const TinyVector<T,ndim>& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a(i) - b(i);
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator- (const TinyVector<T,ndim>& a, const T& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a(i) - b;
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator- (const T a, const TinyVector<T,ndim>& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a - b(i);
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator* (const TinyVector<T,ndim>& a, const TinyVector<T,ndim>& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a(i) * b(i);
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator* (const TinyVector<T,ndim>& a, const T& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a(i) * b;
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator* (const T a, const TinyVector<T,ndim>& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a * b(i);
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator/ (const TinyVector<T,ndim>& a, const TinyVector<T,ndim>& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a(i) / b(i);
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator/ (const TinyVector<T,ndim>& a, const T& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a(i) / b;
  return res;

};

template <typename T,int ndim>
TinyVector<T,ndim> operator/ (const T a, const TinyVector<T,ndim>& b) {

  TinyVector<T,ndim> res;
  for(int i=0;i<ndim;i++)
    res(i) = a / b(i);
  return res;

};

template<int ndim>
double Dot_prod(const TinyVector<double,ndim>& V1, const TinyVector<double,ndim>& V2) {
  double sum(V1(0)*V2(0));
  for (int i=1;i<ndim;i++)
    sum += V1(i)*V2(i);
  return sum;
}

template<int ndim>
double norm2(const TinyVector<double,ndim>& V) {
  return sqrt(Dot_prod(V,V));
}



#endif
