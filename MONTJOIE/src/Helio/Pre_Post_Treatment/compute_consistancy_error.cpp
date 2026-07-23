#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>

#include <blitz/array.h>

#include <complex>
#include <algorithm>
#include <vector>

#include "utils.hpp"

using namespace std;

// consistancy error is computed on a cartesian grid
// This is why we need to exclude points at a distance lesser than a parameter dist from the boundary
// We also exclude points near the Dirac source

void compute_Helio_scalar_consistancy(CField&, Field&, double&, Field&, Field&, Field&, Field&, Field&, double&, const Vertex&, geom&);

int main (int argc, char** argv) {

  if (argc<2) {
    cout << "Usage : compute_consistancy_error.x solution_file" << endl;
    abort();
  }

  string solution_file(argv[1]);

  geom G(solution_file);

  CField U;                       // tested solution
  Field  rho,sigma,mu,alpha,beta; // background
  Field  err;
  cout << "Read data" << endl;
  Read_VTK(U    ,solution_file      ,G);  
  Read_VTK(rho  ,string("rho.vtk")  ,G);  
  Read_VTK(sigma,string("sigma.vtk"),G);  
  Read_VTK(mu   ,string("mu.vtk")   ,G);  
  Read_VTK(alpha,string("alpha.vtk"),G);  
  Read_VTK(beta ,string("beta.vtk") ,G);  

  double freq = 3.e-3*2.e0*M_PI*698.e6;
  double d = 0.1e0;

  compute_Helio_scalar_consistancy(U,err,freq,rho,sigma,mu,alpha,beta,d,Vertex(0.e0,0.8e0),G);

  return 0;

}


void compute_Helio_scalar_consistancy (CField& U0  , Field& err, double& omega,
                                       Field&  rho0 , Field& sigma0,
                                       Field&  mu0  , Field& alpha0,
                                       Field&  beta0, double& dist, 
                                       const Vertex& s, geom& g0 ) {


  // Reduce the fields to x>0 only

  geom g(g0);
  g.O(0) = 0.e0;
  g.N    = Index(g0.N(0)/2,g0.N(1),g0.N(2));
  g.h    = Vertex(2.e0/(double)g0.N(0),2.e0/(double)g0.N(1),1.e0);

  Index lower(g0.N(0)/2,0);
  Index upper(g0.N(0)-1,g0.N(1)-1);
  RectDomain<NDIM> dom(lower,upper);

  lower = Index(0,0);
  upper = Index(g0.N(0)/2-1,g0.N(1)-1);
  RectDomain<NDIM> dom2(lower,upper);
 
  CField U;  U.resize(g.N);
  U(dom2) = U0(dom);
 
  Field alpha(U.shape()); alpha(dom2) = alpha0(dom);
  Field beta (U.shape()); beta(dom2)  = beta0 (dom);
  Field sigma(U.shape()); sigma(dom2) = sigma0(dom);
  Field rho  (U.shape()); rho(dom2)   = rho0  (dom);
  Field mu   (U.shape()); mu(dom2)    = mu0   (dom);

  CField source(U.shape()); source = Complex(0.e0,0.e0);

// Test computation with some analytical solution
/*

  int N = 250;
  U.resize(N/2,N);
  err.resize(N,N);
  alpha.resize(U.shape()); alpha = 1.e0;
  beta .resize(U.shape()); beta  = 1.e0;
  rho  .resize(U.shape()); rho   = 1.e0;
  sigma.resize(U.shape()); sigma = 0.e0;
  mu   .resize(U.shape()); mu    = 1.e0;
  omega = 1.e0;
  dist = 0.e0;
  g.O(0) = 0.e0;
  g.h = Vertex(2.e0/(double)N,2.e0/(double)N,1.e0);
  g.N = Index(N/2,N,N);

  for(CField::iterator it=source.begin(); it!=source.end(); ++it) {

     Index p = it.position();

     Vertex coords;
     for (int dim=0;dim<NDIM;dim++) {
       coords(dim) = g.O(dim) + p(dim)*g.h(dim);
     }

     double TP = 2.e0*M_PI;
     double x = coords(0);
     double y = coords(1);
     if (x!=0.e0)
       source(p) = Complex(sin(TP*y)*(TP*sin(TP*x)/x + 2.e0*TP*TP*cos(TP*x)-cos(TP*x))-3.e0,0.e0);
     else
       source(p) = Complex(0.e0,0.e0);
     U(p)        = Complex(3.e0+cos(TP*x)*sin(TP*y),0.e0);
  }

*/

  // Build the mask of points where the error will be computed
  
  cout << "Build mask" << endl;
  IField mask  (U.shape());
  Field  radius(U.shape());
  mask = 1;
  int offset = 4;

  for(IField::iterator it=mask.begin(); it!=mask.end(); ++it) {

     Index p = it.position();

     Vertex coords;
     for (int dim=0;dim<NDIM;dim++) {
       coords(dim) = g.O(dim) + p(dim)*g.h(dim);
     }

     double d1 = 1.e0 - norm2(coords); // distance to the boundary
     double d2 = norm2(coords-s);      // distance to the Dirac

     // Remove points close to the source
     //if (d2<dist) mask(p) = 0;
     // Remove points close to the boundary
     //if (d1 < 5.e0*sqrt(2.e0)*g.h(0)) mask(p) = 0;
     // Remove points close to the border of the domain (also, r=0 axis)
     if ((p(1)<offset) || (p(1) > g.N(1)-1-offset) || (p(0)<offset) || (p(0) >g.N(0)-1-offset) ) mask(p) = 0;
     radius(p) = abs(coords(0));

  }
  

  // Compute left hand side : (\omega^2\rho +i\omega\sigma)U

  cout << "Compute LHS" << endl; 
  CField LHS(U.shape());
  LHS  = -omega*omega*rho;
  LHS -= Complex(0.e0,1.e0)*omega*sigma;
  LHS *= U;

  // Compute right hand side : \beta div ( \mu grad(\alpha U))

  cout << "Compute RHS" << endl; 
  CField RHS (U.shape());
  CField tmp1(U.shape());
  CField tmp2(U.shape());

  tmp1  = alpha*U; 
  //tmp1  = U; 
  tmp1  = dx(tmp1,g);
  tmp1 *= mu*radius;
  tmp1  = dx(tmp1,g);
  tmp1 *= where(radius!=0.e0,beta/radius,0.e0);
  //tmp1 *= where(radius!=0.e0,1.e0/radius,0.e0);

  tmp2  = alpha*U; 
  //tmp2  = U; 
  tmp2  = dy(tmp2,g);
  tmp2 *= mu;
  tmp2  = dy(tmp2,g);
  tmp2 *= beta;

  RHS = -tmp1-tmp2;

  // Compute error


  Field MRHS(U.shape()); MRHS = abs(RHS);
  Field MLHS(U.shape()); MLHS = abs(LHS);
  double thr = 1.e-7*max(MLHS);

  RHS = where(MLHS>thr,RHS,0.e0); 
  LHS = where(MLHS>thr,LHS,0.e0); 


  cout << "Compute error" << endl;
  RHS*=mask;
  LHS*=mask;
  source*=mask;

  Field diff(U.shape()); diff = abs(RHS+LHS-source);
  Field val (U.shape()); val  = abs(LHS);

  // Erreur médiane

  vector<double> errs;
  for(Field::iterator it=diff.begin();it!=diff.end();it++) {
    errs.push_back(*it);
  }
  sort(errs.begin(),errs.end());
  double limit = errs[floor(0.8*errs.size())];

  diff = where(diff<=limit,diff,0.e0);
  val  = where(diff<=limit,val ,0.e0);


  err.resize(U.shape());
  err= where(val!=0.e0,diff/val,0.e0);

  cout << "Relative cons. error in L2 norm : " << L2_norm(diff,g)/L2_norm(LHS,g)   << endl;
  //cout << "abs L2 norm of cons. error : " << L2_norm(diff,g)  << endl;
  cout << "Write err.vtk" << endl;
  Write_VTK(err,string("err.vtk"),string("Relative_error"),g);
  cout << "Write sol.vtk (without points in the header!)" << endl;
  Write_VTK(U,string("sol.vtk"),string("U"),g);

};


