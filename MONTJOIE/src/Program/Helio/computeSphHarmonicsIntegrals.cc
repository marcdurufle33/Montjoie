#define MONTJOIE_WITH_TWO_DIM
//#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

#include "Elliptic/Helmholtz/AxiSymHelmholtz.hxx"
#include "Elliptic/Helmholtz/AxiSymHelmholtzInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/AxiSymHelmholtz.cxx"
#endif

#include "Helio/Kernel1D.cxx"
#include <time.h>
using namespace Montjoie;


// ======================================================================================================

int main(int argc, char **argv)
{  
  //InitMontjoie(argc, argv);

  string toCompute = argv[1];
  double prec = 16;

#ifdef SELDON_WITH_MPI
  MPI_Init(&argc, &argv);
#endif

  if (!toCompute.compare("--help") )
    {
      cout << "./computeKernel1D.x keyword options" << endl;
      cout << "\t keyword can be SINGLE_Ylm, SINGLE_Plm, Ylm, Plm, AllTables." << endl;
      cout << "\t ./computeKernel1D.x keyword --help can be used to know the options." << endl;
    }

  else if (!toCompute.compare("SINGLE_Ylm") )
    {
      string str = argv[2];
      if (!str.compare("--help"))
	{
	  cout << "./computeKernel1D.x SINGLE_Ylm l1 l2 l3 m1 m2 m3" << endl;
	  cout << "\t Compute the Gaunt integral int Y_l1^m1 Y_l2^{m2} Y_l3^m3." << endl;
	}
      else
	{     
	  int l1 = to_num<int>(argv[2]);
	  int l2 = to_num<int>(argv[3]);
	  int l3 = to_num<int>(argv[4]);
	  int m1  = to_num<int>(argv[5]);
	  int m2  = to_num<int>(argv[6]);
	  int m3  = to_num<int>(argv[7]);
 
	  Real_wp g = Gaunt(l1,l2,l3,m1,m2,m3);
	  if (argc > 8)
	    prec = atof(argv[8]);
	  cout.precision(prec);

	  cout << g << endl; 
	}   
    }
  else if (!toCompute.compare("3j") )
    {
      string str = argv[2];
      if (!str.compare("--help"))
	{
	  cout << "./computeKernel1D.x 3j l1Max l2 l3 m1 m2 m3" << endl;
	  cout << "\t Returns the Wigner 3j symbol up to l1Max." << endl;
	}
      else
	{     
	  int l1Max = to_num<int>(argv[2]);
	  int l2 = to_num<int>(argv[3]);
	  int l3 = to_num<int>(argv[4]);
	  int m1  = to_num<int>(argv[5]);
	  int m2  = to_num<int>(argv[6]);
	  int m3  = to_num<int>(argv[7]);
 
	  Vector<Real_wp> g = Wigner3j_full(l2,l3,m1,m2,m3,l1Max);
	  if (argc > 8)
	    prec = atof(argv[8]);
	  cout.precision(prec);
	  for (int i=0; i < l1Max; i++)
	    cout << g(i) << endl; 
	}   
    }
  else if (!toCompute.compare("SINGLE_3j") )
    {
      string str = argv[2];
      if (!str.compare("--help"))
	{
	  cout << "./computeKernel1D.x SINGLE_3j l1 l2 l3 m1 m2 m3" << endl;
	  cout << "\t Returns the Wigner 3j symbol." << endl;
	}
      else
	{     
	  int l1 = to_num<int>(argv[2]);
	  int l2 = to_num<int>(argv[3]);
	  int l3 = to_num<int>(argv[4]);
	  int m1  = to_num<int>(argv[5]);
	  int m2  = to_num<int>(argv[6]);
	  int m3  = to_num<int>(argv[7]);
 
	  Real_wp g = Wigner3j(l1,l2,l3,m1,m2,m3);
	  if (argc > 8)
	    prec = atof(argv[8]);
	  cout.precision(prec);

	  cout << g << endl; 
	}   
    }
  else if (!toCompute.compare("Ylm"))
    {  
      string str = argv[2];
      if (!str.compare("--help"))
	{
	  cout << "./computeKernel1D.x Ylm l1Max l2Max l3Min l3Max m1 m2 m3 directory" << endl;
	  cout << "\t Compute the integral int Y_l1^m1 Y_l2^m2 Y_l3^m3 when m1+m2+m3=0." << endl;
	  cout << "\t The computation is made for l1 from 0 to l1Max, l2 from 0 to l2Max and l3 from l3Min to l3Max." << endl; 
	  cout << "\t The results are written in directory." << endl;
	}
      else
	{ 
	  int l1Max = to_num<int>(argv[2]);
	  int l2Max = to_num<int>(argv[3]);
	  int l3Min = to_num<int>(argv[4]);
	  int l3Max = to_num<int>(argv[5]);
	  int m1  = to_num<int>(argv[6]);
	  int m2  = to_num<int>(argv[7]);
	  int m3  = to_num<int>(argv[8]);
	  string directory = argv[9];

	  if (m1+m2+m3 != 0)
	    {
	      cout << "The sum of the m must be equal to 0." << endl;
	      abort();
	    }
	  fullGaunt(l1Max,l2Max,l3Min,l3Max,m1,m2,m3,directory); 
	} 
    }
  else if (!toCompute.compare("SINGLE_Plm"))
    {
      string str = argv[2];
      if (!str.compare("--help"))
	{
	  cout << "./computeKernel1D.x SINGLE_Plm l1 l2 l3 m1 m2 m3" << endl;
	  cout << "\t Compute the integral int P_l1^m1 P_l2^m2 P_l3^m3." << endl;
	  cout << "\t It is implemented for m1+m2+m3=1 or -1." << endl;
	}
      else
	{ 
	  int l1 = to_num<int>(argv[2]);
	  int l2 = to_num<int>(argv[3]);
	  int l3 = to_num<int>(argv[4]);
	  int m1  = to_num<int>(argv[5]);
	  int m2  = to_num<int>(argv[6]);
	  int m3  = to_num<int>(argv[7]);

	  if (argc > 8)
	    prec = atof(argv[8]);
	  cout.precision(prec);

	  Real_wp g = generalGauntSingle(l1,l2,l3,m1,m2,m3);
	  cout << g << endl;
	}
    }
  else if (!toCompute.compare("Plm"))
    {
      string str = argv[2];
      if (!str.compare("--help"))
	{
	  cout << "./computeKernel1D.x Plm l1Max l2Max l3Min l3Max m1 m2 m3 directory" << endl;
	  cout << "\t Compute the integral int P_l1^m1 P_l2^m2 P_l3^m3 when m1+m2+m3=1 or -1." << endl;
	  cout << "\t The computation is made for l1 from 0 to l1Max, l2 from 0 to l2Max and l3 from l3Min to l3Max." << endl; 
	  cout << "\t The results are written in directory." << endl;
	}
      else
	{ 
	  int l1Max = to_num<int>(argv[2]);
	  int l2Max = to_num<int>(argv[3]);
	  int l3Min = to_num<int>(argv[4]);
	  int l3Max = to_num<int>(argv[5]);
	  int m1  = to_num<int>(argv[6]);
	  int m2  = to_num<int>(argv[7]);
	  int m3  = to_num<int>(argv[8]);
	  string directory;
	  if (argc > 9)
	    directory = argv[9];
	  else
	    directory = string();
	  generalGaunt(l1Max,l2Max,l3Min,l3Max,m1,m2,m3,directory);
 
	  /*
	    for (int i1=0; i1<l1Max+1;i1++)
	    for (int i2=0; i2<l2Max+1;i2++)
	    for (int i3=l3Min; i3<l3Max+1;i3++)
	    cout << "gaunt(" << i3 << "," << i1 << "," << i2 << ") =" << g(i3-l3Min)(i1,i2) << endl;
	  */
	}
    }

  else if (!toCompute.compare("AllTables"))
    {
      string str = argv[2];
      if (!str.compare("--help"))
	{
	  cout << "./computeKernel1D.x AllTables l1Max l2Max l3Min l3Max m1 m2 m3 directory" << endl;
	  cout << "\t Compute the integral int Y_l1^m Y_l2^-m Y_l3^0, int P_l1^m P_l2^1-m P_l3^0 and int P_l1^m P_l2^-1-m P_l3^0." << endl;
	  cout << "\t The computation is made for l1 from 0 to l1Max, l2 from 0 to l2Max, l3 from l3Min to l3Max and m from 0 to l1Max." << endl; 
	  cout << "\t The results are written in directory." << endl;
	}
      else
	{
	  int l1Max = to_num<int>(argv[2]);
	  int l2Max = to_num<int>(argv[3]);
	  int l3Max = to_num<int>(argv[4]);
	  int shift  = to_num<int>(argv[5]);
	  string directory = argv[6];
	  for (int m=0; m<l1Max; m++)
	    {
	      if (shift == 0)
		fullGaunt(l1Max,l2Max,0,l3Max,m,-m,0,directory); 
	      else
		generalGaunt(l1Max,l2Max,0,l3Max,m,-m+shift,0,directory);
	    }
	}
    }
  else if (!toCompute.compare("Kernel"))
    { 
  	  Kernel1D kernel = Kernel1D();
	  int indFreqMin = 0;
	  if (argc > 3)
	    indFreqMin = to_num<int>(argv[3]);
	  kernel.ReadInputFile(argv[2], indFreqMin);
	  kernel.compute(0);
    }
  else if (!toCompute.compare("Kernel3D"))
    { 
  	  Kernel1D kernel = Kernel1D();
	  int indFreqMin = 0;
	  if (argc > 3)
	    indFreqMin = to_num<int>(argv[3]);
	  kernel.ReadInputFile(argv[2], indFreqMin);
	  kernel.computeKernel3D();
    }

  else if (!toCompute.compare("Xu"))
    { 
      int lMax = to_num<int>(argv[2]);
      int lb = to_num<int>(argv[3]);
      int mb  = to_num<int>(argv[4]);
      //int mu  = to_num<int>(argv[5]);
      

      clock_t begin_time = clock();
      for (int l=0; l < lMax; l++)
	{
	  for  (int m=-l; m <=l ; m++)
	    {
	      Vector<Real_wp> g = Xu_getGaunt(lb, l, mb, m);
	    }
	}
      Real_wp secondsXu = 1.*(clock() - begin_time) / CLOCKS_PER_SEC;
      begin_time = clock();
      for (int l=0; l < lMax; l++)
	{
	  for  (int m=-l; m <=l ; m++)
	    {
	      Vector<Real_wp> g0 = GauntVect(lMax, l, lb, -m-mb, m, mb);
	    }
	}
      Real_wp secondsG = 1.*(clock() - begin_time) / CLOCKS_PER_SEC;


      /*
      //for (int i=0; i < g.GetM(); i++)
      //	cout << g(i) << endl;

      // Compare with old version
      ///Real_wp cl1m1 = sqrt((2.*n+1)/(4.*pi_wp) / Xu_RatioFactorial(n,m)); //factorial(n-m)/factorial(n+m));
      //Real_wp cl2m2 = sqrt((2.*nu+1)/(4.*pi_wp) / Xu_RatioFactorial(nu,mu)); //factorial(nu-mu)/factorial(nu+mu));
      int qmax = min(min(n,nu),(n+nu-abs(m+mu))/2);
      for (int p=0; p <= qmax; p++)
	{
	  int l3 = n+nu-2*p;
	  //Real_wp cl3m3 = sqrt((2.*l3+1)/(4.*pi_wp) / Xu_RatioFactorial(l3,m+mu)); // factorial(l3-(m+mu))/factorial(l3+(m+mu)));
	  Real_wp gTest = Gaunt(n,nu,l3,m,mu,-m-mu);
	  //gTest /= (2.*pi_wp*cl1m1*cl2m2*cl3m3) * pow(-1.,m+mu);
	  //gTest /= (2.*pi_wp) * pow(-1.,m+mu);
	  //gTest *= (2.*l3+1.)/2. / Xu_RatioFactorial(l3,m+mu); //factorial(l3-m-mu) / factorial(l3+m+mu);
	  //cout << "g = " << gTest << ", new g = " << g(p) << " , ratio = " << (gTest / g(p)) << endl;
	}
      */

      cout << "Time Xu = " << secondsXu << ", time W3j = " << secondsG << endl;
    }

  else if (!toCompute.compare("Plm2"))
    {
      string str = argv[2];
      if (!str.compare("--help"))
	{
	  cout << "./computeKernel1D.x Plm l1Max l2Max l3Min l3Max m1 m2 m3 directory" << endl;
	  cout << "\t Compute the integral int P_l1^m1 P_l2^m2 P_l3^m3 when m1+m2+m3=1 or -1." << endl;
	  cout << "\t The computation is made for l1 from 0 to l1Max, l2 from 0 to l2Max and l3 from l3Min to l3Max." << endl; 
	  cout << "\t The results are written in directory." << endl;
	}
      else
	{ 
	  int l1Max = to_num<int>(argv[2]);
	  int l2Max = to_num<int>(argv[3]);
	  int l3Min = to_num<int>(argv[4]);
	  int l3Max = to_num<int>(argv[5]);
	  int m3  = to_num<int>(argv[6]);
	  Vector<int> l3s(l3Max-l3Min);
	  for (int i=0; i < l3Max-l3Min; i++)
	    l3s(i) = l3Min+i;
	  generalGauntNew(l1Max,l2Max,l3s,m3,1);
 
	  /*
	    for (int i1=0; i1<l1Max+1;i1++)
	    for (int i2=0; i2<l2Max+1;i2++)
	    for (int i3=l3Min; i3<l3Max+1;i3++)
	    cout << "gaunt(" << i3 << "," << i1 << "," << i2 << ") =" << g(i3-l3Min)(i1,i2) << endl;
	  */
	}
    }

  else
    {
      cout << "This code needs a data file and a file of list of sources in argument" << endl;
      cout << "helmholtz_axi.x configurationFileName srcListFileName" << endl;
      cout << "is a good syntax" << endl;
    }
  //return FinalizeMontjoie();

#ifdef SELDON_WITH_MPI
    MPI_Finalize();
#endif

  return 0;
}
