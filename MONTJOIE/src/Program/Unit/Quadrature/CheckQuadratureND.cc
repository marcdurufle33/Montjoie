#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Quadrature/MontjoieQuadrature.hxx"

using namespace Montjoie;

Real_wp threshold;

int CheckTriangleQuadrature(int pmax, const VectR2& Points2D, const VectReal_wp& Weights2D)
{
  // testing integral of x^i y^j
  for (int p = 0; p <= pmax; p++)
    for (int i = 0; i <= p; i++)
      {
	int j = p-i;
        double val_exact = tgamma(i+1)*tgamma(j+1)/tgamma(i+j+3);
        double val = 0;
        for (int m = 0; m < Points2D.GetM(); m++)
          val += Weights2D(m)*pow(Points2D(m)(0), i)*pow(Points2D(m)(1), j);
        
        if (abs(val-val_exact)/val > threshold)
          {
	    return p-1;
            //cout << " Order " << i+j << " not exactly integrated for order p = " << p <<endl;
            //DISP(i); DISP(j); DISP(val); DISP(val_exact);
          }
      }  
  
  return pmax;
}

int CheckQuadrangleQuadrature(int pmax, const VectR2& Points2D, const VectReal_wp& Weights2D)
{
  // testing integral of x^i y^j
  for (int p = 0; p <= pmax; p++)
    for (int m1 = 0; m1 <= 2*p; m1++)
      {
	int i = p, j = p;
	if (m1 < p)
	  j = m1;
	else
	  i = m1-p;
	
        double val_exact = Real_wp(1)/((i+1)*(j+1));
        double val = 0;
        for (int m = 0; m < Points2D.GetM(); m++)
          val += Weights2D(m)*pow(Points2D(m)(0), i)*pow(Points2D(m)(1), j);
        
	if (abs(val-val_exact)/val > threshold)
          {
	    return p-1;
            //cout << " Order " << i+j << " not exactly integrated for order p = " << p <<endl;
            //DISP(i); DISP(j); DISP(val); DISP(val_exact);
          }
      }  
  
  return pmax;
}

int CheckTetrahedronQuadrature(int pmax, const VectR3& Points3D, VectReal_wp& Weights3D)
{
  int nb_points = Points3D.GetM();
  // testing integral of x^i y^j z^k
  for (int p = 0; p <= pmax; p++)
    for (int i = 0; i <= p; i++)
      for (int j = 0; j <= p-i; j++)
        {
	  int k = p-i-j;
          double val_exact = tgamma(i+1)*tgamma(j+1)*tgamma(k+1)/tgamma(i+j+k+4);
          double val = 0;
          for (int m = 0; m < nb_points; m++)
            val += Weights3D(m)*pow(Points3D(m)(0), i)*pow(Points3D(m)(1), j)*pow(Points3D(m)(2), k);
          
          if (abs(val-val_exact)/val > threshold)
            {
	      return p-1;
              //cout << " Order " << i+j+k << " not exactly integrated for order p = " << p <<endl;
              //DISP(i); DISP(j); DISP(k); DISP(val); DISP(val_exact);
            }
        }
  
  return pmax;
}

int CheckPyramidQuadrature(int pmax, const VectR3& Points3D, const VectReal_wp& Weights3D, bool poly = false)
{
  // integration de C_p
  for (int k = 0; k <= pmax; k++)
    for (int i = 0; i <= k; i++)
      {
	int jmax = k;
	if (poly)
	  jmax = k-i;
	
	for (int j = 0; j <= jmax; j++)
	  {
	    Real_wp val_ex = 4.0/Real_wp(k+3);
	    if (i%2 == 1)
	      val_ex = 0.0;
	    else
	      val_ex *= 1.0/Real_wp(i+1);
	    
	    if (j%2 == 1)
	      val_ex = 0.0;
	    else
	      val_ex *= 1.0/Real_wp(j+1);
	    
	    Real_wp val_num = 0;
	    for (int n = 0; n < Points3D.GetM(); n++)
	      val_num += Weights3D(n)*pow(Points3D(n)(0), i)*pow(Points3D(n)(1), j)*pow(1.0-Points3D(n)(2), k-i-j);
	    
	    if ((abs(val_num - val_ex) > threshold) || isnan(val_num) || isnan(val_ex))
	      {
		return k-1;
		//cout << "Integrale de x^" << i << " y^" << j << " (1-z)^" << k << " differe" << endl;
		//DISP(val_ex); DISP(val_num);
		//test_integral = false;
	      }
	  }
      }
  
  return pmax;
}

int CheckHexahedronQuadrature(int pmax, const VectR3& Points3D, const VectReal_wp& Weights3D)
{
  // testing integral of x^i y^j z^k
  for (int p = 0; p <= pmax; p++)
    for (int m1 = 0; m1 <= p; m1++) 
      for (int m2 = 0; m2 <= p; m2++)
	for (int m3 = 0; m3 <= 2; m3++)
	  {	    
	    int i = p, j = p, k = p;
	    switch(m3)
	      {
	      case 0: i = m1; j = m2; break;
	      case 1: i = m1; k = m2; break;
	      case 2: j = m1; k = m2; break;
	      }
	    
	    double val_exact = Real_wp(1)/((i+1)*(j+1)*(k+1));
	    double val = 0;
	    for (int m = 0; m < Points3D.GetM(); m++)
	      val += Weights3D(m)*pow(Points3D(m)(0), i)*pow(Points3D(m)(1), j)*pow(Points3D(m)(2), k);
	    
	    if (abs(val-val_exact)/val > threshold)
	      {
		return p-1;
		//cout << " Order " << i+j << " not exactly integrated for order p = " << p <<endl;
		//DISP(i); DISP(j); DISP(val); DISP(val_exact);
	      }
	  }  
  
  return pmax;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 1e-13;
  
  // testing TriangleQuadrature
  VectR2 Points2D; VectReal_wp Weights2D;
  int order_int;
  for (int p = 1; p <= 20; p++)
    {
      TriangleQuadrature::ConstructQuadrature(p, Points2D, Weights2D);
      order_int = CheckTriangleQuadrature(p+2, Points2D, Weights2D);
      if (order_int != p)
	{
	  DISP(order_int); DISP(p);
	  cout << "ConstructQuadrature(triangle) incorrect" << endl;
	  abort();
	}
    }

  for (int p = 1; p <= 20; p++)
    {
      TriangleQuadrature::ConstructQuadrature(p, Points2D, Weights2D, TriangleQuadrature::QUADRATURE_TENSOR);
      order_int = CheckTriangleQuadrature(p+2, Points2D, Weights2D);
      if (order_int != 2*(p/2)+1)
	{
	  DISP(p); DISP(order_int);
	  cout << "ConstructQuadrature(triangle tensor) incorrect" << endl;
	  abort();
	}
    }
  
  // testing QuadrangleQuadrature
  for (int p = 1; p <= 20; p++)
    {
      QuadrangleQuadrature::ConstructQuadrature(p, Points2D, Weights2D);
      order_int = CheckQuadrangleQuadrature(p+2, Points2D, Weights2D);
      if (order_int != 2*(p/2)+1)
	{
	  DISP(p); DISP(order_int);
	  cout << "ConstructQuadrature(quadrangle) incorrect" << endl;
	}
    }
  
  // testing TetrahedronQuadrature
  threshold = 1e-11;
  VectR3 Points3D; VectReal_wp Weights3D;
  for (int p = 1; p <= 21; p++)
    {
      TetrahedronQuadrature::ConstructQuadrature(p, Points3D, Weights3D);
      order_int = CheckTetrahedronQuadrature(p+2, Points3D, Weights3D);
      int pref = p;
      if (p >= 10)
	pref = 2*(p/2)+1;
      
      if (order_int != pref)
	{
	  cout << "ConstructQuadrature(tetrahedron) incorrect" << endl;
	  abort();
	}
    }

  for (int p = 1; p <= 21; p++)
    {
      TetrahedronQuadrature::ConstructQuadrature(p, Points3D, Weights3D, TetrahedronQuadrature::QUADRATURE_TENSOR);
      order_int = CheckTetrahedronQuadrature(p+2, Points3D, Weights3D);
      int pref = 2*(p/2)+1;
      
      if (order_int != pref)
	{
	  cout << "ConstructQuadrature(tetrahedron_tensor) incorrect" << endl;
	  abort();
	}
    }
  
  // testing PyramidQuadrature
  threshold = 1.0e-12;
  for (int p = 1; p <= 9; p++)
    {
      PyramidQuadrature::ConstructPolynomialRule(p, Points3D, Weights3D);
      order_int = CheckPyramidQuadrature(p+2, Points3D, Weights3D, true);
      if (order_int != p)
	{
	  cout << "ConstructQuadrature(pyramid polynomial) incorrect" << endl;
	  abort();
	}
    }
  
  threshold = 0.5e-13;
  VectReal_wp points1d, weights1d;
  for (int p = 1; p <= 10; p++)
    {
      PyramidQuadrature::ConstructQuadrature(p, Points3D, Weights3D, 
					     points1d, weights1d, PyramidQuadrature::QUADRATURE_GAUSS);
      order_int = CheckPyramidQuadrature(2*p+3, Points3D, Weights3D, false);
      if (order_int != 2*p-1)
	{
	  cout << "ConstructQuadrature(pyramid gauss) incorrect" << endl;
	  abort();	  
	}

      PyramidQuadrature::ConstructQuadrature(p, Points3D, Weights3D, 
					     points1d, weights1d, PyramidQuadrature::QUADRATURE_JACOBI1);
      order_int = CheckPyramidQuadrature(2*p+3, Points3D, Weights3D, false);
      if (order_int != 2*p)
	{
	  cout << "ConstructQuadrature(pyramid jacobi1) incorrect" << endl;
	  abort();	  
	}

      PyramidQuadrature::ConstructQuadrature(p, Points3D, Weights3D, 
					     points1d, weights1d, PyramidQuadrature::QUADRATURE_JACOBI2);
      order_int = CheckPyramidQuadrature(2*p+3, Points3D, Weights3D, false);
      if (order_int != 2*p+1)
	{
	  cout << "ConstructQuadrature(pyramid jacobi2) incorrect" << endl;
	  abort();	  
	}
    }
  
  // testing HexahedronQuadrature
  for (int p = 1; p <= 20; p++)
    {
      HexahedronQuadrature::ConstructQuadrature(p, Points3D, Weights3D);
      order_int = CheckHexahedronQuadrature(p+2, Points3D, Weights3D);
      if (order_int != 2*(p/2)+1)
	{
	  DISP(p); DISP(order_int);
	  cout << "ConstructQuadrature(hexahedron) incorrect" << endl;
	}
    }  

  cout << "All tests passed successfully" << endl;    
  
  return FinalizeMontjoie();
}
