#include "Solver/MontjoieSolver.hxx"
#include "Output/MontjoieOutput.hxx"
#include "Instationary/MontjoieTime.hxx"

//#define USE_DERIVATIVE_CFL

using namespace Montjoie;

class RknOrderCondition
{
protected :
  int order, nb_eq;
  bool bbar_from_b; bool cond_a_ci2;
  
public :
  void SetOrder(int r, bool cond_bbar, bool cond_a)
  {
    order = r;
    nb_eq = 0;
    bbar_from_b = cond_bbar;
    cond_a_ci2 = cond_a;
    if (cond_bbar)
      {
	if (cond_a)
	  {
	    switch(r)
	      {
	      case 1 : nb_eq = 1; break;
	      case 2 : nb_eq = 2; break;
	      case 3 : nb_eq = 3; break;
	      case 4 : nb_eq = 5; break;
	      case 5 : nb_eq = 8; break;
	      case 6 : nb_eq = 13; break;
	      case 7 : nb_eq = 22; break;
	      case 8 : nb_eq = 37; break;
	      case 9 : nb_eq = 63; break;
	      case 10 : nb_eq = 111; break;
	      default : cout << "not implemented" << endl; abort();
	      }
	  }
	else
	  {
	    switch(r)
	      {
	      case 1 : nb_eq = 1; break;
	      case 2 : nb_eq = 3; break;
	      case 3 : nb_eq = 5; break;
	      case 4 : nb_eq = 8; break;
	      case 5 : nb_eq = 14; break;
	      case 6 : nb_eq = 24; break;
	      case 7 : nb_eq = 44; break;
	      case 8 : nb_eq = 80; break;
	      default : cout << "not implemented" << endl; abort();
	      }
	  }	
      }
    else
      {
	if (cond_a)
	  {
	    cout << "Condition sum bar{a}_{i, j} = c_i^2 / 2 "
		 << "cannot be taken if bar{b}_i ne b_i (1 - c_i)" << endl;
	    
	    abort();
	  }
	
	switch(r)
	  {
	  case 1 : nb_eq = 1; break;
	  case 2 : nb_eq = 3; break;
	  case 3 : nb_eq = 6; break;
	  case 4 : nb_eq = 11; break;
	  case 5 : nb_eq = 20; break;
	  case 6 : nb_eq = 36; break;
	  case 7 : nb_eq = 66; break;
	  case 8 : nb_eq = 122; break;
	  default : cout << "not implemented" << endl; abort();
	  }
      }
	
  }
  
  int GetM() const
  {
    return nb_eq;
  }
  
  void GetCoefficients(const VectReal_wp& param, Matrix<Real_wp>& Abar, VectReal_wp& B,
		       VectReal_wp& Bbar, VectReal_wp& C) const
  {
  }
  
  void EvaluateFunction(const VectReal_wp& param, 
			VectReal_wp& eval) const
  {
    Matrix<Real_wp> Abar;
    VectReal_wp B, Bbar, C;
    GetCoefficients(param, Abar, B, Bbar, C);
    
    EvaluateFunction(Abar, B, Bbar, C, eval);
  }

  void EvaluateFunction(const Matrix<Real_wp>& Abar, const VectReal_wp& B,
			const VectReal_wp& Bbar, const VectReal_wp& C, 
			VectReal_wp& eval) const
  {
    eval.Fill(0);
    
    int n = C.GetM();
    Real_wp one(1);
    int nb = 0;
    if (order >= 1)
      {
	eval(nb) = one;
	for (int i = 0; i < n; i++)
	  eval(nb) -= B(i);
	
	nb++;
      }

    if (order >= 2)
      {
	eval(nb) = one/2;
	for (int i = 0; i < n; i++)
	  eval(nb) -= B(i)*C(i);
	
	nb++;
	
	if (!bbar_from_b)
	  {
	    eval(nb) = one/2;
	    for (int i = 0; i < n; i++)
	      eval(nb) -= Bbar(i);
	    
	    nb++;
	  }
      }
    
    if (order >= 3)
      {	
	eval(nb) = one/3;
	for (int i = 0; i < n; i++)
	  eval(nb) -= B(i)*C(i)*C(i);
	
	nb++;
	
	if (!cond_a_ci2)
	  {
	    eval(nb) = one/6;
	    for (int i = 0; i < n; i++)
	      for (int j = 0; j < i; j++)
		eval(nb) -= Abar(i, j)*B(i);  

	    nb++;
	  }

	if (!bbar_from_b)
	  {
	    eval(nb) = one/6;
	    for (int i = 0; i < n; i++)
	      eval(nb) -= Bbar(i)*C(i);
	    
	    nb++;
	  }
      }
    
    if (order >= 4)
      {
	eval(nb) = one/4; eval(nb+1) = one/24;
	for (int i = 0; i < n; i++)
	  {
	    eval(nb) -= B(i)*pow(C(i), 3);
	    for (int j = 0; j < i; j++)
	      eval(nb+1) -= B(i)*Abar(i, j)*C(j);
	  }
	
	nb += 2;
	
	if (!cond_a_ci2)
	  {
	    eval(nb) = one/8;
	    for (int i = 0; i < n; i++)
	      for (int j = 0; j < i; j++)
		eval(nb) -= B(i)*Abar(i, j)*C(i);
	    
	    nb++;
	  }
	
	if (!bbar_from_b)
	  {
	    eval(nb) = one/12; eval(nb+1) = one/24;
	    for (int i = 0; i < n; i++)
	      {	    
		eval(nb) -= Bbar(i)*C(i)*C(i);
		for (int j = 0; j < i; j++)
		  eval(nb+1) -= Bbar(i)*Abar(i, j);
	      }
	    
	    nb += 2;
	  }
      }
    
    if (order >= 5)
      {
	eval(nb) = one/5; eval(nb+1) = one/30; eval(nb+2) = one/60; 
	for (int i = 0; i < n; i++)
	  {
	    eval(nb) -= B(i)*pow(C(i), 4);
	    for (int j = 0; j < i; j++)
	      {
		eval(nb+1) -= B(i)*Abar(i, j)*C(i)*C(j);
		eval(nb+2) -= B(i)*Abar(i, j)*C(j)*C(j);
	      }
	  }
	
	nb += 3;

	if (!cond_a_ci2)
	  {
	    eval(nb) = one/10; eval(nb+1) = one/120; eval(nb+2) = one/20;
	    for (int i = 0; i < n; i++)
	      {
		for (int j = 0; j < i; j++)
		  {
		    eval(nb) -= B(i)*Abar(i, j)*C(i)*C(i);
		    for (int k = 0; k < j; k++)
		      eval(nb+1) -= B(i)*Abar(i, j)*Abar(j, k);
		    
		    for (int k = 0; k < i; k++)
		      eval(nb+2) -= B(i)*Abar(i, j)*Abar(i, k);
		  }
	      }
	    
	    nb += 3;
	  }
	
	if (!bbar_from_b)
	  {
	    eval(nb) = one/20; eval(nb+1) = one/40; eval(nb+2) = one/120;
	    for (int i = 0; i < n; i++)
	      {
		eval(nb) -= Bbar(i)*pow(C(i), 3);
		for (int j = 0; j < i; j++)
		  {
		    eval(nb+1) -= Bbar(i)*Abar(i, j)*C(i);
		    eval(nb+2) -= Bbar(i)*Abar(i, j)*C(j);
		  }
	      }
	    
	    nb += 3;
	  }
      }
    
    if (order >= 6)
      {
	eval(nb) = one/6; eval(nb+1) = one/36; eval(nb+2) = one/72;
	eval(nb+3) = one/120; eval(nb+4) = one/720;
	for (int i = 0; i < n; i++)
	  {
	    eval(nb) -= B(i)*pow(C(i), 5);
	    for (int j = 0; j < i; j++)
	      {
		eval(nb+1) -= B(i)*Abar(i, j)*C(i)*C(i)*C(j);
		eval(nb+2) -= B(i)*Abar(i, j)*C(i)*C(j)*C(j);
		eval(nb+3) -= B(i)*Abar(i, j)*pow(C(j), 3);
		for (int k = 0; k < j; k++)
		  eval(nb+4) -= B(i)*Abar(i, j)*Abar(j, k)*C(k);
	      }
	  }
	
	nb += 5;
	
	if (!cond_a_ci2)
	  {
	    eval(nb) = one/12; eval(nb+1) = one/24; eval(nb+2) = one/144;
	    eval(nb+3) = one/72; eval(nb+4) = one/240;
	    for (int i = 0; i < n; i++)
	      for (int j = 0; j < i; j++)
		{
		  eval(nb) -= B(i)*Abar(i, j)*pow(C(i), 3);
		  for (int k = 0; k < i; k++)
		    {
		      eval(nb+1) -= B(i)*Abar(i, j)*Abar(i, k)*C(i);
		      eval(nb+3) -= B(i)*Abar(i, j)*Abar(i, k)*C(j);
		    }
		  
		  for (int k = 0; k < j; k++)
		    {
		      eval(nb+2) -= B(i)*Abar(i, j)*Abar(j, k)*C(i);
		      eval(nb+4) -= B(i)*Abar(i, j)*Abar(j, k)*C(j);
		    }
		}
	    
	    nb += 5;
	  }
	
	if (!bbar_from_b)
	  {
	    eval(nb) = one/30; eval(nb+1) = one/60; eval(nb+2) = one/180;
	    eval(nb+3) = one/360; eval(nb+4) = one/720; eval(nb+5) = one/120;
	    for (int i = 0; i < n; i++)
	      {
		eval(nb) -= Bbar(i)*pow(C(i), 4);
		for (int j = 0; j < i; j++)
		  {
		    eval(nb+1) -= Bbar(i)*Abar(i, j)*C(i)*C(i);
		    eval(nb+2) -= Bbar(i)*Abar(i, j)*C(i)*C(j);
		    eval(nb+3) -= Bbar(i)*Abar(i, j)*C(j)*C(j);
		    for (int k = 0; k < j; k++)
		      eval(nb+4) -= Bbar(i)*Abar(i, j)*Abar(j, k);
		    
		    for (int k = 0; k < i; k++)
		      eval(nb+5) -= Bbar(i)*Abar(i, j)*Abar(i, k);
		  }
	      }
	    
	    nb += 6;
	  }
      }
    
    if (order >= 7)
      {
	eval(nb) = one/7; eval(nb+1) = one/42; eval(nb+2) = one/252;
	eval(nb+3) = one/140; eval(nb+4) = one/84; eval(nb+5) = one/840;
	eval(nb+6) = one/210; eval(nb+7) = one/1260; eval(nb+8) = one/2520;
	for (int i = 0; i < n; i++)
	  {
	    eval(nb) -= B(i)*pow(C(i), 6);
	    for (int j = 0; j < i; j++)
	      {
		eval(nb+1) -= B(i)*pow(C(i), 3)*Abar(i, j)*C(j);
		eval(nb+3) -= B(i)*C(i)*Abar(i, j)*pow(C(j), 3);
		eval(nb+4) -= B(i)*C(i)*C(i)*Abar(i, j)*C(j)*C(j);
		eval(nb+6) -= B(i)*Abar(i, j)*pow(C(j), 4);
		for (int k = 0; k < i; k++)
		  eval(nb+2) -= B(i)*Abar(i, j)*Abar(i, k)*C(j)*C(k);
		
		for (int k = 0; k < j; k++)
		  {
		    eval(nb+5) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*C(k);
		    eval(nb+7) -= B(i)*Abar(i, j)*C(j)*Abar(j, k)*C(k);
		    eval(nb+8) -= B(i)*Abar(i, j)*Abar(j, k)*C(k)*C(k);
		  }
	      }
	  }
	
	nb += 9;
	
	if (!cond_a_ci2)
	  {
	    eval(nb) = one/14; eval(nb+1) = one/28; eval(nb+2) = one/56;
	    eval(nb+3) = one/84; eval(nb+4) = one/168; eval(nb+5) = one/336;
	    eval(nb+6) = one/168; eval(nb+7) = one/280; eval(nb+8) = one/420;
	    eval(nb+9) = one/5040; eval(nb+10) = one/840;
	    for (int i = 0; i < n; i++)
	      for (int j = 0; j < i; j++)
		{
		  eval(nb) -= B(i)*pow(C(i), 4)*Abar(i, j);
		  for (int k = 0; k < i; k++)
		    {
		      eval(nb+1) -= B(i)*C(i)*C(i)*Abar(i, j)*Abar(i, k);
		      eval(nb+3) -= B(i)*C(i)*Abar(i, j)*Abar(i, k)*C(j);
		      eval(nb+6) -= B(i)*Abar(i, j)*Abar(i, k)*C(j)*C(j);
		      for (int l = 0; l < i; l++)
			eval(nb+2) -= B(i)*Abar(i, j)*Abar(i, k)*Abar(i, l);
		      
		      for (int l = 0; l < k; l++)
			eval(nb+5) -= B(i)*Abar(i, j)*Abar(i, k)*Abar(k, l);
		    }
		  
		  for (int k = 0; k < j; k++)
		    {
		      eval(nb+4) -= B(i)*C(i)*C(i)*Abar(i, j)*Abar(j, k);
		      eval(nb+7) -= B(i)*C(i)*Abar(i, j)*C(j)*Abar(j, k);
		      eval(nb+8) -= B(i)*Abar(i, j)*C(j)*C(j)*Abar(j, k);
		      for (int l = 0; l < k; l++)
			eval(nb+9) -= B(i)*Abar(i, j)*Abar(j, k)*Abar(k, l);

		      for (int l = 0; l < j; l++)
			eval(nb+10) -= B(i)*Abar(i, j)*Abar(j, k)*Abar(j, l);
		    }
		}
	    
	    nb += 11;
	  }
	
	if (!bbar_from_b)
	  {
	    eval(nb) = one/42; eval(nb+1) = one/252; eval(nb+2) = one/504;
	    eval(nb+3) = one/840; eval(nb+4) = one/5040;
	    for (int i = 0; i < n; i++)
	      {
		eval(nb) -= Bbar(i)*pow(C(i), 5);
		for (int j = 0; j < i; j++)
		  {
		    eval(nb+1) -= Bbar(i)*Abar(i, j)*C(i)*C(i)*C(j);
		    eval(nb+2) -= Bbar(i)*Abar(i, j)*C(i)*C(j)*C(j);
		    eval(nb+3) -= Bbar(i)*Abar(i, j)*pow(C(j), 3);
		    for (int k = 0; k < j; k++)
		      eval(nb+4) -= Bbar(i)*Abar(i, j)*Abar(j, k)*C(k);
		  }
	      }
	    
	    nb += 5;
	    
	    eval(nb) = one/84; eval(nb+1) = one/168; eval(nb+2) = one/1008;
	    eval(nb+3) = one/504; eval(nb+4) = one/1680;
	    for (int i = 0; i < n; i++)
	      for (int j = 0; j < i; j++)
		{
		  eval(nb) -= Bbar(i)*Abar(i, j)*pow(C(i), 3);
		  for (int k = 0; k < i; k++)
		    {
		      eval(nb+1) -= Bbar(i)*Abar(i, j)*Abar(i, k)*C(i);
		      eval(nb+3) -= Bbar(i)*Abar(i, j)*Abar(i, k)*C(j);
		    }
		  
		  for (int k = 0; k < j; k++)
		    {
		      eval(nb+2) -= Bbar(i)*Abar(i, j)*Abar(j, k)*C(i);
		      eval(nb+4) -= Bbar(i)*Abar(i, j)*Abar(j, k)*C(j);
		    }
		}
	    
	    nb += 5;
	  }
      }
    
    if (order >= 8)
      {
	eval(nb) = one/8; eval(nb+1) = one/48; eval(nb+2) = one/96;
	eval(nb+3) = one/288; eval(nb+4) = one/960; eval(nb+5) = one/576;
	eval(nb+6) = one/160; eval(nb+7) = one/240; eval(nb+8) = one/1440;
	eval(nb+9) = one/2880; eval(nb+10) = one/336; eval(nb+11) = one/4032;
	eval(nb+12) = one/2016; eval(nb+13) = one/6720; eval(nb+14) = one/40320;
	for (int i = 0; i < n; i++)
	  {
	    eval(nb) -= B(i)*pow(C(i), 7);
	    for (int j = 0; j < i; j++)
	      {
		eval(nb+1) -= B(i)*pow(C(i), 4)*Abar(i, j)*C(j);
		eval(nb+2) -= B(i)*pow(C(i), 3)*Abar(i, j)*C(j)*C(j);
		eval(nb+6) -= B(i)*C(i)*C(i)*Abar(i, j)*pow(C(j), 3);
		eval(nb+7) -= B(i)*C(i)*Abar(i, j)*pow(C(j), 4);
		eval(nb+10) -= B(i)*Abar(i, j)*pow(C(j), 5);
		for (int k = 0; k < i; k++)
		  {
		    eval(nb+3) -= B(i)*C(i)*Abar(i, j)*C(j)*Abar(i, k)*C(k);
		    eval(nb+5) -= B(i)*Abar(i, j)*Abar(i, k)*C(j)*C(k)*C(k);
		  }
		
		for (int k = 0; k < j; k++)
		  {
		    eval(nb+4) -= B(i)*C(i)*C(i)*Abar(i, j)*Abar(j, k)*C(k);
		    eval(nb+8) -= B(i)*C(i)*Abar(i, j)*C(j)*Abar(j, k)*C(k);
		    eval(nb+9) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*C(k)*C(k);
		    eval(nb+11) -= B(i)*Abar(i, j)*C(j)*Abar(j, k)*C(k)*C(k);
		    eval(nb+12) -= B(i)*Abar(i, j)*C(j)*C(j)*Abar(j, k)*C(k);
		    eval(nb+13) -= B(i)*Abar(i, j)*Abar(j, k)*pow(C(k), 3);
		    for (int l = 0; l < k; l++)
		      eval(nb+14) -= B(i)*Abar(i, j)*Abar(j, k)*Abar(k, l)*C(l);
		  }
	      }	   
	  }
	
	nb += 15;
	
	if (!cond_a_ci2)
	  {
	    eval(nb) = one/16; eval(nb+1) = one/32; eval(nb+2) = one/64;
	    eval(nb+3) = one/96; eval(nb+4) = one/192; eval(nb+5) = one/192;
	    eval(nb+6) = one/192; eval(nb+7) = one/384; eval(nb+8) = one/320;
	    eval(nb+9) = one/1152; eval(nb+10) = one/320; eval(nb+11) = one/640;
	    eval(nb+12) = one/1920; eval(nb+13) = one/480; eval(nb+14) = one/5760;
	    eval(nb+15) = one/960; eval(nb+16) = one/672; eval(nb+17) = one/1344;
	    eval(nb+18) = one/8064; eval(nb+19) = one/4032; eval(nb+20) = one/13440;
	    for (int i = 0; i < n; i++)
	      for (int j = 0; j < i; j++)
		{
		  eval(nb) -= B(i)*Abar(i, j)*pow(C(i), 5);
		  for (int k = 0; k < i; k++)
		    {
		      eval(nb+1) -= B(i)*pow(C(i), 3)*Abar(i, j)*Abar(i, k);
		      eval(nb+3) -= B(i)*C(i)*C(i)*Abar(i, j)*C(j)*Abar(i, k);
		      eval(nb+6) -= B(i)*C(i)*Abar(i, j)*C(j)*C(j)*Abar(i, k);
		      eval(nb+10) -= B(i)*Abar(i, j)*pow(C(j), 3)*Abar(i, k);
		      for (int l = 0; l < i; l++)
			{
			  eval(nb+2) -= B(i)*C(i)*Abar(i, j)*Abar(i, k)*Abar(i, l);
			  eval(nb+5) -= B(i)*Abar(i, j)*C(j)*Abar(i, k)*Abar(i, l);
			}
		      
		      for (int l = 0; l < j; l++)
			{
			  eval(nb+9) -= B(i)*Abar(i, j)*Abar(i, k)*C(k)*Abar(j, l);
			  eval(nb+11) -= B(i)*Abar(i, j)*Abar(i, k)*Abar(j, l)*C(j);
			  eval(nb+12) -= B(i)*Abar(i, j)*Abar(i, k)*Abar(j, l)*C(l);
			}
		    }
		  
		  for (int k = 0; k < j; k++)
		    {
		      eval(nb+4) -= B(i)*pow(C(i), 3)*Abar(i, j)*Abar(j, k);
		      eval(nb+8) -= B(i)*C(i)*C(i)*Abar(i, j)*C(j)*Abar(j, k);
		      eval(nb+13) -= B(i)*C(i)*Abar(i, j)*C(j)*C(j)*Abar(j, k);
		      eval(nb+16) -= B(i)*Abar(i, j)*pow(C(j), 3)*Abar(j, k);
		      for (int l = 0; l < i; l++)
			eval(nb+7) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*Abar(i, l);
		      
		      for (int l = 0; l < j; l++)
			{
			  eval(nb+15) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*Abar(j, l);
			  eval(nb+17) -= B(i)*Abar(i, j)*C(j)*Abar(j, k)*Abar(j, l);
			  eval(nb+19) -= B(i)*Abar(i, j)*Abar(j, k)*C(k)*Abar(j, l);
			}
		      
		      for (int l = 0; l < k; l++)
			{
			  eval(nb+14) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*Abar(k, l);
			  eval(nb+18) -= B(i)*Abar(i, j)*C(j)*Abar(j, k)*Abar(k, l);
			  eval(nb+20) -= B(i)*Abar(i, j)*Abar(j, k)*C(k)*Abar(k, l);
			}
		    }
		}
	    
	    nb += 21;
	  }
	
	if (!bbar_from_b)
	  {
	    eval(nb) = one/56; eval(nb+1) = one/336; eval(nb+2) = one/2016;
	    eval(nb+3) = one/1120; eval(nb+4) = one/672; eval(nb+5) = one/6720;
	    eval(nb+6) = one/1680; eval(nb+7) = one/10080; eval(nb+8) = one/20160;
	    for (int i = 0; i < n; i++)
	      {
		eval(nb) -= Bbar(i)*pow(C(i), 6);
		for (int j = 0; j < i; j++)
		  {
		    eval(nb+1) -= Bbar(i)*pow(C(i), 3)*Abar(i, j)*C(j);
		    eval(nb+3) -= Bbar(i)*C(i)*Abar(i, j)*pow(C(j), 3);
		    eval(nb+4) -= Bbar(i)*C(i)*C(i)*Abar(i, j)*C(j)*C(j);
		    eval(nb+6) -= Bbar(i)*Abar(i, j)*pow(C(j), 4);
		    for (int k = 0; k < i; k++)
		      eval(nb+2) -= Bbar(i)*Abar(i, j)*Abar(i, k)*C(j)*C(k);
		    
		    for (int k = 0; k < j; k++)
		      {
			eval(nb+5) -= Bbar(i)*C(i)*Abar(i, j)*Abar(j, k)*C(k);
			eval(nb+7) -= Bbar(i)*Abar(i, j)*C(j)*Abar(j, k)*C(k);
			eval(nb+8) -= Bbar(i)*Abar(i, j)*Abar(j, k)*C(k)*C(k);
		      }
		  }
	      }
	    
	    nb += 9;
	    
	    eval(nb) = one/112; eval(nb+1) = one/224; eval(nb+2) = one/448;
	    eval(nb+3) = one/672; eval(nb+4) = one/1344; eval(nb+5) = one/2688;
	    eval(nb+6) = one/1344; eval(nb+7) = one/2240; eval(nb+8) = one/3360;
	    eval(nb+9) = one/40320; eval(nb+10) = one/6720;
	    for (int i = 0; i < n; i++)
	      for (int j = 0; j < i; j++)
		{
		  eval(nb) -= Bbar(i)*pow(C(i), 4)*Abar(i, j);
		  for (int k = 0; k < i; k++)
		    {
		      eval(nb+1) -= Bbar(i)*C(i)*C(i)*Abar(i, j)*Abar(i, k);
		      eval(nb+3) -= Bbar(i)*C(i)*Abar(i, j)*Abar(i, k)*C(j);
		      eval(nb+6) -= Bbar(i)*Abar(i, j)*Abar(i, k)*C(j)*C(j);
		      for (int l = 0; l < i; l++)
			eval(nb+2) -= Bbar(i)*Abar(i, j)*Abar(i, k)*Abar(i, l);
		      
		      for (int l = 0; l < k; l++)
			eval(nb+5) -= Bbar(i)*Abar(i, j)*Abar(i, k)*Abar(k, l);
		    }
		  
		  for (int k = 0; k < j; k++)
		    {
		      eval(nb+4) -= Bbar(i)*C(i)*C(i)*Abar(i, j)*Abar(j, k);
		      eval(nb+7) -= Bbar(i)*C(i)*Abar(i, j)*C(j)*Abar(j, k);
		      eval(nb+8) -= Bbar(i)*Abar(i, j)*C(j)*C(j)*Abar(j, k);
		      for (int l = 0; l < k; l++)
			eval(nb+9) -= Bbar(i)*Abar(i, j)*Abar(j, k)*Abar(k, l);
		      
		      for (int l = 0; l < j; l++)
			eval(nb+10) -= Bbar(i)*Abar(i, j)*Abar(j, k)*Abar(j, l);
		    }
		}
	    
	    nb += 11;
	  }
      }
    
    if (order >= 9)
      {
	eval(nb) = one/9; eval(nb+1) = one/54; eval(nb+2) = one/108; eval(nb+3) = one/324;
	eval(nb+4) = one/180; eval(nb+5) = one/1080; eval(nb+6) = one/648; eval(nb+7) = one/270;
	eval(nb+8) = one/1620; eval(nb+9) = one/3240; eval(nb+10) = one/1296;
	eval(nb+11) = one/1080; eval(nb+12) = one/378; eval(nb+13) = one/2268;
	eval(nb+14) = one/4536; eval(nb+15) = one/7560; eval(nb+16) = one/45360;
	eval(nb+17) = one/504; eval(nb+18) = one/3024; eval(nb+19) = one/18144;
	eval(nb+20) = one/10080; eval(nb+21) = one/6048; eval(nb+22) = one/60480;
	eval(nb+23) = one/15120; eval(nb+24) = one/90720; eval(nb+25) = one/181440;
	for (int i = 0; i < n; i++)
	  {
	    eval(nb) -= B(i)*pow(C(i), 8);
	    for (int j = 0; j < i; j++)
	      {
		eval(nb+1) -= B(i)*pow(C(i), 5)*Abar(i, j)*C(j);
		eval(nb+2) -= B(i)*pow(C(i), 4)*Abar(i, j)*C(j)*C(j);
		eval(nb+4) -= B(i)*pow(C(i), 3)*Abar(i, j)*pow(C(j), 3);
		eval(nb+7) -= B(i)*C(i)*C(i)*Abar(i, j)*pow(C(j), 4);
		eval(nb+12) -= B(i)*C(i)*Abar(i, j)*pow(C(j), 5);
		eval(nb+17) -= B(i)*Abar(i, j)*pow(C(j), 6);
		for (int k = 0; k < i; k++)
		  {
		    eval(nb+3) -= B(i)*C(i)*C(i)*Abar(i, j)*Abar(i, k)*C(j)*C(k);
		    eval(nb+6) -= B(i)*C(i)*Abar(i, j)*Abar(i, k)*C(j)*C(j)*C(k);
		    eval(nb+10) -= B(i)*Abar(i, j)*C(j)*C(j)*Abar(i, k)*C(k)*C(k);
		    eval(nb+11) -= B(i)*Abar(i, j)*pow(C(j), 3)*Abar(i, k)*C(k);
		  }
		
		for (int k = 0; k < j; k++)
		  {
		    eval(nb+5) -= B(i)*pow(C(i), 3)*Abar(i, j)*Abar(j, k)*C(k);
		    eval(nb+8) -= B(i)*C(i)*C(i)*Abar(i, j)*C(j)*Abar(j, k)*C(k);
		    eval(nb+9) -= B(i)*C(i)*C(i)*Abar(i, j)*Abar(j, k)*C(k)*C(k);
		    eval(nb+13) -= B(i)*C(i)*Abar(i, j)*C(j)*C(j)*Abar(j, k)*C(k);
		    eval(nb+14) -= B(i)*C(i)*Abar(i, j)*C(j)*Abar(j, k)*C(k)*C(k);
		    eval(nb+15) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*pow(C(k), 3);
		    eval(nb+18) -= B(i)*Abar(i, j)*pow(C(j), 3)*Abar(j, k)*C(k);
		    eval(nb+20) -= B(i)*Abar(i, j)*C(j)*Abar(j, k)*pow(C(k), 3);
		    eval(nb+21) -= B(i)*Abar(i, j)*C(j)*C(j)*Abar(j, k)*C(k)*C(k);
		    eval(nb+23) -= B(i)*Abar(i, j)*Abar(j, k)*pow(C(k), 4);
		    for (int l = 0; l < k; l++)
		      {
			eval(nb+16) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*Abar(k, l)*C(l);
			eval(nb+22) -= B(i)*Abar(i, j)*C(j)*Abar(j, k)*Abar(k, l)*C(l);
			eval(nb+24) -= B(i)*Abar(i, j)*Abar(j, k)*C(k)*Abar(k, l)*C(l);
			eval(nb+25) -= B(i)*Abar(i, j)*Abar(j, k)*Abar(k, l)*C(l)*C(l);
		      }
		    
		    for (int l = 0; l < j; l++)
		      eval(nb+19) -= B(i)*Abar(i, j)*Abar(j, k)*Abar(j, l)*C(k)*C(l);
		  }
	      }
	  }
	
	nb += 26;
      }

    if (order >= 10)
      {
	eval(nb) = one/10; eval(nb+1) = one/60; eval(nb+2) = one/120; eval(nb+3) = one/360;
	eval(nb+4) = one/200; eval(nb+5) = one/1200; eval(nb+6) = one/720; eval(nb+7) = one/300;
	eval(nb+8) = one/1800; eval(nb+9) = one/3600; eval(nb+10) = one/2160;
	eval(nb+11) = one/1200; eval(nb+12) = one/7200; eval(nb+13) = one/1440;
	eval(nb+14) = one/420; eval(nb+15) = one/2520; eval(nb+16) = one/5040;
	eval(nb+17) = one/8400; eval(nb+18) = one/50400; eval(nb+19) = one/1800;
	eval(nb+20) = one/10800; eval(nb+21) = one/21600; eval(nb+22) = one/2400;
	eval(nb+23) = one/14400; eval(nb+24) = one/560; eval(nb+25) = one/3360;
	eval(nb+26) = one/20160; eval(nb+27) = one/11200; eval(nb+28) = one/67200;
	eval(nb+29) = one/16800; eval(nb+30) = one/100800; eval(nb+31) = one/6720;
	eval(nb+32) = one/201600; eval(nb+33) = one/720; eval(nb+34) = one/4320;
	eval(nb+35) = one/8640; eval(nb+36) = one/25920; eval(nb+37) = one/14400;
	eval(nb+38) = one/86400; eval(nb+39) = one/51840; eval(nb+40) = one/21600;
	eval(nb+41) = one/129600; eval(nb+42) = one/259200; eval(nb+43) = one/30240;
	eval(nb+44) = one/181440; eval(nb+45) = one/362880; eval(nb+46) = one/604800;
	eval(nb+47) = one/3628800;
	for (int i = 0; i < n; i++)
	  {
	    eval(nb) -= B(i)*pow(C(i), 9);
	    for (int j = 0; j < i; j++)
	      {
		eval(nb+1) -= B(i)*pow(C(i), 6)*Abar(i, j)*C(j);
		eval(nb+2) -= B(i)*pow(C(i), 5)*Abar(i, j)*C(j)*C(j);
		eval(nb+4) -= B(i)*pow(C(i), 4)*Abar(i, j)*pow(C(j), 3);
		eval(nb+7) -= B(i)*pow(C(i), 3)*Abar(i, j)*pow(C(j), 4);
		eval(nb+14) -= B(i)*C(i)*C(i)*Abar(i, j)*pow(C(j), 5);
		eval(nb+24) -= B(i)*C(i)*Abar(i, j)*pow(C(j), 6);
		eval(nb+33) -= B(i)*Abar(i, j)*pow(C(j), 7);
		for (int k = 0; k < i; k++)
		  {
		    eval(nb+3) -= B(i)*pow(C(i), 3)*Abar(i, j)*Abar(i, k)*C(j)*C(k);
		    eval(nb+6) -= B(i)*C(i)*C(i)*Abar(i, j)*Abar(i, k)*C(j)*C(j)*C(k);
		    eval(nb+11) -= B(i)*C(i)*Abar(i, j)*Abar(i, k)*pow(C(j), 3)*C(k);
		    eval(nb+13) -= B(i)*C(i)*Abar(i, j)*Abar(i, k)*C(j)*C(j)*C(k)*C(k);
		    eval(nb+19) -= B(i)*Abar(i, j)*Abar(i, k)*pow(C(j), 4)*C(k);
		    eval(nb+22) -= B(i)*Abar(i, j)*Abar(i, k)*pow(C(j), 3)*C(k)*C(k);
		    for (int l = 0; l < i; l++)
		      eval(nb+10) -= B(i)*Abar(i, j)*Abar(i, k)*Abar(i, l)*C(j)*C(k)*C(l);
		    
		    for (int l = 0; l < j; l++)
		      {
			eval(nb+12) -= B(i)*C(i)*Abar(i, j)*Abar(i, k)*C(k)*Abar(j, l)*C(l);
			eval(nb+20) -= B(i)*Abar(i, j)*Abar(i, k)*C(j)*C(k)*Abar(j, l)*C(l);
			eval(nb+21) -= B(i)*Abar(i, j)*Abar(i, k)*C(k)*Abar(j, l)*C(l)*C(l);
			eval(nb+23) -= B(i)*Abar(i, j)*Abar(i, k)*C(k)*C(k)*Abar(j, l)*C(l);
		      }
		  }
		
		for (int k = 0; k < j; k++)
		  {
		    eval(nb+5) -= B(i)*pow(C(i), 4)*Abar(i, j)*Abar(j, k)*C(k);
		    eval(nb+8) -= B(i)*pow(C(i), 3)*Abar(i, j)*C(j)*Abar(j, k)*C(k);
		    eval(nb+9) -= B(i)*pow(C(i), 3)*Abar(i, j)*Abar(j, k)*C(k)*C(k);
		    eval(nb+15) -= B(i)*C(i)*C(i)*Abar(i, j)*C(j)*C(j)*Abar(j, k)*C(k);
		    eval(nb+16) -= B(i)*C(i)*C(i)*Abar(i, j)*C(j)*Abar(j, k)*C(k)*C(k);
		    eval(nb+17) -= B(i)*C(i)*C(i)*Abar(i, j)*Abar(j, k)*pow(C(k), 3);
		    eval(nb+25) -= B(i)*C(i)*Abar(i, j)*pow(C(j), 3)*Abar(j, k)*C(k);
		    eval(nb+27) -= B(i)*C(i)*Abar(i, j)*C(j)*Abar(j, k)*pow(C(k), 3);
		    eval(nb+29) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*pow(C(k), 4);
		    eval(nb+31) -= B(i)*C(i)*Abar(i, j)*C(j)*C(j)*Abar(j, k)*C(k)*C(k);
		    eval(nb+34) -= B(i)*Abar(i, j)*pow(C(j), 4)*Abar(j, k)*C(k);
		    eval(nb+35) -= B(i)*Abar(i, j)*pow(C(j), 3)*Abar(j, k)*C(k)*C(k);
		    eval(nb+37) -= B(i)*Abar(i, j)*C(j)*C(j)*Abar(j, k)*pow(C(k), 3);
		    eval(nb+40) -= B(i)*Abar(i, j)*C(j)*Abar(j, k)*pow(C(k), 4);
		    eval(nb+43) -= B(i)*Abar(i, j)*Abar(j, k)*pow(C(k), 5);
		    for (int l = 0; l < k; l++)
		      {
			eval(nb+18) -= B(i)*C(i)*C(i)*Abar(i, j)*Abar(j, k)*Abar(k, l)*C(l);
			eval(nb+28) -= B(i)*C(i)*Abar(i, j)*C(j)*Abar(j, k)*Abar(k, l)*C(l);
			eval(nb+30) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*C(k)*Abar(k, l)*C(l);
			eval(nb+32) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*Abar(k, l)*C(l)*C(l);
			eval(nb+38) -= B(i)*Abar(i, j)*C(j)*C(j)*Abar(j, k)*Abar(k, l)*C(l);
			eval(nb+41) -= B(i)*Abar(i, j)*C(j)*Abar(j, k)*C(k)*Abar(k, l)*C(l);
			eval(nb+42) -= B(i)*Abar(i, j)*C(j)*Abar(j, k)*Abar(k, l)*C(l)*C(l);
			eval(nb+44) -= B(i)*Abar(i, j)*Abar(j, k)*C(k)*C(k)*Abar(k, l)*C(l);
			eval(nb+45) -= B(i)*Abar(i, j)*Abar(j, k)*C(k)*Abar(k, l)*C(l)*C(l);
			eval(nb+46) -= B(i)*Abar(i, j)*Abar(j, k)*Abar(k, l)*pow(C(l), 3);
			for (int m = 0; m < l; m++)
			  eval(nb+47) -= B(i)*Abar(i, j)*Abar(j, k)*Abar(k, l)*Abar(l, m)*C(m);
		      }
		    
		    for (int l = 0; l < j; l++)
		      {
			eval(nb+26) -= B(i)*C(i)*Abar(i, j)*Abar(j, k)*Abar(j, l)*C(k)*C(l);
			eval(nb+36) -= B(i)*Abar(i, j)*C(j)*Abar(j, k)*Abar(j, l)*C(k)*C(l);
			eval(nb+39) -= B(i)*Abar(i, j)*Abar(j, k)*Abar(j, l)*C(k)*C(k)*C(l);
		      }
		  }
	      }
	  }
	
	nb += 48;
      }
  }

  void CheckCondition(const Matrix<Real_wp>& Abar, const VectReal_wp& B,
		      const VectReal_wp& Bbar, const VectReal_wp& C) const
  {
    Real_wp err, one(1);
    cout << "Erreur sur bar{b}_i = b_i (1 - c_i)" << endl;
    for (int i = 0; i < Bbar.GetM(); i++)
      {
	err = abs(Bbar(i) - B(i)*(one-C(i)));
	cout << "i = " << i << ", errB = " << err << endl;
      }
    
    cout << endl;
    cout << "Erreur sur c_i^2/2 = sum a_{i, j}" << endl;
    for (int i = 0; i < Bbar.GetM(); i++)
      {
	Real_wp sumA(0);
	for (int j = 0; j < i; j++)
	  sumA += Abar(i, j);
	
	err = abs(C(i)*C(i)/2 - sumA);
	cout << "i = " << i << ", errA = " << err << endl;
      }    

    cout << endl;
    cout << "Erreur sur c_i^3/6 = sum a_{i, j} c_j" << endl;
    for (int i = 0; i < Bbar.GetM(); i++)
      {
	Real_wp sumA(0);
	for (int j = 0; j < i; j++)
	  sumA += Abar(i, j)*C(j);
	
	err = abs(C(i)*C(i)*C(i)/6 - sumA);
	cout << "i = " << i << ", B(i) = " << B(i) << ", errC = " << err << endl;
      }    

    cout << endl;
    cout << "Erreur sur c_i^4/12 = sum a_{i, j} c_j^2" << endl;
    for (int i = 0; i < Bbar.GetM(); i++)
      {
	Real_wp sumA(0);
	for (int j = 0; j < i; j++)
	  sumA += Abar(i, j)*C(j)*C(j);
	
	err = abs(C(i)*C(i)*C(i)*C(i)/12 - sumA);
	cout << "i = " << i << ", B(i) = " << B(i) << ", errC = " << err << endl;
      }    

    cout << endl;
    cout << "Erreur sur c_i^5/20 = sum a_{i, j} c_j^3" << endl;
    for (int i = 0; i < Bbar.GetM(); i++)
      {
	Real_wp sumA(0);
	for (int j = 0; j < i; j++)
	  sumA += Abar(i, j)*C(j)*C(j)*C(j);
	
	err = abs(C(i)*C(i)*C(i)*C(i)*C(i)/20 - sumA);
	cout << "i = " << i << ", B(i) = " << B(i) << ", errC = " << err << endl;
      }    

    cout << endl;
    cout << "Erreur sur c_i^6/30 = sum a_{i, j} c_j^4" << endl;
    for (int i = 0; i < Bbar.GetM(); i++)
      {
	Real_wp sumA(0);
	for (int j = 0; j < i; j++)
	  sumA += Abar(i, j)*square(C(j)*C(j));
	
	err = abs(pow(C(i), 6)/30 - sumA);
	cout << "i = " << i << ", B(i) = " << B(i) << ", errC = " << err << endl;
      }    

  }
  
};

bool GetLU_success(Matrix<Real_wp>& A, IVect& pivot)
{
#ifdef SELDON_LAPACK_CHECK_INFO
  try
    {
#endif
      GetLU(A, pivot);
#ifdef SELDON_LAPACK_CHECK_INFO
    }
  catch (LapackError)
    {
      return false;
    }
#endif
  
  return true;
}


void GetMatrixDz(const Matrix<Real_wp>& Abar, const VectReal_wp& B,
		 const VectReal_wp& Bbar, const VectReal_wp& C, const Real_wp& z0,
		 TinyMatrix<UnivariatePolynomial<Real_wp>, General, 2, 2>& D)
{
  UnivariatePolynomial<Real_wp> one, z, zero;
  one.SetOrder(0); one(0) = Real_wp(1);
  zero.SetOrder(0); zero(0) = Real_wp(0);
  z.SetOrder(1); z(1) = Real_wp(1); z(0) = Real_wp(0);
  
  D(0, 0) = one; D(0, 1) = z;
  D(1, 0) = zero; D(1, 1) = one;

  Vector<UnivariatePolynomial<Real_wp> > kw(Abar.GetM()), ky(Abar.GetM());
  for (int i = 0; i < Abar.GetM(); i++)
    {
      kw(i) = C(i)*z;
      ky(i) = one;
      for (int j = 0; j < i; j++)
	{
	  kw(i) += Abar(i, j)*z*kw(j);
	  ky(i) += Abar(i, j)*z*ky(j);
	}
      
      D(0, 0) += Bbar(i)*ky(i)*z;
      D(0, 1) += Bbar(i)*kw(i)*z;
      D(1, 0) += B(i)*ky(i);
      D(1, 1) += B(i)*kw(i);
    }

  TinyMatrix<Real_wp, General, 2, 2> Dz;
  Dz(0, 0) = D(0, 0).Evaluate(z0);
  Dz(0, 1) = D(0, 1).Evaluate(z0);
  Dz(1, 0) = D(1, 0).Evaluate(z0);
  Dz(1, 1) = D(1, 1).Evaluate(z0);

  DISP(Dz);
}

Real_wp GetAmplificationRKN(const Matrix<Real_wp>& Abar, const VectReal_wp& B,
			    const VectReal_wp& Bbar, const VectReal_wp& C,
			    const Real_wp& z, Complex_wp& L1, Complex_wp& L2)
{
  Real_wp one(1);
  TinyMatrix<Real_wp, General, 2, 2> D;
  D(0, 0) = one; D(0, 1) = z;
  D(1, 0) = 0; D(1, 1) = one;
  VectReal_wp kw(Abar.GetM()), ky(Abar.GetM());
  for (int i = 0; i < Abar.GetM(); i++)
    {
      kw(i) = C(i)*z*z;
      ky(i) = z;
      for (int j = 0; j < i; j++)
	{
	  kw(i) += Abar(i, j)*z*kw(j);
	  ky(i) += Abar(i, j)*z*ky(j);
	}
      
      D(0, 0) += Bbar(i)*ky(i);
      D(0, 1) += Bbar(i)*kw(i);
      D(1, 0) += one/z*B(i)*ky(i);
      D(1, 1) += one/z*B(i)*kw(i);
    }
  DISP(z); DISP(D);
  
  TinyVector<Real_wp, 2> Lr, Li;
  GetEigenvalues(D, Lr, Li);
  
  L1 = Complex_wp(Lr(0), Li(0));
  L2 = Complex_wp(Lr(1), Li(1));
  Real_wp G = max(abs(L1), abs(L2));
  return G;
}

Real_wp GetAmplificationRKN(const Matrix<Real_wp>& Abar, const VectReal_wp& B,
			    const VectReal_wp& Bbar, const VectReal_wp& C,
			    const Real_wp& z)
{
  Complex_wp L1, L2;
  return GetAmplificationRKN(Abar, B, Bbar, C, z, L1, L2);
}


void PermuteCloseRoots(const Complex_wp& L1_old, const Complex_wp& L2_old,
		       Complex_wp& L1_new, Complex_wp& L2_new)
{
  Complex_wp L1 = L1_new, L2 = L2_new;
  if (abs(L1 - L1_old) < abs(L1 - L2_old))
    {
      L1_new = L1;
      L2_new = L2;
    }
  else
    {
      L1_new = L2;
      L2_new = L1;
    }
}

Real_wp GetDifferenceAngle(const Complex_wp& z1, const Complex_wp& z2)
{
  Real_wp angle1 = arg(z1);
  Real_wp angle2 = arg(z2);
  if (angle2 < angle1-pi_wp)
    angle2 += 2*pi_wp;
  
  if (angle2 > angle1 + pi_wp)
    angle2 -= 2*pi_wp;
  
  return angle2 - angle1;
}


Real_wp GetIntersectionZ(const Complex_wp& pt1, const Complex_wp& pt2,
                         const Real_wp& v1, const Real_wp& v2)
{
  if (v1 == v2)
    return Real_wp(1);
  
  Real_wp angle1 = arg(pt1);
  Real_wp angle2 = arg(pt2);
  Real_wp zmin = 1e300;
  Real_wp zi = (angle2 - angle1) / (v1 - v2);
  if (zi > 0)
    zmin = zi;
  
  zi = (angle2 + 2*pi_wp - angle1) / (v1 - v2);
  if (zi > 0)
    zmin = min(zmin, zi);
  
  zi = (angle2 - 2*pi_wp - angle1) / (v1 - v2);
  if (zi > 0)
    zmin = min(zmin, zi);          
  
  return zmin;
}

Real_wp GetCFL_RungeKuttaNystrom(const Matrix<Real_wp>& Abar, const VectReal_wp& B,
				 const VectReal_wp& Bbar, const VectReal_wp& C, bool display = false)
{
  Real_wp one(1), zero(0);
  Real_wp dz_min(1e-5), dz_max(1);
  Real_wp z0 = -dz_min;
  if (epsilon_machine < 1e-30)
    {
      z0 = -Real_wp(1e-12);
      dz_min = Real_wp(1e-7);
    }
  
  Real_wp target = one + 1e-12;
  Complex_wp L1, L2, L1_new, L2_new;
  Real_wp G0 = GetAmplificationRKN(Abar, B, Bbar, C, z0, L1, L2);
  Real_wp dz_target = 1e-12; int nb_iter_max(50);
  Real_wp dz_loc = 1e-12, target_loc = one;
  
  if (epsilon_machine < 1e-30)
    {
      target = one + 1e-16;  
      dz_loc = 1e-14;
      target_loc = one -1e-18;
    }

  if (G0 > target)
    {
      // instable algorithm
      return zero;
    }
  
  // searching first point with |G| > 1
  Real_wp z1 = z0, G1 = G0, Gprev = G0, zprev = z0;
  int nb_eval = 0;
  while (G1 < target)
    {
      zprev = z0; z0 = z1; Gprev = G0; G0 = G1;
      
      // roots at z-epsilon are evaluated to estimate the velocity angle for each root
      GetAmplificationRKN(Abar, B, Bbar, C, z0-dz_min/2, L1_new, L2_new);
      PermuteCloseRoots(L1, L2, L1_new, L2_new);
      
      Real_wp v1 = GetDifferenceAngle(L1, L1_new) / dz_min;
      Real_wp v2 = GetDifferenceAngle(L2, L2_new) / dz_min;
      
      // setting the new dz by computing the angle where an eventual intersection may occur
      Real_wp dist(0);
      if ((v1 == 0) && (v2 == 0))
	{
	  // case where the angle does not change
	  // the amplitude defines the velocity
	  v1 = (abs(L1_new) - abs(L1)) / dz_min;
	  v2 = (abs(L2_new) - abs(L2)) / dz_min;
	  dist = (abs(L1) - abs(L2)) / (v2 - v1);
	}
      else
	dist = GetIntersectionZ(L1, L2, v1, v2);
      
      //Real_wp dz_prev = dz;
      Real_wp dz = abs(dist)/2; 
      if (dz < dz_min)
        dz = dz_min;
      else if (dz > dz_max)
        dz = dz_max;
      
      //DISP(v1); DISP(v2); DISP(L1); DISP(L2); DISP(dz);
      // evaluating the roots for the next point and amplification factor
      z1 -= dz;      
      G1 = GetAmplificationRKN(Abar, B, Bbar, C, z1, L1, L2);
      //DISP(z1); DISP(G1);
      
      // detecting a local maxima
      if ((G0 > G1) && (G0 > Gprev) && (z0 < -one))
        {
          if (display)
            cout << "local maxima at z = " << z0 << endl;

          // the local maximum is found by dichotomy
          Real_wp z_min = z1, z_max = zprev, z_med = z0, Gmed = G0;
          int nb_iter = 0;
          while ((abs(z_max-z_min) > dz_loc) && (nb_iter < nb_iter_max))
            {
              Real_wp za = 0.5*(z_min+z_med), zb = 0.5*(z_max+z_med);
              Real_wp Ga = GetAmplificationRKN(Abar, B, Bbar, C, za);
              Real_wp Gb = GetAmplificationRKN(Abar, B, Bbar, C, zb);
	      nb_eval += 2;
              if (Gmed > max(Ga, Gb))
                {
                  z_min = za;
                  z_max = zb;
                }
              else if (Ga > max(Gmed, Gb))
                {
                  z_max = z_med;
                  z_med = za;
                  Gmed = Ga;
                }
              else
                {
                  z_min = z_med;
                  z_med = zb;
                  Gmed = Gb;
                }
              
              nb_iter++;
            }
          
          if (display)
            cout << "value at z = " << z_med << " : " << Gmed << endl;
          
          if (Gmed >= target_loc)
            {
              G1 = Gmed;
              z0 = zprev;
              z1 = z_med;
	      if (Gmed < target)
		{
		  if (display)
		    cout << "Value at z = " << z1 << " : " << G1 << " Number evaluations : " << nb_eval << endl;
		  
		  return z_med;
		}
            }
        }
      
      // maximal searched value
      nb_eval += 2;
      if (z1 < -Real_wp(1000))
	return -Real_wp(1000);
    }
  
  if (display)
    cout << "Value at z = " << z1 << " : " << G1 << " Number evaluations : " << nb_eval << endl;

  // then dichotomy method to find solution of |G| = 1
  Real_wp z, G; int nb_iter = 0;
  while ((abs(z1-z0) > dz_target) && (nb_iter < nb_iter_max))
    {
      z = 0.5*(z0 + z1);
      G = GetAmplificationRKN(Abar, B, Bbar, C, z);
      if (G > target)
	z1 = z;
      else
	z0 = z;
      
      nb_iter++;
    }
  
  z = 0.5*(z0 + z1);
  if (abs(z) > 9.61)
    {
      //DISP(z); DISP(G);
      //abort();
    }

  return z;
}

//! class for maximizing the CFL with free parameters
/*!
  This class implements the computation of the coefficients a_ij, b_i, bbar_i and c_i
  with free parameters
*/
template<class T>
class RknMinimizationCFL : public VirtualMinimizedFunction<double>
{
  // pour changer r, appeler SetOrder
  int r;
  int permut;
  int type_init;
  enum {RANDOM, USER};
  Vector<double> x_init;
 
public : 
  //! default constructor
  RknMinimizationCFL()
  {
    r = 5;
    type_init = RANDOM;
    this->n = 2;
    //this->type_algo = VirtualMinimizedFunction<T>::SIMPLEX;
    this->type_algo = VirtualMinimizedFunction<T>::SIMPLEX2_RAND;
    //this->type_algo = VirtualMinimizedFunction<T>::CG;
    //this->type_algo = VirtualMinimizedFunction<T>::CG_FR;
    //this->type_algo = VirtualMinimizedFunction<T>::STEEPEST_DESCENT;
    //this->type_algo = VirtualMinimizedFunction<T>::BFGS;
    //this->type_algo = VirtualMinimizedFunction<T>::BFGS2;

    this->step_size = 0.2;
  }
  
  void SetInitialGuess(const Vector<double>& x)
  {
    type_init = USER;
    x_init = x;
  }

  //! initial guess for free parameters
  void FindInitGuess(Vector<double>& x)
  {
    if (type_init == USER)
      {
	x = x_init;
	return;
      }

    x.Reallocate(this->n);
    //x(0) = 2.0*T(rand())/RAND_MAX; x(1) = 2.0*T(rand())/RAND_MAX;
    for (int i = 0; i < this->n; i++)
      x(i) = double(rand())/RAND_MAX;
    //x(i) = 3.0*T(rand())/RAND_MAX - 1.5;
    
  }
  
  //! computes coefficients B and Bbar that satisfy the relations \sum_i b_i c_i^k = 1/(k+1)
  //! and bbar_i = b_i (1-c_i)
  bool FindCoefB_FromC(const Vector<T>& C, Matrix<T>& powC,
		       Vector<T>& B, Vector<T>& Bbar)
  {
    int n = C.GetM();
    T one(1);
    
    // storing powers of c_j
    powC.Reallocate(C.GetM(), C.GetM());
    powC.Fill(one);
    for (int i = 1; i < n; i++)
      for (int j = 0; j < C.GetM(); j++)
	powC(j, i) = powC(j, i-1)*C(j);
    
    // vandermonde matrix
    // to invert \sum B_i C_i^j = 1/(j+1)
    Matrix<T> Vdm(n, n);
    Vdm.Fill(one);
    for (int i = 1; i < n; i++)
      for (int j = 0; j < n; j++)
	Vdm(i, j) = powC(j, i);
    
    //DISP(Vdm); //Vdm.WriteText("vdm.dat");
    IVect pivot(n);
    bool success_lu = GetLU_success(Vdm, pivot);
    
    if (success_lu)
      {
	for (int i = 0; i < n; i++)
	  B(i) = one/(i+1);
	
	SolveLU(Vdm, pivot, B);
      }
    else
      {
	B.Zero();
	Bbar.Reallocate(n);
	Bbar.Zero();
	return false;
      }
    
    // then Bbar_i = B_i (1 - C_i)
    Bbar.Reallocate(n);
    for (int i = 0; i < n; i++)
      Bbar(i) = B(i)*(one - C(i));
    
    return true;
  }
  
  //! computes coefficients a_ij from b_i and c_i
  bool FindCoefA_FromBC(const Vector<T>& x, const Vector<T>& B, const Vector<T>& C,
			Matrix<T>& powC, Matrix<T>& A)
  {
    T zero(0), one(1);
    int nb_coef = (A.GetM()*(A.GetM()-1))/2;
    Matrix<T> sys(nb_coef, nb_coef);
    sys.Fill(zero);
    
    Vector<T> rhs(nb_coef); rhs.Fill(0);
    int N = A.GetM(); Matrix<int> num(N, N);
    num.Fill(-1);
    int nb = 0;
    for (int i = 0; i < N; i++)
      for (int j = 0; j < i; j++)
	num(i, j) = nb++;
    
    nb = 0;
    int kmin = 0;
    
    if (r==6)
      {
	kmin = 0;
	// equation \sum_{i=2}^6 b_i c_i^l \sum_{j=0}^{i-1} a_{i, j} c_j^k = 1/[(k+l+3)(k+2)(k+1)]
	for (int k = kmin; k <= r-3; k++)
	  {
	    
	    int lmax = -1;
	    lmax = r-3-k;
	    
	    for (int l = 0; l <= lmax; l++)
	      {
		for (int i = 0; i < r-1; i++)
		  for (int j = 0; j < i; j++)
		    sys(nb, num(i, j)) = B(i)*powC(i, l)*powC(j, k);
		
		rhs(nb) = one/((k+l+3)*(k+2)*(k+1));
		nb++;
	      }
	  }

      }
    else
      {
	// r >=7
	kmin = 2;
	// equation \sum a_{i, j} = c_i^2/2, i=1..r-1
	for (int i = 1; i <= r-1; i++)
	  {
	    for (int j = 0; j < i; j++)
	      sys(nb, num(i, j)) = one;
	    
	    rhs(nb) = powC(i, 2)/2;
	    nb++;	
	  }
	
	// equation \sum a_{i, j} c_j = c_i^3/6, i=2..r-1
	for (int i = 2; i <= r-1; i++)
	  {
	    for (int j = 0; j < i; j++)
	      sys(nb, num(i, j)) = C(j);
	    
	    rhs(nb) = powC(i, 3)/6;
	    nb++;	
	  }
	
	if (r == 8)
	  {
	    kmin = 3;
	    // equation \sum a_{i, j} c_j^2 = c_i^4/12, i=3..r-1
	    for (int i = 3; i <= r-1; i++)
	      {
		for (int j = 0; j < i; j++)
		  sys(nb, num(i, j)) = C(j)*C(j);
		
		rhs(nb) = powC(i, 4)/12;
		nb++;	
	      }
	  }
	
	// equation \sum_{i=2}^6 b_i c_i^l \sum_{j=0}^{i-1} a_{i, j} c_j^k = 1/[(k+l+3)(k+2)(k+1)]
	for (int k = kmin; k <= r-3; k++)
	  {
	    int lmax = -1;
	    lmax = r-2-k;
	    if (k >= 4)
	      lmax = r-3-k;
	    
	    for (int l = 0; l <= lmax; l++)
	      {
		for (int i = 2; i <= r-1; i++)
		  for (int j = 0; j < i; j++)
		    sys(nb, num(i, j)) = B(i)*powC(i, l)*powC(j, k);
		
		rhs(nb) = one/((k+l+3)*(k+2)*(k+1));
		nb++;
	      }
	  }
	
	// equation \sum_{i=2}^6 b_i c_i^l a_{i,1} = 0
	for (int l = 0; l <= r-6; l++)
	  {
	    for (int i = 2; i <= r-1; i++)
	      sys(nb, num(i, 1)) = B(i)*powC(i, l);
	    
	    rhs(nb) = zero;
	    nb++;
	  }
	
      }
    
    //DISP(sys);
    //DISP(rhs);
    //DISP(nb);
    if(r == 10)
      {
	// Step 1
	A(1,0) = powC(1,2)/2.0;
	A(2,1) = powC(2,3)/(6*C(1));
	A(2,0) = powC(2,2)/2.0 - A(2,1);

	//A(3,2) = (powC(3,3)/3 - C(1)*powC(3,2))/powC(2,2);
	A(3, 2) = powC(3, 3)*(C(3) - 2*C(1)) / (12*(powC(2, 2) - C(1)*C(2)));
	A(3, 1) = powC(3,3)/(6*C(1)) - 2*A(3,2);
	A(3, 0) = powC(3,2)/2.0 - A(3,1) - A(3,2);

	A(4, 3) = (powC(4,4) - 2*powC(4,3)*C(2))/(12*(powC(3,2)-C(3)*C(2)));
	A(4, 2) = (powC(4,3)/6 - C(3)*A(4,3))/C(2);
	A(4, 1) = 0;
	A(4, 0) = powC(4,2)/2 - A(4,2) - A(4,3);

	// Step 2
	A(5,1) = 0;
	A(6,1) = 0;
	A(10,1) = 0;
	A(7,1) = -(B(3)*A(3,1))/B(7);
	A(8,1) = -(B(2)*A(2,1))/B(8);

	// Step 3
	A(10,9) = (B(9)/B(10))*(powC(9,2)/2 - C(9) + 0.5);

	A(10,8) = (B(8)*(powC(8,3)/6 - C(8)/2 + one/3)) / (C(10)*B(10));
	A(9,8) = (B(8)*(C(8)*C(8)/2 - C(8) + one/2) - B(10)*A(10,8)) / B(9);

	T alpha = B(7)*(powC(7,2)/2 - C(7) + one/2);
	T beta = B(7)*(powC(7,3)/6 - C(7)/2 + one/3);
	T gamma = B(7)*(powC(7,4)/12 - C(7)/3 + one/4);
	T d8 = B(8)*C(8); T d9 = B(9)*C(9); T d10 = B(10)*C(10);
	T e8 = B(8)*powC(8,2); T e9 = B(9)*powC(9,2); T e10 = B(10)*powC(10,2);
	A(8,7) =  -(B(9)*beta*e10 - B(10)*beta*e9 - B(9)*d10*gamma + B(10)*d9*gamma - alpha*d9*e10 + alpha*d10*e9)/(B(8)*d9*e10 - B(8)*d10*e9 - B(9)*d8*e10 + B(9)*d10*e8 + B(10)*d8*e9 - B(10)*d9*e8);
	A(9,7) =   (B(8)*beta*e10 - B(10)*beta*e8 - B(8)*d10*gamma + B(10)*d8*gamma - alpha*d8*e10 + alpha*d10*e8)/(B(8)*d9*e10 - B(8)*d10*e9 - B(9)*d8*e10 + B(9)*d10*e8 + B(10)*d8*e9 - B(10)*d9*e8);
	A(10,7) = -(B(8)*beta*e9 - B(9)*beta*e8 - B(8)*d9*gamma + B(9)*d8*gamma - alpha*d8*e9 + alpha*d9*e8)/(B(8)*d9*e10 - B(8)*d10*e9 - B(9)*d8*e10 + B(9)*d10*e8 + B(10)*d8*e9 - B(10)*d9*e8);
	
	// Step 4 (Index start at 1 in the comment, instead 0 for the array)

	// Computation of X_i^3 = \sum_{j=1}^{i-1} a_{ij} C_j^3, i=2,3,...,s-1
	// for i = 2,3,4,5,10 we have (Xpow3(0) = X_1^3 = 0; Xpow3(1) = X_2^3 ... Xpow3(9) = X_{10}^3, Xpow3(10) = X_{11}^3)
	Vector<T> Xpow3(11); Xpow3.Fill(0);
	for(int i = 1; i < 5; i++)
	  for(int j = 0; j < i; j++)
	    Xpow3(i) += A(i,j)*powC(j,3);

	// From (2.k) \sum_{j=1}^{i-1} a_[ij}C_j^3 = \frac{C_i^5}{20}, i=5,10
	Xpow3(9) = powC(9,5)/20;

	// We solve  \sum_{i=1}^{s-1}b_i C_i^k X_i^3 = \frac{1}{(k+6)*20}, k=0,1,2,3,4 for X_i^3, i=6,7,8,9,11
	Vector<T> rhsXpow3(5); rhsXpow3.Fill(0);
	Matrix<T> matXpow3(5, 5);
	matXpow3.Fill(one);
	for (int k = 0; k < 5; k++)
	  {
	    for (int j = 0; j < 5; j++)
	      if(j==4) 
		matXpow3(k, j) = B(10)*powC(10, k);
	      else 
		matXpow3(k, j) = B(j+5)*powC(j+5, k);
	    
	    rhsXpow3(k) = one/((k+6)*20);
	    for(int i=1; i <=4; i++)
	      rhsXpow3(k) -= B(i)*powC(i,k)*Xpow3(i);

	    rhsXpow3(k) -= B(9)*powC(9,k)*Xpow3(9);
	  }
	
	IVect pivotXpow3;
	GetLU(matXpow3,pivotXpow3);
	SolveLU(matXpow3,pivotXpow3, rhsXpow3);
	for(int j=0; j<5; j++)
	  if(j==4)
	    Xpow3(10) = rhsXpow3(j);
	  else
	    Xpow3(j+5) = rhsXpow3(j);
	 
	// now we compute $a_{6,1}, \, a_{6,3}, \, a_{6,4}$ and $a_{6,5}$
	Vector<T> rhsA5(4); rhsA5.Fill(0);
        Matrix<T> VdmA5(4, 4); VdmA5.Fill(one);
	rhsA5(0) = powC(5,2)/2;
	for(int i=1; i <4; i++)
	  {
	    for(int j=0; j<4; j++)
	      {
		if(j==0)
		  VdmA5(i,j) = zero;
		else
		  VdmA5(i,j) = powC(j+1,i);
	      }

	    if(i==3)
	      rhsA5(i) = Xpow3(5);
	    else
	      rhsA5(i) = powC(5,i+2)/((i+2)*(i+1));
	  }
	
	IVect pivotA5;
        GetLU(VdmA5,pivotA5);
        SolveLU(VdmA5,pivotA5, rhsA5);
	A(5,0) = rhsA5(0);
	A(5,2) = rhsA5(1);
	A(5,3) = rhsA5(2);
	A(5,4) = rhsA5(3);

	// step 5 (Index start at 1 in the comment, instead 0 for the array)
	// Computation of X_i^4 = \sum_{j=1}^{i-1} a_{ij} C_j^4, i=2,3,...,s-1   
        // for i = 2,3,4,5,10 we have (Xpow4(0) = X_1^4 = 0; Xpow4(1) = X_2^4 ... Xpow4(9) = X_{10}^4, Xpow4(10) = X_{11}^4)
	Vector<T> Xpow4(11); Xpow4.Fill(0);
	for(int i = 1; i < 6; i++)
	  for(int j = 0; j < i; j++)
	    Xpow4(i) += A(i,j)*powC(j,4);

	// From (2.m) \sum_{j=1}^{i-1} a_[ij}C_j^4 = 0, i=10
	Xpow4(9) = zero;

	// We solve  \sum_{i=1}^{s-1}b_i C_i^k X_i^4 = \frac{1}{(k+7)*30}, k=0,1,2,3 for X_i^4, i=7,8,9,11
	Vector<T> rhsXpow4(4); rhsXpow4.Fill(0);
	Matrix<T> matXpow4(4, 4); matXpow4.Fill(one);
	for (int k = 0; k < 4; k++)
	  {
	    for (int j = 0; j < 4; j++)
	      if(j==3) 
		matXpow4(k, j) = B(10)*powC(10, k);
	      else 
		matXpow4(k, j) = B(j+6)*powC(j+6, k);
	    
	    rhsXpow4(k) = one/((k+7)*30);
	    for(int i=1; i <=5; i++)
	      rhsXpow4(k) -= B(i)*powC(i,k)*Xpow4(i);

	    rhsXpow4(k) -= B(9)*powC(9,k)*Xpow4(9);
	  }
	
	IVect pivotXpow4;
        GetLU(matXpow4,pivotXpow4);
        SolveLU(matXpow4,pivotXpow4, rhsXpow4);
        for(int j=0; j<4; j++)
          if(j==3)
            Xpow4(10) = rhsXpow4(j);
          else
            Xpow4(j+6) = rhsXpow4(j);

        // now we compute $a_{7,1}, \,_{7,3}, \, a_{7,4}, \, a_{7,5}$ and $a_{7,6}$ 
	Vector<T> rhsA6(5); rhsA6.Fill(0);
        Matrix<T> VdmA6(5, 5); VdmA6.Fill(one);
	rhsA6(0) = powC(6,2)/2;
	
	for(int i=1; i < 5; i++)
	  {
	    for(int j=0; j< 5; j++)
	      {
		if(j==0)
		  VdmA6(i,j) = zero;
		else
		  VdmA6(i,j) = powC(j+1,i);
	      }
	    if(i < 3)
	      rhsA6(i) = powC(6,i+2)/((i+2)*(i+1));
	  }
	rhsA6(3) = Xpow3(6);
	rhsA6(4) = Xpow4(6);

	IVect pivotA6;
        GetLU(VdmA6,pivotA6);
        SolveLU(VdmA6,pivotA6, rhsA6);

	A(6,0) = rhsA6(0);
	A(6,2) = rhsA6(1);
	A(6,3) = rhsA6(2);
	A(6,4) = rhsA6(3);
	A(6,5) = rhsA6(4);

	// step 6 (Index start at 1 in the comment, instead 0 for the array)
	// Computation of X_i^4 = \sum_{j=1}^{i-1} a_{ij} C_j^4, i=2,3,...,s-1   
        // we know Xpow5 for i <= 7  (Xpow5(0) = X_1^4 = 0; Xpow4(1) = X_2^5 ... Xpow4(9) = X_{10}^5, Xpow4(10) = X_{11}^5)
	Vector<T> Xpow5(11); Xpow5.Fill(0);
	for(int i = 1; i < 7; i++)
	  for(int j = 0; j < i; j++)
	    Xpow5(i) += A(i,j)*powC(j,5);

	// We solve  \sum_{i=1}^{s-1}b_i C_i^k X_i^5 = \frac{1}{(k+8)*42}, k=0,1,2,3 for X_i^5, i=8,9,10,11
	Vector<T> rhsXpow5(4); rhsXpow5.Fill(0);
	Matrix<T> matXpow5(4, 4); matXpow5.Fill(one);
	for (int k = 0; k < 4; k++)
	  {
	    for (int i = 0; i < 4; i++) 
	      matXpow5(k, i) = B(i+7)*powC(i+7, k);
	    
	    rhsXpow5(k) = one/((k+8)*42);
	    if (k == 3)
	      rhsXpow5(k) = x(3);
	    
	    for(int i=1; i <=6; i++)
	      rhsXpow5(k) -= B(i)*powC(i,k)*Xpow5(i);
	  }
	
	IVect pivotXpow5;
        GetLU(matXpow5,pivotXpow5);
        SolveLU(matXpow5,pivotXpow5, rhsXpow5);
        for(int j=0; j<4; j++)
          Xpow5(j+7) = rhsXpow5(j);
		
	// now we compute $a_{i,1}, \,_{i,3}, \, a_{i,4}, \, a_{i,5}, \, a_{i,6}$ and $a_{i,7}$ for $i=8,9,10,11$
	for(int i=7; i < 11; i++)
	  {
	    Vector<T> rhsA(6); rhsA.Fill(0);
	    Matrix<T> VdmA(6, 6); VdmA.Fill(one);
	    rhsA(0) = powC(i,2)/2;
	    
	    for(int k=1; k < 6; k++)
	      {
		for(int j=0; j< 6; j++)
		  {
		    if(j==0)
		      VdmA(k,j) = zero;
		    else
		      VdmA(k,j) = powC(j+1,k);
		  }
		
		if(k < 3)
		  rhsA(k) = powC(i,k+2)/((k+2)*(k+1));
	      }
	    
	    rhsA(3) = Xpow3(i);
	    rhsA(4) = Xpow4(i);
	    rhsA(5) = Xpow5(i);

	    if (i == 7)
	      for (int k = 0; k < 6; k++)
		rhsA(k) -= A(7, 1)*powC(1, k);
	    else if (i == 8)
	      for (int k = 0; k < 6; k++)
		rhsA(k) -= A(8, 1)*powC(1, k) + A(8, 7)*powC(7, k);
	    else if (i == 9)
	      for (int k = 0; k < 6; k++)
		rhsA(k) -= A(9, 7)*powC(7, k) + A(9, 8)*powC(8, k);
	    else if (i == 10)
	      for (int k = 0; k < 6; k++)
		rhsA(k) -= A(10, 7)*powC(7, k) + A(10, 8)*powC(8, k) + A(10, 9)*powC(9, k);
	    
	    IVect pivotA;
	    GetLU(VdmA,pivotA);
	    SolveLU(VdmA,pivotA, rhsA);

	    A(i,0) = rhsA(0);
	    A(i,2) = rhsA(1);
	    A(i,3) = rhsA(2);
	    A(i,4) = rhsA(3);
	    A(i,5) = rhsA(4);
	    A(i,6) = rhsA(5);
	  }
	
	// conditions
	/* cout << "Condition 2.5" << endl;
	   for (int k = 0; k < 5; k++)
	   {
	   Real_wp sum(0);
	   for (int i = 0; i < B.GetM(); i++)
	   for (int j = 0; j < i; j++)
	   sum += B(i)*powC(i, k)*A(i, j)*powC(j, 3);
	    
	   Real_wp sum_r = one / (20*(k+6));
	   DISP(sum-sum_r);
	   }

	   cout << "Condition 2.6" << endl;
	   for (int k = 0; k < 4; k++)
	   {
	   Real_wp sum(0);
	   for (int i = 0; i < B.GetM(); i++)
	   for (int j = 0; j < i; j++)
	   sum += B(i)*powC(i, k)*A(i, j)*powC(j, 4);
	    
	   Real_wp sum_r = one / (30*(k+7));
	   DISP(sum-sum_r);
	   }

	   cout << "Condition 2.7" << endl;
	   for (int k = 0; k < 3; k++)
	   {
	   Real_wp sum(0);
	   for (int i = 0; i < B.GetM(); i++)
	   for (int j = 0; j < i; j++)
	   sum += B(i)*powC(i, k)*A(i, j)*powC(j, 5);
	    
	   Real_wp sum_r = one / (42*(k+8));
	   DISP(sum-sum_r);
	   }

	   cout << "Condition 2.8" << endl;
	   for (int i = 1; i < B.GetM(); i++)
	   {
	   Real_wp sum(0);
	   for (int j = 0; j < i; j++)
	   sum += A(i, j);
	    
	   DISP(i+1); DISP(sum - powC(i, 2)/2);
	   }

	   cout << "Condition 2.9" << endl;
	   for (int i = 2; i < B.GetM(); i++)
	   {
	   Real_wp sum(0);
	   for (int j = 0; j < i; j++)
	   sum += A(i, j)*C(j);
	    
	   DISP(i+1); DISP(sum - powC(i, 3)/6);
	   }

	   cout << "Condition 2.10" << endl;
	   for (int i = 2; i < B.GetM(); i++)
	   {
	   Real_wp sum(0);
	   for (int j = 0; j < i; j++)
	   sum += A(i, j)*powC(j, 2);
	    
	   DISP(i+1); DISP(sum - powC(i, 4)/12);
	   }

	   cout << "Condition 2.11" << endl;
	   for (int i = 4; i < 10; i += 5)
	   {
	   Real_wp sum(0);
	   for (int j = 0; j < i; j++)
	   sum += A(i, j)*powC(j, 3);
	    
	   DISP(i+1); DISP(sum - powC(i, 5)/20);
	   }

	   cout << "Condition 2.12" << endl;
	   for (int i = 9; i < 10; i += 5)
	   {
	   Real_wp sum(0);
	   for (int j = 0; j < i; j++)
	   sum += A(i, j)*powC(j, 4);
	    
	   DISP(i+1); DISP(sum);
	   }
	
	   cout << "Condition 2.13" << endl;
	   for (int k = 0; k < 5; k++)
	   {
	   Real_wp sum(0);
	   for (int i = 2; i < B.GetM(); i++)
	   sum += B(i)*powC(i, k)*A(i, 1);
	    
	   DISP(sum);
	   }

	   cout << "Condition 2.14" << endl;
	   for (int j = 7; j <= 9; j++)
	   {
	   Real_wp sum(0);
	   for (int i = j+1; i < B.GetM(); i++)
	   sum += B(i)*A(i, j);
	    
	   Real_wp sum_r = B(j)*(powC(j, 2)/2 - C(j) + one/2);
	   DISP(sum-sum_r);
	   }

	   cout << "Condition 2.15" << endl;
	   for (int j = 7; j <= 8; j++)
	   {
	   Real_wp sum(0);
	   for (int i = j+1; i < B.GetM(); i++)
	   sum += B(i)*C(i)*A(i, j);
	    
	   Real_wp sum_r = B(j)*(powC(j, 3)/6 - C(j)/2 + one/3);
	   DISP(sum-sum_r);
	   }

	   cout << "Condition 2.16" << endl;
	   for (int j = 7; j <= 7; j++)
	   {
	   Real_wp sum(0);
	   for (int i = j+1; i < B.GetM(); i++)
	   sum += B(i)*C(i)*C(i)*A(i, j);
	    
	   Real_wp sum_r = B(j)*(powC(j, 4)/12 - C(j)/3 + one/4);
	   DISP(sum-sum_r);
	   }
	*/
	  
	//exit(0);
      }
    
    if(r != 10)
      {
	IVect pivot;
	bool success = GetLU_success(sys, pivot);
	if (success)
	  {
	    SolveLU(sys, pivot, rhs);
	    
	    nb = 0;
	    for (int i = 0; i < N; i++)
	      for (int j = 0; j < i; j++)
		A(i, j) = rhs(nb++);
	  }
	else
	  {
	    A.Zero();
	    return false;
	  }
      }
    //DISP(num);
    //DISP(rhs);
    //DISP(nb);

    return true;
  }

  //! returns true if the coefficients c_i are valid (i.e. distinct and |c_i| <= 2 )
  bool InvalidCoefC(Vector<T>& C, bool check_distinct = true)
  {
    for (int i = 0; i < C.GetM(); i++)
      if (abs(C(i)) > 2.0)
        return true;
    
    if (check_distinct)
      for (int i = 0; i < C.GetM(); i++)
	for (int j = 0; j < i; j++)
	  if (abs(C(i) - C(j)) <= 1e-8)
	    return true;
    
    return false;
  }
  
  //! initializes the number of free parameters
  void SetOrder(int r_, int p = 0)
  {
    r = r_;
    permut = p;
    switch(r)
      {
      case 2: this->n = 0; break;
      case 3: this->n = 1; break;
      case 4: this->n = 1; break;
      case 5: this->n = 2; break;
      case 6: this->n = 2; break;
      case 7: this->n = 4; break;
      case 8: this->n = 4; break;
      case 10: this->n = 4; break;
      default: cout << "Order not implemented" << endl; abort();
      }
  }
  
  //! computes coefficients a_ij, b_i, bbar_i, c_i from free parameters contained in x
  /*!
    This method returns true if the coefficients are valid (with nice c_i)
  */
  bool GetCoefficients(const Vector<T>& x, Matrix<T>& Abar,
		       Vector<T>& B, Vector<T>& Bbar, Vector<T>& C)
  {
    if (r == 2)
      {
	Abar.Reallocate(1, 1); Bbar.Reallocate(1);
	B.Reallocate(1); C.Reallocate(1);
	Abar.Fill(0);
	C(0) = Real_wp(0.5);
	B(0) = Real_wp(1); Bbar(0) = Real_wp(0.5);
      }
    else if (r == 3)
      {
	Abar.Reallocate(2, 2); Bbar.Reallocate(2);
	B.Reallocate(2); C.Reallocate(2);
	Abar.Fill(0);
	
	T one(1), alpha(x(0));
	C(0) = alpha; C(1) = (2.0-3.0*alpha) / (3.0*(one-2*alpha));
	B(0) = (one/2*C(1) - one/3)/(C(0)*(C(1) - C(0))); B(1) = one - B(0);
	Bbar(0) = (one/2*C(1) - one/6)/(C(1) - C(0)); Bbar(1) = one/2 - Bbar(0);
	Abar(1, 0) = one/(6*B(1));
      }
    else if (r == 4)
      {
	Abar.Reallocate(3, 3); Bbar.Reallocate(3);
	B.Reallocate(3); C.Reallocate(3);
	Abar.Fill(0);
	
	T one(1), alpha(x(0));
	C(0) = alpha; C(1) = one/2; C(2) = one-alpha;
	
	Abar(1, 0) = (one - 4*alpha)*(one - 2*alpha) / (Real_wp(8)*(6*alpha*(alpha-one)+one));
	Abar(2, 0) = 2*alpha*(one-2*alpha);
	Abar(2, 1) = (one-2*alpha)*(one-4*alpha)/2;
	
	B(0) = one/(Real_wp(6)*square(one-2*alpha));
	B(1) = Real_wp(4)*(6*alpha*(alpha-one) + one)/(Real_wp(6)*square(one-2*alpha));
	B(2) = B(0);
	
	Bbar(0) = B(0)*(one-C(0));
	Bbar(1) = B(1)*(one-C(1));
	Bbar(2) = B(2)*(one-C(2));	
      }
    else if (r == 5)
      {
	T one(1), zero(0);
	Abar.Reallocate(4, 4); Bbar.Reallocate(4);
	B.Reallocate(4); C.Reallocate(4);
	Abar.Fill(0);
	
	C(1) = x(0); C(3) = x(1);
	C(0) = zero;
	C(2) = (12 - 15*(C(1)+C(3)) + 20*C(1)*C(3)) / (15-20*(C(1)+C(3))+30*C(1)*C(3));
	B(1) = (3-4*(C(2)+C(3))+6*C(2)*C(3)) / (12*C(1)*(C(1)-C(2))*(C(1)-C(3)));
	B(2) = (3-4*(C(1)+C(3))+6*C(1)*C(3)) / (12*C(2)*(C(2)-C(1))*(C(2)-C(3)));
	B(3) = (3-4*(C(1)+C(2))+6*C(1)*C(2)) / (12*C(3)*(C(3)-C(1))*(C(3)-C(2)));
	B(0) = one - B(1) - B(2) - B(3);
	
	Bbar(0) = B(0)*(one-C(0));
	Bbar(1) = B(1)*(one-C(1));
	Bbar(2) = B(2)*(one-C(2));
	Bbar(3) = B(3)*(one-C(3));
	
	Abar(1, 0) = C(1)*C(1)/2;
	Abar(2, 1) = (Real_wp(4) - 5*C(3)) / (Real_wp(120)*B(2)*C(1)*(C(2)-C(3)));
	Abar(2, 0) = C(2)*C(2)/2 - Abar(2, 1);
	Abar(3, 1) = (C(2)*(6 - 5*C(2)) + 5*C(1)*C(3) - 4*C(1) - 2*C(3))
	  / (Real_wp(120)*B(3)*C(1)*(C(3)-C(2))*(C(2)-C(1)));
	Abar(3, 2) = (Real_wp(2) - 5*C(1)) / (Real_wp(120)*B(3)*C(2)*(C(2)-C(1)));
	Abar(3, 0) = C(3)*C(3)/2 - Abar(3, 1) - Abar(3, 2);
	
	/*Real_wp alpha = x(0);
	  C(0) = alpha;
	  C(1) = (Real_wp(4) - sqrt(Real_wp(6)))/10;
	  C(2) = (Real_wp(4) + sqrt(Real_wp(6)))/10;
	  C(3) = one;
	
	  B(0) = zero;
	  B(1) = (Real_wp(16) - sqrt(Real_wp(6)))/36;
	  B(2) = (Real_wp(16) + sqrt(Real_wp(6)))/36;
	  B(3) = one/9;
	
	  Bbar(0) = B(0)*(one-C(0));
	  œBbar(1) = B(1)*(one-C(1));
	  Bbar(2) = B(2)*(one-C(2));
	  Bbar(3) = B(3)*(one-C(3));
	
	  Abar(1, 0) = C(1)*C(1)/2;
	  Abar(2, 1) = (one-5*C(0))/(120*B(2)*(C(1)-C(0))*(one-C(2)));
	  Abar(2, 0) = C(2)*C(2)/2 - Abar(2, 1);
	  Abar(3, 1) = Real_wp(3)/Real_wp(40)
	  *(20-25*C(1)-(one-5*C(0))/((C(1)-C(0))*(one-C(2))));
	
	  Abar(3, 2) = 15*C(1)/8;
	  Abar(3, 0) = one/2 - Abar(3, 1) - Abar(3, 2); */
      }
    else if (r == 6)
      {
	/*T one(1), zero(0);
	  Abar.Reallocate(5, 5); Bbar.Reallocate(5);
	  B.Reallocate(5); C.Reallocate(5);
	  Abar.Fill(zero);
	
	  T alpha = x(0);
	  C(0) = zero; C(1) = alpha; C(2) = one/2; C(3) = one-alpha; C(4) = one;
	  B(0) = -(10*alpha*(alpha-one)+one) / (60*alpha*(one-alpha));
	  B(1) = one / (60*alpha*(one-alpha)*square(one-2*alpha));
	  B(2) = 8*(5*alpha*(alpha-one) + one) / (15*square(one-2*alpha));
	  B(3) = B(1); B(4) = B(0);
	
	  Bbar(0) = B(0)*(one-C(0));
	  Bbar(1) = B(1)*(one-C(1));
	  Bbar(2) = B(2)*(one-C(2));
	  Bbar(3) = B(3)*(one-C(3));
	  Bbar(4) = B(4)*(one-C(4));
	
	  Abar(1, 0) = alpha*alpha/2;
	  Abar(2, 0) = (5*alpha-one)*(6*alpha*(alpha-one)+one) / (48*alpha*(5*alpha*(alpha-one)+one));
	  Abar(2, 1) = (one-3*alpha)*(one-2*alpha) / (48*alpha*(5*alpha*(alpha-one)+one));
	  Abar(3, 0) = (one-alpha)/(6*alpha)*(one+alpha*(-9+alpha*(29-24*alpha)));
	  Abar(3, 1) = -(one-alpha)/(6*alpha)*(one-2*alpha)*(one-6*alpha);
	  Abar(3, 2) = 2*(one-alpha)*(one-2*alpha)*(one-3*alpha)/3;
	
	  Abar(4, 0) = (one-3*alpha)*(one-4*alpha)/(6*alpha*(10*alpha*(alpha-one)+one));
	  Abar(4, 1) = (one-alpha)*(-one+alpha*(9+alpha*(-29+24*alpha)))
	  / (6*alpha*(10*alpha*(alpha-one)+one)*square(one-2*alpha));
	  Abar(4, 2) = 2*(one-alpha)*(one-5*alpha)*(6*alpha*(alpha-one)+one)
	  / (3*(10*alpha*(alpha-one)+one)*square(one-2*alpha));

	  Abar(4, 3) = -alpha*alpha / (2*(10*alpha*(alpha-one)+one)*square(one-2*alpha));*/

	
	T one(1), zero(0);
	Abar.Reallocate(5, 5); Bbar.Reallocate(5);
	B.Reallocate(5); C.Reallocate(5);
	Abar.Fill(zero);
	Matrix<T> powC;

	C(1) = x(0); C(2) = x(1);
	C(0) = zero; C(4) = one;
	C(3) = (one/30 - (one/20)*(C(1)+C(2)) + (one/12)*C(1)*C(2)) / (one/20 - (one/12)*(C(1)+C(2)) + (one/6)*C(1)*C(2));

	bool success = FindCoefB_FromC(C, powC, B, Bbar);
        if (!success)
	  return false;

	//DISP(1111);
        success = FindCoefA_FromBC(x, B, C, powC, Abar);
        if (!success)
	  return false;
      }
    else if (r == 7)
      {
	T one(1), zero(0);
	Abar.Reallocate(7, 7); Bbar.Reallocate(7);
	B.Reallocate(7); C.Reallocate(7);
	Abar.Fill(zero);
	Matrix<T> powC;
	
	// symmetric case
	//C(0) = zero; C(1) = 0.5; C(2) = x(0); C(3) = x(1); C(4) = one-x(1);
	
	// non-symmetric case
	C(0) = zero; C(1) = x(0); C(2) = x(1); C(3) = x(2); C(4) = x(3);
	
	// C(5) found by solving
	// 1/7 - 1/6 s1 + 1/5 s2 - 1/4 s3 + 1/3 s4 -1/2 s5
	T sum_ci = C(2) + C(3) + C(4);
	T sum_cicj = C(2)*C(3) + C(4)*(C(2) + C(3));
	T prod_ci = C(2)*C(3)*C(4);
	
	T s1_cte = sum_ci + one; T s1_var = one;
	T s2_cte = sum_ci + sum_cicj; T s2_var = s1_cte;
	T s3_cte = sum_cicj + prod_ci; T s3_var = s2_cte;
	T s4_cte = prod_ci; T s4_var = s3_cte; T s5_var = s4_cte;
	
	C(5) = (-one/7 + s1_cte/6 - s2_cte/5 + s3_cte/4 - s4_cte/3)
	  / (-s1_var/6 + s2_var/5 - s3_var/4 + s4_var/3 - s5_var/2);
	
	C(6) = one;

	bool success = FindCoefB_FromC(C, powC, B, Bbar);
        if (!success)
	  return false;

	//DISP(1111);
        success = FindCoefA_FromBC(x, B, C, powC, Abar);
        if (!success)
	  return false;
      }
    else if (r == 8)
      {
	T one(1), zero(0);
	Abar.Reallocate(8, 8); Bbar.Reallocate(8);
	B.Reallocate(8); C.Reallocate(8);
	Abar.Fill(zero);
	
	Matrix<T> powC;
	
	// non-symmetric case
	C(0) = zero; C(1) = x(0)/2; C(2) = x(0); C(3) = x(1); C(4) = x(2); C(5) = x(3);
	
	// C(6) found by solving
	// 1/8 - 1/7 s1 + 1/6 s2 - 1/5 s3 + 1/4 s4 -1/3 s5 + 1/2 s6 = 0
	T sum_ci = C(2) + C(3) + C(4) + C(5);
	T sum_cicj = C(2)*(C(3) + C(4) + C(5)) + C(3)*(C(4) + C(5)) + C(4)*C(5);
	T sum_cicjck = C(2)*(C(3)*(C(4)+C(5)) + C(4)*C(5)) + C(3)*C(4)*C(5);
	T prod_ci = C(2)*C(3)*C(4)*C(5);
	
	T s1_cte = sum_ci + one; T s1_var = one;
	T s2_cte = sum_ci + sum_cicj; T s2_var = s1_cte;
	T s3_cte = sum_cicj + sum_cicjck; T s3_var = s2_cte;
	T s4_cte = sum_cicjck + prod_ci; T s4_var = s3_cte;
	T s5_cte = prod_ci; T s5_var = s4_cte; T s6_var = s5_cte;
	
	C(6) = (-one/8 + s1_cte/7 - s2_cte/6 + s3_cte/5 - s4_cte/4 + s5_cte/3)
	  / (-s1_var/7 + s2_var/6 - s3_var/5 + s4_var/4 - s5_var/3 + s6_var/2);
	
	C(7) = one;
        
        if (InvalidCoefC(C))
          return false;

	bool success = FindCoefB_FromC(C, powC, B, Bbar);
        if (!success)
	  return false;

	//DISP(1111);
        success = FindCoefA_FromBC(x, B, C, powC, Abar);
        if (!success)
	  return false;
      }
    else if (r == 10)
      {
	T one(1), zero(0);
	Abar.Reallocate(11, 11); Bbar.Reallocate(11);
	B.Reallocate(11); C.Reallocate(11);
	Abar.Fill(zero); Bbar.Fill(zero);
	B.Fill(zero); C.Fill(zero);
	
	Matrix<T> powC;
	
	// must be permutation of the node of the nodes of Lobatto-quadrature
	Vector<T> gamma_s(4); gamma_s.Fill(zero);
	gamma_s(0) = 0.5 * (one - sqrt( (Real_wp(7)+2*sqrt(Real_wp(7)))/21 ));
	gamma_s(1) = 0.5 * (one - sqrt( (Real_wp(7)-2*sqrt(Real_wp(7)))/21 ));
	gamma_s(2) = one - gamma_s(1);
	gamma_s(3) = one - gamma_s(0);
	//for(int i = 0; i < 4; i++)
	//C(i+3) = gamma_s(i);
	
	int num_perm = permut;
	Vector<int> indices(4); indices.Fill(0);
	Matrix<int> permIndices(24, 4); permIndices.Fill(0);
	for (int i=0; i<4; i++)
	  indices(i) = i;
	
	//generate the 24 pemutations
	genPerm(indices,4,permIndices);
	for(int i = 0; i < 4; i++)
	  C(3+i) = gamma_s(permIndices(num_perm,i)); 

	//DISP(num_perm);
	//DISP(indices);
	//DISP(permIndices);

	// Hairer permutation
	//C(3) = gamma_s(3); C(4) = gamma_s(1); C(5) = gamma_s(0); C(6) = gamma_s(2);
	
	// (2.18) gives c_2 = c_4(3 c_4 - 5 c_3) / (5 c_4 - 10 c_3) 
	C(2) = C(4)*(3*C(4) - 5*C(3)) / (5*C(4) - 10*C(3));

	// (2.1) gives 
	C(0) = zero; C(7) = C(3); C(8) = C(2); C(9) = zero; C(10) = one;
	
	// (2.17) gives c_1 = c_2/2
	C(1) = C(2)/2;
	
        if (InvalidCoefC(C, false))
          return false;
        
	powC.Reallocate(C.GetM(), C.GetM());
	powC.Fill(one);
	for (int i = 1; i < C.GetM(); i++)
	  for (int j = 0; j < C.GetM(); j++)
	    powC(j, i) = powC(j, i-1)*C(j);
	
	//B_0, B_2, B_3 are free parameters
	B(0) = x(0); B(2) = x(1); B(3) = x(2);
	
	//(2.2) gives
	B(1) = zero; 
	B(8) = -B(2);
	
	// (2.20) gives
	//Globatto<Real_wp> lob;
	//lob.ConstructQuadrature(5, lob.QUADRATURE_LOBATTO);
	//DISP(lob.Points()); DISP(lob.Weights());
	
	LegendrePolynomial<Real_wp> Pn(6);
	VectReal_wp EvalPn;
	for(int i = 4; i <= 6; i++)
	  {
	    Pn.EvaluatePn(6, 2*C(i)-1.0, EvalPn);
	    T p5ci = EvalPn(5);
	    B(i) = one / (30*p5ci*p5ci);
	  }
	
	Pn.EvaluatePn(6, 2.0*C(3)-1.0, EvalPn);
	T p5c3 = EvalPn(5);
	B(7) = one/(30*p5c3*p5c3) - B(3);
	
	B(9) = one/30 - B(0);

	B(10) = one/30;
	
	// we compute Bbar
	for (int i = 0; i < B.GetM(); i++)
	  Bbar(i) = B(i)*(Real_wp(1) - C(i));
	
	//DISP(B); DISP(Bbar);
        bool success = FindCoefA_FromBC(x, B, C, powC, Abar); 
        if (!success)
          return false;
	//DISP(Abar);
      }

    return true;
  }

  /* function to swap array elements */

  void swap (Vector<int>& v, int i, int j) 
  {
    int tmp;
    int nb = v.GetM();
    if((i <0 || i > nb) || (j < 0 || j > nb))
      {
	cout << "Index out of bounds in swap function." << endl;
	abort();
      }
    tmp = v(i);
    v(i) = v(j);
    v(j) = tmp;
  }

  /*
   *This function generate the 24 permutation possible for an array of four elements                        
   */
  void genPerm(Vector<int>& v, int nb_elem, Matrix<int>& matPerm)
  {
    int i,j, counter=0;
    int fix_index=0;
    //
    for(fix_index = 0; fix_index < nb_elem; fix_index++)
      {
	//first permutation for the fixed index 
	swap(v,0,fix_index);

	for(i=1; i < nb_elem; i++)
	  {
	    // second fixed index
	    swap(v,1,i);
	    for(j=0; j< nb_elem; j++)
	      matPerm(counter,j) = v(j);
	    counter ++;

	    swap(v,2,3);
	    for(j=0; j< nb_elem; j++)
	      matPerm(counter,j) = v(j);
	    counter ++;

	    //swap back values
	    swap(v,2,3);
	    swap(v,1,i);
	  }
	//swap back firs permutation
	swap(v,0,fix_index);
      }
  }

  //! computes feval = -cfl from free parameters contained in x
  void EvaluateFunction(const Vector<double>& x, double& feval)
  {
    Vector<T> xb(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      xb(i) = x(i);
    
    T fevalb;
    Matrix<T> Abar; Vector<T> B, Bbar, C;
    bool success = GetCoefficients(xb, Abar, B, Bbar, C);
    if (success)
      fevalb = GetCFL_RungeKuttaNystrom(Abar, B, Bbar, C);
    else
      fevalb = T(0);

    feval = toDouble(fevalb);
  }
  
  //! computes the gradient of EvaluateFunction
  void EvaluateFunctionGradient(const Vector<double>& xb, 
				double& feval_d, Vector<double>& fjac)
  {
    Vector<T> x(xb.GetM());
    for (int i = 0; i < xb.GetM(); i++)
      x(i) = xb(i);

    Matrix<T> Abar; Vector<T> B, Bbar, C;
    bool success = GetCoefficients(x, Abar, B, Bbar, C);
    if (!success)
      {
        feval_d = double(0);
        fjac.Zero();
        return;
      }
    
    T feval = GetCFL_RungeKuttaNystrom(Abar, B, Bbar, C);
    feval_d = toDouble(feval);
    
    Vector<T> xp(x), xm(x);
    T h = pow(epsilon_machine, T(1)/3), fm, fp;
    for (int i = 0; i < this->n; i++)
      {
	xp(i) = x(i) + h;
	xm(i) = x(i) - h;

	GetCoefficients(xp, Abar, B, Bbar, C);
	fp = GetCFL_RungeKuttaNystrom(Abar, B, Bbar, C);

	GetCoefficients(xm, Abar, B, Bbar, C);
	fm = GetCFL_RungeKuttaNystrom(Abar, B, Bbar, C);
	
	fjac(i) = toDouble((fp - fm)/(2.0*h));
	
	xp(i) = x(i); xm(i) = x(i);
      }
  }

  //! writes the coefficients in a file coefs.dat in a C++ manner
  void WriteCoefficients(Matrix<T>& Abar,
			 Vector<T>& B, Vector<T>& Bbar, Vector<T>& C)
  {
    ofstream file_out("coefs.dat");
    file_out.precision(35);
    file_out.setf(ios::scientific);
    
    for (int i = 0; i < C.GetM(); i++)
      file_out << "C(" << i << ") = to_num<Real_wp>(\"" << C(i) << "\");" << '\n';
    
    file_out << '\n';
    for (int i = 0; i < B.GetM(); i++)
      for (int j = 0; j < i; j++)
	file_out << "Abar(" << i << ", " << j << ") = to_num<Real_wp>(\"" << Abar(i, j) << "\");" << '\n';
    
    file_out << '\n';
    for (int i = 0; i < Bbar.GetM(); i++)
      file_out << "Bbar(" << i << ") = to_num<Real_wp>(\"" << Bbar(i) << "\");" << '\n';
    
    file_out << '\n';
    for (int i = 0; i < B.GetM(); i++)
      file_out << "B(" << i << ") = to_num<Real_wp>(\"" << B(i) << "\");" << '\n';
    
    file_out.close();
  }
  
};

class SequenceParameter
{
public:
  virtual int GetNbPoints() const = 0;
  virtual void SetParameters(int k, VectReal_wp& x) = 0;
  
};

void ComputeCFLParam(RknMinimizationCFL<Real_wp>& fct,
                     SequenceParameter& seq, Vector<double>& CflP)
{
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#else
  int nb_proc(1), rank_proc(0);
#endif
  
  // points are distributed to the different processors
  int N = seq.GetNbPoints(), nb_pts_proc, offset_proc;
  GetParallelDistributionPoints(nb_proc, rank_proc, N, nb_pts_proc, offset_proc);
  
  // in the root processor, we retrieve the number of points
  // computed on each proc and offset
  Vector<int64_t> xtmp;
  IVect NbPointsPerProc, OffsetPointsProc;
  // nmax is the maximum number of points (among all proc)
  int nmax = N;
  NbPointsPerProc.Reallocate(nb_proc);
  OffsetPointsProc.Reallocate(nb_proc);
  if (nb_proc > 1)
    {
#ifdef SELDON_WITH_MPI
      MPI_Gather(&nb_pts_proc, 1, MPI_INTEGER,
                 NbPointsPerProc.GetData(), 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      
      MPI_Gather(&offset_proc, 1, MPI_INTEGER,
                 OffsetPointsProc.GetData(), 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

      nmax = NbPointsPerProc.GetNormInf();
      
      MPI_Bcast(&nmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#else
      cout << "Compile with MPI" << endl;
      abort();
#endif
    }
  else
    {
      NbPointsPerProc(0) = N;
      OffsetPointsProc(0) = 0;
    }
    
  Matrix<Real_wp> Abar;
  Vector<Real_wp> B, Bbar, C;

  // main loop over points
  Vector<double> Cfl(nmax); Cfl.Zero();
  VectReal_wp xmin(fct.GetM()); xmin.Zero();
  for (int k = 0; k < N; k++)
    {
      if ((k >= offset_proc) && (k < offset_proc+nb_pts_proc))
        {
          seq.SetParameters(k, xmin);
          // we check that the parameters are all different
          bool coef_ok = true;
          for (int i = 0; i < xmin.GetM(); i++)
            for (int j = 0; j < xmin.GetM(); j++)
              if ((i != j) && (abs(xmin(i) - xmin(j)) <= epsilon_machine))
                coef_ok = false;
          
          if (coef_ok)
            coef_ok = fct.GetCoefficients(xmin, Abar, B, Bbar, C);
          
          if (coef_ok)
            {
              Real_wp cfl = GetCFL_RungeKuttaNystrom(Abar, B, Bbar, C);
              Cfl(k-offset_proc) = toDouble(sqrt(abs(cfl)));
            }
        }
    }
  
  Vector<double> AllCfl;
  if (rank_proc == 0)
    AllCfl.Reallocate(nmax*nb_proc);

#ifdef SELDON_WITH_MPI
  MpiGather(MPI_COMM_WORLD, Cfl, xtmp, AllCfl, Cfl.GetM(), 0);
#endif
  
  if (rank_proc == 0)
    {
      CflP.Reallocate(N); CflP.Zero();
      for (int p = 0; p < nb_proc; p++)
        for (int k = 0; k < NbPointsPerProc(p); k++)
          {
            int num = OffsetPointsProc(p) + k;
            CflP(num)  = AllCfl(nmax*p+k);
          }      
    }

#ifdef SELDON_WITH_MPI  
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}


class SequenceGrid4D : public SequenceParameter
{
protected:
  int nx, ny, nz, nt;
  VectReal_wp param_x, param_y, param_z, param_t;
  
public :
  SequenceGrid4D(int nx_, int ny_, int nz_, int nt_,
		 const VectReal_wp& x, const VectReal_wp& y, const VectReal_wp& z, const VectReal_wp& t)
    : nx(nx_), ny(ny_), nz(nz_), nt(nt_), param_x(x), param_y(y), param_z(z), param_t(t)
  {
  }

  int GetNbPoints() const { return nx*ny*nz*nt; }

  void SetParameters(int num, VectReal_wp& x)
  {
    int nzt = nz*nt, nyzt = ny*nzt;
    int i = num/nyzt, j = (num-i*nyzt)/nzt, k = (num-i*nyzt-j*nzt)/nt, l = num%nt;
    x(0) = param_x(i); x(1) = param_y(j); x(2) = param_z(k); x(3) = param_t(l);
  }

};


class SequenceGrid3D : public SequenceParameter
{
protected:
  int nx, ny, nz;
  VectReal_wp param_x, param_y, param_z;
  
public :
  SequenceGrid3D(int nx_, int ny_, int nz_, const VectReal_wp& x, const VectReal_wp& y, const VectReal_wp& z)
    : nx(nx_), ny(ny_), nz(nz_), param_x(x), param_y(y), param_z(z)
  {
  }

  int GetNbPoints() const { return nx*ny*nz; }

  void SetParameters(int num, VectReal_wp& x)
  {
    int nyz = ny*nz;
    int i = num/nyz, j = (num-i*nyz)/nz, k = num%nz;
    x(0) = param_x(i); x(1) = param_y(j); x(2) = param_z(k);
  }

};



class SequenceGrid2D : public SequenceParameter
{
protected:
  int nx, ny;
  VectReal_wp param_x, param_y;
  
public :
  SequenceGrid2D(int nx_, int ny_, const VectReal_wp& x, const VectReal_wp& y)
    : nx(nx_), ny(ny_), param_x(x), param_y(y)
  {
  }

  int GetNbPoints() const { return nx*ny; }

  void SetParameters(int k, VectReal_wp& x)
  {
    int i = k/ny, j = k%ny;
    x(0) = param_x(i); x(1) = param_y(j);
  }

};

void ComputeCFLTwoParameter(RknMinimizationCFL<Real_wp>& fct,
                            int nx, int ny, const string& name)
{
  VectReal_wp PhiParam, ThetaParam;
  Linspace(Real_wp(1e-3), Real_wp(0.999), nx, PhiParam);
  Linspace(Real_wp(1e-3), Real_wp(0.999), ny, ThetaParam);

  SequenceGrid2D grid(nx, ny, PhiParam, ThetaParam);
  
  Vector<double> CflAll;
  ComputeCFLParam(fct, grid, CflAll);

#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0);
#endif

  if (rank_proc == 0)
    {
      Matrix<double> CflMat(nx, ny);
      int num = 0;
      for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
          CflMat(i, j) = CflAll(num++);
      
      CflMat.Write(name);
    }
}


class SequenceGrid1D : public SequenceParameter
{
protected:
  int nx;
  VectReal_wp param_x;
  
public :
  SequenceGrid1D(int nx_, const VectReal_wp& x)
    : nx(nx_), param_x(x)
  {
  }

  int GetNbPoints() const { return nx; }

  void SetParameters(int i, VectReal_wp& x)
  {
    //x(0) = 0; x(1) = -0.1; x(2) = 0.2;
    x(0) = param_x(i);
    // x(1) = Real_wp(1)/2;
    // x(1) = Real_wp(1) - x(0);
  }
  
};

void ComputeCFLOneParameter(RknMinimizationCFL<Real_wp>& fct,
                            int nx, const string& name)
{
  VectReal_wp ThetaParam;
  Linspace(Real_wp(1e-6), Real_wp(0.01), nx, ThetaParam);

  SequenceGrid1D grid(nx, ThetaParam);
  
  Vector<double> CflAll;
  ComputeCFLParam(fct, grid, CflAll);

#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0);
#endif
  
  if (rank_proc == 0)
    CflAll.Write(name);
}


void ComputeCFLOptimal(RknMinimizationCFL<Real_wp>& fct,
                       SequenceParameter& seq, int nb_selection)
{
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#else
  int nb_proc(1), rank_proc(0);
#endif

  // points are distributed to the different processors
  int N = seq.GetNbPoints(), nb_pts_proc, offset_proc;
  GetParallelDistributionPoints(nb_proc, rank_proc, N, nb_pts_proc, offset_proc);
  
  Vector<int64_t> xtmp;
  Matrix<Real_wp> Abar;
  Vector<Real_wp> B, Bbar, C;

  // main loop over points to find best candidates
  Vector<double> Cfl(nb_selection); Cfl.Zero();
  Vector<int> NumParam(nb_selection); NumParam.Fill(-1);
  VectReal_wp xmin(fct.GetM()); xmin.Zero();
  DISP(nb_selection); DISP(N);
  for (int k = 0; k < N; k++)
    {
      if ((k >= offset_proc) && (k < offset_proc+nb_pts_proc))
        {
	  seq.SetParameters(k, xmin);
          
          // we check that the parameters are all different
          bool coef_ok = true;
          for (int i = 0; i < xmin.GetM(); i++)
            for (int j = 0; j < xmin.GetM(); j++)
              if ((i != j) && (abs(xmin(i) - xmin(j)) <= epsilon_machine))
                coef_ok = false;
          
          if (coef_ok)
            coef_ok = fct.GetCoefficients(xmin, Abar, B, Bbar, C);
          
          if (coef_ok)
            {
              Real_wp zSquare = GetCFL_RungeKuttaNystrom(Abar, B, Bbar, C);
              Real_wp cfl = toDouble(sqrt(abs(zSquare)));
              // a quel rang se situe cette cfl ?
              int rank = nb_selection;
              for (int p = nb_selection-1; p >= 0; p--)
                if (cfl > Cfl(p))
                  rank = p;

              if (rank < nb_selection)
                {
                  // on insere ce nouveau optimum
                  for (int p = nb_selection-2; p >= rank; p--)
                    {
                      NumParam(p+1) = NumParam(p);
                      Cfl(p+1) = Cfl(p);
                    }
                  
                  NumParam(rank) = k;
                  Cfl(rank) = toDouble(cfl);
                }
            }
        }
      //DISP(N);
    }
  
  DISP(Cfl); DISP(NumParam);
  cout << "Starting optimization" << endl;
  
  // then optimization is performed for each candidate  
  Vector<double> x_opt;
#ifdef MONTJOIE_WITH_GSL
  double fmin(0), feval; 
  for (int k = 0; k < nb_selection; k++)
    {      
      Vector<double> xsol(xmin.GetM()); xsol.Zero();
      
      VectReal_wp xsol_(xmin.GetM()); xsol_.Zero();
      seq.SetParameters(NumParam(k), xsol_);
      
      for (int k = 0; k < xsol.GetM(); k++)
	xsol(k) = xsol_(k);
      
      fct.SetInitialGuess(xsol);
      
      MinimizeParametersGsl(fct, xsol);

      fct.EvaluateFunction(xsol, feval);
      if (feval < fmin)
	{
	  x_opt = xsol;
	  fmin = feval;
	}
    }
#else
  seq.SetParameters(NumParam(0), xmin);
  x_opt.Reallocate(xmin.GetM());
  for (int k = 0; k < xmin.GetM(); k++)
    x_opt(k) = toDouble(xmin(k));
#endif  

  for (int k = 0; k < xmin.GetM(); k++)
    xmin(k) = x_opt(k);

  // parameters and cfl are gathered
  Vector<double> AllCfl;
  Vector<double> AllParam;
  double cfl_opt(0);
  bool coef_ok = fct.GetCoefficients(xmin, Abar, B, Bbar, C);
          
  if (coef_ok)
    {
      Real_wp zSquare = GetCFL_RungeKuttaNystrom(Abar, B, Bbar, C);
      cfl_opt = toDouble(sqrt(abs(zSquare)));
    }
  
  if (rank_proc == 0)
    {
      AllCfl.Reallocate(nb_proc);
      AllParam.Reallocate(nb_proc*xmin.GetM());
    }
  
#ifdef SELDON_WITH_MPI
  MpiGather(MPI_COMM_WORLD, &cfl_opt, xtmp, AllCfl.GetData(), 1, 0);
  MpiGather(MPI_COMM_WORLD, x_opt, xtmp, AllParam, xmin.GetM(), 0);
#else
  cout << "Recompile with MPI" << endl;
  abort();
#endif
  
  if (rank_proc == 0)
    {
      Vector<Vector<double> > all_xmin(nb_proc);
      int num = 0;
      for (int i = 0; i < nb_proc; i++)
        {
          all_xmin(i).Reallocate(xmin.GetM());
          for (int j = 0; j < all_xmin(i).GetM(); j++)
            all_xmin(i)(j) = AllParam(num++);
        }

      Sort(AllCfl, all_xmin);
      for (int k = nb_proc-1; k >= 0; k--)
        {
          cout << "Cfl = " << AllCfl(k) << endl;
          cout << "For parameters = " << all_xmin(k) << endl;
        }
    }
  
#ifdef SELDON_WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
}


void ComputeCFLOptimal4D(RknMinimizationCFL<Real_wp>& fct, int nx, int ny, int nz, int nt)
{
  VectReal_wp PhiParam, ThetaParam, PsiParam, ZetaParam;
  Linspace(Real_wp(-0.4), Real_wp(0.4), nx, PhiParam);
  Linspace(Real_wp(-0.4), Real_wp(0.4), ny, ThetaParam);
  Linspace(Real_wp(-0.4), Real_wp(0.4), nz, PsiParam);
  Linspace(Real_wp(1e-5), Real_wp(0.01), nz, ZetaParam);

  SequenceGrid4D grid(nx, ny, nz, nz, PhiParam, ThetaParam, PsiParam, ZetaParam);
  
  ComputeCFLOptimal(fct, grid, 4);
}

void ComputeCFLOptimal3D(RknMinimizationCFL<Real_wp>& fct, int nx, int ny, int nz)
{
  VectReal_wp PhiParam, ThetaParam, PsiParam;
  Linspace(Real_wp(1e-3), Real_wp(0.999), nx, PhiParam);
  Linspace(Real_wp(1e-3), Real_wp(0.999), ny, ThetaParam);
  Linspace(Real_wp(1e-3), Real_wp(0.999), nz, PsiParam);

  SequenceGrid3D grid(nx, ny, nz, PhiParam, ThetaParam, PsiParam);
  
  ComputeCFLOptimal(fct, grid, 4);
}


void ComputeCFLOptimal2D(RknMinimizationCFL<Real_wp>& fct, int nx, int ny)
{
  VectReal_wp PhiParam, ThetaParam;
  Linspace(Real_wp(1e-3), Real_wp(0.999), nx, PhiParam);
  Linspace(Real_wp(1e-3), Real_wp(0.999), ny, ThetaParam);

  SequenceGrid2D grid(nx, ny, PhiParam, ThetaParam);
  
  ComputeCFLOptimal(fct, grid, 4);  
}

void ComputeCFLOptimal1D(RknMinimizationCFL<Real_wp>& fct, int nx)
{
  VectReal_wp ThetaParam;
  Linspace(Real_wp(1e-6), Real_wp(0.01), nx, ThetaParam);

  SequenceGrid1D grid(nx, ThetaParam);
  
  ComputeCFLOptimal(fct, grid, 4);  
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc != 3)
    {
      cout << "Usage : ./test.x order permut" << endl;
      cout << "Order or permut has not been given" << endl;
      abort();
    }

#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0);
#endif

  int order = atoi(argv[1]);  
  int permut = atoi(argv[2]);
  time_t nr = time(NULL);
  srand(nr);
  
  // class implementing function to minimize
  Vector<Real_wp> xmin; 

  if (rank_proc == 0)
    {
      cout << "Order = " << order << endl;
      cout << "Random seed = " << nr << endl;
    }
  
  DISP(epsilon_machine);
  
  RknMinimizationCFL<Real_wp> fct;
  fct.SetOrder(order, permut);

  Matrix<Real_wp> A0, Abar;
  Vector<Real_wp> B, Bbar, C;
  
  // ici on verifie que la famille est correcte
  // pour un choix particulier de parametres libres
  if (false)
    {
      DISP(fct.GetM());
      VectReal_wp param(fct.GetM());
      //param.FillRand(); param *= 0.2e-9; param(3) *= 0.01;
      //param(0) = 0.2; param(1) = 0.3; param(2) = 0.14;
      
      // ordre 3
      //param(0) = (Real_wp(3) - sqrt(Real_wp(3)))/6;
      
      // ordre 4
      //param(0) = Real_wp(1) / (Real_wp(4)*(Real_wp(1) + cos(pi_wp/9)));
      
      // ordre 5
      //param(0) = 0.2776745179; param(1) = 0.7366565526;
      
      // ordre 6, choix instable
      //param(0) = 0.0816464646464646; param(1) = 0.968757575757576;
      
      // ordre 6, params optimaux
      //param(0) = 0.229183265; param(1) = 0.5;
      
      // ordre 7, choix instable
      //param(0) = 0.105717824550642; param(1) = 0.166238304366948;
      //param(2) = 0.433062752147395; param(3) = 0.576375379285793;

      // ordre 7, param quasi-optimaux
      //param(0) = 0.110451398065702; param(1) = 0.173816271367107;
      //param(2) = 0.459433163929695; param(3) = 0.652002232653235;
      
      // ordre 8, param quasi-optimaux
      //param(0) = 0.135294127286225; param(1) = 0.24015308384744;
      //param(2) = 0.453046953126355; param(3) = 0.695039606659698;
      
      // ordre 10, choix d'Hairer
      //param(0) = 0; param(1) = -0.1; param(2) = 0.0; param(3) = Real_wp(1)/462;

      // ordre 10, param quasi-optimaux
      // permut : 21
      param(0) = -0.0; param(1) = -0.1; param(2) = 0.0; param(3) = 0.0021632268153138;
      
      DISP(param);
      bool coef_ok = fct.GetCoefficients(param, Abar, B, Bbar, C);
      if (!coef_ok)
	{
          cout << "Bad choice of coefficients, try again, fail again" << endl;
          abort();
	}
      
      RknOrderCondition sys;
      sys.SetOrder(order, true, true);
    
      int n = sys.GetM();
      VectReal_wp eval(n);
      sys.EvaluateFunction(Abar, B, Bbar, C, eval);
      //DISP(Abar); DISP(B); DISP(Bbar); DISP(C);
      sys.CheckCondition(Abar, B, Bbar, C);
      
      DISP(eval); DISP(Norm2(eval));
      for (int i = 0; i < eval.GetM(); i++)
        if (abs(eval(i)) > 1e-14)
          {
            DISP(i); DISP(eval(i));
          }
      
      // computing the amplification factor
      VectReal_wp zFac, G;
      VectComplex_wp L1_F, L2_F; Complex_wp L1, L2;
      Linspace(-Real_wp(1)/100, -Real_wp(10), 4001, zFac);
      G.Reallocate(zFac.GetM());
      L1_F.Reallocate(zFac.GetM());
      L2_F.Reallocate(zFac.GetM());
      for (int p = 0; p < zFac.GetM(); p++)
	{
	  G(p) = GetAmplificationRKN(Abar, B, Bbar, C, -zFac(p)*zFac(p), L1, L2);
	  if (p == 0)
	    {
	      L1_F(p) = L1;
	      L2_F(p) = L2;
	    }
	  else
	    {
	      if (abs(L1 - L1_F(p-1)) < abs(L2 - L1_F(p-1)))
		{
		  L1_F(p) = L1;
		  L2_F(p) = L2;
		}
	      else
		{
		  L1_F(p) = L2;
		  L2_F(p) = L1;
		}
	    }
	}
      
      DISP(G.GetNormInf());
      L1_F.Write("L1.dat");
      L2_F.Write("L2.dat");
      G.Write("G.dat");
      //DISP(
      //DISP(GetAmplificationRKN(Abar, B, Bbar, C, Real_wp(0.24)));
      
      Real_wp cfl = GetCFL_RungeKuttaNystrom(Abar, B, Bbar, C, true);
      DISP(sqrt(abs(cfl))); 
      
      return 0;
    }  
      
  // ici cfl en fonction d'un jeu de parametres
  // ComputeCFLTwoParameter(fct, 400, 400, "cfl2D.dat");
  //ComputeCFLOneParameter(fct, 10000, "cfl1D.dat");
  
  if (rank_proc == 0)
    {
      
      cout << "Number of unknowns = " << fct.GetM() << endl;
      xmin.Reallocate(fct.GetM()); xmin.Zero();
      xmin.FillRand(); xmin *= Real_wp(1)/RAND_MAX;
      DISP(xmin);
      switch (order)
        {
        case 3:
	  // cfl = 2.4986
          xmin(0) = (Real_wp(3) - sqrt(Real_wp(3)))/6;
          break;
        case 4:
	  // cfl = 3.9392
          xmin(0) = Real_wp(1) / (Real_wp(4)*(Real_wp(1) + cos(pi_wp/9)));
          break;
        case 5:
	  // cfl = 2.9090
          xmin(0) = 0.2776745179; xmin(1) = 0.7366565526;
          break;
        case 6 :
	  // cfl = 3.0895
          xmin(0) = 0.229183265; xmin(1) = 0.5;
          break;
	case 7 :
	  // Cfl =  7.0875
	  xmin(0) = 0.110451398065702; xmin(1) = 0.173816271367107;
	  xmin(2) = 0.459433163929695; xmin(3) = 0.652002232653235;
	  break;
	case 8 :
	  // Cfl = 7.8525
	  xmin(0) = 0.135294127286225; xmin(1) = 0.24015308384744;
	  xmin(2) = 0.453046953126355; xmin(3) = 0.695039606659698;
	  break;
	case 10 :
	  // Cfl = 4.7527 (permut = 21)
	  xmin(0) = 0.0; xmin(1) = -0.1; xmin(2) = 0.0; xmin(3) = 0.0021632268153138;
	  break;
        }
      
      bool coef_ok = fct.GetCoefficients(xmin, Abar, B, Bbar, C);
      if (!coef_ok)
        {
          cout << "Bad choice of coefficients" << endl;
          abort();
        }
      
      fct.WriteCoefficients(Abar, B, Bbar, C);
      DISP(Abar); DISP(B); DISP(Bbar); DISP(C);

      TinyMatrix<UnivariatePolynomial<Real_wp>, General, 2, 2> D;
      GetMatrixDz(Abar, B, Bbar, C, Real_wp(0.34), D);
      DISP(D);

      /* Real_wp r_5(0);
      for (int i = 0; i < Abar.GetM(); i++)
	for (int j = 0; j < i; j++)
	  r_5 += B(i)*pow(C(i), 3)*Abar(i, j)*pow(C(j), 5);

	  DISP(r_5); */
      GetAmplificationRKN(Abar, B, Bbar, C, Real_wp(0.34));
      
      //Real_wp cfl = GetCFL_RungeKuttaNystrom(Abar, B, Bbar, C);
      //DISP(sqrt(abs(cfl))); 
      
      RknOrderCondition sys;
      if (order > 8)
	sys.SetOrder(order, true, true);
      else
	sys.SetOrder(order, true, false);
      
      int n = sys.GetM();
      VectReal_wp eval(n);
      sys.EvaluateFunction(Abar, B, Bbar, C, eval);
      sys.CheckCondition(Abar, B, Bbar, C);
      
      DISP(eval); DISP(Norm2(eval));
      if (order <= 8)
	for (int i = 0; i < eval.GetM(); i++)
	  if (abs(eval(i)) > 1e-14)
	    {
	      DISP(i); DISP(eval(i));
	    }
    }
  
  return FinalizeMontjoie();
}
