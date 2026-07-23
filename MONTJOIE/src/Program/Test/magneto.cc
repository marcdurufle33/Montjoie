#include "Montjoie.hxx"

using namespace Montjoie;

void CreatePolygon(VectR3& polygon)
{
  /*polygon.Reallocate(5);
  polygon(0) = R3(0, 0, 0);
  polygon(1) = R3(0, 15, 0.0);
  polygon(2) = R3(20, 15, 0.0);
  polygon(3) = R3(20, 0, 0.0);
  polygon(4) = R3(0, 0);*/

  polygon.Reallocate(13);
  polygon(0) = R3(1, 1, 0);
  polygon(1) = R3(1, 3, 0);
  polygon(2) = R3(1, 3, 0.5);
  polygon(3) = R3(1, 6.505, 0.5);
  polygon(4) = R3(1, 6.505, 1.0);
  polygon(5) = R3(1, 9, 1);
  polygon(6) = R3(9, 9, 1);
  polygon(7) = R3(9, 6.505, 1);
  polygon(8) = R3(9, 6.505, 0.5);
  polygon(9) = R3(9, 3, 0.5);
  polygon(10) = R3(9, 3, 0);
  polygon(11) = R3(9, 1, 0);
  polygon(12) = R3(1, 1, 0);
  
  /*polygon.Reallocate(9);
  polygon(0) = R2(0, 0, 0);
  polygon(1) = R2(0, 15, 1);
  polygon(2) = R2(10, 15, 0.5);
  polygon(3) = R2(10, 0, -0.2);
  polygon(4) = R2(20, 0);
  polygon(5) = R2(20, 15);
  polygon(6) = R2(10, 15);
  polygon(7) = R2(10, 0);
  polygon(8) = R2(0, 0);*/
}

void ComputeFieldE()
{
  //double xmin = -5.0, xmax = 5.0, ymin = -5.0, ymax = 5.0, zmin = -5.0, zmax = 5.0;
  double xmin = 1.0, xmax = 3.0, ymin = 0.0, ymax = 1.0, zmin = 0.0, zmax = 1.0;
  int nbx = 101, nby = 101, nbz = 101;

  GridInterpolationFull<Dimension3> grid;
  R3 pt0(2.0, 0.0, 0.0), ptMin(xmin, ymin, zmin), ptMax(xmax, ymax, zmax);
  grid.SetThreePlanesOutput(pt0, ptMin, ptMax, nbx, nby, nbz);

  VectR3 Points3D;
  grid.GenerateGridPoints(Points3D);

  int r = 12;
  Globatto<Real_wp> lob;
  lob.ConstructQuadrature(r);
  
  Vector<TinyVector<Complex_wp, 3> > fieldE(Points3D.GetM());
  Real_wp alpha = log(1e6), beta = sqrt(alpha/pi_wp);
  R3 ptRp, vec_v, ptR; Real_wp omega = 2.0*pi_wp; TinyVector<Complex_wp, 3> Jc, vloc;
  Matrix3_3sym_Complex_wp hessian_phi, dyadic_G; R3_Complex_wp grad_phi; Complex_wp phi;
  Complex_wp czero(0, 0); Real_wp omega2 = omega*omega; R3 center(0, 0, 1);
  for (int j1 = 0; j1 < lob.GetNbPointsQuad(); j1++)
    for (int j2 = 0; j2 < lob.GetNbPointsQuad(); j2++)
      for (int j3 = 0; j3 < lob.GetNbPointsQuad(); j3++)
	{
	  ptRp.Init(2.0*lob.Points(j1)-1.0, 2.0*lob.Points(j2)-1.0, 2.0*lob.Points(j3)-1.0);
	  Real_wp coef = 8.0*lob.Weights(j1)*lob.Weights(j2)*lob.Weights(j3)
	    *beta*exp(-alpha*DotProd(ptRp, ptRp));

          ptRp += center;
	  Jc.Init(Complex_wp(coef, 0), czero, czero);
	  
	  for (int k = 0; k < Points3D.GetM(); k++)
	    {
              TransparencySolver_Base::ComputeGreenKernel(Points3D(k), ptRp, omega, phi,
                                                          grad_phi, hessian_phi);              
              
	      // G = \phi(x,y) I + 1/k^2 \nabla_y \nabla_y \phi(x,y)
	      Copy(hessian_phi, dyadic_G); Mlt(Complex_wp(1.0/omega2), dyadic_G);
        
	      dyadic_G(0, 0) += phi;
	      dyadic_G(1, 1) += phi;
	      dyadic_G(2, 2) += phi;
              
	      Mlt(dyadic_G, Jc, vloc);
              fieldE(k) += vloc;
              
              ptR.Init(ptRp(0), ptRp(1), -ptRp(2));
              TransparencySolver_Base::ComputeGreenKernel(Points3D(k), ptR, omega, phi,
                                                          grad_phi, hessian_phi);

              Copy(hessian_phi, dyadic_G); Mlt(Complex_wp(1.0/omega2), dyadic_G);
              
	      dyadic_G(0, 0) += phi;
	      dyadic_G(1, 1) += phi;
	      dyadic_G(2, 2) += phi;
              
	      Mlt(dyadic_G, Jc, vloc);

	      fieldE(k) -= vloc;
              
	      //TimesProd(grad_phi, Jc, vloc);
	      //fieldE(k) += vloc;
	    }
	}
  
  VectComplex_wp val(Points3D.GetM());
  for (int k = 0; k < Points3D.GetM(); k++)
    val(k) = fieldE(k)(0);
  
  int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
  WriteMatlab(val, grid, "totoX.dat", double_prec);

  for (int k = 0; k < Points3D.GetM(); k++)
    val(k) = fieldE(k)(1);
  
  WriteMatlab(val, grid, "totoY.dat", double_prec);

  for (int k = 0; k < Points3D.GetM(); k++)
    val(k) = fieldE(k)(2);
  
  WriteMatlab(val, grid, "totoZ.dat", double_prec);

}


void ComputeFieldH(const VectR3& polygon)
{
  double xmin = -2.0, xmax = 12.0, ymin = -2.0, ymax = 12.0;
  int nbx = 201, nby = 201;
  
  GridInterpolationFull<Dimension2> grid;
  grid.SetPlaneOutput(xmin, xmax, ymin, ymax, nbx, nby);
  
  VectR2 Points2D; VectReal_wp TetaPoints;
  grid.GenerateGridPoints(Points2D, TetaPoints);
  
  int r = 40;
  Globatto<Real_wp> lob;
  lob.ConstructQuadrature(r);

  VectComplex_wp val(Points2D.GetM());
  Vector<TinyVector<Complex_wp, 3> > fieldH(Points2D.GetM());
  val.Zero();
  
  R3 vec_u, vec_v, vec_w; R3_Complex_wp vec_wc, grad_phi;
  R3 ptR, ptRp;
  Real_wp coef = 1.0/(4.0*pi_wp);
  Real_wp pos_ecoute = 1.2;
  Real_wp pos_boucle = 0.0;
  Real_wp pos_metal = -0.1; // position du metal
  Real_wp h = pos_ecoute - pos_boucle;
  Real_wp pos_sym = 2*pos_metal - pos_boucle;
  Real_wp h2 = pos_ecoute - pos_sym;
  bool dirichlet = false;
  bool magnetostatic = true;
  Real_wp freq = 1.6e3;
  Real_wp omega = 2.0*pi_wp*freq / PhysicalConstant::speed_light;
  R3 ur, vr, wr;
  bool use_numerical_integration = false;
  bool compute_H = true;

  bool presence_gradins = true;
  int nb_marches = 2;
  Real_wp h_gradin = 1.0;
  VectReal_wp y_div(3);
  y_div(0) = 3.0; y_div(1) = 6.505; y_div(2) = 10.0;
  DISP(y_div);

  VectR3 Points3D(Points2D.GetM());
  for (int k = 0; k < Points2D.GetM(); k++)
    {
      int num = 0;
      double y = Points2D(k)(1);
      while ((num < y_div.GetM()) && (y >= y_div(num)))
        num++;

      if (num == y_div.GetM())
        num = 0;

      Real_wp z = h + Real_wp(num)*h_gradin/nb_marches;      
      Points3D(k).Init(Points2D(k)(0), Points2D(k)(1), z);
    }
  
  for (int i = 0; i < polygon.GetM()-1; i++)
    {
      vec_u = polygon(i+1) - polygon(i);

      Real_wp L = Norm2(vec_u);
      Real_wp invNorm = 1.0/L;
      ur = invNorm * vec_u;
      GetVectorPlane(ur, vr, wr);
      //DISP(i); DISP(ur); DISP(vr); DISP(L); DISP(wr);

      if (use_numerical_integration)
	for (int j = 0; j < lob.GetNbPointsQuad(); j++)
	  {
	    Real_wp zeta = lob.Points(j);
	    ptRp = (1.0-zeta)*polygon(i) + zeta*polygon(i+1);
	    for (int k = 0; k < Points3D.GetM(); k++)
	      {
		ptR = Points3D(k);
                vec_v = ptR - ptRp;
		
		if (magnetostatic)
		  {
                    if (compute_H)
                      {
                        TimesProd(vec_u, vec_v, vec_w);
                        Real_wp vloc = lob.Weights(j)/pow(Norm2(vec_v), 3.0);                      
                        Add(vloc, vec_w, fieldH(k));

                        if (dirichlet)
                          {
                            vec_v.Init(vec_v(0), vec_v(1), h2);
                            TimesProd(vec_u, vec_v, vec_w);
                            Real_wp vloc = -lob.Weights(j)/pow(Norm2(vec_v), 3.0);
                            Add(vloc, vec_w, fieldH(k));
                          }
                      }
                    else
                      {
                        Real_wp vloc = lob.Weights(j) / Norm2(vec_v);
                        Add(vloc, vec_u, fieldH(k));
                      }
		  }
		else
		  {
		    Real_wp T = Norm2(vec_v), inv_T = Real_wp(1)/T;
		    Complex_wp arg = Iwp*omega*T;
		    
		    // evaluation of green function
		    Complex_wp phi = exp(arg) * inv_T;
		    if (compute_H)
                      {
                        Complex_wp phi_div_T = phi * inv_T;
                        
                        Complex_wp alpha0 = (-Iwp*omega + inv_T)*phi_div_T;
                        
                        grad_phi = vec_v;    
                        Mlt( alpha0 , grad_phi);    
                        TimesProd(vec_u, grad_phi, vec_wc);
                      }
                    else
                      vec_wc = phi*vec_u;
                    
		    Real_wp vloc = lob.Weights(j);
		    Add(vloc, vec_wc, fieldH(k));
		  }
	      }
	  }
      else
	{
	  for (int k = 0; k < Points3D.GetM(); k++)
	    {
	      R3 diff = Points3D(k)-polygon(i);
	      Real_wp beta = DotProd(diff, ur);
	      Real_wp gamma = DotProd(diff, vr);
              Real_wp h0 = DotProd(diff, wr);
	      //DISP(beta); DISP(gamma); DISP(h);
	      Real_wp offset = h0*h0 + gamma*gamma;
              Real_wp offset2 = h2*h2 + gamma*gamma;
	      Real_wp valInt;
              if (compute_H)
                {
                  valInt = (L - beta) / sqrt(offset + square(L-beta)) + beta / sqrt(offset + beta*beta);
                  valInt /= offset;
                  fieldH(k) += -h0*valInt*vr + gamma*valInt*wr;
                  
                  if (dirichlet)
                    {
                      valInt = (L - beta) / sqrt(offset2 + square(L-beta)) + beta / sqrt(offset2 + beta*beta);
                      valInt /= offset2;
                      //fieldH(k)(0) += h2*valInt*vr(0);
                      //fieldH(k)(1) += h2*valInt*vr(1);
                      //fieldH(k)(2) -= gamma*valInt;
                      abort();
                    }
                }
              else
                {
                  valInt = -log(sqrt(square(L-beta) + offset) + beta-L) + log(sqrt(offset+beta*beta) + beta);
                  fieldH(k) += valInt*ur;
                  //fieldH(k)(0) += valInt*ur(0);
                  //fieldH(k)(1) += valInt*ur(1);
                  //fieldH(k)(2) = 0.0;
                }
              
	    }
	}
    }
  
  int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
  for (int k = 0; k < Points2D.GetM(); k++)
    fieldH(k) *= coef;

  //DISP(fieldH);
  for (int k = 0; k < Points2D.GetM(); k++)
    val(k) = fieldH(k)(0);
  
  WriteMatlab(val, grid, "totoX.dat", double_prec);

  for (int k = 0; k < Points2D.GetM(); k++)
    val(k) = fieldH(k)(1);
  
  WriteMatlab(val, grid, "totoY.dat", double_prec);

  for (int k = 0; k < Points2D.GetM(); k++)
    val(k) = fieldH(k)(2);

  WriteMatlab(val, grid, "totoZ.dat", double_prec);
  
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  //ComputeFieldE();
  
  VectR3 polygon;
  CreatePolygon(polygon);
  ComputeFieldH(polygon);

  return FinalizeMontjoie();
}
