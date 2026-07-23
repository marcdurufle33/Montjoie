#ifndef MONTJOIE_FILE_KDV_PROBLEM_CXX

namespace Montjoie
{
    
  KdvEquation::KdvEquation() : GenericWaterWaveEquation()
  {
  }
    
    
  template<class GenericPb>
  void KdvEquation::ComputePhysicalProperty(const GenericPb& var, double xmin, double xmax)
  {      
    // DISP(var.mesh.GetOrder());
    // we approximate bottom with higher degree in order to have
    // a correct approximation of Cxxx
    InitOrder(var.mesh.GetOrder() + 3, var);
    int Nlob = lob.GetNbPointsQuad();
    VectReal_wp Xhigh(Nlob*var.mesh.GetNbElt());
    int nb = 0;
    
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      {
	Real_wp h = var.mesh.Vertex(i+1) - var.mesh.Vertex(i);
	for (int j = 0; j < Nlob; j++)
	  Xhigh(nb++) = var.mesh.Vertex(i) + h*lob.Points(j);
      }
      
    VectReal_wp Bhigh, Chigh, C5;
    ComputeBottom(xmin, xmax, Xhigh, Bhigh, Chigh);
    C5 = Chigh;
    for (int i = 0; i < C5.GetM(); i++)
      C5(i) = pow(Chigh(i), 5.0);
    
    // computation of Cx, Cxx and Cxxx
    VectReal_wp CxDof, C5xDof, C5xxDof, C5xxxDof;
    GetGradient(var, Chigh, CxDof);
    
    GetGradient(var, C5, C5xDof);
    GetGradient(var, C5xDof, C5xxDof);
    GetGradient(var, C5xxDof, C5xxxDof);
    
    GetUquadrature(var, Chigh, Cquad);
    GetUquadrature(var, CxDof, CxQuad);
    
    GetUquadrature(var, C5, C5Quad);
    GetUquadrature(var, C5xDof, C5xQuad);
    GetUquadrature(var, C5xxDof, C5xxQuad);
    GetUquadrature(var, C5xxxDof, C5xxxQuad);
    
    invCquad = Cquad;
    for (int i = 0; i < Cquad.GetM(); i++)
      invCquad(i) = 1.0/Cquad(i);
    
  }

  
  KdvStiffnessOperatorDG::KdvStiffnessOperatorDG(KdvEquationDG& var_kdv) : var(var_kdv)
  {
  }
    
    
  template<class Vector1, class Vector2>
  void KdvStiffnessOperatorDG
  ::EvaluateNumericalFlux(int i, int k, const Vector1& Um, const Vector1& Up, Vector2& flux) const
  {
    flux.Zero();
    switch (var.type_model)
      {
      case KdvEquation::GENTLE :
        {
          Real_wp flux_m = (var.Cquad(i) + 0.75*var.epsilon*Um(0))*Um(0);
          Real_wp flux_p = (var.Cquad(i) + 0.75*var.epsilon*Up(0))*Up(0);
          flux(0) = 0.5*(flux_p - flux_m);    
          if (k == 0)
            {
              flux(0) -= 0.5*var.alpha_lax*(Um(0) - Up(0));
              flux(1) = -Up(0) + Um(0);
            }
          else
            {
              flux(0) += 0.5*var.alpha_lax*(Um(0) - Up(0)) + var.mu/6*(Up(2)-Um(2));
              flux(2) = -Up(1) + Um(1);
            }
        }
        break;
      case KdvEquation::STRONG :
        {
          Real_wp flux_m = (var.Cquad(i) + 0.75*var.invCquad(i)*var.epsilon*Um(0))*Um(0);
          Real_wp flux_p = (var.Cquad(i) + 0.75*var.invCquad(i)*var.epsilon*Up(0))*Up(0);
          flux(0) = 0.5*(flux_p - flux_m);    
          if (k == 0)
            {
              flux(0) -= 0.5*var.alpha_lax*(Um(0) - Up(0));
              flux(1) = -Up(0) + Um(0);
            }
          else
            {
              flux(0) += 0.5*var.alpha_lax*(Um(0) - Up(0)) + var.mu/6*var.C5Quad(i)*(Up(2)-Um(2));
              flux(2) = -Up(1) + Um(1);
            }
        }
        break;
      case KdvEquation::ORIGINAL :
        {
          Real_wp flux_m = (var.Cquad(i) + 0.75*var.invCquad(i)*var.epsilon*Um(0))*Um(0);
          Real_wp flux_p = (var.Cquad(i) + 0.75*var.invCquad(i)*var.epsilon*Up(0))*Up(0);
          flux(0) = 0.5*(flux_p - flux_m);    
          if (k == 0)
            {
              flux(0) -= 0.5*var.alpha_lax*(Um(0) - Up(0));
              flux(1) = -Up(0) + Um(0);
            }
          else
            {
              flux(0) += 0.5*var.alpha_lax*(Um(0) - Up(0)) + var.mu/6*var.C5Quad(i)*(Up(2)-Um(2));
              flux(2) = -Up(1) + Um(1);
            }
        }
        break;
      }
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void KdvStiffnessOperatorDG
  ::EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um, const Vector1& Up,
                                    Vector2& flux,
                                    Matrix1& df_dum, Matrix1& df_dup) const
  {
    flux.Zero(); df_dum.Zero(); df_dup.Zero();
    switch(var.type_model)
      {
      case KdvEquation::GENTLE :
        {
          Real_wp flux_m = (var.Cquad(i) + 0.75*var.epsilon*Um(0))*Um(0);
          Real_wp flux_p = (var.Cquad(i) + 0.75*var.epsilon*Up(0))*Up(0);
          flux(0) = 0.5*(flux_p - flux_m);    
          if (k == 0)
            {
              flux(0) -= 0.5*var.alpha_lax*(Um(0) - Up(0));
              flux(1) = -Up(0)+Um(0);
            }
          else
            {
              flux(0) += 0.5*var.alpha_lax*(Um(0) - Up(0)) + var.mu/6*(Up(2)-Um(2));
              flux(2) = -Up(1)+Um(1);
            }
          
          df_dum(0,0) = -0.5*var.Cquad(i) - 0.75*var.epsilon*Um(0);
          df_dup(0,0) = 0.5*var.Cquad(i) + 0.75*var.epsilon*Up(0);
          
          if (k == 0)
            {
              df_dum(0,0) -= 0.5*var.alpha_lax;
              df_dup(0,0) += 0.5*var.alpha_lax;	
              
              df_dum(1,0) = 1.0;
              df_dup(1,0) = -1.0;
            }
          else
            {
              df_dum(0,0) += 0.5*var.alpha_lax;
              df_dup(0,0) -= 0.5*var.alpha_lax;
              
              df_dum(0,2) = -var.mu/6;
              df_dup(0,2) = var.mu/6;
              
              df_dum(2,1) = 1.0;
              df_dup(2,1) = -1.0;
            }
        }
        break;
      case KdvEquation::STRONG :
        {
          Real_wp flux_m = (var.Cquad(i) + 0.75*var.invCquad(i)*var.epsilon*Um(0))*Um(0);
          Real_wp flux_p = (var.Cquad(i) + 0.75*var.invCquad(i)*var.epsilon*Up(0))*Up(0);
          flux(0) = 0.5*(flux_p - flux_m);    
          if (k == 0)
            {
              flux(0) -= 0.5*var.alpha_lax*(Um(0) - Up(0));
              flux(1) = -Up(0)+Um(0);
            }
          else
            {
              flux(0) += 0.5*var.alpha_lax*(Um(0) - Up(0)) + var.mu/6*var.C5Quad(i)*(Up(2)-Um(2));
              flux(2) = -Up(1)+Um(1);
            }
          
          df_dum(0,0) = -0.5*var.Cquad(i) - 0.75*var.invCquad(i)*var.epsilon*Um(0);
          df_dup(0,0) = 0.5*var.Cquad(i) + 0.75*var.invCquad(i)*var.epsilon*Up(0);
          
          if (k == 0)
            {
              df_dum(0,0) -= 0.5*var.alpha_lax;
              df_dup(0,0) += 0.5*var.alpha_lax;
              
              df_dum(1,0) = 1.0;
              df_dup(1,0) = -1.0;
            }
          else
            {
              df_dum(0,0) += 0.5*var.alpha_lax;
              df_dup(0,0) -= 0.5*var.alpha_lax;
              
              df_dum(0,2) = -var.mu/6*var.C5Quad(i);
              df_dup(0,2) = var.mu/6*var.C5Quad(i);
              
              df_dum(2,1) = 1.0;
              df_dup(2,1) = -1.0;
            }
        }
        break;
      case KdvEquation::ORIGINAL :
        {
          Real_wp flux_m = (var.Cquad(i) + 0.75*var.invCquad(i)*var.epsilon*Um(0))*Um(0);
          Real_wp flux_p = (var.Cquad(i) + 0.75*var.invCquad(i)*var.epsilon*Up(0))*Up(0);
          flux(0) = 0.5*(flux_p - flux_m);    
          if (k == 0)
            {
              flux(0) -= 0.5*var.alpha_lax*(Um(0) - Up(0));
              flux(1) = -Up(0)+Um(0);
            }
          else
            {
              flux(0) += 0.5*var.alpha_lax*(Um(0) - Up(0)) + var.mu/6*var.C5Quad(i)*(Up(2)-Um(2));
              flux(2) = -Up(1)+Um(1);
            }
          
          df_dum(0,0) = -0.5*var.Cquad(i) - 0.75*var.invCquad(i)*var.epsilon*Um(0);
          df_dup(0,0) = 0.5*var.Cquad(i) + 0.75*var.invCquad(i)*var.epsilon*Up(0);
          
          if (k == 0)
            {
              df_dum(0,0) -= 0.5*var.alpha_lax;
              df_dup(0,0) += 0.5*var.alpha_lax;
              
              df_dum(1,0) = 1.0;
              df_dup(1,0) = -1.0;
            }
          else
            {
              df_dum(0,0) += 0.5*var.alpha_lax;
              df_dup(0,0) -= 0.5*var.alpha_lax;
              
              df_dum(0,2) = -var.mu/6*var.C5Quad(i);
              df_dup(0,2) = var.mu/6*var.C5Quad(i);
              
              df_dum(2,1) = 1.0;
              df_dup(2,1) = -1.0;
            }
        }
      }
  }
  
  
  template<class Vector1, class Vector2>
  void KdvStiffnessOperatorDG::EvaluateFunction(int i, const Vector1& u, const Vector1& du,
						Vector2& g, Vector2& h) const
  {
    switch (var.type_model)
      {
      case KdvEquation::GENTLE :
        {
          h(0) = du(0)*(var.Cquad(i) + 1.5*var.epsilon*u(0)) 
            + 0.5*var.CxQuad(i)*u(0) + var.mu/6*du(2); 
          h(1) = u(1) - du(0);
          h(2) = u(2) - du(1);
        }
        break;
      case KdvEquation::STRONG :
        {
          Real_wp h_cte = 0.5*var.CxQuad(i) + var.mu/48*var.C5xxxQuad(i);
          h(0) = h_cte*u(0) - 0.5*var.epsilon*square(var.invCquad(i)*u(0))*var.CxQuad(i);
          h(0) += 0.125*var.mu*var.C5xxQuad(i)*u(1) + 0.25*var.mu*var.C5xQuad(i)*u(2);
          h(0) += du(0)*(var.Cquad(i) + 1.5*var.invCquad(i)*var.epsilon*u(0))
            + var.mu/6*var.C5Quad(i)*du(2);
          h(1) = u(1) - du(0);
          h(2) = u(2) - du(1);
        }
        break;
      case KdvEquation::ORIGINAL :
        {
          Real_wp h_cte = 0.5*var.CxQuad(i);
          h(0) = h_cte*u(0) + du(0)*(var.Cquad(i) + 1.5*var.invCquad(i)*var.epsilon*u(0))
            + du(2)*var.mu/6*var.C5Quad(i);
          h(1) = u(1) - du(0);
          h(2) = u(2) - du(1);
        }
        break;
      }
  }
  
  
  template<class Vector1, class Vector2, class Matrix1>
  void KdvStiffnessOperatorDG
  ::EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                       Vector2& g, Vector2& h, Matrix1& dg_du,
                       Matrix1& dh_du, Matrix1& dg_ddu, Matrix1& dh_ddu) const
  {
    switch (var.type_model)
      {
      case KdvEquation::GENTLE :
        {
          h(0) = 0.5*var.CxQuad(i)*u(0) 
            + du(0)*(var.Cquad(i) + 1.5*var.epsilon*u(0)) + var.mu/6*du(2); 
          h(1) = u(1) - du(0);
          h(2) = u(2) - du(1);
          
          dh_du(0, 0) = 0.5*var.CxQuad(i) + 1.5*var.epsilon*du(0);
          dh_du(1, 1) = 1.0;
          dh_du(2, 2) = 1.0;
          
          dh_ddu(0, 0) = var.Cquad(i) + 1.5*var.epsilon*u(0);
          dh_ddu(0, 2) = var.mu/6;
          
          dh_ddu(1, 0) = -1.0;
          dh_ddu(2, 1) = -1.0;
        }
        break;
      case KdvEquation::STRONG :
        {
          Real_wp h_cte = 0.5*var.CxQuad(i) + var.mu/48*var.C5xxxQuad(i);
          h(0) = h_cte*u(0) - 0.5*var.epsilon*square(var.invCquad(i)*u(0))*var.CxQuad(i);
          h(0) += 0.125*var.mu*var.C5xxQuad(i)*u(1) + 0.25*var.mu*var.C5xQuad(i)*u(2);
          h(0) += du(0)*(var.Cquad(i) + 1.5*var.invCquad(i)*var.epsilon*u(0))
            + du(2)*var.mu/6*var.C5Quad(i);
          h(1) = u(1) - du(0);
          h(2) = u(2) - du(1);
          
          dh_du(0, 0) = h_cte - var.epsilon*square(var.invCquad(i))*u(0)*var.CxQuad(i);
          dh_du(0, 0) += 1.5*du(0)*var.invCquad(i)*var.epsilon;
          dh_du(0, 1) = 0.125*var.mu*var.C5xxQuad(i);
          dh_du(0, 2) = 0.25*var.mu*var.C5xQuad(i);
          
          dh_du(1, 1) = 1.0;
          dh_du(2, 2) = 1.0;
          
          dh_ddu(0, 0) = var.Cquad(i) + 1.5*var.invCquad(i)*var.epsilon*u(0);
          dh_ddu(0, 2) = var.mu/6*var.C5Quad(i);
          
          dh_ddu(1, 0) = -1.0;
          dh_ddu(2, 1) = -1.0;
        }
        break;
      case KdvEquation::ORIGINAL :
        {
          Real_wp h_cte = 0.5*var.CxQuad(i);
          h(0) = h_cte*u(0) + du(0)*(var.Cquad(i) + 1.5*var.invCquad(i)*var.epsilon*u(0))
            + du(2)*var.mu/6*var.C5Quad(i);
          h(1) = u(1) - du(0);
          h(2) = u(2) - du(1);	
          
          dh_du(0,0) = h_cte + du(0)*1.5*var.invCquad(i)*var.epsilon;
          dh_du(1, 1) = 1.0;
          dh_du(2, 2) = 1.0;
          
          dh_ddu(0, 0) = var.Cquad(i) + 1.5*var.invCquad(i)*var.epsilon*u(0);
          dh_ddu(0, 2) = var.mu/6*var.C5Quad(i);
          dh_ddu(1, 0) = -1.0;
          dh_ddu(2, 1) = -1.0;
        }
        break;
      }
  }

  
  KdvEquationDG::KdvEquationDG() : var_stiff(*this)
  {
    alpha_lax = 3.0;
  }
  
  
  KdvEquationDG::KdvEquationDG(const KdvEquationDG& var) : KdvEquation(var), var_stiff(*this)
  {
    alpha_lax = var.alpha_lax;
  }
  
  
  KdvStiffnessOperatorFem::KdvStiffnessOperatorFem(const KdvEquationFem& var_eq) : var(var_eq)
  {
  }
  
  
  template<class Vector1, class Vector2>
  void KdvStiffnessOperatorFem::
  EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                   Vector2& g, Vector2& h) const
  {
    h(1) = u(1) - du(0);
    h(2) = u(2) - du(1);

    switch (var.type_model)
      {
      case KdvEquation::GENTLE :
        {
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0)
            + 1.5*var.epsilon*u(0)*du(0) + var.mu/6*du(2);
        }
        break;
      case KdvEquation::STRONG :
        {
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0)
            + 1.5*var.invCquad(i)*var.epsilon*u(0)*du(0) + 0.125*var.mu*var.C5xxQuad(i)*du(0)
            + 0.25*var.mu*var.C5xQuad(i)*du(1) + var.mu/6*var.C5Quad(i)*du(2)
            + var.mu/48*var.C5xxxQuad(i)*u(0) 
            - 0.5*var.epsilon*square(var.invCquad(i)*u(0))*var.CxQuad(i);
          
        }
        break;
      case KdvEquation::ORIGINAL :
        {
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0)
            + 1.5*var.invCquad(i)*var.epsilon*u(0)*du(0) + var.mu/6*var.C5Quad(i)*du(2);          
        }
        break;
      }
  }
   
  
  template<class Vector1, class Vector2, class Matrix1>
  void KdvStiffnessOperatorFem::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                     Vector2& g, Vector2& h, Matrix1& dg_du,
                     Matrix1& dh_du, Matrix1& dg_ddu, Matrix1& dh_ddu) const
  {
    h(1) = u(1) - du(0);
    h(2) = u(2) - du(1);
    
    dh_du(1, 1) = 1.0;
    dh_ddu(1, 0) = -1.0;
    
    dh_du(2, 2) = 1.0;
    dh_ddu(2, 1) = -1.0;
    
    switch (var.type_model)
      {
      case KdvEquation::GENTLE :
        {
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0)
            + 1.5*var.epsilon*u(0)*du(0) + var.mu/6*du(2);

          dh_du(0, 0) = 0.5*var.CxQuad(i) + 1.5*var.epsilon*du(0);
          dh_ddu(0, 0) = var.Cquad(i) + 1.5*var.epsilon*u(0);
          dh_ddu(0, 2) = var.mu/6;
          
        }
        break;
      case KdvEquation::STRONG :
        {
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0)
            + 1.5*var.invCquad(i)*var.epsilon*u(0)*du(0) + 0.125*var.mu*var.C5xxQuad(i)*du(0)
            + 0.25*var.mu*var.C5xQuad(i)*du(1) + var.mu/6*var.C5Quad(i)*du(2)
            + var.mu/48*var.C5xxxQuad(i)*u(0) 
            - 0.5*var.epsilon*square(var.invCquad(i)*u(0))*var.CxQuad(i);
          
          dh_du(0, 0) = 0.5*var.CxQuad(i) + var.mu/48*var.C5xxxQuad(i) 
            + 1.5*var.invCquad(i)*var.epsilon*du(0)
            - var.epsilon*square(var.invCquad(i))*var.CxQuad(i)*u(0);
          
          dh_ddu(0, 0) = var.Cquad(i) + 1.5*var.invCquad(i)*var.epsilon*u(0)
            + 0.125*var.mu*var.C5xxQuad(i);
          
          dh_ddu(0, 1) = 0.25*var.mu;
          dh_ddu(0, 2) = var.mu/6*var.C5Quad(i);
        }
        break;
      case KdvEquation::ORIGINAL :
        {
          h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0)
            + 1.5*var.invCquad(i)*var.epsilon*u(0)*du(0) + var.mu/6*var.C5Quad(i)*du(2);          

          dh_du(0, 0) = 0.5*var.CxQuad(i) + 1.5*var.invCquad(i)*var.epsilon*du(0);
          
          dh_ddu(0, 0) = var.Cquad(i) + 1.5*var.invCquad(i)*var.epsilon*u(0);
          dh_ddu(0, 1) = 0.25*var.mu;
          dh_ddu(0, 2) = var.mu/6*var.C5Quad(i);
        }
        break;
      }
  }
 
  
  KdvEquationFem::KdvEquationFem() : var_stiff(*this)
  {
  }
  
  KdvEquationFem::KdvEquationFem(const KdvEquationFem& A) : KdvEquation(A), var_stiff(*this)
  {
  }

  inline bool KdvEquationFem::LinearMassTerm()
  {
    return false;
  }


  /*************************
   * Spectral KdV Equation *
   *************************/

  
  //! default constructor
  KdvEquationSpectral::KdvEquationSpectral()
  {
    constant_coef = false;
  }
  
  
  //! computation of c and cx 
  void KdvEquationSpectral
  ::ComputePhysicalProperty(FftRealInterface& fft, const VectReal_wp& nu,
			    double xmin, double xmax, int N)
  {      
    VectReal_wp Xdof(N);
    Real_wp dx = (xmax-xmin)/N;
    for (int i = 0; i < N; i++)
      Xdof(i) = xmin + i*dx;
    
    VectReal_wp Bdof;
    ComputeBottom(xmin, xmax, Xdof, Bdof, Cquad);
    
    Cquad.WriteText("bottom_sp.dat");
    C5Quad.Reallocate(N);
    invCquad.Reallocate(N);
    for (int i = 0; i < N; i++)
      {
        C5Quad(i) = pow(Cquad(i), 5.0);
        invCquad(i) = 1.0/Cquad(i);
      }
    
    // computation of Cx, (c^5)_x, (c^5)_xx, (c^5)_xxx
    VectComplex_wp Chat(N/2+1), C5hat(N/2+1);
    VectComplex_wp C5hat_x(N/2+1), C5hat_xx(N/2+1);
    fft.ApplyForward(Cquad, Chat);
    fft.ApplyForward(C5Quad, C5hat);
    
    for (int i = 0; i < nu.GetM(); i++)
      {
        Chat(i) *= Iwp*nu(i);
        C5hat(i) *= Iwp*nu(i);
        C5hat_x(i) = Iwp*nu(i)*C5hat(i);
        C5hat_xx(i) = Iwp*nu(i)*C5hat_x(i);
      }
    
    CxQuad.Reallocate(N);
    C5xQuad.Reallocate(N);
    C5xxQuad.Reallocate(N);
    C5xxxQuad.Reallocate(N);
    fft.ApplyInverse(Chat, CxQuad);
    fft.ApplyInverse(C5hat, C5xQuad);
    fft.ApplyInverse(C5hat_x, C5xxQuad);
    fft.ApplyInverse(C5hat_xx, C5xxxQuad);
    
    ProdReal.Reallocate(N);
    XnReel.Reallocate(N);
    dXreal.Reallocate(N);
    d2Xreal.Reallocate(N);
    d3Xreal.Reallocate(N);
    
    XnHat.Reallocate(N/2+1);
    dXn.Reallocate(N/2+1);
    d2Xn.Reallocate(N/2+1);
    d3Xn.Reallocate(N/2+1);
    
    coef_Gamma1.Reallocate(nu.GetM());
    coef_Gamma2.Reallocate(nu.GetM());
    coef_Gamma3.Reallocate(nu.GetM());
    for (int i = 0; i < nu.GetM(); i++)
      {
        coef_Gamma1(i) = Iwp*nu(i);
        coef_Gamma2(i) = -nu(i)*nu(i);
        coef_Gamma3(i) = -Iwp*pow(nu(i), 3.0);
      }
  }
  
  
  void KdvEquationSpectral
  ::InitSinus(const Real_wp& dt, const VectReal_wp& nu)
  {
    coef_Gamma1.Reallocate(nu.GetM());
    coef_Gamma2.Reallocate(nu.GetM());
    coef_Gamma3.Reallocate(nu.GetM());
    for (int i = 0; i < nu.GetM(); i++)
      {
        coef_Gamma1(i) = Iwp/dt*sin(dt*nu(i));
        coef_Gamma2(i) = -1.0/dt*sin(dt*nu(i)*nu(i));
        coef_Gamma3(i) = -Iwp/dt*sin(dt*pow(nu(i), 3.0));
      }
  }
  
  
  //! applying mass operator
  void KdvEquationSpectral
  ::ApplyMass(const VectReal_wp& X, FftRealInterface& fft, VectReal_wp& nu, VectReal_wp& Y)
  {
    Seldon::Copy(X, Y);
  }
  
  
  //! adding stiffness operator
  void KdvEquationSpectral
  ::AddFunction(const Real_wp& alpha, const VectReal_wp& X,
                FftRealInterface& fft, VectReal_wp& nu,
                VectReal_wp& Y, bool solve_mass)
  {
    VectReal_wp Prod(Y);
    EvaluateFunction(0.0, X, fft, nu, Prod, solve_mass);    
    Add(alpha, Prod, Y);
  }
  
  
  //! evaluation of the operator that evolves in time
  void KdvEquationSpectral
  ::EvaluateFunction(const Real_wp& t, const VectComplex_wp& X,
		     FftRealInterface& fft, VectReal_wp& nu,
		     VectComplex_wp& Y, bool solve_mass)
  {
    fft.ApplyInverse(X, XnReel);
    
    GetFunction(X, XnReel, fft, nu, ProdReal);
    
    fft.ApplyForward(ProdReal, Y);
  }
  
  
  //! evaluation of the operator that evolves in time
  void KdvEquationSpectral
  ::EvaluateFunction(const Real_wp& t, const VectReal_wp& X,
		     FftRealInterface& fft, VectReal_wp& nu,
		     VectReal_wp& Y, bool solve_mass)
  {
    fft.ApplyForward(X, XnHat);
    
    GetFunction(XnHat, X, fft, nu, Y);    
  }    

  
  void KdvEquationSpectral
  ::GetFunction(const VectComplex_wp& Xn, const VectReal_wp& Xreal,
                FftRealInterface& fft, VectReal_wp& nu, VectReal_wp& Y)
  {
    for (int i = 0; i < Xn.GetM(); i++)
      dXn(i) = Iwp*nu(i)*Xn(i);
    
    fft.ApplyInverse(dXn, dXreal);
    
    if (type_model == GENTLE)
      {	
	for (int i = 0; i < Xreal.GetM(); i++)
	  Y(i) = -(Cquad(i)*dXreal(i) + 0.5*CxQuad(i)*Xreal(i) + 1.5*epsilon*Xreal(i)*dXreal(i));

        for (int i = 0; i < XnHat.GetM(); i++)
          d3Xn(i) = coef_Gamma3(i)*Xn(i);
	
	fft.ApplyInverse(d3Xn, d3Xreal);
        
        for (int i = 0; i < Y.GetM(); i++)
          Y(i) -= mu/6.0*d3Xreal(i);
      }
    else
      {
        for (int i = 0; i < Xn.GetM(); i++)
          {
            d2Xn(i) = coef_Gamma2(i)*Xn(i);
            d3Xn(i) = coef_Gamma3(i)*Xn(i);
          }
        
        fft.ApplyInverse(d2Xn, d2Xreal);
        fft.ApplyInverse(d3Xn, d3Xreal);
	
        if (type_model == STRONG)
          {
            for (int i = 0; i < Xreal.GetM(); i++)
              Y(i) = -(Cquad(i)*dXreal(i) + 0.5*CxQuad(i)*Xreal(i)
                       + 1.5*epsilon*invCquad(i)*Xreal(i)*dXreal(i)
                       - 0.5*CxQuad(i)*square(invCquad(i)*Xreal(i)))
                - mu/6.0*(C5Quad(i)*d3Xreal(i) + 1.5*C5xQuad(i)*d2Xreal(i)
                          + 0.75*C5xxQuad(i)*dXreal(i) + 0.125*C5xxxQuad(i)*Xreal(i));
          }
        else
          {
            for (int i = 0; i < Xreal.GetM(); i++)
              Y(i) = -(Cquad(i)*dXreal(i) + 0.5*CxQuad(i)*Xreal(i) 
                       + 1.5*epsilon*invCquad(i)*Xreal(i)*dXreal(i))
                - mu/6.0*C5Quad(i)*d3Xreal(i);
          }
      }
  }
  
}

#define MONTJOIE_FILE_KDV_PROBLEM_CXX
#endif
