#ifndef MONTJOIE_FILE_TIME_MAXWELL3D_CXX

#include "ProdMatVectMaxwell3D.cxx"

namespace Montjoie
{
  
  
  //! initialization before time iterations
  template<class TypeEquation>
  void TimeMaxwell_3D<TypeEquation>::InitTimeIterations()
  {    
    VarInstationary<TypeEquation>::InitTimeIterations();
    
    int nb_dofV = this->GetNbVectorialUnknowns();
    int nb_dofU = this->GetNbScalarUnknowns();
    
    // DISP(nb_dofU); DISP(nb_dofV);
    VectReal_wp En(nb_dofU), Hn(nb_dofV);
    En.Fill(0); Hn.Fill(0);
    
    GetExtrapolationMaxwell3D(En, Hn, -1, true, true);

    
  }
  

  template<class TypeEquation>
  void TimeMaxwell_3D<TypeEquation>::RunTimeIterations()
  {
    VarInstationary<TypeEquation>::RunTimeIterations();
  }
  

  template<class TypeEquation>
  void TimeMaxwell_3D<TypeEquation>
  ::GetExtrapolationMaxwell3D(const VectReal_wp& U, const VectReal_wp& V, int lvl, 
			      bool extrapol_u, bool extrapol_v)
  {
  }
  
  
  /***********************
   * Continuous elements *
   ***********************/
  
  
  void HyperbolicProblem<TimeMaxwellEquation_3D>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X, VectReal_wp& Y,
                             bool invert_mass, bool source)
  {
    EvaluateDerivativeFunction(tn, nb_deriv, 1.0, MatrixVectorProductLevel::ALL_LEVELS,
                               X, Y, invert_mass, source);
  }
  
  
  //! compute Y = G(tn,X) if we put time-maxwell equation to the form dX/dt = G(t,X)
  /*! 
    \param[in] tn input time
    \param[in] nb_deriv derivative number of G,
    if 0 we evaluate G, if 1 dG/dt, if 2 d^2 G /dt^2 ...
    \param[in] X iterate
    \param[out] Y result of G(tn,X)
  */
  void HyperbolicProblem<TimeMaxwellEquation_3D>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                             const Real_wp& alpha, int level,
                             const VectReal_wp& X, VectReal_wp& Y,
                             bool invert_mass, bool source)
  {
    EllipticProblem<StaticMaxwellEquation_3D>& var = this->var_harmonic;

    int nb_dof_E = var.offset_dof_unknown(1);
    int nb_dof_H = var.GetNbDof() - nb_dof_E;
    Real_wp one(1), zero(0);

    VectReal_wp En, Hn, Prod_En, Prod_Hn;
    Real_wp* x_ptr = X.GetData(); Real_wp* y_ptr = Y.GetData();
    En.SetData(nb_dof_E, x_ptr); Prod_En.SetData(nb_dof_E, y_ptr);
    x_ptr += nb_dof_E; y_ptr += nb_dof_E;
    
    Hn.SetData(nb_dof_H, x_ptr); Prod_Hn.SetData(nb_dof_H, y_ptr);
    
    // dE/dt = Dh^{-1} ( Rh^S H + Fh - Sh E)
    this->MltAddStiffnessScalar(one, level, Hn, zero, Prod_En);

    this->ApplyOperatorSh(-one, tn, En, one, Prod_En);
    if (source)
      {
	this->AddScalarSourceAtTime(one, tn, nb_deriv, Prod_En);
	this->SetDirichletConditionSource(tn, nb_deriv+1, Prod_En);
      }
    else
      var.ImposeNullDirichletCondition(Prod_En);

    if (invert_mass)
      this->SolveOperatorDh(Prod_En);
    
    // dH/dt = Bh^{-1} ( Rh^V E - ShV H + FhV)
    this->MltAddStiffnessVectorial(one, level, En, zero, Prod_Hn);
    if (source)
      this->AddVectorialSourceAtTime(one, tn, nb_deriv, Prod_Hn);
    
    this->ApplyOperatorShVectorial(-one, tn, Hn, one, Prod_Hn);
    if (invert_mass)
      this->SolveOperatorBh(Prod_Hn);    
    
    // nullify temporary vectors
    En.Nullify(); Hn.Nullify(); Prod_En.Nullify(); Prod_Hn.Nullify();
  }


  void HyperbolicProblem<TimeMaxwellEquation_3D>
  ::MltAddStiffnessScalar(const Real_wp& alpha, int level, const VectReal_wp& B,
			  const Real_wp& beta, VectReal_wp& C)
  {
    MltAdd_SquareHex_ScalarHcurl3D(alpha, level, *this, B, beta, C);
  }


  void HyperbolicProblem<TimeMaxwellEquation_3D>
  ::MltAddStiffnessVectorial(const Real_wp& alpha, int level, const VectReal_wp& B,
			     const Real_wp& beta, VectReal_wp& C)
  {
    MltAdd_SquareHex_VectorialHcurl3D(alpha, level, *this, B, beta, C);
  }


  void HyperbolicProblem<TimeMaxwellEquation_3D>
  ::AddElementStiffnessTensor(int i, int& offset, const Real_wp& alpha_,
                              const VectReal_wp& Uh, VectReal_wp& Vh,
                              const Real_wp& coef_mu, const Real_wp& coef_sigma)
  {
    const ElementReference<Dimension3, 2>& Fb = this->var_harmonic.GetReferenceElementHcurl(i);
    int ref = this->var_harmonic.mesh.Element(i).GetReference();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    const VectReal_wp& WeightsND = Fb.WeightsND();
    Real_wp coef; R3 sigma, epsilon;
    bool affine = this->var_harmonic.mesh.IsElementAffine(i);
    Real_wp alpha = coef_mu*alpha_;
    if (coef_mu != Real_wp(0))
      for (int j = 0; j < nb_points_quad; j++)
        {
          if (affine)
            coef = this->var_harmonic.Glob_jacobian(i)(0)*WeightsND(j)*alpha;
          else
            coef = this->var_harmonic.Glob_jacobian(i)(j)*alpha;
          
          sigma(0) = Uh(offset);
          sigma(1) = Uh(offset+1);
          sigma(2) = Uh(offset+2);
          
          this->var_harmonic.ref_mu(ref).MltMatrix(this->var_harmonic, i, j, sigma, epsilon);
          
          Vh(offset) += coef*epsilon(0);
          Vh(offset+1) += coef*epsilon(1);
          Vh(offset+2) += coef*epsilon(2);
          
          offset += 3;
        }
    else
      offset += 3*nb_points_quad;
    
    int offset2 = offset - 3*nb_points_quad;
    if (this->var_harmonic.InsidePML(i))
      {
        int i1 = i - this->var_harmonic.mesh.GetNbElt() + this->var_harmonic.GetNbEltPML();
        R3 tau;
        Real_wp mu = this->var_harmonic.ref_mu(ref).GetConstant()(0, 0), beta;
        for (int j = 0; j < nb_points_quad; j++)
          {
            if (affine)
              coef = this->var_harmonic.Glob_jacobian(i)(0)*WeightsND(j)*alpha_;
            else
              coef = this->var_harmonic.Glob_jacobian(i)(j)*alpha_;
            
            alpha = coef_mu*coef;
            Vh(offset) += alpha*(Uh(offset) - Uh(offset2));
            Vh(offset+1) += alpha*(Uh(offset+1) - Uh(offset2+1));
            Vh(offset+2) += alpha*(Uh(offset+2) - Uh(offset2+2));
            
            if (coef_sigma != Real_wp(0))
              {
                alpha = coef_sigma*coef;
                beta = alpha*mu;
                tau = this->var_harmonic.GetTauPML(i1, j); 
                Vh(offset2) += beta*tau(1)*Uh(offset2);
                Vh(offset) += alpha*(tau(2)*Uh(offset) - tau(0)*Uh(offset2));

                Vh(offset2+1) += beta*tau(2)*Uh(offset2+1);
                Vh(offset+1) += alpha*(tau(0)*Uh(offset+1) - tau(1)*Uh(offset2+1));
                
                Vh(offset2+2) += beta*tau(0)*Uh(offset2+2);
                Vh(offset+2) += alpha*(tau(1)*Uh(offset+2) - tau(2)*Uh(offset2+2));
              }
            
            offset2 += 3;
            offset += 3;
          }
      }
  }
    
  
  void HyperbolicProblem<TimeMaxwellEquation_3D>
  ::SolveElementStiffnessTensor(int i, int& offset, VectReal_wp& Vh,
                                const Real_wp& coef_sigma)
  {
    SolveElementStiffnessTensorGen(i, offset, Real_wp(1), Vh, Vh, coef_sigma, false);
  }


  void HyperbolicProblem<TimeMaxwellEquation_3D>
  ::SolveElementStiffnessTensorGen(int i, int& offset, const Real_wp& alpha,
                                   const VectReal_wp& Uh, VectReal_wp& Vh,
                                   const Real_wp& coef_sigma, bool add)
  {
    const ElementReference<Dimension3, 2>& Fb = this->var_harmonic.GetReferenceElementHcurl(i);
    int ref = this->var_harmonic.mesh.Element(i).GetReference();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    const VectReal_wp& WeightsND = Fb.WeightsND();
    Real_wp coef; R3 sigma, epsilon, tau;
    bool affine = this->var_harmonic.mesh.IsElementAffine(i);
    int i1 = i - this->var_harmonic.mesh.GetNbElt() + this->var_harmonic.GetNbEltPML();
    for (int j = 0; j < nb_points_quad; j++)
      {
        if (affine)
          coef = alpha/(this->var_harmonic.Glob_jacobian(i)(0)*WeightsND(j));
        else
          coef = alpha/this->var_harmonic.Glob_jacobian(i)(j);
        
        sigma(0) = Uh(offset);
        sigma(1) = Uh(offset+1);
        sigma(2) = Uh(offset+2);
        
        this->var_harmonic.ref_invMu(ref).MltMatrix(this->var_harmonic, i, j, sigma, epsilon);
        
        if (this->var_harmonic.InsidePML(i))
          {
            int offset2 = offset + 3*nb_points_quad;

            if (coef_sigma == Real_wp(0))
              {
                if (add)
                  { 
                    Vh(offset2) += coef*(Uh(offset2) + epsilon(0));
                    Vh(offset2+1) += coef*(Uh(offset2+1) + epsilon(1));
                    Vh(offset2+2) += coef*(Uh(offset2+2) + epsilon(2));
                  }
                else
                  {
                    Vh(offset2) = coef*(Uh(offset2) + epsilon(0));
                    Vh(offset2+1) = coef*(Uh(offset2+1) + epsilon(1));
                    Vh(offset2+2) = coef*(Uh(offset2+2) + epsilon(2));
                  }
              }
            else
              {
                tau = this->var_harmonic.GetTauPML(i1, j); 
                
                Real_wp d0 = 1.0/(1.0 + coef_sigma*tau(0));
                Real_wp d1 = 1.0/(1.0 + coef_sigma*tau(1));
                Real_wp d2 = 1.0/(1.0 + coef_sigma*tau(2));
                
                epsilon(0) *= d1; 
                epsilon(1) *= d2;
                epsilon(2) *= d0;
                
                if (add)
                  { 
                    Vh(offset2) += coef*d2*(Uh(offset2) + (1.0+tau(0)*coef_sigma)*epsilon(0));
                    Vh(offset2+1) += coef*d0*(Uh(offset2+1) + (1.0+tau(1)*coef_sigma)*epsilon(1));
                    Vh(offset2+2) += coef*d1*(Uh(offset2+2) + (1.0+tau(2)*coef_sigma)*epsilon(2));
                  }
                else
                  { 
                    Vh(offset2) = coef*d2*(Uh(offset2) + (1.0+tau(0)*coef_sigma)*epsilon(0));
                    Vh(offset2+1) = coef*d0*(Uh(offset2+1) + (1.0+tau(1)*coef_sigma)*epsilon(1));
                    Vh(offset2+2) = coef*d1*(Uh(offset2+2) + (1.0+tau(2)*coef_sigma)*epsilon(2));
                  }
              }
          }
        
        if (add)
          {
            Vh(offset) += coef*epsilon(0);
            Vh(offset+1) += coef*epsilon(1);
            Vh(offset+2) += coef*epsilon(2);
          }
        else
          {
            Vh(offset) = coef*epsilon(0);
            Vh(offset+1) = coef*epsilon(1);
            Vh(offset+2) = coef*epsilon(2);
          }
        
        offset += 3;
      }

    if (this->var_harmonic.InsidePML(i))
      offset += 3*nb_points_quad;
  }
  
  
  void HyperbolicProblem<TimeMaxwellEquation_3D>
  ::SolveAddElementStiffnessTensor(int i, int& offset, const Real_wp& alpha,
                                   const VectReal_wp& Uh, VectReal_wp& Vh,
                                   const Real_wp& coef_sigma)
  {
    SolveElementStiffnessTensorGen(i, offset, alpha, Uh, Vh, coef_sigma, true);
  }
  
  
  /**************************
   * Discontinuous elements *
   **************************/
  

  void HyperbolicProblem<TimeMaxwellEquation_3D_DG>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X, VectReal_wp& Y,
                             bool invert_mass, bool source)
  {
    EvaluateDerivativeFunction(tn, nb_deriv, 1.0, MatrixVectorProductLevel::ALL_LEVELS,
                               X, Y, invert_mass, source);
  }
  
    
  //! compute Y = G(tn,X) if we put time-maxwell equation to the form dX/dt = G(t,X)
  /*! 
    \param[in] tn input time
    \param[in] nb_deriv derivative number of G,
    if 0 we evaluate G, if 1 dG/dt, if 2 d^2 G /dt^2 ...
    \param[in] X iterate
    \param[out] Y result of G(tn,X)
  */
  void HyperbolicProblem<TimeMaxwellEquation_3D_DG>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                             const Real_wp& alpha, int level,
                             const VectReal_wp& X, VectReal_wp& Y, bool invert_mass,
                             bool source)
  {
    EllipticProblem<StaticMaxwellEquation_3D_DG>& var = this->var_harmonic;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int nb_unknowns = 3*mesh_num.GetNbDof();
    int nb_unknowns_pml = 3*mesh_num.GetNbDofPML();
    
    VectReal_wp En, Hn, En_star, Hn_star;
    VectReal_wp Prod_En, Prod_Hn, Prod_En_star, Prod_Hn_star;
    
    Real_wp one(1), zero(0);
    
    // X is split in (Uh, Vh, Vh_pml, Vh_diamond, Vh_star)
    // idem for Y
    Real_wp* x_ptr = X.GetData(); Real_wp* y_ptr = Y.GetData();
    En.SetData(nb_unknowns, x_ptr); Prod_En.SetData(nb_unknowns, y_ptr);
    x_ptr += nb_unknowns; y_ptr += nb_unknowns;
    
    Hn.SetData(nb_unknowns, x_ptr); Prod_Hn.SetData(nb_unknowns, y_ptr);
    x_ptr += nb_unknowns; y_ptr += nb_unknowns;
    if (nb_unknowns_pml > 0)
      {
	En_star.SetData(nb_unknowns_pml, x_ptr); Prod_En_star.SetData(nb_unknowns_pml, y_ptr); 
	x_ptr += nb_unknowns_pml; y_ptr += nb_unknowns_pml;
	
	Hn_star.SetData(nb_unknowns_pml,x_ptr); Prod_Hn_star.SetData(nb_unknowns_pml,y_ptr);
	x_ptr += nb_unknowns_pml; y_ptr += nb_unknowns_pml;
      }
    
    GetExtrapolationMaxwell3D(En, Hn, level, true, true);
    
    // dE/dt = Dh^{-1} ( Rh^S H + Fh - Sh E)
    this->ApplyOperatorRhScalar(one, tn, Hn, zero, Prod_En, false);
    if (source)
      this->AddScalarSourceAtTime(one, tn, nb_deriv, Prod_En);
    //this->ApplyOperatorSh(-one, tn, En, one, Prod_En);
    
    // dH/dt = Bh^{-1} ( Rh^V E - ShV H + FhV) 
    this->ApplyOperatorRhVectorial(one, tn, En, zero, Prod_Hn, false);
    if (source)
      this->AddVectorialSourceAtTime(one, tn, nb_deriv, Prod_Hn);
    //this->ApplyOperatorShVectorial(-one, tn, Hn, one, Prod_Hn);

    if (invert_mass)
      {
        this->SolveOperatorDh(Prod_En);
        this->SolveOperatorBh(Prod_Hn);
      }
    
    glob_chrono.Start(VirtualTimer::PML);
    if (nb_unknowns_pml > 0)
      {
        int offsetx_pml = 0; int Nvol = mesh_num.GetNbDof();        
        for (int i = 0; i < var.GetNbEltPML(); i++)
          {
            int i1 = i + var.mesh.GetNbElt() - var.GetNbEltPML();
            int ref = var.mesh.Element(i1).GetReference();
            Real_wp mu = var.ref_mu(ref)(0, 0);
            Real_wp epsilon = var.ref_epsilon(ref)(0, 0);
            Real_wp c = 1.0/sqrt(mu*epsilon);
            int nb_pts_quad = var.GetNbPointsQuadratureInside(i1);
            R3 tau;
            int offset_x = mesh_num.Element(i1).GetNumberDof(0);
            int offset_y = offset_x + Nvol;
            int offset_z = offset_y + Nvol;
            int offsety_pml = offsetx_pml + nb_unknowns_pml/3;
            int offsetz_pml = offsety_pml + nb_unknowns_pml/3;
            for (int j = 0; j < nb_pts_quad; j++)
              {
                //Real_wp invMass = 1.0/(jacob*WeightsND(j));
                tau(0) = c*var.GetTauPML(i, j, 0); 
                tau(1) = c*var.GetTauPML(i, j, 1); 
                tau(2) = c*var.GetTauPML(i, j, 2);                
                
                Prod_En_star(offsetx_pml + j)
                  = Prod_En(offset_x + j) - tau(1)*En_star(offsetx_pml + j);
                Prod_En_star(offsety_pml + j)
                  = Prod_En(offset_y + j) - tau(2)*En_star(offsety_pml + j);
                Prod_En_star(offsetz_pml + j)
                  = Prod_En(offset_z + j) - tau(0)*En_star(offsetz_pml + j);
                
                Prod_Hn_star(offsetx_pml + j)
                  = Prod_Hn(offset_x + j) - tau(1)*Hn_star(offsetx_pml + j);
                Prod_Hn_star(offsety_pml + j)
                  = Prod_Hn(offset_y + j) - tau(2)*Hn_star(offsety_pml + j);
                Prod_Hn_star(offsetz_pml + j)
                  = Prod_Hn(offset_z + j) - tau(0)*Hn_star(offsetz_pml + j);
                
                Prod_En(offset_x + j) = Prod_En_star(offsetx_pml+j)
                  + tau(0)*En_star(offsetx_pml+j) - tau(2)*En(offset_x+j);
                Prod_En(offset_y + j) = Prod_En_star(offsety_pml+j)
                  + tau(1)*En_star(offsety_pml+j) - tau(0)*En(offset_y+j);
                Prod_En(offset_z + j) = Prod_En_star(offsetz_pml+j)
                  + tau(2)*En_star(offsetz_pml+j) - tau(1)*En(offset_z+j);
                
                Prod_Hn(offset_x + j) = Prod_Hn_star(offsetx_pml+j)
                  + tau(0)*Hn_star(offsetx_pml+j) - tau(2)*Hn(offset_x+j);
                Prod_Hn(offset_y + j) = Prod_Hn_star(offsety_pml+j)
                  + tau(1)*Hn_star(offsety_pml+j) - tau(0)*Hn(offset_y+j);
                Prod_Hn(offset_z + j) = Prod_Hn_star(offsetz_pml+j)
                  + tau(2)*Hn_star(offsetz_pml+j) - tau(1)*Hn(offset_z+j);
              }
            
            offsetx_pml += nb_pts_quad;
          }
      }

    glob_chrono.Stop(VirtualTimer::PML);
    
    // nullify temporary vectors
    En.Nullify(); Hn.Nullify(); Prod_En.Nullify(); Prod_Hn.Nullify();
    if (nb_unknowns_pml > 0)
      {
	Hn_star.Nullify(); En_star.Nullify();
	Prod_Hn_star.Nullify(); Prod_En_star.Nullify();
      }
  }
  
  
  void HyperbolicProblem<TimeMaxwellEquation_3D_DG>
  ::ComputeExtrapolationVectorial(const VectReal_wp& En, const VectReal_wp& Hn)
  {
    int level = this->global_level_operator;
    this->GetExtrapolationMaxwell3D(En, Hn, level, false, true);
  }
  
  
  void HyperbolicProblem<TimeMaxwellEquation_3D_DG>
  ::ComputeExtrapolationScalar(const VectReal_wp& En, const VectReal_wp& Hn)
  {
    int level = this->global_level_operator;
    this->GetExtrapolationMaxwell3D(En, Hn, level, true, false);    
  }   


  void HyperbolicProblem<TimeMaxwellEquation_3D_DG>
  ::GetExtrapolationMaxwell3D(const VectReal_wp& U, const VectReal_wp& V, int lvl, 
			      bool extrapol_u, bool extrapol_v)
  {
    Montjoie::GetExtrapolationMaxwell3D(this->GetLeafClass(), U, V, lvl, 
					this->extrapolE, this->extrapolH, this->Eneighbor, this->Hneighbor,
					this->evalEn_quad, this->evalHn_quad, extrapol_u, extrapol_v);
  }

  
  void HyperbolicProblem<TimeMaxwellEquation_3D_DG>
  ::MltAddStiffnessScalar(const Real_wp& alpha, int level, const VectReal_wp& B,
			  const Real_wp& beta, VectReal_wp& C)
  {
    MltAddStiffnessScalar_MaxwellDG(alpha, level, *this, B, beta, C);
  }

  
  void HyperbolicProblem<TimeMaxwellEquation_3D_DG>
  ::MltAddStiffnessVectorial(const Real_wp& alpha, int level, const VectReal_wp& B,
			     const Real_wp& beta, VectReal_wp& C)
  {
    MltAddStiffnessVectorial_MaxwellDG(alpha, level, *this, B, beta, C);
  }


  void HyperbolicProblem<TimeMaxwellEquation_3D_DG>
  ::GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp& rho, Real_wp& sigma)
  {
    // isotropic case only
    if (num <= 2)
      {
	rho = this->var_harmonic.ref_epsilon(ref).GetCoefficient(this->var_harmonic, i, j, num, num);
	sigma = this->var_harmonic.ref_sigma(ref).GetCoefficient(this->var_harmonic, i, j, num, num);
      }
    else
      {
	rho = this->var_harmonic.ref_mu(ref).GetCoefficient(this->var_harmonic, i, j, num-3, num-3);
	sigma = 0.0;
      }
  }


  /*****************************
   * TimeMaxwellEquationHdg_3D *
   *****************************/
  
  
  void HyperbolicProblem<TimeMaxwellEquationHdg_3D>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X, VectReal_wp& Y,
                             bool invert_mass, bool source)
  {
    EllipticProblem<StaticMaxwellEquationHdg_3D>& var = this->var_harmonic;

    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int nb_dof_E = 0, nb_dof_H = 0, N = 0;
    int Nvol = var.GetNbMainUnknownDof();
    nb_dof_E = Nvol;
    nb_dof_H = Nvol;
    N = nb_dof_E + nb_dof_H;
    
    Real_wp one(1);    
    VectReal_wp En, Hn;
    VectReal_wp Prod_En, Prod_Hn;
    
    if ((N > X.GetM()) || (N > Y.GetM()))
      {
        cout << "X or Y is not large enough" << endl;
        DISP(N);
        DISP(X.GetM()); DISP(Y.GetM());
        abort();
      }

    // we retrieve E and H
    Real_wp* x_ptr = X.GetData(); Real_wp* y_ptr = Y.GetData();
    En.SetData(nb_dof_E, x_ptr); Prod_En.SetData(nb_dof_E, y_ptr);
    x_ptr += nb_dof_E; y_ptr += nb_dof_E;
    
    Hn.SetData(nb_dof_H, x_ptr); Prod_Hn.SetData(nb_dof_H, y_ptr);
    int level = MatrixVectorProductLevel::ALL_LEVELS;

    /* cout << "Je suis dans EvaluateDerivative" << endl;
    VectReal_wp& Xb = const_cast<VectReal_wp&>(X);
    Xb.FillRand(); Xb *= 1e-9;

    Xb.WriteText("X.dat"); */
    
    // HDG case
    int nb_dof_L = mesh_num.GetNbDof();
    VectReal_wp& Lambda = evalLambda_n;
    Lambda.Reallocate(nb_dof_L);	
    Y.Zero();
    
    if (this->Glob_mat_RhS.GetM() > 0)
      MltAdd(-Real_wp(1), this->Glob_mat_RhS, X, Real_wp(0), Lambda);
    else
      this->GetExtrapolationLambda(En, Hn, level, Lambda);

    //Lambda.WriteText("ProdLambda.dat");
    
    if (source)
      this->AddTimeSourceHDG(Real_wp(1), tn, nb_deriv, Y, Lambda);

    this->Glob_mat_Dh->SolveOperatorCh(Lambda);

    //Lambda.WriteText("Lambda.dat");
    //Y.WriteText("Y.dat");
    
    // testing MltAddStiffnessScalar/Vectorial
    if (this->Glob_mat_Ch_Lambda.GetM() > 0)
      {
	MltAdd(Real_wp(-1), this->Glob_mat_Ch_Lambda, Lambda, one, Y);
	MltAdd(Real_wp(-1), this->Glob_mat_RhV, X, one, Y);
      }
    else
      this->MltAddStiffness(one, level, level, En, Hn, Lambda, one, Prod_En, Prod_Hn);

    //Prod_En.WriteText("ProdE.dat");
    //Prod_Hn.WriteText("ProdH.dat");
    //Y.WriteText("ProdY.dat");
	
    this->ApplyOperatorSh(-one, tn, En, one, Prod_En);

    //Prod_En.WriteText("ProdE2.dat");
	
    if (invert_mass)
      this->SolveMassMatrix(Y);

    //Y.WriteText("AY.dat");
    //int test_input; cout << "waiting" << endl; cin >> test_input;
    
    // nullify temporary vectors
    En.Nullify(); Hn.Nullify(); Prod_En.Nullify(); Prod_Hn.Nullify();
  }
  
  
  //! compute Y = G(tn,X) if we put time-maxwell equation to the form dX/dt = G(t,X)
  /*! 
    \param[in] tn input time
    \param[in] nb_deriv derivative number of G,
    if 0 we evaluate G, if 1 dG/dt, if 2 d^2 G /dt^2 ...
    \param[in] X iterate
    \param[out] Y result of G(tn,X)
  */
  void HyperbolicProblem<TimeMaxwellEquationHdg_3D>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                             const Real_wp& alpha, int level,
                             const VectReal_wp& X, VectReal_wp& Y, bool invert_mass,
                             bool source)
  {
    cout << "Not implemented" << endl;
    abort();
  }


  void HyperbolicProblem<TimeMaxwellEquationHdg_3D>
  ::GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp& rho, Real_wp& sigma)
  {
    if (num == 0)
      {
	rho = this->var_harmonic.ref_epsilon(ref).GetCoefficient(this->var_harmonic, i, j, 0, 0);
	sigma = this->var_harmonic.ref_sigma(ref).GetCoefficient(this->var_harmonic, i, j, 0, 0);
      }
    else
      {
	rho = this->var_harmonic.ref_mu(ref).GetCoefficient(this->var_harmonic, i, j, 0, 0);
	sigma = 0.0;
      }
  }
  
  
  Real_wp HyperbolicProblem<TimeMaxwellEquationHdg_3D>
  ::GetCoefficientTauHDG(int num_face) const
  {
    return this->var_harmonic.alpha_penalization*this->var_harmonic.Glob_CoefPenalDG(num_face);
  }

  
  void HyperbolicProblem<TimeMaxwellEquationHdg_3D>
  ::AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t, int nb_deriv, VectReal_wp& b_src)
  {
    if ((this->var_harmonic.FormulationDG() == ElementReference_Base::HDG)
	&& (!this->SourceDoesNotDependOnTime()))
      {
	if (t >= this->tlimit_source)
	  return;

	int nb_dof_L = this->var_harmonic.GetNbMeshDof();
	VectReal_wp Lambda(nb_dof_L);
	Lambda.Zero();
	this->AddTimeSourceHDG(alpha, t, nb_deriv, b_src, Lambda);

	this->Glob_mat_Dh->SolveOperatorCh(Lambda);

	if (this->Glob_mat_Ch_Lambda.GetM() > 0)
	  MltAdd(Real_wp(-1), this->Glob_mat_Ch_Lambda, Lambda, Real_wp(1), b_src);
	else
	  {
	    cout << "not implemented" << endl;
	    abort();
	    //MltAddLambdaAcousticHDG(Real_wp(1), -1, *this, Lambda, b_src);
	  }
      }
    else
      TimeMaxwell_3D<TimeMaxwellEquationHdg_3D>::AddPrimitiveTimeSource(alpha, t, nb_deriv, b_src);
  }

  
  void HyperbolicProblem<TimeMaxwellEquationHdg_3D>
  ::MltAddStiffness(const Real_wp& alpha, int level, int level2,
		    const VectReal_wp& Uh, const VectReal_wp& Vh,
		    const VectReal_wp& Lambda, const Real_wp& beta,
		    VectReal_wp& ProdUh, VectReal_wp& ProdVh)
  {
    MltAddStiffnessMaxwellHDG(alpha, level, level, *this, Uh, Vh, Lambda,
			      this->evalEn_quad, this->evalHn_quad,
			      this->extrapolE, this->extrapolH, beta, ProdUh, ProdVh);
  }
    

  void HyperbolicProblem<TimeMaxwellEquationHdg_3D>
  ::GetExtrapolationLambda(const VectReal_wp& Uh, const VectReal_wp& Vh, int level, VectReal_wp& Lambda)
  {
    GetExtrapolationMaxwellHDG(*this, level, Uh, Vh, this->evalEn_quad, this->evalHn_quad,
    				this->extrapolE, this->extrapolH, Lambda);
  }
  

  void HyperbolicProblem<TimeMaxwellEquationHdg_3D>
  ::AddTimeSourceHDG(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& Y, VectReal_wp& Lambda)
  {
    EllipticProblem<StaticMaxwellEquationHdg_3D>& var = this->var_harmonic;
    if (t >= this->tlimit_source)
      return;
    
    if (this->SourceDoesNotDependOnTime())
      this->AddPrimitiveSourceAtTime(alpha, t, n, Y);
    else
      {
	VectReal_wp b_src(var.GetNbDof());
	b_src.Zero();
	
	Vector<VectReal_wp> b_vec;
        b_vec.SetData(1, &b_src);
        
	VirtualSourceFEM<Real_wp, Dimension>* f;
	f = var.GetNewSourceEquationObject(0);
	
	f->Init(t, this->GetTimeStep(), var.print_level, n, true);
	
	Vector<VirtualSourceFEM<Real_wp, Dimension>* > f_vec(1);
	f_vec(0) = f;
	var.AddSurfaceSource(Real_wp(1), b_vec, f_vec);
        b_vec.Nullify();
	
	int nb_dof_L = var.nb_unknowns_hdg*var.GetNbMeshDof();
	int nb_dof_uv = var.GetNbDof() - nb_dof_L;
	
	for (int i = 0; i < nb_dof_L; i++)
	  Lambda(i) += alpha*b_src(i);

	for (int i = 0; i < nb_dof_uv; i++)
	  Y(i) += alpha*b_src(nb_dof_L + i);
      }
  }


  
  VirtualMassMatrix* HyperbolicProblem<TimeMaxwellEquationHdg_3D>::GetNewMassMatrix()
  {
    return new DiscontinuousMassMatrixHdgMaxwell3D(this->GetLeafClass());
  }


  /***************************************
   * DiscontinuousMassMatrixHdgMaxwell3D *
   ***************************************/

  
  //! constructor with a given problem
  template<class TypeEquation>
  DiscontinuousMassMatrixHdgMaxwell3D
  ::DiscontinuousMassMatrixHdgMaxwell3D(HyperbolicProblem<TypeEquation>& var)
    : DiscontinuousUnsteadyMassMatrix<Dimension3>(var),
    var_problem(var.var_harmonic), var_time(var), var_boundary(var.var_harmonic),
    var_maxwell(var.var_harmonic)
  {
    type_matrix = 0;
  }
    

  //! initialisation of mass matrix
  void DiscontinuousMassMatrixHdgMaxwell3D::Init(bool compute_time)
  {
    Vector<bool> diag_elt;
    int type = var_problem.GetMassMatrixType(diag_elt);    
    if (type == FemMassMatrix::BLOCK_DIAGONAL_UNSYM)
      type = FemMassMatrix::BLOCK_DIAGONAL;

    type_matrix = type;
    ComputeOperatorCh();
    ComputeOperatorBhGeom();

    // checking which matrices have to be computed
    bool compute_Dh, compute_invDh, compute_DhMinusdtSh,
      compute_invDhPlusdtSh, compute_Sh, compute_Bh, compute_invBh,
      compute_BhMinusdtSh, compute_invBhPlusdtSh, compute_ShVec;
    
    this->FindMatricesToCompute(compute_Dh, compute_invDh, compute_DhMinusdtSh,
				compute_invDhPlusdtSh, compute_Sh);
    
    this->FindMatricesToComputeVec(compute_Bh, compute_invBh, compute_BhMinusdtSh,
				   compute_invBhPlusdtSh, compute_ShVec);

    /* int Nvol = var_problem.GetNbMainUnknownDof();
    VectReal_wp X(Nvol), Y(Nvol);
    X.FillRand(); X *= 1e-9; X.Write("X.dat"); */

    const Mesh<Dimension3>& mesh = var_problem.mesh;
    VectReal_wp coef_block(mesh.GetNbElt());
    if (compute_Dh || compute_invDh)
      {
	this->Dh = this->GetNewMassMatrix(FemMassMatrix::BLOCK_ELT);
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    int ref = var_problem.mesh.Element(i).GetReference();
	    coef_block(i) = var_maxwell.ref_epsilon(ref).GetConstant()(0, 0);
	  }
	
	this->Dh->CopyBlockCoef(coef_block, 0);

	// testing ApplyOperatorDh
	/*
	this->ApplyOperatorDh(Real_wp(1), Real_wp(0), X, Real_wp(0), Y);
	
	Y.Write("Y.dat");
	
	Y.Zero(); Y = X;
	this->SolveOperatorDh(Y);
	
	Y.Write("Z.dat"); */
      }

    if (compute_Sh)
      {
	this->Sh = this->GetNewMassMatrix(FemMassMatrix::BLOCK_ELT);
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    int ref = var_problem.mesh.Element(i).GetReference();
	    coef_block(i) = var_maxwell.ref_sigma(ref).GetConstant()(0, 0);
	  }
	
	this->Sh->CopyBlockCoef(coef_block, 0);

	// testing ApplyOperatorSh
	/*this->ApplyOperatorSh(Real_wp(1), Real_wp(0), X, Real_wp(0), Y);
	
	  Y.Write("Y.dat"); */
      }

    if (compute_DhMinusdtSh || compute_invDhPlusdtSh ||
	compute_BhMinusdtSh || compute_invBhPlusdtSh )
      {
	cout << "Scheme not adapted for HDG" << endl;
	abort();
      }
    
    if (compute_Bh || compute_invBh)
      {
	this->Bh = this->GetNewMassMatrix(FemMassMatrix::BLOCK_ELT);
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    int ref = var_problem.mesh.Element(i).GetReference();
	    coef_block(i) = var_maxwell.ref_mu(ref).GetConstant()(0, 0);
	  }
	
	this->Bh->CopyBlockCoef(coef_block, 1);

	// testing ApplyOperatorDh
	/* this->ApplyOperatorBh(Real_wp(1), Real_wp(0), X, Real_wp(0), Y);
	
	Y.Write("Y.dat");
	
	Y.Zero(); Y = X;
	this->SolveOperatorBh(Y);
	
	Y.Write("Z.dat"); */
      }
  }


  //! computes operator Bh
  void DiscontinuousMassMatrixHdgMaxwell3D::ComputeOperatorBhGeom()
  {
    const Mesh<Dimension3>& mesh = var_problem.mesh;
    this->GeomMass = GetNewMassMatrix(FemMassMatrix::BLOCK_ELT);
    this->GeomMass->ReallocateStoredMass();

    // checking that the indices are constant and isotropic
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	int ref = var_problem.mesh.Element(i).GetReference();
	if (var_maxwell.ref_epsilon(ref).IsVarying() ||
	    (var_maxwell.ref_epsilon(ref).GetAnisotropy() != TensorPhysicalIndice<Dimension3, 3, Real_wp>::ISOTROPE) || var_maxwell.ref_sigma(ref).IsVarying() ||
	    (var_maxwell.ref_sigma(ref).GetAnisotropy() != TensorPhysicalIndice<Dimension3, 3, Real_wp>::ISOTROPE) || var_maxwell.ref_mu(ref).IsVarying() ||
	    (var_maxwell.ref_mu(ref).GetAnisotropy() != TensorPhysicalIndice<Dimension3, 3, Real_wp>::ISOTROPE))
	  {
	    cout << "Case not handled" << endl;
	    abort();
	  }	
      }
    
    DiscontinuousBlockMassMatrix* block;
    Vector<int> num_row; VectReal_wp diagonal;
    Vector<TinyMatrix<Real_wp, Symmetric, 3, 3> > Bh_b;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        HarmonicMaxwell3D_PhysGeomInfo<Real_wp>& mass = var_maxwell.Glob_matMass_elem(i);
        const ElementReference<Dimension3, 2>& Fb = var_problem.GetReferenceElementHcurl(i);
        const HexahedronHcurlLobatto& Fb_hex = dynamic_cast<const HexahedronHcurlLobatto&>(Fb);
	
        int ref = var_problem.mesh.Element(i).GetReference();
        Real_wp invEpsilon = Real_wp(1)/var_maxwell.ref_epsilon(ref).GetConstant()(0, 0);
	
        int nb_points_quad = Fb.GetNbPointsQuadratureInside();
        if (mass.IsDiagonalMass())
          {
            diagonal.Reallocate(Fb.GetNbDof());
            for (int j = 0; j < Fb.GetNbDof(); j++)
              diagonal(j) = invEpsilon*mass.GetDh(j);

            block = new DiscontinuousDiagonalMassMatrix(diagonal);
          }
        else
          {
            num_row.Reallocate(3*nb_points_quad);
            Bh_b.Reallocate(nb_points_quad);
            for (int j = 0; j < nb_points_quad; j++)
              {
                int jx, jy, jz;
                Fb_hex.GetDofNumber_FromPointNode(j, jx, jy, jz);
                
                num_row(3*j) = jx; num_row(3*j+1) = jy; num_row(3*j+2) = jz;
                Bh_b(j) = invEpsilon*mass.GetBh(j);
              }
            
            block = new DiscontinuousBlockDiagMassMatrix<3>(num_row, Bh_b);
          }
        
        this->GeomMass->SetLocalStoredBlock(i, block);
      }
  }
  

  //! computes operator Ch
  void DiscontinuousMassMatrixHdgMaxwell3D::ComputeOperatorCh()
  {
    this->InitSurface(false);
    
    const Mesh<Dimension3>& mesh = var_problem.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    
    if (type_matrix == FemMassMatrix::DIAGONAL)
      {
	// extracting Ch (mass matrix for lambda)
	VectReal_wp diag_Ch;
	diag_Ch.Reallocate(mesh_num.GetNbDof());
	diag_Ch.Zero();
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    HarmonicMaxwell3D_PhysGeomInfo<Real_wp>& mass = var_maxwell.Glob_matMass_elem(i);
	    int offsetS = 0, offsetL = 0;
	    for (int num_loc = 0; num_loc < 6; num_loc++)
	      {
		int num_face = mesh.Element(i).numBoundary(num_loc);
		int ref_boundary = mesh.Boundary(num_face).GetReference();
		Real_wp coef(1);
		if ((ref_boundary != 0) && (mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_ABSORBING))
		  coef = Real_wp(2);
		
		const ElementReference<Dimension2, 2>& Fb_s = var_problem.GetSurfaceFiniteElementHcurl(num_face);
		const QuadrangleHcurlLobatto& Fb_quad = dynamic_cast<const QuadrangleHcurlLobatto&>(Fb_s);
		int nb_quad = Fb_s.GetNbPointsQuadratureInside();
		const IVect& num_dof = mesh_num.Element(i).GetNodle();
		for (int k = 0; k < nb_quad; k++)
		  {
		    const Matrix2_2sym& Bh = mass.GetBhSurf(offsetS+k);
		    int kx = Fb_quad.GetXdofNumber(k);
		    int ky = Fb_quad.GetYdofNumber(k);
		    diag_Ch(num_dof(offsetL+kx)) -= coef*Bh(0, 0);
		    diag_Ch(num_dof(offsetL+ky)) -= coef*Bh(1, 1);
		  }
		
		offsetL += Fb_s.GetNbDof();
		offsetS += nb_quad;
	      }	    
	  }

	// the diagonal is assembled
	var_problem.AddDomains(diag_Ch, 1);
	
	// we fill SurfMass
	DiscontinuousBlockMassMatrix* block;
	VectReal_wp diagonal;
	this->SurfMass->ReallocateStoredMass();
	int offsetL = 0;
	for (int i = 0; i < mesh.GetNbBoundary(); i++)
	  {
	    const ElementReference<Dimension2, 2>& Fb_s = var_problem.GetSurfaceFiniteElementHcurl(i);

	    diagonal.Reallocate(Fb_s.GetNbDof());
	    for (int k = 0; k < Fb_s.GetNbDof(); k++)
	      diagonal(k) = diag_Ch(offsetL + k);
	    
	    block = new DiscontinuousDiagonalMassMatrix(diagonal);
	    this->SurfMass->SetLocalStoredBlock(i, block);
	    offsetL += Fb_s.GetNbDof();
	  }
	
	// then ChSurf
	VectReal_wp coef_block(mesh.GetNbBoundary());
	coef_block.Fill(Real_wp(1));
	this->ChSurf->CopyBlockCoef(coef_block, 0);
	
	//diag_Ch.Write("Ch.dat");
      }
    else if (type_matrix == FemMassMatrix::BLOCK_DIAGONAL)
      {
	// Ch stored as a block-diagonal matrix with blocks 2x2
	Vector<int> dof_block_Ch;
	Vector<Matrix2_2sym> block_Ch;
	
	int nb_blocks = mesh_num.GetNbDof()/2;
	block_Ch.Reallocate(nb_blocks);
	dof_block_Ch.Reallocate(2*nb_blocks);
	Vector<int> IndexBlock(2*nb_blocks);
	IndexBlock.Fill(-1);

	nb_blocks = 0;
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    int offsetL = 0;
	    for (int num_loc = 0; num_loc < 6; num_loc++)
	      {
		int num_face = mesh.Element(i).numBoundary(num_loc);
		const ElementReference<Dimension2, 2>& Fb_s = var_problem.GetSurfaceFiniteElementHcurl(num_face);
		const QuadrangleHcurlLobatto& Fb_quad = dynamic_cast<const QuadrangleHcurlLobatto&>(Fb_s);
	      	int nb_quad = Fb_s.GetNbPointsQuadratureInside();
		const IVect& num_dof = mesh_num.Element(i).GetNodle();
		for (int k = 0; k < nb_quad; k++)
		  {
		    int kx = Fb_quad.GetXdofNumber(k);
		    int ky = Fb_quad.GetYdofNumber(k);
		    if (IndexBlock(num_dof(offsetL+kx)) == -1)
		      {
			dof_block_Ch(2*nb_blocks) = num_dof(offsetL + kx);
			dof_block_Ch(2*nb_blocks+1) = num_dof(offsetL + ky);
			IndexBlock(num_dof(offsetL+kx)) = nb_blocks;
			IndexBlock(num_dof(offsetL+ky)) = nb_blocks;
			nb_blocks++;
		      }
		  }
		
		offsetL += Fb_s.GetNbDof();
	      }
	  }
	
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    HarmonicMaxwell3D_PhysGeomInfo<Real_wp>& mass = var_maxwell.Glob_matMass_elem(i);
	    int offsetS = 0, offsetL = 0;
	    const IVect& negative_dof = mesh_num.Element(i).GetNegativeDofNumber();
	    for (int num_loc = 0; num_loc < 6; num_loc++)
	      {
		int num_face = mesh.Element(i).numBoundary(num_loc);
		int ref_boundary = mesh.Boundary(num_face).GetReference();
		Real_wp coef(1);
		if ((ref_boundary != 0) && (mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_ABSORBING))
		  coef = Real_wp(2);
		
		const ElementReference<Dimension2, 2>& Fb_s = var_problem.GetSurfaceFiniteElementHcurl(num_face);
		const QuadrangleHcurlLobatto& Fb_quad = dynamic_cast<const QuadrangleHcurlLobatto&>(Fb_s);
		int nb_quad = Fb_s.GetNbPointsQuadratureInside();
		const IVect& num_dof = mesh_num.Element(i).GetNodle();
		Vector<bool> SignDof(Fb_s.GetNbDof());
		SignDof.Fill(true);
		int order = Fb_quad.GetOrder();
		for (int k = 2*(order+1); k < 4*(order+1); k++)
		  SignDof(k) = false;

		for (int k = 0; k < negative_dof.GetM(); k++)
		  if ((negative_dof(k) >= offsetL) && (negative_dof(k) < offsetL+Fb_s.GetNbDof()))
		    SignDof(negative_dof(k) - offsetL) = !SignDof(negative_dof(k) - offsetL);
		
		for (int k = 0; k < nb_quad; k++)
		  {
		    const Matrix2_2sym& Bh = mass.GetBhSurf(offsetS+k);
		    int kx = Fb_quad.GetXdofNumber(k);
		    int ky = Fb_quad.GetYdofNumber(k);
		    int nblock = IndexBlock(num_dof(offsetL+kx));
		    int posx = 0, posy = 1;
		    if (num_dof(offsetL+kx) == dof_block_Ch(2*nblock+1))
		      { posx = 1; posy = 0; }

		    block_Ch(nblock)(posx, posx) += coef*Bh(0, 0);
		    block_Ch(nblock)(posy, posy) += coef*Bh(1, 1);
		    if (SignDof(kx) != SignDof(ky))
		      block_Ch(nblock)(posx, posy) -= coef*Bh(0, 1);
		    else
		      block_Ch(nblock)(posx, posy) += coef*Bh(0, 1);
		  }
		
		offsetL += Fb_s.GetNbDof();
		offsetS += nb_quad;
	      }	    
	  }

#ifdef SELDON_WITH_MPI
	// diagonal is assembled
	VectReal_wp diag(mesh_num.GetNbDof());
	diag.Zero();
	for (int i = 0; i < nb_blocks; i++)
	  {
	    diag(dof_block_Ch(2*i)) = block_Ch(i)(0, 0);
	    diag(dof_block_Ch(2*i+1)) = block_Ch(i)(1, 1);
	  }

	var_problem.AddDomains(diag, 1);
	for (int i = 0; i < nb_blocks; i++)
	  {
	    block_Ch(i)(0, 0) = diag(dof_block_Ch(2*i));
	    block_Ch(i)(1, 1) = diag(dof_block_Ch(2*i+1));

	    diag(dof_block_Ch(2*i)) = block_Ch(i)(0, 1);
	    diag(dof_block_Ch(2*i+1)) = block_Ch(i)(0, 1);
	  }

	// and extra-diagonal
	var_problem.AddDomains(diag, 1);
	for (int i = 0; i < nb_blocks; i++)
	  block_Ch(i)(0, 1) = diag(dof_block_Ch(2*i));	
#endif

	// we fill SurfMass
	Vector<bool> BoundaryUsed(mesh.GetNbBoundary());
	BoundaryUsed.Fill(false);
	DiscontinuousBlockMassMatrix* block;
	Vector<int> num_row;
	Vector<TinyMatrix<Real_wp, Symmetric, 2, 2> > Ch_b;
	this->SurfMass->ReallocateStoredMass();
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    const IVect& num_dof = mesh_num.Element(i).GetNodle();
	    int offsetL = 0;
	    for (int num_loc = 0; num_loc < 6; num_loc++)
	      {
		int num_face = mesh.Element(i).numBoundary(num_loc);
		const ElementReference<Dimension2, 2>& Fb_s = var_problem.GetSurfaceFiniteElementHcurl(num_face);
		if (!BoundaryUsed(num_face))
		  {
		    BoundaryUsed(num_face) = true;		  
		    const QuadrangleHcurlLobatto& Fb_quad = dynamic_cast<const QuadrangleHcurlLobatto&>(Fb_s);
		    int nb_points_quad = Fb_s.GetNbPointsQuadratureInside();
		    num_row.Reallocate(2*nb_points_quad);
		    Ch_b.Reallocate(nb_points_quad);
		    int n0 = num_dof(offsetL);
		    for (int k = 1; k < Fb_s.GetNbDof(); k++)
		      if (num_dof(offsetL+k) < n0)
			n0 = num_dof(offsetL+k);
		    
		    for (int k = 0; k < nb_points_quad; k++)
		      {
			int kx = Fb_quad.GetXdofNumber(k);
			int ky = Fb_quad.GetYdofNumber(k);
			int nblock = IndexBlock(num_dof(offsetL+kx));
			
			num_row(2*k) = num_dof(offsetL+kx) - n0;
			num_row(2*k+1) = num_dof(offsetL+ky) - n0;
			Ch_b(k) = block_Ch(nblock);
		      }
		    
		    block = new DiscontinuousBlockDiagMassMatrix<2>(num_row, Ch_b);
		    this->SurfMass->SetLocalStoredBlock(num_face, block);
		  }
		
		offsetL += Fb_s.GetNbDof();
	      }
	  }

	// then ChSurf
	VectReal_wp coef_block(mesh.GetNbBoundary());
	coef_block.Fill(Real_wp(-1));
	this->ChSurf->CopyBlockCoef(coef_block, 0);

	/* ofstream file_out("Ch.dat");
	file_out.precision(15);
	for (int i = 0; i < nb_blocks; i++)
	  {
	    int nx = dof_block_Ch(2*i);
	    int ny = dof_block_Ch(2*i+1);
	    file_out << nx+1 << " " << nx+1 << " " << block_Ch(i)(0, 0) << '\n';
	    file_out << nx+1 << " " << ny+1 << " " << block_Ch(i)(0, 1) << '\n';
	    file_out << ny+1 << " " << nx+1 << " " << block_Ch(i)(1, 0) << '\n';
	    file_out << ny+1 << " " << ny+1 << " " << block_Ch(i)(1, 1) << '\n';
	  }

	file_out.close();
	*/
	
      }
    else
      {
	cout << "Storage not implemented " << endl;
	abort();
      }

    // testing SolveOperatorCh
    /* VectReal_wp X(mesh_num.GetNbDof()), Y(mesh_num.GetNbDof());
    X.FillRand(); X *= 1e-9; X.Write("X.dat");

    this->SolveOperatorCh(X);

    X.Write("Y.dat");
    */
  }

  
  //! solves by operator Ch (Y is overwritten with C_h^{-1} Y)
  void DiscontinuousMassMatrixHdgMaxwell3D::SolveOperatorCh(VectReal_wp& Y)
  {
    DiscontinuousUnsteadyMassMatrix<Dimension3>::SolveOperatorCh(Y);
    
    /*var_problem.AddDomains(Y);
    if (type_matrix == FemMassMatrix::DIAGONAL)
      Y = Y*diag_Ch;
    else if (type_matrix == FemMassMatrix::BLOCK_DIAGONAL)
      {
	int nb_blocks = block_Ch.GetM(); int nb = 0;
	R2 vec_u, vec_v;
	for (int i = 0; i < nb_blocks; i++)
	  {
	    vec_u.Init(Y(dof_block_Ch(nb)), Y(dof_block_Ch(nb+1)));
	    Mlt(block_Ch(i), vec_u, vec_v);
	    Y(dof_block_Ch(nb)) = vec_v(0);
	    Y(dof_block_Ch(nb+1)) = vec_v(1);
	    nb += 2;
	  }
      }
    */
  }
     
} // end namespace

#define MONTJOIE_FILE_TIME_MAXWELL3D_CXX
#endif
  
