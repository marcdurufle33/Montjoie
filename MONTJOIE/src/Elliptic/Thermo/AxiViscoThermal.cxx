#ifndef MONTJOIE_FILE_AXI_VISCO_THERMAL_CXX

#include "AxiViscoThermal.hxx"

namespace Montjoie
{
  
  template<class T>
  template<class TypeEquation, class NatureMat, class MatStiff>
  void ViscoThermalEquationAxi<T>::GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                                     int num_elem, int jloc, const NatureMat& nat_mat,
                                                     int ref, MatStiff& Cgrad_grad)
  {
    Cgrad_grad.Zero();
    /*Real_wp r = vars.Glob_radius(num_elem)(jloc);
    Cgrad_grad(0,0)(0,0) = 2.0*M_PI*(2.0*vars.mu + vars.zeta)*r;
    Cgrad_grad(1,1)(1,1) = 2.0*M_PI*(2.0*vars.mu + vars.zeta)*r;
    Cgrad_grad(1,1)(0,0) = 2.0*M_PI*vars.mu*r;
    Cgrad_grad(0,0)(1,1) = 2.0*M_PI*vars.mu*r;
    Cgrad_grad(1,0)(0,1) = 2.0*M_PI*vars.mu*r;
    Cgrad_grad(0,1)(1,0) = 2.0*M_PI*vars.mu*r;
    Cgrad_grad(2,2)(0,0) = 2.0*M_PI*(-vars.lambda*r)/vars.T0;
    Cgrad_grad(2,2)(1,1) = 2.0*M_PI*(-vars.lambda*r)/vars.T0;
    Cgrad_grad(1,0)(1,0) = 2.0*M_PI*vars.zeta*r;
    Cgrad_grad(0,1)(0,1) = 2.0*M_PI*vars.zeta*r;
    Cgrad_grad(2,2)(2,2) = 2.0*M_PI*(-vars.lambda)/(vars.T0*r);*/
    /*Cgrad_grad(2, 2)(0, 0) = 1.0;
      Cgrad_grad(2, 2)(1, 1) = 1.0;*/
    //DISP(Cgrad_grad);
  }

  
  template<class T>
  template<class TypeEquation, class T0, class MatStiff>
  void ViscoThermalEquationAxi<T>
  ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                     int num_elem, int jloc,
                     const GlobalGenericMatrix<T0>& nat_mat, int ref,
                     MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    Dgrad_phi.Zero();
    Ephi_grad.Zero();
    //DISP(Dgrad_phi); DISP(Ephi_grad);
    Real_wp r = vars.Glob_radius(num_elem)(jloc);
    if (vars.alternative_ipp == false)
      {
	Dgrad_phi(0,3)(0) = -2.0*M_PI*r;
	Dgrad_phi(1,3)(1) = -2.0*M_PI*r;
	//Dgrad_phi(0,0)(0) = vars.zeta*2.0*M_PI;
	//Dgrad_phi(1,0)(1) = vars.zeta*2.0*M_PI;
	Ephi_grad(3,1)(1) = -2.0*M_PI*r;
	Ephi_grad(3,0)(0) = -2.0*M_PI*r;
	//Ephi_grad(0,0)(0) = vars.zeta*2.0*M_PI;
	//Ephi_grad(0,1)(1) = vars.zeta*2.0*M_PI;
	}
    else
      {
	Dgrad_phi(3,0)(0) = 2.0*M_PI*r;
	Dgrad_phi(3,1)(1) = 2.0*M_PI*r;
	//Dgrad_phi(0,0)(0) = vars.zeta*2*M_PI;
        //Dgrad_phi(1,0)(1) = vars.zeta*2*M_PI;
	Ephi_grad(1,3)(1) = 2.0*M_PI*r;
	Ephi_grad(0,3)(0) = 2.0*M_PI*r;
	//Ephi_grad(0,0)(0) = vars.zeta*2.0*M_PI;
	//Ephi_grad(0,1)(1) = vars.zeta*2.0*M_PI;
	}
    /*Dgrad_phi(3, 0)(0) = -r;
    Dgrad_phi(3, 1)(1) = -r;
    Ephi_grad(0, 3)(0) = -r;
    Ephi_grad(1, 3)(1) = -r;*/
    //DISP(Dgrad_phi); DISP(Ephi_grad);
  }
  

  template<class T>
  template<class TypeEquation, class T0, class MatMass>
  void ViscoThermalEquationAxi<T>
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                  int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    mass.Zero();
    Real_wp r = vars.Glob_radius(i)(j);
    T0 m_iomega; vars.GetMiomega(m_iomega);
    mass(0, 0) = -2.0*M_PI*m_iomega*vars.rho0*r; //+ (4.0*M_PI*vars.mu)/r + (2.0*M_PI*vars.zeta)/r;
    mass(1, 1) = -2.0*M_PI*m_iomega*vars.rho0*r;
    mass(2, 2) = (2.0*M_PI*m_iomega*vars.rho0*vars.c_p*r)/vars.T0;
    mass(2, 3) = (2.0*M_PI*(-m_iomega)*r)/vars.T0;
    mass(3, 2) = (2.0*M_PI*(-m_iomega)*r)/vars.T0;
    mass(3, 3) = (2.0*M_PI*m_iomega*r)/vars.P0;
    //DISP(vars.T0);
    //mass(3,0) = -2.0*M_PI;
    //mass(0,3) = -2.0*M_PI;
    //DISP(mass);
    if (!vars.alternative_ipp)
      {
	mass(3,0) = -2.0*M_PI;
	mass(0,3) = -2.0*M_PI;
      }
    /*mass(0, 0) = m_iomega*r;
    mass(1, 1) = m_iomega*r;
    mass(2, 2) = m_iomega*r;
    mass(3, 3) = -m_iomega*r;*/
  }


  template<class TypeEquation>
  VarViscoThermal_Axi<TypeEquation>::VarViscoThermal_Axi(EllipticProblem<TypeEquation>& var)
    : VarAxisymProblem(var), var_problem(var)
  {
  }
  
  
  template<class TypeEquation>
  void VarViscoThermal_Axi<TypeEquation>::AllocateMassMatrices()
  {
    Glob_radius.Reallocate(var_problem.mesh.GetNbElt());    
  }


  template<class TypeEquation>
  void VarViscoThermal_Axi<TypeEquation>
  ::ComputeLocalMassMatrix(int i, int N, bool linear_sparse,
                           SetPoints<Dimension2>& PointsElem,
                           SetMatrices<Dimension2>& MatricesElem,
                           IVect& OrderFace, const ElementGeomReference<Dimension2>& Fb)
  {
    VarHarmonic<TypeEquation>::ComputeLocalMassMatrix(i, N, linear_sparse, PointsElem,
                                                      MatricesElem, OrderFace, Fb);
    
    this->Glob_radius(i).Reallocate(N);
    for (int j = 0; j < N; j++)
      this->Glob_radius(i)(j) = var_problem.Glob_PointsQuadrature(i)(j)(0);    
  }

  
  template<class TypeEquation>
  void VarViscoThermal_Axi<TypeEquation>::UpdateDirichlet(int n)
  {
  }
    

  template<class TypeEquation>
  void VarViscoThermal_Axi<TypeEquation>::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (description_field == "AlternativeIPP")
      {
	if (parameters(0) == "YES")
	  {
	    alternative_ipp = true;
	  }
	else
	  {
	    alternative_ipp = false;
	  }
      }
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
  }

  
  template<class TypeEquation>
  void VarViscoThermal_Axi<TypeEquation>::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    if (name_media == "cp")                                                                                                                                                                                    
      c_p = to_num<Real_wp>(parameters(1));                                                                                                                                                                    
    else if (name_media == "T0")                                                                                                                                                                               
      T0 = to_num<Real_wp>(parameters(1));                                                                                                                                                                     
    else if (name_media == "P0")                                                                                                                                                                               
      P0 = to_num<Real_wp>(parameters(1));                                                                                                                                                                     
    else if (name_media == "rho0")                                                                                                                                                                             
      rho0 = to_num<Real_wp>(parameters(1));                                                                                                                                                                   
    else if (name_media ==  "mu")                                                                                                                                                                              
      mu = to_num<Real_wp>(parameters(1));                                                                                                                                                                     
    else if (name_media == "lambda")                                                                                                                                                                           
      lambda = to_num<Real_wp>(parameters(1));                                                                                                                                                                 
    else if (name_media == "zeta")                                                                                                                                                                             
      zeta = to_num<Real_wp>(parameters(1));
    //DISP(c_p); DISP(T0); DISP(P0); DISP(rho0); DISP(mu); DISP(lambda); DISP(zeta);
  }
  

  template<class TypeEquation>
  string VarViscoThermal_Axi<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return string("");
  }

  template<class TypeEquation>
  VolumetricSource<ViscoThermalEquationAxi<Complex_wp> >
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
                     const Vector<VectString>& param)
    : VirtualSourceFEM<Complex_wp, Dimension2>(var), var_problem(var)
  {
    fsrc = NULL;

    GaussianSourceField<Complex_wp, Dimension3>* f_gauss = NULL;

    Vector<Complex_wp> polar;
    var.GetPolarization(polar);
    unif_source_surf.Reallocate(var.mesh.GetNbReferences()+1);
    gaussian_source_surf.Reallocate(var.mesh.GetNbReferences()+1);
    coef_source_surf.Reallocate(var_problem.mesh.GetNbReferences()+1);
    coef_source_surf.Fill(polar);
    unif_source_surf.Fill(false); gaussian_source_surf.Fill(false);
    
    for (int num = 0; num < param.GetM(); num++)
      {
	if (param(num).GetM() > 1)
	  if (param(num)(0) == "SRC_SURFACE")
            {
              for (int k = 1; k < param(num).GetM(); k++)
                {
                  if (param(num)(k) == "Polarization")
                    {
                      for (int p = 0; p < polar.GetM(); p++)
                        polar(p) = to_num<Complex_wp>(param(num)(k+p+1));
                      
                      break;
                    }
                }
              
              if (param(num)(2) == "UNIFORM")
                {
                  int ref = to_num<int>(param(num)(1));
                  unif_source_surf(ref) = true;
                  coef_source_surf(ref) = polar;
                }
              
              if (param(num)(2) == "GAUSSIAN")
                {
                  int ref = to_num<int>(param(num)(1));
                  gaussian_source_surf(ref) = true;
                  
                  int nb = 3;
                  f_gauss = new GaussianSourceField<Complex_wp, Dimension3>();                
                  var.InitGaussianParameter(*f_gauss, param(num), nb);
                  f_gauss->SetPolarization(polar);
                  fsrc = f_gauss;
                }
            }
      }
  }
  

  VolumetricSource<ViscoThermalEquationAxi<Complex_wp> >::~VolumetricSource()
  {
    if (fsrc != NULL)
      delete fsrc;
  }


  void VolumetricSource<ViscoThermalEquationAxi<Complex_wp> >
  ::SetSurfaceSource(int ref, VirtualSourceField<Complex_wp, Dimension3>* f)
  {
    if (fsrc != NULL)
      delete fsrc;
    
    fsrc = f;
  }

  
  bool VolumetricSource<ViscoThermalEquationAxi<Complex_wp> >
  ::IsNonNullSurfacicSource(int ref)
  {
    if (unif_source_surf(ref))
      return true;

    if (gaussian_source_surf(ref))
      return true;

    return false;
  }
  

  void VolumetricSource<ViscoThermalEquationAxi<Complex_wp> >
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
                           const SetMatrices<Dimension2>& MatricesElem, Vector<Complex_wp>& f)
  {
    R2 x = PointsElem.GetPointQuadratureBoundary(k);
    int ref = this->ref_boundary_;
    if (unif_source_surf(ref))
      {
        for (int k = 0; k < coef_source_surf(ref).GetM(); k++)
          f(k) = x(0)*coef_source_surf(ref)(k);
      }

    if (gaussian_source_surf(ref))
      {
        EvaluateMode(x, f);
        f *= x(0);
      }
  }
  
  void VolumetricSource<ViscoThermalEquationAxi<Complex_wp> >
  ::EvaluateFunction(int i, int j, const R2& x, Vector<Complex_wp>& f)
  {
    f.Zero();
    int ref = this->ref_boundary_;
    if (unif_source_surf(ref))
      {
        for (int k = 0; k < coef_source_surf(ref).GetM(); k++)
          f(k) = coef_source_surf(ref)(k);
      }

    if (gaussian_source_surf(ref))
      {
        EvaluateMode(x, f);
      }
  }

  void VolumetricSource<ViscoThermalEquationAxi<Complex_wp> >
  ::EvaluateMode(const R2& x, Vector<Complex_wp>& f)
  {
    R3 pt3D(x(0), 0, x(1));
    fsrc->EvaluateFunction(pt3D, f);
  }

  void EllipticProblem<ViscoThermalEquationAxi<Complex_wp> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>&,
                            const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Case impossible" << endl;
    abort();
  }
  
  
  void EllipticProblem<ViscoThermalEquationAxi<Complex_wp> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
                          const GlobalGenericMatrix<Real_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
  }

  
  void EllipticProblem<ViscoThermalEquationAxi<Complex_wp> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
                            const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
                                      this->GetReferenceElementH1(i));
  }

  
  void EllipticProblem<ViscoThermalEquationAxi<Complex_wp> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
                          const GlobalGenericMatrix<Complex_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
  }
  
}

#define MONTJOIE_FILE_AXI_VISCO_THERMAL_CXX
#endif
