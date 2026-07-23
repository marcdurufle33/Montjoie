#ifndef MONTJOIE_FILE_VISCO_THERMAL_EQUATION_HXX

#include "ViscoThermalEquation.hxx"

namespace Montjoie
{

  /************************
   * ViscoThermalEquation *
   ************************/

  // providing C in Cgrad_grad
  template<class Dimension, class T>
  template<class TypeEquation, class NatureMat, class MatStiff>
  void ViscoThermalEquation<Dimension, T>
  ::GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                      int num_elem, int jloc, const NatureMat& nat_mat,
                      int ref, MatStiff& Cgrad)
  {
    Cgrad.Zero();
    //DISP(vars.c_p);     DISP(vars.T0);
    
    if (Dimension::dim_N == 2)
      {
        Cgrad(0,0)(0,0) = 2*vars.mu+vars.zeta;
        Cgrad(1,1)(1,1) = 2*vars.mu+vars.zeta;
        Cgrad(1,1)(0,0) = vars.mu;
        Cgrad(0,0)(1,1) = vars.mu;
        Cgrad(1,0)(0,1) = vars.mu;
        Cgrad(0,1)(1,0) = vars.mu;
        Cgrad(1,0)(1,0) = vars.zeta;
        Cgrad(0,1)(0,1) = vars.zeta;
        Cgrad(2,2)(0,0) = -vars.lambda/vars.T0;
        Cgrad(2,2)(1,1) = -vars.lambda/vars.T0;
      }
    else
      {
        Cgrad(0,0)(0,0) = 2*vars.mu+vars.zeta;
        Cgrad(1,1)(1,1) = 2*vars.mu+vars.zeta;
        Cgrad(2,2)(2,2) = 2*vars.mu+vars.zeta;
        Cgrad(0,0)(1,1) = vars.mu;
        Cgrad(1,0)(0,1) = vars.mu;
        Cgrad(0,1)(1,0) = vars.mu;
        Cgrad(1,1)(0,0) = vars.mu;
        Cgrad(0,0)(2,2) = vars.mu;
        Cgrad(2,0)(0,2) = vars.mu;
        Cgrad(0,2)(2,0) = vars.mu;
        Cgrad(2,2)(0,0) = vars.mu;
        Cgrad(1,1)(2,2) = vars.mu;
        Cgrad(2,1)(1,2) = vars.mu;
        Cgrad(1,2)(2,1) = vars.mu;
	Cgrad(2,2)(1,1) = vars.mu;
        Cgrad(1,0)(1,0) = vars.zeta;
        Cgrad(2,0)(2,0) = vars.zeta;
        Cgrad(0,1)(0,1) = vars.zeta;
        Cgrad(2,1)(2,1) = vars.zeta;
        Cgrad(0,2)(0,2) = vars.zeta;
        Cgrad(1,2)(1,2) = vars.zeta;
        Cgrad(3,3)(0,0) = -vars.lambda/vars.T0;
        Cgrad(3,3)(1,1) = -vars.lambda/vars.T0;
        Cgrad(3,3)(2,2) = -vars.lambda/vars.T0;
      }
   }



  template<class Dimension, class T>
  template<class TypeEquation, class T0, class MatStiff>
  void ViscoThermalEquation<Dimension, T>::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                                            int num_elem, int jloc,
                                                            const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                                            MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
   {
     if (Dimension::dim_N == 2)
       {
         if (vars.alternative_ipp == false)
           {
             Dgrad_phi.Zero(); Ephi_grad.Zero();
             Dgrad_phi(0,3)(0) = -1.;
             Dgrad_phi(1,3)(1) = -1.;
             Ephi_grad(3,0)(0) = -1.;
             Ephi_grad(3,1)(1) = -1.;
           }
         else
           {
             Dgrad_phi.Zero(); Ephi_grad.Zero();
             Dgrad_phi(3,0)(0) = 1.;
             Dgrad_phi(3,1)(1) = 1.;
             Ephi_grad(0,3)(0) = 1.;
             Ephi_grad(1,3)(1) = 1.;
           }
       }
     else
       {
         Dgrad_phi.Zero(); Ephi_grad.Zero();
         if (vars.alternative_ipp ==false)
           {
             Dgrad_phi(0,4)(0) = -1.;
             Dgrad_phi(1,4)(1) = -1.;
             Dgrad_phi(2,4)(2) = -1.;
             Ephi_grad(4,0)(0) = -1.;
             Ephi_grad(4,1)(1) = -1.;
             Ephi_grad(4,2)(2) = -1.; 
           }
         else
           {
             Dgrad_phi(4,0)(0) = 1.;
             Dgrad_phi(4,1)(1) = 1.;
             Dgrad_phi(4,2)(2) = 1.;
             Ephi_grad(0,4)(0) = 1.;
             Ephi_grad(1,4)(1) = 1.;
             Ephi_grad(2,4)(2) = 1.;
           }
       }
  }

  template<class Dimension, class T>
  template<class TypeEquation, class T0, class MatMass>
  void ViscoThermalEquation<Dimension, T>::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                                                         int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    mass.Zero(); 
    Complex_wp m_iomega; vars.GetMiomega(m_iomega);

    if (Dimension::dim_N == 2)
      {
        mass(0,0) = -m_iomega*vars.rho0;
        mass(1,1) = -m_iomega*vars.rho0;
        mass(2,2) = (m_iomega*vars.rho0*vars.c_p)/vars.T0;
        mass(2,3) = (-m_iomega)/vars.T0;
        mass(3,2) = (-m_iomega)/vars.T0;
        mass(3,3) = m_iomega/vars.P0;
      }
    else
      {
        mass(0,0) = -m_iomega*vars.rho0;
        mass(1,1) = -m_iomega*vars.rho0;
        mass(2,2) = -m_iomega*vars.rho0;
        mass(3,3) = (m_iomega*vars.rho0*vars.c_p)/vars.T0;
        mass(4,3) = (-m_iomega)/vars.T0;
        mass(3,4) = (-m_iomega)/vars.T0;
        mass(4,4) = m_iomega/vars.P0;
      }

    //DISP(mass);
    
    /*Complex_wp m_iomega; vars.GetMiomega(m_iomega);
    Complex_wp m_omega2; vars.GetMomega2(m_omega2);
    mass(0, 0) = m_omega2*nat_mat.GetCoefMass();
    mass(1, 1) = m_omega2*nat_mat.GetCoefMass();
    mass(2, 2) = m_iomega*nat_mat.GetCoefDamping(); */
  }

  /*****************                                                                                                                                                                                              
  ViscoThermalIndex                                                                                                                                                                                               
  ******************/

  void ViscoThermalIndex::UpdateTemp(Real_wp Temp, Real_wp h_v, Real_wp x_c)
  {
    T0 = 293.15;
    zeta = 1.3e-5;
    Real_wp puis;
    puis = 5.21899-5.8294*T0/Temp-1.0252*pow(T0/Temp,2);
    if (puis>0) // possible probleme de puissance negative
      x_v = h_v*pow(10, puis);
    else
      x_v = h_v/pow(10, -puis);

    c_p = 1012.25*(1 + 0.5438*(x_v-x_v0) + 0.638*pow(x_v-x_v0,2) - 0.1594*(x_c-x_c0) + 0.075*pow(x_c-x_c0,2) + 0.00952*(Temp/T0-1) + 0.0406*pow(Temp/T0-1,2) + 0.3976*(x_c-x_c0)*(Temp/T0-1));
    mu = 1.8206e-5*(1 + 0.77013*(Temp/T0-1));
    lambda = 2.5562e-2*(1 + 0.8490*(Temp/T0-1));
    rho0 = 1.19930*(T0/Temp)*(1 - 0.3767*(x_v-x_v0) + 0.4162*(x_c-x_c0) - 0.0029*(Temp/T0-1));
    //c_0 = 343.986*sqrt(Temp/T0)*sqrt(1 + 0.314*(x_v-x_v0) - 0.520*(x_c-x_c0) + 0.25*pow(x_c-x_c0,2) - 0.16*(x_c-x_c0)*(Temp/T0-1));
    //gamma = 1.40108*(1 - 0.060*(x_v-x_v0) - 0.104*(x_c-x_c0) - 0.0087*(Temp/T0-1) - 0.154*(x_c-x_c0)*(Temp/T0-1));
    //P0 = (pow(c_0,2)*rho0)/gamma;
    P0 = 101325;
  }

  void ViscoThermalIndex::SetInputData(const string& description_field, const VectString& parameters)
  {
    bool updatetemp = false;
    if (description_field == "Temperature")
      {
        updatetemp = true;
        Temp = to_num<Real_wp>(parameters(0));
      }
    else if (description_field == "Humidity")
      {
        updatetemp = true;
        h_v = to_num<Real_wp>(parameters(0));
      }
    else if (description_field == "CO2mixing")
      {
        updatetemp = true;
        x_c = to_num<Real_wp>(parameters(0));
      }
    else
      {
        updatetemp = false;
      }
    if (updatetemp)
      {
	ViscoThermalIndex::UpdateTemp(Temp,h_v,x_c);
	DISP(c_p); DISP(T0); DISP(P0); DISP(rho0); DISP(mu); DISP(lambda); DISP(zeta);
      }
  }



  /**********************
   * VarViscoThermal_Eq *
   **********************/


  template<class TypeEquation>
  void VarViscoThermal_Eq<TypeEquation>::SetInputData(const string& description_field, const VectString& parameters)
  {
    ViscoThermalIndex::SetInputData(description_field, parameters);
    if (description_field == "AlternativeIPP")
      {
        if (parameters(0) == "YES")
          alternative_ipp = true;
        else
          alternative_ipp = false;
      }
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
  }

  template<class TypeEquation>
  void VarViscoThermal_Eq<TypeEquation>::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
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
  }

  template<class TypeEquation>
  string VarViscoThermal_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return string("");
  }

  

  /**********************
   * EllipticProblem *
   **********************/


  template<class Dimension>
  void EllipticProblem<ViscoThermalEquation<Dimension, Complex_wp> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>&,
                            const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl; abort();
  }

  template<class Dimension>
  void EllipticProblem<ViscoThermalEquation<Dimension, Complex_wp> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
                          const GlobalGenericMatrix<Real_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
    cout << "Not possible" << endl; abort();
  }


  template<class Dimension>
  void EllipticProblem<ViscoThermalEquation<Dimension, Complex_wp> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
                            const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
                                      this->GetReferenceElementH1(i));
  }

  template<class Dimension>
  void EllipticProblem<ViscoThermalEquation<Dimension, Complex_wp> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
                          const GlobalGenericMatrix<Complex_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }


}

#define MONTJOIE_FILE_VISCO_THERMAL_EQUATION_HXX
#endif
