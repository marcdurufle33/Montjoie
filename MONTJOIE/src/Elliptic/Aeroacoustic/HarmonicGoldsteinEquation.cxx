#ifndef MONTJOIE_FILE_HARMONIC_GOLDSTEIN_EQUATION_CXX
namespace Montjoie
{
  

  /*********************
   * GoldsteinEquation *
   *********************/
  
  template<class T, class Dim>
  template<class TypeEquation, class T0, class Vector1>
  void GoldsteinEquation<T,Dim>::GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
                                                     const GlobalGenericMatrix<T0>& nat_mat,
                                                     Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    unknown_to_derive.Fill(true);
    fct_test_to_derive.Fill(true); 
  }
  
  
  template<class T, class Dim>
  template<class TypeEquation, class NatureMat, class MatStiff>
  void GoldsteinEquation<T,Dim>::GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                                   int num_elem, int jloc, const NatureMat& nat_mat,
                                                   int ref, MatStiff& Cgrad_grad)
  {
    FillZero(Cgrad_grad); 
    //Real_wp rho = vars.eval_rho0(num_elem)(jloc);
    // Convected helmholtz
    //Cgrad_grad(0,0).SetIdentity();
    //Cgrad_grad(0,0) *= rho* nat_mat.GetCoefStiffness();
    // Transport part
    Real_wp d = vars.penal_supg;
    R_N v0 = vars.eval_flow(num_elem)(jloc);
    for(int m=0;m<nb_unknowns;m++)
      {
        for(int i=0;i<Dimension::dim_N;i++)
          {
            for(int j=0;j<Dimension::dim_N;j++)
              {
                Cgrad_grad(m,m)(i,j) = v0(i)*v0(j);
                Cgrad_grad(m,m)(i,j) *= d*nat_mat.GetCoefStiffness(); 
              }
          }
      }
  }
  
  
  template<class T, class Dim>
  template<class TypeEquation, class T0, class MatMass>
  void GoldsteinEquation<T,Dim>::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                                               int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    FillZero(mass);
    MatrixN_N grad_M = vars.grad_flow(num_elem)(jloc);
    R_N M = vars.eval_flow(num_elem)(jloc);
    T omeg;
    vars.GetMiomega(omeg);
    for(int m=0;m<nb_unknowns;m++)
      {
        mass(m,m) = omeg;
        for(int n=0;n<nb_unknowns;n++)
          mass(m,n) += grad_M(n,m);
      }
    mass *= nat_mat.GetCoefStiffness();
    
  }
  
  
  template<class T, class Dim>
  template<class TypeEquation, class T0, class MatStiff>
  void GoldsteinEquation<T,Dim>
  ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                     int num_elem, int jloc,
                     const GlobalGenericMatrix<T0>& nat_mat, int ref,
                     MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
    
    MatrixN_N grad_M = vars.grad_flow(num_elem)(jloc);
    R_N M = vars.eval_flow(num_elem)(jloc);
    Real_wp d = vars.penal_supg;
    T omeg;
    vars.GetMiomega(omeg);
    
    for(int m=0;m<nb_unknowns;m++)
      {
        //Dgrad_phi(m,m) = omeg*M(m);
        //Ephi_grad(m,m) = omeg;
        
        for(int i=0;i<Dimension::dim_N;i++)
          {
            Dgrad_phi(m,m)(i) = omeg*M(i);
            Ephi_grad(m,m)(i) = M(i);
          }
        
        for(int n=0;n<nb_unknowns;n++)
          {
            for(int i=0;i<Dimension::dim_N;i++)
              {
                Dgrad_phi(m,n)(i) += grad_M(n,m)*M(i);
                //Ephi_grad(m,n)(i) += grad_M(m,n);
              }
          }
      }
    Dgrad_phi *= d*nat_mat.GetCoefStiffness();
    Ephi_grad *= nat_mat.GetCoefStiffness();
  }
  
  template<class T, class Dim>
  void GoldsteinEquation<T, Dim>::SetIndexToCompute(VarGalbrunIndex_Base<Dimension>& var)
  {
    var.compute_grad_flow = true;
    var.compute_source_momentum = false;
  }
  
  
  /***********************
   * GoldsteinEquationDG *
   ***********************/
  
  template<class T, class Dim>
  template<class TypeEquation, class T0, class Vector1>
  void GoldsteinEquationDG<T, Dim>::GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
                                                        const GlobalGenericMatrix<T0>& nat_mat,
                                                        Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    unknown_to_derive.Fill(true);
    fct_test_to_derive.Fill(true); 
  }
  
  
  //! computes mass = -i omega I - grad M
  template<class T, class Dim>
  template<class TypeEquation, class T0, class MatMass>
  void GoldsteinEquationDG<T, Dim>::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                                                int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    FillZero(mass);    
    MatrixN_N grad_M = vars.grad_flow(num_elem)(jloc);
    R_N M = vars.eval_flow(num_elem)(jloc);
    T omeg;
    vars.GetMiomega(omeg);
    for(int m = 0; m < nb_unknowns; m++)
      {
        mass(m, m) = omeg;
        for(int n = 0; n < nb_unknowns; n++)
          mass(m, n) -= grad_M(n, m);
      }
    
    mass *= nat_mat.GetCoefStiffness();    
  }
  
  
  //! 
  template<class T, class Dim>
  template<class TypeEquation, class T0, class MatStiff>
  void GoldsteinEquationDG<T,Dim>::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                                    int num_elem, int jloc,
                                                    const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                                    MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
    
    MatrixN_N grad_M = vars.grad_flow(num_elem)(jloc);
    R_N M = vars.eval_flow(num_elem)(jloc);
    T omeg;
    vars.GetMiomega(omeg);
    
    for (int m = 0; m < nb_unknowns; m++)
      {
        //Dgrad_phi(m,m) = omeg*M(m);
        //Ephi_grad(m,m) = omeg;

        for(int i = 0; i < Dimension::dim_N; i++)
          Ephi_grad(m,m)(i) = M(i);
      }
    
    Ephi_grad *= nat_mat.GetCoefStiffness();
  }
  
  
  template<class T, class Dim>
  void GoldsteinEquationDG<T, Dim>::SetIndexToCompute(VarGalbrunIndex_Base<Dimension>& var)
  {
    var.compute_grad_flow = true;
    var.compute_source_momentum = false;
  }
  
  
  template<class T, class Dim>
  template<class Matrix1, class GenericPb, class T0>
  void GoldsteinEquationDG<T, Dim>
  ::GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int k,
               int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
               const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    T0 v0_dot_n = DotProd(normale, vars.eval_flow(iquad)(k));      
    T0 coef = -abs(v0_dot_n)*nat_mat.GetCoefDamping();
    Nabc.SetDiagonal(coef);
  }

  
  template<class T, class Dim>
  template<class Matrix1, class GenericPb, class T0>
  void GoldsteinEquationDG<T, Dim>
  ::GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
            int ref, int iquad, int k,
            const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
            const GenericPb& vars, const ElementReference<Dimension, 1>&)
  {
    T0 v0_dot_n = DotProd(normale, vars.eval_flow(iquad)(k));      
    T0 coef = abs(v0_dot_n)*nat_mat.GetCoefDamping();
    Nabc.SetDiagonal(coef);
  }
  
  
  /*******************
   * VarGoldstein_Eq *
   *******************/

  //template<class Complexe, class Dimension> template<class TypeEquation>
  template<class TypeEquation>
  VarGoldstein_Eq<TypeEquation>::VarGoldstein_Eq(EllipticProblem<TypeEquation>& var) : VarGalbrunIndex_Base<typename TypeEquation::Dimension>(var)
  {
    //    adjustment_neumann = false;
    //    compute_div_flow = false;
    //    compute_grad_rho = false;
    //    compute_grad_sigma = false;
    //    compute_grad_c0 = false;
    //    compute_hessian_flow = false;
    //    compute_hessian_p0 = false;
    //    compute_gamma = false;
    //    compute_grad_gamma = false;
    //    store_grad_rho0_c0 = false;
  }
  
  //  template<class TypeEquation>
  //    VarGoldstein_Eq<TypeEquation>::VarGoldstein_Eq() 
  //    {
  //
  //    }
  template<class TypeEquation>
  VarGoldstein_Eq<TypeEquation>::VarGoldstein_Eq()
    : VarGalbrunIndex_Base<typename TypeEquation::Dimension>(this->GetLeafClass()),VarHarmonic<TypeEquation>()
  {
    
  }
  
  template<class TypeEquation>
  void VarGoldstein_Eq<TypeEquation>::ComputePhysicalCoefficients()
  {
    TypeEquation::SetIndexToCompute(*this);
    
    VarHarmonic<TypeEquation>::ComputePhysicalCoefficients();
    VarGalbrunIndex_Base<Dimension>::ComputePhysicalCoefficients(); 
  }

  
  template<class TypeEquation>
  void VarGoldstein_Eq<TypeEquation>::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("PenalSUPG"))
      {
        if (parameters.GetM() <= 0)
          {
            std::cout << "In SetInputData of PenalSUPG :" << std::endl;
            std::cout << "Error ! PenalSUPG need one argument" << std::endl;
            abort();
          }
        
        penal_supg = to_num<Real_wp>(parameters(0));
      }
    else
      {
        VarGalbrunIndex_Base<Dimension>::SetInputData(description_field,parameters);
        VarHarmonic<TypeEquation>::SetInputData(description_field,parameters); 
      }
  }
  
  
  /*******************
   * EllipticProblem *
   *******************/  
  

  template<class Dimension>
  void EllipticProblem<GoldsteinEquation<Complex_wp,Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
                            const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i,num_dof,mat_elem,nat_mat,*this,this->GetReferenceElementH1(i));
  }
  

  template<class Dimension>
  void EllipticProblem<GoldsteinEquationDG<Complex_wp,Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
                            const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this, this->GetReferenceElementH1(i));
  }

  template<class Dimension>
  void EllipticProblem<GoldsteinEquation<Complex_wp,Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
                          const GlobalGenericMatrix<Complex_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
  }
  
  template<class Dimension>
  void EllipticProblem<GoldsteinEquationDG<Complex_wp,Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
                          const GlobalGenericMatrix<Complex_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);   
  }

}

#define MONTJOIE_FILE_HARMONIC_GOLDSTEIN_EQUATION_CXX
#endif
