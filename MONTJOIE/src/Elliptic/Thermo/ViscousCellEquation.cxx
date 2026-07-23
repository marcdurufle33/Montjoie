#ifndef MONTJOIE_FILE_VISCO_THERMAL_EQUATION_HXX

#include "ViscousCellEquation.hxx"

// #define OTHER_FORMULATION_VCELL

namespace Montjoie
{

  /************************
   * ViscousCellEquation *
   ************************/

  // providing C in Cgrad_grad
  template<class T>
  template<class TypeEquation, class NatureMat, class MatStiff>
  void ViscousCellEquation<T>::GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                                  int num_elem, int jloc, const NatureMat& nat_mat,
                                                  int ref, MatStiff& Cgrad)
  {
    FillZero(Cgrad);
    T coef = -vars.viscosity*nat_mat.GetCoefStiffness();
    Cgrad(0, 0)(0, 0) = coef;
    Cgrad(0, 0)(1, 1) = coef;
    Cgrad(1, 1)(0, 0) = coef;
    Cgrad(1, 1)(1, 1) = coef;
  }


  template<class T>
  template<class TypeEquation, class T0, class MatStiff>
  void ViscousCellEquation<T>::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                                 int num_elem, int jloc,
                                                 const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    Dgrad_phi.Zero(); Ephi_grad.Zero();
#ifndef OTHER_FORMULATION_VCELL
    // Formulation div-div
    Dgrad_phi(0, 2)(0) = 1;
    Dgrad_phi(1, 2)(1) = 1;
    Ephi_grad(2, 0)(0) = 1;
    Ephi_grad(2, 1)(1) = 1;
#else
    // Formulation grad-grad
    Ephi_grad(0, 2)(0) = -1;
    Ephi_grad(1, 2)(1) = -1;
    Dgrad_phi(2, 0)(0) = -1;
    Dgrad_phi(2, 1)(1) = -1;
#endif
  }

  template<class T>
  template<class TypeEquation, class T0, class MatMass>
  void ViscousCellEquation<T>::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                                              int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    mass.Zero();
    Complex_wp m_iomega; vars.GetMiomega(m_iomega); // m_iomega = - i omega
    // Complex_wp m_omega2; vars.GetMomega2(m_omega2);
    T coef = m_iomega*nat_mat.GetCoefMass() * vars.density;
    // T coef = -m_iomega*density;
    mass(0, 0) = coef;
    mass(1, 1) = coef;
    Real_wp epsilon = 1e-6; // Avoid having a singular mass matrix
    mass(2, 2) = epsilon*coef;
  }


  /**********************
   * VarViscousCell_Eq *
   **********************/


  template<class TypeEquation>
  void VarViscousCell_Eq<TypeEquation>::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
  }

  template<class TypeEquation>
  void VarViscousCell_Eq<TypeEquation>::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    if (name_media == "viscosity")
    {
      viscosity = to_num<Real_wp>(parameters(0));
      cout << "VarViscousCell_Eq: Set viscosity to " << viscosity << endl;
    }
    else if (name_media == "density")
    {
      density = to_num<Real_wp>(parameters(0));
      cout << "VarViscousCell_Eq: Set density to " << density << endl;
    }
  }

  template<class TypeEquation>
  string VarViscousCell_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return string("");
  }


  /**********************
   * EllipticProblem *
   **********************/
  //
  //
  void EllipticProblem<ViscousCellEquation<Complex_wp> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
                            const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
                                      this->GetReferenceElementH1(i));
  }
  //
  void EllipticProblem<ViscousCellEquation<Complex_wp> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
                          const GlobalGenericMatrix<Complex_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }

}

#define MONTJOIE_FILE_VISCO_THERMAL_EQUATION_HXX
#endif
