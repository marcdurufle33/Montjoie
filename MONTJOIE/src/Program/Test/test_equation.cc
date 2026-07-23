#include "Montjoie.hxx"

using namespace Montjoie;

#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"

namespace Montjoie
{

  // visco-thermal equation
  template<class T>
  class MyEquation : public GenericEquation<T>
  {
  public:
    typedef Dimension2 Dimension;

    enum{nb_unknowns = 2, nb_components_en = 1, nb_components_hn = 1,
         nb_unknowns_scal = 2, nb_unknowns_vec = 0 , nb_unknowns_hdg = 0};

    static const bool FormulationDG = false;
    static const bool FirstOrderFormulation = false;
    static const bool DiscontinuousDiMatrix = false;
    static const bool TensorStiffnessSymmetric = false;

    inline static bool SymmetricGlobalMatrix() { return false; } 
    inline static bool SymmetricElementaryMatrix() { return false; } 

    static inline IVect GetOtherElementType() { IVect num(1); num(0) = 1; return num; }
    
    // providing C in Cgrad_grad
    template<class TypeEquation, class NatureMat, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                  int num_elem, int jloc, const NatureMat& nat_mat,
                                  int ref, MatStiff& Cgrad_grad)
    {
      Cgrad_grad.Zero();

      Cgrad_grad(0, 0)(0, 0) = -vars.lambda/vars.T0;
      Cgrad_grad(0, 0)(1, 1) = -vars.lambda/vars.T0;

      Complex_wp m_iomega; vars.GetMiomega(m_iomega);
      Complex_wp coef;
      coef = Real_wp(1) / (m_iomega*vars.rho0);
      Cgrad_grad(1, 1)(0, 0) = coef;
      Cgrad_grad(1, 1)(1, 1) = coef;
      
    }

    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc,
                                 const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
    {

      Dgrad_phi.Zero(); Ephi_grad.Zero();
      
      // equations des ondes en formulation mixte
      /*
      Dgrad_phi(0, 1)(0) = nat_mat.GetCoefStiffness();
      Dgrad_phi(0, 2)(1) = nat_mat.GetCoefStiffness();

      Ephi_grad(1, 0)(0) = -nat_mat.GetCoefStiffness();
      Ephi_grad(2, 0)(1) = -nat_mat.GetCoefStiffness(); */
    }

    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                              int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
    {
      Complex_wp m_iomega; vars.GetMiomega(m_iomega);
      mass.Zero();
      mass(0, 0) = m_iomega*vars.rho0*vars.cp/vars.T0;
      mass(0, 1) = -m_iomega/vars.T0;
      mass(1, 0) = -m_iomega/vars.T0;
      mass(1, 1) = m_iomega/vars.P0;
      
      /*
        equations des ondes en formulation mixte
      Complex_wp m_iomega; vars.GetMiomega(m_iomega);
      T0 coef = m_iomega*nat_mat.GetCoefMass();
      mass(0,0) = -coef/vars.P0;
      mass(1,1) = -coef*vars.rho0;
      mass(2,2) = -coef*vars.rho0;
      */
    }
    
  };


  //!
  template<class TypeEquation>
  class VarMyEquation_Eq : public VarHarmonic<TypeEquation>
  {
    typedef Complex_wp Complexe;
    typedef Dimension2 Dimension;

  public:
    Real_wp rho0, P0, cp, lambda, T0;

    VarMyEquation_Eq()
    { rho0 = 1.0; P0 = 1.0; cp = 1.0; lambda = 1.0; T0 = 1.0; }

    void SetInputData(const string& description_field, const VectString& parameters)
    {
      VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
    }

    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
    {
      if (name_media == "P0")
        P0 = to_num<Real_wp>(parameters(1));
      else if (name_media == "rho0")
        rho0 = to_num<Real_wp>(parameters(1));
      else if (name_media == "T0")
        T0 = to_num<Real_wp>(parameters(1));
      else if (name_media == "lambda")
        lambda = to_num<Real_wp>(parameters(1));
      else if (name_media == "cp")
        cp = to_num<Real_wp>(parameters(1));
    }
    
    string GetPhysicalIndexName(int m) const
    {
      return string("");
    }

    inline void InitIndices(int n) { }
    inline int GetNbPhysicalIndices() const{ return 100;}
    inline void SetIndices(int i, const VectString& parameters) {}
    inline bool IsVaryingMedia(int i) const { return false; }
    inline bool IsVaryingMedia(int m, int i) const { return false; }
    inline Real_wp GetVelocityOfMedia(int ref) const{ return 1.0;}
    inline Real_wp GetVelocityOfInfinity() const{return 1.0;}

    inline void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
                                  IVect& num_ref, IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                                  Vector<bool>& compute_hess) {}

  };
  


  //---------------------------
  // EllipticProblem
  //---------------------------
  template<>
  class EllipticProblem<MyEquation<Complex_wp> >
    : public VarMyEquation_Eq<MyEquation<Complex_wp> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                                 CondensationBlockSolver_Base<Complex_wp>&,
                                 const GlobalGenericMatrix<Complex_wp>& nat_mat)
    {
      Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
                                        this->GetReferenceElementH1(i));
    }
    
    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0)
    {
    }
    
  };

}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      string input_file(argv[1]), type_element, type_equation;
      getElement_Equation(input_file, type_element, type_equation);
      
      EllipticProblem<MyEquation<Complex_wp> > vars;
      vars.RunAll(input_file, type_element, type_equation);
    }
  else
    {
      std::cout << "Please provide a data file" << std::endl;
      abort();
    }
  
  return FinalizeMontjoie();
}
