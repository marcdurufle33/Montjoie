#ifndef MONTJOIE_FILE_HARMONIC_GOLDSTEIN_EQUATION_HXX

namespace Montjoie
{

  template<class Complexe, class Dimension>
  class VarGoldstein_Dim;
  
  //---------------------------
  //GoldsteinEquation
  //---------------------------
  template<class T, class Dim>
  class GoldsteinEquation : public GenericEquation<T>
  {
  public:
    typedef Dim Dimension;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename Dimension::R_N R_N;
    
    enum{nb_unknowns = Dimension::dim_N, nb_components_en = Dimension::dim_N, nb_components_hn = Dimension::dim_N,
         nb_unknowns_scal = Dimension::dim_N, nb_unknowns_vec = 0 , nb_unknowns_hdg=0};
    
    
    static const bool FormulationDG = false;
    static const bool FirstOrderFormulation = true;
    static const bool DiscontinuousDiMatrix = false;
    static const bool TensorStiffnessSymmetric = false;
    
    
    // providing C in Cgrad_grad
    template<class TypeEquation, class NatureMat, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                  int num_elem, int jloc, const NatureMat& nat_mat,
                                  int ref, MatStiff& Cgrad_grad);
    
    template<class TypeEquation, class T0, class Vector1>
    static void GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
                                    const GlobalGenericMatrix<T0>& nat_mat,
                                    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc,
                                 const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                              int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    static void SetIndexToCompute(VarGalbrunIndex_Base<Dimension>& var);
    
  };


  // ---------------------------
  // GoldsteinEquationDG
  // ---------------------------
  
  template<class T, class Dim>
  class GoldsteinEquationDG : public GenericEquation<T>
  {
  public:
    typedef Dim Dimension;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename Dimension::R_N R_N;
    
    enum {nb_unknowns = Dimension::dim_N, nb_components_en = Dimension::dim_N, nb_components_hn = Dimension::dim_N,
          nb_unknowns_scal = Dimension::dim_N, nb_unknowns_vec = 0 , nb_unknowns_hdg=0};
    
    static const bool FormulationDG = true;
    static const bool FirstOrderFormulation = true;
    static const bool DiscontinuousDiMatrix = false;
    static const bool TensorStiffnessSymmetric = false;
    
    template<class TypeEquation, class T0, class Vector1>
    static void GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
                                    const GlobalGenericMatrix<T0>& nat_mat,
                                    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc,
                                 const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                              int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int k,
                           int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);

    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
			int ref, int iquad, int k,
			const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
                        const GenericPb& vars, const ElementReference<Dimension, 1>&);
    
    static void SetIndexToCompute(VarGalbrunIndex_Base<Dimension>& var);
    
  };
  

  //---------------------------
  //VarGoldStein_eq
  //---------------------------


  
  template<class TypeEquation>
  class VarGoldstein_Eq : public VarGalbrunIndex_Base<typename TypeEquation::Dimension>, public VarHarmonic<TypeEquation>
  {
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeEquation::Dimension Dimension;
    
  public:
    VarGoldstein_Eq();
    VarGoldstein_Eq(EllipticProblem<TypeEquation>& var);
    
    Real_wp penal_supg = 0.0;
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    inline void InitIndices(int n){VarGalbrunIndex_Base<Dimension>::InitIndices(n); }
    inline int GetNbPhysicalIndices() const{return VarGalbrunIndex_Base<Dimension>::GetNbPhysicalIndices();}
    inline void SetIndices(int i, const VectString& parameters){VarGalbrunIndex_Base<Dimension>::SetIndices(i,parameters);}
    inline void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters){VarGalbrunIndex_Base<Dimension>::SetPhysicalIndex(name_media,i,parameters);}
    inline string GetPhysicalIndexName(int m) const{return VarGalbrunIndex_Base<Dimension>::GetPhysicalIndexName(m);}
    inline bool IsVaryingMedia(int i) const{return VarGalbrunIndex_Base<Dimension>::IsVaryingMedia(i);}
    inline bool IsVaryingMedia(int m, int i) const{return VarGalbrunIndex_Base<Dimension>::IsVaryingMedia(i);}
    inline Real_wp GetVelocityOfMedia(int ref) const{return VarGalbrunIndex_Base<Dimension>::GetVelocityOfMedia(ref);}
    inline Real_wp GetVelocityOfInfinity() const{return VarGalbrunIndex_Base<Dimension>::GetVelocityOfInfinity();}

    inline void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
                                  IVect& num_ref, IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                                  Vector<bool>& compute_hess) {VarGalbrunIndex_Base<Dimension>::GetVaryingIndices(rho_real, num_ref, num_index, num_component, compute_grad, compute_hess);}    
    
    void ComputePhysicalCoefficients();
  };
  
  

  //---------------------------
  // EllipticProblem
  //---------------------------

  template<class Dimension >
  class EllipticProblem<GoldsteinEquation<Complex_wp,Dimension> >
    : public VarGoldstein_Eq<GoldsteinEquation<Complex_wp,Dimension> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                                 CondensationBlockSolver_Base<Complex_wp>&,
                                 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
  };


  template<class Dimension >
  class EllipticProblem<GoldsteinEquationDG<Complex_wp, Dimension> >
    : public VarGoldstein_Eq<GoldsteinEquationDG<Complex_wp, Dimension> >
  {
  public:
    //EllipticProblem() { DISP((GoldsteinEquationDG<Complex_wp, Dimension>::FirstOrderFormulation)); DISP(this->dg_formulation); }
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                                 CondensationBlockSolver_Base<Complex_wp>&,
                                 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
  };
  
  
  
  
}

#define MONTJOIE_FILE_HARMONIC_GOLDSTEIN_EQUATION_HXX
#endif
