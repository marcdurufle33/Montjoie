#ifndef MONTJOIE_FILE_VISCO_THERMAL_EQUATION_HXX

namespace Montjoie
{

  // visco-thermal equation
  template<class Dim, class T>
  class ViscoThermalEquation : public GenericEquation<T>
  {
  public:
    typedef Dim Dimension;

    enum{nb_unknowns = 2+Dimension::dim_N, nb_components_en = 1, nb_components_hn = 1,
         nb_unknowns_scal = 2+Dimension::dim_N, nb_unknowns_vec = 0 , nb_unknowns_hdg = 0};

    static const bool FormulationDG = false;
    static const bool FirstOrderFormulation = false;
    static const bool DiscontinuousDiMatrix = false;
    static const bool TensorStiffnessSymmetric = false;

    static inline IVect GetOtherElementType() { IVect num(1); num(0) = 1; return num; }

    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return true; }

    // providing C in Cgrad_grad
    template<class TypeEquation, class NatureMat, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                  int num_elem, int jloc, const NatureMat& nat_mat,
                                  int ref, MatStiff& Cgrad_grad);

    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc,
                                 const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);

    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                              int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);

  };


  //------------------------
  // ViscoThermalIndex
  //------------------------
  
  class ViscoThermalIndex
  {
  public:
    Real_wp c_p, mu, lambda, T0, rho0, P0, zeta; // lambda = kappa
    Real_wp Temp, h_v, x_c, x_v;
    Real_wp h_v0, x_c0, x_v0, gamma, c_0;
    
    inline ViscoThermalIndex() {h_v0 = 0.5; x_c0 = 4.2e-4; x_v0 = 1.1571e-2; T0 = 293.15; Temp = 1.0; h_v = 1.0; x_c = 1.0; x_v = 1.0; };
 
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void UpdateTemp(Real_wp Temp, Real_wp h_v, Real_wp x_c);
  };


  //!
  template<class TypeEquation>
  class VarViscoThermal_Eq:public ViscoThermalIndex, public VarHarmonic<TypeEquation>
  {
    typedef Complex_wp Complexe;
    typedef typename TypeEquation::Dimension Dimension;

  public:
    // Real_wp c_p, T0, P0, rho0, mu, lambda, zeta;
    bool alternative_ipp;

  public:
    inline VarViscoThermal_Eq() { c_p = 1.0; P0 = 1.0; rho0 = 1.0; mu = 1.0; lambda = 1.0; zeta = 1.0; alternative_ipp = false; }


    void SetInputData(const string& description_field, const VectString& parameters);

    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;

    inline void InitIndices(int n) {}
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
  template<class Dimension>
  class EllipticProblem<ViscoThermalEquation<Dimension, Complex_wp> >
    : public VarViscoThermal_Eq<ViscoThermalEquation<Dimension, Complex_wp> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                                 CondensationBlockSolver_Base<Real_wp>&,
                                 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                                 CondensationBlockSolver_Base<Complex_wp>&,
                                 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
  };



  
}

#define MONTJOIE_FILE_VISCO_THERMAL_EQUATION_HXX
#endif
