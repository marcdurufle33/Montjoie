#ifndef MONTJOIE_FILE_VISCOUS_CELL_EQUATION_HXX

#include "Montjoie.hxx"

namespace Montjoie
{

  // visco-thermal equation
  template<class T>
  class ViscousCellEquation : public GenericEquation<T>
  {
  public:
    typedef Dimension2 Dimension;

    enum{nb_unknowns = 3, nb_components_en = 3, nb_components_hn = 3,
         nb_unknowns_scal = 3, nb_unknowns_vec = 0, nb_unknowns_hdg = 0};
    
    static const bool FormulationDG = false;
    static const bool FirstOrderFormulation = false;
    static const bool DiscontinuousDiMatrix = false;
    static const bool TensorStiffnessSymmetric = false;

    inline static bool SymmetricGlobalMatrix() { return true; } 
    inline static bool SymmetricElementaryMatrix() { return true; } 

    static inline IVect GetOtherElementType() { IVect num(1); num(0) = 1; return num; }
    
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


  //!
  template<class TypeEquation>
  class VarViscousCell_Eq : public VarHarmonic<TypeEquation>
  {
    typedef Complex_wp Complexe;
    typedef Dimension2 Dimension;

  public:
    Real_wp viscosity, density;

  public:
    inline VarViscousCell_Eq() { viscosity = 1.0; density = 1.0; }


    void SetInputData(const string& description_field, const VectString& parameters);

    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;

    inline void InitIndices(int n) {}
    // inline int GetNbPhysicalIndices() const{ return 2; }
    int GetNbPhysicalIndices() const{ return 2; }
    inline void SetIndices(int i, const VectString& parameters) {}
    inline bool IsVaryingMedia(int i) const { return false; }
    inline bool IsVaryingMedia(int m, int i) const { return false; }
    // inline Real_wp GetVelocityOfMedia(int ref) const{ return density;} // FIXME
    // inline Real_wp GetVelocityOfInfinity() const{return rho;}

    inline void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
                                  IVect& num_ref, IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                                  Vector<bool>& compute_hess) {}

  };



  //---------------------------
  // EllipticProblem
  //---------------------------
  template<>
  class EllipticProblem<ViscousCellEquation<Complex_wp> >
    : public VarViscousCell_Eq<ViscousCellEquation<Complex_wp> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                                 CondensationBlockSolver_Base<Complex_wp>&,
                                 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
  };

}

#define MONTJOIE_FILE_VISCOUS_CELL_EQUATION_HXX
#endif
