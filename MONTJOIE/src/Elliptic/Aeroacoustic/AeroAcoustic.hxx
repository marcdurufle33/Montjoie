#ifndef MONTJOIE_FILE_AERO_ACOUSTIC_HXX

namespace Montjoie
{
  // file used to define aero-acoustic classes for the resolution
  // in time-harmonic (or stationary)
  // definition of methods and functions are in file AeroAcoustic.cxx
  
  template<class Dimension> class AeroAcoustic_Base;
  
  //! class to specify resolution of linearized Euler equations
  /*!
    (-i omega + sigma + M \cdot \nabla) \rho + rho div(M) + div u = 0
    (-i omega + sigma + M \cdot \nabla) u - M \cdot \nabla rho_0 / rho_0 u + \nabla p + \nabla M (u + rho M) = 0
    (-i omega + sigma + M \cdot \nabla) p +  div(c_0^2 u) + gamma div(M) p - (gamma-1)/rho0 u \cdot \nabla p0
                                          - p0 / rho_0 \nabla \gamma \cdot u  = 0
    
    rho0, p0, c0, M are given coefficients (stationary flow)
    gamma is given by gamma = rho0 c0^2 / p0    
    unknowns are (p, \rho, u) = (p', \rho', \rho' u')
    where p', \rho', u' are the classical perturbations of Euler system
   */
  template<class T, class Dim>
  class LinearizedEulerEquation_Base : public GenericEquation<T>
  {
  public :
    typedef Dim Dimension;
    
    static const bool FirstOrderFormulation = true;
    
    enum {nb_unknowns = 2+Dimension::dim_N, nb_unknowns_hdg=0,
	  nb_components_en = 1, nb_components_hn = 1,
          nb_unknowns_scal = 2+Dimension::dim_N, nb_unknowns_vec = 0};

    static inline bool SymmetricGlobalMatrix() { return false; }
    static inline bool SymmetricElementaryMatrix() { return false; }

    static void SetIndexToCompute(AeroAcoustic_Base<Dimension>& var);
        
    // for compatbility purpose
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int num_elem, const ElementReference_Dim<Dimension>& Fb);
        
    // For a detailed description of the following methods
    // look at the file GenericEquation.hxx (class GenericEquation_Base)
    template<class GenericPb, class T0, class Vector1>
    static void GetNeededDerivative(const GenericPb& vars, const GlobalGenericMatrix<T0>& nat_mat,
				    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class T1, class T0>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                                     int ref, const TinyVector<TinyVector<T0, 2>, 4>& dU,
                                     TinyVector<T0, 4>& V);

    template<class TypeEquation, class T1, class T0>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                                     int ref, const TinyVector<TinyVector<T0, 3>, 5>& dU,
                                     TinyVector<T0, 5>& V);

    template<class TypeEquation, class T1, class T0>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                                     int ref, const TinyVector<T0, 4>& U,
                                     TinyVector<TinyVector<T0, 2>, 4>& dV);

    template<class TypeEquation, class T1, class T0>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                                     int ref, const TinyVector<T0, 5>& U,
                                     TinyVector<TinyVector<T0, 3>, 5>& dV);
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T1, class T0>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T1>& nat_mat, int ref,
                                const TinyVector<T0, 4>& U, TinyVector<T0, 4>& V);

    template<class TypeEquation, class T1, class T0>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T1>& nat_mat, int ref,
                                const TinyVector<T0, 5>& U, TinyVector<T0, 5>& V);

    template<class T0>
    static void GetAbsoluteD(TinyMatrix<T0, General, 4, 4>& Dtest, const R2& normale,
                             const Real_wp& c0);

    template<class T0>
    static void GetAbsoluteD(TinyMatrix<T0, General, 5, 5>& Dtest, const R3& normale,
                             const Real_wp& c0);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
			int ref, int iquad, int k,
			const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
                        const GenericPb& vars, const ElementReference<Dimension, 1>&);
    
    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(typename Dimension::R_N& normale, int ref, const Vector1& Vn, Vector1& Un,
                        int num_elem1, int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
                        const EllipticProblem<TypeEquation>& vars,
			const ElementReference<Dimension, 1>& Fb);

    template<class Matrix1, class GenericPb, class T0>
    static void GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int k,
                           int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class T0, class GenericPb, class T1>
    static void MltPenalDG(const typename Dimension::R_N& normale,
                           const TinyVector<T0, 4>& Vn, TinyVector<T0, 4>& Un,
                           int i, int k, int nf, const GlobalGenericMatrix<T1>& nat_mat,
                           int ref, int ref2, const GenericPb& vars,
			   const ElementReference<Dimension, 1>& Fb);

    template<class T0, class GenericPb, class T1>
    static void MltPenalDG(const typename Dimension::R_N& normale,
                           const TinyVector<T0, 5>& Vn, TinyVector<T0, 5>& Un,
                           int i, int k, int nf, const GlobalGenericMatrix<T1>& nat_mat,
                           int ref, int ref2, const GenericPb& vars,
			   const ElementReference<Dimension, 1>& Fb);
    
  };
  
  
  //! stationary aeroacoustic equation for uniform flows
  template<class Dim>
  class StationaryLinearizedEulerEquation
    : public LinearizedEulerEquation_Base<Real_wp, Dim>
  {
  };

  //! time-harmonic aeroacoustic equation for uniform flows
  template<class Dim>
  class HarmonicLinearizedEulerEquation
    : public LinearizedEulerEquation_Base<Complex_wp, Dim>
  {
  };


  template<class Dimension>
  class FemMatrixFreeClass<Real_wp, StationaryLinearizedEulerEquation<Dimension> >
    : public FemMatrixFreeClass_Eq<Real_wp, StationaryLinearizedEulerEquation<Dimension> >
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<StationaryLinearizedEulerEquation<Dimension> >& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };


  //! class to specify resolution of aeroacoustic equation for any flow
  /*!
    Implementation of approximate models of Linearized Euler Equations
    Simplified model :
    (-i omega + sigma + M \cdot \nabla) p + div(c0^2 u) = 0
    (-i omega + sigma + M \cdot \nabla) u + \nabla p  + (\nabla M) u = 0
    
    Bogey-Bailly-Juve :
    (-i omega + sigma + M \cdot \nabla) p + div(c0^2 u) = 0
    (-i omega + sigma + M \cdot \nabla) u + \nabla p = 0
    
    Conservative model :
    rho0 (-i omega + sigma + M \cdot \nabla) p + div(rho0 c0 v) = 0
    rho0 (-i omega + sigma + M \cdot \nabla) u + rho0 c0 \nabla p = 0
    
    Galbrun model
    \rho (-i omega + sigma + M \cdot \nabla) p + (\rho c)^2 div u = 0
    \rho (-i omega + sigma + M \cdot \nabla) u + \nabla p + ((div u) - (\nabla u)^T) \nabla p_0 / (-i omega + \sigma) = 0
    
    Bogey-Bailly-Juve and conservative model avoid Kelvin-Helmholtz instabilities
    while the simplified model should exhibit instabilities
   */
  template<class T, class Dim>
  class AeroStationaryEquation_Base : public GenericEquation<T>
  {
  public :
    typedef Dim Dimension;

    static const bool FirstOrderFormulation = true;
    
    enum {nb_unknowns = 1+Dimension::dim_N, nb_unknowns_hdg=0,
	  nb_components_en = 1, nb_components_hn = 1,
          nb_unknowns_scal = 1+Dimension::dim_N, nb_unknowns_vec = 0};

    static inline bool SymmetricGlobalMatrix() { return false; }
    static inline bool SymmetricElementaryMatrix() { return false; }

    static void SetIndexToCompute(AeroAcoustic_Base<Dimension>& var);
    
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int num_elem, const ElementReference_Dim<Dimension>& Fb);
    
    // For a detailed description of the following methods
    // look at the file GenericEquation.hxx (class GenericEquation_Base)
    template<class GenericPb, class T0, class Vector1>
    static void GetNeededDerivative(const GenericPb& vars,
                                    const GlobalGenericMatrix<T0>& nat_mat,
				    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class T1, class T0>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                                     int ref, const TinyVector<TinyVector<T0, 2>, 3>& dU, TinyVector<T0, 3>& V);

    template<class TypeEquation, class T1, class T0>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                                     int ref, const TinyVector<TinyVector<T0, 3>, 4>& dU, TinyVector<T0, 4>& V);

    template<class TypeEquation, class T1, class T0>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                                     int ref, const TinyVector<T0, 3>& U, TinyVector<TinyVector<T0, 2>, 3>& dV);

    template<class TypeEquation, class T1, class T0>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                                     int ref, const TinyVector<T0, 4>& U, TinyVector<TinyVector<T0, 3>, 4>& dV);

    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T1, class T0>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T1>& nat_mat, int ref,
                                const TinyVector<T0, 3>& U, TinyVector<T0, 3>& V);

    template<class TypeEquation, class T1, class T0>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T1>& nat_mat, int ref,
                                const TinyVector<T0, 4>& U, TinyVector<T0, 4>& V);

    template<class T0>
    static void GetAbsoluteD(TinyMatrix<T0, General, 3, 3>& Dtest, const R2& normale,
                             const Real_wp& c0, bool conservative, bool galbrun, const T0& gamma_);
    
    template<class T0>
    static void GetAbsoluteD(TinyMatrix<T0, General, 4, 4>& Dtest, const R3& normale,
                             const Real_wp& c0, bool conservative, bool galbrun,
                             const T0& gamma_);

    template<class Matrix1, class GenericPb, class T0>
    static void GetAbsoluteD(Matrix1& Nabc, const typename Dimension::R_N& normale, int iquad, int k,
                             const GlobalGenericMatrix<T0>& nat_mat, const GenericPb& vars);

    template<class Matrix1, class GenericPb, class T0>
    static void GetAbsoluteMatrixD(Matrix1& Nabc, const typename Dimension::R_N& normale, int iquad, int k,
                                   const GlobalGenericMatrix<T0>& nat_mat, const GenericPb& vars);
        
    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
			int ref, int iquad, int k,
			const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
                        const GenericPb& vars, const ElementReference<Dimension, 1>&);
    
    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(typename Dimension::R_N& normale, int ref, const Vector1& Vn, Vector1& Un,
                        int num_elem1, int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
                        const EllipticProblem<TypeEquation>& vars,
			const ElementReference<Dimension, 1>& Fb);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int k,
                           int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class T0, class GenericPb, class T1>
    static void MltPenalDG(const typename Dimension::R_N& normale,
                           const TinyVector<T0, 3>& Vn, TinyVector<T0, 3>& Un,
                           int i, int k, int nf, const GlobalGenericMatrix<T1>& nat_mat,
                           int ref, int ref2, const GenericPb& vars,
			   const ElementReference<Dimension, 1>& Fb);

    template<class T0, class GenericPb, class T1>
    static void MltPenalDG(const typename Dimension::R_N& normale,
                           const TinyVector<T0, 4>& Vn, TinyVector<T0, 4>& Un,
                           int i, int k, int nf, const GlobalGenericMatrix<T1>& nat_mat,
                           int ref, int ref2, const GenericPb& vars,
			   const ElementReference<Dimension, 1>& Fb);
    
  };
  
  
  //! stationary aeroacoustic equation for any flow
  template<class Dim>
  class AeroStationaryEquation : public AeroStationaryEquation_Base<Real_wp, Dim>
  {
  public :
  };

  //! time-harmonic aeroacoustic equation for any flow
  template<class Dim>
  class HarmonicAeroEquation : public AeroStationaryEquation_Base<Complex_wp, Dim>
  {
  public :
  };


  template<class Dimension>
  class FemMatrixFreeClass<Real_wp, AeroStationaryEquation<Dimension> >
    : public FemMatrixFreeClass_Eq<Real_wp, AeroStationaryEquation<Dimension> >
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<AeroStationaryEquation<Dimension> >& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };


  //! stationary or time-harmonic base class for aeroacoustic
  template<class Dimension>
  class AeroAcoustic_Base : public VarGalbrunIndex_Base<Dimension>
  {
  public :
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    
    enum {LEE_MODEL, SIMPLIFIED_LEE, BOGEY_BAILLY_JUVE, CONSERVATIVE, GALBRUN};
    int type_model;
    
  private:
    DistributedProblem<Dimension>& var_problem;
    
  public:
    template<class TypeEquation>
    AeroAcoustic_Base(EllipticProblem<TypeEquation>&);
    
    // additional parameters of the data file
    void SetInputData(const string& description_field, const VectString& parameters);
        
    template<class T0>
    void ModifyVolumetricSource(int i, int j, const R_N& x,
				const VirtualSourceField<T0, Dimension>& fsrc,
                                Vector<T0>& f) const;
    
  };


  template<class TypeEquation>
  class AeroAcoustic_Eq : public VarHarmonic<TypeEquation>,
			  public AeroAcoustic_Base<typename TypeEquation::Dimension>
  {
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeEquation::Dimension Dimension;
    
  public:
    AeroAcoustic_Eq();
   
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void InitIndices(int n);
    int GetNbPhysicalIndices() const;
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
			   IVect& num_ref, IVect& num_index, IVect& num_component,
			   Vector<bool>& compute_grad, Vector<bool>& compute_hess);

    void ComputePhysicalCoefficients();

    void PerformOtherInitializations();

    template<class T, class Dim>
    void ModifyVolumetricSourceGen(int i, int j, const typename Dim::R_N&,
                                   const VirtualSourceField<T, Dim>&,
                                   Vector<T>&) const;

    void ModifyVolumetricSourceGen(int i, int j, const typename Dimension::R_N&,
                                   const VirtualSourceField<Complexe, Dimension>&,
                                   Vector<Complexe>&) const;
    
    void ModifyVolumetricSource(int i, int j, const R2&,
                                const VirtualSourceField<Real_wp, Dimension2>&,
                                Vector<Real_wp>&) const;
    
    void ModifyVolumetricSource(int i, int j, const R2&,
                                const VirtualSourceField<Complex_wp, Dimension2>&,
                                Vector<Complex_wp>&) const;

    void ModifyVolumetricSource(int i, int j, const R3&,
                                const VirtualSourceField<Real_wp, Dimension3>&,
                                Vector<Real_wp>&) const;
    
    void ModifyVolumetricSource(int i, int j, const R3&,
                                const VirtualSourceField<Complex_wp, Dimension3>&,
                                Vector<Complex_wp>&) const;
    
  };


  //! class used to solve harmonic aeroacoustic equation with DG method
  template<class Dimension>
  class EllipticProblem<HarmonicLinearizedEulerEquation<Dimension> >
    : public AeroAcoustic_Eq<HarmonicLinearizedEulerEquation<Dimension> >
  {
  public:
    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Real_wp>&,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>&);

    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Complex_wp>&,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>&);

  };

  //! class used to solve real aeroacoustic equation with DG method
  template<class Dimension>
  class EllipticProblem<StationaryLinearizedEulerEquation<Dimension> >
    : public AeroAcoustic_Eq<StationaryLinearizedEulerEquation<Dimension> >
  {
  public:
  };

  //! class used to solve harmonic aeroacoustic equation with DG method
  template<class Dimension>
  class EllipticProblem<HarmonicAeroEquation<Dimension> >
    : public AeroAcoustic_Eq<HarmonicAeroEquation<Dimension> >
  {
  public:
    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Real_wp>&,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>&);

    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Complex_wp>&,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>&);
    
  };

  //! class used to solve real aeroacoustic equation with DG method
  template<class Dimension>
  class EllipticProblem<AeroStationaryEquation<Dimension> >
    : public AeroAcoustic_Eq<AeroStationaryEquation<Dimension> >
  {
  public:
    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Real_wp>&,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>&);

    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Complex_wp>&,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>&);

  };
  
}

#define MONTJOIE_FILE_AERO_ACOUSTIC_HXX
#endif
