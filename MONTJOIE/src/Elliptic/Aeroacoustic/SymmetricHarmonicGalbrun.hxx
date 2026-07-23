#ifndef MONTJOIE_FILE_SYMMETRIC_HARMONIC_GALBRUN_HXX

namespace Montjoie
{
  
  //! class to specify resolution of harmonic Galbrun equation for any flow
  /*!
    Galbrun's equations is given as :
    rho (-i omega + sigma + M \cdot \nabla)^2 u - \nabla(rho c^2 div u) + (div u) \nabla p_0 - (\nabla u)^T \nabla p_0 = f
    where (\nabla u)^T = | du_x/dx du_y/dx |
                         | du_x/dy du_y/dy |
                         
    rho, c, M, p_0 and sigma are given coefficients (stationary flow)

    Dirichlet condition : u.n = 0
    Neumann condition p = 0 (i.e. div u = 0)
   */
  template<class T, class Dim>
  class GalbrunEquationSipg_Base : public GenericEquation<T>
  {
  public :
    typedef Dim Dimension;
    
    enum {nb_unknowns = Dimension::dim_N, nb_components_u = 1,
          nb_components_grad = Dimension::dim_N, nb_unknowns_hdg=0,
	  nb_components_en = 1, nb_components_hn = 1,
	  nb_unknowns_scal = Dimension::dim_N, nb_unknowns_vec = 0};

    static inline bool SymmetricGlobalMatrix() { return false; }
    static inline bool SymmetricElementaryMatrix() { return false; }

    static void SetIndexToCompute(VarGalbrunIndex_Base<Dimension>& var);
    
    // for compatbility purpose
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int num_elem, const ElementReference_Dim<Dimension>& Fb);
    
    // For a detailed description of the following methods
    // look at the file GenericEquation.hxx (class GenericEquation_Base)
    template<class GenericPb, class T0, class Vector1>
    static void GetNeededDerivative(const GenericPb& vars,
                                    const GlobalGenericMatrix<T0>& nat_mat,
				    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var,
                                     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& dU, Vector1& dV);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                  int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                                  int ref, MatStiff& Cgrad_grad);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& Vn, Vector2& Un);
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& Cj);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& Un, Vector1& Vn);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
			int ref, int iquad, int npoint,
			const GlobalGenericMatrix<T0>& nat_mat, int ref2, 
			const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(typename Dimension::R_N& normale, int ref,
                        const Vector1& Vn, Vector1& Un, int num_elem1,
                        int npoint, const GlobalGenericMatrix<T0>& nat_mat, int ref2, 
                        const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& );
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int npoint,
			   int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class Vector2, class GenericPb, class T0>
    static void MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector2& Un,
			   int iquad, int npoint, int nf, const GlobalGenericMatrix<T0>& nat_mat,
                           int ref, int ref2, const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
  };
  
  
  //! stationary Galbrun's equation for any flow and SIPG
  template<class Dimension>
  class GalbrunStationaryEquationSipg : public GalbrunEquationSipg_Base<Real_wp, Dimension>
  {
  public :
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V);
    
  };
  
  
  //! time-harmonic Galbrun's equation for any flow and SIPG
  template<class Dimension>
  class HarmonicGalbrunEquationSipg : public GalbrunEquationSipg_Base<Complex_wp, Dimension>
  {
  public :
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V);
    
  };


  //! class for absorbing boundary condition for Galbrun's equation and SIPG formulation
  template<class T, class Dimension>
  class ImpedanceABC<T, HarmonicGalbrunEquationSipg<Dimension> >
    : public ImpedanceFunction_Base<T, Dimension>
  {
  private:
    const VarProblem<Dimension>& var_problem;
    const VarGalbrunIndex_Base<Dimension>& var_galbrun;

  public :
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);

    void EvaluateImpedancePhi_H1(int i, int iquad, int num_edge, int num_loc, int k,
                                 const GlobalGenericMatrix<T>& nat_mat, int ref_domain,
                                 const SetPoints<Dimension>&, const SetMatrices<Dimension>&);

    void ApplyImpedancePhi_H1(int m, int j, int, const TinyVector<Real_wp, 1>& phi,
                              const typename Dimension::R_N& grad_phi, Vector<T>& f_phi);
    
  };
  
  
  //! class used to solve time-harmonic Galbrun's equation with SIPG method
  template<class Dimension>
  class EllipticProblem<HarmonicGalbrunEquationSipg<Dimension> >
    : public VarGalbrun_Eq<HarmonicGalbrunEquationSipg<Dimension> >
  {
  public:
    inline EllipticProblem() { this->dg_formulation = ElementReference_Base::DISCONTINUOUS;
      this->alpha_penalization = -Real_wp(1); this->delta_penalization = -Real_wp(1);}
    
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
  
  
  //! class used to solve real Galbrun's equation with SIPG method
  template<class Dimension>
  class EllipticProblem<GalbrunStationaryEquationSipg<Dimension> >
    : public VarGalbrun_Eq<GalbrunStationaryEquationSipg<Dimension> >
  {
  public:
    inline EllipticProblem() { this->dg_formulation = ElementReference_Base::DISCONTINUOUS;
      this->alpha_penalization = -Real_wp(1); this->delta_penalization = -Real_wp(1);}
  };


  //! class to specify resolution of harmonic Galbrun equation for any flow
  /*!
    Galbrun's equations is given as :
    rho (-i omega + sigma + M \cdot \nabla)^2 u - \nabla(rho c^2 div u) + (div u) \nabla p_0 - (\nabla u)^T \nabla p_0 = f
    where (\nabla u)^T = | du_x/dx du_y/dx |
                         | du_x/dy du_y/dy |
                         
    rho, c, M, p_0 and sigma are given coefficients (stationary flow)

    Dirichlet condition : u.n = 0
    Neumann condition p = 0 (i.e. div u = 0)
   */
  template<class T>
  class GalbrunEquationH1_2D : public GenericEquation<T>
  {
  public :
    typedef Dimension2 Dimension;
    
    // for H1(DG), set nb_unknowns=4    
    enum {nb_unknowns = 3, nb_components_u = 1,
          nb_components_grad = 2, nb_unknowns_hdg=0,
	  nb_components_en = 1, nb_components_hn = 1,
	  nb_unknowns_scal = 3, nb_unknowns_vec = 0};

    static inline bool SymmetricGlobalMatrix() { return false; }
    static inline bool SymmetricElementaryMatrix() { return false; }

    static void SetIndexToCompute(VarGalbrunIndex_Base<Dimension2>& var);
    
    // for compatbility purpose
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int num_elem, const ElementReference_Dim<Dimension>& Fb);
        
    // For a detailed description of the following methods
    // look at the file GenericEquation.hxx (class GenericEquation_Base)
    template<class GenericPb, class T0, class Vector1>
    static void GetNeededDerivative(const GenericPb& vars,
                                    const GlobalGenericMatrix<T0>& nat_mat, Vector1& At, Vector1& A);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var,
                                     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& dU, Vector1& dV);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                  int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                                  int ref, MatStiff& Cgrad_grad);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& dU, Vector2& V);

    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
			int ref, int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
			const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(typename Dimension::R_N& normale, int ref,
                        const Vector1& Vn, Vector1& Un, int num_elem1,
                        int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
                        const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& );
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int k,
			   int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class Vector2, class GenericPb, class T0>
    static void MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector2& Un,
			   int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
                           int ref, int ref2, const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
  };
  
  
  //! stationary Galbrun's equation for any flow and SIPG
  class GalbrunStationaryEquationH1 : public GalbrunEquationH1_2D<Real_wp>
  {
  public :
  };
  
  
  //! time-harmonic Galbrun's equation for any flow and H1
  class HarmonicGalbrunEquationH1 : public GalbrunEquationH1_2D<Complex_wp>
  {
  public :
  };


  //! volumetric source
  template<class T>
  class VolumetricSource_GalbrunH1 : public VirtualSourceFEM<T, Dimension2>
  {
  protected:
    GaussianSource<Dimension2> fsrc;
    T coef_vol;
    
  private:
    const DistributedProblem<Dimension2>& var_problem;
    const VarGalbrunIndex_Base<Dimension2>& var_galbrun;
    const VarGalbrun_Base& var_galbrun_base;
    
  public :
    
    template<class TypeEquation>
    VolumetricSource_GalbrunH1(const EllipticProblem<TypeEquation>& var,
			       const Vector<VectString>& param);

    bool IsNonNullVolumetricSource(const VectR2& s);
    void EvaluateVolumetricSource(int i, int j, const R2& x, Vector<T>& f);
    
  };
  
  
  //! volumetric source
  template<>
  class VolumetricSource<HarmonicGalbrunEquationH1> : public VolumetricSource_GalbrunH1<Complex_wp>
  {
  public :
    template<class TypeEquation>
    VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param);
    
  };


  //! class for absorbing boundary condition for Galbrun's equation and H1 formulation
  template<class T>
  class ImpedanceABC<T, HarmonicGalbrunEquationH1>
    : public ImpedanceABC<T, HarmonicGalbrunEquationSipg<Dimension2> >
  {
  public:
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);
    
  };


  //! class used to solve time-harmonic Galbrun's equation with H^1 method
  template<>
  class EllipticProblem<HarmonicGalbrunEquationH1>
    : public VarGalbrun_Eq<HarmonicGalbrunEquationH1>
  {
  public:
    inline EllipticProblem()
    {
      if (this->nb_unknowns == 4)
	{
	  this->dg_formulation = ElementReference_Base::DISCONTINUOUS;
	  this->alpha_penalization = -Real_wp(1); this->delta_penalization = -Real_wp(1);
	}
    }

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
  
  
  //! class used to solve real Galbrun's equation with H^1 method
  template<>
  class EllipticProblem<GalbrunStationaryEquationH1 >
    : public VarGalbrun_Eq<GalbrunStationaryEquationH1>
  {
  public:
    inline EllipticProblem()
    {
      if (this->nb_unknowns == 4)
	{
	  this->dg_formulation = ElementReference_Base::DISCONTINUOUS;
	  this->alpha_penalization = -Real_wp(1); this->delta_penalization = -Real_wp(1);
	}
    }

  };
  
}

#define MONTJOIE_FILE_SYMMETRIC_HARMONIC_GALBRUN_HXX
#endif
