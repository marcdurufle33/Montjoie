#ifndef MONTJOIE_FILE_STATIC_MAXWELL_AXI_HXX

namespace Montjoie
{
  
  //! parameters for axisymmetric Maxwell equations
  class StaticMaxwellEquation_HcurlAxi : public GenericEquation<Real_wp>
  {
  public :
    typedef Dimension2 Dimension; 

    enum{nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_hdg=0,
	 nb_components_en = 2, nb_components_hn = 2, nb_unknowns_vec = 3,
    type_element = 2};

    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();

    static inline IVect GetOtherElementType() { IVect num(1); num(0) = 1; return num; }

    static Real_wp GetCoefficientMassMatrix(const VarProblem_Base& vars);
    static Real_wp GetCoefficientStiffnessMatrix(const VarProblem_Base& vars);
    static Real_wp GetCoefMode(const VarProblem_Base& vars);
    
  };
  
  
  //! class to solve static Maxwell equations in axisymmetric domains 
  //! with first family on quadrilaterals
  template<>
  class EllipticProblem<StaticMaxwellEquation_HcurlAxi>
    : public HarmonicMaxwellAxi_Hcurl<StaticMaxwellEquation_HcurlAxi>
  {
  public :
    void GetDof_Eteta(VectBool& IsDofTeta) const;
    
  };


  //! volumetric source
  template<>
  class VolumetricSource<StaticMaxwellEquation_HcurlAxi>
    : public VolumetricSource_MaxwellAxi<Real_wp>
  {
  public :
    template<class TypeEquation>
    VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param);
    
  };
  
  
  //! base class to solve static Maxwell equation on axisymmetric domain with DG methods
  template<class TypeEquation>
  class StaticMaxwell_Axi_DG : public HarmonicMaxwellAxi_Eq<TypeEquation>
  {
  public :
    typedef typename TypeEquation::Complexe Complexe; //!< real or complex
    
    void GetDof_Eteta(VectBool& IsDofTeta) const;
    
    void EvaluateH_MixedFormulation(const Vector<Complexe>& U0,
                                    Vector<Complexe>& EvalH_Nodal) const;
    
  };
  
  
  //! class to specify the resolution of stationary Maxwell equations with axisymmetric geometry
  class StaticMaxwellEquation_Axi_DG
    : public GenericEquation<Real_wp>
  {
  public :
    typedef Dimension2 Dimension; 

    static const bool FirstOrderFormulation = true;
    
    enum {nb_unknowns = 6, nb_unknowns_hdg=0,
	  nb_components_en = 2, nb_components_hn = 2, nb_unknowns_scal = 3, nb_unknowns_vec = 3};

    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();
    
    template<class TypeEquation, class T0, class Vector1>
    static void GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
                                    const GlobalGenericMatrix<T0>& nat_mat,
                                    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);

    template<class Vector1, class Matrix1, class T0, class TypeEquation>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& vars,
                                     int i,int  k, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, const Matrix1& Vn, Vector1& Un);

    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& vars,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& Un, Vector2& Vn);

    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
				 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
				 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);

    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, R2& normale, int ref, int iquad, int k,
			const GlobalGenericMatrix<T0>& nat_mat, int ref2,
                        const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb);

    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(R2& normale, int ref, const Vector1& Vn, Vector1& Un, int num_elem1, int k,
                        const GlobalGenericMatrix<T0>& nat_mat, int ref2, 
                        const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension2, 1>& Fb);

    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& Un, Vector1& Vn);

    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& Cj);
    
  };


  //! class used to solve static Maxwell equation with DG Gauss hexahedra
  template<>
  class EllipticProblem<StaticMaxwellEquation_Axi_DG>
    : public StaticMaxwell_Axi_DG<StaticMaxwellEquation_Axi_DG>
  {    
  };
    

  //! volumetric source
  /* template<>
  class VolumetricSource<StaticMaxwellEquation_Axi_DG>
    : public VolumetricSource_MaxwellAxi<Real_wp>
  {
    
  public :
    //! constructor with given problem
    VolumetricSource(EllipticProblem<StaticMaxwellEquation_Axi_DG>& var) 
      : VolumetricSource_MaxwellAxi<Real_wp>(var)
    { }
    
    };*/
  
}

#define MONTJOIE_FILE_STATIC_MAXWELL_AXI_HXX
#endif

