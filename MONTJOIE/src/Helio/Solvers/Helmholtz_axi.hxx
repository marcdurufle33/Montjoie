#ifndef MONTJOIE_FILE_HELIO_AXISYM_HELMHOLTZ_HXX

namespace Montjoie
{

  //! Scalar Helmholtz equation for helioseismology with specific diffusion term
  //! -(\omega^2\rho + i\omega\sigma)u -2i\omega v\cdot\nabla u -2\dfrac{m}{r}\omega v_\theta u -\beta\nabla\cdot(\mu\nabla(\alpha u)) =s
  class HELIO_HelmholtzAxisymEquation : public GenericEquation<Complex_wp, 1>
  {
    public :

      typedef Dimension2          Dimension;

      enum {nb_unknowns      = 1,
            nb_unknowns_scal = 1,
            nb_components_en = 1,
            nb_components_hn = 1,
            nb_unknowns_vec  = 3};

    static inline bool SymmetricGlobalMatrix() { return false; }
    static inline bool SymmetricElementaryMatrix() { return false; }
    
  };


  //! Class to solve Helioseismology Helmholtz equation on axisymmetric domains
  template<class TypeEquation>
  class HELIO_VarHelmholtzAxisym : public VarHelmholtzAxi_Eq<TypeEquation>
  {
    
  public :
    typedef typename TypeEquation::Complexe     Complexe;
    
    // Additional physical coefficients (rho, sigma and mu are already present in VarHelmholtzAxi
    // =====================
    Vector<ScalarPhysicalIndice<Dimension2, Complexe> > ref_alpha, ref_beta; //!< Physical medium properties in term \beta \nabla\cdot(\mu\nabla(\alpha f))
    Complexe alpha0,beta0;                                                //!< Values of alpha and beta at infinity
    
    Real_wp  R_sun; //!< R_sun is needed to normalize the results of \grad\alpha and \grad\beta
    
    HELIO_VarHelmholtzAxisym();
    
    // Overloaded methods
    // ==================
    
    // Additionnal input data
    // ======================
    void SetInputData(const string& description_field, const VectString& parameters);
    
    // Physical properties
    // ===================
    
    // Allocation of indices
    void InitIndices(int n);
    
    // Sets physical indices of media with reference i
    void SetIndices(int i, const VectString& parameters);
    
    // Sets physical indices named in the data file with PhysicalMedia=...
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    // Filling arrays with physical properties
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
                           IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);

    // Sets values at infinity with the maximum of the coefficients
    void SetPhysicalIndexAtInfinity(const Vector<bool>&);


    void ComputeElementaryMatrix(int ID_elt, IVect& num_dof, VirtualMatrix<Complexe>& mat,
				 const GlobalGenericMatrix<Complexe>& nat_mat);
    
  };


  //! Class for absorbing boundary condition for Helioseismology Helmholtz equation in axisymmetric domain
  template<>
  class ImpedanceABC<Complex_wp, HELIO_HelmholtzAxisymEquation>
    : public ImpedanceFunction_Base<Complex_wp, Dimension2, 1>
  {
  private:
    const EllipticProblem<HELIO_HelmholtzAxisymEquation>& var_helm;
  
  public :
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>&);
    
    // Evaluation of impedance coefficient
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
			      const SetPoints<Dimension2>& Pts,
			      const SetMatrices<Dimension2>& Mat);
    
  };


  //! Class to solve Helioseismology Helmholtz equation on axisymmetric domains
  template<>
  class EllipticProblem<HELIO_HelmholtzAxisymEquation>
    : public HELIO_VarHelmholtzAxisym<HELIO_HelmholtzAxisymEquation>
  {
  };
  
}

#define MONTJOIE_FILE_HELIO_AXISYM_HELMHOLTZ_HXX
#endif
