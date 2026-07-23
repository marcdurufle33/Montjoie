#ifndef MONTJOIE_FILE_VLASOV_MAXWELL_3D_HXX

namespace Montjoie
{
  
  //! class to solve Maxwell-Vlasov equations in 3-D
  class VlasovMaxwell<Dimension3> : public VlasovMaxwell_Base<Dimension3>
  {
  public :
    //! true if a uniform magnetic field is set
    bool uniform_magnetic_field;
    //! value of the uniform magnetic field if present
    R3 value_uniform_magnetic;
    
    VectReal_wp WeightsBoundary;

  public:    
    VlasovMaxwell_3D();
    
    //! reading of input file
    void SetInputData(const string& description_field, const VectString& parameters);
    int GetOffset_PosParticle() const;
    
    // computation of the localization grid for particles
    void GeneratePeriodicDirections();
    void PreComputeCloudsGrid(VectR3& pts, const Real_wp& radius,
                              Vector<IVect>& Cloud_elem, GridInterpolation<Dimension3>& grid);
    
    // time iterations
    void InitTimeIterations();
    void ComputeProjectors(VectR3& PtsInterpol3D, VectR3& PtsQuad3D, VectReal_wp& Weights3D);
    void RunTimeIterations();

    //! definition of the operator A for the ode dU/dt = A U
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                                    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert_mass = true, bool source = true);
    
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                                    const Real_wp& alpha, int level,
                                    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert_mass = true, bool source = false);
    
    Real_wp GetNormeSolution(const VectReal_wp& Uh) const;
    void SetUniform_MagneticField(const R3& valH, VectReal_wp& Y, GhostIf<false>& dg_form);
    void SetUniform_MagneticField(const R3& valH, VectReal_wp& Y, GhostIf<true>& dg_form);
    
    // creation of particles
    // preload
    void SetInitialVector(const Real_wp& t_begin, VectReal_wp& Y);
    // or for beams and space charge limited emission
    void GiveIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y);
    void AddParticles(int nb_part, Real_wp& q0, Real_wp& m0,
		      VectR3& p0, VectR3& pts, VectReal_wp& poids, int offset, VectReal_wp& Y);
    
    // localization of particles
    void LocalizeParticles(const VectReal_wp& X, IVect& NumBoxGrid_Point);
    void MovePoints_Periodicity(VectR3& Points);
    void CreateParticlesBreakdownField(const VectReal_wp& En, const VectReal_wp& Hn, int nb_iter,
				       VectReal_wp& Y, int offset,
                                       VectReal_wp& WeightParticle, VectR3& CoorParticle);
    
    // methods used for boris correction
    void EvaluateRho(const VectReal_wp& Y, VectReal_wp& bsrc,
                     CurrentSource_Vlasov<TypeElement,TypeEquation>& Jsrc_rho);
    void EvaluateSourceBoris(VectReal_wp& b_rhs, VectReal_wp& evalRho);
    void EvaluateSourceBoris(VectReal_wp& b_rhs, VectReal_wp& evalRho,
                             VectR3& evalEn, VectR3& evalEnBoundary);
    void ComputeValue_PhiNodal(VectReal_wp& phi, const VectReal_wp& evalRho);
    void AddSurfacicIntegral_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                                      VectReal_wp& Prod_En, GhostIf<false>&);
    void AddSurfacicIntegral_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                                      VectReal_wp& Prod_En, GhostIf<true>&);
    void AddVolumetricIntegral_PhiNodal(const Real_wp&, const VectReal_wp& phi,
                                        VectReal_wp& evalE, GhostIf<false>&);
    void AddVolumetricIntegral_PhiNodal(const Real_wp&, const VectReal_wp& phi,
                                        VectReal_wp& evalE, GhostIf<true>&);
    void AddVolumetricProjection_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                                          VectReal_wp& Prod_En, GhostIf<true>&);
    void AddVolumetricProjection_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                                          VectReal_wp& Prod_En, GhostIf<false>&);

    // treatment of the coupling matrix Ch
    void EvaluateE_H_Quadrature(const VectReal_wp& En, const VectReal_wp& Hn,
				VectR3& Equad, VectR3& Hquad, bool, GhostIf<true>& dg_form);
    void EvaluateE_H_Quadrature(const VectReal_wp& En, const VectReal_wp& Hn,
				VectR3& Equad, VectR3& Hquad, bool, GhostIf<false>& dg_form);
    void EvaluateE_H_Interpolate(const VectR3& En, const VectR3& Hn,
				 VectR3& Equad, VectR3& Hquad);
    void AddIntegralCurrentJ(const Real_wp& alpha, const VectR3& Jn_quad,
			     VectReal_wp& Prod_En, GhostIf<true>& dg_form);
    void AddIntegralCurrentJ(const Real_wp& alpha, const VectR3& Jn_quad,
			     VectReal_wp& Prod_En, GhostIf<false>& dg_form);
    void EvaluateCurrentQuadrature(const VectR3& Jn_interp, VectR3& Jn_quad);
    

    void EvaluateE_H_QuadratureBoundary(const VectReal_wp& En, const VectReal_wp& Hn,
					VectR3& Equad, VectR3& Hquad, GhostIf<true>& dg_form);
    void EvaluateE_H_QuadratureBoundary(const VectReal_wp& En, const VectReal_wp& Hn,
					VectR3& Equad, VectR3& Hquad, GhostIf<false>& dg_form);
    
    //! writing snapshots of the solution
    void WriteSnapshot(int nb_iter, const Real_wp& t, const VectReal_wp& Uh);
    
  };
  

  //! empty class
  class VlasovMaxwellEquation3D : public TimeMaxwellEquation_3D
  {
  public :
    
  };


  //! empty class
  class VlasovMaxwellEquation_3D_DG : public TimeMaxwellEquation_3D_DG
  {
  public :
    
  };
  
  
  //! class to solve Vlasov-Maxwell equations with edge finite elements
  template<>
  class HyperbolicProblem<VlasovMaxwellEquation3D>
    : public VlasovMaxwell<Dimension3>, public TimeMaxwell_3D<VlasovMaxwellEquation3D>
  {
  public :
    //! default constructor
    HyperbolicProblem();
    
  };
  
  
  //! class to solve Vlasov-Maxwell equations with DG
  template<>
  class HyperbolicProblem<VlasovMaxwellEquation_3D_DG>
    : public VlasovMaxwell<Dimension3>, public TimeMaxwell_3D<VlasovMaxwellEquation_3D_DG>
  {
  public :
    //! default constructor
    HyperbolicProblem();
    
  };
  
}

#define MONTJOIE_FILE_VLASOV_MAXWELL_3D_HXX
#endif
