#ifndef MONTJOIE_FILE_VLASOV_MAXWELL_2D_HXX

namespace Montjoie
{
  
  //! class to solve Maxwell-Vlasov equations in 2-D
  class VlasovMaxwell<Dimension2> : public VlasovMaxwell_Base<Dimension2>
  {
  public :
    //! true if a uniform magnetic field is set
    bool uniform_magnetic_field;
    //! value of the uniform magnetic field if present
    Real_wp value_uniform_magnetic;
    
    bool mirror_right_side;
    VectReal_wp WeightsBoundary;
    Globatto<Real_wp> lob_boundary;

  private:
    
  public:
    template<class TypeEquation>
    VlasovMaxwell(HyperbolicProblem<TypeEquation>&);
    
    //! reading of input file
    void SetInputData(const string& description_field, const VectString& parameters);
    int GetOffset_PosParticle() const;
    
    // computation of the localization grid for particles
    void GeneratePeriodicDirections();
    void PreComputeCloudsGrid(VectR2& pts, const Real_wp& radius,
                              Vector<IVect>& Cloud_elem, GridInterpolation<Dimension2>& grid);
    
    // time iterations
    void ComputeProjectors(VectR2& PtsInterpol2D, VectR2& PtsQuad2D, VectReal_wp& Weights2D);
    void RunTimeIterations();

    //! definition of the operator A for the ode dU/dt = A U
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                                    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert_mass = true, bool source = true);
    
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp&, int level,
                                    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert_mass = true, bool source = false);
    
    Real_wp GetNormeSolution(const VectReal_wp& Uh) const;
    void SetUniform_MagneticField(const Real_wp& valH, VectReal_wp& Y, GhostIf<false>& dg_form);
    void SetUniform_MagneticField(const Real_wp& valH, VectReal_wp& Y, GhostIf<true>& dg_form);
    
    // creation of particles
    // preload
    void SetInitialVector(const Real_wp& t_begin, VectReal_wp& Y);
    // or for beams and space charge limited emission
    void GiveIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y);
    void AddParticles(int nb_part, Real_wp& q0, Real_wp& m0,
		      VectR2& p0, VectR2& pts, VectReal_wp& poids, int offset, VectReal_wp& Y);
    void ComputeWeightsBeam(Vlasov_BeamCurrent<Dimension2>& current,
                            VectR2& NewCoor, Real_wp& poids, VectReal_wp& NewWeight,
			    CurrentSource_Vlasov<TypeElement,TypeEquation>& Jsrc);
    void ComputePoints_ValEqualConstant(const Globatto<Real_wp>& lob, const VectReal_wp& value,
					const Real_wp& cte, VectReal_wp& pts);
    
    // localization of particles
    void LocalizeParticles(const VectReal_wp& X, IVect& NumBoxGrid_Point);
    void MovePoints_Periodicity(VectR2& Points);
    void CreateParticlesBreakdownField(const VectReal_wp& En, const VectReal_wp& Hn, int nb_iter,
				       VectReal_wp& Y, int offset,
                                       VectReal_wp& WeightParticle, VectR2& CoorParticle);
    
    // methods used for boris correction
    void EvaluateRho(const VectReal_wp& Y, VectReal_wp& bsrc,
                     CurrentSource_Vlasov<TypeElement,TypeEquation>& Jsrc_rho);
    void EvaluateSourceBoris(VectReal_wp& b_rhs, VectReal_wp& evalRho);
    void EvaluateSourceBoris(VectReal_wp& b_rhs, VectReal_wp& evalRho,
                             VectR2& evalEn, VectR2& evalEnBoundary);
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
				VectR2& Equad, VectReal_wp& Hquad, bool, GhostIf<true>& dg_form);
    void EvaluateE_H_Quadrature(const VectReal_wp& En, const VectReal_wp& Hn,
				VectR2& Equad, VectReal_wp& Hquad, bool, GhostIf<false>& dg_form);
    void EvaluateE_H_Interpolate(const VectR2& En, const VectReal_wp& Hn,
				 VectR2& Equad, VectReal_wp& Hquad);
    void AddIntegralCurrentJ(const Real_wp& alpha, const VectR2& Jn_quad,
			     VectReal_wp& Prod_En, GhostIf<true>& dg_form);
    void AddIntegralCurrentJ(const Real_wp& alpha, const VectR2& Jn_quad,
			     VectReal_wp& Prod_En, GhostIf<false>& dg_form);
    void EvaluateCurrentQuadrature(const VectR2& Jn_interp, VectR2& Jn_quad);
    

    void EvaluateE_H_QuadratureBoundary(const VectReal_wp& En, const VectReal_wp& Hn,
					VectR2& Equad, VectReal_wp& Hquad, GhostIf<true>& dg_form);
    void EvaluateE_H_QuadratureBoundary(const VectReal_wp& En, const VectReal_wp& Hn,
					VectR2& Equad, VectReal_wp& Hquad,
                                        GhostIf<false>& dg_form);
    
    //! writing snapshots of the solution
    void WriteSnapshot(int nb_iter, const Real_wp& t, const VectReal_wp& Uh);
    
  };
  
  //! empty class
  class VlasovMaxwellEquation2D : public TimeMaxwellEquation_2D
  {
  public :
    
  };


  //! empty class
  class VlasovMaxwellEquation_2D_DG : public TimeMaxwellEquation_2D_DG
  {
  public :
    
  };
  
  
  //! class to solve Vlasov-Maxwell equations with edge finite elements
  template<>
  class HyperbolicProblem<VlasovMaxwellEquation2D>
    : public VlasovMaxwell<Dimension2>, public TimeMaxwell_2D<VlasovMaxwellEquation2D>
  {    
  public:
    HyperbolicProblem();
    
  };
  
  
  //! class to solve Vlasov-Maxwell equations with DG
  template<>
  class HyperbolicProblem<VlasovMaxwellEquation_2D_DG>
    : public VlasovMaxwell<Dimension2>, public TimeMaxwell_2D<VlasovMaxwellEquation_2D_DG>
  {
  public :
    HyperbolicProblem();
    
  };
  
}

#define MONTJOIE_FILE_VLASOV_MAXWELL_2D_HXX
#endif
