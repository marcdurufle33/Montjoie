#ifndef MONTJOIE_FILE_TIME_REVERSAL_HXX

namespace Montjoie
{

  template<class Dimension>
  class TimeAcousticReversal;
  
  //! class used to handle time-reversal method
  /*!
    class used for direct simulation (storing values on a surface)
    and for reverse simulation (source on a surface)
   */
  template<class Dimension>
  class TimeReversalSource : public VirtualSourceFEM<Real_wp, Dimension>
  {
  protected :    
    //! order of interpolation
    int order;
    //! interpolation data
    Globatto<Real_wp> lob;
    //! history of v dot n at different times
    Vector<VectReal_wp> ValueVn;
    VectReal_wp ValueTime; //!< time t for each value v dot n
    VectReal_wp InterpolVn; //!< value v dot n (related to ValueTime)
    //! time to reach in order to write the next copy
    Real_wp t_current;
    Real_wp t0_buffer; //!< buffer [t0,t1] 
    Real_wp t1_buffer; //!< buffer [t0,t1] 
    int n0_buffer; //!< iteration number related to t0_buffer
    int n1_buffer; //!< iteration number related to t1_buffer
    
  private:
    TimeAcousticReversal<Dimension>& time_reversal;

  public:
    //! constructor with given problem
    template<class TypeEquation>
    TimeReversalSource(const EllipticProblem<TypeEquation>& var,
		       TimeAcousticReversal<Dimension>& time_rev);

    void SetCurrentTime(const Real_wp& t);
    const Real_wp& GetCurrentTime() const;
    
    void InitInterpolation();
    
    void ComputeSource(const Real_wp& t, const Real_wp&, int n);
    void InitSourceSurface(const Real_wp& t, const Real_wp&);

    bool IsNonNullSurfacicSource(int ref);
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				const SetMatrices<Dimension>& MatricesElem, VectReal_wp& f);
    
  };

  
  //! Base class for time-reversal with acoustic
  template<class Dimension>
  class TimeAcousticReversal
  {
    template<class Dim> friend class TimeReversalSource;
    
  protected: 
    //! time-reversal experiment (direct or reverse)
    bool presence_time_reversal;
    //! direct or inverse simulation ?
    int type_simulation;
    //! time interval [t0,t1]
    Real_wp tbegin, tend;
    //! deltat, time step between each copy on the disk
    Real_wp deltat;
    //! object used to compute finite element right hand side
    TimeReversalSource<Dimension> src_reverse;

    IVect OffsetQuad; //!< offset for quadrature points of the surface Gamma
    Vector<int> IndexSurface;
    
    Vector<int> is_ref_on_mode;    
    MeshInterpolationFEM<Dimension> mesh_interp;
    string file_mesh, file_out_dtn;
    
  private:
    VarProblem<Dimension>& var_problem;
    VarSourceProblem_Dim<Dimension>& var_source;
    VarInstationary_Dim<Dimension>& var_time;

  public:
    typedef float value_type;
    enum {DIRECT, INVERSE}; //!< direct or inverse
    
    template<class TypeEquation>
    TimeAcousticReversal(HyperbolicProblem<TypeEquation>& var);

    // Inline methods
    bool IsPresent() const;
    int GetSimulationType() const;
    Real_wp GetInitialTime() const;
    Real_wp GetFinalTime() const;
    
    // other methods
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void InitComputation();
    void WriteOutput(int nb_iter, const Real_wp& t, const VectReal_wp Vh);
   
    void GetVn(const VectReal_wp& Vn, Vector<value_type>& du) const;
    void ApplyFilterVn(VectReal_wp& dU_dn, const Real_wp& t);
    
    // computation of source
    bool SourceDoesNotDependOnTime() const;
    void AddScalarSourceAtTime(const Real_wp&, const Real_wp&, int nb_deriv, VectReal_wp& b_src);
    
  };

}

#define MONTJOIE_FILE_TIME_REVERSAL_HXX
#endif
