#ifndef MONTJOIE_FILE_USER_SOURCE_HXX

namespace Montjoie
{
  
  //! class to define source, that can be changed by the user
  /*!
    definition of methods is in file UserSource.cxx
   */
  template<class T, class Dimension>
  class UserDefinedSource : public VirtualSourceFEM<T, Dimension>
  {
  public :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    
    const Vector<string>& source_space_param;
   
    template<class TypeEquation>
    UserDefinedSource(const EllipticProblem<TypeEquation>& var,
		      const Vector<VectString>& param);

    void EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f);
    
    bool IsNonNullVolumetricSource(const VectR_N& s);
    void EvaluateVolumetricSource(int i, int j, const R_N& x, Vector<T>& f);
    
    bool IsNonNullGradientSource(const VectR_N& s);
    void EvaluateGradientSource(int i, int j, const R_N& x, Vector<T>& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				const SetMatrices<Dimension>& MatricesElem, Vector<T>& f);
    
    bool IsNonNullSurfacicSourceGradient(int cond_ref);
    void EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension>& PointsElem,
					const SetMatrices<Dimension>& MatricesElem, Vector<T>& f);
    
  };

  
  //! class that enables the user to change initial conditions for fields
  template<class T, class Dimension>
  class InitialUserFunction : public VirtualProjectorFEM<T, Dimension>
  {
  public :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    
    const Vector<string>& param_initial_condition;
    
    template<class TypeEquation>
    InitialUserFunction(const EllipticProblem<TypeEquation>&,
			const Vector<string>& param);
    
    void InitElement(int i, const typename Dimension::VectR_N& s);
    void EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f);
    
  };
  

  //! Base class used to specify a time source (function that depends only on t)
  /*!
    The user should not modify this class
  */
  template<class T>
  class VirtualTimeSource
  {
  public:
    virtual ~VirtualTimeSource();
    
    // returns the function f (should always be overloaded)
    virtual T Evaluate(const Real_wp&) = 0;
    // derivative is optional
    virtual T EvaluateDerivative(const Real_wp&);
    
  };
  
  
  //! source in time, which can be modified by the user in file UserSource.cxx
  class TimeUserSource : public VirtualTimeSource<Real_wp>
  {
  public :
    TimeUserSource(const Real_wp& freq, Real_wp& t_final, const VectString& param);
    
    Real_wp Evaluate(const Real_wp&);
    
  };

  
  // definition of variable media (rho, mu, epsilon, etc)
  template<class Dimension, class T>
  void ComputeVariableUserIndex(const VarGeometryProblem<Dimension>& var,
                                const IVect& ElementRho, int ref, int num_index, int num_component,
                                const Vector<typename Dimension::VectR_N>& PointsQuadrature,
				bool compute_grad, bool compute_hess,
                                Vector<Vector<T> >& rho,
				Vector<Vector<TinyVector<T, Dimension::dim_N> > >& grad_rho,
				Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N,
				Dimension::dim_N> > >& hess_rho,
				const T& offset, const T& amplitude);
  
} // namespace Montjoie

#define MONTJOIE_FILE_USER_SOURCE_HXX
#endif
