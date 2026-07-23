#ifndef MONTJOIE_FILE_MULTI_DIM_PHYSICAL_INDEX_HXX

namespace Montjoie
{

  //! base class for multi-dimensional parameters
  template<class T, class Dimension>
  class VariableParameterND_Base
  {
  protected:
    typedef typename Dimension::R_N R_N;
    
    // number of degrees of freedom for the parameter
    int nb_points_grid_c;
    // position of points used to project the function
    Vector<R_N> all_points;
    // order of approximation
    int order;
    
  public:
    
    VariableParameterND_Base();

    virtual ~VariableParameterND_Base();

    int GetNbDof() const;
    const typename Dimension::VectR_N& GetInterpolationPoints() const;
    
    static VariableParameterND_Base<T, Dimension>* GetNewApproximation(DistributedProblem<Dimension>& var,
								       const Vector<string>& param);
    virtual T Evaluate(const R_N& pt_glob, int ne, const R_N& pt_loc) = 0;
    
    virtual void ProjectFunction(T (&f)(const R_N&), Vector<T>& projF);

    virtual void ComputeValuesPhiRef(int i, int j, const R_N& point, 
				     const R_N& point_loc, Vector<Real_wp, VectSparse>& phi) = 0;
    
    virtual void Init(const Vector<T>&) = 0;
    
    virtual void GetWeightFunction(VectReal_wp&);
    virtual void KeepElements(const Vector<bool>& num);
    
  };
  
  
  //! multi-dimensional parameter based on the same mesh as the computation of the solution
  template<class T, class Dimension>
  class SameMeshIndex : public VariableParameterND_Base<T, Dimension>
  {
  protected:
    typedef typename Dimension::R_N R_N;
    Vector<T> val_cn;
    Mesh<Dimension> mesh;
    Vector<int> ElemIndex;
    bool discontinuous;
    Vector<IVect> Nodle;
    IVect offset;
    
  public:
    SameMeshIndex();

    void Init(const Vector<T>& cn);    
    void SetApproximation(int order, bool discont, DistributedProblem<Dimension>& var);

    T Evaluate(const R_N& pt_glob, int ne, const R_N& pt_loc);

    void ComputeValuesPhiRef(int i, int jloc, const R_N& point, 
			     const R_N& point_loc, Vector<Real_wp, VectSparse>& phi);
    
    IVect GetDofsOnBoundary(int ne, int num_loc);
    
    void GetWeightFunction(VectReal_wp&);
    void KeepElements(const Vector<bool>& num);

  };


  //! constant multi-dimensional parameter 
  template<class T, class Dimension>
  class UniformParameterIndex : public VariableParameterND_Base<T, Dimension>
  {
  protected:
    typedef typename Dimension::R_N R_N;
    T val_cn;
    
  public:
    UniformParameterIndex();

    void Init(const Vector<T>& cn);    

    T Evaluate(const R_N& pt_glob, int ne, const R_N& pt_loc);

    void ComputeValuesPhiRef(int i, int jloc, const R_N& point, 
			     const R_N& point_loc, Vector<Real_wp, VectSparse>& phi);
    
    void GetWeightFunction(VectReal_wp&);
    
  };
  
  
  //! multi-dimensional parameter defined as tensor product of B-splines
  template<class T, class Dimension>
  class RectangularBSplineIndex : public VariableParameterND_Base<T, Dimension>
  {
  protected:
    typedef typename Dimension::R_N R_N;
    BSplineInterpolation<T> spline_x, spline_y;
    mutable Vector<Real_wp> value_phi_x, value_phi_y;
    Vector<T> val_cn;
    
  public:
    void Init(const Vector<T>& cn);

    void SetApproximation(const VectReal_wp& points_x, const VectReal_wp& points_y, int order_x, int order_y);
    
    T Evaluate(const R_N& pt_glob, int ne, const R_N& pt_loc);
    void ProjectFunction(T (&f)(const R_N&), Vector<T>& projF);
    
    void ComputeValuesPhiRef(int i, int j, const R_N& point, 
			     const R_N& point_loc, Vector<Real_wp, VectSparse>& phi);
    
  };
  
}

#define MONTJOIE_FILE_MULTI_DIM_PHYSICAL_INDEX_HXX
#endif

