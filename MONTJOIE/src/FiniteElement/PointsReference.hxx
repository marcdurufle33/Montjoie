#ifndef MONTJOIE_FILE_POINTS_REFERENCE_HXX

namespace Montjoie
{
  //! base class for "reference points", where dof, nodal and quadrature points are different
  template<class Dimension>
  class SetPoints
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    typedef typename Dimension::VectR_N VectR_N;  //!< vector of points
    //! nodal points after transformation F_i
    VectR_N PointsNodal;
    VectR_N PointsBoundary; //!< quadrature points on boundary after transformation F_i
    VectR_N PointsQuadrature; //!< quadrature points after transformation F_i
    VectR_N PointsDof; //!< dof points after transformation F_i
    VectR_N PointsDofBoundary; //!< dof points on the surface
    
  public :
    SetPoints();
    
    // Inline methods
    void ReallocatePointsDof(int N);
    void ReallocatePointsQuadrature(int N);
    void ReallocatePointsNodal(int N);
    void ReallocatePointsQuadratureBoundary(int N);
    void ReallocatePointsDofBoundary(int N);
    
    void SetPointDof(int i, const R_N& mat);
    void SetPointQuadrature(int i, const R_N& mat);
    void SetPointNodal(int i, const R_N& mat);
    void SetPointQuadratureBoundary(int k, const R_N& mat);
    void SetPointDofBoundary(int k, const R_N& mat);
    
    int GetNbPointsNodal() const;
    int GetNbPointsQuadrature() const;
    int GetNbPointsDof() const;
    int GetNbPointsQuadratureBoundary() const;
    int GetNbPointsDofBoundary() const;
    
    const R_N& GetPointDof(int i) const;
    const R_N& GetPointQuadrature(int i) const;
    const R_N& GetPointNodal(int i) const;
    const Vector<R_N>& GetPointNodal() const;
    const R_N& GetPointQuadratureBoundary(int k) const;
    const R_N& GetPointDofBoundary(int k) const;
    
    void CopyNodalToQuadrature();
    void CopyNodalToDof();
    void CopyQuadratureToDof();
    
    // other methods
    void RotatePoints(const Real_wp& teta);
    void TranslatePoints(const R_N& vec_u);
    
    template<class Dim>
    friend ostream& operator <<(ostream& out,const SetPoints<Dim>& e);
    
  };

  template<class Dim>
  ostream& operator <<(ostream& out,const SetPoints<Dim>& e);
  
  //! base class for jacobian matrices, where dof, nodal and quadrature points are different
  template<class Dimension>
  class SetMatrices
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector or points
    typedef typename Dimension::MatrixN_N MatrixN_N; //!< 2x2 or 3x3 matrix
    
    VectReal_wp DsQuadrature; //!< surfacic element in integration formula
    VectReal_wp K1_curve; //!< curvature K1
    VectReal_wp K2_curve; //!< curvature K2
    VectR_N NormaleQuadrature; //!< normale on each quadrature point
    //! jacobian matrices on nodal points DF_i
    Vector<MatrixN_N> MatricesNodal;
    //! jacobian matrices on boundary quadrature points DF_i
    Vector<MatrixN_N> MatricesBoundary; 
    Vector<MatrixN_N> MatricesNodalBoundary, invMatricesNodalBoundary;
    Vector<MatrixN_N> MatricesQuadrature; //!< jacobian matrices on quadrature points DF_i
    Vector<MatrixN_N> MatricesDof; //!< jacobian matrices on dof points DF_i
    Vector<MatrixN_N> MatricesDofBoundary; //!< jacobian matrices on boundary dof points
    VectReal_wp Weights_ND; //!< weight of integration
    SetPoints<Dimension>* points_elem;
    
  public :
    
    SetMatrices();

    // Inline methods
    void ReallocatePointsDof(int N);
    void ReallocatePointsQuadrature(int N);
    void ReallocatePointsNodal(int N);
    void ReallocatePointsDofBoundary(int N);
    
    int GetNbPointsQuadratureBoundary() const;
    int GetNbPointsDofBoundary() const;
    
    void SetPointDof(int i,const MatrixN_N& mat);
    void SetPointQuadrature(int i,const MatrixN_N& mat);
    void SetPointNodal(int i,const MatrixN_N& mat);
    void SetPointQuadratureBoundary(int k,const MatrixN_N& mat);
    void SetPointDofBoundary(int k,const MatrixN_N& mat);
    void SetPointNodalBoundary(const Vector<MatrixN_N>& mat);
    void SetInversePointNodalBoundary(const Vector<MatrixN_N>& mat);

    void SetNormaleQuadratureBoundary(int k,const R_N & normale_fj);
    void SetDsQuadratureBoundary(int k,const Real_wp& dsj);
    void SetK1QuadratureBoundary(int k,const Real_wp& h_courbure);
    void SetK2QuadratureBoundary(int k,const Real_wp& h_courbure);

    void FillNodal(const MatrixN_N& mat);
    void FillQuadrature(const MatrixN_N& mat);
    void FillDof(const MatrixN_N& mat);

    const MatrixN_N& GetPointQuadrature(int i) const;
    const MatrixN_N& GetPointDof(int i) const;
    const MatrixN_N& GetPointNodal(int i) const;
    const MatrixN_N& GetPointQuadratureBoundary(int k) const;
    const MatrixN_N& GetPointDofBoundary(int k) const;
    const MatrixN_N& GetPointNodalBoundary(int i) const;
    const MatrixN_N& GetInversePointNodalBoundary(int i) const;
    
    const R_N& GetNormaleQuadratureBoundary(int k) const ;
    const Real_wp& GetDsQuadratureBoundary(int k) const;
    const Real_wp& GetK1QuadratureBoundary(int k) const;
    const Real_wp& GetK2QuadratureBoundary(int k) const;
    
    void CopyNodalToQuadrature();
    void CopyNodalToDof();
    void CopyQuadratureToDof();
    
    // other methods
    void ReallocatePointsQuadratureBoundary(int N);

    void ComputeNodalBoundary(int num_loc, const ElementReference_Dim<Dimension>& Fb);
    
    void RotateNormale(const Real_wp& teta);

    void InitSetPoints(const SetPoints<Dimension>& pts);
    const SetPoints<Dimension>& GetSetPoints() const;
    
    //! displays informations about class SetMatrices
    template<class Dim>
    friend ostream& operator <<(ostream& out,const SetMatrices<Dim>& e);
    
  };

  template<class Dim>
  ostream& operator <<(ostream& out,const SetMatrices<Dim>& e);
  
  //! generic class to compute the inverse of transformation Fi
  template<class Dimension>
  class FjInverseProblem : public NonLinearEquations_Base<Real_wp, typename Dimension::R_N,
							  typename Dimension::MatrixN_N>,
			   public NonLinearEquations_Newton<Real_wp, typename Dimension::R_N,
							    typename Dimension::R_N>,
			   public VirtualLeastSquaresFunction<Real_wp, typename Dimension::R_N,
							      typename Dimension::R_N,
							      typename Dimension::MatrixN_N,
							      typename Dimension::MatrixN_Nsym>
  {
  protected :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    
    const ElementGeomReference_Base<Dimension>& Fb; //!< finite element class
    const Mesh<Dimension>& mesh; //!< given mesh    
    const int num_elem; //!< element number
    
    typename Dimension::VectR_N Vertices; //!< vertices of the element
    SetPoints<Dimension> PointsElem; //!< points after transformation F_i
    SetMatrices<Dimension> MatricesElem;
    typename Dimension::R_N point_global; //!< we want to solve F_i(x) = point_global
 
    bool affine, curved;
    MatrixN_N invDF;
    R_N F0, x_eval;
    VectReal_wp phi;
    Real_wp Rmax, norm_DF;

    NewtonSolver<Real_wp, R_N, R_N> newton_solver;

  public :
    
    static int non_linear_solver;
    enum{NEWTON_SOLVER, MINPACK_SOLVER, LVM_SOLVER};
    //! maximum number of iterations to invert \f$ F_i \f$
    static int nb_max_iterations;
    
    static int64_t nb_iter_all;
    //! stopping criterion
    static Real_wp threshold_newton, coef_safety_solver;
    
    FjInverseProblem(const Mesh<Dimension>& mesh, int nquad);
    
    typename Dimension::VectR_N& GetVertices();
    SetPoints<Dimension>& GetSetPoints();

    static void SetInputData(const string&, const VectString&);
    
    void FindInitGuess(R_N& x);    
    void EvaluateFunction(const R_N& x, R_N& fvec);    
    void EvaluateJacobian(const R_N& x, MatrixN_N& fjac);

    void ComputeScheme(const R_N&, R_N&);
    void ComputeAndFactoriseDiff(const R_N&, const R_N&);
    void SolveDifferential(const R_N&, R_N&);
    Real_wp GetNorm2Vector(const R_N&);
    
    void EvaluateF(const R_N& x, R_N& feval);
    void EvaluateJacobian(const R_N& x, R_N& feval, MatrixN_N& fjac);

    bool Solve(const R_N& point, R_N& res, int print_level = 0);
    
  };
    
}
  
#define MONTJOIE_FILE_POINTS_REFERENCE_HXX  
#endif
