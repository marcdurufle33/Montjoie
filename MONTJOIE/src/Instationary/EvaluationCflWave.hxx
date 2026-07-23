#ifndef MONTJOIE_FILE_EVALUATION_CFL_WAVE_HXX

namespace Montjoie
{

  class EigenvalueComparisonStabFunc : public EigenvalueComparisonClass<Real_wp>
  {
  private:
    Real_wp dt_nom;
    UnivariatePolynomial<Real_wp> Pol;
    
  public:
    void SetFunction(const Real_wp& dt, const UnivariatePolynomial<Real_wp>& P);

    int CompareEigenvalue(const Real_wp& Lr, const Real_wp& Li, const Real_wp& Lr2, const Real_wp& Li2);
    
    Real_wp FindTimeStep(const Real_wp& Lr, const Real_wp& Li) const;

  };
  

  //! class to evaluate large eigenvalues of M^-1 K
  class EigenProblemTimeMontjoie
    : public VirtualEigenProblem<Real_wp, Real_wp, Real_wp>      
  {
  private :
    VarProblem_Base& var_problem;
    VarBoundaryCondition_Base& var_boundary;
    VarInstationary_Base& var_time;
    DistributedProblem_Base& var_comm;
    
    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    EigenProblemTimeMontjoie(HyperbolicProblem<TypeEquation>& var);
    
    void FactorizeCholeskyMass();    
    void ComputeDiagonalMass();    

    void MltStiffness(const Real_wp&, const Real_wp&,
                      const Vector<Real_wp>& X, Vector<Real_wp>& Y);

    void MltStiffness(const Vector<Real_wp>& X, Vector<Real_wp>& Y);
    void MltStiffness(const SeldonTranspose&, const Vector<Real_wp>& X, Vector<Real_wp>& Y);

    void MltStiffness(const Real_wp&, const Real_wp&,
                      const Vector<Complex_wp>& X, Vector<Complex_wp>& Y);

    void MltStiffness(const Vector<Complex_wp>& X, Vector<Complex_wp>& Y);
    void MltStiffness(const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y);
    
    bool IsSymmetricProblem() const;
    bool IsHermitianProblem() const;
    
    void CheckStiffnessMatrix();
    
  };
  
  
  template<class Dimension>
  class MeshGeometryCFL
  {
  protected:
    Vector<int> IndexVertex;
    Vector<bool> ElementOnSubdomain;

    Vector<int> IndexNeighbor;
    typename Dimension::VectR_N PointsNeighbor;
    Vector<IVect> ElementNeighbor; 
    Vector<int> ReferenceNeighbor;
    
    Vector<int> NeighIndexVertex;
    Vector<bool> NeighElementOnSubdomain;
    
  public :
    void Init(const DistributedProblem<Dimension>& var);
    
    void ExtractSubMesh(const Mesh<Dimension>& mesh, const MeshNumbering<Dimension>&,
                        Mesh<Dimension>& sub_mesh, MeshNumbering<Dimension>&,
                        const IVect& num, const IVect& num_vertex,
                        const IVect& num_neigh, const IVect& num_vertex_neigh);
    
    void CreateSmallMesh(int i, const Mesh<Dimension>& mesh, const MeshNumbering<Dimension>& mesh_num,
                         Mesh<Dimension>& small_mesh, MeshNumbering<Dimension>& small_mesh_num);
    
  };
  
  
  template<class TypeEquation>
  Real_wp EvaluateCFL(HyperbolicProblem<TypeEquation>& var_time,
                      const UnivariatePolynomial<Real_wp>& P,
                      bool init_computation = true);
  
  Real_wp EvaluateCFL_Unsteady(VarComputationProblem& var, VarInstationary_Base& var_time,
			       EigenProblemTimeMontjoie& pb_eig,
			       const UnivariatePolynomial<Real_wp>&, bool init_computation = true);
  
  Real_wp EvaluateCFL_Generic(EigenProblemTimeMontjoie&, const UnivariatePolynomial<Real_wp>&,
                              const Real_wp& dt, bool first_order);
  
  template<class Dimension>
  void ComputeLocalTimeStep(VarInstationary_Dim<Dimension>& var_time,
                            VectReal_wp& dt);
  
}

#define MONTJOIE_FILE_EVALUATION_CFL_WAVE_HXX
#endif

