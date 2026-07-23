#ifndef MONTJOIE_FILE_DISCONTINUOUS_UNSTEADY_MASS_MATRIX_HXX

namespace Montjoie
{
  
  //! class storing the mass matrix, damping matrix or combinations for time schemes 
  //! and discontinuous formulation
  /*!
    This class is able to store Dh (mass matrix), Sh (damping matrix)
    and combinations Dh - dt/2 Sh, Dh + dt/2 Sh
    For first-order formulations, it can also store Bh (mass matrix of additional unknowns), ShV
    and combinations Bh - dt/2 ShV, Bh + dt/2 ShV
   */
  class DiscontinuousUnsteadyMassMatrix_Base : public VirtualMassMatrix
  {
  protected :
    
    // mass matrix with only jacobian
    DiscontinuousMassMatrix_Base* GeomMass;
    
    // different mass matrices we need for scalar unknown
    DiscontinuousMassMatrix_Base* Dh, *DhMinusdtSh, *DhPlusdtSh, *Sh;
    
    // and mass matrices for vectorial unknown
    DiscontinuousMassMatrix_Base* Bh, *BhMinusdtSh, *BhPlusdtSh, *ShVec;

    // surface mass matrix (for HDG)
    DiscontinuousMassMatrix_Base* SurfMass, *ChSurf;
    
  private:
    VarComputationProblem& var_computation;
    VarProblem_Base& var_problem;
    VarBoundaryCondition_Base& var_boundary;
    VarInstationary_Base& var_time;
    DistributedProblem_Base& var_comm;
    
    void InitDefaultValues();
    
  public :    
    
    template<class TypeEquation>
    DiscontinuousUnsteadyMassMatrix_Base(HyperbolicProblem<TypeEquation>& vars);

    ~DiscontinuousUnsteadyMassMatrix_Base();
    
    void Clear();
    
    // computation of the mass matrix
    void Init(bool compute_time = true);
    void InitSurface(bool compute_mass = true);

    virtual size_t GetMemorySize() const;
    
    // inline methods
    const VectReal_wp& GetDiagonalDh() const;
    const VectReal_wp& GetInverseDiagonalDh() const;
    const VectReal_wp& GetDiagonalSh() const;
    const VectReal_wp& GetDiagonalDhMinusdtSh() const;
    const VectReal_wp& GetInverseDiagonalDhPlusdtSh() const;
    const VectReal_wp& GetDiagonalBh() const;
    const VectReal_wp& GetDiagonalBhMinusdtSh() const;
    const VectReal_wp& GetInverseDiagonalBhPlusdtSh() const;

    DiscontinuousMassMatrix_Base& GetOperatorDh();
    DiscontinuousMassMatrix_Base& GetOperatorBh();
    DiscontinuousMassMatrix_Base& GetOperatorSh();
    DiscontinuousMassMatrix_Base& GetOperatorShVec();
    DiscontinuousMassMatrix_Base& GetOperatorCh();
    
    virtual DiscontinuousMassMatrix_Base* GetNewMassMatrix(int type_mat) = 0;
    virtual DiscontinuousMassMatrix_Base* GetNewSurfaceMassMatrix(int type_mat, bool geom) = 0;

    // other methods
    void SolveMassMatrix(VectReal_wp& Y);    
    void SolveMassMatrix(VectReal_wp& Y, int m1, int m2);

    void ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t,
			 const VectReal_wp& x, const Real_wp& beta, VectReal_wp& Y);    
    
    void ApplyOperatorSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			 const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void ApplyOperatorShVectorial(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
				  const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			 const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
				    const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void SolveOperatorDhPlusdtSh(const Real_wp&, const VectReal_wp& X, VectReal_wp& Y);
    void SolveOperatorDhPlusdtSh(VectReal_wp& Y);
    
    void SolveOperatorDh(VectReal_wp& Y);
    
    void SolveCholeskyDh(const SeldonTranspose& transA, VectReal_wp& Y);
    
    void SolveOperatorDh(const Real_wp&, const VectReal_wp& X, VectReal_wp& Y);
    
    void SolveOperatorBh(VectReal_wp& Vh);
    void SolveOperatorBh(const Real_wp&, const VectReal_wp& Uh, VectReal_wp& Vh);

    void ApplyOperatorBh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
                         const Real_wp& beta, VectReal_wp& Vh);
    
    void ApplyOperatorBhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
                                  const Real_wp& beta, VectReal_wp& Vh);
    
    void SolveOperatorBhPlusdtSh(const Real_wp&, const VectReal_wp& X, VectReal_wp& Y);
    void SolveOperatorBhPlusdtSh(VectReal_wp& Y);

    void SolveOperatorCh(VectReal_wp& Y);
    
  };
  
  
  template<class Dimension>  
  class DiscontinuousUnsteadyMassMatrix : public DiscontinuousUnsteadyMassMatrix_Base
  {
  private:
    VarProblem<Dimension>& var_problem;
    VarInstationary_Fem<Dimension>& var_time;
    
  public :
    template<class TypeEquation>
    DiscontinuousUnsteadyMassMatrix(HyperbolicProblem<TypeEquation>& var);

    DiscontinuousMassMatrix_Base* GetNewMassMatrix(int type_mat);
    DiscontinuousMassMatrix_Base* GetNewSurfaceMassMatrix(int type_mat, bool geom);
    
  };

} // end namespace

#define MONTJOIE_FILE_DISCONTINUOUS_UNSTEADY_MASS_MATRIX_HXX
#endif
