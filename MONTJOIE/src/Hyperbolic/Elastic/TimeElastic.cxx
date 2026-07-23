#ifndef MONTJOIE_FILE_TIME_ELASTIC_CXX

#include "ProdMatVectElastic.cxx"

namespace Montjoie
{
  
  /********************
   * TimeElastic_Base *
   ********************/
  
  
  //! default constructor
  template<class TypeEquation>
  TimeElastic_Base<TypeEquation>::TimeElastic_Base() : VarInstationary<TypeEquation>()
  {
  }
  

  template<class TypeEquation>
  bool TimeElastic_Base<TypeEquation>::SplitSystem() const
  {
    return true;
  }
  
  
  // initialization for some variables
  template<class TypeEquation>
  bool TimeElastic_Base<TypeEquation>::IsDampedMedia(int ref)
  {
    return false;
  }
  

  //! initialization before time iterations
  template<class TypeEquation>
  void TimeElastic_Base<TypeEquation>::InitTimeIterations()
  {    
    // use of split PML => jacobian matrices DFi must be diagonal
    this->CheckIdentityJacobianInPML();
    
    VarInstationary<TypeEquation>::InitTimeIterations();
  }
  

  //! compute Y = G(tn,X) if we put elastodynamic equation with pml to the form dX/dt = G(t,X)
  /*! 
    \param[in] tn input time
    \param[in] nb_deriv derivative number of G, if 0 we evaluate G,
    if 1 dG/dt, if 2 d^2 G /dt^2 ...
    \param[in] X iterate
    \param[out] Y result of G(tn,X)
  */
  template<class TypeEquation>
  void TimeElastic_Base<TypeEquation>
  ::EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                               const VectReal_wp& X, VectReal_wp& Y, bool invert_mass, bool source)
  {
    const MeshNumbering<Dimension>& mesh_num = this->var_harmonic.GetMeshNumbering(0);
    VectReal_wp Uh, Vh;
    VectReal_wp Prod_Uh, Prod_Vh;
    int Nvol = mesh_num.GetNbDof();
    int Npml = mesh_num.GetNbDofPML();
    int N = this->var_harmonic.GetNbDof();
    int d = Dimension::dim_N;
    int nb_dof_u = d*Nvol + d*(d-1)*Npml;
    int nb_dof_v = N - nb_dof_u;

    Real_wp one(1), zero(0);
    
    if ((N > X.GetM()) || (N > Y.GetM()))
      {
        cout << "X or Y is not large enough" << endl;
        DISP(N);
        DISP(X.GetM()); DISP(Y.GetM());
        abort();
      }
    
    Real_wp* x_ptr = X.GetData(); Real_wp* y_ptr = Y.GetData();
    Uh.SetData(nb_dof_u, x_ptr); Prod_Uh.SetData(nb_dof_u, y_ptr);
    x_ptr += nb_dof_u; y_ptr += nb_dof_u;
    
    Vh.SetData(nb_dof_v, x_ptr); Prod_Vh.SetData(nb_dof_v, y_ptr);
    
    // dU/dt = Dh^{-1} ( Rh^S V + Fh - Sh U)
    this->ApplyOperatorRhScalar(one, tn, Vh, zero, Prod_Uh);
    if (source)
      this->AddScalarSourceAtTime(one, tn, nb_deriv, Prod_Uh);
    
    this->ApplyOperatorSh(-one, tn, Uh, one, Prod_Uh);
    
    if (invert_mass)
      this->SolveOperatorDh(Prod_Uh);
    
    // dV/dt = Bh^{-1} ( Rh^V U - ShV V + FhV)  (for elements outside PML)
    this->ApplyOperatorRhVectorial(one, tn, Uh, zero, Prod_Vh);
    if (source)
      this->AddVectorialSourceAtTime(one, tn, nb_deriv, Prod_Vh);
    
    this->ApplyOperatorShVectorial(-one, tn, Vh, one, Prod_Vh);
    if (invert_mass)
      {
        this->SolveOperatorBh(Prod_Vh);
        this->SetDirichletCondition(tn, nb_deriv+1, Prod_Uh);
      }
    
    // nullify temporary vectors
    Uh.Nullify(); Vh.Nullify(); Prod_Uh.Nullify(); Prod_Vh.Nullify();
  }


  //! compute Y = G(tn,X) if we put elastodynamic equation with pml to the form dX/dt = G(t,X)
  /*! 
    \param[in] tn input time
    \param[in] nb_deriv derivative number of G, if 0 we evaluate G,
    if 1 dG/dt, if 2 d^2 G /dt^2 ...
    \param[in] X iterate
    \param[out] Y result of G(tn,X)
  */
  template<class TypeEquation>
  void TimeElastic_Base<TypeEquation>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
                             const VectReal_wp& X, VectReal_wp& Y, bool invert_mass, bool source)
  {
    EvaluateDerivativeFunction(tn, nb_deriv, X, Y, invert_mass, source);
  }


  template<class TypeEquation>
  VirtualMassMatrix* TimeElastic_Base<TypeEquation>::GetNewMassMatrix()
  {
    return new ContinuousMassMatrixElastic<Dimension>(this->GetLeafClass());
  }


  template<class Dimension>
  void HyperbolicProblem<TimeElasticEquation<Dimension> >
  ::MltAddStiffnessScalar(const Real_wp& alpha, int level, const VectReal_wp& B,
			  const Real_wp& beta, VectReal_wp& C)
  {    
    if (this->Glob_mat_Kh == NULL)
      {
        cout << "Stiffness matrix not allocated" << endl;
        abort();
      }

    FemMatrixFreeClass<Real_wp, ElasticEquation<Dimension> >&
      Kh = static_cast<FemMatrixFreeClass<Real_wp, ElasticEquation<Dimension> >& >(*this->Glob_mat_Kh);

    MltAdd_SquareElasHex_ScalarH1(alpha, level, Kh, B, beta, C);
  }

  template<class Dimension>
  void HyperbolicProblem<TimeElasticEquation<Dimension> >
  ::MltAddStiffnessVectorial(const Real_wp& alpha, int level, const VectReal_wp& B,
			     const Real_wp& beta, VectReal_wp& C)
  {
    if (this->Glob_mat_Kh == NULL)
      {
        cout << "Stiffness matrix not allocated" << endl;
        abort();
      }

    FemMatrixFreeClass<Real_wp, ElasticEquation<Dimension> >&
      Kh = static_cast<FemMatrixFreeClass<Real_wp, ElasticEquation<Dimension> >& >(*this->Glob_mat_Kh);

    MltAdd_SquareElasHex_VectorialH1(alpha, level, Kh, B, beta, C);
  }

  

  /********************************
   * ContinuousMassMatrixElastic *
   ********************************/
  
  
  //! constructor
  template<class Dimension> template<class TypeEquation>
  ContinuousMassMatrixElastic<Dimension>
  ::ContinuousMassMatrixElastic(HyperbolicProblem<TypeEquation>& var)
    : ContinuousUnsteadyMassMatrix<Dimension>(var),
      var_problem(var.var_harmonic), var_boundary(var.var_harmonic), var_time(var)
  {
  }
  
  
  //! Computes the mass matrices
  template<class Dimension>
  void ContinuousMassMatrixElastic<Dimension>::Init(bool compute_time)
  {
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    int Nscal = var_problem.offset_dof_unknown(1);
    int Npml = mesh_num.GetNbDofPML();
    int N = var_problem.GetNbDof();
    int d = Dimension::dim_N;
    int nb_dof_u = d*(Nvol + (d-1)*Npml);
    int nb_dof_v = N - nb_dof_u;
    Matrix<Real_wp, Symmetric, DiagonalRow> Dh(N, N);    

    int Npml_all = Npml;
    int rank_proc = 0;
    
#ifdef SELDON_WITH_MPI
    MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    MPI_Allreduce(&Npml, &Npml_all, 1, MPI_INTEGER, MPI_MAX, var_problem.comm_group_mode);
#endif

    Vector<bool> diag_elt;
    int type = var_problem.GetMassMatrixType(diag_elt);
    bool diag_matrices = true;
    if (type != FemMassMatrix::DIAGONAL)
      diag_matrices = false;
    
#ifdef MONTJOIE_WITH_TRANSMISSION
    // non-diagonal mass matrix for transmission conditions
    if (var_boundary.GetTransmissionProblem().GetNbInterfaces() > 0)
      diag_matrices = false;
#endif

    bool compute_Dh, compute_invDh, compute_DhMinusdtSh,
      compute_invDhPlusdtSh, compute_Sh;
    
    this->FindMatricesToCompute(compute_Dh, compute_invDh, compute_DhMinusdtSh,
				compute_invDhPlusdtSh, compute_Sh);
    
    bool compute_Bh, compute_invBh, compute_BhMinusdtSh,
      compute_invBhPlusdtSh, compute_ShVec;
    
    this->FindMatricesToComputeVec(compute_Bh, compute_invBh, compute_BhMinusdtSh,
				   compute_invBhPlusdtSh, compute_ShVec);
    
    Real_wp deltat = var_time.GetTimeStep();
    if (!diag_matrices)
      {
        if (Npml > 0)
          {
            cout << "not implemented" << endl;
            abort();
          }
        else
          {
	    GlobalGenericMatrix<Real_wp> nat_mat;
            nat_mat.SetCoefMass(1.0);
            nat_mat.SetCoefDamping(0.0);
            nat_mat.SetCoefStiffness(0.0);
            this->Dh = this->GetNewMassMatrix(FemMassMatrix::MATRIX_SPARSE);
            var_problem.AddMatrixWithBC(this->Dh->GetSparseMatrix(), nat_mat);
	    
            nat_mat.SetCoefMass(0.0);
            nat_mat.SetCoefDamping(1.0);
            nat_mat.SetCoefStiffness(0.0);
	    this->Sh = this->GetNewMassMatrix(FemMassMatrix::MATRIX_SPARSE);
            var_problem.AddMatrixWithBC(this->Sh->GetSparseMatrix(), nat_mat);
                        
            if (compute_DhMinusdtSh)
              {
                this->DhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::MATRIX_SPARSE);
                this->DhMinusdtSh->GetSparseMatrix() = this->Dh->GetSparseMatrix();
                Add(-0.5*deltat, this->Sh->GetSparseMatrix(), this->DhMinusdtSh->GetSparseMatrix());
              }
          }
        
        return;
      }
    
    // evaluation of mass matrix
    glob_chrono.Reset(VirtualTimer::MASS);
    glob_chrono.Start(VirtualTimer::MASS);
    
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(1.0);
    nat_mat.SetCoefDamping(0.0);
    nat_mat.SetCoefStiffness(0.0);
    var_problem.AddMatrixWithBC(Dh, nat_mat);
    
    glob_chrono.Stop(VirtualTimer::MASS);
    if (rank_proc == 0)
      cout << "Time to compute Dh = " << glob_chrono.GetSeconds(VirtualTimer::MASS) << " seconds" << endl;
    
    if ((compute_Dh) || (compute_invDh))
      this->Dh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);

    if (compute_Dh)
      {
	VectReal_wp& diagonal = this->Dh->GetDiagonal();
        diagonal.Reallocate(nb_dof_u);
        for (int i = 0; i < nb_dof_u; i++)
          diagonal(i) = Dh(i, i);
      }

    if (compute_invDh)
      {
	VectReal_wp& invDiagonal = this->Dh->GetInverseDiagonal();
        invDiagonal.Reallocate(nb_dof_u);
        for (int i = 0; i < nb_dof_u; i++)
          invDiagonal(i) = Dh(i, i);
        
        var_problem.AddDomains(invDiagonal);
        for (int i = 0; i < nb_dof_u; i++)
          invDiagonal(i) = 1.0/invDiagonal(i);
      }
    
    if (Npml_all == 0)
      {
        glob_chrono.Reset(VirtualTimer::MASS);
        glob_chrono.Start(VirtualTimer::MASS);

        Matrix<Real_wp, Symmetric, DiagonalRow> Sh(N, N);            
        //DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> Sh(N, N);
        // in that case, mass matrix is diagonal
        nat_mat.SetCoefMass(0.0);
        nat_mat.SetCoefDamping(1.0);
        nat_mat.SetCoefStiffness(0.0);
        var_problem.AddMatrixWithBC(Sh, nat_mat);
        //Sh.WriteText("Sh.dat");
        
        glob_chrono.Stop(VirtualTimer::MASS);

        if (rank_proc == 0)
          cout << "Time to compute Sh = " << glob_chrono.GetSeconds(VirtualTimer::MASS) << " seconds" << endl;
        
        if (compute_DhMinusdtSh)
          {
            this->DhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    
	    VectReal_wp& diagonal = this->DhMinusdtSh->GetDiagonal();
            diagonal.Reallocate(nb_dof_u);
            for (int i = 0; i < nb_dof_u; i++)
              diagonal(i) = Dh(i, i) - 0.5*deltat*Sh(i, i);
          }
        
        if (compute_invDhPlusdtSh)
          {
            this->DhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& invDiagonal = this->DhPlusdtSh->GetInverseDiagonal();
	    invDiagonal.Reallocate(nb_dof_u);
            for (int i = 0; i < nb_dof_u; i++)
              invDiagonal(i) = Dh(i, i) + 0.5*deltat*Sh(i, i);
            
            var_problem.AddDomains(invDiagonal);
            for (int i = 0; i < nb_dof_u; i++)
              invDiagonal(i) = 1.0/invDiagonal(i);
          }
	
        if (compute_Sh)
          {
	    this->Sh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
            VectReal_wp& diagonal = this->Sh->GetDiagonal();
	    diagonal.Reallocate(nb_dof_u);
            for (int i = 0; i < nb_dof_u; i++)
              diagonal(i) = Sh(i, i);
          }
        
        if ((compute_Bh) || (compute_invBh))
	  this->Bh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	
	if (compute_Bh)
	  {
	    VectReal_wp& diagonal = this->Bh->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Dh(nb_dof_u+i, nb_dof_u+i);          
          }
	
        if (compute_invBh)
          {
	    VectReal_wp& invDiagonal = this->Bh->GetInverseDiagonal();
            invDiagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              invDiagonal(i) = 1.0/Dh(nb_dof_u+i, nb_dof_u+i);          
          }
        
        if (compute_ShVec)
          {
	    this->ShVec = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& diagonal = this->ShVec->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Sh(nb_dof_u+i, nb_dof_u+i);          
          }
	
        if (compute_BhMinusdtSh)
          {
	    this->BhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& diagonal = this->BhMinusdtSh->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Dh(nb_dof_u+i, nb_dof_u+i) - 0.5*deltat*Sh(nb_dof_u+i, nb_dof_u+i);          
          }
	
        if (compute_invBhPlusdtSh)
          {
	    this->BhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& invDiagonal = this->BhPlusdtSh->GetInverseDiagonal();
	    invDiagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              invDiagonal(i)
                = 1.0/(Dh(nb_dof_u+i, nb_dof_u+i) + 0.5*deltat*Sh(nb_dof_u+i, nb_dof_u+i));
          }
      }
    else
      {
        Vector<IVect> pattern(d*Nvol + nb_dof_v);
        DistributedMatrix<Real_wp, General, BlockDiagRow> Sh;
        
        for (int i = 0; i < Nvol; i++)
          {
            int num = mesh_num.GetDofPML(i);
            if (num >= 0)
              {
                pattern(i).Reallocate(d);
                pattern(i)(0) = i;
                pattern(i)(1) = Nvol + num;
                if (d == 3)
                  pattern(i)(2) = Nvol + Npml + num;
                
                pattern(i+Nvol).Reallocate(d);
                pattern(i+Nvol)(0) = i + Nscal;
                pattern(i+Nvol)(1) = Nscal + Nvol + num;
                if (d == 3)
                  {
                    pattern(i+Nvol)(2) = Nscal + Nvol + Npml + num;

                    pattern(i+2*Nvol).Reallocate(d);
                    pattern(i+2*Nvol)(0) = i + 2*Nscal;
                    pattern(i+2*Nvol)(1) = 2*Nscal + Nvol + num;
                    pattern(i+2*Nvol)(2) = 2*Nscal + Nvol + Npml + num;
                  }                
              }
            else
              {
                pattern(i).Reallocate(1);
                pattern(i)(0) = i;
                pattern(Nvol+i).Reallocate(1);
                pattern(Nvol+i)(0) = Nscal+i;
                if (d == 3)
                  {
                    pattern(2*Nvol+i).Reallocate(1);
                    pattern(2*Nvol+i)(0) = 2*Nscal+i; 
                  }
              }
          }
        
        for (int i = 0; i < nb_dof_v; i++)
          {
            pattern(d*Nvol+i).Reallocate(1);
            pattern(d*Nvol+i)(0) = nb_dof_u+i;
          }

        var_problem.InitDistributedMatrix(Sh);        
        Sh.SetPattern(pattern);
        pattern.Resize(d*Nvol);
        
        glob_chrono.Reset(VirtualTimer::MASS);
        glob_chrono.Start(VirtualTimer::MASS);
        
        nat_mat.SetCoefMass(0.0);
        nat_mat.SetCoefDamping(1.0);
        nat_mat.SetCoefStiffness(0.0);
        var_problem.AddMatrixWithBC(Sh, nat_mat);
	
        glob_chrono.Stop(VirtualTimer::MASS);
        if (rank_proc == 0)
          cout << "Time to compute Sh = " << glob_chrono.GetSeconds(VirtualTimer::MASS) << " seconds" << endl;
        
	this->Sh = this->GetNewMassMatrix(FemMassMatrix::BLOCK_DIAGONAL_UNSYM);
	DistributedMatrix<Real_wp, General, BlockDiagRow>&
	  Sh_block_diagonal = this->Sh->GetUnsymmetricBlockDiagonal();

        var_problem.InitDistributedMatrix(Sh_block_diagonal);	
        Sh_block_diagonal.SetPattern(pattern);
        for (int i = 0; i < d*Nvol; i++)
          {
            int n = Sh_block_diagonal.GetBlockSize(i);
            for (int j = 0; j < n; j++)
              for (int k = 0; k < n; k++)
                Sh_block_diagonal.Value(i, j, k) = Sh.Value(i, j, k);
          }
        
        if (compute_DhMinusdtSh)
          {
	    this->DhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::BLOCK_DIAGONAL_UNSYM);
	    DistributedMatrix<Real_wp, General, BlockDiagRow>&
	      block_diagonal = this->DhMinusdtSh->GetUnsymmetricBlockDiagonal();

            block_diagonal = Sh_block_diagonal;            
            Mlt(-0.5*deltat, block_diagonal);
            for (int i = 0; i < nb_dof_u; i++)
	      block_diagonal.AddInteraction(i, i, Dh(i, i));          
          }

        if (compute_invDhPlusdtSh)
          {
	    this->DhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::BLOCK_DIAGONAL_UNSYM);
	    DistributedMatrix<Real_wp, General, BlockDiagRow>&
	      block_diagonal = this->DhPlusdtSh->GetUnsymmetricBlockDiagonal();
	    
            block_diagonal = Sh_block_diagonal;            
            Mlt(0.5*deltat, block_diagonal);
            for (int i = 0; i < nb_dof_u; i++)
	      block_diagonal.AddInteraction(i, i, Dh(i, i));          
	    
            this->DhPlusdtSh->Invert(false);
          }
        
        if (!compute_Sh)
	  Sh_block_diagonal.Clear();
	
        if ((compute_Bh) || (compute_invBh))
	  this->Bh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	
        if (compute_Bh)
          {
	    VectReal_wp& diagonal = this->Bh->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Dh(nb_dof_u+i, nb_dof_u+i);          
          }
        
        if (compute_invBh)
          {
	    VectReal_wp& invDiagonal = this->Bh->GetInverseDiagonal();
	    invDiagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              invDiagonal(i) = 1.0/Dh(nb_dof_u+i, nb_dof_u+i);          
          }
        
        if (compute_ShVec)
          {
	    this->ShVec = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& diagonal = this->ShVec->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Sh(nb_dof_u+i, nb_dof_u+i);          
          }
        
        if (compute_BhMinusdtSh)
          {
	    this->BhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& diagonal = this->BhMinusdtSh->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Dh(nb_dof_u+i, nb_dof_u+i) - 0.5*deltat*Sh(nb_dof_u+i, nb_dof_u+i);          
          }
        
        if (compute_invBhPlusdtSh)
          {
	    this->BhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& invDiagonal = this->BhPlusdtSh->GetInverseDiagonal();
            invDiagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              invDiagonal(i) = 1.0/
		(Dh(nb_dof_u+i, nb_dof_u+i) + 0.5*deltat*Sh(nb_dof_u+i, nb_dof_u+i));
          }
      }
  }
  
}

#define MONTJOIE_FILE_TIME_ELASTIC_CXX
#endif
