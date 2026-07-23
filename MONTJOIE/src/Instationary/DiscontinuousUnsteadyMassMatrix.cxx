#ifndef MONTJOIE_FILE_DISCONTINUOUS_UNSTEADY_MASS_MATRIX_CXX

namespace Montjoie
{
  
  void DiscontinuousUnsteadyMassMatrix_Base::InitDefaultValues()
  {
    GeomMass = NULL;
    Dh = NULL;
    DhMinusdtSh = NULL;
    DhPlusdtSh = NULL;
    Sh = NULL;
    Bh = NULL;
    BhMinusdtSh = NULL;
    BhPlusdtSh = NULL;
    ShVec = NULL;
    SurfMass = NULL;
    ChSurf = NULL;
  }

  
  DiscontinuousUnsteadyMassMatrix_Base::~DiscontinuousUnsteadyMassMatrix_Base()
  {
    Clear();
  }


  void DiscontinuousUnsteadyMassMatrix_Base::Clear()
  {
    if (Dh != NULL)
      delete Dh;
  
    if (DhMinusdtSh != NULL)
      delete DhMinusdtSh;

    if (DhPlusdtSh != NULL)
      delete DhPlusdtSh;

    if (Sh != NULL)
      delete Sh;

    if (Bh != NULL)
      delete Bh;
  
    if (BhMinusdtSh != NULL)
      delete BhMinusdtSh;

    if (BhPlusdtSh != NULL)
      delete BhPlusdtSh;

    if (ShVec != NULL)
      delete ShVec;

    if (GeomMass != NULL)
      delete GeomMass;

    if (SurfMass != NULL)
      delete SurfMass;

    if (ChSurf != NULL)
      delete ChSurf;

    InitDefaultValues();
  }

  
  //! computation of needed mass matrices Dh, Sh, Bh etc for the chosen time scheme
  void DiscontinuousUnsteadyMassMatrix_Base::Init(bool compute_time)
  {    
    // checking which matrices have to be computed
    bool compute_Dh, compute_invDh, compute_DhMinusdtSh,
      compute_invDhPlusdtSh, compute_Sh, compute_Bh, compute_invBh,
      compute_BhMinusdtSh, compute_invBhPlusdtSh, compute_ShVec;
    
    this->FindMatricesToCompute(compute_Dh, compute_invDh, compute_DhMinusdtSh,
				compute_invDhPlusdtSh, compute_Sh);
    
    this->FindMatricesToComputeVec(compute_Bh, compute_invBh, compute_BhMinusdtSh,
				   compute_invBhPlusdtSh, compute_ShVec);
    
    // allocating matrices
    int type_matrix = FemMassMatrix::BLOCK_ELT;
    GeomMass = GetNewMassMatrix(type_matrix);
    if (compute_invDh || compute_Dh)
      {
	Dh = GetNewMassMatrix(type_matrix);
	Dh->SetCoefficient(1.0, 0.0);	
      }
    
    if (compute_DhMinusdtSh)
      {
	DhMinusdtSh = GetNewMassMatrix(type_matrix);
	DhMinusdtSh->SetCoefficient(1.0, -0.5*var_time.GetTimeStep());
      }

    if (compute_invDhPlusdtSh)
      {
	DhPlusdtSh = GetNewMassMatrix(type_matrix);
	DhPlusdtSh->SetCoefficient(1.0, 0.5*var_time.GetTimeStep());
      }

    if (compute_Sh)
      {
	Sh = GetNewMassMatrix(type_matrix);
	Sh->SetCoefficient(0.0, 1.0);	
      }

    if (compute_invBh || compute_Bh)
      {
	Bh = GetNewMassMatrix(type_matrix);
	Bh->SetCoefficient(1.0, 0.0);	
      }
    
    if (compute_BhMinusdtSh)
      {
	BhMinusdtSh = GetNewMassMatrix(type_matrix);
	BhMinusdtSh->SetCoefficient(1.0, -0.5*var_time.GetTimeStep());
      }

    if (compute_invBhPlusdtSh)
      {
	BhPlusdtSh = GetNewMassMatrix(type_matrix);
	BhPlusdtSh->SetCoefficient(1.0, 0.5*var_time.GetTimeStep());
      }

    if (compute_ShVec)
      {
	ShVec = GetNewMassMatrix(type_matrix);
	ShVec->SetCoefficient(0.0, 1.0);	
      }
    
    // computation of mass matrices for explicit schemes
    int jeton = glob_chrono.GetNumber();
    if (compute_time)
      glob_chrono.Start(jeton);
    
    GeomMass->ComputeMass(true);
    
    if (compute_time)
      {
	glob_chrono.Stop(jeton);    
	glob_chrono.DisplayTime(jeton);
      }
    
    glob_chrono.ReleaseNumber(jeton);
    
    // considering only isotropic coefficients
    int ms = var_problem.nb_unknowns_scal;
    int mv = var_problem.nb_unknowns_scal + var_problem.nb_unknowns_vec;
    if ((compute_Dh)||(compute_invDh))
      Dh->ComputeMass(false, 0, ms);
    
    if (compute_DhMinusdtSh)
      DhMinusdtSh->ComputeMass(false, 0, ms);
    
    if (compute_invDhPlusdtSh)
      DhPlusdtSh->ComputeMass(false, 0, ms);
    
    if (compute_Sh)
      Sh->ComputeMass(false, 0, ms);
          
    if ((compute_Bh)||(compute_invBh))
      Bh->ComputeMass(false, ms, mv);
    
    if (compute_BhMinusdtSh)
      BhMinusdtSh->ComputeMass(false, ms, mv);
    
    if (compute_invBhPlusdtSh)
      BhPlusdtSh->ComputeMass(false, ms, mv);

    if (compute_ShVec)
      ShVec->ComputeMass(false, ms, mv);

    InitSurface();
  }


  //! computation of mass matrices for surface integrals
  void DiscontinuousUnsteadyMassMatrix_Base::InitSurface(bool compute_mass)
  {
    int type_matrix = FemMassMatrix::BLOCK_ELT;
    // computing surface mass matrix
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	SurfMass = GetNewSurfaceMassMatrix(type_matrix, true);
	SurfMass->SetCoefficient(Real_wp(1), Real_wp(0));
	if (compute_mass)
	  SurfMass->ComputeMass(true);

	ChSurf = GetNewSurfaceMassMatrix(type_matrix, false);
	ChSurf->SetCoefficient(Real_wp(1), Real_wp(0));
	if (compute_mass)
	  ChSurf->ComputeMass(false);
      }
  }
  

  //! returns the size used by the object in bytes
  size_t DiscontinuousUnsteadyMassMatrix_Base::GetMemorySize() const
  {
    size_t taille = 0;

    if (GeomMass != NULL)
      taille += GeomMass->GetMemorySize();

    if (Dh != NULL)
      taille += Dh->GetMemorySize();

    if (DhMinusdtSh != NULL)
      taille += DhMinusdtSh->GetMemorySize();
    
    if (DhPlusdtSh != NULL)
      taille += DhPlusdtSh->GetMemorySize();
    
    if (Sh != NULL)
      taille += Sh->GetMemorySize();

    if (Bh != NULL)
      taille += Bh->GetMemorySize();

    if (BhMinusdtSh != NULL)
      taille += BhMinusdtSh->GetMemorySize();
    
    if (BhPlusdtSh != NULL)
      taille += BhPlusdtSh->GetMemorySize();
    
    if (ShVec != NULL)
      taille += ShVec->GetMemorySize();

    if (SurfMass != NULL)
      taille += SurfMass->GetMemorySize();
    
    if (ChSurf != NULL)
      taille += ChSurf->GetMemorySize();

    return taille;
  }

  
  //! solving by mass matrix, Y = M^{-1} Y, where M regroups Dh and Bh
  void DiscontinuousUnsteadyMassMatrix_Base
  ::SolveMassMatrix(VectReal_wp& Y)
  {
    SolveMassMatrix(Y, 0, var_problem.nb_unknowns);
  }

  
  //! applies both Dh and Bh
  void DiscontinuousUnsteadyMassMatrix_Base
  ::ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t,
		    const VectReal_wp& X, const Real_wp& beta, VectReal_wp& Y)
  {
    if (beta == Real_wp(0))
      Y.Zero();
    else if (beta != Real_wp(1))
      Y *= beta;

    // glob_chrono.Start(VirtualTimer::MASS);

    int offset = 0;
    if (Dh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Dh->MltMass(offset, alpha, X, Y);
    else
      {
	for (int m = 0; m < var_problem.nb_unknowns_scal; m++)
	  {
	    Dh->SetUnknownNumber(m);
	    offset = m*var_comm.GetNbMainUnknownDof();
	    Dh->MltMass(offset, alpha, X, Y);
	  }
      }
    
    int nodl_mesh = var_problem.GetNbMeshDof();
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      nodl_mesh = var_comm.GetNbMainUnknownDof();

    int offset_u = var_problem.nb_unknowns_scal*var_comm.GetNbMainUnknownDof();
    if (Bh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Bh->MltMass(offset_u, alpha, X, Y);
    else
      {
	for (int m = var_problem.nb_unknowns_scal; m < var_problem.nb_unknowns; m++)
	  {
	    Bh->SetUnknownNumber(m);
	    offset = offset_u + (m-var_problem.nb_unknowns_scal)*nodl_mesh;
	    Bh->MltMass(offset, alpha, X, Y);
	  }
      }

    // glob_chrono.Stop(VirtualTimer::MASS);    
  }


  //! solving by mass matrix for a subset of unknowns
  void DiscontinuousUnsteadyMassMatrix_Base
  ::SolveMassMatrix(VectReal_wp& Y, int m1, int m2)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    int offset = 0;
    if (Dh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Dh->SolveMass(offset, Y);
    else
      {
	for (int m = m1; m < var_problem.nb_unknowns_scal; m++)
	  {
	    Dh->SetUnknownNumber(m);
	    offset = m*var_comm.GetNbMainUnknownDof();
	    Dh->SolveMass(offset, Y);
	  }
      }
    
    int nodl_mesh = var_problem.GetNbMeshDof();
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      nodl_mesh = var_comm.GetNbMainUnknownDof();
    
    int offset_u = var_problem.nb_unknowns_scal*var_comm.GetNbMainUnknownDof();
    if (Bh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Bh->SolveMass(offset_u, Y);
    else
      {
	for (int m = var_problem.nb_unknowns_scal; m < m2; m++)
	  {
	    Bh->SetUnknownNumber(m);
	    offset = offset_u + (m-var_problem.nb_unknowns_scal)*nodl_mesh;
	    Bh->SolveMass(offset, Y);
	  }
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! applies operator Sh, Prod_Uh = beta Prod_Uh + alpha Sh Uh  
  void DiscontinuousUnsteadyMassMatrix_Base
  ::ApplyOperatorSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
		    const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    if (beta == Real_wp(0))
      Prod_Uh.Zero();
    else if (beta != Real_wp(1))
      Mlt(beta, Prod_Uh);

    if (Sh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Sh->MltMass(0, alpha, Uh, Prod_Uh);
    else
      {
        int offset = 0;
        for (int m = 0; m < var_problem.nb_unknowns_scal; m++)
          {
            Sh->SetUnknownNumber(m);
            offset = var_problem.offset_dof_unknown(m);
            Sh->MltMass(offset, alpha, Uh, Prod_Uh);
          }
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! applies operator Dh, Prod_Uh = beta Prod_Uh + alpha Dh Uh
  void DiscontinuousUnsteadyMassMatrix_Base
  ::ApplyOperatorDh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
		  const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    if (beta == Real_wp(0))
      Prod_Uh.Zero();
    else if (beta != Real_wp(1))
      Mlt(beta, Prod_Uh);
    
    int offset = 0;
    if (Dh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Dh->MltMass(0, alpha, Uh, Prod_Uh);
    else
      {
        for (int m = 0; m < var_problem.nb_unknowns_scal; m++)
          {
            Dh->SetUnknownNumber(m);
            offset = var_problem.offset_dof_unknown(m);
            Dh->MltMass(offset, alpha, Uh, Prod_Uh);
          }
        
        // glob_chrono.Stop(VirtualTimer::MASS);
      }
  }
  

  //! applies operator Dh - dt/2 Sh, Prod_Uh = beta Prod_Uh + alpha (Dh - dt/2 Sh) Uh  
  void DiscontinuousUnsteadyMassMatrix_Base::
  ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
                           const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    if (beta == Real_wp(0))
      Prod_Uh.Zero();
    else if (beta != Real_wp(1))
      Mlt(beta, Prod_Uh);

    int offset = 0;
    if (DhMinusdtSh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      DhMinusdtSh->MltMass(0, alpha, Uh, Prod_Uh);
    else
      {
        for (int m = 0; m < var_problem.nb_unknowns_scal; m++)
          {
            DhMinusdtSh->SetUnknownNumber(m);
            offset = var_problem.offset_dof_unknown(m);
            DhMinusdtSh->MltMass(offset, alpha, Uh, Prod_Uh);
          }    
      }

    // glob_chrono.Stop(VirtualTimer::MASS);
  }

  
  //! solves by operator Dh + dt/2 Sh, Y = Y + alpha (Dh + dt/2)^{-1} X
  void DiscontinuousUnsteadyMassMatrix_Base::
  SolveOperatorDhPlusdtSh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    int offset = 0;
    if (DhPlusdtSh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      DhPlusdtSh->SolveMass(0, alpha, X, Y);
    else
      {
        for (int m = 0; m < var_problem.nb_unknowns_scal; m++)
          {
            DhPlusdtSh->SetUnknownNumber(m);
            offset = var_problem.offset_dof_unknown(m);
            DhPlusdtSh->SolveMass(offset, alpha, X, Y);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! solves by operator Dh + dt/2 Sh, Y = (Dh + dt/2)^{-1} Y
  void DiscontinuousUnsteadyMassMatrix_Base
  ::SolveOperatorDhPlusdtSh(VectReal_wp& Y)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    int offset = 0;
    if (DhPlusdtSh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      DhPlusdtSh->SolveMass(0, Y);
    else
      {
        for (int m = 0; m < var_problem.nb_unknowns_scal; m++)
          {
            DhPlusdtSh->SetUnknownNumber(m);
            offset = var_problem.offset_dof_unknown(m);
            DhPlusdtSh->SolveMass(offset, Y);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! solves by operator Dh, Y = Dh^{-1} Y  
  void DiscontinuousUnsteadyMassMatrix_Base::
  SolveOperatorDh(VectReal_wp& Y)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    int offset = 0;
    if (Dh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Dh->SolveMass(0, Y);
    else
      {
        for (int m = 0; m < var_problem.nb_unknowns_scal; m++)
          {
            Dh->SetUnknownNumber(m);
            offset = var_problem.offset_dof_unknown(m);
            Dh->SolveMass(offset, Y);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! solves by Cholesky factor of Dh, Y = L^{-1} Y or Y = L^{-T} Y where Dh = L L^T  
  void DiscontinuousUnsteadyMassMatrix_Base::
  SolveCholeskyDh(const SeldonTranspose& transA, VectReal_wp& Y)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    int offset = 0;
    if (Dh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Dh->SolveCholesky(0, transA, Y);
    else
      {
        for (int m = 0; m < var_problem.nb_unknowns_scal; m++)
          {
            Dh->SetUnknownNumber(m);
            offset = var_problem.offset_dof_unknown(m);
            Dh->SolveCholesky(offset, transA, Y);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! solves by operator Dh, Y = Y + alpha Dh^{-1} X  
  void DiscontinuousUnsteadyMassMatrix_Base::
  SolveOperatorDh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    int offset = 0;
    if (Dh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Dh->SolveMass(0, alpha, X, Y);
    else
      {
        for (int m = 0; m < var_problem.nb_unknowns_scal; m++)
          {
            Dh->SetUnknownNumber(m);
            offset = var_problem.offset_dof_unknown(m);
            Dh->SolveMass(offset, alpha, X, Y);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! applies operator Bh, Vh = beta Vh + alpha Bh Uh      
  void DiscontinuousUnsteadyMassMatrix_Base::
  ApplyOperatorBh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
                  const Real_wp& beta, VectReal_wp& Vh)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    if (beta == Real_wp(0))
      Vh.Zero();
    else if (beta != Real_wp(1))
      Mlt(beta, Vh);
   
    int offset = 0;
    if (Bh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Bh->MltMass(0, alpha, Uh, Vh);
    else
      {
        int nodl_mesh = var_problem.GetNbMeshDof();
        if (var_problem.FormulationDG() == ElementReference_Base::HDG)
          nodl_mesh = var_comm.GetNbMainUnknownDof();

        for (int m = 0; m < var_problem.nb_unknowns_vec; m++)
          {
            Bh->SetUnknownNumber(m + var_problem.nb_unknowns_scal);
            offset = m*nodl_mesh;
            Bh->MltMass(offset, alpha, Uh, Vh);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }


  //! applies operator (Bh - dt/2 ShV), Vh = beta Vh + alpha (Bh - dt/2 ShV) Uh      
  void DiscontinuousUnsteadyMassMatrix_Base::
  ApplyOperatorBhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
                           const Real_wp& beta, VectReal_wp& Vh)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    if (beta == Real_wp(0))
      Vh.Zero();
    else if (beta != Real_wp(1))
      Mlt(beta, Vh);
   
    int offset = 0;
    if (BhMinusdtSh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      BhMinusdtSh->MltMass(0, alpha, Uh, Vh);
    else
      {
        int nodl_mesh = var_problem.GetNbMeshDof();
        if (var_problem.FormulationDG() == ElementReference_Base::HDG)
          nodl_mesh = var_comm.GetNbMainUnknownDof();

        for (int m = 0; m < var_problem.nb_unknowns_vec; m++)
          {
            BhMinusdtSh->SetUnknownNumber(m + var_problem.nb_unknowns_scal);
            offset = m*nodl_mesh;
            BhMinusdtSh->MltMass(offset, alpha, Uh, Vh);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! solves operator (Bh + dt/2 ShV), Vh = (Bh + dt/2 ShV)^{-1} Vh    
  void DiscontinuousUnsteadyMassMatrix_Base::
  SolveOperatorBhPlusdtSh(const Real_wp& alpha, const VectReal_wp& Uh, VectReal_wp& Vh)
  {
    // glob_chrono.Start(VirtualTimer::MASS);

    int offset = 0;
    if (BhPlusdtSh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      BhPlusdtSh->SolveMass(0, alpha, Uh, Vh);
    else
      {
        int nodl_mesh = var_problem.GetNbMeshDof();
        if (var_problem.FormulationDG() == ElementReference_Base::HDG)
          nodl_mesh = var_comm.GetNbMainUnknownDof();
        
        for (int m = 0; m < var_problem.nb_unknowns_vec; m++)
          {
            BhPlusdtSh->SetUnknownNumber(m + var_problem.nb_unknowns_scal);
            offset = m*nodl_mesh;
            BhPlusdtSh->SolveMass(offset, alpha, Uh, Vh);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! solves operator (Bh + dt/2 ShV), Vh = (Bh + dt/2 ShV)^{-1} Vh
  void DiscontinuousUnsteadyMassMatrix_Base::
  SolveOperatorBhPlusdtSh(VectReal_wp& Vh)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    int offset = 0;
    if (BhPlusdtSh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      BhPlusdtSh->SolveMass(0, Vh);
    else
      {
        int nodl_mesh = var_problem.GetNbMeshDof();
        if (var_problem.FormulationDG() == ElementReference_Base::HDG)
          nodl_mesh = var_comm.GetNbMainUnknownDof();

        for (int m = 0; m < var_problem.nb_unknowns_vec; m++)
          {
            BhPlusdtSh->SetUnknownNumber(m + var_problem.nb_unknowns_scal);
            offset = m*nodl_mesh;
            BhPlusdtSh->SolveMass(offset, Vh);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! solves operator Bh, Vh = Bh^{-1} Vh
  void DiscontinuousUnsteadyMassMatrix_Base
  ::SolveOperatorBh(VectReal_wp& Vh)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    int offset = 0;
    if (Bh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Bh->SolveMass(0, Vh);
    else
      {
        int nodl_mesh = var_problem.GetNbMeshDof();
        if (var_problem.FormulationDG() == ElementReference_Base::HDG)
          nodl_mesh = var_comm.GetNbMainUnknownDof();

        for (int m = 0; m < var_problem.nb_unknowns_vec; m++)
          {
            Bh->SetUnknownNumber(m + var_problem.nb_unknowns_scal);
            offset = m*nodl_mesh;
            Bh->SolveMass(offset, Vh);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! solves operator Bh, Vh = Vh + alpha Bh^{-1} Uh  
  void DiscontinuousUnsteadyMassMatrix_Base::
  SolveOperatorBh(const Real_wp& alpha, const VectReal_wp& Uh, VectReal_wp& Vh)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    int offset = 0;
    if (Bh->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      Bh->SolveMass(0, alpha, Uh, Vh);
    else
      {
        int nodl_mesh = var_problem.GetNbMeshDof();
        if (var_problem.FormulationDG() == ElementReference_Base::HDG)
          nodl_mesh = var_comm.GetNbMainUnknownDof();
        
        for (int m = 0; m < var_problem.nb_unknowns_vec; m++)
          {
            Bh->SetUnknownNumber(m + var_problem.nb_unknowns_scal);
            offset = m*nodl_mesh;
            Bh->SolveMass(offset, alpha, Uh, Vh);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! applies operator ShV, Prod_Uh = beta Prod_Uh + alpha ShV Uh      
  void DiscontinuousUnsteadyMassMatrix_Base::
  ApplyOperatorShVectorial(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			   const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    // glob_chrono.Start(VirtualTimer::MASS);
    
    if (beta == Real_wp(0))
      Prod_Uh.Fill(0);
    else if (beta != Real_wp(1))
      Mlt(beta, Prod_Uh);

    int offset = 0;
    if (ShVec->GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      ShVec->MltMass(0, alpha, Uh, Prod_Uh);
    else
      {
        int nodl_mesh = var_problem.GetNbMeshDof();
        if (var_problem.FormulationDG() == ElementReference_Base::HDG)
          nodl_mesh = var_comm.GetNbMainUnknownDof();

        for (int m = 0; m < var_problem.nb_unknowns_vec; m++)
          {
            ShVec->SetUnknownNumber(m + var_problem.nb_unknowns_scal);
            offset = m*nodl_mesh;
            ShVec->MltMass(offset, alpha, Uh, Prod_Uh);
          }    
      }
    
    // glob_chrono.Stop(VirtualTimer::MASS);
  }


  //! solves mass matrix due to lambda (HDG)
  void DiscontinuousUnsteadyMassMatrix_Base::SolveOperatorCh(VectReal_wp& Y)
  {
    if (ChSurf == NULL)
      {
	cout << "Ch not present" << endl;
	abort();
      }

    var_comm.AddDomains(Y);
    ChSurf->SolveMass(0, Y);
  }


  template<class Dimension>
  DiscontinuousMassMatrix_Base* DiscontinuousUnsteadyMassMatrix<Dimension>
  ::GetNewMassMatrix(int type_mat)
  {
    return new DiscontinuousMassMatrix<Dimension>(var_problem, var_time,
                                                  this->GeomMass, type_mat);
  }


  template<class Dimension>
  DiscontinuousMassMatrix_Base* DiscontinuousUnsteadyMassMatrix<Dimension>
  ::GetNewSurfaceMassMatrix(int type_mat, bool geom)
  {
    if (geom)
      return new DiscontinuousSurfaceMassMatrix<Dimension>(var_problem, var_time,
                                                           this->SurfMass, type_mat);
    else
      {
	DiscontinuousSurfaceMassMatrix<Dimension>& mass
	  = static_cast<DiscontinuousSurfaceMassMatrix<Dimension>& >(*this->SurfMass);
	
	return new DiscontinuousSurfaceMassMatrix<Dimension>(var_problem, var_time,
                                                             this->SurfMass, type_mat,
                                                             mass.GetDs(), mass.GetOffsetDs());
      }
  }

} // end namespace

#define MONTJOIE_FILE_DISCONTINUOUS_UNSTEADY_MASS_MATRIX_CXX
#endif

