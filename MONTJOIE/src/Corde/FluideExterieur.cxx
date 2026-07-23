#ifndef MONTJOIE_FILE_FLUIDE_EXTERIEUR_CXX

namespace Montjoie
{
  //! default constructor
  FluideExterieur::FluideExterieur() : Kh(var_vibro.var_volume.var_harmonic)
  {
    rho = 1.0;
    c2 = 1.0;
    energie = 0.0;
  }
  
  
  //! setting parameters of the class with a line of the data file
  /*!
    \param[in] keyword keyword associated with the line of the data file
    \param[in] param parameters associated with the keyword
   */
  void FluideExterieur::SetInputData(const string & keyword, const VectString & param)
  {
    var_vibro.SetInputData(keyword,param);
    if(!keyword.compare("Density"))
      {
	if (param.GetM() > 0)
	  rho = to_num<Real_wp>(param(0));
	else
	  {
	    cout << "In SetInputData of FluideExterieur" << endl;
	    cout << "Enter a value after density, for instance :" << endl;
	    cout << "Density = 2.0" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }
      }
    if(!keyword.compare("Velocity"))
      {
        if (param.GetM() > 0)
	  c2 = square(to_num<Real_wp>(param(0)));
	else
	  {
	    cout << "In SetInputData of FluideExterieur" << endl;
	    cout << "Enter a value after Velocity, for instance :" << endl;
	    cout << "Velocity = 2.0" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }
      }
  }
  
  
  //! Initialising physical indexes of plate and fluid
  void FluideExterieur::InitIndices(int nplate, int nvol)
  {
    var_vibro.var_plate.var_harmonic.InitIndices(nplate);
    var_vibro.var_volume.var_harmonic.InitIndices(nvol);
  }
  
  
  //! constructing arrays needed for time iterations
  void FluideExterieur::ConstructAcoustics()
  {
    // un alias vers var_vibro.var_volume.var_harmonic
    EllipticProblem<LaplaceEquation<Dimension3> > &
      var_short = var_vibro.var_volume.var_harmonic;

    const MeshNumbering<Dimension3>& mesh_num = var_short.GetMeshNumbering(0);
    const MeshNumbering<Dimension2>& mesh_num_plate = var_vibro.var_plate.var_harmonic.GetMeshNumbering(0);
        
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_vibro.var_volume.var_harmonic.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif

#ifdef SELDON_WITH_MPI
    chrono.SetCommunicator(var_vibro.var_volume.var_harmonic.comm_group_mode);
#endif
    
    // prend le maillage 3D et extrait le maillage de plaque avec des dll
    // qui se correspondent.
    if (rank_proc == 0)
      cout << " Computing 3-D mesh and extracting plate ..." << endl;
    
    string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
    string prefix = var_vibro.var_volume.path_reprise;
    
    if (rank_proc == 0)
      {
        if (var_vibro.var_volume.load_reprise)
          {
            // pour une reprise on utilise le Epart stocke
            if (nb_proc > 1)
              {                
		IVect epart_vol, epart_surf; 
		epart_vol.Read(prefix + "EpartVolPianoSave" + suffix);
                var_vibro.var_volume.var_harmonic.SetEpartSplitting(epart_vol);
                
		epart_surf.Read(prefix + "EpartSurfPianoSave" + suffix);
                var_vibro.var_plate.var_harmonic.SetEpartSplitting(epart_surf);
              }
          }
        else
          {
            if (var_vibro.var_volume.save_reprise)
              {
                if (nb_proc > 1)
                  {
                    // on informe qu'il faut sauver Epart
                    var_vibro.var_volume.var_harmonic
		      .SaveEpartSplitting(prefix + "EpartVolPianoSave" + suffix);
                    
                    var_vibro.var_plate.var_harmonic
		      .SaveEpartSplitting(prefix + "EpartSurfPianoSave" + suffix);
                  }
              }
          }
      }
    
    // donner rho et mu a montjoie tels que 
    // rho d2u/dt2 - div(mu \nabla u) = 0
    TinyMatrix<Real_wp, Symmetric, 3, 3> Mu;
    Mu.SetIdentity();
    Mu *= (1.0/rho);
    for (int m = 0; m < var_vibro.var_volume.var_harmonic.ref_rho.GetM(); m++)
      {
        var_short.ref_rho(m).SetConstant(1.0/(rho*c2));
        var_short.ref_mu(m).SetConstant(Mu);
      }

    chrono.Start(VirtualTimer::ALL);
    var_vibro.ComputeMeshAndFiniteElement();
        
    chrono.Stop(VirtualTimer::ALL);
    chrono.DisplayTime(VirtualTimer::ALL, "construct 2-D and 3-D meshes");
    
    if (rank_proc == 0)
      {
        if (nb_proc > 1)
          {
            cout << "Total number of degrees of freedom in the volume : "
                 << var_short.GetNbGlobalMeshDof() << endl;
            
            cout << "Total number of degrees of freedom for Uz on the surface : "
                 << var_vibro.var_plate.var_harmonic.GetNbGlobalMeshDof() << endl;
          }
        else
          {
            cout << "The mesh contains " << var_short.mesh.GetNbElt() << " elements " << endl;
            cout << "The mesh contains " << var_short.GetNbEltPML() << " elements inside the PML" << endl;
            cout << "The mesh contains " << mesh_num.GetNbDofPML() << " dofs inside the PML " << endl;
            
            cout << "Total number of degrees of freedom in the volume : "
                 << mesh_num.GetNbDof() << endl;
            
            cout << "Total number of degrees of freedom for Uz on the surface : "
                 << mesh_num_plate.GetNbDof() << endl;
          }
      }
    
    // on construit les quantités géométriques : la transfo F, DF et le jacobien 
    if (rank_proc == 0)
      cout << "Computing 3-D jacobian matrices ... " << endl;
    
    chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);    
    var_short.ComputeMassMatrix(); // son nom n'a rien a voir avec ce qu'il fait
    
    if (rank_proc == 0)
      cout << "Computing 2-D jacobian matrices ... " << endl;

    var_vibro.var_plate.var_harmonic.ComputeMassMatrix();
    
    chrono.Stop(VirtualTimer::ALL);
    chrono.DisplayTime(VirtualTimer::ALL, "compute 2-D and 3-D jacobian matrices");
    
    // pas besoin de calculer les matrices de masse de l'instationnaire, car on le fait apres
    if (var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
      {
	if (rank_proc == 0)
	  cout << "Computing 3-D mass matrices ... " << endl;
	
	chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);    	
	var_vibro.var_volume.ComputeMassMatrix();
	
	chrono.Stop(VirtualTimer::ALL);
	chrono.DisplayTime(VirtualTimer::ALL, "compute mass matrices");
      }
    
    if (rank_proc == 0)
      cout << "Initalisation before time iterations... " << endl;
    
    chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);
    var_vibro.var_volume.InitTimeIterations(); 
    
    chrono.Stop(VirtualTimer::ALL);
    chrono.DisplayTime(VirtualTimer::ALL, "initializations");    
    
    // la matrice de masse    
    if (!var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
      {
	if (rank_proc == 0)
	  cout << "Computing mass matrices for acoustics..." << endl;
	
	chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);
	
        GlobalGenericMatrix<Real_wp> nat_mat;
        FemMatrixFreeClass<Real_wp ,LaplaceEquation<Dimension3> > Khb(var_short);
        nat_mat.SetCoefMass(1.0);
        nat_mat.SetCoefStiffness(0.0); nat_mat.SetCoefDamping(0.0);
        var_short.ComputeDiagonalMatrix(Mh, Khb, nat_mat, false);
        
        Mh_unassembled = Mh;
        var_short.AddDomains(Mh); 
        
        // la matrice d'amortissement
        nat_mat.SetCoefMass(0.0); nat_mat.SetCoefStiffness(0.0); nat_mat.SetCoefDamping(1.0);
        var_short.ComputeDiagonalMatrix(Sh, Khb, nat_mat);
        
	ContinuousMassMatrixAcoustic<Dimension3>&
	  Glob_mat_Dh = static_cast<ContinuousMassMatrixAcoustic<Dimension3>& >(*var_vibro.var_volume.Glob_mat_Dh);
        
        Glob_mat_Dh.InitDiagonalMass();
	Glob_mat_Dh.SetDiagonalDh(Mh);
        Glob_mat_Dh.SetDiagonalSh(Sh);
        
        VectReal_wp invMhPlusDtSh(Mh);
        for (int i = 0; i < Mh.GetM(); i++)
          invMhPlusDtSh(i) = 1.0/(Mh(i) + 0.5*var_vibro.var_volume.GetTimeStep()*Sh(i));
        
        Glob_mat_Dh.SetInverseDiagonalDhPlusdtSh(invMhPlusDtSh);
        
        // la matrice de rigidité
        nat_mat.SetCoefMass(0.0); nat_mat.SetCoefStiffness(1.0); nat_mat.SetCoefDamping(0.0);
        var_short.AddMatrixWithBC(static_cast<FemMatrixFreeClass_Base<Real_wp>& >(Kh), nat_mat);
        
	chrono.Stop(VirtualTimer::ALL);
	chrono.DisplayTime(VirtualTimer::ALL, "compute acoustic matrices");
      }
    else
      var_vibro.var_volume.ComputeStiffnessMatrix();

    if (rank_proc == 0)
      cout << "Initialisation of vibroacoustic stuff..." << endl;

    chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);    
    var_vibro.InitTimeIterations(); // Pour avoir MhSurf assemble sur var_vibro.MassVolSurf    
    var_vibro.var_plate.InitTimeIterations(); 
    
    chrono.Stop(VirtualTimer::ALL);
    chrono.DisplayTime(VirtualTimer::ALL, "vibroacoustic initializations");    
    
    //Mh.WriteText("MhVol.dat");
    //Sh.WriteText("ShVol.dat");
    //Kh.WriteText("KhVol.dat");
    int Nvol = mesh_num.GetNbDof() + 2*mesh_num.GetNbDofPML(); 
    if (!var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
      {
        InvMhPlusDtSh.Reallocate(Nvol); 
        for(int i = 0; i < Nvol; i++)
          InvMhPlusDtSh(i) = 1.0/( Mh(i) + 0.5*var_vibro.var_volume.GetTimeStep()*Sh(i));
      }
    
    // calcul de la matrice Cbarre
    if (rank_proc == 0)
      cout << "Computing Ch..." << endl;
    
    chrono.Reset(VirtualTimer::ALL); chrono.Start(VirtualTimer::ALL);    
    ComputeCbarre();

    chrono.Stop(VirtualTimer::ALL);
    chrono.DisplayTime(VirtualTimer::ALL, "compute Ch");    

    // allocation des vecteurs Pn et Pnm1
    Pn.Reallocate(Nvol); Pnm1.Reallocate(Nvol);
    KhPn.Reallocate(Nvol);
    Pn.Fill(0.0); Pnm1.Fill(0.0); 
    KhPn.Fill(0.0);
    int Nvec = 0;
    if (var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
      {
        Nvec = var_vibro.var_volume.GetNbVectorialUnknowns();
        Vh_n.Reallocate(Nvec);
        Vh_nm1.Reallocate(Nvec);
        BhVh_np1.Reallocate(Nvec);
        Vh_np1.Reallocate(Nvec);
        Vh_n.Fill(0);
        Vh_nm1.Fill(0);
        Vh_np1.Fill(0);
        BhVh_np1.Fill(0);
      }
    else
      {
        KhPnm1_unassembled.Reallocate(Nvol); 
        KhPn_unassembled.Reallocate(Nvol); 
        KhPnm1_unassembled.Fill(0.0);  KhPn_unassembled.Fill(0.0); 
      }
    
    // on remet a zero le chronometre
    chrono.Reset(VirtualTimer::ALL);

    if (var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
      {
        IsDofOnDomain.Reallocate(Nvol);
        IsDofVOnDomain.Reallocate(Nvec);
        IsDofOnDomain.Fill(false);
        IsDofVOnDomain.Fill(false);
        for (int i = 0; i < var_short.mesh.GetNbElt(); i++)
          if (!var_short.InsidePML(i))
            {
              for (int j = 0; j < mesh_num.Element(i).GetNbDof(); j++)
                IsDofOnDomain(mesh_num.Element(i).GetNumberDof(j)) = true;
              
              for (int j = var_short.GetOffsetDofV(i); j < var_short.GetOffsetDofV(i+1); j++)
                IsDofVOnDomain(j) = true;
            }
        
      }

  }
  
  
  //! computes the coupling matrix Ch
  void FluideExterieur::ComputeCbarre()
  {    
    Matrix<Real_wp, General, ArrayRowSparse> Ch;
    var_vibro.ComputeCbarre(Ch);
    Copy(Ch, Cbarre);
  }
  
  
  //! returns the energy of the fluid
  Real_wp FluideExterieur::GetEnergy()
  {
    return(energie);
  }
  
  
  //! Computation of P^(n+1/2) and V^(n+1) (for order 1)  P^(n+1) (for order 2)
  /*!
    \param[in] Lambda  P Lambda^n+1/2
    \param[in] temps time t^n
    \param[out] ChLambda_Minus Ch^T P Lambda^{n+1/2} for dofs on Gamma-
    \param[out] ChLambda_Plus Ch^T P Lambda^{n+1/2} for dofs on Gamma+
   */
  void FluideExterieur::Advance(const VectReal_wp & Lambda, const Real_wp& temps,
                                VectReal_wp& ChLambda_Minus, VectReal_wp& ChLambda_Plus)
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_vibro.var_volume.var_harmonic.comm_group_mode;
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    bool maitre = (rank_proc == 0);
    Real_wp dt = var_vibro.var_volume.GetTimeStep();
    const MeshNumbering<Dimension3>& mesh_num = var_vibro.var_volume.var_harmonic.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    
    VectReal_wp FP; 
    chrono.Start(VirtualTimer::FLUX);
    // FP = Ch^T P Lambda^{n+1/2}
    FP.Reallocate(var_vibro.nodl_acous);
    FP.Fill(0);
    var_vibro.MltAddChTranspose(1.0, Lambda, FP);
    chrono.Stop(VirtualTimer::FLUX); chrono.Start(VirtualTimer::COMM);
    var_vibro.AssembleAcoustic(FP);
    
    // sending surface values to volume dofs
    var_vibro.ExchangeSurfaceToVolume(FP, ChLambda_Minus, ChLambda_Plus);

    chrono.Stop(VirtualTimer::COMM);
    chrono.Start(VirtualTimer::SCHEME);
    // adding -Ch^T P Lambda^{n+1/2} to KhPn
    for (int i = 0; i < var_vibro.ddl_vol_minus.GetM(); i++)
      KhPn(var_vibro.ddl_vol_minus(i)) -= ChLambda_Minus(i);
    
    for (int i = 0; i < var_vibro.ddl_vol_plus.GetM(); i++)
      KhPn(var_vibro.ddl_vol_plus(i)) += ChLambda_Plus(i);
    
    if (maitre)
      {    
        if (abs(temps - var_vibro.var_volume.GetFinalTimeSource()) <= dt)
          cout << "La source est sur le point de s'arreter " << endl;        
      }
    
    energie = 0.0;    
    chrono.Stop(VirtualTimer::SCHEME);
    
    if (var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
      {
        chrono.Start(VirtualTimer::MASS);
        // on inverse par (Mh + dt/2 Sh)
        var_vibro.var_volume.Glob_mat_Dh->SolveOperatorDhPlusdtSh(KhPn);
        
        // on obtient ainsi P^n+1/2, on stocke P^n-1/2 aussi
        Copy(Pn, Pnm1);
        Copy(KhPn, Pn);
        chrono.Stop(VirtualTimer::MASS); chrono.Start(VirtualTimer::STIFFNESS);
        
        // on calcule V^n+1 = (Bh + dt/2 Sh^V)^-1 [  (Bh - dt/2 Sh^V) V^n  - dt R^V P^{n+1/2} ]
        var_vibro.var_volume.ApplyOperatorRhVectorial(-dt, temps+0.5*dt, Pn, 0.0, Vh_np1);
	chrono.Stop(VirtualTimer::STIFFNESS); chrono.Start(VirtualTimer::MASS);
        
        const VectReal_wp& BhMinusDtSh = var_vibro.var_volume.Glob_mat_Dh->GetDiagonalBhMinusdtSh();
        const VectReal_wp& invBhPlusDtSh = var_vibro.var_volume.Glob_mat_Dh->GetInverseDiagonalBhPlusdtSh();
        Real_wp tmp, prod_scal = 0.0;
        for (int i = 0; i < Vh_np1.GetM(); i++)
          {
            tmp = invBhPlusDtSh(i)*(Vh_np1(i) + BhMinusDtSh(i)*Vh_n(i));
            if (IsDofVOnDomain(i))
              prod_scal += BhMinusDtSh(i)*Vh_n(i)*tmp;
            
            Vh_n(i) = tmp;
          }
        
        //var_vibro.var_volume.ApplyOperatorBhMinusdtSh(1.0, temps+0.5*dt, Vh_n, 1.0, Vh_np1);
	//var_vibro.var_volume.SolveOperatorBhPlusdtSh(Vh_np1);

        chrono.Stop(VirtualTimer::MASS); chrono.Start(VirtualTimer::EXTRAPOL);

        // calcul de l'energie
        var_vibro.var_volume.ApplyOperatorDh(1.0, temps, Pn, 0.0, KhPn);
        const VectReal_wp& Dh = var_vibro.var_volume.Glob_mat_Dh->GetDiagonalDh();
        for (int i = 0; i < Pn.GetM(); i++)
          if (IsDofOnDomain(i))
            prod_scal += Dh(i)*Pn(i)*Pn(i);
        
        energie = 0.5*prod_scal;
        
        chrono.Stop(VirtualTimer::EXTRAPOL);
      }
    else
      {
        // on inverse masse + dt/2 condition absorbante
        Real_wp tmp;
        
        // on stocke (Pn-Pnm1) et (Pn+Pnm1) pour l'énergie
        Real_wp SommeP, DiffP, KhSommeP, KhDiffP;

        chrono.Start(VirtualTimer::SCHEME);
        
        for (int i = 0; i < Nvol; i++)
          {	
            // attention KhPn = - dt Kh P^n
            tmp = KhPn(i)*dt + Mh(i)*(2*Pn(i)-Pnm1(i)) + dt/2.0*Sh(i)*Pnm1(i);
            
            SommeP = (Pn(i)+Pnm1(i));
            DiffP = Pn(i) - Pnm1(i);
            KhSommeP = KhPn_unassembled(i) + KhPnm1_unassembled(i);
            KhDiffP = KhPn_unassembled(i) - KhPnm1_unassembled(i);
            tmp *= InvMhPlusDtSh(i);
            Pnm1(i) = Pn(i);
            Pn(i) = tmp;
            
            // profiter de la boucle pour calculer l'énergie de l'air
            energie += ( square(DiffP)*Mh_unassembled(i)/dt - dt*DiffP*KhDiffP*0.25
                         + SommeP*KhSommeP*0.25*dt ) *0.5/dt;
          }
        
        chrono.Stop(VirtualTimer::SCHEME);
      }
    
    var_vibro.var_volume.var_harmonic.ImposeNullDirichletCondition(Pn);
    var_vibro.var_volume.var_harmonic.ImposeNullDirichletCondition(Pnm1);
    
#ifdef SELDON_WITH_MPI
    Real_wp energie_loc = energie; Vector<int64_t> xtmp;
    MpiReduce(comm, &energie_loc, xtmp, &energie, 1, MPI_SUM, 0); 
#endif  
  }
  
  
  //! Computation of Gvol = P^{n-1/2} - dt/2 Mh^-1 (Rh^S V^n  - F^n - Ch^T P Lambda^{n-1/2} )
  //!   in assembled form
  /*!
    \param[in] invMhSurfMinus Mh^-1 on the surface Gamma^-
    \param[in] invMhSurfPlus Mh^-1 on the surface Gamma^+
    \param[in] GvolSurf_minus values of Gvol on degrees of freedom of Gamma^-
    \param[in] GvolSurf_plus values of Gvol on degrees of freedom of Gamma^+
    \param[in] temps time used for the computation of right hand side F^n
    For second-order formulation, Gvol = (P^n - P^{n-1})/ dt
    - dt/2 (Kh P^n - F^n - Ch^T P Lambda^{n-1/2})
   */
  void FluideExterieur
  ::CalculeGvol(const VectReal_wp& invMhSurfMinus, const VectReal_wp& invMhSurfPlus,
                const VectReal_wp& ChPu_Lambda_Minus, const VectReal_wp& ChPu_Lambda_Plus,
		VectReal_wp & GvolSurf_minus, VectReal_wp& GvolSurf_plus, Real_wp temps)
  {   
#ifdef SELDON_WITH_MPI
    //MPI_Comm& comm = var_vibro.var_volume.var_harmonic.comm_group_mode;
    //int nb_proc; MPI_Comm_size(comm, &nb_proc);
    //int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    //int nb_proc(1), rank_proc(0);
#endif

    Real_wp dt = var_vibro.var_volume.GetTimeStep();
    Real_wp invdt = 1.0/dt;
    chrono.Start(VirtualTimer::STIFFNESS);
    if (var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
      {
        // calcul de KhPn = Rh^S V^n
        var_vibro.var_volume.ApplyOperatorRhScalar(1.0, temps, Vh_n, 0.0, KhPn);
      }
    else
      {
        // calcul de KhPn = Kh P^n et on stocke le precedent
        Copy(KhPn_unassembled, KhPnm1_unassembled);
        Kh.MltVector(SeldonNoTrans, Pn, KhPn_unassembled, false);
        Copy(KhPn_unassembled, KhPn);
      }
    chrono.Stop(VirtualTimer::STIFFNESS); chrono.Start(VirtualTimer::SCHEME);
    
    if (temps < var_vibro.var_volume.GetFinalTimeSource())
      {
        // on retranche les termes de source
        if (var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
          var_vibro.var_volume.AddPrimitiveSourceAtTime(-1.0, temps, 0, KhPn);
        else
          var_vibro.var_volume.AddPrimitiveSourceAtTime(-1.0, temps, 0, KhPn);
      }
    
    chrono.Stop(VirtualTimer::SCHEME); chrono.Start(VirtualTimer::MASS);
    // on multiplie par -dt et on rajoute (Mh - dt/2 Sh) P^{n-/1/2} pour l'ordre 1
    Mlt(-dt, KhPn);
    if (var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
      var_vibro.var_volume.ApplyOperatorDhMinusdtSh(1.0, temps, Pn, 1.0, KhPn);

    chrono.Stop(VirtualTimer::MASS); chrono.Start(VirtualTimer::COMM);
    // on assemble KhPn entre les differents processeurs
    var_vibro.var_volume.Assemble(KhPn);
    
    // string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
    
    // on rajoute + Ch^T P Lambda^{n-1/2}
    for (int i = 0; i < var_vibro.ddl_vol_minus.GetM(); i++)
      KhPn(var_vibro.ddl_vol_minus(i)) += ChPu_Lambda_Minus(i);
    
    for (int i = 0; i < var_vibro.ddl_vol_plus.GetM(); i++)
      KhPn(var_vibro.ddl_vol_plus(i)) -= ChPu_Lambda_Plus(i);

    // on recupere Kh Pn+ et Kh Pn- sur la surface
    VectReal_wp KhPn_minus(var_vibro.ddl_vol_minus.GetM()); 
    VectReal_wp KhPn_plus(var_vibro.ddl_vol_plus.GetM());
    
    for(int i = 0; i < var_vibro.ddl_vol_minus.GetM(); i++)
      KhPn_minus(i) = KhPn(var_vibro.ddl_vol_minus(i));
    
    for(int i = 0; i < var_vibro.ddl_vol_plus.GetM(); i++)
      KhPn_plus(i) = KhPn(var_vibro.ddl_vol_plus(i));
    
    VectReal_wp KhPnSurf_minus(var_vibro.nodl_acous), KhPnSurf_plus(var_vibro.nodl_acous);
    KhPnSurf_minus.Fill(0); KhPnSurf_plus.Fill(0);
    var_vibro.ExchangeVolumeToSurface(KhPn_minus, KhPn_plus, KhPnSurf_minus, KhPnSurf_plus);
    
    // et on fait pareil avec (P^n - P^(n-1)) / Delta t ou P^n-1/2
    VectReal_wp Gvol_minus(var_vibro.ddl_vol_minus.GetM()); 
    VectReal_wp Gvol_plus(var_vibro.ddl_vol_plus.GetM());
    if (var_vibro.var_volume.var_harmonic.FirstOrderFormulation())
      {
        // pour l'ordre 1, c'est 1/2 P^{n-1/2} au lieu de (P^n - P^(n-1)) / Delta t
        for(int i = 0; i < var_vibro.ddl_vol_minus.GetM(); i++)
          Gvol_minus(i) = 0.5*Pn(var_vibro.ddl_vol_minus(i));
        
        for(int i = 0; i < var_vibro.ddl_vol_plus.GetM(); i++)
          Gvol_plus(i) = 0.5*Pn(var_vibro.ddl_vol_plus(i));
      }
    else
      {
        for(int i = 0; i < var_vibro.ddl_vol_minus.GetM(); i++)
          Gvol_minus(i)
            = (Pn(var_vibro.ddl_vol_minus(i)) - Pnm1(var_vibro.ddl_vol_minus(i)))*invdt;
        
        for(int i = 0; i < var_vibro.ddl_vol_plus.GetM(); i++)
          Gvol_plus(i) = (Pn(var_vibro.ddl_vol_plus(i)) - Pnm1(var_vibro.ddl_vol_plus(i)))*invdt;
      }

    GvolSurf_minus.Reallocate(var_vibro.nodl_acous);
    GvolSurf_plus.Reallocate(var_vibro.nodl_acous);
    GvolSurf_minus.Fill(0); GvolSurf_plus.Fill(0);
    var_vibro.ExchangeVolumeToSurface(Gvol_minus, Gvol_plus, GvolSurf_minus, GvolSurf_plus);
    
    chrono.Stop(VirtualTimer::COMM);     chrono.Start(VirtualTimer::SCHEME);
    
    // calcul de Gvol = (P^n - P^(n-1))/dt - 0.5 dt Mh^-1 (Kh P^n - Ch^T P Lambda^{n-1/2})
    // (pour l'ordre 2)
    //  et Gvol = P^n-1/2 - 0.5 dt Mh^-1 (Rh^S V^n - F^n - Ch^T P Lambda^{n-1/2})  (pour l'ordre 1)
    if (GvolSurf_minus.GetM() != var_vibro.nodl_acous)
      {
        cout << "size incorrect" << endl;
        abort();
      }

    if (GvolSurf_plus.GetM() != var_vibro.nodl_acous)
      {
        cout << "size incorrect" << endl;
        abort();
      }
    
    for (int i = 0; i < var_vibro.nodl_acous; i++)
      {
	Real_wp coef_minus = 0.5*invMhSurfMinus(i);
	Real_wp coef_plus = 0.5*invMhSurfPlus(i);
        GvolSurf_minus(i) += coef_minus*KhPnSurf_minus(i);
	GvolSurf_plus(i) += coef_plus*KhPnSurf_plus(i);
      }

    chrono.Stop(VirtualTimer::SCHEME);
  }
  
}  


#define MONTJOIE_FILE_FLUIDE_EXTERIEUR_CXX
#endif
