#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Hyperbolic/Elastic/MontjoieTimeElastic.hxx"
#include "Elliptic/Helmholtz/MontjoieLaplace.hxx"

using namespace Montjoie;

namespace Montjoie
{
  template<class T>
  class PreconditioningDirect : public Preconditioner_Base<T>
  {
  protected:
    All_MatrixLU<T>& mat_lu;
    
  public:
    PreconditioningDirect(All_MatrixLU<T>& mat) : mat_lu(mat) {}
    
    void Solve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z)
    {
      z = r;
      mat_lu.SetPrintLevel(-1);
      mat_lu.Solve(z);
    }
    
    void TransSolve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z)
    {
      z = r;
      mat_lu.SetPrintLevel(-1);
      mat_lu.TransSolve(z);
    }
  
  };
  
}

//! base class for piezo-electric system (coupling between Laplace equation and elastodynamics)
template<class TypeEquation>
class PiezoInteraction_Base : public InputDataProblem_Base,
			      public VirtualOdeSystem<Real_wp>
{
public:
  //! object for time-domain elastic equation
  HyperbolicProblem<TimeElasticEquation<Dimension2> > time_elas;
  //! object for elastic equation
  EllipticProblem<TypeEquation>& var_elas;
  //! object for Poisson's equation
  EllipticProblem<LaplaceEquation<Dimension2> > var_laplace;
  //! if true, the time-harmonic problem is considered
  bool complex_solve, solve_laplace;

  //! tensor d of piezoelectricity
  Vector<TinySymmetricVecTensor<Real_wp, 2> > ref_media_d;

  //! references where there is a piezo-electrical material
  Vector<int> ref_piezo;
  //! surface references correspoding to the cathodes 
  Vector<int> ref_cathode;
  //! properties associated with each cathode
  VectReal_wp resistance_cathode, potentiel_cathode;
  
  //! coefficients associated with Poisson's problem
  GlobalGenericMatrix<Real_wp> nat_mat_lap;
  //! linear solver used for Poisson's problem
  All_LinearSolver laplace_solver;

  //! operator b
  Matrix<Real_wp, General, ArrayRowSparse> Bh;
  //! columns b(phi_i, :)
  Vector<Vector<Real_wp, VectSparse> > BhCol;
  //! coefficients c(phi_i, phi_i)
  VectReal_wp CoefCii;
  
  //! range of frequency
  VectReal_wp omega;

  bool use_iterative_solver;
  int nb_iter_max_eff_iterative_solver;
  
  //! output class for lambda
  ParamOutputClass sismo_lambda;
  
  //! arrays used for parallel computation
  IVect row_num, overlap_num, proc_num, MatchingProc;
  Vector<IVect> MatchingDofNumber;
  int nodl_global, nodl_local;
  
  //! constructor with an elastic problem
  PiezoInteraction_Base(EllipticProblem<TypeEquation>& var) : var_elas(var), laplace_solver(var_laplace)
  {
    complex_solve = false;
    solve_laplace = true;
    nodl_global = 0; nodl_local = 0;
    use_iterative_solver = false;
    nb_iter_max_eff_iterative_solver = 20;
  }
  
  //! the meshes are constructed (both for elasto and Poisson's problem
  void ConstructPiezoMesh(IVect& NumEltPiezo)
  {
    Mesh<Dimension2>& mesh = var_elas.mesh;    
    if (var_elas.GetRankProcMode() == 0)
      {
        // mesh for elastodynamics is constructed
        mesh.ConstructMesh(0, var_elas.mesh_data(0));
        
        int ref_max = 1;
        for (int i = 0; i < mesh.GetNbElt(); i++)
          ref_max = max(ref_max, mesh.Element(i).GetReference());
        
        // is_piezo will be true if the reference i is associated with a piezo material
        Vector<bool> is_piezo(ref_max+1);
        is_piezo.Fill(false);
        for (int i = 0; i < ref_piezo.GetM(); i++)
          if (ref_piezo(i) <= ref_max)
            is_piezo(ref_piezo(i)) = true;
        
        // we store boundary condition for Poisson's problem
        IVect laplace_boundary_condition = var_laplace.mesh.GetBoundaryCondition();

        // we copy the mesh
        var_laplace.mesh = mesh;
        for (int ref = 1; ref < laplace_boundary_condition.GetM(); ref++)
          var_laplace.mesh.SetBoundaryCondition(ref, laplace_boundary_condition(ref));
        
        // NumEltPiezo will contains all the elements associated with a piezo-electric material
        // in var_laplace.mesh, we keep only these elements, other elements are removed
        int nb_elt = 0;
        NumEltPiezo.Reallocate(mesh.GetNbElt());
        for (int i = 0; i < mesh.GetNbElt(); i++)
          {
            int ref = var_laplace.mesh.Element(i).GetReference();
            if (is_piezo(ref))
              {
                var_laplace.mesh.Element(nb_elt) = mesh.Element(i);
                NumEltPiezo(nb_elt) = i;
                nb_elt++;
              }
          }
        
        NumEltPiezo.Resize(nb_elt);
        var_laplace.mesh.ResizeElements(nb_elt);
        
        // unused vertices and edges are removed
        var_laplace.mesh.ForceCoherenceMesh();
        var_laplace.mesh.ReorientElements();
        var_laplace.mesh.FindConnectivity();
        
        // the mesh associated with piezo-electric material is written
        // then it will be read by the master processor for Poisson's problem
        var_laplace.mesh.Write("piezo.msh");
      }
  }
  
  //! constructors a MPI communicator for Poisson's problem
  //! involving only processors that contain some piezo elements
  //! NumEltPiezo : list of elements that are piezo
  bool ConstructCommunicatorLaplace(IVect& NumEltPiezo)
  {
    // it will return true if the current processor does contain piezo elements
    // otherwise it will return false
    bool solve_lap = true;
    
#ifdef SELDON_WITH_MPI
    if (var_elas.GetNbProcPerMode() > 1)
      {
        int my_key = -1; MPI_Status status;
        IVect EpartLap; int my_send_epart, nb_elt_local_laplace;
        if (var_elas.GetRankProcMode() == 0)
          {
            int nb_process = var_elas.GetNbProcPerMode();
            IVect Epart; Epart.Read("Epart.dat");
            // NumLoc(i) contains the local number of the element on the local processor
            Vector<int> NumLoc(Epart.GetM());
            Vector<int> NumElem(nb_process);
            NumElem.Zero();
            for (int i = 0; i < Epart.GetM(); i++)
              {
                int p = Epart(i);
                NumLoc(i) = NumElem(p);
                NumElem(p)++;
              }
            
            // which processors are involved in Poisson problem ?
            // we store also for each processor the number of elements
            Vector<bool> ProcUsed(nb_process);
            ProcUsed.Fill(false); NumElem.Zero();
            int nb_proc_used = 0;
            for (int i = 0; i < NumEltPiezo.GetM(); i++)
              {
                int ne = NumEltPiezo(i);
                int p = Epart(ne);
                if (!ProcUsed(p))
                  {
                    ProcUsed(p) = true;
                    nb_proc_used++;
                  }
                
                NumElem(p)++;
              }
            
            // NumProcLaplace will contain the list of processors involved in Poisson problem
            // first_proc is the first element of NumProcLaplace
            IVect NumProcLaplace(nb_proc_used), InvNumProc(nb_process);
            nb_proc_used = 0; int first_proc = -1;
            for (int i = 0; i < ProcUsed.GetM(); i++)
              {
                int key = -1;
                if (ProcUsed(i))
                  {
                    InvNumProc(i) = nb_proc_used;
                    key = nb_proc_used;
                    if (nb_proc_used == 0)
                      first_proc = i;
                    
                    NumProcLaplace(nb_proc_used++) = i;                
                  }
                
                if (i == 0)
                  my_key = key;
                else
                  MPI_Send(&key, 1, MPI_INTEGER, i, 98, var_elas.comm_group_mode);
                
                int send_epart = 0;
                if (first_proc == i)
                  send_epart = var_laplace.mesh.GetNbElt();
                
                int nb_elt_loc = NumElem(i);
                
                if (i == 0)
                  {
                    my_send_epart = send_epart;
                    nb_elt_local_laplace = nb_elt_loc;
                  }
                else
                  {
                    MPI_Send(&send_epart, 1, MPI_INTEGER, i, 99, var_elas.comm_group_mode);
                    MPI_Send(&nb_elt_loc, 1, MPI_INTEGER, i, 101, var_elas.comm_group_mode);
                  }
              }
            
            // generating Epart for var_laplace
            EpartLap.Reallocate(var_laplace.mesh.GetNbElt());
            for (int i = 0; i < NumEltPiezo.GetM(); i++)
              {
                int ne = NumEltPiezo(i);
                int p = Epart(ne);
                EpartLap(i) = InvNumProc(p);
              }                        
            
            if (first_proc != 0)
              MPI_Send(EpartLap.GetData(), EpartLap.GetM(), MPI_INTEGER, first_proc, 100, var_elas.comm_group_mode);

            // generating NumEltPiezo for each proc
            Vector<IVect> NumEltPiezoLoc(nb_process);
            for (int i = 0; i < nb_process; i++)
              NumEltPiezoLoc(i).Reallocate(NumElem(i));
            
            NumElem.Zero();
            for (int i = 0; i < NumEltPiezo.GetM(); i++)
              {
                int ne = NumEltPiezo(i);
                int p = Epart(ne);
                NumEltPiezoLoc(p)(NumElem(p)) = NumLoc(ne);
                NumElem(p)++;
              }
            
            for (int i = 0; i < nb_process; i++)
              {
                if (i == 0)
                  NumEltPiezo = NumEltPiezoLoc(0);
                else
                  {
                    if (NumElem(i) > 0)
                      MPI_Send(NumEltPiezoLoc(i).GetData(), NumElem(i), 
                               MPI_INTEGER, i, 102, var_elas.comm_group_mode);
                  }
              }
          }        
        else
          {
            MPI_Recv(&my_key, 1, MPI_INTEGER, 0, 98, var_elas.comm_group_mode, &status);
            MPI_Recv(&my_send_epart, 1, MPI_INTEGER, 0, 99, var_elas.comm_group_mode, &status);
            MPI_Recv(&nb_elt_local_laplace, 1, MPI_INTEGER, 0, 101, var_elas.comm_group_mode, &status);
            if (my_send_epart > 0)
              {
                EpartLap.Reallocate(my_send_epart);
                MPI_Recv(EpartLap.GetData(), EpartLap.GetM(),
                         MPI_INTEGER, 0, 100, var_elas.comm_group_mode, &status);

                var_laplace.SetEpartSplitting(EpartLap);
              }

            if (nb_elt_local_laplace > 0)
              {
                NumEltPiezo.Reallocate(nb_elt_local_laplace);
                MPI_Recv(NumEltPiezo.GetData(), NumEltPiezo.GetM(),
                         MPI_INTEGER, 0, 102, var_elas.comm_group_mode, &status);
              }
            else
              NumEltPiezo.Clear();            
          }
        
        // creating MPI communicator for var_laplace
        MPI_Comm comm_lap; int color = MPI_UNDEFINED;
        if (my_key >= 0)
          color = 0;
        
        MPI_Comm_split(var_elas.comm_group_mode, color, my_key, &comm_lap);    
        var_laplace.comm_group_mode = comm_lap;
        
        if (comm_lap == MPI_COMM_NULL)
          solve_lap = false;
        else
          {
            int rank_laplace; MPI_Comm_rank(comm_lap, &rank_laplace);
            if (rank_laplace == 0)
              {
                // for master processor of Poisson's problem, we
                // read the piezo mesh, and set the partioning, such that 
                // the partioning between Poisson and elastodynamic coincides
                var_laplace.mesh.Read("piezo.msh");
                var_laplace.SetEpartSplitting(EpartLap);
              }
          }
      }
#endif
    
    return solve_lap;
  }

#ifdef SELDON_WITH_MPI
  //! computes global row numbers of the coupled problem and inits the distributed matrix A
  void InitRowNumbers(DistributedMatrix_Base<Complex_wp>& A)
  {
    if (var_elas.GetNbProcPerMode() <= 1)
      return;
    
    row_num = var_elas.GetGlobalDofNumber();
    IVect row_lap;
    if (solve_laplace)
      {
        row_lap = var_laplace.GetGlobalDofNumber();
        int offset = var_elas.GetNbGlobalDof();
        for (int i = 0; i < row_lap.GetM(); i++)
          row_lap(i) += offset;
        
        row_num.PushBack(row_lap);
        
        IVect row_cathode;
        row_cathode.Reallocate(ref_cathode.GetM());
        offset += var_laplace.GetNbGlobalDof();
        for (int i = 0; i < ref_cathode.GetM(); i++)
          row_cathode(i) = offset + i;
        
        row_num.PushBack(row_cathode);
      }
    
    A.Init(row_num, overlap_num, proc_num, MatchingProc, MatchingDofNumber, var_elas.comm_group_mode);    
    nodl_global = A.GetGlobalM();
    nodl_local = A.GetNodlScalar();
  }


  void InitDistributedMatrix(DistributedMatrix_Base<Complex_wp>& A)
  {
    A.Init(nodl_global, &row_num, &overlap_num, &proc_num, nodl_local, 1,
           &MatchingProc, &MatchingDofNumber, var_elas.comm_group_mode);
  }
#endif

  //! Computes operator b
  void ConstructBh(const IVect& NumEltPiezo)
  {
    Bh.Reallocate(var_laplace.GetNbDof(), 2*var_elas.offset_dof_unknown(1));
    TinyVector<bool, 4> null_term; null_term.Fill(true);
    null_term(1) = false;
    Matrix<Real_wp> mat_elem;
    int Nvol = var_elas.offset_dof_unknown(1);    
    Vector<int> col_interac; VectReal_wp val_interac;
    for (int i = 0; i < var_laplace.mesh.GetNbElt(); i++)
      {
	int ne = NumEltPiezo(i);
	int ref = var_elas.mesh.Element(ne).GetReference();
	TinySymmetricVecTensor<Real_wp, 2> D = ref_media_d(ref);
	const ElementReference<Dimension2, 1>& Fb = var_elas.GetReferenceElementH1(ne);
	
	IVect num_ddl_laplace = var_laplace.GetDofNumberOnElement(i);
	IVect num_ddl_elas = var_elas.GetDofNumberOnElement(ne);
	int nb_dof_loc = num_ddl_elas.GetM();
	num_ddl_elas.Resize(2*nb_dof_loc);
	for (int j = 0; j < nb_dof_loc; j++)
	  num_ddl_elas(nb_dof_loc+j) = Nvol + num_ddl_elas(j);
	
	mat_elem.Reallocate(nb_dof_loc, 2*nb_dof_loc); mat_elem.Zero();
	bool affine = var_laplace.mesh.IsElementAffine(i);
	Real_wp jacob; Matrix2_2 Dgeom, dfjm1;
	if (affine)
	  {
	    Matrix2_2 Cquad0, Cquad1;
	    jacob = var_elas.Glob_jacobian(ne)(0);
	    dfjm1 = var_elas.Glob_DFjm1(ne)(0);
	    Mlt(1.0 / jacob, dfjm1);

	    GetMatrixD_dfj(0, D, dfjm1, Dgeom);
	    Cquad0 = jacob*Dgeom;
	    GetMatrixD_dfj(1, D, dfjm1, Dgeom);
	    Cquad1 = jacob*Dgeom;

	    Real_wp mass; R2 Dquad, Equad;
	    Fb.AddConstantElemMatrix(0, 0, mass, Cquad0, Dquad, Equad,
				     null_term, mat_elem);

	    Fb.AddConstantElemMatrix(0, nb_dof_loc, mass, Cquad1, Dquad, Equad,
				     null_term, mat_elem);

	  }
	else
	  {
	    VectReal_wp mass; VectR2 Dquad, Equad;
	    int Nquad = Fb.GetNbPointsQuadratureInside();
	    Vector<Matrix2_2> Cquad0(Nquad), Cquad1(Nquad);	    
	    for (int k = 0; k < Nquad; k++)
	      {
		jacob = var_elas.Glob_jacobian(ne)(k);
		dfjm1 = var_elas.Glob_DFjm1(ne)(k);
		Mlt(Fb.WeightsND(k)/jacob, dfjm1);

		GetMatrixD_dfj(0, D, dfjm1, Dgeom);
		Cquad0(k) = jacob*Dgeom;
		GetMatrixD_dfj(1, D, dfjm1, Dgeom);
		Cquad1(k) = jacob*Dgeom;
	      }

	    Fb.AddVariableElemMatrix(0, 0, mass, Cquad0, Dquad, Equad,
				     null_term, mat_elem);

	    Fb.AddVariableElemMatrix(0, nb_dof_loc, mass, Cquad1, Dquad, Equad,
				     null_term, mat_elem);
	  }

	col_interac.Reallocate(mat_elem.GetN());
	val_interac.Reallocate(mat_elem.GetN());
	col_interac.Fill(-1); val_interac.Zero();
	for (int j = 0; j < nb_dof_loc; j++)
	  {
	    int num_dof = num_ddl_laplace(j);
	    if (num_dof >= 0)
	      {
		int c_interac = 0;
		for (int k = 0; k < mat_elem.GetN(); k++)
		  if (mat_elem(j, k) != Real_wp(0))
		    {
		      col_interac(c_interac) = num_ddl_elas(k);
		      val_interac(c_interac) = mat_elem(j, k);
		      c_interac++;
		    }

		Bh.AddInteractionRow(num_dof, c_interac, col_interac, val_interac);
	      }
	  }
      }
  }
  
  //! Computes potential phi_i
  void ConstructPotentialPhiI()
  {
    // Laplacian matrix stored in mat_laplace
    FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension2> > mat_laplace(var_laplace);
    var_laplace.SetStorageFiniteElementMatrix(var_laplace.MATRIX_FREE);
    var_laplace.AddMatrixWithBC(mat_laplace, nat_mat_lap);
    mat_laplace.IgnoreDirichletDof();
    
    VectReal_wp phi_i(var_laplace.GetNbDof()), rhs_i(var_laplace.GetNbDof());
    int Ns = 2*var_elas.offset_dof_unknown(1);
    VectReal_wp val_interac(Ns);

    // loop over cathodes
    for (int i = 0; i < ref_cathode.GetM(); i++)
      {
	Vector<int> DofsCathode;
	IVect ref_cathode_i(1); ref_cathode_i(0) = ref_cathode(i);
	var_laplace.FindDofsOnReference(ref_cathode_i, DofsCathode);
	rhs_i.Zero();
        // phi_i = 1 on the specified cathode 
	for (int j = 0; j < DofsCathode.GetM(); j++)
	  rhs_i(DofsCathode(j)) = Real_wp(1);
        
        // we solve Poisson's problem to find the potential phi_i
	phi_i = rhs_i;
	laplace_solver.ComputeSolution(phi_i, nat_mat_lap);	
	
	// we compute c(phi_i, phi_i)
	rhs_i.Zero();
	mat_laplace.MltVector(phi_i, rhs_i);
        
        DistributedVector<Real_wp>* rhs_dist_i = var_laplace.AllocateDistributedVector(rhs_i);
        DistributedVector<Real_wp>* phi_dist_i = var_laplace.AllocateDistributedVector(phi_i);

	CoefCii(i) = DotProd(*rhs_dist_i, *phi_dist_i);

        var_laplace.NullifyDistributedVector(rhs_dist_i);
        var_laplace.NullifyDistributedVector(phi_dist_i);
        
	// we compute b(phi_i, v), these columns are stored in BhCol(i)
	Mlt(SeldonTrans, Bh, phi_i, val_interac);
	Copy(val_interac, BhCol(i));
      }
  }
    
  // a line of the data file is read
  void SetInputData(const string& keyword, const Vector<string>& params)
  {
    if (keyword == "ReferencePiezo")
      {
        // list of references associated with a piezo-electric material
	ref_piezo.Reallocate(params.GetM());
	for (int i = 0; i < params.GetM(); i++)
	  ref_piezo(i) = to_num<int>(params(i));
      }
    else if (keyword == "CathodePiezo")
      {
        // for each cathode the input file must contain a line:
        // CathodePiezo = ref resistance potential
	ref_cathode.PushBack(to_num<int>(params(0)));
	resistance_cathode.PushBack(to_num<Real_wp>(params(1)));
	potentiel_cathode.PushBack(to_num<Real_wp>(params(2)));
      }
    else if (keyword == "PhysicalMedia")
      {
        // PhysicalMedia = ref type_anisotropy coefs
        // in this field the user gives the value of the tensor d
	int ref = to_num<int>(params(0));
	if (params(1) == "6M2")
	  {
	    Real_wp coef = to_num<Real_wp>(params(2));
	    ref_media_d(ref)(0, 0, 0) = coef;
	    ref_media_d(ref)(0, 1, 1) = -coef;
	    ref_media_d(ref)(1, 0, 1) = -coef;
	  }
	else if (params(1) == "ANISOTROPE")
	  {
	    ref_media_d(ref)(0, 0, 0) = to_num<Real_wp>(params(2));
	    ref_media_d(ref)(0, 0, 1) = to_num<Real_wp>(params(3));
	    ref_media_d(ref)(0, 1, 1) = to_num<Real_wp>(params(4));
	    ref_media_d(ref)(1, 0, 0) = to_num<Real_wp>(params(5));
	    ref_media_d(ref)(1, 0, 1) = to_num<Real_wp>(params(6));
	    ref_media_d(ref)(1, 1, 1) = to_num<Real_wp>(params(7));
	  }
	else
	  {
	    cout << "Unknown material" << endl;
	    abort();
	  }
      }
    else if (keyword == "RangeFrequency")
      {
        Real_wp omega_min = 2*pi_wp*to_num<Real_wp>(params(0));
        Real_wp omega_max = 2*pi_wp*to_num<Real_wp>(params(1));
        int nb_omega = to_num<int>(params(2));
	Linspace(omega_min, omega_max, nb_omega, omega);
      }
    else if (keyword == "UseIterativeSolverFreq")
      {
	if (params(0) == "YES")
          {
            use_iterative_solver = true;
            if (params.GetM() > 1)
              nb_iter_max_eff_iterative_solver = to_num<int>(params(1));
          }
        else
          use_iterative_solver = false;
      }
    else if (keyword == "FileOutputLambda")
      {
        sismo_lambda.SetDiffractedFieldFile(params(0));
        sismo_lambda.SetTotalFieldFile(params(1));
      }
    else if (keyword == "ParametersOutputLambda")
      {
        sismo_lambda.SetInputData(params, true);
      }
    else if (keyword.find("Solid") == 0)
      {
        // lines beginning with Solid are dedicated to elastodynamics
        string keyword2 = keyword.substr(5);
        if (this->complex_solve)
          var_elas.SetInputData(keyword2, params);	
        else
          time_elas.SetInputData(keyword2, params);	
      }
    else if (keyword.find("Piezo") == 0)
      {
        // lines beginning with Piezo are dedicated to Poisson's problem
	string keyword2 = keyword.substr(5);
        var_laplace.SetInputData(keyword2, params);
      }    
    else
      {
        // other lines concern both problems
	var_laplace.SetInputData(keyword, params);
	if (this->complex_solve)
          var_elas.SetInputData(keyword, params);
        else
          time_elas.SetInputData(keyword, params);
      }
  }
  
  // intermediate function for computing Bh
  void GetMatrixD_dfj(int i, const TinySymmetricVecTensor<Real_wp, 2>& D,
		      const Matrix2_2& dfjm1, Matrix2_2& Dgeom)
  {
    Dgeom(0, 0) = D(0, i, 0)*dfjm1(0, 0)*dfjm1(0, 0)
      + D(1, i, 0)*dfjm1(0, 1)*dfjm1(0, 0)
      + D(0, i, 1)*dfjm1(0, 0)*dfjm1(0, 1)
      + D(1, i, 1)*dfjm1(0, 1)*dfjm1(0, 1);

    Dgeom(1, 0) = D(0, i, 0)*dfjm1(1, 0)*dfjm1(0, 0)
      + D(1, i, 0)*dfjm1(1, 1)*dfjm1(0, 0)
      + D(0, i, 1)*dfjm1(1, 0)*dfjm1(0, 1)
      + D(1, i, 1)*dfjm1(1, 1)*dfjm1(0, 1);

    Dgeom(0, 1) = D(0, i, 0)*dfjm1(0, 0)*dfjm1(1, 0)
      + D(1, i, 0)*dfjm1(0, 1)*dfjm1(1, 0)
      + D(0, i, 1)*dfjm1(0, 0)*dfjm1(1, 1)
      + D(1, i, 1)*dfjm1(0, 1)*dfjm1(1, 1);
    
    Dgeom(1, 1) = D(0, i, 0)*dfjm1(1, 0)*dfjm1(1, 0)
      + D(1, i, 0)*dfjm1(1, 1)*dfjm1(1, 0)
      + D(0, i, 1)*dfjm1(1, 0)*dfjm1(1, 1)
      + D(1, i, 1)*dfjm1(1, 1)*dfjm1(1, 1);
  }

  // main function for running a simulation
  void RunAll(const string& input_file, const string& type_element,
	      const string& type_equation)
  {    
    // the data file is read
    Vector<string> lines_data_file;
    if (input_file.compare("NONE"))
      {    
#ifdef SELDON_WITH_MPI
        ReadLinesFile(input_file, lines_data_file, MPI_COMM_WORLD);
#else
        ReadLinesFile(input_file, lines_data_file);
#endif
      }

    // the maximal number of indices
    ref_media_d.Reallocate(PhysicalConstant::nb_max_indices);
    var_elas.InitIndices(PhysicalConstant::nb_max_indices);
    var_laplace.InitIndices(PhysicalConstant::nb_max_indices);

    ReadInputFile(lines_data_file, *this);

    // first-order formulation for elastodynamics if a first-order formulation is selected
    if ((!complex_solve) && (time_elas.FirstOrderScheme()))
      var_elas.SetFirstOrderFormulation(true);

    // one frequency ?
    if (omega.GetM() <= 1)
      {
        if (omega.GetM() == 0)
          {
            omega.Reallocate(1);
            omega(0) = var_elas.GetOmega();
          }
      }
    
    // we choose the first frequency of omega
    var_elas.SetOmega(omega(0));
    
    // General mesh is constructed
    IVect NumEltPiezo;
    ConstructPiezoMesh(NumEltPiezo);

    var_elas.mesh_data.Clear();
    var_laplace.mesh_data.Clear();
    
    // mesh for elastodynamics is split into all processors
    var_elas.SaveEpartSplitting("Epart.dat");
    var_elas.ComputeMeshAndFiniteElement(type_element);
    
    // we construct the communicator for Poisson's problem
    // solve_laplace is true if the current processor does contain piezo elements
    solve_laplace = ConstructCommunicatorLaplace(NumEltPiezo);
    
    if (solve_laplace)
      {
        // mesh for Poisson's problem is split into involved processors
        var_laplace.ComputeMeshAndFiniteElement(type_element);
        if (var_laplace.GetRankProcMode() == 0)
          cout << "Laplace done" << endl;
        
        // Needed quantities to compute Laplacian matrix are computed
        var_laplace.PerformOtherInitializations();
        var_laplace.ComputeMassMatrix();
        var_laplace.ComputeQuasiPeriodicPhase();
        
        ReadInputFile(lines_data_file, laplace_solver);
      }
    
    // needed quantites to compute finite element matrix for elasticity are computed
    var_elas.PerformOtherInitializations();
    var_elas.ComputeMassMatrix();
    var_elas.ComputeQuasiPeriodicPhase();

    // we construct operator B
    if (solve_laplace)
      ConstructBh(NumEltPiezo);

    //Bh.WriteText("Bh.dat");
    // Ns : number of dofs for elastodynamics, Np : number of dofs for piezo,
    // Nb : number of dofs for cathodes (lambda_i unknowns)
    int Nvol = var_elas.offset_dof_unknown(1);
    int Ns = 2*Nvol, Np = var_laplace.GetNbDof(), Nb = ref_cathode.GetM();
    if (!solve_laplace)
      {
        Np = 0;
        Nb = 0;
      }
    
    // we find functions phi_i, and compute B^T phi_i
    nat_mat_lap.SetCoefStiffness(Real_wp(1));
    nat_mat_lap.SetCoefMass(Real_wp(0));
    nat_mat_lap.SetCoefDamping(Real_wp(0));

    if (solve_laplace)
      {
        laplace_solver.PerformFactorizationStep(nat_mat_lap);
                
        // we construct functions phi_i
        CoefCii.Reallocate(ref_cathode.GetM());
        BhCol.Reallocate(ref_cathode.GetM());
        
        ConstructPotentialPhiI();
      }
    
    if (complex_solve)
      {
	// we consider that time-harmonic simulation is asked
	if (solve_laplace)
          laplace_solver.Clear();
	
	Matrix<Real_wp, General, ArrayRowSparse> Bh_trans;
	if (solve_laplace)
          Transpose(Bh, Bh_trans);
	
        // general matrix is symmetric
	//DistributedMatrix<Complex_wp, General, ArrayRowSparse> mat_glob;
	DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse> mat_glob, A;
	mat_glob.Reallocate(Ns+Np+Nb, Ns+Np+Nb);
        
	DistributedMatrix<Complex_wp, Symmetric, RowSymSparse> Acsr;

#ifdef SELDON_WITH_MPI
        InitRowNumbers(mat_glob);
#endif        
        
	// part with Bh
        int rank_laplace = -1;
	if (solve_laplace)
          {
            MPI_Comm_rank(var_laplace.comm_group_mode, &rank_laplace);
            for (int i = 0; i < Ns; i++)
              {
                int size_row = Bh_trans.GetRowSize(i);
                mat_glob.ReallocateRow(i, size_row);
                for (int j = 0; j < size_row; j++)
                  {
                    mat_glob.Index(i, j) = Ns + Bh_trans.Index(i, j);
                    mat_glob.Value(i, j) = Bh_trans.Value(i, j);
                  }
              }
            
            if (!mat_glob.IsSymmetric())
              for (int i = 0; i < Np; i++)
                {
                  int size_row = Bh.GetRowSize(i);
                  mat_glob.ReallocateRow(Ns+i, size_row);
                  for (int j = 0; j < size_row; j++)
                    {
                      mat_glob.Index(Ns+i, j) = Bh.Index(i, j);
                      mat_glob.Value(Ns+i, j) = Bh.Value(i, j);
                    }
                }
            
            Bh.Clear();
            Bh_trans.Clear();

            // part with b(phi_i)
            for (int i = 0; i < Nb; i++)
              for (int j = 0; j < BhCol(i).GetM(); j++)
                {
                  int num = BhCol(i).Index(j);
                  mat_glob.AddInteraction(num, Np+Ns+i, BhCol(i).Value(j));
                  mat_glob.AddInteraction(Np+Ns+i, num, BhCol(i).Value(j));
                }
            
            BhCol.Clear();

            // part of Poisson equation -div( epsilon \grad phi)
            GlobalGenericMatrix<Complex_wp> nat_mat;
            nat_mat.SetCoefMass(Complex_wp(0, 0));
            nat_mat.SetCoefDamping(Complex_wp(0, 0));
            nat_mat.SetCoefStiffness(Complex_wp(-1, 0));
            var_laplace.AddMatrixWithBC(mat_glob, nat_mat, 2*Nvol, 2*Nvol);
          }

        VectComplex_wp b;
        b.Reallocate(Ns+Np+Nb);
        b.Zero();
        // var_elas.ComputeRightHandSide(b);	
        
        Matrix<Complex_wp> ValLambda;
        if (rank_laplace == 0)
          {
            ValLambda.Reallocate(omega.GetM(), Nb);
            ValLambda.Zero();
          }

        All_MatrixLU<Complex_wp> mat_lu;
        mat_lu.SetPrintLevel(4);
        
        PreconditioningDirect<Complex_wp> prec(mat_lu);
        
        All_IterativeSolver<Complex_wp> iterative_solver;
        iterative_solver.SetPrintLevel(6);
        iterative_solver.SetSolverType(iterative_solver.GMRES);
        iterative_solver.SetMaxNumberIteration(1000);
        iterative_solver.SetRestart(30);
        iterative_solver.SetStoppingCriterion(1e-6);
        iterative_solver.SetPreconditioning(prec);
        
        VectComplex_wp old_solution;
        bool refacto_matrix = true;
        for (int n_freq = 0; n_freq < omega.GetM(); n_freq++)
          {
            //cout << "calcul geom" << endl;
            var_elas.SetOmega(omega(n_freq));
            var_elas.ComputeMassMatrix();
            var_elas.ComputeQuasiPeriodicPhase();
           
            //cout << "copie matrice" << endl;
            A = mat_glob;
            if (n_freq == omega.GetM()-1)
              mat_glob.Clear();
            
            //cout << "on rajoute la partie diagonale" << endl;
            // diagonal part for cathode/anode
            if (solve_laplace)
              for (int i = 0; i < Nb; i++)
                {              
                  // part (Y_i / (i omega) - c(phi_i, phi_i)) lambda_i
                  if (rank_laplace == 0)
                    {
                      Complex_wp vloc = Real_wp(1) / (Iwp*var_elas.GetOmega()*resistance_cathode(i)) - CoefCii(i);
                      A.AddInteraction(Np+Ns+i, Np+Ns+i, vloc);
                    }
                }    
            
            //cout << "on rajoute l'elasto" << endl;
            // part of elastodynamics
            GlobalGenericMatrix<Complex_wp> nat_mat;
            var_elas.AddMatrixWithBC(A, nat_mat, 0, 0);
	
            //cout << "ini distributed" << endl;
#ifdef SELDON_WITH_MPI
            // using distributed arrays of the piezo
            InitDistributedMatrix(A);
#endif
            
            //cout << "conversion matrice" << endl;
            Copy(A, Acsr); A.Clear();
            
            //cout << "calcul second membre" << endl;
            // potential V_i modifies the right hand side (Yi / (i omega)) Vi
            b.Zero();
            if (rank_laplace == 0)
              for (int i = 0; i < Nb; i++)
                b(Np+Ns+i) = potentiel_cathode(i) / (Iwp*var_elas.GetOmega()*resistance_cathode(i));
            
            if (!refacto_matrix)
              {
                DistributedVector<Complex_wp> source(overlap_num, var_elas.comm_group_mode);
                DistributedVector<Complex_wp> solution(overlap_num, var_elas.comm_group_mode);
                
                solution.SetData(b.GetM(), b.GetData());
                source = solution;
                for (int i = 0; i < b.GetM(); i++)
                  solution(i) = old_solution(i);
                
                solution.Zero();
                //cout << "solveur iteratif" << endl;
                Acsr.SetReadyForMltAdd(true);
                if (var_elas.GetRankProcMode() != 0)
                  iterative_solver.SetPrintLevel(-1);
                
                iterative_solver.Solve(Acsr, solution, source);
                
                solution.Nullify();
                if (iterative_solver.GetNumberIteration() >= nb_iter_max_eff_iterative_solver)
                  refacto_matrix = true;
              }
            
            // finite element matrix is solved
            if (refacto_matrix)
              {
                //cout << "factorization" << endl;
                if (n_freq == 0)
                  mat_lu.PerformAnalysis(Acsr);
                
                mat_lu.PerformFactorization(Acsr);
                        
                //cout << "resolution" << endl;
                // linear system is solved
                mat_lu.Solve(b);

                if (use_iterative_solver)
                  refacto_matrix = false;
              }
            
            old_solution = b;
            if (rank_laplace == 0)
              for (int i = 0; i < Nb; i++)
                cout << "Lambda_i = " << b(Np+Ns+i) << endl;
            
            // outputs are performed
            if (omega.GetM() <= 1)
              {
                // the values of u, phi_0 and lambda_i are extracted
                VectComplex_wp x_elas(Ns), x_lap(Np);
                if (solve_laplace)
                  for (int i = Ns; i < Np+Ns; i++)
                    x_lap(i-Ns) = b(i);
                
                for (int i = 0; i < Ns; i++)
                  x_elas(i) = b(i);
                
                var_elas.WriteDatas(x_elas);
                if (solve_laplace)
                  var_laplace.WriteDatas(x_lap);
              }
            else
              {
                // only lambda is kept
                if (rank_laplace == 0)
                  for (int i = 0; i < Nb; i++)
                    ValLambda(n_freq, i) = b(Np+Ns+i);
              }
          }
        
        if (rank_laplace == 0)
          ValLambda.Write("ValLambda.dat");
      }
    else
      {
	// time-domain simulations
	time_elas.ComputeRightHandSide();
	time_elas.ComputeStiffnessMatrix();
	time_elas.ComputeMassMatrix();
	time_elas.InitTimeIterations();	
        
        sismo_lambda.InitTime(0.0);
        sismo_lambda.InitSismo(string(""), sismo_lambda.GetTotalFieldFile(),
                                sismo_lambda.GetDiffractedFieldFile(), 0, true);
                
        if (solve_laplace)
          laplace_solver.SetPrintLevel(-1);
	
	if (time_elas.FirstOrderScheme())
	  {
            // case of an explicit scheme for the first-order formulation
            // (needed with PML)
	    VectReal_wp U0(Np + Nb + var_elas.GetNbDof());
	    U0.Zero();
	    RunTimeScheme(time_elas.GetInitialTime(), time_elas.GetFinalTime(), time_elas.GetTimeStep(),
			  U0, *this, time_elas.var_time_scheme.GetTimeScheme());
            
	    return;
	  }

        // otherwise only second-order leap-frog is implemented
	if ((time_elas.var_time_scheme.GetTimeSchemeType() != TimeSchemeEnum::LEAP_FROG_ORDER2)
	    || (time_elas.var_time_scheme.GetOrder() != 2))
	  {
	    DISP(time_elas.var_time_scheme.GetTimeSchemeType());
	    DISP(time_elas.var_time_scheme.GetOrder());
	    DISP(TimeSchemeEnum::LEAP_FROG_ORDER2);
	    cout << "Time scheme not implemented" << endl;
	    abort();
	  }
        
	Real_wp dt = time_elas.GetTimeStep(), dt2 = dt*dt;
	int nb_iter = toInteger(ceil((time_elas.GetFinalTime() - time_elas.GetInitialTime())/dt));

	VectReal_wp Un_prev(Ns), Un(Ns), Un_next(Ns), KhUn(Ns), Un_diff(Ns);
	VectReal_wp BhUn(Np), x_sol(Np), KhUn_tmp(Ns);
	VectReal_wp Lambda_prev(Nb), Lambda_n(Nb), Lambda_next(Nb);
	VectReal_wp Zstat(Nb);
        for (int i = 0; i < ref_cathode.GetM(); i++)
          {
            Un.Zero();
            if (solve_laplace)
              ScatterSparseEntry(BhCol(i), Un);
            
            var_elas.AddDomains(Un);
            time_elas.SolveOperatorDh(Un);
            
            // Schur complement is computed (Zstat)
            if (solve_laplace)
              {
                Real_wp scal = DotProd(BhCol(i), Un);            
                Real_wp scal_tmp(scal); Vector<int64_t> xtmp;
                MpiAllreduce(var_laplace.comm_group_mode, &scal_tmp,
                             xtmp, &scal, 1, MPI_SUM);
                
                Zstat(i) = 1.0/resistance_cathode(i) + dt/2*scal;
              }
          }
        
	Un_prev.Zero(); Un.Zero(); Un_next.Zero(); KhUn.Zero();
	BhUn.Zero(); x_sol.Zero(); KhUn_tmp.Zero();
	Lambda_prev.Zero(); Lambda_n.Zero(); Lambda_next.Zero();
	for (int n = 0; n <= nb_iter; n++)
	  {
	    Real_wp tn = n*dt;
	    time_elas.WriteSnapshot(n, tn, Un);

            if (solve_laplace)
              {
                Mlt(Bh, Un, BhUn);
                var_laplace.ImposeNullDirichletCondition(BhUn);
		x_sol = BhUn;
                laplace_solver.ComputeSolution(x_sol, nat_mat_lap);
                Mlt(SeldonTrans, Bh, x_sol, KhUn);
              }
            else
              KhUn.Zero();
            
	    time_elas.ApplyOperatorKh(dt2, tn, Un, -dt2, KhUn);
	    Copy(KhUn, KhUn_tmp);
            var_elas.AddDomains(KhUn);
	    time_elas.SolveOperatorDh(KhUn);

            // first prediction for U^{n+1}
	    Un_next = 2.0*Un - Un_prev + KhUn;
	    var_elas.ImposeNullDirichletCondition(Un_next);	    
	    
            if (solve_laplace)
              {
                // then lambda^{n+1} is computed by using Schur complement
                Real_wp pulse; time_elas.SourceOnlyTime(tn, 1, pulse);
                VectReal_wp potentiel(ref_cathode.GetM());
                for (int i = 0; i < ref_cathode.GetM(); i++)
                  {
                    Un_diff = Un_next - Un_prev;
                    potentiel(i) = DotProd(Un_diff, BhCol(i))/(2*dt);
                  }

                if (var_laplace.GetNbProcPerMode() > 1)
                  {
                    VectReal_wp pot_tmp(potentiel);
                    Vector<int64_t> xtmp;
                    MpiAllreduce(var_laplace.comm_group_mode, pot_tmp, xtmp,
                                 potentiel, potentiel.GetM(), MPI_SUM);
                  }
                
                for (int i = 0; i < ref_cathode.GetM(); i++)
                  {
                    Real_wp Vi = pulse*potentiel_cathode(i);
                    Real_wp rhs = potentiel(i) + Vi/resistance_cathode(i);
                    
                    rhs = 2*dt*rhs + (CoefCii(i) - dt/2*Zstat(i))*Lambda_prev(i)
                      - dt *Zstat(i)*Lambda_n(i);
                    
                    Lambda_next(i) = rhs / (CoefCii(i) + dt/2*Zstat(i));
                    
                    Real_wp coef_lambda = (Lambda_next(i) + 2*Lambda_n(i) + Lambda_prev(i))/4;
                    for (int j = 0; j < BhCol(i).GetM(); j++)
                      KhUn_tmp(BhCol(i).Index(j)) -= dt2*BhCol(i).Value(j)*coef_lambda;
                  }
              }

            // correction for U^{n+1}
            var_elas.AddDomains(KhUn_tmp);
	    time_elas.SolveOperatorDh(KhUn_tmp);
	    Un_next = 2.0*Un - Un_prev + KhUn_tmp;
	    var_elas.ImposeNullDirichletCondition(Un_next);
	    
	    // next iterate
	    Un_prev = Un; Un = Un_next;
	    Lambda_prev = Lambda_n; Lambda_n = Lambda_next;
	  }
      }
  }


  void GiveIterate(int n, const Real_wp& tn, VectReal_wp& Y)
  {
    VectReal_wp Un;
    Un.SetData(var_elas.offset_dof_unknown(1)*2, Y.GetData());
    
    time_elas.WriteSnapshot(n, tn, Un);
    
    if (solve_laplace)
      {
        if (var_laplace.GetRankProcMode() == 0)
          if (sismo_lambda.SnapshotToStore(tn))
            {
              VectReal_wp val_lambda, val_lambda_diff;
              val_lambda.Reallocate(ref_cathode.GetM()+1);
              val_lambda(0) = tn;
              int offset = var_elas.GetNbDof() + var_laplace.GetNbDof();
              for (int i = 0; i < ref_cathode.GetM(); i++)
                val_lambda(1+i) = Y(offset + i);
              
              sismo_lambda.AddVectorSismo(val_lambda, val_lambda_diff);
              sismo_lambda.IncrementSnapshot();
            }
      }
    
    Un.Nullify();
  }

  void GiveFinalIterate(int n, const Real_wp& tn, VectReal_wp& Y)
  {
    GiveIterate(n, tn, Y);
  }

  // function y = g(t, x) involved in the evolution dx/dt = g(t, x)
  void EvaluateFunction(const Real_wp& tn, const VectReal_wp& x, VectReal_wp& y,
			bool invert_mass = true, bool source = true)
  {
    VectReal_wp X_elas, Phi, ProdX_elas, ProdPhi, Lambda, ProdLambda;
    VectReal_wp Un_elas, ProdUn, ProdVn;
    Real_wp* x_ptr = x.GetData(), *y_ptr = y.GetData();
    X_elas.SetData(var_elas.GetNbDof(), x_ptr);
    Un_elas.SetData(var_elas.offset_dof_unknown(1)*2, x_ptr); x_ptr += var_elas.GetNbDof();
    bool solve_laplace = (var_laplace.GetNbDof() > 0);
    if (solve_laplace)
      {
        Phi.SetData(var_laplace.GetNbDof(), x_ptr); x_ptr += var_laplace.GetNbDof();
        Lambda.SetData(ref_cathode.GetM(), x_ptr); x_ptr += ref_cathode.GetM();
      }
    
    ProdX_elas.SetData(var_elas.GetNbDof(), y_ptr);
    ProdUn.SetData(var_elas.offset_dof_unknown(1)*2, y_ptr); y_ptr += ProdUn.GetM();
    ProdVn.SetData(var_elas.GetNbDof() - ProdUn.GetM(), y_ptr);
    y_ptr += var_elas.GetNbDof() - ProdUn.GetM();
    if (solve_laplace)
      {
        ProdPhi.SetData(var_laplace.GetNbDof(), y_ptr); y_ptr += var_laplace.GetNbDof();
        ProdLambda.SetData(ref_cathode.GetM(), y_ptr); y_ptr += ref_cathode.GetM();
      }
        
    time_elas.EvaluateFunction(tn, X_elas, ProdX_elas, false, false);

    if (solve_laplace)
      MltAdd(Real_wp(-1), SeldonTrans, Bh, Phi, Real_wp(1), ProdUn);
    
    if (solve_laplace)
      {
        VectReal_wp BhUn(var_laplace.GetNbDof());
        VectReal_wp potentiel(ref_cathode.GetM());
        Mlt(Bh, Un_elas, BhUn);
        if (invert_mass)
          {
            var_laplace.ImposeNullDirichletCondition(BhUn);
	    ProdPhi = BhUn;
            laplace_solver.ComputeSolution(ProdPhi, nat_mat_lap);
          }
        
        Real_wp pulse; time_elas.SourceOnlyTime(tn, 0, pulse);
        for (int i = 0; i < ref_cathode.GetM(); i++)
          {
            for (int j = 0; j < BhCol(i).GetM(); j++)
              ProdX_elas(BhCol(i).Index(j)) -= BhCol(i).Value(j)*Lambda(i);
            
            potentiel(i) = DotProd(BhCol(i), Un_elas);
          }
        
        if (var_laplace.GetNbProcPerMode() > 1)
          {
            VectReal_wp pot_tmp(potentiel);
            Vector<int64_t> xtmp;
            MpiAllreduce(var_laplace.comm_group_mode, pot_tmp, xtmp,
                         potentiel, potentiel.GetM(), MPI_SUM);
          }
        
        for (int i = 0; i < ref_cathode.GetM(); i++)
          {
            Real_wp Vi = pulse*potentiel_cathode(i);
            Real_wp rhs = potentiel(i) + Vi/resistance_cathode(i);
            
            ProdLambda(i) = -Lambda(i)/resistance_cathode(i) + rhs;
            if (invert_mass)
              ProdLambda(i) /= CoefCii(i);
          }
      }
    
    var_elas.ImposeNullDirichletCondition(ProdUn);
	    
    if (invert_mass)
      {
	time_elas.SolveOperatorDh(ProdUn);
	time_elas.SolveOperatorBh(ProdVn);
      }
    
    X_elas.Nullify();
    Un_elas.Nullify(); ProdUn.Nullify(); ProdVn.Nullify();
    ProdX_elas.Nullify();
    
    if (solve_laplace)
      {
        Phi.Nullify(); Lambda.Nullify();
        ProdPhi.Nullify(); ProdLambda.Nullify();
      }

    //int test_input; cout << "Waiting..." << endl; cin >> test_input;
  }

  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp&, Real_wp alpha = 1.0) {}
  
};


template<class T>
class PiezoInteraction{};

//! class used to solve time-harmonic problem
template<>
class PiezoInteraction<Complex_wp> : public PiezoInteraction_Base<HarmonicElasticEquation<Dimension2> >
{
public:
  EllipticProblem<HarmonicElasticEquation<Dimension2> > var_harmonic;

  PiezoInteraction() : PiezoInteraction_Base<HarmonicElasticEquation<Dimension2> >(var_harmonic)
  {
    this->complex_solve = true;
  }
  
};

//! class used to solve time-domain problem
template<>
class PiezoInteraction<Real_wp> : public PiezoInteraction_Base<ElasticEquation<Dimension2> >
{
public:

  PiezoInteraction() : PiezoInteraction_Base<ElasticEquation<Dimension2> >(this->time_elas.var_harmonic)
  { this->complex_solve = false; }
  
};

// main function
int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  if (argc > 1)
    {

      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      
      if (type_equation == "HARMONIC_PIEZO")
        {
          PiezoInteraction<Complex_wp> var;      
          var.RunAll(file_name_data, type_element, type_equation);
        }
      else
        {
          PiezoInteraction<Real_wp> var;      
          var.RunAll(file_name_data, type_element, type_equation);
        }
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"piezo.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  
  return FinalizeMontjoie();
}
