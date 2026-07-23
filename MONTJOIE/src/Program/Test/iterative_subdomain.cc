#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_HP_H1

#define SELDON_WITH_PRECONDITIONING

#include "Montjoie.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

#include "Quadrature/SingularIntegration2D.cxx"
#include "Quadrature/SingularIntegration3D.cxx"

using namespace Montjoie;

int rank_proc(0), nb_proc(1);

template<class T, int p>
void WriteNoBrackets(const Vector<TinyVector<T, p> >& P, const string& nom)
{
  ofstream file_out(nom.data());
  file_out.precision(cout.precision());
  for (int i = 0; i < P.GetM(); i++)
    {
      PrintNoBrackets(file_out, P(i));
      file_out << '\n';
    }
  
  file_out.close();
}

//! Parameters for the resolution of Helmholtz equation with domain decomposition method
/*!
  The considered Helmholtz equation is equal to :
  - omega^2 rho_i u - div( mu_i grad u) = f
*/
class SubdomainSolverParameter : public EuclidianDistanceClass_Base<Dimension2>,
				 public EuclidianDistanceClass_Base<Dimension3>
{
public :
  //! distances used for the cut-off function
  /*!
    if r < dist_tronc_min, the cut-off function is equal to one
    if r > dist_tronc, the cut-off function is equal to zero
  */
  Real_wp dist_tronc, dist_tronc_min;
  //! type of transmission conditions (operator T is different)
  int type_operator;
  //! available choices
  /*!
    the boundary condition is equal to :
    mu_i du/dn - i omega sqrt(rho_i mu_i) z u - i z c mu_i T u = ...
    Despres solver : T = 0 and z = 1
    Operator solver : 
    T(u, v) = \int_\Gamma \int_\Gamma |x-y|^\beta \curl_\Gamma u(x) \cdot \curl_\Gamma v(y) dx dy
    
    Operator_grad solver :
    T(u, v) = \int_\Gamma \int_\Gamma |x-y|^\beta \grad_\Gamma u(x) \cdot \grad_\Gamma v(y) dx dy
    
    \beta is equal to 1/2 in 2-D, -1/2 in 3-D
  */
  enum {DESPRES, OPERATOR, OPERATOR_GRAD, OPERATOR_DIFF, OPERATOR_T, OPERATOR_LOCAL};
  //! coefficients of the local operator
  Complex_wp CoefC0, CoefD0; VectComplex_wp CoefAl, CoefBl;
  Real_wp epsilon_pade;
  //! coefficients z for the two interfaces
  /*!
    We assume here that each processor possesses an interface with one or two processors,
    no more processors can be involved. The coefficients z are stored for the two potential interfaces
  */
  TinyVector<Complex_wp, 2> coef_z;
  //! processor number for the two interfaces
  /*!
    For the first interface, ref_neighbor_domain(0) is the number of the processor
    sharing the interface with the current processor
    For the second interface, the number is stored in ref_neighbor_domain(1)
  */
  TinyVector<int, 2> ref_neighbor_domain;  
  //! coefficients c for the two interfaces
  TinyVector<Real_wp, 2> coef_c;
  //! relaxation parameter used for Jacobi or Gauss-Seidel algorithm
  Real_wp omega_relax;
  //! if true, the iteration matrix involved in Jacobi solver is computed and stored
  bool compute_jacobi_iteration_matrix;
  //! if true, the unknowns g are computed on quadrature points instead of degrees of freedom
  bool compute_g_on_quadrature_points;
  //! if true, the interfaces are separated in operator T
  bool separated_interfaces;
  //! type of iterative solver used to obtain the solution of the linear system
  int type_resolution;
  //! available iterative solvers
  /*!
    JACOBI_SOLVER : Jacobi algorithm on the linear system with unknowns g
    GAUSS_SEIDEL_SOLVER : Gauss-Seidel algorithm on the linear system with unknowns g
    GMRES_SOLVER : Gmres solver on the linear system with unknowns g
    PRECOND_SOLVER : preconditioned COCG solver on the original linear system    
  */
  enum{JACOBI_SOLVER, GAUSS_SEIDEL_SOLVER, GMRES_SOLVER, PRECOND_SOLVER};
  //! if true the residual is computed by comparison with the exact solution
  //! if false, the residual is computed as || A x - b|| where A is the original matrix
  /*!
    This distinction is made for Jacobi or Gauss-Seidel solvers.
    For other solvers (Gmres, preconditioned), the residual is directly provided by
    the iterative solver
  */
  bool compare_solution;  

  // for periodic condition
  Real_wp L_periodic; bool x_periodic;
  Real_wp omega;
  
  SubdomainSolverParameter()
  {
    dist_tronc = 100.0;
    dist_tronc_min = 25.0;
    type_operator = OPERATOR;
    CoefC0 = 0; CoefD0 = 0;
    coef_z.Fill(1.0);
    ref_neighbor_domain(0) = rank_proc-1;
    ref_neighbor_domain(1) = rank_proc+1;
    coef_c.Fill(0.0);
    omega_relax = 0.5;
    compute_jacobi_iteration_matrix = false;
    compute_g_on_quadrature_points = false;
    separated_interfaces = true;
    type_resolution = JACOBI_SOLVER;
    compare_solution = false;
    
    L_periodic = 1.0; x_periodic = false;
  }

  void SetOmega(const Real_wp& w) { omega = w;  }
  
  bool LocalOperator() const
  {
    if (type_operator == DESPRES)
      return true;
    
    if (type_operator == OPERATOR_LOCAL)
      return true;
    
    return false;
  }
  
  void SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (keyword == "SubdomainSolver")
      {
        if ( (param(0) == "DESPRES") || (param(0) == "OPERATOR") || (param(0) == "OPERATOR_GRAD")
	     || (param(0) == "OPERATOR_DIFF") || (param(0) == "OPERATOR_T") || (param(0) == "OPERATOR_LOCAL"))
          {
            // type d'operateur
            if (param(0) == "DESPRES")
              type_operator = DESPRES;
	    else if (param(0) == "OPERATOR_GRAD")
	      type_operator = OPERATOR_GRAD;
	    else if (param(0) == "OPERATOR_DIFF")
	      type_operator = OPERATOR_DIFF;
	    else if (param(0) == "OPERATOR_T")
	      type_operator = OPERATOR_T;
            else if (param(0) == "OPERATOR_LOCAL")
              {
                if (param.GetM() <= 2)
                  {
                    cout << "Entrez plus de parametres dans SubdomainSolver" << endl;
                    cout << "Par exemple SubdomainSolver = OPERATOR_LOCAL PADE teta n epsilon omega" << endl;
                    cout << "Parametres actuels : " << endl << param << endl;
                    abort();
                  }
                
		type_operator = OPERATOR_LOCAL;
		int nb = 1;
		if (param(nb) == "PADE")
		  {
		    if (param.GetM() <= 5)
		      {
			cout << "Entrez plus de parametres dans SubdomainSolver" << endl;
			cout << "Par exemple SubdomainSolver = OPERATOR_LOCAL PADE teta n epsilon omega" << endl;
			cout << "Parametres actuels : " << endl << param << endl;
			abort();
		      }
		    
		    nb++;
		    Real_wp alpha = to_num<Real_wp>(param(nb++))*pi_wp/180;
		    int n = to_num<int>(param(nb++));
		    epsilon_pade = to_num<Real_wp>(param(nb++));
		    CoefD0 = 0.0;
		    ComputePadeCoefficientsSqrt(alpha, n, CoefC0, CoefAl, CoefBl);
		  }
                else if (param(nb) == "EMDA")
                  {
                    CoefC0 = Complex_wp(1, 0.5); nb++;
                    if (param.GetM() > 3)
                      {
                        CoefC0 = Complex_wp(1, to_num<Real_wp>(param(nb)));
                        nb++;
                      }
                  }
                else if (param(nb) == "OO2")
                  {
                    if (param.GetM() <= 5)
		      {
			cout << "Entrez plus de parametres dans SubdomainSolver" << endl;
			cout << "Par exemple SubdomainSolver = OPERATOR_LOCAL OO2 coef_k kmax dk omega" << endl;
			cout << "Parametres actuels : " << endl << param << endl;
			abort();
		      }
                    
                    nb++;
                    // coef_k such that k = omega coef_k on the interface 
                    //    therefore is should contain \sqrt{rho/mu} on the interface
                    // kmax = pi/h where h is the mesh size
                    // dk = pi/L where L is the size of the domain (for example radius of interface)
                    Real_wp coef_k = to_num<Real_wp>(param(nb++)); 
                    Real_wp kmax = to_num<Real_wp>(param(nb++)); 
                    Real_wp dk = to_num<Real_wp>(param(nb++));
                    Complex_wp k = omega*coef_k, kmin(0);
                    Complex_wp km = k - dk, kp = k + dk;
                    Complex_wp alpha_opt = Iwp*pow((k*k - km*km)*(k*k - kmin*kmin), 0.25);
                    Complex_wp beta_opt = pow((kp*kp-k*k)*(kmax*kmax-k*k), 0.25);
                    CoefC0 = (alpha_opt*beta_opt - k*k) / (alpha_opt+beta_opt);
                    CoefD0 = 1.0/(alpha_opt+beta_opt);
                    CoefC0 /= Iwp*k;
                    CoefD0 /= Iwp/k;                    
                  }
                else if (param(nb) == "ABC2")
                  {
                    nb++;
                    CoefC0 = 1.0;
                    CoefD0 = 0.5;
                  }

		omega_relax = to_num<Real_wp>(param(nb)); nb++;
                DISP(CoefC0); DISP(CoefD0);
		return;
	      }
            else
              type_operator = OPERATOR;
            
            if (param.GetM() <= nb_proc)
              {
                cout << "Entrez plus de parametres dans SubdomainSolver" << endl;
                cout << "Par exemple SubdomainSolver = OPERATOR z0 z1 ... zn c0 .. cn omega" << endl;
                cout << "Parametres actuels : " << endl << param << endl;
                abort();
              }
            
            // on lit les impedances z
            if (rank_proc > 0)
              coef_z(0) = conj(to_num<Complex_wp>(param(rank_proc)));
            
            if (rank_proc < nb_proc-1)
              coef_z(1) = to_num<Complex_wp>(param(1+rank_proc));
            
            int nb = nb_proc;
            if (!LocalOperator())
              {
                if (param.GetM() < 2*nb_proc)
                  {
                    cout << "Entrez plus de parametres dans SubdomainSolver" << endl;
                    cout << "Par exemple SubdomainSolver = OPERATOR z0 z1 ... zn c0 .. cn omega" << endl;
                    cout << "Parametres actuels : " << endl << param << endl;
                    abort();
                  }
                
                if (rank_proc > 0)
                  coef_c(0) = to_num<Real_wp>(param(nb+rank_proc-1));
            
                if (rank_proc < nb_proc-1)
                  coef_c(1) = to_num<Real_wp>(param(nb+rank_proc));
                
                nb += nb_proc-1;
              }
            
            omega_relax = to_num<Real_wp>(param(nb)); nb++;
          }
      }
    else if (keyword == "NeighboringProcessors")
      {
	if (param.GetM() <= 2)
	  {
	    cout << "Entrez plus de parametres dans NeighboringProcessors" << endl;
	    cout << "Par exemple NeighboringProcessors = p proc0 proc1" << endl;
	    cout << "Parametres actuels : " << endl << param << endl;
	    abort();
	  }	
	
	int p = to_num<int>(param(0));
	int p0 = to_num<int>(param(1));
	int p1 = to_num<int>(param(2));
	if (p == rank_proc)
	  {
	    ref_neighbor_domain(0) = p0;
	    ref_neighbor_domain(1) = p1;
	  }
      }
    else if (keyword == "CoefficientOperatorTransmission")
      {
	if (param.GetM() <= 3)
	  {
	    cout << "Entrez plus de parametres dans CoefficientOperatorTransmission" << endl;
	    cout << "Par exemple CoefficientOperatorTransmission = proc0 proc1 z c" << endl;
	    cout << "Parametres actuels : " << endl << param << endl;
	    abort();
	  }
	
	int p0 = to_num<int>(param(0));
	int p1 = to_num<int>(param(1));
	Complex_wp z = to_num<Complex_wp>(param(2));
	Real_wp c = to_num<Real_wp>(param(3));
	if ((ref_neighbor_domain(0) == p0) && (rank_proc == p1))
	  {
	    coef_c(0) = c;
	    coef_z(0) = conj(z);
	  }
	else if ((ref_neighbor_domain(0) == p1) && (rank_proc == p0))
	  {
	    coef_c(0) = c;
	    coef_z(0) = z;
	  }

	if ((ref_neighbor_domain(1) == p0) && (rank_proc == p1))
	  {
	    coef_c(1) = c;
	    coef_z(1) = conj(z);
	  }
	else if ((ref_neighbor_domain(1) == p1) && (rank_proc == p0))
	  {
	    coef_c(1) = c;
	    coef_z(1) = z;
	  }	
      }
    else if (keyword == "DistanceCutOff")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "Entrez plus de parametres dans DistanceCutOff" << endl;
	    cout << "Par exemple DistanceCutOff = d" << endl;
	    cout << "Parametres actuels : " << endl << param << endl;
	    abort();
	  }

        dist_tronc = to_num<Real_wp>(param(0));
        dist_tronc_min = 0.25*dist_tronc;
	if (param.GetM() >= 2)
	  {
	    if (param(1) == "GLOBAL")
	      separated_interfaces = false;
	    else
	      separated_interfaces = true;
	  }
      }
    else if (keyword == "StoreJacobiMatrixIteration")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "Entrez plus de parametres dans StoreJacobiMatrixIteration" << endl;
	    cout << "Par exemple StoreJacobiMatrixIteration = YES" << endl;
	    cout << "Parametres actuels : " << endl << param << endl;
	    abort();
	  }

        if (param(0) == "YES")
          compute_jacobi_iteration_matrix = true;
        else
          compute_jacobi_iteration_matrix = false;
      }
    else if (keyword == "ComputeGinOnQuadraturePoints")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "Entrez plus de parametres dans ComputeGinOnQuadraturePoints" << endl;
	    cout << "Par exemple ComputeGinOnQuadraturePoints = YES" << endl;
	    cout << "Parametres actuels : " << endl << param << endl;
	    abort();
	  }

        if (param(0) == "YES")
          {
            compute_g_on_quadrature_points = true;
            cout << "no longer correct" << endl;
            abort();
          }
        else
          compute_g_on_quadrature_points = false;
      }
    else if (keyword == "SubdomainIterativeAlgorithm")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "Entrez plus de parametres dans SubdomainIterativeAlgorithm" << endl;
	    cout << "Par exemple SubdomainIterativeAlgorithm = Gmres" << endl;
	    cout << "Parametres actuels : " << endl << param << endl;
	    abort();
	  }
	
        if (param(0) == "GaussSeidel")
	  type_resolution = GAUSS_SEIDEL_SOLVER;
	else if (param(0) == "Jacobi")
	  type_resolution = JACOBI_SOLVER;
	else if (param(0) == "Gmres")
	  type_resolution = GMRES_SOLVER;
	else if (param(0) == "Precond")
	  type_resolution = PRECOND_SOLVER;
        else
          {
	    cout << "Unknown solver" << endl;
	    abort();
	  }
      }
    else if (keyword == "FiniteElementQuadrature")
      {
	if (param(0) == "Gauss")
	  QuadrangleGauss::type_quadrature_default
	    = Globatto<Real_wp>::QUADRATURE_GAUSS;
	else if (param(0) == "GaussBlended")
	  {
	    QuadrangleGauss::type_quadrature_default
	      = Globatto<Real_wp>::QUADRATURE_GAUSS_BLENDED;
	    
	    if (param.GetM() >= 2)
	      Globatto<Real_wp>::blending_default = to_num<Real_wp>(param(1));
	  }
	else if (param(1) == "Lobatto")
	  QuadrangleGauss::type_quadrature_default
	    = Globatto<Real_wp>::QUADRATURE_LOBATTO;
	else
	  {
	    cout << "Unknown quadrature" << endl;
	    abort();
	  }
      }
    else if (keyword == "ResidualNorm")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "Entrez plus de parametres dans ResidualNorm" << endl;
	    cout << "Par exemple ResidualNorm = Inverse" << endl;
	    cout << "Parametres actuels : " << endl << param << endl;
	    abort();
	  }

        if (param(0) == "Inverse")
          compare_solution = true;
        else
          compare_solution = false;
      }
    else if (keyword == "PeriodicityX")
      {
	if (param.GetM() <= 1)
	  {
	    cout << "Entrez plus de parametres dans PeriodicityX" << endl;
	    cout << "Par exemple PeriodicityX = YES L" << endl;
	    cout << "Parametres actuels : " << endl << param << endl;
	    abort();
	  }	
	
	if (param(0) == "YES")
	  x_periodic = true;
	else
	  x_periodic = false;
	
	L_periodic = to_num<Real_wp>(param(1));
      }
  }
  
  
  //! returns distance between x and y
  Real_wp GetDistance(const R2& x, const R2& y) const
  {
    Real_wp dist = x.Distance(y);
    if (x_periodic)
      {
	dist = min(dist, x.Distance(R2(y(0)+L_periodic, y(1))));
	dist = min(dist, x.Distance(R2(y(0)-L_periodic, y(1))));
      }
    
    return dist;
  }

  Real_wp GetDistance(const R3& x, const R3& y) const
  {
    return x.Distance(y);
  }
  

  //! returns difference v-u
  R2 GetDiff(const R2& u, const R2& v) const
  {
    R2 diff(v-u);
    if (x_periodic)
      {
	Real_wp L = Norm2(diff);
	
	R2 diff_n(diff);
	diff_n(0) -= L_periodic;
	Real_wp Ln = Norm2(diff_n);

	R2 diff_p(diff);
	diff_p(0) += L_periodic;
	Real_wp Lp = Norm2(diff_p);
	
	if (Ln < min(L, Lp))
	  diff = diff_n;
	else if (Lp < min(L, Ln))
	  diff = diff_p;
      }
    
    //DISP(u); DISP(v); DISP(diff);
    
    return diff;
  }

  R3 GetDiff(const R3& u, const R3& v) const
  {
    return R3(v-u);
  }
  
  
  //! returns cut-off function xsi(x, y)
  template<int d>
  Real_wp GetCutOffFunction(const TinyVector<Real_wp, d>& x,
			    const TinyVector<Real_wp, d>& y) const
  {
    Real_wp r = GetDistance(x, y);
    if (r > dist_tronc)
      return Real_wp(0);
    else if (r < dist_tronc_min)
      {
        return Real_wp(1);
      }
    else 
      {
        Real_wp A0 = 0.5;
        Real_wp A1 = -1.875;
        Real_wp A3 = 5.0; 
        Real_wp A5 = -6.0;
        Real_wp arg = (r - dist_tronc_min) / (dist_tronc - dist_tronc_min) - 0.5;
        Real_wp arg2 = arg*arg;
        return A0 + arg*(A1 + arg2*(A3 + arg2*A5));
      }
  }
  
};


//! fills RefDof, IndexDof, IsDofOnProc0 on a local boundary
template<class TypeEquation, class Dimension>
void FindLocalSurfaceDof(EllipticProblem<TypeEquation>& var,
                         SubdomainSolverParameter& input_var,
                         int i, int num_face, int ne, int num_loc,
                         const Mesh<Dimension>& mesh, int& nquad,
                         IVect& IndexDof, IVect& RefDof, Vector<bool>& IsDofOnProc0,
                         int& nb_ddl_surf, const ElementReference<Dimension, 1>& Fb)
{
  IVect Nodle = var.GetDofNumberOnElement(ne);
  nquad = Fb.GetNbQuadBoundary(num_loc);
  for (int j = 0; j < Fb.GetNbDofBoundary(num_loc); j++)
    {
      int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
      int num_dof = Nodle(num_dof_loc);      
      if (num_dof >= 0)
        {
          if (IndexDof(num_dof) == -1)
            {
              int proc = var.GetProcessorNeighboringFace(num_face);
              if (proc == input_var.ref_neighbor_domain(0))
                IsDofOnProc0(num_dof) = true;
              else if (proc == input_var.ref_neighbor_domain(1))
                IsDofOnProc0(num_dof) = false;
              else
                {
                  cout << "Unknown neighboring processor" << endl;
                  abort();
                }
              
              RefDof(num_dof) = var.mesh.Element(ne).GetReference();
              IndexDof(num_dof) = nb_ddl_surf++;
            }
        }
    }
}
  

//! fills arrays IndexDof, RefDof, IsLocalDofOnProc0
/*!
  ref_cond : Gamma consists of edges/faces such that ref_cond(ref) = 1
  where ref is the reference of each edge/face
  IndexDof(i) is the local dof number on Gamma of the global dof i
  RefDof(i) is the reference of the domain associated with local dof i on Gamma
  IsLocalDofOnProc0 : if IsLocalDofOnProc0(i) is true, the local dof i on Gamma
  belongs to interface 0, otherwise to interface 1
  nb_ddl_surf : number of dofs on Gamma
*/
template<class TypeEquation, class Dimension>
void RetrieveSurfaceNumber(EllipticProblem<TypeEquation>& var,
                           SubdomainSolverParameter& input_var,
                           const IVect& ref_cond, IVect& IndexDof, IVect& RefDof,
                           Vector<bool>& IsLocalDofOnProc0, MeshInterpolationFEM<Dimension>& var_interp,
			   Vector<bool>& IsLocalQuadOnProc0, int& nb_ddl_surf)
{
  Mesh<Dimension>& mesh = var.mesh;
  //MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    
  IVect RefDofGlobal(var.GetNbDof());
  Vector<bool> IsDofOnProc0(var.GetNbDof());
  IsDofOnProc0.Fill(false);
  IndexDof.Reallocate(var.GetNbDof());
  IndexDof.Fill(-1);
  nb_ddl_surf = 0;
  IsLocalQuadOnProc0.Reallocate(var_interp.GetNbPointsQuadrature());
  int offset_quad = 0;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = mesh.BoundaryRef(i).GetReference();
      if ((ref >= 0) && (ref_cond(ref) == 1))
	{
	  int num_face = i;
	  int ne = mesh.Boundary(num_face).numElement(0);
	  if (ne >= 0)
	    {
	      int num_loc = mesh.Element(ne).GetPositionBoundary(num_face);
	      int nb_quad_edge(0);
	      FindLocalSurfaceDof(var, input_var, i, num_face, ne, num_loc, mesh, nb_quad_edge,
				  IndexDof, RefDofGlobal, IsDofOnProc0, nb_ddl_surf, var.GetReferenceElementH1(ne));
	      
	      int proc = var.GetProcessorNeighboringFace(num_face);
	      if (proc == input_var.ref_neighbor_domain(0))
		for (int j = 0; j < nb_quad_edge; j++)
		  IsLocalQuadOnProc0(offset_quad + j) = true;
	      else if (proc == input_var.ref_neighbor_domain(1))
		for (int j = 0; j < nb_quad_edge; j++)
		  IsLocalQuadOnProc0(offset_quad + j) = false;
		
	      offset_quad += nb_quad_edge;
	    }
	}
    }
  
  // renumbering IndexDof
  // replacing global dofs by local dofs for RefDof, IsLocalDofOnProc0
  IsLocalDofOnProc0.Reallocate(nb_ddl_surf);
  RefDof.Reallocate(nb_ddl_surf);
  nb_ddl_surf = 0;
  for (int i = 0; i < IndexDof.GetM(); i++)
    if (IndexDof(i) >= 0)
      {
        IsLocalDofOnProc0(nb_ddl_surf) = IsDofOnProc0(i);
        RefDof(nb_ddl_surf) = RefDofGlobal(i);
        IndexDof(i) = nb_ddl_surf++;  
      }
}

//! main class for solving Helmholtz equation with DDM method with an integral operator (or local operator)
template<class TypeEquation>
class SubdomainMatrixFreeOperator : public InputDataProblem_Base, public VirtualMatrix<Complex_wp>,
				    public Preconditioner_Base<Complex_wp>
{
public :
  typedef typename TypeEquation::Dimension Dimension;
  typedef typename Dimension::DimensionBoundary DimensionB;
  typedef typename TypeEquation::Complexe Complexe;
  typedef Vector<Complexe> VectComplexe;
  
  //! parametres d'entree pour la methode avec operateur
  SubdomainSolverParameter input_var;
  // objet pour calculer u sur les points d'integration
  MeshInterpolationFEM<Dimension> var_interp;
  
  //! surface mesh of Gamma
  Mesh<Dimension> mesh_surf;
  
  // alias vers le probleme element fini
  EllipticProblem<TypeEquation>& var;
  // IndexDofGamma(i) : numero local dans Gamma du ddl global i
  // RefDofGamma(i) : reference du domaine ou est le ddl local i
  IVect IndexDofGamma, RefDofGamma;
  // nombre de ddls dans Gamma
  int nb_ddl_gamma;
  // liste des ddls globaux qui sont sur Gamma
  IVect ListeDof;
  // position de tous les ddls de Gamma
  Vector<typename Dimension::R_N> CoorDof;
  // IsDofOnProc0(i) = true si le ddl local i est sur l'interface 0
  Vector<bool> IsDofOnProc0, IsQuadPointOnProc0;
  Vector<Vector<typename Dimension::R_N> > NormaleNodal;
  Vector<Vector<TinyMatrix<Real_wp, General, 3, 2> > > MatDfjm1Nodal, MatDfjm1Quad;

  // nombre de ddls pour l'inconnue phi, et taille pour gin
  int nb_ddl_phi, nb_ddl_gin;
  // matrices L et M
  Matrix<Real_wp, General, ArrayRowSparse> matL, matM;
  Matrix<Complex_wp, General, ArrayRowSparse> matLstiff;

  // factorisation de la matrice de masse M
  mutable All_MatrixLU<Real_wp> mat_mass_lu;
  // factorisation de la matrice locale elements finis (comprenant la partie avec T)
  mutable All_MatrixLU<Complexe> mat_lu;

  // u et du/dn sur les points d'integration
  mutable Vector<Complexe> trace_g_in_quad;
  // vecteurs intermediaires
  mutable Vector<Complexe> g_in_quad, u_gamma, Mu_gamma;
  mutable Vector<Complexe> Lphi_sol, g_in_vol;
  
  // coefficients d'impedance pour l'interface 0 et l'interface 1
  Complex_wp coef_imped0, coef_imped1, coef_c0, coef_c1;

  // solution sur tous les ddls du domaine
  mutable Vector<Complexe> x_sol_glob;

  // operateur pour le cas du COCG preconditionne
  VectReal_wp OperatorP;

  // objet pour evaluer les integrales singulieres
  SingularDoubleQuadratureGalerkin_Base<Real_wp, DimensionB> singular_int;
  
  //! constructor with finite element problem and data lines of the input file
  SubdomainMatrixFreeOperator(EllipticProblem<TypeEquation>& var_helm,
			      Vector<string> lines_data_file) : var(var_helm)
  {
    input_var.SetOmega(var.GetOmega());
    ReadInputFile(lines_data_file, *this);
    
    nb_ddl_phi = 0;
    nb_ddl_gin = 0;
    nb_ddl_gamma = 0;
  }

  void SetInputData(const string& keyword, const Vector<string>& param)
  {
    MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    input_var.SetInputData(keyword, param);
    
    if (keyword == "SingularIntegrationMethod")
      {
        if (param.GetM() < 3)
          {
            cout << "Entrez plus de parametres dans SingularIntegrationMethod" << endl;
            cout << "Par exemple SingularIntegrationMethod = ANALYTICAL eta AUTO" << endl;
            cout << "Parametres actuels : " << endl << param << endl;
            abort();
          }
        
	Real_wp eta = to_num<Real_wp>(param(1));
	if (param(0) == "ANALYTICAL")
	  {
	    int r_proche, r_joint;
	    if (param(2) == "AUTO")
	      {
		int r = mesh_num.GetOrder();
		r_proche = r+2;
		r_joint = r+11;
	      }
	    else
	      {
		if (param.GetM() < 4)
		  {
		    cout << "Entrez plus de parametres dans SingularIntegrationMethod" << endl;
		    cout << "Par exemple SingularIntegrationMethod = ANALYTICAL eta rp rj" << endl;
		    cout << "Parametres actuels : " << endl << param << endl;
		    abort();
		  }
		
		r_proche = to_num<int>(param(2));
		r_joint = to_num<int>(param(3));
	      }
	    
	    singular_int.SetAnalyticalIntegration(eta, r_proche, r_joint);
	  }
	else
	  {
	    int r_proche, r_joint, r_joint_ext, n_int, n_ext;
	    if (param(2) == "AUTO")
	      {
		int r = mesh_num.GetOrder();
		// best choice  (r/2r for straight elements)
		n_int = 2*r;
		if (input_var.type_operator == input_var.OPERATOR_DIFF)
		  n_int = 4*r;
		
		n_ext = r+10;
		r_proche = r+2;
		r_joint = r+11;
		r_joint_ext = 2*r+16;
	      }
	    else
	      {
		if (param.GetM() < 7)
		  {
		    cout << "Entrez plus de parametres dans SingularIntegrationMethod" << endl;
		    cout << "Par exemple SingularIntegrationMethod = NUMERICAL eta rp rj rj_ext" << endl;
		    cout << "Parametres actuels : " << endl << param << endl;
		    abort();
		  }
		
		r_proche = to_num<int>(param(2));
		r_joint = to_num<int>(param(3));
		r_joint_ext = to_num<int>(param(4));
		n_int = to_num<int>(param(5));
		n_ext = to_num<int>(param(6));
	      }
	    
	    singular_int.SetNumericalIntegration(eta, r_proche, r_joint, r_joint_ext, n_int, n_ext);  
	  }
      }
  }

  void ComputeJacobianMatrices(Dimension2& dim)
  {
    MatDfjm1Nodal.Reallocate(var_interp.GetNbBoundary());
    MatDfjm1Quad.Reallocate(var_interp.GetNbBoundary());
    NormaleNodal.Reallocate(var_interp.GetNbBoundary());
  }

  void ComputeLocalJacobianMatrices(int iloc, const ElementReference<Dimension2, 1>& elt,
				    bool is_curved, const VectR3& pts_nodal,
				    Vector<TinyMatrix<Real_wp, General, 3, 2> >& mat,
				    Vector<TinyMatrix<Real_wp, General, 3, 2> >& matQ)
  {
    int nb_nodes = elt.GetNbPointsNodalElt();
    mat.Reallocate(nb_nodes);
    
    VectR3 dF_dx, dF_dy, dF_dx_quad, dF_dy_quad;
    elt.GetGeometricElement().DFjElemNodal(pts_nodal, is_curved, dF_dx, dF_dy);
    elt.GetGeometricElement().DFjElemQuadrature(pts_nodal, is_curved, dF_dx, dF_dy,
						dF_dx_quad, dF_dy_quad);
    
    const Vector<R3>& NormaleQuadI = var_interp.NormaleQuadrature(iloc);
    R3 vec_u, Col1, Col2;
    NormaleNodal(iloc).Reallocate(nb_nodes);
    for (int i = 0; i < nb_nodes; i++)
      {
	TimesProd(dF_dx(i), dF_dy(i), vec_u);
	Real_wp coef = 1.0/DotProd(vec_u, vec_u);
	
	TimesProd(dF_dy(i), vec_u, Col1);
	TimesProd(vec_u, dF_dx(i), Col2);
	mat(i)(0, 0) = coef*Col1(0); mat(i)(1, 0) = coef*Col1(1); mat(i)(2, 0) = coef*Col1(2);
	mat(i)(0, 1) = coef*Col2(0); mat(i)(1, 1) = coef*Col2(1); mat(i)(2, 1) = coef*Col2(2);
	
	Mlt(1.0/Norm2(vec_u), vec_u);
	NormaleNodal(iloc)(i) = vec_u;
        if (DotProd(NormaleNodal(iloc)(i), NormaleQuadI(0)) < 0)
          NormaleNodal(iloc)(i) *= -1.0;
      }
    
    int nb_quad = elt.GetNbPointsQuadratureInside();
    matQ.Reallocate(nb_quad);
    for (int i = 0; i < nb_quad; i++)
      {
	TimesProd(dF_dx_quad(i), dF_dy_quad(i), vec_u);
	Real_wp coef = 1.0/DotProd(vec_u, vec_u);
	
	TimesProd(dF_dy_quad(i), vec_u, Col1);
	TimesProd(vec_u, dF_dx_quad(i), Col2);
	matQ(i)(0, 0) = coef*Col1(0); matQ(i)(1, 0) = coef*Col1(1); matQ(i)(2, 0) = coef*Col1(2);
	matQ(i)(0, 1) = coef*Col2(0); matQ(i)(1, 1) = coef*Col2(1); matQ(i)(2, 1) = coef*Col2(2);
      }
  }
  
  
  void ComputeJacobianMatrices(Dimension3& dim)
  {
    MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int order = mesh_num.GetOrder();
    Vector<const ElementReference_Dim<Dimension3>* > elt_vol;
    var.GetReferenceElement(elt_vol);

    const ElementReference<Dimension2, 1>* Fb_tri, *Fb_quad;
    Fb_tri = NULL;
    Fb_quad = NULL;

    for (int n = 0; n < elt_vol.GetM(); n++)
      {
	if (elt_vol(n)->GetHybridType() != 0)
	  Fb_quad = &dynamic_cast<const ElementReference<Dimension, 1>* >(elt_vol(n))->GetQuadrangularSurfaceFiniteElement();
	
	if (elt_vol(n)->GetHybridType() != 3)
	  Fb_tri = &dynamic_cast<const ElementReference<Dimension, 1>* >(elt_vol(n))->GetTriangularSurfaceFiniteElement();
      }
    
    MatDfjm1Nodal.Reallocate(var_interp.GetNbBoundary());
    MatDfjm1Quad.Reallocate(var_interp.GetNbBoundary());
    NormaleNodal.Reallocate(var_interp.GetNbBoundary());
    for (int i = 0; i < var_interp.GetNbBoundary(); i++)
      {
	int nb_vert = mesh_surf.BoundaryRef(i).GetNbVertices();
	bool is_curved = mesh_surf.IsBoundaryCurved(i);
	if (nb_vert == 3)
	  ComputeLocalJacobianMatrices(i, *Fb_tri,
				       is_curved, var_interp.PointsNodal(i),
				       MatDfjm1Nodal(i), MatDfjm1Quad(i));
	else
	  ComputeLocalJacobianMatrices(i, *Fb_quad,
				       is_curved, var_interp.PointsNodal(i),
				       MatDfjm1Nodal(i), MatDfjm1Quad(i));
      }
  }
  
  //! constructs surface Gamma, and arrays associated (IndexDofGamma, RefDofGamma, IsDofOnProc0)
  void ComputeSurfaceDofs()
  {
    MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int order = mesh_num.GetOrder();
    int type_quad = var.GetReferenceElement(0).GetTypeIntegrationEdge();
    var_interp.SetGaussQuadrature(order, type_quad);
    
    var_interp.InitProjectionSurface(var.mesh);
    IVect ref_cond(var.mesh.GetNbReferences()+1);
    ref_cond.Fill(0);
    for (int i = 0; i <= var.mesh.GetNbReferences(); i++)
      if (var.mesh.GetBoundaryCondition(i) == BoundaryConditionEnum::LINE_NEIGHBOR)
	ref_cond(i) = 1;
    
    var_interp.ComputeSurfaceMesh(ref_cond, var.mesh, mesh_surf, var);
    
    // on ne regroupe pas les points de quadrature avec les sous-domaines
    var_interp.GatherQuadraturePoints(var.comm_group_mode, false);
    
    RetrieveSurfaceNumber(var, input_var, ref_cond, IndexDofGamma,
                          RefDofGamma, IsDofOnProc0,
			  var_interp, IsQuadPointOnProc0, nb_ddl_gamma);
    
    cout << "Number of dofs in Gamma for processor " << rank_proc << " : " << nb_ddl_gamma << endl;
    
    // on construit la liste des ddls surfaciques
    if (nb_ddl_gamma > 0)
      {
	ListeDof.Reallocate(nb_ddl_gamma);  
	for (int i = 0; i < IndexDofGamma.GetM(); i++)
	  if (IndexDofGamma(i) >= 0)
	    ListeDof(IndexDofGamma(i)) = i;
	
	// on calcule la position des ddls
	Vector<typename Dimension::R_N> GlobalDofs;
	var.ComputeDofCoordinates(GlobalDofs);

	CoorDof.Reallocate(nb_ddl_gamma);
	for (int i = 0; i < ListeDof.GetM(); i++)
	  CoorDof(i) = GlobalDofs(ListeDof(i));
      }
    
    // en 3-D, calcul de DF^*-1 pour reconstruire le gradient
    Dimension dim;
    ComputeJacobianMatrices(dim);    
  }

  
  //! adds local part \int_{e_{iloc}} \int_{e_{jloc}} L(u, v) dx dy to matrix matL
  void AddLocalContributionL(int iloc, int jloc,
			     int num_elem_i, int num_loc_i, int num_elem_j, int num_loc_j,
			     VirtualMatrix<Real_wp>& oper_L, VirtualMatrix<Real_wp>& oper_M,
			     VirtualMatrix<Complex_wp>& oper_stiff,
			     const ElementReference<Dimension, 1>& Fb1,
			     const ElementReference<Dimension, 1>& Fb2)
  {
    int nb_dof = Fb1.GetNbDofBoundary(num_loc_i);
    //int nb_pts_quad = Fb1.GetNbQuadBoundary(num_loc_i);
    
    // elementary matrix
    Matrix<Real_wp> mat_elem(nb_dof, nb_dof), mat_mass(nb_dof, nb_dof);
    Matrix<Real_wp> mat_elem_diagI(nb_dof, nb_dof), mat_elem_diagJ(nb_dof, nb_dof);
    Matrix<Complex_wp> mat_stiff(nb_dof, nb_dof);
    
    // vertices of edges or faces
    int num_face_i = var.mesh.Element(num_elem_i).numBoundary(num_loc_i);
    int num_face_j = var.mesh.Element(num_elem_j).numBoundary(num_loc_j);
    Vector<typename Dimension::R_N> PtsAi(var.mesh.Boundary(num_face_i).GetNbVertices());
    Vector<typename Dimension::R_N> PtsAj(var.mesh.Boundary(num_face_j).GetNbVertices());
    for (int i = 0; i < var.mesh.Boundary(num_face_i).GetNbVertices(); i++)
      {
	int nv = MeshNumbering<Dimension>::
	  GetLocalVertexBoundary(var.mesh.Element(num_elem_i).GetHybridType(), num_loc_i, i);
	
	PtsAi(i) = var.mesh.Vertex(var.mesh.Element(num_elem_i).numVertex(nv));
      }
    
    for (int i = 0; i < var.mesh.Boundary(num_face_j).GetNbVertices(); i++)
      {
	int nv = MeshNumbering<Dimension>::
	  GetLocalVertexBoundary(var.mesh.Element(num_elem_j).GetHybridType(), num_loc_j, i);
	
	PtsAj(i) = var.mesh.Vertex(var.mesh.Element(num_elem_j).numVertex(nv));
      }
    
    // are edges/faces curved ?
    bool curved_edges = false;
    if (var.mesh.IsBoundaryCurved(num_face_i) || var.mesh.IsBoundaryCurved(num_face_j))
      curved_edges = true;
    
    // we retrieve ds on quadrature points
    const VectReal_wp& Dsi = var_interp.DsQuadrature(iloc);
    const VectReal_wp& Dsj = var_interp.DsQuadrature(jloc);
    
    // and quadrature points
    const Vector<typename Dimension::R_N>& PointsQuadI = var_interp.PointsQuadrature(iloc);
    const Vector<typename Dimension::R_N>& PointsQuadJ = var_interp.PointsQuadrature(jloc);
    const Vector<typename Dimension::R_N>& NormaleQuadI = var_interp.NormaleQuadrature(iloc);
    const Vector<typename Dimension::R_N>& NormaleQuadJ = var_interp.NormaleQuadrature(jloc);
    const Vector<typename Dimension::R_N>& PointsNodalI = var_interp.PointsNodal(iloc);
    const Vector<typename Dimension::R_N>& PointsNodalJ = var_interp.PointsNodal(jloc);
    
    // then we evaluate integral \int_0^1 \int_0^1 |X(s)-X(t)|^\beta f(s, t) J(s) J(t) ds dt (in 2-D)
    if (!input_var.LocalOperator())
      singular_int.ComputeElemMatrix(PtsAi, PtsAj, curved_edges, Dsi, Dsj,
				     PointsNodalI, PointsNodalJ,
				     PointsQuadI, PointsQuadJ, NormaleNodal(iloc),
				     NormaleNodal(jloc), NormaleQuadI, NormaleQuadJ,
				     MatDfjm1Nodal(iloc), MatDfjm1Nodal(jloc),
				     MatDfjm1Quad(iloc), MatDfjm1Quad(jloc),
				     input_var, mat_elem, mat_elem_diagI, mat_elem_diagJ);
    else
      {
	mat_elem.Fill(0);
	mat_stiff.Fill(0);
	// for local operator, evaluating stiffness matrix \int 1/k_\eps^2 du/ds dphi/ds ds
	if ((iloc == jloc) && (input_var.type_operator == input_var.OPERATOR_LOCAL))
	  {
	    int ref = var.mesh.Element(num_elem_i).GetReference();
	    Complex_wp rho = var.ref_rho(ref).GetConstant();
	    Complex_wp mu = var.ref_mu(ref).GetConstant()(0, 0);
            Complex_wp k_inf = var.GetOmega()*sqrt(rho/mu);
	    
	    singular_int.ComputeStiffnessMatrix(PtsAi, curved_edges, Dsi, NormaleQuadI,
						MatDfjm1Quad(iloc), k_inf,
						input_var.epsilon_pade, mat_stiff);
	  }
      }
    
    // mass matrix part
    if (iloc == jloc)
      singular_int.ComputeMassMatrix(PtsAi, curved_edges, Dsi, mat_mass);
    
    // elementary matrix is added into the global matrix oper_L
    MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    Real_wp coef_normalization = sqrt(2.0/pi_wp);
    for (int i = 0; i < nb_dof; i++)
      for (int j = 0; j < nb_dof; j++)
	{
	  int num_i = Fb1.GetLocalNumber(num_loc_i, i);
	  int num_glob_i = mesh_num.Element(num_elem_i).GetNumberDof(num_i);
	  int num_dof_i = IndexDofGamma(num_glob_i);

          int num_j = Fb2.GetLocalNumber(num_loc_j, j);
          int num_glob_j = mesh_num.Element(num_elem_j).GetNumberDof(num_j);
          int num_dof_j = IndexDofGamma(num_glob_j);
	  
	  if (input_var.LocalOperator())
	    {
	      if ((iloc == jloc) && (input_var.type_operator == input_var.OPERATOR_LOCAL))
                oper_stiff.AddInteraction(num_dof_i, num_dof_j, mat_stiff(i, j));
	    }
	  else
	    {
	      mat_elem(i, j) *= coef_normalization;
	      oper_L.AddInteraction(num_dof_i, num_dof_j, mat_elem(i, j));
	    }
	  
	  if (iloc == jloc)
	    oper_M.AddInteraction(num_dof_i, num_dof_j, mat_mass(i, j));
	  else
	    {
	      if (input_var.type_operator == input_var.OPERATOR_DIFF)
		{
		  if (num_dof_i != num_dof_j)
		    oper_L.AddInteraction(num_dof_j, num_dof_i, mat_elem(i, j));
		  
		  int num_i2 = Fb1.GetLocalNumber(num_loc_i, j);
		  int num_glob_i2 = mesh_num.Element(num_elem_i).GetNumberDof(num_i2);
		  int num_dof_i2 = IndexDofGamma(num_glob_i2);
		  
		  int num_j2 = Fb2.GetLocalNumber(num_loc_j, i);
		  int num_glob_j2 = mesh_num.Element(num_elem_j).GetNumberDof(num_j2);
		  int num_dof_j2 = IndexDofGamma(num_glob_j2);
		  
		  mat_elem_diagI(i, j) *= coef_normalization;
		  oper_L.AddInteraction(num_dof_i, num_dof_i2, mat_elem_diagI(i, j));

		  mat_elem_diagJ(i, j) *= coef_normalization;
		  oper_L.AddInteraction(num_dof_j, num_dof_j2, mat_elem_diagJ(i, j));
		}
	    }
	}
  }
  

  void ConstructSingularIntegrationObject(Dimension2& dim, int r)
  {
    singular_int.SetFiniteElement(var.GetReferenceElementH1(0).GetSurfaceFiniteElement(0));
  }
  
  void ConstructSingularIntegrationObject(Dimension3& dim, int r)
  {
    Vector<const ElementReference<Dimension3, 1>* > elt_vol;
    var.GetReferenceElement(elt_vol);

    const ElementReference<Dimension2, 1>* Fb_tri, *Fb_quad;
    Fb_tri = NULL;
    Fb_quad = NULL;

    for (int n = 0; n < elt_vol.GetM(); n++)
      {
	if (elt_vol(n)->GetHybridType() != 0)
	  Fb_quad = &elt_vol(n)->GetQuadrangularSurfaceFiniteElement();
	
	if (elt_vol(n)->GetHybridType() != 3)
	  Fb_tri = &elt_vol(n)->GetTriangularSurfaceFiniteElement();
      }

    if (Fb_tri != NULL)
      singular_int.SetTriangularFiniteElement(*Fb_tri,
					      var_interp.PointsReferenceSurface()(0),
					      var_interp.WeightsReferenceSurface(0));
    
    if (Fb_quad != NULL)
      singular_int.SetFiniteElement2(*Fb_quad,
				     var_interp.PointsReferenceSurface()(1),
				     var_interp.WeightsReferenceSurface(1));
  }
  
  
  //! Computes needed operators (L and M)
  void ComputeNeededOperators()
  {
    if (!input_var.LocalOperator())      
      {
	nb_ddl_phi = nb_ddl_gamma;
        int d = 0;
        if (Dimension::dim_N == 3)
          d = 1;
        
	if (input_var.type_operator == input_var.OPERATOR_GRAD)
	  singular_int.SetOperator(singular_int.OPERATOR_GRAD, 0.5-d);
	else if (input_var.type_operator == input_var.OPERATOR_DIFF)
	  singular_int.SetOperator(singular_int.OPERATOR_DIFF, -1.5-d);
	else if (input_var.type_operator == input_var.OPERATOR_T)
	  {
	    nb_ddl_phi = 0;
	    singular_int.SetOperator(singular_int.OPERATOR_DIFF, -2.0-d);
	  }        
	else
	  singular_int.SetOperator(singular_int.OPERATOR_DS, 0.5-d);        
      }
    else
      {
	nb_ddl_phi = 0;
	if (input_var.type_operator == input_var.OPERATOR_LOCAL)
	  nb_ddl_phi = input_var.CoefAl.GetM()*nb_ddl_gamma;
      }

    MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int r = mesh_num.GetOrder();
    Dimension dim;
    ConstructSingularIntegrationObject(dim, r);
    
    //DISP(input_var.coef_z); DISP(input_var.coef_c);

    var.var_chrono.SetMessage("IntegralOperator", "compute the integral operator");
    var.var_chrono.Start("IntegralOperator");

    matL.Reallocate(nb_ddl_gamma, nb_ddl_gamma);
    matM.Reallocate(nb_ddl_gamma, nb_ddl_gamma);  
    if (input_var.type_operator == input_var.OPERATOR_LOCAL)
      matLstiff.Reallocate(nb_ddl_gamma, nb_ddl_gamma);
    
    for (int i = 0; i < var_interp.GetNbBoundary(); i++)
      for (int j = 0; j < var_interp.GetNbBoundary(); j++)
        {
	  Real_wp distance(1e100);
          for (int iloc = 0; iloc < mesh_surf.BoundaryRef(i).GetNbVertices(); iloc++)
	    for (int jloc = 0; jloc < mesh_surf.BoundaryRef(j).GetNbVertices(); jloc++)
	      {
		int ni = mesh_surf.BoundaryRef(i).numVertex(iloc);
		int nj = mesh_surf.BoundaryRef(j).numVertex(jloc);
		distance = min(distance, input_var.GetDistance(mesh_surf.Vertex(ni),
							       mesh_surf.Vertex(nj)));
	      }
	  
          int num_elem1 = var_interp.GetElementNumberOfSurface(i);
          int num_loc1 = var_interp.GetLocalPositionOfSurface(i);
          int num_face1 = var.mesh.Element(num_elem1).numBoundary(num_loc1);
          int proc1 = var.GetProcessorNeighboringFace(num_face1);
	  
          int num_elem2 = var_interp.GetElementNumberOfSurface(j);
          int num_loc2 = var_interp.GetLocalPositionOfSurface(j);
          int num_face2 = var.mesh.Element(num_elem2).numBoundary(num_loc2);
          int proc2 = var.GetProcessorNeighboringFace(num_face2);
	  
          if (input_var.LocalOperator())
            {
              if (i != j)
                distance = input_var.dist_tronc*2;
            }
          
          if ((distance < input_var.dist_tronc)
	      && ( (!input_var.separated_interfaces) || (proc1 == proc2)))
            {
              int num_elem_i = var_interp.GetElementNumberOfSurface(i);
              int num_loc_i = var_interp.GetLocalPositionOfSurface(i);
              int num_elem_j = var_interp.GetElementNumberOfSurface(j);
              int num_loc_j = var_interp.GetLocalPositionOfSurface(j);
              
	      AddLocalContributionL(i, j, num_elem_i, num_loc_i, num_elem_j, num_loc_j,
				    matL, matM, matLstiff,
				    var.GetReferenceElementH1(num_elem_i),
				    var.GetReferenceElementH1(num_elem_j));              
            }
        }
    
    // on applique la fonction de troncature a l'operateur L
    if (!input_var.LocalOperator())      
      {
        for (int i = 0; i < matL.GetM(); i++)
          for (int j2 = 0; j2 < matL.GetRowSize(i); j2++)
            {
              int j = matL.Index(i, j2);
              Real_wp coef_tronc = input_var.GetCutOffFunction(CoorDof(i), CoorDof(j));
              matL.Value(i, j2) *= coef_tronc;
            }

        if (var.print_level >= 8)
	  {
            WriteNoBrackets(CoorDof, "CoorDof"+to_str(rank_proc)+".dat");
	    matM.WriteText("M"+to_str(rank_proc)+".dat");
	    matLstiff.WriteText("Lstiff"+to_str(rank_proc)+".dat");
	    matL.WriteText("L"+to_str(rank_proc)+".dat");      
	  }
      }
    else
      matL.Clear();
    
    var.var_chrono.Stop("IntegralOperator");
  }
  
  
  // calcul et factorisation des matrices pour chaque sous-domaine
  template<class MatrixSparse>
  void ComputeLocalMatrix(FemMatrixFreeClass<Complex_wp, TypeEquation>& free_mat,
                          GlobalGenericMatrix<Complex_wp>& nat_mat,
                          MatrixSparse& mat_direct, int sizeA)
  {
    mat_direct.Reallocate(sizeA, sizeA);
    
    // surface integrals
    ImpedanceABC<Complex_wp, TypeEquation> fimped(var);
    
    // absorbing condition on the extern boundary only
    var.AddMatrixImpedanceBoundary(Complex_wp(1, 0), var.mesh.GetBoundaryCondition(),
				   BoundaryConditionEnum::LINE_ABSORBING,
				   nat_mat, mat_direct, 0, 0, fimped);
    
    // impedance condition
    ImpedanceGeneric<Complex_wp, TypeEquation> fimped_gen(var);
    var.AddMatrixImpedanceBoundary(Complex_wp(1, 0), var.mesh.GetBoundaryCondition(),
				   BoundaryConditionEnum::LINE_IMPEDANCE,
				   nat_mat, mat_direct, 0, 0, fimped_gen);
    
    // volume integrals
    var.AddMatrixFEM(free_mat, mat_direct, nat_mat, 0, 0);
    
    mat_direct.Init(mat_direct.GetM(), NULL, NULL, NULL,
                    mat_direct.GetM(), 1, NULL, NULL, MPI_COMM_SELF);
    
    
    int Nvol = var.GetNbDof();  
    int nb_ddl_gamma = matM.GetM();
    // on construit la liste des ddls surfaciques
    IVect ListeDof(nb_ddl_gamma);  
    for (int i = 0; i < IndexDofGamma.GetM(); i++)
      if (IndexDofGamma(i) >= 0)
        ListeDof(IndexDofGamma(i)) = i;
    
    // partie impedance commune a Despres et l'operateur
    // coef0 = - i omega sqrt(rho mu) z 
    // et similairement pour coef1
    Complex_wp k0(1, 0), k1(1, 0), mu0(1, 0), mu1(1, 0);
    coef_imped0 = 0, coef_imped1 = 0;
    int ref_domain0 = -1, ref_domain1 = -1;
    for (int i = 0; i < ListeDof.GetM(); i++)
      {
        // on a mu du/dn = mu (i k z u)
        // avec k = omega sqrt(rho/mu)
        // donc mu du/dn = i omega z sqrt(rho mu)
        // var.coefficient_impedance_absorbing = sqrt(rho mu)
        Complex_wp coef(0);
        if (IsDofOnProc0(i))
          {
            if (ref_domain0 == -1)
              {
                ref_domain0 = RefDofGamma(i);
		mu0 = var.ref_mu(ref_domain0).GetConstant()(0, 0);
                k0 = var.GetOmega()*var.coefficient_impedance_absorbing(ref_domain0);
                coef_imped0 = -Iwp*k0*input_var.coef_z(0);
              }
            else
              {
                if (input_var.separated_interfaces)
		  if (RefDofGamma(i) != ref_domain0)
		    {
		      cout << "Hetereogeneous interface not treated " << endl;
		      abort();
		    }
              }
            
            coef = coef_imped0;
          }
        else
          {
            if (ref_domain1 == -1)
              {
                ref_domain1 = RefDofGamma(i);
		mu1 = var.ref_mu(ref_domain1).GetConstant()(0, 0);
                k1 = var.GetOmega()*var.coefficient_impedance_absorbing(ref_domain1);
                coef_imped1 = -Iwp*k1*input_var.coef_z(1);
              }
            else
              {
                if (input_var.separated_interfaces)
		  if (RefDofGamma(i) != ref_domain1)
		    {
		      cout << "Hetereogeneous interface not treated " << endl;
		      abort();
		    }
              }
            
            coef = coef_imped1;
          }
        
	if (input_var.type_operator == input_var.OPERATOR_LOCAL)
	  {
	    // part -ik (C0 + \sum A_l/B_l) \int u v
            Complex_wp C0_tilde = input_var.CoefC0;
            for (int l = 0; l < input_var.CoefAl.GetM(); l++)
              C0_tilde += input_var.CoefAl(l) / input_var.CoefBl(l);
            
	    for (int j = 0; j < matM.GetRowSize(i); j++)
	      mat_direct.AddInteraction(ListeDof(i), ListeDof(matM.Index(i, j)),
					coef*C0_tilde*matM.Value(i, j));
	    
	    // part ik D0 \int 1/k_eps^2 du/ds dv/ds
	    if (input_var.CoefD0 != Complex_wp(0, 0))
	      for (int j = 0; j < matLstiff.GetRowSize(i); j++)
		mat_direct.AddInteraction(ListeDof(i), ListeDof(matLstiff.Index(i, j)),
                                          -coef*input_var.CoefD0*matLstiff.Value(i, j));
	    
	    for (int l = 0; l < input_var.CoefAl.GetM(); l++)
	      {
		// part ik A_l/B_l \int psi^l v
                Complex_wp coef_tilde = input_var.CoefAl(l) / input_var.CoefBl(l);
		for (int j = 0; j < matM.GetRowSize(i); j++)
		  mat_direct.AddInteraction(ListeDof(i), Nvol + l*nb_ddl_gamma + matM.Index(i, j),
					    -coef*matM.Value(i, j)*coef_tilde);
		
		// part ik A_l /B_l \int u v_l - ik A_l/B_l \int \psi_l v_l                
		for (int j = 0; j < matM.GetRowSize(i); j++)
		  {
		    mat_direct.AddInteraction(Nvol + l*nb_ddl_gamma + i,
					      Nvol + l*nb_ddl_gamma + matM.Index(i, j),
					      coef*coef_tilde*matM.Value(i, j));

		    mat_direct.AddInteraction(Nvol + l*nb_ddl_gamma + i,
					      ListeDof(matM.Index(i, j)),
					      -coef*coef_tilde*matM.Value(i, j));
		  }
		
		//  part + ik A_l \int 1/k_eps^2 dpsi^l/ds dv_l/ds
		for (int j = 0; j < matLstiff.GetRowSize(i); j++)
		  mat_direct.AddInteraction(Nvol + l*nb_ddl_gamma + i,
					    Nvol + l*nb_ddl_gamma + matLstiff.Index(i, j),
					    -coef*matLstiff.Value(i, j)*input_var.CoefAl(l));
	      }
	  }
	else
	  {
	    int nb_interac = matM.GetRowSize(i);
	    for (int j = 0; j < nb_interac; j++)
	      mat_direct.AddInteraction(ListeDof(i), ListeDof(matM.Index(i, j)), matM.Value(i, j)*coef);
	  }
      }
    
    coef_c0 = input_var.coef_c(0)*mu0/k0;
    coef_c1 = input_var.coef_c(1)*mu1/k1;
    
    mat_direct.WriteText("mat_before" + to_str(rank_proc) + ".dat");
    
    // coefficient a appliquer a l'operateur L
    Complexe coef_L0 = coef_c0*coef_imped0;
    Complexe coef_L1 = coef_c1*coef_imped1;
    //DISP(coef_imped0); DISP(coef_imped1);
    //DISP(input_var.coef_c); DISP(mu0); DISP(mu1); DISP(k0); DISP(k1);
    //DISP(coef_L0); DISP(coef_L1);
    
    if (!input_var.LocalOperator())
      {
        // on rajoute coef_L L et coef_L L^T a la matrice
        IVect index(nb_ddl_gamma); Vector<Complexe> val(nb_ddl_gamma);

	if (input_var.type_operator == input_var.OPERATOR_T)
	  {
	    for (int i = 0; i < ListeDof.GetM(); i++)
	      {
		int nb_interac = matL.GetRowSize(i);
		for (int j = 0; j < nb_interac; j++)
		  {
		    index(j) = matL.Index(i, j);
		    if (IsDofOnProc0(i))
		      {
			if (input_var.separated_interfaces)
			  if (!IsDofOnProc0(matL.Index(i, j)))
			    {
			      cout << "This case should be impossible" << endl;
			      abort();
			    }
			
			val(j) = matL.Value(i, j)*coef_L0;
		      }
		    else
		      {
			if (input_var.separated_interfaces)
			  if (IsDofOnProc0(matL.Index(i, j)))
			    {
			      cout << "This case should be impossible" << endl;
			      abort();
			    }
			
			val(j) = matL.Value(i, j)*coef_L1;
		      }
		  }
		
		mat_direct.AddInteractionRow(ListeDof(i), nb_interac, index, val);
	      }
	  }
	else
	  {	    
	    for (int i = 0; i < ListeDof.GetM(); i++)
	      {
		int nb_interac = matL.GetRowSize(i);
		for (int j = 0; j < nb_interac; j++)
		  {
		    index(j) = Nvol + matL.Index(i, j);
		    if (IsDofOnProc0(i))
		      {
			if (input_var.separated_interfaces)
			  if (!IsDofOnProc0(matL.Index(i, j)))
			    {
			      cout << "This case should be impossible" << endl;
			      abort();
			    }
			
			val(j) = matL.Value(i, j)*coef_L0;
		      }
		    else
		      {
			if (input_var.separated_interfaces)
			  if (IsDofOnProc0(matL.Index(i, j)))
			    {
			      cout << "This case should be impossible" << endl;
			      abort();
			    }
			
			val(j) = matL.Value(i, j)*coef_L1;
		      }
		  }
		
		mat_direct.AddInteractionRow(ListeDof(i), nb_interac, index, val);
		
		// partie transposee
		mat_direct.AddInteractionColumn(ListeDof(i), nb_interac, index, val);
	      }
	    
	    // on rajoute M a la matrice
	    for (int i = 0; i < nb_ddl_gamma; i++)
	      {
		int nb_interac = matM.GetRowSize(i);
		for (int j = 0; j < nb_interac; j++)
		  {
		    index(j) = Nvol + matM.Index(i, j);
		    if (IsDofOnProc0(i))
		      val(j) = -coef_L0*matM.Value(i, j);
		    else
		      val(j) = -coef_L1*matM.Value(i, j);
		  }
		
		mat_direct.AddInteractionRow(Nvol + i, nb_interac, index, val);
	      }
	  }
      }
    
    mat_direct.WriteText("mat_after" + to_str(rank_proc) + ".dat");
    mat_lu.Factorize(mat_direct);
    
    int ierr;
    int type = mat_lu.GetInfoFactorization(ierr);
    PrintFactorizationFailed(type, ierr);
  }


  //! computation and factorization of finite element matrix (including operator T)
  void FactorizeOperators(FemMatrixFreeClass<Complex_wp, TypeEquation>& free_mat,
			  GlobalGenericMatrix<Complex_wp>& nat_mat)  
  {
    // resolution des problemes locaux
    var.var_chrono.SetMessage("LocalProblem", "factorize local problems");
    var.var_chrono.Start("LocalProblem");
    
    VectString param_direct_solver;
    
    mat_lu.SelectDirectSolver(param_direct_solver);
    int Nvol = var.GetNbDof();
    int sizeA = Nvol + nb_ddl_phi;
    if (var.IsSymmetricProblem())
      {
	cout << "Symmetric problem " << endl;
	DistributedMatrix<Complexe, Symmetric, ArrayRowSymSparse> mat_direct;      
	ComputeLocalMatrix(free_mat, nat_mat, mat_direct, sizeA);
      }
    else
      {      
	DistributedMatrix<Complexe, General, ArrayRowSparse> mat_direct;
	ComputeLocalMatrix(free_mat, nat_mat, mat_direct, sizeA);
      }
    
    var.var_chrono.Stop("LocalProblem");
    
    // factorisation de la matrice de masse
    if ((input_var.type_operator != input_var.DESPRES) || input_var.compute_jacobi_iteration_matrix)
      {
	if ((input_var.compute_g_on_quadrature_points) || input_var.compute_jacobi_iteration_matrix)
	  {
	    mat_mass_lu.Factorize(matM, true);
	    
	    int ierr;
	    int type = mat_lu.GetInfoFactorization(ierr);
	    PrintFactorizationFailed(type, ierr);
	  }
      }
  }
  
  //! returns the size of finite element matrix factorized
  int GetNbDof() const
  {
    int Nvol = var.GetNbDof();
    int sizeA = Nvol + nb_ddl_phi;
    return sizeA;
  }  
  
  //! returns the number of degrees of freedom for phi
  int GetNbDdlPhi() const
  {
    return nb_ddl_phi;
  }
  
  //! returns the size of g^in / g^out
  int GetNbDofGin() const
  {
    return nb_ddl_gin;
  }
  
  //! returns the size of linear system to solve for Gmres
  int GetM() const
  {
    return nb_ddl_gin;
  }
  

  //! allocation of intermediary vectors
  void InitVectors()
  {  
    // calcul de la solution (algo de Jacobi) sur les g
    int N = var_interp.GetNbPointsQuadrature(), Nquad = N;
    if (!input_var.compute_g_on_quadrature_points)
      {
	N = nb_ddl_gamma;
	trace_g_in_quad.Reallocate(Nquad);
	g_in_quad.SetData(Nquad, 
			  reinterpret_cast<Complex_wp*>(trace_g_in_quad.GetData()));
	
	u_gamma.Reallocate(nb_ddl_gamma); Mu_gamma.Reallocate(nb_ddl_gamma);
	u_gamma.Fill(0); Mu_gamma.Fill(0);
      }

    g_in_quad.Fill(0);

    // vecteurs intermediaires qu'on alloue une seule fois
    MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    Lphi_sol.Reallocate(Nvol); g_in_vol.Reallocate(Nvol);
    Lphi_sol.Fill(0); g_in_vol.Fill(0);      
    
    nb_ddl_gin = N;
    this->m_ = nb_ddl_gin;
    this->n_ = nb_ddl_gin;

    if (input_var.type_resolution == input_var.GMRES_SOLVER)
      {
	x_sol_glob.Reallocate(Nvol+nb_ddl_phi);
	x_sol_glob.Fill(0);
      }

    if (input_var.type_resolution == input_var.PRECOND_SOLVER)
      {
	x_sol_glob.Reallocate(Nvol+nb_ddl_phi);
	x_sol_glob.Fill(0);

	// operator P
	OperatorP.Reallocate(Nvol);
	OperatorP.Fill(1.0);
	var.AddDomains(OperatorP);
	for (int i = 0; i < var.GetNbDof(); i++)
	  OperatorP(i) = 1.0/OperatorP(i);
      }
  }

  //! computes source due to g_in
  /*!
    The boundary condition reads
    mu du/dn - i omega sqrt(rho mu) z u - i z c mu T(u) = g_in
    then the source is equal to x_sol(i) = \int_Gamma g_in \varphi_i ds
   */
  void ComputeSourceGin(const Vector<Complexe>& g_in, Vector<Complexe>& x_sol) const
  {
    Vector<Complexe> trace_du_dn;
    x_sol.Fill(0);
    if (input_var.compute_g_on_quadrature_points)
      var_interp.AddSourceBoundary(var, g_in, x_sol);
    else
      {
	// etape intermediaire pour interpoler g des ddls
	// aux points d'integration : on appelle ComputeEnHnOnBoundary
	if (input_var.type_operator == input_var.DESPRES)
	  {
	    for (int i = 0; i < nb_ddl_gamma; i++)
	      g_in_vol(ListeDof(i)) = g_in(i);
	    
	    var_interp.ComputeEnHnOnBoundary(var, g_in_vol, trace_g_in_quad, trace_du_dn, false, false);
            
	    // on rajoute la partie surfacique
	    var_interp.AddSourceBoundary(var, g_in_quad, x_sol);
	  }
	else
	  {
	    // cas ou on a M g_in au lieu de g_in
	    // => on le rajoute directement a la source
	    for (int i = 0; i < nb_ddl_gamma; i++)
	      x_sol(ListeDof(i)) += g_in(i);
	  }
      }
  }
  
  //! on input the source b, on output the solution x of linear system A x = b
  //! where A is the finite element matrix with additional terms due to the operator
  void ComputeLocalSolution(Vector<Complexe>& x_sol) const
  {
    mat_lu.Solve(x_sol);
  }
  
  //! we compute Lphi = L phi on degrees of freedom of g_in (dofs or quadrature points)
  void ComputeLphi(const Vector<Complexe>& x_sol, Vector<Complexe>& phi,
		   Vector<Complexe>& trace_Lphi, Vector<Complexe>& Lphi) const
  {
    MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    if (!input_var.LocalOperator())
      {
	if (input_var.type_operator == input_var.OPERATOR_T)
	  {
	    for (int i = 0; i < nb_ddl_gamma; i++)
	      phi(i) = x_sol(ListeDof(i));
	    
	    Mlt(matL, phi, Lphi);
	  }
	else
	  {	
	    // on extrait phi de x_sol
	    for (int i = 0; i < nb_ddl_phi; i++)
	      phi(i) = x_sol(Nvol + i);
	    
	    // on calcule L phi
	    Mlt(matL, phi, Lphi);
	  }    
	
	if (input_var.compute_g_on_quadrature_points)
	  {
	    // calcul de phi sur les points de quadrature
	    Vector<Complexe> trace_du_dn;
	    
	    // on inverse par la masse
	    mat_mass_lu.Solve(Lphi);
	    
	    for (int i = 0; i < nb_ddl_phi; i++)
	      Lphi_sol(ListeDof(i)) = Lphi(i);
	    
	    // on calcule Lphi sur les points d'integration
	    var_interp.ComputeEnHnOnBoundary(var, Lphi_sol, trace_Lphi, trace_du_dn, false, false);
	  }
      }
    else if (input_var.type_operator == input_var.OPERATOR_LOCAL)
      {
	// on extrait phi et u de x_sol
	Vector<Complexe> u(nb_ddl_gamma);
	for (int i = 0; i < nb_ddl_gamma; i++)
	  u(i) = x_sol(ListeDof(i));

	for (int i = 0; i < nb_ddl_phi; i++)
	  phi(i) = x_sol(Nvol + i);
	
	// on calcule L phi
	Mlt(matM, u, Lphi); Mlt(input_var.CoefC0, Lphi);
	if (input_var.CoefD0 != Complexe(0))
	  MltAdd(-input_var.CoefD0, matLstiff, u, Complexe(1), Lphi);
	
	Vector<Complexe> phi_l;
	for (int l = 0; l < input_var.CoefAl.GetM(); l++)
	  {
	    phi_l.SetData(nb_ddl_gamma, &phi(l*nb_ddl_gamma));
	    MltAdd(-input_var.CoefAl(l), matLstiff, phi_l, Complexe(1), Lphi); 
	    phi_l.Nullify();
	  }

	if (input_var.compute_g_on_quadrature_points)
	  {
	    // calcul de phi sur les points de quadrature
	    Vector<Complexe> trace_du_dn;
	    
	    // on inverse par la masse
	    mat_mass_lu.Solve(Lphi);
	    
	    for (int i = 0; i < nb_ddl_gamma; i++)
	      Lphi_sol(ListeDof(i)) = Lphi(i);
	    
	    // on calcule Lphi sur les points d'integration
	    var_interp.ComputeEnHnOnBoundary(var, Lphi_sol, trace_Lphi, trace_du_dn, false, false);
	  }	
      }
  }

  
  //! computation of g_out from g_in
  void UpdateValueGin(Complexe& g_in, int i, const Complexe& u, const Complexe& Lphi, bool domain0) const
  {
    if (domain0)
      g_in = -g_in + 2.0*Iwp*imag(coef_imped0)*(u + coef_c0*Lphi);
    else
      g_in = -g_in + 2.0*Iwp*imag(coef_imped1)*(u + coef_c1*Lphi);
  }


  //! computation of g_out from g_in
  void UpdateValueGinLoc(Complexe& g_in, int i, const Complexe& Lphi, bool domain0) const
  {
    if (domain0)
      g_in = -g_in + 2.0*coef_imped0*Lphi;
    else
      g_in = -g_in + 2.0*coef_imped1*Lphi;
  }
  
  
  //! computation of g_out from g_in, g_in is replaced by g_out
  void UpdateGin(Vector<Complexe>& g_in, const Vector<Complexe>& x_sol,
		 const Vector<Complexe>& trace_Lphi, const Vector<Complexe>& Lphi) const
  {
    Vector<Complexe> trace_u, trace_du_dn;
    int N = g_in.GetM(); Complexe zero;
    SetComplexZero(zero);
    if (input_var.compute_g_on_quadrature_points)
      {
	var_interp.ComputeEnHnOnBoundary(var, x_sol, trace_u, trace_du_dn, false, false);
	
	if (input_var.type_operator == input_var.DESPRES)
	  for (int i = 0; i < N; i++)
            UpdateValueGin(g_in(i), i, trace_u(i), zero, IsQuadPointOnProc0(i));
	else if (input_var.type_operator == input_var.OPERATOR_LOCAL)
	  for (int i = 0; i < N; i++)
            UpdateValueGinLoc(g_in(i), i, trace_Lphi(i), IsQuadPointOnProc0(i));
	else	  
	  for (int i = 0; i < N; i++)
            UpdateValueGin(g_in(i), i, trace_u(i), trace_Lphi(i), IsQuadPointOnProc0(i));
      }
    else
      {
	for (int i = 0; i < nb_ddl_gamma; i++)
	  u_gamma(i) = x_sol(ListeDof(i));
              
	if (input_var.type_operator == input_var.DESPRES)
	  for (int i = 0; i < N; i++)
            UpdateValueGin(g_in(i), i, u_gamma(i), zero, IsDofOnProc0(i));
	else if (input_var.type_operator == input_var.OPERATOR_LOCAL)
	  for (int i = 0; i < N; i++)
            UpdateValueGinLoc(g_in(i), i, Lphi(i), IsDofOnProc0(i));
	else
	  {
	    // produit par la masse              
	    Mlt(matM, u_gamma, Mu_gamma);
            
	    for (int i = 0; i < N; i++)
              UpdateValueGin(g_in(i), i, Mu_gamma(i), Lphi(i), IsDofOnProc0(i));
	  }
      }
  }
  
  
  //! the values are echanged and relaxed in the same time
  void ExchangeRelaxGinGout(Vector<Complexe>& g_in,
			    const Real_wp& omega, int proc) const
  {
    // on echange les valeurs entre sous-domaines      
    if (input_var.compute_g_on_quadrature_points)
      {
	// on echange les valeurs sur les points de quadrature
	var.ExchangeQuadRelaxDomains(g_in, omega, proc);
      }
    else
      {
	// on echange les valeurs sur les ddls
	for (int i = 0; i < nb_ddl_gamma; i++)
	  g_in_vol(ListeDof(i)) = g_in(i);
              
	var.ExchangeRelaxDomains(g_in_vol, omega, proc);
	
	for (int i = 0; i < nb_ddl_gamma; i++)
	  g_in(i) = g_in_vol(ListeDof(i));
      }
  }

  
  //! the values are echanged between interfaces
  void ExchangeGinGout(Vector<Complexe>& g_in) const
  {
    // on echange les valeurs entre sous-domaines      
    if (input_var.compute_g_on_quadrature_points)
      {
	Vector<Vector<Complexe> > g_send, g_recv;
	Vector<Vector<int64_t> > g_send_tmp, g_recv_tmp;
	Vector<MPI_Request> request; int tag = 23;

	// on echange les valeurs sur les points de quadrature
	var.ExchangeUfaceDomains(g_in, g_send, g_send_tmp, g_recv, g_recv_tmp, request, tag);
	var.GetUfaceDomains(g_in, g_send, g_send_tmp, g_recv, g_recv_tmp, request, tag);
      }
    else
      {
	// on echange les valeurs sur les ddls
	for (int i = 0; i < nb_ddl_gamma; i++)
	  g_in_vol(ListeDof(i)) = g_in(i);
              
	var.ExchangeDomains(g_in_vol);
	
	for (int i = 0; i < nb_ddl_gamma; i++)
	  g_in(i) = g_in_vol(ListeDof(i));
      }
  }

  
  //! computes iteration matrix for Jacobi algorithm
  void ComputeJacobiIterationMatrix(Vector<Complexe>& g_in, Vector<Complexe>& x_sol)
  {
    if (!input_var.compute_jacobi_iteration_matrix)
      return;
    
    // si on veut calculer la matrice d'iteration de Jacobi
    // on met source_rhs = 0, et on boucle sur les colonnes de la matrice
    ofstream file_out_jacobi;
    
    string name = "mat_jacobi_P"+to_str(rank_proc)+".dat";
    file_out_jacobi.open(name.data());
    file_out_jacobi.precision(15);
    
    Vector<Complexe> phi(nb_ddl_phi);
    Vector<Complexe> Lphi(nb_ddl_phi);
    if (input_var.type_operator == input_var.OPERATOR_T)
      {
	phi.Reallocate(nb_ddl_gamma);
	Lphi.Reallocate(nb_ddl_gamma);
      }
    else if (input_var.type_operator == input_var.OPERATOR_LOCAL)
      Lphi.Reallocate(nb_ddl_gamma);
    
    // traces de u et phi sur l'interface sigma (valeurs sur les points de quadrature)
    Vector<Complexe> trace_Lphi;

    for (int j = 0; j < g_in.GetM(); j++)
      {
	// pour calculer la matrice d'iteration
	// on fait le produit avec e_i
	g_in.Fill(0);
	g_in(j) = 1.0;
	
	// on part de g_in      
	// on actualise la source qu'on met dans x_sol
	// quelque part on resout -omega^2 u - Laplacien u = 0
	// et du/dn - i omega z T(u) = g_in 
	  
	// partie surfacique (due a g_in)
	ComputeSourceGin(g_in, x_sol);
	
	// resolution des problemes locaux, et calcul de Lphi
	ComputeLocalSolution(x_sol);
	ComputeLphi(x_sol, phi, trace_Lphi, Lphi);
          
	// on calcule g_out = -g_in - 2 ik u
	// g_in est ecrase par g_out
	UpdateGin(g_in, x_sol, trace_Lphi, Lphi);

	// on stocke la colonne de la matrice d'iteration
	for (int i = 0; i < g_in.GetM(); i++)
	  if (g_in(i) != Complex_wp(0, 0))
	    file_out_jacobi << i+1 << " " << j+1 << " " << g_in(i) << '\n';
      }

    file_out_jacobi.close();
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (nb_proc > 2)
      {
        // in that case, we write dof numbers for the two interfaces
        IVect Dof0(nb_ddl_gamma), Dof1(nb_ddl_gamma);
        int nb0 = 0, nb1 = 0;
        for (int i = 0; i < nb_ddl_gamma; i++)
          {
            if (IsDofOnProc0(i))
              Dof0(nb0++) = i;
            else
              Dof1(nb1++) = i;
          }
        
        Dof0.Resize(nb0);
        Dof1.Resize(nb1);
        
        Dof0.WriteText("num_gamma0_P"+to_str(rank_proc)+".dat");
        Dof1.WriteText("num_gamma1_P"+to_str(rank_proc)+".dat");
        
        return;
      }
    
    if (rank_proc == 0)
      {
	// computation of eigenvalues of L_0 L_1
	Matrix<Complex_wp> L0, L1;
	Matrix<Complex_wp, General, ArrayRowSparse> L0_sp, L1_sp;
	
	L0_sp.ReadText("mat_jacobi_P0.dat");
	L1_sp.ReadText("mat_jacobi_P1.dat");
	
	Copy(L0_sp, L0);
	Copy(L1_sp, L1);
	
	Matrix<Complex_wp> A(L0.GetM(), L0.GetM()), V;
	VectComplex_wp lambda;
	Mlt(L0, L1, A);
	
	GetEigenvaluesEigenvectors(A, lambda, V);
	
	lambda.Write("eigenval.dat");
	
	// estimating the maximal eigenvalue with relaxation
	VectComplex_wp mode(L0.GetM());
	Real_wp Tmax = 0.0; Complex_wp l;
	for (int i = 0; i < lambda.GetM(); i++)
	  {
	    for (int j = 0; j < mode.GetM(); j++)
	      mode(j) = V(j, i);
	    
	    mat_mass_lu.Solve(mode);
	    for (int j = 0; j < mode.GetM(); j++)
	      V(j, i) = mode(j);

	    l = sqrt(lambda(i));
	    Tmax = max(Tmax, abs(1.0-input_var.omega_relax + input_var.omega_relax*l));
	    Tmax = max(Tmax, abs(1.0-input_var.omega_relax - input_var.omega_relax*l));
	  }
	
	cout << "Taux attendu = " << Tmax << endl;
	
	// computing the angle theta for each dof
	VectReal_wp Theta(ListeDof.GetM());
	IVect Index(ListeDof.GetM());
	Index.Fill(); Real_wp r;
	for (int i = 0; i < ListeDof.GetM(); i++)
	  Theta(i) = CoorDof(i)(0);
	  //CartesianToPolar(CoorDof(i)(0), CoorDof(i)(1), r, Theta(i));
	
	Sort(Theta, Index);
	
	Theta.WriteText("theta.dat");
	
	// computing the absolute value of each value with relaxation
	// freq_mode will contain the mode number (based on FFT)
	VectReal_wp taux(lambda.GetM()), freq_mode(lambda.GetM());
	int N = mode.GetM();
	FftInterface<Complex_wp> fft;
	fft.Init(N);
	for (int i = 0; i < lambda.GetM(); i++)
	  {
	    l = sqrt(lambda(i));
	    Tmax = abs(1.0-input_var.omega_relax + input_var.omega_relax*l);
	    Tmax = max(Tmax, abs(1.0-input_var.omega_relax - input_var.omega_relax*l));
	    taux(i) = Tmax;
	    
	    for (int j = 0; j < N; j++)
	      mode(j) = V(Index(j), i);
	    
	    fft.ApplyForward(mode);
	    freq_mode(i) = 0;
	    Real_wp amax(0);
	    //freq_mode(i) += j*(abs(mode(j)) + abs(mode(N-j)));
	    for (int j = 1; j < N/2; j++)	    
	      if (amax < max(abs(mode(j)), abs(mode(N-j))))
		{
		  amax = max(abs(mode(j)), abs(mode(N-j)));
		  freq_mode(i) = j;
		}
	    
	    if (amax < abs(mode(0)))
	      {
		amax = abs(mode(0));
		freq_mode(i) = 0;
	      }	    
	  }
	
	// sorting the modes
	IVect IndexF(freq_mode.GetM());
	IndexF.Fill();
	Sort(freq_mode, IndexF);
	DISP(freq_mode);
	
	// the modes are written
	ofstream taux_out("taux.dat");
	taux_out.precision(15);
	for (int i = 0; i < lambda.GetM(); i++)
	  {
	    int n = IndexF(i);
	    taux_out << i << " " << taux(n) << '\n';

	    if (var.print_level >= 8)
              {
                for (int j = 0; j < N; j++)
                  mode(j) = V(Index(j), n);
                
                string name = "ModeN" + to_str(i) + ".dat";
                ofstream file_mode(name.data());
                file_mode.precision(15);
                for (int j = 0; j < N; j++)
                  file_mode << Theta(j) << " " << real(mode(j)) << " " << imag(mode(j)) << '\n';
                
                file_mode.close();
              }
	  }
	
	taux_out.close();
      }

    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Finalize();
    exit(0);
  }
  
  //! performs one iteration of Jacobi step
  void AdvanceJacobi(VectComplexe& g_in, VectComplexe& x_sol, VectComplexe& source_rhs,
		     VectComplexe& g_in_prev)
  {
    Complexe one;
    SetComplexOne(one);
    
    // vecteurs contenant phi et Lphi
    Vector<Complexe> phi(nb_ddl_phi);
    Vector<Complexe> Lphi(nb_ddl_phi);
    if (input_var.type_operator == input_var.OPERATOR_T)
      {
	phi.Reallocate(nb_ddl_gamma);
	Lphi.Reallocate(nb_ddl_gamma);
      }
    else if (input_var.type_operator == input_var.OPERATOR_LOCAL)
      Lphi.Reallocate(nb_ddl_gamma);
        
    // traces de u et phi sur l'interface sigma (valeurs sur les points de quadrature)
    Vector<Complexe> trace_Lphi;
    
    // on part de g_in      
    // on actualise la source qu'on met dans x_sol
    // quelque part on resout -omega^2 u - Laplacien u = f
    // et du/dn - i omega z T(u) = g_in 
    
    // partie surfacique (due a g_in)
    ComputeSourceGin(g_in, x_sol);
    
    // rajout de la partie volumique (due a f)
    Add(one, source_rhs, x_sol);
    
    // resolution des problemes locaux, et calcul de Lphi
    ComputeLocalSolution(x_sol);
    ComputeLphi(x_sol, phi, trace_Lphi, Lphi);
          
    // on stocke la valeur precedent de g_in pour la relaxation
    Copy(g_in, g_in_prev);
    
    // on calcule g_out = -g_in - 2 ik u
    // g_in est ecrase par g_out
    UpdateGin(g_in, x_sol, trace_Lphi, Lphi);
    
    ExchangeGinGout(g_in);
    
    // relaxation
    Real_wp omega_relax = input_var.omega_relax;
    for (int i = 0; i < g_in.GetM(); i++)
      g_in(i) = omega_relax*g_in(i) + (1.0-omega_relax)*g_in_prev(i);      
  }

  //! performs one iteration of Gauss-Seidel algorithm
  void AdvanceGaussSeidel(VectComplexe& g_in, VectComplexe& x_sol, VectComplexe& source_rhs,
			  VectComplexe& g_in_prev)
  {
    Complexe one;
    SetComplexOne(one);
    
    // vecteurs contenant phi et Lphi
    Vector<Complexe> phi(nb_ddl_phi);
    Vector<Complexe> Lphi(nb_ddl_phi);
    if (input_var.type_operator == input_var.OPERATOR_T)
      {
	phi.Reallocate(nb_ddl_gamma);
	Lphi.Reallocate(nb_ddl_gamma);
      }
    else if (input_var.type_operator == input_var.OPERATOR_LOCAL)
      Lphi.Reallocate(nb_ddl_gamma);

    // traces de u et phi sur l'interface sigma (valeurs sur les points de quadrature)
    Vector<Complexe> trace_Lphi;
    Real_wp omega_relax = input_var.omega_relax;
    //if (rank_proc == 1)
    //omega_relax = 1.0;
    
    // boucle sur les processeurs
    for (int n = 0; n < nb_proc; n++)
      {	
	if (rank_proc == n)
	  {
	    // on part de g_in      
	    // on actualise la source qu'on met dans x_sol
	    // quelque part on resout -omega^2 u - Laplacien u = f
	    // et du/dn - i omega z T(u) = g_in 
	    
	    // partie surfacique (due a g_in)
	    ComputeSourceGin(g_in, x_sol);
	  
	    // rajout de la partie volumique (due a f)
	    Add(one, source_rhs, x_sol);
	    
	    // resolution des problemes locaux, et calcul de Lphi
	    ComputeLocalSolution(x_sol);
	    ComputeLphi(x_sol, phi, trace_Lphi, Lphi);
	    
	    // on calcule g_out = -g_in - 2 ik u
	    // g_in est ecrase par g_out
	    // on stocke la valeur precedent de g_in pour la relaxation
	    Copy(g_in, g_in_prev);
	    UpdateGin(g_in, x_sol, trace_Lphi, Lphi);	    
	  }
	
	// on envoie g_out aux autres domaines
	// la relaxation est faite dans la foulee
	ExchangeRelaxGinGout(g_in, omega_relax, n);
	
	// on remplace g_out par le bon g_in
	// pour que la relaxation marche a l'iteration d'apres
	if (rank_proc == n)
	  Copy(g_in_prev, g_in);
      }
  }

  //! applies linear system y = A x, where A is the linear system
  //! to be solved in order to obtain the unknowns g
  template<class Vector1>
  void MltOperator(const Vector1& x, Vector1& y) const
  {
    Complexe one;
    SetComplexOne(one);
    
    // vecteurs contenant phi et Lphi
    Vector<Complexe> phi(nb_ddl_phi);
    Vector<Complexe> Lphi(nb_ddl_phi);
    if (input_var.type_operator == input_var.OPERATOR_T)
      {
	phi.Reallocate(nb_ddl_gamma);
	Lphi.Reallocate(nb_ddl_gamma);
      }
    else if (input_var.type_operator == input_var.OPERATOR_LOCAL)
      Lphi.Reallocate(nb_ddl_gamma);
    
    // traces de u et phi sur l'interface sigma (valeurs sur les points de quadrature)
    Vector<Complexe> trace_Lphi;
    
    // on part de g_in      
    // on actualise la source qu'on met dans x_sol
    // quelque part on resout -omega^2 u - Laplacien u = f
    // et du/dn - i omega z T(u) = g_in 
    
    // partie surfacique (due a g_in)
    ComputeSourceGin(x, x_sol_glob);
    
    // resolution des problemes locaux, et calcul de Lphi
    ComputeLocalSolution(x_sol_glob);
    ComputeLphi(x_sol_glob, phi, trace_Lphi, Lphi);
    
    // on calcule g_out = -g_in - 2 ik u
    // g_in est ecrase par g_out
    Copy(x, y);
    UpdateGin(y, x_sol_glob, trace_Lphi, Lphi);
    
    ExchangeGinGout(y);
    for (int i = 0; i < y.GetM(); i++)
      y(i) = x(i) - y(i);
  }
  
  
  //! applies preconditioning
  void Solve(const VirtualMatrix<Complex_wp>& A, const VectComplex_wp& b, VectComplex_wp& x)
  {
    MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    for (int i = 0; i < Nvol; i++)
      x_sol_glob(i) = OperatorP(i)*b(i);
    
    for (int i = Nvol; i < x_sol_glob.GetM(); i++)
      x_sol_glob(i) = 0;
    
    mat_lu.Solve(x_sol_glob);
    
    var.AddDomains(x_sol_glob);
    for (int i = 0; i < Nvol; i++)
      x(i) = OperatorP(i)*x_sol_glob(i);
  }


  void MltAddVector(const Real_wp& alpha, const Vector<Real_wp>& x,
		    const Real_wp& beta, VectReal_wp& y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }

  void MltAddVector(const Complex_wp& alpha, const Vector<Complex_wp>& x,
		    const Complex_wp& beta, VectComplex_wp& y) const
  {
    VectComplex_wp y2(x);
    this->MltOperator(x, y2);
    
    if (beta == 0)
      y.Fill(0);
    else
      Mlt(beta, y);
    
    Add(alpha, y2, y);
  }

  
  void MltAddVector(const Real_wp& alpha, const SeldonTranspose&, const Vector<Real_wp>& x,
		    const Real_wp& beta, VectReal_wp& y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }

  
  void MltAddVector(const Complex_wp& alpha, const SeldonTranspose&, const Vector<Complex_wp>& x,
		    const Complex_wp& beta, VectComplex_wp& y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }

  
  void MltVector(const VectReal_wp& x, VectReal_wp& y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }
    

  void MltVector(const VectComplex_wp& x, VectComplex_wp& y) const
  {
    this->MltOperator(x, y);
  }

  
  void MltVector(const SeldonTranspose&, const VectReal_wp& x, VectReal_wp& y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }
    

  void MltVector(const SeldonTranspose&, const VectComplex_wp& x, VectComplex_wp& y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
};


template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& input_file, const string& name_element)
{
  typedef typename TypeEquation::Complexe Complexe;
  typedef typename TypeEquation::Dimension Dimension;
  typedef Vector<Complexe> VectComplexe;
  var.InitIndices(100);
      
  // lecture du fichier de donnees
  Vector<string> lines_data_file;
#ifdef SELDON_WITH_MPI
  ReadLinesFile(input_file, lines_data_file, MPI_COMM_WORLD);
#else
  ReadLinesFile(input_file, lines_data_file);
#endif
  
  ReadInputFile(lines_data_file, var);

  // objet pour recuperer le nombre d'iterations max et la tolerance
  All_LinearSolver glob_solver(var);
  ReadInputFile(lines_data_file, glob_solver);
    
  // operateur utilise dans les iterations
  SubdomainMatrixFreeOperator<TypeEquation> glob_mat(var, lines_data_file);
  SubdomainSolverParameter& input_var = glob_mat.input_var;
  
  var.var_chrono.SetMessage("GlobalComputation", "perform the overall simulation");
  var.var_chrono.Start("GlobalComputation");
  
  var.var_chrono.SetMessage("MeshGeneration", "construct and number the mesh");
  var.var_chrono.Start("MeshGeneration");
  
  // we read mesh and construct reference element (for finite element method)
  bool split_mesh = true;
  var.ComputeMeshAndFiniteElement(name_element, split_mesh);
   
  // adding dofs for other variables
  var.PerformOtherInitializations();

  MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
  if (var.print_level >= 0)
    if (rank_proc == 0)
      {
        if (nb_proc > 1)
          {
            cout << "Number of scalar dofs " << var.GetNbGlobalMeshDof() << endl;
            cout << "Number of dofs for the overall problem : " << var.GetNbGlobalDof() << endl;
          }
        else
          {
            cout << "Number of scalar dofs " << mesh_num.GetNbDof() << endl;
            cout << "Number of dofs for the overall problem : " << var.GetNbDof() << endl;
          }
      }
  
  var.var_chrono.Stop("MeshGeneration");
  
  var.var_chrono.SetMessage("JacobianMatricesComputation", "compute jacobian matrices");
  var.var_chrono.Start("JacobianMatricesComputation");

  // computation of Ji, DFi on the mesh
  var.ComputeMassMatrix();

  var.var_chrono.Stop("JacobianMatricesComputation");
  
  // source et solution
  VectComplexe source_rhs, x_sol, r_vec;
  
  // computation of phase for quasi-periodic condition
  var.ComputeQuasiPeriodicPhase();
  
  // calcul de la source
  var.var_chrono.SetMessage("RightHandSideComputation", "compute the right hand side");
  var.var_chrono.Start("RightHandSideComputation");
  
  var.ComputeRightHandSide(source_rhs, false);

  var.var_chrono.Stop("RightHandSideComputation");
    
  // matrice iterative
  var.var_chrono.SetMessage("IterativeMatrix", "compute the iterative matrix");
  var.var_chrono.Start("IterativeMatrix");
  
  GlobalGenericMatrix<Complex_wp> nat_mat;  
  FemMatrixFreeClass<Complex_wp, TypeEquation> free_mat(var), mat_iterative(var);
  var.AddMatrixWithBC(static_cast<FemMatrixFreeClass_Base<Complex_wp>& >(mat_iterative), nat_mat);
  
  // relevement pour la condition de Dirichlet
  var.ApplyDirichletCondition(SeldonNoTrans, mat_iterative, source_rhs);

  var.var_chrono.Stop("IterativeMatrix");
  
  // calcul des points de quadrature sur la surface
  // et la liste des ddls, ListeDof, nb_ddl_gamma
  glob_mat.ComputeSurfaceDofs();
  
  // calcul de l'operateur L et la matrice de masse M
  glob_mat.ComputeNeededOperators();
  
  // resolution du probleme global si on veut comparer avec
  // la solution exacte
  VectComplexe x_sol_ref;
  if (input_var.compare_solution)
    {
      var.var_chrono.SetMessage("GlobalProblem", "factorize global problem");
      var.var_chrono.Start("GlobalProblem");
      
      if (var.IsSymmetricProblem())
        {
          DistributedMatrix<Complexe, Symmetric, ArrayRowSymSparse> mat_direct;
          All_MatrixLU<Complexe> mat_glob_lu;
          var.AddMatrixWithBC(free_mat, mat_direct, nat_mat);
          mat_glob_lu.SetPrintLevel(2);
          mat_glob_lu.Factorize(mat_direct);
          
          int ierr;
          int type = mat_glob_lu.GetInfoFactorization(ierr);
          PrintFactorizationFailed(type, ierr);
          
          x_sol_ref = source_rhs;
          var.AddDomains(x_sol_ref);
          mat_glob_lu.Solve(x_sol_ref);
        }
      else
        {
          DistributedMatrix<Complexe, General, ArrayRowSparse> mat_direct;
          All_MatrixLU<Complexe> mat_glob_lu;
          var.AddMatrixWithBC(free_mat, mat_direct, nat_mat);
          mat_glob_lu.SetPrintLevel(2);
          mat_glob_lu.Factorize(mat_direct);
          
          int ierr;
          int type = mat_glob_lu.GetInfoFactorization(ierr);
          PrintFactorizationFailed(type, ierr);
          
          x_sol_ref = source_rhs;
          var.AddDomains(x_sol_ref);
          mat_glob_lu.Solve(x_sol_ref);
        }
      
      var.var_chrono.Stop("GlobalProblem");
    }
  
  // resolution des problemes locaux
  // et factorisation de la matrice de masse si necessaire
  glob_mat.FactorizeOperators(free_mat, nat_mat);

  // allocation des vecteurs necessaires
  glob_mat.InitVectors();
  int sizeA = glob_mat.GetNbDof();
  int size_gin = glob_mat.GetNbDofGin();
  int Nvol = mesh_num.GetNbDof();
  if (input_var.type_resolution == input_var.PRECOND_SOLVER)
    sizeA = Nvol;
  
  // vecteur d'iteration et le precedent (pour relaxation)
  Vector<Complexe> g_in(size_gin), g_in_prev(size_gin);
  g_in.Fill(0); g_in_prev.Fill(0);
  
  // vecteur distribue solution, residu et la source
  x_sol.Reallocate(sizeA); r_vec.Reallocate(sizeA);
  x_sol.Fill(0); r_vec.Fill(0);
  source_rhs.Resize(sizeA);
  for (int i = Nvol; i < sizeA; i++)
    source_rhs(i) = 0.0;
  
  Vector<Complexe> source_rhs_assembled(source_rhs);
  var.AddDomains(source_rhs_assembled);
  
  DistributedVector<Complexe> b(var.GetOverlappedDofNumber(), var.comm_group_mode);
  DistributedVector<Complexe> r(var.GetOverlappedDofNumber(), var.comm_group_mode);
  
  b.SetData(source_rhs_assembled); r.SetData(r_vec);  
  
  Complexe one;
  SetComplexOne(one);

  // normes absolues des second membres
  // utilisees ensuite pour calculer les erreurs relatives
  int nb_iter = 0, nb_iter_max = glob_solver.GetMaxNumberIteration();
  Real_wp tol = glob_solver.GetStoppingCriterion();
  Real_wp b_norm = Norm2(b), residu = 1.0, x_norm(1);
  if (input_var.compare_solution)
    {
      r.Fill(0);
      for (int i = 0; i < x_sol_ref.GetM(); i++)
        r(i) = x_sol_ref(i);
      
      x_norm = Norm2(r);
    }
  
  ofstream history_residu;
  if (0 == rank_proc)
    {
      history_residu.open("residu.txt");
      history_residu << nb_iter << " 1.0 " << endl;
    }
  
  var.var_chrono.SetMessage("SolutionComputation", "compute the solution of linear system");
  var.var_chrono.SetMessage("Output", "write the solution on asked outputs");
  var.var_chrono.Start("SolutionComputation");
  
  // calcul de la matrice d'iteration si demande
  glob_mat.ComputeJacobiIterationMatrix(g_in, x_sol);
  
  // boucle principale pour Jacobi ou Gauss-Seidel
  if (input_var.type_resolution == input_var.GMRES_SOLVER)
    {
      input_var.omega_relax = 1.0;
      Vector<Complexe> rhs_g_in(g_in.GetM());
      rhs_g_in.Fill(0);
      
      // declaration des vecteurs distribues pour le Gmres
      IVect rows;
      DistributedVector<Complexe> g_in_dist(rows, MPI_COMM_WORLD);
      DistributedVector<Complexe> rhs_g_in_dist(rows, MPI_COMM_WORLD);
      g_in_dist.SetData(g_in);
      rhs_g_in_dist.SetData(rhs_g_in);

      // calcul du second membre avec une iteration de Jacobi
      glob_mat.AdvanceJacobi(rhs_g_in, x_sol, source_rhs, g_in_prev);
            
      /* string name = "mat_gmres_P"+to_str(rank_proc)+".dat";
      ofstream file_out_gmres(name.data());
      file_out_gmres.precision(15);
      
      for (int j = 0; j < 2*g_in.GetM(); j++)
	{
	  // pour calculer la matrice d'iteration
	  // on fait le produit avec e_i
	  g_in.Fill(0);
	  if (j >= g_in.GetM())
	    {
	      if (rank_proc == 1)
		g_in(j-g_in.GetM()) = 1.0;	      
	    }
	  else
	    {
	      if (rank_proc == 0)
		g_in(j) = 1.0;
	    }
	  
	  Mlt(glob_mat, g_in, g_in_prev);
	  
	  // on stocke la colonne de la matrice d'iteration
	  for (int i = 0; i < g_in_prev.GetM(); i++)
	    if (g_in_prev(i) != Complex_wp(0, 0))
	      {
		if (rank_proc == 0)
		  file_out_gmres << i+1 << " " << j+1 << " " << g_in_prev(i) << '\n';
		else
		  file_out_gmres << i+g_in.GetM()+1 << " " << j+1 << " " << g_in_prev(i) << '\n';
	      }
	}

	file_out_gmres.close(); */
      
      Iteration<Real_wp> iter(nb_iter_max, tol);
      if (rank_proc == 0)
	{
	  iter.ShowFullHistory();
	  iter.SaveFullHistory("residu_gmres.txt");
	}
      else
	iter.HideMessages();
      
      iter.SetInitGuess(true);
      iter.SetRestart(glob_solver.GetRestart());
      Preconditioner_Base<Complex_wp> precond; // preconditioneur identite
      Gmres(glob_mat, g_in_dist, rhs_g_in_dist, precond, iter);
      
      // on remplit x_sol avec un coup de Jacobi
      glob_mat.AdvanceJacobi(g_in, x_sol, source_rhs, g_in_prev);

      if (input_var.compare_solution)
	{
	  r.Fill(0);
	  for (int i = 0; i < x_sol_ref.GetM(); i++)
	    r(i) = x_sol(i) - x_sol_ref(i);
	  
          Real_wp err_inv = Norm2(r) / x_norm;
	  if (0 == rank_proc)
            cout << "Erreur inverse = " <<  err_inv << endl;
	}
      
      g_in_dist.Nullify(); rhs_g_in_dist.Nullify();
    }
  else if (input_var.type_resolution == input_var.PRECOND_SOLVER)
    {
      DistributedVector<Complexe> x_sol_dist(var.GetOverlappedDofNumber(), var.comm_group_mode);      
      x_sol_dist.SetData(x_sol);

      x_sol.Resize(Nvol);
      b.Resize(Nvol);

      Iteration<Real_wp> iter(nb_iter_max, tol);
      if (rank_proc == 0)
	{
	  iter.ShowFullHistory();
	  //iter.SaveFullHistory("residu_precond.txt");
	}
      else
	iter.HideMessages();
      
      iter.SetInitGuess(true);
      iter.SetRestart(glob_solver.GetRestart());
      
      CoCg(mat_iterative, x_sol_dist, b, glob_mat, iter);
      
      if (input_var.compare_solution)
	{
	  r.Fill(0);
	  for (int i = 0; i < x_sol_ref.GetM(); i++)
	    r(i) = x_sol(i) - x_sol_ref(i);
	  
	  cout << "Erreur inverse = " << Norm2(r) / x_norm << endl;
	}
      
      x_sol_dist.Nullify();
    }
  else
    {
      Real_wp residu_prec = 1.0;
      while ( (nb_iter < nb_iter_max) && (residu > tol))
	{
	  // on calcule le nouvel itere
	  if (input_var.type_resolution == input_var.GAUSS_SEIDEL_SOLVER)
	    glob_mat.AdvanceGaussSeidel(g_in, x_sol, source_rhs, g_in_prev);
	  else
	    glob_mat.AdvanceJacobi(g_in, x_sol, source_rhs, g_in_prev);
	  
	  // calcul du residu r = b - A x
	  if (input_var.compare_solution)
	    {
	      // erreur inverse : on calcule x - x_ref
	      r.Fill(0);
	      for (int i = 0; i < x_sol_ref.GetM(); i++)
		r(i) = x_sol(i) - x_sol_ref(i);
	      
	      residu_prec = residu;
	      residu = Norm2(r) / x_norm;
	    }
	  else
	    {
	      // erreur directe : on calcule r = b - A x
	      Copy(b, r);
	      //x_sol.Write("xsol"+to_str(rank_proc)+".dat");
	      mat_iterative.MltAddVector(-one, x_sol, one, r_vec);
	      //r_vec.Write("rvec"+to_str(rank_proc)+".dat");
	      residu_prec = residu;
	      residu = Norm2(r)/b_norm;      
	      //DISP(r_vec);
	    }
	  
	  if (rank_proc == 0)
	    {
	      cout << "Residu at iteration " << nb_iter << " = " << residu << endl;      
	      history_residu << nb_iter+1 << " " << residu << endl;
	      
	      cout << "Taux = " << residu/residu_prec << endl;
	    }              
	  
	  nb_iter++;
	}
      
      //g_in.Write("Gin_ref_P"+to_str(rank_proc)+".dat");
    }
  
  if (rank_proc == 0)
    history_residu.close();
  
  var.var_chrono.Stop("SolutionComputation");

  // on remet Dirichlet
  mat_iterative.ImposeDirichletCondition(SeldonNoTrans, x_sol);
  x_sol.Write("sol_sub_P"+to_str(rank_proc)+".dat");
  
  // on ecrit la solution sur le disque
  var.var_chrono.Start("Output");
  var.WriteDatas(x_sol);
  
  var.var_chrono.Stop("Output");  
  var.var_chrono.Stop("GlobalComputation");
  var.var_chrono.DisplayAll();
  
  b.Nullify(); r.Nullify();
  glob_mat.g_in_quad.Nullify();
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
#ifdef SELDON_WITH_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
      MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#endif
      
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
            
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      if (rank_proc == 0)
        cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;
      
      if (rank_proc == 0)
        cout<<"Helmholtz Solver with Lobatto quadrilaterals "<<endl; 
      
      EllipticProblem<HelmholtzEquation<Dimension2> > Vars;
      //EllipticProblem<HelmholtzEquation<Dimension3> > Vars;
      RunAll(Vars, file_name_data, type_element);
      
      if (rank_proc == 0)
        cout<<" we destroy the variables "<<endl;

    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"helmholtz2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  if (rank_proc == 0)
    cout<<"End of the program"<<endl; 
  
  return FinalizeMontjoie();
}
