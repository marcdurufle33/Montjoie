#ifndef MONTJOIE_FILE_NEWTON_SOLVER_CXX

namespace Montjoie
{

  //! default constructor
  template<class T, class VectorSol, class VectorScale>
  NewtonSolver<T, VectorSol, VectorScale>::NewtonSolver()
  {
    tol = 1e4*GetPrecisionMachine(Treal(1));
    tol_convergence = 10.0*tol;
    nbitermax = 5;
    erreur = 0.0;
    scale_newton_test = false;
    force_reevaluation_jacobian = false;
    threshold_norme_scheme = 1e-5;
    print_level = 0;
    backtracking = false;
    decreasing_residual_required = true;
    convergence_required = true;
  }


  //! Returns the memory used by the object
  template<class T, class VectorSol, class VectorScale>
  size_t NewtonSolver<T, VectorSol, VectorScale>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille += x_memory.GetMemorySize() + x_prec.GetMemorySize();
    taille += scheme.GetMemorySize() + save_scheme.GetMemorySize();
    taille += scale_newton.GetMemorySize();
    return taille;
  }
  
  
  //! modification of parameters of the data file
  template<class T, class VectorSol, class VectorScale>
  void NewtonSolver<T, VectorSol, VectorScale>
  ::SetInputData(const string & description_field, const Vector<string> &parameters)
  {
    if (!description_field.compare("NewtonParam"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of NewtonSolver" << endl;
	    cout << "NewtonParam needs at least two parameters, for instance :" << endl;
	    cout << "NewtonParam = tol nb_iter_max" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	tol = to_num<Treal>(parameters(0));
	nbitermax = to_num<int>(parameters(1));
	if (parameters.GetM() < 3)
	  tol_convergence = tol;
	else
	  tol_convergence = to_num<Treal>(parameters(2));
	
	if (parameters.GetM() < 4)
	  scale_newton_test = true;
	else
	  {
	    if (parameters(3) == "NO_SCALE")
	      scale_newton_test = false;
	  }
	
	if (parameters.GetM() >= 5)
	  threshold_norme_scheme = to_num<Treal>(parameters(4));
      }
    else if (!description_field.compare("ForceReevaluationJacobian"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of MegaNewtonSolver" << endl;
	    cout << "ForceReevaluationJacobian needs one parameter, for instance :" << endl;
	    cout << "ForceReevaluationJacobian = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	if (parameters(0) == "YES")
	  force_reevaluation_jacobian = true;
	else
	  force_reevaluation_jacobian = false;
      }
    
  }
  
  
  //! initializing the scheme
  template<class T, class VectorSol, class VectorScale>
  void NewtonSolver<T, VectorSol, VectorScale>
  ::Init(NonLinearEquations_Newton<T, VectorSol, VectorScale>& var, VectorSol& x_init)
  {
    scheme = x_init;
    x_memory = x_init;
    x_prec = x_init;
    save_scheme = x_init;
    
    scheme.Zero();
    x_memory.Zero();
    x_prec.Zero();
    save_scheme.Zero();

    // first jacobian matrix is computed here
    var.ComputeAndFactoriseDiff(x_init, scale_newton);
    
    // fichier stockant les residus
    if (print_level >= 10)
      file_out_residual_newton.open("residu_newton.txt");
  }
  
  
  //! Solves F(x) = 0 and returns the number of iterations
  template<class T, class VectorSol, class VectorScale>
  int NewtonSolver<T, VectorSol, VectorScale>
  ::Solve(NonLinearEquations_Newton<T, VectorSol, VectorScale>& var, VectorSol& x)
  {      
    erreur = 1.0;
    int nb = 0;	  
    
    // F(x) is evaluated and stored in scheme
    var.ComputeScheme(x, scheme);
    if (scale_newton_test)
      scheme = scheme*scale_newton;
        
    // norme_ref is used to compute the relative residue
    Real_wp norme_ref = max(threshold_norme_scheme, var.GetNorm2Vector(scheme));

    // previous values of x and scheme are stored
    x_memory = x;
    save_scheme = scheme;
    
    // a first try is made to compute the solution with Newton solver
    Newton(var, x, nb, norme_ref, false, true);
    
    // if nb = nbitermax, the algorithm failed to converge, we try again 
    // by recomputing the differential matrix
    if (nb == nbitermax)
      {
	x = x_memory;
	scheme = save_scheme;
        
	var.ComputeAndFactoriseDiff(x, scale_newton);
        
	nb = 0;
	Newton(var, x, nb, norme_ref);
      }      
    
    // the residue is printed
    if (print_level >= 10)
      file_out_residual_newton << nb << " " << erreur << " " << norme_ref << endl;
    
    if ((erreur > tol_convergence*norme_ref) && (convergence_required))
      {
	cout << "Newton failed to give the solution"<<endl;
	DISP(norme_ref); DISP(erreur/norme_ref);
	DISP(tol_convergence);
	DISP(nbitermax);
	abort();
      }
    
    return nb;
  }
  
  
  //! Basic Newton algorithm with backtracking
  template<class T, class VectorSol, class VectorScale>
  void NewtonSolver<T, VectorSol, VectorScale>
  ::Newton(NonLinearEquations_Newton<T, VectorSol, VectorScale>& var,
	   VectorSol& x, int &nb, const Treal& norme_ref,
	   bool reevaluate_jacobian, bool scheme_already_evaluated)
  {
    if (!scheme_already_evaluated)
      {
	var.ComputeScheme(x, scheme);
	if (scale_newton_test)
	  scheme = scheme*scale_newton;
      }
    
    // the residue is computed
    erreur = var.GetNorm2Vector(scheme);
    
    // previous residue is initialized to a larger value
    Treal err_prec = erreur*2.0;
    if (print_level >= 2)
      {
	cout << "Residue in Newton solver at iteration " << nb << " = " << erreur << endl;
      }
    
    // on verifie que le newton ne stagne pas, ou que le residu est assez petit
    while ( ((!decreasing_residual_required) 
	     || (erreur < err_prec)) && (erreur > tol*norme_ref)  && (nb < nbitermax))
      {
	nb++;
	err_prec = erreur;
        
	// we compute the new step DF^-1 F(x)
	if (force_reevaluation_jacobian)
	  reevaluate_jacobian = true;
	
	// nb > 1 because ComputeAndFactoriseDiff has been called in Init
	if ((reevaluate_jacobian) && (nb > 1))
	  var.ComputeAndFactoriseDiff(x, scale_newton);
	
	x_prec = scheme;
	var.SolveDifferential(x_prec, scheme);
	x_prec = x;
        
	// x is updated
	x -= scheme;
	
	// F(x) is evaluated to test the convergence and for the next iteration
	var.ComputeScheme(x, scheme);
	if (scale_newton_test)
	  scheme = scheme*scale_newton;
	
	erreur = var.GetNorm2Vector(scheme);
	
	if ((erreur > err_prec) && (backtracking))
	  {
	    // trying backtracking method to ensure a decreasing residue
	    int lvl_max = 14, lvl = 0;
	    Treal alpha = 0.5;
	    save_scheme = x;
	    while ( (erreur > err_prec) && (lvl < lvl_max))
	      {
		x = (Treal(1) - alpha) * x_prec + alpha*save_scheme;
		
		var.ComputeScheme(x, scheme);
		if (scale_newton_test)
		  scheme = scheme*scale_newton;
		
		erreur = var.GetNorm2Vector(scheme);
		alpha *= 0.5; lvl++;
	      }
	  }
	
	if (print_level >= 2)
	  cout << "Residue in Newton solver at iteration " << nb << " = " << erreur << endl;
      }
    
    // If the error increased, we take the previous step
    if (erreur > err_prec)
      x = x_prec;
  }
  
}

#define MONTJOIE_FILE_NEWTON_SOLVER_CXX
#endif
