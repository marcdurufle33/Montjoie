#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

class InputDispersionVariables : public InputDataProblem_Base
{
public :
  
  // maximal k used to determine Taylor expansion of dispersion error
  Real_wp kmax;
  
  // threshold used to replace floats by fractions in Taylor expansion
  Real_wp threshold_rational;
  
  // two wave numbers used to determine the order of approximation
  Real_wp kh1_order, kh2_order;
  
  // ponderation between Gauss (exact integration) and Gauss-Lobatto
  Real_wp alpha_lumping;
  
  // two extremities of the dispersion curve
  Real_wp kmin_courbe, kmax_courbe;
  
  // order of approximation
  int order;
  
  // number of points in the dispersion curve
  int nb_points_courbe;
  
  // verbose level
  int print_level;
  
  // LOBATTO or GAUSS
  int type_quadrature;
  
  // if true, there is a ponderation between Gauss and Gauss-Lobatto
  bool ponderation_lumping;
  
  // if true, Taylor expansion is searched
  bool search_taylor_dvt;
  
  // if true, the order of dispersion error is searched
  bool find_order;
  
  // if true, we try to compute a dispersion curve
  bool courbe_dispersion;
  
  // if true wave equation, otherwise advection equation
  bool wave_equation;
  
  // if true, local refinement is applied to avoid intersections between eigenvalues
  bool local_refinement;
  
  // if true we consider a Discontinuous Galerkin formulation
  bool dg_formulation;
  
  // order of dispersion found
  int order_dispersion;
  
  // maximal order for Taylor expansion
  int order_dvt_taylor;
  
  // maximal wave number for which no eigenvalues are intersecting
  Real_wp kmin_no_spurious;
  
  // penalization parameters
  Real_wp alpha_dg, delta_dg;
  
  InputDispersionVariables()
  {
    order = 1;
    nb_points_courbe = 2;
    print_level = 1;
    type_quadrature = Globatto<Real_wp>::QUADRATURE_GAUSS;
    ponderation_lumping = false;
    search_taylor_dvt = false;
    find_order = true;
    courbe_dispersion = false;
    wave_equation = false;
    local_refinement = false;
    dg_formulation = true;
    kmin_no_spurious = 1;
    alpha_dg = 0;
    delta_dg = 0;
  }
  
  // setting parameters with input file
  void SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("OrderDiscretization"))
      {
	to_num(parameters(0), order);
      }
    else if (!description_field.compare("FormulationDG"))
      {
	if (!parameters(0).compare("YES"))
	  dg_formulation = true;
	else
	  dg_formulation = false;
      }
    else if (!description_field.compare("PenalizationDG"))
      {
	alpha_dg = to_num<Real_wp>(parameters(0));
	delta_dg = to_num<Real_wp>(parameters(1));
      }
    else if (!description_field.compare("WaveEquation"))
      {
	if (!parameters(0).compare("YES"))
	  wave_equation = true;
	else
	  wave_equation = false;
      }
    else if (!description_field.compare("LocalRefinement"))
      {
	if (!parameters(0).compare("YES"))
	  local_refinement = true;
	else
	  local_refinement = false;
      }
    else if (!description_field.compare("TypeQuadrature"))
      {
	if (!parameters(0).compare("LOBATTO"))
	  type_quadrature = Globatto<Real_wp>::QUADRATURE_LOBATTO;
	else
	  type_quadrature = Globatto<Real_wp>::QUADRATURE_GAUSS;
	
	// DISP(type_quadrature);
      }
    else if (!description_field.compare("MaxWaveNumber"))
      {
	to_num(parameters(0), kmax);
      }
    else if (!description_field.compare("ThresholdRational"))
      {
	to_num(parameters(0), threshold_rational);
	// DISP(threshold_rational);
      }
    else if (!description_field.compare("DisplayGraph"))
      {
	if (!parameters(0).compare("YES"))
	  {
	    courbe_dispersion = true;
	    to_num(parameters(1), kmin_courbe);
	    to_num(parameters(2), kmax_courbe);
	    to_num(parameters(3), nb_points_courbe);
	  }
	else
	  courbe_dispersion = false;
      }
    else if (!description_field.compare("TaylorDevelopment"))
      {
	if (!parameters(0).compare("YES"))
	  {
	    search_taylor_dvt = true;
	    to_num(parameters(1), order_dvt_taylor);
	  }
	else
	  search_taylor_dvt = false;
      }
    else if (!description_field.compare("PrintLevel"))
      {
	to_num(parameters(0), print_level);
      }
    else if (!description_field.compare("FindOrder"))
      {
	if (!parameters(0).compare("YES"))
	  {
	    find_order = true;
	    to_num(parameters(1), kh1_order);
	    to_num(parameters(2), kh2_order);
	  }
	else
	  find_order = false;
      }
    else if (!description_field.compare("Ponderation"))
      {
	if (!parameters(0).compare("YES"))
	  {
	    ponderation_lumping = true;
	    int num, denom;
	    to_num(parameters(1), num);
	    to_num(parameters(2), denom);
	    alpha_lumping = Real_wp(num)/Real_wp(denom);
	  }
	else
	  ponderation_lumping = false;
      }
  }

};


// basic class to store intermediary arrays needed for evaluate the dispersion error
class DispersionSolver : public InputDispersionVariables
{
public :
  //! Relation between periodic dofs and originals dofs
  /*!
    PeriodicityDof(i, 0) = number of the original dof obtained by translation of periodic dof i
    PeriodicityDof(i, 1:dim_N) = translation vector between the two dofs
  */
  Matrix<int> PeriodicityDof;
  //! list of original dof numbers
  IVect Num_OriginalDofs;
  //! true is the dof is periodic
  VectBool IsDofPeriodic;
  
  VectComplex_wp phase;
  Matrix<Complex_wp, General, RowMajor> Kh_per, Dh_per, Ah_per;


  // constructs informations on periodicity of dofs of the mesh
  void GetPeriodicityDofs()
  {
    int r = order;
    
    if (dg_formulation)
      {
	if (wave_equation)
	  {
	    // 2(r+1) degrees of freedom on the central element
	    int nb_original_dof = 2*(r+1);
	    Num_OriginalDofs.Reallocate(2*(r+1)); Num_OriginalDofs.Fill();
	    // we consider three elements (left, central, and right)
	    IsDofPeriodic.Reallocate(3*nb_original_dof); 
	    IsDofPeriodic.Fill(false);
	    PeriodicityDof.Reallocate(6*(r+1),2);
	    for (int i = 0; i <= order; i++)
	      {
		// first we number central element
		PeriodicityDof(i, 0) = i; PeriodicityDof(i, 1) = 0;
		PeriodicityDof(r+1+i, 0) = r+1+i; PeriodicityDof(r+1+i, 1) = 0;
		
		// then left element
		PeriodicityDof(2*(r+1)+i, 0) = i; PeriodicityDof(2*(r+1)+i,1) = -1;
		PeriodicityDof(3*(r+1)+i, 0) = r+1+i; PeriodicityDof(3*(r+1)+i,1) = -1;
		
		// then right element
		PeriodicityDof(4*(r+1)+i,0) = i; PeriodicityDof(4*(r+1)+i,1) = 1;
		PeriodicityDof(5*(r+1)+i,0) = r+1+i; PeriodicityDof(5*(r+1)+i,1) = 1;
		
		// only central element is in the original pattern
		IsDofPeriodic(2*(r+1)+i) = true; IsDofPeriodic(3*(r+1)+i) = true;
		IsDofPeriodic(4*(r+1)+i) = true; IsDofPeriodic(5*(r+1)+i) = true;
	      }
	  }
	else
	  {
	    int nb_original_dof = r+1;
	    Num_OriginalDofs.Reallocate(r+1);
	    Num_OriginalDofs.Fill();
	    IsDofPeriodic.Reallocate(3*nb_original_dof);
	    IsDofPeriodic.Fill(false);
	    PeriodicityDof.Reallocate(3*(r+1),2);
	    for (int i = 0; i <= order; i++)
	      {
		PeriodicityDof(i, 0) = i; PeriodicityDof(i,1) = 0;
		PeriodicityDof(r+1+i, 0) = i; PeriodicityDof(r+1+i,1) = -1;
		PeriodicityDof(2*(r+1)+i,0) = i; PeriodicityDof(2*(r+1)+i,1) = 1;
		IsDofPeriodic(r+1+i) = true;
		IsDofPeriodic(2*(r+1)+i) = true;
	      }
	  }
      }
    else
      {
	// continuous elements
	// we first number vertices
	Num_OriginalDofs.Reallocate(r); Num_OriginalDofs(0) = 1;
	IsDofPeriodic.Reallocate(2*r+1); IsDofPeriodic.Fill(false);
	PeriodicityDof.Reallocate(2*r+1,2); 
	PeriodicityDof(0,0) = 0; PeriodicityDof(0,1) = -1;
	PeriodicityDof(1,0) = 0; PeriodicityDof(1,1) = 0;
	PeriodicityDof(2,0) = 0; PeriodicityDof(2,1) = 1;
	IsDofPeriodic(0) = true; IsDofPeriodic(2) = true;
	// then dofs inside the central element and left element
	for (int i = 1; i < r; i++)
	  {
	    PeriodicityDof(2+i,0) = i; PeriodicityDof(2+i,1) = -1;
	    PeriodicityDof(1+i+r,0) = i; PeriodicityDof(1+i+r,1) = 0; 
	    Num_OriginalDofs(i) = 1+i+r; IsDofPeriodic(2+i) = true;
	  }
      }
    
    DISP(PeriodicityDof); DISP(Num_OriginalDofs); DISP(IsDofPeriodic);
  }
  
  
  // kwave : wave number k
  // Dh, Kh : mass and stiffness matrix (for all the elements)
  // error_dispersion : error of dispersion (omega - kwave)
  // omega_numer : omega found
  // eigen_values, eigen_vectors : eigenvalues and eigenvectors
  void GetDispersionError(const Real_wp& kwave, Matrix<Real_wp, General, ArrayRowSparse>& Dh,
			  Matrix<Real_wp, General, ArrayRowSparse>& Kh,			
			  Real_wp& error_dispersion, Real_wp& omega_numer,
			  VectComplex_wp& eigen_values, Matrix<Complex_wp>& eigen_vectors)
  {
    VectComplex_wp Xh;
    IVect ipivot;
    Real_wp omega, arg, coef, err_minimal;
    omega = abs(kwave);
    if ((wave_equation)&&(!dg_formulation))
      omega *= omega;
    
    // phase for periodic dofs
    int nb_original_dof = Num_OriginalDofs.GetM();
    phase.Reallocate(IsDofPeriodic.GetM());
    phase.Fill(Complex_wp(1,0));
    for (int i = 0; i < phase.GetM(); i++)
      if (IsDofPeriodic(i))
	{
	  arg = PeriodicityDof(i,1)*kwave;
	  phase(i) = Complex_wp(cos(arg), sin(arg));
	}
    
    // periodic matrices
    Dh_per.Reallocate(nb_original_dof, nb_original_dof);
    Dh_per.Fill(Complex_wp(0,0));
    for (int i1 = 0; i1 < nb_original_dof; i1++)
      {
	int i = Num_OriginalDofs(i1);
	for (int j1 = 0; j1 < Dh.GetRowSize(i); j1++)
	  {
	    int j = Dh.Index(i, j1);
	    int num_dof = PeriodicityDof(j,0);
	    Dh_per(i1, num_dof) += phase(j)*Dh.Value(i,j1);
	  }
      }
    
    Kh_per.Reallocate(nb_original_dof, nb_original_dof);
    Kh_per.Fill(Complex_wp(0,0));
    for (int i1 = 0; i1 < nb_original_dof; i1++)
      {
	int i = Num_OriginalDofs(i1);
	for (int j1 = 0; j1 < Kh.GetRowSize(i); j1++)
	  {
	    int j = Kh.Index(i, j1);
	    int num_dof = PeriodicityDof(j,0);
	    Kh_per(i1, num_dof) += phase(j)*Kh.Value(i,j1);
	  }
      }
    
    // now we compute Dh_per^{-1} Kh_per
    Xh.Reallocate(nb_original_dof);
    Ah_per.Reallocate(nb_original_dof, nb_original_dof);
    GetLU(Dh_per, ipivot); // DISP(Dh_per);
    for (int i = 0; i < nb_original_dof; i++)
      {
	for (int j = 0; j < nb_original_dof; j++)
	  Xh(j) = Kh_per(j,i);
	
	SolveLU(Dh_per, ipivot, Xh);
	
	for (int j = 0; j < nb_original_dof; j++)
	  Ah_per(j, i) = Xh(j);
      }
    
    // we get eigenvalues and eigen-vectors
    GetEigenvaluesEigenvectors(Ah_per, eigen_values, eigen_vectors);
    DISP(eigen_values);
    err_minimal = Real_wp(100); error_dispersion = Real_wp(1);
    omega_numer = Real_wp(0);
    // we search the closest eigenvalue to thre relation dispersion
    for (int i = 0; i < nb_original_dof; i++)
      {
	arg = abs(eigen_values(i));
	coef = abs(omega-arg);
	if (coef < err_minimal)
	  {
	    error_dispersion = arg-omega;
	    err_minimal = coef;
	    if ((wave_equation)&&(!dg_formulation))
	      omega_numer = sqrt(arg);
	    else
	      omega_numer = arg;
	  }
	// DISP(arg/omega2); DISP(error_dispersion);
      }
  }
  
  
  // main method
  template<class Matrix1>
  void GetDispersionRelation(Matrix1& Dh, Matrix1& Kh, VectReal_wp& xdof)
  {
    // Real_wp arg, coef;
    int nodl = 2*order+1;
    if (dg_formulation)
      {
	nodl = 3*(order+1);
	if (wave_equation)
	  nodl *= 2;
      }
    
    Real_wp kwave;
    int nb_original_dof = Num_OriginalDofs.GetM();
    Matrix<Complex_wp> eigen_vectors(nb_original_dof, nb_original_dof);
    
    Real_wp error_dispersion;
    VectComplex_wp eigen_values(nb_original_dof);
    eigen_values.Fill(Complex_wp(0,0));
    cout<<" Number of original dofs " << nb_original_dof << endl;
    Real_wp err1(0), err2(0), omega_numer;
    if (find_order)
      {
	// searching order by estimating slope of the log-log curve
	// with only two points kh1_order and kh2_order
	kwave = kh1_order;
	GetDispersionError(kwave, Dh, Kh, err1, omega_numer,
			   eigen_values, eigen_vectors);
	
	err1 = abs(err1); DISP(kwave); DISP(err1);
	kwave = kh2_order;
	GetDispersionError(kwave, Dh, Kh, err2, omega_numer,
			   eigen_values, eigen_vectors);
	
	err2 = abs(err2); DISP(kwave); DISP(err2);
	error_dispersion = (log(err2)-log(err1))/(log(kh2_order)-log(kh1_order));
	order_dispersion = toInteger(round(error_dispersion));
	cout<<"We found a dispersion order of "<<(error_dispersion-1)<<endl;
	
	Real_wp coef = err2/(pow(kwave,error_dispersion));
	DISP(Real_wp(1)/Real_wp(coef));
      }
    else
      {
	order_dispersion = 2*(order+1);
      }
  
    if (courbe_dispersion)
      {
	// computation of a dispersion curve
	// we display all the eigenvalues (even spurious) in file root.dat
	VectReal_wp coef_relation(nb_original_dof), val_propre(nb_original_dof);
	ofstream file_out("dispersion.dat"); file_out.precision(15);
	ofstream file_err("error.dat"); file_err.precision(15);
	ofstream multi_root("root.dat"); multi_root.precision(15);
	
	// extrapolation used to find an estimate of the next point of the curve
	int order_extrapolation = 2;
	Globatto<Real_wp> extrapol;
	extrapol.ConstructQuadrature(order_extrapolation, extrapol.QUADRATURE_LOBATTO);
	Matrix<Real_wp> Vh(order_extrapolation+1, nb_original_dof);
	Vh.Fill(Real_wp(0));
	VectReal_wp xsi_history(Vh.GetM());
	xsi_history.Fill(Real_wp(0));
	// estimation of the next eigenvalue on the curve
	VectReal_wp coef_estimation(order_extrapolation+1), val_estimate(nb_original_dof);
	// errors on the eigenvectors
	VectReal_wp err_propre(nb_original_dof), old_err_propre(nb_original_dof);
	
	Real_wp dx_global = Real_wp(1)/Real_wp(nb_points_courbe-1), xsi(0);
	Real_wp dx = dx_global, dx_min = 1e-10; int nb_iter = 0;
	while (xsi < 1)
	  {
	    kwave = (Real_wp(1)-xsi)*kmin_courbe + xsi*kmax_courbe ;
	    VectComplex_wp onde_inc(nb_original_dof);
	    for (int j = 0; j < nb_original_dof; j++)
	      onde_inc(j) = exp(-Iwp*kwave*xdof(j));
	    
	    // computation of numerical eigenvalues
	    GetDispersionError(kwave, Dh, Kh, err1, omega_numer, 
			       eigen_values, eigen_vectors);
	    
	    multi_root<<kwave<<" ";
	    file_err<<kwave<<" ";
	    if ((wave_equation)&&(!dg_formulation))
	      for (int j = 0; j < eigen_values.GetM(); j++)
		coef_relation(j) = sqrt(abs(eigen_values(j)));
	    else
	      for (int j = 0; j < eigen_values.GetM(); j++)
		coef_relation(j) = -imag(eigen_values(j));
	    
	    // coef_relation is equal to the pulsation for all eigenvalues
	    
	    // we compute the error made on the physical eigenvector exp(ikx), for each numerical eigenvector
	    for (int i = 0; i < eigen_values.GetM(); i++)
	      {
		old_err_propre(i) = 0.0;
		for (int j = 0; j < nb_original_dof; j++)
		  old_err_propre(i) += square(abs(onde_inc(j) - eigen_vectors(i,j)/eigen_vectors(i,0)));
		
		old_err_propre(i) = sqrt(old_err_propre(i));
	      }
	    // DISP(old_err_propre); // exit(0);
	    
	    file_out<<kwave<<" "<<omega_numer<<endl;
	  
	    if (nb_iter > 0)
	      {
		if ((nb_iter > order_extrapolation)&&(kwave > kmin_no_spurious))
		  {
		    // on calcule l'extrapolee de chaque valeur propre
		    extrapol.AffectPoints(xsi_history);
		    val_estimate.Fill(Real_wp(0));
		    // DISP(xsi_history); DISP(xsi); DISP(extrapol.Cte_of_G); DISP(extrapol.Points);
		    for (int k = 0; k <= order_extrapolation; k++)
		      coef_estimation(k) = extrapol.EvaluatePhi(k, xsi);
		    
		    // DISP(coef_estimation);
		    for (int j = 0; j < nb_original_dof; j++)
		      {
			for (int k = 0; k <= order_extrapolation; k++)
			  val_estimate(j) += coef_estimation(k)*Vh(k, j);
		      }
		  }
		else
		  for (int j = 0; j < nb_original_dof; j++)
		    val_estimate(j) = Vh(order_extrapolation,j);
		
		// val_estimate is the estimation of the next eigenvalues on each curve
		// DISP(xsi_history); DISP(xsi); DISP(Vh); DISP(val_estimate);
		
		VectBool EstimateUsed(nb_original_dof), RelationUsed(nb_original_dof);
		EstimateUsed.Fill(false); RelationUsed.Fill(false);
		for (int p = 0; p < nb_original_dof; p++)
		  {
		    // pour chaque valeur, on trouve le point le plus proche de l'extrapolation
		    VectReal_wp proximite_valeur(nb_original_dof); IVect index_val(nb_original_dof);
		    proximite_valeur.Fill(1e300); index_val.Fill(-1);
		    for (int j = 0; j < nb_original_dof; j++)
		      if (!EstimateUsed(j))
			{
			  for (int k = 0; k < nb_original_dof; k++)
			    if (!RelationUsed(k))
			      {
				Real_wp dist = abs(coef_relation(k) - val_estimate(j));
				// if ((dist < err_min)&&(!IndexUsed(k)))
				if (dist < proximite_valeur(j))
				  {
				    proximite_valeur(j) = dist;
				    index_val(j) = k;
				  }
			      }
			}
		    
		    // DISP(p); DISP(proximite_valeur); DISP(index_val);
		    // on boucle sur toutes les distances pour choisir la plus petite
		    Real_wp dmin(1e300); int index = -1;
		    for (int k = 0; k < nb_original_dof; k++)
		      if (proximite_valeur(k) < dmin)
			{
			  dmin = proximite_valeur(k);
			  index = k;
			}
		    
		    EstimateUsed(index) = true;
		    RelationUsed(index_val(index)) = true;
		    val_propre(index) = coef_relation(index_val(index));
		    err_propre(index) = old_err_propre(index_val(index));
		  }
		
		if (local_refinement)
		  {
		    // maintenant, on adapte le pas, si des valeurs propres sont trop proches
		    Real_wp dmin(100);
		    for (int j = 0; j < nb_original_dof; j++)
		      for (int k = 0; k < nb_original_dof; k++)
			if (j != k)
			  dmin = min(dmin, 0.5*abs(val_propre(j)-val_propre(k)));
		    
		    dx = min(dmin, dx_global);
		    
		    // plus sioux, on adapte le pas en prévision des intersections
		    int N = Vh.GetM()-1;
		    Real_wp x0, x1, y0, y1, z0, z1, xinter, delta;
		    for (int j = 0; j < nb_original_dof; j++)
		      for (int k = 0; k < nb_original_dof; k++)
			if (j != k)
			  {
			    x0 = xsi_history(N); x1 = xsi;
			    y0 = Vh(N,j); y1 = val_propre(j);
			    z0 = Vh(N,k); z1 = val_propre(k);
			    delta = z1-z0-y1+y0;
			    if (abs(delta) > 1e-10)
			      {
				// possible intersection
				xinter = x0 + (x1-x0)*(y0-z0)/delta;
				if (xinter > x1)
				  dx = min(dx, 0.5*abs(xinter-x1)/(kmax_courbe-kmin_courbe));
			      }
			  }
		    
		    dx = max(dx, dx_min);
		    // DISP(dx);
		  }
	      }
	    else
	      {
		val_propre = coef_relation;
		err_propre = old_err_propre;
	      }
	    
	    // on met a jour Vh
	    for (int j = 0; j < (Vh.GetM()-1); j++)
	      {
		xsi_history(j) = xsi_history(j+1);
		for (int k = 0; k < nb_original_dof; k++)
		  Vh(j,k) = Vh(j+1,k);
	      }
	    
	    // le dernier instantane
	    xsi_history(Vh.GetM()-1) = xsi;
	    for (int k = 0; k < nb_original_dof; k++)
	      Vh(Vh.GetM()-1,k) = val_propre(k);
	    
	    for (int j = 0; j < nb_original_dof; j++)
	      multi_root<<val_propre(j)<<" ";
	    
	    for (int j = 0; j < nb_original_dof; j++)
	      file_err<<err_propre(j)<<" ";
	    
	    multi_root<<endl; file_err<<endl;
	    xsi += dx; nb_iter++;
	  }
	multi_root.close();
	file_err.close();
	file_out.close();
      }
    
    if (search_taylor_dvt)
      {
	// searching Taylor expansion
	int Np = order_dvt_taylor;
	Matrix<Real_wp> VanDerMonde(Np+1, Np+1);
	VectReal_wp coef_relation(Np+1);
	for (int i = 0; i <= Np; i++)
	  {
	    Real_wp xsi = kmax*Real_wp(i)/Real_wp(Np);
	    for (int j = 0; j <= Np; j++)
	      VanDerMonde(i,j) = pow(xsi, Real_wp(j));
	    
	    GetDispersionError(xsi, Dh, Kh, err1, omega_numer, eigen_values, eigen_vectors);
	    
	    coef_relation(i) = err1;
	  }
	
	IVect ipivot_vdm(VanDerMonde.GetM()); GetLU(VanDerMonde, ipivot_vdm);
	SolveLU(VanDerMonde, ipivot_vdm, coef_relation);
	DISP(coef_relation);
	cout<<"we compute rational approximation of each coefficient"<<endl;
	//PrintRational(coef_relation, threshold_rational);
	
      }
    
  }
  
  void ConstructMatrices(Matrix<Real_wp, General, ArrayRowSparse>& Dh,
			 Matrix<Real_wp, General, ArrayRowSparse>& Kh, VectReal_wp& xdof)
  {
    Globatto<Real_wp> lob, gauss;
    gauss.ConstructQuadrature(order, type_quadrature);
    lob.ConstructQuadrature(order, lob.QUADRATURE_LOBATTO);
    Matrix<int> Nodle;
    if (!dg_formulation)
      {
	Nodle.Reallocate(2,order+1);
	Nodle(0,0) = 0; Nodle(0,order) = 1;
	Nodle(1,0) = 1; Nodle(1,order) = 2;
	for (int i = 1; i < order; i++)
	  {
	    Nodle(0,i) = 2+i;
	    Nodle(1,i) = 1+order+i;
	  }
	
	xdof.Reallocate(order);
	for (int i = 0; i < order; i++)
	  xdof(i) = lob.Points(i);
	
      }
    else
      {
	xdof.Reallocate(order+1);
	for (int i = 0; i <= order; i++)
	  xdof(i) = lob.Points(i);
      }
    
    Matrix<Real_wp> Kh_loc(order+1, order+1), Mh_loc(order+1,order+1);
    Matrix<Real_wp> Rh_loc(order+1, order+1);
    Kh_loc.Fill(Real_wp(0));
    Mh_loc.Fill(Real_wp(0));
    Rh_loc.Fill(Real_wp(0));
    Real_wp phi_j, phi_k, dphi_j, dphi_k;
    for (int i = 0; i < gauss.GetNbPointsQuad(); i++)
      {
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    {
	      phi_j = lob.EvaluatePhi(j, gauss.Points(i));
	      phi_k = lob.EvaluatePhi(k, gauss.Points(i));
	      dphi_j = lob.EvaluatePhiGrad(j, gauss.Points(i));
	      dphi_k = lob.EvaluatePhiGrad(k, gauss.Points(i));
	      Mh_loc(j,k) += phi_j*phi_k*gauss.Weights(i);
	      Kh_loc(j,k) += dphi_j*dphi_k*gauss.Weights(i);
	      Rh_loc(j,k) += dphi_j*phi_k*gauss.Weights(i);
	    }
      }
    DISP(Mh_loc); DISP(Kh_loc); DISP(Rh_loc);
    if (ponderation_lumping)
      {
	cout << "not implemented" << endl;
	abort();
      }
    
    // assemblage
    if (dg_formulation)
      {
	int nodl = 3*(order+1);
	if (wave_equation)
	  nodl *= 2;
	
	Dh.Reallocate(nodl, nodl); Kh.Reallocate(nodl, nodl);
	// on ne calcule que les lignes associees au premier element
	if (wave_equation)
	  {
	    for (int i = 0; i <= order; i++)
	      for (int j = 0; j <= order; j++)
		{
		  Dh.AddInteraction(i, j, Mh_loc(i,j));
		  Dh.AddInteraction(order+1+i, order+1+j, Mh_loc(i,j));
		  Kh.AddInteraction(i, order+1+j, -Rh_loc(i,j));
		  Kh.AddInteraction(order+1+i, j, -Rh_loc(i,j));
		}
	    
	    // termes de flux
	    Kh.AddInteraction(0, order+1, -0.5);
	    Kh.AddInteraction(0, 3*(order+1)+order, -0.5);
	    Kh.AddInteraction(order+1, 0, -0.5);
	    Kh.AddInteraction(order+1, 2*(order+1)+order, -0.5);
	    
	    Kh.AddInteraction(order, (order+1)+order, 0.5);
	    Kh.AddInteraction(order, 5*(order+1), 0.5);
	    Kh.AddInteraction((order+1)+order, order, 0.5);
	    Kh.AddInteraction((order+1)+order, 4*(order+1), 0.5);
	    
	    // penalisation
	    Kh.AddInteraction(0, 0, alpha_dg);
	    Kh.AddInteraction(0, 2*(order+1)+order, -alpha_dg);
	    Kh.AddInteraction(order, order, -alpha_dg);
	    Kh.AddInteraction(order, 4*(order+1), alpha_dg);

	    Kh.AddInteraction((order+1), (order+1), delta_dg);
	    Kh.AddInteraction((order+1), 3*(order+1)+order, -delta_dg);
	    Kh.AddInteraction((order+1)+order, (order+1)+order, -delta_dg);
	    Kh.AddInteraction((order+1)+order, 5*(order+1), delta_dg);
	  }
	else
	  {
	    for (int i = 0; i <= order; i++)
	      for (int j = 0; j <= order; j++)
		{
		  Dh.AddInteraction(i, j, Mh_loc(i,j));
		  Kh.AddInteraction(i, j, Rh_loc(i,j));
		}
	    
	    // termes de flux
	    Kh.AddInteraction(0, 0, 0.5);
	    Kh.AddInteraction(0, 2*order+1, 0.5);
	    Kh.AddInteraction(order, order, -0.5);
	    Kh.AddInteraction(order, 2*order+2, -0.5);
	  }
      }
    else
      {
	int nodl = 2*order+1;
	Dh.Reallocate(nodl, nodl); Kh.Reallocate(nodl, nodl);
	for (int num_elem = 0; num_elem < 2; num_elem++)
	  {
	    for (int iloc = 0; iloc <= order; iloc++)
	      for (int jloc = 0; jloc <= order; jloc++)
		{
		  int iglob = Nodle(num_elem,iloc);
		  int jglob = Nodle(num_elem,jloc);
		  Dh.AddInteraction(iglob, jglob, Mh_loc(iloc,jloc));
		  Kh.AddInteraction(iglob, jglob, Kh_loc(iloc,jloc));
		}
	  }
      }
    DISP(Dh); DISP(Kh);
  }
  
};

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      cout<<"entrez le nom du fichier de donnees"<<endl;
      return -1;
    }
  
  InitMontjoie(argc, argv);
  
  DispersionSolver disp;
  ReadInputFile(string(argv[1]), disp);
  
  cout.precision(16);
  
  Real_wp epsilon_machine = GetPrecisionMachine(Real_wp(1));
  DISP(epsilon_machine);
  
  Matrix<Real_wp, General, ArrayRowSparse>  Dh, Kh; VectReal_wp xdof;
  disp.ConstructMatrices(Dh, Kh, xdof);  
  disp.GetPeriodicityDofs();
  
  cout<<"Periodicity treated"<<endl;
  disp.GetDispersionRelation(Dh, Kh, xdof);
  
  return 0;
}
