#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell3D.hxx"

using namespace Montjoie;

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x)
{
  int type = rand()%3;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}

template<class T>
class EigenvalueComparisonMode : public EigenvalueComparisonClass<T>
{
  VectComplex_wp poles;
  Real_wp r0; Complex_wp shift; Real_wp omega;

public:
  void SetPoles(const VectComplex_wp& P, const Real_wp& r) { poles = P; r0 = r; }
  void SetShift(const Real_wp& Om, const Complex_wp& s) { omega = Om; shift = s; }
  
  int CompareEigenvalue(const T& L, const T&, const T& L2, const T&)
  {
    T om; to_complex(omega*(shift + Real_wp(1)/L), om);
    T om2; to_complex(omega*(shift + Real_wp(1)/L2), om2);
    Real_wp dist = abs(L), dist2 = abs(L2);
    for (int i = 0; i < poles.GetM(); i++)
      {
        if (abs(om - poles(i)) <= r0)
          dist *= absSquare(om-poles(i))/r0;

        if (abs(om2 - poles(i)) <= r0)
          dist2 *= absSquare(om2-poles(i))/r0;
      }

    if (dist > dist2)
      return -1;

    return 1;
  }

};

class InputVariables : public InputDataProblem_Base
{
public :
  // parametres pour calculer les valeurs propres (shifts et nombre de vps pour chaque shift)
  Vector<int> number_eigenval;
  Vector<Complex_wp> shift_cplx;
  VectReal_wp radius_shift;
  VectReal_wp ratio_shift, angle_shift;
  Real_wp tolerance_eigenvalue;
  int nb_max_iterations;

  // spectre a selectionner
  Real_wp Lr_max, Li_max, Li_min, spectrum_radius_min;
  Vector<int> number_eigenval_select;

  // solveur a choisir
  int type_eigensolver;
  bool use_rational_eigensolver;
  
  // doit on rajouter des PML au maillage ?
  bool add_pml;

  // fichier contenant les valeurs propres
  string file_name_eigen;

  // type de relevement (pour calculer la solution a l'exterieur de la cavite)
  int type_relevement;
  enum { STANDARD};

  // methode de projection pour calculer les alpha_m
  int type_projection; Complex_wp pole_muljarov;
  enum {PROJ_INTEGRAL, PROJ_BIORTHO, PROJ_MARSEILLE, PROJ_MULJAROV};
  int use_p_unknown;

  // gamme de frequences pour le probleme source
  Real_wp fmin, fmax; int nb_freq, first_freq;

  // pour les PMLs non-dispersives
  bool dispersive_pml;
  bool manual_choice_pml_modes;
  bool compute_fem_sol;
  bool load_fem_sol, write_fem_sol;
  string path_fem_sol;
  
  // fichiers ou on ecrit la solution QNM (et FEM)
  string file_qnm_sol, file_fem_sol;

  // stocke-t-on les matrices Mh et Kh ?
  bool explicit_matrix;

  // valeurs propres de reference (pour dissocier les PMLs des QNM)
  Real_wp threshold_ratio_qnm;
  bool extractQNM;

  // est-ce que les shifts sont adaptes a/remplissent la Spectrum Box qui est definie
  int lambda_par_shift;
  Complex_wp shift_init;
  bool fill_box;
  Real_wp shift_step;
  int nb_extrema_values_to_remove;

  // Threshold_EigenValue
  Real_wp threshold_eigen, threshold_eigen_vector, threshold_eigenval_equal;

  bool source_wei;
  bool compute_err_grad;
  int reference_error;  
  
  bool check_orthogonality;
  bool save_eigen_vec, load_eigen_vec;
  Vector<string> file_name_eigenvectors;
  bool display_eigenvectors;
  string name_output_eigenmode_file;

  bool merge_eigenvalues;
  VectReal_wp exclude_eigenval;
  int rank_proc;
#ifdef SELDON_WITH_MPI
  MPI_Comm comm;
#endif

  Vector<int> ref_cavity;

  bool estimate_number_eigenval;
  int type_integration_feast, nb_points_quad_feast;
  bool drop_negative_modes;

  bool avoid_poles, avoid_zeros; Real_wp radius_poles;
  
  // ordre d'interpolation et nombre d'intervalles
  int order_interpolation, nb_interval_interpolation;
  int nb_levels_interpolation, type_interpolation;
  enum {TCHEBY, CLENSHAW, LEJA};
  
  
  InputVariables()
  {    
    number_eigenval.Reallocate(1);
    number_eigenval(0) = 100;
    tolerance_eigenvalue = 1e-12;
    nb_max_iterations = 1000;

    Lr_max = Real_wp(0.6); Li_max = Real_wp(30); Li_min = Real_wp(0);
    spectrum_radius_min = Real_wp(0);

    merge_eigenvalues = false;
    type_eigensolver = TypeEigenvalueSolver::ARPACK;
    use_rational_eigensolver = false;
    
    add_pml = false;
    file_name_eigen = string("Lambda");
    
    type_relevement = STANDARD;
    type_projection = PROJ_BIORTHO;
    
    fmin = 0.01; fmax = 2.0; nb_freq = 0; first_freq = 0;

    dispersive_pml = true;
    manual_choice_pml_modes = true;

    compute_fem_sol = true;
    file_qnm_sol = "SolQnm";
    file_fem_sol = "SolFem";
    load_fem_sol = false;
    write_fem_sol = false;
    path_fem_sol = ".";

    explicit_matrix = false;

    threshold_ratio_qnm = 10.0;    
    extractQNM = false;

    threshold_eigen=1e-6;
    threshold_eigenval_equal = 1e-6;
    fill_box = false;
    nb_extrema_values_to_remove = 1;

    R2::threshold = 1e-6;
    threshold_eigen_vector = 1e-10;
    
    compute_err_grad = false;
    reference_error = -1;
    
    source_wei = false;
    check_orthogonality = false;
    save_eigen_vec = false;
    load_eigen_vec = false;
    display_eigenvectors = false;
    name_output_eigenmode_file = "ModeU";
    
    estimate_number_eigenval = false;
    rank_proc = 0;
    drop_negative_modes = false;

    avoid_poles = false;
    avoid_zeros = false;
    radius_poles = 0.1;

    ref_cavity.Reallocate(1); ref_cavity(0) = 2;

    type_integration_feast = -1;
    nb_points_quad_feast = 0;

    order_interpolation = 0;
    nb_interval_interpolation = 0;
    nb_levels_interpolation = 1;
    type_interpolation = TCHEBY;
  }
  
  bool SelectEigenvalue(const Complex_wp& z)
  {
    if ((realpart(z) <= Lr_max)
	&& (abs(imagpart(z)) <= Li_max)
	&& (abs(imagpart(z)) >= Li_min)
	&& (abs(z) >= spectrum_radius_min))
      {
        for (int k = 0; k < exclude_eigenval.GetM(); k++)
          if (abs(exclude_eigenval(k) - z) < threshold_eigen)
            return false;
        
        return true;
      }
    
    return false;
  }

  bool SelectEigenvector(const VectComplex_wp& x, int N)
  {
    Real_wp normX = 0, normY;
    for (int i = 0; i < N; i++)
      normX += absSquare(x(i));

    normY = square(Norm2(x));
    R2 sum_tmp(normX, normY), sum_res;
    MPI_Allreduce(&sum_tmp(0), &sum_res(0), 2, MPI_DOUBLE, MPI_SUM, comm);
    normX = sum_res(0); normY = sum_res(1);
    normX = sqrt(normX); normY = sqrt(normY);
    if (normX / normY > threshold_eigen_vector)
      return true;

    return false;
  }
  
  void SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (keyword == "SpectrumBox")
      {
	Lr_max = to_num<Real_wp>(param(0));
	Li_max = to_num<Real_wp>(param(1));
        if (param.GetM() > 2)
          Li_min = to_num<Real_wp>(param(2));
      }
    else if (keyword == "SpectrumNumber")
      number_eigenval_select.ReadText(param(0));
    else if (keyword == "SpectrumRadius")
      spectrum_radius_min = to_num<Real_wp>(param(0));
    else if (keyword == "EigenvalueTolerance")
      tolerance_eigenvalue = to_num<Real_wp>(param(0));
    else if (keyword == "EigenvalueMaxNumberIterations")
      nb_max_iterations = to_num<int>(param(0));
    else if (keyword == "NbEigenvaluesToCompute")
      {
	number_eigenval.Reallocate(param.GetM());
	for (int i = 0; i < param.GetM(); i++)
	  number_eigenval(i) = to_num<int>(param(i));
      }
    else if (keyword == "ShiftEigenvalue")
      {
	shift_cplx.Reallocate(param.GetM());
	for (int i = 0; i < param.GetM(); i++)
	  {
            if (param(i) == "Large")
              shift_cplx(i) = 1e31;
            else
              shift_cplx(i) = to_num<Complex_wp>(param(i));
          }
      }
    else if (keyword == "AvoidPoles")
      {
        if (param(0) == "YES")
          avoid_poles = true;
        else if (param(0) == "ALL")
          {
            avoid_poles = true;
            avoid_zeros = true;
          }
        else
          avoid_poles = false;
        
        radius_poles = to_num<Real_wp>(param(1));
      }
    else if (keyword == "RadiusShift")
      {
        radius_shift.Reallocate(param.GetM());
        for (int i = 0; i < param.GetM(); i++)
          radius_shift(i) = to_num<Real_wp>(param(i));
      }
    else if (keyword == "RatioShift")
      {
        ratio_shift.Reallocate(param.GetM());
        for (int i = 0; i < param.GetM(); i++)
          ratio_shift(i) = to_num<Real_wp>(param(i));
      }
    else if (keyword == "AngleShift")
      {
        angle_shift.Reallocate(param.GetM());
        for (int i = 0; i < param.GetM(); i++)
          angle_shift(i) = to_num<Real_wp>(param(i));
      }
    else if (keyword == "IntegrationFeast")
      {
        type_integration_feast = to_num<int>(param(0));
        nb_points_quad_feast = to_num<int>(param(1));
      }
    else if (keyword == "EstimateNumberEigenvalues")
      {
        if (param(0) == "YES")
          estimate_number_eigenval = true;
        else
          estimate_number_eigenval = false;
      }
    else if (keyword == "EigenvalueSolver")
      {
        if (param(0) == "Arpack")
          type_eigensolver = TypeEigenvalueSolver::ARPACK;
        else if (param(0) == "Slepc")
          type_eigensolver = TypeEigenvalueSolver::SLEPC;
        else if (param(0) == "Feast")
          type_eigensolver = TypeEigenvalueSolver::FEAST;
        else if (param(0) == "Lapack")
          type_eigensolver = -1;
        else if (param(0) == "Rational")
          use_rational_eigensolver = true;
      }
    else if (keyword == "DropNegativeModes")
      {
        if (param(0) == "YES")
          drop_negative_modes = true;
        else
          drop_negative_modes = false;
      }
    else if (keyword == "AddIntermediateLayerPML")
      {
	if (param(0) == "YES")
	  add_pml = true;
      }
    else if (keyword == "ReferenceCavity")
      {
        ref_cavity.Reallocate(param.GetM());
        for (int i = 0; i < param.GetM(); i++)
          ref_cavity(i) = to_num<int>(param(i));
      }
    else if (keyword == "FileEigenvalue")
      file_name_eigen = param(0);
    else if (keyword == "ChoiceModesPML")
      {
	if (param(0) == "AUTO")
	  manual_choice_pml_modes = false;
      }
    else if (keyword == "RangeWaveLength")
      {
	Real_wp Lmin = to_num<Real_wp>(param(0));
	Real_wp Lmax = to_num<Real_wp>(param(1));
	nb_freq = to_num<int>(param(2));

        fmin = Real_wp(1) / Lmax;
        fmax = Real_wp(1) / Lmin;
        if (param.GetM() >= 4)
          first_freq = to_num<int>(param(3));
      }
    else if (keyword == "RangeFrequency")
      {
	fmin = to_num<Real_wp>(param(0));
	fmax = to_num<Real_wp>(param(1));
	nb_freq = to_num<int>(param(2));

        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
	  {
            fmin *= Real_wp(1) / PhysicalConstant::speed_light;
            fmax *= Real_wp(1) / PhysicalConstant::speed_light;
          }

        if (param.GetM() >= 4)
          first_freq = to_num<int>(param(3));
      }
    else if (keyword == "ProjectionQNM")
      {
	if (param(0) == "Integral")
	  type_projection = PROJ_INTEGRAL;
	else if (param(0) == "Biortho")
	  type_projection = PROJ_BIORTHO;
        else if (param(0) == "Marseille")
          type_projection = PROJ_MARSEILLE;
        else if (param(0) == "Muljarov")
          type_projection = PROJ_MULJAROV;
        else
	  {
	    cout << "Unknown projection " << param(0) << endl;
	    abort();
	  }

        use_p_unknown = -1;
        if (param.GetM() > 1)
          {
            if (param(1) == "P")
              use_p_unknown = 0;
            else if (param(1) == "Q")
              use_p_unknown = 1;
            else
              use_p_unknown = to_num<int>(param(1));
          }
      }
    else if (keyword == "InterpolationQNM")
      {
        nb_interval_interpolation = to_num<int>(param(0));
        order_interpolation = to_num<int>(param(1));
        if (param.GetM() > 2)
          {
            nb_levels_interpolation = to_num<int>(param(3));
            if (param(2) == "LEJA")
              type_interpolation = LEJA;
            else if (param(2) == "CLENSHAW")
              type_interpolation = CLENSHAW;
            else
              type_interpolation = TCHEBY;
          }
      }
    else if (keyword == "DispersivePML")
      {
	if (param(0) == "YES")
	  dispersive_pml = true;
	else if (param(0) == "NO")
          dispersive_pml = false;
	else
	  {
	    cout << "Invalid PML" << endl;
	    abort();
	  }
      }
    else if (keyword == "FileSolutionQNM")
      file_qnm_sol = param(0);
    else if (keyword == "FileSolutionFEM")
      {
        compute_fem_sol = true;
        load_fem_sol = false;
        write_fem_sol = false;
        if (param(0) == "NO")
          compute_fem_sol = false;
        else if (param(0) == "Load")
          {
            load_fem_sol = true;
            path_fem_sol = param(1);
            file_fem_sol = param(2);
          }
        else if (param(0) == "Save")
          {
            write_fem_sol = true;
            path_fem_sol = param(1);
            file_fem_sol = param(2);
          }
        else
          file_fem_sol = param(0);
      }
    else if (keyword == "ExplicitMatrixEigen")
      {
	if (param(0) == "YES")
	  explicit_matrix = true;
	else
	  explicit_matrix = false;
      }
    else if (keyword == "ExtractQNM")
      {
	if (param(0) == "YES")
	  {
	    extractQNM = true;
	    if (param.GetM() > 1)
	      threshold_ratio_qnm = to_num<Real_wp>(param(1));
	  }
	else
	  extractQNM = false;
      }
    else if (keyword == "FillBox")
      {
	if (param(0) == "YES")
	  {
	    fill_box = true;
	    if (param.GetM() > 1)
	      {
		lambda_par_shift = to_num<int>(param(1));
		shift_init = to_num<Complex_wp>(param(2));
		shift_step = to_num<Real_wp>(param(3));
                if (param.GetM() > 4)
                  nb_extrema_values_to_remove = to_num<int>(param(4));
	      }	    
	  }
	else
	  fill_box = false;
      }
    else if (keyword == "ThresholdEigenValue")
      {
	threshold_eigen = to_num<Real_wp>(param(0));
      }
    else if (keyword == "ThresholdEigenvalueEqual")
      threshold_eigenval_equal = to_num<Real_wp>(param(0));
    else if (keyword == "UseWeiSource")
      {
        if (param(0) == "YES")
          source_wei = true;
        else
          source_wei = false;
      }
    else if (keyword == "ComputeError")
      {
        if (param(0) == "GRAD")
          compute_err_grad = true;
        else
          compute_err_grad = false;

        if (param.GetM() > 1)
          reference_error = to_num<int>(param(1));
      }
    else if (keyword == "LoadEigenvectors")
      {
        if (param(0) == "YES")
          {
            load_eigen_vec = true;
            file_name_eigenvectors.Reallocate(1);
            file_name_eigenvectors(0) = param(1);
          }
      }
    else if (keyword == "SaveEigenvectors")
      {
        if (param(0) == "YES")
          {
            save_eigen_vec = true;
            file_name_eigenvectors.Reallocate(1);
            file_name_eigenvectors(0) = param(1);
          }
      }
    else if (keyword == "MergeEigenvectors")
      {
        merge_eigenvalues = true;
        file_name_eigenvectors = param;
      }
    else if (keyword == "CheckOrthogonality")
      {
        if (param(0) == "YES")
          check_orthogonality = true;
        else
          check_orthogonality = false;
      }
    else if (keyword == "DisplayEigenvectors")
      {
        if (param(0) == "YES")
          {
            display_eigenvectors = true;
            if (param.GetM() > 1)
              name_output_eigenmode_file = param(1);
          }
        else
          display_eigenvectors = false;
      }
  }
  
  void UpdateAdim(Real_wp L)
  {
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      {
        fmin *= L; fmax *= L;
      }
  }
  
};

void ExtractEigenvector(const Matrix<Real_wp>& new_eigenvec, int n, int type,
			VectComplex_wp& x)
{
  switch (type)
    {
    case 0 :
      {
	for (int p = 0; p < x.GetM(); p++)
	  x(p) = new_eigenvec(p, n);
      }
      break;
    case 1 :
      {
	for (int p = 0; p < x.GetM(); p++)
	  x(p) = Complex_wp(new_eigenvec(p, n), new_eigenvec(p, n+1));	
      }
      break;
    case 2 :
      {
	for (int p = 0; p < x.GetM(); p++)
	  x(p) = Complex_wp(new_eigenvec(p, n-1), -new_eigenvec(p, n));	
      }
      break;
    }  
}

void FillEigenvector(const VectComplex_wp& x, int n, int type,
		     Matrix<Real_wp>& new_eigenvec)
{
  switch (type)
    {
    case 0 :
      {
	for (int p = 0; p < x.GetM(); p++)
	  new_eigenvec(p, n) = realpart(x(p));
      }
      break;
    case 1 :
      {
	for (int p = 0; p < x.GetM(); p++)
	  {
	    new_eigenvec(p, n) = realpart(x(p));
	    new_eigenvec(p, n+1) = imagpart(x(p));
	  }
      }
      break;
    case 2 :
      {
	for (int p = 0; p < x.GetM(); p++)
	  {
	    new_eigenvec(p, n-1) = realpart(x(p));
	    new_eigenvec(p, n) = -imagpart(x(p));
	  }
      }
      break;
    }  
}

void FindEigenvalueNumber(const VectComplex_wp& Lambda, const VectComplex_wp& lambda,
			  const Vector<int>& imag_eigenvec, DistributedVector<Complex_wp>& x,
			  bool duplicate_allowed, IVect& num_mode, Matrix<Real_wp>& new_eigenvec)
{
  VectR2 old_eigenval(Lambda.GetM());
  for (int i = 0; i < Lambda.GetM(); i++)
    old_eigenval(i).Init(realpart(Lambda(i)), imagpart(Lambda(i)));
  
  VectR2 new_eigenval(lambda.GetM());
  for (int i = 0; i < lambda.GetM(); i++)
    new_eigenval(i).Init(realpart(lambda(i)), imagpart(lambda(i)));
  
  IVect permut_old(Lambda.GetM()); permut_old.Fill();
  Sort(old_eigenval, permut_old);
  
  IVect permut_new(lambda.GetM()); permut_new.Fill();
  Sort(new_eigenval, permut_new);
  
  num_mode.Reallocate(lambda.GetM());
  num_mode.Fill(-1);
  int k = 0;
  for (int i = 0; i < old_eigenval.GetM(); i++)
    {
      if ((i > 0) && (old_eigenval(i) == old_eigenval(i-1)))
	continue;
      
      while ((k < lambda.GetM()) && (new_eigenval(k) < old_eigenval(i)))
	k++;
      
      if ((k < lambda.GetM()) && (new_eigenval(k) == old_eigenval(i)))
	{
	  // we found the corresponding eigenvalue
	  num_mode(permut_new(k)) = permut_old(i);
	}

      if (duplicate_allowed)
	{
	  while ((k < lambda.GetM()) && (new_eigenval(k) == old_eigenval(i)))
	    {
	      num_mode(permut_new(k)) = permut_old(i);
	      k++;
	    }
	}
      else
	{
	  if ((k < lambda.GetM()) && (new_eigenval(k) == old_eigenval(i)))
	    {
	      // we found the corresponding eigenvalue
	      num_mode(permut_new(k)) = permut_old(i);
	    }
	  
	  if ((k+1 < lambda.GetM()) && (new_eigenval(k+1) == old_eigenval(i)))
	    {
	      // double eigenvalue not handled
	      cout << "Multiple eigenvalue not treated" << endl;
	      abort();
	    }
	}
    }

  // among new eigenvalues, we check that there is no "duplicate"
  k = 0;
  DistributedVector<Complex_wp> y(x);
  while (k < new_eigenval.GetM())
    {
      int n = 0;
      while ((k+n < lambda.GetM()) && (new_eigenval(k+n) == new_eigenval(k)))
	n++;
      
      if (new_eigenval(k)(1) > 0)
	{
	  k += n;
	  continue;
	}
      
      if (n > 1)
	{
	  // on cree une matrice de masse associee a la base des vecteurs propres
	  Matrix<Complex_wp, Hermitian, RowHermPacked> mass(n, n); VectReal_wp lambda_mass;
	  Matrix<Complex_wp> eigenvec_mass;
	  for (int i = 0; i < n; i++)
	    for (int j = i; j < n; j++)
	      {
		ExtractEigenvector(new_eigenvec, permut_new(k+i),
				   imag_eigenvec(permut_new(k+i)), x);
		
		ExtractEigenvector(new_eigenvec, permut_new(k+j),
				   imag_eigenvec(permut_new(k+j)), y);
		
		mass.Set(i, j, DotProdConj(x, y));
	      }
          
	  GetEigenvaluesEigenvectors(mass, lambda_mass, eigenvec_mass);

	  int dim_kernel = 0;
	  for (int j = 0; j < lambda_mass.GetM(); j++)
	    if (abs(lambda_mass(j)) <= 1e-12)
	      dim_kernel++;
	  
	  if (dim_kernel > 0)
	    {
	      Vector<VectComplex_wp> vec_E(n);
	      for (int i = 0; i < n; i++)
		{
		  vec_E(i).Reallocate(x.GetM());
		  ExtractEigenvector(new_eigenvec, permut_new(k+i),
				     imag_eigenvec(permut_new(k+i)), vec_E(i));
		}
	      
	      // on garde les vecteurs propres qui sont independants
	      int num = 0;
	      for (int j = 0; j < lambda_mass.GetM(); j++)
		if (abs(lambda_mass(j)) > 1e-12)
		  {
		    x.Zero();
		    for (int i = 0; i < n; i++)
		      Add(eigenvec_mass(i, j), vec_E(i), x);

		    FillEigenvector(x, permut_new(k+num),
				    imag_eigenvec(permut_new(k+num)), new_eigenvec);
		    
		    num++;
		  }

	      for (int j = num; j < n; j++)
		{
		  // on met num_mode a zero pour les vecteurs propres restants		  
		  num_mode(permut_new(k+j)) = 0;
		  if (imag_eigenvec(permut_new(k+j)) == 1)
		    num_mode(permut_new(k+j) + 1) = 0;
		  else if (imag_eigenvec(permut_new(k+j)) == 2)
		    num_mode(permut_new(k+j) - 1) = 0;
		}
	    }
	}
      
      k += n;
    }
}


void FindEigenvalueNumber(const VectComplex_wp& Lambda, const VectComplex_wp& lambda,
			  DistributedVector<Complex_wp>& x,
			  bool duplicate_allowed, IVect& num_mode, Matrix<Complex_wp>& new_eigenvec)
{
  VectR2 old_eigenval(Lambda.GetM());
  for (int i = 0; i < Lambda.GetM(); i++)
    old_eigenval(i).Init(realpart(Lambda(i)), imagpart(Lambda(i)));
  
  VectR2 new_eigenval(lambda.GetM());
  for (int i = 0; i < lambda.GetM(); i++)
    new_eigenval(i).Init(realpart(lambda(i)), imagpart(lambda(i)));
  
  IVect permut_old(Lambda.GetM()); permut_old.Fill();
  Sort(old_eigenval, permut_old);
  
  IVect permut_new(lambda.GetM()); permut_new.Fill();
  Sort(new_eigenval, permut_new);
  
  num_mode.Reallocate(lambda.GetM());
  num_mode.Fill(-1);
  int k = 0;
  for (int i = 0; i < old_eigenval.GetM(); i++)
    {
      if ((i > 0) && (old_eigenval(i) == old_eigenval(i-1)))
	continue;
      
      while ((k < lambda.GetM()) && (new_eigenval(k) < old_eigenval(i)))
	k++;
      
      if ((k < lambda.GetM()) && (new_eigenval(k) == old_eigenval(i)))
	{
	  // we found the corresponding eigenvalue
	  num_mode(permut_new(k)) = permut_old(i);
	}

      if (duplicate_allowed)
	{
	  while ((k < lambda.GetM()) && (new_eigenval(k) == old_eigenval(i)))
	    {
	      num_mode(permut_new(k)) = permut_old(i);
	      k++;
	    }
	}
      else
	{
	  if ((k < lambda.GetM()) && (new_eigenval(k) == old_eigenval(i)))
	    {
	      // we found the corresponding eigenvalue
	      num_mode(permut_new(k)) = permut_old(i);
	    }
	  
	  if ((k+1 < lambda.GetM()) && (new_eigenval(k+1) == old_eigenval(i)))
	    {
	      // double eigenvalue not handled
	      cout << "Multiple eigenvalue not treated" << endl;
	      abort();
	    }
	}
    }

  // among new eigenvalues, we check that there is no "duplicate"
  k = 0;
  DistributedVector<Complex_wp> y(x);
  while (k < new_eigenval.GetM())
    {
      int n = 0;
      while ((k+n < lambda.GetM()) && (new_eigenval(k+n) == new_eigenval(k)))
	n++;
      
      if (n > 1)
	{
	  // on cree une matrice de masse associee a la base des vecteurs propres
	  Matrix<Complex_wp, Hermitian, RowHermPacked> mass(n, n); VectReal_wp lambda_mass;
	  Matrix<Complex_wp> eigenvec_mass;
	  for (int i = 0; i < n; i++)
	    for (int j = i; j < n; j++)
	      {
                GetCol(new_eigenvec, permut_new(k+i), x);
                GetCol(new_eigenvec, permut_new(k+j), y);
		mass.Set(i, j, DotProdConj(x, y));
	      }
          
	  GetEigenvaluesEigenvectors(mass, lambda_mass, eigenvec_mass);

	  int dim_kernel = 0;
	  for (int j = 0; j < lambda_mass.GetM(); j++)
	    if (abs(lambda_mass(j)) <= 1e-12)
	      dim_kernel++;
	  
	  if (dim_kernel > 0)
	    {
	      Vector<VectComplex_wp> vec_E(n);
	      for (int i = 0; i < n; i++)
		{
		  vec_E(i).Reallocate(x.GetM());
                  GetCol(new_eigenvec, permut_new(k+i), vec_E(i));
		}
              
	      // on garde les vecteurs propres qui sont independants
	      int num = 0;
	      for (int j = 0; j < lambda_mass.GetM(); j++)
		if (abs(lambda_mass(j)) > 1e-12)
		  {
		    x.Zero();
		    for (int i = 0; i < n; i++)
		      Add(eigenvec_mass(i, j), vec_E(i), x);

                    SetCol(x, permut_new(k+num), new_eigenvec);
		    num++;
		  }

	      for (int j = num; j < n; j++)
		{
		  // on met num_mode a zero pour les vecteurs propres restants
		  num_mode(permut_new(k+j)) = 0;
		}
	    }
	}
      
      k += n;
    }
}

void GetComplexEigenvalues(const VectReal_wp& Lr, const VectReal_wp& Li, VectComplex_wp& L)
{
  int nev = Lr.GetM();
  L.Reallocate(nev);
  for (int i = 0; i < nev; i++)
    L(i) = Complex_wp(Lr(i), Li(i));
}

void MergeEigenvalues(const VectComplex_wp& new_lambda, const VectComplex_wp& new_lambda_imag,
		      Matrix<Complex_wp>& new_eigenvec, VectComplex_wp& lambda,
		      VectComplex_wp& lambda_imag, Matrix<Complex_wp>& eigen_vec,
                      DistributedVector<Complex_wp>& x, bool drop_negative_modes,
                      bool drop_real_negative = false)
{
  IVect num_mode;
  FindEigenvalueNumber(lambda, new_lambda, x, true, num_mode, new_eigenvec);

  if (drop_negative_modes)
    for (int i = 0; i < new_lambda.GetM(); i++)
      if (imagpart(new_lambda(i)) < 0)
        num_mode(i) = 0;

  if (drop_real_negative)
    for (int i = 0; i < new_lambda.GetM(); i++)
      if (realpart(new_lambda(i)) < 0)
        num_mode(i) = 0;
  
  int nb_new = 0;
  for (int k = 0; k < num_mode.GetM(); k++)
    if (num_mode(k) == -1)
      nb_new++;
  
  int nev = lambda.GetM();
  int N = new_eigenvec.GetM();
  lambda.Resize(nev + nb_new);
  lambda_imag.Resize(nev + nb_new);
  lambda_imag.Zero();
  eigen_vec.Resize(N, nev + nb_new);

  nb_new = nev;
  for (int k = 0; k < new_lambda.GetM(); k++)
    if (num_mode(k) == -1)
      {
        lambda(nb_new) = new_lambda(k);
        for (int j = 0; j < N; j++)
          eigen_vec(j, nb_new) = new_eigenvec(j, k);
        
        nb_new++;
      }
}


void MergeEigenvalues(const VectReal_wp& new_lambda, const VectReal_wp& new_lambda_imag,
		      Matrix<Real_wp>& new_eigenvec, VectReal_wp& lambda,
		      VectReal_wp& lambda_imag, Matrix<Real_wp>& eigen_vec,
                      DistributedVector<Complex_wp>& x, bool drop_negative_modes,
                      bool drop_real_negative = false)
{
  VectComplex_wp Lambda(lambda.GetM()), NewLambda(new_lambda.GetM());
  GetComplexEigenvalues(lambda, lambda_imag, Lambda);
  GetComplexEigenvalues(new_lambda, new_lambda_imag, NewLambda);

  Vector<int> imag_eigenvec(new_lambda.GetM()); imag_eigenvec.Fill(0);
  int k = 0;
  while (k < new_lambda.GetM())
    {
      if (new_lambda_imag(k) == Real_wp(0))
	k++;
      else
	{
	  imag_eigenvec(k) = 1;
	  imag_eigenvec(k+1) = 2;
	  k += 2;
	}
    }
    
  IVect num_mode;
  FindEigenvalueNumber(Lambda, NewLambda, imag_eigenvec, x, true, num_mode, new_eigenvec);
  
  int nb_new = 0;
  k = 0;
  while (k < new_lambda.GetM())
    {
      if (new_lambda_imag(k) == Real_wp(0))
	{
	  if (num_mode(k) == -1)
	    nb_new++;
	  
	  k++;
	}
      else
	{
	  if ((num_mode(k) == -1) && (num_mode(k+1) == -1))
	    nb_new += 2;
	  
	  k += 2;
	}
    }

  int nev = lambda.GetM();
  int N = new_eigenvec.GetM();
  lambda.Resize(nev + nb_new);
  lambda_imag.Resize(nev + nb_new);
  eigen_vec.Resize(N, nev + nb_new);

  nb_new = nev;
  k = 0;
  while (k < new_lambda.GetM())
    {
      if (new_lambda_imag(k) == Real_wp(0))
	{
	  if (num_mode(k) == -1)
	    {
	      lambda(nb_new) = new_lambda(k);
	      lambda_imag(nb_new) = new_lambda_imag(k);
	      for (int j = 0; j < N; j++)
		eigen_vec(j, nb_new) = new_eigenvec(j, k);
	      
	      nb_new++;
	    }
	  
	  k++;
	}
      else
	{
	  if ((num_mode(k) == -1) && (num_mode(k+1) == -1))
	    {
	      if (new_lambda_imag(k) > 0)
		{
		  lambda(nb_new) = new_lambda(k);
		  lambda_imag(nb_new) = -new_lambda_imag(k);
		  lambda(nb_new+1) = new_lambda(k);
		  lambda_imag(nb_new+1) = new_lambda_imag(k);
		  
		  for (int j = 0; j < N; j++)
		    {
		      eigen_vec(j, nb_new) = new_eigenvec(j, k);
		      eigen_vec(j, nb_new+1) = -new_eigenvec(j, k+1);
		    }
		}
	      else
		{
		  lambda(nb_new) = new_lambda(k);
		  lambda_imag(nb_new) = new_lambda_imag(k);
		  lambda(nb_new+1) = new_lambda(k);
		  lambda_imag(nb_new+1) = -new_lambda_imag(k);
		  
		  for (int j = 0; j < N; j++)
		    {
		      eigen_vec(j, nb_new) = new_eigenvec(j, k);
		      eigen_vec(j, nb_new+1) = new_eigenvec(j, k+1);
		    }
		}
	      
	      nb_new += 2;
	    }
	  
	  k += 2;
	}
    }
}


void UpdateEigenValues(const FemMatrixFreeClass_Base<Real_wp>& mat_mass,
                       const FemMatrixFreeClass_Base<Real_wp>& mat_stiff,
		       VectReal_wp& lambda0_real, VectReal_wp& lambda0_imag, Matrix<Real_wp>& eigen_vec0,
		       Real_wp threshold)
{
  int N = mat_stiff.GetM();
  Vector<Real_wp> lambda(lambda0_real),lambda_imag(lambda0_imag);
  Matrix<Real_wp> eigen_vec(eigen_vec0);

  const DistributedProblem_Base& var_comm = mat_stiff.GetDistributedProblem();
  DistributedVector<Real_wp> X(var_comm.GetOverlappedDofNumber(), var_comm.comm_group_mode);
  X.Reallocate(N);

  DistributedVector<Complex_wp> Diff(var_comm.GetOverlappedDofNumber(), var_comm.comm_group_mode);
  Diff.Reallocate(N); Diff.Zero();
  
  int rank_proc = var_comm.GetRankProcMode();
  DistributedVector<Real_wp> Xi(X), Y(X), Yi(X), Mx(X), Mxi(X);
  X.Fill(0); Xi.Fill(0);
  Y.Fill(0); Yi.Fill(0);
  Mx.Zero(); Mxi.Zero();
  
  int i = 0;
  int num = 0; Real_wp err_inf = 0;
  while (i < lambda.GetM())
    {
      bool eigen_pair = false;
      if (i < lambda.GetM()-1)
        {
          if ( (lambda(i) == lambda(i+1)) && (lambda_imag(i) == -lambda_imag(i+1)))
            eigen_pair = true;
        }
      
      if ((lambda_imag(i) != 0) && (!eigen_pair))
	break;
	
      Real_wp err = 0, normeX(1);
      if (eigen_pair)
        {
          for (int j = 0; j < N; j++)
            {
              X(j) = eigen_vec(j, i);
              Xi(j) = eigen_vec(j, i+1);
            }
          
          mat_stiff.MltVector(X, Y); 
          mat_stiff.MltVector(Xi, Yi);
          mat_mass.MltVector(X, Mx);
          mat_mass.MltVector(Xi, Mxi);
          normeX = sqrt(DotProd(Mx, Mx) + DotProd(Mxi, Mxi));
          for (int j = 0; j < N; j++)
            Diff(j) = Complex_wp(Y(j), Yi(j))
              - Complex_wp(lambda(i), lambda_imag(i))*Complex_wp(Mx(j), Mxi(j));
                    
          err = Norm2(Diff);
        }
      else
        {
          for (int j = 0; j < N; j++)
            X(j) = eigen_vec(j, i);
          
          mat_stiff.MltVector(X, Y);
	  mat_mass.MltVector(X, Mx);
          normeX = sqrt(DotProd(Mx, Mx));
          for (int j = 0; j < N; j++)
            Xi(j) = Y(j) - lambda(i)*Mx(j);
          
          err = Norm2(Xi);
        }
      
      if (err > threshold*normeX)
        {
          if (rank_proc == 0)
            {
              cout << "Error on eigenvalue " << lambda(i) << " " << lambda_imag(i) << endl;
              cout << "Error = " << err/normeX << endl;
            }
          //abort();
        }
      else
	{
          err_inf = max(err_inf, err/normeX);
	  if (eigen_pair)
	    {
	      lambda0_real(num) = lambda(i);
	      lambda0_imag(num) = lambda_imag(i);
	
	      for (int j = 0; j < N; j++)
		eigen_vec0(j,num) = eigen_vec(j,i);
	      
	      num++;      
	      lambda0_real(num) = lambda(i+1);
	      lambda0_imag(num) = lambda_imag(i+1);
	      
	      for (int j = 0; j < N; j++)
		eigen_vec0(j,num) = eigen_vec(j,i+1);
	      
	      num++;     
	    }
	  else
	    {
	      lambda0_real(num) = lambda(i);
	      lambda0_imag(num) = lambda_imag(i);
	      for (int j = 0; j < N; j++)
		eigen_vec0(j,num) = eigen_vec(j,i);
	      
	      num++;}
	}
      
      if (eigen_pair)
        i += 2;
      else
        i++;
    }
  
  if (rank_proc == 0)
    cout << "Error on accepted eigenvectors = " << err_inf << endl;
  
  lambda0_real.Resize(num);
  lambda0_imag.Resize(num);
  eigen_vec0.Resize(N, num);
}

void UpdateEigenValues(const FemMatrixFreeClass_Base<Complex_wp>& mat_mass,
                       const FemMatrixFreeClass_Base<Complex_wp>& mat_stiff,
		       VectComplex_wp& lambda0_real, VectComplex_wp& lambda0_imag, Matrix<Complex_wp>& eigen_vec0,
		       Real_wp threshold)
{
  int N = mat_stiff.GetM();
  Vector<Complex_wp> lambda(lambda0_real);
  Matrix<Complex_wp> eigen_vec(eigen_vec0);
  
  const DistributedProblem_Base& var_comm = mat_stiff.GetDistributedProblem();
  DistributedVector<Complex_wp> X(var_comm.GetOverlappedDofNumber(), var_comm.comm_group_mode);
  X.Reallocate(N);

  int rank_proc = var_comm.GetRankProcMode();
  DistributedVector<Complex_wp> Y(X), Mx(X), Diff(X);
  X.Zero(); Y.Zero(); Mx.Zero(); Diff.Zero();
  
  int num = 0;
  for (int i = 0; i < lambda0_real.GetM(); i++)
    {
      Real_wp err = 0, normeX(1);
      for (int j = 0; j < N; j++)
	X(j) = eigen_vec(j, i);
          
      mat_stiff.MltVector(X, Y);
      mat_mass.MltVector(X, Mx);
      normeX = Norm2(Mx);
      for (int j = 0; j < N; j++)
	Diff(j) = Y(j) - lambda(i)*Mx(j);
      
      err = Norm2(Diff);
      
      if (err > threshold*normeX)
        {
          if (rank_proc == 0)
            {
              cout << "Error on eigenvalue " << lambda(i) << endl;
              cout << "Error = " << err/normeX << endl;
              //abort();
            }
        }
      else
	{
	  lambda0_real(num) = lambda(i);
	  for (int j = 0; j < N; j++)
	    eigen_vec0(j, num) = eigen_vec(j,i);
	  
	  num++;
	}
    }

  lambda0_real.Resize(num);
  lambda0_imag.Resize(num);
  eigen_vec0.Resize(N, num);
}


void UpdateSecondEigenval(const FemMatrixFreeClass_Base<Real_wp>& mat_mass,
                          const FemMatrixFreeClass_Base<Real_wp>& mat_damp,
                          const FemMatrixFreeClass_Base<Real_wp>& mat_stiff,
                          VectReal_wp& lambda0_real, VectReal_wp& lambda0_imag,
                          Matrix<Real_wp>& eigen_vec0, Real_wp threshold)
{
  cout << "Not implemented" << endl;
  abort();
}

void UpdateSecondEigenval(const FemMatrixFreeClass_Base<Complex_wp>& mat_mass,
                          const FemMatrixFreeClass_Base<Complex_wp>& mat_damp,
                          const FemMatrixFreeClass_Base<Complex_wp>& mat_stiff,
                          VectComplex_wp& lambda0_real, VectComplex_wp& lambda0_imag,
                          Matrix<Complex_wp>& eigen_vec0, Real_wp threshold)
{
  int N = mat_stiff.GetM();
  Vector<Complex_wp> lambda(lambda0_real);
  Matrix<Complex_wp> eigen_vec(N, eigen_vec0.GetN());
  for (int i = 0; i < N; i++)
    for (int j = 0; j < eigen_vec0.GetN(); j++)
      eigen_vec(i, j) = eigen_vec0(i, j);
  
  const DistributedProblem_Base& var_comm = mat_stiff.GetDistributedProblem();
  DistributedVector<Complex_wp> X(var_comm.GetOverlappedDofNumber(), var_comm.comm_group_mode);
  X.Reallocate(N);

  int rank_proc = var_comm.GetRankProcMode();
  DistributedVector<Complex_wp> Y(X), Mx(X), Sx(X), Diff(X);
  X.Zero(); Y.Zero(); Mx.Zero(); Sx.Zero();
  
  int num = 0;
  for (int i = 0; i < lambda0_real.GetM(); i++)
    {
      Real_wp err = 0, normeX(1);
      for (int j = 0; j < N; j++)
	X(j) = eigen_vec(j, i);
          
      mat_stiff.MltVector(X, Y);
      mat_damp.MltVector(X, Sx);
      mat_mass.MltVector(X, Mx);
      normeX = Norm2(Mx);
      for (int j = 0; j < N; j++)
	Diff(j) = Y(j) + lambda(i)*lambda(i)*Mx(j) + lambda(i)*Sx(j);
      
      err = Norm2(Diff);
      
      if (err > threshold*normeX)
        {
          if (rank_proc == 0)
            {
              cout << "Error on eigenvalue " << lambda(i) << endl;
              cout << "Error = " << err/normeX << endl;
              //abort();
            }
        }
      else
	{
	  lambda0_real(num) = lambda(i);
	  for (int j = 0; j < N; j++)
	    eigen_vec0(j, num) = eigen_vec(j,i);
	  
	  num++;
	}
    }

  lambda0_real.Resize(num);
  lambda0_imag.Resize(num);
  eigen_vec0.Resize(N, num);
}

template<class T>
void RetrieveEigenvectors(VirtualEigenProblem<T, T, T>& eigen_solver,
                          const VirtualMatrix<T>& K, Matrix<T>& eigen_vec)
{
  try
    {
      EigenProblemMontjoie<T>& solver = dynamic_cast<EigenProblemMontjoie<T>& >(eigen_solver);
      
      int N = K.GetM();
      int nev = eigen_vec.GetN();
      Matrix<T> eigen_vec0(eigen_vec);
      eigen_vec.Reallocate(N, nev);
      Vector<T> col0(eigen_vec0.GetM()), col(N);
      for (int j = 0; j < nev; j++)
        {
          GetCol(eigen_vec0, j, col0);
          solver.ExpandVector(col0, col, true);
          
          SetCol(col, j, eigen_vec);
        }
      
    }
  catch(...)
    {
    }
}

template<class T>
void ComputeEigenvaluesShift(VirtualEigenProblem<T, T, T>& eigen_solver,
			     const FemMatrixFreeClass_Base<T>& M, const FemMatrixFreeClass_Base<T>& K, int N,
			     InputVariables& data, Vector<T>& lambda, Vector<T>& lambda_imag,
			     Matrix<T>& eigen_vec, string file_eigenvalue, int rank_proc, bool symM, bool drop_negative, 
                             int type_solver = TypeEigenvalueSolver::DEFAULT)
{
  //int lambda_par_shift, bool fill_box
  Vector<T> new_lambda, new_lambda_imag;
  Matrix<T> new_eigenvec;
  VectComplex_wp tabshift;
  const DistributedProblem_Base& var_comm = M.GetDistributedProblem();
  DistributedVector<Complex_wp> x(var_comm.GetOverlappedDofNumber(), var_comm.comm_group_mode);
  x.Reallocate(M.GetM());
  x.Zero();

  lambda.Clear(); lambda_imag.Clear();
  eigen_vec.Clear();
  
  eigen_solver.SetStoppingCriterion(data.tolerance_eigenvalue);
  eigen_solver.SetNbMaximumIterations(data.nb_max_iterations);
  
  if (symM)
    eigen_solver.SetComputationalMode(eigen_solver.SHIFTED_MODE);
  else
    eigen_solver.SetComputationalMode(eigen_solver.INVERT_MODE);
  
  if (type_solver == TypeEigenvalueSolver::SLEPC)
    eigen_solver.SetComputationalMode(eigen_solver.INVERT_MODE);

  if (type_solver == TypeEigenvalueSolver::FEAST)
    eigen_solver.SetComputationalMode(eigen_solver.REGULAR_MODE);
  
  if (!data.fill_box )
    {
      
      for (int p = 0; p < data.shift_cplx.GetM(); p++)
	{
          if (symM)
            eigen_solver.SetComputationalMode(eigen_solver.SHIFTED_MODE);
          else
            eigen_solver.SetComputationalMode(eigen_solver.INVERT_MODE);
          
          if (type_solver == TypeEigenvalueSolver::SLEPC)
            eigen_solver.SetComputationalMode(eigen_solver.INVERT_MODE);
          
          eigen_solver.SetNbAskedEigenvalues(data.number_eigenval(p));
          eigen_solver.SetNbAdditionalEigenvalues(20);
          if (abs(data.shift_cplx(p)) > 1e30)
            {
              if (symM)
                eigen_solver.SetComputationalMode(eigen_solver.REGULAR_MODE);
              else
                eigen_solver.SetComputationalMode(eigen_solver.INVERT_MODE);

              eigen_solver.SetTypeSpectrum(eigen_solver.LARGE_EIGENVALUES,
                                           Complex_wp(0, 0),
                                           eigen_solver.SORTED_MODULUS);
            }
          else
            {
              if (data.rank_proc == 0)
                {
                  DISP(Iwp*data.shift_cplx(p)); DISP(data.number_eigenval(p));              
                }
              
              eigen_solver.SetTypeSpectrum(eigen_solver.CENTERED_EIGENVALUES,
                                           Iwp*data.shift_cplx(p),
                                           eigen_solver.SORTED_MODULUS);
            }
                    
          if (type_solver == TypeEigenvalueSolver::FEAST)
            {
              if (p >= data.radius_shift.GetM())
                {
                  cout << "Missing Radius (RadiusShift) " << endl;
                  DISP(p); DISP(data.radius_shift.GetM());
                }

              FeastParam& param = eigen_solver.GetFeastParameters();
              param.SetCircleSpectrum(Iwp*data.shift_cplx(p), data.radius_shift(p));
              if (data.ratio_shift.GetM() > 0)
                param.SetEllipseSpectrum(Iwp*data.shift_cplx(p), data.radius_shift(p),
                                                data.ratio_shift(p), data.angle_shift(p));
              
              //eigen_solver.SetGlobalPrintLevel(1);
              eigen_solver.SetPrintLevel(1);
              if (data.rank_proc != 0)
                eigen_solver.SetPrintLevel(-1);

              if (data.estimate_number_eigenval)
                eigen_solver.GetFeastParameters().EnableEstimateNumberEigenval();
              else
                eigen_solver.GetFeastParameters().EnableEstimateNumberEigenval(false);

              eigen_solver.GetFeastParameters().SetNumOfQuadraturePoints(data.nb_points_quad_feast);
              eigen_solver.GetFeastParameters().SetTypeIntegration(data.type_integration_feast);
            }
          else
            eigen_solver.SetPrintLevel(3);
          
	  eigen_solver.Init(N);
	  
	  // calcul des valeurs propres
	  if (data.number_eigenval(p) > 0)
            GetEigenvaluesEigenvectors(eigen_solver, new_lambda, new_lambda_imag, new_eigenvec, type_solver);
          else
            {
              new_lambda.Clear();
              new_lambda_imag.Clear();
              new_eigenvec.Clear();
            }
          
          if (data.estimate_number_eigenval)
            exit(0);

          RetrieveEigenvectors(eigen_solver, K, new_eigenvec);
	  UpdateEigenValues(M, K, new_lambda, new_lambda_imag, new_eigenvec, data.threshold_eigen);

	  if (rank_proc == 0)
	    {
	      new_lambda.WriteText(file_eigenvalue + to_str(p) + "_real.dat");
	      new_lambda_imag.WriteText(file_eigenvalue + to_str(p) + "_imag.dat");
              cout << "Number of linear solves = " << eigen_solver.GetNbMatrixVectorProducts() << endl;
	    }
          
	  MergeEigenvalues(new_lambda, new_lambda_imag, new_eigenvec,
			   lambda, lambda_imag, eigen_vec, x, drop_negative);
	}
    }
  else

    {
      VectReal_wp box_real, box_imag;
      Matrix<bool> shift_points_included;
      int nb_real,nb_imag;
      Complex_wp shift_grid = 0;
      
      nb_real = ceil((data.Lr_max-real(Iwp*shift_grid))/data.shift_step);
      nb_imag = ceil((data.Li_max-imag(Iwp*shift_grid))/data.shift_step);

      box_real.Reallocate(nb_real);
      box_imag.Reallocate(nb_imag);
      shift_points_included.Reallocate(nb_real,nb_imag);
      shift_points_included.Fill(false);
      shift_points_included(0,0) = true;
      int nbre_shift_included = 1;

      
      for (int j = 0; j < nb_real; j++)
	box_real(j) = real(Iwp*data.shift_init) + j*data.shift_step;
      for (int k = 0; k < nb_imag; k++)
	box_imag(k) = imag(Iwp*data.shift_init) + k*data.shift_step;
      
      //DISP(box_real);
      //DISP(box_imag);

      Complex_wp shift = data.shift_init;
      int cpt = 0;
      
      while (nbre_shift_included < nb_real * nb_imag)
	{
	  tabshift.PushBack(Iwp*shift);
	  DISP(Iwp*shift); DISP(data.lambda_par_shift);
	  eigen_solver.SetNbAskedEigenvalues(data.lambda_par_shift);
	  eigen_solver.SetNbAdditionalEigenvalues(20);
	  eigen_solver.SetTypeSpectrum(eigen_solver.CENTERED_EIGENVALUES,
				       Iwp*shift,
				       eigen_solver.SORTED_MODULUS);
	  
	  eigen_solver.Init(N);
	  eigen_solver.SetPrintLevel(1);
	  
	  // calcul des valeurs propres
	  GetEigenvaluesEigenvectors(eigen_solver, new_lambda, new_lambda_imag, new_eigenvec);
          RetrieveEigenvectors(eigen_solver, K, new_eigenvec);
	  UpdateEigenValues(M, K, new_lambda, new_lambda_imag, new_eigenvec, data.threshold_eigen);
	  
          // Calcul du rayon
	  Real_wp R = 0;
	  
          VectReal_wp radius(new_lambda.GetM()); int nb = 0;
	  for (int i = 0; i < new_lambda.GetM(); i++)
	    if (imagpart(new_lambda(i)+Iwp*new_lambda_imag(i)) > 0)
              {
                radius(nb) = abs(new_lambda(i)+Iwp*new_lambda_imag(i) - Iwp*shift);
                nb++;
              }

          radius.Resize(nb);
          Sort(radius);
          
          // on enleve les extremes
          R = radius(nb-1-data.nb_extrema_values_to_remove);
          
	  // Detection des shifts inclus	  
	  Real_wp x0,x1,y0,y1;

	  x0 = real(Iwp*shift)-sqrt(2.0)/2.0*R;
	  x1 = real(Iwp*shift)+sqrt(2.0)/2.0*R;
	  y0 = imag(Iwp*shift)-sqrt(2.0)/2.0*R;
	  y1 = imag(Iwp*shift)+sqrt(2.0)/2.0*R;

	  int i0 = ceil((x0-real(Iwp*shift_grid))/data.shift_step);
	  int i1 = floor((x1-real(Iwp*shift_grid))/data.shift_step);
	  int j0 = ceil((y0-imag(Iwp*shift_grid))/data.shift_step);
	  int j1 = floor((y1-imag(Iwp*shift_grid))/data.shift_step);


	  //DISP(i0); DISP(i1); DISP(j0); DISP(j1);
	  
	  for (int i = i0; i <= i1; i++)
	    for (int j = j0; j <= j1; j++)
	      if ((i >= 0) && (i < nb_real) && (j >= 0) && (j < nb_imag))
		{
		  if (!shift_points_included(i,j))
		    {
		      //DISP(i); DISP(j);
		      shift_points_included(i,j) = true;
		      nbre_shift_included++;
		    }
		}
	  
	  //Recherche du nouveau shift
	  //Complex_wp prev_shift = shift;
	  
	  for (int i = 0; i < nb_real; i++)
	    {
	      bool test_arret = false;
	      for (int j = 0; j < nb_imag; j++)
		if (!shift_points_included(i,j))
		  {
		    //cout << "nouveau point" << endl; //DISP(i); DISP(j);
		    shift = -Iwp*(box_real(i) + Iwp*box_imag(j));
		    test_arret = true;
		    break;
		  }
	      
	      if (test_arret)
		break;
	    }
	  
	  //DISP(nbre_shift_included);
	  //DISP(R);
	  //DISP(shift);
	  //DISP(prev_shift);

	  
	  if (rank_proc == 0)
	    {
	      new_lambda.WriteText(file_eigenvalue + to_str(cpt) + "_real.dat");
	      new_lambda_imag.WriteText(file_eigenvalue + to_str(cpt) + "_imag.dat");
	      cpt ++;
	    }
	  
	  MergeEigenvalues(new_lambda, new_lambda_imag, new_eigenvec,
                           lambda, lambda_imag, eigen_vec, x, drop_negative); 

	  //int test; cout << "waiting"<< endl; cin>>test;

	}
      tabshift.Write("liste_shift.dat");
      
    }
  
  if (rank_proc == 0)
    {
      lambda.WriteText(file_eigenvalue + "M_real.dat");
      lambda_imag.WriteText(file_eigenvalue + "M_imag.dat");
    }
  
  // eigen_vec.WriteText(file_eigenvalue + "M_vec.dat");
}


void ComputeEigenvaluesScalapack(DistributedMatrix<Real_wp, General, ArrayRowSparse>& A,
                                 VectReal_wp& lambda, VectReal_wp& lambda_imag,
                                 Matrix<Real_wp>& eigen_vec)
{
  cout << "Coucou je passe la" << endl;
#ifdef SELDON_WITH_SCALAPACK
  // distributed eigensolver
  DistributedMatrix<Real_wp, General, ColMajor> Ah;
  int n = A.GetM();
  int nb = 64;
  Ah.Init(global_blacs_handle, n, n, nb, nb);
  Ah.Zero();
  for (int i = 0; i < A.GetM(); i++)
    {
      for (int j = 0; j < A.GetRowSize(i); j++)
        Ah.SetGlobal(i, A.Index(i, j), A.Value(i, j));   
    }
  
  GetEigenvalues(Ah, lambda, lambda_imag);
#else
  cout << "Recompile with Scalapack" << endl;
  abort();
#endif
}

void ComputeEigenvaluesLapack(DistributedMatrix<Real_wp, General, ArrayRowSparse>& A,
                              VectReal_wp& lambda, VectReal_wp& lambda_imag,
                              Matrix<Real_wp>& eigen_vec)
{
#ifdef SELDON_WITH_MPI
  int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
  if (nb_proc > 1)
    return ComputeEigenvaluesScalapack(A, lambda, lambda_imag, eigen_vec);
#endif
  
  Matrix<Real_wp> Adense;
  Copy(A, Adense);

  GetEigenvaluesEigenvectors(Adense, lambda, lambda_imag, eigen_vec);
}


void ComputeEigenvaluesLapack(DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A,
                              VectComplex_wp& lambda, VectComplex_wp& lambda_imag,
                              Matrix<Complex_wp>& eigen_vec)
{
  Matrix<Complex_wp> Adense;
  Copy(A, Adense);

  GetEigenvaluesEigenvectors(Adense, lambda, eigen_vec);
}


void ComputeEigenvaluesLapack(Matrix<Real_wp>& Adense,
                              VectReal_wp& lambda, VectReal_wp& lambda_imag,
                              Matrix<Real_wp>& eigen_vec)
{
  GetEigenvaluesEigenvectors(Adense, lambda, lambda_imag, eigen_vec);
}


void ComputeEigenvaluesLapack(Matrix<Complex_wp>& Adense,
                              VectComplex_wp& lambda, VectComplex_wp& lambda_imag,
                              Matrix<Complex_wp>& eigen_vec)
{
  GetEigenvaluesEigenvectors(Adense, lambda, eigen_vec);
}

template<class TypeEquation>
class QnmSolver
{
protected:
  typedef typename TypeEquation::Dimension Dimension;
  typedef typename TypeEquation::Complexe Complexe;
  
  // class discretizing the considered equation (Helmholtz or Maxwell)
  EllipticProblem<TypeEquation> vars;
  InputVariables& input_var;

  // linear solver used to solve alpha M_h X + beta K_h X = F
  All_LinearSolver* glob_solver;

  // eigenvalues, eigenvectors and biorthogonal eigenvectors
  VectComplex_wp AllLambda; Vector<bool> eigen_vec_complex;
  Vector<VectComplex_wp> AllEigenVec;

  VectReal_wp diag_Mh, mass_Mh, extra_diag_Mh;
  Vector<int> num_extra_Mh; Complex_wp coef_Mh, coef_Sh;
  FemMatrixFreeClass_Base<Complexe>* Mh_sparse, *Sh_sparse;
  Matrix<Real_wp, Symmetric, ArrayRowSymSparse> sparse_Mh;
  
  Vector<int> num_Es;
  
  // references of the cavity (case of different resonators)
  Vector<int> ref_cavity;

  // pour l'interpolation
  Vector<VectComplex_wp> diff_points_interpolation;
  Vector<SubdivGlobatto> interp;
  Vector<Vector<int> > num_interp;
  
  // matrices pour PML non-diagonales (Maxwell)
  Matrix<Real_wp, Symmetric, ArrayRowSymSparse> Dh_pml, T231_pml, T312_pml;
  
public:
  
  QnmSolver(InputVariables& data) :
    input_var(data)
  { ref_cavity = input_var.ref_cavity; Mh_sparse = NULL; Sh_sparse = NULL;
    SetComplexOne(coef_Mh); SetComplexOne(coef_Sh); }

  ~QnmSolver()
  {
    if (Mh_sparse != NULL)
      delete Mh_sparse;

    if (Sh_sparse != NULL)
      delete Sh_sparse;
  }

  int GetOffsetV();

  void AddIncidentWave(const Real_wp& alpha, const Real_wp& omega, VectComplex_wp& X);

  void ReconstructEigenvector(const VirtualMatrix<Real_wp>& mat_stiff,
                              const Complex_wp& L, VectComplex_wp& x);
  
  bool DiagonalMassMatrix() const;
  bool BlockDiagonalMassMatrix() const;

  Complex_wp GetRatioDeltaEpsilon(const Real_wp& om, const Complex_wp& om_j, int k = 0);

  void LoadEigenvectors(const string& file_name, VectComplex_wp& lambda,
                        Vector<bool>& eigen_val_complex,
                        Vector<VectComplex_wp>& eigenvec)
  {
    string root = GetBaseString(file_name);
    Vector<bool> eigen_vec_cplx;
    VectComplex_wp all_lambda;
    string suffix = "_P" + to_str(input_var.rank_proc) + ".dat";
    eigen_vec_cplx.Read(root + "_cplx" + suffix);
    all_lambda.Read(root + "_lambda" + suffix);

    int num = 0;
    if (input_var.number_eigenval_select.GetM() > 0)
      num = input_var.number_eigenval_select.GetM();
    else
      for (int i = 0; i < all_lambda.GetM(); i++)
        if (input_var.SelectEigenvalue(all_lambda(i)))
          num++;
    
    lambda.Reallocate(num);
    eigen_val_complex.Reallocate(num);
    eigenvec.Reallocate(num);

    num = 0;
    if (input_var.number_eigenval_select.GetM() > 0)
      {
        for (num = 0; num < input_var.number_eigenval_select.GetM(); num++)
          {
            int i = input_var.number_eigenval_select(num); 
            lambda(num) = all_lambda(i);
            eigen_val_complex(num) = eigen_vec_cplx(i);
            eigenvec(num).Read(root + to_str(i) + suffix);
          }
      }
    else
      {
        for (int i = 0; i < all_lambda.GetM(); i++)
          if (input_var.SelectEigenvalue(all_lambda(i)))
            {
              lambda(num) = all_lambda(i);
              eigen_val_complex(num) = eigen_vec_cplx(i);
              eigenvec(num).Read(root + to_str(i) + suffix);
              num++;
            }
      }
  }
  
  void SaveEigenvectors(const string& file_name)
  {
    string root = GetBaseString(file_name);
    string suffix = "_P" + to_str(input_var.rank_proc) + ".dat";
    eigen_vec_complex.Write(root + "_cplx" + suffix);
    AllLambda.Write(root + "_lambda" + suffix);

    int offset_v = this->GetOffsetV();
    VectComplex_wp x(offset_v);
    for (int i = 0; i < AllEigenVec.GetM(); i++)
      {
        for (int j = 0; j < offset_v; j++)
          x(j) = AllEigenVec(i)(j);

        x.Write(root + to_str(i) + suffix);
      }
  }
  
  void ConstructProblem(const Vector<string>& lines_data_file, const string& name_element,
			const string& name_equation)
  {
    // on lit le fichier de donnees
    int rank_proc(0);
#ifdef SELDON_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);    
    input_var.rank_proc = rank_proc;
    input_var.comm = MPI_COMM_WORLD;
#endif
    
    // initialisation des parametres du probleme
    vars.SetTypeEquation(name_equation);
    vars.InitIndices(PhysicalConstant::nb_max_indices);
    ReadInputFile(lines_data_file, vars);
    input_var.UpdateAdim(vars.GetWaveLengthAdim());
    
    if (input_var.dispersive_pml)
      {
        vars.SetFirstOrderFormulation(true);
        vars.linearize_drude = false;
      }
    else
      {
        vars.SetFirstOrderFormulation(false);
        vars.linearize_drude = true;
      }
    
    //DISP(input_var.dispersive_pml);
    //DISP(vars.FirstOrderFormulation());

#ifdef SELDON_WITH_MPI
    if ((input_var.type_eigensolver == TypeEigenvalueSolver::FEAST)
        || (input_var.type_eigensolver == -1))
      vars.comm_group_mode = MPI_COMM_SELF;
#endif
    
    if ((rank_proc != 0) && (input_var.type_eigensolver == TypeEigenvalueSolver::FEAST))
      {
        vars.print_level = -1;
        vars.mesh.print_level = -1;
      }
    
    // initialisation solveur
    glob_solver = vars.GetNewLinearSolver();
    ReadInputFile(lines_data_file, *glob_solver);
    
    bool split_mesh = true;
    if (vars.GetNbProcPerMode() == 1)
      split_mesh = false;

    // on rajoute des elements si demande
    if (input_var.add_pml && (rank_proc == 0))
      {
	vars.mesh.ConstructMesh(0, vars.mesh_data(0));
	vars.mesh_data.Clear();
	for (int i = 0; i < vars.mesh.GetNbElt(); i++)
	  vars.mesh.Element(i).UnsetPML();
	
	// on rerajoute des PML
	for (int i = 0; i < vars.mesh.GetNbPmlAreas(); i++)
	  vars.mesh.GetPmlArea(i).AddPML(i, vars.mesh);		
      }

    // on calcule le maillage et element fini
    vars.ComputeMeshAndFiniteElement(name_element, split_mesh);
    
    if ((vars.GetNbProcPerMode() == 1) && (rank_proc == 0))
      vars.mesh.Write("test.mesh");
    
    // autres initialisations
    vars.PerformOtherInitializations();

    // selection du solveur
    typename TypeEquation::Dimension dim; Symmetric property;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    glob_solver->SelectOptimalLinearSolver(mesh_num.GetOrder(), vars.GetNbDof(), dim, property);

    vars.ComputeMassMatrix(true, true);
    vars.ComputeQuasiPeriodicPhase();
  }  


  // calcule ||u_{pml} || / || u_{all} ||
  template<class T>
  double GetPmlNorm(Vector<T>& eigen_mode)
  {
    double sum_pml = 0, sum_ext = 0;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    for (int i = 0; i < vars.mesh.GetNbElt(); i++)
      {
	const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
	int N   = Fb.GetNbPointsQuadratureInside();

	bool affine = vars.mesh.IsElementAffine(i);
	for (int j = 0; j < N ; j++)
	  {
	    Real_wp jacob = vars.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
	    int n = mesh_num.Element(i).GetNumberDof(j);
	    if (vars.InsidePML(i))
	      sum_pml += jacob*absSquare(eigen_mode(n));
	    else
	      sum_ext += jacob*absSquare(eigen_mode(n));
	  }
      }
  
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(vars.comm_group_mode, &nb_proc);
    if (nb_proc > 1)
      {
        R2 sum_tmp(sum_pml, sum_ext), sum_res;
        MPI_Allreduce(&sum_tmp(0), &sum_res(0), 2, GetMpiDataType(sum_tmp(0)), MPI_MAX, vars.comm_group_mode);
        sum_pml = sum_res(0);
        sum_ext = sum_res(1);
      }
#endif
    
    return sqrt(sum_pml/sum_ext);
  }

  void PerformInitializations(FemMatrixFreeClass_Base<Real_wp>& Mh0, FemMatrixFreeClass_Base<Real_wp>& Kh0);

  void PerformInitializations(FemMatrixFreeClass_Base<Complex_wp>& Mh0, FemMatrixFreeClass_Base<Complex_wp>& Kh0) {}

  // calcule le vecteur w=biortho a partir du vecteur x=eigen_vec
  void ComputeBiorthogonal(const Complex_wp& lambda, const Vector<Complex_wp>& eigen_vec,
                           Vector<Complex_wp>& biortho, bool compute_PQ = true,
                           bool compute_pml = true);
  
  // on traite la matrice de masse comme une matrice diagonale
  void ComputeMassMatrix();
  
  // on traite la matrice de masse comme une matrice diagonale
  void ComputePhysicalMass();

  void MltMass(const VectComplex_wp& x, VectComplex_wp& y)
  {
    y.Reallocate(x.GetM()); y.Zero();
    if (this->DiagonalMassMatrix())
      {
        for (int j = 0; j < diag_Mh.GetM(); j++)
          y(j) = diag_Mh(j)*x(j);
      }
    else if (diag_Mh.GetM() > 0)
      {
        for (int j = 0; j < diag_Mh.GetM(); j++)
          {
            y(j) = diag_Mh(j)*x(j);
            if (num_extra_Mh(j) != -1)
              y(j) += extra_diag_Mh(j)*x(num_extra_Mh(j));
          }
      }
    else
      {
        if (Mh_sparse == NULL)
          {
            cout << "Mh_sparse not available" << endl;
            abort();
          }
        
        Mh_sparse->MltVector(x, y);
        Mlt(coef_Mh, y);
      }    
  }  


  void MltAddDamping(const Complex_wp& alpha, const VectComplex_wp& x, VectComplex_wp& y)
  {
    if (Sh_sparse == NULL)
      {
        cout << "Sh_sparse not available" << endl;
        abort();
      }
    
    Sh_sparse->MltAddVector(alpha*coef_Sh, x, Complex_wp(1, 0), y);
  }
  
  template<class T>
  void CheckOperator(EigenProblemMontjoie<T>* eigen_solver)
  {
    int N = eigen_solver->GetM(); 

    srand(0);

    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

    Vector<Treal> x, y;
    Vector<Tcplx> xc, yc;
    
    GenerateRandomVector(x, N);
    GenerateRandomVector(xc, N);

    y.Reallocate(N);
    yc.Reallocate(N);
    
    eigen_solver->ComputeMassMatrix();

    eigen_solver->MltMass(x, y);
    eigen_solver->MltMass(xc, yc);
    
    x.Write("x.dat"); xc.Write("xc.dat");
    y.Write("y.dat"); yc.Write("yc.dat");

    eigen_solver->MltMass(SeldonTrans, x, y);
    eigen_solver->MltMass(SeldonTrans, xc, yc);

    y.Write("yt.dat"); yc.Write("yt_c.dat");
    
    eigen_solver->ComputeStiffnessMatrix();
    
    eigen_solver->MltStiffness(x, y);
    eigen_solver->MltStiffness(xc, yc);
    
    y.Write("Ky.dat"); yc.Write("Kyc.dat");

    eigen_solver->MltStiffness(SeldonTrans, x, y);
    eigen_solver->MltStiffness(SeldonTrans, xc, yc);
    
    y.Write("Kyt.dat"); yc.Write("Kyt_c.dat");
    
    T a = 0.74, b = -0.27;
    eigen_solver->ComputeStiffnessMatrix(a, b);
    
    eigen_solver->MltStiffness(a, b, x, y);
    eigen_solver->MltStiffness(a, b, xc, yc);
    
    y.Write("Ky_m.dat"); yc.Write("Kyc_m.dat");    

    eigen_solver->ComputeAndFactorizeStiffnessMatrix(a, b);
    
    eigen_solver->ComputeSolution(x, y);

    y.Write("sol_r.dat");

    eigen_solver->ComputeSolution(SeldonTrans, x, y);

    y.Write("solt_r.dat");

    Tcplx ac(0.62,2.1), bc(-0.31,0.93);
    eigen_solver->ComputeAndFactorizeStiffnessMatrix(ac, bc);
    
    eigen_solver->ComputeSolution(xc, yc);

    yc.Write("sol_c.dat");

    eigen_solver->ComputeSolution(SeldonTrans, xc, yc);

    yc.Write("solt_c.dat");

    
  }

  
  void ReconstructEigenvectors()
  {
    FemMatrixFreeClass_Base<Real_wp>* Kh0;
    
    Kh0 = vars.GetNewIterativeMatrix(Real_wp(0));
    
    GlobalGenericMatrix<Real_wp> nat_mat0;
    nat_mat0.SetCoefMass(0.0);
    nat_mat0.SetCoefStiffness(1.0);
    nat_mat0.SetCoefDamping(1.0);
    vars.SetCoefficientDirichlet(Real_wp(0));
    vars.AddMatrixWithBC(*Kh0, nat_mat0);
    
    for (int i = 0; i < AllEigenVec.GetM(); i++)
      this->ReconstructEigenvector(*Kh0, AllLambda(i), AllEigenVec(i));
    
    delete Kh0;
  }

  void ComputeRationalProblem(Vector<DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse> >& vec_Ai,
                              Vector<Vector<Complex_wp> >& coef_num, Vector<Vector<Complex_wp> >& coef_denom)
  {
    cout << "Not implemented in the general case " << endl;
    abort();
  }


  void ComputeRationalProblem(Vector<DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> >& vec_Ai,
                              Vector<Vector<Real_wp> >& coef_num, Vector<Vector<Real_wp> >& coef_denom)
  {
    cout << "Not implemented in the general case " << endl;
    abort();
  }

  template<class T>
  void SolveRationalEigenProblem(const FemMatrixFreeClass_Base<T>& Mh,
                                 const FemMatrixFreeClass_Base<T>& Sh,
                                 const FemMatrixFreeClass_Base<T>& Kh,
                                 Vector<T>& lambda, Vector<T>& lambda_imag,
                                 Matrix<T>& eigen_vec)
  {
    cout << "not possible" << endl;
    abort();
  }

  void SolveRationalEigenProblem(const FemMatrixFreeClass_Base<Complexe>& Mh0,
                                 const FemMatrixFreeClass_Base<Complexe>& Sh0,
                                 const FemMatrixFreeClass_Base<Complexe>& Kh0,
                                 Vector<Complexe>& lambda,
                                 Vector<Complexe>& lambda_imag,
                                 Matrix<Complexe>& eigen_vec)
  {
    cout << "not implemented" << endl;
    abort();    
  }
  
  template<class T>
  void SolvePolynomialEigenProblem(const FemMatrixFreeClass_Base<T>& Mh,
                                   const FemMatrixFreeClass_Base<T>& Sh,
                                   const FemMatrixFreeClass_Base<T>& Kh,
                                   Vector<T>& lambda, Vector<T>& lambda_imag,
                                   Matrix<T>& eigen_vec)
  {
    cout << "not possible" << endl;
    abort();
  }

  void SolvePolynomialEigenProblem(const FemMatrixFreeClass_Base<Complexe>& Mh0,
                                   const FemMatrixFreeClass_Base<Complexe>& Sh0,
                                   const FemMatrixFreeClass_Base<Complexe>& Kh0,
                                   Vector<Complexe>& lambda,
                                   Vector<Complexe>& lambda_imag,
                                   Matrix<Complexe>& eigen_vec)
  {
    PolynomialEigenProblemMontjoie<Complexe>* eigen_solver = NULL;
    eigen_solver = vars.GetNewPolynomialEigenSolver(*glob_solver);
    
    EigenvalueComparisonMode<Complexe> compar;
    if (input_var.avoid_poles)
      {
        VectComplex_wp poles;
        for (int ref = 1; ref <= 10; ref++)
          if (vars.ref_drude(ref).IsEnabled())
            {
              VectComplex_wp P = vars.ref_drude(ref).GetPoleOmega();
              poles.PushBack(P);
              if (input_var.avoid_zeros)
                {
                  P = vars.ref_drude(ref).GetZeroOmega();
                  poles.PushBack(P);                  
                }
            }        
        
        DISP(poles);
        compar.SetPoles(poles, input_var.radius_poles);
      }
    
    eigen_solver->SetStoppingCriterion(input_var.tolerance_eigenvalue);
    eigen_solver->SetNbMaximumIterations(input_var.nb_max_iterations);
    //SlepcParamPep& param = eigen_solver->GetSlepcParameters();
    //param.SetEigensolverType(SlepcParamPep::QARNOLDI);
    
    const DistributedProblem_Base& var_comm = Mh0.GetDistributedProblem();
    DistributedVector<Complex_wp> x(var_comm.GetOverlappedDofNumber(), var_comm.comm_group_mode);
    x.Reallocate(Mh0.GetM());
    x.Zero();
    
    eigen_solver->SetSpectralTransformation();
    Matrix<Complexe> new_eigenvec;
    Vector<Complexe> new_lambda, new_lambda_imag;
    for (int i = 0; i < input_var.shift_cplx.GetM(); i++)
      {
        Complexe shift; to_complex(input_var.shift_cplx(i)/vars.GetOmega(), shift);
        eigen_solver->SetNbAskedEigenvalues(input_var.number_eigenval(i));
        
        if (input_var.avoid_poles)
          {
            compar.SetShift(vars.GetOmega(), shift);
            eigen_solver->SetUserComparisonClass(&compar);
            eigen_solver->SetTypeSpectrum(eigen_solver->CENTERED_EIGENVALUES,
                                          shift, eigen_solver->SORTED_USER);
          }
        else
          eigen_solver->SetTypeSpectrum(eigen_solver->CENTERED_EIGENVALUES,
                                        shift, eigen_solver->SORTED_MODULUS);
        
        eigen_solver->SetPrintLevel(3);
        eigen_solver->InitComputation();
        
        Matrix<Complexe, General, ColMajor> V;
#ifdef SELDON_WITH_SLEPC
        FindEigenvaluesSlepc(*eigen_solver, new_lambda, new_lambda_imag, V);
#else
        cout << "Recompile Montjoie with Slepc" << endl;
        abort();
#endif
        
        DISP(new_lambda); DISP(new_lambda_imag);
        int nev = input_var.number_eigenval(i);
        new_eigenvec.Reallocate(Kh0.GetM(), nev);
        Vector<Complexe> col0(V.GetM()), col(Kh0.GetM());
        for (int j = 0; j < nev; j++)
          {
            GetCol(V, j, col0);
            eigen_solver->ExpandVector(col0, col, true);
            
            SetCol(col, j, new_eigenvec);
          }
        
        UpdateSecondEigenval(Mh0, Sh0, Kh0, new_lambda, new_lambda_imag, new_eigenvec,
                             input_var.threshold_eigen);

        if (input_var.rank_proc == 0)
          new_lambda.Write(input_var.file_name_eigen + to_str(i) + "_tmp.dat");
        
        MergeEigenvalues(new_lambda, new_lambda_imag, new_eigenvec,
                         lambda, lambda_imag, eigen_vec, x, false,
                         input_var.drop_negative_modes);
      }
  }

  template<class T>
  void ComputeEigenvectors(Vector<T>& lambda, Vector<T>& lambda_imag,
                           Matrix<T>& eigen_vec)
  {
    FemMatrixFreeClass_Base<T>* Kh, *Kh0, *Sh0;
    FemMatrixFreeClass_Base<T>* Mh, *Mh0;
    
    Kh = vars.GetNewIterativeMatrix(T(0));
    Mh = vars.GetNewIterativeMatrix(T(0));
    Kh0 = vars.GetNewIterativeMatrix(T(0));
    Mh0 = vars.GetNewIterativeMatrix(T(0));
    Sh0 = vars.GetNewIterativeMatrix(T(0));
    
    GlobalGenericMatrix<T> nat_mat0;
    nat_mat0.SetCoefMass(1.0);
    nat_mat0.SetCoefStiffness(0.0);
    nat_mat0.SetCoefDamping(0.0);
    vars.SetCoefficientDirichlet(Real_wp(1));
    vars.AddMatrixWithBC(*Mh0, nat_mat0);
	
    nat_mat0.SetCoefMass(0.0);
    nat_mat0.SetCoefStiffness(1.0);
    if (vars.FirstOrderFormulation())
      nat_mat0.SetCoefDamping(1.0);
    else
      nat_mat0.SetCoefDamping(0.0);
    
    vars.SetCoefficientDirichlet(Real_wp(0));
    vars.AddMatrixWithBC(*Kh0, nat_mat0);
    
    //DISP(vars.FirstOrderFormulation());

    if (!vars.FirstOrderFormulation())
      {
        nat_mat0.SetCoefMass(0.0);
        nat_mat0.SetCoefStiffness(0.0);
        nat_mat0.SetCoefDamping(1.0);
        
        vars.SetCoefficientDirichlet(Real_wp(0));
        vars.AddMatrixWithBC(*Sh0, nat_mat0);
      }
    
    ComputeMassMatrix();
    
    this->PerformInitializations(*Mh0, *Kh0);
    
    //DISP(vars.FirstOrderFormulation());
    if (input_var.explicit_matrix)
      {
	// cas ou on stocke les matrices
	DistributedMatrix<T, General, ArrayRowSparse> Md;
        DistributedMatrix<T, General, ArrayRowSparse> K;
        DistributedMatrix<T, General, ArrayRowSparse> S;
        
	// calcul matrice de masse
	GlobalGenericMatrix<T> nat_mat;
	nat_mat.SetCoefMass(1.0);
	nat_mat.SetCoefStiffness(0.0);
	nat_mat.SetCoefDamping(0.0);
	vars.SetCoefficientDirichlet(Real_wp(1));
	
	Md.Clear();
	vars.AddMatrixWithBC(Md, nat_mat);
	
	Md.WriteText("Mh.dat");
	
	// calcul matrice de rigidite
	nat_mat.SetCoefMass(0.0);
	nat_mat.SetCoefStiffness(1.0);
	if (vars.FirstOrderFormulation())
          nat_mat.SetCoefDamping(1.0);
        else
          nat_mat.SetCoefDamping(0.0);
	vars.SetCoefficientDirichlet(Real_wp(0));
	
	K.Clear();
	vars.AddMatrixWithBC(K, nat_mat);  
        
	K.WriteText("Kh.dat");

	if (!vars.FirstOrderFormulation())
          {
            nat_mat.SetCoefMass(0.0);
            nat_mat.SetCoefStiffness(0.0);
            nat_mat.SetCoefDamping(1.0);
            vars.SetCoefficientDirichlet(Real_wp(0));
	
            S.Clear();
            vars.AddMatrixWithBC(S, nat_mat);  
        
            S.WriteText("Sh.dat");            
          }

        if ((!vars.FirstOrderFormulation()) && (input_var.type_eigensolver != -1))
          {
#ifdef SELDON_WITH_SLEPC
            Vector<DistributedMatrix<T, Symmetric, ArrayRowSymSparse> > vec_Ai;
            Vector<Vector<T> > coef_num, coef_denom; 
            ComputeRationalProblem(vec_Ai, coef_num, coef_denom);

            for (int i = 0; i < vec_Ai.GetM(); i++)
              {
                vec_Ai(i).WriteText("Ai" + to_str(i) + ".dat");
                DISP(i); DISP(coef_num(i)); DISP(coef_denom(i));
              }

            SplitSparseNonLinearEigenProblem<T, Symmetric, ArrayRowSymSparse> var_eig;
            var_eig.InitMatrix(vec_Ai, coef_num, coef_denom);

            var_eig.SetStoppingCriterion(1e-12);
            var_eig.SetNbAskedEigenvalues(10);

            T one; to_complex(Complex_wp(0.4, -0.1), one);
            var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, one);
            SlepcParamNep& param = var_eig.GetSlepcParameters();
            param.SetEigensolverType(param.NLEIGS);

            param.SetIntervalRegion(0.3, 0.6, -0.3, -1e-4);
            param.EnableCommandLineOptions();
            
            Matrix<T, General, ColMajor> eigen_vec2;
            FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec2);
            DISP(lambda);
#else
            cout << "Recompile Montjoie with Slepc" <<endl;
            abort();
#endif
            
            exit(0);
          }
                
        DistributedMatrix<T, General, ArrayRowSparse> Ah;
        bool generalized_eigenvalue_pb = false;
        if (this->DiagonalMassMatrix())
          {
            // on suppose ici la matrice de masse diagonale
            Vector<T> M(Md.GetM());
            for (int i = 0; i < Md.GetM(); i++)
              M(i) = Real_wp(1) / Md(i, i);
            
            //Md.WriteText("Mh.dat");
        
            //int N = K.GetM();
            Ah = K;
            ScaleLeftMatrix(Ah, M);
          }
        else if (this->BlockDiagonalMassMatrix())
          {
            // cas block-diagonal
            DistributedMatrix<T, General, BlockDiagRow> Mb;
            
            ConvertToBlockDiagonal(Md, Mb);
            GetInverse(Mb);
            ConvertToSparse(Mb, Md);
            Md.WriteText("Mb.dat");

            if (!vars.FirstOrderFormulation())
              {
                DistributedMatrix<T, General, ArrayRowSparse> B, C;
                B.Init(K);
                Mlt(Md, K, B);
                C.Init(S);
                Mlt(Md, S, C);

                int N = K.GetM();
                Ah.Reallocate(2*N, 2*N);
                for (int i = 0; i < N; i++)
                  {
                    Ah.ReallocateRow(i, 1);
                    Ah.Index(i, 0) = i+N;
                    Ah.Value(i, 0) = 1.0;

                    int size_row = B.GetRowSize(i) + C.GetRowSize(i);
                    Ah.ReallocateRow(i+N, size_row);
                    for (int j = 0; j < B.GetRowSize(i); j++)
                      {
                        Ah.Index(i+N, j) = B.Index(i, j);
                        Ah.Value(i+N, j) = -B.Value(i, j);
                      }
                    
                    int pos = B.GetRowSize(i);
                    for (int j = 0; j < C.GetRowSize(i); j++)
                      {
                        Ah.Index(i+N, pos+j) = N+C.Index(i, j);
                        Ah.Value(i+N, pos+j) = -C.Value(i, j);
                      }
                  }
              }
            else
              {
                Ah.Init(K);
                Mlt(Md, K, Ah);
              }
          }
        else
          { 
            generalized_eigenvalue_pb = true;
          }
        
	if (!generalized_eigenvalue_pb)
          Ah.WriteText("Ah.dat");
        
	if (input_var.shift_cplx.GetM() == 0)
	  {
	    input_var.shift_cplx.Reallocate(1);
	    input_var.shift_cplx(0) = Iwp*vars.GetOmega();
	  }
        
        //DISP(generalized_eigenvalue_pb); DISP(K.GetM());
        if (input_var.type_eigensolver == -1)
          {
            if (generalized_eigenvalue_pb)
              {
                SparseDistributedSolver<T> mat_lu;
                mat_lu.Factorize(Md);
                
                Vector<T, VectSparse> col_sparse;
                DistributedMatrix<T, General, ArrayRowSparse> A;
                int N = K.GetM();
                if (vars.FirstOrderFormulation())
                  A.Reallocate(N, N);
                else
                  {
                    A.Reallocate(2*N, 2*N);
                    for (int i = 0; i < N; i++)
                      {
                        A.ReallocateRow(i, 1);
                        A.Index(i, 0) = i+N;
                        A.Value(i, 0) = 1.0;
                      }
                  }
                
                Matrix<T, General, ColMajor> x_sol;
                Vector<T> col(K.GetM());
                for (int i = 0; i < K.GetM(); i += 100)
                  {
                    int size = 100;
                    if (i+size > K.GetM())
                      size = K.GetM() - i;
                    
                    x_sol.Reallocate(K.GetM(), size);
                    for (int k = 0; k < size; k++)
                      {
                        GetCol(K, i+k, col_sparse);
                        col.Zero();
                        for (int j = 0; j < col_sparse.GetM(); j++)
                          col(col_sparse.Index(j)) = col_sparse.Value(j);
                        
                        SetCol(col, k, x_sol);
                      }

                    mat_lu.Solve(x_sol);
                    
                    for (int k = 0; k < size; k++)
                      {
                        GetCol(x_sol, k, col);
                        for (int j = 0; j < K.GetM(); j++)
                          if (abs(col(j)) > 1e-30)
                            {
                              if (!vars.FirstOrderFormulation())
                                A.AddInteraction(j+N, i+k, -col(j));
                              else
                                A.AddInteraction(j, i+k, col(j));
                            }
                      }
                    
                    if (!vars.FirstOrderFormulation())
                      {
                        // part with -S
                        for (int k = 0; k < size; k++)
                          {
                            GetCol(S, i+k, col_sparse);
                            col.Zero();
                            for (int j = 0; j < col_sparse.GetM(); j++)
                              col(col_sparse.Index(j)) = col_sparse.Value(j);
                            
                            SetCol(col, k, x_sol);
                          }
                        
                        mat_lu.Solve(x_sol);
                        
                        for (int k = 0; k < size; k++)
                          {
                            GetCol(x_sol, k, col);
                            for (int j = 0; j < N; j++)
                              if (abs(col(j)) > 1e-30)
                                A.AddInteraction(j+N, N+i+k, -col(j));
                          }
                      }
                  }
                
                cout << "computing eigenvalues of matrix of size " << A.GetM() << endl;
                //A.WriteText("Ah.dat");
                
                ComputeEigenvaluesLapack(A, lambda, lambda_imag, eigen_vec);
              }
            else
              {
                cout << "computing eigenvalues of matrix of size " << Ah.GetM() << endl;
                ComputeEigenvaluesLapack(Ah, lambda, lambda_imag, eigen_vec);
              }
            
            lambda.Write("LambdaLapack.dat");
            lambda_imag.Write("LambdaImagLapack.dat");
            
            if (vars.FirstOrderFormulation())
              UpdateEigenValues(*Mh0, *Kh0, lambda, lambda_imag, eigen_vec, input_var.threshold_eigen);
            else
              UpdateSecondEigenval(*Mh0, *Sh0, *Kh0, lambda, lambda_imag, eigen_vec, input_var.threshold_eigen);
          }
        else
          {
            // on cherche les valeurs propres proches de sigma = shift_cplx
            // les valeurs propres sont triees par module croissant            
            if (generalized_eigenvalue_pb)
              {
                cout << "Not available" << endl; abort();
                //SparseEigenProblem<T, DistributedMatrix<T, General, ArrayRowSparse>,
                //                   DistributedMatrix<T, General, ArrayRowSparse> > eigen_solver_sp;

                //eigen_solver_sp.InitMatrix(K, Md);	
                
                //ComputeEigenvaluesShift(eigen_solver_sp, *Mh0, *Kh0, K.GetM(), input_var,
                //                        lambda, lambda_imag, eigen_vec, input_var.file_name_eigen,
                //                        vars.GetRankProcMode(), false, input_var.type_eigensolver);
              }
            else if (!vars.FirstOrderFormulation())
              {
                cout << "Not implemented " << endl;
                abort();
              }
            else
              {
                SparseEigenProblem<T, DistributedMatrix<T, General, ArrayRowSparse>,
                                   DistributedMatrix<T, Symmetric, ArrayRowSymSparse> > eigen_solver_sp;
                
                eigen_solver_sp.InitMatrix(Ah);	
                
                ComputeEigenvaluesShift(eigen_solver_sp, *Mh0, *Kh0, K.GetM(), input_var,
                                        lambda, lambda_imag, eigen_vec, input_var.file_name_eigen,
                                        vars.GetRankProcMode(), true, input_var.drop_negative_modes, input_var.type_eigensolver);
              }
          }
      }
    else if (!vars.FirstOrderFormulation())
      {
        if (input_var.use_rational_eigensolver)
          SolveRationalEigenProblem(*Mh0, *Sh0, *Kh0, lambda, lambda_imag, eigen_vec);
        else
          SolvePolynomialEigenProblem(*Mh0, *Sh0, *Kh0, lambda, lambda_imag, eigen_vec);
      }
    else
      {
	// initialisation du probleme aux valeurs propres
	EigenProblemMontjoie<T>* eigen_solver = NULL;
	eigen_solver = vars.GetNewEigenSolver(*glob_solver, T());
        
	eigen_solver->InitMatrix(*Kh, *Mh);
        
	if (input_var.rank_proc == 0)
          cout << "On appelle ComputeEigenvaluesShift" << endl;
        
        if (input_var.type_eigensolver == -1)
          {
            cout << "Lapack solver must be selected with ExplicitMatrixEigen = YES" << endl;
            abort();
          }
        
        int print_level = vars.print_level;
        if (input_var.type_eigensolver == TypeEigenvalueSolver::FEAST)
          vars.print_level = -1;
        
        bool symM = this->DiagonalMassMatrix();
        
	ComputeEigenvaluesShift(*eigen_solver, *Mh0, *Kh0, eigen_solver->GetM(), input_var,
				lambda, lambda_imag, eigen_vec, input_var.file_name_eigen,
				vars.GetRankProcMode(), symM, input_var.drop_negative_modes, input_var.type_eigensolver);

        vars.print_level = print_level;
      }
    
    delete Mh; delete Kh; delete Kh0; delete Mh0; delete Sh0;
  }

  
  void ExtractEigenvectors(const VectReal_wp& lambda, const VectReal_wp& lambda_imag,
                           const Matrix<Real_wp>& eigen_vec)
  {
    //DISP(lambda); DISP(lambda_imag);
    //lambda.Write("lambda_real.dat");
    //lambda_imag.Write("lambda_imag.dat");

    int nev = lambda.GetM();
    //int N = vars.GetNbDof();
    
    // we count the number of eigenvalues that are kept
    int num = 0, nb = 0;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    int Nu = mesh_num.GetNbDof();
    VectComplex_wp x(eigen_vec.GetM());
    while (nb < nev)
      {
	if (lambda_imag(nb) == Real_wp(0))
	  {
            for (int j = 0; j < x.GetM(); j++)
              x(j) = Complex_wp(eigen_vec(j, nb), 0);
            
            if (input_var.SelectEigenvalue(Complex_wp(lambda(nb), 0))
                && input_var.SelectEigenvector(x, Nu))
              num++;
            
	    nb++;
	  }
	else if ((nb < nev-1) && input_var.
		 SelectEigenvalue(Complex_wp(lambda(nb), lambda_imag(nb))))
	  {
            for (int j = 0; j < x.GetM(); j++)
              x(j) = Complex_wp(eigen_vec(j, nb), eigen_vec(j, nb+1));
            
	    if (input_var.SelectEigenvector(x, Nu))
              num++;
            
	    nb += 2;
	  }
	else
	  nb += 2;
      }

    AllLambda.Reallocate(num);
    AllEigenVec.Reallocate(num);
    eigen_vec_complex.Reallocate(num);

    // boucle sur les valeurs propres (considerees complexes conjuguees)
    nb = 0; num = 0;
    while (nb < nev)
      {
	if (lambda_imag(nb) == Real_wp(0))
	  {
	    // real mode
            if (input_var.SelectEigenvalue(Complex_wp(lambda(nb), 0)))
              {
                for (int j = 0; j < x.GetM(); j++)
                  x(j) = Complex_wp(eigen_vec(j, nb), 0);

                if (input_var.SelectEigenvector(x, Nu))
                  {
                    eigen_vec_complex(num) = false;                
                    AllLambda(num) = Complex_wp(lambda(nb), 0);
                    AllEigenVec(num) = x;
                    num++;
                  }
              }
            
	    nb++;
	  }
	else if ((nb < nev-1) && input_var.
		 SelectEigenvalue(Complex_wp(lambda(nb), lambda_imag(nb))))
	  {
            for (int j = 0; j < x.GetM(); j++)
              x(j) = Complex_wp(eigen_vec(j, nb), eigen_vec(j, nb+1));

            if (input_var.SelectEigenvector(x, Nu))
              {
                // we store only an eigenvalue (the other being complex conjugate)
                AllLambda(num) = Complex_wp(lambda(nb), lambda_imag(nb));
                eigen_vec_complex(num) = true;
                AllEigenVec(num) = x;
                
                // next eigenvalue
                num++;
              }
            
            if (lambda_imag(nb) != -lambda_imag(nb+1))
              {
                cout << "Eigenvalues not complex conjugate ?" << endl;
                abort();
              }
            
	    nb += 2;
	  }
	else
	  {
	    // modes outside spectrum are not stored
	    nb += 2;
	  }
      }
  }


  void ExtractEigenvectors(const VectComplex_wp& lambda, const VectComplex_wp& lambda_imag,
                           const Matrix<Complex_wp>& eigen_vec)
  {
    int nev = lambda.GetM();
    //int N = vars.GetNbDof();
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    int Nu = mesh_num.GetNbDof();
    VectComplex_wp x(eigen_vec.GetM());
    
    // we count the number of eigenvalues that are kept
    int num = 0;
    for (int i = 0; i < nev; i++)
      if (input_var.SelectEigenvalue(lambda(i)))
        {
          GetCol(eigen_vec, i, x);
          if (input_var.SelectEigenvector(x, Nu))
            num++;
        }

    AllLambda.Reallocate(num);
    AllEigenVec.Reallocate(num);
    eigen_vec_complex.Reallocate(num);
    eigen_vec_complex.Fill(true);
    num = 0;
    for (int i = 0; i < nev; i++)
      if (input_var.SelectEigenvalue(lambda(i)))
        {
          GetCol(eigen_vec, i, x);
          if (input_var.SelectEigenvector(x, Nu))
            {
              AllLambda(num) = lambda(i);
              AllEigenVec(num) = x;
              if (abs(imagpart(lambda(i))) <= 1e-12)
                eigen_vec_complex(num) = false;
              
              num++;
            }
        }
  }
  
  void ExtractQnmModes()
  {
    Vector<bool> is_mode_qnm(AllLambda.GetM());
    VectReal_wp ratio_qnm(AllLambda.GetM());
    is_mode_qnm.Fill(false);
    int Nb_QNM = 0;
    for (int i = 0; i < AllLambda.GetM(); i++)
      {
	double ratio = GetPmlNorm(AllEigenVec(i));
	ratio_qnm(i) = ratio;
	if (ratio < input_var.threshold_ratio_qnm)
	  {
	    is_mode_qnm(i) = true;
	    Nb_QNM++;
	    cout << "MODE " << i << " est a priori QNM" << endl;
	  }
      }
    
    ratio_qnm.WriteText("ratio.dat");
    if (input_var.extractQNM)
      {
	VectComplex_wp oldlambda(AllLambda);
	Vector<VectComplex_wp> oldEigenvec(AllEigenVec);
        Vector<bool> old_eigen_vec_cplx(eigen_vec_complex);
	AllLambda.Reallocate(Nb_QNM);
	AllEigenVec.Reallocate(Nb_QNM);
        eigen_vec_complex.Reallocate(Nb_QNM);
	Nb_QNM = 0; 

	for (int i = 0; i < oldlambda.GetM(); i++)
	  if (is_mode_qnm(i))
	    {
	      AllLambda(Nb_QNM) = oldlambda(i);
	      AllEigenVec(Nb_QNM) = oldEigenvec(i);
              eigen_vec_complex(Nb_QNM) = old_eigen_vec_cplx(i);
	      Nb_QNM ++;
	    }
        
	if ((vars.GetRankProcMode() == 0) && (input_var.rank_proc == 0))
	  AllLambda.Write("LambdaQNM.dat");
      }
  }
  
  void ComputeEigenmodes()
  {
    vars.SetHomogeneousDirichlet(true);
    if (input_var.type_eigensolver == -1)
      input_var.explicit_matrix = true;
    
    if (input_var.dispersive_pml && (input_var.type_eigensolver != TypeEigenvalueSolver::SLEPC))
      {
        VectReal_wp lambda, lambda_imag;
        Matrix<Real_wp> eigen_vec;
        
        ComputeEigenvectors(lambda, lambda_imag, eigen_vec);
        ExtractEigenvectors(lambda, lambda_imag, eigen_vec);
      }
    else
      {
        VectComplex_wp lambda, lambda_imag;
        Matrix<Complex_wp> eigen_vec;

        ComputeEigenvectors(lambda, lambda_imag, eigen_vec);
        if (!input_var.dispersive_pml)
          {
            Mlt(Iwp*vars.GetOmega(), lambda);
            coef_Mh = -1.0/vars.GetSquareOmega();
            coef_Sh = 1.0/(-Iwp*vars.GetOmega());
          }
        
        ExtractEigenvectors(lambda, lambda_imag, eigen_vec);
      }

    // on trie les valeurs propres
    {
      VectR2 lambda_vec(AllLambda.GetM());
      for (int i = 0; i < AllLambda.GetM(); i++)
        lambda_vec(i).Init(realpart(AllLambda(i)), imagpart(AllLambda(i)));

      VectComplex_wp oldlambda(AllLambda);
      Vector<VectComplex_wp> oldEigenvec(AllEigenVec);
      Vector<bool> oldEigenVecComplex(eigen_vec_complex);

      Vector<int> permut(AllLambda.GetM());
      permut.Fill();
      Real_wp threshold = R2::threshold;
      R2::threshold = input_var.threshold_eigenval_equal;
      Sort(lambda_vec, permut);
      R2::threshold = threshold;

      for (int i = 0; i < oldlambda.GetM(); i++)
        {
          AllLambda(i) = oldlambda(permut(i));
          AllEigenVec(i) = oldEigenvec(permut(i));
          eigen_vec_complex(i) = oldEigenVecComplex(permut(i));
        }
    }
    
    if ((vars.GetRankProcMode() == 0) && (input_var.rank_proc == 0))
      AllLambda.Write(input_var.file_name_eigen + "_sorted.dat");    
    
    for (int i = 0; i < AllLambda.GetM(); i++)
      {      
        //AllEigenVec(i).Write("XmodeI" + to_str(i) + ".dat");
        
	// Pour afficher le ratio pour chaque mode
	//double ratio = GetPmlNorm(AllEigenVec(i));
	//DISP(i);
	//DISP(ratio);
	
      }

    if (input_var.extractQNM)
      ExtractQnmModes();
    
    // Gram-schmidt orthogonalization
    DistributedVector<Complex_wp> biortho(vars.GetOverlappedDofNumber(), vars.comm_group_mode);
    DistributedVector<Complex_wp> y(vars.GetOverlappedDofNumber(), vars.comm_group_mode);
    int num = 0; int nmax = 0; int sum_multiple = 0;
    VectReal_wp err_multiple;
    while (num < AllLambda.GetM())
      {
	int n = 0;
	while ((num+n < AllLambda.GetM()) && (abs(AllLambda(num)-AllLambda(num+n)) <= input_var.threshold_eigenval_equal))
	  n++;
        
	if (n > 1)
	  {
            Real_wp errL(0);
            for (int j = 0; j < n; j++)
              for (int k = 0; k < n; k++)
                errL = max(errL, abs(AllLambda(num+j) - AllLambda(num+k)));
            
            err_multiple.PushBack(errL);
            nmax = max(nmax, n);
            sum_multiple += n;
            vector<DistributedVector<Complex_wp> > vec_biortho(n, y);
	    // multiple eigenvalues, we orthogonalize with Gram-Schmidt procedure
            for (int i = 0; i < n; i++)
              {
                VectComplex_wp x = AllEigenVec(num+i);
                MltMass(x, y);
                if (!vars.FirstOrderFormulation())
                  {
                    Mlt(Real_wp(2), y);
                    MltAddDamping(Real_wp(1) / AllLambda(num+i), x, y);
                  }
                
                for (int j = 0; j < i; j++)
                  {
                    Complex_wp alpha = DotProd(y, vec_biortho[j]);
                    Add(-alpha, AllEigenVec(num+j), x);
                  }
                
                ComputeBiorthogonal(AllLambda(num+i), x, vec_biortho[i]);
                Complex_wp scal = DotProd(y, vec_biortho[i]);
                
                if (scal == Complex_wp(0, 0))
                  {
                    x.Zero();
                    vec_biortho[i].Zero();
                  }
                else
                  {
                    Mlt(Real_wp(1)/sqrt(scal), x);
                    Mlt(Real_wp(1)/sqrt(scal), vec_biortho[i]);
                  }

                AllEigenVec(num+i) = x;
              }
	  }

	num += n;
      }

    if (input_var.rank_proc == 0)
      {
        cout << "Maximum multiplicity of eigenvalues = " << nmax << endl;
        cout << "Total number of multiple eigenvalues = " << sum_multiple << endl;
        DISP(err_multiple); DISP(err_multiple.GetNormInf());
      }
        
    // normalization
    if (input_var.rank_proc == 0)
      DISP(AllLambda.GetM());
    
    for (int i = 0; i < AllLambda.GetM(); i++)
      {
	// biorthogonal vector is computed
        ComputeBiorthogonal(AllLambda(i), AllEigenVec(i), biortho);
        
        MltMass(AllEigenVec(i), y);
        if (!vars.FirstOrderFormulation())
          {
            Mlt(Real_wp(2), y);
            MltAddDamping(Real_wp(-1) / AllLambda(i), AllEigenVec(i), y);
          }
        
	Complex_wp scal = DotProd(y, biortho);
        
	// normalisation du mode pour avoir x_j \cdot w_j = 1
        if (scal == Complex_wp(0, 0))
          {
            AllEigenVec(i).Zero();
            biortho.Zero();
          }
        else
          {
            scal = 1.0/sqrt(scal);
            Mlt(scal, AllEigenVec(i));
            Mlt(scal, biortho);
          }

        //AllEigenVec(i).Write("Xmode" + to_str(i) + ".dat");
	//biortho.Write("Xortho"+ to_str(i) + ".dat");
      }

    // on teste ReconstructEigenvector
    /*    if (false)
      {
        FemMatrixFreeClass_Base<Real_wp>* Kh0;
        
        Kh0 = vars.GetNewIterativeMatrix(Real_wp(0));
        
        GlobalGenericMatrix<Real_wp> nat_mat0;
        nat_mat0.SetCoefMass(0.0);
        nat_mat0.SetCoefStiffness(1.0);
        nat_mat0.SetCoefDamping(1.0);
        vars.SetCoefficientDirichlet(Real_wp(0));
        vars.AddMatrixWithBC(*Kh0, nat_mat0);
        
        Real_wp err_max = 0;
        for (int i = 0; i < AllLambda.GetM(); i++)
          {
            y = AllEigenVec(i);
            y.Resize(this->GetOffsetV());
            
            this->ReconstructEigenvector(*Kh0, AllLambda(i), y);
            Real_wp erri = 0, normi = 0;
            for (int j = 0; j < y.GetM(); j++)
              {
                erri = max(erri, abs(y(j) - AllEigenVec(i)(j)));
                normi = max(normi, abs(AllEigenVec(i)(j)));     
                if (erri > 1e-4)
                  {
                    DISP(j); DISP(erri); DISP(y(j)); DISP(AllEigenVec(i)(j));
                    abort();
                  }
              }

            err_max = max(err_max, erri / normi);
          }

        cout << "Error on reconstruction = " << err_max << endl;
        delete Kh0;
        }*/
 
    return;
    
  }
  
  void DisplayEigenmodes()
  {
    if (!input_var.display_eigenvectors)
      return;
    
    // on ecrit les valeurs et vecteurs propres
    for (int i = 0; i < AllLambda.GetM(); i++)
      WriteSolution(AllEigenVec(i), input_var.name_output_eigenmode_file + NumberToString(i));
  }

  
  void WriteSolution(const VectComplex_wp& u, const string& root)
  {
    Vector<VectComplex_wp> eigen_mode_v;
    eigen_mode_v.SetData(1, &const_cast<VectComplex_wp&>(u));
    
    // on change les noms de fichier
    for (int j = 0; j < vars.output_mesh_param.GetM(); j++)
      vars.output_mesh_param(j).
	SetFileName(2, root + "G" + to_str(j) + string(".bb"));
    
    for (int j = 0; j < vars.output_grid_param.GetM(); j++)
      vars.output_grid_param(j).
	SetFileName(2, root + "G" + to_str(j) + string(".dat"));
    
    vars.WriteOutputFile(eigen_mode_v, 2);

    // pour ecrire P
    /*if (Dimension::dim_N == 2)
      {
        int nPole = vars.ref_drude(ref_cavity).gamma.GetM(); DISP(nPole);
        Vector<VectComplex_wp> P(nPole);
        for (int k = 0; k < nPole; k++)
          P(k).Reallocate(vars.mesh_num.GetNbDof());
        
        int offset_P = vars.mesh_num.GetNbDof() + vars.mesh_num.GetNbDofPML();
        for (int i = 0; i < vars.mesh_num.GetNbDof(); i++)
          {
            int n2 = vars.GetDrudeDof(i);
            if (n2 >= 0)
              {
                for (int k = 0; k < nPole; k++)
                  P(k)(i) = u(offset_P + n2 + k);
              }
            else
              {
                for (int k = 0; k < nPole; k++)
                  P(k)(i) = 0.0;
              }
          }

        for (int j = 0; j < vars.output_grid_param.GetM(); j++)
          vars.output_grid_param(j).
            SetFileName(2, root + "PolarG" + to_str(j) + string(".dat"));
        
        vars.WriteOutputFile(P, 2);
      }
    */
    eigen_mode_v.Nullify();
  }
  

  void CheckOrthogonality()
  {
    int nev = AllLambda.GetM();
    //Matrix<Complex_wp> mass(nev, nev);
    Real_wp err_extra, err_diag;
    DistributedVector<Complex_wp> biortho(vars.GetOverlappedDofNumber(), vars.comm_group_mode);
    DistributedVector<Complex_wp> y(vars.GetOverlappedDofNumber(), vars.comm_group_mode);
    for (int i = 0; i < nev; i++)
      {
        ComputeBiorthogonal(AllLambda(i), AllEigenVec(i), biortho);
        
	for (int j = 0; j < nev; j++)
	  {
            MltMass(AllEigenVec(j), y);
            if (!vars.FirstOrderFormulation())
              {
                Mlt(Real_wp(2), y);
                MltAddDamping(Real_wp(-2) / (AllLambda(i) + AllLambda(j)), AllEigenVec(j), y);
              }

            Complex_wp scal = DotProd(y, biortho);
            
	    if (j == i)
	      { DISP(i); DISP(scal); DISP(Norm2(biortho));
                err_diag = max(err_diag, abs(Real_wp(1)-scal));   }

	    if (j != i)
              {
                err_extra = max(err_extra, abs(scal));
                if (abs(scal) > 2e-4)
                  {
                    cout << "Modes not orthogonal " << endl;
                    DISP(i); DISP(j); DISP(scal); DISP(AllLambda(i)); DISP(AllLambda(j));
                    //abort();
                  }
              }
	  }
      }

    cout << "Error on diagonal = " << err_diag << endl;
    cout << "Error on non-diagonal elements = " << err_extra << endl;
    //mass.Write("mass.dat");
  }

  void ComputeSourceOmega(const Real_wp& omega, Vector<VectComplex_wp>& Fcplx);

  void ComputeProjectionQNM(const Vector<VectComplex_wp>& vecF, const VectComplex_wp& F,
			    const Real_wp& omega, VectComplex_wp& decomp)
  {
    VectComplex_wp Fcplx(F);
    if (AllEigenVec.GetM() <= 0)
      return;
    
    int N = AllEigenVec(0).GetM();
    Fcplx.Resize(N);
    
    VectComplex_wp y;
    if (input_var.type_projection == input_var.PROJ_BIORTHO)
      {
	decomp.Reallocate(AllLambda.GetM()*2); decomp.Zero();
	for (int i = 0; i < AllLambda.GetM(); i++)
	  {
            if (input_var.dispersive_pml)
              ComputeBiorthogonal(AllLambda(i), AllEigenVec(i), y, input_var.source_wei, false);           
            else
              y = AllEigenVec(i);
            
            if (input_var.use_p_unknown >= 0)
              {
                for (int k = 0; k < vecF.GetM(); k++)
                  {
                    Complex_wp vloc = DotProd(y, vecF(k));
                    Complex_wp alpha = vloc / (-Iwp*omega + AllLambda(i));
                    alpha *= this->GetRatioDeltaEpsilon(omega, -Iwp*AllLambda(i), k);
                    if (input_var.use_p_unknown == 1)
                      alpha *= (-Iwp*AllLambda(i)) / omega;
                    
                    decomp(2*i) += alpha;
                  }
              }
            else
              {
                Complex_wp vloc = DotProd(y, Fcplx);
                decomp(2*i) = vloc / (-Iwp*omega + AllLambda(i));
              }
            
	    if (eigen_vec_complex(i) && (input_var.dispersive_pml || input_var.drop_negative_modes))
              {
                if (input_var.use_p_unknown >= 0)
                  {
                    for (int k = 0; k < vecF.GetM(); k++)
                      {                       
                        Complex_wp vloc = DotProdConj(y, vecF(k));              
                        Complex_wp alpha = vloc / (-Iwp*omega + conjugate(AllLambda(i)));
                        alpha *= this->GetRatioDeltaEpsilon(omega, -Iwp*conjugate(AllLambda(i)), k);
                        if (input_var.use_p_unknown == 1)
                          alpha *= (-Iwp*conjugate(AllLambda(i))) / omega;
                        
                        decomp(2*i+1) += alpha;
                      }
                  }
                else
                  {                     
                    Complex_wp vloc = DotProdConj(y, Fcplx);              
                    decomp(2*i+1) = vloc / (-Iwp*omega + conjugate(AllLambda(i)));                    
                  }
              }
	  }
      }
    else if (input_var.type_projection == input_var.PROJ_MARSEILLE)
      {
        decomp.Reallocate(AllLambda.GetM()*2); decomp.Zero();
	for (int i = 0; i < AllLambda.GetM(); i++)
	  {
	    Complex_wp vloc = DotProd(AllEigenVec(i), Fcplx);
	    decomp(2*i) = vloc / (-Iwp*omega + AllLambda(i)) * (Iwp*omega / AllLambda(i));
            if (input_var.use_p_unknown == 0)
              decomp(2*i) *= this->GetRatioDeltaEpsilon(omega, -Iwp*AllLambda(i));
            else if (input_var.use_p_unknown == 1)
              decomp(2*i) = vloc / (-Iwp*omega + AllLambda(i)) * AllLambda(i) / (Iwp*omega);
            
	    if (eigen_vec_complex(i) && (input_var.dispersive_pml || input_var.drop_negative_modes))
              {
                vloc = DotProdConj(AllEigenVec(i), Fcplx);
                decomp(2*i+1) = vloc / (-Iwp*omega + conjugate(AllLambda(i)))
                  * (Iwp*omega / conjugate(AllLambda(i)));

                if (input_var.use_p_unknown == 0)
                  decomp(2*i+1) *= this->GetRatioDeltaEpsilon(omega, -Iwp*conjugate(AllLambda(i)));
                else if (input_var.use_p_unknown == 1)
                  decomp(2*i+1) = vloc / (-Iwp*omega + conjugate(AllLambda(i)))
                    * conjugate(AllLambda(i)) / (Iwp*omega);
              }
	  }
      }
    else if (input_var.type_projection == input_var.PROJ_MULJAROV)
      {
        decomp.Reallocate(AllLambda.GetM()*2);
	for (int i = 0; i < AllLambda.GetM(); i++)
	  {
	    Complex_wp vloc = DotProd(AllEigenVec(i), Fcplx);
            Complex_wp Om = input_var.pole_muljarov;
	    decomp(2*i) = vloc / (-Iwp*omega + AllLambda(i)) * (omega - Om) / (-Iwp*AllLambda(i) - Om);
            
	    if (eigen_vec_complex(i) && (input_var.dispersive_pml || input_var.drop_negative_modes))
              {
                vloc = DotProdConj(AllEigenVec(i), Fcplx);
                decomp(2*i+1) = vloc / (-Iwp*omega + conjugate(AllLambda(i)))
                  * (omega - Om) / (-Iwp*conjugate(AllLambda(i)) - Om);
              }
	  }
      }
    else
      {
	cout << "Projection not implemented" << endl;
      }

    // on somme les contributions (cas d'un seul resonateur)
    VectComplex_wp dec(decomp);
    MPI_Allreduce(dec.GetData(), decomp.GetData(), 2*dec.GetM(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  void ComputeLinearCombination_Relevement(const VectComplex_wp& decomp,
					   const Real_wp& omega, VectComplex_wp& X)
  {
    int nev = AllLambda.GetM();
    X.Reallocate(vars.GetNbDof());
    X.Zero();
    for (int i = 0; i < nev; i++)
      {
        if (eigen_vec_complex(i) && (input_var.dispersive_pml || input_var.drop_negative_modes))
          for (int j = 0; j < AllEigenVec(i).GetM(); j++)
            X(j) += decomp(2*i)*AllEigenVec(i)(j) + decomp(2*i+1)*conjugate(AllEigenVec(i)(j));
        else
          for (int j = 0; j < AllEigenVec(i).GetM(); j++)
            X(j) += decomp(2*i)*AllEigenVec(i)(j);
      }
  }
  
  void ConstructExcludeEigenvalPML()
  {
    if (!vars.FirstOrderFormulation())
      return;
    
    TinyVector<Real_wp, 1>::threshold = 1e-12;
    Vector<TinyVector<Real_wp, 1> > sigma;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    num_Es.Reallocate(mesh_num.GetNbDof());
    num_Es.Fill();
    for (int i = 0; i < vars.mesh.GetNbElt(); i++)
      {
        const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
        if (vars.InsidePML(i))
          for (int j = 0; j < Fb.GetNbDof(); j++)
            {
              int num_dof = mesh_num.Element(i).GetNumberDof(j);
              int npml = mesh_num.GetDofPML(num_dof);
              if (npml >= 0)
                num_Es(num_dof) = mesh_num.GetNbDof() + npml;
            }
        
	int N   = vars.GetNbPointsQuadratureInside(i);
	int i1 = i - vars.mesh.GetNbElt() + vars.GetNbEltPML();
        if (vars.InsidePML(i))
          {
            Vector<TinyVector<Real_wp, 1> > sig(3*N);
            for (int j = 0; j < N ; j++)
              {
                TinyVector<Complexe, Dimension::dim_N> tau = vars.GetTauPML(i1, j);
                sig(3*j) = TinyVector<Real_wp, 1>(realpart(tau(0)));
                sig(3*j+1) = TinyVector<Real_wp, 1>(realpart(tau(1)));
                if (Dimension::dim_N == 3)
                  sig(3*j+2) = TinyVector<Real_wp, 1>(realpart(tau(2)));
              }

            RemoveDuplicate(sig);
            sigma.PushBack(sig);
            if (i%10 == 0)
              RemoveDuplicate(sigma);
          }
      }
    
    //num_Es.WriteText("num_Es.dat");
    
    RemoveDuplicate(sigma);

    VectReal_wp sigma_(sigma.GetM());
    for (int i = 0; i < sigma_.GetM(); i++)
      sigma_(i) = sigma(i)(0);

    if (input_var.rank_proc == 0)
      cout << "Values taken by damping sigma" << sigma_ << endl;
    
    input_var.exclude_eigenval = sigma_;
  }

  void GetPhysicalDofs(Vector<bool>& is_dof_physical)
  {
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    is_dof_physical.Reallocate(mesh_num.GetNbDof());
    is_dof_physical.Fill(false);
    for (int i = 0; i < vars.mesh.GetNbElt(); i++)
      {
        int ref = vars.mesh.Element(i).GetReference();
	const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
	int N   = Fb.GetNbDof();
	//bool affine = vars.mesh.IsElementAffine(i);
        if (!vars.InsidePML(i))
          if ((input_var.reference_error==-1) || (ref == input_var.reference_error))
            for (int j = 0; j < N ; j++)
              {
                int n = mesh_num.Element(i).GetNumberDof(j);
                
                is_dof_physical(n) = true;
              }
      }
  }


  Real_wp ComputeError(const VectComplex_wp& X, const VectComplex_wp& Y, Vector<bool>& is_dof_physical)
  {
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    Real_wp err = 0, norm = 0;
    if (Dimension::dim_N == 2)
      {
        for (int i = 0; i < mesh_num.GetNbDof(); i++)
          if (is_dof_physical(i))
            {
              err += mass_Mh(i)*absSquare(X(i) - Y(i));
              norm += mass_Mh(i)*absSquare(Y(i));
            }
      }
    else
      {
        int N = mesh_num.GetNbDof();
        VectComplex_wp Mxy(N), My(N), diff(N), Ysub(N);
        Mxy.Zero(); My.Zero();
        for (int i = 0; i < N; i++)
          diff(i) = X(i) - Y(i);

        for (int i = 0; i < N; i++)
          Ysub(i) = Y(i);
        
        Mlt(sparse_Mh, diff, Mxy);
        Mlt(sparse_Mh, Y, My);

        err = abs(DotProdConj(Mxy, diff));
        norm = abs(DotProdConj(My, Ysub));
      }
    
    VectReal_wp vec0(2), vec(2);
    vec0(0) = err; vec0(1) = norm;
    MPI_Allreduce(vec0.GetData(), vec.GetData(), 2, MPI_DOUBLE, MPI_SUM, vars.comm_group_mode);
    err = vec(0); norm = vec(1);

    return sqrt(err / norm);
  }

  Real_wp ComputeErrorGrad(const Vector<Vector<VectComplex_wp> >& U)
  {
    
    Complex_wp dx1, dx2, dx3, dy1, dy2, dy3;
    Real_wp err = 0, norm = 0;
    for (int i = 0; i < vars.mesh.GetNbElt(); i++)
      {
	const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
	int N = Fb.GetNbPointsQuadratureInside();
	bool affine = vars.mesh.IsElementAffine(i);
        int ref = vars.mesh.Element(i).GetReference();
        if (!vars.InsidePML(i))
          if ((input_var.reference_error == -1) || (input_var.reference_error == ref))
            for (int j = 0; j < N ; j++)
              {
                Real_wp jacob = vars.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
                
                if (Dimension::dim_N == 3)
                  {
                    dx1 = U(i)(0)(j);
                    dx2 = U(i)(1)(j);
                    dx3 = U(i)(2)(j);
                    
                    dy1 = U(i)(3)(j);
                    dy2 = U(i)(4)(j);
                    dy3 = U(i)(5)(j);                  
                    
                    err += jacob*(absSquare(dx1-dy1) + absSquare(dx2-dy2) + absSquare(dx3-dy3));
                    norm  += jacob*(absSquare(dy1) + absSquare(dy2) + absSquare(dy3));
                  }
                else
                  {
                    dx1 = U(i)(0)(2*j);
                    dx2 = U(i)(0)(2*j+1);
                    dy1 = U(i)(1)(2*j);
                    dy2 = U(i)(1)(2*j+1);
                    
                    err += jacob*(absSquare(dx1-dy1) + absSquare(dx2-dy2));
                    norm  += jacob*(absSquare(dy1) + absSquare(dy2));
                  }
              }
      }

    VectReal_wp vec0(2), vec(2);
    vec0(0) = err; vec0(1) = norm;
    MPI_Allreduce(vec0.GetData(), vec.GetData(), 2, MPI_DOUBLE, MPI_SUM, vars.comm_group_mode);
    err = vec(0); norm = vec(1);

    return sqrt(err / norm);
  }


  void FindEigenvectorNumber(const VectComplex_wp& Lambda1, const VectComplex_wp& Lambda2,
                             IVect& num_mode, Vector<VectComplex_wp>& X1, 
                             Vector<VectComplex_wp>& X2)
  {
    VectR2 old_eigenval(Lambda1.GetM());
    for (int i = 0; i < Lambda1.GetM(); i++)
      old_eigenval(i).Init(realpart(Lambda1(i)), imagpart(Lambda1(i)));
    
    VectR2 new_eigenval(Lambda2.GetM());
    for (int i = 0; i < Lambda2.GetM(); i++)
      new_eigenval(i).Init(realpart(Lambda2(i)), imagpart(Lambda2(i)));
    
    IVect permut_old(Lambda1.GetM()); permut_old.Fill();
    Sort(old_eigenval, permut_old);
    
    IVect permut_new(Lambda2.GetM()); permut_new.Fill();
    Sort(new_eigenval, permut_new);
    
    num_mode.Reallocate(Lambda2.GetM());
    num_mode.Fill(-1);
    int k = 0; VectComplex_wp x, y;
    for (int i = 0; i < old_eigenval.GetM(); i++)
      {
        if ((i > 0) && (old_eigenval(i) == old_eigenval(i-1)))
          continue;
        
        // we search multiple eigenvalues for lambda1
        int n = 0;
        while ((i+n < Lambda1.GetM()) && (old_eigenval(i+n) == old_eigenval(i)))
          n++;
        
        while ((k < Lambda2.GetM()) && (new_eigenval(k) < old_eigenval(i)))
          k++;
        
        // multiple eigenvalues for lambda2
        int m = 0;
        while ((k+m < Lambda2.GetM()) && (new_eigenval(k+m) == old_eigenval(i)))
          {
            // we found the corresponding eigenvalue
            num_mode(permut_new(k+m)) = permut_old(i);
            m++;
          }
        
        // treating the case where the eigenvalues coincide
        if (m > 0)
          {
            // on cree une matrice de masse associee a la base des vecteurs propres
            Matrix<Complex_wp, Hermitian, RowHermPacked> mass(m+n, m+n);
            VectReal_wp lambda_mass;
            Matrix<Complex_wp> eigenvec_mass;
            for (int i0 = 0; i0 < m+n; i0++)
              for (int j0 = i0; j0 < m+n; j0++)
                {
                  if (i0 < n)
                    x = X1(permut_old(i+i0));
                  else
                    x = X2(permut_new(k+i0-n));
                  
                  if (j0 < n)
                    y = X1(permut_old(i+j0));
                  else
                    y = X2(permut_new(k+j0-n));
                  
                  mass.Set(i0, j0, DotProdConj(x, y));
                }
            
            GetEigenvaluesEigenvectors(mass, lambda_mass, eigenvec_mass);
            
            int dim_kernel = 0;
            for (int j = 0; j < lambda_mass.GetM(); j++)
              if (abs(lambda_mass(j)) <= 1e-12)
                dim_kernel++;
            
            int rank_mass = m+n - dim_kernel;
            if (dim_kernel == 0)
              {
                // on garde tous les vecteurs propres
                for (int j = 0; j < m; j++)
                  num_mode(permut_new(k+j)) = -1;
              }
            else if (rank_mass > n)
              {
                // on garde les vecteurs propres qui sont independants
                int num = 0;
                for (int j = 0; j < lambda_mass.GetM(); j++)
                  if (abs(lambda_mass(j)) > 1e-12)
                    {
                      x.Zero();
                      for (int i0 = 0; i0 < m+n; i0++)
                        {
                          if (i0 < n)
                            Add(eigenvec_mass(i0, j), X1(permut_old(i+i0)), x);
                          else
                            Add(eigenvec_mass(i0, j), X2(permut_new(k+i0-n)), x);
                        }
                      
                      
                      if (j < n)
                        X1(permut_old(i+j)) = x;
                      else
                        {
                          X2(permut_new(k+j-n)) = x;
                          num_mode(permut_new(k+j-n)) = -1;
                        }
                      
                      num++;
                    }
                
                for (int j = num; j < n+m; j++)
                  {
                    // on met num_mode a zero pour les vecteurs propres restants
                    num_mode(permut_new(k+j)) = 0;
                  }
              }
          }
        
        k += m;
      }
  }
  
  
  void MergeEigenvectors(VectComplex_wp& Lambda1, VectComplex_wp& Lambda2,
                         Vector<bool>& eigen_vec_cplx1, Vector<bool>& eigen_vec_cplx2, 
                         Vector<VectComplex_wp>& X1, Vector<VectComplex_wp>& X2)
  {
    IVect num_mode;
    FindEigenvectorNumber(Lambda1, Lambda2, num_mode, X1, X2);

    int nb_new = 0;
    for (int i = 0; i < num_mode.GetM(); i++)
      if (num_mode(i) == -1)
        nb_new++;
    
    if (nb_new > 0)
      {
        int nev = Lambda1.GetM();
        Lambda1.Resize(nev + nb_new);
        eigen_vec_cplx1.Resize(nev + nb_new);
        X1.Resize(nev + nb_new);
        int num = nev;
        for (int i = 0; i < num_mode.GetM(); i++)
          if (num_mode(i) == -1)
            {
              Lambda1(num) = Lambda2(i);
              eigen_vec_cplx1(num) = eigen_vec_cplx2(i);
              X1(num) = X2(i);
              num++;
            }
      }    
  }

  Real_wp ComputeDisplayFieldE(int m, const Real_wp& omega, const VectComplex_wp& Y,
                               const Vector<Vector<VectComplex_wp> >& GradQuad, VectComplex_wp& rhs);


  void ComputeSolutionOnInterpolationPoints(All_LinearSolver* glob_solver)
  {      
    if (input_var.nb_interval_interpolation <= 0)
      return;

    VectReal_wp leja_points; 
    string points("0.5	0	1	0.211324865405187	0.829353297207782	0.0803729132191221	0.935003574854083	0.652806664558611	0.339146193942521	0.0285104089150469	0.976336635615583	0.739706164461324	0.143680679821208	0.422020317760199	0.887436170755216	0.0102611906570433	0.580582634266648	0.991663154768928	0.269314698789435	0.0540535895403992	0.785948542059528");
    istringstream input(points);
    leja_points.ReadText(input);

    // premiere etape : on calcule les points d'interpolation
    int Ni = input_var.order_interpolation;
    interp.Reallocate(input_var.nb_levels_interpolation);
    int Ni_next = 2*((Ni+1)/2)+1;
    if (input_var.type_interpolation == input_var.TCHEBY)
      {
        Ni_next = 1;
        while (Ni_next <= Ni)
          Ni_next *= 3;
      }

    if (input_var.type_interpolation == input_var.CLENSHAW)
      {
        Ni_next = 2;
        while (Ni_next+1 <= Ni)
          Ni_next *= 2;
      }
    
    for (int lvl = 0; lvl < input_var.nb_levels_interpolation; lvl++)
      {
        VectReal_wp pts(Ni);
        if (input_var.type_interpolation == input_var.TCHEBY)
          {
            for (int i = 0; i < Ni; i++)
              pts(i) = 0.5 + 0.5*cos((2.0*i+1.0) / (2.0*Ni) * pi_wp);
                                                
            Ni = Ni_next; Ni_next *= 3;
          }
        else if (input_var.type_interpolation == input_var.CLENSHAW)
          {
            if (Ni == 1)
              pts(0) = 0.5;
            else
              for (int i = 0; i < Ni; i++)
                pts(i) = 0.5 + 0.5*cos(Real_wp(i) / (Ni-1) * pi_wp);
            
            Ni = Ni_next + 1; Ni_next *= 2;
          }
        else if (input_var.type_interpolation == input_var.LEJA)
          {
            for (int k = 0; k < Ni; k++)
              pts(k) = leja_points(k);

            Ni = Ni_next; Ni_next += 2;
          }

        if (input_var.rank_proc == 0)
          {
            DISP(lvl); DISP(pts);
          }
        
        Sort(pts);
        interp(lvl).InitPoints(1, pts); 
      }

    // on calcule les numeros des niveaux intermediaires
    int lvl_max = input_var.nb_levels_interpolation-1;
    VectReal_wp pts_interp = interp(lvl_max).Points();
    if (input_var.rank_proc == 0)
      DISP(pts_interp);
    num_interp.Reallocate(lvl_max+1);
    num_interp(lvl_max).Reallocate(pts_interp.GetM());
    num_interp(lvl_max).Fill();
    for (int lvl = 0; lvl < lvl_max; lvl++)
      {
        VectReal_wp pts = interp(lvl).Points();
        int num = 0;
        num_interp(lvl).Reallocate(pts.GetM());
        num_interp(lvl).Fill(-1);
        for (int k = 0; k < pts.GetM(); k++)
          {
            while ((num < pts_interp.GetM()) && (pts_interp(num) < pts(k)))
              num++;

            if ((num < pts_interp.GetM()) && (abs(pts_interp(num)-pts(k)) < 1e-12))
              num_interp(lvl)(k) = num;
            else
              {
                cout << "Non-nested sequence of points" << endl;
                DISP(lvl); DISP(pts);
                abort();
              }
          }
      }
    
    // on calcule la solution sur les points d'interpolation du dernier niveau
    // boucle sur les points d'interpolation
    int N = vars.GetNbDof();
    VectComplex_wp Fcplx(N), Xfem(N), decomp, Xqnm(N);
    Vector<VectComplex_wp> vecF;
    diff_points_interpolation.Reallocate(pts_interp.GetM());
    for (int i = 0; i < pts_interp.GetM(); i++)
      {
        // pulsation associee au point i
        Real_wp f = input_var.fmin + pts_interp(i)*(input_var.fmax - input_var.fmin);
        Real_wp omega = 2.0*pi_wp*f;
        
        // calcul de la source
        this->ComputeSourceOmega(omega, vecF);
        Fcplx = vecF(0);
        for (int k = 1; k < vecF.GetM(); k++)
          Fcplx += vecF(k);
        
        // calcul de la solution FEM
        GlobalGenericMatrix<Complex_wp> nat_mat;
        if (input_var.dispersive_pml)
          {
            nat_mat.SetCoefMass(-Iwp*omega);
            nat_mat.SetCoefStiffness(Complex_wp(1, 0));
            nat_mat.SetCoefDamping(Complex_wp(1, 0));
          }
        else
          {
            Real_wp coef = omega*omega / vars.GetSquareOmega();
            nat_mat.SetCoefMass(coef / (-Iwp*omega));
            coef = omega / vars.GetOmega();
            nat_mat.SetCoefDamping(coef / (-Iwp*omega)); 
            coef = Real_wp(1);
            nat_mat.SetCoefStiffness(coef / (-Iwp*omega));
          }

        //glob_solver->EnableMatrixStorage(true, "mat_test.dat");
        glob_solver->EnableStaticCondensation(true, true, true);
        glob_solver->PerformFactorizationStep(nat_mat, true);
        Xfem = Fcplx; glob_solver->ComputeSolution(Xfem);
        
        // calcul de la solution QNM
        //if (!vars.GetSymmetrizationUse())
        //vars.ModifySourceSymmetry(Fcplx);
        
        // on projette sur les QNMs (calcul des alpha_i dans decomp)
	ComputeProjectionQNM(vecF, Fcplx, omega, decomp);
        
	// on calcule la combinaison lineaire \sum \alpha_i qnm_i
	ComputeLinearCombination_Relevement(decomp, omega, Xqnm);
                       
        // on stocke la difference Xfem - Xqnm
        diff_points_interpolation(i) = Xfem - Xqnm;
      }
  }

  void ReconstructInterpolationQNM(int lvl, const Real_wp& omega, VectComplex_wp& X)
  {
    if (input_var.nb_interval_interpolation <= 0)
      return;
    
    // on calcule le x_loc a partir de omega
    Real_wp f = omega / (2.0*pi_wp);
    Real_wp x_loc = (f - input_var.fmin) / (input_var.fmax - input_var.fmin);
    //DISP(x_loc);

    // on rajoute l'interpolation de la difference
    for (int k = 0; k <= interp(lvl).GetOrder(); k++)
      {
        Real_wp phi = interp(lvl).EvaluatePhi(k, x_loc); //DISP(k); DISP(phi);
        Add(phi, diff_points_interpolation(num_interp(lvl)(k)), X);
      }
  }

  void RunAll(const Vector<string>& lines_data_file, const string& name_element,
	      const string& name_equation)
  {    
    glob_chrono.SetMessage("AllSimulation", "complete the simulation");
    glob_chrono.Start("AllSimulation");

    // calcul des matrices Ah, Mh et Kh
    this->ConstructProblem(lines_data_file, name_element, name_equation);
    
    ConstructExcludeEigenvalPML();

    MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    if (input_var.rank_proc == 0)
      {
        cout << "Number of dofs for E = " << mesh_num.GetNbDof() << endl;
        cout << "Size of eigenvalue problem = " << vars.GetNbDof() << endl;
      }
    
    if (input_var.merge_eigenvalues)
      {
        this->LoadEigenvectors(input_var.file_name_eigenvectors(0),
                               AllLambda, eigen_vec_complex, AllEigenVec);
        
        int N = input_var.file_name_eigenvectors.GetM();
        if (N < 3)
          {
            cout << "Provide at least three file names"  << endl;
            abort();
          }
        
        VectComplex_wp LambdaNew;
        Vector<bool> eigen_vec_cplx;
        Vector<VectComplex_wp> EigenVecNew;
        for (int j = 1; j < N-1; j++)
          {
            this->LoadEigenvectors(input_var.file_name_eigenvectors(j),
                                   LambdaNew, eigen_vec_cplx, EigenVecNew);
            
            MergeEigenvectors(AllLambda, LambdaNew, eigen_vec_complex, eigen_vec_cplx,
                              AllEigenVec, EigenVecNew);
          }
        
        this->SaveEigenvectors(input_var.file_name_eigenvectors(N-1));
        return;
      }
    
    if (input_var.load_eigen_vec)
      {
        cout << "Loading eigenvectors..." << endl;
        
        this->LoadEigenvectors(input_var.file_name_eigenvectors(0),
                               AllLambda, eigen_vec_complex, AllEigenVec);
        
        cout << "Reconstructing eigenvectors..." << endl;

        //if (input_var.source_wei)
        this->ReconstructEigenvectors();
        
        if (input_var.extractQNM)
          ExtractQnmModes();
        else
          AllLambda.Write("Lambda.dat");
      }
    else
      {
        // calcul des valeurs propres qu'on conserve (modes QNM)
        if (input_var.rank_proc == 0)
          cout << "On calcule les modes" << endl;
        
        glob_chrono.SetMessage("EigenModes", "compute eigenmodes");
        glob_chrono.Start("EigenModes");
        this->ComputeEigenmodes();

        glob_chrono.Stop("EigenModes");
        glob_chrono.DisplayTime("EigenModes");
        
        // on affiche les modes
        if (input_var.rank_proc == 0)
          cout << "On affiche les modes" << endl;
        
        this->DisplayEigenmodes();
        
        // on verifie la relation d'orthogonalite
        if (input_var.rank_proc == 0)
          cout << "On verifie l'orthogonalite" << endl;
        
        if (input_var.check_orthogonality)
          this->CheckOrthogonality();
        
        if (input_var.save_eigen_vec)
          this->SaveEigenvectors(input_var.file_name_eigenvectors(0));
      }
    
    if (input_var.nb_freq <= 0)
      return;
    
    int N = vars.GetNbDof();
    VectComplex_wp Fcplx(N), X(N), decomp;
    
    // boucle sur les frequences
    if (input_var.rank_proc == 0)
      cout << "Boucle en frequence" << endl;
    
    bool tm_mode = true;
    for (int ref = 0; ref < vars.GetNbPhysicalIndices(); ref++)
      if (vars.ref_drude(ref).IsEnabled() && vars.ref_drude(ref).IsModeTE())
        tm_mode = false;

    vars.SetCoefficientDirichlet(Real_wp(1));
    All_LinearSolver* glob_solver_fem;
    glob_solver_fem = vars.GetNewLinearSolver();
    glob_solver_fem->SetDirectSolver();

    Vector<bool> is_dof_physical;
    GetPhysicalDofs(is_dof_physical);
    ComputePhysicalMass();
    
    //if ((input_var.type_projection == input_var.PROJ_MULJAROV) && (input_var.rank_proc == 0))
    if (input_var.rank_proc == 0)
      {
        VectComplex_wp all_poles = vars.ref_drude(ref_cavity(0)).GetPoleOmega();
        VectComplex_wp all_zeros = vars.ref_drude(ref_cavity(0)).GetZeroOmega();
        cout << "Displaying roots of epsilon" << endl;
        for (int i = 0; i < all_zeros.GetM(); i++)
          {
            DISP(all_zeros(i));
            DISP(vars.ref_drude(ref_cavity(0)).eps_inf + vars.ref_drude(ref_cavity(0)).GetDeltaEpsilon(all_zeros(i)));
          }

        cout << endl << "Displaying roots of epsilon = -1" << endl;
        VectComplex_wp all_mones = vars.ref_drude(ref_cavity(0)).GetZeroOmega(Real_wp(-1));
        for (int i = 0; i < all_mones.GetM(); i++)
          {
            DISP(all_mones(i));
            DISP(vars.ref_drude(ref_cavity(0)).eps_inf + vars.ref_drude(ref_cavity(0)).GetDeltaEpsilon(all_mones(i)));
          }

        cout << endl;
        DISP(vars.ref_drude(ref_cavity(0)).eps_inf);
        DISP(all_poles);
        if (input_var.type_projection == input_var.PROJ_MULJAROV)
          input_var.pole_muljarov = all_poles(input_var.use_p_unknown);
        DISP(input_var.pole_muljarov);
      }

    // on calcule la solution sur les points d'interpolation
    ComputeSolutionOnInterpolationPoints(glob_solver_fem);
    
    VectComplex_wp Y; VectReal_wp err(input_var.nb_freq);
    Matrix<Real_wp> errI(input_var.nb_freq, input_var.nb_levels_interpolation);
    err.Zero(); errI.Zero();
    VectReal_wp errE(input_var.nb_freq); errE.Zero();
    Vector<VectComplex_wp> vecF;
    for (int m = input_var.first_freq; m < input_var.nb_freq; m++)
      {
	Real_wp f = input_var.fmin + m*(input_var.fmax-input_var.fmin)/(input_var.nb_freq-1);
	Real_wp omega = 2*pi_wp*f;
	if (input_var.rank_proc == 0)
          {
            DISP(omega);
            DISP(vars.ref_drude(ref_cavity(0)).GetDeltaEpsilon(omega));
          }
        
	// on remplit la source Fcplx
	this->ComputeSourceOmega(omega, vecF);
        Fcplx = vecF(0);
        for (int k = 1; k < vecF.GetM(); k++)
          Fcplx += vecF(k);
        //Fcplx.Write("rhs.dat");

        if (input_var.compute_fem_sol)
          {
            if (input_var.load_fem_sol)
              X.Read(input_var.path_fem_sol + "/XsolRef"+to_str(m) + "_P" + to_str(input_var.rank_proc) + ".dat");
            else
              {
                GlobalGenericMatrix<Complex_wp> nat_mat;
                if (input_var.dispersive_pml)
                  {
                    nat_mat.SetCoefMass(-Iwp*omega);
                    nat_mat.SetCoefStiffness(Complex_wp(1, 0));
                    nat_mat.SetCoefDamping(Complex_wp(1, 0));
                  }
                else
                  {
                    Real_wp coef = omega*omega / vars.GetSquareOmega();
                    nat_mat.SetCoefMass(coef / (-Iwp*omega));
                    coef = omega / vars.GetOmega();
                    nat_mat.SetCoefDamping(coef / (-Iwp*omega)); 
                    coef = Real_wp(1);
                    nat_mat.SetCoefStiffness(coef / (-Iwp*omega));
                  }
                
                //glob_solver_fem->EnableMatrixStorage(true, "mat_test.dat");
                glob_solver_fem->EnableStaticCondensation(true, true, true);
                glob_solver_fem->PerformFactorizationStep(nat_mat, true);
                X = Fcplx; glob_solver_fem->ComputeSolution(X);
                if (input_var.write_fem_sol)
                  X.Write(input_var.path_fem_sol + "/XsolRef" + to_str(m) + "_P" + to_str(input_var.rank_proc) + ".dat");
              }
            
            //X.Zero(); X(95) = 1.0;
            WriteSolution(X, input_var.file_fem_sol + to_str(m));
            Y = X;
          }
        else
          X = Fcplx;
        
        //X.Write("sol.dat");
        
        X.Zero();
        VectComplex_wp rhs(Fcplx);
        
        // on projette sur les QNMs (calcul des alpha_i dans decomp)
	ComputeProjectionQNM(vecF, Fcplx, omega, decomp);
	//decomp.Write("alpha"+ to_str(m) + ".dat");	
	
	// on calcule la combinaison lineaire \sum \alpha_i qnm_i
	ComputeLinearCombination_Relevement(decomp, omega, X);
        
	//decomp.Write("beta"+ to_str(m) + ".dat");
        
        // on reconstruit la solution en rajoutant l'interpolation de la difference
        VectComplex_wp Xinit = X;
        for (int k = 0; k < input_var.nb_levels_interpolation; k++)
          {
            X = Xinit;
            ReconstructInterpolationQNM(k, omega, X);
            if ((input_var.use_p_unknown >= 0) && input_var.source_wei)
              this->AddIncidentWave(Real_wp(-1), omega, X);
            
            //X.Write("sol_qnm.dat");
            if (k == input_var.nb_levels_interpolation-1)
              WriteSolution(X, input_var.file_qnm_sol + to_str(m));
            
            if (input_var.compute_fem_sol)
              {
                if (input_var.compute_err_grad)
                  {
                    Vector<VectComplex_wp> U(2);
                    Vector<Vector<VectComplex_wp> > Uquad, GradQuad;
                    U(0) = X; U(1) = Y;
                    vars.ComputeQuadratureUgradU(U, Uquad, GradQuad, false, true, true, 2);
                    
                    if (!tm_mode)
                      errE(m) = this->ComputeDisplayFieldE(m, omega, Y, GradQuad, rhs);
                    
                    err(m) = ComputeErrorGrad(GradQuad);
                  }
                else
                  err(m) = ComputeError(X, Y, is_dof_physical);
                
                if (input_var.rank_proc == 0)
                  {
                    DISP(m); DISP(err(m)); DISP(errE(m));
                    err.WriteText("err" + to_str(input_var.first_freq) + ".dat");
                  }

                errI(m, k) = err(m);
              }
            //cout << "Waiting" << endl; int test_input; cin >> test_input;
          }
      }
    
    if (input_var.rank_proc == 0)
      {
        DISP(AllLambda.GetM());
        DISP(err); DISP(err.GetNormInf());
        DISP(errE); DISP(errE.GetNormInf());
        errI.WriteText("ErreurInterpolation" + to_str(input_var.first_freq) + ".dat");
        err.WriteText("err" + to_str(input_var.first_freq) + ".dat");
      }

    glob_chrono.Stop("AllSimulation");
    glob_chrono.DisplayTime("AllSimulation");
  }

};


template<>
int QnmSolver<StaticMaxwellEquation_3D>::GetOffsetV()
{
  const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
  return mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
}

template<>
int QnmSolver<HarmonicMaxwellEquation_3D>::GetOffsetV()
{
  return vars.GetMeshNumbering(0).GetNbDof();
}

template<>
int QnmSolver<LaplaceEquation<Dimension2> >::GetOffsetV()
{
  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  //return vars.GetNbDof();
  return mesh_num.GetNbDof();
  
  //int offset_P = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
  //int offset_Q = offset_P + vars.GetNbDofDrude();
  //int offset_v = offset_Q + vars.GetNbDofDrude();
  //return offset_v;
}

template<>
int QnmSolver<HelmholtzEquation<Dimension2> >::GetOffsetV()
{
  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  return mesh_num.GetNbDof();
  //return vars.GetNbDof();
}

template<>
void QnmSolver<HelmholtzEquation<Dimension2> >
::ComputeRationalProblem(Vector<DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse> >& vec_Ai,
                         Vector<Vector<Complex_wp> >& coef_num, Vector<Vector<Complex_wp> >& coef_denom)
{
  Vector<int> IndexDrude(vars.GetNbPhysicalIndices());
  IndexDrude.Fill(-1);
  int nb_drude = 0;
  for (int i = 1; i < vars.GetNbPhysicalIndices(); i++)
    if (vars.ref_drude(i).IsEnabled())
      {
        IndexDrude(i) = nb_drude;
        nb_drude += vars.ref_drude(i).gamma.GetM();
      }
  
  vec_Ai.Reallocate(2+nb_drude);
  coef_num.Reallocate(2+nb_drude);
  coef_denom.Reallocate(2+nb_drude);
  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  int N = mesh_num.GetNbDof();
  for (int i = 0; i < vec_Ai.GetM(); i++)
    {
      vec_Ai(i).Clear();
      vec_Ai(i).Reallocate(N, N);
    }

  Vector<Complex_wp> coef(3);
  // premier terme en lambda^2
  coef.Zero(); coef(0) = Complex_wp(1, 0);
  coef_num(0) = coef;
  coef.Clear(); coef_denom(0) = coef;

  // second terme constant
  coef.Reallocate(1); coef(0) = Complex_wp(1, 0);
  coef_num(1) = coef;
  coef.Clear(); coef_denom(1) = coef;

  // termes rationnels pour Drude
  for (int i = 1; i < vars.GetNbPhysicalIndices(); i++)
    if (vars.ref_drude(i).IsEnabled())
      {
        int num = IndexDrude(i);
        for (int k = 0; k < vars.ref_drude(i).gamma.GetM(); k++)
          {
            coef.Reallocate(4); coef.Zero();
            coef(0) = Iwp*vars.GetOmega()*vars.ref_drude(i).eps_sigma(k);
            coef(1) = -vars.ref_drude(i).eps_omega_p2(k);
            coef_num(2 + num + k) = coef;
            
            coef.Reallocate(3);
            coef(0) = vars.GetSquareOmega();
            coef(1) = Iwp*vars.GetOmega()*vars.ref_drude(i).gamma(k);
            coef(2) = -vars.ref_drude(i).omega_02(k);            
            coef_denom(2 + num + k) = coef;
          }
      }
  
  // boucle sur les elements pour calculer les matrices
  Matrix<Complex_wp> mat_elem;
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      if (!vars.UseNumericalIntegration(i))
        {
          cout << "Not implemented" << endl;
          abort();
        }
      
      int ref = vars.mesh.Element(i).GetReference();
      int num = IndexDrude(ref);
      const ElementReference<Dimension2, 1>& Fb = vars.GetReferenceElementH1(i);
      int nb_points_quad = Fb.GetNbPointsQuadratureInside();
      Vector<TinyMatrix<Complex_wp, General, 2, 2> > Cgrad(nb_points_quad);
      Vector<TinyVector<Complex_wp, 2> > Dvec, Evec;
      Vector<Complex_wp> Amass;
      for (int j = 0; j < nb_points_quad; j++)
        {
          int num_dof = mesh_num.Element(i).GetNumberDof(j);
          Complex_wp vloc = vars.Glob_matMass_Dh(i)(j);
          vec_Ai(0).AddInteraction(num_dof, num_dof, vloc);
          if (num >= 0)
            {
              vloc /= vars.ref_drude(ref).eps_inf;
              vec_Ai(2+num).AddInteraction(num_dof, num_dof, vloc);
            }
          
          Cgrad(j) = vars.Glob_matMass_Bh(i)(j);
        }

      TinyVector<bool, 4> null_term(true, false, true, true);
      mat_elem.Reallocate(Fb.GetNbDof(), Fb.GetNbDof());
      mat_elem.Zero();
      Fb.AddVariableElemMatrix(0, 0, Amass, Cgrad, Dvec, Evec, null_term, mat_elem);

      Vector<int> nodle = mesh_num.Element(i).GetNodle();
      for (int j = 0; j < Fb.GetNbDof(); j++)
        for (int k = 0; k < Fb.GetNbDof(); k++)
          vec_Ai(1).AddInteraction(nodle(j), nodle(k), mat_elem(j, k));
    }

  for (int i = 1; i < vars.GetNbPhysicalIndices(); i++)
    if (vars.ref_drude(i).IsEnabled())
      {
        int num = IndexDrude(i);
        for (int k = 1; k < vars.ref_drude(i).gamma.GetM(); k++)
          vec_Ai(2+num+k) = vec_Ai(2+num);
      }
  
  for (int i = 0; i < vec_Ai.GetM(); i++)
    {
      EraseRow(vars.GetDirichletDofNumber(), vec_Ai(i));
      EraseCol(vars.GetDirichletDofNumber(), vec_Ai(i));
    }

  for (int i = 0; i < vars.GetNbDirichletDof(); i++)
    {
      int num_dof = vars.GetDirichletDofNumber(i);
      vec_Ai(0).AddInteraction(num_dof, num_dof, Real_wp(1));
    }
}

template<>
void QnmSolver<StaticMaxwellEquation_3D>
::ReconstructEigenvector(const VirtualMatrix<Real_wp>& mat_stiff,
                         const Complex_wp& L, VectComplex_wp& X)
{
  int nb = this->GetOffsetV();
  X.Resize(vars.GetNbDof());
  for (int i = nb; i < X.GetM(); i++)
    SetComplexZero(X(i));
  
  VectComplex_wp Y(vars.GetNbDof());
  Y.Zero();
  mat_stiff.MltVector(X, Y); 

  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      int ref = vars.mesh.Element(i).GetReference();
      bool affine = vars.mesh.IsElementAffine(i);
      const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
      Real_wp mu = vars.ref_mu(ref).GetConstant()(0, 0);
      Complex_wp mu_cplx(mu, 0);
      if (vars.ref_drude(ref).IsEnabled() && vars.ref_drude(ref).IsModeTE())
        {
          for (int kp = 0; kp < vars.ref_drude(ref).gamma.GetM(); kp++)
            {
              Real_wp gamma = vars.ref_drude(ref).gamma(kp);
              Real_wp coef_omega2 = vars.ref_drude(ref).omega_02(kp);
              Real_wp coef_omegap = vars.ref_drude(ref).eps_omega_p2(kp);
              Real_wp sig = vars.ref_drude(ref).eps_sigma(kp);
              
              mu_cplx -= (coef_omegap - L*sig) / (-L*L - coef_omega2 + L*gamma);
            }
        }
      
      int N = Fb.GetNbPointsQuadratureInside();
      int i1 = i - vars.mesh.GetNbElt() + vars.GetNbEltPML();
      
      // loop on quadrature points
      Real_wp jacobian(1); Complex_wp coef; R3 tau;
      int offset = this->GetOffsetV() + vars.GetOffsetDofV(i);
      int offset_Q = vars.GetNbDof() - vars.GetNbVectorialDofDrude();
      int offset_P = offset_Q - vars.GetNbVectorialDofDrude();
      offset_P += vars.OffsetDofDrudeV(i);
      offset_Q += vars.OffsetDofDrudeV(i);
      
      for (int j = 0; j < N; j++)
        {
          int num_dofHx = offset + 3*j;
          int num_dofHy = num_dofHx + 1;
          int num_dofHz = num_dofHy + 1;
          
          if (affine)
            jacobian = vars.Glob_jacobian(i)(0)*Fb.WeightsND(j);
          else
            jacobian = vars.Glob_jacobian(i)(j);
          
          if (vars.InsidePML(i))
            {
              int num_dofHx_s = offset + 3*N + 3*j;
              int num_dofHy_s = num_dofHx_s + 1, num_dofHz_s = num_dofHx_s+2;
              
              tau = vars.GetTauPML(i1, j);
              X(num_dofHx) = Y(num_dofHx) / (jacobian*mu*(L - tau(1)));
              X(num_dofHy) = Y(num_dofHy) / (jacobian*mu*(L - tau(2)));
              X(num_dofHz) = Y(num_dofHz) / (jacobian*mu*(L - tau(0)));

              X(num_dofHx_s) = X(num_dofHx)*(L-tau(0)) / (L-tau(2));
              X(num_dofHy_s) = X(num_dofHy)*(L-tau(1)) / (L-tau(0));
              X(num_dofHz_s) = X(num_dofHz)*(L-tau(2)) / (L-tau(1));
            }
          else
            {
              coef = Real_wp(1) / (jacobian*mu_cplx*L);
              
              X(num_dofHx) = Y(num_dofHx)*coef;
              X(num_dofHy) = Y(num_dofHy)*coef;
              X(num_dofHz) = Y(num_dofHz)*coef;
            }
          
          if (vars.ref_drude(ref).IsEnabled())
            {
              for (int kp = 0; kp < vars.ref_drude(ref).gamma.GetM(); kp++)
                {
                  int num_dofPx = offset_P + 3*(N*kp + j);
                  int num_dofQx = offset_Q + 3*(N*kp + j);
                  
                  Real_wp gamma = vars.ref_drude(ref).gamma(kp);
                  Real_wp coef_omega2 = vars.ref_drude(ref).omega_02(kp);
                  Real_wp coef_omegap = vars.ref_drude(ref).eps_omega_p2(kp);
                  Real_wp sig = vars.ref_drude(ref).eps_sigma(kp);
                  
                  coef = (coef_omegap-L*sig) / (-L*L - coef_omega2 + L*gamma);
                  if (!vars.ref_drude(ref).IsModeTE())
                    {
                      coef /= (jacobian*coef_omegap);
                      
                      X(num_dofPx) = coef*Y(num_dofQx);
                      X(num_dofPx+1) = coef*Y(num_dofQx+1);
                      X(num_dofPx+2) = coef*Y(num_dofQx+2);
                    }
                  else
                    {
                      X(num_dofPx) = -coef*X(num_dofHx);
                      X(num_dofPx+1) = -coef*X(num_dofHy);
                      X(num_dofPx+2) = -coef*X(num_dofHz);
                    }
                  
                  X(num_dofQx) = -L*X(num_dofPx);
                  X(num_dofQx+1) = -L*X(num_dofPx+1);
                  X(num_dofQx+2) = -L*X(num_dofPx+2);
                }
            }
        }
    }
}

template<>
void QnmSolver<HarmonicMaxwellEquation_3D>
::ReconstructEigenvector(const VirtualMatrix<Real_wp>& mat_stiff,
                         const Complex_wp& L, VectComplex_wp& X)
{
  return;
}

template<>
void QnmSolver<LaplaceEquation<Dimension2> >
::ReconstructEigenvector(const VirtualMatrix<Real_wp>& mat_stiff,
                         const Complex_wp& L, VectComplex_wp& x)
{
  return;
  
  cout << "Not implemented" << endl;
  abort();
}

template<>
void QnmSolver<HelmholtzEquation<Dimension2> >
::ReconstructEigenvector(const VirtualMatrix<Real_wp>& mat_stiff,
                         const Complex_wp& L, VectComplex_wp& x)
{
  return;
}

template<>
void QnmSolver<StaticMaxwellEquation_3D>::ComputeSourceOmega(const Real_wp& omega, Vector<VectComplex_wp>& vecF)
{
  vecF.Reallocate(ref_cavity.GetM());
  for (int k = 0; k < vecF.GetM(); k++)
    {
      vecF(k).Reallocate(vars.GetNbDof());
      vecF(k).Zero();
    }
  
  // on calcule la source
  bool incident_source = false;
  R3 k, polar;
  GaussianSourceField<Real_wp, Dimension> gaussian;
  if (vars.GetSourceType(0) == vars.SRC_DIFFRACTED_FIELD)
    {
      incident_source = true;
      k = vars.GetWaveVector();
      Mlt(1.0/Norm2(k), k);
    }
  else if (vars.GetSourceType(0) == vars.SRC_VOLUME)
    {
      VectString param = vars.GetParameterSource(0, 0);
      int nb = -1;
      for (int k = 0; k < param.GetM(); k++)
        if (param(k) == "GAUSSIAN")
          {
            nb = k+1;
            vars.InitGaussianParameter(gaussian, param, nb);
          }
        
      if (nb == -1)
        {
          cout << "Case not implemented" << endl;
          abort();
        }
    }
  else
    {
      cout << "Case not implemented" << endl;
      abort();
    }
  
  vars.GetPolarization(polar);

  const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
  VectR3 s; VectComplex_wp feval, contrib;
  SetPoints<Dimension3> PointsElem; Complex_wp u_inc(0, 0);
  Complex_wp coef_eps(0, 0), coef_mu(0, 0);
  TinyVector<Complex_wp, 3> vec_u, vec_v;
  int offset_H = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
      
      int ref = vars.mesh.Element(i).GetReference();
      int num_f = -1;
      for (int j = 0; j < ref_cavity.GetM(); j++)
        if (ref == ref_cavity(j))
          num_f = j;
      
      if (num_f == -1)
        continue;

      VectComplex_wp& Fcplx = vecF(num_f);
      bool affine = vars.mesh.IsElementAffine(i);
      TinyMatrix<Complex_wp, Symmetric, 3, 3> epsilon, sigma, mu;
      epsilon = vars.ref_epsilon(ref).GetConstant();
      sigma = vars.ref_sigma(ref).GetConstant();
      mu = vars.ref_mu(ref).GetConstant();
      
      Real_wp eps_inf = vars.ref_drude(ref).eps_inf;
      vars.ref_drude(ref).ModifyCoefficientMaxwell(omega, epsilon, sigma, mu);
      
      bool presence_coef_eps = false;
      if (abs(epsilon(0, 0) - vars.epsilon0) > 1e-10)
        {
          presence_coef_eps = true;
          if (input_var.source_wei && vars.ref_drude(ref).IsEnabled())
            coef_eps = Iwp*omega*(eps_inf - vars.epsilon0);
          else
            coef_eps = Iwp*omega*(epsilon(0, 0) - vars.epsilon0);
        }

      bool presence_coef_mu = false;
      if (abs(mu(0, 0) - vars.mu0) > 1e-10)
        {
          presence_coef_mu = true;
          if (input_var.source_wei && vars.ref_drude(ref).IsEnabled())
            coef_mu = Iwp*omega*(eps_inf - vars.mu0);
          else
            coef_mu = Iwp*omega*(mu(0, 0) - vars.mu0);
        }
      
      bool presence_source = true;
      if (incident_source)
        {
          presence_source = false;
          if (presence_coef_eps)
            presence_source = true;

          if (presence_coef_mu)
            presence_source = true;

          if (presence_coef_eps && presence_coef_mu)
            {
              cout << "Case not implemented" << endl;
              abort();
            }
        }
      
      int offset = offset_H + vars.GetOffsetDofV(i);
      int offset_Q = vars.GetNbDof() - vars.GetNbVectorialDofDrude();
      if (vars.ref_drude(ref).IsEnabled())
        offset_Q += vars.OffsetDofDrudeV(i);
      
      if (presence_source)
        {
          vars.mesh.GetVerticesElement(i, s);
          Fb.FjElemQuadrature(s, PointsElem, vars.mesh, i);
          
          Real_wp jacob, jacob_weighted; Matrix3_3 dfjm1;	  
          int Nquad = Fb.GetNbPointsQuadratureInside();
          feval.Reallocate(3*Nquad);
          for (int j = 0; j < Nquad; j++)
            {
              vars.GetInverseJacobianMatrix(i, j, affine, Fb.GetGeometricElement(),
                                            dfjm1, jacob, jacob_weighted);
              
              typename Dimension::R_N x = PointsElem.GetPointQuadrature(j);
              Complex_wp srcE(0, 0);
              if (incident_source)
                {
                  if (presence_source)
                    {
                      Real_wp arg = omega*DotProd(k, x);
                      u_inc = Complex_wp(cos(arg), sin(arg));
                    }
                  
                  if (presence_coef_eps)
                    srcE = jacob_weighted*coef_eps*u_inc;

                  if (presence_coef_mu)
                    srcE = jacob_weighted*coef_mu*u_inc;
                  
                  if (input_var.source_wei && vars.ref_drude(ref).IsEnabled())
                    {
                      for (int kp = 0; kp < vars.ref_drude(ref).gamma.GetM(); kp++)
                        {
                          Complex_wp coef = vars.ref_drude(ref).eps_omega_p2(kp) - Iwp*omega*vars.ref_drude(ref).eps_sigma(kp);
                          
                          if (presence_coef_mu || presence_coef_eps)
                            {
                              Fcplx(offset_Q + 3*(kp*Nquad + j)) = coef*jacob_weighted*u_inc*polar(0);
                              Fcplx(offset_Q + 3*(kp*Nquad + j)+1) = coef*jacob_weighted*u_inc*polar(1);
                              Fcplx(offset_Q + 3*(kp*Nquad + j)+2) = coef*jacob_weighted*u_inc*polar(2);
                            }
                        }
                    }
                }
              else
                srcE = jacob_weighted*gaussian.GetAmplitude(x);
              
              if (presence_coef_eps)
                {
                  vec_u = srcE*polar;
                  Mlt(dfjm1, vec_u, vec_v);
                  feval(3*j) = vec_v(0); feval(3*j+1) = vec_v(1); feval(3*j+2) = vec_v(2);
                }
              else
                {
                  vec_u = srcE*polar;
                  Fcplx(offset + 3*j) = vec_u(0); Fcplx(offset + 3*j+1) = vec_u(1); 
                  Fcplx(offset + 3*j+2) = vec_u(2); 
                }
            }
          
          if (presence_coef_eps)
            {
              contrib.Reallocate(Fb.GetNbDof());
              Fb.ApplyCh(feval, contrib);
          
              mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, contrib, i);
              
              for (int j = 0; j < Fb.GetNbDof(); j++)
                Fcplx(mesh_num.Element(i).GetNumberDof(j)) += contrib(j);
            }
        }
    }

  for (int n = 0; n < vecF.GetM(); n++)
    vars.ImposeNullDirichletCondition(vecF(n));
}


template<>
void QnmSolver<HarmonicMaxwellEquation_3D>::ComputeSourceOmega(const Real_wp& omega, Vector<VectComplex_wp>& vecF)
{
  vecF.Reallocate(ref_cavity.GetM());
  for (int k = 0; k < vecF.GetM(); k++)
    {
      vecF(k).Reallocate(vars.GetNbDof());
      vecF(k).Zero();
    }
  
  // on calcule la source
  bool incident_source = false;
  R3 k, polar;
  GaussianSourceField<Real_wp, Dimension> gaussian;
  if (vars.GetSourceType(0) == vars.SRC_DIFFRACTED_FIELD)
    {
      incident_source = true;
      k = vars.GetWaveVector();
      Mlt(1.0/Norm2(k), k);
    }
  else if (vars.GetSourceType(0) == vars.SRC_VOLUME)
    {
      VectString param = vars.GetParameterSource(0, 0);
      int nb = -1;
      for (int k = 0; k < param.GetM(); k++)
        if (param(k) == "GAUSSIAN")
          {
            nb = k+1;
            vars.InitGaussianParameter(gaussian, param, nb);
          }
        
      if (nb == -1)
        {
          cout << "Case not implemented" << endl;
          abort();
        }
    }
  else
    {
      cout << "Case not implemented" << endl;
      abort();
    }
  
  vars.GetPolarization(polar);

  const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
  VectR3 s; VectComplex_wp feval, contrib;
  SetPoints<Dimension3> PointsElem; Complex_wp u_inc(0, 0);
  Complex_wp coef_eps(0, 0);
  TinyVector<Complex_wp, 3> vec_u, vec_v;
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
      
      int ref = vars.mesh.Element(i).GetReference();
      int num_f = -1;
      for (int j = 0; j < ref_cavity.GetM(); j++)
        if (ref == ref_cavity(j))
          num_f = j;
      
      if (num_f == -1)
        continue;

      VectComplex_wp& Fcplx = vecF(num_f);
      bool affine = vars.mesh.IsElementAffine(i);
      TinyMatrix<Complex_wp, Symmetric, 3, 3> epsilon, sigma, mu;
      epsilon = vars.ref_epsilon(ref).GetConstant();
      sigma = vars.ref_sigma(ref).GetConstant();
      mu = vars.ref_mu(ref).GetConstant();
      
      //Real_wp eps_inf = vars.ref_drude(ref).eps_inf;
      vars.ref_drude(ref).ModifyCoefficientMaxwell(omega, epsilon, sigma, mu);
      
      bool presence_coef_eps = false;
      if (abs(epsilon(0, 0) - vars.epsilon0) > 1e-10)
        {
          presence_coef_eps = true;
          coef_eps = Iwp*omega*(epsilon(0, 0) - vars.epsilon0);
        }
      
      bool presence_source = true;
      if (incident_source)
        {
          presence_source = false;
          if (presence_coef_eps)
            presence_source = true;
        }
      
      int offset_P = vars.GetNbDof() - vars.GetNbVectorialDofDrude();
      if (vars.ref_drude(ref).IsEnabled())
        offset_P += vars.OffsetDofDrudeV(i);
      
      if (presence_source)
        {
          vars.mesh.GetVerticesElement(i, s);
          Fb.FjElemQuadrature(s, PointsElem, vars.mesh, i);
          
          Real_wp jacob, jacob_weighted; Matrix3_3 dfjm1;	  
          int Nquad = Fb.GetNbPointsQuadratureInside();
          feval.Reallocate(3*Nquad);
          for (int j = 0; j < Nquad; j++)
            {
              vars.GetInverseJacobianMatrix(i, j, affine, Fb.GetGeometricElement(),
                                            dfjm1, jacob, jacob_weighted);
              
              typename Dimension::R_N x = PointsElem.GetPointQuadrature(j);
              Complex_wp srcE(0, 0);
              if (incident_source)
                {
                  if (presence_source)
                    {
                      Real_wp arg = omega*DotProd(k, x);
                      u_inc = Complex_wp(cos(arg), sin(arg));
                    }
                  
                  if (presence_coef_eps)
                    srcE = jacob_weighted*coef_eps*u_inc;
                }
              else
                srcE = jacob_weighted*gaussian.GetAmplitude(x);
              
              vec_u = srcE*polar;
              Mlt(dfjm1, vec_u, vec_v);
              feval(3*j) = vec_v(0); feval(3*j+1) = vec_v(1); feval(3*j+2) = vec_v(2);
             }
          
          contrib.Reallocate(Fb.GetNbDof());
          Fb.ApplyCh(feval, contrib);
          
          mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, contrib, i);
              
          for (int j = 0; j < Fb.GetNbDof(); j++)
            Fcplx(mesh_num.Element(i).GetNumberDof(j)) += contrib(j);
        }
    }

  for (int n = 0; n < vecF.GetM(); n++)
    vars.ImposeNullDirichletCondition(vecF(n));
}

template<>
void QnmSolver<LaplaceEquation<Dimension2> >::ComputeSourceOmega(const Real_wp& omega, Vector<VectComplex_wp>& vecF)
{
  vecF.Reallocate(ref_cavity.GetM());
  for (int k = 0; k < vecF.GetM(); k++)
    {
      vecF(k).Reallocate(vars.GetNbDof());
      vecF(k).Zero();
    }
  
  // on calcule la source
  bool incident_source = false;
  GaussianSourceField<Real_wp, Dimension> gaussian;
  IncidentWaveField<Complex_wp, Dimension>* incident_wave = NULL;
  if (vars.GetSourceType(0) == vars.SRC_DIFFRACTED_FIELD)
    {
      typename Dimension::R_N k, origin;
      incident_source = true;
      k = vars.GetWaveVector();
      origin = vars.GetPhaseOrigin();
      Mlt(omega/Norm2(k), k);
      VectString param = vars.GetParameterSource(0, 0);
      int ref = 2;
      Complex_wp rho = vars.ref_rho(ref).GetConstant();
      TinyMatrix<Complex_wp, Symmetric, 2, 2> mu;
      for (int m = 0; m < 2; m++)
        for (int n = 0; n < 2; n++)
          mu(m, n) = vars.ref_mu(ref).GetConstant()(m, n);
      
      Complex_wp sigma = vars.ref_sigma(ref).GetConstant();
      vars.ref_drude(ref).ModifyCoefficientHelmholtz(omega, rho, sigma, mu);

      if (vars.GetIncidentFieldType(0) == vars.INCIDENT_PLANE_WAVE)
        incident_wave = new PlaneWaveIncidentField<Complex_wp, Dimension2>(origin, k);
      else
        {
          param(3) = to_str(rho);
          param(4) = to_str(mu);            
          incident_wave = new LayeredPlaneWaveIncidentField<Complex_wp, Dimension2>(origin, k, omega, 1.0, 1.0, param);
        }
    }
  else if (vars.GetSourceType(0) == vars.SRC_VOLUME)
    {
      VectString param = vars.GetParameterSource(0, 0);
      int nb = -1;
      for (int p = 0; p < param.GetM(); p++)
        if (param(p) == "GAUSSIAN")
          {
            nb = p+1;
            vars.InitGaussianParameter(gaussian, param, nb);
          }
        
      if (nb == -1)
        {
          cout << "Case not implemented" << endl;
          abort();
        }
    }
  else
    {
      cout << "Case not implemented" << endl;
      abort();
    }

  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  typename Dimension::VectR_N s;
  SetPoints<Dimension> PointsElem; Complex_wp u_inc(0, 0);
  TinyVector<Complex_wp, 2> grad_uinc;
  Complex_wp coef_rho(0, 0), coef_mu(0, 0);
  int offset_v = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension, 1>& Fb = vars.GetReferenceElementH1(i);
	
      int offset_Q = vars.GetNbDof() - vars.GetNbVectorialDofDrude();
      int offset_P = offset_Q - vars.GetNbVectorialDofDrude();
      offset_P += vars.OffsetDofDrudeV(i);
      offset_Q += vars.OffsetDofDrudeV(i);
      
      int ref = vars.mesh.Element(i).GetReference();
      int num_f = -1;
      for (int j = 0; j < ref_cavity.GetM(); j++)
        if (ref == ref_cavity(j))
          num_f = j;
      
      if (num_f == -1)
        continue;
      
      VectComplex_wp& Fcplx = vecF(num_f);
      bool affine = vars.mesh.IsElementAffine(i);
      Complex_wp rho = vars.ref_rho(ref).GetConstant();
      TinyMatrix<Complex_wp, Symmetric, Dimension::dim_N, Dimension::dim_N> mu;
      for (int m = 0; m < Dimension::dim_N; m++)
        for (int n = 0; n < Dimension::dim_N; n++)
          mu(m, n) = vars.ref_mu(ref).GetConstant()(m, n);
      
      Complex_wp sigma = vars.ref_sigma(ref).GetConstant();
      Real_wp rho_inf = vars.ref_drude(ref).eps_inf;
      vars.ref_drude(ref).ModifyCoefficientHelmholtz(omega, rho, sigma, mu);
      
      Complex_wp rho0 = vars.rho0, mu0 = vars.mu0;
      vars.mesh.GetVerticesElement(i, s);
      if (incident_source)
        {
          incident_wave->InitElement(i, s);
          incident_wave->UpdateCoefAB(mu0, rho0);
        }
      
      bool presence_coef_rho = false;
      if (abs(rho - rho0) > 1e-10)
        {
          presence_coef_rho = true;
          if (input_var.source_wei && vars.ref_drude(ref).IsEnabled())
            coef_rho = Iwp*omega*(rho_inf - rho0);
          else
            coef_rho = Iwp*omega*(rho - rho0);
        }
      
      bool presence_coef_mu = false;
      if (abs(mu(0, 0) - mu0) > 1e-10)
        {
          presence_coef_mu = true;
          coef_mu = -mu0*(Real_wp(1) / mu(0,0) - Real_wp(1)/mu0);
        }
      
      bool presence_source = true;
      if (incident_source)
        {
          presence_source = false;
          if (presence_coef_rho || presence_coef_mu)
            presence_source = true;
        }
      
      if (presence_source)
        {

          Fb.FjElemQuadrature(s, PointsElem, vars.mesh, i);
	  
          int offset = offset_v + vars.GetOffsetDofV(i);
          int Nquad = Fb.GetNbPointsQuadratureInside();
          for (int j = 0; j < Nquad; j++)
            {
              Real_wp jacob = vars.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
              typename Dimension::R_N x = PointsElem.GetPointQuadrature(j);
              if (incident_source)
                {
                  if (presence_source)
                    incident_wave->EvaluateFunctionGradient(x, u_inc, grad_uinc);
                  
                  if (presence_coef_rho)
                    {
                      if (!Fb.LumpedMassMatrix()) { cout << "Not implemented" << endl; abort(); }
                      int n = mesh_num.Element(i).GetNumberDof(j);
                      Fcplx(n) += jacob*coef_rho*u_inc;
                    }
                  
                  if (presence_coef_mu)
                    {
                      Fcplx(offset + 2*j) += jacob*coef_mu*grad_uinc(0);
                      Fcplx(offset + 2*j+1) += jacob*coef_mu*grad_uinc(1);
                    }
                  
                  if (input_var.source_wei && vars.ref_drude(ref).IsEnabled() && presence_coef_rho)
                    {
                      int nPole = vars.ref_drude(ref).gamma.GetM();
                      for (int kp = 0; kp < nPole; kp++)
                        {
                          Complex_wp coef = vars.ref_drude(ref).eps_omega_p2(kp) - Iwp*omega*vars.ref_drude(ref).eps_sigma(kp);
                          Fcplx(offset_Q + j*nPole + kp) += jacob*coef*u_inc;
                        }
                    }
                }
              else
                {
                  if (!Fb.LumpedMassMatrix()) { cout << "Not implemented" << endl; abort(); }
                  int n = mesh_num.Element(i).GetNumberDof(j);
                  Fcplx(n) += jacob*gaussian.GetAmplitude(x);
                }
            }
        }
    }

  if (vars.GetSourceType(0) == vars.SRC_DIFFRACTED_FIELD)
    delete incident_wave;
  
  for (int n = 0; n < vecF.GetM(); n++)
    vars.ImposeNullDirichletCondition(vecF(n));
}

template<>
void QnmSolver<HelmholtzEquation<Dimension2> >::ComputeSourceOmega(const Real_wp& omega, Vector<VectComplex_wp>& vecF)
{
  vecF.Reallocate(ref_cavity.GetM());
  for (int k = 0; k < vecF.GetM(); k++)
    {
      vecF(k).Reallocate(vars.GetNbDof());
      vecF(k).Zero();
    }
  
  // on calcule la source
  bool incident_source = false;
  typename Dimension::R_N k;
  GaussianSourceField<Real_wp, Dimension> gaussian;
  if (vars.GetSourceType(0) == vars.SRC_DIFFRACTED_FIELD)
    {
      incident_source = true;
      k = vars.GetWaveVector();
      Mlt(1.0/Norm2(k), k);
    }
  else if (vars.GetSourceType(0) == vars.SRC_VOLUME)
    {
      VectString param = vars.GetParameterSource(0, 0);
      int nb = -1;
      for (int k = 0; k < param.GetM(); k++)
        if (param(k) == "GAUSSIAN")
          {
            nb = k+1;
            vars.InitGaussianParameter(gaussian, param, nb);
          }
        
      if (nb == -1)
        {
          cout << "Case not implemented" << endl;
          abort();
        }
    }
  else
    {
      cout << "Case not implemented" << endl;
      abort();
    }

  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  typename Dimension::VectR_N s;
  SetPoints<Dimension> PointsElem; Complex_wp u_inc(0, 0);
  Complex_wp coef_rho(0, 0);
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension, 1>& Fb = vars.GetReferenceElementH1(i);
	
      int offset_P = vars.GetNbDof() - vars.GetNbVectorialDofDrude();
      offset_P += vars.OffsetDofDrudeV(i);
      
      int ref = vars.mesh.Element(i).GetReference();
      int num_f = -1;
      for (int j = 0; j < ref_cavity.GetM(); j++)
        if (ref == ref_cavity(j))
          num_f = j;
      
      if (num_f == -1)
        continue;
      
      VectComplex_wp& Fcplx = vecF(num_f);
      bool affine = vars.mesh.IsElementAffine(i);
      Complex_wp rho = vars.ref_rho(ref).GetConstant();
      TinyMatrix<Complex_wp, Symmetric, Dimension::dim_N, Dimension::dim_N> mu;
      for (int m = 0; m < Dimension::dim_N; m++)
        for (int n = 0; n < Dimension::dim_N; n++)
          mu(m, n) = vars.ref_mu(ref).GetConstant()(m, n);
      
      Complex_wp sigma = vars.ref_sigma(ref).GetConstant();
      //Real_wp rho_inf = vars.ref_drude(ref).eps_inf;
      vars.ref_drude(ref).ModifyCoefficientHelmholtz(omega, rho, sigma, mu);
      
      bool presence_coef_rho = false;
      if (abs(rho - vars.rho0) > 1e-10)
        {
          presence_coef_rho = true;
          coef_rho = Iwp*omega*(rho - vars.rho0);
        }

      bool presence_source = true;
      if (incident_source)
        {
          presence_source = false;
          if (presence_coef_rho)
            presence_source = true;
        }
      
      if (presence_source)
        {
          vars.mesh.GetVerticesElement(i, s);
          Fb.FjElemQuadrature(s, PointsElem, vars.mesh, i);
	  
          int Nquad = Fb.GetNbPointsQuadratureInside();
          for (int j = 0; j < Nquad; j++)
            {
              Real_wp jacob = vars.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
              typename Dimension::R_N x = PointsElem.GetPointQuadrature(j);
              if (incident_source)
                {
                  if (presence_source)
                    {
                      Real_wp arg = omega*DotProd(k, x);
                      u_inc = Complex_wp(cos(arg), sin(arg));
                    }
                  
                  if (presence_coef_rho)
                    {
                      if (!Fb.LumpedMassMatrix()) { cout << "Not implemented" << endl; abort(); }
                      int n = mesh_num.Element(i).GetNumberDof(j);
                      Fcplx(n) += jacob*coef_rho*u_inc;
                    }
                }
              else
                {
                  if (!Fb.LumpedMassMatrix()) { cout << "Not implemented" << endl; abort(); }
                  int n = mesh_num.Element(i).GetNumberDof(j);
                  Fcplx(n) += jacob*gaussian.GetAmplitude(x);
                }
            }
        }
    }
  
  for (int n = 0; n < vecF.GetM(); n++)
    vars.ImposeNullDirichletCondition(vecF(n));
}

template<>
Real_wp QnmSolver<StaticMaxwellEquation_3D>
::ComputeDisplayFieldE(int m, const Real_wp& omega, const VectComplex_wp& Y,
                       const Vector<Vector<VectComplex_wp> >& GradQuad, VectComplex_wp& rhs)
{
  const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
  Real_wp err = 0, norme_err = 0;
  Vector<Vector<TinyVector<Complex_wp, 3> > > Equad(this->vars.mesh.GetNbElt());
  Vector<Vector<TinyVector<Complex_wp, 3> > > Efem_quad(this->vars.mesh.GetNbElt());
  int offset_H = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
  for (int i = 0; i < this->vars.mesh.GetNbElt(); i++)
    {
      int offset = offset_H + vars.GetOffsetDofV(i);
      const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
      
      bool affine = vars.mesh.IsElementAffine(i);
      int Nquad = Fb.GetNbPointsQuadratureInside();
      Equad(i).Reallocate(Nquad);
      Efem_quad(i).Reallocate(Nquad);
      
      int ref = vars.mesh.Element(i).GetReference();
      Complex_wp eps = vars.ref_epsilon(ref).GetConstant()(0, 0);
      Real_wp eps_inf = 1.0;
      if (vars.ref_drude(ref).IsEnabled())
        {
          eps = vars.ref_drude(ref).GetEpsilon(omega);
          eps_inf = vars.ref_drude(ref).eps_inf;
        }
      
      Complex_wp coef = Real_wp(1) / (-Iwp*omega*eps);
      for (int j = 0; j < Nquad; j++)
        {
          Real_wp jacob = vars.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
          Real_wp invPoids = Real_wp(1) / jacob;;
          Complex_wp coefQ = invPoids;
          if (input_var.source_wei && vars.ref_drude(ref).IsEnabled())
            coefQ *= (eps - vars.mu0) / (eps_inf - vars.mu0);
          
          Equad(i)(j)(0) = coef*(-GradQuad(i)(0)(j) + rhs(offset + 3*j)*coefQ);
          Equad(i)(j)(1) = coef*(-GradQuad(i)(1)(j) + rhs(offset + 3*j+1)*coefQ);
          Equad(i)(j)(2) = coef*(-GradQuad(i)(2)(j) + rhs(offset + 3*j+2)*coefQ);

          Efem_quad(i)(j)(0) = coef*(-GradQuad(i)(3)(j) + rhs(offset + 3*j)*coefQ);
          Efem_quad(i)(j)(1) = coef*(-GradQuad(i)(4)(j) + rhs(offset + 3*j+1)*coefQ);
          Efem_quad(i)(j)(2) = coef*(-GradQuad(i)(5)(j) + rhs(offset + 3*j+2)*coefQ);
          
          if (!vars.InsidePML(i))
            {
              err += jacob*AbsSquare(Equad(i)(j) - Efem_quad(i)(j));
              norme_err += jacob*AbsSquare(Efem_quad(i)(j));
            }
        }
    }
  
  err = sqrt(err / norme_err);
  
  if (vars.var_grid.GetM() <= 0)
    return err;
  
  GridInterpolationFull<Dimension3>& var_gr = vars.var_grid(0);
  GridInterpolation<Dimension3>& var_interp = vars.all_points_display;
  const IVect& list_points = var_gr.GetPointNumber();
  const Vector<IVect>& list_points_proc = var_gr.GetPointNumberAllProc();

  int nnz = list_points.GetM();
  int Nglob = var_gr.GetNbGlobalPoints();
  TinyVector<VectComplex_wp, 3> interpE, exactE;
  for (int j = 0; j < 3; j++)
    {
      interpE(j).Reallocate(Nglob);
      exactE(j).Reallocate(Nglob);
      interpE(j).Zero(); exactE(j).Zero();
    }

  VectReal_wp phi;
  for (int i1 = 0; i1 < nnz; i1++)
    {
      int i = list_points(i1);
      if (i < 0)
        continue;
      
      //R3 pt_glob = var_interp.GetGlobalCoordinate(i);
      int iglob = list_points_proc(0)(i1);
      int num_elem = var_interp.GetElementNumber(i);
      R3 pt_loc = var_interp.GetLocalCoordinate(i);
      vars.mesh.ComputeValuesPhiNodalRef(num_elem, pt_loc, phi);
      
      TinyVector<Complex_wp, 3> vec_E, vec_Efem;
      if (num_elem >= 0)
        {
          for (int j = 0; j < phi.GetM(); j++)
            {
              Add(phi(j), Equad(num_elem)(j), vec_E);
              Add(phi(j), Efem_quad(num_elem)(j), vec_Efem);
            }
          
          interpE(0)(iglob) = vec_E(0);
          interpE(1)(iglob) = vec_E(1);
          interpE(2)(iglob) = vec_E(2);
          
          exactE(0)(iglob) = vec_Efem(0);
          exactE(1)(iglob) = vec_Efem(1);
          exactE(2)(iglob) = vec_Efem(2);
        }
    }
  
  WriteMatlab(interpE(0), var_gr, "FieldE" + to_str(m) + "_U0.dat", OutputTypeEnum::DOUBLE_PRECISION, false);
  WriteMatlab(interpE(1), var_gr, "FieldE"  + to_str(m) + "_U1.dat", OutputTypeEnum::DOUBLE_PRECISION, false);
  WriteMatlab(interpE(2), var_gr, "FieldE"  + to_str(m) + "_U2.dat", OutputTypeEnum::DOUBLE_PRECISION, false);

  WriteMatlab(exactE(0), var_gr, "FieldEfem"  + to_str(m) + "_U0.dat", OutputTypeEnum::DOUBLE_PRECISION, false);
  WriteMatlab(exactE(1), var_gr, "FieldEfem" + to_str(m) + "_U1.dat", OutputTypeEnum::DOUBLE_PRECISION, false);
  WriteMatlab(exactE(2), var_gr, "FieldEfem" + to_str(m) + "_U2.dat", OutputTypeEnum::DOUBLE_PRECISION, false);
  
  return err;
}

template<>
Real_wp QnmSolver<HarmonicMaxwellEquation_3D>
::ComputeDisplayFieldE(int m, const Real_wp& omega, const VectComplex_wp& Y,
                       const Vector<Vector<VectComplex_wp> >& GradQuad, VectComplex_wp& rhs)
{
  return Real_wp(0);
}

template<>
Real_wp QnmSolver<LaplaceEquation<Dimension2> >
::ComputeDisplayFieldE(int m, const Real_wp& omega, const VectComplex_wp& Y,
                       const Vector<Vector<VectComplex_wp> >& GradQuad, VectComplex_wp& rhs)
{
  return Real_wp(0);
}

template<>
Real_wp QnmSolver<HelmholtzEquation<Dimension2> >
::ComputeDisplayFieldE(int m, const Real_wp& omega, const VectComplex_wp& Y,
                       const Vector<Vector<VectComplex_wp> >& GradQuad, VectComplex_wp& rhs)
{
  return Real_wp(0);
}

template<>
void QnmSolver<StaticMaxwellEquation_3D>
::ComputeBiorthogonal(const Complex_wp& lambda, const Vector<Complex_wp>& eigen_vec,
                      Vector<Complex_wp>& biortho, bool compute_Q, bool compute_pml)
{
  const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
  int offset_Es = mesh_num.GetNbDof();
  int offset_H = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();

  biortho.Reallocate(eigen_vec.GetM());
  biortho.Zero();

  // a = E^*
  for (int i = 0; i < num_Es.GetM(); i++)
    biortho(i) = eigen_vec(num_Es(i));
  
  if (!compute_pml && !compute_Q)
    {
      for (int i = 0; i < vars.mesh.GetNbElt(); i++)
        if (!vars.InsidePML(i))
          {
            int offset = offset_H + vars.GetOffsetDofV(i);
            const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
            int Nquad = Fb.GetNbPointsQuadratureInside();
            // b = -H^*
            for (int j = 0; j < 3*Nquad; j++)
              biortho(offset + j) = -eigen_vec(offset + j);          
          }

      return;
    }
  
  if ((Dh_pml.GetM() > 0) && (compute_pml))
    {
      // on calcule c dans le cas non-diagonal
      // on factorise (-lambda Mh + T_{3,1,2})
      Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
      Copy(T312_pml, A);
      Add(-lambda, Dh_pml, A);
      
      All_MatrixLU<Complex_wp> mat_lu;
      mat_lu.Factorize(A);

      int nEs = offset_H - offset_Es;
      VectComplex_wp x(nEs), y(nEs);
      x.Zero(); y.Zero();
      for (int i = 0; i < offset_Es; i++)
        if (mesh_num.GetDofPML(i) >= 0)
          x(mesh_num.GetDofPML(i)) = eigen_vec(i);

      Mlt(T231_pml, x, y);
      MltAdd(-lambda, Dh_pml, x, Complex_wp(1, 0), y);
      
      mat_lu.Solve(y);
      for (int i = 0; i < nEs; i++)
        biortho(offset_Es + i) = y(i);
    }
  
  int offset = 0;
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      offset = offset_H + vars.GetOffsetDofV(i);
      const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
      int Nquad = Fb.GetNbPointsQuadratureInside();
      int offset_s = offset;
      if (vars.InsidePML(i))
        offset_s += 3*Nquad;
      
      // b = -H^*
      for (int j = 0; j < 3*Nquad; j++)
        biortho(offset + j) = -eigen_vec(offset_s + j);

      if (vars.InsidePML(i) && compute_pml)
        {
          HarmonicMaxwell3D_PhysGeomInfo<Real_wp>& mass = vars.Glob_matMass_elem(i);
          bool diag_mass = mass.IsDiagonalMass();
          
          int i1 = i - vars.mesh.GetNbElt() + vars.GetNbEltPML();
          int ref = vars.mesh.Element(i).GetReference();
          Real_wp eps = vars.ref_epsilon(ref).GetConstant()(0, 0);
          Real_wp mu = vars.ref_mu(ref).GetConstant()(0, 0);

          // c = epsilon (1 + (tau_{2, 3, 1} - tau_{3, 1, 2}) / (-lambda + T_{3, 1, 2}) ) E
          if (diag_mass && (Dh_pml.GetM() <= 0))
            for (int j = 0; j < Fb.GetNbDof(); j++)
              {
                int num_dof = mesh_num.Element(i).GetNumberDof(j);
                if (mesh_num.GetDofPML(num_dof) >= 0)
                  {
                    int num_dof_Es = mesh_num.GetDofPML(num_dof) + offset_Es;
                    TinyVector<Real_wp, 3> tau = mass.GetTauPML(j);
                    int n = Fb.GetCoordinateDof(j);
                    Complex_wp vloc(0, 0);
                    switch (n)
                      {
                      case 0: vloc = (tau(1)-tau(2)) / (-lambda + tau(2)); break;
                      case 1: vloc = (tau(2)-tau(0)) / (-lambda + tau(0)); break;
                      case 2: vloc = (tau(0)-tau(1)) / (-lambda + tau(1)); break;
                      }
                    
                    biortho(num_dof_Es) = eps*(Real_wp(1) + vloc)*eigen_vec(num_dof);              
                  }
              }
          
          // d = -mu (1 + (tau_{2, 3, 1} - tau_{3, 1, 2}) / (-lambda + T_{3, 1, 2}) ) H
          for (int j = 0; j < Nquad; j++)
            {
              TinyVector<Real_wp, 3> tau = vars.GetTauPML(i1, j);
              Complex_wp vloc = (tau(1)-tau(2)) / (-lambda + tau(2));
              biortho(offset_s + 3*j) = -mu*(Real_wp(1) + vloc)*eigen_vec(offset+3*j);
              
              vloc = (tau(2)-tau(0)) / (-lambda + tau(0));
              biortho(offset_s + 3*j+1) = -mu*(Real_wp(1) + vloc)*eigen_vec(offset+3*j+1);

              vloc = (tau(0)-tau(1)) / (-lambda + tau(1));
              biortho(offset_s + 3*j+2) = -mu*(Real_wp(1) + vloc)*eigen_vec(offset+3*j+2);
            }
        }

      int ref = vars.mesh.Element(i).GetReference();
      if (vars.ref_drude(ref).IsEnabled() && compute_Q)
        for (int kp = 0; kp < vars.ref_drude(ref).gamma.GetM(); kp++)
          {
            Complex_wp coef1(1, 0), coef2(1, 0);
            if (!vars.use_symm_drude)
              {
                Real_wp coef_omega2 = vars.ref_drude(ref).omega_02(kp);
                Real_wp coef_eps_inf = vars.ref_drude(ref).eps_omega_p2(kp);
                Real_wp sig = vars.ref_drude(ref).eps_sigma(kp);
                Complex_wp denom = coef_eps_inf - sig*lambda;
                if (abs(denom) <= 1e-12)
                  {
                    DISP(vars.ref_drude(ref).gamma);
                    DISP(vars.ref_drude(ref).eps_omega_p2);
                    DISP(lambda); DISP(sig*lambda); DISP(coef_eps_inf); DISP(denom);
                    cout << "Null denominator for biortho" << endl;
                    abort();
                  }
                
                coef1 = Real_wp(1)/denom;
                coef2 = coef_omega2 * coef1;
              }
            
            if (vars.ref_drude(ref).IsModeTE())
              {
                coef1 = -coef1;
                coef2 = -coef2;
              }
            
            // biorthogonal for p and q
            int offset_Q = vars.GetNbDof() - vars.GetNbVectorialDofDrude();
            int offset_P = offset_Q - vars.GetNbVectorialDofDrude();
            offset_P += vars.OffsetDofDrudeV(i) + 3*Nquad*kp;
            offset_Q += vars.OffsetDofDrudeV(i) + 3*Nquad*kp;
            for (int j = 0; j < 3*Nquad; j++)
              biortho(offset_Q+j) = -coef1*eigen_vec(offset_Q+j);

            for (int j = 0; j < 3*Nquad; j++)
              biortho(offset_P+j) = coef2*eigen_vec(offset_P+j);
          }
    }
  
}

template<>
void QnmSolver<HarmonicMaxwellEquation_3D>
::ComputeBiorthogonal(const Complex_wp& lambda, const Vector<Complex_wp>& eigen_vec,
                      Vector<Complex_wp>& biortho, bool compute_Q, bool compute_pml)
{
  // E is identical
  biortho = eigen_vec;

  // we modify P
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      int offset_P = vars.GetNbDof() - vars.GetNbVectorialDofDrude();
      offset_P += vars.OffsetDofDrudeV(i);
      int ref = vars.mesh.Element(i).GetReference();
      if (vars.ref_drude(ref).IsEnabled() && compute_Q)
        {
          int nPole = vars.ref_drude(ref).gamma.GetM();
          const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
          for (int kp = 0; kp < nPole; kp++)
            {
              Real_wp coef_eps_inf = vars.ref_drude(ref).eps_omega_p2(kp);
              Real_wp sig = vars.ref_drude(ref).eps_sigma(kp);
              Complex_wp coef = -lambda*lambda / (coef_eps_inf - lambda*sig);
              if (vars.ref_drude(ref).IsModeTM())
                {
                  int Nquad = Fb.GetNbPointsQuadratureInside();
                  int offset = offset_P + 3*Nquad*kp;
                  for (int j = 0; j < 3*Nquad; j++)
                    biortho(offset+j) = coef*eigen_vec(offset+j);
                }
            }
        }
    }
}

template<>
void QnmSolver<LaplaceEquation<Dimension2> >
::ComputeBiorthogonal(const Complex_wp& lambda, const Vector<Complex_wp>& eigen_vec,
                      Vector<Complex_wp>& biortho, bool compute_Q, bool compute_pml)
{
  // without dispersive PMLs we have the equality
  biortho = eigen_vec;
  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  Complex_wp one; SetComplexOne(one);
  int Ndof = mesh_num.GetNbDof();
  Vector<bool> DofTreated(Ndof);
  DofTreated.Fill(false);
  
  // we compute u_ortho = (1.0 - (sigma_x + sigma_y) / (2 lambda)) u
  // we reconstruct u_ortho^* = - (sigma_x - sigma_y) / (-2 lambda + sigma_x + sigma_y) u_ortho
  if (compute_pml)
    for (int i = 0; i < vars.mesh.GetNbElt(); i++)
      {
        const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
        if (!Fb.LumpedMassMatrix())
	  {
	    cout << "Not coherent for elements without mass lumping" << endl;
	    abort();
	  }
        
        int N   = Fb.GetNbPointsQuadratureInside();
        int i1 = i - vars.mesh.GetNbElt() + vars.GetNbEltPML();
        for (int j = 0; j < N ; j++)
          {
            int num_dof = mesh_num.Element(i).GetNumberDof(j);
            if (vars.InsidePML(i))
              {
                if (!DofTreated(num_dof))
                  {
                    DofTreated(num_dof) = true;
                    R2 tau = vars.GetTauPML(i1, j);
                    biortho(num_dof) =
                      (one - (tau(0) + tau(1))/(Real_wp(2)*lambda)) * eigen_vec(num_dof);
                    
                    int num_dof_pml = mesh_num.GetDofPML(num_dof);
                    //biortho(Ndof + num_dof_pml) = -(tau(0) - tau(1)) / (-Real_wp(2)*lambda + tau(0) + tau(1))*biortho(num_dof);
                    biortho(Ndof + num_dof_pml) = (tau(0) - tau(1)) / (Real_wp(2)*lambda) * eigen_vec(num_dof);
                  }
              }
          }
      }
  
  // we reconstruct v_ortho
  int offset_v = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      int offset_Q = vars.GetNbDof() - vars.GetNbVectorialDofDrude();
      int offset_P = offset_Q - vars.GetNbVectorialDofDrude();
      offset_P += vars.OffsetDofDrudeV(i);
      offset_Q += vars.OffsetDofDrudeV(i);
      int offset = offset_v + vars.GetOffsetDofV(i);
      bool affine = vars.mesh.IsElementAffine(i);
      const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
      if (vars.InsidePML(i) && compute_pml)
        {
          int i1 = i - vars.mesh.GetNbElt() + vars.GetNbEltPML();
          Vector<Complex_wp> u_loc(Fb.GetNbDof()), us_loc(Fb.GetNbDof());
          const IVect& Nodle = mesh_num.Element(i).GetNodle();
          for (int j = 0; j < Fb.GetNbDof(); j++)
            {
              u_loc(j) = biortho(Nodle(j));
              us_loc(j) = biortho(Ndof + mesh_num.GetDofPML(Nodle(j)));
            }
          
          Vector<Complex_wp> grad_u_loc(2*Fb.GetNbDof()), grad_us_loc(2*Fb.GetNbDof());
          Fb.ApplyRhTranspose(u_loc, grad_u_loc);
          Fb.ApplyRhTranspose(us_loc, grad_us_loc);
	  
          TinyVector<Complex_wp, 2> vec_u, vec_v, vec_us, vec_vs;
          for (int j = 0; j < Fb.GetNbDof(); j++)
            {
              vec_u.Init(grad_u_loc(2*j), grad_u_loc(2*j+1));
              vec_us.Init(grad_us_loc(2*j), grad_us_loc(2*j+1));
              Real_wp jacob(1);
              if (affine)
                {
                  jacob = vars.Glob_jacobian(i)(0);
                  MltTrans(vars.Glob_DFjm1(i)(0), vec_u, vec_v);
                  MltTrans(vars.Glob_DFjm1(i)(0), vec_us, vec_vs);
                }
              else
                {
                  jacob = vars.Glob_jacobian(i)(j) / Fb.WeightsND(j);
                  MltTrans(vars.Glob_DFjm1(i)(j), vec_u, vec_v);
                  MltTrans(vars.Glob_DFjm1(i)(j), vec_us, vec_vs);
                }
              
              R2 tau = vars.GetTauPML(i1, j);
              
              biortho(offset+2*j) = -(vec_v(0) + vec_vs(0)) / (jacob*(-lambda + tau(0)));
              biortho(offset+2*j+1) = -(vec_v(1) - vec_vs(1)) / (jacob*(-lambda + tau(1)));
            }
        }
      else
        {
          // -H
          for (int j = 0; j < Fb.GetNbDof(); j++)
            {
              biortho(offset+2*j) = -biortho(offset+2*j);
              biortho(offset+2*j+1) = -biortho(offset+2*j+1);
            }

          int ref = vars.mesh.Element(i).GetReference();
          if (vars.ref_drude(ref).IsEnabled() && compute_Q)
            {
              int nPole = vars.ref_drude(ref).gamma.GetM();
              for (int kp = 0; kp < nPole; kp++)
                {
                  Real_wp coef_omega2 = vars.ref_drude(ref).omega_02(kp);
                  Real_wp coef_eps_inf = vars.ref_drude(ref).eps_omega_p2(kp);
                  Real_wp sig = vars.ref_drude(ref).eps_sigma(kp);
                  Complex_wp denom = coef_eps_inf - sig*lambda;
                  if (abs(denom) <= 1e-12)
                    {
                      cout << "Null denominator for biortho" << endl;
                      abort();
                    }
                  
                  Complex_wp coef1 = Real_wp(1)/denom;
                  Complex_wp coef2 = coef_omega2 * coef1;
                  if (vars.ref_drude(ref).IsModeTE())
                    {
                      coef1 = -coef1;
                      coef2 = -coef2;
                    }
                  
                  // biorthogonal for p and q
                  if (vars.ref_drude(ref).IsModeTE())
                    {
                      for (int j = 0; j < Fb.GetNbDof(); j++)
                        {
                          int dof_Px = offset_P + 2*j*nPole + kp, dof_Py = dof_Px + nPole;
                          biortho(dof_Px) = coef2*eigen_vec(dof_Px);
                          biortho(dof_Py) = coef2*eigen_vec(dof_Py);
                          
                          int dof_Qx = offset_Q + 2*j*nPole + kp, dof_Qy = dof_Qx + nPole;
                          biortho(dof_Qx) = -coef1*eigen_vec(dof_Qx);
                          biortho(dof_Qy) = -coef1*eigen_vec(dof_Qy);
                        }
                    }
                  else
                    {
                      for (int j = 0; j < Fb.GetNbDof(); j++)
                        {
                          int dof_P = offset_P + j*nPole + kp;
                          int dof_Q = offset_Q + j*nPole + kp;
                          biortho(dof_P) = coef2*eigen_vec(dof_P);
                          biortho(dof_Q) = -coef1*eigen_vec(dof_Q);
                        }
                    }
                }
            }
        }
    }
}

template<>
void QnmSolver<HelmholtzEquation<Dimension2> >
::ComputeBiorthogonal(const Complex_wp& lambda, const Vector<Complex_wp>& eigen_vec,
                      Vector<Complex_wp>& biortho, bool compute_Q, bool compute_pml)
{
  // E is identical
  biortho = eigen_vec;

  // we modify P
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      int offset_P = vars.GetNbDof() - vars.GetNbVectorialDofDrude();
      offset_P += vars.OffsetDofDrudeV(i);
      int ref = vars.mesh.Element(i).GetReference();
      if (vars.ref_drude(ref).IsEnabled() && compute_Q)
        {
          int nPole = vars.ref_drude(ref).gamma.GetM();
          const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
          for (int kp = 0; kp < nPole; kp++)
            {
              Real_wp coef_eps_inf = vars.ref_drude(ref).eps_omega_p2(kp);
              Real_wp sig = vars.ref_drude(ref).eps_sigma(kp);
              Complex_wp coef = -lambda*lambda / (coef_eps_inf - lambda*sig);
              if (vars.ref_drude(ref).IsModeTM())
                {
                  for (int j = 0; j < Fb.GetNbDof(); j++)
                    {
                      int dof_P = offset_P + j*nPole + kp;
                      biortho(dof_P) = coef*eigen_vec(dof_P);
                    }
                }
            }
        }
    }
}

template<>
bool QnmSolver<StaticMaxwellEquation_3D>::DiagonalMassMatrix() const { return false; }

template<>
bool QnmSolver<StaticMaxwellEquation_3D>::BlockDiagonalMassMatrix() const
{
  return false;
  
  bool diag_mass = true;
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    if (!vars.Glob_matMass_elem(i).IsDiagonalMass())
      diag_mass = false;
  
  return diag_mass; 
}

template<>
bool QnmSolver<HarmonicMaxwellEquation_3D>::DiagonalMassMatrix() const { return false; }

template<>
bool QnmSolver<HarmonicMaxwellEquation_3D>::BlockDiagonalMassMatrix() const
{
  return false;
}
template<>
bool QnmSolver<LaplaceEquation<Dimension2> >::DiagonalMassMatrix() const { return false; }

template<>
bool QnmSolver<LaplaceEquation<Dimension2> >::BlockDiagonalMassMatrix() const { return true; }

template<>
bool QnmSolver<HelmholtzEquation<Dimension2> >::DiagonalMassMatrix() const { return false; }

template<>
bool QnmSolver<HelmholtzEquation<Dimension2> >::BlockDiagonalMassMatrix() const { return true; }

template<>
void QnmSolver<StaticMaxwellEquation_3D>::PerformInitializations(FemMatrixFreeClass_Base<Real_wp>& Mh0, FemMatrixFreeClass_Base<Real_wp>& Kh0)
{
  bool extract_matrix = false;
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    if (vars.InsidePML(i))
      {
          HarmonicMaxwell3D_PhysGeomInfo<Real_wp>& mass = vars.Glob_matMass_elem(i);
          bool diag_mass = mass.IsDiagonalMass();
          if (!diag_mass)
            extract_matrix = true;
      }

  if (!extract_matrix)
    return;
  
  const Matrix<Real_wp, General, ArrayRowSparse>& Mh = Mh0.mat_iterative_unsym;
  const Matrix<Real_wp, General, ArrayRowSparse>& Kh = Kh0.mat_iterative_unsym;
  
  //Mh.WriteText("Mh2.dat");
  //Kh.WriteText("Kh2.dat");
  const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);  
  int offset_Es = mesh_num.GetNbDof();
  int offset_H = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
  
  // on stocke les sous-matrices 
  GetSubMatrix(Mh, offset_Es, offset_H, offset_Es, offset_H, Dh_pml);
  GetSubMatrix(Kh, offset_Es, offset_H, offset_Es, offset_H, T312_pml);
  
  int N = mesh_num.GetNbDofPML();
  T231_pml.Reallocate(N, N);
  //ofstream file_out("num_pml.dat");
  for (int i = 0; i < offset_Es; i++)
    if (mesh_num.GetDofPML(i) >= 0)
      {
        //file_out << i << '\n';
        int iloc = mesh_num.GetDofPML(i);
        for (int j = 0; j < Kh.GetRowSize(i); j++)
          if ((Kh.Index(i, j) < offset_Es) && (mesh_num.GetDofPML(Kh.Index(i, j)) >= 0))
            {
              int jloc = mesh_num.GetDofPML(Kh.Index(i, j));
              T231_pml.AddInteraction(iloc, jloc, Kh.Value(i, j));
            }
      }

  //file_out.close();
  //Dh_pml.WriteText("Dh.dat");
  //T312_pml.WriteText("T312.dat");
  //T231_pml.WriteText("T231.dat");
}

template<>
void QnmSolver<HarmonicMaxwellEquation_3D>::PerformInitializations(FemMatrixFreeClass_Base<Real_wp>& Mh0, FemMatrixFreeClass_Base<Real_wp>& Kh0)
{
}

template<>
void QnmSolver<LaplaceEquation<Dimension2> >::PerformInitializations(FemMatrixFreeClass_Base<Real_wp>& Mh0, FemMatrixFreeClass_Base<Real_wp>& Kh0)
{
}

template<>
void QnmSolver<HelmholtzEquation<Dimension2> >
::PerformInitializations(FemMatrixFreeClass_Base<Real_wp>& Mh0, FemMatrixFreeClass_Base<Real_wp>& Kh0)
{
}

template<>
void QnmSolver<StaticMaxwellEquation_3D>
::ComputeMassMatrix()
{
  Mh_sparse = vars.GetNewIterativeMatrix(Real_wp(0));
  
  GlobalGenericMatrix<Real_wp> nat_mat0;
  nat_mat0.SetCoefMass(1.0);
  nat_mat0.SetCoefStiffness(0.0);
  nat_mat0.SetCoefDamping(0.0);
  vars.SetCoefficientDirichlet(Real_wp(1));
  vars.AddMatrixWithBC(*Mh_sparse, nat_mat0);
  
  vars.SetCoefficientDirichlet(Real_wp(0));
}

template<>
void QnmSolver<HarmonicMaxwellEquation_3D>
::ComputeMassMatrix()
{
  Mh_sparse = vars.GetNewIterativeMatrix(Complex_wp(0, 0));
  
  GlobalGenericMatrix<Complex_wp> nat_mat0;
  nat_mat0.SetCoefMass(1.0);
  nat_mat0.SetCoefStiffness(0.0);
  nat_mat0.SetCoefDamping(0.0);
  vars.SetCoefficientDirichlet(Real_wp(1));
  vars.AddMatrixWithBC(*Mh_sparse, nat_mat0);
  
  vars.SetCoefficientDirichlet(Real_wp(0));

  Sh_sparse = vars.GetNewIterativeMatrix(Complex_wp(0, 0));

  nat_mat0.SetCoefMass(0.0);
  nat_mat0.SetCoefStiffness(0.0);
  nat_mat0.SetCoefDamping(1.0);
  vars.SetCoefficientDirichlet(Real_wp(0));
  vars.AddMatrixWithBC(*Sh_sparse, nat_mat0);
}

// on traite la matrice de masse comme une matrice diagonale
template<>
void QnmSolver<LaplaceEquation<Dimension2> >::ComputeMassMatrix()
{
  diag_Mh.Reallocate(vars.GetNbDof());
  diag_Mh.Zero();

  extra_diag_Mh.Reallocate(vars.GetNbDof());
  extra_diag_Mh.Zero();

  num_extra_Mh.Reallocate(vars.GetNbDof());
  num_extra_Mh.Fill(-1);

  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  int N = mesh_num.GetNbDof();
  int offset_v = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
  
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      int offset_Q = vars.GetNbDof() - vars.GetNbVectorialDofDrude();
      int offset_P = offset_Q - vars.GetNbVectorialDofDrude();
      offset_P += vars.OffsetDofDrudeV(i);
      offset_Q += vars.OffsetDofDrudeV(i);

      int ref = vars.mesh.Element(i).GetReference();
      const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
      bool affine = vars.mesh.IsElementAffine(i);
      int offset = offset_v + vars.GetOffsetDofV(i);
      
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          Real_wp jacob = vars.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
          int num_dof = mesh_num.Element(i).GetNumberDof(j);
          diag_Mh(num_dof) += vars.ref_rho(ref).GetConstant() * jacob;
          if (vars.InsidePML(i))
            {
              int num_dof_pml = mesh_num.GetDofPML(num_dof);
              diag_Mh(N+num_dof_pml) += vars.ref_rho(ref).GetConstant() * jacob;
            }
          
          diag_Mh(offset + 2*j) += vars.ref_invMu(ref).GetConstant()(0, 0)*jacob;
          diag_Mh(offset + 2*j+1) += vars.ref_invMu(ref).GetConstant()(1, 1)*jacob;
          
          if (vars.ref_drude(ref).IsEnabled())
            {
              int nPole = vars.ref_drude(ref).gamma.GetM();
              for (int kp = 0; kp < nPole; kp++)
                {
                  Real_wp sig = vars.ref_drude(ref).eps_sigma(kp);
                  
                  if (vars.ref_drude(ref).IsModeTE())
                    {
                      int dof_Px = offset_P + 2*j*nPole + kp, dof_Py = dof_Px + nPole;
                      
                      diag_Mh(dof_Px) = jacob;
                      diag_Mh(dof_Py) = jacob;
                      
                      int dof_Qx = offset_Q + 2*j*nPole + kp, dof_Qy = dof_Qx + nPole;
                      
                      diag_Mh(dof_Qx) = jacob;
                      diag_Mh(dof_Qy) = jacob;

                      extra_diag_Mh(dof_Qx) = -sig*jacob;
                      extra_diag_Mh(dof_Qy) = -sig*jacob;

                      num_extra_Mh(dof_Qx) = offset + 2*j;
                      num_extra_Mh(dof_Qy) = offset + 2*j+1;
                    }
                  else
                    {
                      int num_dof_P = offset_P + j*nPole + kp;
                      int num_dof_Q = offset_Q + j*nPole + kp;
                      diag_Mh(num_dof_P) += jacob;
                      diag_Mh(num_dof_Q) += jacob;
                      extra_diag_Mh(num_dof_Q) -= sig*jacob;
                      num_extra_Mh(num_dof_Q) = num_dof;
                    }
                }      
            }    
        }
    }
  
  //diag_Mh.WriteText("MhDiag.dat");
  //extra_diag_Mh.WriteText("ExtraMhDiag.dat");
  //num_extra_Mh.WriteText("NumMh.dat");
}

template<>
void QnmSolver<HelmholtzEquation<Dimension2> >::ComputeMassMatrix()
{
  Mh_sparse = vars.GetNewIterativeMatrix(Complex_wp(0, 0));
  
  GlobalGenericMatrix<Complex_wp> nat_mat0;
  nat_mat0.SetCoefMass(1.0);
  nat_mat0.SetCoefStiffness(0.0);
  nat_mat0.SetCoefDamping(0.0);
  vars.SetCoefficientDirichlet(Real_wp(1));
  vars.AddMatrixWithBC(*Mh_sparse, nat_mat0);

  Sh_sparse = vars.GetNewIterativeMatrix(Complex_wp(0, 0));

  nat_mat0.SetCoefMass(0.0);
  nat_mat0.SetCoefStiffness(0.0);
  nat_mat0.SetCoefDamping(1.0);
  vars.SetCoefficientDirichlet(Real_wp(0));
  vars.AddMatrixWithBC(*Sh_sparse, nat_mat0);

}

template<>
void QnmSolver<StaticMaxwellEquation_3D>
::ComputePhysicalMass()
{
  const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
  sparse_Mh.Reallocate(mesh_num.GetNbDof(), mesh_num.GetNbDof());
  sparse_Mh.Zero();
  
  Matrix<Real_wp> mat_elem;
  Real_wp threshold = 1e-14;
  Vector<int> col; VectReal_wp val;
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
      bool affine = vars.mesh.IsElementAffine(i);
      mat_elem.Reallocate(Fb.GetNbDof(), Fb.GetNbDof());
      mat_elem.Zero();
      
      Real_wp jacob, jacob_weighted; Matrix3_3 dfjm1, Amass;
      int N = Fb.GetNbPointsQuadratureInside();
      IVect Nodle = mesh_num.number_map.GetDofNumberOnElement(mesh_num, i);
      
      if (!vars.InsidePML(i))
        {     
          if (affine)
            N = 1;
          
          Vector<Matrix3_3> C(N);
          for (int j = 0; j < N; j++)
            {
              vars.GetInverseJacobianMatrix(i, j, affine, Fb.GetGeometricElement(),
                                            dfjm1, jacob, jacob_weighted);
                
              MltTrans(dfjm1, dfjm1, Amass);
              if (affine)
                Amass *= jacob;
              else
                Amass *= jacob_weighted;
              
              C(j) = Amass;
            }
          
          if (affine)
            Fb.AddConstantMassMatrix(0, 0, C(0), mat_elem);
          else
            Fb.AddVariableMassMatrix(0, 0, C, mat_elem);
          
          mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat_elem, i, 1);
          mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat_elem, i, 1);
          
          Vector<int> permut(Fb.GetNbDof());
          permut.Fill();
          Sort(Nodle, permut);
          
          col.Reallocate(Fb.GetNbDof());
          val.Reallocate(Fb.GetNbDof()); int nb;
          for (int j = 0; j < N; j++)
            {
              nb = 0;
              for (int k = 0; k < N; k++)
                if (abs(mat_elem(permut(j), permut(k))) > threshold)
                  {
                    col(nb) = Nodle(k);
                    val(nb) = mat_elem(permut(j), permut(k));
                    nb++;
                  }
              
              sparse_Mh.AddInteractionRow(Nodle(j), nb, col, val, true);
            }
        }
    }
}

template<>
void QnmSolver<HarmonicMaxwellEquation_3D>
::ComputePhysicalMass()
{
  const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
  sparse_Mh.Reallocate(mesh_num.GetNbDof(), mesh_num.GetNbDof());
  sparse_Mh.Zero();
  
  Matrix<Real_wp> mat_elem;
  Real_wp threshold = 1e-14;
  Vector<int> col; VectReal_wp val;
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
      bool affine = vars.mesh.IsElementAffine(i);
      mat_elem.Reallocate(Fb.GetNbDof(), Fb.GetNbDof());
      mat_elem.Zero();
      
      Real_wp jacob, jacob_weighted; Matrix3_3 dfjm1, Amass;
      int N = Fb.GetNbPointsQuadratureInside();
      IVect Nodle = mesh_num.number_map.GetDofNumberOnElement(mesh_num, i);
      
      if (!vars.InsidePML(i))
        {     
          if (affine)
            N = 1;
          
          Vector<Matrix3_3> C(N);
          for (int j = 0; j < N; j++)
            {
              vars.GetInverseJacobianMatrix(i, j, affine, Fb.GetGeometricElement(),
                                            dfjm1, jacob, jacob_weighted);
                
              MltTrans(dfjm1, dfjm1, Amass);
              if (affine)
                Amass *= jacob;
              else
                Amass *= jacob_weighted;
              
              C(j) = Amass;
            }
          
          if (affine)
            Fb.AddConstantMassMatrix(0, 0, C(0), mat_elem);
          else
            Fb.AddVariableMassMatrix(0, 0, C, mat_elem);
          
          mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat_elem, i, 1);
          mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat_elem, i, 1);
          
          Vector<int> permut(Fb.GetNbDof());
          permut.Fill();
          Sort(Nodle, permut);
          
          col.Reallocate(Fb.GetNbDof());
          val.Reallocate(Fb.GetNbDof()); int nb;
          for (int j = 0; j < N; j++)
            {
              nb = 0;
              for (int k = 0; k < N; k++)
                if (abs(mat_elem(permut(j), permut(k))) > threshold)
                  {
                    col(nb) = Nodle(k);
                    val(nb) = mat_elem(permut(j), permut(k));
                    nb++;
                  }
              
              sparse_Mh.AddInteractionRow(Nodle(j), nb, col, val, true);
            }
        }
    }
}

template<>
void QnmSolver<LaplaceEquation<Dimension2> >:: ComputePhysicalMass()
{
  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  mass_Mh.Reallocate(mesh_num.GetNbDof());
  mass_Mh.Zero();
  
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
      bool affine = vars.mesh.IsElementAffine(i);
      int ref = vars.mesh.Element(i).GetReference();
      
      if (!vars.InsidePML(i))
        {
          if ((input_var.reference_error == -1) || (input_var.reference_error == ref))
            for (int j = 0; j < Fb.GetNbDof(); j++)
              {
                Real_wp jacob = vars.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
                int num_dof = mesh_num.Element(i).GetNumberDof(j);
                mass_Mh(num_dof) += jacob;
              }
        }
    }
}

template<>
void QnmSolver<HelmholtzEquation<Dimension2> >:: ComputePhysicalMass()
{
  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  mass_Mh.Reallocate(mesh_num.GetNbDof());
  mass_Mh.Zero();
  
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      const ElementReference_Dim<Dimension>& Fb = vars.GetReferenceElement(i);
      bool affine = vars.mesh.IsElementAffine(i);
      int ref = vars.mesh.Element(i).GetReference();
      
      if (!vars.InsidePML(i))
        {
          if ((input_var.reference_error == -1) || (input_var.reference_error == ref))
            for (int j = 0; j < Fb.GetNbDof(); j++)
              {
                Real_wp jacob = vars.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
                int num_dof = mesh_num.Element(i).GetNumberDof(j);
                mass_Mh(num_dof) += jacob;
              }
        }
    }
}

template<>
Complex_wp  QnmSolver<LaplaceEquation<Dimension2> >::GetRatioDeltaEpsilon(const Real_wp& om, const Complex_wp& om_j, int k)
{
  int ref = ref_cavity(k);
  if (vars.ref_drude(ref).IsEnabled())
    return vars.ref_drude(ref).GetDeltaEpsilon(om_j) / vars.ref_drude(ref).GetDeltaEpsilon(om);
  
  return Complex_wp(1, 0);
}

template<>
Complex_wp  QnmSolver<HelmholtzEquation<Dimension2> >::GetRatioDeltaEpsilon(const Real_wp& om, const Complex_wp& om_j, int k)
{
  int ref = ref_cavity(k);
  if (vars.ref_drude(ref).IsEnabled())
    return vars.ref_drude(ref).GetDeltaEpsilon(om_j) / vars.ref_drude(ref).GetDeltaEpsilon(om);
  
  return Complex_wp(1, 0);
}

template<>
Complex_wp  QnmSolver<StaticMaxwellEquation_3D>::GetRatioDeltaEpsilon(const Real_wp& om, const Complex_wp& om_j, int k)
{
  int ref = ref_cavity(k);
  if (vars.ref_drude(ref).IsEnabled())
    return vars.ref_drude(ref).GetDeltaEpsilon(om_j) / vars.ref_drude(ref).GetDeltaEpsilon(om);
  
  return Complex_wp(1, 0);
}

template<>
Complex_wp  QnmSolver<HarmonicMaxwellEquation_3D>::GetRatioDeltaEpsilon(const Real_wp& om, const Complex_wp& om_j, int k)
{
  int ref = ref_cavity(k);
  if (vars.ref_drude(ref).IsEnabled())
    return vars.ref_drude(ref).GetDeltaEpsilon(om_j) / vars.ref_drude(ref).GetDeltaEpsilon(om);
  
  return Complex_wp(1, 0);
}

template<>
void QnmSolver<StaticMaxwellEquation_3D>::AddIncidentWave(const Real_wp& alpha, const Real_wp& omega, VectComplex_wp& X)
{
  const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
  R3 kwave = vars.GetWaveVector(), polar;
  Mlt(omega/Norm2(kwave), kwave);
  vars.GetPolarization(polar);
  
  VectR3 s; SetPoints<Dimension3> PointsElem;
  SetMatrices<Dimension3> MatricesElem;
  Vector<VectComplex_wp> feval(1), contrib(1);
  Vector<bool> DofUsed(mesh_num.GetNbDof());
  DofUsed.Fill(false);
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      vars.mesh.GetVerticesElement(i, s);
      const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
      Fb.FjElemDof(s, PointsElem, vars.mesh, i);
      Fb.DFjElemDof(s, PointsElem, MatricesElem, vars.mesh, i);
      
      feval(0).Reallocate(3*Fb.GetNbPointsDofInside()); feval(0).Zero();
      contrib(0).Reallocate(Fb.GetNbDof()); contrib(0).Zero();
      for (int j = 0; j < Fb.GetNbPointsDofInside(); j++)
        {
          R3 pt = PointsElem.GetPointDof(j);
          Complex_wp u_inc = exp(Iwp*DotProd(kwave, pt));
          feval(0)(3*j) = u_inc*polar(0);
          feval(0)(3*j+1) = u_inc*polar(1);
          feval(0)(3*j+2) = u_inc*polar(2);
        }
      
      Fb.ComputeProjectionDof(MatricesElem, feval, contrib, mesh_num, i);
      
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          int num_dof = mesh_num.Element(i).GetNumberDof(j);
          if (!DofUsed(num_dof))
            {
              DofUsed(num_dof) = true;
              X(num_dof) += alpha*contrib(0)(j);
            }
        }          
    }
}

template<>
void QnmSolver<HarmonicMaxwellEquation_3D>::AddIncidentWave(const Real_wp& alpha, const Real_wp& omega, VectComplex_wp& X)
{
  const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
  R3 kwave = vars.GetWaveVector(), polar;
  Mlt(omega/Norm2(kwave), kwave);
  vars.GetPolarization(polar);
  
  VectR3 s; SetPoints<Dimension3> PointsElem;
  SetMatrices<Dimension3> MatricesElem;
  Vector<VectComplex_wp> feval(1), contrib(1);
  Vector<bool> DofUsed(mesh_num.GetNbDof());
  DofUsed.Fill(false);
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      vars.mesh.GetVerticesElement(i, s);
      const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
      Fb.FjElemDof(s, PointsElem, vars.mesh, i);
      Fb.DFjElemDof(s, PointsElem, MatricesElem, vars.mesh, i);
      
      feval(0).Reallocate(3*Fb.GetNbPointsDofInside()); feval(0).Zero();
      contrib(0).Reallocate(Fb.GetNbDof()); contrib(0).Zero();
      for (int j = 0; j < Fb.GetNbPointsDofInside(); j++)
        {
          R3 pt = PointsElem.GetPointDof(j);
          Complex_wp u_inc = exp(Iwp*DotProd(kwave, pt));
          feval(0)(3*j) = u_inc*polar(0);
          feval(0)(3*j+1) = u_inc*polar(1);
          feval(0)(3*j+2) = u_inc*polar(2);
        }
      
      Fb.ComputeProjectionDof(MatricesElem, feval, contrib, mesh_num, i);
      
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          int num_dof = mesh_num.Element(i).GetNumberDof(j);
          if (!DofUsed(num_dof))
            {
              DofUsed(num_dof) = true;
              X(num_dof) += alpha*contrib(0)(j);
            }
        }          
    }
}

template<>
void QnmSolver<LaplaceEquation<Dimension2> >::AddIncidentWave(const Real_wp& alpha, const Real_wp& omega, VectComplex_wp& X)
{
  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  R2 kwave = vars.GetWaveVector();
  Mlt(omega/Norm2(kwave), kwave);
  
  VectR2 s; SetPoints<Dimension2> PointsElem;
  Vector<bool> DofUsed(mesh_num.GetNbDof());
  DofUsed.Fill(false);
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      vars.mesh.GetVerticesElement(i, s);
      const ElementReference<Dimension2, 1>& Fb = vars.GetReferenceElementH1(i);
      Fb.FjElemDof(s, PointsElem, vars.mesh, i);
      
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          int num_dof = mesh_num.Element(i).GetNumberDof(j);
          if (!DofUsed(num_dof))
            {
              DofUsed(num_dof) = true;
              R2 pt = PointsElem.GetPointDof(j);
              Complex_wp u_inc = exp(Iwp*DotProd(kwave, pt));
              
              X(num_dof) += alpha*u_inc;
            }
        }          
    }
}


template<>
void QnmSolver<HelmholtzEquation<Dimension2> >::AddIncidentWave(const Real_wp& alpha, const Real_wp& omega, VectComplex_wp& X)
{
  R2 kwave = vars.GetWaveVector();
  Mlt(omega/Norm2(kwave), kwave);

  const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
  VectR2 s; SetPoints<Dimension2> PointsElem;
  Vector<bool> DofUsed(mesh_num.GetNbDof());
  DofUsed.Fill(false);
  for (int i = 0; i < vars.mesh.GetNbElt(); i++)
    {
      vars.mesh.GetVerticesElement(i, s);
      const ElementReference<Dimension2, 1>& Fb = vars.GetReferenceElementH1(i);
      Fb.FjElemDof(s, PointsElem, vars.mesh, i);
      
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          int num_dof = mesh_num.Element(i).GetNumberDof(j);
          if (!DofUsed(num_dof))
            {
              DofUsed(num_dof) = true;
              R2 pt = PointsElem.GetPointDof(j);
              Complex_wp u_inc = exp(Iwp*DotProd(kwave, pt));
              
              X(num_dof) += alpha*u_inc;
            }
        }          
    }
}

int main(int argc, char **argv) 
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
      
      int rank_proc(0);
      Vector<string> lines_data_file;
      
      // lecture du fichier de donnees
#ifdef SELDON_WITH_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
      ReadLinesFile(file_name_data, lines_data_file, MPI_COMM_WORLD);
#else
      ReadLinesFile(file_name_data, lines_data_file);
#endif
      
      InputVariables input_var;
      ReadInputFile(lines_data_file, input_var);
      
      if (rank_proc == 0)
        cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;
      
      if (type_equation == "MAXWELL")
        {
          if (input_var.dispersive_pml)
            {
              QnmSolver<StaticMaxwellEquation_3D> solver_qnm(input_var);
              solver_qnm.RunAll(lines_data_file, type_element, type_equation);
              
              if (rank_proc == 0)
                cout<<" we destroy the variables "<<endl;
            }
          else
            {
              QnmSolver<HarmonicMaxwellEquation_3D> solver_qnm(input_var);
              solver_qnm.RunAll(lines_data_file, type_element, type_equation);
              
              if (rank_proc == 0)
                cout<<" we destroy the variables "<<endl;
            }
        }
      else
        {
          if (input_var.dispersive_pml)
            {
              QnmSolver<LaplaceEquation<Dimension2> > solver_qnm(input_var);
              solver_qnm.RunAll(lines_data_file, type_element, type_equation);
              
              if (rank_proc == 0)
                cout<<" we destroy the variables "<<endl;
              
            }
          else
            {
              QnmSolver<HelmholtzEquation<Dimension2> > solver_qnm(input_var);
              solver_qnm.RunAll(lines_data_file, type_element, type_equation);
              
              if (rank_proc == 0)
                cout<<" we destroy the variables "<<endl;
            }
        }
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"mode1D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
  
}
