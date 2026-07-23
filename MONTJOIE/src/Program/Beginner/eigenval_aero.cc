#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_NODAL_HCURL

#include "Elliptic/Helmholtz/MontjoieLaplace.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell3D.hxx"

#include "Elliptic/Aeroacoustic/MontjoieAeroacoustic.hxx"

using namespace Montjoie;

class InputVariables : public InputDataProblem_Base
{
  
public :
  // parametres pour calculer les valeurs propres (shifts et nombre de vps pour chaque shift)
  Vector<int> number_eigenval;
  Vector<Complex_wp> shift_cplx;
  VectReal_wp radius_shift;
  Real_wp tolerance_eigenvalue;
  int nb_max_iterations;
  
  // spectre a selectionner
  Real_wp Lr_max, Li_max, Li_min, spectrum_radius_min;
  
  // solveur a choisir
  int type_eigensolver;
  
  // fichier contenant les valeurs propres
  string file_name_eigen;

  // stocke-t-on les matrices Mh et Kh ?
  bool explicit_matrix;

  // Threshold_EigenValue
  Real_wp threshold_eigen, threshold_eigen_vector, threshold_eigenval_equal;
  
  // fichier ou on ecrit les vecteurs propres
  string file_name_eigenvectors;
  
  // affiche-t-on les vecteurs propres
  bool display_eigenvectors;
  
  // rang du processeur
  int rank_proc;

  // mode a choisir
  int computational_mode;  

  InputVariables()
  {    
    number_eigenval.Reallocate(1);
    number_eigenval(0) = 100;
    tolerance_eigenvalue = 1e-12;
    nb_max_iterations = 1000;

    Lr_max = Real_wp(0.6); Li_max = Real_wp(30); Li_min = Real_wp(0);
    spectrum_radius_min = Real_wp(0);

    type_eigensolver = TypeEigenvalueSolver::ARPACK;
    
    file_name_eigen = string("Lambda");
    
    explicit_matrix = false;

    threshold_eigen = 1e-6;
    threshold_eigenval_equal = 1e-6;
    
    R2::threshold = 1e-6;
    threshold_eigen_vector = 1e-10;

    display_eigenvectors = false;

    rank_proc = 0;
    computational_mode = -1;
  }
  
  bool SelectEigenvalue(const Complex_wp& z)
  {
    if ((realpart(z) <= Lr_max)
	&& (abs(imagpart(z)) <= Li_max)
	&& (abs(imagpart(z)) >= Li_min)
	&& (abs(z) >= spectrum_radius_min))
      {
        return true;
      }
    
    return false;
  }

  bool SelectEigenvector(const VectComplex_wp& x, int N)
  {
    return true;
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
    else if (keyword == "SpectrumRadius")
      spectrum_radius_min = to_num<Real_wp>(param(0));
    else if (keyword == "EigenvalueTolerance")
      tolerance_eigenvalue = to_num<Real_wp>(param(0));
    else if (keyword == "EigenvalueMaxNumberIterations")
      nb_max_iterations = to_num<int>(param(0));
    else if (keyword == "ComputationalMode")
      {
        if (param(0) == "SHIFTED")
          computational_mode = EigenProblem_Base<Real_wp>::SHIFTED_MODE;
        else if (param(0) == "REGULAR")
          computational_mode = EigenProblem_Base<Real_wp>::REGULAR_MODE;
        else if (param(0) == "COMPLEX_SHIFT")
          computational_mode = EigenProblem_Base<Real_wp>::IMAG_SHIFTED_MODE;
        else if (param(0) == "BUCKLING")
          computational_mode = EigenProblem_Base<Real_wp>::BUCKLING_MODE;
        else if (param(0) == "CAYLEY")
          computational_mode = EigenProblem_Base<Real_wp>::CAYLEY_MODE;
        else if (param(0) == "INVERT")
          computational_mode = EigenProblem_Base<Real_wp>::INVERT_MODE;
        else
          {
            cout << "Unknown computational mode" << endl;
            abort();
          }
      }
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
    else if (keyword == "RadiusShift")
      {
        radius_shift.Reallocate(param.GetM());
        for (int i = 0; i < param.GetM(); i++)
          radius_shift(i) = to_num<Real_wp>(param(i));
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
      }
    else if (keyword == "FileEigenvalue")
      file_name_eigen = param(0);
    else if (keyword == "ExplicitMatrixEigen")
      {
	if (param(0) == "YES")
	  explicit_matrix = true;
	else
	  explicit_matrix = false;
      }
    else if (keyword == "ThresholdEigenValue")
      {
	threshold_eigen = to_num<Real_wp>(param(0));
      }
    else if (keyword == "ThresholdEigenvalueEqual")
      threshold_eigenval_equal = to_num<Real_wp>(param(0));
    else if (keyword == "DisplayEigenvectors")
      {
        if (param(0) == "YES")
          {
            display_eigenvectors = true;
            file_name_eigenvectors = param(1);
          }
        else
          display_eigenvectors = false;
      }
  }

};


void GetComplexEigenvalues(const VectReal_wp& Lr, const VectReal_wp& Li, VectComplex_wp& L)
{
  int nev = Lr.GetM();
  L.Reallocate(nev);
  for (int i = 0; i < nev; i++)
    L(i) = Complex_wp(Lr(i), Li(i));
}

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
			  const Vector<int>& imag_eigenvec,
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
  Vector<Complex_wp> x(new_eigenvec.GetM()), y(new_eigenvec.GetM());
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

void MergeEigenvalues(const VectComplex_wp& new_lambda, const VectComplex_wp& new_lambda_imag,
		      Matrix<Complex_wp>& new_eigenvec, VectComplex_wp& lambda,
		      VectComplex_wp& lambda_imag, Matrix<Complex_wp>& eigen_vec)
{
  lambda = new_lambda; lambda_imag = new_lambda_imag; eigen_vec = new_eigenvec;
  //cout << "Not implemented" << endl;
  //abort();
}


void MergeEigenvalues(const VectReal_wp& new_lambda, const VectReal_wp& new_lambda_imag,
		      Matrix<Real_wp>& new_eigenvec, VectReal_wp& lambda,
		      VectReal_wp& lambda_imag, Matrix<Real_wp>& eigen_vec)
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
  FindEigenvalueNumber(Lambda, NewLambda, imag_eigenvec, true, num_mode, new_eigenvec);
  
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

void UpdateEigenValues(const VirtualMatrix<Real_wp>& mat_mass, const VirtualMatrix<Real_wp>& mat_stiff,
		       VectReal_wp& lambda0_real, VectReal_wp& lambda0_imag, Matrix<Real_wp>& eigen_vec0,
		       Real_wp threshold)
{
  int N = mat_stiff.GetM();
  Vector<Real_wp> lambda(lambda0_real),lambda_imag(lambda0_imag);
  Matrix<Real_wp> eigen_vec(eigen_vec0);
  
  Vector<Real_wp> X(N), Xi(N), Y(N), Yi(N), Mx(N), Mxi(N);
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
          
          //X.WriteText("X.dat"); Xi.WriteText("Xi.dat");
	  mat_stiff.MltVector(X, Y); 
          mat_stiff.MltVector(Xi, Yi);
          mat_mass.MltVector(X, Mx);
          mat_mass.MltVector(Xi, Mxi);
          normeX = sqrt(DotProd(Mx, Mx) + DotProd(Mxi, Mxi));
          for (int j = 0; j < N; j++)
            err += pow(abs(Complex_wp(Y(j), Yi(j))
                           - Complex_wp(lambda(i), lambda_imag(i))
                           *Complex_wp(Mx(j), Mxi(j)) ), 2.0);
          
          err = sqrt(err);
        }
      else
        {
          for (int j = 0; j < N; j++)
            X(j) = eigen_vec(j, i);
          
          mat_stiff.MltVector(X, Y);
	  mat_mass.MltVector(X, Mx);
          normeX = sqrt(DotProd(Mx, Mx));
          for (int j = 0; j < N; j++)
            err += pow(Y(j) - lambda(i)*Mx(j), 2.0);
          
          err = sqrt(err);
        }
      
      if (err > threshold*normeX)
        {
          cout << "Error on eigenvalue " << lambda(i) << " " << lambda_imag(i) << endl;
          cout << "Error = " << err/normeX << endl;
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

  cout << "Error on accepted eigenvectors = " << err_inf << endl;
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
          solver.ExpandVector(col0, col);
          
          SetCol(col, j, eigen_vec);
        }
      
    }
  catch(...)
    {
    }
}

template<class T>
void ComputeEigenvaluesShift(VirtualEigenProblem<T, T, T>& eigen_solver,
			     const VirtualMatrix<T>& M, const VirtualMatrix<T>& K, int N,
			     InputVariables& data, Vector<T>& lambda, Vector<T>& lambda_imag,
			     Matrix<T>& eigen_vec, string file_eigenvalue, bool symM, bool second_order)
{
  Vector<T> new_lambda, new_lambda_imag;
  Matrix<T> new_eigenvec;

  lambda.Clear(); lambda_imag.Clear();
  eigen_vec.Clear();

  int type_solver = data.type_eigensolver;
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
  
  for (int p = 0; p < data.shift_cplx.GetM(); p++)
    {
      if (symM)
        eigen_solver.SetComputationalMode(eigen_solver.SHIFTED_MODE);
      else
        eigen_solver.SetComputationalMode(eigen_solver.INVERT_MODE);
      
      if (type_solver == TypeEigenvalueSolver::SLEPC)
        eigen_solver.SetComputationalMode(eigen_solver.INVERT_MODE);
      
      if (data.computational_mode >= 0)
        eigen_solver.SetComputationalMode(data.computational_mode);
      
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
          if (second_order)
            {
              if (data.rank_proc == 0)
                {
                  DISP(square(data.shift_cplx(p))); DISP(data.number_eigenval(p));              
                }
              
              eigen_solver.SetTypeSpectrum(eigen_solver.CENTERED_EIGENVALUES,
                                           square(data.shift_cplx(p)),
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
        }
      
      if (type_solver == TypeEigenvalueSolver::FEAST)
        {
          if (p >= data.radius_shift.GetM())
            {
              cout << "Missing Radius (RadiusShift) " << endl;
              DISP(p); DISP(data.radius_shift.GetM());
            }
          
          //eigen_solver.SetCircleSpectrum(Iwp*data.shift_cplx(p), data.radius_shift(p));
          //eigen_solver.SetGlobalPrintLevel(1);
          eigen_solver.SetPrintLevel(1);
          if (data.rank_proc != 0)
            eigen_solver.SetPrintLevel(-1);
        }
      else
        eigen_solver.SetPrintLevel(3);
      
      eigen_solver.Init(N);
      
      // calcul des valeurs propres
      GetEigenvaluesEigenvectors(eigen_solver, new_lambda, new_lambda_imag, new_eigenvec, type_solver);
      if (data.rank_proc != 0)
        return;
      
      RetrieveEigenvectors(eigen_solver, K, new_eigenvec);
      //UpdateEigenValues(M, K, new_lambda, new_lambda_imag, new_eigenvec, data.threshold_eigen);
      
      if (data.rank_proc == 0)
        {
          new_lambda.WriteText(file_eigenvalue + to_str(p) + "_real.dat");
          new_lambda_imag.WriteText(file_eigenvalue + to_str(p) + "_imag.dat");
        }
      
      MergeEigenvalues(new_lambda, new_lambda_imag, new_eigenvec,
                       lambda, lambda_imag, eigen_vec);
    }
}

void ExtractEigenvectors(VarProblem_Base& vars,
                         InputVariables& input_var,
                         const VectReal_wp& lambda, const VectReal_wp& lambda_imag,
                         const Matrix<Real_wp>& eigen_vec,
                         VectComplex_wp& AllLambda, Vector<VectComplex_wp>& AllEigenVec)
{
  int nev = lambda.GetM();
  //int N = vars.GetNbDof();
  
  // we count the number of eigenvalues that are kept
  int num = 0, nb = 0;
  int Nu = vars.GetNbMeshDof();
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
  
  DISP(num);
  AllLambda.Reallocate(num);
  AllEigenVec.Reallocate(num);
  
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


void ExtractEigenvectors(VarProblem_Base& vars,
                         InputVariables& input_var,
                         const VectComplex_wp& lambda, const VectComplex_wp& lambda_imag,
                         const Matrix<Complex_wp>& eigen_vec,
                         VectComplex_wp& AllLambda, Vector<VectComplex_wp>& AllEigenVec)
{
  int nev = lambda.GetM();
  //int N = vars.GetNbDof();
  
  int Nu = vars.GetNbMeshDof();
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
  num = 0;
  for (int i = 0; i < nev; i++)
    if (input_var.SelectEigenvalue(lambda(i)))
      {
        GetCol(eigen_vec, i, x);
        if (input_var.SelectEigenvector(x, Nu))
          {
            AllLambda(num) = lambda(i);
            AllEigenVec(num) = x;
            num++;
          }
      }
}

template<class T>
void WriteSolution(VarOutputProblem_Base& vars,
                   const Vector<T>& u, const string& root)
{
  Vector<Vector<T> > eigen_mode_v;
  eigen_mode_v.SetData(1, &const_cast<Vector<T>& >(u));
  
  // on change les noms de fichier
  for (int j = 0; j < vars.output_mesh_param.GetM(); j++)
    vars.output_mesh_param(j).
      SetFileName(2, root + "G" + to_str(j) + string(".bb"));
  
  for (int j = 0; j < vars.output_grid_param.GetM(); j++)
    vars.output_grid_param(j).
      SetFileName(2, root + "G" + to_str(j) + string(".dat"));
  
  vars.WriteOutputFile(eigen_mode_v, 2);
  eigen_mode_v.Nullify();
}

void WriteEigenvectors(VarOutputProblem_Base& var, 
                       InputVariables& input_var,
                       Vector<VectComplex_wp>& eigen_vec)
{
  for (int j = 0; j < eigen_vec.GetM(); j++)
    {
      string name = input_var.file_name_eigenvectors + "_N" + to_str(j) + "_";      
      WriteSolution(var, eigen_vec(j), name);
    }                  
}


void WriteRealEigenvectors(VarOutputProblem_Base& var, 
                           InputVariables& input_var,
                           Matrix<Real_wp>& eigen_vec)
{
  VectReal_wp u(eigen_vec.GetM());
  for (int j = 0; j < eigen_vec.GetN(); j++)
    {
      string name = input_var.file_name_eigenvectors + "_N" + to_str(j) + "_";      
      GetCol(eigen_vec, j, u);
      
      WriteSolution(var, u, name);
    }
}

// fonction principale
template<class Complexe, class Dimension>
void RunAll(VarHarmonic_Base<Complexe, Dimension>& var,
            InputVariables& input_var, const string& name_file,
            const string& name_element, const string& name_equation)
{
  // initialisation indices physiques
#ifdef SELDON_WITH_MPI
  var.comm_group_mode = MPI_COMM_WORLD;
#endif

  var.InitIndices(PhysicalConstant::nb_max_indices);
  
  // lecture du fichier de donnees
  Vector<string> lines_data_file;
  
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
  ReadLinesFile(name_file, lines_data_file, var.comm_group_mode);
#else
  int rank_proc(0), nb_proc(1);
  ReadLinesFile(name_file, lines_data_file);
#endif

  ReadInputFile(lines_data_file, var);
  ReadInputFile(lines_data_file, input_var);
  
  // dans le cas du solveur Feast, chaque proc calcule la meme matrice
#ifdef SELDON_WITH_MPI
  if (input_var.type_eigensolver == TypeEigenvalueSolver::FEAST)
    {
      var.comm_group_mode = MPI_COMM_SELF;
      input_var.rank_proc = rank_proc;
    }
#endif
  
  if (rank_proc != 0)
    {
      var.print_level = -1;
      var.mesh.print_level = -1;
    }
  
  // construction maillage et elements finis
  MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
  mesh_num.SetSameNumberPeriodicDofs();
  var.ComputeMeshAndFiniteElement(name_element);
  var.PerformOtherInitializations();
  
  // quantites geometriques et autres ddls
  var.ComputeMassMatrix();
  var.ComputeQuasiPeriodicPhase();

  var.SetCurrentModeNumber(0);
  
  // (-omega^2 M_h - i omega S_h + K_h) X = 0 (si FirstOrderFormulation () = false)
  // (-i omega M_h + K_h) X = 0 (si FirstOrderFormulation () = true)
  if (var.FirstOrderFormulation())
    cout << "Formulation d'ordre 1 " << endl;
  else
    cout << "Formulation d'ordre 2 " << endl;
  
  // on stocke les matrices Mh, Kh et Sh
  var.SetHomogeneousDirichlet(true); // pour les noeuds en Dirichlet
  
  if (!input_var.explicit_matrix)
    {
      All_LinearSolver* glob_solver;
      glob_solver = var.GetNewLinearSolver();
      ReadInputFile(lines_data_file, *glob_solver);

      FemMatrixFreeClass_Base<Complexe>* Kh, *Kh0;
      FemMatrixFreeClass_Base<Complexe>* Mh, *Mh0;
      
      Kh = var.GetNewIterativeMatrix(Complexe(0));
      Mh = var.GetNewIterativeMatrix(Complexe(0));
      
      Kh0 = var.GetNewIterativeMatrix(Complexe(0));
      Mh0 = var.GetNewIterativeMatrix(Complexe(0));

      GlobalGenericMatrix<Complexe> nat_mat0;
      nat_mat0.SetCoefMass(1.0);
      nat_mat0.SetCoefStiffness(0.0);
      nat_mat0.SetCoefDamping(0.0);
      var.SetCoefficientDirichlet(Real_wp(1));
      var.AddMatrixWithBC(*Mh0, nat_mat0);
      
      nat_mat0.SetCoefMass(0.0);
      nat_mat0.SetCoefStiffness(1.0);
      nat_mat0.SetCoefDamping(1.0);
      var.SetCoefficientDirichlet(Real_wp(0));
      var.AddMatrixWithBC(*Kh0, nat_mat0);

      EigenProblemMontjoie<Complexe>* eigen_solver = NULL;
      eigen_solver = var.GetNewEigenSolver(*glob_solver, Complexe());
        
      eigen_solver->InitMatrix(*Kh, *Mh);
      
      int print_level = var.print_level;
      var.print_level = -1;
        
      Vector<Complexe> lambda, lambda_imag; Matrix<Complexe> eigen_vec;

      ComputeEigenvaluesShift(*eigen_solver, *Mh0, *Kh0, eigen_solver->GetM(), input_var,
                              lambda, lambda_imag, eigen_vec, input_var.file_name_eigen,
                              true, !var.FirstOrderFormulation());
      
      var.print_level = print_level;

      if (input_var.rank_proc != 0)
        return;
      
      VectComplex_wp AllLambda;
      Vector<VectComplex_wp> AllEigenVec;
      ExtractEigenvectors(var, input_var, lambda, lambda_imag, 
                          eigen_vec, AllLambda, AllEigenVec);
      
      AllLambda.Write(input_var.file_name_eigen + ".dat");
      if (input_var.display_eigenvectors)
        WriteEigenvectors(var, input_var, AllEigenVec);
      
      delete Mh; delete Kh; delete Kh0; delete Mh0;
      delete glob_solver;

      return;
    }
  
  DistributedMatrix<Real_wp, General, ArrayRowSparse> Kh, Sh;
  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> Mh;
  GlobalGenericMatrix<Real_wp> nat_mat;
  
  nat_mat.SetCoefMass(1.0);
  nat_mat.SetCoefDamping(0.0);
  nat_mat.SetCoefStiffness(0.0);
  var.SetCoefficientDirichlet(Real_wp(1));
  var.AddMatrixWithBC(Mh, nat_mat);
  
  nat_mat.SetCoefMass(0.0);
  var.SetCoefficientDirichlet(Real_wp(0));
  if (var.FirstOrderFormulation())
    {
      nat_mat.SetCoefDamping(1.0);
    }
  else
    {
      nat_mat.SetCoefDamping(1.0);
      var.AddMatrixWithBC(Sh, nat_mat);
      
      nat_mat.SetCoefDamping(0.0);
    }
  
  nat_mat.SetCoefStiffness(1.0);
  var.AddMatrixWithBC(Kh, nat_mat);  
  
  if (var.print_level > 4)
    {
      Mh.WriteText("Mh.dat");
      Sh.WriteText("Sh.dat");
      Kh.WriteText("Kh.dat");
    }
  
  DistributedMatrix<Real_wp, General, ArrayRowSparse>* mat_stiff;
  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>* mat_mass;
  
  // linearisation : on se ramene a un systeme aux valeurs propres lineaires
  // (-i omega M + K) X = 0 (dans le cas ordre 1)
  DistributedMatrix<Real_wp, General, ArrayRowSparse> K;
  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> M;
  
  int N = var.GetNbDof();
  if (var.FirstOrderFormulation())
    {
      mat_mass = &Mh;
      mat_stiff = &Kh;
    }
  else
    {
      mat_mass = &M;
      mat_stiff = &K;
      if (Sh.GetNonZeros() == 0)
	{
          // cas symetrique (ou K et M seraient symetriques)
          // (- omega^2 Mh + Kh) X = 0  => valeur propre = omega^2
          cout << "Attention on considere que K est symmetrique " << endl;
          
	  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> Ksym;	  
          VectReal_wp lambda, lambda_imag; Matrix<Real_wp> eigen_vec;
	  Copy(Kh, Ksym);
          if (var.print_level > 4)
            Ksym.WriteText("Ksym.dat");
          
          SparseEigenProblem<Real_wp, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>,
                             DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> > eigen_solver;
          
          eigen_solver.InitMatrix(Ksym, Mh);
          
          ComputeEigenvaluesShift(eigen_solver, Mh, Ksym, Ksym.GetM(),
                                  input_var, lambda, lambda_imag, eigen_vec, input_var.file_name_eigen,
                                  true, true);
                    
          
          lambda.WriteText(input_var.file_name_eigen + ".dat");
          WriteRealEigenvectors(var, input_var, eigen_vec);
          
          return;
	}
      
      // cas non-symetrique
      // on rajoute une inconnue intermediaire
      K.Reallocate(2*N, 2*N);
      M.Reallocate(2*N, 2*N);
      
      IVect col; VectReal_wp val;
      for (int i = 0; i < N; i++)
        {
          M.AddInteraction(i, i, 1.0);
          int nb = Mh.GetRowSize(i);
          col.Reallocate(nb);
          val.Reallocate(nb);
          for (int j = 0; j < Mh.GetRowSize(i); j++)
            {
              col(j) = Mh.Index(i, j) + N;
              val(j) = Mh.Value(i, j);
            }
          
          M.AddInteractionRow(i+N, nb, col, val);
          
          K.AddInteraction(i, i+N, 1.0);
          nb = Kh.GetRowSize(i) + Sh.GetRowSize(i);
          col.Reallocate(nb);
          val.Reallocate(nb);
          nb = 0;
          for (int j = 0; j < Kh.GetRowSize(i); j++)
            {
              col(nb) = Kh.Index(i, j);
              val(nb) = -Kh.Value(i, j);
              nb++;
            }

          for (int j = 0; j < Sh.GetRowSize(i); j++)
            {
              col(nb) = Sh.Index(i, j) + N;
              val(nb) = -Sh.Value(i, j);
              nb++;
            }
          
          K.AddInteractionRow(i+N, nb, col, val);
        }
      
      if (var.print_level > 7)
	{
	  K.WriteText("Stiff.dat");
	  M.WriteText("Mass.dat");
	}
    }


  VectReal_wp lambda, lambda_imag;
  if (input_var.type_eigensolver == -1)
    {
      // cas ou on veut toutes les valeurs propres
      if (nb_proc == 1)
        {
          // dense eigensolver
          Matrix<Real_wp> Ah(Kh.GetM(), Kh.GetM());
          Ah.Zero();
          for (int i = 0; i < Kh.GetM(); i++)
            {
              if (Mh.GetRowSize(i) != 1)
                {
                  cout << "Only diagonal mass matrices handled" << endl;
                  abort();
                }
              
              if (Mh.Index(i, 0) != i)
                {
                  cout << "Only diagonal mass matrices handled" << endl;
                  abort();
                }
              
              Real_wp invMi = Real_wp(1) / Mh.Value(i, 0);
              for (int j = 0; j < Kh.GetRowSize(i); j++)
                Ah(i, Kh.Index(i, j)) = Kh.Value(i, j)*invMi;          
            }
          
          GetEigenvalues(Ah, lambda, lambda_imag);
        }
      else
        {
#ifdef SELDON_WITH_SCALAPACK
          // distributed eigensolver
          DistributedMatrix<Real_wp, General, ColMajor> Ah;
          int n = Kh.GetM();
          int nb = 64;
          Ah.Init(global_blacs_handle, n, n, nb, nb);
          Ah.Zero();
          for (int i = 0; i < Kh.GetM(); i++)
            {
              if (Mh.GetRowSize(i) != 1)
                {
                  cout << "Only diagonal mass matrices handled" << endl;
                  abort();
                }
              
              if (Mh.Index(i, 0) != i)
                {
                  cout << "Only diagonal mass matrices handled" << endl;
                  abort();
                }
              
              Real_wp invMi = Real_wp(1) / Mh.Value(i, 0);
              for (int j = 0; j < Kh.GetRowSize(i); j++)
                Ah.SetGlobal(i, Kh.Index(i, j), Kh.Value(i, j)*invMi);   
            }
          
          GetEigenvalues(Ah, lambda, lambda_imag);
#else
          cout << "Recompile with Scalapack" << endl;
          abort();
#endif
        }
      
      VectComplex_wp AllLambda;
      GetComplexEigenvalues(lambda, lambda_imag, AllLambda);
      AllLambda.Write(input_var.file_name_eigen + ".dat");
      if (input_var.display_eigenvectors)
        {
          cout << "Not implemented" << endl;
          abort();
        }
    }
  else
    {      
      Matrix<Real_wp> eigen_vec;
      
      // sparse eigensolver
      SparseEigenProblem<Real_wp, DistributedMatrix<Real_wp, General, ArrayRowSparse>,
                         DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> > eigen_solver;
      
      eigen_solver.InitMatrix(*mat_stiff, *mat_mass);
      
      ComputeEigenvaluesShift(eigen_solver, *mat_mass, *mat_stiff, mat_stiff->GetM(),
                              input_var, lambda, lambda_imag, eigen_vec, input_var.file_name_eigen,
                              false, true);
      
      if (input_var.rank_proc != 0)
        return;
      
      VectComplex_wp AllLambda;
      Vector<VectComplex_wp> AllEigenVec;
      ExtractEigenvectors(var, input_var, lambda, lambda_imag, 
                          eigen_vec, AllLambda, AllEigenVec);
      
      
      AllLambda.Write(input_var.file_name_eigen + ".dat");
      if (input_var.display_eigenvectors)
        WriteEigenvectors(var, input_var, AllEigenVec);
    }
    
}

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc == 2)
    {
#ifdef SELDON_WITH_MPI
      int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
      int rank_proc(0);
#endif
      
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

      InputVariables data;
      
      if (rank_proc == 0)
	cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;
      
      if ((type_equation == "TIME_MAXWELL3D_DG")
          || (type_equation == "STATIC_MAXWELL3D_DG"))
        {
          EllipticProblem<StaticMaxwellEquation_3D_DG> Vars;
          RunAll(Vars, data, file_name_data, type_element, type_equation);
        }
      else if (type_equation == "LAPLACE")
        {
	  EllipticProblem<LaplaceEquation<Dimension2> > Vars;
          RunAll(Vars, data, file_name_data, type_element, type_equation);
        }
      else if (type_equation == "HARMONIC_GALBRUN")
        {
	  EllipticProblem<HarmonicGalbrunEquation<Dimension2> > Vars;
          RunAll(Vars, data, file_name_data, type_element, type_equation);
        }
      else
        {
          cout << "Unknown equation" << endl;
          abort();
        }
      
      if (rank_proc == 0)
	cout<<" we destroy the variables "<<endl;
      
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"main.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
}
