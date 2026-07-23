#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Elliptic/Helmholtz/MontjoieLaplace.hxx"

using namespace Montjoie;

template<class GenericPb>
bool PointClosePML(const R2& ptA, GenericPb& var)
{
  if ( (ptA(0) < var.GetXmin()+R2::threshold) || 
       (ptA(0) > var.GetXmax()-R2::threshold) ||
       (ptA(1) < var.GetYmin()+R2::threshold) || 
       (ptA(1) > var.GetYmax()-R2::threshold) )
    return true;
  
  return false;
}

template<class GenericPb>
bool PointClosePML(const R3& ptA, GenericPb& var)
{
  if ( (ptA(0) < var.GetXmin()+R2::threshold) || 
       (ptA(0) > var.GetXmax()-R2::threshold) ||
       (ptA(1) < var.GetYmin()+R2::threshold) || 
       (ptA(1) > var.GetYmax()-R2::threshold) ||
       (ptA(2) < var.GetZmin()+R2::threshold) || 
       (ptA(2) > var.GetZmax()-R2::threshold) )
    return true;
  
  return false;
}

// calcule ||u_{pml} || / || u_{all} ||
template<class TypeEquation, class Vector1>
double GetPmlNorm(EllipticProblem<TypeEquation>& var,
                  Vector1& eigen_mode)
{
  Vector<bool> BoundaryVertex(var.mesh.GetNbVertices());
  BoundaryVertex.Fill(false);
  R2 ptA;
  for (int i = 0; i < var.mesh.GetNbVertices(); i++)
    {
      ptA = var.mesh.Vertex(i);
      if (PointClosePML(ptA, var))
        BoundaryVertex(i) = true;
    }
  
  double sum_pml = 0, sum_all = 0;
  Vector<bool> dof_domain(var.mesh_num.GetNbDof());
  dof_domain.Fill(false);
  for (int i = 0; i < var.mesh.GetNbElt(); i++)
    if (!var.InsidePML(i))
      for (int j = 0; j < var.mesh_num.Element(i).GetNbDof(); j++)
	dof_domain(var.mesh_num.Element(i).GetNumberDof(j)) = true;

  for (int i = 0; i < var.mesh_num.GetNbDof(); i++)
    {
      sum_all += abs(eigen_mode(i));
      if (!dof_domain(i))
	sum_pml += abs(eigen_mode(i));
    }
  
#ifdef SELDON_WITH_MPI
  int nb_proc; MPI_Comm_size(var.comm_group_mode, &nb_proc);
  if (nb_proc > 1)
    {
      R2 sum_tmp(sum_pml, sum_all), sum_res;
      MPI_Allreduce(&sum_tmp(0), &sum_res, 2, GetMpiDataType(sum_tmp(0)), MPI_MAX, var.comm_group_mode);
      sum_pml = sum_res(0);
      sum_all = sum_res(1);
    }
#endif
  
  return sum_pml/sum_all;
}

// regeneration du maillage pour le pas de maillage lc donne
// dans le maillage, on doit avoir deux lignes du type :
// lc = 0.5;
// lambda = 0.3;
// lambda est la longueur d'onde, on essaiera de construire le maillage
// en mettant un ecart entre l'objet et le bord du domaine egal a un multiple
// de lambda (par exemple 2*lambda ou 3*lambda) 
void ModifyMeshFile(const Real_wp& wavelength, const Real_wp& lc, const string& file_name,
                    const string& file_base_geo)
{
  string file_geo = GetBaseString(file_name) + ".geo";
  ifstream file_in(file_base_geo.data());
  ofstream file_out(file_geo.data());
  
  if (!file_in.is_open())
    {
      cout << "Mesh file not found" << endl;
      abort();
    }
  
  string ligne;
  while (file_in.good())
    {
      getline(file_in, ligne);
      if (ligne.substr(0, 4) == "lc =")
        file_out << "lc = " << lc << ";" << endl;          
      else if (ligne.substr(0, 8) == "lambda =")
        file_out << "lambda = " << wavelength << ";" << endl;
      else
        file_out << ligne << endl;
    }
  
  file_out.close();
  
  string command_line = "geo2mesh " + file_geo + " >& sort" + file_geo;
  system(command_line.data());  
}

// classe contenant les donnees pour trouver les valeurs propre sur un intervalle
class InputVariablesEigenvaluePML : public InputDataProblem_Base
{
public :
  // intervalle ou on va chercher kr
  Real_wp omega_min_spectrum, omega_max_spectrum;
  // taux utilise pour raffiner le maillage au fur et a mesure que lambda diminue
  // on prend la regle h = coef * lambda
  Real_wp coef_refinement;
  // nom du fichier contenant les valeurs propres
  string file_name_eigenvalue;
  // si true, toutes les valeurs propres sont stockees meme celles en dehors de l'intervalle
  bool storage_all_eigenvalues;
  // nombre de valeurs propres calculees lors de la premiere iteration, et lors des iterations suivantes
  int first_number_eigenval, second_number_eigenval, add_number_eigenval;
  // fichier .geo du maillage
  string file_base_geo;
  // si true, les valeurs propres sont rajoutees a la fin du fichier
  bool append_eigenvalue;
  // si true, on ne modifie pas le maillage
  bool mesh_not_modified;
  // tolerance utilisee pour le calcul des valeurs propres
  double tolerance_eigenvalue;
  // non-vide si on veut utiliser differents maillages lors du calcul
  // si omega < omega_subdiv(0), on utilise le maillage contenu dans file_mesh(0)
  // si omega < omega_subdiv(1), on utilise le maillage contenu dans file_mesh(1)
  // ...
  VectReal_wp omega_subdiv;
  Vector<string> file_mesh;
  
  InputVariablesEigenvaluePML()
  {
    omega_min_spectrum = 1.0; omega_max_spectrum = 10.0;
    coef_refinement = 1.0;
    storage_all_eigenvalues = true;
    first_number_eigenval = 50;
    second_number_eigenval = 100;
    add_number_eigenval = 0;
    append_eigenvalue = false;
    mesh_not_modified = true;
    tolerance_eigenvalue = 1e-12;
  }
  
  void SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (keyword == "MinimalFrequencySpectrum")
      {
        omega_min_spectrum = 2.0*pi_wp*to_num<Real_wp>(param(0)) + to_num<Real_wp>(param(1));
      }
    else if (keyword == "MaximalFrequencySpectrum")
      {
        omega_max_spectrum = 2.0*pi_wp*to_num<Real_wp>(param(0)) + to_num<Real_wp>(param(1));
      }
    else if (keyword == "StoreAllEigenvalues")
      {
        if (param(0) == "YES")
          storage_all_eigenvalues = true;
        else
          storage_all_eigenvalues = false;
      }
    else if (keyword == "FileEigenvalue")
      {
        file_name_eigenvalue = param(0);
      }
    else if (keyword == "CoefRefinement")
      {
        coef_refinement = to_num<Real_wp>(param(0));
      }
    else if (keyword == "NbEigenvaluesToCompute")
      {
        first_number_eigenval = to_num<int>(param(0));
        second_number_eigenval = to_num<int>(param(1));
        if (param.GetM() >= 3)
          add_number_eigenval = to_num<int>(param(2));
      }
    else if (keyword == "FileBaseGeo")
      file_base_geo = param(0);
    else if (keyword == "AppendEigenvalue")
      {
        if (param(0) == "YES")
          append_eigenvalue = true;
        else
          append_eigenvalue = false;
      }
    else if (keyword == "KeepMesh")
      {
        if (param(0) == "YES")
          mesh_not_modified = true;
        else
          mesh_not_modified = false;
      }
    else if (keyword == "ListMeshesToUse")
      {
        int n = param.GetM()/2;
        omega_subdiv.Reallocate(n);
        file_mesh.Reallocate(n+1);
        
        if (param.GetM() >= 1)
          file_mesh(0) = param(0);
        
        for (int i = 0; i < n; i++)
          {
            omega_subdiv(i) = to_num<Real_wp>(param(2*i+1));
            file_mesh(i+1) = param(2*i+2);
          }
      }
    else if (keyword == "EigenvalueTolerance")
      {
        tolerance_eigenvalue = to_num<Real_wp>(param(0));
      }
  }
  
};

// fonction generale pour calculer tous les modes
template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& input_file,
	    const string& name_element)
{
  // on lit le fichier de donnees
  Vector<string> lines_data_file;
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(var.comm_group_mode, &rank_proc);
  ReadLinesFile(input_file, lines_data_file, var.comm_group_mode);
#else
  int rank_proc(0);
  ReadLinesFile(input_file, lines_data_file);
#endif
  
  // initialisation des parametres du probleme
  var.InitIndices(PhysicalConstant::nb_max_indices);
  
  ReadInputFile(lines_data_file, var);

  InputVariablesEigenvaluePML data;
  ReadInputFile(lines_data_file, data);
  
  All_LinearSolver glob_solver(var);  
  ReadInputFile(lines_data_file, glob_solver);

  //VarComputationRCS<TypeElement, TypeEquation> rcs_param(var);
  //ReadInputFile(lines_data_file, rcs_param);
  
  // test de sortie de la boucle
  // tant que all_eigenvalues_not_computed, on continue a les chercher
  bool all_eigenvalues_not_computed = true;
  // numero de la valeur propre en cours
  int num_eigen = 0;
  // shift initial utilise
  Real_wp omega = var.GetOmega(); 
    
  // fichier qui va contenir toutes les valeurs propres
  ofstream file_eigenvalue;
  file_eigenvalue.precision(15);
  
  // omega_cur : pulsation de la valeur propre courante
  Real_wp omega_cur = 0;
  if (!data.append_eigenvalue)
    {
      // on cree un nouveau fichier
      if (rank_proc == 0)
        {
          file_eigenvalue.open(data.file_name_eigenvalue.data());
          file_eigenvalue.close();
        }
    }
  else
    {
      // les valeurs propres sont rajoutees a la suite du fichier
      // on recupere les anciennes valeurs propres pour initialiser num_eigen et omega_cur
      Matrix<Real_wp> old_eigen;
      old_eigen.ReadText(data.file_name_eigenvalue);
      num_eigen = 2*old_eigen.GetM();
      if (rank_proc == 0)
        cout << "Restarting computation with eigenvalue " << num_eigen << endl;
      
      omega_cur = old_eigen(old_eigen.GetM()-1, 1);
    }
  
  Real_wp omega0 = omega_cur;  

  FemMatrixFreeClass_Base<Real_wp>* Kh;
  FemMatrixFreeClass_Base<Real_wp>* Mh;

  Kh = var.GetNewIterativeMatrix(Real_wp(0));
  Mh = var.GetNewIterativeMatrix(Real_wp(0));

  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> M;
  DistributedMatrix<Real_wp, General, ArrayRowSparse> K;

  VectComplex_wp AllLambda;
  Vector<VectComplex_wp> AllEigenVec;
  
  // boucle principale
  while (all_eigenvalues_not_computed)
    {
      var.Restart();
      
      // longueur d'onde lambda
      Real_wp wavelength = 2.0*pi_wp/omega;
      
      // pas de maillage
      Real_wp lc = data.coef_refinement*wavelength;
      
      if ((data.omega_min_spectrum < data.omega_max_spectrum) && (!data.mesh_not_modified))
        {
          if (data.file_mesh.GetM() > 0)
            {
              int n = data.file_mesh.GetM()-1;
              for (int k = n-1; k >= 0; k--)
                {
                  if (omega < data.omega_subdiv(k))
                    {
                      n = k;
                      break;
                    }
                }
              
              // maillage a utiliser
              var.mesh_data.Reallocate(1);
              var.mesh_data(0) = data.file_mesh(n);              
            }
          else
            {
              if (rank_proc == 0)
                ModifyMeshFile(wavelength, lc, var.mesh_data(0)(0), data.file_base_geo);
            }
        }
      
      if (rank_proc == 0)
        cout << "Lambda = " << wavelength << endl;
      
      // on regle l'epaisseur de la PML
      var.mesh.GetPmlArea(0).SetThicknessPML(wavelength);
      
      // on construit le maillage et les elements finis
      var.ComputeMeshAndFiniteElement(name_element);
      
      // on ecrit le maillage
      if (rank_proc == 0)
	var.mesh.Write("test"+to_str(num_eigen)+".mesh");
      
      // on rajoute les ddls en plus
      var.PerformOtherInitializations();
      
      // calcul des coefficients geometriques : DF_i, J_i, etc
      var.ComputeMassMatrix();
      
      GlobalGenericMatrix<Real_wp> nat_mat;
      nat_mat.SetCoefMass(1.0);
      nat_mat.SetCoefStiffness(0.0);
      nat_mat.SetCoefDamping(0.0);
      
      cout << "calcul de la matrice de masse " << endl;
      FemMatrixFreeClass<Real_wp, TypeEquation> free_mat(var);

      M.Clear();
      var.AddMatrixWithBC(free_mat, M, nat_mat);  
      M.WriteText("Mh.dat");
      
      cout << endl << endl << endl << endl;
      cout << "calcul de la matrice de rigidite " << endl;
      
      nat_mat.SetCoefMass(0.0);
      nat_mat.SetCoefStiffness(1.0);
      nat_mat.SetCoefDamping(1.0);
      K.Clear();
      var.AddMatrixWithBC(free_mat, K, nat_mat);  
      K.WriteText("Kh.dat");

      SparseEigenProblem<Real_wp, Matrix<Real_wp, General, ArrayRowSparse> > eigen_solver;
      //cout << endl << endl << endl << endl;
      //cout << "calcul de la matrice ponderee " << endl; */
      
      // initialisation du probleme aux valeurs propres
      //EigenProblemMontjoie<Real_wp> eigen_solver(var, glob_solver);
      //eigen_solver.InitMatrix(*Kh, *Mh);
      
      eigen_solver.SetStoppingCriterion(data.tolerance_eigenvalue);
      eigen_solver.SetNbMaximumIterations(1000);
      if (num_eigen == 0)
        eigen_solver.SetNbAskedEigenvalues(data.first_number_eigenval);
      else
        eigen_solver.SetNbAskedEigenvalues(data.second_number_eigenval);
      
      eigen_solver.SetNbAdditionalEigenvalues(data.add_number_eigenval);
      
      eigen_solver.InitMatrix(K, M);
      eigen_solver.SetDiagonalMass(true);
      DISP(eigen_solver.DiagonalMass());
      DISP(omega); DISP(eigen_solver.GetNbAskedEigenvalues());
      DISP(eigen_solver.GetNbArnoldiVectors());
      // on cherche les valeurs propres proches de sigma = (-0.12, omega+0.1)
      // les valeurs propres sont triees par module croissant
      eigen_solver.SetComputationalMode(eigen_solver.SHIFTED_MODE);
      eigen_solver.SetTypeSpectrum(eigen_solver.CENTERED_EIGENVALUES,
				   complex<double>(-0.12, omega+0.1), eigen_solver.SORTED_MODULUS);
      
      eigen_solver.SetPrintLevel(1);
      DISP(eigen_solver.DiagonalMass());
      // calcul des valeurs propres
      Vector<double> lambda, lambda_imag;
      Matrix<double> eigen_vec;
      GetEigenvaluesEigenvectors(eigen_solver, lambda, lambda_imag, eigen_vec);
      
      lambda.WriteText("eigenvalue_real.dat");
      lambda_imag.WriteText("eigenvalue_imag.dat");
      eigen_vec.WriteText("eigenvec.dat");
      if (rank_proc == 0)
        {
          DISP(lambda); DISP(lambda_imag);
        }
      
      // on ecrit les valeurs et vecteurs propres
      int n = eigen_vec.GetM();
      Vector<double> eigen_mode(n), eigen_mode_imag(n);
      VectComplex_wp eigen_mode_cplx(n);
      Vector<Vector<double> > eigen_mode_v, eigen_mode_imag_v;
      eigen_mode_v.SetData(1, &eigen_mode);
      eigen_mode_imag_v.SetData(1, &eigen_mode_imag);
      
      // seuil pour determiner si un mode est PML ou pas
      double seuil = 0.99;
      
      if (rank_proc == 0)
        file_eigenvalue.open(data.file_name_eigenvalue.data(), ios::app);
      
      int neig = 0;
      if (rank_proc == 0)
        cout << "On prend les valeurs propres plus grandes que " << omega_cur << endl;

      // boucle sur les valeurs propres (considerees complexes conjuguees
      for (int i = 0; i < eigen_solver.GetNbAskedEigenvalues()-4; i += 2)
        if (abs(lambda_imag(i)) >= omega_cur)
          if ( data.storage_all_eigenvalues  || (abs(lambda_imag(i)) < data.omega_max_spectrum) )
            {
	      if (lambda_imag(i) != -lambda_imag(i+1))
		{
		  cout << "Eigenvalues not complex conjugate ?" << endl;
		  abort();
		}
	      
              if (neig == 0)
                omega0 = abs(lambda_imag(i));
              
              for (int j = 0; j < n; j++)
                eigen_mode(j) = eigen_vec(j, i);

	      for (int j = 0; j < n; j++)
                eigen_mode_imag(j) = eigen_vec(j, i+1);

	      for (int j = 0; j < n; j++)
		eigen_mode_cplx(j) = complex<double>(eigen_vec(j, i), eigen_vec(j, i+1));
	      
              // si la norme de la solution est trop grande pres de la PML
              // on considere que c'est un mode PML
              double NormeExt = GetPmlNorm(var, eigen_mode_cplx);
              if ((NormeExt < seuil) || (data.storage_all_eigenvalues))
                {
                  // on ecrit la valeur propre
                  if (rank_proc == 0)
		    file_eigenvalue << num_eigen/2 + 1 << " " << abs(lambda_imag(i)) << " " << lambda(i) << " " << endl;
		  
                  // on ecrit le vecteur propre
                  string numero = NumberToString(num_eigen);
                  for (int j = 0; j < var.output_mesh_param.GetM(); j++)
                    var.output_mesh_param(j).
                      SetFileName(2, GetBaseString(var.output_mesh_param(j).GetTotalFieldFile())
                                  + "." + numero.substr(numero.size()-3,3) + string(".bb"));
                  
                  for (int j = 0; j < var.output_grid_param.GetM(); j++)
                    var.output_grid_param(j).
                      SetFileName(2, GetBaseString(var.output_grid_param(j).GetTotalFieldFile())
                                  + numero + string(".dat"));
                  
                  var.WriteOutputFile(eigen_mode_v, 2);

		  numero = NumberToString(num_eigen+1);
                  for (int j = 0; j < var.output_mesh_param.GetM(); j++)
                    var.output_mesh_param(j).
                      SetFileName(2, GetBaseString(var.output_mesh_param(j).GetTotalFieldFile())
                                  + "." + numero.substr(numero.size()-3,3) + string(".bb"));
                  
                  for (int j = 0; j < var.output_grid_param.GetM(); j++)
                    var.output_grid_param(j).
                      SetFileName(2, GetBaseString(var.output_grid_param(j).GetTotalFieldFile())
                                  + numero + string(".dat"));
                  
                  var.WriteOutputFile(eigen_mode_imag_v, 2);
		  
                  num_eigen += 2;
                  //cout << "mode " << i << " interesting " << endl;
                  //DISP(lambda(i));
                  //DISP(lambda_imag(i));
                  //DISP(NormeExt);
                }
              else
                {
                  cout << "Mode " << i << " ignored " << endl;
                }
              
              omega_cur = abs(lambda_imag(i));
              neig++;
            }
      
      if (rank_proc == 0)
        file_eigenvalue.close();
      
      // on passe aux valeurs propres suivantes
      omega = omega_cur + 0.35*abs((omega_cur - omega0));
      //if (omega_cur > data.omega_max_spectrum)
      all_eigenvalues_not_computed = false;
      
      eigen_mode_v.Nullify();
      eigen_mode_imag_v.Nullify();
    }
  
  delete Kh; delete Mh;
  
  // computation successful, we exit the function
  return;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc != 2)
    {
      cout << "Donnez un fichier de donnees " << endl;
      abort();
    }
  
  EllipticProblem<LaplaceEquation<Dimension2> > vars;
  string input_file(argv[1]);
  RunAll(vars, input_file, string("TRIANGLE_CLASSICAL"));
  
  return FinalizeMontjoie();
}

