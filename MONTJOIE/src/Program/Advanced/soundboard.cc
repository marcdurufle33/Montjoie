#define MONTJOIE_WITH_TWO_DIM
// #define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#include "Hyperbolic/Elastic/MontjoieTimeElastic.hxx"

using namespace Montjoie;

class InputEigenProblem : public InputDataProblem_Base
{
public :
  int nb_eigenvalues_per_run;

  InputEigenProblem() { nb_eigenvalues_per_run = 0; }
  
  void SetInputData(const string& keyword, const VectString& parameters)
  {
    if (keyword == "NbEigenvaluesToComputePerRun")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of Soundboard" << endl;
	    cout << "NbEigenvaluesToComputePerRun needs one parameter, for instance :" << endl;
	    cout << "NbEigenvaluesToComputePerRun = N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        nb_eigenvalues_per_run = to_num<int>(parameters(0));
      }    
  }
  
};


template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& input_file,
	    const string& name_element)
{
  typedef Montjoie::Real_wp Real_wp;
  
  // initialization of physical properties
  var.InitIndices(PhysicalConstant::nb_max_indices);
  var.SetTypeEquation("REISSNER_MINDLIN");

  Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(var.comm_group_mode, &rank_proc);
  int nb_proc; MPI_Comm_size(var.comm_group_mode, &nb_proc);
  ReadLinesFile(input_file, lines_data_file, var.comm_group_mode);
#else
  int rank_proc(0), nb_proc(1);
  ReadLinesFile(input_file, lines_data_file);
#endif

  // lecture fichier de donnees
  ReadInputFile(lines_data_file, var);
  
  All_LinearSolver glob_solver(var);
  ReadInputFile(lines_data_file,  glob_solver);

  // construction du numero de dossier
  // si le dossier n'est pas connu
  int num = -1;
  if (var.DOSSIER_output.size() == 0)
    {
      if (rank_proc == 0)
        {
          var.DOSSIER_output = "[STIFFOUT]/"; 
          EcritDossier(var.DOSSIER_output, input_file, num);
        }

#ifdef SELDON_WITH_MPI
      MPI_Bcast_string(var.DOSSIER_output, 0, var.comm_group_mode);
#endif
    }
  
  // calcul maillage et element fini
  var.ComputeMeshAndFiniteElement(name_element);
  
  if (nb_proc == 1)
    var.mesh.Write("test.mesh");
  
  // autres initialisations
  var.PerformOtherInitializations();
  
  // objet declare le pb aux valeurs propres
  EigenProblemMontjoie<Real_wp> eigen_solver(var,  glob_solver);
  
  // lecture du fichier de donnees
  InputEigenProblem input_var;
  ReadInputFile(lines_data_file, eigen_solver);
  ReadInputFile(lines_data_file, input_var);

  int nev = eigen_solver.GetNbAskedEigenvalues();
  int nb_eigenvalues_per_run = input_var.nb_eigenvalues_per_run;
  if (nev > 0)
    {
      var.ComputeMassMatrix();

      if (var.print_level >= 6)
        {
          DistributedMatrix<Real_wp, General, ArrayRowSparse> Mh, Kh, Kh2;
          GlobalGenericMatrix<Real_wp> nat_mat;
          
          nat_mat.SetCoefStiffness(1.0);
          nat_mat.SetCoefDamping(0.0);
          nat_mat.SetCoefMass(0.0);
          var.AddMatrixWithBC(Kh, nat_mat);
          
	  nat_mat.SetCoefStiffness(0.0);
          nat_mat.SetCoefDamping(0.0);
          nat_mat.SetCoefMass(1.0);
          var.AddMatrixWithBC(Mh, nat_mat);
          
          /*FemMatrixFreeClass<TypeElement, TypeEquation> Ah(var);
          nat_mat.SetCoefStiffness(1.0);
          nat_mat.SetCoefDamping(0.0);
          nat_mat.SetCoefMass(0.0);
          var.AddMatrixWithBC(Ah, nat_mat);
          
          int N = var.GetNbDof();
          Kh2.Reallocate(N, N);
          VectReal_wp Ones(N), KhOnes(N);
          Ones.Fill(0); KhOnes.Fill(0);
          for (int j = 0; j < Ah.GetM(); j++)
            {
              Ones.Fill(0); Ones(j) = 1.0;
              MltAdd(1.0, Ah, Ones, 0.0, KhOnes);
              
              for (int i = 0; i < N; i++)
                if (abs(KhOnes(i)) > 1e-12)
                  Kh2.AddInteraction(i, j, KhOnes(i));
            }
          */
          //Mh.WriteText("Mh.dat");
          //Kh.WriteText("Kh.dat");
          //Kh2.WriteText("Kh2.dat");
        }
      
      if (nb_eigenvalues_per_run > 0)
        eigen_solver.SetNbAskedEigenvalues(nb_eigenvalues_per_run);
      else
        eigen_solver.SetNbAskedEigenvalues(nev+4);

      FemMatrixFreeClass_Base<Real_wp>* Kh;
      FemMatrixFreeClass_Base<Real_wp>* Mh;
      
      Mh = var.GetNewIterativeMatrix(Real_wp(0));
      Kh = var.GetNewIterativeMatrix(Real_wp(0));
      
      eigen_solver.InitMatrix(*Kh, *Mh);
      
      // on alloue les tableaux resultats
      VectReal_wp eigen_values;
      Matrix<Real_wp> eigen_vectors;
      eigen_values.Reallocate(nev);
      eigen_vectors.Reallocate(eigen_solver.GetM(), nev);

      VectReal_wp lambda_r, lambda_i;
      Matrix<Real_wp, General, ColMajor> eigen_vec;
      Real_wp omega_cur = 0;
      Real_wp shift = eigen_solver.GetShiftValue(), shift_imag = 0;
      int n = 0; Real_wp coef = 1.0+1e-8;
      while (n < nev)
        {
          DISP(n); DISP(nev);
          // on calcule les valeurs propres pour le shift choisi
          var.UpdateShiftAdimensionalization(shift, shift_imag);
          eigen_solver.SetShiftValue(shift); eigen_solver.SetImagShiftValue(shift_imag);
          
          GetEigenvaluesEigenvectors(eigen_solver, lambda_r, lambda_i, eigen_vec);
          
          var.UpdateEigenvaluesAdimensionalization(lambda_r,
                                                   lambda_i, eigen_vec);
          
          if (rank_proc == 0)
            cout << "Temporary eigenvalues :" << endl << lambda_r << endl << endl;
                    
          // on stocke toutes les nouvelles valeurs propres
          for (int i = 0; i < lambda_r.GetM()-4; i++)
            {
              if ((lambda_r(i) > coef*omega_cur) && (n < nev))
                {
                  eigen_values(n) = lambda_r(i);
                  for (int j = 0; j < eigen_vec.GetM(); j++)
                    eigen_vectors(j, n) = eigen_vec(j, i);
                  
                  omega_cur = lambda_r(i);
                  n++;
                }
            }
          
          // on passe aux valeurs propres suivantes
          shift = omega_cur + 0.35*abs(omega_cur - lambda_r(0));
        }

      delete Kh; delete Mh;
 
      VectReal_wp eigen_mode(var.GetNbDof());
      Vector<Vector<Real_wp> > eigen_mode_v;
      eigen_mode_v.SetData(1, &eigen_mode);
      
      const IVect& local_col_numbers = eigen_solver.GetLocalColumnNumbers();
      for (int j = 0; j < nev; j++)
        {
          eigen_mode.Fill(0);
          for (int i = 0; i < eigen_vectors.GetM(); i++)
            eigen_mode(local_col_numbers(i)) = eigen_vectors(i, j);
          
          if (nb_proc > 1)
            var.AddDomains(eigen_mode);
          
          string numero = NumberToString(j);
          for (int k = 0; k < var.output_mesh_param.GetM(); k++)
            var.output_mesh_param(k).
              SetFileName(2, GetBaseString(var.output_mesh_param(k)
                                           .GetTotalFieldFile())
                          + "." + numero.substr(numero.size()-3,3) + string(".vtk"));
          
          for (int k = 0; k < var.output_grid_param.GetM(); k++)
            var.output_grid_param(k).
              SetFileName(2, GetBaseString(var.output_grid_param(k)
                                           .GetTotalFieldFile())
                          + numero + string(".dat"));
          
          var.WriteOutputFile(eigen_mode_v, 2);
        }
      
      eigen_values.WriteText(var.DOSSIER_output+"eigen_values.dat");
      for (int i = 0; i < eigen_values.GetM(); i++)
        eigen_values(i) = sqrt(eigen_values(i))/(2.0*pi_wp);
      
      eigen_values.WriteText(var.DOSSIER_output+"frequency.dat");
      
      eigen_mode_v.Nullify();
    }
}



int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      // on recupere le nom du fichier de donnees
      string file_name_data;
      file_name_data = string(argv[1]);
      DISP(file_name_data);

      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      
      EllipticProblem<ReissnerMindlinEquation> Vars;

      RunAll(Vars, file_name_data, type_element);
      cout<<" we destroy the variables "<<endl;
      
      
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"main.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  cout<<"End of the program"<<endl; 

  return FinalizeMontjoie();
}
