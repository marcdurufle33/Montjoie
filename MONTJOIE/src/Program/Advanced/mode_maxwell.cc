// #define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
// #define MONTJOIE_WITH_THREE_DIM

#define SELDON_WITH_PRECONDITIONING

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_HP_HCURL

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/MontjoieMaxwell3D.hxx"
#else
#include "Elliptic/Helmholtz/MontjoieLaplace.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell2D.hxx"
#include "Elliptic/Maxwell/3D/ComputationModeMaxwell.cxx"
#endif



using namespace Montjoie;

// main function
template<class TypeEquation>
int RunAll(EllipticProblem<TypeEquation>& var, string input_file, string name_element, string name_equation)
{
  typedef typename TypeEquation::Complexe Complexe;
  // initialization of physical properties
  var.InitIndices(PhysicalConstant::nb_max_indices);
  var.SetTypeEquation(name_equation);
  
  Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
  ReadLinesFile(input_file, lines_data_file, var.comm_group_mode);
#else
  ReadLinesFile(input_file, lines_data_file);
#endif

  // input file is treated
  ReadInputFile(lines_data_file, var);
    
  All_LinearSolver glob_solver(var);
  ReadInputFile(lines_data_file, glob_solver);
  
  bool split_mesh = true;
  if (var.GetNbProcPerMode() == 1)
    split_mesh = false;

  // mesh and finite element are constructed
  var.ComputeMeshAndFiniteElement(name_element, split_mesh);
  var.PerformOtherInitializations();
    
  // direct solver
  EigenProblemMontjoie<typename TypeEquation::Complexe> eigen_solver(var, glob_solver);
  ReadInputFile(lines_data_file, eigen_solver);
  
  Dimension2 dim; Symmetric property;
  MeshNumbering<Dimension2>& mesh_num = var.GetMeshNumbering(0);
  glob_solver.SelectOptimalLinearSolver(mesh_num.GetOrder(), var.GetNbDof(), dim, property);
  
  // computation of geometric coefficients
  var.ComputeMassMatrix();
    
  // computation of eigenvalues and eigenmodes
  if (var.GetRankProcMode() == 0)
    var.mesh.Write("test.mesh");

  /*DistributedMatrix<Complexe, General, ArrayRowSparse> Dh, Kh;
  GlobalGenericMatrix<Complexe> nat_mat;
  nat_mat.SetCoefMass(1.0);
  nat_mat.SetCoefDamping(0.0);
  nat_mat.SetCoefStiffness(0.0);
  var.AddMatrixWithBC(Dh, nat_mat);
  
  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefDamping(0.0);
  nat_mat.SetCoefStiffness(1.0);
  var.AddMatrixWithBC(Kh, nat_mat);
  
  Dh.WriteText("DhEig.dat");
  Kh.WriteText("KhEig.dat");

  exit(0); */
  
  eigen_solver.ComputeEigenModes();
  
  return 0;
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

      if (type_equation == "MODE_MAXWELL")
        {
          EllipticProblem<ModeMaxwellEquation<Real_wp> > Vars;
          if (Vars.GetRankProcMode() == 0)
            cout<<"Maxwell Solver with triangles of Nedelec's first family "<<endl; 
          
          RunAll(Vars, file_name_data, type_element, type_equation);
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
        }
      else
        {
          EllipticProblem<ModeMaxwellEquation<Complex_wp> > Vars;
          if (Vars.GetRankProcMode() == 0)
            cout<<"Maxwell Solver with triangles of Nedelec's first family "<<endl; 
          
          RunAll(Vars, file_name_data, type_element, type_equation);
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
        }
    }
  
  return FinalizeMontjoie();
}
