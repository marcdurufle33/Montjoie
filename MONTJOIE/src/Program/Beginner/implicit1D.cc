#define MONTJOIE_WITH_ONE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Harmonic/MontjoieHarmonic.hxx"
#include "Instationary/MontjoieTime.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/Helmholtz1D.cxx"
#endif

using namespace Montjoie;


void ConstructAll(EllipticProblem<LaplaceEquation1D>& var_laplace,
		  const string& name_Mh, const string& name_Kh)
{
  var_laplace.ComputeMeshAndFiniteElement(string("EDGE_LOBATTO"));

  var_laplace.ComputePhysicalCoefficients();  
  var_laplace.ComputeMassMatrix();
  
  VectReal_wp mat_mass;
  var_laplace.GetMassMatrix(mat_mass);

  mat_mass.WriteText(name_Mh);
  
  GlobalGenericMatrix<Real_wp> nat_mat;

  Matrix<Real_wp, General, BandedCol> mat_stiff;  
  int N = var_laplace.GetNbDof();
  int r = var_laplace.mesh.GetOrder();
  mat_stiff.Reallocate(N, N, r, r);
  mat_stiff.Zero();
  
  nat_mat.SetCoefMass(0.0);
  var_laplace.AddMatrixFEM(mat_stiff, nat_mat);
  var_laplace.AddBoundaryTerms(mat_stiff, nat_mat);

  mat_stiff.WriteText(name_Kh);
}


// probleme temporel
void RunTimeProblem(const string& input_file)
{
  EllipticProblem<LaplaceEquation1D> var_laplace;
  EllipticProblem<LaplaceEquation1D> var_laplace_2;
  
  ReadInputFile(input_file, var_laplace);
  ReadInputFile(input_file, var_laplace_2);

  Vector<VectString> mesh_data(1);
  mesh_data(0) = var_laplace.GetMeshData(1);
  var_laplace_2.SetMeshData(mesh_data);

  ConstructAll(var_laplace, "Mh1.dat", "Kh1.dat");
  ConstructAll(var_laplace_2, "Mh2.dat", "Kh2.dat");

  var_laplace.GetCoordinateDof().WriteText("Points1.dat");
  var_laplace_2.GetCoordinateDof().WriteText("Points2.dat");
  
}


int main(int argc, char** argv)
{
  if (argc < 2)
    {
      cout << "Donnez un fichier de donnees" << endl;
      abort();
    }
  
  InitMontjoie(argc, argv);
  
  string input_file(argv[1]);
  RunTimeProblem(input_file);
  
  return FinalizeMontjoie();
}
