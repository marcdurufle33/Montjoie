#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Hyperbolic/Elastic/MontjoieTimeElastic.hxx"
#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"

#include "Hyperbolic/Elastic/FluidStructureInteraction.hxx"
#include "Hyperbolic/Elastic/FluidStructureInteractionInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Elastic/FluidStructureInteraction.cxx"
#endif

using namespace Montjoie;

template<class Dimension>
void CheckHelmholtzSolution(AcousticAcousticInteraction<Dimension>& var,
			    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& M,
			    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& K)
{
  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> A(K);
  Add(-var.var_fluid.var_harmonic.GetSquareOmega(), M, A);
  
  All_MatrixLU<Real_wp> mat_lu;
  
  //A.WriteText("Mat.dat");
  mat_lu.ShowMessages();
  mat_lu.Factorize(A);
  
  int N = var.GetNbDof();
  VectReal_wp rhs(N);
  rhs.Fill(0);
  
  int Nu = var.var_solid.var_harmonic.GetNbDof();
  int Np = var.var_fluid.var_harmonic.GetNbDof();
  VectReal_wp rhs_elas(Nu);
  VectReal_wp rhs_acous(Np);  
        
  var.var_fluid.var_harmonic.ComputeRightHandSide(rhs_acous);
  var.var_solid.var_harmonic.ComputeRightHandSide(rhs_elas);
        
  for (int i = 0; i < Nu; i++)
    rhs(i) = rhs_elas(i);

  for (int i = 0; i < Np; i++)
    rhs(Nu+i) = rhs_acous(i);

  mat_lu.Solve(rhs);

  for (int i = 0; i < Nu; i++)
    rhs_elas(i) = rhs(i);

  for (int i = 0; i < Np; i++)
    rhs_acous(i) = rhs(Nu+i);
  
  var.var_fluid.var_harmonic.WriteDatas(rhs_acous);
  var.var_solid.var_harmonic.WriteDatas(rhs_elas);
}

template<class Dimension>
void RunAll(AcousticAcousticInteraction<Dimension>& var, const string& data_file,
	    const string& name_element)
{
  typedef Montjoie::Real_wp Real_wp;
  
  var.ConstructAll(data_file, name_element);

  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> M, K;
  var.ComputeMatrices(M, K);

  // pour verifier que Helmholtz est bien resolu
  // CheckHelmholtzSolution(var, M, K);
  
  M.WriteText("mat_mass.dat");
  K.WriteText("mat_stiff.dat");

  // on peut extraire Mh1, Mh2, Kh1 et Kh2
  int N1 = var.var_solid.var_harmonic.GetNbDof();
  int N2 = var.var_fluid.var_harmonic.GetNbDof();
  int Nl = var.DdlLambda.GetM();
  VectReal_wp Mh1(N1), Mh2(N2);
  Matrix<Real_wp, Symmetric, ArrayRowSymSparse> Kh1(N1, N1), Kh2(N2, N2);
  Matrix<Real_wp, General, ArrayRowSparse> Ch1(N1, Nl), Ch2(N2, Nl);
  
  int offset = N1+N2;
  for (int i = 0; i < N1; i++)
    {
      Mh1(i) = M(i, i);
      for (int j = 0; j < K.GetRowSize(i); j++)
	{
	  int jcol = K.Index(i, j);
	  if (jcol < N1)
	    Kh1.AddInteraction(i, jcol, K.Value(i, j));
	  else
	    Ch1.AddInteraction(i, jcol-offset, K.Value(i, j));
	}
    }

  for (int i = 0; i < N2; i++)
    {
      Mh2(i) = M(N1+i, N1+i);
      for (int j = 0; j < K.GetRowSize(N1+i); j++)
	{
	  int jcol = K.Index(N1+i, j);
	  if (jcol < offset)
	    Kh2.AddInteraction(i, jcol-N1, K.Value(N1+i, j));
	  else
	    Ch2.AddInteraction(i, jcol-offset, K.Value(N1+i, j));
	}
    }
  
  Mh1.WriteText("Mh1.dat");   Mh2.WriteText("Mh2.dat");
  Kh1.WriteText("Kh1.dat");   Kh2.WriteText("Kh2.dat");
  Ch1.WriteText("Ch1.dat");   Ch2.WriteText("Ch2.dat");
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  // premier argument -> fichier de donnees
  if (argc != 2)
    {
      cout << "Cette commande demande un argument " << endl;
      cout << "Par exemple, ./implicit_ND.x EXECUTION/data_file.ini " << endl;
      abort();
    }
  
  string name_data_file(argv[1]);
  
  // we get the type of element selected by the user, and type of equation
  string type_element, type_equation;
  getElement_Equation(name_data_file, type_element, type_equation);
  
  int dim_N = 3;
  if ((type_element.find("TRIANGLE") == 0) || (type_element.find("QUADRANGLE") == 0))
    dim_N = 2;
  
  if (dim_N == 3)
    {
      AcousticAcousticInteraction<Dimension3> var;
      RunAll(var, name_data_file, type_element);
    }
  else
    {
      AcousticAcousticInteraction<Dimension2> var;
      RunAll(var, name_data_file, type_element);
    }
  
  return FinalizeMontjoie();
}

