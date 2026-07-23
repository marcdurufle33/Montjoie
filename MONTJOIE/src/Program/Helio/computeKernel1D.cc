#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

//#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

//#include "Elliptic/Helmholtz/AxiSymHelmholtz.hxx"
//#include "Elliptic/Helmholtz/AxiSymHelmholtzInline.cxx"

//#ifndef SELDON_WITH_COMPILED_LIBRARY
//#include "Elliptic/Helmholtz/AxiSymHelmholtz.cxx"
//#endif

#include "Harmonic/MontjoieHarmonic.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
#include "Elliptic/Helmholtz/HelmholtzRadial.hxx"
//#include "Helio/Solvers/HelmholtzRadial.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D_Inline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/Helmholtz1D.cxx"
#include "Elliptic/Helmholtz/HelmholtzRadial.cxx"
//#include "Helio/Solvers/HelmholtzRadial.cxx"
#endif

#include "Helio/Kernel1D.cxx"
#include <time.h>
using namespace Montjoie;


// ======================================================================================================

void getGreen2(string input_file)
{
  typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;

  MontjoieTimer var_chrono;

  var_chrono.SetMessage("GlobalComputation", "perform the overall simulation");
  var_chrono.Start("GlobalComputation");
  
  var_chrono.SetMessage("ConstructAll", "initialize the computation");
  var_chrono.Start("ConstructAll");

  string name_elt;// = "EDGE_LOBATTO";
  string type_equation;// = "HELMHOLTZ";
  getElement_Equation(input_file, name_elt, type_equation);
  EllipticProblem<HelmholtzEquationRadial> var;
  
  // initialisation of computations
  var.ConstructAll(input_file, name_elt, var_chrono);
  
  var_chrono.Stop("ConstructAll");

  // computation of the right hand sides
  var_chrono.SetMessage("ComputationRhs", "compute the right hand side");
  var_chrono.Start("ComputationRhs");

  Vector<Vector<Vector<Complex_wp> > > rhs;
  var.ComputeRightHandSide(rhs, var_chrono);
  
  var_chrono.Stop("ComputationRhs");

  var_chrono.SetMessage("FactorisationMatrix", "compute and factorize finite element matrix");
  var_chrono.SetMessage("ComputationSolution", "solves linear systems");

  Matrix<Complex_wp, General, BandedCol> mat_lu;
  Vector<int> pivot;
  int lMax = 31;
  int nR = 100;

  Vector<Real_wp> x(nR);
  Vector<int>     eltS(nR); eltS.Zero();
  Vector<Real_wp> locPos(nR); locPos.Zero();

  for (int i =0; i < nR; i++)
    {
      x(i) = 1.*i/nR;

      while ((eltS(i)<var.mesh.GetNbVertices()) && (var.mesh.Vertex(eltS(i))<x(i)))
	eltS(i)++;

      // Point at a mesh boundary or very close
      if (eltS(i)==0) {
	if (abs(var.mesh.Vertex(0)-x(i))<=10*epsilon_machine)
	  locPos(i) = 0.e0;
	else {
	  cout << "point " << x(i) << " not in the mesh " << endl;
	  abort();
	}}
      else if (eltS(i)==var.mesh.GetNbVertices()) {
	if (abs(var.mesh.Vertex(eltS(i)-1)-x(i)) <= 10*epsilon_machine) {
	  eltS(i)    = var.mesh.GetNbElt()-1;
	  locPos(i) = 1.e0;
	}
	else {
	  cout << "point " << x(i) << " not in the mesh " << endl;
	  abort();
	}}
      // point inside the mesh
      else {
	eltS(i)--;
	locPos(i) = (x(i)-var.mesh.Vertex(eltS(i)))/(var.mesh.Vertex(eltS(i)+1)-var.mesh.Vertex(eltS(i)));
      }
    }


  //Vector<Complex_wp> rhs;
  Matrix<Complex_wp> U(lMax, nR);
  Matrix<Complex_wp> drU(lMax, nR);
  U.Zero(); drU.Zero();

  for (int l=0; l < lMax; l++)
    {
      var_chrono.Start("FactorisationMatrix");
      
      // computation and factorisation of finite element matrix
      var.PerformFactorizationStep(mat_lu, pivot, l);
      
      var_chrono.Stop("FactorisationMatrix");
      
      var_chrono.Start("ComputationSolution");

      if (rhs(l)(0).GetM() > 0)
          {
            var.ComputeSolution(mat_lu, pivot, rhs(l)(0));
          }
	var_chrono.Stop("ComputationSolution");

	// 2. Reconstruct the solutions
	VectReal_wp phi;

	for (int i =0; i < nR; i++)
	  {
	    const ElementReference<Dimension1,1>& Fb = var.GetReferenceElement(eltS(i));
	    Fb.ComputeValuesPhiRef(locPos(i),phi);
	    for (int IDloc=0; IDloc<phi.GetM(); IDloc++) 
	      {
		U(l,i) += rhs(l)(0)(var.mesh.GetNumberDof(eltS(i),IDloc)) * phi(IDloc);
	      }
	  }
	/*
	  int nb_elt = var.GetNbElt();
	DISP(rhs(l)(0).GetM());
	const ElementReference<Dimension1, 1>& Fb0 = var.GetReferenceElement(0);
	DISP((nb_elt-1) * Fb0.GetNbPointsQuadratureInside() + Fb0.GetNbPointsQuadratureInside()-1);
	for (int ID_elt = 0; ID_elt < nb_elt; ID_elt++)
	  {	
	    const ElementReference<Dimension1, 1>& Fb = var.GetReferenceElement(ID_elt);
	    for (int i =0; i < nR; i++)
	      {
		Fb.ComputeValuesPhiRef(x(i),phi);
		if (Norm2(rhs(l)(0)) > 1.e6)
		  cout << l << ", " << i << endl;
		int N   = Fb.GetNbPointsQuadratureInside();
		for (int ID_pt = 0; ID_pt < N; ID_pt++)
		  {
		    int cpt = ID_elt*N +ID_pt;
		    //if (i==0)
		    //  cout << cpt << ", " << rhs(l)(0)(cpt) << endl;
		    U(l,i) += rhs(l)(0)(cpt) * phi(ID_pt);
		  }
	      }
	  }
	*/
    }
  U.Write("Ul.dat");


}


int main(int argc, char **argv)
{  
  InitMontjoie(argc, argv);


  //#ifdef SELDON_WITH_MPI
  //MPI_Init(&argc, &argv);
  //#endif

  if (argc > 1)
    {
      string input_file = string(argv[1]); 
      //getGreen(input_file);
      // we get the type of element selected by the user, and type of equation
      //string name_elt;// = "EDGE_LOBATTO";
      //string type_equation;// = "HELMHOLTZ";
      //getElement_Equation(input_file, name_elt, type_equation);
      //EllipticProblem<HelmholtzEquationRadial> var;
      //RunAll(var, input_file, name_elt);

      string name_elt = "EDGE_LOBATTO";
      string type_equation = "HELMHOLTZ";

      Kernel1D kernel = Kernel1D();
      int indFreqMin = 0;
      if (argc > 2)
	indFreqMin = to_num<int>(argv[2]);
      kernel.ReadInputFile(input_file, indFreqMin);
      //Matrix<Complex_wp> G(kernel.lMaxG, kernel.nr);
      //Matrix<Complex_wp> drG(kernel.lMaxG, kernel.nr);
      //G.Zero(); drG.Zero();
      //cout << "Computation of the Green's function ... ";
      //kernel.getGreen(input_file);
      //cout << " done" << endl;
      if (kernel.kernelDimension == "2D")
	kernel.compute(0);
      else
	kernel.computeKernel3D();
    }
  else
    {
      cout << "./computeKernel1D.x configFile" << endl;
      cout << "Compute the Born sensitivity kernel using the decompostion in spherical harmonics from (Fournier et al. 2018)" << endl;
      abort();
    }
  //return FinalizeMontjoie();

#ifdef SELDON_WITH_MPI
    MPI_Finalize();
#endif

  return 0;
}
