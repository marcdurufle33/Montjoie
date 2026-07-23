#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Montjoie.hxx"

using namespace Montjoie;

class InputVar : public InputDataProblem_Base
{
public:
  Vector<string> mesh_data;
  Real_wp omega, tau;
  int order;

  InputVar()
  {
    omega = 2.0*pi_wp; tau = 0.0;
    order = 1;
  }

  void SetInputData(const string& description_field, const VectString& parameters)
  {
    if (description_field == "FileMesh")
      {
	mesh_data = parameters;
      }
    else if (description_field == "Frequency")
      {
	omega = 2.0*pi_wp*to_num<Real_wp>(parameters(0));
      }
    else if (description_field == "OrderDiscretization")
      {
	order = to_num<int>(parameters(0));
      }
    else if (description_field == "PenalizationDG")
      {
	tau = to_num<Real_wp>(parameters(0));
      }

  }
  
};

void RunAll(const string& name_element, const string& name_equation, const string& input_file)
{
  Mesh<Dimension2> mesh;
  InputVar var;

  ReadInputFile(input_file, var);

  mesh.SetGeometryOrder(var.order);
  mesh.ConstructMesh(Mesh<Dimension2>::QUADRILATERAL_MESH, var.mesh_data);
  
  QuadrangleLobatto Fb;
  Fb.ConstructFiniteElement(var.order);
  
  int nb_dof = (var.order+1)*mesh.GetNbBoundary();
  int Nvol = mesh.GetNbElt()*Fb.GetNbDof();
  nb_dof += 3*Nvol;
  Matrix<Complex_wp, General, ArrayRowSparse> A;

  //for (int i = 0; i < mesh.GetNbBoundary(); i++)
  //{ DISP(i); DISP(mesh.Boundary(i)); }
  
  A.Reallocate(nb_dof, nb_dof);
  VectR2 s; int offset = 0;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      mesh.GetVerticesElement(i, s);
      
      Real_wp h = s(1)(0) - s(0)(0);
      
      IVect num_u(Fb.GetNbDof()), num_vx(Fb.GetNbDof()), num_vy(Fb.GetNbDof());
      for (int j = 0; j < Fb.GetNbDof(); j++)
	{
	  num_u(j) = offset + j;
	  num_vx(j) = Nvol + offset + j;
	  num_vy(j) = 2*Nvol + offset + j;
	}
      
      // mass part
      for (int j = 0; j < Fb.GetNbDof(); j++)
	{
	  Real_wp poids = Fb.WeightsND(j)*h*h;
	  Complex_wp vloc = -Iwp*var.omega*poids;
	  A.AddInteraction(num_u(j), num_u(j), vloc);
	  A.AddInteraction(num_vx(j), num_vx(j), vloc);
	  A.AddInteraction(num_vy(j), num_vy(j), vloc);
	}

      // stiffness part
      VectReal_wp phi(Fb.GetNbDof()), grad_phi(2*Fb.GetNbDof());
      for (int j = 0; j < Fb.GetNbDof(); j++)
	{
	  phi.Zero(); phi(j) = 1.0;
	  Fb.ApplyRhTranspose(phi, grad_phi);
	  for (int k = 0; k < Fb.GetNbDof(); k++)
	    {
	      Complex_wp vloc = h*Fb.WeightsND(k)*grad_phi(2*k);
	      if (abs(vloc) > 1e-12)
		{
		  A.AddInteraction(num_u(j), num_vx(k), vloc);
		  A.AddInteraction(num_vx(k), num_u(j), -vloc);
		}

	      vloc = h*Fb.WeightsND(k)*grad_phi(2*k+1);
	      if (abs(vloc) > 1e-12)
		{
		  A.AddInteraction(num_u(j), num_vy(k), vloc);
		  A.AddInteraction(num_vy(k), num_u(j), -vloc);
		}
	    }
	}
      
      // flux part
      for (int num_loc = 0; num_loc < 4; num_loc++)
	{
	  R2 normale = Fb.NormaleLoc(num_loc);
	  IVect num_L(var.order+1);
	  int nf = mesh.Element(i).numBoundary(num_loc);
	  int offset_n = 3*Nvol + nf*(var.order+1);
	  if (mesh.Element(i).GetOrientationBoundary(num_loc))
	    for (int j = 0; j <= var.order; j++)
	      num_L(j) = offset_n + j;
	  else
	    for (int j = 0; j <= var.order; j++)
	      num_L(j) = offset_n + var.order-j;
	  
	  for (int j = 0; j <= var.order; j++)
	    {
	      Real_wp poids = h*Fb.Weights1D(j);
	      int num_dof = Fb.GetLocalNumber(num_loc, j);
	      
	      // part v.n phi et psi.n (u-lambda)
	      Complex_wp vloc = normale(0)*poids;
	      A.AddInteraction(num_u(num_dof), num_vx(num_dof), -vloc);
	      A.AddInteraction(num_vx(num_dof), num_u(num_dof), vloc);
	      A.AddInteraction(num_vx(num_dof), num_L(j), -vloc);
	      A.AddInteraction(num_L(j), num_vx(num_dof), vloc);
	      
	      vloc = normale(1)*poids;
	      A.AddInteraction(num_u(num_dof), num_vy(num_dof), -vloc);
	      A.AddInteraction(num_vy(num_dof), num_u(num_dof), vloc);
	      A.AddInteraction(num_vy(num_dof), num_L(j), -vloc);
	      A.AddInteraction(num_L(j), num_vy(num_dof), vloc);
	      
	      // part tau (u -lambda)
	      vloc = var.tau*poids;
	      A.AddInteraction(num_u(num_dof), num_u(num_dof), -vloc);
	      A.AddInteraction(num_u(num_dof), num_L(j), vloc);
	      A.AddInteraction(num_L(j), num_u(num_dof), vloc);
	      A.AddInteraction(num_L(j), num_L(j), -vloc);
	    }
	}
      
      offset += Fb.GetNbDof();
    }

  //A.WriteText("A.dat");
  
  // right hand side
  SetPoints<Dimension2> Pts;
  offset = 0;
  Real_wp alpha = -log(1e-6);
  Real_wp beta = sqrt(alpha / pi_wp);
  VectComplex_wp b(A.GetM()); b.Zero();
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      mesh.GetVerticesElement(i, s);
      mesh.FjElemNodal(s, Pts, mesh, i);
      
      Real_wp h = s(1)(0) - s(0)(0);
      for (int j = 0; j < Fb.GetNbDof(); j++)
	{
	  Real_wp r2 = AbsSquare(Pts.GetPointNodal(j));
	  Real_wp poids = h*h*Fb.WeightsND(j);
	  b(offset+j) = beta*exp(-alpha*r2)*poids;
	  //b(offset+j) = beta*exp(-alpha*r2);
	}

      offset += Fb.GetNbDof();
    }

  // solution
  All_MatrixLU<Complex_wp> mat_lu;
  mat_lu.SetPrintLevel(2);
  mat_lu.Factorize(A);

  mat_lu.Solve(b);

  GridInterpolation<Dimension2> grid;
  grid.SetXmin(mesh.GetXmin()); grid.SetXmax(mesh.GetXmax());
  grid.SetYmin(mesh.GetYmin()); grid.SetYmax(mesh.GetYmax());
  
  GridInterpolationFull<Dimension2> grid_plane;
  grid_plane.SetPlaneOutput(mesh.GetXmin(), mesh.GetXmax(),
			    mesh.GetYmin(), mesh.GetYmax(), 200, 200);
  
  IVect Index, NumPartie; VectReal_wp TetaPoints; VectR2 Points2D;
  grid_plane.GenerateGridPoints(grid, Index, Points2D, TetaPoints, NumPartie);

  // interpolation
  grid.InitInterpolationGrid(mesh);
  grid.Append(Points2D, TetaPoints, NumPartie, Index);   DISP(grid.GetNbGlobalPoints());
  grid.LocalizePoints(mesh);

  DISP(grid.GetNbGlobalPoints());
  VectComplex_wp val(grid.GetNbGlobalPoints());
  val.Zero();
  for (int i = 0; i < grid.GetNbGlobalPoints(); i++)
    {
      int ne = grid.GetElementNumber(i);
      if (ne >= 0)
	{
	  offset = ne*Fb.GetNbDof();
	  R2 pt_loc = grid.GetLocalCoordinate(i);
	  VectReal_wp phi;
	  Fb.ComputeValuesPhiRef(pt_loc, phi);
	  for (int j = 0; j < Fb.GetNbDof(); j++)
	    val(i) += phi(j)*b(offset + j);
	}
    }

  WriteMatlab(val, grid_plane, "test.dat", true, false);
  
}

int main(int argc, char** argv)
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

      RunAll(type_element, type_equation, file_name_data);
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"acous2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();  
}
