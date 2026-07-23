#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_WIRES

#include "Montjoie.hxx"

using namespace Montjoie;

class PocklingtonIntegral
{
public :
  
  // definition of the wire and associated mesh
  Wire<Dimension3> wire;
  Real_wp frequency, omega;
  R3 kwave;
  
  PocklingtonIntegral()
  {
    frequency = 1.0; omega = 2.0*pi_wp*frequency;
    kwave.Zero(); kwave(0) = omega;
  }
  
  void SetInputData(const VectString& parameters)
  {
    wire.SetInputData(parameters);
  }
  
  template<class Matrix1>
  void ComputeMatrixIntegral(Matrix1& mat_integral)
  {
    int nodl = wire.GetNbDof(), nb_dof_loc = wire.GetNbDof_Element(0);
    mat_integral.Reallocate(nodl, nodl);
    Real_wp eps2 = wire.GetEps()*wire.GetEps();
    MatrixSymFullComplex_wp mat_loc(nb_dof_loc, nb_dof_loc);
    Real_wp xi, xj, xj_p, nu_loc, xsi_sing, nup_loc, z_m_zp2;
    Complex_wp evalG; Real_wp phi, wn, r;
    
    int order = nb_dof_loc-1;
    Globatto<Real_wp> gauss_singZ, gauss_phi;
    gauss_singZ.ConstructQuadrature(2*order+1); int Np = gauss_singZ.nb_points_quadrature;
    VectReal_wp PointsZp(2*Np), WeightsZp(2*Np);
    
    gauss_phi.ConstructQuadrature(2*order+1); Np = gauss_phi.nb_points_quadrature;
    VectReal_wp PointsPhi(Np), WeightsPhi(Np);
    for (int i = 0; i < Np; i++)
      {
	xi = gauss_phi.Points(i);
	PointsPhi(i) = pi_wp*xi*xi;
	WeightsPhi(i) = xi*gauss_phi.Weights(i);
      }
    
    // double loop on edges
    for (int i1 = 0; i1 < wire.GetNbElements(); i1++)
      for (int i2 = 0; i2 < wire.GetNbElements(); i2++)
	{
	  // segments voisins ?
	  int a1 = wire.GetNumberVertex(i1, 0);
	  int b1 = wire.GetNumberVertex(i1, 1);
	  int a2 = wire.GetNumberVertex(i2, 0);
	  int b2 = wire.GetNumberVertex(i2, 1);
	  
	  mat_loc.Zero();
	  Real_wp ds1_ds2 = wire.GetLength(a1,b1)*wire.GetLength(a2,b2);
	  if ((a1==a2)||(a1==b2)||(b1==a2)||(b1=b2))
	    {
	      
	      // integration singuliere
	      for (int j = 0; j < gauss_singZ.nb_points_quadrature; j++)
		{
		  nu_loc = gauss_singZ.Points(j);
		  wire.Fj(i1, nu_loc, xj);
		  
		  // point singulier ?
		  if (i1==i2)
		    xsi_sing = nu_loc;
		  else if ((a2 == a1)||(a2==b1))
		    xsi_sing = 0.0;
		  else
		    xsi_sing = 1.0;
		  
		  for (int k = 0; k < gauss_singZ.nb_points_quadrature; k++)
		    {
		      xi = gauss_singZ.Points(k);
		      PointsZp(2*k) = xsi_sing*(1.0-xi*xi);
		      WeightsZp(2*k) = 2.0*xsi_sing*gauss_singZ.Weights(k)*xi;
		      PointsZp(2*k+1) = xsi_sing*(1.0-xi*xi) + xi*xi;
		      WeightsZp(2*k+1) = 2.0*(1.0-xsi_sing)*gauss_singZ.Weights(k)*xi;
		    }
		  
		  for (int jp = 0; jp < PointsZp.GetM(); jp++)
		    if (WeightsZp(jp) != 0)
		      {
			nup_loc = PointsZp(jp);
			wire.Fj(i2, nup_loc, xj_p);
			z_m_zp2 = square(xj-xj_p);
			evalG = 0.0;
			for (int k = 0; k < PointsPhi.GetM(); k++)
			  {
			    phi = PointsPhi(k);
			    wn = WeightsPhi(k);
			    r = sqrt(2.0*eps2*(1.0-cos(phi)) + z_m_zp2);
			    // la constante 1/(4pi) est dans les poids de phi
			    evalG += exp(Iwp*omega*r)/r*wn;
			  }
			evalG *= WeightsZp(jp)*gauss_singZ.Weights(j)*ds1_ds2;
			
			for (int m = 0; m < nb_dof_loc; m++)
			  for (int n = m; n < nb_dof_loc; n++)
			    mat_loc(m,n) += evalG*wire.GetValue_Phi(m, nu_loc)*wire.GetValue_Phi(n, nup_loc);
		      }
		}
	    }
	  else
	    {
	      // integration reguliere
	    }
	  
	  // assemblage
	  for (int m = 0; m < nb_dof_loc; m++)
	    for (int n = m; n < nb_dof_loc; n++)
	      mat_integral(wire.GetNumberDof(i1,m),wire.GetNumberDof(i2,n)) += mat_loc(m,n);
	}
  }
  
  void ComputeRightHandSide(VectComplex_wp& b)
  {
    Real_wp kortho = sqrt(kwave(0)*kwave(0) + kwave(1)*kwave(1));
    Real_wp kz = kwave(2);
    int nodl = wire.GetNbDof(), nb_dof_loc = wire.GetNbDof_Element(0); 
    b.Reallocate(nodl); b.Zero();
    
    Globatto<Real_wp> gauss;
    int order = nb_dof_loc-1; gauss.ConstructQuadrature(order);
    
    MatrixFullReal_wp ValuePhi(nb_dof_loc, gauss.nb_points_quadrature);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < gauss.nb_points_quadrature; j++)
	ValuePhi(i,j) = wire.GetValue_Phi(i, gauss.Points(j));
    
    for (int i1 = 0; i1 < wire.GetNbElements(); i1++)
      {
	int a1 = wire.GetNumberVertex(i1, 0);
	int b1 = wire.GetNumberVertex(i1, 1);
	Real_wp dsj = wire.GetLength(a1,b1);
	for (int j = 0; j < gauss.nb_points_quadrature; j++)
	  {
	    Real_wp nu_loc = gauss.Points(j);
	    Real_wp pt_glob; wire.Fj(i1, nu_loc, pt_glob);
	    
	    Complex_wp u_inc = exp(Iwp*kz*pt_glob)*j0(kortho*wire.GetEps());
	    u_inc *= dsj*gauss.Weights(j);
	    for (int k = 0; k < nb_dof_loc; k++)
	      b(wire.GetNumberDof(i1,k)) -= u_inc*ValuePhi(k,j);
	  }
      }
  }
  
  void ComputeIntegralRepresentation(const VectComplex_wp& xsol, const VectR3& Vertices, VectComplex_wp& Urep)
  {
    int nb_dof_loc = wire.GetNbDof_Element(0);
    Real_wp eps = wire.GetEps(), xi, phi, R;
    Real_wp nu_loc, xj;
    
    int order = nb_dof_loc-1;
    Globatto<Real_wp> gaussZ, gauss_phi;
    gaussZ.ConstructQuadrature(2*order+1); int Np = gaussZ.nb_points_quadrature;
    
    gauss_phi.ConstructQuadrature(2*order+1); Np = gauss_phi.nb_points_quadrature;
    VectReal_wp PointsPhi(Np), WeightsPhi(Np);
    for (int i = 0; i < Np; i++)
      {
	xi = gauss_phi.Points(i);
	PointsPhi(i) = pi_wp*xi*xi;
	WeightsPhi(i) = xi*gauss_phi.Weights(i);
      }
    
    MatrixFullReal_wp ValuePhi(nb_dof_loc, gaussZ.nb_points_quadrature);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < gaussZ.nb_points_quadrature; j++)
	ValuePhi(i,j) = wire.GetValue_Phi(i, gaussZ.Points(j));
    
    int Npts = Vertices.GetM(); Urep.Reallocate(Npts); Urep.Fill(0);
    // loop on edges
    for (int i1 = 0; i1 < wire.GetNbElements(); i1++)
      {
	int a1 = wire.GetNumberVertex(i1, 0);
	int b1 = wire.GetNumberVertex(i1, 1);
	
	Real_wp dsj = wire.GetLength(a1,b1);
	
	// integration en z
	for (int j = 0; j < gaussZ.nb_points_quadrature; j++)
	  {
	    nu_loc = gaussZ.Points(j);
	    wire.Fj(i1, nu_loc, xj);
	    Real_wp poids = dsj*gaussZ.Weights(j);
	    Complex_wp Jcur = 0.0;
	    for (int n = 0; n < nb_dof_loc; n++)
	      Jcur += ValuePhi(n,j)*xsol(wire.GetNumberDof(i1,n));
	    
	    // loop on all points
	    for (int k = 0; k < Npts; k++)
	      {
		Real_wp z = Vertices(k)(2); Real_wp r = sqrt(square(Vertices(k)(0)) + square(Vertices(k)(1)));
		Real_wp dist = (z-xj)*(z-xj) + (r-eps)*(r-eps);
		Real_wp coef = 2.0*r*eps; Complex_wp evalG = 0.0;
		// integration en teta
		for (int m = 0; m < gauss_phi.nb_points_quadrature; m++)
		  {
		    phi = PointsPhi(m);
		    R = sqrt(dist+coef*(1.0-cos(phi)));
		    evalG += exp(Iwp*omega*R)/R*WeightsPhi(m);
		  }
		
		Urep(k) += poids*Jcur*evalG;
	      }
	  }
      }
  }
  
  void WriteOutput_RegularGrid(const VectComplex_wp& U0, const GridInterpolationFull<Dimension3>& grid, const string& file_name)
  {
    VectComplex_wp Urep;
    ComputeIntegralRepresentation(U0, grid.var_interp->GlobalCoord, Urep);
    ofstream file_out(file_name.data()); int nb_true = 1;
    file_out.write(reinterpret_cast<char*>(const_cast<int*>(&nb_true)),sizeof(int));
    double df(0); grid.Write(file_out, df);
    Urep.Write(file_out);
  }
};

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  PocklingtonIntegral var;
  VectString parameters(10);
  parameters(0) = "0"; parameters(1) = "0"; parameters(2) = "-0.5";
  parameters(3) = "0"; parameters(4) = "0"; parameters(5) = "0.5";
  parameters(6) = "0.01"; parameters(7) = "7"; parameters(8) = "10"; parameters(9) = "10";
  
  var.SetInputData(parameters);
  MatrixSymFullComplex_wp mat_integral;
  // computation of the matrix
  var.ComputeMatrixIntegral(mat_integral);
  // LU factorization
  IVect pivot(mat_integral.GetM());
  GetLU(mat_integral, pivot);
  
  // computation of the rhs
  VectComplex_wp b_source;
  var.ComputeRightHandSide(b_source);
  // we solve linear system
  SolveLU(mat_integral, pivot, b_source);
  
  // grid
  int n = 100;
  GridInterpolation<Dimension3> pts_grid; GridInterpolationFull<Dimension3> grid;
  grid.type_output_file = grid.THREE_PLANES;
  grid.nbPoints_x = n; grid.nbPoints_y = n; grid.nbPoints_z = n;
  grid.center.Zero();
  grid.xmin = -1.0; grid.ymin = -1.0; grid.zmin = -1.0;
  grid.xmax = 1.0; grid.ymax = 1.0; grid.zmax = 1.0;
  grid.InitGrid(pts_grid);
  
  // ecriture sur un fichier de sortie
  string file_name = "Uref.dat";
  var.WriteOutput_RegularGrid(b_source, grid, file_name);
  
  return FinalizeMontjoie();
}
