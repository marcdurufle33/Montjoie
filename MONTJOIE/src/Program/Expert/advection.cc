#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Harmonic/MontjoieHarmonic.hxx"

using namespace Montjoie;

namespace Montjoie
{

  template<class Dim>
  class AdvectionEquationRD : public GenericEquation<Real_wp, 1>
  {
  public :
    typedef Dim Dimension;
    
    enum {nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_vec = 0,
	  nb_components_en = 1, nb_components_hn = 1};
    
    static bool ComputeDFjm1() { return true; }
        
    template<class TypeEquation>
    static void ApplyBoundaryFlux(const EllipticProblem<TypeEquation>& var,
				  int i, int j, int ref, const typename Dimension::R_N& normale,
				  const Real_wp& val_u, Real_wp& flux)
    {
      // DISP(var.ref_lambda(ref).GetCoefficient(var, i, j, var.PointsQuadrature(i)(j)));
      flux = DotProd(normale, var.ref_lambda(ref).GetCoefficient(var, i, j, var.PointsQuadrature(i)(j)) )*val_u;
      // DISP(var.PointsQuadrature(i)(j)); DISP(var.ref_lambda(ref).GetCoefficient(var, i, j, var.PointsQuadrature(i)(j)));
    }
    
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var, int i,
				  const ElementReference<Dimension, 1>& Fb)
    {
      abort();
      /* int Np = Fb.GetNbPointsQuadrature();
      var.PointsQuadrature(i).Reallocate(Np);
      for (int j = 0; j < Np; j++)
      var.PointsQuadrature(i)(j) = Points.GetPointQuadrature(j);*/
    }

  };
  
  template<int d>
  class VarPhysicalVelocity
  {
  public :
    template<class TypeEquation>
    TinyVector<Real_wp, d> GetCoefficient(const EllipticProblem<TypeEquation>& var,
					  int i, int j, const TinyVector<Real_wp, d>& point) const
    {
      // return R2(1, 0);
      return R2(-point(1), point(0));
    } 
    
    Real_wp GetMaximumNorm()
    {
      return Real_wp(1);
    }
    
  };

  class VarPhysicalViscosity
  {
  public :
    template<class TypeEquation, int d>
    Real_wp GetCoefficient(const EllipticProblem<TypeEquation>& var,
			   int i, int j, const TinyVector<Real_wp, d>& point) const
    {
      return Real_wp(0);
    } 
    
  };
  

  template<class TypeEquation>
  class VarAdvectionRD_Base : public VarHarmonic<TypeEquation>
  {
  public :
    typedef typename TypeEquation::Dimension Dimension;
    
    Vector<VarPhysicalVelocity<Dimension::dim_N> > ref_lambda;
    Vector<VarPhysicalViscosity> ref_nu;
    Vector<Vector<typename Dimension::R_N> > PointsQuadrature;
    
    VarAdvectionRD_Base() : VarHarmonic<TypeEquation>()
    {
      this->delta_penalization = 0.1;
    }
    
    int GetNbPhysicalIndices() const { return ref_lambda.GetM(); }
    void InitIndices(int n)
    {
      ref_lambda.Reallocate(n);
    }
    
    void SetIndices(int i, const VectString& parameters)
    {
    }

    void SetPhysicalIndex(const string&, int, const Vector<string>&) {}
    string GetPhysicalIndexName(int m) const { return string(); }
    bool IsVaryingMedia(int ref) const { return true; }
    bool IsVaryingMedia(int m, int ref) const { return true; }

    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
			   IVect& num_ref, IVect& num_index, IVect& num_component,
			   Vector<bool>& compute_grad, Vector<bool>& compute_hess)
    {
    }
    
    void AllocateMassMatrices()
    {
      PointsQuadrature.Reallocate(this->mesh.GetNbElt());
    }
  };
  
  template<class Dimension>
  class EllipticProblem<AdvectionEquationRD<Dimension> > 
    : public VarAdvectionRD_Base<AdvectionEquationRD<Dimension> > 
  {
    
  };


  template<class TypeEquation>
  class AdvectionSourceDirichlet
    : public VirtualSourceFEM<Real_wp, typename TypeEquation::Dimension>
  {
  public :
    typedef typename TypeEquation::Dimension Dimension;
    typedef typename Dimension::R_N R_N;
    
    enum {nb_unknowns = 1, nb_components = 1};
    
    AdvectionSourceDirichlet(const EllipticProblem<TypeEquation>& var)
      : VirtualSourceFEM<Real_wp, typename TypeEquation::Dimension>(var) {}
    
    void InitElement(int num_elem, const typename Dimension::VectR_N& s) {}
    
    void InitSurface(int i, int num_boundary, int num_elem, int num_loc)
    {
    }
    
    bool IsNonNullDirichletSource(int cond)
    {
      // DISP(cond); DISP( (EllipticProblem<TypeElement, TypeEquation>::LINE_DIRICHLET) );
      if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
	return true;
      
      return false;
    }
    
    
    void EvaluateFunction(int num_elem, int j, const R_N& point, VectReal_wp& f)
    {
      Real_wp x = point(0), y = point(1);
      // f(0) = 1.0-y*y;
      f(0) = x*(1.0-x);
      // f(0) = 1.0;
      // DISP(y); DISP(f);
      // if (abs(y) <= 1e-10)
      // f(0) = sqrt(x*x+y*y);
    } 
    
  };
  
}

template<class TypeEquation>
void RunAllRD(EllipticProblem<TypeEquation>& var, const string& file_name,
	      const string& name_element)
{
  typedef typename TypeEquation::Dimension Dimension;
  var.InitIndices(PhysicalConstant::nb_max_indices);
  ReadInputFile(file_name, var);
  
  var.ComputeMeshAndFiniteElement(name_element);
  
  var.TreatDirichletCondition();
  var.ComputeMassMatrix();
  
  // Uh -> solution
  int N = var.GetNbDof();
  VectReal_wp Uh(N), AhU(N), b_rhs(N);
  Uh.Fill(0);
  
  // inhomogeneous dirichlet condition
  AdvectionSourceDirichlet<TypeEquation> fsrc(var);
  b_rhs.Fill(0);
  //var.SetSurfacicProjection(b_rhs, fsrc);
  abort();
  //var.InitDirichletCondition(b_rhs);  
  //var.ImposeDirichletCondition(Uh);
  
  Real_wp delta = 0.1, teta = 0.5;
  Mesh<Dimension>& mesh = var.mesh;
  mesh.Write("test.mesh");
  
  // evaluation of alpha for each triangle
  VectReal_wp alpha(mesh.GetNbElt());
  Real_wp max_grad(0); abort();
  //Real_wp max_grad = var.GetReferenceElement(0).GetMaxIntegralGradient(); // DISP(max_grad);
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      Real_wp hmin = 100;
      int ref = mesh.Element(i).GetReference();
      for (int j = 0; j < mesh.Element(i).GetNbEdges(); j++)
	{
	  int ne = mesh.Element(i).numEdge(j);
	  int n1 = mesh.GetEdge(ne).numVertex(0);
	  int n2 = mesh.GetEdge(ne).numVertex(1);
	  hmin = min(hmin, mesh.Vertex(n1).Distance(mesh.Vertex(n2)));
	}
      
      alpha(i) = 5.0*hmin*var.ref_lambda(ref).GetMaximumNorm()*max_grad;
    }
  // DISP(alpha);

  typename Dimension::R_N normale, lambda, df_lambda;
  int num_dof;
  VectReal_wp Uloc, Uboundary, phiL, penalU, dU; Real_wp flux(0), beta, phiH, vloc;
  // loop to obtain the solution
  bool test_loop = true; int nb_iter = 1;
  // we inject the solution to check residuals
  // Uh.Fill(0);
  //var.AddVolumicProjection(1.0, Uh, fsrc);
  DISP(Uh);
  // var.WriteDatas(Uh);
  while (test_loop)
    {
      // vector AhU will contain \sum phi_i
      AhU.Fill(0);
      
      // loop on elements       
      for (int num_elem = 0; num_elem < mesh.GetNbElt(); num_elem++)
	{
	  const ElementReference<Dimension, 1>& Fb = var.GetReferenceElement(0);
	  int nb_dof_loc = Fb.GetNbDof();
	  int nb_quad_boundaries = Fb.GetNbPointsQuadBoundaries();
	  int ref = mesh.Element(num_elem).GetReference();
	  int nb_quad = Fb.GetNbPointsQuadratureInside();
	  Uloc.Reallocate(nb_dof_loc);
	  penalU.Reallocate(nb_dof_loc);
	  dU.Reallocate(nb_quad*Dimension::dim_N);
	  phiL.Reallocate(nb_dof_loc);
	  Uboundary.Reallocate(nb_quad_boundaries);
	  
	  // we get u on the element
	  for (int j = 0; j < Fb.GetNbDof(); j++)
	    {
	      num_dof = var.mesh_num.Element(num_elem).GetNumberDof(j);
	      Uloc(j) = Uh(num_dof);
	    }
	  
	  // we compute values on the boundary
	  abort();
	  //Fb.ApplySh(Uloc, Uboundary);
	  
	  // DISP(num_elem); DISP(Uloc); DISP(Uboundary);
	  // we compute \int_{\partial T} lambda n u ds
	  Real_wp residual = 0;
	  for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	    {
	      int num_face = mesh.Element(num_elem).numBoundary(num_loc);
	      bool new_face = var.IsNewFace(num_elem)(num_loc);
	      int nb_quad_face = Fb.GetNbQuadBoundary(num_loc);
	      int offset_face = Fb.GetOffsetSh(num_loc);
	      // DISP(nb_quad_face);
	      for (int k = 0; k < nb_quad_face; k++)
		{
		  normale = var.Glob_normale(num_face)(k);
		  if (!new_face)
		    Mlt(-1.0, normale);
		  
		  // DISP(k); DISP(normale);
		  TypeEquation::ApplyBoundaryFlux(var, num_elem, Fb.GetQuadNumber(num_loc, k),
						  ref, normale, Uboundary(offset_face+k), flux);
		  
		  // DISP(flux);
		  // DISP(k); DISP(Uboundary(offset_face+k)); DISP(normale); DISP(flux); DISP(Fb.GetWeightBoundary(k, num_loc));
		  residual += Fb.WeightsQuadratureBoundary(k, num_loc)*flux;
		}
	    }
	  
	  // stabilization term
	  Fb.ApplyRhTranspose(Uloc, dU);
	  
	  // DISP(num_elem); DISP(Uloc); DISP(dU); DISP(residual);
	  
	  for (int i = 0; i < nb_quad; i++)
	    {
	      lambda = var.ref_lambda(ref).GetCoefficient(var, num_elem, i, var.PointsQuadrature(num_elem)(i));
	      if (var.Glob_DFjm1(num_elem).GetM() > 1)
		Mlt(var.Glob_DFjm1(num_elem)(i), lambda, df_lambda);
	      else
		Mlt(var.Glob_DFjm1(num_elem)(0), lambda, df_lambda);
	      
	      vloc = 0;
	      //DISP(df_lambda);
	      // vloc = J_i DF_i^{-1} lambda \cdot grad(u)
	      for (int k = 0; k < Dimension::dim_N; k++)
		vloc += dU(Dimension::dim_N*i + k)*df_lambda(k);
	      
	      if (var.Glob_jacobian(num_elem).GetM() > 1)
		vloc *= teta*Fb.WeightsND(i)*Fb.WeightsND(i)/var.Glob_jacobian(num_elem)(i);
	      else
		vloc *= teta*Fb.WeightsND(i)/var.Glob_jacobian(num_elem)(0);
	      //DISP(vloc);
	      for (int k = 0; k < Dimension::dim_N; k++)
		dU(Dimension::dim_N*i + k) = df_lambda(k)*vloc;
	    }
	  
	  // DISP(dU);
	  Fb.ApplyRh(dU, penalU);
	  // penalU.Fill(0);
	  
	  // DISP(residual);
	  if (abs(residual) > 1e-10)
	    {	    
	      // computation of phi_i^L with Friedrichs scheme
	      Real_wp sum_phiLp = 0;
	      for (int i = 0; i < nb_dof_loc; i++)
		{
		  phiL(i) = 0;
		  for (int j = 0; j < nb_dof_loc; j++)
		    phiL(i) += Uloc(i)-Uloc(j);
		  
		  // DISP(i); DISP(phiL(i));
		  phiL(i) = (residual + alpha(num_elem)*phiL(i))/nb_dof_loc;
		  phiL(i) = max(Real_wp(0), phiL(i)/residual);
		  sum_phiLp += phiL(i);
		}
	      // DISP(phiL); DISP(sum_phiLp);
	      
	      // computation of phi_i^H
	      for (int i = 0; i < nb_dof_loc; i++)
		{
		  if (abs(sum_phiLp) > 1e-10)
		    beta = phiL(i)/sum_phiLp;
		  else
		    beta = 0;
		  phiH = beta*(residual+ penalU(i));
		  // DISP(i); DISP(beta); DISP(phiH);
		  num_dof = var.mesh_num.Element(num_elem).GetNumberDof(i);
		  AhU(num_dof) += phiH;
		}
	    }
	  	  
	}
      
      // DISP(AhU);
      // DISP(Norm2(AhU));
      // u_i^{n+1} = u_i^n - delta \sum phi_i^T
      Add(-delta, AhU, Uh);
      
      // DISP(delta);
      // var.WriteDatas(Uh);
      
      // int test_input; cout<<"we wait"<<endl; cin>>test_input;
      
      DISP(Norm2(AhU));
      // Dirichlet condition
      //var.ImposeDirichletCondition(Uh);
      abort();
      
      Real_wp erreur = Norm2(AhU)/Norm2(Uh);
      if ((nb_iter > 10000) || (erreur <= 1e-6))
	test_loop = false;
      
      nb_iter++;
    }
  
  DISP(nb_iter);
  DISP(Norm2(AhU));
  var.WriteDatas(Uh);
  
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  string file_name_data(argv[1]);
  
  string type_element, type_equation;
  getElement_Equation(file_name_data, type_element, type_equation);
  
  EllipticProblem<AdvectionEquationRD<Dimension2> > vars;
  RunAllRD(vars, file_name_data, type_element);
  
  return FinalizeMontjoie();
}
