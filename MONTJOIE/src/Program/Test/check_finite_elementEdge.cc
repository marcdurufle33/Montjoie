#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_HP_HCURL

#include "Harmonic/MontjoieHarmonic.hxx"

#ifdef MONTJOIE_WITH_FAST_LIBRARY
// includes for fast compilation
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#endif

#define CHECK_ABORT abort()
// #define CHECK_ABORT cout<<"ERROR"<<endl

using namespace Montjoie;

template<class T>
bool EqualVector(const Vector<T>& x, const Vector<T>& y)
{
  if (x.GetM() != y.GetM())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    if (abs(x(i) - y(i)) > 1e-12)
      return false;
  
  return true;
}

Real_wp FuncPowS(const R2& x, int m, int n)
{
  return pow(x(0)+0.5*x(1), m)*pow(0.2+0.9*x(0)+x(1), n);
}


namespace Montjoie
{
  template<class Dim>
  class MyOwnEquation : public GenericEquation<Real_wp>
  {
  public :
    typedef Dim Dimension;
    enum{nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_vec = 0,
	 nb_components_en = 1, nb_components_hn = 1, type_element = 2, nb_unknowns_hdg = 0};
  };
  
  
  template<class Dimension>
  class EllipticProblem<MyOwnEquation<Dimension> > : public VarHarmonic<MyOwnEquation<Dimension> >
  {
  public :
    typedef typename Dimension::VectR_N VectR_N;
    
    virtual void InitIndices(int n) {}
    virtual int GetNbPhysicalIndices() const { return 0; }
    virtual void SetIndices(int ref, const Vector<string>&) {}
    virtual void SetPhysicalIndex(const string&, int, const Vector<string>&) {}
    virtual bool IsVaryingMedia(int iquad) const { return false; }
    virtual bool IsVaryingMedia(int m, int ref) const { return false; }
    virtual string GetPhysicalIndexName(int m) const { return string("rho"); }
    
    virtual void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
				   Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
				   IVect& num_ref, IVect& num_index, IVect& num_component,
				   Vector<bool>& compute_grad, Vector<bool>& compute_hess) {}
    
    virtual void ComputeStoreCoefficientsPML(int i1, int, const VectR_N& AllPoints) {}
    virtual void FinalizeComputationVaryingIndices() {}
    
    virtual void InitCyclicDomain() {}
    virtual void TreatDirichletCondition() {}
    virtual void InitVarGrid() {}
    virtual void ComputeVarGrid() {}
    
    virtual void ComputeArraySpecificEquation() {}
    virtual void FindElementsInsidePML() {}
    
    virtual bool IsSymmetricProblem(bool eigen = false) const { return true; }
    
    virtual void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Real_wp>&,
					 CondensationBlockSolver_Base<Real_wp>&,
					 const GlobalGenericMatrix<Real_wp>&) {}
    
    virtual void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Complex_wp>&,
					 CondensationBlockSolver_Base<Complex_wp>&,
					 const GlobalGenericMatrix<Complex_wp>&) {}
    
  };
  
} // end namespace Montjoie

void CheckSurfaceElement(const ElementReference<Dimension2, 2>& Fb)
{
  const ElementReference<Dimension1, 1>& edge = Fb.GetSurfaceFiniteElement(0);
  int r = Fb.GetOrder();
  if ( (edge.GetGeometryOrder() != r) || 
       (edge.GetNbPointsQuadratureInside() != Fb.GetNbQuadBoundary(0)) || 
       (edge.GetNbDof() != Fb.GetNbDofBoundary(0)) )
    {
      cout << " GetSurfaceFiniteElement incorrect " << endl;
      CHECK_ABORT;
    }
  
  if (!EqualVector(edge.Points(), Fb.Points1D()))
    {
      cout << " GetSurfaceFiniteElement incorrect (1-D points different) " << endl;
      DISP(edge.Points()); DISP(Fb.Points1D());
      CHECK_ABORT;
    }
  
  if (!EqualVector(edge.Weights(), Fb.Weights1D()))
    {
      cout << " GetSurfaceFiniteElement incorrect (1-D weights different) " << endl;
      DISP(edge.Weights()); DISP(Fb.Weights1D());
      CHECK_ABORT;
    }

  if (!EqualVector(edge.PointsDof(), Fb.PointsDofBoundary(0)))
    {
      cout << " GetSurfaceFiniteElement incorrect (1-D dof points different) " << endl;
      DISP(edge.PointsDof()); DISP(Fb.PointsDofBoundary(0));
      //CHECK_ABORT;
    }  

  R2 pt_loc; VectReal_wp phi_surf;
  VectR2 phi;
  for (int j = 0; j <= r; j++)
    {
      Fb.GetLocalCoordOnBoundary(0, edge.Points(j), pt_loc);
      edge.ComputeValuesPhiRef(edge.Points(j), phi_surf);
      Fb.ComputeValuesPhiRef(pt_loc, phi);

      for (int i = 0; i < edge.GetNbDof(); i++)
	{
	  int n = Fb.GetLocalNumber(0, i);
	  if (abs(phi_surf(i)- phi(n)(0)) > 1e-12)
	    {
	      DISP(pt_loc); DISP(i); DISP(n); DISP(phi_surf(i)); DISP(phi(n)(0));
	      cout << "Value of basis functions differ in GetSurfaceFiniteElement" << endl;
	      abort();
	    }
	}
    }
}

void CheckHcurlFiniteElement(ElementReference<Dimension2, 2>& Fb,
			     const string& name_elt, Mesh<Dimension2>& mesh, VectR2& s,
			     bool check_surf = true, bool check_stiff = true)
{
  // checking Hcurl finite elements
  Real_wp aire = 1.0;
  if (s.GetM() == 3)
    aire = 0.5;
  
  Real_wp threshold = 300.0*epsilon_machine;
  MeshNumbering<Dimension2> mesh_num(mesh);
  R2::threshold = threshold;
  
  int num_loc = 0;
  // main loop on orders of approximations
  for (int r = 1; r <= 4; r++)
    {      
      cout << "Checking order " << r << endl;
      // we can use a different order for geometry
      int rgeom = min(r, 3);
      int rquad = r;
      int type_quad = -1;
      Fb.ConstructFiniteElement(r, rgeom, rquad, type_quad);
      Fb.SetMesh(mesh);
      
      if (Fb.GetOrder() != r)
        {
          cout << "Problem when setting order of finite element " << endl;
          CHECK_ABORT;
        }
      
      if (Fb.GetGeometryOrder() != rgeom)
        {
          cout << "Problem when setting geometry order of finite element " << endl;
          CHECK_ABORT;
        }
      
      int Nquad = Fb.GetNbPointsQuadratureInside();
      int nb_dof = 0, nb_quad = 0;
      for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
        {
          if (Fb.GetOffsetSh(num_loc) != nb_quad)
            {
              cout << "number of quadrature points = " << nb_quad << " and offset = " << Fb.GetOffsetSh(num_loc) << endl;
              cout<<"Warning : offset for Sh array may have problems " <<endl;
            }
          
	  if (Fb.GetNbQuadBoundary(num_loc) == 0)
	    {
              cout << "No quadrature points on boundary " << num_loc << endl;
	      CHECK_ABORT;
	    }
	  
          nb_quad += Fb.GetNbQuadBoundary(num_loc);
          nb_dof += Fb.GetNbDofBoundary(num_loc);          
        }
      
      if (nb_quad != Fb.GetNbPointsQuadBoundaries())
        {
          cout << "Number of quadrature points on boundaries not coherent " << endl;
          CHECK_ABORT;
        }
      
      if (nb_dof < Fb.GetNbDofBoundaries())
	{
          cout << "Number of dofs on boundaries not coherent " << endl;
	  CHECK_ABORT;
	}
      
      nb_dof = Fb.GetNbDof();
      if (Fb.GetNbDof() == 0)
        {
          cout << "Element with no degree of freedom " << endl;
          CHECK_ABORT;
        }      
      
      if (Fb.GetNbPointsUsedForSh() != nb_quad)
        cout<<"Warning : number of points different for Sh and boundary integrals " << endl;
      
      // checking nodal points
      VectReal_wp xi, omega;
      ComputeGaussLobatto(xi, omega, rgeom);
      VectReal_wp points1d = Fb.PointsNodal1D();
      for (int k = 0; k <= rgeom; k++)
        if (abs(points1d(k) - xi(k)) > threshold)
          {
            cout << "1-D nodal points not on Lobatto " << endl;
            CHECK_ABORT;
          }
      
      Vector<R2> points_nodal = Fb.PointsNodalND();
      for (int k = 0; k < Fb.GetNbVertices(); k++)
        if (points_nodal(k) != s(k))
          {
            cout << "first nodal points have to be vertices " << endl;
            CHECK_ABORT;
          }
      
      R2 ptA, ptB;
      
      int nb_nodes_bounds = Fb.GetNbVertices()*rgeom;
      
      for (int k = nb_nodes_bounds; k < Fb.GetNbPointsNodalElt(); k++)
        if (Fb.GetDistanceToBoundary(points_nodal(k)) <= 1e-12)
          {
            cout << "Other nodal points have to be inside element " << endl;
            CHECK_ABORT;
          }
      
      // checking quadrature points
      for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
        {
          for (int k = 0; k < Fb.GetNbQuadBoundary(num_loc); k++)
            {
              Real_wp tloc = Fb.PointsQuadratureBoundary(k, num_loc);
              Fb.GetLocalCoordOnBoundary(num_loc, tloc, ptA);
              if (ptA != Fb.PointsND(Fb.GetQuadNumber(num_loc, k)))
                {
                  cout << "Problems on quadrature points of boundaries " << endl;
                  CHECK_ABORT;
                }
            }
        }
      
      Real_wp aire_num = 0;
      for (int k = 0; k < Nquad; k++)
        aire_num += Fb.WeightsND(k);
      
      if (abs(aire - aire_num) > 1e-12)
        {
          cout << "Sum of weights is not equal to the volume of element " << endl;
          CHECK_ABORT;
        }

      Real_wp int_x2 = 0;
      points1d = Fb.Points1D();
      VectReal_wp weights1d = Fb.Weights1D();
      for (int k = 0; k < points1d.GetM(); k++)
        int_x2 += weights1d(k)*square(points1d(k));      
      
      if (rquad > 1)
	if (abs(int_x2 - Real_wp(1)/3) > 1e-12)
	  {
            cout << "1-D quadrature points not correct " << endl;
	    CHECK_ABORT;
	  }
      
      // checking if equality between points have been correctly notified
      bool quad_equal_nodal = true;
      if (Fb.GetNbPointsQuadrature() == Fb.GetNbPointsNodalElt())
	{
	  for (int i = 0; i < Fb.GetNbPointsQuadrature(); i++)
	    if (Fb.PointsND(i) != Fb.PointsNodalND(i))
	      quad_equal_nodal = false;
	}
      else
	quad_equal_nodal = false;
      
      if (quad_equal_nodal)
	{
	  if (!Fb.QuadratureEqualNodal())
	    cout <<"Warning : quadrature and nodal points are equal, but it is not noticed " << endl;
	}
      else
	{
	  if (Fb.QuadratureEqualNodal())
	    {
              cout << "quadrature points = nodal points, but it is not the case " <<endl;
	      CHECK_ABORT;
	    }
	}
      
      bool dof_equal_nodal = true;
      if (Fb.GetNbPointsDof() == Fb.GetNbPointsNodalElt())
	{
	  for (int i = 0; i < Fb.GetNbPointsDof(); i++)
	    if (Fb.PointsDofND(i) != Fb.PointsNodalND(i))
	      dof_equal_nodal = false;
	}
      else
	dof_equal_nodal = false;
      
      if (dof_equal_nodal)
	{
	  if (!Fb.DofEqualNodal())
	    cout <<"Warning : dof and nodal points are equal, but it is not noticed " << endl;
	}
      else
	{
	  if (Fb.DofEqualNodal())
	    {
              cout << "dof points = nodal points, but it is not the case " <<endl;
	      CHECK_ABORT;
	    }
	}
      
      bool dof_equal_quad = true;
      if (Fb.GetNbPointsQuadrature() == Fb.GetNbPointsDof())
	{
	  for (int i = 0; i < Fb.GetNbPointsQuadrature(); i++)
	    if (Fb.PointsND(i) != Fb.PointsDofND(i))
	      dof_equal_quad = false;
	}
      else
	dof_equal_quad = false;
      
      if (dof_equal_quad)
	{
	  if (!Fb.DofEqualQuadrature())
	    cout <<"Warning : quadrature and dof points are equal, but it is not noticed " << endl;
	}
      else
	{
	  if (Fb.DofEqualQuadrature())
	    {
              cout << "quadrature points = dof points, but it is not the case " <<endl;
	      CHECK_ABORT;
	    }
	}
      
      // checking if nodal shape functions satisfy delta_i,j
      VectReal_wp phi_node, phip, phim;
      for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
        {
          Fb.ComputeValuesPhiNodalRef(Fb.PointsNodalND(i), phi_node);
          for (int j = 0; j < Fb.GetNbPointsNodalElt(); j++)
            {
              Real_wp delta = 0;
              if (i == j)
                delta = 1.0;
              
              if (abs(delta - phi_node(j)) > threshold)
                {
                  cout << "nodal shape functions are not satisfying delta_ij " << endl;
                  CHECK_ABORT;
                }
            }
        }

      // checking gradient of nodal shape functions
      VectR2 grad_phi; Real_wp h = pow(epsilon_machine, 1.0/3);
      for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
        {
          Fb.ComputeGradientPhiNodalRef(Fb.PointsNodalND(i), grad_phi);
          Real_wp dphi_num; 
          for (int k = 0; k < 2; k++)
            {
              ptA = Fb.PointsNodalND(i);
              ptA(k) += h;
              Fb.ComputeValuesPhiNodalRef(ptA, phip);

              ptA(k) -= 2.0*h;
              Fb.ComputeValuesPhiNodalRef(ptA, phim);
              
              for (int j = 0; j < Fb.GetNbPointsNodalElt(); j++)
                {
                  dphi_num = (phip(j) - phim(j))/(2.0*h);
                  if (abs(dphi_num - grad_phi(j)(k)) > sqrt(threshold))
                    {
                      cout << "Gradient of nodal shape functions not correct " <<endl;
                      DISP(dphi_num); DISP(grad_phi(j)(k));
                      CHECK_ABORT;
                    }
                }
            }
        }

      // checking dofs      
      if (Fb.GetNbDofBoundaries() > Fb.GetNbDof())
        {
          cout << "Problem for dofs on boundaries, the are greater than total number of dofs " <<endl;
          CHECK_ABORT;
        }
            
      // checking transformation Fi for linear element
      mesh.Element(0).UnsetCurved();
      VectR2 vert_elt; R2 pt_loc, ptC; Matrix2_2 dfj, dfjm1;
      mesh.GetVerticesElement(0, vert_elt);
      SetPoints<Dimension2> PointsElem;
      for (int i = 0; i < s.GetM(); i++)
	{
	  Fb.Fj(vert_elt, PointsElem, s(i), ptA, mesh, 0);
	  int n = mesh.Element(0).numVertex(i);
	  if (mesh.Vertex(n) != ptA)
	    {
              cout << "tranformation Fi not correct " << endl;
	      CHECK_ABORT;
	    }
	  
	  Fb.DFj(vert_elt, PointsElem, s(i), dfj, mesh, 0);
	  for (int k = 0; k < 2; k++)
	    {
	      pt_loc = s(i); pt_loc(k) += h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptC, mesh, 0);
	      
	      pt_loc = s(i); pt_loc(k) -= h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptB, mesh, 0);
	      
	      for (int m = 0; m < 2; m++)
		{
		  Real_wp df_num = (ptC(m) - ptB(m))/(2.0*h);
		  if (abs(df_num - dfj(m, k)) > sqrt(threshold))
                    {
                      cout << "DFj not correctly computed " << endl;
                      CHECK_ABORT;
                    }
		}
	    }
	  
	  pt_loc.Fill(0);
	  FjInverseProblem<Dimension2> inverseFj(mesh, 0);
	  inverseFj.Solve(ptA, pt_loc);

	  if (pt_loc != s(i))
	    {
              cout << "Fj^-1 does not work " << endl;
	      CHECK_ABORT;
	    }
	}
      
      Fb.FjElemNodal(vert_elt, PointsElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
	{
	  Fb.FjLinear(vert_elt, Fb.PointsNodalND(i), ptA);
	  if (ptA != PointsElem.GetPointNodal(i))
	    {
              cout << "Fj not correct for nodal points " << endl;
	      CHECK_ABORT;
	    }
	}

      Fb.FjElemQuadrature(vert_elt, PointsElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsQuadrature(); i++)
	{
	  Fb.FjLinear(vert_elt, Fb.PointsND(i), ptA);
	  if (ptA != PointsElem.GetPointQuadrature(i))
	    {
              cout << "Fj not correct for quadrature points " << endl;
	      CHECK_ABORT;
	    }
	}

      Fb.FjElemDof(vert_elt, PointsElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsDof(); i++)
	{
	  Fb.FjLinear(vert_elt, Fb.PointsDofND(i), ptA);
	  if (ptA != PointsElem.GetPointDof(i))
	    {
              cout << "Fj not correct for dof points " << endl;
	      CHECK_ABORT;
	    }
	}
      
      
      SetMatrices<Dimension2> MatricesElem;
      Fb.DFjElemNodal(vert_elt, PointsElem, MatricesElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
	{
	  Fb.DFj(vert_elt, PointsElem, Fb.PointsNodalND(i), dfj, mesh, 0);
	  if (dfj != MatricesElem.GetPointNodal(i))
	    {
              cout << "DFj not correct for nodal points " << endl;
	      CHECK_ABORT;
	    }
	}

      Fb.DFjElemQuadrature(vert_elt, PointsElem, MatricesElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsQuadrature(); i++)
	{
	  Fb.DFj(vert_elt, PointsElem, Fb.PointsND(i), dfj, mesh, 0);
	  if (dfj != MatricesElem.GetPointQuadrature(i))
	    {
              cout << "DFj not correct for quadrature points " << endl;
	      CHECK_ABORT;
	    }
	}

      Fb.DFjElemDof(vert_elt, PointsElem, MatricesElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsDof(); i++)
	{
	  Fb.DFj(vert_elt, PointsElem, Fb.PointsDofND(i), dfj, mesh, 0);
	  if (dfj != MatricesElem.GetPointDof(i))
	    {
              cout << "DFj not correct for dof points " << endl;
	      CHECK_ABORT;
	    }
	}
      
      // and for curved element
      mesh.Element(0).SetCurved();
      
      // checking transformation Fi for curved element
      mesh.GetVerticesElement(0, vert_elt);
      Fb.FjElemNodal(vert_elt, PointsElem, mesh, 0); 
      for (int i = 0; i < s.GetM(); i++)
	{
	  Fb.Fj(vert_elt, PointsElem, s(i), ptA, mesh, 0);
	  int n = mesh.Element(0).numVertex(i);
	  if (mesh.Vertex(n) != ptA)
	    {
              cout << "Curved Fj does not coincide with vertices " << endl;
	      CHECK_ABORT;
	    }
	  
	  Fb.DFj(vert_elt, PointsElem, s(i), dfj, mesh, 0);
	  for (int k = 0; k < 2; k++)
	    {
	      pt_loc = s(i); pt_loc(k) += h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptC, mesh, 0);
	      
	      pt_loc = s(i); pt_loc(k) -= h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptB, mesh, 0);
	      
	      for (int m = 0; m < 2; m++)
		{
		  Real_wp df_num = (ptC(m) - ptB(m))/(2.0*h);
		  if (abs(df_num - dfj(m, k)) > sqrt(threshold))
		    {
                      cout << "Curved DFj not correct " << endl;
		      CHECK_ABORT;
		    }
		}
	    }
	  
	  pt_loc.Fill(0);
	  FjInverseProblem<Dimension2> inverseFj(mesh, 0);
	  inverseFj.Solve(ptA, pt_loc);
	  
	  if (pt_loc != s(i))
	    {
              //cout << "Curved Fj^-1 not correct " << endl;
              //R_N ptTest;
              //Fb.Fj(vert_elt, PointsElem, pt_loc, ptTest, mesh, 0);
              //DISP(ptTest); DISP(ptA); DISP(pt_loc); DISP(s(i));
	      //CHECK_ABORT;
	    }
	}
      
      Fb.FjElemNodal(vert_elt, PointsElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
	{
	  Fb.Fj(vert_elt, PointsElem, Fb.PointsNodalND(i), ptA, mesh, 0);
	  if (ptA != PointsElem.GetPointNodal(i))
	    {
              cout << "Curved Fj not correct on nodal points " << endl;
	      CHECK_ABORT;
	    }
	}

      Fb.FjElemQuadrature(vert_elt, PointsElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsQuadrature(); i++)
	{
	  Fb.Fj(vert_elt, PointsElem, Fb.PointsND(i), ptA, mesh, 0);
	  if (ptA != PointsElem.GetPointQuadrature(i))
	    {
              cout << "Curved Fj not correct on quadrature points " << endl;
	      CHECK_ABORT;
	    }
	}

      Fb.FjElemDof(vert_elt, PointsElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsDof(); i++)
	{
	  Fb.Fj(vert_elt, PointsElem, Fb.PointsDofND(i), ptA, mesh, 0);
	  if (ptA != PointsElem.GetPointDof(i))
	    {
              cout << "Curved Fj not correct on dof points " << endl;
	      CHECK_ABORT;
	    }
	}
      
      
      Fb.DFjElemNodal(vert_elt, PointsElem, MatricesElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
	{
	  Fb.DFj(vert_elt, PointsElem, Fb.PointsNodalND(i), dfj, mesh, 0);
	  if (dfj != MatricesElem.GetPointNodal(i))
	    {
              cout << "Curved DFj not correct on nodal points " << endl;
	      CHECK_ABORT;
	    }
	}

      Fb.DFjElemQuadrature(vert_elt, PointsElem, MatricesElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsQuadrature(); i++)
	{
	  Fb.DFj(vert_elt, PointsElem, Fb.PointsND(i), dfj, mesh, 0);
	  if (dfj != MatricesElem.GetPointQuadrature(i))
	    {
              cout << "Curved DFj not correct on quadrature points " << endl;
              DISP(dfj); DISP(i); DISP(MatricesElem.GetPointQuadrature(i));
	      CHECK_ABORT;
	    }
	}

      Fb.DFjElemDof(vert_elt, PointsElem, MatricesElem, mesh, 0);
      for (int i = 0; i < Fb.GetNbPointsDof(); i++)
	{
	  Fb.DFj(vert_elt, PointsElem, Fb.PointsDofND(i), dfj, mesh, 0);
	  if (dfj != MatricesElem.GetPointDof(i))
	    {
              cout << "Curved DFj not correct on dof points " << endl;
	      CHECK_ABORT;
	    }
	}
      
      num_loc = 1; Real_wp dsj; R2 normale_fj;
      Fb.FjSurfaceElem(vert_elt, PointsElem, mesh, 0, num_loc);
      Fb.DFjSurfaceElem(vert_elt, PointsElem, MatricesElem, mesh, 0, num_loc);
      for (int i = 0; i < Fb.GetNbQuadBoundary(num_loc); i++)
        {
          Fb.GetLocalCoordOnBoundary(num_loc, Fb.PointsQuadratureBoundary(i, num_loc), pt_loc);
          Fb.Fj(vert_elt, PointsElem, pt_loc, ptC, mesh, 0);
          
          if (ptC != PointsElem.GetPointQuadratureBoundary(i))
            {
              cout << "FjSurfaceElem not correct " << endl;
              CHECK_ABORT;
            }
          
          Fb.DFj(vert_elt, PointsElem, pt_loc, dfj, mesh, 0);
          if (dfj != MatricesElem.GetPointQuadratureBoundary(i))
            {
              cout << "DFjSurfaceElem not correct " << endl;
              CHECK_ABORT;
            }
          
          GetInverse(dfj, dfjm1);
          Fb.GetNormale(dfjm1, normale_fj, dsj, num_loc);
          if (normale_fj != MatricesElem.GetNormaleQuadratureBoundary(i))
            {
              cout << "problem when computing normales " << endl;
              CHECK_ABORT;
            }
          
          if (abs(dsj - MatricesElem.GetDsQuadratureBoundary(i)) > 1e-12)
            {
              cout << "problem when computing surface element " << endl;
              CHECK_ABORT;
            }
        }
      
      // checking continuity 
      {
	EllipticProblem<MyOwnEquation<Dimension2> > var;
        MeshNumbering<Dimension2>& mesh_num = var.GetMeshNumbering(0);
        
        int type_mesh = Mesh<Dimension2>::QUADRILATERAL_MESH;
        if (s.GetM() == 3)
          type_mesh = Mesh<Dimension2>::TRIANGULAR_MESH;
        
        mesh_num.SetOrder(r); var.mesh.SetGeometryOrder(rgeom);
        TinyVector<int, 4> ref_bound; ref_bound.Fill(1);
        var.mesh.CreateRegularMesh(R2(-2,-2), R2(2,2), TinyVector<int, 2>(3,3), 1,
                                   ref_bound, type_mesh);
        
        var.mesh.ClearConnectivity();
        // random permutation of local vertex numbers
        srand(time(NULL));
        //srand(0);
        for (int i = 0; i < var.mesh.GetNbElt(); i++)
          {
            int nb_vert = var.mesh.Element(i).GetNbVertices();
            int ref = var.mesh.Element(i).GetReference();
            IVect num(nb_vert), numb;
            for (int j = 0; j < nb_vert; j++)
              num(j) = var.mesh.Element(i).numVertex(j);
            
            numb = num;            
            
            switch (nb_vert)
              {
              case 3 :
                {
                  int perm = rand()%3;
                  switch (perm)
                    {
                    case 0 :
                      num = numb; break;
                    case 1 :
                      num(0) = numb(1); num(1) = numb(2); num(2) = numb(0); break;
                    case 2 :
                      num(0) = numb(2); num(1) = numb(0); num(2) = numb(1); break;
                    }
                }
                break;
              case 4 :
                {
                  int perm = rand()%4;
                  switch (perm)
                    {
                    case 0 :
                      num = numb; break;
                    case 1 :
                      num(0) = numb(1); num(1) = numb(2); num(2) = numb(3); num(3) = numb(0); break;
                    case 2 :
                      num(0) = numb(2); num(1) = numb(3); num(2) = numb(0); num(3) = numb(1); break;
                    case 3 :
                      num(0) = numb(3); num(1) = numb(0); num(2) = numb(1); num(3) = numb(2); break;
                    }
                }
                break;
              }
            
            var.mesh.Element(i).Init(num, ref);
          }
        
        var.mesh.ReorientElements();
        var.mesh.FindConnectivity();

	mesh_num.SetOrder(r);
	var.ConstructFiniteElement(name_elt);
        
	mesh_num.NumberMesh();
        var.mesh.Write("test.mesh");
        
        var.CheckContinuity();
      }
      
      cout << "Check of continuity successful" << endl;

      // checking surface finite elements
      if (check_surf)
	CheckSurfaceElement(Fb);
      
      // checking mass matrix
      Matrix<R2> ValuePhi(Fb.GetNbDof(), Fb.GetNbPointsQuadrature());
      Real_wp max_cij = 0; VectR2 phi; VectReal_wp curl_phi;
      for (int k = 0; k < Fb.GetNbPointsQuadrature(); k++)
        {
          Fb.ComputeValuesPhiRef(Fb.PointsND(k), phi);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              ValuePhi(i, k) = phi(i);
              max_cij = max(Norm2(ValuePhi(i, k)), max_cij);
            }
        }
      
      Matrix<Real_wp, Symmetric, RowSymPacked> Mh(Fb.GetNbDof(), Fb.GetNbDof()), InvMh; Mh.Fill(0);
      Real_wp max_mij = 0;
      for (int i = 0; i < Fb.GetNbDof(); i++)
        for (int j = i; j < Fb.GetNbDof(); j++)
          {
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              Mh(i, j) += Fb.WeightsND(k)*DotProd(ValuePhi(i, k), ValuePhi(j, k));
                        
            max_mij = max(max_mij, abs(Mh(i, j)));
          }

      Matrix<Real_wp, Symmetric, RowSymPacked> mass_matrix = Mh;
      bool mass_diagonal = true;
      for (int i = 0; i < Fb.GetNbDof(); i++)
        for (int j = i+1; j < Fb.GetNbDof(); j++)
          if (abs(Mh(i, j))/max_mij > threshold)
            mass_diagonal = false;
      
      if (Fb.LumpedMassMatrix())
        {
          if (!mass_diagonal)
            {
              cout << "Mass matrix not lumped " << endl;
              CHECK_ABORT;
            }
        }
      else
        {
          if (mass_diagonal)
            {
              cout <<"Warning : mass matrix declared as not diagonal, but is diagonal " << endl;
              // CHECK_ABORT;
            }
        }
      
      InvMh = Mh;
      GetInverse(InvMh);
      Real_wp max_invMij = 0;
      for (int i = 0; i < Fb.GetNbDof(); i++)
        for (int j = i; j < Fb.GetNbDof(); j++)
          max_invMij = max(max_invMij, abs(InvMh(i, j)));
            
      VectReal_wp Ones(Fb.GetNbDof()), MhOnes(Fb.GetNbDof());
      Ones.Fill(0); MhOnes.Fill(0);
      bool mlt_mass_matrix_correct = true;
      bool solve_mass_matrix_correct = true;
      bool solve_cholesky_correct = true;
      for (int i = 0; i < Fb.GetNbDof(); i++)
        {
          Ones(i) = 1.0; MhOnes = Ones;
          Fb.MltMassMatrix(MhOnes);
          for (int j = 0; j < Fb.GetNbDof(); j++)
            if (abs(MhOnes(j) - Mh(j, i))/max_mij > 100*threshold)
              {
                cout << "MltMassMatrix not correct " << endl;
                mlt_mass_matrix_correct = false;
				DISP(i); DISP(j); DISP(MhOnes(j)); DISP(Mh(j,i));
                CHECK_ABORT;
              }
          
          MhOnes = Ones;
          Fb.SolveCholesky(SeldonNoTrans, MhOnes);
          Fb.SolveCholesky(SeldonTrans, MhOnes);
          for (int j = 0; j < Fb.GetNbDof(); j++)
            if (abs(InvMh(j, i) - MhOnes(j))/max_invMij > 100*threshold)
              {
                cout << "SolveCholesky not correct " << endl;
                solve_cholesky_correct = false;
                CHECK_ABORT;
              }
          
          MhOnes = Ones;
          Fb.SolveMassMatrix(MhOnes);
          for (int j = 0; j < Fb.GetNbDof(); j++)
            if (abs(InvMh(j, i) - MhOnes(j))/max_invMij > 100*threshold)
              {
                cout << "SolveMassMatrix not correct " << endl;
                solve_mass_matrix_correct = false;
                CHECK_ABORT;
              }
          
          Ones(i) = 0.0;
        }
      
      if (!mlt_mass_matrix_correct)
        cout << "Warning : MltMassMatrix not correct " << endl;

      if (!solve_mass_matrix_correct)
        cout << "Warning : SolveMassMatrix not correct " << endl;

      if (!solve_cholesky_correct)
        cout << "Warning : SolveCholesky not correct " << endl;
      
      VectReal_wp jacob(Fb.GetNbPointsQuadratureInside());
      jacob.FillRand(); Mlt(1e-9, jacob); Mh.Fill(0);
      for (int i = 0; i < Fb.GetNbDof(); i++)
        for (int j = i; j < Fb.GetNbDof(); j++)
          {
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              Mh(i, j) += Fb.WeightsND(k)*jacob(k)*DotProd(ValuePhi(i, k), ValuePhi(j, k));
            
          }

      bool sparse_mat = false;
      for (int i = 0; i < Fb.GetNbDof(); i++)
        for (int j = 0; j < Fb.GetNbDof(); j++)
          if (abs(Mh(i, j))/ max_mij <= threshold)
            sparse_mat = true;
      
      if (Fb.SparseMassMatrix())
        {
          if (!sparse_mat)
            {
              cout << "Mass matrix declared sparse but is dense " << endl;
              CHECK_ABORT;
            }
        }
      else
        {
          if (!Fb.LumpedMassMatrix())
            if (sparse_mat)
              cout << "Warning : mass matrix sparse but declared as dense " << endl;
        }
      
      if (Fb.OptimizedComputationMassMatrix())
        {
          for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
            jacob(k) *= Fb.WeightsND(k);
          
          Fb.IntegrateMassMatrix(InvMh, jacob);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            for (int j = 0; j < Fb.GetNbDof(); j++)
              if (abs(Mh(i, j) - InvMh(i, j))/ max_mij > threshold)
                {
                  cout << "IntegrateMassMatrix not correct " << endl;
                  CHECK_ABORT;
                }
        }
      
      Fb.DFjElemQuadrature(vert_elt, PointsElem, MatricesElem, mesh, 0);
      for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
        jacob(k) = Det(MatricesElem.GetPointQuadrature(k));
      
      Mh.Fill(0);
      for (int i = 0; i < Fb.GetNbDof(); i++)
        for (int j = i; j < Fb.GetNbDof(); j++)
          {
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              Mh(i, j) += Fb.WeightsND(k)*jacob(k)*DotProd(ValuePhi(i, k), ValuePhi(j, k));            
          }
      
      if (Fb.LinearSparseMassMatrix())
        {
          Fb.ComputeCoefJacobian(vert_elt, jacob);
          Fb.ComputeMassMatrix(InvMh, jacob);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            for (int j = 0; j < Fb.GetNbDof(); j++)
              if (abs(Mh(i, j) - InvMh(i, j))/ max_mij > threshold)
                {
                  cout << "ComputeMassMatrix not correct " << endl;
                  CHECK_ABORT;
                }
        }
      
      // checking ApplyCh/ComputeIntegralRef (same methods)
      Ones.Reallocate(2*Fb.GetNbPointsQuadratureInside()); Ones.Fill(0);
      for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
        for (int n = 0; n < 2; n++)
          {
            Ones(2*j+n) = 1.0;
            Fb.ComputeIntegralRef(Ones, MhOnes);
            for (int i = 0; i < Fb.GetNbDof(); i++)
              if (abs(MhOnes(i) - ValuePhi(i, j)(n))/max_cij > threshold)
                {
                  cout << "ComputeIntegralRef not correct " << endl;
                  DISP(MhOnes(i)); 
                  DISP(ValuePhi(i, j)(n));
                  CHECK_ABORT;
                }
            
            Fb.ApplyCh(Ones, MhOnes);
            for (int i = 0; i < Fb.GetNbDof(); i++)
              if (abs(MhOnes(i) - ValuePhi(i, j)(n))/max_cij > threshold)
                {
                  cout << "ApplyCh not correct " << endl;
                  CHECK_ABORT;
                }
            
            Ones(2*j+n) = 0.0;
          }
      
      // and ApplyChTranspose
      MhOnes.Fill(0);
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          MhOnes(j) = 1.0;
          Fb.ApplyChTranspose(MhOnes, Ones);
          for (int i = 0; i < Fb.GetNbPointsQuadratureInside(); i++)
            for (int n = 0; n < 2; n++)
              if (abs(Ones(2*i+n) - ValuePhi(j, i)(n))/max_cij > threshold)
                {
                  cout << "ApplyChTranspose not correct " << endl;
                  CHECK_ABORT;
                }
          
          MhOnes(j) = 0.0;
        }
      
      // checking ProjectQuadratureToDofRef
      VectReal_wp ProjOnes(Fb.GetNbDof());
      ProjOnes.Fill(0); Ones.Fill(0);
      bool proj_quadrature_to_dof_correct = true;
      for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
        for (int n = 0; n < 2; n++)
          {
            Ones(2*j+n) = 1.0;
            Fb.ProjectQuadratureToDofRef(Ones, MhOnes);
            Ones(2*j+n) *= Fb.WeightsND(j);
            Fb.ComputeIntegralRef(Ones, ProjOnes);
            Fb.SolveMassMatrix(ProjOnes);
            for (int i = 0; i < Fb.GetNbDof(); i++)
              if (abs(MhOnes(i) - ProjOnes(i)) > threshold)
                {
                  //cout << "ProjectQuadratureToDof is not correct " << endl;
                  proj_quadrature_to_dof_correct = false;
                  //CHECK_ABORT;
                }
            
            Ones(2*j+n) = 0.0;
          }

      if (!proj_quadrature_to_dof_correct)
        cout << "Warning : ProjectQuadratureToDofRef not correct " << endl;
      
      // checking ComputeProjectionDofRef
      VectReal_wp feval(2*Fb.GetNbPointsDofInside()); feval.Fill(0);
      for (int comp = 0; comp < 2; comp++)
        {
          for (int m = 0; m <= r-1; m++)
            for (int n = 0; n <= r-1-m; n++)
              {
		feval.Zero();
                for (int j = 0; j < Fb.GetNbPointsDofInside(); j++)
                  {
                    Real_wp x = Fb.PointsDofND(j)(0);
                    Real_wp y = Fb.PointsDofND(j)(1);
                    feval(2*j+comp) = pow(x, m)*pow(y, n);
                  }                  
		
		int comp2 = 1-comp;
		Fb.ComputeProjectionDofRef(feval, ProjOnes);
		Fb.ApplyChTranspose(ProjOnes, Ones);
		for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
		  {
		    Real_wp x = Fb.PointsND(j)(0);
		    Real_wp y = Fb.PointsND(j)(1);
		    Real_wp val_exact = pow(x, m)*pow(y, n);
		    if ((abs(val_exact - Ones(2*j+comp)) > threshold)
			&& (abs(Ones(2*j+comp2)) > threshold))
		      {
			cout << "ComputeProjectionDofRef not correct " << endl;
			DISP(j); DISP(comp); DISP(val_exact); DISP(Ones(2*j+comp));
			CHECK_ABORT;
		      }
		  }
	      }
        }
      
      // checking ComputeNodalValuesRef and ComputeNodalGradientRef
      Matrix<R2> ValuePhiNodal(Fb.GetNbDof(), Fb.GetNbPointsNodalElt());
      Matrix<R2> GradPhiNodal(Fb.GetNbPointsNodalElt(), Fb.GetNbPointsNodalElt());
      for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
        {
          Fb.ComputeValuesPhiRef(Fb.PointsNodalND(i), phi);
          Fb.ComputeGradientPhiNodalRef(Fb.PointsNodalND(i), grad_phi);          
          for (int j = 0; j < Fb.GetNbDof(); j++)
            ValuePhiNodal(j, i) = phi(j);
          
          for (int j = 0; j < Fb.GetNbPointsNodalElt(); j++)
            GradPhiNodal(j, i) = grad_phi(j);
        }
      
      Ones.Reallocate(Fb.GetNbDof()); Ones.Fill(0);
      MhOnes.Reallocate(2*Fb.GetNbPointsNodalElt()); MhOnes.Fill(0);
      VectReal_wp Ex(Fb.GetNbPointsNodalElt());
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          Ones(j) = 1.0;
          Fb.ComputeNodalValuesRef(Ones, MhOnes);
          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            for (int n = 0; n < 2; n++)
              if (abs(ValuePhiNodal(j, i)(n) - MhOnes(2*i+n)) > threshold)
                {
                  cout << "ComputeNodalValuesRef not correct " << endl;
                  DISP(MhOnes(2*i+n)); DISP(ValuePhiNodal(j, i)(n));
                  CHECK_ABORT;
                }
          
          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            Ex(i) = MhOnes(2*i);
          
          Fb.ComputeNodalGradientRef(Ex, grad_phi);          
          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            {
              R2 grad;
              for (int k = 0; k < Fb.GetNbPointsNodalElt(); k++)
                Add(MhOnes(2*k), GradPhiNodal(k, i), grad);
              
              if (grad != grad_phi(i))
                {
                  cout << "ComputeNodalGradientRef not correct " << endl;
                  CHECK_ABORT;
                }
            }

          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            Ex(i) = MhOnes(2*i+1);
          
          Fb.ComputeNodalGradientRef(Ex, grad_phi);          
          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            {
              R2 grad;
              for (int k = 0; k < Fb.GetNbPointsNodalElt(); k++)
                Add(MhOnes(2*k+1), GradPhiNodal(k, i), grad);
              
              if (grad != grad_phi(i))
                {
                  cout << "ComputeNodalGradientRef not correct " << endl;
                  CHECK_ABORT;
                }
            }
          
          Ones(j) = 0.0;
        }
      
      
      // checking curl of basis functions
      VectR2 phi_p, phi_m;
      for (int i = 0; i < Fb.GetNbPointsDof(); i++)
        {
          Fb.ComputeCurlPhiRef(Fb.PointsDofND(i), curl_phi);
          Real_wp val_max = 0;
          for (int j = 0; j < Fb.GetNbDof(); j++)
            val_max = max(val_max, abs(curl_phi(j)));

          VectReal_wp curlphi_num(Fb.GetNbDof()); 
          FillZero(curlphi_num);
          
          ptA = Fb.PointsDofND(i);
          ptA(0) += h;
          Fb.ComputeValuesPhiRef(ptA, phi_p);
          
          ptA(0) -= 2.0*h;
          Fb.ComputeValuesPhiRef(ptA, phi_m);
          
          for (int j = 0; j < Fb.GetNbDof(); j++)
            curlphi_num(j) += (phi_p(j)(1) - phi_m(j)(1))/(2.0*h);
                      
          ptA = Fb.PointsDofND(i);
          ptA(1) += h;
          Fb.ComputeValuesPhiRef(ptA, phi_p);
          
          ptA(1) -= 2.0*h;
          Fb.ComputeValuesPhiRef(ptA, phi_m);
          
          for (int j = 0; j < Fb.GetNbDof(); j++)
            curlphi_num(j) -= (phi_p(j)(0) - phi_m(j)(0))/(2.0*h);
                     
          for (int j = 0; j < Fb.GetNbDof(); j++)
            {
              if (abs(curlphi_num(j) - curl_phi(j)) > val_max*sqrt(threshold))
                {
                  DISP(Fb.PointsDofND(i)); DISP(j); DISP(curlphi_num(j)); DISP(curl_phi(j));
                  cout << "Curl of basis functions not correct " <<endl;
                  CHECK_ABORT;
                }
            }
        }

      // checking GetValuePhiOnQuadraturePoint (and Gradient)
      Matrix<Real_wp> CurlPhi(Fb.GetNbDof(), Fb.GetNbPointsQuadrature());
      Real_wp max_rij = 0;
      for (int k = 0; k < Fb.GetNbPointsQuadrature(); k++)
        {
          Fb.ComputeCurlPhiRef(Fb.PointsND(k), curl_phi);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              CurlPhi(i, k) = curl_phi(i);
              max_rij = max(max_rij, abs(CurlPhi(i, k)));
            }
        }
      
      for (int k = 0; k < Fb.GetNbPointsQuadrature(); k++)
        {
          Fb.GetValuePhiOnQuadraturePoint(k, phi);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (phi(i) != ValuePhi(i, k))
              {
                cout << "GetValuePhiOnQuadraturePoint not correct " << endl;
                DISP(phi(i)); DISP(ValuePhi(i, k));
                CHECK_ABORT;
              }
          
          Fb.GetCurlPhiOnQuadraturePoint(k, curl_phi);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(curl_phi(i) - CurlPhi(i, k)) > threshold)
              {
                cout << "GetCurlPhiOnQuadraturePoint not correct " << endl;
                DISP(curl_phi(i)); DISP(CurlPhi(i, k));
                CHECK_ABORT;
              }
        }
      
      // checking ComputeIntegralGradientRef / ApplyRh
      Ones.Reallocate(Fb.GetNbPointsQuadratureInside());
      MhOnes.Reallocate(Fb.GetNbDof()); MhOnes.Fill(0);
      Ones.Fill(0);
      for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
        {
          Ones(j) = 1.0;
          Fb.ComputeIntegralCurlRef(Ones, MhOnes);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(MhOnes(i) - CurlPhi(i, j))/max_rij > threshold)
              {
                cout << "ComputeIntegralCurlRef not correct " << endl;
                CHECK_ABORT;
              }

          Fb.ApplyRh(Ones, MhOnes);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(MhOnes(i) - CurlPhi(i, j))/max_rij > threshold)
              {
                cout << "ApplyRh not correct " << endl;
                CHECK_ABORT;
              }
          
          Ones(j) = 0.0;
        }
      
      // checking ApplyRhTranspose
      MhOnes.Fill(0);
      for (int i = 0; i < Fb.GetNbDof(); i++)
        {
          MhOnes(i) = 1.0;
          
          Fb.ApplyRhTranspose(MhOnes, Ones);
          for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
            if (abs(Ones(j) - CurlPhi(i, j))/max_rij > threshold)
              {
                cout << "ApplyRhTranspose not correct " << endl;
                CHECK_ABORT;
              }
          
          MhOnes(i) = 0.0;
        }
      
      // checking ComputeIntegralSurfaceRef/ApplySh
      ComputeGaussLegendre(points1d, weights1d, rquad);

      Matrix<bool> IsDofOnBoundary(Fb.GetNbBoundaries(), Fb.GetNbDof());
      IsDofOnBoundary.Fill(false);
      
      for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
        {
          int nb_pts_quad = Fb.GetNbQuadBoundary(num_loc);
          Ones.Reallocate(2*nb_pts_quad); Ones.Fill(0);
          for (int k = 0; k < nb_pts_quad; k++)
            for (int n = 0; n < 2; n++)
              {
                Ones(2*k+n) = 1.0; MhOnes.Fill(1.0);
                int k2 = Fb.GetQuadNumber(num_loc, k);
                
                Fb.ComputeIntegralSurfaceRef(Ones, MhOnes, num_loc);
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  if (abs(MhOnes(i) - ValuePhi(i, k2)(n))/max_cij > threshold)
                    {
                      DISP(MhOnes(i)); DISP(ValuePhi(i, k2)(n)); DISP(i); DISP(num_loc); DISP(k2); DISP(n);
                      cout << "ComputeIntegralSurfaceRef not correct " << endl;
                      CHECK_ABORT;
                    }
                
                Ones(2*k+n) = 0.0;
              }
          
          Ones.Fill(0);
          for (int k = 0; k < nb_pts_quad; k++)
            {
              Ones(k) = 1.0;
              int k2 = Fb.GetQuadNumber(num_loc, k);
              Fb.ComputeIntegralSurfaceCurlRef(Ones, MhOnes, num_loc);
              for (int i = 0; i < Fb.GetNbDof(); i++)
                if (abs(MhOnes(i) - CurlPhi(i, k2))/max_rij > threshold)
                  {
                    cout << "ComputeIntegralSurfaceCurlRef not correct " << endl;
                    DISP(MhOnes(i)); DISP(CurlPhi(i, k2)); DISP(num_loc);
                    CHECK_ABORT;
                  }
              
              Ones(k) = 0.0;
            }
          
          nb_pts_quad = points1d.GetM();
          
          Ones.Reallocate(2*nb_pts_quad); Ones.Fill(0);
          for (int k = 0; k < nb_pts_quad; k++)
            for (int n = 0; n < 2; n++)
              {
                Ones(2*k+n) = 1.0; MhOnes.Fill(1.0);
                
                Fb.GetLocalCoordOnBoundary(num_loc, points1d(k), pt_loc);
                Fb.ComputeValuesPhiRef(pt_loc, phi);
              
                Fb.ComputeGaussIntegralSurfaceRef(Ones, MhOnes, num_loc);
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  if (abs(MhOnes(i) - phi(i)(n))/max_cij > threshold)
                    {
                      cout << "ComputeGaussIntegralSurfaceRef not correct " << endl;
                      CHECK_ABORT;
                    }
                
                Ones(2*k+n) = 0.0;
		
		// detecting surface dofs
		for (int i = 0; i < Fb.GetNbDof(); i++)
		  {
		    Real_wp phi_times_n = phi(i)(0)*Fb.NormaleLoc(num_loc)(1)
		      - phi(i)(1)*Fb.NormaleLoc(num_loc)(0);
		    
		    if (abs(phi_times_n) > threshold)
		      IsDofOnBoundary(num_loc, i) = true;
		  }
              }
	  
        }
      
      {
	int nb_dof_boundaries_ref = 0;
	for (int i = 0; i < Fb.GetNbDof(); i++)
	  {
	    bool dof_present = false;
	    for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	      if (IsDofOnBoundary(num_loc, i))
		dof_present = true;
	    
	    if (dof_present)
	      nb_dof_boundaries_ref++;
	  }
	
	DISP(nb_dof_boundaries_ref);
	if (nb_dof_boundaries_ref != Fb.GetNbDofBoundaries())
	  {
	    cout << "GetNbDofBoundaries not correct " << endl;
	    DISP(nb_dof_boundaries_ref); DISP(Fb.GetNbDofBoundaries());
	    CHECK_ABORT;
	  }
	
	for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	  {
	    nb_dof_boundaries_ref = 0;
	    for (int i = 0; i < Fb.GetNbDof(); i++)
	      {
		if (IsDofOnBoundary(num_loc, i))
		  nb_dof_boundaries_ref++;
		
		if (IsDofOnBoundary(num_loc, i) != Fb.IsTangentialDof(i, num_loc))
		  {
		    cout << "IsTangentialDof incorrect" << endl;
		    DISP(i); DISP(num_loc);
		    DISP(IsDofOnBoundary(num_loc, i));
		    DISP(Fb.IsTangentialDof(i, num_loc));
		    CHECK_ABORT;
		  }
	      }
	    
	    if (Fb.GetNbDofBoundary(num_loc) != nb_dof_boundaries_ref)
	      {
		cout << "GetNbDofBoundary incorrect" << endl;
		CHECK_ABORT;
	      }
	  }
      }
      
      // checking ComputeProjectionSurfaceDofRef
      if (!Fb.DiscontinuousElement())
	for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
          {
	    if ((Fb.GetNbPointsDofInside() == 0) || (Fb.GetNbPointsDofInside() < Fb.GetNbDof()/2))
	      {
		cout << "GetNbPointsDofInside() incorrect" << endl;
		DISP(Fb.GetNbPointsDofInside());
		abort();
	      }

	    if ((Fb.GetNbPointsDofSurface(num_loc) == 0)
		|| (Fb.GetNbPointsDofSurface(num_loc) < Fb.GetNbDofBoundary(num_loc)))
	      {
		cout << "GetNbPointsDofSurface incorrect" << endl;
		DISP(Fb.GetNbPointsDofSurface(num_loc));
		abort();
	      }

	    const VectReal_wp& pts_face = Fb.PointsDofBoundary(num_loc);
	    if (pts_face.GetM() != Fb.GetNbPointsDofSurface(num_loc))
	      {
		cout << "PointsDofBoundary incorrect" << endl;
		abort();
	      }

	    for (int i = 0; i < Fb.GetNbPointsDofSurface(num_loc); i++)
	      {
		R2 pt_loc;
		Fb.GetLocalCoordOnBoundary(num_loc, pts_face(i), pt_loc);
		if (pt_loc != Fb.PointsDofND(Fb.GetPointDofNumber(num_loc, i)))
		  {
		    DISP(num_loc); DISP(i);
		    DISP(pt_loc); DISP(Fb.PointsDofND(Fb.GetPointDofNumber(num_loc, i)));
		    cout << "PointsDofBoundary incorrect" << endl;
		    abort();
		  }
	      }
	    
	    for (int p = 0; p < 2; p++)
	      for (int m = 0; m <= r-1; m++)
		for (int n = 0; n <= r-1-m; n++)
		  {
		    Vector<Real_wp> feval_surf(2*Fb.GetNbPointsDofSurface(num_loc));
		    feval_surf.Zero();
		    for (int i = 0; i < Fb.GetNbPointsDofSurface(num_loc); i++)
		      feval_surf(2*i+p) = FuncPowS(Fb.PointsDofND(Fb.GetPointDofNumber(num_loc, i)), m, n);
		    
		    VectReal_wp contrib_face(Fb.GetNbDofBoundary(num_loc));
		    Fb.ComputeProjectionSurfaceDofRef(feval_surf, contrib_face, num_loc);
		    
		    VectReal_wp feval_vol(2*Fb.GetNbPointsDofInside());
		    feval_vol.Zero();
		    for (int i = 0; i < Fb.GetNbPointsDofInside(); i++)
		      feval_vol(2*i+p) = FuncPowS(Fb.PointsDofND(i), m, n);
		    
		    Fb.ComputeProjectionDofRef(feval_vol, ProjOnes);

		    for (int i = 0; i < Fb.GetNbDofBoundary(num_loc); i++)
		      if (abs(contrib_face(i) - ProjOnes(Fb.GetLocalNumber(num_loc, i))) > 10*r*threshold)
			{
			  cout << "ComputeProjectionSurfaceDofRef not correct " << endl;
			  DISP(num_loc);
			  DISP(contrib_face(i));
			  DISP(ProjOnes(Fb.GetLocalNumber(num_loc, i)));
			  CHECK_ABORT;
			}
		  }
          }
      
      // checking FluxWeight
      for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
        {
          VectReal_wp poids = Fb.GetFluxWeight(num_loc);
          for (int k = 0; k < Fb.GetNbQuadBoundary(num_loc); k++)
            if (abs(poids(k) - 0.5*Fb.WeightsQuadratureBoundary(k, num_loc)) > threshold)
              {
                cout << "PoidsFlux not correct" << endl;
                CHECK_ABORT;
              }
        }
      
      // checking addition of elementary matrices
      int nb_dof_loc = Fb.GetNbDof(); int off_row = 12, off_col = 30;
      Matrix2_2 Bmass; Bmass.FillRand(); Mlt(1e-9, Bmass);
      Matrix<Real_wp> AhOrig(off_row+nb_dof_loc, off_col+nb_dof_loc), AhNew(off_row+nb_dof_loc, off_col+nb_dof_loc);
      AhOrig.FillRand(); Mlt(1e-9*max_mij, AhOrig); AhNew = AhOrig;
      Fb.AddConstantMassMatrix(off_row, off_col, Bmass, AhNew);
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val_exact = 0;
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              val_exact += DotProd(ValuePhi(i, k), dot(Bmass, ValuePhi(j, k)))*Fb.WeightsND(k);
            
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_mij > threshold)
              {
                DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j)); DISP(val_exact);
                cout << "AddConstantMassMatrix not correct " << endl;
                CHECK_ABORT;
              }
          }
      
      AhNew = AhOrig;
      Real_wp Cgrad  = 0.87;
      Fb.AddConstantStiffnessMatrix(off_row, off_col, Cgrad, AhNew);
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val_exact = 0;
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              val_exact += Cgrad*CurlPhi(i, k)*CurlPhi(j, k)*Fb.WeightsND(k);
            
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_rij > threshold)
              {
                cout << "AddConstantStiffnessMatrix not correct " << endl;
                DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j)); DISP(val_exact);
                CHECK_ABORT;
              }
          }


      AhNew = AhOrig;
      Vector<TinyMatrix<Real_wp, General, 2, 2> > Bvar(Fb.GetNbPointsQuadratureInside());
      for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
	{
	  Bvar(k).FillRand();
	  Bvar(k) *= 1e-9;	 
	  Bvar(k) *= Fb.WeightsND(k);
	}
      
      Fb.AddVariableMassMatrix(off_row, off_col, Bvar, AhNew);
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val_exact = 0;
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              val_exact += DotProd(ValuePhi(i, k), dot(Bvar(k), ValuePhi(j, k)));
            
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_mij > threshold)
              {
                DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j)); DISP(val_exact);
                cout << "AddVariableMassMatrix not correct " << endl;
                CHECK_ABORT;
              }
          }
      
      AhNew = AhOrig;
      VectReal_wp Cvar(Fb.GetNbPointsQuadratureInside());
      Cvar.FillRand(); Mlt(1e-9, Cvar);
      for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
	Cvar(k) *= Fb.WeightsND(k);
      
      Fb.AddVariableStiffnessMatrix(off_row, off_col, Cvar, AhNew);
      if (check_stiff)
	for (int i = 0; i < nb_dof_loc; i++)
	  for (int j = 0; j < nb_dof_loc; j++)
	    {
	      Real_wp val_exact = 0;
	      for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
		val_exact += Cvar(k)*CurlPhi(i, k)*CurlPhi(j, k);
	      
	      if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_rij > threshold)
		{
		  cout << "AddVariableStiffnessMatrix not correct " << endl;
		  DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j)); DISP(val_exact);
		  CHECK_ABORT;
		}
	    }
      
      // checking ValueBoundaryRef, and GradientBoundaryRef
      Ones.Reallocate(Fb.GetNbDof()); Ones.Fill(0);
      MhOnes.Reallocate(2*Fb.GetNbPointsNodalElt()); MhOnes.Fill(0);
      VectR2 dEx, dEy; VectReal_wp feval_surf, Ey;
      Ex.Reallocate(Fb.GetNbPointsNodalElt());
      Ey.Reallocate(Fb.GetNbPointsNodalElt());
            
      dEx.Reallocate(Fb.GetNbPointsNodalElt());
      dEy.Reallocate(Fb.GetNbPointsNodalElt());
      if (r == rgeom)
        for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
          {
            feval_surf.Reallocate(2*Fb.GetNbNodalBoundary(num_loc));
            feval.Reallocate(2*Fb.GetNbNodalBoundary(num_loc));
            for (int i = 0; i < Fb.GetNbDof(); i++)
              {
                Ones(i) = 1.0;
                Fb.ComputeNodalValuesRef(Ones, MhOnes);
                for (int j = 0; j < Fb.GetNbPointsNodalElt(); j++)
                  {
                    Ex(j) = MhOnes(2*j);
                    Ey(j) = MhOnes(2*j+1);
                  }
                
                Fb.ComputeNodalGradientRef(Ex, dEx);
                Fb.ComputeNodalGradientRef(Ey, dEy);
                                
                Fb.ComputeValueBoundaryRef(Ones, feval_surf, num_loc);
                Fb.ComputeCurlBoundaryRef(Ones, feval, num_loc);
                
                for (int k = 0; k < Fb.GetNbNodalBoundary(num_loc); k++)
                  {
                    int j = Fb.GetNodalNumber(num_loc, k);
                    if ((abs(feval_surf(2*k) - MhOnes(2*j)) > threshold) 
                        || (abs(feval_surf(2*k+1) - MhOnes(2*j+1)) > threshold))
                      {
                        cout << "ComputeValueBoundaryRef not correct " << endl;
                        CHECK_ABORT;
                      }
                    
                  }
                
                Ones(i) = 0.0;
              }
          }
    }
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  // 2-D mesh with only one triangle
  Mesh<Dimension2> mesh2d;
  mesh2d.SetGeometryOrder(4);
  
  mesh2d.ReallocateVertices(3);
  mesh2d.Vertex(0).Init(-0.2, 0.2);
  mesh2d.Vertex(1).Init(1, 0);
  mesh2d.Vertex(2).Init(0.5, 0.5*sqrt(3.0));
  
  mesh2d.ReallocateBoundariesRef(3);
  mesh2d.BoundaryRef(0).Init(0, 1, 1);
  mesh2d.BoundaryRef(1).Init(1, 2, 2);
  mesh2d.BoundaryRef(2).Init(0, 2, 1);
  
  mesh2d.ReallocateElements(1);
  mesh2d.Element(0).InitTriangular(0, 1, 2, 1);
  
  mesh2d.SetCurveType(2, mesh2d.CURVE_CIRCLE);
  VectReal_wp param(3);
  param(0) = 0; param(1) = 0; param(2) = 1.0;
  mesh2d.SetCurveParameter(2, param);
  
  mesh2d.ReorientElements();
  mesh2d.FindConnectivity();
  mesh2d.ProjectPointsOnCurves();

  // checking triangular finite elements
  Vector<R2> vertex_tri(3);
  vertex_tri(0).Init(0, 0); vertex_tri(1).Init(1, 0); vertex_tri(2).Init(0, 1);

  {
    cout << "Testing TriangleHcurlFirstFamily " << endl;
    TriangleHcurlFirstFamily tri; string name_elt("TRIANGLE_FIRST_FAMILY");
    CheckHcurlFiniteElement(tri, name_elt, mesh2d, vertex_tri);
  }

  {
    cout << "Testing TriangleHcurlOptimalFirstFamily " << endl;
    TriangleHcurlOptimalFirstFamily tri; string name_elt("TRIANGLE_OPTIMAL_FIRST_FAMILY");
    CheckHcurlFiniteElement(tri, name_elt, mesh2d, vertex_tri);
  }

  {
    cout << "Testing TriangleHcurlHierarchic " << endl;
    TriangleHcurlHierarchic tri; string name_elt("TRIANGLE_HP_SECOND_FAMILY");
    CheckHcurlFiniteElement(tri, name_elt, mesh2d, vertex_tri, false);
  }
    
  {
    cout << "Testing TriangleHcurlSecondFamily " << endl;
    TriangleHcurlSecondFamily tri;  string name_elt("TRIANGLE_SECOND_FAMILY");
    CheckHcurlFiniteElement(tri, name_elt, mesh2d, vertex_tri);
  }

  {
    cout << "Testing TriangleHcurlOptimalHpFirstFamily " << endl;
    TriangleHcurlOptimalHpFirstFamily tri;  string name_elt("TRIANGLE_OPTIMAL_HP_FIRST_FAMILY");
    CheckHcurlFiniteElement(tri, name_elt, mesh2d, vertex_tri);
  }


  // checking quadrilateral finite elements
  mesh2d.ResizeVertices(4);
  mesh2d.Vertex(3).Init(-0.3, 1.2);
  
  mesh2d.ResizeBoundariesRef(4);
  mesh2d.BoundaryRef(2).Init(2, 3, 1);
  mesh2d.BoundaryRef(3).Init(0, 3, 1);
    
  mesh2d.Element(0).InitQuadrangular(0, 1, 2, 3, 1);
  mesh2d.ReorientElements();
  mesh2d.FindConnectivity();
  mesh2d.ProjectPointsOnCurves();
  
  Vector<R2> vertex_quad(4);
  vertex_quad(0).Init(0, 0); vertex_quad(1).Init(1, 0);
  vertex_quad(2).Init(1, 1); vertex_quad(3).Init(0, 1);
  
  {
    cout << endl << endl << "Testing QuadrangleHcurlFirstFamily " << endl;
    QuadrangleHcurlFirstFamily quad;  string name_elt("QUADRANGLE_FIRST_FAMILY");
    CheckHcurlFiniteElement(quad, name_elt, mesh2d, vertex_quad, true, false);
  }

  {
    cout << endl << endl << "Testing QuadrangleHcurlGaussFirstFamily " << endl;
    QuadrangleHcurlGaussFirstFamily quad; string name_elt("QUADRANGLE_GAUSS_FIRST_FAMILY");
    CheckHcurlFiniteElement(quad, name_elt, mesh2d, vertex_quad);
  }

  {
    cout << endl << endl << "Testing QuadrangleHcurlOptimalFirstFamily " << endl;
    QuadrangleHcurlOptimalFirstFamily quad;  string name_elt("TRIANGLE_OPTIMAL_FIRST_FAMILY");
    CheckHcurlFiniteElement(quad, name_elt, mesh2d, vertex_quad);
  }
  
  {
    cout << endl << endl << "Testing QuadrangleHcurlLobatto " << endl;
    QuadrangleHcurlLobatto quad; string name_elt("QUADRANGLE_HCURL_LOBATTO");
    CheckHcurlFiniteElement(quad, name_elt, mesh2d, vertex_quad);
  }

  {
    cout << endl << endl << "Testing QuadrangleHcurlHpFirstFamily " << endl;
    QuadrangleHcurlHpFirstFamily quad; string name_elt("TRIANGLE_HP_FIRST_FAMILY");
    CheckHcurlFiniteElement(quad, name_elt, mesh2d, vertex_quad);
  }
  
  {
    cout << endl << endl << "Testing QuadrangleHcurlOptimalHpFirstFamily " << endl;
    QuadrangleHcurlOptimalHpFirstFamily quad; string name_elt("TRIANGLE_OPTIMAL_HP_FIRST_FAMILY");
    CheckHcurlFiniteElement(quad, name_elt, mesh2d, vertex_quad);
  }

  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
