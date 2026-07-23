#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_NODAL_H1

#define MONTJOIE_WITH_HP_HCURL

#include "Harmonic/MontjoieHarmonic.hxx"

//#include "Elliptic/Maxwell/MontjoieMaxwell3D.hxx"
#include "FiniteElement/Pyramid/PyramidHcurlOther.hxx"
#include "FiniteElement/Pyramid/PyramidHcurlOther.cxx"

#ifdef MONTJOIE_WITH_FAST_LIBRARY
// includes for fast compilation
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#endif

#define CHECK_ABORT abort()
//#define CHECK_ABORT cout<<"ERROR"<<endl

using namespace Montjoie;

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
bool EqualVector(const Vector<T>& x, const Vector<T>& y)
{
  if (x.GetM() != y.GetM())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    if (abs(x(i)-y(i)) > 1e-12)
      return false;
  
  return true;
}

template<class T, int p>
bool EqualVector(const Vector<TinyVector<T, p> >& x, const Vector<TinyVector<T, p> >& y)
{
  if (x.GetM() != y.GetM())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    if (x(i) != y(i))
      return false;
  
  return true;
}

Real_wp FuncPowS(const R3& x, int m, int n)
{
  return pow(x(0)+0.5*x(1)-0.3*x(2), m)*pow(0.2+0.9*x(0)+x(1)+0.8*x(2), n);
}

namespace Montjoie
{
  template<class Dim>
  class MyOwnEquation : public GenericEquation<Real_wp>
  {
  public :
    typedef Dim Dimension;
    enum{nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_vec = 0,
	 nb_components_en = 1, nb_components_hn = 1, nb_unknowns_hdg = 0, type_element = 2};
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
    
    virtual void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Real_wp>&, CondensationBlockSolver_Base<Real_wp>&,
					 const GlobalGenericMatrix<Real_wp>&) {}
    
    virtual void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Complex_wp>&, CondensationBlockSolver_Base<Complex_wp>&,
					 const GlobalGenericMatrix<Complex_wp>&) {}
        
  };
  
} // end namespace Montjoie

void CreateRegularMesh(Mesh<Dimension3>& mesh, const VectR3& s, int r, int rgeom)
{
  int type_mesh = Mesh<Dimension3>::HEXAHEDRAL_MESH;
  switch (s.GetM())
    {
    case 4 : type_mesh = Mesh<Dimension3>::TETRAHEDRAL_MESH; break;
      //case 4 : type_mesh = Mesh<Dimension3>::HYBRID_MESH; break;
    case 5 : type_mesh = Mesh<Dimension3>::PYRAMID_MESH; break;
    case 6 : type_mesh = Mesh<Dimension3>::WEDGE_MESH; break;
    }
  
  mesh.SetGeometryOrder(rgeom);
  TinyVector<int, 6> ref_bound; ref_bound.Fill(1);
  mesh.CreateRegularMesh(R3(-2,-2,-2), R3(2,2,2), TinyVector<int, 3>(3,3,3), 1,
			 ref_bound, type_mesh);
  
  mesh.ClearConnectivity(); //return;
  // random permutation of local vertex numbers
  //srand(time(NULL));
  srand(0);
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      int nb_vert = mesh.Element(i).GetNbVertices();
      int ref = mesh.Element(i).GetReference();
      IVect num(nb_vert), numb;
      for (int j = 0; j < nb_vert; j++)
	num(j) = mesh.Element(i).numVertex(j);
      
      numb = num;            
      
      switch (nb_vert)
	{
	case 4 :
	  {
	    int perm = rand()%12;
	    switch (perm)
	      {
	      case 0 :
		num = numb; break;
	      case 1 :
		num(0) = numb(1); num(1) = numb(2); num(2) = numb(0); num(3) = numb(3); break;
	      case 2 :
		num(0) = numb(2); num(1) = numb(0); num(2) = numb(1); num(3) = numb(3); break;
	      case 3 :
		num(0) = numb(1); num(1) = numb(3); num(2) = numb(2); num(3) = numb(0); break;
	      case 4 :
		num(0) = numb(2); num(1) = numb(1); num(2) = numb(3); num(3) = numb(0); break;
	      case 5 :
		num(0) = numb(3); num(1) = numb(2); num(2) = numb(1); num(3) = numb(0); break;
	      case 6 :
		num(0) = numb(0); num(1) = numb(2); num(2) = numb(3); num(3) = numb(1); break;
	      case 7 :
		num(0) = numb(3); num(1) = numb(0); num(2) = numb(2); num(3) = numb(1); break;
	      case 8 :
		num(0) = numb(2); num(1) = numb(3); num(2) = numb(0); num(3) = numb(1); break;
	      case 9 :
		num(0) = numb(0); num(1) = numb(3); num(2) = numb(1); num(3) = numb(2); break;
	      case 10 :
		num(0) = numb(1); num(1) = numb(0); num(2) = numb(3); num(3) = numb(2); break;
	      case 11 :
		num(0) = numb(3); num(1) = numb(1); num(2) = numb(0); num(3) = numb(2); break;
	      }
	  }
	  break;
	case 5 :
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
	case 6 :
	  {
	    int perm = rand()%6;
	    switch (perm)
	      {
	      case 0 :
		num = numb; break;
	      case 1 :
		num(0) = numb(1); num(1) = numb(2); num(2) = numb(0);
		num(3) = numb(4); num(4) = numb(5); num(5) = numb(3); break;
	      case 2 :
		num(0) = numb(2); num(1) = numb(0); num(2) = numb(1);
		num(3) = numb(5); num(4) = numb(3); num(5) = numb(4); break;
	      case 3 :
		num(0) = numb(3); num(1) = numb(5); num(2) = numb(4);
		num(3) = numb(0); num(4) = numb(2); num(5) = numb(1); break;
	      case 4 :
		num(0) = numb(4); num(1) = numb(3); num(2) = numb(5);
		num(3) = numb(1); num(4) = numb(0); num(5) = numb(2); break;
	      case 5 :
		num(0) = numb(5); num(1) = numb(4); num(2) = numb(3);
		num(3) = numb(2); num(4) = numb(1); num(5) = numb(0); break;
	      }
	  }
	  break;
	case 8 :
	  {
	    int perm = rand()%24;
	    switch (perm)
	      {
	      case 0 :
		num = numb; break;
	      case 1 :
		num(0) = numb(1); num(1) = numb(2); num(2) = numb(3); num(3) = numb(0);
		num(4) = numb(5); num(5) = numb(6); num(6) = numb(7); num(7) = numb(4); break;
	      case 2 :
		num(0) = numb(2); num(1) = numb(3); num(2) = numb(0); num(3) = numb(1);
		num(4) = numb(6); num(5) = numb(7); num(6) = numb(4); num(7) = numb(5); break;
	      case 3 :
		num(0) = numb(3); num(1) = numb(0); num(2) = numb(1); num(3) = numb(2);
		num(4) = numb(7); num(5) = numb(4); num(6) = numb(5); num(7) = numb(6); break;
	      case 4 :
		num(0) = numb(5); num(1) = numb(6); num(2) = numb(7); num(3) = numb(4);
		num(4) = numb(1); num(5) = numb(2); num(6) = numb(3); num(7) = numb(0); break;
	      case 5 :
		num(0) = numb(6); num(1) = numb(7); num(2) = numb(4); num(3) = numb(5);
		num(4) = numb(2); num(5) = numb(3); num(6) = numb(0); num(7) = numb(1); break;
	      case 6 :
		num(0) = numb(7); num(1) = numb(4); num(2) = numb(5); num(3) = numb(6);
		num(4) = numb(3); num(5) = numb(0); num(6) = numb(1); num(7) = numb(2); break;
	      case 7 :
		num(0) = numb(4); num(1) = numb(5); num(2) = numb(6); num(3) = numb(7);
		num(4) = numb(0); num(5) = numb(1); num(6) = numb(2); num(7) = numb(3); break;
	      case 8 :
		num(0) = numb(1); num(1) = numb(2); num(2) = numb(6); num(3) = numb(5);
		num(4) = numb(0); num(5) = numb(3); num(6) = numb(7); num(7) = numb(4); break;
	      case 9 :
		num(0) = numb(2); num(1) = numb(6); num(2) = numb(5); num(3) = numb(1);
		num(4) = numb(3); num(5) = numb(7); num(6) = numb(4); num(7) = numb(0); break;
	      case 10 :
		num(0) = numb(6); num(1) = numb(5); num(2) = numb(1); num(3) = numb(2);
		num(4) = numb(7); num(5) = numb(4); num(6) = numb(0); num(7) = numb(3); break;
	      case 11 :
		num(0) = numb(5); num(1) = numb(1); num(2) = numb(2); num(3) = numb(6);
		num(4) = numb(4); num(5) = numb(0); num(6) = numb(3); num(7) = numb(7); break;
	      case 12 :
		num(0) = numb(0); num(1) = numb(3); num(2) = numb(7); num(3) = numb(4);
		num(4) = numb(1); num(5) = numb(2); num(6) = numb(6); num(7) = numb(5); break;
	      case 13 :
		num(0) = numb(3); num(1) = numb(7); num(2) = numb(4); num(3) = numb(0);
		num(4) = numb(2); num(5) = numb(6); num(6) = numb(5); num(7) = numb(1); break;
	      case 14 :
		num(0) = numb(7); num(1) = numb(4); num(2) = numb(0); num(3) = numb(3);
		num(4) = numb(6); num(5) = numb(5); num(6) = numb(1); num(7) = numb(2); break;
	      case 15 :
		num(0) = numb(4); num(1) = numb(0); num(2) = numb(3); num(3) = numb(7);
		num(4) = numb(5); num(5) = numb(1); num(6) = numb(2); num(7) = numb(6); break;
	      case 16 :
		num(0) = numb(3); num(1) = numb(2); num(2) = numb(6); num(3) = numb(7);
		num(4) = numb(0); num(5) = numb(1); num(6) = numb(5); num(7) = numb(4); break;
	      case 17 :
		num(0) = numb(2); num(1) = numb(6); num(2) = numb(7); num(3) = numb(3);
		num(4) = numb(1); num(5) = numb(5); num(6) = numb(4); num(7) = numb(0); break;
	      case 18 :
		num(0) = numb(6); num(1) = numb(7); num(2) = numb(3); num(3) = numb(2);
		num(4) = numb(5); num(5) = numb(4); num(6) = numb(0); num(7) = numb(1); break;
	      case 19 :
		num(0) = numb(7); num(1) = numb(3); num(2) = numb(2); num(3) = numb(6);
		num(4) = numb(4); num(5) = numb(0); num(6) = numb(1); num(7) = numb(5); break;
	      case 20 :
		num(0) = numb(0); num(1) = numb(1); num(2) = numb(5); num(3) = numb(4);
		num(4) = numb(3); num(5) = numb(2); num(6) = numb(6); num(7) = numb(7); break;
	      case 21 :
		num(0) = numb(1); num(1) = numb(5); num(2) = numb(4); num(3) = numb(0);
		num(4) = numb(2); num(5) = numb(6); num(6) = numb(7); num(7) = numb(3); break;
	      case 22 :
		num(0) = numb(5); num(1) = numb(4); num(2) = numb(0); num(3) = numb(1);
		num(4) = numb(6); num(5) = numb(7); num(6) = numb(3); num(7) = numb(2); break;
	      case 23 :
		num(0) = numb(4); num(1) = numb(0); num(2) = numb(1); num(3) = numb(5);
		num(4) = numb(7); num(5) = numb(3); num(6) = numb(2); num(7) = numb(6); break;
	      }
	  }
	}
      
      mesh.Element(i).Init(num, ref);
    }
  
}

void CheckSurfaceElement(const ElementReference<Dimension3, 2>& Fb)
{
  bool presence_tri = false, presence_quad = false;
  if (Fb.GetNbVertices() != 8)
    presence_tri = true;

  if (Fb.GetNbVertices() != 4)
    presence_quad = true;
  
  int r = Fb.GetOrder();
  int type_elt = Fb.GetHybridType();
  
  // small mesh with one element
  Mesh<Dimension3> mesh;
  mesh.ReallocateVertices(Fb.GetNbVertices());
  
  switch (type_elt)
    {
    case 0 :
      {
	mesh.Vertex(0).Init(0, 0, 0);
	mesh.Vertex(1).Init(1, 0, 0);
	mesh.Vertex(2).Init(0, 1, 0);
	mesh.Vertex(3).Init(0, 0, 1);
      }
      break;
    case 1 :
      {
	mesh.Vertex(0).Init(-1, -1, 0);
	mesh.Vertex(1).Init(1, -1, 0);
	mesh.Vertex(2).Init(1, 1, 0);
	mesh.Vertex(3).Init(-1, 1, 0);
	mesh.Vertex(4).Init(0, 0, 1);
      }
      break;
    case 2 :
      {
	mesh.Vertex(0).Init(0, 0, 0);
	mesh.Vertex(1).Init(1, 0, 0);
	mesh.Vertex(2).Init(0, 1, 0);
	mesh.Vertex(3).Init(0, 0, 1);
	mesh.Vertex(4).Init(1, 0, 1);
	mesh.Vertex(5).Init(0, 1, 1);
      }
      break;
    case 3 :
      {
	mesh.Vertex(0).Init(0, 0, 0);
	mesh.Vertex(1).Init(1, 0, 0);
	mesh.Vertex(2).Init(1, 1, 0);
	mesh.Vertex(3).Init(0, 1, 0);
	mesh.Vertex(4).Init(0, 0, 1);
	mesh.Vertex(5).Init(1, 0, 1);
	mesh.Vertex(6).Init(1, 1, 1);
	mesh.Vertex(7).Init(0, 1, 1);
      }
    }
  
  mesh.ReallocateElements(1);
  int ref = 1; IVect num(Fb.GetNbVertices());
  num.Fill();
  mesh.Element(0).Init(num, ref);
  mesh.FindConnectivity();
  
  VectR3 sFace; Matrix3_3 DF2, invDF2;
  TinyVector<R3, 3> vec_DF;
  
  if (presence_tri)
    {
      const ElementReference<Dimension2, 2>& tri = Fb.GetTriangularSurfaceFiniteElement();
        
      for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	if (!MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(type_elt, num_loc))
	  {
	    
	    if ( (tri.GetGeometryOrder() != r) || 
		 (tri.GetNbPointsQuadratureInside() != Fb.GetNbQuadBoundary(num_loc)) || 
		 (tri.GetNbDof() != Fb.GetNbDofBoundary(num_loc)) )
	      {
		cout << " GetTriangularSurfaceFiniteElement incorrect " << endl;
		CHECK_ABORT;
	      }
	    
	    /* if (!EqualVector(tri.PointsQuadInsideND(), Fb.Points2D_tri()))
	      {
		cout << " GetTriangularSurfaceFiniteElement incorrect (2-D points different) " << endl;
		DISP(tri.PointsQuadInsideND()); DISP(Fb.Points2D_tri());
		CHECK_ABORT;
	      }
            
	    if (!EqualVector(tri.WeightsND(), Fb.Weights2D_tri()))
	      {
		cout << " GetTriangularSurfaceFiniteElement incorrect (2-D weights different) " << endl;
		DISP(tri.WeightsND()); DISP(Fb.Weights2D_tri());
		CHECK_ABORT;
	      }
            */
            
	    VectR2 points_dof_tri(tri.GetNbPointsDofInside());
	    for (int i = 0; i < points_dof_tri.GetM(); i++)
	      points_dof_tri(i) = tri.PointsDofND(i);
	    
	    if (!EqualVector(points_dof_tri, Fb.PointsDofBoundary(num_loc)))
	      {
		cout << " GetTriangularSurfaceFiniteElement incorrect (2-D points different) " << endl;
		DISP(points_dof_tri); DISP(Fb.PointsDofBoundary(num_loc));
		CHECK_ABORT;
	      }

	    sFace.Reallocate(3);
	    for (int j = 0; j < 3; j++)
	      {
		int nv = MeshNumbering<Dimension3>::
		  GetVertexNumberOfFace(type_elt, num_loc, j);
				
		sFace(j) = mesh.Vertex(nv);
	      }
	    
	    // computation of DF*-1
	    vec_DF(0) = sFace(1) - sFace(0);
	    vec_DF(1) = sFace(2) - sFace(0);
	    TimesProd(vec_DF(0), vec_DF(1), vec_DF(2));
	    DF2 = Matrix3_3(vec_DF);
	    GetInverse(DF2, invDF2);
	    
	    R3 normale_fj = Fb.NormaleLoc(num_loc);
	    R3 pt_loc; VectR2 phi_surf;
	    VectR3 phi; R3 phi_tri, phi3d, vec_u, vec_v;
	    for (int j = 0; j < tri.GetNbPointsQuadratureInside(); j++)
	      {
		Real_wp x = tri.PointsND(j)(0);
		Real_wp y = tri.PointsND(j)(1);
		pt_loc = (1.0-x-y)*sFace(0) + x*sFace(1) + y*sFace(2);

		tri.ComputeValuesPhiRef(tri.PointsND(j), phi_surf);
		Fb.ComputeValuesPhiRef(pt_loc, phi);
		
		for (int i = 0; i < tri.GetNbDof(); i++)
		  {
		    vec_u.Init(phi_surf(i)(0), phi_surf(i)(1), 0);
		    MltTrans(invDF2, vec_u, vec_v);
		    TimesProd(normale_fj, vec_v, phi_tri);

		    int n = Fb.GetLocalNumber(num_loc, i);
		    TimesProd(normale_fj, phi(n), phi3d);
		    
		    // signs may differ...
		    if ( (phi_tri.Distance(phi3d) > 1e-11) && (phi_tri.Distance(-phi3d) > 1e-11))
		      {
			DISP(pt_loc); DISP(i); DISP(num_loc); DISP(j); DISP(phi_tri); DISP(phi3d);
			cout << "Value of basis functions differ in GetTriangularSurfaceFiniteElement" << endl;
			abort();
		      }
		  }
	      }
	  }
    }
  
  if (presence_quad)
    {
      const ElementReference<Dimension2, 2>& quad = Fb.GetQuadrangularSurfaceFiniteElement();
      
      for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	if (MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(type_elt, num_loc))
	  {
	    
	    if ( (quad.GetGeometryOrder() != r) || 
		 (quad.GetNbPointsQuadratureInside() != Fb.GetNbQuadBoundary(num_loc)) || 
		 (quad.GetNbDof() != Fb.GetNbDofBoundary(num_loc)) )
	      {
		cout << " GetQuadrangularSurfaceFiniteElement incorrect " << endl;
		CHECK_ABORT;
	      }

            /*
	    if (!EqualVector(quad.PointsQuadInsideND(), Fb.Points2D_quad()))
	      {
		cout << " GetQuadrangularSurfaceFiniteElement incorrect (2-D points different) " << endl;
		DISP(quad.PointsQuadInsideND()); DISP(Fb.Points2D_quad());
		CHECK_ABORT;
	      }
	    
	    if (!EqualVector(quad.WeightsND(), Fb.Weights2D_quad()))
	      {
		cout << " GetQuadrangularSurfaceFiniteElement incorrect (2-D weights different) " << endl;
		DISP(quad.WeightsND()); DISP(Fb.Weights2D_quad());
		CHECK_ABORT;
	      }
            */
            
	    VectR2 points_dof_quad(quad.GetNbPointsDofInside());
	    for (int i = 0; i < points_dof_quad.GetM(); i++)
	      points_dof_quad(i) = quad.PointsDofND(i);
	    
	    if (!EqualVector(points_dof_quad, Fb.PointsDofBoundary(num_loc)))
	      {
		cout << " GetQuadrangularSurfaceFiniteElement incorrect (2-D points different) " << endl;
		DISP(points_dof_quad); DISP(Fb.PointsDofBoundary(num_loc));
		CHECK_ABORT;
	      }
	    
	    sFace.Reallocate(4);
	    for (int j = 0; j < 4; j++)
	      {
		int nv = MeshNumbering<Dimension3>::
		  GetVertexNumberOfFace(type_elt, num_loc, j);
				
		sFace(j) = mesh.Vertex(nv);
	      }
	    
	    // computation of DF*-1
	    vec_DF(0) = sFace(1) - sFace(0);
	    vec_DF(1) = sFace(3) - sFace(0);
	    TimesProd(vec_DF(0), vec_DF(1), vec_DF(2));
	    DF2 = Matrix3_3(vec_DF);
	    GetInverse(DF2, invDF2);
	    
	    R3 normale_fj = Fb.NormaleLoc(num_loc);
	    R3 pt_loc; VectR2 phi_surf;
	    VectR3 phi; R3 phi_quad, phi3d, vec_u, vec_v;
	    for (int j = 0; j < quad.GetNbPointsQuadratureInside(); j++)
	      {
		Real_wp x = quad.PointsND(j)(0);
		Real_wp y = quad.PointsND(j)(1);
		pt_loc = (1.0-x)*(1.0-y)*sFace(0) + x*(1.0-y)*sFace(1)
		  + x*y*sFace(2) + (1.0-x)*y*sFace(3);

		quad.ComputeValuesPhiRef(quad.PointsND(j), phi_surf);
		Fb.ComputeValuesPhiRef(pt_loc, phi);
		
		for (int i = 0; i < quad.GetNbDof(); i++)
		  {
		    vec_u.Init(phi_surf(i)(0), phi_surf(i)(1), 0);
		    MltTrans(invDF2, vec_u, vec_v);
		    TimesProd(normale_fj, vec_v, phi_quad);

		    int n = Fb.GetLocalNumber(num_loc, i);
		    TimesProd(normale_fj, phi(n), phi3d);
		    
		    // signs may differ...
		    if ( (phi_quad.Distance(phi3d) > 1e-11) && (phi_quad.Distance(-phi3d) > 1e-11))
		      {
			DISP(pt_loc); DISP(i); DISP(num_loc); DISP(j); DISP(phi_quad); DISP(phi3d);
			cout << "Value of basis functions differ in GetQuadrangularSurfaceFiniteElement" << endl;
			abort();
		      }
		  }
	      }
	  }
    }
}

void CheckHcurlFiniteElement(ElementReference<Dimension3, 2>& Fb, const string& name_elt,
			     Mesh<Dimension3>& mesh, VectR3& s, bool test_mass_matrix = false,
			     bool check_surf = true)
{
  // checking Hcurl finite elements
  Real_wp aire = 1.0;
  switch (s.GetM())
    {
    case 4 : aire = 1.0/6; break;
    case 5 : aire = 4.0/3; break;
    case 6 : aire = 0.5; break;
    }
  
  Real_wp threshold = 300.0*epsilon_machine;
  R3::threshold = threshold;
  MeshNumbering<Dimension3> mesh_num(mesh);
  
  int num_loc = 0;
  // main loop on orders of approximations
  for (int r = 1; r <= 4; r++)
    {      
      cout << "Checking order " << r << endl;
      // we can use a different order for geometry
      int rgeom = min(r, 4);

      Fb.ConstructFiniteElement(r, rgeom);
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
      
      Vector<R3> points_nodal = Fb.PointsNodalND();
      for (int k = 0; k < Fb.GetNbVertices(); k++)
        if (points_nodal(k) != s(k))
          {
            cout << "first nodal points have to be vertices " << endl;
            CHECK_ABORT;
          }
      
      R3 ptA, ptB;
      
      int nb_nodes_bounds = (Fb.GetNbVertices()-2)*rgeom*rgeom + 2;
      
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
              TinyVector<Real_wp, 2> tloc = Fb.PointsQuadratureBoundary(k, num_loc);
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
      
      if (r > 1)
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
      VectR3 grad_phi; Real_wp h = pow(epsilon_machine, 1.0/3);
      for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
        {
          Fb.ComputeGradientPhiNodalRef(Fb.PointsNodalND(i), grad_phi);
          Real_wp dphi_num; 
          for (int k = 0; k < 3; k++)
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
      VectR3 vert_elt; R3 pt_loc, ptC; Matrix3_3 dfj, dfjm1;
      mesh.GetVerticesElement(0, vert_elt);
      SetPoints<Dimension3> PointsElem;
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
	  for (int k = 0; k < 3; k++)
	    {
	      pt_loc = s(i); pt_loc(k) += h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptC, mesh, 0);
	      
	      pt_loc = s(i); pt_loc(k) -= h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptB, mesh, 0);
	      
	      for (int m = 0; m < 3; m++)
		{
		  Real_wp df_num = (ptC(m) - ptB(m))/(2.0*h);
		  if (abs(df_num - dfj(m, k)) > sqrt(threshold))
		    {
                      if ((Fb.GetNbVertices() == 5) && (i == 4))
                        {
                          if ((m == 0)&&(k==0))
                            cout << "Warning : DFj not correctly computed on apex of pyramid" << endl;
                        }
                      else
                        {
                          cout << "DFj not correctly computed " << endl;
                          CHECK_ABORT;
                        }
		    }
		}
	    }
	  
	  pt_loc.Fill(0);
	  FjInverseProblem<Dimension3> inverseFj(mesh, 0);
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
      
      SetMatrices<Dimension3> MatricesElem;
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
	  for (int k = 0; k < 3; k++)
	    {
	      pt_loc = s(i); pt_loc(k) += h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptC, mesh, 0);
	      
	      pt_loc = s(i); pt_loc(k) -= h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptB, mesh, 0);
	      
	      for (int m = 0; m < 3; m++)
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
	  FjInverseProblem<Dimension3> inverseFj(mesh, 0);
	  inverseFj.Solve(ptA, pt_loc);
	  
	  if (i != 4)
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
      
      num_loc = 1; Real_wp dsj; R3 normale_fj;
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
	EllipticProblem<MyOwnEquation<Dimension3> > var;
	MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
        
        CreateRegularMesh(var.mesh, s, r, rgeom);
        
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

      //exit(0);
      // checking mass matrix
      Matrix<R3> ValuePhi(Fb.GetNbDof(), Fb.GetNbPointsQuadrature());
      Real_wp max_cij = 0; VectR3 phi, curl_phi;
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
              //CHECK_ABORT;
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
      if (test_mass_matrix)
        {
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              Ones(i) = 1.0; MhOnes = Ones;
              Fb.MltMassMatrix(MhOnes);
              for (int j = 0; j < Fb.GetNbDof(); j++)
                if (abs(MhOnes(j) - Mh(j, i))/max_mij > 100*threshold)
                  {
                    //cout << "MltMassMatrix not correct " << endl;
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
                    //cout << "SolveCholesky not correct " << endl;
                    solve_cholesky_correct = false;
                    CHECK_ABORT;
                  }
              
              MhOnes = Ones;
              Fb.SolveMassMatrix(MhOnes);
              for (int j = 0; j < Fb.GetNbDof(); j++)
                if (abs(InvMh(j, i) - MhOnes(j))/max_invMij > 100*threshold)
                  {
                    //cout << "SolveMassMatrix not correct " << endl;
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
          
	  mesh.Element(0).UnsetCurved();

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
        }
          
      // checking ApplyCh/ComputeIntegralRef (same methods)
      Ones.Reallocate(3*Fb.GetNbPointsQuadratureInside()); Ones.Fill(0);
      for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
        for (int n = 0; n < 3; n++)
          {
            Ones(3*j+n) = 1.0;
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
            
            Ones(3*j+n) = 0.0;
          }
      
      // and ApplyChTranspose
      MhOnes.Fill(0);
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          MhOnes(j) = 1.0;
          Fb.ApplyChTranspose(MhOnes, Ones);
          for (int i = 0; i < Fb.GetNbPointsQuadratureInside(); i++)
            for (int n = 0; n < 3; n++)
              if (abs(Ones(3*i+n) - ValuePhi(j, i)(n))/max_cij > threshold)
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
      for (int comp = 0; comp < 3; comp++)
        {
          for (int m = 0; m <= r-1; m++)
            for (int n = 0; n <= r-1-m; n++)
              for (int k = 0; k <= r-1-m-n; k++)
                {
                  Ones.Fill(0);
                  for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
                    {
                      Real_wp x = Fb.PointsND(j)(0);
                      Real_wp y = Fb.PointsND(j)(1);
                      Real_wp z = Fb.PointsND(j)(2);
                      Ones(3*j+comp) = pow(x, m)*pow(y, n)*pow(z, k);
                    }
                  
                  Fb.ProjectQuadratureToDofRef(Ones, MhOnes);
                  
                  for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
                    Ones(3*j+comp) *= Fb.WeightsND(j);
                  
                  Fb.ComputeIntegralRef(Ones, ProjOnes);
                  Fb.SolveMassMatrix(ProjOnes);
                  for (int i = 0; i < Fb.GetNbDof(); i++)
                    if (abs(MhOnes(i) - ProjOnes(i)) > threshold)
                      {
                        proj_quadrature_to_dof_correct = false;
                        cout << "ProjectQuadratureToDofRef incorrect" << endl;
                        DISP(i); DISP(MhOnes(i)); DISP(ProjOnes(i));
                        CHECK_ABORT;
                      }
                }
        }
      
      if (!proj_quadrature_to_dof_correct)
        cout << "Warning : ProjectQuadratureToDofRef not correct " << endl;
      
      // checking ComputeProjectionDofRef
      VectReal_wp feval(3*Fb.GetNbPointsDofInside());
      VectComplex_wp feval_cplx(3*Fb.GetNbPointsDofInside());
      VectComplex_wp ProjOnes_cplx(Fb.GetNbDof()), Ones_cplx(3*Fb.GetNbPointsQuadratureInside());
      ProjOnes_cplx.Fill(0); Ones_cplx.Fill(0);
      for (int comp = 0; comp < 3; comp++)
        {
          for (int m = 0; m <= r-1; m++)
            for (int n = 0; n <= r-1-m; n++)
              for (int k = 0; k <= r-1-m-n; k++)
                {
                  feval.Fill(0); feval_cplx.Fill(0);
                  for (int j = 0; j < Fb.GetNbPointsDofInside(); j++)
                    {
                      Real_wp x = Fb.PointsDofND(j)(0);
                      Real_wp y = Fb.PointsDofND(j)(1);
                      Real_wp z = Fb.PointsDofND(j)(2);
                      feval(3*j+comp) = pow(x, m)*pow(y, n)*pow(z, k);
                      feval_cplx(3*j+comp) = feval(3*j+comp);
                    }
		  
		  int comp1 = 0, comp2 = 1;
		  if (comp == 0)
		    { comp1 = 1; comp2 = 2; }
		  else if (comp == 1)
		    { comp1 = 0; comp2 = 2; }
		  
                  Fb.ComputeProjectionDofRef(feval, ProjOnes);
                  Fb.ApplyChTranspose(ProjOnes, Ones);
                  for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
                    {
                      Real_wp x = Fb.PointsND(j)(0);
                      Real_wp y = Fb.PointsND(j)(1);
                      Real_wp z = Fb.PointsND(j)(2);
                      Real_wp val_exact = pow(x, m)*pow(y, n)*pow(z, k);
                      if ((abs(val_exact - Ones(3*j+comp)) > 10.0*r*r*r*threshold)
			|| (abs(Ones(3*j+comp1)) > 10.0*r*r*r*threshold)
			  || (abs(Ones(3*j+comp2)) > 10.0*r*r*r*threshold) )
                        {
                          cout << "ComputeProjectionDofRef not correct " << endl;
                          DISP(j); DISP(m); DISP(n); DISP(k);
                          DISP(comp); DISP(val_exact); DISP(Ones(3*j+comp));
                          CHECK_ABORT;
                        }
                    }

                  Fb.ComputeProjectionDofRef(feval_cplx, ProjOnes_cplx);
                  Fb.ApplyChTranspose(ProjOnes_cplx, Ones_cplx);
                  for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
                    {
                      Real_wp x = Fb.PointsND(j)(0);
                      Real_wp y = Fb.PointsND(j)(1);
                      Real_wp z = Fb.PointsND(j)(2);
                      Real_wp val_exact = pow(x, m)*pow(y, n)*pow(z, k);
                      if ((abs(val_exact - Ones_cplx(3*j+comp)) > 10.0*r*r*r*threshold)
			|| (abs(Ones_cplx(3*j+comp1)) > 10.0*r*r*r*threshold)
			  || (abs(Ones_cplx(3*j+comp2)) > 10.0*r*r*r*threshold) )
                        {
                          cout << "ComputeProjectionDofRef not correct " << endl;
                          DISP(j); DISP(m); DISP(n); DISP(k);
                          DISP(comp); DISP(val_exact); DISP(Ones_cplx(3*j+comp));
                          CHECK_ABORT;
                        }
                    }

                }
        }
      
      // checking ComputeNodalValuesRef and ComputeNodalGradientRef
      Matrix<R3> ValuePhiNodal(Fb.GetNbDof(), Fb.GetNbPointsNodalElt());
      Matrix<R3> GradPhiNodal(Fb.GetNbPointsNodalElt(), Fb.GetNbPointsNodalElt());
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
      MhOnes.Reallocate(3*Fb.GetNbPointsNodalElt()); MhOnes.Fill(0);
      VectReal_wp Ex(Fb.GetNbPointsNodalElt());
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          Ones(j) = 1.0;
          Fb.ComputeNodalValuesRef(Ones, MhOnes);
          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            if (i != 4)
	      for (int n = 0; n < 3; n++)
		if (abs(ValuePhiNodal(j, i)(n) - MhOnes(3*i+n)) > threshold)
		  {
		    cout << "ComputeNodalValuesRef not correct " << endl;
		    DISP(MhOnes(3*i+n)); DISP(ValuePhiNodal(j, i)(n));
		    CHECK_ABORT;
		  }
          
          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            Ex(i) = MhOnes(3*i);
          
          Fb.ComputeNodalGradientRef(Ex, grad_phi);          
          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            if (i != 4)
	      {
		R3 grad;
		for (int k = 0; k < Fb.GetNbPointsNodalElt(); k++)
		  Add(MhOnes(3*k), GradPhiNodal(k, i), grad);
		
		if (grad != grad_phi(i))
		  {
		    cout << "ComputeNodalGradientRef not correct " << endl;
		    CHECK_ABORT;
		  }
	      }

          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            Ex(i) = MhOnes(3*i+1);
          
          Fb.ComputeNodalGradientRef(Ex, grad_phi);          
          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            if (i != 4)
	      {
		R3 grad;
		for (int k = 0; k < Fb.GetNbPointsNodalElt(); k++)
		  Add(MhOnes(3*k+1), GradPhiNodal(k, i), grad);
		
		if (grad != grad_phi(i))
		  {
		    cout << "ComputeNodalGradientRef not correct " << endl;
		    CHECK_ABORT;
		  }
	      }

          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            Ex(i) = MhOnes(3*i+2);
          
          Fb.ComputeNodalGradientRef(Ex, grad_phi);          
          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            if (i != 4)
	      {
		R3 grad;
		for (int k = 0; k < Fb.GetNbPointsNodalElt(); k++)
		  Add(MhOnes(3*k+2), GradPhiNodal(k, i), grad);
		
		if (grad != grad_phi(i))
		  {
		    cout << "ComputeNodalGradientRef not correct " << endl;
		    CHECK_ABORT;
		  }
	      }
          
          Ones(j) = 0.0;
        }
      
      
      // checking curl of basis functions
      VectR3 phi_p, phi_m;
      for (int i = 0; i < Fb.GetNbPointsDof(); i++)
        {
          Fb.ComputeCurlPhiRef(Fb.PointsDofND(i), curl_phi);
          Real_wp val_max = 0;
          for (int j = 0; j < Fb.GetNbDof(); j++)
            val_max = max(val_max, Norm2(curl_phi(j)));

          VectR3 curlphi_num(Fb.GetNbDof()); 
          FillZero(curlphi_num);
          
          ptA = Fb.PointsDofND(i);
          ptA(0) += h;
          Fb.ComputeValuesPhiRef(ptA, phi_p);
          
          ptA(0) -= 2.0*h;
          Fb.ComputeValuesPhiRef(ptA, phi_m);
          
          for (int j = 0; j < Fb.GetNbDof(); j++)
            {
              curlphi_num(j)(1) -= (phi_p(j)(2) - phi_m(j)(2))/(2.0*h);
              curlphi_num(j)(2) += (phi_p(j)(1) - phi_m(j)(1))/(2.0*h);
            }
          
          ptA = Fb.PointsDofND(i);
          ptA(1) += h;
          Fb.ComputeValuesPhiRef(ptA, phi_p);
          
          ptA(1) -= 2.0*h;
          Fb.ComputeValuesPhiRef(ptA, phi_m);
          
          for (int j = 0; j < Fb.GetNbDof(); j++)
            {
              curlphi_num(j)(0) += (phi_p(j)(2) - phi_m(j)(2))/(2.0*h);
              curlphi_num(j)(2) -= (phi_p(j)(0) - phi_m(j)(0))/(2.0*h);
            }
          
          ptA = Fb.PointsDofND(i);
          ptA(2) += h;
          Fb.ComputeValuesPhiRef(ptA, phi_p);
          
          ptA(2) -= 2.0*h;
          Fb.ComputeValuesPhiRef(ptA, phi_m);
          
          for (int j = 0; j < Fb.GetNbDof(); j++)
            {
              curlphi_num(j)(0) -= (phi_p(j)(1) - phi_m(j)(1))/(2.0*h);
              curlphi_num(j)(1) += (phi_p(j)(0) - phi_m(j)(0))/(2.0*h);
            }
          
          for (int j = 0; j < Fb.GetNbDof(); j++)
            {
              if ( curlphi_num(j).Distance(curl_phi(j)) > val_max*sqrt(threshold))
                {
                  DISP(Fb.PointsDofND(i)); DISP(j); DISP(curlphi_num(j)); DISP(curl_phi(j));
                  cout << "Curl of basis functions not correct " <<endl;
                  CHECK_ABORT;
                }
            }
        }

      // checking GetValuePhiOnQuadraturePoint (and Gradient)
      Matrix<R3> CurlPhi(Fb.GetNbDof(), Fb.GetNbPointsQuadrature());
      Real_wp max_rij = 0;
      for (int k = 0; k < Fb.GetNbPointsQuadrature(); k++)
        {
          Fb.ComputeCurlPhiRef(Fb.PointsND(k), curl_phi);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              CurlPhi(i, k) = curl_phi(i);
              max_rij = max(max_rij, Norm2(CurlPhi(i, k)));
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
            if (curl_phi(i) != CurlPhi(i, k))
              {
                cout << "GetCurlPhiOnQuadraturePoint not correct " << endl;
                DISP(curl_phi(i)); DISP(CurlPhi(i, k));
                CHECK_ABORT;
              }
        }
      
      // checking ComputeIntegralGradientRef / ApplyRh
      Ones.Reallocate(3*Fb.GetNbPointsQuadratureInside());
      MhOnes.Reallocate(Fb.GetNbDof()); MhOnes.Fill(0);
      Ones.Fill(0);
      for (int j = 0; j < 3*Fb.GetNbPointsQuadratureInside(); j++)
        {
          int j2 = j/3;
          int m = j%3;
          Ones(j) = 1.0;
          Fb.ComputeIntegralCurlRef(Ones, MhOnes);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(MhOnes(i) - CurlPhi(i, j2)(m))/max_rij > threshold)
              {
                cout << "ComputeIntegralCurlRef not correct " << endl;
                CHECK_ABORT;
              }

          Fb.ApplyRh(Ones, MhOnes);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(MhOnes(i) - CurlPhi(i, j2)(m))/max_rij > threshold)
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
            for (int m = 0; m < 3; m++)
              if (abs(Ones(j*3+m) - CurlPhi(i, j)(m))/max_rij > threshold)
                {
                  cout << "ApplyRhTranspose not correct " << endl;
                  CHECK_ABORT;
                }
          
          MhOnes(i) = 0.0;
        }
      
      // checking ComputeIntegralSurfaceRef/ApplySh
      ComputeGaussLegendre(points1d, weights1d, r);
      VectR2 points2d_tri, points2d_quad;
      VectReal_wp weights2d_tri, weights2d_quad;
      TriangleQuadrature::ConstructQuadrature(2*r, points2d_tri, weights2d_tri);
      Matrix<int> NumNodes2D_quad, coor;
      MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(r, NumNodes2D_quad, coor);
      weights2d_quad.Reallocate((r+1)*(r+1));
      points2d_quad.Reallocate((r+1)*(r+1));
      for (int i = 0; i <= r; i++)
        for (int j = 0; j <= r; j++)
          {
            points2d_quad(NumNodes2D_quad(i, j)).Init(points1d(i), points1d(j));
            weights2d_quad(NumNodes2D_quad(i, j)) = weights1d(i)*weights1d(j);
          }
      
      Matrix<bool> IsDofOnBoundary(Fb.GetNbBoundaries(), Fb.GetNbDof());
      IsDofOnBoundary.Fill(false);

      // IVect num_point_quad;
      for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
        {
          int nb_pts_quad = Fb.GetNbQuadBoundary(num_loc);
          Ones.Reallocate(3*nb_pts_quad); Ones.Fill(0);
          for (int k = 0; k < nb_pts_quad; k++)
            for (int n = 0; n < 3; n++)
              {
                Ones(3*k+n) = 1.0; MhOnes.Fill(1.0);
                int k2 = Fb.GetQuadNumber(num_loc, k);
                
                Fb.ComputeIntegralSurfaceRef(Ones, MhOnes, num_loc);
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  if (abs(MhOnes(i) - ValuePhi(i, k2)(n))/max_cij > threshold)
                    {
                      DISP(MhOnes(i)); DISP(ValuePhi(i, k2)(n)); DISP(i); DISP(num_loc); DISP(k2); DISP(n);
                      cout << "ComputeIntegralSurfaceRef not correct " << endl;
                      CHECK_ABORT;
                    }
                
                MhOnes.Fill(1.5);
                Fb.ApplySh(2.5, num_loc, Ones, MhOnes);
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  if (abs(MhOnes(i) - 2.5*ValuePhi(i, k2)(n) - 1.5)/max_cij > threshold)
                    {
                      DISP(MhOnes(i)); DISP(2.5*ValuePhi(i, k2)(n)+1.5);
                      DISP(i); DISP(num_loc); DISP(k2); DISP(n);
                      cout << "ApplySh not correct " << endl;
                      CHECK_ABORT;
                    }
                
                Ones(3*k+n) = 0.0;
              }
          
          Ones.Fill(0);
          for (int k = 0; k < nb_pts_quad; k++)
            for (int m = 0; m < 3; m++)
              {
                Ones(k*3+m) = 1.0;
                int k2 = Fb.GetQuadNumber(num_loc, k);
                Fb.ComputeIntegralSurfaceCurlRef(Ones, MhOnes, num_loc);
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  if (abs(MhOnes(i) - CurlPhi(i, k2)(m))/max_rij > threshold)
                    {
                      cout << "ComputeIntegralSurfaceCurlRef not correct " << endl;
                      DISP(MhOnes(i)); DISP(CurlPhi(i, k2)(m)); DISP(num_loc);
                      CHECK_ABORT;
                    }
                
                
                MhOnes.Fill(1.5);
                Fb.ApplyNablaSh(2.5, num_loc, Ones, MhOnes);
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  if (abs(MhOnes(i) - 2.5*CurlPhi(i, k2)(m) - 1.5)/max_rij > threshold)
                    {
                      cout << "ApplyNablaSh not correct " << endl;
                      DISP(MhOnes(i)); DISP(2.5*CurlPhi(i, k2)(m) + 1.5); DISP(num_loc);
                      CHECK_ABORT;
                    }
                
                Ones(k*3+m) = 0.0;
              }
          
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              MhOnes.Fill(0); MhOnes(i) = 1.0;
                              
              Fb.ApplyShTranspose(num_loc, MhOnes, Ones);
              for (int k = 0; k < nb_pts_quad; k++)
                {
                  int k2 = Fb.GetQuadNumber(num_loc, k);
                  for (int n = 0; n < 3; n++)
                    if (abs(Ones(3*k+n) - ValuePhi(i, k2)(n))/max_cij > threshold)
                      {
                        DISP(Ones(3*k+n)); DISP(ValuePhi(i, k2)(n));
                        DISP(i); DISP(num_loc); DISP(k2); DISP(n);
                        cout << "ApplyShTranspose not correct " << endl;
                        CHECK_ABORT;
                      }
                }
          
              Fb.ApplyNablaShTranspose(num_loc, MhOnes, Ones);
              for (int k = 0; k < nb_pts_quad; k++)
                {
                  int k2 = Fb.GetQuadNumber(num_loc, k);
                  for (int n = 0; n < 3; n++)
                    if (abs(Ones(3*k+n) - CurlPhi(i, k2)(n))/max_rij > threshold)
                      {
                        cout << "ApplyNablaShTranspose not correct " << endl;
                        DISP(i); DISP(k); DISP(n);
                        DISP(Ones(3*k+n)); DISP(CurlPhi(i, k2)(n)); DISP(num_loc);
                        CHECK_ABORT;
                      }
                }
            }
          
          nb_pts_quad = points1d.GetM();
          VectR2 pts_integration; VectReal_wp poids_integration;
          if (MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(Fb.GetHybridType(), num_loc))
            {
              pts_integration = points2d_quad;
              poids_integration = weights2d_quad;
              nb_pts_quad = points2d_quad.GetM();
            }
          else
            {
              pts_integration = points2d_tri;
              poids_integration = weights2d_tri;
              nb_pts_quad = points2d_tri.GetM();
            }
          
          Ones.Reallocate(3*nb_pts_quad); Ones.Fill(0);
          for (int k = 0; k < nb_pts_quad; k++)
            for (int n = 0; n < 3; n++)
              {
                Ones(3*k+n) = 1.0; MhOnes.Fill(1.0);
                TinyVector<Real_wp, 2> xi;
              
                xi(0) = pts_integration(k)(0);
                xi(1) = pts_integration(k)(1);
                
                Fb.GetLocalCoordOnBoundary(num_loc, xi, pt_loc);
                Fb.ComputeValuesPhiRef(pt_loc, phi);
              
                Fb.ComputeGaussIntegralSurfaceRef(Ones, MhOnes, num_loc);
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  if (abs(MhOnes(i) - phi(i)(n))/max_cij > threshold)
                    {
                      cout << "ComputeGaussIntegralSurfaceRef not correct " << endl;
                      CHECK_ABORT;
                    }
                
                Ones(3*k+n) = 0.0;

		// detecting surface dofs
		for (int i = 0; i < Fb.GetNbDof(); i++)
		  {
		    R3 phi_times_n;
		    TimesProd(phi(i), Fb.NormaleLoc(num_loc), phi_times_n);
		    
		    if (Norm2(phi_times_n) > 10*threshold)
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
	    if ((Fb.GetNbPointsDofInside() == 0) || (Fb.GetNbPointsDofInside() < Fb.GetNbDof()/3))
	      {
		cout << "GetNbPointsDofInside() incorrect" << endl;
		DISP(Fb.GetNbPointsDofInside());
		abort();
	      }

	    if ((Fb.GetNbPointsDofSurface(num_loc) == 0)
		|| (Fb.GetNbPointsDofSurface(num_loc) < Fb.GetNbDofBoundary(num_loc)/2))
	      {
		cout << "GetNbPointsDofSurface incorrect" << endl;
		DISP(Fb.GetNbPointsDofSurface(num_loc));
		abort();
	      }
	    
	    const VectR2& pts_face = Fb.PointsDofBoundary(num_loc);
	    if (pts_face.GetM() != Fb.GetNbPointsDofSurface(num_loc))
	      {
		cout << "PointsDofBoundary incorrect" << endl;
		abort();
	      }
	    
	    for (int i = 0; i < Fb.GetNbPointsDofSurface(num_loc); i++)
	      {
		R3 pt_loc;
		Fb.GetLocalCoordOnBoundary(num_loc, pts_face(i), pt_loc);
		if (pt_loc != Fb.PointsDofND(Fb.GetPointDofNumber(num_loc, i)))
		  {
		    DISP(num_loc); DISP(i);
		    DISP(pt_loc); DISP(Fb.PointsDofND(Fb.GetPointDofNumber(num_loc, i)));
		    cout << "PointsDofBoundary incorrect" << endl;
		    abort();
		  }
	      }
	    
	    for (int p = 0; p < 3; p++)
	      for (int m = 0; m <= r-1; m++)
		for (int n = 0; n <= r-1-m; n++)
		  {
		    Vector<Real_wp> feval_surf(3*Fb.GetNbPointsDofSurface(num_loc));
		    feval_surf.Zero();
		    for (int i = 0; i < Fb.GetNbPointsDofSurface(num_loc); i++)
		      feval_surf(3*i+p) = FuncPowS(Fb.PointsDofND(Fb.GetPointDofNumber(num_loc, i)), m, n);
		    
		    VectReal_wp contrib_face(Fb.GetNbDofBoundary(num_loc));
		    Fb.ComputeProjectionSurfaceDofRef(feval_surf, contrib_face, num_loc);
		    
		    VectReal_wp feval_vol(3*Fb.GetNbPointsDofInside());
		    feval_vol.Zero();
		    for (int i = 0; i < feval_vol.GetM()/3; i++)
		      feval_vol(3*i+p) = FuncPowS(Fb.PointsDofND(i), m, n);
		    
		    Fb.ComputeProjectionDofRef(feval_vol, ProjOnes);
		    for (int i = 0; i < Fb.GetNbDofBoundary(num_loc); i++)
		      if (abs(contrib_face(i) - ProjOnes(Fb.GetLocalNumber(num_loc, i))) > 100*r*r*r*threshold)
			{
			  cout << "ComputeProjectionSurfaceDofRef not correct " << endl;
			  DISP(num_loc); DISP(i); DISP(m); DISP(n); DISP(p);
			  DISP(contrib_face(i));
			  DISP(ProjOnes(Fb.GetLocalNumber(num_loc, i)));
			  DISP(abs(contrib_face(i) - ProjOnes(Fb.GetLocalNumber(num_loc, i))));
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
      
      // checking ApplyRhQuadrature/ApplyRhQuadratureTranspose (when tensorized points are used)
      if (Fb.UseQuadraturePointsForRh())
        {
          int nb_dof_loc = Fb.GetNbDof();
          int nb_points_quadrature_inside = Fb.GetNbPointsQuadratureInside();
          VectReal_wp Vh(3*nb_points_quadrature_inside);
          VectReal_wp AhUh(nb_dof_loc);
          VectReal_wp Uh(3*nb_points_quadrature_inside);
          VectReal_wp AhUh2(nb_dof_loc);
          for (int i = 0; i < Uh.GetM(); i++)
            {
              Uh.Fill(0); Uh(i) = 1.0;
              Fb.ApplyRh(Uh, AhUh2);
              Fb.ApplyRhQuadrature(Uh, Vh);
              Fb.ApplyCh(Vh, AhUh);
              for (int j = 0; j < AhUh.GetM(); j++)
                if (abs(AhUh(j) - AhUh2(j))/max_rij > threshold)
                  {
                    DISP(AhUh(j));
                    DISP(AhUh2(j));
                    cout << "ApplyRhQuadrature not correct " << endl;
                    CHECK_ABORT;
                  }
            }
                    
          AhUh.Reallocate(3*nb_points_quadrature_inside);
          Uh.Reallocate(nb_dof_loc);
          AhUh2.Reallocate(3*nb_points_quadrature_inside);
          for (int i = 0; i < Uh.GetM(); i++)
            {
              Uh.Fill(0); Uh(i) = 1.0;
              Fb.ApplyRhTranspose(Uh, AhUh);
              Fb.ApplyChTranspose(Uh, Vh);
              Fb.ApplyRhQuadratureTranspose(Vh, AhUh2);
              for (int j = 0; j < AhUh.GetM(); j++)
                if (abs(AhUh(j) - AhUh2(j))/max_rij > threshold)
                  {
                    DISP(AhUh(j));
                    DISP(AhUh2(j));
                    cout << "ApplyRhQuadratureTranspose not correct " << endl;
                    CHECK_ABORT;
                  }
            }
          
        }
      
      // checking addition of elementary matrices
      int nb_dof_loc = Fb.GetNbDof(); int off_row = 12, off_col = 30;
      Matrix3_3 Bmass; Bmass(1,1) = 1.2; Bmass.FillRand(); Mlt(1e-9, Bmass);
      Bmass(1, 0) = Bmass(0, 1); Bmass(2, 0) = Bmass(0, 2); Bmass(2, 1) = Bmass(1, 2);
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
      Matrix3_3 Cgrad;      
      Cgrad.FillRand(); Mlt(1e-9, Cgrad);
      Cgrad(1, 0) = Cgrad(0, 1); Cgrad(2, 0) = Cgrad(0, 2); Cgrad(2, 1) = Cgrad(1, 2);
      Fb.AddConstantStiffnessMatrix(off_row, off_col, Cgrad, AhNew);
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val_exact = 0;
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              val_exact += DotProd(CurlPhi(i, k), dot(Cgrad, CurlPhi(j, k)))*Fb.WeightsND(k);
            
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_rij > threshold)
              {
                cout << "AddConstantStiffnessMatrix not correct " << endl;
                DISP(i); DISP(j); DISP(AhNew(i, j) - AhOrig(i, j)); DISP(val_exact);
                CHECK_ABORT;
              }
          }
            
      // checking ValueBoundaryRef, and GradientBoundaryRef
      Ones.Reallocate(Fb.GetNbDof()); Ones.Fill(0);
      MhOnes.Reallocate(3*Fb.GetNbPointsNodalElt()); MhOnes.Fill(0);
      
      VectReal_wp feval_surf;
      
      bool compute_gradient_boundary = true;
      bool compute_value_boundary = true;
      for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
        {
          feval_surf.Reallocate(3*Fb.GetNbNodalBoundary(num_loc));
          feval.Reallocate(3*Fb.GetNbNodalBoundary(num_loc));
          int nb_pts_nodal = Fb.GetNbNodalBoundary(num_loc);
          Matrix<R3> E_nodal(Fb.GetNbDof(), nb_pts_nodal), curlE_nodal(Fb.GetNbDof(), nb_pts_nodal);
          for (int k = 0; k < nb_pts_nodal; k++)
            {
              int j = Fb.GetNodalNumber(num_loc, k);
              Fb.ComputeValuesPhiRef(Fb.PointsNodalND(j), phi_p);
              Fb.ComputeCurlPhiRef(Fb.PointsNodalND(j), curl_phi);
              for (int i = 0; i < Fb.GetNbDof(); i++)
                {
                  E_nodal(i, k) = phi_p(i);
                  curlE_nodal(i, k) = curl_phi(i);
                }
            }
          
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              Ones(i) = 1.0;
              
              Fb.ComputeValueBoundaryRef(Ones, feval_surf, num_loc);
              Fb.ComputeCurlBoundaryRef(Ones, feval, num_loc);
                
              for (int k = 0; k < nb_pts_nodal; k++)
                {
                  if ((abs(feval_surf(3*k) - E_nodal(i, k)(0)) > threshold) 
                      || (abs(feval_surf(3*k+1) - E_nodal(i, k)(1)) > threshold) 
                      || (abs(feval_surf(3*k+2) - E_nodal(i, k)(2)) > threshold))
                    {
                      compute_value_boundary = false;
                      //cout << "ComputeValueBoundaryRef not correct " << endl;
		      //DISP(feval_surf(3*k)); DISP(feval_surf(3*k+1)); DISP(feval_surf(3*k+2));
		      //DISP(E_nodal(i, k)(0)); DISP(E_nodal(i, k)(1)); DISP(E_nodal(i, k)(2));
                      //CHECK_ABORT;
                    }
                  
                  if ((abs(feval(3*k) - curlE_nodal(i, k)(0)) > 20*threshold) 
                      || (abs(feval(3*k+1) - curlE_nodal(i, k)(1)) > 20*threshold) 
                      || (abs(feval(3*k+2) - curlE_nodal(i, k)(2)) > 20*threshold))
                    {
                      compute_gradient_boundary = false;
                      //cout << "ComputeCurlBoundaryRef not correct " << endl;
		      //DISP(feval_surf(3*k)); DISP(feval_surf(3*k+1)); DISP(feval_surf(3*k+2));
		      //DISP(curlE_nodal(i, k)(0)); DISP(curlE_nodal(i, k)(1));
		      //DISP(curlE_nodal(i, k)(2));
                      //CHECK_ABORT;
                    }
                }
              Ones(i) = 0.0;
            }
        }
      
      if (!compute_value_boundary)
        {
          cout << "Warning : ComputeValueBoundaryRef incorrect" << endl;
        }

      if (!compute_gradient_boundary)
        {
          cout << "Warning : ComputeGradientBoundaryRef incorrect" << endl;
        }
            
    }
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  // 3-D mesh with only one tetrahedron
  Mesh<Dimension3> mesh3d;
  mesh3d.SetGeometryOrder(4);
  
  mesh3d.ReallocateVertices(4);
  mesh3d.Vertex(0).Init(-0.2, 0.2, 0.1);
  mesh3d.Vertex(1).Init(1, 0, 0);
  mesh3d.Vertex(2).Init(0.5, 0.5*sqrt(3.0), 0);
  mesh3d.Vertex(3).Init(0.2, 0.6, sqrt(0.6));
  
  mesh3d.ReallocateBoundariesRef(4);
  mesh3d.BoundaryRef(0).InitTriangular(0, 1, 2, 1);
  mesh3d.BoundaryRef(1).InitTriangular(1, 2, 3, 2);
  mesh3d.BoundaryRef(2).InitTriangular(0, 1, 3, 1);
  mesh3d.BoundaryRef(3).InitTriangular(0, 2, 3, 1);
  
  mesh3d.ReallocateElements(1);
  mesh3d.Element(0).InitTetrahedral(0, 1, 2, 3, 1);
  
  mesh3d.SetCurveType(2, mesh3d.CURVE_SPHERE);
  VectReal_wp param; param.Reallocate(4);
  param(0) = 0; param(1) = 0; param(2) = 0.0; param(3) = 1.0;
  mesh3d.SetCurveParameter(2, param);
  
  mesh3d.ReorientElements();
  mesh3d.FindConnectivity();
  mesh3d.ProjectPointsOnCurves();
  
  Vector<R3> vertex_tet(4);
  vertex_tet(0).Init(0, 0, 0);
  vertex_tet(1).Init(1, 0, 0);
  vertex_tet(2).Init(0, 1, 0);
  vertex_tet(3).Init(0, 0, 1);
  
  {
    cout << "Testing TetrahedronHcurlFirstFamily " << endl;
    TetrahedronHcurlFirstFamily tet; string name_elt("TETRAHEDRON_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_tet);
  }

  {
    cout << "Testing TetrahedronHcurlOptimalFirstFamily " << endl;
    TetrahedronHcurlOptimalFirstFamily tet; string name_elt("TETRAHEDRON_OPTIMAL_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_tet);
  }
	
  {
    cout << "Testing TetrahedronHcurlOptimalHpFirstFamily " << endl;
    TetrahedronHcurlOptimalHpFirstFamily tet; string name_elt("TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_tet);
  }
  
  mesh3d.Clear();
  mesh3d.ReallocateVertices(5);
  mesh3d.Vertex(0).Init(1, 0, 0);
  mesh3d.Vertex(1).Init(0.5, 0.5*sqrt(3.0), 0);
  mesh3d.Vertex(2).Init(0.8, 0.2, sqrt(0.68));
  mesh3d.Vertex(3).Init(0.7, -0.3, sqrt(0.58));
  mesh3d.Vertex(4).Init(-0.2, 0.2, 0.1);

  mesh3d.ReallocateBoundariesRef(5);
  mesh3d.BoundaryRef(0).InitQuadrangular(0, 1, 2, 3, 2);
  mesh3d.BoundaryRef(1).InitTriangular(0, 1, 4, 1);
  mesh3d.BoundaryRef(2).InitTriangular(1, 2, 4, 1);
  mesh3d.BoundaryRef(3).InitTriangular(2, 3, 4, 1);
  mesh3d.BoundaryRef(4).InitTriangular(0, 3, 4, 1);
  
  mesh3d.ReallocateElements(1);
  mesh3d.Element(0).InitPyramidal(0, 1, 2, 3, 4, 1);
  
  mesh3d.ReorientElements();
  mesh3d.FindConnectivity();
  mesh3d.ProjectPointsOnCurves();

  Vector<R3> vertex_pyr(5);
  vertex_pyr(0).Init(-1, -1, 0);
  vertex_pyr(1).Init(1, -1, 0);
  vertex_pyr(2).Init(1, 1, 0);
  vertex_pyr(3).Init(-1, 1, 0);
  vertex_pyr(4).Init(0, 0, 1);
  
  {
    cout << "Testing PyramidHcurlFirstFamily " << endl;
    PyramidHcurlFirstFamily tet; string name_elt("TETRAHEDRON_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_pyr);
  }
  
  {
    cout << "Testing PyramidHcurlOptimalFirstFamily " << endl;
    PyramidHcurlOptimalFirstFamily tet; string name_elt("TETRAHEDRON_OPTIMAL_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_pyr);
  }

  {
    cout << "Testing PyramidHcurlOptimalHpFirstFamily " << endl;
    PyramidHcurlOptimalHpFirstFamily tet; string name_elt("TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_pyr);
  }

  {
    cout << "Testing PyramidHcurlHpFirstFamily " << endl;
    PyramidHcurlHpFirstFamily tet; string name_elt("TETRAHEDRON_HP_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_pyr);
  }
  
  {
    //cout << "Testing PyramidHcurlOther " << endl;
    //PyramidHcurlOther tet;
    //CheckHcurlFiniteElement(tet, mesh3d, vertex_pyr);
  }
  
  mesh3d.Clear();
  mesh3d.ReallocateVertices(6);
  mesh3d.Vertex(0).Init(0.0, 0.0, 0.0);
  mesh3d.Vertex(1).Init(1.0, 0.0, 0.0);
  mesh3d.Vertex(2).Init(0.0, 1.0, 0.0);
  //mesh3d.Vertex(3).Init(0.0, 0.0, 1.0);
  //mesh3d.Vertex(4).Init(1.0, 0.0, 1.0);
  //mesh3d.Vertex(5).Init(0.0, 1.0, 1.0);
  mesh3d.Vertex(3).Init(0.0, 0.0, 0.32);
  mesh3d.Vertex(4).Init(sqrt(1.0-0.32*0.32), 0, 0.32);
  mesh3d.Vertex(5).Init(0.0, sqrt(1.0-0.32*0.32), 0.32);
  
  mesh3d.ReallocateBoundariesRef(5);
  mesh3d.BoundaryRef(0).InitTriangular(0, 1, 2, 1);
  mesh3d.BoundaryRef(1).InitQuadrangular(0, 1, 4, 3, 1);
  mesh3d.BoundaryRef(2).InitQuadrangular(1, 2, 5, 4, 2);
  mesh3d.BoundaryRef(3).InitQuadrangular(0, 2, 5, 3, 1);
  mesh3d.BoundaryRef(4).InitTriangular(3, 4, 5, 1);
  
  mesh3d.ReallocateElements(1);
  mesh3d.Element(0).InitWedge(0, 1, 2, 3, 4, 5, 1);
  
  mesh3d.ReorientElements();
  mesh3d.FindConnectivity();
  mesh3d.ProjectPointsOnCurves();
  
  Vector<R3> vertex_wed(6);
  vertex_wed(0).Init(0, 0, 0);
  vertex_wed(1).Init(1, 0, 0);
  vertex_wed(2).Init(0, 1, 0);
  vertex_wed(3).Init(0, 0, 1);
  vertex_wed(4).Init(1, 0, 1);
  vertex_wed(5).Init(0, 1, 1);
  
  {
    cout << "Testing WedgeHcurlFirstFamily " << endl;
    WedgeHcurlFirstFamily tet; string name_elt("TETRAHEDRON_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_wed);
  }

  {
    cout << "Testing WedgeHcurlOptimalFirstFamily " << endl;
    WedgeHcurlOptimalFirstFamily tet; string name_elt("TETRAHEDRON_OPTIMAL_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_wed);
  }
  
  {
    cout << "Testing WedgeHcurlHpFirstFamily " << endl;
    WedgeHcurlHpFirstFamily tet; string name_elt("TETRAHEDRON_HP_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_wed);
  }

  {
    cout << "Testing WedgeHcurlOptimalHpFirstFamily " << endl;
    WedgeHcurlOptimalHpFirstFamily tet; string name_elt("TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_wed);
  }
  
  // one hexahedron  
  mesh3d.ReallocateVertices(8);
  mesh3d.Vertex(0).Init(1, 0, 0);
  mesh3d.Vertex(1).Init(0.5, 0.5*sqrt(3.0), 0);
  mesh3d.Vertex(2).Init(0.8, 0.2, sqrt(0.68));
  mesh3d.Vertex(3).Init(0.7, -0.3, sqrt(0.58));
  mesh3d.Vertex(4).Init(-0.2, 0.2, 0.1);
  mesh3d.Vertex(5).Init(-0.1, 0.72, -0.13);
  mesh3d.Vertex(6).Init(0.14, 0.81, 0.73);
  mesh3d.Vertex(7).Init(0.08, -0.17, 0.69);
  
  mesh3d.ReallocateBoundariesRef(6);
  mesh3d.BoundaryRef(0).InitQuadrangular(0, 1, 2, 3, 2);
  mesh3d.BoundaryRef(1).InitQuadrangular(0, 1, 5, 4, 1);
  mesh3d.BoundaryRef(2).InitQuadrangular(0, 3, 7, 4, 1);
  mesh3d.BoundaryRef(3).InitQuadrangular(1, 2, 6, 5, 1);
  mesh3d.BoundaryRef(4).InitQuadrangular(2, 3, 7, 6, 1);
  mesh3d.BoundaryRef(5).InitQuadrangular(4, 5, 6, 7, 1);
  
  mesh3d.ReallocateElements(1);
  mesh3d.Element(0).InitHexahedral(0, 1, 2, 3, 4, 5, 6, 7, 1);
  
  mesh3d.ReorientElements();
  mesh3d.FindConnectivity();
  mesh3d.ProjectPointsOnCurves();

  Vector<R3> vertex_hex(8);
  vertex_hex(0).Init(0, 0, 0);
  vertex_hex(1).Init(1, 0, 0);
  vertex_hex(2).Init(1, 1, 0);
  vertex_hex(3).Init(0, 1, 0);
  vertex_hex(4).Init(0, 0, 1);
  vertex_hex(5).Init(1, 0, 1);
  vertex_hex(6).Init(1, 1, 1);
  vertex_hex(7).Init(0, 1, 1);
  
  {
    cout << "Testing HexahedronHcurlFirstFamily " << endl;
    HexahedronHcurlFirstFamily tet; string name_elt("TETRAHEDRON_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_hex);
  }

  {
    cout << "Testing HexahedronHcurlOptimalFirstFamily " << endl;
    HexahedronHcurlOptimalFirstFamily tet; string name_elt("TETRAHEDRON_OPTIMAL_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_hex);
  }
  
  {  
    cout << "Testing HexahedronHcurlHpFirstFamily " << endl;
    HexahedronHcurlHpFirstFamily tet; string name_elt("TETRAHEDRON_HP_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_hex, true);
  } 
  
  {  
    cout << "Testing HexahedronHcurlOptimalHpFirstFamily " << endl;
    HexahedronHcurlOptimalHpFirstFamily tet; string name_elt("TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_hex, true);
  }
  
  {
    cout << "Testing HexahedronHcurlLobatto " << endl;
    HexahedronHcurlLobatto tet; string name_elt("HEXAHEDRON_HCURL_LOBATTO");
    CheckHcurlFiniteElement(tet, name_elt, mesh3d, vertex_hex);
  }

  cout << "All tests passed successfully" << endl;

  return FinalizeMontjoie();
}
