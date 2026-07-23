#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_ORTHO_DG
#define MONTJOIE_WITH_LEGENDRE_DG
#define MONTJOIE_WITH_HP_H1

#include "Harmonic/MontjoieHarmonic.hxx"

#ifdef MONTJOIE_WITH_FAST_LIBRARY
// includes for fast compilation
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#endif

#define CHECK_ABORT abort()
//#define CHECK_ABORT cout<<"ERROR"<<endl

int r_over_quad = 0;

using namespace Montjoie;

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}

Real_wp FuncPowS(const R2& x, int m, int n)
{
  return pow(x(0)+0.5*x(1), m)*pow(0.2+0.9*x(0)+x(1), n);
}

Real_wp FuncPowS(const R3& x, int m, int n)
{
  return pow(x(0)+0.5*x(1)-0.3*x(2), m)*pow(0.2+0.9*x(0)+x(1)+0.8*x(2), n);
  //return pow(x(0), m);
}

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

namespace Montjoie
{
  template<class Dim>
  class MyOwnEquation : public GenericEquation<Real_wp>
  {
  public :
    typedef Dim Dimension;
    enum{nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_vec = 0,
	 nb_components_en = 1, nb_components_hn = 1, nb_unknowns_hdg=0};
  };

  template<class Dimension>
  class TestDim
  {};

  template<>
  class TestDim<Dimension2>
  {
  public:
    template<int t>
    static void GetLocalCoord(const ElementReference<Dimension2, t>& Fb,
                              int num_loc, Real_wp tloc, R2& ptA)
    {
      Fb.GetLocalCoordOnBoundary(num_loc, tloc, ptA);
    }

    template<int t>
    static void GetLocalCoord(const ElementReference<Dimension2, t>& Fb,
                              int num_loc, TinyVector<Real_wp, 1> tloc, R2& ptA)
    {
      Fb.GetLocalCoordOnBoundary(num_loc, tloc(0), ptA);
    }

    template<int t>
    static void GetLocalCoord(const ElementReference<Dimension2, t>& Fb,
                       int num_loc, R2 tloc, R2& ptA)
    {
    }
    
  };

  template<>
  class TestDim<Dimension3>
  {
  public:
    template<int t>
    static void GetLocalCoord(const ElementReference<Dimension3, t>& Fb,
                       int num_loc, Real_wp tloc, R3& ptA)
    {
    }

    template<int t>
    static void GetLocalCoord(const ElementReference<Dimension3, t>& Fb,
                              int num_loc, R2 tloc, R3& ptA)
    {
      Fb.GetLocalCoordOnBoundary(num_loc, tloc, ptA);
    }
    
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
    virtual bool IsVaryingMedia(int ref) const { return false; }
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

#ifdef MONTJOIE_WITH_THREE_DIM
void CreateRegularMesh(Mesh<Dimension3>& mesh, const VectR3& s, int r, int rgeom)
{
  int type_mesh = Mesh<Dimension3>::HEXAHEDRAL_MESH;
  switch (s.GetM())
    {
    case 4 : type_mesh = Mesh<Dimension3>::TETRAHEDRAL_MESH; break;
    case 5 : type_mesh = Mesh<Dimension3>::PYRAMID_MESH; break;
    case 6 : type_mesh = Mesh<Dimension3>::WEDGE_MESH; break;
    }
  
  mesh.SetGeometryOrder(rgeom);
  TinyVector<int, 6> ref_bound; ref_bound.Fill(1);
  mesh.CreateRegularMesh(R3(-2,-2,-2), R3(2,2,2), TinyVector<int, 3>(3,3,3), 1,
			 ref_bound, type_mesh);
  
  mesh.ClearConnectivity();
  // random permutation of local vertex numbers
  //srand(time(NULL));
  //srand(0);
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
#endif

void CreateRegularMesh(Mesh<Dimension2>& mesh, const VectR2& s, int r, int rgeom)
{
  int type_mesh = Mesh<Dimension2>::QUADRILATERAL_MESH;
  if (s.GetM() == 3)
    type_mesh = Mesh<Dimension2>::TRIANGULAR_MESH;
  
  DISP(type_mesh);
  mesh.SetGeometryOrder(rgeom);
  TinyVector<int, 4> ref_bound; ref_bound.Fill(1);
  mesh.CreateRegularMesh(R2(-2,-2), R2(2,2), TinyVector<int, 2>(3,3), 1,
			 ref_bound, type_mesh);
  
  mesh.ClearConnectivity();
  // random permutation of local vertex numbers
  //srand(time(NULL));
  //srand(0);
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
	case 3 :
	  {
	    int perm = rand()%6;
	    switch (perm)
	      {
	      case 0 :
		num = numb; break;
	      case 1 :
		num(0) = numb(1); num(1) = numb(2); num(2) = numb(0); break;
	      case 2 :
		num(0) = numb(2); num(1) = numb(0); num(2) = numb(1); break;
	      case 3 :
		num(0) = numb(0); num(1) = numb(2); num(2) = numb(1); break;
	      case 4 :
		num(0) = numb(2); num(1) = numb(1); num(2) = numb(0); break;
	      case 5 :
		num(0) = numb(1); num(1) = numb(0); num(2) = numb(2); break;
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
      
      mesh.Element(i).Init(num, ref);
    }
  
}


void CheckRhSplit(const ElementReference<Dimension2, 1>& Fb,
		  VectReal_wp& Ones, VectReal_wp& MhX, VectReal_wp& MhY, 
		  VectReal_wp& MhZ, Matrix<R2>& GradientPhi, int m, int j2, const Real_wp& threshold)
{
  Fb.ApplyRhSplit(Ones, MhX, MhY);
  for (int i = 0; i < Fb.GetNbDof(); i++)
    {
      Real_wp val = GradientPhi(i, j2)(m);
      if ( ( (m == 0) && (abs(MhX(i)-val) > 100.0*threshold) && (abs(MhY(i)) > 100.0*threshold))
	   || ( (m == 1) && (abs(MhX(i)) > 100.0*threshold) && (abs(MhY(i)-val) > 100.0*threshold)) )
	{
	  cout << "ApplyRhSplit not correct " << endl;
	  CHECK_ABORT;
	}
    }
}


#ifdef MONTJOIE_WITH_THREE_DIM
void CheckRhSplit(const ElementReference<Dimension3, 1>& Fb,
		  VectReal_wp& Ones, VectReal_wp& MhX, VectReal_wp& MhY, 
		  VectReal_wp& MhZ, Matrix<R3>& GradientPhi, int m, int j2, const Real_wp& threshold)
{
  Fb.ApplyRhSplit(Ones, MhX, MhY, MhZ);
  for (int i = 0; i < Fb.GetNbDof(); i++)
    {
      Real_wp val = GradientPhi(i, j2)(m);
      if ( ( (m == 0) && (abs(MhX(i)-val) > 100.0*threshold) && (abs(MhY(i)) > 100.0*threshold) && (abs(MhZ(i)) > 100.0*threshold))
	   || ( (m == 1) && (abs(MhX(i)) > 100.0*threshold) && (abs(MhY(i)-val) > 100.0*threshold) && (abs(MhZ(i)) > 100.0*threshold)) 
	   || ( (m == 2) && (abs(MhX(i)) > 100.0*threshold) && (abs(MhY(i)) > 100.0*threshold) && (abs(MhZ(i)-val) > 100.0*threshold)) )
	{
	  cout << "ApplyRhSplit not correct " << endl;
	  CHECK_ABORT;
	}
    }
}
#endif

void CheckConstantRhSplit(const ElementReference<Dimension2, 1>& Fb,
			  VectReal_wp& Ones, VectReal_wp& MhX, VectReal_wp& MhY,
			  VectReal_wp& Mhz, Matrix<Real_wp>& grad_matrix, Matrix<R2>& GradientPhi,
			  int i, int m, const Real_wp& threshold)
{
  Fb.ApplyConstantRhSplit(Ones, MhX, MhY);
  for (int k = 0; k < Fb.GetNbDof(); k++)
    {
      Real_wp val = grad_matrix(k, 2*i+m);
      if ( ( (m == 0) && ((abs(MhX(k)-val) > 100.0*threshold) || (abs(MhY(k)) > 100.0*threshold)))
	   || ( (m == 1) && ((abs(MhX(k)) > 100.0*threshold) || (abs(MhY(k)-val) > 100.0*threshold))) )
	{
	  cout << "ApplyRhSplit not correct " << endl;                        
	  CHECK_ABORT;
	}
    }
}

#ifdef MONTJOIE_WITH_THREE_DIM
void CheckConstantRhSplit(const ElementReference<Dimension3, 1>& Fb,
			  VectReal_wp& Ones, VectReal_wp& MhX, VectReal_wp& MhY,
			  VectReal_wp& MhZ, Matrix<Real_wp>& grad_matrix, Matrix<R3>& GradientPhi,
			  int i, int m, const Real_wp& threshold)
{
  Fb.ApplyConstantRhSplit(Ones, MhX, MhY, MhZ);
  for (int k = 0; k < Fb.GetNbDof(); k++)
    {
      Real_wp val = grad_matrix(k, 3*i+m);
      if ( ( (m == 0) && ((abs(MhX(k)-val) > 100.0*threshold) || (abs(MhY(k)) > 100.0*threshold) || (abs(MhZ(k)) > 100.0*threshold) ))
	   || ( (m == 1) && ((abs(MhX(k)) > 100.0*threshold) || (abs(MhY(k)-val) > 100.0*threshold) || (abs(MhZ(k)) > 100.0*threshold) ))
	   || ( (m == 2) && ((abs(MhX(k)) > 100.0*threshold) || (abs(MhY(k)) > 100.0*threshold) || (abs(MhZ(k)-val) > 100.0*threshold) )) )
	{
	  cout << "ApplyConstantRhSplit not correct " << endl;
	  CHECK_ABORT;
	}
    }
}
#endif

void CheckSurfaceElement(const ElementReference<Dimension2, 1>& Fb)
{
  const ElementReference<Dimension1, 1>& edge = Fb.GetSurfaceFiniteElement(0);
  int r = Fb.GetOrder();
  if ((edge.GetOrder() != r) || (edge.GetGeometryOrder() != r) || 
      edge.GetNbPointsQuadratureInside() != Fb.GetNbQuadBoundary(0) || edge.GetNbDof() != (r+1))
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
      CHECK_ABORT;
    }  

  R2 pt_loc; VectR2 grad_phi;
  VectReal_wp phi_surf, phi, dphi_surf;
  for (int j = 0; j <= r; j++)
    {
      Fb.GetLocalCoordOnBoundary(0, edge.Points(j), pt_loc);
      edge.ComputeValuesPhiRef(edge.Points(j), phi_surf);
      Fb.ComputeValuesPhiRef(pt_loc, phi);

      edge.ComputeGradientPhiRef(edge.Points(j), dphi_surf);
      Fb.ComputeGradientPhiRef(pt_loc, grad_phi);
      for (int i = 0; i <= r; i++)
	{
	  int n = Fb.GetLocalNumber(0, i);
	  if (abs(phi_surf(i)- phi(n)) > 1e-12)
	    {
	      DISP(pt_loc); DISP(i); DISP(n); DISP(phi_surf(i)); DISP(phi(n));
	      cout << "Value of basis functions differ in GetSurfaceFiniteElement" << endl;
	      abort();
	    }

	  if (abs(dphi_surf(i) - grad_phi(n)(0)) > 1e-12)
	    {
	      cout << "Gradient of basis functions differ in GetSurfaceFiniteElement" << endl;
	      abort();
	    }
	}
    }
}

#ifdef MONTJOIE_WITH_THREE_DIM
void CheckSurfaceElement(const ElementReference<Dimension3, 1>& Fb)
{
  bool presence_tri = false, presence_quad = false;
  if (Fb.GetNbVertices() != 8)
    presence_tri = true;

  if (Fb.GetNbVertices() != 4)
    presence_quad = true;

  int r = Fb.GetOrder();
  int type_elt = Fb.GetHybridType();
  
  if (presence_tri)
    {
      const ElementReference<Dimension2, 1>& tri = Fb.GetTriangularSurfaceFiniteElement();

      for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	if (!MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(type_elt, num_loc))
	  {
	    
	    if ((tri.GetOrder() != r) || (tri.GetGeometryOrder() != r) || 
		tri.GetNbPointsQuadratureInside() != Fb.GetNbQuadBoundary(num_loc)
		|| tri.GetNbDof() != Fb.GetNbDofBoundary(num_loc))
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
	    
	    R3 pt_loc;
	    VectReal_wp phi_surf, phi;
	    for (int j = 0; j < tri.WeightsND().GetM(); j++)
	      {
		Fb.GetLocalCoordOnBoundary(num_loc, tri.PointsND(j), pt_loc);
		tri.ComputeValuesPhiRef(tri.PointsND(j), phi_surf);
		Fb.ModifySignProjectionSurface(phi_surf, num_loc);
		Fb.ComputeValuesPhiRef(pt_loc, phi);
		for (int i = 0; i < tri.GetNbDof(); i++)
		  {
		    int n = Fb.GetLocalNumber(num_loc, i);
		    if (abs(phi_surf(i)- phi(n)) > 1e-12)
		      {
			DISP(pt_loc); DISP(i); DISP(n); DISP(phi_surf(i)); DISP(phi(n));
			cout << "Value of basis functions differ in GetTriangularSurfaceFiniteElement" << endl;
			abort();
		      }
		  }
	      }
	  }
    }

  if (presence_quad)
    {
      const ElementReference<Dimension2, 1>& quad = Fb.GetQuadrangularSurfaceFiniteElement();
      
      for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	if (MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(type_elt, num_loc))
	  {
	    
	    if ((quad.GetOrder() != r) || (quad.GetGeometryOrder() != r) || 
		quad.GetNbPointsQuadratureInside() != Fb.GetNbQuadBoundary(num_loc)
		|| quad.GetNbDof() != Fb.GetNbDofBoundary(num_loc))
	      {
		cout << " GetQuadrangularSurfaceFiniteElement incorrect " << endl;
		CHECK_ABORT;
	      }
	    
	    /* if (!EqualVector(quad.PointsQuadInsideND(), Fb.Points2D_quad()))
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
	    
	    R3 pt_loc;
	    VectReal_wp phi_surf, phi;
	    for (int j = 0; j < quad.WeightsND().GetM(); j++)
	      {
		Fb.GetLocalCoordOnBoundary(num_loc, quad.PointsND(j), pt_loc);
		quad.ComputeValuesPhiRef(quad.PointsND(j), phi_surf);
		Fb.ModifySignProjectionSurface(phi_surf, num_loc);
		Fb.ComputeValuesPhiRef(pt_loc, phi);
		for (int i = 0; i < quad.GetNbDof(); i++)
		  {
		    int n = Fb.GetLocalNumber(num_loc, i);
		    if (abs(phi_surf(i)- phi(n)) > 1e-12)
		      {
			DISP(pt_loc); DISP(i); DISP(n); DISP(phi_surf(i)); DISP(phi(n));
			cout << "Value of basis functions differ in GetQuadrangularSurfaceFiniteElement" << endl;
			abort();
		      }
		  }
	      }
	  }
    }
}
#endif

template<class Dimension>
void CheckContinuousFiniteElement(ElementReference<Dimension, 1>& Fb,
				  const string& name_elt,
				  Mesh<Dimension>& mesh, typename Dimension::VectR_N& s,
                                  bool check_continuity = true, bool check_surf = true)
{
  typedef typename Dimension::R_N R_N;
  typedef typename Dimension::VectR_N VectR_N;
  typedef typename Dimension::MatrixN_N MatrixN_N;

  MeshNumbering<Dimension> mesh_num(mesh);
  // checking continuous finite elements (H1 or DG)
  Real_wp aire = 1.0;
  if (Dimension::dim_N == 2)
    {
      if (s.GetM() == 3)
        aire = 0.5;
    }
  else
    {
      switch (s.GetM())
        {
        case 4 : aire = 1.0/6; break;
        case 5 : aire = 4.0/3; break;
        case 6 : aire = 0.5; break;
        }
    }
  
  Real_wp threshold = 300.0*epsilon_machine;
  R_N::threshold = threshold;
  
  int num_loc = 0;
  for (int r = 1; r <= 4; r++)
    {      
      cout << "Checking order " << r << endl;
      int rgeom = min(r, 4);
      int rquad = r +  r_over_quad;
      Fb.ConstructFiniteElement(r, rgeom, rquad);
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
      
      if (Fb.GetQuadratureOrder() != rquad)
        {
          cout << "Problem when setting quadrature order of finite element " << endl;
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
      
      Vector<R_N> points_nodal = Fb.PointsNodalND();
      for (int k = 0; k < Fb.GetNbVertices(); k++)
        if (points_nodal(k) != s(k))
          {
            cout << "first nodal points have to be vertices " << endl;
            CHECK_ABORT;
          }
      
      R_N ptA, ptB; int offset = Fb.GetNbVertices();
      for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
        {
          if (Dimension::dim_N == 2)
            for (int k = 0; k < rgeom-1; k++)
              {
                Real_wp tloc;
                tloc = points1d(k+1);
                TestDim<Dimension>::GetLocalCoord(Fb, num_loc, tloc, ptA);
                if (ptA != points_nodal(offset+k))
                  {
                    cout << "nodal points on edges have to be Lobatto points " << endl;
                    CHECK_ABORT;
                  }
                
                if (Fb.GetNodalNumber(num_loc, k+1) != offset+k)
                  {
                    cout << "Problem for nodes numbering " << endl;
                    CHECK_ABORT;
                  }
              }
          
          offset += rgeom-1;
        }
      
      int nb_nodes_bounds = Fb.GetNbVertices()*rgeom;
      if (Dimension::dim_N == 3)
        nb_nodes_bounds = (Fb.GetNbVertices()-2)*rgeom*rgeom + 2;
      
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
	      typedef typename Dimension::DimensionBoundary DimensionBoundary;
              typename DimensionBoundary::R_N tloc = Fb.PointsQuadratureBoundary(k, num_loc);
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
          cout << "Sum of weights is not equal to the volume/area of element " << endl;
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
      VectReal_wp phi, phip, phim;
      for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
        {
          Fb.ComputeValuesPhiNodalRef(Fb.PointsNodalND(i), phi);
          for (int j = 0; j < Fb.GetNbPointsNodalElt(); j++)
            {
              Real_wp delta = 0;
              if (i == j)
                delta = 1.0;
              
              if (abs(delta - phi(j)) > threshold)
                {
                  cout << "nodal shape functions are not satisfying delta_ij " << endl;
                  CHECK_ABORT;
                }
            }
        }

      // checking gradient of nodal shape functions
      VectR_N grad_phi; Real_wp h = pow(epsilon_machine, 1.0/3);
      for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
        {
          Fb.ComputeGradientPhiNodalRef(Fb.PointsNodalND(i), grad_phi);
          Real_wp dphi_num; 
          for (int k = 0; k < Dimension::dim_N; k++)
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
          cout << "Problem for dofs on boundaries, they are greater than total number of dofs " <<endl;
          CHECK_ABORT;
        }

      // checking gradient of basis functions
      for (int i = 0; i < Fb.GetNbPointsDof(); i++)
        if (i != 4)
          {
            Real_wp val_max = 0;
            Fb.ComputeGradientPhiRef(Fb.PointsDofND(i), grad_phi);
            for (int j = 0; j < Fb.GetNbDof(); j++)
              val_max = max(val_max, Norm2(grad_phi(j)));
            
            Real_wp dphi_num; 
            for (int k = 0; k < Dimension::dim_N; k++)
              {
                ptA = Fb.PointsDofND(i);
                ptA(k) += h;
                Fb.ComputeValuesPhiRef(ptA, phip);
                
                ptA(k) -= 2.0*h;
                Fb.ComputeValuesPhiRef(ptA, phim);
                
                for (int j = 0; j < Fb.GetNbDof(); j++)
                  {
                    dphi_num = (phip(j) - phim(j))/(2.0*h);
                    if (abs(dphi_num - grad_phi(j)(k))/val_max > sqrt(threshold))
                      {
                        DISP(i); DISP(j); DISP(k); DISP(dphi_num); DISP(grad_phi(j)(k));
                        cout << "Gradient of basis functions not correct " <<endl;
                        CHECK_ABORT;
                      }
                  }
              }
          }
      
      // checking transformation Fi for linear element
      mesh.Element(0).UnsetCurved();
      VectR_N vert_elt; R_N pt_loc, ptC; MatrixN_N dfj, dfjm1;
      mesh.GetVerticesElement(0, vert_elt);
      SetPoints<Dimension> PointsElem;
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
	  for (int k = 0; k < Dimension::dim_N; k++)
	    {
	      pt_loc = s(i); pt_loc(k) += h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptC, mesh, 0);
	      
	      pt_loc = s(i); pt_loc(k) -= h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptB, mesh, 0);
	      
	      for (int m = 0; m < Dimension::dim_N; m++)
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
	  FjInverseProblem<Dimension> inverseFj(mesh, 0);
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
      
      
      SetMatrices<Dimension> MatricesElem;
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
	  for (int k = 0; k < Dimension::dim_N; k++)
	    {
	      pt_loc = s(i); pt_loc(k) += h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptC, mesh, 0);
	      
	      pt_loc = s(i); pt_loc(k) -= h;
	      Fb.Fj(vert_elt, PointsElem, pt_loc, ptB, mesh, 0);
	      
	      for (int m = 0; m < Dimension::dim_N; m++)
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
	  FjInverseProblem<Dimension> inverseFj(mesh, 0);
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
      
      num_loc = 1; Real_wp dsj; R_N normale_fj;
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
      if ((!Fb.DiscontinuousElement()) && (check_continuity))
	{
	  EllipticProblem<MyOwnEquation<Dimension> > var;
          MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
	  
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
      
      cout << "Surface element ok" << endl;
      // checking mass matrix
      Matrix<Real_wp> ValuePhi(Fb.GetNbDof(), Fb.GetNbPointsQuadrature());
      Real_wp max_cij = 0;
      for (int k = 0; k < Fb.GetNbPointsQuadrature(); k++)
        {
          Fb.ComputeValuesPhiRef(Fb.PointsND(k), phi);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              ValuePhi(i, k) = phi(i);
              max_cij = max(ValuePhi(i, k), max_cij);
            }
        }
      
      Matrix<Real_wp, Symmetric, RowSymPacked> Mh(Fb.GetNbDof(), Fb.GetNbDof()), InvMh; Mh.Fill(0);
      Real_wp max_mij = 0;
      for (int i = 0; i < Fb.GetNbDof(); i++)
        for (int j = i; j < Fb.GetNbDof(); j++)
          {
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              Mh(i, j) += Fb.WeightsND(k)*ValuePhi(i, k)*ValuePhi(j, k);
                        
            max_mij = max(max_mij, abs(Mh(i, j)));
          }

      Matrix<Real_wp, Symmetric, RowSymPacked> mass_matrix = Mh;
          
      InvMh = Mh;
      GetInverse(InvMh);
      Real_wp max_invMij = 0;
      for (int i = 0; i < Fb.GetNbDof(); i++)
        for (int j = i; j < Fb.GetNbDof(); j++)
          max_invMij = max(max_invMij, abs(InvMh(i, j)));
      
      VectReal_wp Ones(Fb.GetNbDof()), MhOnes(Fb.GetNbDof());
      Ones.Fill(0); MhOnes.Fill(0);
      for (int i = 0; i < Fb.GetNbDof(); i++)
        {
          Ones(i) = 1.0; MhOnes = Ones;
          Fb.MltMassMatrix(MhOnes);
          for (int j = 0; j < Fb.GetNbDof(); j++)
            if (abs(MhOnes(j) - Mh(j, i))/max_mij > threshold)
              {
                cout << "MltMassMatrix not correct " << endl;
                CHECK_ABORT;
              }
          
	  MhOnes = Ones;
          Fb.SolveCholesky(SeldonNoTrans, MhOnes);
          Fb.SolveCholesky(SeldonTrans, MhOnes);
          for (int j = 0; j < Fb.GetNbDof(); j++)
            if (abs(InvMh(j, i) - MhOnes(j))/max_invMij > 100.0*threshold)
              {
                cout << "SolveCholesky not correct " << endl;
                CHECK_ABORT;
              }
          
          MhOnes = Ones;
          Fb.SolveMassMatrix(MhOnes);
          for (int j = 0; j < Fb.GetNbDof(); j++)
            if (abs(InvMh(j, i) - MhOnes(j))/max_invMij > 100.0*threshold)
              {
                cout << "SolveMassMatrix not correct " << endl;
                CHECK_ABORT;
              }
	  
          Ones(i) = 0.0;
        }
      
      VectReal_wp jacob(Fb.GetNbPointsQuadratureInside());
      jacob.FillRand(); Mlt(1e-9, jacob); Mh.Fill(0);
      max_mij = 0;
      for (int i = 0; i < Fb.GetNbDof(); i++)
        for (int j = i; j < Fb.GetNbDof(); j++)
          {
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              Mh(i, j) += Fb.WeightsND(k)*jacob(k)*ValuePhi(i, k)*ValuePhi(j, k);
            
            max_mij = max(max_mij, abs(Mh(i, j)));
          }

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
      
      if (Fb.LumpedMassMatrix() != Fb.DiagonalMassMatrix())
        {
          cout << "For H1 elements, mass lumping => diagonal mass matrix " << endl;
          CHECK_ABORT;
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
      
      max_mij = 0; Mh.Fill(0);
      for (int i = 0; i < Fb.GetNbDof(); i++)
        for (int j = i; j < Fb.GetNbDof(); j++)
          {
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              Mh(i, j) += Fb.WeightsND(k)*jacob(k)*ValuePhi(i, k)*ValuePhi(j, k);
            
            max_mij = max(max_mij, abs(Mh(i, j)));
          }
      
      if (Fb.LinearSparseMassMatrix())
        {
          Fb.ComputeCoefJacobian(vert_elt, jacob);
          Fb.ComputeMassMatrix(InvMh, jacob);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            for (int j = 0; j < Fb.GetNbDof(); j++)
              if (abs(Mh(i, j) - InvMh(i, j))/ max_mij > threshold)
                {
                  DISP(i); DISP(j);
                  DISP(Mh(i, j)); DISP(InvMh(i, j));
                  cout << "ComputeMassMatrix not correct " << endl;
                  CHECK_ABORT;
                }
        }
      
      // checking ApplyCh/ComputeIntegralRef (same methods)
      Ones.Reallocate(Fb.GetNbPointsQuadratureInside()); Ones.Fill(0);
      for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
        {
          Ones(j) = 1.0;
          Fb.ComputeIntegralRef(Ones, MhOnes);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(MhOnes(i) - ValuePhi(i, j))/max_cij > threshold)
              {
                cout << "ComputeIntegralRef not correct " << endl;
                DISP(i); DISP(j); DISP(MhOnes(i)); 
                DISP(ValuePhi(i, j));
                CHECK_ABORT;
              }
          
          Fb.ApplyCh(Ones, MhOnes);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(MhOnes(i) - ValuePhi(i, j))/max_cij > threshold)
              {
                cout << "ApplyCh not correct " << endl;
                DISP(MhOnes(i)); DISP(ValuePhi(i, j));
                CHECK_ABORT;
              }
          
          Ones(j) = 0.0;
        }
      
      // and ApplyChTranspose
      MhOnes.Fill(0);
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          MhOnes(j) = 1.0;
          Fb.ApplyChTranspose(MhOnes, Ones);
          for (int i = 0; i < Fb.GetNbPointsQuadratureInside(); i++)
            if (abs(Ones(i) - ValuePhi(j, i))/max_cij > threshold)
              {
                cout << "ApplyChTranspose not correct " << endl;
                DISP(i); DISP(j); DISP(Ones(i)); DISP(ValuePhi(j, i));
                CHECK_ABORT;
              }
          
          MhOnes(j) = 0.0;
        }
      
      // checking ProjectQuadratureToDofRef
      VectReal_wp ProjOnes(Fb.GetNbDof());
      ProjOnes.Fill(0); Ones.Fill(0);
      for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
        {
          Ones(j) = 1.0;
          Fb.ProjectQuadratureToDofRef(Ones, MhOnes);
          Ones(j) *= Fb.WeightsND(j);
          Fb.ComputeIntegralRef(Ones, ProjOnes);
          Fb.SolveMassMatrix(ProjOnes);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(MhOnes(i) - ProjOnes(i)) > threshold)
              {
                cout << "ProjectQuadratureToDofRef not correct " << endl;
                CHECK_ABORT;
              }
          
          Ones(j) = 0.0;
        }
      
      // checking ComputeProjectionDofRef
      VectReal_wp feval(Fb.GetNbPointsDofInside()); feval.Fill(0);
      if (Dimension::dim_N == 2)
	{
	  for (int m = 0; m <= r; m++)
	    for (int n = 0; n <= r-m; n++)
	      {
		for (int j = 0; j < Fb.GetNbPointsDofInside(); j++)
		  {
		    Real_wp x = Fb.PointsDofND(j)(0);
		    Real_wp y = Fb.PointsDofND(j)(1);
		    feval(j) = pow(x, m)*pow(y, n);
		  }
		
		Fb.ComputeProjectionDofRef(feval, ProjOnes);
                Fb.ApplyChTranspose(ProjOnes, Ones);
		for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
		  {
		    Real_wp x = Fb.PointsND(j)(0);
		    Real_wp y = Fb.PointsND(j)(1);
		    Real_wp val_exact = pow(x, m)*pow(y, n);
		    if (abs(val_exact - Ones(j)) > 10*r*r*threshold)
		      {
			cout << "ComputeProjectionDofRef not correct " << endl;
			DISP(m); DISP(n); DISP(j); DISP(val_exact); DISP(Ones(j));
			CHECK_ABORT;
		      }
		  }
	      }
	}
      else
	{
	  for (int m = 0; m <= r; m++)
	    for (int n = 0; n <= r-m; n++)
	      for (int p = 0; p <= r-m-n; p++)
		{
		  for (int j = 0; j < Fb.GetNbPointsDofInside(); j++)
		    {
		      Real_wp x = Fb.PointsDofND(j)(0);
		      Real_wp y = Fb.PointsDofND(j)(1);
		      Real_wp z = Fb.PointsDofND(j)(2);
		      feval(j) = pow(x, m)*pow(y, n)*pow(z, p);
		    }
		  
		  Fb.ComputeProjectionDofRef(feval, ProjOnes);
		  Fb.ApplyChTranspose(ProjOnes, Ones);
		  for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
		    {
		      Real_wp x = Fb.PointsND(j)(0);
		      Real_wp y = Fb.PointsND(j)(1);
		      Real_wp z = Fb.PointsND(j)(2);
		      Real_wp val_exact = pow(x, m)*pow(y, n)*pow(z, p);
		      if (abs(val_exact - Ones(j)) > 10*r*r*r*threshold)
			{
			  cout << "ComputeProjectionDofRef not correct " << endl;
			  DISP(m); DISP(n); DISP(p); DISP(j); DISP(val_exact); DISP(Ones(j));
			  CHECK_ABORT;
			}
		    }
		}
	}      
      
      // checking ComputeNodalValuesRef and ComputeNodalGradientRef
      Matrix<Real_wp> ValuePhiNodal(Fb.GetNbDof(), Fb.GetNbPointsNodalElt());
      Matrix<R_N> GradPhiNodal(Fb.GetNbPointsNodalElt(), Fb.GetNbPointsNodalElt());
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
      MhOnes.Reallocate(Fb.GetNbPointsNodalElt()); MhOnes.Fill(0);
      for (int j = 0; j < Fb.GetNbDof(); j++)
        {
          Ones(j) = 1.0;
          Fb.ComputeNodalValuesRef(Ones, MhOnes);
          Fb.ComputeNodalGradientRef(MhOnes, grad_phi);
          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            if (abs(ValuePhiNodal(j, i) - MhOnes(i)) > threshold)
              {
                cout << "ComputeNodalValuesRef not correct " << endl;
                DISP(MhOnes(i)); DISP(ValuePhiNodal(j, i));
                CHECK_ABORT;
              }

          for (int i = 0; i < Fb.GetNbPointsNodalElt(); i++)
            {
              R_N grad;
              for (int k = 0; k < Fb.GetNbPointsNodalElt(); k++)
                Add(MhOnes(k), GradPhiNodal(k, i), grad);
              
              if (grad != grad_phi(i))
                {
                  cout << "ComputeNodalGradientRef not correct " << endl;
                  CHECK_ABORT;
                }
            }
          
          Ones(j) = 0.0;
        }
      
      // checking GetValuePhiOnQuadraturePoint (and Gradient)
      Matrix<R_N> GradientPhi(Fb.GetNbDof(), Fb.GetNbPointsQuadrature());
      Real_wp max_rij = 0;
      for (int k = 0; k < Fb.GetNbPointsQuadrature(); k++)
        {
          Fb.ComputeGradientPhiRef(Fb.PointsND(k), grad_phi);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              GradientPhi(i, k) = grad_phi(i);
              max_rij = max(max_rij, Norm2(GradientPhi(i, k)));
            }
        }
      
      for (int k = 0; k < Fb.GetNbPointsQuadrature(); k++)
        {
          Fb.GetValuePhiOnQuadraturePoint(k, phi);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(phi(i) - ValuePhi(i, k))/max_cij > threshold)
              {
                cout << "GetValuePhiOnQuadraturePoint not correct " << endl;
                DISP(phi(i)); DISP(ValuePhi(i, k));
                CHECK_ABORT;
              }
          
          Fb.GetGradientPhiOnQuadraturePoint(k, grad_phi);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (grad_phi(i).Distance(GradientPhi(i, k))/max_rij > threshold)
              {
                cout << "GetGradientPhiOnQuadraturePoint not correct " << endl;
                DISP(grad_phi(i)); DISP(GradientPhi(i, k));
                CHECK_ABORT;
              }
        }
      
      phi.Reallocate(Fb.GetNbPointsQuadratureInside());
      grad_phi.Reallocate(Dimension::dim_N*phi.GetM());
      for (int i = 0; i < Fb.GetNbDof(); i++)
	{	  
	  Fb.GetValueSinglePhiQuadrature(i, phi);

          for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
            if (abs(phi(k) - ValuePhi(i, k))/max_cij > threshold)
              {
                cout << "GetValueSinglePhiQuadrature not correct " << endl;
                DISP(phi(k)); DISP(ValuePhi(i, k));
                CHECK_ABORT;
              }
	  
	  Fb.GetGradientSinglePhiQuadrature(i, phi, grad_phi);
	  
	  for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
	    {
	      if ( (abs(phi(k) - ValuePhi(i, k))/max_cij > threshold)
		   || (grad_phi(k).Distance(GradientPhi(i, k))/max_rij > threshold) )
		{
		  cout << "GetGradientSinglePhiQuadrature not correct " << endl;
		  DISP(phi(k)); DISP(ValuePhi(i, k));
		  DISP(grad_phi(k)); DISP(GradientPhi(i, k));
		  CHECK_ABORT;		  
		}
	    }		
	}
      
      // checking ComputeIntegralGradientRef / ApplyRh
      Ones.Reallocate(Dimension::dim_N*Fb.GetNbPointsQuadratureInside());
      MhOnes.Reallocate(Fb.GetNbDof()); MhOnes.Fill(0);
      VectReal_wp MhX(Fb.GetNbDof()), MhY(Fb.GetNbDof()), MhZ(Fb.GetNbDof());
      Ones.Fill(0);
      for (int j = 0; j < Dimension::dim_N*Fb.GetNbPointsQuadratureInside(); j++)
        {
          int j2 = j/Dimension::dim_N;
          int m = j%Dimension::dim_N;
          Ones.Fill(0); Ones(j) = 1.0;
          Fb.ComputeIntegralGradientRef(Ones, MhOnes);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(MhOnes(i) - GradientPhi(i, j2)(m))/max_rij > 100.0*threshold)
              {
                cout << "ComputeIntegralGradientRef not correct " << endl;
		DISP(i); DISP(j2); DISP(m);
		DISP(MhOnes(i)); DISP(GradientPhi(i, j2)(m));
                CHECK_ABORT;
              }

	  Ones.Fill(0); Ones(j) = 1.0;
          Fb.ApplyRh(Ones, MhOnes);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            if (abs(MhOnes(i) - GradientPhi(i, j2)(m))/max_rij > 100.0*threshold)
              {
                cout << "ApplyRh not correct " << endl;
                CHECK_ABORT;
              }

	  CheckRhSplit(Fb, Ones, MhX, MhY, MhZ, GradientPhi, m, j2, threshold);
          
	  Ones(j) = 0.0;
        }
      
      // checking ApplyRhTranspose
      MhOnes.Fill(0);
      for (int i = 0; i < Fb.GetNbDof(); i++)
        {
          MhOnes(i) = 1.0;
          
          Fb.ApplyRhTranspose(MhOnes, Ones);
          for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
            for (int m = 0; m < Dimension::dim_N; m++)
              if (abs(Ones(j*Dimension::dim_N+m) - GradientPhi(i, j)(m))/max_rij > 100.0*threshold)
                {
                  cout << "ApplyRhTranspose not correct " << endl;
                  CHECK_ABORT;
                }
          
          MhOnes(i) = 0.0;
        }
      
      // checking ComputeIntegralSurfaceRef/ApplySh
      ComputeGaussLegendre(points1d, weights1d, rquad);
      VectR2 points2d_tri, points2d_quad;
      VectReal_wp weights2d_tri, weights2d_quad;
      TriangleQuadrature::ConstructQuadrature(2*rquad, points2d_tri, weights2d_tri);
      Matrix<int> NumNodes2D_quad, coor;
      MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(rquad, NumNodes2D_quad, coor);
      weights2d_quad.Reallocate((rquad+1)*(rquad+1));
      points2d_quad.Reallocate((rquad+1)*(rquad+1));
      for (int i = 0; i <= rquad; i++)
        for (int j = 0; j <= rquad; j++)
          {
            points2d_quad(NumNodes2D_quad(i, j)).Init(points1d(i), points1d(j));
            weights2d_quad(NumNodes2D_quad(i, j)) = weights1d(i)*weights1d(j);
          }
      
      IVect num_point_quad;
      Matrix<bool> IsDofOnBoundary(Fb.GetNbBoundaries(), Fb.GetNbDof());
      IsDofOnBoundary.Fill(false);
      for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
        {
          int nb_pts_quad = Fb.GetNbQuadBoundary(num_loc);
          Ones.Reallocate(nb_pts_quad); Ones.Fill(0);
          for (int k = 0; k < nb_pts_quad; k++)
            {
              Ones(k) = 1.0; MhOnes.Fill(1.0);
              int k2 = Fb.GetQuadNumber(num_loc, k);
              num_point_quad.PushBack(k2);
              Fb.ApplySh(2.5, num_loc, Ones, MhOnes);
              for (int i = 0; i < Fb.GetNbDof(); i++)
                if (abs(MhOnes(i) - 2.5*ValuePhi(i, k2) - 1.0)/max_cij > threshold)
                  {
                    cout << "ApplySh not correct " << endl;
                    DISP(MhOnes(i) - 1.0);
                    DISP(2.5*ValuePhi(i, k2));
                    DISP(num_loc);
                    CHECK_ABORT;
                  }

              Fb.ComputeIntegralSurfaceRef(Ones, MhOnes, num_loc);
              for (int i = 0; i < Fb.GetNbDof(); i++)
                if (abs(MhOnes(i) - ValuePhi(i, k2))/max_cij > threshold)
                  {
                    cout << "ComputeIntegralSurfaceRef not correct " << endl;
                    CHECK_ABORT;
                  }
                            
              Ones(k) = 0.0;
            }
          
          // checking ApplyShTranspose
          Ones.Fill(0); MhOnes.Fill(0);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              MhOnes(i) = 1.0;
              Fb.ApplyShTranspose(num_loc, MhOnes, Ones);
              for (int k = 0; k < nb_pts_quad; k++)
                {
                  int k2 = Fb.GetQuadNumber(num_loc, k);
                  if (abs(Ones(k) - ValuePhi(i, k2))/max_cij > threshold)
                    {
                      cout << "ApplyShTranspose not correct " << endl;
                      DISP(Ones(k)); DISP(ValuePhi(i, k2));
                      DISP(num_loc);
                      CHECK_ABORT;
                    }
                }
              
              MhOnes(i) = 0.0;
            }
          
          Ones.Reallocate(nb_pts_quad*Dimension::dim_N);
          Ones.Fill(0);
          for (int k = 0; k < nb_pts_quad; k++)
            for (int m = 0; m < Dimension::dim_N; m++)
              {
                Ones(k*Dimension::dim_N+m) = 1.0;
                int k2 = Fb.GetQuadNumber(num_loc, k);
                Fb.ComputeIntegralSurfaceGradientRef(Ones, MhOnes, num_loc);
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  if (abs(MhOnes(i) - GradientPhi(i, k2)(m))/max_rij > threshold)
                    {
                      cout << "ComputeIntegralSurfaceGradientRef not correct " << endl;
                      DISP(MhOnes(i)); DISP(GradientPhi(i, k2)(m)); DISP(num_loc);
                      CHECK_ABORT;
                    }
                
                Ones(k*Dimension::dim_N+m) = 0.0;
              }
          
          Ones.Fill(0);
          for (int k = 0; k < nb_pts_quad; k++)
            for (int m = 0; m < Dimension::dim_N; m++)
              {
                Ones(k*Dimension::dim_N+m) = 1.0;
                int k2 = Fb.GetQuadNumber(num_loc, k);
                MhOnes.Fill(0);
                Fb.ApplyNablaSh(1.0, num_loc, Ones, MhOnes);
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  if (abs(MhOnes(i) - GradientPhi(i, k2)(m))/max_rij > 100.0*threshold)
                    {
                      cout << "ApplyNablaSh not correct " << endl;
                      DISP(MhOnes(i)); DISP(GradientPhi(i, k2)(m)); DISP(num_loc);
                      CHECK_ABORT;
                    }
                
                Ones(k*Dimension::dim_N+m) = 0.0;
              }

          MhOnes.Fill(0);
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              MhOnes(i) = 1.0;
              Fb.ApplyNablaShTranspose(num_loc, MhOnes, Ones);
              for (int k = 0; k < nb_pts_quad; k++)
                for (int m = 0; m < Dimension::dim_N; m++)
                  {
                    int k2 = Fb.GetQuadNumber(num_loc, k);
                    if (abs(Ones(k*Dimension::dim_N+m) - GradientPhi(i, k2)(m))/max_rij > 100.0*threshold)
                      {
                        cout << "ApplyNablaShTranspose not correct " << endl;
                        DISP(Ones(k*Dimension::dim_N+m)); DISP(GradientPhi(i, k2)(m)); DISP(num_loc);
                        CHECK_ABORT;
                      }
                  }
              
              MhOnes(i) = 0;
            }
          

          nb_pts_quad = points1d.GetM();
          VectR2 pts_integration; VectReal_wp poids_integration;
#ifdef MONTJOIE_WITH_THREE_DIM
          if (Dimension::dim_N == 3)
            {
              if (MeshNumbering<Dimension3>::
		  IsLocalFaceQuadrilateral(Fb.GetHybridType(), num_loc))
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
            }
#endif
              
          Ones.Reallocate(nb_pts_quad); Ones.Fill(0);
          for (int k = 0; k < Ones.GetM(); k++)
            {
              Ones(k) = 1.0; MhOnes.Fill(1.0);
              TinyVector<Real_wp, Dimension::dim_N-1> xi;
              
              if (Dimension::dim_N == 3)
                {
                  xi(0) = pts_integration(k)(0);
                  xi(1) = pts_integration(k)(1);
                }
              else
                {
                  xi(0) = points1d(k);
                }
              
              TestDim<Dimension>::GetLocalCoord(Fb, num_loc, xi, pt_loc);
              Fb.ComputeValuesPhiRef(pt_loc, phi);
              
	      Fb.ComputeGaussIntegralSurfaceRef(Ones, MhOnes, num_loc);
              for (int i = 0; i < Fb.GetNbDof(); i++)
                if (abs(MhOnes(i) - phi(i))/max_cij > threshold)
                  {
                    cout << "ComputeGaussIntegralSurfaceRef not correct " << endl;
                    CHECK_ABORT;
                  }
	      
              Ones(k) = 0.0;

	      // detecting surface dofs
	      for (int i = 0; i < Fb.GetNbDof(); i++)
		if (abs(phi(i)) > threshold)
		  {
		    IsDofOnBoundary(num_loc, i) = true;
		  }
            }
        }
    
      // checking dof boundaries
      if (!Fb.DiscontinuousElement())
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
      
      // checking Apply, ApplyShTranspose for variable order
      //if (Fb.DiscontinuousElement())
      if (false)
	for (int rv = 1; rv <= 8; rv++)
	  {	    
	    Mesh<Dimension> mesh_elt;
	    MeshNumbering<Dimension> mesh_num(mesh_elt);
	    mesh_elt.ReallocateVertices(2*Fb.GetNbVertices());
	    for (int i = 0; i < s.GetM(); i++)
	      {
		mesh_elt.Vertex(i) = s(i);
		mesh_elt.Vertex(i+s.GetM()) = s(i);
	      }
	    
	    mesh_elt.ReallocateElements(2);
	    IVect num(s.GetM()); num.Fill();
	    mesh_elt.Element(0).Init(num, 1);
	    for (int k = 0; k < num.GetM(); k++)
	      num(k) += num.GetM();
	    
	    mesh_elt.Element(1).Init(num, 1);
	    
	    mesh_elt.FindConnectivity();
	    mesh_num.SetOrderElement(0, r);
	    mesh_num.SetOrderElement(1, rv);
	    if (r != rv)
              for (int i = 0; i < mesh_elt.GetNbBoundary(); i++)
                mesh_num.SetOrderQuadrature(i, rv);
	    
	    Fb.SetVariableOrder(mesh_elt, mesh_num);
	    int type_quad = Fb.GetTypeIntegrationEdge();
	    int type_quad_tri = Fb.GetTypeIntegrationTriangle();
	    int type_quad_quad = Fb.GetTypeIntegrationQuadrangle();
	    if (Dimension::dim_N == 2)
	      {
		Globatto<Real_wp> gauss;
		gauss.ConstructQuadrature(rv, type_quad);
		points1d = gauss.Points();
		weights1d = gauss.Weights();
	      }
	    else
	      {
		Globatto<Real_wp> gauss;
		gauss.ConstructQuadrature(rv, type_quad_quad);
		points1d = gauss.Points();
		weights1d = gauss.Weights();
		TriangleQuadrature::ConstructQuadrature(2*rv, points2d_tri, weights2d_tri, type_quad_tri);
		MeshNumbering<Dimension2>::
		  ConstructQuadrilateralNumbering(rv, NumNodes2D_quad, coor);
		
		weights2d_quad.Reallocate((rv+1)*(rv+1));
		points2d_quad.Reallocate((rv+1)*(rv+1));
		for (int i = 0; i <= rv; i++)
		  for (int j = 0; j <= rv; j++)
		    {
		      points2d_quad(NumNodes2D_quad(i, j)).Init(points1d(i), points1d(j));
		      weights2d_quad(NumNodes2D_quad(i, j)) = weights1d(i)*weights1d(j);
		    }
		
	      }
	    
	    VectReal_wp MhQuad(Fb.GetNbPointsQuadratureInside());
	    for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	      {
		int nb_pts_quad = points1d.GetM();
		VectR2 pts_integration; VectReal_wp poids_integration;
#ifdef MONTJOIE_WITH_THREE_DIM
		if (Dimension::dim_N == 3)
		  {
		    if (MeshNumbering<Dimension3>::
			IsLocalFaceQuadrilateral(Fb.GetHybridType(), num_loc))
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
		  }
#endif
		
		Ones.Reallocate(nb_pts_quad); Ones.Fill(0);
                VectReal_wp Gones(nb_pts_quad*Dimension::dim_N); Gones.Fill(0);
		Matrix<Real_wp> ValPhiR(Fb.GetNbDof(), nb_pts_quad);
                Matrix<R_N> GradPhiR(Fb.GetNbDof(), nb_pts_quad);
		for (int k = 0; k < nb_pts_quad; k++)
		  {
		    TinyVector<Real_wp, Dimension::dim_N-1> xi;
		    
		    if (Dimension::dim_N == 3)
		      {
			xi(0) = pts_integration(k)(0);
			xi(1) = pts_integration(k)(1);
		      }
		    else
		      {
			xi(0) = points1d(k);
		      }
		    
		    TestDim<Dimension>::GetLocalCoord(Fb, num_loc, xi, pt_loc);
		    Fb.ComputeValuesPhiRef(pt_loc, phi);
		    Fb.ComputeGradientPhiRef(pt_loc, grad_phi);
		    
		    Ones.Fill(0); Ones(k) = 1.0; MhOnes.Fill(1.0); MhQuad.Fill(0.0);
		    if (Fb.UseQuadraturePointsForSh())
		      {
			Fb.ApplyShQuadrature(2.5, num_loc, Ones, MhQuad, rv);
			Fb.ApplyCh(MhQuad, MhOnes);
			for (int n = 0; n < MhOnes.GetM(); n++)
			  MhOnes(n) += 1.0;
		      }
		    else
		      Fb.ApplySh(2.5, num_loc, Ones, MhOnes, rv);
		    
		    for (int i = 0; i < Fb.GetNbDof(); i++)
		      {
			ValPhiR(i, k) = phi(i);
			if (abs(MhOnes(i) - 2.5*phi(i) - 1.0)/max_cij > threshold)
			  {
			    cout << "ApplySh (variable order) not correct " << endl;
			    DISP(MhOnes(i) - 1.0);
			    DISP(2.5*phi(i));
			    DISP(num_loc);
			    CHECK_ABORT;
			  }
		      }
                    
                    for (int m = 0; m < Dimension::dim_N; m++)
                      {
                        Gones(k*Dimension::dim_N+m) = 1.0;
                        MhOnes.Fill(1.0); MhQuad.Fill(0.0);
                        if (Fb.UseQuadraturePointsForSh())
                          {
                            Fb.ApplyNablaShQuadrature(2.5, num_loc, Gones, MhQuad, rv);
                            Fb.ApplyCh(MhQuad, MhOnes);
                            for (int n = 0; n < MhOnes.GetM(); n++)
                              MhOnes(n) += 1.0;
                          }
                        else
                          Fb.ApplyNablaSh(2.5, num_loc, Gones, MhOnes, rv);
                        
                        for (int i = 0; i < Fb.GetNbDof(); i++)
                          {
                            GradPhiR(i, k)(m) = grad_phi(i)(m);
                            if (abs(MhOnes(i) - 2.5*grad_phi(i)(m) - 1.0)/max_rij > 10.0*r*rv*threshold)
                              {
                                cout << "ApplyNablaSh (variable order) not correct " << endl;
                                DISP(MhOnes(i) - 1.0);
                                DISP(2.5*grad_phi(i)(m));
                                DISP(i); DISP(m); DISP(num_loc);
                                CHECK_ABORT;
                              }
                          }
                        
                        Gones(k*Dimension::dim_N+m) = 0.0;
                      }
                  }
                
		Ones.Fill(0); MhOnes.Fill(0);
		for (int i = 0; i < Fb.GetNbDof(); i++)
		  {
		    MhOnes(i) = 1.0;
		    if (Fb.UseQuadraturePointsForSh())
		      {
			Fb.ApplyChTranspose(MhOnes, MhQuad);
			Fb.ApplyShQuadratureTranspose(num_loc, MhQuad, Ones, rv);
		      }
		    else
		      Fb.ApplyShTranspose(num_loc, MhOnes, Ones, rv);
		    
		    for (int k = 0; k < nb_pts_quad; k++)
		      {
			if (abs(Ones(k) - ValPhiR(i, k))/max_cij > threshold)
			  {
			    cout << "ApplyShTranspose (variable order) not correct " << endl;
			    DISP(Ones(k)); DISP(ValPhiR(i, k));
			    DISP(num_loc);
			    CHECK_ABORT;
			  }
		      }

                    if (Fb.UseQuadraturePointsForSh())
		      {
			Fb.ApplyChTranspose(MhOnes, MhQuad);
			Fb.ApplyNablaShQuadratureTranspose(num_loc, MhQuad, Gones, rv);
		      }
		    else
		      Fb.ApplyNablaShTranspose(num_loc, MhOnes, Gones, rv);
		    
                    for (int m = 0; m < Dimension::dim_N; m++)
                      for (int k = 0; k < nb_pts_quad; k++)
                        {
                          if (abs(Gones(k*Dimension::dim_N+m) - GradPhiR(i, k)(m))/max_rij > 10.0*r*rv*threshold)
                            {
                              cout << "ApplyNablaShTranspose (variable order) not correct " << endl;
                              DISP(Gones(k*Dimension::dim_N+m)); DISP(GradPhiR(i, k)(m));
                              DISP(num_loc);
                              CHECK_ABORT;
                            }
                        }
                    
		    MhOnes(i) = 0.0;
		  }
	      }          
	    
	  }      
      
      // checking ComputeProjectionSurfaceDofRef
      if (!Fb.DiscontinuousElement())
	for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
          {
	    if ((Fb.GetNbPointsDofInside() == 0) || (Fb.GetNbPointsDofInside() < Fb.GetNbDof()))
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

	    const typename Dimension::DimensionBoundary::VectR_N& pts_face = Fb.PointsDofBoundary(num_loc);
	    if (pts_face.GetM() != Fb.GetNbPointsDofSurface(num_loc))
	      {
		cout << "PointsDofBoundary incorrect" << endl;
		abort();
	      }
	    
	    for (int i = 0; i < Fb.GetNbPointsDofSurface(num_loc); i++)
	      {
		typename Dimension::R_N pt_loc;
		Fb.GetLocalCoordOnBoundary(num_loc, pts_face(i), pt_loc);
		if (pt_loc != Fb.PointsDofND(Fb.GetPointDofNumber(num_loc, i)))
		  {
		    DISP(num_loc); DISP(i);
		    DISP(pt_loc); DISP(Fb.PointsDofND(Fb.GetPointDofNumber(num_loc, i)));
		    cout << "PointsDofBoundary incorrect" << endl;
		    abort();
		  }
	      }
	    
            for (int m = 0; m <= r; m++)
              for (int n = 0; n <= r-m; n++)
                {
                  Vector<Real_wp> feval_surf(Fb.GetNbPointsDofSurface(num_loc));
                  for (int i = 0; i < feval_surf.GetM(); i++)
		    feval_surf(i) = FuncPowS(Fb.PointsDofND(Fb.GetPointDofNumber(num_loc, i)), m, n);
                  
		  VectReal_wp contrib_face(Fb.GetNbDofBoundary(num_loc));
                  Fb.ComputeProjectionSurfaceDofRef(feval_surf, contrib_face, num_loc);
                  VectReal_wp feval_vol(Fb.GetNbPointsDofInside());
                  for (int i = 0; i < feval_vol.GetM(); i++)
                    feval_vol(i) = FuncPowS(Fb.PointsDofND(i), m, n);
                  
                  Fb.ComputeProjectionDofRef(feval_vol, ProjOnes);
                  for (int i = 0; i < Fb.GetNbDofBoundary(num_loc); i++)
                    if (abs(contrib_face(i) - ProjOnes(Fb.GetLocalNumber(num_loc, i))) > 20*r*r*r*threshold)
                      {			
                        cout << "ComputeProjectionSurfaceDofRef not correct " << endl;
			DISP(i); DISP(num_loc); DISP(m); DISP(n); 
                        DISP(contrib_face(i));
                        DISP(ProjOnes(Fb.GetLocalNumber(num_loc, i)));

			DISP(contrib_face);
			for (int i2 = 0; i2 < Fb.GetNbDofBoundary(num_loc); i2++)
			  {
			    DISP(i2);
			    DISP(ProjOnes(Fb.GetLocalNumber(num_loc, i2)));
			  }
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
          VectReal_wp Vh(nb_points_quadrature_inside);
          VectReal_wp AhUh(nb_dof_loc);
          VectReal_wp Uh(Dimension::dim_N*nb_points_quadrature_inside);
          VectReal_wp AhUh2(nb_dof_loc);
          for (int i = 0; i < Uh.GetM(); i++)
            {
              Uh.Fill(0); Uh(i) = 1.0;
              Fb.ApplyRh(Uh, AhUh2);
	      Uh.Fill(0); Uh(i) = 1.0;
              Fb.ApplyRhQuadrature(Uh, Vh);
              Fb.ApplyCh(Vh, AhUh);
              for (int j = 0; j < AhUh.GetM(); j++)
                if (abs(AhUh(j) - AhUh2(j))/max_rij > 10.0*threshold)
                  {
                    DISP(AhUh(j));
                    DISP(AhUh2(j));
                    cout << "ApplyRhQuadrature not correct " << endl;
                    CHECK_ABORT;
                  }
            }
                    
          AhUh.Reallocate(Dimension::dim_N*nb_points_quadrature_inside);
          Uh.Reallocate(nb_dof_loc); AhUh2.Reallocate(Dimension::dim_N*nb_points_quadrature_inside);
          for (int i = 0; i < Uh.GetM(); i++)
            {
              Uh.Fill(0); Uh(i) = 1.0;
              Fb.ApplyRhTranspose(Uh, AhUh);
              Fb.ApplyChTranspose(Uh, Vh);
              Fb.ApplyRhQuadratureTranspose(Vh, AhUh2);
              for (int j = 0; j < AhUh.GetM(); j++)
                if (abs(AhUh(j) - AhUh2(j))/max_rij > 10.0*threshold)
                  {
                    DISP(AhUh(j));
                    DISP(AhUh2(j));
                    cout << "ApplyRhQuadratureTranspose not correct " << endl;
                    CHECK_ABORT;
                  }
            }
          
        }
      
      if (Fb.UseQuadraturePointsForSh())
        {
          int nb_dof_loc = Fb.GetNbDof();
          int nb_points_quadrature_inside = Fb.GetNbPointsQuadratureInside();
          VectReal_wp Vh(nb_points_quadrature_inside);
          VectReal_wp AhUh, Uh, AhUh2;
          
          for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
            {
              int nb_points_quadrature_boundary = Fb.GetNbQuadBoundary(num_loc);
              Vh.Reallocate(nb_points_quadrature_inside);
              AhUh.Reallocate(nb_dof_loc);
              Uh.Reallocate(nb_points_quadrature_boundary);
              AhUh2.Reallocate(nb_dof_loc);
              for (int i = 0; i < Uh.GetM(); i++)
                {
                  Uh.Fill(0); Uh(i) = 1.0;
                  AhUh.Fill(0);
                  Fb.ApplySh(1.0, num_loc, Uh, AhUh);
                  AhUh2.Fill(0); Vh.Fill(0);
                  Fb.ApplyShQuadrature(1.0, num_loc, Uh, Vh);
                  Fb.ApplyCh(Vh, AhUh2);
                  for (int j = 0; j < AhUh.GetM(); j++)
                    if (abs(AhUh(j) - AhUh2(j))/max_cij > threshold)
                      {
                        DISP(AhUh(j));
                        DISP(AhUh2(j));
                        cout << "ApplyShQuadrature not correct " << endl;
                        CHECK_ABORT;
                      }
                }
              
              Uh.Reallocate(Dimension::dim_N*nb_points_quadrature_boundary);
              for (int i = 0; i < Uh.GetM(); i++)
                {
                  Uh.Fill(0); Uh(i) = 1.0;
                  AhUh.Fill(0);
                  Fb.ApplyNablaSh(1.0, num_loc, Uh, AhUh);
                  AhUh2.Fill(0); Vh.Fill(0);
                  Fb.ApplyNablaShQuadrature(1.0, num_loc, Uh, Vh);
                  Fb.ApplyCh(Vh, AhUh2);
                  for (int j = 0; j < AhUh.GetM(); j++)
                    if (abs(AhUh(j) - AhUh2(j))/max_rij > threshold)
                      {
                        DISP(num_loc);
                        DISP(i); DISP(AhUh(j)); DISP(AhUh2(j));
                        cout << "ApplyNablaShQuadrature not correct " << endl;
                        CHECK_ABORT;
                      }
                }
              
              
              AhUh.Reallocate(nb_points_quadrature_boundary);
              Uh.Reallocate(nb_dof_loc);
              AhUh2.Reallocate(nb_points_quadrature_boundary);
              for (int i = 0; i < Uh.GetM(); i++)
                {
                  Uh.Fill(0); Uh(i) = 1.0;
                  Fb.ApplyShTranspose(num_loc, Uh, AhUh);
                  AhUh2.Fill(0);
                  Fb.ApplyChTranspose(Uh, Vh);
                  Fb.ApplyShQuadratureTranspose(num_loc, Vh, AhUh2);
                  for (int j = 0; j < AhUh.GetM(); j++)
                    if (abs(AhUh(j) - AhUh2(j))/max_cij > threshold)
                      {
                        DISP(i); DISP(j); DISP(AhUh(j));
                        DISP(AhUh2(j)); DISP(num_loc);
                        cout << "ApplyShQuadratureTranspose not correct " << endl;
                        CHECK_ABORT;
                      }
                }

              AhUh.Reallocate(Dimension::dim_N*nb_points_quadrature_boundary);
              Uh.Reallocate(nb_dof_loc);
              AhUh2.Reallocate(Dimension::dim_N*nb_points_quadrature_boundary);
              for (int i = 0; i < Uh.GetM(); i++)
                {
                  Uh.Fill(0); Uh(i) = 1.0;
                  Fb.ApplyNablaShTranspose(num_loc, Uh, AhUh);
                  AhUh2.Fill(0);
                  Fb.ApplyChTranspose(Uh, Vh);
                  Fb.ApplyNablaShQuadratureTranspose(num_loc, Vh, AhUh2);
                  for (int j = 0; j < AhUh.GetM(); j++)
                    if (abs(AhUh(j) - AhUh2(j))/max_rij > threshold)
                      {
                        DISP(i); DISP(j); DISP(AhUh(j));
                        DISP(AhUh2(j)); DISP(num_loc);
                        cout << "ApplyNablaShQuadratureTranspose not correct " << endl;
                        CHECK_ABORT;
                      }
                }
            }
          
        } 
      
      // computing gradient and stiffness matrix
      int nb_dof_loc = Fb.GetNbDof();
      Matrix<Real_wp> grad_matrix(nb_dof_loc, Dimension::dim_N*nb_dof_loc);
      Matrix<Real_wp> stiff_matrix(Dimension::dim_N*nb_dof_loc, Dimension::dim_N*nb_dof_loc);
      grad_matrix.Fill(0); stiff_matrix.Fill(0);
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              for (int m = 0; m < Dimension::dim_N; m++)
                {
                  grad_matrix(i, Dimension::dim_N*j + m) += Fb.WeightsND(k)*GradientPhi(i, k)(m)*ValuePhi(j, k);
                  for (int n = 0; n < Dimension::dim_N; n++)
                    stiff_matrix(i*Dimension::dim_N + m, j*Dimension::dim_N+n) += Fb.WeightsND(k)*GradientPhi(i, k)(m)*GradientPhi(j, k)(n);
                }
          }
      
      // checking addition of elementary matrices
      int off_row = 12, off_col = 25;
      Matrix<Real_wp> AhOrig(off_row+nb_dof_loc, off_col+nb_dof_loc), AhNew(off_row + nb_dof_loc, off_col + nb_dof_loc);
      Matrix<Real_wp> AhSum(nb_dof_loc, nb_dof_loc);
      AhNew.Zero();
      
      Real_wp Amass; GetRandNumber(Amass);
      AhOrig.FillRand(); Mlt(1e-9*max_mij, AhOrig); AhNew = AhOrig;
      Fb.AddConstantMassMatrix(off_row, off_col, Amass, AhNew);
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
	    AhSum(i, j) = AhOrig(off_row+i, off_col+j) + Amass*mass_matrix(i, j);
            if (abs(AhNew(off_row+i, off_col + j) - AhOrig(off_row+i, off_col+j) - Amass*mass_matrix(i, j))/max_mij > threshold)
              {
                DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j)); DISP(Amass*mass_matrix(i, j));
                cout << "AddConstantMassMatrix not correct " << endl;
                CHECK_ABORT;
              }
          }
      
      AhNew = AhOrig;
      
      R_N Dgrad, Egrad; MatrixN_N Cgrad;
      Cgrad.FillRand(); Mlt(1e-9, Cgrad);
      Dgrad.FillRand(); Mlt(1e-9, Dgrad);
      Egrad.FillRand(); Mlt(1e-9, Egrad);
      Fb.AddConstantElemMatrix(off_row, off_col, Amass, Cgrad, Dgrad, Egrad,
			       TinyVector<bool, 4>(true, true, true, false), AhNew);
      
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val_exact = 0;
            for (int m = 0; m < Dimension::dim_N; m++)
              val_exact += Egrad(m)*grad_matrix(j, Dimension::dim_N*i+m);
	    
	    AhSum(i, j) += val_exact;
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_rij > threshold)
              {
                cout << "AddConstantElemMatrix not correct " << endl;
                CHECK_ABORT;
              }
          }
      
      AhNew = AhOrig;
      Fb.AddConstantElemMatrix(off_row, off_col, Amass, Cgrad, Dgrad, Egrad,
			       TinyVector<bool, 4>(true, true, false, true), AhNew);
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val_exact = 0;
            for (int m = 0; m < Dimension::dim_N; m++)
              val_exact += Dgrad(m)*grad_matrix(i, Dimension::dim_N*j+m);
            
	    AhSum(i, j) += val_exact;
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_rij > threshold)
              {
                cout << "AddConstantElemMatrix not correct " << endl;
                CHECK_ABORT;
              }
          }
      
      AhNew = AhOrig;      
      Fb.AddConstantElemMatrix(off_row, off_col, Amass, Cgrad, Dgrad, Egrad,
			       TinyVector<bool, 4>(true, false, true, true), AhNew);

      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val_exact = 0;
            for (int m = 0; m < Dimension::dim_N; m++)
              for (int n = 0; n < Dimension::dim_N; n++)
                val_exact += Cgrad(m, n)*stiff_matrix(i*Dimension::dim_N+m, j*Dimension::dim_N+n);
            
	    AhSum(i, j) += val_exact;
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_rij > threshold)
              {
                cout << "AddConstantElemMatrix not correct " << endl;
                DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j)); DISP(val_exact);
                CHECK_ABORT;
              }
          }

      AhNew = AhOrig;
      
      Fb.AddConstantElemMatrix(off_row, off_col, Amass, Cgrad, Dgrad, Egrad,
			       TinyVector<bool, 4>(false, false, false, false), AhNew);

      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
	    if ((abs(AhNew(off_row+i, off_col+j) - AhSum(i, j))/max_rij > threshold)  || isnan(abs(AhNew(i, j) - AhSum(i, j))))
	      {
		cout << "AddConstantElemMatrix not correct" << endl;
		abort();
	      }
	  }

      // variable matrices
      AhNew = AhOrig;
      VectReal_wp coefA;
      GenerateRandomVector(coefA, Fb.GetNbPointsQuadratureInside());
      Fb.AddVariableMassMatrix(off_row, off_col, coefA, AhNew);
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val = 0;
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              val += coefA(k)*ValuePhi(i, k)*ValuePhi(j, k);
            
	    AhSum(i, j) = AhOrig(off_row+i, off_col+j) + val;
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val)/max_mij > 100.0*threshold)
              {
                DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j)); DISP(val);
                cout << "AddVariableMassMatrix not correct " << endl;
                CHECK_ABORT;
              }
          }
      
      AhNew = AhOrig;
      Vector<R_N> coefD(Fb.GetNbPointsQuadratureInside());
      Vector<R_N> coefE(Fb.GetNbPointsQuadratureInside());
      Vector<MatrixN_N> coefC(Fb.GetNbPointsQuadratureInside());
      for (int i = 0; i < Fb.GetNbPointsQuadratureInside(); i++)
        {
          coefC(i).FillRand(); Mlt(1e-9, coefC(i));
          coefD(i).FillRand(); Mlt(1e-9, coefD(i));
          coefE(i).FillRand(); Mlt(1e-9, coefE(i));
        }
      
      Fb.AddVariableElemMatrix(off_row, off_col, coefA, coefC, coefD, coefE,
			       TinyVector<bool, 4>(true, true, true, false), AhNew);
      
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val_exact = 0;
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              val_exact += DotProd(coefE(k), GradientPhi(j, k))*ValuePhi(i, k);
	    
	    AhSum(i, j) += val_exact;
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_rij > 100.0*threshold)
              {
                DISP(AhNew(off_row+i, off_col+j)); DISP(AhOrig(off_row+i, off_col+j) + val_exact);
                DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact);
                cout << "AddVariableElemMatrix not correct " << endl;
                CHECK_ABORT;
              }
          }
      
      AhNew = AhOrig; R_N vec_u;
      Fb.AddVariableElemMatrix(off_row, off_col, coefA, coefC, coefD, coefE,
			       TinyVector<bool, 4>(true, true, false, true), AhNew);

      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val_exact = 0;
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              val_exact += DotProd(coefD(k), GradientPhi(i, k))*ValuePhi(j, k);
            
	    AhSum(i, j) += val_exact;
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_rij > 100.0*threshold)
              {
                DISP(AhNew(off_row+i, off_col+j)); DISP(AhOrig(off_row+i, off_col+j) + val_exact);
                DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact);
                cout << "AddVariableElemMatrix not correct " << endl;
                CHECK_ABORT;
              }
          }
      
      AhNew = AhOrig;      
      Fb.AddVariableElemMatrix(off_row, off_col, coefA, coefC, coefD, coefE,
			       TinyVector<bool, 4>(true, false, true, true), AhNew);

      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
            Real_wp val_exact = 0;
            for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
              {
                Mlt(coefC(k), GradientPhi(j, k), vec_u);
                val_exact += DotProd(vec_u, GradientPhi(i, k));
              }
            
	    AhSum(i, j) += val_exact;
            if (abs(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact)/max_rij > 100.0*threshold)
              {
                DISP(AhNew(off_row+i, off_col+j)); DISP(AhOrig(off_row+i, off_col+j) + val_exact);
                DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j) - val_exact);
                cout << "AddVariableElemMatrix not correct " << endl;
                DISP(AhNew(off_row+i, off_col+j) - AhOrig(off_row+i, off_col+j)); DISP(val_exact);
                CHECK_ABORT;
              }
          }
      
      AhNew = AhOrig;
      Fb.AddVariableElemMatrix(off_row, off_col, coefA, coefC, coefD, coefE,
			       TinyVector<bool, 4>(false, false, false, false), AhNew);
      
      for (int i = 0; i < nb_dof_loc; i++)
        for (int j = 0; j < nb_dof_loc; j++)
          {
	    if ((abs(AhNew(off_row+i, off_col+j) - AhSum(i, j))/max_rij > 100.0*r*threshold)  || isnan(abs(AhNew(off_row+i, off_col+j) - AhSum(i, j))))
	      {
                DISP(AhNew(off_row+i, off_col+j)); DISP(AhSum(i, j));
                DISP(AhNew(off_row+i, off_col+j) - AhSum(i, j));
		cout << "AddVariableElemMatrix not correct" << endl;
		abort();
	      }
	  }
      
      // checking ApplyConstantRh/ApplyConstantRhTranspose
      Ones.Reallocate(Dimension::dim_N*nb_dof_loc);
      MhOnes.Reallocate(nb_dof_loc);
      Ones.Fill(0); MhOnes.Fill(0);
      for (int i = 0; i < nb_dof_loc; i++)
        for (int m = 0; m < Dimension::dim_N; m++)
          {
            Ones(Dimension::dim_N*i+m) = 1.0;
            Fb.ApplyConstantRh(Ones, MhOnes);
            for (int k = 0; k < nb_dof_loc; k++)
              if (abs(MhOnes(k) - grad_matrix(k, Dimension::dim_N*i+m))/max_rij > threshold)
                {
                  cout << "ApplyConstantRh not correct " << endl;
                  DISP(MhOnes(k)); DISP(grad_matrix(k, Dimension::dim_N*i+m));
                  CHECK_ABORT;
                }

	    CheckConstantRhSplit(Fb, Ones, MhX, MhY, MhZ, grad_matrix, GradientPhi,
				 i, m, threshold);
	                
            Ones(Dimension::dim_N*i+m) = 0.0;
          }

      Ones.Fill(0); MhOnes.Fill(0);
      for (int k = 0; k < nb_dof_loc; k++)
        {
          MhOnes(k) = 1.0;
          Fb.ApplyConstantRhTranspose(MhOnes, Ones);
          for (int i = 0; i < nb_dof_loc; i++)
            for (int m = 0; m < Dimension::dim_N; m++)
              {
                if (abs(Ones(i*Dimension::dim_N+m) - grad_matrix(k, Dimension::dim_N*i+m))/max_rij > threshold)
                  {
                    cout << "ApplyConstantRhTranspose not correct " << endl;
                    DISP(Ones(i*Dimension::dim_N+m)); DISP(grad_matrix(k, Dimension::dim_N*i+m));
                    CHECK_ABORT;
                  }
                
              }
          
          MhOnes(k) = 0.0;
        }
      
      // checking ValueBoundaryRef, and GradientBoundaryRef
      Ones.Reallocate(Fb.GetNbDof()); Ones.Fill(0);
      MhOnes.Reallocate(Fb.GetNbPointsNodalElt()); MhOnes.Fill(0);
      grad_phi.Reallocate(Fb.GetNbPointsNodalElt());
      for (num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
        {
          VectReal_wp feval_surf(Fb.GetNbNodalBoundary(num_loc));
          feval.Reallocate(Dimension::dim_N*Fb.GetNbNodalBoundary(num_loc));
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              Ones(i) = 1.0;
              Fb.ComputeNodalValuesRef(Ones, MhOnes);
              Fb.ComputeNodalGradientRef(MhOnes, grad_phi);
              Fb.ComputeValueBoundaryRef(Ones, feval_surf, num_loc);
              Fb.ComputeGradientBoundaryRef(Ones, feval, num_loc);
              
              for (int k = 0; k < feval_surf.GetM(); k++)
                {
                  int j = Fb.GetNodalNumber(num_loc, k);
                  if (abs(feval_surf(k) - MhOnes(j)) > threshold)
                    {
                      cout << "ComputeValueBoundaryRef not correct " << endl;
                      DISP(feval_surf(k)); DISP(MhOnes(j)); DISP(num_loc);
                      CHECK_ABORT;
                    }
                  
                  if ((r == rgeom) && (j != 4))
                    for (int m = 0; m < Dimension::dim_N; m++)
                      if (abs(feval(Dimension::dim_N*k+m) - grad_phi(j)(m))/max_rij > threshold)
                        {
                          cout << "ComputeGradientBoundaryRef not correct " << endl;
                          DISP(feval(Dimension::dim_N*k+m)); DISP(grad_phi(j)(m)); DISP(num_loc);
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

  if (argc > 1)
    r_over_quad = atoi(argv[1]);
  
  // 2-D mesh with only one triangle
  Mesh<Dimension2> mesh2d;
  mesh2d.SetGeometryOrder(4);
  
  mesh2d.ReallocateVertices(3);
  mesh2d.Vertex(0).Init(-0.2, 0.2);
  mesh2d.Vertex(1).Init(1, 0);
  mesh2d.Vertex(2).Init(0.5, 0.5*sqrt(Real_wp(3)));
  
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
    cout << "Testing TriangleClassical " << endl;
    TriangleClassical tri;  string name_elt("TRIANGLE_CLASSICAL");
    CheckContinuousFiniteElement(tri, name_elt, mesh2d, vertex_tri);
  }
  
  {
    cout << "Testing TriangleLobatto " << endl;
    TriangleLobatto tri;  string name_elt("TRIANGLE_LOBATTO");
    CheckContinuousFiniteElement(tri, name_elt, mesh2d, vertex_tri);
  }
  
  {
    cout << "Testing TriangleHierarchic " << endl;
    TriangleHierarchic tri; string name_elt("TRIANGLE_HIERARCHIC");
    CheckContinuousFiniteElement(tri, name_elt, mesh2d, vertex_tri);
  }

  if (r_over_quad == 0)
  {
    cout << "Testing TriangleDgOrtho " << endl;
    TriangleDgOrtho tri; string name_elt("TRIANGLE_DG_ORTHO");
    CheckContinuousFiniteElement(tri, name_elt, mesh2d, vertex_tri, false, false);
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
    cout << endl << endl << "Testing QuadrangleGauss " << endl;
    QuadrangleGauss quad;  string name_elt("TRIANGLE_CLASSICAL");
    CheckContinuousFiniteElement(quad, name_elt, mesh2d, vertex_quad);
  }
  
  {
    cout << endl << endl << "Testing QuadrangleLobatto " << endl;
    QuadrangleLobatto quad; string name_elt("TRIANGLE_LOBATTO");
    CheckContinuousFiniteElement(quad, name_elt, mesh2d, vertex_quad);
  }

  if (r_over_quad == 0)
  {
    cout << endl << endl << "Testing QuadrangleRadau " << endl;
    QuadrangleRadau quad; string name_elt("TRIANGLE_RADAU");
    // quad.SetRadauInterpolationPoints();    
    CheckContinuousFiniteElement(quad, name_elt, mesh2d, vertex_quad, false);
    //CheckContinuousFiniteElement(quad, name_elt, mesh2d, vertex_quad);
  }
 
  {
    cout << endl << endl << "Testing QuadrangleHierarchic " << endl;
    QuadrangleHierarchic quad; string name_elt("TRIANGLE_HIERARCHIC");
    CheckContinuousFiniteElement(quad, name_elt, mesh2d, vertex_quad);
  }
    
  {
    cout << endl << endl << "Testing QuadrangleDgGauss " << endl;
    QuadrangleDgGauss quad; string name_elt("TRIANGLE_DG_CLASSICAL");
    CheckContinuousFiniteElement(quad, name_elt, mesh2d, vertex_quad, false, false);
  }

#ifdef MONTJOIE_WITH_THREE_DIM
  // 3-D mesh with only one tetrahedron
  Mesh<Dimension3> mesh3d;
  mesh3d.SetGeometryOrder(4);
  
  mesh3d.ReallocateVertices(4);
  mesh3d.Vertex(0).Init(-0.2, 0.2, 0.1);
  mesh3d.Vertex(1).Init(1, 0, 0);
  mesh3d.Vertex(2).Init(0.5, 0.5*sqrt(Real_wp(3.0)), 0);
  mesh3d.Vertex(3).Init(0.2, 0.6, sqrt(Real_wp(0.6)));
  
  mesh3d.ReallocateBoundariesRef(4);
  mesh3d.BoundaryRef(0).InitTriangular(0, 1, 2, 1);
  mesh3d.BoundaryRef(1).InitTriangular(1, 2, 3, 2);
  mesh3d.BoundaryRef(2).InitTriangular(0, 1, 3, 1);
  mesh3d.BoundaryRef(3).InitTriangular(0, 2, 3, 1);
  
  mesh3d.ReallocateElements(1);
  mesh3d.Element(0).InitTetrahedral(0, 1, 2, 3, 1);
  
  mesh3d.SetCurveType(2, mesh3d.CURVE_SPHERE);
  param.Reallocate(4);
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
    cout << "Testing TetrahedronClassical " << endl;
    TetrahedronClassical tet; string name_elt("TETRAHEDRON_CLASSICAL");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_tet);
  }

  {
    cout << "Testing TetrahedronHierarchic " << endl;
    TetrahedronHierarchic tet;  string name_elt("TETRAHEDRON_HIERARCHIC");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_tet);
  }

  {
    cout << "Testing TetrahedronDgOrtho " << endl;
    TetrahedronDgOrtho tet; string name_elt("TETRAHEDRON_DG_ORTHO");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_tet, false, false);
  }
  
  mesh3d.Clear();
  mesh3d.ReallocateVertices(5);
  mesh3d.Vertex(0).Init(1, 0, 0);
  mesh3d.Vertex(1).Init(0.5, 0.5*sqrt(Real_wp(3.0)), 0);
  mesh3d.Vertex(2).Init(0.8, 0.2, sqrt(Real_wp(0.68)));
  mesh3d.Vertex(3).Init(0.7, -0.3, sqrt(Real_wp(0.58)));
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
    cout << "Testing PyramidLobatto " << endl;
    PyramidLobatto tet; string name_elt("TETRAHEDRON_LOBATTO");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_pyr);
  }
  
  {
    cout << "Testing PyramidClassical " << endl;
    PyramidClassical tet; string name_elt("TETRAHEDRON_CLASSICAL");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_pyr);
  }
  
  {
    cout << "Testing PyramidHierarchic " << endl;
    PyramidHierarchic tet; string name_elt("TETRAHEDRON_HIERARCHIC");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_pyr);
  }

  {
    cout << "Testing PyramidDgOrtho " << endl;
    PyramidDgOrtho tet; string name_elt("TETRAHEDRON_DG_ORTHO");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_pyr, false, false);
  }
  
  {
    cout << "Testing PyramidDgLegendre " << endl;
    PyramidDgLegendre tet; string name_elt("TETRAHEDRON_DG_LEGENDRE");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_pyr, false, false);
  }

  mesh3d.Clear(); mesh3d.ClearCurves();
  mesh3d.ReallocateVertices(6);
  mesh3d.Vertex(0).Init(0.0, 0.0, 0.0);
  mesh3d.Vertex(1).Init(1.0, 0.0, 0.0);
  mesh3d.Vertex(2).Init(0.0, 1.0, 0.0);
  //mesh3d.Vertex(3).Init(0.12, 0.03, 0.98);
  //mesh3d.Vertex(4).Init(1.32, 0.07, 0.96);
  //mesh3d.Vertex(5).Init(0.1, 1.2, 0.8);
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
    cout << "Testing WedgeClassical " << endl;
    WedgeClassical tet; string name_elt("TETRAHEDRON_CLASSICAL");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_wed);
  }

  {
    cout << "Testing WedgeLobatto " << endl;
    WedgeLobatto tet; string name_elt("TETRAHEDRON_LOBATTO");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_wed);
  }

  {
    cout << "Testing WedgeHierarchic " << endl;
    WedgeHierarchic tet; string name_elt("TETRAHEDRON_HIERARCHIC");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_wed, true, false);
  }

  {
    cout << "Testing WedgeDgClassical " << endl;
    WedgeDgClassical tet;  string name_elt("TETRAHEDRON_DG_CLASSICAL");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_wed, false, false);
  }

  {
    cout << "Testing WedgeDgOrtho " << endl;
    WedgeDgOrtho tet; string name_elt("TETRAHEDRON_DG_ORTHO");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_wed, false, false);
  }

  {
    cout << "Testing WedgeDgLegendre " << endl;
    WedgeDgLegendre tet;  string name_elt("TETRAHEDRON_DG_LEGENDRE");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_wed, false, false);
  }
  
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
    cout << "Testing HexahedronGauss " << endl;
    HexahedronGauss tet; string name_elt("TETRAHEDRON_CLASSICAL");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_hex);
  }

  {
    cout << "Testing HexahedronHierarchic " << endl;
    HexahedronHierarchic tet; string name_elt("TETRAHEDRON_HIERARCHIC");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_hex, true, false);
  }

  {
    cout << "Testing HexahedronDgGauss " << endl;
    HexahedronDgGauss tet;  string name_elt("TETRAHEDRON_DG_CLASSICAL");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_hex, false, false);
  }

  {
    cout << "Testing HexahedronLobatto " << endl;
    HexahedronLobatto tet; string name_elt("TETRAHEDRON_LOBATTO");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_hex);
  }

  {
    cout << "Testing HexahedronDgLegendre " << endl;
    HexahedronDgLegendre tet; string name_elt("TETRAHEDRON_DG_LEGENDRE");
    CheckContinuousFiniteElement(tet, name_elt, mesh3d, vertex_hex, false, false);
  }
#endif
  
  return FinalizeMontjoie();
}
