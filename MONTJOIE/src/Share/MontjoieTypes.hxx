#ifndef MONTJOIE_FILE_MONTJOIE_TYPES_HXX

// some gsl files
#ifdef MONTJOIE_WITH_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>
#endif

namespace Montjoie
{

  // prior declaration of edge, face, volume (entities of a mesh)
  template<class Dimension> class Edge {};
  template<class Dimension> class Face {};
  class Volume;
  class ElementNumbering;
  template<class Dimension> class ElementGeomReference {};
  template<class Dimension> class ElementGeomReference_Base;
  template<class Dimension> class ElementReference_Dim;
  template<class Dimension, int type> class ElementReference;
  
  class EdgeLobatto;
  class TriangleGeomReference;
  class QuadrangleGeomReference;
  
#ifdef MONTJOIE_WITH_THREE_DIM
  class TetrahedronGeomReference;
  class PyramidGeomReference;
  class WedgeGeomReference;
  class HexahedronGeomReference;
#endif
  
  template<int t = 0> class TriangleReference;
  template<int t = 0> class QuadrangleReference;
  template<int t = 0> class FaceReference;

#ifdef MONTJOIE_WITH_THREE_DIM
  template<int t = 0> class VolumeReference;
  template<int t = 0> class TetrahedronReference;
  template<int t = 0> class PyramidReference;
  template<int t = 0> class WedgeReference;
  template<int t = 0> class HexahedronReference;
#endif

  class FiniteElementInterpolator;
  
  // VectReal_wp = real vector, VectComplex_wp = complex vector
  // VectBool = boolean bector
  //! vector of real numbers
  typedef Vector<Real_wp> VectReal_wp;
  //! vector of complex numbers
  typedef Vector<Complex_wp> VectComplex_wp;
  //! vector of boolean
  typedef Vector<bool> VectBool;
  //! vector of strings
  typedef Vector<string> VectString;
    
  // small vectors in R^2
  //! real vector of length 2
  typedef TinyVector<Real_wp, 2> R2;
  //! complex vector of length 2
  typedef TinyVector<Complex_wp, 2> R2_Complex_wp;
  
  // VectR2 : array of points in R^2
  //! vector of points in R2
  typedef Vector<R2> VectR2;
  //! complex vector of points in R2
  typedef Vector<R2_Complex_wp> VectR2_Complex_wp;
  
  //! real matrices 2x2
  typedef TinyMatrix<Real_wp, General, 2, 2> Matrix2_2;
  //! real symmetric matrices 2x2
  typedef TinyMatrix<Real_wp, Symmetric, 2, 2> Matrix2_2sym;
  //! complex matrices 2x2
  typedef TinyMatrix<Complex_wp, General, 2, 2> Matrix2_2_Complex_wp;
  //! complex symmetric matrices 2x2
  typedef TinyMatrix<Complex_wp, Symmetric, 2, 2> Matrix2_2sym_Complex_wp;
  
  // similar things in 3-D
  //! real vector of length 3
  typedef TinyVector<Real_wp, 3> R3;
  //! complex vector of length 3
  typedef TinyVector<Complex_wp, 3> R3_Complex_wp;
  
  //! vector of points in R3
  typedef Vector<R3> VectR3;
  //! complex vector of points in R3
  typedef Vector<R3_Complex_wp> VectR3_Complex_wp;
    
  //! real matrices 3x3
  typedef TinyMatrix<Real_wp, General, 3, 3> Matrix3_3;
  //! real symmetric matrices 3x3
  typedef TinyMatrix<Real_wp, Symmetric, 3, 3> Matrix3_3sym;
  //! complex matrices 3x3
  typedef TinyMatrix<Complex_wp, General, 3, 3> Matrix3_3_Complex_wp;
  //! complex symmetric matrices 3x3
  typedef TinyMatrix<Complex_wp, Symmetric, 3, 3> Matrix3_3sym_Complex_wp;
  
  //! base class for 1-D
  /*!
    This class is intended to be a parameter of template classes
    It defines basic objects (Rn, elements, etc) in 1-D.
   */
  class Dimension1
  {
  public :
    enum{ dim_N = 1}; //!< dimension
    
    // d = dimension = 1 here
    typedef Real_wp R_N; //!< alias for float vector of size d 
    typedef Complex_wp R_N_Complex_wp; //!< alias for complex vector of size d
    typedef VectReal_wp VectR_N; //!< alias for vector of points in R^d
    //! alias for vector of complex points in R^d
    typedef VectComplex_wp VectR_N_Complex_wp;
    
    typedef Real_wp MatrixN_N; //!< alias for d x d matrix
    typedef Real_wp MatrixN_Nsym;  //!< alias for d x d symmetric matrix
    //! alias for d x d complex matrix 
    typedef Complex_wp MatrixN_N_Complex_wp; 
    //! alias for d x d symmetric complex matrix
    typedef Complex_wp MatrixN_Nsym_Complex_wp;
    //! alias for vector of d x d matrices
    typedef VectReal_wp VectMatrixN_N;
    //! alias for vector of symmetric d x d matrices
    typedef VectReal_wp VectMatrixN_Nsym;
    
    //! alias for elements in dimension d
    typedef Edge<Dimension1> Element_Rn;
    typedef Dimension1 DimensionBoundary; //!< dimension d-1
  };
  
  
  //! base class for 2-D
  /*!
    This class is intended to be a parameter of template classes
    It defines basic objects (Rn, elements, etc) in 2-D.
   */
  class Dimension2
  {
  public :
    enum{ dim_N = 2}; //!< dimension
    
    // d = dimension = 2 here
    typedef R2 R_N; //!< alias for float vector of size d 
    typedef R2_Complex_wp R_N_Complex_wp; //!< alias for complex vector of size d
    typedef VectR2 VectR_N; //!< alias for vector of points in R^d
    //! alias for vector of complex points in R^d
    typedef VectR2_Complex_wp VectR_N_Complex_wp;

    typedef Matrix2_2 MatrixN_N; //!< alias for d x d matrix
    typedef Matrix2_2sym MatrixN_Nsym; //!< alias for d x d symmetric matrix
    //! alias for d x d complex matrix 
    typedef Matrix2_2_Complex_wp MatrixN_N_Complex_wp;
    //! alias for d x d symmetric complex matrix
    typedef Matrix2_2sym_Complex_wp MatrixN_Nsym_Complex_wp;
    
    //! alias for elements in dimension d
    typedef Face<Dimension2> Element_Rn;
    typedef Dimension1 DimensionBoundary; //!< dimension d-1
  };
  
  
  //! base class for 3-D
  /*!
    This class is intended to be a parameter of template classes
    It defines basic objects (Rn, elements, etc) in 3-D.
   */
  class Dimension3
  {
  public :
    enum{ dim_N = 3}; //!< dimension
    
    // d = dimension = 2 here
    typedef R3 R_N; //!< alias for float vector of size d 
    typedef R3_Complex_wp R_N_Complex_wp; //!< alias for complex vector of size d
    typedef VectR3 VectR_N; //!< alias for vector of points in R^d
    //! alias for vector of complex points in R^d
    typedef VectR3_Complex_wp VectR_N_Complex_wp;

    typedef Matrix3_3 MatrixN_N; //!< alias for d x d matrix
    typedef Matrix3_3sym MatrixN_Nsym; //!< alias for d x d symmetric matrix
    //! alias for d x d complex matrix 
    typedef Matrix3_3_Complex_wp MatrixN_N_Complex_wp;
    //! alias for d x d symmetric complex matrix
    typedef Matrix3_3sym_Complex_wp MatrixN_Nsym_Complex_wp;
    
    //! alias for elements in dimension d
    typedef Volume Element_Rn;
    typedef Dimension2 DimensionBoundary; //!< dimension d-1
  };
  
  template<class T> class GlobalGenericMatrix;
  
  template<class Prop>
  class GenericStorage {};
  
  template<>
  class GenericStorage<General>
  {
  public :
    typedef BlockDiagRow BlockDiagonal;
    typedef RowMajor DenseStorage;
  };

  template<>
  class GenericStorage<Symmetric>
  {
  public :
    typedef BlockDiagRowSym BlockDiagonal;
    typedef RowSymPacked DenseStorage;
  };

  template<>
  class GenericStorage<Real_wp>
  {
  public :
    typedef DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> SparseSymMatrix;
    typedef DistributedMatrix<Real_wp, General, ArrayRowSparse> SparseMatrix;
    typedef DistributedMatrix<Real_wp, Symmetric, RowSymSparse> SparseSymMatrixCSR;
    typedef DistributedMatrix<Real_wp, General, RowSparse> SparseMatrixCSR;
  };

  template<>
  class GenericStorage<Complex_wp>
  {
  public :
    typedef DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> SparseSymMatrix;
    typedef DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse> SparseMatrix;
    typedef DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse> SparseSymMatrixCSR;
    typedef DistributedMatrix<Complex_wp, General, RowComplexSparse> SparseMatrixCSR;
  };

  template<>
  class GenericStorage<ArrayRowSparse>
  {
  public :
    typedef RowSparse StorageCSR;
  };

  template<>
  class GenericStorage<ArrayRowSymSparse>
  {
  public :
    typedef RowSymSparse StorageCSR;
  };
  
  /********************************************
   * Prior declaration of some useful classes *
   ********************************************/
  
  // mesh classes
  template<class Dimension, class T = Real_wp> class Mesh {};
  template<class Dimension, class T = Real_wp> class MeshNumbering {};
  template<class Dimension, class T = Real_wp> class MeshBoundaries {};
  template<class Dimension, class T = Real_wp> class SurfacicMesh {};
  template<class Dimension, class T = Real_wp> class SurfacicMeshIntegration {};
  template<class Dimension> class GridInterpolation {};
  template<class Dimension> class GridInterpolationFull {};
  
  // class defining a varying scalar field
  template<class Dimension, class T>
  class PhysicalVaryingMedia;
  
  // declarations of base classes for EllipticProblem
  class VarProblem_Base;
  class DistributedProblem_Base;
  class VarComputationProblem_Base;
  class VarComputationProblem;
  class VarBoundaryCondition_Base;
  class VarSourceProblem_Base;
  class VarOutputProblem_Base;

  class VarInstationary_Base;
  
  template<class Dimension> class VarGeometryProblem;
  template<class Dimension> class DistributedProblem;

  template<class Dimension> class VarProblem;
  template<class Dimension> class VarOutputProblem;
  template<class Dimension> class VarInstationary_Fem;
  
  template<class T, class Dimension> class VarSourceProblem;
  template<class T, class Dimension> class VarBoundaryCondition;
  
  template<class TypeEquation> class VarHarmonic;

#ifdef MONTJOIE_WITH_TRANSMISSION
  template<class Dimension>
  class VarTransmission_Base;
#endif

  class VarGeneralizedImpedance_Base;

  template<class Dimension> class VarComputationRCS_Base;
  template<class TypeEquation> class VarComputationRCS;
  
  //! empty class that will be specialized for each equation
  /*!
    EllipticProblem is the main class when solving time-harmonic or static
    problems. Since, there is no "generic" equation implemented in Montjoie, 
    this class is specialized for each equation solved by Montjoie. When
    you specify a new equation, you will have to specialize this class 
    for your own class defining the equation.
   */
  template<class TypeEquation> class EllipticProblem
  {
  };


  //! empty class that will be specialized for each equation
  /*!
    HyperbolicProblem is the main class when solving unsteady (time-dependent)
    problems. Since, there is no "generic" equation implemented in Montjoie, 
    this class is specialized for each equation solved by Montjoie. When
    you specify a new equation, you will have to specialize this class 
    for your own class defining the equation.
   */
  template<class TypeEquation> class HyperbolicProblem
  {
  };
  
  
  // declaration of class defining global and mass matrix for finite element
  template<class T> class All_Preconditioner_Base;
  template<class T> class FemMatrixFreeClass_Base;
  template<class T, class TypeEquation> class FemMatrixFreeClass;
    
#ifdef SELDON_WITH_PRECONDITIONING
  template<class T, class TypeEquation> class MultigridPreconditioning;
#endif
  
  // Polynomial classes
  template<class T, class Allocator>
  class UnivariatePolynomial;
  
  template<class T, class Allocator>
  class MultivariatePolynomial;
  
  // objects containing points and jacobian matrices
  template<class Dimension>
  class SetPoints;
  
  template<class Dimension>
  class SetMatrices;
  
  // projection on tensorized elements
  template<class Dimension>
  class TensorizedProjector {};
  
  //! class specialized for Dimension2 and Dimension3 
  template<class Dimension>
  class FiniteElementHcurl {};
  
  // physical property
  template<class Dimension, class T>
  class GenericPhysicalIndice;
  
  //! class used to avoid if structures when value is known at compilation
  /*!
    instead of writing :
    if (condition)
       statement 1
    else
       statement 2
    
    we sometimes write :
    
    GhostIf<condition> cond_true;
    CallOverloadedFunction(args, ..., cond_true);

    and we have two definitions of this function :
    
    void CallOverloadedFunction(args, ..., GhostIf<true>& )
    {
      statement 1
    }

    void CallOverloadedFunction(args, ..., GhostIf<false>& )
    {
      statement 2
    }
   */
  template<bool>
  class GhostIf {};
  
  // class for definition of source, by default null source is set
  template<class T, class Dimension>
  class VirtualSourceFEM;

  template<class T, class Dimension>
  class IncidentWaveField;

  template<class T> class VirtualTimeSource;
  
  //! definition source when diffracted field is computed
  /*!
    Since the definition of this source depends on the solved equation (Helmholtz/Maxwell)
    the generic class is empty, you have to specialize this class if you 
    know how to define source induced by the computation of diffracted field for your own equation
  */
  template<class TypeEquation>
  class DiffractedWaveSource : public VirtualSourceFEM<typename TypeEquation::Complexe,
						       typename TypeEquation::Dimension>
  {
  public :
    
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& var,
			 IncidentWaveField<typename TypeEquation::Complexe,
			 typename TypeEquation::Dimension>& fsrc_);
    
  };
  
  //! definition source when total field is computed
  /*!
    Since the definition of this source depends on the solved equation (Helmholtz/Maxwell)
    the generic class is empty, you have to specialize this class if you 
    know how to define source induced by the computation of total field for your own equation
  */
  template<class TypeEquation>
  class TotalWaveSource  : public VirtualSourceFEM<typename TypeEquation::Complexe,
						   typename TypeEquation::Dimension>
  {
  public :  
    
    inline TotalWaveSource(const EllipticProblem<TypeEquation>& var,
			   IncidentWaveField<typename TypeEquation::Complexe,
			   typename TypeEquation::Dimension>& fsrc_);
    
  };
  
  
#ifdef SELDON_WITH_SCALAPACK
  extern BlacsHandle global_blacs_handle;
#endif
  
  void InitMontjoie(int argc, char** argv, int lvl = 1);
  
  int FinalizeMontjoie();
  
} // namespace Montjoie


namespace Seldon
{
  double Norm2(const complex<double>& c);
  double Norm2(const double& c);

#ifdef MONTJOIE_WITH_MULTIPLE
  Real_wp Norm2(const Complex_wp& c);
  Real_wp Norm2(const Real_wp& c);
#endif
  
}

#define MONTJOIE_FILE_MONTJOIE_TYPES_HXX
#endif

