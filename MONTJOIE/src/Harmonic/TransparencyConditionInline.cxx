#ifndef MONTJOIE_FILE_TRANSPARENCY_CONDITION_INLINE_CXX

namespace Montjoie
{
  
  /**************************
   * VarComputationRCS_Base *
   **************************/
  

  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline VarComputationRCS_Base<Dimension>
  ::VarComputationRCS_Base(EllipticProblem<TypeEquation>& var) : var_problem(var)
  {
    InitDefaultValues();
  }
  
  
  template<class Dimension>
  inline bool VarComputationRCS_Base<Dimension>::RcsToBeComputed() const
  {
    return rcs_to_be_computed;
  }
  
  
  template<class Dimension>
  inline int VarComputationRCS_Base<Dimension>::GetNbPointsOutside() const
  {
    return PointsOutside.GetM();
  }


  template<class Dimension>
  inline int VarComputationRCS_Base<Dimension>::GetNbAngles() const
  {
    return nb_angles_RCS; 
  }
  
  
  template<class Dimension>
  inline int VarComputationRCS_Base<Dimension>::GetRcsType() const
  {
    return type_rcs;
  }
  
  
  template<class Dimension>
  inline const MeshInterpolationFEM<Dimension>& VarComputationRCS_Base<Dimension>
  ::GetInterpolationMesh() const
  {
    return var_mesh;
  }
  
  
  template<class Dimension>
  inline const typename Dimension::R_N& VarComputationRCS_Base<Dimension>
  ::GetOutsidePoint(int i) const
  {
    return PointsOutside(i);
  }
  
  
  template<class Dimension>
  inline const typename Dimension::VectR_N& VarComputationRCS_Base<Dimension>
  ::GetOutsidePoint() const
  {
    return PointsOutside;
  }
  
  
  template<class Dimension>
  inline void VarComputationRCS_Base<Dimension>
  ::SetOutsidePoints(const typename Dimension::VectR_N& pts)
  {
    PointsOutside = pts;
  }
  
  
  template<class Dimension>
  inline void VarComputationRCS_Base<Dimension>::SetTimeStep(const Real_wp& dt)
  {
    deltat = dt;
  }

  
  template<class Dimension>
  inline void VarComputationRCS_Base<Dimension>
  ::ComputeIntegralRepresentation(const VectReal_wp& trace_En, const VectReal_wp& trace_Hn,
				  const MeshInterpolationFEM<Dimension>& mesh_interp,
				  const typename Dimension::R_N& pointX, VectReal_wp&) const
  {
    abort();
  }
  
  
  template<class Dimension>
  inline void VarComputationRCS_Base<Dimension>
  ::ComputeIntegralRepresentation(const VectComplex_wp& trace_En, const VectComplex_wp& trace_Hn,
				  const MeshInterpolationFEM<Dimension>& mesh_interp,
				  const typename Dimension::R_N& pointX, VectComplex_wp&) const
  {
    abort(); 
  }


  template<class TypeEquation>
  inline VarComputationRCS<TypeEquation>
  ::VarComputationRCS(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_Base<typename TypeEquation::Dimension>(var) 
  {
  }
    

  template<class TypeEquation>
  inline void VarComputationRCS<TypeEquation>::ComputeRCS(const VectReal_wp& U0)
  {
    abort();
  }
  
  
  template<class TypeEquation>
  inline void VarComputationRCS<TypeEquation>::ComputeRCS(const VectComplex_wp& U0)
  {
    abort();
  }

  
  /***************************
   * TransparencySolver_Base *
   ***************************/


  //! default constructor
  template<class TypeEquation>
  inline TransparencySolver_Base
  ::TransparencySolver_Base(EllipticProblem<TypeEquation>& var, All_LinearSolver& solver)
    : glob_solver(solver), var_problem(var), var_comm(var)
  {
    InitDefaultValues();
  }
  
  
  inline bool TransparencySolver_Base::UseTransparentCondition() const
  {
    return use_transparency_condition;
  }  

  
  inline void TransparencySolver_Base
  ::Solve(VectReal_wp & x_sol, VectReal_wp & b_source)
  {
  }


  inline void TransparencySolver_Base::Reallocate(int, int)
  {
    abort();
  }
  
  
  inline size_t TransparencySolver_Base::GetMemorySize() const
  {
    return 0;
  }
  
  
  inline void TransparencySolver_Base::Clear() 
  {
  }
  
  
  //! no iterations history displayed
  inline void TransparencySolver_Base::HideMessages()
  {
    display_messages = false;
  }
  
  
  //! iterations history displayed
  inline void TransparencySolver_Base::ShowMessages()
  {
    display_messages = true;
  }
  
  
  inline void TransparencySolver_Base::ComputeSolution(VectComplex_wp& rhs, VectComplex_wp& sol) const
  {
    sol = rhs;
    glob_solver.ComputeSolution(sol, nat_mat);
  }
  

  //! Constructor
  template<class Dimension> template<class TypeEquation>
  inline TransparencySolver_Dim<Dimension>
  ::TransparencySolver_Dim(EllipticProblem<TypeEquation>& var, All_LinearSolver& solver)
    : TransparencySolver_Base(var, solver), var_problem(var), var_source(var)
  {
  }
  
  
  template<class Dimension, int nb_en, int nb_hn> template<class TypeEquation>
  inline TransparencySolver_Fem<Dimension, nb_en, nb_hn>
  ::TransparencySolver_Fem(EllipticProblem<TypeEquation>& var, All_LinearSolver& solver) :
    TransparencySolver_Dim<Dimension>(var, solver)
  {
  }
  

  template<class TypeEquation>
  inline TransparencySolver<TypeEquation>
  ::TransparencySolver(EllipticProblem<TypeEquation>& var, All_LinearSolver& solver) :
    TransparencySolver_Fem<Dimension, nb_en, nb_hn>(var, solver) {}
  
  
  template<class TypeEquation>
  inline void TransparencySolver<TypeEquation>
  ::GetSource(const VectComplex_wp&, const VectComplex_wp&, int n,
	      const Real_wp&, const R_N&, const R_N&, Vector<VectComplex_wp>&, int j) const
  { 
    abort();
  }

  
  template<class TypeEquation>
  inline void TransparencySolver<TypeEquation>
  ::ComputeIntegralRepresentation(const VectComplex_wp&, const VectComplex_wp&,
				  const MeshInterpolationFEM<Dimension>&,
				  const R_N&, const R_N&, TinyVector<Complex_wp, nb_en>&,
				  TinyVector<Complex_wp, nb_hn>&) const
  {
    cout << "ComputeIntegralRepresentation not implemented for this equation" << endl;
    abort();
  }
  
}

#define MONTJOIE_FILE_TRANSPARENCY_CONDITION_INLINE_CXX
#endif
