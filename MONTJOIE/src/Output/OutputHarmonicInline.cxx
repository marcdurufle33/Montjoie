#ifndef MONTJOIE_FILE_OUTPUT_HARMONIC_INLINE_CXX

namespace Montjoie
{
 
  //! default constructor
  template<class TypeEquation>
  inline VarOutputProblem_Base::VarOutputProblem_Base(EllipticProblem<TypeEquation>& var)
    : var_problem(var), var_comm(var), var_boundary(var), var_source(var)
  {
    InitDefaultValues();
  }

  
  //! returns true if output files are written in double precision
  inline bool VarOutputProblem_Base::OutputWrittenInDoublePrecision()
  {
    return precision_output_file == OutputTypeEnum::DOUBLE_PRECISION;
  }


  //! initialisation of outputs
  inline void VarOutputProblem_Base::InitOutput(const Real_wp& t0)
  {
    // initializations of outputs
    for (int i = 0; i < output_grid_param.GetM(); i++)
      output_grid_param(i).InitTime(t0);
    
    for (int i = 0; i < output_mesh_param.GetM(); i++)
      output_mesh_param(i).InitTime(t0);
  }


  inline void VarOutputProblem_Base::SetOutputComponent(int t)
  {
    choice_field_output = t;
  }
  
  
  inline int VarOutputProblem_Base::GetOutputComponent() const
  {
    return choice_field_output; 
  }


  inline void VarOutputProblem_Base::KeepGridLocalizationArrays()
  {
    keep_grid_localization_array = true;
  }

  template<class T>
  inline void VarOutputProblem_Base::AddContribInterpol(const T& coef, const TinyVector<Real_wp, 1>& phi, int p, Vector<T>& val)
  {
    val(p) += coef*phi(0);
  }

  template<class T>
  inline void VarOutputProblem_Base::AddContribInterpol(const T& coef, const R2& phi, int p, Vector<T>& val)
  {
    int offset = 2*p;
    val(offset) += coef*phi(0);
    val(offset+1) += coef*phi(1);
  }

  template<class T>
  inline void VarOutputProblem_Base::AddContribInterpol(const T& coef, const R3& phi, int p, Vector<T>& val)
  {
    int offset = 3*p;
    val(offset) += coef*phi(0);
    val(offset+1) += coef*phi(1);
    val(offset+2) += coef*phi(2);
  }
  
  /************************
   * VarOutputProblem_Dim *
   ************************/
  
  
  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline VarOutputProblem_Dim<Dimension>::VarOutputProblem_Dim(EllipticProblem<TypeEquation>& var)
    : VarOutputProblem_Base(var), var_boundary(var), var_source(var),
      var_source_base(var), var_problem(var), output_rcs_param(var.output_rcs_param)
  {
  }

  
  template<class Dimension>
  inline VarComputationRCS_Base<Dimension>& VarOutputProblem_Dim<Dimension>::GetParameterOutputRCS()
  {
    return output_rcs_param;
  }
  
    
  /********************
   * VarOutputProblem *
   ********************/
  
  
  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline VarOutputProblem<Dimension>::VarOutputProblem(EllipticProblem<TypeEquation>& var)
    : VarOutputProblem_Dim<Dimension>(var), var_problem(var), var_source(var), var_boundary(var)
  {
  }


  //! the user can change the values of the solution before output on files
  template<class Dimension>
  inline void VarOutputProblem<Dimension>::
  ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u, int i,
		      const GridInterpolation<Dimension>& var_interp,
		      int iquad, bool compute_grad) const
  {
  }


  //! the user can change the values of the solution before output on files
  template<class Dimension>
  inline void VarOutputProblem<Dimension>::
  ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u, int i,
		      const GridInterpolation<Dimension>& var_interp,
		      int iquad, bool compute_grad) const
  {
  }


  template<class Dimension>
  inline void VarOutputProblem<Dimension>
  ::ModifyOutputUnknown(Vector<VectReal_wp>&, Vector<VectReal_wp>&,
			int, bool, bool) const 
  {
  }
  
  
  template<class Dimension>
  inline void VarOutputProblem<Dimension>
  ::ModifyOutputUnknown(Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
			int, bool, bool) const 
  {
  }
  

  template<class Dimension>
  inline void VarOutputProblem<Dimension>
  ::ComputeInterpolationUloc(const Vector<Vector<Real_wp> > & U0,
			     const GridInterpolation<Dimension> & var_interp,
			     Vector<Real_wp>& trace_vec, Vector<Real_wp>& trace_grad_vec,
			     const IVect& list_points, int nnz, bool compute_grad) const
  {
    ComputeInterpolationUlocGen(U0, var_interp, trace_vec, trace_grad_vec, list_points, nnz, compute_grad);
  }
  
  
  template<class Dimension>
  inline void VarOutputProblem<Dimension>
  ::ComputeInterpolationUloc(const Vector<Vector<Complex_wp> > & U0,
			     const GridInterpolation<Dimension> & var_interp,
			     Vector<Complex_wp>& trace_vec, Vector<Complex_wp>& trace_grad_vec,
			     const IVect& list_points, int nnz, bool compute_grad) const
  {
    ComputeInterpolationUlocGen(U0, var_interp, trace_vec, trace_grad_vec, list_points, nnz, compute_grad);
  }

  
  template<class Dimension>
  inline void VarOutputProblem<Dimension>
  ::ComputeNodalUgradU(const Vector<VectReal_wp>& U0,
		       Vector<Vector<VectReal_wp> >& Unodal,
		       Vector<Vector<VectReal_wp> >& GradNodal,
		       bool u_component, bool grad_component,
		       bool hess_component, bool true_unknown, int nb_u) const
  {
    ComputeNodalUgradU_Gen(U0, Unodal, GradNodal, u_component,
			   grad_component, hess_component, true_unknown, nb_u);
  }
  
  
  template<class Dimension>
  inline void VarOutputProblem<Dimension>
  ::ComputeNodalUgradU(const Vector<VectComplex_wp>& U0,
		       Vector<Vector<VectComplex_wp> >& Unodal,
		       Vector<Vector<VectComplex_wp> >& GradNodal,
		       bool u_component, bool grad_component,
		       bool hess_component, bool true_unknown, int nb_u) const
  {
    ComputeNodalUgradU_Gen(U0, Unodal, GradNodal, u_component,
			   grad_component, hess_component, true_unknown, nb_u);
  }

  
  template<class Dimension>
  inline void VarOutputProblem<Dimension>
  ::ComputeQuadratureUgradU(const Vector<VectReal_wp>& U0,
			    Vector<Vector<VectReal_wp> >& Uquad,
			    Vector<Vector<VectReal_wp> >& GradQuad,
			    bool u_component, bool grad_component,
			    bool true_unknown, int nb_u) const
  {
    ComputeQuadratureUgradU_Gen(U0, Uquad, GradQuad, u_component, grad_component,
				true_unknown, nb_u);
  }
  
  
  template<class Dimension>
  inline void VarOutputProblem<Dimension>
  ::ComputeQuadratureUgradU(const Vector<VectComplex_wp>& U0,
			    Vector<Vector<VectComplex_wp> >& Uquad,
			    Vector<Vector<VectComplex_wp> >& GradQuad,
			    bool u_component, bool grad_component,
			    bool true_unknown, int nb_u) const
  {
    ComputeQuadratureUgradU_Gen(U0, Uquad, GradQuad, u_component, grad_component,
				true_unknown, nb_u);
  }
    
}

#define MONTJOIE_FILE_OUTPUT_HARMONIC_INLINE_CXX
#endif
