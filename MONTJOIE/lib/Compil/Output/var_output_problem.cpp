#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Output/OutputHarmonic.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template void VarOutputProblem_Base::WriteDatas(VectReal_wp&);
  SELDON_EXTERN template void VarOutputProblem_Base::WriteDatas(VectComplex_wp&);
  
  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFileReshaped(const Vector<VectReal_wp>&, const string&) const;
  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFileReshaped(const Vector<VectComplex_wp>&, const string&) const;

  SELDON_EXTERN template void VarOutputProblem_Base::ReshapeVector(const Vector<Real_wp>&, Vector<VectReal_wp>&, int ) const;
  SELDON_EXTERN template void VarOutputProblem_Base::ReshapeVector(const Vector<Complex_wp>&, Vector<VectComplex_wp>&, int) const;

  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFile(const VectReal_wp&, const string&, int) const;
  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFile(const VectComplex_wp&, const string&, int) const;
  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFile(const Vector<VectReal_wp>&, int, int) const;
  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFile(const Vector<VectComplex_wp>&, int, int) const;
  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFile(const VectReal_wp&, int, int, const Vector<int>&) const;
  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFile(const VectComplex_wp&, int, int, const Vector<int>&) const;
  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFile(const VectReal_wp&, const string&, int, const Vector<int>&) const;
  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFile(const VectComplex_wp&, const string&, int, const Vector<int>&) const;

  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFile(const VectReal_wp&, int, int) const;
  SELDON_EXTERN template void VarOutputProblem_Base::WriteOutputFile(const VectComplex_wp&, int, int) const;

  SELDON_EXTERN template class VarOutputProblem_Dim<Dimension2>;
  SELDON_EXTERN template class VarOutputProblem_Dim<Dimension3>;
  
  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension2>::
  ComputeInterpolationU_GradU(const Vector<VectReal_wp> &, const GridInterpolation<Dimension2>&, const GridInterpolationFull<Dimension2>&,
			      Vector<VectReal_wp>&, Vector<VectReal_wp>&, int, const Real_wp&,
			      Vector<VectReal_wp>& trace_u_diff, Vector<VectReal_wp>& trace_grad_diff, bool) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension2>::
  ComputeInterpolationU_GradU(const Vector<VectComplex_wp> &, const GridInterpolation<Dimension2>&, const GridInterpolationFull<Dimension2>&,
			      Vector<VectComplex_wp>&, Vector<VectComplex_wp>&, int, const Real_wp&,
			      Vector<VectComplex_wp>& trace_u_diff, Vector<VectComplex_wp>& trace_grad_diff, bool) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension3>::
  ComputeInterpolationU_GradU(const Vector<VectReal_wp> &, const GridInterpolation<Dimension3>&, const GridInterpolationFull<Dimension3>&,
			      Vector<VectReal_wp>&, Vector<VectReal_wp>&, int, const Real_wp&,
			      Vector<VectReal_wp>& trace_u_diff, Vector<VectReal_wp>& trace_grad_diff, bool) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension3>::
  ComputeInterpolationU_GradU(const Vector<VectComplex_wp> &, const GridInterpolation<Dimension3>&, const GridInterpolationFull<Dimension3>&,
			      Vector<VectComplex_wp>&, Vector<VectComplex_wp>&, int, const Real_wp&,
			      Vector<VectComplex_wp>& trace_u_diff, Vector<VectComplex_wp>& trace_grad_diff, bool) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension2>::
  ComputeInterpolationU_GradU(const Vector<Real_wp> & U0, int num_grid, bool compute_grad,
                              Vector<Vector<Real_wp> >& trace_u, Vector<Vector<Real_wp> >& trace_grad_u, int) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension2>::
  ComputeInterpolationU_GradU(const Vector<Complex_wp> & U0, int num_grid, bool compute_grad,
                              Vector<Vector<Complex_wp> >& trace_u, Vector<Vector<Complex_wp> >& trace_grad_u, int) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension3>::
  ComputeInterpolationU_GradU(const Vector<Real_wp> & U0, int num_grid, bool compute_grad,
                              Vector<Vector<Real_wp> >& trace_u, Vector<Vector<Real_wp> >& trace_grad_u, int) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension3>::
  ComputeInterpolationU_GradU(const Vector<Complex_wp> & U0, int num_grid, bool compute_grad,
                              Vector<Vector<Complex_wp> >& trace_u, Vector<Vector<Complex_wp> >& trace_grad_u, int) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension2>::
  ComputeInterpolationNodalU(const Vector<Vector<VectComplex_wp> >&,
			     const Vector<Vector<VectComplex_wp> >&,
			     const GridInterpolationFull<Dimension2>&, int, Vector<VectComplex_wp>&) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension2>::
  ComputeInterpolationNodalU(const Vector<Vector<VectReal_wp> >&,
			     const Vector<Vector<VectReal_wp> >&,
			     const GridInterpolationFull<Dimension2>&, int, Vector<VectReal_wp>&) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension3>::
  ComputeInterpolationNodalU(const Vector<Vector<VectComplex_wp> >&,
			     const Vector<Vector<VectComplex_wp> >&,
			     const GridInterpolationFull<Dimension3>&, int, Vector<VectComplex_wp>&) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension3>::
  ComputeInterpolationNodalU(const Vector<Vector<VectReal_wp> >&,
			     const Vector<Vector<VectReal_wp> >&,
			     const GridInterpolationFull<Dimension3>&, int, Vector<VectReal_wp>&) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension2>::
  WriteSnapshot(int, const Real_wp&, Vector<Real_wp>&, const Vector<MeshInterpolation<Dimension2> >&, bool) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension2>::
  WriteSnapshot(int, const Real_wp&, Vector<Complex_wp>&, const Vector<MeshInterpolation<Dimension2> >&, bool) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension3>::
  WriteSnapshot(int, const Real_wp&, Vector<Real_wp>&, const Vector<MeshInterpolation<Dimension3> >&, bool) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension3>::
  WriteSnapshot(int, const Real_wp&, Vector<Complex_wp>&, const Vector<MeshInterpolation<Dimension3> >&, bool) const;
  
  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension2>::
  WriteSnapshot(int, const Real_wp&, Vector<Real_wp>&, const ParamOutputClass&, const Vector<GridInterpolationFull<Dimension2> >&, bool) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension2>::
  WriteSnapshot(int, const Real_wp&, Vector<Complex_wp>&, const ParamOutputClass&, const Vector<GridInterpolationFull<Dimension2> >&, bool) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension3>::
  WriteSnapshot(int, const Real_wp&, Vector<Real_wp>&, const ParamOutputClass&, const Vector<GridInterpolationFull<Dimension3> >&, bool) const;

  SELDON_EXTERN template void VarOutputProblem_Dim<Dimension3>::
  WriteSnapshot(int, const Real_wp&, Vector<Complex_wp>&, const ParamOutputClass&, const Vector<GridInterpolationFull<Dimension3> >&, bool) const;
  

  SELDON_EXTERN template class VarOutputProblem<Dimension2>;
  SELDON_EXTERN template class VarOutputProblem<Dimension3>;

  SELDON_EXTERN template void VarOutputProblem<Dimension2>
  ::ComputeInterpolationUlocGen(const Vector<VectReal_wp>&, const GridInterpolation<Dimension2>&, VectReal_wp&, VectReal_wp&, const IVect&, int, bool) const;

  SELDON_EXTERN template void VarOutputProblem<Dimension2>
  ::ComputeInterpolationUlocGen(const Vector<VectComplex_wp>&, const GridInterpolation<Dimension2>&, VectComplex_wp&, VectComplex_wp&, const IVect&, int, bool) const;


  SELDON_EXTERN template void VarOutputProblem<Dimension2>
  ::ComputeNodalUgradU_Gen(const Vector<VectReal_wp>&, Vector<Vector<VectReal_wp> >&,
			   Vector<Vector<VectReal_wp> >&, bool , bool , bool, bool, int) const;

  SELDON_EXTERN template void VarOutputProblem<Dimension2>
  ::ComputeNodalUgradU_Gen(const Vector<VectComplex_wp>&, Vector<Vector<VectComplex_wp> >&,
			   Vector<Vector<VectComplex_wp> >&, bool , bool , bool, bool, int) const;

  SELDON_EXTERN template void VarOutputProblem<Dimension2>
  ::ComputeQuadratureUgradU_Gen(const Vector<VectReal_wp>&, Vector<Vector<VectReal_wp> >&,
				Vector<Vector<VectReal_wp> >&, bool, bool, bool, int) const;

  SELDON_EXTERN template void VarOutputProblem<Dimension2>
  ::ComputeQuadratureUgradU_Gen(const Vector<VectComplex_wp>&, Vector<Vector<VectComplex_wp> >&,
				Vector<Vector<VectComplex_wp> >&, bool, bool, bool, int) const;


  SELDON_EXTERN template void VarOutputProblem<Dimension3>
  ::ComputeInterpolationUlocGen(const Vector<VectReal_wp>&, const GridInterpolation<Dimension3>&, VectReal_wp&, VectReal_wp&, const IVect&, int, bool) const;

  SELDON_EXTERN template void VarOutputProblem<Dimension3>
  ::ComputeInterpolationUlocGen(const Vector<VectComplex_wp>&, const GridInterpolation<Dimension3>&, VectComplex_wp&, VectComplex_wp&, const IVect&, int, bool) const;


  SELDON_EXTERN template void VarOutputProblem<Dimension3>
  ::ComputeNodalUgradU_Gen(const Vector<VectReal_wp>&, Vector<Vector<VectReal_wp> >&,
			   Vector<Vector<VectReal_wp> >&, bool , bool , bool, bool, int) const;

  SELDON_EXTERN template void VarOutputProblem<Dimension3>
  ::ComputeNodalUgradU_Gen(const Vector<VectComplex_wp>&, Vector<Vector<VectComplex_wp> >&,
			   Vector<Vector<VectComplex_wp> >&, bool , bool , bool, bool, int) const;

  SELDON_EXTERN template void VarOutputProblem<Dimension3>
  ::ComputeQuadratureUgradU_Gen(const Vector<VectReal_wp>&, Vector<Vector<VectReal_wp> >&,
				Vector<Vector<VectReal_wp> >&, bool, bool, bool, int) const;

  SELDON_EXTERN template void VarOutputProblem<Dimension3>
  ::ComputeQuadratureUgradU_Gen(const Vector<VectComplex_wp>&, Vector<Vector<VectComplex_wp> >&,
				Vector<Vector<VectComplex_wp> >&, bool, bool, bool, int) const;

}
