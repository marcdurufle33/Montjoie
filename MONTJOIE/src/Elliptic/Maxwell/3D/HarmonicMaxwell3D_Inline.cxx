#ifndef MONTJOIE_FILE_HARMONIC_MAXWELL_3D_INLINE_CXX

namespace Montjoie
{

  //! default constructor
  template<class TypeEquation>
  inline HarmonicMaxwell_3D_Base
  ::HarmonicMaxwell_3D_Base(EllipticProblem<TypeEquation>& var)
    : var_problem(var), var_computation(var), var_boundary(var), var_output(var)
#ifdef MONTJOIE_WITH_TRANSMISSION
    , var_transmission_base(var.var_transmission)
#endif
  {
    nb_dof_drude_vec = 0; nb_dof_drude_vec_all = 0;
    this->use_symm_drude = false; this->linearize_drude = false;
    var.GetMeshNumbering(0).drop_interface_pml_dof = true;
  }
  

  template<class Complexe>
  inline void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::SetAh(int i, const TinyMatrix<Complexe, Symmetric, 3, 3>& M)
  {
    Ah_Bh(offset_Ah + i) = M;
  }


  template<class Complexe>
  inline void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::SetBh(int i, const TinyMatrix<Complexe, Symmetric, 3, 3>& M)
  {
    Ah_Bh(offset_Bh + i) = M;
  }


  template<class Complexe>
  inline void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::SetBhSigma(int i, const TinyMatrix<Complexe, Symmetric, 3, 3>& M)
  {
    Ah_Bh(offset_BhSigma + i) = M;
  }

  
  template<class Complexe>
  inline void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::SetBhSurf(int i, const TinyMatrix<Complexe, Symmetric, 2, 2>& M)
  {
    BhSurf(i) = M;
  }


  template<class Complexe>
  inline const TinyMatrix<Complexe, Symmetric, 3, 3>& HarmonicMaxwell3D_PhysGeomInfo<Complexe>::GetBh(int i) const
  {
    return Ah_Bh(offset_Bh + i);
  }
    

  template<class Complexe>
  inline const TinyMatrix<Complexe, Symmetric, 3, 3>& HarmonicMaxwell3D_PhysGeomInfo<Complexe>::GetBhSigma(int i) const
  {
    return Ah_Bh(offset_BhSigma + i);
  }
  

  template<class Complexe>
  inline const TinyMatrix<Complexe, Symmetric, 3, 3>& HarmonicMaxwell3D_PhysGeomInfo<Complexe>::GetAh(int i) const
  {
    return Ah_Bh(offset_Ah + i);
  }


  template<class Complexe>
  inline const TinyMatrix<Complexe, Symmetric, 2, 2>& HarmonicMaxwell3D_PhysGeomInfo<Complexe>::GetBhSurf(int i) const
  {
    return BhSurf(i);
  }


  template<class Complexe>
  inline bool HarmonicMaxwell3D_PhysGeomInfo<Complexe>::IsDiagonalMass() const
  {
    if (diag_Dh.GetM() <= 0)
      return false;

    return true;
  }
  
  
  template<class Complexe>
  inline void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::SetDh(int i, const Complexe& d)
  {
    diag_Dh(i) = d;
  }


  template<class Complexe>
  inline void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::SetDhSigma(int i, const Complexe& d)
  {
    diag_Dh(offset_diagDhSigma + i) = d;
  }


  template<class Complexe>
  inline void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::SetDhStiff(int i, const Complexe& d)
  {
    diag_Dh(offset_diagDhStiff + i) = d;
  }


  template<class Complexe>
  inline void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::SetTauPML(int i, const TinyVector<Real_wp, 3>& t)
  {
    tauPML(i) = t;
  }


  template<class Complexe>
  inline const Complexe& HarmonicMaxwell3D_PhysGeomInfo<Complexe>::GetDh(int i) const
  {
    return diag_Dh(i);
  }


  template<class Complexe>
  inline const Complexe& HarmonicMaxwell3D_PhysGeomInfo<Complexe>::GetDhSigma(int i) const
  {
    return diag_Dh(offset_diagDhSigma + i);
  }


  template<class Complexe>
  inline const Complexe& HarmonicMaxwell3D_PhysGeomInfo<Complexe>::GetDhStiff(int i) const
  {
    return diag_Dh(offset_diagDhStiff + i);
  }


  template<class Complexe>
  inline const TinyVector<Real_wp, 3>& HarmonicMaxwell3D_PhysGeomInfo<Complexe>::GetTauPML(int i) const
  {
    return tauPML(i);
  }


  //! default constructor
  template<class Complexe> template<class TypeEquation>
  inline HarmonicMaxwell_3D<Complexe>
  ::HarmonicMaxwell_3D(EllipticProblem<TypeEquation>& var)
    : HarmonicMaxwell_3D_Base(var),
      var_problem(var), var_computation(var), var_output(var), var_boundary(var),
      var_source_base(var), var_source(var)
  {
    epsilon0 = 1.0; mu0 = 1.0; invMu0 = 1.0;
  }
  

  template<class TypeEquation>
  inline HarmonicMaxwell_3D_Eq<TypeEquation>::HarmonicMaxwell_3D_Eq()
    : VarHarmonic<TypeEquation>(), HarmonicMaxwell_3D<Complexe>(this->GetLeafClass())
  {
  }
  

  //! Constructor with the given Maxwell problem
  template<class T, class Complexe> template<class TypeEquation>
  inline ImpedanceABC_Maxwell3D<T, Complexe>
  ::ImpedanceABC_Maxwell3D(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<T, Dimension3>(var),
      var_maxwell(var), var_boundary(var), var_problem(var)
  {
  }
  

  template<class TypeEquation>
  inline ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_3D>
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceABC_Maxwell3D<Complex_wp, Complex_wp>(var) 
  {
  }

  
  template<class TypeEquation>
  inline TransparencySolver_Maxwell3D::
  TransparencySolver_Maxwell3D(EllipticProblem<TypeEquation>& var,
			       All_LinearSolver& solver)
    : TransparencySolver_Fem<Dimension3, 3, 3>(var,solver), var_problem(var)
  {
  }
  
  
  template<class T>
  inline CondensationBlockSolver<T, HarmonicMaxwellEquation_3D>
  ::CondensationBlockSolver(EllipticProblem<HarmonicMaxwellEquation_3D>& var)
    : CondensationBlockSolver_Maxwell3D<T, HarmonicMaxwellEquation_3D>(var)
  {
  }
  
  
  template<class T, class Complexe> template<class TypeEquation>
  inline CondensationBlockSolver_MaxwellHdg3D<T, Complexe>
  ::CondensationBlockSolver_MaxwellHdg3D(EllipticProblem<TypeEquation>& var)
    : CondensationBlockSolver_Fem<T>(var) , var_problem(var)
  {
  }
  

  template<class T>
  inline CondensationBlockSolver<T, HarmonicMaxwellEquationHdg_3D>
  ::CondensationBlockSolver(EllipticProblem<HarmonicMaxwellEquationHdg_3D>& var)
    : CondensationBlockSolver_MaxwellHdg3D<T, Complex_wp>(var)
  {
  }

  
  template<class TypeEquation>
  inline VarComputationRCS_Maxwell3D::VarComputationRCS_Maxwell3D(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_Base<Dimension3>(var) , var_problem(var)
  {
  }
  
  
  template<class TypeEquation>
  inline VarComputationRCS<HarmonicMaxwellEquation_3D>::VarComputationRCS(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_Maxwell3D(var) 
  {
  }
  

  template<class TypeEquation>
  inline VarComputationRCS<HarmonicMaxwellEquation_3D_DG>::VarComputationRCS(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_Maxwell3D(var)
  {
  }


  template<class TypeEquation>
  inline VarComputationRCS<HarmonicMaxwellEquationHdg_3D>::VarComputationRCS(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_Maxwell3D(var)
  {
  }

  inline FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquation_3D>
  ::FemMatrixFreeClass(const EllipticProblem<HarmonicMaxwellEquation_3D>& var_)
    : FemMatrixFreeClass_Eq<Complex_wp, HarmonicMaxwellEquation_3D>(var_)
  {
  }

  inline FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquationHdg_3D>
  ::FemMatrixFreeClass(const EllipticProblem<HarmonicMaxwellEquationHdg_3D>& var_)
    : FemMatrixFreeClass_Eq<Complex_wp, HarmonicMaxwellEquationHdg_3D>(var_)
  {
  }

}

#define MONTJOIE_FILE_HARMONIC_MAXWELL_3D_INLINE_CXX
#endif

