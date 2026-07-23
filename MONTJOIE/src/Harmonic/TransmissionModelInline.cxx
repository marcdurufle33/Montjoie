#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_INLINE_CXX

namespace Montjoie
{

  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline VarTransmission_Base<Dimension>::VarTransmission_Base(EllipticProblem<TypeEquation>& var)
    : mesh_plate_num(mesh_plate), var_problem(var), var_boundary(var)
  {
    InitDefaultValues();
  }


  template<class Dimension>
  inline void VarTransmission_Base<Dimension>::InitTransmission()
  {
  }
    

  template<class Dimension>
  inline void VarTransmission_Base<Dimension>::UpdateNumberOfDofs(int& nodl, int& nodl_scalar)
  {
  }

  
  template<class Dimension>
  inline void VarTransmission_Base<Dimension>::UpdateGlobalDofs()
  {
  }
  
  template<class Dimension>
  inline bool VarTransmission_Base<Dimension>::IsSymmetricCondition() const
  {
    return true;
  }

  
  template<class Dimension>
  inline const typename Dimension::MatrixN_N& VarTransmission_Base<Dimension>
  ::GetDfjMatrixQuadraturePoint(int i, int n, int k) const
  {
    return MatDfjNodes(i)(n)(k); 
  }
  
  
  template<class Dimension>
  inline const IVect& VarTransmission_Base<Dimension>::GetDofNumberOmegaMinus() const
  {
    return DdlVolMinus; 
  }
  
  
  template<class Dimension>
  inline const IVect& VarTransmission_Base<Dimension>::GetDofNumberOmegaPlus() const
  {
    return DdlVolPlus;
  }

  
  template<class Dimension> 
  inline int VarTransmission_Base<Dimension>::GetNbVolumeDofOnSurface() const
  {
    return nb_ddl_volume_on_surface;
  }


  template<class Dimension> 
  inline int VarTransmission_Base<Dimension>::GetNbInterfaces() const
  {
    return transmission_references.GetM();
  }
    
  
  template<class Dimension> 
  inline void VarTransmission_Base<Dimension>::SetReferences(const IVect& ref)
  {
    transmission_references = ref;
  }
  
  
  template<class Dimension> 
  inline const Mesh<Dimension2>& VarTransmission_Base<Dimension>::GetMeshPlate() const
  {
    return mesh_plate; 
  }
  
  
  template<class Dimension> 
  inline void VarTransmission_Base<Dimension>::SetMeshPlate(const Mesh<Dimension2>& mesh_p)
  {
    gamma_is_plate = true;
    mesh_plate = mesh_p;
  }
  
  
  template<class Dimension> 
  inline void VarTransmission_Base<Dimension>
  ::SetMeshPlateNumbering(const MeshNumbering<Dimension2>& mesh_n)
  {
    mesh_plate_num = mesh_n;
  }
  
  
  template<class Dimension> 
  inline void VarTransmission_Base<Dimension>::DuplicateExtremities(bool ext)
  {
    duplicate_extremity = ext;
  }

  
  template<class Dimension> 
  inline Vector<IVect>& VarTransmission_Base<Dimension>::GetNodleAcous()
  {
    return NodleAcous;
  }
  
  
  template<class Dimension> 
  inline IVect& VarTransmission_Base<Dimension>::GetDirichletAcous()
  {
    return DirichletAcous; 
  }
  
  
  template<class Dimension> 
  inline IVect& VarTransmission_Base<Dimension>::GetProcVolMinus()
  {
    return ProcVolMinus;
  }
  
  
  template<class Dimension> 
  inline IVect& VarTransmission_Base<Dimension>::GetProcVolPlus()
  {
    return ProcVolPlus;
  }
  
  
  template<class Dimension> 
  inline Vector<VectReal_wp>& VarTransmission_Base<Dimension>::GetJacobianAcous()
  {
    return JacobianAcous;
  }
  
  
  template<class Dimension> 
  inline Vector<IVect>& VarTransmission_Base<Dimension>::GetNodleVolMinus()
  {
    return NodleVolMinus;
  }
  
  
  template<class Dimension> 
  inline Vector<IVect>& VarTransmission_Base<Dimension>::GetNodleVolPlus()
  {
    return NodleVolPlus;
  }

  
  template<class Dimension> 
  inline void VarTransmission_Base<Dimension>
  ::AddTransmissionTerms(const Real_wp& alpha,
			 const GlobalGenericMatrix<Real_wp>& nat_mat,
			 VirtualMatrix<Real_wp>& mat_sp, int, int) 
  {
  }
  
  
  template<class Dimension> 
  inline void VarTransmission_Base<Dimension>
  ::AddTransmissionTerms(const Complex_wp& alpha,
			 const GlobalGenericMatrix<Complex_wp>& nat_mat,
			 VirtualMatrix<Complex_wp>& mat_sp, int, int) 
  {
  }
  
  
  template<class TypeEquation>
  inline VarTransmission<TypeEquation>
  ::VarTransmission(EllipticProblem<TypeEquation>& var)
    : VarTransmission_Base<typename TypeEquation::Dimension>(var) 
  {
  }
  
}

#define MONTJOIE_FILE_TRANSMISSION_MODEL_INLINE_CXX
#endif
