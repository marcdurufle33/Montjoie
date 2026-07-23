#ifndef MONTJOIE_FILE_VAR_GEOMETRY_PROBLEM_INLINE_CXX

namespace Montjoie
{

  /**********************
   * VarGeometryProblem *
   **********************/
  
  
  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline VarGeometryProblem<Dimension>::VarGeometryProblem(EllipticProblem<TypeEquation>& var)
    : VarProblem_Base(var), mesh(), mesh_num(mesh), var_boundary(var), var_comm(var), var_output(var)
  {
    InitDefaultValues();
  }
  

  //! returns the wave vector k of the incident field, eg a plane wave
  /*!
    The wave vector is used for the definition of the incident field
    when the user wants to compute the scattering of an object
   */
  template<class Dimension>
  inline const typename Dimension::R_N& VarGeometryProblem<Dimension>::GetWaveVector() const
  {
    return this->kwave;
  }


  //! sets the wave vector k of the incident field, eg a plane wave
  /*!
    The wave vector is used for the definition of the incident field
    when the user wants to compute the scattering of an object
   */
  template<class Dimension>
  inline void VarGeometryProblem<Dimension> ::SetWaveVector(const typename Dimension::R_N& k)
  {
    this->kwave = k;
  }
  
  
  //! returns the origin of the incident field
  /*!
    For example, the plane wave will be equal to exp( ik (x-x0))
    x0 is the origin of the phase
   */
  template<class Dimension>
  inline const typename Dimension::R_N& VarGeometryProblem<Dimension>::GetPhaseOrigin() const
  {
    return this->origine_phase;
  }


  //! returns the polarization of the incident field or of the source
  template<class Dimension>
  inline const VectReal_wp& VarGeometryProblem<Dimension>::GetPolarization() const
  {
    return this->polarization;
  }


  //! returns the polarization of the incident field or of the source
  template<class Dimension>
  inline const VectReal_wp& VarGeometryProblem<Dimension>::GetPolarizationGrad() const
  {
    return this->polarization_grad;
  }


  //! returns the polarization of the incident field or of the source
  template<class Dimension> template<class T>
  inline void VarGeometryProblem<Dimension>::GetPolarization(Vector<T>& P) const
  {
    P.Reallocate(this->polarization.GetM());
    for (int i = 0; i < P.GetM(); i++)
      P(i) = this->polarization(i);    
  }


  //! returns the polarization for gradient sources
  template<class Dimension> template<class T>
  inline void VarGeometryProblem<Dimension>::GetPolarizationGrad(Vector<T>& P) const
  {
    P.Reallocate(this->polarization_grad.GetM());
    for (int i = 0; i < P.GetM(); i++)
      P(i) = this->polarization_grad(i);    
  }


  //! returns the polarization of the incident field or of the source
  template<class Dimension> template<class T, int p>
  inline void VarGeometryProblem<Dimension>::GetPolarization(TinyVector<T, p>& P) const
  {
    for (int i = 0; i < p; i++)
      P(i) = this->polarization(i);
  }


  //! sets the polarization of the incident field or of the source
  template<class Dimension> template<int p>
  inline void VarGeometryProblem<Dimension>::SetPolarization(const TinyVector<Real_wp, p>& P)
  {
    for (int k = 0; k < p; k++)
      this->polarization(k) = P(k);
  }
  

  //! method called once the physical indices have been computed
  template<class Dimension>
  inline void VarGeometryProblem<Dimension>::FinalizeComputationVaryingIndices()
  {
    // default equation : nothing is done
  }

  
  //! method called to allocate arrays needed to the computation of the elementary matrix
  template<class Dimension>
  inline void VarGeometryProblem<Dimension>::AllocateMassMatrices()
  {
  }


  template<class Dimension>
  inline void VarGeometryProblem<Dimension>::DoNotComputeGrid()
  {
    var_output.grid_to_be_computed = false;
  }

  
  //! returns the dimension
  template<class Dimension>
  inline int VarGeometryProblem<Dimension>::GetDimension() const
  {
    return Dimension::dim_N;
  }  

  
  //! retrieves components of u on the element i
  /*!
    \param[in] B components of u on all the degrees of freedom of the mesh
    \param[in] i element number
    \param[out] Eloc components of u on degrees of freedom of element i
    Somehow the method performs Eloc = B(num_dof)
    where num_dof would be the dof numbers of element i
   */
  template<class Dimension>
  template<class Vector1, class Vector2>
  inline void VarGeometryProblem<Dimension>
  ::GetLocalUnknownVector(const Vector1& B, int i, Vector2& Eloc) const
  {
    this->mesh_num.number_map.GetLocalUnknownVector(this->mesh_num, B, i, Eloc);
  }
  
  
  //! adds a local vector to the global vector
  /*!
    \param[in] alpha coefficient
    \param[in] Eloc local vector to add to C
    \param[in] i element number
    \param[inout] C components on all the degrees of freedom of the mesh    
    Somehow the method performs C(num_dof) = C(num_dof) + alpha Eloc
    where num_dof would be the dof numbers of element i
   */
  template<class Dimension>
  template<class T0, class T1, int nb_u, class Vector2>
  inline void VarGeometryProblem<Dimension>
  ::AddLocalUnknownVector(const T0& alpha, TinyVector<Vector<T1>, nb_u>& Eloc,
			  int i, Vector2& C) const
  {
    this->mesh_num.number_map.AddLocalUnknownVector(this->mesh_num, alpha, Eloc, i, C);
  }
  
  
  //! Modification of a local vector so that it can be directly added to the global one  
  template<class Dimension> template<class Vector1>
  inline void VarGeometryProblem<Dimension>
  ::ModifyLocalUnknownVector(Vector1& U, int iquad) const
  {
    this->mesh_num.number_map.ModifyLocalUnknownVector(this->mesh_num, U, iquad);
  }
  

  //! Modification of a local vector so that it can be directly added to the global one  
  template<class Dimension> template<class Vector1>
  inline void VarGeometryProblem<Dimension>
  ::ModifyLocalComponentVector(Vector1& U, int iquad) const
  {
    this->mesh_num.number_map.ModifyLocalComponentVector(this->mesh_num, U, iquad);
  }
  

  //! projection from a global vector to a local vector  
  template<class Dimension> template<class Vector1>
  inline void VarGeometryProblem<Dimension>::GetGlobalUnknownVector(Vector1& U, int iquad) const
  {
    this->mesh_num.number_map.GetGlobalUnknownVector(this->mesh_num, U, iquad);
  }


  //! returns true if element i is inside PML
  template<class Dimension>
  inline bool VarGeometryProblem<Dimension>::InsidePML(int i) const
  {
    return this->mesh.Element(i).IsPML();
  }


  //! returns true if element i is inside PML
  template<class Dimension>
  inline bool VarGeometryProblem<Dimension>::ElementInsidePML(int i) const
  {
    return this->mesh.Element(i).IsPML();
  }

  
  //! writes the mesh on a file
  template<class Dimension>
  inline void VarGeometryProblem<Dimension>::WriteMesh(const string& file_name)
  {
    this->mesh.Write(file_name);
  }
  
  
  //! sets periodic dofs with the same dof numbers as original dofs
  template<class Dimension>
  inline void VarGeometryProblem<Dimension>::SetSameNumberPeriodicDofs()
  {
    this->mesh_num.SetSameNumberPeriodicDofs();
    for (int i = 0; i < this->other_mesh_num.GetM(); i++)
      this->other_mesh_num(i)->SetSameNumberPeriodicDofs();
  }

}

#define MONTJOIE_FILE_VAR_GEOMETRY_PROBLEM_INLINE_CXX
#endif

