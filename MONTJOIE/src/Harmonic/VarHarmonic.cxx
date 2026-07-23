#ifndef MONTJOIE_FILE_VAR_HARMONIC_CXX

namespace Montjoie
{
  
  //! constructor of VarHarmonicClass, set default values
  template<class TypeEquation>
  VarHarmonic<TypeEquation>::VarHarmonic()
    : VarHarmonic_Base<Complexe, Dimension>(GetLeafClass()),
      fct_impedance_absorbing(GetLeafClass()), fct_impedance_generic(GetLeafClass())
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL    
    , fct_impedance_high_conduc(GetLeafClass())
#endif
    , output_rcs_param(this->GetLeafClass())
#ifdef MONTJOIE_WITH_TRANSMISSION
    , var_transmission(this->GetLeafClass())
#endif
    , var_gibc(this->GetLeafClass())
  {                    
  }

  
  template<class TypeEquation>
  void VarHarmonic<TypeEquation>::GetMemoryUsed(map<string, size_t>& var) const
  {
    VarHarmonic_Base<Complexe, Dimension>::GetMemoryUsed(var);
    size_t n = var_gibc.GetMemorySize();
    n += output_rcs_param.GetMemorySize();
#ifdef MONTJOIE_WITH_TRANSMISSION
    n += var_transmission.GetMemorySize();
#endif
    var["TransmissionData"] = n;
  }


  //! returns a new pointer with the description of the source
  template<class TypeEquation>
  VirtualSourceFEM<typename TypeEquation::Complexe, typename TypeEquation::Dimension>* 
  VarHarmonic<TypeEquation>::GetNewSourceEquationObject(int n)
  {
    const Vector<VectString>& param = this->source_space_param(n);
    switch (this->type_source(n))
      {
      case VarSourceProblem_Base::SRC_DIFFRACTED_FIELD:
	return new DiffractedWaveSource<TypeEquation>(this->GetLeafClass(),
						      *this->incident_wave_field(n));
      case VarSourceProblem_Base::SRC_TOTAL_FIELD:
	return new TotalWaveSource<TypeEquation>(this->GetLeafClass(),
						 *this->incident_wave_field(n));
      case VarSourceProblem_Base::SRC_USER:
	return new UserDefinedSource<Complexe, Dimension>(this->GetLeafClass(), param);
      case VarSourceProblem_Base::SRC_VOLUME:
	return new VolumetricSource<TypeEquation>(this->GetLeafClass(), param);
      case VarSourceProblem_Base::SRC_DIRAC:
	return new DiracSource<TypeEquation>(this->GetLeafClass(), param);
      }
    
    return NULL;
  }


  //! returns a new pointer with a Dirac Source
  template<class TypeEquation>
  VirtualSourceFEM<typename TypeEquation::Complexe, typename TypeEquation::Dimension>* 
  VarHarmonic<TypeEquation>::GetNewDiracSourceEquationObject(const Vector<VectString>& param) const
  {
    return new DiracSource<TypeEquation>(this->GetLeafClass(), param);
  }


  //! we add the mass matrix alpha*Mh
  template<class TypeEquation>
  void VarHarmonic<TypeEquation>::GetMassMatrix(VectComplexe& Dh, bool assemble) const
  {
    const EllipticProblem<TypeEquation>& var_leaf = this->GetLeafClass();
     GlobalGenericMatrix<Complexe> nat_mat;
    Complexe one, zero;
    SetComplexOne(one); SetComplexZero(zero);
    nat_mat.SetCoefStiffness(zero);
    nat_mat.SetCoefMass(one);
    nat_mat.SetCoefDamping(zero);
    Dh.Reallocate(var_leaf.offset_dof_unknown(this->nb_unknowns_scal));
    Dh.Zero();
    int nb_mesh_num = this->GetNbMeshNumberings();
    Vector<IVect> Nodle(nb_mesh_num);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	const ElementReference_Dim<Dimension>& Fb = this->GetReferenceElement(i);
	bool variable = this->UseNumericalIntegration(i);
	bool affine = this->mesh.IsElementAffine(i);
	int nb_dof_loc = Fb.GetNbDof();
	int ref = this->mesh.Element(i).GetReference();
	TinyMatrix<Complexe, General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> mass;
        for (int nm = 0; nm < nb_mesh_num; nm++)
          Nodle(nm) = this->GetMeshNumbering(nm).Element(i).GetNodle();
        
	if (variable)
	  {
	    for (int j = 0; j < nb_dof_loc; j++)
	      {
		TypeEquation::GetMassPhiDFiTensor(var_leaf, i, j, nat_mat, ref, mass,
						  variable, affine, Fb);
                
		for (int m = 0; m < TypeEquation::nb_unknowns; m++)
		  {
                    int nm = this->mesh_num_unknown(m);
                    Dh(Nodle(nm)(j) + this->offset_dof_unknown(m)) += mass(m, m);
                  }
	      }
	  }
	else
	  {
	    TypeEquation::GetMassPhiDFiTensor(var_leaf, i, 0, nat_mat, ref, mass,
					      variable, affine, Fb);
	    
	    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	      {
                int nm = this->mesh_num_unknown(m);
                for (int j = 0; j < nb_dof_loc; j++)
                  Dh(Nodle(nm)(j) + this->offset_dof_unknown(m)) += mass(m, m)*Fb.WeightsND(j);
              }
	  }
      }
    
    // adding terms of other domains (parallel implementation)
    if (assemble)
      if (this->FormulationDG() == ElementReference_Base::CONTINUOUS)
	this->AddDomains(Dh);
  }


  //! copying physical properties of another similar problem
  template<class TypeEquation> template<class TypeEq2>
  void VarHarmonic<TypeEquation>::
  CopyIndices(const EllipticProblem<TypeEq2>& var)
  {
    cout << "This method should be overloaded in derived class " << endl;
    abort();
  }
  
  
  //! computes arrays specific to the problem solved
  template<class TypeEquation>
  void VarHarmonic<TypeEquation>::ComputeArraySpecificEquation()
  {
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	TypeEquation::ComputeMassMatrix(this->GetLeafClass(), i,
					this->GetReferenceElement(i));
      }    
  }
  
}

#define MONTJOIE_FILE_VAR_HARMONIC_CXX
#endif
