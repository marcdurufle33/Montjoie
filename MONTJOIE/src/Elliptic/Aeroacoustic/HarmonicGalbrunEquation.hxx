#ifndef MONTJOIE_FILE_HARMONIC_GALBRUN_EQUATION_HXX

namespace Montjoie
{

  template<class TypeEquation>
  class VarGalbrun_Eq : public VarHarmonic<TypeEquation>,
			public VarGalbrun_Dim<typename TypeEquation::Complexe,
					      typename TypeEquation::Dimension>
  {
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeEquation::Dimension Dimension;
    
  public:
    VarGalbrun_Eq();
   
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void InitIndices(int n);
    int GetNbPhysicalIndices() const;
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
			   IVect& num_ref, IVect& num_index, IVect& num_component,
			   Vector<bool>& compute_grad, Vector<bool>& compute_hess);

    void ComputeMeshAndFiniteElement(const string & name_elt, bool split_mesh = true);    
    void ComputePhysicalCoefficients();

    void PerformOtherInitializations();

    template<class T, class Dim>
    void ModifyVolumetricSourceGen(int i, int j, const typename Dim::R_N&,
                                   const VirtualSourceField<T, Dim>&,
                                   Vector<T>&) const;
    
    void ModifyVolumetricSourceGen(int i, int j, const typename Dimension::R_N&,
                                   const VirtualSourceField<Complexe, Dimension>&,
                                   Vector<Complexe>&) const;
    
    void ModifyVolumetricSource(int i, int j, const R2&,
                                const VirtualSourceField<Real_wp, Dimension2>&,
                                Vector<Real_wp>&) const;
    
    void ModifyVolumetricSource(int i, int j, const R2&,
                                const VirtualSourceField<Complex_wp, Dimension2>&,
                                Vector<Complex_wp>&) const;

    void ModifyVolumetricSource(int i, int j, const R3&,
                                const VirtualSourceField<Real_wp, Dimension3>&,
                                Vector<Real_wp>&) const;
    
    void ModifyVolumetricSource(int i, int j, const R3&,
                                const VirtualSourceField<Complex_wp, Dimension3>&,
                                Vector<Complex_wp>&) const;

  };
    
}

#define MONTJOIE_FILE_HARMONIC_GALBRUN_EQUATION_HXX
#endif
