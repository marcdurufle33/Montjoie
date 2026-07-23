#ifndef MONTJOIE_FILE_BRIDGE_HXX

namespace Montjoie
{
  
  class Bridge
  {
  protected :

    int nb_lagrange_mult;
    Real_wp scal_a, scal_b, scal_c, scal_d;
    Real_wp coef_couplage_V, coef_couplage_H, coef_couplage_any;
    string inertia_file_name, damp_file_name, stiff_file_name;
    
    Vector<Real_wp> vect_scal_a;
    Vector<Real_wp> vect_scal_b;	
    Vector<Real_wp> vect_couplage_plaque;
    Vector<Real_wp> sum_LaMu;
    
    int nev;

  public :
    
    int offset_bridge, offset_extra_bridge;
	  int size_bridge;
	  Real_wp phi_control, height_of_bridge;
    Real_wp vect_inertia, vect_stiff, vect_damp;
    bool extra_bridge, matchNLaMuToStringDOF;
    
    Bridge();

    void SetInputData(const string & keyword, const Vector<string> &param);
    template <class TypeEquation>
    void ConstructAll(MultiString_Base& multi, TypeEquation &var_eq, const string&);
    void ReadBridgeImpedance(const string& DOSSIER);    
    int GetSize();
    int GetNbLagrangeMultipliers() const;
    void SetNbLagrangeMultipliers(int n);

    void AddScheme(MultiString_Base& multi,
                   const VectReal_wp& somme_corde_0, const VectReal_wp& somme_corde_2,
                   Matrix<Real_wp>&, VectReal_wp& scheme);
    
    void AddRightHandSide(MultiString_Base& multi,
                          const VectReal_wp& somme_corde_0,
                          Matrix<Real_wp>& val2, VectReal_wp& scheme);
    
    void AddDiff(MultiString_Base& multi, VirtualMatrix<Real_wp>& DiffMatrix,
                  Matrix<Real_wp>& val2);
    
    Real_wp Compute_coef_couplage_V(MultiString_Base& multi, Matrix<Real_wp>& val2);
    Real_wp Compute_coef_couplage_H(MultiString_Base& multi, Matrix<Real_wp>& val2);
    Real_wp Compute_coef_couplage_any(int g2, MultiString_Base& multi, Matrix<Real_wp>& val2);
    Real_wp Compute_coef_couplage_Theta(MultiString_Base& multi, Matrix<Real_wp>& val2);
    Real_wp ComputeEnergy(MultiString_Base & multi, Matrix<Real_wp>& val2);
    
    inline void SetVecScal_a(int i, Real_wp scal_a_)
    {
      vect_scal_a(i) = scal_a_;
    }
    inline void SetVecScal_b(int i, Real_wp scal_b_)
    {
      vect_scal_b(i) = scal_b_;
    }    
    
            
    inline void SetScal_a(Real_wp scal_a_)
    {
      scal_a = scal_a_;
    }
    inline void SetScal_b(Real_wp scal_b_)
    {
      scal_b = scal_b_;
    }
    inline Real_wp GetScal_a()
    {
      return(scal_a);
    }
    inline Real_wp GetScal_b()
    {
      return(scal_b);
    }

    inline void SetScal_c(Real_wp scal_c_)
    {
      scal_c = scal_c_;
    }
    inline void SetScal_d(Real_wp scal_d_)
    {
      scal_d = scal_d_;
    }
    inline Real_wp GetScal_c()
    {
      return(scal_c);
    }
    inline Real_wp GetScal_d()
    {
      return(scal_d);
    }

  };


}

#define MONTJOIE_FILE_BRIDGE_HXX
#endif

