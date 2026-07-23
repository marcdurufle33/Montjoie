#ifndef MONTJOIE_FILE_MODEL_WIRES_HXX
 
namespace Montjoie
{

  ///////////
  // WIRES //
  ///////////
  class VarModelWires
  {
  protected :
    Vector<Wire<Dimension> > Tab_Fil; //!< array of wires
    IVect num_ddl_fil; //!< dof numbers for each wire
    IVect num_ddl_fctfil; //!< dof numbers for additional function, associated with wires
    string file_courant_fil;
    //!< file name for output of Lagrange mutliplier associated with the wire
    string file_courant_eff;
    //!< file name for output of additional function associated with the wire
    bool holland_wire; //!< use of Holland model for wires ?
    bool lagrange_wire; //!< ...
    bool holland_delta_app_dirac; //!< should be false
    
    // Attribut utilisés uniquement pour des tests sur les fils
    int DofCentral;  //!< internal variable
    int NumDofAdj; //!< internal variable
    R_N NoeudCentral; //!< internal variable
    VectR_N NoeudAdj; //!< internal variable
    VectComplexe x_test; //!< internal variable
    VectComplexe b_test; //!< internal variable
    
#ifdef MONTJOIE_WITH_TWO_DIM 
    void IncrementNumberDofs_Wires(Dimension2 & dim);
    template<class MatrixSparse> 
    void TreatWires(MatrixSparse& mat_sp, Dimension2& dim);
    //! to be done
    void AddComponent_AdditionalFunction(const VectComplexe& U0, int num_fil,
                                         const R2& pt_glob, VectComplexe& val_u) const {}
#endif

#ifdef MONTJOIE_WITH_THREE_DIM 
    void IncrementNumberDofs_Wires(Dimension3 & dim);
    template<class MatrixSparse> 
    void TreatWires(MatrixSparse& mat_sp, Dimension3& dim);
    //! scalar case, no scalar product
    template<class Vector1, class Vector2>
    Real_wp DotProd_WireDirection(const Vector1& val, const Vector2& axis) { return val(0); }
    template<class MatrixSparse>
    void TreatWire_Holland(MatrixSparse& mat_sp);
    template<class VectMatrix>
    void ComputeHollandOperator(VectMatrix& Glob_mat_Bholland);
    template<class VectMatrix>
    void ComputeOperator_FictitiousDomain(VectMatrix& Cfict, VectMatrix& Badd_stiff, 
					  VectMatrix& Badd_mass, VectMatrix& beta_stiff,
                                          VectMatrix& beta_mass);
    
    void AddComponent_AdditionalFunction(const VectComplexe& U0, int num_fil,
                                         const R3& pt_glob, VectComplexe& val_u) const;
#endif

    // calcul pour les fils
    static const int FCT_WIRE_AUCUNE = 0; //!< no additional function
    static const int FCT_WIRE_GAUSSIENNE = 1; //!< gaussian cut-off
    static const int FCT_WIRE_PLATEAU_H2 = 2; //!< H^2 cut-off
    static const int FCT_WIRE_PLATEAU_H3 = 3 ; //!< H^3 cut-off

#ifdef MONTJOIE_WITH_TWO_DIM
    template<class TypeSource>
    void SetWireSource(VectComplexe& b_source, TypeSource* f, Dimension2& dim, Complex_wp a);
    //! overloaded
    template<class TypeSource>
    void SetWireSource(VectComplexe& b_source, TypeSource* f, Dimension2& dim, Real_wp a) {}
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
    template<class TypeSource>
    void SetWireSource(VectComplexe& b_source, TypeSource* f, Dimension3& dim, Complex_wp a);
    //! overloaded
    template<class TypeSource>
    void SetWireSource(VectComplexe& b_source, TypeSource* f, Dimension3& dim, Real_wp a) {}
#endif
    
    template<class TypeSource>
    void SetPonctualSource(VectComplexe& b_source, TypeSource* f);    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_MODEL_WIRES_HXX
#endif
