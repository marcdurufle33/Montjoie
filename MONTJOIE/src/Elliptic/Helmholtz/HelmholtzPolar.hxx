#ifndef MONTJOIE_FILE_HELMHOLTZ_POLAR_HXX

namespace Montjoie
{
  class HelmholtzEquationPolar : public HelmholtzEquation1D
  {
  };

  
  //! base class for solving Helmholtz equation in polar coordinates
  template<class Complexe>
  class VarHelmholtz_Polar : public VarHelmholtz_1D<Complexe>
  {
  public :
    typedef Vector<Complexe> VectComplexe; //!< "complex" vector
    
    //! maximal mode number (l)
    int Lmax;
    //! threshold used to stop the computation
    Real_wp threshold_mode;
    //! if true the number of involved modes is automatically computed
    bool number_mode_to_be_computed;
    Vector<ParamOutputClass> output_grid_param; //!< output parameters
    //! predefined display grids
    Vector<GridInterpolationFull<Dimension1> > var_grid;
    //! diagonal matrix for \int mu \varphi_j \varphi_i
    Vector<Complexe> mass_matrix;
    //! matrix for -\omega^2 \int \rho r^2 \varphi_j \varphi_i
    //! + \int \mu r^2 d/dr \varphi_j d/dr \varphi_i
    Matrix<Complex_wp, General, BandedCol> stiffness_matrix;
    Matrix<Complex_wp, General, ArrayRowSparse> stiffness_matrix_sp;
    //! type of right hand side
    int type_source;
    enum{SRC_DIFFRACTED_FIELD, SRC_TOTAL_FIELD, SRC_DIRAC, SRC_GAUSSIAN, SRC_GAUSSIAN_R};
    //! direction of plane wave
    R2 wave_vector;
    //! origin of Dirac
    R2 origin_dirac; Real_wp radius_gaussian_src, radius_gaussian_cut_off;
    //! quadrature rules used to perform an integration in theta-coordinate
    Globatto<Real_wp> gauss_teta;
    //! impedance for transparent condition
    Vector<Complexe> coef_impedance_transparent;
    
    bool display_du_dr;
    
    VarHelmholtz_Polar();
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void ConstructAll(const string&, const string&, MontjoieTimer&);
    void ComputeVarGrid();
    
    void PerformFactorizationStep(Matrix<Complex_wp, General, BandedCol>& mat_lu,
				  All_MatrixLU<Complex_wp>& mat_lu_sp,
                                  IVect& pivot, int l);
    
    void ComputeSolution(Matrix<Complex_wp, General, BandedCol>& mat_lu,
			 All_MatrixLU<Complex_wp>& mat_lu_sp,
                         IVect& pivot, Vector<Complex_wp>& x);
    
    void ComputeRightHandSide(Vector<Vector<Complexe> >& rhs, MontjoieTimer&);
    
    void WriteDatas(Vector<Vector<Complexe> >& sol);
    
    // Inline methods
    int GetLmax() const;
    
    int GetIndexOutputFiles(int type);
    
  };
  
  
  template<>
  class EllipticProblem<HelmholtzEquationPolar> : public VarHelmholtz_Polar<Complex_wp>
  {
  };
  
}

#define MONTJOIE_FILE_HELMHOLTZ_POLAR_HXX
#endif
