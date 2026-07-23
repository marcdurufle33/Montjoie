#ifndef MONTJOIE_FILE_VAR_MIGRATION_HXX

namespace Montjoie
{
  
  template<class Dimension>
  class VarMigration_Base : public InputDataProblem_Base
  {
  public:
    // data file for the experience and the simulation
    string input_file_exp, input_file_simu;
    
  protected:
    // references where the data is measured and reversed
    Vector<int> ref_measure_exp;
    Vector<int> ref_measure_simu;

    bool write_intermediate_output;
    
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    
    // description of the sources
    Vector<Vector<string> > param_source;
    
    // mesh used to interpolate the solution on quadrature points of the referenced surfaces
    MeshInterpolationFEM<Dimension> mesh_interp;

    // weights for measurements
    VectReal_wp weights_measure;
    VectR_N points_measure;
    int nb_points_global_measure;
    
#ifdef SELDON_WITH_MPI
    // arrays for MPI communications
    Vector<IVect> num_points_to_recv, num_points_to_send;
#endif
    
    // operator for interpolating the experience solution at the quadrature points
    Matrix<Real_wp, General, ArrayRowSparse> proj_exp, proj_simu;
    
    // parameters for time impulsion
    int nb_points_time;
    Real_wp tau_impulse, Tmax, Tsinus;
    
  public:
    bool axisymmetric_problem;

    VarMigration_Base();

    int GetNbPointsQuadrature() const;
    void GetWeightsQuadrature(VectReal_wp& weights) const;
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    int RemoveDuplicateLocalization(VarProblem<Dimension>& var, IVect& NumElement, IVect& ProcLocal);
    
    // the quadrature points are constructed and localized for the experimental data
    void ComputeQuadrature(VarProblem<Dimension>& var,
			   VarProblem<Dimension>& var_exp);

    template<class T>
    void ProjectExperimentData(VarHarmonic_Base<T, Dimension>& var_exp,
			       const Matrix<T, General, ColMajor>& rhs,
			       Matrix<T, General, ColMajor>& data_exp);
    
    template<class T>
    void ProjectSimulationData(VarHarmonic_Base<T, Dimension>& var_simu,
			       const Matrix<T, General, ColMajor>& rhs,
			       Matrix<T, General, ColMajor>& data_simu);
    
    void ComputeExperiment(VarHarmonic_Base<Complex_wp, Dimension>& var_exp,
			   All_LinearSolver* solver_exp, Matrix<Complex_wp, General, ColMajor>& data_exp);

    void ComputeRightHandSide(VarHarmonic_Base<Complex_wp, Dimension>& var_simu,
			      Matrix<Complex_wp, General, ColMajor>& rhs);

    void ComputeSimulation(VarHarmonic_Base<Complex_wp, Dimension>& var_simu,
			   All_LinearSolver* solver_simu,
			   Matrix<Complex_wp, General, ColMajor>& data_simu,
			   Matrix<Complex_wp, General, ColMajor>& rhs);

    void ComputeIntegralResidu(VarProblem<Dimension>& var_simu,
			       Matrix<Complex_wp, General, ColMajor>& residu,
			       Matrix<Complex_wp, General, ColMajor>& sol_residu);
    
    void ComputeAdjoint(VarProblem<Dimension>& var_simu, All_LinearSolver* solver_simu,
			Matrix<Complex_wp, General, ColMajor>& residu,
			Matrix<Complex_wp, General, ColMajor>& sol_residu);
    
  };

  
  template<class TypeEquation>
  class VarMigration : public VarMigration_Base<typename TypeEquation::Dimension>
  {
  };

}

#define MONTJOIE_FILE_VAR_MIGRATION_HXX
#endif

