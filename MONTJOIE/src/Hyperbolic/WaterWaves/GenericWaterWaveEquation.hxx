#ifndef MONTJOIE_FILE_GENERIC_WATER_WAVE_EQUATION_HXX

namespace Montjoie
{
    
  class GenericWaterWaveEquation
  {
  public :
    double a; // amplitude
    double h0; // height
    double epsilon; // epsilon = a/h0
    double mu; 
    double lambda;
    double beta;
    double alpha, slope_bottom, yintercept;
    double c1, nu; // for soliton solution
    double center_init, freq_init;
    int type_bottom; double freq_bottom;
    enum {FLAT, SINUS, LINEAR, HALF_SINUS, PERIODIC_SINUS};
    int type_initial;
    enum {GAUSSIAN, SOLITON, CRENEAU, INIT_SINUS, BREAKING, SOLITON_GN};
    Globatto<Real_wp> lob;
    Matrix<Real_wp> ValPhi, GradPhi;
    int type_model;
    enum {GENTLE, STRONG, ORIGINAL};
    
    GenericWaterWaveEquation();
    
    template<class TypeEq>
    void Copy(const TypeEq& var);
    
    void InitEpsilon(double eps, double, double, double);
    void SetInitialCondition(int, int, double, double);
    void ComputeBottom(double, double, const VectReal_wp& x, VectReal_wp& b, VectReal_wp& c);
    
    void ComputeInitialCondition(double, double, const VectReal_wp& x, VectReal_wp& z0);
    void ComputeInitialCondition(double, double, const VectReal_wp& x, VectReal_wp& z0,
                                 VectReal_wp& z0_x, VectReal_wp& z0_xx);
    
    template<class GenericPb>
    void InitOrder(int order, const GenericPb& var);
    
    template<class GenericPb>
    void GetGradient(const GenericPb& var, const VectReal_wp& Un, VectReal_wp& dUn);
    
    template<class GenericPb>
    void GetUquadrature(const GenericPb& var, const VectReal_wp& Un, VectReal_wp& UnQuad);
    
  };


  class InputVariablesWaterWaves : public InputDataProblem_Base
  {
  public :
    double frequency, center_source;
    int type_initial_condition, type_model, type_bottom;
    double x0, xN, xg_0, xg_N;
    int nb_elt, order, nb_points_fdtd;
    double dt;
    double Tfinal;
    int number_grid_points; // number of points where we want to know the solution
    int delta_snapshot; 

    double log_eps_min, log_eps_max;
    int nb_points_log_eps;
    bool log_log_curve;
    
    int time_order;
    int type_time_scheme;
    
    double mu, cte_epsilon, power_epsilon, cte_alpha, power_alpha, cte_beta, power_beta;
    int n1_order, n2_order;

    // formulation used to solve the equation
    int type_formulation;
    // available formulations
    // SPECTRAL : Spectral Method
    // FDTD : Finite-Fifference Method
    // LDG : Local Discontinuous Galerkin
    // FEM : Finite Element Method
    enum {SPECTRAL, FDTD, LDG, FEM};
    
    string file_output_history;
    
    InputVariablesWaterWaves();
    
    void SetInputData(const string& keyword,
		      const Vector<string>& param);
  };
  
}

#define MONTJOIE_FILE_GENERIC_WATER_WAVE_EQUATION_HXX
#endif
