#ifndef MONTJOIE_FILE_MULTI_FREQUENCY_PROBLEM_HXX

namespace Montjoie
{
  class MultiFrequencyParameter;

  //! main class for solving multi-frequency problems
  class MultiFrequencyProblem : public InputDataProblem_Base
  {
  protected:
    //! explored frequencies
    VectReal_wp omega;
    Real_wp wavelength_adim;
    
    //! number of other parameters
    int nb_glob_parameters;
    
    // how to solve the different frequencies
    int nb_parallel_proc_freq;
    bool same_pattern_LU, use_iterative_solver, do_not_store_facto;
    int nb_iter_max_eff_iterative_solver;
    
    //! if true, both parameters and frequencies are distributed
    //! if false, only frequencies are distributed
    bool distribute_param;
    
    // each processor treats frequencies between offset_omega and
    // offset_omega + omega.GetM()
    int offset_omega;    

    // each processor treats parameters between offset_param and
    // offset_param + nb_local_parameters
    int offset_param, nb_local_parameters;
    
#ifdef SELDON_WITH_MPI
    // comm_freq : global communicator between frequencies
    // comm_group_mode : communicator for each (frequency, parameter)
    // comm_param : communicator between parameters (for a given interval of frequencies)
    MPI_Comm comm_freq, comm_group_mode, comm_param, comm_first_param;
#endif
    
    IVect ref_impedance; string file_output_impedance;

  public:
    MultiFrequencyProblem();
    virtual inline ~MultiFrequencyProblem(){}
    
#ifdef SELDON_WITH_MPI
    inline const MPI_Comm& GetFrequencyCommunicator() const { return comm_freq; }
    inline const MPI_Comm& GetLocalCommunicator() const { return comm_group_mode; }
#endif

    inline void SetWaveLengthAdim(const Real_wp& L) { wavelength_adim = L; }
    void SetInputData(const string&, const Vector<string>& param);
    
    void DistributeFrequency(VectReal_wp& coef_input_frequency);
    void SolveAllFreq(VarComputationProblem& var_exp, All_LinearSolver*,
                      MultiFrequencyParameter& param);    

    template <class T>
    void GatherVectorFrequency(Vector<T>& y);
    
  };

  
  //! intermediate class specifying the parameters for each frequency
  class MultiFrequencyParameter
  {
  public:
    virtual void SetPulsation(int m, int mg, const Real_wp& omega,
                              VarComputationProblem&, MultiFrequencyProblem&);
    
    virtual void SetParameters(int k, int kg, int m, int mg,
                               VarComputationProblem&, MultiFrequencyProblem&);
    
    virtual void SetSourceFrequency(int k, int kg, int m, int mg,
                                    const VarComputationProblem&,
                                    Matrix<Complex_wp, General, ColMajor>& sol) = 0;
    
    virtual void GiveSolutionFrequency(int k, int kg, int m, int mg,
                                       const VarComputationProblem&,
                                       Matrix<Complex_wp, General, ColMajor>& sol_exp) = 0;
    
  };
  
}

#define MONTJOIE_FILE_MULTI_FREQUENCY_PROBLEM_HXX
#endif
