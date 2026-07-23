#ifndef MONTJOIE_FILE_PARAMETER_OUTPUT_HXX

namespace Montjoie
{
  //! base class to specify file ouputs, begin time and final time for output, etc
  class ParamOutputClass
  {
  protected :
    Real_wp t_begin; //!< initial time 
    Real_wp t_end; //!< final time
    mutable Real_wp t_current; //!< current time
    Real_wp deltat; //!< time step between two snapshots
    Real_wp epsilon_time; //!< threshold to distinguish times
    
    int nature; //!< type of output
    int component; //!< which components to write ?
    mutable int num_instantane; //!< snapshot number
    int dimension; //!< dimension of the output
    
    // output files
    VectString name_file;
    
    // is a sismo on a point ?
    bool sismo_point; //!< seismogramm ?
    
    // result vectors
    mutable Vector<Vector<Real_wp> > real_result;
    mutable Vector<Vector<Complex_wp> > cplx_result;
    WriteOnTheGoWithBuffer<Real_wp> output_sismo, output_sismo_diff;
    int size_buffer_sismo;
    bool first_output_sismo;
    //! if true, sismos are written in binary
    bool seismogramm_to_write_in_binary;

 
  public :
    ParamOutputClass();
    
    int GetComponent() const;
    int GetNature() const;
    void SetNature(int);
    void SetComponent(int);
    int IncrementSnapshot() const;
    int GetSnapshotNumber() const;
    bool GradientToCompute(int Nmax) const;
    bool IsSeismogramOutput() const;
    
    void SetInterval(const Real_wp& t0, const Real_wp& tf, const Real_wp& dt,
                     Real_wp eps = -1);
    
    void SetBinaryOutput(bool b = true);
    void SetBufferSize(int s);

    static int GetTrueType(int type, int dim);
    static int GetFalseType(int type, int dim, int dim_N);
    
    const string& GetTotalFieldFile() const;
    const string& GetDiffractedFieldFile() const;
    const string& GetFileName(int i) const;

    int GetNbFile() const;
    void ResizeNbFile(int);
    void AdimTime(const Real_wp& t0);
    
    const Real_wp& GetCurrentTime() const;
    void SetThresholdTime(const Real_wp& eps);
    
    void SetTotalFieldFile(const string&);
    void SetDiffractedFieldFile(const string&);
    void SetFileName(int i, const string&);

    size_t GetMemorySize() const;    

    void SetInputData(const VectString& param, bool sismo = false);
    
    bool SnapshotToStore(const Real_wp& t) const;

    void InitTime(const Real_wp& t0);
    void ChangeTime(int n, const Real_wp& t);
    
    void InitResult(const Vector<Real_wp>& x, int n) const;
    void InitResult(const Vector<Complex_wp>& x, int n) const;

    void UpdateResult(int num_file, int m, Vector<Real_wp>& x) const;
    void UpdateResult(int num_file, int m, Vector<Complex_wp>& x) const;
    
    void InitSismo(const string& DOSSIER_output, const string& name_file,
		   const string& name_file_diff, int add_total_field, bool double_prec);

    void AddVectorSismo(VectReal_wp& val, VectReal_wp& val_diff);
    void AddVectorSismo(VectComplex_wp& val, VectComplex_wp& val_diff);
    
    friend ostream& operator <<(ostream& out, const ParamOutputClass& p);
    
  };
      
} // namespace Montjoie

#define MONTJOIE_FILE_PARAMETER_OUTPUT_HXX
#endif

