#ifndef MONTJOIE_FILE_PARAMETER_OUTPUT_REPRISE_HXX

namespace Montjoie
{
  class ParameterOutputReprise
  {
  protected :
    ofstream file_out;
    ifstream file_in;
    int master_proc, rank_proc;
    int nb_procs_slave;
    
  public :    
    ParameterOutputReprise();
    
    void RegroupWritingOnMaster(int nb_proc_per_node);
    
    void StartReading(const string& input_file);
        
    void Read(VectReal_wp& x, bool master_only = false);
    void Read(Vector<int>& x, bool master_only = false);
    
    template<int p>
    void Read(Vector<TinyVector<Real_wp, p> >& x,
              bool master_only = false);
    
    void Read(Matrix<Real_wp>& A, bool master_only = false);    
    void Read(int& n, bool master_only = false);
    
    void FinalizeReading();
    
    void StartWriting(const string& output_file);

    void Write(const VectReal_wp& x, bool master_only = false);
    void Write(const Vector<int>& x, bool master_only = false);
    
    template<int p>
    void Write(const Vector<TinyVector<Real_wp, p> >& x,
               bool master_only = false);
    
    void Write(const Matrix<Real_wp>& A, bool master_only = false);
    void Write(int n, bool master_only = false);
    
    void FinalizeWriting();
    
  };
  
}

#define MONTJOIE_FILE_PARAMETER_OUTPUT_REPRISE_HXX
#endif
