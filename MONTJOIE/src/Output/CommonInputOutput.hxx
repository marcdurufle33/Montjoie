#ifndef MONTJOIE_FILE_COMMON_INPUT_OUTPUT_HXX

namespace Montjoie
{
  //! class grouping enumerations for outputs
  class OutputTypeEnum
  {
  public :
    
    // output on Matlab/Medit or Vtk
    enum { ASCII_MATLAB, BINARY_MATLAB, ASCII_VTK, BINARY_VTK};
    // output in single/double precision
    enum { SINGLE_PRECISION, DOUBLE_PRECISION, QUADRUPLE_PRECISION, MULTIPLE_PRECISION};
    // type of appends, NO : result is directly written, YES : result is added
    //                 FINAL : result is added and written on files
    enum {APPEND_NO, APPEND_YES, APPEND_FINAL};
  };
  
  
  void GetParallelDistributionPoints(int nb_proc, int rank_proc, int nbPointsRCS,
                                     int& nb_points_effective, int& num_angle0, int type_algo = 0);
  
  
  /***********************
   * Medit input/outputs *
   ************************/

  
  template<class T>
  void ReadMedit(Vector<complex<T> >& x_sol, const string& file_name);
  
  template<class T>
  void ReadMedit(Vector<T>& x_sol, const string& file_name);
  
  template<class T>
  void WriteMedit(const Vector<complex<T> >& x_sol,
                  const string& file_name, const Dimension2& dim,
		  bool double_prec = false);
  
  template<class T>
  void WriteMedit(const Vector<T>& x_sol,
                  const string& file_name, const Dimension2& dim,
		  bool double_prec = false);
  
  template<class T>
  void WriteMedit(const Vector<complex<T> >& x_sol,
                  const string& file_name, const Dimension3& dim,
		  bool double_prec = false);
  
  template<class T>
  void WriteMedit(const Vector<T>& x_sol,
                  const string& file_name, const Dimension3& dim,
		  bool double_prec = false);

  
  /*************************
   * Montjoie output files *
   *************************/

  
  template<class T>
  void ReadMeshData(Vector<Vector<T> >& nu, const string& file_name);
  
  template<class T>
  void ReadMeshData(Vector<Vector<complex<T> > >& nu, const string& file_name);
  
  template<class T>
  void WriteMeshData(const Vector<Vector<T> >& nu, const string& file_name);
  
  template<class T>
  void WriteMeshData(const Vector<Vector<complex<T> > >& nu,
                     const string& file_name);

  void CopyVector(const Vector<float>& x, VectReal_wp& y);
  void CopyVector(const Vector<double>& x, VectReal_wp& y);
  void CopyVector(const Vector<complex<float> >& x, VectReal_wp& y);
  void CopyVector(const Vector<complex<double> >& x, VectReal_wp& y);

  void CopyVector(const Vector<float>& x, VectComplex_wp& y);
  void CopyVector(const Vector<double>& x, VectComplex_wp& y);
  void CopyVector(const Vector<complex<float> >& x, VectComplex_wp& y);
  void CopyVector(const Vector<complex<double> >& x, VectComplex_wp& y);

  /*************
   * Vtk files *
   *************/
  

#ifdef MONTJOIE_WITH_MPFR
  mpreal swapEndian(mpreal d);
#endif

  double swapEndian(double d);
  float swapEndian(float d);
  int swapEndian(int d);
  
  class VtkDataTypeEnum
  {
  public :
    // type of datas in vtk files (we put only types handled by Montjoie)
    enum {SCALAR, VECTOR, TENSOR};
  };

  template<class T>
  void WriteVtk(const Vector<T>& x_sol, const string& name,
		ostream&, bool double_prec = false, bool ascii = false);
  
  template<class T>
  void WriteVtk(const Vector<complex<T> >& x_sol, const string& name,
		ostream&, bool double_prec = false, bool ascii = false);

  template<class T>
  void WriteVtk(const Vector<TinyVector<T, 3> >& x_sol, const string& name,
		ostream&, bool double_prec = false, bool ascii = false);
  
  template<class T>
  void WriteVtk(const Vector<TinyVector<complex<T>, 3> >& x_sol, const string& name,
		ostream&, bool double_prec = false, bool ascii = false);
  
  template<class T>
  void WriteVtk(const Vector<TinyMatrix<T, Symmetric, 3, 3> >& x_sol,
		const string& name,
		ostream&, bool double_prec = false, bool ascii = false);
		
  template<class T>
  void WriteVtk(const Vector<TinyMatrix<complex<T>, Symmetric, 3, 3> >& x_sol,
		const string& name,
		ostream&, bool double_prec = false, bool ascii = false);
  
  template<class T>
  void WriteVtk(const Vector<Vector<complex<T> > >& x_sol,
		const IVect& type, const Vector<IVect>& decomp,
		const Vector<string>& name,
		ostream&, bool double_prec = false, bool ascii = false);
  
  template<class T>
  void WriteVtk(const Vector<Vector<T> >& x_sol,
		const IVect& type,
		const Vector<IVect>& decomp, const Vector<string>& name,
		ostream&, bool double_prec = false, bool ascii = false);

  template<class T>
  void WriteComplexVectorVtk(const Vector<complex<T> > & output_vector,
			     ostream& file_out, bool double_prec = false,
			     bool with_size = false, bool swap = false);


  void WriteBinaryInteger(IVect&, ostream& file_out,
			  bool with_size = true, bool swap = false);

  template<class T, class Allocator>
  void WriteBinaryDoubleOrFloat(const Vector<T, VectFull, Allocator>& output_vector,
				ostream& file_out, bool double_prec = false,
				bool with_size = true, bool swap = false);
  
  template<class reel, class Allocator>
  void WriteBinaryDoubleOrFloat(const Vector<complex<reel>, VectFull,
				 Allocator>& output_vector,
				ostream& file_out, bool double_prec = false,
				bool with_size = true, bool swap = false);

  template<class real, int N>
  void WriteBinaryDoubleOrFloat(const TinyVector<real, N>& output_vector,
				ostream& file_out, bool double_prec = false,
				bool with_size = true, bool swap = false);

  template<class reel, int N>
  void WriteBinaryDoubleOrFloat(const TinyVector<complex<reel>, N>& output_vector,
				ostream& file_out, bool double_prec = false,
				bool with_size = true, bool swap = false);
  

  template<class T>
  void ReadVtk(Vector<complex<T> >& x_sol, string& name, istream&, bool ascii);
  
  template<class T>
  void ReadVtk(Vector<T>& x_sol, string& name, istream&, bool ascii);

  template<class T>
  void ReadVtk(Vector<TinyVector<complex<T>, 3> >& x_sol,
	       string& name, istream&, bool ascii);
  
  template<class T>
  void ReadVtk(Vector<TinyVector<T, 3> >& x_sol,
	       string& name, istream&, bool ascii);

  template<class T>
  void ReadVtk(Vector<TinyMatrix<complex<T>, Symmetric, 3, 3> >& x_sol,
	       string& name, istream&, bool ascii);
  
  template<class T>
  void ReadVtk(Vector<TinyMatrix<T, Symmetric, 3, 3> >& x_sol,
	       string& name, istream&, bool ascii);

  void ReadBinaryInteger(IVect& num, istream& file_in,
			 bool with_size = true, bool swap = false);
  
  template<class real, class Allocator>
  void ReadBinaryDoubleOrFloat(Vector<real, VectFull, Allocator>& output_vector,
			       istream& file_in, bool double_prec = false,
			       bool with_size = true, bool swap = false);
  
  template<class real, class Allocator>
  void ReadBinaryDoubleOrFloat(Vector<complex<real>, VectFull,
				      Allocator>& output_vector, istream& file_in,
			       bool double_prec = false,
			       bool with_size = true, bool swap = false);

  // functions to read input files
  istream& getlineM(istream& file_in, string& line);
  
#ifdef SELDON_WITH_MPI
  void ReadLinesFile(const string& file_name, Vector<string>& all_lines,
                     const MPI_Comm& comm = MPI_COMM_SELF);
#else
  void ReadLinesFile(const string& file_name, Vector<string>& all_lines);
#endif
  
  //! base class to derive in order to read a data file with the function ReadInputFile
  class InputDataProblem_Base
  {
  public :
    // in the derived classes overload the function SetInputData
    virtual void SetInputData(const string&, const Vector<string>&) = 0;
    
  };
  
  void ReadInputFile(const string& file_name, InputDataProblem_Base& vars);
  void ReadInputFile(const Vector<string>& all_lines, InputDataProblem_Base& vars);
  
  void EcritDossier(string & DOSSIER, string name_data_file, int num = -1);
  
  
  
  /***********************************
   * Writing on files with buffering *
   ***********************************/
  
  
  // inline functions
  void WriteOnTheGo(const string & data_file, Real_wp value);
  
  template<class T, int nb>
  void WriteOnTheGo(const string & data_file, TinyVector<T, nb> & value);
  
  template<class T>
  void WriteOnTheGo(const string & data_file, Vector<T> & value);
  
  //! class for writing seismogramms on files with delayed writing
  /*!
    The output file will consists of n columns :
    x0 x1 x2 ... x_n
    v0 v1 v2 ... v_n
    ...
    In order to initialise the output file (n must be known at compilation time) :
    
    WriteOnTheGoWithTinyBuffer<double, n> output;
    output.Init("fichier.dat", 100);
    
    And to add lines to write :
    
    TinyVector<double, n> line;
    output.AddTinyVect(line);
    
    Here the writing in the file will be performed each 100 lines
   */
  template<class T, int nb>
  class WriteOnTheGoWithTinyBuffer
  {
  protected:
    int size_buffer, iter;
    Vector< TinyVector<T,nb> > Buffer;
    string name_file;
    bool is_binary, is_double;
    
  public:
    WriteOnTheGoWithTinyBuffer();
    
    size_t GetMemorySize() const;
    
    void Init(const string &data_file, int size,
              bool remove_file = true, bool write_prec = true);
    
    void AddTinyVect(TinyVector<T,nb> &data);    
    void SetBinary(bool i_want_binary = true);
    void SetDoublePrecision(bool i_want_double = true);

    void WritePrecision();
    void CloseBuffer();
    
    ~WriteOnTheGoWithTinyBuffer();
    
  };
  
  
  //! class for writing seismogramms on files with delayed writing
  /*!
    The output file will consists of n columns :
    x0 x1 x2 ... x_n
    v0 v1 v2 ... v_n
    ...
    In order to initialise the output file :
    
    WriteOnTheGoWithTinyBuffer<double> output;
    output.Init("fichier.dat", 100);
    
    And to add lines to write :
    
    Vector<double> line(n);
    output.AddTinyVect(line);
    
    Here the writing in the file will be performed each 100 lines
   */
  template<class T>
  class WriteOnTheGoWithBuffer
  {
  protected:
    int size_buffer, iter;
    Vector< Vector<T> > Buffer;
    string name_file;
    bool is_binary, is_double;
    
  public:
    WriteOnTheGoWithBuffer();

    size_t GetMemorySize() const;
    
    void Init(const string &data_file, int size,
              bool remove_file = true, bool write_prec = true);
    
    void AddVect(Vector<T> &data);
    void AddVect(Vector<complex<T> > &data);
    void SetBinary(bool i_want_binary = true);
    void SetDoublePrecision(bool i_want_double = true);
    
    void WritePrecision();
    void CloseBuffer();
    
    ~WriteOnTheGoWithBuffer();
    
  };
    
}

#define MONTJOIE_FILE_COMMON_INPUT_OUTPUT_HXX
#endif
  
