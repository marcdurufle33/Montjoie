#ifndef MONTJOIE_FILE_COMMON_INPUT_OUTPUT_CXX

#include "CommonInputOutput.hxx"

namespace Montjoie
{
  
  //! distributes to each processor a number of points to treat
  /*!
    \param[in] nbPointsRCS total number of points to treat
    \param[out] nb_points_effective for each processor the number of points to treat
    \param[out] num_angle0 for each processor the number of the first point to treat
    Assuming that you have 100 points to treat (e.g. computation of radar cross section
    for 100 different angle), each computation associated with a point is independent.
    You wish that the processor 0 treats point 0, 1, ..., 10, processor 1 treats 
    11, 12, ..., 20, processor 3, 21, 22, etc
    By calling this function, you obtain for each processor, the number of points 
    to treat and the first point number 
  */
  void GetParallelDistributionPoints(int nb_proc, int rank_proc, int nbPointsRCS,
                                     int& nb_points_effective, int& num_angle0, int type_algo)
  {
    if (nb_proc > 1)
      {
        int nb_points_div = nbPointsRCS/nb_proc;
        if (nbPointsRCS == nb_points_div*nb_proc)
          {
            // each processor has exactly the same number of points
            num_angle0 = nb_points_div*rank_proc;
            nb_points_effective = nb_points_div;
          }
        else
          {
	    if (type_algo == 0)
	      {
		// last processor has less points
		nb_points_div = toInteger(ceil(double(nbPointsRCS)/nb_proc));
		num_angle0 = nb_points_div*rank_proc;
		nb_points_effective = nb_points_div;
		if (rank_proc == (nb_proc-1))
		  nb_points_effective = nbPointsRCS - num_angle0;
	      }
	    else
	      {
		// better : the first p processors have n points
		//          and other processors n-1 points
		int nb_points_div = toInteger(ceil(double(nbPointsRCS)/nb_proc));
		int p = nbPointsRCS - (nb_points_div-1)*nb_proc;
		if (rank_proc < p)
		  {
		    nb_points_effective = nb_points_div;
		    num_angle0 = nb_points_div*rank_proc;
		  }
		else
		  {
		    nb_points_effective = nb_points_div-1;
		    num_angle0 = p*nb_points_div + (nb_points_div-1)*(rank_proc-p);
		  }
	      }
          }
      }
    else
      {
        nb_points_effective = nbPointsRCS;
        num_angle0 = 0;
      }
  }
  
  
  /***********************
   * Medit input/outputs *
   ************************/
  
  
  //! reading of a bb file (medit format)
  /*!
    \param[out] x_sol solution on the vertices of the mesh
    \param[in] file_name name of the bb file
    This functions extracts a solution from a bb file
    For complex numbers, it is assumed that there are two associated files
    toto_real.bb and toto_imag.bb containing the real and imaginary part
  */
  template<class T>
  void ReadMedit(Vector<complex<T> >& x_sol, const string& file_name)
  {
    string file_name_real, file_name_imag;
    file_name_real = file_name + string("_real.bb");
    file_name_imag = file_name + string("_imag.bb");
    ifstream file_out_real(file_name_real.data());
    ifstream file_out_imag(file_name_imag.data());
    int size_vec, tmp;
    if ((!file_out_real.is_open())||(!file_out_imag.is_open()))
      {
	cout << "Unable to read files" << endl;
	cout << file_name_real << endl;
	cout << file_name_imag << endl;
	abort();
      }
    
    file_out_real >> tmp >> tmp >> size_vec >> tmp;
    file_out_imag >> tmp >> tmp >> tmp >> tmp;
    x_sol.Reallocate(size_vec); Real_wp val_real, val_imag;
    for (int i = 0; i < x_sol.GetM(); i++)
      {
	file_out_real >> val_real; 
	file_out_imag >> val_imag;
	x_sol(i) = Complex_wp(val_real, val_imag);
      }
    
    file_out_real.close();
    file_out_imag.close();
  }
  
  
  //! reading of a bb file (medit format)
  /*!
    \param[out] x_sol solution on the vertices of the mesh
    \param[in] file_name name of the bb file
    This functions extracts a solution from a bb file
  */
  template<class T>
  void ReadMedit(Vector<T>& x_sol, const string& file_name)
  {
    ifstream file_out(file_name.data());
    int size_vec, tmp;
    if (!file_out.is_open())
      {
	cout << "Unable to read files " << endl;
	cout << file_name << endl;
	abort();
      }
    
    file_out>>tmp>>tmp>>size_vec>>tmp;
    x_sol.Reallocate(size_vec);
    for (int i = 0; i < x_sol.GetM(); i++)
      file_out >> x_sol(i); 
    
    file_out.close();
  }


  //! writing of a bb file (medit format)
  /*!
    \param[in] x_sol solution to write
    \param[in] file_name name of the bb file
    \param[in] dim Dimension of the associated mesh
    \param[in] double_prec true if the solution is written in double precision
    This function writes a solution in a bb file.
    For a complex solution, three bb files will be created,
    if file_name is "toto.bb", the files "toto_real.bb", "toto_imag.bb" and
    "toto_abs.bb" are created and will contain the real, imaginary part and modulus
    of the solution
  */
  template<class T>
  void WriteMedit(const Vector<complex<T> >& x_sol,
                  const string& file_name, const Dimension2& dim, bool double_prec)
  {
    string file_name_real, file_name_imag, file_name_abs;
    string root = GetBaseString(file_name);
    
    // symbolic links to the mesh file
    // and we write separately real, imaginary and modulus of the solution
    file_name_real = root + string("_real.bb");
    MakeLink(root+".mesh", root+"_real.mesh");
    file_name_imag = root + string("_imag.bb");
    MakeLink(root+".mesh", root+"_imag.mesh");
    file_name_abs = root+ string("_abs.bb");
    MakeLink(root+".mesh", root+"_abs.mesh");
    
    ofstream file_out_real(file_name_real.data());
    ofstream file_out_imag(file_name_imag.data());
    ofstream file_out_abs(file_name_abs.data()); 
    
    file_out_real<<"2 1 "<<x_sol.GetM()<<" 2"<<endl;
    file_out_imag<<"2 1 "<<x_sol.GetM()<<" 2"<<endl;
    file_out_abs<<"2 1 "<<x_sol.GetM()<<" 2"<<endl;
    
    int nb_digits = 8;
    if (double_prec)
      nb_digits = 16;
    
    file_out_real.precision(nb_digits);
    file_out_imag.precision(nb_digits);
    file_out_abs.precision(nb_digits);
    
    for (int i = 0; i < x_sol.GetM(); i++)
      {
	file_out_real<<real(x_sol(i))<<endl;
	file_out_imag<<imag(x_sol(i))<<endl;
	file_out_abs<<abs(x_sol(i))<<endl;
      }
    
    file_out_real.close(); file_out_imag.close();  file_out_abs.close();
  }
  
  
  //! writing of a bb file (medit format)
  /*!
    \param[in] x_sol solution to write
    \param[in] file_name name of the bb file
    \param[in] dim Dimension of the associated mesh
    \param[in] double_prec true if the solution is written in double precision
    This function writes a solution in a bb file.
  */
  template<class T>
  void WriteMedit(const Vector<T>& x_sol,
                  const string& file_name, const Dimension2& dim, bool double_prec)
  {
    string file_name_true;
    file_name_true = GetBaseString(file_name) + string(".bb");
    ofstream file_out(file_name_true.data());
    file_out<<"2 1 "<<x_sol.GetM()<<" 2"<<endl;
    if (double_prec)
      file_out.precision(16);
    else
      file_out.precision(8);
    
    for (int i=0 ; i<x_sol.GetM() ; i++)
      file_out<<x_sol(i)<<endl;
    
    file_out.close();
  }
  
  
  //! writing of a bb file (medit format)
  /*!
    \param[in] x_sol solution to write
    \param[in] file_name name of the bb file
    \param[in] dim Dimension of the associated mesh
    \param[in] double_prec true if the solution is written in double precision
    This function writes a solution in a bb file.
    For a complex solution, three bb files will be created,
    if file_name is "toto.bb", the files "toto_real.bb", "toto_imag.bb" and
    "toto_abs.bb" are created and will contain the real, imaginary part and modulus
    of the solution
  */
  template<class T>
  void WriteMedit(const Vector<complex<T> >& x_sol,
                  const string& file_name, const Dimension3 & dim, bool double_prec)
  {
    string file_name_real,file_name_imag, file_name_abs;
    string root = GetBaseString(file_name);

    // symbolic links to the mesh file
    // and we write separately real, imaginary and modulus of the solution
    file_name_real = root + string("_real.bb");
    MakeLink(root+".mesh", root+"_real.mesh");
    file_name_imag = root + string("_imag.bb");
    MakeLink(root+".mesh", root+"_imag.mesh");
    file_name_abs = root+ string("_abs.bb");
    MakeLink(root+".mesh", root+"_abs.mesh");

    ofstream file_out_real(file_name_real.data());
    ofstream file_out_imag(file_name_imag.data());
    ofstream file_out_abs(file_name_abs.data()); 
    file_out_real<<"3 1 "<<x_sol.GetM()<<" 2"<<endl;
    file_out_imag<<"3 1 "<<x_sol.GetM()<<" 2"<<endl;
    file_out_abs<<"3 1 "<<x_sol.GetM()<<" 2"<<endl;
    
    int nb_digits = 8;
    if (double_prec)
      nb_digits = 16;
    
    file_out_real.precision(nb_digits);
    file_out_imag.precision(nb_digits);
    file_out_abs.precision(nb_digits);
    
    for (int i = 0; i < x_sol.GetM(); i++)
      {
	file_out_real<<real(x_sol(i))<<endl;
	file_out_imag<<imag(x_sol(i))<<endl;
	file_out_abs<<abs(x_sol(i))<<endl;
      }
    file_out_real.close();
    file_out_imag.close();
    file_out_abs.close();
  }
  
  
  //! writing of a bb file (medit format)
  /*!
    \param[in] x_sol solution to write
    \param[in] file_name name of the bb file
    \param[in] dim Dimension of the associated mesh
    \param[in] double_prec true if the solution is written in double precision
    This function writes a solution in a bb file.
  */
  template<class T>
  void WriteMedit(const Vector<T>& x_sol,
                  const string& file_name, const Dimension3& dim, bool double_prec)
  {
    string file_name_true;
    if (file_name.find(".bb",0) == string::npos)
      file_name_true = file_name + string(".bb");
    else
      file_name_true = file_name;
    
    ofstream file_out(file_name_true.data());
    file_out<<"3 1 "<<x_sol.GetM()<<" 2"<<endl;
    if (double_prec)
      file_out.precision(16);
    else
      file_out.precision(8);
	
    for (int i = 0; i < x_sol.GetM(); i++)
      file_out<<x_sol(i)<<endl;
    
    file_out.close();
  }
  

  /*************************
   * Montjoie output files *
   *************************/

  
  //! reads a discontinuous index from a binary file
  /*!
    \param[in] nu discontinuous index
    \param[in] file_name name of the binary file
    This function extracts a discontinuous index that has been written
    by calling function WriteMeshData. 
  */
  template<class T>
  void ReadMeshData(Vector<Vector<T> >& nu, const string& file_name)
  {
    ifstream file_in(file_name.data());
    int nb_elt, type_data;
    file_in.read(reinterpret_cast<char*>(&type_data), sizeof(int));
    if (type_data != 0) 
      {
	cout << "Error while reading the file " << file_name << endl;
        cout << "Reading a complex index into a real one" << endl;
        cout << "Please provide a real index" << endl;
        abort();
      }
    
    file_in.read(reinterpret_cast<char*>(&nb_elt), sizeof(int));
    
    IVect offset(nb_elt+1);
    file_in.read(reinterpret_cast<char*>(offset.GetData()),
                 (nb_elt+1)*sizeof(int));
    
    nu.Reallocate(nb_elt);
    for (int i = 0; i < nb_elt; i++)
      {
        int nb_nodes = offset(i+1) - offset(i);
        nu(i).Reallocate(nb_nodes);
        file_in.read(reinterpret_cast<char*>(nu(i).GetData()),
                     nb_nodes*sizeof(T));
      }
    
    file_in.close();
  }
  

  //! reads a discontinuous index from a binary file
  /*!
    \param[in] nu discontinuous index
    \param[in] file_name name of the binary file
    This function extracts a discontinuous index that has been written
    by calling function WriteMeshData. 
  */  
  template<class T>
  void ReadMeshData(Vector<Vector<complex<T> > >& nu, const string& file_name)
  {
    ifstream file_in(file_name.data());
    int nb_elt, type_data;
    file_in.read(reinterpret_cast<char*>(&type_data), sizeof(int));
    if (type_data != 1) 
      {
	cout << "Error while reading the file " << file_name << endl;
        cout << "Reading a real index into a complex one" << endl;
        cout << "Please provide a complex index" << endl;
	cout << " in ReadMeshData" << endl;
        abort();
      }

    file_in.read(reinterpret_cast<char*>(&nb_elt), sizeof(int));
    IVect offset(nb_elt+1);
    file_in.read(reinterpret_cast<char*>(offset.GetData()),
                 (nb_elt+1)*sizeof(int));
    
    nu.Reallocate(nb_elt);
    for (int i = 0; i < nb_elt; i++)
      {
        int nb_nodes = offset(i+1) - offset(i);
        nu(i).Reallocate(nb_nodes);
        file_in.read(reinterpret_cast<char*>(nu(i).GetData()),
                     2*nb_nodes*sizeof(T));
      }
    
    file_in.close();
  }
  
  
  //! writes a discontinuous index in a binary file
  /*!
    \param[in] nu discontinuous index
    \param[in] file_name name of the binary file
    This function writes a discontinuous index file in a binary file.
    We assume that this index nu is given by its values on quadrature points
    (or nodal points) for each element of the mesh. Therefore nu is
    an array of array.
  */
  template<class T>
  void WriteMeshData(const Vector<Vector<T> >& nu, const string& file_name)
  {
    int nb_elt = nu.GetM();
    IVect offset(nb_elt + 1);
    offset(0) = 0;
    for (int i = 0; i < nb_elt; i++)
      offset(i+1) = offset(i) + nu(i).GetM();
    
    int type_data = 0; // real numbers    
    ofstream file_out(file_name.data());
    file_out.write(reinterpret_cast<char*>(&type_data), sizeof(int));
    file_out.write(reinterpret_cast<char*>(&nb_elt), sizeof(int));
    file_out.write(reinterpret_cast<char*>(offset.GetData()),
                   (nb_elt+1)*sizeof(int));
    
    for (int i = 0; i < nb_elt; i++)
      file_out.write(reinterpret_cast<char*>(nu(i).GetData()),
                     (offset(i+1) - offset(i))*sizeof(T));
    
    file_out.close();
  }
  

  //! writes a discontinuous index in a binary file
  /*!
    \param[in] nu discontinuous index
    \param[in] file_name name of the binary file
    This function writes a discontinuous index file in a binary file.
    We assume that this index nu is given by its values on quadrature points
    (or nodal points) for each element of the mesh. Therefore nu is
    an array of array.
  */  
  template<class T>
  void WriteMeshData(const Vector<Vector<complex<T> > >& nu,
                     const string& file_name)
  {
    int nb_elt = nu.GetM();
    IVect offset(nb_elt + 1);
    offset(0) = 0;
    for (int i = 0; i < nb_elt; i++)
      offset(i+1) = offset(i) + nu(i).GetM();
    
    int type_data = 1; // complex numbers    
    ofstream file_out(file_name.data());
    file_out.write(reinterpret_cast<char*>(&type_data), sizeof(int));
    file_out.write(reinterpret_cast<char*>(&nb_elt), sizeof(int));
    file_out.write(reinterpret_cast<char*>(offset.GetData()),
                   (nb_elt+1)*sizeof(int));
    
    for (int i = 0; i < nb_elt; i++)
      {
        file_out.write(reinterpret_cast<char*>(nu(i).GetData()),
                       2*(offset(i+1) - offset(i))*sizeof(T));
      }
    
    file_out.close();
  }
  
  
  //! conversion from float vector to double vector
  void CopyVector(const Vector<float>& x, VectReal_wp& y)
  {
    y.Reallocate(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      y(i) = x(i);
  }
  
  
  //! conversion from double vector to double vector
  void CopyVector(const Vector<double>& x, VectReal_wp& y)
  {
    y.Reallocate(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      y(i) = x(i);
  }
  
  
  //! conversion from complex float vector to double vector (we take the real part)
  void CopyVector(const Vector<complex<float> >& x, VectReal_wp& y)
  {
    y.Reallocate(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      y(i) = real(x(i));
  }
  
  
  //! conversion from complex double vector to double vector (we take the real part)
  void CopyVector(const Vector<complex<double> >& x, VectReal_wp& y)
  {
    y.Reallocate(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      y(i) = real(x(i));
  }
  
  
  //! conversion from float vector to complex double vector
  void CopyVector(const Vector<float>& x, VectComplex_wp& y)
  {
    y.Reallocate(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      y(i) = x(i);
  }
  
  
  //! conversion from double vector to complex double vector
  void CopyVector(const Vector<double>& x, VectComplex_wp& y)
  {
    y.Reallocate(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      y(i) = x(i);
  }
  
  
  //! conversion from complex float vector to complex double vector
  void CopyVector(const Vector<complex<float> >& x, VectComplex_wp& y)
  {
    y.Reallocate(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      y(i) = x(i);
  }
  
  
  //! conversion from complex double vector to complex double vector
  void CopyVector(const Vector<complex<double> >& x, VectComplex_wp& y)
  {
    y.Reallocate(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      y(i) = x(i);
  }
  
  
  /*************
   * Vtk files *
   *************/
  
  
#ifdef MONTJOIE_WITH_MPFR
  mpreal swapEndian(mpreal d)
  {
    return d;
  }
#endif
  
  
  //! swapping octets of a real number (conversion big endian <-> little endian)
  double swapEndian(double d)
  {
    double a;
    char *dst = reinterpret_cast<char*>(&a);
    char *src = reinterpret_cast<char*>(&d);
    
    dst[0] = src[7];
    dst[1] = src[6];
    dst[2] = src[5];
    dst[3] = src[4];
    dst[4] = src[3];
    dst[5] = src[2];
    dst[6] = src[1];
    dst[7] = src[0];
    
    return a;
  }


  //! swapping octets of a real number (conversion big endian <-> little endian)
  float swapEndian(float d)
  {
    float a;
    char *dst = reinterpret_cast<char*>(&a);
    char *src = reinterpret_cast<char*>(&d);
    
    dst[0] = src[3];
    dst[1] = src[2];
    dst[2] = src[1];
    dst[3] = src[0];
    
    return a;
  }


  //! swapping octets of an integer (conversion big endian <-> little endian)
  int swapEndian(int d)
  {
    int a;
    char *dst = reinterpret_cast<char*>(&a);
    char *src = reinterpret_cast<char*>(&d);
    
    dst[0] = src[3];
    dst[1] = src[2];
    dst[2] = src[1];
    dst[3] = src[0];
    
    return a;
  }
  
  
  //! writes a scalar solution in vtk format to an output stream
  /*!
    \param[in] x_sol solution
    \param[in] name name of the data
    \param[inout] file_out output stream
    \param[in] double_prec if true the solution is written in double precision
    otherwise in single precision
    \param[out] ascii if true the solution is written in ascii otherwise in binary
    The solution is considered as a scalar data (only one component)
  */
  template<class T>
  void WriteVtk(const Vector<T>& x_sol, const string& name,
		ostream& file_out, bool double_prec, bool ascii)
  {
    if (double_prec)
      {
	file_out.precision(15);
	file_out << "SCALARS " << name <<" double 1" << endl;
      }
    else
      {
	file_out << "SCALARS " << name <<" float 1" << endl;
	file_out.precision(7);
      }
    
    file_out << "LOOKUP_TABLE default" << endl;

    if (ascii)
      for (int i = 0; i < x_sol.GetM(); i++)
	file_out << x_sol(i) << '\n';
    else
      {
	// values are swapped, problem of endian
	WriteBinaryDoubleOrFloat(x_sol, file_out, double_prec, false, true);
	file_out << '\n';
      }    
  }

  
  //! writes a complex solution in vtk format to an output stream
  /*!
    \param[in] x_sol complex solution
    \param[in] name name of the data
    \param[inout] file_out output stream
    \param[in] double_prec if true the solution is written in double precision
    otherwise in single precision
    \param[out] ascii if true the solution is written in ascii otherwise in binary
    The solution is considered as a vector data,
    x-component of the vector will contain the real part
    y-component of the vector will contain the imaginary part
  */
  template<class T>
  void WriteVtk(const Vector<complex<T> >& x_sol, const string& name,
		ostream& file_out, bool double_prec, bool ascii)
  {
    if (double_prec)
      {
	file_out.precision(15);
	file_out << "VECTORS " << name <<" double" << endl;
      }
    else
      {
	file_out << "VECTORS " << name <<" float" << endl;
	file_out.precision(7);
      }

    if (ascii)
      for (int i = 0; i < x_sol.GetM(); i++)
	file_out << real(x_sol(i)) << " " << imag(x_sol(i)) << " 0.0 " << '\n';
    else
      {
	// values are swapped to obtain data with big endian convention
	WriteComplexVectorVtk(x_sol, file_out, double_prec, false, true);
	file_out << '\n';
      }    
  }


  //! writes a vectorial solution in vtk format to an output stream
  /*!
    \param[in] x_sol solution
    \param[in] name name of the data
    \param[inout] file_out output stream
    \param[in] double_prec if true the solution is written in double precision
    otherwise in single precision
    \param[out] ascii if true the solution is written in ascii otherwise in binary
    The solution is considered as a vector data (three components)
  */
  template<class T>
  void WriteVtk(const Vector<TinyVector<T, 3> >& x_sol, const string& name,
		ostream& file_out, bool double_prec, bool ascii)
  {
    if (double_prec)
      {
	file_out.precision(15);
	file_out << "VECTORS " << name <<" double" << endl;
      }
    else
      {
	file_out << "VECTORS " << name <<" float" << endl;
	file_out.precision(7);
      }

    if (ascii)
      for (int i = 0; i < x_sol.GetM(); i++)
	file_out << x_sol(i)(0) << " " << x_sol(i)(1) << " " << x_sol(i)(2) << '\n';
    else
      {
	// values are swapped to satisfy big endian convention
	Vector<T> xvec;
	xvec.SetData(3*x_sol.GetM(), reinterpret_cast<T*>(x_sol.GetData()));
	WriteBinaryDoubleOrFloat(xvec, file_out, double_prec, false, true);
	xvec.Nullify();
	file_out << '\n';
      }
  }

  
  //! writes a vectorial solution in vtk format to an output stream
  /*!
    \param[in] x_sol solution
    \param[in] name name of the data
    \param[inout] file_out output stream
    \param[in] double_prec if true the solution is written in double precision
    otherwise in single precision
    \param[out] ascii if true the solution is written in ascii otherwise in binary
    The solution is considered as a vector data (three components)
    the real part and imaginary part are written as two different vectors
  */
  template<class T>
  void WriteVtk(const Vector<TinyVector<complex<T>, 3> >& x_sol, const string& name,
		ostream& file_out, bool double_prec, bool ascii)
  {
    if (double_prec)
      {
	file_out.precision(15);
	file_out << "VECTORS " << name <<"_real " << " double" << endl;
      }
    else
      {
	file_out << "VECTORS " << name << "_real " << " float" << endl;
	file_out.precision(7);
      }

    int N = x_sol.GetM();
    Vector<T> xvec;
    if (ascii)
      for (int i = 0; i < x_sol.GetM(); i++)
	file_out << real(x_sol(i)(0)) << " " << real(x_sol(i)(1)) << " "
		 << real(x_sol(i)(2)) << '\n';
    else
      {
	xvec.Reallocate(3*N);
	for (int i = 0; i < N; i++)
	  {
	    xvec(3*i) = real(x_sol(i)(0));
	    xvec(3*i+1) = real(x_sol(i)(1));
	    xvec(3*i+2) = real(x_sol(i)(2));
	  }
	
	// values are swapped to obtain data with big endian convention
	WriteBinaryDoubleOrFloat(xvec, file_out, double_prec, false, true);
	file_out << '\n';
      }    

    if (double_prec)
      file_out << "VECTORS " << name <<"_imag " << " double" << endl;
    else
      file_out << "VECTORS " << name << "_imag " << " float" << endl;

    if (ascii)
      for (int i = 0; i < x_sol.GetM(); i++)
	file_out << imag(x_sol(i)(0)) << " " << imag(x_sol(i)(1)) << " "
		 << imag(x_sol(i)(2)) << '\n';
    else
      {
	for (int i = 0; i < N; i++)
	  {
	    xvec(3*i) = imag(x_sol(i)(0));
	    xvec(3*i+1) = imag(x_sol(i)(1));
	    xvec(3*i+2) = imag(x_sol(i)(2));
	  }
	
	// values are swapped to obtain data with big endian convention
	WriteBinaryDoubleOrFloat(xvec, file_out, double_prec, false, true);
	file_out << '\n';
      }    
    
  }

  
  //! writes a tensorial solution in vtk format to an output stream
  /*!
    \param[in] x_sol solution
    \param[in] name name of the data
    \param[inout] file_out output stream
    \param[in] double_prec if true the solution is written in double precision
    otherwise in single precision
    \param[out] ascii if true the solution is written in ascii otherwise in binary
    The solution is considered as a tensor data (3x3 symmetric matrix)
  */
  template<class T>
  void WriteVtk(const Vector<TinyMatrix<T, Symmetric, 3, 3> >& x_sol,
                const string& name,
                ostream& file_out, bool double_prec, bool ascii)
  {
    if (double_prec)
      {
	file_out.precision(15);
	file_out << "TENSORS " << name <<" double" << endl;
      }
    else
      {
	file_out << "TENSORS " << name <<" float" << endl;
	file_out.precision(7);
      }

    if (ascii)
      for (int i = 0; i < x_sol.GetM(); i++)
	{
	  file_out << x_sol(i)(0, 0) << " " << x_sol(i)(0, 1) << " " << x_sol(i)(0, 2) << '\n';
	  file_out << x_sol(i)(1, 0) << " " << x_sol(i)(1, 1) << " " << x_sol(i)(1, 2) << '\n';
	  file_out << x_sol(i)(2, 0) << " " << x_sol(i)(2, 1) << " " << x_sol(i)(2, 2) << '\n';
	  if (i < x_sol.GetM() - 1)
	    file_out << '\n';
	}
    else
      {
	// values are swapped to satisfy big endian convention
	int N = x_sol.GetM();
	Vector<T> xvec(9*N);
	for (int i = 0; i < N; i++)
	  {
	    xvec(9*i) = x_sol(i)(0, 0); xvec(9*i+1) = x_sol(i)(0, 1);
	    xvec(9*i+2) = x_sol(i)(0, 2); xvec(9*i+3) = x_sol(i)(1, 0);
	    xvec(9*i+4) = x_sol(i)(1, 1); xvec(9*i+5) = x_sol(i)(1, 2);
	    xvec(9*i+6) = x_sol(i)(2, 0); xvec(9*i+7) = x_sol(i)(2, 1);
	    xvec(9*i+8) = x_sol(i)(2, 2);
	  }
	
	WriteBinaryDoubleOrFloat(xvec, file_out, double_prec, false, true);
	file_out << '\n';
      }
  }


  //! writes a tensorial solution in vtk format to an output stream
  /*!
    \param[in] x_sol solution
    \param[in] name name of the data
    \param[inout] file_out output stream
    \param[in] double_prec if true the solution is written in double precision
    otherwise in single precision
    \param[out] ascii if true the solution is written in ascii otherwise in binary
    The solution is considered as a tensor data (3x3 symmetric matrix)
  */
  template<class T>
  void WriteVtk(const Vector<TinyMatrix<complex<T>, Symmetric, 3, 3> >& x_sol, const string& name,
		ostream& file_out, bool double_prec, bool ascii)
  {
    if (double_prec)
      {
	file_out.precision(15);
	file_out << "TENSORS " << name << "_real double" << endl;
      }
    else
      {
	file_out << "TENSORS " << name <<"_real float" << endl;
	file_out.precision(7);
      }

    if (ascii)
      for (int i = 0; i < x_sol.GetM(); i++)
	{
	  file_out << real(x_sol(i)(0, 0)) << " " << real(x_sol(i)(0, 1)) 
                   << " " << real(x_sol(i)(0, 2)) << '\n';
	  file_out << real(x_sol(i)(1, 0)) << " " << real(x_sol(i)(1, 1)) 
                   << " " << real(x_sol(i)(1, 2)) << '\n';
	  file_out << real(x_sol(i)(2, 0)) << " " << real(x_sol(i)(2, 1)) 
                   << " " << real(x_sol(i)(2, 2)) << '\n';
	  if (i < x_sol.GetM() - 1)
	    file_out << '\n';
	}
    else
      {
	// values are swapped to satisfy big endian convention
	int N = x_sol.GetM();
	Vector<T> xvec(9*N);
	for (int i = 0; i < N; i++)
	  {
	    xvec(9*i) = real(x_sol(i)(0, 0)); xvec(9*i+1) = real(x_sol(i)(0, 1));
	    xvec(9*i+2) = real(x_sol(i)(0, 2)); xvec(9*i+3) = real(x_sol(i)(1, 0));
	    xvec(9*i+4) = real(x_sol(i)(1, 1)); xvec(9*i+5) = real(x_sol(i)(1, 2));
	    xvec(9*i+6) = real(x_sol(i)(2, 0)); xvec(9*i+7) = real(x_sol(i)(2, 1));
	    xvec(9*i+8) = real(x_sol(i)(2, 2));
	  }
	
	WriteBinaryDoubleOrFloat(xvec, file_out, double_prec, false, true);
	file_out << '\n';
      }

    if (double_prec)
      file_out << "TENSORS " << name << "_imag double" << endl;
    else
      file_out << "TENSORS " << name <<"_imag float" << endl;

    if (ascii)
      for (int i = 0; i < x_sol.GetM(); i++)
	{
	  file_out << imag(x_sol(i)(0, 0)) << " " << imag(x_sol(i)(0, 1)) << " " 
                   << imag(x_sol(i)(0, 2)) << '\n';
	  file_out << imag(x_sol(i)(1, 0)) << " " << imag(x_sol(i)(1, 1)) 
                   << " " << imag(x_sol(i)(1, 2)) << '\n';
	  file_out << imag(x_sol(i)(2, 0)) << " " << imag(x_sol(i)(2, 1)) 
                   << " " << imag(x_sol(i)(2, 2)) << '\n';
	  if (i < x_sol.GetM() - 1)
	    file_out << '\n';
	}
    else
      {
	// values are swapped to satisfy big endian convention
	int N = x_sol.GetM();
	Vector<T> xvec(9*N);
	for (int i = 0; i < N; i++)
	  {
	    xvec(9*i) = imag(x_sol(i)(0, 0)); xvec(9*i+1) = imag(x_sol(i)(0, 1));
	    xvec(9*i+2) = imag(x_sol(i)(0, 2)); xvec(9*i+3) = imag(x_sol(i)(1, 0));
	    xvec(9*i+4) = imag(x_sol(i)(1, 1)); xvec(9*i+5) = imag(x_sol(i)(1, 2));
	    xvec(9*i+6) = imag(x_sol(i)(2, 0)); xvec(9*i+7) = imag(x_sol(i)(2, 1));
	    xvec(9*i+8) = imag(x_sol(i)(2, 2));
	  }
	
	WriteBinaryDoubleOrFloat(xvec, file_out, double_prec, false, true);
	file_out << '\n';
      }
  }


  //! writes a vectorial solution in vtk format to an output stream
  /*!
    \param[in] x_sol solution
    \param[in] type the different data sets to write
    \param[in] decomp for each data set, the component numbers
    \param[in] name associated name with each data set
    \param[in] name name of the data
    \param[inout] file_out output stream
    \param[in] double_prec if true the solution is written in double precision
    otherwise in single precision
    \param[out] ascii if true the solution is written in ascii otherwise in binary
    The solution contains p unknowns, these unknowns are grouped into several 
    data sets (i.e. scalars, vectors or tensors). For each data set, the vtk type
    is given in type(i) (equal to VtkDataTypeEnum::SCALAR, VECTOR or TENSOR),
    the component numbers are given in decomp(i) and the name of the data set is given
    in name(i).
  */
  template<class T>
  void WriteVtk(const Vector<Vector<complex<T> > >& x_sol,
		const IVect& type, const Vector<IVect>& decomp,
		const Vector<string>& name,
		ostream& file_out, bool double_prec, bool ascii)
  {
    if (x_sol.GetM() <= 0)
      return;

    Vector<T> xscal;
    Vector<TinyVector<T, 3> > xvec;
    Vector<TinyMatrix<T, Symmetric, 3, 3> > xtens;
    int N = x_sol(0).GetM();
    // loop over each data set
    for (int k = 0; k < name.GetM(); k++)
      {
	switch (type(k))
	  {
	  case VtkDataTypeEnum::SCALAR :
	    {
	      xscal.Reallocate(N);
	      int m = decomp(k)(0)/2;
	      if (decomp(k)(0)%2 == 0)
		for (int i = 0; i < N; i++)
		  xscal(i) = real(x_sol(m)(i));
	      else
		for (int i = 0; i < N; i++)
		  xscal(i) = imag(x_sol(m)(i));
	      
	      WriteVtk(xscal, name(k), file_out, double_prec, ascii);
	    }
	    break;
	  case VtkDataTypeEnum::VECTOR :
	    {
	      xvec.Reallocate(N);
	      int m = decomp(k)(0)/2;
	      if (m >= 0)
		{
		  if (decomp(k)(0)%2 == 0)
		    for (int i = 0; i < N; i++)
		      xvec(i)(0) = real(x_sol(m)(i));
		  else
		    for (int i = 0; i < N; i++)
		      xvec(i)(0) = imag(x_sol(m)(i));
		}

	      m = decomp(k)(1)/2;
	      if (m >= 0)
		{
		  if (decomp(k)(1)%2 == 0)
		    for (int i = 0; i < N; i++)
		      xvec(i)(1) = real(x_sol(m)(i));
		  else
		    for (int i = 0; i < N; i++)
		      xvec(i)(1) = imag(x_sol(m)(i));
		}
	      
	      m = decomp(k)(2)/2;
	      if (m >= 0)
		{
		  if (decomp(k)(2)%2 == 0)
		    for (int i = 0; i < N; i++)
		      xvec(i)(2) = real(x_sol(m)(i));
		  else
		    for (int i = 0; i < N; i++)
		      xvec(i)(2) = imag(x_sol(m)(i));
		}
	      
	      WriteVtk(xvec, name(k), file_out, double_prec, ascii);
	    }
	    break;
	  case VtkDataTypeEnum::TENSOR :
	    {
	      xtens.Reallocate(N);
	      int m = decomp(k)(0)/2;
	      if (m >= 0)
		{
		  if (decomp(k)(0)%2 == 0)
		    for (int i = 0; i < N; i++)
		      xtens(i)(0, 0) = real(x_sol(m)(i));
		  else
		    for (int i = 0; i < N; i++)
		      xtens(i)(0, 0) = imag(x_sol(m)(i));
		}
	      
	      m = decomp(k)(1)/2;
	      if (m >= 0)
		{	
		  if (decomp(k)(1)%2 == 0)
		    for (int i = 0; i < N; i++)
		      xtens(i)(0, 1) = real(x_sol(m)(i));
		  else
		    for (int i = 0; i < N; i++)
		      xtens(i)(0, 1) = imag(x_sol(m)(i));
		}
	      
	      m = decomp(k)(2)/2;
	      if (m >= 0)
		{
		  if (decomp(k)(2)%2 == 0)
		    for (int i = 0; i < N; i++)
		      xtens(i)(0, 2) = real(x_sol(m)(i));
		  else
		    for (int i = 0; i < N; i++)
		      xtens(i)(0, 2) = imag(x_sol(m)(i));
		}
	      
	      m = decomp(k)(3)/2;
	      if (m >= 0)
		{
		  if (decomp(k)(3)%2 == 0)
		    for (int i = 0; i < N; i++)
		      xtens(i)(1, 1) = real(x_sol(m)(i));
		  else
		    for (int i = 0; i < N; i++)
		      xtens(i)(1, 1) = imag(x_sol(m)(i));
		}
	      
	      m = decomp(k)(4)/2;
	      if (m >= 0)
		{
		  if (decomp(k)(4)%2 == 0)
		    for (int i = 0; i < N; i++)
		      xtens(i)(1, 2) = real(x_sol(m)(i));
		  else
		    for (int i = 0; i < N; i++)
		      xtens(i)(1, 2) = imag(x_sol(m)(i));
		}
	      
	      m = decomp(k)(5)/2;
	      if (m >= 0)
		{
		  if (decomp(k)(5)%2 == 0)
		    for (int i = 0; i < N; i++)
		      xtens(i)(2, 2) = real(x_sol(m)(i));
		  else
		    for (int i = 0; i < N; i++)
		      xtens(i)(2, 2) = imag(x_sol(m)(i));
		}
	      
	      WriteVtk(xtens, name(k), file_out, double_prec, ascii);

	    }
	    break;	    
	  }
      }
  }
  

  //! writes a vectorial solution in vtk format to an output stream
  /*!
    \param[in] x_sol solution
    \param[in] type the different data sets to write
    \param[in] decomp for each data set, the component numbers
    \param[in] name associated name with each data set
    \param[in] name name of the data
    \param[inout] file_out output stream
    \param[in] double_prec if true the solution is written in double precision
    otherwise in single precision
    \param[out] ascii if true the solution is written in ascii otherwise in binary
    The solution contains p unknowns, these unknowns are grouped into several 
    data sets (i.e. scalars, vectors or tensors). For each data set, the vtk type
    is given in type(i) (equal to VtkDataTypeEnum::SCALAR, VECTOR or TENSOR),
    the component numbers are given in decomp(i) and the name of the data set is given
    in name(i).
  */  
  template<class T>
  void WriteVtk(const Vector<Vector<T> >& x_sol,
		const IVect& type,
		const Vector<IVect>& decomp, const Vector<string>& name,
		ostream& file_out, bool double_prec, bool ascii)
  {
    if (x_sol.GetM() <= 0)
      return;
    
    Vector<T> xscal;
    Vector<TinyVector<T, 3> > xvec;
    Vector<TinyMatrix<T, Symmetric, 3, 3> > xtens;
    int N = x_sol(0).GetM();
    // loop over each data set
    for (int k = 0; k < name.GetM(); k++)
      {
	switch (type(k))
	  {
	  case VtkDataTypeEnum::SCALAR :
	    {
	      xscal.Reallocate(N);
	      int m = decomp(k)(0);
	      for (int i = 0; i < N; i++)
		xscal(i) = x_sol(m)(i);
	      
	      WriteVtk(xscal, name(k), file_out, double_prec, ascii);
	    }
	    break;
	  case VtkDataTypeEnum::VECTOR :
	    {
	      xvec.Reallocate(N);
	      int m = decomp(k)(0);
	      if (m >= 0)
		for (int i = 0; i < N; i++)
		  xvec(i)(0) = x_sol(m)(i);
	      
	      m = decomp(k)(1);
	      if (m >= 0)
		for (int i = 0; i < N; i++)
		  xvec(i)(1) = x_sol(m)(i);
	      
	      m = decomp(k)(2);
	      if (m >= 0)
		for (int i = 0; i < N; i++)
		  xvec(i)(2) = x_sol(m)(i);
	      
	      WriteVtk(xvec, name(k), file_out, double_prec, ascii);
	    }
	    break;
	  case VtkDataTypeEnum::TENSOR :
	    {
	      xtens.Reallocate(N);
	      int m = decomp(k)(0);
	      if (m >= 0)
		for (int i = 0; i < N; i++)
		  xtens(i)(0, 0) = x_sol(m)(i);
	      
	      m = decomp(k)(1);
	      if (m >= 0)
		for (int i = 0; i < N; i++)
		  xtens(i)(0, 1) = x_sol(m)(i);
	      
	      m = decomp(k)(2);
	      if (m >= 0)
		for (int i = 0; i < N; i++)
		  xtens(i)(0, 2) = x_sol(m)(i);

	      m = decomp(k)(3);
	      if (m >= 0)
		for (int i = 0; i < N; i++)
		  xtens(i)(1, 1) = x_sol(m)(i);
	      
	      m = decomp(k)(4);
	      if (m >= 0)
		for (int i = 0; i < N; i++)
		  xtens(i)(1, 2) = x_sol(m)(i);

	      m = decomp(k)(5);
	      if (m >= 0)
		for (int i = 0; i < N; i++)
		  xtens(i)(2, 2) = x_sol(m)(i);
	      
	      WriteVtk(xtens, name(k), file_out, double_prec, ascii);
	    }
	    break;	    
	  }
      }
  }
  

  //! writing a complex vector as a three-component vector (x = real, y = imag, z = 0)
  template<class T>
  void WriteComplexVectorVtk(const Vector<complex<T> >& output_vector,
			     ostream& file_out, bool double_prec,
			     bool with_size, bool swap)
  {
    int N = output_vector.GetM();
    if (double_prec)
      {
	Vector<double> output(3*N);
	if (swap)
	  for (int i = 0; i < N; i++)
	    {
	      output(3*i) = swapEndian(double(real(output_vector(i))));
	      output(3*i+1) = swapEndian(double(imag(output_vector(i))));
	      output(3*i+2) = swapEndian(double(0.0));
	    }
	else
	  for (int i = 0; i < N; i++)
	    {
	      output(3*i) = toDouble(realpart(output_vector(i)));
	      output(3*i+1) = toDouble(imagpart(output_vector(i)));
	      output(3*i+2) = 0.0;
	    }
	
	output.Write(file_out, with_size);
      }
    else
      {
	Vector<float> output(3*N);
	if (swap)
	  for (int i = 0; i < N; i++)
	    {
	      output(3*i) = swapEndian(float(real(output_vector(i))));
	      output(3*i+1) = swapEndian(float(imag(output_vector(i))));
	      output(3*i+2) = swapEndian(float(0.0));
	    }
	else
	  for (int i = 0; i < N; i++)
	    {
	      output(3*i) = toDouble(real(output_vector(i)));
	      output(3*i+1) = toDouble(imag(output_vector(i)));
	      output(3*i+2) = 0;
	    }
	
	output.Write(file_out, with_size);
      }
  }

  
  //! writes a vector of integers in binary with eventual swapping
  void WriteBinaryInteger(IVect& output_vector,
			  ostream& file_out, bool with_size, bool swap)
  {
    int N = output_vector.GetM();
    if (swap)
      for (int i = 0; i < N; i++)
	output_vector(i) = swapEndian(output_vector(i));
    
    output_vector.Write(file_out, with_size);
  }
  
  
  //! writing a real vector in single or double precision
  template<class real, class Allocator>
  void WriteBinaryDoubleOrFloat(const Vector<real, VectFull, Allocator>& output_vector,
				ostream& file_out, bool double_prec, bool with_size, bool swap)
  {
    int N = output_vector.GetM();
    if (double_prec)
      {
	Vector<double> output(N);
	if (swap)
	  for (int i = 0; i < N; i++)
	    output(i) = swapEndian(double(output_vector(i)));
	else
	  for (int i = 0; i < N; i++)
	    output(i) = double(output_vector(i));
	
	output.Write(file_out, with_size);
      }
    else
      {
	Vector<float> output(N);
	if (swap)
	  for (int i = 0; i < N; i++)
	    output(i) = swapEndian(float(output_vector(i)));
	else
	  for (int i = 0; i < N; i++)
	    output(i) = float(output_vector(i));
	
	output.Write(file_out, with_size);
      }
  }
  

  //! writing a complex vector in single or double precision
  template<class reel, class Allocator>
  void WriteBinaryDoubleOrFloat(const Vector<complex<reel>, VectFull, Allocator>& output_vector,
				ostream& file_out, bool double_prec, bool with_size, bool swap)
  {
    int N = output_vector.GetM();
    if (double_prec)
      {
	Vector<complex<double> > output(N);
	for (int i = 0; i < N; i++)
	  output(i) = complex<double>(double(real(output_vector(i))), double(imag(output_vector(i))));
	
	output.Write(file_out, with_size);
      }
    else
      {
	Vector<complex<float> > output(N);
	for (int i = 0; i < N; i++)
	  output(i) = complex<float>(float(real(output_vector(i))), float(imag(output_vector(i))));
	
	output.Write(file_out, with_size);
      }
  }
  
  
  //! writing a real vector in single or double precision
  template<class real, int N>
  void WriteBinaryDoubleOrFloat(const TinyVector<real, N>& output_vector,
				ostream& file_out, bool double_prec, bool with_size, bool swap)
  {
    if (double_prec)
      {
	Vector<double> output(N);
	if (swap)
	  for (int i = 0; i < N; i++)
	    output(i) = swapEndian(double(output_vector(i)));
	else
	  for (int i = 0; i < N; i++)
	    output(i) = double(output_vector(i));
	
	output.Write(file_out, with_size);
      }
    else
      {
	Vector<float> output(N);
	if (swap)
	  for (int i = 0; i < N; i++)
	    output(i) = swapEndian(float(output_vector(i)));
	else
	  for (int i = 0; i < N; i++)
	    output(i) = float(output_vector(i));
	
	output.Write(file_out, with_size);
      }
  }
  
  
  //! writing a complex vector in single or double precision
  template<class real, int N>
  void WriteBinaryDoubleOrFloat(const TinyVector<complex<real>, N>& output_vector,
				ostream& file_out, bool double_prec, bool with_size, bool swap)
  {
    if (double_prec)
      {
	Vector<complex<double> > output(N);
	for (int i = 0; i < N; i++)
	  output(i) = complex<double>(double(real(output_vector(i))), double(imag(output_vector(i))));
	
	output.Write(file_out, with_size);
      }
    else
      {
	Vector<complex<float> > output(N);
	for (int i = 0; i < N; i++)
	  output(i) = complex<float>(float(real(output_vector(i))), float(imag(output_vector(i))));
	
	output.Write(file_out, with_size);
      }
  }
  
  
  //! reads an array of complex numbers from vtk file
  template<class T>
  void ReadVtk(Vector<complex<T> >& x_sol, string& name, istream& file_in, bool ascii)
  {
    string keyword, dataType, ligne;
    file_in >> keyword >> name >> dataType;
    if (keyword != "VECTORS")
      {
	cout << "Incorrect type of data" << endl;
	abort();
      }

    bool double_prec = false;
    if (dataType == "float")
      double_prec = false;
    else if (dataType == "double")
      double_prec = true;
    else
      {
	cout << "Type of data must be float or double" << endl;
	abort();
      }
    
    if (ascii)
      {
	T vx, vy, vz;
	for (int i = 0; i < x_sol.GetM(); i++)
	  {
	    file_in >> vx >> vy >> vz;
	    x_sol(i) = complex<T>(vx, vy);
	  }	
      }
    else
      {
	// datas written in big endian convention
	getline(file_in, ligne);
	Vector<T> xvec(3*x_sol.GetM());
	ReadBinaryDoubleOrFloat(xvec, file_in, double_prec, false, true);
	
	for (int i = 0; i < x_sol.GetM(); i++)
	  x_sol(i) = complex<T>(xvec(3*i), xvec(3*i+1));
      }
  }
  
  
  //! reads an array of real number from vtk file
  template<class T>
  void ReadVtk(Vector<T>& x_sol, string& name, istream& file_in, bool ascii)
  {
    string keyword, dataType, ligne; int ncomp;
    file_in >> keyword >> name >> dataType >> ncomp;
    if (keyword != "SCALARS")
      {
	cout << "Incorrect type of data" << endl;
	abort();
      }
    
    // reading LOOKUP_TABLE tableName
    file_in >> keyword >> ligne;
    
    // precision ?
    bool double_prec = false;
    if (dataType == "float")
      double_prec = false;
    else if (dataType == "double")
      double_prec = true;
    else
      {
	cout << "Type of data must be float or double" << endl;
	abort();
      }
    
    // reading values
    if (ascii)
      {
	for (int i = 0; i < x_sol.GetM(); i++)
	  file_in >> x_sol(i);

      }
    else
      {
	// datas written in big endian convention
	getline(file_in, ligne);
	ReadBinaryDoubleOrFloat(x_sol, file_in, double_prec, false, true);
      }
  }
  
  
  //! reads an array of small complex vectors from a vtk file
  template<class T>
  void ReadVtk(Vector<TinyVector<complex<T>, 3> >& x_sol,
	       string& name, istream& file_in, bool ascii)
  {
    Vector<TinyVector<T, 3> > xreal(x_sol.GetM()), ximag(x_sol.GetM());
    ReadVtk(xreal, name, file_in, ascii);
    ReadVtk(ximag, name, file_in, ascii);
    
    for (int i = 0; i < x_sol.GetM(); i++)
      x_sol(i).Init(complex<T>(xreal(i)(0), ximag(i)(0)),
		    complex<T>(xreal(i)(1), ximag(i)(1)),
		    complex<T>(xreal(i)(2), ximag(i)(2)));
  }
  
  
  //! reads an array of tiny vectors from a vtk file
  template<class T>
  void ReadVtk(Vector<TinyVector<T, 3> >& x_sol,
	       string& name, istream& file_in, bool ascii)
  {
    string keyword, dataType, ligne;
    file_in >> keyword >> name >> dataType;
    if (keyword != "VECTORS")
      {
	cout << "Incorrect type of data" << endl;
	abort();
      }
    
    // precision ?
    bool double_prec = false;
    if (dataType == "float")
      double_prec = false;
    else if (dataType == "double")
      double_prec = true;
    else
      {
	cout << "Type of data must be float or double" << endl;
	abort();
      }
    
    // reading values
    if (ascii)
      {
	for (int i = 0; i < x_sol.GetM(); i++)
	  file_in >> x_sol(i)(0) >> x_sol(i)(1) >> x_sol(i)(2);	
      }
    else
      {
	// datas written in big endian convention
	getline(file_in, ligne);
	Vector<T> xvec(3*x_sol.GetM());
	ReadBinaryDoubleOrFloat(xvec, file_in, double_prec, false, true);
	
	for (int i = 0; i < x_sol.GetM(); i++)
	  x_sol(i).Init(xvec(3*i), xvec(3*i+1), xvec(3*i+2));
      }
  }
  
  
  //! Reads an array of tiny complex matrices from a vtk file
  template<class T>
  void ReadVtk(Vector<TinyMatrix<complex<T>, Symmetric, 3, 3> >& x_sol,
	       string& name, istream& file_in, bool ascii)
  {
    Vector<TinyMatrix<T, Symmetric, 3, 3> > xreal(x_sol.GetM()), ximag(x_sol.GetM());
    ReadVtk(xreal, name, file_in, ascii);
    ReadVtk(ximag, name, file_in, ascii);
    
    for (int i = 0; i < x_sol.GetM(); i++)
      {
	x_sol(i)(0, 0) = complex<T>(xreal(i)(0, 0), ximag(i)(0, 0));
	x_sol(i)(0, 1) = complex<T>(xreal(i)(0, 1), ximag(i)(0, 1));
	x_sol(i)(0, 2) = complex<T>(xreal(i)(0, 2), ximag(i)(0, 2));
	x_sol(i)(1, 1) = complex<T>(xreal(i)(1, 1), ximag(i)(1, 1));
	x_sol(i)(1, 2) = complex<T>(xreal(i)(1, 2), ximag(i)(1, 2));
	x_sol(i)(2, 2) = complex<T>(xreal(i)(2, 2), ximag(i)(2, 2));
      }
  }
  
  
  //! Reads an array of tiny matrices from a vtk file
  template<class T>
  void ReadVtk(Vector<TinyMatrix<T, Symmetric, 3, 3> >& x_sol,
	       string& name, istream& file_in, bool ascii)
  {
    string keyword, dataType, ligne;
    file_in >> keyword >> name >> dataType;
    if (keyword != "TENSORS")
      {
	cout << "Incorrect type of data" << endl;
	abort();
      }
    
    // precision ?
    bool double_prec = false;
    if (dataType == "float")
      double_prec = false;
    else if (dataType == "double")
      double_prec = true;
    else
      {
	cout << "Type of data must be float or double" << endl;
	abort();
      }
    
    // reading values
    if (ascii)
      {
	for (int i = 0; i < x_sol.GetM(); i++)
	  {
	    file_in >> x_sol(i)(0, 0) >> x_sol(i)(0, 1) >> x_sol(i)(0, 2);
	    file_in >> x_sol(i)(1, 0) >> x_sol(i)(1, 1) >> x_sol(i)(1, 2);
	    file_in >> x_sol(i)(2, 0) >> x_sol(i)(2, 1) >> x_sol(i)(2, 2);
	  }
      }
    else
      {
	// datas written in big endian convention
	getline(file_in, ligne);
	Vector<T> xvec(9*x_sol.GetM());
	ReadBinaryDoubleOrFloat(xvec, file_in, double_prec, false, true);
	
	for (int i = 0; i < x_sol.GetM(); i++)
	  {
	    int j = 9*i;
	    x_sol(i)(0, 0) = xvec(j);
	    x_sol(i)(0, 1) = xvec(j+1);
	    x_sol(i)(0, 2) = xvec(j+2);
	    x_sol(i)(1, 1) = xvec(j+4);
	    x_sol(i)(1, 2) = xvec(j+5);
	    x_sol(i)(2, 2) = xvec(j+8);
	  }
      }
  }

  
  //! reading of a vector of integers
  void ReadBinaryInteger(IVect& output_vector, istream& file_in,
			 bool with_size, bool swap)
  {
    long N = output_vector.GetM();
    Vector<int> output;
    if (!with_size)
      output.Reallocate(N);
    
    output.Read(file_in, with_size);
    if (swap)
      for (int i = 0; i < min(N, output.GetM()); i++)
	output_vector(i) = swapEndian(output(i));
  }

  
  //! reading of any vector in single or double precision
  template<class real, class Allocator>
  void ReadBinaryDoubleOrFloat(Vector<real, Vect_Full, Allocator>& output_vector,
			       istream& file_in, bool double_prec, bool with_size, bool swap)
  {
    long N = output_vector.GetM();
    if (double_prec)
      {
	Vector<double> output;
	if (!with_size)
	  output.Reallocate(N);
	
	output.Read(file_in, with_size);
	if (swap)
	  for (int i = 0; i < min(N, output.GetM()); i++)
	    output_vector(i) = swapEndian(output(i));
	else
	  for (int i = 0; i < min(N, output.GetM()); i++)
            output_vector(i) = output(i);
      }
    else
      {
	Vector<float> output;
	if (!with_size)
	  output.Reallocate(N);
	
	output.Read(file_in, with_size);
	if (swap)
	  for (int i = 0; i < min(N, output.GetM()); i++)
	    output_vector(i) = swapEndian(output(i));
	else
	  for (int i = 0; i < min(N, output.GetM()); i++)
	    output_vector(i) = output(i);
      }
  }
  
  
  //! reading of any complex vector in single or double precision
  template<class real, class Alloc>
  void ReadBinaryDoubleOrFloat(Vector<complex<real>, Vect_Full, Alloc>& output_vector,
			       istream& file_in, bool double_prec, bool with_size, bool swap)
  {
    long N = output_vector.GetM();
    if (double_prec)
      {
	Vector<complex<double> > output;
	if (!with_size)
	  output.Reallocate(N);
	
	output.Read(file_in, with_size);
	for (int i = 0; i < min(N, output.GetM()); i++)
	  output_vector(i) = output(i);
      }
    else
      {
	Vector<complex<float> > output;
	if (!with_size)
	  output.Reallocate(N);
	
	output.Read(file_in, with_size);
	for (int i = 0; i < min(N, output.GetM()); i++)
	  output_vector(i) = output(i);
      }
  }
  
  
  //! Reads a line by taking care of \r\n problem
  istream& getlineM(istream& file_in, string& line)
  {
    // reading a line of the file
    getline(file_in, line, '\n');
    
    // removing the character 13 to handle dos files
    if (line.size() > 0)
      {
	if (line[line.size()-1] == 13)
	  line = line.substr(0, line.size()-1);
      }
    
    return file_in;
  }
  

  //! the data file is read and stored in the array of strings all_lines
  /*!
    \param[in] file_name name of the data file to read
    \param[out] all_lines all the lines of the data file
    \param[in] comm MPI communicator for processor that will share the data file
    Only the 0-ranked processor reads the data file, and sends its content 
    to the other processors of the communicator
  */
#ifdef SELDON_WITH_MPI
  void ReadLinesFile(const string& file_name, Vector<string>& all_lines,
			    const MPI_Comm& comm)
#else
  void ReadLinesFile(const string& file_name, Vector<string>& all_lines)
#endif
  {
#ifdef SELDON_WITH_MPI
    // we read all the lines of the file in the root processor
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (rank_proc == 0)
      {
#endif
        // attempt to open the input file
        ifstream file_in(file_name.data());
        if (!file_in.is_open())
          {
            cout << "File with datas not found" << endl;
            abort();
          }
    
        string line;  all_lines.Clear();
        // reading the file until the end
        while (!file_in.eof())
          {
	    getlineM(file_in, line);
            all_lines.PushBack(line);
          }
    
#ifdef SELDON_WITH_MPI
      }
#endif
    
#ifdef SELDON_WITH_MPI
    // broadcasting the lines to the other processors
    if (nb_proc > 1)
      {
        int n = all_lines.GetM();
        MPI_Bcast(&n, 1, MPI_INTEGER, 0, comm); 

        if (rank_proc != 0)
          all_lines.Reallocate(n);
        
        for (int i = 0; i < n; i++)
          MPI_Bcast_string(all_lines(i), 0, comm);
      }
#endif
    
  }

  
  //! reading of the input file
  /*!
    \param[in] file_name name of the input file
    \param[in,out] vars generic problem to solve
    \return 0 if successful
   */
  void ReadInputFile(const string& file_name, InputDataProblem_Base& vars)
  {
    // we read all the lines
    Vector<string> all_lines;
    ReadLinesFile(file_name, all_lines);
    
    // then treating them
    ReadInputFile(all_lines, vars);
  }
  
  
  //! analysing the lines of the data file with the object vars
  void ReadInputFile(const Vector<string>& all_lines,
		     InputDataProblem_Base& vars)
  {
    string line;
    string description_field;
    VectString parameters, words;
   
    // reading the file until the end
    int n = 0;
    while (n < all_lines.GetM())
      {
        line = all_lines(n); n++;
        
	// if a line begins with #, it is a comment (ignored)
	if ((line.size() > 0)&&(line[0]!='#'))
	  {
            description_field = string("");
            
	    // the line is split in words
            StringTokenize(line, words, string(" \t"));
	    
	    // first word of the line -> description_field
            if (words.GetM() > 0)
              description_field = words(0);
	    
	    // all the words after '=' are retrieved and stored in parameters
	    // For example, if a line of the input file is equal to
	    // GridBoundingBox = -1 1 -1 1
	    // description_field will be equal to GridBoundingBox
	    // and parameters(0) = -1, parameters(1) = 1 ...
            parameters.Clear();
            if (words.GetM() > 2)
              {
                parameters.Reallocate(words.GetM()-2);
                for (int i = 2; i < words.GetM(); i++)
                  parameters(i-2) = words(i);
                
                int Nw = words.GetM()-1;
		// if the line ends with \, we look at the following line
                while ((words.GetM() > 0) && (words(Nw)[words(Nw).size()-1] == '\\'))
                  {
                    int nb_old = parameters.GetM()-1;
                    if (words(Nw).size() > 1)
                      {
                        parameters(nb_old) = words(Nw).substr(0, words(Nw).size()-1);
                        nb_old++;
                      }
                    
                    line = all_lines(n); n++;
                    StringTokenize(line, words, string(" \t"));
                    
                    int nb_param = nb_old + words.GetM();
                    parameters.Resize(nb_param);
                    
                    for (int i = 0; i < words.GetM(); i++)
                      parameters(nb_old+i) = words(i);
                    
                    Nw = words.GetM() - 1;
                  }
	      }
	    
	    // keyword and associated values are treated by the
	    // method SetInputData that can be overloaded in derived classes
            if (parameters.GetM() > 0)
	      vars.SetInputData(description_field, parameters);
          }
        
      } // end of lines
  }

  
  //! Ecriture d'un nouveau dossier dans [STIFFOUT]
  void EcritDossier(string & DOSSIER, string name_data_file, int num)
  {    
    if (num >= 1)
      {
        DOSSIER += to_str(num) + "/";
        return;
      }
    
    // trouver à quel numéro on est
    int success_read = system("wc -l stiffindex.dat > tmp.dat");
    int intnblignes = 1;
    if (success_read != 0)
      {
	cout << "Unable to count lines " << endl;
      }
    else
      {
        string nblignes;
        ifstream tmp2("tmp.dat");
        while (!tmp2.eof())
          {
            char ligne[256];
            string s;
            tmp2.getline (ligne, sizeof (ligne));
            s = ligne;    
            // DISP(ligne);
            int nom = s.find("index");
            if (nom >0)
              {
                nblignes = s.substr(0, nom-1);
                intnblignes = to_num<int>(nblignes);
                intnblignes += 1;
              }
          }
        
        tmp2.close();
      }
    
    // faire un nouveau dossier avec le bon numéro
    ofstream f("stiffindex.dat", ifstream::app);
    ostringstream nbl;
    nbl << intnblignes;
    string tmpp = nbl.str() + " : " + DOSSIER + "\n";
    f << tmpp.data();
    f.close();
    DISP(tmpp);

    string dos1 = DOSSIER ;
    string ordre="mkdir " +dos1 ;
    DISP(ordre);
        success_read = system(ordre.data());
     if (success_read != 0)
      {
	cout << "Unable to create the directory" << dos1 << endl;
	//abort();
      }
    
    // string DOSSIER1;
    DOSSIER = DOSSIER +  nbl.str() + "/";
    ordre="mkdir " + DOSSIER;
    success_read = system(ordre.data());
    if (success_read != 0)
      {
	cout << "Unable to create the directory" << DOSSIER << endl;
	//abort();
      }
    
    size_t num_slash = name_data_file.find_last_of("/");
    string s_slash = name_data_file;
    if (num_slash != string::npos)
      s_slash = s_slash.substr(num_slash+1);
    
    ordre = "cp " + name_data_file + " " + DOSSIER + s_slash;
    success_read = system(ordre.data());    
    if (success_read != 0)
      {
	cout << "Unable to copy the file " << name_data_file << endl;
	abort();
      }
  }
  
  
  /**************************
   * Outputs with buffering *
   **************************/

  
  //! default constructor
  template<class T, int nb>
  WriteOnTheGoWithTinyBuffer<T, nb>::WriteOnTheGoWithTinyBuffer()
  {
    size_buffer = 100;
    iter = 0;
    is_binary = false;
    is_double = true;
  }
  
  
  template<class T, int nb>
  size_t WriteOnTheGoWithTinyBuffer<T, nb>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille += Buffer.GetMemorySize();
    taille += Seldon::GetMemorySize(name_file);
    return taille;
  }


  //! initialisation of the buffer
  /*!
    \param[in] data_file the destination file of the buffer
    \param[in] size size of the buffer (consecutive values stored before writing them in the file)
    \param[in] remove_file if true, the file is erased before rewriting it,
    otherwise the new values are simply appended at the end of the file
  */
  template<class T, int nb>
  void WriteOnTheGoWithTinyBuffer<T, nb>::Init(const string &data_file, int size,
                                               bool remove_file, bool write_prec)
  {
    name_file = data_file;
    if (remove_file)
      std::remove(name_file.data());
    
    size_buffer = size;
    iter = 0;
    Buffer.Reallocate(size_buffer);
    // precision is written for binary files
    if (write_prec && is_binary)
      WritePrecision();
  }
  
  
  //! Pushes a tiny vector into the buffer
  /*!
    When the buffer is full, contents of the buffer are appended to the output file,
    and buffer is emptied, ready for the next push.
  */
  template<class T, int nb>
  void WriteOnTheGoWithTinyBuffer<T, nb>
  ::AddTinyVect(TinyVector<T, nb> &data)
  {
    Buffer(iter) = data;
    iter++;
    if (iter == size_buffer)
      {
	ofstream file_out(name_file.data(), ios::app);
        file_out.setf(ios::scientific);
        if (is_binary)
          {
            for (int p = 0; p < size_buffer; p++)
              WriteBinaryDoubleOrFloat(Buffer(p), file_out, is_double, true, false);
          }
        else
          {
            if (is_double)
              file_out.precision(15);
            else
              file_out.precision(7);
            
            for (int p = 0 ; p < size_buffer; p++)
              {
                PrintNoBrackets(file_out, Buffer(p));
                file_out << '\n';
              }
          }
        
	file_out.close();
	iter = 0;
      }    
  }
                                
  
  //! values are written in binary or in ascii format
  template<class T, int nb>
  void WriteOnTheGoWithTinyBuffer<T, nb>::SetBinary(bool i_want_binary)
  {
    is_binary = i_want_binary;
  }
  
  
  //! values are written in single or double precision
  template<class T, int nb>
  void WriteOnTheGoWithTinyBuffer<T, nb>::SetDoublePrecision(bool i_want_double)
  { 
    is_double = i_want_double;
  }
  
  
  //! type of data is written in binary file
  template<class T, int nb>
  void WriteOnTheGoWithTinyBuffer<T, nb>::WritePrecision()
  {
    int itmp = 0;
    if (is_double)
      itmp = 1;
    
    ofstream file_out(name_file.data(), ios::app);
    file_out.write(reinterpret_cast<char*>(const_cast<int*>(&itmp)), sizeof(int));
    file_out.close();
  }
                                
                                
  //! Buffer is emptied, writing its contents into the output file
  template<class T, int nb>
  void WriteOnTheGoWithTinyBuffer<T, nb>::CloseBuffer()
  {
    if (iter == 0)
      return;
    
    ofstream file_out(name_file.data(), ios::app);
    file_out.setf(ios::scientific);
    if (is_binary)
      {
        for (int p = 0; p < iter; p++)
          WriteBinaryDoubleOrFloat(Buffer(p), file_out, is_double, true, false);
      }
    else 
      {
        if (is_double)
          file_out.precision(15);
        else
          file_out.precision(7);

        for (int p = 0 ; p < iter; p++)
          {
            PrintNoBrackets(file_out, Buffer(p));
            file_out << '\n';
          }
      }
    
    file_out.close();         
    
    iter = 0;
  }
  
  
  //! destructor
  template<class T, int nb>
  WriteOnTheGoWithTinyBuffer<T, nb>::~WriteOnTheGoWithTinyBuffer()
  {
    CloseBuffer();
  }
  
 
  //! default constructor
  template<class T>
  WriteOnTheGoWithBuffer<T>::WriteOnTheGoWithBuffer()
  {
    size_buffer = 100;
    iter = 0;
    is_binary = false;
    is_double = true;
  }
  

  template<class T>
  size_t WriteOnTheGoWithBuffer<T>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille += Seldon::GetMemorySize(Buffer);
    taille += Seldon::GetMemorySize(name_file);
    return taille;
  }


  //! initialisation of the buffer
  /*!
    \param[in] data_file the destination file of the buffer
    \param[in] size size of the buffer (consecutive values stored before writing them in the file)
    \param[in] remove_file if true, the file is erased before rewriting it,
    otherwise the new values are simply appended at the end of the file
  */  
  template<class T>
  void WriteOnTheGoWithBuffer<T>::Init(const string &data_file, int size,
                                       bool remove_file, bool write_prec)
  {
    name_file	  = data_file;
    if (remove_file)
      std::remove(name_file.data());
    
    size_buffer = size;
    iter = 0;
    Buffer.Reallocate(size_buffer);
    // precision is written for binary files
    if (is_binary && write_prec)
      WritePrecision();
  }


  //! Pushes a vector of values into the buffer
  /*!
    When the buffer is full, contents of the buffer are appended to the output file,
    and buffer is emptied, ready for the next push.
  */    
  template<class T>
  void WriteOnTheGoWithBuffer<T>::AddVect(Vector<T> &data)
  {
    Buffer(iter) = data;
    iter++;
    if (iter == size_buffer)
      {
	ofstream file_out(name_file.data(), ios::app);
        file_out.setf(ios::scientific);
	if (is_binary)
          {
            for (int p = 0; p < size_buffer; p++)
              WriteBinaryDoubleOrFloat(Buffer(p), file_out, is_double, true, false);
          }
        else 
          {
            if (is_double)
              file_out.precision(15);
            else
              file_out.precision(7);
            
            for (int p = 0 ; p < size_buffer; p++)
              {
                file_out << Buffer(p);
                file_out << '\n';
              }
          }  
        
        file_out.close(); 
	iter = 0;
      }
    
  }
  

  //! Pushes a vector of values into the buffer
  /*!
    When the buffer is full, contents of the buffer are appended to the output file,
    and buffer is emptied, ready for the next push.
  */    
  template<class T>
  void WriteOnTheGoWithBuffer<T>::AddVect(Vector<complex<T> > &data)
  {
    cout << "not implemented" << endl;
    abort();
  }
  
  
  //! values are written in binary or in ascii format
  template<class T>
  void WriteOnTheGoWithBuffer<T>::SetBinary(bool i_want_binary)
  {
    is_binary = i_want_binary;
  }
  
  
  //! values are written in single or double precision
  template<class T>
  void WriteOnTheGoWithBuffer<T>::SetDoublePrecision(bool i_want_double)
  { 
    is_double = i_want_double;
  }
  
  
  //! type of data is written in binary file
  template<class T>
  void WriteOnTheGoWithBuffer<T>::WritePrecision()
  {
    int itmp = 0;
    if (is_double)
      itmp = 1;
    
    ofstream file_out(name_file.data(), ios::app);
    file_out.write(reinterpret_cast<char*>(const_cast<int*>(&itmp)), sizeof(int));
    file_out.close();
  }
  
  
  //! Buffer is emptied, writing its contents into the output file  
  template<class T>
  void WriteOnTheGoWithBuffer<T>::CloseBuffer()
  {
    if (iter == 0)
      return;
    
    ofstream file_out(name_file.data(), ios::app);
    file_out.setf(ios::scientific);
    if (is_binary)
      {
        for (int p = 0; p < iter; p++)
          WriteBinaryDoubleOrFloat(Buffer(p), file_out, is_double, true, false);
      }
    else 
      {
        if (is_double)
          file_out.precision(15);
        else
          file_out.precision(7);

        for (int p = 0 ; p < iter; p++)
          {
            file_out << Buffer(p);
            file_out << '\n';
          }
      }    
    
    file_out.close(); 
    
    iter = 0;
  }
  
  
  //! destructor
  template<class T>	
  WriteOnTheGoWithBuffer<T>::~WriteOnTheGoWithBuffer()
  {
    CloseBuffer();
  }
  
}

#define MONTJOIE_FILE_COMMON_INPUT_OUTPUT_CXX
#endif
  
