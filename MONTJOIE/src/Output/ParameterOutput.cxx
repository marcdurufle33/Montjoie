#ifndef MONTJOIE_FILE_PARAMETER_OUTPUT_CXX

namespace Montjoie
{
  //! default constructor
  ParamOutputClass::ParamOutputClass()
  {
    t_begin = 0.0;
    t_end = 0.0;
    t_current = 0.0;
    deltat = 0.0;
    nature = 0;
    component = 0;
    num_instantane = 0;
    name_file.Reallocate(2);
    sismo_point = false;
    epsilon_time = 1e-10;
    size_buffer_sismo = 100;
    first_output_sismo = true;
    seismogramm_to_write_in_binary = false;
  }

  
  void ParamOutputClass::SetInterval(const Real_wp& t0, const Real_wp& tf, const Real_wp& dt,
                                     Real_wp eps)
  {
    t_begin = t0;
    t_end = tf;
    deltat = dt;    
    if (eps < 0)
      epsilon_time = abs(tf - t0)*1e-14;
    else
      epsilon_time = eps;
  }
    
  
  void ParamOutputClass::SetBinaryOutput(bool b)
  {
    seismogramm_to_write_in_binary = b;
  }
  
  
  void ParamOutputClass::SetBufferSize(int s)
  {
    size_buffer_sismo = s;
  }
  

  size_t ParamOutputClass::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + Seldon::GetMemorySize(name_file)
      + Seldon::GetMemorySize(real_result) + Seldon::GetMemorySize(cplx_result)
      + output_sismo.GetMemorySize() + output_sismo_diff.GetMemorySize();
    
    return taille;
  }


  //! changing parameters of the ouput (data file)
  /*!
    \param[in] name output name
    \param[in] param values of the line of the data file
    \param[in] nb number of values
  */
  void ParamOutputClass::SetInputData(const VectString& param, bool sismo)
  {
    if (param.GetM() <= 2)
      {
	cout << "In SetInputData of ParamOutputClass" << endl;
	cout << "SismoPoint needs more parameters, for instance :" << endl;
	cout << "SismoPoint = t0 tf dt" << endl;
	cout << "Current parameters are : " << endl << param << endl;
	abort();
      }

    t_begin = to_num<Real_wp>(param(0));
    t_end = to_num<Real_wp>(param(1));
    deltat = to_num<Real_wp>(param(2));
    epsilon_time = abs(t_end-t_begin)*1e-14;
    if (sismo)
      {
	sismo_point = true;
	if (param.GetM() > 3)
	  {
	    if (param(3) == "BINARY")
	      seismogramm_to_write_in_binary = true;
	    else
	      seismogramm_to_write_in_binary = false;
	  }
	
	if (param.GetM() > 4)
	  size_buffer_sismo = to_num<int>(param(4));
      }
    else
      sismo_point = false;
  }
  
  
  //! Time to store a snapshot ?
  bool ParamOutputClass::SnapshotToStore(const Real_wp& t) const
  {
    if (this->deltat == 0.0)
      return false;
    
    if ((t >= (this->t_current-epsilon_time))&&(t <= (this->t_end+epsilon_time)))
      return true;
    
    return false;
  }
  
  
  void ParamOutputClass::InitResult(const Vector<Real_wp>& x, int n) const
  {
    if (real_result.GetM() != n*name_file.GetM())
      real_result.Reallocate(n*name_file.GetM());
    
    for (int i = 0; i < real_result.GetM(); i++)
      {
        if (real_result(i).GetM() != x.GetM())
          {
            real_result(i).Reallocate(x.GetM());
            real_result(i).Fill(0);
          }
      }
  }
  

  void ParamOutputClass::InitResult(const Vector<Complex_wp>& x, int n) const
  {
    if (cplx_result.GetM() != n*name_file.GetM())
      cplx_result.Reallocate(n*name_file.GetM());
    
    for (int i = 0; i < cplx_result.GetM(); i++)
      {
        if (cplx_result(i).GetM() != x.GetM())
          {
            cplx_result(i).Reallocate(x.GetM());
            cplx_result(i).Fill(0);
          }
      }
  }


  void ParamOutputClass::UpdateResult(int num_file, int m, Vector<Real_wp>& x) const
  {
    int p = m*name_file.GetM() + num_file;
    for (int i = 0; i < x.GetM(); i++)
      {
        real_result(p)(i) += x(i);
        x(i) = real_result(p)(i);
      }
  }
  
  
  void ParamOutputClass::UpdateResult(int num_file, int m, Vector<Complex_wp>& x) const
  {
    int p = m*name_file.GetM() + num_file;
    for (int i = 0; i < x.GetM(); i++)
      {
        cplx_result(p)(i) += x(i);
        x(i) = cplx_result(p)(i);
      }
  }
  

  //! initialization before time iterations
  void ParamOutputClass::InitTime(const Real_wp& t0)
  {
    t_current = t_begin;
    if (t_current < t0)
      {
	t_begin = t0;
	t_current = t0;
      }

    num_instantane = 0;
    if (sismo_point)
      {
	// ofstream file_out(file_total_field.data());
	// DISP(file_total_field);
	remove(name_file(0).data());
	remove(name_file(1).data());
	// file_out<<tinitial<<" "<<endl;
	// file_out.close();
      }
  }
  
  
  void ParamOutputClass::ChangeTime(int nb_iter, const Real_wp& t)
  {
    if (this->deltat == 0.0)
      num_instantane = 0;
    else
      num_instantane = toInteger(ceil((t-this->t_begin)/this->deltat));
    
    this->t_current = this->t_begin + num_instantane*this->deltat;
  }
    
  
  bool ParamOutputClass::GradientToCompute(int Nmax) const
  {
    if (component == -1)
      return true;
    
    if (component >= Nmax)
      return true;
    
    return false;
  }
  

  void ParamOutputClass::InitSismo(const string& DOSSIER_output, const string& name_file,
				   const string& name_file_diff, int add_total_field, bool double_prec)
  {
    if (first_output_sismo)
      {
        sismo_point = true;
	if (double_prec)
	  output_sismo.SetDoublePrecision();
	else
	  output_sismo.SetDoublePrecision(false);
		
	if (seismogramm_to_write_in_binary)
	  output_sismo.SetBinary();
	else
	  output_sismo.SetBinary(false);
	
	output_sismo.Init(DOSSIER_output+name_file, size_buffer_sismo, false);
	
	if (add_total_field != 0)
	  {
	    if (double_prec)
	      output_sismo_diff.SetDoublePrecision();
	    else
	      output_sismo_diff.SetDoublePrecision(false);
	    
	    if (seismogramm_to_write_in_binary)
	      output_sismo_diff.SetBinary();
	    else
	      output_sismo_diff.SetBinary(false);
	    
	    output_sismo_diff.Init(DOSSIER_output+name_file_diff, size_buffer_sismo, false);
	  }
	
	first_output_sismo = false;
      }
  }            

    
  void ParamOutputClass::AddVectorSismo(VectReal_wp& val, VectReal_wp& val_diff)
  {
    output_sismo.AddVect(val);
    
    if (val_diff.GetM() > 0)
      output_sismo_diff.AddVect(val_diff);
  }


  void ParamOutputClass::AddVectorSismo(VectComplex_wp& val, VectComplex_wp& val_diff)
  {
    cout << "Seismograms with complex values are not implemented" << endl;
    abort();
  }

  
  //! display informations about class ParamOutputClass
  ostream& operator <<(ostream& out, const ParamOutputClass& p)
  {
    out<<"sismos between "<<p.t_begin<<" and "
       <<p.t_end<<" separated of "<<p.deltat<<endl;
    out<<"nature "<<p.nature<<" component "<<p.component<<endl;
    out<<"output files "<<p.name_file<<endl;
    
    return out;
  }
    
} // namespace Montjoie

#define MONTJOIE_FILE_PARAMETER_OUTPUT_CXX
#endif

