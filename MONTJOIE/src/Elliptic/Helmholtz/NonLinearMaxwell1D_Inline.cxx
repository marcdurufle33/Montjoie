#ifndef MONTJOIE_FILE_NON_LINEAR_MAXWELL_1D_INLINE_CXX

namespace Montjoie
{

  //! returns the final time
  inline Real_wp MaxwellProblem1D::GetTimeStep() const
  {
    return dt; 
  }
  
  
  //! returns the initial time
  inline Real_wp MaxwellProblem1D::GetInitialTime() const
  {
    return initial_time;
  }
  
  
  //! returns the final time
  inline Real_wp MaxwellProblem1D::GetFinalTime() const
  {
    return final_time;
  }
  
  
  //! returns the number of degrees of freedom
  inline int MaxwellProblem1D::GetNbDof() const
  {
    return nb_points_z; 
  }
  
  
  //! returns E(t^n)
  inline Vector<Complex_wp>& MaxwellProblem1D::GetIterate()
  {
    return E_current;
  }
  
  
  inline void MaxwellProblem1D
  ::SetDirichletCondition(const Real_wp& t, int n, VectComplex_wp& Y,
			  Real_wp alpha)
  {
  }
  
  
  //! returns the file name used for writing Fourier transform of E
  inline string MaxwellProblem1D::GetFourierOutputFile() const
  {
    return var_laplace.GetOutputParameters().GetDiffractedFieldFile();
  }
  
  
  //! finalizing output of seismogramms if needed
  inline void MaxwellProblem1D::CloseBuffers()
  {
    for (int i = 0; i < En_sismo.GetM(); i++)
      if (En_sismo(i).GetM() > 0)
	WriteSismo(i);
  }
  
}

#define MONTJOIE_FILE_NON_LINEAR_MAXWELL_1D_INLINE_CXX
#endif
